#!/usr/bin/env python3
"""Migrate frozen significance CSV files to Verification schema v2.

The migration is deliberately byte-preserving for every pre-existing CSV
field token.  It parses RFC4180-style records at the byte level so quoted
commas, escaped quotes, CRLF, and embedded newlines are retained exactly.
Only the existing VerificationClass token is replaced; all new provenance
fields are appended.

The output directory is an immutable publication target: it must not exist.
All manifest inputs are verified before a sibling staging directory is
created, and the completed staging directory is renamed into place once.
Files with unresolved provenance or incompatible class pairs receive audit
records but no canonical migrated CSV.
"""

from __future__ import annotations

import argparse
import ctypes
import csv
import errno
import hashlib
import io
import json
import os
import re
import secrets
import shlex
import shutil
import stat
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple


SCHEMA_VERSION = 2

APPENDED_HEADERS: Tuple[str, ...] = (
    "VerificationSchemaVersion",
    "VerificationClass_V1_Deprecated",
    "LabelFirstSupport",
    "ClusterFirstSupport",
    "WithinHPSupport",
    "DispersionWarning",
    "EvidencePath",
    "EvidenceDerivation",
    "LOH_Subtype_LegacyVC",
)

REQUIRED_IDENTITY_HEADERS: Tuple[str, ...] = (
    "RegionID",
    "Chr",
    "Pos",
    "Ref",
    "Alt",
)

V1_CURRENT_TO_PATH: Mapping[str, str] = {
    "LOH-Structure": "LOH_STRUCTURE",
    "MultiGroupNoLabel": "WITHIN_HP_MULTIGROUP",
    "LabelShift": "LABEL_SHIFT",
    "PermanovaLocation": "PERMANOVA_LOCATION",
    "StructureNoLabel": "HP_AUC_STRUCTURE_NO_LABEL",
    "DispersionStructure": "DISPERSION_STRUCTURE",
    "Noise_Uniform": "NOISE_UNIFORM",
    "Noise_Chaotic": "NOISE_CHAOTIC",
    "Noise_Uncorrelated": "NOISE_UNCORRELATED",
}

LEGACY_CLASSES = frozenset(("Strong", "Subclone", "Weak", "Noise"))
LOH_SUBTYPES = frozenset(("None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"))
LOH_BY_LEGACY: Mapping[str, str] = {
    "Strong": "LOH_Strong",
    "Subclone": "LOH_Subclone",
    "Weak": "LOH_Weak",
    "Noise": "LOH_Noise",
}

MANIFEST_BASE_HEADERS = ("file", "rows", "sha256")
MANIFEST_COUNT_HEADERS = (
    "input_final_strong",
    "legacy_strong_to_final_strong",
    "legacy_subclone_to_final_strong",
    "legacy_strong_or_subclone_exceptions",
)

SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


class MigrationError(RuntimeError):
    """Fatal contract or I/O error."""


class CsvFormatError(MigrationError):
    """Malformed byte-level CSV input."""


class RowConflict(ValueError):
    """A row cannot be mapped under the authoritative v2 truth table."""


@dataclass(frozen=True)
class RawCsvRecord:
    fields: Tuple[bytes, ...]
    newline: bytes


@dataclass(frozen=True)
class ManifestEntry:
    filename: str
    rows: int
    sha256: str
    expected_counts: Mapping[str, int]


@dataclass(frozen=True)
class VerificationDecision:
    verification_class_v2: str
    evidence_path: str
    label_first_support: str
    cluster_first_support: str


@dataclass
class ValidationResult:
    status: str
    reason: str
    input_rows: int = 0
    before_current: Counter = field(default_factory=Counter)
    before_legacy: Counter = field(default_factory=Counter)
    before_loh: Counter = field(default_factory=Counter)
    after_current: Counter = field(default_factory=Counter)
    after_evidence: Counter = field(default_factory=Counter)
    observed_manifest_counts: Counter = field(default_factory=Counter)
    conflicts: List[Dict[str, str]] = field(default_factory=list)
    stable_key_uniqueness: str = "NOT_RUN"


@dataclass
class FileReport:
    filename: str
    status: str
    reason: str
    input_rows: int
    output_rows: int
    input_sha256: str
    output_sha256: str
    schema_version: str
    unmapped_count: int
    before_current: Mapping[str, int]
    before_legacy: Mapping[str, int]
    before_loh: Mapping[str, int]
    after_current: Mapping[str, int]
    after_evidence: Mapping[str, int]
    raw_token_preservation: str
    significant_stable_key_invariant: str
    stable_key_uniqueness: str


@dataclass(frozen=True)
class StagingIdentity:
    path: Path
    resolved_path: Path
    resolved_parent: Path
    name_prefix: str
    device: int
    inode: int
    owner_uid: int
    parent_device: int
    parent_inode: int


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _call_libc_renameat2_noreplace(source: Path, destination: Path) -> None:
    """Call Linux renameat2(RENAME_NOREPLACE) without a lossy fallback."""

    if not sys.platform.startswith("linux"):
        raise OSError(errno.ENOSYS, "Linux renameat2 is unavailable on this platform")

    libc = ctypes.CDLL(None, use_errno=True)
    try:
        renameat2 = libc.renameat2
    except AttributeError as error:
        raise OSError(errno.ENOSYS, "libc does not expose renameat2") from error

    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    ctypes.set_errno(0)
    result = renameat2(
        -100,  # AT_FDCWD
        os.fsencode(source),
        -100,  # AT_FDCWD
        os.fsencode(destination),
        1,  # RENAME_NOREPLACE
    )
    if result != 0:
        error_number = ctypes.get_errno() or errno.EIO
        raise OSError(error_number, os.strerror(error_number), str(destination))


def atomic_publish_noreplace(staging: Path, output_dir: Path) -> None:
    """Atomically publish staging, failing if no-replace is unavailable or destination exists."""

    try:
        _call_libc_renameat2_noreplace(staging, output_dir)
    except OSError as error:
        if error.errno in (errno.EEXIST, errno.ENOTEMPTY):
            raise MigrationError(
                "output directory appeared during migration; atomic no-replace publication refused: "
                + str(output_dir)
            ) from error
        unsupported = {
            errno.ENOSYS,
            errno.EINVAL,
            errno.EOPNOTSUPP,
            getattr(errno, "ENOTSUP", errno.EOPNOTSUPP),
        }
        if error.errno in unsupported:
            raise MigrationError(
                "atomic no-replace publication is unsupported by libc, kernel, or filesystem; "
                "refusing non-atomic fallback for "
                + str(output_dir)
            ) from error
        raise MigrationError(
            "atomic no-replace publication failed for %s: errno=%s (%s)"
            % (output_dir, error.errno, error.strerror or "unknown error")
        ) from error


def capture_staging_identity(
    staging: Path,
    expected_parent: Path,
    name_prefix: str,
) -> StagingIdentity:
    """Capture the exact mkdtemp directory identity used by fail-safe cleanup."""

    staging = staging.absolute()
    resolved_parent = expected_parent.resolve(strict=True)
    resolved_staging = staging.resolve(strict=True)
    metadata = staging.lstat()
    parent_metadata = resolved_parent.stat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise MigrationError("staging path is not a real directory: " + str(staging))
    if resolved_staging.parent != resolved_parent:
        raise MigrationError("staging directory resolved outside its expected parent")
    if not staging.name.startswith(name_prefix):
        raise MigrationError("staging directory does not have the expected prefix")
    if metadata.st_uid != os.geteuid():
        raise MigrationError("staging directory is not owned by the effective user")
    return StagingIdentity(
        path=staging,
        resolved_path=resolved_staging,
        resolved_parent=resolved_parent,
        name_prefix=name_prefix,
        device=metadata.st_dev,
        inode=metadata.st_ino,
        owner_uid=metadata.st_uid,
        parent_device=parent_metadata.st_dev,
        parent_inode=parent_metadata.st_ino,
    )


def cleanup_owned_staging(identity: StagingIdentity) -> None:
    """Remove only the unchanged mkdtemp directory captured in identity."""

    try:
        metadata = identity.path.lstat()
    except FileNotFoundError:
        return

    resolved_parent = identity.path.parent.resolve(strict=True)
    parent_metadata = resolved_parent.stat()
    if resolved_parent != identity.resolved_parent:
        raise MigrationError("refusing staging cleanup: resolved parent changed")
    if (
        parent_metadata.st_dev != identity.parent_device
        or parent_metadata.st_ino != identity.parent_inode
    ):
        raise MigrationError("refusing staging cleanup: parent identity changed")
    if not identity.path.name.startswith(identity.name_prefix):
        raise MigrationError("refusing staging cleanup: prefix mismatch")
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise MigrationError("refusing staging cleanup: path is not the created directory")
    if metadata.st_uid != identity.owner_uid or metadata.st_uid != os.geteuid():
        raise MigrationError("refusing staging cleanup: ownership changed")
    if metadata.st_dev != identity.device or metadata.st_ino != identity.inode:
        raise MigrationError("refusing staging cleanup: directory identity changed")
    if identity.path.resolve(strict=True) != identity.resolved_path:
        raise MigrationError("refusing staging cleanup: resolved path changed")

    shutil.rmtree(identity.path)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _newline_at(data: bytes, index: int) -> Optional[Tuple[bytes, int]]:
    byte = data[index]
    if byte == 10:
        return b"\n", 1
    if byte == 13:
        if index + 1 < len(data) and data[index + 1] == 10:
            return b"\r\n", 2
        return b"\r", 1
    return None


def iter_csv_records(data: bytes) -> Iterator[RawCsvRecord]:
    """Yield byte-exact CSV records while honoring quoted embedded newlines."""

    fields: List[bytes] = []
    field_start = 0
    index = 0
    in_quotes = False
    closed_quote = False

    while index < len(data):
        byte = data[index]

        if in_quotes:
            if byte == 34:
                if index + 1 < len(data) and data[index + 1] == 34:
                    index += 2
                    continue
                in_quotes = False
                closed_quote = True
            index += 1
            continue

        newline = _newline_at(data, index)

        if closed_quote:
            if byte == 44:
                fields.append(data[field_start:index])
                field_start = index + 1
                closed_quote = False
                index += 1
                continue
            if newline is not None:
                terminator, length = newline
                fields.append(data[field_start:index])
                yield RawCsvRecord(tuple(fields), terminator)
                fields = []
                index += length
                field_start = index
                closed_quote = False
                continue
            raise CsvFormatError(
                "unexpected byte after closing CSV quote at byte offset " + str(index)
            )

        if byte == 34:
            if index != field_start:
                raise CsvFormatError(
                    "quote appears inside an unquoted CSV field at byte offset " + str(index)
                )
            in_quotes = True
            index += 1
            continue

        if byte == 44:
            fields.append(data[field_start:index])
            field_start = index + 1
            index += 1
            continue

        if newline is not None:
            terminator, length = newline
            fields.append(data[field_start:index])
            yield RawCsvRecord(tuple(fields), terminator)
            fields = []
            index += length
            field_start = index
            continue

        index += 1

    if in_quotes:
        raise CsvFormatError("unterminated quoted CSV field at end of file")

    if field_start < len(data) or fields or closed_quote:
        fields.append(data[field_start:len(data)])
        yield RawCsvRecord(tuple(fields), b"")


def decode_csv_token(token: bytes) -> str:
    if token.startswith(b'"'):
        if len(token) < 2 or not token.endswith(b'"'):
            raise CsvFormatError("quoted CSV token does not end with a quote")
        payload = token[1:-1].replace(b'""', b'"')
    else:
        if b'"' in token:
            raise CsvFormatError("unquoted CSV token contains a quote")
        payload = token
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise CsvFormatError("CSV token is not valid UTF-8") from error


def encode_csv_token(value: str, force_quote: bool = False) -> bytes:
    encoded = value.encode("utf-8")
    quote = force_quote or any(marker in encoded for marker in (b",", b'"', b"\r", b"\n"))
    if not quote:
        return encoded
    return b'"' + encoded.replace(b'"', b'""') + b'"'


def header_names(record: RawCsvRecord) -> List[str]:
    names = [decode_csv_token(token) for token in record.fields]
    if names and names[0].startswith("\ufeff"):
        names[0] = names[0].lstrip("\ufeff")
    return names


def count_csv_data_rows(data: bytes) -> int:
    """Count logical data records without allocating every field token."""

    if not data:
        raise CsvFormatError("CSV file is empty and has no header")
    index = 0
    records = 0
    in_quotes = False
    closed_quote = False
    at_field_start = True
    last_record_end = 0

    while index < len(data):
        byte = data[index]
        if in_quotes:
            if byte == 34:
                if index + 1 < len(data) and data[index + 1] == 34:
                    index += 2
                    continue
                in_quotes = False
                closed_quote = True
            index += 1
            continue

        newline = _newline_at(data, index)
        if closed_quote:
            if byte == 44:
                closed_quote = False
                at_field_start = True
                index += 1
                continue
            if newline is not None:
                _terminator, length = newline
                records += 1
                index += length
                last_record_end = index
                closed_quote = False
                at_field_start = True
                continue
            raise CsvFormatError(
                "unexpected byte after closing CSV quote at byte offset " + str(index)
            )

        if byte == 34:
            if not at_field_start:
                raise CsvFormatError(
                    "quote appears inside an unquoted CSV field at byte offset " + str(index)
                )
            in_quotes = True
            at_field_start = False
            index += 1
            continue
        if byte == 44:
            at_field_start = True
            index += 1
            continue
        if newline is not None:
            _terminator, length = newline
            records += 1
            index += length
            last_record_end = index
            at_field_start = True
            continue
        at_field_start = False
        index += 1

    if in_quotes:
        raise CsvFormatError("unterminated quoted CSV field at end of file")
    if last_record_end < len(data):
        records += 1
    if records == 0:
        raise CsvFormatError("CSV file is empty and has no header")
    return records - 1


def classify_offline(v1_current: str, legacy: str) -> VerificationDecision:
    if legacy not in LEGACY_CLASSES:
        raise RowConflict("UNKNOWN_LEGACY_CLASS:" + legacy)

    label_support = "true" if legacy in ("Strong", "Weak") else "false"
    cluster_support = "true" if legacy in ("Strong", "Subclone") else "false"

    if v1_current == "Strong":
        if legacy == "Strong":
            return VerificationDecision(
                "Strong_Bidirectional", "BIDIRECTIONAL", label_support, cluster_support
            )
        if legacy == "Subclone":
            return VerificationDecision(
                "ClusterFirstOnly", "CLUSTER_FIRST_ONLY", label_support, cluster_support
            )
        raise RowConflict("INCOMPATIBLE_STRONG_LEGACY_PAIR:" + legacy)

    evidence_path = V1_CURRENT_TO_PATH.get(v1_current)
    if evidence_path is None:
        raise RowConflict("UNKNOWN_V1_CURRENT_CLASS:" + v1_current)
    if legacy not in ("Weak", "Noise"):
        raise RowConflict(
            "INCOMPATIBLE_NONSTRONG_LEGACY_PAIR:" + v1_current + ":" + legacy
        )
    return VerificationDecision(v1_current, evidence_path, label_support, cluster_support)


def expected_loh_subtype(potential_loh: str, legacy: str) -> str:
    if potential_loh == "false":
        return "None"
    if potential_loh != "true":
        raise RowConflict("INVALID_POTENTIAL_LOH:" + potential_loh)
    if legacy not in LOH_BY_LEGACY:
        raise RowConflict("UNKNOWN_LEGACY_FOR_LOH:" + legacy)
    return LOH_BY_LEGACY[legacy]


def _parse_nonnegative_int(value: str, context: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise MigrationError(context + " is not an integer: " + repr(value)) from error
    if parsed < 0:
        raise MigrationError(context + " must be non-negative")
    return parsed


def read_manifest(path: Path) -> List[ManifestEntry]:
    if not path.is_file():
        raise MigrationError("manifest is not a regular file: " + str(path))
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise MigrationError("manifest has no header")
            missing = [name for name in MANIFEST_BASE_HEADERS if name not in reader.fieldnames]
            if missing:
                raise MigrationError("manifest missing required columns: " + ",".join(missing))
            entries: List[ManifestEntry] = []
            seen = set()
            for line_number, row in enumerate(reader, start=2):
                filename = (row.get("file") or "").strip()
                if not filename:
                    raise MigrationError("manifest line %d has an empty file" % line_number)
                if Path(filename).name != filename or filename in (".", ".."):
                    raise MigrationError("manifest file must be a basename: " + filename)
                if filename in seen:
                    raise MigrationError("manifest contains duplicate file: " + filename)
                seen.add(filename)
                rows = _parse_nonnegative_int(
                    (row.get("rows") or "").strip(), "manifest rows for " + filename
                )
                sha256 = (row.get("sha256") or "").strip().lower()
                if not SHA256_RE.fullmatch(sha256):
                    raise MigrationError("manifest sha256 is invalid for " + filename)
                expected_counts: Dict[str, int] = {}
                for name in MANIFEST_COUNT_HEADERS:
                    raw = row.get(name)
                    if raw is not None and raw.strip() != "":
                        expected_counts[name] = _parse_nonnegative_int(
                            raw.strip(), "manifest %s for %s" % (name, filename)
                        )
                entries.append(ManifestEntry(filename, rows, sha256, expected_counts))
    except UnicodeDecodeError as error:
        raise MigrationError("manifest is not valid UTF-8") from error
    except csv.Error as error:
        raise MigrationError("manifest is not valid TSV: " + str(error)) from error
    if not entries:
        raise MigrationError("manifest contains no data rows")
    return entries


def resolve_input_path(input_root: Path, filename: str) -> Path:
    root = input_root.resolve()
    path = (root / filename).resolve()
    if path.parent != root:
        raise MigrationError("manifest input escapes input root: " + filename)
    if not path.is_file():
        raise MigrationError("manifest input is not a regular file: " + str(path))
    return path


def preflight_inputs(entries: Sequence[ManifestEntry], input_root: Path) -> None:
    for entry in entries:
        path = resolve_input_path(input_root, entry.filename)
        data = path.read_bytes()
        actual_sha = sha256_bytes(data)
        if actual_sha != entry.sha256:
            raise MigrationError(
                "manifest SHA-256 mismatch for %s: expected=%s actual=%s"
                % (entry.filename, entry.sha256, actual_sha)
            )
        actual_rows = count_csv_data_rows(data)
        if actual_rows != entry.rows:
            raise MigrationError(
                "manifest row-count mismatch for %s: expected=%d actual=%d"
                % (entry.filename, entry.rows, actual_rows)
            )


def _conflict_row(
    filename: str,
    row_number: int,
    names: Sequence[str],
    values: Sequence[str],
    reason: str,
) -> Dict[str, str]:
    lookup = dict(zip(names, values))

    conflict = {
        "file": filename,
        "row_number": str(row_number),
        "RegionID": lookup.get("RegionID", ""),
        "Chr": lookup.get("Chr", ""),
        "Pos": lookup.get("Pos", ""),
        "Ref": lookup.get("Ref", ""),
        "Alt": lookup.get("Alt", ""),
        "VerificationClass": lookup.get("VerificationClass", ""),
        "VerificationClass_Legacy": lookup.get("VerificationClass_Legacy", ""),
        "Potential_LOH": lookup.get("Potential_LOH", ""),
        "LOH_Subtype": lookup.get("LOH_Subtype", ""),
        "reason": reason,
    }
    return conflict


def _validate_headers(names: Sequence[str]) -> Tuple[Optional[str], Optional[str]]:
    duplicate = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicate:
        return "UNRESOLVED_SCHEMA", "DUPLICATE_HEADERS:" + ",".join(duplicate)
    already_present = sorted(name for name in APPENDED_HEADERS if name in names)
    if already_present:
        return "UNRESOLVED_SCHEMA", "OUTPUT_FIELDS_ALREADY_PRESENT:" + ",".join(already_present)
    if "VerificationClass_Legacy" not in names:
        return "UNRESOLVED_LEGACY_PROVENANCE", "MISSING_VERIFICATION_CLASS_LEGACY"
    missing_loh = [name for name in ("Potential_LOH", "LOH_Subtype") if name not in names]
    if missing_loh:
        return "UNRESOLVED_LOH_PROVENANCE", "MISSING_LOH_FIELDS:" + ",".join(missing_loh)
    required = list(REQUIRED_IDENTITY_HEADERS) + ["VerificationClass", "Significant"]
    missing = [name for name in required if name not in names]
    if missing:
        return "UNRESOLVED_SCHEMA", "MISSING_REQUIRED_FIELDS:" + ",".join(missing)
    return None, None


def validate_input_file(entry: ManifestEntry, path: Path) -> ValidationResult:
    data = path.read_bytes()
    if sha256_bytes(data) != entry.sha256:
        raise MigrationError("input changed after preflight: " + entry.filename)

    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise CsvFormatError("CSV input is not valid UTF-8: " + entry.filename) from error
    logical_rows = csv.reader(io.StringIO(text, newline=""), strict=True)
    try:
        names = next(logical_rows)
    except StopIteration as error:
        raise CsvFormatError("CSV has no header: " + entry.filename) from error
    except csv.Error as error:
        raise CsvFormatError("invalid CSV header in %s: %s" % (entry.filename, error)) from error
    if names and names[0].startswith("\ufeff"):
        names[0] = names[0].lstrip("\ufeff")
    status, reason = _validate_headers(names)
    if status is not None:
        return ValidationResult(
            status=status,
            reason=reason or "UNRESOLVED",
            input_rows=entry.rows,
        )

    indices = {name: index for index, name in enumerate(names)}
    result = ValidationResult(
        status="VALID",
        reason="READY",
        stable_key_uniqueness="PASS",
    )
    stable_keys = set()

    try:
        for row_number, values in enumerate(logical_rows, start=2):
            result.input_rows += 1
            if len(values) != len(names):
                result.stable_key_uniqueness = "FAIL"
                result.conflicts.append(
                    _conflict_row(
                        entry.filename,
                        row_number,
                        names,
                        values,
                        "FIELD_COUNT_MISMATCH:expected=%d:actual=%d"
                        % (len(names), len(values)),
                    )
                )
                continue

            key = tuple(values[indices[name]] for name in REQUIRED_IDENTITY_HEADERS)
            if any(value == "" for value in key):
                result.stable_key_uniqueness = "FAIL"
                result.conflicts.append(
                    _conflict_row(entry.filename, row_number, names, values, "EMPTY_STABLE_KEY")
                )
                continue
            if key in stable_keys:
                result.stable_key_uniqueness = "FAIL"
                result.conflicts.append(
                    _conflict_row(entry.filename, row_number, names, values, "DUPLICATE_STABLE_KEY")
                )
                continue
            stable_keys.add(key)

            current = values[indices["VerificationClass"]]
            legacy = values[indices["VerificationClass_Legacy"]]
            potential_loh = values[indices["Potential_LOH"]]
            loh_subtype = values[indices["LOH_Subtype"]]

            result.before_current[current] += 1
            result.before_legacy[legacy] += 1
            result.before_loh[loh_subtype] += 1
            if current == "Strong":
                result.observed_manifest_counts["input_final_strong"] += 1
            if current == "Strong" and legacy == "Strong":
                result.observed_manifest_counts["legacy_strong_to_final_strong"] += 1
            if current == "Strong" and legacy == "Subclone":
                result.observed_manifest_counts["legacy_subclone_to_final_strong"] += 1
            if (current == "Strong") != (legacy in ("Strong", "Subclone")):
                result.observed_manifest_counts["legacy_strong_or_subclone_exceptions"] += 1

            try:
                decision = classify_offline(current, legacy)
                if loh_subtype not in LOH_SUBTYPES:
                    raise RowConflict("UNKNOWN_LOH_SUBTYPE:" + loh_subtype)
                expected_loh = expected_loh_subtype(potential_loh, legacy)
                if loh_subtype != expected_loh:
                    raise RowConflict(
                        "LOH_PROVENANCE_MISMATCH:expected=%s:actual=%s"
                        % (expected_loh, loh_subtype)
                    )
            except RowConflict as error:
                result.conflicts.append(
                    _conflict_row(entry.filename, row_number, names, values, str(error))
                )
                continue

            result.after_current[decision.verification_class_v2] += 1
            result.after_evidence[decision.evidence_path] += 1
    except csv.Error as error:
        raise CsvFormatError("invalid CSV data in %s: %s" % (entry.filename, error)) from error

    if result.input_rows != entry.rows:
        raise MigrationError(
            "row count changed after preflight for %s: expected=%d actual=%d"
            % (entry.filename, entry.rows, result.input_rows)
        )

    manifest_mismatches = []
    for name, expected in entry.expected_counts.items():
        actual = result.observed_manifest_counts[name]
        if actual != expected:
            manifest_mismatches.append("%s:expected=%d:actual=%d" % (name, expected, actual))
    if manifest_mismatches:
        result.status = "FAILED_MANIFEST_COUNTS"
        result.reason = ";".join(manifest_mismatches)
        return result

    if result.conflicts:
        result.status = "FAILED_CONFLICT"
        result.reason = "UNMAPPED_OR_INCOMPATIBLE_ROWS"
        return result
    result.reason = "MAPPABLE"
    return result


def _row_decision(record: RawCsvRecord, indices: Mapping[str, int]) -> VerificationDecision:
    current = decode_csv_token(record.fields[indices["VerificationClass"]])
    legacy = decode_csv_token(record.fields[indices["VerificationClass_Legacy"]])
    return classify_offline(current, legacy)


def _migrated_record(
    record: RawCsvRecord,
    indices: Mapping[str, int],
    is_header: bool,
) -> bytes:
    fields = list(record.fields)
    if is_header:
        fields.extend(encode_csv_token(name) for name in APPENDED_HEADERS)
        return b",".join(fields) + record.newline

    original_current = fields[indices["VerificationClass"]]
    original_loh = fields[indices["LOH_Subtype"]]
    decision = _row_decision(record, indices)
    force_quote = original_current.startswith(b'"')
    fields[indices["VerificationClass"]] = encode_csv_token(
        decision.verification_class_v2, force_quote=force_quote
    )
    fields.extend(
        (
            b"2",
            original_current,
            decision.label_first_support.encode("ascii"),
            decision.cluster_first_support.encode("ascii"),
            b"NA",
            b"NA",
            decision.evidence_path.encode("ascii"),
            b"LEGACY_CLASS",
            original_loh,
        )
    )
    return b",".join(fields) + record.newline


def audit_raw_preservation(input_data: bytes, output_data: bytes) -> int:
    input_records = iter(iter_csv_records(input_data))
    output_records = iter(iter_csv_records(output_data))
    row_count = 0
    input_header: Optional[RawCsvRecord] = None
    indices: Dict[str, int] = {}
    input_stable_keys = set()
    output_stable_keys = set()

    while True:
        try:
            input_record = next(input_records)
        except StopIteration:
            input_record = None
        try:
            output_record = next(output_records)
        except StopIteration:
            output_record = None
        if input_record is None or output_record is None:
            if input_record is not output_record:
                raise MigrationError("raw preservation audit found a record-count mismatch")
            break

        if input_header is None:
            input_header = input_record
            names = header_names(input_record)
            indices = {name: index for index, name in enumerate(names)}
            expected_headers = tuple(encode_csv_token(name) for name in APPENDED_HEADERS)
            if output_record.fields[: len(input_record.fields)] != input_record.fields:
                raise MigrationError("raw preservation audit found a modified existing header token")
            if output_record.fields[len(input_record.fields) :] != expected_headers:
                raise MigrationError("raw preservation audit found incorrect appended headers")
        else:
            row_count += 1
            if len(output_record.fields) != len(input_record.fields) + len(APPENDED_HEADERS):
                raise MigrationError("raw preservation audit found a field-count mismatch")
            verification_index = indices["VerificationClass"]
            input_key = tuple(
                decode_csv_token(input_record.fields[indices[name]])
                for name in REQUIRED_IDENTITY_HEADERS
            )
            output_key = tuple(
                decode_csv_token(output_record.fields[indices[name]])
                for name in REQUIRED_IDENTITY_HEADERS
            )
            if any(value == "" for value in input_key) or input_key != output_key:
                raise MigrationError("Significant stable-key audit found a missing or changed key")
            if input_key in input_stable_keys or output_key in output_stable_keys:
                raise MigrationError("stable-key uniqueness audit found a duplicate key")
            input_stable_keys.add(input_key)
            output_stable_keys.add(output_key)
            significant_index = indices["Significant"]
            if output_record.fields[significant_index] != input_record.fields[significant_index]:
                raise MigrationError("Significant stable-key invariant audit failed")
            for index, token in enumerate(input_record.fields):
                if index == verification_index:
                    continue
                if output_record.fields[index] != token:
                    raise MigrationError(
                        "raw preservation audit found a changed existing token at data row %d field %d"
                        % (row_count + 1, index + 1)
                    )
            decision = _row_decision(input_record, indices)
            if decode_csv_token(output_record.fields[verification_index]) != decision.verification_class_v2:
                raise MigrationError("raw preservation audit found an incorrect v2 class")
            appended = output_record.fields[len(input_record.fields) :]
            expected_appended = {
                0: b"2",
                2: decision.label_first_support.encode("ascii"),
                3: decision.cluster_first_support.encode("ascii"),
                4: b"NA",
                5: b"NA",
                6: decision.evidence_path.encode("ascii"),
                7: b"LEGACY_CLASS",
            }
            for appended_index, expected_token in expected_appended.items():
                if appended[appended_index] != expected_token:
                    raise MigrationError(
                        "post-write audit found an incorrect appended token at data row %d field %s"
                        % (row_count + 1, APPENDED_HEADERS[appended_index])
                    )
            if appended[1] != input_record.fields[verification_index]:
                raise MigrationError("v1 deprecated class is not an exact raw-token copy")
            if appended[8] != input_record.fields[indices["LOH_Subtype"]]:
                raise MigrationError("LOH legacy alias is not an exact raw-token copy")
        if output_record.newline != input_record.newline:
            raise MigrationError("raw preservation audit found a changed record newline")
    return row_count


def write_migrated_file(input_path: Path, output_path: Path, entry: ManifestEntry) -> Tuple[int, str]:
    input_data = input_path.read_bytes()
    if sha256_bytes(input_data) != entry.sha256:
        raise MigrationError("input changed before write: " + entry.filename)
    records = iter_csv_records(input_data)
    try:
        header = next(records)
    except StopIteration as error:
        raise MigrationError("input has no header during write: " + entry.filename) from error
    names = header_names(header)
    indices = {name: index for index, name in enumerate(names)}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.parent / (
        "." + output_path.name + ".tmp." + secrets.token_hex(8)
    )
    digest = hashlib.sha256()
    rows = 0
    with temporary.open("xb") as handle:
        migrated = _migrated_record(header, indices, is_header=True)
        handle.write(migrated)
        digest.update(migrated)
        for record in records:
            migrated = _migrated_record(record, indices, is_header=False)
            handle.write(migrated)
            digest.update(migrated)
            rows += 1
        handle.flush()
        os.fsync(handle.fileno())

    output_data = temporary.read_bytes()
    audited_rows = audit_raw_preservation(input_data, output_data)
    if audited_rows != rows or rows != entry.rows:
        raise MigrationError(
            "post-write row audit failed for %s: manifest=%d written=%d audited=%d"
            % (entry.filename, entry.rows, rows, audited_rows)
        )
    if sha256_bytes(output_data) != digest.hexdigest():
        raise MigrationError("post-write SHA-256 audit failed for " + entry.filename)
    os.replace(temporary, output_path)
    return rows, digest.hexdigest()


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / ("." + path.name + ".tmp." + secrets.token_hex(8))
    with temporary.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _counter_dict(counter: Counter) -> Dict[str, int]:
    return {key: counter[key] for key in sorted(counter)}


def _json_cell(mapping: Mapping[str, int]) -> str:
    return json.dumps(dict(mapping), sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _tsv_bytes(headers: Sequence[str], rows: Iterable[Mapping[str, object]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=list(headers),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({name: row.get(name, "") for name in headers})
    return stream.getvalue().encode("utf-8")


def write_reports(
    staging: Path,
    reports: Sequence[FileReport],
    conflicts: Sequence[Mapping[str, str]],
    unresolved: Sequence[Mapping[str, str]],
    command: str,
    manifest_path: Path,
    manifest_sha256: str,
    input_root: Path,
    output_dir: Path,
    generated_at: str,
) -> Mapping[str, object]:
    file_headers = (
        "file",
        "status",
        "reason",
        "input_rows",
        "output_rows",
        "input_sha256",
        "output_sha256",
        "schema_version",
        "unmapped_count",
        "raw_token_preservation",
        "significant_stable_key_invariant",
        "stable_key_uniqueness",
        "before_current_counts",
        "before_legacy_counts",
        "before_loh_counts",
        "after_current_counts",
        "after_evidence_counts",
    )
    file_rows = []
    for report in reports:
        file_rows.append(
            {
                "file": report.filename,
                "status": report.status,
                "reason": report.reason,
                "input_rows": report.input_rows,
                "output_rows": report.output_rows,
                "input_sha256": report.input_sha256,
                "output_sha256": report.output_sha256,
                "schema_version": report.schema_version,
                "unmapped_count": report.unmapped_count,
                "raw_token_preservation": report.raw_token_preservation,
                "significant_stable_key_invariant": report.significant_stable_key_invariant,
                "stable_key_uniqueness": report.stable_key_uniqueness,
                "before_current_counts": _json_cell(report.before_current),
                "before_legacy_counts": _json_cell(report.before_legacy),
                "before_loh_counts": _json_cell(report.before_loh),
                "after_current_counts": _json_cell(report.after_current),
                "after_evidence_counts": _json_cell(report.after_evidence),
            }
        )
    atomic_write_bytes(staging / "migration_file_report.tsv", _tsv_bytes(file_headers, file_rows))

    conflict_headers = (
        "file",
        "row_number",
        "RegionID",
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "VerificationClass",
        "VerificationClass_Legacy",
        "Potential_LOH",
        "LOH_Subtype",
        "reason",
    )
    atomic_write_bytes(
        staging / "unmapped_conflicts.tsv", _tsv_bytes(conflict_headers, conflicts)
    )
    atomic_write_bytes(
        staging / "unresolved_files.tsv",
        _tsv_bytes(("file", "status", "reason"), unresolved),
    )

    output_manifest_rows = (
        {
            "file": report.filename,
            "rows": report.output_rows,
            "sha256": report.output_sha256,
            "schema_version": report.schema_version,
        }
        for report in reports
        if report.status == "VALID"
    )
    atomic_write_bytes(
        staging / "migrated_outputs_manifest.tsv",
        _tsv_bytes(("file", "rows", "sha256", "schema_version"), output_manifest_rows),
    )

    aggregate_before_current: Counter = Counter()
    aggregate_before_legacy: Counter = Counter()
    aggregate_before_loh: Counter = Counter()
    aggregate_after_current: Counter = Counter()
    aggregate_after_evidence: Counter = Counter()
    for report in reports:
        aggregate_before_current.update(report.before_current)
        aggregate_before_legacy.update(report.before_legacy)
        aggregate_before_loh.update(report.before_loh)
        if report.status == "VALID":
            aggregate_after_current.update(report.after_current)
            aggregate_after_evidence.update(report.after_evidence)

    failed = sum(report.status != "VALID" for report in reports)
    status = "VALID" if failed == 0 else "FAILED"
    reason = "ALL_FILES_MIGRATED" if failed == 0 else "FILE_MIGRATION_FAILURES"
    input_rows = sum(report.input_rows for report in reports)
    output_rows = sum(report.output_rows for report in reports)
    raw_token_preservation = (
        "PASS" if all(report.raw_token_preservation == "PASS" for report in reports) else "FAIL"
    )
    significant_stable_key_invariant = (
        "PASS"
        if all(report.significant_stable_key_invariant == "PASS" for report in reports)
        else "FAIL"
    )
    stable_key_uniqueness = (
        "PASS" if all(report.stable_key_uniqueness == "PASS" for report in reports) else "FAIL"
    )
    summary: Dict[str, object] = {
        "status": status,
        "reason": reason,
        "schema_version": SCHEMA_VERSION,
        "generated_at": generated_at,
        "command": command,
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": manifest_sha256,
        "input_root": str(input_root.resolve()),
        "output_dir": str(output_dir.resolve()),
        "total_files": len(reports),
        "valid_files": len(reports) - failed,
        "failed_files": failed,
        "input_rows": input_rows,
        "output_rows": output_rows,
        "unmapped_rows": len(conflicts),
        "raw_token_preservation": raw_token_preservation,
        "significant_stable_key_invariant": significant_stable_key_invariant,
        "stable_key_uniqueness": stable_key_uniqueness,
        "before_current_counts": _counter_dict(aggregate_before_current),
        "before_legacy_counts": _counter_dict(aggregate_before_legacy),
        "before_loh_counts": _counter_dict(aggregate_before_loh),
        "after_current_counts": _counter_dict(aggregate_after_current),
        "after_evidence_counts": _counter_dict(aggregate_after_evidence),
    }
    atomic_write_bytes(
        staging / "migration_summary.json",
        (json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    atomic_write_bytes(staging / "migration_command.txt", (command + "\n").encode("utf-8"))
    status_row = {
        "status": status,
        "reason": reason,
        "schema_version": SCHEMA_VERSION,
        "total_files": len(reports),
        "valid_files": len(reports) - failed,
        "failed_files": failed,
        "input_rows": input_rows,
        "output_rows": output_rows,
        "unmapped_rows": len(conflicts),
        "raw_token_preservation": raw_token_preservation,
        "significant_stable_key_invariant": significant_stable_key_invariant,
        "stable_key_uniqueness": stable_key_uniqueness,
        "generated_at": generated_at,
    }
    atomic_write_bytes(
        staging / "migration_status.tsv",
        _tsv_bytes(tuple(status_row.keys()), (status_row,)),
    )
    return summary


def run_migration(
    manifest_path: Path,
    input_root: Path,
    output_dir: Path,
    command: str,
) -> Mapping[str, object]:
    manifest_path = manifest_path.resolve()
    input_root = input_root.resolve()
    output_dir = output_dir.absolute()
    output_dir = output_dir.parent.resolve() / output_dir.name

    if not input_root.is_dir():
        raise MigrationError("input root is not a directory: " + str(input_root))
    if os.path.lexists(output_dir):
        raise MigrationError("output directory already exists: " + str(output_dir))
    if output_dir == input_root:
        raise MigrationError("output directory must not equal input root")
    try:
        output_dir.relative_to(input_root)
    except ValueError:
        pass
    else:
        raise MigrationError("output directory must not be inside input root")

    manifest_sha256 = sha256_file(manifest_path)
    entries = read_manifest(manifest_path)
    if sha256_file(manifest_path) != manifest_sha256:
        raise MigrationError("manifest changed while it was being read")
    preflight_inputs(entries, input_root)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_prefix = "." + output_dir.name + ".staging."
    staging = Path(tempfile.mkdtemp(prefix=staging_prefix, dir=str(output_dir.parent)))
    staging_identity = capture_staging_identity(staging, output_dir.parent, staging_prefix)

    try:
        generated_at = utc_now()
        reports: List[FileReport] = []
        all_conflicts: List[Mapping[str, str]] = []
        unresolved: List[Mapping[str, str]] = []

        for entry in entries:
            input_path = resolve_input_path(input_root, entry.filename)
            validation = validate_input_file(entry, input_path)
            all_conflicts.extend(validation.conflicts)
            output_rows = 0
            output_sha = ""
            raw_token_preservation = "NOT_RUN"
            significant_stable_key_invariant = "NOT_RUN"
            stable_key_uniqueness = validation.stable_key_uniqueness
            if validation.status == "VALID":
                output_rows, output_sha = write_migrated_file(
                    input_path, staging / entry.filename, entry
                )
                raw_token_preservation = "PASS"
                significant_stable_key_invariant = "PASS"
                stable_key_uniqueness = "PASS"
            else:
                unresolved.append(
                    {
                        "file": entry.filename,
                        "status": validation.status,
                        "reason": validation.reason,
                    }
                )
            reports.append(
                FileReport(
                    filename=entry.filename,
                    status=validation.status,
                    reason=validation.reason,
                    input_rows=validation.input_rows,
                    output_rows=output_rows,
                    input_sha256=entry.sha256,
                    output_sha256=output_sha,
                    schema_version=str(SCHEMA_VERSION) if validation.status == "VALID" else "",
                    unmapped_count=len(validation.conflicts),
                    before_current=_counter_dict(validation.before_current),
                    before_legacy=_counter_dict(validation.before_legacy),
                    before_loh=_counter_dict(validation.before_loh),
                    after_current=_counter_dict(validation.after_current),
                    after_evidence=_counter_dict(validation.after_evidence),
                    raw_token_preservation=raw_token_preservation,
                    significant_stable_key_invariant=significant_stable_key_invariant,
                    stable_key_uniqueness=stable_key_uniqueness,
                )
            )

        if sha256_file(manifest_path) != manifest_sha256:
            raise MigrationError("manifest changed during migration")
        summary = write_reports(
            staging,
            reports,
            all_conflicts,
            unresolved,
            command,
            manifest_path,
            manifest_sha256,
            input_root,
            output_dir,
            generated_at,
        )
        if sha256_file(manifest_path) != manifest_sha256:
            raise MigrationError("manifest changed before publication")
        if os.path.lexists(output_dir):
            raise MigrationError("output directory appeared during migration: " + str(output_dir))
        atomic_publish_noreplace(staging, output_dir)
    except BaseException as migration_error:
        try:
            cleanup_owned_staging(staging_identity)
        except BaseException as cleanup_error:
            raise MigrationError(
                "migration failed (%s); safe staging cleanup failed or was refused (%s)"
                % (migration_error, cleanup_error)
            ) from cleanup_error
        raise
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Byte-preserving offline migration to Verification schema v2."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--require-output-dir-absent",
        required=True,
        action="store_true",
        help="Required safety acknowledgement; migration never merges or overwrites output.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    raw_args = list(argv) if argv is not None else sys.argv[1:]
    command = shlex.join([sys.executable, str(Path(__file__).resolve()), *raw_args])
    try:
        summary = run_migration(args.manifest, args.input_root, args.output_dir, command)
    except (MigrationError, OSError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 2
    print(
        "migration status=%s valid_files=%s failed_files=%s input_rows=%s output_rows=%s output=%s"
        % (
            summary["status"],
            summary["valid_files"],
            summary["failed_files"],
            summary["input_rows"],
            summary["output_rows"],
            args.output_dir,
        )
    )
    return 0 if summary["status"] == "VALID" else 1


if __name__ == "__main__":
    raise SystemExit(main())
