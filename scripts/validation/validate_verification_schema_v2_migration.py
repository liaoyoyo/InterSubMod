#!/usr/bin/env python3
"""Independently validate a published Verification schema-v2 migration.

This program is deliberately read-only and does not import the migration
implementation.  It audits the frozen source corpus, the published CSVs, and
all seven provenance reports.  Every mismatch is fatal; no report is created
or repaired by this validator.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import mmap
import os
import re
import stat
import sys
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple


SCHEMA_VERSION = 2
PROVENANCE_REPORTS: Tuple[str, ...] = (
    "migration_summary.json",
    "migration_status.tsv",
    "migration_file_report.tsv",
    "migrated_outputs_manifest.tsv",
    "unmapped_conflicts.tsv",
    "unresolved_files.tsv",
    "migration_command.txt",
)
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
IDENTITY_HEADERS: Tuple[str, ...] = ("RegionID", "Chr", "Pos", "Ref", "Alt")
REQUIRED_SOURCE_HEADERS: Tuple[str, ...] = IDENTITY_HEADERS + (
    "VerificationClass",
    "VerificationClass_Legacy",
    "Potential_LOH",
    "LOH_Subtype",
    "Significant",
)
MANIFEST_HEADERS: Tuple[str, ...] = (
    "file",
    "rows",
    "sha256",
    "input_final_strong",
    "legacy_strong_to_final_strong",
    "legacy_subclone_to_final_strong",
    "legacy_strong_or_subclone_exceptions",
)
FILE_REPORT_HEADERS: Tuple[str, ...] = (
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
STATUS_HEADERS: Tuple[str, ...] = (
    "status",
    "reason",
    "schema_version",
    "total_files",
    "valid_files",
    "failed_files",
    "input_rows",
    "output_rows",
    "unmapped_rows",
    "raw_token_preservation",
    "significant_stable_key_invariant",
    "stable_key_uniqueness",
    "generated_at",
)
OUTPUT_MANIFEST_HEADERS: Tuple[str, ...] = ("file", "rows", "sha256", "schema_version")
CONFLICT_HEADERS: Tuple[str, ...] = (
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
UNRESOLVED_HEADERS: Tuple[str, ...] = ("file", "status", "reason")
PASS_FLAGS: Tuple[str, ...] = (
    "raw_token_preservation",
    "significant_stable_key_invariant",
    "stable_key_uniqueness",
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
UTC_SECONDS_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")

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


class ValidationError(RuntimeError):
    """The published migration does not satisfy the frozen contract."""


@dataclass(frozen=True)
class ExpectedCorpus:
    total_files: int = 14
    total_rows: int = 328_697
    input_final_strong: int = 66_841
    strong_bidirectional: int = 59_910
    cluster_first_only: int = 6_931
    exceptions: int = 0
    hcc1395_file: str = "significance_summary_HCC1395.csv"
    hcc1395_input_strong: int = 9_228
    hcc1395_cluster_first_only: int = 1_516


@dataclass(frozen=True)
class ManifestEntry:
    filename: str
    rows: int
    sha256: str
    input_final_strong: int
    legacy_strong_to_final_strong: int
    legacy_subclone_to_final_strong: int
    exceptions: int


@dataclass(frozen=True)
class RawRecord:
    fields: Tuple[bytes, ...]
    newline: bytes


@dataclass(frozen=True)
class Decision:
    current_v2: str
    evidence_path: str
    label_support: str
    cluster_support: str


@dataclass
class FileAudit:
    rows: int
    input_sha256: str
    output_sha256: str
    before_current: Counter
    before_legacy: Counter
    before_loh: Counter
    after_current: Counter
    after_evidence: Counter


@dataclass(frozen=True)
class ValidationSummary:
    total_files: int
    total_rows: int
    input_final_strong: int
    strong_bidirectional: int
    cluster_first_only: int
    unmapped_rows: int


def _fail(message: str) -> None:
    raise ValidationError(message)


def _parse_int(value: object, context: str) -> int:
    if isinstance(value, bool):
        _fail(context + " must be an integer, not boolean")
    try:
        parsed = int(str(value))
    except (TypeError, ValueError):
        _fail(context + " is not an integer: " + repr(value))
    if parsed < 0:
        _fail(context + " must be non-negative")
    return parsed


def _require_equal(actual: object, expected: object, context: str) -> None:
    if actual != expected:
        _fail("%s mismatch: expected=%r actual=%r" % (context, expected, actual))


def _regular_path(path: Path, context: str) -> Path:
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        _fail(context + " is missing: " + str(path))
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        _fail(context + " must be a non-symlink regular file: " + str(path))
    return path


@contextmanager
def _open_regular(path: Path, context: str) -> Iterator[BinaryIO]:
    _regular_path(path, context)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise ValidationError(context + " cannot be opened safely: " + str(error)) from error
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            _fail(context + " changed to a non-regular file while opening")
        with os.fdopen(descriptor, "rb", closefd=False) as handle:
            yield handle
    finally:
        os.close(descriptor)


def _read_regular_bytes(path: Path, context: str, max_bytes: int = 16 * 1024 * 1024) -> bytes:
    with _open_regular(path, context) as handle:
        metadata = os.fstat(handle.fileno())
        if metadata.st_size > max_bytes:
            _fail("%s exceeds bounded report size (%d bytes)" % (context, max_bytes))
        data = handle.read(max_bytes + 1)
    if len(data) > max_bytes:
        _fail("%s exceeds bounded report size (%d bytes)" % (context, max_bytes))
    return data


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _decode_token(token: bytes) -> str:
    if token.startswith(b'"'):
        if len(token) < 2 or not token.endswith(b'"'):
            _fail("quoted CSV token does not end with a quote")
        payload = token[1:-1].replace(b'""', b'"')
    else:
        if b'"' in token:
            _fail("unquoted CSV token contains a quote")
        payload = token
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValidationError("CSV token is not valid UTF-8") from error


def _encode_token(value: str, force_quote: bool = False) -> bytes:
    encoded = value.encode("utf-8")
    if not force_quote and not any(marker in encoded for marker in (b",", b'"', b"\r", b"\n")):
        return encoded
    return b'"' + encoded.replace(b'"', b'""') + b'"'


def _newline_at(data: mmap.mmap, index: int) -> Optional[Tuple[bytes, int]]:
    byte = data[index]
    if byte == 10:
        return b"\n", 1
    if byte == 13:
        if index + 1 < len(data) and data[index + 1] == 10:
            return b"\r\n", 2
        return b"\r", 1
    return None


def _iter_records(data: mmap.mmap) -> Iterator[RawRecord]:
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
                yield RawRecord(tuple(fields), terminator)
                fields = []
                index += length
                field_start = index
                closed_quote = False
                continue
            _fail("unexpected byte after closing CSV quote at offset %d" % index)

        if byte == 34:
            if index != field_start:
                _fail("quote appears inside unquoted CSV field at offset %d" % index)
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
            yield RawRecord(tuple(fields), terminator)
            fields = []
            index += length
            field_start = index
            continue
        index += 1

    if in_quotes:
        _fail("unterminated quoted CSV field")
    if field_start < len(data) or fields or closed_quote:
        fields.append(data[field_start:len(data)])
        yield RawRecord(tuple(fields), b"")


def _headers(record: RawRecord) -> List[str]:
    names = [_decode_token(token) for token in record.fields]
    if names and names[0].startswith("\ufeff"):
        names[0] = names[0].lstrip("\ufeff")
    duplicate = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicate:
        _fail("duplicate CSV headers: " + ",".join(duplicate))
    return names


def _classify(current: str, legacy: str) -> Decision:
    label = "true" if legacy in ("Strong", "Weak") else "false"
    cluster = "true" if legacy in ("Strong", "Subclone") else "false"
    if legacy not in ("Strong", "Subclone", "Weak", "Noise"):
        _fail("unknown legacy class: " + legacy)
    if current == "Strong":
        if legacy == "Strong":
            return Decision("Strong_Bidirectional", "BIDIRECTIONAL", label, cluster)
        if legacy == "Subclone":
            return Decision("ClusterFirstOnly", "CLUSTER_FIRST_ONLY", label, cluster)
        _fail("incompatible Strong/legacy pair: " + legacy)
    evidence = V1_CURRENT_TO_PATH.get(current)
    if evidence is None:
        _fail("unknown v1 current class: " + current)
    if legacy not in ("Weak", "Noise"):
        _fail("incompatible non-Strong/legacy pair: %s/%s" % (current, legacy))
    return Decision(current, evidence, label, cluster)


def _read_tsv(path: Path, expected_headers: Sequence[str], context: str) -> List[Dict[str, str]]:
    data = _read_regular_bytes(path, context)
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValidationError(context + " is not UTF-8") from error
    try:
        reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
        _require_equal(tuple(reader.fieldnames or ()), tuple(expected_headers), context + " headers")
        rows = list(reader)
    except csv.Error as error:
        raise ValidationError(context + " is malformed TSV: " + str(error)) from error
    if any(None in row for row in rows):
        _fail(context + " contains rows with extra fields")
    return rows


def _read_json_object(path: Path, context: str) -> Dict[str, object]:
    data = _read_regular_bytes(path, context)
    try:
        parsed = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValidationError(context + " is not valid UTF-8 JSON: " + str(error)) from error
    if not isinstance(parsed, dict):
        _fail(context + " must contain one JSON object")
    return parsed


def _index_unique(rows: Iterable[Mapping[str, str]], field: str, context: str) -> Dict[str, Mapping[str, str]]:
    indexed: Dict[str, Mapping[str, str]] = {}
    for row in rows:
        key = row.get(field, "")
        if not key or Path(key).name != key or key in (".", ".."):
            _fail(context + " contains an invalid basename: " + repr(key))
        if key in indexed:
            _fail(context + " contains duplicate " + key)
        indexed[key] = row
    return indexed


def _read_manifest(path: Path) -> Tuple[List[ManifestEntry], str]:
    data = _read_regular_bytes(path, "frozen manifest")
    manifest_sha = _sha256_bytes(data)
    try:
        reader = csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""), delimiter="\t", strict=True)
        _require_equal(tuple(reader.fieldnames or ()), MANIFEST_HEADERS, "frozen manifest headers")
        raw_rows = list(reader)
    except (UnicodeDecodeError, csv.Error) as error:
        raise ValidationError("frozen manifest is malformed: " + str(error)) from error
    if not raw_rows:
        _fail("frozen manifest has no entries")
    entries: List[ManifestEntry] = []
    seen = set()
    for row in raw_rows:
        if None in row:
            _fail("frozen manifest contains an extra field")
        filename = row["file"]
        if Path(filename).name != filename or filename in ("", ".", "..") or filename in seen:
            _fail("frozen manifest contains an invalid or duplicate basename: " + repr(filename))
        seen.add(filename)
        sha256 = row["sha256"]
        if not SHA256_RE.fullmatch(sha256):
            _fail("frozen manifest has an invalid SHA-256 for " + filename)
        entries.append(
            ManifestEntry(
                filename=filename,
                rows=_parse_int(row["rows"], "manifest rows for " + filename),
                sha256=sha256,
                input_final_strong=_parse_int(row["input_final_strong"], "input_final_strong for " + filename),
                legacy_strong_to_final_strong=_parse_int(
                    row["legacy_strong_to_final_strong"], "legacy strong count for " + filename
                ),
                legacy_subclone_to_final_strong=_parse_int(
                    row["legacy_subclone_to_final_strong"], "legacy subclone count for " + filename
                ),
                exceptions=_parse_int(
                    row["legacy_strong_or_subclone_exceptions"], "exception count for " + filename
                ),
            )
        )
    return entries, manifest_sha


def _parse_counter(value: str, context: str) -> Counter:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise ValidationError(context + " is not valid JSON") from error
    if not isinstance(parsed, dict):
        _fail(context + " must be a JSON object")
    result: Counter = Counter()
    for key, raw_count in parsed.items():
        if not isinstance(key, str):
            _fail(context + " has a non-string key")
        result[key] = _parse_int(raw_count, context + "[" + key + "]")
    return result


def _audit_csv_pair(input_path: Path, output_path: Path, entry: ManifestEntry) -> FileAudit:
    with _open_regular(input_path, "frozen input " + entry.filename) as input_handle:
        with _open_regular(output_path, "migrated output " + entry.filename) as output_handle:
            if os.fstat(input_handle.fileno()).st_size == 0 or os.fstat(output_handle.fileno()).st_size == 0:
                _fail("input/output CSV must not be empty for " + entry.filename)
            with mmap.mmap(input_handle.fileno(), 0, access=mmap.ACCESS_READ) as input_data:
                with mmap.mmap(output_handle.fileno(), 0, access=mmap.ACCESS_READ) as output_data:
                    input_sha = hashlib.sha256(input_data).hexdigest()
                    output_sha = hashlib.sha256(output_data).hexdigest()
                    _require_equal(input_sha, entry.sha256, "frozen input SHA-256 for " + entry.filename)

                    input_records = iter(_iter_records(input_data))
                    output_records = iter(_iter_records(output_data))
                    try:
                        input_header = next(input_records)
                        output_header = next(output_records)
                    except StopIteration:
                        _fail("input/output CSV lacks a header for " + entry.filename)
                    input_names = _headers(input_header)
                    missing = [name for name in REQUIRED_SOURCE_HEADERS if name not in input_names]
                    if missing:
                        _fail("frozen input %s lacks headers: %s" % (entry.filename, ",".join(missing)))
                    indices = {name: index for index, name in enumerate(input_names)}
                    expected_appended_headers = tuple(_encode_token(name) for name in APPENDED_HEADERS)
                    _require_equal(
                        output_header.fields[: len(input_header.fields)],
                        input_header.fields,
                        "existing header tokens for " + entry.filename,
                    )
                    _require_equal(
                        output_header.fields[len(input_header.fields) :],
                        expected_appended_headers,
                        "append-only headers for " + entry.filename,
                    )
                    _require_equal(output_header.newline, input_header.newline, "header newline for " + entry.filename)

                    rows = 0
                    input_keys = set()
                    output_keys = set()
                    before_current: Counter = Counter()
                    before_legacy: Counter = Counter()
                    before_loh: Counter = Counter()
                    after_current: Counter = Counter()
                    after_evidence: Counter = Counter()
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
                                _fail("record count differs for " + entry.filename)
                            break
                        rows += 1
                        if len(input_record.fields) != len(input_header.fields):
                            _fail("input field count differs at %s row %d" % (entry.filename, rows + 1))
                        if len(output_record.fields) != len(input_record.fields) + len(APPENDED_HEADERS):
                            _fail("output field count differs at %s row %d" % (entry.filename, rows + 1))

                        input_key = tuple(_decode_token(input_record.fields[indices[name]]) for name in IDENTITY_HEADERS)
                        output_key = tuple(_decode_token(output_record.fields[indices[name]]) for name in IDENTITY_HEADERS)
                        if any(value == "" for value in input_key + output_key):
                            _fail("empty stable key at %s row %d" % (entry.filename, rows + 1))
                        if input_key in input_keys or output_key in output_keys:
                            _fail("duplicate stable key at %s row %d" % (entry.filename, rows + 1))
                        input_keys.add(input_key)
                        output_keys.add(output_key)
                        _require_equal(output_key, input_key, "stable-key/order at %s row %d" % (entry.filename, rows + 1))

                        verification_index = indices["VerificationClass"]
                        significant_index = indices["Significant"]
                        _require_equal(
                            output_record.fields[significant_index],
                            input_record.fields[significant_index],
                            "Significant raw token for key " + repr(input_key),
                        )
                        for index, token in enumerate(input_record.fields):
                            if index != verification_index:
                                _require_equal(
                                    output_record.fields[index],
                                    token,
                                    "existing raw token at %s row %d field %d" % (entry.filename, rows + 1, index + 1),
                                )
                        _require_equal(output_record.newline, input_record.newline, "record newline at %s row %d" % (entry.filename, rows + 1))

                        current = _decode_token(input_record.fields[verification_index])
                        legacy = _decode_token(input_record.fields[indices["VerificationClass_Legacy"]])
                        loh_subtype = _decode_token(input_record.fields[indices["LOH_Subtype"]])
                        decision = _classify(current, legacy)
                        expected_current_token = _encode_token(
                            decision.current_v2, force_quote=input_record.fields[verification_index].startswith(b'"')
                        )
                        _require_equal(
                            output_record.fields[verification_index],
                            expected_current_token,
                            "VerificationClass v2 token at %s row %d" % (entry.filename, rows + 1),
                        )
                        appended = output_record.fields[len(input_record.fields) :]
                        expected_appended = (
                            b"2",
                            input_record.fields[verification_index],
                            decision.label_support.encode("ascii"),
                            decision.cluster_support.encode("ascii"),
                            b"NA",
                            b"NA",
                            decision.evidence_path.encode("ascii"),
                            b"LEGACY_CLASS",
                            input_record.fields[indices["LOH_Subtype"]],
                        )
                        _require_equal(
                            appended,
                            expected_appended,
                            "appended provenance contract at %s row %d" % (entry.filename, rows + 1),
                        )
                        before_current[current] += 1
                        before_legacy[legacy] += 1
                        before_loh[loh_subtype] += 1
                        after_current[decision.current_v2] += 1
                        after_evidence[decision.evidence_path] += 1

                    _require_equal(rows, entry.rows, "CSV row count for " + entry.filename)
                    _require_equal(input_keys, output_keys, "stable-key set for " + entry.filename)
                    _require_equal(hashlib.sha256(input_data).hexdigest(), input_sha, "input stability for " + entry.filename)
                    _require_equal(hashlib.sha256(output_data).hexdigest(), output_sha, "output stability for " + entry.filename)
    return FileAudit(
        rows=rows,
        input_sha256=input_sha,
        output_sha256=output_sha,
        before_current=before_current,
        before_legacy=before_legacy,
        before_loh=before_loh,
        after_current=after_current,
        after_evidence=after_evidence,
    )


def validate_migration(
    manifest_path: Path,
    input_root: Path,
    output_dir: Path,
    expected: ExpectedCorpus = ExpectedCorpus(),
) -> ValidationSummary:
    """Validate one immutable publication without writing any files."""

    manifest_path = manifest_path.resolve(strict=True)
    input_root = input_root.resolve(strict=True)
    output_dir = output_dir.resolve(strict=True)
    if not input_root.is_dir() or input_root.is_symlink():
        _fail("input root must resolve to a real directory: " + str(input_root))
    if not output_dir.is_dir() or output_dir.is_symlink():
        _fail("output directory must resolve to a real directory: " + str(output_dir))
    if output_dir == input_root:
        _fail("output directory must differ from input root")

    entries, manifest_sha = _read_manifest(manifest_path)
    _require_equal(len(entries), expected.total_files, "frozen manifest file count")
    _require_equal(sum(entry.rows for entry in entries), expected.total_rows, "frozen manifest row total")
    _require_equal(
        sum(entry.input_final_strong for entry in entries), expected.input_final_strong, "manifest input Strong total"
    )
    _require_equal(
        sum(entry.legacy_strong_to_final_strong for entry in entries),
        expected.strong_bidirectional,
        "manifest legacy Strong-to-Strong total",
    )
    _require_equal(
        sum(entry.legacy_subclone_to_final_strong for entry in entries),
        expected.cluster_first_only,
        "manifest legacy Subclone-to-Strong total",
    )
    _require_equal(sum(entry.exceptions for entry in entries), expected.exceptions, "manifest exception total")
    entry_by_name = {entry.filename: entry for entry in entries}
    if expected.hcc1395_file not in entry_by_name:
        _fail("authoritative HCC1395 file is missing from manifest: " + expected.hcc1395_file)
    hcc_entry = entry_by_name[expected.hcc1395_file]
    _require_equal(hcc_entry.input_final_strong, expected.hcc1395_input_strong, "HCC1395 manifest Strong")
    _require_equal(
        hcc_entry.legacy_subclone_to_final_strong,
        expected.hcc1395_cluster_first_only,
        "HCC1395 manifest ClusterFirstOnly",
    )

    report_bytes = {
        name: _read_regular_bytes(output_dir / name, "provenance report " + name)
        for name in PROVENANCE_REPORTS
    }
    report_hashes = {name: _sha256_bytes(data) for name, data in report_bytes.items()}
    summary = _read_json_object(output_dir / "migration_summary.json", "migration summary")
    status_rows = _read_tsv(output_dir / "migration_status.tsv", STATUS_HEADERS, "migration status")
    file_rows = _read_tsv(output_dir / "migration_file_report.tsv", FILE_REPORT_HEADERS, "file report")
    output_manifest_rows = _read_tsv(
        output_dir / "migrated_outputs_manifest.tsv", OUTPUT_MANIFEST_HEADERS, "output manifest"
    )
    conflict_rows = _read_tsv(output_dir / "unmapped_conflicts.tsv", CONFLICT_HEADERS, "unmapped conflicts")
    unresolved_rows = _read_tsv(output_dir / "unresolved_files.tsv", UNRESOLVED_HEADERS, "unresolved files")
    _require_equal(len(status_rows), 1, "migration status row count")
    _require_equal(len(conflict_rows), 0, "unmapped conflict row count")
    _require_equal(len(unresolved_rows), 0, "unresolved file row count")

    status_row = status_rows[0]
    file_by_name = _index_unique(file_rows, "file", "file report")
    output_manifest_by_name = _index_unique(output_manifest_rows, "file", "output manifest")
    expected_names = set(entry_by_name)
    _require_equal(set(file_by_name), expected_names, "file-report filenames")
    _require_equal(set(output_manifest_by_name), expected_names, "output-manifest filenames")

    expected_directory_entries = expected_names | set(PROVENANCE_REPORTS)
    actual_directory_entries = {item.name for item in os.scandir(output_dir)}
    _require_equal(actual_directory_entries, expected_directory_entries, "published output directory entries")

    aggregate_before_current: Counter = Counter()
    aggregate_before_legacy: Counter = Counter()
    aggregate_before_loh: Counter = Counter()
    aggregate_after_current: Counter = Counter()
    aggregate_after_evidence: Counter = Counter()
    audits: Dict[str, FileAudit] = {}
    for entry in entries:
        input_path = input_root / entry.filename
        output_path = output_dir / entry.filename
        if input_path.resolve().parent != input_root or output_path.resolve().parent != output_dir:
            _fail("bounded-path check failed for " + entry.filename)
        report = file_by_name[entry.filename]
        output_manifest = output_manifest_by_name[entry.filename]
        audit = _audit_csv_pair(input_path, output_path, entry)
        audits[entry.filename] = audit
        _require_equal(output_manifest["sha256"], audit.output_sha256, "output-manifest SHA for " + entry.filename)
        _require_equal(_parse_int(output_manifest["rows"], "output rows"), audit.rows, "output-manifest rows for " + entry.filename)
        _require_equal(output_manifest["schema_version"], str(SCHEMA_VERSION), "output schema for " + entry.filename)
        _require_equal(report["status"], "VALID", "file status for " + entry.filename)
        _require_equal(report["reason"], "MAPPABLE", "file reason for " + entry.filename)
        _require_equal(_parse_int(report["input_rows"], "file input rows"), audit.rows, "file input rows for " + entry.filename)
        _require_equal(_parse_int(report["output_rows"], "file output rows"), audit.rows, "file output rows for " + entry.filename)
        _require_equal(report["input_sha256"], audit.input_sha256, "file input SHA for " + entry.filename)
        _require_equal(report["output_sha256"], audit.output_sha256, "file output SHA for " + entry.filename)
        _require_equal(report["schema_version"], str(SCHEMA_VERSION), "file schema for " + entry.filename)
        _require_equal(_parse_int(report["unmapped_count"], "file unmapped count"), 0, "file unmapped count for " + entry.filename)
        for flag in PASS_FLAGS:
            _require_equal(report[flag], "PASS", "%s for %s" % (flag, entry.filename))
        counter_contracts = (
            ("before_current_counts", audit.before_current),
            ("before_legacy_counts", audit.before_legacy),
            ("before_loh_counts", audit.before_loh),
            ("after_current_counts", audit.after_current),
            ("after_evidence_counts", audit.after_evidence),
        )
        for field, actual_counter in counter_contracts:
            _require_equal(_parse_counter(report[field], "%s/%s" % (entry.filename, field)), actual_counter, "%s for %s" % (field, entry.filename))
        aggregate_before_current.update(audit.before_current)
        aggregate_before_legacy.update(audit.before_legacy)
        aggregate_before_loh.update(audit.before_loh)
        aggregate_after_current.update(audit.after_current)
        aggregate_after_evidence.update(audit.after_evidence)

    _require_equal(aggregate_before_current["Strong"], expected.input_final_strong, "independent input Strong total")
    _require_equal(aggregate_after_current["Strong_Bidirectional"], expected.strong_bidirectional, "independent Strong_Bidirectional total")
    _require_equal(aggregate_after_current["ClusterFirstOnly"], expected.cluster_first_only, "independent ClusterFirstOnly total")
    _require_equal(
        aggregate_after_current["Strong_Bidirectional"] + aggregate_after_current["ClusterFirstOnly"],
        aggregate_before_current["Strong"],
        "Strong split conservation",
    )
    hcc_audit = audits[expected.hcc1395_file]
    _require_equal(hcc_audit.before_current["Strong"], expected.hcc1395_input_strong, "HCC1395 input Strong")
    _require_equal(
        hcc_audit.after_current["ClusterFirstOnly"],
        expected.hcc1395_cluster_first_only,
        "HCC1395 ClusterFirstOnly",
    )

    _require_equal(summary.get("status"), "VALID", "summary status")
    _require_equal(summary.get("reason"), "ALL_FILES_MIGRATED", "summary reason")
    _require_equal(_parse_int(summary.get("schema_version"), "summary schema version"), SCHEMA_VERSION, "summary schema version")
    _require_equal(summary.get("manifest"), str(manifest_path), "summary manifest path")
    _require_equal(summary.get("manifest_sha256"), manifest_sha, "summary manifest SHA")
    _require_equal(summary.get("input_root"), str(input_root), "summary input root")
    _require_equal(summary.get("output_dir"), str(output_dir), "summary output directory")
    _require_equal(_parse_int(summary.get("total_files"), "summary total files"), expected.total_files, "summary total files")
    _require_equal(_parse_int(summary.get("valid_files"), "summary valid files"), expected.total_files, "summary valid files")
    _require_equal(_parse_int(summary.get("failed_files"), "summary failed files"), 0, "summary failed files")
    _require_equal(_parse_int(summary.get("input_rows"), "summary input rows"), expected.total_rows, "summary input rows")
    _require_equal(_parse_int(summary.get("output_rows"), "summary output rows"), expected.total_rows, "summary output rows")
    _require_equal(_parse_int(summary.get("unmapped_rows"), "summary unmapped rows"), expected.exceptions, "summary unmapped rows")
    for flag in PASS_FLAGS:
        _require_equal(summary.get(flag), "PASS", "summary " + flag)
    summary_counters = (
        ("before_current_counts", aggregate_before_current),
        ("before_legacy_counts", aggregate_before_legacy),
        ("before_loh_counts", aggregate_before_loh),
        ("after_current_counts", aggregate_after_current),
        ("after_evidence_counts", aggregate_after_evidence),
    )
    for field, counter in summary_counters:
        value = summary.get(field)
        if not isinstance(value, dict):
            _fail("summary %s must be an object" % field)
        parsed_counter = Counter({str(key): _parse_int(count, "summary %s/%s" % (field, key)) for key, count in value.items()})
        _require_equal(parsed_counter, counter, "summary " + field)

    _require_equal(status_row["status"], "VALID", "status TSV status")
    _require_equal(status_row["reason"], "ALL_FILES_MIGRATED", "status TSV reason")
    status_integer_contracts = {
        "schema_version": SCHEMA_VERSION,
        "total_files": expected.total_files,
        "valid_files": expected.total_files,
        "failed_files": 0,
        "input_rows": expected.total_rows,
        "output_rows": expected.total_rows,
        "unmapped_rows": expected.exceptions,
    }
    for field, expected_value in status_integer_contracts.items():
        _require_equal(_parse_int(status_row[field], "status " + field), expected_value, "status " + field)
    for flag in PASS_FLAGS:
        _require_equal(status_row[flag], "PASS", "status " + flag)
    generated_at = summary.get("generated_at")
    if not isinstance(generated_at, str) or not UTC_SECONDS_RE.fullmatch(generated_at):
        _fail("summary generated_at is not UTC second precision")
    _require_equal(status_row["generated_at"], generated_at, "status generated_at")
    command = summary.get("command")
    if not isinstance(command, str) or not command.strip():
        _fail("summary command is empty")
    _require_equal(report_bytes["migration_command.txt"], (command + "\n").encode("utf-8"), "migration command report")

    _require_equal(_sha256_bytes(_read_regular_bytes(manifest_path, "frozen manifest final read")), manifest_sha, "manifest stability")
    final_entries = {item.name for item in os.scandir(output_dir)}
    _require_equal(final_entries, expected_directory_entries, "final published output directory entries")
    for name, expected_hash in report_hashes.items():
        _require_equal(
            _sha256_bytes(_read_regular_bytes(output_dir / name, "final provenance report " + name)),
            expected_hash,
            "provenance report stability for " + name,
        )

    return ValidationSummary(
        total_files=expected.total_files,
        total_rows=expected.total_rows,
        input_final_strong=expected.input_final_strong,
        strong_bidirectional=expected.strong_bidirectional,
        cluster_first_only=expected.cluster_first_only,
        unmapped_rows=expected.exceptions,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read-only fail-closed validator for the frozen Verification schema-v2 migration."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = validate_migration(args.manifest, args.input_root, args.output_dir)
    except (ValidationError, OSError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 2
    print(
        "validator status=VALID files=%d rows=%d input_strong=%d "
        "strong_bidirectional=%d cluster_first_only=%d unmapped=%d"
        % (
            result.total_files,
            result.total_rows,
            result.input_final_strong,
            result.strong_bidirectional,
            result.cluster_first_only,
            result.unmapped_rows,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
