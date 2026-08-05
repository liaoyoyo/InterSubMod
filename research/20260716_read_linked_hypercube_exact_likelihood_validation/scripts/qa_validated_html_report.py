#!/usr/bin/env python3
"""Browser QA and immutable verification for the professor-facing HTML report.

The historical invocation remains supported::

    qa_validated_html_report.py --html REPORT --output-dir QA --expect-status final

After QA, the exact artifact can be revalidated without launching a browser::

    qa_validated_html_report.py verify-only --html REPORT --output-dir QA

The trust model is checksum provenance under a non-hostile same-UID account.  A
successful QA directory is an exact five-file, read-only artifact.  Failed
work is preserved under a visibly named diagnostic directory and is never
deleted or authenticated as a successful QA receipt.
"""

from __future__ import annotations

import argparse
import ctypes
import datetime as dt
import errno
import hashlib
import json
import os
import stat
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import unquote, urlparse

from playwright.sync_api import sync_playwright


SCHEMA_NAME = "intersubmod.validated_html_browser_qa"
SCHEMA_VERSION = "2.0.0"
ARTIFACT_CLASS = "IMMUTABLE_BROWSER_QA_PASS"
RECEIPT_NAME = "browser_qa_receipt.json"
SIDECAR_NAME = f"{RECEIPT_NAME}.sha256"
OUTPUT_FILENAMES = {
    "desktop_screenshot": "desktop_full.png",
    "mobile_screenshot": "mobile_full.png",
    "print_pdf": "print_qa.pdf",
}
EXACT_QA_FILES = tuple((*OUTPUT_FILENAMES.values(), RECEIPT_NAME, SIDECAR_NAME))
IDENTITY_KEYS = {"path", "size_bytes", "sha256", "mode_octal", "st_nlink"}
CHECK_NAMES = (
    "title_present",
    "lang_is_zh_hant",
    "one_h1",
    "minimum_sections",
    "minimum_svg_figures",
    "collapsible_definitions_present",
    "details_toggle_works",
    "unique_ids",
    "all_internal_anchors_resolve",
    "all_local_link_targets_exist",
    "no_external_links",
    "no_remote_requests",
    "no_failed_requests",
    "no_console_errors",
    "no_page_errors",
    "desktop_no_body_overflow",
    "mobile_no_body_overflow",
    "status_matches_expected",
    "html_unchanged_during_qa",
    "qa_producer_unchanged_during_qa",
)
METRIC_NAMES = (
    "title",
    "status_text",
    "h1_count",
    "section_count",
    "svg_count",
    "details_count",
    "desktop_widths",
    "mobile_widths",
    "missing_anchor_targets",
    "missing_local_files",
    "external_links",
    "remote_requests",
    "failed_requests",
    "console_errors",
    "page_errors",
)
ASSURANCE = {
    "scientific_evidence_authority": False,
    "content_binding": "absolute path, byte size, and full SHA-256 for HTML, QA producer, and three rendered outputs",
    "publication": "same-filesystem renameat2(RENAME_NOREPLACE); receipt and sidecar created with O_EXCL",
    "filesystem": "HTML and five QA files mode 0444 with st_nlink=1; exact QA directory mode 0555",
    "threat_model": "checksum provenance; no hostile same-UID actor",
}


class QaArtifactError(RuntimeError):
    """Fail-closed QA artifact construction or verification error."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise QaArtifactError(message)


def _absolute(path: Path) -> Path:
    """Return an absolute lexical path without silently accepting a symlink."""

    return Path(os.path.abspath(os.fspath(path)))


def _stable_read(path: Path, label: str) -> tuple[bytes, os.stat_result]:
    """Read one non-aliased regular file while detecting concurrent changes."""

    absolute = _absolute(path)
    try:
        before_path = os.lstat(absolute)
    except OSError as exc:
        raise QaArtifactError(f"{label} is unavailable: {absolute}: {exc}") from exc
    require(stat.S_ISREG(before_path.st_mode), f"{label} is not a regular file: {absolute}")
    require(not stat.S_ISLNK(before_path.st_mode), f"{label} is a symlink: {absolute}")
    require(before_path.st_nlink == 1, f"{label} hardlink count is not one: {absolute}")

    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(absolute, flags)
    except OSError as exc:
        raise QaArtifactError(f"cannot open {label}: {absolute}: {exc}") from exc
    try:
        before_fd = os.fstat(descriptor)
        require(stat.S_ISREG(before_fd.st_mode), f"opened {label} is not regular: {absolute}")
        require(before_fd.st_nlink == 1, f"opened {label} hardlink count is not one: {absolute}")
        require(
            (before_path.st_dev, before_path.st_ino) == (before_fd.st_dev, before_fd.st_ino),
            f"{label} path changed before read: {absolute}",
        )
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        after_fd = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = os.lstat(absolute)
    stable_fields = ("st_dev", "st_ino", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
    require(
        all(getattr(before_fd, key) == getattr(after_fd, key) for key in stable_fields)
        and all(getattr(after_fd, key) == getattr(after_path, key) for key in stable_fields),
        f"{label} changed during read: {absolute}",
    )
    require(after_fd.st_nlink == 1, f"{label} gained a hardlink during read: {absolute}")
    payload = b"".join(chunks)
    require(len(payload) == after_fd.st_size, f"{label} byte count changed during read: {absolute}")
    return payload, after_fd


def artifact_identity(path: Path, *, recorded_path: Path | None = None) -> dict[str, object]:
    """Return an exact path/content/mode/link identity for a regular file."""

    payload, observed = _stable_read(path, "artifact")
    logical = _absolute(path if recorded_path is None else recorded_path)
    return {
        "path": str(logical),
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "mode_octal": f"{stat.S_IMODE(observed.st_mode):04o}",
        "st_nlink": int(observed.st_nlink),
    }


def file_identity(path: Path) -> dict[str, object]:
    """Compatibility helper returning the original path/size/SHA identity."""

    identity = artifact_identity(path)
    return {key: identity[key] for key in ("path", "size_bytes", "sha256")}


def file_identity_matches(path: Path, expected: Mapping[str, object]) -> bool:
    """Fail closed when legacy identity fields are absent, malformed, or stale."""

    if set(expected) != {"path", "size_bytes", "sha256"}:
        return False
    recorded_path = expected.get("path")
    size_bytes = expected.get("size_bytes")
    digest = expected.get("sha256")
    if not isinstance(recorded_path, str) or not Path(recorded_path).is_absolute():
        return False
    if str(_absolute(Path(recorded_path))) != recorded_path:
        return False
    if not isinstance(size_bytes, int) or isinstance(size_bytes, bool) or size_bytes < 0:
        return False
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or digest != digest.lower()
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        return False
    try:
        return file_identity(path) == dict(expected)
    except (OSError, QaArtifactError):
        return False


def artifact_identity_matches(path: Path, expected: Mapping[str, object]) -> bool:
    if set(expected) != IDENTITY_KEYS:
        return False
    try:
        _validate_identity_record(expected, "expected artifact identity")
        return artifact_identity(path) == dict(expected)
    except (OSError, QaArtifactError):
        return False


def _write_new_file(path: Path, payload: bytes, mode: int = 0o444) -> None:
    """Create one file exclusively; never replace or truncate an existing path."""

    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    descriptor = os.open(path, flags, mode)
    try:
        offset = 0
        while offset < len(payload):
            offset += os.write(descriptor, payload[offset:])
        os.fsync(descriptor)
        os.fchmod(descriptor, mode)
        observed = os.fstat(descriptor)
        require(stat.S_ISREG(observed.st_mode), f"new path is not a regular file: {path}")
        require(observed.st_nlink == 1, f"new path is aliased: {path}")
    finally:
        os.close(descriptor)


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _rename_noreplace(source: Path, destination: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    require(renameat2 is not None, "renameat2(RENAME_NOREPLACE) is unavailable")
    renameat2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    renameat2.restype = ctypes.c_int
    result = renameat2(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
    if result != 0:
        observed_errno = ctypes.get_errno()
        if observed_errno in {errno.EEXIST, errno.ENOTEMPTY}:
            raise QaArtifactError(f"refusing to overwrite QA directory: {destination}")
        raise QaArtifactError(f"atomic QA publication failed: {os.strerror(observed_errno)}")


def _json_payload(document: Mapping[str, object]) -> bytes:
    return json.dumps(document, ensure_ascii=False, allow_nan=False, indent=2).encode("utf-8") + b"\n"


def _strict_json_load(payload: bytes, label: str) -> object:
    def reject_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise QaArtifactError(f"duplicate JSON key in {label}: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> object:
        raise QaArtifactError(f"non-finite JSON value in {label}: {value}")

    try:
        return json.loads(payload.decode("utf-8"), object_pairs_hook=reject_duplicates, parse_constant=reject_constant)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise QaArtifactError(f"invalid JSON in {label}: {exc}") from exc


def _exact_keys(value: object, expected: set[str], label: str) -> Mapping[str, object]:
    require(isinstance(value, dict), f"{label} is not an object")
    require(set(value) == expected, f"{label} keys differ: expected {sorted(expected)}, got {sorted(value)}")
    return value


def _validate_identity_record(
    value: object,
    label: str,
    *,
    required_mode: str | None = None,
) -> Mapping[str, object]:
    row = _exact_keys(value, IDENTITY_KEYS, label)
    path = row["path"]
    size_bytes = row["size_bytes"]
    digest = row["sha256"]
    mode = row["mode_octal"]
    nlink = row["st_nlink"]
    require(isinstance(path, str) and Path(path).is_absolute(), f"{label}.path is not absolute")
    require(str(_absolute(Path(path))) == path, f"{label}.path is not lexically canonical")
    require(isinstance(size_bytes, int) and not isinstance(size_bytes, bool) and size_bytes >= 0, f"{label}.size_bytes invalid")
    require(isinstance(digest, str) and len(digest) == 64 and digest == digest.lower(), f"{label}.sha256 invalid")
    require(all(character in "0123456789abcdef" for character in digest), f"{label}.sha256 is not hexadecimal")
    require(isinstance(mode, str) and len(mode) == 4 and mode.isdigit(), f"{label}.mode_octal invalid")
    require(all(character in "01234567" for character in mode), f"{label}.mode_octal is not octal")
    require(nlink == 1 and not isinstance(nlink, bool), f"{label}.st_nlink is not one")
    if required_mode is not None:
        require(mode == required_mode, f"{label}.mode_octal is not {required_mode}")
    return row


def build_receipt(
    *,
    html_path: Path,
    html_input: Mapping[str, object],
    qa_producer: Mapping[str, object],
    expected_status: str,
    checks: Mapping[str, bool],
    metrics: Mapping[str, object],
    outputs: Mapping[str, Mapping[str, object]],
    output_dir: Path,
    created_at_utc: str | None = None,
) -> dict[str, object]:
    """Build the schema-2, content-bound browser-QA success receipt."""

    require(expected_status in {"partial", "final"}, "expected_status is invalid")
    require(set(checks) == set(CHECK_NAMES), "browser QA check allowlist mismatch")
    require(all(value is True for value in checks.values()), "cannot build a success receipt from failing checks")
    require(set(metrics) == set(METRIC_NAMES), "browser QA metric allowlist mismatch")
    require(set(outputs) == set(OUTPUT_FILENAMES), "browser QA output role allowlist mismatch")
    root = _absolute(output_dir)
    timestamp = created_at_utc or dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")
    return {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "artifact_class": ARTIFACT_CLASS,
        "created_at_utc": timestamp,
        # Kept for the legacy reader while inputs.html is the authoritative identity.
        "html": str(_absolute(html_path)),
        "inputs": {"html": dict(html_input), "qa_producer": dict(qa_producer)},
        "expected_status": expected_status,
        "checks": dict(checks),
        "metrics": dict(metrics),
        "all_pass": True,
        "outputs": {role: dict(outputs[role]) for role in OUTPUT_FILENAMES},
        "layout": {
            "qa_root": str(root),
            "receipt_path": str(root / RECEIPT_NAME),
            "sidecar_path": str(root / SIDECAR_NAME),
            "directory_mode_octal": "0555",
            "exact_files": list(EXACT_QA_FILES),
            "file_mode_octal": "0444",
            "file_st_nlink": 1,
        },
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": SIDECAR_NAME,
            "covers": RECEIPT_NAME,
        },
        "assurance": dict(ASSURANCE),
    }


def _validate_receipt_document(
    document: object,
    *,
    output_dir: Path,
    html_path: Path,
    producer_path: Path,
    expected_status: str | None,
) -> Mapping[str, object]:
    top = _exact_keys(
        document,
        {
            "schema_name", "schema_version", "task_type", "artifact_class", "created_at_utc",
            "html", "inputs", "expected_status", "checks", "metrics", "all_pass", "outputs",
            "layout", "receipt_integrity", "assurance",
        },
        "browser QA receipt",
    )
    require(top["schema_name"] == SCHEMA_NAME and top["schema_version"] == SCHEMA_VERSION, "browser QA schema mismatch")
    require(top["task_type"] == "B_COMPREHENSIVE_VALIDATION", "browser QA task type mismatch")
    require(top["artifact_class"] == ARTIFACT_CLASS, "browser QA artifact class mismatch")
    try:
        created = dt.datetime.fromisoformat(str(top["created_at_utc"]).replace("Z", "+00:00"))
    except ValueError as exc:
        raise QaArtifactError("browser QA creation time is invalid") from exc
    require(created.tzinfo is not None, "browser QA creation time lacks timezone")
    status_value = top["expected_status"]
    require(status_value in {"partial", "final"}, "browser QA expected_status invalid")
    if expected_status is not None:
        require(status_value == expected_status, "browser QA expected_status differs from requested status")

    inputs = _exact_keys(top["inputs"], {"html", "qa_producer"}, "browser QA inputs")
    html_row = _validate_identity_record(inputs["html"], "inputs.html", required_mode="0444")
    producer_row = _validate_identity_record(inputs["qa_producer"], "inputs.qa_producer")
    require(html_row["path"] == str(_absolute(html_path)), "receipt HTML path differs from requested HTML")
    require(producer_row["path"] == str(_absolute(producer_path)), "receipt producer path differs from executing verifier")
    require(top["html"] == html_row["path"], "legacy HTML path differs from inputs.html")

    checks = _exact_keys(top["checks"], set(CHECK_NAMES), "browser QA checks")
    require(all(checks[name] is True for name in CHECK_NAMES), "browser QA receipt contains a failing check")
    metrics = _exact_keys(top["metrics"], set(METRIC_NAMES), "browser QA metrics")
    require(top["all_pass"] is True, "browser QA all_pass is not true")
    require(isinstance(metrics["title"], str) and isinstance(metrics["status_text"], str), "browser QA text metrics invalid")

    root = _absolute(output_dir)
    outputs = _exact_keys(top["outputs"], set(OUTPUT_FILENAMES), "browser QA outputs")
    for role, filename in OUTPUT_FILENAMES.items():
        row = _validate_identity_record(outputs[role], f"outputs.{role}", required_mode="0444")
        require(row["path"] == str(root / filename), f"browser QA output path mismatch: {role}")

    layout = _exact_keys(
        top["layout"],
        {
            "qa_root", "receipt_path", "sidecar_path", "directory_mode_octal",
            "exact_files", "file_mode_octal", "file_st_nlink",
        },
        "browser QA layout",
    )
    require(layout == {
        "qa_root": str(root),
        "receipt_path": str(root / RECEIPT_NAME),
        "sidecar_path": str(root / SIDECAR_NAME),
        "directory_mode_octal": "0555",
        "exact_files": list(EXACT_QA_FILES),
        "file_mode_octal": "0444",
        "file_st_nlink": 1,
    }, "browser QA exact layout contract mismatch")
    require(top["receipt_integrity"] == {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": SIDECAR_NAME,
        "covers": RECEIPT_NAME,
    }, "browser QA sidecar contract mismatch")
    require(top["assurance"] == ASSURANCE, "browser QA assurance contract mismatch")
    return top


def local_link_checks(html_path: Path, hrefs: list[str]) -> tuple[list[str], list[str]]:
    missing_files: list[str] = []
    external_links: list[str] = []
    for href in hrefs:
        if not href or href.startswith("#"):
            continue
        parsed = urlparse(href)
        if parsed.scheme in {"http", "https"}:
            external_links.append(href)
            continue
        if parsed.scheme not in {"", "file"}:
            continue
        raw_path = unquote(parsed.path)
        target = Path(raw_path) if parsed.scheme == "file" else html_path.parent / raw_path
        if not target.resolve().exists():
            missing_files.append(href)
    return sorted(set(missing_files)), sorted(set(external_links))


def _candidate_output_identities(candidate: Path, output_dir: Path) -> dict[str, dict[str, object]]:
    identities: dict[str, dict[str, object]] = {}
    seen_inodes: set[tuple[int, int]] = set()
    for role, filename in OUTPUT_FILENAMES.items():
        path = candidate / filename
        observed = os.lstat(path)
        require(stat.S_ISREG(observed.st_mode) and not stat.S_ISLNK(observed.st_mode), f"QA output is not regular: {path}")
        require(observed.st_nlink == 1, f"QA output hardlink count is not one: {path}")
        inode = (observed.st_dev, observed.st_ino)
        require(inode not in seen_inodes, f"QA output aliases another role: {path}")
        seen_inodes.add(inode)
        os.chmod(path, 0o444)
        identities[role] = artifact_identity(path, recorded_path=output_dir / filename)
    return identities


def _set_candidate_readonly(candidate: Path) -> None:
    for child in candidate.iterdir():
        observed = os.lstat(child)
        if stat.S_ISREG(observed.st_mode) and not stat.S_ISLNK(observed.st_mode):
            os.chmod(child, 0o444)
    os.chmod(candidate, 0o555)
    _fsync_directory(candidate)


def _failure_destination(output_dir: Path, kind: str) -> Path:
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return output_dir.parent / f"{output_dir.name}.{kind}.{stamp}.{os.getpid()}"


def _invalidate_success_names(root: Path) -> None:
    """Preserve receipt bytes while removing every valid-success filename."""

    for name in (RECEIPT_NAME, SIDECAR_NAME):
        source = root / name
        if os.path.lexists(source):
            destination = root / f"UNAUTHENTICATED_{name}"
            _rename_noreplace(source, destination)


def _preserve_failed_staging(candidate: Path, output_dir: Path, message: str) -> Path:
    """Preserve a never-published candidate with an explicit failure marker."""

    if not os.path.lexists(candidate):
        return candidate
    try:
        _invalidate_success_names(candidate)
        marker = candidate / "QA_NOT_PUBLISHED.txt"
        if not os.path.lexists(marker):
            _write_new_file(marker, (message + "\n").encode("utf-8"))
        _set_candidate_readonly(candidate)
        destination = _failure_destination(output_dir, "failed-staging")
        _rename_noreplace(candidate, destination)
        _fsync_directory(output_dir.parent)
        return destination
    except Exception:
        # The original candidate remains intact and is still visibly a staging path.
        return candidate


def _archive_failed_publication(output_dir: Path, message: str) -> Path:
    """Make a failed post-rename tree unambiguously non-verifiable, then preserve it."""

    if not os.path.lexists(output_dir):
        return output_dir
    try:
        os.chmod(output_dir, 0o755)
        _invalidate_success_names(output_dir)
        marker = output_dir / "PUBLICATION_VERIFICATION_FAILED.txt"
        if not os.path.lexists(marker):
            _write_new_file(marker, (message + "\n").encode("utf-8"))
        _set_candidate_readonly(output_dir)
        destination = _failure_destination(output_dir, "failed-publication")
        _rename_noreplace(output_dir, destination)
        _fsync_directory(output_dir.parent)
        return destination
    except Exception:
        return output_dir


def seal_qa_artifact(
    candidate: Path,
    output_dir: Path,
    *,
    html_path: Path,
    producer_path: Path,
    expected_status: str,
    checks: Mapping[str, bool],
    metrics: Mapping[str, object],
    pre_html_identity: Mapping[str, object],
    pre_producer_identity: Mapping[str, object],
    created_at_utc: str | None = None,
) -> dict[str, object]:
    """Seal and atomically publish an already rendered three-output candidate."""

    candidate = _absolute(candidate)
    output_dir = _absolute(output_dir)
    html_path = _absolute(html_path)
    producer_path = _absolute(producer_path)
    require(not os.path.lexists(output_dir), f"output directory already exists: {output_dir}")
    candidate_stat = os.lstat(candidate)
    require(stat.S_ISDIR(candidate_stat.st_mode) and not stat.S_ISLNK(candidate_stat.st_mode), "QA candidate is not a real directory")
    require({child.name for child in candidate.iterdir()} == set(OUTPUT_FILENAMES.values()), "QA candidate output set differs from allowlist")
    require(file_identity_matches(html_path, pre_html_identity), "HTML changed before QA sealing")
    require(artifact_identity_matches(producer_path, pre_producer_identity), "QA producer changed before QA sealing")
    require(set(checks) == set(CHECK_NAMES) and all(checks.values()), "refusing to seal failing QA checks")

    html_state = os.lstat(html_path)
    require(stat.S_ISREG(html_state.st_mode) and not stat.S_ISLNK(html_state.st_mode), "HTML is not a regular file")
    require(html_state.st_nlink == 1, "HTML hardlink count is not one")
    os.chmod(html_path, 0o444)
    html_identity = artifact_identity(html_path)
    require(html_identity["mode_octal"] == "0444", "HTML mode did not become 0444")
    require(
        {key: html_identity[key] for key in ("path", "size_bytes", "sha256")} == dict(pre_html_identity),
        "HTML changed between QA check and sealing identity capture",
    )
    producer_identity = artifact_identity(producer_path)
    require(
        producer_identity == dict(pre_producer_identity),
        "QA producer changed between QA check and sealing identity capture",
    )
    outputs = _candidate_output_identities(candidate, output_dir)
    receipt = build_receipt(
        html_path=html_path,
        html_input=html_identity,
        qa_producer=producer_identity,
        expected_status=expected_status,
        checks=checks,
        metrics=metrics,
        outputs=outputs,
        output_dir=output_dir,
        created_at_utc=created_at_utc,
    )
    _validate_receipt_document(
        receipt,
        output_dir=output_dir,
        html_path=html_path,
        producer_path=producer_path,
        expected_status=expected_status,
    )
    receipt_payload = _json_payload(receipt)
    _write_new_file(candidate / RECEIPT_NAME, receipt_payload)
    require(_stable_read(candidate / RECEIPT_NAME, "staged browser QA receipt")[0] == receipt_payload, "staged receipt bytes differ")
    require(not os.path.lexists(candidate / SIDECAR_NAME), "staged sidecar unexpectedly exists")

    published = False
    try:
        require(not os.path.lexists(output_dir), f"output directory appeared before publication: {output_dir}")
        _rename_noreplace(candidate, output_dir)
        published = True
        digest = hashlib.sha256(receipt_payload).hexdigest()
        _write_new_file(output_dir / SIDECAR_NAME, f"{digest}  {RECEIPT_NAME}\n".encode("ascii"))
        for name in EXACT_QA_FILES:
            os.chmod(output_dir / name, 0o444)
        os.chmod(output_dir, 0o555)
        _fsync_directory(output_dir)
        _fsync_directory(output_dir.parent)
        verification = verify_qa_artifact(
            output_dir,
            html_path=html_path,
            producer_path=producer_path,
            expected_status=expected_status,
        )
        require(verification["all_pass"] is True, "post-publication QA verification failed")
        return receipt
    except Exception as exc:
        if published:
            preserved = _archive_failed_publication(output_dir, str(exc))
            raise QaArtifactError(f"QA publication failed verification; diagnostics preserved at {preserved}: {exc}") from exc
        raise


def verify_qa_artifact(
    output_dir: Path,
    *,
    html_path: Path,
    producer_path: Path | None = None,
    expected_status: str | None = None,
) -> dict[str, object]:
    """Independently revalidate the exact immutable QA artifact."""

    root = _absolute(output_dir)
    html = _absolute(html_path)
    producer = _absolute(Path(__file__) if producer_path is None else producer_path)
    try:
        root_state = os.lstat(root)
    except OSError as exc:
        raise QaArtifactError(f"QA directory is unavailable: {root}: {exc}") from exc
    require(stat.S_ISDIR(root_state.st_mode) and not stat.S_ISLNK(root_state.st_mode), "QA root is not a real directory")
    require(stat.S_IMODE(root_state.st_mode) == 0o555, "QA root mode is not 0555")

    entries = list(root.iterdir())
    require({entry.name for entry in entries} == set(EXACT_QA_FILES), "QA directory differs from exact five-file layout")
    inode_set: set[tuple[int, int]] = set()
    for entry in entries:
        observed = os.lstat(entry)
        require(stat.S_ISREG(observed.st_mode) and not stat.S_ISLNK(observed.st_mode), f"QA entry is not a regular file: {entry}")
        require(stat.S_IMODE(observed.st_mode) == 0o444, f"QA entry mode is not 0444: {entry}")
        require(observed.st_nlink == 1, f"QA entry hardlink count is not one: {entry}")
        inode = (observed.st_dev, observed.st_ino)
        require(inode not in inode_set, f"QA entry aliases another file: {entry}")
        inode_set.add(inode)

    receipt_payload, _ = _stable_read(root / RECEIPT_NAME, "browser QA receipt")
    sidecar_payload, _ = _stable_read(root / SIDECAR_NAME, "browser QA sidecar")
    digest = hashlib.sha256(receipt_payload).hexdigest()
    require(sidecar_payload == f"{digest}  {RECEIPT_NAME}\n".encode("ascii"), "browser QA sidecar mismatch")
    document = _strict_json_load(receipt_payload, str(root / RECEIPT_NAME))
    receipt = _validate_receipt_document(
        document,
        output_dir=root,
        html_path=html,
        producer_path=producer,
        expected_status=expected_status,
    )

    html_row = receipt["inputs"]["html"]
    producer_row = receipt["inputs"]["qa_producer"]
    require(artifact_identity_matches(html, html_row), "bound HTML identity changed after QA")
    require(artifact_identity_matches(producer, producer_row), "bound QA producer identity changed after QA")
    external_inodes: set[tuple[int, int]] = set(inode_set)
    for label, path in (("HTML", html), ("QA producer", producer)):
        observed = os.lstat(path)
        inode = (observed.st_dev, observed.st_ino)
        require(inode not in external_inodes, f"{label} aliases a QA artifact file")
        external_inodes.add(inode)
    require(html != producer, "HTML and QA producer paths alias")

    outputs = receipt["outputs"]
    for role, filename in OUTPUT_FILENAMES.items():
        observed_identity = artifact_identity(root / filename)
        require(observed_identity == outputs[role], f"QA output identity changed: {role}")
    return {
        "all_pass": True,
        "artifact_class": ARTIFACT_CLASS,
        "qa_root": str(root),
        "receipt": str(root / RECEIPT_NAME),
        "receipt_sha256": digest,
        "html_sha256": html_row["sha256"],
        "producer_sha256": producer_row["sha256"],
        "outputs_verified": len(OUTPUT_FILENAMES),
        "exact_layout_verified": True,
    }


def _run_browser_qa(html_path: Path, candidate: Path, expected_status: str) -> tuple[dict[str, bool], dict[str, object]]:
    console_errors: list[str] = []
    page_errors: list[str] = []
    failed_requests: list[str] = []
    remote_requests: list[str] = []

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000}, device_scale_factor=1)
        page.on("console", lambda message: console_errors.append(message.text) if message.type == "error" else None)
        page.on("pageerror", lambda error: page_errors.append(str(error)))
        page.on("requestfailed", lambda request: failed_requests.append(request.url))
        page.on(
            "request",
            lambda request: remote_requests.append(request.url)
            if urlparse(request.url).scheme in {"http", "https"}
            else None,
        )
        page.goto(html_path.as_uri(), wait_until="networkidle")

        status_text = page.locator(".status-ribbon").inner_text()
        title = page.title()
        html_lang = page.locator("html").get_attribute("lang")
        h1_count = page.locator("h1").count()
        section_count = page.locator("section").count()
        svg_count = page.locator("svg").count()
        details_count = page.locator("details").count()
        ids = page.locator("[id]").evaluate_all("els => els.map(el => el.id)")
        hrefs = page.locator("a[href]").evaluate_all("els => els.map(el => el.getAttribute('href'))")
        internal_targets = page.locator("a[href^='#']").evaluate_all(
            "els => els.map(el => el.getAttribute('href').slice(1))"
        )
        missing_anchor_targets = sorted({target for target in internal_targets if target and target not in ids})
        desktop_widths = page.evaluate(
            "() => ({scroll: document.documentElement.scrollWidth, client: document.documentElement.clientWidth})"
        )
        page.screenshot(path=str(candidate / OUTPUT_FILENAMES["desktop_screenshot"]), full_page=True)

        details_toggle_ok = True
        if details_count:
            first_details = page.locator("details").first
            before = first_details.get_attribute("open") is not None
            first_details.locator("summary").click()
            after = first_details.get_attribute("open") is not None
            details_toggle_ok = before != after

        page.set_viewport_size({"width": 390, "height": 844})
        page.wait_for_timeout(100)
        narrow_widths = page.evaluate(
            "() => ({scroll: document.documentElement.scrollWidth, client: document.documentElement.clientWidth})"
        )
        page.screenshot(path=str(candidate / OUTPUT_FILENAMES["mobile_screenshot"]), full_page=True)

        page.emulate_media(media="print")
        page.pdf(path=str(candidate / OUTPUT_FILENAMES["print_pdf"]), format="A4", print_background=True)
        browser.close()

    missing_files, external_links = local_link_checks(html_path, hrefs)
    expected_status_ok = (
        "PARTIAL PREVIEW" in status_text
        if expected_status == "partial"
        else "FINAL" in status_text and "PARTIAL" not in status_text
    )
    checks: dict[str, bool] = {
        "title_present": bool(title.strip()),
        "lang_is_zh_hant": html_lang == "zh-Hant",
        "one_h1": h1_count == 1,
        "minimum_sections": section_count >= 8,
        "minimum_svg_figures": svg_count >= (5 if expected_status == "partial" else 6),
        "collapsible_definitions_present": details_count >= 5,
        "details_toggle_works": details_toggle_ok,
        "unique_ids": len(ids) == len(set(ids)),
        "all_internal_anchors_resolve": not missing_anchor_targets,
        "all_local_link_targets_exist": not missing_files,
        "no_external_links": not external_links,
        "no_remote_requests": not remote_requests,
        "no_failed_requests": not failed_requests,
        "no_console_errors": not console_errors,
        "no_page_errors": not page_errors,
        "desktop_no_body_overflow": desktop_widths["scroll"] <= desktop_widths["client"],
        "mobile_no_body_overflow": narrow_widths["scroll"] <= narrow_widths["client"],
        "status_matches_expected": expected_status_ok,
        # Filled by run_browser_qa at the same boundary as producer verification.
        "html_unchanged_during_qa": False,
        "qa_producer_unchanged_during_qa": False,
    }
    metrics: dict[str, object] = {
        "title": title,
        "status_text": status_text,
        "h1_count": h1_count,
        "section_count": section_count,
        "svg_count": svg_count,
        "details_count": details_count,
        "desktop_widths": desktop_widths,
        "mobile_widths": narrow_widths,
        "missing_anchor_targets": missing_anchor_targets,
        "missing_local_files": missing_files,
        "external_links": external_links,
        "remote_requests": sorted(set(remote_requests)),
        "failed_requests": sorted(set(failed_requests)),
        "console_errors": console_errors,
        "page_errors": page_errors,
    }
    return checks, metrics


def run_browser_qa(html_path: Path, output_dir: Path, expected_status: str) -> dict[str, object]:
    """Render, validate, seal, and publish one immutable browser QA artifact."""

    html = _absolute(html_path)
    output = _absolute(output_dir)
    producer = _absolute(Path(__file__))
    require(expected_status in {"partial", "final"}, "expect-status must be partial or final")
    require(not os.path.lexists(output), f"refusing to overwrite QA directory: {output}")
    pre_html = file_identity(html)
    pre_producer = artifact_identity(producer)
    output.parent.mkdir(parents=True, exist_ok=True)
    parent_state = os.lstat(output.parent)
    require(stat.S_ISDIR(parent_state.st_mode) and not stat.S_ISLNK(parent_state.st_mode), "QA output parent is not a real directory")
    require(not os.path.lexists(output), f"QA output appeared before staging: {output}")
    candidate = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=str(output.parent)))
    try:
        checks, metrics = _run_browser_qa(html, candidate, expected_status)
        checks["html_unchanged_during_qa"] = file_identity_matches(html, pre_html)
        checks["qa_producer_unchanged_during_qa"] = artifact_identity_matches(producer, pre_producer)
        if not all(checks.values()):
            diagnostic = {
                "schema_name": "intersubmod.validated_html_browser_qa_diagnostic",
                "schema_version": "1.0.0",
                "all_pass": False,
                "expected_status": expected_status,
                "checks": checks,
                "metrics": metrics,
                "note": "NOT A VALIDATION RECEIPT; no authentication sidecar is issued",
            }
            _write_new_file(candidate / "browser_qa_diagnostic.json", _json_payload(diagnostic))
            raise QaArtifactError("browser QA checks failed")
        return seal_qa_artifact(
            candidate,
            output,
            html_path=html,
            producer_path=producer,
            expected_status=expected_status,
            checks=checks,
            metrics=metrics,
            pre_html_identity=pre_html,
            pre_producer_identity=pre_producer,
        )
    except Exception as exc:
        if os.path.lexists(candidate):
            preserved = _preserve_failed_staging(candidate, output, str(exc))
            raise QaArtifactError(f"browser QA failed; diagnostics preserved at {preserved}: {exc}") from exc
        raise


def _build_legacy_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expect-status", choices=("partial", "final"), required=True)
    return parser


def _build_verify_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Independently verify an immutable browser QA artifact")
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expect-status", choices=("partial", "final"))
    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    tokens = list(sys.argv[1:] if argv is None else argv)
    if tokens and tokens[0] == "verify-only":
        args = _build_verify_parser().parse_args(tokens[1:])
        args.command = "verify-only"
        return args
    if tokens and tokens[0] == "run":
        tokens = tokens[1:]
    args = _build_legacy_run_parser().parse_args(tokens)
    args.command = "run"
    return args


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        if args.command == "verify-only":
            result = verify_qa_artifact(
                args.output_dir,
                html_path=args.html,
                expected_status=args.expect_status,
            )
        else:
            result = run_browser_qa(args.html, args.output_dir, args.expect_status)
    except (QaArtifactError, OSError, ValueError) as exc:
        print(json.dumps({"all_pass": False, "error": str(exc)}, ensure_ascii=False, indent=2))
        raise SystemExit(2)
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
