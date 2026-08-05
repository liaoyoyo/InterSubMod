#!/usr/bin/env python3
"""Freeze and verify the final presentation-layer code and evidence bindings.

This contract is intentionally separate from the scientific exact-11 release.
It snapshots exactly three presentation programs, records the exact identity of
the twelve files consumed by the final report, and records the Python,
Playwright, and Chromium runtime.  Evidence inputs are identity-bound but are
not copied; ``verify-only`` therefore fails if any bound input later changes.

The snapshot is published only after a complete staged verification, using a
same-filesystem ``renameat2(RENAME_NOREPLACE)``.  All files are 0444 and all
directories are 0555.  This is provenance for a presentation artifact, not a
new scientific-evidence authority.
"""

from __future__ import annotations

import argparse
import ctypes
import dataclasses
import datetime as dt
import errno
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import stat
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence


RESEARCH_RELATIVE = Path("research/20260716_read_linked_hypercube_exact_likelihood_validation")
SCRIPTS_RELATIVE = RESEARCH_RELATIVE / "scripts"
SOURCE_ALLOWLIST: tuple[tuple[str, Path], ...] = (
    ("html_builder", SCRIPTS_RELATIVE / "build_validated_html_report.py"),
    ("html_browser_qa", SCRIPTS_RELATIVE / "qa_validated_html_report.py"),
    ("numeric_summary_builder", SCRIPTS_RELATIVE / "build_final_numeric_summary.py"),
)
EVIDENCE_ROLES = (
    "canonical_json",
    "funnel_receipt",
    "funnel_verification_receipt",
    "m0_receipt",
    "m0_verification_receipt",
    "pilot_receipt",
    "method_audit",
    "literature_audit",
    "m2_extraction_receipt",
    "m2_ranking_receipt",
    "m2_verification_receipt",
    "final_numeric_summary",
)

SCHEMA_NAME = "intersubmod.presentation_snapshot"
SCHEMA_VERSION = "1.0.0"
MANIFEST_NAME = "presentation_snapshot_manifest.json"
SIDECAR_NAME = f"{MANIFEST_NAME}.sha256"
ARTIFACT_CLASS = "PRESENTATION_PROVENANCE_NOT_SCIENTIFIC_EVIDENCE"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
CHECK_NAMES = (
    "exact_three_program_allowlist_frozen",
    "all_program_sources_stable_regular_nonaliased",
    "all_program_copies_match_sources",
    "exact_twelve_evidence_inputs_identity_bound",
    "all_evidence_inputs_stable_regular_nonaliased",
    "runtime_identity_frozen",
    "staged_snapshot_verified_before_atomic_publish",
)


class PresentationSnapshotError(RuntimeError):
    """Fail-closed presentation snapshot or identity violation."""


@dataclasses.dataclass(frozen=True)
class _ReadResult:
    payload: bytes
    persistent_identity: dict[str, Any]
    volatile_identity: tuple[int, int, int, int, int, int]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PresentationSnapshotError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256_bytes(payload)


def _exact_keys(value: Any, expected: set[str], label: str) -> Mapping[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    require(set(value) == expected, f"{label} schema mismatch: {sorted(set(value) ^ expected)}")
    return value


def _reject_duplicate_key(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PresentationSnapshotError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _strict_json_load(payload: bytes, label: str) -> Any:
    try:
        return json.loads(
            payload.decode("utf-8", errors="strict"),
            object_pairs_hook=_reject_duplicate_key,
            parse_constant=lambda value: (_ for _ in ()).throw(
                PresentationSnapshotError(f"non-finite JSON number in {label}: {value}")
            ),
        )
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise PresentationSnapshotError(f"invalid JSON {label}: {exc}") from exc


def _volatile_stat(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_nlink),
        int(value.st_size),
        int(value.st_mtime_ns),
        int(value.st_ctime_ns),
    )


def _stable_read(path: Path, label: str, *, require_single_link: bool = True) -> _ReadResult:
    """Read a regular non-symlink file while detecting path/inode drift."""
    absolute = path.absolute()
    try:
        before_path = os.lstat(absolute)
    except OSError as exc:
        raise PresentationSnapshotError(f"cannot lstat {label}: {absolute}: {exc}") from exc
    require(not stat.S_ISLNK(before_path.st_mode), f"{label} is a symlink: {absolute}")
    require(stat.S_ISREG(before_path.st_mode), f"{label} is not a regular file: {absolute}")
    if require_single_link:
        require(before_path.st_nlink == 1, f"{label} hardlink count is not one: {absolute}")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(absolute, flags)
    except OSError as exc:
        raise PresentationSnapshotError(f"cannot safely open {label}: {absolute}: {exc}") from exc
    try:
        before_fd = os.fstat(descriptor)
        require(_volatile_stat(before_fd) == _volatile_stat(before_path), f"{label} path changed before read")
        pieces: list[bytes] = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            pieces.append(block)
        after_fd = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        after_path = os.lstat(absolute)
    except OSError as exc:
        raise PresentationSnapshotError(f"cannot lstat {label} after read: {absolute}: {exc}") from exc
    require(
        _volatile_stat(before_path) == _volatile_stat(before_fd)
        == _volatile_stat(after_fd) == _volatile_stat(after_path),
        f"{label} changed while reading: {absolute}",
    )
    payload = b"".join(pieces)
    require(len(payload) == before_path.st_size, f"{label} size changed while reading: {absolute}")
    return _ReadResult(
        payload=payload,
        persistent_identity={
            "path": str(absolute),
            "size_bytes": len(payload),
            "sha256": sha256_bytes(payload),
        },
        volatile_identity=_volatile_stat(before_path),
    )


def _validate_file_identity(value: Any, label: str) -> dict[str, Any]:
    row = _exact_keys(value, {"path", "size_bytes", "sha256"}, label)
    require(isinstance(row["path"], str) and Path(row["path"]).is_absolute(), f"{label}.path is not absolute")
    require(
        isinstance(row["size_bytes"], int)
        and not isinstance(row["size_bytes"], bool)
        and row["size_bytes"] >= 0,
        f"{label}.size_bytes is invalid",
    )
    require(
        isinstance(row["sha256"], str) and SHA256_RE.fullmatch(row["sha256"]) is not None,
        f"{label}.sha256 is invalid",
    )
    return dict(row)


def _runtime_identity() -> dict[str, Any]:
    try:
        playwright_version = importlib.metadata.version("playwright")
        import playwright  # type: ignore
        from playwright.sync_api import sync_playwright  # type: ignore
    except (ImportError, importlib.metadata.PackageNotFoundError) as exc:
        raise PresentationSnapshotError(f"Playwright runtime is unavailable: {exc}") from exc
    module_path = Path(playwright.__file__).resolve(strict=True)
    python_path = Path(sys.executable).resolve(strict=True)
    try:
        with sync_playwright() as engine:
            chromium_path = Path(engine.chromium.executable_path).resolve(strict=True)
    except Exception as exc:  # Playwright raises implementation-specific errors.
        raise PresentationSnapshotError(f"cannot resolve Playwright Chromium runtime: {exc}") from exc
    return {
        "python": {
            "executable": _stable_read(python_path, "Python executable").persistent_identity,
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "playwright": {
            "version": playwright_version,
            "module": _stable_read(module_path, "Playwright module").persistent_identity,
        },
        "chromium": {
            "executable": _stable_read(chromium_path, "Playwright Chromium executable").persistent_identity,
        },
        "platform": platform.platform(),
    }


def _validate_runtime(value: Any, label: str) -> dict[str, Any]:
    runtime = _exact_keys(value, {"python", "playwright", "chromium", "platform"}, label)
    python = _exact_keys(runtime["python"], {"executable", "version", "implementation"}, f"{label}.python")
    playwright = _exact_keys(runtime["playwright"], {"version", "module"}, f"{label}.playwright")
    chromium = _exact_keys(runtime["chromium"], {"executable"}, f"{label}.chromium")
    for field, item in (
        ("python.version", python["version"]),
        ("python.implementation", python["implementation"]),
        ("playwright.version", playwright["version"]),
        ("platform", runtime["platform"]),
    ):
        require(isinstance(item, str) and item, f"{label}.{field} must be a non-empty string")
    return {
        "python": {
            "executable": _validate_file_identity(python["executable"], f"{label}.python.executable"),
            "version": python["version"],
            "implementation": python["implementation"],
        },
        "playwright": {
            "version": playwright["version"],
            "module": _validate_file_identity(playwright["module"], f"{label}.playwright.module"),
        },
        "chromium": {
            "executable": _validate_file_identity(chromium["executable"], f"{label}.chromium.executable"),
        },
        "platform": runtime["platform"],
    }


def _discover_repo_root() -> Path:
    script = Path(__file__).resolve()
    for parent in script.parents:
        if all((parent / relative).is_file() for _, relative in SOURCE_ALLOWLIST):
            return parent
    raise PresentationSnapshotError("cannot discover InterSubMod repository root")


def _write_new_file(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    descriptor = os.open(path, flags, 0o444)
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.fchmod(descriptor, 0o444)
        require(os.fstat(descriptor).st_nlink == 1, f"new file is aliased: {path}")
    finally:
        os.close(descriptor)


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _set_readonly_tree(root: Path) -> None:
    directories = [path for path in root.rglob("*") if path.is_dir()]
    for directory in sorted(directories, key=lambda item: len(item.parts), reverse=True):
        os.chmod(directory, 0o555)
        _fsync_directory(directory)
    os.chmod(root, 0o555)
    _fsync_directory(root)


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
            raise PresentationSnapshotError(f"refusing to overwrite presentation snapshot: {destination}")
        raise PresentationSnapshotError(f"atomic presentation snapshot publish failed: {os.strerror(observed_errno)}")


def _preserve_failed_staging(staging: Path, output_root: Path) -> tuple[Path, str | None]:
    """Preserve, never delete, a failed staging tree for audit and recovery."""
    if not staging.exists():
        return staging, "staging path disappeared before failure preservation"
    try:
        _set_readonly_tree(staging)
    except (OSError, PresentationSnapshotError) as exc:
        return staging, f"could not enforce readonly staging modes: {exc}"
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    archived = output_root.parent / f"{output_root.name}.failed-staging.{stamp}.{os.getpid()}"
    try:
        _rename_noreplace(staging, archived)
        _fsync_directory(output_root.parent)
        return archived, None
    except (OSError, PresentationSnapshotError) as exc:
        # The original staging directory remains intact when no-replace rename
        # itself fails.  Its exact path is returned to the caller.
        return staging, f"failed-staging atomic rename was unavailable: {exc}"


def _manifest_payload(document: Mapping[str, Any]) -> bytes:
    return json.dumps(document, ensure_ascii=False, allow_nan=False, indent=2).encode("utf-8") + b"\n"


def _authenticate_manifest(root: Path) -> tuple[dict[str, Any], str]:
    manifest = _stable_read(root / MANIFEST_NAME, "presentation snapshot manifest")
    sidecar = _stable_read(root / SIDECAR_NAME, "presentation snapshot manifest sidecar")
    expected = f"{manifest.persistent_identity['sha256']}  {MANIFEST_NAME}\n".encode("ascii")
    require(sidecar.payload == expected, "presentation snapshot manifest sidecar mismatch")
    document = _strict_json_load(manifest.payload, str(root / MANIFEST_NAME))
    require(isinstance(document, dict), "presentation snapshot manifest is not an object")
    return document, str(manifest.persistent_identity["sha256"])


def _normalize_evidence_inputs(evidence_inputs: Mapping[str, Path]) -> dict[str, Path]:
    require(set(evidence_inputs) == set(EVIDENCE_ROLES), "evidence role allowlist mismatch")
    return {role: Path(evidence_inputs[role]).absolute() for role in EVIDENCE_ROLES}


def freeze_presentation_snapshot(
    evidence_inputs: Mapping[str, Path],
    output_root: Path,
    *,
    _repo_root: Path | None = None,
    _runtime: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Create a new immutable presentation snapshot; never overwrite a path."""
    evidence_paths = _normalize_evidence_inputs(evidence_inputs)
    output_root = output_root.absolute()
    require(not os.path.lexists(output_root), f"output root must not exist: {output_root}")
    repo_root = (_discover_repo_root() if _repo_root is None else Path(_repo_root)).resolve(strict=True)

    source_reads: dict[str, _ReadResult] = {}
    source_inodes: set[tuple[int, int]] = set()
    for role, relative in SOURCE_ALLOWLIST:
        observed = _stable_read(repo_root / relative, f"program source/{role}")
        inode = observed.volatile_identity[:2]
        require(inode not in source_inodes, f"program source aliases another role: {role}")
        source_inodes.add(inode)
        source_reads[role] = observed

    evidence_reads: dict[str, _ReadResult] = {}
    evidence_inodes: set[tuple[int, int]] = set()
    for role, path in evidence_paths.items():
        observed = _stable_read(path, f"evidence/{role}")
        inode = observed.volatile_identity[:2]
        require(inode not in evidence_inodes, f"evidence file aliases another role: {role}")
        require(inode not in source_inodes, f"evidence file aliases a program source: {role}")
        evidence_inodes.add(inode)
        evidence_reads[role] = observed

    runtime = _validate_runtime(
        _runtime_identity() if _runtime is None else dict(_runtime), "freeze runtime"
    )
    output_root.parent.mkdir(parents=True, exist_ok=True)
    require(output_root.parent.is_dir() and not output_root.parent.is_symlink(), "output parent is not a real directory")
    require(not os.path.lexists(output_root), f"output root appeared before staging: {output_root}")
    staging = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.staging.", dir=str(output_root.parent)))
    published = False
    try:
        snapshot_rows: list[dict[str, Any]] = []
        for role, relative in SOURCE_ALLOWLIST:
            snapshot_relative = Path("code_snapshot") / relative.name
            observed = source_reads[role]
            _write_new_file(staging / snapshot_relative, observed.payload)
            copied = _stable_read(staging / snapshot_relative, f"staged program/{role}")
            require(copied.payload == observed.payload, f"staged program differs from source: {role}")
            snapshot_rows.append(
                {
                    "role": role,
                    "repo_relative_path": relative.as_posix(),
                    "source": observed.persistent_identity,
                    "snapshot": {
                        "path": snapshot_relative.as_posix(),
                        "size_bytes": len(copied.payload),
                        "sha256": copied.persistent_identity["sha256"],
                        "mode_octal": "0444",
                    },
                }
            )

        evidence_rows = [
            {"role": role, **evidence_reads[role].persistent_identity}
            for role in EVIDENCE_ROLES
        ]
        producer = _stable_read(Path(__file__), "presentation snapshot freezer").persistent_identity
        entrypoints = {
            role: (Path("code_snapshot") / relative.name).as_posix()
            for role, relative in SOURCE_ALLOWLIST
        }
        checks = {name: True for name in CHECK_NAMES}
        document: dict[str, Any] = {
            "schema_name": SCHEMA_NAME,
            "schema_version": SCHEMA_VERSION,
            "task_type": "B_COMPREHENSIVE_VALIDATION",
            "artifact_class": ARTIFACT_CLASS,
            "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z"),
            "entrypoints": entrypoints,
            "source_snapshot": {
                "repo_root": str(repo_root),
                "entries": snapshot_rows,
                "semantic_sha256": canonical_json_sha256(snapshot_rows),
            },
            "evidence_inputs": {
                "entries": evidence_rows,
                "semantic_sha256": canonical_json_sha256(evidence_rows),
            },
            "runtime": runtime,
            "runtime_semantic_sha256": canonical_json_sha256(runtime),
            "producer": producer,
            "assurance": {
                "scientific_evidence_authority": False,
                "code": "exact three-program byte snapshot",
                "evidence": "external files identity-bound by absolute path, size, and full SHA-256; not copied",
                "runtime": "Python executable, Playwright module/version, and Chromium executable full SHA-256",
                "publication": "same-filesystem renameat2(RENAME_NOREPLACE) after staged verification",
                "threat_model": "checksum provenance; no hostile same-UID actor",
            },
            "checks": checks,
            "all_pass": True,
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": SIDECAR_NAME,
                "covers": MANIFEST_NAME,
            },
        }
        _write_new_file(staging / MANIFEST_NAME, _manifest_payload(document))
        digest = sha256_bytes(_manifest_payload(document))
        _write_new_file(staging / SIDECAR_NAME, f"{digest}  {MANIFEST_NAME}\n".encode("ascii"))
        _set_readonly_tree(staging)

        # Recheck every origin and evidence file immediately before publication.
        for role, relative in SOURCE_ALLOWLIST:
            observed = _stable_read(repo_root / relative, f"publish-boundary source/{role}")
            require(observed == source_reads[role], f"program source changed before publish: {role}")
        for role, path in evidence_paths.items():
            observed = _stable_read(path, f"publish-boundary evidence/{role}")
            require(observed == evidence_reads[role], f"evidence changed before publish: {role}")
        # Production recomputes the real runtime at the publication boundary.
        # The private injection exists only so isolated synthetic tests do not
        # depend on a local Playwright/browser installation.
        staged_verification = verify_presentation_snapshot(
            staging, _runtime=runtime if _runtime is not None else None
        )
        require(staged_verification["all_pass"] is True, "staged presentation snapshot verification failed")
        require(not os.path.lexists(output_root), f"output root appeared before atomic publish: {output_root}")
        _rename_noreplace(staging, output_root)
        published = True
        _fsync_directory(output_root.parent)
        return {
            "all_pass": True,
            "artifact_class": ARTIFACT_CLASS,
            "snapshot_root": str(output_root),
            "manifest": str(output_root / MANIFEST_NAME),
            "manifest_sha256": digest,
            "programs_frozen": len(SOURCE_ALLOWLIST),
            "evidence_inputs_bound": len(EVIDENCE_ROLES),
        }
    except Exception as exc:
        if published:
            raise PresentationSnapshotError(
                f"presentation snapshot was atomically published at {output_root}, "
                f"but post-publish durability confirmation failed: {exc}"
            ) from exc
        preserved, preservation_issue = _preserve_failed_staging(staging, output_root)
        detail = f"; preservation note: {preservation_issue}" if preservation_issue else ""
        raise PresentationSnapshotError(
            f"presentation snapshot freeze failed: {exc}; failed staging preserved at {preserved}{detail}"
        ) from exc


def verify_presentation_snapshot(
    snapshot_root: Path,
    *,
    _runtime: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Verify the immutable code snapshot, external evidence, and runtime."""
    root = snapshot_root.absolute()
    try:
        root_state = os.lstat(root)
    except OSError as exc:
        raise PresentationSnapshotError(f"snapshot root is unavailable: {root}: {exc}") from exc
    require(stat.S_ISDIR(root_state.st_mode) and not stat.S_ISLNK(root_state.st_mode), "snapshot root is not a real directory")
    require(stat.S_IMODE(root_state.st_mode) == 0o555, "snapshot root mode is not 0555")
    document, manifest_sha256 = _authenticate_manifest(root)
    expected_top_keys = {
        "schema_name", "schema_version", "task_type", "artifact_class", "created_at_utc",
        "entrypoints", "source_snapshot", "evidence_inputs", "runtime",
        "runtime_semantic_sha256", "producer", "assurance", "checks", "all_pass",
        "receipt_integrity",
    }
    _exact_keys(document, expected_top_keys, "presentation snapshot manifest")
    require(document["schema_name"] == SCHEMA_NAME and document["schema_version"] == SCHEMA_VERSION, "manifest schema mismatch")
    require(document["task_type"] == "B_COMPREHENSIVE_VALIDATION", "manifest task type mismatch")
    require(document["artifact_class"] == ARTIFACT_CLASS, "manifest artifact class mismatch")
    try:
        created = dt.datetime.fromisoformat(str(document["created_at_utc"]).replace("Z", "+00:00"))
    except ValueError as exc:
        raise PresentationSnapshotError("manifest creation time is invalid") from exc
    require(created.tzinfo is not None, "manifest creation time lacks timezone")
    require(document["all_pass"] is True, "manifest all_pass is not true")
    checks = _exact_keys(document["checks"], set(CHECK_NAMES), "manifest checks")
    require(all(checks[name] is True for name in CHECK_NAMES), "manifest contains a failing named check")
    require(document["receipt_integrity"] == {
        "scheme": "external_sha256_sidecar_v1", "sidecar_name": SIDECAR_NAME, "covers": MANIFEST_NAME,
    }, "manifest receipt integrity contract mismatch")
    require(document["assurance"] == {
        "scientific_evidence_authority": False,
        "code": "exact three-program byte snapshot",
        "evidence": "external files identity-bound by absolute path, size, and full SHA-256; not copied",
        "runtime": "Python executable, Playwright module/version, and Chromium executable full SHA-256",
        "publication": "same-filesystem renameat2(RENAME_NOREPLACE) after staged verification",
        "threat_model": "checksum provenance; no hostile same-UID actor",
    }, "manifest assurance contract mismatch")

    expected_entrypoints = {
        role: (Path("code_snapshot") / relative.name).as_posix()
        for role, relative in SOURCE_ALLOWLIST
    }
    require(document["entrypoints"] == expected_entrypoints, "entrypoint allowlist mismatch")
    source_snapshot = _exact_keys(document["source_snapshot"], {"repo_root", "entries", "semantic_sha256"}, "source_snapshot")
    require(isinstance(source_snapshot["repo_root"], str) and Path(source_snapshot["repo_root"]).is_absolute(), "source repo root is invalid")
    entries = source_snapshot["entries"]
    require(isinstance(entries, list) and len(entries) == len(SOURCE_ALLOWLIST), "source snapshot entry count mismatch")
    require(source_snapshot["semantic_sha256"] == canonical_json_sha256(entries), "source snapshot semantic digest mismatch")
    expected_by_role = dict(SOURCE_ALLOWLIST)
    seen_roles: set[str] = set()
    expected_files = {root / MANIFEST_NAME, root / SIDECAR_NAME}
    snapshot_inodes: set[tuple[int, int]] = set()
    for row in entries:
        _exact_keys(row, {"role", "repo_relative_path", "source", "snapshot"}, "source snapshot entry")
        role = row["role"]
        require(isinstance(role, str) and role in expected_by_role and role not in seen_roles, f"invalid or duplicate source role: {role}")
        seen_roles.add(role)
        require(row["repo_relative_path"] == expected_by_role[role].as_posix(), f"source relative path mismatch: {role}")
        source_identity = _validate_file_identity(row["source"], f"source/{role}")
        snapshot = _exact_keys(row["snapshot"], {"path", "size_bytes", "sha256", "mode_octal"}, f"snapshot/{role}")
        expected_relative = expected_entrypoints[role]
        require(snapshot["path"] == expected_relative and snapshot["mode_octal"] == "0444", f"snapshot path/mode record mismatch: {role}")
        require(snapshot["size_bytes"] == source_identity["size_bytes"] and snapshot["sha256"] == source_identity["sha256"], f"snapshot/source identity mismatch: {role}")
        copy_path = root / expected_relative
        copy = _stable_read(copy_path, f"snapshot program/{role}")
        require(copy.persistent_identity["size_bytes"] == snapshot["size_bytes"] and copy.persistent_identity["sha256"] == snapshot["sha256"], f"snapshot bytes mismatch: {role}")
        require(stat.S_IMODE(os.lstat(copy_path).st_mode) == 0o444, f"snapshot file mode is not 0444: {role}")
        inode = copy.volatile_identity[:2]
        require(inode not in snapshot_inodes, f"snapshot program aliases another role: {role}")
        snapshot_inodes.add(inode)
        expected_files.add(copy_path)
    require(seen_roles == set(expected_by_role), "source role set mismatch")

    evidence = _exact_keys(document["evidence_inputs"], {"entries", "semantic_sha256"}, "evidence_inputs")
    evidence_entries = evidence["entries"]
    require(isinstance(evidence_entries, list) and len(evidence_entries) == len(EVIDENCE_ROLES), "evidence entry count mismatch")
    require(evidence["semantic_sha256"] == canonical_json_sha256(evidence_entries), "evidence semantic digest mismatch")
    seen_evidence: set[str] = set()
    evidence_inodes: set[tuple[int, int]] = set()
    for row in evidence_entries:
        _exact_keys(row, {"role", "path", "size_bytes", "sha256"}, "evidence entry")
        role = row["role"]
        require(isinstance(role, str) and role in EVIDENCE_ROLES and role not in seen_evidence, f"invalid or duplicate evidence role: {role}")
        seen_evidence.add(role)
        expected_identity = _validate_file_identity(
            {key: row[key] for key in ("path", "size_bytes", "sha256")}, f"evidence/{role}"
        )
        observed = _stable_read(Path(expected_identity["path"]), f"evidence/{role}")
        require(observed.persistent_identity == expected_identity, f"bound evidence identity changed: {role}")
        inode = observed.volatile_identity[:2]
        require(inode not in evidence_inodes, f"bound evidence aliases another role: {role}")
        evidence_inodes.add(inode)
    require(seen_evidence == set(EVIDENCE_ROLES), "evidence role set mismatch")

    frozen_runtime = _validate_runtime(document["runtime"], "manifest runtime")
    require(document["runtime_semantic_sha256"] == canonical_json_sha256(frozen_runtime), "runtime semantic digest mismatch")
    current_runtime = _validate_runtime(
        _runtime_identity() if _runtime is None else dict(_runtime), "current runtime"
    )
    require(current_runtime == frozen_runtime, "current runtime differs from frozen presentation runtime")
    _validate_file_identity(document["producer"], "producer")

    actual_files: set[Path] = set()
    actual_directories: set[Path] = {root}
    for path in root.rglob("*"):
        observed = os.lstat(path)
        require(not stat.S_ISLNK(observed.st_mode), f"snapshot contains a symlink: {path}")
        if stat.S_ISDIR(observed.st_mode):
            require(stat.S_IMODE(observed.st_mode) == 0o555, f"snapshot directory mode is not 0555: {path}")
            actual_directories.add(path)
        elif stat.S_ISREG(observed.st_mode):
            require(stat.S_IMODE(observed.st_mode) == 0o444, f"snapshot file mode is not 0444: {path}")
            actual_files.add(path)
        else:
            raise PresentationSnapshotError(f"snapshot contains unsupported filesystem entry: {path}")
    require(actual_directories == {root, root / "code_snapshot"}, "snapshot directory set differs from exact layout")
    require(actual_files == expected_files, "snapshot file set differs from exact allowlist")
    return {
        "all_pass": True,
        "artifact_class": ARTIFACT_CLASS,
        "snapshot_root": str(root),
        "manifest_sha256": manifest_sha256,
        "programs_verified": len(seen_roles),
        "evidence_inputs_verified": len(seen_evidence),
        "runtime_match": True,
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    freeze = subparsers.add_parser("freeze", help="create a new immutable presentation snapshot")
    freeze.add_argument("--output-root", required=True, type=Path)
    for role in EVIDENCE_ROLES:
        freeze.add_argument(f"--{role.replace('_', '-')}", required=True, type=Path)
    verify = subparsers.add_parser("verify-only", help="verify without modifying the snapshot")
    verify.add_argument("--snapshot-root", required=True, type=Path)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    try:
        if args.command == "freeze":
            evidence = {role: getattr(args, role) for role in EVIDENCE_ROLES}
            result = freeze_presentation_snapshot(evidence, args.output_root)
        else:
            result = verify_presentation_snapshot(args.snapshot_root)
    except (PresentationSnapshotError, OSError, ValueError) as exc:
        print(json.dumps({"all_pass": False, "error": str(exc)}, ensure_ascii=False, indent=2))
        raise SystemExit(2)
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
