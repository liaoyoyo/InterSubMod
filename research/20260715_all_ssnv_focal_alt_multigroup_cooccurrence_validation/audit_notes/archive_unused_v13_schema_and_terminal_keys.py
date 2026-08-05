#!/usr/bin/env python3
"""Archive the unused v13 authority and terminal keys after pre-authority rejection."""

from __future__ import annotations

import ctypes
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Any


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
RESULT_ROOT = TOPIC_ROOT / "results"
REVIEW_ROOT = TOPIC_ROOT / "reviews"
EVIDENCE_ROOT = (
    AUDIT_ROOT
    / "rejected_pre_authority_reviews"
    / "20260723_v13_candidate_nash_medium_findings"
)
LEDGER = EVIDENCE_ROOT / "key_archive_ledger.v1.jsonl"
STATUS_FILENAME = "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
AT_FDCWD = -100
RENAME_NOREPLACE = 1

ENTRIES = (
    {
        "role": "schema_recovery_authority_v13",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_verifier_schema_recovery/20260723_v13"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_verifier_schema_recovery/archive/"
            "20260723_v13_unused_pre_authority_review_rejection_01"
        ),
        "public_filename": "ed25519_public.pem",
        "private_filename": "ed25519_private_one_time.pem",
        "public_sha256": (
            "99e8836c17adfa0559ad92b8605f2c8b1fd7e0a2064de207e34862b4f48ee56c"
        ),
    },
    {
        "role": "downstream_terminal_v3",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_downstream_continuation/20260723_m2v5_terminal_v3"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_downstream_continuation/archive/"
            "20260723_m2v5_terminal_v3_unused_pre_authority_review_rejection_01"
        ),
        "public_filename": "ed25519_public.pem",
        "private_filename": "ed25519_private_one_time_resident.pem",
        "public_sha256": (
            "a1bbafabe577ae3a05dffc0e61566d51ad7436baf980590f9fdac51179ca0d94"
        ),
    },
)

FORMAL_V13_SLOTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v13.json.ed25519.sig",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_mendel.v13.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_nash.v13.json",
    REVIEW_ROOT / "20260723_tumor_ref_schema_recovery_external_claude_opus.v13.json",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v13.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v13.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v13.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v13.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v13.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v13.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v13.json",
)


class ArchiveError(RuntimeError):
    """Raised when the append-only archive contract cannot be satisfied."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mode(path: Path) -> str:
    return oct(stat.S_IMODE(os.lstat(path).st_mode))


def fsync_directory(path: Path) -> None:
    descriptor = os.open(
        path,
        os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0),
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def rename_no_replace(source: Path, target: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise ArchiveError("renameat2 is required")
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    if (
        renameat2(
            AT_FDCWD,
            os.fsencode(source),
            AT_FDCWD,
            os.fsencode(target),
            RENAME_NOREPLACE,
        )
        != 0
    ):
        error_number = ctypes.get_errno()
        raise ArchiveError(
            f"No-replace archive move failed: {os.strerror(error_number)}"
        )


def encode_json(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, sort_keys=True) + "\n"
    ).encode("utf-8")


def write_exclusive_readonly(path: Path, value: Any) -> None:
    data = encode_json(value)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o400)
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise ArchiveError(f"Short write: {path}")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def append_event(descriptor: int, value: Any) -> None:
    data = encode_json(value)
    offset = 0
    while offset < len(data):
        written = os.write(descriptor, data[offset:])
        if written <= 0:
            raise ArchiveError("Short append to key archive ledger")
        offset += written
    os.fsync(descriptor)


def preflight(entry: dict[str, Any]) -> dict[str, Any]:
    source = entry["source"]
    archive = entry["archive"]
    if not source.is_dir() or source.is_symlink() or os.path.lexists(archive):
        raise ArchiveError(f"Unsafe source/archive state for {entry['role']}")
    observed_names = sorted(path.name for path in source.iterdir())
    expected_names = sorted(
        (entry["public_filename"], entry["private_filename"])
    )
    if observed_names != expected_names:
        raise ArchiveError(f"Unexpected key-directory contents: {observed_names}")
    public_key = source / entry["public_filename"]
    private_key = source / entry["private_filename"]
    if (
        mode(source) != "0o700"
        or mode(public_key) != "0o444"
        or mode(private_key) != "0o400"
        or sha256(public_key) != entry["public_sha256"]
    ):
        raise ArchiveError(f"Key identity/mode drift for {entry['role']}")
    return {
        "role": entry["role"],
        "source": str(source),
        "archive": str(archive),
        "reason": "REJECTED_PRE_AUTHORITY_BY_INDEPENDENT_MEDIUM_FINDINGS",
        "signature_created": False,
        "must_never_be_used": True,
        "public_key": {
            "path": str(public_key),
            "mode": mode(public_key),
            "size_bytes": public_key.stat().st_size,
            "sha256": sha256(public_key),
        },
        "private_key": {
            "path": str(private_key),
            "mode": mode(private_key),
            "size_bytes": private_key.stat().st_size,
            "sha256": sha256(private_key),
        },
    }


def main() -> None:
    if os.path.lexists(LEDGER):
        raise ArchiveError(f"Archive ledger already exists: {LEDGER}")
    occupied = [str(path) for path in FORMAL_V13_SLOTS if os.path.lexists(path)]
    if occupied:
        raise ArchiveError(f"Formal v13 slots are occupied: {occupied}")
    if not EVIDENCE_ROOT.is_dir() or EVIDENCE_ROOT.is_symlink():
        raise ArchiveError(f"Evidence root is missing or unsafe: {EVIDENCE_ROOT}")
    records = [preflight(dict(entry)) for entry in ENTRIES]
    for entry in ENTRIES:
        parent = entry["archive"].parent
        if os.path.lexists(parent):
            if not parent.is_dir() or parent.is_symlink():
                raise ArchiveError(f"Unsafe archive parent: {parent}")
        else:
            os.mkdir(parent, 0o700)
            fsync_directory(parent.parent)

    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    ledger_fd = os.open(LEDGER, flags, 0o600)
    completed: list[dict[str, Any]] = []
    try:
        append_event(
            ledger_fd,
            {
                "event": "PREFLIGHT_COMPLETE",
                "formal_v13_slot_count": len(FORMAL_V13_SLOTS),
                "formal_v13_slots_all_absent": True,
                "entry_count": len(records),
            },
        )
        for entry, record in zip(ENTRIES, records, strict=True):
            source = entry["source"]
            archive = entry["archive"]
            write_exclusive_readonly(source / STATUS_FILENAME, record)
            fsync_directory(source)
            append_event(
                ledger_fd,
                {"event": "ARCHIVE_MOVE_START", "role": record["role"]},
            )
            rename_no_replace(source, archive)
            fsync_directory(archive)
            fsync_directory(archive.parent)
            fsync_directory(source.parent)
            if os.path.lexists(source) or not archive.is_dir():
                raise ArchiveError(f"Post-move path drift: {record['role']}")
            complete = {
                **record,
                "archive_record": str(archive / STATUS_FILENAME),
                "archive_directory_mode": mode(archive),
                "status": "ARCHIVED_UNUSED_NEVER_SIGNED",
            }
            completed.append(complete)
            append_event(
                ledger_fd,
                {"event": "ARCHIVE_MOVE_COMPLETE", **complete},
            )
        append_event(
            ledger_fd,
            {
                "event": "ARCHIVE_BATCH_COMPLETE",
                "entry_count": len(completed),
                "status": "PASS",
            },
        )
    finally:
        os.fchmod(ledger_fd, 0o444)
        os.fsync(ledger_fd)
        os.close(ledger_fd)
        fsync_directory(EVIDENCE_ROOT)
    print(
        json.dumps(
            {
                "pass": len(completed) == len(ENTRIES),
                "ledger": str(LEDGER),
                "archived_roles": [record["role"] for record in completed],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
