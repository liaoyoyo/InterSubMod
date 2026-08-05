#!/usr/bin/env python3
"""Archive four unused v5 signer key directories without deletion or overwrite."""

from __future__ import annotations

import ctypes
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Any


HERE = Path(__file__).resolve().parent
LEDGER = HERE / "20260719_unused_v5_signer_key_archive.v1.jsonl"
PRIOR_ARCHIVE_LEDGER = HERE / "20260718_unused_v4_signer_key_archive.v1.jsonl"
PRIOR_ARCHIVE_LEDGER_SHA256 = (
    "b7720d1ab680f9f87b23f238a2f9295df1f78acf9419a1133fb44f492028a9a9"
)
PUBLIC_FILENAME = "ed25519_public.pem"
PRIVATE_FILENAME = "ed25519_private_encrypted_unrecoverable_after_signing.pem"
STATUS_FILENAME = "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
AT_FDCWD = -100
RENAME_NOREPLACE = 1

ENTRIES = (
    {
        "role": "source",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
            "20260718_all_ssnv_v5_summary_hotfix"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/archive/"
            "20260718_all_ssnv_v5_summary_hotfix_unused_v5_terminal_tuple_race_01"
        ),
        "target": Path(
            "/big7_disk/liaoyoyo2001/InterSubMod/docs/provenance/"
            "source_authorities/"
            "20260718_all_ssnv_focal_alt_release_source_authority.v5.approval.json"
        ),
        "public_sha256": (
            "677f73c56f8e3fb998ca4f1e360dd58dd73b90c6389f2f3d3d17350c2f5ec698"
        ),
    },
    {
        "role": "result",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "20260717_all_ssnv_result_v2"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "archive/20260717_all_ssnv_result_v2_unused_v5_terminal_tuple_race_01"
        ),
        "target": Path(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/"
            "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
            "results/task_b_final_dataset_release_receipt.v1.json"
        ),
        "public_sha256": (
            "1c7abc3a62da417ed6c94ab47ec28e58f3a9c0c4249f413a5997cadc53e366c3"
        ),
    },
    {
        "role": "report",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "20260717_all_ssnv_report_v2"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "archive/20260717_all_ssnv_report_v2_unused_v5_terminal_tuple_race_01"
        ),
        "target": Path(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/"
            "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
            "results/task_b_final_report_release_receipt.v1.json"
        ),
        "public_sha256": (
            "0dcf55f8a532ac9f506e1aba703fc3caee1a0cde5de602bfd5dd5a3d35e40a88"
        ),
    },
    {
        "role": "supplemental",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_supplemental_report_integrity/"
            "20260718_positional_singleton_v1"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_supplemental_report_integrity/archive/"
            "20260718_positional_singleton_v1_unused_v5_terminal_tuple_race_01"
        ),
        "target": Path(
            "/big7_disk/liaoyoyo2001/InterSubMod/research/"
            "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
            "results/positional_singleton_supplemental_release_receipt.v1.json"
        ),
        "public_sha256": (
            "7413426516497930be3f73c42bac451a5e0c85021f1cd7660e1c7fc7a5e1cf61"
        ),
    },
)


class ArchiveError(RuntimeError):
    """Raised when the no-delete archive contract cannot be satisfied."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mode(path: Path) -> str:
    return oct(stat.S_IMODE(path.stat().st_mode))


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
    result = renameat2(
        AT_FDCWD,
        os.fsencode(source),
        AT_FDCWD,
        os.fsencode(target),
        RENAME_NOREPLACE,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise ArchiveError(
            f"No-replace archive move failed: {os.strerror(error_number)}"
        )


def write_exclusive_readonly_json(path: Path, payload: dict[str, Any]) -> None:
    encoded = (
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o444)
    try:
        offset = 0
        while offset < len(encoded):
            written = os.write(descriptor, encoded[offset:])
            if written <= 0:
                raise ArchiveError(f"Unable to write archive record: {path}")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
    finally:
        os.close(descriptor)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(
        path,
        os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0),
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def preflight(entry: dict[str, Any]) -> dict[str, Any]:
    source = entry["source"]
    archive = entry["archive"]
    target = entry["target"]
    if not source.is_dir() or source.is_symlink():
        raise ArchiveError(f"Source key directory is missing or unsafe: {source}")
    if os.path.lexists(archive):
        raise ArchiveError(f"Archive destination already exists: {archive}")
    for formal_path in (
        target,
        Path(str(target) + ".ed25519.sig"),
        Path(str(target) + ".ed25519.sig.pending"),
    ):
        if os.path.lexists(formal_path):
            raise ArchiveError(f"Formal slot is not empty: {formal_path}")
    observed_names = sorted(path.name for path in source.iterdir())
    expected_names = sorted((PRIVATE_FILENAME, PUBLIC_FILENAME))
    if observed_names != expected_names:
        raise ArchiveError(
            f"Unexpected key-directory contents for {entry['role']}: {observed_names}"
        )
    public_key = source / PUBLIC_FILENAME
    private_key = source / PRIVATE_FILENAME
    if (
        mode(source) != "0o700"
        or mode(public_key) != "0o444"
        or mode(private_key) != "0o400"
    ):
        raise ArchiveError(f"Unexpected key modes for {entry['role']}")
    observed_public_sha256 = sha256(public_key)
    if observed_public_sha256 != entry["public_sha256"]:
        raise ArchiveError(f"Public-key identity drifted for {entry['role']}")
    return {
        "role": entry["role"],
        "source": str(source),
        "archive": str(archive),
        "target": str(target),
        "signer_version": "v5",
        "signature_created": False,
        "formal_slots_absent": True,
        "reason": (
            "UNUSED_AFTER_INTERNAL_REVIEW_FOUND_POST_RETIREMENT_"
            "TERMINAL_SIGNED_TUPLE_FALSE_SUCCESS_WINDOW"
        ),
        "source_directory_mode_before_archive": mode(source),
        "public_key": {
            "filename": PUBLIC_FILENAME,
            "mode": mode(public_key),
            "size_bytes": public_key.stat().st_size,
            "sha256": observed_public_sha256,
        },
        "private_key": {
            "filename": PRIVATE_FILENAME,
            "mode": mode(private_key),
            "size_bytes": private_key.stat().st_size,
            "sha256": sha256(private_key),
            "retired_by_signing": False,
            "must_never_be_used": True,
        },
    }


def append_ledger(descriptor: int, payload: dict[str, Any]) -> None:
    encoded = (
        json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n"
    ).encode("utf-8")
    offset = 0
    while offset < len(encoded):
        written = os.write(descriptor, encoded[offset:])
        if written <= 0:
            raise ArchiveError("Unable to append archive ledger")
        offset += written
    os.fsync(descriptor)


def main() -> None:
    if os.path.lexists(LEDGER):
        raise ArchiveError(f"Ledger already exists: {LEDGER}")
    if (
        not PRIOR_ARCHIVE_LEDGER.is_file()
        or mode(PRIOR_ARCHIVE_LEDGER) != "0o444"
        or sha256(PRIOR_ARCHIVE_LEDGER) != PRIOR_ARCHIVE_LEDGER_SHA256
    ):
        raise ArchiveError("Prior archive ledger identity drifted")
    records = [preflight(dict(entry)) for entry in ENTRIES]
    if len(records) != len(ENTRIES):
        raise ArchiveError("Preflight record cardinality mismatch")
    for entry in ENTRIES:
        archive_parent = entry["archive"].parent
        if os.path.lexists(archive_parent):
            if not archive_parent.is_dir() or archive_parent.is_symlink():
                raise ArchiveError(f"Unsafe archive parent: {archive_parent}")
        else:
            os.mkdir(archive_parent, 0o700)
            fsync_directory(archive_parent.parent)

    ledger_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        ledger_flags |= os.O_NOFOLLOW
    ledger_fd = os.open(LEDGER, ledger_flags, 0o600)
    completed: list[dict[str, Any]] = []
    try:
        append_ledger(
            ledger_fd,
            {
                "event": "PREFLIGHT_COMPLETE",
                "entry_count": len(records),
                "roles": [record["role"] for record in records],
                "prior_archive_batch": {
                    "path": str(PRIOR_ARCHIVE_LEDGER),
                    "sha256": PRIOR_ARCHIVE_LEDGER_SHA256,
                    "status": "PASS",
                    "signer_version": "v4",
                },
            },
        )
        for entry, record in zip(ENTRIES, records):
            source = entry["source"]
            archive = entry["archive"]
            write_exclusive_readonly_json(source / STATUS_FILENAME, record)
            fsync_directory(source)
            append_ledger(
                ledger_fd,
                {
                    "event": "ARCHIVE_MOVE_START",
                    "role": record["role"],
                    "source": str(source),
                    "archive": str(archive),
                },
            )
            rename_no_replace(source, archive)
            os.chmod(archive, 0o700)
            fsync_directory(archive)
            fsync_directory(archive.parent)
            fsync_directory(source.parent)
            if os.path.lexists(source) or not archive.is_dir():
                raise ArchiveError(f"Post-move path verification failed: {record['role']}")
            if mode(archive) != "0o700" or mode(archive / STATUS_FILENAME) != "0o444":
                raise ArchiveError(f"Post-move mode verification failed: {record['role']}")
            completed_record = {
                **record,
                "archive_directory_mode": mode(archive),
                "archive_record": str(archive / STATUS_FILENAME),
                "status": "ARCHIVED_UNUSED_NEVER_SIGNED",
            }
            completed.append(completed_record)
            append_ledger(
                ledger_fd,
                {
                    "event": "ARCHIVE_MOVE_COMPLETE",
                    **completed_record,
                },
            )
        append_ledger(
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
        fsync_directory(LEDGER.parent)
    print(
        json.dumps(
            {
                "pass": len(completed) == len(ENTRIES),
                "ledger": str(LEDGER),
                "archived": completed,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
