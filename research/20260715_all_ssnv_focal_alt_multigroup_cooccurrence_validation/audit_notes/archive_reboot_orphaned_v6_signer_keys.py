#!/usr/bin/env python3
"""Archive reboot-orphaned unused signer keys without deletion or overwrite."""

from __future__ import annotations

import ctypes
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Any


HERE = Path(__file__).resolve().parent
LEDGER = HERE / "20260719_reboot_orphaned_v6_signer_key_archive.v1.jsonl"
PUBLIC_FILENAME = "ed25519_public.pem"
PRIVATE_FILENAME = "ed25519_private_encrypted_unrecoverable_after_signing.pem"
STATUS_FILENAME = "UNUSED_KEY_ARCHIVE_RECORD.v1.json"
AT_FDCWD = -100
RENAME_NOREPLACE = 1

TOPIC = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)

ENTRIES = (
    {
        "role": "source",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
            "20260719_all_ssnv_v8_command_parity_bootstrap"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
            "archive/20260719_all_ssnv_v8_command_parity_bootstrap_"
            "reboot_orphaned_unused_01"
        ),
        "target": Path(
            "/big7_disk/liaoyoyo2001/InterSubMod/docs/provenance/"
            "source_authorities/"
            "20260719_all_ssnv_focal_alt_release_source_authority.v6.approval.json"
        ),
        "public_sha256": (
            "f19961cd55c144e125c2e72622e87b77942346267fc62c4cc21e3a857e0292fa"
        ),
    },
    {
        "role": "result",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_result_release_authority/"
            "20260719_all_ssnv_result_v4_bound_bootstrap"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_result_release_authority/archive/"
            "20260719_all_ssnv_result_v4_bound_bootstrap_"
            "reboot_orphaned_unused_01"
        ),
        "target": TOPIC / "results/task_b_final_dataset_release_receipt.v1.json",
        "public_sha256": (
            "3779ee723bb960a4a1ca684dc92edb9da32b7b0faddc64258031abc0077ec116"
        ),
    },
    {
        "role": "report",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_report_release_authority/"
            "20260719_all_ssnv_report_v4_bound_bootstrap"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_report_release_authority/archive/"
            "20260719_all_ssnv_report_v4_bound_bootstrap_"
            "reboot_orphaned_unused_01"
        ),
        "target": TOPIC / "results/task_b_final_report_release_receipt.v1.json",
        "public_sha256": (
            "a43458292c04f131c872e31276a181305657b22c6bdf23e4272687e5328bf9bc"
        ),
    },
    {
        "role": "supplemental",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_supplemental_report_integrity/"
            "20260719_positional_singleton_v3_bound_bootstrap"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_supplemental_report_integrity/archive/"
            "20260719_positional_singleton_v3_bound_bootstrap_"
            "reboot_orphaned_unused_01"
        ),
        "target": (
            TOPIC
            / "results/positional_singleton_supplemental_release_receipt.v1.json"
        ),
        "public_sha256": (
            "21337cfc649e7fd4e344af6585cbe22b59059c7127f60bbf1876ac60d83ff63d"
        ),
    },
    {
        "role": "reviewer_A",
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_external_review_attestation/"
            "20260719_reviewer_A_v20_command_parity_bootstrap"
        ),
        "archive": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_external_review_attestation/archive/"
            "20260719_reviewer_A_v20_command_parity_bootstrap_"
            "reboot_orphaned_unused_01"
        ),
        "target": (
            TOPIC
            / "audit_notes/"
            "20260719_external_claude_reviewer_A_v19_"
            "command_parity_process_attestation.v1.json"
        ),
        "public_sha256": (
            "522928c52f1dad8556a9d49f4de89e01eb22eace6e49378b6159aeed0eeb5bb0"
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


def boot_time_epoch() -> int:
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith("btime "):
            return int(line.split()[1])
    raise ArchiveError("Unable to determine host boot time")


def matching_processes(key_directory: Path) -> list[int]:
    needle = os.fsencode(str(key_directory))
    matches: list[int] = []
    for process_dir in Path("/proc").iterdir():
        if not process_dir.name.isdigit():
            continue
        try:
            command_line = (process_dir / "cmdline").read_bytes()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if needle in command_line:
            matches.append(int(process_dir.name))
    return sorted(matches)


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


def write_exclusive_json(path: Path, payload: dict[str, Any], file_mode: int) -> None:
    encoded = (
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, file_mode)
    try:
        offset = 0
        while offset < len(encoded):
            written = os.write(descriptor, encoded[offset:])
            if written <= 0:
                raise ArchiveError(f"Unable to write {path}")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, file_mode)
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


def preflight(entry: dict[str, Any], boot_epoch: int) -> dict[str, Any]:
    source = entry["source"]
    archive = entry["archive"]
    target = entry["target"]
    if not source.is_dir() or source.is_symlink():
        raise ArchiveError(f"Missing or unsafe source key directory: {source}")
    if os.path.lexists(archive):
        raise ArchiveError(f"Archive destination already exists: {archive}")
    for formal_path in (
        target,
        Path(str(target) + ".ed25519.sig"),
        Path(str(target) + ".ed25519.sig.pending"),
    ):
        if os.path.lexists(formal_path):
            raise ArchiveError(f"Formal slot is not empty: {formal_path}")
    expected_names = sorted((PRIVATE_FILENAME, PUBLIC_FILENAME))
    observed_names = sorted(path.name for path in source.iterdir())
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
    if sha256(public_key) != entry["public_sha256"]:
        raise ArchiveError(f"Public-key identity drifted for {entry['role']}")
    if (
        public_key.stat().st_mtime >= boot_epoch
        or private_key.stat().st_mtime >= boot_epoch
    ):
        raise ArchiveError(f"Key material does not predate boot: {entry['role']}")
    live_processes = matching_processes(source)
    if live_processes:
        raise ArchiveError(
            f"Signer process still references {entry['role']}: {live_processes}"
        )
    return {
        "role": entry["role"],
        "source": str(source),
        "archive": str(archive),
        "target": str(target),
        "signer_version": "v6",
        "signature_created": False,
        "formal_slots_absent": True,
        "host_boot_time_epoch": boot_epoch,
        "key_material_predates_boot": True,
        "live_signer_processes": [],
        "reason": "HOST_REBOOT_ORPHANED_IN_MEMORY_PASSPHRASE",
        "status": "ARCHIVED_UNUSED_NEVER_SIGNED",
        "source_directory_mode_before_archive": mode(source),
        "public_key": {
            "filename": PUBLIC_FILENAME,
            "mode": mode(public_key),
            "size_bytes": public_key.stat().st_size,
            "sha256": sha256(public_key),
        },
        "private_key": {
            "filename": PRIVATE_FILENAME,
            "mode": mode(private_key),
            "size_bytes": private_key.stat().st_size,
            "sha256": sha256(private_key),
            "retired_by_signing": False,
            "passphrase_recoverable": False,
            "must_never_be_used": True,
        },
    }


def append_jsonl(descriptor: int, payload: dict[str, Any]) -> None:
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
    boot_epoch = boot_time_epoch()
    records = [preflight(dict(entry), boot_epoch) for entry in ENTRIES]
    for entry in ENTRIES:
        archive_parent = entry["archive"].parent
        if os.path.lexists(archive_parent):
            if not archive_parent.is_dir() or archive_parent.is_symlink():
                raise ArchiveError(f"Unsafe archive parent: {archive_parent}")
        else:
            os.mkdir(archive_parent, 0o700)
            fsync_directory(archive_parent.parent)

    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    ledger_fd = os.open(LEDGER, flags, 0o600)
    completed: list[str] = []
    try:
        append_jsonl(
            ledger_fd,
            {
                "event": "PREFLIGHT_COMPLETE",
                "entry_count": len(records),
                "host_boot_time_epoch": boot_epoch,
                "roles": [record["role"] for record in records],
            },
        )
        for entry, record in zip(ENTRIES, records):
            source = entry["source"]
            archive = entry["archive"]
            write_exclusive_json(source / STATUS_FILENAME, record, 0o444)
            fsync_directory(source)
            append_jsonl(
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
                raise ArchiveError(f"Post-move verification failed: {record['role']}")
            if mode(archive) != "0o700":
                raise ArchiveError(f"Archive mode drifted: {record['role']}")
            completed.append(record["role"])
            append_jsonl(
                ledger_fd,
                {
                    "event": "ARCHIVE_MOVE_COMPLETE",
                    "role": record["role"],
                    "archive": str(archive),
                    "record_sha256": sha256(archive / STATUS_FILENAME),
                    "status": "ARCHIVED_UNUSED_NEVER_SIGNED",
                },
            )
        if len(completed) != len(ENTRIES):
            raise ArchiveError("Archive completion cardinality mismatch")
        append_jsonl(
            ledger_fd,
            {
                "event": "ARCHIVE_BATCH_COMPLETE",
                "completed_roles": completed,
                "entry_count": len(completed),
                "status": "PASS",
            },
        )
        os.fchmod(ledger_fd, 0o444)
        os.fsync(ledger_fd)
    finally:
        os.close(ledger_fd)
    print(
        json.dumps(
            {
                "status": "PASS",
                "ledger": str(LEDGER),
                "roles": completed,
                "boot_time_epoch": boot_epoch,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
