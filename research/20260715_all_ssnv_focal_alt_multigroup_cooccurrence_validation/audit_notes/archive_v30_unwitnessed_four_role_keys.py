#!/usr/bin/env python3
"""Archive unused v30 keys before a witnessed four-role rebootstrap."""

from __future__ import annotations

from datetime import datetime, timezone
import ctypes
import fcntl
import hashlib
import json
import os
from pathlib import Path
import secrets
import signal
import stat
import sys
from typing import Any, Mapping


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
RESULT_ROOT = TOPIC_ROOT / "results"
REVIEW_ROOT = TOPIC_ROOT / "reviews"
PREPARE = AUDIT_ROOT / "20260724_v30_four_role_key_rebootstrap_archive_prepare.json"
RECEIPT = AUDIT_ROOT / "20260724_v30_four_role_key_rebootstrap_archive_receipt.json"
PROGRESS = AUDIT_ROOT / "20260724_v30_four_role_key_rebootstrap_archive_progress.jsonl"
PARTIAL = AUDIT_ROOT / "20260724_v30_four_role_key_rebootstrap_archive_partial.json"
OLD_BOOTSTRAP_ARTIFACTS = {
    AUDIT_ROOT / "20260724_v30_authority_terminal_key_bootstrap_prepare.json": (
        2_392,
        "db0e9db0c9d7024341ef2caa06b629358eef66958afddbcb0b7a8e7bde609d2a",
    ),
    AUDIT_ROOT / "20260724_v30_authority_terminal_key_bootstrap_receipt.json": (
        3_907,
        "7362ca6d77a9b924b0275d2683457543d667e80ae6b15e0861fee6e6dc12fbec",
    ),
    AUDIT_ROOT / "20260724_v30_authority_terminal_key_bootstrap_success.json": (
        2_620,
        "ab65a4dc0076300c474ed788c023230f76edb58ccd0ac0d3710a47aa94213257",
    ),
}
EXPECTED_SIGNER_PIDS = (3_580_365, 3_580_366)
AT_FDCWD = -100
RENAME_NOREPLACE = 1
LEASE_BREAK_EVENTS: list[int] = []

KEY_SPECS: dict[str, dict[str, Any]] = {
    "authority_v30_initial": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "20260724_v30"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
            "archive/20260724_v30_unused_pre_authority_four_role_rebootstrap_01"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "private_size": 119,
        "public_sha256": "208e7ed2f119e080806bbe6807ede8dd83556fcd16e9b0b4193043e75caf1236",
    },
    "terminal_v20_initial": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260724_m2v5_terminal_v20"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v20_unused_pre_authority_four_role_rebootstrap_01"
        ),
        "private_name": "ed25519_private_one_time_resident.pem",
        "private_size": 119,
        "public_sha256": "a69177e76b7943f43349a7c1c727e559b7103527427bcbc58d76ec3ca58a8c29",
    },
    "result_v7_initial": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "20260724_all_ssnv_result_v7_v30_recovery"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "archive/20260724_all_ssnv_result_v7_unused_pre_authority_four_role_rebootstrap_01"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "private_size": 290,
        "public_sha256": "6390a84edc65876be3c30f8432bc6efd4a1ee04f6ceaa963315d9d144e7a2430",
    },
    "report_v7_initial": {
        "source": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "20260724_all_ssnv_report_v7_v30_recovery"
        ),
        "destination": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "archive/20260724_all_ssnv_report_v7_unused_pre_authority_four_role_rebootstrap_01"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "private_size": 290,
        "public_sha256": "ca74b50e459f3a4c0bcdc0e6efaab97d753a40c74316e064b26c8cc208b7cf23",
    },
}

PROTECTED_OUTPUTS = (
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.bundle",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.json",
    RESULT_ROOT / "tumor_ref_promotion_schema_recovery_authority.v30.json.ed25519.sig",
    RESULT_ROOT / "tumor_ref_source_receipt_promotion_verification.recovery.v30.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.json",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.log",
    RESULT_ROOT / "m2v5_runner_only_gate_replay.recovery.v30.success_witness.json",
    RESULT_ROOT / "m2v5_downstream_continuation.recovery.v30.json",
    RESULT_ROOT / "m2v5_downstream_continuation_exit_attestation.recovery.v30.json",
    RESULT_ROOT
    / "m2v5_downstream_continuation_exit_attestation.recovery.v30.json.ed25519.sig",
    RESULT_ROOT / "m2v5_downstream_continuation_supervisor_success.recovery.v30.json",
    RESULT_ROOT / "m2v5_downstream_continuation_incident.recovery.v30.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json.ed25519.sig",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
    RESULT_ROOT / "task_b_final_report_release_receipt.v1.json.ed25519.sig",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v30.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v30.json",
    REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_external_claude_opus.v30.json",
)


class ArchiveError(RuntimeError):
    """Raised before any overwrite or ambiguous key retirement is accepted."""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_record(path: Path, *, expected_mode: int = 0o444) -> dict[str, Any]:
    descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        before = os.fstat(descriptor)
        data = b"".join(
            os.pread(descriptor, min(1024 * 1024, before.st_size - offset), offset)
            for offset in range(0, before.st_size, 1024 * 1024)
        )
        after = os.fstat(descriptor)
        live = os.lstat(path)
        identity = lambda item: (
            item.st_dev,
            item.st_ino,
            item.st_mode,
            item.st_size,
            item.st_mtime_ns,
            item.st_ctime_ns,
            item.st_nlink,
        )
        if (
            not stat.S_ISREG(before.st_mode)
            or len(data) != before.st_size
            or identity(before) != identity(after)
            or identity(after) != identity(live)
            or stat.S_IMODE(after.st_mode) != expected_mode
            or after.st_nlink != 1
            or path.resolve(strict=True) != path
        ):
            raise ArchiveError(f"Artifact identity/mode drift: {path}")
        return {
            "path": str(path),
            "size_bytes": len(data),
            "sha256": sha256_bytes(data),
            "mode": oct(expected_mode),
            "link_count": 1,
        }
    finally:
        os.close(descriptor)


def private_metadata(path: Path, expected_size: int) -> dict[str, Any]:
    descriptor = os.open(
        path, os.O_PATH | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        observed = os.fstat(descriptor)
        live = os.lstat(path)
        if (
            not stat.S_ISREG(observed.st_mode)
            or (observed.st_dev, observed.st_ino, observed.st_mode, observed.st_size)
            != (live.st_dev, live.st_ino, live.st_mode, live.st_size)
            or stat.S_IMODE(observed.st_mode) != 0o400
            or observed.st_size != expected_size
            or observed.st_nlink != 1
            or path.resolve(strict=True) != path
        ):
            raise ArchiveError(f"Private-key metadata drift: {path}")
        return {
            "path": str(path),
            "size_bytes": observed.st_size,
            "mode": "0o400",
            "link_count": observed.st_nlink,
        }
    finally:
        os.close(descriptor)


def write_exclusive(path: Path, data: bytes) -> dict[str, Any]:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
        0o600,
    )
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise ArchiveError(f"Short write: {path}")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return file_record(path)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def strict_json_bytes(data: bytes, label: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise ArchiveError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    value = json.loads(
        data,
        object_pairs_hook=no_duplicates,
        parse_constant=lambda item: (_ for _ in ()).throw(
            ArchiveError(f"Non-finite JSON in {label}: {item}")
        ),
    )
    if not isinstance(value, dict):
        raise ArchiveError(f"JSON root is not an object in {label}")
    return value


def publish_json_atomic(
    path: Path,
    payload: dict[str, Any],
    *,
    transaction_id: str,
    stage: str,
) -> dict[str, Any]:
    staging = path.with_name(f".{path.name}.staging.{transaction_id}.{stage}")
    if os.path.lexists(path) or os.path.lexists(staging):
        raise ArchiveError(f"Archive publication path occupied: {path}")
    data = json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True).encode(
        "ascii"
    ) + b"\n"
    staging_record = write_exclusive(staging, data)
    if strict_json_bytes(staging.read_bytes(), f"staged {stage}") != payload:
        raise ArchiveError(f"Staged archive JSON round-trip drift: {stage}")
    parent_fd = os.open(
        path.parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
    )
    try:
        parent_before = os.fstat(parent_fd)
        move_no_replace(staging.name, path.name, parent_fd, parent_fd)
        os.fsync(parent_fd)
        parent_after = os.fstat(parent_fd)
        if (
            parent_before.st_dev,
            parent_before.st_ino,
            parent_before.st_mode,
        ) != (
            parent_after.st_dev,
            parent_after.st_ino,
            parent_after.st_mode,
        ):
            raise ArchiveError("Archive publication parent identity changed")
    finally:
        os.close(parent_fd)
    final_record = file_record(path)
    if (
        final_record["size_bytes"] != staging_record["size_bytes"]
        or final_record["sha256"] != staging_record["sha256"]
        or strict_json_bytes(path.read_bytes(), f"published {stage}") != payload
    ):
        raise ArchiveError(f"Published archive JSON round-trip drift: {stage}")
    return final_record


class ProgressLedger:
    def __init__(self, path: Path, transaction_id: str) -> None:
        self.path = path
        self.transaction_id = transaction_id
        self.fd = os.open(
            path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
            0o600,
        )
        self.closed = False

    def append(self, event: str, **details: Any) -> None:
        if self.closed:
            raise ArchiveError("Progress ledger is already closed")
        record = {
            "created_at_utc": now_utc(),
            "transaction_id": self.transaction_id,
            "event": event,
            **details,
        }
        data = json.dumps(record, ensure_ascii=True, sort_keys=True).encode("ascii") + b"\n"
        offset = 0
        while offset < len(data):
            written = os.write(self.fd, data[offset:])
            if written <= 0:
                raise ArchiveError("Short progress-ledger write")
            offset += written
        os.fsync(self.fd)

    def finish(self) -> dict[str, Any]:
        if not self.closed:
            os.fchmod(self.fd, 0o444)
            os.fsync(self.fd)
            os.close(self.fd)
            self.closed = True
            fsync_directory(self.path.parent)
        return file_record(self.path)

    def close_best_effort(self) -> None:
        if self.closed:
            return
        try:
            os.fchmod(self.fd, 0o444)
            os.fsync(self.fd)
        finally:
            os.close(self.fd)
            self.closed = True


def move_no_replace(
    source_name: str,
    destination_name: str,
    source_parent_fd: int,
    destination_parent_fd: int,
) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise ArchiveError("renameat2 is required for no-replace archival")
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        source_parent_fd,
        os.fsencode(source_name),
        destination_parent_fd,
        os.fsencode(destination_name),
        RENAME_NOREPLACE,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise ArchiveError(
            f"No-replace archive rename failed: {source_name} -> {destination_name}: "
            f"{os.strerror(error_number)}"
        )


def live_key_roots() -> dict[str, Path]:
    roots: dict[str, Path] = {}
    for role, spec in KEY_SPECS.items():
        source = Path(spec["source"])
        destination = Path(spec["destination"])
        present = [path for path in (source, destination) if os.path.lexists(path)]
        if len(present) != 1:
            raise ArchiveError(
                f"Expected exactly one live source/archive root for {role}: {present}"
            )
        roots[role] = present[0]
    return roots


def record_restricted_procfs_fd_visibility(process: Path) -> dict[str, Any]:
    pid = int(process.name)
    comm = (process / "comm").read_bytes()
    cmdline = (process / "cmdline").read_bytes()
    cgroup = (process / "cgroup").read_bytes()
    stat_line = (process / "stat").read_text(encoding="ascii").strip()
    try:
        _, stat_tail = stat_line.rsplit(") ", 1)
        stat_fields = stat_tail.split()
        parent_pid = int(stat_fields[1])
        starttime_ticks = int(stat_fields[19])
    except (ValueError, IndexError) as error:
        raise ArchiveError(f"Cannot parse restricted PID stat identity: {pid}") from error
    fd_directory = process / "fd"
    process_stat = process.stat()
    fd_stat = os.lstat(fd_directory)
    live_roots = live_key_roots()
    key_needles = [os.fsencode(path) for path in live_roots.values()]
    if (
        process_stat.st_uid != os.getuid()
        or not stat.S_ISDIR(fd_stat.st_mode)
        or any(needle in cmdline for needle in key_needles)
    ):
        raise ArchiveError(f"Restricted-procfs process evidence drift: {pid}")
    return {
        "pid": pid,
        "process_uid": process_stat.st_uid,
        "comm_sha256": sha256_bytes(comm),
        "cmdline_sha256": sha256_bytes(cmdline),
        "cgroup_sha256": sha256_bytes(cgroup),
        "starttime_ticks": starttime_ticks,
        "parent_pid": parent_pid,
        "fd_directory_mode": oct(stat.S_IMODE(fd_stat.st_mode)),
        "fd_directory_uid": fd_stat.st_uid,
        "cmdline_contains_no_live_key_root": True,
        "fd_enumeration_available": False,
        "kernel_write_lease_required_before_key_move": True,
    }


def active_key_processes() -> list[dict[str, Any]]:
    needles = {
        role: os.fsencode(path) for role, path in live_key_roots().items()
    }
    found: list[dict[str, Any]] = []
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        if int(process.name) == os.getpid():
            continue
        try:
            process_stat = process.stat()
            if process_stat.st_uid != os.getuid():
                continue
            cmdline = (process / "cmdline").read_bytes()
        except (FileNotFoundError, ProcessLookupError):
            continue
        except PermissionError as error:
            raise ArchiveError(
                f"Cannot inspect same-UID process command line: {process.name}"
            ) from error
        matches = {
            role: os.fsdecode(value)
            for role, value in needles.items()
            if value in cmdline
        }
        if matches:
            found.append({"pid": int(process.name), "matched_roots": matches})
    return found


def active_key_descriptor_scan() -> dict[str, Any]:
    protected: dict[tuple[int, int], str] = {}
    for role, root in live_key_roots().items():
        spec = KEY_SPECS[role]
        for path in (
            root,
            root / "ed25519_public.pem",
            root / str(spec["private_name"]),
        ):
            opened = os.lstat(path)
            protected[(opened.st_dev, opened.st_ino)] = str(path)
    found: list[dict[str, Any]] = []
    direct_scan_pids: list[int] = []
    restricted_procfs_visibility: list[dict[str, Any]] = []
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        if int(process.name) == os.getpid():
            continue
        descriptor_root = process / "fd"
        try:
            process_stat = process.stat()
            if process_stat.st_uid != os.getuid():
                continue
            descriptors = list(descriptor_root.iterdir())
        except (FileNotFoundError, ProcessLookupError):
            continue
        except PermissionError:
            restricted_procfs_visibility.append(
                record_restricted_procfs_fd_visibility(process)
            )
            continue
        direct_scan_pids.append(int(process.name))
        for descriptor in descriptors:
            try:
                opened = descriptor.stat()
            except (FileNotFoundError, ProcessLookupError):
                continue
            except PermissionError as error:
                raise ArchiveError(
                    f"Cannot inspect same-UID descriptor: {process.name}/{descriptor.name}"
                ) from error
            path = protected.get((opened.st_dev, opened.st_ino))
            if path is not None:
                found.append(
                    {"pid": int(process.name), "fd": descriptor.name, "path": path}
                )
    return {
        "matches": found,
        "direct_scan_pids": sorted(direct_scan_pids),
        "restricted_procfs_visibility": restricted_procfs_visibility,
        "procfs_fd_visibility_complete": not restricted_procfs_visibility,
        "kernel_write_lease_gate_required_before_each_key_move": bool(
            restricted_procfs_visibility
        ),
        "pass": not found,
    }


def active_key_descriptors() -> list[dict[str, Any]]:
    return active_key_descriptor_scan()["matches"]


def lease_break_handler(signum: int, _frame: Any) -> None:
    LEASE_BREAK_EVENTS.append(signum)


def acquire_key_write_leases(
    role: str,
    descriptors: Mapping[str, tuple[int, Path]],
) -> dict[str, Any]:
    event_offset = len(LEASE_BREAK_EVENTS)
    acquired: list[str] = []
    records: dict[str, Any] = {}
    try:
        for label, (descriptor, path) in descriptors.items():
            try:
                fcntl.fcntl(descriptor, fcntl.F_SETLEASE, fcntl.F_WRLCK)
            except OSError as error:
                raise ArchiveError(
                    f"Kernel write lease rejected an existing opener: {role}/{label}: "
                    f"errno={error.errno}"
                ) from error
            acquired.append(label)
            observed = os.fstat(descriptor)
            lease_state = fcntl.fcntl(descriptor, fcntl.F_GETLEASE)
            if lease_state != fcntl.F_WRLCK:
                raise ArchiveError(f"Kernel write lease state drift: {role}/{label}")
            records[label] = {
                "path_before_move": str(path),
                "device": observed.st_dev,
                "inode": observed.st_ino,
                "mode": oct(stat.S_IMODE(observed.st_mode)),
                "lease_state": "F_WRLCK",
            }
        if len(LEASE_BREAK_EVENTS) != event_offset:
            raise ArchiveError(f"Kernel write lease break requested during acquire: {role}")
        return {
            "role": role,
            "event_offset": event_offset,
            "files": records,
            "existing_external_openers_rejected_by_kernel": True,
            "new_external_opens_blocked_until_post_move_verification": True,
            "pass": True,
        }
    except BaseException:
        for label in reversed(acquired):
            descriptor = descriptors[label][0]
            try:
                fcntl.fcntl(descriptor, fcntl.F_SETLEASE, fcntl.F_UNLCK)
            except OSError:
                pass
        raise


def verify_key_write_leases(
    role: str,
    descriptors: Mapping[str, tuple[int, Path]],
    evidence: Mapping[str, Any],
    *,
    stage: str,
) -> None:
    if len(LEASE_BREAK_EVENTS) != evidence.get("event_offset"):
        raise ArchiveError(f"Kernel write lease break requested {stage}: {role}")
    for label, (descriptor, _path) in descriptors.items():
        observed = os.fstat(descriptor)
        expected = evidence["files"][label]
        if (
            fcntl.fcntl(descriptor, fcntl.F_GETLEASE) != fcntl.F_WRLCK
            or observed.st_dev != expected["device"]
            or observed.st_ino != expected["inode"]
        ):
            raise ArchiveError(f"Kernel write lease/inode drift {stage}: {role}/{label}")


def release_key_write_leases(
    descriptors: Mapping[str, tuple[int, Path]],
) -> None:
    for descriptor, _path in descriptors.values():
        fcntl.fcntl(descriptor, fcntl.F_SETLEASE, fcntl.F_UNLCK)


def preflight() -> dict[str, Any]:
    if any(Path(f"/proc/{pid}").exists() for pid in EXPECTED_SIGNER_PIDS):
        raise ArchiveError("Result/report signer PID remains active")
    occupied = [
        str(path)
        for path in (*PROTECTED_OUTPUTS, PREPARE, RECEIPT, PROGRESS, PARTIAL)
        if os.path.lexists(path)
    ]
    if occupied:
        raise ArchiveError(f"Protected output already exists: {occupied}")
    process_matches = active_key_processes()
    descriptor_scan = active_key_descriptor_scan()
    descriptor_matches = descriptor_scan["matches"]
    if process_matches or descriptor_matches:
        raise ArchiveError(
            f"Key root remains held: cmdline={process_matches} fds={descriptor_matches}"
        )
    old_bootstrap = {}
    for path, (expected_size, expected_sha256) in OLD_BOOTSTRAP_ARTIFACTS.items():
        record = file_record(path)
        if record["size_bytes"] != expected_size or record["sha256"] != expected_sha256:
            raise ArchiveError(f"Initial bootstrap artifact drift: {path}")
        old_bootstrap[str(path)] = record
    public_hashes = set()
    key_inputs = {}
    for role, spec in KEY_SPECS.items():
        source = Path(spec["source"])
        destination = Path(spec["destination"])
        if (
            not source.is_dir()
            or source.resolve(strict=True) != source
            or stat.S_IMODE(os.lstat(source).st_mode) != 0o700
            or os.path.lexists(destination)
            or destination.parent.resolve(strict=True) != destination.parent
        ):
            raise ArchiveError(f"Invalid key source/destination: {role}")
        public = file_record(source / "ed25519_public.pem")
        private = private_metadata(
            source / str(spec["private_name"]), int(spec["private_size"])
        )
        if public["sha256"] != spec["public_sha256"] or public["size_bytes"] != 113:
            raise ArchiveError(f"Public-key identity drift: {role}")
        public_hashes.add(public["sha256"])
        key_inputs[role] = {"public_key": public, "private_key_metadata_only": private}
    if len(public_hashes) != 4:
        raise ArchiveError("Initial four role public keys are not pairwise distinct")
    return {
        "expected_signer_pids": list(EXPECTED_SIGNER_PIDS),
        "expected_signer_pids_absent": True,
        "active_key_processes": process_matches,
        "active_key_descriptors": descriptor_matches,
        "active_key_descriptor_scan": descriptor_scan,
        "old_bootstrap_artifacts": old_bootstrap,
        "key_inputs": key_inputs,
        "protected_outputs_absent": True,
        "private_key_bytes_read": False,
        "pass": True,
    }


def archive_key(
    role: str,
    spec: Mapping[str, Any],
    *,
    transaction_id: str,
    ledger: ProgressLedger,
) -> dict[str, Any]:
    source = Path(spec["source"])
    destination = Path(spec["destination"])
    source_parent_fd = destination_parent_fd = root_fd = -1
    public_fd = private_fd = -1
    lease_descriptors: dict[str, tuple[int, Path]] = {}
    leases_held = False
    try:
        source_parent_fd = os.open(
            source.parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        )
        destination_parent_fd = os.open(
            destination.parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        )
        root_fd = os.open(
            source.name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=source_parent_fd,
        )
        public_fd = os.open(
            "ed25519_public.pem",
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=root_fd,
        )
        private_fd = os.open(
            str(spec["private_name"]),
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=root_fd,
        )
        root_opened = os.fstat(root_fd)
        source_from_parent = os.stat(
            source.name, dir_fd=source_parent_fd, follow_symlinks=False
        )
        if (
            root_opened.st_dev,
            root_opened.st_ino,
            root_opened.st_mode,
            root_opened.st_nlink,
        ) != (
            source_from_parent.st_dev,
            source_from_parent.st_ino,
            source_from_parent.st_mode,
            source_from_parent.st_nlink,
        ):
            raise ArchiveError(f"Bound source root identity drift: {role}")
        public = file_record(source / "ed25519_public.pem")
        private = private_metadata(
            source / str(spec["private_name"]), int(spec["private_size"])
        )
        ledger.append(
            "BEFORE_MUTATION",
            role=role,
            source=str(source),
            destination=str(destination),
            source_root_device=root_opened.st_dev,
            source_root_inode=root_opened.st_ino,
        )
        record = {
            "schema_name": "intersubmod.pre_authority_key_archive",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "transaction_id": transaction_id,
            "role": role,
            "state": "UNUSED_PRE_AUTHORITY_SUPERSEDED_BY_WITNESSED_FOUR_ROLE_REBOOTSTRAP",
            "source_directory": str(source),
            "archive_directory": str(destination),
            "public_key": public,
            "private_key_metadata_only": private,
            "private_key_bytes_read": False,
            "key_reuse_forbidden": True,
            "pass": True,
        }
        record_data = json.dumps(
            record, ensure_ascii=True, indent=2, sort_keys=True
        ).encode("ascii") + b"\n"
        archive_fd = os.open(
            "FAILED_KEY_ARCHIVE_RECORD.v1.json",
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
            0o600,
            dir_fd=root_fd,
        )
        try:
            offset = 0
            while offset < len(record_data):
                written = os.write(archive_fd, record_data[offset:])
                if written <= 0:
                    raise ArchiveError(f"Short key archive-record write: {role}")
                offset += written
            os.fsync(archive_fd)
            os.fchmod(archive_fd, 0o444)
            os.fsync(archive_fd)
            archive_stat = os.fstat(archive_fd)
        finally:
            os.close(archive_fd)
        os.fsync(root_fd)
        archive_record = {
            "path": str(destination / "FAILED_KEY_ARCHIVE_RECORD.v1.json"),
            "size_bytes": archive_stat.st_size,
            "sha256": sha256_bytes(record_data),
            "mode": "0o444",
            "link_count": archive_stat.st_nlink,
        }
        ledger.append(
            "ARCHIVE_RECORD_WRITTEN",
            role=role,
            archive_record=archive_record,
        )
        process_matches = active_key_processes()
        descriptor_matches = active_key_descriptors()
        if process_matches or descriptor_matches:
            raise ArchiveError(
                f"Key became held before move {role}: "
                f"cmdline={process_matches} fds={descriptor_matches}"
            )
        lease_descriptors = {
            "public_key": (public_fd, source / "ed25519_public.pem"),
            "private_key": (
                private_fd,
                source / str(spec["private_name"]),
            ),
        }
        lease_evidence = acquire_key_write_leases(role, lease_descriptors)
        leases_held = True
        verify_key_write_leases(
            role,
            lease_descriptors,
            lease_evidence,
            stage="before move",
        )
        ledger.append("KERNEL_WRITE_LEASE_GATE_ACQUIRED", evidence=lease_evidence)
        ledger.append("BEFORE_MOVE", role=role, kernel_write_leases_verified=True)
        move_no_replace(
            source.name,
            destination.name,
            source_parent_fd,
            destination_parent_fd,
        )
        os.fsync(source_parent_fd)
        os.fsync(destination_parent_fd)
        moved_from_parent = os.stat(
            destination.name,
            dir_fd=destination_parent_fd,
            follow_symlinks=False,
        )
        if (
            root_opened.st_dev,
            root_opened.st_ino,
            root_opened.st_mode,
            root_opened.st_nlink,
        ) != (
            moved_from_parent.st_dev,
            moved_from_parent.st_ino,
            moved_from_parent.st_mode,
            moved_from_parent.st_nlink,
        ) or os.path.lexists(source):
            raise ArchiveError(f"Post-archive path/inode state drift: {role}")
        verify_key_write_leases(
            role,
            lease_descriptors,
            lease_evidence,
            stage="after move",
        )
        release_key_write_leases(lease_descriptors)
        leases_held = False
        lease_evidence = {
            **lease_evidence,
            "verified_immediately_before_move": True,
            "verified_immediately_after_move": True,
            "released_only_after_post_move_verification": True,
        }
        ledger.append("KERNEL_WRITE_LEASE_GATE_RELEASED", evidence=lease_evidence)
        output = {
            "role": role,
            "archive_root": str(destination),
            "archive_record": archive_record,
            "public_key": {
                **public,
                "path": str(destination / "ed25519_public.pem"),
            },
            "private_key_metadata_only": {
                **private,
                "path": str(destination / str(spec["private_name"])),
            },
            "kernel_write_lease_gate": lease_evidence,
            "key_reuse_forbidden": True,
            "pass": True,
        }
        ledger.append("MOVED", role=role, archive=output)
        return output
    finally:
        if leases_held:
            try:
                release_key_write_leases(lease_descriptors)
            except OSError:
                pass
        for descriptor in (
            private_fd,
            public_fd,
            root_fd,
            destination_parent_fd,
            source_parent_fd,
        ):
            if descriptor >= 0:
                os.close(descriptor)


def validate_final_archives(
    archives: Mapping[str, Any], transaction_id: str
) -> dict[str, Any]:
    if set(archives) != set(KEY_SPECS):
        raise ArchiveError("Final archive role set is incomplete")
    validated: dict[str, Any] = {}
    for role, spec in KEY_SPECS.items():
        source = Path(spec["source"])
        destination = Path(spec["destination"])
        if (
            os.path.lexists(source)
            or not destination.is_dir()
            or destination.resolve(strict=True) != destination
            or stat.S_IMODE(os.lstat(destination).st_mode) != 0o700
            or set(path.name for path in destination.iterdir())
            != {
                "ed25519_public.pem",
                str(spec["private_name"]),
                "FAILED_KEY_ARCHIVE_RECORD.v1.json",
            }
        ):
            raise ArchiveError(f"Final key archive root state drift: {role}")
        archive_record = file_record(destination / "FAILED_KEY_ARCHIVE_RECORD.v1.json")
        public = file_record(destination / "ed25519_public.pem")
        private = private_metadata(
            destination / str(spec["private_name"]), int(spec["private_size"])
        )
        archive_payload = strict_json_bytes(
            (destination / "FAILED_KEY_ARCHIVE_RECORD.v1.json").read_bytes(),
            f"archive record {role}",
        )
        expected = archives[role]
        lease = expected.get("kernel_write_lease_gate", {})
        lease_files = lease.get("files", {}) if isinstance(lease, dict) else {}
        public_stat = os.lstat(destination / "ed25519_public.pem")
        private_stat = os.lstat(destination / str(spec["private_name"]))
        if (
            archive_payload.get("transaction_id") != transaction_id
            or archive_payload.get("role") != role
            or archive_payload.get("archive_directory") != str(destination)
            or archive_payload.get("key_reuse_forbidden") is not True
            or archive_payload.get("private_key_bytes_read") is not False
            or archive_payload.get("pass") is not True
            or archive_record != expected.get("archive_record")
            or public != expected.get("public_key")
            or private != expected.get("private_key_metadata_only")
            or set(lease_files) != {"public_key", "private_key"}
            or lease.get("existing_external_openers_rejected_by_kernel") is not True
            or lease.get("new_external_opens_blocked_until_post_move_verification")
            is not True
            or lease.get("verified_immediately_before_move") is not True
            or lease.get("verified_immediately_after_move") is not True
            or lease.get("released_only_after_post_move_verification") is not True
            or lease.get("pass") is not True
            or lease_files["public_key"].get("device") != public_stat.st_dev
            or lease_files["public_key"].get("inode") != public_stat.st_ino
            or lease_files["private_key"].get("device") != private_stat.st_dev
            or lease_files["private_key"].get("inode") != private_stat.st_ino
        ):
            raise ArchiveError(f"Final key archive evidence drift: {role}")
        validated[role] = {
            "archive_root": str(destination),
            "archive_record": archive_record,
            "public_key": public,
            "private_key_metadata_only": private,
            "kernel_write_lease_gate": lease,
            "key_reuse_forbidden": True,
            "pass": True,
        }
    return validated


def main() -> int:
    source = file_record(Path(__file__).resolve(strict=True))
    preflight_record = preflight()
    transaction_id = secrets.token_hex(16)
    prepare_payload = {
        "schema_name": "intersubmod.v30_four_role_key_rebootstrap_archive_prepare",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v30",
        "transaction_id": transaction_id,
        "source": source,
        "preflight": preflight_record,
        "intended_archives": {
            role: str(spec["destination"]) for role, spec in KEY_SPECS.items()
        },
        "status": "PREPARED_BEFORE_NO_REPLACE_ARCHIVAL",
        "release_authority_granted": False,
        "pass": False,
    }
    prepare_record = publish_json_atomic(
        PREPARE,
        prepare_payload,
        transaction_id=transaction_id,
        stage="prepare",
    )
    ledger = ProgressLedger(PROGRESS, transaction_id)
    archives: dict[str, Any] = {}
    progress_record: dict[str, Any] | None = None
    receipt_record: dict[str, Any] | None = None
    committed = False
    previous_mask: set[signal.Signals] | None = None
    previous_sigio_handler: Any = None
    blocked_signals = {
        signal.SIGINT,
        signal.SIGTERM,
        signal.SIGHUP,
        signal.SIGQUIT,
        signal.SIGUSR1,
        signal.SIGUSR2,
    }
    try:
        ledger.append("TRANSACTION_STARTED", prepare=prepare_record)
        current_mask = signal.pthread_sigmask(signal.SIG_BLOCK, set())
        if signal.SIGIO in current_mask:
            raise ArchiveError("SIGIO is blocked; kernel lease-break detection unavailable")
        LEASE_BREAK_EVENTS.clear()
        previous_sigio_handler = signal.getsignal(signal.SIGIO)
        signal.signal(signal.SIGIO, lease_break_handler)
        previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, blocked_signals)
        ledger.append(
            "SIGNAL_FENCE_ENTERED",
            blocked_signals=sorted(item.name for item in blocked_signals),
            sigio_lease_break_handler_installed=True,
        )
        for role, spec in KEY_SPECS.items():
            archives[role] = archive_key(
                role,
                spec,
                transaction_id=transaction_id,
                ledger=ledger,
            )
        postflight_processes = active_key_processes()
        postflight_descriptors = active_key_descriptors()
        if postflight_processes or postflight_descriptors:
            raise ArchiveError(
                "Archived key remains held after moves: "
                f"cmdline={postflight_processes} fds={postflight_descriptors}"
            )
        if LEASE_BREAK_EVENTS:
            raise ArchiveError(
                f"Kernel lease-break events observed: {LEASE_BREAK_EVENTS}"
            )
        validated_archives = validate_final_archives(archives, transaction_id)
        all_original_roots_absent = all(
            not os.path.lexists(spec["source"]) for spec in KEY_SPECS.values()
        )
        all_archive_roots_present = all(
            os.path.lexists(spec["destination"]) for spec in KEY_SPECS.values()
        )
        if (
            not all_original_roots_absent
            or not all_archive_roots_present
            or len(validated_archives) != 4
        ):
            raise ArchiveError("Final four-role archive presence gate failed")
        ledger.append(
            "TRANSACTION_COMMITTED",
            validated_archive_roles=sorted(validated_archives),
        )
        progress_record = ledger.finish()
        receipt_payload = {
            "schema_name": "intersubmod.v30_four_role_key_rebootstrap_archive",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "transaction_id": transaction_id,
            "source": source,
            "prepare": prepare_record,
            "progress": progress_record,
            "archives": validated_archives,
            "archive_count": len(validated_archives),
            "all_original_roots_absent": all_original_roots_absent,
            "all_archive_roots_present": all_archive_roots_present,
            "all_archive_contents_roundtrip_verified": True,
            "kernel_write_lease_gate_verified_for_all_archives": True,
            "kernel_lease_break_events": list(LEASE_BREAK_EVENTS),
            "signal_fence_covered_all_mutations_and_final_validation": True,
            "key_reuse_forbidden": True,
            "private_key_bytes_read": False,
            "release_authority_granted": False,
            "pass": True,
        }
        receipt_record = publish_json_atomic(
            RECEIPT,
            receipt_payload,
            transaction_id=transaction_id,
            stage="receipt",
        )
        committed = True
    except BaseException as error:
        if not ledger.closed:
            try:
                ledger.append(
                    "TRANSACTION_FAILED",
                    error_type=type(error).__name__,
                    error_message=str(error),
                    moved_roles=sorted(archives),
                )
            finally:
                progress_record = ledger.finish()
        if not committed and not os.path.lexists(RECEIPT) and not os.path.lexists(PARTIAL):
            role_states = {
                role: {
                    "source": str(spec["source"]),
                    "source_present": os.path.lexists(spec["source"]),
                    "destination": str(spec["destination"]),
                    "destination_present": os.path.lexists(spec["destination"]),
                }
                for role, spec in KEY_SPECS.items()
            }
            partial_payload = {
                "schema_name": "intersubmod.v30_four_role_key_rebootstrap_archive_partial",
                "schema_version": "1.0.0",
                "created_at_utc": now_utc(),
                "generation": "v30",
                "transaction_id": transaction_id,
                "source": source,
                "prepare": prepare_record,
                "progress": progress_record,
                "moved_roles": sorted(archives),
                "role_states": role_states,
                "error": {"type": type(error).__name__, "message": str(error)},
                "kernel_lease_break_events": list(LEASE_BREAK_EVENTS),
                "rebootstrap_forbidden_until_manual_reconciliation": True,
                "release_authority_granted": False,
                "pass": False,
            }
            publish_json_atomic(
                PARTIAL,
                partial_payload,
                transaction_id=transaction_id,
                stage="partial",
            )
        raise
    finally:
        if not ledger.closed:
            ledger.close_best_effort()
        if previous_mask is not None:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        if previous_sigio_handler is not None:
            signal.signal(signal.SIGIO, previous_sigio_handler)
    if not committed or receipt_record is None or progress_record is None:
        raise ArchiveError("Archive transaction returned without a committed receipt")
    print(
        json.dumps(
            {
                "prepare": prepare_record,
                "progress": progress_record,
                "receipt": receipt_record,
                "archive_count": len(archives),
                "private_key_bytes_read": False,
                "pass": True,
            },
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except BaseException as error:
        print(
            json.dumps(
                {
                    "error": {"type": type(error).__name__, "message": str(error)},
                    "pass": False,
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        exit_code = 1
    raise SystemExit(exit_code)
