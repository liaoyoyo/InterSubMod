#!/usr/bin/env python3
"""Archive the failed v8 result signer key after its pre-ready contract rejection."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import secrets
import signal
import stat
import sys
import types
from typing import Any


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
AUDIT_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "audit_notes"
)
RESULT_ROOT = AUDIT_ROOT.parent / "results"
HELPER = AUDIT_ROOT / "archive_v30_unwitnessed_four_role_keys.py"
EXPECTED_HELPER_SHA256 = (
    "793d4a34a15341b0c78d4865d103dc818b0995f2e885650f8f559dca3b9fb6ef"
)
SOURCE_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260724_all_ssnv_result_v8_v30_recovery"
)
DESTINATION_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "archive/20260724_all_ssnv_result_v8_failed_pre_ready_parent_nlink_contract_01"
)
PRIVATE_NAME = "ed25519_private_one_time.pem"
PREVIOUS_PREPARE = AUDIT_ROOT / "20260724_v30_result_v8_signer_prepare.json"
PREVIOUS_READY = AUDIT_ROOT / "20260724_v30_result_v8_signer_ready.json"
PREVIOUS_PARTIAL = AUDIT_ROOT / "20260724_v30_result_v8_signer_partial.json"
EXPECTED_PREPARE_SIZE = 1_898
EXPECTED_PREPARE_SHA256 = (
    "0b7b50f6228571fedc3f44ee8a2525d54980dbf9719d36f620aa1adc734db1ae"
)
EXPECTED_PUBLIC_SHA256 = (
    "04d1775e1dacb2cb0222816fb25244f8cdd13cf92eab86fa92d0b172898d918c"
)
TARGET = RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json"
SIGNATURE = Path(str(TARGET) + ".ed25519.sig")
PREPARE = AUDIT_ROOT / "20260724_v30_result_v8_failed_pre_ready_archive_prepare.json"
PROGRESS = AUDIT_ROOT / "20260724_v30_result_v8_failed_pre_ready_archive_progress.jsonl"
RECEIPT = AUDIT_ROOT / "20260724_v30_result_v8_failed_pre_ready_archive_receipt.json"
PARTIAL = AUDIT_ROOT / "20260724_v30_result_v8_failed_pre_ready_archive_partial.json"
ARCHIVE_RECORD_NAME = "FAILED_SIGNER_PRE_READY_RECORD.v1.json"


class V8ArchiveError(RuntimeError):
    """Raised when the failed v8 key cannot be retired unambiguously."""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def stat_identity(observed: os.stat_result) -> tuple[int, ...]:
    return (
        observed.st_dev,
        observed.st_ino,
        observed.st_mode,
        observed.st_size,
        observed.st_mtime_ns,
        observed.st_ctime_ns,
        observed.st_nlink,
    )


def load_helper():
    descriptor = os.open(HELPER, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        before = os.fstat(descriptor)
        live_before = os.lstat(HELPER)
        data = b"".join(
            os.pread(descriptor, min(1024 * 1024, before.st_size - offset), offset)
            for offset in range(0, before.st_size, 1024 * 1024)
        )
        digest = hashlib.sha256(data).hexdigest()
        if (
            not stat.S_ISREG(before.st_mode)
            or stat_identity(before) != stat_identity(live_before)
            or len(data) != before.st_size
            or stat.S_IMODE(before.st_mode) != 0o444
            or before.st_nlink != 1
            or HELPER.resolve(strict=True) != HELPER
            or digest != EXPECTED_HELPER_SHA256
        ):
            raise V8ArchiveError("Frozen archive-helper pre-execution identity drift")
        helper_record = {
            "path": str(HELPER),
            "size_bytes": len(data),
            "sha256": digest,
            "mode": "0o444",
            "link_count": 1,
        }
        module = types.ModuleType("v8_archive_helper")
        module.__file__ = str(HELPER)
        module.__package__ = None
        exec(compile(data, str(HELPER), "exec"), module.__dict__)
        after = os.fstat(descriptor)
        live_after = os.lstat(HELPER)
        if (
            stat_identity(after) != stat_identity(before)
            or stat_identity(live_after) != stat_identity(before)
        ):
            raise V8ArchiveError("Frozen archive-helper changed during bound execution")
        return module, helper_record
    finally:
        os.close(descriptor)


def publish_json(module, path: Path, payload: dict[str, Any], transaction_id: str):
    return module.publish_json_atomic(
        path,
        payload,
        transaction_id=transaction_id,
        stage=path.stem,
    )


def verify_exact_published_json(
    module,
    path: Path,
    expected_payload: dict[str, Any],
) -> dict[str, Any]:
    expected_data = json.dumps(
        expected_payload, ensure_ascii=True, indent=2, sort_keys=True
    ).encode("ascii") + b"\n"
    descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        before = os.fstat(descriptor)
        data = b"".join(
            os.pread(descriptor, min(1024 * 1024, before.st_size - offset), offset)
            for offset in range(0, before.st_size, 1024 * 1024)
        )
        after = os.fstat(descriptor)
        live = os.lstat(path)
        if (
            not stat.S_ISREG(before.st_mode)
            or stat_identity(before) != stat_identity(after)
            or stat_identity(after) != stat_identity(live)
            or stat.S_IMODE(after.st_mode) != 0o444
            or after.st_nlink != 1
            or path.resolve(strict=True) != path
            or len(data) != after.st_size
            or data != expected_data
            or module.strict_json_bytes(data, f"recovered terminal {path.name}")
            != expected_payload
        ):
            raise V8ArchiveError(f"Published terminal artifact drift: {path}")
        parent_fd = os.open(
            path.parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        )
        try:
            parent_before = os.fstat(parent_fd)
            os.fsync(parent_fd)
            parent_after = os.fstat(parent_fd)
        finally:
            os.close(parent_fd)
        after_parent_fsync = os.fstat(descriptor)
        live_after_parent_fsync = os.lstat(path)
        if (
            not stat.S_ISDIR(parent_before.st_mode)
            or (
                parent_before.st_dev,
                parent_before.st_ino,
                parent_before.st_mode,
            )
            != (
                parent_after.st_dev,
                parent_after.st_ino,
                parent_after.st_mode,
            )
            or stat_identity(after_parent_fsync) != stat_identity(before)
            or stat_identity(live_after_parent_fsync) != stat_identity(before)
        ):
            raise V8ArchiveError(
                f"Published terminal identity drift after parent fsync: {path}"
            )
        return {
            "path": str(path),
            "size_bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "mode": "0o444",
            "link_count": 1,
        }
    finally:
        os.close(descriptor)


def assert_absent_artifacts(stage: str) -> dict[str, Any]:
    paths = (PREVIOUS_READY, PREVIOUS_PARTIAL, TARGET, SIGNATURE)
    present = [str(path) for path in paths if os.path.lexists(path)]
    if present:
        raise V8ArchiveError(f"Unexpected v8 artifact at {stage}: {present}")
    return {
        "stage": stage,
        "observed_at_utc": now_utc(),
        "paths": [str(path) for path in paths],
        "all_absent": True,
    }


def external_holder_scan(module) -> dict[str, Any]:
    original_specs = module.KEY_SPECS
    module.KEY_SPECS = {
        "result_v8_failed": {
            "source": SOURCE_ROOT,
            "destination": DESTINATION_ROOT,
            "private_name": PRIVATE_NAME,
            "private_size": 290,
            "public_sha256": EXPECTED_PUBLIC_SHA256,
        }
    }
    try:
        processes = module.active_key_processes()
        descriptor_scan = module.active_key_descriptor_scan()
    finally:
        module.KEY_SPECS = original_specs
    if processes or descriptor_scan.get("matches"):
        raise V8ArchiveError(
            "Failed-v8 key remains externally held: "
            f"cmdline={processes} fds={descriptor_scan.get('matches')}"
        )
    return {
        "active_key_processes": processes,
        "active_key_descriptors": descriptor_scan.get("matches"),
        "active_key_descriptor_scan": descriptor_scan,
        "pass": True,
    }


def bound_key_inventory(
    *,
    root_path: Path,
    root_fd: int,
    public_fd: int,
    private_fd: int,
    expected_entries: set[str],
) -> dict[str, Any]:
    root_opened = os.fstat(root_fd)
    root_live = os.lstat(root_path)
    public_opened_before = os.fstat(public_fd)
    private_opened = os.fstat(private_fd)
    public_live = os.lstat(root_path / "ed25519_public.pem")
    private_live = os.lstat(root_path / PRIVATE_NAME)
    public_data = b"".join(
        os.pread(
            public_fd,
            min(1024 * 1024, public_opened_before.st_size - offset),
            offset,
        )
        for offset in range(0, public_opened_before.st_size, 1024 * 1024)
    )
    public_opened_after = os.fstat(public_fd)
    entries = set(os.listdir(root_fd))
    if (
        not stat.S_ISDIR(root_opened.st_mode)
        or stat_identity(root_opened) != stat_identity(root_live)
        or stat.S_IMODE(root_opened.st_mode) != 0o700
        or root_opened.st_uid != os.getuid()
        or root_path.resolve(strict=True) != root_path
        or entries != expected_entries
        or not stat.S_ISREG(public_opened_before.st_mode)
        or stat_identity(public_opened_before) != stat_identity(public_opened_after)
        or stat_identity(public_opened_after) != stat_identity(public_live)
        or stat.S_IMODE(public_opened_after.st_mode) != 0o444
        or public_opened_after.st_size != 113
        or public_opened_after.st_nlink != 1
        or len(public_data) != 113
        or hashlib.sha256(public_data).hexdigest() != EXPECTED_PUBLIC_SHA256
        or not stat.S_ISREG(private_opened.st_mode)
        or stat_identity(private_opened) != stat_identity(private_live)
        or stat.S_IMODE(private_opened.st_mode) != 0o400
        or private_opened.st_size != 290
        or private_opened.st_nlink != 1
    ):
        raise V8ArchiveError(f"Failed-v8 bound key inventory drift: {root_path}")
    return {
        "root": {
            "path": str(root_path),
            "device": root_opened.st_dev,
            "inode": root_opened.st_ino,
            "mode": "0o700",
            "owner_uid": root_opened.st_uid,
            "link_count": root_opened.st_nlink,
            "entries": sorted(entries),
        },
        "public_key": {
            "path": str(root_path / "ed25519_public.pem"),
            "size_bytes": len(public_data),
            "sha256": EXPECTED_PUBLIC_SHA256,
            "mode": "0o444",
            "link_count": 1,
            "device": public_opened_after.st_dev,
            "inode": public_opened_after.st_ino,
        },
        "private_key_metadata_only": {
            "path": str(root_path / PRIVATE_NAME),
            "size_bytes": private_opened.st_size,
            "mode": "0o400",
            "link_count": private_opened.st_nlink,
            "device": private_opened.st_dev,
            "inode": private_opened.st_ino,
        },
        "private_key_bytes_read": False,
        "pass": True,
    }


def preflight(module) -> dict[str, Any]:
    occupied = [
        str(path)
        for path in (DESTINATION_ROOT, PREPARE, PROGRESS, RECEIPT, PARTIAL)
        if os.path.lexists(path)
    ]
    if occupied:
        raise V8ArchiveError(f"Failed-v8 archive output occupied: {occupied}")
    absence = assert_absent_artifacts("preflight")
    prepare = module.file_record(PREVIOUS_PREPARE)
    prepare_payload = module.strict_json_bytes(
        PREVIOUS_PREPARE.read_bytes(), "v8 signer prepare"
    )
    if (
        prepare["size_bytes"] != EXPECTED_PREPARE_SIZE
        or prepare["sha256"] != EXPECTED_PREPARE_SHA256
        or prepare_payload.get("schema_name")
        != "intersubmod.attested_one_time_signer_prepare"
        or prepare_payload.get("role") != "result"
        or prepare_payload.get("key_root") != str(SOURCE_ROOT)
        or prepare_payload.get("status") != "PREPARED_NOT_READY"
        or prepare_payload.get("release_authority_granted") is not False
        or prepare_payload.get("pass") is not False
    ):
        raise V8ArchiveError("Failed-v8 prepare receipt drift")
    if (
        not SOURCE_ROOT.is_dir()
        or SOURCE_ROOT.resolve(strict=True) != SOURCE_ROOT
        or stat.S_IMODE(os.lstat(SOURCE_ROOT).st_mode) != 0o700
        or set(path.name for path in SOURCE_ROOT.iterdir())
        != {"ed25519_public.pem", PRIVATE_NAME}
    ):
        raise V8ArchiveError("Failed-v8 key-root state drift")
    public = module.file_record(SOURCE_ROOT / "ed25519_public.pem")
    private = module.private_metadata(SOURCE_ROOT / PRIVATE_NAME, 290)
    if public["size_bytes"] != 113 or public["sha256"] != EXPECTED_PUBLIC_SHA256:
        raise V8ArchiveError("Failed-v8 public-key identity drift")
    holder_scan = external_holder_scan(module)
    return {
        "previous_prepare": prepare,
        "absence": absence,
        "external_holder_scan": holder_scan,
        "public_key": public,
        "private_key_metadata_only": private,
        "private_key_bytes_read": False,
        "pass": True,
    }


def main() -> int:
    module, helper_record = load_helper()
    source_record = module.file_record(Path(__file__).resolve(strict=True))
    preflight_record = preflight(module)
    transaction_id = secrets.token_hex(16)
    prepare_payload = {
        "schema_name": "intersubmod.v8_failed_pre_ready_archive_prepare",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "transaction_id": transaction_id,
        "source": source_record,
        "frozen_helper": helper_record,
        "preflight": preflight_record,
        "source_root": str(SOURCE_ROOT),
        "destination_root": str(DESTINATION_ROOT),
        "status": "PREPARED_BEFORE_FAILED_KEY_ARCHIVE",
        "release_authority_granted": False,
        "pass": False,
    }
    prepare_record = publish_json(module, PREPARE, prepare_payload, transaction_id)
    ledger = module.ProgressLedger(PROGRESS, transaction_id)
    committed = False
    progress_record = None
    previous_mask = None
    previous_sigio_handler = None
    root_fd = public_fd = private_fd = source_parent_fd = destination_parent_fd = -1
    leases_held = False
    lease_descriptors: dict[str, tuple[int, Path]] = {}
    lease: dict[str, Any] | None = None
    source_inventory_before_record: dict[str, Any] | None = None
    source_inventory_after_record: dict[str, Any] | None = None
    destination_inventory: dict[str, Any] | None = None
    pre_move_holder_scan: dict[str, Any] | None = None
    post_move_holder_scan: dict[str, Any] | None = None
    absence_before_move: dict[str, Any] | None = None
    absence_after_move: dict[str, Any] | None = None
    absence_before_receipt: dict[str, Any] | None = None
    resources_settled = False
    runtime_cleanup: dict[str, Any] | None = None
    signal_fence_restored = False
    receipt_payload: dict[str, Any] | None = None
    receipt_publication_recovered = False

    def settle_runtime_resources() -> dict[str, Any]:
        nonlocal leases_held
        nonlocal private_fd, public_fd, root_fd
        nonlocal destination_parent_fd, source_parent_fd
        nonlocal resources_settled, runtime_cleanup
        if resources_settled:
            assert runtime_cleanup is not None
            return runtime_cleanup
        errors: list[dict[str, str]] = []
        if leases_held:
            try:
                module.release_key_write_leases(lease_descriptors)
            except BaseException as error:
                errors.append(
                    {
                        "stage": "release_write_leases",
                        "type": type(error).__name__,
                        "message": str(error),
                    }
                )
            leases_held = False
        descriptors = {
            "private_key": private_fd,
            "public_key": public_fd,
            "key_root": root_fd,
            "destination_parent": destination_parent_fd,
            "source_parent": source_parent_fd,
        }
        descriptor_states: dict[str, str] = {}
        for label, descriptor in descriptors.items():
            if descriptor < 0:
                descriptor_states[label] = "NOT_OPENED"
                continue
            try:
                os.close(descriptor)
            except BaseException as error:
                descriptor_states[label] = "CLOSE_FAILED_UNKNOWN_STATE"
                errors.append(
                    {
                        "stage": f"close_{label}",
                        "type": type(error).__name__,
                        "message": str(error),
                    }
                )
            else:
                descriptor_states[label] = "CLOSED"
        private_fd = public_fd = root_fd = -1
        destination_parent_fd = source_parent_fd = -1
        runtime_cleanup = {
            "completed_before_final_receipt_or_partial": not errors,
            "descriptors": descriptor_states,
            "write_leases_released": not any(
                item["stage"] == "release_write_leases" for item in errors
            ),
            "signal_fence_restoration_deferred_until_after_terminal_publication": True,
            "errors": errors,
        }
        resources_settled = True
        return runtime_cleanup

    def restore_signal_fence_after_terminal() -> None:
        nonlocal previous_mask, previous_sigio_handler, signal_fence_restored
        if signal_fence_restored:
            return
        if previous_sigio_handler is not None:
            signal.signal(signal.SIGIO, previous_sigio_handler)
            previous_sigio_handler = None
        if previous_mask is not None:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
            previous_mask = None
        signal_fence_restored = True

    try:
        ledger.append("TRANSACTION_STARTED", prepare=prepare_record)
        blocked = {
            signal.SIGINT,
            signal.SIGTERM,
            signal.SIGHUP,
            signal.SIGQUIT,
            signal.SIGUSR1,
            signal.SIGUSR2,
        }
        if signal.SIGIO in signal.pthread_sigmask(signal.SIG_BLOCK, set()):
            raise V8ArchiveError("SIGIO is blocked")
        module.LEASE_BREAK_EVENTS.clear()
        previous_sigio_handler = signal.getsignal(signal.SIGIO)
        signal.signal(signal.SIGIO, module.lease_break_handler)
        previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, blocked)
        ledger.append("SIGNAL_FENCE_ENTERED", sigio_handler_installed=True)

        source_parent_fd = os.open(
            SOURCE_ROOT.parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        )
        destination_parent_fd = os.open(
            DESTINATION_ROOT.parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        )
        root_fd = os.open(
            SOURCE_ROOT.name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=source_parent_fd,
        )
        public_fd = os.open(
            "ed25519_public.pem",
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=root_fd,
        )
        private_fd = os.open(
            PRIVATE_NAME,
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=root_fd,
        )
        root_opened = os.fstat(root_fd)
        root_from_parent = os.stat(
            SOURCE_ROOT.name, dir_fd=source_parent_fd, follow_symlinks=False
        )
        public_opened = os.fstat(public_fd)
        private_opened = os.fstat(private_fd)
        if (
            stat_identity(root_opened) != stat_identity(root_from_parent)
            or not stat.S_ISREG(public_opened.st_mode)
            or not stat.S_ISREG(private_opened.st_mode)
        ):
            raise V8ArchiveError("Failed-v8 bound inode identity drift")
        source_inventory_before_record = bound_key_inventory(
            root_path=SOURCE_ROOT,
            root_fd=root_fd,
            public_fd=public_fd,
            private_fd=private_fd,
            expected_entries={"ed25519_public.pem", PRIVATE_NAME},
        )
        if (
            source_inventory_before_record["public_key"]["sha256"]
            != preflight_record["public_key"]["sha256"]
            or source_inventory_before_record["private_key_metadata_only"][
                "size_bytes"
            ]
            != preflight_record["private_key_metadata_only"]["size_bytes"]
        ):
            raise V8ArchiveError("Failed-v8 preflight/bound inventory mismatch")
        failure_record = {
            "schema_name": "intersubmod.failed_signer_pre_ready",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "role_generation": "result_v8",
            "state": "FAILED_PRE_READY_PARENT_NLINK_CONTRACT_BUG",
            "observed_session_exit_code": 3,
            "observed_error": "Signer ready contract drift: result",
            "root_cause": (
                "wrapper required key-parent st_nlink stability although signer mkdir "
                "correctly increments it by exactly one"
            ),
            "previous_prepare": preflight_record["previous_prepare"],
            "ready_receipt_published": False,
            "target_or_signature_published": False,
            "private_key_bytes_read": False,
            "key_reuse_forbidden": True,
            "release_authority_granted": False,
            "pass": False,
        }
        failure_data = json.dumps(
            failure_record, ensure_ascii=True, indent=2, sort_keys=True
        ).encode("ascii") + b"\n"
        record_fd = os.open(
            ARCHIVE_RECORD_NAME,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
            0o600,
            dir_fd=root_fd,
        )
        try:
            offset = 0
            while offset < len(failure_data):
                written = os.write(record_fd, failure_data[offset:])
                if written <= 0:
                    raise V8ArchiveError("Short failed-v8 archive-record write")
                offset += written
            os.fsync(record_fd)
            os.fchmod(record_fd, 0o444)
            os.fsync(record_fd)
        finally:
            os.close(record_fd)
        os.fsync(root_fd)
        failure_artifact_before_move = module.file_record(
            SOURCE_ROOT / ARCHIVE_RECORD_NAME
        )
        if (
            failure_artifact_before_move["size_bytes"] != len(failure_data)
            or failure_artifact_before_move["sha256"]
            != hashlib.sha256(failure_data).hexdigest()
            or module.strict_json_bytes(
                (SOURCE_ROOT / ARCHIVE_RECORD_NAME).read_bytes(),
                "failed-v8 archive record before move",
            )
            != failure_record
        ):
            raise V8ArchiveError("Failed-v8 archive-record round-trip drift")
        source_inventory_after_record = bound_key_inventory(
            root_path=SOURCE_ROOT,
            root_fd=root_fd,
            public_fd=public_fd,
            private_fd=private_fd,
            expected_entries={
                "ed25519_public.pem",
                PRIVATE_NAME,
                ARCHIVE_RECORD_NAME,
            },
        )
        ledger.append(
            "FAILURE_RECORD_WRITTEN",
            artifact=failure_artifact_before_move,
            exact_inventory=source_inventory_after_record,
        )

        pre_move_holder_scan = external_holder_scan(module)
        absence_before_move = assert_absent_artifacts("immediately_before_lease")
        ledger.append(
            "PRE_MOVE_EXTERNAL_HOLDER_GATE",
            holder_scan=pre_move_holder_scan,
            absence=absence_before_move,
        )

        lease_descriptors = {
            "public_key": (public_fd, SOURCE_ROOT / "ed25519_public.pem"),
            "private_key": (private_fd, SOURCE_ROOT / PRIVATE_NAME),
        }
        lease = module.acquire_key_write_leases(
            "result_v8_failed", lease_descriptors
        )
        leases_held = True
        module.verify_key_write_leases(
            "result_v8_failed", lease_descriptors, lease, stage="before move"
        )
        absence_before_move = assert_absent_artifacts(
            "immediately_before_no_replace_move"
        )
        ledger.append("KERNEL_WRITE_LEASE_GATE_ACQUIRED", evidence=lease)
        module.move_no_replace(
            SOURCE_ROOT.name,
            DESTINATION_ROOT.name,
            source_parent_fd,
            destination_parent_fd,
        )
        os.fsync(source_parent_fd)
        os.fsync(destination_parent_fd)
        moved_from_parent = os.stat(
            DESTINATION_ROOT.name,
            dir_fd=destination_parent_fd,
            follow_symlinks=False,
        )
        if (
            moved_from_parent.st_dev,
            moved_from_parent.st_ino,
            moved_from_parent.st_mode,
        ) != (
            root_opened.st_dev,
            root_opened.st_ino,
            root_opened.st_mode,
        ) or os.path.lexists(SOURCE_ROOT):
            raise V8ArchiveError("Failed-v8 no-replace move identity drift")
        module.verify_key_write_leases(
            "result_v8_failed", lease_descriptors, lease, stage="after move"
        )
        destination_inventory = bound_key_inventory(
            root_path=DESTINATION_ROOT,
            root_fd=root_fd,
            public_fd=public_fd,
            private_fd=private_fd,
            expected_entries={
                "ed25519_public.pem",
                PRIVATE_NAME,
                ARCHIVE_RECORD_NAME,
            },
        )
        absence_after_move = assert_absent_artifacts(
            "after_move_while_write_leases_held"
        )
        if module.LEASE_BREAK_EVENTS:
            raise V8ArchiveError(
                "Kernel lease-break event observed before release: "
                f"{module.LEASE_BREAK_EVENTS}"
            )
        module.release_key_write_leases(lease_descriptors)
        leases_held = False
        if module.LEASE_BREAK_EVENTS:
            raise V8ArchiveError(
                "Kernel lease-break event observed during release: "
                f"{module.LEASE_BREAK_EVENTS}"
            )
        lease = {
            **lease,
            "verified_immediately_before_move": True,
            "verified_immediately_after_move": True,
            "released_only_after_post_move_verification": True,
        }
        ledger.append("KEY_ROOT_MOVED_NO_REPLACE", lease=lease)

        if (
            os.path.lexists(SOURCE_ROOT)
            or not DESTINATION_ROOT.is_dir()
            or destination_inventory is None
            or destination_inventory["root"]["device"] != root_opened.st_dev
            or destination_inventory["root"]["inode"] != root_opened.st_ino
            or destination_inventory["public_key"]["device"]
            != lease["files"]["public_key"]["device"]
            or destination_inventory["public_key"]["inode"]
            != lease["files"]["public_key"]["inode"]
            or destination_inventory["private_key_metadata_only"]["device"]
            != lease["files"]["private_key"]["device"]
            or destination_inventory["private_key_metadata_only"]["inode"]
            != lease["files"]["private_key"]["inode"]
        ):
            raise V8ArchiveError("Failed-v8 post-move state drift")
        public = destination_inventory["public_key"]
        private = destination_inventory["private_key_metadata_only"]
        failure_artifact = module.file_record(DESTINATION_ROOT / ARCHIVE_RECORD_NAME)
        if (
            failure_artifact != {
                **failure_artifact_before_move,
                "path": str(DESTINATION_ROOT / ARCHIVE_RECORD_NAME),
            }
            or module.strict_json_bytes(
                (DESTINATION_ROOT / ARCHIVE_RECORD_NAME).read_bytes(),
                "failed-v8 archive record after move",
            )
            != failure_record
        ):
            raise V8ArchiveError("Failed-v8 final archive-record identity drift")
        runtime_cleanup = settle_runtime_resources()
        if runtime_cleanup["completed_before_final_receipt_or_partial"] is not True:
            raise V8ArchiveError(
                f"Failed-v8 runtime cleanup failed: {runtime_cleanup['errors']}"
            )
        post_move_holder_scan = external_holder_scan(module)
        if module.LEASE_BREAK_EVENTS:
            raise V8ArchiveError(
                "Kernel lease-break event observed after postflight scan: "
                f"{module.LEASE_BREAK_EVENTS}"
            )
        absence_before_receipt = assert_absent_artifacts(
            "immediately_before_receipt_publication"
        )
        ledger.append(
            "READY_TO_PUBLISH_RECEIPT",
            final_inventory=destination_inventory,
            post_move_holder_scan=post_move_holder_scan,
            absence=absence_before_receipt,
            kernel_lease_break_events=list(module.LEASE_BREAK_EVENTS),
        )
        progress_record = ledger.finish()
        receipt_payload = {
            "schema_name": "intersubmod.v8_failed_pre_ready_archive",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "transaction_id": transaction_id,
            "source": source_record,
            "frozen_helper": helper_record,
            "prepare": prepare_record,
            "progress": progress_record,
            "previous_prepare": preflight_record["previous_prepare"],
            "archive_root": str(DESTINATION_ROOT),
            "source_inventory_before_record": source_inventory_before_record,
            "source_inventory_after_record": source_inventory_after_record,
            "destination_inventory": destination_inventory,
            "failure_record": failure_artifact,
            "public_key": public,
            "private_key_metadata_only": private,
            "kernel_write_lease_gate": lease,
            "pre_move_external_holder_scan": pre_move_holder_scan,
            "post_move_external_holder_scan": post_move_holder_scan,
            "runtime_cleanup": runtime_cleanup,
            "signal_fence_active_during_terminal_publication": (
                previous_mask is not None and previous_sigio_handler is not None
            ),
            "absence_before_move": absence_before_move,
            "absence_after_move": absence_after_move,
            "absence_before_receipt": absence_before_receipt,
            "kernel_lease_break_events": list(module.LEASE_BREAK_EVENTS),
            "source_root_absent": True,
            "ready_target_signature_absent_at_final_gate": True,
            "private_key_bytes_read": False,
            "key_reuse_forbidden": True,
            "status": "ARCHIVE_COMMITTED_BY_RECEIPT_PUBLICATION",
            "release_authority_granted": False,
            "pass": True,
        }
        try:
            receipt_record = publish_json(
                module, RECEIPT, receipt_payload, transaction_id
            )
        except BaseException:
            if not os.path.lexists(RECEIPT):
                raise
            receipt_record = verify_exact_published_json(
                module, RECEIPT, receipt_payload
            )
            receipt_publication_recovered = True
        committed = True
        print(
            json.dumps(
                {
                    "prepare": prepare_record,
                    "progress": progress_record,
                    "receipt": receipt_record,
                    "receipt_publication_recovered_after_atomic_rename": (
                        receipt_publication_recovered
                    ),
                    "archive_root": str(DESTINATION_ROOT),
                    "pass": True,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    except BaseException as error:
        runtime_cleanup = settle_runtime_resources()
        failure_handling_errors: list[dict[str, str]] = []
        if not ledger.closed:
            try:
                ledger.append(
                    "TRANSACTION_FAILED",
                    error_type=type(error).__name__,
                    error_message=str(error),
                )
            except BaseException as ledger_error:
                failure_handling_errors.append(
                    {
                        "stage": "append_transaction_failed",
                        "type": type(ledger_error).__name__,
                        "message": str(ledger_error),
                    }
                )
            try:
                progress_record = ledger.finish()
            except BaseException as ledger_error:
                failure_handling_errors.append(
                    {
                        "stage": "finish_progress_ledger",
                        "type": type(ledger_error).__name__,
                        "message": str(ledger_error),
                    }
                )
                if not ledger.closed:
                    try:
                        ledger.close_best_effort()
                    except BaseException as close_error:
                        failure_handling_errors.append(
                            {
                                "stage": "close_progress_ledger",
                                "type": type(close_error).__name__,
                                "message": str(close_error),
                            }
                        )
        if (
            not committed
            and not os.path.lexists(RECEIPT)
            and not os.path.lexists(PARTIAL)
        ):
            partial_payload = {
                "schema_name": "intersubmod.v8_failed_pre_ready_archive_partial",
                "schema_version": "1.0.0",
                "created_at_utc": now_utc(),
                "transaction_id": transaction_id,
                "source": source_record,
                "prepare": prepare_record,
                "progress": progress_record,
                "runtime_cleanup": runtime_cleanup,
                "signal_fence_active_during_terminal_publication": (
                    previous_mask is not None
                    and previous_sigio_handler is not None
                ),
                "failure_handling_errors": failure_handling_errors,
                "source_root_present": os.path.lexists(SOURCE_ROOT),
                "destination_root_present": os.path.lexists(DESTINATION_ROOT),
                "error": {"type": type(error).__name__, "message": str(error)},
                "rebootstrap_forbidden_until_manual_reconciliation": True,
                "release_authority_granted": False,
                "pass": False,
            }
            try:
                publish_json(module, PARTIAL, partial_payload, transaction_id)
            except BaseException as publication_error:
                raise V8ArchiveError(
                    "Failed-v8 partial publication failed after runtime cleanup: "
                    f"original={type(error).__name__}:{error}; "
                    f"partial={type(publication_error).__name__}:{publication_error}"
                ) from publication_error
        raise
    finally:
        if not resources_settled:
            settle_runtime_resources()
        if not ledger.closed:
            ledger.close_best_effort()
        restore_signal_fence_after_terminal()


if __name__ == "__main__":
    try:
        exit_code = main()
    except BaseException as error:
        print(
            json.dumps(
                {"error": {"type": type(error).__name__, "message": str(error)}, "pass": False},
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        exit_code = 1
    raise SystemExit(exit_code)
