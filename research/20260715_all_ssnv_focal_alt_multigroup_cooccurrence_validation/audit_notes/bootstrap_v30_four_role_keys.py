#!/usr/bin/env python3
"""Create and cross-bind all four fresh v30 release-role keypairs."""

from __future__ import annotations

from datetime import datetime, timezone
import ctypes
import hashlib
import json
import os
from pathlib import Path
import secrets
import signal
import stat
import subprocess
import sys
from typing import Any


OPENSSL = Path("/usr/bin/openssl")
EXPECTED_OPENSSL_SHA256 = (
    "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
)
AUTHORITY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/"
    "20260724_v30_four_role_rebootstrap"
)
AUTHORITY_PRIVATE = AUTHORITY_ROOT / "ed25519_private_one_time.pem"
AUTHORITY_PUBLIC = AUTHORITY_ROOT / "ed25519_public.pem"
TERMINAL_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260724_m2v5_terminal_v21"
)
TERMINAL_PRIVATE = TERMINAL_ROOT / "ed25519_private_one_time_resident.pem"
TERMINAL_PUBLIC = TERMINAL_ROOT / "ed25519_public.pem"
EXPECTED_PARENT_MODES = {
    AUTHORITY_ROOT.parent: 0o775,
    TERMINAL_ROOT.parent: 0o700,
}
RESULT_PUBLIC = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260724_all_ssnv_result_v9_v30_recovery/ed25519_public.pem"
)
RESULT_PRIVATE = RESULT_PUBLIC.with_name("ed25519_private_one_time.pem")
REPORT_PUBLIC = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
    "20260724_all_ssnv_report_v9_v30_recovery/ed25519_public.pem"
)
REPORT_PRIVATE = REPORT_PUBLIC.with_name("ed25519_private_one_time.pem")
RESULT_READY_RECEIPT = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/"
    "20260724_v30_result_v9_signer_ready.json"
)
RESULT_PREPARE_RECEIPT = RESULT_READY_RECEIPT.with_name(
    "20260724_v30_result_v9_signer_prepare.json"
)
REPORT_READY_RECEIPT = RESULT_READY_RECEIPT.with_name(
    "20260724_v30_report_v9_signer_ready.json"
)
REPORT_PREPARE_RECEIPT = RESULT_READY_RECEIPT.with_name(
    "20260724_v30_report_v9_signer_prepare.json"
)
RESULT_TARGET = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/"
    "task_b_final_dataset_release_receipt.v1.json"
)
REPORT_TARGET = RESULT_TARGET.with_name("task_b_final_report_release_receipt.v1.json")
RECEIPT = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/"
    "20260724_v30_four_role_key_bootstrap_receipt.json"
)
PREPARE_RECEIPT = RECEIPT.with_name(
    "20260724_v30_four_role_key_bootstrap_prepare.json"
)
SUCCESS_WITNESS = RECEIPT.with_name(
    "20260724_v30_four_role_key_bootstrap_success.json"
)
PROGRESS_RECEIPT = RECEIPT.with_name(
    "20260724_v30_four_role_key_bootstrap_progress.jsonl"
)
PARTIAL_RECEIPT = RECEIPT.with_name(
    "20260724_v30_four_role_key_bootstrap_partial.json"
)
ARCHIVE_PREPARE = RECEIPT.with_name(
    "20260724_v30_four_role_key_rebootstrap_archive_prepare.json"
)
ARCHIVE_RECEIPT = RECEIPT.with_name(
    "20260724_v30_four_role_key_rebootstrap_archive_receipt.json"
)
ARCHIVE_PROGRESS = RECEIPT.with_name(
    "20260724_v30_four_role_key_rebootstrap_archive_progress.jsonl"
)
ARCHIVE_PARTIAL = RECEIPT.with_name(
    "20260724_v30_four_role_key_rebootstrap_archive_partial.json"
)
FAILED_V8_ARCHIVE_PREPARE = RECEIPT.with_name(
    "20260724_v30_result_v8_failed_pre_ready_archive_prepare.json"
)
FAILED_V8_ARCHIVE_PROGRESS = RECEIPT.with_name(
    "20260724_v30_result_v8_failed_pre_ready_archive_progress.jsonl"
)
FAILED_V8_ARCHIVE_RECEIPT = RECEIPT.with_name(
    "20260724_v30_result_v8_failed_pre_ready_archive_receipt.json"
)
FAILED_V8_ARCHIVE_PARTIAL = RECEIPT.with_name(
    "20260724_v30_result_v8_failed_pre_ready_archive_partial.json"
)
FAILED_V8_PREPARE = RECEIPT.with_name(
    "20260724_v30_result_v8_signer_prepare.json"
)
FAILED_V8_SOURCE_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260724_all_ssnv_result_v8_v30_recovery"
)
FAILED_V8_ARCHIVE_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "archive/20260724_all_ssnv_result_v8_failed_pre_ready_parent_nlink_contract_01"
)
FAILED_V8_ARCHIVE_TOOL = RECEIPT.with_name("archive_failed_v8_result_pre_ready.py")
FAILED_V8_PRIVATE_NAME = "ed25519_private_one_time.pem"
FAILED_V8_RECORD_NAME = "FAILED_SIGNER_PRE_READY_RECORD.v1.json"
EXPECTED_FAILED_V8_ARCHIVE_RECEIPT_SIZE = 80_393
EXPECTED_FAILED_V8_ARCHIVE_RECEIPT_SHA256 = (
    "19aa68f7ed9d43b86d85d027377857423a03e4c6c4c8d8df744c37a91fd8e6c5"
)
EXPECTED_FAILED_V8_ARCHIVE_TOOL = {
    "path": str(FAILED_V8_ARCHIVE_TOOL),
    "size_bytes": 36_283,
    "sha256": "87a3abb23c6f88168697c8c948278151ce08d1e6d36dd1c1f0baf707f3777100",
    "mode": "0o444",
    "link_count": 1,
}
EXPECTED_FAILED_V8_HELPER = {
    "path": str(RECEIPT.with_name("archive_v30_unwitnessed_four_role_keys.py")),
    "size_bytes": 44_381,
    "sha256": "793d4a34a15341b0c78d4865d103dc818b0995f2e885650f8f559dca3b9fb6ef",
    "mode": "0o444",
    "link_count": 1,
}
EXPECTED_FAILED_V8_PREPARE = {
    "path": str(FAILED_V8_PREPARE),
    "size_bytes": 1_898,
    "sha256": "0b7b50f6228571fedc3f44ee8a2525d54980dbf9719d36f620aa1adc734db1ae",
    "mode": "0o444",
    "link_count": 1,
}
EXPECTED_FAILED_V8_PUBLIC_SHA256 = (
    "04d1775e1dacb2cb0222816fb25244f8cdd13cf92eab86fa92d0b172898d918c"
)
EXPECTED_ARCHIVE_RECEIPT_SIZE = 11_606
EXPECTED_ARCHIVE_RECEIPT_SHA256 = (
    "bf167b420eab0dbf3794bc260d292a4797263beb16a26132cd14f9c1c29ac30c"
)
ARCHIVED_KEY_ROOTS = {
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
    },
}
AT_FDCWD = -100
RENAME_NOREPLACE = 1
ENVIRONMENT = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/tmp",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
}
EXPECTED_SIGNER_OPENSSL = {
    "path": "/usr/bin/openssl",
    "size_bytes": 1_001_272,
    "sha256": EXPECTED_OPENSSL_SHA256,
    "mode": "0o755",
}
EXPECTED_SIGNING_PROTOCOL = "canonical SIGN JSON with exact path,size_bytes,sha256,mode"
EXPECTED_PUBLICATION_SEMANTICS = (
    "exit 0 publishes ceremony output but grants no release authority; "
    "an independent consumer must reopen and verify every live path"
)
EXPECTED_WRAPPER_SOURCE = {
    "path": (
        "/big7_disk/liaoyoyo2001/InterSubMod/research/"
        "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/"
        "run_attested_one_time_ed25519_signer_v7.py"
    ),
    "size_bytes": 48_702,
    "sha256": "41faf8202ac0d671cb23042c883e602484247e0ec64aef330759ad0482813058",
    "mode": "0o444",
    "link_count": 1,
}
EXPECTED_SIGNER_SOURCE = {
    "path": (
        "/big7_disk/liaoyoyo2001/InterSubMod/research/"
        "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/"
        "run_one_time_ed25519_signer_v6.py"
    ),
    "size_bytes": 34_458,
    "sha256": "c9c32222d4e8ca6daad8f6b3f8a91746fc99dd200d24cf3d77e500234c0f8bee",
    "mode": "0o444",
    "link_count": 1,
}
EXPECTED_SIGNER_PYTHON = {
    "path": "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11",
    "size_bytes": 25_409_784,
    "sha256": "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e",
    "mode": "0o775",
    "link_count": 1,
}


class BootstrapError(RuntimeError):
    """Raised before any overwrite or ambiguous key state is accepted."""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rename_no_replace(source: Path, destination: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise BootstrapError("renameat2 is required for no-replace publication")
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
        os.fsencode(destination),
        RENAME_NOREPLACE,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise BootstrapError(
            f"No-replace publication failed: {source} -> {destination}: "
            f"{os.strerror(error_number)}"
        )


def read_regular(path: Path, expected_mode: int) -> tuple[bytes, os.stat_result]:
    descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise BootstrapError(f"Not a regular file: {path}")
        chunks: list[bytes] = []
        offset = 0
        while offset < before.st_size:
            chunk = os.pread(descriptor, before.st_size - offset, offset)
            if not chunk:
                break
            chunks.append(chunk)
            offset += len(chunk)
        data = b"".join(chunks)
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
            len(data) != before.st_size
            or identity(before) != identity(after)
            or identity(after) != identity(live)
            or stat.S_IMODE(after.st_mode) != expected_mode
            or after.st_nlink != 1
        ):
            raise BootstrapError(f"File identity/mode drift: {path}")
        return data, after
    finally:
        os.close(descriptor)


def public_record(path: Path) -> dict[str, Any]:
    data, observed = read_regular(path, 0o444)
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
        "mode": oct(stat.S_IMODE(observed.st_mode)),
        "link_count": observed.st_nlink,
    }


def source_record(path: Path) -> dict[str, Any]:
    data, observed = read_regular(path, 0o444)
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
        "mode": oct(stat.S_IMODE(observed.st_mode)),
        "link_count": observed.st_nlink,
    }


def private_metadata(path: Path) -> dict[str, Any]:
    observed = os.lstat(path)
    if (
        not stat.S_ISREG(observed.st_mode)
        or path.resolve(strict=True) != path
        or stat.S_IMODE(observed.st_mode) != 0o400
        or observed.st_uid != os.getuid()
        or observed.st_nlink != 1
        or observed.st_size <= 0
    ):
        raise BootstrapError(f"Private-key metadata drift: {path}")
    return {
        "path": str(path),
        "size_bytes": observed.st_size,
        "mode": oct(stat.S_IMODE(observed.st_mode)),
        "link_count": observed.st_nlink,
    }


def write_exclusive(path: Path, data: bytes, mode: int) -> dict[str, Any]:
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
                raise BootstrapError(f"Short write: {path}")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
        observed = os.fstat(descriptor)
        if (
            observed.st_size != len(data)
            or stat.S_IMODE(observed.st_mode) != mode
            or observed.st_nlink != 1
        ):
            raise BootstrapError(f"Published key state drift: {path}")
        return {
            "path": str(path),
            "size_bytes": observed.st_size,
            "mode": oct(stat.S_IMODE(observed.st_mode)),
            "link_count": observed.st_nlink,
        }
    finally:
        os.close(descriptor)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def enforce_and_record_key_directory(path: Path) -> dict[str, Any]:
    descriptor = os.open(
        path,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
    )
    try:
        os.fchmod(descriptor, 0o700)
        os.fsync(descriptor)
        observed = os.fstat(descriptor)
        live = os.lstat(path)
        if (
            not stat.S_ISDIR(observed.st_mode)
            or (observed.st_dev, observed.st_ino, observed.st_mode, observed.st_nlink)
            != (live.st_dev, live.st_ino, live.st_mode, live.st_nlink)
            or stat.S_IMODE(observed.st_mode) != 0o700
            or observed.st_uid != os.getuid()
        ):
            raise BootstrapError(f"Key directory identity/permission drift: {path}")
        return {
            "path": str(path),
            "mode": oct(stat.S_IMODE(observed.st_mode)),
            "owner_uid": observed.st_uid,
            "link_count": observed.st_nlink,
        }
    finally:
        os.close(descriptor)


def open_bound_openssl() -> tuple[int, dict[str, Any]]:
    descriptor = os.open(OPENSSL, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    observed = os.fstat(descriptor)
    live = os.lstat(OPENSSL)
    if (
        not stat.S_ISREG(observed.st_mode)
        or (observed.st_dev, observed.st_ino, observed.st_mode, observed.st_size)
        != (live.st_dev, live.st_ino, live.st_mode, live.st_size)
        or stat.S_IMODE(observed.st_mode) != 0o755
        or observed.st_nlink != 1
    ):
        os.close(descriptor)
        raise BootstrapError("OpenSSL descriptor/path identity drift")
    data = b"".join(
        os.pread(descriptor, min(1024 * 1024, observed.st_size - offset), offset)
        for offset in range(0, observed.st_size, 1024 * 1024)
    )
    if len(data) != observed.st_size or sha256_bytes(data) != EXPECTED_OPENSSL_SHA256:
        os.close(descriptor)
        raise BootstrapError("OpenSSL content identity drift")
    return descriptor, {
        "path": str(OPENSSL),
        "size_bytes": observed.st_size,
        "sha256": sha256_bytes(data),
        "mode": oct(stat.S_IMODE(observed.st_mode)),
        "link_count": observed.st_nlink,
        "execution_path": f"/proc/self/fd/{descriptor}",
    }


def openssl(
    executable_fd: int,
    arguments: list[str],
    input_data: bytes | None = None,
) -> bytes:
    completed = subprocess.run(
        [f"/proc/self/fd/{executable_fd}", *arguments],
        input=input_data,
        capture_output=True,
        check=False,
        env=ENVIRONMENT,
        pass_fds=(executable_fd,),
        timeout=60,
    )
    if completed.returncode != 0 or completed.stderr:
        raise BootstrapError(
            f"OpenSSL failed: command={arguments} exit={completed.returncode} "
            f"stderr_sha256={sha256_bytes(completed.stderr)}"
        )
    return completed.stdout


def create_staged_keypair(
    executable_fd: int,
    root: Path,
    private_filename: str,
) -> dict[str, Any]:
    if os.path.lexists(root):
        raise BootstrapError(f"Staging key root already exists: {root}")
    parent = root.parent.resolve(strict=True)
    parent_stat = os.lstat(parent)
    expected_parent_mode = EXPECTED_PARENT_MODES.get(parent)
    if (
        parent != root.parent
        or expected_parent_mode is None
        or stat.S_IMODE(parent_stat.st_mode) != expected_parent_mode
        or parent_stat.st_uid != os.getuid()
        or stat.S_IMODE(parent_stat.st_mode) & 0o002
    ):
        raise BootstrapError(f"Key parent identity/permission drift: {parent}")
    os.mkdir(root, mode=0o700)
    root_record = enforce_and_record_key_directory(root)
    private_path = root / private_filename
    public_path = root / "ed25519_public.pem"
    private_pem = openssl(
        executable_fd, ["genpkey", "-algorithm", "Ed25519"]
    )
    public_pem = openssl(
        executable_fd, ["pkey", "-pubout"], input_data=private_pem
    )
    if (
        not private_pem.startswith(b"-----BEGIN PRIVATE KEY-----\n")
        or not public_pem.startswith(b"-----BEGIN PUBLIC KEY-----\n")
        or len(private_pem) != 119
        or len(public_pem) != 113
    ):
        raise BootstrapError("Unexpected Ed25519 PEM encoding")
    private_written = write_exclusive(private_path, private_pem, 0o400)
    public_written = write_exclusive(public_path, public_pem, 0o444)
    fsync_directory(root)
    return {
        "root": str(root),
        "root_metadata": root_record,
        "parent": {
            "path": str(parent),
            "mode": oct(stat.S_IMODE(parent_stat.st_mode)),
            "owner_uid": parent_stat.st_uid,
            "world_writable": False,
        },
        "private_key_metadata_only": private_written,
        "private_key_bytes_recorded_in_receipt": False,
        "public_key": {
            **public_written,
            "sha256": sha256_bytes(public_pem),
        },
    }


def strict_json_bytes(data: bytes, label: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise BootstrapError(f"Duplicate JSON key in {label}: {key}")
            output[key] = value
        return output

    try:
        value = json.loads(
            data,
            object_pairs_hook=no_duplicates,
            parse_constant=lambda item: (_ for _ in ()).throw(
                BootstrapError(f"Non-finite JSON in {label}: {item}")
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise BootstrapError(f"Invalid JSON in {label}") from error
    if not isinstance(value, dict):
        raise BootstrapError(f"JSON root is not an object: {label}")
    return value


def json_record(path: Path) -> dict[str, Any]:
    data, observed = read_regular(path, 0o444)
    value = strict_json_bytes(data, str(path))
    return {
        "artifact": {
            "path": str(path),
            "size_bytes": len(data),
            "sha256": sha256_bytes(data),
            "mode": oct(stat.S_IMODE(observed.st_mode)),
            "link_count": observed.st_nlink,
        },
        "payload": value,
    }


def exact_json_record(
    path: Path,
    expected_payload: dict[str, Any],
    *,
    recover_parent_durability: bool = False,
) -> dict[str, Any]:
    expected_data = json.dumps(
        expected_payload, ensure_ascii=True, indent=2, sort_keys=True
    ).encode("ascii") + b"\n"
    data, observed = read_regular(path, 0o444)
    if (
        data != expected_data
        or strict_json_bytes(data, f"exact JSON {path}") != expected_payload
    ):
        raise BootstrapError(f"Canonical JSON bytes drift: {path}")
    if recover_parent_durability:
        fsync_directory(path.parent)
        recovered_data, recovered = read_regular(path, 0o444)
        if (
            recovered_data != expected_data
            or (
                observed.st_dev,
                observed.st_ino,
                observed.st_mode,
                observed.st_size,
                observed.st_mtime_ns,
                observed.st_ctime_ns,
                observed.st_nlink,
            )
            != (
                recovered.st_dev,
                recovered.st_ino,
                recovered.st_mode,
                recovered.st_size,
                recovered.st_mtime_ns,
                recovered.st_ctime_ns,
                recovered.st_nlink,
            )
        ):
            raise BootstrapError(f"Recovered canonical JSON identity drift: {path}")
    return {
        "artifact": {
            "path": str(path),
            "size_bytes": len(data),
            "sha256": sha256_bytes(data),
            "mode": "0o444",
            "link_count": 1,
        },
        "payload": expected_payload,
    }


def publish_json_atomic(
    path: Path,
    payload: dict[str, Any],
    *,
    transaction_id: str,
    stage_role: str,
) -> dict[str, Any]:
    staging = path.with_name(f".{path.name}.staging.{transaction_id}.{stage_role}")
    if os.path.lexists(path) or os.path.lexists(staging):
        raise BootstrapError(f"Bootstrap publication path occupied: {path}")
    data = json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True).encode(
        "ascii"
    ) + b"\n"
    write_exclusive(staging, data, 0o444)
    staged = exact_json_record(staging, payload)
    rename_no_replace(staging, path)
    fsync_directory(path.parent)
    published = exact_json_record(path, payload)
    if (
        published["payload"] != payload
        or published["artifact"]["sha256"] != staged["artifact"]["sha256"]
        or published["artifact"]["size_bytes"] != staged["artifact"]["size_bytes"]
    ):
        raise BootstrapError(f"Published bootstrap JSON drift: {stage_role}")
    return published["artifact"]


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
            raise BootstrapError("Bootstrap progress ledger is closed")
        payload = {
            "created_at_utc": now_utc(),
            "transaction_id": self.transaction_id,
            "event": event,
            **details,
        }
        data = json.dumps(payload, ensure_ascii=True, sort_keys=True).encode("ascii") + b"\n"
        offset = 0
        while offset < len(data):
            written = os.write(self.fd, data[offset:])
            if written <= 0:
                raise BootstrapError("Short bootstrap progress-ledger write")
            offset += written
        os.fsync(self.fd)

    def finish(self) -> dict[str, Any]:
        if not self.closed:
            os.fchmod(self.fd, 0o444)
            os.fsync(self.fd)
            os.close(self.fd)
            self.closed = True
            fsync_directory(self.path.parent)
        return source_record(self.path)

    def close_best_effort(self) -> None:
        if self.closed:
            return
        try:
            os.fchmod(self.fd, 0o444)
            os.fsync(self.fd)
        finally:
            os.close(self.fd)
            self.closed = True


def validate_pre_authority_archive() -> dict[str, Any]:
    if (
        EXPECTED_ARCHIVE_RECEIPT_SIZE <= 0
        or len(EXPECTED_ARCHIVE_RECEIPT_SHA256) != 64
    ):
        raise BootstrapError("Archive receipt identity has not been pinned")
    if os.path.lexists(ARCHIVE_PARTIAL):
        raise BootstrapError("Partial archive witness forbids four-role bootstrap")
    prepare = json_record(ARCHIVE_PREPARE)
    progress = source_record(ARCHIVE_PROGRESS)
    receipt = json_record(ARCHIVE_RECEIPT)
    payload = receipt["payload"]
    if (
        receipt["artifact"]["size_bytes"] != EXPECTED_ARCHIVE_RECEIPT_SIZE
        or receipt["artifact"]["sha256"] != EXPECTED_ARCHIVE_RECEIPT_SHA256
        or payload.get("schema_name")
        != "intersubmod.v30_four_role_key_rebootstrap_archive"
        or payload.get("schema_version") != "1.0.0"
        or payload.get("generation") != "v30"
        or payload.get("prepare") != prepare["artifact"]
        or payload.get("progress") != progress
        or payload.get("archive_count") != 4
        or payload.get("all_original_roots_absent") is not True
        or payload.get("all_archive_roots_present") is not True
        or payload.get("all_archive_contents_roundtrip_verified") is not True
        or payload.get("kernel_write_lease_gate_verified_for_all_archives")
        is not True
        or payload.get("kernel_lease_break_events") != []
        or payload.get("signal_fence_covered_all_mutations_and_final_validation")
        is not True
        or payload.get("key_reuse_forbidden") is not True
        or payload.get("private_key_bytes_read") is not False
        or payload.get("release_authority_granted") is not False
        or payload.get("pass") is not True
        or set(payload.get("archives", {})) != set(ARCHIVED_KEY_ROOTS)
    ):
        raise BootstrapError("Pre-authority four-role archive receipt drift")
    validated: dict[str, Any] = {}
    for role, spec in ARCHIVED_KEY_ROOTS.items():
        source = Path(spec["source"])
        destination = Path(spec["destination"])
        private_name = str(spec["private_name"])
        observed = payload["archives"][role]
        if (
            os.path.lexists(source)
            or not destination.is_dir()
            or destination.resolve(strict=True) != destination
            or stat.S_IMODE(os.lstat(destination).st_mode) != 0o700
            or set(path.name for path in destination.iterdir())
            != {
                "ed25519_public.pem",
                private_name,
                "FAILED_KEY_ARCHIVE_RECORD.v1.json",
            }
        ):
            raise BootstrapError(f"Archived key root state drift: {role}")
        archive_record = source_record(
            destination / "FAILED_KEY_ARCHIVE_RECORD.v1.json"
        )
        public = public_record(destination / "ed25519_public.pem")
        private = private_metadata(destination / private_name)
        lease = observed.get("kernel_write_lease_gate", {})
        lease_files = lease.get("files", {}) if isinstance(lease, dict) else {}
        public_stat = os.lstat(destination / "ed25519_public.pem")
        private_stat = os.lstat(destination / private_name)
        if (
            observed.get("archive_root") != str(destination)
            or observed.get("archive_record") != archive_record
            or observed.get("public_key") != public
            or observed.get("private_key_metadata_only") != private
            or set(lease_files) != {"public_key", "private_key"}
            or lease.get("existing_external_openers_rejected_by_kernel") is not True
            or lease.get("verified_immediately_before_move") is not True
            or lease.get("verified_immediately_after_move") is not True
            or lease.get("released_only_after_post_move_verification") is not True
            or lease_files["public_key"].get("device") != public_stat.st_dev
            or lease_files["public_key"].get("inode") != public_stat.st_ino
            or lease_files["private_key"].get("device") != private_stat.st_dev
            or lease_files["private_key"].get("inode") != private_stat.st_ino
            or observed.get("key_reuse_forbidden") is not True
            or observed.get("pass") is not True
        ):
            raise BootstrapError(f"Archived key evidence drift: {role}")
        validated[role] = observed
    return {
        "prepare": prepare["artifact"],
        "progress": progress,
        "receipt": receipt["artifact"],
        "archives": validated,
        "partial_witness_absent": True,
        "all_original_roots_absent": True,
        "key_reuse_forbidden": True,
        "private_key_bytes_read": False,
        "pass": True,
    }


def validate_failed_v8_archive() -> dict[str, Any]:
    if (
        EXPECTED_FAILED_V8_ARCHIVE_RECEIPT_SIZE <= 0
        or len(EXPECTED_FAILED_V8_ARCHIVE_RECEIPT_SHA256) != 64
    ):
        raise BootstrapError("Failed-v8 archive receipt identity is not pinned")
    failed_v8_ready = FAILED_V8_PREPARE.with_name(
        "20260724_v30_result_v8_signer_ready.json"
    )
    failed_v8_partial = FAILED_V8_PREPARE.with_name(
        "20260724_v30_result_v8_signer_partial.json"
    )
    protected_absence = (
        FAILED_V8_SOURCE_ROOT,
        FAILED_V8_ARCHIVE_PARTIAL,
        failed_v8_ready,
        failed_v8_partial,
        RESULT_TARGET,
        Path(str(RESULT_TARGET) + ".ed25519.sig"),
    )
    occupied = [str(path) for path in protected_absence if os.path.lexists(path)]
    if occupied:
        raise BootstrapError(f"Failed-v8 archive terminal state drift: {occupied}")

    prepare = json_record(FAILED_V8_ARCHIVE_PREPARE)
    progress = source_record(FAILED_V8_ARCHIVE_PROGRESS)
    receipt = json_record(FAILED_V8_ARCHIVE_RECEIPT)
    previous_prepare = source_record(FAILED_V8_PREPARE)
    archive_tool = source_record(FAILED_V8_ARCHIVE_TOOL)
    helper = source_record(Path(EXPECTED_FAILED_V8_HELPER["path"]))
    payload = receipt["payload"]
    expected_keys = {
        "absence_after_move",
        "absence_before_move",
        "absence_before_receipt",
        "archive_root",
        "created_at_utc",
        "destination_inventory",
        "failure_record",
        "frozen_helper",
        "kernel_lease_break_events",
        "kernel_write_lease_gate",
        "key_reuse_forbidden",
        "pass",
        "post_move_external_holder_scan",
        "pre_move_external_holder_scan",
        "prepare",
        "previous_prepare",
        "private_key_bytes_read",
        "private_key_metadata_only",
        "progress",
        "public_key",
        "ready_target_signature_absent_at_final_gate",
        "release_authority_granted",
        "runtime_cleanup",
        "schema_name",
        "schema_version",
        "signal_fence_active_during_terminal_publication",
        "source",
        "source_inventory_after_record",
        "source_inventory_before_record",
        "source_root_absent",
        "status",
        "transaction_id",
    }
    if (
        receipt["artifact"]["size_bytes"]
        != EXPECTED_FAILED_V8_ARCHIVE_RECEIPT_SIZE
        or receipt["artifact"]["sha256"]
        != EXPECTED_FAILED_V8_ARCHIVE_RECEIPT_SHA256
        or set(payload) != expected_keys
        or payload.get("schema_name")
        != "intersubmod.v8_failed_pre_ready_archive"
        or payload.get("schema_version") != "1.0.0"
        or payload.get("status") != "ARCHIVE_COMMITTED_BY_RECEIPT_PUBLICATION"
        or payload.get("prepare") != prepare["artifact"]
        or payload.get("progress") != progress
        or payload.get("previous_prepare") != previous_prepare
        or payload.get("source") != archive_tool
        or payload.get("frozen_helper") != helper
        or archive_tool != EXPECTED_FAILED_V8_ARCHIVE_TOOL
        or helper != EXPECTED_FAILED_V8_HELPER
        or previous_prepare != EXPECTED_FAILED_V8_PREPARE
        or payload.get("archive_root") != str(FAILED_V8_ARCHIVE_ROOT)
        or payload.get("source_root_absent") is not True
        or payload.get("ready_target_signature_absent_at_final_gate") is not True
        or payload.get("private_key_bytes_read") is not False
        or payload.get("key_reuse_forbidden") is not True
        or payload.get("signal_fence_active_during_terminal_publication")
        is not True
        or payload.get("kernel_lease_break_events") != []
        or payload.get("release_authority_granted") is not False
        or payload.get("pass") is not True
    ):
        raise BootstrapError("Failed-v8 archive receipt contract drift")

    expected_entries = {
        "ed25519_public.pem",
        FAILED_V8_PRIVATE_NAME,
        FAILED_V8_RECORD_NAME,
    }
    if (
        not FAILED_V8_ARCHIVE_ROOT.is_dir()
        or FAILED_V8_ARCHIVE_ROOT.resolve(strict=True) != FAILED_V8_ARCHIVE_ROOT
        or stat.S_IMODE(os.lstat(FAILED_V8_ARCHIVE_ROOT).st_mode) != 0o700
        or set(path.name for path in FAILED_V8_ARCHIVE_ROOT.iterdir())
        != expected_entries
    ):
        raise BootstrapError("Failed-v8 archive root inventory drift")
    root_stat = os.lstat(FAILED_V8_ARCHIVE_ROOT)
    public_path = FAILED_V8_ARCHIVE_ROOT / "ed25519_public.pem"
    private_path = FAILED_V8_ARCHIVE_ROOT / FAILED_V8_PRIVATE_NAME
    failure_path = FAILED_V8_ARCHIVE_ROOT / FAILED_V8_RECORD_NAME
    public = public_record(public_path)
    private = private_metadata(private_path)
    failure = json_record(failure_path)
    public_stat = os.lstat(public_path)
    private_stat = os.lstat(private_path)
    public_enriched = {
        **public,
        "device": public_stat.st_dev,
        "inode": public_stat.st_ino,
    }
    private_enriched = {
        **private,
        "device": private_stat.st_dev,
        "inode": private_stat.st_ino,
    }
    destination_inventory = payload.get("destination_inventory", {})
    destination_root = destination_inventory.get("root", {})
    if (
        public["size_bytes"] != 113
        or public["sha256"] != EXPECTED_FAILED_V8_PUBLIC_SHA256
        or private["size_bytes"] != 290
        or payload.get("public_key") != public_enriched
        or payload.get("private_key_metadata_only") != private_enriched
        or payload.get("failure_record") != failure["artifact"]
        or failure["payload"].get("schema_name")
        != "intersubmod.failed_signer_pre_ready"
        or failure["payload"].get("state")
        != "FAILED_PRE_READY_PARENT_NLINK_CONTRACT_BUG"
        or failure["payload"].get("private_key_bytes_read") is not False
        or failure["payload"].get("key_reuse_forbidden") is not True
        or failure["payload"].get("release_authority_granted") is not False
        or failure["payload"].get("pass") is not False
        or destination_inventory.get("public_key") != public_enriched
        or destination_inventory.get("private_key_metadata_only")
        != private_enriched
        or destination_inventory.get("private_key_bytes_read") is not False
        or destination_inventory.get("pass") is not True
        or destination_root.get("path") != str(FAILED_V8_ARCHIVE_ROOT)
        or destination_root.get("device") != root_stat.st_dev
        or destination_root.get("inode") != root_stat.st_ino
        or destination_root.get("mode") != "0o700"
        or destination_root.get("owner_uid") != os.getuid()
        or destination_root.get("entries") != sorted(expected_entries)
    ):
        raise BootstrapError("Failed-v8 final archive evidence drift")

    lease = payload.get("kernel_write_lease_gate", {})
    lease_files = lease.get("files", {}) if isinstance(lease, dict) else {}
    expected_lease_keys = {
        "event_offset",
        "existing_external_openers_rejected_by_kernel",
        "files",
        "new_external_opens_blocked_until_post_move_verification",
        "pass",
        "released_only_after_post_move_verification",
        "role",
        "verified_immediately_after_move",
        "verified_immediately_before_move",
    }
    if (
        set(lease) != expected_lease_keys
        or set(lease_files) != {"public_key", "private_key"}
        or lease.get("role") != "result_v8_failed"
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
        raise BootstrapError("Failed-v8 kernel write-lease evidence drift")

    for label in ("pre_move_external_holder_scan", "post_move_external_holder_scan"):
        scan = payload.get(label, {})
        descriptor_scan = scan.get("active_key_descriptor_scan", {})
        if (
            scan.get("active_key_processes") != []
            or scan.get("active_key_descriptors") != []
            or scan.get("pass") is not True
            or descriptor_scan.get("matches") != []
            or descriptor_scan.get("pass") is not True
        ):
            raise BootstrapError(f"Failed-v8 external-holder evidence drift: {label}")
    absence_stages = {
        "absence_before_move": "immediately_before_no_replace_move",
        "absence_after_move": "after_move_while_write_leases_held",
        "absence_before_receipt": "immediately_before_receipt_publication",
    }
    for label, expected_stage in absence_stages.items():
        absence = payload.get(label, {})
        if absence.get("stage") != expected_stage or absence.get("all_absent") is not True:
            raise BootstrapError(f"Failed-v8 protected-absence evidence drift: {label}")
    cleanup = payload.get("runtime_cleanup", {})
    if (
        cleanup.get("completed_before_final_receipt_or_partial") is not True
        or cleanup.get("write_leases_released") is not True
        or cleanup.get(
            "signal_fence_restoration_deferred_until_after_terminal_publication"
        )
        is not True
        or cleanup.get("errors") != []
        or set(cleanup.get("descriptors", {}).values()) != {"CLOSED"}
    ):
        raise BootstrapError("Failed-v8 runtime-cleanup evidence drift")

    progress_data, _ = read_regular(FAILED_V8_ARCHIVE_PROGRESS, 0o444)
    progress_rows = [json.loads(line) for line in progress_data.splitlines()]
    if (
        not progress_rows
        or any(row.get("transaction_id") != payload.get("transaction_id") for row in progress_rows)
        or progress_rows[-1].get("event") != "READY_TO_PUBLISH_RECEIPT"
        or any(row.get("event") == "TRANSACTION_COMMITTED" for row in progress_rows)
    ):
        raise BootstrapError("Failed-v8 progress ledger drift")
    return {
        "prepare": prepare["artifact"],
        "progress": progress,
        "receipt": receipt["artifact"],
        "previous_prepare": previous_prepare,
        "archive_root": str(FAILED_V8_ARCHIVE_ROOT),
        "failure_record": failure["artifact"],
        "public_key": public_enriched,
        "private_key_metadata_only": private_enriched,
        "kernel_write_lease_gate": lease,
        "partial_witness_absent": True,
        "source_root_absent": True,
        "private_key_bytes_read": False,
        "key_reuse_forbidden": True,
        "pass": True,
    }


def validate_attested_signer_ready(
    *,
    role: str,
    prepare_path: Path,
    ready_path: Path,
    public_path: Path,
    private_path: Path,
    expected_target: Path,
) -> dict[str, Any]:
    prepare = json_record(prepare_path)
    ready = json_record(ready_path)
    prepare_payload = prepare["payload"]
    ready_payload = ready["payload"]
    live_public = public_record(public_path)
    live_private = private_metadata(private_path)
    creation = ready_payload.get("creation_contract")
    key_material = ready_payload.get("key_material")
    absence = (
        creation.get("key_root_absence_observation", {})
        if isinstance(creation, dict)
        else {}
    )
    parent_transition = (
        creation.get("key_parent_identity_transition_verified", {})
        if isinstance(creation, dict)
        else {}
    )
    prepare_parent = prepare_payload.get("key_parent", {})
    root_inventory = (
        creation.get("root_exact_inventory_verified", {})
        if isinstance(creation, dict)
        else {}
    )
    child_ready = (
        key_material.get("ready_payload", {})
        if isinstance(key_material, dict)
        else {}
    )
    child_bootstrap = (
        child_ready.get("bootstrap", {}) if isinstance(child_ready, dict) else {}
    )
    signer_source = ready_payload.get("signer_source")
    wrapper_source = ready_payload.get("wrapper_source")
    python_record = ready_payload.get("python")
    signer_artifact = (
        {
            key: signer_source.get(key)
            for key in ("path", "size_bytes", "sha256", "mode")
        }
        if isinstance(signer_source, dict)
        else None
    )
    python_artifact = (
        {
            key: python_record.get(key)
            for key in ("path", "size_bytes", "sha256", "mode")
        }
        if isinstance(python_record, dict)
        else None
    )
    if (
        prepare_payload.get("schema_name")
        != "intersubmod.attested_one_time_signer_prepare"
        or prepare_payload.get("schema_version") != "1.0.0"
        or prepare_payload.get("generation") != "v30"
        or prepare_payload.get("role") != role
        or prepare_payload.get("key_root") != str(public_path.parent)
        or prepare_payload.get("key_root_absent_before_child_launch") is not True
        or prepare_payload.get("expected_target") != str(expected_target)
        or prepare_payload.get("expected_target_and_signature_absent") is not True
        or prepare_payload.get("status") != "PREPARED_NOT_READY"
        or prepare_payload.get("release_authority_granted") is not False
        or prepare_payload.get("pass") is not False
        or ready_payload.get("schema_name")
        != "intersubmod.attested_one_time_signer_ready"
        or ready_payload.get("schema_version") != "1.0.0"
        or ready_payload.get("generation") != "v30"
        or ready_payload.get("role") != role
        or ready_payload.get("status") != "READY_AWAITING_EXACT_SIGN_TARGET"
        or ready_payload.get("prepare") != prepare["artifact"]
        or ready_payload.get("transaction_id") != prepare_payload.get("transaction_id")
        or prepare_payload.get("wrapper_source") != wrapper_source
        or prepare_payload.get("signer_source") != signer_source
        or prepare_payload.get("python") != python_record
        or ready_payload.get("expected_target") != str(expected_target)
        or ready_payload.get("expected_target_and_signature_absent_at_ready")
        is not True
        or ready_payload.get("release_authority_granted") is not False
        or ready_payload.get("pass") is not True
        or not isinstance(creation, dict)
        or set(creation)
        != {
            "key_root_absence_observation",
            "key_parent_identity_transition_verified",
            "signer_exclusive_creation_source_bound",
            "root_exact_inventory_verified",
        }
        or set(absence)
        != {
            "parent_device",
            "parent_inode",
            "root_name",
            "present",
            "observed_before_child_launch",
        }
        or not isinstance(absence.get("parent_device"), int)
        or not isinstance(absence.get("parent_inode"), int)
        or absence.get("root_name") != public_path.parent.name
        or absence.get("present") is not False
        or absence.get("observed_before_child_launch") is not True
        or set(parent_transition)
        != {
            "path",
            "device",
            "inode",
            "mode",
            "link_count_before",
            "link_count_after",
            "expected_link_count_delta_for_one_new_subdirectory",
        }
        or parent_transition.get("path") != str(public_path.parent.parent)
        or parent_transition.get("device") != absence.get("parent_device")
        or parent_transition.get("inode") != absence.get("parent_inode")
        or parent_transition.get("mode") != prepare_parent.get("mode")
        or parent_transition.get("link_count_before")
        != prepare_parent.get("link_count")
        or not isinstance(parent_transition.get("link_count_before"), int)
        or parent_transition.get("link_count_after")
        != parent_transition.get("link_count_before") + 1
        or parent_transition.get(
            "expected_link_count_delta_for_one_new_subdirectory"
        )
        != 1
        or creation.get("signer_exclusive_creation_source_bound") != signer_source
        or not isinstance(key_material, dict)
        or key_material.get("private_key_bytes_read") is not False
        or key_material.get("root") != str(public_path.parent)
        or root_inventory != key_material.get("root_metadata")
        or key_material.get("public_key") != live_public
        or key_material.get("private_key_metadata_only") != live_private
        or key_material.get("pass") is not True
        or child_ready.get("event") != "SIGNER_READY"
        or child_ready.get("openssl") != EXPECTED_SIGNER_OPENSSL
        or child_ready.get("expected_target") != str(expected_target)
        or child_ready.get("signing_protocol") != EXPECTED_SIGNING_PROTOCOL
        or child_ready.get("publication_semantics")
        != EXPECTED_PUBLICATION_SEMANTICS
        or child_bootstrap.get("script") != signer_artifact
        or child_bootstrap.get("python") != python_artifact
        or child_bootstrap.get("isolated_mode") is not True
        or child_bootstrap.get("no_user_site") is not True
        or wrapper_source != EXPECTED_WRAPPER_SOURCE
        or signer_source != EXPECTED_SIGNER_SOURCE
        or python_record != EXPECTED_SIGNER_PYTHON
        or os.path.lexists(expected_target)
        or os.path.lexists(Path(str(expected_target) + ".ed25519.sig"))
    ):
        raise BootstrapError(f"Attested {role} signer-ready contract drift")
    return {
        "prepare_receipt": prepare["artifact"],
        "ready_receipt": ready["artifact"],
        "key_root": key_material["root"],
        "root_metadata": key_material["root_metadata"],
        "private_key_metadata_only": live_private,
        "private_key_bytes_read": False,
        "public_key": live_public,
        "wrapper_source": wrapper_source,
        "signer_source": signer_source,
        "python": python_record,
        "creation_contract": dict(creation),
        "pass": True,
    }


def main() -> int:
    occupied = [
        str(path)
        for path in (
            AUTHORITY_ROOT,
            TERMINAL_ROOT,
            PREPARE_RECEIPT,
            RECEIPT,
            SUCCESS_WITNESS,
            PROGRESS_RECEIPT,
            PARTIAL_RECEIPT,
        )
        if os.path.lexists(path)
    ]
    if occupied:
        raise BootstrapError(f"Bootstrap output already exists: {occupied}")
    archive_attestation = validate_pre_authority_archive()
    failed_v8_archive_attestation = validate_failed_v8_archive()
    source = source_record(Path(__file__).resolve(strict=True))
    executable_fd, openssl_record = open_bound_openssl()
    result_attestation = validate_attested_signer_ready(
        role="result",
        prepare_path=RESULT_PREPARE_RECEIPT,
        ready_path=RESULT_READY_RECEIPT,
        public_path=RESULT_PUBLIC,
        private_path=RESULT_PRIVATE,
        expected_target=RESULT_TARGET,
    )
    report_attestation = validate_attested_signer_ready(
        role="report",
        prepare_path=REPORT_PREPARE_RECEIPT,
        ready_path=REPORT_READY_RECEIPT,
        public_path=REPORT_PUBLIC,
        private_path=REPORT_PRIVATE,
        expected_target=REPORT_TARGET,
    )
    result_public = result_attestation["public_key"]
    report_public = report_attestation["public_key"]
    transaction_id = secrets.token_hex(16)
    authority_staging = AUTHORITY_ROOT.with_name(
        f".{AUTHORITY_ROOT.name}.staging.{transaction_id}"
    )
    terminal_staging = TERMINAL_ROOT.with_name(
        f".{TERMINAL_ROOT.name}.staging.{transaction_id}"
    )
    receipt_staging = RECEIPT.with_name(
        f".{RECEIPT.name}.staging.{transaction_id}"
    )
    witness_staging = SUCCESS_WITNESS.with_name(
        f".{SUCCESS_WITNESS.name}.staging.{transaction_id}"
    )
    prepare = {
        "schema_name": "intersubmod.v30_four_role_key_bootstrap_prepare",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "generation": "v30",
        "transaction_id": transaction_id,
        "source": source,
        "openssl": openssl_record,
        "intended_outputs": {
            "authority_root": str(AUTHORITY_ROOT),
            "terminal_root": str(TERMINAL_ROOT),
            "receipt": str(RECEIPT),
            "success_witness": str(SUCCESS_WITNESS),
            "progress": str(PROGRESS_RECEIPT),
            "partial_on_failure": str(PARTIAL_RECEIPT),
        },
        "staging_paths": {
            "authority_root": str(authority_staging),
            "terminal_root": str(terminal_staging),
            "receipt": str(receipt_staging),
            "success_witness": str(witness_staging),
        },
        "attested_signer_ready_inputs": {
            "result": {
                "root": str(RESULT_PUBLIC.parent),
                "ready_receipt": result_attestation["ready_receipt"],
            },
            "report": {
                "root": str(REPORT_PUBLIC.parent),
                "ready_receipt": report_attestation["ready_receipt"],
            },
        },
        "pre_authority_archive": archive_attestation,
        "failed_v8_pre_ready_archive": failed_v8_archive_attestation,
        "status": "PREPARED_NOT_COMMITTED",
        "release_authority_granted": False,
        "pass": False,
    }
    write_exclusive(
        PREPARE_RECEIPT,
        json.dumps(prepare, ensure_ascii=True, indent=2, sort_keys=True).encode(
            "ascii"
        )
        + b"\n",
        0o444,
    )
    fsync_directory(PREPARE_RECEIPT.parent)
    prepare_record = exact_json_record(PREPARE_RECEIPT, prepare)

    ledger: ProgressLedger | None = None
    progress_record: dict[str, Any] | None = None
    committed = False
    receipt_publication_recovered = False
    witness_publication_recovered = False
    previous_mask: set[signal.Signals] | None = None
    blocked_signals = {
        signal.SIGINT,
        signal.SIGTERM,
        signal.SIGHUP,
        signal.SIGQUIT,
        signal.SIGUSR1,
        signal.SIGUSR2,
    }
    try:
        previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, blocked_signals)
        ledger = ProgressLedger(PROGRESS_RECEIPT, transaction_id)
        ledger.append("TRANSACTION_STARTED", prepare=prepare_record["artifact"])
        ledger.append(
            "SIGNAL_FENCE_ENTERED",
            blocked_signals=sorted(item.name for item in blocked_signals),
        )
        create_staged_keypair(
            executable_fd,
            authority_staging,
            AUTHORITY_PRIVATE.name,
        )
        ledger.append("AUTHORITY_STAGING_KEY_CREATED", path=str(authority_staging))
        create_staged_keypair(
            executable_fd,
            terminal_staging,
            TERMINAL_PRIVATE.name,
        )
        ledger.append("TERMINAL_STAGING_KEY_CREATED", path=str(terminal_staging))
        rename_no_replace(authority_staging, AUTHORITY_ROOT)
        fsync_directory(AUTHORITY_ROOT.parent)
        ledger.append("AUTHORITY_KEY_ROOT_PUBLISHED", path=str(AUTHORITY_ROOT))
        rename_no_replace(terminal_staging, TERMINAL_ROOT)
        fsync_directory(TERMINAL_ROOT.parent)
        ledger.append("TERMINAL_KEY_ROOT_PUBLISHED", path=str(TERMINAL_ROOT))

        authority_root_record = enforce_and_record_key_directory(AUTHORITY_ROOT)
        terminal_root_record = enforce_and_record_key_directory(TERMINAL_ROOT)
        authority_public = public_record(AUTHORITY_PUBLIC)
        terminal_public = public_record(TERMINAL_PUBLIC)
        authority_private = private_metadata(AUTHORITY_PRIVATE)
        terminal_private = private_metadata(TERMINAL_PRIVATE)
        public_hashes = {
            authority_public["sha256"],
            terminal_public["sha256"],
            result_public["sha256"],
            report_public["sha256"],
        }
        if len(public_hashes) != 4:
            raise BootstrapError("The four v30 role public keys are not pairwise distinct")

        receipt = {
            "schema_name": "intersubmod.v30_four_role_key_bootstrap",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "transaction_id": transaction_id,
            "status": "FRESH_KEYS_CREATED_AWAITING_AUTHORITY",
            "prepare_receipt": prepare_record["artifact"],
            "source": source,
            "openssl": openssl_record,
            "pre_authority_archive": archive_attestation,
            "failed_v8_pre_ready_archive": failed_v8_archive_attestation,
            "authority": {
                "root": str(AUTHORITY_ROOT),
                "root_metadata": authority_root_record,
                "private_key_metadata_only": authority_private,
                "private_key_bytes_recorded_in_receipt": False,
                "public_key": authority_public,
            },
            "terminal": {
                "root": str(TERMINAL_ROOT),
                "root_metadata": terminal_root_record,
                "private_key_metadata_only": terminal_private,
                "private_key_bytes_recorded_in_receipt": False,
                "public_key": terminal_public,
            },
            "result": result_attestation,
            "report": report_attestation,
            "four_role_public_keys_pairwise_distinct": True,
            "private_key_bytes_recorded_in_receipt": False,
            "release_authority_granted": False,
            "pass": True,
        }
        receipt_data = json.dumps(
            receipt, ensure_ascii=True, indent=2, sort_keys=True
        ).encode("ascii") + b"\n"
        try:
            write_exclusive(receipt_staging, receipt_data, 0o444)
            fsync_directory(receipt_staging.parent)
            exact_json_record(receipt_staging, receipt)
            rename_no_replace(receipt_staging, RECEIPT)
            fsync_directory(RECEIPT.parent)
            receipt_record = exact_json_record(RECEIPT, receipt)
        except BaseException:
            if not os.path.lexists(RECEIPT):
                raise
            receipt_record = exact_json_record(
                RECEIPT,
                receipt,
                recover_parent_durability=True,
            )
            receipt_publication_recovered = True
        ledger.append(
            "BOOTSTRAP_RECEIPT_PUBLISHED_AWAITING_SUCCESS_WITNESS",
            receipt=receipt_record["artifact"],
        )
        progress_record = ledger.finish()

        witness = {
            "schema_name": "intersubmod.v30_four_role_key_bootstrap_success",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "transaction_id": transaction_id,
            "status": "BOOTSTRAP_COMMITTED",
            "prepare_receipt": prepare_record["artifact"],
            "receipt": receipt_record["artifact"],
            "progress": progress_record,
            "pre_authority_archive": archive_attestation,
            "failed_v8_pre_ready_archive": failed_v8_archive_attestation,
            "authority_public_key": authority_public,
            "terminal_public_key": terminal_public,
            "result": result_attestation,
            "report": report_attestation,
            "checks": {
                "both_final_key_roots_published_no_replace": True,
                "both_private_keys_mode_0400_link_count_one": True,
                "all_public_keys_mode_0444_link_count_one": True,
                "four_role_public_keys_pairwise_distinct": True,
                "result_and_report_signer_ready_receipts_verified": True,
                "pre_authority_four_role_archive_verified": True,
                "failed_v8_pre_ready_archive_verified": True,
                "all_four_private_keys_mode_0400_link_count_one": True,
                "receipt_round_trip_before_witness": True,
                "progress_ledger_fsynced_before_witness": True,
                "private_key_bytes_not_recorded": True,
            },
            "release_authority_granted": False,
            "pass": True,
        }
        witness_data = json.dumps(
            witness, ensure_ascii=True, indent=2, sort_keys=True
        ).encode("ascii") + b"\n"
        try:
            write_exclusive(witness_staging, witness_data, 0o444)
            fsync_directory(witness_staging.parent)
            exact_json_record(witness_staging, witness)
            rename_no_replace(witness_staging, SUCCESS_WITNESS)
            fsync_directory(SUCCESS_WITNESS.parent)
            final_witness_record = exact_json_record(SUCCESS_WITNESS, witness)
        except BaseException:
            if not os.path.lexists(SUCCESS_WITNESS):
                raise
            final_witness_record = exact_json_record(
                SUCCESS_WITNESS,
                witness,
                recover_parent_durability=True,
            )
            witness_publication_recovered = True
        committed = True
        print(
            json.dumps(
                {
                    "prepare_receipt": prepare_record["artifact"],
                    "receipt": receipt_record["artifact"],
                    "progress": progress_record,
                    "success_witness": final_witness_record["artifact"],
                    "receipt_publication_recovered_after_atomic_rename": (
                        receipt_publication_recovered
                    ),
                    "witness_publication_recovered_after_atomic_rename": (
                        witness_publication_recovered
                    ),
                    "pre_authority_archive": archive_attestation,
                    "failed_v8_pre_ready_archive": failed_v8_archive_attestation,
                    "authority_public_key": authority_public,
                    "terminal_public_key": terminal_public,
                    "result": result_attestation,
                    "report": report_attestation,
                    "pass": True,
                },
                ensure_ascii=True,
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    except BaseException as error:
        failure_handling_errors: list[dict[str, str]] = []
        if ledger is not None and not ledger.closed:
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
            and not os.path.lexists(SUCCESS_WITNESS)
            and not os.path.lexists(PARTIAL_RECEIPT)
        ):
            partial_payload = {
                "schema_name": "intersubmod.v30_four_role_key_bootstrap_partial",
                "schema_version": "1.0.0",
                "created_at_utc": now_utc(),
                "generation": "v30",
                "transaction_id": transaction_id,
                "source": source,
                "prepare": prepare_record["artifact"],
                "progress": progress_record,
                "failure_handling_errors": failure_handling_errors,
                "pre_authority_archive": archive_attestation,
                "failed_v8_pre_ready_archive": failed_v8_archive_attestation,
                "path_states": {
                    str(path): os.path.lexists(path)
                    for path in (
                        authority_staging,
                        terminal_staging,
                        AUTHORITY_ROOT,
                        TERMINAL_ROOT,
                        receipt_staging,
                        RECEIPT,
                        witness_staging,
                        SUCCESS_WITNESS,
                    )
                },
                "error": {"type": type(error).__name__, "message": str(error)},
                "rebootstrap_forbidden_until_manual_reconciliation": True,
                "release_authority_granted": False,
                "pass": False,
            }
            publish_json_atomic(
                PARTIAL_RECEIPT,
                partial_payload,
                transaction_id=transaction_id,
                stage_role="partial",
            )
        raise
    finally:
        if ledger is not None and not ledger.closed:
            ledger.close_best_effort()
        try:
            os.close(executable_fd)
        finally:
            if previous_mask is not None:
                signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            json.dumps(
                {
                    "error": {
                        "type": type(error).__name__,
                        "message": str(error),
                    },
                    "pass": False,
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)
