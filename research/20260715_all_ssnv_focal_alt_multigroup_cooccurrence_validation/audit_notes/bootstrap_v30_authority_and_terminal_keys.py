#!/usr/bin/env python3
"""Create fresh v30 authority and terminal Ed25519 keys with no-clobber IO."""

from __future__ import annotations

from datetime import datetime, timezone
import ctypes
import hashlib
import json
import os
from pathlib import Path
import secrets
import stat
import subprocess
import sys
from typing import Any


OPENSSL = Path("/usr/bin/openssl")
EXPECTED_OPENSSL_SHA256 = (
    "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
)
AUTHORITY_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v30"
)
AUTHORITY_PRIVATE = AUTHORITY_ROOT / "ed25519_private_one_time.pem"
AUTHORITY_PUBLIC = AUTHORITY_ROOT / "ed25519_public.pem"
TERMINAL_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
    "20260724_m2v5_terminal_v20"
)
TERMINAL_PRIVATE = TERMINAL_ROOT / "ed25519_private_one_time_resident.pem"
TERMINAL_PUBLIC = TERMINAL_ROOT / "ed25519_public.pem"
EXPECTED_PARENT_MODES = {
    AUTHORITY_ROOT.parent: 0o775,
    TERMINAL_ROOT.parent: 0o700,
}
RESULT_PUBLIC = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
    "20260724_all_ssnv_result_v7_v30_recovery/ed25519_public.pem"
)
REPORT_PUBLIC = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
    "20260724_all_ssnv_report_v7_v30_recovery/ed25519_public.pem"
)
RECEIPT = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/"
    "20260724_v30_authority_terminal_key_bootstrap_receipt.json"
)
PREPARE_RECEIPT = RECEIPT.with_name(
    "20260724_v30_authority_terminal_key_bootstrap_prepare.json"
)
SUCCESS_WITNESS = RECEIPT.with_name(
    "20260724_v30_authority_terminal_key_bootstrap_success.json"
)
AT_FDCWD = -100
RENAME_NOREPLACE = 1
ENVIRONMENT = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/tmp",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
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
        or observed.st_nlink != 1
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


def json_record(path: Path) -> dict[str, Any]:
    data, observed = read_regular(path, 0o444)
    value = json.loads(data)
    if not isinstance(value, dict):
        raise BootstrapError(f"JSON root is not an object: {path}")
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


def main() -> int:
    occupied = [
        str(path)
        for path in (
            AUTHORITY_ROOT,
            TERMINAL_ROOT,
            PREPARE_RECEIPT,
            RECEIPT,
            SUCCESS_WITNESS,
        )
        if os.path.lexists(path)
    ]
    if occupied:
        raise BootstrapError(f"Bootstrap output already exists: {occupied}")
    source = source_record(Path(__file__).resolve(strict=True))
    executable_fd, openssl_record = open_bound_openssl()
    result_public = public_record(RESULT_PUBLIC)
    report_public = public_record(REPORT_PUBLIC)
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
        "schema_name": "intersubmod.v30_authority_terminal_key_bootstrap_prepare",
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
        },
        "staging_paths": {
            "authority_root": str(authority_staging),
            "terminal_root": str(terminal_staging),
            "receipt": str(receipt_staging),
            "success_witness": str(witness_staging),
        },
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
    prepare_record = json_record(PREPARE_RECEIPT)
    if prepare_record["payload"] != prepare:
        raise BootstrapError("Bootstrap prepare receipt round-trip drift")

    try:
        create_staged_keypair(
            executable_fd,
            authority_staging,
            AUTHORITY_PRIVATE.name,
        )
        create_staged_keypair(
            executable_fd,
            terminal_staging,
            TERMINAL_PRIVATE.name,
        )
        rename_no_replace(authority_staging, AUTHORITY_ROOT)
        fsync_directory(AUTHORITY_ROOT.parent)
        rename_no_replace(terminal_staging, TERMINAL_ROOT)
        fsync_directory(TERMINAL_ROOT.parent)

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
            "schema_name": "intersubmod.v30_authority_terminal_key_bootstrap",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "transaction_id": transaction_id,
            "status": "FRESH_KEYS_CREATED_AWAITING_AUTHORITY",
            "prepare_receipt": prepare_record["artifact"],
            "source": source,
            "openssl": openssl_record,
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
            "result_public_key": result_public,
            "report_public_key": report_public,
            "four_role_public_keys_pairwise_distinct": True,
            "private_key_bytes_recorded_in_receipt": False,
            "release_authority_granted": False,
            "pass": True,
        }
        receipt_data = json.dumps(
            receipt, ensure_ascii=True, indent=2, sort_keys=True
        ).encode("ascii") + b"\n"
        write_exclusive(receipt_staging, receipt_data, 0o444)
        fsync_directory(receipt_staging.parent)
        staged_receipt_record = json_record(receipt_staging)
        if staged_receipt_record["payload"] != receipt:
            raise BootstrapError("Bootstrap receipt staging round-trip drift")
        rename_no_replace(receipt_staging, RECEIPT)
        fsync_directory(RECEIPT.parent)
        receipt_record = json_record(RECEIPT)
        if receipt_record["payload"] != receipt:
            raise BootstrapError("Published bootstrap receipt round-trip drift")

        witness = {
            "schema_name": "intersubmod.v30_authority_terminal_key_bootstrap_success",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "transaction_id": transaction_id,
            "status": "BOOTSTRAP_COMMITTED",
            "prepare_receipt": prepare_record["artifact"],
            "receipt": receipt_record["artifact"],
            "authority_public_key": authority_public,
            "terminal_public_key": terminal_public,
            "result_public_key": result_public,
            "report_public_key": report_public,
            "checks": {
                "both_final_key_roots_published_no_replace": True,
                "both_private_keys_mode_0400_link_count_one": True,
                "all_public_keys_mode_0444_link_count_one": True,
                "four_role_public_keys_pairwise_distinct": True,
                "receipt_round_trip_before_witness": True,
                "private_key_bytes_not_recorded": True,
            },
            "release_authority_granted": False,
            "pass": True,
        }
        write_exclusive(
            witness_staging,
            json.dumps(witness, ensure_ascii=True, indent=2, sort_keys=True).encode(
                "ascii"
            )
            + b"\n",
            0o444,
        )
        fsync_directory(witness_staging.parent)
        staged_witness_record = json_record(witness_staging)
        if staged_witness_record["payload"] != witness:
            raise BootstrapError("Bootstrap success witness staging round-trip drift")
        rename_no_replace(witness_staging, SUCCESS_WITNESS)
        fsync_directory(SUCCESS_WITNESS.parent)
        final_witness_record = json_record(SUCCESS_WITNESS)
        if final_witness_record["payload"] != witness:
            raise BootstrapError("Published bootstrap success witness round-trip drift")
        print(
            json.dumps(
                {
                    "prepare_receipt": prepare_record["artifact"],
                    "receipt": receipt_record["artifact"],
                    "success_witness": final_witness_record["artifact"],
                    "authority_public_key": authority_public,
                    "terminal_public_key": terminal_public,
                    "result_public_key": result_public,
                    "report_public_key": report_public,
                    "pass": True,
                },
                ensure_ascii=True,
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    finally:
        os.close(executable_fd)


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
