#!/usr/bin/env python3
"""Run a bound one-time signer and publish its immutable ready attestation."""

from __future__ import annotations

from collections import deque
from datetime import datetime, timezone
import ctypes
import errno
import hashlib
import json
import os
from pathlib import Path
import secrets
import selectors
import stat
import subprocess
import sys
import time
from typing import Any


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
RESULT_ROOT = TOPIC_ROOT / "results"
WRAPPER = AUDIT_ROOT / "run_attested_one_time_ed25519_signer_v7.py"
SIGNER = AUDIT_ROOT / "run_one_time_ed25519_signer_v6.py"
EXPECTED_SIGNER = {
    "size_bytes": 34_458,
    "sha256": "c9c32222d4e8ca6daad8f6b3f8a91746fc99dd200d24cf3d77e500234c0f8bee",
    "mode": "0o444",
}
PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
PYTHON_TARGET = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
)
EXPECTED_PYTHON = {
    "size_bytes": 25_409_784,
    "sha256": "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e",
    "mode": "0o775",
}
EXPECTED_OPENSSL = {
    "path": "/usr/bin/openssl",
    "size_bytes": 1_001_272,
    "sha256": "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08",
    "mode": "0o755",
}
EXPECTED_SIGNING_PROTOCOL = "canonical SIGN JSON with exact path,size_bytes,sha256,mode"
EXPECTED_PUBLICATION_SEMANTICS = (
    "exit 0 publishes ceremony output but grants no release authority; "
    "an independent consumer must reopen and verify every live path"
)
SIGNER_ENVIRONMENT = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/tmp",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PYTHONHASHSEED": "0",
    "PYTHONNOUSERSITE": "1",
}
WRAPPER_ENVIRONMENT = {
    **SIGNER_ENVIRONMENT,
    "PYTHONDONTWRITEBYTECODE": "1",
}
AT_FDCWD = -100
RENAME_NOREPLACE = 1
MAX_PROTOCOL_LINE_BYTES = 1_048_576
INITIAL_READY_TIMEOUT_SECONDS = 30.0
COMMAND_RESPONSE_TIMEOUT_SECONDS = 120.0

ROLE_SPECS = {
    "result": {
        "key_root": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/"
            "20260724_all_ssnv_result_v9_v30_recovery"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "target": RESULT_ROOT / "task_b_final_dataset_release_receipt.v1.json",
        "prepare": AUDIT_ROOT / "20260724_v30_result_v9_signer_prepare.json",
        "ready": AUDIT_ROOT / "20260724_v30_result_v9_signer_ready.json",
        "partial": AUDIT_ROOT / "20260724_v30_result_v9_signer_partial.json",
    },
    "report": {
        "key_root": Path(
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/"
            "20260724_all_ssnv_report_v9_v30_recovery"
        ),
        "private_name": "ed25519_private_one_time.pem",
        "target": RESULT_ROOT / "task_b_final_report_release_receipt.v1.json",
        "prepare": AUDIT_ROOT / "20260724_v30_report_v9_signer_prepare.json",
        "ready": AUDIT_ROOT / "20260724_v30_report_v9_signer_ready.json",
        "partial": AUDIT_ROOT / "20260724_v30_report_v9_signer_partial.json",
    },
}


class AttestedSignerError(RuntimeError):
    """Raised before an unbound or unwitnessed signer state is accepted."""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_bound_file(
    descriptor: int,
    path: Path,
    *,
    expected_mode: int,
) -> tuple[bytes, dict[str, Any]]:
    opened = os.fstat(descriptor)
    live = os.lstat(path)
    data = b"".join(
        os.pread(descriptor, min(1024 * 1024, opened.st_size - offset), offset)
        for offset in range(0, opened.st_size, 1024 * 1024)
    )
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
        not stat.S_ISREG(opened.st_mode)
        or len(data) != opened.st_size
        or identity(opened) != identity(live)
        or stat.S_IMODE(opened.st_mode) != expected_mode
        or opened.st_nlink != 1
        or path.resolve(strict=True) != path
    ):
        raise AttestedSignerError(f"Bound file identity/mode drift: {path}")
    return data, {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
        "mode": oct(expected_mode),
        "link_count": 1,
    }


def read_regular(path: Path, *, expected_mode: int) -> tuple[bytes, dict[str, Any]]:
    descriptor = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        return read_bound_file(descriptor, path, expected_mode=expected_mode)
    finally:
        os.close(descriptor)


def private_metadata(path: Path) -> dict[str, Any]:
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
            or observed.st_size != 290
            or observed.st_nlink != 1
            or path.resolve(strict=True) != path
        ):
            raise AttestedSignerError(f"Private-key metadata drift: {path}")
        return {
            "path": str(path),
            "size_bytes": observed.st_size,
            "mode": "0o400",
            "link_count": observed.st_nlink,
        }
    finally:
        os.close(descriptor)


def write_exclusive(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    data = json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True).encode(
        "ascii"
    ) + b"\n"
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
                raise AttestedSignerError(f"Short write: {path}")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    _, record = read_regular(path, expected_mode=0o444)
    return record


def strict_json_bytes(data: bytes, label: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise AttestedSignerError(f"Duplicate JSON key in {label}: {key}")
            value[key] = item
        return value

    try:
        payload = json.loads(
            data,
            object_pairs_hook=no_duplicates,
            parse_constant=lambda value: (_ for _ in ()).throw(
                AttestedSignerError(f"Non-finite JSON in {label}: {value}")
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise AttestedSignerError(f"Invalid JSON in {label}") from error
    if not isinstance(payload, dict):
        raise AttestedSignerError(f"JSON root is not an object in {label}")
    return payload


def rename_no_replace(
    source_name: str,
    destination_name: str,
    parent_fd: int,
) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise AttestedSignerError("renameat2 is required for no-replace publication")
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        parent_fd,
        os.fsencode(source_name),
        parent_fd,
        os.fsencode(destination_name),
        RENAME_NOREPLACE,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise AttestedSignerError(
            f"No-replace publication failed: {source_name} -> {destination_name}: "
            f"{os.strerror(error_number)}"
        )


def publish_json_atomic(
    path: Path,
    payload: dict[str, Any],
    *,
    transaction_id: str,
    stage_role: str,
) -> dict[str, Any]:
    staging = path.with_name(f".{path.name}.staging.{transaction_id}.{stage_role}")
    if os.path.lexists(path) or os.path.lexists(staging):
        raise AttestedSignerError(f"Attestation publication path occupied: {path}")
    staged_record = write_exclusive(staging, payload)
    staged_data, _ = read_regular(staging, expected_mode=0o444)
    if strict_json_bytes(staged_data, f"staged {stage_role}") != payload:
        raise AttestedSignerError(f"Staged attestation round-trip drift: {stage_role}")
    parent_fd = os.open(
        path.parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
    )
    try:
        parent_before = os.fstat(parent_fd)
        rename_no_replace(staging.name, path.name, parent_fd)
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
            raise AttestedSignerError("Attestation parent identity changed")
    finally:
        os.close(parent_fd)
    final_data, final_record = read_regular(path, expected_mode=0o444)
    if (
        strict_json_bytes(final_data, f"published {stage_role}") != payload
        or final_record["sha256"] != staged_record["sha256"]
        or final_record["size_bytes"] != staged_record["size_bytes"]
    ):
        raise AttestedSignerError(f"Published attestation round-trip drift: {stage_role}")
    return final_record


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def ensure_bound_self_execution(role: str) -> None:
    binding_keys = {
        "INTERSUBMOD_ATTESTED_SIGNER_SOURCE_FD",
        "INTERSUBMOD_ATTESTED_SIGNER_PYTHON_FD",
        "INTERSUBMOD_ATTESTED_SIGNER_CANONICAL_PATH",
        "INTERSUBMOD_ATTESTED_SIGNER_PYTHON_ALIAS",
    }
    present = binding_keys.intersection(os.environ)
    if present:
        if present != binding_keys:
            raise AttestedSignerError("Incomplete wrapper binding environment")
        return
    if (
        Path(__file__).resolve(strict=True) != WRAPPER
        or Path(sys.argv[0]).resolve(strict=True) != WRAPPER
        or role not in ROLE_SPECS
    ):
        raise AttestedSignerError("Wrapper bootstrap requires canonical source and role")
    python_fd = os.open(PYTHON_TARGET, os.O_RDONLY)
    source_fd = os.open(WRAPPER, os.O_RDONLY)
    os.set_inheritable(python_fd, True)
    os.set_inheritable(source_fd, True)
    environment = {
        **WRAPPER_ENVIRONMENT,
        "INTERSUBMOD_ATTESTED_SIGNER_SOURCE_FD": str(source_fd),
        "INTERSUBMOD_ATTESTED_SIGNER_PYTHON_FD": str(python_fd),
        "INTERSUBMOD_ATTESTED_SIGNER_CANONICAL_PATH": str(WRAPPER),
        "INTERSUBMOD_ATTESTED_SIGNER_PYTHON_ALIAS": str(PYTHON),
    }
    os.execve(
        f"/proc/self/fd/{python_fd}",
        [str(PYTHON), "-I", "-B", f"/proc/self/fd/{source_fd}", role],
        environment,
    )


def verify_self_binding(role: str) -> tuple[dict[str, Any], int, dict[str, Any]]:
    source_fd = int(os.environ["INTERSUBMOD_ATTESTED_SIGNER_SOURCE_FD"])
    python_fd = int(os.environ["INTERSUBMOD_ATTESTED_SIGNER_PYTHON_FD"])
    canonical = Path(os.environ["INTERSUBMOD_ATTESTED_SIGNER_CANONICAL_PATH"])
    python_alias = os.environ["INTERSUBMOD_ATTESTED_SIGNER_PYTHON_ALIAS"]
    _, source_record = read_bound_file(source_fd, canonical, expected_mode=0o444)
    _, python_record = read_bound_file(
        python_fd, PYTHON_TARGET, expected_mode=0o775
    )
    expected_environment = {
        **WRAPPER_ENVIRONMENT,
        "INTERSUBMOD_ATTESTED_SIGNER_SOURCE_FD": str(source_fd),
        "INTERSUBMOD_ATTESTED_SIGNER_PYTHON_FD": str(python_fd),
        "INTERSUBMOD_ATTESTED_SIGNER_CANONICAL_PATH": str(WRAPPER),
        "INTERSUBMOD_ATTESTED_SIGNER_PYTHON_ALIAS": str(PYTHON),
    }
    expected_argv = [f"/proc/self/fd/{source_fd}", role]
    expected_cmdline = b"\0".join(
        os.fsencode(value)
        for value in (
            str(PYTHON),
            "-I",
            "-B",
            f"/proc/self/fd/{source_fd}",
            role,
        )
    ) + b"\0"
    observed_cmdline = Path("/proc/self/cmdline").read_bytes()
    if (
        canonical != WRAPPER
        or python_alias != str(PYTHON)
        or sys.executable != str(PYTHON)
        or Path(__file__).as_posix() != f"/proc/self/fd/{source_fd}"
        or dict(os.environ) != expected_environment
        or sys.flags.isolated != 1
        or sys.flags.no_user_site != 1
        or sys.flags.dont_write_bytecode != 1
        or sys.argv != expected_argv
        or observed_cmdline != expected_cmdline
        or python_record["size_bytes"] != EXPECTED_PYTHON["size_bytes"]
        or python_record["sha256"] != EXPECTED_PYTHON["sha256"]
    ):
        raise AttestedSignerError("Wrapper bound execution identity drift")
    return source_record, python_fd, python_record


def validate_ready(
    role: str,
    spec: dict[str, Any],
    ready: dict[str, Any],
    signer_record: dict[str, Any],
    python_record: dict[str, Any],
    parent_fd: int,
    parent_opened: os.stat_result,
) -> dict[str, Any]:
    key_root = Path(spec["key_root"])
    public_path = key_root / "ed25519_public.pem"
    private_path = key_root / str(spec["private_name"])
    public_data, public = read_regular(public_path, expected_mode=0o444)
    private = private_metadata(private_path)
    root_fd = os.open(
        key_root.name,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        dir_fd=parent_fd,
    )
    root_stat = os.fstat(root_fd)
    root_from_parent = os.stat(
        key_root.name, dir_fd=parent_fd, follow_symlinks=False
    )
    parent_after = os.fstat(parent_fd)
    bootstrap = ready.get("bootstrap", {})
    signer_artifact_keys = ("path", "size_bytes", "sha256", "mode")
    signer_public = {key: public[key] for key in signer_artifact_keys}
    signer_source = {key: signer_record[key] for key in signer_artifact_keys}
    signer_python = {key: python_record[key] for key in signer_artifact_keys}
    ready_keys = {
        "event",
        "openssl",
        "public_key",
        "private_key_mode",
        "key_directory_mode",
        "expected_target",
        "signing_protocol",
        "publication_semantics",
        "bootstrap",
    }
    bootstrap_keys = {
        "schema_name",
        "schema_version",
        "argv0",
        "script",
        "python",
        "isolated_mode",
        "no_user_site",
        "environment",
    }
    try:
        root_entries = set(os.listdir(root_fd))
        root_identity_ok = (
            root_stat.st_dev,
            root_stat.st_ino,
            root_stat.st_mode,
            root_stat.st_nlink,
        ) == (
            root_from_parent.st_dev,
            root_from_parent.st_ino,
            root_from_parent.st_mode,
            root_from_parent.st_nlink,
        )
    finally:
        os.close(root_fd)
    if (
        set(ready) != ready_keys
        or ready.get("event") != "SIGNER_READY"
        or ready.get("openssl") != EXPECTED_OPENSSL
        or ready.get("expected_target") != str(spec["target"])
        or ready.get("key_directory_mode") != "0o700"
        or ready.get("private_key_mode") != "0o400"
        or ready.get("public_key") != signer_public
        or ready.get("signing_protocol") != EXPECTED_SIGNING_PROTOCOL
        or ready.get("publication_semantics") != EXPECTED_PUBLICATION_SEMANTICS
        or len(public_data) != 113
        or not stat.S_ISDIR(root_stat.st_mode)
        or stat.S_IMODE(root_stat.st_mode) != 0o700
        or root_stat.st_uid != os.getuid()
        or not root_identity_ok
        or root_entries != {str(spec["private_name"]), "ed25519_public.pem"}
        or (
            parent_after.st_dev,
            parent_after.st_ino,
            parent_after.st_mode,
        )
        != (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        )
        or parent_after.st_nlink != parent_opened.st_nlink + 1
        or set(bootstrap) != bootstrap_keys
        or bootstrap.get("schema_name") != "intersubmod.bound_python_entrypoint"
        or bootstrap.get("schema_version") != "1.0.0"
        or not isinstance(bootstrap.get("argv0"), str)
        or not bootstrap["argv0"].startswith("/proc/self/fd/")
        or not bootstrap["argv0"].removeprefix("/proc/self/fd/").isdigit()
        or bootstrap.get("script") != signer_source
        or bootstrap.get("python") != signer_python
        or bootstrap.get("isolated_mode") is not True
        or bootstrap.get("no_user_site") is not True
        or bootstrap.get("environment") != SIGNER_ENVIRONMENT
    ):
        raise AttestedSignerError(f"Signer ready contract drift: {role}")
    return {
        "root": str(key_root),
        "root_metadata": {
            "path": str(key_root),
            "mode": "0o700",
            "owner_uid": root_stat.st_uid,
            "link_count": root_stat.st_nlink,
            "entries": sorted(root_entries),
        },
        "private_key_metadata_only": private,
        "private_key_bytes_read": False,
        "public_key": public,
        "ready_payload": ready,
        "pass": True,
    }


def read_initial_ready(
    child: subprocess.Popen[bytes],
) -> dict[str, Any]:
    assert child.stdout is not None and child.stderr is not None
    os.set_blocking(child.stdout.fileno(), False)
    os.set_blocking(child.stderr.fileno(), False)
    selector = selectors.DefaultSelector()
    selector.register(child.stdout, selectors.EVENT_READ, "stdout")
    selector.register(child.stderr, selectors.EVENT_READ, "stderr")
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    deadline = time.monotonic() + INITIAL_READY_TIMEOUT_SECONDS
    try:
        while time.monotonic() < deadline:
            for key, _ in selector.select(max(0.0, deadline - time.monotonic())):
                stream = key.fileobj
                label = key.data
                chunk = os.read(stream.fileno(), 65_536)
                if not chunk:
                    selector.unregister(stream)
                    continue
                buffers[label].extend(chunk)
                if len(buffers[label]) > MAX_PROTOCOL_LINE_BYTES:
                    raise AttestedSignerError(f"Signer initial {label} exceeded size cap")
            if b"\n" in buffers["stdout"]:
                line, remainder = bytes(buffers["stdout"]).split(b"\n", 1)
                if remainder or buffers["stderr"]:
                    raise AttestedSignerError("Signer emitted extra output before ready")
                return strict_json_bytes(line, "signer ready line")
            if child.poll() is not None and not selector.get_map():
                break
    finally:
        selector.close()
    raise AttestedSignerError(
        "Signer child returned no ready record before deadline: "
        f"exit={child.poll()} stderr_sha256={sha256_bytes(bytes(buffers['stderr']))}"
    )


def validate_terminal_line(
    line: bytes,
    *,
    stream: str,
    expected_target: Path,
) -> str:
    if stream == "stdout" and line == b"WAITING_FOR_EXACT_SIGN_TARGET":
        return "waiting"
    if stream == "stderr" and line == b"SIGNER_INPUT_CLOSED_WITHOUT_SIGNATURE":
        return "input_closed"
    payload = strict_json_bytes(line, f"signer terminal {stream}")
    event = payload.get("event")
    if stream == "stdout" and event == "CEREMONY_OUTPUT_AVAILABLE_REQUIRES_INDEPENDENT_VERIFICATION":
        if (
            set(payload)
            != {
                "event",
                "target",
                "signature",
                "private_key_mode",
                "point_in_time_checks",
                "consumer_verification_required",
                "release_authority_granted",
            }
            or payload.get("private_key_mode") != "0o0"
            or payload.get("point_in_time_checks")
            != "PASS_BEFORE_DESCRIPTOR_TEARDOWN"
            or payload.get("consumer_verification_required") is not True
            or payload.get("release_authority_granted") is not False
            or not isinstance(payload.get("target"), dict)
            or payload["target"].get("path") != str(expected_target)
            or payload["target"].get("mode") != "0o444"
            or not isinstance(payload.get("signature"), dict)
            or payload["signature"].get("path")
            != str(expected_target) + ".ed25519.sig"
            or payload["signature"].get("size_bytes") != 64
            or payload["signature"].get("mode") != "0o444"
        ):
            raise AttestedSignerError("Signer success response contract drift")
        return "success"
    if stream == "stderr" and event == "SIGN_FAILED":
        if (
            set(payload)
            != {
                "event",
                "error_type",
                "error",
                "signature_published_by_signer",
                "signature_path_exists",
                "private_retirement_started",
                "private_key_mode",
                "private_key_state",
            }
            or not isinstance(payload.get("error_type"), str)
            or not isinstance(payload.get("error"), str)
            or not isinstance(payload.get("signature_published_by_signer"), bool)
            or not isinstance(payload.get("signature_path_exists"), bool)
            or not isinstance(payload.get("private_retirement_started"), bool)
            or payload.get("private_key_mode") not in {None, "0o0", "0o400"}
            or not isinstance(payload.get("private_key_state"), str)
        ):
            raise AttestedSignerError("Signer failure response contract drift")
        return "failure"
    raise AttestedSignerError(f"Unexpected signer protocol line on {stream}")


def proxy_signer_protocol(
    child: subprocess.Popen[bytes],
    *,
    expected_target: Path,
) -> int:
    assert child.stdin is not None and child.stdout is not None and child.stderr is not None
    selector = selectors.DefaultSelector()
    parent_input = sys.stdin.buffer
    parent_fd = parent_input.fileno()
    child_stdin_fd = child.stdin.fileno()
    os.set_blocking(parent_fd, False)
    os.set_blocking(child_stdin_fd, False)
    os.set_blocking(child.stdout.fileno(), False)
    os.set_blocking(child.stderr.fileno(), False)
    manual_parent_read = stat.S_ISREG(os.fstat(parent_fd).st_mode)
    parent_registered = False
    if not manual_parent_read:
        try:
            selector.register(parent_input, selectors.EVENT_READ, "stdin")
            parent_registered = True
        except OSError as error:
            if error.errno != errno.EPERM:
                raise
            manual_parent_read = True
    selector.register(child.stdout, selectors.EVENT_READ, "stdout")
    selector.register(child.stderr, selectors.EVENT_READ, "stderr")
    output_buffers = {"stdout": bytearray(), "stderr": bytearray()}
    input_queue = bytearray()
    command_deadlines: deque[float] = deque()
    parent_partial_line_bytes = 0
    eof_deadline: float | None = None
    terminal_event: str | None = None
    parent_eof = False
    child_input_closed = False
    child_write_registered = False

    def disable_parent_read() -> None:
        nonlocal parent_registered
        if parent_registered:
            selector.unregister(parent_input)
            parent_registered = False

    def enable_child_write() -> None:
        nonlocal child_write_registered
        if not child_write_registered:
            selector.register(child.stdin, selectors.EVENT_WRITE, "child_stdin")
            child_write_registered = True

    def disable_child_write() -> None:
        nonlocal child_write_registered
        if child_write_registered:
            selector.unregister(child.stdin)
            child_write_registered = False

    def close_child_input_after_drain() -> None:
        nonlocal child_input_closed, eof_deadline
        if parent_eof and not input_queue and not child_input_closed:
            disable_child_write()
            child.stdin.close()
            child_input_closed = True
            if eof_deadline is None:
                eof_deadline = time.monotonic() + COMMAND_RESPONSE_TIMEOUT_SECONDS

    def accept_parent_chunk(chunk: bytes) -> None:
        nonlocal parent_eof, parent_partial_line_bytes, eof_deadline
        if not chunk:
            observed_at = time.monotonic()
            if parent_partial_line_bytes:
                command_deadlines.append(
                    observed_at + COMMAND_RESPONSE_TIMEOUT_SECONDS
                )
                parent_partial_line_bytes = 0
            if eof_deadline is None:
                eof_deadline = observed_at + COMMAND_RESPONSE_TIMEOUT_SECONDS
            parent_eof = True
            disable_parent_read()
            close_child_input_after_drain()
            return
        if len(input_queue) + len(chunk) > MAX_PROTOCOL_LINE_BYTES:
            raise AttestedSignerError("Parent signer input exceeded bounded queue")
        pieces = chunk.split(b"\n")
        if len(pieces) == 1:
            parent_partial_line_bytes += len(chunk)
            if parent_partial_line_bytes > MAX_PROTOCOL_LINE_BYTES:
                raise AttestedSignerError("Parent signer protocol line exceeded size cap")
        else:
            if parent_partial_line_bytes + len(pieces[0]) > MAX_PROTOCOL_LINE_BYTES:
                raise AttestedSignerError("Parent signer protocol line exceeded size cap")
            if any(len(piece) > MAX_PROTOCOL_LINE_BYTES for piece in pieces[1:-1]):
                raise AttestedSignerError("Parent signer protocol line exceeded size cap")
            parent_partial_line_bytes = len(pieces[-1])
        input_queue.extend(chunk)
        observed_at = time.monotonic()
        for _ in range(chunk.count(b"\n")):
            command_deadlines.append(observed_at + COMMAND_RESPONSE_TIMEOUT_SECONDS)
        enable_child_write()

    try:
        while True:
            if manual_parent_read and not parent_eof:
                try:
                    accept_parent_chunk(os.read(parent_fd, 65_536))
                except BlockingIOError:
                    pass
            deadline_candidates = list(command_deadlines)
            if eof_deadline is not None:
                deadline_candidates.append(eof_deadline)
            oldest_deadline = min(deadline_candidates, default=None)
            if oldest_deadline is not None and time.monotonic() >= oldest_deadline:
                raise AttestedSignerError("Signer command response timed out")
            timeout = (
                1.0
                if oldest_deadline is None
                else max(0.0, min(1.0, oldest_deadline - time.monotonic()))
            )
            for key, _ in selector.select(timeout):
                stream = key.fileobj
                label = key.data
                if label == "stdin":
                    try:
                        accept_parent_chunk(os.read(parent_fd, 65_536))
                    except BlockingIOError:
                        pass
                    continue
                if label == "child_stdin":
                    try:
                        written = os.write(child_stdin_fd, input_queue)
                    except BlockingIOError:
                        continue
                    except BrokenPipeError as error:
                        raise AttestedSignerError(
                            "Signer child closed stdin with queued input"
                        ) from error
                    if written <= 0:
                        raise AttestedSignerError("Short nonblocking signer-input write")
                    del input_queue[:written]
                    if not input_queue:
                        disable_child_write()
                        close_child_input_after_drain()
                    continue
                chunk = os.read(stream.fileno(), 65_536)
                if not chunk:
                    selector.unregister(stream)
                    continue
                destination = sys.stdout.buffer if label == "stdout" else sys.stderr.buffer
                destination.write(chunk)
                destination.flush()
                output_buffers[label].extend(chunk)
                if len(output_buffers[label]) > MAX_PROTOCOL_LINE_BYTES:
                    raise AttestedSignerError(f"Signer {label} protocol line exceeded size cap")
                while b"\n" in output_buffers[label]:
                    raw_line, remainder = bytes(output_buffers[label]).split(b"\n", 1)
                    output_buffers[label] = bytearray(remainder)
                    event = validate_terminal_line(
                        raw_line,
                        stream=label,
                        expected_target=expected_target,
                    )
                    if event in {"waiting", "success", "failure"}:
                        if command_deadlines:
                            command_deadlines.popleft()
                        elif not (event == "waiting" and parent_eof):
                            raise AttestedSignerError(
                                "Signer emitted a response without a complete command"
                            )
                    if event in {"success", "failure", "input_closed"}:
                        if terminal_event is not None:
                            raise AttestedSignerError("Signer emitted multiple terminal events")
                        if event == "input_closed" and not parent_eof:
                            raise AttestedSignerError(
                                "Signer reported input closure before parent EOF"
                            )
                        terminal_event = event
            exit_code = child.poll()
            child_streams_open = any(
                key.data in {"stdout", "stderr"} for key in selector.get_map().values()
            )
            if exit_code is not None and not child_streams_open:
                if any(output_buffers.values()):
                    raise AttestedSignerError("Signer ended with unterminated protocol output")
                expected_event = {0: "success", 7: "input_closed", 8: "failure"}.get(
                    exit_code
                )
                if expected_event is None or terminal_event != expected_event:
                    raise AttestedSignerError(
                        f"Signer exit/event mismatch: exit={exit_code} event={terminal_event}"
                    )
                return exit_code
    finally:
        selector.close()


def stop_and_reap(child: subprocess.Popen[bytes]) -> None:
    if child.poll() is not None:
        child.wait()
        return
    child.terminate()
    try:
        child.wait(timeout=5)
    except subprocess.TimeoutExpired:
        child.kill()
        child.wait(timeout=5)


def exception_record(stage: str, error: BaseException) -> dict[str, str]:
    return {
        "stage": stage,
        "type": type(error).__name__,
        "message": str(error),
    }


def reap_child_for_attestation(
    child: subprocess.Popen[bytes] | None,
) -> tuple[dict[str, Any], list[tuple[str, BaseException]]]:
    failures: list[tuple[str, BaseException]] = []
    if child is None:
        return {
            "started": False,
            "pid": None,
            "reaped": True,
            "return_code": None,
            "recovery_kill_attempted": False,
        }, failures

    recovery_kill_attempted = False
    try:
        stop_and_reap(child)
    except BaseException as error:
        failures.append(("stop_and_reap", error))

    if child.returncode is None:
        recovery_kill_attempted = True
        try:
            child.kill()
        except BaseException as error:
            failures.append(("recovery_kill", error))
        try:
            child.wait(timeout=5)
        except BaseException as error:
            failures.append(("recovery_wait", error))

    if child.returncode is None:
        failures.append(
            (
                "verify_child_reaped",
                AttestedSignerError(
                    "Signer child was not confirmed reaped after forced recovery"
                ),
            )
        )

    return {
        "started": True,
        "pid": child.pid,
        "reaped": child.returncode is not None,
        "return_code": child.returncode,
        "recovery_kill_attempted": recovery_kill_attempted,
    }, failures


def close_attestation_descriptors(
    descriptors: dict[str, int],
) -> tuple[dict[str, str], list[tuple[str, BaseException]]]:
    states: dict[str, str] = {}
    failures: list[tuple[str, BaseException]] = []
    for label, descriptor in descriptors.items():
        if descriptor < 0:
            states[label] = "NOT_OPENED"
            continue
        try:
            os.close(descriptor)
        except BaseException as error:
            states[label] = "CLOSE_FAILED_UNKNOWN_STATE"
            failures.append((f"close_{label}", error))
        else:
            states[label] = "CLOSED"
    return states, failures


def sample_key_state_after_cleanup(
    key_root: Path,
    private_name: str,
) -> dict[str, Any]:
    key_state: dict[str, Any] = {
        "root": str(key_root),
        "root_present": os.path.lexists(key_root),
        "sampled_after_cleanup": True,
        "private_key_bytes_read": False,
    }
    if key_root.is_dir():
        root_stat = os.lstat(key_root)
        key_state["root_metadata"] = {
            "mode": oct(stat.S_IMODE(root_stat.st_mode)),
            "owner_uid": root_stat.st_uid,
            "link_count": root_stat.st_nlink,
            "entries": sorted(path.name for path in key_root.iterdir()),
        }
        public_path = key_root / "ed25519_public.pem"
        private_path = key_root / private_name
        if public_path.is_file():
            _, key_state["public_key"] = read_regular(
                public_path, expected_mode=0o444
            )
        if private_path.is_file():
            key_state["private_key_metadata_only"] = private_metadata(private_path)
    return key_state


def run(role: str) -> int:
    ensure_bound_self_execution(role)
    source_record, python_fd, python_record = verify_self_binding(role)
    spec = ROLE_SPECS[role]
    key_root = Path(spec["key_root"])
    prepare_path = Path(spec["prepare"])
    ready_path = Path(spec["ready"])
    partial_path = Path(spec["partial"])
    target = Path(spec["target"])
    signature = Path(str(target) + ".ed25519.sig")
    parent = key_root.parent.resolve(strict=True)
    signer_fd = -1
    parent_fd = -1
    child: subprocess.Popen[bytes] | None = None
    prepare_record: dict[str, Any] | None = None
    ready_record: dict[str, Any] | None = None
    transaction_id: str | None = None
    try:
        if any(
            os.path.lexists(path)
            for path in (
                key_root,
                prepare_path,
                ready_path,
                partial_path,
                target,
                signature,
            )
        ):
            raise AttestedSignerError(f"Attested signer output already exists: {role}")
        parent_fd = os.open(
            parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
        )
        parent_opened = os.fstat(parent_fd)
        parent_live = os.lstat(parent)
        if (
            parent != key_root.parent
            or not stat.S_ISDIR(parent_opened.st_mode)
            or (
                parent_opened.st_dev,
                parent_opened.st_ino,
                parent_opened.st_mode,
                parent_opened.st_nlink,
            )
            != (
                parent_live.st_dev,
                parent_live.st_ino,
                parent_live.st_mode,
                parent_live.st_nlink,
            )
            or parent_opened.st_uid != os.getuid()
            or stat.S_IMODE(parent_opened.st_mode) & 0o002
        ):
            raise AttestedSignerError(f"Signer key parent is unsafe: {parent}")
        try:
            os.stat(key_root.name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise AttestedSignerError(f"Signer key root is occupied: {key_root}")
        signer_fd = os.open(SIGNER, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        _, signer_record = read_bound_file(signer_fd, SIGNER, expected_mode=0o444)
        if any(signer_record.get(key) != value for key, value in EXPECTED_SIGNER.items()):
            raise AttestedSignerError("Frozen signer source identity drift")
        transaction_id = secrets.token_hex(16)
        prepare_payload = {
            "schema_name": "intersubmod.attested_one_time_signer_prepare",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "role": role,
            "transaction_id": transaction_id,
            "wrapper_source": source_record,
            "signer_source": signer_record,
            "python": python_record,
            "key_root": str(key_root),
            "key_parent": {
                "path": str(parent),
                "device": parent_opened.st_dev,
                "inode": parent_opened.st_ino,
                "mode": oct(stat.S_IMODE(parent_opened.st_mode)),
                "owner_uid": parent_opened.st_uid,
                "link_count": parent_opened.st_nlink,
            },
            "key_root_absent_before_child_launch": True,
            "expected_target": str(target),
            "expected_target_and_signature_absent": True,
            "status": "PREPARED_NOT_READY",
            "release_authority_granted": False,
            "pass": False,
        }
        prepare_record = publish_json_atomic(
            prepare_path,
            prepare_payload,
            transaction_id=transaction_id,
            stage_role="prepare",
        )
        try:
            os.stat(key_root.name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise AttestedSignerError("Key root appeared before child launch")
        child = subprocess.Popen(
            [
                str(PYTHON),
                "-I",
                f"/proc/self/fd/{signer_fd}",
                str(key_root),
                str(spec["private_name"]),
                str(target),
            ],
            executable=f"/proc/self/fd/{python_fd}",
            cwd=REPO_ROOT,
            env=SIGNER_ENVIRONMENT,
            pass_fds=(python_fd, signer_fd),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if child.stdin is None or child.stdout is None or child.stderr is None:
            raise AttestedSignerError("Signer child pipe setup failed")
        ready = read_initial_ready(child)
        validated = validate_ready(
            role,
            spec,
            ready,
            signer_record,
            python_record,
            parent_fd,
            parent_opened,
        )
        parent_after_ready = os.fstat(parent_fd)
        if (
            parent_after_ready.st_dev,
            parent_after_ready.st_ino,
            parent_after_ready.st_mode,
        ) != (
            parent_opened.st_dev,
            parent_opened.st_ino,
            parent_opened.st_mode,
        ) or parent_after_ready.st_nlink != parent_opened.st_nlink + 1:
            raise AttestedSignerError("Key parent changed during signer setup")
        ready_payload = {
            "schema_name": "intersubmod.attested_one_time_signer_ready",
            "schema_version": "1.0.0",
            "created_at_utc": now_utc(),
            "generation": "v30",
            "role": role,
            "transaction_id": transaction_id,
            "status": "READY_AWAITING_EXACT_SIGN_TARGET",
            "prepare": prepare_record,
            "wrapper_source": source_record,
            "signer_source": signer_record,
            "python": python_record,
            "child_pid_at_ready": child.pid,
            "creation_contract": {
                "key_root_absence_observation": {
                    "parent_device": parent_opened.st_dev,
                    "parent_inode": parent_opened.st_ino,
                    "root_name": key_root.name,
                    "present": False,
                    "observed_before_child_launch": True,
                },
                "key_parent_identity_transition_verified": {
                    "path": str(parent),
                    "device": parent_opened.st_dev,
                    "inode": parent_opened.st_ino,
                    "mode": oct(stat.S_IMODE(parent_opened.st_mode)),
                    "link_count_before": parent_opened.st_nlink,
                    "link_count_after": parent_after_ready.st_nlink,
                    "expected_link_count_delta_for_one_new_subdirectory": 1,
                },
                "signer_exclusive_creation_source_bound": signer_record,
                "root_exact_inventory_verified": validated["root_metadata"],
            },
            "key_material": validated,
            "expected_target": str(target),
            "expected_target_and_signature_absent_at_ready": (
                not os.path.lexists(target) and not os.path.lexists(signature)
            ),
            "release_authority_granted": False,
            "pass": True,
        }
        if ready_payload["expected_target_and_signature_absent_at_ready"] is not True:
            raise AttestedSignerError("Target appeared before signer ready publication")
        ready_record = publish_json_atomic(
            ready_path,
            ready_payload,
            transaction_id=transaction_id,
            stage_role="ready",
        )
        print(
            json.dumps(
                {
                    "event": "ATTESTED_SIGNER_READY",
                    "role": role,
                    "ready_receipt": ready_record,
                    "public_key": validated["public_key"],
                    "private_key_metadata_only": validated[
                        "private_key_metadata_only"
                    ],
                    "release_authority_granted": False,
                    "pass": True,
                },
                ensure_ascii=True,
                sort_keys=True,
            ),
            flush=True,
        )
        return proxy_signer_protocol(child, expected_target=target)
    finally:
        active_error = sys.exc_info()[1]
        child_state, cleanup_failures = reap_child_for_attestation(child)
        descriptor_states, descriptor_failures = close_attestation_descriptors(
            {
                "signer_source": signer_fd,
                "key_parent": parent_fd,
                "python_runtime": python_fd,
            }
        )
        cleanup_failures.extend(descriptor_failures)
        descriptors_settled = all(
            state in {"CLOSED", "NOT_OPENED"}
            for state in descriptor_states.values()
        )
        cleanup_record = {
            "completed_before_key_state_sample": (
                child_state["reaped"] is True and descriptors_settled
            ),
            "child": child_state,
            "descriptors": descriptor_states,
            "errors": [
                exception_record(stage, error)
                for stage, error in cleanup_failures
            ],
        }
        failure_requires_partial = (
            active_error is not None or bool(cleanup_failures)
        )
        if prepare_record is not None and failure_requires_partial:
            attestation_failures: list[dict[str, str]] = []
            if cleanup_record["completed_before_key_state_sample"] is True:
                try:
                    key_state = sample_key_state_after_cleanup(
                        key_root, str(spec["private_name"])
                    )
                except BaseException as error:
                    attestation_failures.append(
                        exception_record("sample_key_state_after_cleanup", error)
                    )
                    key_state = {
                        "root": str(key_root),
                        "sampling_performed": False,
                        "sampling_attempted_after_cleanup": True,
                        "private_key_bytes_read": False,
                    }
            else:
                key_state = {
                    "root": str(key_root),
                    "sampling_performed": False,
                    "sampling_attempted_after_cleanup": False,
                    "reason": "CLEANUP_NOT_CONFIRMED_COMPLETE",
                    "private_key_bytes_read": False,
                }
            primary_error = (
                exception_record("run", active_error)
                if active_error is not None
                else cleanup_record["errors"][0]
            )
            partial_payload = {
                "schema_name": "intersubmod.attested_one_time_signer_partial",
                "schema_version": "1.0.0",
                "created_at_utc": now_utc(),
                "generation": "v30",
                "role": role,
                "prepare": prepare_record,
                "ready": ready_record,
                "wrapper_source": source_record,
                "cleanup": cleanup_record,
                "key_state": key_state,
                "attestation_errors": attestation_failures,
                "expected_target": str(target),
                "target_present": os.path.lexists(target),
                "signature_present": os.path.lexists(signature),
                "error": primary_error,
                "rebootstrap_same_key_root_forbidden": True,
                "release_authority_granted": False,
                "pass": False,
            }
            try:
                if transaction_id is None:
                    raise AttestedSignerError(
                        "Prepared signer failure has no transaction identifier"
                    )
                publish_json_atomic(
                    partial_path,
                    partial_payload,
                    transaction_id=transaction_id,
                    stage_role="partial",
                )
            except BaseException as publication_error:
                original = (
                    exception_record("run", active_error)
                    if active_error is not None
                    else primary_error
                )
                raise AttestedSignerError(
                    "Partial signer attestation publication failed after cleanup: "
                    f"original={original} publication="
                    f"{exception_record('publish_partial', publication_error)}"
                ) from publication_error
        if active_error is None and cleanup_failures:
            first_stage, first_error = cleanup_failures[0]
            raise AttestedSignerError(
                f"Signer cleanup failed after protocol: {first_stage}: {first_error}"
            ) from first_error


def main() -> int:
    if len(sys.argv) != 2 or sys.argv[1] not in ROLE_SPECS:
        raise AttestedSignerError("Usage: wrapper result|report")
    return run(sys.argv[1])


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
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
        raise SystemExit(3)
