#!/usr/bin/env python3
"""Publish the reviewed recovery v30 authority as one atomic directory bundle."""

from __future__ import annotations

import ctypes
from datetime import datetime, timezone
import errno
import fcntl
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import secrets
import signal
import stat
import subprocess
import sys
from types import ModuleType
from typing import Any, Callable, Mapping


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
CEREMONY_BUILDER = AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v30.py"
VALIDATOR = AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v30.py"
EXPECTED_VALIDATOR_SHA256 = (
    "74c3e445fec567cf81184611b983d260b1b60344b8def63fd2feac099004daf1"
)

MEMBERS = ["authority.ed25519.sig", "authority.json", "commit.json"]
RENAME_NOREPLACE = 1
LINUX_MFD_CLOEXEC = getattr(os, "MFD_CLOEXEC", 0x0001)
LINUX_MFD_ALLOW_SEALING = getattr(os, "MFD_ALLOW_SEALING", 0x0002)
LINUX_F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
LINUX_F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
LINUX_F_SEAL_SEAL = getattr(fcntl, "F_SEAL_SEAL", 0x0001)
LINUX_F_SEAL_SHRINK = getattr(fcntl, "F_SEAL_SHRINK", 0x0002)
LINUX_F_SEAL_GROW = getattr(fcntl, "F_SEAL_GROW", 0x0004)
LINUX_F_SEAL_WRITE = getattr(fcntl, "F_SEAL_WRITE", 0x0008)

_RENAME_LIBC = ctypes.PyDLL(None, use_errno=True)
_RENAMEAT2 = getattr(_RENAME_LIBC, "renameat2", None)
if _RENAMEAT2 is not None:
    _RENAMEAT2.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    _RENAMEAT2.restype = ctypes.c_int

_MEMFD_CREATE = getattr(_RENAME_LIBC, "memfd_create", None)
if _MEMFD_CREATE is not None:
    _MEMFD_CREATE.argtypes = [ctypes.c_char_p, ctypes.c_uint]
    _MEMFD_CREATE.restype = ctypes.c_int

_INOTIFY_INIT1 = getattr(_RENAME_LIBC, "inotify_init1", None)
_INOTIFY_ADD_WATCH = getattr(_RENAME_LIBC, "inotify_add_watch", None)
if _INOTIFY_INIT1 is not None:
    _INOTIFY_INIT1.argtypes = [ctypes.c_int]
    _INOTIFY_INIT1.restype = ctypes.c_int
if _INOTIFY_ADD_WATCH is not None:
    _INOTIFY_ADD_WATCH.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32]
    _INOTIFY_ADD_WATCH.restype = ctypes.c_int

INOTIFY_MUTATION_MASK = (
    0x00000002  # IN_MODIFY
    | 0x00000004  # IN_ATTRIB
    | 0x00000008  # IN_CLOSE_WRITE
    | 0x00000040  # IN_MOVED_FROM
    | 0x00000080  # IN_MOVED_TO
    | 0x00000100  # IN_CREATE
    | 0x00000200  # IN_DELETE
    | 0x00000400  # IN_DELETE_SELF
    | 0x00000800  # IN_MOVE_SELF
    | 0x00004000  # IN_Q_OVERFLOW
    | 0x00008000  # IN_IGNORED
)


def stable_json(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    ).encode("ascii")


def read_fd(descriptor: int, size: int) -> bytes:
    return b"".join(
        os.pread(descriptor, min(8 * 1024 * 1024, size - offset), offset)
        for offset in range(0, size, 8 * 1024 * 1024)
    )


def full_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


class BoundInput:
    def __init__(
        self,
        path: Path,
        descriptor: int,
        opened: os.stat_result,
        data: bytes,
    ) -> None:
        self.path = path
        self.descriptor = descriptor
        self.opened = opened
        self.data = data
        self.sha256 = hashlib.sha256(data).hexdigest()

    @classmethod
    def open(cls, path: Path, *, expected_mode: str) -> "BoundInput":
        if not path.is_absolute() or path.resolve(strict=True) != path:
            raise RuntimeError(f"Bound input path is not canonical: {path}")
        before = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or oct(stat.S_IMODE(before.st_mode)) != expected_mode
        ):
            raise RuntimeError(f"Bound input metadata drift: {path}")
        flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        data = read_fd(descriptor, opened.st_size)
        after = os.stat(path, follow_symlinks=False)
        if (
            full_identity(before) != full_identity(opened)
            or full_identity(opened) != full_identity(after)
            or len(data) != opened.st_size
        ):
            os.close(descriptor)
            raise RuntimeError(f"Bound input changed while opening: {path}")
        return cls(path, descriptor, opened, data)

    def record(self, *, mode: str | None = None) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "size_bytes": len(self.data),
            "sha256": self.sha256,
            "mode": mode or oct(stat.S_IMODE(self.opened.st_mode)),
        }

    def recheck(self, *, expected_mode: str, allow_retirement: bool = False) -> None:
        opened = os.fstat(self.descriptor)
        live = os.stat(self.path, follow_symlinks=False)
        if allow_retirement:
            stable = lambda value: (
                value.st_dev,
                value.st_ino,
                value.st_nlink,
                value.st_size,
                value.st_mtime_ns,
            )
            valid_identity = (
                stable(opened) == stable(self.opened)
                and stable(live) == stable(self.opened)
                and opened.st_ctime_ns >= self.opened.st_ctime_ns
                and live.st_ctime_ns == opened.st_ctime_ns
            )
        else:
            valid_identity = (
                full_identity(opened) == full_identity(self.opened)
                and full_identity(live) == full_identity(self.opened)
            )
        data = read_fd(self.descriptor, opened.st_size)
        if (
            not valid_identity
            or oct(stat.S_IMODE(opened.st_mode)) != expected_mode
            or oct(stat.S_IMODE(live.st_mode)) != expected_mode
            or len(data) != len(self.data)
            or hashlib.sha256(data).hexdigest() != self.sha256
        ):
            raise RuntimeError(f"Bound input terminal recheck failed: {self.path}")

    def close(self) -> None:
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1


def bind_executing_builder() -> tuple[BoundInput, dict[str, Any]]:
    canonical = str(CEREMONY_BUILDER)
    source_path = Path(__file__).resolve(strict=True)
    argv_path = Path(sys.argv[0]).resolve(strict=True)
    cmdline = Path("/proc/self/cmdline").read_bytes().split(b"\0")
    canonical_token = os.fsencode(canonical)
    if (
        source_path != CEREMONY_BUILDER
        or argv_path != CEREMONY_BUILDER
        or sys.argv[0] != canonical
        or cmdline.count(canonical_token) != 1
    ):
        raise RuntimeError(
            "Ceremony builder must execute from its exact canonical reviewed path"
        )
    bound = BoundInput.open(CEREMONY_BUILDER, expected_mode="0o444")
    live = os.stat(source_path, follow_symlinks=False)
    if full_identity(live) != full_identity(bound.opened):
        bound.close()
        raise RuntimeError("Executing ceremony builder source/path identity drift")
    return bound, {
        "source": bound.record(),
        "canonical_file": canonical,
        "sys_argv0": sys.argv[0],
        "proc_cmdline_contains_exact_canonical_script_once": True,
        "bound_before_validator_load": True,
        "same_path_inode_as_bound_source": True,
        "threat_model": "trusted_same_uid_account_no_malicious_runtime_code_injection",
        "pass": True,
    }


class WriterLease:
    def __init__(self, bound: BoundInput) -> None:
        self.bound = bound

    @classmethod
    def acquire(cls, module: ModuleType) -> "WriterLease":
        bound = BoundInput.open(module.PUBLIC_KEY, expected_mode="0o444")
        try:
            if bound.sha256 != module.EXPECTED_PUBLIC_KEY_SHA256:
                raise RuntimeError("Recovery writer-lock public-key identity drift")
            try:
                fcntl.flock(bound.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RuntimeError(
                    "Another cooperating recovery output writer holds the v30 lease"
                ) from error
            bound.recheck(expected_mode="0o444")
            return cls(bound)
        except BaseException:
            bound.close()
            raise

    def recheck(self) -> None:
        self.bound.recheck(expected_mode="0o444")

    def close(self) -> None:
        if self.bound.descriptor >= 0:
            fcntl.flock(self.bound.descriptor, fcntl.LOCK_UN)
            self.bound.close()


class DirectoryLease:
    def __init__(self, path: Path) -> None:
        if not path.is_absolute() or path.resolve(strict=True) != path:
            raise RuntimeError(f"Directory path is not canonical: {path}")
        before = os.stat(path, follow_symlinks=False)
        flags = os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        self.descriptor = os.open(path, flags)
        self.opened = os.fstat(self.descriptor)
        self.path = path
        if (
            not stat.S_ISDIR(before.st_mode)
            or (before.st_dev, before.st_ino, before.st_mode)
            != (self.opened.st_dev, self.opened.st_ino, self.opened.st_mode)
        ):
            self.close()
            raise RuntimeError(f"Directory changed while binding: {path}")

    def recheck(self) -> None:
        opened = os.fstat(self.descriptor)
        live = os.stat(self.path, follow_symlinks=False)
        identity = lambda value: (value.st_dev, value.st_ino, value.st_mode)
        if identity(opened) != identity(self.opened) or identity(live) != identity(
            self.opened
        ):
            raise RuntimeError(f"Directory lease drift: {self.path}")

    def generation(self) -> tuple[int, ...]:
        opened = os.fstat(self.descriptor)
        live = os.stat(self.path, follow_symlinks=False)
        if full_identity(opened) != full_identity(live):
            raise RuntimeError(f"Directory generation changed while sampled: {self.path}")
        return full_identity(opened)

    def close(self) -> None:
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1


class DirectoryMutationWatch:
    def __init__(self, leases: Mapping[Path, DirectoryLease]) -> None:
        if _INOTIFY_INIT1 is None or _INOTIFY_ADD_WATCH is None:
            raise RuntimeError("inotify mutation-watch primitives are unavailable")
        descriptor = _INOTIFY_INIT1(os.O_NONBLOCK | os.O_CLOEXEC)
        if descriptor < 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error), "inotify_init1")
        self.descriptor = descriptor
        self.paths = frozenset(leases)
        try:
            for path in sorted(leases, key=str):
                watch = _INOTIFY_ADD_WATCH(
                    descriptor, os.fsencode(path), INOTIFY_MUTATION_MASK
                )
                if watch < 0:
                    error = ctypes.get_errno()
                    raise OSError(error, os.strerror(error), str(path))
                leases[path].recheck()
        except BaseException:
            self.close()
            raise

    def assert_quiet(self) -> None:
        observed = bytearray()
        while True:
            try:
                block = os.read(self.descriptor, 64 * 1024)
            except BlockingIOError:
                break
            if not block:
                break
            observed.extend(block)
        if observed:
            raise RuntimeError(
                "Ceremony namespace mutation event observed during absence scan"
            )

    def close(self) -> None:
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1


def open_absence_directory_leases(module: ModuleType) -> dict[Path, DirectoryLease]:
    parents = (module.RESULT_ROOT, module.REVIEW_ROOT, module.WORKSPACE_ROOT)
    leases: dict[Path, DirectoryLease] = {}
    try:
        for path in parents:
            leases[path] = DirectoryLease(path)
    except BaseException:
        for lease in leases.values():
            lease.close()
        raise
    return leases


def require_ceremony_absence(
    module: ModuleType,
    leases: Mapping[Path, DirectoryLease],
    *,
    allowed_stage_name: str | None = None,
    mutation_watch: DirectoryMutationWatch | None = None,
) -> dict[str, Any]:
    expected_parents = {module.RESULT_ROOT, module.REVIEW_ROOT, module.WORKSPACE_ROOT}
    contract = module._validate_ceremony_absence_contract()
    if (
        set(leases) != expected_parents
        or contract.get("forbidden_output_slot_count")
        != module.EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT
    ):
        raise RuntimeError("Ceremony absence lease/contract drift")

    ordered_parents = tuple(sorted(leases, key=str))

    def generations() -> dict[Path, tuple[int, ...]]:
        return {path: leases[path].generation() for path in ordered_parents}

    def scan_namespace() -> tuple[list[str], list[str]]:
        occupied: list[str] = []
        for path in module.CEREMONY_FORBIDDEN_OUTPUT_SLOTS:
            lease = leases.get(path.parent)
            if lease is None:
                raise RuntimeError(f"Unleased forbidden-output parent: {path}")
            try:
                os.stat(path.name, dir_fd=lease.descriptor, follow_symlinks=False)
            except FileNotFoundError as error:
                if error.errno != errno.ENOENT:
                    raise
            else:
                occupied.append(str(path))
        result_lease = leases[module.RESULT_ROOT]
        unexpected_staging = sorted(
            name
            for name in os.listdir(result_lease.descriptor)
            if any(
                fnmatch.fnmatchcase(name, pattern)
                for pattern in module.CEREMONY_STAGING_PATTERNS
            )
            and name != allowed_stage_name
        )
        return occupied, unexpected_staging

    owns_watch = mutation_watch is None
    watch = mutation_watch or DirectoryMutationWatch(leases)
    if watch.descriptor < 0 or watch.paths != frozenset(leases):
        raise RuntimeError("Ceremony mutation-watch lease set drift")
    try:
        generation_t0 = generations()
        occupied_first, staging_first = scan_namespace()
        generation_t1 = generations()
        if generation_t1 != generation_t0:
            raise RuntimeError("Ceremony namespace generation drift during first scan")
        occupied_second, staging_second = scan_namespace()
        generation_t2 = generations()
        if generation_t2 != generation_t1:
            raise RuntimeError("Ceremony namespace generation drift during second scan")
        watch.assert_quiet()

        occupied = sorted(set(occupied_first) | set(occupied_second))
        unexpected_staging = sorted(set(staging_first) | set(staging_second))
        if occupied or unexpected_staging:
            raise RuntimeError(
                "Ceremony forbidden namespace became occupied: "
                f"slots={occupied}, staging={unexpected_staging}"
            )
        return {
            "forbidden_output_slots_checked": len(
                module.CEREMONY_FORBIDDEN_OUTPUT_SLOTS
            ),
            "staging_patterns_checked": len(module.CEREMONY_STAGING_PATTERNS),
            "allowed_stage_name": allowed_stage_name,
            "directory_generation_tokens": {
                str(path): list(generation_t2[path]) for path in ordered_parents
            },
            "scan_count": 2,
            "inotify_mutation_watch": "quiet",
            "pass": True,
        }
    finally:
        if owns_watch:
            watch.close()


def load_validator() -> tuple[ModuleType, BoundInput]:
    lease = BoundInput.open(VALIDATOR, expected_mode="0o444")
    if lease.sha256 != EXPECTED_VALIDATOR_SHA256:
        lease.close()
        raise RuntimeError("Frozen recovery validator SHA-256 drift")
    module = ModuleType("bound_tumor_ref_schema_recovery_authority_builder_v30")
    module.__file__ = str(VALIDATOR)
    module.__package__ = ""
    code = compile(lease.data, str(VALIDATOR), "exec", dont_inherit=True)
    exec(code, module.__dict__)
    return module, lease


def require_final_slots_absent(module: ModuleType) -> None:
    occupied = [
        str(path)
        for path in (module.AUTHORITY_BUNDLE,)
        if os.path.lexists(path)
    ]
    if occupied:
        raise RuntimeError(f"Recovery authority final slot is occupied: {occupied}")


def require_no_staging(module: ModuleType) -> None:
    pattern = f".{module.AUTHORITY_BUNDLE.name}.staging.*"
    staging = sorted(str(path) for path in module.RESULT_ROOT.glob(pattern))
    if staging:
        raise RuntimeError(f"Recovery authority staging evidence exists: {staging}")


def run_probe(module: ModuleType) -> tuple[dict[str, Any], dict[str, Any]]:
    probe = BoundInput.open(module.READONLY_PROBE, expected_mode="0o444")
    python_target_path = module.PYTHON.resolve(strict=True)
    python_target = BoundInput.open(python_target_path, expected_mode="0o775")
    environment = {
        "PATH": "/usr/bin:/bin",
        "HOME": "/tmp",
        "LANG": "C.UTF-8",
        "LC_ALL": "C.UTF-8",
        "PYTHONHASHSEED": "0",
        "PYTHONNOUSERSITE": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "INTERSUBMOD_PROBE_SOURCE_FD": str(probe.descriptor),
        "INTERSUBMOD_PROBE_CANONICAL_PATH": str(module.READONLY_PROBE),
        "INTERSUBMOD_PROBE_PYTHON_FD": str(python_target.descriptor),
        "INTERSUBMOD_PROBE_PYTHON_ALIAS": str(module.PYTHON),
    }
    launch_script = 'exec -a "$1" "/proc/self/fd/$2" -I -B "/proc/self/fd/$3"'
    command = [
        "/bin/bash",
        "-c",
        launch_script,
        "bound-probe-launcher",
        str(module.PYTHON),
        str(python_target.descriptor),
        str(probe.descriptor),
    ]
    try:
        result = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=environment,
            pass_fds=(python_target.descriptor, probe.descriptor),
            capture_output=True,
            check=False,
            timeout=900,
        )
        probe.recheck(expected_mode="0o444")
        python_target.recheck(expected_mode="0o775")
    finally:
        python_target.close()
        probe.close()
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(f"Recovery read-only probe failed: exit={result.returncode}")
    payload = json.loads(result.stdout)
    if (
        payload.get("pass") is not True
        or payload.get("no_output_writes") is not True
        or payload.get("prior_inputs_verified") is not True
        or payload.get("prior_recorded_verifier_recheck", {}).get("pass") is not True
        or payload.get("static_contract", {}).get("pass") is not True
        or payload.get("forbidden_output_slots_checked")
        != module.EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT
        or payload.get("regression_tests", {}).get("summary")
        != module.EXPECTED_REGRESSION_SUMMARY
        or payload.get("review_evidence_state") != "all_present"
        or not isinstance(payload.get("sources"), Mapping)
        or set(payload["sources"]) != set(module.RECOVERY_SOURCE_PATHS)
        or payload.get("self_execution_binding", {}).get("pass") is not True
        or payload.get("regression_tests", {}).get("source_execution_binding", {}).get(
            "pass"
        )
        is not True
    ):
        raise RuntimeError("Recovery read-only probe contract drift")
    execution = {
        "probe_source": payload["self_execution_binding"]["source"],
        "python_runtime": payload["self_execution_binding"]["python_runtime"],
        "launch": "bound_python_fd_exec_a_canonical_alias_bound_probe_source_fd",
        "regression_test_source": payload["regression_tests"][
            "source_execution_binding"
        ],
        "pass": True,
    }
    return payload, execution


def required_module_paths(module: ModuleType) -> set[Path]:
    return (
        (set(module.RECOVERY_SOURCE_PATHS.values()) - {module.VALIDATOR})
        | set(module.LEGACY_SOURCE_PATHS.values())
        | set(module.ORIGINAL_CHAIN_PATHS.values())
        | set(module.PRIOR_RECOVERY_CHAIN_PATHS.values())
        | set(module.REVIEW_PATHS.values())
        | {
            module.PUBLIC_KEY,
            module.PRIOR_RECOVERY_PRIVATE_KEY,
            module.REJECTED_V2_EVIDENCE,
            module.REJECTED_V2_PRIVATE_KEY,
            module.REJECTED_V3_EVIDENCE,
            module.REJECTED_V3_PRIVATE_KEY,
            module.REJECTED_V4_EVIDENCE,
            module.REJECTED_V4_PRIVATE_KEY,
            module.REJECTED_V5_EVIDENCE,
            module.REJECTED_V5_PRIVATE_KEY,
            module.REJECTED_V6_EVIDENCE,
            module.REJECTED_V6_PRIVATE_KEY,
            *module.REJECTED_V6_ARCHIVED_REVIEWS.values(),
            module.REJECTED_V7_EVIDENCE,
            module.REJECTED_V7_PRIVATE_KEY,
            module.REJECTED_V8_EVIDENCE,
            module.REJECTED_V8_PRIVATE_KEY,
            module.REJECTED_V19_ROUND1_EVIDENCE,
            *module.REJECTED_V19_ROUND1_REVIEW_PATHS.values(),
            module.REJECTED_V19_ROUND1_EXTERNAL_ENVELOPE,
            module.REJECTED_V19_ROUND1_EXTERNAL_STDERR,
            *module.REJECTED_V19_ROUND1_TRANSPORT_PATHS.values(),
            *module.REJECTED_V19_ROUND1_SOURCE_PATHS.values(),
            module.FAILED_V9_EVIDENCE,
            module.FAILED_V9_PUBLIC_KEY,
            module.FAILED_V9_PRIVATE_KEY,
            *module.FAILED_V9_ARCHIVED_ARTIFACT_PATHS.values(),
            module.FAILED_V10_EVIDENCE,
            module.FAILED_V10_PUBLIC_KEY,
            module.FAILED_V10_PRIVATE_KEY,
            *module.FAILED_V10_SOURCE_PATHS.values(),
            *module.FAILED_V10_ARCHIVED_ARTIFACT_PATHS.values(),
            module.FAILED_V11_EVIDENCE,
            module.FAILED_V11_PUBLIC_KEY,
            module.FAILED_V11_PRIVATE_KEY,
            module.FAILED_V11_CONTINUATION_PUBLIC_KEY,
            module.FAILED_V11_CONTINUATION_PRIVATE_KEY,
            *module.FAILED_V11_SOURCE_PATHS.values(),
            *module.FAILED_V11_ARCHIVED_ARTIFACT_PATHS.values(),
            module.FAILED_V12_EVIDENCE,
            module.FAILED_V12_PUBLIC_KEY,
            module.FAILED_V12_PRIVATE_KEY,
            module.FAILED_V12_CONTINUATION_PUBLIC_KEY,
            module.FAILED_V12_CONTINUATION_PRIVATE_KEY,
            *module.FAILED_V12_SOURCE_PATHS.values(),
            *module.FAILED_V12_ARCHIVED_ARTIFACT_PATHS.values(),
            module.REJECTED_V20_EVIDENCE,
            *module.REJECTED_V20_SOURCE_PATHS.values(),
            *module.REJECTED_V20_INTERNAL_REVIEW_PATHS.values(),
            module.REJECTED_V20_EXTERNAL_ENVELOPE,
            module.REJECTED_V20_EXTERNAL_STDERR,
            *module.REJECTED_V20_TRANSPORT_PATHS.values(),
            module.REJECTED_V20_PUBLIC_KEY,
            module.REJECTED_V20_PRIVATE_KEY,
            module.REJECTED_V20_CONTINUATION_PUBLIC_KEY,
            module.REJECTED_V20_CONTINUATION_PRIVATE_KEY,
            module.FAILED_V21_EVIDENCE,
            module.FAILED_V21_PUBLIC_KEY,
            module.FAILED_V21_PRIVATE_KEY,
            module.FAILED_V21_CONTINUATION_PUBLIC_KEY,
            module.FAILED_V21_CONTINUATION_PRIVATE_KEY,
            *module.FAILED_V21_SOURCE_PATHS.values(),
            *module.FAILED_V21_ARCHIVED_ARTIFACT_PATHS.values(),
            *module.FAILED_V21_REVIEW_TRANSPORT_PATHS.values(),
            module.REJECTED_V22_EVIDENCE,
            *module.REJECTED_V22_SOURCE_PATHS.values(),
            *module.REJECTED_V22_REVIEW_PATHS.values(),
            *module.REJECTED_V22_TRANSPORT_PATHS.values(),
            module.REJECTED_V22_PUBLIC_KEY,
            module.REJECTED_V22_PRIVATE_KEY,
            module.REJECTED_V22_CONTINUATION_PUBLIC_KEY,
            module.REJECTED_V22_CONTINUATION_PRIVATE_KEY,
            module.REJECTED_V23_EVIDENCE,
            *module.REJECTED_V23_SOURCE_PATHS.values(),
            *module.REJECTED_V23_REVIEW_PATHS.values(),
            *module.REJECTED_V23_TRANSPORT_PATHS.values(),
            module.REJECTED_V23_PUBLIC_KEY,
            module.REJECTED_V23_PRIVATE_KEY,
            module.REJECTED_V23_CONTINUATION_PUBLIC_KEY,
            module.REJECTED_V23_CONTINUATION_PRIVATE_KEY,
            module.REJECTED_V24_EVIDENCE,
            *module.REJECTED_V24_SOURCE_PATHS.values(),
            *module.REJECTED_V24_REVIEW_PATHS.values(),
            *module.REJECTED_V24_TRANSPORT_PATHS.values(),
            module.REJECTED_V24_PUBLIC_KEY,
            module.REJECTED_V24_PRIVATE_KEY,
            module.REJECTED_V24_CONTINUATION_PUBLIC_KEY,
            module.REJECTED_V24_CONTINUATION_PRIVATE_KEY,
            module.AUTHORIZED_CONTINUATION_PUBLIC_KEY,
            module.AUTHORIZED_CONTINUATION_PRIVATE_KEY,
            module.RECOVERY_CONTINUATION_PUBLIC_KEY,
            module.RECOVERY_CONTINUATION_PRIVATE_KEY,
        }
    )


def build_payloads(
    module: ModuleType,
    builder: BoundInput,
    validator: BoundInput,
    transaction_id: str,
    *,
    builder_binding: Mapping[str, Any],
    probe_execution: Mapping[str, Any],
    expected_probe_sources: Mapping[str, Any],
) -> tuple[dict[str, dict[str, Any]], bytes, int]:
    recovery_paths = {
        role: path
        for role, path in module.RECOVERY_SOURCE_PATHS.items()
        if role not in {"authority_validator", "ceremony_builder"}
    }
    module_builder_record = module._records(
        {"ceremony_builder": module.CEREMONY_BUILDER}
    )["ceremony_builder"]
    recovery_sources = {
        "authority_validator": validator.record(),
        "ceremony_builder": builder.record(),
        **module._records(recovery_paths),
    }
    if (
        not module._strict_equal(module_builder_record, builder.record())
        or not module._strict_equal(recovery_sources, expected_probe_sources)
        or not module._strict_equal(builder_binding.get("source"), builder.record())
    ):
        raise RuntimeError(
            "Probe source inventory differs from authority-bound recovery sources"
        )
    ceremony_execution_bindings = {
        "builder": dict(builder_binding),
        "readonly_probe_and_regression_tests": dict(probe_execution),
        "review_attribution_semantics": module.REVIEW_ATTRIBUTION_SEMANTICS,
        "pass": True,
    }
    module._validate_ceremony_execution_bindings(
        ceremony_execution_bindings, recovery_sources
    )
    legacy_sources = module._records(module.LEGACY_SOURCE_PATHS)
    original_chain = module._records(module.ORIGINAL_CHAIN_PATHS)
    module._validate_pinned_records(legacy_sources, original_chain)
    prior_records = module._records(module.PRIOR_RECOVERY_CHAIN_PATHS)
    prior_recovery_chain = module._validate_prior_recovery_chain(prior_records)
    rejected_generations = module._validate_rejected_generations()
    prior_failed_signed_recovery = module._validate_prior_failed_signed_recoveries()
    fresh_key_bootstrap = module._validate_v30_key_bootstrap()
    terminal_key_rotation, _ = module._validate_terminal_key_rotation(
        expected_recovery_private_mode="0o400"
    )

    public_data, public_key, public_fd = module._open_bound(module.PUBLIC_KEY)
    if (
        len(public_data) != 113
        or public_key["sha256"] != module.EXPECTED_PUBLIC_KEY_SHA256
    ):
        raise RuntimeError("Recovery public-key identity drift")

    review_records: dict[str, dict[str, Any]] = {}
    for role, path in module.REVIEW_PATHS.items():
        review_data, review_record, _ = module._open_bound(path)
        payload = module._load_json(review_data, f"pre-existing recovery review {role}")
        module._validate_review(
            role,
            payload,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
        review_records[role] = review_record

    authority = {
        "schema_name": module.SCHEMA_NAME,
        "schema_version": module.SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "authority_id": module.AUTHORITY_ID,
        "task_type": "B_comprehensive_validation",
        "status": module.AUTHORITY_STATUS,
        "authority_bundle": {
            "path": str(module.AUTHORITY_BUNDLE),
            "publication": "atomic_directory_rename_noreplace",
            "transaction_id": transaction_id,
            "members": MEMBERS,
        },
        "scope": module.RECOVERY_SCOPE,
        "terminal_key_rotation": terminal_key_rotation,
        "public_key": public_key,
        "retired_private_key": {"path": str(module.PRIVATE_KEY), "mode": "0o0"},
        "original_signed_chain": original_chain,
        "prior_recovery_chain": prior_recovery_chain,
        "prior_failed_signed_recovery": prior_failed_signed_recovery,
        "fresh_key_bootstrap": fresh_key_bootstrap,
        "rejected_unsigned_generations": rejected_generations,
        "legacy_sources": legacy_sources,
        "recovery_sources": recovery_sources,
        "legacy_commands": module.LEGACY_COMMANDS,
        "recovery_commands": module.RECOVERY_COMMANDS,
        "schemas": module.SCHEMAS,
        "receipt_paths": module.RECEIPT_PATHS,
        "review_artifacts": review_records,
        "review_attribution_semantics": module.REVIEW_ATTRIBUTION_SEMANTICS,
        "ceremony_execution_bindings": ceremony_execution_bindings,
        "checks": module.AUTHORITY_CHECKS,
        "pass": True,
        "pass_semantics": module.AUTHORITY_PASS_SEMANTICS,
    }
    module._validate_authority_payload(
        authority,
        transaction_id=transaction_id,
        public_record=public_key,
        expected_private={"path": str(module.PRIVATE_KEY), "mode": "0o0"},
        original_chain=original_chain,
        prior_recovery_chain=prior_recovery_chain,
        prior_failed_signed_recovery=prior_failed_signed_recovery,
        fresh_key_bootstrap=fresh_key_bootstrap,
        rejected_generations=rejected_generations,
        legacy_sources=legacy_sources,
        recovery_sources=recovery_sources,
        review_records=review_records,
        terminal_key_rotation=terminal_key_rotation,
    )
    return review_records, stable_json(authority), public_fd


def create_memfd(name: str, flags: int) -> int:
    if hasattr(os, "memfd_create"):
        return os.memfd_create(name, flags=flags)
    if _MEMFD_CREATE is None:
        raise RuntimeError("memfd_create is unavailable")
    descriptor = _MEMFD_CREATE(os.fsencode(name), flags)
    if descriptor < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), name)
    return descriptor


def sealed_memfd(name: str, data: bytes) -> int:
    flags = LINUX_MFD_CLOEXEC | LINUX_MFD_ALLOW_SEALING
    descriptor = create_memfd(name, flags)
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise RuntimeError(f"Short memfd write: {name}")
            offset += written
        os.fchmod(descriptor, 0o400)
        fcntl.fcntl(
            descriptor,
            LINUX_F_ADD_SEALS,
            LINUX_F_SEAL_SEAL
            | LINUX_F_SEAL_SHRINK
            | LINUX_F_SEAL_GROW
            | LINUX_F_SEAL_WRITE,
        )
        os.lseek(descriptor, 0, os.SEEK_SET)
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def stage_signing_key(private: BoundInput) -> int:
    private.recheck(expected_mode="0o400")
    descriptor = sealed_memfd("schema_recovery_signing_key_v30", private.data)
    try:
        staged = os.fstat(descriptor)
        expected_seals = (
            LINUX_F_SEAL_SEAL
            | LINUX_F_SEAL_SHRINK
            | LINUX_F_SEAL_GROW
            | LINUX_F_SEAL_WRITE
        )
        if (
            not stat.S_ISREG(staged.st_mode)
            or stat.S_IMODE(staged.st_mode) != 0o400
            or staged.st_nlink != 0
            or staged.st_size != len(private.data)
            or fcntl.fcntl(descriptor, LINUX_F_GET_SEALS) != expected_seals
            or read_fd(descriptor, staged.st_size) != private.data
        ):
            raise RuntimeError("Sealed recovery signing-key memfd drift")
        private.recheck(expected_mode="0o400")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def verify_signature_fds(
    module: ModuleType, data_fd: int, public_fd: int, signature_fd: int
) -> None:
    result = subprocess.run(
        [
            str(module.OPENSSL),
            "pkeyutl",
            "-verify",
            "-pubin",
            "-inkey",
            f"/proc/self/fd/{public_fd}",
            "-rawin",
            "-in",
            f"/proc/self/fd/{data_fd}",
            "-sigfile",
            f"/proc/self/fd/{signature_fd}",
        ],
        env={"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8"},
        pass_fds=(data_fd, public_fd, signature_fd),
        capture_output=True,
        check=False,
        timeout=60,
    )
    if result.returncode != 0:
        raise RuntimeError("Recovery authority signature verification failed")


def sign_and_verify(
    module: ModuleType,
    authority_data: bytes,
    signing_key_fd: int,
    public_fd: int,
) -> bytes:
    authority_fd = sealed_memfd("schema_recovery_authority_v30", authority_data)
    signature_fd = create_memfd(
        "schema_recovery_signature_v30",
        LINUX_MFD_CLOEXEC | LINUX_MFD_ALLOW_SEALING,
    )
    try:
        result = subprocess.run(
            [
                str(module.OPENSSL),
                "pkeyutl",
                "-sign",
                "-rawin",
                "-inkey",
                f"/proc/self/fd/{signing_key_fd}",
                "-in",
                f"/proc/self/fd/{authority_fd}",
                "-out",
                f"/proc/self/fd/{signature_fd}",
            ],
            env={"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8"},
            pass_fds=(signing_key_fd, authority_fd, signature_fd),
            capture_output=True,
            check=False,
            timeout=60,
        )
        if result.returncode != 0:
            raise RuntimeError("Recovery authority signing failed")
        signature = os.pread(signature_fd, 128, 0)
        if len(signature) != 64:
            raise RuntimeError("Recovery authority signature is not 64 bytes")
        fcntl.fcntl(
            signature_fd,
            LINUX_F_ADD_SEALS,
            LINUX_F_SEAL_SEAL
            | LINUX_F_SEAL_SHRINK
            | LINUX_F_SEAL_GROW
            | LINUX_F_SEAL_WRITE,
        )
        verify_signature_fds(module, authority_fd, public_fd, signature_fd)
        return signature
    finally:
        os.close(signature_fd)
        os.close(authority_fd)


def member_record(name: str, data: bytes) -> dict[str, Any]:
    return {
        "name": name,
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": "0o444",
    }


def published_record(path: Path, data: bytes) -> dict[str, Any]:
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": "0o444",
    }


def commit_data(
    module: ModuleType,
    transaction_id: str,
    authority_data: bytes,
    signature: bytes,
) -> bytes:
    return stable_json(
        {
            "schema_name": "intersubmod.tumor_ref_schema_recovery_atomic_commit",
            "schema_version": module.SCHEMA_VERSION,
            "transaction_id": transaction_id,
            "members": MEMBERS,
            "authority": member_record("authority.json", authority_data),
            "signature": member_record("authority.ed25519.sig", signature),
            "retired_private_key": {"path": str(module.PRIVATE_KEY), "mode": "0o0"},
            "pass": True,
        }
    )


def write_member(directory_fd: int, name: str, data: bytes) -> int:
    flags = os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(name, flags, 0o400, dir_fd=directory_fd)
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise RuntimeError(f"Short atomic-bundle member write: {name}")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
        observed = os.fstat(descriptor)
        if (
            observed.st_nlink != 1
            or stat.S_IMODE(observed.st_mode) != 0o444
            or read_fd(descriptor, observed.st_size) != data
        ):
            raise RuntimeError(f"Atomic-bundle member verification failed: {name}")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def stage_bundle(
    module: ModuleType,
    parent: DirectoryLease,
    transaction_id: str,
    authority_data: bytes,
    signature: bytes,
    public_fd: int,
) -> tuple[str, int, dict[str, int]]:
    stage_name = f".{module.AUTHORITY_BUNDLE.name}.staging.{transaction_id}"
    os.mkdir(stage_name, 0o700, dir_fd=parent.descriptor)
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    stage_fd = os.open(stage_name, flags, dir_fd=parent.descriptor)
    os.fchmod(stage_fd, 0o700)
    descriptors: dict[str, int] = {}
    try:
        descriptors["authority.json"] = write_member(
            stage_fd, "authority.json", authority_data
        )
        descriptors["authority.ed25519.sig"] = write_member(
            stage_fd, "authority.ed25519.sig", signature
        )
        commit = commit_data(module, transaction_id, authority_data, signature)
        descriptors["commit.json"] = write_member(stage_fd, "commit.json", commit)
        if set(os.listdir(stage_fd)) != set(MEMBERS):
            raise RuntimeError("Atomic-bundle staged member set drift")
        verify_signature_fds(
            module,
            descriptors["authority.json"],
            public_fd,
            descriptors["authority.ed25519.sig"],
        )
        parsed = module._validate_bundle_commit(
            commit,
            {
                "path": str(module.AUTHORITY),
                **member_record("authority.json", authority_data),
            },
            {
                "path": str(module.AUTHORITY_SIGNATURE),
                **member_record("authority.ed25519.sig", signature),
            },
            {"path": str(module.PRIVATE_KEY), "mode": "0o0"},
        )
        if parsed["transaction_id"] != transaction_id:
            raise RuntimeError("Atomic-bundle transaction ID drift")
        os.fsync(stage_fd)
        return stage_name, stage_fd, descriptors
    except BaseException:
        for descriptor in descriptors.values():
            os.close(descriptor)
        os.close(stage_fd)
        raise


def recheck_staged_bundle(
    module: ModuleType,
    parent: DirectoryLease,
    stage_name: str,
    stage_fd: int,
    member_fds: Mapping[str, int],
    transaction_id: str,
    authority_data: bytes,
    signature: bytes,
    public_fd: int,
) -> dict[str, Any]:
    terminal_stage_identity_recheck(
        module,
        parent,
        stage_name,
        stage_fd,
        member_fds,
        transaction_id,
        authority_data,
        signature,
    )
    expected_data = {
        "authority.json": authority_data,
        "authority.ed25519.sig": signature,
        "commit.json": commit_data(module, transaction_id, authority_data, signature),
    }
    verify_signature_fds(
        module,
        member_fds["authority.json"],
        public_fd,
        member_fds["authority.ed25519.sig"],
    )
    module._validate_bundle_commit(
        expected_data["commit.json"],
        {
            "path": str(module.AUTHORITY),
            **member_record("authority.json", authority_data),
        },
        {
            "path": str(module.AUTHORITY_SIGNATURE),
            **member_record("authority.ed25519.sig", signature),
        },
        {"path": str(module.PRIVATE_KEY), "mode": "0o0"},
    )
    return terminal_stage_identity_recheck(
        module,
        parent,
        stage_name,
        stage_fd,
        member_fds,
        transaction_id,
        authority_data,
        signature,
    )


def terminal_stage_identity_recheck(
    module: ModuleType,
    parent: DirectoryLease,
    stage_name: str,
    stage_fd: int,
    member_fds: Mapping[str, int],
    transaction_id: str,
    authority_data: bytes,
    signature: bytes,
    *,
    expected_tokens: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    parent.recheck()
    staged = os.fstat(stage_fd)
    live_stage = os.stat(
        stage_name, dir_fd=parent.descriptor, follow_symlinks=False
    )
    if (
        not stat.S_ISDIR(staged.st_mode)
        or stat.S_IMODE(staged.st_mode) != 0o700
        or staged.st_nlink != 2
        or full_identity(staged) != full_identity(live_stage)
        or set(os.listdir(stage_fd)) != set(MEMBERS)
        or set(member_fds) != set(MEMBERS)
    ):
        raise RuntimeError("Atomic-bundle staging directory identity drift")

    expected_data = {
        "authority.json": authority_data,
        "authority.ed25519.sig": signature,
        "commit.json": commit_data(module, transaction_id, authority_data, signature),
    }
    member_tokens: dict[str, tuple[int, ...]] = {}
    for name, expected in expected_data.items():
        descriptor = member_fds[name]
        opened = os.fstat(descriptor)
        live = os.stat(name, dir_fd=stage_fd, follow_symlinks=False)
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or stat.S_IMODE(opened.st_mode) != 0o444
            or full_identity(opened) != full_identity(live)
            or read_fd(descriptor, opened.st_size) != expected
        ):
            raise RuntimeError(f"Atomic-bundle staged member drift: {name}")
        member_tokens[name] = full_identity(opened)

    parent_opened = os.fstat(parent.descriptor)
    parent_live = os.stat(parent.path, follow_symlinks=False)
    staged_terminal = os.fstat(stage_fd)
    live_stage_terminal = os.stat(
        stage_name, dir_fd=parent.descriptor, follow_symlinks=False
    )
    if (
        full_identity(parent_opened) != full_identity(parent_live)
        or full_identity(staged_terminal) != full_identity(live_stage_terminal)
        or full_identity(staged_terminal) != full_identity(staged)
    ):
        raise RuntimeError("Atomic-bundle terminal directory generation drift")
    tokens: dict[str, Any] = {
        "parent": full_identity(parent_opened),
        "stage": full_identity(staged_terminal),
        "members": member_tokens,
    }
    if expected_tokens is not None and tokens != expected_tokens:
        raise RuntimeError("Atomic-bundle terminal identity tokens changed")
    return tokens


def rename_noreplace(parent_fd: int, source: str, destination: str) -> None:
    if _RENAMEAT2 is None:
        raise RuntimeError("renameat2(RENAME_NOREPLACE) is unavailable")
    result = _RENAMEAT2(
        parent_fd,
        os.fsencode(source),
        parent_fd,
        os.fsencode(destination),
        RENAME_NOREPLACE,
    )
    if result != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), destination)


def atomic_publish_after_terminal_contract(
    parent: DirectoryLease,
    absence_parents: Mapping[Path, DirectoryLease],
    mutation_watch: DirectoryMutationWatch,
    expected_absence_generations: Mapping[str, Any],
    stage_name: str,
    stage_fd: int,
    member_fds: Mapping[str, int],
    destination: str,
    expected_tokens: Mapping[str, Any],
    post_publish: Callable[[], None],
) -> None:
    if _RENAMEAT2 is None or not hasattr(signal, "pthread_sigmask"):
        raise RuntimeError("Terminal no-replace publication primitives are unavailable")
    source_bytes = os.fsencode(stage_name)
    destination_bytes = os.fsencode(destination)
    member_names = tuple(sorted(member_fds))
    expected_parent = tuple(expected_tokens["parent"])
    expected_stage = tuple(expected_tokens["stage"])
    expected_members = {
        name: tuple(expected_tokens["members"][name]) for name in member_names
    }
    expected_absence_paths = {str(path) for path in absence_parents}
    if (
        mutation_watch.descriptor < 0
        or mutation_watch.paths != frozenset(absence_parents)
        or set(expected_absence_generations) != expected_absence_paths
        or str(parent.path) not in expected_absence_generations
        or tuple(expected_absence_generations[str(parent.path)]) != expected_parent
    ):
        raise RuntimeError("Terminal absence-parent token contract drift")
    absence_parent_checks = tuple(
        (
            str(path),
            absence_parents[path].descriptor,
            path,
            tuple(expected_absence_generations[str(path)]),
        )
        for path in sorted(absence_parents, key=str)
    )
    blocked_signals = set(signal.valid_signals()) - {signal.SIGKILL, signal.SIGSTOP}
    previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, blocked_signals)
    try:
        if len(os.listdir("/proc/self/task")) != 1:
            raise RuntimeError("Terminal publication requires a single-threaded process")
        if set(os.listdir(stage_fd)) != set(MEMBERS) or set(member_fds) != set(MEMBERS):
            raise RuntimeError("Atomic-bundle terminal member-set drift")
        for name in member_names:
            opened = os.fstat(member_fds[name])
            live = os.stat(name, dir_fd=stage_fd, follow_symlinks=False)
            opened_identity = (
                opened.st_dev,
                opened.st_ino,
                opened.st_mode,
                opened.st_nlink,
                opened.st_size,
                opened.st_mtime_ns,
                opened.st_ctime_ns,
            )
            live_identity = (
                live.st_dev,
                live.st_ino,
                live.st_mode,
                live.st_nlink,
                live.st_size,
                live.st_mtime_ns,
                live.st_ctime_ns,
            )
            if opened_identity != expected_members[name] or live_identity != expected_members[name]:
                raise RuntimeError(f"Atomic-bundle terminal member token drift: {name}")
        staged = os.fstat(stage_fd)
        live_stage = os.stat(stage_name, dir_fd=parent.descriptor, follow_symlinks=False)
        staged_identity = (
            staged.st_dev,
            staged.st_ino,
            staged.st_mode,
            staged.st_nlink,
            staged.st_size,
            staged.st_mtime_ns,
            staged.st_ctime_ns,
        )
        live_stage_identity = (
            live_stage.st_dev,
            live_stage.st_ino,
            live_stage.st_mode,
            live_stage.st_nlink,
            live_stage.st_size,
            live_stage.st_mtime_ns,
            live_stage.st_ctime_ns,
        )
        if staged_identity != expected_stage or live_stage_identity != expected_stage:
            raise RuntimeError("Atomic-bundle terminal stage token drift")
        for path_text, descriptor, path, expected_generation in absence_parent_checks:
            opened_parent = os.fstat(descriptor)
            live_parent = os.stat(path, follow_symlinks=False)
            opened_parent_identity = (
                opened_parent.st_dev,
                opened_parent.st_ino,
                opened_parent.st_mode,
                opened_parent.st_nlink,
                opened_parent.st_size,
                opened_parent.st_mtime_ns,
                opened_parent.st_ctime_ns,
            )
            live_parent_identity = (
                live_parent.st_dev,
                live_parent.st_ino,
                live_parent.st_mode,
                live_parent.st_nlink,
                live_parent.st_size,
                live_parent.st_mtime_ns,
                live_parent.st_ctime_ns,
            )
            if (
                opened_parent_identity != expected_generation
                or live_parent_identity != expected_generation
            ):
                raise RuntimeError(
                    f"Atomic-bundle terminal absence-parent token drift: {path_text}"
                )
        mutation_watch.assert_quiet()
        result = _RENAMEAT2(
            parent.descriptor,
            source_bytes,
            parent.descriptor,
            destination_bytes,
            RENAME_NOREPLACE,
        )
        if result != 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error), destination)
        post_publish()
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)


def publish_staged_bundle(
    module: ModuleType,
    parent: DirectoryLease,
    absence_parents: Mapping[Path, DirectoryLease],
    stage_name: str,
    stage_fd: int,
    member_fds: Mapping[str, int],
    transaction_id: str,
    authority_data: bytes,
    signature: bytes,
    public_fd: int,
) -> None:
    require_final_slots_absent(module)
    post_external_tokens = recheck_staged_bundle(
        module,
        parent,
        stage_name,
        stage_fd,
        member_fds,
        transaction_id,
        authority_data,
        signature,
        public_fd,
    )
    mutation_watch = DirectoryMutationWatch(absence_parents)
    try:
        final_absence = require_ceremony_absence(
            module,
            absence_parents,
            allowed_stage_name=stage_name,
            mutation_watch=mutation_watch,
        )
        terminal_tokens = terminal_stage_identity_recheck(
            module,
            parent,
            stage_name,
            stage_fd,
            member_fds,
            transaction_id,
            authority_data,
            signature,
            expected_tokens=post_external_tokens,
        )

        def post_publish() -> None:
            os.fsync(parent.descriptor)
            final = os.stat(module.AUTHORITY_BUNDLE, follow_symlinks=False)
            staged = os.fstat(stage_fd)
            if (
                (
                    final.st_dev,
                    final.st_ino,
                    stat.S_IMODE(final.st_mode),
                    final.st_nlink,
                )
                != (staged.st_dev, staged.st_ino, 0o700, 2)
                or set(os.listdir(stage_fd)) != set(MEMBERS)
            ):
                raise RuntimeError("Atomic-bundle final identity drift")
            terminal_stage_identity_recheck(
                module,
                parent,
                module.AUTHORITY_BUNDLE.name,
                stage_fd,
                member_fds,
                transaction_id,
                authority_data,
                signature,
            )

        atomic_publish_after_terminal_contract(
            parent,
            absence_parents,
            mutation_watch,
            final_absence["directory_generation_tokens"],
            stage_name,
            stage_fd,
            member_fds,
            module.AUTHORITY_BUNDLE.name,
            terminal_tokens,
            post_publish,
        )
    finally:
        mutation_watch.close()


def require_rename_noreplace_capability() -> None:
    if _RENAMEAT2 is None or not hasattr(signal, "pthread_sigmask"):
        raise RuntimeError("renameat2(RENAME_NOREPLACE) is unavailable")


def terminal_recheck(
    module: ModuleType,
    builder: BoundInput,
    validator: BoundInput,
    private: BoundInput,
    writer: WriterLease,
    parent: DirectoryLease,
    absence_parents: Mapping[Path, DirectoryLease],
    *,
    private_mode: str,
    allowed_stage_name: str | None = None,
) -> None:
    builder.recheck(expected_mode="0o444")
    validator.recheck(expected_mode="0o444")
    writer.recheck()
    module._require_leases(required_paths=required_module_paths(module))
    private.recheck(
        expected_mode=private_mode,
        allow_retirement=private_mode == "0o0",
    )
    parent.recheck()
    require_ceremony_absence(
        module, absence_parents, allowed_stage_name=allowed_stage_name
    )


def retire_private_key(private: BoundInput) -> None:
    os.fchmod(private.descriptor, 0o000)
    os.fsync(private.descriptor)
    private.recheck(expected_mode="0o0", allow_retirement=True)


def preflight(
    module: ModuleType,
    builder: BoundInput,
    validator: BoundInput,
    builder_binding: Mapping[str, Any],
) -> dict[str, Any]:
    require_rename_noreplace_capability()
    writer = WriterLease.acquire(module)
    absence_parents: dict[Path, DirectoryLease] = {}
    private: BoundInput | None = None
    try:
        require_final_slots_absent(module)
        require_no_staging(module)
        absence_parents = open_absence_directory_leases(module)
        parent = absence_parents[module.RESULT_ROOT]
        private = BoundInput.open(module.PRIVATE_KEY, expected_mode="0o400")
        require_ceremony_absence(module, absence_parents)
        probe, probe_execution = run_probe(module)
        transaction_id = secrets.token_hex(16)
        reviews, authority_data, public_fd = build_payloads(
            module,
            builder,
            validator,
            transaction_id,
            builder_binding=builder_binding,
            probe_execution=probe_execution,
            expected_probe_sources=probe["sources"],
        )
        terminal_recheck(
            module,
            builder,
            validator,
            private,
            writer,
            parent,
            absence_parents,
            private_mode="0o400",
        )
        require_final_slots_absent(module)
        return {
            "reviews": reviews,
            "authority_data": authority_data,
            "probe": probe,
            "transaction_id": transaction_id,
            "public_fd": public_fd,
            "private": private,
            "writer": writer,
            "parent": parent,
            "absence_parents": absence_parents,
        }
    except BaseException:
        if private is not None:
            private.close()
        for lease in absence_parents.values():
            lease.close()
        writer.close()
        raise


def close_state(state: Mapping[str, Any]) -> None:
    state["private"].close()
    for lease in state["absence_parents"].values():
        lease.close()
    state["writer"].close()


def publish(
    module: ModuleType,
    builder: BoundInput,
    validator: BoundInput,
    builder_binding: Mapping[str, Any],
) -> dict[str, Any]:
    state = preflight(module, builder, validator, builder_binding)
    private: BoundInput = state["private"]
    writer: WriterLease = state["writer"]
    parent: DirectoryLease = state["parent"]
    absence_parents: Mapping[Path, DirectoryLease] = state["absence_parents"]
    retirement_required = False
    signing_key_fd = -1
    stage_fd = -1
    member_fds: dict[str, int] = {}
    try:
        terminal_recheck(
            module,
            builder,
            validator,
            private,
            writer,
            parent,
            absence_parents,
            private_mode="0o400",
        )
        signing_key_fd = stage_signing_key(private)
        terminal_recheck(
            module,
            builder,
            validator,
            private,
            writer,
            parent,
            absence_parents,
            private_mode="0o400",
        )
        retirement_required = True
        retire_private_key(private)
        terminal_recheck(
            module,
            builder,
            validator,
            private,
            writer,
            parent,
            absence_parents,
            private_mode="0o0",
        )
        signature = sign_and_verify(
            module,
            state["authority_data"],
            signing_key_fd,
            state["public_fd"],
        )
        os.close(signing_key_fd)
        signing_key_fd = -1
        terminal_recheck(
            module,
            builder,
            validator,
            private,
            writer,
            parent,
            absence_parents,
            private_mode="0o0",
        )
        stage_name, stage_fd, member_fds = stage_bundle(
            module,
            parent,
            state["transaction_id"],
            state["authority_data"],
            signature,
            state["public_fd"],
        )
        terminal_recheck(
            module,
            builder,
            validator,
            private,
            writer,
            parent,
            absence_parents,
            private_mode="0o0",
            allowed_stage_name=stage_name,
        )
        publish_staged_bundle(
            module,
            parent,
            absence_parents,
            stage_name,
            stage_fd,
            member_fds,
            state["transaction_id"],
            state["authority_data"],
            signature,
            state["public_fd"],
        )
        commit = commit_data(
            module, state["transaction_id"], state["authority_data"], signature
        )
        return {
            "schema_name": "intersubmod.tumor_ref_schema_recovery_authority_build",
            "schema_version": module.SCHEMA_VERSION,
            "authority_bundle": str(module.AUTHORITY_BUNDLE),
            "transaction_id": state["transaction_id"],
            "authority": published_record(module.AUTHORITY, state["authority_data"]),
            "signature": published_record(module.AUTHORITY_SIGNATURE, signature),
            "commit": published_record(module.AUTHORITY_COMMIT, commit),
            "review_count": len(state["reviews"]),
            "private_key_mode": "0o0",
            "probe_pass": state["probe"]["pass"],
            "publication_validation": {
                "basis": "retained_stage_directory_and_member_fds",
                "validator_path_reopened": False,
                "terminal_stage_recheck_after_external_verification": True,
                "pass": True,
            },
            "pass": True,
        }
    except BaseException:
        if retirement_required:
            try:
                retire_private_key(private)
            except BaseException as retirement_error:
                raise RuntimeError(
                    "Signed v30 ceremony failed and private-key retirement also failed"
                ) from retirement_error
        raise
    finally:
        if signing_key_fd >= 0:
            os.close(signing_key_fd)
        for descriptor in member_fds.values():
            os.close(descriptor)
        if stage_fd >= 0:
            os.close(stage_fd)
        close_state(state)


def main() -> int:
    builder, builder_binding = bind_executing_builder()
    validator: BoundInput | None = None
    try:
        module, validator = load_validator()
        if sys.argv[1:] == ["--preflight"]:
            state = preflight(module, builder, validator, builder_binding)
            try:
                result = {
                    "mode": "preflight",
                    "review_count": len(state["reviews"]),
                    "predicted_authority_size_bytes": len(state["authority_data"]),
                    "predicted_authority_sha256": hashlib.sha256(
                        state["authority_data"]
                    ).hexdigest(),
                    "probe_pass": state["probe"]["pass"],
                    "no_output_writes": state["probe"]["no_output_writes"],
                    "private_key_mode": "0o400",
                    "pass": True,
                }
            finally:
                close_state(state)
        elif sys.argv[1:] == ["--publish-after-three-approvals"]:
            result = publish(module, builder, validator, builder_binding)
        else:
            raise RuntimeError(
                "Usage: builder --preflight | --publish-after-three-approvals"
            )
        print(json.dumps(result, ensure_ascii=True, indent=2, sort_keys=True))
        return 0
    finally:
        if validator is not None:
            validator.close()
        builder.close()


if __name__ == "__main__":
    raise SystemExit(main())
