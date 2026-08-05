#!/usr/bin/env python3
"""Fail-closed lifecycle and frozen-provenance primitives for layered v3 runs.

This module is deliberately independent from the live v2 runner.  It provides a
small integration surface for a future v3 launcher:

* a per-run lock file held with ``flock`` in the run parent;
* a private staging directory which is published only after preflight artifacts
  have been frozen and cross-bound by a launch receipt;
* an atomic, hash-chained state machine and heartbeat;
* a source bundle that stores executable bytes, not only working-tree paths;
* an environment lock with a strict ``SM_*`` allowlist;
* tracked process groups with fail-fast sibling cancellation; and
* a final ``_SUCCESS`` marker that is created only after verification.

The module uses only the Python standard library so it can be imported by a
bootstrap runner before project-specific dependencies are trusted.
"""

from __future__ import annotations

import argparse
import atexit
import datetime as dt
import fcntl
import hashlib
import importlib.metadata
import json
import locale
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import socket
import stat
import subprocess
import sys
import threading
import time
import uuid
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union


PathLike = Union[str, os.PathLike]
ZERO_DIGEST = "0" * 64
RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
SAFE_LABEL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$")


class ContractError(RuntimeError):
    """Base exception carrying a stable fail-closed error code."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


class RunLockError(ContractError):
    pass


class StateTransitionError(ContractError):
    pass


class SourceBundleError(ContractError):
    pass


class EnvironmentContractError(ContractError):
    pass


class HeartbeatError(ContractError):
    pass


class ChildProcessFailure(ContractError):
    pass


class LifecycleSignal(RuntimeError):
    """Raised by signal handlers so the context manager can run one trap path."""

    def __init__(self, signum: int):
        self.signum = signum
        try:
            self.signal_name = signal.Signals(signum).name
        except ValueError:
            self.signal_name = str(signum)
        super().__init__(f"received signal {self.signal_name}")


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def canonical_json_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    fd = os.open(str(path), flags)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _regular_file_bytes(path: PathLike, reject_symlink: bool = False) -> Tuple[bytes, os.stat_result, Path]:
    requested = Path(path)
    try:
        lst = requested.lstat()
    except FileNotFoundError as exc:
        raise ContractError("E_FILE_MISSING", f"file does not exist: {requested}") from exc
    if reject_symlink and stat.S_ISLNK(lst.st_mode):
        raise ContractError("E_SYMLINK_FORBIDDEN", f"symlink is not accepted: {requested}")
    resolved = requested.resolve(strict=True)
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        fd = os.open(str(resolved), flags)
    except OSError as exc:
        raise ContractError("E_FILE_OPEN", f"cannot open regular file: {resolved}: {exc}") from exc
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode):
            raise ContractError("E_NOT_REGULAR_FILE", f"not a regular file: {resolved}")
        chunks: List[bytes] = []
        while True:
            chunk = os.read(fd, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after:
        raise ContractError("E_FILE_CHANGED_DURING_READ", f"file changed while being read: {resolved}")
    payload = b"".join(chunks)
    if len(payload) != before.st_size:
        raise ContractError("E_FILE_SHORT_READ", f"expected {before.st_size} bytes, read {len(payload)}: {resolved}")
    return payload, before, resolved


def sha256_file(path: PathLike, reject_symlink: bool = False) -> str:
    payload, _, _ = _regular_file_bytes(path, reject_symlink=reject_symlink)
    return sha256_bytes(payload)


def atomic_write_bytes(
    path: PathLike,
    payload: bytes,
    mode: int = 0o600,
    no_overwrite: bool = False,
) -> str:
    """Write bytes using same-directory temp, fsync, replace, and directory fsync."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if no_overwrite and target.exists():
        raise ContractError("E_ATOMIC_TARGET_EXISTS", f"refusing to overwrite: {target}")
    temp = target.parent / f".{target.name}.partial.{os.getpid()}.{uuid.uuid4().hex}"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    fd = os.open(str(temp), flags, mode)
    try:
        with os.fdopen(fd, "wb", closefd=True) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
            os.fchmod(handle.fileno(), mode)
        if no_overwrite and target.exists():
            raise ContractError("E_ATOMIC_TARGET_EXISTS", f"refusing to overwrite: {target}")
        os.replace(str(temp), str(target))
        _fsync_directory(target.parent)
    except BaseException:
        # A partial file is intentionally retained for failure forensics.
        raise
    return sha256_bytes(payload)


def atomic_write_json(
    path: PathLike,
    payload: Mapping[str, Any],
    mode: int = 0o600,
    no_overwrite: bool = False,
) -> str:
    encoded = canonical_json_bytes(payload) + b"\n"
    # Validate the exact bytes before they acquire the formal filename.
    parsed = json.loads(encoded.decode("utf-8"))
    if parsed != payload:
        raise ContractError("E_JSON_ROUNDTRIP", f"JSON payload did not round-trip for {path}")
    return atomic_write_bytes(path, encoded, mode=mode, no_overwrite=no_overwrite)


def read_json_with_sha256(path: PathLike) -> Tuple[Dict[str, Any], str]:
    payload, _, resolved = _regular_file_bytes(path, reject_symlink=True)
    try:
        document = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("E_JSON_INVALID", f"invalid JSON: {resolved}: {exc}") from exc
    if not isinstance(document, dict):
        raise ContractError("E_JSON_NOT_OBJECT", f"expected a JSON object: {resolved}")
    return document, sha256_bytes(payload)


class ParentRunLock:
    """A non-blocking per-run flock stored in the run parent directory."""

    def __init__(self, run_parent: PathLike, run_id: str):
        validate_run_id(run_id)
        self.run_parent = Path(run_parent).resolve()
        self.run_parent.mkdir(parents=True, exist_ok=True)
        if not self.run_parent.is_dir():
            raise RunLockError("E_RUN_PARENT", f"run parent is not a directory: {self.run_parent}")
        self.path = self.run_parent / f".layered_v3.{run_id}.launch.lock"
        self._fd: Optional[int] = None

    def acquire(self) -> None:
        if self._fd is not None:
            return
        fd = os.open(str(self.path), os.O_RDWR | os.O_CREAT, 0o600)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            os.close(fd)
            raise RunLockError("E_RUN_LOCKED", f"another launcher holds {self.path}") from exc
        os.ftruncate(fd, 0)
        os.write(fd, f"pid={os.getpid()} acquired_at_utc={utc_now()}\n".encode("utf-8"))
        os.fsync(fd)
        self._fd = fd

    def release(self) -> None:
        if self._fd is None:
            return
        try:
            fcntl.flock(self._fd, fcntl.LOCK_UN)
        finally:
            os.close(self._fd)
            self._fd = None

    @property
    def held(self) -> bool:
        return self._fd is not None


def validate_run_id(run_id: str) -> None:
    if not RUN_ID_RE.fullmatch(run_id) or run_id in {".", ".."}:
        raise ContractError("E_RUN_ID_INVALID", f"unsafe run ID: {run_id!r}")


def _safe_label(value: str, field: str) -> None:
    if not SAFE_LABEL_RE.fullmatch(value):
        raise ContractError("E_LABEL_INVALID", f"unsafe {field}: {value!r}")


def _sanitize_detail(value: str, limit: int = 4096) -> str:
    """Keep trap/state details single-line and bounded without hiding delimiters."""

    sanitized = value.replace("\r", "\\r").replace("\n", "\\n").replace("\t", "\\t")
    return sanitized if len(sanitized) <= limit else sanitized[: limit - 3] + "..."


@dataclass(frozen=True)
class SourceBundleResult:
    root: Path
    manifest_path: Path
    manifest_sha256: str
    content_sha256: str
    role_paths: Mapping[str, Path]


def _source_content_digest(entries: Sequence[Mapping[str, Any]]) -> str:
    identity = [
        {
            "role": item["role"],
            "bundled_path": item["bundled_path"],
            "mode": item["mode"],
            "size": item["size"],
            "sha256": item["sha256"],
        }
        for item in entries
    ]
    return sha256_bytes(canonical_json_bytes(identity))


def build_source_bundle(
    bundle_root: PathLike,
    runner: PathLike,
    validator: PathLike,
    verifier: PathLike,
    imported_scripts: Sequence[PathLike],
) -> SourceBundleResult:
    """Copy required executable source bytes and write a deterministic manifest."""

    root = Path(bundle_root)
    if root.exists():
        raise SourceBundleError("E_SOURCE_BUNDLE_EXISTS", f"bundle root already exists: {root}")
    if not imported_scripts:
        raise SourceBundleError("E_SOURCE_BUNDLE_IMPORTS_EMPTY", "at least one imported local script is required")
    root.mkdir(parents=True, mode=0o700)

    specifications: List[Tuple[str, PathLike, str]] = [
        ("runner", runner, "files/core/runner.py"),
        ("validator", validator, "files/core/validator.py"),
        ("verifier", verifier, "files/core/verifier.py"),
    ]
    for index, imported in enumerate(imported_scripts):
        imported_name = Path(imported).name
        safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", imported_name)
        specifications.append(
            (f"imported:{index:03d}:{safe_name}", imported, f"files/imported/{index:03d}_{safe_name}")
        )

    entries: List[Dict[str, Any]] = []
    role_paths: Dict[str, Path] = {}
    observed_sources: Set[Tuple[int, int]] = set()
    for role, source, relative in specifications:
        payload, source_stat, resolved = _regular_file_bytes(source, reject_symlink=True)
        source_identity = (source_stat.st_dev, source_stat.st_ino)
        if source_identity in observed_sources:
            raise SourceBundleError(
                "E_SOURCE_BUNDLE_DUPLICATE_SOURCE",
                f"source appears in more than one role: {resolved}",
            )
        observed_sources.add(source_identity)
        mode = stat.S_IMODE(source_stat.st_mode)
        bundled = root / relative
        atomic_write_bytes(bundled, payload, mode=mode)
        copied, copied_stat, _ = _regular_file_bytes(bundled, reject_symlink=True)
        digest = sha256_bytes(payload)
        if copied != payload or sha256_bytes(copied) != digest:
            raise SourceBundleError("E_SOURCE_BUNDLE_COPY_MISMATCH", f"copy readback mismatch: {bundled}")
        entry = {
            "role": role,
            "source_path": str(resolved),
            "bundled_path": relative,
            "mode": mode,
            "size": len(payload),
            "sha256": digest,
            "source_device": source_stat.st_dev,
            "source_inode": source_stat.st_ino,
            "source_mtime_ns": source_stat.st_mtime_ns,
            "bundled_device": copied_stat.st_dev,
            "bundled_inode": copied_stat.st_ino,
        }
        entries.append(entry)
        role_paths[role] = bundled

    content_digest = _source_content_digest(entries)
    manifest = {
        "schema_name": "intersubmod.layered_source_bundle",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "file_count": len(entries),
        "bundle_content_sha256": content_digest,
        "files": entries,
    }
    manifest_path = root / "source_bundle_manifest.json"
    manifest_sha = atomic_write_json(manifest_path, manifest)
    result = SourceBundleResult(root, manifest_path, manifest_sha, content_digest, role_paths)
    verify_source_bundle(
        manifest_path,
        expected_manifest_sha256=manifest_sha,
        expected_content_sha256=content_digest,
    )
    return result


def verify_source_bundle(
    manifest_path: PathLike,
    expected_manifest_sha256: Optional[str] = None,
    expected_content_sha256: Optional[str] = None,
) -> SourceBundleResult:
    manifest_file = Path(manifest_path)
    try:
        manifest, manifest_sha = read_json_with_sha256(manifest_file)
    except ContractError as exc:
        raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", str(exc)) from exc
    if expected_manifest_sha256 and manifest_sha != expected_manifest_sha256:
        raise SourceBundleError(
            "E_SOURCE_BUNDLE_MISMATCH",
            f"manifest digest mismatch: expected {expected_manifest_sha256}, observed {manifest_sha}",
        )
    required_top = {
        "schema_name",
        "schema_version",
        "created_at_utc",
        "file_count",
        "bundle_content_sha256",
        "files",
    }
    if set(manifest) != required_top:
        raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", "source bundle manifest has missing or unknown fields")
    if manifest["schema_name"] != "intersubmod.layered_source_bundle" or manifest["schema_version"] != "1.0.0":
        raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", "unsupported source bundle schema")
    files = manifest["files"]
    if not isinstance(files, list) or manifest["file_count"] != len(files):
        raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", "file_count does not match files array")
    roles = [item.get("role") for item in files if isinstance(item, dict)]
    if not {"runner", "validator", "verifier"}.issubset(set(roles)) or not any(
        isinstance(role, str) and role.startswith("imported:") for role in roles
    ):
        raise SourceBundleError(
            "E_SOURCE_BUNDLE_MISMATCH",
            "required runner/validator/verifier/imported roles are absent",
        )
    if len(roles) != len(set(roles)):
        raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", "duplicate source bundle role")

    root = manifest_file.parent.resolve()
    role_paths: Dict[str, Path] = {}
    expected_paths: Set[Path] = {manifest_file.resolve()}
    for item in files:
        required_item = {
            "role",
            "source_path",
            "bundled_path",
            "mode",
            "size",
            "sha256",
            "source_device",
            "source_inode",
            "source_mtime_ns",
            "bundled_device",
            "bundled_inode",
        }
        if not isinstance(item, dict) or set(item) != required_item:
            raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", "source bundle file entry is malformed")
        relative = Path(item["bundled_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", f"unsafe bundled path: {relative}")
        bundled = (root / relative).resolve()
        if bundled == root or root not in bundled.parents:
            raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", f"bundled path escapes root: {relative}")
        try:
            payload, observed_stat, _ = _regular_file_bytes(bundled, reject_symlink=True)
        except ContractError as exc:
            raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", str(exc)) from exc
        observed = (len(payload), sha256_bytes(payload), stat.S_IMODE(observed_stat.st_mode))
        expected = (item["size"], item["sha256"], item["mode"])
        if observed != expected:
            raise SourceBundleError(
                "E_SOURCE_BUNDLE_MISMATCH",
                f"bundled bytes/mode mismatch for {relative}: expected={expected}, observed={observed}",
            )
        expected_paths.add(bundled)
        role_paths[item["role"]] = bundled

    observed_paths = {path.resolve() for path in root.rglob("*") if path.is_file()}
    if observed_paths != expected_paths:
        unexpected = sorted(str(path.relative_to(root)) for path in observed_paths - expected_paths)
        missing = sorted(str(path.relative_to(root)) for path in expected_paths - observed_paths)
        raise SourceBundleError(
            "E_SOURCE_BUNDLE_MISMATCH",
            f"unmanifested or missing bundle files: unexpected={unexpected}, missing={missing}",
        )
    content_digest = _source_content_digest(files)
    if content_digest != manifest["bundle_content_sha256"]:
        raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", "bundle content digest does not match manifest")
    if expected_content_sha256 and content_digest != expected_content_sha256:
        raise SourceBundleError(
            "E_SOURCE_BUNDLE_MISMATCH",
            f"bundle content digest mismatch: expected {expected_content_sha256}, observed {content_digest}",
        )
    return SourceBundleResult(root, manifest_file, manifest_sha, content_digest, role_paths)


@dataclass(frozen=True)
class ToolSpec:
    name: str
    executable: str
    version_args: Tuple[str, ...] = ("--version",)


def capture_environment_lock(
    output_path: PathLike,
    sm_allowlist: Iterable[str],
    sm_defaults: Optional[Mapping[str, Optional[str]]] = None,
    required_sm: Iterable[str] = (),
    environment: Optional[Mapping[str, str]] = None,
    distributions: Sequence[str] = (),
    tools: Sequence[ToolSpec] = (),
    storage_path: Optional[PathLike] = None,
    determinism_required: Optional[Mapping[str, str]] = None,
) -> str:
    """Capture an environment lock without serializing unrelated secrets."""

    env = dict(os.environ if environment is None else environment)
    allowlist = set(sm_allowlist)
    defaults = dict(sm_defaults or {})
    required = set(required_sm)
    for name in allowlist | set(defaults) | required:
        if not re.fullmatch(r"SM_[A-Z0-9_]+", name):
            raise EnvironmentContractError("E_ENVIRONMENT_MISMATCH", f"invalid SM allowlist name: {name}")
    if not set(defaults).issubset(allowlist) or not required.issubset(allowlist):
        raise EnvironmentContractError("E_ENVIRONMENT_MISMATCH", "defaults/required variables must be allowlisted")
    unknown = sorted(name for name in env if name.startswith("SM_") and name not in allowlist)
    if unknown:
        raise EnvironmentContractError("E_UNKNOWN_SM_ENV", f"unknown SM_* variables: {unknown}")

    sm_values: Dict[str, Dict[str, Optional[str]]] = {}
    for name in sorted(allowlist):
        if name in env:
            value: Optional[str] = env[name]
            origin = "environment"
        elif name in defaults:
            value = defaults[name]
            origin = "default"
        else:
            value = None
            origin = "unset"
        if name in required and (value is None or value == ""):
            raise EnvironmentContractError("E_ENVIRONMENT_MISMATCH", f"required environment variable is unset: {name}")
        sm_values[name] = {"value": value, "origin": origin}

    determinism = dict(
        {"LC_ALL": "C", "TZ": "UTC", "PYTHONHASHSEED": "0"}
        if determinism_required is None
        else determinism_required
    )
    effective_determinism: Dict[str, str] = {}
    for name, expected in sorted(determinism.items()):
        observed = env.get(name)
        if observed != expected:
            raise EnvironmentContractError(
                "E_ENVIRONMENT_MISMATCH",
                f"determinism variable {name} expected {expected!r}, observed {observed!r}",
            )
        effective_determinism[name] = observed

    python_path = Path(sys.executable).resolve(strict=True)
    python_record = {
        "requested_path": sys.executable,
        "realpath": str(python_path),
        "sha256": sha256_file(python_path),
        "version": platform.python_version(),
        "implementation": platform.python_implementation(),
    }
    distribution_records: List[Dict[str, str]] = []
    for name in sorted(set(distributions)):
        try:
            version = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError as exc:
            raise EnvironmentContractError("E_ENVIRONMENT_MISMATCH", f"Python distribution is missing: {name}") from exc
        distribution_records.append({"name": name, "version": version})

    tool_records: List[Dict[str, Any]] = []
    tool_names: Set[str] = set()
    for spec in tools:
        _safe_label(spec.name, "tool name")
        if spec.name in tool_names:
            raise EnvironmentContractError("E_ENVIRONMENT_MISMATCH", f"duplicate tool name: {spec.name}")
        tool_names.add(spec.name)
        requested = spec.executable
        found = requested if os.path.isabs(requested) else shutil.which(requested, path=env.get("PATH"))
        if not found:
            raise EnvironmentContractError("E_ENVIRONMENT_MISMATCH", f"tool not found: {requested}")
        realpath = Path(found).resolve(strict=True)
        try:
            completed = subprocess.run(
                [str(realpath), *spec.version_args],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
                timeout=10,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise EnvironmentContractError(
                "E_ENVIRONMENT_MISMATCH",
                f"tool version probe failed: {spec.name}: {exc}",
            ) from exc
        if completed.returncode != 0:
            raise EnvironmentContractError(
                "E_ENVIRONMENT_MISMATCH",
                f"tool version probe returned {completed.returncode}: {spec.name}",
            )
        tool_records.append(
            {
                "name": spec.name,
                "requested": requested,
                "realpath": str(realpath),
                "sha256": sha256_file(realpath),
                "version_command": [str(realpath), *spec.version_args],
                "version_output": completed.stdout.strip(),
            }
        )

    storage_root = Path(storage_path if storage_path is not None else Path(output_path).parent).resolve()
    storage_stat = storage_root.stat()
    lock = {
        "schema_name": "intersubmod.layered_environment_lock",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "python": python_record,
        "python_distributions": distribution_records,
        "tools": tool_records,
        "system": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "kernel": platform.release(),
            "glibc": list(platform.libc_ver()),
            "locale": locale.setlocale(locale.LC_ALL),
            "timezone_env": env.get("TZ"),
            "timezone_names": list(time.tzname),
            "storage_path": str(storage_root),
            "storage_device": storage_stat.st_dev,
        },
        "determinism": effective_determinism,
        "sm_variables": sm_values,
    }
    return atomic_write_json(output_path, lock)


def _pid_start_time(pid: int) -> str:
    try:
        raw = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        closing = raw.rfind(")")
        fields = raw[closing + 2 :].split()
        # fields[0] is proc stat field 3 (state); field 22 is therefore index 19.
        return fields[19]
    except (OSError, IndexError) as exc:
        raise ContractError("E_PID_IDENTITY", f"cannot read Linux process start time for PID {pid}") from exc


@dataclass
class TrackedChild:
    label: str
    argv: Tuple[str, ...]
    command_sha256: str
    process: subprocess.Popen
    pid_start_time: str
    launched_at_utc: str
    returncode: Optional[int] = None


class TrackedChildRegistry:
    """Launch each child in a new process group and retain exact PID identity."""

    def __init__(self, snapshot_path: Callable[[], Path]):
        self._snapshot_path = snapshot_path
        self._children: Dict[int, TrackedChild] = {}
        self._lock = threading.RLock()

    def launch(self, argv: Sequence[str], label: str, **popen_kwargs: Any) -> TrackedChild:
        _safe_label(label, "child label")
        if not argv or any(not isinstance(arg, str) or "\x00" in arg for arg in argv):
            raise ChildProcessFailure("E_CHILD_COMMAND", "child argv must be non-empty strings without NUL")
        if popen_kwargs.pop("shell", False):
            raise ChildProcessFailure("E_CHILD_COMMAND", "shell=True is forbidden")
        if "start_new_session" in popen_kwargs:
            raise ChildProcessFailure("E_CHILD_COMMAND", "start_new_session is controlled by the registry")
        command_digest = sha256_bytes(canonical_json_bytes(list(argv)))
        process = subprocess.Popen(list(argv), start_new_session=True, **popen_kwargs)
        try:
            start_time = _pid_start_time(process.pid)
        except BaseException:
            process.terminate()
            process.wait()
            raise
        child = TrackedChild(
            label=label,
            argv=tuple(argv),
            command_sha256=command_digest,
            process=process,
            pid_start_time=start_time,
            launched_at_utc=utc_now(),
        )
        with self._lock:
            self._children[process.pid] = child
            self._persist()
        return child

    def _refresh(self) -> None:
        for child in self._children.values():
            if child.returncode is None:
                code = child.process.poll()
                if code is not None:
                    child.returncode = code

    def _persist(self) -> None:
        self._refresh()
        payload = {
            "schema_name": "intersubmod.layered_child_registry",
            "schema_version": "1.0.0",
            "updated_at_utc": utc_now(),
            "children": [
                {
                    "label": child.label,
                    "pid": child.process.pid,
                    "pid_start_time": child.pid_start_time,
                    "command_sha256": child.command_sha256,
                    "argv": list(child.argv),
                    "launched_at_utc": child.launched_at_utc,
                    "returncode": child.returncode,
                }
                for child in sorted(self._children.values(), key=lambda item: item.process.pid)
            ],
        }
        atomic_write_json(self._snapshot_path(), payload)

    def heartbeat_pid_start_times(self) -> Dict[str, str]:
        with self._lock:
            self._refresh()
            return {
                str(pid): child.pid_start_time
                for pid, child in sorted(self._children.items())
                if child.returncode is None
            }

    def active_count(self) -> int:
        with self._lock:
            self._refresh()
            return sum(child.returncode is None for child in self._children.values())

    def all_reaped(self) -> bool:
        return self.active_count() == 0

    def _same_process(self, child: TrackedChild) -> bool:
        if child.process.poll() is not None:
            child.returncode = child.process.returncode
            return False
        try:
            return _pid_start_time(child.process.pid) == child.pid_start_time
        except ContractError:
            return False

    def terminate_all(self, grace_seconds: float = 2.0) -> None:
        with self._lock:
            self._refresh()
            active = [child for child in self._children.values() if child.returncode is None]
            identity_errors: List[int] = []
            for child in active:
                if not self._same_process(child):
                    identity_errors.append(child.process.pid)
                    continue
                try:
                    os.killpg(child.process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            deadline = time.monotonic() + max(0.0, grace_seconds)
            while time.monotonic() < deadline and any(child.process.poll() is None for child in active):
                time.sleep(min(0.02, max(0.0, deadline - time.monotonic())))
            for child in active:
                if child.process.poll() is None:
                    if not self._same_process(child):
                        identity_errors.append(child.process.pid)
                        continue
                    try:
                        os.killpg(child.process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
            for child in active:
                try:
                    child.returncode = child.process.wait(timeout=2)
                except subprocess.TimeoutExpired as exc:
                    raise ChildProcessFailure(
                        "E_CHILD_NOT_REAPED",
                        f"child did not terminate: PID {child.process.pid}",
                    ) from exc
            self._persist()
            if identity_errors:
                raise ChildProcessFailure(
                    "E_PID_IDENTITY",
                    f"refused to signal PID identities: {sorted(set(identity_errors))}",
                )

    def wait_all_fail_fast(
        self,
        poll_seconds: float = 0.02,
        grace_seconds: float = 2.0,
    ) -> Mapping[str, int]:
        while True:
            with self._lock:
                self._refresh()
                failures = [child for child in self._children.values() if child.returncode not in (None, 0)]
                active = [child for child in self._children.values() if child.returncode is None]
            if failures:
                self.terminate_all(grace_seconds=grace_seconds)
                first = sorted(failures, key=lambda item: item.process.pid)[0]
                raise ChildProcessFailure(
                    "E_CHILD_FAILED",
                    f"child {first.label} PID {first.process.pid} exited {first.returncode}",
                )
            if not active:
                with self._lock:
                    self._persist()
                    return {child.label: int(child.returncode or 0) for child in self._children.values()}
            time.sleep(max(0.001, poll_seconds))


class HeartbeatWriter:
    def __init__(self, lifecycle: "RunLifecycle", interval_seconds: float):
        if interval_seconds <= 0:
            raise HeartbeatError("E_HEARTBEAT_INTERVAL", "heartbeat interval must be positive")
        self.lifecycle = lifecycle
        self.interval_seconds = interval_seconds
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._seq = 0
        self.error: Optional[str] = None

    def write_once(self, active_samples: Optional[Sequence[str]] = None) -> Mapping[str, Any]:
        self._seq += 1
        samples = list(active_samples if active_samples is not None else self.lifecycle.active_samples)
        payload = {
            "schema_name": "intersubmod.layered_heartbeat",
            "schema_version": "1.0.0",
            "seq": self._seq,
            "wall_time_utc": utc_now(),
            "monotonic_seconds": time.monotonic() - self.lifecycle.started_monotonic,
            "host": socket.gethostname(),
            "launcher_pid": os.getpid(),
            "launcher_pid_start_time": self.lifecycle.launcher_pid_start_time,
            "state": self.lifecycle.state,
            "stage": self.lifecycle.stage,
            "active_samples": samples,
            "child_pid_start_times": self.lifecycle.children.heartbeat_pid_start_times(),
            "last_event_seq": self.lifecycle.sequence,
            "frozen_lock_sha256": self.lifecycle.frozen_lock_sha256,
            "source_bundle_sha256": self.lifecycle.source_bundle_content_sha256,
        }
        atomic_write_json(self.lifecycle.root / "heartbeat.json", payload)
        return payload

    def _run(self) -> None:
        while not self._stop.wait(self.interval_seconds):
            try:
                self.write_once()
            except BaseException as exc:  # Preserve failure for the launcher to gate on.
                self.error = repr(exc)
                self._stop.set()

    def start(self) -> None:
        if self._thread is not None:
            raise HeartbeatError("E_HEARTBEAT_ALREADY_RUNNING", "heartbeat is already running")
        self.write_once()
        self._thread = threading.Thread(target=self._run, name="layered-v3-heartbeat", daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.interval_seconds * 2))
            if self._thread.is_alive():
                raise HeartbeatError("E_HEARTBEAT_STOP", "heartbeat thread did not stop")


def validate_heartbeat(
    path: PathLike,
    max_age_seconds: float = 180.0,
    expected_launcher_pid_start_time: Optional[str] = None,
    now: Optional[dt.datetime] = None,
) -> Mapping[str, Any]:
    payload, _ = read_json_with_sha256(path)
    required = {
        "schema_name",
        "schema_version",
        "seq",
        "wall_time_utc",
        "monotonic_seconds",
        "host",
        "launcher_pid",
        "launcher_pid_start_time",
        "state",
        "stage",
        "active_samples",
        "child_pid_start_times",
        "last_event_seq",
        "frozen_lock_sha256",
        "source_bundle_sha256",
    }
    if set(payload) != required:
        raise HeartbeatError("E_HEARTBEAT_INVALID", "heartbeat has missing or unknown fields")
    if payload["schema_name"] != "intersubmod.layered_heartbeat" or payload["schema_version"] != "1.0.0":
        raise HeartbeatError("E_HEARTBEAT_INVALID", "heartbeat schema mismatch")
    try:
        observed_time = dt.datetime.fromisoformat(str(payload["wall_time_utc"]).replace("Z", "+00:00"))
    except ValueError as exc:
        raise HeartbeatError("E_HEARTBEAT_INVALID", "heartbeat timestamp is invalid") from exc
    reference = now or dt.datetime.now(dt.timezone.utc)
    age = (reference - observed_time).total_seconds()
    if age > max_age_seconds or age < -5:
        raise HeartbeatError("E_HEARTBEAT_STALE", f"heartbeat age is {age:.3f} seconds")
    if (
        expected_launcher_pid_start_time is not None
        and payload["launcher_pid_start_time"] != expected_launcher_pid_start_time
    ):
        raise HeartbeatError("E_HEARTBEAT_STALE", "launcher PID start time does not match")
    return payload


LEGAL_TRANSITIONS: Mapping[str, Set[str]] = {
    "CREATING": {"PREFLIGHT", "FAILED", "ABORTED"},
    "PREFLIGHT": {"READY", "FAILED", "ABORTED"},
    "READY": {"RUNNING", "FAILED", "ABORTED"},
    "RUNNING": {"VERIFYING", "FAILED", "ABORTED"},
    "VERIFYING": {"SUCCEEDED", "FAILED", "ABORTED"},
    "SUCCEEDED": set(),
    "FAILED": set(),
    "ABORTED": set(),
}
TERMINAL_STATES = {"SUCCEEDED", "FAILED", "ABORTED"}


class RunLifecycle:
    """Own the v3 staging/publish boundary and all launcher lifecycle state."""

    def __init__(self, run_parent: PathLike, run_id: str, install_traps: bool = True):
        validate_run_id(run_id)
        self.run_id = run_id
        self.started_monotonic = time.monotonic()
        self.launcher_pid_start_time = _pid_start_time(os.getpid())
        self.lock = ParentRunLock(run_parent, run_id)
        self.lock.acquire()
        self.run_parent = self.lock.run_parent
        self.published_root = self.run_parent / run_id
        self.staging_root = self.run_parent / f".{run_id}.staging.{os.getpid()}.{uuid.uuid4().hex[:12]}"
        self._root = self.staging_root
        self._published = False
        self._closed = False
        self._terminalizing = False
        self._old_signal_handlers: Dict[int, Any] = {}
        self._atexit_registered = False
        self.state = "CREATING"
        self.stage = "bootstrap"
        self.sequence = 0
        self._state_digest = ZERO_DIGEST
        self.active_samples: List[str] = []
        self.frozen_lock_sha256: Optional[str] = None
        self.environment_lock_sha256: Optional[str] = None
        self.source_bundle_manifest_sha256: Optional[str] = None
        self.source_bundle_content_sha256: Optional[str] = None
        self.launch_receipt_sha256: Optional[str] = None
        self._source_bundle: Optional[SourceBundleResult] = None
        self._heartbeat: Optional[HeartbeatWriter] = None
        try:
            if self.published_root.exists():
                raise ContractError("E_RUN_ROOT_EXISTS", f"published run root already exists: {self.published_root}")
            self.staging_root.mkdir(mode=0o700)
            (self.staging_root / "state_events").mkdir(mode=0o700)
            self.children = TrackedChildRegistry(lambda: self.root / "children.json")
            self._write_state("CREATING", stage="bootstrap", reason="lifecycle initialized", error_code=None)
            if install_traps:
                self.install_traps()
        except BaseException:
            try:
                self._restore_traps()
            except BaseException:
                pass
            self.lock.release()
            raise

    @property
    def root(self) -> Path:
        return self._root

    @property
    def published(self) -> bool:
        return self._published

    def _state_record(
        self,
        new_state: str,
        stage: str,
        reason: str,
        error_code: Optional[str],
        sample: Optional[str],
    ) -> Dict[str, Any]:
        sequence = self.sequence + 1
        record: Dict[str, Any] = {
            "schema_name": "intersubmod.layered_run_state",
            "schema_version": "1.0.0",
            "run_id": self.run_id,
            "sequence": sequence,
            "previous_state_sha256": self._state_digest,
            "timestamp_utc": utc_now(),
            "state": new_state,
            "stage": stage,
            "sample": sample,
            "launcher_pid": os.getpid(),
            "launcher_pid_start_time": self.launcher_pid_start_time,
            "reason": reason,
            "error_code": error_code,
        }
        record["state_sha256"] = sha256_bytes(canonical_json_bytes(record))
        return record

    def _write_state(
        self,
        new_state: str,
        stage: str,
        reason: str,
        error_code: Optional[str],
        sample: Optional[str] = None,
    ) -> None:
        record = self._state_record(new_state, stage, reason, error_code, sample)
        event_path = self.root / "state_events" / f"{record['sequence']:06d}_{new_state}.json"
        atomic_write_json(event_path, record, no_overwrite=True)
        atomic_write_json(self.root / "run_state.json", record)
        self.sequence = int(record["sequence"])
        self._state_digest = str(record["state_sha256"])
        self.state = new_state
        self.stage = stage

    def transition(
        self,
        new_state: str,
        stage: str,
        reason: str,
        error_code: Optional[str] = None,
        sample: Optional[str] = None,
    ) -> None:
        if new_state not in LEGAL_TRANSITIONS.get(self.state, set()):
            raise StateTransitionError("E_STATE_TRANSITION", f"illegal transition {self.state} -> {new_state}")
        if "\n" in reason or "\r" in reason or "\t" in reason:
            raise StateTransitionError("E_STATE_DETAIL", "state reason may not contain tabs or newlines")
        _safe_label(stage, "stage")
        if sample is not None:
            _safe_label(sample, "sample")
        self._write_state(new_state, stage, reason, error_code, sample)

    def begin_preflight(self) -> None:
        self.transition("PREFLIGHT", "preflight", "preflight started")

    def write_frozen_lock(self, payload: Mapping[str, Any]) -> str:
        if self.state != "PREFLIGHT":
            raise StateTransitionError("E_STATE_TRANSITION", "frozen lock can only be written during PREFLIGHT")
        if self.frozen_lock_sha256 is not None:
            raise ContractError("E_FROZEN_LOCK_EXISTS", "frozen lock has already been written")
        digest = atomic_write_json(self.root / "frozen_input_lock.json", payload, no_overwrite=True)
        self.frozen_lock_sha256 = digest
        return digest

    def build_source_bundle(
        self,
        runner: PathLike,
        validator: PathLike,
        verifier: PathLike,
        imported_scripts: Sequence[PathLike],
    ) -> SourceBundleResult:
        if self.state != "PREFLIGHT":
            raise StateTransitionError("E_STATE_TRANSITION", "source bundle can only be built during PREFLIGHT")
        if self._source_bundle is not None:
            raise SourceBundleError("E_SOURCE_BUNDLE_EXISTS", "source bundle has already been built")
        bundle = build_source_bundle(self.root / "source_bundle", runner, validator, verifier, imported_scripts)
        self._source_bundle = bundle
        self.source_bundle_manifest_sha256 = bundle.manifest_sha256
        self.source_bundle_content_sha256 = bundle.content_sha256
        return bundle

    def capture_environment_lock(self, **kwargs: Any) -> str:
        if self.state != "PREFLIGHT":
            raise StateTransitionError("E_STATE_TRANSITION", "environment lock can only be written during PREFLIGHT")
        if self.environment_lock_sha256 is not None:
            raise EnvironmentContractError("E_ENVIRONMENT_LOCK_EXISTS", "environment lock has already been written")
        digest = capture_environment_lock(self.root / "environment_lock.json", **kwargs)
        self.environment_lock_sha256 = digest
        return digest

    def _verify_preflight_artifacts(self) -> None:
        required = {
            "frozen_lock_sha256": self.frozen_lock_sha256,
            "environment_lock_sha256": self.environment_lock_sha256,
            "source_bundle_manifest_sha256": self.source_bundle_manifest_sha256,
            "source_bundle_content_sha256": self.source_bundle_content_sha256,
        }
        missing = sorted(name for name, value in required.items() if not value)
        if missing:
            raise ContractError("E_PREFLIGHT_INCOMPLETE", f"missing preflight artifacts: {missing}")
        if sha256_file(self.root / "frozen_input_lock.json", reject_symlink=True) != self.frozen_lock_sha256:
            raise ContractError("E_FROZEN_LOCK_MISMATCH", "frozen lock changed before publish")
        if sha256_file(self.root / "environment_lock.json", reject_symlink=True) != self.environment_lock_sha256:
            raise EnvironmentContractError("E_ENVIRONMENT_MISMATCH", "environment lock changed before publish")
        if self._source_bundle is None:
            raise SourceBundleError("E_SOURCE_BUNDLE_MISMATCH", "source bundle metadata is absent")
        verify_source_bundle(
            self._source_bundle.manifest_path,
            expected_manifest_sha256=self.source_bundle_manifest_sha256,
            expected_content_sha256=self.source_bundle_content_sha256,
        )

    def publish_ready(self, receipt_extra: Optional[Mapping[str, Any]] = None) -> Path:
        if self.state != "PREFLIGHT":
            raise StateTransitionError("E_STATE_TRANSITION", "publish requires PREFLIGHT state")
        if self.published_root.exists():
            raise ContractError("E_RUN_ROOT_EXISTS", f"published run root already exists: {self.published_root}")
        self._verify_preflight_artifacts()
        receipt = {
            "schema_name": "intersubmod.layered_run_receipt",
            "schema_version": "1.0.0",
            "run_id": self.run_id,
            "created_at_utc": utc_now(),
            "launcher_pid": os.getpid(),
            "launcher_pid_start_time": self.launcher_pid_start_time,
            "frozen_lock_sha256": self.frozen_lock_sha256,
            "environment_lock_sha256": self.environment_lock_sha256,
            "source_bundle_manifest_sha256": self.source_bundle_manifest_sha256,
            "source_bundle_content_sha256": self.source_bundle_content_sha256,
            "extra": dict(receipt_extra or {}),
        }
        self.launch_receipt_sha256 = atomic_write_json(
            self.root / "launch_receipt.json", receipt, no_overwrite=True
        )
        self.transition("READY", "ready", "preflight artifacts frozen and launch receipt written")
        _fsync_directory(self.staging_root)
        os.rename(str(self.staging_root), str(self.published_root))
        _fsync_directory(self.run_parent)
        self._root = self.published_root
        self._published = True
        if self._source_bundle is not None:
            old_bundle = self._source_bundle
            new_bundle_root = self.root / "source_bundle"
            relocated_roles = {
                role: new_bundle_root / path.relative_to(old_bundle.root)
                for role, path in old_bundle.role_paths.items()
            }
            self._source_bundle = SourceBundleResult(
                new_bundle_root,
                new_bundle_root / old_bundle.manifest_path.name,
                old_bundle.manifest_sha256,
                old_bundle.content_sha256,
                relocated_roles,
            )
        return self.root

    def begin_running(self) -> None:
        if not self._published:
            raise StateTransitionError("E_STATE_TRANSITION", "RUNNING requires a published run root")
        self.transition("RUNNING", "running", "workers may now start")

    def begin_verifying(self) -> None:
        if self.children.active_count() != 0:
            raise ChildProcessFailure("E_CHILD_ACTIVE", "cannot verify while tracked children are active")
        self.transition("VERIFYING", "verifying", "independent verification started")

    def set_active_samples(self, samples: Sequence[str]) -> None:
        for sample in samples:
            _safe_label(sample, "sample")
        self.active_samples = list(samples)

    def start_heartbeat(self, interval_seconds: float = 60.0) -> HeartbeatWriter:
        if self.state not in {"READY", "RUNNING", "VERIFYING"} or not self._published:
            raise HeartbeatError("E_HEARTBEAT_STATE", "heartbeat requires a published READY/RUNNING/VERIFYING run")
        if self._heartbeat is not None:
            raise HeartbeatError("E_HEARTBEAT_ALREADY_RUNNING", "heartbeat already exists")
        heartbeat = HeartbeatWriter(self, interval_seconds)
        heartbeat.start()
        self._heartbeat = heartbeat
        return heartbeat

    def stop_heartbeat(self) -> None:
        if self._heartbeat is None:
            return
        self._heartbeat.stop()
        if self._heartbeat.error:
            raise HeartbeatError("E_HEARTBEAT_WRITE", self._heartbeat.error)
        self._heartbeat = None

    def _stop_runtime(self, grace_seconds: float) -> None:
        heartbeat_error: Optional[BaseException] = None
        try:
            self.stop_heartbeat()
        except BaseException as exc:
            heartbeat_error = exc
        child_error: Optional[BaseException] = None
        try:
            if self.children.active_count():
                self.children.terminate_all(grace_seconds=grace_seconds)
        except BaseException as exc:
            child_error = exc
        if heartbeat_error:
            raise heartbeat_error
        if child_error:
            raise child_error

    def fail(self, error_code: str, reason: str, grace_seconds: float = 2.0) -> None:
        if self.state in TERMINAL_STATES or self._terminalizing:
            return
        self._terminalizing = True
        runtime_note = ""
        try:
            try:
                self._stop_runtime(grace_seconds)
            except BaseException as exc:
                runtime_note = f"; runtime cleanup error={exc!r}"
            self.transition(
                "FAILED",
                "failed",
                _sanitize_detail(reason + runtime_note),
                error_code=error_code,
            )
        finally:
            self._terminalizing = False

    def abort(self, signum: int, grace_seconds: float = 2.0) -> None:
        if self.state in TERMINAL_STATES or self._terminalizing:
            return
        self._terminalizing = True
        try:
            try:
                self._stop_runtime(grace_seconds)
                cleanup_note = ""
            except BaseException as exc:
                cleanup_note = f"; runtime cleanup error={exc!r}"
            try:
                signal_name = signal.Signals(signum).name
            except ValueError:
                signal_name = str(signum)
            self.transition(
                "ABORTED",
                "aborted",
                _sanitize_detail(f"received {signal_name}{cleanup_note}"),
                error_code=f"E_SIGNAL_{signal_name}",
            )
        finally:
            self._terminalizing = False

    def succeed(self, verification_artifact: PathLike, success_extra: Optional[Mapping[str, Any]] = None) -> Path:
        if self.state != "VERIFYING":
            raise StateTransitionError("E_STATE_TRANSITION", "success requires VERIFYING state")
        self.stop_heartbeat()
        if self.children.active_count() != 0:
            raise ChildProcessFailure("E_CHILD_ACTIVE", "cannot succeed while tracked children are active")
        self._verify_preflight_artifacts()
        if not self.launch_receipt_sha256 or sha256_file(
            self.root / "launch_receipt.json", reject_symlink=True
        ) != self.launch_receipt_sha256:
            raise ContractError("E_LAUNCH_RECEIPT_MISMATCH", "launch receipt changed after publish")
        verification_path = Path(verification_artifact)
        verification_sha = sha256_file(verification_path, reject_symlink=True)
        if (self.root / "_SUCCESS").exists():
            raise ContractError("E_SUCCESS_EXISTS", "_SUCCESS already exists")
        self.transition("SUCCEEDED", "succeeded", "verification and provenance readback passed")
        marker = {
            "schema_name": "intersubmod.layered_success_marker",
            "schema_version": "1.0.0",
            "run_id": self.run_id,
            "created_at_utc": utc_now(),
            "state_sha256": self._state_digest,
            "launch_receipt_sha256": self.launch_receipt_sha256,
            "frozen_lock_sha256": self.frozen_lock_sha256,
            "environment_lock_sha256": self.environment_lock_sha256,
            "source_bundle_manifest_sha256": self.source_bundle_manifest_sha256,
            "source_bundle_content_sha256": self.source_bundle_content_sha256,
            "verification_path": str(verification_path.resolve()),
            "verification_sha256": verification_sha,
            "extra": dict(success_extra or {}),
        }
        marker_path = self.root / "_SUCCESS"
        atomic_write_json(marker_path, marker, no_overwrite=True)
        return marker_path

    def handle_signal(self, signum: int, _frame: Any = None) -> None:
        raise LifecycleSignal(signum)

    def install_traps(self) -> None:
        if self._old_signal_handlers:
            return
        if threading.current_thread() is not threading.main_thread():
            raise ContractError("E_SIGNAL_THREAD", "signal traps must be installed by the main thread")
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            self._old_signal_handlers[signum] = signal.getsignal(signum)
            signal.signal(signum, self.handle_signal)
        atexit.register(self._atexit_guard)
        self._atexit_registered = True

    def _restore_traps(self) -> None:
        for signum, handler in self._old_signal_handlers.items():
            signal.signal(signum, handler)
        self._old_signal_handlers.clear()
        if self._atexit_registered:
            try:
                atexit.unregister(self._atexit_guard)
            except AttributeError:
                pass
            self._atexit_registered = False

    def _atexit_guard(self) -> None:
        if self._closed:
            return
        try:
            if self.state not in TERMINAL_STATES:
                self.fail("E_PROCESS_EXIT", "process exited without terminal success", grace_seconds=0.2)
        finally:
            self.lock.release()

    def close(self) -> None:
        if self._closed:
            return
        if self.state not in TERMINAL_STATES:
            self.fail("E_INCOMPLETE_EXIT", "lifecycle closed without terminal state", grace_seconds=0.2)
        self._restore_traps()
        self.lock.release()
        self._closed = True

    def __enter__(self) -> "RunLifecycle":
        return self

    def __exit__(self, exc_type: Any, exc: Any, _traceback: Any) -> bool:
        try:
            if isinstance(exc, LifecycleSignal):
                self.abort(exc.signum, grace_seconds=0.2)
            elif exc is not None and self.state not in TERMINAL_STATES:
                code = exc.code if isinstance(exc, ContractError) else "E_UNHANDLED_EXCEPTION"
                self.fail(code, str(exc), grace_seconds=0.2)
            elif exc is None and self.state not in TERMINAL_STATES:
                self.fail("E_INCOMPLETE_EXIT", "context exited without success", grace_seconds=0.2)
        finally:
            self.close()
        return False


def _main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Layered v3 lifecycle contract utilities")
    subparsers = parser.add_subparsers(dest="command", required=True)
    source_parser = subparsers.add_parser("verify-source-bundle")
    source_parser.add_argument("manifest")
    source_parser.add_argument("--expected-manifest-sha256")
    source_parser.add_argument("--expected-content-sha256")
    heartbeat_parser = subparsers.add_parser("validate-heartbeat")
    heartbeat_parser.add_argument("heartbeat")
    heartbeat_parser.add_argument("--max-age-seconds", type=float, default=180.0)
    args = parser.parse_args(argv)
    try:
        if args.command == "verify-source-bundle":
            result = verify_source_bundle(
                args.manifest,
                expected_manifest_sha256=args.expected_manifest_sha256,
                expected_content_sha256=args.expected_content_sha256,
            )
            output = {
                "pass": True,
                "manifest_sha256": result.manifest_sha256,
                "content_sha256": result.content_sha256,
            }
        else:
            payload = validate_heartbeat(args.heartbeat, max_age_seconds=args.max_age_seconds)
            output = {"pass": True, "seq": payload["seq"], "state": payload["state"]}
        print(json.dumps(output, sort_keys=True))
        return 0
    except ContractError as exc:
        print(
            json.dumps({"pass": False, "error_code": exc.code, "message": exc.message}, sort_keys=True),
            file=sys.stderr,
        )
        return 6


if __name__ == "__main__":
    raise SystemExit(_main())
