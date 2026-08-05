#!/usr/bin/env python3
"""Layered-v3 fail-closed launcher and frozen-source worker adapter.

This file is the only new execution entry point for layered v3.  It does not
modify or source the live v2 shell runner.  The bootstrap process owns the
parent flock and PREFLIGHT state, freezes all executable bytes, runs the
bundled validator, and publishes the run root only after every launch gate has
passed.  Scientific workers execute this bundled runner and load only bundled
scientific modules; sample inputs are resolved exclusively from
``frozen_input_lock.json``.

The production command is intentionally not started by this module's tests.
Use ``--preflight-only`` for a terminal, clearly labelled preflight receipt
with zero scientific workers.  A normal invocation performs the full
7-dataset/6-biological-sample chr1-22 workflow and calls the bundled verifier
before ``RunLifecycle.succeed`` creates ``_SUCCESS``.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import hashlib
import importlib.abc
import importlib.util
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Any, Callable, Mapping, Optional, Sequence


def _load_lifecycle_module() -> Any:
    """Import lifecycle from the worktree or its prefixed frozen-bundle path."""

    try:
        import layered_v3_lifecycle as module

        return module
    except ImportError:
        here = Path(__file__).resolve()
        candidates = sorted((here.parent.parent / "imported").glob("*_layered_v3_lifecycle.py"))
        if len(candidates) != 1:
            raise RuntimeError(
                f"E_LIFECYCLE_IMPORT: expected one bundled lifecycle module; found {candidates}"
            )
        spec = importlib.util.spec_from_file_location("layered_v3_lifecycle", candidates[0])
        if spec is None or spec.loader is None:
            raise RuntimeError(f"E_LIFECYCLE_IMPORT: cannot load {candidates[0]}")
        module = importlib.util.module_from_spec(spec)
        sys.modules["layered_v3_lifecycle"] = module
        spec.loader.exec_module(module)
        return module


L = _load_lifecycle_module()
ContractError = L.ContractError
RunLifecycle = L.RunLifecycle
ToolSpec = L.ToolSpec

REPO_ROOT = Path(__file__).resolve().parents[1]
METHOD_ROOT = REPO_ROOT / "docs/methodology/_assets/20260627_subclone_4axis_teaching"
METHOD_SCRIPTS = METHOD_ROOT / "scripts"
SCHEMA_ROOT = METHOD_ROOT / "schemas"

EXPECTED_BINDING = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
CANONICAL_TREE_INPUT_CONTRACT = "longphase_s_recalibrated_FILTER_PASS"
SENSITIVITY_TREE_INPUT_CONTRACT = "clairs_FILTER_PASS_sensitivity"
TREE_PROFILES = {
    CANONICAL_TREE_INPUT_CONTRACT: {
        "task_type": "comprehensive_validation",
        "backbone_role": "longphase_s_recalibrated_filter_pass",
        "ledger_tree_contract": "longphase_recalibrated_PASS",
    },
    SENSITIVITY_TREE_INPUT_CONTRACT: {
        "task_type": "backbone_sensitivity",
        "backbone_role": "clairs_filter_pass_sensitivity",
        "ledger_tree_contract": "clairs_PASS_input",
    },
}
PARTS = {
    1: ("chr1", "chr6", "chr11", "chr16", "chr21"),
    2: ("chr2", "chr7", "chr12", "chr17", "chr22"),
    3: ("chr3", "chr8", "chr13", "chr18"),
    4: ("chr4", "chr9", "chr14", "chr19"),
    5: ("chr5", "chr10", "chr15", "chr20"),
}
SCIENCE_BASENAMES = (
    "sm_linkage_genomewide.py",
    "sm_multilocus_combinations.py",
    "tree_enumeration_solver.py",
    "layered_tree_reconstruction.py",
    "build_region_view.py",
    "build_ssnv_site_ledger.py",
)
SCHEMA_BASENAMES = (
    "layered_input_manifest_v3.schema.json",
    "layered_input_lock_v1.schema.json",
    "longphase_raw_all_capture_receipt_v2.schema.json",
)
SCIENTIFIC_ROLES = tuple(f"mlhp_part_{number}" for number in range(1, 6)) + (
    "layered_reconstruction",
    "layered_region_view",
    "site_ledger",
    "site_ledger_summary",
)
PROVENANCE_FIELDS = (
    "frozen_lock_sha256",
    "launch_receipt_sha256",
    "environment_lock_sha256",
    "source_bundle_manifest_sha256",
    "source_bundle_content_sha256",
    "input_set_sha256",
)
FROZEN_CONTRACT_SHA256 = {
    "layered_input_manifest_v3.schema.json": "47f54d86159beae66b68719e034c7755eba59b2d1719d4546153fa0ae74ff19f",
    "layered_input_lock_v1.schema.json": "1a40107f696e19e375836247587c728fec781a1f6c69c34ad60cbec8d72a4fdb",
    "longphase_raw_all_capture_receipt_v2.schema.json": "93bee1231a646f66ead1d31c0f60df5253df671a914606f4018be39d22616afe",
    "validate_layered_v3_inputs.py": "d4e69586be124b337dd9d687fef1c544f7bdbe6288bdc74cf8fe4d2146c04c79",
    "verify_layered_v3.py": "ef7a4893b76ea219df7dce33f58ec0278ec50d3cb5703b16d987f68aead816bd",
}
SM_ALLOWLIST = frozenset(
    {
        "SM_INPUT_MANIFEST",
        "SM_RUN_ID",
        "SM_RUN_PARENT",
        "SM_PARALLEL_SAMPLES",
        "SM_VERIFY_EVERY",
        "SM_ANALYSIS_TREE_CAP",
        "SM_DISPLAY_TREE_CAP",
        "SM_MINREAD",
        "SM_MAX_SNV",
        "SM_TIER_R",
        "SM_MAPQ_MIN",
        "SM_BASEQ_MIN",
        "SM_HEARTBEAT_INTERVAL",
        "SM_MIN_LOGICAL_CPUS",
        "SM_MIN_AVAILABLE_MEMORY_GIB",
        "SM_MIN_FREE_DISK_GIB",
        "SM_MAX_LOAD_PER_CPU",
        "SM_MAX_IOWAIT_PERCENT",
        "SM_RESOURCE_SAMPLE_SECONDS",
        "SM_MAX_NFS_READ_MBPS",
        "SM_NFS_MOUNTPOINT",
    }
)
GLOBAL_SCOPE_LOCK_NAME = ".layered_chr1_22_7dataset_full.lock"
CONFLICTING_PROCESS_BASENAMES = frozenset(
    {
        "run_layered_7samples_newbb.sh",
        "run_layered_v3.py",
        "sm_multilocus_combinations.py",
        "layered_tree_reconstruction.py",
        "build_region_view.py",
        "build_ssnv_site_ledger.py",
        "longphase-s",
        "run_longphase_raw_all_production_sidecars.sh",
        "capture_longphase_tagged_bam_sidecar.py",
        "validate_streamed_longphase_sidecar.py",
        "finalize_longphase_raw_all_capture_receipts.py",
        "build_longphase_raw_all_capture_receipt_v2.py",
    }
)
PRODUCTION_RESOURCE_POLICY = {
    "min_logical_cpus": 8,
    "min_available_memory_gib": 128.0,
    "min_free_disk_gib": 500.0,
    "max_load_per_cpu": 1.25,
    "max_iowait_percent": 20.0,
    "resource_sample_seconds": 300.0,
    "max_nfs_read_mbps": 80.0,
    "nfs_mountpoint": "/big8_disk",
}


class GlobalScopeLock:
    """One cooperative full-scope lock shared by every run ID."""

    def __init__(self, run_parent: Path, run_id: str):
        L.validate_run_id(run_id)
        self.run_parent = run_parent.resolve()
        self.run_parent.mkdir(parents=True, exist_ok=True)
        self.path = self.run_parent / GLOBAL_SCOPE_LOCK_NAME
        self.run_id = run_id
        self._fd: Optional[int] = None

    def acquire(self) -> None:
        if self._fd is not None:
            return
        if self.path.is_symlink():
            raise ContractError("E_GLOBAL_RUN_LOCKED", f"global lock path may not be a symlink: {self.path}")
        flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        try:
            fd = os.open(self.path, flags, 0o600)
        except OSError as exc:
            raise ContractError("E_GLOBAL_RUN_LOCKED", f"cannot safely open global lock {self.path}: {exc}") from exc
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            os.lseek(fd, 0, os.SEEK_SET)
            owner = os.read(fd, 4096).decode("utf-8", errors="replace").strip()
            os.close(fd)
            raise ContractError(
                "E_GLOBAL_RUN_LOCKED",
                f"another chr1-22/7-dataset launcher holds {self.path}; owner={owner!r}",
            ) from exc
        payload = (
            f"run_id={self.run_id} pid={os.getpid()} "
            f"pid_start_time={L._pid_start_time(os.getpid())} acquired_at_utc={L.utc_now()}\n"
        ).encode()
        os.ftruncate(fd, 0)
        os.write(fd, payload)
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

    def __enter__(self) -> "GlobalScopeLock":
        self.acquire()
        return self

    def __exit__(self, _kind: Any, _value: Any, _traceback: Any) -> bool:
        self.release()
        return False


@dataclass(frozen=True)
class SourceInventory:
    runner: Path
    lifecycle: Path
    validator: Path
    verifier: Path
    scientific: tuple[Path, ...]
    schemas: tuple[Path, ...]
    distributions: tuple[str, ...] = ("pysam",)
    tools: tuple[ToolSpec, ...] = ()
    expected_sha256: tuple[tuple[str, str], ...] = ()

    @property
    def imported(self) -> tuple[Path, ...]:
        return (self.lifecycle, *self.scientific, *self.schemas)


@dataclass(frozen=True)
class RunConfig:
    manifest: Path
    run_parent: Path
    run_id: str
    parallel_samples: int = 2
    verify_every: int = 1
    analysis_tree_cap: int = 0
    display_tree_cap: int = 32
    minread: int = 3
    max_snv: int = 8
    tier_r: int = 50000
    mapq_min: int = 20
    baseq_min: int = 0
    heartbeat_interval: float = 60.0
    preflight_only: bool = False
    min_logical_cpus: int = 8
    min_available_memory_gib: float = 128.0
    min_free_disk_gib: float = 500.0
    max_load_per_cpu: float = 1.25
    max_iowait_percent: float = 20.0
    resource_sample_seconds: float = 300.0
    max_nfs_read_mbps: float = 80.0
    nfs_mountpoint: Path = Path("/big8_disk")

    def params(self) -> dict[str, Any]:
        return {
            "scope": "chr1-22",
            "contigs": list(AUTOSOMES),
            "dataset_count": 7,
            "biological_sample_count": 6,
            "parallel_samples": self.parallel_samples,
            "parallel_parts_per_sample": 1,
            "VERIFY_EVERY": self.verify_every,
            "ANALYSIS_TREE_CAP": self.analysis_tree_cap,
            "DISPLAY_TREE_CAP": self.display_tree_cap,
            "MINREAD": self.minread,
            "MAX_SNV": self.max_snv,
            "TIER_R": self.tier_r,
            "MAPQ_MIN": self.mapq_min,
            "BASEQ_MIN": self.baseq_min,
            "resource_policy": {
                "min_logical_cpus": self.min_logical_cpus,
                "min_available_memory_gib": self.min_available_memory_gib,
                "min_free_disk_gib": self.min_free_disk_gib,
                "max_load_per_cpu": self.max_load_per_cpu,
                "max_iowait_percent": self.max_iowait_percent,
                "sample_seconds": self.resource_sample_seconds,
                "max_nfs_read_mbps": self.max_nfs_read_mbps,
                "nfs_mountpoint": str(self.nfs_mountpoint),
            },
        }


def production_inventory() -> SourceInventory:
    return SourceInventory(
        runner=Path(__file__).resolve(),
        lifecycle=(REPO_ROOT / "scripts/layered_v3_lifecycle.py").resolve(),
        validator=(METHOD_SCRIPTS / "validate_layered_v3_inputs.py").resolve(),
        verifier=(REPO_ROOT / "scripts/verify_layered_v3.py").resolve(),
        scientific=tuple((METHOD_SCRIPTS / name).resolve() for name in SCIENCE_BASENAMES),
        schemas=tuple((SCHEMA_ROOT / name).resolve() for name in SCHEMA_BASENAMES),
        expected_sha256=tuple(sorted(FROZEN_CONTRACT_SHA256.items())),
    )


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path, *, reject_symlink: bool = False) -> str:
    return L.sha256_file(path, reject_symlink=reject_symlink)


def _no_duplicate_pairs(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError("E_JSON_DUPLICATE_KEY", f"duplicate key {key!r}")
        result[key] = value
    return result


def load_json_strict(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(
            path.read_text(encoding="utf-8", errors="strict"),
            object_pairs_hook=_no_duplicate_pairs,
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-finite JSON constant: {value}")
            ),
        )
    except ContractError:
        raise
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as exc:
        raise ContractError("E_JSON_INVALID", f"cannot parse strict JSON {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise ContractError("E_JSON_NOT_OBJECT", f"expected JSON object: {path}")
    return document


def _validate_config(config: RunConfig) -> None:
    L.validate_run_id(config.run_id)
    if not config.manifest.is_file() or config.manifest.is_symlink():
        raise ContractError("E_FILE_MISSING", f"manifest is missing/non-regular: {config.manifest}")
    if not 1 <= config.parallel_samples <= 4:
        raise ContractError("E_PARAM_RANGE", "parallel_samples must be between 1 and 4")
    if config.verify_every != 1 or config.analysis_tree_cap != 0:
        raise ContractError(
            "E_PARAM_PROFILE",
            "comprehensive mode requires VERIFY_EVERY=1 and ANALYSIS_TREE_CAP=0",
        )
    integer_bounds = {
        "DISPLAY_TREE_CAP": (config.display_tree_cap, 1, 1000000),
        "MINREAD": (config.minread, 1, 1000000),
        "MAX_SNV": (config.max_snv, 2, 1000),
        "TIER_R": (config.tier_r, 1, 1000000000),
        "MAPQ_MIN": (config.mapq_min, 0, 255),
        "BASEQ_MIN": (config.baseq_min, 0, 255),
    }
    for name, (value, lower, upper) in integer_bounds.items():
        if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
            raise ContractError("E_PARAM_RANGE", f"{name} must be in [{lower}, {upper}]")
    if config.heartbeat_interval <= 0:
        raise ContractError("E_PARAM_RANGE", "heartbeat interval must be positive")
    if config.min_logical_cpus < 1:
        raise ContractError("E_PARAM_RANGE", "min_logical_cpus must be positive")
    for name, value in (
        ("min_available_memory_gib", config.min_available_memory_gib),
        ("min_free_disk_gib", config.min_free_disk_gib),
        ("max_load_per_cpu", config.max_load_per_cpu),
        ("resource_sample_seconds", config.resource_sample_seconds),
        ("max_nfs_read_mbps", config.max_nfs_read_mbps),
    ):
        if value < 0:
            raise ContractError("E_PARAM_RANGE", f"{name} must be non-negative")
    if not 0 <= config.max_iowait_percent <= 100:
        raise ContractError("E_PARAM_RANGE", "max_iowait_percent must be in [0, 100]")
    if not config.nfs_mountpoint.is_absolute():
        raise ContractError("E_PARAM_RANGE", "nfs_mountpoint must be absolute")


def _validate_production_resource_policy(config: RunConfig, inventory: SourceInventory) -> None:
    if not inventory.expected_sha256:
        return
    observed = {
        "min_logical_cpus": config.min_logical_cpus,
        "min_available_memory_gib": config.min_available_memory_gib,
        "min_free_disk_gib": config.min_free_disk_gib,
        "max_load_per_cpu": config.max_load_per_cpu,
        "max_iowait_percent": config.max_iowait_percent,
        "resource_sample_seconds": config.resource_sample_seconds,
        "max_nfs_read_mbps": config.max_nfs_read_mbps,
        "nfs_mountpoint": str(config.nfs_mountpoint),
    }
    if observed != PRODUCTION_RESOURCE_POLICY:
        raise ContractError(
            "E_RESOURCE_POLICY_DRIFT",
            f"production resource policy is frozen; observed={observed} expected={PRODUCTION_RESOURCE_POLICY}",
        )


def _validate_source_inventory(inventory: SourceInventory) -> None:
    paths = (
        inventory.runner,
        inventory.lifecycle,
        inventory.validator,
        inventory.verifier,
        *inventory.scientific,
        *inventory.schemas,
    )
    by_name: dict[str, Path] = {}
    for path in paths:
        if path.is_symlink() or not path.is_file():
            raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"source is missing/symlink: {path}")
        if path.name in by_name:
            raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"duplicate source basename: {path.name}")
        by_name[path.name] = path
    if {path.name for path in inventory.scientific} != set(SCIENCE_BASENAMES):
        raise ContractError("E_SOURCE_BUNDLE_MISMATCH", "scientific source inventory is incomplete/extra")
    if {path.name for path in inventory.schemas} != set(SCHEMA_BASENAMES):
        raise ContractError("E_SOURCE_BUNDLE_MISMATCH", "schema source inventory is incomplete/extra")
    for basename, expected in inventory.expected_sha256:
        if basename not in by_name:
            raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"pinned source absent: {basename}")
        observed = file_sha256(by_name[basename], reject_symlink=True)
        if observed != expected:
            raise ContractError(
                "E_FROZEN_CONTRACT_DRIFT",
                f"reviewed contract changed: {basename}; expected={expected} observed={observed}",
            )


def _deterministic_environment(source: Optional[Mapping[str, str]] = None) -> dict[str, str]:
    env = dict(os.environ if source is None else source)
    required = {"LC_ALL": "C", "TZ": "UTC", "PYTHONHASHSEED": "0"}
    mismatches = {name: env.get(name) for name, expected in required.items() if env.get(name) != expected}
    if mismatches:
        raise ContractError(
            "E_ENVIRONMENT_MISMATCH",
            f"launch with LC_ALL=C TZ=UTC PYTHONHASHSEED=0; observed {mismatches}",
        )
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def _controlled_child_environment(base: Optional[Mapping[str, str]] = None) -> dict[str, str]:
    env = _deterministic_environment(base)
    for name in list(env):
        if name.startswith("SM_"):
            del env[name]
    return env


def _cpu_times() -> tuple[int, int]:
    fields = Path("/proc/stat").read_text(encoding="utf-8").splitlines()[0].split()
    if not fields or fields[0] != "cpu" or len(fields) < 6:
        raise ContractError("E_RESOURCE_GATE", "cannot parse aggregate /proc/stat CPU row")
    values = [int(value) for value in fields[1:]]
    return sum(values), values[4]


def _memory_available_bytes() -> int:
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if line.startswith("MemAvailable:"):
            fields = line.split()
            return int(fields[1]) * 1024
    raise ContractError("E_RESOURCE_GATE", "MemAvailable is absent from /proc/meminfo")


def _nfs_read_bytes(mountpoint: Path, mountstats: Path = Path("/proc/self/mountstats")) -> dict[str, Any]:
    """Read NFS per-op READ bytes_recv for one exact mountpoint."""

    try:
        lines = mountstats.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ContractError("E_NFS_BASELINE_UNAVAILABLE", f"cannot read {mountstats}: {exc}") from exc
    matches: list[dict[str, Any]] = []
    current: Optional[dict[str, Any]] = None
    for line in lines:
        if line.startswith("device "):
            try:
                device, remainder = line[len("device ") :].split(" mounted on ", 1)
                mounted, remainder = remainder.split(" with fstype ", 1)
                fstype = remainder.split()[0]
            except (ValueError, IndexError) as exc:
                raise ContractError("E_NFS_BASELINE_UNAVAILABLE", f"malformed mountstats header: {line}") from exc
            current = {
                "device": device,
                "mountpoint": mounted,
                "fstype": fstype,
                "read_bytes_recv": None,
            }
            if mounted == str(mountpoint):
                matches.append(current)
            continue
        if current is None or current.get("mountpoint") != str(mountpoint):
            continue
        stripped = line.strip()
        if stripped.startswith("READ:"):
            fields = stripped.split()
            try:
                values = [int(value) for value in fields[1:]]
            except ValueError as exc:
                raise ContractError("E_NFS_BASELINE_UNAVAILABLE", f"non-integer NFS READ row: {line}") from exc
            # Linux mountstats v1.x per-op fields:
            # ops, transmissions, timeouts, bytes_sent, bytes_recv,
            # queue, rtt, execute, errors.
            if len(values) < 9:
                raise ContractError("E_NFS_BASELINE_UNAVAILABLE", f"short NFS READ row: {line}")
            current["read_bytes_recv"] = values[4]
    if len(matches) != 1:
        raise ContractError(
            "E_NFS_BASELINE_UNAVAILABLE",
            f"expected one mountstats stanza for {mountpoint}; observed {len(matches)}",
        )
    result = matches[0]
    if not str(result["fstype"]).startswith("nfs") or result["read_bytes_recv"] is None:
        raise ContractError(
            "E_NFS_BASELINE_UNAVAILABLE",
            f"mount {mountpoint} is not NFS or lacks per-op READ bytes_recv: {result}",
        )
    result["counter_source"] = f"{mountstats}: per-op READ bytes_recv field"
    return result


def _default_process_snapshot() -> list[dict[str, Any]]:
    snapshot: list[dict[str, Any]] = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        if pid == os.getpid():
            continue
        try:
            raw = (entry / "cmdline").read_bytes()
            if not raw:
                continue
            argv = [part.decode("utf-8", errors="replace") for part in raw.split(b"\0") if part]
            snapshot.append({"pid": pid, "pid_start_time": L._pid_start_time(pid), "argv": argv})
        except (OSError, ContractError):
            continue
    return snapshot


def _bounded_sleep(seconds: float) -> None:
    deadline = time.monotonic() + seconds
    while True:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return
        time.sleep(min(30.0, remaining))


def _classify_process_snapshot(
    snapshot: Sequence[Mapping[str, Any]], observation_point: str
) -> dict[str, Any]:
    normalized: list[dict[str, Any]] = []
    conflicts: list[dict[str, Any]] = []
    for process in snapshot:
        try:
            pid = int(process["pid"])
            start_time = str(process["pid_start_time"])
            argv = list(process["argv"])
            if not argv or not all(isinstance(token, str) and token for token in argv):
                raise ValueError("invalid argv")
        except (KeyError, TypeError, ValueError) as exc:
            raise ContractError("E_PROCESS_OBSERVATION", f"malformed process snapshot entry: {process}") from exc
        record = {
            "pid": pid,
            "pid_start_time": start_time,
            "argv": argv,
            "argv_sha256": canonical_sha256(argv),
        }
        normalized.append(record)
        matched = sorted(CONFLICTING_PROCESS_BASENAMES & {Path(token).name for token in argv})
        if matched:
            conflicts.append(
                {
                    **record,
                    "matched_basenames": matched,
                    "observation_points": [observation_point],
                }
            )
    normalized.sort(key=lambda item: (item["pid"], item["pid_start_time"]))
    conflicts.sort(key=lambda item: (item["pid"], item["pid_start_time"]))
    return {
        "captured_at_utc": L.utc_now(),
        "observation_point": observation_point,
        "processes_inspected": len(normalized),
        "snapshot_sha256": canonical_sha256(normalized),
        "conflicts": conflicts,
    }


def _merge_process_conflicts(*observations: Mapping[str, Any]) -> list[dict[str, Any]]:
    merged: dict[tuple[int, str, str], dict[str, Any]] = {}
    for observation in observations:
        for raw in observation.get("conflicts", []):
            key = (int(raw["pid"]), str(raw["pid_start_time"]), str(raw["argv_sha256"]))
            if key not in merged:
                merged[key] = dict(raw)
            else:
                points = set(merged[key]["observation_points"]) | set(raw["observation_points"])
                merged[key]["observation_points"] = sorted(points)
    return sorted(merged.values(), key=lambda item: (item["pid"], item["pid_start_time"]))


def observe_conflicting_processes(
    global_lock: GlobalScopeLock,
    config: RunConfig,
    process_snapshot_provider: Optional[Callable[[], Sequence[Mapping[str, Any]]]] = None,
) -> dict[str, Any]:
    """Record a bounded /proc command observation and reject legacy/unlocked workers."""

    provider = process_snapshot_provider or _default_process_snapshot
    process_start = _classify_process_snapshot(list(provider()), "baseline_start")
    cpu_count = os.cpu_count() or 0
    memory_available = _memory_available_bytes()
    disk_free = shutil.disk_usage(global_lock.run_parent).free
    load_1m = os.getloadavg()[0]
    load_per_cpu = load_1m / cpu_count if cpu_count else float("inf")
    nfs_start = _nfs_read_bytes(config.nfs_mountpoint)
    before_total, before_iowait = _cpu_times()
    sample_started_at_utc = L.utc_now()
    sample_started_monotonic = time.monotonic()
    _bounded_sleep(config.resource_sample_seconds)
    sample_seconds = time.monotonic() - sample_started_monotonic
    sample_ended_at_utc = L.utc_now()
    after_total, after_iowait = _cpu_times()
    nfs_end = _nfs_read_bytes(config.nfs_mountpoint)
    process_end = _classify_process_snapshot(list(provider()), "baseline_end")
    matches = _merge_process_conflicts(process_start, process_end)
    if (
        nfs_start["device"] != nfs_end["device"]
        or nfs_start["mountpoint"] != nfs_end["mountpoint"]
        or nfs_start["fstype"] != nfs_end["fstype"]
    ):
        raise ContractError("E_NFS_BASELINE_UNAVAILABLE", "NFS mount identity changed during baseline")
    nfs_delta = int(nfs_end["read_bytes_recv"]) - int(nfs_start["read_bytes_recv"])
    if nfs_delta < 0 or sample_seconds <= 0:
        raise ContractError("E_NFS_BASELINE_UNAVAILABLE", "NFS READ counter reset or sample duration is zero")
    nfs_read_mbps = nfs_delta / sample_seconds / 1_000_000.0
    delta_total = after_total - before_total
    iowait_percent = (
        100.0 * (after_iowait - before_iowait) / delta_total if delta_total > 0 else 100.0
    )
    gib = 1024**3
    thresholds = {
        "min_logical_cpus": config.min_logical_cpus,
        "min_available_memory_bytes": int(config.min_available_memory_gib * gib),
        "min_free_disk_bytes": int(config.min_free_disk_gib * gib),
        "max_load_per_cpu": config.max_load_per_cpu,
        "max_iowait_percent": config.max_iowait_percent,
        "sample_seconds": config.resource_sample_seconds,
        "max_nfs_read_mbps": config.max_nfs_read_mbps,
        "nfs_mountpoint": str(config.nfs_mountpoint),
    }
    observed_resources = {
        "logical_cpus": cpu_count,
        "available_memory_bytes": memory_available,
        "free_disk_bytes": disk_free,
        "load_1m": load_1m,
        "load_per_cpu": load_per_cpu,
        "iowait_percent": iowait_percent,
        "nfs_read_baseline": {
            "mountpoint": nfs_start["mountpoint"],
            "device": nfs_start["device"],
            "fstype": nfs_start["fstype"],
            "counter_source": nfs_start["counter_source"],
            "start_read_bytes_recv": nfs_start["read_bytes_recv"],
            "end_read_bytes_recv": nfs_end["read_bytes_recv"],
            "read_bytes_recv_delta": nfs_delta,
            "sample_started_at_utc": sample_started_at_utc,
            "sample_ended_at_utc": sample_ended_at_utc,
            "sample_seconds": sample_seconds,
            "read_mbps_decimal": nfs_read_mbps,
            "threshold_max_mbps_decimal": config.max_nfs_read_mbps,
        },
    }
    resource_checks = {
        "logical_cpus": cpu_count >= thresholds["min_logical_cpus"],
        "available_memory": memory_available >= thresholds["min_available_memory_bytes"],
        "free_disk": disk_free >= thresholds["min_free_disk_bytes"],
        "load_per_cpu": load_per_cpu <= thresholds["max_load_per_cpu"],
        "iowait_percent": iowait_percent <= thresholds["max_iowait_percent"],
        "nfs_read_mbps": nfs_read_mbps < thresholds["max_nfs_read_mbps"],
    }
    resources = {
        "policy": thresholds,
        "observed": observed_resources,
        "checks": resource_checks,
        "pass": all(resource_checks.values()),
    }
    evidence = {
        "schema_name": "intersubmod.layered_process_observation",
        "schema_version": "1.0.0",
        "created_at_utc": L.utc_now(),
        "observer_pid": os.getpid(),
        "observer_pid_start_time": L._pid_start_time(os.getpid()),
        "global_scope_lock": {
            "path": str(global_lock.path),
            "held": global_lock._fd is not None,
            "run_id": global_lock.run_id,
        },
        "match_policy": sorted(CONFLICTING_PROCESS_BASENAMES),
        "processes_inspected": process_start["processes_inspected"] + process_end["processes_inspected"],
        "process_observations": {"start": process_start, "end": process_end},
        "conflicts": matches,
        "resources": resources,
        "pass": not matches and resources["pass"],
    }
    return evidence


def _bundle_source_map(bundle_manifest: Path) -> dict[str, Path]:
    document = load_json_strict(bundle_manifest)
    root = bundle_manifest.parent.resolve()
    result: dict[str, Path] = {}
    for raw in document.get("files", []):
        if not isinstance(raw, dict):
            raise ContractError("E_SOURCE_BUNDLE_MISMATCH", "malformed bundle entry")
        source_name = Path(str(raw.get("source_path", ""))).name
        bundled = (root / str(raw.get("bundled_path", ""))).resolve(strict=True)
        if source_name in result:
            raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"duplicate source basename: {source_name}")
        result[source_name] = bundled
    return result


def _require_bundle_source(bundle_manifest: Path, basename: str) -> Path:
    try:
        return _bundle_source_map(bundle_manifest)[basename]
    except KeyError as exc:
        raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"bundled source absent: {basename}") from exc


class _BundleFinder(importlib.abc.MetaPathFinder):
    def __init__(self, source_map: Mapping[str, Path]):
        self.modules = {
            Path(name).stem: path for name, path in source_map.items() if name.endswith(".py")
        }

    def find_spec(
        self,
        fullname: str,
        path: Optional[Sequence[str]] = None,
        target: Any = None,
    ) -> Optional[importlib.machinery.ModuleSpec]:
        if "." in fullname or fullname not in self.modules:
            return None
        return importlib.util.spec_from_file_location(fullname, self.modules[fullname])


def _load_module_from_bundle(bundle_manifest: Path, basename: str, module_name: str) -> Any:
    L.verify_source_bundle(bundle_manifest)
    source = _require_bundle_source(bundle_manifest, basename)
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:
        raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"cannot import bundled {basename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _exec_bundled_script(bundle_manifest: Path, basename: str, argv: Sequence[str]) -> int:
    L.verify_source_bundle(bundle_manifest)
    sources = _bundle_source_map(bundle_manifest)
    source = _require_bundle_source(bundle_manifest, basename)
    finder = _BundleFinder(sources)
    old_argv = sys.argv
    sys.meta_path.insert(0, finder)
    try:
        sys.argv = [str(source), *argv]
        try:
            runpy.run_path(str(source), run_name="__main__")
            return 0
        except SystemExit as exc:
            if exc.code is None:
                return 0
            if isinstance(exc.code, int):
                return exc.code
            print(exc.code, file=sys.stderr)
            return 1
    finally:
        sys.argv = old_argv
        with contextlib.suppress(ValueError):
            sys.meta_path.remove(finder)


def _run_bundled_validator(
    bundle_manifest: Path,
    manifest: Path,
    schema_dir: Path,
    frozen_candidate: Path,
    failure_report: Path,
) -> int:
    """Execute exact bundled validator bytes while routing bundled schemas explicitly."""

    validator_path = _require_bundle_source(bundle_manifest, "validate_layered_v3_inputs.py")
    validator = _load_module_from_bundle(
        bundle_manifest,
        "validate_layered_v3_inputs.py",
        "_layered_v3_frozen_validator",
    )
    manifest_schema = schema_dir / "layered_input_manifest_v3.schema.json"
    lock_schema = schema_dir / "layered_input_lock_v1.schema.json"
    # The validator resolves the producer-receipt schema from __file__.parents[1].
    # Point only that resource lookup at the frozen schema copy; executable bytes
    # and the recorded validator identity remain the exact core bundle path.
    validator.__file__ = str(schema_dir.parent / "scripts/validate_layered_v3_inputs.py")
    manifest_path = manifest.resolve(strict=True)
    repo_roots = [parent for parent in manifest_path.parents if parent.name == "InterSubMod"]
    if len(repo_roots) > 1:
        raise ContractError(
            "E_SOURCE_BUNDLE_MISMATCH",
            f"manifest has ambiguous InterSubMod repo ancestors: {[str(path) for path in repo_roots]}",
        )
    validator.set_method_receipt_repo_root(repo_roots[0] if repo_roots else None)
    try:
        lock = validator.validate_and_build_lock(
            manifest,
            manifest_schema,
            validator_path=validator_path,
        )
        validator.atomic_publish_frozen_lock(frozen_candidate, lock, lock_schema)
        return 0
    except validator.ContractError as exc:
        validator.atomic_write_json(
            failure_report,
            validator.failure_document(exc, manifest, frozen_candidate),
        )
        print(f"{exc.code}: {exc}", file=sys.stderr)
        return int(exc.exit_code)
    except BaseException as exc:
        internal = validator.ContractError("E_INTERNAL", repr(exc), exit_code=70, stage="internal")
        validator.atomic_write_json(
            failure_report,
            validator.failure_document(internal, manifest, frozen_candidate),
        )
        print(f"{internal.code}: {internal}", file=sys.stderr)
        return 70


def _copy_schema_bundle(bundle_manifest: Path, destination: Path) -> None:
    destination.mkdir(parents=True, mode=0o700)
    for basename in SCHEMA_BASENAMES:
        source = _require_bundle_source(bundle_manifest, basename)
        payload = source.read_bytes()
        L.atomic_write_bytes(destination / basename, payload, mode=0o444, no_overwrite=True)
        if file_sha256(destination / basename) != file_sha256(source):
            raise ContractError("E_SOURCE_BUNDLE_MISMATCH", f"schema copy drift: {basename}")


def _relocate_validator_paths(lock: dict[str, Any], staging: Path, published: Path) -> dict[str, Any]:
    relocated = json.loads(json.dumps(lock))
    validator = relocated.get("validator")
    if not isinstance(validator, dict):
        raise ContractError("E_FROZEN_LOCK_MISMATCH", "validator identity is absent from frozen lock")
    for field in ("path", "schema_path"):
        value = Path(str(validator.get(field, "")))
        try:
            relative = value.resolve(strict=False).relative_to(staging.resolve())
        except ValueError as exc:
            raise ContractError(
                "E_FROZEN_LOCK_MISMATCH",
                f"validator {field} is not under preflight staging root: {value}",
            ) from exc
        validator[field] = str(published / relative)
    digest_payload = dict(relocated)
    digest_payload.pop("lock_sha256", None)
    relocated["lock_sha256"] = canonical_sha256(digest_payload)
    return relocated


def _full_identity_ok(spec: Any, label: str) -> Path:
    if not isinstance(spec, dict) or set(spec) != {"path", "identity"}:
        raise ContractError("E_REQUIRED_ROLE", f"{label} must be a full-hash artifact")
    identity = spec["identity"]
    if not isinstance(identity, dict) or set(identity) != {"policy", "size_bytes", "sha256"}:
        raise ContractError("E_REQUIRED_ROLE", f"{label}.identity is malformed")
    if identity["policy"] != "full_sha256":
        raise ContractError("E_REQUIRED_ROLE", f"{label} is not full_sha256")
    path = Path(str(spec["path"]))
    if not path.is_file():
        raise ContractError("E_FILE_MISSING", f"{label} is missing: {path}")
    observed = (path.stat().st_size, file_sha256(path))
    expected = (identity["size_bytes"], identity["sha256"])
    if observed != expected:
        raise ContractError("E_HASH_MISMATCH", f"{label} identity changed")
    return path


def _indexed_artifact_ok(spec: Any, label: str) -> Path:
    if not isinstance(spec, dict) or set(spec) != {"path", "identity", "index"}:
        raise ContractError("E_REQUIRED_ROLE", f"{label} must be an indexed full-hash artifact")
    path = _full_identity_ok({"path": spec["path"], "identity": spec["identity"]}, label)
    _full_identity_ok(spec["index"], f"{label}.index")
    return path


def _assert_tree_vcf_pass_only(path: Path) -> None:
    try:
        import pysam
    except ImportError as exc:
        raise ContractError("E_ENVIRONMENT_MISMATCH", "pysam is required for tree VCF gate") from exc
    count = 0
    try:
        with pysam.VariantFile(str(path)) as handle:
            for record in handle:
                count += 1
                if set(record.filter.keys()) != {"PASS"}:
                    raise ContractError(
                        "E_TREE_VCF_FILTER",
                        f"tree_vcf contains non-PASS record {record.contig}:{record.pos}",
                    )
    except ContractError:
        raise
    except (OSError, ValueError) as exc:
        raise ContractError("E_TREE_VCF_FILTER", f"cannot scan tree_vcf {path}: {exc}") from exc
    if count == 0:
        raise ContractError("E_TREE_VCF_FILTER", f"tree_vcf is empty: {path}")


def _validate_copy_number_contract(copy_number: Any, sample: str) -> None:
    expected_keys = {
        "availability",
        "source",
        "semantics",
        "coordinate_system",
        "unlisted_position_semantics",
        "allowed_states",
        "overlap_policy",
        "reason",
        "cn_bed",
        "cn_int_gain",
        "cn_int_loss",
        "integration_json",
    }
    if not isinstance(copy_number, dict) or set(copy_number) != expected_keys:
        raise ContractError("E_CN_CONTRACT", f"{sample} frozen CN fields are incomplete/extra")
    availability = copy_number.get("availability")
    if availability == "measured":
        scalar_ok = (
            isinstance(copy_number.get("source"), str)
            and bool(copy_number["source"])
            and isinstance(copy_number.get("semantics"), str)
            and bool(copy_number["semantics"])
            and copy_number.get("coordinate_system") == "0_based_half_open"
            and copy_number.get("unlisted_position_semantics") == "neutral"
            and copy_number.get("allowed_states") == ["gain", "loss", "loh", "neutral"]
            and copy_number.get("overlap_policy") == "forbid"
            and copy_number.get("reason") is None
        )
        if not scalar_ok or copy_number.get("cn_bed") is None:
            raise ContractError("E_CN_CONTRACT", f"{sample} measured CN semantics/artifacts are ambiguous")
        for role in ("cn_bed", "cn_int_gain", "cn_int_loss", "integration_json"):
            if copy_number.get(role) is not None:
                _full_identity_ok(copy_number[role], f"{sample}.copy_number.{role}")
        return
    unavailable = {
        "availability": "unavailable",
        "source": "unavailable",
        "semantics": "missing; never interpreted neutral",
        "coordinate_system": None,
        "unlisted_position_semantics": "unavailable",
        "allowed_states": [],
        "overlap_policy": "not_applicable",
        "reason": "No reviewed measured CN source is available",
        "cn_bed": None,
        "cn_int_gain": None,
        "cn_int_loss": None,
        "integration_json": None,
    }
    if copy_number != unavailable:
        raise ContractError("E_CN_CONTRACT", f"{sample} unavailable CN contract is ambiguous")


def validate_execution_lock(lock: Mapping[str, Any], *, scan_tree_vcf: bool = True) -> dict[str, Any]:
    if (
        lock.get("schema_name") != "intersubmod.layered_input_lock"
        or lock.get("schema_version") != "1.0.0"
        or lock.get("all_pass") is not True
    ):
        raise ContractError("E_SCHEMA_VERSION", "validator did not emit layered input lock 1.0.0 PASS")
    samples = lock.get("samples")
    if not isinstance(samples, list) or len(samples) != 7:
        raise ContractError("E_DATASET_SET_MISMATCH", "frozen lock must contain exactly seven datasets")
    analysis = lock.get("analysis_contract")
    if not isinstance(analysis, dict):
        raise ContractError("E_SCHEMA_INVALID", "frozen analysis contract is missing")
    tree_contract = analysis.get("tree_input_contract")
    profile = TREE_PROFILES.get(tree_contract)
    if profile is None or analysis.get("task_type") != profile["task_type"]:
        raise ContractError("E_TREE_VCF_ROLE", f"unsupported frozen tree contract: {tree_contract}")
    observed: dict[str, str] = {}
    for item in samples:
        if not isinstance(item, dict) or not isinstance(item.get("sample"), str):
            raise ContractError("E_SCHEMA_INVALID", "malformed frozen sample entry")
        sample = item["sample"]
        if sample in observed:
            raise ContractError("E_SAMPLE_DUPLICATE", f"duplicate frozen sample: {sample}")
        observed[sample] = item.get("biological_id")
        if item.get("pass") is not True:
            raise ContractError("E_STATE_INVALID", f"sample lock is not PASS: {sample}")
        alignment = item.get("alignment_payload")
        if not isinstance(alignment, dict) or set(alignment) != {"path", "storage_identity_v1"}:
            raise ContractError("E_REQUIRED_ROLE", f"{sample} alignment_payload is not frozen")
        if not Path(str(alignment["path"])).is_file():
            raise ContractError("E_FILE_MISSING", f"{sample} BAM is missing")
        somatic = item.get("somatic")
        expected_roles = {
            "caller_raw_vcf",
            "longphase_input_vcf",
            "caller_pass_baseline_vcf",
            "longphase_recalibrated_all_vcf",
            "longphase_recalibrated_pass_vcf",
            "tree_vcf",
        }
        if not isinstance(somatic, dict) or set(somatic) != expected_roles:
            raise ContractError(
                "E_REQUIRED_ROLE",
                f"{sample} somatic roles must be exact {sorted(expected_roles)}; observed {sorted(somatic or {})}",
            )
        paths = {role: _indexed_artifact_ok(somatic[role], f"{sample}.{role}") for role in expected_roles}
        if paths["longphase_recalibrated_pass_vcf"].resolve() == paths["longphase_recalibrated_all_vcf"].resolve():
            raise ContractError("E_TREE_VCF_ROLE", f"{sample} LongPhase-S PASS must differ from recalibrated-all")
        selected_role = (
            "longphase_recalibrated_pass_vcf"
            if tree_contract == CANONICAL_TREE_INPUT_CONTRACT
            else "caller_pass_baseline_vcf"
        )
        if paths["tree_vcf"].resolve() != paths[selected_role].resolve():
            raise ContractError("E_TREE_VCF_ROLE", f"{sample} selected tree does not match {tree_contract}")
        if scan_tree_vcf:
            _assert_tree_vcf_pass_only(paths["tree_vcf"])
            _assert_tree_vcf_pass_only(paths["longphase_recalibrated_pass_vcf"])
        tags = item.get("read_tags")
        required_tag_roles = {
            "sidecar",
            "index",
            "validation",
            "producer_capture_receipt",
            "duplicate_identity_policy",
            "subject_binding",
            "producer_policy",
            "producer_evidence",
        }
        if not isinstance(tags, dict) or set(tags) != required_tag_roles:
            raise ContractError("E_REQUIRED_ROLE", f"{sample} read-tag lock roles are incomplete/extra")
        for role in ("sidecar", "index", "validation", "producer_capture_receipt"):
            _full_identity_ok(tags[role], f"{sample}.read_tags.{role}")
        binding = tags["subject_binding"]
        if (
            not isinstance(binding, dict)
            or binding.get("sample") != sample
            or tags.get("duplicate_identity_policy")
            != "collapse_redundant_rows_with_identical_HP_PS"
            or binding.get("duplicate_identity_policy")
            != tags.get("duplicate_identity_policy")
            or binding.get("mapped_alignment_count")
            != binding.get("identity_unique_count", 0) + binding.get("duplicate_count", -1)
            or binding.get("conflict_count") != 0
            or binding.get("longphase_recalibrated_pass_vcf_sha256")
            != somatic["longphase_recalibrated_pass_vcf"]["identity"]["sha256"]
            or binding.get("longphase_input_vcf_sha256")
            != somatic["longphase_input_vcf"]["identity"]["sha256"]
            or binding.get("caller_pass_baseline_vcf_sha256")
            != somatic["caller_pass_baseline_vcf"]["identity"]["sha256"]
        ):
            raise ContractError("E_SIDECAR_SUBJECT_MISMATCH", f"{sample} tree/input roles are not subject-bound")
        _validate_copy_number_contract(item.get("copy_number"), sample)
    if observed != EXPECTED_BINDING:
        raise ContractError(
            "E_DATASET_SET_MISMATCH",
            f"expected canonical 7/6 binding; observed {observed}",
        )
    scope = analysis.get("scope") if isinstance(analysis, dict) else None
    if (
        not isinstance(scope, dict)
        or scope.get("name") != "whole_autosomes_chr1_22"
        or scope.get("contigs") != list(AUTOSOMES)
        or analysis.get("longphase_input_contract") != "normalized_ClairS_raw_all"
        or analysis.get("tree_input_contract") not in TREE_PROFILES
    ):
        raise ContractError("E_SCOPE_NOT_COMPREHENSIVE", "frozen scope is not ordered chr1-22")
    return {item["sample"]: item for item in samples}


def _sample_from_lock(lock: Mapping[str, Any], sample: str) -> Mapping[str, Any]:
    samples = lock.get("samples")
    if not isinstance(samples, list) or len(samples) != 7:
        raise ContractError("E_DATASET_SET_MISMATCH", "runtime lock no longer contains seven datasets")
    by_id = {
        item.get("sample"): item
        for item in samples
        if isinstance(item, dict) and isinstance(item.get("sample"), str)
    }
    if set(by_id) != set(EXPECTED_BINDING):
        raise ContractError("E_DATASET_SET_MISMATCH", "runtime dataset set differs from frozen 7/6 contract")
    try:
        return by_id[sample]
    except KeyError as exc:
        raise ContractError("E_DATASET_SET_MISMATCH", f"sample is absent from lock: {sample}") from exc


def _input_set_sha256(sample_meta: Mapping[str, Any]) -> str:
    return canonical_sha256(
        {
            "sample": sample_meta["sample"],
            "biological_id": sample_meta["biological_id"],
            "alignment_payload": sample_meta["alignment_payload"],
            "somatic": sample_meta["somatic"],
            "read_tags": sample_meta["read_tags"],
            "copy_number": sample_meta["copy_number"],
        }
    )


def _read_runtime_authorities(run_root: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, str]]:
    lock_path = run_root / "frozen_input_lock.json"
    receipt_path = run_root / "launch_receipt.json"
    environment_path = run_root / "environment_lock.json"
    bundle_manifest = run_root / "source_bundle/source_bundle_manifest.json"
    lock = load_json_strict(lock_path)
    receipt = load_json_strict(receipt_path)
    bundle = L.verify_source_bundle(bundle_manifest)
    bindings = {
        "frozen_lock_sha256": file_sha256(lock_path),
        "launch_receipt_sha256": file_sha256(receipt_path),
        "environment_lock_sha256": file_sha256(environment_path),
        "source_bundle_manifest_sha256": bundle.manifest_sha256,
        "source_bundle_content_sha256": bundle.content_sha256,
    }
    for field in (
        "frozen_lock_sha256",
        "environment_lock_sha256",
        "source_bundle_manifest_sha256",
        "source_bundle_content_sha256",
    ):
        if receipt.get(field) != bindings[field]:
            raise ContractError("E_LAUNCH_RECEIPT_MISMATCH", f"receipt binding changed: {field}")
    return lock, receipt, bindings


def _atomic_json_with_adapter_fields(path: Path, additions: Mapping[str, Any]) -> dict[str, Any]:
    document = load_json_strict(path)
    overlap = set(document) & set(additions)
    if overlap:
        raise ContractError("E_OUTPUT_SCHEMA", f"adapter refuses to overwrite fields {sorted(overlap)} in {path}")
    document.update(additions)
    L.atomic_write_json(path, document)  # os.replace is atomic even when target exists.
    return document


def _run_script_process(
    bundle_manifest: Path,
    basename: str,
    argv: Sequence[str],
    env: Mapping[str, str],
    log_path: Path,
) -> None:
    runner = _require_bundle_source(bundle_manifest, "run_layered_v3.py")
    with log_path.open("ab") as log:
        completed = subprocess.run(
            [
                sys.executable,
                str(runner),
                "_exec-script",
                "--bundle-manifest",
                str(bundle_manifest),
                "--source",
                basename,
                "--",
                *argv,
            ],
            env=dict(env),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        raise ContractError(
            "E_SCIENTIFIC_WORKER",
            f"bundled {basename} exited {completed.returncode}; log={log_path}",
        )


def _output_artifact(role: str, path: Path, sample_root: Path) -> dict[str, Any]:
    try:
        relative = path.resolve(strict=True).relative_to(sample_root.resolve(strict=True))
    except (OSError, ValueError) as exc:
        raise ContractError("E_PATH_ESCAPE", f"output {role} escaped sample root: {path}") from exc
    return {
        "role": role,
        "path": str(relative),
        "size_bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def _somatic_roles(meta: Mapping[str, Any], tree_contract: str) -> dict[str, str]:
    somatic = meta["somatic"]
    profile = TREE_PROFILES[tree_contract]
    return {
        "longphase_input_role": "normalized_clairs_raw_all",
        "longphase_input_vcf_sha256": somatic["longphase_input_vcf"]["identity"]["sha256"],
        "caller_pass_baseline_role": (
            "clairs_filter_pass_selected_sensitivity_tree"
            if tree_contract == SENSITIVITY_TREE_INPUT_CONTRACT
            else "clairs_filter_pass_sensitivity_only"
        ),
        "caller_pass_baseline_vcf_sha256": somatic["caller_pass_baseline_vcf"]["identity"]["sha256"],
        "tree_input_contract": tree_contract,
        "tree_backbone_role": profile["backbone_role"],
        "tree_vcf_sha256": somatic["tree_vcf"]["identity"]["sha256"],
        "ledger_role": "clairs_raw",
        "caller_raw_vcf_sha256": somatic["caller_raw_vcf"]["identity"]["sha256"],
        "longphase_recalibrated_all_vcf_sha256": somatic["longphase_recalibrated_all_vcf"]["identity"]["sha256"],
        "longphase_recalibrated_pass_vcf_sha256": somatic[
            "longphase_recalibrated_pass_vcf"
        ]["identity"]["sha256"],
    }


def _run_sample_worker(run_root: Path, sample: str) -> int:
    lock, receipt, base_provenance = _read_runtime_authorities(run_root)
    meta = _sample_from_lock(lock, sample)
    params = receipt.get("extra", {}).get("analysis_params")
    if not isinstance(params, dict):
        raise ContractError("E_LAUNCH_RECEIPT_MISMATCH", "analysis_params absent from receipt")
    expected_profile = {
        "scope": "chr1-22",
        "contigs": list(AUTOSOMES),
        "dataset_count": 7,
        "biological_sample_count": 6,
    }
    if any(params.get(key) != value for key, value in expected_profile.items()):
        raise ContractError("E_PARAM_PROFILE", "receipt is not comprehensive chr1-22/7/6")
    if params.get("VERIFY_EVERY") != 1 or params.get("ANALYSIS_TREE_CAP") != 0:
        raise ContractError("E_PARAM_PROFILE", "receipt disabled exhaustive verification")

    bundle_manifest = run_root / "source_bundle/source_bundle_manifest.json"
    sample_root = run_root / "samples" / sample
    sample_root.mkdir(parents=True, exist_ok=False)
    somatic = meta["somatic"]
    tags = meta["read_tags"]
    copy_number = meta["copy_number"]
    bam = str(meta["alignment_payload"]["path"])
    tree_vcf = str(somatic["tree_vcf"]["path"])
    longphase_input = str(somatic["longphase_input_vcf"]["path"])
    raw_vcf = str(somatic["caller_raw_vcf"]["path"])
    recalibrated_all = str(somatic["longphase_recalibrated_all_vcf"]["path"])
    sidecar = str(tags["sidecar"]["path"])
    cn_bed = ""
    cn_gain = ""
    cn_loss = ""
    integration = ""
    cn_source = "unavailable"
    if copy_number.get("availability") == "measured":
        if not copy_number.get("source") or not copy_number.get("semantics"):
            raise ContractError("E_CN_CONTRACT", f"{sample} measured CN source/semantics are not frozen")
        cn_source = str(copy_number["source"])
        cn_bed = str(copy_number["cn_bed"]["path"]) if copy_number.get("cn_bed") else ""
        cn_gain = str(copy_number["cn_int_gain"]["path"]) if copy_number.get("cn_int_gain") else ""
        cn_loss = str(copy_number["cn_int_loss"]["path"]) if copy_number.get("cn_int_loss") else ""
        integration = str(copy_number["integration_json"]["path"]) if copy_number.get("integration_json") else ""

    provenance = {**base_provenance, "input_set_sha256": _input_set_sha256(meta)}
    tree_contract = lock["analysis_contract"]["tree_input_contract"]
    tree_profile = TREE_PROFILES[tree_contract]
    roles = _somatic_roles(meta, tree_contract)
    base_env = _controlled_child_environment()
    base_env.update(
        {
            "SM_WORKDIR": str(sample_root),
            "SM_TBAM": bam,
            "SM_VD": str(Path(tree_vcf).parent),
            "SM_SOMATIC_VCF": tree_vcf,
            "SM_CNBED": cn_bed,
            "SM_CN_SOURCE": cn_source,
            "SM_READ_TAG_SIDECAR": sidecar,
            "SM_MINREAD": str(params["MINREAD"]),
            "SM_MAX_SNV": str(params["MAX_SNV"]),
            "SM_TIER_R": str(params["TIER_R"]),
            "SM_MAPQ_MIN": str(params["MAPQ_MIN"]),
            "SM_BASEQ_MIN": str(params["BASEQ_MIN"]),
        }
    )
    outputs: list[dict[str, Any]] = []
    part_paths: list[Path] = []
    tag_binding = tags["subject_binding"]
    for part, chromosomes in PARTS.items():
        target = sample_root / f"mlhp_part_{part}.json"
        temporary = sample_root / f".mlhp_part_{part}.json.partial.{os.getpid()}"
        _run_script_process(
            bundle_manifest,
            "sm_multilocus_combinations.py",
            [",".join(chromosomes), str(temporary)],
            base_env,
            sample_root / f"mlhp_part_{part}.log",
        )
        _atomic_json_with_adapter_fields(
            temporary,
            {
                "sample": sample,
                "part": part,
                "chromosomes": list(chromosomes),
                "somatic_roles": roles,
                "provenance": provenance,
            },
        )
        part_document = load_json_strict(temporary)
        census = part_document.get("read_tag_census")
        if not isinstance(census, dict):
            raise ContractError("E_OUTPUT_SCHEMA", f"{sample} part {part} lacks read_tag_census")
        if census.get("sidecar_missing") != 0 or census.get("sidecar_conflicts") != 0:
            raise ContractError("E_JOIN_MISSING", f"{sample} part {part} sidecar join is not exact")
        if census.get("alignment_identity_allele_conflicts") != 0:
            raise ContractError(
                "E_JOIN_CONFLICT", f"{sample} part {part} has conflicting alleles for one alignment identity"
            )
        if census.get("sidecar_exact_matches") != census.get("alignment_group_exposures"):
            raise ContractError("E_JOIN_MULTIPLICITY", f"{sample} part {part} exposure mismatch")
        census_additions = {
            "identity_schema": "coordinate_join_v1",
            "sidecar_sha256": tags["sidecar"]["identity"]["sha256"],
            "sidecar_index_sha256": tags["index"]["identity"]["sha256"],
            "alignment_payload_identity_sha256": meta["alignment_payload"]["storage_identity_v1"]["identity_sha256"],
            "sidecar_duplicates": int(tag_binding["duplicate_count"]),
            "duplicate_identity_policy": tags["duplicate_identity_policy"],
            "sidecar_extra": 0,
            "sidecar_malformed": 0,
        }
        incompatible = {
            key: census[key]
            for key, expected in census_additions.items()
            if key in census and census[key] != expected
        }
        if incompatible:
            raise ContractError(
                "E_JOIN_MULTIPLICITY",
                f"{sample} part {part} emitted diagnostics conflicting with frozen evidence: {incompatible}",
            )
        census.update(census_additions)
        if any(census[name] != 0 for name in ("sidecar_extra", "sidecar_malformed")):
            raise ContractError("E_JOIN_MULTIPLICITY", f"{sample} frozen sidecar census is not clean")
        L.atomic_write_json(temporary, part_document)
        os.replace(temporary, target)
        part_paths.append(target)
        outputs.append(_output_artifact(f"mlhp_part_{part}", target, sample_root))

    layered = sample_root / f"layered_reconstruction_{sample}.json"
    layered_temp = sample_root / f".{layered.name}.partial.{os.getpid()}"
    layered_env = dict(base_env)
    layered_env.update(
        {
            "SM_ML_GLOB": str(sample_root / "mlhp_part_*.json"),
            "SM_OUT": str(layered_temp),
            "SM_VERIFY_EVERY": str(params["VERIFY_EVERY"]),
            "SM_ANALYSIS_TREE_CAP": str(params["ANALYSIS_TREE_CAP"]),
            "SM_DISPLAY_TREE_CAP": str(params["DISPLAY_TREE_CAP"]),
            "SM_CN_INT_GAIN": cn_gain,
            "SM_CN_INT_LOSS": cn_loss,
        }
    )
    _run_script_process(
        bundle_manifest,
        "layered_tree_reconstruction.py",
        [],
        layered_env,
        sample_root / "layered.log",
    )
    _atomic_json_with_adapter_fields(
        layered_temp,
        {"sample": sample, "provenance": provenance},
    )
    os.replace(layered_temp, layered)
    outputs.append(_output_artifact("layered_reconstruction", layered, sample_root))

    region = sample_root / f"layered_region_view_{sample}.json"
    region_temp = sample_root / f".{region.name}.partial.{os.getpid()}"
    region_env = dict(base_env)
    region_env.update(
        {
            "SM_LAYERED": str(layered),
            "SM_OUT": str(region_temp),
            "SM_SAMPLE": sample,
            "SM_ML_GLOB": str(sample_root / "mlhp_part_*.json"),
            "SM_SOMATIC_VCF": tree_vcf,
            "SM_INTEGRATION": integration,
            "SM_BACKBONE_SOURCE": tree_profile["backbone_role"],
        }
    )
    _run_script_process(
        bundle_manifest,
        "build_region_view.py",
        [],
        region_env,
        sample_root / "region_view.log",
    )
    copy_number_contract = json.loads(json.dumps(copy_number, sort_keys=True))
    _atomic_json_with_adapter_fields(
        region_temp,
        {"provenance": provenance, "copy_number_contract": copy_number_contract},
    )
    os.replace(region_temp, region)
    outputs.append(_output_artifact("layered_region_view", region, sample_root))

    ledger = sample_root / f"ssnv_site_ledger_{sample}.tsv.gz"
    ledger_summary = sample_root / f"ssnv_site_ledger_{sample}.summary.json"
    ledger_stage = sample_root / f".site_ledger.stage.{os.getpid()}"
    ledger_stage.mkdir(mode=0o700)
    staged_ledger = ledger_stage / ledger.name
    staged_summary = ledger_stage / ledger_summary.name
    _run_script_process(
        bundle_manifest,
        "build_ssnv_site_ledger.py",
        [
            "--sample",
            sample,
            "--caller-raw-vcf",
            raw_vcf,
            "--longphase-input-vcf",
            longphase_input,
            "--tree-input-vcf",
            tree_vcf,
            "--tree-contract",
            tree_profile["ledger_tree_contract"],
            "--longphase-input-contract",
            "clairs_raw_all",
            "--recalibrated-vcf",
            recalibrated_all,
            "--mlhp-glob",
            str(sample_root / "mlhp_part_*.json"),
            "--output-tsv-gz",
            str(staged_ledger),
            "--output-summary",
            str(staged_summary),
        ],
        base_env,
        sample_root / "site_ledger.log",
    )
    staged_index = Path(f"{staged_ledger}.tbi")
    if not staged_ledger.is_file() or not staged_index.is_file() or not staged_summary.is_file():
        raise ContractError("E_REQUIRED_OUTPUT", f"{sample} site ledger artifact set is incomplete")
    _atomic_json_with_adapter_fields(staged_summary, {"provenance": provenance})
    os.replace(staged_ledger, ledger)
    os.replace(staged_index, Path(f"{ledger}.tbi"))
    os.replace(staged_summary, ledger_summary)
    outputs.append(_output_artifact("site_ledger", ledger, sample_root))
    outputs.append(_output_artifact("site_ledger_summary", ledger_summary, sample_root))

    output_manifest = {
        "schema_name": "intersubmod.layered_sample_output_manifest",
        "schema_version": "1.0.0",
        "sample": sample,
        "biological_id": meta["biological_id"],
        "run_id": run_root.name,
        **provenance,
        "somatic_roles": roles,
        "copy_number_contract": copy_number_contract,
        "outputs": outputs,
    }
    L.atomic_write_json(sample_root / "output_manifest.json", output_manifest, no_overwrite=True)
    return 0


def _launch_full_run(
    lifecycle: RunLifecycle,
    config: RunConfig,
    bundle_manifest: Path,
) -> Path:
    lifecycle.begin_running()
    lifecycle.start_heartbeat(interval_seconds=config.heartbeat_interval)
    runner = _require_bundle_source(bundle_manifest, "run_layered_v3.py")
    samples = sorted(EXPECTED_BINDING)
    child_env = _controlled_child_environment()
    for offset in range(0, len(samples), config.parallel_samples):
        batch = samples[offset : offset + config.parallel_samples]
        lifecycle.set_active_samples(batch)
        handles = []
        try:
            for sample in batch:
                log_path = lifecycle.root / f"worker_{sample}.log"
                handle = log_path.open("ab")
                handles.append(handle)
                lifecycle.children.launch(
                    [
                        sys.executable,
                        str(runner),
                        "_worker",
                        "--run-root",
                        str(lifecycle.root),
                        "--sample",
                        sample,
                    ],
                    label=f"worker:{sample}",
                    env=child_env,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                )
            lifecycle.children.wait_all_fail_fast(poll_seconds=0.05, grace_seconds=2.0)
        finally:
            for handle in handles:
                handle.close()
    lifecycle.set_active_samples([])

    manifest_index = []
    for sample in samples:
        path = lifecycle.root / "samples" / sample / "output_manifest.json"
        if not path.is_file():
            raise ContractError("E_REQUIRED_OUTPUT", f"missing sample output manifest: {sample}")
        manifest_index.append(
            {"sample": sample, "path": str(path.relative_to(lifecycle.root)), "sha256": file_sha256(path)}
        )
    L.atomic_write_json(
        lifecycle.root / "output_manifests.json",
        {"dataset_count": 7, "manifests": manifest_index},
        no_overwrite=True,
    )

    lifecycle.begin_verifying()
    lifecycle.set_active_samples(["VERIFY"])
    verifier = _require_bundle_source(bundle_manifest, "verify_layered_v3.py")
    verification = lifecycle.root / "verification_summary.json"
    verify_log = (lifecycle.root / "verifier.log").open("ab")
    try:
        lifecycle.children.launch(
            [
                sys.executable,
                str(verifier),
                "--run-root",
                str(lifecycle.root),
                "--frozen-lock",
                str(lifecycle.root / "frozen_input_lock.json"),
                "--launch-receipt",
                str(lifecycle.root / "launch_receipt.json"),
                "--output",
                str(verification),
            ],
            label="verifier",
            env=child_env,
            stdout=verify_log,
            stderr=subprocess.STDOUT,
        )
        lifecycle.children.wait_all_fail_fast(poll_seconds=0.05, grace_seconds=2.0)
    finally:
        verify_log.close()
    summary = load_json_strict(verification)
    if summary.get("all_pass") is not True or summary.get("dataset_count") != 7:
        raise ContractError("E_VERIFICATION_FAILED", "bundled verifier did not return exact 7/7 PASS")
    lifecycle.set_active_samples([])
    return lifecycle.succeed(
        verification,
        success_extra={
            "mode": "full",
            "dataset_count": 7,
            "biological_sample_count": 6,
            "scope": "chr1-22",
        },
    )


def run_pipeline(
    config: RunConfig,
    inventory: Optional[SourceInventory] = None,
    *,
    after_bundle_hook: Optional[Callable[[Path], None]] = None,
    process_snapshot_provider: Optional[Callable[[], Sequence[Mapping[str, Any]]]] = None,
) -> Path:
    """Hold the fixed full-scope lock across preflight, runtime and terminal state."""

    with GlobalScopeLock(config.run_parent, config.run_id) as global_lock:
        return _run_pipeline_locked(
            config,
            inventory,
            global_lock=global_lock,
            after_bundle_hook=after_bundle_hook,
            process_snapshot_provider=process_snapshot_provider,
        )


def _run_pipeline_locked(
    config: RunConfig,
    inventory: Optional[SourceInventory] = None,
    *,
    global_lock: GlobalScopeLock,
    after_bundle_hook: Optional[Callable[[Path], None]] = None,
    process_snapshot_provider: Optional[Callable[[], Sequence[Mapping[str, Any]]]] = None,
) -> Path:
    """Run preflight or the complete frozen v3 pipeline.

    ``after_bundle_hook`` exists solely for deterministic mutation-injection
    tests.  The production CLI never exposes it.
    """

    _validate_config(config)
    environment = _deterministic_environment()
    inventory = inventory or production_inventory()
    _validate_source_inventory(inventory)
    _validate_production_resource_policy(config, inventory)
    with RunLifecycle(config.run_parent, config.run_id) as lifecycle:
        lifecycle.begin_preflight()
        process_observation_path = lifecycle.root / "process_observation.json"
        try:
            process_observation = observe_conflicting_processes(
                global_lock,
                config,
                process_snapshot_provider=process_snapshot_provider,
            )
        except ContractError as exc:
            process_observation = {
                "schema_name": "intersubmod.layered_process_observation",
                "schema_version": "1.0.0",
                "created_at_utc": L.utc_now(),
                "observer_pid": os.getpid(),
                "observer_pid_start_time": L._pid_start_time(os.getpid()),
                "global_scope_lock": {
                    "path": str(global_lock.path),
                    "held": global_lock._fd is not None,
                    "run_id": global_lock.run_id,
                },
                "match_policy": sorted(CONFLICTING_PROCESS_BASENAMES),
                "processes_inspected": 0,
                "process_observations": None,
                "conflicts": [],
                "resources": {
                    "policy": config.params()["resource_policy"],
                    "observed": None,
                    "checks": {"evidence_available": False},
                    "pass": False,
                    "error": {"code": exc.code, "message": exc.message},
                },
                "pass": False,
            }
            L.atomic_write_json(
                process_observation_path,
                process_observation,
                no_overwrite=True,
            )
            evidence_sha = file_sha256(process_observation_path)
            raise ContractError(
                exc.code,
                f"resource/process evidence unavailable; evidence={process_observation_path} sha256={evidence_sha}; {exc.message}",
            ) from exc
        L.atomic_write_json(
            process_observation_path,
            process_observation,
            no_overwrite=True,
        )
        process_observation_sha256 = file_sha256(process_observation_path)
        if process_observation["conflicts"]:
            raise ContractError(
                "E_CONFLICTING_FULL_RUN",
                f"unlocked/legacy layered workers are active; evidence={process_observation_path}",
            )
        if not process_observation["resources"]["pass"]:
            failed_checks = sorted(
                name
                for name, passed in process_observation["resources"]["checks"].items()
                if not passed
            )
            raise ContractError(
                "E_RESOURCE_GATE",
                f"host resource gate failed {failed_checks}; evidence={process_observation_path}",
            )
        bundle = lifecycle.build_source_bundle(
            inventory.runner,
            inventory.validator,
            inventory.verifier,
            inventory.imported,
        )
        if after_bundle_hook is not None:
            after_bundle_hook(bundle.manifest_path)
        sm_defaults = {
            "SM_INPUT_MANIFEST": str(config.manifest.resolve()),
            "SM_RUN_ID": config.run_id,
            "SM_RUN_PARENT": str(config.run_parent.resolve()),
            "SM_PARALLEL_SAMPLES": str(config.parallel_samples),
            "SM_VERIFY_EVERY": str(config.verify_every),
            "SM_ANALYSIS_TREE_CAP": str(config.analysis_tree_cap),
            "SM_DISPLAY_TREE_CAP": str(config.display_tree_cap),
            "SM_MINREAD": str(config.minread),
            "SM_MAX_SNV": str(config.max_snv),
            "SM_TIER_R": str(config.tier_r),
            "SM_MAPQ_MIN": str(config.mapq_min),
            "SM_BASEQ_MIN": str(config.baseq_min),
            "SM_HEARTBEAT_INTERVAL": str(config.heartbeat_interval),
            "SM_MIN_LOGICAL_CPUS": str(config.min_logical_cpus),
            "SM_MIN_AVAILABLE_MEMORY_GIB": str(config.min_available_memory_gib),
            "SM_MIN_FREE_DISK_GIB": str(config.min_free_disk_gib),
            "SM_MAX_LOAD_PER_CPU": str(config.max_load_per_cpu),
            "SM_MAX_IOWAIT_PERCENT": str(config.max_iowait_percent),
            "SM_RESOURCE_SAMPLE_SECONDS": str(config.resource_sample_seconds),
            "SM_MAX_NFS_READ_MBPS": str(config.max_nfs_read_mbps),
            "SM_NFS_MOUNTPOINT": str(config.nfs_mountpoint),
        }
        capture_environment = dict(environment)
        # CLI/environment values have already been resolved into ``config``.
        # Remove recognised SM_* names so the lock records their exact effective
        # values from sm_defaults, while leaving unknown SM_* names visible to
        # the lifecycle's fail-closed rejection.
        for name in SM_ALLOWLIST:
            capture_environment.pop(name, None)
        lifecycle.capture_environment_lock(
            sm_allowlist=SM_ALLOWLIST,
            sm_defaults=sm_defaults,
            environment=capture_environment,
            distributions=inventory.distributions,
            tools=inventory.tools,
            storage_path=config.run_parent,
        )
        schema_dir = lifecycle.root / "validator_runtime/schemas"
        _copy_schema_bundle(bundle.manifest_path, schema_dir)
        candidate = lifecycle.root / "validator_runtime/validated_input_lock.json"
        failure = lifecycle.root / "validator_runtime/validation_failure.json"
        runner = bundle.role_paths["runner"]
        validator_log = lifecycle.root / "validator_runtime/validator.log"
        with validator_log.open("ab") as log:
            completed = subprocess.run(
                [
                    sys.executable,
                    str(runner),
                    "_run-validator",
                    "--bundle-manifest",
                    str(bundle.manifest_path),
                    "--manifest",
                    str(config.manifest.resolve()),
                    "--schema-dir",
                    str(schema_dir),
                    "--frozen-candidate",
                    str(candidate),
                    "--failure-report",
                    str(failure),
                ],
                env=_controlled_child_environment(environment),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode != 0 or not candidate.is_file():
            code = "E_VALIDATOR_FAILED"
            if failure.is_file():
                failure_document = load_json_strict(failure)
                code = str(failure_document.get("error", {}).get("code") or code)
            raise ContractError(code, f"bundled validator rejected manifest; report={failure}")
        lock = _relocate_validator_paths(
            load_json_strict(candidate),
            lifecycle.staging_root,
            lifecycle.published_root,
        )
        validate_execution_lock(lock, scan_tree_vcf=True)
        source_sha = lock.get("source_manifest", {}).get("byte_sha256")
        if source_sha != file_sha256(config.manifest.resolve()):
            raise ContractError("E_SOURCE_MANIFEST_MISMATCH", "manifest changed after bundled validation")
        manifest_snapshot = lifecycle.root / "input_manifest.snapshot.json"
        L.atomic_write_bytes(
            manifest_snapshot,
            config.manifest.resolve().read_bytes(),
            mode=0o444,
            no_overwrite=True,
        )
        if file_sha256(manifest_snapshot) != source_sha:
            raise ContractError("E_SOURCE_MANIFEST_MISMATCH", "manifest snapshot differs from validator input")
        lifecycle.write_frozen_lock(lock)
        tree_contract = lock["analysis_contract"]["tree_input_contract"]
        tree_profile = TREE_PROFILES[tree_contract]
        mode = "preflight_only" if config.preflight_only else "full"
        lifecycle.publish_ready(
            receipt_extra={
                "mode": mode,
                "task_type": tree_profile["task_type"],
                "analysis_params": config.params(),
                "source_manifest_snapshot": {
                    "path": "input_manifest.snapshot.json",
                    "sha256": source_sha,
                },
                "process_observation": {
                    "path": "process_observation.json",
                    "sha256": process_observation_sha256,
                    "conflict_count": 0,
                },
                "validator_execution": "bundled_core_via_runner_adapter",
                "worker_input_authority": "frozen_input_lock.json",
                "tree_backbone_role": tree_profile["backbone_role"],
                "ledger_roles": [
                    "caller_raw_vcf",
                    "longphase_input_vcf",
                    "caller_pass_baseline_vcf",
                    "longphase_recalibrated_all_vcf",
                    "longphase_recalibrated_pass_vcf",
                    "tree_vcf",
                ],
            }
        )
        published_bundle = lifecycle.root / "source_bundle/source_bundle_manifest.json"
        L.verify_source_bundle(
            published_bundle,
            expected_manifest_sha256=lifecycle.source_bundle_manifest_sha256,
            expected_content_sha256=lifecycle.source_bundle_content_sha256,
        )
        if config.preflight_only:
            lifecycle.begin_running()
            lifecycle.begin_verifying()
            report = lifecycle.root / "preflight_verification.json"
            L.atomic_write_json(
                report,
                {
                    "schema_name": "intersubmod.layered_preflight_verification",
                    "schema_version": "1.0.0",
                    "mode": "preflight_only",
                    "pass": True,
                    "scientific_workers_started": 0,
                    "dataset_count": 7,
                    "biological_sample_count": 6,
                    "scope": list(AUTOSOMES),
                },
                no_overwrite=True,
            )
            lifecycle.succeed(
                report,
                success_extra={
                    "mode": "preflight_only",
                    "scientific_workers_started": 0,
                    "production_result": False,
                },
            )
            return lifecycle.root
        _launch_full_run(lifecycle, config, published_bundle)
        return lifecycle.root


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--run-parent", type=Path)
    parser.add_argument("--run-id")
    parser.add_argument("--parallel-samples", type=int)
    parser.add_argument("--verify-every", type=int)
    parser.add_argument("--analysis-tree-cap", type=int)
    parser.add_argument("--display-tree-cap", type=int)
    parser.add_argument("--minread", type=int)
    parser.add_argument("--max-snv", type=int)
    parser.add_argument("--tier-r", type=int)
    parser.add_argument("--mapq-min", type=int)
    parser.add_argument("--baseq-min", type=int)
    parser.add_argument("--heartbeat-interval", type=float)
    parser.add_argument("--min-logical-cpus", type=int)
    parser.add_argument("--min-available-memory-gib", type=float)
    parser.add_argument("--min-free-disk-gib", type=float)
    parser.add_argument("--max-load-per-cpu", type=float)
    parser.add_argument("--max-iowait-percent", type=float)
    parser.add_argument("--resource-sample-seconds", type=float)
    parser.add_argument("--max-nfs-read-mbps", type=float)
    parser.add_argument("--nfs-mountpoint", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    subparsers = parser.add_subparsers(dest="internal_command")
    validator = subparsers.add_parser("_run-validator", help=argparse.SUPPRESS)
    validator.add_argument("--bundle-manifest", required=True, type=Path)
    validator.add_argument("--manifest", required=True, type=Path)
    validator.add_argument("--schema-dir", required=True, type=Path)
    validator.add_argument("--frozen-candidate", required=True, type=Path)
    validator.add_argument("--failure-report", required=True, type=Path)
    executor = subparsers.add_parser("_exec-script", help=argparse.SUPPRESS)
    executor.add_argument("--bundle-manifest", required=True, type=Path)
    executor.add_argument("--source", required=True)
    executor.add_argument("argv", nargs=argparse.REMAINDER)
    worker = subparsers.add_parser("_worker", help=argparse.SUPPRESS)
    worker.add_argument("--run-root", required=True, type=Path)
    worker.add_argument("--sample", required=True)
    return parser


def _argument_or_env(value: Any, name: str, default: Any, converter: Callable[[str], Any]) -> Any:
    if value is not None:
        return value
    raw = os.environ.get(name)
    return default if raw is None else converter(raw)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.internal_command == "_run-validator":
            return _run_bundled_validator(
                args.bundle_manifest,
                args.manifest,
                args.schema_dir,
                args.frozen_candidate,
                args.failure_report,
            )
        if args.internal_command == "_exec-script":
            script_argv = list(args.argv)
            if script_argv and script_argv[0] == "--":
                script_argv = script_argv[1:]
            return _exec_bundled_script(args.bundle_manifest, args.source, script_argv)
        if args.internal_command == "_worker":
            return _run_sample_worker(args.run_root.resolve(strict=True), args.sample)

        manifest = _argument_or_env(
            args.manifest,
            "SM_INPUT_MANIFEST",
            METHOD_ROOT / "data/layered_input_manifest_v3.json",
            Path,
        )
        run_parent = _argument_or_env(
            args.run_parent,
            "SM_RUN_PARENT",
            Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds"),
            Path,
        )
        run_id = _argument_or_env(
            args.run_id,
            "SM_RUN_ID",
            time.strftime("%Y%m%d_%H%M%S_layered_reconstruction_v3"),
            str,
        )
        config = RunConfig(
            manifest=Path(manifest),
            run_parent=Path(run_parent),
            run_id=str(run_id),
            parallel_samples=_argument_or_env(args.parallel_samples, "SM_PARALLEL_SAMPLES", 2, int),
            verify_every=_argument_or_env(args.verify_every, "SM_VERIFY_EVERY", 1, int),
            analysis_tree_cap=_argument_or_env(args.analysis_tree_cap, "SM_ANALYSIS_TREE_CAP", 0, int),
            display_tree_cap=_argument_or_env(args.display_tree_cap, "SM_DISPLAY_TREE_CAP", 32, int),
            minread=_argument_or_env(args.minread, "SM_MINREAD", 3, int),
            max_snv=_argument_or_env(args.max_snv, "SM_MAX_SNV", 8, int),
            tier_r=_argument_or_env(args.tier_r, "SM_TIER_R", 50000, int),
            mapq_min=_argument_or_env(args.mapq_min, "SM_MAPQ_MIN", 20, int),
            baseq_min=_argument_or_env(args.baseq_min, "SM_BASEQ_MIN", 0, int),
            heartbeat_interval=_argument_or_env(
                args.heartbeat_interval,
                "SM_HEARTBEAT_INTERVAL",
                60.0,
                float,
            ),
            min_logical_cpus=_argument_or_env(
                args.min_logical_cpus, "SM_MIN_LOGICAL_CPUS", 8, int
            ),
            min_available_memory_gib=_argument_or_env(
                args.min_available_memory_gib, "SM_MIN_AVAILABLE_MEMORY_GIB", 128.0, float
            ),
            min_free_disk_gib=_argument_or_env(
                args.min_free_disk_gib, "SM_MIN_FREE_DISK_GIB", 500.0, float
            ),
            max_load_per_cpu=_argument_or_env(
                args.max_load_per_cpu, "SM_MAX_LOAD_PER_CPU", 1.25, float
            ),
            max_iowait_percent=_argument_or_env(
                args.max_iowait_percent, "SM_MAX_IOWAIT_PERCENT", 20.0, float
            ),
            resource_sample_seconds=_argument_or_env(
                args.resource_sample_seconds, "SM_RESOURCE_SAMPLE_SECONDS", 300.0, float
            ),
            max_nfs_read_mbps=_argument_or_env(
                args.max_nfs_read_mbps, "SM_MAX_NFS_READ_MBPS", 80.0, float
            ),
            nfs_mountpoint=Path(
                _argument_or_env(args.nfs_mountpoint, "SM_NFS_MOUNTPOINT", Path("/big8_disk"), Path)
            ),
            preflight_only=args.preflight_only,
        )
        root = run_pipeline(config)
        print(f"LAYERED V3 {'PREFLIGHT' if config.preflight_only else 'FULL'} PASS -> {root}")
        return 0
    except ContractError as exc:
        print(json.dumps({"pass": False, "error_code": exc.code, "message": exc.message}), file=sys.stderr)
        return 7
    except BaseException as exc:
        print(json.dumps({"pass": False, "error_code": "E_INTERNAL", "message": repr(exc)}), file=sys.stderr)
        return 70


if __name__ == "__main__":
    raise SystemExit(main())
