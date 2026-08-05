#!/usr/bin/env python3
"""Fail-closed one-shot handoff from LongPhase production to layered-v3.

The default mode is a read-only dry run.  Production execution additionally
requires both ``--execute-reviewed-plan`` and an out-of-band SHA-256 of the
reviewed launch-plan bytes.  No command is retried and no output is
overwritten.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
METHOD_ROOT = REPO_ROOT / "docs/methodology/_assets/20260627_subclone_4axis_teaching"
METHOD_SCRIPTS = METHOD_ROOT / "scripts"
METHOD_SCHEMAS = METHOD_ROOT / "schemas"
PROBE_ROOT = (
    REPO_ROOT
    / "research/20260710_layered_reconstruction_v2/probes/"
    "20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1"
)

EXPECTED_PRODUCTION_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260711_longphase_s_production_sidecars_v1"
)
EXPECTED_RUN_PARENT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds"
)
EXPECTED_REVIEWED_PRODUCER_MANIFEST = (
    METHOD_ROOT / "data/longphase_production_sidecar_manifest.json"
)
EXPECTED_BASE_MANIFEST = METHOD_ROOT / "data/layered_v2_input_manifest.json"
EXPECTED_REAL_RECEIPT = PROBE_ROOT / "equivalence_probe.json"
EXPECTED_SYNTHETIC_RECEIPT = PROBE_ROOT / "synthetic_contract_receipt.json"
EXPECTED_LONGPHASE_BINARY = Path(
    "/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s"
)
EXPECTED_LONGPHASE_SHA256 = (
    "07cbd0aa0c9f33ed59c5baff45fbe3554ef96d55457635de7348c4501b283f54"
)
EXPECTED_PS_BINARY = Path("/usr/bin/ps")

EXPECTED_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)

REQUIRED_INPUT_ROLES = frozenset(
    {
        "producer_manifest_reviewed",
        "producer_manifest_snapshot",
        "base_manifest",
        "real_method_receipt",
        "synthetic_method_receipt",
    }
)
REQUIRED_EXECUTABLE_ROLES = frozenset(
    {
        "python",
        "handoff_supervisor",
        "longphase_binary",
        "producer_finalizer",
        "receipt_normalizer",
        "v3_preparer",
        "v3_runner",
        "v3_lifecycle",
        "v3_validator",
        "v3_verifier",
        "schema_manifest",
        "schema_lock",
        "schema_receipt",
        "science_sm_linkage_genomewide",
        "science_sm_multilocus_combinations",
        "science_tree_enumeration_solver",
        "science_layered_tree_reconstruction",
        "science_build_region_view",
        "science_build_ssnv_site_ledger",
        "tool_stat",
        "tool_samtools",
        "tool_bcftools",
        "tool_bgzip",
        "tool_tabix",
        "tool_ps",
    }
)

RUNNER_OPTIONS = (
    "--parallel-samples", "2",
    "--verify-every", "1",
    "--analysis-tree-cap", "0",
    "--display-tree-cap", "32",
    "--minread", "3",
    "--max-snv", "8",
    "--tier-r", "50000",
    "--mapq-min", "20",
    "--baseq-min", "0",
    "--heartbeat-interval", "60",
    "--min-logical-cpus", "8",
    "--min-available-memory-gib", "128",
    "--min-free-disk-gib", "500",
    "--max-load-per-cpu", "1.25",
    "--max-iowait-percent", "20",
    "--resource-sample-seconds", "300",
    "--max-nfs-read-mbps", "80",
    "--nfs-mountpoint", "/big8_disk",
)

SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")


class HandoffError(RuntimeError):
    def __init__(self, code: str, message: str, *, exit_code: int = 7):
        super().__init__(message)
        self.code = code
        self.message = message
        self.exit_code = exit_code


class ProducerPending(HandoffError):
    def __init__(self, message: str):
        super().__init__("E_PRODUCER_PENDING", message, exit_code=7)


def fail(code: str, message: str, *, exit_code: int = 7) -> None:
    raise HandoffError(code, message, exit_code=exit_code)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def strict_json_load(path: Path) -> Any:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                fail("E_JSON_DUPLICATE_KEY", f"duplicate JSON key {key!r}: {path}")
            result[key] = value
        return result

    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=no_duplicates)
    except HandoffError:
        raise
    except Exception as exc:
        fail("E_JSON_INVALID", f"cannot parse {path}: {exc}")


def atomic_write_json(path: Path, value: Any, *, mode: int = 0o444) -> None:
    if path.exists() or path.is_symlink():
        fail("E_OUTPUT_EXISTS", f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode()
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError:
        fail("E_OUTPUT_EXISTS", f"refusing to overwrite {path}")
    with os.fdopen(fd, "wb") as handle:
        handle.write(data)
        handle.flush()
        os.fchmod(handle.fileno(), mode)
        os.fsync(handle.fileno())
    directory_fd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def exact_path(raw: Any, label: str) -> Path:
    if not isinstance(raw, str) or not raw:
        fail("E_PLAN_SCHEMA", f"{label} must be a non-empty absolute path")
    path = Path(raw)
    if not path.is_absolute():
        fail("E_PLAN_SCHEMA", f"{label} must be absolute: {path}")
    return path


@dataclass(frozen=True)
class ArtifactPin:
    role: str
    path: Path
    resolved_path: Path
    sha256: str

    @classmethod
    def from_json(cls, role: str, value: Any) -> "ArtifactPin":
        if not isinstance(value, dict) or set(value) != {"path", "resolved_path", "sha256"}:
            fail("E_PLAN_SCHEMA", f"pin {role} requires path/resolved_path/sha256 only")
        path = exact_path(value["path"], f"pins.{role}.path")
        resolved = exact_path(value["resolved_path"], f"pins.{role}.resolved_path")
        digest = value["sha256"]
        if not isinstance(digest, str) or not HEX64.fullmatch(digest):
            fail("E_PLAN_SCHEMA", f"pin {role} has invalid sha256")
        return cls(role, path, resolved, digest)

    def verify(self) -> None:
        if not self.path.exists() or not self.path.is_file():
            fail("E_PIN_MISSING", f"pinned artifact missing: {self.role}={self.path}")
        observed_resolved = self.path.resolve(strict=True)
        if observed_resolved != self.resolved_path:
            fail(
                "E_PIN_PATH_DRIFT",
                f"{self.role} resolved path drift: {observed_resolved} != {self.resolved_path}",
            )
        observed = sha256_file(observed_resolved)
        if observed != self.sha256:
            fail("E_PIN_HASH_DRIFT", f"{self.role} hash drift: {observed} != {self.sha256}")


@dataclass(frozen=True)
class LaunchPlan:
    path: Path
    byte_sha256: str
    plan_id: str
    execute_authorized: bool
    reviewed_by: str | None
    reviewed_at_utc: str | None
    production_root: Path
    handoff_lock: Path
    handoff_workspace: Path
    source_manifest: Path
    source_manifest_failure: Path
    manifest_id: str
    run_parent: Path
    run_id: str
    inputs: Mapping[str, ArtifactPin]
    executables: Mapping[str, ArtifactPin]

    @property
    def run_root(self) -> Path:
        return self.run_parent / self.run_id

    @property
    def all_pins(self) -> tuple[ArtifactPin, ...]:
        return tuple(self.inputs.values()) + tuple(self.executables.values())


def _pin_map(raw: Any, required: frozenset[str], label: str) -> dict[str, ArtifactPin]:
    if not isinstance(raw, dict) or set(raw) != set(required):
        missing = sorted(set(required) - set(raw or {})) if isinstance(raw, dict) else sorted(required)
        extra = sorted(set(raw or {}) - set(required)) if isinstance(raw, dict) else []
        fail("E_PLAN_SCHEMA", f"{label} roles differ; missing={missing} extra={extra}")
    return {role: ArtifactPin.from_json(role, raw[role]) for role in sorted(required)}


def load_launch_plan(path: Path) -> LaunchPlan:
    if not path.is_file() or path.is_symlink():
        fail("E_PLAN_MISSING", f"launch plan must be a regular non-symlink file: {path}")
    byte_sha = sha256_file(path)
    raw = strict_json_load(path)
    required_top = {
        "schema_name", "schema_version", "plan_id", "execution_authorization",
        "production_root", "handoff", "layered_run", "inputs", "executables",
    }
    if not isinstance(raw, dict) or set(raw) != required_top:
        fail("E_PLAN_SCHEMA", "launch plan top-level fields differ from the frozen contract")
    if raw["schema_name"] != "intersubmod.layered_v3_handoff_launch_plan" or raw["schema_version"] != "1.0.0":
        fail("E_PLAN_SCHEMA", "unsupported handoff launch-plan schema")
    plan_id = raw["plan_id"]
    if not isinstance(plan_id, str) or not SAFE_ID.fullmatch(plan_id):
        fail("E_PLAN_SCHEMA", "plan_id is not a safe identifier")
    auth = raw["execution_authorization"]
    if not isinstance(auth, dict) or set(auth) != {"execute", "reviewed_by", "reviewed_at_utc"}:
        fail("E_PLAN_SCHEMA", "execution_authorization fields differ")
    if not isinstance(auth["execute"], bool):
        fail("E_PLAN_SCHEMA", "execution_authorization.execute must be boolean")
    review_values = (auth["reviewed_by"], auth["reviewed_at_utc"])
    if auth["execute"]:
        if not all(isinstance(value, str) and value for value in review_values):
            fail("E_PLAN_SCHEMA", "authorized execution requires non-empty review metadata")
    elif review_values != (None, None):
        fail("E_PLAN_SCHEMA", "unreviewed execute=false plan requires null review metadata")
    production_root = exact_path(raw["production_root"], "production_root")
    handoff = raw["handoff"]
    if not isinstance(handoff, dict) or set(handoff) != {
        "lock", "workspace", "source_manifest", "source_manifest_failure", "manifest_id"
    }:
        fail("E_PLAN_SCHEMA", "handoff fields differ")
    layered = raw["layered_run"]
    if not isinstance(layered, dict) or set(layered) != {"run_parent", "run_id"}:
        fail("E_PLAN_SCHEMA", "layered_run fields differ")
    run_id = layered["run_id"]
    manifest_id = handoff["manifest_id"]
    if not isinstance(run_id, str) or not SAFE_ID.fullmatch(run_id):
        fail("E_PLAN_SCHEMA", "run_id is not a safe identifier")
    if not isinstance(manifest_id, str) or not SAFE_ID.fullmatch(manifest_id):
        fail("E_PLAN_SCHEMA", "manifest_id is not a safe identifier")
    plan = LaunchPlan(
        path=path.resolve(),
        byte_sha256=byte_sha,
        plan_id=plan_id,
        execute_authorized=auth["execute"],
        reviewed_by=auth["reviewed_by"],
        reviewed_at_utc=auth["reviewed_at_utc"],
        production_root=production_root,
        handoff_lock=exact_path(handoff["lock"], "handoff.lock"),
        handoff_workspace=exact_path(handoff["workspace"], "handoff.workspace"),
        source_manifest=exact_path(handoff["source_manifest"], "handoff.source_manifest"),
        source_manifest_failure=exact_path(
            handoff["source_manifest_failure"], "handoff.source_manifest_failure"
        ),
        manifest_id=manifest_id,
        run_parent=exact_path(layered["run_parent"], "layered_run.run_parent"),
        run_id=run_id,
        inputs=_pin_map(raw["inputs"], REQUIRED_INPUT_ROLES, "inputs"),
        executables=_pin_map(raw["executables"], REQUIRED_EXECUTABLE_ROLES, "executables"),
    )
    validate_plan_relationships(plan)
    return plan


def validate_plan_relationships(plan: LaunchPlan) -> None:
    if plan.production_root != EXPECTED_PRODUCTION_ROOT:
        fail("E_PLAN_PATH", f"production_root must be {EXPECTED_PRODUCTION_ROOT}")
    if plan.run_parent != EXPECTED_RUN_PARENT:
        fail("E_PLAN_PATH", f"run_parent must be {EXPECTED_RUN_PARENT}")
    expected_snapshot = plan.production_root / "input_manifest.json"
    fixed = {
        "producer_manifest_snapshot": expected_snapshot,
        "producer_manifest_reviewed": EXPECTED_REVIEWED_PRODUCER_MANIFEST,
        "base_manifest": EXPECTED_BASE_MANIFEST,
        "real_method_receipt": EXPECTED_REAL_RECEIPT,
        "synthetic_method_receipt": EXPECTED_SYNTHETIC_RECEIPT,
        "longphase_binary": EXPECTED_LONGPHASE_BINARY,
        "handoff_supervisor": Path(__file__).resolve(),
        "producer_finalizer": METHOD_SCRIPTS / "finalize_longphase_production_sidecars.py",
        "receipt_normalizer": METHOD_SCRIPTS / "build_longphase_production_capture_receipt_v3.py",
        "v3_preparer": METHOD_SCRIPTS / "prepare_clean_layered_manifest_v3.py",
        "v3_runner": REPO_ROOT / "scripts/run_layered_v3.py",
        "v3_lifecycle": REPO_ROOT / "scripts/layered_v3_lifecycle.py",
        "v3_validator": METHOD_SCRIPTS / "validate_layered_v3_inputs.py",
        "v3_verifier": REPO_ROOT / "scripts/verify_layered_v3.py",
        "schema_manifest": METHOD_SCHEMAS / "layered_input_manifest_v3.schema.json",
        "schema_lock": METHOD_SCHEMAS / "layered_input_lock_v1.schema.json",
        "schema_receipt": METHOD_SCHEMAS / "longphase_production_capture_receipt_v1.schema.json",
        "science_sm_linkage_genomewide": METHOD_SCRIPTS / "sm_linkage_genomewide.py",
        "science_sm_multilocus_combinations": METHOD_SCRIPTS / "sm_multilocus_combinations.py",
        "science_tree_enumeration_solver": METHOD_SCRIPTS / "tree_enumeration_solver.py",
        "science_layered_tree_reconstruction": METHOD_SCRIPTS / "layered_tree_reconstruction.py",
        "science_build_region_view": METHOD_SCRIPTS / "build_region_view.py",
        "science_build_ssnv_site_ledger": METHOD_SCRIPTS / "build_ssnv_site_ledger.py",
        "tool_ps": EXPECTED_PS_BINARY,
    }
    for role, expected in fixed.items():
        collection = plan.inputs if role in plan.inputs else plan.executables
        if collection[role].path != expected:
            fail("E_PLAN_PATH", f"{role} must be {expected}")
    longphase = plan.executables["longphase_binary"]
    if longphase.sha256 != EXPECTED_LONGPHASE_SHA256:
        fail("E_PLAN_HASH", "LongPhase binary digest is not the reviewed production digest")
    expected_lock = plan.production_root.parent / f".{plan.production_root.name}.handoff.lock"
    if plan.handoff_lock != expected_lock:
        fail("E_PLAN_PATH", f"handoff lock must be {expected_lock}")
    expected_handoff_parent = plan.run_parent / ".layered_v3_handoff"
    if plan.handoff_workspace.parent != expected_handoff_parent:
        fail("E_PLAN_PATH", f"handoff workspace must be directly under {expected_handoff_parent}")
    if plan.handoff_workspace.name != plan.plan_id:
        fail("E_PLAN_PATH", "handoff workspace basename must equal plan_id")
    if plan.source_manifest != plan.handoff_workspace / "layered_input_manifest_v3.json":
        fail("E_PLAN_PATH", "source manifest must use the fixed workspace filename")
    if plan.source_manifest_failure != plan.handoff_workspace / "layered_input_manifest_v3.failure.json":
        fail("E_PLAN_PATH", "source-manifest failure path must use the fixed workspace filename")
    if plan.run_root == plan.production_root or plan.handoff_workspace in plan.run_root.parents:
        fail("E_PLAN_PATH", "producer, handoff, and layered run roots must be distinct")


def verify_pins(plan: LaunchPlan) -> None:
    for pin in plan.all_pins:
        pin.verify()
    reviewed = plan.inputs["producer_manifest_reviewed"]
    snapshot = plan.inputs["producer_manifest_snapshot"]
    if reviewed.sha256 != snapshot.sha256:
        fail("E_PRODUCER_MANIFEST_DRIFT", "reviewed producer manifest and frozen snapshot hashes differ")
    if reviewed.path.read_bytes() != snapshot.path.read_bytes():
        fail("E_PRODUCER_MANIFEST_DRIFT", "reviewed producer manifest and frozen snapshot bytes differ")


def read_sha256_manifest(path: Path) -> dict[Path, str]:
    if not path.is_file() or path.is_symlink():
        fail("E_HASH_MANIFEST", f"hash manifest missing/non-regular: {path}")
    result: dict[Path, str] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})\s+\*?(.+)", line)
        if not match:
            fail("E_HASH_MANIFEST", f"malformed line {number}: {path}")
        artifact = Path(match.group(2))
        if not artifact.is_absolute():
            fail("E_HASH_MANIFEST", f"relative artifact on line {number}: {path}")
        resolved = artifact.resolve(strict=False)
        if resolved in result:
            fail("E_HASH_MANIFEST", f"duplicate artifact: {artifact}")
        if not artifact.is_file():
            fail("E_HASH_MISSING", f"hashed artifact missing: {artifact}")
        observed = sha256_file(artifact)
        if observed != match.group(1):
            fail("E_HASH_DRIFT", f"hash mismatch: {artifact}")
        result[resolved] = observed
    if not result:
        fail("E_HASH_MANIFEST", f"empty hash manifest: {path}")
    return result


def load_expected_samples(plan: LaunchPlan) -> tuple[str, ...]:
    manifest = strict_json_load(plan.inputs["producer_manifest_reviewed"].path)
    samples = tuple(item.get("sample") for item in manifest.get("samples", []))
    if samples != EXPECTED_SAMPLES:
        fail("E_DATASET_SET", f"reviewed producer manifest order/set differs: {samples}")
    if manifest.get("dataset_count") != 7 or manifest.get("biological_sample_count") != 6:
        fail("E_DATASET_SET", "reviewed producer manifest is not 7 datasets / 6 biological samples")
    return samples


def inspect_status(path: Path, samples: Sequence[str]) -> tuple[bool, bytes]:
    if not path.is_file():
        raise ProducerPending(f"run_status missing: {path}")
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        fail("E_STATUS_SCHEMA", f"run_status is not UTF-8: {exc}")
    lines = text.splitlines()
    exact_header = "timestamp\tsample\tstage\tstatus\tdetail"
    if not lines or lines[0] != exact_header:
        fail("E_STATUS_SCHEMA", "run_status header differs")
    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames != ["timestamp", "sample", "stage", "status", "detail"]:
        fail("E_STATUS_SCHEMA", "run_status contains duplicate/reordered header columns")
    rows = list(reader)
    states = {sample: "NOT_STARTED" for sample in samples}
    all_seen = False
    for index, row in enumerate(rows):
        sample, stage, status = row["sample"], row["stage"], row["status"]
        if status == "FAIL":
            fail("E_PRODUCER_FAILED", f"run_status contains FAIL at row {index + 2}: {row}")
        if sample == "ALL":
            if (stage, status) != ("verify", "PASS"):
                fail("E_STATUS_UNEXPECTED", f"unexpected ALL status row: {row}")
            if all_seen or index != len(rows) - 1 or any(value != "PASS" for value in states.values()):
                fail("E_STATUS_ORDER", "ALL verify PASS must be unique, last, and follow all sample PASS rows")
            all_seen = True
            continue
        if sample not in states or stage != "production_tagging" or status not in {"START", "PASS"}:
            fail("E_STATUS_UNEXPECTED", f"unexpected status row: {row}")
        if status == "START":
            if states[sample] != "NOT_STARTED":
                fail("E_STATUS_ORDER", f"duplicate/out-of-order START for {sample}")
            states[sample] = "START"
        else:
            if states[sample] != "START":
                fail("E_STATUS_ORDER", f"PASS does not immediately follow a prior START for {sample}")
            states[sample] = "PASS"
    complete = all_seen and all(value == "PASS" for value in states.values())
    return complete, raw


def validate_status(path: Path, samples: Sequence[str]) -> None:
    complete, _ = inspect_status(path, samples)
    if not complete:
        fail("E_STATUS_INCOMPLETE", "run_status has legal progress but is not exact terminal 7/7 + ALL PASS")


def validate_verification(path: Path, production_root: Path, samples: Sequence[str]) -> None:
    if not path.is_file():
        fail("E_PRODUCER_NOT_READY", f"verification missing: {path}")
    value = strict_json_load(path)
    if (
        value.get("dataset_count") != 7
        or value.get("n_pass") != 7
        or value.get("all_pass") is not True
        or len(value.get("samples", [])) != 7
    ):
        fail("E_VERIFICATION_INCOMPLETE", "producer verification is not exact 7/7 all-pass")
    observed: list[str] = []
    for item in value["samples"]:
        if item.get("pass") is not True:
            fail("E_VERIFICATION_INCOMPLETE", "producer verification contains a non-pass sample")
        sidecar = exact_path(item.get("sidecar"), "verification.samples.sidecar")
        try:
            relative = sidecar.relative_to(production_root / "samples")
        except ValueError:
            fail("E_VERIFICATION_SUBJECT", f"verification sidecar escapes production root: {sidecar}")
        if len(relative.parts) < 2:
            fail("E_VERIFICATION_SUBJECT", f"verification sidecar has no sample subject: {sidecar}")
        observed.append(relative.parts[0])
    if Counter(observed) != Counter(samples):
        fail("E_VERIFICATION_SUBJECT", f"verification sample set differs: {observed}")


def system_process_snapshot(ps_binary: Path) -> list[dict[str, Any]]:
    proc = subprocess.run(
        [str(ps_binary), "-eo", "pid=,ppid=,args="], text=True, capture_output=True, check=False
    )
    if proc.returncode != 0:
        fail("E_PROCESS_OBSERVATION", f"ps failed with exit {proc.returncode}")
    result = []
    for line in proc.stdout.splitlines():
        fields = line.strip().split(None, 2)
        if len(fields) == 3:
            result.append({"pid": int(fields[0]), "ppid": int(fields[1]), "argv": fields[2]})
    return result


def process_conflicts(rows: Sequence[Mapping[str, Any]], production_root: Path) -> list[dict[str, Any]]:
    conflicts: list[dict[str, Any]] = []
    root = str(production_root)
    producer_names = {
        "run_longphase_production_sidecars.sh",
        "capture_longphase_tagged_bam_sidecar.py",
        "longphase-s",
    }
    layered_names = {"run_layered_7samples_newbb.sh", "run_layered_v3.py"}
    for row in rows:
        argv = str(row.get("argv", ""))
        tokens = argv.split()
        basenames = {Path(token).name for token in tokens if token and not token.startswith("-")}
        reason = None
        if basenames & producer_names and (
            root in argv or "run_longphase_production_sidecars.sh" in basenames
        ):
            reason = "producer_active"
        elif basenames & layered_names:
            reason = "overlapping_layered_launcher"
        if reason:
            conflicts.append({"pid": row.get("pid"), "ppid": row.get("ppid"), "reason": reason, "argv": argv})
    return conflicts


def validate_producer_common(
    plan: LaunchPlan,
    process_provider: Callable[[], Sequence[Mapping[str, Any]]],
) -> tuple[str, ...]:
    verify_pins(plan)
    samples = load_expected_samples(plan)
    root = plan.production_root
    if not root.is_dir() or root.is_symlink():
        fail("E_PRODUCER_ROOT", f"production root missing/non-directory: {root}")
    params = strict_json_load(root / "params.json")
    if params.get("truth_flags") is not False or params.get("output_mode") != "read_tag_sidecar":
        fail("E_PRODUCTION_POLICY", "producer params do not prove truth-free read-tag sidecar mode")
    read_sha256_manifest(root / "code.sha256")
    status_complete, _ = inspect_status(root / "run_status.tsv", samples)
    conflicts = process_conflicts(process_provider(), root)
    layered = [item for item in conflicts if item["reason"] == "overlapping_layered_launcher"]
    producer = [item for item in conflicts if item["reason"] == "producer_active"]
    if layered:
        fail("E_PROCESS_ACTIVE", f"overlapping layered launcher: {layered}")
    if not status_complete:
        if producer:
            raise ProducerPending("producer has legal non-terminal status and matching active process")
        fail("E_PRODUCER_ABANDONED", "producer status is incomplete but no wrapper/LongPhase/capture process remains")
    validate_verification(root / "verification_summary.json", root, samples)
    if producer:
        raise ProducerPending("producer is 7/7 verified but wrapper/LongPhase/capture has not exited")
    return samples


def validate_closeout(plan: LaunchPlan, samples: Sequence[str]) -> bool:
    root = plan.production_root
    closeout_dir = root / "closeout"
    success_path = root / "_SUCCESS"
    if not closeout_dir.exists() and not success_path.exists():
        return False
    if not closeout_dir.is_dir() or not success_path.is_file():
        fail("E_CLOSEOUT_PARTIAL", "closeout and _SUCCESS must appear atomically as a complete pair")
    receipt_path = closeout_dir / "production_closeout.json"
    hashes_path = closeout_dir / "artifacts.final.sha256"
    if not receipt_path.is_file() or not hashes_path.is_file():
        fail("E_CLOSEOUT_PARTIAL", "closeout receipt/hash manifest missing")
    receipt = strict_json_load(receipt_path)
    receipt_samples = tuple(item.get("sample") for item in receipt.get("samples", []))
    if (
        receipt.get("status") != "PASS"
        or receipt.get("dataset_count") != 7
        or receipt.get("n_pass") != 7
        or receipt.get("truth_flags") is not False
        or receipt.get("binary_sha256") != EXPECTED_LONGPHASE_SHA256
        or Counter(receipt_samples) != Counter(samples)
    ):
        fail("E_CLOSEOUT_INVALID", "production closeout does not prove the exact reviewed contract")
    success = strict_json_load(success_path)
    if (
        success.get("status") != "SUCCESS"
        or Path(success.get("closeout_receipt", "")).resolve() != receipt_path.resolve()
        or success.get("closeout_receipt_sha256") != sha256_file(receipt_path)
        or Path(success.get("artifacts_manifest", "")).resolve() != hashes_path.resolve()
        or success.get("artifacts_manifest_sha256") != sha256_file(hashes_path)
    ):
        fail("E_CLOSEOUT_INVALID", "producer _SUCCESS does not bind closeout bytes")
    read_sha256_manifest(root / "code.sha256")
    read_sha256_manifest(hashes_path)
    return True


def receipt_paths(plan: LaunchPlan, sample: str) -> tuple[Path, Path]:
    return (
        plan.production_root / "samples" / sample / "producer_capture_receipt_v3.json",
        plan.handoff_workspace / "receipt_failures" / f"{sample}.failure.json",
    )


def assert_receipt_state_new(plan: LaunchPlan, samples: Sequence[str]) -> None:
    states = []
    for sample in samples:
        output, failure_path = receipt_paths(plan, sample)
        states.append((sample, output.exists(), failure_path.exists()))
    if any(output or failure for _, output, failure in states):
        fail("E_PARTIAL_RECEIPTS", f"one-shot supervisor requires zero prior receipt/failure outputs: {states}")


def validate_exact_receipts(plan: LaunchPlan, samples: Sequence[str]) -> dict[str, str]:
    receipts: dict[str, str] = {}
    sample_root = plan.production_root / "samples"
    observed = []
    for child in sample_root.iterdir():
        if child.is_dir() and (child / "producer_capture_receipt_v3.json").exists():
            observed.append(child.name)
    if Counter(observed) != Counter(samples):
        fail("E_RECEIPT_COUNT", f"receipt subjects are not exact seven: {observed}")
    for sample in samples:
        output, failure_path = receipt_paths(plan, sample)
        if failure_path.exists() or not output.is_file() or output.is_symlink():
            fail("E_RECEIPT_INVALID", f"receipt/failure state invalid for {sample}")
        value = strict_json_load(output)
        if (
            value.get("schema_name") != "intersubmod.longphase_production_capture_receipt"
            or value.get("schema_version") != "1.0.0"
            or value.get("sample") != sample
        ):
            fail("E_RECEIPT_INVALID", f"receipt identity invalid for {sample}")
        receipts[sample] = sha256_file(output)
    return receipts


class RealExecutor:
    def run(self, argv: Sequence[str], log_path: Path, env: Mapping[str, str]) -> int:
        if log_path.exists():
            fail("E_OUTPUT_EXISTS", f"command log exists: {log_path}")
        with log_path.open("xb") as log:
            completed = subprocess.run(list(argv), stdout=log, stderr=subprocess.STDOUT, env=dict(env), check=False)
        return int(completed.returncode)

    def popen(self, argv: Sequence[str], log_path: Path, env: Mapping[str, str]) -> tuple[Any, Any]:
        if log_path.exists():
            fail("E_OUTPUT_EXISTS", f"runner log exists: {log_path}")
        log = log_path.open("xb")
        try:
            process = subprocess.Popen(
                list(argv), stdout=log, stderr=subprocess.STDOUT, env=dict(env), start_new_session=False
            )
        except BaseException:
            log.close()
            raise
        return process, log


def controlled_environment(plan: LaunchPlan) -> tuple[dict[str, str], dict[str, str]]:
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("SM_"):
            del env[key]
    tool_roles = {
        "samtools": "tool_samtools",
        "bgzip": "tool_bgzip",
        "tabix": "tool_tabix",
        "stat": "tool_stat",
        "bcftools": "tool_bcftools",
    }
    directories: list[str] = []
    for role in tool_roles.values():
        parent = str(plan.executables[role].resolved_path.parent)
        if parent not in directories:
            directories.append(parent)
    env.update({
        "LC_ALL": "C",
        "TZ": "UTC",
        "PYTHONHASHSEED": "0",
        "PYTHONDONTWRITEBYTECODE": "1",
        "PATH": os.pathsep.join(directories),
    })
    resolutions: dict[str, str] = {}
    for command, role in tool_roles.items():
        observed = shutil.which(command, path=env["PATH"])
        if observed is None or Path(observed).resolve() != plan.executables[role].resolved_path:
            fail(
                "E_TOOL_RESOLUTION",
                f"controlled PATH does not resolve {command} to pin {plan.executables[role].resolved_path}; observed={observed}",
            )
        resolutions[command] = str(Path(observed).resolve())
    return env, resolutions


def command_paths(plan: LaunchPlan) -> dict[str, str]:
    return {role: str(pin.resolved_path) for role, pin in plan.executables.items()}


def validate_command_policy(commands: Sequence[tuple[str, Sequence[str]]]) -> None:
    forbidden = ("--truth-vcf", "--truth-bed")
    for stage, argv in commands:
        for token in argv:
            if any(token == option or token.startswith(f"{option}=") for option in forbidden):
                fail("E_COMMAND_POLICY", f"{stage} contains forbidden production option {token}")


def build_commands(plan: LaunchPlan, samples: Sequence[str]) -> list[tuple[str, list[str]]]:
    exe = command_paths(plan)
    commands: list[tuple[str, list[str]]] = []
    if not (plan.production_root / "_SUCCESS").exists():
        commands.append(("producer_closeout", [
            exe["python"], exe["producer_finalizer"],
            "--run-root", str(plan.production_root),
            "--expected-manifest", str(plan.inputs["producer_manifest_reviewed"].path),
            "--longphase-binary", exe["longphase_binary"],
            "--expected-binary-sha256", EXPECTED_LONGPHASE_SHA256,
        ]))
    for sample in samples:
        output, failure_path = receipt_paths(plan, sample)
        commands.append((f"receipt:{sample}", [
            exe["python"], exe["receipt_normalizer"],
            "--sample-dir", str(plan.production_root / "samples" / sample),
            "--production-manifest", str(plan.inputs["producer_manifest_snapshot"].path),
            "--sample", sample,
            "--run-root", str(plan.production_root),
            "--longphase-binary", exe["longphase_binary"],
            "--output", str(output),
            "--failure-report", str(failure_path),
        ]))
    commands.append(("prepare_v3_manifest", [
        exe["python"], exe["v3_preparer"],
        "--base-manifest", str(plan.inputs["base_manifest"].path),
        "--longphase-manifest", str(plan.inputs["producer_manifest_reviewed"].path),
        "--production-root", str(plan.production_root),
        "--real-data-receipt", str(plan.inputs["real_method_receipt"].path),
        "--synthetic-receipt", str(plan.inputs["synthetic_method_receipt"].path),
        "--manifest-id", plan.manifest_id,
        "--output", str(plan.source_manifest),
        "--failure-report", str(plan.source_manifest_failure),
    ]))
    commands.append(("run_layered_v3", [
        exe["python"], exe["v3_runner"],
        "--manifest", str(plan.source_manifest),
        "--run-parent", str(plan.run_parent),
        "--run-id", plan.run_id,
        *RUNNER_OPTIONS,
    ]))
    validate_command_policy(commands)
    return commands


def validate_prepared_manifest(plan: LaunchPlan, samples: Sequence[str]) -> str:
    if plan.source_manifest_failure.exists() or not plan.source_manifest.is_file():
        fail("E_MANIFEST_BUILD", "v3 source manifest absent or failure report present")
    value = strict_json_load(plan.source_manifest)
    observed = tuple(item.get("sample") for item in value.get("samples", []))
    if (
        value.get("schema_name") != "intersubmod.layered_input_manifest"
        or value.get("schema_version") != "3.0.0"
        or value.get("manifest_id") != plan.manifest_id
        or value.get("dataset_count") != 7
        or value.get("biological_sample_count") != 6
        or observed != tuple(samples)
    ):
        fail("E_MANIFEST_BUILD", "prepared v3 manifest identity/scope readback failed")
    return sha256_file(plan.source_manifest)


def validate_layered_success(plan: LaunchPlan) -> dict[str, Any]:
    marker_path = plan.run_root / "_SUCCESS"
    if not marker_path.is_file() or marker_path.is_symlink():
        fail("E_LAYERED_NO_SUCCESS", f"layered runner exited without _SUCCESS: {marker_path}")
    marker = strict_json_load(marker_path)
    verification_path = Path(marker.get("verification_path", ""))
    if (
        marker.get("schema_name") != "intersubmod.layered_success_marker"
        or marker.get("schema_version") != "1.0.0"
        or marker.get("run_id") != plan.run_id
        or marker.get("extra", {}).get("mode") != "full"
        or marker.get("extra", {}).get("dataset_count") != 7
        or marker.get("extra", {}).get("scope") != "chr1-22"
    ):
        fail("E_LAYERED_NO_SUCCESS", "layered _SUCCESS identity/scope is invalid")
    try:
        verification_path.resolve().relative_to(plan.run_root.resolve())
    except ValueError:
        fail("E_LAYERED_NO_SUCCESS", "layered verification path escapes run root")
    if not verification_path.is_file() or sha256_file(verification_path) != marker.get("verification_sha256"):
        fail("E_LAYERED_NO_SUCCESS", "layered verification hash binding failed")
    verification = strict_json_load(verification_path)
    if verification.get("all_pass") is not True or verification.get("dataset_count") != 7:
        fail("E_LAYERED_NO_SUCCESS", "layered verification is not exact 7/7 pass")
    return {"success_sha256": sha256_file(marker_path), "verification_sha256": sha256_file(verification_path)}


class HandoffSupervisor:
    def __init__(
        self,
        plan: LaunchPlan,
        *,
        process_provider: Callable[[], Sequence[Mapping[str, Any]]] | None = None,
        executor: Any | None = None,
        sleeper: Callable[[float], None] = time.sleep,
        poll_seconds: float = 30.0,
    ):
        if poll_seconds < 0 or poll_seconds > 60:
            fail("E_POLL_INTERVAL", "poll interval must be between 0 and 60 seconds")
        self.plan = plan
        self.process_provider = process_provider or (
            lambda: system_process_snapshot(plan.executables["tool_ps"].resolved_path)
        )
        self.executor = executor or RealExecutor()
        self.sleeper = sleeper
        self.poll_seconds = poll_seconds

    def _lock(self) -> Any:
        self.plan.handoff_lock.parent.mkdir(parents=True, exist_ok=True)
        handle = self.plan.handoff_lock.open("a+b")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            handle.close()
            fail("E_HANDOFF_LOCKED", f"another handoff supervisor holds {self.plan.handoff_lock}")
        return handle

    def _static_preflight(self) -> tuple[str, ...]:
        samples = validate_producer_common(self.plan, self.process_provider)
        if self.plan.run_root.exists() or self.plan.run_root.is_symlink():
            fail("E_RUN_ROOT_EXISTS", f"layered run-id already exists: {self.plan.run_root}")
        if self.plan.handoff_workspace.exists() or self.plan.handoff_workspace.is_symlink():
            fail("E_HANDOFF_EXISTS", f"one-shot handoff workspace already exists: {self.plan.handoff_workspace}")
        assert_receipt_state_new(self.plan, samples)
        validate_closeout(self.plan, samples)
        return samples

    def _preflight_or_wait(self, wait_for_producer: bool) -> tuple[str, ...]:
        previous_status: bytes | None = None
        while True:
            status_path = self.plan.production_root / "run_status.tsv"
            current_status = status_path.read_bytes() if status_path.is_file() else b""
            if previous_status is not None and not current_status.startswith(previous_status):
                fail("E_STATUS_REGRESSION", "run_status bytes were truncated or rewritten while waiting")
            previous_status = current_status
            try:
                return self._static_preflight()
            except ProducerPending:
                if not wait_for_producer:
                    raise
                self.sleeper(self.poll_seconds)

    def dry_run(self, *, wait_for_producer: bool = False) -> dict[str, Any]:
        with self._lock():
            samples = self._preflight_or_wait(wait_for_producer)
            commands = build_commands(self.plan, samples)
            return {
                "schema_name": "intersubmod.layered_v3_handoff_dry_run",
                "schema_version": "1.0.0",
                "pass": True,
                "mode": "dry_run",
                "plan_sha256": self.plan.byte_sha256,
                "dataset_count": len(samples),
                "subprocesses_started": 0,
                "commands": [{"stage": stage, "argv": argv} for stage, argv in commands],
            }

    def execute(self, *, wait_for_producer: bool = False) -> dict[str, Any]:
        with self._lock():
            samples = self._preflight_or_wait(wait_for_producer)
            self.plan.handoff_workspace.parent.mkdir(parents=True, exist_ok=True)
            self.plan.handoff_workspace.mkdir(mode=0o700)
            atomic_write_json(
                self.plan.handoff_workspace / "reviewed_plan_binding.json",
                {
                    "schema_name": "intersubmod.reviewed_handoff_plan_binding",
                    "schema_version": "1.0.0",
                    "plan_path": str(self.plan.path),
                    "plan_sha256": self.plan.byte_sha256,
                    "reviewed_by": self.plan.reviewed_by,
                    "reviewed_at_utc": self.plan.reviewed_at_utc,
                },
            )
            env, tool_resolution = controlled_environment(self.plan)
            commands = dict(build_commands(self.plan, samples))
            closeout_exists = validate_closeout(self.plan, samples)
            if not closeout_exists:
                self._run_stage("producer_closeout", commands["producer_closeout"], env)
                validate_producer_common(self.plan, self.process_provider)
                if not validate_closeout(self.plan, samples):
                    fail("E_CLOSEOUT_MISSING", "finalizer exited zero without immutable closeout")
            for sample in samples:
                validate_status(self.plan.production_root / "run_status.tsv", samples)
                validate_closeout(self.plan, samples)
                verify_pins(self.plan)
                self._run_stage(f"receipt_{sample}", commands[f"receipt:{sample}"], env)
                output, failure_path = receipt_paths(self.plan, sample)
                if failure_path.exists() or not output.is_file():
                    fail("E_RECEIPT_INVALID", f"normalizer did not publish a clean receipt for {sample}")
            receipt_hashes = validate_exact_receipts(self.plan, samples)
            validate_closeout(self.plan, samples)
            verify_pins(self.plan)
            self._run_stage("prepare_v3_manifest", commands["prepare_v3_manifest"], env)
            manifest_sha = validate_prepared_manifest(self.plan, samples)

            # Final launch gate: repeat every mutable/operational observation.
            validate_producer_common(self.plan, self.process_provider)
            validate_closeout(self.plan, samples)
            validate_exact_receipts(self.plan, samples)
            verify_pins(self.plan)
            if sha256_file(self.plan.source_manifest) != manifest_sha:
                fail("E_MANIFEST_DRIFT", "v3 source manifest changed after preparation")
            if self.plan.run_root.exists() or self.plan.run_root.is_symlink():
                fail("E_RUN_ROOT_EXISTS", f"layered run-id appeared before launch: {self.plan.run_root}")

            runner_argv = commands["run_layered_v3"]
            runner_log = self.plan.handoff_workspace / "run_layered_v3.foreground.log"
            atomic_write_json(
                self.plan.handoff_workspace / "runner.planned.json",
                {
                    "schema_name": "intersubmod.layered_v3_runner_plan",
                    "schema_version": "1.0.0",
                    "argv": runner_argv,
                    "environment": {
                        key: env[key]
                        for key in ("LC_ALL", "TZ", "PYTHONHASHSEED", "PYTHONDONTWRITEBYTECODE", "PATH")
                    },
                    "tool_resolution": tool_resolution,
                    "log": str(runner_log),
                    "source_manifest_sha256": manifest_sha,
                    "receipt_sha256": receipt_hashes,
                },
            )
            process, log_handle = self.executor.popen(runner_argv, runner_log, env)
            atomic_write_json(
                self.plan.handoff_workspace / "runner.started.json",
                {
                    "schema_name": "intersubmod.layered_v3_runner_started",
                    "schema_version": "1.0.0",
                    "pid": int(process.pid),
                    "argv": runner_argv,
                    "log": str(runner_log),
                },
            )
            try:
                runner_exit = int(process.wait())
            finally:
                log_handle.close()
            atomic_write_json(
                self.plan.handoff_workspace / "runner.exit.json",
                {
                    "schema_name": "intersubmod.layered_v3_runner_exit",
                    "schema_version": "1.0.0",
                    "pid": int(process.pid),
                    "exit_code": runner_exit,
                    "argv": runner_argv,
                    "log": str(runner_log),
                },
            )
            if runner_exit != 0:
                fail("E_LAYERED_RUNNER_EXIT", f"layered runner exited {runner_exit}", exit_code=runner_exit)
            success = validate_layered_success(self.plan)
            result = {
                "schema_name": "intersubmod.layered_v3_handoff_success",
                "schema_version": "1.0.0",
                "pass": True,
                "plan_sha256": self.plan.byte_sha256,
                "production_root": str(self.plan.production_root),
                "receipt_sha256": receipt_hashes,
                "source_manifest": str(self.plan.source_manifest),
                "source_manifest_sha256": manifest_sha,
                "layered_run_root": str(self.plan.run_root),
                **success,
            }
            atomic_write_json(self.plan.handoff_workspace / "handoff.success.json", result)
            return result

    def _run_stage(self, label: str, argv: Sequence[str], env: Mapping[str, str]) -> None:
        verify_pins(self.plan)
        log_path = self.plan.handoff_workspace / "logs" / f"{label}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        rc = int(self.executor.run(argv, log_path, env))
        atomic_write_json(
            self.plan.handoff_workspace / "stage_receipts" / f"{label}.json",
            {
                "schema_name": "intersubmod.layered_v3_handoff_stage",
                "schema_version": "1.0.0",
                "stage": label,
                "argv": list(argv),
                "exit_code": rc,
                "log": str(log_path),
                "log_sha256": sha256_file(log_path) if log_path.is_file() else None,
            },
        )
        if rc != 0:
            fail("E_STAGE_EXIT", f"{label} exited {rc}", exit_code=rc)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--launch-plan", required=True, type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="read-only default; start no subprocess")
    mode.add_argument(
        "--execute-reviewed-plan",
        action="store_true",
        help="execute only when JSON authorization and out-of-band plan hash also match",
    )
    parser.add_argument(
        "--wait-for-producer",
        action="store_true",
        help="hold the handoff flock and poll only exact producer authority files/processes every 30 seconds",
    )
    parser.add_argument("--expected-launch-plan-sha256")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        plan = load_launch_plan(args.launch_plan)
        execute = bool(args.execute_reviewed_plan)
        if args.expected_launch_plan_sha256 is not None:
            if not HEX64.fullmatch(args.expected_launch_plan_sha256):
                fail("E_PLAN_HASH", "--expected-launch-plan-sha256 must be 64 lowercase hex")
            if args.expected_launch_plan_sha256 != plan.byte_sha256:
                fail("E_PLAN_HASH", "out-of-band launch-plan hash does not match plan bytes")
        if execute:
            if args.expected_launch_plan_sha256 is None:
                fail("E_EXECUTION_GATE", "execution requires --expected-launch-plan-sha256")
            if not plan.execute_authorized or not plan.reviewed_by or not plan.reviewed_at_utc:
                fail("E_EXECUTION_GATE", "reviewed JSON authorization is absent/incomplete")
        supervisor = HandoffSupervisor(plan)
        result = (
            supervisor.execute(wait_for_producer=args.wait_for_producer)
            if execute
            else supervisor.dry_run(wait_for_producer=args.wait_for_producer)
        )
        print(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2))
        return 0
    except HandoffError as exc:
        print(
            json.dumps({"pass": False, "error_code": exc.code, "message": exc.message}, ensure_ascii=False),
            file=sys.stderr,
        )
        return exc.exit_code if 1 <= exc.exit_code <= 125 else 70
    except BaseException as exc:
        print(
            json.dumps({"pass": False, "error_code": "E_INTERNAL", "message": repr(exc)}, ensure_ascii=False),
            file=sys.stderr,
        )
        return 70


if __name__ == "__main__":
    raise SystemExit(main())
