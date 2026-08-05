#!/usr/bin/env python3
"""Run the HCC1395 chr1-22 exact-PS x HP k<=12 PARTIAL pilot.

This runner is intentionally fail closed.  It verifies the frozen Big7 input
binding, compiles the independent C++ kernel once, executes every requested
chromosome through extraction, Python partitioning, normalized-input
decompression, C++ parity, and comparison, then writes a single auditable
receipt.  It never resumes or overwrites a non-empty output root.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Mapping, Sequence


TOPIC_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TOPIC_ROOT.parents[1]
DATASET = "HCC1395"
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))
SCHEMA_NAME = "intersubmod.hcc1395_exact_ps_k12_partial_pilot"
SCHEMA_VERSION = "1.0.0"
INPUT_BINDING_SCHEMA_NAME = "intersubmod.big7_hcc1395_input_binding_receipt"
INPUT_BINDING_SCHEMA_VERSION = "1.0.0"
EXTRACTION_SCHEMA_NAME = "intersubmod.lossless_read_linkage_chromosome_receipt"
EXTRACTION_SCHEMA_VERSION = "1.3.0"

DEFAULT_VERIFIER = TOPIC_ROOT / "scripts" / "verify_big7_input_binding.py"
DEFAULT_EXTRACTOR = (
    REPO_ROOT
    / "research"
    / "20260718_k_gt8_read_supported_segmentation"
    / "scripts"
    / "extract_lossless_read_linkage_collapsing.py"
)
DEFAULT_PARTITIONER = TOPIC_ROOT / "scripts" / "exact_ps_k12_partition.py"
DEFAULT_CPP_SOURCE = TOPIC_ROOT / "cpp" / "exact_ps_k12_partition.cpp"
DEFAULT_COMPARATOR = TOPIC_ROOT / "scripts" / "compare_python_cpp_partitions.py"

EXTRACTION_SUFFIXES = (
    "targets.bed",
    "molecule_sparse_calls.tsv.gz",
    "site_catalog.tsv.gz",
    "cut_support.tsv.gz",
    "components.tsv.gz",
    "site_component_membership.tsv.gz",
)
GZIP_EXTRACTION_SUFFIXES = frozenset(EXTRACTION_SUFFIXES[1:])

SUMMARY_FIELDS = (
    "chrom",
    "status",
    "S",
    "units",
    "unit_memberships",
    "unique_sites",
    "k_eq_1_units",
    "k_2_to_8_units",
    "k_9_to_12_units",
    "k_gt_12_units",
    "blocks",
    "constraints",
    "pattern_count_total",
    "retained_constraints",
    "cut_constraints",
    "unavoidable_constraints",
    "molecule_weight_total",
    "molecule_weight_retained",
    "molecule_weight_cut",
    "molecule_weight_unavoidable",
    "cross_ps_violations",
    "cross_hp_violations",
    "python_cpp_mismatch_count",
    "historical_artifacts_compared",
    "historical_semantic_mismatches",
    "historical_physical_only_differences",
    "stage_receipt",
    "stage_receipt_sha256",
)


class RunnerError(RuntimeError):
    """A fail-closed orchestration or evidence-contract failure."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise RunnerError(f"not a regular file: {path}")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def strict_json_load(path: Path) -> Any:
    def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise RunnerError(f"duplicate JSON key in {path}: {key}")
            result[key] = value
        return result

    try:
        return json.loads(
            path.read_text(encoding="utf-8", errors="strict"),
            object_pairs_hook=pairs_hook,
            parse_constant=lambda value: (_ for _ in ()).throw(
                RunnerError(f"non-finite JSON constant in {path}: {value}")
            ),
        )
    except RunnerError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RunnerError(f"cannot read strict JSON {path}: {exc}") from exc


def write_json_exclusive(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            indent=2,
        )
        + "\n"
    )
    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise RunnerError(f"refusing to overwrite output: {path}") from exc


def write_sha256_sidecar(path: Path) -> Path:
    sidecar = path.with_name(f"{path.name}.sha256")
    try:
        with sidecar.open("x", encoding="ascii") as handle:
            handle.write(f"{sha256_path(path)}  {path.name}\n")
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise RunnerError(f"refusing to overwrite checksum: {sidecar}") from exc
    return sidecar


def verify_sha256_sidecar(path: Path) -> dict[str, Any]:
    sidecar = path.with_name(f"{path.name}.sha256")
    if not sidecar.is_file():
        raise RunnerError(f"missing SHA-256 sidecar: {sidecar}")
    fields = sidecar.read_text(encoding="ascii").strip().split()
    observed = sha256_path(path)
    if len(fields) != 2 or fields[1] != path.name or fields[0] != observed:
        raise RunnerError(f"invalid SHA-256 sidecar: {sidecar}")
    return file_identity(sidecar)


def parse_chromosomes(value: str) -> tuple[str, ...]:
    tokens = tuple(token.strip() for token in value.split(",") if token.strip())
    if not tokens:
        raise argparse.ArgumentTypeError("chromosome list must not be empty")
    invalid = [chrom for chrom in tokens if chrom not in AUTOSOMES]
    if invalid:
        raise argparse.ArgumentTypeError(
            f"only chr1-chr22 are allowed; invalid={','.join(invalid)}"
        )
    if len(tokens) != len(set(tokens)):
        raise argparse.ArgumentTypeError("chromosome list contains duplicates")
    return tokens


def resolve_executable(value: str | Path, label: str) -> Path:
    raw = str(value)
    candidate = Path(raw)
    resolved_text = shutil.which(raw) if candidate.name == raw else None
    resolved = Path(resolved_text) if resolved_text else candidate.expanduser().resolve()
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise RunnerError(f"{label} is not an executable file: {resolved}")
    return resolved


def require_source(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise RunnerError(f"{label} is missing: {resolved}")
    return resolved


def prepare_output_root(path: Path) -> Path:
    output_root = path.expanduser().absolute()
    if output_root.is_symlink():
        raise RunnerError(f"output root must not be a symlink: {output_root}")
    if output_root.exists():
        if not output_root.is_dir():
            raise RunnerError(f"output root is not a directory: {output_root}")
        if next(output_root.iterdir(), None) is not None:
            raise RunnerError(f"output root must be new or empty: {output_root}")
    else:
        output_root.mkdir(parents=True, exist_ok=False)
    return output_root.resolve()


def _log_identity(path: Path) -> dict[str, Any]:
    return file_identity(path)


def run_command(
    stage: str,
    command: Sequence[str],
    stdout_log: Path,
    stderr_log: Path,
) -> dict[str, Any]:
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    started_utc = utc_now()
    started = time.perf_counter()
    exit_code = 127
    launch_error: str | None = None
    with stdout_log.open("xb") as stdout, stderr_log.open("xb") as stderr:
        try:
            completed = subprocess.run(
                list(command),
                cwd=str(REPO_ROOT),
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
            exit_code = completed.returncode
        except OSError as exc:
            launch_error = f"{type(exc).__name__}: {exc}"
            stderr.write((launch_error + "\n").encode("utf-8"))
    return {
        "stage": stage,
        "started_at_utc": started_utc,
        "ended_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - started,
        "command": list(command),
        "command_shell": shlex.join(command),
        "exit_code": exit_code,
        "launch_error": launch_error,
        "logs": {
            "stdout": _log_identity(stdout_log),
            "stderr": _log_identity(stderr_log),
        },
    }


def _identity_matches(spec: Mapping[str, Any], path: Path) -> bool:
    observed = file_identity(path)
    size = spec.get("size_bytes", spec.get("bytes"))
    return (
        str(Path(str(spec.get("path", ""))).resolve()) == observed["path"]
        and size == observed["size_bytes"]
        and spec.get("sha256") == observed["sha256"]
    )


def validate_declared_outputs(
    document: Mapping[str, Any], expected_parent: Path
) -> dict[str, Any]:
    declared = document.get("outputs")
    if declared is None:
        return {"declared": False, "verified": 0}
    if not isinstance(declared, Mapping):
        raise RunnerError("child outputs field is not an object")
    verified = 0
    parent = expected_parent.resolve()
    for name, raw in declared.items():
        if not isinstance(raw, Mapping):
            raise RunnerError(f"declared output identity is not an object: {name}")
        path = Path(str(raw.get("path", "")))
        if not path.is_file() or path.resolve().parent != parent:
            raise RunnerError(f"declared output escaped or is missing: {name} -> {path}")
        if not _identity_matches(raw, path):
            raise RunnerError(f"declared output identity mismatch: {path}")
        verified += 1
    return {"declared": True, "verified": verified}


def child_receipt_info(
    path: Path,
    label: str,
    *,
    expected_output_parent: Path | None = None,
    require_sidecar: bool = False,
) -> dict[str, Any]:
    if not path.is_file():
        raise RunnerError(f"{label} receipt is missing: {path}")
    document = strict_json_load(path)
    if not isinstance(document, Mapping):
        raise RunnerError(f"{label} receipt root is not an object")
    sidecar = None
    if require_sidecar:
        sidecar = verify_sha256_sidecar(path)
    declared = None
    if expected_output_parent is not None:
        declared = validate_declared_outputs(document, expected_output_parent)
    return {
        "identity": file_identity(path),
        "sha256_sidecar": sidecar,
        "schema_name": document.get("schema_name"),
        "schema_version": document.get("schema_version"),
        "all_pass": document.get("all_pass") is True,
        "declared_outputs": declared,
        "document": document,
    }


def validate_input_binding_receipt(
    path: Path,
    *,
    expected_manifest: Mapping[str, Any],
    expected_samtools: Path,
) -> dict[str, Any]:
    document = strict_json_load(path)
    if not isinstance(document, Mapping):
        raise RunnerError("input binding receipt root is not an object")
    expected_claim = {
        "schema_name": INPUT_BINDING_SCHEMA_NAME,
        "schema_version": INPUT_BINDING_SCHEMA_VERSION,
        "task_type": "exploratory_pilot",
        "claim_status": "PARTIAL",
        "validation_evidence_eligible": False,
        "sample": DATASET,
        "all_pass": True,
        "failure": None,
    }
    observed_claim = {key: document.get(key) for key in expected_claim}
    if observed_claim != expected_claim:
        raise RunnerError(
            "input binding receipt schema/claim mismatch: "
            f"expected={expected_claim!r} observed={observed_claim!r}"
        )

    integrity = document.get("receipt_integrity")
    if not isinstance(integrity, Mapping) or set(integrity) != {"policy", "sha256"}:
        raise RunnerError("input binding receipt integrity record is malformed")
    if integrity.get("policy") != "semantic_json_sha256_without_receipt_integrity":
        raise RunnerError("input binding receipt integrity policy is unsupported")
    integrity_payload = copy.deepcopy(dict(document))
    del integrity_payload["receipt_integrity"]
    observed_integrity = canonical_json_sha256(integrity_payload)
    if integrity.get("sha256") != observed_integrity:
        raise RunnerError("input binding receipt semantic integrity SHA-256 mismatch")

    declared_manifest = document.get("manifest")
    if not isinstance(declared_manifest, Mapping):
        raise RunnerError("input binding receipt manifest identity is missing")
    for key in ("path", "size_bytes", "sha256"):
        if declared_manifest.get(key) != expected_manifest.get(key):
            raise RunnerError(
                f"input binding receipt manifest {key} mismatch: "
                f"expected={expected_manifest.get(key)!r} "
                f"observed={declared_manifest.get(key)!r}"
            )

    checks = document.get("checks")
    if (
        not isinstance(checks, list)
        or not checks
        or any(not isinstance(row, Mapping) or row.get("pass") is not True for row in checks)
    ):
        raise RunnerError("input binding receipt does not contain a non-empty all-PASS check set")

    observed = document.get("observed")
    quickcheck = observed.get("samtools_quickcheck") if isinstance(observed, Mapping) else None
    if not isinstance(quickcheck, Mapping):
        raise RunnerError("input binding receipt samtools quickcheck record is missing")
    expected_executable = str(expected_samtools.resolve())
    if (
        quickcheck.get("executable") != expected_executable
        or quickcheck.get("returncode") != 0
        or not isinstance(quickcheck.get("command"), list)
        or not quickcheck["command"]
        or quickcheck["command"][0] != expected_executable
    ):
        raise RunnerError("input binding receipt used a different samtools executable")

    snapshot = {
        "manifest": document.get("manifest"),
        "source_authority": document.get("source_authority"),
        "binding_contract": document.get("binding_contract"),
        "observed": document.get("observed"),
        "checks": document.get("checks"),
    }
    return {
        "schema_name": document["schema_name"],
        "schema_version": document["schema_version"],
        "sample": document["sample"],
        "check_count": len(checks),
        "receipt_integrity_sha256": observed_integrity,
        "binding_snapshot_sha256": canonical_json_sha256(snapshot),
        "samtools_executable": expected_executable,
    }


def attach_input_binding_validation(
    stage: dict[str, Any],
    receipt_path: Path,
    *,
    expected_manifest: Mapping[str, Any],
    expected_samtools: Path,
) -> dict[str, Any]:
    try:
        validation = validate_input_binding_receipt(
            receipt_path,
            expected_manifest=expected_manifest,
            expected_samtools=expected_samtools,
        )
    except Exception as exc:
        stage["input_binding_validation"] = None
        stage["all_pass"] = False
        stage["failure"] = f"{type(exc).__name__}: {exc}"
        return stage
    stage["input_binding_validation"] = validation
    return stage


def validate_extraction_receipt(
    path: Path,
    *,
    chrom: str,
    expected_manifest: Mapping[str, Any],
    expected_samtools: Path,
) -> dict[str, Any]:
    document = strict_json_load(path)
    if not isinstance(document, Mapping):
        raise RunnerError("extraction receipt root is not an object")
    if (
        document.get("schema_name") != EXTRACTION_SCHEMA_NAME
        or document.get("schema_version") != EXTRACTION_SCHEMA_VERSION
        or document.get("all_pass") is not True
        or document.get("scope", {}).get("dataset") != DATASET
        or document.get("scope", {}).get("chrom") != chrom
    ):
        raise RunnerError(f"{chrom} extraction receipt schema/scope mismatch")
    provenance_manifest = document.get("provenance", {}).get("manifest")
    if not isinstance(provenance_manifest, Mapping):
        raise RunnerError(f"{chrom} extraction manifest provenance is missing")
    for key in ("path", "sha256"):
        if provenance_manifest.get(key) != expected_manifest.get(key):
            raise RunnerError(f"{chrom} extraction manifest {key} mismatch")
    expected_executable = str(expected_samtools.resolve())
    command = document.get("command")
    parameters = document.get("parameters")
    if (
        not isinstance(command, list)
        or not command
        or command[0] != expected_executable
        or not isinstance(parameters, Mapping)
        or parameters.get("samtools_executable") != expected_executable
    ):
        raise RunnerError(f"{chrom} extraction used a different samtools executable")
    return {
        "schema_name": document["schema_name"],
        "schema_version": document["schema_version"],
        "sample": DATASET,
        "chrom": chrom,
        "manifest_sha256": provenance_manifest["sha256"],
        "samtools_executable": expected_executable,
    }


def external_child_stage(
    *,
    stage: str,
    command: Sequence[str],
    logs_dir: Path,
    receipt_path: Path,
    required_outputs: Sequence[Path] = (),
    expected_output_parent: Path | None = None,
    require_receipt_sidecar: bool = False,
) -> dict[str, Any]:
    record = run_command(
        stage,
        command,
        logs_dir / f"{stage}.stdout.log",
        logs_dir / f"{stage}.stderr.log",
    )
    validation_error = None
    child = None
    required: dict[str, Any] = {}
    try:
        child = child_receipt_info(
            receipt_path,
            stage,
            expected_output_parent=expected_output_parent,
            require_sidecar=require_receipt_sidecar,
        )
        for path in required_outputs:
            required[path.name] = file_identity(path)
    except Exception as exc:
        validation_error = f"{type(exc).__name__}: {exc}"
    record["child_receipt"] = (
        None
        if child is None
        else {key: value for key, value in child.items() if key != "document"}
    )
    record["required_outputs"] = required
    record["validation_error"] = validation_error
    record["all_pass"] = (
        record["exit_code"] == 0
        and validation_error is None
        and child is not None
        and child["all_pass"]
    )
    if not record["all_pass"]:
        if child is not None and not child["all_pass"]:
            record["failure"] = f"{stage} child receipt all_pass=false"
        elif record["exit_code"] != 0:
            record["failure"] = f"{stage} command exited {record['exit_code']}"
        else:
            record["failure"] = validation_error
    return record


def simple_external_stage(
    *,
    stage: str,
    command: Sequence[str],
    logs_dir: Path,
    required_outputs: Sequence[Path] = (),
) -> dict[str, Any]:
    record = run_command(
        stage,
        command,
        logs_dir / f"{stage}.stdout.log",
        logs_dir / f"{stage}.stderr.log",
    )
    validation_error = None
    required: dict[str, Any] = {}
    try:
        for path in required_outputs:
            required[path.name] = file_identity(path)
    except Exception as exc:
        validation_error = f"{type(exc).__name__}: {exc}"
    record["required_outputs"] = required
    record["validation_error"] = validation_error
    record["all_pass"] = record["exit_code"] == 0 and validation_error is None
    if not record["all_pass"]:
        record["failure"] = (
            f"{stage} command exited {record['exit_code']}"
            if record["exit_code"] != 0
            else validation_error
        )
    return record


def semantic_file_identity(path: Path, *, gzip_encoded: bool) -> dict[str, Any]:
    if not path.is_file():
        raise RunnerError(f"comparison artifact is missing: {path}")
    digest = hashlib.sha256()
    semantic_bytes = 0
    opener = gzip.open if gzip_encoded else Path.open
    try:
        with opener(path, "rb") as handle:  # type: ignore[arg-type]
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
                semantic_bytes += len(block)
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise RunnerError(f"cannot read semantic bytes from {path}: {exc}") from exc
    physical = file_identity(path)
    return {
        "path": physical["path"],
        "physical_bytes": physical["size_bytes"],
        "physical_sha256": physical["sha256"],
        "semantic_encoding": "gzip_decompressed_bytes" if gzip_encoded else "direct_bytes",
        "semantic_bytes": semantic_bytes,
        "semantic_sha256": digest.hexdigest(),
    }


def compare_semantic_artifacts(
    current: Path, historical: Path, *, gzip_encoded: bool
) -> dict[str, Any]:
    current_identity = semantic_file_identity(current, gzip_encoded=gzip_encoded)
    historical_identity = semantic_file_identity(
        historical, gzip_encoded=gzip_encoded
    )
    semantic_equal = (
        current_identity["semantic_bytes"] == historical_identity["semantic_bytes"]
        and current_identity["semantic_sha256"]
        == historical_identity["semantic_sha256"]
    )
    physical_equal = (
        current_identity["physical_bytes"] == historical_identity["physical_bytes"]
        and current_identity["physical_sha256"]
        == historical_identity["physical_sha256"]
    )
    return {
        "current": current_identity,
        "historical": historical_identity,
        "semantic_equal": semantic_equal,
        "physical_equal": physical_equal,
        "physical_only_difference": semantic_equal and not physical_equal,
    }


def find_historical_artifact(root: Path, chrom: str, filename: str) -> Path:
    candidates = (
        root / filename,
        root / chrom / filename,
        root / chrom / "extract" / filename,
        root / chrom / "extraction" / filename,
        root / "chromosomes" / chrom / "extraction" / filename,
        root / f"{DATASET}_{chrom}" / filename,
        root / f"{DATASET}_{chrom}" / "extraction_v2" / filename,
    )
    direct = {path.resolve() for path in candidates if path.is_file()}
    if not direct:
        direct = {path.resolve() for path in root.rglob(filename) if path.is_file()}
    if len(direct) != 1:
        raise RunnerError(
            f"historical artifact must resolve uniquely: {filename}; matches={len(direct)}"
        )
    return next(iter(direct))


def historical_comparison_stage(
    *,
    chrom: str,
    extraction_dir: Path,
    historical_root: Path,
    logs_dir: Path,
) -> dict[str, Any]:
    stage = "historical_semantic_compare"
    stdout_path = logs_dir / f"{stage}.stdout.log"
    stderr_path = logs_dir / f"{stage}.stderr.log"
    logs_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    started_utc = utc_now()
    comparisons: dict[str, Any] = {}
    error = None
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        try:
            for suffix in EXTRACTION_SUFFIXES:
                filename = f"{DATASET}.{chrom}.{suffix}"
                historical = find_historical_artifact(
                    historical_root, chrom, filename
                )
                comparison = compare_semantic_artifacts(
                    extraction_dir / filename,
                    historical,
                    gzip_encoded=suffix in GZIP_EXTRACTION_SUFFIXES,
                )
                comparisons[suffix] = comparison
                stdout.write(
                    (
                        f"{suffix}\tsemantic_equal={comparison['semantic_equal']}\t"
                        f"physical_equal={comparison['physical_equal']}\n"
                    ).encode("utf-8")
                )
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
            stderr.write((error + "\n").encode("utf-8"))
    semantic_mismatch_count = sum(
        not value["semantic_equal"] for value in comparisons.values()
    )
    scientific_equivalence_all = (
        len(comparisons) == len(EXTRACTION_SUFFIXES)
        and semantic_mismatch_count == 0
    )
    if error is None and not scientific_equivalence_all:
        error = (
            "historical semantic comparison failed: "
            f"compared={len(comparisons)}/{len(EXTRACTION_SUFFIXES)}, "
            f"semantic_mismatches={semantic_mismatch_count}"
        )
        with stderr_path.open("ab") as stderr:
            stderr.write((error + "\n").encode("utf-8"))
    return {
        "stage": stage,
        "started_at_utc": started_utc,
        "ended_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - started,
        "command": [
            "internal:semantic-artifact-compare",
            str(extraction_dir.resolve()),
            str(historical_root.resolve()),
        ],
        "command_shell": "internal semantic comparison (gzip inputs use decompressed bytes)",
        "exit_code": 0 if error is None else 1,
        "all_pass": error is None and scientific_equivalence_all,
        "failure": error,
        "logs": {
            "stdout": file_identity(stdout_path),
            "stderr": file_identity(stderr_path),
        },
        "artifact_count": len(comparisons),
        "semantic_mismatch_count": semantic_mismatch_count,
        "physical_only_difference_count": sum(
            value["physical_only_difference"] for value in comparisons.values()
        ),
        "scientific_equivalence_all": scientific_equivalence_all,
        "comparisons": comparisons,
    }


def decompress_normalized_stage(
    python_dir: Path, normalized_dir: Path, logs_dir: Path
) -> dict[str, Any]:
    stage = "normalize_for_cpp"
    stdout_path = logs_dir / f"{stage}.stdout.log"
    stderr_path = logs_dir / f"{stage}.stderr.log"
    logs_dir.mkdir(parents=True, exist_ok=True)
    started_utc = utc_now()
    started = time.perf_counter()
    error = None
    inputs: dict[str, Any] = {}
    outputs: dict[str, Any] = {}
    command: list[str] = ["internal:gzip-stream-decompress"]
    try:
        if normalized_dir.exists() or normalized_dir.is_symlink():
            raise RunnerError(f"normalized directory already exists: {normalized_dir}")
        normalized_dir.mkdir(parents=True, exist_ok=False)
        for name in ("units", "constraints"):
            source = python_dir / f"{name}.tsv.gz"
            target = normalized_dir / f"{name}.tsv"
            inputs[source.name] = file_identity(source)
            command.extend([str(source.resolve()), str(target.resolve())])
            with gzip.open(source, "rb") as reader, target.open("xb") as writer:
                shutil.copyfileobj(reader, writer, length=1 << 20)
                writer.flush()
                os.fsync(writer.fileno())
            outputs[target.name] = file_identity(target)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    with stdout_path.open("xb") as stdout:
        if error is None:
            stdout.write(b"units.tsv and constraints.tsv decompressed\n")
    with stderr_path.open("xb") as stderr:
        if error is not None:
            stderr.write((error + "\n").encode("utf-8"))
    return {
        "stage": stage,
        "started_at_utc": started_utc,
        "ended_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - started,
        "command": command,
        "command_shell": shlex.join(command),
        "exit_code": 0 if error is None else 1,
        "all_pass": error is None,
        "failure": error,
        "inputs": inputs,
        "outputs": outputs,
        "logs": {
            "stdout": file_identity(stdout_path),
            "stderr": file_identity(stderr_path),
        },
    }


def read_tsv(path: Path, *, compressed: bool) -> list[dict[str, str]]:
    opener = gzip.open if compressed else Path.open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:  # type: ignore[arg-type]
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise RunnerError(f"TSV has no header: {path}")
            return list(reader)
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeError, csv.Error) as exc:
        raise RunnerError(f"cannot read TSV {path}: {exc}") from exc


def _nonnegative(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise RunnerError(f"{label} is not an integer: {value!r}") from exc
    if parsed < 0:
        raise RunnerError(f"{label} is negative: {parsed}")
    return parsed


def _require_columns(
    rows: Sequence[Mapping[str, str]], required: Iterable[str], label: str
) -> None:
    if not rows:
        return
    missing = set(required) - set(rows[0])
    if missing:
        raise RunnerError(f"{label} missing columns: {sorted(missing)}")


def chromosome_metrics(
    *,
    chrom: str,
    extraction_dir: Path,
    python_dir: Path,
    comparison_path: Path,
    historical: Mapping[str, Any] | None,
) -> dict[str, int]:
    prefix = f"{DATASET}.{chrom}"
    sites = read_tsv(extraction_dir / f"{prefix}.site_catalog.tsv.gz", compressed=True)
    units = read_tsv(python_dir / "units.tsv.gz", compressed=True)
    constraints = read_tsv(python_dir / "constraints.tsv.gz", compressed=True)
    blocks = read_tsv(python_dir / "blocks.tsv.gz", compressed=True)
    membership = read_tsv(python_dir / "membership.tsv.gz", compressed=True)
    dispositions = read_tsv(python_dir / "dispositions.tsv.gz", compressed=True)
    _require_columns(sites, ("site_index", "chrom", "pos1"), "site catalog")
    _require_columns(
        units, ("unit_id", "hp_family", "phase_set", "k"), "units"
    )
    _require_columns(
        constraints,
        (
            "unit_id",
            "constraint_id",
            "hp_family",
            "phase_set",
            "molecule_weight",
            "pattern_count",
        ),
        "constraints",
    )
    _require_columns(
        dispositions,
        (
            "unit_id",
            "constraint_id",
            "hp_family",
            "phase_set",
            "molecule_weight",
            "pattern_count",
            "disposition",
        ),
        "dispositions",
    )
    _require_columns(
        blocks, ("unit_id", "hp_family", "phase_set"), "blocks"
    )
    _require_columns(
        membership,
        ("unit_id", "hp_family", "phase_set", "site_index", "chrom", "pos1"),
        "membership",
    )

    catalog_by_index: dict[str, tuple[str, int]] = {}
    catalog_loci: set[tuple[str, int]] = set()
    for row in sites:
        site_index = row["site_index"]
        site_chrom = row["chrom"]
        pos1 = _nonnegative(row["pos1"], f"{chrom} catalog pos1")
        if pos1 == 0 or site_chrom != chrom:
            raise RunnerError(f"invalid catalog locus in {chrom}: {row}")
        locus = (site_chrom, pos1)
        if site_index in catalog_by_index:
            raise RunnerError(f"duplicate catalog site_index in {chrom}: {site_index}")
        if locus in catalog_loci:
            raise RunnerError(f"duplicate catalog (chrom,pos1) in {chrom}: {locus}")
        catalog_by_index[site_index] = locus
        catalog_loci.add(locus)

    active_loci: set[tuple[str, int]] = set()
    for row in membership:
        site_index = row["site_index"]
        locus = catalog_by_index.get(site_index)
        if locus is None:
            raise RunnerError(f"{chrom} membership references unknown site_index: {site_index}")
        observed_locus = (
            row["chrom"],
            _nonnegative(row["pos1"], f"{chrom} membership pos1"),
        )
        if observed_locus != locus:
            raise RunnerError(
                f"{chrom} membership/catalog locus mismatch for site_index={site_index}"
            )
        active_loci.add(locus)

    unit_by_id: dict[str, tuple[str, str]] = {}
    k_bins = {"eq1": 0, "2to8": 0, "9to12": 0, "gt12": 0}
    for row in units:
        unit_id = row["unit_id"]
        if unit_id in unit_by_id:
            raise RunnerError(f"duplicate unit_id in {chrom}: {unit_id}")
        unit_by_id[unit_id] = (row["hp_family"], row["phase_set"])
        k = _nonnegative(row["k"], f"{chrom} unit k")
        if k == 1:
            k_bins["eq1"] += 1
        elif k <= 8:
            k_bins["2to8"] += 1
        elif k <= 12:
            k_bins["9to12"] += 1
        else:
            k_bins["gt12"] += 1

    cross_ps = 0
    cross_hp = 0
    for table in (constraints, blocks, membership, dispositions):
        for row in table:
            route = unit_by_id.get(row.get("unit_id", ""))
            if route is None:
                raise RunnerError(
                    f"{chrom} output references unknown unit: {row.get('unit_id')}"
                )
            cross_hp += row.get("hp_family") != route[0]
            cross_ps += row.get("phase_set") != route[1]

    constraint_values: dict[tuple[str, str], tuple[int, int]] = {}
    for row in constraints:
        key = (row["unit_id"], row["constraint_id"])
        if key in constraint_values:
            raise RunnerError(f"duplicate constraint in {chrom}: {key}")
        constraint_values[key] = (
            _nonnegative(row["molecule_weight"], f"{chrom} molecule weight"),
            _nonnegative(row["pattern_count"], f"{chrom} pattern count"),
        )

    statuses = {
        "retained": [0, 0, 0],
        "cut": [0, 0, 0],
        "unavoidable_span_gt_max_block_size": [0, 0, 0],
    }
    disposition_keys: set[tuple[str, str]] = set()
    for row in dispositions:
        key = (row["unit_id"], row["constraint_id"])
        if key in disposition_keys or key not in constraint_values:
            raise RunnerError(f"invalid disposition constraint key in {chrom}: {key}")
        disposition_keys.add(key)
        weight = _nonnegative(row["molecule_weight"], f"{chrom} disposition weight")
        pattern_count = _nonnegative(
            row["pattern_count"], f"{chrom} disposition pattern count"
        )
        if (weight, pattern_count) != constraint_values[key]:
            raise RunnerError(f"constraint/disposition mass mismatch in {chrom}: {key}")
        status = row["disposition"]
        if status not in statuses:
            raise RunnerError(f"unknown disposition in {chrom}: {status}")
        statuses[status][0] += 1
        statuses[status][1] += weight
        statuses[status][2] += pattern_count
    if disposition_keys != set(constraint_values):
        raise RunnerError(f"constraint dispositions do not conserve keys in {chrom}")

    total_weight = sum(value[0] for value in constraint_values.values())
    total_patterns = sum(value[1] for value in constraint_values.values())
    if total_weight != sum(value[1] for value in statuses.values()):
        raise RunnerError(f"molecule weight is not conserved in {chrom}")
    if total_patterns != sum(value[2] for value in statuses.values()):
        raise RunnerError(f"pattern count is not conserved in {chrom}")

    comparison = strict_json_load(comparison_path)
    if not isinstance(comparison, Mapping) or comparison.get("all_pass") is not True:
        raise RunnerError(f"comparison receipt is not PASS: {comparison_path}")
    mismatch_count = _nonnegative(
        str(comparison.get("mismatch_count")), f"{chrom} mismatch_count"
    )

    historical_artifacts = 0
    historical_mismatches = 0
    historical_physical_only = 0
    if historical is not None:
        historical_artifacts = int(historical["artifact_count"])
        historical_mismatches = int(historical["semantic_mismatch_count"])
        historical_physical_only = int(
            historical["physical_only_difference_count"]
        )

    return {
        "S": len(catalog_loci),
        "units": len(units),
        "unit_memberships": len(membership),
        "unique_sites": len(active_loci),
        "k_eq_1_units": k_bins["eq1"],
        "k_2_to_8_units": k_bins["2to8"],
        "k_9_to_12_units": k_bins["9to12"],
        "k_gt_12_units": k_bins["gt12"],
        "blocks": len(blocks),
        "constraints": len(constraints),
        "pattern_count_total": total_patterns,
        "retained_constraints": statuses["retained"][0],
        "cut_constraints": statuses["cut"][0],
        "unavoidable_constraints": statuses[
            "unavoidable_span_gt_max_block_size"
        ][0],
        "molecule_weight_total": total_weight,
        "molecule_weight_retained": statuses["retained"][1],
        "molecule_weight_cut": statuses["cut"][1],
        "molecule_weight_unavoidable": statuses[
            "unavoidable_span_gt_max_block_size"
        ][1],
        "cross_ps_violations": cross_ps,
        "cross_hp_violations": cross_hp,
        "python_cpp_mismatch_count": mismatch_count,
        "historical_artifacts_compared": historical_artifacts,
        "historical_semantic_mismatches": historical_mismatches,
        "historical_physical_only_differences": historical_physical_only,
    }


def _fail_if_stage_failed(record: Mapping[str, Any]) -> None:
    if record.get("all_pass") is not True:
        raise RunnerError(str(record.get("failure") or "stage failed"))


def extraction_artifacts(extraction_dir: Path, chrom: str) -> tuple[Path, ...]:
    return tuple(
        extraction_dir / f"{DATASET}.{chrom}.{suffix}"
        for suffix in EXTRACTION_SUFFIXES
    )


def process_chromosome(
    *,
    chrom: str,
    args: argparse.Namespace,
    cpp_binary: Path,
    manifest_identity: Mapping[str, Any],
    tool_identities: Mapping[str, Any],
) -> dict[str, Any]:
    chrom_root = args.output_root / "chromosomes" / chrom
    chrom_root.mkdir(parents=True, exist_ok=False)
    logs_dir = chrom_root / "logs"
    stages: list[dict[str, Any]] = []
    metrics: dict[str, int] | None = None
    historical: dict[str, Any] | None = None
    manifest_checks: dict[str, Any] = {}
    failure = None
    try:
        manifest_before = file_identity(args.manifest)
        if manifest_before != dict(manifest_identity):
            raise RunnerError(f"{chrom} manifest identity changed before extraction")
        manifest_checks["before_extraction"] = manifest_before
        extraction_dir = chrom_root / "extraction"
        extract_command = [
            str(args.python),
            str(args.extractor),
            "--manifest",
            str(args.manifest),
            "--sample",
            DATASET,
            "--chrom",
            chrom,
            "--output-dir",
            str(extraction_dir),
            "--mapq-min",
            "20",
            "--baseq-min",
            "20",
            "--bridge-thresholds",
            "1,2,3,5",
            "--samtools-threads",
            "1",
            "--samtools",
            str(args.samtools),
        ]
        extraction = external_child_stage(
            stage="extract",
            command=extract_command,
            logs_dir=logs_dir,
            receipt_path=extraction_dir / "receipt.json",
            required_outputs=extraction_artifacts(extraction_dir, chrom),
            expected_output_parent=extraction_dir,
            require_receipt_sidecar=True,
        )
        try:
            extraction["extraction_contract_validation"] = validate_extraction_receipt(
                extraction_dir / "receipt.json",
                chrom=chrom,
                expected_manifest=manifest_identity,
                expected_samtools=args.samtools,
            )
        except Exception as exc:
            extraction["extraction_contract_validation"] = None
            extraction["all_pass"] = False
            extraction["failure"] = f"{type(exc).__name__}: {exc}"
        stages.append(extraction)
        _fail_if_stage_failed(extraction)

        manifest_after = file_identity(args.manifest)
        if manifest_after != dict(manifest_identity):
            raise RunnerError(f"{chrom} manifest identity changed during extraction")
        manifest_checks["after_extraction"] = manifest_after

        if args.historical_extraction_root is not None:
            historical = historical_comparison_stage(
                chrom=chrom,
                extraction_dir=extraction_dir,
                historical_root=args.historical_extraction_root,
                logs_dir=logs_dir,
            )
            stages.append(historical)
            _fail_if_stage_failed(historical)

        prefix = f"{DATASET}.{chrom}"
        python_dir = chrom_root / "python_partition"
        partition_command = [
            str(args.python),
            str(args.partitioner),
            "--dataset",
            DATASET,
            "--chrom",
            chrom,
            "--site-catalog",
            str(extraction_dir / f"{prefix}.site_catalog.tsv.gz"),
            "--site-component-membership",
            str(extraction_dir / f"{prefix}.site_component_membership.tsv.gz"),
            "--molecule-calls",
            str(extraction_dir / f"{prefix}.molecule_sparse_calls.tsv.gz"),
            "--output-dir",
            str(python_dir),
            "--threshold",
            "3",
            "--max-block-size",
            "12",
        ]
        partition = external_child_stage(
            stage="python_partition",
            command=partition_command,
            logs_dir=logs_dir,
            receipt_path=python_dir / "receipt.json",
            required_outputs=tuple(
                python_dir / name
                for name in (
                    "units.tsv.gz",
                    "constraints.tsv.gz",
                    "blocks.tsv.gz",
                    "membership.tsv.gz",
                    "dispositions.tsv.gz",
                )
            ),
            expected_output_parent=python_dir,
        )
        stages.append(partition)
        _fail_if_stage_failed(partition)

        normalized_dir = chrom_root / "normalized_cpp_inputs"
        normalized = decompress_normalized_stage(
            python_dir, normalized_dir, logs_dir
        )
        stages.append(normalized)
        _fail_if_stage_failed(normalized)

        cpp_dir = chrom_root / "cpp_partition"
        cpp_command = [
            str(cpp_binary),
            "--units",
            str(normalized_dir / "units.tsv"),
            "--constraints",
            str(normalized_dir / "constraints.tsv"),
            "--output-dir",
            str(cpp_dir),
            "--max-block-size",
            "12",
        ]
        cpp_stage = external_child_stage(
            stage="cpp_partition",
            command=cpp_command,
            logs_dir=logs_dir,
            receipt_path=cpp_dir / "summary.json",
            required_outputs=tuple(
                cpp_dir / name
                for name in (
                    "blocks.tsv",
                    "membership.tsv",
                    "dispositions.tsv",
                    "summary.json",
                )
            ),
        )
        stages.append(cpp_stage)
        _fail_if_stage_failed(cpp_stage)

        comparison_path = chrom_root / "comparison.json"
        compare_command = [
            str(args.python),
            str(args.comparator),
            "--python-dir",
            str(python_dir),
            "--cpp-input-units",
            str(normalized_dir / "units.tsv"),
            "--cpp-input-constraints",
            str(normalized_dir / "constraints.tsv"),
            "--cpp-dir",
            str(cpp_dir),
            "--output",
            str(comparison_path),
        ]
        comparison = external_child_stage(
            stage="compare_python_cpp",
            command=compare_command,
            logs_dir=logs_dir,
            receipt_path=comparison_path,
            required_outputs=(comparison_path,),
        )
        stages.append(comparison)
        _fail_if_stage_failed(comparison)

        metrics = chromosome_metrics(
            chrom=chrom,
            extraction_dir=extraction_dir,
            python_dir=python_dir,
            comparison_path=comparison_path,
            historical=historical,
        )
        if (
            metrics["cross_ps_violations"] != 0
            or metrics["cross_hp_violations"] != 0
            or metrics["python_cpp_mismatch_count"] != 0
        ):
            raise RunnerError(f"{chrom} final invariant count is nonzero: {metrics}")
    except Exception as exc:
        failure = {"type": type(exc).__name__, "message": str(exc)}

    all_pass = failure is None and metrics is not None and all(
        stage.get("all_pass") is True for stage in stages
    )
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.chromosome_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "task_type": "exploratory_pilot",
        "claim_status": "PARTIAL",
        "validation_evidence_eligible": False,
        "all_pass": all_pass,
        "sample": DATASET,
        "chrom": chrom,
        "manifest": dict(manifest_identity),
        "tools": dict(tool_identities),
        "stages": stages,
        "metrics": metrics,
        "historical_semantic_comparison": historical,
        "manifest_identity_checks": manifest_checks,
        "failure": failure,
        "claim_ceiling": (
            "Endpoint is extraction -> exact-PS segmentation -> C++ parity for "
            "HCC1395 local engineering evidence only. Topology inference, VAF "
            "ranking, partial-read likelihood adaptation, and unique cellular "
            "clone/evolutionary-tree claims are explicitly out of scope."
        ),
    }
    receipt_path = chrom_root / "stage_receipt.json"
    write_json_exclusive(receipt_path, receipt)
    write_sha256_sidecar(receipt_path)
    return {
        "chrom": chrom,
        "all_pass": all_pass,
        "metrics": metrics,
        "failure": failure,
        "receipt": file_identity(receipt_path),
    }


def write_chromosome_summary(
    path: Path, chromosomes: Sequence[str], results: Mapping[str, Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for chrom in chromosomes:
            result = results.get(chrom)
            row: dict[str, Any] = {field: "" for field in SUMMARY_FIELDS}
            row["chrom"] = chrom
            row["status"] = "PASS" if result and result["all_pass"] else "FAIL"
            if result and result.get("metrics"):
                row.update(result["metrics"])
            if result and result.get("receipt"):
                row["stage_receipt"] = result["receipt"]["path"]
                row["stage_receipt_sha256"] = result["receipt"]["sha256"]
            writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def aggregate_results(
    chromosomes: Sequence[str], results: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    metric_fields = SUMMARY_FIELDS[2:-2]
    totals = {field: 0 for field in metric_fields}
    results_with_metrics = [
        result for result in results.values() if result.get("metrics") is not None
    ]
    passed_results_with_metrics = [
        result for result in results_with_metrics if result.get("all_pass") is True
    ]
    for result in passed_results_with_metrics:
        metrics = result.get("metrics")
        assert metrics is not None
        for field in metric_fields:
            totals[field] += int(metrics[field])
    return {
        "chromosomes_requested": len(chromosomes),
        "chromosomes_passed": sum(
            bool(results.get(chrom, {}).get("all_pass")) for chrom in chromosomes
        ),
        "chromosomes_failed": sum(
            not bool(results.get(chrom, {}).get("all_pass"))
            for chrom in chromosomes
        ),
        "chromosomes_with_metrics": len(results_with_metrics),
        "chromosomes_included_in_totals": len(passed_results_with_metrics),
        "failed_chromosomes_with_metrics_excluded": (
            len(results_with_metrics) - len(passed_results_with_metrics)
        ),
        "aggregation_policy": "PASS_chromosomes_with_complete_metrics_only",
        "S": totals["S"],
        "units": totals["units"],
        "unit_memberships": totals["unit_memberships"],
        "unique_sites": totals["unique_sites"],
        "k_bins": {
            "k=1": totals["k_eq_1_units"],
            "k=2-8": totals["k_2_to_8_units"],
            "k=9-12": totals["k_9_to_12_units"],
            "k>12": totals["k_gt_12_units"],
        },
        "blocks": totals["blocks"],
        "constraints": totals["constraints"],
        "pattern_count_total": totals["pattern_count_total"],
        "constraint_dispositions": {
            "retained": totals["retained_constraints"],
            "cut": totals["cut_constraints"],
            "unavoidable": totals["unavoidable_constraints"],
        },
        "molecule_weights": {
            "total": totals["molecule_weight_total"],
            "retained": totals["molecule_weight_retained"],
            "cut": totals["molecule_weight_cut"],
            "unavoidable": totals["molecule_weight_unavoidable"],
        },
        "cross_ps_violations": totals["cross_ps_violations"],
        "cross_hp_violations": totals["cross_hp_violations"],
        "python_cpp_mismatch_count": totals["python_cpp_mismatch_count"],
        "historical_comparison": {
            "artifacts_compared": totals["historical_artifacts_compared"],
            "semantic_mismatches": totals["historical_semantic_mismatches"],
            "physical_only_differences": totals[
                "historical_physical_only_differences"
            ],
        },
    }


def normalize_args(args: argparse.Namespace) -> argparse.Namespace:
    args.manifest = require_source(args.manifest, "manifest")
    args.input_verifier = require_source(args.input_verifier, "input verifier")
    args.extractor = require_source(args.extractor, "collapsing extractor")
    if "collapsing" not in args.extractor.name and args.extractor == DEFAULT_EXTRACTOR.resolve():
        raise RunnerError("default extractor must be the identity-safe collapsing extractor")
    args.partitioner = require_source(args.partitioner, "Python partitioner")
    args.cpp_source = require_source(args.cpp_source, "C++ source")
    args.comparator = require_source(args.comparator, "comparator")
    args.python = resolve_executable(args.python, "Python")
    args.cxx = resolve_executable(args.cxx, "C++ compiler")
    args.samtools = resolve_executable(args.samtools, "samtools")
    if args.historical_extraction_root is not None:
        args.historical_extraction_root = (
            args.historical_extraction_root.expanduser().resolve()
        )
        if not args.historical_extraction_root.is_dir():
            raise RunnerError(
                "historical extraction root is not a directory: "
                f"{args.historical_extraction_root}"
            )
    return args


def execute(args: argparse.Namespace) -> dict[str, Any]:
    args = normalize_args(args)
    tool_identities = {
        "runner": file_identity(Path(__file__).resolve()),
        "input_verifier": file_identity(args.input_verifier),
        "collapsing_extractor": file_identity(args.extractor),
        "python_partitioner": file_identity(args.partitioner),
        "cpp_source": file_identity(args.cpp_source),
        "comparator": file_identity(args.comparator),
        "python": file_identity(args.python),
        "cxx": file_identity(args.cxx),
        "samtools": file_identity(args.samtools),
    }
    manifest_identity = file_identity(args.manifest)
    args.output_root = prepare_output_root(args.output_root)
    started_at = utc_now()
    global_stages: list[dict[str, Any]] = []
    chromosome_results: dict[str, dict[str, Any]] = {}
    failure = None
    pre_binding_validation: dict[str, Any] | None = None
    cpp_binary = args.output_root / "bin" / "exact_ps_k12_partition"

    try:
        verifier_receipt = args.output_root / "input_binding_receipt.pre.json"
        verifier_command = [
            str(args.python),
            str(args.input_verifier),
            "--manifest",
            str(args.manifest),
            "--output",
            str(verifier_receipt),
            "--samtools",
            str(args.samtools),
        ]
        input_stage = external_child_stage(
            stage="input_verifier_pre",
            command=verifier_command,
            logs_dir=args.output_root / "logs",
            receipt_path=verifier_receipt,
            required_outputs=(verifier_receipt,),
        )
        attach_input_binding_validation(
            input_stage,
            verifier_receipt,
            expected_manifest=manifest_identity,
            expected_samtools=args.samtools,
        )
        global_stages.append(input_stage)
        _fail_if_stage_failed(input_stage)
        pre_binding_validation = input_stage["input_binding_validation"]

        cpp_binary.parent.mkdir(parents=True, exist_ok=False)
        compile_command = [
            str(args.cxx),
            "-std=c++17",
            "-O2",
            "-Wall",
            "-Wextra",
            "-Werror",
            str(args.cpp_source),
            "-o",
            str(cpp_binary),
        ]
        compile_stage = simple_external_stage(
            stage="compile_cpp",
            command=compile_command,
            logs_dir=args.output_root / "logs",
            required_outputs=(cpp_binary,),
        )
        global_stages.append(compile_stage)
        _fail_if_stage_failed(compile_stage)

        if args.workers == 1:
            for chrom in args.chromosomes:
                chromosome_results[chrom] = process_chromosome(
                    chrom=chrom,
                    args=args,
                    cpp_binary=cpp_binary,
                    manifest_identity=manifest_identity,
                    tool_identities=tool_identities,
                )
        else:
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                future_to_chrom = {
                    executor.submit(
                        process_chromosome,
                        chrom=chrom,
                        args=args,
                        cpp_binary=cpp_binary,
                        manifest_identity=manifest_identity,
                        tool_identities=tool_identities,
                    ): chrom
                    for chrom in args.chromosomes
                }
                for future in as_completed(future_to_chrom):
                    chrom = future_to_chrom[future]
                    chromosome_results[chrom] = future.result()

        post_verifier_receipt = args.output_root / "input_binding_receipt.post.json"
        post_verifier_command = [
            str(args.python),
            str(args.input_verifier),
            "--manifest",
            str(args.manifest),
            "--output",
            str(post_verifier_receipt),
            "--samtools",
            str(args.samtools),
        ]
        post_input_stage = external_child_stage(
            stage="input_verifier_post",
            command=post_verifier_command,
            logs_dir=args.output_root / "logs",
            receipt_path=post_verifier_receipt,
            required_outputs=(post_verifier_receipt,),
        )
        attach_input_binding_validation(
            post_input_stage,
            post_verifier_receipt,
            expected_manifest=manifest_identity,
            expected_samtools=args.samtools,
        )
        post_validation = post_input_stage.get("input_binding_validation")
        if (
            post_input_stage.get("all_pass") is True
            and pre_binding_validation is not None
            and isinstance(post_validation, Mapping)
            and post_validation.get("binding_snapshot_sha256")
            != pre_binding_validation.get("binding_snapshot_sha256")
        ):
            post_input_stage["all_pass"] = False
            post_input_stage["failure"] = (
                "input binding snapshot changed between pre-run and post-run verification"
            )
        global_stages.append(post_input_stage)
        _fail_if_stage_failed(post_input_stage)
    except (Exception, KeyboardInterrupt) as exc:
        failure = {"type": type(exc).__name__, "message": str(exc)}

    summary_path = args.output_root / "chromosome_summary.tsv"
    write_chromosome_summary(summary_path, args.chromosomes, chromosome_results)
    aggregate = aggregate_results(args.chromosomes, chromosome_results)
    all_pass = (
        failure is None
        and len(chromosome_results) == len(args.chromosomes)
        and all(result["all_pass"] for result in chromosome_results.values())
        and all(stage.get("all_pass") is True for stage in global_stages)
    )
    if not all_pass and failure is None:
        failed = [
            chrom
            for chrom in args.chromosomes
            if not chromosome_results.get(chrom, {}).get("all_pass")
        ]
        failure = {
            "type": "ChromosomeStageFailure",
            "message": f"chromosome stage failed: {','.join(failed)}",
        }

    invocation = getattr(
        args,
        "invocation",
        [str(Path(__file__).resolve()), *sys.argv[1:]],
    )
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.run_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "started_at_utc": started_at,
        "task_type": "exploratory_pilot",
        "claim_status": "PARTIAL",
        "validation_evidence_eligible": False,
        "all_pass": all_pass,
        "sample": DATASET,
        "scope": {
            "chromosomes": list(args.chromosomes),
            "hcc1395_only": True,
            "autosomes_only": True,
            "workers": args.workers,
        },
        "invocation": list(invocation),
        "manifest": manifest_identity,
        "sources": tool_identities,
        "parameters": {
            "threshold": 3,
            "max_block_size": 12,
            "mapq_min": 20,
            "baseq_min": 20,
            "bridge_thresholds": [1, 2, 3, 5],
            "samtools_threads_per_extractor": 1,
            "resume": False,
            "historical_extraction_root": (
                None
                if args.historical_extraction_root is None
                else str(args.historical_extraction_root)
            ),
            "historical_scientific_comparison": (
                "BED direct bytes; .gz artifacts decompressed semantic bytes; "
                "physical gzip SHA-256 is diagnostic only"
            ),
            "pipeline_endpoint": (
                "extraction -> exact-PS segmentation -> C++ parity; excludes "
                "topology, VAF, and partial-read likelihood adapter"
            ),
        },
        "global_stages": global_stages,
        "chromosomes": [
            chromosome_results[chrom]
            for chrom in args.chromosomes
            if chrom in chromosome_results
        ],
        "aggregate": aggregate,
        "outputs": {"chromosome_summary": file_identity(summary_path)},
        "failure": failure,
        "claim_ceiling": (
            "Endpoint is extraction -> exact-PS segmentation -> C++ parity for "
            "the requested HCC1395 autosomes. This PARTIAL receipt excludes "
            "topology inference, VAF ranking, and the partial-read likelihood "
            "adapter; it is neither seven-sample validation nor proof of a "
            "unique true clone tree."
        ),
    }
    receipt_path = args.output_root / "run_receipt.json"
    write_json_exclusive(receipt_path, receipt)
    write_sha256_sidecar(receipt_path)

    if all_pass:
        marker = {
            "schema_name": f"{SCHEMA_NAME}.success_marker",
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": utc_now(),
            "all_pass": True,
            "claim_status": "PARTIAL",
            "sample": DATASET,
            "chromosomes": list(args.chromosomes),
            "run_receipt": file_identity(receipt_path),
        }
        write_json_exclusive(args.output_root / "_SUCCESS", marker)
    return receipt


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--chromosomes",
        type=parse_chromosomes,
        default=AUTOSOMES,
        help="comma-separated subset of chr1-chr22; default is chr1-chr22",
    )
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1)
    parser.add_argument("--historical-extraction-root", type=Path)

    # Explicit overrides keep the orchestration regression-testable without
    # weakening the fixed scientific parameters or HCC1395-only scope.
    parser.add_argument("--input-verifier", type=Path, default=DEFAULT_VERIFIER)
    parser.add_argument("--extractor", type=Path, default=DEFAULT_EXTRACTOR)
    parser.add_argument("--partitioner", type=Path, default=DEFAULT_PARTITIONER)
    parser.add_argument("--cpp-source", type=Path, default=DEFAULT_CPP_SOURCE)
    parser.add_argument("--comparator", type=Path, default=DEFAULT_COMPARATOR)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args(argv)
    args.invocation = [str(Path(__file__).resolve()), *(argv or sys.argv[1:])]
    return args


def main(argv: Sequence[str] | None = None) -> int:
    try:
        receipt = execute(parse_args(argv))
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    summary = {
        "all_pass": receipt["all_pass"],
        "sample": receipt["sample"],
        "chromosomes": receipt["scope"]["chromosomes"],
        "aggregate": receipt["aggregate"],
    }
    stream = sys.stdout if receipt["all_pass"] else sys.stderr
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True), file=stream)
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
