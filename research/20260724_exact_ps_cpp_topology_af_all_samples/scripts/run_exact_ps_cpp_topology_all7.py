#!/usr/bin/env python3
"""Run the auditable exact-PS k<=12 partition-to-topology cohort pipeline.

The runner reuses immutable extraction and strict-region artifacts.  It creates
new per-chromosome partition outputs, verifies Python/C++ partition parity,
builds the compact MLHP input, and invokes an external C++ topology binary.

Resume is fail closed: an existing chromosome is accepted only when its PASS
receipt still binds every declared artifact identity.  A partial chromosome
directory is never overwritten or silently repaired.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
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
AUTOSOMES = tuple(f"chr{index}" for index in range(1, 23))
CANONICAL_DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)

SCHEMA_NAME = "intersubmod.exact_ps_cpp_topology_cohort"
SCHEMA_VERSION = "1.0.0"
MANIFEST_SCHEMA = "intersubmod.exact_ps_cpp_topology_input_manifest"
MANIFEST_VERSION = "1.0.0"
EXTRACTION_SCHEMA = "intersubmod.lossless_read_linkage_chromosome_receipt"
EXTRACTION_VERSION = "1.3.0"
STRICT_SCHEMA = "intersubmod.strict_exact_ps_hp_regions"
STRICT_VERSION = "1.1.0"
PARTITION_SCHEMA = "intersubmod.exact_ps_k12_partition"
PARTITION_VERSION = "0.1.0"
CPP_PARTITION_SCHEMA = "intersubmod.exact_ps_k12_cpp_partition_summary"
CPP_PARTITION_VERSION = "1.0.0"
COMPARISON_SCHEMA = "intersubmod.exact_ps_k12_python_cpp_comparison"
COMPARISON_VERSION = "1.0.0"
MLHP_RECEIPT_SCHEMA = "intersubmod.exact_ps_layered_tree_input.receipt"

DEFAULT_MANIFEST = TOPIC_ROOT / "input_contract" / "all7_exact_ps_inputs.local.json"
DEFAULT_PARTITIONER = (
    REPO_ROOT
    / "research"
    / "20260722_exact_ps_k12_hcc1395_pilot"
    / "scripts"
    / "exact_ps_k12_partition.py"
)
DEFAULT_COMPARATOR = (
    REPO_ROOT
    / "research"
    / "20260722_exact_ps_k12_hcc1395_pilot"
    / "scripts"
    / "compare_python_cpp_partitions.py"
)
DEFAULT_ADAPTER = (
    REPO_ROOT
    / "docs"
    / "methodology"
    / "_assets"
    / "20260627_subclone_4axis_teaching"
    / "scripts"
    / "exact_ps_partition_to_mlhp.py"
)
DEFAULT_PARTITION_BINARY = REPO_ROOT / "build" / "bin" / "exact_ps_partition"

PYTHON_OUTPUTS = (
    "units.tsv.gz",
    "constraints.tsv.gz",
    "blocks.tsv.gz",
    "membership.tsv.gz",
    "dispositions.tsv.gz",
)
CPP_OUTPUTS = (
    "blocks.tsv",
    "membership.tsv",
    "dispositions.tsv",
    "summary.json",
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


def identity(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise RunnerError(f"not a regular file: {path}")
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise RunnerError(f"missing JSON file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RunnerError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RunnerError(f"JSON root is not an object: {path}")
    return value


def write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        indent=2,
    ) + "\n"
    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise RunnerError(f"refusing to overwrite output: {path}") from exc


def write_sidecar_exclusive(path: Path) -> Path:
    sidecar = path.with_name(f"{path.name}.sha256")
    try:
        with sidecar.open("x", encoding="ascii") as handle:
            handle.write(f"{sha256_path(path)}  {path.name}\n")
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise RunnerError(f"refusing to overwrite checksum: {sidecar}") from exc
    return sidecar


def verify_sidecar(path: Path) -> dict[str, Any]:
    sidecar = path.with_name(f"{path.name}.sha256")
    if not sidecar.is_file():
        raise RunnerError(f"missing SHA-256 sidecar: {sidecar}")
    fields = sidecar.read_text(encoding="ascii").strip().split()
    if (
        len(fields) != 2
        or fields[1] != path.name
        or fields[0] != sha256_path(path)
    ):
        raise RunnerError(f"invalid SHA-256 sidecar: {sidecar}")
    return identity(sidecar)


def identity_matches(spec: Mapping[str, Any], path: Path) -> bool:
    try:
        observed = identity(path)
    except RunnerError:
        return False
    size = spec.get("size_bytes", spec.get("bytes"))
    return (
        str(Path(str(spec.get("path", ""))).resolve()) == observed["path"]
        and size == observed["size_bytes"]
        and spec.get("sha256") == observed["sha256"]
    )


def verify_identity(spec: Mapping[str, Any], label: str) -> Path:
    path = Path(str(spec.get("path", "")))
    if not identity_matches(spec, path):
        raise RunnerError(f"{label} identity mismatch: {path}")
    return path.resolve()


def require_all_true(checks: object, label: str) -> int:
    if not isinstance(checks, Mapping) or not checks:
        raise RunnerError(f"{label} checks must be a non-empty object")
    failed = sorted(key for key, value in checks.items() if value is not True)
    if failed:
        raise RunnerError(f"{label} failed checks: {','.join(failed)}")
    return len(checks)


def parse_csv(value: str, allowed: Sequence[str], label: str) -> tuple[str, ...]:
    tokens = tuple(token.strip() for token in value.split(",") if token.strip())
    if not tokens:
        raise argparse.ArgumentTypeError(f"{label} must not be empty")
    if len(tokens) != len(set(tokens)):
        raise argparse.ArgumentTypeError(f"{label} contains duplicates")
    invalid = [token for token in tokens if token not in allowed]
    if invalid:
        raise argparse.ArgumentTypeError(
            f"invalid {label}: {','.join(invalid)}"
        )
    return tokens


def resolve_executable(value: str | Path, label: str) -> Path:
    raw = str(value)
    candidate = Path(raw)
    if candidate.name == raw:
        discovered = shutil.which(raw)
        resolved = Path(discovered).resolve() if discovered else candidate.resolve()
    else:
        resolved = candidate.expanduser().resolve()
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise RunnerError(f"{label} is not executable: {resolved}")
    return resolved


def resolve_regular_file(value: str | Path, label: str) -> Path:
    resolved = Path(value).expanduser().resolve()
    if not resolved.is_file():
        raise RunnerError(f"{label} is missing: {resolved}")
    return resolved


def load_manifest(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    document = read_json(path)
    if (
        document.get("schema_name") != MANIFEST_SCHEMA
        or document.get("schema_version") != MANIFEST_VERSION
    ):
        raise RunnerError("input manifest schema mismatch")
    order = document.get("dataset_order")
    samples = document.get("samples")
    if tuple(order or ()) != CANONICAL_DATASETS:
        raise RunnerError("input manifest dataset_order is not the canonical seven")
    if not isinstance(samples, Mapping) or set(samples) != set(CANONICAL_DATASETS):
        raise RunnerError("input manifest must define exactly the canonical seven datasets")
    normalized: dict[str, Any] = {}
    for sample in CANONICAL_DATASETS:
        row = samples[sample]
        if not isinstance(row, Mapping):
            raise RunnerError(f"manifest sample row is not an object: {sample}")
        extraction_root = Path(str(row.get("extraction_root", ""))).expanduser()
        strict_root = Path(str(row.get("strict_root", ""))).expanduser()
        normalized[sample] = {
            "extraction_root": extraction_root,
            "strict_root": strict_root,
        }
    return document, normalized


def declared_output(
    receipt: Mapping[str, Any],
    key: str,
    expected_path: Path,
    label: str,
) -> dict[str, Any]:
    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping):
        raise RunnerError(f"{label} outputs are missing")
    spec = outputs.get(key)
    if not isinstance(spec, Mapping) or not identity_matches(spec, expected_path):
        raise RunnerError(f"{label} output identity mismatch: {key}")
    return identity(expected_path)


def validate_upstream(
    sample: str,
    chrom: str,
    extraction_dir: Path,
    strict_dir: Path,
) -> dict[str, Any]:
    extraction_receipt_path = extraction_dir / "receipt.json"
    extraction_receipt = read_json(extraction_receipt_path)
    if (
        extraction_receipt.get("schema_name") != EXTRACTION_SCHEMA
        or extraction_receipt.get("schema_version") != EXTRACTION_VERSION
        or extraction_receipt.get("all_pass") is not True
    ):
        raise RunnerError(f"extraction receipt contract mismatch: {extraction_receipt_path}")
    extraction_scope = extraction_receipt.get("scope")
    if not isinstance(extraction_scope, Mapping) or (
        extraction_scope.get("dataset"),
        extraction_scope.get("chrom"),
    ) != (sample, chrom):
        raise RunnerError(f"extraction receipt scope mismatch: {extraction_receipt_path}")
    extraction_checks = require_all_true(
        extraction_receipt.get("checks"), "extraction receipt"
    )
    extraction_sidecar = verify_sidecar(extraction_receipt_path)

    strict_receipt_path = strict_dir / "receipt.json"
    strict_receipt = read_json(strict_receipt_path)
    if (
        strict_receipt.get("schema_name") != STRICT_SCHEMA
        or strict_receipt.get("schema_version") != STRICT_VERSION
        or strict_receipt.get("all_pass") is not True
    ):
        raise RunnerError(f"strict receipt contract mismatch: {strict_receipt_path}")
    strict_scope = strict_receipt.get("scope")
    if not isinstance(strict_scope, Mapping) or (
        strict_scope.get("dataset"),
        strict_scope.get("chrom"),
    ) != (sample, chrom):
        raise RunnerError(f"strict receipt scope mismatch: {strict_receipt_path}")
    strict_checks = require_all_true(strict_receipt.get("checks"), "strict receipt")
    strict_sidecar = verify_sidecar(strict_receipt_path)

    site_catalog = extraction_dir / f"{sample}.{chrom}.site_catalog.tsv.gz"
    molecule_calls = extraction_dir / f"{sample}.{chrom}.molecule_sparse_calls.tsv.gz"
    extraction_membership = (
        extraction_dir / f"{sample}.{chrom}.site_component_membership.tsv.gz"
    )
    strict_membership = (
        strict_dir / f"{sample}.{chrom}.site_component_membership.tsv.gz"
    )
    for path in (site_catalog, molecule_calls, extraction_membership, strict_membership):
        if not path.is_file():
            raise RunnerError(f"required upstream artifact is missing: {path}")

    declared_output(
        extraction_receipt,
        site_catalog.name,
        site_catalog,
        "extraction receipt",
    )
    declared_output(
        extraction_receipt,
        molecule_calls.name,
        molecule_calls,
        "extraction receipt",
    )
    declared_output(
        extraction_receipt,
        extraction_membership.name,
        extraction_membership,
        "extraction receipt",
    )
    declared_output(
        strict_receipt,
        "membership",
        strict_membership,
        "strict receipt",
    )
    return {
        "extraction_receipt": identity(extraction_receipt_path),
        "extraction_receipt_sidecar": extraction_sidecar,
        "strict_receipt": identity(strict_receipt_path),
        "strict_receipt_sidecar": strict_sidecar,
        "site_catalog": identity(site_catalog),
        "molecule_calls": identity(molecule_calls),
        "strict_membership": identity(strict_membership),
        "checks": extraction_checks + strict_checks,
    }


def ensure_symlink(path: Path, target: Path, *, resume: bool) -> None:
    expected = target.resolve()
    if os.path.lexists(path):
        if not resume or not path.is_symlink() or path.resolve() != expected:
            raise RunnerError(f"existing symlink target is not resumable: {path}")
        return
    os.symlink(str(expected), path)


def run_command(
    *,
    stage: str,
    command: Sequence[str],
    logs_dir: Path,
) -> dict[str, Any]:
    logs_dir.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while True:
        stem = stage if attempt == 1 else f"{stage}.attempt{attempt}"
        stdout_path = logs_dir / f"{stem}.stdout.log"
        stderr_path = logs_dir / f"{stem}.stderr.log"
        stage_record_path = logs_dir / f"{stem}.stage.json"
        if not any(
            path.exists() for path in (stdout_path, stderr_path, stage_record_path)
        ):
            break
        attempt += 1
    started_utc = utc_now()
    started = time.perf_counter()
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        completed = subprocess.run(
            list(command),
            cwd=str(REPO_ROOT),
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            check=False,
        )
    record = {
        "stage": stage,
        "attempt": attempt,
        "started_at_utc": started_utc,
        "ended_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - started,
        "command": list(command),
        "command_shell": shlex.join(command),
        "exit_code": completed.returncode,
        "logs": {
            "stdout": identity(stdout_path),
            "stderr": identity(stderr_path),
        },
    }
    write_json_exclusive(stage_record_path, record)
    if completed.returncode != 0:
        raise RunnerError(
            f"stage failed with exit code {completed.returncode}: {stage}; "
            f"stderr={stderr_path}"
        )
    return record


def decompress_exclusive(source: Path, target: Path) -> dict[str, Any]:
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        with gzip.open(source, "rb") as input_handle, target.open("xb") as output:
            shutil.copyfileobj(input_handle, output, length=1 << 20)
            output.flush()
            os.fsync(output.fileno())
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise RunnerError(f"cannot decompress {source}: {exc}") from exc
    return identity(target)


def validate_partition_receipt(
    path: Path,
    sample: str,
    chrom: str,
) -> tuple[dict[str, Any], int]:
    receipt = read_json(path)
    if (
        receipt.get("schema_name") != PARTITION_SCHEMA
        or receipt.get("schema_version") != PARTITION_VERSION
        or receipt.get("all_pass") is not True
    ):
        raise RunnerError(f"Python partition receipt contract mismatch: {path}")
    scope = receipt.get("scope")
    if not isinstance(scope, Mapping) or (
        scope.get("dataset"),
        scope.get("chrom"),
    ) != (sample, chrom):
        raise RunnerError(f"Python partition scope mismatch: {path}")
    checks = require_all_true(receipt.get("checks"), "Python partition")
    for name in PYTHON_OUTPUTS:
        declared_output(receipt, name, path.parent / name, "Python partition")
    return receipt, checks


def validate_cpp_summary(path: Path) -> tuple[dict[str, Any], int]:
    summary = read_json(path)
    if (
        summary.get("schema_name") != CPP_PARTITION_SCHEMA
        or summary.get("schema_version") != CPP_PARTITION_VERSION
        or summary.get("all_pass") is not True
    ):
        raise RunnerError(f"C++ partition summary contract mismatch: {path}")
    return summary, require_all_true(summary.get("checks"), "C++ partition")


def validate_comparison(path: Path) -> tuple[dict[str, Any], int]:
    comparison = read_json(path)
    if (
        comparison.get("schema_name") != COMPARISON_SCHEMA
        or comparison.get("schema_version") != COMPARISON_VERSION
        or comparison.get("all_pass") is not True
        or comparison.get("mismatch_count") != 0
    ):
        raise RunnerError(f"Python/C++ comparison did not pass: {path}")
    return comparison, require_all_true(
        comparison.get("checks"), "Python/C++ comparison"
    )


def bound_artifacts(paths: Iterable[Path]) -> list[dict[str, Any]]:
    return [identity(path) for path in paths]


def validate_chromosome_resume(
    receipt_path: Path,
    *,
    sample: str,
    chrom: str,
    extraction_dir: Path,
    strict_dir: Path,
) -> dict[str, Any]:
    receipt = read_json(receipt_path)
    if (
        receipt.get("schema_name") != f"{SCHEMA_NAME}.chromosome_receipt"
        or receipt.get("schema_version") != SCHEMA_VERSION
        or receipt.get("all_pass") is not True
        or receipt.get("sample") != sample
        or receipt.get("chrom") != chrom
    ):
        raise RunnerError(f"existing chromosome receipt is not resumable: {receipt_path}")
    require_all_true(receipt.get("checks"), "chromosome resume")
    artifacts = receipt.get("bound_artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise RunnerError("chromosome resume receipt lacks bound artifact identities")
    for index, spec in enumerate(artifacts):
        if not isinstance(spec, Mapping):
            raise RunnerError(f"invalid bound artifact identity at index {index}")
        verify_identity(spec, f"chromosome bound artifact {index}")
    chrom_root = receipt_path.parent
    for name, target in (
        ("extraction", extraction_dir),
        ("strict_regions", strict_dir),
    ):
        link = chrom_root / name
        if not link.is_symlink() or link.resolve() != target.resolve():
            raise RunnerError(f"chromosome source symlink mismatch: {link}")
    verify_sidecar(receipt_path)
    return receipt


def process_chromosome(
    *,
    sample: str,
    chrom: str,
    source: Mapping[str, Path],
    sample_root: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    extraction_dir = (
        Path(source["extraction_root"]) / "chromosomes" / chrom / "extraction"
    ).resolve()
    strict_dir = (
        Path(source["strict_root"]) / "chromosomes" / chrom
    ).resolve()
    upstream = validate_upstream(sample, chrom, extraction_dir, strict_dir)
    chrom_root = sample_root / "chromosomes" / chrom
    receipt_path = chrom_root / "stage_receipt.json"
    if chrom_root.exists():
        if not args.resume:
            raise RunnerError(f"chromosome output already exists: {chrom_root}")
        return validate_chromosome_resume(
            receipt_path,
            sample=sample,
            chrom=chrom,
            extraction_dir=extraction_dir,
            strict_dir=strict_dir,
        )

    chrom_root.mkdir(parents=True, exist_ok=False)
    ensure_symlink(chrom_root / "extraction", extraction_dir, resume=False)
    ensure_symlink(chrom_root / "strict_regions", strict_dir, resume=False)
    logs_dir = chrom_root / "logs"
    stages: list[dict[str, Any]] = []

    prefix = f"{sample}.{chrom}"
    python_dir = chrom_root / "python_partition"
    partition_command = [
        str(args.python),
        str(args.partitioner),
        "--dataset",
        sample,
        "--chrom",
        chrom,
        "--site-catalog",
        str(extraction_dir / f"{prefix}.site_catalog.tsv.gz"),
        "--site-component-membership",
        str(strict_dir / f"{prefix}.site_component_membership.tsv.gz"),
        "--molecule-calls",
        str(extraction_dir / f"{prefix}.molecule_sparse_calls.tsv.gz"),
        "--output-dir",
        str(python_dir),
        "--threshold",
        "3",
        "--max-block-size",
        "12",
    ]
    stages.append(
        run_command(
            stage="python_partition",
            command=partition_command,
            logs_dir=logs_dir,
        )
    )
    partition_receipt, python_checks = validate_partition_receipt(
        python_dir / "receipt.json", sample, chrom
    )
    partition_membership = (partition_receipt.get("inputs") or {}).get(
        "site_component_membership"
    )
    if not isinstance(partition_membership, Mapping) or not identity_matches(
        partition_membership,
        strict_dir / f"{prefix}.site_component_membership.tsv.gz",
    ):
        raise RunnerError("partition is not bound to strict endpoint membership")

    normalized_dir = chrom_root / "normalized_cpp_inputs"
    units_plain = decompress_exclusive(
        python_dir / "units.tsv.gz", normalized_dir / "units.tsv"
    )
    constraints_plain = decompress_exclusive(
        python_dir / "constraints.tsv.gz", normalized_dir / "constraints.tsv"
    )
    cpp_dir = chrom_root / "cpp_partition"
    cpp_command = [
        str(args.partition_binary),
        "--units",
        str(normalized_dir / "units.tsv"),
        "--constraints",
        str(normalized_dir / "constraints.tsv"),
        "--output-dir",
        str(cpp_dir),
        "--max-block-size",
        "12",
    ]
    stages.append(
        run_command(
            stage="cpp_partition",
            command=cpp_command,
            logs_dir=logs_dir,
        )
    )
    _, cpp_checks = validate_cpp_summary(cpp_dir / "summary.json")

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
    stages.append(
        run_command(
            stage="compare_python_cpp",
            command=compare_command,
            logs_dir=logs_dir,
        )
    )
    comparison, comparison_checks = validate_comparison(comparison_path)

    output_paths = [
        python_dir / "receipt.json",
        *(python_dir / name for name in PYTHON_OUTPUTS),
        Path(units_plain["path"]),
        Path(constraints_plain["path"]),
        *(cpp_dir / name for name in CPP_OUTPUTS),
        comparison_path,
    ]
    source_paths = [
        Path(value["path"])
        for key, value in upstream.items()
        if key != "checks" and isinstance(value, Mapping) and "path" in value
    ]
    all_checks = upstream["checks"] + python_checks + cpp_checks + comparison_checks
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.chromosome_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "sample": sample,
        "chrom": chrom,
        "parameters": {
            "threshold": 3,
            "max_active_bits": 12,
            "partition_membership_source": "strict_endpoint_regions",
        },
        "checks": {
            "upstream_receipts_pass": True,
            "partition_uses_strict_membership": True,
            "python_partition_pass": True,
            "cpp_partition_pass": True,
            "python_cpp_mismatch_zero": comparison.get("mismatch_count") == 0,
            "all_child_checks_pass": True,
        },
        "check_count": all_checks,
        "stages": stages,
        "bound_artifacts": bound_artifacts([*source_paths, *output_paths]),
        "source_links": {
            "extraction": str(extraction_dir),
            "strict_regions": str(strict_dir),
        },
        "counts": {
            "python": comparison.get("counts", {}).get("python_rows", {}),
            "cpp": comparison.get("counts", {}).get("cpp_rows", {}),
        },
    }
    write_json_exclusive(receipt_path, receipt)
    write_sidecar_exclusive(receipt_path)
    return receipt


def validate_sample_partition_receipt(
    path: Path,
    *,
    sample: str,
    expected_chroms: Sequence[str],
) -> dict[str, Any]:
    receipt = read_json(path)
    if (
        receipt.get("schema_name") != f"{SCHEMA_NAME}.sample_partition_receipt"
        or receipt.get("schema_version") != SCHEMA_VERSION
        or receipt.get("all_pass") is not True
        or receipt.get("sample") != sample
    ):
        raise RunnerError(f"sample partition receipt mismatch: {path}")
    scope = receipt.get("scope")
    if not isinstance(scope, Mapping) or scope.get("chromosomes") != list(expected_chroms):
        raise RunnerError(f"sample partition receipt scope mismatch: {path}")
    require_all_true(receipt.get("checks"), "sample partition")
    chrom_receipts = receipt.get("chromosome_receipts")
    if not isinstance(chrom_receipts, Mapping) or set(chrom_receipts) != set(expected_chroms):
        raise RunnerError("sample partition chromosome receipt grid mismatch")
    for chrom, spec in chrom_receipts.items():
        if not isinstance(spec, Mapping):
            raise RunnerError(f"invalid chromosome receipt identity: {chrom}")
        verify_identity(spec, f"sample chromosome receipt {chrom}")
    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping) or not inputs:
        raise RunnerError("sample partition receipt input identities are missing")
    for name, spec in inputs.items():
        if not isinstance(spec, Mapping):
            raise RunnerError(f"invalid sample input identity: {name}")
        verify_identity(spec, f"sample partition input {name}")
    verify_sidecar(path)
    return receipt


def create_sample_partition_receipt(
    sample_root: Path,
    *,
    sample: str,
    chroms: Sequence[str],
    chrom_receipts: Mapping[str, Mapping[str, Any]],
    manifest_path: Path,
    args: argparse.Namespace,
    full_scope: bool,
) -> tuple[Path, dict[str, Any]]:
    receipt_name = "run_receipt.json" if full_scope else "partial_partition_receipt.json"
    receipt_path = sample_root / receipt_name
    if receipt_path.exists():
        if not args.resume:
            raise RunnerError(f"sample receipt already exists: {receipt_path}")
        return receipt_path, validate_sample_partition_receipt(
            receipt_path, sample=sample, expected_chroms=chroms
        )
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.sample_partition_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "task_type": (
            "B_comprehensive_validation" if full_scope else "A_exploratory_pilot"
        ),
        "claim_status": "TECHNICAL_PARTITION_ONLY",
        "validation_evidence_eligible": full_scope,
        "all_pass": True,
        "sample": sample,
        "scope": {
            "chromosomes": list(chroms),
            "autosomes_complete": full_scope,
        },
        "parameters": {
            "threshold": 3,
            "max_active_bits": 12,
            "partition_membership_source": "strict_endpoint_regions",
        },
        "checks": {
            "all_requested_chromosomes_present": len(chrom_receipts) == len(chroms),
            "all_chromosome_receipts_pass": all(
                value.get("all_pass") is True for value in chrom_receipts.values()
            ),
            "all_python_cpp_comparisons_pass": all(
                value.get("checks", {}).get("python_cpp_mismatch_zero") is True
                for value in chrom_receipts.values()
            ),
            "full_receipt_requires_chr1_22": (not full_scope)
            or tuple(chroms) == AUTOSOMES,
        },
        "counts": {
            "chromosomes_requested": len(chroms),
            "chromosomes_pass": sum(
                value.get("all_pass") is True for value in chrom_receipts.values()
            ),
            "child_checks": sum(
                int(value.get("check_count", 0)) for value in chrom_receipts.values()
            ),
        },
        "chromosome_receipts": {
            chrom: identity(sample_root / "chromosomes" / chrom / "stage_receipt.json")
            for chrom in chroms
        },
        "inputs": {
            "manifest": identity(manifest_path),
            "partitioner": identity(args.partitioner),
            "partition_binary": identity(args.partition_binary),
            "comparator": identity(args.comparator),
        },
        "claim_ceiling": (
            "This receipt proves strict-membership exact-PS k<=12 partitioning and "
            "Python/C++ partition parity. It does not prove topology completeness, "
            "a unique true tree, or cellular clone/subclone identity."
        ),
    }
    require_all_true(receipt["checks"], "sample partition construction")
    write_json_exclusive(receipt_path, receipt)
    write_sidecar_exclusive(receipt_path)
    return receipt_path, receipt


def validate_mlhp_output(output_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    receipt_path = output_path.with_suffix(output_path.suffix + ".receipt.json")
    receipt = read_json(receipt_path)
    if (
        receipt.get("schema_name") != MLHP_RECEIPT_SCHEMA
        or receipt.get("all_pass") is not True
    ):
        raise RunnerError(f"MLHP adapter receipt mismatch: {receipt_path}")
    output_spec = receipt.get("output")
    if not isinstance(output_spec, Mapping) or not identity_matches(output_spec, output_path):
        raise RunnerError("MLHP receipt output identity mismatch")
    document = read_json(output_path)
    if document.get("sample") is None or not isinstance(document.get("groups"), list):
        raise RunnerError("MLHP output contract mismatch")
    return receipt, document


def validate_topology_receipt(
    receipt_path: Path,
    *,
    input_path: Path,
    output_path: Path,
    max_family_size: int,
    max_search_nodes: int,
) -> dict[str, Any]:
    receipt = read_json(receipt_path)
    if receipt.get("all_pass") is not True:
        raise RunnerError(f"topology technical receipt is not PASS: {receipt_path}")
    input_spec = receipt.get("input")
    output_spec = receipt.get("output")
    if not isinstance(input_spec, Mapping) or not identity_matches(input_spec, input_path):
        raise RunnerError("topology receipt input identity mismatch")
    if not isinstance(output_spec, Mapping) or not identity_matches(output_spec, output_path):
        raise RunnerError("topology receipt output identity mismatch")
    completeness = receipt.get("all_mutation_bearing_families_complete")
    if not isinstance(completeness, bool):
        raise RunnerError("topology receipt lacks explicit family completeness")
    objective_certified = receipt.get("all_units_objective_certified")
    if not isinstance(objective_certified, bool):
        raise RunnerError("topology receipt lacks explicit objective certification")
    parameters = receipt.get("parameters")
    if not isinstance(parameters, Mapping) or (
        parameters.get("max_family_size"),
        parameters.get("max_search_nodes"),
    ) != (max_family_size, max_search_nodes):
        raise RunnerError("topology receipt guard parameters mismatch")
    census = receipt.get("status_census")
    required_census = {
        "unit_status",
        "objective_status",
        "family_status",
        "read_af_status",
    }
    if (
        not isinstance(census, Mapping)
        or set(census) != required_census
        or any(
            not isinstance(subcensus, Mapping)
            or any(
                isinstance(value, bool) or not isinstance(value, int) or value < 0
                for value in subcensus.values()
            )
            for subcensus in census.values()
        )
    ):
        raise RunnerError("topology receipt status_census is invalid")
    counts = receipt.get("counts")
    if (
        not isinstance(counts, Mapping)
        or isinstance(counts.get("total_units"), bool)
        or not isinstance(counts.get("total_units"), int)
        or counts["total_units"] < 0
    ):
        raise RunnerError("topology receipt total-unit count is invalid")
    rows = 0
    with output_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                raise RunnerError(f"blank JSONL row at line {line_number}")
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RunnerError(
                    f"invalid topology JSONL at line {line_number}: {exc}"
                ) from exc
            if not isinstance(value, dict):
                raise RunnerError(f"topology JSONL row is not an object: {line_number}")
            rows += 1
    if rows == 0:
        raise RunnerError("topology JSONL is empty")
    if rows != counts["total_units"]:
        raise RunnerError("topology JSONL row count differs from receipt")
    receipt["_verified_jsonl_rows"] = rows
    return receipt


def run_full_sample_downstream(
    *,
    sample: str,
    sample_root: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    mlhp_path = sample_root / f"{sample}.exact_ps_mlhp.json"
    mlhp_receipt_path = mlhp_path.with_suffix(mlhp_path.suffix + ".receipt.json")
    adapter_stage: dict[str, Any]
    if mlhp_path.exists() or mlhp_receipt_path.exists():
        if not args.resume or not mlhp_path.is_file() or not mlhp_receipt_path.is_file():
            raise RunnerError(f"non-resumable MLHP output state: {mlhp_path}")
        mlhp_receipt, mlhp = validate_mlhp_output(mlhp_path)
        adapter_stage = {
            "stage": "exact_ps_partition_to_mlhp",
            "resumed": True,
            "wall_seconds": 0.0,
        }
    else:
        command = [
            str(args.python),
            str(args.adapter),
            "--partition-root",
            str(sample_root),
            "--output",
            str(mlhp_path),
            "--sample",
            sample,
            "--chroms",
            ",".join(AUTOSOMES),
            "--min-read",
            "3",
            "--require-strict-endpoint-receipt",
        ]
        adapter_stage = run_command(
            stage="exact_ps_partition_to_mlhp",
            command=command,
            logs_dir=sample_root / "logs",
        )
        mlhp_receipt, mlhp = validate_mlhp_output(mlhp_path)
    if mlhp.get("sample") != sample:
        raise RunnerError("MLHP sample mismatch")

    topology_path = sample_root / f"{sample}.topology.jsonl"
    topology_receipt_path = sample_root / f"{sample}.topology.receipt.json"
    topology_stage: dict[str, Any]
    if topology_path.exists() or topology_receipt_path.exists():
        if (
            not args.resume
            or not topology_path.is_file()
            or not topology_receipt_path.is_file()
        ):
            raise RunnerError(f"non-resumable topology output state: {topology_path}")
        topology_receipt = validate_topology_receipt(
            topology_receipt_path,
            input_path=mlhp_path,
            output_path=topology_path,
            max_family_size=args.max_family_size,
            max_search_nodes=args.max_search_nodes,
        )
        topology_stage = {
            "stage": "cpp_topology",
            "resumed": True,
            "wall_seconds": 0.0,
        }
    else:
        command = [
            str(args.topology_binary),
            "--input",
            str(mlhp_path),
            "--output",
            str(topology_path),
            "--receipt",
            str(topology_receipt_path),
            "--max-family-size",
            str(args.max_family_size),
            "--max-search-nodes",
            str(args.max_search_nodes),
        ]
        topology_stage = run_command(
            stage="cpp_topology",
            command=command,
            logs_dir=sample_root / "logs",
        )
        topology_receipt = validate_topology_receipt(
            topology_receipt_path,
            input_path=mlhp_path,
            output_path=topology_path,
            max_family_size=args.max_family_size,
            max_search_nodes=args.max_search_nodes,
        )

    pipeline_path = sample_root / "pipeline_receipt.json"
    if pipeline_path.exists():
        if not args.resume:
            raise RunnerError(f"pipeline receipt already exists: {pipeline_path}")
        pipeline = read_json(pipeline_path)
        if (
            pipeline.get("schema_name") != f"{SCHEMA_NAME}.sample_pipeline_receipt"
            or pipeline.get("all_pass") is not True
            or pipeline.get("sample") != sample
        ):
            raise RunnerError(f"pipeline receipt is not resumable: {pipeline_path}")
        expected_parameters = {
            "max_active_bits": 12,
            "max_family_size": args.max_family_size,
            "max_search_nodes": args.max_search_nodes,
        }
        if pipeline.get("parameters") != expected_parameters:
            raise RunnerError(
                "pipeline receipt guard parameters differ from this invocation"
            )
        for spec in pipeline.get("bound_artifacts", []):
            if not isinstance(spec, Mapping):
                raise RunnerError("invalid pipeline bound artifact")
            verify_identity(spec, "pipeline bound artifact")
        verify_sidecar(pipeline_path)
        return pipeline

    complete = topology_receipt["all_mutation_bearing_families_complete"]
    objective_certified = topology_receipt["all_units_objective_certified"]
    pipeline = {
        "schema_name": f"{SCHEMA_NAME}.sample_pipeline_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "technical_all_pass": True,
        "all_mutation_bearing_families_complete": complete,
        "all_units_objective_certified": objective_certified,
        "sample": sample,
        "parameters": {
            "max_active_bits": 12,
            "max_family_size": args.max_family_size,
            "max_search_nodes": args.max_search_nodes,
        },
        "stages": [adapter_stage, topology_stage],
        "status_census": topology_receipt["status_census"],
        "counts": {
            "mlhp_groups": len(mlhp["groups"]),
            "topology_jsonl_rows": topology_receipt["_verified_jsonl_rows"],
        },
        "bound_artifacts": bound_artifacts(
            [
                sample_root / "run_receipt.json",
                mlhp_path,
                mlhp_receipt_path,
                topology_path,
                topology_receipt_path,
            ]
        ),
        "claim_ceiling": (
            "Technical PASS means the pipeline and identity gates completed. "
            + (
                "All mutation-bearing families completed and all units carry objective "
                "certificates under the configured guards; "
                if complete and objective_certified
                else (
                    "At least one family ABSTAINED, hit a configured guard, or "
                    "lacks an objective certificate; "
                )
            )
            + "the output remains a minimum recurrence-allowed model ranked by "
            "read-AF, not proof of a unique biological clone tree."
        ),
    }
    write_json_exclusive(pipeline_path, pipeline)
    write_sidecar_exclusive(pipeline_path)
    return pipeline


def prepare_output_root(path: Path, resume: bool) -> Path:
    root = path.expanduser().absolute()
    if root.is_symlink():
        raise RunnerError(f"output root must not be a symlink: {root}")
    if root.exists():
        if not root.is_dir():
            raise RunnerError(f"output root is not a directory: {root}")
        if not resume and next(root.iterdir(), None) is not None:
            raise RunnerError(f"output root must be new or empty: {root}")
    else:
        root.mkdir(parents=True, exist_ok=False)
    return root.resolve()


def execute(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    args.manifest = resolve_regular_file(args.manifest, "input manifest")
    _, sources = load_manifest(args.manifest)
    args.partitioner = resolve_regular_file(args.partitioner, "Python partitioner")
    args.comparator = resolve_regular_file(args.comparator, "partition comparator")
    args.adapter = resolve_regular_file(args.adapter, "MLHP adapter")
    args.python = resolve_executable(args.python, "Python")
    args.partition_binary = resolve_executable(
        args.partition_binary, "C++ partition binary"
    )
    full_scope = tuple(args.chroms) == AUTOSOMES
    if full_scope:
        if args.topology_binary is None:
            raise RunnerError("--topology-binary is required for chr1-22 runs")
        args.topology_binary = resolve_executable(
            args.topology_binary, "C++ topology binary"
        )
    args.output_root = prepare_output_root(args.output_root, args.resume)

    sample_results: dict[str, dict[str, Any]] = {}
    for sample in args.samples:
        sample_root = args.output_root / "samples" / sample
        if sample_root.exists():
            if not args.resume or not sample_root.is_dir():
                raise RunnerError(f"sample output already exists: {sample_root}")
        else:
            sample_root.mkdir(parents=True, exist_ok=False)
        chrom_receipts: dict[str, dict[str, Any]] = {}
        if args.workers == 1 or len(args.chroms) == 1:
            for chrom in args.chroms:
                chrom_receipts[chrom] = process_chromosome(
                    sample=sample,
                    chrom=chrom,
                    source=sources[sample],
                    sample_root=sample_root,
                    args=args,
                )
        else:
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                future_to_chrom = {
                    executor.submit(
                        process_chromosome,
                        sample=sample,
                        chrom=chrom,
                        source=sources[sample],
                        sample_root=sample_root,
                        args=args,
                    ): chrom
                    for chrom in args.chroms
                }
                for future in as_completed(future_to_chrom):
                    chrom = future_to_chrom[future]
                    chrom_receipts[chrom] = future.result()
            chrom_receipts = {
                chrom: chrom_receipts[chrom] for chrom in args.chroms
            }
        _, partition_receipt = create_sample_partition_receipt(
            sample_root,
            sample=sample,
            chroms=args.chroms,
            chrom_receipts=chrom_receipts,
            manifest_path=args.manifest,
            args=args,
            full_scope=full_scope,
        )
        if full_scope:
            pipeline = run_full_sample_downstream(
                sample=sample,
                sample_root=sample_root,
                args=args,
            )
            sample_results[sample] = {
                "partition": partition_receipt,
                "pipeline": pipeline,
            }
        else:
            sample_results[sample] = {
                "partition": partition_receipt,
                "pipeline": None,
            }

    receipt_path = (
        args.output_root / "cohort_receipt.json"
        if full_scope
        else args.output_root / "partial_cohort_receipt.json"
    )
    if receipt_path.exists():
        if not args.resume:
            raise RunnerError(f"cohort receipt already exists: {receipt_path}")
        receipt = read_json(receipt_path)
        scope = receipt.get("scope")
        parameters = receipt.get("parameters")
        if (
            receipt.get("schema_name") != f"{SCHEMA_NAME}.cohort_receipt"
            or receipt.get("schema_version") != SCHEMA_VERSION
            or receipt.get("all_pass") is not True
            or not isinstance(scope, Mapping)
            or scope.get("samples") != list(args.samples)
            or scope.get("chromosomes") != list(args.chroms)
            or not isinstance(parameters, Mapping)
            or parameters.get("max_family_size") != args.max_family_size
            or parameters.get("max_search_nodes") != args.max_search_nodes
        ):
            raise RunnerError("existing cohort receipt is not PASS")
        sample_specs = receipt.get("sample_receipts")
        if not isinstance(sample_specs, Mapping):
            raise RunnerError("cohort receipt sample identities are missing")
        for sample, stage_specs in sample_specs.items():
            if not isinstance(stage_specs, Mapping):
                raise RunnerError(f"invalid cohort sample identity block: {sample}")
            for stage, spec in stage_specs.items():
                if spec is None:
                    continue
                if not isinstance(spec, Mapping):
                    raise RunnerError(f"invalid cohort identity: {sample}/{stage}")
                verify_identity(spec, f"cohort identity {sample}/{stage}")
        verify_sidecar(receipt_path)
        return receipt

    technical_all_pass = all(
        result["partition"]["all_pass"] is True
        and (
            result["pipeline"] is None
            or result["pipeline"]["technical_all_pass"] is True
        )
        for result in sample_results.values()
    )
    family_complete = full_scope and all(
        result["pipeline"]["all_mutation_bearing_families_complete"] is True
        for result in sample_results.values()
    )
    objective_complete = full_scope and all(
        result["pipeline"]["all_units_objective_certified"] is True
        for result in sample_results.values()
    )
    all7 = tuple(args.samples) == CANONICAL_DATASETS
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.cohort_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "task_type": (
            "B_comprehensive_validation"
            if full_scope and all7
            else "A_exploratory_pilot"
        ),
        "all_pass": technical_all_pass,
        "technical_all_pass": technical_all_pass,
        "all_mutation_bearing_families_complete": family_complete,
        "all_units_objective_certified": objective_complete,
        "validation_evidence_eligible": (
            full_scope and all7 and family_complete and objective_complete
        ),
        "scope": {
            "samples": list(args.samples),
            "chromosomes": list(args.chroms),
            "canonical_all7": all7,
            "autosomes_complete": full_scope,
        },
        "parameters": {
            "threshold": 3,
            "max_active_bits": 12,
            "max_family_size": args.max_family_size,
            "max_search_nodes": args.max_search_nodes,
            "chromosome_workers": args.workers,
            "resume": args.resume,
        },
        "wall_seconds": time.perf_counter() - started,
        "sample_receipts": {
            sample: {
                "partition": identity(
                    args.output_root / "samples" / sample
                    / (
                        "run_receipt.json"
                        if full_scope
                        else "partial_partition_receipt.json"
                    )
                ),
                "pipeline": (
                    identity(
                        args.output_root
                        / "samples"
                        / sample
                        / "pipeline_receipt.json"
                    )
                    if full_scope
                    else None
                ),
            }
            for sample in args.samples
        },
        "claim_ceiling": (
            "Technical PASS is reported separately from family completeness. "
            "A cohort with any guard-triggered ABSTAIN is not topology-complete "
            "and cannot support a claim that every family has a resolved tree."
        ),
    }
    write_json_exclusive(receipt_path, receipt)
    write_sidecar_exclusive(receipt_path)
    return receipt


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--samples",
        default=",".join(CANONICAL_DATASETS),
        help="comma-separated subset in canonical order",
    )
    parser.add_argument(
        "--chroms",
        default=",".join(AUTOSOMES),
        help="comma-separated subset of chr1-chr22",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--workers",
        type=int,
        choices=tuple(range(1, 9)),
        default=4,
        help="independent chromosome workers per sample (1-8)",
    )
    parser.add_argument("--partitioner", type=Path, default=DEFAULT_PARTITIONER)
    parser.add_argument("--comparator", type=Path, default=DEFAULT_COMPARATOR)
    parser.add_argument("--adapter", type=Path, default=DEFAULT_ADAPTER)
    parser.add_argument(
        "--partition-binary", type=Path, default=DEFAULT_PARTITION_BINARY
    )
    parser.add_argument("--topology-binary", type=Path)
    parser.add_argument("--max-family-size", type=int, default=100_000)
    parser.add_argument("--max-search-nodes", type=int, default=100_000_000)
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args(argv)
    try:
        args.samples = parse_csv(args.samples, CANONICAL_DATASETS, "samples")
        args.chroms = parse_csv(args.chroms, AUTOSOMES, "chroms")
    except argparse.ArgumentTypeError as exc:
        parser.error(str(exc))
    if args.max_family_size < 0 or args.max_search_nodes < 0:
        parser.error("max-family-size and max-search-nodes must be nonnegative")
    if tuple(args.samples) != tuple(
        sample for sample in CANONICAL_DATASETS if sample in args.samples
    ):
        parser.error("samples must follow canonical dataset order")
    if tuple(args.chroms) != tuple(chrom for chrom in AUTOSOMES if chrom in args.chroms):
        parser.error("chroms must follow chr1-chr22 order")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    try:
        receipt = execute(parse_args(argv))
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "technical_all_pass": receipt["technical_all_pass"],
                "all_mutation_bearing_families_complete": receipt[
                    "all_mutation_bearing_families_complete"
                ],
                "scope": receipt["scope"],
                "wall_seconds": receipt["wall_seconds"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if receipt["technical_all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
