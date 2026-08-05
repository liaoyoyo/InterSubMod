#!/usr/bin/env python3
"""Production launcher for strict exact-PS/HP read-linked mutation regions.

The v4 scientific unit is deliberately narrower than a phase set:

* container: dataset x chromosome x exact non-missing PS x HP1/HP2;
* node: sSNV with a fixed R/A observation in the container;
* edge: at least three distinct canonical molecules fixed R/A at both ends;
* region: a connected component of those endpoint edges;
* tree input: only multisite regions (k > 1).  Unlinked singletons are retained
  in the regional funnel as abstentions and can never enter partition/topology.

The launcher defaults to all seven canonical datasets and chr1-chr22.  A
subset is possible only with ``--allow-partial-scope`` and is labelled
ineligible for production validation.  Every output root must be new; child
artifacts, receipts and SHA-256 sidecars are never overwritten.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
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
from typing import Any, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
CANONICAL_DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
STAGES = {"regions": 1, "partition": 2, "topology": 3}
SCHEMA_NAME = "intersubmod.layered_v4_strict_production_run"
SCHEMA_VERSION = "1.0.0"
PRIMARY_THRESHOLD = 3
MAX_BLOCK_SIZE = 12

DEFAULT_EXTRACTOR = (
    REPO_ROOT
    / "research/20260718_k_gt8_read_supported_segmentation/scripts/"
    "extract_lossless_read_linkage_collapsing.py"
)
DEFAULT_REGION_BUILDER = REPO_ROOT / "scripts/build_strict_ps_hp_regions.py"
DEFAULT_PARTITIONER = (
    REPO_ROOT
    / "research/20260722_exact_ps_k12_hcc1395_pilot/scripts/"
    "exact_ps_k12_partition.py"
)
DEFAULT_PARTITION_CPP = (
    REPO_ROOT
    / "research/20260722_exact_ps_k12_hcc1395_pilot/cpp/"
    "exact_ps_k12_partition.cpp"
)
DEFAULT_PARTITION_COMPARATOR = (
    REPO_ROOT
    / "research/20260722_exact_ps_k12_hcc1395_pilot/scripts/"
    "compare_python_cpp_partitions.py"
)
DEFAULT_STRICT_CPP = (
    REPO_ROOT
    / "research/20260723_production_exact_ps_strict_read_linkage/cpp/"
    "strict_endpoint_graph_verify.cpp"
)
DEFAULT_STRICT_CPP_HEADER = DEFAULT_STRICT_CPP.with_name("strict_endpoint_graph.hpp")
DEFAULT_ADAPTER = (
    REPO_ROOT
    / "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/"
    "exact_ps_partition_to_mlhp.py"
)
DEFAULT_LAYERED = (
    REPO_ROOT
    / "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/"
    "layered_tree_reconstruction.py"
)
DEFAULT_REGION_VIEW = (
    REPO_ROOT
    / "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/"
    "build_region_view.py"
)


class RunnerError(RuntimeError):
    """A fail-closed production contract violation."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    path = path.resolve(strict=True)
    if not path.is_file():
        raise RunnerError(f"not a regular file: {path}")
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": stat.st_size,
        "sha256": sha256_path(path),
    }


def read_json_strict(path: Path) -> Any:
    def reject_duplicate(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise RunnerError(f"{path}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=reject_duplicate)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RunnerError(f"cannot read strict JSON {path}: {exc}") from exc


def write_json_exclusive(path: Path, document: Any) -> None:
    payload = json.dumps(
        document, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False
    ) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def write_sha256_sidecar(path: Path) -> Path:
    sidecar = path.with_name(path.name + ".sha256")
    with sidecar.open("x", encoding="ascii") as handle:
        handle.write(f"{sha256_path(path)}  {path.name}\n")
        handle.flush()
        os.fsync(handle.fileno())
    return sidecar


def write_failure_marker(run_root: Path, exc: BaseException) -> None:
    failure_path = run_root / "_FAILED"
    if failure_path.exists():
        return
    write_json_exclusive(
        failure_path,
        {
            "schema_name": f"{SCHEMA_NAME}.failure",
            "schema_version": SCHEMA_VERSION,
            "all_pass": False,
            "failed_at_utc": utc_now(),
            "error_type": type(exc).__name__,
            "error": str(exc),
        },
    )


def validate_sha256_sidecar(path: Path) -> None:
    sidecar = path.with_name(path.name + ".sha256")
    try:
        tokens = sidecar.read_text(encoding="ascii").strip().split()
    except OSError as exc:
        raise RunnerError(f"missing SHA-256 sidecar for {path}") from exc
    if len(tokens) != 2 or tokens[1] != path.name or tokens[0] != sha256_path(path):
        raise RunnerError(f"invalid SHA-256 sidecar: {sidecar}")


def parse_unique_csv(value: str, allowed: Sequence[str], label: str) -> tuple[str, ...]:
    values = tuple(token for token in value.split(",") if token)
    if not values or len(values) != len(set(values)):
        raise argparse.ArgumentTypeError(f"{label} must be a non-empty unique CSV")
    unknown = sorted(set(values) - set(allowed))
    if unknown:
        raise argparse.ArgumentTypeError(f"unknown {label}: {','.join(unknown)}")
    return values


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--datasets",
        type=lambda value: parse_unique_csv(value, CANONICAL_DATASETS, "datasets"),
        default=CANONICAL_DATASETS,
    )
    parser.add_argument(
        "--chromosomes",
        type=lambda value: parse_unique_csv(value, AUTOSOMES, "chromosomes"),
        default=AUTOSOMES,
    )
    parser.add_argument(
        "--allow-partial-scope",
        action="store_true",
        help="Allow an explicitly PARTIAL, validation-ineligible dataset/chromosome subset.",
    )
    parser.add_argument(
        "--stage-through",
        choices=tuple(STAGES),
        default="partition",
        help=(
            "Default partition is the currently runtime-validated production boundary; "
            "topology additionally invokes the strict-receipt adapter and tree renderers."
        ),
    )
    parser.add_argument(
        "--extraction-cache-pattern",
        help=(
            "Read-only cache directory pattern containing {chrom} and optionally {sample}; "
            "the validated directory is symlinked into the new run."
        ),
    )
    parser.add_argument(
        "--python",
        default=None,
        help="Python >=3.10 executable; auto-detected when omitted.",
    )
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--gzip", default="gzip")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--samtools-threads", type=int, default=1)
    parser.add_argument("--extractor", type=Path, default=DEFAULT_EXTRACTOR)
    parser.add_argument("--region-builder", type=Path, default=DEFAULT_REGION_BUILDER)
    parser.add_argument("--partitioner", type=Path, default=DEFAULT_PARTITIONER)
    parser.add_argument("--partition-cpp", type=Path, default=DEFAULT_PARTITION_CPP)
    parser.add_argument("--partition-comparator", type=Path, default=DEFAULT_PARTITION_COMPARATOR)
    parser.add_argument("--strict-cpp", type=Path, default=DEFAULT_STRICT_CPP)
    parser.add_argument("--strict-cpp-header", type=Path, default=DEFAULT_STRICT_CPP_HEADER)
    parser.add_argument("--adapter", type=Path, default=DEFAULT_ADAPTER)
    parser.add_argument("--layered-reconstruction", type=Path, default=DEFAULT_LAYERED)
    parser.add_argument("--region-view", type=Path, default=DEFAULT_REGION_VIEW)
    args = parser.parse_args(argv)
    args.invocation = [str(Path(__file__).resolve()), *(argv or sys.argv[1:])]
    return args


def validate_scope(args: argparse.Namespace) -> bool:
    full = tuple(args.datasets) == CANONICAL_DATASETS and tuple(args.chromosomes) == AUTOSOMES
    if not full and not args.allow_partial_scope:
        raise RunnerError(
            "dataset/chromosome subset requires --allow-partial-scope; production defaults "
            "are all seven datasets and chr1-chr22"
        )
    if (
        args.extraction_cache_pattern is not None
        and "{chrom}" not in args.extraction_cache_pattern
    ):
        raise RunnerError("--extraction-cache-pattern must contain {chrom}; {sample} is optional")
    if args.samtools_threads < 1:
        raise RunnerError("--samtools-threads must be positive")
    return full


def load_manifest(
    path: Path, *, allow_partial: bool = False
) -> tuple[dict[str, Mapping[str, Any]], dict[str, Any]]:
    document = read_json_strict(path)
    accepted_schemas = {"intersubmod.layered_input_manifest"}
    if allow_partial:
        accepted_schemas.add("intersubmod.big7_hcc1395_partial_pilot_manifest")
    if not isinstance(document, dict) or document.get("schema_name") not in accepted_schemas:
        raise RunnerError("manifest is not an InterSubMod layered input manifest")
    samples = document.get("samples")
    if not isinstance(samples, list):
        raise RunnerError("manifest samples must be a list")
    entries: dict[str, Mapping[str, Any]] = {}
    for item in samples:
        if not isinstance(item, dict) or not isinstance(item.get("sample"), str):
            raise RunnerError("manifest contains a malformed sample entry")
        sample = item["sample"]
        if sample in entries:
            raise RunnerError(f"manifest contains duplicate sample {sample}")
        entries[sample] = item
    if not allow_partial and set(entries) != set(CANONICAL_DATASETS):
        raise RunnerError(
            f"canonical manifest dataset set mismatch: {sorted(entries)} != {sorted(CANONICAL_DATASETS)}"
        )
    if not entries or not set(entries).issubset(CANONICAL_DATASETS):
        raise RunnerError(f"manifest contains an unsupported dataset set: {sorted(entries)}")
    return entries, document


def require_tool(path: Path, label: str) -> Path:
    path = path.resolve(strict=True)
    if not path.is_file():
        raise RunnerError(f"{label} is not a regular file: {path}")
    return path


def require_executable(value: str, label: str) -> str:
    candidate = shutil.which(value)
    if candidate is None:
        raise RunnerError(f"cannot resolve {label} executable: {value}")
    return str(Path(candidate).resolve())


def require_python(value: str | None) -> str:
    candidates: list[str] = []
    if value is not None:
        candidates.append(value)
    else:
        candidates.extend([sys.executable, "python3", "/usr/bin/python3"])
    seen: set[str] = set()
    failures: list[str] = []
    for candidate in candidates:
        resolved = shutil.which(candidate) if not Path(candidate).is_absolute() else candidate
        if not resolved:
            failures.append(f"{candidate}:not-found")
            continue
        resolved = str(Path(resolved).resolve())
        if resolved in seen:
            continue
        seen.add(resolved)
        completed = subprocess.run(
            [resolved, "-c", "import sys; raise SystemExit(0 if sys.version_info >= (3,10) else 1)"],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if completed.returncode == 0:
            return resolved
        failures.append(f"{resolved}:requires-Python>=3.10")
    raise RunnerError("no compatible Python interpreter: " + "; ".join(failures))


def prepare_output_root(path: Path) -> Path:
    path = path.expanduser().absolute()
    if path.exists() or path.is_symlink():
        raise RunnerError(f"output root must not already exist: {path}")
    path.mkdir(parents=True, exist_ok=False)
    return path.resolve()


def run_command(
    *,
    stage: str,
    command: Sequence[str],
    logs_dir: Path,
    env: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / f"{stage}.stdout.log"
    stderr_path = logs_dir / f"{stage}.stderr.log"
    started = time.perf_counter()
    started_at = utc_now()
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        completed = subprocess.run(
            list(command),
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            env=dict(env) if env is not None else None,
            check=False,
        )
    record = {
        "stage": stage,
        "command": list(command),
        "command_shell": shlex.join(command),
        "started_at_utc": started_at,
        "ended_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - started,
        "exit_code": completed.returncode,
        "all_pass": completed.returncode == 0,
        "logs": {"stdout": file_identity(stdout_path), "stderr": file_identity(stderr_path)},
    }
    if completed.returncode != 0:
        tail = stderr_path.read_text(encoding="utf-8", errors="replace")[-4000:]
        raise RunnerError(f"stage {stage} failed with exit {completed.returncode}: {tail}")
    return record


def run_streamed_strict_cpp(
    *,
    gzip_executable: str,
    binary: Path,
    molecule_calls: Path,
    output_dir: Path,
    logs_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=False)
    stdout_path = logs_dir / "strict_cpp.stdout.log"
    stderr_path = logs_dir / "strict_cpp.stderr.log"
    producer_stderr_path = logs_dir / "strict_cpp_gzip.stderr.log"
    command_gzip = [gzip_executable, "-dc", str(molecule_calls)]
    command_cpp = [
        str(binary),
        "--input",
        "-",
        "--threshold",
        str(PRIMARY_THRESHOLD),
        "--edges-output",
        str(output_dir / "edges.tsv"),
        "--components-output",
        str(output_dir / "components.tsv"),
        "--receipt-output",
        str(output_dir / "receipt.json"),
    ]
    started = time.perf_counter()
    started_at = utc_now()
    logs_dir.mkdir(parents=True, exist_ok=True)
    with (
        stdout_path.open("xb") as stdout,
        stderr_path.open("xb") as stderr,
        producer_stderr_path.open("xb") as producer_stderr,
    ):
        producer = subprocess.Popen(
            command_gzip,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=producer_stderr,
        )
        assert producer.stdout is not None
        consumer = subprocess.Popen(
            command_cpp,
            stdin=producer.stdout,
            stdout=stdout,
            stderr=stderr,
        )
        producer.stdout.close()
        consumer_code = consumer.wait()
        producer_code = producer.wait()
    if producer_code != 0 or consumer_code != 0:
        raise RunnerError(
            f"strict C++ stream failed: gzip={producer_code}, verifier={consumer_code}"
        )
    return {
        "stage": "strict_cpp",
        "commands": [command_gzip, command_cpp],
        "command_shell": f"{shlex.join(command_gzip)} | {shlex.join(command_cpp)}",
        "started_at_utc": started_at,
        "ended_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - started,
        "exit_code": 0,
        "all_pass": True,
        "logs": {
            "stdout": file_identity(stdout_path),
            "stderr": file_identity(stderr_path),
            "gzip_stderr": file_identity(producer_stderr_path),
        },
    }


def expected_extraction_inputs(entry: Mapping[str, Any]) -> dict[str, Any]:
    try:
        return {
            "bam": entry["alignment_payload"],
            "tree_vcf": entry["somatic"]["tree_vcf"],
            "read_tag_sidecar": entry["read_tags"]["sidecar"],
            "read_tag_sidecar_index": entry["read_tags"]["index"],
        }
    except (KeyError, TypeError) as exc:
        raise RunnerError("manifest sample lacks extractor input identities") from exc


def validate_output_identity(declared: Mapping[str, Any], path: Path, label: str) -> None:
    observed = file_identity(path)
    if declared.get("sha256") != observed["sha256"] or declared.get("size_bytes") != observed["size_bytes"]:
        raise RunnerError(f"{label} identity drift: {path}")


def validate_extraction(
    directory: Path,
    *,
    sample: str,
    chrom: str,
    entry: Mapping[str, Any],
) -> dict[str, Any]:
    receipt_path = directory / "receipt.json"
    validate_sha256_sidecar(receipt_path)
    receipt = read_json_strict(receipt_path)
    if (
        receipt.get("schema_name") != "intersubmod.lossless_read_linkage_chromosome_receipt"
        or receipt.get("schema_version") != "1.3.0"
        or receipt.get("all_pass") is not True
        or receipt.get("scope", {}).get("dataset") != sample
        or receipt.get("scope", {}).get("chrom") != chrom
    ):
        raise RunnerError(f"invalid extraction receipt scope/schema: {receipt_path}")
    parameters = receipt.get("parameters") or {}
    if (
        parameters.get("mapq_min") != 20
        or parameters.get("baseq_min") != 20
        or parameters.get("bridge_thresholds") != [1, 2, 3, 5]
    ):
        raise RunnerError(f"extraction parameter drift: {receipt_path}")
    if receipt.get("inputs") != expected_extraction_inputs(entry):
        raise RunnerError(f"extraction input identity does not match current manifest: {receipt_path}")
    checks = receipt.get("checks") or {}
    if not checks or any(value is not True for value in checks.values()):
        raise RunnerError(f"extraction receipt contains a failed check: {receipt_path}")
    prefix = f"{sample}.{chrom}"
    needed = {
        f"{prefix}.site_catalog.tsv.gz": directory / f"{prefix}.site_catalog.tsv.gz",
        f"{prefix}.molecule_sparse_calls.tsv.gz": directory / f"{prefix}.molecule_sparse_calls.tsv.gz",
    }
    declared_outputs = receipt.get("outputs") or {}
    for name, path in needed.items():
        if name not in declared_outputs:
            raise RunnerError(f"extraction receipt does not bind {name}")
        validate_output_identity(declared_outputs[name], path, f"extraction output {name}")
    return receipt


def materialize_extraction(
    *,
    args: argparse.Namespace,
    sample: str,
    chrom: str,
    entry: Mapping[str, Any],
    chrom_root: Path,
    logs_dir: Path,
) -> tuple[Path, dict[str, Any]]:
    destination = chrom_root / "extraction"
    if args.extraction_cache_pattern:
        cache = Path(
            args.extraction_cache_pattern.format(sample=sample, chrom=chrom)
        ).expanduser().resolve(strict=True)
        if not cache.is_dir() or cache.is_symlink():
            raise RunnerError(f"cache must resolve to a real directory: {cache}")
        receipt = validate_extraction(cache, sample=sample, chrom=chrom, entry=entry)
        os.symlink(cache, destination, target_is_directory=True)
        return destination, {
            "stage": "extraction_cache",
            "all_pass": True,
            "read_only": True,
            "cache": str(cache),
            "symlink": str(destination),
            "receipt": file_identity(cache / "receipt.json"),
        }

    command = [
        args.python,
        str(args.extractor),
        "--manifest",
        str(args.manifest),
        "--sample",
        sample,
        "--chrom",
        chrom,
        "--output-dir",
        str(destination),
        "--mapq-min",
        "20",
        "--baseq-min",
        "20",
        "--bridge-thresholds",
        "1,2,3,5",
        "--samtools-threads",
        str(args.samtools_threads),
        "--samtools",
        args.samtools,
    ]
    stage = run_command(stage="extract", command=command, logs_dir=logs_dir)
    validate_extraction(destination, sample=sample, chrom=chrom, entry=entry)
    stage["receipt"] = file_identity(destination / "receipt.json")
    return destination, stage


def read_tsv(path: Path, *, compressed: bool) -> list[dict[str, str]]:
    opener = gzip.open if compressed else open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise RunnerError(f"TSV has no header: {path}")
            return list(reader)
    except (OSError, EOFError, UnicodeError, csv.Error) as exc:
        raise RunnerError(f"cannot read TSV {path}: {exc}") from exc


def validate_strict_membership(
    membership_path: Path,
    *,
    sample: str,
    chrom: str,
    threshold: int = PRIMARY_THRESHOLD,
) -> dict[str, int]:
    """Reject legacy membership and prove that no singleton enters PRIMARY."""

    required = {
        "dataset",
        "chrom",
        "linkage_basis",
        "phase_set",
        "phase_set_status",
        "inference_role",
        "component_class",
        "tree_eligible",
        "threshold",
        "site_index",
        "component_id",
        "linkage_rule",
    }
    rows = read_tsv(membership_path, compressed=True)
    if rows:
        missing = sorted(required - set(rows[0]))
        if missing:
            raise RunnerError(f"legacy/non-strict membership is forbidden; missing {missing}")
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen_owner: set[tuple[str, str, str]] = set()
    for line_number, row in enumerate(rows, start=2):
        if row["dataset"] != sample or row["chrom"] != chrom:
            raise RunnerError(f"membership scope mismatch at line {line_number}")
        if row["linkage_rule"] != "strict_fixed_ra_endpoint_pair":
            raise RunnerError(f"legacy linkage rule at line {line_number}")
        if row["linkage_basis"] not in {"PS_HP1", "PS_HP2"} or not row["phase_set"]:
            raise RunnerError(f"non-exact PS/HP membership at line {line_number}")
        key = (row["linkage_basis"], row["phase_set"], row["site_index"], row["threshold"])
        if key in seen_owner:
            raise RunnerError(f"duplicate site ownership at line {line_number}")
        seen_owner.add(key)
        if int(row["threshold"]) == threshold:
            groups[row["component_id"]].append(row)

    counts: Counter[str] = Counter()
    for component_id, component in groups.items():
        k = len(component)
        roles = {row["inference_role"] for row in component}
        eligibility = {row["tree_eligible"] for row in component}
        classes = {row["component_class"] for row in component}
        if len(roles) != 1 or len(eligibility) != 1 or len(classes) != 1:
            raise RunnerError(f"component contract is internally inconsistent: {component_id}")
        role = next(iter(roles))
        if k == 1:
            if (
                role != "ABSTAIN_SINGLETON_UNLINKED"
                or eligibility != {"false"}
                or classes != {"UNLINKED_SINGLETON_COMPONENT"}
            ):
                raise RunnerError(f"singleton incorrectly enters tree inference: {component_id}")
            counts["singleton_abstain_regions"] += 1
        else:
            if (
                role != "PRIMARY_PS_AWARE"
                or eligibility != {"true"}
                or classes != {"READ_LINKED_MULTISITE_REGION"}
            ):
                raise RunnerError(f"multisite strict component is not PRIMARY: {component_id}")
            counts["primary_multisite_regions"] += 1
            counts["primary_site_memberships"] += k
    counts["threshold_regions_all"] = len(groups)
    if counts["threshold_regions_all"] != (
        counts["singleton_abstain_regions"] + counts["primary_multisite_regions"]
    ):
        raise RunnerError("strict component classification mass is not conserved")
    return dict(counts)


def validate_strict_receipt(
    directory: Path,
    *,
    sample: str,
    chrom: str,
) -> tuple[dict[str, Any], dict[str, int]]:
    receipt_path = directory / "receipt.json"
    validate_sha256_sidecar(receipt_path)
    receipt = read_json_strict(receipt_path)
    if (
        receipt.get("schema_name") != "intersubmod.strict_exact_ps_hp_regions"
        or receipt.get("schema_version") != "1.1.0"
        or receipt.get("all_pass") is not True
        or receipt.get("scope", {}).get("dataset") != sample
        or receipt.get("scope", {}).get("chrom") != chrom
    ):
        raise RunnerError(f"strict region receipt schema/scope failure: {receipt_path}")
    contract = receipt.get("contract") or {}
    parameters = receipt.get("parameters") or {}
    checks = receipt.get("checks") or {}
    if (
        contract.get("container") != "chromosome x exact nonmissing PS x HP1/HP2"
        or contract.get("region") != "connected component of endpoint edges with support >= threshold"
        or parameters.get("primary_threshold") != PRIMARY_THRESHOLD
        or parameters.get("linkage_rule") != "strict_fixed_ra_endpoint_pair"
        or not checks
        or any(value is not True for value in checks.values())
    ):
        raise RunnerError(f"strict region contract failure: {receipt_path}")
    membership_path = directory / f"{sample}.{chrom}.site_component_membership.tsv.gz"
    declared = (receipt.get("outputs") or {}).get("membership")
    if not isinstance(declared, dict):
        raise RunnerError("strict receipt does not bind membership output")
    validate_output_identity(declared, membership_path, "strict membership")
    return receipt, validate_strict_membership(
        membership_path, sample=sample, chrom=chrom
    )


def compare_strict_python_cpp(
    *,
    strict_dir: Path,
    cpp_dir: Path,
    sample: str,
    chrom: str,
) -> dict[str, Any]:
    python_edges = read_tsv(strict_dir / f"{sample}.{chrom}.endpoint_edges.tsv.gz", compressed=True)
    cpp_edges = read_tsv(cpp_dir / "edges.tsv", compressed=False)

    def py_edge(row: Mapping[str, str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
        key = (
            row["dataset"], row["chrom"], row["phase_set"], row["hp_family"],
            row["site_i_index"], row["site_j_index"], row["pos_i1"], row["pos_j1"],
        )
        value = (
            row["support_total"], row["support_RR"], row["support_RA"],
            row["support_AR"], row["support_AA"],
            "1" if row["passes_primary_threshold"] == "true" else "0",
        )
        return key, value

    def cpp_edge(row: Mapping[str, str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
        hp = {"HP1": "1", "HP2": "2"}.get(row["hp_family"], row["hp_family"])
        key = (
            row["dataset"], row["chrom"], row["phase_set"], hp,
            row["left_site_index"], row["right_site_index"], row["left_pos1"], row["right_pos1"],
        )
        value = (
            row["support_total"], row["support_RR"], row["support_RA"],
            row["support_AR"], row["support_AA"], row["qualifies"],
        )
        return key, value

    py_edge_map = dict(py_edge(row) for row in python_edges)
    cpp_edge_map = dict(cpp_edge(row) for row in cpp_edges)
    if len(py_edge_map) != len(python_edges) or len(cpp_edge_map) != len(cpp_edges):
        raise RunnerError("duplicate endpoint edge key in strict parity inputs")

    membership = read_tsv(
        strict_dir / f"{sample}.{chrom}.site_component_membership.tsv.gz",
        compressed=True,
    )
    component_meta = {
        row["component_id"]: row
        for row in read_tsv(
            strict_dir / f"{sample}.{chrom}.components.tsv.gz", compressed=True
        )
        if int(row["threshold"]) == PRIMARY_THRESHOLD
    }
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in membership:
        if int(row["threshold"]) == PRIMARY_THRESHOLD:
            grouped[row["component_id"]].append(row)
    py_components: dict[tuple[str, ...], tuple[str, ...]] = {}
    for component_id, rows in grouped.items():
        rows.sort(key=lambda row: int(row["site_index"]))
        first = rows[0]
        meta = component_meta[component_id]
        hp = "1" if first["linkage_basis"] == "PS_HP1" else "2"
        key = (
            sample, chrom, first["phase_set"], hp,
            ",".join(row["site_index"] for row in rows),
            ",".join(row["pos1"] for row in rows),
        )
        py_components[key] = (
            str(len(rows)), meta["retained_edge_count"], meta["component_class"],
            meta["tree_eligible"], meta["inference_role"],
        )
    cpp_components: dict[tuple[str, ...], tuple[str, ...]] = {}
    for row in read_tsv(cpp_dir / "components.tsv", compressed=False):
        hp = {"HP1": "1", "HP2": "2"}.get(row["hp_family"], row["hp_family"])
        key = (
            row["dataset"], row["chrom"], row["phase_set"], hp,
            row["site_indices"], row["positions1"],
        )
        cpp_components[key] = (
            row["k"], row["qualifying_edge_count"], row["component_class"],
            row["tree_eligible"], row["inference_role"],
        )
    edge_mismatch = py_edge_map != cpp_edge_map
    component_mismatch = py_components != cpp_components
    receipt = read_json_strict(cpp_dir / "receipt.json")
    invariants = receipt.get("invariants") or {}
    cpp_receipt_pass = (
        receipt.get("schema_name") == "intersubmod.strict_endpoint_graph_receipt"
        and receipt.get("schema_version") == "1.1.0"
        and receipt.get("threshold") == PRIMARY_THRESHOLD
        and invariants
        and all(value is True for value in invariants.values())
    )
    result = {
        "schema_name": "intersubmod.strict_endpoint_graph_python_cpp_comparison",
        "schema_version": "1.0.0",
        "all_pass": cpp_receipt_pass and not edge_mismatch and not component_mismatch,
        "mismatch_count": int(edge_mismatch) + int(component_mismatch),
        "checks": {
            "cpp_receipt_contract_pass": cpp_receipt_pass,
            "endpoint_edges_exact": not edge_mismatch,
            "components_exact": not component_mismatch,
            "singleton_abstention_exact": not component_mismatch,
        },
        "counts": {
            "python_edges": len(py_edge_map),
            "cpp_edges": len(cpp_edge_map),
            "python_components": len(py_components),
            "cpp_components": len(cpp_components),
        },
        "inputs": {
            "python_edges": file_identity(
                strict_dir / f"{sample}.{chrom}.endpoint_edges.tsv.gz"
            ),
            "python_components": file_identity(
                strict_dir / f"{sample}.{chrom}.components.tsv.gz"
            ),
            "python_membership": file_identity(
                strict_dir / f"{sample}.{chrom}.site_component_membership.tsv.gz"
            ),
            "cpp_edges": file_identity(cpp_dir / "edges.tsv"),
            "cpp_components": file_identity(cpp_dir / "components.tsv"),
            "cpp_receipt": file_identity(cpp_dir / "receipt.json"),
        },
    }
    if not result["all_pass"]:
        raise RunnerError(f"strict Python/C++ semantic parity failed: {result}")
    return result


def decompress_partition_inputs(python_dir: Path, normalized_dir: Path) -> None:
    normalized_dir.mkdir(parents=True, exist_ok=False)
    for name in ("units", "constraints"):
        source = python_dir / f"{name}.tsv.gz"
        target = normalized_dir / f"{name}.tsv"
        with gzip.open(source, "rb") as reader, target.open("xb") as writer:
            shutil.copyfileobj(reader, writer, length=1 << 20)
            writer.flush()
            os.fsync(writer.fileno())


def validate_partition_receipt(
    path: Path, *, sample: str, chrom: str, strict_membership: Path
) -> dict[str, Any]:
    receipt = read_json_strict(path)
    if (
        receipt.get("schema_name") != "intersubmod.exact_ps_k12_partition"
        or receipt.get("schema_version") != "0.1.0"
        or receipt.get("all_pass") is not True
        or receipt.get("scope") != {"dataset": sample, "chrom": chrom}
        or receipt.get("parameters", {}).get("accepted_inference_role") != "PRIMARY_PS_AWARE"
    ):
        raise RunnerError(f"partition receipt contract failed: {path}")
    checks = receipt.get("checks") or {}
    if not checks or any(value is not True for value in checks.values()):
        raise RunnerError(f"partition receipt has failed check: {path}")
    declared = receipt.get("inputs", {}).get("site_component_membership") or {}
    observed = file_identity(strict_membership)
    declared_size = declared.get("size_bytes", declared.get("bytes"))
    if declared.get("sha256") != observed["sha256"] or declared_size != observed["size_bytes"]:
        raise RunnerError("partition did not consume the strict membership artifact")
    return receipt


def compile_binaries(args: argparse.Namespace, run_root: Path) -> tuple[Path, Path, list[dict[str, Any]]]:
    bin_dir = run_root / "bin"
    bin_dir.mkdir(parents=True, exist_ok=False)
    logs = run_root / "logs"
    strict_binary = bin_dir / "strict_endpoint_graph_verify"
    partition_binary = bin_dir / "exact_ps_k12_partition"
    records = [
        run_command(
            stage="compile_strict_cpp",
            command=[
                args.cxx, "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
                "-I", str(args.strict_cpp.parent), str(args.strict_cpp), "-o", str(strict_binary),
            ],
            logs_dir=logs,
        )
    ]
    if STAGES[args.stage_through] >= STAGES["partition"]:
        records.append(
            run_command(
                stage="compile_partition_cpp",
                command=[
                    args.cxx, "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
                    str(args.partition_cpp), "-o", str(partition_binary),
                ],
                logs_dir=logs,
            )
        )
    return strict_binary, partition_binary, records


def process_chromosome(
    *,
    args: argparse.Namespace,
    sample: str,
    chrom: str,
    entry: Mapping[str, Any],
    sample_root: Path,
    strict_binary: Path,
    partition_binary: Path,
) -> dict[str, Any]:
    chrom_root = sample_root / "chromosomes" / chrom
    chrom_root.mkdir(parents=True, exist_ok=False)
    logs_dir = chrom_root / "logs"
    stages: list[dict[str, Any]] = []
    extraction, extraction_stage = materialize_extraction(
        args=args,
        sample=sample,
        chrom=chrom,
        entry=entry,
        chrom_root=chrom_root,
        logs_dir=logs_dir,
    )
    stages.append(extraction_stage)
    prefix = f"{sample}.{chrom}"
    site_catalog = extraction / f"{prefix}.site_catalog.tsv.gz"
    molecule_calls = extraction / f"{prefix}.molecule_sparse_calls.tsv.gz"

    strict_dir = chrom_root / "strict_regions"
    stages.append(
        run_command(
            stage="strict_regions",
            command=[
                args.python,
                str(args.region_builder),
                "--dataset", sample,
                "--chrom", chrom,
                "--site-catalog", str(site_catalog),
                "--molecule-calls", str(molecule_calls),
                "--output-dir", str(strict_dir),
                "--thresholds", "1,2,3,5",
                "--primary-threshold", str(PRIMARY_THRESHOLD),
            ],
            logs_dir=logs_dir,
        )
    )
    strict_receipt, strict_counts = validate_strict_receipt(
        strict_dir, sample=sample, chrom=chrom
    )

    strict_cpp_dir = chrom_root / "strict_cpp"
    stages.append(
        run_streamed_strict_cpp(
            gzip_executable=args.gzip,
            binary=strict_binary,
            molecule_calls=molecule_calls,
            output_dir=strict_cpp_dir,
            logs_dir=logs_dir,
        )
    )
    strict_comparison = compare_strict_python_cpp(
        strict_dir=strict_dir,
        cpp_dir=strict_cpp_dir,
        sample=sample,
        chrom=chrom,
    )
    strict_comparison_path = chrom_root / "strict_comparison.json"
    write_json_exclusive(strict_comparison_path, strict_comparison)
    write_sha256_sidecar(strict_comparison_path)

    metrics: Counter[str] = Counter(strict_counts)
    metrics["S"] = int(strict_receipt["scope"]["candidate_loci_S"])
    primary_summary = strict_receipt.get("summary_by_threshold", {}).get(str(PRIMARY_THRESHOLD), {})
    metrics["strict_active_node_memberships"] = int(primary_summary.get("active_node_memberships", 0))
    metrics["strict_retained_endpoint_edges"] = int(primary_summary.get("retained_endpoint_edges", 0))
    metrics["python_cpp_strict_mismatch_count"] = strict_comparison["mismatch_count"]
    partition_comparison_identity = None

    if STAGES[args.stage_through] >= STAGES["partition"]:
        strict_membership = strict_dir / f"{prefix}.site_component_membership.tsv.gz"
        python_dir = chrom_root / "python_partition"
        stages.append(
            run_command(
                stage="python_partition",
                command=[
                    args.python, str(args.partitioner),
                    "--dataset", sample,
                    "--chrom", chrom,
                    "--site-catalog", str(site_catalog),
                    "--site-component-membership", str(strict_membership),
                    "--molecule-calls", str(molecule_calls),
                    "--output-dir", str(python_dir),
                    "--threshold", str(PRIMARY_THRESHOLD),
                    "--max-block-size", str(MAX_BLOCK_SIZE),
                ],
                logs_dir=logs_dir,
            )
        )
        partition_receipt = validate_partition_receipt(
            python_dir / "receipt.json",
            sample=sample,
            chrom=chrom,
            strict_membership=strict_membership,
        )
        normalized_dir = chrom_root / "normalized_cpp_inputs"
        decompress_partition_inputs(python_dir, normalized_dir)
        cpp_dir = chrom_root / "cpp_partition"
        stages.append(
            run_command(
                stage="cpp_partition",
                command=[
                    str(partition_binary),
                    "--units", str(normalized_dir / "units.tsv"),
                    "--constraints", str(normalized_dir / "constraints.tsv"),
                    "--output-dir", str(cpp_dir),
                    "--max-block-size", str(MAX_BLOCK_SIZE),
                ],
                logs_dir=logs_dir,
            )
        )
        comparison_path = chrom_root / "comparison.json"
        stages.append(
            run_command(
                stage="partition_python_cpp_comparison",
                command=[
                    args.python, str(args.partition_comparator),
                    "--python-dir", str(python_dir),
                    "--cpp-input-units", str(normalized_dir / "units.tsv"),
                    "--cpp-input-constraints", str(normalized_dir / "constraints.tsv"),
                    "--cpp-dir", str(cpp_dir),
                    "--output", str(comparison_path),
                ],
                logs_dir=logs_dir,
            )
        )
        comparison = read_json_strict(comparison_path)
        if comparison.get("all_pass") is not True or comparison.get("mismatch_count") != 0:
            raise RunnerError(f"partition Python/C++ comparison failed: {comparison_path}")
        partition_counts = partition_receipt.get("counts") or {}
        metrics["units"] = int(partition_counts.get("eligible_units", 0))
        metrics["unit_memberships"] = int(partition_counts.get("eligible_unit_sites", 0))
        metrics["blocks"] = int(partition_counts.get("blocks", 0))
        metrics["python_cpp_mismatch_count"] = int(comparison["mismatch_count"])
        partition_comparison_identity = file_identity(comparison_path)
        if metrics["units"] != metrics["primary_multisite_regions"]:
            raise RunnerError("partition unit count differs from strict multisite region count")

    receipt = {
        "schema_name": f"{SCHEMA_NAME}.chromosome_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "sample": sample,
        "chrom": chrom,
        "stage_through": args.stage_through,
        "contract": {
            "container": "chromosome x exact nonmissing PS x HP1/HP2",
            "region": "connected component of fixed-R/A endpoint edges supported by >=3 distinct molecules",
            "tree_input": "PRIMARY_PS_AWARE multisite strict components only (k>1)",
            "singleton": "ABSTAIN_SINGLETON_UNLINKED; excluded before partition",
        },
        "metrics": dict(metrics),
        "stages": stages,
        "strict_receipt": file_identity(strict_dir / "receipt.json"),
        "strict_comparison": file_identity(strict_comparison_path),
        "partition_comparison": partition_comparison_identity,
    }
    receipt_path = chrom_root / "stage_receipt.json"
    write_json_exclusive(receipt_path, receipt)
    write_sha256_sidecar(receipt_path)
    return receipt


def optional_manifest_path(value: Any) -> Path | None:
    if isinstance(value, dict) and isinstance(value.get("path"), str):
        candidate = Path(value["path"])
        return candidate if candidate.is_file() else None
    return None


def run_topology(
    *,
    args: argparse.Namespace,
    sample: str,
    entry: Mapping[str, Any],
    sample_root: Path,
) -> dict[str, Any]:
    topology_dir = sample_root / "topology"
    topology_dir.mkdir(parents=True, exist_ok=False)
    logs_dir = topology_dir / "logs"
    exact_input = topology_dir / f"{sample}.strict_exact_ps_mlhp.json"
    adapter_command = [
        args.python, str(args.adapter),
        "--partition-root", str(sample_root),
        "--output", str(exact_input),
        "--sample", sample,
        "--chroms", ",".join(args.chromosomes),
        "--min-read", "3",
        "--require-strict-endpoint-receipt",
    ]
    copy_number = entry.get("copy_number") if isinstance(entry, dict) else None
    cn_bed = optional_manifest_path(copy_number.get("cn_bed")) if isinstance(copy_number, dict) else None
    if cn_bed is not None:
        adapter_command.extend(["--cn-bed", str(cn_bed), "--cn-source", str(copy_number.get("source", "manifest"))])
    stages = [run_command(stage="exact_ps_adapter", command=adapter_command, logs_dir=logs_dir)]
    adapter_receipt = read_json_strict(exact_input.with_suffix(exact_input.suffix + ".receipt.json"))
    if adapter_receipt.get("all_pass") is not True:
        raise RunnerError("exact-PS adapter receipt is not PASS")

    layered_output = topology_dir / f"{sample}.layered_reconstruction.json"
    env = dict(os.environ)
    env.update(
        {
            "SM_ML": str(exact_input),
            "SM_OUT": str(layered_output),
            "SM_VERIFY_EVERY": "1",
            "SM_ANALYSIS_TREE_CAP": "0",
            "SM_DISPLAY_TREE_CAP": "32",
            # Explicit empties prevent the legacy HCC1395 defaults leaking into other samples.
            "SM_CN_INT_GAIN": "",
            "SM_CN_INT_LOSS": "",
        }
    )
    if isinstance(copy_number, dict):
        gain = optional_manifest_path(copy_number.get("cn_int_gain"))
        loss = optional_manifest_path(copy_number.get("cn_int_loss"))
        if gain is not None:
            env["SM_CN_INT_GAIN"] = str(gain)
        if loss is not None:
            env["SM_CN_INT_LOSS"] = str(loss)
    stages.append(
        run_command(
            stage="layered_reconstruction",
            command=[args.python, str(args.layered_reconstruction)],
            logs_dir=logs_dir,
            env=env,
        )
    )
    layered_document = read_json_strict(layered_output)
    if not isinstance(layered_document, dict) or not isinstance(layered_document.get("detail"), list):
        raise RunnerError("layered reconstruction output lacks detail list")

    region_output = topology_dir / f"{sample}.layered_region_view.json"
    region_env = dict(os.environ)
    tree_vcf = entry["somatic"]["tree_vcf"]["path"]
    region_env.update(
        {
            "SM_LAYERED": str(layered_output),
            "SM_OUT": str(region_output),
            "SM_SAMPLE": sample,
            "SM_ML_GLOB": str(exact_input),
            "SM_SOMATIC_VCF": str(tree_vcf),
            "SM_BACKBONE_SOURCE": "LongPhase-S recalibrated FILTER=PASS",
        }
    )
    stages.append(
        run_command(
            stage="region_view",
            command=[args.python, str(args.region_view)],
            logs_dir=logs_dir,
            env=region_env,
        )
    )
    region_document = read_json_strict(region_output)
    if not isinstance(region_document, dict) or not isinstance(region_document.get("regions"), list):
        raise RunnerError("region view output lacks regions list")
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.topology_receipt",
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "sample": sample,
        "stages": stages,
        "outputs": {
            "exact_tree_input": file_identity(exact_input),
            "layered_reconstruction": file_identity(layered_output),
            "region_view": file_identity(region_output),
        },
        "counts": {
            "adapter_groups": int(read_json_strict(exact_input).get("n_groups_analyzed", 0)),
            "layered_units": len(layered_document["detail"]),
            "region_view_regions": len(region_document["regions"]),
        },
        "claim_ceiling": "Candidate regional mutation-state trees; not a unique cellular clone tree.",
    }
    receipt_path = topology_dir / "receipt.json"
    write_json_exclusive(receipt_path, receipt)
    write_sha256_sidecar(receipt_path)
    return receipt


def execute(args: argparse.Namespace) -> dict[str, Any]:
    full_scope = validate_scope(args)
    args.manifest = require_tool(args.manifest, "manifest")
    entries, _ = load_manifest(args.manifest, allow_partial=not full_scope)
    missing_requested = sorted(set(args.datasets) - set(entries))
    if missing_requested:
        raise RunnerError(f"requested datasets absent from manifest: {missing_requested}")
    args.python = require_python(args.python)
    args.cxx = require_executable(args.cxx, "C++ compiler")
    args.gzip = require_executable(args.gzip, "gzip")
    args.samtools = require_executable(args.samtools, "samtools")
    for attribute in (
        "extractor", "region_builder", "partitioner", "partition_cpp",
        "partition_comparator", "strict_cpp", "strict_cpp_header", "adapter",
        "layered_reconstruction", "region_view",
    ):
        setattr(args, attribute, require_tool(getattr(args, attribute), attribute))
    run_root = prepare_output_root(args.output_root)
    manifest_before = file_identity(args.manifest)
    source_paths = {
        name: getattr(args, name)
        for name in (
            "extractor", "region_builder", "partitioner", "partition_cpp",
            "partition_comparator", "strict_cpp", "strict_cpp_header", "adapter",
            "layered_reconstruction", "region_view",
        )
    }
    source_paths["launcher"] = Path(__file__).resolve()
    source_before = {name: file_identity(path) for name, path in source_paths.items()}
    try:
        strict_binary, partition_binary, compile_records = compile_binaries(args, run_root)
    except Exception as exc:
        write_failure_marker(run_root, exc)
        raise
    samples_root = run_root / "samples"
    sample_receipts: dict[str, Any] = {}
    aggregate: Counter[str] = Counter()
    try:
        samples_root.mkdir(parents=True, exist_ok=False)
        for sample in args.datasets:
            sample_root = samples_root / sample
            sample_root.mkdir(parents=True, exist_ok=False)
            chromosome_receipts: dict[str, Any] = {}
            sample_counts: Counter[str] = Counter()
            for chrom in args.chromosomes:
                receipt = process_chromosome(
                    args=args,
                    sample=sample,
                    chrom=chrom,
                    entry=entries[sample],
                    sample_root=sample_root,
                    strict_binary=strict_binary,
                    partition_binary=partition_binary,
                )
                chromosome_receipts[chrom] = file_identity(
                    sample_root / "chromosomes" / chrom / "stage_receipt.json"
                )
                sample_counts.update(receipt["metrics"])
            partition_receipt = {
                "schema_name": f"{SCHEMA_NAME}.sample_partition_receipt",
                "schema_version": SCHEMA_VERSION,
                "created_at_utc": utc_now(),
                "all_pass": True,
                "sample": sample,
                "scope": {
                    "chromosomes": list(args.chromosomes),
                    "scope_label": "chr1-22" if tuple(args.chromosomes) == AUTOSOMES else "PARTIAL",
                },
                "stage_through": args.stage_through,
                "metrics": dict(sample_counts),
                "chromosome_receipts": chromosome_receipts,
            }
            sample_run_receipt_path = sample_root / "run_receipt.json"
            write_json_exclusive(sample_run_receipt_path, partition_receipt)
            write_sha256_sidecar(sample_run_receipt_path)
            topology_receipt = None
            if args.stage_through == "topology":
                topology_receipt = run_topology(
                    args=args,
                    sample=sample,
                    entry=entries[sample],
                    sample_root=sample_root,
                )
            sample_receipts[sample] = {
                "partition": file_identity(sample_run_receipt_path),
                "topology": (
                    file_identity(sample_root / "topology" / "receipt.json")
                    if topology_receipt is not None
                    else None
                ),
            }
            aggregate.update(sample_counts)

        if file_identity(args.manifest) != manifest_before:
            raise RunnerError("manifest identity changed during the run")
        source_after = {name: file_identity(path) for name, path in source_paths.items()}
        if source_after != source_before:
            raise RunnerError("executable source identity changed during the run")
        receipt = {
            "schema_name": SCHEMA_NAME,
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": utc_now(),
            "all_pass": True,
            "task_type": "production_deployment" if full_scope else "exploratory_partial",
            "claim_status": "FULL" if full_scope else "PARTIAL",
            "validation_evidence_eligible": full_scope,
            "stage_through": args.stage_through,
            "scope": {"datasets": list(args.datasets), "chromosomes": list(args.chromosomes)},
            "parameters": {
                "mapq_min": 20,
                "baseq_min": 20,
                "endpoint_edge_min_distinct_molecules": PRIMARY_THRESHOLD,
                "max_block_size": MAX_BLOCK_SIZE,
            },
            "contract": {
                "container": "chromosome x exact nonmissing PS x HP1/HP2",
                "node": "sSNV with fixed R/A observation in the exact container",
                "region": "connected component of retained fixed-R/A endpoint edges",
                "tree_input": "strict multisite component only (k>1)",
                "singleton": "funnel-only abstention; never a tree region",
            },
            "manifest": manifest_before,
            "sources": source_before,
            "compile_stages": compile_records,
            "sample_receipts": sample_receipts,
            "aggregate": dict(aggregate),
            "invocation": args.invocation,
            "claim_ceiling": (
                "The workflow produces candidate local mutation-state trees. Exact PS and read "
                "connectivity do not by themselves prove a unique cellular clone tree."
            ),
        }
        receipt_path = run_root / "run_receipt.json"
        write_json_exclusive(receipt_path, receipt)
        write_sha256_sidecar(receipt_path)
        success = {
            "schema_name": f"{SCHEMA_NAME}.success",
            "schema_version": SCHEMA_VERSION,
            "all_pass": True,
            "run_receipt": file_identity(receipt_path),
        }
        write_json_exclusive(run_root / "_SUCCESS", success)
        return receipt
    except Exception as exc:
        write_failure_marker(run_root, exc)
        raise


def main(argv: Sequence[str] | None = None) -> int:
    try:
        receipt = execute(parse_args(argv))
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "all_pass": receipt["all_pass"],
                "scope": receipt["scope"],
                "stage_through": receipt["stage_through"],
                "aggregate": receipt["aggregate"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
