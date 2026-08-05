#!/usr/bin/env python3
"""Run a cached physical-span sensitivity grid after HCC1395 full completion.

The source HCC1395 full run is immutable input.  Before creating an output
root, this runner authenticates the source success marker, aggregate receipt,
run contract, every requested chromosome extraction stage/child receipt, and
every extraction child output identity.  It never invokes the BAM extractor.

For each physical-span cap and chromosome with a k>8 target, the independent
``build_k_gt8_partitions_span_cap.py`` CLI is run under ``/usr/bin/time -v``.
The zero-target chromosome (canonical chr21) receives an explicit skip receipt
and time file.  Fresh and resume modes are fail-closed: incomplete, altered,
or unauthenticated artifacts are never deleted or overwritten.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import run_hcc1395_full_segmentation as full_contract


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PARTITIONER = SCRIPT_DIR / "build_k_gt8_partitions_span_cap.py"
TIME_BINARY = Path("/usr/bin/time")
AUTOSOMES = tuple(f"chr{value}" for value in range(1, 23))
CANONICAL_SPAN_CAPS = (50_000, 100_000, 200_000)
DATASET = "HCC1395"
SCHEMA_NAME = "intersubmod.hcc1395_cached_span_sensitivity"
SCHEMA_VERSION = "1.0.0"
SOURCE_CONTRACT_SCHEMA = (
    "intersubmod.hcc1395_full_k_gt8_segmentation.run_contract"
)
SOURCE_RECEIPT_SCHEMA = (
    "intersubmod.hcc1395_full_k_gt8_segmentation.run_receipt"
)
SOURCE_MARKER_SCHEMA = (
    "intersubmod.hcc1395_full_k_gt8_segmentation.success_marker"
)
SOURCE_STAGE_SCHEMA = (
    "intersubmod.hcc1395_full_k_gt8_segmentation.stage_receipt"
)
SOURCE_CHILD_SCHEMA = "intersubmod.lossless_read_linkage_chromosome_receipt"
SPAN_CHILD_SCHEMA = (
    "intersubmod.k_gt8_read_supported_segmentation_span_cap"
)
CANONICAL_TOTALS = {
    "chromosomes": 22,
    "ssnv_sites": 79_687,
    "k_gt8_components": 408,
    "k_gt8_sites": 47_570,
    "k_gt8_max_k": 3_574,
    "partitioned_chromosomes": 21,
    "zero_target_skipped_chromosomes": 1,
}

ContractError = full_contract.ContractError
file_identity = full_contract.file_identity
strict_json_load = full_contract.strict_json_load
verify_sha256_sidecar = full_contract.verify_sha256_sidecar
verify_output_identities = full_contract.verify_output_identities
derive_site_catalog_inventory = full_contract.derive_site_catalog_inventory
write_json_atomic = full_contract.write_json_atomic
write_sha256_sidecar = full_contract.write_sha256_sidecar
sha256_path = full_contract.sha256_path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def require_regular_file(
    path: Path,
    label: str,
    *,
    allow_symlink: bool = False,
) -> None:
    if not path.is_file() or (path.is_symlink() and not allow_symlink):
        raise ContractError(f"{label} is missing, symlinked, or not regular: {path}")


def verify_sha256_sidecar_strict(path: Path, label: str) -> str:
    require_regular_file(path, label)
    require_regular_file(
        path.with_name(f"{path.name}.sha256"),
        f"{label} SHA-256 sidecar",
    )
    return verify_sha256_sidecar(path)


def verify_identity_record(
    raw: Any,
    *,
    expected_path: Path,
    label: str,
) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ContractError(f"{label} identity is not an object")
    if set(raw) != {"path", "size_bytes", "sha256"}:
        raise ContractError(f"{label} identity keys differ")
    require_regular_file(expected_path, label)
    if Path(str(raw["path"])).resolve() != expected_path.resolve():
        raise ContractError(f"{label} identity path differs")
    observed = file_identity(expected_path)
    if observed != dict(raw):
        raise ContractError(
            f"{label} identity drift: expected={dict(raw)} observed={observed}"
        )
    return observed


def verify_sidecar_bound_identity(
    raw: Any,
    *,
    expected_path: Path,
    label: str,
) -> dict[str, Any]:
    observed = verify_identity_record(
        raw,
        expected_path=expected_path,
        label=label,
    )
    sidecar_sha = verify_sha256_sidecar_strict(expected_path, label)
    if sidecar_sha != observed["sha256"]:
        raise ContractError(f"{label} sidecar and identity SHA differ")
    return observed


def require_true_checks(receipt: Mapping[str, Any], label: str) -> None:
    checks = receipt.get("checks")
    if not isinstance(checks, Mapping) or not checks:
        raise ContractError(f"{label} has no checks")
    failed = [name for name, value in checks.items() if value is not True]
    if failed:
        raise ContractError(f"{label} failed checks: {failed}")


def parse_chromosomes(value: str) -> tuple[str, ...]:
    chromosomes = tuple(
        token.strip() for token in value.split(",") if token.strip()
    )
    if not chromosomes or len(chromosomes) != len(set(chromosomes)):
        raise argparse.ArgumentTypeError(
            "chromosomes must be nonempty and unique"
        )
    invalid = [chrom for chrom in chromosomes if chrom not in AUTOSOMES]
    if invalid:
        raise argparse.ArgumentTypeError(f"invalid autosomes: {invalid}")
    return chromosomes


def parse_span_caps(value: str) -> tuple[int, ...]:
    try:
        caps = tuple(
            int(token.strip())
            for token in value.split(",")
            if token.strip()
        )
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "span caps must be comma-separated integers"
        ) from exc
    if not caps or len(caps) != len(set(caps)):
        raise argparse.ArgumentTypeError(
            "span caps must be nonempty and unique"
        )
    if any(cap < 0 for cap in caps):
        raise argparse.ArgumentTypeError("span caps must be non-negative")
    return caps


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--source-full-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--span-caps",
        type=parse_span_caps,
        default=CANONICAL_SPAN_CAPS,
        help=(
            "comma-separated inclusive block end-start caps in bp; "
            "production scope is fixed to 50000,100000,200000"
        ),
    )
    parser.add_argument(
        "--chromosomes",
        type=parse_chromosomes,
        default=AUTOSOMES,
        help="comma-separated autosomes; subset allowed only in --test-mode",
    )
    parser.add_argument("--partitioner", type=Path, default=DEFAULT_PARTITIONER)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--test-mode",
        action="store_true",
        help=(
            "accept an authenticated source _TEST_SUCCESS run and chromosome "
            "subset; writes _TEST_SUCCESS, never production _SUCCESS"
        ),
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.test_mode and args.chromosomes != AUTOSOMES:
        raise ContractError(
            "production cached sensitivity scope must be exact chr1-22"
        )
    if not args.test_mode and args.span_caps != CANONICAL_SPAN_CAPS:
        raise ContractError(
            "production cached sensitivity caps must be exact "
            "50000,100000,200000 in canonical order"
        )
    if not args.source_full_root.is_dir() or args.source_full_root.is_symlink():
        raise ContractError(
            f"source full root must be a real directory: {args.source_full_root}"
        )
    require_regular_file(args.partitioner, "span-cap partitioner")
    require_regular_file(
        args.python,
        "Python executable",
        allow_symlink=True,
    )
    require_regular_file(TIME_BINARY, "/usr/bin/time")
    require_regular_file(
        Path(full_contract.__file__).resolve(),
        "source contract verifier module",
    )
    if args.output_root.exists() and args.output_root.is_symlink():
        raise ContractError("output root must not be a symlink")
    source = args.source_full_root.resolve()
    # strict=False still resolves any existing symlinked parent components,
    # preventing a seemingly external output path from writing into source.
    output = args.output_root.resolve(strict=False)
    try:
        output.relative_to(source)
    except ValueError:
        pass
    else:
        raise ContractError("output root must not be inside source full root")
    try:
        source.relative_to(output)
    except ValueError:
        pass
    else:
        raise ContractError("source full root must not be inside output root")


def verify_source_extraction(
    *,
    source_root: Path,
    chrom: str,
    inventory: Mapping[str, Any],
    source_contract: Mapping[str, Any],
    source_contract_sha256: str,
    aggregate_entry: Mapping[str, Any],
) -> dict[str, Any]:
    chrom_root = source_root / "chromosomes" / chrom
    extraction_dir = chrom_root / "extract"
    stage_path = chrom_root / "extract.stage_receipt.json"
    stage_identity = verify_sidecar_bound_identity(
        aggregate_entry.get("extraction"),
        expected_path=stage_path,
        label=f"{chrom} extraction stage receipt",
    )
    stage = strict_json_load(stage_path)
    if not isinstance(stage, Mapping):
        raise ContractError(f"{chrom} extraction stage receipt is not an object")
    expected_stage_scalars = {
        "schema_name": SOURCE_STAGE_SCHEMA,
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "stage": "extract",
        "status": "COMPLETED",
        "exit_code": 0,
        "contract_sha256": source_contract_sha256,
    }
    for key, expected in expected_stage_scalars.items():
        if stage.get(key) != expected:
            raise ContractError(
                f"{chrom} extraction stage {key} mismatch: "
                f"expected={expected} observed={stage.get(key)}"
            )
    extractor_identity = source_contract.get("tools", {}).get("extractor")
    if stage.get("tool") != extractor_identity:
        raise ContractError(f"{chrom} extraction stage tool identity drift")
    source_parameters = source_contract["parameters"]
    observed_command = stage.get("command")
    if not isinstance(observed_command, list) or not observed_command:
        raise ContractError(f"{chrom} extraction stage command is missing")
    python_invocation = Path(str(observed_command[0]))
    if (
        not python_invocation.is_file()
        or python_invocation.resolve()
        != Path(source_contract["tools"]["python"]["path"]).resolve()
    ):
        raise ContractError(
            f"{chrom} extraction Python invocation differs from source tool"
        )
    expected_command = [
        str(python_invocation),
        source_contract["tools"]["extractor"]["path"],
        "--manifest",
        source_contract["frozen_manifest"]["path"],
        "--sample",
        DATASET,
        "--chrom",
        chrom,
        "--output-dir",
        str(extraction_dir.resolve()),
        "--mapq-min",
        str(source_parameters["mapq_min"]),
        "--baseq-min",
        str(source_parameters["baseq_min"]),
        "--bridge-thresholds",
        str(source_parameters["bridge_thresholds"]),
        "--samtools-threads",
        str(source_parameters["samtools_threads"]),
    ]
    if observed_command != expected_command:
        raise ContractError(f"{chrom} extraction stage command drift")
    expected_log_paths = {
        "stdout": chrom_root / "extract.stdout.log",
        "stderr": chrom_root / "extract.stderr.log",
        "resource": chrom_root / "extract.time_v.txt",
    }
    verified_logs = {}
    for name, path in expected_log_paths.items():
        verified_logs[name] = verify_identity_record(
            stage.get("logs", {}).get(name),
            expected_path=path,
            label=f"{chrom} extraction {name}",
        )
    expected_timed_command = [
        source_contract["tools"]["time"]["path"],
        "-v",
        "-o",
        str(expected_log_paths["resource"].resolve()),
        "--",
        *expected_command,
    ]
    if stage.get("timed_command") != expected_timed_command:
        raise ContractError(f"{chrom} extraction timed command drift")

    child_path = extraction_dir / "receipt.json"
    child_identity = verify_sidecar_bound_identity(
        stage.get("child_receipt"),
        expected_path=child_path,
        label=f"{chrom} extraction child receipt",
    )
    child = strict_json_load(child_path)
    if (
        not isinstance(child, Mapping)
        or child.get("schema_name") != SOURCE_CHILD_SCHEMA
        or child.get("all_pass") is not True
    ):
        raise ContractError(f"{chrom} extraction child receipt is not all_pass")
    expected_scope = {
        "dataset": DATASET,
        "chrom": chrom,
        "n_sSNV": inventory["ssnv_sites"],
    }
    if child.get("scope") != expected_scope:
        raise ContractError(f"{chrom} extraction child scope mismatch")
    expected_provenance = {
        "extractor": {
            "path": source_contract["tools"]["extractor"]["path"],
            "sha256": source_contract["tools"]["extractor"]["sha256"],
        },
        "manifest": {
            "path": source_contract["frozen_manifest"]["path"],
            "sha256": source_contract["frozen_manifest"]["sha256"],
        },
    }
    if child.get("provenance") != expected_provenance:
        raise ContractError(f"{chrom} extraction child provenance drift")
    require_true_checks(child, f"{chrom} extraction child")
    child_outputs = verify_output_identities(
        child.get("outputs", {}),
        extraction_dir,
    )
    prefix = f"{DATASET}.{chrom}"
    site_name = f"{prefix}.site_catalog.tsv.gz"
    calls_name = f"{prefix}.molecule_sparse_calls.tsv.gz"
    for required_name in (site_name, calls_name):
        if required_name not in child_outputs:
            raise ContractError(
                f"{chrom} extraction child lacks required output {required_name}"
            )
    site_catalog = extraction_dir / site_name
    molecule_calls = extraction_dir / calls_name
    extracted_inventory = derive_site_catalog_inventory(
        site_catalog,
        expected_chrom=chrom,
        legacy_gap_bp=int(source_contract["parameters"]["legacy_gap_bp"]),
        max_block_size=int(source_contract["parameters"]["max_block_size"]),
    )
    inventory_fields = (
        "ssnv_sites",
        "positional_components_all",
        "k_gt8_components",
        "k_gt8_sites",
        "k_gt8_max_k",
        "positions_sha256",
    )
    if extracted_inventory != {
        key: inventory[key] for key in inventory_fields
    }:
        raise ContractError(
            f"{chrom} extraction site inventory differs from source contract"
        )
    if aggregate_entry.get("inventory") != dict(inventory):
        raise ContractError(
            f"{chrom} aggregate inventory differs from source contract"
        )
    return {
        "inventory": dict(inventory),
        "stage_receipt": stage_identity,
        "child_receipt": child_identity,
        "logs": verified_logs,
        "outputs": child_outputs,
        "site_catalog": child_outputs[site_name],
        "molecule_calls": child_outputs[calls_name],
    }


def verify_source_full(
    args: argparse.Namespace,
) -> dict[str, Any]:
    source_root = args.source_full_root.resolve()
    marker_name = "_TEST_SUCCESS" if args.test_mode else "_SUCCESS"
    marker_path = source_root / marker_name
    require_regular_file(marker_path, f"source {marker_name}")
    marker = strict_json_load(marker_path)
    if not isinstance(marker, Mapping):
        raise ContractError("source success marker is not an object")
    if (
        marker.get("schema_name") != SOURCE_MARKER_SCHEMA
        or marker.get("all_pass") is not True
        or marker.get("sample") != DATASET
    ):
        raise ContractError("source success marker is not all_pass HCC1395")
    expected_comprehensive = not args.test_mode
    if marker.get("comprehensive") is not expected_comprehensive:
        raise ContractError("source success marker comprehensive flag mismatch")

    receipt_path = source_root / "receipt.json"
    receipt_sha = verify_sha256_sidecar_strict(
        receipt_path,
        "source run receipt",
    )
    receipt = strict_json_load(receipt_path)
    if not isinstance(receipt, Mapping):
        raise ContractError("source run receipt is not an object")
    expected_receipt_scalars = {
        "schema_name": SOURCE_RECEIPT_SCHEMA,
        "all_pass": True,
        "comprehensive_all_pass": expected_comprehensive,
        "sample": DATASET,
    }
    for key, expected in expected_receipt_scalars.items():
        if receipt.get(key) != expected:
            raise ContractError(f"source run receipt {key} mismatch")
    expected_marker_receipt = {
        "path": str(receipt_path.resolve()),
        "sha256": receipt_sha,
    }
    if marker.get("run_receipt") != expected_marker_receipt:
        raise ContractError("source success marker does not bind run receipt")

    contract_path = source_root / "run_contract.json"
    contract_sha = verify_sha256_sidecar_strict(
        contract_path,
        "source run contract",
    )
    contract_identity = verify_identity_record(
        receipt.get("contract"),
        expected_path=contract_path,
        label="source run contract",
    )
    if contract_identity["sha256"] != contract_sha:
        raise ContractError("source contract sidecar and receipt identity differ")
    contract = strict_json_load(contract_path)
    if (
        not isinstance(contract, Mapping)
        or contract.get("schema_name") != SOURCE_CONTRACT_SCHEMA
        or contract.get("sample") != DATASET
    ):
        raise ContractError("unexpected source run contract")

    source_scope = tuple(contract.get("scope", {}).get("chromosomes", ()))
    receipt_scope = tuple(receipt.get("scope", {}).get("chromosomes", ()))
    marker_scope = tuple(marker.get("scope", ()))
    if not source_scope or not (
        source_scope == receipt_scope == marker_scope == args.chromosomes
    ):
        raise ContractError(
            "source contract, receipt, marker, and requested chromosome "
            "scopes must match exactly"
        )
    if bool(contract.get("scope", {}).get("test_mode")) is not args.test_mode:
        raise ContractError("source contract test_mode mismatch")
    if bool(receipt.get("scope", {}).get("test_mode")) is not args.test_mode:
        raise ContractError("source receipt test_mode mismatch")
    if not args.test_mode:
        if source_scope != AUTOSOMES:
            raise ContractError("production source scope is not chr1-22")
        if contract.get("scope", {}).get("comprehensive_chr1_22") is not True:
            raise ContractError("production source is not comprehensive")
        if receipt.get("totals") != CANONICAL_TOTALS:
            raise ContractError("production source canonical totals mismatch")

    inventory_rows = contract.get("vcf_derived_inventory")
    aggregate_rows = receipt.get("chromosomes")
    if not isinstance(inventory_rows, list) or not isinstance(
        aggregate_rows, list
    ):
        raise ContractError("source contract/receipt chromosome rows missing")
    if tuple(row.get("chrom") for row in inventory_rows) != source_scope:
        raise ContractError("source contract inventory order/scope mismatch")
    if tuple(row.get("chrom") for row in aggregate_rows) != source_scope:
        raise ContractError("source aggregate chromosome order/scope mismatch")
    derived_source_totals = {
        "chromosomes": len(inventory_rows),
        "ssnv_sites": sum(int(row["ssnv_sites"]) for row in inventory_rows),
        "k_gt8_components": sum(
            int(row["k_gt8_components"]) for row in inventory_rows
        ),
        "k_gt8_sites": sum(
            int(row["k_gt8_sites"]) for row in inventory_rows
        ),
        "k_gt8_max_k": max(
            (int(row["k_gt8_max_k"]) for row in inventory_rows),
            default=0,
        ),
        "partitioned_chromosomes": sum(
            int(row["k_gt8_components"]) > 0 for row in inventory_rows
        ),
        "zero_target_skipped_chromosomes": sum(
            int(row["k_gt8_components"]) == 0 for row in inventory_rows
        ),
    }
    if receipt.get("totals") != derived_source_totals:
        raise ContractError("source aggregate totals differ from inventory")
    aggregate_outputs = verify_output_identities(
        receipt.get("outputs", {}),
        source_root,
    )

    verified_chromosomes = {}
    for inventory, aggregate in zip(inventory_rows, aggregate_rows):
        if not isinstance(inventory, Mapping) or not isinstance(
            aggregate, Mapping
        ):
            raise ContractError("source chromosome entry is not an object")
        chrom = str(inventory["chrom"])
        if not args.test_mode:
            expected = full_contract.CANONICAL_EXPECTED[chrom]
            observed = (
                inventory.get("ssnv_sites"),
                inventory.get("k_gt8_components"),
                inventory.get("k_gt8_sites"),
                inventory.get("k_gt8_max_k"),
            )
            if observed != expected:
                raise ContractError(
                    f"source canonical inventory mismatch for {chrom}"
                )
        verified_chromosomes[chrom] = verify_source_extraction(
            source_root=source_root,
            chrom=chrom,
            inventory=inventory,
            source_contract=contract,
            source_contract_sha256=contract_sha,
            aggregate_entry=aggregate,
        )

    return {
        "root": str(source_root),
        "success_marker": file_identity(marker_path),
        "run_receipt": file_identity(receipt_path),
        "run_contract": contract_identity,
        "aggregate_outputs": aggregate_outputs,
        "parameters": {
            "legacy_gap_bp": int(contract["parameters"]["legacy_gap_bp"]),
            "max_block_size": int(contract["parameters"]["max_block_size"]),
        },
        "chromosomes": verified_chromosomes,
    }


def establish_bytes(path: Path, payload: bytes, *, resume: bool) -> None:
    if path.exists() or path.is_symlink():
        if not resume:
            raise ContractError(f"output exists without --resume: {path}")
        verify_sha256_sidecar_strict(path, "resumed byte output")
        if path.read_bytes() != payload:
            raise ContractError(f"existing output differs on resume: {path}")
        return
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    write_sha256_sidecar(path)


def establish_json_stable(
    path: Path,
    value: Mapping[str, Any],
    *,
    resume: bool,
    stable_keys: Sequence[str],
) -> tuple[dict[str, Any], dict[str, Any]]:
    if path.exists() or path.is_symlink():
        if not resume:
            raise ContractError(f"JSON output exists without --resume: {path}")
        sha = verify_sha256_sidecar_strict(path, "resumed JSON output")
        observed = strict_json_load(path)
        if not isinstance(observed, Mapping):
            raise ContractError(f"existing JSON is not an object: {path}")
        if {key: observed.get(key) for key in stable_keys} != {
            key: value.get(key) for key in stable_keys
        }:
            raise ContractError(f"existing JSON differs on resume: {path}")
        return dict(observed), {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha,
        }
    write_json_atomic(path, value)
    sha = sha256_path(path)
    write_sha256_sidecar(path)
    return dict(value), {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha,
    }


def build_run_contract(
    args: argparse.Namespace,
    source: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema_name": f"{SCHEMA_NAME}.run_contract",
        "schema_version": SCHEMA_VERSION,
        "sample": DATASET,
        "scope": {
            "chromosomes": list(args.chromosomes),
            "span_caps_bp": list(args.span_caps),
            "test_mode": args.test_mode,
            "comprehensive_chr1_22": (
                not args.test_mode and args.chromosomes == AUTOSOMES
            ),
        },
        "source_full_validation": source,
        "tools": {
            "runner": file_identity(Path(__file__).resolve()),
            "source_contract_verifier": file_identity(
                Path(full_contract.__file__).resolve()
            ),
            "partitioner": file_identity(args.partitioner),
            "python": file_identity(args.python),
            "time": file_identity(TIME_BINARY),
        },
        "parameters": {
            "legacy_gap_bp": source["parameters"]["legacy_gap_bp"],
            "max_block_size": source["parameters"]["max_block_size"],
            "cached_extraction_only": True,
            "sequential_cap_chromosome_order": True,
            "scheduling_policy": "inherited; apply nice/ionice externally",
            "zero_target_policy": (
                "write explicit skip receipt and time file; never invoke "
                "partitioner"
            ),
        },
    }


def establish_run_contract(
    args: argparse.Namespace,
    contract: Mapping[str, Any],
) -> tuple[Path, str]:
    output_root = args.output_root
    contract_path = output_root / "run_contract.json"
    if output_root.exists():
        if not args.resume:
            raise ContractError(
                f"output root exists; use --resume only for exact contract: "
                f"{output_root}"
            )
        if not output_root.is_dir() or output_root.is_symlink():
            raise ContractError("resume output root is not a real directory")
        sha = verify_sha256_sidecar_strict(
            contract_path,
            "cached sensitivity run contract",
        )
        observed = strict_json_load(contract_path)
        if observed != dict(contract):
            raise ContractError("resume run contract differs")
        return contract_path, sha
    if args.resume:
        raise ContractError(
            f"--resume output root does not exist: {output_root}"
        )
    output_root.mkdir(parents=True, exist_ok=False)
    write_json_atomic(contract_path, contract)
    sha = sha256_path(contract_path)
    write_sha256_sidecar(contract_path)
    return contract_path, sha


def cap_label(cap: int) -> str:
    return f"span_{cap}bp"


def stage_paths(chrom_root: Path) -> dict[str, Path]:
    return {
        "output": chrom_root / "partition",
        "stdout": chrom_root / "partition.stdout.log",
        "stderr": chrom_root / "partition.stderr.log",
        "resource": chrom_root / "partition.time_v.txt",
        "stage_receipt": chrom_root / "partition.stage_receipt.json",
        "skip_receipt": chrom_root / "partition.skip_receipt.json",
        "failure": chrom_root / "partition.failure.json",
    }


def expected_partition_command(
    args: argparse.Namespace,
    *,
    cap: int,
    chrom: str,
    output_dir: Path,
    source_chrom: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> list[str]:
    inventory = source_chrom["inventory"]
    return [
        str(args.python.resolve()),
        str(args.partitioner.resolve()),
        "--dataset",
        DATASET,
        "--chrom",
        chrom,
        "--site-catalog",
        source_chrom["site_catalog"]["path"],
        "--molecule-calls",
        source_chrom["molecule_calls"]["path"],
        "--output-dir",
        str(output_dir.resolve()),
        "--legacy-gap-bp",
        str(contract["parameters"]["legacy_gap_bp"]),
        "--max-block-size",
        str(contract["parameters"]["max_block_size"]),
        "--max-block-span-bp",
        str(cap),
        "--expected-target-components",
        str(inventory["k_gt8_components"]),
        "--expected-target-sites",
        str(inventory["k_gt8_sites"]),
    ]


def verify_span_child(
    *,
    cap: int,
    chrom: str,
    output_dir: Path,
    command: Sequence[str],
    source_chrom: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    receipt_path = output_dir / "receipt.json"
    receipt_sha = verify_sha256_sidecar_strict(
        receipt_path,
        f"span child receipt {cap} {chrom}",
    )
    receipt = strict_json_load(receipt_path)
    if (
        not isinstance(receipt, Mapping)
        or receipt.get("schema_name") != SPAN_CHILD_SCHEMA
        or receipt.get("all_pass") is not True
    ):
        raise ContractError(
            f"span partition child is not all_pass: {cap} {chrom}"
        )
    scope = receipt.get("scope", {})
    inventory = source_chrom["inventory"]
    expected_scope = {
        "dataset": DATASET,
        "chrom": chrom,
        "site_catalog_sites": inventory["ssnv_sites"],
    }
    for key, expected in expected_scope.items():
        if scope.get(key) != expected:
            raise ContractError(
                f"span child scope mismatch: {cap} {chrom} {key}"
            )
    counts = receipt.get("counts", {})
    for key, expected in {
        "target_components": inventory["k_gt8_components"],
        "target_sites": inventory["k_gt8_sites"],
    }.items():
        if counts.get(key, 0) != expected:
            raise ContractError(
                f"span child count mismatch: {cap} {chrom} {key}"
            )
    parameters = receipt.get("parameters", {})
    if parameters.get("max_block_size") != contract["parameters"][
        "max_block_size"
    ]:
        raise ContractError("span child max_block_size drift")
    if parameters.get("max_block_span_bp") != cap:
        raise ContractError("span child max_block_span_bp drift")
    require_true_checks(receipt, f"span child {cap} {chrom}")
    if receipt.get("command") != list(command):
        raise ContractError(f"span child command drift: {cap} {chrom}")
    expected_inputs = {
        "site_catalog": source_chrom["site_catalog"],
        "molecule_calls": source_chrom["molecule_calls"],
    }
    if receipt.get("inputs") != expected_inputs:
        raise ContractError(
            f"span child cached input identity drift: {cap} {chrom}"
        )
    outputs = verify_output_identities(
        receipt.get("outputs", {}),
        output_dir,
    )
    return {
        "receipt": dict(receipt),
        "receipt_identity": {
            "path": str(receipt_path.resolve()),
            "size_bytes": receipt_path.stat().st_size,
            "sha256": receipt_sha,
        },
        "outputs": outputs,
    }


def source_binding(source_chrom: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "extraction_stage_receipt": source_chrom["stage_receipt"],
        "extraction_child_receipt": source_chrom["child_receipt"],
        "site_catalog": source_chrom["site_catalog"],
        "molecule_calls": source_chrom["molecule_calls"],
    }


def ensure_partition_stage_absent(paths: Mapping[str, Path]) -> None:
    occupied = [
        path
        for path in paths.values()
        if path.exists() or path.is_symlink()
    ]
    if occupied:
        raise ContractError(
            "refusing to overwrite partition stage artifacts: "
            + ", ".join(map(str, occupied))
        )


def validate_stage_receipt(
    *,
    args: argparse.Namespace,
    cap: int,
    chrom: str,
    paths: Mapping[str, Path],
    command: Sequence[str],
    source_chrom: Mapping[str, Any],
    contract: Mapping[str, Any],
    contract_sha256: str,
) -> dict[str, Any]:
    stage_sha = verify_sha256_sidecar_strict(
        paths["stage_receipt"],
        f"cached partition stage receipt {cap} {chrom}",
    )
    stage = strict_json_load(paths["stage_receipt"])
    expected_scalars = {
        "schema_name": f"{SCHEMA_NAME}.stage_receipt",
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "span_cap_bp": cap,
        "stage": "cached_partition",
        "status": "COMPLETED",
        "exit_code": 0,
        "contract_sha256": contract_sha256,
    }
    for key, expected in expected_scalars.items():
        if stage.get(key) != expected:
            raise ContractError(
                f"stage receipt mismatch: {cap} {chrom} {key}"
            )
    if stage.get("command") != list(command):
        raise ContractError(f"stage command drift: {cap} {chrom}")
    if stage.get("tool") != contract["tools"]["partitioner"]:
        raise ContractError(f"stage tool drift: {cap} {chrom}")
    if stage.get("source_binding") != source_binding(source_chrom):
        raise ContractError(f"stage source binding drift: {cap} {chrom}")
    for name in ("stdout", "stderr", "resource"):
        verify_identity_record(
            stage.get("logs", {}).get(name),
            expected_path=paths[name],
            label=f"{cap} {chrom} partition {name}",
        )
    child = verify_span_child(
        cap=cap,
        chrom=chrom,
        output_dir=paths["output"],
        command=command,
        source_chrom=source_chrom,
        contract=contract,
    )
    if stage.get("child_receipt") != child["receipt_identity"]:
        raise ContractError(f"stage child binding drift: {cap} {chrom}")
    return {
        "stage_receipt": dict(stage),
        "stage_receipt_identity": {
            "path": str(paths["stage_receipt"].resolve()),
            "size_bytes": paths["stage_receipt"].stat().st_size,
            "sha256": stage_sha,
        },
        "child": child,
        "resumed": True,
    }


def run_partition_stage(
    *,
    args: argparse.Namespace,
    cap: int,
    chrom: str,
    paths: Mapping[str, Path],
    command: Sequence[str],
    source_chrom: Mapping[str, Any],
    contract: Mapping[str, Any],
    contract_sha256: str,
) -> dict[str, Any]:
    ensure_partition_stage_absent(paths)
    paths["stdout"].parent.mkdir(parents=True, exist_ok=True)
    timed_command = [
        str(TIME_BINARY),
        "-v",
        "-o",
        str(paths["resource"]),
        "--",
        *command,
    ]
    started = time.perf_counter()
    with paths["stdout"].open("xb") as stdout, paths["stderr"].open(
        "xb"
    ) as stderr:
        completed = subprocess.run(
            timed_command,
            stdout=stdout,
            stderr=stderr,
            check=False,
        )
    wall_seconds = time.perf_counter() - started
    if completed.returncode != 0:
        failure = {
            "schema_name": f"{SCHEMA_NAME}.stage_failure",
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": utc_now(),
            "all_pass": False,
            "sample": DATASET,
            "chrom": chrom,
            "span_cap_bp": cap,
            "stage": "cached_partition",
            "exit_code": completed.returncode,
            "wall_seconds": wall_seconds,
            "command": list(command),
            "timed_command": timed_command,
            "contract_sha256": contract_sha256,
            "source_binding": source_binding(source_chrom),
            "logs": {
                name: file_identity(paths[name])
                for name in ("stdout", "stderr", "resource")
                if paths[name].is_file()
            },
        }
        write_json_atomic(paths["failure"], failure)
        write_sha256_sidecar(paths["failure"])
        raise ContractError(
            f"span partition failed: cap={cap} chrom={chrom} "
            f"exit={completed.returncode}"
        )
    try:
        child = verify_span_child(
            cap=cap,
            chrom=chrom,
            output_dir=paths["output"],
            command=command,
            source_chrom=source_chrom,
            contract=contract,
        )
    except Exception as exc:
        failure = {
            "schema_name": f"{SCHEMA_NAME}.stage_failure",
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": utc_now(),
            "all_pass": False,
            "sample": DATASET,
            "chrom": chrom,
            "span_cap_bp": cap,
            "stage": "cached_partition_validation",
            "exit_code": completed.returncode,
            "wall_seconds": wall_seconds,
            "command": list(command),
            "timed_command": timed_command,
            "contract_sha256": contract_sha256,
            "source_binding": source_binding(source_chrom),
            "validation_error": str(exc),
            "logs": {
                name: file_identity(paths[name])
                for name in ("stdout", "stderr", "resource")
                if paths[name].is_file()
            },
        }
        write_json_atomic(paths["failure"], failure)
        write_sha256_sidecar(paths["failure"])
        raise
    stage = {
        "schema_name": f"{SCHEMA_NAME}.stage_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "span_cap_bp": cap,
        "stage": "cached_partition",
        "status": "COMPLETED",
        "exit_code": completed.returncode,
        "wall_seconds": wall_seconds,
        "command": list(command),
        "timed_command": timed_command,
        "scheduling_policy": "inherited_from_outer_process",
        "contract_sha256": contract_sha256,
        "tool": contract["tools"]["partitioner"],
        "source_binding": source_binding(source_chrom),
        "child_receipt": child["receipt_identity"],
        "logs": {
            name: file_identity(paths[name])
            for name in ("stdout", "stderr", "resource")
        },
    }
    write_json_atomic(paths["stage_receipt"], stage)
    stage_sha = sha256_path(paths["stage_receipt"])
    write_sha256_sidecar(paths["stage_receipt"])
    return {
        "stage_receipt": stage,
        "stage_receipt_identity": {
            "path": str(paths["stage_receipt"].resolve()),
            "size_bytes": paths["stage_receipt"].stat().st_size,
            "sha256": stage_sha,
        },
        "child": child,
        "resumed": False,
    }


def get_or_run_partition_stage(
    *,
    args: argparse.Namespace,
    cap: int,
    chrom: str,
    paths: Mapping[str, Path],
    command: Sequence[str],
    source_chrom: Mapping[str, Any],
    contract: Mapping[str, Any],
    contract_sha256: str,
) -> dict[str, Any]:
    any_present = any(
        path.exists() or path.is_symlink() for path in paths.values()
    )
    if not any_present:
        return run_partition_stage(
            args=args,
            cap=cap,
            chrom=chrom,
            paths=paths,
            command=command,
            source_chrom=source_chrom,
            contract=contract,
            contract_sha256=contract_sha256,
        )
    if not args.resume:
        raise ContractError(
            f"partition artifacts exist without --resume: {cap} {chrom}"
        )
    if paths["failure"].exists() or paths["failure"].is_symlink():
        raise ContractError(
            f"prior partition failure must be archived: {paths['failure']}"
        )
    unexpected = [
        name
        for name in ("skip_receipt",)
        if paths[name].exists() or paths[name].is_symlink()
    ]
    if unexpected:
        raise ContractError(
            f"completed target has skip artifacts: {cap} {chrom}"
        )
    return validate_stage_receipt(
        args=args,
        cap=cap,
        chrom=chrom,
        paths=paths,
        command=command,
        source_chrom=source_chrom,
        contract=contract,
        contract_sha256=contract_sha256,
    )


def get_or_write_skip(
    *,
    args: argparse.Namespace,
    cap: int,
    chrom: str,
    paths: Mapping[str, Path],
    source_chrom: Mapping[str, Any],
    contract_sha256: str,
) -> dict[str, Any]:
    forbidden = [
        paths[name]
        for name in ("output", "stdout", "stderr", "stage_receipt", "failure")
        if paths[name].exists() or paths[name].is_symlink()
    ]
    if forbidden:
        raise ContractError(
            "zero-target chromosome has unexpected partition artifacts: "
            + ", ".join(map(str, forbidden))
        )
    resource_existed = paths["resource"].exists() or paths[
        "resource"
    ].is_symlink()
    receipt_existed = paths["skip_receipt"].exists() or paths[
        "skip_receipt"
    ].is_symlink()
    if resource_existed != receipt_existed:
        raise ContractError(
            f"incomplete skip artifacts cannot be resumed: cap={cap} "
            f"chrom={chrom}"
        )
    time_payload = (
        'Command being timed: "SKIP_NO_K_GT8_TARGET"\n'
        "Cached partition command executed: no\n"
        "Elapsed (wall clock) time (seconds): 0.000000\n"
        "Exit status: 0\n"
    ).encode("utf-8")
    if resource_existed:
        if not args.resume:
            raise ContractError(
                f"skip time file exists without --resume: {paths['resource']}"
            )
        if paths["resource"].read_bytes() != time_payload:
            raise ContractError(f"skip time file drift: {paths['resource']}")
    else:
        paths["resource"].parent.mkdir(parents=True, exist_ok=True)
        with paths["resource"].open("xb") as handle:
            handle.write(time_payload)
            handle.flush()
            os.fsync(handle.fileno())

    expected = {
        "schema_name": f"{SCHEMA_NAME}.skip_receipt",
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "span_cap_bp": cap,
        "stage": "cached_partition",
        "status": "SKIP_NO_K_GT8_TARGET",
        "reason": (
            "authenticated source extraction inventory contains zero legacy "
            "positional k>8 targets"
        ),
        "wall_seconds": 0.0,
        "contract_sha256": contract_sha256,
        "source_binding": source_binding(source_chrom),
        "inventory": source_chrom["inventory"],
        "resource": file_identity(paths["resource"]),
    }
    value = {**expected, "created_at_utc": utc_now()}
    observed, identity = establish_json_stable(
        paths["skip_receipt"],
        value,
        resume=args.resume,
        stable_keys=tuple(expected),
    )
    return {
        "stage_receipt": observed,
        "stage_receipt_identity": identity,
        "child": None,
        "resumed": receipt_existed,
    }


def tsv_payload(
    rows: Sequence[Mapping[str, Any]],
    fieldnames: Sequence[str],
) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def child_count(child: Mapping[str, Any] | None, key: str) -> Any:
    if child is None:
        return 0
    return child["receipt"].get("counts", {}).get(key, 0)


def run(args: argparse.Namespace) -> dict[str, Any]:
    validate_args(args)
    source_validation_started = time.perf_counter()
    source = verify_source_full(args)
    source_validation_seconds = time.perf_counter() - source_validation_started
    contract = build_run_contract(args, source)
    contract_path, contract_sha = establish_run_contract(args, contract)

    all_rows = []
    cap_receipts = []
    for cap in args.span_caps:
        cap_rows = []
        cap_chromosomes = []
        for chrom in args.chromosomes:
            source_chrom = source["chromosomes"][chrom]
            inventory = source_chrom["inventory"]
            chrom_root = (
                args.output_root
                / "caps"
                / cap_label(cap)
                / "chromosomes"
                / chrom
            )
            paths = stage_paths(chrom_root)
            if inventory["k_gt8_components"] == 0:
                stage = get_or_write_skip(
                    args=args,
                    cap=cap,
                    chrom=chrom,
                    paths=paths,
                    source_chrom=source_chrom,
                    contract_sha256=contract_sha,
                )
            else:
                command = expected_partition_command(
                    args,
                    cap=cap,
                    chrom=chrom,
                    output_dir=paths["output"],
                    source_chrom=source_chrom,
                    contract=contract,
                )
                stage = get_or_run_partition_stage(
                    args=args,
                    cap=cap,
                    chrom=chrom,
                    paths=paths,
                    command=command,
                    source_chrom=source_chrom,
                    contract=contract,
                    contract_sha256=contract_sha,
                )
            child = stage["child"]
            row = {
                "span_cap_bp": cap,
                "chrom": chrom,
                "status": stage["stage_receipt"]["status"],
                "ssnv_sites": inventory["ssnv_sites"],
                "k_gt8_components": inventory["k_gt8_components"],
                "k_gt8_sites": inventory["k_gt8_sites"],
                "k_gt8_max_k": inventory["k_gt8_max_k"],
                "wall_seconds": stage["stage_receipt"].get(
                    "wall_seconds", 0.0
                ),
                "new_blocks": child_count(child, "new_blocks"),
                "exact_patterns": child_count(child, "exact_patterns"),
                "raw_total_molecule_weight": child_count(
                    child, "raw_total_molecule_weight"
                ),
                "raw_retained_molecule_weight": child_count(
                    child, "raw_retained_molecule_weight"
                ),
                "raw_lost_molecule_weight": child_count(
                    child, "raw_lost_molecule_weight"
                ),
                "unavoidable_patterns": child_count(
                    child, "unavoidable_patterns"
                ),
                "unavoidable_size_patterns": child_count(
                    child, "unavoidable_size_patterns"
                ),
                "unavoidable_span_cap_patterns": child_count(
                    child, "unavoidable_span_cap_patterns"
                ),
                "unavoidable_both_limits_patterns": child_count(
                    child, "unavoidable_both_limits_patterns"
                ),
                "stage_receipt": stage["stage_receipt_identity"]["path"],
                "stage_receipt_sha256": stage[
                    "stage_receipt_identity"
                ]["sha256"],
                "time_v": str(paths["resource"].resolve()),
                "time_v_sha256": sha256_path(paths["resource"]),
            }
            cap_rows.append(row)
            all_rows.append(row)
            cap_chromosomes.append(
                {
                    "chrom": chrom,
                    "inventory": inventory,
                    "status": row["status"],
                    "stage_receipt": stage["stage_receipt_identity"],
                    "child_receipt": (
                        None if child is None else child["receipt_identity"]
                    ),
                }
            )

        cap_root = args.output_root / "caps" / cap_label(cap)
        summary_path = cap_root / "chromosome_summary.tsv"
        summary_fields = tuple(cap_rows[0])
        establish_bytes(
            summary_path,
            tsv_payload(cap_rows, summary_fields),
            resume=args.resume,
        )
        totals = {
            "chromosomes": len(cap_rows),
            "completed_partition_chromosomes": sum(
                row["status"] == "COMPLETED" for row in cap_rows
            ),
            "zero_target_skipped_chromosomes": sum(
                row["status"] == "SKIP_NO_K_GT8_TARGET"
                for row in cap_rows
            ),
            "k_gt8_components": sum(
                int(row["k_gt8_components"]) for row in cap_rows
            ),
            "k_gt8_sites": sum(int(row["k_gt8_sites"]) for row in cap_rows),
            "cached_partition_wall_seconds": sum(
                float(row["wall_seconds"])
                for row in cap_rows
                if row["status"] == "COMPLETED"
            ),
            "new_blocks": sum(int(row["new_blocks"]) for row in cap_rows),
            "exact_patterns": sum(
                int(row["exact_patterns"]) for row in cap_rows
            ),
            "raw_total_molecule_weight": sum(
                int(row["raw_total_molecule_weight"]) for row in cap_rows
            ),
            "raw_retained_molecule_weight": sum(
                int(row["raw_retained_molecule_weight"])
                for row in cap_rows
            ),
            "raw_lost_molecule_weight": sum(
                int(row["raw_lost_molecule_weight"]) for row in cap_rows
            ),
            "unavoidable_patterns": sum(
                int(row["unavoidable_patterns"]) for row in cap_rows
            ),
            "unavoidable_size_patterns": sum(
                int(row["unavoidable_size_patterns"]) for row in cap_rows
            ),
            "unavoidable_span_cap_patterns": sum(
                int(row["unavoidable_span_cap_patterns"])
                for row in cap_rows
            ),
            "unavoidable_both_limits_patterns": sum(
                int(row["unavoidable_both_limits_patterns"])
                for row in cap_rows
            ),
        }
        cap_receipt_value = {
            "schema_name": f"{SCHEMA_NAME}.cap_receipt",
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": utc_now(),
            "all_pass": True,
            "sample": DATASET,
            "span_cap_bp": cap,
            "scope": {
                "chromosomes": list(args.chromosomes),
                "test_mode": args.test_mode,
            },
            "contract": {
                "path": str(contract_path.resolve()),
                "size_bytes": contract_path.stat().st_size,
                "sha256": contract_sha,
            },
            "totals": totals,
            "chromosomes": cap_chromosomes,
            "outputs": {
                "chromosome_summary": file_identity(summary_path),
            },
        }
        cap_receipt_path = cap_root / "receipt.json"
        cap_receipt, cap_receipt_identity = establish_json_stable(
            cap_receipt_path,
            cap_receipt_value,
            resume=args.resume,
            stable_keys=(
                "schema_name",
                "schema_version",
                "all_pass",
                "sample",
                "span_cap_bp",
                "scope",
                "contract",
                "totals",
                "chromosomes",
                "outputs",
            ),
        )
        cap_receipts.append(
            {
                "span_cap_bp": cap,
                "receipt": cap_receipt_identity,
                "totals": cap_receipt["totals"],
            }
        )

    summary_path = args.output_root / "summary.tsv"
    establish_bytes(
        summary_path,
        tsv_payload(all_rows, tuple(all_rows[0])),
        resume=args.resume,
    )
    aggregate_totals = {
        "span_caps": len(args.span_caps),
        "chromosome_cap_tasks": len(all_rows),
        "completed_partition_tasks": sum(
            row["status"] == "COMPLETED" for row in all_rows
        ),
        "zero_target_skipped_tasks": sum(
            row["status"] == "SKIP_NO_K_GT8_TARGET"
            for row in all_rows
        ),
        "cached_partition_wall_seconds": sum(
            float(row["wall_seconds"])
            for row in all_rows
            if row["status"] == "COMPLETED"
        ),
    }
    receipt_value = {
        "schema_name": f"{SCHEMA_NAME}.run_receipt",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "comprehensive_all_pass": (
            not args.test_mode and args.chromosomes == AUTOSOMES
        ),
        "sample": DATASET,
        "scope": {
            "chromosomes": list(args.chromosomes),
            "span_caps_bp": list(args.span_caps),
            "test_mode": args.test_mode,
        },
        "contract": {
            "path": str(contract_path.resolve()),
            "size_bytes": contract_path.stat().st_size,
            "sha256": contract_sha,
        },
        "source_validation_seconds_this_invocation": (
            source_validation_seconds
        ),
        "timing_scope": (
            "cached_partition_wall_seconds excludes source identity "
            "validation and excludes upstream BAM extraction"
        ),
        "totals": aggregate_totals,
        "caps": cap_receipts,
        "outputs": {
            "summary": file_identity(summary_path),
        },
        "claim_ceiling": (
            "engineering sensitivity of deterministic local partitions to a "
            "hard physical span cap; not a unique true evolutionary tree"
        ),
    }
    receipt_path = args.output_root / "receipt.json"
    final_receipt, final_identity = establish_json_stable(
        receipt_path,
        receipt_value,
        resume=args.resume,
        stable_keys=(
            "schema_name",
            "schema_version",
            "all_pass",
            "comprehensive_all_pass",
            "sample",
            "scope",
            "contract",
            "timing_scope",
            "totals",
            "caps",
            "outputs",
            "claim_ceiling",
        ),
    )
    marker_name = "_TEST_SUCCESS" if args.test_mode else "_SUCCESS"
    marker_path = args.output_root / marker_name
    marker_value = {
        "schema_name": f"{SCHEMA_NAME}.success_marker",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": True,
        "comprehensive": not args.test_mode,
        "sample": DATASET,
        "scope": {
            "chromosomes": list(args.chromosomes),
            "span_caps_bp": list(args.span_caps),
        },
        "run_receipt": final_identity,
    }
    marker, marker_identity = establish_json_stable(
        marker_path,
        marker_value,
        resume=args.resume,
        stable_keys=(
            "schema_name",
            "schema_version",
            "all_pass",
            "comprehensive",
            "sample",
            "scope",
            "run_receipt",
        ),
    )
    print(
        json.dumps(
            {
                "all_pass": marker["all_pass"],
                "comprehensive_all_pass": final_receipt[
                    "comprehensive_all_pass"
                ],
                "totals": final_receipt["totals"],
                "receipt": final_identity,
                "marker": marker_identity,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return final_receipt


def main() -> None:
    try:
        args = parse_args()
        run(args)
    except (
        ContractError,
        OSError,
        ValueError,
        subprocess.SubprocessError,
    ) as exc:
        print(f"FAIL-CLOSED: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
