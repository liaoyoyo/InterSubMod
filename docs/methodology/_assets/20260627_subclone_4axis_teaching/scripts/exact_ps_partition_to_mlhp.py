#!/usr/bin/env python3
"""Convert validated exact-PS partition artifacts into layered-tree input.

The adapter is deliberately fail closed.  Every output group represents one
``exact PS x HP x read-linked component x bounded block``.  Only constraints
classified as ``retained`` are projected into that block.  Different phase
sets are never pooled, even when their HP family labels are equal.

Tree patterns are rebuilt from the extractor's authoritative per-molecule
``R/A/O/D/S/L/X`` rows; the segmentation constraints are used only to validate
the block boundary and disposition mass.  Non-R/A observations are marginalized
to X instead of being copied into multiple completions.  The current tree solver
still consumes pattern presence, not a weighted read likelihood, so projected
HP-specific coverage is not reported as caller VAF.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Iterable, Mapping


SCHEMA_NAME = "intersubmod.exact_ps_layered_tree_input"
SCHEMA_VERSION = "1.0.0"
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
RETAINED = "retained"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chroms", default=",".join(AUTOSOMES))
    parser.add_argument("--min-read", type=int, default=3)
    parser.add_argument("--cn-bed", type=Path)
    parser.add_argument("--cn-source", default="unavailable")
    parser.add_argument(
        "--allow-missing-cpp-parity",
        action="store_true",
        help="Research-only escape hatch; production validation must not use it.",
    )
    parser.add_argument(
        "--require-strict-endpoint-receipt",
        action="store_true",
        help=(
            "Production gate: require each chromosome's strict_regions/receipt.json, "
            "verify its exact-PS endpoint contract, and bind its membership SHA-256 "
            "to the partition input."
        ),
    )
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def receipt_size_bytes(spec: Mapping[str, object]) -> int | None:
    """Normalize legacy ``bytes`` and current ``size_bytes`` receipt fields."""

    value = spec.get("size_bytes", spec.get("bytes"))
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        return None
    return value


def read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return value


def read_tsv_gz(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(required) - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{path}: missing columns {','.join(missing)}")
        return [dict(row) for row in reader]


def parse_positions(value: str, *, label: str) -> tuple[int, ...]:
    if not value:
        raise ValueError(f"{label}: empty positions")
    try:
        values = tuple(int(token) for token in value.split(","))
    except ValueError as exc:
        raise ValueError(f"{label}: non-integer position") from exc
    if values != tuple(sorted(set(values))):
        raise ValueError(f"{label}: positions must be strictly increasing and unique")
    return values


def load_cn(path: Path | None) -> dict[str, list[tuple[int, int, str]]] | None:
    if path is None:
        return None
    if not path.is_file():
        raise FileNotFoundError(path)
    result: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected >=4 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid interval") from exc
            result[fields[0]].append((start, end, fields[3].lower()))
    for values in result.values():
        values.sort()
    return result


def cn_state(
    intervals: Mapping[str, list[tuple[int, int, str]]] | None,
    chrom: str,
    pos1: int,
) -> str:
    if intervals is None:
        return "unavailable"
    for start, end, label in intervals.get(chrom, []):
        if start <= pos1 <= end:
            return label
    return "neutral"


def atomic_write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ) + "\n"
    descriptor, temp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".pending", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_name, path)
    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def locate_chrom_dir(root: Path, chrom: str) -> Path:
    candidate = root / "chromosomes" / chrom
    if candidate.is_dir():
        return candidate
    candidate = root / chrom
    if candidate.is_dir():
        return candidate
    raise FileNotFoundError(f"no chromosome directory for {chrom} under {root}")


def build_chrom_groups(
    chrom_dir: Path,
    *,
    sample: str,
    chrom: str,
    min_read: int,
    cn_intervals: Mapping[str, list[tuple[int, int, str]]] | None,
    require_cpp_parity: bool,
    require_strict_endpoint_receipt: bool = False,
) -> tuple[list[dict[str, object]], Counter[str], dict[str, object]]:
    comparison_path = chrom_dir / "comparison.json"
    comparison = read_json(comparison_path)
    if require_cpp_parity and (
        comparison.get("all_pass") is not True
        or int(comparison.get("mismatch_count", 0)) != 0
    ):
        raise ValueError(f"{comparison_path}: Python/C++ parity did not pass")

    partition_dir = chrom_dir / "python_partition"
    receipt_path = partition_dir / "receipt.json"
    receipt = read_json(receipt_path)
    if receipt.get("all_pass") is not True:
        raise ValueError(f"{receipt_path}: partition receipt is not all_pass")
    scope = receipt.get("scope", {})
    if scope.get("dataset") != sample or scope.get("chrom") != chrom:
        raise ValueError(f"{receipt_path}: scope mismatch")
    checks = receipt.get("checks", {})
    for check in ("cross_ps_zero", "cross_hp_zero", "constraint_mass_conserved"):
        if checks.get(check) is not True:
            raise ValueError(f"{receipt_path}: required check failed: {check}")

    strict_receipt_path = chrom_dir / "strict_regions" / "receipt.json"
    strict_receipt = None
    strict_funnel: dict[str, int] = {}
    if require_strict_endpoint_receipt:
        strict_receipt = read_json(strict_receipt_path)
        if (
            strict_receipt.get("schema_name")
            != "intersubmod.strict_exact_ps_hp_regions"
            or strict_receipt.get("schema_version") != "1.1.0"
            or strict_receipt.get("all_pass") is not True
        ):
            raise ValueError(f"{strict_receipt_path}: strict endpoint receipt is not PASS")
        strict_scope = strict_receipt.get("scope", {})
        strict_contract = strict_receipt.get("contract", {})
        strict_checks = strict_receipt.get("checks", {})
        if (
            strict_scope.get("dataset") != sample
            or strict_scope.get("chrom") != chrom
            or strict_contract.get("container")
            != "chromosome x exact nonmissing PS x HP1/HP2"
            or strict_contract.get("region")
            != "connected component of endpoint edges with support >= threshold"
            or strict_contract.get("tree_region")
            != (
                "multisite connected component (k>1); unlinked singleton components "
                "are retained only for funnel conservation and never enter tree inference"
            )
            or strict_checks.get("all_reported_edges_have_direct_fixed_endpoint_support")
            is not True
            or strict_checks.get("every_multisite_component_connected_by_retained_edges")
            is not True
        ):
            raise ValueError(f"{strict_receipt_path}: strict endpoint contract mismatch")
        strict_membership = (strict_receipt.get("outputs") or {}).get("membership")
        partition_membership = (receipt.get("inputs") or {}).get(
            "site_component_membership"
        )
        if (
            not isinstance(strict_membership, dict)
            or not isinstance(partition_membership, dict)
            or strict_membership.get("sha256") != partition_membership.get("sha256")
            or receipt_size_bytes(strict_membership)
            != receipt_size_bytes(partition_membership)
            or receipt_size_bytes(strict_membership) is None
        ):
            raise ValueError(
                f"{strict_receipt_path}: strict membership is not the partition input"
            )
        strict_summary = (strict_receipt.get("summary_by_threshold") or {}).get("3")
        if not isinstance(strict_summary, dict):
            raise ValueError(f"{strict_receipt_path}: threshold-3 summary is missing")
        strict_membership_rows = read_tsv_gz(
            Path(str(strict_membership["path"])),
            ("threshold", "chrom", "pos1", "inference_role", "tree_eligible"),
        )
        threshold_rows = [
            row for row in strict_membership_rows if int(row["threshold"]) == 3
        ]
        active_unique_loci = len(
            {(row["chrom"], int(row["pos1"])) for row in threshold_rows}
        )
        strict_funnel = {
            "strict_candidate_loci_S": int(strict_scope["candidate_loci_S"]),
            "strict_active_unique_loci": active_unique_loci,
            "strict_exact_ps_hp_containers": int(strict_summary["containers"]),
            "strict_active_node_memberships": int(
                strict_summary["active_node_memberships"]
            ),
            "strict_components_all": int(strict_summary["components_all"]),
            "strict_singleton_unlinked_components": int(
                strict_summary["singleton_unlinked_components"]
            ),
            "strict_tree_eligible_regions": int(
                strict_summary["tree_eligible_regions"]
            ),
        }
        if (
            len(threshold_rows) != strict_funnel["strict_active_node_memberships"]
            or strict_funnel["strict_components_all"]
            != strict_funnel["strict_singleton_unlinked_components"]
            + strict_funnel["strict_tree_eligible_regions"]
            or any(
                (row["inference_role"] == "PRIMARY_PS_AWARE")
                != (row["tree_eligible"] == "true")
                for row in threshold_rows
            )
        ):
            raise ValueError(f"{strict_receipt_path}: strict singleton/tree funnel mismatch")
    stage_receipt_path = chrom_dir / "stage_receipt.json"
    if stage_receipt_path.is_file():
        stage_receipt = read_json(stage_receipt_path)
        if stage_receipt.get("all_pass") is not True:
            raise ValueError(f"{stage_receipt_path}: stage receipt is not all_pass")
        stage_metrics = stage_receipt.get("metrics", {})
    else:
        stage_receipt = None
        stage_metrics = {}

    block_path = partition_dir / "blocks.tsv.gz"
    disposition_path = partition_dir / "dispositions.tsv.gz"
    blocks = read_tsv_gz(
        block_path,
        (
            "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
            "hp_family", "phase_set", "block_id", "block_index", "k", "positions1",
        ),
    )
    dispositions = read_tsv_gz(
        disposition_path,
        (
            "dataset", "chrom", "unit_id", "constraint_id", "hp_family",
            "phase_set", "positions1", "call_codes", "molecule_weight",
            "disposition", "retained_block_index",
        ),
    )

    block_by_key: dict[tuple[str, str], dict[str, object]] = {}
    patterns: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    projected_molecule_rows: Counter[tuple[str, str]] = Counter()
    retained_constraint_weight: Counter[tuple[str, str]] = Counter()
    metrics: Counter[str] = Counter()
    metrics.update(strict_funnel)
    for row in blocks:
        if row["dataset"] != sample or row["chrom"] != chrom:
            raise ValueError(f"{block_path}: dataset/chrom mismatch")
        if row["hp_family"] not in {"1", "2"} or not row["phase_set"]:
            raise ValueError(f"{block_path}: non-primary HP or missing PS")
        positions = parse_positions(row["positions1"], label=row["block_id"])
        if len(positions) != int(row["k"]):
            raise ValueError(f"{block_path}: k/positions mismatch for {row['block_id']}")
        key = (row["unit_id"], row["block_index"])
        if key in block_by_key:
            raise ValueError(f"{block_path}: duplicate unit/block key {key}")
        block_by_key[key] = {**row, "positions": positions}
        metrics["blocks_total"] += 1
        metrics[f"blocks_k_{min(len(positions), 13)}"] += 1

    partition_counts = receipt.get("counts", {})
    metrics["candidate_loci_S"] = int(
        strict_funnel.get(
            "strict_candidate_loci_S",
            stage_metrics.get("S", partition_counts.get("site_catalog_rows", 0)),
        )
    )
    metrics["primary_unique_loci"] = int(
        strict_funnel.get(
            "strict_active_unique_loci", stage_metrics.get("unique_sites", 0)
        )
    )
    if metrics["primary_unique_loci"] == 0:
        metrics["primary_unique_loci"] = len(
            {position for block in block_by_key.values() for position in block["positions"]}
        )
    metrics["exact_ps_hp_units"] = int(
        stage_metrics.get("units", partition_counts.get("eligible_units", 0))
    )
    metrics["unit_memberships"] = int(
        stage_metrics.get(
            "unit_memberships", partition_counts.get("eligible_unit_sites", 0)
        )
    )
    metrics["bounded_blocks"] = int(
        stage_metrics.get("blocks", partition_counts.get("blocks", len(blocks)))
    )
    metrics["python_cpp_mismatch_count"] = int(
        stage_metrics.get("python_cpp_mismatch_count", comparison.get("mismatch_count", 0))
    )

    for row in dispositions:
        if row["dataset"] != sample or row["chrom"] != chrom:
            raise ValueError(f"{disposition_path}: dataset/chrom mismatch")
        status = row["disposition"]
        weight = int(row["molecule_weight"])
        metrics["constraint_rows_total"] += 1
        metrics["constraint_weight_total"] += weight
        metrics[f"constraint_rows_{status}"] += 1
        metrics[f"constraint_weight_{status}"] += weight
        if status != RETAINED:
            continue
        block_index = row["retained_block_index"]
        key = (row["unit_id"], block_index)
        block = block_by_key.get(key)
        if block is None:
            raise ValueError(f"{disposition_path}: retained constraint has no block {key}")
        if row["hp_family"] != block["hp_family"] or row["phase_set"] != block["phase_set"]:
            raise ValueError(f"{disposition_path}: retained constraint crosses HP/PS")
        observed = parse_positions(row["positions1"], label=row["constraint_id"])
        codes = row["call_codes"]
        if len(observed) != len(codes) or set(codes) - {"R", "A"}:
            raise ValueError(f"{disposition_path}: invalid R/A constraint")
        block_positions = block["positions"]
        if not set(observed).issubset(block_positions):
            raise ValueError(f"{disposition_path}: constraint leaves retained block")
        retained_constraint_weight[key] += weight

    molecule_candidates = sorted(
        (chrom_dir / "extraction").glob("*.molecule_sparse_calls.tsv.gz")
    )
    if len(molecule_candidates) != 1:
        raise ValueError(
            f"{chrom_dir}: expected exactly one molecule_sparse_calls.tsv.gz, "
            f"found {len(molecule_candidates)}"
        )
    molecule_path = molecule_candidates[0]
    molecule_rows = read_tsv_gz(
        molecule_path,
        (
            "dataset", "chrom", "molecule_id", "hp_family", "phase_set",
            "positions1", "call_codes",
        ),
    )
    route: dict[tuple[str, str, int], tuple[str, str]] = {}
    for key, block in block_by_key.items():
        for position in block["positions"]:
            route_key = (str(block["hp_family"]), str(block["phase_set"]), position)
            previous = route.get(route_key)
            if previous is not None:
                raise ValueError(
                    f"{block_path}: route-site belongs to multiple blocks: {route_key}"
                )
            route[route_key] = key

    seen_molecules: set[str] = set()
    for row in molecule_rows:
        metrics["molecule_rows_total"] += 1
        if row["dataset"] != sample or row["chrom"] != chrom:
            raise ValueError(f"{molecule_path}: dataset/chrom mismatch")
        molecule_id = row["molecule_id"]
        if not molecule_id or molecule_id in seen_molecules:
            raise ValueError(f"{molecule_path}: empty or duplicate molecule_id")
        seen_molecules.add(molecule_id)
        hp = row["hp_family"]
        phase_set = row["phase_set"]
        if hp not in {"1", "2"} or not phase_set:
            metrics["molecule_rows_nonprimary_or_missing_ps"] += 1
            continue
        observed = parse_positions(row["positions1"], label=molecule_id)
        codes = row["call_codes"]
        if len(observed) != len(codes) or set(codes) - {"R", "A", "O", "D", "S", "L", "X"}:
            raise ValueError(f"{molecule_path}: invalid molecule call vector")
        projected: dict[tuple[str, str], dict[int, str]] = defaultdict(dict)
        for position, code in zip(observed, codes):
            key = route.get((hp, phase_set, position))
            if key is not None:
                projected[key][position] = code if code in {"R", "A"} else "X"
        for key, calls in projected.items():
            block_positions = block_by_key[key]["positions"]
            vector = "".join(calls.get(position, "X") for position in block_positions)
            if set(vector) == {"X"}:
                metrics["projected_molecule_blocks_without_fixed_ra"] += 1
                continue
            patterns[key][vector] += 1
            projected_molecule_rows[key] += 1
            metrics["projected_molecule_block_incidences"] += 1

    groups: list[dict[str, object]] = []
    stable_ids: set[str] = set()
    for key, block in sorted(
        block_by_key.items(),
        key=lambda item: (item[1]["chrom"], item[1]["positions"], item[1]["block_id"]),
    ):
        positions = block["positions"]
        if len(positions) < 2:
            metrics["blocks_k1_not_tree_eligible"] += 1
            continue
        supported = Counter(
            {vector: weight for vector, weight in patterns.get(key, {}).items() if weight >= min_read}
        )
        if not supported:
            metrics["blocks_pattern_unsupported"] += 1
            metrics["memberships_pattern_unsupported"] += len(positions)
            continue
        full = Counter({vector: weight for vector, weight in supported.items() if "X" not in vector})
        partial = Counter({vector: weight for vector, weight in supported.items() if "X" in vector})
        hp = str(block["hp_family"])
        phase_set = str(block["phase_set"])
        region_id = f"{chrom}|PS={phase_set}|HP={hp}|{block['block_id']}"
        if region_id in stable_ids:
            raise ValueError(f"duplicate stable region_id: {region_id}")
        stable_ids.add(region_id)
        coverage = {str(position): [0, 0] for position in positions}
        for vector, weight in supported.items():
            for position, code in zip(positions, vector):
                if code == "R":
                    coverage[str(position)][0] += weight
                elif code == "A":
                    coverage[str(position)][1] += weight
        projected_weight = sum(patterns.get(key, {}).values())
        supported_weight = sum(supported.values())
        group = {
            "region_id": region_id,
            "analysis_unit": "exact_ps_hp_component_bounded_block",
            "dataset": sample,
            "chrom": chrom,
            "start": positions[0],
            "end": positions[-1],
            "span": positions[-1] - positions[0],
            "positions": list(positions),
            "n_sSNV": len(positions),
            "unit_id": block["unit_id"],
            "component_id": block["component_id"],
            "block_id": block["block_id"],
            "block_index": int(block["block_index"]),
            "linkage_basis": block["linkage_basis"],
            "hp_family": hp,
            "phase_set": phase_set,
            "phase_set_status": "KNOWN_PS_PRIMARY",
            "populations_by_hp": {hp: dict(sorted(full.items()))} if full else {},
            "subread_groups_by_hp": {hp: dict(sorted(partial.items()))} if partial else {},
            "col_coverage_by_hp": {hp: coverage},
            "reads_by_hp": {hp: supported_weight},
            "n_full_cov_reads": sum(full.values()),
            "retained_segmentation_constraint_weight": retained_constraint_weight.get(key, 0),
            "projected_molecule_block_incidences": projected_weight,
            "tree_supported_molecule_block_incidences": supported_weight,
            "projected_molecule_rows": projected_molecule_rows.get(key, 0),
            "tree_supported_pattern_count": len(supported),
            "min_read": min_read,
            "vaf_eligible": False,
            "coverage_interpretation": (
                "exact-PS projected molecule support after O/D/S/L/X->X marginalization; "
                "HP-specific and not caller VAF"
            ),
            "cn": cn_state(cn_intervals, chrom, positions[len(positions) // 2]),
            "cn_data_available": cn_intervals is not None,
        }
        groups.append(group)
        metrics["groups_tree_input"] += 1
        metrics["memberships_tree_input"] += len(positions)
        metrics["tree_supported_pattern_count"] += len(supported)
        metrics["tree_supported_molecule_block_incidences"] += supported_weight

    inputs = {
        "comparison": identity(comparison_path),
        "partition_receipt": identity(receipt_path),
        "blocks": identity(block_path),
        "dispositions": identity(disposition_path),
        "molecule_sparse_calls": identity(molecule_path),
    }
    if strict_receipt is not None:
        inputs["strict_endpoint_receipt"] = identity(strict_receipt_path)
    if stage_receipt is not None:
        inputs["stage_receipt"] = identity(stage_receipt_path)
    return groups, metrics, inputs


def main() -> int:
    args = parse_args()
    if args.min_read < 1:
        raise ValueError("min-read must be positive")
    chroms = tuple(token for token in args.chroms.split(",") if token)
    if not chroms or len(chroms) != len(set(chroms)):
        raise ValueError("chroms must be a non-empty unique comma-separated list")
    if any(chrom not in AUTOSOMES for chrom in chroms):
        raise ValueError("this adapter currently accepts chr1-chr22 only")
    run_receipt_path = args.partition_root / "run_receipt.json"
    run_receipt = read_json(run_receipt_path)
    if run_receipt.get("all_pass") is not True:
        raise ValueError("partition run_receipt is not all_pass")
    if run_receipt.get("sample") != args.sample:
        raise ValueError("partition run_receipt sample mismatch")

    cn_intervals = load_cn(args.cn_bed)
    groups: list[dict[str, object]] = []
    totals: Counter[str] = Counter()
    chrom_inputs: dict[str, object] = {}
    for chrom in chroms:
        chrom_dir = locate_chrom_dir(args.partition_root, chrom)
        chrom_groups, metrics, inputs = build_chrom_groups(
            chrom_dir,
            sample=args.sample,
            chrom=chrom,
            min_read=args.min_read,
            cn_intervals=cn_intervals,
            require_cpp_parity=not args.allow_missing_cpp_parity,
            require_strict_endpoint_receipt=args.require_strict_endpoint_receipt,
        )
        groups.extend(chrom_groups)
        totals.update(metrics)
        chrom_inputs[chrom] = inputs

    region_ids = [str(group["region_id"]) for group in groups]
    if len(region_ids) != len(set(region_ids)):
        raise AssertionError("stable region IDs are not globally unique")
    if any(not group["phase_set"] or group["hp_family"] not in {"1", "2"} for group in groups):
        raise AssertionError("non-primary exact-PS group escaped validation")

    funnel = {
        "funnel_contract": (
            "strict_endpoint_exact_ps_membership_v2"
            if args.require_strict_endpoint_receipt
            else "exact_ps_membership_v1"
        ),
        "candidate_loci_S": totals["candidate_loci_S"],
        "active_unique_loci": totals["primary_unique_loci"],
        "primary_unique_loci": totals["primary_unique_loci"],
        "strict_exact_ps_hp_containers": totals["strict_exact_ps_hp_containers"],
        "strict_active_node_memberships": totals["strict_active_node_memberships"],
        "strict_components_all": totals["strict_components_all"],
        "singleton_unlinked_components_abstain": totals[
            "strict_singleton_unlinked_components"
        ],
        "tree_eligible_read_linked_regions": totals["strict_tree_eligible_regions"],
        "exact_ps_hp_units": totals["exact_ps_hp_units"],
        "unit_memberships": totals["unit_memberships"],
        "bounded_blocks": totals["bounded_blocks"],
        "tree_input_groups": totals["groups_tree_input"],
        "tree_input_memberships": totals["memberships_tree_input"],
        "k1_blocks_not_tree_eligible": totals["blocks_k1_not_tree_eligible"],
        "pattern_unsupported_blocks": totals["blocks_pattern_unsupported"],
        "pattern_unsupported_memberships": totals["memberships_pattern_unsupported"],
        "constraint_weight_total": totals["constraint_weight_total"],
        "constraint_weight_retained": totals["constraint_weight_retained"],
        "constraint_weight_cut": totals["constraint_weight_cut"],
        "constraint_weight_unavoidable": sum(
            value
            for key, value in totals.items()
            if key.startswith("constraint_weight_unavoidable")
        ),
        "projected_molecule_block_incidences": totals["projected_molecule_block_incidences"],
        "tree_supported_molecule_block_incidences": totals[
            "tree_supported_molecule_block_incidences"
        ],
        "check_constraint_weight_conservation": totals["constraint_weight_total"]
        == totals["constraint_weight_retained"]
        + totals["constraint_weight_cut"]
        + sum(
            value
            for key, value in totals.items()
            if key.startswith("constraint_weight_unavoidable")
        ),
        "check_stable_region_ids_unique": len(region_ids) == len(set(region_ids)),
        "check_cross_ps_zero": True,
        "check_cross_hp_zero": True,
        "check_cpp_parity_required": not args.allow_missing_cpp_parity,
        "check_strict_component_conservation": (
            not args.require_strict_endpoint_receipt
            or totals["strict_components_all"]
            == totals["strict_singleton_unlinked_components"]
            + totals["strict_tree_eligible_regions"]
        ),
        "check_singletons_excluded_before_partition": (
            not args.require_strict_endpoint_receipt
            or totals["exact_ps_hp_units"] == totals["strict_tree_eligible_regions"]
        ),
        "check_strict_membership_funnel_conserved": (
            not args.require_strict_endpoint_receipt
            or totals["strict_active_node_memberships"]
            == totals["unit_memberships"]
            + totals["strict_singleton_unlinked_components"]
        ),
    }
    document = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "sample": args.sample,
        "scope": {"chromosomes": list(chroms), "scope_label": "chr1-22"},
        "params": {
            "MINREAD": args.min_read,
            "MAX_SNV": 12,
            "input_mode": "exact_ps_partition_adapter",
            "phase_set_role": "exact_nonmissing_primary_evidence_boundary",
            "analysis_unit": "exact PS x HP x read-linked component x bounded block",
            "region_linkage_rule": (
                "strict fixed-R/A endpoint-pair support connected component"
                if args.require_strict_endpoint_receipt
                else "upstream exact-PS read-linked component (strict receipt not required)"
            ),
            "strict_endpoint_receipt_required": args.require_strict_endpoint_receipt,
            "pattern_weight_role": (
                "per-molecule exact-PS block projection for MINREAD; solver uses pattern presence"
            ),
            "partial_read_policy": "direct R/A projection; O/D/S/L/X and out-of-block positions marginalize to X",
            "vaf_role": "not_available_from_sparse_partition_constraints",
            "cn_source": args.cn_source if cn_intervals is not None else "unavailable",
        },
        "input_funnel": funnel,
        "read_tag_census": {
            "phase_set_policy": "exact_nonmissing_fail_closed",
            "cross_ps_violations": 0,
            "cross_hp_violations": 0,
            "python_cpp_mismatch_count": totals["python_cpp_mismatch_count"],
        },
        "n_groups_analyzed": len(groups),
        "groups": groups,
        "inputs": {
            "partition_run_receipt": identity(run_receipt_path),
            "chromosomes": chrom_inputs,
            "cn_bed": identity(args.cn_bed) if args.cn_bed else None,
        },
        "claim_ceiling": (
            "Exact-PS regional mutation-state tree input rebuilt from authoritative per-molecule "
            "R/A/O/D/S/L/X rows. Counts support MINREAD but not a likelihood objective; no VAF "
            "ranking or cellular clone claim."
        ),
    }
    atomic_write_json(args.output, document)
    receipt_path = args.output.with_suffix(args.output.suffix + ".receipt.json")
    receipt = {
        "schema_name": f"{SCHEMA_NAME}.receipt",
        "schema_version": SCHEMA_VERSION,
        "all_pass": all(
            funnel[key]
            for key in (
                "check_constraint_weight_conservation",
                "check_stable_region_ids_unique",
                "check_cross_ps_zero",
                "check_cross_hp_zero",
                "check_cpp_parity_required",
                "check_strict_component_conservation",
                "check_singletons_excluded_before_partition",
                "check_strict_membership_funnel_conserved",
            )
        ),
        "output": identity(args.output),
        "counts": dict(sorted(totals.items())),
        "funnel": funnel,
    }
    atomic_write_json(receipt_path, receipt)
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "receipt": str(receipt_path.resolve()),
                "groups": len(groups),
                "cross_ps": 0,
                "cross_hp": 0,
                "all_pass": receipt["all_pass"],
            },
            sort_keys=True,
        )
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
