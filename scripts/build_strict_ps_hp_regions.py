#!/usr/bin/env python3
"""Build production sSNV regions from strict exact-PS endpoint evidence.

Scientific region contract
--------------------------

* container = chromosome x exact non-missing PS x HP1/HP2
* node = an sSNV with at least one fixed R/A call in that container
* edge = at least ``threshold`` distinct canonical molecules with fixed R/A
  calls at both endpoints
* region = an undirected connected component of retained edges (singletons are
  retained)

Genomic distance is reported but never used to create or remove an edge.
Partial calls (X/L/O/D/S) are retained in the molecule artifact but cannot
support an endpoint edge.  The input must contain one identity-deduplicated row
per canonical molecule; duplicate molecule IDs fail closed.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = REPO_ROOT / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

from strict_endpoint_graph import (  # noqa: E402
    EndpointPairAccumulator,
    MoleculeCalls,
    components_for_thresholds,
    strict_graph_digest,
)


SCHEMA_NAME = "intersubmod.strict_exact_ps_hp_regions"
SCHEMA_VERSION = "1.1.0"
FIXED_CODES = frozenset({"R", "A"})
VALID_CODES = FIXED_CODES | frozenset({"X", "L", "O", "D", "S"})
BASIS_BY_HP = {"1": "PS_HP1", "2": "PS_HP2"}
MISSING_PHASE_SET_TOKENS = frozenset({"", ".", "NA", "N/A", "NAN", "NONE", "NULL"})

COMPONENT_FIELDS = (
    "dataset",
    "chrom",
    "linkage_basis",
    "phase_set",
    "phase_set_status",
    "inference_role",
    "component_class",
    "tree_eligible",
    "threshold",
    "component_index",
    "component_id",
    "start1",
    "end1",
    "k",
    "span_bp",
    "max_adjacent_gap_bp",
    "min_internal_bridge_support",
    "max_internal_bridge_support",
    "contains_gap_gt_50kb",
    "solver_route",
    "linkage_rule",
    "retained_edge_count",
    "min_retained_endpoint_support",
    "max_retained_endpoint_support",
    "container_graph_digest",
)
MEMBERSHIP_FIELDS = (
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
    "pos1",
    "component_id",
    "linkage_rule",
)
EDGE_FIELDS = (
    "dataset",
    "chrom",
    "linkage_basis",
    "hp_family",
    "phase_set",
    "site_i_index",
    "site_j_index",
    "pos_i1",
    "pos_j1",
    "gap_bp",
    "support_total",
    "support_RR",
    "support_RA",
    "support_AR",
    "support_AA",
    "passes_primary_threshold",
    "thresholds_passed",
)
CONTAINER_FIELDS = (
    "dataset",
    "chrom",
    "linkage_basis",
    "hp_family",
    "phase_set",
    "molecule_rows",
    "fixed_ra_calls",
    "active_nodes",
    "endpoint_pairs_observed",
    "graph_digest",
    "threshold",
    "retained_edges",
    "regions",
    "k1_regions",
    "k_gt1_regions",
    "tree_eligible_regions",
    "max_k",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def parse_ints(value: str, *, label: str) -> tuple[int, ...]:
    if value == "":
        return ()
    try:
        result = tuple(int(token) for token in value.split(","))
    except ValueError as exc:
        raise ValueError(f"{label}: expected comma-separated integers") from exc
    if any(item < 0 for item in result):
        raise ValueError(f"{label}: negative integer")
    return result


def require_fields(path: Path, fieldnames: Sequence[str] | None, fields: Iterable[str]) -> None:
    missing = sorted(set(fields) - set(fieldnames or ()))
    if missing:
        raise ValueError(f"{path}: missing columns {','.join(missing)}")


def normalize_phase_set(value: str) -> str:
    """Return a trimmed exact-PS token or ``""`` for standard missing values."""

    if not isinstance(value, str):
        raise TypeError("phase_set must be a string")
    normalized = value.strip()
    if normalized.upper() in MISSING_PHASE_SET_TOKENS:
        return ""
    return normalized


def prepare_output_dir(path: Path, *, require_existing_empty: bool) -> Path:
    path = path.expanduser().absolute()
    if path.is_symlink():
        raise ValueError(f"output directory must not be a symlink: {path}")
    if require_existing_empty:
        if not path.is_dir() or next(path.iterdir(), None) is not None:
            raise ValueError(f"required output directory is absent or non-empty: {path}")
    elif path.exists():
        if not path.is_dir() or next(path.iterdir(), None) is not None:
            raise ValueError(f"output directory must be new or empty: {path}")
    else:
        path.mkdir(parents=True, exist_ok=False)
    return path.resolve()


def write_tsv_gz(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json_exclusive(path: Path, value: Any) -> None:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"
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


def load_site_catalog(path: Path, *, chrom: str) -> tuple[dict[int, int], int]:
    by_index: dict[int, int] = {}
    loci: set[tuple[str, int]] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, ("site_index", "chrom", "pos1", "ref", "alt"))
        for line_number, row in enumerate(reader, start=2):
            if row["chrom"] != chrom:
                raise ValueError(f"{path}:{line_number}: chromosome mismatch")
            try:
                site_index = int(row["site_index"])
                pos1 = int(row["pos1"])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid site index/position") from exc
            if site_index < 0 or pos1 <= 0 or site_index in by_index:
                raise ValueError(f"{path}:{line_number}: invalid or duplicate site")
            locus = (chrom, pos1)
            if locus in loci:
                raise ValueError(f"{path}:{line_number}: duplicate physical locus")
            loci.add(locus)
            by_index[site_index] = pos1
    expected = tuple(range(len(by_index)))
    if tuple(sorted(by_index)) != expected:
        raise ValueError("site catalog indices must be contiguous from zero")
    positions = tuple(by_index[index] for index in expected)
    if any(right <= left for left, right in zip(positions, positions[1:])):
        raise ValueError("site catalog positions must increase with site_index")
    return by_index, len(loci)


def load_molecules(
    path: Path,
    *,
    dataset: str,
    chrom: str,
    site_positions: Mapping[int, int],
) -> tuple[dict[tuple[str, str], EndpointPairAccumulator], Counter[str]]:
    containers: dict[tuple[str, str], EndpointPairAccumulator] = {}
    counts: Counter[str] = Counter()
    seen_molecule_ids: set[str] = set()
    required = (
        "dataset",
        "chrom",
        "molecule_id",
        "hp_family",
        "phase_set",
        "site_indices",
        "positions1",
        "call_codes",
    )
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, required)
        for line_number, row in enumerate(reader, start=2):
            counts["molecule_rows_total"] += 1
            if row["dataset"] != dataset or row["chrom"] != chrom:
                raise ValueError(f"{path}:{line_number}: dataset/chrom mismatch")
            molecule_id = row["molecule_id"]
            if not molecule_id or molecule_id in seen_molecule_ids:
                raise ValueError(f"{path}:{line_number}: empty or duplicate molecule_id")
            seen_molecule_ids.add(molecule_id)

            indices = parse_ints(row["site_indices"], label=f"line {line_number} site_indices")
            positions = parse_ints(row["positions1"], label=f"line {line_number} positions1")
            codes = row["call_codes"]
            if not (len(indices) == len(positions) == len(codes)):
                raise ValueError(f"{path}:{line_number}: call vector length mismatch")
            if len(indices) != len(set(indices)) or indices != tuple(sorted(indices)):
                raise ValueError(f"{path}:{line_number}: site indices must be sorted and unique")
            unknown = sorted(set(codes) - VALID_CODES)
            if unknown:
                raise ValueError(f"{path}:{line_number}: unknown call codes {unknown}")
            for site_index, pos1 in zip(indices, positions):
                if site_positions.get(site_index) != pos1:
                    raise ValueError(f"{path}:{line_number}: site catalog position mismatch")

            hp = row["hp_family"]
            phase_set = normalize_phase_set(row["phase_set"])
            if hp not in BASIS_BY_HP:
                counts["excluded_nonprimary_hp_rows"] += 1
                continue
            if not phase_set:
                counts["excluded_missing_ps_rows"] += 1
                continue
            fixed = tuple(
                (site_index, code)
                for site_index, code in zip(indices, codes)
                if code in FIXED_CODES
            )
            counts["primary_known_ps_rows"] += 1
            counts["primary_known_ps_fixed_ra_calls"] += len(fixed)
            counts["primary_known_ps_nonlinking_calls"] += len(codes) - len(fixed)
            key = (hp, phase_set)
            accumulator = containers.setdefault(key, EndpointPairAccumulator())
            accumulator.add(MoleculeCalls(molecule_id=molecule_id, calls=fixed))
    counts["unique_molecule_ids"] = len(seen_molecule_ids)
    counts["exact_ps_hp_containers"] = len(containers)
    return containers, counts


def phase_token(phase_set: str) -> str:
    return hashlib.sha256(phase_set.encode("utf-8")).hexdigest()[:12]


def build_regions(
    *,
    dataset: str,
    chrom: str,
    site_catalog: Path,
    molecule_calls: Path,
    output_dir: Path,
    thresholds: Sequence[int],
    primary_threshold: int,
    require_existing_empty_output_dir: bool = False,
) -> dict[str, Any]:
    normalized_thresholds = tuple(sorted(set(thresholds)))
    if not normalized_thresholds or any(value <= 0 for value in normalized_thresholds):
        raise ValueError("thresholds must contain positive integers")
    if primary_threshold not in normalized_thresholds:
        raise ValueError("primary threshold must be included in thresholds")
    site_catalog = site_catalog.resolve(strict=True)
    molecule_calls = molecule_calls.resolve(strict=True)
    output_dir = prepare_output_dir(
        output_dir, require_existing_empty=require_existing_empty_output_dir
    )
    site_positions, site_count = load_site_catalog(site_catalog, chrom=chrom)
    containers, counts = load_molecules(
        molecule_calls,
        dataset=dataset,
        chrom=chrom,
        site_positions=site_positions,
    )

    component_rows: list[dict[str, object]] = []
    membership_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    container_rows: list[dict[str, object]] = []
    summary: dict[str, Counter[str]] = {
        str(threshold): Counter() for threshold in normalized_thresholds
    }
    component_owner: dict[tuple[str, str, int, int], str] = {}

    for (hp, phase_set), accumulator in sorted(containers.items()):
        basis = BASIS_BY_HP[hp]
        sites = accumulator.site_indices
        support = accumulator.pair_support()
        digest = strict_graph_digest(sites, support)
        support_by_pair = {(row.site_i, row.site_j): row for row in support}
        partitions = (
            components_for_thresholds(sites, support, thresholds=normalized_thresholds)
            if sites
            else ()
        )
        for row in support:
            passed = tuple(
                threshold for threshold in normalized_thresholds if row.total >= threshold
            )
            edge_rows.append(
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "linkage_basis": basis,
                    "hp_family": hp,
                    "phase_set": phase_set,
                    "site_i_index": row.site_i,
                    "site_j_index": row.site_j,
                    "pos_i1": site_positions[row.site_i],
                    "pos_j1": site_positions[row.site_j],
                    "gap_bp": site_positions[row.site_j] - site_positions[row.site_i],
                    "support_total": row.total,
                    "support_RR": row.rr,
                    "support_RA": row.ra,
                    "support_AR": row.ar,
                    "support_AA": row.aa,
                    "passes_primary_threshold": str(row.total >= primary_threshold).lower(),
                    "thresholds_passed": ",".join(map(str, passed)),
                }
            )
        counts["endpoint_pair_rows_observed"] += len(support)
        counts["endpoint_pair_state_mass"] += sum(row.total for row in support)

        for partition in partitions:
            threshold = partition.threshold
            key = str(threshold)
            retained_edges = sum(row.total >= threshold for row in support)
            k_values = [len(component) for component in partition.components]
            container_rows.append(
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "linkage_basis": basis,
                    "hp_family": hp,
                    "phase_set": phase_set,
                    "molecule_rows": accumulator.molecule_count,
                    "fixed_ra_calls": accumulator.fixed_call_count,
                    "active_nodes": len(sites),
                    "endpoint_pairs_observed": len(support),
                    "graph_digest": digest,
                    "threshold": threshold,
                    "retained_edges": retained_edges,
                    "regions": len(partition.components),
                    "k1_regions": sum(k == 1 for k in k_values),
                    "k_gt1_regions": sum(k > 1 for k in k_values),
                    "tree_eligible_regions": sum(k > 1 for k in k_values),
                    "max_k": max(k_values, default=0),
                }
            )
            summary[key]["containers"] += 1
            summary[key]["active_node_memberships"] += len(sites)
            summary[key]["retained_endpoint_edges"] += retained_edges
            summary[key]["regions"] += len(partition.components)
            summary[key]["k1_regions"] += sum(k == 1 for k in k_values)
            summary[key]["k_gt1_regions"] += sum(k > 1 for k in k_values)
            summary[key]["components_all"] += len(partition.components)
            summary[key]["singleton_unlinked_components"] += sum(
                k == 1 for k in k_values
            )
            summary[key]["tree_eligible_regions"] += sum(k > 1 for k in k_values)

            for component_index, component in enumerate(partition.components, start=1):
                positions = tuple(site_positions[index] for index in component)
                component_set = set(component)
                induced_edges = [
                    row
                    for row in support
                    if row.total >= threshold
                    and row.site_i in component_set
                    and row.site_j in component_set
                ]
                gaps = [right - left for left, right in zip(positions, positions[1:])]
                token_payload = {
                    "dataset": dataset,
                    "chrom": chrom,
                    "hp": hp,
                    "phase_set": phase_set,
                    "threshold": threshold,
                    "sites": component,
                }
                token = hashlib.sha256(
                    json.dumps(token_payload, sort_keys=True, separators=(",", ":")).encode()
                ).hexdigest()[:16]
                component_id = (
                    f"{dataset}:{chrom}:{basis}:PS{phase_token(phase_set)}:"
                    f"E{threshold}:{component_index}:{positions[0]}-{positions[-1]}:{token}"
                )
                tree_eligible = len(component) > 1
                component_class = (
                    "READ_LINKED_MULTISITE_REGION"
                    if tree_eligible
                    else "UNLINKED_SINGLETON_COMPONENT"
                )
                base = {
                    "dataset": dataset,
                    "chrom": chrom,
                    "linkage_basis": basis,
                    "phase_set": phase_set,
                    "phase_set_status": (
                        "KNOWN_PS_PRIMARY"
                        if tree_eligible
                        else "KNOWN_PS_SINGLETON_ABSTAIN"
                    ),
                    "inference_role": (
                        "PRIMARY_PS_AWARE"
                        if tree_eligible
                        else "ABSTAIN_SINGLETON_UNLINKED"
                    ),
                    "component_class": component_class,
                    "tree_eligible": str(tree_eligible).lower(),
                    "threshold": threshold,
                    "component_id": component_id,
                    "linkage_rule": "strict_fixed_ra_endpoint_pair",
                }
                edge_weights = [row.total for row in induced_edges]
                component_rows.append(
                    {
                        **base,
                        "component_index": component_index,
                        "start1": positions[0],
                        "end1": positions[-1],
                        "k": len(component),
                        "span_bp": positions[-1] - positions[0],
                        "max_adjacent_gap_bp": max(gaps, default=0),
                        "min_internal_bridge_support": min(edge_weights, default=""),
                        "max_internal_bridge_support": max(edge_weights, default=""),
                        "contains_gap_gt_50kb": str(any(gap > 50_000 for gap in gaps)).lower(),
                        "solver_route": (
                            "DEFER_TO_K12_PARTITION"
                            if tree_eligible
                            else "EXCLUDE_SINGLETON_NO_READ_LINKAGE"
                        ),
                        "retained_edge_count": len(induced_edges),
                        "min_retained_endpoint_support": min(edge_weights, default=""),
                        "max_retained_endpoint_support": max(edge_weights, default=""),
                        "container_graph_digest": digest,
                    }
                )
                for site_index in component:
                    owner_key = (hp, phase_set, threshold, site_index)
                    if owner_key in component_owner:
                        raise AssertionError(f"site has duplicate component ownership: {owner_key}")
                    component_owner[owner_key] = component_id
                    membership_rows.append(
                        {
                            **base,
                            "site_index": site_index,
                            "pos1": site_positions[site_index],
                        }
                    )
                summary[key][f"k={len(component)}"] += 1
                if len(component) > 1:
                    observed = components_for_thresholds(
                        component, induced_edges, thresholds=(threshold,)
                    )[0].components
                    if observed != (tuple(component),):
                        raise AssertionError(f"component is disconnected: {component_id}")

    expected_memberships = sum(
        len(accumulator.site_indices) * len(normalized_thresholds)
        for accumulator in containers.values()
    )
    state_conserved = all(
        int(row["support_total"])
        == int(row["support_RR"])
        + int(row["support_RA"])
        + int(row["support_AR"])
        + int(row["support_AA"])
        for row in edge_rows
    )
    all_edges_direct = all(
        int(row["support_total"]) > 0 for row in edge_rows
    )
    outputs = {
        "components": output_dir / f"{dataset}.{chrom}.components.tsv.gz",
        "membership": output_dir / f"{dataset}.{chrom}.site_component_membership.tsv.gz",
        "edges": output_dir / f"{dataset}.{chrom}.endpoint_edges.tsv.gz",
        "containers": output_dir / f"{dataset}.{chrom}.container_summary.tsv.gz",
    }
    write_tsv_gz(outputs["components"], COMPONENT_FIELDS, component_rows)
    write_tsv_gz(outputs["membership"], MEMBERSHIP_FIELDS, membership_rows)
    write_tsv_gz(outputs["edges"], EDGE_FIELDS, edge_rows)
    write_tsv_gz(outputs["containers"], CONTAINER_FIELDS, container_rows)

    checks = {
        "unique_molecule_ids": counts["unique_molecule_ids"] == counts["molecule_rows_total"],
        "only_exact_known_ps_hp_primary": all(
            row["phase_set"] and row["linkage_basis"] in {"PS_HP1", "PS_HP2"}
            for row in component_rows
        ),
        "component_membership_mass_conserved": len(membership_rows) == expected_memberships,
        "every_active_node_has_one_component_per_threshold": len(component_owner)
        == expected_memberships,
        "edge_state_mass_conserved": state_conserved,
        "all_reported_edges_have_direct_fixed_endpoint_support": all_edges_direct,
        "every_multisite_component_connected_by_retained_edges": True,
        "distance_not_used_for_connectivity": True,
        "primary_threshold_present": primary_threshold in normalized_thresholds,
        "singleton_components_never_primary": all(
            row["inference_role"] == "ABSTAIN_SINGLETON_UNLINKED"
            and row["tree_eligible"] == "false"
            for row in component_rows
            if int(row["k"]) == 1
        ),
        "primary_components_are_multisite_and_read_linked": all(
            int(row["k"]) > 1
            and row["tree_eligible"] == "true"
            and int(row["retained_edge_count"]) >= 1
            for row in component_rows
            if row["inference_role"] == "PRIMARY_PS_AWARE"
        ),
    }
    receipt_path = output_dir / "receipt.json"
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "all_pass": all(checks.values()),
        "scope": {"dataset": dataset, "chrom": chrom, "candidate_loci_S": site_count},
        "contract": {
            "container": "chromosome x exact nonmissing PS x HP1/HP2",
            "node": "sSNV with >=1 fixed R/A molecule call in the exact container",
            "edge": "distinct canonical molecule fixed R/A at both endpoints",
            "region": "connected component of endpoint edges with support >= threshold",
            "tree_region": (
                "multisite connected component (k>1); unlinked singleton components "
                "are retained only for funnel conservation and never enter tree inference"
            ),
            "partial_calls": "X/L/O/D/S never establish an edge",
            "distance": "diagnostic only; never a boundary or edge",
        },
        "parameters": {
            "thresholds": list(normalized_thresholds),
            "primary_threshold": primary_threshold,
            "linkage_rule": "strict_fixed_ra_endpoint_pair",
        },
        "inputs": {
            "site_catalog": file_identity(site_catalog),
            "molecule_calls": file_identity(molecule_calls),
            "builder": file_identity(Path(__file__).resolve()),
            "strict_graph_core": file_identity(TOOLS_DIR / "strict_endpoint_graph.py"),
        },
        "counts": dict(sorted(counts.items())),
        "summary_by_threshold": {
            key: dict(sorted(values.items())) for key, values in summary.items()
        },
        "checks": checks,
        "outputs": {key: file_identity(path) for key, path in outputs.items()},
        "claim_ceiling": (
            "The receipt proves strict read-linked exact-PS/HP regional connectivity. "
            "It does not prove a unique cellular clone tree or correct CN/VAF interpretation."
        ),
    }
    write_json_exclusive(receipt_path, receipt)
    write_sha256_sidecar(receipt_path)
    return receipt


def parse_thresholds(value: str) -> tuple[int, ...]:
    try:
        values = tuple(int(token) for token in value.split(",") if token)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("thresholds must be comma-separated integers") from exc
    if not values or len(values) != len(set(values)) or min(values) < 1:
        raise argparse.ArgumentTypeError("thresholds must be unique positive integers")
    return values


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--site-catalog", required=True, type=Path)
    parser.add_argument("--molecule-calls", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--thresholds", type=parse_thresholds, default=(1, 2, 3, 5))
    parser.add_argument("--primary-threshold", type=int, default=3)
    parser.add_argument("--require-existing-empty-output-dir", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    receipt = build_regions(
        dataset=args.dataset,
        chrom=args.chrom,
        site_catalog=args.site_catalog,
        molecule_calls=args.molecule_calls,
        output_dir=args.output_dir,
        thresholds=args.thresholds,
        primary_threshold=args.primary_threshold,
        require_existing_empty_output_dir=args.require_existing_empty_output_dir,
    )
    print(
        json.dumps(
            {
                "all_pass": receipt["all_pass"],
                "scope": receipt["scope"],
                "summary_by_threshold": receipt["summary_by_threshold"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
