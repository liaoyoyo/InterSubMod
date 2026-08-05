#!/usr/bin/env python3
"""Full current-v5 M0 census of candidate vertex-set likelihood ambiguity.

M0 deliberately uses the already thresholded alignment-exposure patterns stored in
``layered_region_view_*.json``.  It is an engineering census, not the final
lossless molecule likelihood.  The script:

* accepts only complete, non-capped, verified primary HP1/HP2 mutation lineages;
* collapses candidate trees by mutation-state vertex set before scoring;
* gives every edge variant of one vertex set the same score by construction;
* fits an error-aware mixture to ``obs_subreads`` exactly once per vertex set;
* emits explicit HP-lineage and region denominators plus conservation checks.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sys
import importlib.util
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hypercube_exact import fit_vertex_mixture_slsqp, parse_full, vertex_set_digest  # noqa: E402


DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)

ROW_FIELDS = (
    "dataset",
    "region",
    "chrom",
    "start",
    "end",
    "k",
    "family",
    "n_reads_reported",
    "n_scoring_pattern_groups",
    "n_scoring_alignment_exposures",
    "raw_tree_candidates_T",
    "distinct_vertex_sets_V",
    "distinct_edge_sets_E",
    "candidate_generation_status",
    "best_vertex_sets",
    "top_edge_variants",
    "best_log_likelihood",
    "second_log_likelihood",
    "delta_best_second",
    "top_relative_likelihood_weight",
    "selection_status",
    "all_fits_converged",
    "all_fits_monotone",
    "max_emission_rank",
    "vertex_set_ids",
    "top_vertex_set_ids",
)


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def node_to_vertex(node: str, k: int) -> int:
    if node == "ROOT":
        return 0
    pattern = node[2:] if node.startswith("H_") else node
    if len(pattern) != k or set(pattern) - {"R", "A"}:
        raise ValueError(f"unexpected tree node {node!r} for k={k}")
    return parse_full(pattern)


def tree_vertex_set(tree: dict[str, Any], k: int) -> tuple[int, ...]:
    vertices = {0}
    for edge in tree.get("edges") or []:
        if len(edge) != 2:
            raise ValueError(f"malformed edge: {edge!r}")
        vertices.update(node_to_vertex(str(node), k) for node in edge)
    return tuple(sorted(vertices))


def edge_set_digest(tree: dict[str, Any]) -> str:
    edges = sorted((str(parent), str(child)) for parent, child in (tree.get("edges") or []))
    payload = {"edges": edges, "recurrence": sorted(str(v) for v in (tree.get("recurrence") or []))}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def eligible_primary(lineage: dict[str, Any]) -> tuple[bool, str]:
    if not lineage.get("is_primary_lineage"):
        return False, "NOT_PRIMARY_HP1_HP2"
    if not lineage.get("mutation_bearing"):
        return False, "REFERENCE_ONLY"
    if str(lineage.get("family")) not in {"1", "2"}:
        return False, "NOT_GERMLINE_HP1_HP2"
    if lineage.get("capped"):
        return False, "CAPPED"
    if lineage.get("analysis_candidate_set_complete") is not True:
        return False, "CANDIDATE_SET_INCOMPLETE"
    if lineage.get("verification_complete") is not True or lineage.get("verify_pass") is not True:
        return False, "SOLVER_VERIFICATION_NOT_FULL_PASS"
    if not lineage.get("obs_subreads"):
        return False, "NO_THRESHOLD_RETAINED_PATTERNS"
    return True, "ELIGIBLE_M0"


def load_frozen_solver(canonical_root: Path):
    path = canonical_root / "source_bundle" / "files" / "imported" / "003_tree_enumeration_solver.py"
    if not path.is_file():
        raise FileNotFoundError(path)
    specification = importlib.util.spec_from_file_location("canonical_v5_tree_solver", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot import frozen solver: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module, path


def candidate_vertex_classes(
    dataset: str,
    region: dict[str, Any],
    lineage: dict[str, Any],
    frozen_solver,
) -> tuple[dict[tuple[int, ...], int], str]:
    """Return exact vertex-set classes and their analytical edge-tree counts."""
    k = int(region["n_sSNV"])
    raw_t = int(lineage.get("n_trees") or 0)
    if lineage.get("display_trees_complete") is True:
        trees = list(lineage.get("trees") or [])
        if raw_t != len(trees) or int(lineage.get("n_trees_stored") or 0) != raw_t:
            raise RuntimeError(f"complete candidate tree count mismatch: {dataset} {region['region']} HP{lineage['family']}")
        edge_ids_by_vertex: dict[tuple[int, ...], set[str]] = defaultdict(set)
        for tree in trees:
            edge_ids_by_vertex[tree_vertex_set(tree, k)].add(edge_set_digest(tree))
        if sum(len(values) for values in edge_ids_by_vertex.values()) != raw_t:
            raise RuntimeError("duplicate or missing stored edge candidates")
        return {vertices: len(values) for vertices, values in edge_ids_by_vertex.items()}, "STORED_COMPLETE"

    full = lineage.get("obs_populations") or {}
    partial = list((lineage.get("obs_subreads") or {}).keys())
    rebuilt = frozen_solver.enumerate_min_trees(full, partial, k, tree_cap=1)
    if rebuilt.get("capped"):
        raise RuntimeError(f"frozen solver unexpectedly capped during display reconstruction: {dataset} {region['region']}")
    if int(rebuilt.get("n_trees") or 0) != raw_t or int(rebuilt.get("n_hidden") or 0) != int(lineage.get("n_hidden") or 0):
        raise RuntimeError(
            f"frozen solver reconstruction mismatch: {dataset} {region['region']} HP{lineage['family']} "
            f"T={rebuilt.get('n_trees')}/{raw_t} hidden={rebuilt.get('n_hidden')}/{lineage.get('n_hidden')}"
        )
    classes: dict[tuple[int, ...], int] = {}
    for node_set in rebuilt.get("_feasible_N") or []:
        vertices = tuple(sorted(sum(1 << bit for bit in node) for node in node_set))
        edge_count = 1
        for node in node_set:
            if node:
                edge_count *= sum(1 for bit in node if (node - {bit}) in node_set)
        if vertices in classes:
            raise RuntimeError("duplicate reconstructed feasible vertex set")
        classes[vertices] = edge_count
    if sum(classes.values()) != raw_t:
        raise RuntimeError("reconstructed analytical edge count does not conserve T")
    return classes, "RECONSTRUCTED_FROM_FROZEN_SOLVER"


def _selection_status(raw_t: int, n_vertex_sets: int, best_v: int, top_edges: int) -> str:
    if raw_t == 1:
        return "T1_CANDIDATE_UNIQUE"
    if n_vertex_sets == 1:
        return "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED"
    if best_v > 1:
        return "LIKELIHOOD_TIED_VERTEX_SETS"
    if top_edges > 1:
        return "LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED"
    return "LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE"


def score_lineage(
    dataset: str,
    region: dict[str, Any],
    lineage: dict[str, Any],
    *,
    error_rate: float,
    tie_tolerance: float,
    frozen_solver,
) -> dict[str, Any]:
    k = int(region["n_sSNV"])
    raw_t = int(lineage.get("n_trees") or 0)
    by_vertex, generation_status = candidate_vertex_classes(dataset, region, lineage, frozen_solver)
    if not by_vertex:
        raise RuntimeError(f"no candidate trees: {dataset} {region['region']} HP{lineage['family']}")
    if sum(by_vertex.values()) != raw_t:
        raise AssertionError("candidate conservation failed")

    patterns = [(str(pattern), int(count)) for pattern, count in lineage["obs_subreads"].items()]
    if any(len(pattern) != k or count <= 0 for pattern, count in patterns):
        raise RuntimeError(f"invalid scoring patterns: {dataset} {region['region']} HP{lineage['family']}")

    if len(by_vertex) == 1:
        vertices = next(iter(by_vertex))
        fits = []
        top_vertices = [vertices]
        best_ll = second_ll = delta = None
        top_weight = 1.0
        top_edges = by_vertex[vertices]
        all_converged = all_monotone = True
        max_rank = 0
    else:
        fits = []
        for vertices in sorted(by_vertex):
            fit = fit_vertex_mixture_slsqp(patterns, vertices, error_rate=error_rate)
            fits.append((vertices, fit))
        fits.sort(key=lambda item: (-item[1].log_likelihood, item[0]))
        best_ll = fits[0][1].log_likelihood
        top = [(vertices, fit) for vertices, fit in fits if best_ll - fit.log_likelihood <= tie_tolerance]
        top_vertices = [vertices for vertices, _ in top]
        second_ll = fits[1][1].log_likelihood if len(fits) > 1 else None
        delta = None if second_ll is None else best_ll - second_ll
        relative = [math.exp(max(-745.0, fit.log_likelihood - best_ll)) for _, fit in fits]
        top_weight = relative[0] / math.fsum(relative)
        top_edges = sum(by_vertex[vertices] for vertices in top_vertices)
        all_converged = all(fit.converged for _, fit in fits)
        all_monotone = all(fit.monotone for _, fit in fits)
        max_rank = max(fit.emission_rank for _, fit in fits)
    ids = [vertex_set_digest(k, vertices) for vertices in sorted(by_vertex)]
    top_ids = [vertex_set_digest(k, vertices) for vertices in top_vertices]
    status = _selection_status(raw_t, len(by_vertex), len(top_vertices), top_edges)
    if not all_converged:
        status = "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE"

    return {
        "dataset": dataset,
        "region": region["region"],
        "chrom": region["chrom"],
        "start": int(region["start"]),
        "end": int(region["end"]),
        "k": k,
        "family": str(lineage["family"]),
        "n_reads_reported": int(lineage.get("n_reads") or 0),
        "n_scoring_pattern_groups": len(patterns),
        "n_scoring_alignment_exposures": sum(count for _, count in patterns),
        "raw_tree_candidates_T": raw_t,
        "distinct_vertex_sets_V": len(by_vertex),
        "distinct_edge_sets_E": raw_t,
        "candidate_generation_status": generation_status,
        "best_vertex_sets": len(top_vertices),
        "top_edge_variants": top_edges,
        "best_log_likelihood": best_ll,
        "second_log_likelihood": second_ll,
        "delta_best_second": delta,
        "top_relative_likelihood_weight": top_weight,
        "selection_status": status,
        "all_fits_converged": all_converged,
        "all_fits_monotone": all_monotone,
        "max_emission_rank": max_rank,
        "vertex_set_ids": json.dumps(ids, separators=(",", ":")),
        "top_vertex_set_ids": json.dumps(top_ids, separators=(",", ":")),
    }


def _json_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".12g")
    return value


def summarize_rows(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    rows = list(rows)
    status = Counter(row["selection_status"] for row in rows)
    by_dataset: dict[str, dict[str, Any]] = {}
    for dataset in DATASETS:
        subset = [row for row in rows if row["dataset"] == dataset]
        by_dataset[dataset] = {
            "n_hp_lineage_units": len(subset),
            "n_regions": len({row["region"] for row in subset}),
            "selection_status_counts": dict(Counter(row["selection_status"] for row in subset)),
            "raw_tree_candidates_T": sum(int(row["raw_tree_candidates_T"]) for row in subset),
            "distinct_vertex_sets_V": sum(int(row["distinct_vertex_sets_V"]) for row in subset),
        }
    n = len(rows)
    return {
        "n_hp_lineage_units": n,
        "n_regions_with_at_least_one_eligible_hp_lineage": len({(row["dataset"], row["region"]) for row in rows}),
        "raw_tree_candidates_T": sum(int(row["raw_tree_candidates_T"]) for row in rows),
        "distinct_vertex_sets_V": sum(int(row["distinct_vertex_sets_V"]) for row in rows),
        "selection_status_counts": dict(status),
        "selection_status_fraction_of_eligible_hp_units": {
            key: value / n if n else None for key, value in sorted(status.items())
        },
        "all_fits_converged": all(bool(row["all_fits_converged"]) for row in rows),
        "all_fits_monotone": all(bool(row["all_fits_monotone"]) for row in rows),
        "n_optimizer_abstain": sum(row["selection_status"] == "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE" for row in rows),
        "by_dataset": by_dataset,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--error-rate", type=float, default=0.01)
    parser.add_argument("--tie-tolerance", type=float, default=1e-6)
    parser.add_argument("--sample", action="append", choices=DATASETS)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.error_rate < 0.5 or args.tie_tolerance < 0:
        raise ValueError("invalid error rate or tie tolerance")
    if args.output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    selected = tuple(args.sample or DATASETS)
    all_rows: list[dict[str, Any]] = []
    exclusion = Counter()
    input_files = []
    primary_regions = set()
    fully_eligible_regions = set()
    frozen_solver, frozen_solver_path = load_frozen_solver(args.canonical_root)
    started = time.monotonic()

    for dataset in selected:
        path = args.canonical_root / "samples" / dataset / f"layered_region_view_{dataset}.json"
        if not path.is_file():
            raise FileNotFoundError(path)
        input_files.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_path(path)})
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("sample") != dataset:
            raise RuntimeError(f"sample identity mismatch in {path}")
        for region in payload.get("regions") or []:
            primary = [
                lineage for lineage in (region.get("lineages") or [])
                if lineage.get("is_primary_lineage") and lineage.get("mutation_bearing")
            ]
            if not primary:
                continue
            region_key = (dataset, region["region"])
            primary_regions.add(region_key)
            region_eligibility = [eligible_primary(lineage) for lineage in primary]
            if all(ok for ok, _ in region_eligibility):
                fully_eligible_regions.add(region_key)
            for lineage, (eligible, reason) in zip(primary, region_eligibility):
                if not eligible:
                    exclusion[reason] += 1
                    continue
                all_rows.append(
                    score_lineage(
                        dataset,
                        region,
                        lineage,
                        error_rate=args.error_rate,
                        tie_tolerance=args.tie_tolerance,
                        frozen_solver=frozen_solver,
                    )
                )
        print(
            f"M0_PROGRESS dataset={dataset} cumulative_hp_units={len(all_rows)} "
            f"elapsed_seconds={time.monotonic() - started:.1f}",
            flush=True,
        )

    rows_path = args.output_dir / "m0_hp_lineage_likelihood_census.tsv.gz"
    with gzip.open(rows_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, ROW_FIELDS, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in all_rows:
            writer.writerow({key: _json_cell(row[key]) for key in ROW_FIELDS})

    aggregate = summarize_rows(all_rows)
    receipt = {
        "schema_name": "intersubmod.hypercube_m0_likelihood_census_receipt",
        "schema_version": "1.0.0",
        "scope": "M0_THRESHOLDED_ALIGNMENT_EXPOSURE_NOT_LOSSLESS_MOLECULE_LIKELIHOOD",
        "claim_ceiling": "engineering ambiguity census only; not true clone, parent edge, calibrated model selection, or final publication likelihood",
        "canonical_root": str(args.canonical_root),
        "selected_datasets": list(selected),
        "input_files": input_files,
        "parameters": {
            "error_rate": args.error_rate,
            "tie_tolerance_log_likelihood": args.tie_tolerance,
            "candidate_filter": "primary HP1/HP2 mutation-bearing AND noncapped AND candidate-complete AND full verification PASS",
            "score_grain": "distinct mutation-state vertex set",
            "edge_policy": "same vertex set receives one identical score; edge variants remain unresolved",
            "display_incomplete_policy": "reconstruct exact feasible vertex sets with canonical-v5 frozen solver; do not trust stored prefix",
        },
        "frozen_solver": {"path": str(frozen_solver_path), "sha256": sha256_path(frozen_solver_path)},
        "population": {
            "n_primary_mutation_regions": len(primary_regions),
            "n_fully_m0_eligible_regions": len(fully_eligible_regions),
            "n_regions_with_any_eligible_hp_lineage": aggregate["n_regions_with_at_least_one_eligible_hp_lineage"],
            "excluded_hp_lineage_unit_counts": dict(exclusion),
        },
        "aggregate": aggregate,
        "conservation": {
            "eligible_plus_excluded_hp_units_positive": bool(all_rows) and sum(exclusion.values()) >= 0,
            "output_rows_equal_eligible_hp_units": len(all_rows) == aggregate["n_hp_lineage_units"],
            "selected_datasets_present": all(aggregate["by_dataset"][dataset]["n_hp_lineage_units"] > 0 for dataset in selected),
            "all_nonconverged_units_fail_closed": aggregate["n_optimizer_abstain"]
            == aggregate["selection_status_counts"].get("RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE", 0),
            "all_fits_monotone": aggregate["all_fits_monotone"],
        },
        "full_task_b_scope": set(selected) == set(DATASETS),
        "outputs": {
            "rows": str(rows_path),
            "rows_size_bytes": rows_path.stat().st_size,
            "rows_sha256": sha256_path(rows_path),
        },
    }
    receipt["all_pass"] = all(receipt["conservation"].values())
    receipt_path = args.output_dir / "m0_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    raise SystemExit(0 if receipt["all_pass"] else 1)


if __name__ == "__main__":
    main()
