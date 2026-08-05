#!/usr/bin/env python3
"""Re-fit only M0 optimizer-abstain units with a global KKT certificate.

This is a bounded numerical-method audit.  It never reads BAM and never adds
VAF.  The likelihood and candidate vertex sets are identical to M0; only the
optimizer stopping contract changes from ``scipy.success`` to a concavity-based
global log-likelihood gap certificate.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import scipy

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_m0_likelihood_census as m0  # noqa: E402
from hypercube_exact import fit_vertex_mixture_slsqp, vertex_set_digest  # noqa: E402


UNIT_FIELDS = (
    "dataset", "region", "chrom", "family", "k", "raw_tree_candidates_T",
    "distinct_vertex_sets_V", "n_candidate_fits", "n_slsqp_success",
    "n_slsqp_status_failure", "slsqp_status_counts_json", "n_robust_certified",
    "n_rank_deficient_weight_models", "max_warm_start_global_ll_gap_bound",
    "max_final_global_ll_gap_bound", "max_abs_ll_refinement_change",
    "new_best_log_likelihood", "new_second_log_likelihood", "new_delta_best_second",
    "new_best_vertex_sets", "new_top_vertex_set_ids_json", "old_top_vertex_set_ids_json",
    "top_vertex_set_ids_changed", "counterfactual_selection_status",
    "all_candidates_certified", "all_simplex_residuals_pass", "all_refinements_monotone",
)

CANDIDATE_FIELDS = (
    "dataset", "region", "family", "k", "vertex_set_id", "n_vertices",
    "emission_rank", "augmented_emission_rank", "mixture_weights_identifiable",
    "slsqp_success", "slsqp_status", "slsqp_message", "slsqp_iterations",
    "warm_start_log_likelihood", "warm_start_global_ll_gap_bound",
    "refinement_iterations", "final_log_likelihood", "ll_refinement_change",
    "final_global_ll_gap_bound", "simplex_residual", "optimizer_status",
    "robust_certified", "monotone",
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


def read_abstain_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["selection_status"] == "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE"
        ]
    if len(rows) != len({(row["dataset"], row["region"], row["family"]) for row in rows}):
        raise RuntimeError("duplicate M0 optimizer-abstain unit keys")
    return rows


def json_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".12g")
    return value


def quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"min": None, "p50": None, "p90": None, "p99": None, "max": None}
    array = np.asarray(values, dtype=float)
    return {
        "min": float(np.min(array)),
        "p50": float(np.quantile(array, 0.50)),
        "p90": float(np.quantile(array, 0.90)),
        "p99": float(np.quantile(array, 0.99)),
        "max": float(np.max(array)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-root", required=True, type=Path)
    parser.add_argument("--m0-output-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--error-rate", type=float, default=0.01)
    parser.add_argument("--tie-tolerance", type=float, default=1e-6)
    parser.add_argument("--certificate-tolerance", type=float, default=1e-8)
    parser.add_argument("--sample", action="append", choices=m0.DATASETS)
    parser.add_argument("--max-units", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {args.output_dir}")
    if not 0 < args.error_rate < 0.5 or args.tie_tolerance < 0 or args.certificate_tolerance <= 0:
        raise ValueError("invalid optimizer audit parameters")
    m0_rows_path = args.m0_output_dir / "m0_hp_lineage_likelihood_census.tsv.gz"
    m0_receipt_path = args.m0_output_dir / "m0_receipt.json"
    m0_receipt = json.loads(m0_receipt_path.read_text(encoding="utf-8"))
    all_abstain = read_abstain_rows(m0_rows_path)
    selected_samples = set(args.sample or m0.DATASETS)
    selected = [row for row in all_abstain if row["dataset"] in selected_samples]
    if args.max_units is not None:
        if args.max_units <= 0:
            raise ValueError("max-units must be positive")
        selected = selected[: args.max_units]
    targets = {(row["dataset"], row["region"], row["family"]): row for row in selected}
    args.output_dir.mkdir(parents=True, exist_ok=False)
    frozen_solver, frozen_solver_path = m0.load_frozen_solver(args.canonical_root)
    started = time.monotonic()
    unit_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    found: set[tuple[str, str, str]] = set()

    for dataset in m0.DATASETS:
        dataset_targets = {key for key in targets if key[0] == dataset}
        if not dataset_targets:
            continue
        input_path = args.canonical_root / "samples" / dataset / f"layered_region_view_{dataset}.json"
        payload = json.loads(input_path.read_text(encoding="utf-8"))
        for region in payload.get("regions") or []:
            region_keys = {key for key in dataset_targets if key[1] == region["region"]}
            if not region_keys:
                continue
            for lineage in region.get("lineages") or []:
                key = (dataset, region["region"], str(lineage.get("family")))
                if key not in targets:
                    continue
                found.add(key)
                old = targets[key]
                by_vertex, _ = m0.candidate_vertex_classes(
                    dataset, region, lineage, frozen_solver
                )
                patterns = [
                    (str(pattern), int(count))
                    for pattern, count in lineage["obs_subreads"].items()
                ]
                fits = []
                local_candidates = []
                for vertices in sorted(by_vertex):
                    fit = fit_vertex_mixture_slsqp(
                        patterns,
                        vertices,
                        error_rate=args.error_rate,
                    )
                    vertex_id = vertex_set_digest(int(region["n_sSNV"]), vertices)
                    fits.append((vertices, fit, vertex_id))
                    row = {
                        "dataset": dataset,
                        "region": region["region"],
                        "family": str(lineage["family"]),
                        "k": int(region["n_sSNV"]),
                        "vertex_set_id": vertex_id,
                        "n_vertices": len(vertices),
                        "emission_rank": fit.emission_rank,
                        "augmented_emission_rank": fit.augmented_emission_rank,
                        "mixture_weights_identifiable": fit.mixture_weights_identifiable,
                        "slsqp_success": fit.slsqp_success,
                        "slsqp_status": fit.slsqp_status,
                        "slsqp_message": fit.slsqp_message,
                        "slsqp_iterations": fit.iterations - fit.refinement_iterations,
                        "warm_start_log_likelihood": fit.warm_start_log_likelihood,
                        "warm_start_global_ll_gap_bound": fit.warm_start_global_log_likelihood_gap_bound,
                        "refinement_iterations": fit.refinement_iterations,
                        "final_log_likelihood": fit.log_likelihood,
                        "ll_refinement_change": fit.log_likelihood - fit.warm_start_log_likelihood,
                        "final_global_ll_gap_bound": fit.global_log_likelihood_gap_bound,
                        "simplex_residual": fit.simplex_residual,
                        "optimizer_status": fit.optimizer_status,
                        "robust_certified": fit.converged,
                        "monotone": fit.monotone,
                    }
                    candidate_rows.append(row)
                    local_candidates.append(row)

                fits.sort(key=lambda item: (-item[1].log_likelihood, item[0]))
                best_ll = fits[0][1].log_likelihood
                top = [item for item in fits if best_ll - item[1].log_likelihood <= args.tie_tolerance]
                top_vertices = [item[0] for item in top]
                top_ids = [item[2] for item in top]
                second_ll = fits[1][1].log_likelihood if len(fits) > 1 else None
                top_edges = sum(by_vertex[vertices] for vertices in top_vertices)
                counterfactual = m0._selection_status(
                    int(old["raw_tree_candidates_T"]),
                    int(old["distinct_vertex_sets_V"]),
                    len(top_vertices),
                    top_edges,
                )
                old_top_ids = json.loads(old["top_vertex_set_ids"])
                status_counts = Counter(str(row["slsqp_status"]) for row in local_candidates)
                unit_rows.append({
                    "dataset": dataset,
                    "region": region["region"],
                    "chrom": region["chrom"],
                    "family": str(lineage["family"]),
                    "k": int(region["n_sSNV"]),
                    "raw_tree_candidates_T": int(old["raw_tree_candidates_T"]),
                    "distinct_vertex_sets_V": int(old["distinct_vertex_sets_V"]),
                    "n_candidate_fits": len(local_candidates),
                    "n_slsqp_success": sum(bool(row["slsqp_success"]) for row in local_candidates),
                    "n_slsqp_status_failure": sum(not bool(row["slsqp_success"]) for row in local_candidates),
                    "slsqp_status_counts_json": json.dumps(status_counts, sort_keys=True, separators=(",", ":")),
                    "n_robust_certified": sum(bool(row["robust_certified"]) for row in local_candidates),
                    "n_rank_deficient_weight_models": sum(
                        not bool(row["mixture_weights_identifiable"]) for row in local_candidates
                    ),
                    "max_warm_start_global_ll_gap_bound": max(
                        float(row["warm_start_global_ll_gap_bound"]) for row in local_candidates
                    ),
                    "max_final_global_ll_gap_bound": max(
                        float(row["final_global_ll_gap_bound"]) for row in local_candidates
                    ),
                    "max_abs_ll_refinement_change": max(
                        abs(float(row["ll_refinement_change"])) for row in local_candidates
                    ),
                    "new_best_log_likelihood": best_ll,
                    "new_second_log_likelihood": second_ll,
                    "new_delta_best_second": None if second_ll is None else best_ll - second_ll,
                    "new_best_vertex_sets": len(top),
                    "new_top_vertex_set_ids_json": json.dumps(top_ids, separators=(",", ":")),
                    "old_top_vertex_set_ids_json": json.dumps(old_top_ids, separators=(",", ":")),
                    "top_vertex_set_ids_changed": set(top_ids) != set(old_top_ids),
                    "counterfactual_selection_status": counterfactual,
                    "all_candidates_certified": all(bool(row["robust_certified"]) for row in local_candidates),
                    "all_simplex_residuals_pass": all(
                        float(row["simplex_residual"]) <= 1e-12 for row in local_candidates
                    ),
                    "all_refinements_monotone": all(bool(row["monotone"]) for row in local_candidates),
                })
        print(
            f"OPTIMIZER_AUDIT_PROGRESS dataset={dataset} units={len(unit_rows)}/{len(selected)} "
            f"candidate_fits={len(candidate_rows)} elapsed_seconds={time.monotonic()-started:.1f}",
            flush=True,
        )

    if found != set(targets):
        raise RuntimeError(f"target lookup mismatch missing={len(set(targets)-found)} extra={len(found-set(targets))}")
    unit_path = args.output_dir / "m0_optimizer_nonconvergence_units.tsv.gz"
    candidate_path = args.output_dir / "m0_optimizer_candidate_diagnostics.tsv.gz"
    for path, fields, rows in (
        (unit_path, UNIT_FIELDS, unit_rows),
        (candidate_path, CANDIDATE_FIELDS, candidate_rows),
    ):
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="raise")
            writer.writeheader()
            for row in rows:
                writer.writerow({field: json_cell(row[field]) for field in fields})

    counterfactual = Counter(row["counterfactual_selection_status"] for row in unit_rows)
    slsqp_statuses = Counter(str(row["slsqp_status"]) for row in candidate_rows)
    full_scope = (
        args.max_units is None
        and selected_samples == set(m0.DATASETS)
        and len(selected) == int(m0_receipt["aggregate"]["n_optimizer_abstain"])
    )
    checks = {
        "selected_units_equal_written_units": len(selected) == len(unit_rows),
        "candidate_count_conserves_distinct_vertex_sets": len(candidate_rows)
        == sum(int(row["distinct_vertex_sets_V"]) for row in unit_rows),
        "all_candidates_globally_kkt_certified": all(
            bool(row["robust_certified"]) for row in candidate_rows
        ),
        "all_candidate_simplex_residuals_pass": all(
            float(row["simplex_residual"]) <= 1e-12 for row in candidate_rows
        ),
        "all_pairwise_refinements_monotone": all(bool(row["monotone"]) for row in candidate_rows),
        "all_final_ll_not_below_warm_start_beyond_roundoff": all(
            float(row["ll_refinement_change"]) >= -1e-8 for row in candidate_rows
        ),
        "source_m0_rows_hash_matches_receipt": sha256_path(m0_rows_path)
        == m0_receipt["outputs"]["rows_sha256"],
    }
    receipt = {
        "schema_name": "intersubmod.m0_optimizer_nonconvergence_kkt_audit",
        "schema_version": "1.0.0",
        "scope": {
            "task_type": "B_COMPREHENSIVE_VALIDATION" if full_scope else "A_BOUNDED_DIAGNOSTIC",
            "selected_samples": sorted(selected_samples),
            "n_original_optimizer_abstain_units": len(all_abstain),
            "n_selected_optimizer_abstain_units": len(selected),
            "full_original_abstain_scope": full_scope,
        },
        "method": {
            "objective_unchanged": "sum_pattern count * log(sum_state pi_state * emission(pattern|state))",
            "same_read_vaf_added": False,
            "optimizer": "SLSQP uniform warm start + deterministic monotone pairwise mass transfer",
            "certificate": "Frank-Wolfe/KKT global log-likelihood gap upper bound",
            "concavity": "log of positive affine mixture, nonnegative count sum; any KKT point is global",
            "weight_identifiability": "augmented rank([Q;1]) == number of states; otherwise pi is not uniquely interpretable",
            "claim_ceiling": (
                "certifies numerical optimum of the declared read-pattern mixture only; mixture pi is a latent "
                "state-exposure proportion and is not a cellular clone fraction, true tree, or calibrated model probability"
            ),
        },
        "parameters": {
            "error_rate": args.error_rate,
            "tie_tolerance_log_likelihood": args.tie_tolerance,
            "certificate_tolerance_global_log_likelihood": args.certificate_tolerance,
        },
        "inputs": {
            "m0_receipt": {"path": str(m0_receipt_path), "sha256": sha256_path(m0_receipt_path)},
            "m0_rows": {"path": str(m0_rows_path), "sha256": sha256_path(m0_rows_path)},
            "canonical_root": str(args.canonical_root),
            "frozen_solver": {"path": str(frozen_solver_path), "sha256": sha256_path(frozen_solver_path)},
        },
        "aggregate": {
            "n_units": len(unit_rows),
            "n_candidate_fits": len(candidate_rows),
            "n_slsqp_status_failure_candidates": sum(
                not bool(row["slsqp_success"]) for row in candidate_rows
            ),
            "n_slsqp_status_success_candidates": sum(
                bool(row["slsqp_success"]) for row in candidate_rows
            ),
            "slsqp_status_counts": dict(sorted(slsqp_statuses.items())),
            "n_robust_certified_candidates": sum(
                bool(row["robust_certified"]) for row in candidate_rows
            ),
            "n_fully_rescued_units": sum(
                bool(row["all_candidates_certified"]) for row in unit_rows
            ),
            "n_units_with_changed_top_vertex_set_ids": sum(
                bool(row["top_vertex_set_ids_changed"]) for row in unit_rows
            ),
            "counterfactual_selection_status_counts": dict(sorted(counterfactual.items())),
            "n_rank_deficient_weight_candidates": sum(
                not bool(row["mixture_weights_identifiable"]) for row in candidate_rows
            ),
            "warm_start_global_ll_gap_quantiles": quantiles([
                float(row["warm_start_global_ll_gap_bound"]) for row in candidate_rows
            ]),
            "final_global_ll_gap_quantiles": quantiles([
                float(row["final_global_ll_gap_bound"]) for row in candidate_rows
            ]),
            "absolute_ll_refinement_change_quantiles": quantiles([
                abs(float(row["ll_refinement_change"])) for row in candidate_rows
            ]),
            "refinement_iteration_quantiles": quantiles([
                float(row["refinement_iterations"]) for row in candidate_rows
            ]),
        },
        "checks": checks,
        "outputs": {
            unit_path.name: {"path": str(unit_path), "size_bytes": unit_path.stat().st_size, "sha256": sha256_path(unit_path)},
            candidate_path.name: {"path": str(candidate_path), "size_bytes": candidate_path.stat().st_size, "sha256": sha256_path(candidate_path)},
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "elapsed_seconds": time.monotonic() - started,
        },
    }
    receipt["all_pass"] = all(checks.values())
    receipt_path = args.output_dir / "optimizer_audit_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    raise SystemExit(0 if receipt["all_pass"] else 1)


if __name__ == "__main__":
    main()
