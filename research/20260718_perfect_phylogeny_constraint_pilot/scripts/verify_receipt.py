#!/usr/bin/env python3
"""Independent invariant checks for the bounded perfect-phylogeny pilot."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def rounded_objective(result: dict[str, Any]) -> int | None:
    value = result.get("objective")
    return None if value is None else int(round(float(value)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "receipt",
        nargs="?",
        type=Path,
        default=PROJECT_ROOT / "results" / "pilot_receipt.json",
    )
    args = parser.parse_args()
    with args.receipt.open(encoding="utf-8") as handle:
        receipt = json.load(handle)

    checks = 0

    def require(condition: bool, message: str) -> None:
        nonlocal checks
        checks += 1
        if not condition:
            raise AssertionError(message)

    require(
        receipt["schema"] == "intersubmod.perfect_phylogeny_constraint_pilot.v1",
        "Unexpected receipt schema",
    )
    require(
        receipt["scope_flag"] == "PARTIAL_NOT_COMPREHENSIVE_VALIDATION",
        "Pilot must retain the partial-scope flag",
    )
    require(len(receipt["synthetic_results"]) == 3, "Expected three synthetic cases")
    require(len(receipt["real_results"]) == 8, "Expected eight real cases")

    synthetic = {row["case_key"]: row for row in receipt["synthetic_results"]}
    sisters = synthetic["SYN_sisters_10_01"]
    require(
        sisters["recurrence_allowed"]["status"] == "OPTIMAL"
        and sisters["strict_perfect"]["status"] == "OPTIMAL"
        and rounded_objective(sisters["recurrence_allowed"]) == 0
        and rounded_objective(sisters["strict_perfect"]) == 0,
        "The 10/01 sister control must be feasible with h*=0 in both models",
    )
    recurrence = synthetic["SYN_recurrence_10_01_11"]
    require(
        recurrence["recurrence_allowed"]["status"] == "OPTIMAL"
        and rounded_objective(recurrence["recurrence_allowed"]) == 0
        and recurrence["strict_perfect"]["status"] == "INFEASIBLE",
        "The 10/01/11 recurrence counterexample must separate the models",
    )
    require(
        recurrence["strict_proven_changes_feasibility_or_objective"],
        "The synthetic recurrence case must record a proven model difference",
    )
    partial = synthetic["SYN_partial_AX_XA"]
    require(
        partial["recurrence_allowed"]["status"] == "OPTIMAL"
        and partial["strict_perfect"]["status"] == "OPTIMAL"
        and rounded_objective(partial["recurrence_allowed"]) == 2
        and rounded_objective(partial["strict_perfect"]) == 2,
        "The AX/XA partial-state control must have h*=2 in both models",
    )

    both_certified = 0
    cross_certified = 0
    for row in receipt["real_results"]:
        model_a = row["recurrence_allowed"]
        model_b = row["strict_perfect"]
        require(
            not model_a["selected_vertices"] or model_a["structural_checker_pass"],
            f"Model A incumbent failed structural checks: {row['case_key']}",
        )
        require(
            not model_b["selected_vertices"]
            or (model_b["checker"] is not None and model_b["checker"]["pass"]),
            f"Model B incumbent failed strict checks: {row['case_key']}",
        )
        for label, result in (("A", model_a), ("B", model_b)):
            if result["status"] == "OPTIMAL":
                require(
                    result["objective"] is not None
                    and result["objective_bound"] is not None
                    and math.isclose(
                        float(result["objective"]),
                        float(result["objective_bound"]),
                        abs_tol=1e-6,
                    )
                    and math.isclose(float(result["mip_gap"]), 0.0, abs_tol=1e-9),
                    f"Certified optimum lacks a matching bound: {label} {row['case_key']}",
                )
            if result["status"] == "LIMIT_REACHED":
                require(
                    result["objective"] is not None
                    and result["objective_bound"] is not None
                    and float(result["objective_bound"])
                    <= float(result["objective"]) + 1e-6,
                    f"Limited solve has an invalid incumbent/bound: {label} {row['case_key']}",
                )
        if model_a["status"] == "OPTIMAL" and model_b["status"] == "OPTIMAL":
            both_certified += 1
            require(
                rounded_objective(model_a) == rounded_objective(model_b),
                f"Certified objectives differ: {row['case_key']}",
            )
        if row["strict_optimum_cross_certified_by_model_a"]:
            cross_certified += 1
            require(
                model_a["status"] == "OPTIMAL"
                and model_a["strict_perfect_checker"]["pass"]
                and row["strict_optimum_cross_certified_objective"]
                == rounded_objective(model_a),
                f"Invalid cross-model optimum certificate: {row['case_key']}",
            )
        require(
            not row["strict_proven_changes_feasibility_or_objective"],
            f"Pilot must not claim a proven strict-model difference for real case: "
            f"{row['case_key']}",
        )

    require(
        both_certified >= 4,
        "At least the four M=31/63 controls should be certified by both models",
    )
    require(
        receipt["summary"]["real_both_models_certified_optimal"] == both_certified
        and receipt["summary"]["real_both_certified_same_objective"] == both_certified,
        "Summary disagrees with certified real-case objectives",
    )
    require(
        receipt["summary"][
            "real_strict_proven_changes_feasibility_or_objective"
        ]
        == 0,
        "No real case in this pilot proves a strict-model objective/feasibility change",
    )
    require(
        cross_certified == receipt["summary"]["real_model_a_optimal"]
        and receipt["summary"][
            "real_strict_optimum_cross_certified_by_model_a"
        ]
        == cross_certified,
        "Every certified Model A optimum in this pilot should also certify strict h*",
    )

    require(len(receipt["strict_phase_b"]) == 2, "Expected two bounded Phase B cases")
    for row in receipt["strict_phase_b"]:
        require(
            row["status"] == "INCOMPLETE"
            and not row["complete"]
            and row["n_vertex_sets"] > 0
            and row["parent_tree_count_lower_bound"] >= row["n_vertex_sets"],
            f"Phase B must be explicitly incomplete with a lower bound: {row['case_key']}",
        )

    print(f"PASS checks={checks}")
    print(f"INPUT {args.receipt.resolve()}")
    print(
        "SUMMARY "
        f"synthetic={len(receipt['synthetic_results'])} "
        f"real={len(receipt['real_results'])} "
        f"both_certified={both_certified} "
        f"strict_cross_certified={cross_certified} "
        "real_proven_model_changes=0 "
        "phase_b_complete=0/2"
    )


if __name__ == "__main__":
    main()
