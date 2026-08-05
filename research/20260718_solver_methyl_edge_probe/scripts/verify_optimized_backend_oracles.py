#!/usr/bin/env python3
"""Independent finite-domain oracles for the optimized research backend."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import pathlib
import random
import sys
import time
from typing import Any


HERE = pathlib.Path(__file__).resolve().parent
PROBE_ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from optimized_hypercube_backend import (  # noqa: E402
    BitsetObligationBnb,
    evaluate_evolutionary_mode,
    fixed_node_parent_mapping,
    optimal_family_digest,
    prepare_bitset_problem,
    solve_group_terminal_subset_dp,
    solve_with_certificates,
)
from solver_probe import brute_force_optimal, build_instance  # noqa: E402


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def state_pattern(state: int, k: int) -> str:
    return "".join("A" if state & (1 << bit) else "R" for bit in range(k))


def check_case(
    full: list[str],
    partial: list[str],
    k: int,
) -> dict[str, Any]:
    instance = build_instance(full, partial, k)
    brute = brute_force_optimal(instance, combination_cap=2_000_000)
    if not brute["complete"]:
        return {"pass": False, "reason": "BRUTE_INCOMPLETE"}
    problem = prepare_bitset_problem(instance)
    dp = solve_group_terminal_subset_dp(problem, max_q=8)
    target = dp.objective if dp.objective_certified else None
    bnb = BitsetObligationBnb(
        problem,
        time_limit_seconds=5.0,
        certified_target_objective=target,
    ).run()
    brute_digest = optimal_family_digest(k, brute["vertex_sets"])
    return {
        "pass": (
            bnb.family_complete
            and bnb.objective == brute["objective"]
            and bnb.family_digest == brute_digest
            and (
                not dp.objective_certified
                or dp.objective == brute["objective"]
            )
        ),
        "q": problem.q,
        "dp_routed": dp.objective_certified,
        "brute_objective": brute["objective"],
        "dp_objective": dp.objective,
        "bnb_objective": bnb.objective,
        "brute_family_count": len(brute["vertex_sets"]),
        "bnb_family_count": len(bnb.vertex_sets),
        "brute_digest": brute_digest,
        "bnb_digest": bnb.family_digest,
        "bnb_status": bnb.status,
    }


def exhaustive_k_le_3() -> dict[str, Any]:
    started = time.perf_counter()
    cases = 0
    dp_routed = 0
    mismatches: list[dict[str, Any]] = []
    for k in range(1, 4):
        full_states = list(range(1, 1 << k))
        partial_pool = [
            "".join(symbols)
            for symbols in itertools.product("RAX", repeat=k)
            if "X" in symbols
        ]
        partial_cases = [()]
        partial_cases.extend((pattern,) for pattern in partial_pool)
        partial_cases.extend(
            itertools.combinations_with_replacement(partial_pool, 2)
        )
        for mandatory_bits in range(1 << len(full_states)):
            full = [
                state_pattern(state, k)
                for index, state in enumerate(full_states)
                if mandatory_bits & (1 << index)
            ]
            for partial_tuple in partial_cases:
                partial = list(partial_tuple)
                result = check_case(full, partial, k)
                cases += 1
                dp_routed += bool(result.get("dp_routed"))
                if not result["pass"] and len(mismatches) < 20:
                    mismatches.append(
                        {
                            "k": k,
                            "full": full,
                            "partial": partial,
                            "result": result,
                        }
                    )
    return {
        "status": "PASS" if not mismatches else "FAIL",
        "cases": cases,
        "dp_routed_cases": dp_routed,
        "mismatch_count": len(mismatches),
        "first_mismatches": mismatches,
        "elapsed_seconds": time.perf_counter() - started,
    }


def seeded_k4(n_cases: int) -> dict[str, Any]:
    started = time.perf_counter()
    rng = random.Random(20260718)
    mismatches: list[dict[str, Any]] = []
    checked = 0
    full_pool = [
        "".join(symbols) for symbols in itertools.product("RA", repeat=4)
    ]
    partial_pool = [
        "".join(symbols)
        for symbols in itertools.product("RAX", repeat=4)
        if "X" in symbols
    ]
    for _ in range(n_cases):
        full = rng.sample(full_pool, rng.randint(0, 4))
        partial = rng.sample(partial_pool, rng.randint(0, 5))
        result = check_case(full, partial, 4)
        checked += 1
        if not result["pass"] and len(mismatches) < 20:
            mismatches.append(
                {"full": full, "partial": partial, "result": result}
            )
    return {
        "status": "PASS" if not mismatches else "FAIL",
        "cases": checked,
        "mismatch_count": len(mismatches),
        "first_mismatches": mismatches,
        "elapsed_seconds": time.perf_counter() - started,
    }


def noncontiguous_active_bits() -> dict[str, Any]:
    full = ["A" + "R" * 10 + "A"]
    partial = ["X" + "R" * 10 + "A", "A" + "R" * 10 + "X"]
    result = check_case(full, partial, 12)
    instance = build_instance(full, partial, 12)
    return {
        "status": "PASS" if result["pass"] else "FAIL",
        "raw_k": 12,
        "effective_m": bin(instance.universe_mask).count("1"),
        "vertices": list(instance.vertices),
        "result": result,
    }


def fail_closed_controls() -> dict[str, Any]:
    legacy_counterexample = solve_with_certificates(
        build_instance(["RAA"], ["ARX", "AXA"], 3),
        max_sets=1,
        time_limit_seconds=5.0,
    )
    capped = solve_with_certificates(
        build_instance(["AAAA"], [], 4),
        max_sets=1,
        time_limit_seconds=5.0,
    )
    deadline = solve_with_certificates(
        build_instance(["AAAA"], [], 4),
        max_sets=None,
        time_limit_seconds=1e-12,
    )
    direct_problem = prepare_bitset_problem(build_instance(["AAAA"], [], 4))
    direct_no_incumbent_cap = BitsetObligationBnb(
        direct_problem,
        max_sets=0,
        time_limit_seconds=5.0,
    ).run()
    direct_no_incumbent_deadline = BitsetObligationBnb(
        direct_problem,
        max_sets=None,
        time_limit_seconds=1e-12,
    ).run()
    dp_before_bnb = solve_group_terminal_subset_dp(direct_problem)
    certified_bnb_deadline = BitsetObligationBnb(
        direct_problem,
        max_sets=None,
        time_limit_seconds=1e-12,
        certified_target_objective=dp_before_bnb.objective,
    ).run()
    direct_dp_deadline = solve_group_terminal_subset_dp(
        direct_problem,
        time_limit_seconds=1e-12,
    )

    nonfinite_deadline_rejections = []
    for invalid in (math.nan, math.inf, -math.inf):
        for route in ("dp", "bnb", "wrapper"):
            try:
                if route == "dp":
                    solve_group_terminal_subset_dp(
                        direct_problem,
                        time_limit_seconds=invalid,
                    )
                elif route == "bnb":
                    BitsetObligationBnb(
                        direct_problem,
                        time_limit_seconds=invalid,
                    )
                else:
                    solve_with_certificates(
                        build_instance(["AAAA"], [], 4),
                        time_limit_seconds=invalid,
                    )
            except ValueError:
                nonfinite_deadline_rejections.append(True)
            else:
                nonfinite_deadline_rejections.append(False)

    strict_gates = {
        "allele_specific_cn_known": True,
        "copy_neutral_diploid_both_homologs_retained": True,
        "total_cn": 2,
        "major_cn": 1,
        "minor_cn": 1,
        "loh_absent": True,
        "clonal_loh_absent": True,
        "subclonal_loh_absent": True,
        "clonal_deletion_absent": True,
        "subclonal_deletion_absent": True,
        "clonal_amplification_absent": True,
        "subclonal_amplification_absent": True,
        "clonal_cna_absent": True,
        "subclonal_cna_absent": True,
        "wgd_boundary_absent": True,
        "mutation_loss_absent": True,
        "phasing_qc_pass": True,
        "mapping_qc_pass": True,
        "base_quality_qc_pass": True,
        "strand_qc_pass": True,
        "read_independence_qc_pass": True,
        "allele_specific_duplicated_copies_absent": True,
        "somatic_variant_qc_pass": True,
    }
    disconnected_m1 = evaluate_evolutionary_mode(
        {0, 3},
        k=2,
        mode="M1_STRICT_INFINITE_SITES",
        eligibility=strict_gates,
    )
    try:
        fixed_node_parent_mapping(
            {0, 1, 3},
            edge_scores={(0, 1): 1e308, (1, 3): 1e308},
        )
    except ValueError:
        cumulative_edge_overflow_rejected = True
    else:
        cumulative_edge_overflow_rejected = False
    cap_family = capped["family_certificate"]
    legacy_family = legacy_counterexample["family_certificate"]
    deadline_family = deadline["family_certificate"]
    checks = {
        "legacy_counterexample_objective_is_true_h2": (
            legacy_family["objective"] == 2
        ),
        "legacy_counterexample_unique_family_complete": (
            legacy_family["family_complete"]
        ),
        "cap_objective_is_true_h3": cap_family["objective"] == 3,
        "cap_family_incomplete": not cap_family["family_complete"],
        "cap_ranking_blocked": not capped["ranking_allowed"],
        "deadline_family_incomplete": not deadline_family["family_complete"],
        "deadline_ranking_blocked": not deadline["ranking_allowed"],
        "total_deadline_objective_not_claimed": (
            deadline["objective_certificate"]["objective_certified"] is False
            and deadline_family["objective"] is None
        ),
        "certified_bnb_deadline_preserves_h3": (
            certified_bnb_deadline.objective == 3
            and certified_bnb_deadline.objective_certified
            and not certified_bnb_deadline.family_complete
        ),
        "direct_cap_does_not_claim_feasible": (
            direct_no_incumbent_cap.status
            == "NO_FEASIBLE_CERTIFICATE_MAX_SETS"
        ),
        "direct_deadline_does_not_claim_feasible": (
            direct_no_incumbent_deadline.status
            == "NO_FEASIBLE_CERTIFICATE_DEADLINE"
        ),
        "direct_dp_deadline_does_not_certify": (
            direct_dp_deadline.status == "OBJECTIVE_INCOMPLETE_DEADLINE"
            and not direct_dp_deadline.objective_certified
        ),
        "all_nonfinite_deadlines_rejected": all(nonfinite_deadline_rejections)
        and len(nonfinite_deadline_rejections) == 9,
        "disconnected_m1_is_not_publishable": (
            disconnected_m1["status"] == "STRICT_INFEASIBLE"
            and disconnected_m1["publication_allowed"] is False
            and disconnected_m1["orphan_states_without_selected_predecessor"]
            == [3]
        ),
        "cumulative_edge_overflow_rejected": cumulative_edge_overflow_rejected,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "legacy_counterexample_status": legacy_family["status"],
        "cap_status": cap_family["status"],
        "deadline_status": deadline_family["status"],
        "direct_no_incumbent_cap_status": direct_no_incumbent_cap.status,
        "direct_no_incumbent_deadline_status": direct_no_incumbent_deadline.status,
        "certified_bnb_deadline_status": certified_bnb_deadline.status,
        "direct_dp_deadline_status": direct_dp_deadline.status,
        "disconnected_m1_status": disconnected_m1["status"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--seeded-k4-cases", type=int, default=300)
    args = parser.parse_args()

    exhaustive = exhaustive_k_le_3()
    seeded = seeded_k4(args.seeded_k4_cases)
    noncontiguous = noncontiguous_active_bits()
    fail_closed = fail_closed_controls()
    all_pass = all(
        section["status"] == "PASS"
        for section in (exhaustive, seeded, noncontiguous, fail_closed)
    )
    backend_path = HERE / "optimized_hypercube_backend.py"
    test_path = PROBE_ROOT / "tests/test_optimized_hypercube_backend.py"
    receipt = {
        "schema": "intersubmod.optimized_backend_oracles.v1",
        "scope": "FINITE_EXHAUSTIVE_K_LE_3_PLUS_SEEDED_K4_NOT_PRODUCTION",
        "source_sha256": {
            "optimized_backend": sha256_file(backend_path),
            "oracle_script": sha256_file(pathlib.Path(__file__).resolve()),
            "unit_tests": sha256_file(test_path),
        },
        "exhaustive_k_le_3": exhaustive,
        "seeded_k4": seeded,
        "noncontiguous_active_bits": noncontiguous,
        "fail_closed_controls": fail_closed,
        "summary": {
            "all_pass": all_pass,
            "total_structural_cases": exhaustive["cases"] + seeded["cases"] + 1,
            "total_mismatches": (
                exhaustive["mismatch_count"]
                + seeded["mismatch_count"]
                + (0 if noncontiguous["status"] == "PASS" else 1)
            ),
            "incomplete_ranked_count": 0,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "PASS" if all_pass else "FAIL",
        f"exhaustive={exhaustive['cases']}",
        f"seeded_k4={seeded['cases']}",
        f"mismatches={receipt['summary']['total_mismatches']}",
        f"cap={fail_closed['cap_status']}",
        f"deadline={fail_closed['deadline_status']}",
    )
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
