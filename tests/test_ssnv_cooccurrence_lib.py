from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from scipy.stats import fisher_exact


LIB_PATH = Path(__file__).resolve().parents[1] / "scripts" / "ssnv_cooccurrence_lib.py"
SPEC = importlib.util.spec_from_file_location("test_shared_ssnv_cooccurrence_lib", LIB_PATH)
assert SPEC is not None and SPEC.loader is not None
LIB = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = LIB
SPEC.loader.exec_module(LIB)


def test_marker_association_pearson_path_is_scipy_tuple_compatible() -> None:
    labels = ["g1"] * 10 + ["g2"] * 10 + ["g3"] * 10
    calls = (
        ["R"] * 8
        + ["A"] * 2
        + ["R"] * 5
        + ["A"] * 5
        + ["R"] * 2
        + ["A"] * 8
    )

    observed = LIB.methyl_group_marker_association(labels, calls)
    expected = LIB.chi2_contingency(observed["table"], correction=False)

    assert observed["testable"] is True
    assert observed["analytic_test"] == "pearson_chi_square_descriptive_only"
    assert observed["p_analytic"] == pytest.approx(float(expected[1]))


def test_fixed_margins_exact_matches_fisher_for_2x2() -> None:
    table = [[8, 2], [1, 9]]
    result = LIB.fisher_freeman_halton_kx2(table)
    assert result["identifiable"] is True
    assert result["state_space_status"] == "EXACT_ENUMERATED"
    assert result["p_value"] == pytest.approx(float(fisher_exact(table)[1]), rel=1e-12)


def test_sparse_kx2_uses_exact_probability_ordering() -> None:
    result = LIB.fisher_freeman_halton_kx2([[3, 0], [0, 3], [2, 1]])
    assert result["state_space_size"] == 12
    assert result["p_value"] == pytest.approx(1.0 / 7.0)
    assert result["method"] == "fisher_freeman_halton_fixed_margins_probability_ordering"


def test_exact_state_space_ceiling_fails_closed_without_chi_square_fallback() -> None:
    result = LIB.fisher_freeman_halton_kx2(
        [[3, 3], [3, 3], [3, 3]], state_space_ceiling=2
    )
    assert result["p_value"] is None
    assert result["identifiable"] is False
    assert result["state_space_status"] == "NOT_IDENTIFIABLE_STATE_SPACE_LIMIT"
    assert result["state_space_lower_bound"] == 3
    assert result["fallback"] == "fail_closed_not_identifiable"


def test_conditional_permutation_is_deterministic_for_sparse_kx2() -> None:
    labels = ["g1"] * 4 + ["g2"] * 4 + ["g3"] * 4
    categories = ["R", "R", "R", "A"] + ["R", "A", "A", "A"] + ["R", "R", "A", "A"]
    strata = ["HP1|PS=10|strand=+"] * 6 + ["HP2|PS=10|strand=-"] * 6
    first = LIB.conditional_label_permutation(labels, categories, strata, seed=17, permutations=99)
    second = LIB.conditional_label_permutation(labels, categories, strata, seed=17, permutations=99)
    assert first == second
    assert first["status"] == "PERMUTABLE"
    assert first["permutations"] == 99


def test_conditional_permutation_nonexchangeable_is_not_identifiable() -> None:
    result = LIB.conditional_label_permutation(
        ["g1"] * 6 + ["g2"] * 6,
        ["R"] * 6 + ["A"] * 6,
        ["HP1|PS=10|strand=+"] * 6 + ["HP2|PS=10|strand=+"] * 6,
        seed=3,
        permutations=999,
    )
    assert result["p_conditional_perm"] is None
    assert result["permutable"] is False
    assert result["status"] == "NOT_IDENTIFIABLE_NONEXCHANGEABLE"


def test_low_depth_25_percent_forbidden_state_is_not_accepted() -> None:
    focal = ["R"] * 3 + ["A"] * 3 + ["R"] + ["A"] * 3
    partner = ["R"] * 3 + ["R"] * 3 + ["A"] + ["A"] * 3
    result = LIB.four_state_summary(focal, partner)
    focal_gate = result["focal_ancestor_violation"]
    assert result["state_counts_with_noncallable"] == {
        "RR": 3,
        "AR": 3,
        "RA": 1,
        "AA": 3,
        "O": 0,
        "X": 0,
    }
    assert focal_gate["rate"] == pytest.approx(0.25)
    assert focal_gate["upper_bound"] > focal_gate["threshold"]
    assert focal_gate["status"] == "NOT_IDENTIFIABLE_LOW_PRECISION"
    assert result["relation_compatibility"].startswith("NOT_IDENTIFIABLE")
    assert result["minimum_zero_violation_depth"] == 203


def test_zero_forbidden_reads_requires_fixed_ceiling_depth() -> None:
    focal = ["R"] * 50 + ["A"] * 50 + ["A"] * 203
    partner = ["R"] * 50 + ["R"] * 50 + ["A"] * 203
    result = LIB.four_state_summary(focal, partner)
    gate = result["focal_ancestor_violation"]
    assert gate["violations"] == 0
    assert gate["denominator"] == 203
    assert gate["upper_bound"] <= 0.02
    assert gate["status"] == "COMPATIBLE_WITH_FIXED_ERROR_CEILING"
    assert result["relation_compatibility"] == "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    assert result["compatible_relation_models"] == ["FOCAL_ANCESTOR"]
    assert result["n_compatible_relation_models"] == 1


def test_multiple_compatible_order_models_are_not_hidden_by_priority() -> None:
    rr, ar, ra, aa = 3, 40_000, 3, 500
    focal = ["R"] * rr + ["A"] * ar + ["R"] * ra + ["A"] * aa
    partner = ["R"] * rr + ["R"] * ar + ["A"] * ra + ["A"] * aa
    result = LIB.four_state_summary(focal, partner)
    assert result["compatible_relation_models"] == ["FOCAL_ANCESTOR", "BRANCHING"]
    assert result["n_compatible_relation_models"] == 2
    assert result["relation_compatibility"] == (
        "MULTIPLE_MUTATION_ORDER_MODELS_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    )


def test_spatial_separation_is_not_named_statistical_independence() -> None:
    assert LIB.spatially_separated_positions([100, 105, 120, 141], 20) == [100, 120, 141]
    assert LIB.count_spatially_separated_positions([100, 105, 120, 141], 20) == 3
    assert "Deprecated alias" in (LIB.count_independent_positions.__doc__ or "")
