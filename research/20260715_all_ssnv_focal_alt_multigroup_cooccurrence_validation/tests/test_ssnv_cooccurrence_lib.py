from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pysam
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "ssnv_cooccurrence_lib.py"
SPEC = importlib.util.spec_from_file_location("ssnv_cooccurrence_lib", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_user_two_descendant_state_scenario_is_not_overcalled_at_low_depth() -> None:
    # Wild type RR, clone 1 AR, clone 2 AA; forbidden RA is absent.
    focal = ["R"] * 6 + ["A"] * 6 + ["A"] * 6
    partner = ["R"] * 6 + ["R"] * 6 + ["A"] * 6
    summary = MODULE.four_state_summary(focal, partner)
    assert summary["state_counts"] == {"RR": 6, "AR": 6, "RA": 0, "AA": 6}
    assert summary["minimum_zero_violation_depth"] == 203
    assert summary["familywise_confidence"] == pytest.approx(0.95)
    assert summary["confidence"] == pytest.approx(1.0 - 0.05 / 3.0)
    assert summary["multiplicity_method"] == "bonferroni_three_relation_models"
    assert summary["relation_compatibility"] == "NOT_IDENTIFIABLE_FIXED_ERROR_CEILING"

    labels = ["clone1"] * 6 + ["clone2"] * 6
    association = MODULE.methyl_group_marker_association(labels, partner[6:])
    assert association["testable"] is True
    assert association["cramers_v"] == pytest.approx(1.0)
    assert association["delta_alt_fraction"] == pytest.approx(1.0)


def test_user_two_descendant_state_scenario_is_compatible_at_prespecified_depth() -> None:
    # The same RR -> AR -> AA geometry becomes compatible only after the
    # zero-violation denominator clears the fixed 2% simultaneous 95% familywise gate.
    focal = ["R"] * 6 + ["A"] * 6 + ["A"] * 203
    partner = ["R"] * 6 + ["R"] * 6 + ["A"] * 203
    summary = MODULE.four_state_summary(focal, partner)
    assert summary["state_counts"] == {"RR": 6, "AR": 6, "RA": 0, "AA": 203}
    assert summary["focal_ancestor_violation"]["upper_bound"] <= 0.02
    assert summary["relation_compatibility"] == (
        "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    )
    assert summary["compatible_relation_models"] == ["FOCAL_ANCESTOR"]


def test_alt_only_observation_cannot_identify_ancestry() -> None:
    summary = MODULE.four_state_summary(["A"] * 12, ["R"] * 6 + ["A"] * 6)
    assert summary["complete_four_state_testable"] is False
    assert summary["relation_compatibility"] == (
        "NOT_IDENTIFIABLE_INSUFFICIENT_FOUR_STATE_DEPTH"
    )


def test_conditional_permutation_and_bh_are_reproducible() -> None:
    labels = ["g1"] * 6 + ["g2"] * 6
    calls = ["R"] * 6 + ["A"] * 6
    strata = ["HP1:+"] * 4 + ["HP2:-"] * 2 + ["HP1:+"] * 2 + ["HP2:-"] * 4
    first = MODULE.conditional_label_permutation(labels, calls, strata, seed=7, permutations=99)
    second = MODULE.conditional_label_permutation(labels, calls, strata, seed=7, permutations=99)
    assert first == second
    assert first["permutable"] is True
    assert MODULE.benjamini_hochberg([0.01, 0.04, 0.03, None]) == pytest.approx(
        [0.03, 0.04, 0.04, None]
    )


def test_snv_call_respects_base_quality() -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [1000])
    alignment = pysam.AlignedSegment(header)
    alignment.query_name = "read-a"
    alignment.query_sequence = "ACGT"
    alignment.query_qualities = [30, 30, 10, 30]
    alignment.reference_id = 0
    alignment.reference_start = 99
    alignment.cigarstring = "4M"
    assert MODULE.call_snv(alignment, MODULE.Variant("chr1", 100, "A", "T")) == "R"
    assert MODULE.call_snv(alignment, MODULE.Variant("chr1", 101, "A", "C")) == "A"
    assert MODULE.call_snv(alignment, MODULE.Variant("chr1", 102, "C", "G")) == "X"
    assert MODULE.call_snv(alignment, MODULE.Variant("chr1", 103, "A", "C")) == "O"


def test_independent_marker_count_uses_fixed_spacing() -> None:
    assert MODULE.count_independent_positions([100, 105, 120, 141], minimum_separation=20) == 3
