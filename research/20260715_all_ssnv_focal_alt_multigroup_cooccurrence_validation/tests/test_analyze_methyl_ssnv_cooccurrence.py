from __future__ import annotations

import csv
import importlib.util
import gzip
import json
import subprocess
import sys
from array import array
from dataclasses import replace
from pathlib import Path

import pysam
import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
SCRIPT = SCRIPT_DIR / "analyze_methyl_ssnv_cooccurrence.py"
SPEC = importlib.util.spec_from_file_location("analyze_methyl_ssnv_cooccurrence", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


@pytest.fixture(autouse=True)
def test_only_release_source_authority(monkeypatch: pytest.MonkeyPatch) -> None:
    authority = {"authority_id": "TEST_ONLY_UNSIGNED_AUTHORITY", "pass": True}
    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY,
        "validate_release_source_authority",
        lambda *_args, **_kwargs: authority,
    )


def make_alignment(
    header: pysam.AlignmentHeader,
    name: str,
    *,
    start: int = 100,
    flag: int = 0,
    mapq: int = 60,
    length: int = 1_200,
    add_methyl_tags: bool = True,
) -> pysam.AlignedSegment:
    alignment = pysam.AlignedSegment(header)
    alignment.query_name = name
    alignment.query_sequence = "C" + "A" * (length - 1)
    alignment.query_qualities = [35] * length
    alignment.flag = flag
    alignment.reference_id = 0
    alignment.reference_start = start
    alignment.mapping_quality = mapq
    alignment.cigarstring = f"{length}M"
    if add_methyl_tags:
        alignment.set_tag("MM", "C+m?,0;")
        alignment.set_tag("ML", array("B", [255]))
    return alignment


def write_bam(path: Path, alignments: list[pysam.AlignedSegment]) -> Path:
    header = alignments[0].header.to_dict()
    with pysam.AlignmentFile(path, "wb", header=header) as handle:
        for alignment in alignments:
            handle.write(alignment)
    pysam.index(str(path))
    return path


def projection(name: str, *, start: int = 100, end: int = 1_300, mapq: int = 60) -> tuple:
    return MODULE.TAGS.projection_key(name, "chr1", start, end, mapq, "+")


def populate_empty_raw_identity_audit(audit: dict, expected: int) -> None:
    audit.update(
        {
            "equivalence_contract": MODULE.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
            "allowed_differing_auxiliary_tags": ["RG"],
            "missing_projection_policy": MODULE.RAW_IDENTITY_MISSING_POLICY,
            "conflicting_analysis_payload_policy": MODULE.RAW_IDENTITY_CONFLICT_POLICY,
            "expected_projections": expected,
            "matched_records": expected,
            "duplicate_projections_collapsed": 0,
            "duplicate_extra_records_collapsed": 0,
            "exact_duplicate_projections_collapsed": 0,
            "rg_only_duplicate_projections_collapsed": 0,
            "duplicate_projection_sha256": MODULE.projection_digest([]),
            "alignment_payload_sha256": "0" * 64,
            "analysis_scope_identity_contract": MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT,
            "recovered_payload_sha256": "1" * 64,
            "analysis_site_payload_sha256": "2" * 64,
            "duplicate_projection_examples": [],
            "duplicate_projection_records": [],
        }
    )


def test_raw_bam_identity_resolves_one_eligible_primary_alignment(tmp_path: Path) -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    primary = make_alignment(header, "target")
    secondary_same_identity = make_alignment(header, "target", flag=0x100)
    duplicate_same_identity = make_alignment(header, "target", flag=0x400)
    low_mapq = make_alignment(header, "low-mapq", mapq=19)
    too_short = make_alignment(header, "too-short", length=999)
    no_mm_ml = make_alignment(header, "no-tags", add_methyl_tags=False)
    bam_path = write_bam(
        tmp_path / "unique.bam",
        [primary, secondary_same_identity, duplicate_same_identity, low_mapq, too_short, no_mm_ml],
    )
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        matches = MODULE.resolve_unique_eligible_alignments(bam, "chr1", 101, [projection("target")])
    assert list(matches) == [projection("target")]
    assert matches[projection("target")].flag == 0
    assert MODULE.cigar_query_length(primary) == 1_200
    assert MODULE.eligible_alignment(primary) is True
    assert MODULE.eligible_alignment(secondary_same_identity) is False
    assert MODULE.eligible_alignment(low_mapq) is False
    assert MODULE.eligible_alignment(too_short) is False
    assert MODULE.eligible_alignment(no_mm_ml) is False


def test_raw_bam_identity_exact_duplicate_collapses_and_missing_fails(tmp_path: Path) -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    first = make_alignment(header, "duplicate-projection")
    second = make_alignment(header, "duplicate-projection")
    bam_path = write_bam(tmp_path / "ambiguous.bam", [first, second])
    audit = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        matches = MODULE.resolve_unique_eligible_alignments(
            bam,
            "chr1",
            101,
            [projection("duplicate-projection")],
            audit,
        )
    assert list(matches) == [projection("duplicate-projection")]
    assert audit["duplicate_projections_collapsed"] == 1
    assert audit["duplicate_extra_records_collapsed"] == 1
    assert audit["exact_duplicate_projections_collapsed"] == 1
    assert audit["rg_only_duplicate_projections_collapsed"] == 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with pytest.raises(RuntimeError, match="missing"):
            MODULE.resolve_unique_eligible_alignments(bam, "chr1", 101, [projection("absent")])


def test_raw_bam_identity_rg_only_duplicate_collapses_with_audit(tmp_path: Path) -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    first = make_alignment(header, "rg-duplicate")
    second = make_alignment(header, "rg-duplicate")
    first.set_tag("RG", "run-a")
    second.set_tag("RG", "run-b")
    bam_path = write_bam(tmp_path / "rg-only.bam", [first, second])
    audit = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        matches = MODULE.resolve_unique_eligible_alignments(
            bam, "chr1", 101, [projection("rg-duplicate")], audit
        )
    assert list(matches) == [projection("rg-duplicate")]
    assert audit["equivalence_contract"] == MODULE.RAW_DUPLICATE_EQUIVALENCE_CONTRACT
    assert audit["allowed_differing_auxiliary_tags"] == ["RG"]
    assert audit["exact_duplicate_projections_collapsed"] == 0
    assert audit["rg_only_duplicate_projections_collapsed"] == 1


def test_raw_bam_identity_three_rg_duplicates_collapse_with_exact_extra_count(
    tmp_path: Path,
) -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    alignments = [make_alignment(header, "triple-rg") for _ in range(3)]
    alignments[0].set_tag("RG", "run-a")
    alignments[1].set_tag("RG", "run-b")
    alignments[2].set_tag("RG", "run-b")
    bam_path = write_bam(tmp_path / "triple-rg.bam", alignments)
    audit = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        MODULE.resolve_unique_eligible_alignments(
            bam, "chr1", 101, [projection("triple-rg")], audit
        )
    assert audit["duplicate_projections_collapsed"] == 1
    assert audit["duplicate_extra_records_collapsed"] == 2
    assert audit["rg_only_duplicate_projections_collapsed"] == 1
    assert audit["duplicate_projection_records"][0]["record_count"] == 3


def test_raw_bam_identity_rg_present_vs_absent_is_allowed_and_audited(
    tmp_path: Path,
) -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    first = make_alignment(header, "rg-present-absent")
    second = make_alignment(header, "rg-present-absent")
    second.set_tag("RG", "run-b")
    bam_path = write_bam(tmp_path / "rg-present-absent.bam", [first, second])
    audit = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        MODULE.resolve_unique_eligible_alignments(
            bam, "chr1", 101, [projection("rg-present-absent")], audit
        )
    assert audit["duplicate_projection_records"][0]["classification"] == (
        "RG_ONLY_DUPLICATE"
    )
    assert audit["duplicate_projection_records"][0]["differing_auxiliary_tags"] == [
        "RG"
    ]


@pytest.mark.parametrize(
    "conflict",
    ["sequence", "quality", "cigar", "flag", "MM", "ML", "HP", "PS", "AS"],
)
def test_raw_bam_identity_material_duplicate_conflicts_fail_closed(
    tmp_path: Path, conflict: str
) -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    first = make_alignment(header, "conflicting-duplicate")
    second = make_alignment(header, "conflicting-duplicate")
    if conflict == "sequence":
        second.query_sequence = "G" + "A" * 1_199
        second.query_qualities = [35] * 1_200
    elif conflict == "quality":
        second.query_qualities = [34] * 1_200
    elif conflict == "cigar":
        second.cigarstring = "600M1I1D599M"
    elif conflict == "flag":
        second.flag = 0x2
    elif conflict == "MM":
        second.set_tag("MM", "C+m?,1;")
    elif conflict == "ML":
        second.set_tag("ML", array("B", [254]))
    elif conflict == "HP":
        first.set_tag("HP", 1)
        second.set_tag("HP", 2)
    elif conflict == "PS":
        first.set_tag("PS", 10)
        second.set_tag("PS", 11)
    elif conflict == "AS":
        first.set_tag("AS", 1)
        second.set_tag("AS", 2)
    bam_path = write_bam(tmp_path / f"conflict-{conflict}.bam", [first, second])
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with pytest.raises(RuntimeError, match="conflicting_analysis_payload"):
            MODULE.resolve_unique_eligible_alignments(
                bam, "chr1", 101, [projection("conflicting-duplicate")]
            )


def test_raw_bam_identity_preserves_b_array_subtype_after_bam_roundtrip(
    tmp_path: Path,
) -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    signed = make_alignment(header, "typed-ml")
    unsigned = make_alignment(header, "typed-ml")
    signed.set_tag("ML", array("b", [1]))
    unsigned.set_tag("ML", array("B", [1]))
    bam_path = write_bam(tmp_path / "typed-ml.bam", [signed, unsigned])
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        observed = list(bam.fetch("chr1", 100, 101))
        assert observed[0].get_tag("ML").typecode == "b"
        assert observed[1].get_tag("ML").typecode == "B"
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with pytest.raises(RuntimeError, match="conflicting_analysis_payload"):
            MODULE.resolve_unique_eligible_alignments(
                bam, "chr1", 101, [projection("typed-ml")]
            )


def test_analysis_scope_alignment_digest_ignores_rg_but_binds_sequence_and_ml() -> None:
    header = pysam.AlignmentHeader.from_references(["chr1"], [10_000])
    baseline = make_alignment(header, "payload")
    key = projection("payload")
    baseline_digest = MODULE.alignment_payload_digest({key: baseline})

    rg_changed = make_alignment(header, "payload")
    rg_changed.set_tag("RG", "another-run")
    assert MODULE.alignment_payload_digest({key: rg_changed}) == baseline_digest

    sequence_changed = make_alignment(header, "payload")
    sequence_changed.query_sequence = "G" + "A" * 1_199
    sequence_changed.query_qualities = [35] * 1_200
    assert MODULE.alignment_payload_digest({key: sequence_changed}) != baseline_digest

    ml_changed = make_alignment(header, "payload")
    ml_changed.set_tag("ML", array("B", [254]))
    assert MODULE.alignment_payload_digest({key: ml_changed}) != baseline_digest


def test_analysis_scope_recovered_digest_binds_latest_tags_label_and_partner_calls() -> None:
    baseline = recovered_read(1, "g1", "A")
    baseline_digest = MODULE.canonical_payload_sha256(
        MODULE.recovered_read_analysis_payload(baseline)
    )
    mutations = (
        replace(baseline, latest_hp="1-1"),
        replace(baseline, latest_ps=202),
        replace(baseline, label="g2"),
        replace(baseline, partner_calls={200: "R"}),
    )
    assert all(
        MODULE.canonical_payload_sha256(
            MODULE.recovered_read_analysis_payload(mutated)
        )
        != baseline_digest
        for mutated in mutations
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("matched_records", 1),
        ("duplicate_extra_records_collapsed", 0),
        ("exact_duplicate_projections_collapsed", 1),
        ("duplicate_projections_collapsed", 2),
        ("duplicate_projection_sha256", "not-a-sha256"),
    ],
)
def test_raw_identity_audit_invariants_fail_closed(field: str, value: object) -> None:
    key = projection("audit")
    audit = {
        "equivalence_contract": MODULE.RAW_DUPLICATE_EQUIVALENCE_CONTRACT,
        "allowed_differing_auxiliary_tags": ["RG"],
        "missing_projection_policy": MODULE.RAW_IDENTITY_MISSING_POLICY,
        "conflicting_analysis_payload_policy": MODULE.RAW_IDENTITY_CONFLICT_POLICY,
        "expected_projections": 1,
        "matched_records": 2,
        "duplicate_projections_collapsed": 1,
        "duplicate_extra_records_collapsed": 1,
        "exact_duplicate_projections_collapsed": 0,
        "rg_only_duplicate_projections_collapsed": 1,
        "duplicate_projection_sha256": MODULE.projection_digest([key]),
        "alignment_payload_sha256": "0" * 64,
        "duplicate_projection_examples": [list(key)],
        "duplicate_projection_records": [
            {
                "projection": list(key),
                "record_count": 2,
                "classification": "RG_ONLY_DUPLICATE",
                "differing_auxiliary_tags": ["RG"],
            }
        ],
    }
    audit[field] = value
    with pytest.raises(RuntimeError):
        MODULE.validate_raw_identity_audit(audit, expected_recovered_count=1)


def test_raw_identity_audit_missing_field_and_empty_digest_mismatch_fail_closed() -> None:
    audit = {}
    populate_empty_raw_identity_audit(audit, 1)
    del audit["matched_records"]
    with pytest.raises(RuntimeError, match="omitted"):
        MODULE.validate_raw_identity_audit(audit, expected_recovered_count=1)

    populate_empty_raw_identity_audit(audit, 1)
    audit["duplicate_projection_sha256"] = "a" * 64
    with pytest.raises(RuntimeError, match="digest"):
        MODULE.validate_raw_identity_audit(audit, expected_recovered_count=1)


def recovered_read(index: int, label: str | None, partner_call: str, focal_call: str = "A"):
    return MODULE.RecoveredRead(
        read_id=f"read-{index}",
        read_name=f"read-{index}",
        expected_focal_call=focal_call,
        focal_call=focal_call,
        latest_hp="1-1" if index % 2 == 0 else "2-1",
        latest_ps=101,
        strand="+" if index % 3 else "-",
        label=label,
        partner_calls={200: partner_call},
    )


def test_endpoint_a_keeps_r_a_o_x_and_never_converts_o_x_to_reference() -> None:
    calls = ["R"] * 5 + ["O"] * 3 + ["X"] * 2 + ["A"] * 5 + ["O"] * 3 + ["X"] * 2
    labels = ["g1"] * 10 + ["g2"] * 10
    reads = [recovered_read(index, label, call) for index, (label, call) in enumerate(zip(labels, calls))]
    row = MODULE._pair_result(
        {"sample": "S", "biological_id": "B"},
        MODULE.LIB.Variant("chr1", 100, "A", "C"),
        MODULE.LIB.Variant("chr1", 200, "A", "G"),
        reads,
    )
    assert row["core_partner_call_counts"] == {"A": 5, "O": 6, "R": 5, "X": 4}
    assert row["endpoint_a_n_informative"] == 10
    assert row["endpoint_a_table"] == [[5, 0], [0, 5]]
    assert row["endpoint_a_testable"] is True
    assert row["endpoint_a_cramers_v"] == pytest.approx(1.0)
    assert row["endpoint_a_delta_alt_fraction"] == pytest.approx(1.0)
    assert row["callability_table"] == [[5, 5], [5, 5]]
    assert set(row).difference({"_permutation_payload"}) == set(MODULE.PAIR_COLUMNS)


def test_ffh_fixed_margins_golden_kx2_probability_ordering() -> None:
    result = MODULE.LIB.fisher_freeman_halton_kx2([[3, 0], [3, 0], [2, 1], [0, 3]])
    assert result["identifiable"] is True
    assert result["state_space_status"] == "EXACT_ENUMERATED"
    assert result["p_value"] == pytest.approx(0.0727272727, rel=1e-9)


def test_ffh_state_space_limit_fails_closed_without_analytic_fallback() -> None:
    result = MODULE.LIB.fisher_freeman_halton_kx2(
        [[3, 3], [3, 3], [3, 3]], state_space_ceiling=2
    )
    assert result["p_value"] is None
    assert result["identifiable"] is False
    assert result["state_space_status"] == "NOT_IDENTIFIABLE_STATE_SPACE_LIMIT"
    assert result["state_space_lower_bound"] == 3
    assert result["fallback"] == "fail_closed_not_identifiable"


def test_four_state_endpoint_reports_all_states_and_alt_only_is_not_identifiable() -> None:
    focal = ["R"] * 6 + ["A"] * 6 + ["R"] * 3 + ["A"] * 3
    partner = ["R"] * 6 + ["R"] * 6 + ["A"] * 3 + ["A"] * 3
    summary = MODULE.LIB.four_state_summary(focal, partner)
    assert summary["state_counts"] == {"RR": 6, "AR": 6, "RA": 3, "AA": 3}
    assert summary["complete_four_state_testable"] is True
    assert MODULE.endpoint_b_status(summary) == summary["relation_compatibility"]

    alt_only = MODULE.LIB.four_state_summary(["A"] * 12, ["R"] * 6 + ["A"] * 6)
    assert MODULE.endpoint_b_status(alt_only) == "NOT_IDENTIFIABLE_NO_FOCAL_REF"
    assert "not proof" in alt_only["claim_guardrail"]

    low_depth_forbidden = MODULE.LIB.four_state_summary(
        ["R"] * 3 + ["A"] * 3 + ["R"] + ["A"] * 3,
        ["R"] * 3 + ["R"] * 3 + ["A"] + ["A"] * 3,
    )
    assert low_depth_forbidden["state_counts"] == {"RR": 3, "AR": 3, "RA": 1, "AA": 3}
    assert low_depth_forbidden["focal_ancestor_violation"]["rate"] == pytest.approx(0.25)
    assert low_depth_forbidden["focal_ancestor_violation"]["status"] == (
        "NOT_IDENTIFIABLE_LOW_PRECISION"
    )
    assert low_depth_forbidden["relation_compatibility"].startswith("NOT_IDENTIFIABLE")


def test_global_exact_bh_by_uses_only_m2_family_and_separates_effect_gate() -> None:
    rows = [
        {
            "endpoint_a_testable": True,
            "endpoint_a_p_analytic": 1e-12,
            "endpoint_a_p_fixed_margins_exact": 0.001,
            "endpoint_a_global_fdr_family_status": "ELIGIBLE_M2_EXACT_FAMILY",
            "endpoint_a_cramers_v": 0.31,
            "endpoint_a_delta_alt_fraction": 0.50,
            "callability_p_analytic": 0.20,
            "callability_testable": False,
            "callability_cramers_v": None,
            "core_partner_call_counts": {"R": 10, "A": 10},
        },
        {
            "endpoint_a_testable": True,
            "endpoint_a_p_analytic": 1e-12,
            "endpoint_a_p_fixed_margins_exact": 0.040,
            "endpoint_a_global_fdr_family_status": "ELIGIBLE_M2_EXACT_FAMILY",
            "endpoint_a_cramers_v": 0.29,
            "endpoint_a_delta_alt_fraction": 1.00,
            "callability_p_analytic": None,
            "callability_testable": False,
            "callability_cramers_v": None,
            "core_partner_call_counts": {"R": 10, "A": 10},
        },
        {
            "endpoint_a_testable": True,
            "endpoint_a_p_analytic": 1e-20,
            "endpoint_a_p_fixed_margins_exact": 1e-20,
            "endpoint_a_global_fdr_family_status": "INELIGIBLE_M2_SCREEN",
            "endpoint_a_cramers_v": 1.0,
            "endpoint_a_delta_alt_fraction": 1.0,
            "callability_p_analytic": 0.01,
            "callability_testable": True,
            "callability_cramers_v": 0.8,
            "core_partner_call_counts": {"R": 8, "A": 8, "X": 4},
        },
    ]
    MODULE.apply_endpoint_a_bh_and_gate(rows)
    assert rows[0]["endpoint_a_q_global_bh"] == pytest.approx(0.002)
    assert rows[1]["endpoint_a_q_global_bh"] == pytest.approx(0.040)
    assert rows[0]["endpoint_a_q_global_by"] == pytest.approx(0.003)
    assert rows[1]["endpoint_a_q_global_by"] == pytest.approx(0.060)
    assert rows[2]["endpoint_a_q_global_bh"] is None
    assert rows[2]["endpoint_a_q_global_by"] is None
    assert rows[2]["endpoint_a_global_fdr_family_status"] == "INELIGIBLE_M2_SCREEN"
    assert rows[0]["endpoint_a_exact_by_discovery"] is True
    assert rows[0]["endpoint_a_effect_gate_pass"] is True
    assert rows[0]["endpoint_a_pre_candidate"] is True
    assert rows[0]["endpoint_a_pre_candidate_by_sensitivity"] is True
    assert rows[1]["endpoint_a_exact_bh_discovery"] is True
    assert rows[1]["endpoint_a_exact_by_discovery"] is False
    assert rows[1]["endpoint_a_effect_gate_pass"] is False
    assert rows[1]["endpoint_a_pre_candidate"] is False
    assert rows[2]["endpoint_a_pre_candidate"] is False
    assert rows[2]["callability_q_global_bh"] == pytest.approx(0.02)
    assert rows[0]["callability_gate_status"] == "PASS_ALL_CORE_READS_CALLABLE"
    assert rows[0]["callability_gate_pass"] is True
    assert rows[2]["callability_gate_status"] == "FAIL_DIFFERENTIAL_CALLABILITY_SIGNAL"
    assert rows[2]["callability_gate_pass"] is False


def test_summary_declares_execution_pass_without_scientific_upgrade() -> None:
    summary = MODULE.build_summary(
        entries=[],
        site_rows=[],
        pair_rows=[],
        dedup={},
        replication={},
        geometry={},
        stable_partner_audit={},
        conditional_sensitivity_audit={},
        joint_signature_global_fdr_audit={},
        site_pair_reconciliation={},
        full_scope=False,
    )

    assert summary["status"] == "EXECUTION_PASS"
    assert summary["pass"] is True
    assert summary["pass_semantics"] == (
        "execution_integrity_only_not_scientific_confirmation"
    )
    assert summary["n_multi_marker_molecular_haplotype_base_candidates"] == 0


def test_m2_axis_status_summary_conserves_pre_axis_not_evaluable_site() -> None:
    complete_statuses = {
        axis: {"status": "NOT_ALIGNED_EFFECT_BELOW_THRESHOLD_WITH_ADEQUATE_POWER"}
        for axis, _kind, _effect_suffix, _threshold in MODULE.M2_GATE.AXIS_SPECS
    }
    pre_axis_status = (
        "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM:"
        "observed=11:maximum=10"
    )
    rows = [
        {
            "sample": "S",
            "chrom": "chr1",
            "pos": 100,
            "m2_screen_eligibility_status": "INELIGIBLE_SCREEN_HP_CONFOUND",
            "m2_axis_statuses": complete_statuses,
        },
        {
            "sample": "S",
            "chrom": "chr1",
            "pos": 200,
            "m2_screen_eligibility_status": pre_axis_status,
            "m2_axis_statuses": {},
        },
    ]

    counts = MODULE.summarize_m2_axis_status_counts(rows)

    pre_axis_census_status = f"{MODULE.M2_PRE_AXIS_SUMMARY_STATUS}:{pre_axis_status}"
    assert set(counts) == {
        axis for axis, _kind, _effect_suffix, _threshold in MODULE.M2_GATE.AXIS_SPECS
    }
    assert all(sum(axis_counts.values()) == 2 for axis_counts in counts.values())
    assert all(axis_counts[pre_axis_census_status] == 1 for axis_counts in counts.values())


def test_m2_axis_status_summary_rejects_partial_axis_schema() -> None:
    rows = [
        {
            "sample": "S",
            "chrom": "chr1",
            "pos": 100,
            "m2_screen_eligibility_status": "NOT_EVALUABLE_M2_AXIS_INDETERMINATE",
            "m2_axis_statuses": {"hp_exact": {"status": "INDETERMINATE"}},
        }
    ]

    with pytest.raises(RuntimeError, match="axis status schema drift"):
        MODULE.summarize_m2_axis_status_counts(rows)


def test_m2_axis_status_summary_rejects_extra_axis_schema() -> None:
    complete_statuses = {
        axis: {"status": "NOT_ALIGNED_EFFECT_BELOW_THRESHOLD_WITH_ADEQUATE_POWER"}
        for axis, _kind, _effect_suffix, _threshold in MODULE.M2_GATE.AXIS_SPECS
    }
    complete_statuses["unexpected_axis"] = {"status": "INDETERMINATE"}
    rows = [
        {
            "sample": "S",
            "chrom": "chr1",
            "pos": 100,
            "m2_screen_eligibility_status": "NOT_EVALUABLE_M2_AXIS_INDETERMINATE",
            "m2_axis_statuses": complete_statuses,
        }
    ]

    with pytest.raises(RuntimeError, match="axis status schema drift"):
        MODULE.summarize_m2_axis_status_counts(rows)


def test_m2_axis_status_summary_rejects_post_axis_status_with_empty_map() -> None:
    rows = [
        {
            "sample": "HCC1954",
            "chrom": "chr5",
            "pos": 751076,
            "m2_screen_eligibility_status": (
                "NOT_EVALUABLE_M2_AXIS_INDETERMINATE:hp_exact"
            ),
            "m2_axis_statuses": {},
        }
    ]
    with pytest.raises(RuntimeError, match="allowlisted pre-axis"):
        MODULE.summarize_m2_axis_status_counts(rows)


def test_exact_by_effect_and_conditional_must_pass_on_the_same_formal_pair() -> None:
    labels = ["g1"] * 10 + ["g2"] * 10
    calls = ["R"] * 10 + ["A"] * 10
    reads = [recovered_read(index, label, call) for index, (label, call) in enumerate(zip(labels, calls))]
    row = MODULE._pair_result(
        {"sample": "S", "biological_id": "B"},
        MODULE.LIB.Variant("chr1", 100, "A", "C"),
        MODULE.LIB.Variant("chr1", 200, "A", "G"),
        reads,
        m2_eligible=True,
    )
    MODULE.apply_endpoint_a_bh_and_gate([row])
    assert row["endpoint_a_exact_by_discovery"] is True
    assert row["endpoint_a_effect_gate_pass"] is True
    assert row["endpoint_a_formal_pair_by_confirmed"] is False

    audit = MODULE.run_conditional_permutations([row], permutations=999)
    assert audit["n_exact_by_discoveries_submitted"] == 1
    assert audit["global_fdr_calibrated"] is False
    assert row["endpoint_a_permutable"] is True
    assert row["endpoint_a_p_conditional_perm"] <= 0.05
    assert row["endpoint_a_conditional_sensitivity_pass"] is True
    assert row["endpoint_a_formal_pair_by_confirmed"] is True


def test_differential_callability_blocks_formal_g1_pair() -> None:
    labels = ["g1"] * 10 + ["g2"] * 10
    calls = ["R"] * 10 + ["A"] * 6 + ["X"] * 4
    reads = [
        recovered_read(index, label, call)
        for index, (label, call) in enumerate(zip(labels, calls))
    ]
    row = MODULE._pair_result(
        {"sample": "S", "biological_id": "B"},
        MODULE.LIB.Variant("chr1", 100, "A", "C"),
        MODULE.LIB.Variant("chr1", 200, "A", "G"),
        reads,
        m2_eligible=True,
    )
    MODULE.apply_endpoint_a_bh_and_gate([row])
    assert row["callability_gate_pass"] is False
    assert row["callability_gate_status"] in {
        "FAIL_DIFFERENTIAL_CALLABILITY_SIGNAL",
        "NOT_IDENTIFIABLE_DIFFERENTIAL_CALLABILITY",
    }
    MODULE.run_conditional_permutations([row], permutations=999)
    assert row["endpoint_a_conditional_sensitivity_pass"] is True
    assert row["endpoint_a_formal_pair_by_confirmed"] is False


def test_conditional_permutation_nonexchangeable_is_not_identifiable() -> None:
    result = MODULE.LIB.conditional_label_permutation(
        ["g1"] * 6 + ["g2"] * 6,
        ["R"] * 6 + ["A"] * 6,
        ["HP1|PS=10|strand=+"] * 6 + ["HP2|PS=10|strand=+"] * 6,
        seed=11,
        permutations=999,
    )
    assert result["p_conditional_perm"] is None
    assert result["permutable"] is False
    assert result["status"] == "NOT_IDENTIFIABLE_NONEXCHANGEABLE"


def test_geometry_counts_every_dense_window_edge_not_only_new_right_edges() -> None:
    dense = MODULE.count_partner_geometry({"chr1": [100, 101, 102, 103]})
    assert dense == {
        "n_sites": 4,
        "n_sites_with_at_least_one_partner": 4,
        "n_unordered_pairs": 6,
    }
    mixed = MODULE.count_partner_geometry(
        {"chr1": [100, 200, 300, 5_100, 5_200], "chr2": [1]}
    )
    assert mixed == {
        "n_sites": 6,
        "n_sites_with_at_least_one_partner": 5,
        "n_unordered_pairs": 9,
    }
    assert sum(MODULE.EXPECTED_GEOMETRY["unordered_pairs_by_sample"].values()) == 864_255


def test_stable_focal_partner_oracle_hard_fails_a_missing_pair() -> None:
    oracle = {
        ("S", "chr1", 100): {
            ("chr1", 200, "A", "G"),
            ("chr1", 300, "C", "T"),
        }
    }
    rows = [
        {
            "sample": "S",
            "focal_chrom": "chr1",
            "focal_pos": 100,
            "partner_chrom": "chr1",
            "partner_pos": position,
            "partner_ref": ref,
            "partner_alt": alt,
        }
        for position, ref, alt in [(200, "A", "G"), (300, "C", "T")]
    ]
    audit = MODULE.validate_stable_partner_rows(oracle, rows)
    assert audit["n_expected_directed_pair_rows"] == 2
    assert audit["pass"] is True
    with pytest.raises(RuntimeError, match="enumeration mismatch"):
        MODULE.validate_stable_partner_rows(oracle, rows[:1])


def test_stable_site_and_assignment_denominators_must_match_both_directions(tmp_path: Path) -> None:
    site_path = tmp_path / "sites.tsv.gz"
    with gzip.open(site_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "sample\tchrom\tpos\tstable_null_multigroup\tref\talt\talt_readset_sha256\tscreen_contract\n"
        )
        handle.write(f"S\tchr1\t100\ttrue\tA\tC\tdigest-100\t{MODULE.SCREEN_CONTRACT}\n")
        handle.write(f"S\tchr1\t200\ttrue\tG\tT\tdigest-200\t{MODULE.SCREEN_CONTRACT}\n")
    assignments = {
        ("S", "chr1", 100): {
            "sample": "S",
            "chrom": "chr1",
            "pos": 100,
            "screen_contract": MODULE.SCREEN_CONTRACT,
        }
    }
    with pytest.raises(RuntimeError, match="missing read assignments"):
        MODULE.load_stable_site_rows(site_path, assignments, {"S"})


def test_candidate_screen_loader_requires_cluster_sizes_before_worker_launch(
    tmp_path: Path,
) -> None:
    site_path = tmp_path / "sites.tsv.gz"
    fieldnames = [
        "sample",
        "chrom",
        "pos",
        "stable_null_multigroup",
        "ref",
        "alt",
        "alt_readset_sha256",
        "screen_contract",
        "modal_assignment_ari_min",
        "hp_axis_confound",
        "technical_axis_confound",
        "residual_unexplained_multigroup",
        "phase_anchored_robust_epigenetic_candidate",
        *MODULE.M2_GATE.AXIS_SOURCE_FIELDS,
    ]
    with gzip.open(site_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
    with pytest.raises(RuntimeError, match="cluster_sizes"):
        MODULE.load_stable_site_rows(
            site_path,
            {},
            {"S"},
            require_candidate_screen_fields=True,
        )


def test_worker_task_strips_truth_ledger_and_assignment_posthoc_fields() -> None:
    key = ("S", "chr1", 100)
    assignment = {
        "sample": "S",
        "chrom": "chr1",
        "pos": 100,
        "region_dir": "/tmp/region",
        "screen_contract": MODULE.SCREEN_CONTRACT,
        "posthoc": {"truth_label": "FP", "ssnv_branch": "retained"},
    }
    site_row = {
        "region_dir": "/tmp/region",
        "alt_readset_sha256": "digest",
        "ref": "A",
        "alt": "C",
        "screen_contract": MODULE.SCREEN_CONTRACT,
        "cluster_sizes": '{"g1": 2, "g2": 1}',
        "truth_label": "FP",
        "ssnv_branch": "retained",
        "component_id": "component",
        "selected_group_id": "group",
    }
    tasks = MODULE.build_tasks(
        [{"sample": "S"}],
        {key: assignment},
        {key: site_row},
        top_markers=3,
    )
    assert "posthoc" not in tasks[0]["assignment"]
    assert tasks[0]["screen_row"] == {
        "region_dir": "/tmp/region",
        "alt_readset_sha256": "digest",
        "ref": "A",
        "alt": "C",
        "screen_contract": MODULE.SCREEN_CONTRACT,
        "cluster_sizes": '{"g1": 2, "g2": 1}',
    }


def topology_region(*, capped: bool = False, n_trees: int = 1) -> dict:
    return {
        "region": "chr1:100-300",
        "positions": [100, 200, 300],
        "lineages": [
            {
                "is_primary_lineage": True,
                "obs_col_coverage": {"100": [5, 5], "200": [5, 5], "300": [5, 5]},
                "capped": capped,
                "analysis_candidate_set_complete": True,
                "n_trees": n_trees,
                "trees": [
                    {
                        "edges": [["ROOT", "H_ARR"], ["H_ARR", "H_AAR"]],
                        "recurrence": [],
                    }
                ],
            }
        ],
    }


def test_topology_direct_scope_requires_retained_and_both_region_positions() -> None:
    region = topology_region()
    assert (
        MODULE.topology_direct_scope("retained", "retained", 100, 200, region)
        == "DIRECT_RETAINED_SHARED_REGION"
    )
    assert (
        MODULE.topology_direct_scope("max_snv_excluded", "retained", 100, 200, region)
        == "INDIRECT_MAX_SNV_EXCLUDED"
    )
    assert (
        MODULE.topology_direct_scope("retained", "positional_singleton", 100, 200, region)
        == "INDIRECT_NON_RETAINED"
    )
    assert (
        MODULE.topology_direct_scope("retained", "retained", 100, 400, region)
        == "OUTSIDE_REGION_POSITIONS"
    )
    assert MODULE.infer_topology_order(region, 100, 200) == "FOCAL_BEFORE_PARTNER"
    assert MODULE.infer_topology_order(topology_region(capped=True), 100, 200) == "NOT_INFERABLE_CAPPED"
    assert MODULE.infer_topology_order(topology_region(n_trees=2), 100, 200) == "NOT_INFERABLE_MULTI_TREE"


def minimal_site(sample: str, chrom: str, pos: int, digest: str = "unique") -> dict:
    return {
        "sample": sample,
        "biological_id": "HCC1395" if sample == "HCC1395_DORADO" else sample,
        "chrom": chrom,
        "pos": pos,
        "ref": "A",
        "alt": "C",
        "component_id": f"component-{pos}",
        "alt_readset_sha256": digest,
        "alt_readset_representative": False,
    }


def minimal_pair(sample: str, chrom: str, focal_pos: int, partner_pos: int) -> dict:
    return {
        "sample": sample,
        "focal_chrom": chrom,
        "focal_pos": focal_pos,
        "partner_pos": partner_pos,
    }


def test_oracle_extraction_schema_and_shared_readset_dedup() -> None:
    sites = [
        minimal_site("H2009", "chr3", 193_395_128),
        minimal_site("HCC1395_DORADO", "chr5", 750_311),
        minimal_site("COLO829", "chr4", 92_183_865),
        minimal_site("COLO829", "chr20", 52_761_564, "same-readset"),
        minimal_site("COLO829", "chr20", 52_761_565, "same-readset"),
    ]
    dedup = MODULE.apply_dedup_statistics(sites)
    pairs = [
        minimal_pair("HCC1395_DORADO", "chr5", 750_311, position)
        for position in [748_455, 750_429, 750_926, 751_297]
    ] + [minimal_pair("COLO829", "chr4", 92_183_865, 92_182_834)]
    payload = MODULE.extract_oracle_cases(sites, pairs)
    MODULE.validate_oracle_cases(payload, full_scope=True)
    assert payload["schema_name"].endswith(".oracle_cases")
    assert all(case["partner_window_oracle_pass"] for case in payload["focal_cases"])
    assert payload["shared_readset_case"]["same_alt_readset"] is True
    assert payload["shared_readset_case"]["one_alt_readset_representative"] is True
    assert dedup == {
        "n_dataset_sites": 5,
        "n_unique_biological_sites": 5,
        "n_unique_components": 5,
        "n_unique_alt_readsets": 4,
    }


def replication_row(
    sample: str,
    partner_alt: str = "A",
    *,
    relation: str = "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
    read_names: list[str] | None = None,
) -> dict:
    return {
        "sample": sample,
        "focal_chrom": "chr5",
        "focal_pos": 750_311,
        "focal_ref": "A",
        "focal_alt": "C",
        "partner_chrom": "chr5",
        "partner_pos": 750_429,
        "partner_ref": "C",
        "partner_alt": partner_alt,
        "endpoint_a_formal_pair_by_confirmed": True,
        "endpoint_a_effect_direction": "PARTNER_ALT_ENRICHED_IN_DOMINANT_METHYL_GROUP",
        "endpoint_a_cramers_v": 0.80,
        "endpoint_a_delta_alt_fraction": 0.90,
        "endpoint_b_relation_compatibility": relation,
        "cross_platform_exact_pair_present": False,
        "cross_platform_replication_status": "NOT_APPLICABLE",
        "_permutation_payload": {"core_raw_read_names": read_names or []},
    }


def test_cross_platform_replication_requires_exact_identity_and_reports_read_overlap() -> None:
    hcc = replication_row("HCC1395", read_names=["ont-only", "shared"])
    dorado = replication_row("HCC1395_DORADO", read_names=["dorado-only", "shared"])
    nonmatching = replication_row("HCC1395_DORADO", "G", read_names=["other"])
    summary = MODULE.apply_cross_platform_replication([hcc, dorado, nonmatching])
    assert summary["n_exact_pairs_present_both"] == 1
    assert summary["n_exact_pairs_conditional_by_effect_direction_replicated"] == 1
    assert summary["n_exact_pairs_with_core_raw_read_name_overlap"] == 1
    assert summary["biological_n"] == 1
    assert summary["independent_technical_library_claim_allowed"] is False
    assert hcc["cross_platform_core_read_name_intersection_n"] == 1
    assert hcc["cross_platform_core_read_name_union_n"] == 3
    assert hcc["cross_platform_core_read_name_jaccard"] == pytest.approx(1.0 / 3.0)
    assert hcc["cross_platform_core_read_name_overlap_present"] is True
    assert hcc["cross_platform_biological_n"] == 1
    assert hcc["cross_platform_four_state_compatible"] is True
    assert hcc["cross_platform_replication_status"] == (
        "SAME_MOLECULE_REPROCESSING_CONDITIONAL_BY_CONCORDANCE"
    )
    assert dorado["cross_platform_replication_status"] == (
        "SAME_MOLECULE_REPROCESSING_CONDITIONAL_BY_CONCORDANCE"
    )
    assert nonmatching["cross_platform_replication_status"] == "ONE_PLATFORM_ONLY"
    assert "_permutation_payload" not in hcc
    assert "_permutation_payload" not in dorado


def test_cross_platform_incompatible_four_state_is_not_mislabelled_compatible() -> None:
    relation = "INCOMPATIBLE_OR_COMPLEX_UNDER_FIXED_ERROR_MODEL"
    hcc = replication_row("HCC1395", relation=relation, read_names=["hcc"])
    dorado = replication_row("HCC1395_DORADO", relation=relation, read_names=["dorado"])
    MODULE.apply_cross_platform_replication([hcc, dorado])
    assert hcc["cross_platform_four_state_compatible"] is False
    assert dorado["cross_platform_four_state_compatible"] is False
    assert hcc["cross_platform_conditional_by_effect_direction_replication_pass"] is True
    assert dorado["cross_platform_conditional_by_effect_direction_replication_pass"] is True


def test_joint_signature_uses_complete_reads_and_spaced_markers_control_formal_candidate() -> None:
    partners = [
        MODULE.LIB.Variant("chr1", 1_100, "A", "G"),
        MODULE.LIB.Variant("chr1", 1_120, "C", "T"),
    ]
    reads = []
    for index in range(24):
        first_group = index < 12
        reads.append(
            MODULE.RecoveredRead(
                read_id=f"joint-{index}",
                read_name=f"joint-{index}",
                expected_focal_call="A",
                focal_call="A",
                latest_hp="1-1" if index % 2 == 0 else "2-1",
                latest_ps=101,
                strand="+",
                label="g1" if first_group else "g2",
                partner_calls={
                    1_100: "R" if first_group else "A",
                    1_120: "R" if first_group else "A",
                },
            )
        )
    reads.append(
        MODULE.RecoveredRead(
            read_id="incomplete",
            read_name="incomplete",
            expected_focal_call="A",
            focal_call="A",
            latest_hp="1-1",
            latest_ps=101,
            strand="+",
            label="g1",
            partner_calls={1_100: "R", 1_120: "X"},
        )
    )
    signature = MODULE._signature_result(reads, partners, seed=29, permutations=999)
    assert signature["joint_signature_n_complete_reads"] == 24
    assert signature["joint_signature_testable"] is True
    assert signature["joint_signature_permutable"] is True
    assert signature["joint_signature_sensitivity_pass"] is True
    assert signature["joint_signature_postselection_fdr_calibrated"] is False
    assert signature["joint_signature_n_complete_marker_effect_supported"] == 2

    site = {
        "sample": "S",
        "chrom": "chr1",
        "pos": 1_000,
        "n_partner_markers": 3,
        "m2_screen_eligible": True,
        "top_marker_positions": [1_100, 1_110, 1_120],
        "joint_signature_complete_marker_effect_supported_positions": [1_100, 1_120],
        "joint_signature_sensitivity_pass": signature["joint_signature_sensitivity_pass"],
        "joint_signature_global_by_discovery": True,
    }
    pairs = [
        {
            "sample": "S",
            "focal_chrom": "chr1",
            "focal_pos": 1_000,
            "partner_pos": position,
            "endpoint_a_testable": True,
            "endpoint_a_global_fdr_family_status": "ELIGIBLE_M2_EXACT_FAMILY",
            "endpoint_a_p_fixed_margins_exact": 0.001,
            "endpoint_a_exact_bh_discovery": True,
            "endpoint_a_exact_by_discovery": True,
            "endpoint_a_formal_pair_by_confirmed": True,
            "endpoint_a_permutable": True,
        }
        for position in [1_100, 1_110, 1_120]
    ]
    MODULE.finalize_site_inference([site], pairs)
    assert site["n_pair_by_confirmed"] == 3
    assert site["n_spatially_separated_pair_by_20bp"] == 2
    assert site["spatially_separated_pair_by_positions_20bp"] == [1_100, 1_120]
    assert site["n_top_marker_pair_by_confirmed"] == 3
    assert site["multi_marker_molecular_haplotype_base_candidate"] is True
    assert site[
        "spatially_separated_same_complete_read_effect_supported_top_pair_by_positions_20bp"
    ] == [1_100, 1_120]
    assert site["spaced_marker_20bp_contract"] == (
        "genomic_distance_at_least_20bp_only_not_statistical_independence"
    )


def test_joint_signature_cannot_be_driven_by_only_one_marker_on_complete_reads() -> None:
    partners = [
        MODULE.LIB.Variant("chr1", 1_100, "A", "G"),
        MODULE.LIB.Variant("chr1", 1_120, "C", "T"),
    ]
    reads = [
        MODULE.RecoveredRead(
            read_id=f"disjoint-{index}",
            read_name=f"disjoint-{index}",
            expected_focal_call="A",
            focal_call="A",
            latest_hp="1-1" if index % 2 == 0 else "2-1",
            latest_ps=101,
            strand="+",
            label="g1" if index < 10 else "g2",
            partner_calls={
                1_100: "R" if index < 10 else "A",
                1_120: "R",
            },
        )
        for index in range(20)
    ]
    signature = MODULE._signature_result(reads, partners, seed=29, permutations=999)
    assert signature["joint_signature_sensitivity_pass"] is True
    assert signature["joint_signature_n_complete_marker_effect_supported"] == 1
    assert signature["joint_signature_complete_marker_effect_supported_positions"] == [1_100]

    site = {
        "sample": "S",
        "chrom": "chr1",
        "pos": 1_000,
        "n_partner_markers": 2,
        "m2_screen_eligible": True,
        "top_marker_positions": [1_100, 1_120],
        "joint_signature_complete_marker_effect_supported_positions": [1_100],
        "joint_signature_sensitivity_pass": True,
    }
    pairs = [
        {
            "sample": "S",
            "focal_chrom": "chr1",
            "focal_pos": 1_000,
            "partner_pos": position,
            "endpoint_a_testable": True,
            "endpoint_a_global_fdr_family_status": "ELIGIBLE_M2_EXACT_FAMILY",
            "endpoint_a_p_fixed_margins_exact": 0.001,
            "endpoint_a_exact_bh_discovery": True,
            "endpoint_a_exact_by_discovery": True,
            "endpoint_a_formal_pair_by_confirmed": True,
            "endpoint_a_permutable": True,
        }
        for position in [1_100, 1_120]
    ]
    MODULE.finalize_site_inference([site], pairs)
    assert site["n_spatially_separated_pair_by_20bp"] == 2
    assert site["multi_marker_molecular_haplotype_base_candidate"] is False


def test_formal_claim_names_do_not_call_spaced_markers_independent_or_subclonal() -> None:
    formal_claim_fields = {
        "endpoint_a_formal_pair_by_confirmed",
        "n_spatially_separated_pair_by_20bp",
        "spatially_separated_pair_by_positions_20bp",
        "spaced_marker_20bp_contract",
        "multi_marker_molecular_haplotype_base_candidate",
    }
    assert formal_claim_fields.issubset(set(MODULE.PAIR_COLUMNS) | set(MODULE.SITE_COLUMNS))
    assert all("independent" not in field.lower() for field in formal_claim_fields)
    assert all("subclone" not in field.lower() for field in formal_claim_fields)
    assert MODULE.SPACED_MARKER_MIN_SEPARATION_BP == 20


def test_m2_categorical_level_counts_uses_stable_assignment_fields() -> None:
    assignment = {
        "core_reads": [
            {"label": "g1", "latest_hp": "1-1", "strand": "+"},
            {"label": "g1", "latest_hp": "1-2", "strand": "-"},
            {"label": "g2", "latest_hp": "2-1", "strand": "+"},
        ]
    }
    screen_row = {"cluster_sizes": '{"g1": 2, "g2": 1}'}

    assert MODULE.m2_categorical_level_counts(assignment, screen_row) == {
        "hp_exact": 3,
        "hp_family": 2,
        "strand": 2,
    }


def test_site_task_ranks_top_markers_by_coverage_not_p_value(monkeypatch: pytest.MonkeyPatch) -> None:
    focal = MODULE.LIB.Variant("chr1", 100, "A", "C")
    partners = [
        MODULE.LIB.Variant("chr1", 200, "A", "G"),
        MODULE.LIB.Variant("chr1", 300, "C", "T"),
        MODULE.LIB.Variant("chr1", 400, "G", "A"),
    ]
    reads = []
    for index in range(20):
        reads.append(
            MODULE.RecoveredRead(
                read_id=f"r{index}",
                read_name=f"r{index}",
                expected_focal_call="A",
                focal_call="A",
                latest_hp="1-1" if index < 10 else "2-1",
                latest_ps=10,
                strand="+",
                label="g1" if index < 10 else "g2",
                partner_calls={
                    200: "R" if index < 10 else "A",
                    300: "R" if index < 5 else "A" if index < 10 else "X",
                    400: "X",
                },
                readset_identity=f"r{index}|chr1|0|1200|+",
            )
        )
    digest = MODULE.focal_alt_readset_digest(reads)
    monkeypatch.setattr(MODULE, "get_vcf", lambda *_: object())
    monkeypatch.setattr(MODULE, "discover_focal_and_partners", lambda *_: (focal, partners))

    def recover_with_audit(*args):
        populate_empty_raw_identity_audit(args[-1], len(reads))
        return reads

    monkeypatch.setattr(MODULE, "recover_site_reads", recover_with_audit)
    monkeypatch.setattr(
        MODULE,
        "m2_categorical_level_counts",
        lambda *_: {"hp_exact": 2, "hp_family": 2, "strand": 1},
    )
    result = MODULE.analyze_site_task(
        {
            "entry": {
                "sample": "S",
                "biological_id": "B",
                "all_ssnv_vcf": {"path": "/unused.vcf.gz"},
                "all_ssnv_vcf_index": {"path": "/unused.vcf.gz.csi"},
            },
            "assignment": {"sample": "S", "chrom": "chr1", "pos": 100, "region_dir": "/unused"},
            "screen_row": {
                "ref": "A",
                "alt": "C",
                "truth_label": "FP",
                "ssnv_branch": "retained",
                "component_id": "component",
                "selected_group_id": "group",
                "alt_readset_sha256": digest,
                "screen_contract": MODULE.SCREEN_CONTRACT,
            },
            "top_markers": 2,
        }
    )
    assert result["site"]["top_marker_positions"] == [200, 300]
    assert set(result["site"]) == set(MODULE.SITE_COLUMNS)
    assert result["raw_identity_duplicates"] == []
    ranks = {row["partner_pos"]: row["top_coverage_rank"] for row in result["pairs"]}
    assert ranks == {200: 1, 300: 2, 400: None}


def test_output_directory_is_never_reused(tmp_path: Path) -> None:
    output = tmp_path / "new-output"
    MODULE.create_output_dir(output)
    assert output.is_dir()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.create_output_dir(output)


def test_manifest_frozen_identity_policy_is_explicit_and_fail_closed(
    tmp_path: Path,
) -> None:
    artifact_path = tmp_path / "artifact.bin"
    artifact_path.write_bytes(b"frozen-input\n")
    metadata_record = {
        "path": str(artifact_path),
        "size_bytes": artifact_path.stat().st_size,
        "mtime_ns": artifact_path.stat().st_mtime_ns,
    }

    with pytest.raises(RuntimeError, match="required SHA-256"):
        MODULE.verify_manifest_frozen_field("all_ssnv_vcf", metadata_record)

    large_identity = MODULE.verify_manifest_frozen_field(
        "raw_alignment", metadata_record
    )
    assert large_identity["identity_mode"] == (
        "explicit_large_artifact_size_mtime_path_v1"
    )
    assert large_identity["sha256"] is None

    hashed_record = {**metadata_record, "sha256": MODULE.sha256(artifact_path)}
    hashed_identity = MODULE.verify_manifest_frozen_field(
        "all_ssnv_vcf", hashed_record
    )
    assert hashed_identity["identity_mode"] == "sha256_size_path_v1"

    with pytest.raises(RuntimeError, match="lacks an identity policy"):
        MODULE.verify_manifest_frozen_field("unexpected_role", hashed_record)


def test_formal_full_scope_cli_is_directly_executable_and_uses_production_defaults(
    tmp_path: Path,
) -> None:
    completed = subprocess.run(
        [sys.executable, "-I", str(SCRIPT), "--help"],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert "--exact-state-space-ceiling" in completed.stdout
    assert "--primary-artifact-audit-pre" in completed.stdout
    argv = [
        "--manifest",
        str(tmp_path / "manifest.json"),
        "--assignments",
        str(tmp_path / "assignments.jsonl.gz"),
        "--sites",
        str(tmp_path / "sites.tsv.gz"),
        "--primary-artifact-audit-pre",
        str(tmp_path / "stable_primary_artifact_audit.v1_pre_downstream.json"),
        "--independent-m2-audit",
        str(tmp_path / "independent_m2.json"),
        "--preflight-receipt",
        str(tmp_path / "cooccurrence_preflight.json"),
        "--intersubmod-root",
        str(tmp_path / "runs"),
        "--output-dir",
        str(tmp_path / "output"),
    ]
    args = MODULE.parse_args(argv)
    assert args.sample == []
    assert args.allow_partial_scope is False
    assert args.permutations == 999
    assert args.exact_state_space_ceiling == MODULE.EXACT_STATE_SPACE_CEILING
    assert args.max_pending_chunks >= args.workers
    expected_command = MODULE.canonical_producer_command(
        manifest=Path(argv[1]),
        assignments=Path(argv[3]),
        sites=Path(argv[5]),
        primary_artifact_audit_pre=Path(argv[7]),
        independent_m2_audit=Path(argv[9]),
        preflight_receipt=Path(argv[11]),
        intersubmod_root=Path(argv[13]),
        output_dir=Path(argv[15]),
    )
    expected_cache = tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME
    assert expected_command[:5] == [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={expected_cache}",
        str(SCRIPT.resolve()),
    ]
    with pytest.raises(ValueError, match=r"--permutations must equal canonical value 999"):
        MODULE.main([*argv, "--permutations", "998"])
    with pytest.raises(ValueError, match=r"--top-markers must equal canonical value 3"):
        MODULE.main([*argv, "--top-markers", "2"])
    with pytest.raises(ValueError, match=r"--exact-state-space-ceiling must equal canonical value"):
        MODULE.main([*argv, "--exact-state-space-ceiling", "10"])


def test_formal_release_command_rejects_imported_main_and_requires_proc_command(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    expected = [sys.executable, "-I", str(SCRIPT.resolve()), "--help"]
    formal_authority = {"authority_id": MODULE.SOURCE_AUTHORITY.AUTHORITY_ID}
    with pytest.raises(ValueError, match="direct-CLI only"):
        MODULE.resolve_release_command(
            argv=["--help"],
            script_path=SCRIPT,
            expected_command=expected,
            source_authority=formal_authority,
            role="producer",
        )

    monkeypatch.setattr(MODULE, "observed_process_command", lambda: list(expected))
    assert MODULE.resolve_release_command(
        argv=None,
        script_path=SCRIPT,
        expected_command=expected,
        source_authority=formal_authority,
        role="producer",
    ) == expected

    monkeypatch.setattr(
        MODULE,
        "observed_process_command",
        lambda: [sys.executable, str(SCRIPT.resolve()), "--help"],
    )
    with pytest.raises(ValueError, match="canonical Task-B command"):
        MODULE.resolve_release_command(
            argv=None,
            script_path=SCRIPT,
            expected_command=expected,
            source_authority=formal_authority,
            role="producer",
        )


def test_formal_python_cache_contract_requires_exact_nonsymlink_mode0700(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "cooccurrence"
    cache_root = tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME
    cache_root.mkdir(mode=0o700)
    monkeypatch.setattr(MODULE.sys, "pycache_prefix", str(cache_root))
    monkeypatch.setattr(
        MODULE.sys,
        "_xoptions",
        {"pycache_prefix": str(cache_root)},
    )
    assert MODULE.capture_python_cache_contract(output_dir) == {
        "path": str(cache_root),
        "mode": "0o700",
        "sys_pycache_prefix": str(cache_root),
        "xoption_pycache_prefix": str(cache_root),
        "contract": "canonical_nonsymlink_mode0700_python_xoption_v1",
    }

    cache_root.chmod(0o755)
    with pytest.raises(RuntimeError, match="mode-0700"):
        MODULE.capture_python_cache_contract(output_dir)

    cache_root.chmod(0o700)
    monkeypatch.setattr(MODULE.sys, "pycache_prefix", str(tmp_path / "wrong"))
    with pytest.raises(RuntimeError, match="canonical Task-B path"):
        MODULE.capture_python_cache_contract(output_dir)


def make_valid_preflight_gate_fixture(
    tmp_path: Path,
) -> tuple[Path, dict[str, Path], dict[str, dict]]:
    inputs = {
        name: tmp_path / filename
        for name, filename in (
            ("manifest", "manifest.json"),
            ("assignments", "assignments.jsonl.gz"),
            ("sites", "sites.tsv.gz"),
            ("primary_artifact_audit_pre", "primary_pre.json"),
            ("independent_m2_audit", "independent_m2.json"),
        )
    }
    for name, path in inputs.items():
        path.write_text(f"{name}\n", encoding="utf-8")
    current_code = MODULE.capture_code_source_identity()
    preflight_sources = {
        "preflight": MODULE.artifact(
            SCRIPT_DIR / "audit_cooccurrence_task_contract_preflight.py"
        ),
        **current_code,
    }
    modes = {role: "0o444" for role in preflight_sources}
    intersubmod_root = tmp_path / "runs"
    intersubmod_root.mkdir()
    path = tmp_path / "preflight.json"
    input_identity = {
        name: MODULE.artifact(input_path) for name, input_path in inputs.items()
    }
    runtime_identity = MODULE.capture_runtime_environment_identity()
    payload = {
        "schema_name": MODULE.PREFLIGHT_SCHEMA_NAME,
        "schema_version": MODULE.PREFLIGHT_SCHEMA_VERSION,
        "task_type": "B_comprehensive_validation",
        "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
        "source_authority": (
            MODULE.SOURCE_AUTHORITY.validate_release_source_authority()
        ),
        "started_at_utc": "2026-07-15T00:00:00+00:00",
        "finished_at_utc": "2026-07-15T00:01:00+00:00",
        "created_at_utc": "2026-07-15T00:01:00+00:00",
        "command": MODULE.canonical_preflight_command(
            manifest=inputs["manifest"],
            assignments=inputs["assignments"],
            sites=inputs["sites"],
            intersubmod_root=intersubmod_root,
            independent_m2_audit=inputs["independent_m2_audit"],
            primary_artifact_audit_pre=inputs["primary_artifact_audit_pre"],
            output=path,
        ),
        "inputs": input_identity,
        "input_lock": {
            "identity_before": input_identity,
            "identity_after": input_identity,
            "all_primary_inputs_unchanged": True,
        },
        "runtime_environment_lock": {
            "identity_before": runtime_identity,
            "identity_after": runtime_identity,
            "direct_runtime_unchanged": True,
        },
        "parameters": dict(MODULE.PREFLIGHT_CANONICAL_PARAMETERS),
        "code": {
            "source_identity_before": preflight_sources,
            "source_identity_after": preflight_sources,
            "source_modes_before": modes,
            "source_modes_after": modes,
        },
        "observed": {
            "task_count": MODULE.EXPECTED_STABLE_SITE_TASKS,
            "runtime_statistical_api_probe": MODULE.runtime_statistical_api_probe(),
            "raw_bam_identity_recovery": {
                "aggregate": {"site_tasks": MODULE.EXPECTED_STABLE_SITE_TASKS},
                "all_result_rows_passed_invariant_validation": True,
                "missing_projection_policy": MODULE.RAW_IDENTITY_MISSING_POLICY,
                "conflicting_analysis_payload_policy": (
                    MODULE.RAW_IDENTITY_CONFLICT_POLICY
                ),
                "failure_counts_materialized": False,
                "analysis_scope_identity_contract": (
                    MODULE.ANALYSIS_SCOPE_IDENTITY_CONTRACT
                ),
            },
        },
        "checks": {
            key: True for key in MODULE.PREFLIGHT_CANONICAL_CHECK_KEYS
        },
        "pass": True,
        "pass_semantics": MODULE.PREFLIGHT_PASS_SEMANTICS,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    path.chmod(0o444)
    return path, inputs, current_code


def test_preflight_gate_binds_inputs_and_exact_source_roles(tmp_path: Path) -> None:
    path, inputs, current_code = make_valid_preflight_gate_fixture(tmp_path)
    observed = MODULE.validate_preflight_gate(
        path,
        manifest=inputs["manifest"],
        assignments=inputs["assignments"],
        sites=inputs["sites"],
        intersubmod_root=tmp_path / "runs",
        independent_m2_audit=inputs["independent_m2_audit"],
        primary_artifact_audit_pre=inputs["primary_artifact_audit_pre"],
        current_code_identity=current_code,
        current_runtime_identity=MODULE.capture_runtime_environment_identity(),
        current_runtime_probe=MODULE.runtime_statistical_api_probe(),
        producer_started_at_utc="2026-07-15T00:02:00+00:00",
    )
    assert observed["pass"] is True


def test_preflight_gate_rejects_source_role_superset(tmp_path: Path) -> None:
    path, inputs, current_code = make_valid_preflight_gate_fixture(tmp_path)
    payload = json.loads(path.read_text())
    for field_name in ("source_identity_before", "source_identity_after"):
        payload["code"][field_name]["unexpected"] = MODULE.artifact(
            inputs["manifest"]
        )
    for field_name in ("source_modes_before", "source_modes_after"):
        payload["code"][field_name]["unexpected"] = "0o444"
    path.chmod(0o644)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    path.chmod(0o444)
    with pytest.raises(RuntimeError, match="source identity map drifted"):
        MODULE.validate_preflight_gate(
            path,
            manifest=inputs["manifest"],
            assignments=inputs["assignments"],
            sites=inputs["sites"],
            intersubmod_root=tmp_path / "runs",
            independent_m2_audit=inputs["independent_m2_audit"],
            primary_artifact_audit_pre=inputs["primary_artifact_audit_pre"],
            current_code_identity=current_code,
            current_runtime_identity=MODULE.capture_runtime_environment_identity(),
            current_runtime_probe=MODULE.runtime_statistical_api_probe(),
            producer_started_at_utc="2026-07-15T00:02:00+00:00",
        )


def test_run_receipt_binds_pre_audit_before_completion_timestamp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = tmp_path / "manifest.json"
    assignments = tmp_path / "assignments.jsonl.gz"
    sites = tmp_path / "sites.tsv.gz"
    pre_audit = tmp_path / "stable_primary_artifact_audit.v1_pre_downstream.json"
    independent_m2 = tmp_path / "independent_m2.json"
    preflight = tmp_path / "cooccurrence_preflight.json"
    for path in (manifest, assignments, sites, pre_audit, independent_m2, preflight):
        path.write_text(f"{path.name}\n", encoding="utf-8")
    runs = tmp_path / "runs"
    runs.mkdir()
    output = tmp_path / "output"

    monkeypatch.setattr(
        MODULE,
        "load_manifest",
        lambda *_args, **_kwargs: ({"samples": []}, [], {"samples": []}),
    )
    monkeypatch.setattr(
        MODULE,
        "audit_partner_geometry",
        lambda *_args, **_kwargs: {"aggregate": {"n_unordered_pairs": 0}},
    )
    monkeypatch.setattr(MODULE, "load_assignments", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(MODULE, "load_stable_site_rows", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(MODULE, "build_stable_partner_oracle", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(MODULE, "build_tasks", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(MODULE, "bounded_process_map", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(
        MODULE,
        "validate_stable_partner_rows",
        lambda *_args, **_kwargs: {"n_expected_directed_pair_rows": 0},
    )
    monkeypatch.setattr(
        MODULE,
        "run_conditional_permutations",
        lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(
        MODULE,
        "finalize_site_inference",
        lambda *_args, **_kwargs: {
            "n_m2_global_fdr_eligible_sites": 0,
            "n_m2_exact_global_fdr_family_pairs": 0,
        },
    )
    monkeypatch.setattr(MODULE, "attach_posthoc_annotations", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(MODULE, "apply_dedup_statistics", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(MODULE, "apply_cross_platform_replication", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(MODULE, "extract_oracle_cases", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(MODULE, "validate_oracle_cases", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        MODULE,
        "validate_preflight_gate",
        lambda *_args, **_kwargs: {
            "schema_name": MODULE.PREFLIGHT_SCHEMA_NAME,
            "schema_version": MODULE.PREFLIGHT_SCHEMA_VERSION,
            "observed": {"task_count": MODULE.EXPECTED_STABLE_SITE_TASKS},
            "pass": True,
            "raw_identity_release_contract": MODULE.RAW_IDENTITY_RELEASE_CONTRACT,
        },
    )
    monkeypatch.setattr(
        MODULE,
        "build_summary",
        lambda *_args, **_kwargs: {
            "pass": True,
            "frozen_input_identity_policy": {"artifacts": []},
                "m2_axis_statistic_provenance": {
                    "downstream_raw_read_axis_statistic_recomputed": False
                },
                "raw_bam_identity_recovery_audit": {
                    "n_expected_projection_occurrences": 0,
                    "n_duplicate_projection_occurrences_collapsed": 0,
                    "n_duplicate_extra_record_occurrences_collapsed": 0,
                    "all_site_results_passed_invariant_validation": True,
                },
            },
        )

    timestamps = [
        "2026-07-15T00:00:00+00:00",
        "2026-07-15T00:01:00+00:00",
    ]
    clock_calls: list[str] = []

    def fake_now_utc() -> str:
        value = timestamps[len(clock_calls)]
        clock_calls.append(value)
        return value

    original_artifact = MODULE.artifact
    artifact_clock_counts: list[int] = []

    def tracked_artifact(path: Path, include_hash: bool = True) -> dict:
        artifact_clock_counts.append(len(clock_calls))
        return original_artifact(path, include_hash=include_hash)

    monkeypatch.setattr(MODULE, "now_utc", fake_now_utc)
    monkeypatch.setattr(MODULE, "artifact", tracked_artifact)
    monkeypatch.setattr(
        MODULE,
        "capture_code_source_modes",
        lambda: {role: "0o444" for role in MODULE.code_source_paths()},
    )
    test_cache = {
        "path": str(tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME),
        "mode": "0o700",
        "sys_pycache_prefix": str(
            tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME
        ),
        "xoption_pycache_prefix": str(
            tmp_path / MODULE.CANONICAL_PYTHON_CACHE_DIRNAME
        ),
        "contract": "canonical_nonsymlink_mode0700_python_xoption_v1",
    }
    monkeypatch.setattr(
        MODULE,
        "resolve_release_command",
        lambda **kwargs: list(kwargs["expected_command"]),
    )
    monkeypatch.setattr(
        MODULE,
        "capture_python_cache_contract",
        lambda _output_dir: dict(test_cache),
    )

    argv = [
        "--manifest",
        str(manifest),
        "--assignments",
        str(assignments),
        "--sites",
        str(sites),
        "--primary-artifact-audit-pre",
        str(pre_audit),
        "--independent-m2-audit",
        str(independent_m2),
        "--preflight-receipt",
        str(preflight),
        "--intersubmod-root",
        str(runs),
        "--output-dir",
        str(output),
        "--workers",
        "40",
        "--chunk-size",
        "8",
        "--max-pending-chunks",
        "80",
        "--permutations",
        "999",
        "--top-markers",
        "3",
        "--exact-state-space-ceiling",
        "250000",
    ]
    assert MODULE.main(argv) == 0

    receipt = json.loads((output / MODULE.RECEIPT_OUTPUT_NAME).read_text(encoding="utf-8"))
    assert receipt["inputs"]["primary_artifact_audit_pre"] == original_artifact(pre_audit)
    assert set(receipt["inputs"]["primary_artifact_audit_pre"]) == {
        "path",
        "size_bytes",
        "sha256",
    }
    assert receipt["started_at_utc"] <= receipt["finished_at_utc"]
    assert receipt["created_at_utc"] == receipt["finished_at_utc"]
    assert clock_calls == timestamps
    assert artifact_clock_counts
    assert receipt["code"]["m2_screen_gate"] == original_artifact(
        Path(MODULE.M2_GATE.__file__).resolve()
    )
    assert set(artifact_clock_counts) == {1}
    assert (output / MODULE.RAW_IDENTITY_DUPLICATE_OUTPUT_NAME).is_file()
    assert receipt["source_lock"]["all_sources_read_only_and_unchanged"] is True
