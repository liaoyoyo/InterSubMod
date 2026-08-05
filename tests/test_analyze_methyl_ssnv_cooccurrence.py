from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "scripts"
    / "analyze_methyl_ssnv_cooccurrence.py"
)
SPEC = importlib.util.spec_from_file_location("test_analyze_methyl_ssnv_cooccurrence_v2", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_categorical_association_uses_scipy_tuple_compatible_p_value() -> None:
    labels = ["g1"] * 10 + ["g2"] * 10 + ["g3"] * 10
    categories = (
        ["R"] * 8
        + ["A"] * 2
        + ["R"] * 5
        + ["A"] * 5
        + ["R"] * 2
        + ["A"] * 8
    )

    observed = MODULE.categorical_association(labels, categories)
    expected = MODULE.chi2_contingency(observed["table"], correction=False)

    assert observed["testable"] is True
    assert observed["analytic_test"] == "pearson_chi_square"
    assert observed["p_analytic"] == pytest.approx(float(expected[1]))


def recovered_reads(partner_calls: dict[int, list[str]]) -> list[MODULE.RecoveredRead]:
    count = len(next(iter(partner_calls.values())))
    assert all(len(calls) == count for calls in partner_calls.values())
    reads = []
    for index in range(count):
        reads.append(
            MODULE.RecoveredRead(
                read_id=f"read-{index}",
                read_name=f"raw-{index}",
                expected_focal_call="A",
                focal_call="A",
                latest_hp="1-1" if index % 2 == 0 else "2-1",
                latest_ps=101,
                strand="+",
                label="g1" if index < count // 2 else "g2",
                partner_calls={position: calls[index] for position, calls in partner_calls.items()},
                readset_identity=f"raw-{index}|chr1|0|1200|+",
            )
        )
    return reads


def test_global_exact_family_excludes_screen_confounded_pairs() -> None:
    calls = ["R"] * 10 + ["A"] * 10
    reads = recovered_reads({200: calls})
    focal = MODULE.LIB.Variant("chr1", 100, "A", "C")
    partner = MODULE.LIB.Variant("chr1", 200, "A", "G")
    eligible = MODULE._pair_result(
        {"sample": "eligible", "biological_id": "B"},
        focal,
        partner,
        reads,
        m2_eligible=True,
    )
    confounded = MODULE._pair_result(
        {"sample": "confounded", "biological_id": "C"},
        focal,
        partner,
        reads,
        m2_eligible=False,
    )
    MODULE.apply_endpoint_a_bh_and_gate([eligible, confounded])
    assert eligible["endpoint_a_q_global_by"] is not None
    assert eligible["endpoint_a_exact_by_discovery"] is True
    assert confounded["endpoint_a_p_fixed_margins_exact"] == pytest.approx(
        eligible["endpoint_a_p_fixed_margins_exact"]
    )
    assert confounded["endpoint_a_q_global_bh"] is None
    assert confounded["endpoint_a_q_global_by"] is None
    assert confounded["endpoint_a_global_fdr_family_status"] == "INELIGIBLE_M2_SCREEN"

    audit = MODULE.run_conditional_permutations([eligible, confounded], permutations=99)
    assert audit["n_exact_by_discoveries_submitted"] == 1
    assert eligible["endpoint_a_conditional_status"] == "PERMUTABLE"
    assert eligible["endpoint_a_formal_pair_by_confirmed"] is True
    assert confounded["endpoint_a_conditional_status"] == "NOT_RUN_NOT_EXACT_BY_DISCOVERY"


def test_analytic_p_cannot_enter_formal_global_family() -> None:
    row = {
        "endpoint_a_testable": True,
        "endpoint_a_p_analytic": 1e-12,
        "endpoint_a_cramers_v": 1.0,
        "endpoint_a_delta_alt_fraction": 1.0,
        "callability_p_analytic": None,
    }
    MODULE.apply_endpoint_a_bh_and_gate([row])
    assert row["endpoint_a_q_global_bh"] is None
    assert row["endpoint_a_q_global_by"] is None
    assert row["endpoint_a_exact_bh_discovery"] is False
    assert row["endpoint_a_exact_by_discovery"] is False


def test_pair_conditional_payload_uses_hp_family_ps_and_strand() -> None:
    reads = recovered_reads({200: ["R"] * 10 + ["A"] * 10})
    row = MODULE._pair_result(
        {"sample": "S", "biological_id": "B"},
        MODULE.LIB.Variant("chr1", 100, "A", "C"),
        MODULE.LIB.Variant("chr1", 200, "A", "G"),
        reads,
        m2_eligible=True,
    )
    strata = row["_permutation_payload"]["strata"]
    assert row["endpoint_a_conditional_strata"] == "latest_HP_family_x_PS_x_strand"
    assert all("PS=101" in stratum and "strand=+" in stratum for stratum in strata)


def test_different_reads_cannot_be_stitched_into_joint_signature() -> None:
    reads = recovered_reads(
        {
            200: ["R"] * 6 + ["X"] * 6,
            300: ["X"] * 6 + ["A"] * 6,
        }
    )
    result = MODULE._signature_result(
        reads,
        [
            MODULE.LIB.Variant("chr1", 200, "A", "G"),
            MODULE.LIB.Variant("chr1", 300, "C", "T"),
        ],
        seed=5,
        permutations=99,
    )
    assert result["joint_signature_n_complete_reads"] == 0
    assert result["joint_signature_testable"] is False
    assert result["joint_signature_sensitivity_pass"] is False
    assert result["joint_signature_postselection_fdr_calibrated"] is False
    assert result["signature_p_analytic"] is None
    assert result["signature_q_global_by"] is None


def formal_pair(position: int, *, confirmed: bool = True, top: bool = True) -> dict[str, object]:
    return {
        "sample": "S",
        "focal_chrom": "chr1",
        "focal_pos": 1_000,
        "partner_pos": position,
        "endpoint_a_global_fdr_family_status": "ELIGIBLE_M2_EXACT_FAMILY",
        "endpoint_a_testable": True,
        "endpoint_a_p_fixed_margins_exact": 0.001,
        "endpoint_a_exact_bh_discovery": confirmed,
        "endpoint_a_exact_by_discovery": confirmed,
        "endpoint_a_formal_pair_by_confirmed": confirmed,
        "endpoint_a_permutable": confirmed,
        "top_coverage_marker": top,
    }


def formal_site(
    *,
    joint_pass: bool = True,
    complete_effect_positions: list[int] | None = None,
) -> dict[str, object]:
    if complete_effect_positions is None:
        complete_effect_positions = [1_100, 1_120]
    return {
        "sample": "S",
        "chrom": "chr1",
        "pos": 1_000,
        "n_partner_markers": 3,
        "m2_screen_eligible": True,
        "top_marker_positions": [1_100, 1_120],
        "joint_signature_complete_marker_effect_supported_positions": complete_effect_positions,
        "joint_signature_sensitivity_pass": joint_pass,
        "joint_signature_global_by_discovery": joint_pass,
    }


def test_multi_marker_candidate_requires_spaced_top_by_hits_and_joint_gate() -> None:
    pairs = [formal_pair(1_100), formal_pair(1_120), formal_pair(1_300, confirmed=False)]
    site = formal_site()
    audit = MODULE.finalize_site_inference([site], pairs)
    assert audit["pass"] is True
    assert site["n_pair_by_confirmed"] == 2
    assert site["pair_by_confirmed_positions"] == [1_100, 1_120]
    assert site["n_spatially_separated_pair_by_20bp"] == 2
    assert site["n_top_marker_pair_by_confirmed"] == 2
    assert site["multi_marker_molecular_haplotype_base_candidate"] is True
    assert "not_statistical_independence" in str(site["spaced_marker_20bp_contract"])

    blocked = formal_site(joint_pass=False)
    MODULE.finalize_site_inference([blocked], pairs)
    assert blocked["multi_marker_molecular_haplotype_base_candidate"] is False

    incomplete_effect_support = formal_site(complete_effect_positions=[1_100])
    MODULE.finalize_site_inference([incomplete_effect_support], pairs)
    assert incomplete_effect_support[
        "multi_marker_molecular_haplotype_base_candidate"
    ] is False


def replication_row(sample: str, read_names: list[str]) -> dict[str, object]:
    return {
        "sample": sample,
        "focal_chrom": "chr5",
        "focal_pos": 750_311,
        "focal_ref": "A",
        "focal_alt": "C",
        "partner_chrom": "chr5",
        "partner_pos": 750_429,
        "partner_ref": "C",
        "partner_alt": "A",
        "endpoint_a_formal_pair_by_confirmed": True,
        "endpoint_a_effect_direction": "PARTNER_ALT_ENRICHED_IN_DOMINANT_METHYL_GROUP",
        "endpoint_a_cramers_v": 0.8,
        "endpoint_a_delta_alt_fraction": 0.9,
        "endpoint_b_relation_compatibility": "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL",
        "cross_platform_exact_pair_present": False,
        "_permutation_payload": {"core_raw_read_names": read_names},
    }


def test_cross_platform_replication_reports_same_raw_molecule_overlap() -> None:
    hcc = replication_row("HCC1395", ["raw-a", "raw-shared"])
    dorado = replication_row("HCC1395_DORADO", ["raw-shared", "raw-b"])
    summary = MODULE.apply_cross_platform_replication([hcc, dorado])
    assert summary["biological_n"] == 1
    assert summary["n_exact_pairs_present_both"] == 1
    assert summary["n_exact_pairs_with_core_raw_read_name_overlap"] == 1
    assert hcc["cross_platform_core_read_name_intersection_n"] == 1
    assert hcc["cross_platform_core_read_name_union_n"] == 3
    assert hcc["cross_platform_core_read_name_jaccard"] == pytest.approx(1.0 / 3.0)
    assert hcc["cross_platform_core_read_name_overlap_present"] is True
    assert hcc["cross_platform_replication_status"] == (
        "SAME_MOLECULE_REPROCESSING_CONDITIONAL_BY_CONCORDANCE"
    )
    assert "_permutation_payload" not in hcc
    assert "_permutation_payload" not in dorado


def test_m2_screen_gate_is_fail_closed_and_truth_blind() -> None:
    clean = {
        "stable_null_multigroup": "true",
        "modal_assignment_ari_min": "0.9",
        "hp_axis_confound": "false",
        "technical_axis_confound": "false",
        "residual_unexplained_multigroup": "true",
        "phase_anchored_robust_epigenetic_candidate": "true",
        "cluster_sizes": '{"g1":100,"g2":100}',
        "m2_categorical_level_counts": {
            "hp_exact": 2,
            "hp_family": 2,
            "strand": 2,
        },
    }
    for prefix in ("hp_exact", "hp_family", "strand"):
        clean.update(
            {
                f"{prefix}_v": "0.1",
                f"{prefix}_p_perm": "0.5",
                f"{prefix}_n": "200",
                f"{prefix}_aligned": "false",
            }
        )
    for prefix in ("start", "end", "length", "mapq", "cpg_called"):
        clean.update(
            {
                f"{prefix}_eta2": "0.1",
                f"{prefix}_p_perm": "0.5",
                f"{prefix}_n": "200",
                f"{prefix}_aligned": "false",
            }
        )
    gate = MODULE.m2_screen_eligibility(clean)
    assert gate["eligible"] is True
    assert gate["evaluable"] is True
    assert gate["status"] == "ELIGIBLE_M2_RESIDUAL_UNEXPLAINED_AND_AXES_DETERMINATE"

    confounded = {
        **clean,
        "strand_v": "0.4",
        "strand_p_perm": "0.01",
        "strand_aligned": "true",
        "technical_axis_confound": "true",
        "residual_unexplained_multigroup": "false",
        "phase_anchored_robust_epigenetic_candidate": "false",
    }
    confounded_gate = MODULE.m2_screen_eligibility(confounded)
    assert confounded_gate["eligible"] is False
    assert confounded_gate["evaluable"] is True

    indeterminate = {
        **clean,
        "mapq_eta2": "0.2",
        "mapq_p_perm": "0.2",
    }
    indeterminate_gate = MODULE.m2_screen_eligibility(indeterminate)
    assert indeterminate_gate["eligible"] is False
    assert indeterminate_gate["evaluable"] is False
    assert indeterminate_gate["indeterminate_axes"] == ["mapq"]

    constant = {
        **clean,
        "hp_exact_v": "",
        "hp_exact_p_perm": "",
        "m2_categorical_level_counts": {
            "hp_exact": 1,
            "hp_family": 2,
            "strand": 2,
        },
    }
    constant_gate = MODULE.m2_screen_eligibility(constant)
    assert constant_gate["eligible"] is True
    assert constant_gate["constant_axes"] == ["hp_exact"]

    aligned_low_power = {
        **clean,
        "cluster_sizes": '{"g1":10,"g2":10}',
        "strand_v": "0.4",
        "strand_p_perm": "0.01",
        "strand_aligned": "true",
        "technical_axis_confound": "true",
        "residual_unexplained_multigroup": "false",
        "phase_anchored_robust_epigenetic_candidate": "false",
    }
    for prefix in ("hp_exact", "hp_family", "strand", "start", "end", "length", "mapq", "cpg_called"):
        aligned_low_power[f"{prefix}_n"] = "20"
    aligned_low_power_gate = MODULE.m2_screen_eligibility(aligned_low_power)
    assert aligned_low_power_gate["evaluable"] is False
    assert aligned_low_power_gate[
        "aligned_below_negative_evaluability_power_axes"
    ] == ["strand"]
    assert "strand" not in aligned_low_power_gate["low_power_axes"]

    low_power = {**clean, "cluster_sizes": '{"g1":10,"g2":10}'}
    for prefix in ("hp_exact", "hp_family", "strand", "start", "end", "length", "mapq", "cpg_called"):
        low_power[f"{prefix}_n"] = "20"
    low_power_gate = MODULE.m2_screen_eligibility(low_power)
    assert low_power_gate["eligible"] is False
    assert low_power_gate["evaluable"] is False
    assert set(low_power_gate["low_power_axes"]) == {
        "hp_exact", "hp_family", "strand", "start", "end", "length", "mapq", "cpg_called"
    }

    eleven_group_sizes = {f"g{index}": 10 for index in range(1, 12)}
    eleven_groups = {
        **clean,
        "cluster_sizes": json.dumps(eleven_group_sizes),
    }
    for prefix in (
        "hp_exact",
        "hp_family",
        "strand",
        "start",
        "end",
        "length",
        "mapq",
        "cpg_called",
    ):
        eleven_groups[f"{prefix}_n"] = "110"
    eleven_group_gate = MODULE.m2_screen_eligibility(eleven_groups)
    assert eleven_group_gate == {
        "contract": MODULE.M2_GATE.GATE_CONTRACT,
        "eligible": False,
        "evaluable": False,
        "status": (
            "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM:"
            "observed=11:maximum=10"
        ),
        "axis_statuses": {},
        "indeterminate_axes": [],
        "low_power_axes": [],
        "aligned_axes": [],
        "constant_axes": [],
        "aligned_below_negative_evaluability_power_axes": [],
        "categorical_level_counts": {
            "hp_exact": 2,
            "hp_family": 2,
            "strand": 2,
        },
    }
    eleven_groups_as_mapping = {
        **eleven_groups,
        "cluster_sizes": eleven_group_sizes,
    }
    assert MODULE.m2_screen_eligibility(eleven_groups_as_mapping) == eleven_group_gate

    zero_n = {**clean, "hp_exact_n": "0"}
    with pytest.raises(MODULE.M2_GATE.M2GateError, match="count/core-cluster drift"):
        MODULE.m2_screen_eligibility(zero_n)

    missing_gate = MODULE.m2_screen_eligibility({})
    assert missing_gate["evaluable"] is False
    assert missing_gate["status"].startswith("NOT_EVALUABLE_M2_SCREEN_FIELDS_MISSING")


def test_m2_positive_confound_is_evaluable_below_negative_power_threshold() -> None:
    row = {
        "stable_null_multigroup": "true",
        "modal_assignment_ari_min": "0.9",
        "hp_axis_confound": "true",
        "technical_axis_confound": "false",
        "residual_unexplained_multigroup": "false",
        "phase_anchored_robust_epigenetic_candidate": "false",
        "cluster_sizes": '{"g1":70,"g2":70}',
        "m2_categorical_level_counts": {
            "hp_exact": 2,
            "hp_family": 2,
            "strand": 2,
        },
    }
    for prefix in ("hp_exact", "hp_family", "strand"):
        row.update(
            {
                f"{prefix}_v": "0.1",
                f"{prefix}_p_perm": "0.5",
                f"{prefix}_n": "140",
                f"{prefix}_aligned": "false",
            }
        )
    for prefix in ("start", "end", "length", "mapq", "cpg_called"):
        row.update(
            {
                f"{prefix}_eta2": "0.1",
                f"{prefix}_p_perm": "0.5",
                f"{prefix}_n": "140",
                f"{prefix}_aligned": "false",
            }
        )
    row.update(
        {
            "hp_exact_v": "0.4",
            "hp_exact_p_perm": "0.01",
            "hp_exact_aligned": "true",
        }
    )

    gate = MODULE.m2_screen_eligibility(row)

    assert MODULE.M2_GATE.minimum_n_for_target_power(
        "hp_exact", "categorical", 2, 0.30
    ) == 152
    assert gate["eligible"] is False
    assert gate["evaluable"] is True
    assert gate["status"] == "INELIGIBLE_SCREEN_HP_CONFOUND"
    assert gate["indeterminate_axes"] == []
    assert gate["low_power_axes"] == []
    assert gate["aligned_below_negative_evaluability_power_axes"] == ["hp_exact"]
    assert gate["axis_statuses"]["hp_exact"]["positive_alignment_overrides_power"] is True
    assert all(
        gate["axis_statuses"][axis]["adequate_information_for_non_alignment"] is True
        for axis in ("hp_family", "strand", "start", "end", "length", "mapq", "cpg_called")
    )


def test_site_worker_entry_excludes_truth_ledger_and_topology_inputs() -> None:
    key = ("S", "chr1", 100)
    safe_fields = {
        "all_ssnv_vcf": {"path": "/safe/all.vcf.gz"},
        "all_ssnv_vcf_index": {"path": "/safe/all.vcf.gz.csi"},
        "latest_read_tag_sidecar": {"path": "/safe/tags.tsv.gz"},
        "latest_read_tag_sidecar_index": {"path": "/safe/tags.tsv.gz.tbi"},
        "raw_alignment": {"path": "/safe/raw.bam"},
        "raw_alignment_index": {"path": "/safe/raw.bam.bai"},
    }
    entry = {
        "sample": "S",
        "biological_id": "B",
        **safe_fields,
        "truth_fp_vcf": {"path": "/forbidden/fp.vcf.gz"},
        "truth_tp_vcf": {"path": "/forbidden/tp.vcf.gz"},
        "ledger": {"path": "/forbidden/ledger.tsv.gz"},
        "layered_region_view": {"path": "/forbidden/regions.jsonl"},
    }
    assignment = {
        "sample": "S",
        "chrom": "chr1",
        "pos": 100,
        "region_dir": "/safe/region",
        "posthoc": {"truth_label": "FP"},
    }
    site = {
        "region_dir": "/safe/region",
        "alt_readset_sha256": "digest",
        "ref": "A",
        "alt": "C",
        "screen_contract": MODULE.SCREEN_CONTRACT,
    }
    task = MODULE.build_tasks([entry], {key: assignment}, {key: site}, top_markers=3)[0]
    assert set(task["entry"]) == {"sample", "biological_id", *safe_fields}
    assert "posthoc" not in task["assignment"]


def test_primary_artifact_audit_pre_is_required_by_cli(tmp_path: Path) -> None:
    argv = [
        "--manifest",
        str(tmp_path / "manifest.json"),
        "--assignments",
        str(tmp_path / "assignments.jsonl.gz"),
        "--sites",
        str(tmp_path / "sites.tsv.gz"),
        "--intersubmod-root",
        str(tmp_path / "runs"),
        "--output-dir",
        str(tmp_path / "output"),
    ]
    with pytest.raises(SystemExit) as error:
        MODULE.parse_args(argv)
    assert error.value.code == 2

    args = MODULE.parse_args(
        [
            *argv,
            "--primary-artifact-audit-pre",
            str(tmp_path / "stable_primary_artifact_audit.v1_pre_downstream.json"),
            "--independent-m2-audit",
            str(tmp_path / "independent_m2_gate_recount.v3.json"),
            "--preflight-receipt",
            str(tmp_path / "cooccurrence_task_contract_preflight.json"),
        ]
    )
    assert args.primary_artifact_audit_pre == (
        tmp_path / "stable_primary_artifact_audit.v1_pre_downstream.json"
    )
    assert args.preflight_receipt == (
        tmp_path / "cooccurrence_task_contract_preflight.json"
    )
