from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "build_all_ssnv_report_artifact.py"
)
SPEC = importlib.util.spec_from_file_location("build_all_ssnv_report_artifact", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


SAMPLE_TRUTH_COUNTS = {
    "HCC1395": {"TP": 79_687, "FP": 0, "UNASSESSED": 0},
    "HCC1395_DORADO": {"TP": 79_739, "FP": 0, "UNASSESSED": 0},
    "COLO829": {"TP": 37_788, "FP": 0, "UNASSESSED": 0},
    "H1437": {"TP": 77_080, "FP": 0, "UNASSESSED": 0},
    "H2009": {"TP": 61_002, "FP": 0, "UNASSESSED": 93_463},
    "HCC1937": {"TP": 0, "FP": 7_745, "UNASSESSED": 10_945},
    "HCC1954": {"TP": 0, "FP": 0, "UNASSESSED": 22_400},
}

M1_PASS = {
    "HCC1395": 3,
    "HCC1395_DORADO": 2,
    "COLO829": 1,
    "H1437": 1,
    "H2009": 1,
    "HCC1937": 1,
    "HCC1954": 1,
}
M2_PASS = {
    "HCC1395": 2,
    "HCC1395_DORADO": 1,
    "COLO829": 1,
    "H1437": 0,
    "H2009": 1,
    "HCC1937": 0,
    "HCC1954": 1,
}
G1_PASS = {
    "HCC1395": 1,
    "HCC1395_DORADO": 1,
    "COLO829": 0,
    "H1437": 0,
    "H2009": 1,
    "HCC1937": 0,
    "HCC1954": 0,
}
STATUSES = ("PASS", "FAIL", "NOT_EVALUABLE", "NOT_RUN")


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_final_tsv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: MODULE.serialized_tsv_value(row.get(field))
                    for field in fields
                }
            )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def clear_tsv_rows(path: Path) -> None:
    header = path.read_text(encoding="utf-8").splitlines()[0]
    path.write_text(header + "\n", encoding="utf-8")


def artifact(path: Path) -> dict:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": MODULE.sha256(path),
    }


def refresh_final_receipt(inputs: MODULE.ReportInputs) -> None:
    receipt = read_json(inputs.final_receipt)
    receipt["inputs"] = read_json(inputs.final_dataset)["input_artifacts"]
    receipt["outputs"]["final_report_dataset"] = artifact(inputs.final_dataset)
    write_json(inputs.final_receipt, receipt)


def refresh_final_input_artifact(
    inputs: MODULE.ReportInputs, role: str, path: Path
) -> None:
    final = read_json(inputs.final_dataset)
    final["input_artifacts"][role] = artifact(path)
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)


def metric_from_statuses(counts: dict[str, int], claim_id: str) -> dict:
    numerator = counts["PASS"]
    denominator = numerator + counts["FAIL"]
    return {
        "numerator": numerator,
        "denominator": denominator,
        "ratio": numerator / denominator if denominator else None,
        "not_evaluable": counts["NOT_EVALUABLE"],
        "not_run": counts["NOT_RUN"],
        "population": sum(counts.values()),
        "denominator_definition": (
            MODULE.M1_DENOMINATOR_DEFINITION
            if claim_id == "M1"
            else MODULE.M2_DENOMINATOR_DEFINITION
            if claim_id == "M2"
            else "claim-specific evaluable sites (PASS + FAIL)"
        ),
    }


def add_status_counts(rows: list[dict[str, int]]) -> dict[str, int]:
    return {status: sum(row[status] for row in rows) for status in STATUSES}


def sample_status_counts(
    sample: str,
    *,
    with_candidate: bool,
    eleven_group_site: bool = False,
) -> dict[str, dict[str, int]]:
    population = sum(SAMPLE_TRUTH_COUNTS[sample].values())
    m1_pass = M1_PASS[sample]
    m1_fail = population - m1_pass
    group_limit_not_evaluable = int(eleven_group_site and sample == "HCC1954")
    original_m2_pass = M2_PASS[sample]
    m2_pass = original_m2_pass - group_limit_not_evaluable
    g1_pass = G1_PASS[sample]
    g2_pass = int(with_candidate and sample == "HCC1395")
    candidate_row = int(with_candidate and sample == "HCC1395")
    return {
        "M1": {
            "PASS": m1_pass,
            "FAIL": m1_fail,
            "NOT_EVALUABLE": 0,
            "NOT_RUN": 0,
        },
        "M2": {
            "PASS": m2_pass,
            "FAIL": m1_pass - original_m2_pass,
            "NOT_EVALUABLE": group_limit_not_evaluable,
            "NOT_RUN": population - m1_pass,
        },
        "G1": {
            "PASS": g1_pass,
            "FAIL": m2_pass - g1_pass,
            "NOT_EVALUABLE": 0,
            "NOT_RUN": population - m2_pass,
        },
        "G2": {
            "PASS": g2_pass,
            "FAIL": m2_pass - g2_pass,
            "NOT_EVALUABLE": 0,
            "NOT_RUN": population - m2_pass,
        },
        "R1": {
            "PASS": candidate_row,
            "FAIL": 0,
            "NOT_EVALUABLE": 0,
            "NOT_RUN": population - candidate_row,
        },
        "B1": {
            "PASS": candidate_row,
            "FAIL": 0,
            "NOT_EVALUABLE": 0,
            "NOT_RUN": population - candidate_row,
        },
        "C1": {
            "PASS": 0,
            "FAIL": 0,
            "NOT_EVALUABLE": candidate_row,
            "NOT_RUN": population - candidate_row,
        },
        "L1": {
            "PASS": 0,
            "FAIL": 0,
            "NOT_EVALUABLE": candidate_row,
            "NOT_RUN": population - candidate_row,
        },
        "L2": {
            "PASS": 0,
            "FAIL": 0,
            "NOT_EVALUABLE": 0,
            "NOT_RUN": population,
        },
    }


def split_status_counts(
    counts: dict[str, int], truth_populations: dict[str, int]
) -> dict[str, dict[str, int]]:
    remaining = dict(counts)
    result: dict[str, dict[str, int]] = {}
    for truth in MODULE.TRUTH_LABELS:
        capacity = truth_populations[truth]
        cell = {status: 0 for status in STATUSES}
        for status in STATUSES:
            take = min(capacity, remaining[status])
            cell[status] = take
            remaining[status] -= take
            capacity -= take
        assert capacity == 0
        result[truth] = cell
    assert all(value == 0 for value in remaining.values())
    return result


def build_funnel(
    with_candidate: bool,
    *,
    eleven_group_site: bool = False,
) -> tuple[dict, dict[str, dict[str, int]]]:
    sample_counts = {
        sample: sample_status_counts(
            sample,
            with_candidate=with_candidate,
            eleven_group_site=eleven_group_site,
        )
        for sample in MODULE.EXPECTED_DATASETS
    }
    per_sample = {
        sample: {
            claim_id: metric_from_statuses(sample_counts[sample][claim_id], claim_id)
            for claim_id in MODULE.CLAIM_IDS
        }
        for sample in MODULE.EXPECTED_DATASETS
    }
    split_cells: dict[str, dict[str, dict[str, int]]] = {}
    for sample in MODULE.EXPECTED_DATASETS:
        for claim_id in MODULE.CLAIM_IDS:
            split = split_status_counts(
                sample_counts[sample][claim_id], SAMPLE_TRUTH_COUNTS[sample]
            )
            for truth in MODULE.TRUTH_LABELS:
                split_cells.setdefault(f"{sample}|{truth}", {})[claim_id] = split[truth]
    sample_by_truth = {
        key: {
            claim_id: metric_from_statuses(split_cells[key][claim_id], claim_id)
            for claim_id in MODULE.CLAIM_IDS
        }
        for key in (
            f"{sample}|{truth}"
            for sample in MODULE.EXPECTED_DATASETS
            for truth in MODULE.TRUTH_LABELS
        )
    }
    truth_strata = {
        truth: {
            claim_id: metric_from_statuses(
                add_status_counts(
                    [split_cells[f"{sample}|{truth}"][claim_id] for sample in MODULE.EXPECTED_DATASETS]
                ),
                claim_id,
            )
            for claim_id in MODULE.CLAIM_IDS
        }
        for truth in MODULE.TRUTH_LABELS
    }
    pooled_statuses = {
        claim_id: add_status_counts(
            [sample_counts[sample][claim_id] for sample in MODULE.EXPECTED_DATASETS]
        )
        for claim_id in MODULE.CLAIM_IDS
    }
    pooled = {
        claim_id: metric_from_statuses(pooled_statuses[claim_id], claim_id)
        for claim_id in MODULE.CLAIM_IDS
    }
    return {
        "pooled": pooled,
        "per_sample": per_sample,
        "truth_strata": truth_strata,
        "sample_by_truth": sample_by_truth,
    }, pooled_statuses


def claim_ladder(funnel: dict) -> list[dict]:
    rows = []
    dorado = funnel["per_sample"]["HCC1395_DORADO"]
    for claim_id in MODULE.CLAIM_IDS:
        metric = funnel["pooled"][claim_id]
        biological_numerator = metric["numerator"] - dorado[claim_id]["numerator"]
        biological_denominator = metric["denominator"] - dorado[claim_id]["denominator"]
        rows.append(
            {
                "claim_id": claim_id,
                "claim_name": MODULE.CLAIM_NAMES[claim_id],
                "dataset_numerator": metric["numerator"],
                "dataset_denominator": metric["denominator"],
                "dataset_ratio": metric["ratio"],
                "dataset_not_evaluable": metric["not_evaluable"],
                "dataset_not_run": metric["not_run"],
                "biological_numerator": biological_numerator,
                "biological_denominator": biological_denominator,
                "biological_ratio": (
                    biological_numerator / biological_denominator
                    if biological_denominator
                    else None
                ),
                "status": MODULE.aggregate_status(metric),
                "denominator_definition": metric["denominator_definition"],
                "guardrail": "No automatic claim promotion.",
                "automatic_upgrade_prohibited": True,
            }
        )
    return rows


def pair_row(
    partner_pos: int,
    *,
    exact_p: float,
    effect: float,
    delta: float,
    formal: bool,
) -> dict:
    return {
        "sample": "HCC1395",
        "focal_chrom": "chr1",
        "focal_pos": 100,
        "focal_ref": "A",
        "focal_alt": "C",
        "partner_chrom": "chr1",
        "partner_pos": partner_pos,
        "partner_ref": "G",
        "partner_alt": "T",
        "endpoint_a_groups": json.dumps(["methyl-A", "methyl-B"]),
        "endpoint_a_table": json.dumps([[8, 2], [1, 9]]),
        "endpoint_a_testable": "true",
        "endpoint_a_p_fixed_margins_exact": exact_p,
        "endpoint_a_q_global_by": 0.04 if formal else 0.20,
        "endpoint_a_cramers_v": effect,
        "endpoint_a_delta_alt_fraction": delta,
        "endpoint_a_p_conditional_perm": 0.01 if formal else 0.20,
        "endpoint_a_permutations": 999,
        "endpoint_a_permutable": "true",
        "endpoint_a_conditional_status": "PERMUTABLE",
        "endpoint_a_conditional_sensitivity_pass": str(formal).lower(),
        "endpoint_a_formal_pair_by_confirmed": str(formal).lower(),
        "callability_testable": "false",
        "callability_q_global_by": "",
        "callability_cramers_v": "",
        "callability_noncallable_core_reads": 0,
        "callability_gate_status": "PASS_ALL_CORE_READS_CALLABLE",
        "callability_gate_pass": "true",
        "endpoint_b_state_counts": json.dumps(
            {"RR": 8, "AR": 2, "RA": 1, "AA": 9, "O": 0, "X": 0}
        ),
        "endpoint_b_n_called_depth": 20,
        "endpoint_b_error_ceiling": 0.02,
        "endpoint_b_error_model_confidence": 0.9833333333333333,
        "endpoint_b_familywise_confidence": 0.95,
        "endpoint_b_relation_family_size": 3,
        "endpoint_b_multiplicity_method": "bonferroni_three_relation_models",
        "endpoint_b_minimum_zero_violation_depth": 203,
        "endpoint_b_focal_ancestor_violation_p_exact": 0.016,
        "endpoint_b_focal_ancestor_violation_upper_bound": 0.394,
        "endpoint_b_focal_ancestor_violation_threshold": 0.02,
        "endpoint_b_focal_ancestor_violation_status": "VIOLATES_FIXED_ERROR_CEILING",
        "endpoint_b_partner_ancestor_violation_p_exact": 0.016,
        "endpoint_b_partner_ancestor_violation_upper_bound": 0.394,
        "endpoint_b_partner_ancestor_violation_threshold": 0.02,
        "endpoint_b_partner_ancestor_violation_status": "VIOLATES_FIXED_ERROR_CEILING",
        "endpoint_b_branching_violation_p_exact": 0.001,
        "endpoint_b_branching_violation_upper_bound": 0.895,
        "endpoint_b_branching_violation_threshold": 0.02,
        "endpoint_b_branching_violation_status": "VIOLATES_FIXED_ERROR_CEILING",
        "endpoint_b_complete_four_state_testable": "true",
        "endpoint_b_relation_compatibility": (
            "NOT_IDENTIFIABLE_FIXED_ERROR_CEILING"
        ),
        "endpoint_b_compatible_relation_models": json.dumps([]),
        "endpoint_b_n_compatible_relation_models": 0,
        "focal_truth_label": "TP",
        "partner_truth_label": "TP",
        "focal_ssnv_branch": "retained",
        "partner_ssnv_branch": "retained",
        "focal_component_id": "component-focal",
        "partner_component_id": "component-partner",
        "topology_scope": "SHARED_RETAINED_REGION",
        "topology_region": "chr1:1-1000",
        "topology_order_status": "PAIRWISE_ORDER_UNRESOLVED",
        "topology_claim_guardrail": (
            "Topology is posthoc context; not proof of cellular ancestry."
        ),
    }


def build_fixture(
    tmp_path: Path,
    *,
    with_candidate: bool,
    with_history: bool = False,
    eleven_group_site: bool = False,
) -> tuple[MODULE.ReportInputs, Path]:
    root = tmp_path / "machine-local-inputs"
    repo_root = tmp_path / "portable-repo-root"
    root.mkdir(parents=True)
    repo_root.mkdir(parents=True)

    claim_contract_path = root / "claim-contract-v3.md"
    claim_contract_path.write_text(
        (SCRIPT.parents[1] / "claim-contract-v3.md").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    final_builder_path = root / "build_all_ssnv_final_report_dataset.py"
    final_builder_path.write_text(
        "#!/usr/bin/env python3\n# Frozen test producer identity.\n",
        encoding="utf-8",
    )

    layered_samples = []
    for sample in MODULE.EXPECTED_DATASETS:
        producer_receipt_path = root / "producer-receipts" / sample / "producer_capture_receipt_v2.json"
        write_json(
            producer_receipt_path,
            {
                "schema_name": "intersubmod.longphase_raw_all_capture_receipt",
                "schema_version": "2.0.0",
                "sample": sample,
                "bam_output_policy": {
                    "transport": "named_fifo",
                    "persisted_bam": False,
                    "regular_bam_count": 0,
                    "consumed_fifo_path": str(
                        (producer_receipt_path.parent / "consumed_tagged_bam.fifo").resolve()
                    ),
                    "is_fifo_at_closeout": True,
                },
            },
        )
        producer_identity = artifact(producer_receipt_path)
        layered_samples.append(
            {
                "sample": sample,
                "read_tags": {
                    "producer_capture_receipt": {
                        "path": producer_identity["path"],
                        "identity": {
                            "policy": "full_sha256",
                            "size_bytes": producer_identity["size_bytes"],
                            "sha256": producer_identity["sha256"],
                        },
                    }
                },
            }
        )

    layered_path = root / "layered_manifest.json"
    write_json(
        layered_path,
        {
            "schema_name": MODULE.LAYERED_SCHEMA,
            "schema_version": MODULE.LAYERED_SCHEMA_VERSION,
            "dataset_count": 7,
            "biological_sample_count": 6,
            "production_summary": {"all_pass": True, "passed_dataset_count": 7},
            "analysis_contract": {
                "longphase_input_contract": "normalized_ClairS_raw_all",
                "tree_input_contract": "longphase_s_recalibrated_FILTER_PASS",
                "read_tag_mode": "external_sidecar",
                "embedded_tag_policy": "ignore",
                "require_exact_join": True,
                "tagging_semantics": "latest LongPhase-S HP/PS projection",
            },
            "samples": layered_samples,
        },
    )

    manifest_path = root / "all_ssnv_input_manifest.json"
    samples = []
    truth_totals = {truth: 0 for truth in MODULE.TRUTH_LABELS}
    for sample in MODULE.EXPECTED_DATASETS:
        truth_counts = SAMPLE_TRUTH_COUNTS[sample]
        for truth in MODULE.TRUTH_LABELS:
            truth_totals[truth] += truth_counts[truth]
        samples.append(
            {
                "sample": sample,
                "counts": {
                    "all_ssnv": sum(truth_counts.values()),
                    "truth_tp": truth_counts["TP"],
                    "truth_fp": truth_counts["FP"],
                    "truth_unassessed": truth_counts["UNASSESSED"],
                },
            }
        )
    write_json(
        manifest_path,
        {
            "schema_name": MODULE.MANIFEST_SCHEMA,
            "schema_version": MODULE.MANIFEST_SCHEMA_VERSION,
            "pass": True,
            "samples": samples,
            "totals": {
                "all_ssnv": MODULE.EXPECTED_SITES,
                "truth_tp": truth_totals["TP"],
                "truth_fp": truth_totals["FP"],
                "truth_unassessed": truth_totals["UNASSESSED"],
            },
            "layered_root": artifact(layered_path),
        },
    )

    screen_path = root / "all_ssnv_summary.json"
    terminal_audit = {
        "authoritative_tag_source": "same_run_LongPhase_S_external_HP_PS_sidecar",
        "embedded_reads_tsv_hp_used_for_analysis": False,
        "join_occurs_before_focal_ALT_selection": True,
        "n_sites": MODULE.EXPECTED_SITES,
        "join_status_counts": {"PASS": MODULE.EXPECTED_SITES},
        "all_sites_pass": True,
        "n_reads_tsv_site_rows": 1234,
        "n_exact_hp_ps_site_read_joins": 1234,
        "every_reads_tsv_row_joined": True,
        "all_projection_identities_unique": True,
        "n_projection_multimatch_site_reads": 0,
        "pass": True,
    }
    write_json(
        screen_path,
        {
            "schema_name": MODULE.SCREEN_SCHEMA,
            "schema_version": MODULE.SCREEN_SCHEMA_VERSION,
            "pass": True,
            "scope": {
                "full_469849": True,
                "selected_datasets": list(MODULE.EXPECTED_DATASETS),
                "expected_sites": MODULE.EXPECTED_SITES,
                "processed_sites": MODULE.EXPECTED_SITES,
            },
            "pooled_site_weighted": {"n_sites": MODULE.EXPECTED_SITES},
            "latest_hp_ps_terminal_join_audit": terminal_audit,
        },
    )
    screen_sites_path = root / "all_ssnv_site_results.tsv.gz"
    screen_assignments_path = root / "all_ssnv_stable_assignments.jsonl.gz"
    screen_sites_path.write_text("fixture-screen-sites\n", encoding="utf-8")
    screen_assignments_path.write_text("fixture-screen-assignments\n", encoding="utf-8")
    screen_receipt_path = root / "screen_run_manifest.json"
    write_json(
        screen_receipt_path,
        {
            "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
            "schema_version": "1.2.0",
            "status": "EXECUTION_PASS",
            "pass": True,
            "started_at_utc": "2026-07-15T00:00:00+00:00",
            "finished_at_utc": "2026-07-15T00:02:00+00:00",
            "execution": {"recovery_merge": False},
        },
    )

    site_path = root / "methyl_ssnv_site_results.tsv"
    write_tsv(
        site_path,
        [
            {
                "sample": "HCC1395",
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "C",
                "truth_label": "TP",
                "ssnv_branch": "retained",
                "component_id": "component-focal",
                "modal_assignment_ari_min": 0.92,
                "hp_axis_confound": "false",
                "technical_axis_confound": "false",
                "residual_unexplained_multigroup": "true",
                "joint_signature_n_complete_reads": 20,
                "joint_signature_testable": "true",
                "joint_signature_groups": json.dumps(["methyl-A", "methyl-B"]),
                "joint_signature_categories": json.dumps(["R|R", "R|A"]),
                "joint_signature_table": json.dumps([[8, 2], [1, 9]]),
                "joint_signature_conditional_status": (
                    "PERMUTABLE"
                ),
                "joint_signature_p_conditional_perm": (
                    0.01 if with_candidate else 0.20
                ),
                "joint_signature_permutations": 999,
                "joint_signature_permutable": "true",
                "joint_signature_sensitivity_pass": str(with_candidate).lower(),
                "joint_signature_global_fdr_family_status": (
                    "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
                ),
                "joint_signature_q_global_by": 0.04 if with_candidate else 0.20,
                "joint_signature_global_by_discovery": str(with_candidate).lower(),
                "multi_marker_molecular_haplotype_base_candidate": str(
                    with_candidate
                ).lower(),
            }
        ],
    )
    pair_path = root / "methyl_ssnv_pair_results.tsv"
    write_tsv(
        pair_path,
        [
            pair_row(
                120,
                exact_p=0.02,
                effect=0.70,
                delta=0.60,
                formal=with_candidate,
            ),
            pair_row(
                150,
                exact_p=0.02,
                effect=0.80,
                delta=0.50,
                formal=False,
            ),
        ],
    )
    cooccurrence_summary_path = root / "cooccurrence_summary.json"
    write_json(
        cooccurrence_summary_path,
        {
            "schema_name": MODULE.COOCCURRENCE_SCHEMA,
            "schema_version": MODULE.COOCCURRENCE_SCHEMA_VERSION,
            "pass": True,
            "samples": list(MODULE.EXPECTED_DATASETS),
            "site_pair_count_reconciliation": {"pass": True},
        },
    )
    cooccurrence_receipt_path = root / "cooccurrence_run_receipt.json"
    write_json(
        cooccurrence_receipt_path,
        {
            "schema_name": "intersubmod.methyl_ssnv_cooccurrence_validation.run_receipt",
            "schema_version": "2.0.0",
            "started_at_utc": "2026-07-15T00:02:10+00:00",
            "finished_at_utc": "2026-07-15T00:02:20+00:00",
            "pass": True,
        },
    )
    tumor_ref_sites_path = root / "tumor_ref_sites.tsv.gz"
    tumor_ref_sites_path.write_text("fixture-tumor-ref-sites\n", encoding="utf-8")
    tumor_ref_summary_path = root / "tumor_ref_summary.json"
    tumor_ref_receipt_path = root / "tumor_ref_run_manifest.json"
    write_json(tumor_ref_summary_path, {"schema_name": "fixture.tumor_ref", "pass": True})
    write_json(
        tumor_ref_receipt_path,
        {
            "schema_name": "intersubmod.all_ssnv_tumor_ref_controls.run_manifest",
            "schema_version": "1.0.0",
            "started_at_utc": "2026-07-15T00:02:10+00:00",
            "finished_at_utc": "2026-07-15T00:02:20+00:00",
            "pass": True,
        },
    )
    primary_pre_path = root / "stable_primary.pre.json"
    primary_post_path = root / "stable_primary.post.json"
    write_json(
        primary_pre_path,
        {
            "schema_name": "intersubmod.stable_primary_artifact_audit",
            "started_at_utc": "2026-07-15T00:02:00+00:00",
            "finished_at_utc": "2026-07-15T00:02:05+00:00",
            "pass": True,
        },
    )
    write_json(
        primary_post_path,
        {
            "schema_name": "intersubmod.stable_primary_artifact_audit",
            "started_at_utc": "2026-07-15T00:02:30+00:00",
            "finished_at_utc": "2026-07-15T00:02:35+00:00",
            "pass": True,
        },
    )

    reconciliation_path = root / "output_reconciliation.json"
    exact_totals = {
        field: MODULE.EXPECTED_SITES
        for field in ("expected_sites", "reads_keys", "methylation_keys", "bernoulli_keys")
    }
    zero_totals = {
        field: 0
        for field in (
            "region_failures",
            "reads_missing",
            "reads_extra",
            "reads_duplicates",
            "reads_empty",
            "methylation_missing",
            "methylation_extra",
            "methylation_duplicates",
            "methylation_empty",
            "bernoulli_missing",
            "bernoulli_extra",
            "bernoulli_duplicates",
            "bernoulli_empty",
        )
    }
    write_json(
        reconciliation_path,
        {
            "schema_name": MODULE.RECONCILIATION_SCHEMA,
            "schema_version": MODULE.RECONCILIATION_SCHEMA_VERSION,
            "pass": True,
            "status": "EXECUTION_PASS",
            "totals": {**exact_totals, **zero_totals},
            "datasets": [
                {"sample": sample, "pass": True}
                for sample in MODULE.EXPECTED_DATASETS
            ],
        },
    )

    immutability_path = root / "frozen_input_immutability.post.json"
    write_json(
        immutability_path,
        {
            "schema_name": MODULE.IMMUTABILITY_SCHEMA,
            "schema_version": MODULE.IMMUTABILITY_SCHEMA_VERSION,
            "pass": True,
            "audit_phase": "post_run",
            "created_at_utc": "2026-07-15T00:03:00+00:00",
            "totals": {
                "n_samples": 7,
                "n_sample_pass": 7,
                "n_artifacts": 77,
                "n_artifact_pass": 77,
            },
        },
    )

    tree_path = root / "latest_tree_input_contract_audit.json"
    write_json(
        tree_path,
        {
            "schema_name": MODULE.TREE_AUDIT_SCHEMA,
            "schema_version": MODULE.TREE_AUDIT_SCHEMA_VERSION,
            "pass": True,
            "totals": {"all_ssnv": MODULE.EXPECTED_SITES},
            "top_level_checks": {
                "same_run_longphase_s": True,
                "recalibrated_filter_pass": True,
            },
            "samples": [
                {"sample": sample, "pass": True}
                for sample in MODULE.EXPECTED_DATASETS
            ],
        },
    )

    reference_audit_path = root / "extraction_reference_identity_audit.v1.json"
    write_json(
        reference_audit_path,
        {
            "schema_name": MODULE.REFERENCE_AUDIT_SCHEMA,
            "schema_version": MODULE.REFERENCE_AUDIT_SCHEMA_VERSION,
            "pass": True,
            "task_type": "B_comprehensive_validation",
            "pass_semantics": MODULE.REFERENCE_AUDIT_PASS_SEMANTICS,
            "scope": "7 datasets; 469,849 extraction dataset-sites; test reference",
            "reference": {
                "path": "/fixture/reference.fa",
                "size_bytes": 1_000,
                "full_sha256": "1" * 64,
                "fai": {
                    "path": "/fixture/reference.fa.fai",
                    "size_bytes": 100,
                    "full_sha256": "2" * 64,
                },
            },
            "sample_receipts": [
                {
                    "sample": sample,
                    "path": f"/fixture/{sample}/run_receipt.json",
                    "sha256": "3" * 64,
                    "receipt_pass": True,
                    "exit_code": 0,
                    "reference_path_equal": True,
                    "command_reference_path_equal": True,
                }
                for sample in MODULE.EXPECTED_DATASETS
            ],
            "checks": {
                "seven_extraction_receipts_present": True,
                "seven_extraction_receipts_pass": True,
                "current_reference_full_sha256_captured": True,
            },
            "limitations": [
                "Fixture preserves the independent post-extraction binding limitation."
            ],
        },
    )

    funnel, pooled_statuses = build_funnel(
        with_candidate,
        eleven_group_site=eleven_group_site,
    )
    candidate_key = json.dumps(
        ["HCC1395", "chr1", 100, "A", "C"], separators=(",", ":")
    )
    witness_pair_key = json.dumps(
        ["HCC1395", "chr1", 100, "A", "C", "chr1", 200, "G", "T"],
        separators=(",", ":"),
    )
    candidate_catalog = (
        [
            {
                "candidate_key": candidate_key,
                "sample": "HCC1395",
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "C",
                "g2_status": "PASS",
                "m2_screen_gate_contract": (
                    "m2-measured-axis-v4_asymmetric-positive-confound-and-observed-categorical-levels"
                ),
                "m2_screen_evaluable": True,
                "m2_screen_eligibility_status": (
                    "ELIGIBLE_M2_RESIDUAL_UNEXPLAINED_AND_AXES_DETERMINATE"
                ),
                "m2_indeterminate_axes": [],
                "m2_low_power_axes": [],
                "joint_signature_global_fdr_family_status": (
                    "ELIGIBLE_M2_JOINT_SIGNATURE_GLOBAL_FDR_FAMILY"
                ),
                "joint_signature_q_global_bh": 0.02,
                "joint_signature_q_global_by": 0.04,
                "joint_signature_global_bh_discovery": True,
                "joint_signature_global_by_discovery": True,
                "n_same_pair_four_state_witnesses": 1,
                "b1_prespecified_pair_key": witness_pair_key,
                "b1_uses_posthoc_compatible_pair_search": False,
                "b1_prespecified_pair_is_witness": True,
                "b1_status": "PASS",
                "tumor_ref_status": "PASS",
                "normal_control_status": "PASS",
                "normal_called_reads": 12,
                "normal_alt_reads": 0,
            }
        ]
        if with_candidate
        else []
    )
    candidate_witness_pairs = (
        [
            {
                "witness_pair_key": witness_pair_key,
                "sample": "HCC1395",
                "focal_chrom": "chr1",
                "focal_pos": 100,
                "focal_ref": "A",
                "focal_alt": "C",
                "partner_chrom": "chr1",
                "partner_pos": 200,
                "partner_ref": "G",
                "partner_alt": "T",
                "global_by_q": 0.01,
                "four_state_familywise_confidence": 0.95,
                "four_state_relation_family_size": 3,
                "four_state_multiplicity_method": "bonferroni_three_relation_models",
                "four_state_minimum_zero_violation_depth": 203,
                "same_pair_four_state_witness": True,
                "b1_prespecified_pair": True,
            }
        ]
        if with_candidate
        else []
    )
    final_path = root / "final_report_dataset.json"
    final = {
        "schema_name": MODULE.FINAL_SCHEMA,
        "schema_version": MODULE.FINAL_SCHEMA_VERSION,
        "pass": True,
        "task_type": "B_comprehensive_validation",
        "pass_semantics": MODULE.PASS_SEMANTICS,
        "scope": {
            "datasets": list(MODULE.EXPECTED_DATASETS),
            "dataset_count": 7,
            "biological_sample_count": 6,
            "expected_screen_sites": MODULE.EXPECTED_SITES,
            "observed_screen_sites": MODULE.EXPECTED_SITES,
            "technical_replicate": {
                "primary": "HCC1395",
                "replicate": "HCC1395_DORADO",
                "counts_as_independent_biological_n": False,
            },
        },
        "input_artifacts": {
            "manifest": artifact(manifest_path),
            "screen_sites": artifact(screen_sites_path),
            "screen_assignments": artifact(screen_assignments_path),
            "screen_summary": artifact(screen_path),
            "screen_receipt": artifact(screen_receipt_path),
            "cooccurrence_sites": artifact(site_path),
            "cooccurrence_pairs": artifact(pair_path),
            "cooccurrence_summary": artifact(cooccurrence_summary_path),
            "cooccurrence_receipt": artifact(cooccurrence_receipt_path),
            "tumor_ref_sites": artifact(tumor_ref_sites_path),
            "tumor_ref_summary": artifact(tumor_ref_summary_path),
            "tumor_ref_receipt": artifact(tumor_ref_receipt_path),
            "primary_artifact_audit_pre": artifact(primary_pre_path),
            "primary_artifact_audit_post": artifact(primary_post_path),
        },
        "funnel_metrics": funnel,
        "m1_operational_screen": {
            "status_semantics": "FLAGGED_VS_NOT_FLAGGED_OPERATIONAL_SCREEN_ONLY",
            "denominator_definition": MODULE.M1_DENOMINATOR_DEFINITION,
            "n_all_dataset_sites": MODULE.EXPECTED_SITES,
            "n_screen_evaluable": MODULE.EXPECTED_SITES - 2,
            "n_screen_not_evaluable": 2,
            "n_flagged_stable_multigroup": sum(M1_PASS.values()),
            "n_not_flagged_all": MODULE.EXPECTED_SITES - sum(M1_PASS.values()),
            "flag_yield": sum(M1_PASS.values()) / MODULE.EXPECTED_SITES,
            "flag_yield_among_screen_evaluable": (
                sum(M1_PASS.values()) / (MODULE.EXPECTED_SITES - 2)
            ),
            "n_evaluable_not_flagged": (
                MODULE.EXPECTED_SITES - 2 - sum(M1_PASS.values())
            ),
            "global_null_validity_exported_for_nonstable_sites": False,
            "nonflagged_scientific_interpretation": (
                "NOT_IDENTIFIABLE_AS_TRUE_NEGATIVE_VS_NULL_INVALID_FROM_SITE_TSV"
            ),
            "biological_prevalence_estimate": None,
        },
        "m2_evaluability_contract": {
            "gate_contract": MODULE.M2_GATE_CONTRACT,
            "denominator_definition": MODULE.M2_DENOMINATOR_DEFINITION,
            "minimum_supported_methyl_groups": 2,
            "maximum_supported_methyl_groups": 10,
            "categorical_planning_level_ceilings": dict(
                MODULE.M2_CATEGORICAL_PLANNING_LEVEL_CEILINGS
            ),
            "assignment_observed_levels_role": (
                MODULE.M2_ASSIGNMENT_OBSERVED_LEVELS_ROLE
            ),
            "n_m1_pass": sum(M1_PASS.values()),
            "n_m2_evaluable": sum(M1_PASS.values()) - int(eleven_group_site),
            "n_m2_not_evaluable": int(eleven_group_site),
            "not_evaluable_reason_counts": (
                {
                    "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM": 1
                }
                if eleven_group_site
                else {}
            ),
            "n_group_count_exceeds_planning_model_maximum": int(eleven_group_site),
            "group_count_exceeds_planning_model_examples": (
                [
                    {
                        "dataset": "HCC1954",
                        "chrom": "chr5",
                        "pos": 751076,
                        "ref": "G",
                        "alt": "A",
                        "observed_methyl_groups": 11,
                        "maximum_supported_methyl_groups": 10,
                        "m1_status": "PASS",
                        "m2_status": "NOT_EVALUABLE",
                        "g1_status": "NOT_RUN",
                        "g2_status": "NOT_RUN",
                        "b1_status": "NOT_RUN",
                        "reason": (
                            "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_"
                            "MAXIMUM:observed=11:maximum=10"
                        ),
                    }
                ]
                if eleven_group_site
                else []
            ),
            "group_count_exceeds_examples_complete": True,
            "group_count_exceeds_claim_behavior": {
                "M1": "PASS retained",
                "M2": "NOT_EVALUABLE excluded from PASS/FAIL denominator",
                "G1": "NOT_RUN",
                "G2": "NOT_RUN",
                "B1": "NOT_RUN",
            },
        },
        "tumor_ref_source_identity_attestation": {
            "status": "NOT_INCLUDED_INTERMEDIATE_TERMINAL_BUILD",
            "release_gate_pass": False,
            "publishable_task_b_release": False,
            "interpretation": (
                "Computational intermediate only; a passing bounded retrospective source-file "
                "identity receipt is required for the final release bundle."
            ),
        },
        "background_control_replication_gate": {
            "contract": MODULE.BACKGROUND_CONTROL_REPLICATION_GATE_CONTRACT,
            "applies_to": ["tumor_REF", "matched_normal_REF"],
            "required_conditions": [
                "coarse_ng>=2",
                "modal_fraction>=0.7_via_unstable_false",
            ],
            "membership_ari_minimum_required": False,
            "relation_to_primary_m1_replication_flags": (
                MODULE.BACKGROUND_CONTROL_RELATION_TO_PRIMARY_M1
            ),
            "b1_pass_direction": "requires_no_lenient_background_replication",
            "false_positive_direction": (
                "cannot_increase_B1_passes_vs_ARI_qualified_predicate_on_same_background_payload"
            ),
            "false_negative_direction": (
                "may_conservatively_reduce_B1_passes_when_K_is_stable_but_membership_is_not"
            ),
            "scientific_interpretation": (
                "background nonreplication guardrail, not an exact primary-M1 replay"
            ),
        },
        "claim_ladder": claim_ladder(funnel),
        "candidate_catalog": candidate_catalog,
        "candidate_witness_pairs": candidate_witness_pairs,
        "m2_axis_statistic_provenance": {
            "axis_effect_and_permutation_p_source": (
                "source_locked_focal_alt_multigroup_screen_producer"
            ),
            "screen_recovery_source_validation": {
                "mode": "source_locked_prefix_plus_seed_parallel_replacement",
                "pinned_analyzer_sha256": "a" * 64,
                "prefix_source_lock_pass": True,
                "serial_parallel_exact_equivalence_pass": True,
                "real_nested_equivalence_fixture": {},
            },
            "downstream_raw_read_axis_statistic_recomputed": False,
            "downstream_axis_classification_recomputed": True,
            "downstream_checks": [
                "axis sample-size reconciliation",
                "499-permutation add-one p-value grid",
                "effect threshold classification",
                "80-percent planning-power evaluability",
                "assignment-derived categorical constant-axis proof",
                "asymmetric positive-confound versus negative-evaluability power decision",
            ],
            "categorical_planning_level_ceilings": dict(
                MODULE.M2_CATEGORICAL_PLANNING_LEVEL_CEILINGS
            ),
            "assignment_observed_levels_role": (
                MODULE.M2_ASSIGNMENT_OBSERVED_LEVELS_ROLE
            ),
            "claim_guardrail": (
                "M2 effect/p are producer-derived; terminal validation independently "
                "reclassifies the frozen statistics."
            ),
        },
        "focal_partner_truth_matrix": [
            {
                "focal_truth_label": focal_truth,
                "partner_truth_label": partner_truth,
                "n_all_pair_rows": int(focal_truth == "TP" and partner_truth == "TP") * 2,
                "n_g1_formal_pair_rows": (
                    int(with_candidate)
                    if focal_truth == "TP" and partner_truth == "TP"
                    else 0
                ),
            }
            for focal_truth in MODULE.TRUTH_LABELS
            for partner_truth in MODULE.TRUTH_LABELS
        ],
        "technical_replication": {
            "status": "ANY_CONCORDANT_EXACT_PAIR_OBSERVED",
            "numerator": 1,
            "denominator": 1,
            "ratio": 1.0,
            "not_evaluable_one_platform_only": 1,
            "denominator_definition": "exact shared focal and partner opportunities",
            "biological_n": 1,
            "independent_biological_replication_n": 0,
            "replication_claim_status": "NOT_EVALUABLE_BIOLOGICAL_N1",
            "inferential_confidence_interval": None,
            "pair_independence_assumption_met": False,
            "required_for_b1": False,
        },
        "biological_replication": {
            "status": "NOT_RUN",
            "independent_biological_replication_n": 0,
        },
        "counts": {
            "cooccurrence_pair_rows": 2,
            "g2_candidates": int(with_candidate),
            "candidate_catalog_rows": int(with_candidate),
            "candidate_witness_pair_rows": int(with_candidate),
            "claim_status_counts": pooled_statuses,
        },
    }
    write_json(final_path, final)
    candidate_catalog_path = root / "candidate_catalog.tsv"
    candidate_witness_pairs_path = root / "candidate_witness_pairs.tsv"
    write_final_tsv(
        candidate_catalog_path,
        candidate_catalog,
        MODULE.FINAL_CANDIDATE_TSV_REQUIRED_FIELDS,
    )
    write_final_tsv(
        candidate_witness_pairs_path,
        candidate_witness_pairs,
        MODULE.FINAL_WITNESS_TSV_REQUIRED_FIELDS,
    )
    final_receipt_path = root / "final_report_dataset.run_receipt.json"
    write_json(
        final_receipt_path,
        {
            "schema_name": MODULE.FINAL_RECEIPT_SCHEMA,
            "schema_version": MODULE.FINAL_RECEIPT_SCHEMA_VERSION,
            "pass": True,
            "task_type": "B_comprehensive_validation",
            "pass_semantics": MODULE.PASS_SEMANTICS,
            "command": [
                "python",
                str(final_builder_path.resolve()),
                "--output-dir",
                str(final_path.parent.resolve()),
            ],
            "outputs": {
                "final_report_dataset": artifact(final_path),
                "candidate_catalog": artifact(candidate_catalog_path),
                "candidate_witness_pairs": artifact(candidate_witness_pairs_path),
            },
            "inputs": final["input_artifacts"],
            "code": {
                "final_report_dataset_builder": artifact(final_builder_path),
            },
        },
    )

    history_path = None
    if with_history:
        history_path = root / "earlier_fp_only_report.md"
        history_path.write_text("# Earlier FP-only specificity analysis\n", encoding="utf-8")

    return (
        MODULE.ReportInputs(
            final_dataset=final_path,
            final_receipt=final_receipt_path,
            candidate_catalog=candidate_catalog_path,
            candidate_witness_pairs=candidate_witness_pairs_path,
            claim_contract=claim_contract_path,
            manifest=manifest_path,
            screen_summary=screen_path,
            cooccurrence_sites=site_path,
            cooccurrence_pairs=pair_path,
            cooccurrence_summary=cooccurrence_summary_path,
            output_reconciliation=reconciliation_path,
            post_immutability_audit=immutability_path,
            tree_input_audit=tree_path,
            reference_identity_audit=reference_audit_path,
            earlier_fp_report=history_path,
        ),
        repo_root,
    )


def build_report(
    tmp_path: Path,
    *,
    with_candidate: bool,
    with_history: bool = False,
    eleven_group_site: bool = False,
) -> tuple[dict[str, Path], MODULE.ReportInputs, Path]:
    inputs, repo_root = build_fixture(
        tmp_path,
        with_candidate=with_candidate,
        with_history=with_history,
        eleven_group_site=eleven_group_site,
    )
    output_dir = tmp_path / "formal-output"
    outputs = MODULE.build_outputs(
        inputs,
        output_dir,
        repo_root=repo_root,
        generated_at="2026-07-15T00:00:00+00:00",
    )
    return outputs, inputs, repo_root


def test_nonzero_candidate_builds_three_canonical_outputs(tmp_path: Path) -> None:
    outputs, inputs, _ = build_report(
        tmp_path,
        with_candidate=True,
        with_history=True,
    )
    assert set(outputs) == {
        "report.md",
        "artifact.json",
        "report_build_receipt.json",
    }
    assert {path.name for path in outputs["report.md"].parent.iterdir()} == {
        "report.md",
        "artifact.json",
        "report_build_receipt.json",
    }

    report = outputs["report.md"].read_text(encoding="utf-8")
    artifact_payload = read_json(outputs["artifact.json"])
    build_receipt = read_json(outputs["report_build_receipt.json"])
    manifest = artifact_payload["manifest"]
    datasets = artifact_payload["snapshot"]["datasets"]
    assert "目前最強的分子證據結論" in report
    assert "multi-marker molecular-haplotype base candidate" in report
    assert "先把所有 M2 exact-testable focal-partner pairs 納入同一全域 BY family" in report
    assert "Permutation 可執行性不是 global FDR family membership" in report
    assert "共同 ancestral ALT 模型的分子預測可被檢驗" in report
    assert "不表示共同 ancestral ALT 已被直接觀察" in report
    assert "屬結構性 `NOT_EVALUABLE`；不是 CN/CCF 否證 0 個候選" in report
    assert "沒有產生或覆寫原本 tagged BAM" in report
    assert "在同一 background payload 上" in report
    assert "這不比較 ALT 與 REF 的實際 flag 集合" in report
    assert "只可能保守減少候選" in report
    assert "仍保守列為可評估的 M2 FAIL" in report
    assert "HP-exact/family/strand=7/5/2" in report
    assert "specificity，不是全 sSNV prevalence" in report
    assert str(inputs.final_dataset.resolve()) in report
    assert artifact_payload["surface"] == "report"
    assert isinstance(datasets, dict)
    assert all(isinstance(rows, list) for rows in datasets.values())
    assert sum(map(len, datasets.values())) <= MODULE.MAX_SNAPSHOT_ROWS
    assert len(datasets["claim_ladder"]) == 9
    assert len(datasets["overview_case_chart"]) == 12
    assert {row["panel"] for row in datasets["overview_case_chart"]} == {
        "claim ladder",
        "actual case",
    }
    assert datasets["case_summary"][0]["selection_mode"] == "G2_PASS_BASE_CANDIDATE"
    assert datasets["case_summary"][0]["negative_result"] is False
    registry_by_role = {
        row["view_role"]: row for row in datasets["case_selection_registry"]
    }
    canonical_view = registry_by_role["canonical_pre_registered_oracle"]
    assert canonical_view["selection_status"] == "TARGET_UNAVAILABLE_NO_SUBSTITUTION"
    assert canonical_view["selected_focal_site"] == "N/A"
    assert canonical_view["selected_partner_site"] == "N/A"
    assert registry_by_role["evaluable_statistical_negative"]["negative_result"] is True
    view_by_role = {row["view_role"]: row for row in datasets["case_view_summary"]}
    statistical_negative = view_by_role["evaluable_statistical_negative"]
    assert statistical_negative["conditional_permutations"] == 999
    assert statistical_negative["conditional_sensitivity_pass"] is False
    assert statistical_negative["negative_result_scope"].startswith(
        "evaluated conditional endpoint FAIL"
    )
    assert view_by_role["primary_report_case"]["normal_partner_status"] == (
        "NOT_EVALUATED_BY_DESIGN"
    )
    assert len(datasets["case_four_state_models"]) == 3
    assert manifest["blocks"][0]["body"] == f"# {manifest['title']}"
    assert manifest["charts"]
    assert any(block["type"] == "chart" for block in manifest["blocks"])
    assert "|---|" in report
    assert all(
        "|---|" not in str(block.get("body", ""))
        for block in manifest["blocks"]
        if block["type"] == "markdown"
    )
    assert build_receipt["schema_name"] == MODULE.REPORT_RECEIPT_SCHEMA
    assert build_receipt["schema_version"] == MODULE.REPORT_RECEIPT_SCHEMA_VERSION
    assert build_receipt["pass"] is True
    assert build_receipt["snapshot_rows"] == sum(map(len, datasets.values()))
    assert build_receipt["claim_contract"]["sha256"] == MODULE.sha256(
        inputs.claim_contract
    )
    assert build_receipt["outputs"]["report_md"]["sha256"] == MODULE.sha256(
        outputs["report.md"]
    )
    assert build_receipt["outputs"]["artifact_json"]["sha256"] == MODULE.sha256(
        outputs["artifact.json"]
    )
    assert build_receipt["artifact_presentation_scope"].startswith(
        "complete_portable_technical_report"
    )
    inventory = {row["role"]: row for row in datasets["source_inventory"]}
    assert inventory[MODULE.CLAIM_CONTRACT_SOURCE_ROLE]["sha256"] == MODULE.sha256(
        inputs.claim_contract
    )
    assert "final dataset producer receipt" in inventory
    assert "extraction reference identity audit" in inventory
    assert "current terminal screen receipt" in inventory
    assert "stable-primary pre-consumer audit" in inventory
    assert "stable-primary post-consumer audit" in inventory
    assert all(
        f"producer BAM-output receipt {sample}" in inventory
        for sample in MODULE.EXPECTED_DATASETS
    )
    assert len(datasets["producer_bam_output"]) == 7
    assert all(row["persisted_bam"] is False for row in datasets["producer_bam_output"])
    for table in manifest["tables"]:
        sort = table["defaultSort"]
        assert sort["direction"] in {"asc", "desc"}
        assert sort["field"] in {column["field"] for column in table["columns"]}
        assert len(table["columns"]) <= 8
    serialized_artifact = json.dumps(artifact_payload, ensure_ascii=False)
    assert str(tmp_path.resolve()) not in serialized_artifact
    for source in artifact_payload["sources"]:
        assert not Path(source["path"]).is_absolute()
        assert ".." not in Path(source["path"]).parts
        assert all(not Path(path).is_absolute() for path in source["query"]["tables_used"])
        assert all(".." not in Path(path).parts for path in source["query"]["tables_used"])


def test_report_exposes_eleven_group_m1_pass_m2_not_evaluable_boundary(
    tmp_path: Path,
) -> None:
    outputs, _, _ = build_report(
        tmp_path,
        with_candidate=True,
        eleven_group_site=True,
    )
    report = outputs["report.md"].read_text(encoding="utf-8")
    artifact_payload = read_json(outputs["artifact.json"])
    datasets = artifact_payload["snapshot"]["datasets"]
    assert "觀察 2-10 個 methyl groups" in report
    assert "本次超過上限者為 1 個" in report
    assert "HCC1954 chr5:751076 G>A (11 groups)" in report
    m2_rows = datasets["m2_evaluability"]
    assert m2_rows[0]["count"] == 1
    assert m2_rows[1] == {
        "record_type": "observed_example",
        "site": "HCC1954 chr5:751076 G>A",
        "count": 1,
        "observed_methyl_groups": 11,
        "maximum_supported_methyl_groups": 10,
        "m1_status": "PASS",
        "m2_status": "NOT_EVALUABLE",
        "g1_status": "NOT_RUN",
        "g2_status": "NOT_RUN",
        "b1_status": "NOT_RUN",
        "claim_behavior": "M1 PASS; M2 NE; G1/G2/B1 NOT_RUN",
        "reason": (
            "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_"
            "MAXIMUM:observed=11:maximum=10"
        ),
    }
    m2_claim = next(
        row for row in datasets["claim_ladder"] if row["claim_id"] == "M2"
    )
    assert m2_claim["dataset_not_evaluable"] == 1
    assert m2_claim["denominator_definition"] == MODULE.M2_DENOMINATOR_DEFINITION
    hcc1954 = [
        row
        for row in datasets["stratum_claim_metrics"]
        if row["stratum_type"] == "sample" and row["stratum"] == "HCC1954"
    ]
    by_claim = {row["claim_id"]: row for row in hcc1954}
    assert by_claim["M2"]["not_evaluable"] == 1
    hcc1954_population = sum(SAMPLE_TRUTH_COUNTS["HCC1954"].values())
    assert by_claim["G1"]["not_run"] == hcc1954_population
    assert by_claim["G2"]["not_run"] == hcc1954_population
    assert any(
        block.get("tableId") == "table-m2-evaluability"
        for block in artifact_payload["manifest"]["blocks"]
    )


def test_zero_candidate_uses_ranked_nonconfirming_fallback(tmp_path: Path) -> None:
    outputs, _, _ = build_report(tmp_path, with_candidate=False)
    report = outputs["report.md"].read_text(encoding="utf-8")
    artifact_payload = read_json(outputs["artifact.json"])
    case = artifact_payload["snapshot"]["datasets"]["case_summary"][0]
    assert case["selection_mode"] == "NON_CONFIRMING_FALLBACK"
    assert case["partner_site"] == "chr1:150 G>T"
    assert case["exact_p"] == pytest.approx(0.02)
    assert case["cramers_v"] == pytest.approx(0.80)
    assert case["g2_status"] == "FAIL"
    assert "未確認個案（non-confirming）" in report
    assert "不得計為候選" in report
    assert "甲基異質性與局部共分離" in report
    assert artifact_payload["snapshot"]["datasets"]["claim_ladder"][3][
        "dataset_numerator"
    ] == 0


def test_zero_g1_and_g2_fallback_does_not_invent_local_cosegregation() -> None:
    answer = MODULE.biological_answer_for_claims(
        g1={"dataset_numerator": 0},
        g2={"dataset_numerator": 0},
    )
    assert "G1 也為 0" in answer
    assert "沒有 read-level partner co-segregation 證據" in answer
    assert "甲基異質性與局部共分離" not in answer


def test_zero_sortable_pair_table_is_legal_structured_not_applicable(
    tmp_path: Path,
) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=False)
    clear_tsv_rows(inputs.cooccurrence_pairs)
    refresh_final_input_artifact(inputs, "cooccurrence_pairs", inputs.cooccurrence_pairs)

    outputs = MODULE.build_outputs(
        inputs,
        tmp_path / "formal-output",
        repo_root=repo_root,
        generated_at="2026-07-15T00:00:00+00:00",
    )
    report = outputs["report.md"].read_text(encoding="utf-8")
    artifact = read_json(outputs["artifact.json"])
    datasets = artifact["snapshot"]["datasets"]
    case = datasets["case_summary"][0]
    assert case["selection_mode"] == "NOT_APPLICABLE_NO_ELIGIBLE_PAIR"
    assert case["evidence_status"] == "NOT_APPLICABLE"
    assert case["negative_result"] is False
    assert case["focal_site"] == "N/A"
    assert datasets["case_group_partner"] == [
        {
            "row_order": 1,
            "methyl_group": "N/A",
            "partner_call": "N/A",
            "count": None,
            "status": "NOT_APPLICABLE_NO_ELIGIBLE_PAIR",
        }
    ]
    assert len(datasets["overview_case_chart"]) == 8
    assert {row["panel"] for row in datasets["overview_case_chart"]} == {
        "claim ladder"
    }
    assert all(
        row["negative_result"] is False
        for row in datasets["case_selection_registry"]
    )
    assert "不是 negative biological result" in report
    portable_summary = next(
        block["body"]
        for block in artifact["manifest"]["blocks"]
        if block["id"] == "block-portable-summary"
    )
    assert "actual case=N/A" in portable_summary
    assert "chart-case-group-partner" not in {
        chart["id"] for chart in artifact["manifest"]["charts"]
    }


def test_four_view_registry_prefers_pre_registered_oracle_independently(
    tmp_path: Path,
) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=False)
    rows = read_tsv(inputs.cooccurrence_pairs)
    canonical = pair_row(
        750_411,
        exact_p=0.03,
        effect=0.65,
        delta=0.55,
        formal=False,
    )
    canonical.update(
        {
            "sample": "HCC1395_DORADO",
            "focal_chrom": "chr5",
            "focal_pos": 750_311,
            "focal_ref": "C",
            "focal_alt": "T",
            "partner_chrom": "chr5",
        }
    )
    rows.append({field: str(value) for field, value in canonical.items()})
    write_tsv(inputs.cooccurrence_pairs, rows)
    site_rows = read_tsv(inputs.cooccurrence_sites)
    canonical_site = dict(site_rows[0])
    canonical_site.update(
        {
            "sample": "HCC1395_DORADO",
            "chrom": "chr5",
            "pos": "750311",
            "ref": "C",
            "alt": "T",
            "component_id": "component-canonical",
            "multi_marker_molecular_haplotype_base_candidate": "false",
            "joint_signature_sensitivity_pass": "false",
        }
    )
    site_rows.append(canonical_site)
    write_tsv(inputs.cooccurrence_sites, site_rows)
    refresh_final_input_artifact(inputs, "cooccurrence_pairs", inputs.cooccurrence_pairs)
    refresh_final_input_artifact(inputs, "cooccurrence_sites", inputs.cooccurrence_sites)

    outputs = MODULE.build_outputs(
        inputs,
        tmp_path / "formal-output",
        repo_root=repo_root,
        generated_at="2026-07-15T00:00:00+00:00",
    )
    registry = read_json(outputs["artifact.json"])["snapshot"]["datasets"][
        "case_selection_registry"
    ]
    by_role = {row["view_role"]: row for row in registry}
    assert set(by_role) == {
        "aggregate",
        "canonical_pre_registered_oracle",
        "extreme_global_exact_effect",
        "well_explained",
        "evaluable_statistical_negative",
    }
    assert len({row["selection_definition"] for row in registry}) == 5
    assert by_role["aggregate"]["selected_focal_site"] == "N/A"
    oracle = by_role["canonical_pre_registered_oracle"]
    assert oracle["selection_status"] == "SELECTED_PRE_REGISTERED_ORACLE"
    assert oracle["selected_focal_site"] == "HCC1395_DORADO chr5:750311 C>T"
    assert by_role["extreme_global_exact_effect"]["selected_focal_site"].startswith(
        "HCC1395 chr1:100"
    )


def test_artifact_has_complete_sections_and_exact_source_mapping(tmp_path: Path) -> None:
    outputs, _, _ = build_report(tmp_path, with_candidate=True)
    report = outputs["report.md"].read_text(encoding="utf-8")
    artifact_payload = read_json(outputs["artifact.json"])
    manifest = artifact_payload["manifest"]
    tables = {table["id"]: table for table in manifest["tables"]}
    charts = {chart["id"]: chart for chart in manifest["charts"]}
    blocks = {block["id"]: block for block in manifest["blocks"]}
    sources = {source["id"]: source for source in artifact_payload["sources"]}

    assert tables["table-flow"]["sourceId"] == "src-flow-view"
    assert tables["table-audits"]["sourceId"] == "src-audit-view"
    assert tables["table-case-selection-registry"]["sourceId"] == (
        "src-case-selection-view"
    )
    assert tables["table-case-summary"]["sourceId"] == "src-case-summary-view"
    assert tables["table-case-group-partner"]["sourceId"] == "src-case-group-view"
    assert tables["table-case-four-state"]["sourceId"] == "src-case-four-state-view"
    assert tables["table-case-joint-signature"]["sourceId"] == "src-case-joint-view"
    assert charts["chart-case-group-partner"]["sourceId"] == "src-case-group-view"
    assert "non-confirming witness" in charts["chart-case-group-partner"]["subtitle"]
    assert "cellular clone" in charts["chart-case-group-partner"]["subtitle"]
    assert charts["chart-overview-case"]["sourceId"] == "src-overview-case-view"
    assert "M1=FLAGGED/all dataset-sites" in charts["chart-overview-case"]["subtitle"]
    assert (
        "M2-G2=PASS/claim-specific evaluable denominator"
        in charts["chart-overview-case"]["subtitle"]
    )
    assert "partner R/A read" in charts["chart-overview-case"]["subtitle"]
    overview_rows = artifact_payload["snapshot"]["datasets"]["overview_case_chart"]
    m1_rows = [
        row
        for row in overview_rows
        if row["panel"] == "claim ladder" and row["category"] == "M1"
    ]
    assert {row["series"] for row in m1_rows} == {
        "M1 FLAGGED / all dataset-sites",
        "M1 descriptive biological-site aggregation",
    }
    assert all("not biological prevalence" in row["interpretation"] for row in m1_rows)
    overview_metric_definitions = sources["src-overview-case-view"]["query"][
        "metric_definitions"
    ]
    assert "M1 dataset-site proportion = FLAGGED / all dataset-sites" in (
        overview_metric_definitions
    )
    assert "M2-G2 claim proportion = PASS / claim-specific evaluable denominator" in (
        overview_metric_definitions
    )
    assert all(
        "claim proportion = PASS / evaluable" != definition
        for definition in overview_metric_definitions
    )
    assert len(sources["src-flow-view"]["query"]["tables_used"]) == 9
    assert len(sources["src-audit-view"]["query"]["tables_used"]) == 12
    assert len(sources["src-case-summary-view"]["query"]["tables_used"]) == 2
    assert len(sources["src-case-group-view"]["query"]["tables_used"]) == 1
    assert len(sources["src-overview-case-view"]["query"]["tables_used"]) == 2
    assert len(sources["src-case-four-state-view"]["query"]["tables_used"]) == 1
    assert len(sources["src-case-joint-view"]["query"]["tables_used"]) == 1

    required_blocks = {
        "block-title",
        "block-portable-summary",
        "block-conclusion",
        "block-overview-case-chart",
        "block-claim-table",
        "block-strata-table",
        "block-audit-table",
        "block-case-selection",
        "block-case-view-summary",
        "block-case-four-state-models",
        "block-case-joint",
        "block-ancestral-hypothesis",
        "block-methods",
        "block-limits",
        "block-sources",
        "block-source-table",
    }
    assert required_blocks.issubset(blocks)
    assert len(blocks) >= 30
    assert blocks["block-title"]["body"] == f"# {MODULE.PORTABLE_TITLE}"
    assert "469,849" in blocks["block-portable-summary"]["body"]
    assert "cellular subclone" in blocks["block-portable-summary"]["body"]
    assert blocks["block-overview-case-chart"]["chartId"] == "chart-overview-case"
    assert "外部方法對照" in report
    assert "建議的下一步" in report
    assert "尚待回答的問題" in report
    assert "dataset-site union" in report
    assert "unique biological loci" in report


@pytest.mark.parametrize(
    ("mutation", "pattern"),
    (
        (lambda receipt: receipt.update(pass_=False), "did not PASS"),
        (
            lambda receipt: receipt.update(schema_version="9.9.9"),
            "schema version",
        ),
        (
            lambda receipt: receipt.update(task_type="A_exploratory_pilot"),
            "task type",
        ),
        (
            lambda receipt: receipt.update(pass_semantics="scientific_confirmation"),
            "pass semantics",
        ),
        (lambda receipt: receipt.update(command=[]), "command must be a non-empty"),
        (
            lambda receipt: receipt["outputs"]["final_report_dataset"].update(
                sha256="0" * 64
            ),
            "SHA-256 mismatch",
        ),
        (
            lambda receipt: receipt["code"]["final_report_dataset_builder"].update(
                sha256="0" * 64
            ),
            "SHA-256 mismatch",
        ),
    ),
)
def test_final_receipt_contract_fails_closed(
    tmp_path: Path, mutation, pattern: str
) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    receipt = read_json(inputs.final_receipt)
    mutation(receipt)
    if "pass_" in receipt:
        receipt["pass"] = receipt.pop("pass_")
    write_json(inputs.final_receipt, receipt)
    with pytest.raises(MODULE.ReportContractError, match=pattern):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_claim_contract_vocabulary_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    contract = inputs.claim_contract.read_text(encoding="utf-8")
    inputs.claim_contract.write_text(
        contract.replace(MODULE.CLAIM_NAMES["G2"], "removed G2 name", 1),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        MODULE, "TERMINAL_CLAIM_CONTRACT_SHA256", MODULE.sha256(inputs.claim_contract)
    )
    with pytest.raises(MODULE.ReportContractError, match="complete terminal claim vocabulary"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_claim_contract_terminal_semantics_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=False)
    contract = inputs.claim_contract.read_text(encoding="utf-8")
    inputs.claim_contract.write_text(
        contract.replace("bonferroni_three_relation_models", "legacy_single_model"),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        MODULE, "TERMINAL_CLAIM_CONTRACT_SHA256", MODULE.sha256(inputs.claim_contract)
    )
    with pytest.raises(MODULE.ReportContractError, match="missing_fragments"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_claim_contract_any_semantic_drift_fails_exact_sha(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=False)
    contract = inputs.claim_contract.read_text(encoding="utf-8")
    inputs.claim_contract.write_text(
        contract.replace("Cramer's V>=0.30", "Cramer's V>=0.00"),
        encoding="utf-8",
    )
    with pytest.raises(MODULE.ReportContractError, match="SHA-256 mismatch"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_source_attested_release_requires_exact_v5_claim_contract(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    source_receipt = tmp_path / "source_identity_receipt.json"
    source_snapshot = tmp_path / "source_identity_snapshot.json"
    write_json(source_receipt, {"pass": True, "schema_version": "1.1.0"})
    write_json(source_snapshot, {"pass": True, "schema_version": "1.1.0"})
    final = read_json(inputs.final_dataset)
    final["input_artifacts"]["tumor_ref_source_identity_receipt"] = artifact(
        source_receipt
    )
    final["tumor_ref_source_identity_attestation"] = {
        "status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
        "release_gate_pass": True,
        "publishable_task_b_release": True,
        "audit_class": "bounded_retrospective_source_file_identity",
        "receipt": artifact(source_receipt),
        "snapshot": artifact(source_snapshot),
        "source_roles": ["analyzer", "focal_alt_cluster_lib"],
        "source_sha256": {
            "analyzer": "a" * 64,
            "focal_alt_cluster_lib": "b" * 64,
        },
        "limitation": (
            "Retrospective bounded source-file attestation only; not a prelaunch lock or a "
            "complete environment attestation."
        ),
    }
    final["m2_evaluability_contract"]["independent_logic_audit"] = {
        "status": "PASS_LOGIC_INDEPENDENT_RECOUNT",
        "pass": True,
        "production_gate_imported": False,
        "production_gate_functions_called": False,
        "counts": {
            "all_rows": MODULE.EXPECTED_SITES,
            "m1_stable_rows": sum(M1_PASS.values()),
            "eligible": sum(M2_PASS.values()),
            "evaluable_ineligible": sum(M1_PASS.values()) - sum(M2_PASS.values()),
            "not_evaluable_axis_indeterminate": 0,
            "not_evaluable_group_count_gt10": 0,
        },
        "n_evaluable_sites_with_aligned_axis_below_negative_evaluability_power": 0,
    }
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)

    with pytest.raises(
        MODULE.ReportContractError,
        match="Claim contract SHA-256 mismatch for release state",
    ):
        MODULE.build_outputs(
            inputs,
            tmp_path / "release-with-v3",
            repo_root=repo_root,
        )
    assert not (tmp_path / "release-with-v3").exists()

    inputs.claim_contract.write_text(
        (SCRIPT.parents[1] / "claim-contract-v5.md").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    outputs = MODULE.build_outputs(
        inputs,
        tmp_path / "release-with-v5",
        repo_root=repo_root,
    )
    receipt = read_json(outputs["report_build_receipt.json"])
    assert receipt["validations"]["claim_contract_release_state_semantics"] == (
        "PASS_RELEASE_V5"
    )
    assert receipt["validations"][
        "tumor_ref_bounded_source_identity_release_gate"
    ] == "PASS"
    assert receipt["claim_contract"]["sha256"] == MODULE.RELEASE_CLAIM_CONTRACT_SHA256


def test_final_candidate_tsv_hash_tamper_fails_closed(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    inputs.candidate_catalog.write_text(
        inputs.candidate_catalog.read_text(encoding="utf-8").replace("0.04", "0.05"),
        encoding="utf-8",
    )
    with pytest.raises(MODULE.ReportContractError, match="candidate_catalog.*SHA-256"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_final_candidate_tsv_resigned_semantic_drift_fails_closed(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    rows = read_tsv(inputs.candidate_catalog)
    rows[0]["joint_signature_q_global_by"] = "0.05"
    write_final_tsv(
        inputs.candidate_catalog,
        rows,
        MODULE.FINAL_CANDIDATE_TSV_REQUIRED_FIELDS,
    )
    receipt = read_json(inputs.final_receipt)
    receipt["outputs"]["candidate_catalog"] = artifact(inputs.candidate_catalog)
    write_json(inputs.final_receipt, receipt)
    with pytest.raises(MODULE.ReportContractError, match="TSV/JSON drift"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_final_witness_tsv_resigned_semantic_drift_fails_closed(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    rows = read_tsv(inputs.candidate_witness_pairs)
    rows[0]["four_state_relation_family_size"] = "1"
    write_final_tsv(
        inputs.candidate_witness_pairs,
        rows,
        MODULE.FINAL_WITNESS_TSV_REQUIRED_FIELDS,
    )
    receipt = read_json(inputs.final_receipt)
    receipt["outputs"]["candidate_witness_pairs"] = artifact(
        inputs.candidate_witness_pairs
    )
    write_json(inputs.final_receipt, receipt)
    with pytest.raises(MODULE.ReportContractError, match="TSV/JSON drift"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_cli_requires_final_receipt_and_claim_contract() -> None:
    actions = {action.dest: action for action in MODULE.build_parser()._actions}
    assert actions["final_receipt"].required is True
    assert actions["final_release_receipt"].required is True
    assert actions["final_release_signature"].required is True
    assert actions["final_release_public_key"].required is True
    assert actions["candidate_catalog"].required is True
    assert actions["candidate_witness_pairs"].required is True
    assert actions["claim_contract"].required is True
    assert actions["reference_identity_audit"].required is True


@pytest.mark.parametrize(
    "forbidden",
    (
        "subclone-compatible",
        "independent marker",
        "C1-C5",
        "C2 stable methyl",
        "C3 co-segregation",
        "C4 compatible clone",
        "C5 lineage",
    ),
)
def test_generated_outputs_exclude_legacy_terms(tmp_path: Path, forbidden: str) -> None:
    outputs, _, _ = build_report(tmp_path, with_candidate=True)
    rendered = "\n".join(path.read_text(encoding="utf-8") for path in outputs.values())
    assert forbidden.lower() not in rendered.lower()


def test_current_terminal_tag_audit_fails_closed(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    screen = read_json(inputs.screen_summary)
    screen["latest_hp_ps_terminal_join_audit"]["pass"] = False
    write_json(inputs.screen_summary, screen)
    refresh_final_input_artifact(inputs, "screen_summary", inputs.screen_summary)
    with pytest.raises(MODULE.ReportContractError, match="terminal HP/PS audit failed closed"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_post_immutability_must_follow_terminal_screen(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    immutability = read_json(inputs.post_immutability_audit)
    immutability["created_at_utc"] = "2026-07-14T23:59:59+00:00"
    write_json(inputs.post_immutability_audit, immutability)
    with pytest.raises(
        MODULE.ReportContractError,
        match="predates the selected terminal screen",
    ):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)


def test_final_receipt_input_inventory_must_match_final_dataset(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    receipt = read_json(inputs.final_receipt)
    receipt["inputs"].pop("primary_artifact_audit_post")
    write_json(inputs.final_receipt, receipt)
    with pytest.raises(
        MODULE.ReportContractError,
        match="input_artifacts differ",
    ):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)


def test_persisted_producer_bam_policy_fails_closed(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    manifest = read_json(inputs.manifest)
    layered_path = Path(manifest["layered_root"]["path"])
    layered = read_json(layered_path)
    declared = layered["samples"][0]["read_tags"]["producer_capture_receipt"]
    receipt_path = Path(declared["path"])
    receipt = read_json(receipt_path)
    receipt["bam_output_policy"]["persisted_bam"] = True
    write_json(receipt_path, receipt)
    receipt_identity = artifact(receipt_path)
    declared["identity"].update(
        {
            "size_bytes": receipt_identity["size_bytes"],
            "sha256": receipt_identity["sha256"],
        }
    )
    write_json(layered_path, layered)
    manifest["layered_root"] = artifact(layered_path)
    write_json(inputs.manifest, manifest)
    refresh_final_input_artifact(inputs, "manifest", inputs.manifest)

    with pytest.raises(MODULE.ReportContractError, match="Producer BAM-output policy failed"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_ratio_and_denominator_reconciliation_fail_before_output(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    final["funnel_metrics"]["pooled"]["M1"]["ratio"] = 0.12345
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)
    with pytest.raises(MODULE.ReportContractError, match="Ratio does not reconcile"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_m1_all_site_operational_denominator_fails_closed_on_drift(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    final["m1_operational_screen"]["n_all_dataset_sites"] -= 1
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)
    with pytest.raises(
        MODULE.ReportContractError,
        match="M1 screen evaluability counts do not reconcile",
    ):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_background_control_gate_disclosure_fails_closed_on_ari_drift(
    tmp_path: Path,
) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    final["background_control_replication_gate"][
        "membership_ari_minimum_required"
    ] = True
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)
    with pytest.raises(
        MODULE.ReportContractError,
        match="must disclose that membership ARI is not required",
    ):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_report_validates_verified_tumor_ref_source_identity_release_gate(
    tmp_path: Path,
) -> None:
    receipt = tmp_path / "source_identity_receipt.json"
    snapshot = tmp_path / "source_identity_snapshot.json"
    write_json(receipt, {"pass": True})
    write_json(snapshot, {"pass": True})
    final = {
        "tumor_ref_source_identity_attestation": {
            "status": "VERIFIED_BOUNDED_RETROSPECTIVE_SOURCE_IDENTITY",
            "release_gate_pass": True,
            "publishable_task_b_release": True,
            "audit_class": "bounded_retrospective_source_file_identity",
            "receipt": artifact(receipt),
            "snapshot": artifact(snapshot),
            "source_roles": ["analyzer", "focal_alt_cluster_lib"],
            "source_sha256": {
                "analyzer": "a" * 64,
                "focal_alt_cluster_lib": "b" * 64,
            },
            "limitation": (
                "Retrospective bounded source-file attestation only; not a prelaunch lock or a "
                "complete environment attestation."
            ),
        }
    }
    observed = MODULE.validate_tumor_ref_source_identity_attestation(final)
    assert observed["release_gate_pass"] is True
    assert observed["_receipt_path"] == receipt.resolve()
    assert observed["_snapshot_path"] == snapshot.resolve()

    final["tumor_ref_source_identity_attestation"]["source_sha256"]["analyzer"] = "bad"
    with pytest.raises(MODULE.ReportContractError, match="SHA is malformed"):
        MODULE.validate_tumor_ref_source_identity_attestation(final)


def test_report_rejects_candidate_without_global_by_confirmation(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    final["candidate_catalog"][0]["joint_signature_q_global_by"] = 0.90
    write_json(inputs.final_dataset, final)
    candidate_rows = read_tsv(inputs.candidate_catalog)
    candidate_rows[0]["joint_signature_q_global_by"] = "0.9"
    write_final_tsv(
        inputs.candidate_catalog,
        candidate_rows,
        MODULE.FINAL_CANDIDATE_TSV_REQUIRED_FIELDS,
    )
    refresh_final_receipt(inputs)
    receipt = read_json(inputs.final_receipt)
    receipt["outputs"]["candidate_catalog"] = artifact(inputs.candidate_catalog)
    write_json(inputs.final_receipt, receipt)
    with pytest.raises(MODULE.ReportContractError, match="Candidate G2 global-BY gate drift"):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_report_rejects_false_raw_read_recomputation_claim(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    final["m2_axis_statistic_provenance"][
        "downstream_raw_read_axis_statistic_recomputed"
    ] = True
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)
    with pytest.raises(
        MODULE.ReportContractError, match="M2 axis-statistic provenance contract drift"
    ):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_report_requires_complete_asymmetric_m2_downstream_checks(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    final["m2_axis_statistic_provenance"]["downstream_checks"].remove(
        "asymmetric positive-confound versus negative-evaluability power decision"
    )
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)

    with pytest.raises(
        MODULE.ReportContractError, match="M2 axis-statistic provenance contract drift"
    ):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_report_rejects_categorical_planning_level_drift(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    final["m2_evaluability_contract"]["categorical_planning_level_ceilings"][
        "hp_exact"
    ] = 2
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)

    with pytest.raises(
        MODULE.ReportContractError, match="M2 categorical planning-level contract drift"
    ):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


@pytest.mark.parametrize(
    ("field", "pattern"),
    (
        ("output_reconciliation", "output reconciliation did not PASS"),
        ("post_immutability_audit", "post-run immutability audit did not PASS"),
        ("tree_input_audit", "tree-input audit did not PASS"),
        (
            "reference_identity_audit",
            "extraction reference identity audit did not PASS",
        ),
    ),
)
def test_release_blocking_audits_fail_closed(
    tmp_path: Path,
    field: str,
    pattern: str,
) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    audit_path = getattr(inputs, field)
    payload = read_json(audit_path)
    payload["pass"] = False
    write_json(audit_path, payload)
    with pytest.raises(MODULE.ReportContractError, match=pattern):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()


def test_existing_output_directory_is_never_overwritten(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    output_dir = tmp_path / "formal-output"
    output_dir.mkdir()
    sentinel = output_dir / "sentinel.txt"
    sentinel.write_text("keep\n", encoding="utf-8")
    with pytest.raises(MODULE.ReportContractError, match="overwrite prohibited"):
        MODULE.build_outputs(inputs, output_dir, repo_root=repo_root)
    assert sentinel.read_text(encoding="utf-8") == "keep\n"


def test_l1_l2_cannot_pass_without_orthogonal_evidence(tmp_path: Path) -> None:
    inputs, repo_root = build_fixture(tmp_path, with_candidate=True)
    final = read_json(inputs.final_dataset)
    l1 = final["funnel_metrics"]["pooled"]["L1"]
    l1.update(
        {
            "numerator": 1,
            "denominator": 1,
            "ratio": 1.0,
            "not_evaluable": 0,
            "not_run": MODULE.EXPECTED_SITES - 1,
        }
    )
    final["counts"]["claim_status_counts"]["L1"] = {
        "PASS": 1,
        "FAIL": 0,
        "NOT_EVALUABLE": 0,
        "NOT_RUN": MODULE.EXPECTED_SITES - 1,
    }
    final["claim_ladder"][7].update(
        {
            "dataset_numerator": 1,
            "dataset_denominator": 1,
            "dataset_ratio": 1.0,
            "dataset_not_evaluable": 0,
            "dataset_not_run": MODULE.EXPECTED_SITES - 1,
            "status": "PASS",
        }
    )
    write_json(inputs.final_dataset, final)
    refresh_final_receipt(inputs)
    with pytest.raises(MODULE.ReportContractError):
        MODULE.build_outputs(inputs, tmp_path / "formal-output", repo_root=repo_root)
    assert not (tmp_path / "formal-output").exists()
