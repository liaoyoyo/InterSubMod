#!/usr/bin/env python3
"""Build the single source dataset used by the final Markdown, figures, and HTML."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scipy.stats import fisher_exact


TOPIC_ROOT = Path(__file__).resolve().parents[1]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def parse_json_cell(value: str, default: Any) -> Any:
    if value in {"", "NA", "None"}:
        return default
    return json.loads(value)


def parse_float(value: str) -> float | None:
    if value in {"", "NA", "None"}:
        return None
    return float(value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def endpoint(row: dict[str, str], name: str) -> bool:
    if name == "forced_silhouette_multigroup":
        return row.get("forced_k", "") not in {"", "0", "0.0", "NA", "None"}
    return parse_bool(row.get(name, ""))


def classify_hp(label: str) -> str:
    if label in {"1-1", "2-1"}:
        return "first_somatic_haplotype"
    if label in {"1-2", "2-2"}:
        return "second_somatic_haplotype"
    if label == "3":
        return "ambiguous"
    return "untagged_or_other"


def aggregate_hp(rows: list[dict[str, str]]) -> dict[str, Any]:
    exact: Counter[str] = Counter()
    family: Counter[str] = Counter()
    expected = 0
    for row in rows:
        counts = parse_json_cell(row.get("alt_hp_counts", ""), {})
        for label, count in counts.items():
            count = int(count)
            exact[str(label)] += count
            family[classify_hp(str(label))] += count
        expected += int(row.get("n_alt_raw", 0) or 0)
    observed = sum(exact.values())
    if expected != observed:
        raise RuntimeError(f"ALT HP count mismatch: n_alt_raw={expected}, parsed_tags={observed}")
    return {
        "n_alt_reads": observed,
        "exact_counts": dict(sorted(exact.items())),
        "family_counts": dict(sorted(family.items())),
        "second_somatic_haplotype_fraction": rate(family["second_somatic_haplotype"], observed),
    }


def coverage_bands(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    bands = [
        ("6-9", 6, 9),
        ("10-19", 10, 19),
        ("20-39", 20, 39),
        (">=40", 40, None),
    ]
    result = []
    for label, low, high in bands:
        selected = [
            row
            for row in rows
            if int(row["n_alt_after_peel"]) >= low
            and (high is None or int(row["n_alt_after_peel"]) <= high)
        ]
        stable = sum(endpoint(row, "stable_null_multigroup") for row in selected)
        high_threshold = sum(
            endpoint(row, "phase_anchored_robust_epigenetic_candidate") for row in selected
        )
        result.append(
            {
                "band": label,
                "n_evaluable": len(selected),
                "n_stable": stable,
                "stable_fraction": rate(stable, len(selected)),
                "n_hp_tag_covered_high_threshold": high_threshold,
                "hp_tag_covered_high_threshold_fraction": rate(high_threshold, len(selected)),
            }
        )
    return result


def shared_read_components(assignments_path: Path) -> dict[str, Any]:
    rows = [json.loads(line) for line in assignments_path.read_text(encoding="utf-8").splitlines() if line]
    parent = list(range(len(rows)))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    first_site: dict[tuple[str, str], int] = {}
    read_counts: Counter[tuple[str, str]] = Counter()
    n_assignments = 0
    for index, row in enumerate(rows):
        for read_name in row["read_names"]:
            key = (row["sample"], read_name)
            n_assignments += 1
            read_counts[key] += 1
            if key in first_site:
                union(index, first_site[key])
            else:
                first_site[key] = index

    component_sizes = Counter(find(index) for index in range(len(rows)))
    return {
        "n_stable_sites": len(rows),
        "n_read_assignments": n_assignments,
        "n_unique_dataset_reads": len(read_counts),
        "n_reused_dataset_reads": sum(count > 1 for count in read_counts.values()),
        "n_shared_read_components": len(component_sizes),
        "largest_component_size_sites": max(component_sizes.values()),
        "n_components_with_multiple_sites": sum(size > 1 for size in component_sizes.values()),
    }


def selected_case_specs() -> list[dict[str, Any]]:
    return [
        {"case_id": "replicated_chr1_hcc", "sample": "HCC1395", "chrom": "chr1", "pos": 175089892,
         "role": "technical_replication_positive"},
        {"case_id": "replicated_chr1_dorado", "sample": "HCC1395_DORADO", "chrom": "chr1", "pos": 175089892,
         "role": "technical_replication_positive"},
        {"case_id": "replicated_chr9_hcc", "sample": "HCC1395", "chrom": "chr9", "pos": 39585394,
         "role": "replicated_but_cluster_count_changes"},
        {"case_id": "replicated_chr9_dorado", "sample": "HCC1395_DORADO", "chrom": "chr9", "pos": 39585394,
         "role": "replicated_but_cluster_count_changes"},
        {"case_id": "linear_context_not_lineage", "sample": "HCC1937", "chrom": "chr9", "pos": 23323391,
         "role": "high_threshold_linear_regional_context"},
        {"case_id": "branching_context_same_endpoint", "sample": "HCC1937", "chrom": "chr14", "pos": 58925240,
         "role": "high_threshold_branching_regional_context"},
        {"case_id": "hp_axis_counterexample", "sample": "HCC1937", "chrom": "chr2", "pos": 143925487,
         "role": "hp_axis_confound"},
        {"case_id": "technical_axis_counterexample", "sample": "H2009", "chrom": "chr1", "pos": 244024747,
         "role": "technical_axis_confound"},
        {"case_id": "forced_only_counterexample", "sample": "H2009", "chrom": "chr18", "pos": 60839476,
         "role": "forced_silhouette_but_null_rejected"},
        {"case_id": "strict_branching_candidate", "sample": "H1437", "chrom": "chr1", "pos": 85105364,
         "role": "strict_followup_with_branching_context"},
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=TOPIC_ROOT / "results" / "report_dataset_v1")
    args = parser.parse_args()

    paths = {
        "preflight": TOPIC_ROOT / "results" / "latest_input_preflight.json",
        "fp_summary": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "latest_full_v3_frozen" / "latest_summary.json",
        "fp_sites": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "latest_full_v3_frozen" / "latest_site_results_with_topology.tsv",
        "fp_assignments": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "latest_full_v3_frozen" / "latest_stable_multigroup_read_assignments.jsonl",
        "tp_summary": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "matched_tp_v2_frozen" / "matched_tp_summary.json",
        "comparison": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "fp_matched_tp_comparison_v1" / "fp_matched_tp_comparison_summary.json",
        "robustness": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "fp_matched_tp_robustness_v1" / "fp_matched_tp_robustness_summary.json",
        "ref_background": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "ref_background_v2" / "ref_background_summary.json",
        "strict_summary": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "strict_survival_v1" / "strict_null_assignment_summary.json",
        "strict_sites": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "strict_survival_v1" / "strict_null_assignment_site_results.tsv",
        "topology": TOPIC_ROOT / "results" / "focal_alt_multicluster" / "latest_full_v3_frozen" / "latest_topology_context_summary.json",
        "replication": TOPIC_ROOT / "results" / "hcc1395_cross_dataset_replication.v3.json",
        "compatibility": TOPIC_ROOT / "results" / "latest_canonical_compatibility_audit.json",
        "coverage_audit": TOPIC_ROOT / "results" / "significance_summary_coverage_audit.json",
        "materialization_retry1": TOPIC_ROOT / "results" / "latest_tagged_subset_materialization.retry1.json",
        "materialization_retry2": TOPIC_ROOT / "results" / "latest_tagged_subset_materialization.retry2.json",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing report inputs: {missing}")

    preflight = load_json(paths["preflight"])
    fp_summary = load_json(paths["fp_summary"])
    tp_summary = load_json(paths["tp_summary"])
    comparison = load_json(paths["comparison"])
    robustness = load_json(paths["robustness"])
    ref_background = load_json(paths["ref_background"])
    strict_summary = load_json(paths["strict_summary"])
    topology = load_json(paths["topology"])
    replication = load_json(paths["replication"])
    compatibility = load_json(paths["compatibility"])
    coverage_audit = load_json(paths["coverage_audit"])
    fp_rows = load_tsv(paths["fp_sites"])
    strict_rows = load_tsv(paths["strict_sites"])

    if len(fp_rows) != fp_summary["pooled_site_weighted"]["n_sites"]:
        raise RuntimeError("FP site table and summary denominator differ")
    indexed = {(row["sample"], row["chrom"], int(row["pos"])): row for row in fp_rows}
    if len(indexed) != len(fp_rows):
        raise RuntimeError("Duplicate FP site key")

    evaluable_rows = [row for row in fp_rows if row["analysis_status"] == "evaluable"]
    stable_rows = [row for row in fp_rows if endpoint(row, "stable_null_multigroup")]
    residual_rows = [row for row in fp_rows if endpoint(row, "residual_unexplained_multigroup")]
    high_threshold_rows = [
        row for row in fp_rows if endpoint(row, "phase_anchored_robust_epigenetic_candidate")
    ]
    forced_rows = [row for row in fp_rows if endpoint(row, "forced_silhouette_multigroup")]
    strict_keys = {
        tuple(site.rsplit(":", 2)) for site in strict_summary["pooled"]["strict_sites"]
    }
    strict_keys = {(sample, chrom, int(pos)) for sample, chrom, pos in strict_keys}
    strict_fp_rows = [indexed[key] for key in strict_keys]

    forced_keys = {(row["sample"], row["chrom"], row["pos"]) for row in forced_rows}
    stable_keys = {(row["sample"], row["chrom"], row["pos"]) for row in stable_rows}
    forced_stable_both = len(forced_keys & stable_keys)

    e4_normal_af_ge_001 = sum((parse_float(row["normal_af"]) or 0) >= 0.01 for row in high_threshold_rows)
    e4_normal_af_ge_005 = sum((parse_float(row["normal_af"]) or 0) >= 0.05 for row in high_threshold_rows)
    e4_normal_alt_supported = 0
    for row in high_threshold_rows:
        normal_ad = parse_json_cell(row["normal_ad"], [])
        if len(normal_ad) >= 2 and int(normal_ad[1]) > 0:
            e4_normal_alt_supported += 1
    e4_tumor_af_ge_08 = sum((parse_float(row["caller_af"]) or 0) >= 0.8 for row in high_threshold_rows)

    linear_context = "PRIMARY_ALL_STORED_TREES_LINEAR"
    e4_linear = sum(row["layered_topology_context"] == linear_context for row in high_threshold_rows)
    non_e4 = [row for row in fp_rows if row not in high_threshold_rows]
    non_e4_linear = sum(row["layered_topology_context"] == linear_context for row in non_e4)
    topology_odds, topology_p = fisher_exact(
        [[e4_linear, len(high_threshold_rows) - e4_linear],
         [non_e4_linear, len(non_e4) - non_e4_linear]],
        alternative="two-sided",
    )

    stable_ari_robust = sum(parse_bool(row["modal_assignment_all_pairs_ari_ge_0_8"]) for row in stable_rows)
    high_threshold_ari_robust = sum(
        parse_bool(row["modal_assignment_all_pairs_ari_ge_0_8"]) for row in high_threshold_rows
    )

    strict_unique_readsets = len({(row["sample"], row["alt_readset_sha256"]) for row in strict_fp_rows})
    strict_layered_components = Counter(
        (row["sample"], row["component_id"]) for row in strict_fp_rows
    )
    strict_topology = Counter(row["layered_topology_context"] for row in strict_fp_rows)

    samples = sorted(fp_summary["per_sample"])
    truth_bed_by_sample = {
        row["sample"]: row.get("truth_bed") for row in preflight["samples"]
    }
    per_sample_rows: list[dict[str, Any]] = []
    for sample in samples:
        fp = fp_summary["per_sample"][sample]
        tp = tp_summary["per_sample"][sample]
        matched = comparison["per_sample"][sample]["all_pairs"]
        per_sample_rows.append(
            {
                "sample": sample,
                "biological_sample": "HCC1395" if sample == "HCC1395_DORADO" else sample,
                "truth_confident_region_bed_available": truth_bed_by_sample[sample] is not None,
                "fp_sites": fp["n_sites"],
                "fp_evaluable": fp["n_evaluable_alt_ge6_after_peel"],
                "fp_evaluable_fraction": fp["evaluable_fraction_all"],
                "fp_stable_n": fp["n_stable_null_multigroup"],
                "fp_stable_fraction_evaluable": fp["stable_null_fraction_evaluable"],
                "fp_residual_n": fp["n_residual_unexplained_multigroup"],
                "fp_residual_fraction_evaluable": fp["residual_fraction_evaluable"],
                "fp_hp_tag_covered_high_threshold_n": fp["n_phase_anchored_robust_epigenetic_candidate"],
                "fp_hp_tag_covered_high_threshold_fraction_evaluable": fp["phase_anchored_robust_fraction_evaluable"],
                "strict_followup_n": strict_summary["per_sample"][sample]["n_strict_epigenetic_followup_candidates"],
                "tp_control_evaluable": tp["n_evaluable_alt_ge6_after_peel"],
                "matched_stable_risk_difference": matched["stable_null_multigroup"]["paired_risk_difference_fp_minus_tp"],
                "matched_stable_mcnemar_p": matched["stable_null_multigroup"]["exact_mcnemar_p"],
                "matched_high_threshold_risk_difference": matched["phase_anchored_robust_epigenetic_candidate"]["paired_risk_difference_fp_minus_tp"],
                "matched_high_threshold_mcnemar_p": matched["phase_anchored_robust_epigenetic_candidate"]["exact_mcnemar_p"],
            }
        )

    case_rows: list[dict[str, Any]] = []
    for spec in selected_case_specs():
        key = (spec["sample"], spec["chrom"], spec["pos"])
        if key not in indexed:
            raise RuntimeError(f"Missing selected case: {key}")
        row = indexed[key]
        case_rows.append(
            {
                **spec,
                "ref": row["ref"],
                "alt": row["alt"],
                "n_alt_raw": int(row["n_alt_raw"]),
                "n_alt_after_peel": int(row["n_alt_after_peel"]),
                "cluster_sizes": parse_json_cell(row["cluster_sizes"], {}),
                "alt_hp_counts": parse_json_cell(row["alt_hp_counts"], {}),
                "caller_af": parse_float(row["caller_af"]),
                "normal_af": parse_float(row["normal_af"]),
                "stable": endpoint(row, "stable_null_multigroup"),
                "residual": endpoint(row, "residual_unexplained_multigroup"),
                "hp_tag_covered_high_threshold": endpoint(
                    row, "phase_anchored_robust_epigenetic_candidate"
                ),
                "hp_axis_confound": endpoint(row, "hp_axis_confound"),
                "technical_axis_confound": endpoint(row, "technical_axis_confound"),
                "strict_followup": key in strict_keys,
                "topology_context": row["layered_topology_context"],
                "region_dir": row["region_dir"],
            }
        )

    receipt_by_sample: dict[str, dict[str, Any]] = {}
    for receipt_path in (paths["materialization_retry1"], paths["materialization_retry2"]):
        for receipt in load_json(receipt_path)["receipts"]:
            receipt_by_sample[receipt["sample"]] = receipt
    if sorted(receipt_by_sample) != samples:
        raise RuntimeError("Final tagged materialization receipts do not cover all samples")
    materialization_totals = {
        "n_samples": len(receipt_by_sample),
        "written_unique_alignments": sum(
            row["diagnostics"]["written_unique_alignments"] for row in receipt_by_sample.values()
        ),
        "with_MM_ML": sum(row["diagnostics"]["with_MM_ML"] for row in receipt_by_sample.values()),
        "duplicate_identity_collapsed": sum(
            row["diagnostics"]["duplicate_identity_collapsed"] for row in receipt_by_sample.values()
        ),
        "output_size_bytes": sum(row["output_size_bytes"] for row in receipt_by_sample.values()),
        "all_pass": all(row["pass"] for row in receipt_by_sample.values()),
        "all_sidecar_missing_zero": all(
            row["diagnostics"]["sidecar_missing"] == 0 for row in receipt_by_sample.values()
        ),
        "all_outputs_outside_canonical": all(
            "/output/canonical/" not in row["output_bam"] for row in receipt_by_sample.values()
        ),
        "output_bams": {sample: row["output_bam"] for sample, row in sorted(receipt_by_sample.items())},
    }

    all_pairs = comparison["pooled_site_weighted"]["all_pairs"]
    dual_eval = comparison["pooled_site_weighted"]["both_evaluable"]
    biological_robustness = robustness["results"]["all_pairs"]
    biological_robustness_dual = robustness["results"]["both_evaluable"]

    report = {
        "schema_name": "intersubmod.single_fp_alt_multicluster_subclone_validation_report_dataset",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "research_goals": ["G3", "G4", "G5"],
        "verdict": {
            "subclone_claim": "NOT_SUPPORTED",
            "linear_evolution_claim": "NOT_IDENTIFIABLE_FROM_ONE_SSNV",
            "supported_claim": "A small subset shows reproducible focal-ALT read-level methylation heterogeneity worthy of orthogonal follow-up.",
            "recommended_term": "strict epigenetic follow-up candidate",
        },
        "scope_and_provenance": {
            "datasets": preflight["dataset_count"],
            "biological_samples": preflight["biological_sample_count"],
            "scope": preflight["scope"],
            "latest_manifest": preflight["latest_manifest"],
            "latest_fp": preflight["totals"]["latest_fp"],
            "latest_tp": preflight["totals"]["latest_tp"],
            "latest_fp_overlap_legacy": preflight["totals"]["latest_fp_overlap_canonical"],
            "latest_fp_new_or_previously_unmaterialized": preflight["totals"]["latest_fp_promoted_or_previously_unmaterialized"],
            "colo829_truth_bed_available": truth_bed_by_sample["COLO829"] is not None,
            "site_result_sha256": sha256(paths["fp_sites"]),
            "analysis_script_sha256": sha256(TOPIC_ROOT / "scripts" / "analyze_focal_alt_multicluster.py"),
            "cluster_library_sha256": sha256(TOPIC_ROOT / "scripts" / "focal_alt_cluster_lib.py"),
            "git_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=TOPIC_ROOT, text=True
            ).strip(),
            "working_tree_dirty": bool(
                subprocess.check_output(["git", "status", "--porcelain"], cwd=TOPIC_ROOT, text=True).strip()
            ),
        },
        "tagged_bam_materialization": materialization_totals,
        "fp_primary": {
            "n_sites": len(fp_rows),
            "n_evaluable": len(evaluable_rows),
            "evaluable_fraction_all": rate(len(evaluable_rows), len(fp_rows)),
            "n_forced_silhouette_multigroup": len(forced_rows),
            "forced_fraction_evaluable": rate(len(forced_rows), len(evaluable_rows)),
            "n_stable": len(stable_rows),
            "stable_fraction_evaluable": rate(len(stable_rows), len(evaluable_rows)),
            "stable_fraction_all": rate(len(stable_rows), len(fp_rows)),
            "n_residual": len(residual_rows),
            "residual_fraction_evaluable": rate(len(residual_rows), len(evaluable_rows)),
            "residual_fraction_all": rate(len(residual_rows), len(fp_rows)),
            "n_hp_tag_covered_high_threshold": len(high_threshold_rows),
            "hp_tag_covered_high_threshold_fraction_evaluable": rate(len(high_threshold_rows), len(evaluable_rows)),
            "hp_tag_covered_high_threshold_fraction_all": rate(len(high_threshold_rows), len(fp_rows)),
            "n_orthogonal_subclone_confirmed": sum(endpoint(row, "orthogonal_subclone_confirmed") for row in fp_rows),
            "forced_stable_overlap": forced_stable_both,
            "forced_stable_jaccard": rate(forced_stable_both, len(forced_keys | stable_keys)),
            "forced_only": len(forced_keys - stable_keys),
            "stable_only": len(stable_keys - forced_keys),
            "stable_assignment_ari_ge_0_8": stable_ari_robust,
            "hp_tag_covered_high_threshold_assignment_ari_ge_0_8": high_threshold_ari_robust,
            "coverage_bands": coverage_bands(evaluable_rows),
            "shared_read_dependence": shared_read_components(paths["fp_assignments"]),
        },
        "hp_tags": {
            "all_fp": aggregate_hp(fp_rows),
            "stable": aggregate_hp(stable_rows),
            "hp_tag_covered_high_threshold": aggregate_hp(high_threshold_rows),
            "strict_followup": aggregate_hp(strict_fp_rows),
            "guardrail": "HP 1-1/2-1 are LongPhase-S first-somatic-haplotype tags; methylation cluster labels with similar text are separate objects.",
        },
        "normal_and_ref_background": {
            "e4_n": len(high_threshold_rows),
            "e4_normal_alt_supported_n": e4_normal_alt_supported,
            "e4_normal_af_ge_0_01_n": e4_normal_af_ge_001,
            "e4_normal_af_ge_0_05_n": e4_normal_af_ge_005,
            "e4_tumor_af_ge_0_8_n": e4_tumor_af_ge_08,
            "tumor_ref_background": ref_background,
            "normal_methylation_background": "NOT_TESTED",
        },
        "strict_sensitivity": {
            **strict_summary["pooled"],
            "n_unique_alt_readsets": strict_unique_readsets,
            "n_unique_layered_regional_components": len(strict_layered_components),
            "layered_regional_component_site_counts": {
                f"{sample}:{component_id}": count
                for (sample, component_id), count in sorted(strict_layered_components.items())
            },
            "topology_context_counts": dict(strict_topology),
            "guardrail": strict_summary["guardrail"],
        },
        "matched_tp_specificity": {
            "n_pairs": comparison["n_pairs"],
            "all_pairs": {
                endpoint_name: all_pairs[endpoint_name]
                for endpoint_name in (
                    "evaluable",
                    "stable_null_multigroup",
                    "residual_unexplained_multigroup",
                    "phase_anchored_robust_epigenetic_candidate",
                )
            },
            "both_evaluable": {
                endpoint_name: dual_eval[endpoint_name]
                for endpoint_name in (
                    "stable_null_multigroup",
                    "residual_unexplained_multigroup",
                    "phase_anchored_robust_epigenetic_candidate",
                )
            },
            "stable_biological_sample_equal_weight": biological_robustness["stable_null_multigroup"]["biological_sample_equal_weight"],
            "high_threshold_biological_sample_equal_weight": biological_robustness["phase_anchored_robust_epigenetic_candidate"]["biological_sample_equal_weight"],
            "both_evaluable_stable_biological_sample_equal_weight": biological_robustness_dual["stable_null_multigroup"]["biological_sample_equal_weight"],
            "both_evaluable_leave_hcc1395_stable": biological_robustness_dual["stable_null_multigroup"]["leave_one_biological_sample_out_pooled"]["HCC1395"],
            "leave_hcc1395_biological_sample_out": {
                "stable": biological_robustness["stable_null_multigroup"]["leave_one_biological_sample_out_pooled"]["HCC1395"],
                "residual": biological_robustness["residual_unexplained_multigroup"]["leave_one_biological_sample_out_pooled"]["HCC1395"],
                "high_threshold": biological_robustness["phase_anchored_robust_epigenetic_candidate"]["leave_one_biological_sample_out_pooled"]["HCC1395"],
            },
        },
        "topology": {
            **topology["counts"],
            "e4_linear_context_fisher_odds_ratio": topology_odds,
            "e4_linear_context_fisher_p": topology_p,
            "guardrail": topology["guardrail"],
        },
        "technical_replication": replication,
        "legacy_sensitivity": compatibility,
        "summary_coverage_audit": coverage_audit,
        "per_sample": per_sample_rows,
        "case_catalog": case_rows,
        "limitations": [
            "No matched-normal methylation background was run in the latest InterSubMod batch.",
            "No copy-number/purity-adjusted cancer-cell fraction or orthogonal single-cell lineage evidence is available.",
            "HCC1395 and HCC1395_DORADO are technical datasets from one biological sample.",
            "COLO829 truth splitting lacks a confident-region BED in the current preflight.",
            "Site-level confidence intervals do not fully account for shared reads and repeated readsets.",
            "Strict survival is a sensitivity gate, not a calibrated genome-wide FDR procedure.",
        ],
        "pass": True,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / "report_dataset.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    per_sample_path = args.output_dir / "per_sample_metrics.tsv"
    with per_sample_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_sample_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(per_sample_rows)

    strict_path = args.output_dir / "strict_followup_candidates.tsv"
    strict_fields = [
        "sample", "chrom", "pos", "ref", "alt", "n_alt_raw", "cluster_sizes", "alt_readset_sha256",
        "caller_af", "normal_af", "longphase_filter_transition", "component_id", "component_size",
        "ssnv_branch", "layered_topology_context", "region_dir",
    ]
    with strict_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=strict_fields, delimiter="\t")
        writer.writeheader()
        for row in sorted(strict_fp_rows, key=lambda value: (value["sample"], value["chrom"], int(value["pos"]))):
            writer.writerow({field: row[field] for field in strict_fields})

    case_path = args.output_dir / "case_catalog.tsv"
    case_fields = list(case_rows[0])
    with case_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=case_fields, delimiter="\t")
        writer.writeheader()
        for row in case_rows:
            serialized = {
                key: json.dumps(value, ensure_ascii=False, separators=(",", ":"))
                if isinstance(value, (dict, list))
                else value
                for key, value in row.items()
            }
            writer.writerow(serialized)

    print(
        json.dumps(
            {
                "report_dataset": str(report_path),
                "per_sample_table": str(per_sample_path),
                "strict_candidate_table": str(strict_path),
                "case_catalog": str(case_path),
                "n_fp_sites": len(fp_rows),
                "n_strict_followup": len(strict_fp_rows),
                "pass": report["pass"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
