#!/usr/bin/env python3
"""Build the canonical Data Analytics artifact for the VAF/CCF validation report.

This script intentionally fails closed when any required PyClone-VI analysis
table is absent, incomplete, or contains a non-PASS fit.  It only authors the
canonical ``artifact.json`` input; the shared portable report builder remains
the sole HTML renderer.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence


SCRIPT_PATH = Path(__file__).resolve()
TOPIC_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]

RAW_REQUIRED = {
    "summary": "data/raw_vaf_summary.tsv",
    "histogram": "data/raw_vaf_histogram_0.02.tsv",
    "bands": "data/raw_vaf_band_composition.tsv",
    "hcc_metrics": "data/hcc1395_pair_metrics.tsv",
    "hcc_records": "data/hcc1395_pair_exact_shared_records.tsv.gz",
    "truth_stats": "data/truth_intersection_stats.tsv",
    "neighbor": "pairwise_distributions/data/hcc1395_distribution_neighbor_summary.tsv",
    "pairwise": "pairwise_distributions/data/pairwise_truth_confirmed_distribution_distances.tsv",
    "provenance": "provenance_and_claim_ceiling.json",
}

PYCLONE_REQUIRED = {
    "analysis_summary": "analysis_summary.json",
    "fit_summaries": "fit_summaries.tsv",
    "cluster_profiles": "cluster_profiles.tsv",
    "joint": "hcc1395_joint_comparisons.tsv",
    "separate": "hcc1395_separate_fit_comparisons.tsv",
    "sensitivity": "sensitivity_comparisons.tsv",
    "all_ready": "all_ready_summary.tsv",
    "blocked": "blocked_samples.tsv",
    "provenance": "provenance.json",
}

REGIONAL_PRECISION_REQUIRED = {
    "summary": "regional_precision_summary.json",
    "metrics": "regional_precision_metrics.tsv",
    "k_counts": "k_region_counts.tsv",
    "checks": "validation_checks.tsv",
    "confusion": "coarse_confusion_by_k.tsv",
}

CLONE_BRIDGE_REQUIRED = {
    "summary": "summary.json",
    "checks": "checks.tsv",
    "join_coverage": "join_coverage.tsv",
    "concordance": "global_vs_regional_cluster_concordance.tsv",
    "patterns": "region_cluster_pattern_summary.tsv",
    "edges": "directed_edge_summary.tsv",
    "flags": "diagnostic_flags.tsv",
}

CAUSE_AUDIT_REQUIRED = {
    "summary": "summary.json",
    "checks": "checks.tsv",
    "strata": "factor_strata_summary.tsv",
    "contrasts": "contrast_effects.tsv",
    "models": "multivariable_logistic_models.tsv",
}

DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

PRIMARY_JOINT_BUNDLE = "hcc1395_pair_primary_joint_main"
PRIMARY_PROFILE_BUNDLES = {
    "HCC1395": PRIMARY_JOINT_BUNDLE,
    "HCC1395_DORADO": PRIMARY_JOINT_BUNDLE,
    "H1437": "H1437_individual_main",
    "H2009": "H2009_individual_main",
    "HCC1954": "HCC1954_individual_main",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=TOPIC_ROOT / "results" / "raw_vaf_validation_v1",
    )
    parser.add_argument(
        "--pyclone-analysis-dir",
        type=Path,
        default=TOPIC_ROOT / "runs" / "pyclone_vi" / "analysis",
    )
    parser.add_argument(
        "--previous-artifact",
        type=Path,
        default=REPO_ROOT
        / "docs"
        / "reports"
        / "in_progress"
        / "2026"
        / "07"
        / "20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01"
        / "artifact.json",
    )
    parser.add_argument(
        "--topology-context",
        type=Path,
        default=TOPIC_ROOT / "results" / "integrated_topology_context_v1.json",
    )
    parser.add_argument(
        "--regional-precision-dir",
        type=Path,
        default=TOPIC_ROOT / "results" / "regional_precision_validation_v1",
    )
    parser.add_argument(
        "--clone-region-bridge-dir",
        type=Path,
        default=TOPIC_ROOT / "results" / "clone_region_bridge_v1",
    )
    parser.add_argument(
        "--cause-audit-dir",
        type=Path,
        default=TOPIC_ROOT / "results" / "cause_decomposition_audit_v1",
    )
    parser.add_argument("--output", type=Path, default=TOPIC_ROOT / "report" / "artifact.json")
    parser.add_argument(
        "--scatter-max-rows",
        type=int,
        default=500,
        help="Deterministic maximum rows for each exact-locus scatter snapshot",
    )
    return parser.parse_args()


def require_files(root: Path, files: Mapping[str, str], label: str) -> dict[str, Path]:
    resolved = {key: root / relative for key, relative in files.items()}
    missing = [str(path) for path in resolved.values() if not path.is_file()]
    if missing:
        joined = "\n  - ".join(missing)
        raise FileNotFoundError(f"{label} is incomplete; required files are missing:\n  - {joined}")
    return resolved


INTEGER_RE = re.compile(r"^[+-]?\d+$")
FLOAT_RE = re.compile(r"^[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+-]?\d+)?$")


def typed_value(value: str) -> Any:
    value = value.strip()
    if value == "":
        return None
    if INTEGER_RE.fullmatch(value):
        return int(value)
    if FLOAT_RE.fullmatch(value):
        converted = float(value)
        if not math.isfinite(converted):
            raise ValueError(f"Non-finite numeric value is not portable JSON: {value}")
        return converted
    return value


def read_tsv(path: Path, required_columns: Iterable[str]) -> list[dict[str, Any]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = set(required_columns).difference(columns)
        if missing:
            raise ValueError(f"{path} lacks required columns: {sorted(missing)}")
        rows = [{key: typed_value(value) for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path} has no data rows")
    return rows


def read_json(path: Path) -> Any:
    with path.open() as handle:
        value = json.load(handle)
    if value is None:
        raise ValueError(f"{path} is empty JSON")
    return value


def pct(value: float, digits: int = 2) -> str:
    return f"{100 * value:.{digits}f}%"


def num(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def first(rows: Sequence[Mapping[str, Any]], **conditions: Any) -> Mapping[str, Any]:
    for row in rows:
        if all(row.get(key) == value for key, value in conditions.items()):
            return row
    raise ValueError(f"No row matches {conditions}")


def repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError as exc:
        raise ValueError(f"Artifact provenance must remain repo-relative: {path}") from exc


def bounded_exact_loci(path: Path, limit: int) -> list[dict[str, Any]]:
    if limit < 100:
        raise ValueError("--scatter-max-rows must be >=100")
    required = {
        "key",
        "af_hcc1395",
        "af_dorado",
        "dp_hcc1395",
        "dp_dorado",
        "truth_confirmed",
        "source",
        "vaf_delta_dorado_minus_hcc1395",
        "vaf_mean",
    }
    rows = read_tsv(path, required)
    eligible = [row for row in rows if row["source"] == "latest_lps_pass"]
    if not eligible:
        raise ValueError("No latest_lps_pass exact-locus rows")
    eligible.sort(key=lambda row: hashlib.sha256(str(row["key"]).encode()).hexdigest())
    sampled = eligible[:limit]
    return [
        {
            "locus": row["key"],
            "HCC1395_VAF": row["af_hcc1395"],
            "HCC1395_DORADO_VAF": row["af_dorado"],
            "mean_VAF": row["vaf_mean"],
            "DORADO_minus_HCC": row["vaf_delta_dorado_minus_hcc1395"],
            "HCC1395_DP": row["dp_hcc1395"],
            "HCC1395_DORADO_DP": row["dp_dorado"],
            "truth_scope": "SEQC2-confirmed" if str(row["truth_confirmed"]).lower() == "true" else "caller-only",
        }
        for row in sampled
    ]


def normalize_raw(raw_paths: Mapping[str, Path], scatter_limit: int) -> dict[str, Any]:
    summary = read_tsv(
        raw_paths["summary"],
        {"sample", "source", "scope", "n", "vaf_mean", "vaf_median", "dp_median"},
    )
    histogram = read_tsv(
        raw_paths["histogram"],
        {"sample", "source", "scope", "bin_left", "bin_right", "count", "fraction"},
    )
    bands = read_tsv(
        raw_paths["bands"],
        {"sample", "source", "scope", "band", "count", "fraction"},
    )
    metrics = read_tsv(
        raw_paths["hcc_metrics"],
        {
            "source",
            "scope",
            "metric",
            "value",
            "ci95_low_1mb_block_bootstrap",
            "ci95_high_1mb_block_bootstrap",
            "n_union",
            "n_shared",
        },
    )
    truth_stats = read_tsv(
        raw_paths["truth_stats"],
        {"sample", "source", "truth_label", "exact_truth_intersection_n"},
    )
    neighbor = read_tsv(
        raw_paths["neighbor"],
        {
            "source",
            "metric",
            "hcc1395_nearest_neighbor",
            "dorado_rank_among_6_peers",
            "dorado_distance",
            "nearest_other_cell_line",
            "nearest_other_over_dorado_distance_ratio",
            "median_other_over_dorado_distance_ratio",
        },
    )
    pairwise = read_tsv(
        raw_paths["pairwise"],
        {"source", "sample_a", "sample_b", "metric", "distance", "n_a", "n_b"},
    )

    latest_summary = [
        row for row in summary if row["source"] == "latest_lps_pass" and row["scope"] == "truth_confirmed"
    ]
    latest_hist = [
        {
            **row,
            "bin_mid": (float(row["bin_left"]) + float(row["bin_right"])) / 2,
        }
        for row in histogram
        if row["source"] == "latest_lps_pass" and row["scope"] == "truth_confirmed"
    ]
    latest_bands = [
        row for row in bands if row["source"] == "latest_lps_pass" and row["scope"] == "truth_confirmed"
    ]
    if {row["sample"] for row in latest_summary} != set(DATASETS):
        raise ValueError("Latest truth-confirmed VAF summary does not cover all seven dataset rows")
    if {row["sample"] for row in latest_hist} != set(DATASETS):
        raise ValueError("Latest truth-confirmed VAF histogram does not cover all seven dataset rows")
    if {row["sample"] for row in latest_bands} != set(DATASETS):
        raise ValueError("Latest truth-confirmed VAF bands do not cover all seven dataset rows")

    exact_rows = [
        row for row in metrics if row["source"] == "latest_lps_pass" and row["scope"] == "exact_caller_union"
    ]
    metric_map = {row["metric"]: row for row in exact_rows}
    required_metrics = {
        "callset_jaccard",
        "vaf_ccc",
        "vaf_spearman",
        "vaf_mae",
        "vaf_median_abs_delta",
        "vaf_within_0.10",
        "vaf_mean_signed_delta_dorado_minus_hcc1395",
        "binomial_noise_rmse_ratio",
        "binomial_95pct_compatible_fraction",
    }
    missing_metrics = required_metrics.difference(metric_map)
    if missing_metrics:
        raise ValueError(f"Latest exact-locus metrics are incomplete: {sorted(missing_metrics)}")

    latest_neighbors = [row for row in neighbor if row["source"] == "latest_lps_pass"]
    if {row["metric"] for row in latest_neighbors} != {
        "js_divergence_50bin_nats",
        "ks_statistic",
        "wasserstein_vaf",
    }:
        raise ValueError("Latest nearest-neighbor table must contain JSD, KS, and Wasserstein metrics")
    for row in latest_neighbors:
        if row["hcc1395_nearest_neighbor"] != "HCC1395_DORADO" or row["dorado_rank_among_6_peers"] != 1:
            raise ValueError("HCC pair is not the latest-spectrum nearest-neighbor pair; report narrative must be revisited")

    return {
        "summary": latest_summary,
        "histogram": latest_hist,
        "bands": latest_bands,
        "metrics": exact_rows,
        "metric_map": metric_map,
        "truth_stats": truth_stats,
        "neighbor": latest_neighbors,
        "pairwise": [row for row in pairwise if row["source"] == "latest_lps_pass"],
        "scatter": bounded_exact_loci(raw_paths["hcc_records"], scatter_limit),
        "provenance": read_json(raw_paths["provenance"]),
    }


def normalize_pyclone(paths: Mapping[str, Path]) -> dict[str, Any]:
    analysis_summary = read_json(paths["analysis_summary"])
    provenance = read_json(paths["provenance"])
    fit_summaries = read_tsv(
        paths["fit_summaries"],
        {
            "bundle_id",
            "samples",
            "mutations",
            "clusters",
            "assignment_prob_mean",
            "assignment_prob_median",
            "assignment_prob_lt_0.8_fraction",
        },
    )
    profiles = read_tsv(
        paths["cluster_profiles"],
        {
            "bundle_id",
            "sample_id",
            "cluster_id",
            "mutation_count",
            "mutation_fraction",
            "cellular_prevalence",
            "cellular_prevalence_std",
            "clone_class",
            "assignment_probability_mean",
            "assignment_probability_median",
            "assignment_probability_lt_0_8_fraction",
        },
    )
    joint = read_tsv(
        paths["joint"],
        {
            "bundle_id",
            "mutations",
            "clusters",
            "weighted_mean_abs_cp_delta",
            "cp_pearson",
            "cp_spearman",
            "mutation_weighted_cluster_cp_profile_jsd_bits",
            "clonal_state_agreement",
            "clonal_jaccard",
        },
    )
    separate = read_tsv(
        paths["separate"],
        {
            "scenario",
            "mutations_intersection",
            "clusters_a",
            "clusters_b",
            "ari",
            "nmi",
            "hungarian_agreement",
            "cluster_mutation_count_profile_jsd_bits",
            "cp_mae",
            "cp_pearson",
            "cp_spearman",
            "clonal_state_agreement",
            "clonal_jaccard",
        },
    )
    sensitivity = read_tsv(
        paths["sensitivity"],
        {
            "comparison",
            "sample_id",
            "mutations_intersection",
            "mutations_main",
            "mutations_sensitivity",
            "clusters_main",
            "clusters_sensitivity",
            "ari",
            "nmi",
            "hungarian_agreement",
            "prevalence_mae",
            "prevalence_spearman",
        },
    )
    all_ready = read_tsv(
        paths["all_ready"],
        {"bundle_id", "status", "mutations", "samples", "clusters", "num_restarts", "seed", "results_path"},
    )
    for row in all_ready:
        if row.get("results_path"):
            row["results_path"] = repo_relative(Path(str(row["results_path"])))
    blocked = read_tsv(paths["blocked"], {"sample", "status", "reason"})

    non_pass = [row for row in all_ready if row["status"] != "PASS"]
    if non_pass:
        states = ", ".join(f"{row['bundle_id']}={row['status']}" for row in non_pass)
        raise RuntimeError(f"PyClone analysis is not final; refusing to build report: {states}")
    if len(all_ready) != 14:
        raise RuntimeError(f"Expected 14/14 PyClone bundles, observed {len(all_ready)}/14")
    expected_bundles = {row["bundle_id"] for row in all_ready}
    if {row["bundle_id"] for row in fit_summaries} != expected_bundles:
        raise RuntimeError("fit_summaries.tsv and all_ready_summary.tsv bundle sets differ")
    primary_joint = first(joint, bundle_id=PRIMARY_JOINT_BUNDLE)
    primary_profiles = [row for row in profiles if row["bundle_id"] == PRIMARY_JOINT_BUNDLE]
    if {row["sample_id"] for row in primary_profiles} != {"HCC1395", "HCC1395_DORADO"}:
        raise RuntimeError("Primary joint cluster profile lacks one HCC technical source")
    primary_separate_candidates = [row for row in separate if str(row["scenario"]).lower() == "main"]
    if len(primary_separate_candidates) != 1:
        raise RuntimeError("Expected exactly one primary/main separate-fit comparison")
    primary_separate = primary_separate_candidates[0]
    blocked_map = {row["sample"]: row for row in blocked}
    if not {"COLO829", "HCC1937"}.issubset(blocked_map):
        raise RuntimeError("blocked_samples.tsv must explicitly contain COLO829 and HCC1937")
    if any(blocked_map[sample]["status"] != "BLOCKED" for sample in ("COLO829", "HCC1937")):
        raise RuntimeError("COLO829 and HCC1937 must remain fail-closed for CN-aware PyClone")

    mode_map = analysis_summary.get("hcc1395_model_mode_comparisons")
    if not isinstance(mode_map, dict) or set(mode_map) != {"HCC1395", "HCC1395_DORADO"}:
        raise RuntimeError("analysis_summary.json lacks the two HCC joint-vs-separate model-mode comparisons")
    model_mode = [mode_map[sample] for sample in ("HCC1395", "HCC1395_DORADO")]
    separate_map = analysis_summary.get("hcc1395_separate_fit_comparisons")
    if not isinstance(separate_map, dict) or any(
        not isinstance(separate_map.get(scenario), dict) for scenario in ("main", "near_integer")
    ):
        raise RuntimeError("analysis_summary.json lacks main/near-integer independent minor-group metrics")
    required_minor = {
        "subclonal_mutation_jaccard",
        "clonal_subclonal_state_kappa",
        "clonal_subclonal_state_mcc",
        "subclonal_state_f1",
    }
    independent_minor = []
    for scenario in ("main", "near_integer"):
        separate_detail = separate_map[scenario]
        both_subclonal = separate_detail.get("both_subclonal_cluster_metrics")
        if not required_minor.issubset(separate_detail) or not isinstance(both_subclonal, dict):
            raise RuntimeError(f"Independent-fit minor-group metrics are incomplete: {scenario}")
        if "ari" not in both_subclonal or "hungarian_agreement" not in both_subclonal:
            raise RuntimeError(f"Both-subclonal partition metrics are incomplete: {scenario}")
        independent_minor.append(
            {
                "scenario": f"independent {scenario} fits; minor-group focus",
                "subclonal_mutation_jaccard": separate_detail["subclonal_mutation_jaccard"],
                "clonal_subclonal_state_kappa": separate_detail["clonal_subclonal_state_kappa"],
                "clonal_subclonal_state_mcc": separate_detail["clonal_subclonal_state_mcc"],
                "subclonal_state_f1": separate_detail["subclonal_state_f1"],
                "both_subclonal_ari": both_subclonal["ari"],
                "both_subclonal_nmi": both_subclonal.get("nmi"),
                "both_subclonal_hungarian_agreement": both_subclonal["hungarian_agreement"],
                "both_subclonal_mutations": separate_detail.get("subclonal_contingency", {}).get("both_subclonal"),
                "either_subclonal_mutations": separate_detail.get("subclonal_contingency", {}).get("either_subclonal"),
                "interpretation": separate_detail.get("imbalance_note", ""),
            }
        )
    sensitivity_map = analysis_summary.get("sensitivity_comparisons")
    if not isinstance(sensitivity_map, dict):
        raise RuntimeError("analysis_summary.json lacks detailed sensitivity comparisons")
    sensitivity_detail = []
    for comparison, comparison_payload in sensitivity_map.items():
        samples_payload = comparison_payload.get("samples") if isinstance(comparison_payload, dict) else None
        if not isinstance(samples_payload, dict):
            raise RuntimeError(f"Sensitivity comparison has no sample mapping: {comparison}")
        for sample_id, values in samples_payload.items():
            row = {"comparison": comparison, "sample_id": sample_id, **values}
            row["sensitivity_retention_fraction"] = row["mutations_sensitivity"] / row["mutations_main"]
            sensitivity_detail.append(row)
    if {(row["comparison"], row["sample_id"]) for row in sensitivity_detail} != {
        (row["comparison"], row["sample_id"]) for row in sensitivity
    }:
        raise RuntimeError("Detailed JSON and sensitivity TSV comparison sets differ")

    fit_by_bundle = {row["bundle_id"]: row for row in fit_summaries}
    all_sample_profiles: list[dict[str, Any]] = []
    for sample in DATASETS:
        if sample in blocked_map:
            all_sample_profiles.append(
                {
                    "dataset": sample,
                    "status": "BLOCKED",
                    "active_clusters": None,
                    "high_CP_mutation_fraction": None,
                    "normalized_cluster_entropy": None,
                    "assignment_probability_mean": None,
                    "assignment_probability_median": None,
                    "assignment_probability_lt_0_8_fraction": None,
                    "reason": blocked_map[sample]["reason"],
                }
            )
            continue
        bundle = PRIMARY_PROFILE_BUNDLES[sample]
        sample_profiles = [
            row for row in profiles if row["bundle_id"] == bundle and row["sample_id"] == sample
        ]
        if not sample_profiles:
            raise RuntimeError(f"No primary PyClone profile for {sample} ({bundle})")
        fractions = [float(row["mutation_fraction"]) for row in sample_profiles]
        if abs(sum(fractions) - 1.0) > 0.005:
            raise RuntimeError(f"Cluster mutation fractions do not conserve to 1 for {sample}: {sum(fractions)}")
        entropy = -sum(value * math.log(value) for value in fractions if value > 0)
        normalized_entropy = entropy / math.log(len(fractions)) if len(fractions) > 1 else 0.0
        high_cp_fraction = sum(
            float(row["mutation_fraction"])
            for row in sample_profiles
            if float(row["cellular_prevalence"]) >= 0.9
        )
        fit = fit_by_bundle[bundle]
        low_assignment_fraction = float(fit["assignment_prob_lt_0.8_fraction"])
        confidence_note = (
            "QA PASS, but cluster assignments are low-confidence; do not treat cluster count as stable truth"
            if low_assignment_fraction >= 0.5
            else "QA PASS; still model-conditional"
        )
        all_sample_profiles.append(
            {
                "dataset": sample,
                "status": "PASS / LOW CLUSTER CONFIDENCE" if low_assignment_fraction >= 0.5 else "PASS",
                "active_clusters": len(sample_profiles),
                "high_CP_mutation_fraction": high_cp_fraction,
                "normalized_cluster_entropy": normalized_entropy,
                "assignment_probability_mean": fit["assignment_prob_mean"],
                "assignment_probability_median": fit["assignment_prob_median"],
                "assignment_probability_lt_0_8_fraction": low_assignment_fraction,
                "reason": confidence_note,
            }
        )

    return {
        "analysis_summary": analysis_summary,
        "fit_summaries": fit_summaries,
        "cluster_profiles": profiles,
        "primary_profiles": primary_profiles,
        "joint": joint,
        "primary_joint": primary_joint,
        "separate": separate,
        "primary_separate": primary_separate,
        "independent_minor": independent_minor,
        "model_mode": model_mode,
        "sensitivity": sensitivity_detail,
        "all_ready": all_ready,
        "blocked": blocked,
        "all_sample_profiles": all_sample_profiles,
        "provenance": provenance,
    }


def normalize_topology(path: Path) -> dict[str, Any]:
    data = read_json(path)
    if data.get("schema_name") != "intersubmod.integrated_topology_context_snapshot":
        raise ValueError(f"Unexpected topology context schema: {data.get('schema_name')}")
    shapes = data.get("final_shape_all_datasets")
    pair = data.get("hcc1395_pair")
    if not isinstance(shapes, list) or {row.get("dataset") for row in shapes} != set(DATASETS):
        raise ValueError("Integrated topology context does not contain all seven final-shape rows")
    if not isinstance(pair, dict) or pair.get("fixed_denominator") != 5720:
        raise ValueError("Integrated topology context must preserve the fixed 5,720 HCC-pair denominator")
    composition = []
    categories = [
        ("single_only", "Single"),
        ("sister_only", "Sister"),
        ("direct_only", "Direct"),
        ("sister_and_direct", "Sister+direct"),
    ]
    for row in shapes:
        for field, label in categories:
            composition.append(
                {
                    "dataset": row["dataset"],
                    "category": label,
                    "count": row[field],
                    "fraction_resolved": row[field] / row["final_single_shape_regions"],
                }
            )
    coarse = pair["coarse_topology"]
    read_layer = pair["read_full_candidate_constraints"]
    vaf_layer = pair["vaf_official_candidate_constraints"]
    unranked = pair["unranked_exact_candidate_tree_space"]
    selected = pair["vaf_selected_tree_and_shape"]
    topology_ladder = [
        {
            "endpoint": "Coarse five-state topology",
            "agreement_n": coarse["agree_n"],
            "fixed_denominator": 5720,
            "fixed_agreement": coarse["agree_n"] / 5720,
            "evaluable_n": 5720,
            "conditional_agreement": coarse["agree_n"] / 5720,
            "claim_boundary": "Class agreement above chromosome-preserving null; not accuracy",
        },
        {
            "endpoint": "Read candidate constraints: strict + true-induced",
            "agreement_n": read_layer["strict_plus_true_induced_n"],
            "fixed_denominator": 5720,
            "fixed_agreement": read_layer["strict_plus_true_induced_n"] / 5720,
            "evaluable_n": read_layer["evaluable_n"],
            "conditional_agreement": read_layer["strict_plus_true_induced_n"] / read_layer["evaluable_n"],
            "claim_boundary": "Candidate-space compatibility; ambiguous and same-read evidence",
        },
        {
            "endpoint": "VAF candidate constraints: strict + true-induced",
            "agreement_n": vaf_layer["strict_plus_true_induced_n"],
            "fixed_denominator": 5720,
            "fixed_agreement": vaf_layer["strict_plus_true_induced_n"] / 5720,
            "evaluable_n": vaf_layer["evaluable_n"],
            "conditional_agreement": vaf_layer["strict_plus_true_induced_n"] / vaf_layer["evaluable_n"],
            "claim_boundary": "Same-read VAF argmax; 1,234 conflicts; not posterior/truth",
        },
        {
            "endpoint": "Unranked exact candidate-tree set (HP-swap tolerant)",
            "agreement_n": unranked["hp_swap_tolerant_agree_n"],
            "fixed_denominator": 5720,
            "fixed_agreement": unranked["hp_swap_tolerant_agree_n"] / 5720,
            "evaluable_n": 5720,
            "conditional_agreement": unranked["hp_swap_tolerant_agree_n"] / 5720,
            "claim_boundary": "Feasible candidate-set identity; no selected true tree",
        },
        {
            "endpoint": "VAF-unique mutation-labeled exact forest (HP-swap tolerant)",
            "agreement_n": selected["unique_mutation_labeled_exact_forest"]["hp_swap_tolerant_agree_n"],
            "fixed_denominator": 5720,
            "fixed_agreement": selected["unique_mutation_labeled_exact_forest"]["hp_swap_tolerant_agree_n"] / 5720,
            "evaluable_n": selected["unique_mutation_labeled_exact_forest"]["n"],
            "conditional_agreement": selected["unique_mutation_labeled_exact_forest"]["hp_swap_tolerant_agreement_pct"] / 100,
            "claim_boundary": "Mutation-labeled heuristic exact forest; selected subset",
        },
        {
            "endpoint": "Structure-first + VAF single shape (HP-swap tolerant)",
            "agreement_n": selected["structure_first_plus_vaf_single_shape"]["hp_swap_tolerant_agree_n"],
            "fixed_denominator": 5720,
            "fixed_agreement": selected["structure_first_plus_vaf_single_shape"]["hp_swap_tolerant_agree_n"] / 5720,
            "evaluable_n": selected["structure_first_plus_vaf_single_shape"]["n"],
            "conditional_agreement": selected["structure_first_plus_vaf_single_shape"]["hp_swap_tolerant_agreement_pct"] / 100,
            "claim_boundary": "Unlabeled branching skeleton; mutation direction removed",
        },
        {
            "endpoint": "Original Topo>1 resolved to one VAF shape (HP-swap tolerant)",
            "agreement_n": selected["original_topogt1_vaf_single_shape"]["hp_swap_tolerant_agree_n"],
            "fixed_denominator": 5720,
            "fixed_agreement": selected["original_topogt1_vaf_single_shape"]["hp_swap_tolerant_agree_n"] / 5720,
            "evaluable_n": selected["original_topogt1_vaf_single_shape"]["n"],
            "conditional_agreement": selected["original_topogt1_vaf_single_shape"]["hp_swap_tolerant_agreement_pct"] / 100,
            "claim_boundary": "Conditional selected-shape reproducibility; not exact tree",
        },
    ]
    return {"raw": data, "shapes": shapes, "composition": composition, "pair": pair, "ladder": topology_ladder}


def normalize_regional_precision(root: Path) -> dict[str, Any]:
    paths = require_files(root, REGIONAL_PRECISION_REQUIRED, "Regional precision validation")
    summary = read_json(paths["summary"])
    if summary.get("schema_name") != "intersubmod.hcc1395_regional_precision_validation":
        raise ValueError(f"Unexpected regional precision schema: {summary.get('schema_name')}")
    if summary.get("status") != "PASS_WITH_CLAIM_CEILING":
        raise RuntimeError(f"Regional precision summary is not validated: {summary.get('status')}")
    scope = summary.get("scope", {})
    if scope.get("fixed_population") != 5720:
        raise RuntimeError(f"Regional precision fixed population is not 5,720: {scope.get('fixed_population')}")

    metrics = read_tsv(
        paths["metrics"],
        {
            "stratum",
            "metric_group",
            "metric_id",
            "metric_label",
            "global_fixed_denominator",
            "stratum_fixed_denominator",
            "evaluable_denominator",
            "numerator",
            "rate_global_fixed",
            "rate_stratum_fixed",
            "rate_evaluable",
            "availability",
            "claim_ceiling",
        },
    )
    metric_keys = [(row["stratum"], row["metric_id"]) for row in metrics]
    if len(metric_keys) != len(set(metric_keys)):
        raise RuntimeError("Regional precision metrics contain duplicate stratum+metric_id keys")
    metric_by_key = {(row["stratum"], row["metric_id"]): row for row in metrics}

    k_counts = read_tsv(
        paths["k_counts"],
        {"k_stratum", "k_basis", "region_count", "rate_fixed_5720", "availability"},
    )
    if {row["k_stratum"] for row in k_counts} != {"k=1", "k=2", "k=3", "k>=4"}:
        raise RuntimeError("Regional k counts must contain exactly k=1, k=2, k=3, and k>=4")
    k_map = {row["k_stratum"]: int(row["region_count"]) for row in k_counts}
    if k_map != {"k=1": 0, "k=2": 3121, "k=3": 1506, "k>=4": 1093}:
        raise RuntimeError(f"Regional k counts changed; narrative review required: {k_map}")
    if sum(k_map.values()) != 5720:
        raise RuntimeError("Regional k counts do not conserve to 5,720")
    if scope.get("k_counts") != k_map:
        raise RuntimeError("Regional summary and k-count TSV disagree")

    checks = read_tsv(paths["checks"], {"check", "expected", "observed", "status", "evidence"})
    failed = [row for row in checks if row["status"] != "PASS"]
    if failed:
        raise RuntimeError("Regional precision validation checks failed: " + ", ".join(row["check"] for row in failed))
    if len(checks) != 32 or summary.get("validation", {}).get("pass") != 32:
        raise RuntimeError("Regional precision validation receipt must contain 32/32 PASS checks")

    confusion = read_tsv(
        paths["confusion"],
        {"stratum", "category_a", "category_b", "regions", "stratum_denominator", "diagonal"},
    )
    unresolved_k4 = first(
        confusion,
        stratum="k>=4",
        category_a="topology_multiple_unresolved",
        category_b="topology_multiple_unresolved",
    )
    if unresolved_k4["regions"] != 888 or unresolved_k4["stratum_denominator"] != 1093:
        raise RuntimeError(f"k>=4 unresolved diagonal changed; narrative review required: {unresolved_k4}")

    chart_contract = [
        ("coarse_five_state_agreement", "Coarse 5-state", "rate_stratum_fixed", "fixed k-stratum denominator"),
        ("read_strict_plus_induced", "Read strict+induced", "rate_evaluable", "read-evaluable denominator"),
        ("vaf_official_strict_plus_induced", "VAF strict+induced", "rate_evaluable", "VAF-evaluable denominator"),
        ("vaf_official_conflict", "VAF conflict", "rate_evaluable", "VAF-evaluable denominator"),
    ]
    chart_rows = []
    for stratum in ("k=2", "k=3", "k>=4"):
        for metric_id, label, value_field, denominator_basis in chart_contract:
            key = (stratum, metric_id)
            if key not in metric_by_key:
                raise RuntimeError(f"Regional precision chart metric is missing: {key}")
            row = metric_by_key[key]
            value = row[value_field]
            if value is None:
                raise RuntimeError(f"Regional precision chart metric has no value: {key}/{value_field}")
            chart_rows.append(
                {
                    "k_stratum": stratum,
                    "endpoint": label,
                    "value": value,
                    "numerator": row["numerator"],
                    "global_fixed_denominator": row["global_fixed_denominator"],
                    "stratum_fixed_denominator": row["stratum_fixed_denominator"],
                    "evaluable_denominator": row["evaluable_denominator"],
                    "denominator_basis": denominator_basis,
                }
            )

    all_metrics = [row for row in metrics if row["stratum"] == "ALL"]
    headline = {(row["metric_id"]): row for row in all_metrics}
    for metric_id in (
        "site_set_equal",
        "all_shared_alleles_identical_region",
        "coarse_five_state_agreement",
        "shape_hp_swap_agreement",
        "unique_exact_forest_hp_swap_agreement",
    ):
        if metric_id not in headline:
            raise RuntimeError(f"Regional headline metric is missing: {metric_id}")
    allele_headline = summary.get("headline", {}).get("shared_allele_identity", {})
    if allele_headline.get("identity_alleles") != 15713 or allele_headline.get("shared_alleles") != 15713:
        raise RuntimeError("Regional shared-allele identity no longer equals 15,713/15,713")

    return {
        "paths": paths,
        "summary": summary,
        "metrics": metrics,
        "k_counts": k_counts,
        "checks": checks,
        "confusion": confusion,
        "k4_unresolved_diagonal": unresolved_k4,
        "chart": chart_rows,
        "headline": headline,
        "allele_headline": allele_headline,
    }


def normalize_clone_bridge(root: Path) -> dict[str, Any]:
    paths = require_files(root, CLONE_BRIDGE_REQUIRED, "Clone-region bridge")
    summary = read_json(paths["summary"])
    if summary.get("schema_name") != "intersubmod.clone_region_bridge" or summary.get("status") != "PASS":
        raise RuntimeError(f"Clone-region bridge is not a stable PASS artifact: {summary.get('schema_name')}/{summary.get('status')}")
    if summary.get("scope") != "HCC1395 vs HCC1395_DORADO; fixed 5,720 exact-coordinate complete-both historical layered-v2 regions":
        raise RuntimeError("Clone-region bridge scope changed; narrative review required")
    checks = read_tsv(paths["checks"], {"check", "pass", "observed", "expected"})
    if len(checks) != 14 or any(row["pass"] != "True" for row in checks):
        raise RuntimeError("Clone-region bridge must have 14/14 hard checks PASS")
    if len(summary.get("checks", [])) != 14 or any(not row.get("pass") for row in summary["checks"]):
        raise RuntimeError("Clone-region bridge JSON check receipt disagrees with checks.tsv")

    coverage = read_tsv(paths["join_coverage"], {"level", "metric", "n", "denominator", "share"})
    coverage_map = {row["metric"]: row for row in coverage}
    concordance = read_tsv(
        paths["concordance"],
        {
            "population", "stratum", "n", "ari", "clonal_state_kappa", "subclonal_intersection",
            "subclonal_union", "subclonal_jaccard",
        },
    )
    regional_all = first(concordance, population="fixed_5720_regional_subset", stratum="all_joined")
    regional_high = first(concordance, population="fixed_5720_regional_subset", stratum="both_assignment_ge_0.8")
    patterns = read_tsv(
        paths["patterns"],
        {
            "endpoint", "stratum", "evaluable_multilocus_regions", "partition_exact",
            "partition_exact_share_evaluable", "both_single_cluster", "both_multicluster_exact_partition",
            "both_multicluster_different_partition", "both_multicluster_evaluable", "both_multicluster_exact_share",
        },
    )
    pattern_all = first(patterns, endpoint="all_joined", stratum="ALL")
    edges = read_tsv(
        paths["edges"],
        {
            "sample", "tolerance", "confidence_stratum", "evaluable_directed_joined_edges", "compatible",
            "conflict", "uninformative_same_cluster", "compatible_share_evaluable", "conflict_share_evaluable",
        },
    )
    edge_rows = {
        (sample, confidence): first(edges, sample=sample, tolerance=0.02, confidence_stratum=confidence)
        for sample in ("HCC1395", "HCC1395_DORADO")
        for confidence in ("all", "both_assignment_ge_0.8")
    }
    flags = read_tsv(paths["flags"], {"flag", "severity", "status", "observed", "threshold", "interpretation"})
    if len(flags) != 3 or any(row["severity"] != "WARN" for row in flags):
        raise RuntimeError("Clone-region bridge must preserve all three explicit WARN diagnostics")

    required_coverage = {
        "input_expanded_alleles": (15713, 15713),
        "joined_pyclone_both": (14369, 15713),
        "both_sources_determinate_direction": (598, 8096),
        "same_direction_among_both_determinate": (596, 598),
        "opposite_direction_among_both_determinate": (2, 598),
    }
    for metric_id, expected in required_coverage.items():
        row = coverage_map.get(metric_id)
        observed = (row.get("n"), row.get("denominator")) if row else None
        if observed != expected:
            raise RuntimeError(f"Clone bridge metric changed: {metric_id} expected={expected} observed={observed}")
    if (regional_all["subclonal_intersection"], regional_all["subclonal_union"]) != (62, 259):
        raise RuntimeError("Regional subclonal intersection/union changed")
    if regional_all["subclonal_jaccard"] != 0.23938223938223938 or regional_all["clonal_state_kappa"] != 0.3793636987708486:
        raise RuntimeError("Regional concordance headline changed")
    if (regional_high["subclonal_intersection"], regional_high["subclonal_union"]) != (1, 1):
        raise RuntimeError("High-confidence degeneracy changed")
    if (
        pattern_all["partition_exact"], pattern_all["evaluable_multilocus_regions"],
        pattern_all["both_single_cluster"], pattern_all["both_multicluster_exact_partition"],
        pattern_all["both_multicluster_evaluable"],
    ) != (5028, 5189, 5007, 21, 34):
        raise RuntimeError("Bridge partition decomposition changed")
    for sample, expected in {"HCC1395": (19, 0, 1098, 1117), "HCC1395_DORADO": (12, 0, 502, 514)}.items():
        row = edge_rows[(sample, "all")]
        observed = (row["compatible"], row["conflict"], row["uninformative_same_cluster"], row["evaluable_directed_joined_edges"])
        if observed != expected:
            raise RuntimeError(f"Bridge directed-edge headline changed for {sample}: {observed}")
    for sample in ("HCC1395", "HCC1395_DORADO"):
        row = edge_rows[(sample, "both_assignment_ge_0.8")]
        if row["compatible"] != 0 or row["conflict"] != 0 or row["uninformative_same_cluster"] != row["evaluable_directed_joined_edges"]:
            raise RuntimeError(f"High-confidence edge degeneracy changed for {sample}")

    bridge_audit = [
        {
            "layer": "join coverage", "endpoint": "PyClone joined in both sources", "population": "regional exact alleles",
            "numerator": 14369, "denominator": 15713, "value": 14369 / 15713, "status": "AVAILABLE",
            "interpretation": "91.45% of regional exact alleles enter the independent-fit bridge",
        },
        {
            "layer": "regional independent fits", "endpoint": "Subclonal mutation Jaccard", "population": "14,369 joined mutations",
            "numerator": 62, "denominator": 259, "value": regional_all["subclonal_jaccard"], "status": "PARTIAL",
            "interpretation": "Minor-group overlap is low; clonal-majority agreement is not substituted",
        },
        {
            "layer": "regional independent fits", "endpoint": "Clonal/subclonal state κ", "population": "14,369 joined mutations",
            "numerator": None, "denominator": 14369, "value": regional_all["clonal_state_kappa"], "status": "PARTIAL",
            "interpretation": "Chance-corrected binary-state agreement remains moderate/low",
        },
        {
            "layer": "assignment>=0.8", "endpoint": "Subclonal mutation Jaccard", "population": "10,965 joined mutations",
            "numerator": 1, "denominator": 1, "value": 1.0, "status": "WARN_VACUOUS",
            "interpretation": "Perfect value is selection-induced because the subclonal union contains one mutation",
        },
        {
            "layer": "regional partitions", "endpoint": "Partition exact", "population": "5,189 evaluable multilocus regions",
            "numerator": 5028, "denominator": 5189, "value": 5028 / 5189, "status": "WARN_DOMINATED",
            "interpretation": "96.9% is dominated by regions that are one cluster in both fits",
        },
        {
            "layer": "regional partitions", "endpoint": "Both single-cluster", "population": "5,189 evaluable multilocus regions",
            "numerator": 5007, "denominator": 5189, "value": 5007 / 5189, "status": "DOMINANT",
            "interpretation": "This majority makes overall partition-exact agreement weak evidence for multi-clone recovery",
        },
        {
            "layer": "regional partitions", "endpoint": "Both-multicluster exact partition", "population": "34 both-multicluster regions",
            "numerator": 21, "denominator": 34, "value": 21 / 34, "status": "INFORMATIVE_SMALL_N",
            "interpretation": "61.76% is the informative partition subset; denominator is only 34",
        },
    ]
    bridge_edge_audit = [
        {
            "sample": "Cross-source", "confidence": "read singleton direction", "endpoint": "Same determinate direction",
            "numerator": 596, "denominator": 598, "rate": 596 / 598, "conflict": 2,
            "interpretation": "Read directions reproduce, but this does not test whether endpoints occupy different clones",
        },
        {
            "sample": "Cross-source", "confidence": "read singleton direction", "endpoint": "Opposite determinate direction",
            "numerator": 2, "denominator": 598, "rate": 2 / 598, "conflict": 2,
            "interpretation": "Only two directions are opposite across sources",
        },
    ]
    for sample in ("HCC1395", "HCC1395_DORADO"):
        all_row = edge_rows[(sample, "all")]
        high_row = edge_rows[(sample, "both_assignment_ge_0.8")]
        bridge_edge_audit.extend(
            [
                {
                    "sample": sample, "confidence": "all assignments", "endpoint": "Compatible distinct-cluster direction",
                    "numerator": all_row["compatible"], "denominator": all_row["evaluable_directed_joined_edges"],
                    "rate": all_row["compatible_share_evaluable"], "conflict": all_row["conflict"],
                    "interpretation": "No CP-order conflict, but very few edges link different inferred clusters",
                },
                {
                    "sample": sample, "confidence": "all assignments", "endpoint": "Uninformative same cluster",
                    "numerator": all_row["uninformative_same_cluster"], "denominator": all_row["evaluable_directed_joined_edges"],
                    "rate": all_row["uninformative_same_cluster"] / all_row["evaluable_directed_joined_edges"],
                    "conflict": all_row["conflict"],
                    "interpretation": "Same-cluster endpoints cannot validate ancestry ordering",
                },
                {
                    "sample": sample, "confidence": "both assignment>=0.8", "endpoint": "Distinct-cluster directional information",
                    "numerator": high_row["compatible"] + high_row["conflict"], "denominator": high_row["evaluable_directed_joined_edges"],
                    "rate": 0.0, "conflict": high_row["conflict"],
                    "interpretation": "All high-confidence joined directed endpoints share one cluster; ancestry test is vacuous",
                },
            ]
        )

    return {
        "paths": paths,
        "summary": summary,
        "checks": checks,
        "flags": flags,
        "audit": bridge_audit,
        "edges": bridge_edge_audit,
    }


def normalize_cause_audit(root: Path) -> dict[str, Any]:
    paths = require_files(root, CAUSE_AUDIT_REQUIRED, "Cause-decomposition audit")
    summary = read_json(paths["summary"])
    if summary.get("task_type") != "B_comprehensive_validation_red_team_audit":
        raise RuntimeError("Cause-decomposition task type changed")
    if summary.get("population", {}).get("regions") != 5720 or summary.get("population", {}).get("hp_count_mismatch_regions") != 703:
        raise RuntimeError("Cause-decomposition population changed")
    if summary.get("checks_passed") != 16 or summary.get("checks_total") != 16:
        raise RuntimeError("Cause-decomposition summary is not 16/16 PASS")
    checks = read_tsv(paths["checks"], {"check", "observed", "expected", "pass"})
    if len(checks) != 16 or any(row["pass"] != "True" for row in checks):
        raise RuntimeError("Cause-decomposition checks.tsv is not 16/16 PASS")
    strata = read_tsv(
        paths["strata"],
        {"factor", "factor_column", "level", "factor_level_n", "endpoint", "endpoint_scope", "denominator", "numerator", "rate"},
    )
    contrasts = read_tsv(
        paths["contrasts"],
        {
            "contrast", "endpoint", "endpoint_scope", "denominator_a", "numerator_a", "rate_a",
            "denominator_b", "numerator_b", "rate_b", "rate_difference_a_minus_b", "odds_ratio_haldane",
            "rate_difference_ci95_low", "rate_difference_ci95_high", "bh_q_exploratory",
        },
    )
    models = read_tsv(
        paths["models"],
        {
            "model", "model_role", "outcome", "n", "predictor", "predictor_interpretation", "odds_ratio",
            "odds_ratio_ci95_low", "odds_ratio_ci95_high", "cluster_robust_p",
            "bh_q_across_all_nonintercept_coefficients", "warning",
        },
    )

    def contrast(contrast_id: str, endpoint: str) -> Mapping[str, Any]:
        return first(contrasts, contrast=contrast_id, endpoint=endpoint)

    def model(model_id: str, predictor: str) -> Mapping[str, Any]:
        return first(models, model=model_id, predictor=predictor)

    hp_coarse = contrast("hp_count_mismatch_vs_matched", "coarse_agree_fixed")
    hp_read_eval = contrast("hp_count_mismatch_vs_matched", "read_evaluable_fixed")
    vaf_delta = contrast("highest_vs_lowest_baseline_abs_vaf_delta_quartile", "vaf_exact_or_induced_given_evaluable")
    vaf_delta_adjusted = model("vaf_strict_plus_candidate_given_evaluable", "x_abs_vaf_delta_per_0_05")
    if vaf_delta["rate_difference_a_minus_b"] != -0.0732923453292989:
        raise RuntimeError("High-vs-low |delta VAF| contrast changed")
    if vaf_delta_adjusted["odds_ratio"] != 0.9211056112519098:
        raise RuntimeError("Adjusted |delta VAF| model changed")
    if (hp_coarse["denominator_a"], hp_coarse["numerator_a"], hp_read_eval["numerator_a"]) != (703, 357, 0):
        raise RuntimeError("HP mismatch fail-closed result changed")

    audit_rows: list[dict[str, Any]] = []

    def add_contrast(factor: str, label: str, row: Mapping[str, Any], interpretation: str) -> None:
        audit_rows.append(
            {
                "factor": factor,
                "evidence": "unadjusted block-bootstrap contrast",
                "endpoint": row["endpoint"],
                "comparison": label,
                "numerator_a": row["numerator_a"], "denominator_a": row["denominator_a"], "rate_a": row["rate_a"],
                "numerator_b": row["numerator_b"], "denominator_b": row["denominator_b"], "rate_b": row["rate_b"],
                "effect": row["rate_difference_a_minus_b"], "odds_ratio": row["odds_ratio_haldane"],
                "ci95_low": row["rate_difference_ci95_low"], "ci95_high": row["rate_difference_ci95_high"],
                "adjusted_q": row["bh_q_exploratory"], "interpretation": interpretation,
            }
        )

    def add_model(factor: str, label: str, row: Mapping[str, Any], interpretation: str) -> None:
        audit_rows.append(
            {
                "factor": factor,
                "evidence": "multivariable cluster-robust logistic model",
                "endpoint": row["outcome"],
                "comparison": label,
                "numerator_a": None, "denominator_a": row["n"], "rate_a": None,
                "numerator_b": None, "denominator_b": None, "rate_b": None,
                "effect": None, "odds_ratio": row["odds_ratio"],
                "ci95_low": row["odds_ratio_ci95_low"], "ci95_high": row["odds_ratio_ci95_high"],
                "adjusted_q": row["bh_q_across_all_nonintercept_coefficients"], "interpretation": interpretation,
            }
        )

    add_contrast("|delta VAF|", "Q4 high minus Q1 low", vaf_delta,
                 "-7.33 pp before adjustment; same-read diagnostic, not validation")
    add_model("|delta VAF|", "per 0.05 after k/depth/HP/site/CN/candidate adjustment", vaf_delta_adjusted,
              "CI crosses 1 and BH q is not small; crude VAF-delta contrast is not stable after complexity adjustment")
    add_contrast("HP mapping", "HP-count mismatch minus matched", hp_coarse,
                 "Coarse agreement is 21.21 pp lower with HP mismatch")
    add_contrast("HP mapping", "HP-count mismatch minus matched", hp_read_eval,
                 "All 703 HP-mismatch regions fail closed for fine read/VAF endpoints")
    add_contrast("read candidate multiplicity", ">20 versus 1", contrast("read_candidate_topology_gt20_vs_1", "read_exact_or_induced_given_evaluable"),
                 "Conditional exact+induced falls by 73.57 pp as read candidate space expands")
    add_contrast("VAF candidate multiplicity", ">=2 versus 1", contrast("vaf_candidate_topology_ge2_vs_1", "vaf_exact_or_induced_given_evaluable"),
                 "Conditional exact+induced falls by 46.71 pp; candidate definition/selection also affects conflict")
    add_contrast("depth", "Q4 high minus Q1 low", contrast("highest_vs_lowest_baseline_depth_quartile", "read_exact_or_induced_given_evaluable"),
                 "Higher-depth strata look worse for fine agreement, but complexity/candidate selection confounds this contrast")
    add_model("depth", "per 2x depth with candidate adjustment", model("read_strict_plus_candidate_given_evaluable", "x_log2_mean_pair_min_depth"),
              "Diagnostic association persists for read strict, but does not imply low depth is beneficial")
    add_model("depth", "per 2x depth with candidate adjustment", model("vaf_strict_plus_candidate_given_evaluable", "x_log2_mean_pair_min_depth"),
              "VAF strict association is borderline after BH adjustment; do not infer a beneficial low-depth mechanism")
    add_contrast("shared CN/LOH", "any LOH versus no LOH", contrast("full_cn_any_loh_vs_no_loh", "vaf_exact_or_induced_given_evaluable"),
                 "Descriptive shared-CN context only; DORADO has no source-specific CN corroboration")
    add_model("read candidate multiplicity", "per 2x candidate topology maximum", model("read_strict_plus_candidate_given_evaluable", "x_read_candidate_log2"),
              "More read candidates independently reduce fine strict agreement")
    add_model("VAF candidate multiplicity", "per 2x candidate topology maximum", model("vaf_strict_plus_candidate_given_evaluable", "x_vaf_candidate_log2"),
              "More VAF candidates strongly reduce strict agreement after context adjustment")
    for predictor, label in (("x_k3_vs_k2", "k=3 vs k=2"), ("x_k4_vs_k2", "k=4 vs k=2"), ("x_kge5_vs_k2", "k>=5 vs k=2")):
        add_model("shared-site k", label, model("vaf_strict_plus_candidate_given_evaluable", predictor),
                  "Higher k strongly lowers VAF strict agreement after candidate/context adjustment")
        add_model("shared-site k", label, model("vaf_conflict_plus_candidate_given_evaluable", predictor),
                  "Higher k strongly raises VAF conflict after candidate/context adjustment")

    return {
        "paths": paths,
        "summary": summary,
        "checks": checks,
        "audit": audit_rows,
        "headline": {
            "vaf_delta": vaf_delta,
            "vaf_delta_adjusted": vaf_delta_adjusted,
            "hp_coarse": hp_coarse,
            "hp_read_eval": hp_read_eval,
            "read_candidate": contrast("read_candidate_topology_gt20_vs_1", "read_exact_or_induced_given_evaluable"),
            "vaf_candidate": contrast("vaf_candidate_topology_ge2_vs_1", "vaf_exact_or_induced_given_evaluable"),
        },
    }


def table_spec(
    table_id: str,
    title: str,
    subtitle: str,
    dataset: str,
    source_id: str,
    columns: Sequence[tuple[str, str, str]],
    *,
    default_sort: tuple[str, str] | None = None,
) -> dict[str, Any]:
    spec: dict[str, Any] = {
        "id": table_id,
        "title": title,
        "subtitle": subtitle,
        "dataset": dataset,
        "sourceId": source_id,
        "density": "dense",
        "columns": [
            {"field": field, "label": label, **({"format": fmt} if fmt != "text" else {"type": "text"})}
            for field, label, fmt in columns
        ],
    }
    if default_sort:
        spec["defaultSort"] = {"field": default_sort[0], "direction": default_sort[1]}
    return spec


def build_artifact(
    raw: Mapping[str, Any],
    pyclone: Mapping[str, Any],
    topology: Mapping[str, Any],
    regional: Mapping[str, Any],
    bridge: Mapping[str, Any],
    cause: Mapping[str, Any],
    previous: Mapping[str, Any],
    paths: Mapping[str, Path],
) -> dict[str, Any]:
    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    metric = raw["metric_map"]
    joint = pyclone["primary_joint"]
    separate = pyclone["primary_separate"]
    independent_minor = pyclone["independent_minor"][0]
    independent_minor_near = pyclone["independent_minor"][1]
    model_mode = pyclone["model_mode"]
    analysis_summary = pyclone["analysis_summary"]
    primary_fit_detail = analysis_summary["fit_summaries"][PRIMARY_JOINT_BUNDLE]
    joint_clonal_fraction = primary_fit_detail["per_sample"]["HCC1395"]["clonal_fraction"]
    truth_count = int(raw["provenance"]["inputs"]["truth_vcfs"]["HCC1395"]["biallelic_autosomal_snv_count"])
    pass_union_joint = first(pyclone["joint"], bundle_id="hcc1395_pair_pass_union_joint_main")
    primary_truth_coverage = float(joint["mutations"]) / truth_count
    pass_union_truth_coverage = float(pass_union_joint["mutations"]) / truth_count
    hcc_shape = first(topology["shapes"], dataset="HCC1395")
    dorado_shape = first(topology["shapes"], dataset="HCC1395_DORADO")
    hcc1954_profile = first(pyclone["all_sample_profiles"], dataset="HCC1954")
    hcc1954_sensitivity = first(
        pyclone["sensitivity"], comparison="HCC1954_main_vs_near_integer", sample_id="HCC1954"
    )
    h1437_sensitivity = first(
        pyclone["sensitivity"], comparison="H1437_main_vs_near_integer", sample_id="H1437"
    )
    regional_headline = regional["headline"]
    regional_k_counts = {row["k_stratum"]: row["region_count"] for row in regional["k_counts"]}
    k4_unresolved = regional["k4_unresolved_diagonal"]
    regional_chart_by_key = {(row["k_stratum"], row["endpoint"]): row for row in regional["chart"]}
    regional_metric_by_key = {(row["stratum"], row["metric_id"]): row for row in regional["metrics"]}
    bridge_summary = bridge["summary"]
    cause_headline = cause["headline"]
    cause_read_candidate_or = first(
        cause["audit"], factor="read candidate multiplicity", evidence="multivariable cluster-robust logistic model"
    )
    cause_vaf_candidate_or = first(
        cause["audit"], factor="VAF candidate multiplicity", evidence="multivariable cluster-robust logistic model"
    )
    cause_k_vaf_or = {
        label: first(
            cause["audit"], factor="shared-site k", evidence="multivariable cluster-robust logistic model",
            endpoint="vaf_strict", comparison=label,
        )
        for label in ("k=3 vs k=2", "k=4 vs k=2", "k>=5 vs k=2")
    }
    cause_k_conflict_or = {
        label: first(
            cause["audit"], factor="shared-site k", evidence="multivariable cluster-robust logistic model",
            endpoint="vaf_conflict", comparison=label,
        )
        for label in ("k=3 vs k=2", "k=4 vs k=2", "k>=5 vs k=2")
    }
    h2009_profile = first(pyclone["all_sample_profiles"], dataset="H2009")
    h2009_sensitivity = first(
        pyclone["sensitivity"], comparison="H2009_main_vs_near_integer", sample_id="H2009"
    )
    old_datasets = previous.get("snapshot", {}).get("datasets", {})
    annotation = old_datasets.get("annotation_compare")
    old_claim_ceiling = old_datasets.get("claim_ceiling")
    if not isinstance(annotation, list) or not annotation:
        raise ValueError("Previous artifact lacks annotation_compare evidence")
    if not isinstance(old_claim_ceiling, list) or not old_claim_ceiling:
        raise ValueError("Previous artifact lacks topology claim_ceiling evidence")

    old_by_layer = {row["layer"]: row for row in old_claim_ceiling}
    for required_layer in ("粗拓撲類別", "Exact T 集合", "VAF-supported selected T／shape"):
        if required_layer not in old_by_layer:
            raise ValueError(f"Previous artifact lacks topology layer: {required_layer}")

    evidence_ladder = [
        {
            "endpoint": "Raw VAF exact-locus",
            "denominator": int(metric["callset_jaccard"]["n_shared"]),
            "agreement": f"CCC={num(metric['vaf_ccc']['value'])}; |delta|≤0.10={pct(metric['vaf_within_0.10']['value'])}",
            "independence": "same cell line, different technical source; shared caller family",
            "claim_ceiling": "strong technical frequency reproducibility; not clone identity",
        },
        {
            "endpoint": "Marginal VAF spectrum",
            "denominator": 7,
            "agreement": "HCC pair is mutual rank-1 nearest neighbor on JSD/KS/Wasserstein",
            "independence": "truth-confirmed sites; cell-line-specific truth sets differ",
            "claim_ceiling": "profile similarity only; marginal distributions can collide across cell lines",
        },
        {
            "endpoint": "Coarse regional topology",
            "denominator": 5720,
            "agreement": old_by_layer["粗拓撲類別"]["evidence"],
            "independence": "historical layered-v2; same method family",
            "claim_ceiling": old_by_layer["粗拓撲類別"]["allowed_claim"],
        },
        {
            "endpoint": "Mutation-labeled exact T set",
            "denominator": 5720,
            "agreement": old_by_layer["Exact T 集合"]["evidence"],
            "independence": "historical layered-v2; candidate-set comparison",
            "claim_ceiling": old_by_layer["Exact T 集合"]["allowed_claim"],
        },
        {
            "endpoint": "VAF-supported selected T/shape",
            "denominator": 5720,
            "agreement": old_by_layer["VAF-supported selected T／shape"]["evidence"],
            "independence": "same-read VAF heuristic; circular as independent validation",
            "claim_ceiling": old_by_layer["VAF-supported selected T／shape"]["allowed_claim"],
        },
        {
            "endpoint": "PyClone-VI joint CCF clusters",
            "denominator": joint["mutations"],
            "agreement": (
                f"shared-label model: clusters={joint['clusters']}; CP Spearman={num(joint['cp_spearman'])}; "
                f"independent minor-group Jaccard={num(independent_minor['subclonal_mutation_jaccard'])}"
            ),
            "independence": "external implementation, but shared HCC CN/purity and overlapping calls",
            "claim_ceiling": "joint labels are shared by construction; independent fits support only moderate minor-group reproducibility; no tree",
        },
    ]

    summary_metric_row = {
        "shared_loci": metric["callset_jaccard"]["n_shared"],
        "callset_jaccard": metric["callset_jaccard"]["value"],
        "vaf_ccc": metric["vaf_ccc"]["value"],
        "within_0_10": metric["vaf_within_0.10"]["value"],
        "noise_rmse_ratio": metric["binomial_noise_rmse_ratio"]["value"],
        "independent_subclone_jaccard": independent_minor["subclonal_mutation_jaccard"],
        "independent_state_kappa": independent_minor["clonal_subclonal_state_kappa"],
        "both_subclonal_ari": independent_minor["both_subclonal_ari"],
    }

    sources = [
        {
            "id": "raw_vaf_validation",
            "label": "chr1-22 raw VAF technical-reproducibility analysis",
            "path": repo_relative(paths["raw_dir"] / "provenance_and_claim_ceiling.json"),
            "query": {
                "engine": "python + bcftools",
                "sql": "SELECT * FROM snapshot.raw_vaf_validation",
                "language": "python",
                "description": "Exact allele joins, caller-count reconstruction, truth intersections, 1 Mb block bootstrap, and all-sample VAF spectra.",
                "tables_used": [repo_relative(paths["raw_dir"] / relative) for relative in RAW_REQUIRED.values()],
                "filters": ["chr1-22", "FILTER=PASS biallelic sSNV", "latest LongPhase-S recalibrated PASS"],
                "metric_definitions": [
                    "ALT count = round-half-up(DP×AF); REF count = DP-ALT; FORMAT/AD[0] is not used",
                    "exact locus = CHROM:POS:REF:ALT",
                    "HCC pair intervals use 500 nonparametric 1 Mb genomic-block bootstrap replicates",
                ],
            },
        },
        {
            "id": "pairwise_vaf_distribution",
            "label": "All-sample truth-confirmed pairwise VAF distance analysis",
            "path": repo_relative(paths["raw_dir"] / "pairwise_distributions" / "pairwise_distribution_analysis.json"),
            "query": {
                "engine": "python",
                "sql": "SELECT * FROM snapshot.pairwise_vaf_distribution",
                "language": "python",
                "description": "Pairwise 50-bin Jensen-Shannon divergence, VAF Wasserstein distance, KS statistic, and nearest-neighbor ranking.",
                "tables_used": [repo_relative(paths["raw_dir"] / RAW_REQUIRED["pairwise"])],
                "metric_definitions": ["Lower distance means more similar marginal raw-VAF spectrum; it does not establish clone identity"],
            },
        },
        {
            "id": "pyclone_analysis",
            "label": "PyClone-VI 0.2.0 CN/purity-aware clustering analysis",
            "path": repo_relative(paths["pyclone_dir"] / "provenance.json"),
            "query": {
                "engine": "PyClone-VI 0.2.0 + Python",
                "sql": "SELECT * FROM snapshot.pyclone_analysis",
                "language": "python",
                "description": "Beta-binomial joint/separate fits, cluster CP profiles, assignment QA, main/near-integer sensitivity, and explicit CN gates.",
                "tables_used": [repo_relative(paths["pyclone_dir"] / relative) for relative in PYCLONE_REQUIRED.values()],
                "filters": ["truth-confirmed PASS sSNVs", "SAVANA allele-specific CN where accepted", "14/14 requested fits PASS"],
                "metric_definitions": [
                    "cellular prevalence is model-conditional on tumor content, major/minor CN, and count model",
                    "PyClone-VI cluster output is not an evolutionary tree",
                    "DORADO uses shared HCC1395 CN as a sensitivity assumption, not an independent CN confirmation",
                ],
            },
        },
        {
            "id": "integrated_topology_context",
            "label": "Seven-dataset final-shape and fixed-denominator HCC topology context",
            "path": repo_relative(paths["topology_context"]),
            "query": {
                "engine": "validated integrated snapshot",
                "sql": "SELECT * FROM snapshot.integrated_topology_context",
                "description": "Seven-dataset final-shape census and fixed 5,720-region coarse/read/VAF/exact/shape evidence ladder.",
                "filters": ["GRCh38 chr1-22", "historical layered-v2 engineering snapshot", "exact-coordinate complete-both HCC regions"],
                "metric_definitions": [
                    "four-class percentages use VAF-resolved final single-shape regions",
                    "unresolved percentage uses all complete regions",
                    "fixed agreement always divides by 5,720; conditional agreement divides by evaluable subset",
                ],
            },
        },
        {
            "id": "regional_precision_analysis",
            "label": "Fixed 5,720-region k-stratified precision validation",
            "path": repo_relative(regional["paths"]["summary"]),
            "query": {
                "engine": "Python fixed-population audit",
                "sql": "SELECT * FROM snapshot.regional_precision_analysis",
                "language": "python",
                "description": "Exact-site, allele, coarse, read-constraint, VAF-constraint, shape, and mutation-labeled forest endpoints stratified by caller_shared_k.",
                "tables_used": [
                    repo_relative(regional["paths"][key])
                    for key in ("summary", "metrics", "k_counts", "checks", "confusion")
                ],
                "filters": [
                    "GRCh38 chr1-22",
                    "5,720 exact-coordinate complete-both HCC regions",
                    "k defined as caller_shared_k",
                    "32/32 conservation checks PASS",
                ],
                "metric_definitions": [
                    "global fixed denominator is always 5,720",
                    "k-stratum fixed denominator includes evaluable and not-evaluable regions within that k stratum",
                    "conditional denominator is endpoint-specific evaluable regions",
                    "k>=4 coarse agreement is inflated by 888 unresolved-unresolved diagonal regions, not by Direct-only recovery",
                ],
            },
        },
        {
            "id": "clone_region_bridge",
            "label": "Independent PyClone mutation-cluster to regional read-edge bridge",
            "path": repo_relative(bridge["paths"]["summary"]),
            "query": {
                "engine": "Python exact-allele bridge audit",
                "sql": "SELECT * FROM snapshot.clone_region_bridge",
                "language": "python",
                "description": "Joins independent separate-fit PyClone assignments to 15,713 regional alleles, region partitions, and singleton read directions without using joint-fit labels as validation.",
                "tables_used": [repo_relative(bridge["paths"][key]) for key in CLONE_BRIDGE_REQUIRED],
                "filters": ["fixed 5,720 regions", "exact mutation key", "14/14 hard checks PASS", "three WARN diagnostics retained"],
                "metric_definitions": [
                    "regional subclonal Jaccard uses the union of independent-fit subclonal mutation labels",
                    "partition exact is label-invariant and must be decomposed into both-single and both-multicluster regions",
                    "directed-edge compatibility uses singleton read directions and cellular-prevalence ordering with tolerance 0.02",
                ],
            },
        },
        {
            "id": "cause_decomposition_audit",
            "label": "Regional topology discordance cause-decomposition red-team audit",
            "path": repo_relative(cause["paths"]["summary"]),
            "query": {
                "engine": "Python block-bootstrap + cluster-robust logistic diagnostics",
                "sql": "SELECT * FROM snapshot.cause_decomposition_audit",
                "language": "python",
                "description": "Separates site sampling, k, HP mapping, candidate multiplicity, depth, VAF delta, release coverage, and shared CN/LOH associations across fixed and evaluable endpoints.",
                "tables_used": [repo_relative(cause["paths"][key]) for key in CAUSE_AUDIT_REQUIRED],
                "filters": ["fixed 5,720 regions", "16/16 checks PASS", "500 chromosome+1Mb block bootstraps", "diagnostic not causal"],
                "metric_definitions": [
                    "unadjusted contrasts use rate A minus rate B and preserve endpoint denominators",
                    "multivariable odds ratios use cluster-robust standard errors by chromosome+1Mb block",
                    "same-read VAF/depth associations and shared CN are diagnostics, not independent validation",
                ],
            },
        },
        {
            "id": "historical_topology",
            "label": "Historical layered-v2 HCC1395 topology and VAF-selection artifact",
            "path": repo_relative(paths["previous_artifact"]),
            "query": {
                "engine": "validated artifact snapshot",
                "sql": "SELECT * FROM snapshot.historical_topology",
                "description": "Imports only the previously validated topology evidence ladder and annotation null-test rows.",
                "filters": ["exact-coordinate complete-both HCC pair regions", "historical layered-v2 engineering snapshot"],
            },
        },
        {
            "id": "pyclone_official",
            "label": "PyClone-VI official repository and input/output contract",
            "href": "https://github.com/Roth-Lab/pyclone-vi",
        },
        {
            "id": "seqc2_hcc1395",
            "label": "SEQC2 HCC1395 reference study (Fang et al.)",
            "href": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/",
        },
    ]

    cards = [
        {
            "id": "raw_pair_card",
            "description": "Latest PASS exact-allele technical comparison.",
            "dataset": "headline_metrics",
            "sourceId": "raw_vaf_validation",
            "metrics": [
                {"label": "Shared exact loci", "field": "shared_loci", "format": "number"},
                {"label": "Callset Jaccard", "field": "callset_jaccard", "format": "percent"},
            ],
        },
        {
            "id": "vaf_pair_card",
            "description": "Site-level frequency agreement across shared exact alleles.",
            "dataset": "headline_metrics",
            "sourceId": "raw_vaf_validation",
            "metrics": [
                {"label": "VAF CCC", "field": "vaf_ccc", "format": "number"},
                {"label": "|delta| ≤ 0.10", "field": "within_0_10", "format": "percent"},
            ],
        },
        {
            "id": "noise_card",
            "description": "Observed paired RMSE relative to a simple independent-binomial lower bound.",
            "dataset": "headline_metrics",
            "sourceId": "raw_vaf_validation",
            "metrics": [{"label": "Noise RMSE ratio", "field": "noise_rmse_ratio", "format": "number"}],
        },
        {
            "id": "pyclone_card",
            "description": "Independent HCC fits, minor-group-focused metrics; avoids the clonal-majority and shared-label inflation.",
            "dataset": "headline_metrics",
            "sourceId": "pyclone_analysis",
            "metrics": [
                {"label": "Subclone Jaccard", "field": "independent_subclone_jaccard", "format": "number"},
                {"label": "State κ", "field": "independent_state_kappa", "format": "number"},
                {"label": "Both-subclonal ARI", "field": "both_subclonal_ari", "format": "number"},
            ],
        },
    ]

    charts = [
        {
            "id": "vaf_distribution_chart",
            "title": "7 datasets 的 truth-confirmed raw VAF spectrum",
            "subtitle": "Latest LongPhase-S PASS；每個 dataset 以 0.02-bin 內比例表示，不把 VAF band 命名為 clone。",
            "type": "line",
            "dataset": "vaf_histogram",
            "sourceId": "raw_vaf_validation",
            "intent": "distribution",
            "question": "HCC pair 的整體頻率形狀是否相近，其他 cell lines 又落在哪裡？",
            "rationale": "同一固定 binning 的 normalized lines 顯示 spectrum 形狀，同時保留不同樣本量的可比性。",
            "comparisonContext": {"denominator": "truth-confirmed PASS sSNVs within each dataset", "grain": "0.02 VAF bin", "unit": "fraction"},
            "encodings": {
                "x": {"field": "bin_mid", "type": "quantitative", "label": "Raw VAF"},
                "y": {"field": "fraction", "type": "quantitative", "label": "Within-dataset fraction", "format": "percent"},
                "color": {"field": "sample", "type": "nominal", "label": "Dataset"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
        {
            "id": "vaf_band_chart",
            "title": "Raw VAF band composition",
            "subtitle": "每個 dataset 內歸一化為 100%；band 只是描述性分箱，不是 clonal/subclonal label。",
            "type": "stackedBar100",
            "dataset": "vaf_bands",
            "sourceId": "raw_vaf_validation",
            "intent": "composition",
            "question": "不同 datasets 的低、中、高 VAF 份額如何分布？",
            "rationale": "100% stacked bars expose composition shifts without being dominated by different site counts.",
            "comparisonContext": {"denominator": "truth-confirmed PASS sSNVs within each dataset", "grain": "VAF band", "normalization": "within dataset", "unit": "share"},
            "encodings": {
                "x": {"field": "sample", "type": "nominal", "label": "Dataset"},
                "y": {"field": "count", "type": "quantitative", "label": "Sites"},
                "color": {"field": "band", "type": "ordinal", "label": "Raw VAF band"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
        {
            "id": "hcc_exact_scatter_chart",
            "title": "HCC pair shared exact loci 的 raw VAF",
            "subtitle": f"Deterministic bounded snapshot n={len(raw['scatter']):,}；完整 effect sizes 與 CI 使用全部 shared loci。兩軸同為 0–1。",
            "type": "scatter",
            "dataset": "hcc_exact_scatter",
            "sourceId": "raw_vaf_validation",
            "intent": "relationship",
            "question": "相同 allele 在兩個技術來源中的 VAF 是否靠近 identity diagonal？",
            "rationale": "A bounded deterministic scatter shows locus-level spread while full-data metrics prevent visual subsampling from defining the effect size.",
            "comparisonContext": {"denominator": "latest PASS shared exact alleles", "grain": "CHROM:POS:REF:ALT", "unit": "raw VAF"},
            "encodings": {
                "x": {"field": "HCC1395_VAF", "type": "quantitative", "label": "HCC1395 VAF"},
                "y": {"field": "HCC1395_DORADO_VAF", "type": "quantitative", "label": "HCC1395_DORADO VAF"},
                "color": {"field": "truth_scope", "type": "nominal", "label": "Truth scope"},
                "tooltip": [
                    {"field": "locus", "type": "text", "label": "Locus"},
                    {"field": "HCC1395_DP", "type": "quantitative", "label": "HCC DP"},
                    {"field": "HCC1395_DORADO_DP", "type": "quantitative", "label": "DORADO DP"},
                ],
            },
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
            "maxRows": len(raw["scatter"]),
        },
        {
            "id": "final_shape_composition_chart",
            "title": "7 dataset rows 的 resolved final-shape composition",
            "subtitle": "Single/Sister/Direct/Sister+direct 在各 dataset 的 final single-shape regions 內歸一化為 100%；unresolved 另以 complete regions 為分母列於表中。",
            "type": "stackedBar100",
            "dataset": "topology_shape_composition",
            "sourceId": "integrated_topology_context",
            "intent": "composition",
            "question": "七個 datasets 的粗 shape 組成如何，HCC pair 的 Direct 主成分是否相近？",
            "rationale": "A resolved-only 100% composition preserves the four-class contract, while the adjacent table keeps unresolved yield visible with its distinct denominator.",
            "comparisonContext": {"denominator": "final single-shape regions per dataset", "grain": "dataset × resolved shape", "normalization": "within dataset", "unit": "share"},
            "encodings": {
                "x": {"field": "dataset", "type": "nominal", "label": "Dataset"},
                "y": {"field": "count", "type": "quantitative", "label": "Resolved regions"},
                "color": {"field": "category", "type": "nominal", "label": "Final shape"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
        {
            "id": "hcc_bland_altman_chart",
            "title": "HCC pair VAF 差值與平均值",
            "subtitle": "DORADO−HCC；虛線為全資料 mean bias，零線表示無系統偏移。",
            "type": "scatter",
            "dataset": "hcc_exact_scatter",
            "sourceId": "raw_vaf_validation",
            "intent": "relationship",
            "question": "VAF 差值是否有系統偏移或隨平均 VAF 改變？",
            "rationale": "Bland–Altman coordinates separate systematic bias from correlation driven by the broad VAF range.",
            "comparisonContext": {"denominator": "latest PASS shared exact alleles", "grain": "CHROM:POS:REF:ALT", "unit": "VAF difference"},
            "encodings": {
                "x": {"field": "mean_VAF", "type": "quantitative", "label": "Mean VAF"},
                "y": {"field": "DORADO_minus_HCC", "type": "quantitative", "label": "DORADO − HCC"},
                "color": {"field": "truth_scope", "type": "nominal", "label": "Truth scope"},
            },
            "referenceLines": [
                {"axis": "y", "value": 0, "label": "No bias", "lineStyle": "solid", "color": "neutral"},
                {
                    "axis": "y",
                    "value": metric["vaf_mean_signed_delta_dorado_minus_hcc1395"]["value"],
                    "label": "Mean bias",
                    "lineStyle": "dashed",
                    "color": "neutral",
                },
            ],
            "palette": {"kind": "categorical"},
            "maxRows": len(raw["scatter"]),
        },
        {
            "id": "pyclone_joint_chart",
            "title": "Primary joint PyClone-VI cluster cellular prevalence",
            "subtitle": "HCC1395 與 DORADO 共用 mutation clusters；DORADO 沿用 HCC1395 CN/purity，故屬 shared-CN sensitivity。",
            "type": "bar",
            "dataset": "pyclone_primary_profiles",
            "sourceId": "pyclone_analysis",
            "intent": "comparison",
            "question": "在同一 joint clustering 下，兩技術來源的 cluster CP 是否相近？",
            "rationale": "Grouped bars preserve cluster identity and compare the two source-specific model estimates directly.",
            "comparisonContext": {"denominator": "truth-confirmed shared primary PASS sSNVs", "grain": "PyClone-VI cluster × technical source", "unit": "cellular prevalence"},
            "encodings": {
                "x": {"field": "cluster_id", "type": "nominal", "label": "Cluster"},
                "y": {"field": "cellular_prevalence", "type": "quantitative", "label": "Cellular prevalence"},
                "color": {"field": "sample_id", "type": "nominal", "label": "Technical source"},
                "tooltip": [
                    {"field": "mutation_count", "type": "quantitative", "label": "Mutations"},
                    {"field": "mutation_fraction", "type": "quantitative", "label": "Mutation fraction"},
                    {"field": "assignment_probability_mean", "type": "quantitative", "label": "Mean assignment probability"},
                    {"field": "assignment_probability_lt_0_8_fraction", "type": "quantitative", "label": "Assignment P<0.8 fraction"},
                ],
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
        {
            "id": "all_sample_cluster_chart",
            "title": "可執行 datasets 的 primary PyClone-VI active clusters",
            "subtitle": "HCC pair 使用 joint primary fit；H1437/H2009/HCC1954 使用 individual main fit；COLO829/HCC1937 因 CN gate 不畫值。",
            "type": "horizontalBar",
            "dataset": "pyclone_ready_profiles",
            "sourceId": "pyclone_analysis",
            "intent": "comparison",
            "question": "在各自主要模型中，哪些 datasets 顯示多少個 active clusters？",
            "rationale": "A compact count view communicates model complexity while the adjacent table preserves CP composition and blocked states.",
            "comparisonContext": {"denominator": "primary fit per dataset", "grain": "dataset", "unit": "active clusters"},
            "encodings": {
                "x": {"field": "dataset", "type": "nominal", "label": "Dataset"},
                "y": {"field": "active_clusters", "type": "quantitative", "label": "Active clusters"},
            },
            "valueFormat": "number",
            "palette": {"kind": "categorical"},
        },
        {
            "id": "regional_k_resolution_chart",
            "title": "多位點 k 增加時，coarse 與細粒度 endpoints 走向分離",
            "subtitle": "Coarse 使用各 k 全部 regions；Read/VAF strict+induced 與 VAF conflict 使用各自 evaluable denominator。k=1 為 0 regions，故不畫。",
            "type": "bar",
            "dataset": "regional_precision_chart",
            "sourceId": "regional_precision_analysis",
            "intent": "comparison",
            "question": "k=2、k=3、k>=4 時，粗分類、read/VAF 子結構與 VAF conflict 如何改變？",
            "rationale": "Grouped bars expose the opposite movement of coarse agreement and fine-grained structural endpoints while retaining denominator type in each row and tooltip.",
            "comparisonContext": {
                "baseline": "k=2",
                "denominator": "coarse: k-stratum fixed; read/VAF: endpoint-evaluable conditional",
                "grain": "caller_shared_k stratum × endpoint",
                "unit": "agreement or conflict rate",
            },
            "encodings": {
                "x": {"field": "k_stratum", "type": "ordinal", "label": "Shared-site complexity"},
                "y": {"field": "value", "type": "quantitative", "label": "Rate", "format": "percent"},
                "color": {"field": "endpoint", "type": "nominal", "label": "Endpoint"},
                "tooltip": [
                    {"field": "numerator", "type": "quantitative", "label": "Numerator"},
                    {"field": "stratum_fixed_denominator", "type": "quantitative", "label": "k fixed denominator"},
                    {"field": "evaluable_denominator", "type": "quantitative", "label": "Evaluable denominator"},
                    {"field": "denominator_basis", "type": "text", "label": "Displayed denominator"},
                ],
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
    ]

    tables = [
        table_spec(
            "vaf_summary_table",
            "7 datasets 的 latest truth-confirmed VAF summary",
            "Raw VAF 與 depth 描述；不能直接把峰值命名為 clone。",
            "vaf_summary",
            "raw_vaf_validation",
            [
                ("sample", "Dataset", "text"),
                ("n", "Truth-confirmed sSNVs", "number"),
                ("vaf_mean", "Mean VAF", "number"),
                ("vaf_median", "Median VAF", "number"),
                ("dp_median", "Median DP", "number"),
            ],
            default_sort=("vaf_median", "desc"),
        ),
        table_spec(
            "hcc_metrics_table",
            "HCC pair full-data exact-locus effect sizes",
            "95% CI 由 500 次 1 Mb genomic-block bootstrap 得到。",
            "hcc_metrics",
            "raw_vaf_validation",
            [
                ("metric", "Metric", "text"),
                ("value", "Estimate", "number"),
                ("ci95_low_1mb_block_bootstrap", "95% CI low", "number"),
                ("ci95_high_1mb_block_bootstrap", "95% CI high", "number"),
                ("n_shared", "Shared loci", "number"),
            ],
        ),
        table_spec(
            "neighbor_table",
            "HCC1395 的全樣本 VAF-spectrum 最近鄰",
            "距離是 marginal distribution 指標；不能用來識別同一 clone。",
            "vaf_neighbors",
            "pairwise_vaf_distribution",
            [
                ("metric", "Distance", "text"),
                ("hcc1395_nearest_neighbor", "Nearest neighbor", "text"),
                ("dorado_distance", "HCC–DORADO distance", "number"),
                ("nearest_other_cell_line", "Nearest other cell line", "text"),
                ("nearest_other_over_dorado_distance_ratio", "Nearest-other / DORADO", "number"),
                ("median_other_over_dorado_distance_ratio", "Median-other / DORADO", "number"),
            ],
        ),
        table_spec(
            "joint_table",
            "HCC pair joint-fit concordance",
            "Cluster-CP profile JSD 是 mutation-weighted CP profile distance，不稱 clone abundance。",
            "pyclone_joint",
            "pyclone_analysis",
            [
                ("bundle_id", "Bundle", "text"),
                ("mutations", "Mutations", "number"),
                ("clusters", "Clusters", "number"),
                ("weighted_mean_abs_cp_delta", "Weighted CP MAE", "number"),
                ("cp_pearson", "CP Pearson", "number"),
                ("cp_spearman", "CP Spearman", "number"),
                ("mutation_weighted_cluster_cp_profile_jsd_bits", "CP-profile JSD (bits)", "number"),
                ("clonal_state_agreement", "Clonal-state agreement", "percent"),
            ],
        ),
        table_spec(
            "cluster_profile_table",
            "Primary joint clusters 的 CP、mutation fraction 與 assignment confidence",
            "Minor clusters 的 mean assignment probability 可明顯低於主群；群數與 CP 必須連同 assignment uncertainty 解讀。",
            "pyclone_primary_profiles",
            "pyclone_analysis",
            [
                ("sample_id", "Technical source", "text"),
                ("cluster_id", "Cluster", "number"),
                ("mutation_count", "Mutations", "number"),
                ("mutation_fraction", "Mutation fraction", "percent"),
                ("cellular_prevalence", "Cellular prevalence", "percent"),
                ("cellular_prevalence_std", "CP SD", "number"),
                ("assignment_probability_mean", "Mean assignment P", "number"),
                ("assignment_probability_median", "Median assignment P", "number"),
                ("assignment_probability_lt_0_8_fraction", "Assignment P<0.8", "percent"),
                ("clone_class", "Model class", "text"),
            ],
            default_sort=("mutation_fraction", "desc"),
        ),
        table_spec(
            "separate_table",
            "HCC pair independent-fit cluster comparison",
            "Separate fits remove joint cluster labels; Hungarian agreement and ARI/NMI quantify post-hoc partition concordance.仍非生物 truth。",
            "pyclone_separate",
            "pyclone_analysis",
            [
                ("scenario", "Scenario", "text"),
                ("mutations_intersection", "Shared mutations", "number"),
                ("clusters_a", "HCC clusters", "number"),
                ("clusters_b", "DORADO clusters", "number"),
                ("ari", "ARI", "number"),
                ("nmi", "NMI", "number"),
                ("hungarian_agreement", "Hungarian agreement", "percent"),
                ("cp_mae", "CP MAE", "number"),
                ("cp_spearman", "CP Spearman", "number"),
                ("cluster_mutation_count_profile_jsd_bits", "Mutation-count profile JSD (bits)", "number"),
            ],
        ),
        table_spec(
            "independent_minor_table",
            "Independent fits 的 minor-group reproducibility",
            "整體 98% label agreement 受 95–98% clonal majority 主導；此表聚焦 subclonal 位點與 both-subclonal partitions。",
            "pyclone_independent_minor",
            "pyclone_analysis",
            [
                ("scenario", "Scenario", "text"),
                ("subclonal_mutation_jaccard", "Subclone Jaccard", "number"),
                ("clonal_subclonal_state_kappa", "Binary state κ", "number"),
                ("clonal_subclonal_state_mcc", "Binary state MCC", "number"),
                ("subclonal_state_f1", "Subclone-state F1", "number"),
                ("both_subclonal_ari", "Both-subclonal ARI", "number"),
                ("both_subclonal_hungarian_agreement", "Both-subclonal Hungarian", "percent"),
                ("both_subclonal_mutations", "Both-subclonal n", "number"),
                ("either_subclonal_mutations", "Either-subclonal n", "number"),
            ],
        ),
        table_spec(
            "model_mode_table",
            "Joint vs separate model mode：subclonal mutation set 並不穩定",
            "Joint fitting 共享 cluster structure；與 separate fit 的差距是 modeling-mode sensitivity，不是 biological change。",
            "pyclone_model_mode",
            "pyclone_analysis",
            [
                ("sample", "Technical source", "text"),
                ("mutations_intersection", "Mutations", "number"),
                ("joint_subclonal_mutations", "Joint subclonal", "number"),
                ("separate_subclonal_mutations", "Separate subclonal", "number"),
                ("both_subclonal_mutations", "Overlap", "number"),
                ("subclonal_mutation_jaccard", "Subclone Jaccard", "number"),
                ("ari", "Partition ARI", "number"),
                ("cellular_prevalence_spearman", "CP Spearman", "number"),
            ],
        ),
        table_spec(
            "sensitivity_table",
            "CN rounding / mutation-universe sensitivity",
            "Main 與 near-integer 或 union scenarios 的 partition/CP 差異；低一致代表模型輸入敏感。",
            "pyclone_sensitivity",
            "pyclone_analysis",
            [
                ("comparison", "Comparison", "text"),
                ("sample_id", "Sample", "text"),
                ("mutations_intersection", "Intersection", "number"),
                ("sensitivity_retention_fraction", "Sensitivity / main mutations", "percent"),
                ("clusters_main", "Main clusters", "number"),
                ("clusters_sensitivity", "Sensitivity clusters", "number"),
                ("ari", "ARI", "number"),
                ("nmi", "NMI", "number"),
                ("subclonal_mutation_jaccard", "Subclone Jaccard", "number"),
                ("clonal_subclonal_state_kappa", "State κ", "number"),
                ("hungarian_agreement", "Hungarian agreement", "percent"),
                ("prevalence_mae", "CP MAE", "number"),
                ("prevalence_spearman", "CP Spearman", "number"),
            ],
        ),
        table_spec(
            "sample_profile_table",
            "全樣本 PyClone primary profile 與資料 gate",
            "High-CP fraction 定義為 CP≥0.90 clusters 的 mutation fraction；normalized entropy 描述 cluster mutation-count dispersion。",
            "pyclone_all_sample_profiles",
            "pyclone_analysis",
            [
                ("dataset", "Dataset", "text"),
                ("status", "Status", "text"),
                ("active_clusters", "Active clusters", "number"),
                ("high_CP_mutation_fraction", "High-CP mutation fraction", "percent"),
                ("normalized_cluster_entropy", "Normalized entropy", "number"),
                ("assignment_probability_mean", "Mean assignment P", "number"),
                ("assignment_probability_median", "Median assignment P", "number"),
                ("assignment_probability_lt_0_8_fraction", "Assignment P<0.8", "percent"),
                ("reason", "Gate / note", "text"),
            ],
        ),
        table_spec(
            "final_shape_table",
            "7 dataset final-shape 完整 census",
            "四類百分比以 final single-shape regions 為分母；unresolved % 以 complete regions 為分母，不能直接混為同一 100%。",
            "topology_final_shapes",
            "integrated_topology_context",
            [
                ("dataset", "Dataset", "text"),
                ("primary_regions", "Primary", "number"),
                ("complete_regions", "Complete", "number"),
                ("final_single_shape_regions", "Resolved final shape", "number"),
                ("single_only", "Single n", "number"),
                ("single_only_pct_final_shape", "Single %", "number"),
                ("sister_only", "Sister n", "number"),
                ("sister_only_pct_final_shape", "Sister %", "number"),
                ("direct_only", "Direct n", "number"),
                ("direct_only_pct_final_shape", "Direct %", "number"),
                ("sister_and_direct", "Mixed n", "number"),
                ("sister_and_direct_pct_final_shape", "Mixed %", "number"),
                ("unresolved_regions", "Unresolved n", "number"),
                ("unresolved_pct_complete", "Unresolved % complete", "number"),
            ],
            default_sort=("dataset", "asc"),
        ),
        table_spec(
            "topology_fixed_ladder_table",
            "HCC pair：固定 5,720 regions 的 coarse/read/VAF/exact/shape ladder",
            "Fixed agreement 一律除以 5,720；conditional agreement 只除以 evaluable subset，兩者不可互換。",
            "topology_fixed_ladder",
            "integrated_topology_context",
            [
                ("endpoint", "Endpoint", "text"),
                ("agreement_n", "Agreement n", "number"),
                ("fixed_denominator", "Fixed denominator", "number"),
                ("fixed_agreement", "Fixed agreement", "percent"),
                ("evaluable_n", "Evaluable n", "number"),
                ("conditional_agreement", "Conditional agreement", "percent"),
                ("claim_boundary", "Claim boundary", "text"),
            ],
        ),
        table_spec(
            "evidence_ladder_table",
            "解析度愈細，能支持的 claim 愈窄",
            "Raw frequency、coarse shape、mutation-labeled tree、CCF clustering 必須分開解讀。",
            "evidence_ladder",
            "historical_topology",
            [
                ("endpoint", "Endpoint", "text"),
                ("denominator", "Denominator", "number"),
                ("agreement", "Observed concordance", "text"),
                ("independence", "Independence / circularity", "text"),
                ("claim_ceiling", "Allowed claim", "text"),
            ],
        ),
        table_spec(
            "annotation_table",
            "Cancer-gene/drug strata 未顯示顯著 enrichment",
            "Conditional null preserves chromosome + region-length decile；所有 p>0.05，只能做 context/face validity。",
            "annotation_compare",
            "historical_topology",
            [
                ("feature", "Annotation", "text"),
                ("present_n", "Present n", "number"),
                ("present_agreement", "Present agreement", "percent"),
                ("absent_agreement", "Absent agreement", "percent"),
                ("present_minus_absent_pp", "Difference (pp)", "number"),
                ("conditional_p_two_sided", "Conditional p", "number"),
                ("interpretation", "Interpretation", "text"),
            ],
            default_sort=("conditional_p_two_sided", "asc"),
        ),
        table_spec(
            "fit_receipts_table",
            "14/14 PyClone-VI fit receipts",
            "Report generator fails closed unless every requested bundle is PASS。",
            "pyclone_all_ready",
            "pyclone_analysis",
            [
                ("bundle_id", "Bundle", "text"),
                ("status", "Status", "text"),
                ("mutations", "Mutations", "number"),
                ("samples", "Samples", "text"),
                ("clusters", "Clusters", "number"),
                ("num_restarts", "Restarts", "number"),
                ("seed", "Seed", "number"),
            ],
        ),
        table_spec(
            "regional_precision_audit_table",
            "多位點 k 與解析度完整 audit table",
            "Global fixed=5,720；stratum fixed 含該 k 全部 regions；conditional 只含 endpoint-evaluable regions。k=1 為 0，不以空白假裝可估。",
            "regional_precision_audit",
            "regional_precision_analysis",
            [
                ("stratum", "k stratum", "text"),
                ("metric_group", "Layer", "text"),
                ("metric_label", "Endpoint", "text"),
                ("numerator", "Numerator", "number"),
                ("global_fixed_denominator", "Global fixed", "number"),
                ("stratum_fixed_denominator", "k fixed", "number"),
                ("evaluable_denominator", "Evaluable", "number"),
                ("rate_global_fixed", "Rate / 5,720", "percent"),
                ("rate_stratum_fixed", "Rate / k fixed", "percent"),
                ("rate_evaluable", "Conditional rate", "percent"),
                ("availability", "Availability", "text"),
                ("claim_ceiling", "Claim boundary", "text"),
            ],
        ),
        table_spec(
            "clone_bridge_concordance_table",
            "Clone-region bridge：join、minor concordance 與 partition selection audit",
            "Independent separate fits only；assignment>=0.8 的 perfect Jaccard 因 subclonal union=1 而為 vacuous selection。",
            "clone_bridge_audit",
            "clone_region_bridge",
            [
                ("layer", "Layer", "text"),
                ("endpoint", "Endpoint", "text"),
                ("population", "Population", "text"),
                ("numerator", "Numerator", "number"),
                ("denominator", "Denominator", "number"),
                ("value", "Value (0–1)", "number"),
                ("status", "Status", "text"),
                ("interpretation", "Interpretation", "text"),
            ],
        ),
        table_spec(
            "clone_bridge_edge_table",
            "Read direction → independent PyClone cluster/CP compatibility",
            "596/598 cross-source read directions相同；但多數joined endpoints同群，無法驗證ancestry。Tolerance=0.02。",
            "clone_bridge_edges",
            "clone_region_bridge",
            [
                ("sample", "Sample/scope", "text"),
                ("confidence", "Confidence", "text"),
                ("endpoint", "Endpoint", "text"),
                ("numerator", "Numerator", "number"),
                ("denominator", "Denominator", "number"),
                ("rate", "Rate", "percent"),
                ("conflict", "Conflicts", "number"),
                ("interpretation", "Interpretation", "text"),
            ],
        ),
        table_spec(
            "cause_decomposition_audit_table",
            "Discordance cause decomposition：unadjusted contrasts 與 multivariable diagnostics",
            "Effect 為 A−B proportion；OR/CI 為 cluster-robust model。全部是同reads、selection-conditional diagnostics，不是因果或accuracy證據。",
            "cause_decomposition_audit",
            "cause_decomposition_audit",
            [
                ("factor", "Factor", "text"),
                ("evidence", "Evidence", "text"),
                ("endpoint", "Endpoint", "text"),
                ("comparison", "Comparison", "text"),
                ("numerator_a", "A numerator", "number"),
                ("denominator_a", "A denominator/model n", "number"),
                ("rate_a", "A rate", "percent"),
                ("numerator_b", "B numerator", "number"),
                ("denominator_b", "B denominator", "number"),
                ("rate_b", "B rate", "percent"),
                ("effect", "A−B", "number"),
                ("odds_ratio", "Odds ratio", "number"),
                ("ci95_low", "CI low", "number"),
                ("ci95_high", "CI high", "number"),
                ("adjusted_q", "BH q", "number"),
                ("interpretation", "Interpretation", "text"),
            ],
        ),
    ]

    technical_summary = (
        "## 結論先講：頻率結構高度再現，真實演化樹尚未被證明\n\n"
        f"- **Raw exact-locus：** latest PASS 共享 {int(metric['callset_jaccard']['n_shared']):,} 個 alleles，"
        f"callset Jaccard={num(metric['callset_jaccard']['value'])}、VAF CCC={num(metric['vaf_ccc']['value'])}、"
        f"|delta|≤0.10 為 {pct(metric['vaf_within_0.10']['value'])}。\n"
        f"- **技術底噪：** 觀察 RMSE 為簡單獨立二項抽樣下限的 {num(metric['binomial_noise_rmse_ratio']['value'])} 倍；"
        f"{pct(metric['binomial_95pct_compatible_fraction']['value'])} 位點差落在簡單 binomial 95% band 內。\n"
        f"- **Regional final shape：** HCC/DORADO resolved regions 的 Direct 為 {hcc_shape['direct_only_pct_final_shape']:.2f}%/"
        f"{dorado_shape['direct_only_pct_final_shape']:.2f}%，但 unresolved/complete 為 {hcc_shape['unresolved_pct_complete']:.2f}%/"
        f"{dorado_shape['unresolved_pct_complete']:.2f}%；dominant coarse component 相近，不代表完整 profile 等價。\n"
        f"- **Regional precision：** 固定 5,720 regions 中 site-set equal={pct(regional_headline['site_set_equal']['rate_stratum_fixed'])}，"
        f"共享 allele identity={int(regional['allele_headline']['identity_alleles']):,}/{int(regional['allele_headline']['shared_alleles']):,}=100%。"
        f"k=2/3/>=4 regions={regional_k_counts['k=2']:,}/{regional_k_counts['k=3']:,}/{regional_k_counts['k>=4']:,}；"
        f"coarse 由 {pct(regional_chart_by_key[('k=2', 'Coarse 5-state')]['value'])} 升至 "
        f"{pct(regional_chart_by_key[('k>=4', 'Coarse 5-state')]['value'])}，但 k>=4 有 "
        f"{int(k4_unresolved['regions']):,}/{int(k4_unresolved['stratum_denominator']):,} unresolved↔unresolved 對角，"
        "不可解讀為 Direct 或 exact tree 更準；同時 read/VAF exact-substructure 與 mutation-labeled forest 隨 k 下降。\n"
        f"- **PyClone-VI：** joint fit 的 {pct(joint_clonal_fraction)} 主群與共享 labels 是 model construction 的結果，不能當獨立再現。"
        f"真正較嚴格的 independent-fit minor metrics 為 subclone Jaccard={num(independent_minor['subclonal_mutation_jaccard'])}、"
        f"state κ={num(independent_minor['clonal_subclonal_state_kappa'])}/MCC={num(independent_minor['clonal_subclonal_state_mcc'])}、"
        f"both-subclonal ARI={num(independent_minor['both_subclonal_ari'])}，只支持中度 minor-group 再現。\n"
        f"- **Clone→region bridge：** 14,369/15,713 regional alleles 可接到兩側 independent fits；regional subclone Jaccard/κ="
        f"{num(bridge_summary['regional_concordance']['subclonal_jaccard'])}/{num(bridge_summary['regional_concordance']['clonal_state_kappa'])}。"
        "Assignment>=0.8 的 perfect Jaccard 是 subclonal union=1 的 vacuous selection；partition exact 96.9% 又由 5,007 both-single regions 主導，"
        "both-multicluster 才是 21/34=61.76%。Read direction 596/598 同向仍因 endpoints 多半同群而不能驗 ancestry。\n"
        f"- **Cause audit：** 全域 binomial-compatible={pct(metric['binomial_95pct_compatible_fraction']['value'])}、noise ratio="
        f"{num(metric['binomial_noise_rmse_ratio']['value'])}x 只支持 site sampling；較強 identifiability drivers 是 k、candidate multiplicity 與 HP mapping。"
        f"高−低 |delta VAF| 的 VAF exact+induced 差為 {100 * cause_headline['vaf_delta']['rate_difference_a_minus_b']:.2f} pp，"
        f"但 complexity-adjusted OR={num(cause_headline['vaf_delta_adjusted']['odds_ratio'])} "
        f"(CI {num(cause_headline['vaf_delta_adjusted']['odds_ratio_ci95_low'])}–{num(cause_headline['vaf_delta_adjusted']['odds_ratio_ci95_high'])})，不穩定。\n"
        f"- **Selection ceiling：** primary {int(joint['mutations']):,} 與 PASS-union complete {int(pass_union_joint['mutations']):,} mutations "
        f"只覆蓋 SEQC2 HC {truth_count:,} 的 {pct(primary_truth_coverage)}/{pct(pass_union_truth_coverage)}；"
        "雙側 count/PASS 選擇可能優先保留高 VAF 主群，所以 95–98% clonal mutation fraction 不是細胞組成估計。\n"
        "- **結論上限：** 支持「同一 cell line 跨技術來源的頻率結構訊號部分再現」；"
        "不支持「每區域已準確回復唯一真實演化樹」或「方法 accuracy 已被證明」。"
    )

    blocks = [
        {"id": "title", "type": "markdown", "body": "# HCC1395 跨來源 VAF／CCF 與 subclone 外部驗證"},
        {
            "id": "framework",
            "type": "markdown",
            "body": "**敘述框架：SCQA + evidence ladder。** 先回答是否再現，再逐層區分 raw VAF、CCF clustering、regional topology 與 biological truth。",
        },
        {
            "id": "partial_scope",
            "type": "markdown",
            "body": (
                "**PARTIAL / COMPREHENSIVE-SCOPE WARNING — chr1–22、7 dataset rows 都已納入 raw VAF，但 CN-aware PyClone-VI 只能對 HCC pair、H1437、H2009、HCC1954 執行。** "
                "COLO829/HCC1937 因缺乏通過 gate 的 allele-specific CN 而 fail closed；DORADO 沒有 source-specific CN，沿用 HCC1395 CN 只是 sensitivity。"
                "Clean layered-v3 topology 未納入本 snapshot，topology 部分明示為 historical layered-v2。"
            ),
            "sourceId": "pyclone_analysis",
        },
        {"id": "technical_summary", "type": "markdown", "body": technical_summary},
        {"id": "headline_metrics", "type": "metric-strip", "cardIds": ["raw_pair_card", "vaf_pair_card", "noise_card", "pyclone_card"]},
        {
            "id": "distribution_heading",
            "type": "markdown",
            "body": (
                "## 全樣本 VAF：HCC pair 互為最近鄰，但 marginal similarity 不等於 clone identity\n\n"
                "Latest truth-confirmed 資料中，HCC1395 與 HCC1395_DORADO 在 JSD、Wasserstein 與 KS 都互為 rank-1 最近鄰。"
                "然而 COLO829–H2009 也可出現相近 marginal spectrum，因此這層只能證明頻率 profile 接近，不能單獨證明同一 subclone。"
            ),
            "sourceId": "pairwise_vaf_distribution",
        },
        {"id": "distribution_chart", "type": "chart", "chartId": "vaf_distribution_chart"},
        {
            "id": "distribution_chart_note",
            "type": "markdown",
            "body": "折線高度是每個 dataset 在該 0.02 VAF bin 的位點比例。峰的位置同時受 purity、CN/LOH、mutation multiplicity、depth 與 caller selection 影響。",
        },
        {"id": "band_chart", "type": "chart", "chartId": "vaf_band_chart"},
        {
            "id": "band_note",
            "type": "markdown",
            "body": "Band composition 用來描述 spectrum 形狀。本報告不把 low/mid/high VAF band 重命名為 subclone；只有納入 CN/purity 的 PyClone cellular prevalence 才使用模型條件下的 CCF 語意。",
        },
        {"id": "vaf_summary", "type": "table", "tableId": "vaf_summary_table"},
        {"id": "neighbor", "type": "table", "tableId": "neighbor_table"},
        {
            "id": "final_shape_heading",
            "type": "markdown",
            "body": (
                "## 7 dataset final shape：Direct 是 HCC pair 的共同主成分，但 unresolved yield 不同\n\n"
                f"HCC1395 與 DORADO 在成功 resolved 的 regions 中，Direct 分別為 {hcc_shape['direct_only_pct_final_shape']:.2f}% 與 "
                f"{dorado_shape['direct_only_pct_final_shape']:.2f}%（差 {abs(hcc_shape['direct_only_pct_final_shape'] - dorado_shape['direct_only_pct_final_shape']):.2f} pp），"
                f"顯示 dominant coarse shape 相近；但 unresolved/complete 為 {hcc_shape['unresolved_pct_complete']:.2f}% 與 "
                f"{dorado_shape['unresolved_pct_complete']:.2f}% ，且 Single/Sister/Mixed 比例也不同。"
                "因此可說 dominant substructure profile 部分再現，不能說整體拓撲分布等價。"
            ),
            "sourceId": "integrated_topology_context",
        },
        {"id": "final_shape_chart", "type": "chart", "chartId": "final_shape_composition_chart"},
        {
            "id": "final_shape_note",
            "type": "markdown",
            "body": "100% stacked chart 只比較已得到單一 final shape 的四類 composition；unresolved 的分母是全部 complete regions，必須看下表。這個 denominator separation 避免 DORADO 較高 unresolved 被圖形正規化隱藏。",
        },
        {"id": "final_shape_table_block", "type": "table", "tableId": "final_shape_table"},
        {
            "id": "exact_heading",
            "type": "markdown",
            "body": (
                "## HCC pair 逐位點：一致訊號強，差異略高於簡單抽樣底噪\n\n"
                f"Full-data CCC={num(metric['vaf_ccc']['value'])}，MAE={num(metric['vaf_mae']['value'])}，median |delta|={num(metric['vaf_median_abs_delta']['value'])}。"
                f"平均 DORADO−HCC bias={num(metric['vaf_mean_signed_delta_dorado_minus_hcc1395']['value'], 4)}。"
                f"RMSE ratio={num(metric['binomial_noise_rmse_ratio']['value'])} 表示觀察差異約比簡單 binomial 下限高 {(metric['binomial_noise_rmse_ratio']['value'] - 1) * 100:.1f}%；"
                "可解釋為輕度額外技術/模型變異，不應立即解釋為生物演化分歧。"
            ),
            "sourceId": "raw_vaf_validation",
        },
        {"id": "exact_scatter", "type": "chart", "chartId": "hcc_exact_scatter_chart"},
        {
            "id": "exact_scatter_note",
            "type": "markdown",
            "body": "Scatter 只用 deterministic bounded rows 避免 HTML 過大；所有報告數值與 95% CI 的計算仍使用完整 shared-locus 集合，不用這個取樣圖回推 effect size。",
        },
        {"id": "bland_altman", "type": "chart", "chartId": "hcc_bland_altman_chart"},
        {
            "id": "bland_note",
            "type": "markdown",
            "body": "Bland–Altman 座標把「整體相關高」與「個別位點差多大」分開。零線周邊的散布主要是 depth/caller 與其他技術變異；沒有 CN/purity/multiplicity 不能把某一點的偏移命名為 clone change。",
        },
        {"id": "hcc_metrics", "type": "table", "tableId": "hcc_metrics_table"},
        {
            "id": "pyclone_heading",
            "type": "markdown",
            "body": (
                "## 外部軟體：PyClone-VI 量化條件式 subclone-profile 再現\n\n"
                f"Primary joint fit 使用 {int(joint['mutations']):,} shared truth-confirmed mutations，得到 {joint['clusters']} clusters；"
                f"weighted CP MAE={num(joint['weighted_mean_abs_cp_delta'])}、CP Pearson={num(joint['cp_pearson'])}、Spearman={num(joint['cp_spearman'])}。"
                "但 joint fitting 由模型建構共享 cluster labels，且 95.29% mutations 落在主群；這個近乎完美的 joint concordance 不能當獨立再現。"
                f"獨立 fits 的 minor-group 指標較嚴格：subclone Jaccard={num(independent_minor['subclonal_mutation_jaccard'])}、"
                f"binary state κ/MCC={num(independent_minor['clonal_subclonal_state_kappa'])}/{num(independent_minor['clonal_subclonal_state_mcc'])}、"
                f"both-subclonal ARI={num(independent_minor['both_subclonal_ari'])}。"
                f"Joint→separate model mode 的 subclonal-set Jaccard 只有 {num(model_mode[0]['subclonal_mutation_jaccard'])}（HCC）與 "
                f"{num(model_mode[1]['subclonal_mutation_jaccard'])}（DORADO）；因此結論是 broad dominant component 再現、minor assignment 中度且 model-sensitive。"
                f"Near-integer CN 子集的 independent subclone Jaccard/κ 為 "
                f"{num(independent_minor_near['subclonal_mutation_jaccard'])}/{num(independent_minor_near['clonal_subclonal_state_kappa'])}，"
                "雖略高於 main，仍未升為高度一致，表示跨來源 minor 差異不只由 CN 四捨五入造成。"
            ),
            "sourceId": "pyclone_analysis",
        },
        {"id": "pyclone_joint_chart_block", "type": "chart", "chartId": "pyclone_joint_chart"},
        {
            "id": "pyclone_joint_note",
            "type": "markdown",
            "body": "Joint chart 比較同一 cluster label 在兩來源的 cellular prevalence。Cluster labels 是 joint model 共享的，不是兩個獨立 fits 後自然對上的群；minor clusters 的 assignment probability 也可能只有約 0.47–0.57。這能回答「同一條件模型是否看到類似頻率結構」，不能回答「獨立分析是否找回同一 minor clone」或「哪個 mutation 在演化樹上是祖先」。",
        },
        {"id": "joint_table_block", "type": "table", "tableId": "joint_table"},
        {"id": "cluster_profile_table_block", "type": "table", "tableId": "cluster_profile_table"},
        {"id": "separate_table_block", "type": "table", "tableId": "separate_table"},
        {"id": "independent_minor_table_block", "type": "table", "tableId": "independent_minor_table"},
        {"id": "model_mode_table_block", "type": "table", "tableId": "model_mode_table"},
        {"id": "sensitivity_table_block", "type": "table", "tableId": "sensitivity_table"},
        {
            "id": "all_sample_heading",
            "type": "markdown",
            "body": (
                "## 全樣本 PyClone 狀況：5 dataset rows 可計算，2 rows 因 CN 品質 fail closed\n\n"
                "不同 cell lines 的 cluster count、high-CP mutation fraction 與 entropy 只用於 profile characterization，不做 mutation-cluster identity matching。"
                "HCC pair 才有 same-cell-line 跨來源的再現性語意。"
                f"HCC1954 雖然 input/fit QA PASS 且主模型有 {hcc1954_profile['active_clusters']} clusters，"
                f"但 assignment probability mean/median 只有 {num(hcc1954_profile['assignment_probability_mean'])}/"
                f"{num(hcc1954_profile['assignment_probability_median'])}，"
                f"{pct(hcc1954_profile['assignment_probability_lt_0_8_fraction'])} mutations 的 assignment P<0.8；"
                f"near-integer 只保留 {int(hcc1954_sensitivity['mutations_sensitivity']):,}/{int(hcc1954_sensitivity['mutations_main']):,}="
                f"{pct(hcc1954_sensitivity['sensitivity_retention_fraction'])} mutations，main→near ARI/NMI="
                f"{num(hcc1954_sensitivity['ari'])}/{num(hcc1954_sensitivity['nmi'])}、subclone Jaccard="
                f"{num(hcc1954_sensitivity['subclonal_mutation_jaccard'])}、state κ={num(hcc1954_sensitivity['clonal_subclonal_state_kappa'])}、"
                f"CP Spearman/MAE={num(hcc1954_sensitivity['prevalence_spearman'])}/{num(hcc1954_sensitivity['prevalence_mae'])}，已顯示明顯不穩。"
                f"相對地 H1437 main→near ARI={num(h1437_sensitivity['ari'])}、subclone Jaccard="
                f"{num(h1437_sensitivity['subclonal_mutation_jaccard'])}、state κ={num(h1437_sensitivity['clonal_subclonal_state_kappa'])}，"
                f"H2009 主模型雖有 {pct(h2009_profile['assignment_probability_lt_0_8_fraction'])} assignment P<0.8，"
                f"main→near ARI/subclone Jaccard={num(h2009_sensitivity['ari'])}/{num(h2009_sensitivity['subclonal_mutation_jaccard'])}；"
                "這表示 broad partition 可對 CN 篩選穩定，但 individual posterior confidence 仍依賴輸入。"
                "整體上 data QA PASS ≠ biological cluster confidence，穩健性必須逐樣本判斷。"
            ),
            "sourceId": "pyclone_analysis",
        },
        {"id": "all_sample_chart", "type": "chart", "chartId": "all_sample_cluster_chart"},
        {
            "id": "all_sample_chart_note",
            "type": "markdown",
            "body": "Active-cluster count 同時受 mutation universe、CN discretization、purity、model restarts 與有效樣本數影響；不可把「相同 cluster count」當成「相同 clones」。尤其 HCC1954 的低 assignment confidence 顯示 cluster count 可以通過檔案/數值 QA，卻仍不具穩定生物詮釋。",
        },
        {"id": "sample_profile", "type": "table", "tableId": "sample_profile_table"},
        {
            "id": "topology_heading",
            "type": "markdown",
            "body": "## 與 InterSubMod topology 對齊：粗結構較能再現，exact mutation tree 明顯較弱\n\n固定 5,720 exact-coordinate complete-both regions，coarse topology 為 3,969/5,720=69.39%（kappa=0.497）；Read strict+true-induced 為 1,599/5,720=27.95%，VAF official 為 1,790/5,720=31.29%。未排名 exact candidate-tree set 的 HP-swap tolerant agreement 為 2,014/5,720=35.21%；VAF-unique mutation-labeled exact forest 只有 949/5,720=16.59% fixed（949/2,543=37.32% conditional）。Structure-first + VAF 的 unlabeled shape 可到 3,667/5,720=64.11% fixed（70.96% conditional），但移除了 mutation labels/direction。",
            "sourceId": "integrated_topology_context",
        },
        {"id": "topology_fixed_ladder", "type": "table", "tableId": "topology_fixed_ladder_table"},
        {
            "id": "topology_ladder_note",
            "type": "markdown",
            "body": "同一列同時顯示 fixed 與 conditional 分母，是為了避免只報 evaluable subset 的高比例。VAF 會增加 decisiveness，但 official VAF candidate constraints 也產生 1,234 conflicts；因此 VAF 修正不是單調改善，更不是獨立真值。",
        },
        {
            "id": "regional_precision_heading",
            "type": "markdown",
            "body": (
                "## 多位點 k 與解析度：位點相同不代表 exact tree 能穩定重建\n\n"
                f"固定 5,720 regions 中，k=1 為 **{regional_k_counts['k=1']:,}**、k=2 為 **{regional_k_counts['k=2']:,}**、"
                f"k=3 為 **{regional_k_counts['k=3']:,}**、k>=4 為 **{regional_k_counts['k>=4']:,}**；k 定義是兩來源 exact region 的 `caller_shared_k`。"
                f"Site-set equal 為 {int(regional_headline['site_set_equal']['numerator']):,}/5,720="
                f"{pct(regional_headline['site_set_equal']['rate_stratum_fixed'])}，共享 alleles 為 "
                f"{int(regional['allele_headline']['identity_alleles']):,}/{int(regional['allele_headline']['shared_alleles']):,}=100% identity。"
                f"然而 coarse agreement 從 k=2 的 {pct(regional_chart_by_key[('k=2', 'Coarse 5-state')]['value'])} 升到 k>=4 的 "
                f"{pct(regional_chart_by_key[('k>=4', 'Coarse 5-state')]['value'])}，主要由 k>=4 的 `Topo>1 未定↔未定` "
                f"{int(k4_unresolved['regions']):,}/{int(k4_unresolved['stratum_denominator']):,}="
                f"{pct(int(k4_unresolved['regions']) / int(k4_unresolved['stratum_denominator']))} 對角堆高，**不是 Direct-only 或 exact tree recovery 增強**。"
                f"相反地，Read strict+induced conditional 由 "
                f"{pct(regional_metric_by_key[('k=2', 'read_strict_plus_induced')]['rate_evaluable'])} 降至 "
                f"{pct(regional_metric_by_key[('k>=4', 'read_strict_plus_induced')]['rate_evaluable'])}，VAF strict+induced 由 "
                f"{pct(regional_metric_by_key[('k=2', 'vaf_official_strict_plus_induced')]['rate_evaluable'])} 降至 "
                f"{pct(regional_metric_by_key[('k>=4', 'vaf_official_strict_plus_induced')]['rate_evaluable'])}，而 VAF conflict 升至 "
                f"{pct(regional_metric_by_key[('k>=4', 'vaf_official_conflict')]['rate_evaluable'])}。"
            ),
            "sourceId": "regional_precision_analysis",
        },
        {"id": "regional_precision_chart_block", "type": "chart", "chartId": "regional_k_resolution_chart"},
        {
            "id": "regional_precision_chart_note",
            "type": "markdown",
            "body": "圖中 coarse bar 的分母是該 k 全部 regions；Read/VAF 三組 bars 的分母是各 endpoint 可評估 subset，因此不可只比柱高而忽略 tooltip 的 numerator、k fixed 與 evaluable denominators。Shape 與 mutation-labeled exact forest 留在下方完整表，避免把不同抽象層級塞進同一圖。",
        },
        {"id": "regional_precision_audit_block", "type": "table", "tableId": "regional_precision_audit_table"},
        {
            "id": "clone_region_bridge_heading",
            "type": "markdown",
            "body": (
                "## Clone→region bridge：可連接不等於可驗證 ancestry\n\n"
                "Regional 15,713 exact alleles 中有 **14,369（91.45%）**能同時接到 HCC 與 DORADO 的 independent separate-fit PyClone labels。"
                f"這個 regional subset 的 subclonal Jaccard={num(bridge_summary['regional_concordance']['subclonal_jaccard'])}、"
                f"binary state κ={num(bridge_summary['regional_concordance']['clonal_state_kappa'])}，比全域結果更弱。"
                "Assignment>=0.8 後表面上 Jaccard=1，但 subclonal union 只有 **1 mutation**，是 selection-induced vacuous perfect score。"
                "Multilocus partition exact 為 5,028/5,189=96.9%，其中 **5,007** 是兩側都只有一個 cluster；"
                "真正兩側都有多群的 subset 只有 34 regions，exact partition 為 **21/34=61.76%**。"
            ),
            "sourceId": "clone_region_bridge",
        },
        {"id": "clone_bridge_concordance_block", "type": "table", "tableId": "clone_bridge_concordance_table"},
        {
            "id": "clone_bridge_edge_heading",
            "type": "markdown",
            "body": (
                "### Read direction 高度重現，但 CP/cluster bridge 大多沒有 ancestry 資訊\n\n"
                "在兩來源都得到 singleton read direction 的 598 relation pairs 中，596 同向、2 反向。"
                "然而接到 independent PyClone 後，HCC1395 只有 **19/1,117**、DORADO 只有 **12/514** 是不同 cluster 且 CP ordering compatible，"
                "兩側 conflict 都是 0；其餘 HCC 1,098 與 DORADO 502 條都落在同一 cluster，對 ancestry ordering 沒有資訊。"
                "在 assignment>=0.8 subset，975/975 與 332/332 更全部同群，所以 0 conflict 不能當成 ancestry validation。"
            ),
            "sourceId": "clone_region_bridge",
        },
        {"id": "clone_bridge_edge_block", "type": "table", "tableId": "clone_bridge_edge_table"},
        {
            "id": "cause_sampling_heading",
            "type": "markdown",
            "body": (
                "## Cause audit 第一層：site sampling 很穩，但只解釋逐位點誤差\n\n"
                f"Exact-locus VAF 有 {pct(metric['binomial_95pct_compatible_fraction']['value'])} differences 落在簡單 binomial 95% band，"
                f"觀察 RMSE 只有理想抽樣下限的 {num(metric['binomial_noise_rmse_ratio']['value'])} 倍。"
                "這支持兩技術來源的 site-level sampling 一致，卻不會自動決定 HP mapping、候選樹數、induced substructure 或 mutation ancestry。"
            ),
            "sourceId": "raw_vaf_validation",
        },
        {
            "id": "cause_decomposition_heading",
            "type": "markdown",
            "body": (
                "### 較強的 identifiability drivers 是 k、candidate multiplicity 與 HP mapping\n\n"
                f"HP-count mismatch 有 **703 regions**：coarse agreement 為 {pct(cause_headline['hp_coarse']['rate_a'])}，"
                f"相對 matched 的 {pct(cause_headline['hp_coarse']['rate_b'])} 低 21.21 pp，且 read/VAF fine endpoints 全部 fail closed（0/703 evaluable）。"
                f"Read candidates >20 vs 1 的 exact+induced conditional 差為 {100 * cause_headline['read_candidate']['rate_difference_a_minus_b']:.2f} pp；"
                f"multivariable 每 2x read candidate OR={num(cause_read_candidate_or['odds_ratio'])}。"
                f"VAF candidates >=2 vs 1 差為 {100 * cause_headline['vaf_candidate']['rate_difference_a_minus_b']:.2f} pp；"
                f"每 2x candidate OR={num(cause_vaf_candidate_or['odds_ratio'])}。"
                f"控制 candidate/context 後，VAF strict 的 k=3/4/>=5 vs k=2 OR="
                f"{num(cause_k_vaf_or['k=3 vs k=2']['odds_ratio'])}/{num(cause_k_vaf_or['k=4 vs k=2']['odds_ratio'])}/"
                f"{num(cause_k_vaf_or['k>=5 vs k=2']['odds_ratio'])}；VAF conflict OR="
                f"{num(cause_k_conflict_or['k=3 vs k=2']['odds_ratio'])}/{num(cause_k_conflict_or['k=4 vs k=2']['odds_ratio'])}/"
                f"{num(cause_k_conflict_or['k>=5 vs k=2']['odds_ratio'])}。"
                f"高−低 |delta VAF| 的粗對比為 {100 * cause_headline['vaf_delta']['rate_difference_a_minus_b']:.2f} pp，"
                f"但加入 complexity/candidate 後 OR={num(cause_headline['vaf_delta_adjusted']['odds_ratio'])}、"
                f"CI={num(cause_headline['vaf_delta_adjusted']['odds_ratio_ci95_low'])}–{num(cause_headline['vaf_delta_adjusted']['odds_ratio_ci95_high'])}，不穩定。"
                "Depth 的反直覺負相關也受 complexity/selection confounding；不可解讀成低 depth 更好。CN/LOH 只來自共享 HCC1395 SAVANA context，"
                "沒有 DORADO source-specific CN，因此不是獨立機制驗證。"
            ),
            "sourceId": "cause_decomposition_audit",
        },
        {"id": "cause_decomposition_audit_block", "type": "table", "tableId": "cause_decomposition_audit_table"},
        {"id": "evidence_ladder", "type": "table", "tableId": "evidence_ladder_table"},
        {
            "id": "annotation_heading",
            "type": "markdown",
            "body": "## 癌症基因/藥物資料是生物脈絡，不是 topology truth\n\nGENCODE gene body、COSMIC CGC、DGIdb interaction/approved drug 與 COSMIC CLP strata 都未在 conditional null 下顯示 p<0.05 enrichment。因此這些資料可作 face validity 與後續優先排序，不能回推 tree accuracy，也不是用藥建議。",
            "sourceId": "historical_topology",
        },
        {"id": "annotation", "type": "table", "tableId": "annotation_table"},
        {
            "id": "method_heading",
            "type": "markdown",
            "body": (
                "## 方法與定義\n\n"
                "1. Raw VAF 主分析使用 chr1–22 latest LongPhase-S recalibrated PASS biallelic sSNVs；各 cell line 再與其 truth/benchmark VCF 取 exact-allele intersection。\n"
                "2. 計數由 DP 與 AF 重建：ALT=round-half-up(DP×AF)、REF=DP−ALT；不使用已確認無效的 AD[0]。\n"
                "3. HCC pair 用 CHROM:POS:REF:ALT exact join；95% CI 以 chromosome+1 Mb block 為抽樣單位、500 replicates。\n"
                "4. PyClone-VI 0.2.0 使用 beta-binomial、seed=20260712、最多 40 components。主模型用 SAVANA allele-specific CN/purity；near-integer CN 與 PASS-union 為 sensitivity。\n"
                "5. `C` 是 read-supported mutation-state groups/有向 rooted Steiner structure 的節點群數，不是 PyClone cluster count；`T` 是候選樹結構；`Topo` 是去除部分標籤後的拓撲形狀。\n"
                "6. PyClone-VI 產生 clusters 與 cellular prevalence，**不產生演化樹**。\n"
                "7. Clone-region bridge 只使用兩側 separate fits；joint-fit shared labels 不作獨立 concordance。Read edge 只有 singleton F/R direction 才進 CP ordering diagnostic。\n"
                "8. Cause audit 的 contrasts 使用 500 次 chromosome+1Mb block bootstrap；multivariable logistic models 使用相同 genomic blocks 的 cluster-robust SE。兩者皆為 diagnostic，不是 causal model。"
            ),
            "sourceId": "pyclone_analysis",
        },
        {"id": "fit_receipts", "type": "table", "tableId": "fit_receipts_table"},
        {
            "id": "limitations",
            "type": "markdown",
            "body": (
                "## 限制與為何一致率沒有想像中高\n\n"
                "- **Identifiability ceiling：** marginal VAF/CCF 可由多種 CN、purity、multiplicity 與 tree ordering 組合產生；一組 reads 通常無法唯一識別 exact tree。\n"
                "- **Technical selection：** basecalling、alignment、HP assignment、depth 與 PASS filter 改變可評估 loci/regions，會同時改變分母與 candidate multiplicity。\n"
                f"- **Mutation-universe selection：** PyClone primary/union-complete 只覆蓋 SEQC2 HC 的 {pct(primary_truth_coverage)}/{pct(pass_union_truth_coverage)}，且要求 caller PASS 與雙側完整 counts；低 VAF subclones 可能被優先漏掉，所以 {pct(joint_clonal_fraction)} clonal mutation fraction 不是真實細胞比例。\n"
                "- **Bridge degeneracy：** assignment>=0.8 的 regional subclone Jaccard=1 只有 1 個 subclonal-union mutation；partition exact 96.9% 又由 5,007 both-single regions 主導，兩者都不能代表高度 multi-clone concordance。\n"
                "- **Cause audit 非因果：** k、candidate multiplicity、depth、VAF delta 與 topology endpoints 復用同批 reads 且互相相關；depth 的反直覺方向不能解讀為低 depth 有益，所有 OR/contrast 都只是 conditional diagnostics。\n"
                "- **CN 非獨立：** DORADO 沒有 source-specific CN；shared HCC1395 SAVANA CN/LOH 可做情境分層，卻不能驗證兩來源有相同 CN 機制，也會高估 CN-corrected independence。\n"
                "- **Same-read circularity：** InterSubMod 的 VAF winner 又使用同批 reads，可做 ranking heuristic，不是 orthogonal validation。\n"
                "- **Truth 空白：** 現在沒有這些 exact mutation trees 的 single-cell/multi-region ground truth；Fang/SEQC2 的 HCC1395 heterogeneity 只表示大方向相容。\n"
                "- **外部工具也是模型：** PyClone-VI 的多 restarts 只提高 optimization 穩定，不把 posterior cluster 變成 biological truth。"
            ),
        },
        {
            "id": "next_steps",
            "type": "markdown",
            "body": (
                "## 最短的補強路徑\n\n"
                "1. 對 HCC1395 與 DORADO 各自從 BAM 重數 union loci，並產生 DORADO source-specific allele-specific CN。\n"
                "2. 在每個技術來源內做 read/molecule split-half，先估計 within-source reproducibility ceiling，再解釋 cross-source gap。\n"
                "3. 使用 MOBSTER 或另一種實作獨立的 CCF mixture 作 orthogonal sensitivity；只在 cluster stability 通過後才接 Pairtree，並報 posterior tree set，不只報最高分單樹。\n"
                "4. 以 single-cell DNA/CNV、multi-region 或 synthetic spike-in truth 驗證 mutation→cluster→edge；這一步才能將「再現」升級為「準確」。\n"
                "5. 對 703 個 HP-count mismatch 與高 candidate-multiplicity regions 分開建立 failure benchmark；不要用 coarse unresolved 對角掩蓋 fine-endpoint fail-closed。\n"
                "6. 在 directed-edge bridge 中預先要求足夠的 different-cluster endpoints；目前高信心 subset 全部同群，無法檢驗 CP ancestry ordering。\n"
                "7. Clean layered-v3 topology release 通過 final receipt 後，以 locked classifier 全數重算 evidence ladder。"
            ),
        },
        {
            "id": "further_questions",
            "type": "markdown",
            "body": "## 尚待回答\n\n- 兩來源的 passage/library DNA 是否完全同批，還是存在 culture drift？\n- 在 source-specific CN 下，PyClone joint/separate 的 cluster CP 結論是否維持？\n- 若用 orthogonal/single-cell clusters 取代目前 independent PyClone labels，regional subJ=0.239 與 both-multicluster 21/34 是否維持？\n- HP-count mismatch 的 703 regions 是 basecalling/phasing 技術差異，還是 region definition 本身不穩？\n- 在 matched k/candidate multiplicity 的設計下，depth 與 |delta VAF| 的效果是否仍存在？\n- 哪些 regional topology edges 能跨不同且高信心的 clusters，真正提供 ancestry-ordering 檢驗？",
        },
    ]

    snapshot = {
        "version": 1,
        "generatedAt": generated_at,
        "status": "partial",
        "datasets": {
            "headline_metrics": [summary_metric_row],
            "vaf_summary": raw["summary"],
            "vaf_histogram": raw["histogram"],
            "vaf_bands": raw["bands"],
            "hcc_metrics": raw["metrics"],
            "hcc_exact_scatter": raw["scatter"],
            "vaf_neighbors": raw["neighbor"],
            "pairwise_vaf_distances": raw["pairwise"],
            "pyclone_primary_profiles": pyclone["primary_profiles"],
            "pyclone_joint": pyclone["joint"],
            "pyclone_separate": pyclone["separate"],
            "pyclone_independent_minor": pyclone["independent_minor"],
            "pyclone_model_mode": pyclone["model_mode"],
            "pyclone_sensitivity": pyclone["sensitivity"],
            "pyclone_all_sample_profiles": pyclone["all_sample_profiles"],
            "pyclone_ready_profiles": [row for row in pyclone["all_sample_profiles"] if row["status"] != "BLOCKED"],
            "pyclone_all_ready": pyclone["all_ready"],
            "pyclone_fit_summaries": pyclone["fit_summaries"],
            "pyclone_blocked": pyclone["blocked"],
            "topology_final_shapes": topology["shapes"],
            "topology_shape_composition": topology["composition"],
            "topology_fixed_ladder": topology["ladder"],
            "regional_precision_chart": regional["chart"],
            "regional_precision_audit": regional["metrics"],
            "regional_precision_k_counts": regional["k_counts"],
            "regional_precision_checks": regional["checks"],
            "clone_bridge_audit": bridge["audit"],
            "clone_bridge_edges": bridge["edges"],
            "clone_bridge_checks": bridge["checks"],
            "clone_bridge_warnings": bridge["flags"],
            "cause_decomposition_audit": cause["audit"],
            "cause_decomposition_checks": cause["checks"],
            "evidence_ladder": evidence_ladder,
            "annotation_compare": annotation,
        },
    }

    manifest = {
        "version": 1,
        "surface": "report",
        "title": "HCC1395 跨來源 VAF／CCF 與 subclone 外部驗證",
        "description": "chr1-22 全七 dataset raw VAF、HCC1395 exact-locus technical concordance、PyClone-VI CN/purity-aware clustering、historical topology evidence ladder 與 gene/drug null evidence。",
        "generatedAt": generated_at,
        "cards": cards,
        "charts": charts,
        "tables": tables,
        "sources": sources,
        "blocks": blocks,
    }
    return {"surface": "report", "manifest": manifest, "snapshot": snapshot, "sources": sources}


def main() -> int:
    args = parse_args()
    raw_paths = require_files(args.raw_dir, RAW_REQUIRED, "Raw VAF analysis")
    pyclone_paths = require_files(args.pyclone_analysis_dir, PYCLONE_REQUIRED, "PyClone analysis")
    if not args.previous_artifact.is_file():
        raise FileNotFoundError(f"Previous validated topology artifact is missing: {args.previous_artifact}")
    if not args.topology_context.is_file():
        raise FileNotFoundError(f"Integrated topology context is missing: {args.topology_context}")

    raw = normalize_raw(raw_paths, args.scatter_max_rows)
    pyclone = normalize_pyclone(pyclone_paths)
    topology = normalize_topology(args.topology_context)
    regional = normalize_regional_precision(args.regional_precision_dir)
    bridge = normalize_clone_bridge(args.clone_region_bridge_dir)
    cause = normalize_cause_audit(args.cause_audit_dir)
    previous = read_json(args.previous_artifact)
    artifact = build_artifact(
        raw,
        pyclone,
        topology,
        regional,
        bridge,
        cause,
        previous,
        {
            "raw_dir": args.raw_dir,
            "pyclone_dir": args.pyclone_analysis_dir,
            "previous_artifact": args.previous_artifact,
            "topology_context": args.topology_context,
            "regional_precision_dir": args.regional_precision_dir,
            "clone_region_bridge_dir": args.clone_region_bridge_dir,
            "cause_audit_dir": args.cause_audit_dir,
        },
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": "PASS",
                "output": str(args.output.resolve()),
                "snapshot_status": artifact["snapshot"]["status"],
                "datasets": len(artifact["snapshot"]["datasets"]),
                "blocks": len(artifact["manifest"]["blocks"]),
                "charts": len(artifact["manifest"]["charts"]),
                "tables": len(artifact["manifest"]["tables"]),
                "sources": len(artifact["sources"]),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
