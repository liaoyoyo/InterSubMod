#!/usr/bin/env python3
"""Build deep-dive workspace for TO residual FP diagnosis."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

from build_to_fp_provenance_analysis import (
    DEFAULT_SAMPLES,
    SAMPLE_CONFIGS,
    load_configs,
    load_paired_feature_maps,
    load_to_feature_maps,
    rule_trigger,
)
from research_common import as_bool, compute_metrics, ensure_dir, markdown_table, to_float, write_tsv_rows


DEFAULT_PROVENANCE_WORKSPACE = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260322_to_fp_provenance_analysis"
)
DEFAULT_OUTPUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260323_to_residual_fp_deep_dive"
)

STAGE_LABELS = {
    "paired_raw_absent": "raw_absent",
    "paired_persistent_final_fp": "persistent",
    "paired_raw_artifact_filtered": "raw_filter",
    "paired_longphase_s_removed": "longphase_s",
    "paired_postprocess_removed": "paired_postprocess",
    "germline_filter": "germline_filter",
    "normal_alt_support": "normal_alt_support",
}

FINAL_FEATURE_FIELDS = [
    "qual",
    "gq",
    "dp",
    "af",
    "ad_ref",
    "ad_alt",
    "pairwise_median_dist",
    "cramers_v",
    "quality_score",
    "allele_delta",
    "verification_class",
    "quality_tier",
    "potential_loh",
    "coverage_category",
    "suggest_filter",
    "hp_assign_rate",
    "allele_assign_rate",
    "agreement_type",
    "class_shift",
    "effect_size",
    "cluster_class",
    "label_class",
]

RESIDUAL_FP_FIELDS = [
    "sample",
    "platform",
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "stage_label",
    "paired_resolution_stage",
    "paired_resolution_reason",
    "paired_normal_resolvable",
    "shared_exact_across_platform",
    "shared_same_stage_with_other_platform",
    "shared_same_pattern_with_other_platform",
    "to_caller_filter",
    "to_pon_hit",
    "to_pon_hit_count",
    "to_verdict_somatic",
    "to_verdict_subclonal",
    "to_verdict_germline",
    "to_has_h_flag",
    "to_multihap_flag",
    "to_noancestry_flag",
    "to_caller_qual",
    "to_caller_gq",
    "to_caller_dp",
    "to_caller_af",
    "to_caller_ad_ref",
    "to_caller_ad_alt",
    "qual",
    "gq",
    "dp",
    "af",
    "ad_ref",
    "ad_alt",
    "pairwise_median_dist",
    "cramers_v",
    "quality_score",
    "allele_delta",
    "verification_class",
    "quality_tier",
    "potential_loh",
    "coverage_category",
    "suggest_filter",
    "hp_assign_rate",
    "allele_assign_rate",
    "agreement_type",
    "class_shift",
    "effect_size",
    "cluster_class",
    "label_class",
    "paired_raw_filter",
    "paired_raw_qual",
    "paired_raw_gq",
    "paired_raw_dp",
    "paired_raw_af",
    "paired_raw_naf",
    "paired_raw_ndp",
    "paired_raw_nad_ref",
    "paired_raw_nad_alt",
    "paired_final_qual",
    "paired_final_gq",
    "paired_final_dp",
    "paired_final_af",
    "paired_final_ad_ref",
    "paired_final_ad_alt",
    "paired_final_pairwise_median_dist",
    "paired_final_cramers_v",
    "paired_final_quality_score",
    "paired_final_allele_delta",
    "paired_final_verification_class",
    "paired_final_hp_assign_rate",
    "paired_final_allele_assign_rate",
    "paired_final_agreement_type",
    "paired_final_class_shift",
]

FINAL_TP_FIELDS = [
    "sample",
    "platform",
    "variant_key",
    "qual",
    "gq",
    "dp",
    "af",
    "ad_ref",
    "ad_alt",
    "pairwise_median_dist",
    "cramers_v",
    "quality_score",
    "allele_delta",
    "verification_class",
    "quality_tier",
    "potential_loh",
    "coverage_category",
    "suggest_filter",
    "hp_assign_rate",
    "allele_assign_rate",
    "agreement_type",
    "class_shift",
    "effect_size",
    "cluster_class",
    "label_class",
]

SAMPLE_SUMMARY_FIELDS = [
    "sample",
    "platform",
    "truth_total",
    "to_final_tp",
    "to_final_fp",
    "to_final_fn",
    "to_final_precision",
    "to_final_recall",
    "to_final_f1",
    "raw_absent_count",
    "raw_filter_count",
    "longphase_s_count",
    "persistent_count",
    "raw_absent_shared_exact_count",
    "persistent_shared_exact_count",
]

STAGE_SUMMARY_FIELDS = [
    "sample",
    "stage_label",
    "count",
    "fraction_of_final_fp",
    "ideal_f1_if_perfect_removal",
    "ideal_delta_f1_if_perfect_removal",
]

FEATURE_GROUP_FIELDS = [
    "sample",
    "group_name",
    "count",
    "median_af",
    "median_gq",
    "median_qual",
    "median_dp",
    "median_ad_alt",
    "median_quality_score",
    "median_pairwise_median_dist",
    "median_allele_delta",
    "median_hp_assign_rate",
    "median_allele_assign_rate",
    "top_verification_class",
    "top_agreement_type",
    "top_class_shift",
]

PATTERN_SUMMARY_FIELDS = [
    "sample",
    "stage_label",
    "verification_class",
    "agreement_type",
    "class_shift",
    "count",
    "fraction_of_stage",
    "tp_match_count",
    "fp_to_tp_ratio",
    "shared_exact_count",
    "median_af",
    "median_gq",
    "median_quality_score",
    "median_pairwise_median_dist",
]

RULE_SCAN_FIELDS = [
    "sample",
    "stage_label",
    "scan_family",
    "rule_rank",
    "rule_id",
    "rule_label",
    "target_stage_total",
    "target_stage_removed",
    "other_fp_removed",
    "fp_removed_total",
    "tp_removed",
    "tp_total_before",
    "fp_total_before",
    "tp_after",
    "fp_after",
    "precision_after",
    "recall_after",
    "f1_after",
    "delta_f1_vs_final",
]

RECURRENCE_SUMMARY_FIELDS = [
    "stage_label",
    "sample",
    "stage_total",
    "shared_exact_count",
    "platform_specific_count",
    "fraction_shared_exact",
    "same_pattern_shared_exact_count",
]

SHARED_HOTSPOT_FIELDS = [
    "stage_label",
    "sample_a",
    "sample_b",
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "same_pattern",
    "sample_a_verification_class",
    "sample_a_agreement_type",
    "sample_a_class_shift",
    "sample_a_af",
    "sample_a_gq",
    "sample_a_qual",
    "sample_a_quality_score",
    "sample_a_pairwise_median_dist",
    "sample_a_hp_assign_rate",
    "sample_b_verification_class",
    "sample_b_agreement_type",
    "sample_b_class_shift",
    "sample_b_af",
    "sample_b_gq",
    "sample_b_qual",
    "sample_b_quality_score",
    "sample_b_pairwise_median_dist",
    "sample_b_hp_assign_rate",
]

PERSISTENT_CANDIDATE_FIELDS = [
    "sample",
    "shared_exact_across_platform",
    "shared_same_pattern_with_other_platform",
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "verification_class",
    "agreement_type",
    "class_shift",
    "af",
    "gq",
    "qual",
    "quality_score",
    "pairwise_median_dist",
    "allele_delta",
    "hp_assign_rate",
    "suggest_filter",
    "potential_loh",
    "paired_raw_filter",
    "paired_raw_af",
    "paired_raw_naf",
    "paired_raw_gq",
    "paired_raw_dp",
    "paired_final_verification_class",
    "paired_final_agreement_type",
    "paired_final_class_shift",
    "paired_final_af",
    "paired_final_gq",
    "paired_final_quality_score",
    "paired_final_pairwise_median_dist",
    "paired_final_hp_assign_rate",
]

SINGLE_CATEGORICAL_FEATURES = [
    "verification_class",
    "agreement_type",
    "class_shift",
    "cluster_class",
    "label_class",
    "suggest_filter",
    "potential_loh",
    "coverage_category",
]

PAIR_CATEGORICAL_FEATURES = [
    ("verification_class", "agreement_type"),
    ("verification_class", "class_shift"),
    ("agreement_type", "class_shift"),
]

NUMERIC_FEATURE_DIRECTIONS = {
    "gq": ("le",),
    "qual": ("le",),
    "af": ("le", "ge"),
    "dp": ("le",),
    "ad_alt": ("le",),
    "quality_score": ("le",),
    "pairwise_median_dist": ("le", "ge"),
    "hp_assign_rate": ("le",),
    "allele_assign_rate": ("le",),
    "allele_delta": ("le", "ge"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", default=None, help="Sample name, repeatable")
    parser.add_argument(
        "--provenance-workspace",
        default=str(DEFAULT_PROVENANCE_WORKSPACE),
        help="Workspace from 20260322 TO FP provenance analysis",
    )
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Output workspace directory")
    return parser.parse_args()


def normalize_category(value: object) -> str:
    text = str(value or "").strip()
    return text if text else "(blank)"


def format_float(value: float) -> str:
    if value is None or math.isnan(value):
        return ""
    return f"{value:.6f}"


def median_or_blank(values: Iterable[float]) -> str:
    clean = sorted(value for value in values if value is not None and not math.isnan(value))
    if not clean:
        return ""
    mid = len(clean) // 2
    if len(clean) % 2:
        return format_float(clean[mid])
    return format_float((clean[mid - 1] + clean[mid]) / 2.0)


def top_or_blank(values: Iterable[str]) -> str:
    counter = Counter(value for value in values if str(value or "").strip())
    if not counter:
        return ""
    return counter.most_common(1)[0][0]


def quantile_thresholds(values: Sequence[float]) -> List[float]:
    clean = sorted(value for value in values if value is not None and not math.isnan(value))
    if not clean:
        return []
    quantiles = [0.05, 0.10, 0.15, 0.20, 0.25, 0.33, 0.50, 0.67, 0.75, 0.80, 0.90, 0.95]
    thresholds: List[float] = []
    for q in quantiles:
        idx = min(max(int(round((len(clean) - 1) * q)), 0), len(clean) - 1)
        thresholds.append(round(clean[idx], 6))
    unique: List[float] = []
    for threshold in thresholds:
        if threshold not in unique:
            unique.append(threshold)
    return unique


def write_gzip_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    ensure_dir(path.parent)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_residual_master_rows(path: Path) -> Dict[str, Dict[str, str]]:
    rows: Dict[str, Dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("primary_class") != "to_residual_final_fp":
                continue
            key = row["variant_key"]
            rows[key] = row
    return rows


def build_final_tp_row(sample: str, platform: str, variant_key: str, feature: Dict[str, object]) -> Dict[str, object]:
    row: Dict[str, object] = {"sample": sample, "platform": platform, "variant_key": variant_key}
    for field in FINAL_FEATURE_FIELDS:
        row[field] = feature.get(field, "")
    return row


def build_residual_fp_row(
    sample: str,
    platform: str,
    master_row: Dict[str, str],
    feature: Dict[str, object],
    paired_feature: Dict[str, object],
) -> Dict[str, object]:
    stage_label = STAGE_LABELS.get(master_row["paired_resolution_reason"], master_row["paired_resolution_reason"])
    row: Dict[str, object] = {
        "sample": sample,
        "platform": platform,
        "variant_key": master_row["variant_key"],
        "chrom": master_row["chrom"],
        "pos": master_row["pos"],
        "ref": master_row["ref"],
        "alt": master_row["alt"],
        "stage_label": stage_label,
        "paired_resolution_stage": master_row["paired_resolution_stage"],
        "paired_resolution_reason": master_row["paired_resolution_reason"],
        "paired_normal_resolvable": master_row["paired_normal_resolvable"],
        "shared_exact_across_platform": False,
        "shared_same_stage_with_other_platform": False,
        "shared_same_pattern_with_other_platform": False,
        "to_caller_filter": master_row["to_caller_filter"],
        "to_pon_hit": master_row["to_pon_hit"],
        "to_pon_hit_count": master_row["to_pon_hit_count"],
        "to_verdict_somatic": master_row["to_verdict_somatic"],
        "to_verdict_subclonal": master_row["to_verdict_subclonal"],
        "to_verdict_germline": master_row["to_verdict_germline"],
        "to_has_h_flag": master_row["to_has_h_flag"],
        "to_multihap_flag": master_row["to_multihap_flag"],
        "to_noancestry_flag": master_row["to_noancestry_flag"],
        "to_caller_qual": master_row["to_qual"],
        "to_caller_gq": master_row["to_gq"],
        "to_caller_dp": master_row["to_dp"],
        "to_caller_af": master_row["to_af"],
        "to_caller_ad_ref": master_row["to_ad_ref"],
        "to_caller_ad_alt": master_row["to_ad_alt"],
        "paired_raw_filter": master_row["paired_raw_filter"],
        "paired_raw_qual": master_row["paired_raw_qual"],
        "paired_raw_gq": master_row["paired_raw_gq"],
        "paired_raw_dp": master_row["paired_raw_dp"],
        "paired_raw_af": master_row["paired_raw_af"],
        "paired_raw_naf": master_row["paired_raw_naf"],
        "paired_raw_ndp": master_row["paired_raw_ndp"],
        "paired_raw_nad_ref": master_row["paired_raw_nad_ref"],
        "paired_raw_nad_alt": master_row["paired_raw_nad_alt"],
        "paired_final_qual": paired_feature.get("qual", ""),
        "paired_final_gq": paired_feature.get("gq", ""),
        "paired_final_dp": paired_feature.get("dp", ""),
        "paired_final_af": paired_feature.get("af", ""),
        "paired_final_ad_ref": paired_feature.get("ad_ref", ""),
        "paired_final_ad_alt": paired_feature.get("ad_alt", ""),
        "paired_final_pairwise_median_dist": paired_feature.get("pairwise_median_dist", ""),
        "paired_final_cramers_v": paired_feature.get("cramers_v", ""),
        "paired_final_quality_score": paired_feature.get("quality_score", ""),
        "paired_final_allele_delta": paired_feature.get("allele_delta", ""),
        "paired_final_verification_class": paired_feature.get("verification_class", ""),
        "paired_final_hp_assign_rate": paired_feature.get("hp_assign_rate", ""),
        "paired_final_allele_assign_rate": paired_feature.get("allele_assign_rate", ""),
        "paired_final_agreement_type": paired_feature.get("agreement_type", ""),
        "paired_final_class_shift": paired_feature.get("class_shift", ""),
    }
    for field in FINAL_FEATURE_FIELDS:
        row[field] = feature.get(field, "")
    return row


def build_group_summary(sample: str, group_name: str, rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    return {
        "sample": sample,
        "group_name": group_name,
        "count": len(rows),
        "median_af": median_or_blank(to_float(row.get("af")) for row in rows),
        "median_gq": median_or_blank(to_float(row.get("gq")) for row in rows),
        "median_qual": median_or_blank(to_float(row.get("qual")) for row in rows),
        "median_dp": median_or_blank(to_float(row.get("dp")) for row in rows),
        "median_ad_alt": median_or_blank(to_float(row.get("ad_alt")) for row in rows),
        "median_quality_score": median_or_blank(to_float(row.get("quality_score")) for row in rows),
        "median_pairwise_median_dist": median_or_blank(to_float(row.get("pairwise_median_dist")) for row in rows),
        "median_allele_delta": median_or_blank(to_float(row.get("allele_delta")) for row in rows),
        "median_hp_assign_rate": median_or_blank(to_float(row.get("hp_assign_rate")) for row in rows),
        "median_allele_assign_rate": median_or_blank(to_float(row.get("allele_assign_rate")) for row in rows),
        "top_verification_class": top_or_blank(str(row.get("verification_class", "")) for row in rows),
        "top_agreement_type": top_or_blank(str(row.get("agreement_type", "")) for row in rows),
        "top_class_shift": top_or_blank(str(row.get("class_shift", "")) for row in rows),
    }


def build_pattern_summary(
    sample: str,
    stage_label: str,
    stage_rows: Sequence[Dict[str, object]],
    tp_rows: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    tp_counter = Counter(
        (
            normalize_category(row.get("verification_class")),
            normalize_category(row.get("agreement_type")),
            normalize_category(row.get("class_shift")),
        )
        for row in tp_rows
    )
    grouped: Dict[Tuple[str, str, str], List[Dict[str, object]]] = defaultdict(list)
    for row in stage_rows:
        key = (
            normalize_category(row.get("verification_class")),
            normalize_category(row.get("agreement_type")),
            normalize_category(row.get("class_shift")),
        )
        grouped[key].append(row)

    out_rows: List[Dict[str, object]] = []
    total = len(stage_rows)
    for key, rows in sorted(grouped.items(), key=lambda item: len(item[1]), reverse=True):
        tp_match_count = tp_counter.get(key, 0)
        out_rows.append(
            {
                "sample": sample,
                "stage_label": stage_label,
                "verification_class": key[0],
                "agreement_type": key[1],
                "class_shift": key[2],
                "count": len(rows),
                "fraction_of_stage": format_float(len(rows) / total if total else math.nan),
                "tp_match_count": tp_match_count,
                "fp_to_tp_ratio": format_float(len(rows) / tp_match_count if tp_match_count else math.nan),
                "shared_exact_count": sum(1 for row in rows if as_bool(row.get("shared_exact_across_platform"))),
                "median_af": median_or_blank(to_float(row.get("af")) for row in rows),
                "median_gq": median_or_blank(to_float(row.get("gq")) for row in rows),
                "median_quality_score": median_or_blank(to_float(row.get("quality_score")) for row in rows),
                "median_pairwise_median_dist": median_or_blank(to_float(row.get("pairwise_median_dist")) for row in rows),
            }
        )
    return out_rows


def evaluate_rule(
    sample: str,
    stage_label: str,
    scan_family: str,
    rule_id: str,
    rule_label: str,
    predicate,
    stage_rows: Sequence[Dict[str, object]],
    fp_rows: Sequence[Dict[str, object]],
    tp_rows: Sequence[Dict[str, object]],
    truth_total: int,
    final_f1: float,
) -> Dict[str, object]:
    target_stage_removed = sum(1 for row in stage_rows if predicate(row))
    fp_removed_total = sum(1 for row in fp_rows if predicate(row))
    tp_removed = sum(1 for row in tp_rows if predicate(row))
    metrics = compute_metrics(len(tp_rows) - tp_removed, len(fp_rows) - fp_removed_total, truth_total)
    return {
        "sample": sample,
        "stage_label": stage_label,
        "scan_family": scan_family,
        "rule_rank": 0,
        "rule_id": rule_id,
        "rule_label": rule_label,
        "target_stage_total": len(stage_rows),
        "target_stage_removed": target_stage_removed,
        "other_fp_removed": fp_removed_total - target_stage_removed,
        "fp_removed_total": fp_removed_total,
        "tp_removed": tp_removed,
        "tp_total_before": len(tp_rows),
        "fp_total_before": len(fp_rows),
        "tp_after": metrics["tp"],
        "fp_after": metrics["fp"],
        "precision_after": format_float(metrics["precision"]),
        "recall_after": format_float(metrics["recall"]),
        "f1_after": format_float(metrics["f1"]),
        "delta_f1_vs_final": format_float(metrics["f1"] - final_f1),
    }


def build_stage_rules(
    sample: str,
    stage_label: str,
    stage_rows: Sequence[Dict[str, object]],
    fp_rows: Sequence[Dict[str, object]],
    tp_rows: Sequence[Dict[str, object]],
    truth_total: int,
    final_f1: float,
) -> List[Dict[str, object]]:
    rules: List[Dict[str, object]] = []

    for feature in SINGLE_CATEGORICAL_FEATURES:
        counter = Counter(normalize_category(row.get(feature)) for row in stage_rows)
        for value, count in counter.items():
            if value == "(blank)" or count < 2:
                continue
            rules.append(
                evaluate_rule(
                    sample,
                    stage_label,
                    "categorical_single",
                    f"{feature}__eq__{value}",
                    f"{feature} == {value}",
                    lambda row, feature=feature, value=value: normalize_category(row.get(feature)) == value,
                    stage_rows,
                    fp_rows,
                    tp_rows,
                    truth_total,
                    final_f1,
                )
            )

    for left, right in PAIR_CATEGORICAL_FEATURES:
        counter = Counter((normalize_category(row.get(left)), normalize_category(row.get(right))) for row in stage_rows)
        for (left_value, right_value), count in counter.items():
            if "(blank)" in {left_value, right_value} or count < 2:
                continue
            rules.append(
                evaluate_rule(
                    sample,
                    stage_label,
                    "categorical_pair",
                    f"{left}__{right}__eq__{left_value}__{right_value}",
                    f"{left} == {left_value} and {right} == {right_value}",
                    lambda row, left=left, right=right, left_value=left_value, right_value=right_value: (
                        normalize_category(row.get(left)) == left_value
                        and normalize_category(row.get(right)) == right_value
                    ),
                    stage_rows,
                    fp_rows,
                    tp_rows,
                    truth_total,
                    final_f1,
                )
            )

    for feature, directions in NUMERIC_FEATURE_DIRECTIONS.items():
        values = [to_float(row.get(feature)) for row in stage_rows]
        thresholds = quantile_thresholds(values)
        for threshold in thresholds:
            for direction in directions:
                comparator = "<=" if direction == "le" else ">="
                rules.append(
                    evaluate_rule(
                        sample,
                        stage_label,
                        f"numeric_{direction}",
                        f"{feature}__{direction}__{threshold:.6f}",
                        f"{feature} {comparator} {threshold:.6f}",
                        lambda row, feature=feature, threshold=threshold, direction=direction: compare_numeric(
                            row.get(feature),
                            threshold,
                            direction,
                        ),
                        stage_rows,
                        fp_rows,
                        tp_rows,
                        truth_total,
                        final_f1,
                    )
                )

    rules.sort(
        key=lambda row: (
            float(row["delta_f1_vs_final"]),
            row["target_stage_removed"],
            -row["tp_removed"],
        ),
        reverse=True,
    )
    for rank, row in enumerate(rules, start=1):
        row["rule_rank"] = rank
    return rules


def compare_numeric(value: object, threshold: float, direction: str) -> bool:
    numeric = to_float(value)
    if math.isnan(numeric):
        return False
    if direction == "le":
        return numeric <= threshold
    if direction == "ge":
        return numeric >= threshold
    raise ValueError(f"Unsupported direction: {direction}")


def analyze_sample(
    config: Dict[str, object],
    provenance_workspace: Path,
) -> Dict[str, object]:
    sample = str(config["sample"])
    platform = str(config["platform"])
    master_path = provenance_workspace / f"{sample.lower()}_to_fp_provenance_master.tsv.gz"
    residual_master = read_residual_master_rows(master_path)

    to_tp_map, to_fp_map = load_to_feature_maps(Path(config["to_round_dir"]))
    _paired_tp_map, paired_fp_map = load_paired_feature_maps(Path(config["paired_dir"]))

    to_rule = config["metrics_to"]["rule"]
    to_post_removed_tp_keys = {key for key, row in to_tp_map.items() if rule_trigger(row, to_rule)}
    to_post_removed_fp_keys = {key for key, row in to_fp_map.items() if rule_trigger(row, to_rule)}

    final_tp_rows = [
        build_final_tp_row(sample, platform, key, feature)
        for key, feature in to_tp_map.items()
        if key not in to_post_removed_tp_keys
    ]
    final_fp_keys = {key for key in to_fp_map if key not in to_post_removed_fp_keys}
    master_keys = set(residual_master)
    if final_fp_keys != master_keys:
        missing_in_master = sorted(final_fp_keys - master_keys)[:5]
        missing_in_feature = sorted(master_keys - final_fp_keys)[:5]
        raise SystemExit(
            f"{sample}: final FP key mismatch between provenance master and feature map; "
            f"missing_in_master={missing_in_master}, missing_in_feature={missing_in_feature}"
        )

    residual_fp_rows = [
        build_residual_fp_row(sample, platform, residual_master[key], to_fp_map[key], paired_fp_map.get(key, {}))
        for key in sorted(final_fp_keys)
    ]

    final_metrics = config["metrics_to"]["filtered"]
    final_f1 = to_float(final_metrics["f1"])
    truth_total = int(config["truth_total"])

    stage_counter = Counter(row["stage_label"] for row in residual_fp_rows)
    stage_summary_rows: List[Dict[str, object]] = []
    for stage_label, count in sorted(stage_counter.items()):
        perfect_metrics = compute_metrics(len(final_tp_rows), len(residual_fp_rows) - count, truth_total)
        stage_summary_rows.append(
            {
                "sample": sample,
                "stage_label": stage_label,
                "count": count,
                "fraction_of_final_fp": format_float(count / len(residual_fp_rows)),
                "ideal_f1_if_perfect_removal": format_float(perfect_metrics["f1"]),
                "ideal_delta_f1_if_perfect_removal": format_float(perfect_metrics["f1"] - final_f1),
            }
        )

    feature_group_rows = [
        build_group_summary(sample, "final_tp", final_tp_rows),
        build_group_summary(sample, "raw_absent_all", [row for row in residual_fp_rows if row["stage_label"] == "raw_absent"]),
        build_group_summary(sample, "raw_filter_all", [row for row in residual_fp_rows if row["stage_label"] == "raw_filter"]),
        build_group_summary(sample, "longphase_s_all", [row for row in residual_fp_rows if row["stage_label"] == "longphase_s"]),
        build_group_summary(sample, "persistent_all", [row for row in residual_fp_rows if row["stage_label"] == "persistent"]),
    ]

    sample_summary = {
        "sample": sample,
        "platform": platform,
        "truth_total": truth_total,
        "to_final_tp": len(final_tp_rows),
        "to_final_fp": len(residual_fp_rows),
        "to_final_fn": truth_total - len(final_tp_rows),
        "to_final_precision": format_float(to_float(final_metrics["precision"])),
        "to_final_recall": format_float(to_float(final_metrics["recall"])),
        "to_final_f1": format_float(final_f1),
        "raw_absent_count": stage_counter.get("raw_absent", 0),
        "raw_filter_count": stage_counter.get("raw_filter", 0) + stage_counter.get("germline_filter", 0) + stage_counter.get("normal_alt_support", 0),
        "longphase_s_count": stage_counter.get("longphase_s", 0),
        "persistent_count": stage_counter.get("persistent", 0),
        "raw_absent_shared_exact_count": 0,
        "persistent_shared_exact_count": 0,
    }

    return {
        "sample": sample,
        "platform": platform,
        "truth_total": truth_total,
        "final_f1": final_f1,
        "final_tp_rows": final_tp_rows,
        "residual_fp_rows": residual_fp_rows,
        "sample_summary": sample_summary,
        "stage_summary_rows": stage_summary_rows,
        "feature_group_rows": feature_group_rows,
    }


def annotate_cross_platform(results: Dict[str, Dict[str, object]]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]]]:
    samples = list(results)
    if len(samples) < 2:
        return [], [], []

    recurrence_rows: List[Dict[str, object]] = []
    shared_hotspot_rows: List[Dict[str, object]] = []
    extra_feature_groups: List[Dict[str, object]] = []

    sample_a, sample_b = samples[:2]
    for stage_label in ("raw_absent", "persistent"):
        sample_a_map = {
            row["variant_key"]: row
            for row in results[sample_a]["residual_fp_rows"]
            if row["stage_label"] == stage_label
        }
        sample_b_map = {
            row["variant_key"]: row
            for row in results[sample_b]["residual_fp_rows"]
            if row["stage_label"] == stage_label
        }
        shared_keys = set(sample_a_map) & set(sample_b_map)
        same_pattern_keys = {
            key
            for key in shared_keys
            if (
                normalize_category(sample_a_map[key].get("verification_class"))
                == normalize_category(sample_b_map[key].get("verification_class"))
                and normalize_category(sample_a_map[key].get("agreement_type"))
                == normalize_category(sample_b_map[key].get("agreement_type"))
                and normalize_category(sample_a_map[key].get("class_shift"))
                == normalize_category(sample_b_map[key].get("class_shift"))
            )
        }

        for sample, row_map in ((sample_a, sample_a_map), (sample_b, sample_b_map)):
            for key, row in row_map.items():
                if key in shared_keys:
                    row["shared_exact_across_platform"] = True
                    row["shared_same_stage_with_other_platform"] = True
                    row["shared_same_pattern_with_other_platform"] = key in same_pattern_keys

        for sample, row_map in ((sample_a, sample_a_map), (sample_b, sample_b_map)):
            stage_total = len(row_map)
            recurrence_rows.append(
                {
                    "stage_label": stage_label,
                    "sample": sample,
                    "stage_total": stage_total,
                    "shared_exact_count": len(shared_keys),
                    "platform_specific_count": stage_total - len(shared_keys),
                    "fraction_shared_exact": format_float(len(shared_keys) / stage_total if stage_total else math.nan),
                    "same_pattern_shared_exact_count": len(same_pattern_keys),
                }
            )

        for key in sorted(shared_keys, key=lambda value: (sample_a_map[value]["chrom"], int(sample_a_map[value]["pos"]))):
            row_a = sample_a_map[key]
            row_b = sample_b_map[key]
            shared_hotspot_rows.append(
                {
                    "stage_label": stage_label,
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "variant_key": key,
                    "chrom": row_a["chrom"],
                    "pos": row_a["pos"],
                    "ref": row_a["ref"],
                    "alt": row_a["alt"],
                    "same_pattern": key in same_pattern_keys,
                    "sample_a_verification_class": row_a.get("verification_class", ""),
                    "sample_a_agreement_type": row_a.get("agreement_type", ""),
                    "sample_a_class_shift": row_a.get("class_shift", ""),
                    "sample_a_af": format_float(to_float(row_a.get("af"))),
                    "sample_a_gq": row_a.get("gq", ""),
                    "sample_a_qual": format_float(to_float(row_a.get("qual"))),
                    "sample_a_quality_score": format_float(to_float(row_a.get("quality_score"))),
                    "sample_a_pairwise_median_dist": format_float(to_float(row_a.get("pairwise_median_dist"))),
                    "sample_a_hp_assign_rate": format_float(to_float(row_a.get("hp_assign_rate"))),
                    "sample_b_verification_class": row_b.get("verification_class", ""),
                    "sample_b_agreement_type": row_b.get("agreement_type", ""),
                    "sample_b_class_shift": row_b.get("class_shift", ""),
                    "sample_b_af": format_float(to_float(row_b.get("af"))),
                    "sample_b_gq": row_b.get("gq", ""),
                    "sample_b_qual": format_float(to_float(row_b.get("qual"))),
                    "sample_b_quality_score": format_float(to_float(row_b.get("quality_score"))),
                    "sample_b_pairwise_median_dist": format_float(to_float(row_b.get("pairwise_median_dist"))),
                    "sample_b_hp_assign_rate": format_float(to_float(row_b.get("hp_assign_rate"))),
                }
            )

        for sample in (sample_a, sample_b):
            stage_rows = [row for row in results[sample]["residual_fp_rows"] if row["stage_label"] == stage_label]
            shared_rows = [row for row in stage_rows if as_bool(row.get("shared_exact_across_platform"))]
            specific_rows = [row for row in stage_rows if not as_bool(row.get("shared_exact_across_platform"))]
            extra_feature_groups.append(build_group_summary(sample, f"{stage_label}_shared_exact", shared_rows))
            extra_feature_groups.append(build_group_summary(sample, f"{stage_label}_platform_specific", specific_rows))

    for sample in samples:
        summary = results[sample]["sample_summary"]
        summary["raw_absent_shared_exact_count"] = sum(
            1
            for row in results[sample]["residual_fp_rows"]
            if row["stage_label"] == "raw_absent" and as_bool(row.get("shared_exact_across_platform"))
        )
        summary["persistent_shared_exact_count"] = sum(
            1
            for row in results[sample]["residual_fp_rows"]
            if row["stage_label"] == "persistent" and as_bool(row.get("shared_exact_across_platform"))
        )

    return recurrence_rows, shared_hotspot_rows, extra_feature_groups


def build_rule_scans(results: Dict[str, Dict[str, object]]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]]]:
    full_scan_rows: List[Dict[str, object]] = []
    best_rows: List[Dict[str, object]] = []
    pattern_rows: List[Dict[str, object]] = []

    for sample, result in results.items():
        tp_rows = result["final_tp_rows"]
        fp_rows = result["residual_fp_rows"]
        truth_total = result["truth_total"]
        final_f1 = result["final_f1"]
        stage_specs = [
            ("raw_absent", lambda row: row["stage_label"] == "raw_absent"),
            (
                "raw_absent_shared_exact",
                lambda row: row["stage_label"] == "raw_absent" and as_bool(row.get("shared_exact_across_platform")),
            ),
            (
                "raw_absent_platform_specific",
                lambda row: row["stage_label"] == "raw_absent" and not as_bool(row.get("shared_exact_across_platform")),
            ),
            ("persistent", lambda row: row["stage_label"] == "persistent"),
            (
                "persistent_shared_exact",
                lambda row: row["stage_label"] == "persistent" and as_bool(row.get("shared_exact_across_platform")),
            ),
            (
                "persistent_platform_specific",
                lambda row: row["stage_label"] == "persistent" and not as_bool(row.get("shared_exact_across_platform")),
            ),
        ]
        for stage_label, predicate in stage_specs:
            stage_rows = [row for row in fp_rows if predicate(row)]
            if not stage_rows:
                continue
            if stage_label in {"raw_absent", "persistent"}:
                pattern_rows.extend(build_pattern_summary(sample, stage_label, stage_rows, tp_rows))
            rules = build_stage_rules(sample, stage_label, stage_rows, fp_rows, tp_rows, truth_total, final_f1)
            full_scan_rows.extend(rules)
            best_rows.extend(rules[:20])
    return full_scan_rows, best_rows, pattern_rows


def build_persistent_candidate_rows(results: Dict[str, Dict[str, object]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for sample, result in results.items():
        persistent_rows = [row for row in result["residual_fp_rows"] if row["stage_label"] == "persistent"]
        for row in sorted(persistent_rows, key=lambda item: (not as_bool(item["shared_exact_across_platform"]), item["chrom"], int(item["pos"]))):
            rows.append(
                {
                    "sample": sample,
                    "shared_exact_across_platform": row["shared_exact_across_platform"],
                    "shared_same_pattern_with_other_platform": row["shared_same_pattern_with_other_platform"],
                    "variant_key": row["variant_key"],
                    "chrom": row["chrom"],
                    "pos": row["pos"],
                    "ref": row["ref"],
                    "alt": row["alt"],
                    "verification_class": row.get("verification_class", ""),
                    "agreement_type": row.get("agreement_type", ""),
                    "class_shift": row.get("class_shift", ""),
                    "af": format_float(to_float(row.get("af"))),
                    "gq": row.get("gq", ""),
                    "qual": format_float(to_float(row.get("qual"))),
                    "quality_score": format_float(to_float(row.get("quality_score"))),
                    "pairwise_median_dist": format_float(to_float(row.get("pairwise_median_dist"))),
                    "allele_delta": format_float(to_float(row.get("allele_delta"))),
                    "hp_assign_rate": format_float(to_float(row.get("hp_assign_rate"))),
                    "suggest_filter": row.get("suggest_filter", ""),
                    "potential_loh": row.get("potential_loh", ""),
                    "paired_raw_filter": row.get("paired_raw_filter", ""),
                    "paired_raw_af": row.get("paired_raw_af", ""),
                    "paired_raw_naf": row.get("paired_raw_naf", ""),
                    "paired_raw_gq": row.get("paired_raw_gq", ""),
                    "paired_raw_dp": row.get("paired_raw_dp", ""),
                    "paired_final_verification_class": row.get("paired_final_verification_class", ""),
                    "paired_final_agreement_type": row.get("paired_final_agreement_type", ""),
                    "paired_final_class_shift": row.get("paired_final_class_shift", ""),
                    "paired_final_af": format_float(to_float(row.get("paired_final_af"))),
                    "paired_final_gq": row.get("paired_final_gq", ""),
                    "paired_final_quality_score": format_float(to_float(row.get("paired_final_quality_score"))),
                    "paired_final_pairwise_median_dist": format_float(to_float(row.get("paired_final_pairwise_median_dist"))),
                    "paired_final_hp_assign_rate": format_float(to_float(row.get("paired_final_hp_assign_rate"))),
                }
            )
    return rows


def build_workspace_report(
    output_dir: Path,
    sample_summary_rows: Sequence[Dict[str, object]],
    stage_summary_rows: Sequence[Dict[str, object]],
    recurrence_rows: Sequence[Dict[str, object]],
    pattern_rows: Sequence[Dict[str, object]],
    best_rule_rows: Sequence[Dict[str, object]],
) -> None:
    top_patterns = [row for row in pattern_rows if row["stage_label"] == "raw_absent"][:8]
    top_persistent_rules = [row for row in best_rule_rows if row["stage_label"] == "persistent"][:6]
    top_raw_absent_rules = [row for row in best_rule_rows if row["stage_label"] == "raw_absent"][:6]

    lines = [
        "# TO residual FP deep dive workspace",
        "",
        "## Sample summary",
        "",
        markdown_table(
            [
                "sample",
                "to_final_tp",
                "to_final_fp",
                "to_final_f1",
                "raw_absent_count",
                "persistent_count",
                "raw_absent_shared_exact_count",
                "persistent_shared_exact_count",
            ],
            sample_summary_rows,
        ),
        "",
        "## Stage ceiling",
        "",
        markdown_table(
            [
                "sample",
                "stage_label",
                "count",
                "fraction_of_final_fp",
                "ideal_delta_f1_if_perfect_removal",
            ],
            stage_summary_rows,
        ),
        "",
        "## Cross-platform recurrence",
        "",
        markdown_table(
            [
                "stage_label",
                "sample",
                "stage_total",
                "shared_exact_count",
                "fraction_shared_exact",
                "same_pattern_shared_exact_count",
            ],
            recurrence_rows,
        ),
        "",
        "## Raw_absent top patterns",
        "",
        markdown_table(
            [
                "sample",
                "verification_class",
                "agreement_type",
                "class_shift",
                "count",
                "tp_match_count",
                "fp_to_tp_ratio",
                "shared_exact_count",
            ],
            top_patterns,
        ),
        "",
        "## Raw_absent best rules",
        "",
        markdown_table(
            [
                "sample",
                "scan_family",
                "rule_label",
                "target_stage_removed",
                "tp_removed",
                "delta_f1_vs_final",
            ],
            top_raw_absent_rules,
        ),
        "",
        "## Persistent best rules",
        "",
        markdown_table(
            [
                "sample",
                "scan_family",
                "rule_label",
                "target_stage_removed",
                "tp_removed",
                "delta_f1_vs_final",
            ],
            top_persistent_rules,
        ),
        "",
        "## Outputs",
        "",
        "- `sample_level_summary.tsv`",
        "- `residual_stage_summary.tsv`",
        "- `feature_group_summary.tsv`",
        "- `stage_pattern_summary.tsv`",
        "- `stage_rule_scan.tsv` / `stage_best_rules.tsv`",
        "- `cross_platform_recurrence_summary.tsv` / `cross_platform_shared_hotspots.tsv`",
        "- `paired_persistent_candidate_regions.tsv`",
    ]
    (output_dir / "analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    samples = args.sample or DEFAULT_SAMPLES
    for sample in samples:
        if sample not in SAMPLE_CONFIGS:
            raise SystemExit(f"Unsupported sample: {sample}")

    provenance_workspace = Path(args.provenance_workspace).resolve()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    sample_manifest_rows: List[Dict[str, object]] = []
    results: Dict[str, Dict[str, object]] = {}
    feature_group_rows: List[Dict[str, object]] = []
    stage_summary_rows: List[Dict[str, object]] = []

    for sample in samples:
        config = load_configs(sample)
        sample_manifest_rows.append(
            {
                "sample": sample,
                "platform": config["platform"],
                "truth_total": config["truth_total"],
                "to_round_dir": str(config["to_round_dir"]),
                "paired_dir": str(config["paired_dir"]),
                "provenance_master": str(provenance_workspace / f"{sample.lower()}_to_fp_provenance_master.tsv.gz"),
            }
        )
        result = analyze_sample(config, provenance_workspace)
        results[sample] = result
        feature_group_rows.extend(result["feature_group_rows"])
        stage_summary_rows.extend(result["stage_summary_rows"])

    recurrence_rows, shared_hotspot_rows, extra_feature_groups = annotate_cross_platform(results)
    feature_group_rows.extend(extra_feature_groups)

    sample_summary_rows = [results[sample]["sample_summary"] for sample in samples]
    full_scan_rows, best_rule_rows, pattern_rows = build_rule_scans(results)
    persistent_candidate_rows = build_persistent_candidate_rows(results)

    sample_manifest_fields = ["sample", "platform", "truth_total", "to_round_dir", "paired_dir", "provenance_master"]
    write_tsv_rows(output_dir / "sample_manifest.tsv", sample_manifest_fields, sample_manifest_rows)
    write_tsv_rows(output_dir / "sample_level_summary.tsv", SAMPLE_SUMMARY_FIELDS, sample_summary_rows)
    write_tsv_rows(output_dir / "residual_stage_summary.tsv", STAGE_SUMMARY_FIELDS, stage_summary_rows)
    write_tsv_rows(output_dir / "feature_group_summary.tsv", FEATURE_GROUP_FIELDS, feature_group_rows)
    write_tsv_rows(output_dir / "stage_pattern_summary.tsv", PATTERN_SUMMARY_FIELDS, pattern_rows)
    write_tsv_rows(output_dir / "stage_rule_scan.tsv", RULE_SCAN_FIELDS, full_scan_rows)
    write_tsv_rows(output_dir / "stage_best_rules.tsv", RULE_SCAN_FIELDS, best_rule_rows)
    write_tsv_rows(output_dir / "cross_platform_recurrence_summary.tsv", RECURRENCE_SUMMARY_FIELDS, recurrence_rows)
    write_tsv_rows(output_dir / "cross_platform_shared_hotspots.tsv", SHARED_HOTSPOT_FIELDS, shared_hotspot_rows)
    write_tsv_rows(output_dir / "paired_persistent_candidate_regions.tsv", PERSISTENT_CANDIDATE_FIELDS, persistent_candidate_rows)

    for sample in samples:
        lower = sample.lower()
        write_gzip_tsv(output_dir / f"{lower}_residual_fp_enriched.tsv.gz", RESIDUAL_FP_FIELDS, results[sample]["residual_fp_rows"])
        write_gzip_tsv(output_dir / f"{lower}_final_tp_enriched.tsv.gz", FINAL_TP_FIELDS, results[sample]["final_tp_rows"])

    build_workspace_report(
        output_dir,
        sample_summary_rows,
        stage_summary_rows,
        recurrence_rows,
        pattern_rows,
        best_rule_rows,
    )
    print(f"[build_to_residual_fp_deep_dive] Wrote workspace to {output_dir}")


if __name__ == "__main__":
    main()
