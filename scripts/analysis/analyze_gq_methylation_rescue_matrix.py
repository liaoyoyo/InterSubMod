#!/usr/bin/env python3
"""Systematic GQ / methylation rescue analysis across candidate-specific datasets."""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence

import pandas as pd

from research_common import compute_metrics, ensure_dir, write_tsv_rows


@dataclass
class DatasetConfig:
    dataset_id: str
    label: str
    raw_pool_tsv: Path
    joined_tsv: Path
    baseline_tp: int
    baseline_fp: int
    truth_total: int


RuleMask = Callable[[pd.DataFrame], pd.Series]


POOL_ELIGIBILITY_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "caller",
    "mode",
    "downstream_status",
    "rows_total",
    "candidate_eligible_rows",
    "candidate_eligible_fraction",
    "excluded_rows",
    "excluded_gq_ge_10",
    "excluded_gq_ge_15",
    "excluded_gq_ge_20",
    "eligible_gq_ge_10",
    "eligible_gq_ge_15",
    "eligible_gq_ge_20",
]

EXCLUSION_REASON_FIELDS = [
    "dataset_id",
    "label",
    "downstream_status",
    "rescue_exclusion_reason",
    "count",
    "fraction",
]

DATASET_OVERVIEW_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "caller",
    "mode",
    "baseline_tp",
    "baseline_fp",
    "baseline_fn",
    "baseline_f1",
    "raw_lost_tp_total",
    "raw_removed_fp_total",
    "joined_candidate_rows",
    "joined_lost_tp_rows",
    "joined_removed_fp_rows",
    "analyzed_rows",
    "analyzed_lost_tp_rows",
    "analyzed_removed_fp_rows",
    "lost_tp_coverage",
    "removed_fp_coverage",
]

GQ_SWEEP_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "caller",
    "mode",
    "threshold_type",
    "threshold_value",
    "rule_id",
    "trigger_count",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "fp_per_tp",
    "meets_safety",
]

METHYL_ONLY_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "caller",
    "mode",
    "feature_family",
    "rule_id",
    "rule_notes",
    "trigger_count",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "fp_per_tp",
    "meets_safety",
]

COMBO_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "caller",
    "mode",
    "gate_id",
    "gate_notes",
    "rule_scope",
    "feature_family",
    "rule_id",
    "rule_notes",
    "trigger_count",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "delta_f1_vs_gate",
    "fp_per_tp",
    "meets_safety",
]

ARTIFACT_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "caller",
    "mode",
    "artifact_rule",
    "role",
    "base_gate",
    "tp_selected",
    "fp_selected",
    "tp_removed_by_veto",
    "fp_removed_by_veto",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "delta_f1_vs_gate",
]

INTERVAL_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "feature",
    "subset_scope",
    "bin_label",
    "tp_count",
    "fp_count",
    "tp_fraction",
    "fp_fraction",
    "tp_to_fp_ratio",
    "enrichment_ratio",
]

GRID_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "grid_id",
    "x_bin",
    "y_bin",
    "tp_count",
    "fp_count",
    "tp_fraction",
    "fp_fraction",
    "tp_to_fp_ratio",
]

ORTHOGONAL_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "feature_a",
    "feature_b",
    "jaccard",
    "a_only_tp",
    "a_only_fp",
    "b_only_tp",
    "b_only_fp",
    "both_tp",
    "both_fp",
    "union_tp",
    "union_fp",
]

SUMMARY_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "summary_type",
    "rule_id",
    "tp_rescued",
    "fp_reintroduced",
    "delta_f1_vs_baseline",
    "delta_f1_vs_gate",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        required=True,
        help="Format: dataset_id|label|raw_candidate_pool_tsv|joined_tsv|baseline_tp|baseline_fp|truth_total",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def parse_dataset(value: str) -> DatasetConfig:
    parts = value.split("|")
    if len(parts) != 7:
        raise ValueError(f"Invalid --dataset value: {value}")
    dataset_id, label, raw_pool_tsv, joined_tsv, baseline_tp, baseline_fp, truth_total = parts
    return DatasetConfig(
        dataset_id=dataset_id,
        label=label,
        raw_pool_tsv=Path(raw_pool_tsv).resolve(),
        joined_tsv=Path(joined_tsv).resolve(),
        baseline_tp=int(baseline_tp),
        baseline_fp=int(baseline_fp),
        truth_total=int(truth_total),
    )


def as_bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).map(lambda value: str(value).strip().lower() in {"1", "true", "yes", "y"})


def load_raw_pool(path: Path) -> pd.DataFrame:
    usecols = [
        "sample",
        "platform",
        "caller",
        "mode",
        "region_key",
        "truth_status",
        "downstream_status",
        "filter",
        "qual",
        "gq",
        "dp",
        "af",
        "ad_ref",
        "ad_alt",
        "candidate_eligible",
        "rescue_exclusion_reason",
    ]
    df = pd.read_csv(path, sep="\t", usecols=usecols)
    for column in ["qual", "gq", "dp", "af", "ad_ref", "ad_alt"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["candidate_eligible"] = as_bool_series(df["candidate_eligible"])
    df["rescue_exclusion_reason"] = df["rescue_exclusion_reason"].fillna("NA").replace("", "NA")
    df["truth_status"] = df["truth_status"].fillna("NA")
    df["downstream_status"] = df["downstream_status"].fillna("NA")
    return df


def load_joined(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    numeric_columns = [
        "qual",
        "gq",
        "dp",
        "af",
        "ad_ref",
        "ad_alt",
        "PairwiseMeanDist",
        "PairwiseMedianDist",
        "AlleleDelta",
        "CramersV",
        "GlobalP",
        "Quality_Score",
        "hp_assign_rate",
        "allele_assign_rate",
    ]
    bool_columns = [
        "candidate_eligible",
        "rescue_filter_eligible",
        "meets_numeric_thresholds",
        "has_h_flag",
        "verdict_somatic",
        "verdict_subclonal",
        "verdict_germline",
        "pon_hit",
        "multihap_flag",
        "noancestry_flag",
        "PassedGating",
    ]
    text_columns = {
        "VerificationClass": "NA",
        "agreement_type": "NA",
        "cluster_class": "NA",
        "label_class": "NA",
        "DominantLabel": "NA",
        "class_shift": "NA",
        "source_scope": "",
        "truth_status": "NA",
        "downstream_status": "NA",
        "filter": "NA",
    }

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    for column in bool_columns:
        if column in df.columns:
            df[column] = as_bool_series(df[column])
        else:
            df[column] = False
    for column, default in text_columns.items():
        if column in df.columns:
            df[column] = df[column].fillna(default).replace("", default)
        else:
            df[column] = default

    df["analysis_available"] = df["source_scope"].isin({"tp", "fp"})
    df["strong_subclone"] = df["VerificationClass"].isin({"Strong", "Subclone"})
    df["agreement_positive"] = df["agreement_type"].isin({"label_upgrade", "consistent_strong", "consistent_subclone"})
    return df


def safe_ratio(num: int, den: int) -> str:
    if den == 0:
        return "0.000000"
    return f"{num / den:.6f}"


def format_float(value: float) -> str:
    if math.isinf(value):
        return "inf"
    if math.isnan(value):
        return "nan"
    return f"{value:.6f}"


def evaluate_mask(df: pd.DataFrame, mask: pd.Series, baseline_tp: int, baseline_fp: int, truth_total: int) -> Dict[str, object]:
    selected = df[mask].copy()
    tp_rescued = int((selected["downstream_status"] == "caller_lost_tp").sum())
    fp_reintroduced = int((selected["downstream_status"] == "caller_removed_fp").sum())
    metrics = compute_metrics(baseline_tp + tp_rescued, baseline_fp + fp_reintroduced, truth_total)
    baseline = compute_metrics(baseline_tp, baseline_fp, truth_total)
    fp_per_tp = fp_reintroduced / tp_rescued if tp_rescued else float("inf")
    return {
        "trigger_count": int(mask.sum()),
        "tp_rescued": tp_rescued,
        "fp_reintroduced": fp_reintroduced,
        "precision": format_float(metrics["precision"]),
        "recall": format_float(metrics["recall"]),
        "f1": format_float(metrics["f1"]),
        "delta_f1_vs_baseline": format_float(metrics["f1"] - baseline["f1"]),
        "fp_per_tp": format_float(fp_per_tp),
        "meets_safety": tp_rescued > 0 and fp_reintroduced <= tp_rescued,
    }


def raw_status_counts(raw_df: pd.DataFrame, status: str) -> pd.DataFrame:
    return raw_df[raw_df["downstream_status"] == status].copy()


def candidate_eval_frame(joined_df: pd.DataFrame) -> pd.DataFrame:
    keep_status = {"caller_lost_tp", "caller_removed_fp"}
    mask = joined_df["candidate_eligible"] & joined_df["downstream_status"].isin(keep_status)
    return joined_df[mask].copy()


def gate_definitions() -> List[Dict[str, object]]:
    thresholds = [5, 8, 10, 12, 15, 18, 20, 25]
    rules: List[Dict[str, object]] = []
    for threshold in thresholds:
        rules.append(
            {
                "threshold_type": "gq_ge",
                "threshold_value": threshold,
                "rule_id": f"gq_ge_{threshold}",
                "func": lambda df, t=threshold: df["gq"] >= t,
            }
        )
    rules.append(
        {
            "threshold_type": "qual_or_gq",
            "threshold_value": 10,
            "rule_id": "qual_ge_10_or_gq_ge_20",
            "func": lambda df: (df["qual"] >= 10.0) | (df["gq"] >= 20),
        }
    )
    return rules


def support_rule_defs() -> List[Dict[str, object]]:
    return [
        {"family": "pairwise_median", "rule_id": "pairwise_ge_015", "notes": "PairwiseMedianDist >= 0.15", "func": lambda df: df["PairwiseMedianDist"] >= 0.15},
        {"family": "pairwise_median", "rule_id": "pairwise_ge_020", "notes": "PairwiseMedianDist >= 0.20", "func": lambda df: df["PairwiseMedianDist"] >= 0.20},
        {"family": "pairwise_median", "rule_id": "pairwise_ge_025", "notes": "PairwiseMedianDist >= 0.25", "func": lambda df: df["PairwiseMedianDist"] >= 0.25},
        {"family": "pairwise_median", "rule_id": "pairwise_le_020", "notes": "PairwiseMedianDist <= 0.20", "func": lambda df: df["PairwiseMedianDist"] <= 0.20},
        {"family": "pairwise_median", "rule_id": "pairwise_le_015", "notes": "PairwiseMedianDist <= 0.15", "func": lambda df: df["PairwiseMedianDist"] <= 0.15},
        {"family": "pairwise_mean", "rule_id": "pairwise_mean_ge_020", "notes": "PairwiseMeanDist >= 0.20", "func": lambda df: df["PairwiseMeanDist"] >= 0.20},
        {"family": "quality_score", "rule_id": "quality_ge_40", "notes": "Quality_Score >= 40", "func": lambda df: df["Quality_Score"] >= 40},
        {"family": "quality_score", "rule_id": "quality_ge_50", "notes": "Quality_Score >= 50", "func": lambda df: df["Quality_Score"] >= 50},
        {"family": "quality_score", "rule_id": "quality_ge_60", "notes": "Quality_Score >= 60", "func": lambda df: df["Quality_Score"] >= 60},
        {"family": "quality_score", "rule_id": "quality_ge_70", "notes": "Quality_Score >= 70", "func": lambda df: df["Quality_Score"] >= 70},
        {"family": "quality_score", "rule_id": "quality_ge_80", "notes": "Quality_Score >= 80", "func": lambda df: df["Quality_Score"] >= 80},
        {"family": "agreement", "rule_id": "agreement_positive", "notes": "agreement_type in label_upgrade/consistent_strong/consistent_subclone", "func": lambda df: df["agreement_positive"]},
        {"family": "verification", "rule_id": "strong_subclone", "notes": "VerificationClass in Strong/Subclone", "func": lambda df: df["strong_subclone"]},
        {"family": "agreement", "rule_id": "label_upgrade", "notes": "agreement_type == label_upgrade", "func": lambda df: df["agreement_type"] == "label_upgrade"},
        {"family": "agreement", "rule_id": "consistent_strong", "notes": "agreement_type == consistent_strong", "func": lambda df: df["agreement_type"] == "consistent_strong"},
        {"family": "phase", "rule_id": "hp_assign_ge_090", "notes": "hp_assign_rate >= 0.90", "func": lambda df: df["hp_assign_rate"] >= 0.90},
        {"family": "phase", "rule_id": "hp_assign_ge_099", "notes": "hp_assign_rate >= 0.99", "func": lambda df: df["hp_assign_rate"] >= 0.99},
        {"family": "phase", "rule_id": "allele_assign_ge_099", "notes": "allele_assign_rate >= 0.99", "func": lambda df: df["allele_assign_rate"] >= 0.99},
        {"family": "global_p", "rule_id": "globalp_le_050", "notes": "GlobalP <= 0.50", "func": lambda df: df["GlobalP"] <= 0.50},
        {"family": "global_p", "rule_id": "globalp_le_005", "notes": "GlobalP <= 0.05", "func": lambda df: df["GlobalP"] <= 0.05},
    ]


def artifact_low_vaf_high_adelta(df: pd.DataFrame) -> pd.Series:
    return (df["af"] < 0.24) & (df["AlleleDelta"] > 0.25)


def artifact_low_vaf_high_adelta_low_cv(df: pd.DataFrame) -> pd.Series:
    return artifact_low_vaf_high_adelta(df) & (df["CramersV"] < 0.05)


def artifact_combined(df: pd.DataFrame) -> pd.Series:
    return (
        artifact_low_vaf_high_adelta_low_cv(df)
        | ((df["cluster_class"] == "Strong") & (df["label_class"] == "Weak"))
        | ((df["PairwiseMedianDist"] < 0.12) & (df["AlleleDelta"] > 0.15))
        | ((df["hp_assign_rate"] < 0.50) & (df["cluster_class"] == "Strong"))
    )


def artifact_defs() -> List[Dict[str, object]]:
    return [
        {
            "rule_id": "lowvaf_highadelta",
            "notes": "AF < 0.24 and AlleleDelta > 0.25",
            "func": artifact_low_vaf_high_adelta,
        },
        {
            "rule_id": "lowvaf_highadelta_lowcv",
            "notes": "AF < 0.24 and AlleleDelta > 0.25 and CramersV < 0.05",
            "func": artifact_low_vaf_high_adelta_low_cv,
        },
        {
            "rule_id": "combined_artifact_veto",
            "notes": "lowVAF/highAlleleDelta/lowCV or cluster Strong label Weak or low Pairwise/high AlleleDelta or low HP assign Strong",
            "func": artifact_combined,
        },
    ]


def gq_combo_gates() -> List[Dict[str, object]]:
    return [
        {"gate_id": "gq_ge_10", "notes": "candidate_eligible & gq>=10", "func": lambda df: df["gq"] >= 10},
        {"gate_id": "gq_ge_15", "notes": "candidate_eligible & gq>=15", "func": lambda df: df["gq"] >= 15},
        {"gate_id": "gq_ge_20", "notes": "candidate_eligible & gq>=20", "func": lambda df: df["gq"] >= 20},
        {
            "gate_id": "qual_ge_10_or_gq_ge_20",
            "notes": "candidate_eligible & (qual>=10 or gq>=20)",
            "func": lambda df: (df["qual"] >= 10.0) | (df["gq"] >= 20),
        },
    ]


def selected_combo_rule_ids() -> Sequence[str]:
    return [
        "pairwise_ge_015",
        "pairwise_ge_020",
        "pairwise_le_020",
        "quality_ge_60",
        "quality_ge_70",
        "agreement_positive",
        "strong_subclone",
        "hp_assign_ge_099",
        "globalp_le_050",
    ]


def feature_bins() -> Dict[str, Sequence[float]]:
    return {
        "gq": [0, 5, 10, 15, 20, 25, math.inf],
        "PairwiseMedianDist": [0.0, 0.10, 0.15, 0.20, 0.25, math.inf],
        "Quality_Score": [0.0, 50.0, 60.0, 70.0, 80.0, math.inf],
        "af": [0.0, 0.05, 0.10, 0.20, 0.30, 0.50, 1.01],
        "AlleleDelta": [0.0, 0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 1.01],
        "CramersV": [0.0, 0.02, 0.05, 0.10, 0.15, 1.01],
        "hp_assign_rate": [0.0, 0.50, 0.70, 0.85, 0.95, 0.99, 1.01],
    }


def bin_labels(edges: Sequence[float]) -> List[str]:
    labels: List[str] = []
    for left, right in zip(edges[:-1], edges[1:]):
        right_text = "inf" if math.isinf(right) else f"{right:g}"
        labels.append(f"[{left:g},{right_text})")
    return labels


def add_pool_eligibility_rows(config: DatasetConfig, raw_df: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    meta = raw_df.iloc[0].to_dict() if not raw_df.empty else {}
    for status in ["caller_lost_tp", "caller_removed_fp"]:
        sub = raw_status_counts(raw_df, status)
        eligible = sub[sub["candidate_eligible"]]
        excluded = sub[~sub["candidate_eligible"]]
        rows.append(
            {
                "dataset_id": config.dataset_id,
                "label": config.label,
                "sample": meta.get("sample", "unknown"),
                "platform": meta.get("platform", "unknown"),
                "caller": meta.get("caller", "unknown"),
                "mode": meta.get("mode", "unknown"),
                "downstream_status": status,
                "rows_total": len(sub),
                "candidate_eligible_rows": len(eligible),
                "candidate_eligible_fraction": safe_ratio(len(eligible), len(sub)),
                "excluded_rows": len(excluded),
                "excluded_gq_ge_10": int((excluded["gq"] >= 10).sum()),
                "excluded_gq_ge_15": int((excluded["gq"] >= 15).sum()),
                "excluded_gq_ge_20": int((excluded["gq"] >= 20).sum()),
                "eligible_gq_ge_10": int((eligible["gq"] >= 10).sum()),
                "eligible_gq_ge_15": int((eligible["gq"] >= 15).sum()),
                "eligible_gq_ge_20": int((eligible["gq"] >= 20).sum()),
            }
        )
    return rows


def add_exclusion_reason_rows(config: DatasetConfig, raw_df: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for status in ["caller_lost_tp", "caller_removed_fp"]:
        excluded = raw_status_counts(raw_df, status)
        excluded = excluded[~excluded["candidate_eligible"]]
        if excluded.empty:
            continue
        counts = excluded["rescue_exclusion_reason"].value_counts()
        total = int(counts.sum())
        for reason, count in counts.items():
            rows.append(
                {
                    "dataset_id": config.dataset_id,
                    "label": config.label,
                    "downstream_status": status,
                    "rescue_exclusion_reason": reason,
                    "count": int(count),
                    "fraction": safe_ratio(int(count), total),
                }
            )
    return rows


def dataset_overview_row(config: DatasetConfig, raw_df: pd.DataFrame, eval_df: pd.DataFrame) -> Dict[str, object]:
    baseline = compute_metrics(config.baseline_tp, config.baseline_fp, config.truth_total)
    meta = raw_df.iloc[0].to_dict() if not raw_df.empty else {}
    raw_lost = raw_status_counts(raw_df, "caller_lost_tp")
    raw_fp = raw_status_counts(raw_df, "caller_removed_fp")
    analyzed = eval_df[eval_df["analysis_available"]]
    analyzed_lost = analyzed[analyzed["downstream_status"] == "caller_lost_tp"]
    analyzed_fp = analyzed[analyzed["downstream_status"] == "caller_removed_fp"]
    joined_lost = eval_df[eval_df["downstream_status"] == "caller_lost_tp"]
    joined_fp = eval_df[eval_df["downstream_status"] == "caller_removed_fp"]
    return {
        "dataset_id": config.dataset_id,
        "label": config.label,
        "sample": meta.get("sample", "unknown"),
        "platform": meta.get("platform", "unknown"),
        "caller": meta.get("caller", "unknown"),
        "mode": meta.get("mode", "unknown"),
        "baseline_tp": config.baseline_tp,
        "baseline_fp": config.baseline_fp,
        "baseline_fn": baseline["fn"],
        "baseline_f1": format_float(baseline["f1"]),
        "raw_lost_tp_total": len(raw_lost),
        "raw_removed_fp_total": len(raw_fp),
        "joined_candidate_rows": len(eval_df),
        "joined_lost_tp_rows": len(joined_lost),
        "joined_removed_fp_rows": len(joined_fp),
        "analyzed_rows": len(analyzed),
        "analyzed_lost_tp_rows": len(analyzed_lost),
        "analyzed_removed_fp_rows": len(analyzed_fp),
        "lost_tp_coverage": safe_ratio(len(analyzed_lost), len(joined_lost)),
        "removed_fp_coverage": safe_ratio(len(analyzed_fp), len(joined_fp)),
    }


def add_gq_sweep_rows(config: DatasetConfig, eval_df: pd.DataFrame) -> List[Dict[str, object]]:
    meta = eval_df.iloc[0].to_dict() if not eval_df.empty else {}
    rows: List[Dict[str, object]] = []
    for gate in gate_definitions():
        metrics = evaluate_mask(eval_df, gate["func"](eval_df), config.baseline_tp, config.baseline_fp, config.truth_total)
        rows.append(
            {
                "dataset_id": config.dataset_id,
                "label": config.label,
                "sample": meta.get("sample", "unknown"),
                "platform": meta.get("platform", "unknown"),
                "caller": meta.get("caller", "unknown"),
                "mode": meta.get("mode", "unknown"),
                "threshold_type": gate["threshold_type"],
                "threshold_value": gate["threshold_value"],
                "rule_id": gate["rule_id"],
                **metrics,
            }
        )
    return rows


def add_methyl_only_rows(config: DatasetConfig, eval_df: pd.DataFrame) -> List[Dict[str, object]]:
    meta = eval_df.iloc[0].to_dict() if not eval_df.empty else {}
    rows: List[Dict[str, object]] = []
    analyzed = eval_df[eval_df["analysis_available"]].copy()
    for rule in support_rule_defs():
        mask = rule["func"](analyzed)
        metrics = evaluate_mask(analyzed, mask, config.baseline_tp, config.baseline_fp, config.truth_total)
        rows.append(
            {
                "dataset_id": config.dataset_id,
                "label": config.label,
                "sample": meta.get("sample", "unknown"),
                "platform": meta.get("platform", "unknown"),
                "caller": meta.get("caller", "unknown"),
                "mode": meta.get("mode", "unknown"),
                "feature_family": rule["family"],
                "rule_id": rule["rule_id"],
                "rule_notes": rule["notes"],
                **metrics,
            }
        )
    return rows


def gate_baselines(config: DatasetConfig, eval_df: pd.DataFrame) -> Dict[str, Dict[str, object]]:
    baselines: Dict[str, Dict[str, object]] = {}
    for gate in gq_combo_gates():
        metrics = evaluate_mask(eval_df, gate["func"](eval_df), config.baseline_tp, config.baseline_fp, config.truth_total)
        baselines[gate["gate_id"]] = metrics
    return baselines


def add_combo_rows(config: DatasetConfig, eval_df: pd.DataFrame) -> List[Dict[str, object]]:
    meta = eval_df.iloc[0].to_dict() if not eval_df.empty else {}
    rows: List[Dict[str, object]] = []
    analyzed = eval_df[eval_df["analysis_available"]].copy()
    support_rules = [rule for rule in support_rule_defs() if rule["rule_id"] in set(selected_combo_rule_ids())]
    gate_metrics = gate_baselines(config, eval_df)

    for gate in gq_combo_gates():
        gate_mask = gate["func"](analyzed)
        for rule in support_rules:
            mask = gate_mask & rule["func"](analyzed)
            metrics = evaluate_mask(analyzed, mask, config.baseline_tp, config.baseline_fp, config.truth_total)
            rows.append(
                {
                    "dataset_id": config.dataset_id,
                    "label": config.label,
                    "sample": meta.get("sample", "unknown"),
                    "platform": meta.get("platform", "unknown"),
                    "caller": meta.get("caller", "unknown"),
                    "mode": meta.get("mode", "unknown"),
                    "gate_id": gate["gate_id"],
                    "gate_notes": gate["notes"],
                    "rule_scope": "single_support",
                    "feature_family": rule["family"],
                    "rule_id": f"{gate['gate_id']}__{rule['rule_id']}",
                    "rule_notes": rule["notes"],
                    **metrics,
                    "delta_f1_vs_gate": format_float(float(metrics["f1"]) - float(gate_metrics[gate["gate_id"]]["f1"])),
                }
            )

        for rule_a, rule_b in itertools.combinations(support_rules, 2):
            mask = gate_mask & rule_a["func"](analyzed) & rule_b["func"](analyzed)
            metrics = evaluate_mask(analyzed, mask, config.baseline_tp, config.baseline_fp, config.truth_total)
            rows.append(
                {
                    "dataset_id": config.dataset_id,
                    "label": config.label,
                    "sample": meta.get("sample", "unknown"),
                    "platform": meta.get("platform", "unknown"),
                    "caller": meta.get("caller", "unknown"),
                    "mode": meta.get("mode", "unknown"),
                    "gate_id": gate["gate_id"],
                    "gate_notes": gate["notes"],
                    "rule_scope": "pair_support",
                    "feature_family": f"{rule_a['family']}+{rule_b['family']}",
                    "rule_id": f"{gate['gate_id']}__{rule_a['rule_id']}__{rule_b['rule_id']}",
                    "rule_notes": f"{rule_a['notes']} && {rule_b['notes']}",
                    **metrics,
                    "delta_f1_vs_gate": format_float(float(metrics["f1"]) - float(gate_metrics[gate["gate_id"]]["f1"])),
                }
            )
    return rows


def add_artifact_rows(config: DatasetConfig, eval_df: pd.DataFrame) -> List[Dict[str, object]]:
    meta = eval_df.iloc[0].to_dict() if not eval_df.empty else {}
    analyzed = eval_df[eval_df["analysis_available"]].copy()
    gate_metrics = gate_baselines(config, eval_df)
    rows: List[Dict[str, object]] = []

    for artifact in artifact_defs():
        art_mask = artifact["func"](analyzed)
        metrics = evaluate_mask(analyzed, art_mask, config.baseline_tp, config.baseline_fp, config.truth_total)
        rows.append(
            {
                "dataset_id": config.dataset_id,
                "label": config.label,
                "sample": meta.get("sample", "unknown"),
                "platform": meta.get("platform", "unknown"),
                "caller": meta.get("caller", "unknown"),
                "mode": meta.get("mode", "unknown"),
                "artifact_rule": artifact["rule_id"],
                "role": "early_selection",
                "base_gate": "none",
                "tp_selected": metrics["tp_rescued"],
                "fp_selected": metrics["fp_reintroduced"],
                "tp_removed_by_veto": 0,
                "fp_removed_by_veto": 0,
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1": metrics["f1"],
                "delta_f1_vs_baseline": metrics["delta_f1_vs_baseline"],
                "delta_f1_vs_gate": "nan",
            }
        )
        for gate in gq_combo_gates():
            gate_mask = gate["func"](analyzed)
            selected_mask = gate_mask & (~art_mask)
            selected_metrics = evaluate_mask(analyzed, selected_mask, config.baseline_tp, config.baseline_fp, config.truth_total)
            veto_mask = gate_mask & art_mask
            tp_removed = int(((analyzed["downstream_status"] == "caller_lost_tp") & veto_mask).sum())
            fp_removed = int(((analyzed["downstream_status"] == "caller_removed_fp") & veto_mask).sum())
            rows.append(
                {
                    "dataset_id": config.dataset_id,
                    "label": config.label,
                    "sample": meta.get("sample", "unknown"),
                    "platform": meta.get("platform", "unknown"),
                    "caller": meta.get("caller", "unknown"),
                    "mode": meta.get("mode", "unknown"),
                    "artifact_rule": artifact["rule_id"],
                    "role": "rescue_veto",
                    "base_gate": gate["gate_id"],
                    "tp_selected": selected_metrics["tp_rescued"],
                    "fp_selected": selected_metrics["fp_reintroduced"],
                    "tp_removed_by_veto": tp_removed,
                    "fp_removed_by_veto": fp_removed,
                    "precision": selected_metrics["precision"],
                    "recall": selected_metrics["recall"],
                    "f1": selected_metrics["f1"],
                    "delta_f1_vs_baseline": selected_metrics["delta_f1_vs_baseline"],
                    "delta_f1_vs_gate": format_float(float(selected_metrics["f1"]) - float(gate_metrics[gate["gate_id"]]["f1"])),
                }
            )
    return rows


def cut_series(series: pd.Series, edges: Sequence[float]) -> pd.Series:
    return pd.cut(series, bins=edges, right=False, labels=bin_labels(edges), include_lowest=True)


def add_interval_rows(config: DatasetConfig, eval_df: pd.DataFrame) -> List[Dict[str, object]]:
    meta = eval_df.iloc[0].to_dict() if not eval_df.empty else {}
    rows: List[Dict[str, object]] = []
    bins = feature_bins()

    for feature, edges in bins.items():
        if feature not in eval_df.columns:
            continue
        subset_scope = "analyzed" if feature != "gq" else "candidate"
        base = eval_df if feature == "gq" else eval_df[eval_df["analysis_available"]].copy()
        base = base[pd.notna(base[feature])].copy()
        if base.empty:
            continue
        base["bin_label"] = cut_series(base[feature], edges)
        total_tp = int((base["downstream_status"] == "caller_lost_tp").sum())
        total_fp = int((base["downstream_status"] == "caller_removed_fp").sum())
        for bin_label, part in base.groupby("bin_label", observed=False):
            tp_count = int((part["downstream_status"] == "caller_lost_tp").sum())
            fp_count = int((part["downstream_status"] == "caller_removed_fp").sum())
            tp_fraction = tp_count / total_tp if total_tp else 0.0
            fp_fraction = fp_count / total_fp if total_fp else 0.0
            ratio = tp_count / fp_count if fp_count else (float("inf") if tp_count else 0.0)
            enrichment = tp_fraction / fp_fraction if fp_fraction else (float("inf") if tp_fraction else 0.0)
            rows.append(
                {
                    "dataset_id": config.dataset_id,
                    "label": config.label,
                    "sample": meta.get("sample", "unknown"),
                    "mode": meta.get("mode", "unknown"),
                    "feature": feature,
                    "subset_scope": subset_scope,
                    "bin_label": str(bin_label),
                    "tp_count": tp_count,
                    "fp_count": fp_count,
                    "tp_fraction": format_float(tp_fraction),
                    "fp_fraction": format_float(fp_fraction),
                    "tp_to_fp_ratio": format_float(ratio),
                    "enrichment_ratio": format_float(enrichment),
                }
            )
    return rows


def add_grid_rows(config: DatasetConfig, eval_df: pd.DataFrame) -> List[Dict[str, object]]:
    meta = eval_df.iloc[0].to_dict() if not eval_df.empty else {}
    rows: List[Dict[str, object]] = []
    analyzed = eval_df[eval_df["analysis_available"]].copy()
    grids = [
        ("gq_vs_pairwise", "gq", "PairwiseMedianDist"),
        ("gq_vs_quality", "gq", "Quality_Score"),
        ("af_vs_allele_delta", "af", "AlleleDelta"),
        ("pairwise_vs_hp_assign", "PairwiseMedianDist", "hp_assign_rate"),
    ]
    bins = feature_bins()
    for grid_id, x_feature, y_feature in grids:
        if x_feature not in analyzed.columns or y_feature not in analyzed.columns:
            continue
        base = analyzed[pd.notna(analyzed[x_feature]) & pd.notna(analyzed[y_feature])].copy()
        if base.empty:
            continue
        base["x_bin"] = cut_series(base[x_feature], bins[x_feature])
        base["y_bin"] = cut_series(base[y_feature], bins[y_feature])
        total_tp = int((base["downstream_status"] == "caller_lost_tp").sum())
        total_fp = int((base["downstream_status"] == "caller_removed_fp").sum())
        grouped = base.groupby(["x_bin", "y_bin"], observed=False)
        for (x_bin, y_bin), part in grouped:
            tp_count = int((part["downstream_status"] == "caller_lost_tp").sum())
            fp_count = int((part["downstream_status"] == "caller_removed_fp").sum())
            ratio = tp_count / fp_count if fp_count else (float("inf") if tp_count else 0.0)
            rows.append(
                {
                    "dataset_id": config.dataset_id,
                    "label": config.label,
                    "sample": meta.get("sample", "unknown"),
                    "mode": meta.get("mode", "unknown"),
                    "grid_id": grid_id,
                    "x_bin": str(x_bin),
                    "y_bin": str(y_bin),
                    "tp_count": tp_count,
                    "fp_count": fp_count,
                    "tp_fraction": format_float(tp_count / total_tp if total_tp else 0.0),
                    "fp_fraction": format_float(fp_count / total_fp if total_fp else 0.0),
                    "tp_to_fp_ratio": format_float(ratio),
                }
            )
    return rows


def orthogonal_feature_masks(df: pd.DataFrame) -> Dict[str, pd.Series]:
    return {
        "pairwise_ge_020": df["PairwiseMedianDist"] >= 0.20,
        "pairwise_le_020": df["PairwiseMedianDist"] <= 0.20,
        "quality_ge_60": df["Quality_Score"] >= 60,
        "agreement_positive": df["agreement_positive"],
        "strong_subclone": df["strong_subclone"],
        "hp_assign_ge_099": df["hp_assign_rate"] >= 0.99,
        "globalp_le_050": df["GlobalP"] <= 0.50,
        "lowvaf_highadelta_lowcv": artifact_low_vaf_high_adelta_low_cv(df),
    }


def add_orthogonal_rows(config: DatasetConfig, eval_df: pd.DataFrame) -> List[Dict[str, object]]:
    meta = eval_df.iloc[0].to_dict() if not eval_df.empty else {}
    analyzed = eval_df[eval_df["analysis_available"]].copy()
    masks = orthogonal_feature_masks(analyzed)
    rows: List[Dict[str, object]] = []
    for name_a, name_b in itertools.combinations(masks.keys(), 2):
        mask_a = masks[name_a]
        mask_b = masks[name_b]
        intersection = mask_a & mask_b
        union = mask_a | mask_b
        a_only = mask_a & (~mask_b)
        b_only = mask_b & (~mask_a)
        union_count = int(union.sum())
        inter_count = int(intersection.sum())
        rows.append(
            {
                "dataset_id": config.dataset_id,
                "label": config.label,
                "sample": meta.get("sample", "unknown"),
                "mode": meta.get("mode", "unknown"),
                "feature_a": name_a,
                "feature_b": name_b,
                "jaccard": format_float(inter_count / union_count if union_count else 0.0),
                "a_only_tp": int(((analyzed["downstream_status"] == "caller_lost_tp") & a_only).sum()),
                "a_only_fp": int(((analyzed["downstream_status"] == "caller_removed_fp") & a_only).sum()),
                "b_only_tp": int(((analyzed["downstream_status"] == "caller_lost_tp") & b_only).sum()),
                "b_only_fp": int(((analyzed["downstream_status"] == "caller_removed_fp") & b_only).sum()),
                "both_tp": int(((analyzed["downstream_status"] == "caller_lost_tp") & intersection).sum()),
                "both_fp": int(((analyzed["downstream_status"] == "caller_removed_fp") & intersection).sum()),
                "union_tp": int(((analyzed["downstream_status"] == "caller_lost_tp") & union).sum()),
                "union_fp": int(((analyzed["downstream_status"] == "caller_removed_fp") & union).sum()),
            }
        )
    return rows


def summarize_best_rows(gq_rows: List[Dict[str, object]], methyl_rows: List[Dict[str, object]], combo_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    def best_by(rows: Iterable[Dict[str, object]], key_field: str = "delta_f1_vs_baseline") -> Dict[str, object] | None:
        items = list(rows)
        if not items:
            return None
        return max(items, key=lambda row: (float(row[key_field]), int(row["tp_rescued"]), -int(row["fp_reintroduced"])))

    summaries: List[Dict[str, object]] = []
    datasets = sorted({row["dataset_id"] for row in gq_rows + methyl_rows + combo_rows})
    for dataset_id in datasets:
        gq_subset = [row for row in gq_rows if row["dataset_id"] == dataset_id]
        methyl_subset = [row for row in methyl_rows if row["dataset_id"] == dataset_id]
        combo_subset = [row for row in combo_rows if row["dataset_id"] == dataset_id]
        best_gq = best_by(gq_subset)
        best_methyl = best_by(methyl_subset)
        best_combo = best_by(combo_subset, "delta_f1_vs_gate")
        for summary_type, row, note in [
            ("best_gq_only", best_gq, "依 baseline F1 最佳"),
            ("best_methyl_only", best_methyl, "只使用分析可得的甲基/標籤/phase 特徵"),
            ("best_combo_vs_gate", best_combo, "依相對同一 caller gate 的 delta F1 最佳"),
        ]:
            if row is None:
                continue
            summaries.append(
                {
                    "dataset_id": row["dataset_id"],
                    "label": row["label"],
                    "sample": row["sample"],
                    "mode": row["mode"],
                    "summary_type": summary_type,
                    "rule_id": row["rule_id"],
                    "tp_rescued": row["tp_rescued"],
                    "fp_reintroduced": row["fp_reintroduced"],
                    "delta_f1_vs_baseline": row["delta_f1_vs_baseline"],
                    "delta_f1_vs_gate": row.get("delta_f1_vs_gate", "nan"),
                    "notes": note,
                }
            )
    return summaries


def write_summary_markdown(
    output_dir: Path,
    pool_rows: List[Dict[str, object]],
    overview_rows: List[Dict[str, object]],
    gq_rows: List[Dict[str, object]],
    methyl_rows: List[Dict[str, object]],
    combo_rows: List[Dict[str, object]],
    artifact_rows: List[Dict[str, object]],
) -> None:
    lines: List[str] = []
    lines.append("# GQ / 甲基 rescue 系統性分析摘要")
    lines.append("")
    lines.append("## 破題結論")
    lines.append("")
    lines.append("1. `GQ` 目前仍是最穩定的 caller-first 主訊號。")
    lines.append("2. `ClairS-TO` 不一開始就放寬，不是因為 `gq>=10` 完全沒用，而是因為大多數被排除的 TO FP 根本不在 rescue-eligible 候選池裡；若全域放寬，會把大量非 somatic / PON 命中噪音一起放回來。")
    lines.append("3. `5kHz TO` 上確實存在有效甲基 support，但它目前仍是第二層訊號，尚未穩定超過最佳 `GQ only`。")
    lines.append("4. `5kHz paired` 與 `DORADO paired` 目前都支持 caller-first；甲基 rescue 尚未呈現跨樣本穩定規則。")
    lines.append("")

    pool_df = pd.DataFrame(pool_rows)
    overview_df = pd.DataFrame(overview_rows)
    gq_df = pd.DataFrame(gq_rows)
    methyl_df = pd.DataFrame(methyl_rows)
    combo_df = pd.DataFrame(combo_rows)
    artifact_df = pd.DataFrame(artifact_rows)

    if not pool_df.empty:
        lines.append("## 為何 ClairS-TO 不一開始就全域放寬 GQ")
        lines.append("")
        to_rows = pool_df[(pool_df["dataset_id"] == "hcc1395_5khz_to") & (pool_df["downstream_status"] == "caller_removed_fp")]
        if not to_rows.empty:
            row = to_rows.iloc[0]
            lines.append(
                f"- `HCC1395 5kHz TO` 的 `caller_removed_fp` 總量是 `{row['rows_total']}`，但 `candidate_eligible` 只有 `{row['candidate_eligible_rows']}`（`{row['candidate_eligible_fraction']}`）。"
            )
            lines.append(
                f"- 在被排除的 FP 中，`GQ>=10 / >=15 / >=20` 仍分別有 `{row['excluded_gq_ge_10']}` / `{row['excluded_gq_ge_15']}` / `{row['excluded_gq_ge_20']}` 筆。"
            )
            lines.append("- 這代表若不先用 caller 的 rescue-eligible 條件縮小候選池，而是直接全域放寬 GQ，會把大量 PON / NonSomatic / 非救援型噪音一起放回來。")
            lines.append("")

    lines.append("## 各資料集最佳規則")
    lines.append("")
    for dataset_id in overview_df["dataset_id"].tolist():
        dataset_overview = overview_df[overview_df["dataset_id"] == dataset_id].iloc[0]
        gq_best = gq_df[gq_df["dataset_id"] == dataset_id].sort_values("delta_f1_vs_baseline", ascending=False).iloc[0]
        methyl_best = methyl_df[methyl_df["dataset_id"] == dataset_id].sort_values("delta_f1_vs_baseline", ascending=False).iloc[0]
        combo_best = combo_df[combo_df["dataset_id"] == dataset_id].sort_values("delta_f1_vs_gate", ascending=False).iloc[0]
        lines.append(f"### {dataset_overview['label']}")
        lines.append("")
        lines.append(f"- baseline F1: `{dataset_overview['baseline_f1']}`")
        lines.append(f"- 最佳 `GQ only`: `{gq_best['rule_id']}`，`{gq_best['tp_rescued']} TP / {gq_best['fp_reintroduced']} FP`，`delta F1={gq_best['delta_f1_vs_baseline']}`")
        lines.append(f"- 最佳 `甲基 only`: `{methyl_best['rule_id']}`，`{methyl_best['tp_rescued']} TP / {methyl_best['fp_reintroduced']} FP`，`delta F1={methyl_best['delta_f1_vs_baseline']}`")
        lines.append(f"- 最佳 `GQ + 甲基`（相對 gate）: `{combo_best['rule_id']}`，`{combo_best['tp_rescued']} TP / {combo_best['fp_reintroduced']} FP`，`delta F1 vs gate={combo_best['delta_f1_vs_gate']}`")
        lines.append("")

    if not artifact_df.empty:
        lines.append("## Artifact 特徵角色")
        lines.append("")
        for dataset_id in artifact_df["dataset_id"].unique():
            subset = artifact_df[(artifact_df["dataset_id"] == dataset_id) & (artifact_df["role"] == "rescue_veto")].copy()
            if subset.empty:
                continue
            best = subset.sort_values("delta_f1_vs_gate", ascending=False).iloc[0]
            lines.append(
                f"- `{dataset_id}` 最佳 artifact veto 組合為 `{best['artifact_rule']} @ {best['base_gate']}`，相對 gate 的 `delta F1={best['delta_f1_vs_gate']}`；若仍為負值，表示它較適合作為後段 triage，而不是前段 TP rescue。"
            )
        lines.append("")

    lines.append("## 主要輸出")
    lines.append("")
    for name in [
        "pool_eligibility_summary.tsv",
        "exclusion_reason_summary.tsv",
        "dataset_overview.tsv",
        "gq_threshold_sweep.tsv",
        "methylation_only_rule_sweep.tsv",
        "gq_plus_methylation_rule_sweep.tsv",
        "artifact_role_comparison.tsv",
        "feature_interval_enrichment.tsv",
        "pairwise_feature_grid.tsv",
        "orthogonal_feature_candidates.tsv",
        "top_rule_summary.tsv",
    ]:
        path = (output_dir / name).resolve()
        lines.append(f"- [{name}]({path})")
    lines.append("")
    (output_dir / "analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    configs = [parse_dataset(value) for value in args.dataset]
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    pool_rows: List[Dict[str, object]] = []
    exclusion_rows: List[Dict[str, object]] = []
    overview_rows: List[Dict[str, object]] = []
    gq_rows: List[Dict[str, object]] = []
    methyl_rows: List[Dict[str, object]] = []
    combo_rows: List[Dict[str, object]] = []
    artifact_rows: List[Dict[str, object]] = []
    interval_rows: List[Dict[str, object]] = []
    grid_rows: List[Dict[str, object]] = []
    orthogonal_rows: List[Dict[str, object]] = []

    for config in configs:
        raw_df = load_raw_pool(config.raw_pool_tsv)
        joined_df = load_joined(config.joined_tsv)
        eval_df = candidate_eval_frame(joined_df)

        pool_rows.extend(add_pool_eligibility_rows(config, raw_df))
        exclusion_rows.extend(add_exclusion_reason_rows(config, raw_df))
        overview_rows.append(dataset_overview_row(config, raw_df, eval_df))
        gq_rows.extend(add_gq_sweep_rows(config, eval_df))
        methyl_rows.extend(add_methyl_only_rows(config, eval_df))
        combo_rows.extend(add_combo_rows(config, eval_df))
        artifact_rows.extend(add_artifact_rows(config, eval_df))
        interval_rows.extend(add_interval_rows(config, eval_df))
        grid_rows.extend(add_grid_rows(config, eval_df))
        orthogonal_rows.extend(add_orthogonal_rows(config, eval_df))

    summary_rows = summarize_best_rows(gq_rows, methyl_rows, combo_rows)

    write_tsv_rows(output_dir / "pool_eligibility_summary.tsv", POOL_ELIGIBILITY_FIELDS, pool_rows)
    write_tsv_rows(output_dir / "exclusion_reason_summary.tsv", EXCLUSION_REASON_FIELDS, exclusion_rows)
    write_tsv_rows(output_dir / "dataset_overview.tsv", DATASET_OVERVIEW_FIELDS, overview_rows)
    write_tsv_rows(output_dir / "gq_threshold_sweep.tsv", GQ_SWEEP_FIELDS, gq_rows)
    write_tsv_rows(output_dir / "methylation_only_rule_sweep.tsv", METHYL_ONLY_FIELDS, methyl_rows)
    write_tsv_rows(output_dir / "gq_plus_methylation_rule_sweep.tsv", COMBO_FIELDS, combo_rows)
    write_tsv_rows(output_dir / "artifact_role_comparison.tsv", ARTIFACT_FIELDS, artifact_rows)
    write_tsv_rows(output_dir / "feature_interval_enrichment.tsv", INTERVAL_FIELDS, interval_rows)
    write_tsv_rows(output_dir / "pairwise_feature_grid.tsv", GRID_FIELDS, grid_rows)
    write_tsv_rows(output_dir / "orthogonal_feature_candidates.tsv", ORTHOGONAL_FIELDS, orthogonal_rows)
    write_tsv_rows(output_dir / "top_rule_summary.tsv", SUMMARY_FIELDS, summary_rows)
    write_summary_markdown(output_dir, pool_rows, overview_rows, gq_rows, methyl_rows, combo_rows, artifact_rows)

    print(f"[analyze_gq_methylation_rescue_matrix] Wrote outputs to {output_dir}")


if __name__ == "__main__":
    main()
