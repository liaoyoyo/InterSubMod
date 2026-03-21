#!/usr/bin/env python3
"""Analyze methylation rescue feature space across candidate-specific joined feature tables."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Sequence

import pandas as pd

from research_common import compute_metrics, ensure_dir, markdown_table, write_tsv_rows


@dataclass
class DatasetConfig:
    dataset_id: str
    label: str
    input_tsv: Path
    baseline_tp: int
    baseline_fp: int
    truth_total: int


GateFunc = Callable[[pd.DataFrame], pd.Series]
PredicateFunc = Callable[[pd.DataFrame], pd.Series]


DATASET_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "caller",
    "mode",
    "total_rows",
    "candidate_rows",
    "candidate_tp",
    "candidate_fp",
    "analyzed_rows",
    "analyzed_tp",
    "analyzed_fp",
    "lost_tp_total",
    "removed_fp_total",
    "lost_tp_analyzed",
    "removed_fp_analyzed",
    "lost_tp_coverage",
    "removed_fp_coverage",
    "baseline_tp",
    "baseline_fp",
    "baseline_fn",
    "baseline_f1",
]

NUMERIC_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "subset",
    "truth_status",
    "feature",
    "n",
    "mean",
    "median",
    "p25",
    "p75",
    "p10",
    "p90",
]

CATEGORICAL_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "subset",
    "truth_status",
    "feature",
    "category",
    "count",
    "fraction",
]

GATE_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "gate_id",
    "gate_notes",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "fp_per_tp",
]

RULE_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "gate_id",
    "gate_notes",
    "feature_type",
    "feature_family",
    "rule_scope",
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
]

TOP_RULE_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "mode",
    "gate_id",
    "ranking",
    "rule_scope",
    "feature_family",
    "rule_id",
    "rule_notes",
    "tp_rescued",
    "fp_reintroduced",
    "delta_f1_vs_baseline",
    "delta_f1_vs_gate",
    "fp_per_tp",
]


NUMERIC_FEATURES = [
    "gq",
    "qual",
    "dp",
    "af",
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

CATEGORICAL_FEATURES = [
    "VerificationClass",
    "agreement_type",
    "cluster_class",
    "label_class",
    "DominantLabel",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        required=True,
        help="Format: dataset_id|label|joined_tsv|baseline_tp|baseline_fp|truth_total",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def parse_dataset(raw: str) -> DatasetConfig:
    parts = raw.split("|")
    if len(parts) != 6:
        raise ValueError(f"Invalid --dataset value: {raw}")
    dataset_id, label, input_tsv, baseline_tp, baseline_fp, truth_total = parts
    return DatasetConfig(
        dataset_id=dataset_id,
        label=label,
        input_tsv=Path(input_tsv).resolve(),
        baseline_tp=int(baseline_tp),
        baseline_fp=int(baseline_fp),
        truth_total=int(truth_total),
    )


def as_bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).map(lambda value: str(value).strip().lower() in {"1", "true", "yes", "y"})


def prepare_frame(path: Path) -> pd.DataFrame:
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
    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
        else:
            df[column] = math.nan

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
    for column in bool_columns:
        if column in df.columns:
            df[column] = as_bool_series(df[column])
        else:
            df[column] = False

    text_fill = {
        "VerificationClass": "NA",
        "agreement_type": "NA",
        "cluster_class": "NA",
        "label_class": "NA",
        "DominantLabel": "NA",
        "class_shift": "NA",
        "source_scope": "",
    }
    for column, default in text_fill.items():
        if column in df.columns:
            df[column] = df[column].fillna(default).replace("", default)
        else:
            df[column] = default

    df["truth_status"] = df["truth_status"].fillna("NA")
    df["downstream_status"] = df["downstream_status"].fillna("NA")
    df["candidate_eligible"] = df["candidate_eligible"].fillna(False)
    analyzed_mask = df["source_scope"].isin({"tp", "fp"})
    df["analysis_available"] = analyzed_mask
    df["agreement_positive"] = df["agreement_type"].isin({"label_upgrade", "consistent_strong", "consistent_subclone"})
    df["strong_subclone"] = df["VerificationClass"].isin({"Strong", "Subclone"})
    return df


def dataset_metadata(df: pd.DataFrame) -> Dict[str, str]:
    if df.empty:
        return {"sample": "unknown", "platform": "unknown", "caller": "unknown", "mode": "unknown"}
    row = df.iloc[0]
    return {
        "sample": str(row.get("sample", "unknown")),
        "platform": str(row.get("platform", "unknown")),
        "caller": str(row.get("caller", "unknown")),
        "mode": str(row.get("mode", "unknown")),
    }


def gate_definitions() -> List[Dict[str, object]]:
    return [
        {
            "gate_id": "candidate_only",
            "notes": "candidate_eligible",
            "func": lambda df: df["candidate_eligible"],
        },
        {
            "gate_id": "gq_ge_10",
            "notes": "candidate_eligible & gq>=10",
            "func": lambda df: df["candidate_eligible"] & (df["gq"] >= 10),
        },
        {
            "gate_id": "gq_ge_15",
            "notes": "candidate_eligible & gq>=15",
            "func": lambda df: df["candidate_eligible"] & (df["gq"] >= 15),
        },
        {
            "gate_id": "gq_ge_20",
            "notes": "candidate_eligible & gq>=20",
            "func": lambda df: df["candidate_eligible"] & (df["gq"] >= 20),
        },
        {
            "gate_id": "qual_ge_10_or_gq_ge_20",
            "notes": "candidate_eligible & (qual>=10 or gq>=20)",
            "func": lambda df: df["candidate_eligible"] & ((df["qual"] >= 10.0) | (df["gq"] >= 20)),
        },
    ]


def predicate_definitions() -> List[Dict[str, object]]:
    analyzed = lambda df: df["analysis_available"]
    return [
        {
            "feature_type": "methylation",
            "family": "pairwise_median",
            "rule_id": "pairwise_ge_015",
            "notes": "PairwiseMedianDist >= 0.15",
            "func": lambda df: analyzed(df) & (df["PairwiseMedianDist"] >= 0.15),
        },
        {
            "feature_type": "methylation",
            "family": "pairwise_median",
            "rule_id": "pairwise_ge_020",
            "notes": "PairwiseMedianDist >= 0.20",
            "func": lambda df: analyzed(df) & (df["PairwiseMedianDist"] >= 0.20),
        },
        {
            "feature_type": "methylation",
            "family": "pairwise_median",
            "rule_id": "pairwise_ge_025",
            "notes": "PairwiseMedianDist >= 0.25",
            "func": lambda df: analyzed(df) & (df["PairwiseMedianDist"] >= 0.25),
        },
        {
            "feature_type": "methylation",
            "family": "pairwise_median",
            "rule_id": "pairwise_le_020",
            "notes": "PairwiseMedianDist <= 0.20",
            "func": lambda df: analyzed(df) & (df["PairwiseMedianDist"] <= 0.20),
        },
        {
            "feature_type": "methylation",
            "family": "pairwise_median",
            "rule_id": "pairwise_le_025",
            "notes": "PairwiseMedianDist <= 0.25",
            "func": lambda df: analyzed(df) & (df["PairwiseMedianDist"] <= 0.25),
        },
        {
            "feature_type": "methylation",
            "family": "pairwise_mean",
            "rule_id": "pmean_ge_020",
            "notes": "PairwiseMeanDist >= 0.20",
            "func": lambda df: analyzed(df) & (df["PairwiseMeanDist"] >= 0.20),
        },
        {
            "feature_type": "methylation",
            "family": "pairwise_mean",
            "rule_id": "pmean_le_020",
            "notes": "PairwiseMeanDist <= 0.20",
            "func": lambda df: analyzed(df) & (df["PairwiseMeanDist"] <= 0.20),
        },
        {
            "feature_type": "methylation",
            "family": "allele_delta",
            "rule_id": "allele_delta_ge_001",
            "notes": "AlleleDelta >= 0.01",
            "func": lambda df: analyzed(df) & (df["AlleleDelta"] >= 0.01),
        },
        {
            "feature_type": "methylation",
            "family": "allele_delta",
            "rule_id": "allele_delta_ge_005",
            "notes": "AlleleDelta >= 0.05",
            "func": lambda df: analyzed(df) & (df["AlleleDelta"] >= 0.05),
        },
        {
            "feature_type": "methylation",
            "family": "allele_delta",
            "rule_id": "allele_delta_ge_015",
            "notes": "AlleleDelta >= 0.15",
            "func": lambda df: analyzed(df) & (df["AlleleDelta"] >= 0.15),
        },
        {
            "feature_type": "methylation",
            "family": "allele_delta",
            "rule_id": "allele_delta_le_001",
            "notes": "AlleleDelta <= 0.01",
            "func": lambda df: analyzed(df) & (df["AlleleDelta"] <= 0.01),
        },
        {
            "feature_type": "methylation",
            "family": "cramers_v",
            "rule_id": "cramersv_le_005",
            "notes": "CramersV <= 0.05",
            "func": lambda df: analyzed(df) & (df["CramersV"] <= 0.05),
        },
        {
            "feature_type": "methylation",
            "family": "global_p",
            "rule_id": "globalp_le_005",
            "notes": "GlobalP <= 0.05",
            "func": lambda df: analyzed(df) & (df["GlobalP"] <= 0.05),
        },
        {
            "feature_type": "methylation",
            "family": "global_p",
            "rule_id": "globalp_le_050",
            "notes": "GlobalP <= 0.50",
            "func": lambda df: analyzed(df) & (df["GlobalP"] <= 0.50),
        },
        {
            "feature_type": "methylation",
            "family": "quality_score",
            "rule_id": "quality_ge_60",
            "notes": "Quality_Score >= 60",
            "func": lambda df: analyzed(df) & (df["Quality_Score"] >= 60),
        },
        {
            "feature_type": "methylation",
            "family": "quality_score",
            "rule_id": "quality_ge_75",
            "notes": "Quality_Score >= 75",
            "func": lambda df: analyzed(df) & (df["Quality_Score"] >= 75),
        },
        {
            "feature_type": "phase_label",
            "family": "hp_assign_rate",
            "rule_id": "hp_assign_ge_090",
            "notes": "hp_assign_rate >= 0.90",
            "func": lambda df: analyzed(df) & (df["hp_assign_rate"] >= 0.90),
        },
        {
            "feature_type": "phase_label",
            "family": "hp_assign_rate",
            "rule_id": "hp_assign_ge_099",
            "notes": "hp_assign_rate >= 0.99",
            "func": lambda df: analyzed(df) & (df["hp_assign_rate"] >= 0.99),
        },
        {
            "feature_type": "phase_label",
            "family": "allele_assign_rate",
            "rule_id": "allele_assign_ge_099",
            "notes": "allele_assign_rate >= 0.99",
            "func": lambda df: analyzed(df) & (df["allele_assign_rate"] >= 0.99),
        },
        {
            "feature_type": "label_cluster",
            "family": "verification_class",
            "rule_id": "strong_subclone",
            "notes": "VerificationClass in Strong/Subclone",
            "func": lambda df: analyzed(df) & df["strong_subclone"],
        },
        {
            "feature_type": "label_cluster",
            "family": "agreement_type",
            "rule_id": "agreement_positive",
            "notes": "agreement_type in label_upgrade/consistent_strong/consistent_subclone",
            "func": lambda df: analyzed(df) & df["agreement_positive"],
        },
        {
            "feature_type": "label_cluster",
            "family": "agreement_type",
            "rule_id": "consistent_strong",
            "notes": "agreement_type == consistent_strong",
            "func": lambda df: analyzed(df) & (df["agreement_type"] == "consistent_strong"),
        },
        {
            "feature_type": "label_cluster",
            "family": "agreement_type",
            "rule_id": "label_upgrade",
            "notes": "agreement_type == label_upgrade",
            "func": lambda df: analyzed(df) & (df["agreement_type"] == "label_upgrade"),
        },
        {
            "feature_type": "label_cluster",
            "family": "passed_gating",
            "rule_id": "passed_gating",
            "notes": "PassedGating == true",
            "func": lambda df: analyzed(df) & df["PassedGating"],
        },
    ]


def top_predicates_for_combinations() -> Sequence[str]:
    return [
        "pairwise_ge_020",
        "pairwise_le_020",
        "globalp_le_050",
        "quality_ge_60",
        "quality_ge_75",
        "hp_assign_ge_099",
        "strong_subclone",
        "agreement_positive",
        "label_upgrade",
        "allele_delta_ge_005",
    ]


def coverage_fraction(num: int, den: int) -> str:
    if den == 0:
        return "0.000000"
    return f"{num / den:.6f}"


def summarize_dataset(config: DatasetConfig, df: pd.DataFrame, meta: Dict[str, str]) -> Dict[str, object]:
    candidate = df[df["candidate_eligible"]].copy()
    analyzed = candidate[candidate["analysis_available"]].copy()
    lost_tp_total = int((candidate["downstream_status"] == "caller_lost_tp").sum())
    removed_fp_total = int((candidate["downstream_status"] == "caller_removed_fp").sum())
    lost_tp_analyzed = int(((analyzed["downstream_status"] == "caller_lost_tp")).sum())
    removed_fp_analyzed = int(((analyzed["downstream_status"] == "caller_removed_fp")).sum())
    baseline = compute_metrics(config.baseline_tp, config.baseline_fp, config.truth_total)
    return {
        "dataset_id": config.dataset_id,
        "label": config.label,
        "sample": meta["sample"],
        "platform": meta["platform"],
        "caller": meta["caller"],
        "mode": meta["mode"],
        "total_rows": len(df),
        "candidate_rows": len(candidate),
        "candidate_tp": int((candidate["truth_status"] == "TP").sum()),
        "candidate_fp": int((candidate["truth_status"] == "FP").sum()),
        "analyzed_rows": len(analyzed),
        "analyzed_tp": int((analyzed["truth_status"] == "TP").sum()),
        "analyzed_fp": int((analyzed["truth_status"] == "FP").sum()),
        "lost_tp_total": lost_tp_total,
        "removed_fp_total": removed_fp_total,
        "lost_tp_analyzed": lost_tp_analyzed,
        "removed_fp_analyzed": removed_fp_analyzed,
        "lost_tp_coverage": coverage_fraction(lost_tp_analyzed, lost_tp_total),
        "removed_fp_coverage": coverage_fraction(removed_fp_analyzed, removed_fp_total),
        "baseline_tp": config.baseline_tp,
        "baseline_fp": config.baseline_fp,
        "baseline_fn": config.truth_total - config.baseline_tp,
        "baseline_f1": f"{baseline['f1']:.6f}",
    }


def numeric_distributions(config: DatasetConfig, df: pd.DataFrame, meta: Dict[str, str]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    subsets = {
        "candidate_analyzed": df[df["candidate_eligible"] & df["analysis_available"]].copy(),
    }
    for subset_name, subset in subsets.items():
        for truth_status in ("TP", "FP"):
            group = subset[subset["truth_status"] == truth_status].copy()
            for feature in NUMERIC_FEATURES:
                values = group[feature].dropna()
                if values.empty:
                    continue
                rows.append(
                    {
                        "dataset_id": config.dataset_id,
                        "label": config.label,
                        "sample": meta["sample"],
                        "mode": meta["mode"],
                        "subset": subset_name,
                        "truth_status": truth_status,
                        "feature": feature,
                        "n": len(values),
                        "mean": f"{values.mean():.6f}",
                        "median": f"{values.median():.6f}",
                        "p25": f"{values.quantile(0.25):.6f}",
                        "p75": f"{values.quantile(0.75):.6f}",
                        "p10": f"{values.quantile(0.10):.6f}",
                        "p90": f"{values.quantile(0.90):.6f}",
                    }
                )
    return rows


def categorical_distributions(config: DatasetConfig, df: pd.DataFrame, meta: Dict[str, str]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    subset = df[df["candidate_eligible"] & df["analysis_available"]].copy()
    for truth_status in ("TP", "FP"):
        group = subset[subset["truth_status"] == truth_status].copy()
        total = len(group)
        if total == 0:
            continue
        for feature in CATEGORICAL_FEATURES:
            counts = group[feature].fillna("NA").value_counts(dropna=False)
            for category, count in counts.items():
                rows.append(
                    {
                        "dataset_id": config.dataset_id,
                        "label": config.label,
                        "sample": meta["sample"],
                        "mode": meta["mode"],
                        "subset": "candidate_analyzed",
                        "truth_status": truth_status,
                        "feature": feature,
                        "category": str(category),
                        "count": int(count),
                        "fraction": f"{count / total:.6f}",
                    }
                )
    return rows


def evaluate_rule(
    df: pd.DataFrame,
    mask: pd.Series,
    baseline_tp: int,
    baseline_fp: int,
    truth_total: int,
    baseline_f1: float,
) -> Dict[str, object]:
    tp_rescued = int(((df["truth_status"] == "TP") & mask).sum())
    fp_reintroduced = int(((df["truth_status"] == "FP") & mask).sum())
    metrics = compute_metrics(baseline_tp + tp_rescued, baseline_fp + fp_reintroduced, truth_total)
    delta = metrics["f1"] - baseline_f1
    fp_per_tp = fp_reintroduced / tp_rescued if tp_rescued else float("inf")
    return {
        "trigger_count": tp_rescued + fp_reintroduced,
        "tp_rescued": tp_rescued,
        "fp_reintroduced": fp_reintroduced,
        "precision": f"{metrics['precision']:.6f}",
        "recall": f"{metrics['recall']:.6f}",
        "f1": f"{metrics['f1']:.6f}",
        "delta_f1_vs_baseline": delta,
        "fp_per_tp": fp_per_tp,
    }


def format_fp_per_tp(value: float) -> str:
    if math.isinf(value):
        return "inf"
    return f"{value:.6f}"


def sweep_rules(
    config: DatasetConfig,
    df: pd.DataFrame,
    meta: Dict[str, str],
) -> tuple[list[dict], list[dict], list[dict]]:
    baseline = compute_metrics(config.baseline_tp, config.baseline_fp, config.truth_total)
    gate_rows: List[Dict[str, object]] = []
    single_rows: List[Dict[str, object]] = []
    combo_rows: List[Dict[str, object]] = []

    gates = gate_definitions()
    predicates = predicate_definitions()
    predicate_lookup = {item["rule_id"]: item for item in predicates}
    selected_for_combo = [predicate_lookup[item] for item in top_predicates_for_combinations()]

    gate_masks: Dict[str, pd.Series] = {}
    gate_metrics: Dict[str, Dict[str, object]] = {}

    for gate in gates:
        gate_mask = gate["func"](df)
        gate_masks[gate["gate_id"]] = gate_mask
        metrics = evaluate_rule(
            df,
            gate_mask,
            config.baseline_tp,
            config.baseline_fp,
            config.truth_total,
            baseline["f1"],
        )
        gate_metrics[gate["gate_id"]] = metrics
        gate_rows.append(
            {
                "dataset_id": config.dataset_id,
                "label": config.label,
                "sample": meta["sample"],
                "mode": meta["mode"],
                "gate_id": gate["gate_id"],
                "gate_notes": gate["notes"],
                "tp_rescued": metrics["tp_rescued"],
                "fp_reintroduced": metrics["fp_reintroduced"],
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1": metrics["f1"],
                "delta_f1_vs_baseline": f"{metrics['delta_f1_vs_baseline']:.6f}",
                "fp_per_tp": format_fp_per_tp(metrics["fp_per_tp"]),
            }
        )

    for gate in gates:
        gate_mask = gate_masks[gate["gate_id"]]
        gate_delta = gate_metrics[gate["gate_id"]]["delta_f1_vs_baseline"]
        for predicate in predicates:
            mask = gate_mask & predicate["func"](df)
            metrics = evaluate_rule(
                df,
                mask,
                config.baseline_tp,
                config.baseline_fp,
                config.truth_total,
                baseline["f1"],
            )
            single_rows.append(
                {
                    "dataset_id": config.dataset_id,
                    "label": config.label,
                    "sample": meta["sample"],
                    "mode": meta["mode"],
                    "gate_id": gate["gate_id"],
                    "gate_notes": gate["notes"],
                    "feature_type": predicate["feature_type"],
                    "feature_family": predicate["family"],
                    "rule_scope": "single_feature",
                    "rule_id": predicate["rule_id"],
                    "rule_notes": predicate["notes"],
                    "trigger_count": metrics["trigger_count"],
                    "tp_rescued": metrics["tp_rescued"],
                    "fp_reintroduced": metrics["fp_reintroduced"],
                    "precision": metrics["precision"],
                    "recall": metrics["recall"],
                    "f1": metrics["f1"],
                    "delta_f1_vs_baseline": f"{metrics['delta_f1_vs_baseline']:.6f}",
                    "delta_f1_vs_gate": f"{(metrics['delta_f1_vs_baseline'] - gate_delta):.6f}",
                    "fp_per_tp": format_fp_per_tp(metrics["fp_per_tp"]),
                }
            )

        for idx, predicate_a in enumerate(selected_for_combo):
            for predicate_b in selected_for_combo[idx + 1 :]:
                if predicate_a["family"] == predicate_b["family"]:
                    continue
                mask = gate_mask & predicate_a["func"](df) & predicate_b["func"](df)
                metrics = evaluate_rule(
                    df,
                    mask,
                    config.baseline_tp,
                    config.baseline_fp,
                    config.truth_total,
                    baseline["f1"],
                )
                combo_rows.append(
                    {
                        "dataset_id": config.dataset_id,
                        "label": config.label,
                        "sample": meta["sample"],
                        "mode": meta["mode"],
                        "gate_id": gate["gate_id"],
                        "gate_notes": gate["notes"],
                        "feature_type": f"{predicate_a['feature_type']}+{predicate_b['feature_type']}",
                        "feature_family": f"{predicate_a['family']}+{predicate_b['family']}",
                        "rule_scope": "feature_combo",
                        "rule_id": f"{predicate_a['rule_id']}__and__{predicate_b['rule_id']}",
                        "rule_notes": f"{predicate_a['notes']} AND {predicate_b['notes']}",
                        "trigger_count": metrics["trigger_count"],
                        "tp_rescued": metrics["tp_rescued"],
                        "fp_reintroduced": metrics["fp_reintroduced"],
                        "precision": metrics["precision"],
                        "recall": metrics["recall"],
                        "f1": metrics["f1"],
                        "delta_f1_vs_baseline": f"{metrics['delta_f1_vs_baseline']:.6f}",
                        "delta_f1_vs_gate": f"{(metrics['delta_f1_vs_baseline'] - gate_delta):.6f}",
                        "fp_per_tp": format_fp_per_tp(metrics["fp_per_tp"]),
                    }
                )

    top_rows: List[Dict[str, object]] = []
    gate_ids = [item["gate_id"] for item in gates]
    combined = pd.DataFrame(single_rows + combo_rows)
    if not combined.empty:
        combined["tp_rescued"] = pd.to_numeric(combined["tp_rescued"], errors="coerce").fillna(0).astype(int)
        combined["fp_reintroduced"] = pd.to_numeric(combined["fp_reintroduced"], errors="coerce").fillna(0).astype(int)
        combined["delta_f1_vs_baseline"] = pd.to_numeric(combined["delta_f1_vs_baseline"], errors="coerce").fillna(0.0)
        combined["delta_f1_vs_gate"] = pd.to_numeric(combined["delta_f1_vs_gate"], errors="coerce").fillna(0.0)
        combined["fp_per_tp_numeric"] = pd.to_numeric(
            combined["fp_per_tp"].replace({"inf": math.inf}), errors="coerce"
        ).fillna(math.inf)
        for gate_id in gate_ids:
            subset = combined[combined["gate_id"] == gate_id].copy()
            if subset.empty:
                continue
            best_overall = subset.sort_values(
                by=["delta_f1_vs_baseline", "tp_rescued", "fp_per_tp_numeric"],
                ascending=[False, False, True],
            ).head(5)
            for _, row in best_overall.iterrows():
                top_rows.append(
                    {
                        "dataset_id": row["dataset_id"],
                        "label": row["label"],
                        "sample": row["sample"],
                        "mode": row["mode"],
                        "gate_id": gate_id,
                        "ranking": "best_overall_delta",
                        "rule_scope": row["rule_scope"],
                        "feature_family": row["feature_family"],
                        "rule_id": row["rule_id"],
                        "rule_notes": row["rule_notes"],
                        "tp_rescued": int(row["tp_rescued"]),
                        "fp_reintroduced": int(row["fp_reintroduced"]),
                        "delta_f1_vs_baseline": f"{row['delta_f1_vs_baseline']:.6f}",
                        "delta_f1_vs_gate": f"{row['delta_f1_vs_gate']:.6f}",
                        "fp_per_tp": row["fp_per_tp"],
                    }
                )
            best_incremental = subset.sort_values(
                by=["delta_f1_vs_gate", "tp_rescued", "fp_per_tp_numeric"],
                ascending=[False, False, True],
            ).head(5)
            for _, row in best_incremental.iterrows():
                top_rows.append(
                    {
                        "dataset_id": row["dataset_id"],
                        "label": row["label"],
                        "sample": row["sample"],
                        "mode": row["mode"],
                        "gate_id": gate_id,
                        "ranking": "best_incremental_vs_gate",
                        "rule_scope": row["rule_scope"],
                        "feature_family": row["feature_family"],
                        "rule_id": row["rule_id"],
                        "rule_notes": row["rule_notes"],
                        "tp_rescued": int(row["tp_rescued"]),
                        "fp_reintroduced": int(row["fp_reintroduced"]),
                        "delta_f1_vs_baseline": f"{row['delta_f1_vs_baseline']:.6f}",
                        "delta_f1_vs_gate": f"{row['delta_f1_vs_gate']:.6f}",
                        "fp_per_tp": row["fp_per_tp"],
                    }
                )
    return gate_rows, single_rows, combo_rows, top_rows


def build_summary_markdown(
    dataset_rows: Sequence[Dict[str, object]],
    gate_rows: Sequence[Dict[str, object]],
    top_rows: Sequence[Dict[str, object]],
    output_dir: Path,
) -> None:
    dataset_df = pd.DataFrame(dataset_rows)
    gate_df = pd.DataFrame(gate_rows)
    top_df = pd.DataFrame(top_rows)

    lines = [
        "# Methylation Rescue Feature Space Summary",
        "",
        "## Dataset coverage",
        "",
    ]
    coverage_columns = [
        "label",
        "mode",
        "candidate_rows",
        "analyzed_rows",
        "lost_tp_analyzed",
        "lost_tp_coverage",
        "removed_fp_analyzed",
        "removed_fp_coverage",
        "baseline_f1",
    ]
    lines.append(markdown_table(coverage_columns, dataset_df[coverage_columns].to_dict("records")))
    lines.extend(["", "## Gate baselines", ""])
    gate_columns = [
        "label",
        "mode",
        "gate_id",
        "tp_rescued",
        "fp_reintroduced",
        "delta_f1_vs_baseline",
        "fp_per_tp",
    ]
    lines.append(markdown_table(gate_columns, gate_df[gate_columns].to_dict("records")))
    lines.extend(["", "## Top feature rules by gate", ""])
    if top_df.empty:
        lines.append("_No top rules_")
    else:
        top_columns = [
            "label",
            "mode",
            "gate_id",
            "ranking",
            "rule_id",
            "tp_rescued",
            "fp_reintroduced",
            "delta_f1_vs_baseline",
            "delta_f1_vs_gate",
        ]
        lines.append(markdown_table(top_columns, top_df[top_columns].to_dict("records")))
    (output_dir / "feature_space_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    configs = [parse_dataset(item) for item in args.dataset]
    dataset_rows: List[Dict[str, object]] = []
    numeric_rows: List[Dict[str, object]] = []
    categorical_rows: List[Dict[str, object]] = []
    gate_rows: List[Dict[str, object]] = []
    single_rows: List[Dict[str, object]] = []
    combo_rows: List[Dict[str, object]] = []
    top_rows: List[Dict[str, object]] = []

    for config in configs:
        df = prepare_frame(config.input_tsv)
        meta = dataset_metadata(df)
        dataset_rows.append(summarize_dataset(config, df, meta))
        numeric_rows.extend(numeric_distributions(config, df, meta))
        categorical_rows.extend(categorical_distributions(config, df, meta))
        gate_part, single_part, combo_part, top_part = sweep_rules(config, df, meta)
        gate_rows.extend(gate_part)
        single_rows.extend(single_part)
        combo_rows.extend(combo_part)
        top_rows.extend(top_part)

    write_tsv_rows(output_dir / "dataset_summary.tsv", DATASET_FIELDS, dataset_rows)
    write_tsv_rows(output_dir / "numeric_feature_distribution.tsv", NUMERIC_FIELDS, numeric_rows)
    write_tsv_rows(output_dir / "categorical_feature_distribution.tsv", CATEGORICAL_FIELDS, categorical_rows)
    write_tsv_rows(output_dir / "gate_baselines.tsv", GATE_FIELDS, gate_rows)
    write_tsv_rows(output_dir / "single_feature_rule_sweep.tsv", RULE_FIELDS, single_rows)
    write_tsv_rows(output_dir / "combination_rule_sweep.tsv", RULE_FIELDS, combo_rows)
    write_tsv_rows(output_dir / "top_rules_by_gate.tsv", TOP_RULE_FIELDS, top_rows)
    build_summary_markdown(dataset_rows, gate_rows, top_rows, output_dir)

    print(f"[analyze_methylation_rescue_feature_space] Wrote {output_dir / 'dataset_summary.tsv'}")


if __name__ == "__main__":
    main()
