#!/usr/bin/env python3
"""Build and evaluate analysis-layer annotations from phase 2 finer interval results."""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd

from research_common import compute_metrics, ensure_dir, markdown_table, write_tsv_rows
from run_phase2_paired_model_feature_analysis import (
    RESCUE_DATASETS,
    candidate_eval_frame,
    evaluate_mask,
    load_joined,
)


PHASE2_OUTPUT_DIR = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
    "20260311_phase2_paired_model_feature_analysis"
)
OUTPUT_ROOT_DEFAULT = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
    "20260311_phase2_annotation_layer"
)


@dataclass(frozen=True)
class DatasetAnnotationConfig:
    dataset_id: str
    caller_primary_gq: int
    caller_strict_gq: int
    support_topbin_features: Tuple[str, ...]
    support_flag_features: Tuple[str, ...]
    weak_support_flag_features: Tuple[str, ...]
    qc_topbin_features: Tuple[str, ...]
    qc_flag_features: Tuple[str, ...]
    observed_only_topbin_features: Tuple[str, ...]
    notes: str


ANNOTATION_CONFIG: Dict[str, DatasetAnnotationConfig] = {
    "hcc1395_5khz_to": DatasetAnnotationConfig(
        dataset_id="hcc1395_5khz_to",
        caller_primary_gq=10,
        caller_strict_gq=15,
        support_topbin_features=("Quality_Score", "PairwiseMedianDist"),
        support_flag_features=("agreement_positive",),
        weak_support_flag_features=("cluster_first_support",),
        qc_topbin_features=(),
        qc_flag_features=("hp_assign_high_095", "hp_assign_veryhigh_099"),
        observed_only_topbin_features=("hp_assign_rate",),
        notes="TO 主線：caller-first + quality/pairwise/agreement support；hp_assign 保留 phase/QC。",
    ),
    "hcc1395_dorado_to": DatasetAnnotationConfig(
        dataset_id="hcc1395_dorado_to",
        caller_primary_gq=10,
        caller_strict_gq=15,
        support_topbin_features=("Quality_Score", "PairwiseMedianDist"),
        support_flag_features=(),
        weak_support_flag_features=("cluster_first_support",),
        qc_topbin_features=(),
        qc_flag_features=("hp_assign_high_095", "hp_assign_veryhigh_099"),
        observed_only_topbin_features=("hp_assign_rate",),
        notes="DORADO TO：quality 與低 pairwise 作 dataset-aware support；hp_assign 只作 QC。",
    ),
    "hcc1395_5khz_paired": DatasetAnnotationConfig(
        dataset_id="hcc1395_5khz_paired",
        caller_primary_gq=15,
        caller_strict_gq=20,
        support_topbin_features=(),
        support_flag_features=(),
        weak_support_flag_features=(),
        qc_topbin_features=("Quality_Score", "PairwiseMedianDist", "hp_assign_rate"),
        qc_flag_features=("hp_assign_high_095", "hp_assign_veryhigh_099"),
        observed_only_topbin_features=(),
        notes="5kHz paired：甲基訊號暫不升級為 support，只保留 annotation/QC。",
    ),
    "hcc1395_dorado_paired": DatasetAnnotationConfig(
        dataset_id="hcc1395_dorado_paired",
        caller_primary_gq=15,
        caller_strict_gq=20,
        support_topbin_features=(),
        support_flag_features=(),
        weak_support_flag_features=("cluster_first_support",),
        qc_topbin_features=("Quality_Score", "PairwiseMedianDist", "hp_assign_rate"),
        qc_flag_features=("hp_assign_high_095", "hp_assign_veryhigh_099"),
        observed_only_topbin_features=(),
        notes="DORADO paired：caller-first 仍主導；高 hp_assign/quality/低 pairwise 只保留弱 support 或 QC。",
    ),
}


ANNOTATION_CONFIG_FIELDS = [
    "dataset_id",
    "verification_support_mode",
    "caller_primary_gq",
    "caller_strict_gq",
    "support_topbin_features",
    "support_flag_features",
    "weak_support_flag_features",
    "qc_topbin_features",
    "qc_flag_features",
    "observed_only_topbin_features",
    "notes",
]

ANNOTATED_FIELDS = [
    "dataset_id",
    "label",
    "sample",
    "platform",
    "mode",
    "caller",
    "region_key",
    "downstream_status",
    "source_scope",
    "candidate_eligible",
    "analysis_available",
    "qual",
    "gq",
    "dp",
    "af",
    "ad_ref",
    "ad_alt",
    "PairwiseMedianDist",
    "Quality_Score",
    "AlleleDelta",
    "CramersV",
    "hp_assign_rate",
    "VerificationClass",
    "label_first_support",
    "cluster_first_support",
    "verification_support_mode",
    "verification_support_source",
    "verification_support_schema_status",
    "agreement_type",
    "class_shift",
    "caller_primary_flag",
    "caller_strict_flag",
    "quality_topbin_flag",
    "pairwise_topbin_flag",
    "hp_assign_topbin_flag",
    "agreement_positive_flag",
    "cluster_first_support_flag",
    "artifact_lowvaf_highad_flag",
    "artifact_lowvaf_highad_lowcv_flag",
    "hp_assign_high_095",
    "hp_assign_veryhigh_099",
    "support_score_primary",
    "support_score_all",
    "qc_score",
    "annotation_support_any",
    "annotation_support_primary_any",
    "annotation_qc_any",
    "annotation_topbin_summary",
    "annotation_support_summary",
    "annotation_qc_summary",
    "annotation_role_summary",
]

SUMMARY_FIELDS = [
    "dataset_id",
    "label",
    "annotation_name",
    "annotation_role",
    "trigger_count",
    "tp_count",
    "fp_count",
    "tp_fraction",
    "fp_fraction",
    "tp_to_fp_ratio",
    "enrichment_ratio",
]

POLICY_FIELDS = [
    "dataset_id",
    "label",
    "policy_id",
    "policy_notes",
    "trigger_count",
    "tp_rescued",
    "fp_reintroduced",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_baseline",
    "delta_f1_vs_primary_gate",
]

BEST_POLICY_FIELDS = [
    "dataset_id",
    "label",
    "best_policy_by_baseline",
    "best_delta_f1_vs_baseline",
    "best_policy_vs_primary_gate",
    "best_delta_f1_vs_primary_gate",
]

OVERLAP_FIELDS = [
    "dataset_id",
    "label",
    "feature_a",
    "feature_b",
    "jaccard",
    "a_only_tp",
    "a_only_fp",
    "b_only_tp",
    "b_only_fp",
    "both_tp",
    "both_fp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase2-dir", default=str(PHASE2_OUTPUT_DIR), help="Directory containing phase 2 outputs")
    parser.add_argument("--output-dir", default=str(OUTPUT_ROOT_DEFAULT), help="Output directory")
    parser.add_argument(
        "--verification-support-mode",
        choices=("evidence-v2", "legacy"),
        required=True,
        help="Explicit source for the historical Strong/Subclone support cohort",
    )
    return parser.parse_args()


def format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf"
    return f"{value:.6f}"


def parse_bin_label(label: str) -> Tuple[float, float]:
    match = re.fullmatch(r"\[(.*?),(.*?)\)", str(label).strip())
    if not match:
        raise ValueError(f"Unsupported bin label: {label}")
    left_text, right_text = match.groups()
    left = float(left_text)
    right = float("inf") if right_text == "inf" else float(right_text)
    return left, right


def in_interval(series: pd.Series, label: str) -> pd.Series:
    left, right = parse_bin_label(label)
    mask = series >= left
    if math.isinf(right):
        return mask
    return mask & (series < right)


def low_vaf_high_ad(df: pd.DataFrame) -> pd.Series:
    return (df["af"] < 0.24) & (df["AlleleDelta"] > 0.25)


def low_vaf_high_ad_low_cv(df: pd.DataFrame) -> pd.Series:
    return low_vaf_high_ad(df) & (df["CramersV"] < 0.05)


def load_top_bins(phase2_dir: Path) -> Dict[Tuple[str, str], str]:
    df = pd.read_csv(phase2_dir / "fine_feature_interval_top_bins.tsv", sep="\t")
    selected = df[df["feature"].isin(["Quality_Score", "PairwiseMedianDist", "hp_assign_rate"])].copy()
    return {(row["dataset_id"], row["feature"]): row["bin_label"] for _, row in selected.iterrows()}


def annotation_mask(df: pd.DataFrame, flag_name: str) -> pd.Series:
    direct_flag_column = f"{flag_name}_flag"
    if direct_flag_column in df.columns:
        return df[direct_flag_column].fillna(False).astype(bool)
    if flag_name in df.columns:
        return df[flag_name].fillna(False).astype(bool)
    if flag_name == "agreement_positive":
        return df["agreement_positive_flag"]
    if flag_name == "cluster_first_support":
        return df["cluster_first_support_flag"]
    if flag_name == "hp_assign_high_095":
        return df["hp_assign_high_095"]
    if flag_name == "hp_assign_veryhigh_099":
        return df["hp_assign_veryhigh_099"]
    raise KeyError(flag_name)


def build_annotation_strings(row: pd.Series, support_names: Sequence[str], qc_names: Sequence[str], observed_names: Sequence[str]) -> Tuple[str, str, str, str]:
    def row_flag(name: str) -> bool:
        if f"{name}_flag" in row.index:
            return bool(row.get(f"{name}_flag", False))
        return bool(row.get(name, False))

    topbin_tokens = [name for name in ["quality_topbin", "pairwise_topbin", "hp_assign_topbin"] if row_flag(name)]
    support_tokens = [name for name in support_names if row_flag(name)]
    qc_tokens = [name for name in qc_names if row_flag(name)]
    observed_tokens = [name for name in observed_names if row_flag(name)]
    topbin_summary = ",".join(topbin_tokens) if topbin_tokens else "none"
    support_summary = ",".join(support_tokens) if support_tokens else "none"
    qc_summary = ",".join(qc_tokens) if qc_tokens else "none"
    role_tokens = []
    if support_tokens:
        role_tokens.append("support")
    if qc_tokens:
        role_tokens.append("qc")
    if observed_tokens:
        role_tokens.append("observed")
    role_summary = ",".join(role_tokens) if role_tokens else "none"
    return topbin_summary, support_summary, qc_summary, role_summary


def annotate_dataset(
    dataset,
    top_bins: Dict[Tuple[str, str], str],
    support_mode: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    config = ANNOTATION_CONFIG[dataset.dataset_id]
    joined_df = load_joined(dataset.joined_tsv, support_mode)
    base_df = joined_df[joined_df["candidate_eligible"]].copy()
    base_df["dataset_id"] = dataset.dataset_id
    base_df["label"] = dataset.label
    base_df["analysis_available"] = base_df["source_scope"].isin({"tp", "fp"})

    quality_bin = top_bins.get((dataset.dataset_id, "Quality_Score"))
    pairwise_bin = top_bins.get((dataset.dataset_id, "PairwiseMedianDist"))
    hp_bin = top_bins.get((dataset.dataset_id, "hp_assign_rate"))

    base_df["caller_primary_flag"] = base_df["gq"] >= config.caller_primary_gq
    base_df["caller_strict_flag"] = base_df["gq"] >= config.caller_strict_gq
    base_df["quality_topbin_flag"] = in_interval(base_df["Quality_Score"], quality_bin) if quality_bin else False
    base_df["pairwise_topbin_flag"] = in_interval(base_df["PairwiseMedianDist"], pairwise_bin) if pairwise_bin else False
    base_df["hp_assign_topbin_flag"] = in_interval(base_df["hp_assign_rate"], hp_bin) if hp_bin else False
    base_df["agreement_positive_flag"] = base_df["agreement_type"].isin(
        ["label_upgrade", "consistent_strong", "consistent_subclone"]
    )
    base_df["cluster_first_support_flag"] = base_df["cluster_first_support"].fillna(False)
    base_df["artifact_lowvaf_highad_flag"] = low_vaf_high_ad(base_df)
    base_df["artifact_lowvaf_highad_lowcv_flag"] = low_vaf_high_ad_low_cv(base_df)
    base_df["hp_assign_high_095"] = base_df["hp_assign_rate"] >= 0.95
    base_df["hp_assign_veryhigh_099"] = base_df["hp_assign_rate"] >= 0.99

    support_flag_names = []
    if "Quality_Score" in config.support_topbin_features:
        support_flag_names.append("quality_topbin")
    if "PairwiseMedianDist" in config.support_topbin_features:
        support_flag_names.append("pairwise_topbin")
    if "hp_assign_rate" in config.support_topbin_features:
        support_flag_names.append("hp_assign_topbin")
    support_flag_names.extend(config.support_flag_features)
    weak_support_flag_names = list(config.weak_support_flag_features)

    qc_flag_names = []
    if "Quality_Score" in config.qc_topbin_features:
        qc_flag_names.append("quality_topbin")
    if "PairwiseMedianDist" in config.qc_topbin_features:
        qc_flag_names.append("pairwise_topbin")
    if "hp_assign_rate" in config.qc_topbin_features:
        qc_flag_names.append("hp_assign_topbin")
    qc_flag_names.extend(config.qc_flag_features)

    observed_only_flag_names = []
    if "Quality_Score" in config.observed_only_topbin_features:
        observed_only_flag_names.append("quality_topbin")
    if "PairwiseMedianDist" in config.observed_only_topbin_features:
        observed_only_flag_names.append("pairwise_topbin")
    if "hp_assign_rate" in config.observed_only_topbin_features:
        observed_only_flag_names.append("hp_assign_topbin")

    primary_support_masks = [annotation_mask(base_df, name) for name in support_flag_names]
    weak_support_masks = [annotation_mask(base_df, name) for name in weak_support_flag_names]
    qc_masks = [annotation_mask(base_df, name) for name in qc_flag_names]

    if primary_support_masks:
        primary_any = primary_support_masks[0].copy()
        for mask in primary_support_masks[1:]:
            primary_any |= mask
    else:
        primary_any = pd.Series(False, index=base_df.index)

    if weak_support_masks:
        weak_any = weak_support_masks[0].copy()
        for mask in weak_support_masks[1:]:
            weak_any |= mask
    else:
        weak_any = pd.Series(False, index=base_df.index)

    if qc_masks:
        qc_any = qc_masks[0].copy()
        for mask in qc_masks[1:]:
            qc_any |= mask
    else:
        qc_any = pd.Series(False, index=base_df.index)

    base_df["support_score_primary"] = 0
    for name in support_flag_names:
        base_df["support_score_primary"] += annotation_mask(base_df, name).astype(int)
    base_df["support_score_all"] = base_df["support_score_primary"]
    for name in weak_support_flag_names:
        base_df["support_score_all"] += annotation_mask(base_df, name).astype(int)
    base_df["qc_score"] = 0
    for name in qc_flag_names:
        base_df["qc_score"] += annotation_mask(base_df, name).astype(int)

    base_df["annotation_support_primary_any"] = primary_any
    base_df["annotation_support_any"] = primary_any | weak_any
    base_df["annotation_qc_any"] = qc_any

    string_rows = base_df.apply(
        lambda row: build_annotation_strings(row, support_flag_names + weak_support_flag_names, qc_flag_names, observed_only_flag_names),
        axis=1,
        result_type="expand",
    )
    string_rows.columns = [
        "annotation_topbin_summary",
        "annotation_support_summary",
        "annotation_qc_summary",
        "annotation_role_summary",
    ]
    base_df = pd.concat([base_df, string_rows], axis=1)

    analyzed_df = candidate_eval_frame(joined_df)
    analyzed_df = analyzed_df[analyzed_df["analysis_available"]].copy()
    analyzed_df = analyzed_df.merge(
        base_df[
            [
                "region_key",
                "caller_primary_flag",
                "caller_strict_flag",
                "quality_topbin_flag",
                "pairwise_topbin_flag",
                "hp_assign_topbin_flag",
                "agreement_positive_flag",
                "cluster_first_support_flag",
                "artifact_lowvaf_highad_flag",
                "artifact_lowvaf_highad_lowcv_flag",
                "hp_assign_high_095",
                "hp_assign_veryhigh_099",
                "support_score_primary",
                "support_score_all",
                "qc_score",
                "annotation_support_primary_any",
                "annotation_support_any",
                "annotation_qc_any",
            ]
        ],
        on="region_key",
        how="left",
    )

    summary_rows: List[Dict[str, object]] = []
    summary_targets = [
        ("caller_primary_flag", "caller_primary", "support"),
        ("caller_strict_flag", "caller_strict", "support"),
        ("quality_topbin_flag", "quality_topbin", "support" if "Quality_Score" in config.support_topbin_features else "annotation"),
        ("pairwise_topbin_flag", "pairwise_topbin", "support" if "PairwiseMedianDist" in config.support_topbin_features else "annotation"),
        ("hp_assign_topbin_flag", "hp_assign_topbin", "qc" if "hp_assign_rate" in config.qc_topbin_features or "hp_assign_rate" in config.observed_only_topbin_features else "annotation"),
        ("agreement_positive_flag", "agreement_positive", "support" if "agreement_positive" in config.support_flag_features else "annotation"),
        ("cluster_first_support_flag", "cluster_first_support", "support" if "cluster_first_support" in config.weak_support_flag_features else "annotation"),
        ("artifact_lowvaf_highad_flag", "artifact_lowvaf_highad", "artifact"),
        ("artifact_lowvaf_highad_lowcv_flag", "artifact_lowvaf_highad_lowcv", "artifact"),
        ("annotation_support_primary_any", "annotation_support_primary_any", "support"),
        ("annotation_support_any", "annotation_support_any", "support"),
        ("annotation_qc_any", "annotation_qc_any", "qc"),
    ]
    total_tp = int((analyzed_df["downstream_status"] == "caller_lost_tp").sum())
    total_fp = int((analyzed_df["downstream_status"] == "caller_removed_fp").sum())
    for field, name, role in summary_targets:
        mask = analyzed_df[field].fillna(False).astype(bool)
        part = analyzed_df[mask]
        tp_count = int((part["downstream_status"] == "caller_lost_tp").sum())
        fp_count = int((part["downstream_status"] == "caller_removed_fp").sum())
        tp_fraction = tp_count / total_tp if total_tp else 0.0
        fp_fraction = fp_count / total_fp if total_fp else 0.0
        summary_rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "label": dataset.label,
                "annotation_name": name,
                "annotation_role": role,
                "trigger_count": int(mask.sum()),
                "tp_count": tp_count,
                "fp_count": fp_count,
                "tp_fraction": format_float(tp_fraction),
                "fp_fraction": format_float(fp_fraction),
                "tp_to_fp_ratio": format_float(tp_count / fp_count if fp_count else (float("inf") if tp_count else 0.0)),
                "enrichment_ratio": format_float(tp_fraction / fp_fraction if fp_fraction else (float("inf") if tp_fraction else 0.0)),
            }
        )
    summary_df = pd.DataFrame(summary_rows)

    def policy_mask(policy_id: str) -> pd.Series:
        if policy_id == "caller_primary_only":
            return analyzed_df["caller_primary_flag"]
        if policy_id == "caller_strict_only":
            return analyzed_df["caller_strict_flag"]
        if policy_id == "strict_or_support_primary":
            return analyzed_df["caller_strict_flag"] | (
                analyzed_df["caller_primary_flag"] & analyzed_df["annotation_support_primary_any"]
            )
        if policy_id == "strict_or_support_any":
            return analyzed_df["caller_strict_flag"] | (
                analyzed_df["caller_primary_flag"] & analyzed_df["annotation_support_any"]
            )
        if policy_id == "strict_or_support_score_ge2":
            return analyzed_df["caller_strict_flag"] | (
                analyzed_df["caller_primary_flag"] & (analyzed_df["support_score_all"] >= 2)
            )
        if policy_id == "strict_or_support_any_not_artifact":
            return (analyzed_df["caller_strict_flag"] | (analyzed_df["caller_primary_flag"] & analyzed_df["annotation_support_any"])) & (
                ~analyzed_df["artifact_lowvaf_highad_lowcv_flag"]
            )
        if policy_id == "primary_plus_qc":
            return analyzed_df["caller_primary_flag"] & analyzed_df["annotation_qc_any"]
        raise KeyError(policy_id)

    policy_notes = {
        "caller_primary_only": f"gq>={config.caller_primary_gq}",
        "caller_strict_only": f"gq>={config.caller_strict_gq}",
        "strict_or_support_primary": f"gq>={config.caller_strict_gq} OR (gq>={config.caller_primary_gq} and primary support annotation)",
        "strict_or_support_any": f"gq>={config.caller_strict_gq} OR (gq>={config.caller_primary_gq} and any support annotation)",
        "strict_or_support_score_ge2": f"gq>={config.caller_strict_gq} OR (gq>={config.caller_primary_gq} and support_score>=2)",
        "strict_or_support_any_not_artifact": f"strict_or_support_any and not lowVAF/highAD/lowCV artifact",
        "primary_plus_qc": f"gq>={config.caller_primary_gq} and QC annotation only",
    }
    primary_metrics = evaluate_mask(
        analyzed_df,
        analyzed_df["caller_primary_flag"],
        dataset.baseline_tp,
        dataset.baseline_fp,
        dataset.truth_total,
    )
    policy_rows: List[Dict[str, object]] = []
    for policy_id, notes in policy_notes.items():
        mask = policy_mask(policy_id)
        metrics = evaluate_mask(analyzed_df, mask, dataset.baseline_tp, dataset.baseline_fp, dataset.truth_total)
        policy_rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "label": dataset.label,
                "policy_id": policy_id,
                "policy_notes": notes,
                "trigger_count": int(mask.sum()),
                "tp_rescued": int(metrics["tp_rescued"]),
                "fp_reintroduced": int(metrics["fp_reintroduced"]),
                "precision": format_float(float(metrics["precision"])),
                "recall": format_float(float(metrics["recall"])),
                "f1": format_float(float(metrics["f1"])),
                "delta_f1_vs_baseline": format_float(float(metrics["delta_f1_vs_baseline"])),
                "delta_f1_vs_primary_gate": format_float(float(metrics["f1"]) - float(primary_metrics["f1"])),
            }
        )
    policy_df = pd.DataFrame(policy_rows)

    overlap_rows: List[Dict[str, object]] = []
    overlap_features = [
        "quality_topbin_flag",
        "pairwise_topbin_flag",
        "agreement_positive_flag",
        "cluster_first_support_flag",
        "hp_assign_high_095",
        "artifact_lowvaf_highad_lowcv_flag",
    ]
    for i, feature_a in enumerate(overlap_features):
        for feature_b in overlap_features[i + 1 :]:
            mask_a = analyzed_df[feature_a].fillna(False).astype(bool)
            mask_b = analyzed_df[feature_b].fillna(False).astype(bool)
            inter = mask_a & mask_b
            union = mask_a | mask_b
            overlap_rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "label": dataset.label,
                    "feature_a": feature_a,
                    "feature_b": feature_b,
                    "jaccard": format_float(inter.sum() / union.sum() if int(union.sum()) else 0.0),
                    "a_only_tp": int(((analyzed_df["downstream_status"] == "caller_lost_tp") & (mask_a & ~mask_b)).sum()),
                    "a_only_fp": int(((analyzed_df["downstream_status"] == "caller_removed_fp") & (mask_a & ~mask_b)).sum()),
                    "b_only_tp": int(((analyzed_df["downstream_status"] == "caller_lost_tp") & (mask_b & ~mask_a)).sum()),
                    "b_only_fp": int(((analyzed_df["downstream_status"] == "caller_removed_fp") & (mask_b & ~mask_a)).sum()),
                    "both_tp": int(((analyzed_df["downstream_status"] == "caller_lost_tp") & inter).sum()),
                    "both_fp": int(((analyzed_df["downstream_status"] == "caller_removed_fp") & inter).sum()),
                }
            )
    overlap_df = pd.DataFrame(overlap_rows)

    return base_df, summary_df, policy_df, overlap_df


def write_dataset_annotations(output_dir: Path, annotated_df: pd.DataFrame) -> None:
    if annotated_df.empty:
        return
    dataset_id = str(annotated_df["dataset_id"].iloc[0])
    path = output_dir / "annotated_candidates" / f"{dataset_id}_annotated_candidates.tsv"
    rows = annotated_df.copy()
    write_tsv_rows(path, ANNOTATED_FIELDS, rows[ANNOTATED_FIELDS].to_dict(orient="records"))


def write_markdown(
    output_dir: Path,
    config_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    policy_df: pd.DataFrame,
    overlap_df: pd.DataFrame,
    best_policy_df: pd.DataFrame,
) -> None:
    report_path = output_dir / "annotation_layer_summary.md"
    lines = [
        "# Phase 2 annotation layer summary",
        "",
        "## 1. Dataset annotation config",
        "",
        markdown_table(ANNOTATION_CONFIG_FIELDS, config_df.to_dict(orient="records")),
        "",
        "## 2. Annotation enrichment summary",
        "",
        markdown_table(
            SUMMARY_FIELDS,
            summary_df.sort_values(["dataset_id", "annotation_role", "annotation_name"]).to_dict(orient="records"),
        ),
        "",
        "解讀重點：",
        "- `annotation_role=support` 代表目前允許回接到 analysis-layer 第二層 support 的 annotation。",
        "- `annotation_role=qc` 代表 phase/QC 訊號，只建議做註記，不建議直接當 keep 規則。",
        "- `annotation_role=artifact` 代表 artifact triage 訊號，預設只用於後段否決。",
        "",
        "## 3. Annotation policy evaluation",
        "",
        markdown_table(
            POLICY_FIELDS,
            policy_df.sort_values(["dataset_id", "delta_f1_vs_primary_gate"], ascending=[True, False]).to_dict(orient="records"),
        ),
        "",
        "解讀重點：",
        "- `caller_primary_only` 是 phase 2 回接前的 caller-first 基線。",
        "- `delta_f1_vs_primary_gate` > 0 才表示 annotation tier 比原 caller 基線更好。",
        "",
        "## 4. Annotation overlap",
        "",
        markdown_table(
            OVERLAP_FIELDS,
            overlap_df.sort_values(["dataset_id", "jaccard"], ascending=[True, False]).head(40).to_dict(orient="records"),
        ),
        "",
        "## 5. Best policy summary",
        "",
        markdown_table(
            BEST_POLICY_FIELDS,
            best_policy_df.to_dict(orient="records"),
        ),
    ]
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    phase2_dir = Path(args.phase2_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    top_bins = load_top_bins(phase2_dir)
    config_rows = []
    for config in ANNOTATION_CONFIG.values():
        config_rows.append(
            {
                "dataset_id": config.dataset_id,
                "verification_support_mode": args.verification_support_mode,
                "caller_primary_gq": config.caller_primary_gq,
                "caller_strict_gq": config.caller_strict_gq,
                "support_topbin_features": ",".join(config.support_topbin_features) or "none",
                "support_flag_features": ",".join(config.support_flag_features) or "none",
                "weak_support_flag_features": ",".join(config.weak_support_flag_features) or "none",
                "qc_topbin_features": ",".join(config.qc_topbin_features) or "none",
                "qc_flag_features": ",".join(config.qc_flag_features) or "none",
                "observed_only_topbin_features": ",".join(config.observed_only_topbin_features) or "none",
                "notes": config.notes,
            }
        )
    config_df = pd.DataFrame(config_rows)

    all_summary = []
    all_policy = []
    all_overlap = []
    annotated_outputs = []
    for dataset in RESCUE_DATASETS:
        annotated_df, summary_df, policy_df, overlap_df = annotate_dataset(
            dataset,
            top_bins,
            args.verification_support_mode,
        )
        annotated_outputs.append(annotated_df)
        all_summary.append(summary_df)
        all_policy.append(policy_df)
        all_overlap.append(overlap_df)

    ensure_dir(output_dir)
    for annotated_df in annotated_outputs:
        write_dataset_annotations(output_dir, annotated_df)

    summary_df = pd.concat(all_summary, ignore_index=True)
    policy_df = pd.concat(all_policy, ignore_index=True)
    overlap_df = pd.concat(all_overlap, ignore_index=True)
    best_policy_rows = []
    for dataset_id, part in policy_df.groupby("dataset_id", sort=False):
        best_baseline = part.sort_values("delta_f1_vs_baseline", ascending=False).iloc[0]
        best_vs_gate = part.sort_values("delta_f1_vs_primary_gate", ascending=False).iloc[0]
        best_policy_rows.append(
            {
                "dataset_id": dataset_id,
                "label": best_baseline["label"],
                "best_policy_by_baseline": best_baseline["policy_id"],
                "best_delta_f1_vs_baseline": best_baseline["delta_f1_vs_baseline"],
                "best_policy_vs_primary_gate": best_vs_gate["policy_id"],
                "best_delta_f1_vs_primary_gate": best_vs_gate["delta_f1_vs_primary_gate"],
            }
        )
    best_policy_df = pd.DataFrame(best_policy_rows)

    write_tsv_rows(output_dir / "annotation_layer_config.tsv", ANNOTATION_CONFIG_FIELDS, config_df.to_dict(orient="records"))
    write_tsv_rows(output_dir / "annotation_presence_summary.tsv", SUMMARY_FIELDS, summary_df.to_dict(orient="records"))
    write_tsv_rows(output_dir / "annotation_policy_evaluation.tsv", POLICY_FIELDS, policy_df.to_dict(orient="records"))
    write_tsv_rows(output_dir / "annotation_overlap_summary.tsv", OVERLAP_FIELDS, overlap_df.to_dict(orient="records"))
    write_tsv_rows(output_dir / "annotation_best_policy.tsv", BEST_POLICY_FIELDS, best_policy_df.to_dict(orient="records"))
    write_markdown(output_dir, config_df, summary_df, policy_df, overlap_df, best_policy_df)


if __name__ == "__main__":
    main()
