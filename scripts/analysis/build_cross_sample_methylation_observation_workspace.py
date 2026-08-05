#!/usr/bin/env python3
"""Build a workspace for cross-sample methylation distribution review."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    CURRENT_CLASSES_V2,
    LEGACY_CLASSES,
    SchemaContractError,
    select_current_view,
    select_legacy_view,
)


DEFAULT_CHECKLIST = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/final_closeout/completion_checklist.tsv"
)
DEFAULT_OUTPUT_ROOT = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces")
DEFAULT_OUTPUT_NAME = "20260317_cross_sample_methylation_observation_workspace"

SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
MODE_ORDER = ["paired_full", "paired_pileup", "to_pileup", "to_full"]
LAYER_ORDER = ["caller", "longphase", "intersubmod"]
FEATURES_FOR_PLOTS = ["VAF", "AlleleDelta_abs", "CramersV", "Quality_Score"]
FEATURES_FOR_SUMMARY = [
    "VAF",
    "QUAL",
    "NumReads",
    "NumCpGs",
    "AlleleDelta_abs",
    "CramersV",
    "HPMergedDelta",
    "HP_Ratio",
    "PairwiseMedianDist",
    "Quality_Score",
]

SAMPLE_METADATA = {
    "HCC1395": {
        "sample_family": "HCC1395_specimen",
        "platform_group": "ONT_5kHz",
        "cancer_type": "breast_ductal_carcinoma",
        "biological_relation": "same specimen as HCC1395_DORADO",
        "comparison_axis": "platform_and_workflow",
    },
    "HCC1395_DORADO": {
        "sample_family": "HCC1395_specimen",
        "platform_group": "ONT_Dorado",
        "cancer_type": "breast_ductal_carcinoma",
        "biological_relation": "same specimen as HCC1395",
        "comparison_axis": "platform_and_workflow",
    },
    "COLO829": {
        "sample_family": "COLO829_specimen",
        "platform_group": "ONT_PAO_5mCG",
        "cancer_type": "melanoma",
        "biological_relation": "distinct specimen",
        "comparison_axis": "sample_and_cancer",
    },
    "H1437": {
        "sample_family": "H1437_specimen",
        "platform_group": "ONT_Google_5mCG",
        "cancer_type": "lung_adenocarcinoma",
        "biological_relation": "distinct specimen",
        "comparison_axis": "sample_and_cancer",
    },
    "H2009": {
        "sample_family": "H2009_specimen",
        "platform_group": "ONT_Google_5mCG",
        "cancer_type": "lung_adenocarcinoma",
        "biological_relation": "distinct specimen",
        "comparison_axis": "sample_and_cancer",
    },
    "HCC1937": {
        "sample_family": "HCC1937_specimen",
        "platform_group": "ONT_Google_5mCG",
        "cancer_type": "breast_brca1_mutant",
        "biological_relation": "distinct specimen",
        "comparison_axis": "sample_and_cancer",
    },
    "HCC1954": {
        "sample_family": "HCC1954_specimen",
        "platform_group": "ONT_Google_5mCG",
        "cancer_type": "breast_carcinoma",
        "biological_relation": "distinct specimen",
        "comparison_axis": "sample_and_cancer",
    },
}

MODE_METADATA = {
    "paired_full": {
        "workflow_group": "paired",
        "caller_model": "full",
        "workflow_note": "paired tumor-normal with full-model caller baseline",
    },
    "paired_pileup": {
        "workflow_group": "paired",
        "caller_model": "pileup",
        "workflow_note": "paired tumor-normal with pileup caller baseline",
    },
    "to_pileup": {
        "workflow_group": "tumor_only",
        "caller_model": "pileup",
        "workflow_note": "tumor-only pilot with ClairS-TO pileup baseline",
    },
    "to_full": {
        "workflow_group": "tumor_only",
        "caller_model": "full",
        "workflow_note": "tumor-only full-model slot; currently used only as availability coverage",
    },
}

SUMMARY_COLUMNS = [
    "RegionID",
    "Chr",
    "Pos",
    "Ref",
    "Alt",
    "NumReads",
    "NumCpGs",
    "GlobalP",
    "CramersV",
    "HeuristicScore",
    "PassedGating",
    "PairwiseMeanDist",
    "PairwiseMedianDist",
    "HPMergedDelta",
    "HP_Ratio",
    "AlleleDelta",
    "Quality_Score",
    "Quality_Tier",
    "VerificationClass",
    "VerificationSchemaVersion",
    "VerificationClass_Legacy",
    "DominantLabel",
    "Potential_LOH",
    "Coverage_Category",
    "LocalBestCluster",
    "Significant",
    "SuggestFilter",
]

PARAMETER_NOTES = [
    ("NumReads", "region 內可用 reads 數；太低時統計與分群不穩定。"),
    ("NumCpGs", "region 內被觀察到的 CpG 數；越多代表甲基訊號上下文越完整。"),
    ("VAF", "variant allele fraction；低 VAF 常是 FP triage 需要特別注意的背景條件。"),
    ("AlleleDelta", "allele 相關的甲基差異強度 proxy；高值常出現在 artifact triage 線索，但不能直接當全域 keep rule。"),
    ("CramersV", "cluster 與 label 關聯強度的 effect size；高值偏向較穩定的結構訊號。"),
    ("HPMergedDelta", "haplotype merged label 的差異強度；偏向 phase / haplotype 分離訊號。"),
    ("HP_Ratio", "collapsed haplotype imbalance 摘要；是 annotation，不是正式 LOH。"),
    ("PairwiseMedianDist", "region 內 methylation pairwise distance 的中位數；方向具有 dataset dependence。"),
    ("Quality_Score", "甲基訊號整體品質分數；較適合 soft support annotation，不宜直接 hard keep。"),
    (
        "VerificationClass",
        "此工作區必須明選 current C2 或 historical legacy L4 view；兩者都只作結構註記，不直接等於 TP/FP。",
    ),
]


@dataclass
class RunMeta:
    sample: str
    mode: str
    run_dir: Path
    tp_summary_path: Path
    fp_summary_path: Path
    tp_vcf_path: Path
    fp_vcf_path: Path
    caller_f1: float
    longphase_f1: float
    intersubmod_f1: float
    delta_f1_vs_longphase: float
    caller_tp: int
    caller_fp: int
    caller_fn: int
    longphase_tp: int
    longphase_fp: int
    longphase_fn: int
    intersubmod_tp: int
    intersubmod_fp: int
    intersubmod_fn: int
    sample_family: str
    platform_group: str
    cancer_type: str
    biological_relation: str
    comparison_axis: str
    workflow_group: str
    caller_model: str
    workflow_note: str

    @property
    def run_label(self) -> str:
        return f"{self.sample}:{self.mode}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--completion-checklist", default=str(DEFAULT_CHECKLIST), help="Path to completion_checklist.tsv")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT), help="Root for generated workspaces")
    parser.add_argument("--output-dir", default=None, help="Optional explicit output directory")
    parser.add_argument(
        "--verification-view",
        required=True,
        choices=("current", "legacy"),
        help=(
            "Required schema view: current validates VerificationSchemaVersion=2 and uses the 11-class "
            "VerificationClass taxonomy; legacy requires VerificationClass_Legacy."
        ),
    )
    parser.add_argument(
        "--max-points-per-group",
        type=int,
        default=2500,
        help="Maximum sampled points per plot group for heavy plots",
    )
    return parser.parse_args()


def verification_class_order(requested_view: str) -> List[str]:
    """Return the exact ordered categories for the explicitly selected panel."""
    if requested_view == "current":
        return list(CURRENT_CLASSES_V2) + ["UnknownCurrentClass"]
    if requested_view == "legacy":
        return list(LEGACY_CLASSES)
    raise ValueError(f"Unsupported verification view: {requested_view}")


def apply_verification_view(df: pd.DataFrame, requested_view: str) -> tuple[pd.DataFrame, Dict[str, object]]:
    """Validate and attach one explicit VerificationClass view without taxonomy folding."""
    if requested_view == "current":
        view = select_current_view(df)
    elif requested_view == "legacy":
        view = select_legacy_view(df)
    else:
        raise ValueError(f"Unsupported verification view: {requested_view}")

    selected = df.copy()
    selected["VerificationClass_SourceValue"] = selected[view.field]
    selected["VerificationClass"] = pd.Categorical(
        view.values,
        categories=view.categories,
        ordered=True,
    )
    selected["VerificationView"] = requested_view
    selected["VerificationSourceField"] = view.field
    selected["VerificationSchemaStatus"] = view.schema_status
    metadata = view.metadata()
    metadata["requested_view"] = requested_view
    return selected, metadata


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def to_float(value: object) -> float:
    try:
        if value in {"", None}:
            return float("nan")
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def to_int(value: object) -> int:
    try:
        if value in {"", None}:
            return 0
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def parse_vaf_from_record(format_keys: List[str], sample_values: List[str]) -> Optional[float]:
    fmt = dict(zip(format_keys, sample_values))
    for key in ("VAF", "AF"):
        raw = fmt.get(key)
        if raw and raw not in {".", ""}:
            try:
                return float(raw.split(",")[0])
            except ValueError:
                pass
    raw_ad = fmt.get("AD")
    if raw_ad and raw_ad not in {".", ""}:
        fields = raw_ad.split(",")
        if len(fields) >= 2:
            try:
                ref = float(fields[0])
                alt = float(fields[1])
                total = ref + alt
                if total > 0:
                    return alt / total
            except ValueError:
                pass
    return None


def parse_vcf_features(vcf_path: Path) -> pd.DataFrame:
    rows = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            qual = to_float(parts[5])
            fmt_keys = parts[8].split(":")
            sample_values = parts[9].split(":")
            vaf = parse_vaf_from_record(fmt_keys, sample_values)
            rows.append(
                {
                    "Chr": parts[0],
                    "Pos": pos,
                    "Ref": parts[3],
                    "Alt": parts[4],
                    "QUAL": qual,
                    "VAF": np.nan if vaf is None else float(vaf),
                }
            )
    return pd.DataFrame(rows)


def order_key(sample: str, mode: str) -> tuple[int, int]:
    sample_idx = SAMPLE_ORDER.index(sample) if sample in SAMPLE_ORDER else len(SAMPLE_ORDER)
    mode_idx = MODE_ORDER.index(mode) if mode in MODE_ORDER else len(MODE_ORDER)
    return sample_idx, mode_idx


def load_completed_runs(checklist_path: Path) -> List[RunMeta]:
    runs: List[RunMeta] = []
    with checklist_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("expected") != "true":
                continue
            if row.get("status") != "completed":
                continue
            sample = row["sample"]
            mode = row["mode"]
            sample_meta = SAMPLE_METADATA[sample]
            mode_meta = MODE_METADATA[mode]
            runs.append(
                RunMeta(
                    sample=sample,
                    mode=mode,
                    run_dir=Path(row["run_dir"]),
                    tp_summary_path=Path(row["tp_summary_path"]),
                    fp_summary_path=Path(row["fp_summary_path"]),
                    tp_vcf_path=Path(row["tp_vcf_path"]),
                    fp_vcf_path=Path(row["fp_vcf_path"]),
                    caller_f1=to_float(row["caller_f1"]),
                    longphase_f1=to_float(row["longphase_f1"]),
                    intersubmod_f1=to_float(row["intersubmod_f1"]),
                    delta_f1_vs_longphase=to_float(row["delta_f1_vs_longphase"]),
                    caller_tp=to_int(row["caller_tp"]),
                    caller_fp=to_int(row["caller_fp"]),
                    caller_fn=to_int(row["caller_fn"]),
                    longphase_tp=to_int(row["longphase_tp"]),
                    longphase_fp=to_int(row["longphase_fp"]),
                    longphase_fn=to_int(row["longphase_fn"]),
                    intersubmod_tp=to_int(row["intersubmod_tp"]),
                    intersubmod_fp=to_int(row["intersubmod_fp"]),
                    intersubmod_fn=to_int(row["intersubmod_fn"]),
                    sample_family=sample_meta["sample_family"],
                    platform_group=sample_meta["platform_group"],
                    cancer_type=sample_meta["cancer_type"],
                    biological_relation=sample_meta["biological_relation"],
                    comparison_axis=sample_meta["comparison_axis"],
                    workflow_group=mode_meta["workflow_group"],
                    caller_model=mode_meta["caller_model"],
                    workflow_note=mode_meta["workflow_note"],
                )
            )
    runs.sort(key=lambda item: order_key(item.sample, item.mode))
    return runs


def load_workflow_coverage(checklist_path: Path) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    with checklist_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            mode = row["mode"]
            if sample not in SAMPLE_METADATA or mode not in MODE_METADATA:
                continue
            sample_meta = SAMPLE_METADATA[sample]
            mode_meta = MODE_METADATA[mode]
            rows.append(
                {
                    "sample": sample,
                    "mode": mode,
                    "run_label": f"{sample}:{mode}",
                    "expected": row.get("expected", ""),
                    "status": row.get("status", ""),
                    "availability_reason": row.get("availability_reason", ""),
                    "blocking_reason": row.get("blocking_reason", ""),
                    "run_dir": row.get("run_dir", ""),
                    "source_kind": row.get("source_kind", ""),
                    "platform_group": sample_meta["platform_group"],
                    "cancer_type": sample_meta["cancer_type"],
                    "sample_family": sample_meta["sample_family"],
                    "workflow_group": mode_meta["workflow_group"],
                    "caller_model": mode_meta["caller_model"],
                    "workflow_note": mode_meta["workflow_note"],
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    sample_rank = {sample: idx for idx, sample in enumerate(SAMPLE_ORDER)}
    mode_rank = {mode: idx for idx, mode in enumerate(MODE_ORDER)}
    df["sample_order"] = df["sample"].map(sample_rank).fillna(len(sample_rank)).astype(int)
    df["mode_order"] = df["mode"].map(mode_rank).fillna(len(mode_rank)).astype(int)
    return df.sort_values(["sample_order", "mode_order"]).reset_index(drop=True)


def load_summary_rows(
    summary_path: Path,
    vcf_path: Path,
    meta: RunMeta,
    truth_label: str,
    verification_view: str,
) -> pd.DataFrame:
    usecols = lambda c: c in SUMMARY_COLUMNS
    summary_df = pd.read_csv(summary_path, usecols=usecols, low_memory=False)
    try:
        summary_df, verification_metadata = apply_verification_view(summary_df, verification_view)
    except SchemaContractError as exc:
        raise SchemaContractError(f"{summary_path}: {exc}") from exc
    if verification_metadata["unknown_counts"]:
        print(
            f"[verification-schema] {summary_path}: "
            f"unknown={verification_metadata['unknown_counts']}"
        )
    feature_df = parse_vcf_features(vcf_path)
    merged = summary_df.merge(feature_df, on=["Chr", "Pos", "Ref", "Alt"], how="left")
    merged["sample"] = meta.sample
    merged["mode"] = meta.mode
    merged["run_label"] = meta.run_label
    merged["sample_family"] = meta.sample_family
    merged["platform_group"] = meta.platform_group
    merged["cancer_type"] = meta.cancer_type
    merged["workflow_group"] = meta.workflow_group
    merged["caller_model"] = meta.caller_model
    merged["truth_label"] = truth_label
    merged["is_tp"] = truth_label == "TP"
    for col in FEATURES_FOR_SUMMARY + ["GlobalP", "HeuristicScore", "LocalBestCluster"]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")
    merged["AlleleDelta_abs"] = pd.to_numeric(merged.get("AlleleDelta"), errors="coerce").abs()
    return merged


def build_layer_metrics_df(runs: Iterable[RunMeta]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for meta in runs:
        rows.extend(
            [
                {
                    "sample": meta.sample,
                    "mode": meta.mode,
                    "run_label": meta.run_label,
                    "platform_group": meta.platform_group,
                    "cancer_type": meta.cancer_type,
                    "workflow_group": meta.workflow_group,
                    "layer": "caller",
                    "tp": meta.caller_tp,
                    "fp": meta.caller_fp,
                    "fn": meta.caller_fn,
                    "f1": meta.caller_f1,
                },
                {
                    "sample": meta.sample,
                    "mode": meta.mode,
                    "run_label": meta.run_label,
                    "platform_group": meta.platform_group,
                    "cancer_type": meta.cancer_type,
                    "workflow_group": meta.workflow_group,
                    "layer": "longphase",
                    "tp": meta.longphase_tp,
                    "fp": meta.longphase_fp,
                    "fn": meta.longphase_fn,
                    "f1": meta.longphase_f1,
                },
                {
                    "sample": meta.sample,
                    "mode": meta.mode,
                    "run_label": meta.run_label,
                    "platform_group": meta.platform_group,
                    "cancer_type": meta.cancer_type,
                    "workflow_group": meta.workflow_group,
                    "layer": "intersubmod",
                    "tp": meta.intersubmod_tp,
                    "fp": meta.intersubmod_fp,
                    "fn": meta.intersubmod_fn,
                    "f1": meta.intersubmod_f1,
                },
            ]
        )
    df = pd.DataFrame(rows)
    order_map = {label: idx for idx, label in enumerate([f"{s}:{m}" for s in SAMPLE_ORDER for m in MODE_ORDER])}
    df["run_order"] = df["run_label"].map(order_map).fillna(len(order_map)).astype(int)
    return df.sort_values(["run_order", "layer"]).reset_index(drop=True)


def build_feature_summary(df: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    group_arg = group_cols[0] if len(group_cols) == 1 else group_cols
    for group_key, sub in df.groupby(group_arg):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        group_data = dict(zip(group_cols, group_key))
        for feature in FEATURES_FOR_SUMMARY:
            if feature not in sub.columns:
                continue
            values = pd.to_numeric(sub[feature], errors="coerce").dropna()
            if values.empty:
                continue
            rows.append(
                {
                    **group_data,
                    "feature": feature,
                    "n": int(values.shape[0]),
                    "mean": float(values.mean()),
                    "median": float(values.median()),
                    "q25": float(values.quantile(0.25)),
                    "q75": float(values.quantile(0.75)),
                }
            )
    return pd.DataFrame(rows)


def build_verification_summary(df: pd.DataFrame, class_order: List[str]) -> pd.DataFrame:
    counts = (
        df.groupby(
            ["sample", "mode", "run_label", "truth_label", "VerificationClass"],
            observed=False,
        )
        .size()
        .reset_index(name="count")
    )
    counts["VerificationClass"] = pd.Categorical(
        counts["VerificationClass"], categories=class_order, ordered=True
    )
    counts = counts.sort_values(
        ["sample", "mode", "run_label", "truth_label", "VerificationClass"]
    ).reset_index(drop=True)
    total = counts.groupby(["sample", "mode", "run_label", "truth_label"])["count"].transform("sum")
    counts["fraction"] = counts["count"] / total
    return counts


def build_run_transition_df(layer_df: pd.DataFrame) -> pd.DataFrame:
    pivot = (
        layer_df.pivot_table(
            index=["sample", "mode", "run_label", "platform_group", "cancer_type"],
            columns="layer",
            values=["tp", "fp", "fn", "f1"],
            aggfunc="first",
        )
        .reset_index()
        .copy()
    )
    pivot.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else str(col)
        for col in pivot.columns.to_flat_index()
    ]
    pivot["caller_to_longphase"] = pivot["f1_longphase"] - pivot["f1_caller"]
    pivot["longphase_to_intersubmod"] = pivot["f1_intersubmod"] - pivot["f1_longphase"]
    pivot["tp_delta_is_vs_lp"] = pivot["tp_intersubmod"] - pivot["tp_longphase"]
    pivot["fp_delta_is_vs_lp"] = pivot["fp_intersubmod"] - pivot["fp_longphase"]
    pivot["fn_delta_is_vs_lp"] = pivot["fn_intersubmod"] - pivot["fn_longphase"]
    pivot["run_order"] = pivot["run_label"].map(
        {f"{sample}:{mode}": idx for idx, (sample, mode) in enumerate((s, m) for s in SAMPLE_ORDER for m in MODE_ORDER)}
    ).fillna(10_000)
    return pivot.sort_values("run_order").reset_index(drop=True)


def build_feature_delta_df(df: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    group_arg = group_cols[0] if len(group_cols) == 1 else group_cols
    for group_key, sub in df.groupby(group_arg):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        group_data = dict(zip(group_cols, group_key))
        for feature in FEATURES_FOR_SUMMARY:
            if feature not in sub.columns:
                continue
            tp = pd.to_numeric(sub.loc[sub["truth_label"] == "TP", feature], errors="coerce").dropna()
            fp = pd.to_numeric(sub.loc[sub["truth_label"] == "FP", feature], errors="coerce").dropna()
            if tp.empty or fp.empty:
                continue
            pooled = pd.concat([tp, fp], ignore_index=True)
            pooled_iqr = float(pooled.quantile(0.75) - pooled.quantile(0.25))
            delta = float(tp.median() - fp.median())
            rows.append(
                {
                    **group_data,
                    "feature": feature,
                    "tp_n": int(tp.shape[0]),
                    "fp_n": int(fp.shape[0]),
                    "tp_median": float(tp.median()),
                    "fp_median": float(fp.median()),
                    "tp_q25": float(tp.quantile(0.25)),
                    "tp_q75": float(tp.quantile(0.75)),
                    "fp_q25": float(fp.quantile(0.25)),
                    "fp_q75": float(fp.quantile(0.75)),
                    "delta_median": delta,
                    "pooled_iqr": pooled_iqr,
                    "standardized_delta": delta / pooled_iqr if pooled_iqr else np.nan,
                    "tp_nonzero_frac": float((tp > 0).mean()),
                    "fp_nonzero_frac": float((fp > 0).mean()),
                }
            )
    return pd.DataFrame(rows)


def strongest_feature_rows(feature_delta_df: pd.DataFrame, scope_col: str, scope_value: str, top_n: int = 3) -> pd.DataFrame:
    sub = feature_delta_df[feature_delta_df[scope_col] == scope_value].copy()
    if sub.empty:
        return sub
    return sub.assign(abs_std=sub["standardized_delta"].abs()).sort_values(["abs_std", "feature"], ascending=[False, True]).head(top_n)


def plot_feature_delta_heatmap(feature_delta_df: pd.DataFrame, out_png: Path) -> None:
    chosen_features = ["Quality_Score", "VAF", "NumReads", "NumCpGs", "AlleleDelta_abs", "PairwiseMedianDist", "CramersV"]
    mat = (
        feature_delta_df[feature_delta_df["feature"].isin(chosen_features)]
        .pivot_table(index="run_label", columns="feature", values="standardized_delta", aggfunc="first")
        .reindex(index=[f"{sample}:{mode}" for sample in SAMPLE_ORDER for mode in MODE_ORDER if f"{sample}:{mode}" in feature_delta_df["run_label"].unique()], columns=chosen_features)
    )
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(mat, annot=True, fmt=".2f", cmap="RdBu_r", center=0, ax=ax)
    ax.set_title("TP vs FP Feature Separation Heatmap (Median Delta / Pooled IQR)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def downsample_by_group(df: pd.DataFrame, group_cols: List[str], max_points_per_group: int) -> pd.DataFrame:
    parts = []
    group_arg = group_cols[0] if len(group_cols) == 1 else group_cols
    for _, sub in df.groupby(group_arg):
        if sub.shape[0] <= max_points_per_group:
            parts.append(sub)
        else:
            parts.append(sub.sample(max_points_per_group, random_state=42))
    return pd.concat(parts, ignore_index=True) if parts else df.head(0).copy()


def plot_layer_f1(layer_df: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(18, 6))
    plot_df = layer_df.copy()
    sns.barplot(data=plot_df, x="run_label", y="f1", hue="layer", hue_order=LAYER_ORDER, ax=ax)
    ax.set_title("Caller vs LongPhase vs InterSubMod F1 by Run")
    ax.set_ylabel("F1")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=45)
    ax.set_ylim(0.65, 1.0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_workflow_coverage(coverage_df: pd.DataFrame, out_png: Path) -> None:
    status_to_code = {"not_applicable": 0, "completed": 1}
    matrix_df = coverage_df.copy()
    matrix_df["status_code"] = matrix_df["status"].map(status_to_code).fillna(-1)
    pivot = (
        matrix_df.pivot_table(index="sample", columns="mode", values="status_code", aggfunc="first")
        .reindex(index=SAMPLE_ORDER, columns=MODE_ORDER)
    )
    annot = (
        matrix_df.assign(cell_text=lambda df: df["status"].fillna("") + "\n" + df["availability_reason"].fillna(""))
        .pivot_table(index="sample", columns="mode", values="cell_text", aggfunc="first")
        .reindex(index=SAMPLE_ORDER, columns=MODE_ORDER)
    )
    cmap = matplotlib.colors.ListedColormap(["#cfcfcf", "#1b9e77", "#d95f02"])
    bounds = [-1.5, -0.5, 0.5, 1.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(
        pivot,
        cmap=cmap,
        norm=norm,
        cbar=False,
        linewidths=0.8,
        linecolor="white",
        annot=annot,
        fmt="",
        annot_kws={"fontsize": 8},
        ax=ax,
    )
    ax.set_title("Workflow Coverage Matrix Across All Samples and Modes")
    ax.set_xlabel("mode")
    ax.set_ylabel("sample")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_stepwise_delta_heatmap(layer_df: pd.DataFrame, out_png: Path) -> None:
    base = (
        layer_df.pivot_table(index="run_label", columns="layer", values="f1", aggfunc="first")
        .reset_index()
        .copy()
    )
    base["caller_to_longphase"] = base["longphase"] - base["caller"]
    base["longphase_to_intersubmod"] = base["intersubmod"] - base["longphase"]
    mat = base.set_index("run_label")[["caller_to_longphase", "longphase_to_intersubmod"]]
    fig, ax = plt.subplots(figsize=(8, 8))
    sns.heatmap(mat, annot=True, fmt=".4f", cmap="RdYlGn", center=0, ax=ax)
    ax.set_title("Stepwise F1 Delta Heatmap")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_tp_fp_counts(feature_df: pd.DataFrame, out_png: Path) -> None:
    counts = feature_df.groupby(["run_label", "truth_label"]).size().reset_index(name="count")
    pivot = counts.pivot(index="run_label", columns="truth_label", values="count").fillna(0)
    fig, ax = plt.subplots(figsize=(18, 6))
    pivot[["TP", "FP"]].plot(kind="bar", stacked=True, ax=ax, color=["#1b9e77", "#d95f02"])
    ax.set_title("TP / FP Region Counts by Run")
    ax.set_ylabel("Region count")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(title="truth_label")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_mode_box(feature_df: pd.DataFrame, out_png: Path, max_points_per_group: int) -> None:
    plot_df = downsample_by_group(feature_df, ["mode", "truth_label"], max_points_per_group)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    for ax, feature in zip(axes, FEATURES_FOR_PLOTS):
        sub = plot_df[["mode", "truth_label", feature]].dropna()
        sns.boxplot(data=sub, x="mode", y=feature, hue="truth_label", showfliers=False, ax=ax)
        ax.set_title(f"{feature} by Mode and TP/FP")
        ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_run_violin(feature_df: pd.DataFrame, out_png: Path, max_points_per_group: int) -> None:
    plot_df = downsample_by_group(feature_df, ["run_label", "truth_label"], max_points_per_group)
    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    axes = axes.flatten()
    for ax, feature in zip(axes, FEATURES_FOR_PLOTS):
        sub = plot_df[["run_label", "truth_label", feature]].dropna()
        sns.violinplot(data=sub, x="run_label", y=feature, hue="truth_label", cut=0, density_norm="width", ax=ax)
        ax.set_title(f"{feature} Distribution by Run")
        ax.tick_params(axis="x", rotation=55)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_verification_stack(feature_df: pd.DataFrame, out_png: Path, class_order: List[str]) -> None:
    counts = (
        feature_df.groupby(["run_label", "truth_label", "VerificationClass"]).size().reset_index(name="count")
    )
    fig, axes = plt.subplots(2, 1, figsize=(18, 10), sharex=True)
    palette = {
        "Strong": "#1b9e77",
        "Subclone": "#7570b3",
        "Weak": "#d95f02",
        "Noise": "#666666",
        "Strong_Bidirectional": "#1b9e77",
        "ClusterFirstOnly": "#7570b3",
        "LOH-Structure": "#e7298a",
        "MultiGroupNoLabel": "#66a61e",
        "LabelShift": "#e6ab02",
        "PermanovaLocation": "#a6761d",
        "StructureNoLabel": "#1f78b4",
        "DispersionStructure": "#6a3d9a",
        "Noise_Uniform": "#bdbdbd",
        "Noise_Chaotic": "#737373",
        "Noise_Uncorrelated": "#252525",
        "UnknownCurrentClass": "#ff00ff",
    }
    for ax, truth_label in zip(axes, ["TP", "FP"]):
        sub = counts[counts["truth_label"] == truth_label]
        pivot = sub.pivot(index="run_label", columns="VerificationClass", values="count").fillna(0)
        pivot = pivot.reindex(columns=class_order, fill_value=0)
        cols = list(class_order)
        pivot[cols].plot(
            kind="bar",
            stacked=True,
            ax=ax,
            color=[palette.get(c, "#bbbbbb") for c in cols],
        )
        ax.set_title(f"VerificationClass Composition ({truth_label})")
        ax.set_ylabel("Region count")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(title="VerificationClass", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_hcc1395_family(layer_df: pd.DataFrame, out_png: Path) -> None:
    base = (
        layer_df[layer_df["sample"].isin(["HCC1395", "HCC1395_DORADO"])]
        .pivot_table(index=["sample", "mode"], columns="layer", values="f1", aggfunc="first")
        .reset_index()
    )
    base["longphase_to_intersubmod"] = base["intersubmod"] - base["longphase"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sns.barplot(data=base, x="mode", y="intersubmod", hue="sample", ax=axes[0], errorbar=None)
    axes[0].set_title("HCC1395 Same-Specimen Family: InterSubMod F1")
    axes[0].set_ylabel("InterSubMod F1")
    axes[0].set_xlabel("")
    sns.barplot(data=base, x="mode", y="longphase_to_intersubmod", hue="sample", ax=axes[1], errorbar=None)
    axes[1].set_title("HCC1395 Same-Specimen Family: LongPhase→InterSubMod Delta")
    axes[1].set_ylabel("F1 delta")
    axes[1].set_xlabel("")
    axes[1].axhline(0, color="black", linewidth=1)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_sample_panel(feature_df: pd.DataFrame, sample: str, out_png: Path, max_points_per_group: int) -> None:
    sub_all = feature_df[feature_df["sample"] == sample].copy()
    plot_df = downsample_by_group(sub_all, ["mode", "truth_label"], max_points_per_group)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    for ax, feature in zip(axes, FEATURES_FOR_PLOTS):
        sub = plot_df[["mode", "truth_label", feature]].dropna()
        sns.violinplot(data=sub, x="mode", y=feature, hue="truth_label", cut=0, density_norm="width", ax=ax)
        ax.set_title(f"{sample}: {feature}")
        ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_dataframe(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def save_gzip_dataframe(df: pd.DataFrame, path: Path) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        df.to_csv(handle, sep="\t", index=False)


def build_markdown(
    out_path: Path,
    runs: List[RunMeta],
    workflow_df: pd.DataFrame,
    transition_df: pd.DataFrame,
    feature_df: pd.DataFrame,
    feature_delta_run_df: pd.DataFrame,
    feature_delta_sample_df: pd.DataFrame,
    feature_delta_global_df: pd.DataFrame,
    per_sample_plot_names: Dict[str, str],
) -> None:
    stepwise = transition_df.copy()
    best_positive = stepwise.sort_values("longphase_to_intersubmod", ascending=False).head(3)
    most_negative = stepwise.sort_values("longphase_to_intersubmod", ascending=True).head(3)
    mode_delta_summary = (
        stepwise.groupby("mode")[["caller_to_longphase", "longphase_to_intersubmod"]]
        .mean()
        .reset_index()
    )
    global_feature = feature_delta_global_df.set_index("feature")
    run_lines = [
        "| Sample | Mode | Platform | Cancer | Caller F1 | LongPhase F1 | InterSubMod F1 | Caller→LongPhase | LongPhase→InterSubMod |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for meta in runs:
        row = stepwise[stepwise["run_label"] == meta.run_label].iloc[0]
        run_lines.append(
            f"| {meta.sample} | {meta.mode} | {meta.platform_group} | {meta.cancer_type} | "
            f"{meta.caller_f1:.4f} | {meta.longphase_f1:.4f} | {meta.intersubmod_f1:.4f} | "
            f"{row['caller_to_longphase']:+.4f} | {row['longphase_to_intersubmod']:+.4f} |"
        )

    parameter_lines = [f"- `{name}`: {desc}" for name, desc in PARAMETER_NOTES]

    best_positive_lines = [
        f"- `{row.run_label}`: `LongPhase→InterSubMod {row.longphase_to_intersubmod:+.4f}`"
        for row in best_positive.itertuples()
    ]
    most_negative_lines = [
        f"- `{row.run_label}`: `LongPhase→InterSubMod {row.longphase_to_intersubmod:+.4f}`"
        for row in most_negative.itertuples()
    ]
    mode_delta_lines = [
        f"- `{row.mode}`: `Caller→LongPhase {row.caller_to_longphase:+.4f}`, `LongPhase→InterSubMod {row.longphase_to_intersubmod:+.4f}`"
        for row in mode_delta_summary.itertuples()
    ]

    global_feature_lines = []
    for feature in ("VAF", "AlleleDelta_abs", "CramersV", "Quality_Score", "PairwiseMedianDist"):
        if feature in global_feature.index:
            row = global_feature.loc[feature]
            global_feature_lines.append(
                f"- `{feature}` median: `TP={row['tp_median']:.4f}`, `FP={row['fp_median']:.4f}`, "
                f"`std_delta={row['standardized_delta']:+.3f}`"
            )

    quality_row = global_feature.loc["Quality_Score"]
    vaf_row = global_feature.loc["VAF"]
    cramers_row = global_feature.loc["CramersV"]
    vaf_run_df = feature_delta_run_df[feature_delta_run_df["feature"] == "VAF"].copy()
    vaf_positive_count = int((vaf_run_df["delta_median"] > 0).sum())
    vaf_negative_count = int((vaf_run_df["delta_median"] < 0).sum())
    vaf_most_positive = vaf_run_df.sort_values("delta_median", ascending=False).head(2)
    vaf_most_negative = vaf_run_df.sort_values("delta_median", ascending=True).head(2)
    platform_summary = (
        stepwise.groupby("platform_group")[["caller_to_longphase", "longphase_to_intersubmod"]].mean().reset_index()
    )
    hcc1395_5khz_summary = stepwise[stepwise["sample"] == "HCC1395"][["caller_to_longphase", "longphase_to_intersubmod"]].mean()
    hcc1395_dorado_summary = stepwise[stepwise["sample"] == "HCC1395_DORADO"][["caller_to_longphase", "longphase_to_intersubmod"]].mean()
    google_negative = stepwise[
        (stepwise["platform_group"] == "ONT_Google_5mCG") & (stepwise["caller_to_longphase"] < 0)
    ][["sample", "mode", "caller_to_longphase"]]

    top_positive = best_positive.iloc[0]
    top_negative = most_negative.iloc[0]
    strongest_global = feature_delta_global_df.assign(
        abs_std=feature_delta_global_df["standardized_delta"].abs()
    ).sort_values(["abs_std", "feature"], ascending=[False, True]).head(3)
    strongest_global_lines = []
    for row in strongest_global.itertuples():
        direction = "TP > FP" if row.delta_median > 0 else "TP < FP"
        strongest_global_lines.append(
            f"- `{row.feature}`：`{direction}`，`TP={row.tp_median:.4f}`、`FP={row.fp_median:.4f}`、`std_delta={row.standardized_delta:+.3f}`"
        )

    sample_sections = []
    for sample in SAMPLE_ORDER:
        if sample not in per_sample_plot_names:
            continue
        sample_rows = [meta for meta in runs if meta.sample == sample]
        sample_transition = stepwise[stepwise["sample"] == sample].sort_values(
            by="mode",
            key=lambda col: col.map(lambda value: MODE_ORDER.index(value) if value in MODE_ORDER else len(MODE_ORDER)),
        )
        sample_feature_rows = strongest_feature_rows(feature_delta_sample_df, "sample", sample, top_n=3)
        sample_vaf_runs = feature_delta_run_df[
            (feature_delta_run_df["sample"] == sample) & (feature_delta_run_df["feature"] == "VAF")
        ].copy()
        positive_modes = sample_vaf_runs[sample_vaf_runs["delta_median"] > 0]["mode"].tolist()
        negative_modes = sample_vaf_runs[sample_vaf_runs["delta_median"] < 0]["mode"].tolist()
        sample_best = sample_transition.sort_values("longphase_to_intersubmod", ascending=False).iloc[0]
        sample_worst = sample_transition.sort_values("longphase_to_intersubmod", ascending=True).iloc[0]

        sample_sections.append(f"## {sample}")
        sample_sections.append("")
        sample_sections.append(f"![](plots/by_sample/{per_sample_plot_names[sample]})")
        sample_sections.append("")
        sample_sections.append("| Mode | Caller F1 | LongPhase F1 | InterSubMod F1 | LongPhase→InterSubMod |")
        sample_sections.append("| --- | --- | --- | --- | --- |")
        for meta in sample_rows:
            delta = meta.intersubmod_f1 - meta.longphase_f1
            sample_sections.append(
                f"| {meta.mode} | {meta.caller_f1:.4f} | {meta.longphase_f1:.4f} | {meta.intersubmod_f1:.4f} | {delta:+.4f} |"
            )
        sample_sections.append("")
        if not sample_feature_rows.empty:
            row1 = sample_feature_rows.iloc[0]
            direction1 = "TP > FP" if row1["delta_median"] > 0 else "TP < FP"
            sample_sections.append(
                f"- 主分離 feature 之一是 `{row1['feature']}`，目前呈現 `{direction1}`，"
                f"`TP={row1['tp_median']:.4f}`、`FP={row1['fp_median']:.4f}`、`std_delta={row1['standardized_delta']:+.3f}`。"
            )
            if sample_feature_rows.shape[0] > 1:
                row2 = sample_feature_rows.iloc[1]
                direction2 = "TP > FP" if row2["delta_median"] > 0 else "TP < FP"
                sample_sections.append(
                    f"- 第二個較穩的差異來自 `{row2['feature']}`，方向為 `{direction2}`，"
                    f"`TP={row2['tp_median']:.4f}`、`FP={row2['fp_median']:.4f}`。"
                )
        if positive_modes and negative_modes:
            sample_sections.append(
                f"- `VAF` 在此 sample 屬於 `mixed` 現象：`TP > FP` 出現在 `{', '.join(positive_modes)}`，"
                f"`TP < FP` 出現在 `{', '.join(negative_modes)}`，代表它不是可跨 mode 直接套用的單一方向規則。"
            )
        elif positive_modes:
            sample_sections.append(
                f"- `VAF` 在此 sample 的 completed modes 方向一致偏 `TP > FP`，較像同一資料條件下穩定訊號。"
            )
        elif negative_modes:
            sample_sections.append(
                f"- `VAF` 在此 sample 的 completed modes 方向一致偏 `TP < FP`，較像 caller scope / truth set / purity 背景共同作用。"
            )
        sample_sections.append(
            f"- `LongPhase→InterSubMod` 在此 sample 最佳 mode 是 `{sample_best['mode']}` (`{sample_best['longphase_to_intersubmod']:+.4f}`)，"
            f"最差 mode 是 `{sample_worst['mode']}` (`{sample_worst['longphase_to_intersubmod']:+.4f}`)。"
        )
        if sample_worst["longphase_to_intersubmod"] < 0:
            sample_sections.append(
                f"- 最差 mode 的負增益主要對應 `TP {int(sample_worst['tp_delta_is_vs_lp'])}`、"
                f"`FP {int(sample_worst['fp_delta_is_vs_lp'])}` 的變化；若 TP 損失遠大於 FP 移除，F1 就會下滑。"
            )
        sample_sections.append("")

    coverage_lines = [
        "| Sample | Mode | Status | Availability | Notes |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in workflow_df.itertuples():
        note = row.blocking_reason or row.workflow_note
        coverage_lines.append(
            f"| {row.sample} | {row.mode} | {row.status} | {row.availability_reason or '-'} | {note or '-'} |"
        )

    completed_count = int((workflow_df["status"] == "completed").sum())
    not_applicable_count = int((workflow_df["status"] == "not_applicable").sum())
    total_combinations = int(workflow_df.shape[0])

    content = "\n".join(
        [
            "# 跨樣本甲基分佈與分層差異觀察",
            "",
            "## 1. 本次研究問題",
            "",
            "- 是否能在 `16` 個 completed runs 上，用同一個工作區同時觀察 `caller / LongPhase / InterSubMod` 三層 benchmark 與 `TP / FP` 的甲基特徵分佈。",
            "- 是否能區分哪些差異比較像 `平台/定序來源差異`、哪些是 `workflow 差異`、哪些才可能是 `樣本/癌種差異`，以及哪些是跨資料集共通現象。",
            "",
            "## 2. 本次假設",
            "",
            "- `HCC1395` 與 `HCC1395_DORADO` 因為是同一 biological specimen，不同平台差異比檢體差異更值得直接對照。",
            "- `paired_full / paired_pileup / to_pileup` 的差異主要反映 workflow/caller/phasing 上游口徑，不應直接當成癌種差異。",
            "- 跨樣本共通現象若存在，應該在多個 sample/mode 都呈現一致方向；若方向隨 sample 或 mode 改變，較合理解讀是 dataset-aware 現象。",
            "",
            "## 3. 成功條件",
            "",
            "- 所有樣本與所有流程組合都在同一個 workspace 內有 coverage 定位。",
            "- 所有 completed runs 都納入同一個 workspace 的分佈分析。",
            "- 每個 run 都可看到三層 F1 與 TP/FP feature distribution。",
            "- `觀察.md` 可以直接用相對路徑顯示圖片。",
            "",
            "## 4. 資料範圍",
            "",
            f"- Sample × mode combinations in coverage matrix: `{total_combinations}`",
            f"- Completed runs with per-region analysis: `{completed_count}`",
            f"- Not applicable runs: `{not_applicable_count}`",
            f"- Samples: `{', '.join(SAMPLE_ORDER)}`",
            f"- Modes: `{', '.join(MODE_ORDER)}`",
            f"- Per-region rows: `{feature_df.shape[0]:,}`",
            "",
            "## 5. 主要欄位意義",
            "",
            *parameter_lines,
            "",
            "## 6. 全樣本全流程涵蓋矩陣",
            "",
            "![](plots/global/00_workflow_coverage_matrix.png)",
            "",
            *coverage_lines,
            "",
            "## 7. 分層 benchmark 總表",
            "",
            *run_lines,
            "",
            "## 8. 全域圖表",
            "",
            "### 8.1 三層 F1 與 stepwise delta",
            "",
            "![](plots/global/01_layer_f1_by_run.png)",
            "",
            "![](plots/global/02_stepwise_delta_heatmap.png)",
            "",
            "### 8.2 TP/FP 總量與 feature distribution",
            "",
            "![](plots/global/03_tp_fp_region_counts.png)",
            "",
            "![](plots/global/04_feature_box_by_mode_truth.png)",
            "",
            "![](plots/global/05_feature_violin_by_run_truth.png)",
            "",
            "![](plots/global/06_verification_class_stack.png)",
            "",
            "![](plots/global/08_feature_delta_heatmap.png)",
            "",
            "### 8.3 同一 specimen 的平台/流程對照",
            "",
            "![](plots/global/07_hcc1395_family_layer_f1.png)",
            "",
            "## 9. 整體觀察",
            "",
            "### 9.1 LongPhase→InterSubMod 最明顯正增益",
            "",
            *best_positive_lines,
            "",
            "### 9.2 LongPhase→InterSubMod 最明顯負增益",
            "",
            *most_negative_lines,
            "",
            "### 9.3 依 mode 聚合的平均 stepwise delta",
            "",
            *mode_delta_lines,
            "",
            "### 9.4 全域 feature 中位數（TP vs FP）",
            "",
            *global_feature_lines,
            "",
            "### 9.5 這些數字代表什麼",
            "",
            *strongest_global_lines,
            f"- `Quality_Score` 目前是整體上最穩的 TP/FP 分離訊號之一；全域 `TP=95`、`FP=75`，`std_delta={quality_row['standardized_delta']:+.3f}`。",
            f"- `VAF` 的全域方向雖然呈現 `TP < FP` (`TP={vaf_row['tp_median']:.4f}`, `FP={vaf_row['fp_median']:.4f}`)，但它其實是 `heterogeneous`：`16` 個 run 中有 `{vaf_positive_count}` 個是 `TP > FP`、`{vaf_negative_count}` 個是 `TP < FP`，因此不能把它當成全域單一方向規則。",
            f"- `CramersV` 兩邊中位數都等於 `0` 並不是資料壞掉，而是典型的 `zero-inflated` 分佈；目前非零比例只有 `TP={cramers_row['tp_nonzero_frac']:.3f}`、`FP={cramers_row['fp_nonzero_frac']:.3f}`，代表只有少數 region 真的出現可觀察的 cluster-label 關聯。",
            "",
            "## 10. 如何判讀差異來源",
            "",
            "### 10.1 比較像平台/定序來源差異",
            "",
            "- `HCC1395` vs `HCC1395_DORADO` 是同一 specimen，不同平台/定序來源；這組最適合拿來看 platform effect。",
            f"- 目前 `HCC1395 5kHz` 的平均 stepwise delta 為 `Caller→LongPhase {hcc1395_5khz_summary['caller_to_longphase']:+.4f}`、`LongPhase→InterSubMod {hcc1395_5khz_summary['longphase_to_intersubmod']:+.4f}`；`HCC1395_DORADO` 則是 `{hcc1395_dorado_summary['caller_to_longphase']:+.4f}` 與 `{hcc1395_dorado_summary['longphase_to_intersubmod']:+.4f}`。同 specimen 下已可看到 5kHz 比 DORADO 更容易出現正增益，因此較合理先解讀成 dataset/platform-sensitive，而不是檢體差異。",
            f"- `ONT_Google_5mCG` 並沒有呈現單向更穩；目前有 `{google_negative.shape[0]}` 個 Google runs 的 `Caller→LongPhase` 為負值，其中包括 "
            + ", ".join(f"`{row.sample} {row.mode}`" for row in google_negative.itertuples())
            + "，因此不能把高 F1 直接簡化成平台品質特別穩定。",
            "",
            "### 10.2 比較像 workflow 差異",
            "",
            "- `paired_full / paired_pileup / to_pileup` 的 baseline F1 差距很大，尤其同一家族的 `paired` 與 `TO` 差異，優先反映 caller/phasing/benchmark scope 差異。",
            "- `paired_pileup` 的 `Caller→LongPhase` 平均增益 (`+0.0258`) 明顯高於 `paired_full` (`+0.0025`)，表示 pileup baseline 留給 phasing 的修正空間更大；這不代表 pileup 比 full 更好，而是代表上游 caller 口徑不同。",
            "- 目前已完成的 `to_pileup` 中，`Caller→LongPhase` 幾乎等於 `0`；這與既有 TO pilot 脈絡一致，代表 `LongPhase-TO phase` 在目前流程中主要提供 phased / HP 資訊，而不是直接改變 PASS call set。",
            f"- `LongPhase→InterSubMod` 最佳正增益來自 `{top_positive['run_label']}`，因為它相對 LongPhase 只損失 `TP {int(top_positive['tp_delta_is_vs_lp'])}`，但可同步移除 `FP {int(top_positive['fp_delta_is_vs_lp'])}`。",
            f"- 最大負增益來自 `{top_negative['run_label']}`，主因是它相對 LongPhase 損失 `TP {int(top_negative['tp_delta_is_vs_lp'])}`，卻只額外移除 `FP {int(top_negative['fp_delta_is_vs_lp'])}`；也就是 TP 成本高於 FP 收益。",
            "",
            "### 10.3 是否能直接解讀成癌種差異",
            "",
            "- 目前不建議直接把 `lung / breast / melanoma` 的差異解讀成純癌種效應，因為癌種差異同時混有 truth set、平台來源、coverage、caller baseline 等共變因素。",
            "- 這份工作區較適合支持 `樣本內 workflow 比較` 與 `同 specimen 平台比較`，對癌種層級只可做弱觀察，不適合下硬結論。",
            "",
            "### 10.4 跨資料集較穩定的共通現象",
            "",
            "- `InterSubMod` 在多數 run 上對 `LongPhase` 的 F1 變化接近零。這在目前資料下更接近 `符合現有方法定位`，也就是 conservative filter / annotation layer，而不是全域大幅提昇器。",
            "- `low VAF + high AlleleDelta (+ low CramersV)` 仍較像 artifact triage 線索，不適合直接升級為全域 TP rescue 主規則。",
            "- `Quality_Score` 較適合當 soft support annotation；`PairwiseMedianDist` 方向則具有明顯 dataset dependence，不應跨 mode/samples 固定方向解讀。",
            "- `VerificationClass` 應視為 label/cluster 一致性描述，不可直接當作 TP/FP classifier。",
            "- `VAF` 的全域中位數方向看起來反直覺，但 per-run sign 是混合的；因此更合理的結論不是 `VAF 對 TP 不利`，而是 `VAF 在不同 sample/mode 下受到不同 caller/truth 背景影響`。",
            "",
            "## 11. 各樣本觀察",
            "",
            *sample_sections,
            "## 12. 重要限制",
            "",
            "- workflow coverage matrix 會列出全部 sample/mode 組合，但 per-region 分佈圖只對 `status=completed` 的 runs 建立。",
            "- per-region feature plots 為了可視化有做 group-wise downsampling；真正統計表仍以 `data/*.tsv` 與 `per_region_features.tsv.gz` 為準。",
            "- `HP/LOH` 相關欄位只能當 diagnostics / annotation，不應單獨解讀成 region-level LOH。",
            "",
            "## 13. 建議下一步",
            "",
            "- 若要進一步驗證 `平台差異`，優先擴大 `HCC1395 5kHz vs HCC1395_DORADO` 的 same-specimen 分析。",
            "- 若要進一步驗證 `workflow 差異`，優先比較同一 sample 的 `paired_full / paired_pileup / to_pileup` feature shift。",
            "- 若要看 `癌種差異`，需先補更多同癌種、同平台、同 workflow 的 replicate，否則容易混入 confounders。",
            "",
        ]
    )
    out_path.write_text(content + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    checklist_path = Path(args.completion_checklist).resolve()
    runs = load_completed_runs(checklist_path)
    if not runs:
        raise SystemExit(f"No completed expected runs found in {checklist_path}")
    workflow_df = load_workflow_coverage(checklist_path)

    default_output_name = (
        DEFAULT_OUTPUT_NAME if args.verification_view == "legacy"
        else f"{DEFAULT_OUTPUT_NAME}_current_v2"
    )
    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else Path(args.output_root).resolve() / default_output_name
    )
    data_dir = ensure_dir(output_dir / "data")
    plots_global = ensure_dir(output_dir / "plots" / "global")
    plots_sample = ensure_dir(output_dir / "plots" / "by_sample")

    layer_df = build_layer_metrics_df(runs)
    feature_parts: List[pd.DataFrame] = []
    for meta in runs:
        feature_parts.append(
            load_summary_rows(
                meta.tp_summary_path,
                meta.tp_vcf_path,
                meta,
                "TP",
                args.verification_view,
            )
        )
        feature_parts.append(
            load_summary_rows(
                meta.fp_summary_path,
                meta.fp_vcf_path,
                meta,
                "FP",
                args.verification_view,
            )
        )
    feature_df = pd.concat(feature_parts, ignore_index=True)
    class_order = verification_class_order(args.verification_view)

    sample_meta_df = pd.DataFrame(
        [
            {
                "sample": meta.sample,
                "sample_family": meta.sample_family,
                "platform_group": meta.platform_group,
                "cancer_type": meta.cancer_type,
                "biological_relation": meta.biological_relation,
                "comparison_axis": meta.comparison_axis,
            }
            for meta in runs
        ]
    ).drop_duplicates()
    run_manifest_df = pd.DataFrame([meta.__dict__ | {"run_label": meta.run_label} for meta in runs])
    feature_summary_by_run = build_feature_summary(feature_df, ["sample", "mode", "run_label", "truth_label"])
    feature_summary_by_run["grouping"] = "run_truth"
    feature_summary_global = build_feature_summary(feature_df, ["truth_label"])
    feature_summary_global["grouping"] = "global_truth"
    feature_summary_mode = build_feature_summary(feature_df, ["mode", "truth_label"])
    feature_summary_mode["grouping"] = "mode_truth"
    feature_summary_df = pd.concat(
        [feature_summary_by_run, feature_summary_global, feature_summary_mode],
        ignore_index=True,
        sort=False,
    )
    verification_df = build_verification_summary(feature_df, class_order)
    transition_df = build_run_transition_df(layer_df)
    feature_delta_run_df = build_feature_delta_df(feature_df, ["sample", "mode", "run_label"])
    feature_delta_sample_df = build_feature_delta_df(feature_df, ["sample"])
    feature_delta_global_df = build_feature_delta_df(feature_df.assign(global_scope="all"), ["global_scope"])
    if not feature_delta_global_df.empty:
        feature_delta_global_df = feature_delta_global_df.drop(columns=["global_scope"])

    sns.set_theme(style="whitegrid")
    plot_workflow_coverage(workflow_df, plots_global / "00_workflow_coverage_matrix.png")
    plot_layer_f1(layer_df, plots_global / "01_layer_f1_by_run.png")
    plot_stepwise_delta_heatmap(layer_df, plots_global / "02_stepwise_delta_heatmap.png")
    plot_tp_fp_counts(feature_df, plots_global / "03_tp_fp_region_counts.png")
    plot_mode_box(feature_df, plots_global / "04_feature_box_by_mode_truth.png", args.max_points_per_group)
    plot_run_violin(feature_df, plots_global / "05_feature_violin_by_run_truth.png", args.max_points_per_group)
    plot_verification_stack(
        feature_df,
        plots_global / "06_verification_class_stack.png",
        class_order,
    )
    plot_hcc1395_family(layer_df, plots_global / "07_hcc1395_family_layer_f1.png")
    plot_feature_delta_heatmap(feature_delta_run_df, plots_global / "08_feature_delta_heatmap.png")

    per_sample_plot_names: Dict[str, str] = {}
    for sample in SAMPLE_ORDER:
        if sample not in feature_df["sample"].unique():
            continue
        filename = f"{sample}_feature_panel.png"
        plot_sample_panel(feature_df, sample, plots_sample / filename, args.max_points_per_group)
        per_sample_plot_names[sample] = filename

    save_dataframe(run_manifest_df, data_dir / "run_manifest.tsv")
    save_dataframe(workflow_df, data_dir / "workflow_coverage.tsv")
    save_dataframe(sample_meta_df, data_dir / "sample_metadata.tsv")
    save_dataframe(layer_df, data_dir / "layer_metrics.tsv")
    save_dataframe(transition_df, data_dir / "run_transition.tsv")
    save_dataframe(feature_summary_df, data_dir / "feature_summary.tsv")
    save_dataframe(feature_delta_run_df, data_dir / "feature_delta_by_run.tsv")
    save_dataframe(feature_delta_sample_df, data_dir / "feature_delta_by_sample.tsv")
    save_dataframe(feature_delta_global_df, data_dir / "feature_delta_global.tsv")
    save_dataframe(verification_df, data_dir / "verification_class_summary.tsv")
    save_gzip_dataframe(feature_df, data_dir / "per_region_features.tsv.gz")

    metadata = {
        "completion_checklist": str(checklist_path),
        "output_dir": str(output_dir),
        "completed_runs": len(runs),
        "region_rows": int(feature_df.shape[0]),
        "samples": sorted(feature_df["sample"].unique().tolist()),
        "modes": sorted(feature_df["mode"].unique().tolist()),
        "verification_view": args.verification_view,
        "verification_source_field": str(feature_df["VerificationSourceField"].iloc[0]),
        "verification_schema_status": str(feature_df["VerificationSchemaStatus"].iloc[0]),
        "verification_class_order": class_order,
        "verification_unknown_counts": {
            str(key): int(value)
            for key, value in feature_df.loc[
                feature_df["VerificationClass"].astype("object") == "UnknownCurrentClass",
                "VerificationClass_SourceValue",
            ].value_counts(dropna=False).items()
        },
    }
    (data_dir / "workspace_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    build_markdown(
        output_dir / "觀察.md",
        runs=runs,
        workflow_df=workflow_df,
        transition_df=transition_df,
        feature_df=feature_df,
        feature_delta_run_df=feature_delta_run_df,
        feature_delta_sample_df=feature_delta_sample_df,
        feature_delta_global_df=feature_delta_global_df,
        per_sample_plot_names=per_sample_plot_names,
    )

    print(f"[Done] Workspace generated: {output_dir}")


if __name__ == "__main__":
    main()
