#!/usr/bin/env python3
"""
Global relaxed-standard evaluation across pure tumor samples.

Input:
- Existing pure_tumor_eval directory that contains:
  - data/per_variant_decision.tsv.gz
  - data/run_manifest.tsv

Output:
- A subdirectory with:
  - data/*.tsv
  - plots/*.png
  - report/README.md
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

RANDOM_STATE = 42


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Global relaxed standard evaluation")
    parser.add_argument(
        "--eval-dir",
        default="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418",
        help="Existing pure_tumor_eval directory",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: <eval-dir>/global_relaxed_standard_<timestamp>)",
    )
    parser.add_argument(
        "--max-scatter-points",
        type=int,
        default=18000,
        help="Max points per sample for AD-VAF scatter",
    )
    return parser.parse_args()


def ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def metric_from_counts(tp: int, fp: int, truth_total: int) -> Tuple[float, float, float]:
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / truth_total if truth_total > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return precision, recall, f1


def build_sample_baseline(manifest: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in manifest.itertuples(index=False):
        with open(r.metrics, "r") as f:
            m = json.load(f)
        rows.append(
            {
                "sample": r.sample,
                "run_name": r.run_name,
                "run_dir": r.run_dir,
                "truth_total": int(m["truth_total"]),
                "tp_baseline": int(m["baseline"]["tp"]),
                "fp_baseline": int(m["baseline"]["fp"]),
                "precision_baseline": float(m["baseline"]["precision"]),
                "recall_baseline": float(m["baseline"]["recall"]),
                "f1_baseline": float(m["baseline"]["f1"]),
                "tp_current": int(m["filtered"]["tp"]),
                "fp_current": int(m["filtered"]["fp"]),
                "precision_current": float(m["filtered"]["precision"]),
                "recall_current": float(m["filtered"]["recall"]),
                "f1_current": float(m["filtered"]["f1"]),
                "f1_delta_current": float(m["improvement"]["f1_delta"]),
            }
        )
    return pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)


def evaluate_threshold(
    df: pd.DataFrame,
    baseline_df: pd.DataFrame,
    ad_thr: float,
    cv_thr: float,
    vaf_thr: float,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    rows = []

    for base in baseline_df.itertuples(index=False):
        sub = df[df["sample"] == base.sample]
        trigger = (
            (sub["AlleleDelta_abs"] > ad_thr)
            & (sub["CramersV"] < cv_thr)
            & (sub["VAF"] < vaf_thr)
        )
        tp_removed = int((trigger & (sub["is_tp"] == True)).sum())
        fp_removed = int((trigger & (sub["is_tp"] == False)).sum())

        tp = int(base.tp_baseline) - tp_removed
        fp = int(base.fp_baseline) - fp_removed
        p, r, f1 = metric_from_counts(tp, fp, int(base.truth_total))
        f1_delta = f1 - float(base.f1_baseline)

        rows.append(
            {
                "sample": base.sample,
                "ad_thr": ad_thr,
                "cv_thr": cv_thr,
                "vaf_thr": vaf_thr,
                "tp_removed": tp_removed,
                "fp_removed": fp_removed,
                "tp_proposed": tp,
                "fp_proposed": fp,
                "precision_proposed": p,
                "recall_proposed": r,
                "f1_proposed": f1,
                "f1_delta_vs_baseline": f1_delta,
                "f1_delta_vs_current": f1 - float(base.f1_current),
                "net_removed_fp_minus_tp": fp_removed - tp_removed,
            }
        )

    per_sample = pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)
    deltas = per_sample["f1_delta_vs_baseline"]
    improved_count = int((deltas > 0).sum())
    non_negative_count = int((deltas >= 0).sum())
    worst_delta = float(deltas.min())
    mean_delta = float(deltas.mean())
    sum_delta = float(deltas.sum())
    median_delta = float(deltas.median())
    total_fp_minus_tp = int(per_sample["net_removed_fp_minus_tp"].sum())

    summary = {
        "ad_thr": ad_thr,
        "cv_thr": cv_thr,
        "vaf_thr": vaf_thr,
        "improved_count": improved_count,
        "non_negative_count": non_negative_count,
        "worst_delta": worst_delta,
        "mean_delta": mean_delta,
        "median_delta": median_delta,
        "sum_delta": sum_delta,
        "total_fp_minus_tp": total_fp_minus_tp,
    }
    return per_sample, summary


def pick_global_best(score_df: pd.DataFrame) -> pd.Series:
    # Priority:
    # 1) Max non-negative samples
    # 2) Max improved samples
    # 3) Max mean delta
    # 4) Max worst-case delta
    # 5) Max sum delta
    ranked = score_df.sort_values(
        by=["non_negative_count", "improved_count", "mean_delta", "worst_delta", "sum_delta"],
        ascending=[False, False, False, False, False],
    ).reset_index(drop=True)
    return ranked.iloc[0]


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def plot_f1_compare(merged: pd.DataFrame, out_png: Path) -> None:
    plot_df = merged.melt(
        id_vars=["sample"],
        value_vars=["f1_baseline", "f1_current", "f1_proposed"],
        var_name="metric",
        value_name="F1",
    )
    metric_map = {
        "f1_baseline": "Baseline",
        "f1_current": "Current Rule",
        "f1_proposed": "Proposed Relaxed Rule",
    }
    plot_df["metric"] = plot_df["metric"].map(metric_map)

    fig, ax = plt.subplots(figsize=(11, 5))
    sns.barplot(data=plot_df, x="sample", y="F1", hue="metric", ax=ax)
    ax.set_title("F1 Comparison by Sample")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_delta_compare(merged: pd.DataFrame, out_png: Path) -> None:
    plot_df = merged[["sample", "f1_delta_current", "f1_delta_vs_baseline"]].copy()
    plot_df = plot_df.rename(columns={"f1_delta_vs_baseline": "f1_delta_proposed"})
    melt = plot_df.melt(id_vars=["sample"], var_name="type", value_name="F1 delta")
    type_map = {
        "f1_delta_current": "Current vs Baseline",
        "f1_delta_proposed": "Proposed vs Baseline",
    }
    melt["type"] = melt["type"].map(type_map)

    fig, ax = plt.subplots(figsize=(11, 5))
    sns.barplot(data=melt, x="sample", y="F1 delta", hue="type", ax=ax)
    ax.axhline(0, color="black", linewidth=1)
    ax.set_title("F1 Delta Comparison")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_vaf_distribution(df: pd.DataFrame, vaf_thr: float, out_png: Path) -> None:
    samples = sorted(df["sample"].unique())
    n = len(samples)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.5 * rows), sharex=True, sharey=False)
    axes = np.array(axes).reshape(rows, cols)

    for i, s in enumerate(samples):
        r = i // cols
        c = i % cols
        ax = axes[r, c]
        sub = df[df["sample"] == s]
        tp = sub[sub["is_tp"] == True]["VAF"].dropna()
        fp = sub[sub["is_tp"] == False]["VAF"].dropna()
        bins = np.linspace(0, 1.0, 40)
        ax.hist(tp, bins=bins, alpha=0.55, label="TP", color="#1b9e77", density=True)
        ax.hist(fp, bins=bins, alpha=0.55, label="FP", color="#d95f02", density=True)
        ax.axvline(vaf_thr, linestyle="--", color="black", linewidth=1)
        ax.set_title(s)
        if i == 0:
            ax.legend(fontsize=8)
    for j in range(n, rows * cols):
        axes[j // cols, j % cols].set_visible(False)

    fig.suptitle(f"VAF Distribution by Sample (threshold={vaf_thr})", y=1.01)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_ad_distribution(df: pd.DataFrame, ad_thr: float, out_png: Path) -> None:
    samples = sorted(df["sample"].unique())
    n = len(samples)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.5 * rows), sharex=True, sharey=False)
    axes = np.array(axes).reshape(rows, cols)

    for i, s in enumerate(samples):
        r = i // cols
        c = i % cols
        ax = axes[r, c]
        sub = df[df["sample"] == s]
        tp = sub[sub["is_tp"] == True]["AlleleDelta_abs"].dropna()
        fp = sub[sub["is_tp"] == False]["AlleleDelta_abs"].dropna()
        bins = np.linspace(0, 1.0, 40)
        ax.hist(tp, bins=bins, alpha=0.55, label="TP", color="#1b9e77", density=True)
        ax.hist(fp, bins=bins, alpha=0.55, label="FP", color="#d95f02", density=True)
        ax.axvline(ad_thr, linestyle="--", color="black", linewidth=1)
        ax.set_title(s)
        if i == 0:
            ax.legend(fontsize=8)
    for j in range(n, rows * cols):
        axes[j // cols, j % cols].set_visible(False)

    fig.suptitle(f"AlleleDelta Distribution by Sample (threshold={ad_thr})", y=1.01)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_cv_distribution(df: pd.DataFrame, cv_thr: float, out_png: Path) -> None:
    samples = sorted(df["sample"].unique())
    n = len(samples)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.5 * rows), sharex=True, sharey=False)
    axes = np.array(axes).reshape(rows, cols)

    for i, s in enumerate(samples):
        r = i // cols
        c = i % cols
        ax = axes[r, c]
        sub = df[df["sample"] == s]
        tp = sub[sub["is_tp"] == True]["CramersV"].dropna()
        fp = sub[sub["is_tp"] == False]["CramersV"].dropna()
        bins = np.linspace(0, max(0.1, float(sub["CramersV"].max())), 40)
        ax.hist(tp, bins=bins, alpha=0.55, label="TP", color="#1b9e77", density=True)
        ax.hist(fp, bins=bins, alpha=0.55, label="FP", color="#d95f02", density=True)
        ax.axvline(cv_thr, linestyle="--", color="black", linewidth=1)
        ax.set_title(s)
        if i == 0:
            ax.legend(fontsize=8)
    for j in range(n, rows * cols):
        axes[j // cols, j % cols].set_visible(False)

    fig.suptitle(f"CramersV Distribution by Sample (threshold={cv_thr})", y=1.01)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_scatter_ad_vaf(
    df: pd.DataFrame,
    ad_thr: float,
    vaf_thr: float,
    max_points: int,
    out_png: Path,
) -> None:
    samples = sorted(df["sample"].unique())
    n = len(samples)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.8 * rows), sharex=True, sharey=True)
    axes = np.array(axes).reshape(rows, cols)

    for i, s in enumerate(samples):
        ax = axes[i // cols, i % cols]
        sub = df[df["sample"] == s][["AlleleDelta_abs", "VAF", "is_tp"]].dropna()
        if len(sub) > max_points:
            sub = sub.sample(max_points, random_state=RANDOM_STATE)
        colors = np.where(sub["is_tp"].values, "#1b9e77", "#d95f02")
        ax.scatter(sub["AlleleDelta_abs"], sub["VAF"], c=colors, s=6, alpha=0.35, linewidths=0)
        # trigger region: AD > ad_thr and VAF < vaf_thr
        ax.axvline(ad_thr, color="black", linestyle="--", linewidth=1)
        ax.axhline(vaf_thr, color="black", linestyle="--", linewidth=1)
        ax.fill_between(
            x=np.array([ad_thr, 1.0]),
            y1=0,
            y2=vaf_thr,
            color="#bbbbbb",
            alpha=0.25,
        )
        ax.set_title(s)
    for j in range(n, rows * cols):
        axes[j // cols, j % cols].set_visible(False)
    fig.suptitle("AD-VAF Scatter by Sample (gray area = potential trigger zone)", y=1.01)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(
    report_path: Path,
    best: pd.Series,
    merged: pd.DataFrame,
    impossible_all_improve: bool,
) -> None:
    improved = int((merged["f1_delta_vs_baseline"] > 0).sum())
    non_negative = int((merged["f1_delta_vs_baseline"] >= 0).sum())
    n = int(len(merged))

    lines: List[str] = []
    lines.append("<!--")
    lines.append(f"建立時間: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("目標: 以單一放寬標準評估所有純 tumor 樣本，輸出TP/FP/F1與分佈觀察圖")
    lines.append("-->")
    lines.append("")
    lines.append("# 全樣本放寬標準評估報告")
    lines.append("")
    lines.append("## 1. 全局最佳單一標準")
    lines.append("")
    lines.append(f"- AlleleDelta > `{best['ad_thr']}`")
    lines.append(f"- CramersV < `{best['cv_thr']}`")
    lines.append(f"- VAF < `{best['vaf_thr']}`")
    lines.append("")
    lines.append("此標準下（相對 baseline）：")
    lines.append(f"- 改善樣本數：`{improved}/{n}`")
    lines.append(f"- 非負樣本數：`{non_negative}/{n}`")
    lines.append(f"- 平均 F1 delta：`{best['mean_delta']:+.6f}`")
    lines.append(f"- 最差 F1 delta：`{best['worst_delta']:+.6f}`")
    lines.append(f"- 全樣本 delta 總和：`{best['sum_delta']:+.6f}`")
    lines.append("")
    if impossible_all_improve:
        lines.append("- 觀察：在本次搜尋範圍內，**不存在**可讓所有樣本同時 F1 提升的單一門檻。")
        lines.append("")

    lines.append("## 2. 完整結果表")
    lines.append("")
    lines.append("- [per_sample_metrics_proposed.tsv](../data/per_sample_metrics_proposed.tsv)")
    lines.append("- [grid_score_summary.tsv](../data/grid_score_summary.tsv)")
    lines.append("- [best_global_standard.tsv](../data/best_global_standard.tsv)")
    lines.append("")
    lines.append("## 3. 甲基分佈觀察圖")
    lines.append("")
    lines.append("- [F1 三組比較（Baseline/Current/Proposed）](../plots/f1_compare_baseline_current_proposed.png)")
    lines.append("- [F1 delta 比較（Current vs Proposed）](../plots/f1_delta_compare_current_vs_proposed.png)")
    lines.append("- [VAF 分佈（各樣本 TP/FP）](../plots/vaf_distribution_by_sample.png)")
    lines.append("- [AlleleDelta 分佈（各樣本 TP/FP）](../plots/allele_delta_distribution_by_sample.png)")
    lines.append("- [CramersV 分佈（各樣本 TP/FP）](../plots/cramersv_distribution_by_sample.png)")
    lines.append("- [AD-VAF 散點與觸發區域](../plots/ad_vaf_scatter_by_sample.png)")

    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    np.random.seed(RANDOM_STATE)
    sns.set_theme(style="whitegrid")

    eval_dir = Path(args.eval_dir)
    data_dir_in = eval_dir / "data"
    per_variant_path = data_dir_in / "per_variant_decision.tsv.gz"
    manifest_path = data_dir_in / "run_manifest.tsv"

    if not per_variant_path.exists() or not manifest_path.exists():
        raise FileNotFoundError("Missing per_variant_decision.tsv.gz or run_manifest.tsv in eval-dir/data.")

    if args.output_dir:
        out_root = Path(args.output_dir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_root = eval_dir / f"global_relaxed_standard_{ts}"

    data_out = ensure(out_root / "data")
    plots_out = ensure(out_root / "plots")
    report_out = ensure(out_root / "report")

    df = pd.read_csv(per_variant_path, sep="\t", low_memory=False)
    manifest = pd.read_csv(manifest_path, sep="\t")
    baseline_df = build_sample_baseline(manifest)

    ad_grid = [0.15, 0.17, 0.20, 0.25, 0.30]
    cv_grid = [0.03, 0.025, 0.02, 0.015, 0.01]
    vaf_grid = [0.15, 0.12, 0.10, 0.08, 0.05]

    score_rows: List[Dict[str, float]] = []
    per_combo: Dict[Tuple[float, float, float], pd.DataFrame] = {}

    for ad_thr, cv_thr, vaf_thr in itertools.product(ad_grid, cv_grid, vaf_grid):
        per_sample, score = evaluate_threshold(df, baseline_df, ad_thr, cv_thr, vaf_thr)
        score_rows.append(score)
        per_combo[(ad_thr, cv_thr, vaf_thr)] = per_sample

    score_df = pd.DataFrame(score_rows)
    best = pick_global_best(score_df)
    key = (float(best["ad_thr"]), float(best["cv_thr"]), float(best["vaf_thr"]))
    best_per_sample = per_combo[key].copy()

    merged = baseline_df.merge(best_per_sample, on="sample", how="left")
    merged = merged.sort_values("sample").reset_index(drop=True)

    impossible_all_improve = not (score_df["improved_count"] == len(merged)).any()

    # Save tables
    save_tsv(score_df.sort_values(
        by=["non_negative_count", "improved_count", "mean_delta", "worst_delta", "sum_delta"],
        ascending=[False, False, False, False, False],
    ), data_out / "grid_score_summary.tsv")
    save_tsv(pd.DataFrame([best]), data_out / "best_global_standard.tsv")
    save_tsv(merged, data_out / "per_sample_metrics_proposed.tsv")

    # Save per-variant with proposed decision for traceability
    trigger_best = (
        (df["AlleleDelta_abs"] > key[0])
        & (df["CramersV"] < key[1])
        & (df["VAF"] < key[2])
    )
    proposed = df.copy()
    proposed["rule_triggered_proposed"] = trigger_best
    proposed["decision_class_proposed"] = np.where(
        proposed["is_tp"] & (~trigger_best),
        "TP_correct_keep",
        np.where(
            proposed["is_tp"] & trigger_best,
            "TP_error_removed",
            np.where(trigger_best, "FP_correct_removed", "FP_error_kept"),
        ),
    )
    with gzip.open(data_out / "per_variant_decision_proposed.tsv.gz", "wt") as gz:
        proposed.to_csv(gz, sep="\t", index=False)

    # Plots
    plot_f1_compare(merged, plots_out / "f1_compare_baseline_current_proposed.png")
    plot_delta_compare(merged, plots_out / "f1_delta_compare_current_vs_proposed.png")
    plot_vaf_distribution(df, key[2], plots_out / "vaf_distribution_by_sample.png")
    plot_ad_distribution(df, key[0], plots_out / "allele_delta_distribution_by_sample.png")
    plot_cv_distribution(df, key[1], plots_out / "cramersv_distribution_by_sample.png")
    plot_scatter_ad_vaf(df, key[0], key[2], args.max_scatter_points, plots_out / "ad_vaf_scatter_by_sample.png")

    # Report
    write_report(report_out / "README.md", best, merged, impossible_all_improve)

    print("[Done] Global relaxed standard evaluation complete")
    print(f"[Done] Output: {out_root}")
    print(
        "[Done] Best standard: "
        f"AD>{key[0]}, CV<{key[1]}, VAF<{key[2]}"
    )


if __name__ == "__main__":
    main()

