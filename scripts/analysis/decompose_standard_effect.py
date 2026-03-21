#!/usr/bin/env python3
"""
Decompose global filtering standard into single/pair/full rule effects.

Given a selected standard (A, C, V thresholds):
  A: AlleleDelta_abs > ad_thr
  C: CramersV < cv_thr
  V: VAF < vaf_thr

This script evaluates all subsets of {A,C,V}, computes per-sample TP/FP/F1,
and estimates each condition's marginal contribution via Shapley values.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


FEATURES = ["A", "C", "V"]
FEATURE_LABEL = {"A": "AD", "C": "CV", "V": "VAF"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Decompose filtering standard effects")
    parser.add_argument(
        "--eval-dir",
        default="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418",
        help="Base pure_tumor_eval directory (contains data/per_variant_decision.tsv.gz)",
    )
    parser.add_argument(
        "--global-dir",
        default="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418/global_relaxed_standard_20260305_152027",
        help="Global standard output directory (contains data/best_global_standard.tsv)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: <global-dir>/standard_decomposition_<timestamp>)",
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


def subset_name(subset: Set[str]) -> str:
    if not subset:
        return "none"
    return "+".join(FEATURE_LABEL[x] for x in sorted(subset))


def subset_mask(df: pd.DataFrame, subset: Set[str], ad_thr: float, cv_thr: float, vaf_thr: float) -> pd.Series:
    if not subset:
        return pd.Series(False, index=df.index)

    mask = pd.Series(True, index=df.index)
    if "A" in subset:
        mask &= df["AlleleDelta_abs"] > ad_thr
    if "C" in subset:
        mask &= df["CramersV"] < cv_thr
    if "V" in subset:
        mask &= df["VAF"] < vaf_thr
    return mask


def all_subsets() -> List[Set[str]]:
    out: List[Set[str]] = []
    for r in range(0, len(FEATURES) + 1):
        for comb in itertools.combinations(FEATURES, r):
            out.append(set(comb))
    return out


def compute_baseline_map(manifest_df: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    m: Dict[str, Dict[str, float]] = {}
    for r in manifest_df.itertuples(index=False):
        with open(r.metrics, "r") as f:
            j = json.load(f)
        truth_total = int(j["truth_total"])
        tp_baseline = int(j["baseline"]["tp"])
        fp_baseline = int(j["baseline"]["fp"])
        _, _, f1_baseline_exact = metric_from_counts(tp_baseline, fp_baseline, truth_total)

        tp_current = int(j["filtered"]["tp"])
        fp_current = int(j["filtered"]["fp"])
        _, _, f1_current_exact = metric_from_counts(tp_current, fp_current, truth_total)
        m[r.sample] = {
            "truth_total": truth_total,
            "tp_baseline": tp_baseline,
            "fp_baseline": fp_baseline,
            "f1_baseline": f1_baseline_exact,
            "f1_current": f1_current_exact,
            "f1_delta_current": f1_current_exact - f1_baseline_exact,
        }
    return m


def evaluate_subsets(
    df: pd.DataFrame,
    baseline_map: Dict[str, Dict[str, float]],
    ad_thr: float,
    cv_thr: float,
    vaf_thr: float,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    subsets = all_subsets()
    samples = sorted(df["sample"].unique())

    for sample in samples:
        sub_df = df[df["sample"] == sample]
        base = baseline_map[sample]
        tp0 = int(base["tp_baseline"])
        fp0 = int(base["fp_baseline"])
        truth = int(base["truth_total"])
        f1_base = float(base["f1_baseline"])

        for subset in subsets:
            mask = subset_mask(sub_df, subset, ad_thr, cv_thr, vaf_thr)
            tp_removed = int((mask & (sub_df["is_tp"] == True)).sum())
            fp_removed = int((mask & (sub_df["is_tp"] == False)).sum())
            tp = tp0 - tp_removed
            fp = fp0 - fp_removed
            p, r, f1 = metric_from_counts(tp, fp, truth)

            rows.append(
                {
                    "sample": sample,
                    "subset_key": "".join(sorted(subset)),
                    "subset_name": subset_name(subset),
                    "subset_size": len(subset),
                    "tp_removed": tp_removed,
                    "fp_removed": fp_removed,
                    "net_removed_fp_minus_tp": fp_removed - tp_removed,
                    "tp": tp,
                    "fp": fp,
                    "precision": p,
                    "recall": r,
                    "f1": f1,
                    "f1_delta_vs_baseline": f1 - f1_base,
                    "ad_thr": ad_thr,
                    "cv_thr": cv_thr,
                    "vaf_thr": vaf_thr,
                }
            )
    return pd.DataFrame(rows)


def build_global_summary(per_subset_df: pd.DataFrame) -> pd.DataFrame:
    g = (
        per_subset_df.groupby(["subset_key", "subset_name", "subset_size"], as_index=False)
        .agg(
            mean_f1_delta=("f1_delta_vs_baseline", "mean"),
            median_f1_delta=("f1_delta_vs_baseline", "median"),
            min_f1_delta=("f1_delta_vs_baseline", "min"),
            max_f1_delta=("f1_delta_vs_baseline", "max"),
            sum_f1_delta=("f1_delta_vs_baseline", "sum"),
            mean_net_removed_fp_minus_tp=("net_removed_fp_minus_tp", "mean"),
            sum_net_removed_fp_minus_tp=("net_removed_fp_minus_tp", "sum"),
        )
        .sort_values(["subset_size", "mean_f1_delta"], ascending=[True, False])
        .reset_index(drop=True)
    )
    return g


def build_shapley(per_subset_df: pd.DataFrame) -> pd.DataFrame:
    subset_value: Dict[str, Dict[str, float]] = {}
    for sample, sdf in per_subset_df.groupby("sample"):
        subset_value[sample] = dict(zip(sdf["subset_key"], sdf["f1_delta_vs_baseline"]))

    perms = list(itertools.permutations(FEATURES))
    rows: List[Dict[str, object]] = []

    def key_of(s: Set[str]) -> str:
        return "".join(sorted(s))

    for sample, vmap in subset_value.items():
        for f in FEATURES:
            contrib = 0.0
            for perm in perms:
                cur: Set[str] = set()
                for p in perm:
                    before = key_of(cur)
                    if p == f:
                        after_set = set(cur)
                        after_set.add(f)
                        after = key_of(after_set)
                        contrib += vmap.get(after, 0.0) - vmap.get(before, 0.0)
                        break
                    cur.add(p)
            contrib /= len(perms)
            rows.append(
                {
                    "sample": sample,
                    "feature": FEATURE_LABEL[f],
                    "shapley_delta": contrib,
                }
            )

    out = pd.DataFrame(rows)
    global_rows = (
        out.groupby("feature", as_index=False)["shapley_delta"]
        .mean()
        .assign(sample="GLOBAL_MEAN")
    )
    out = pd.concat([out, global_rows[["sample", "feature", "shapley_delta"]]], ignore_index=True)
    return out


def build_ablation(per_subset_df: pd.DataFrame) -> pd.DataFrame:
    # Compare full (ACV) vs removing one feature (AC, AV, CV).
    samples = sorted(per_subset_df["sample"].unique())
    rows = []
    for s in samples:
        sdf = per_subset_df[per_subset_df["sample"] == s]
        full = float(sdf[sdf["subset_key"] == "ACV"]["f1_delta_vs_baseline"].iloc[0])
        ac = float(sdf[sdf["subset_key"] == "AC"]["f1_delta_vs_baseline"].iloc[0])
        av = float(sdf[sdf["subset_key"] == "AV"]["f1_delta_vs_baseline"].iloc[0])
        cv = float(sdf[sdf["subset_key"] == "CV"]["f1_delta_vs_baseline"].iloc[0])
        rows.extend(
            [
                {"sample": s, "ablation": "remove_VAF", "f1_delta_drop": full - ac},
                {"sample": s, "ablation": "remove_CV", "f1_delta_drop": full - av},
                {"sample": s, "ablation": "remove_AD", "f1_delta_drop": full - cv},
            ]
        )
    return pd.DataFrame(rows)


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def plot_subset_delta(global_df: pd.DataFrame, out_png: Path) -> None:
    order = (
        global_df.sort_values(["subset_size", "mean_f1_delta"], ascending=[True, False])["subset_name"]
        .drop_duplicates()
        .tolist()
    )
    fig, ax = plt.subplots(figsize=(11, 5))
    sns.barplot(data=global_df, x="subset_name", y="mean_f1_delta", order=order, ax=ax, color="#4C72B0")
    ax.axhline(0, color="black", linewidth=1)
    ax.set_title("Mean F1 Delta by Subset Rule")
    ax.set_xlabel("Rule subset")
    ax.set_ylabel("Mean F1 delta vs baseline")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_shapley_heatmap(shapley_df: pd.DataFrame, out_png: Path) -> None:
    p = shapley_df[shapley_df["sample"] != "GLOBAL_MEAN"].pivot(
        index="sample", columns="feature", values="shapley_delta"
    )
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(p, annot=True, fmt=".6f", cmap="RdYlGn", center=0, ax=ax)
    ax.set_title("Shapley Contribution (F1 delta)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_ablation(ablation_df: pd.DataFrame, out_png: Path) -> None:
    g = (
        ablation_df.groupby("ablation", as_index=False)["f1_delta_drop"]
        .mean()
        .sort_values("f1_delta_drop", ascending=False)
    )
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.barplot(data=g, x="ablation", y="f1_delta_drop", ax=ax, color="#DD8452")
    ax.axhline(0, color="black", linewidth=1)
    ax.set_title("Ablation Impact (mean drop from full rule)")
    ax.set_ylabel("Mean F1 delta drop")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def write_report(
    path: Path,
    ad_thr: float,
    cv_thr: float,
    vaf_thr: float,
    global_df: pd.DataFrame,
    shapley_df: pd.DataFrame,
    ablation_df: pd.DataFrame,
) -> None:
    full_row = global_df[global_df["subset_key"] == "ACV"].iloc[0]
    single_rows = global_df[global_df["subset_size"] == 1].sort_values("mean_f1_delta", ascending=False)
    global_shapley = (
        shapley_df[shapley_df["sample"] == "GLOBAL_MEAN"]
        .sort_values("shapley_delta", ascending=False)
        .reset_index(drop=True)
    )
    abl_mean = (
        ablation_df.groupby("ablation", as_index=False)["f1_delta_drop"]
        .mean()
        .sort_values("f1_delta_drop", ascending=False)
    )

    lines: List[str] = []
    lines.append("<!--")
    lines.append(f"建立時間: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("目標: 拆解單一標準對TP/FP/F1提升的邊際影響並比較")
    lines.append("-->")
    lines.append("")
    lines.append("# 標準拆解分析報告")
    lines.append("")
    lines.append("## 1. 被拆解標準")
    lines.append("")
    lines.append(f"- Full rule: `AD>{ad_thr}` AND `CV<{cv_thr}` AND `VAF<{vaf_thr}`")
    lines.append("")
    lines.append("## 2. 單一標準影響（全樣本平均）")
    lines.append("")
    for r in single_rows.itertuples(index=False):
        lines.append(
            f"- {r.subset_name}: mean F1 delta `{r.mean_f1_delta:+.6f}`, "
            f"worst `{r.min_f1_delta:+.6f}`, sum `{r.sum_f1_delta:+.6f}`"
        )
    lines.append("")
    lines.append("## 3. 完整標準（ACV）表現")
    lines.append("")
    lines.append(
        f"- mean F1 delta `{full_row['mean_f1_delta']:+.6f}`, "
        f"worst `{full_row['min_f1_delta']:+.6f}`, "
        f"sum `{full_row['sum_f1_delta']:+.6f}`"
    )
    lines.append("")
    lines.append("## 4. 邊際貢獻（Shapley, GLOBAL_MEAN）")
    lines.append("")
    for r in global_shapley.itertuples(index=False):
        lines.append(f"- {r.feature}: `{r.shapley_delta:+.6f}`")
    lines.append("")
    lines.append("## 5. Ablation（從 Full 移除單一條件）")
    lines.append("")
    for r in abl_mean.itertuples(index=False):
        lines.append(f"- {r.ablation}: mean drop `{r.f1_delta_drop:+.6f}`")
    lines.append("")
    lines.append("## 6. 產出檔案")
    lines.append("")
    lines.append("- [subset_performance_per_sample.tsv](../data/subset_performance_per_sample.tsv)")
    lines.append("- [subset_global_summary.tsv](../data/subset_global_summary.tsv)")
    lines.append("- [single_standard_comparison.tsv](../data/single_standard_comparison.tsv)")
    lines.append("- [shapley_contribution.tsv](../data/shapley_contribution.tsv)")
    lines.append("- [ablation_impact.tsv](../data/ablation_impact.tsv)")
    lines.append("- [subset_mean_delta_bar.png](../plots/subset_mean_delta_bar.png)")
    lines.append("- [shapley_heatmap.png](../plots/shapley_heatmap.png)")
    lines.append("- [ablation_impact_bar.png](../plots/ablation_impact_bar.png)")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    sns.set_theme(style="whitegrid")

    eval_dir = Path(args.eval_dir)
    global_dir = Path(args.global_dir)
    if args.output_dir:
        out_root = Path(args.output_dir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_root = global_dir / f"standard_decomposition_{ts}"

    data_out = ensure(out_root / "data")
    plots_out = ensure(out_root / "plots")
    report_out = ensure(out_root / "report")

    per_variant = pd.read_csv(eval_dir / "data" / "per_variant_decision.tsv.gz", sep="\t", low_memory=False)
    manifest_df = pd.read_csv(eval_dir / "data" / "run_manifest.tsv", sep="\t")
    best_std = pd.read_csv(global_dir / "data" / "best_global_standard.tsv", sep="\t").iloc[0]

    ad_thr = float(best_std["ad_thr"])
    cv_thr = float(best_std["cv_thr"])
    vaf_thr = float(best_std["vaf_thr"])

    baseline_map = compute_baseline_map(manifest_df)
    per_subset = evaluate_subsets(per_variant, baseline_map, ad_thr, cv_thr, vaf_thr)
    global_summary = build_global_summary(per_subset)
    single_std = global_summary[global_summary["subset_size"] == 1].copy().sort_values("mean_f1_delta", ascending=False)
    shapley = build_shapley(per_subset)
    ablation = build_ablation(per_subset)

    save_tsv(per_subset, data_out / "subset_performance_per_sample.tsv")
    save_tsv(global_summary, data_out / "subset_global_summary.tsv")
    save_tsv(single_std, data_out / "single_standard_comparison.tsv")
    save_tsv(shapley, data_out / "shapley_contribution.tsv")
    save_tsv(ablation, data_out / "ablation_impact.tsv")

    plot_subset_delta(global_summary, plots_out / "subset_mean_delta_bar.png")
    plot_shapley_heatmap(shapley, plots_out / "shapley_heatmap.png")
    plot_ablation(ablation, plots_out / "ablation_impact_bar.png")

    write_report(
        report_out / "README.md",
        ad_thr,
        cv_thr,
        vaf_thr,
        global_summary,
        shapley,
        ablation,
    )

    print("[Done] Standard decomposition complete")
    print(f"[Done] Output: {out_root}")
    print(f"[Done] Standard: AD>{ad_thr}, CV<{cv_thr}, VAF<{vaf_thr}")


if __name__ == "__main__":
    main()
