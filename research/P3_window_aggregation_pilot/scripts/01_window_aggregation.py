#!/usr/bin/env python3
"""P3 pilot: Window-level aggregation 是否突破 region-level AUC ≤0.58？

Goal (ISM 目標 2 / Phase 2B):
  Region-level HPFineNGroups / HPFineF pure methylation 上限 ≤0.58。
  若 gene-level (or window-level) aggregation 能突破，代表 Phase 2B 有潛力。

Method:
  - Window: 1Mb bins (chr + pos // 1e6)
  - Per window aggregate: mean HPFineNGroups, max, median; same for HPFineF, AlleleDelta
  - Window label: majority truth_label (TP if >=50% regions are TP)
  - AUC over windows (TP vs FP)
  - Comparison: region-level AUC of same metrics (same filter)

Output:
  - window_auc_vs_region_auc.tsv
  - window_auc_comparison.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = OUT_DIR / "data"
FIG_DIR = OUT_DIR / "figures"

SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
METRICS = ["HPFineNGroups", "HPFineF", "AlleleDelta", "HPMergedDelta", "PairwiseMedianDist"]
WINDOW_SIZE = 1_000_000  # 1 Mb


def safe_auc(labels, scores):
    try:
        if len(np.unique(labels)) < 2:
            return np.nan
        return roc_auc_score(labels, scores)
    except Exception:
        return np.nan


def analyze_mode(df: pd.DataFrame, mode: str, loh_filter: str) -> pd.DataFrame:
    """Compute region-level vs window-level AUC per sample per metric."""
    sub = df[df["mode"] == mode].copy()
    if loh_filter == "NonLOH":
        sub = sub[sub["to_loh_bed_hit"] == False]
    elif loh_filter == "LOH":
        sub = sub[sub["to_loh_bed_hit"] == True]
    sub = sub[sub["NumReads"] >= 40].copy()
    sub = sub[sub["truth_label"].isin(["TP", "FP"])].copy()
    sub["y"] = (sub["truth_label"] == "TP").astype(int)
    sub["window_id"] = sub["Chr"].astype(str) + "_" + (sub["Pos"] // WINDOW_SIZE).astype(str)

    rows = []
    for sample in SAMPLES:
        s = sub[sub["sample"] == sample]
        if len(s) < 50:
            continue
        region_n = len(s)
        for metric in METRICS:
            if metric not in s.columns:
                continue
            s_m = s[[metric, "y", "window_id"]].dropna()
            if len(s_m) < 50:
                continue
            # Region-level AUC
            region_auc = safe_auc(s_m["y"], s_m[metric])
            # Window-level aggregation
            win_agg = s_m.groupby("window_id").agg(
                n_regions=("y", "size"),
                n_tp=("y", "sum"),
                mean_metric=(metric, "mean"),
                median_metric=(metric, "median"),
                max_metric=(metric, "max"),
            )
            win_agg["y_maj"] = (win_agg["n_tp"] / win_agg["n_regions"] >= 0.5).astype(int)
            # Drop windows with <3 regions (too noisy)
            win_clean = win_agg[win_agg["n_regions"] >= 3]
            if len(win_clean) < 20:
                continue
            window_auc_mean = safe_auc(win_clean["y_maj"], win_clean["mean_metric"])
            window_auc_median = safe_auc(win_clean["y_maj"], win_clean["median_metric"])
            window_auc_max = safe_auc(win_clean["y_maj"], win_clean["max_metric"])
            best_window_auc = max(
                [x for x in [window_auc_mean, window_auc_median, window_auc_max] if not np.isnan(x)],
                default=np.nan,
            )
            rows.append({
                "sample": sample, "mode": mode, "loh_filter": loh_filter, "metric": metric,
                "region_n": region_n, "window_n": len(win_clean),
                "region_auc": region_auc,
                "window_auc_mean": window_auc_mean,
                "window_auc_median": window_auc_median,
                "window_auc_max": window_auc_max,
                "best_window_auc": best_window_auc,
                "delta_vs_region": best_window_auc - region_auc if best_window_auc == best_window_auc else np.nan,
            })
    return pd.DataFrame(rows)


def plot_comparison(result: pd.DataFrame, path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, mode in zip(axes, ["to", "paired"]):
        m_rows = result[(result["mode"] == mode) & (result["loh_filter"] == "NonLOH")]
        if m_rows.empty:
            ax.text(0.5, 0.5, f"{mode}: no data", ha="center", va="center", transform=ax.transAxes)
            continue
        x = m_rows["region_auc"]
        y = m_rows["best_window_auc"]
        colors = ["tab:blue" if s.startswith("HCC1395") else "tab:orange" if s.startswith("HCC") else
                  "tab:green" if s.startswith("H") else "tab:red" for s in m_rows["sample"]]
        ax.scatter(x, y, c=colors, alpha=0.7, s=60)
        for _, r in m_rows.iterrows():
            ax.annotate(f"{r['sample'][:8]}|{r['metric'][:6]}",
                        (r["region_auc"], r["best_window_auc"]), fontsize=7, alpha=0.6)
        lim = [0.4, 0.8]
        ax.plot(lim, lim, "k--", alpha=0.5, label="y=x")
        ax.axhline(0.60, color="red", linestyle=":", label="target=0.60")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("Region-level AUC")
        ax.set_ylabel("Best Window-level AUC")
        ax.set_title(f"{mode} NonLOH NR≥40 (1Mb window)")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close()


def main():
    print("=== P3 pilot: Window-level aggregation AUC ===\n")
    print(f"Loading {MASTER.name}...")
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)
    print(f"  rows: {len(df):,}\n")

    all_results = []
    for mode in ["to", "paired"]:
        for loh in ["NonLOH", "LOH"]:
            print(f"[{mode} {loh}]")
            r = analyze_mode(df, mode, loh)
            if r.empty:
                print(f"  (no data — LOH annotation unavailable in Paired)")
                continue
            all_results.append(r)
            # per-sample best-metric summary
            print(r.groupby("sample")[["region_auc", "best_window_auc", "delta_vs_region"]].max().round(3).to_string())
            print()

    result = pd.concat(all_results, ignore_index=True)
    result.to_csv(DATA_DIR / "window_auc_vs_region_auc.tsv", sep="\t", index=False)

    plot_comparison(result, FIG_DIR / "01_window_vs_region_auc.png")

    # FINAL VERDICT
    print("\n=== FINAL VERDICT ===")
    nonloh_to = result[(result["mode"] == "to") & (result["loh_filter"] == "NonLOH")]
    n_samples_breakthrough = 0
    for sample in SAMPLES:
        s = nonloh_to[nonloh_to["sample"] == sample]
        if s.empty:
            continue
        max_win = s["best_window_auc"].max()
        max_reg = s["region_auc"].max()
        print(f"  {sample}: region_max={max_reg:.3f} → window_max={max_win:.3f} (Δ={max_win-max_reg:+.3f})")
        if max_win >= 0.60 and max_win > max_reg:
            n_samples_breakthrough += 1

    print(f"\n  樣本 with window_AUC ≥0.60 AND > region_AUC: {n_samples_breakthrough}/7")
    if n_samples_breakthrough >= 4:
        print("  → P3 POTENTIAL: 4/7 樣本 window aggregation 突破 region-level")
    elif n_samples_breakthrough >= 2:
        print("  → P3 MIXED: 2-3 樣本突破，需更深 feature engineering")
    else:
        print("  → P3 NEGATIVE: <2 樣本突破，window aggregation 無潛力")


if __name__ == "__main__":
    main()
