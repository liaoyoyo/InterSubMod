#!/usr/bin/env python3
"""P4 pilot: 二次打擊順序在 AF × methylation pattern 聯合分佈上可否區分？

ISM 目標 3 = 推斷二次打擊順序（two-hit hypothesis）

假設邏輯：
  方向 1 (methylation silencing → then LOH):
    - Silencing 先發生在 germline-heterozygous site
    - 之後 LOH 把一個 allele 移除
    - 結果：intermediate AF region 應有 "biallelic methylation pattern"
           即 AlleleDelta 分散（0 到 1 間變異）、HPFineNGroups 高
  方向 2 (LOH → then methylation silencing):
    - LOH 先發生
    - 之後 methylation 在 surviving allele 累積
    - 結果：intermediate AF region 應有 "concentrated methylation"
           即 AlleleDelta 靠近 0 或 1、HPFineNGroups 較低

方法（避免 P3 spatial auto-correlation confound）:
  - 使用 region-level 特徵（無 aggregation）
  - 按 AF bin × LOH flag 分層
  - 比較 AlleleDelta distribution shape + HPFineNGroups mean
  - 核心對比：
      LOH + intermediate AF (0.3-0.7) vs Non-LOH + intermediate AF
      → 若兩者 methylation pattern 分佈型態顯著差異，支持 two-hit order 可區分
      → 若相同，則單 region 特徵無法支撐該目標

小規模驗證範圍:
  - 7 樣本 TO mode (有 LOH annotation)
  - NR ≥ 40 TP+FP regions
  - AF bin: 0.1-0.3 / 0.3-0.5 / 0.5-0.7
  - 指標: AlleleDelta bimodality (Kurtosis / Bimodality Coefficient)
         HPFineNGroups mean + variance

判準:
  - 若 LOH vs Non-LOH AF-matched subset 的 bimodality 差異 effect size >0.3 → POSITIVE（值得深究）
  - 若差異 <0.15 → NEGATIVE（region-level 無二次打擊順序訊號）
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = OUT_DIR / "data"
FIG_DIR = OUT_DIR / "figures"

SAMPLES = ["HCC1395", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
AF_BINS = [(0.1, 0.3), (0.3, 0.5), (0.5, 0.7)]
NR_MIN = 40


def bimodality_coefficient(x: np.ndarray) -> float:
    """SAS bimodality coefficient: BC = (skew^2 + 1) / (kurt + 3*(n-1)^2/((n-2)(n-3)))
    BC > 5/9 ≈ 0.556 suggests bimodal."""
    x = np.asarray(x)
    x = x[~np.isnan(x)]
    if len(x) < 10:
        return np.nan
    n = len(x)
    g1 = stats.skew(x)
    g2 = stats.kurtosis(x, fisher=True)
    return (g1 ** 2 + 1) / (g2 + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a), np.asarray(b)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 5 or len(b) < 5:
        return np.nan
    pooled_sd = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
    if pooled_sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled_sd


def analyze_sample(df: pd.DataFrame, sample: str) -> pd.DataFrame:
    s = df[(df["sample"] == sample) & (df["mode"] == "to")]
    s = s[(s["NumReads"] >= NR_MIN) & s["truth_label"].isin(["TP", "FP"])]
    s = s.dropna(subset=["caller_af", "AlleleDelta", "HPFineNGroups"]).copy()

    rows = []
    for af_lo, af_hi in AF_BINS:
        af_sub = s[(s["caller_af"] >= af_lo) & (s["caller_af"] < af_hi)]
        for loh_flag, loh_name in [(True, "LOH"), (False, "NonLOH")]:
            subset = af_sub[af_sub["to_loh_bed_hit"] == loh_flag]
            if len(subset) < 20:
                continue
            for truth in ["TP", "FP", "ALL"]:
                if truth == "ALL":
                    t_sub = subset
                else:
                    t_sub = subset[subset["truth_label"] == truth]
                if len(t_sub) < 20:
                    continue
                ad = t_sub["AlleleDelta"].values
                ng = t_sub["HPFineNGroups"].values
                rows.append({
                    "sample": sample, "af_bin": f"{af_lo:.1f}-{af_hi:.1f}",
                    "loh": loh_name, "truth": truth,
                    "n": len(t_sub),
                    "ad_mean": np.nanmean(ad),
                    "ad_std": np.nanstd(ad),
                    "ad_bc": bimodality_coefficient(ad),
                    "ng_mean": np.nanmean(ng),
                    "ng_std": np.nanstd(ng),
                    "ng_high_frac": (ng >= 3).mean(),
                })
    return pd.DataFrame(rows)


def compute_loh_vs_nonloh_delta(result: pd.DataFrame) -> pd.DataFrame:
    """For each (sample, af_bin, truth), compute LOH vs NonLOH effect size on AD/NG."""
    rows = []
    for (sample, af_bin, truth), g in result.groupby(["sample", "af_bin", "truth"]):
        loh = g[g["loh"] == "LOH"]
        non = g[g["loh"] == "NonLOH"]
        if loh.empty or non.empty:
            continue
        rows.append({
            "sample": sample, "af_bin": af_bin, "truth": truth,
            "n_loh": loh["n"].values[0], "n_non": non["n"].values[0],
            "delta_ad_bc": loh["ad_bc"].values[0] - non["ad_bc"].values[0],
            "delta_ad_std": loh["ad_std"].values[0] - non["ad_std"].values[0],
            "delta_ng_mean": loh["ng_mean"].values[0] - non["ng_mean"].values[0],
            "delta_ng_high_frac": loh["ng_high_frac"].values[0] - non["ng_high_frac"].values[0],
        })
    return pd.DataFrame(rows)


def plot_joint_heatmap(df: pd.DataFrame, out_path: Path):
    """Plot AlleleDelta × HPFineNGroups heatmap per (sample, loh, af_bin)."""
    df = df[(df["mode"] == "to") & (df["truth_label"].isin(["TP", "FP"])) & (df["NumReads"] >= NR_MIN)]
    df = df.dropna(subset=["caller_af", "AlleleDelta", "HPFineNGroups"])
    fig, axes = plt.subplots(len(SAMPLES), 6, figsize=(24, 4 * len(SAMPLES)), constrained_layout=True)
    for i, sample in enumerate(SAMPLES):
        s = df[df["sample"] == sample]
        for j, (af_lo, af_hi) in enumerate(AF_BINS):
            for k, loh_flag in enumerate([True, False]):
                col = j * 2 + k
                ax = axes[i, col] if len(SAMPLES) > 1 else axes[col]
                sub = s[(s["caller_af"] >= af_lo) & (s["caller_af"] < af_hi) & (s["to_loh_bed_hit"] == loh_flag)]
                if len(sub) < 10:
                    ax.text(0.5, 0.5, f"n={len(sub)}", ha="center", va="center", transform=ax.transAxes)
                    ax.set_title(f"{sample}\nAF{af_lo:.1f}-{af_hi:.1f} {'LOH' if loh_flag else 'Non'}", fontsize=8)
                    continue
                h = ax.hist2d(sub["AlleleDelta"], sub["HPFineNGroups"],
                              bins=[20, 5], cmap="viridis", cmin=1)
                ax.set_title(f"{sample}\nAF{af_lo:.1f}-{af_hi:.1f} {'LOH' if loh_flag else 'Non'}\nn={len(sub)}",
                             fontsize=8)
                ax.set_xlabel("AlleleDelta", fontsize=7)
                ax.set_ylabel("HPFineNGroups", fontsize=7)
                ax.tick_params(labelsize=6)
    fig.savefig(out_path, dpi=120)
    plt.close()


def main():
    print("=== P4 pilot: AF × Methylation Pattern 聯合分佈 (二次打擊順序) ===\n")
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)
    print(f"  master rows: {len(df):,}\n")

    all_rows = []
    for sample in SAMPLES:
        r = analyze_sample(df, sample)
        if r.empty:
            print(f"  {sample}: no qualifying regions")
            continue
        all_rows.append(r)
    result = pd.concat(all_rows, ignore_index=True)
    result.to_csv(DATA_DIR / "af_loh_methylation_joint_stats.tsv", sep="\t", index=False)

    print("=== Per-sample x AF bin x LOH summary (TP+FP pooled) ===")
    all_truth = result[result["truth"] == "ALL"]
    print(all_truth.pivot_table(
        index=["sample", "af_bin"], columns="loh",
        values=["n", "ad_bc", "ad_std", "ng_mean", "ng_high_frac"]
    ).round(3).to_string())

    # LOH vs NonLOH delta
    delta = compute_loh_vs_nonloh_delta(result[result["truth"] == "ALL"])
    delta.to_csv(DATA_DIR / "loh_vs_nonloh_delta.tsv", sep="\t", index=False)

    print("\n=== LOH vs NonLOH deltas (ALL regions, same AF bin) ===")
    print(delta.round(3).to_string(index=False))

    # VERDICT
    print("\n=== VERDICT ===")
    intermediate = delta[delta["af_bin"] == "0.3-0.5"]
    if intermediate.empty:
        print("  (no intermediate AF 0.3-0.5 data)")
    else:
        mean_delta_bc = intermediate["delta_ad_bc"].abs().mean()
        mean_delta_ng = intermediate["delta_ng_mean"].abs().mean()
        print(f"  Intermediate AF (0.3-0.5) LOH vs NonLOH:")
        print(f"    mean |Δ BimodalityCoef|: {mean_delta_bc:.3f}")
        print(f"    mean |Δ NGroups|:         {mean_delta_ng:.3f}")
        if mean_delta_bc >= 0.15 or mean_delta_ng >= 0.5:
            print("  → POSITIVE: 二次打擊順序可能在 region-level 特徵上可區分")
        elif mean_delta_bc >= 0.08 or mean_delta_ng >= 0.3:
            print("  → MIXED: 微弱訊號；需 per-sample consistency 檢查")
        else:
            print("  → NEGATIVE: region-level AF×methylation 聯合分佈 LOH vs NonLOH 差異微小")

    # consistency across samples
    print("\n=== Per-sample consistency (intermediate AF) ===")
    if not intermediate.empty:
        print(intermediate.groupby("sample")[["delta_ad_bc", "delta_ng_mean"]].mean().round(3).to_string())

    # plot
    print("\n=== Plotting joint heatmap ===")
    plot_joint_heatmap(df, FIG_DIR / "01_af_loh_joint_distribution.png")
    print(f"  saved: {FIG_DIR / '01_af_loh_joint_distribution.png'}")


if __name__ == "__main__":
    main()
