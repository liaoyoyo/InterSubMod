#!/usr/bin/env python3
"""
Part B 合併分析：B.1-3 + B.2-5 + B.2-2
==========================================
B.1-3：7/7 effect size per sample — Cohen's d、rank-biserial r
B.2-5：cnLOH (CN≈2) vs deletion-LOH (CN<0.75) 分層
B.2-2：Coverage_Multiple 分佈驗證（bimodality、樣本間一致性）

輸入：
  - all_region_rows.tsv.gz

輸出：
  - data/b1_3_effect_size_per_sample.tsv (HPFineNGroups + LOH AF×NG)
  - data/b2_5_cn_stratified_effects.tsv (cnLOH vs deletion-LOH)
  - data/b2_2_coverage_multiple_distribution.tsv
  - figures/01_b1_3_effect_size_forest.png
  - figures/02_b2_5_cn_stratified_forest.png
  - figures/03_b2_2_cov_multiple_dist.png
"""

import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUTDIR = ROOT / "research/partB_effect_size_cn_stratification"
DATA = OUTDIR / "data"
FIG = OUTDIR / "figures"
DATA.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]


def cohens_d(a, b):
    """Cohen's d for two groups."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2 + 1e-12)
    return (a.mean() - b.mean()) / pooled


def rank_biserial(a, b):
    """Rank-biserial r from Mann-Whitney U."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    u, _ = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u / (len(a) * len(b)) - 1


def load_master():
    print("Loading master dataset...")
    cols = ["sample", "mode", "truth_label", "NumReads", "HPFineNGroups",
            "caller_af", "Coverage_Multiple", "to_loh_bed_hit",
            "AlleleDelta", "HPFineF"]
    df = pd.read_csv(MASTER, sep="\t", usecols=cols, low_memory=False)
    print(f"  rows: {len(df):,}")
    return df


def b1_3_effect_size(df):
    """B.1-3: Per-sample effect size for HPFineNGroups (N=4 vs <4) in TO NonLOH NR>=40"""
    print("\n[B.1-3] HPFineNGroups effect size per sample (TO NonLOH NR>=40)")
    rows = []
    for sample in SAMPLE_ORDER:
        sub = df[(df["sample"] == sample) &
                 (df["mode"] == "to") &
                 (df["to_loh_bed_hit"].astype(str).str.lower() == "false") &
                 (df["NumReads"] >= 40)]
        n4 = sub[sub["HPFineNGroups"] == 4]
        nlt4 = sub[sub["HPFineNGroups"] < 4]
        if len(n4) < 10 or len(nlt4) < 10:
            continue
        # Use TP rate as binary outcome — Cohen's d with variance from Bernoulli
        tp_n4 = (n4["truth_label"] == "TP").astype(int)
        tp_nlt = (nlt4["truth_label"] == "TP").astype(int)
        d = cohens_d(tp_n4, tp_nlt)
        r = rank_biserial(tp_n4, tp_nlt)
        delta_tp = tp_n4.mean() - tp_nlt.mean()
        rows.append(dict(
            sample=sample,
            n_ngroups4=len(n4),
            n_ngroups_lt4=len(nlt4),
            tp_rate_n4=tp_n4.mean(),
            tp_rate_nlt4=tp_nlt.mean(),
            delta_tp_rate=delta_tp,
            cohens_d=d,
            rank_biserial_r=r,
        ))
    df_out = pd.DataFrame(rows)
    df_out.to_csv(DATA / "b1_3_effect_size_per_sample.tsv", sep="\t", index=False)
    print(df_out.to_string(index=False))
    return df_out


def b2_5_cn_stratified(df):
    """B.2-5: LOH AF × Methylation by CN tier (deletion-LOH CN<0.75 vs cnLOH CN 1.5-2.5)"""
    print("\n[B.2-5] LOH AF × HPFineNGroups stratified by CN")

    df_loh = df[(df["to_loh_bed_hit"].astype(str).str.lower() == "true") &
                (df["mode"] == "to") &
                (df["truth_label"] == "TP")].copy()

    # CN tiers
    df_loh["cn_tier"] = pd.cut(df_loh["Coverage_Multiple"],
                               bins=[-np.inf, 0.75, 1.25, 1.75, np.inf],
                               labels=["CN1_deletion", "CN2_near", "CN3_cnLOH_high", "CN4plus"])

    # Intermediate vs Extreme AF
    df_loh["af_class"] = np.where(
        ((df_loh["caller_af"] >= 0.1) & (df_loh["caller_af"] <= 0.4)) |
        ((df_loh["caller_af"] >= 0.6) & (df_loh["caller_af"] <= 0.9)),
        "Intermediate", "Extreme"
    )

    rows = []
    for sample in SAMPLE_ORDER:
        for cn in ["CN1_deletion", "CN2_near", "CN3_cnLOH_high", "CN4plus"]:
            sub = df_loh[(df_loh["sample"] == sample) & (df_loh["cn_tier"] == cn)]
            if len(sub) < 30:
                rows.append(dict(sample=sample, cn_tier=cn,
                                 n_extreme=np.nan, n_inter=np.nan,
                                 ng_extreme=np.nan, ng_inter=np.nan,
                                 delta_ng=np.nan, mw_p=np.nan, r=np.nan,
                                 n_total=len(sub)))
                continue
            ext = sub[sub["af_class"] == "Extreme"]["HPFineNGroups"].dropna()
            inter = sub[sub["af_class"] == "Intermediate"]["HPFineNGroups"].dropna()
            if len(ext) < 10 or len(inter) < 10:
                continue
            try:
                _, p = scipy_stats.mannwhitneyu(inter, ext, alternative="two-sided")
            except ValueError:
                p = np.nan
            r = rank_biserial(inter, ext)
            rows.append(dict(
                sample=sample, cn_tier=cn,
                n_extreme=len(ext), n_inter=len(inter),
                ng_extreme=ext.mean(), ng_inter=inter.mean(),
                delta_ng=inter.mean() - ext.mean(),
                mw_p=p, r=r,
                n_total=len(sub),
            ))
    df_out = pd.DataFrame(rows)
    df_out.to_csv(DATA / "b2_5_cn_stratified_effects.tsv", sep="\t", index=False)
    print(df_out.to_string(index=False))

    # Summary by CN tier
    print("\n  Summary by CN tier (means across samples where |Δ|>0):")
    summary = df_out.dropna(subset=["delta_ng"]).groupby("cn_tier", observed=True).agg(
        mean_delta=("delta_ng", "mean"),
        median_delta=("delta_ng", "median"),
        min_delta=("delta_ng", "min"),
        max_delta=("delta_ng", "max"),
        samples_positive=("delta_ng", lambda x: (x > 0).sum()),
        samples_total=("delta_ng", "count"),
    )
    print(summary.to_string())
    return df_out


def b2_2_coverage_dist(df):
    """B.2-2: Coverage_Multiple distribution analysis per sample."""
    print("\n[B.2-2] Coverage_Multiple 分佈分析")

    df_loh = df[(df["to_loh_bed_hit"].astype(str).str.lower() == "true") &
                (df["mode"] == "to")].copy()

    rows = []
    for sample in SAMPLE_ORDER:
        sub = df_loh[df_loh["sample"] == sample]
        cm = sub["Coverage_Multiple"].dropna()
        if len(cm) < 50:
            continue
        # Bimodality via KDE peaks: coarse — count samples in CN-tier ranges
        cn1 = (cm < 0.75).sum()
        cn2 = ((cm >= 0.75) & (cm < 1.25)).sum()
        cn3 = ((cm >= 1.25) & (cm < 1.75)).sum()
        cn4 = (cm >= 1.75).sum()
        total = len(cm)
        # Check for bimodal structure: compare peak separation
        median = np.median(cm)
        q1, q3 = np.percentile(cm, [25, 75])
        iqr = q3 - q1
        # Dip statistic proxy: count modes via histogram
        hist, edges = np.histogram(cm, bins=20)
        # Find local maxima
        peaks = 0
        for i in range(1, len(hist) - 1):
            if hist[i] > hist[i-1] and hist[i] > hist[i+1] and hist[i] > total * 0.02:
                peaks += 1
        rows.append(dict(
            sample=sample,
            n=total,
            median_cm=median,
            iqr_cm=iqr,
            pct_CN1_deletion=100 * cn1 / total,
            pct_CN2_near=100 * cn2 / total,
            pct_CN3_cnLOH=100 * cn3 / total,
            pct_CN4plus=100 * cn4 / total,
            n_histogram_peaks=peaks,
        ))
    df_out = pd.DataFrame(rows)
    df_out.to_csv(DATA / "b2_2_coverage_multiple_distribution.tsv", sep="\t", index=False)
    print(df_out.to_string(index=False))

    # Plot distributions
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()
    for i, sample in enumerate(SAMPLE_ORDER):
        sub = df_loh[df_loh["sample"] == sample]
        cm = sub["Coverage_Multiple"].dropna()
        if len(cm) < 50:
            axes[i].text(0.5, 0.5, f"{sample}\nn={len(cm)}", transform=axes[i].transAxes, ha="center")
            continue
        axes[i].hist(cm.clip(0, 3), bins=40, color="steelblue", edgecolor="none")
        axes[i].axvline(0.75, color="red", ls="--", lw=1, label="CN1/CN2")
        axes[i].axvline(1.25, color="red", ls="--", lw=1)
        axes[i].axvline(1.75, color="red", ls="--", lw=1)
        axes[i].set_title(f"{sample} (n={len(cm):,})")
        axes[i].set_xlabel("Coverage_Multiple")
    axes[-1].axis("off")
    plt.tight_layout()
    plt.savefig(FIG / "03_b2_2_cov_multiple_dist.png", dpi=150)
    plt.close()
    return df_out


def plot_b1_3_forest(df, out_path):
    fig, ax = plt.subplots(figsize=(10, 5))
    y = np.arange(len(df))
    colors = ["red" if r < 0 or abs(r) < 0.1 else "steelblue" for r in df["rank_biserial_r"]]
    ax.barh(y, df["rank_biserial_r"], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(df["sample"])
    ax.set_xlabel("Rank-biserial r (NGroups=4 vs <4, TP rate)")
    ax.axvline(0, color="black", ls=":", lw=1)
    ax.axvline(0.1, color="gray", ls="--", lw=1, label="small effect threshold")
    ax.axvline(0.3, color="gray", ls="--", lw=1, label="medium effect threshold")
    for i, (r, delta) in enumerate(zip(df["rank_biserial_r"], df["delta_tp_rate"])):
        ax.text(r + 0.01 if r > 0 else r - 0.01, i,
                f"r={r:.2f} (ΔTP={delta*100:+.1f}pp)",
                va="center", ha="left" if r > 0 else "right", fontsize=8)
    ax.legend()
    ax.set_title("B.1-3: Per-sample effect size (HPFineNGroups=4 vs <4, TO NonLOH NR≥40)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def plot_b2_5_forest(df, out_path):
    df = df.dropna(subset=["delta_ng"]).copy()
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, cn in zip(axes, ["CN1_deletion", "CN2_near", "CN3_cnLOH_high"]):
        sub = df[df["cn_tier"] == cn].set_index("sample").reindex(SAMPLE_ORDER).reset_index()
        sub = sub.dropna(subset=["delta_ng"])
        if len(sub) == 0:
            ax.set_title(f"{cn}\n(no data)")
            continue
        y = np.arange(len(sub))
        colors = ["red" if s == "HCC1954" else "steelblue" for s in sub["sample"]]
        ax.barh(y, sub["delta_ng"], color=colors)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["sample"])
        ax.axvline(0, color="black", ls=":", lw=1)
        for i, (d, r, n) in enumerate(zip(sub["delta_ng"], sub["r"], sub["n_total"])):
            ax.text(d, i, f"Δ={d:+.2f} r={r:+.2f} n={n}",
                    va="center", ha="left" if d > 0 else "right", fontsize=7)
        ax.set_title(f"{cn}")
        ax.set_xlabel("Δ mean_NGroups (Inter − Ext AF)")
    fig.suptitle("B.2-5: LOH AF × HPFineNGroups by CN tier (TO, TP only)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def main():
    print("=== Part B 合併：B.1-3 + B.2-5 + B.2-2 ===")

    df = load_master()

    print("\n" + "="*60)
    b1_3 = b1_3_effect_size(df)
    plot_b1_3_forest(b1_3, FIG / "01_b1_3_effect_size_forest.png")

    print("\n" + "="*60)
    b2_5 = b2_5_cn_stratified(df)
    plot_b2_5_forest(b2_5, FIG / "02_b2_5_cn_stratified_forest.png")

    print("\n" + "="*60)
    b2_2 = b2_2_coverage_dist(df)

    # Final verdict per sub-analysis
    print("\n=== FINAL VERDICTS ===")

    # B.1-3
    print("\n[B.1-3]")
    min_r = b1_3["rank_biserial_r"].min()
    max_r = b1_3["rank_biserial_r"].max()
    n_small = (b1_3["rank_biserial_r"].abs() < 0.1).sum()
    n_medium = ((b1_3["rank_biserial_r"].abs() >= 0.1) & (b1_3["rank_biserial_r"].abs() < 0.3)).sum()
    n_large = (b1_3["rank_biserial_r"].abs() >= 0.3).sum()
    print(f"  r range: [{min_r:.3f}, {max_r:.3f}]")
    print(f"  effect size tier: {n_small} negligible / {n_medium} small / {n_large} medium+")
    if n_small >= 3:
        print("  → 「7/7 一致」強度精確化需要：多樣本 negligible (r<0.1)")
    if (b1_3["rank_biserial_r"] < 0).any():
        print("  → 有負向樣本存在")

    # B.2-5
    print("\n[B.2-5]")
    for cn in ["CN1_deletion", "CN2_near", "CN3_cnLOH_high"]:
        sub = b2_5[b2_5["cn_tier"] == cn].dropna(subset=["delta_ng"])
        if len(sub) == 0:
            continue
        n_pos = (sub["delta_ng"] > 0).sum()
        n_total = len(sub)
        median_d = sub["delta_ng"].median()
        print(f"  {cn}: {n_pos}/{n_total} POS, median Δ_NG={median_d:+.3f}")

    # B.2-2
    print("\n[B.2-2]")
    n_bimodal = (b2_2["n_histogram_peaks"] >= 2).sum()
    print(f"  Samples with bimodal-ish Coverage_Multiple dist: {n_bimodal}/{len(b2_2)}")
    print("  Per-sample CN1/CN2/CN3/CN4+ 分佈：")
    for _, r in b2_2.iterrows():
        print(f"    {r['sample']}: CN1={r['pct_CN1_deletion']:.1f}%, "
              f"CN2={r['pct_CN2_near']:.1f}%, CN3={r['pct_CN3_cnLOH']:.1f}%, "
              f"CN4+={r['pct_CN4plus']:.1f}%")


if __name__ == "__main__":
    main()
