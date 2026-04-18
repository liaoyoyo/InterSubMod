#!/usr/bin/env python3
"""
S0.4 / S0.5 / S4.1-S4.4: Self-Phasing Deep Verification + LOH Decomposition
=============================================================================
Continues the phasing causal chain analysis (Phase I = S0.1-S0.3 + S1-S3).

Sections:
  S0.4 - AF-dependent HP imbalance fine-grain analysis
  S0.5 - Self-phasing LOH contribution estimation (statistical simulation)
  S4.1 - LOH.bed vs ISM HP-Ratio LOH concordance (Cohen's kappa)
  S4.2 - LOH flip loci decomposition
  S4.3 - Null model (binomial) vs observed LOH rate
  S4.4 - Self-phasing LOH vs structural LOH separation

Input:
  Master dataset: all_region_rows.tsv.gz (748K rows, 116 cols)
Output:
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
  20260402_phasing_causal_chain/
"""

import os
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.special import comb as _comb

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

# ── Paths ──────────────────────────────────────────────────────────────────
MASTER_DS = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260327_loh_round1_cross_sample_audit/"
    "all_region_rows.tsv.gz"
)
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260402_phasing_causal_chain"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "COLO829",
    "H1437", "H2009", "HCC1937", "HCC1954",
]

# ── Style ──────────────────────────────────────────────────────────────────
COLOR_TP = "#2196F3"
COLOR_FP = "#F44336"
COLOR_PAIRED = "#1565C0"
COLOR_TO = "#E65100"
TRUTH_PALETTE = {"TP": COLOR_TP, "FP": COLOR_FP}

sns.set_theme(style="whitegrid", font_scale=1.05)
plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 13,
})


# ── Helpers ────────────────────────────────────────────────────────────────
def cohens_d(a, b):
    """Cohen's d (pooled SD)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    ma, mb = a.mean(), b.mean()
    sa, sb = a.std(ddof=1), b.std(ddof=1)
    pooled = np.sqrt(((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2))
    return (ma - mb) / pooled if pooled > 0 else 0.0


def save_fig(fig, name, dpi=150):
    p = OUT_DIR / name
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [saved] {p.name} ({p.stat().st_size / 1024:.0f} KB)")


def save_tsv(df, name):
    p = OUT_DIR / name
    df.to_csv(p, sep="\t", index=False)
    print(f"  [saved] {p.name} ({len(df)} rows)")


def map_sample(s):
    """Map sample_label to short sample name if needed."""
    mapping = {
        "HCC1395_HKU_5kHz": "HCC1395",
        "HCC1395_DORADO": "HCC1395_DORADO",
    }
    return mapping.get(s, s)


# ── Load master dataset ───────────────────────────────────────────────────
print("=" * 70)
print("Loading master dataset ...")
print("=" * 70)

COLS = [
    "variant_key", "sample", "sample_label", "mode", "truth_label",
    "HP1FamilyN", "HP2FamilyN", "HP_Ratio", "effective_hp_reads",
    "caller_af", "core_loh_like", "to_loh_bed_hit", "Potential_LOH",
    "effective_hp_bin",
]

df_all = pd.read_csv(MASTER_DS, sep="\t", compression="gzip",
                      usecols=COLS, low_memory=False)
print(f"  Loaded {len(df_all):,} rows x {len(df_all.columns)} cols")

# Normalise sample names
df_all["sample_short"] = df_all["sample_label"].apply(map_sample)
# Verify mapping covers all
for s in df_all["sample_short"].unique():
    if s not in SAMPLE_ORDER:
        print(f"  WARNING: unmapped sample_short={s}")

# Ensure boolean columns
for col in ["core_loh_like", "to_loh_bed_hit", "Potential_LOH"]:
    df_all[col] = df_all[col].astype(bool)

# Derived: HP imbalance
df_all["hp_imbalance"] = (df_all["HP_Ratio"] - 0.5).abs()

print(f"  Samples: {sorted(df_all['sample_short'].unique())}")
print(f"  Modes: {df_all['mode'].value_counts().to_dict()}")
print()

# Sub-frames
df_to = df_all[df_all["mode"] == "to"].copy()
df_paired = df_all[df_all["mode"] == "paired"].copy()


# ═══════════════════════════════════════════════════════════════════════════
# S0.4: AF-Dependent HP Imbalance — Fine-Grain Verification
# ═══════════════════════════════════════════════════════════════════════════
def run_s0_4():
    print("=" * 70)
    print("S0.4: AF-Dependent HP Imbalance (fine-grain, 0.05 step)")
    print("=" * 70)

    # Filter: TO mode, effective_hp_reads >= 30
    dt = df_to[df_to["effective_hp_reads"] >= 30].copy()
    dp = df_paired[df_paired["effective_hp_reads"] >= 30].copy()
    print(f"  TO rows (hp>=30): {len(dt):,}")
    print(f"  Paired rows (hp>=30): {len(dp):,}")

    # AF bins with 0.05 step
    bins_edges = np.arange(0, 1.05, 0.05)
    bins_labels = [f"{bins_edges[i]:.2f}-{bins_edges[i+1]:.2f}" for i in range(len(bins_edges) - 1)]
    dt["af_bin_fine"] = pd.cut(dt["caller_af"], bins=bins_edges, labels=bins_labels, include_lowest=True)
    dp["af_bin_fine"] = pd.cut(dp["caller_af"], bins=bins_edges, labels=bins_labels, include_lowest=True)

    # Aggregate: mean imbalance per (sample, truth_label, af_bin)
    rows = []
    for samp in SAMPLE_ORDER:
        for af_lab in bins_labels:
            af_lo = float(af_lab.split("-")[0])
            af_hi = float(af_lab.split("-")[1])
            af_mid = (af_lo + af_hi) / 2

            to_tp = dt[(dt["sample_short"] == samp) & (dt["truth_label"] == "TP") & (dt["af_bin_fine"] == af_lab)]
            to_fp = dt[(dt["sample_short"] == samp) & (dt["truth_label"] == "FP") & (dt["af_bin_fine"] == af_lab)]
            pa_tp = dp[(dp["sample_short"] == samp) & (dp["truth_label"] == "TP") & (dp["af_bin_fine"] == af_lab)]

            d_tp_vs_fp = cohens_d(to_tp["hp_imbalance"].values, to_fp["hp_imbalance"].values)
            d_to_vs_paired = cohens_d(to_tp["hp_imbalance"].values, pa_tp["hp_imbalance"].values)

            rows.append({
                "sample": samp,
                "af_bin": af_lab,
                "af_mid": af_mid,
                "to_tp_n": len(to_tp),
                "to_tp_mean_imbalance": to_tp["hp_imbalance"].mean() if len(to_tp) > 0 else np.nan,
                "to_fp_n": len(to_fp),
                "to_fp_mean_imbalance": to_fp["hp_imbalance"].mean() if len(to_fp) > 0 else np.nan,
                "paired_tp_n": len(pa_tp),
                "paired_tp_mean_imbalance": pa_tp["hp_imbalance"].mean() if len(pa_tp) > 0 else np.nan,
                "d_TO_TP_vs_FP": d_tp_vs_fp,
                "d_TO_TP_vs_Paired_TP": d_to_vs_paired,
            })

    result = pd.DataFrame(rows)
    save_tsv(result, "s0_4_af_dependence_fine_grain.tsv")

    # ── Plot: 7-sample facet ──
    fig, axes = plt.subplots(2, 4, figsize=(20, 10), sharey=True)
    axes_flat = axes.flatten()

    for idx, samp in enumerate(SAMPLE_ORDER):
        ax = axes_flat[idx]
        sub = result[result["sample"] == samp].copy()

        # TO TP line
        mask_tp = sub["to_tp_n"] >= 5
        ax.plot(sub.loc[mask_tp, "af_mid"], sub.loc[mask_tp, "to_tp_mean_imbalance"],
                "o-", color=COLOR_TP, label="TO TP", markersize=4, linewidth=1.5)
        # TO FP line
        mask_fp = sub["to_fp_n"] >= 5
        ax.plot(sub.loc[mask_fp, "af_mid"], sub.loc[mask_fp, "to_fp_mean_imbalance"],
                "s--", color=COLOR_FP, label="TO FP", markersize=4, linewidth=1.5)
        # Paired TP line
        mask_pa = sub["paired_tp_n"] >= 5
        ax.plot(sub.loc[mask_pa, "af_mid"], sub.loc[mask_pa, "paired_tp_mean_imbalance"],
                "^:", color=COLOR_PAIRED, label="Paired TP", markersize=4, linewidth=1.5)

        ax.set_title(samp, fontweight="bold")
        ax.set_xlabel("AF")
        if idx % 4 == 0:
            ax.set_ylabel("Mean HP Imbalance")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 0.5)
        ax.legend(fontsize=7, loc="upper left")
        ax.axhline(0, color="grey", linewidth=0.5)

    # Hide unused subplot
    axes_flat[7].set_visible(False)

    fig.suptitle("S0.4: HP Imbalance vs AF (0.05 bins, eff_hp >= 30)", fontweight="bold", y=1.01)
    fig.tight_layout()
    save_fig(fig, "s0_4_af_vs_imbalance.png")

    # ── Summary per-AF-bin across all samples ──
    agg_all = []
    for af_lab in bins_labels:
        af_mid = (float(af_lab.split("-")[0]) + float(af_lab.split("-")[1])) / 2
        sub = result[result["af_bin"] == af_lab]
        agg_all.append({
            "af_bin": af_lab,
            "af_mid": af_mid,
            "mean_d_TO_TP_vs_FP": sub["d_TO_TP_vs_FP"].mean(),
            "mean_d_TO_TP_vs_Paired_TP": sub["d_TO_TP_vs_Paired_TP"].mean(),
            "total_to_tp_n": sub["to_tp_n"].sum(),
        })
    agg_df = pd.DataFrame(agg_all)

    print("\n  Cross-sample mean effect sizes by AF bin:")
    print("  " + "-" * 65)
    for _, r in agg_df.iterrows():
        d1 = f"{r['mean_d_TO_TP_vs_FP']:+.3f}" if np.isfinite(r["mean_d_TO_TP_vs_FP"]) else "  N/A"
        d2 = f"{r['mean_d_TO_TP_vs_Paired_TP']:+.3f}" if np.isfinite(r["mean_d_TO_TP_vs_Paired_TP"]) else "  N/A"
        print(f"  AF {r['af_bin']}: d(TO_TP-FP)={d1}  d(TO_TP-PairedTP)={d2}  n={int(r['total_to_tp_n'])}")

    return result


# ═══════════════════════════════════════════════════════════════════════════
# S0.5: Self-Phasing LOH Contribution Estimation
# ═══════════════════════════════════════════════════════════════════════════
def run_s0_5():
    print("\n" + "=" * 70)
    print("S0.5: Self-Phasing LOH Contribution Estimation (Statistical)")
    print("=" * 70)

    # TO TP only, effective_hp_reads >= 30
    dt = df_to[(df_to["truth_label"] == "TP") & (df_to["effective_hp_reads"] >= 30)].copy()
    print(f"  TO TP rows (hp>=30): {len(dt):,}")

    # Compute adjusted HP ratio after removing self-phasing effect
    hp1 = dt["HP1FamilyN"].values.astype(float)
    hp2 = dt["HP2FamilyN"].values.astype(float)
    n_total = hp1 + hp2
    af = dt["caller_af"].values
    n_alt = np.clip(np.round(n_total * af), 0, n_total)

    # Self-phasing hypothesis: all ALT reads are pushed to the majority HP
    majority = np.maximum(hp1, hp2)
    minority = np.minimum(hp1, hp2)

    # Remove ALT reads from majority, add to minority
    adj_majority = np.clip(majority - n_alt, 0, None)
    adj_minority = minority + n_alt
    adj_total = adj_majority + adj_minority

    # Adjusted ratio (majority / total)
    with np.errstate(divide="ignore", invalid="ignore"):
        adj_ratio = np.where(adj_total > 0, adj_majority / adj_total, 0.5)

    # LOH definitions: ratio < 0.1 or > 0.9
    obs_loh = dt["core_loh_like"].values
    adj_loh_like = (adj_ratio < 0.1) | (adj_ratio > 0.9)

    dt["adjusted_hp_ratio"] = adj_ratio
    dt["adjusted_loh_like"] = adj_loh_like
    dt["observed_loh"] = obs_loh

    # LOH disappearance: was LOH, no longer LOH after adjustment
    dt["loh_disappeared"] = obs_loh & (~adj_loh_like)

    # AF bins (coarser for this analysis)
    af_bins_edges = np.arange(0, 1.05, 0.1)
    af_bins_labels = [f"{af_bins_edges[i]:.1f}-{af_bins_edges[i+1]:.1f}" for i in range(len(af_bins_edges) - 1)]
    dt["af_bin_10"] = pd.cut(dt["caller_af"], bins=af_bins_edges, labels=af_bins_labels, include_lowest=True)

    # Per sample x AF bin
    rows = []
    for samp in SAMPLE_ORDER:
        for af_lab in af_bins_labels:
            sub = dt[(dt["sample_short"] == samp) & (dt["af_bin_10"] == af_lab)]
            n = len(sub)
            n_obs_loh = sub["observed_loh"].sum()
            n_adj_loh = sub["adjusted_loh_like"].sum()
            n_disappeared = sub["loh_disappeared"].sum()
            rate = n_disappeared / n_obs_loh if n_obs_loh > 0 else np.nan
            rows.append({
                "sample": samp,
                "af_bin": af_lab,
                "n_total": n,
                "n_observed_loh": int(n_obs_loh),
                "n_adjusted_loh": int(n_adj_loh),
                "n_loh_disappeared": int(n_disappeared),
                "loh_disappearance_rate": rate,
            })

    result = pd.DataFrame(rows)
    save_tsv(result, "s0_5_self_phasing_loh_contribution.tsv")

    # ── Plot: LOH disappearance rate by AF bin per sample ──
    fig, axes = plt.subplots(2, 4, figsize=(20, 10), sharey=True)
    axes_flat = axes.flatten()

    for idx, samp in enumerate(SAMPLE_ORDER):
        ax = axes_flat[idx]
        sub = result[(result["sample"] == samp) & (result["n_observed_loh"] >= 5)]
        af_mids = [(float(x.split("-")[0]) + float(x.split("-")[1])) / 2 for x in sub["af_bin"]]
        ax.bar(af_mids, sub["loh_disappearance_rate"] * 100, width=0.08,
               color=COLOR_TO, alpha=0.7, edgecolor="black", linewidth=0.5)
        ax.axhline(50, color="red", linestyle="--", linewidth=1, label="50% threshold")
        ax.set_title(samp, fontweight="bold")
        ax.set_xlabel("AF")
        if idx % 4 == 0:
            ax.set_ylabel("LOH Disappearance Rate (%)")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 100)
        ax.legend(fontsize=7)

    axes_flat[7].set_visible(False)

    fig.suptitle("S0.5: LOH Disappearance After Removing Self-Phasing (TO TP, hp>=30)",
                 fontweight="bold", y=1.01)
    fig.tight_layout()
    save_fig(fig, "s0_5_loh_disappearance_by_af.png")

    # Summary
    overall = result[result["n_observed_loh"] > 0]
    total_obs = overall["n_observed_loh"].sum()
    total_dis = overall["n_loh_disappeared"].sum()
    print(f"\n  Overall: {total_dis:,}/{total_obs:,} LOH disappeared = "
          f"{100*total_dis/total_obs:.1f}%")

    return result


# ═══════════════════════════════════════════════════════════════════════════
# S4.1: LOH.bed vs ISM HP-Ratio LOH Concordance
# ═══════════════════════════════════════════════════════════════════════════
def run_s4_1():
    print("\n" + "=" * 70)
    print("S4.1: LOH.bed vs ISM core_loh_like Concordance")
    print("=" * 70)

    dt = df_to.copy()
    print(f"  TO rows: {len(dt):,}")

    # Cohen's kappa helper
    def cohens_kappa(a, b):
        """Compute Cohen's kappa between two boolean arrays."""
        a, b = np.asarray(a, dtype=int), np.asarray(b, dtype=int)
        n = len(a)
        if n == 0:
            return np.nan, np.nan, np.nan, np.nan
        # Confusion matrix
        tp = ((a == 1) & (b == 1)).sum()
        tn = ((a == 0) & (b == 0)).sum()
        fp = ((a == 0) & (b == 1)).sum()
        fn = ((a == 1) & (b == 0)).sum()
        po = (tp + tn) / n  # observed agreement
        pe = ((tp + fp) * (tp + fn) + (tn + fn) * (tn + fp)) / (n * n)  # expected
        kappa = (po - pe) / (1 - pe) if (1 - pe) > 0 else 0.0
        return kappa, po, pe, n

    rows = []

    # Overall
    k, po, pe, n = cohens_kappa(dt["to_loh_bed_hit"], dt["core_loh_like"])
    rows.append({
        "sample": "ALL", "stratum": "overall",
        "n": n, "kappa": k, "agreement": po, "expected": pe,
        "loh_bed_rate": dt["to_loh_bed_hit"].mean(),
        "core_loh_rate": dt["core_loh_like"].mean(),
    })

    # By sample
    for samp in SAMPLE_ORDER:
        sub = dt[dt["sample_short"] == samp]
        if len(sub) == 0:
            continue
        k, po, pe, n = cohens_kappa(sub["to_loh_bed_hit"], sub["core_loh_like"])
        rows.append({
            "sample": samp, "stratum": "overall",
            "n": n, "kappa": k, "agreement": po, "expected": pe,
            "loh_bed_rate": sub["to_loh_bed_hit"].mean(),
            "core_loh_rate": sub["core_loh_like"].mean(),
        })

        # By effective_hp_bin
        for eb in sorted(sub["effective_hp_bin"].dropna().unique()):
            sub2 = sub[sub["effective_hp_bin"] == eb]
            if len(sub2) < 10:
                continue
            k2, po2, pe2, n2 = cohens_kappa(sub2["to_loh_bed_hit"], sub2["core_loh_like"])
            rows.append({
                "sample": samp, "stratum": f"hp_bin={eb}",
                "n": n2, "kappa": k2, "agreement": po2, "expected": pe2,
                "loh_bed_rate": sub2["to_loh_bed_hit"].mean(),
                "core_loh_rate": sub2["core_loh_like"].mean(),
            })

    result = pd.DataFrame(rows)
    save_tsv(result, "s4_1_loh_bed_vs_ism_concordance.tsv")

    # Print summary
    print("\n  Concordance summary (overall per sample):")
    for _, r in result[result["stratum"] == "overall"].iterrows():
        print(f"    {r['sample']:20s}  kappa={r['kappa']:.3f}  "
              f"agreement={r['agreement']:.3f}  n={int(r['n']):,}")

    return result


# ═══════════════════════════════════════════════════════════════════════════
# S4.2: LOH Flip Loci Decomposition
# ═══════════════════════════════════════════════════════════════════════════
def run_s4_2():
    print("\n" + "=" * 70)
    print("S4.2: LOH Flip Loci Decomposition")
    print("=" * 70)

    # Build same-locus pairs
    to_cols = ["variant_key", "sample_short", "core_loh_like", "to_loh_bed_hit",
               "caller_af", "effective_hp_reads", "truth_label", "HP_Ratio"]
    pa_cols = ["variant_key", "sample_short", "core_loh_like",
               "caller_af", "effective_hp_reads", "truth_label", "HP_Ratio"]

    dt = df_to[to_cols].rename(columns={
        "core_loh_like": "to_loh", "HP_Ratio": "to_hp_ratio",
        "caller_af": "to_af", "effective_hp_reads": "to_hp_n",
        "truth_label": "to_truth",
    })
    dp = df_paired[pa_cols].rename(columns={
        "core_loh_like": "paired_loh", "HP_Ratio": "paired_hp_ratio",
        "caller_af": "paired_af", "effective_hp_reads": "paired_hp_n",
        "truth_label": "paired_truth",
    })

    # Merge on variant_key + sample_short
    merged = dt.merge(dp, on=["variant_key", "sample_short"], how="inner")
    print(f"  Matched locus pairs: {len(merged):,}")

    # Define flip: TO LOH + Paired non-LOH
    merged["is_flip"] = merged["to_loh"] & (~merged["paired_loh"])

    # Stratification bins
    af_bins = pd.cut(merged["to_af"], bins=[0, 0.2, 0.5, 1.0],
                     labels=["<0.2", "0.2-0.5", ">0.5"], include_lowest=True)
    merged["af_stratum"] = af_bins

    hp_bins = pd.cut(merged["to_hp_n"], bins=[0, 30, 50, 1e6],
                     labels=["<30", "30-49", ">=50"], include_lowest=True)
    merged["hp_stratum"] = hp_bins

    # Cross-tabulation with multiple strata
    rows = []
    for samp in SAMPLE_ORDER:
        ms = merged[merged["sample_short"] == samp]
        if len(ms) == 0:
            continue
        for loh_bed in [True, False]:
            for af_s in ["<0.2", "0.2-0.5", ">0.5"]:
                for hp_s in ["<30", "30-49", ">=50"]:
                    for tl in ["TP", "FP"]:
                        sub = ms[(ms["to_loh_bed_hit"] == loh_bed) &
                                 (ms["af_stratum"] == af_s) &
                                 (ms["hp_stratum"] == hp_s) &
                                 (ms["to_truth"] == tl)]
                        n = len(sub)
                        n_flip = sub["is_flip"].sum()
                        rate = n_flip / n if n > 0 else np.nan
                        rows.append({
                            "sample": samp,
                            "to_loh_bed_hit": loh_bed,
                            "af_stratum": af_s,
                            "hp_stratum": hp_s,
                            "truth_label": tl,
                            "n_total": n,
                            "n_flip": int(n_flip),
                            "flip_rate": rate,
                        })

    result = pd.DataFrame(rows)
    save_tsv(result, "s4_2_loh_flip_decomposition.tsv")

    # Summary: aggregate flip stats
    valid = result[result["n_total"] > 0]
    total_pairs = valid["n_total"].sum()
    total_flips = valid["n_flip"].sum()
    print(f"\n  Total matched pairs in strata: {total_pairs:,}")
    print(f"  Total flips (TO-LOH + Paired-nonLOH): {total_flips:,}")
    print(f"  Overall flip rate: {100*total_flips/total_pairs:.1f}%")

    # Key breakdowns
    print("\n  Flip rate by loh_bed_hit:")
    for lb in [True, False]:
        sub = valid[valid["to_loh_bed_hit"] == lb]
        n, nf = sub["n_total"].sum(), sub["n_flip"].sum()
        print(f"    loh_bed={lb}: {nf:,}/{n:,} = {100*nf/n:.1f}%")

    print("\n  Flip rate by truth_label:")
    for tl in ["TP", "FP"]:
        sub = valid[valid["truth_label"] == tl]
        n, nf = sub["n_total"].sum(), sub["n_flip"].sum()
        print(f"    {tl}: {nf:,}/{n:,} = {100*nf/n:.1f}%")

    return result


# ═══════════════════════════════════════════════════════════════════════════
# S4.3: Null Model (Binomial) vs Observed LOH Rate
# ═══════════════════════════════════════════════════════════════════════════
def run_s4_3():
    print("\n" + "=" * 70)
    print("S4.3: Null Model (Binomial p=0.5) vs Observed LOH Rate")
    print("=" * 70)

    # Theoretical LOH rate under Binomial(n, 0.5)
    # LOH-like = ratio < 0.1 or > 0.9
    # This means X < 0.1*n or X > 0.9*n where X ~ Binomial(n, 0.5)
    ns = np.arange(1, 201)
    null_loh_rate = np.zeros(len(ns))
    for i, n in enumerate(ns):
        lo_thresh = int(np.floor(0.1 * n))
        hi_thresh = int(np.ceil(0.9 * n))
        p_lo = stats.binom.cdf(lo_thresh, n, 0.5)
        p_hi = 1 - stats.binom.cdf(hi_thresh - 1, n, 0.5)
        null_loh_rate[i] = p_lo + p_hi

    null_df = pd.DataFrame({"n": ns, "null_loh_rate": null_loh_rate})

    # Observed LOH rate by effective_hp_reads
    rows = []
    for mode_val, mode_df in [("to", df_to), ("paired", df_paired)]:
        for tl in ["TP", "FP"]:
            sub = mode_df[mode_df["truth_label"] == tl]
            for n in ns:
                s = sub[sub["effective_hp_reads"] == n]
                if len(s) < 5:
                    continue
                obs_rate = s["core_loh_like"].mean()
                rows.append({
                    "n_reads": n,
                    "mode": mode_val,
                    "truth_label": tl,
                    "observed_loh_rate": obs_rate,
                    "count": len(s),
                })

    obs_df = pd.DataFrame(rows)

    # Merge with null
    obs_df = obs_df.merge(null_df, left_on="n_reads", right_on="n", how="left")
    obs_df["excess_loh_rate"] = obs_df["observed_loh_rate"] - obs_df["null_loh_rate"]

    # Also compute binned version for smoother curves
    bin_edges = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200]
    bin_labels = [f"{bin_edges[i]}-{bin_edges[i+1]}" for i in range(len(bin_edges) - 1)]

    rows_binned = []
    for mode_val, mode_df in [("to", df_to), ("paired", df_paired)]:
        for tl in ["TP", "FP"]:
            sub = mode_df[mode_df["truth_label"] == tl]
            for i in range(len(bin_edges) - 1):
                lo, hi = bin_edges[i], bin_edges[i + 1]
                s = sub[(sub["effective_hp_reads"] >= lo) & (sub["effective_hp_reads"] < hi)]
                if len(s) < 20:
                    continue
                obs_rate = s["core_loh_like"].mean()
                n_mid = (lo + hi) / 2
                # Null rate at midpoint
                lo_t = int(np.floor(0.1 * n_mid))
                hi_t = int(np.ceil(0.9 * n_mid))
                p_lo = stats.binom.cdf(lo_t, int(n_mid), 0.5)
                p_hi = 1 - stats.binom.cdf(hi_t - 1, int(n_mid), 0.5)
                null_rate = p_lo + p_hi
                rows_binned.append({
                    "n_bin": bin_labels[i],
                    "n_mid": n_mid,
                    "mode": mode_val,
                    "truth_label": tl,
                    "observed_loh_rate": obs_rate,
                    "null_loh_rate": null_rate,
                    "excess": obs_rate - null_rate,
                    "count": len(s),
                })

    binned_df = pd.DataFrame(rows_binned)
    combined = pd.concat([obs_df.assign(source="per_n"), binned_df.assign(source="binned")])

    save_tsv(obs_df, "s4_3_null_model_loh_rate.tsv")

    # ── Plot: Observed vs Null LOH rate ──
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Panel A: line plot by n_mid (binned)
    ax = axes[0]
    for mode_val, ls in [("to", "-"), ("paired", "--")]:
        for tl, c in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            sub = binned_df[(binned_df["mode"] == mode_val) & (binned_df["truth_label"] == tl)]
            if len(sub) == 0:
                continue
            ax.plot(sub["n_mid"], sub["observed_loh_rate"] * 100, f"o{ls}",
                    color=c, markersize=5, linewidth=1.5,
                    label=f"{mode_val} {tl}",
                    alpha=0.8 if mode_val == "to" else 0.5)

    # Null model curve
    ax.plot(ns, null_loh_rate * 100, "k-", linewidth=2, label="Null (Binomial p=0.5)")
    ax.set_xlabel("effective_hp_reads (bin midpoint)")
    ax.set_ylabel("LOH Rate (%)")
    ax.set_title("Observed vs Null LOH Rate by Read Count")
    ax.legend(fontsize=8)
    ax.set_xlim(10, 200)
    ax.set_ylim(0, None)

    # Panel B: excess LOH (observed - null)
    ax = axes[1]
    for mode_val, ls in [("to", "-"), ("paired", "--")]:
        for tl, c in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            sub = binned_df[(binned_df["mode"] == mode_val) & (binned_df["truth_label"] == tl)]
            if len(sub) == 0:
                continue
            ax.plot(sub["n_mid"], sub["excess"] * 100, f"o{ls}",
                    color=c, markersize=5, linewidth=1.5,
                    label=f"{mode_val} {tl}",
                    alpha=0.8 if mode_val == "to" else 0.5)

    ax.axhline(0, color="black", linewidth=1)
    ax.set_xlabel("effective_hp_reads (bin midpoint)")
    ax.set_ylabel("Excess LOH Rate (obs - null, %)")
    ax.set_title("Excess LOH Rate Above Null")
    ax.legend(fontsize=8)
    ax.set_xlim(10, 200)

    fig.suptitle("S4.3: Null Model (Binomial p=0.5) vs Observed LOH Rate", fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "s4_3_observed_vs_null_loh.png")

    # Self-phasing LOH estimate
    print("\n  Self-phasing LOH estimate (binned):")
    print("  " + "-" * 70)
    for _, r in binned_df.iterrows():
        print(f"    {r['n_bin']:10s}  {r['mode']:8s}  {r['truth_label']}  "
              f"obs={r['observed_loh_rate']:.3f}  null={r['null_loh_rate']:.3f}  "
              f"excess={r['excess']:+.3f}  n={int(r['count']):,}")

    # Compute self-phasing estimate: excess(TO TP) - excess(Paired TP)
    print("\n  Self-phasing LOH contribution = excess(TO_TP) - excess(Paired_TP):")
    to_tp_b = binned_df[(binned_df["mode"] == "to") & (binned_df["truth_label"] == "TP")]
    pa_tp_b = binned_df[(binned_df["mode"] == "paired") & (binned_df["truth_label"] == "TP")]

    if len(to_tp_b) > 0 and len(pa_tp_b) > 0:
        # Merge on n_bin
        sp = to_tp_b[["n_bin", "n_mid", "excess"]].rename(columns={"excess": "to_tp_excess"})
        sp = sp.merge(pa_tp_b[["n_bin", "excess"]].rename(columns={"excess": "pa_tp_excess"}),
                      on="n_bin", how="inner")
        sp["self_phasing_loh"] = sp["to_tp_excess"] - sp["pa_tp_excess"]
        for _, r in sp.iterrows():
            print(f"    n_bin={r['n_bin']:10s}  TO_excess={r['to_tp_excess']:+.4f}  "
                  f"Pa_excess={r['pa_tp_excess']:+.4f}  self_phasing={r['self_phasing_loh']:+.4f}")

    return binned_df


# ═══════════════════════════════════════════════════════════════════════════
# S4.4: Self-Phasing LOH vs Structural LOH Separation
# ═══════════════════════════════════════════════════════════════════════════
def run_s4_4():
    print("\n" + "=" * 70)
    print("S4.4: Self-Phasing LOH vs Structural LOH Separation")
    print("=" * 70)

    # Build same-locus pairs
    to_cols = ["variant_key", "sample_short", "core_loh_like", "truth_label"]
    pa_cols = ["variant_key", "sample_short", "core_loh_like"]

    dt = df_to[to_cols].rename(columns={"core_loh_like": "to_loh"})
    dp = df_paired[pa_cols].rename(columns={"core_loh_like": "paired_loh"})

    merged = dt.merge(dp, on=["variant_key", "sample_short"], how="inner")
    print(f"  Matched locus pairs: {len(merged):,}")

    # Three categories:
    # 1. Structural LOH: both TO and Paired are LOH
    # 2. Self-phasing LOH: TO is LOH but Paired is not
    # 3. Shared non-LOH: neither is LOH (or only Paired is LOH)
    merged["loh_type"] = "shared_non_LOH"
    merged.loc[merged["to_loh"] & merged["paired_loh"], "loh_type"] = "structural_LOH"
    merged.loc[merged["to_loh"] & (~merged["paired_loh"]), "loh_type"] = "self_phasing_LOH"
    merged.loc[(~merged["to_loh"]) & merged["paired_loh"], "loh_type"] = "paired_only_LOH"

    # Per sample x truth_label
    rows = []
    for samp in SAMPLE_ORDER:
        for tl in ["TP", "FP"]:
            sub = merged[(merged["sample_short"] == samp) & (merged["truth_label"] == tl)]
            n = len(sub)
            if n == 0:
                continue

            n_to_loh = sub["to_loh"].sum()
            n_structural = (sub["loh_type"] == "structural_LOH").sum()
            n_self_phasing = (sub["loh_type"] == "self_phasing_LOH").sum()
            n_paired_only = (sub["loh_type"] == "paired_only_LOH").sum()
            n_shared_nonloh = (sub["loh_type"] == "shared_non_LOH").sum()

            sp_pct_of_to_loh = n_self_phasing / n_to_loh if n_to_loh > 0 else np.nan

            rows.append({
                "sample": samp,
                "truth_label": tl,
                "n_pairs": n,
                "n_to_loh": int(n_to_loh),
                "n_structural_LOH": int(n_structural),
                "n_self_phasing_LOH": int(n_self_phasing),
                "n_paired_only_LOH": int(n_paired_only),
                "n_shared_non_LOH": int(n_shared_nonloh),
                "structural_frac": n_structural / n,
                "self_phasing_frac": n_self_phasing / n,
                "paired_only_frac": n_paired_only / n,
                "shared_non_loh_frac": n_shared_nonloh / n,
                "self_phasing_pct_of_to_loh": sp_pct_of_to_loh,
            })

    result = pd.DataFrame(rows)
    save_tsv(result, "s4_4_loh_type_decomposition.tsv")

    # Summary
    print("\n  LOH Type Decomposition Summary:")
    print("  " + "-" * 90)
    for _, r in result.iterrows():
        print(f"    {r['sample']:20s} {r['truth_label']}  "
              f"structural={r['n_structural_LOH']:5d} ({100*r['structural_frac']:.1f}%)  "
              f"self-phase={r['n_self_phasing_LOH']:5d} ({100*r['self_phasing_frac']:.1f}%)  "
              f"paired-only={r['n_paired_only_LOH']:5d} ({100*r['paired_only_frac']:.1f}%)  "
              f"SP/TO_LOH={100*r['self_phasing_pct_of_to_loh']:.1f}%")

    # Overall aggregation
    tp_all = result[result["truth_label"] == "TP"]
    total_to_loh = tp_all["n_to_loh"].sum()
    total_sp = tp_all["n_self_phasing_LOH"].sum()
    total_str = tp_all["n_structural_LOH"].sum()
    print(f"\n  Overall TP: TO LOH={total_to_loh:,}  "
          f"structural={total_str:,} ({100*total_str/total_to_loh:.1f}%)  "
          f"self-phasing={total_sp:,} ({100*total_sp/total_to_loh:.1f}%)")

    fp_all = result[result["truth_label"] == "FP"]
    if len(fp_all) > 0 and fp_all["n_to_loh"].sum() > 0:
        total_to_loh_fp = fp_all["n_to_loh"].sum()
        total_sp_fp = fp_all["n_self_phasing_LOH"].sum()
        total_str_fp = fp_all["n_structural_LOH"].sum()
        print(f"  Overall FP: TO LOH={total_to_loh_fp:,}  "
              f"structural={total_str_fp:,} ({100*total_str_fp/total_to_loh_fp:.1f}%)  "
              f"self-phasing={total_sp_fp:,} ({100*total_sp_fp/total_to_loh_fp:.1f}%)")

    return result


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("S0.4 / S0.5 / S4.1-S4.4: Self-Phasing Deep + LOH Decomposition")
    print(f"Output: {OUT_DIR}")
    print("=" * 70, "\n")

    r04 = run_s0_4()
    r05 = run_s0_5()
    r41 = run_s4_1()
    r42 = run_s4_2()
    r43 = run_s4_3()
    r44 = run_s4_4()

    # ── SUMMARY ────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    # S0.4 summary
    print("\n[S0.4] AF-Dependent HP Imbalance:")
    if r04 is not None:
        # Aggregate across samples for high-AF and low-AF
        low_af = r04[(r04["af_mid"] <= 0.1)]
        mid_af = r04[(r04["af_mid"] > 0.2) & (r04["af_mid"] <= 0.5)]
        high_af = r04[(r04["af_mid"] > 0.5)]
        for label, sub in [("AF<=0.10", low_af), ("AF 0.20-0.50", mid_af), ("AF>0.50", high_af)]:
            d1 = sub["d_TO_TP_vs_FP"].mean()
            d2 = sub["d_TO_TP_vs_Paired_TP"].mean()
            n = sub["to_tp_n"].sum()
            print(f"  {label}: d(TO_TP-FP)={d1:+.3f}  d(TO_TP-PairedTP)={d2:+.3f}  total_n={int(n):,}")

    # S0.5 summary
    print("\n[S0.5] Self-Phasing LOH Contribution:")
    if r05 is not None:
        valid05 = r05[r05["n_observed_loh"] > 0]
        total_obs = valid05["n_observed_loh"].sum()
        total_dis = valid05["n_loh_disappeared"].sum()
        print(f"  LOH disappeared after removing self-phasing: {total_dis:,}/{total_obs:,} = "
              f"{100*total_dis/total_obs:.1f}%")
        if 100 * total_dis / total_obs > 50:
            print("  -> CONFIRMED: Self-phasing is the dominant LOH source (>50% disappearance)")
        else:
            print("  -> Self-phasing is NOT the dominant LOH source (<50% disappearance)")

    # S4.1 summary
    print("\n[S4.1] LOH.bed vs ISM Concordance:")
    if r41 is not None:
        overall = r41[r41["sample"] == "ALL"].iloc[0]
        print(f"  Overall kappa = {overall['kappa']:.3f} (agreement={overall['agreement']:.3f})")
        if overall["kappa"] < 0.4:
            print("  -> Low agreement: LOH.bed and ISM HP-ratio LOH are measuring different things")
        elif overall["kappa"] < 0.6:
            print("  -> Moderate agreement")
        else:
            print("  -> Substantial agreement")

    # S4.2 summary
    print("\n[S4.2] LOH Flip Loci Decomposition:")
    if r42 is not None:
        valid42 = r42[r42["n_total"] > 0]
        total_n = valid42["n_total"].sum()
        total_f = valid42["n_flip"].sum()
        print(f"  Total flip (TO-LOH + Paired-nonLOH): {total_f:,}/{total_n:,} = "
              f"{100*total_f/total_n:.1f}%")

        # TP vs FP flip
        for tl in ["TP", "FP"]:
            sub = valid42[valid42["truth_label"] == tl]
            n, nf = sub["n_total"].sum(), sub["n_flip"].sum()
            if n > 0:
                print(f"  {tl} flip rate: {nf:,}/{n:,} = {100*nf/n:.1f}%")

    # S4.3 summary
    print("\n[S4.3] Null Model LOH Rate:")
    # Already printed in run_s4_3

    # S4.4 summary
    print("\n[S4.4] Self-Phasing vs Structural LOH:")
    if r44 is not None:
        tp44 = r44[r44["truth_label"] == "TP"]
        total_to_loh = tp44["n_to_loh"].sum()
        total_sp = tp44["n_self_phasing_LOH"].sum()
        total_str = tp44["n_structural_LOH"].sum()
        if total_to_loh > 0:
            print(f"  TP: self-phasing = {total_sp:,}/{total_to_loh:,} TO LOH = "
                  f"{100*total_sp/total_to_loh:.1f}%")
            print(f"  TP: structural   = {total_str:,}/{total_to_loh:,} TO LOH = "
                  f"{100*total_str/total_to_loh:.1f}%")

        fp44 = r44[r44["truth_label"] == "FP"]
        total_to_loh_fp = fp44["n_to_loh"].sum()
        total_sp_fp = fp44["n_self_phasing_LOH"].sum()
        if total_to_loh_fp > 0:
            print(f"  FP: self-phasing = {total_sp_fp:,}/{total_to_loh_fp:,} TO LOH = "
                  f"{100*total_sp_fp/total_to_loh_fp:.1f}%")

    print("\n" + "=" * 70)
    print(f"All outputs saved to: {OUT_DIR}")
    print("=" * 70)
