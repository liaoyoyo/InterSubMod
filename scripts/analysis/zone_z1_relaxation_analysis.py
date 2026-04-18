#!/usr/bin/env python3
"""
Zone Z1 Relaxation Analysis — Explore coverage vs TP rate tradeoff
under different Z1 definitions.

Goal: Z1 in its current form (LOH + Intermediate AF + NGroups>=3)
covers only 0.04% of TO regions. Need to find a relaxed definition
that maintains high TP rate while improving coverage.

Tested variants:
  Z1a: LOH + Intermediate AF (drop NGroups requirement)
  Z1b: LOH + Intermediate AF + NGroups >= 2
  Z1c: LOH + NGroups >= 3 (drop AF requirement)
  Z1d: LOH + NGroups >= 2 (drop AF requirement, relax NGroups)
  Z1e: to_loh_bed_hit only (simplest — any LOH.bed overlap in TO)
  Z1f: LOH + NR >= 40 (coverage-based)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/" \
         "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/zone_aware_validation"
FIG_DIR = os.path.join(OUT_DIR, "figures")

SAMPLE_ORDER = ["H2009", "H1437", "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "COLO829"]

# ── Load ───────────────────────────────────────────────────────────
print("Loading...")
cols = ["sample", "mode", "truth_label", "core_loh_like", "to_loh_bed_hit",
        "caller_af", "HPFineNGroups", "NumReads"]
df = pd.read_csv(MASTER, sep='\t', usecols=cols, low_memory=False)
print(f"Loaded {len(df):,} rows")

df["caller_af"] = pd.to_numeric(df["caller_af"], errors="coerce")
df["HPFineNGroups"] = pd.to_numeric(df["HPFineNGroups"], errors="coerce")
df["NumReads"] = pd.to_numeric(df["NumReads"], errors="coerce")
for col in ["core_loh_like", "to_loh_bed_hit"]:
    df[col] = df[col].map({"True": True, "False": False, True: True, False: False}).fillna(False)
df["is_tp"] = df["truth_label"] == "TP"
df["loh_flag"] = np.where(df["mode"] == "paired", df["core_loh_like"], df["to_loh_bed_hit"])

def is_intermediate_af(af):
    if pd.isna(af): return False
    return (0.1 <= af <= 0.4) or (0.6 <= af <= 0.9)

df["af_intermediate"] = df["caller_af"].apply(is_intermediate_af)

# ── Define Z1 Variants ────────────────────────────────────────────
z1_variants = {
    "Z1_original": lambda d: d["loh_flag"] & d["af_intermediate"] & (d["HPFineNGroups"] >= 3),
    "Z1a_LOH+IntAF": lambda d: d["loh_flag"] & d["af_intermediate"],
    "Z1b_LOH+IntAF+NG2": lambda d: d["loh_flag"] & d["af_intermediate"] & (d["HPFineNGroups"] >= 2),
    "Z1c_LOH+NG3": lambda d: d["loh_flag"] & (d["HPFineNGroups"] >= 3),
    "Z1d_LOH+NG2": lambda d: d["loh_flag"] & (d["HPFineNGroups"] >= 2),
    "Z1e_LOH_only": lambda d: d["loh_flag"],
    "Z1f_LOH+NR40": lambda d: d["loh_flag"] & (d["NumReads"] >= 40),
}

# ── Evaluate ──────────────────────────────────────────────────────
results = []
for mode in ["paired", "to"]:
    dm = df[df["mode"] == mode]
    global_tp_rate = dm["is_tp"].mean()
    global_n = len(dm)

    for vname, vfunc in z1_variants.items():
        mask = vfunc(dm)
        n = mask.sum()
        tp = dm.loc[mask, "is_tp"].sum()
        fp = n - tp
        tp_rate = tp / n if n > 0 else np.nan
        pct = n / global_n * 100

        # Per-sample stats
        sample_tp_rates = []
        sample_ns = []
        n_positive = 0
        n_sig = 0
        for sample in SAMPLE_ORDER:
            sm = dm[dm["sample"] == sample]
            if len(sm) == 0: continue
            z_mask = vfunc(sm)
            sn = z_mask.sum()
            if sn == 0:
                sample_tp_rates.append(np.nan)
                sample_ns.append(0)
                continue
            stp = sm.loc[z_mask, "is_tp"].sum()
            sfp = sn - stp
            stp_rate = stp / sn
            sample_tp_rates.append(stp_rate)
            sample_ns.append(sn)

            s_global_tp = sm["is_tp"].mean()
            if stp_rate > s_global_tp:
                n_positive += 1
            # Fisher test
            non_z = sm[~z_mask]
            nzt = non_z["is_tp"].sum()
            nzf = len(non_z) - nzt
            if sn > 0 and len(non_z) > 0:
                _, pval = stats.fisher_exact([[stp, sfp], [nzt, nzf]], alternative='greater')
                if pval < 0.05:
                    n_sig += 1

        mean_sample_tp = np.nanmean(sample_tp_rates)
        min_sample_tp = np.nanmin(sample_tp_rates) if any(not np.isnan(x) for x in sample_tp_rates) else np.nan

        results.append({
            "mode": mode, "variant": vname,
            "n": n, "pct": pct, "tp": tp, "fp": fp,
            "tp_rate": tp_rate, "global_tp_rate": global_tp_rate,
            "delta": tp_rate - global_tp_rate if not np.isnan(tp_rate) else np.nan,
            "mean_sample_tp": mean_sample_tp,
            "min_sample_tp": min_sample_tp,
            "n_positive": n_positive, "n_sig": n_sig,
            "n_samples_with_data": sum(1 for n in sample_ns if n > 0)
        })

rdf = pd.DataFrame(results)
rdf.to_csv(os.path.join(OUT_DIR, "Z1_relaxation_variants.tsv"), sep='\t', index=False)

# ── Print Summary ──────────────────────────────────────────────────
for mode in ["paired", "to"]:
    print(f"\n{'='*80}")
    print(f"  {mode.upper()} MODE")
    print(f"{'='*80}")
    mdf = rdf[rdf["mode"] == mode]
    print(f"{'Variant':<25s} {'N':>7s} {'%':>6s} {'TP Rate':>8s} {'Delta':>8s} "
          f"{'Mean(S)':>8s} {'Min(S)':>8s} {'Pos':>4s} {'Sig':>4s}")
    print("-" * 85)
    for _, row in mdf.iterrows():
        print(f"{row['variant']:<25s} {row['n']:>7,d} {row['pct']:>5.1f}% "
              f"{row['tp_rate']:>7.4f} {row['delta']:>+7.4f} "
              f"{row['mean_sample_tp']:>7.4f} {row['min_sample_tp']:>7.4f} "
              f"{row['n_positive']:>3.0f}/{row['n_samples_with_data']:.0f} "
              f"{row['n_sig']:>3.0f}/{row['n_samples_with_data']:.0f}")

# ── Figure: Coverage vs TP Rate Tradeoff ──────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

for idx, mode in enumerate(["paired", "to"]):
    ax = axes[idx]
    mdf = rdf[rdf["mode"] == mode]

    colors = ['red', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    markers = ['X', 'o', 's', '^', 'D', 'v', 'p']

    for i, (_, row) in enumerate(mdf.iterrows()):
        ax.scatter(row["pct"], row["tp_rate"], s=150,
                  c=colors[i], marker=markers[i], edgecolors='black', zorder=5,
                  label=f'{row["variant"]} (N={row["n"]:,})')

    ax.axhline(y=mdf.iloc[0]["global_tp_rate"], color='gray', linestyle='--', alpha=0.5,
               label=f'Global ({mdf.iloc[0]["global_tp_rate"]:.3f})')

    ax.set_xlabel("Coverage (% of all regions)", fontsize=11)
    ax.set_ylabel("TP Rate", fontsize=11)
    ax.set_title(f"Z1 Variants: Coverage vs TP Rate — {mode.upper()}", fontweight='bold')
    ax.legend(fontsize=7, loc='lower left', framealpha=0.9)
    ax.grid(alpha=0.3)

    # Add annotation arrows for Pareto-optimal points
    if mode == "to":
        ax.annotate("← More selective\n   Higher TP rate", xy=(0.05, 0.95),
                   xycoords='axes fraction', fontsize=8, ha='left', color='gray')
        ax.annotate("More coverage →\n  Lower TP rate", xy=(0.7, 0.55),
                   xycoords='axes fraction', fontsize=8, ha='left', color='gray')

plt.suptitle("Z1 Definition Relaxation: Coverage vs TP Rate Tradeoff",
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Z1_relaxation_tradeoff.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"\n→ Figure: {FIG_DIR}/Z1_relaxation_tradeoff.png")

# ── Per-sample breakdown for top 3 TO variants ────────────────────
print("\n\n=== Per-Sample Detail for Top TO Variants ===")
dm = df[df["mode"] == "to"]

top_variants = ["Z1a_LOH+IntAF", "Z1b_LOH+IntAF+NG2", "Z1d_LOH+NG2"]
for vname in top_variants:
    vfunc = z1_variants[vname]
    print(f"\n--- {vname} ---")
    for sample in SAMPLE_ORDER:
        sm = dm[dm["sample"] == sample]
        if len(sm) == 0: continue
        z_mask = vfunc(sm)
        sn = z_mask.sum()
        if sn == 0:
            print(f"  {sample:16s}: N=0")
            continue
        stp = sm.loc[z_mask, "is_tp"].sum()
        stp_rate = stp / sn
        s_global = sm["is_tp"].mean()
        print(f"  {sample:16s}: N={sn:5d} TP_rate={stp_rate:.4f} Global={s_global:.4f} "
              f"delta={stp_rate - s_global:+.4f}")

print("\nDone.")
