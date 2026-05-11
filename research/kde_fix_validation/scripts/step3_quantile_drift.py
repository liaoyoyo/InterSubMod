#!/usr/bin/env python3
"""Step 3: Cross-sample CovM quantile drift + Coverage_Category migration.

Compares stale-binary master (all TPs, mode in {paired, to}) vs new KDE-fix
master (14 combos TP, mode in {paired_pileup, paired_full}) per sample.

Outputs:
  - quantile_drift_per_sample.tsv  (p5/p25/p50/p75/p95, mean CovM per combo)
  - category_migration_per_sample.tsv  (6x6 migration estimate from distribution)
  - fig5_quantile_drift_per_sample.png  (14 panels: histogram overlay)
  - fig6_category_migration.png  (bar chart: stale vs fixed % per category)

Note: No per-region pairwise join — stale "paired" may not RegionID-align with
new "paired_pileup"/"paired_full" due to re-run differences. Distribution-level
comparison is sufficient for drift quantification.
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STALE = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
         "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
FIXED = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/kde_rerun_B_14combos/all_region_rows_kde_B_tp.tsv.gz"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/kde_fix_validation/outputs/step3_quantile_drift")
FIG_OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260420_kde_fix_acceptance")

SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437", "COLO829"]
CATEGORIES = ["CNV_Loss", "Low", "Normal", "Elevated", "CNV_Gain", "High_Copy"]

STALE_BASELINE = 75.0


def load_stale():
    cols = ["sample", "mode", "truth_label", "Coverage_Multiple", "Coverage_Category"]
    chunks = []
    for chk in pd.read_csv(STALE, sep="\t", usecols=cols, chunksize=200_000, low_memory=False):
        m = (chk["truth_label"] == "TP") & (chk["mode"] == "paired") & chk["sample"].isin(SAMPLES)
        if m.any():
            chunks.append(chk[m].copy())
    df = pd.concat(chunks, ignore_index=True)
    print(f"Stale TP paired: {len(df):,} rows")
    return df


def load_fixed():
    cols = ["sample", "mode", "Coverage_Multiple", "Diploid_Coverage_Used", "Coverage_Category"]
    df = pd.read_csv(FIXED, sep="\t", usecols=cols, low_memory=False)
    print(f"Fixed: {len(df):,} rows")
    return df


def quantiles(series, qs=(0.05, 0.25, 0.5, 0.75, 0.95)):
    return {f"p{int(q*100):02d}": float(np.nanpercentile(series, q * 100)) for q in qs}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    print("Loading masters...")
    stale = load_stale()
    fixed = load_fixed()

    # Per-sample quantiles
    rows = []
    for s in SAMPLES:
        st = stale[stale["sample"] == s]["Coverage_Multiple"].dropna()
        for mode in ["paired_pileup", "paired_full"]:
            fx = fixed[(fixed["sample"] == s) & (fixed["mode"] == mode)]
            fx_covm = fx["Coverage_Multiple"].dropna()
            if fx.empty:
                continue
            new_baseline = fx["Diploid_Coverage_Used"].iloc[0]
            row = {
                "sample": s,
                "fixed_mode": mode,
                "stale_baseline_x": STALE_BASELINE,
                "new_baseline_x": float(new_baseline),
                "ratio_stale_over_new": STALE_BASELINE / new_baseline,
                "n_stale": len(st),
                "n_fixed": len(fx_covm),
            }
            sq = quantiles(st)
            fq = quantiles(fx_covm)
            for k, v in sq.items():
                row[f"stale_{k}"] = v
            for k, v in fq.items():
                row[f"fixed_{k}"] = v
                row[f"delta_{k}"] = fq[k] - sq[k]
            row["stale_mean"] = float(st.mean())
            row["fixed_mean"] = float(fx_covm.mean())
            row["mean_shift"] = row["fixed_mean"] - row["stale_mean"]
            rows.append(row)

    quant = pd.DataFrame(rows)
    quant.to_csv(OUT / "quantile_drift_per_sample.tsv", sep="\t", index=False)
    print("\nQuantile drift (p50 stale vs fixed):")
    disp = quant[["sample", "fixed_mode", "new_baseline_x", "ratio_stale_over_new",
                  "stale_p50", "fixed_p50", "delta_p50"]].copy()
    print(disp.to_string(index=False))

    # Category distribution per combo
    cat_rows = []
    for s in SAMPLES:
        st = stale[stale["sample"] == s]
        st_dist = st["Coverage_Category"].value_counts(normalize=True).to_dict()
        for mode in ["paired_pileup", "paired_full"]:
            fx = fixed[(fixed["sample"] == s) & (fixed["mode"] == mode)]
            if fx.empty:
                continue
            fx_dist = fx["Coverage_Category"].value_counts(normalize=True).to_dict()
            for cat in CATEGORIES:
                cat_rows.append({
                    "sample": s, "fixed_mode": mode, "category": cat,
                    "stale_pct": 100 * st_dist.get(cat, 0.0),
                    "fixed_pct": 100 * fx_dist.get(cat, 0.0),
                    "delta_pp": 100 * (fx_dist.get(cat, 0.0) - st_dist.get(cat, 0.0)),
                })
    cat = pd.DataFrame(cat_rows)
    cat.to_csv(OUT / "category_migration_per_sample.tsv", sep="\t", index=False)

    # Figure 5: 14-panel histogram overlay (2 rows paired_pileup/full × 7 samples)
    fig, axes = plt.subplots(2, 7, figsize=(22, 9), sharex=True)
    for ci, s in enumerate(SAMPLES):
        st = stale[stale["sample"] == s]["Coverage_Multiple"].clip(0, 4).dropna()
        for ri, mode in enumerate(["paired_pileup", "paired_full"]):
            ax = axes[ri, ci]
            fx = fixed[(fixed["sample"] == s) & (fixed["mode"] == mode)]
            if fx.empty:
                ax.set_visible(False); continue
            fx_covm = fx["Coverage_Multiple"].clip(0, 4).dropna()
            nb = fx["Diploid_Coverage_Used"].iloc[0]
            bins = np.linspace(0, 4, 60)
            ax.hist(st.values, bins=bins, alpha=0.5, color="#cc4444",
                    label=f"Stale (75x)", density=True)
            ax.hist(fx_covm.values, bins=bins, alpha=0.5, color="#3366cc",
                    label=f"Fixed ({nb:.0f}x)", density=True)
            ax.axvline(1.0, color="black", linestyle="--", alpha=0.4, linewidth=0.8)
            ax.axvline(st.median(), color="#cc4444", linestyle=":", alpha=0.7, linewidth=1)
            ax.axvline(fx_covm.median(), color="#3366cc", linestyle=":", alpha=0.7, linewidth=1)
            ax.set_title(f"{s}\n{mode}  ratio={STALE_BASELINE/nb:.2f}x", fontsize=9)
            if ci == 0:
                ax.set_ylabel(f"Density ({mode})")
            if ri == 1:
                ax.set_xlabel("CovM")
            ax.grid(alpha=0.2)
            if ri == 0 and ci == 0:
                ax.legend(loc="upper right", fontsize=7)
    fig.suptitle("Cross-sample CovM distribution drift: stale 75x vs KDE-fixed baseline\n"
                 "(large stale/new ratio → strong right-shift; COLO829 ratio=2.59x most extreme)", y=1.00, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_OUT / "fig5_quantile_drift_per_sample.png", dpi=130, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {FIG_OUT / 'fig5_quantile_drift_per_sample.png'}")

    # Figure 6: Category migration bar chart
    fig, axes = plt.subplots(2, 7, figsize=(22, 7), sharey=True)
    x = np.arange(len(CATEGORIES))
    w = 0.38
    for ci, s in enumerate(SAMPLES):
        for ri, mode in enumerate(["paired_pileup", "paired_full"]):
            ax = axes[ri, ci]
            sub = cat[(cat["sample"] == s) & (cat["fixed_mode"] == mode)]
            if sub.empty:
                ax.set_visible(False); continue
            sub = sub.set_index("category").reindex(CATEGORIES)
            ax.bar(x - w/2, sub["stale_pct"], w, color="#cc4444", label="Stale")
            ax.bar(x + w/2, sub["fixed_pct"], w, color="#3366cc", label="Fixed")
            ax.set_xticks(x); ax.set_xticklabels(CATEGORIES, rotation=35, ha="right", fontsize=7)
            ax.set_title(f"{s} {mode}", fontsize=9)
            if ci == 0:
                ax.set_ylabel("% regions")
            if ri == 0 and ci == 0:
                ax.legend(fontsize=7)
            ax.grid(alpha=0.2, axis="y")
    fig.suptitle("Coverage_Category migration: stale 75x vs KDE-fixed per-sample baseline", y=1.00, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_OUT / "fig6_category_migration_per_sample.png", dpi=130, bbox_inches="tight")
    plt.close()
    print(f"Saved {FIG_OUT / 'fig6_category_migration_per_sample.png'}")

    # Summary
    print("\n" + "=" * 70)
    print("KEY DRIFT METRICS")
    print("=" * 70)
    summary = quant.groupby("sample").agg(
        new_baseline=("new_baseline_x", "first"),
        ratio=("ratio_stale_over_new", "first"),
        stale_p50=("stale_p50", "first"),
        fixed_p50_mean=("fixed_p50", "mean"),
        delta_p50_mean=("delta_p50", "mean"),
    ).sort_values("ratio", ascending=False)
    print(summary.to_string())

    # Most severe drift
    worst = summary.iloc[0]
    best = summary.iloc[-1]
    print(f"\nMost severe drift:  {summary.index[0]}  ratio={worst['ratio']:.2f}  Δp50={worst['delta_p50_mean']:.3f}")
    print(f"Least severe drift: {summary.index[-1]}  ratio={best['ratio']:.2f}   Δp50={best['delta_p50_mean']:.3f}")


if __name__ == "__main__":
    main()
