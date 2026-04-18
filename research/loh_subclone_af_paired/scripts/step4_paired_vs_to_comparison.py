#!/usr/bin/env python3
"""
Step 4: Paired vs TO Direct Comparison
========================================
Cross-mode analysis: compare effect sizes, AF concordance, NGroups agreement,
and LOH annotation consistency between paired and TO modes.

Input: all_region_rows.tsv.gz (master dataset, both modes) + LOH.bed + step2 results
Output:
  - p14: Cross-mode AF concordance scatter
  - p15: Cross-mode NGroups concordance
  - p16: Effect size forest plot (TO r vs Paired r)
  - p17: Cross-mode LOH annotation agreement
  - step4_effect_size_comparison.tsv
  - step4_cross_mode_summary.tsv
"""

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats
from utils import (
    MASTER, FIGDIR, DATADIR, SAMPLE_ORDER, SAMPLE_SHORT, SAMPLE_LABELS,
    AF_COLORS, classify_af, cn_tier, load_loh_beds, annotate_loh_bed
)
import warnings
warnings.filterwarnings("ignore")

FIGDIR.mkdir(parents=True, exist_ok=True)
DATADIR.mkdir(parents=True, exist_ok=True)

# TO mode data paths (for effect size comparison)
TO_DATADIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/loh_subclone_af/data"


def load_both_modes():
    """Load master dataset with both paired and TO modes."""
    print("Loading master dataset (both modes)...")
    cols = [
        "RegionID", "Chr", "Pos", "variant_key",
        "sample", "mode", "truth_label", "caller_af",
        "to_loh_bed_hit",
        "HPFineNGroups", "HPFineF", "HPFineP",
        "Coverage_Multiple", "NumReads", "HP_Ratio",
        "AlleleDelta", "Quality_Score",
        "HP1FamilyN", "HP2FamilyN",
    ]
    df = pd.read_csv(MASTER, sep="\t", usecols=cols, low_memory=False)
    print(f"  Total rows: {len(df):,}")

    df_paired = df[df["mode"] == "paired"].copy()
    df_to = df[df["mode"] == "to"].copy()
    print(f"  Paired: {len(df_paired):,}, TO: {len(df_to):,}")

    # Classifications
    for d in [df_paired, df_to]:
        d["af_class"] = d["caller_af"].apply(classify_af)
        d["cn_tier"] = d["Coverage_Multiple"].apply(cn_tier)

    return df_paired, df_to


def match_cross_mode_variants(df_paired, df_to):
    """Match same variants across paired and TO modes using variant_key + sample."""
    print("\nMatching cross-mode variants...")

    # Merge on variant_key + sample
    merged = df_paired.merge(
        df_to,
        on=["variant_key", "sample"],
        suffixes=("_paired", "_to"),
        how="inner"
    )
    print(f"  Matched variants: {len(merged):,}")

    for sample in SAMPLE_ORDER:
        n = (merged["sample"] == sample).sum()
        print(f"    {SAMPLE_SHORT[sample]}: {n:,}")

    return merged


def fig14_af_concordance(merged):
    """p14: Cross-mode AF concordance scatter."""
    print("\n=== Figure p14: AF Concordance ===")

    fig, axes = plt.subplots(2, 4, figsize=(24, 10))
    axes_flat = axes.flatten()

    concordance_records = []

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        sub = merged[merged["sample"] == sample]

        if len(sub) == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            ax.set_title(SAMPLE_SHORT[sample])
            continue

        # Subsample for plotting
        sub_plot = sub.sample(min(5000, len(sub)), random_state=42)

        # Color by truth label
        tp = sub_plot[sub_plot["truth_label_paired"] == "TP"]
        fp = sub_plot[sub_plot["truth_label_paired"] == "FP"]

        if len(tp) > 0:
            ax.scatter(tp["caller_af_to"], tp["caller_af_paired"],
                       alpha=0.15, s=3, c="#4CAF50", label=f"TP ({len(tp):,})")
        if len(fp) > 0:
            ax.scatter(fp["caller_af_to"], fp["caller_af_paired"],
                       alpha=0.15, s=3, c="#F44336", label=f"FP ({len(fp):,})")

        ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
        ax.set_xlabel("TO AF")
        ax.set_ylabel("Paired AF")
        ax.set_title(SAMPLE_SHORT[sample], fontsize=11, fontweight="bold")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=7, loc="lower right")

        # Correlation
        r, p = scipy_stats.pearsonr(sub["caller_af_to"], sub["caller_af_paired"])
        ax.text(0.05, 0.95, f"r={r:.3f}", transform=ax.transAxes,
                fontsize=10, va="top",
                bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.8))

        concordance_records.append({
            "sample": sample,
            "n_matched": len(sub),
            "pearson_r_af": r,
            "mean_af_diff": (sub["caller_af_paired"] - sub["caller_af_to"]).mean(),
            "median_af_diff": (sub["caller_af_paired"] - sub["caller_af_to"]).median(),
        })

    axes_flat[-1].axis("off")

    fig.suptitle("Step 4: AF Concordance — Paired vs TO Mode\n"
                 "Same variant_key matched across modes",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p14_cross_mode_af_concordance.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return pd.DataFrame(concordance_records)


def fig15_ngroups_concordance(merged):
    """p15: Cross-mode NGroups concordance."""
    print("\n=== Figure p15: NGroups Concordance ===")

    fig, axes = plt.subplots(2, 4, figsize=(24, 10))
    axes_flat = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        sub = merged[(merged["sample"] == sample) &
                     (merged["truth_label_paired"] == "TP")]

        if len(sub) == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            ax.set_title(SAMPLE_SHORT[sample])
            continue

        # 2D histogram (NGroups is integer-like)
        ng_to = sub["HPFineNGroups_to"].clip(0, 6)
        ng_paired = sub["HPFineNGroups_paired"].clip(0, 6)

        h, xedges, yedges = np.histogram2d(ng_to, ng_paired,
                                            bins=[np.arange(-0.5, 7.5, 1),
                                                  np.arange(-0.5, 7.5, 1)])
        # Normalize by row (TO NGroups)
        h_norm = h / (h.sum(axis=1, keepdims=True) + 1e-10) * 100

        im = ax.imshow(h_norm.T, origin="lower", aspect="auto",
                       extent=[-0.5, 6.5, -0.5, 6.5],
                       cmap="YlOrRd", vmin=0, vmax=100)

        # Add text counts
        for xi in range(int(h.shape[0])):
            for yi in range(int(h.shape[1])):
                if h[xi, yi] > 0:
                    ax.text(xi, yi, f"{int(h[xi, yi])}", ha="center", va="center",
                            fontsize=5, color="black" if h_norm[xi, yi] < 50 else "white")

        ax.plot([0, 6], [0, 6], "k--", alpha=0.3)
        ax.set_xlabel("TO NGroups")
        ax.set_ylabel("Paired NGroups")
        ax.set_title(SAMPLE_SHORT[sample], fontsize=11, fontweight="bold")

        # Agreement
        agree = (sub["HPFineNGroups_to"] == sub["HPFineNGroups_paired"]).mean()
        ax.text(0.05, 0.95, f"agree={agree:.1%}", transform=ax.transAxes,
                fontsize=9, va="top",
                bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.8))

    axes_flat[-1].axis("off")

    fig.suptitle("Step 4: NGroups Concordance — Paired vs TO Mode (TP only)\n"
                 "Cell values = variant count; color = row-normalized %",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p15_cross_mode_ngroups_concordance.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig16_effect_size_forest(paired_mw_path, to_mw_path):
    """p16: Forest plot comparing effect sizes across modes."""
    print("\n=== Figure p16: Effect Size Forest Plot ===")

    paired_mw = pd.read_csv(paired_mw_path, sep="\t")
    to_mw = pd.read_csv(to_mw_path, sep="\t")

    # Standardize sample names
    paired_mw["sample_short"] = paired_mw["sample"].map(SAMPLE_SHORT)
    to_mw["sample_short"] = to_mw["sample"].map(SAMPLE_SHORT)

    fig, ax = plt.subplots(figsize=(12, 8))

    samples = [SAMPLE_SHORT[s] for s in SAMPLE_ORDER]
    y_pos = np.arange(len(samples))

    # Plot TO r values
    to_r = []
    paired_r = []
    for s in SAMPLE_ORDER:
        s_short = SAMPLE_SHORT[s]
        to_row = to_mw[to_mw["sample_short"] == s_short]
        paired_row = paired_mw[paired_mw["sample_short"] == s_short]

        to_val = abs(to_row["rank_biserial_r"].values[0]) if len(to_row) > 0 else np.nan
        paired_val = abs(paired_row["rank_biserial_r"].values[0]) if len(paired_row) > 0 else np.nan

        to_r.append(to_val)
        paired_r.append(paired_val)

    ax.barh(y_pos + 0.15, to_r, 0.3, label="TO mode", color="#2196F3", alpha=0.8)
    ax.barh(y_pos - 0.15, paired_r, 0.3, label="Paired mode", color="#E91E63", alpha=0.8)

    # Add value labels
    for j, (t, p) in enumerate(zip(to_r, paired_r)):
        if not np.isnan(t):
            ax.text(t + 0.01, j + 0.15, f"{t:.3f}", va="center", fontsize=9)
        if not np.isnan(p):
            ax.text(p + 0.01, j - 0.15, f"{p:.3f}", va="center", fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples)
    ax.set_xlabel("|Rank-biserial r| (Intermediate vs Extreme NGroups, CN1)")
    ax.set_title("Effect Size Comparison: Paired vs TO Mode\n"
                 "Mann-Whitney |r| for NGroups difference (CN1 LOH TP)",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="lower right")
    ax.axvline(0.5, color="gray", ls="--", alpha=0.3)

    # Compute and display medians
    median_to = np.nanmedian(to_r)
    median_paired = np.nanmedian(paired_r)
    ax.text(0.98, 0.02,
            f"Median |r|:\n  TO = {median_to:.3f}\n  Paired = {median_paired:.3f}",
            transform=ax.transAxes, fontsize=11, ha="right", va="bottom",
            bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.9))

    plt.tight_layout()
    out = FIGDIR / "p16_cross_mode_effect_size_comparison.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    # Save comparison table
    comparison = pd.DataFrame({
        "sample": SAMPLE_ORDER,
        "sample_short": samples,
        "to_r": to_r,
        "paired_r": paired_r,
        "paired_minus_to": [p - t for p, t in zip(paired_r, to_r)],
        "direction_match": ["Yes" if (t > 0) == (p > 0) else "No"
                           for t, p in zip(to_r, paired_r)],
    })
    comparison.to_csv(DATADIR / "step4_effect_size_comparison.tsv", sep="\t", index=False)
    print(f"  Saved: {DATADIR / 'step4_effect_size_comparison.tsv'}")

    return comparison, median_to, median_paired


def fig17_loh_annotation_agreement(merged, loh_beds):
    """p17: LOH annotation agreement between to_loh_bed_hit and our LOH.bed annotation."""
    print("\n=== Figure p17: LOH Annotation Agreement ===")

    # Annotate paired variants with LOH.bed
    merged_copy = merged.copy()
    merged_copy["Chr"] = merged_copy.get("Chr", merged_copy.get("Chr_paired", None))
    merged_copy["Pos"] = merged_copy.get("Pos", merged_copy.get("Pos_paired", None))

    # Check which columns exist
    chr_col = "Chr" if "Chr" in merged_copy.columns else None
    pos_col = "Pos" if "Pos" in merged_copy.columns else None

    if chr_col is None:
        # Try to find from suffixed columns
        for c in merged_copy.columns:
            if "Chr" in c:
                chr_col = c
                break
        for c in merged_copy.columns:
            if "Pos" in c and "Potential" not in c:
                pos_col = c
                break

    if chr_col is None or pos_col is None:
        print("  Cannot find Chr/Pos columns, skipping p17")
        return

    # Rename for annotate_loh_bed
    df_for_annot = merged_copy.rename(columns={chr_col: "Chr", pos_col: "Pos"})
    df_for_annot["loh_bed_hit_ours"] = annotate_loh_bed(df_for_annot, loh_beds)

    # TO's own LOH annotation
    to_loh_col = "to_loh_bed_hit_to"
    if to_loh_col not in df_for_annot.columns:
        to_loh_col = "to_loh_bed_hit"
    if to_loh_col not in df_for_annot.columns:
        print("  Cannot find to_loh_bed_hit column, skipping p17")
        return

    df_for_annot["to_loh"] = df_for_annot[to_loh_col].astype(str).str.lower() == "true"

    fig, axes = plt.subplots(2, 4, figsize=(24, 10))
    axes_flat = axes.flatten()

    agreement_records = []

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        sub = df_for_annot[df_for_annot["sample"] == sample]

        if len(sub) == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            ax.set_title(SAMPLE_SHORT[sample])
            continue

        # 2x2 contingency
        both = ((sub["to_loh"]) & (sub["loh_bed_hit_ours"])).sum()
        to_only = ((sub["to_loh"]) & (~sub["loh_bed_hit_ours"])).sum()
        ours_only = ((~sub["to_loh"]) & (sub["loh_bed_hit_ours"])).sum()
        neither = ((~sub["to_loh"]) & (~sub["loh_bed_hit_ours"])).sum()

        matrix = np.array([[both, to_only], [ours_only, neither]])
        total = len(sub)

        # Plot confusion matrix
        im = ax.imshow(matrix, cmap="Blues", aspect="auto")

        for xi in range(2):
            for yi in range(2):
                pct = 100 * matrix[xi, yi] / total
                ax.text(yi, xi, f"{matrix[xi, yi]:,}\n({pct:.1f}%)",
                        ha="center", va="center", fontsize=9,
                        color="white" if matrix[xi, yi] > total * 0.3 else "black")

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["LOH.bed+", "LOH.bed-"])
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["TO LOH+", "TO LOH-"])
        ax.set_title(SAMPLE_SHORT[sample], fontsize=11, fontweight="bold")

        agree_pct = 100 * (both + neither) / total
        agreement_records.append({
            "sample": sample,
            "n_total": total,
            "both_loh": both,
            "to_only": to_only,
            "bed_only": ours_only,
            "neither": neither,
            "agreement_pct": agree_pct,
        })

    axes_flat[-1].axis("off")

    fig.suptitle("Step 4: LOH Annotation Agreement — to_loh_bed_hit vs LOH.bed Annotation\n"
                 "Same variant matched across modes; LOH.bed applied to both",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p17_cross_mode_loh_annotation_agreement.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return pd.DataFrame(agreement_records)


def compute_cross_mode_summary(concordance_df, effect_df, median_to, median_paired):
    """Compute and print cross-mode summary."""
    print("\n" + "=" * 80)
    print("STEP 4 SUMMARY: Paired vs TO Cross-Mode Comparison")
    print("=" * 80)

    print("\n  AF Concordance:")
    for _, row in concordance_df.iterrows():
        print(f"    {SAMPLE_SHORT[row['sample']]}: r={row['pearson_r_af']:.3f}, "
              f"mean diff={row['mean_af_diff']:.4f} (n={row['n_matched']:,})")

    print(f"\n  Effect Size Comparison (|r|, MW Intermediate vs Extreme NGroups, CN1):")
    n_match = 0
    for _, row in effect_df.iterrows():
        delta = row["paired_minus_to"]
        direction = "Paired stronger" if delta > 0 else "TO stronger"
        print(f"    {row['sample_short']}: TO={row['to_r']:.3f}, Paired={row['paired_r']:.3f}, "
              f"delta={delta:+.3f} ({direction})")
        if row["direction_match"] == "Yes":
            n_match += 1

    print(f"\n  Median |r|: TO = {median_to:.3f}, Paired = {median_paired:.3f}")
    print(f"  Direction consistency: {n_match}/{len(effect_df)} samples")

    # H2p verdict
    if median_paired >= median_to:
        print(f"\n  >>> H2p POSITIVE: Paired median |r| ({median_paired:.3f}) >= TO ({median_to:.3f})")
    else:
        print(f"\n  >>> H2p NEGATIVE: Paired median |r| ({median_paired:.3f}) < TO ({median_to:.3f})")

    # H4p verdict
    if n_match >= 5:
        print(f"  >>> H4p POSITIVE: {n_match}/{len(effect_df)} effect directions match")
    else:
        print(f"  >>> H4p NEGATIVE: only {n_match}/{len(effect_df)} directions match")

    # Save summary
    summary = pd.DataFrame({
        "hypothesis": ["H2p", "H4p"],
        "criterion": [
            f"Paired median |r| >= TO median |r|",
            f">=5/7 effect directions match",
        ],
        "result": [
            f"Paired={median_paired:.3f} vs TO={median_to:.3f}",
            f"{n_match}/7 match",
        ],
        "verdict": [
            "POSITIVE" if median_paired >= median_to else "NEGATIVE",
            "POSITIVE" if n_match >= 5 else "NEGATIVE",
        ],
    })
    summary.to_csv(DATADIR / "step4_cross_mode_summary.tsv", sep="\t", index=False)
    print(f"\n  Saved: {DATADIR / 'step4_cross_mode_summary.tsv'}")

    return summary


def main():
    print("=" * 70)
    print("Step 4: Paired vs TO Direct Comparison")
    print("=" * 70)

    # Load both modes
    df_paired, df_to = load_both_modes()

    # Load LOH.bed
    print("\nLoading LOH.bed files...")
    loh_beds = load_loh_beds()

    # Annotate paired mode with LOH.bed
    print(f"\nAnnotating paired variants with LOH.bed...")
    df_paired["loh_bed_hit"] = annotate_loh_bed(df_paired, loh_beds)

    # Match cross-mode variants
    merged = match_cross_mode_variants(df_paired, df_to)

    # Generate figures
    concordance_df = fig14_af_concordance(merged)
    fig15_ngroups_concordance(merged)

    # Effect size forest plot
    paired_mw_path = DATADIR / "step2_paired_ngroups_mw_test_cn1.tsv"
    to_mw_path = f"{TO_DATADIR}/step2_ngroups_mw_test_cn1.tsv"
    effect_df, median_to, median_paired = fig16_effect_size_forest(paired_mw_path, to_mw_path)

    # LOH annotation agreement
    agree_df = fig17_loh_annotation_agreement(merged, loh_beds)

    # Summary
    summary = compute_cross_mode_summary(concordance_df, effect_df, median_to, median_paired)

    print("\n" + "=" * 70)
    print("Step 4 complete!")
    print(f"  Figures: {FIGDIR}/p14-p17_*.png")
    print(f"  Data: {DATADIR}/step4_*.tsv")
    print("=" * 70)


if __name__ == "__main__":
    main()
