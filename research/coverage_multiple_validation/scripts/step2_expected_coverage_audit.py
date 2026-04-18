#!/usr/bin/env python3
"""Step 2: Cross-sample expected_coverage baseline audit.

Context: Step 1 showed HCC1395 Gain recall = 14.6% → KDE-estimated baseline
    is pulled up by high-ploidy regions. Quantify the baseline bias across all
    7 samples using reverse-engineering from master dataset.

    bias[sample] = (estimated_expected_cov - true_neutral_median) / true_neutral_median
       >0  → baseline over-estimated  → CovM under-estimated for real gains
       <0  → baseline under-estimated → CovM over-estimated (false gains)

For HCC1395, "true neutral" is defined by SEQC2 Neutral BED regions.
For other 6 samples (no CNV truth), proxy neutral regions use:
    - chr1p (1-120Mb), chr2q (120-240Mb), chr3p (1-90Mb), chr9q (70-140Mb)
    (common somatic-neutral regions across cancer types; not universally true
     but a reasonable pan-cancer approximation)

Outputs:
  - data/step2_expected_coverage_per_sample.tsv
  - figures/step2_expected_coverage_bias.png
"""
import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

MASTER = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
          "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
ALIGNED = ("/big7_disk/liaoyoyo2001/InterSubMod/research/coverage_multiple_validation/"
           "data/step1_per_region_truth_aligned.tsv.gz")
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/coverage_multiple_validation"
DATA_DIR = os.path.join(OUT_DIR, "data")
FIG_DIR = os.path.join(OUT_DIR, "figures")

# Proxy neutral regions (hg38) for non-HCC1395 samples
# Select regions commonly NOT altered in multiple cancer types
PROXY_NEUTRAL = [
    ("chr1", 1_000_000, 120_000_000),   # chr1p
    ("chr2", 120_000_000, 240_000_000), # chr2q
    ("chr3", 1_000_000, 90_000_000),    # chr3p
    ("chr9", 70_000_000, 140_000_000),  # chr9q
]


def in_proxy_neutral(chrom, pos):
    for c, s, e in PROXY_NEUTRAL:
        if chrom == c and s <= pos < e:
            return True
    return False


def main():
    print("=" * 70)
    print("Step 2: Cross-sample expected_coverage baseline audit")
    print("=" * 70)

    # ── Load master (all 7 samples × modes) for reverse-engineering ───
    print("\nLoading master dataset...")
    COLS = ["Chr", "Pos", "NumReads", "Coverage_Multiple", "sample", "mode"]
    chunks = []
    for chunk in pd.read_csv(MASTER, sep='\t', usecols=COLS,
                             chunksize=200_000, low_memory=False):
        chunks.append(chunk)
    df = pd.concat(chunks, ignore_index=True)
    for col in ["NumReads", "Coverage_Multiple", "Pos"]:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df["Chr"] = df["Chr"].astype(str)
    print(f"  Loaded {len(df):,} rows")

    # ── Reverse-engineer expected_coverage per sample × mode ──────────
    # Coverage_Multiple = NumReads / expected_coverage
    # → expected_coverage = median(NumReads / Coverage_Multiple)
    print("\nReverse-engineering expected_coverage per sample × mode...")
    ec_rows = []
    for (sample, mode), g in df.groupby(["sample", "mode"]):
        valid = g.dropna(subset=["NumReads", "Coverage_Multiple"])
        valid = valid[valid["Coverage_Multiple"] > 0]
        if len(valid) == 0:
            continue
        ratios = valid["NumReads"] / valid["Coverage_Multiple"]
        expected_cov = float(np.median(ratios))
        global_median_nr = float(valid["NumReads"].median())
        global_mean_nr = float(valid["NumReads"].mean())
        ec_rows.append({
            "sample": sample, "mode": mode,
            "n_regions": len(valid),
            "expected_coverage_estimate": expected_cov,
            "global_median_NumReads": global_median_nr,
            "global_mean_NumReads": global_mean_nr,
        })
    ec_df = pd.DataFrame(ec_rows)

    # ── Compute "true neutral" reference per sample × mode ────────────
    # HCC1395: use SEQC2 Neutral BED (from aligned data)
    # Others: use proxy neutral regions
    print("\nComputing true neutral baseline per sample × mode...")
    ref_rows = []

    # HCC1395 from aligned TSV
    aligned = pd.read_csv(ALIGNED, sep='\t')
    for col in ["NumReads", "Coverage_Multiple", "Pos"]:
        aligned[col] = pd.to_numeric(aligned[col], errors='coerce')
    for (sample, mode), g in aligned.groupby(["sample", "mode"]):
        neutral = g[(g["truth_label_seqc2"] == "Neutral")
                    & g["NumReads"].notna()]
        if len(neutral) > 10:
            ref_rows.append({
                "sample": sample, "mode": mode,
                "neutral_source": "SEQC2_Neutral_BED",
                "n_neutral": len(neutral),
                "neutral_median_NumReads": float(neutral["NumReads"].median()),
                "neutral_mean_NumReads": float(neutral["NumReads"].mean()),
                "neutral_median_CovM": float(neutral["Coverage_Multiple"].median()),
            })

    # Other samples: proxy neutral
    other_samples = [s for s in df["sample"].unique()
                     if s not in {"HCC1395", "HCC1395_DORADO"}]
    for sample in other_samples:
        sub_s = df[df["sample"] == sample]
        mask = sub_s.apply(lambda r: in_proxy_neutral(r["Chr"], r["Pos"])
                           if pd.notna(r["Pos"]) else False, axis=1)
        neutral = sub_s[mask & sub_s["NumReads"].notna()]
        for mode, g in neutral.groupby("mode"):
            if len(g) < 10:
                continue
            ref_rows.append({
                "sample": sample, "mode": mode,
                "neutral_source": "proxy_chr1p_2q_3p_9q",
                "n_neutral": len(g),
                "neutral_median_NumReads": float(g["NumReads"].median()),
                "neutral_mean_NumReads": float(g["NumReads"].mean()),
                "neutral_median_CovM": float(g["Coverage_Multiple"].median()),
            })

    ref_df = pd.DataFrame(ref_rows)

    # ── Merge and compute bias ────────────────────────────────────────
    out = ec_df.merge(ref_df, on=["sample", "mode"], how="left")
    out["bias_vs_neutral"] = (out["expected_coverage_estimate"]
                              - out["neutral_median_NumReads"]) / out["neutral_median_NumReads"]
    out["bias_pct"] = out["bias_vs_neutral"] * 100
    # If the baseline is over-estimated vs true neutral → positive bias
    # CovM = NumReads / baseline, so over-estimated baseline → under-estimated CovM
    # This explains HCC1395 Gain recall=14.6%

    out = out.sort_values(["sample", "mode"]).reset_index(drop=True)
    out.to_csv(os.path.join(DATA_DIR, "step2_expected_coverage_per_sample.tsv"),
               sep='\t', index=False)

    print("\n=== Expected Coverage Bias Table ===")
    display_cols = ["sample", "mode", "n_regions",
                    "expected_coverage_estimate", "neutral_median_NumReads",
                    "neutral_median_CovM", "bias_pct", "neutral_source"]
    print(out[display_cols].to_string(index=False))

    # ── Figure: baseline bias per sample ──────────────────────────────
    to_only = out[out["mode"] == "to"].copy()
    paired_only = out[out["mode"] == "paired"].copy()

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    # Panel 1: expected_coverage vs neutral baseline
    for ax, subset, title in [
        (axes[0], to_only, "TO mode"),
        (axes[1], paired_only, "Paired mode"),
    ]:
        x = np.arange(len(subset))
        w = 0.4
        ax.bar(x - w/2, subset["expected_coverage_estimate"], w,
               label='estimated expected_coverage', color='#e67e22', alpha=0.8)
        ax.bar(x + w/2, subset["neutral_median_NumReads"], w,
               label='neutral median NumReads (truth)', color='#27ae60', alpha=0.8)
        for i, row in subset.reset_index().iterrows():
            bias = row["bias_pct"]
            if pd.notna(bias):
                color = 'red' if abs(bias) > 15 else 'black'
                ax.text(i, max(row["expected_coverage_estimate"],
                                row["neutral_median_NumReads"]) * 1.05,
                        f"{bias:+.1f}%", ha='center', fontsize=9, color=color,
                        fontweight='bold' if abs(bias) > 15 else 'normal')
        ax.set_xticks(x)
        ax.set_xticklabels(subset["sample"], rotation=30, ha='right')
        ax.set_ylabel("coverage (reads)")
        ax.set_title(f"{title} — KDE-estimated baseline vs neutral-region truth  "
                     f"(red = bias >15%)")
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "step2_expected_coverage_bias.png"),
                dpi=150, bbox_inches='tight')
    plt.close()

    # ── Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY — H-CN1 判定")
    print("=" * 70)
    print("Bias >15% indicates KDE baseline mis-calibration:")
    for _, row in out.iterrows():
        bias = row["bias_pct"]
        if pd.isna(bias):
            continue
        flag = "⚠ POSITIVE" if abs(bias) > 15 else "   ok    "
        print(f"  {flag}  {row['sample']:<18} {row['mode']:<7} bias = {bias:+6.1f}%  "
              f"(est={row['expected_coverage_estimate']:.1f}, "
              f"neutral={row['neutral_median_NumReads']:.1f})")

    hcc1954_bias = out[(out["sample"] == "HCC1954") & (out["mode"] == "to")]["bias_pct"].values
    if len(hcc1954_bias):
        print(f"\nHCC1954 (TO) bias = {hcc1954_bias[0]:+.1f}%")
    hcc1395_bias = out[(out["sample"] == "HCC1395") & (out["mode"] == "to")]["bias_pct"].values
    if len(hcc1395_bias):
        print(f"HCC1395 (TO) bias = {hcc1395_bias[0]:+.1f}%  (SEQC2-verified neutral)")

    print("=" * 70)


if __name__ == "__main__":
    main()
