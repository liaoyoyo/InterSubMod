#!/usr/bin/env python3
"""Step 2 (KDE-fix rerun): Z3 amplicon blacklist ΔF1 with corrected baselines.

Context: Z3 flag is TO-mode specific (loh_flag + af_extreme + HPFineNGroups<=1).
The 14-combo rerun is paired-only (no TO); but Z3 analysis needs TO data.
Solution: use stale master's TO rows but RECOMPUTE Coverage_Multiple using new
per-sample KDE baselines (NumReads / new_baseline). This isolates S3 strategy's
sensitivity to the baseline fix, while S1/S2 (coordinate-based) are unchanged
by construction (sanity check).

Outputs:
  - step2_summary.tsv  (14 rows: 7 samples × S1/S2/S3 old vs new)
  - step2_compare.md   (human-readable before/after)
"""
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

STALE = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
         "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
AGG = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/kde_rerun_B_14combos/aggregate_summary.tsv"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/kde_fix_validation/outputs/step2_z3_amplicon")

STALE_BASELINE = 75.0
SAMPLE_ORDER = ["H2009", "H1437", "HCC1395", "HCC1395_DORADO",
                "HCC1937", "HCC1954", "COLO829"]

# Use paired_pileup baseline as canonical per-sample KDE baseline
# (paired_full and paired_pileup agree within ±2x per aggregate.tsv)
NEW_BASELINE = {
    "HCC1395": 55.0, "HCC1395_DORADO": 53.0,
    "HCC1937": 91.0, "HCC1954": 61.0,
    "H2009": 79.0, "H1437": 69.0, "COLO829": 29.0,
}

S1_LITERATURE_ARMS = [
    ("chr5",   1_000_000, 180_000_000),
    ("chr8",           0,  45_000_000),
    ("chr8",  120_000_000, 145_000_000),
    ("chr17",          0,  25_000_000),
    ("chr17", 35_000_000,  42_000_000),
]
S2_WHOLE_CHR = [
    ("chr5",  0, 300_000_000),
    ("chr8",  0, 300_000_000),
    ("chr17", 0, 300_000_000),
]


def in_blacklist(df, coords):
    mask = np.zeros(len(df), dtype=bool)
    for chrom, start, end in coords:
        hit = ((df["Chr"] == chrom) & (df["Pos"] >= start) & (df["Pos"] <= end))
        mask |= hit.values
    return mask


def compute_f1(tp, fp, fn):
    if tp + fp == 0 or tp + fn == 0:
        return 0.0, 0.0, 0.0
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f1


def run_strategy(df_sample, coords=None, covm_thresh=None, covm_col="CovM"):
    tp_total = int(df_sample["is_tp"].sum())
    fp_total = int((df_sample["truth_label"] == "FP").sum())
    fn_total = int(df_sample["is_fn"].sum())

    if coords is not None:
        bl_mask = in_blacklist(df_sample, coords)
    elif covm_thresh is not None:
        bl_mask = (df_sample[covm_col] > covm_thresh).fillna(False).values
    else:
        bl_mask = np.zeros(len(df_sample), dtype=bool)

    mask = bl_mask & df_sample["in_z3"].values

    p0, r0, f0 = compute_f1(tp_total, fp_total, fn_total)
    tp_kept = int((df_sample["is_tp"] & (~mask)).sum())
    fp_kept = int(((df_sample["truth_label"] == "FP") & (~mask)).sum())
    tp_lost = int((df_sample["is_tp"] & mask).sum())
    fp_removed = int(((df_sample["truth_label"] == "FP") & mask).sum())
    fn_new = fn_total + tp_lost
    p1, r1, f1n = compute_f1(tp_kept, fp_kept, fn_new)
    return {
        "f1_before": f0, "f1_after": f1n, "f1_delta": f1n - f0,
        "fp_removed": fp_removed, "tp_lost": tp_lost,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    cols = ["sample", "mode", "truth_label", "Chr", "Pos",
            "core_loh_like", "to_loh_bed_hit",
            "caller_af", "HPFineNGroups", "Coverage_Multiple", "NumReads"]
    print(f"Loading {STALE}...")
    df = pd.read_csv(STALE, sep="\t", usecols=cols, low_memory=False)
    for c in ["caller_af", "HPFineNGroups", "Coverage_Multiple", "Pos", "NumReads"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ["core_loh_like", "to_loh_bed_hit"]:
        df[c] = df[c].map({"True": True, "False": False,
                           True: True, False: False}).fillna(False)

    df["is_tp"] = df["truth_label"] == "TP"
    df["is_fn"] = df["truth_label"] == "FN"
    df["loh_flag"] = np.where(df["mode"] == "paired",
                               df["core_loh_like"], df["to_loh_bed_hit"])

    df_to = df[df["mode"] == "to"].copy()
    df_to["af_extreme"] = df_to["caller_af"].apply(
        lambda x: (not pd.isna(x)) and (x < 0.1 or x > 0.9))
    df_to["in_z3"] = (df_to["loh_flag"] & df_to["af_extreme"]
                      & (df_to["HPFineNGroups"] <= 1))

    # Stale CovM column (as in master, baseline=75x)
    df_to["CovM_stale"] = df_to["Coverage_Multiple"]
    # Fixed CovM: NumReads / per-sample KDE baseline
    df_to["CovM_fixed"] = df_to.apply(
        lambda r: r["NumReads"] / NEW_BASELINE.get(r["sample"], STALE_BASELINE)
        if pd.notna(r["NumReads"]) else np.nan, axis=1)

    # Sample-specific 95th percentile cutoffs (non-Z3 baseline), for both
    cut_stale, cut_fixed = {}, {}
    for s in SAMPLE_ORDER:
        base = df_to[(df_to["sample"] == s) & (~df_to["in_z3"])]
        cut_stale[s] = base["CovM_stale"].dropna().quantile(0.95)
        cut_fixed[s] = base["CovM_fixed"].dropna().quantile(0.95)

    print("\n95th percentile cutoffs:")
    print(f"{'sample':<16s}  {'stale':>8s}  {'fixed':>8s}  ratio")
    for s in SAMPLE_ORDER:
        print(f"{s:<16s}  {cut_stale[s]:>8.3f}  {cut_fixed[s]:>8.3f}"
              f"  {cut_stale[s]/cut_fixed[s] if cut_fixed[s] else np.nan:>5.3f}")

    # Run S1, S2, S3(stale), S3(fixed) per sample
    rows = []
    for s in SAMPLE_ORDER:
        ds = df_to[df_to["sample"] == s]
        r1 = run_strategy(ds, coords=S1_LITERATURE_ARMS)
        r2 = run_strategy(ds, coords=S2_WHOLE_CHR)
        r3s = run_strategy(ds, covm_thresh=cut_stale[s], covm_col="CovM_stale")
        r3f = run_strategy(ds, covm_thresh=cut_fixed[s], covm_col="CovM_fixed")
        rows.append({
            "sample": s,
            "new_baseline_x": NEW_BASELINE.get(s, STALE_BASELINE),
            "f1_before": r1["f1_before"],
            "s1_delta_f1": r1["f1_delta"], "s1_fp_rm": r1["fp_removed"], "s1_tp_lost": r1["tp_lost"],
            "s2_delta_f1": r2["f1_delta"], "s2_fp_rm": r2["fp_removed"], "s2_tp_lost": r2["tp_lost"],
            "s3_stale_cutoff": cut_stale[s],
            "s3_stale_delta_f1": r3s["f1_delta"], "s3_stale_fp_rm": r3s["fp_removed"], "s3_stale_tp_lost": r3s["tp_lost"],
            "s3_fixed_cutoff": cut_fixed[s],
            "s3_fixed_delta_f1": r3f["f1_delta"], "s3_fixed_fp_rm": r3f["fp_removed"], "s3_fixed_tp_lost": r3f["tp_lost"],
            "s3_delta_change": r3f["f1_delta"] - r3s["f1_delta"],
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "step2_summary.tsv", sep="\t", index=False)

    print("\n" + "=" * 90)
    print("Strategy comparison (ΔF1 per sample, TO mode, stale vs fixed CovM)")
    print("=" * 90)
    cols_show = ["sample", "f1_before", "s1_delta_f1", "s2_delta_f1",
                 "s3_stale_delta_f1", "s3_fixed_delta_f1", "s3_delta_change"]
    print(summary[cols_show].to_string(index=False))

    # Sanity check: S1/S2 should be unchanged by CovM (both re-derived from stale data only)
    print("\nSanity: S1/S2 are coordinate-based only, so fixing CovM has 0 effect on them.")
    print(f"  S1 Σ|Δ| = {summary['s1_delta_f1'].abs().sum():.4f}")
    print(f"  S2 Σ|Δ| = {summary['s2_delta_f1'].abs().sum():.4f}")
    print(f"  S3 stale→fixed change Σ|Δ| = {summary['s3_delta_change'].abs().sum():.4f}")

    # Markdown table
    md = OUT / "step2_compare.md"
    with open(md, "w") as f:
        f.write("# Step 2: Z3 amplicon blacklist — stale vs KDE-fixed ΔF1\n\n")
        f.write("## Method\n\n")
        f.write("- Stale CovM = NumReads / 75.0 (from master `Coverage_Multiple`)\n")
        f.write("- Fixed CovM = NumReads / per-sample KDE baseline\n")
        f.write("- Z3 flag (TO mode, loh_flag + af_extreme + HPFineNGroups<=1) unchanged\n")
        f.write("- S1 (literature arms) + S2 (whole-chr) coordinate-based → CovM-independent (sanity ∼0)\n")
        f.write("- S3 (CovM > 95%ile non-Z3) → per-sample cutoff changes with baseline\n\n")
        f.write("## Per-sample ΔF1\n\n")
        f.write("| Sample | Baseline new | F1 base | S1 Δ | S2 Δ | S3 stale Δ | S3 fixed Δ | S3 change |\n")
        f.write("|--------|-------------:|--------:|-----:|-----:|-----------:|-----------:|----------:|\n")
        for _, r in summary.iterrows():
            f.write(f"| {r['sample']} | {r['new_baseline_x']:.0f}x | "
                    f"{r['f1_before']:.4f} | "
                    f"{r['s1_delta_f1']:+.4f} | {r['s2_delta_f1']:+.4f} | "
                    f"{r['s3_stale_delta_f1']:+.4f} | {r['s3_fixed_delta_f1']:+.4f} | "
                    f"{r['s3_delta_change']:+.4f} |\n")
        f.write("\n## 95th percentile cutoffs\n\n")
        f.write("| Sample | Stale cutoff | Fixed cutoff | ratio stale/fixed |\n")
        f.write("|--------|-------------:|-------------:|------------------:|\n")
        for s in SAMPLE_ORDER:
            ratio = cut_stale[s] / cut_fixed[s] if cut_fixed[s] else float("nan")
            f.write(f"| {s} | {cut_stale[s]:.3f} | {cut_fixed[s]:.3f} | {ratio:.3f} |\n")
    print(f"\nWrote {md}")
    print(f"Wrote {OUT}/step2_summary.tsv")


if __name__ == "__main__":
    main()
