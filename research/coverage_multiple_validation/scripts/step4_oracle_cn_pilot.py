#!/usr/bin/env python3
"""Step 4: Oracle CN pilot using corrected baselines and SEQC2 truth CN.

Context: Step 1 showed CovM Gain recall=14.6%. Step 2 discovered that
    expected_coverage is hardcoded at 75.0 across all 7 samples → CovM is
    cross-sample incomparable. This Step tests two corrections:

    Variant A — Re-center: CovM_v2 = NumReads / neutral_median_NumReads
                (per-sample, using SEQC2 for HCC1395 / proxy for others)
    Variant B — Oracle CN (HCC1395 only): replace CovM with SEQC2 truth CN
                directly, rerun Z3 amplicon blacklist (S2 style)

Outputs:
  - data/step4_covm_v2_per_sample.tsv     (re-centered Coverage_Multiple stats)
  - data/step4_oracle_pilot_results.tsv   (rerun S2 with various CN versions)
  - figures/step4_oracle_vs_covm_deltaf1.png
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
STEP2 = ("/big7_disk/liaoyoyo2001/InterSubMod/research/coverage_multiple_validation/"
         "data/step2_expected_coverage_per_sample.tsv")
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/coverage_multiple_validation"
DATA_DIR = os.path.join(OUT_DIR, "data")
FIG_DIR = os.path.join(OUT_DIR, "figures")

COLS = ["Chr", "Pos", "NumReads", "Coverage_Multiple",
        "sample", "mode", "truth_label",
        "core_loh_like", "to_loh_bed_hit",
        "caller_af", "HPFineNGroups"]

# HCC1954 cytogenetic amplicon bed (from Z3 pilot S2)
S2_WHOLE_CHR = [
    ("chr5", 0, 300_000_000),
    ("chr8", 0, 300_000_000),
    ("chr17", 0, 300_000_000),
]


def f1_from_mask(ds, mask, is_tp, fp_mask):
    """Compute F1 after removing rows where mask==True."""
    tp_kept = int((is_tp & (~mask)).sum())
    fp_kept = int((fp_mask & (~mask)).sum())
    tp_lost = int((is_tp & mask).sum())
    fp_lost = int((fp_mask & mask).sum())
    p = tp_kept / max(tp_kept + fp_kept, 1)
    r = tp_kept / max(tp_kept + tp_lost, 1)  # FN_orig = 0 for all samples
    f1 = 2 * p * r / max(p + r, 1e-12)
    return f1, p, r, tp_kept, fp_kept, tp_lost, fp_lost


def in_bed(chrom, pos, bed):
    for c, s, e in bed:
        if chrom == c and s <= pos < e:
            return True
    return False


def main():
    print("=" * 70)
    print("Step 4: Oracle CN pilot — corrected baselines vs SEQC2 truth")
    print("=" * 70)

    # ── Load Step 2 neutral baselines ─────────────────────────────────
    print("\nLoading Step 2 baselines...")
    ec = pd.read_csv(STEP2, sep='\t')
    # per (sample, mode) baseline map
    baseline_map = {(r["sample"], r["mode"]): r["neutral_median_NumReads"]
                    for _, r in ec.iterrows()}

    # ── Load master dataset (7 samples) for re-centering ──────────────
    print("\nLoading master dataset (all 7 samples, TO mode)...")
    chunks = []
    for chunk in pd.read_csv(MASTER, sep='\t', usecols=COLS,
                             chunksize=200_000, low_memory=False):
        chunks.append(chunk[chunk["mode"] == "to"].copy())
    df = pd.concat(chunks, ignore_index=True)
    for col in ["NumReads", "Coverage_Multiple", "Pos", "caller_af", "HPFineNGroups"]:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df["Chr"] = df["Chr"].astype(str)
    for col in ["core_loh_like", "to_loh_bed_hit"]:
        df[col] = df[col].map({"True": True, "False": False,
                               True: True, False: False}).fillna(False)
    print(f"  Loaded {len(df):,} TO rows (7 samples)")

    # ── Recompute Coverage_Multiple_v2 per sample ─────────────────────
    print("\nRecomputing Coverage_Multiple_v2 (per-sample neutral baseline)...")
    df["baseline_v2"] = df.apply(
        lambda r: baseline_map.get((r["sample"], r["mode"]), np.nan), axis=1)
    df["CovM_v2"] = df["NumReads"] / df["baseline_v2"]

    covm_stats = df.groupby("sample").agg(
        n=("CovM_v2", "size"),
        covm_old_median=("Coverage_Multiple", "median"),
        covm_v2_median=("CovM_v2", "median"),
        covm_old_p95=("Coverage_Multiple", lambda x: np.nanpercentile(x, 95)),
        covm_v2_p95=("CovM_v2", lambda x: np.nanpercentile(x, 95)),
    ).reset_index()
    covm_stats.to_csv(os.path.join(DATA_DIR, "step4_covm_v2_per_sample.tsv"),
                      sep='\t', index=False)
    print(covm_stats.to_string(index=False))

    # ── Compute zone flags ─────────────────────────────────────────────
    df["loh_flag"] = np.where(df["mode"] == "paired",
                              df["core_loh_like"], df["to_loh_bed_hit"])
    df["af_extreme"] = df["caller_af"].apply(
        lambda x: (not pd.isna(x)) and (x < 0.1 or x > 0.9))
    df["in_z3"] = df["loh_flag"] & df["af_extreme"] & (df["HPFineNGroups"] <= 1)
    df["is_tp"] = df["truth_label"] == "TP"
    df["is_fp"] = df["truth_label"] == "FP"

    # ── Variant A — S2 ∩ Z3 with new CovM_v2 for high-CN threshold ───
    # Strategy: reject Z3 rows where CovM_v2 > per-sample (non-Z3) 95th percentile
    print("\n=== Variant A — CovM_v2 > non-Z3 95%ile ∩ Z3 ===")

    rows = []
    for sample, g in df.groupby("sample"):
        is_tp = g["is_tp"].values
        is_fp = g["is_fp"].values
        in_z3 = g["in_z3"].values

        # Baseline F1
        tp_total = int(is_tp.sum())
        fp_total = int(is_fp.sum())
        f1_before = 2 * tp_total / (2 * tp_total + fp_total) if (tp_total + fp_total) else 0

        # Cutoff: non-Z3 CovM_v2 95th percentile
        nonz3_covm = g.loc[~g["in_z3"], "CovM_v2"].dropna().values
        if len(nonz3_covm) < 20:
            continue
        cutoff_v2 = float(np.nanpercentile(nonz3_covm, 95))

        # Cutoff (old): non-Z3 CovM_old 95th percentile
        nonz3_covm_old = g.loc[~g["in_z3"], "Coverage_Multiple"].dropna().values
        cutoff_old = float(np.nanpercentile(nonz3_covm_old, 95))

        for variant_name, cm_col, cutoff in [
            ("CovM_old_p95_capZ3", "Coverage_Multiple", cutoff_old),
            ("CovM_v2_p95_capZ3",  "CovM_v2",            cutoff_v2),
            ("S2_whole_chr_capZ3", None,                 None),
        ]:
            if variant_name == "S2_whole_chr_capZ3":
                mask = g.apply(lambda r: in_bed(r["Chr"], r["Pos"], S2_WHOLE_CHR)
                               if pd.notna(r["Pos"]) else False, axis=1).values
            else:
                mask = (g[cm_col].values > cutoff)
            mask = mask & in_z3

            f1, p, r_val, tk, fk, tl, fl = f1_from_mask(
                g, mask, is_tp, is_fp)
            rows.append({
                "sample": sample, "variant": variant_name,
                "cutoff": cutoff if cutoff else None,
                "f1_before": f1_before, "f1_after": f1,
                "delta_f1": f1 - f1_before,
                "tp_lost": tl, "fp_removed": fk - fp_total + int(mask.sum() & is_fp.sum()),  # recomputed below
            })
    res = pd.DataFrame(rows)

    # Recompute fp_removed properly
    rows2 = []
    for sample, g in df.groupby("sample"):
        is_tp = g["is_tp"].values
        is_fp = g["is_fp"].values
        in_z3 = g["in_z3"].values
        tp_total = int(is_tp.sum())
        fp_total = int(is_fp.sum())
        f1_before = 2 * tp_total / (2 * tp_total + fp_total) if (tp_total + fp_total) else 0

        nonz3_covm = g.loc[~g["in_z3"], "CovM_v2"].dropna().values
        nonz3_covm_old = g.loc[~g["in_z3"], "Coverage_Multiple"].dropna().values
        if len(nonz3_covm) < 20 or len(nonz3_covm_old) < 20:
            continue
        cutoff_v2 = float(np.nanpercentile(nonz3_covm, 95))
        cutoff_old = float(np.nanpercentile(nonz3_covm_old, 95))

        pos_series = g["Pos"].values
        chr_series = g["Chr"].values
        s2_mask = np.array([
            in_bed(c, p, S2_WHOLE_CHR) if pd.notna(p) else False
            for c, p in zip(chr_series, pos_series)
        ])

        for variant_name, mask in [
            ("CovM_old_p95_capZ3", (g["Coverage_Multiple"].values > cutoff_old) & in_z3),
            ("CovM_v2_p95_capZ3",  (g["CovM_v2"].values > cutoff_v2) & in_z3),
            ("S2_whole_chr_capZ3", s2_mask & in_z3),
        ]:
            f1, p, r_val, tk, fk, tl, fl = f1_from_mask(g, mask, is_tp, is_fp)
            rows2.append({
                "sample": sample, "variant": variant_name,
                "f1_before": f1_before, "f1_after": f1,
                "delta_f1": f1 - f1_before,
                "tp_lost": tl, "fp_removed": fl,
                "n_z3": int(in_z3.sum()),
            })
    res2 = pd.DataFrame(rows2)
    res2.to_csv(os.path.join(DATA_DIR, "step4_variantA_rerun.tsv"), sep='\t', index=False)
    print(res2.to_string(index=False))

    # ── Variant B — HCC1395 oracle CN (using SEQC2 truth) ──────────────
    print("\n=== Variant B — HCC1395 oracle CN rerun ===")
    aligned = pd.read_csv(ALIGNED, sep='\t')
    for col in ["NumReads", "Coverage_Multiple", "Pos", "truth_cn"]:
        aligned[col] = pd.to_numeric(aligned[col], errors='coerce')
    al_to = aligned[aligned["mode"] == "to"].copy()

    # Need to bring in zone flags — re-merge from df
    # Key on (sample, Chr, Pos)
    al_to["Chr"] = al_to["Chr"].astype(str)
    merge_cols = ["sample", "Chr", "Pos", "in_z3", "is_tp", "is_fp"]
    df_to_z = df[["sample", "Chr", "Pos", "in_z3", "is_tp", "is_fp"]].copy()
    al_merged = al_to.merge(df_to_z, on=["sample", "Chr", "Pos"], how="inner")
    print(f"  merged HCC1395 TO rows with zone flags: {len(al_merged):,}")

    rows3 = []
    for sample in ["HCC1395", "HCC1395_DORADO"]:
        g = al_merged[al_merged["sample"] == sample]
        if len(g) == 0:
            continue
        is_tp = g["is_tp"].values
        is_fp = g["is_fp"].values
        in_z3 = g["in_z3"].values
        tp_total = int(is_tp.sum())
        fp_total = int(is_fp.sum())
        f1_before = 2 * tp_total / (2 * tp_total + fp_total) if (tp_total + fp_total) else 0

        # Oracle: truth_cn ≥ 3 ∩ Z3
        for cn_cut, label in [(2.5, "truth_cn≥2.5"),
                              (3.5, "truth_cn≥3.5"),
                              (4.5, "truth_cn≥4.5")]:
            mask = (g["truth_cn"].values >= cn_cut) & in_z3
            f1, p, r_val, tk, fk, tl, fl = f1_from_mask(g, mask, is_tp, is_fp)
            rows3.append({
                "sample": sample, "variant": f"oracle_{label}_capZ3",
                "f1_before": f1_before, "f1_after": f1,
                "delta_f1": f1 - f1_before,
                "tp_lost": tl, "fp_removed": fl,
                "n_z3": int(in_z3.sum()),
            })
    res3 = pd.DataFrame(rows3)
    res3.to_csv(os.path.join(DATA_DIR, "step4_variantB_oracle.tsv"), sep='\t', index=False)
    print(res3.to_string(index=False))

    # ── Combined summary ──────────────────────────────────────────────
    all_results = pd.concat([res2, res3], ignore_index=True)
    all_results.to_csv(os.path.join(DATA_DIR, "step4_oracle_pilot_results.tsv"),
                       sep='\t', index=False)

    # ── Figure: ΔF1 per sample × variant ──────────────────────────────
    pivot = all_results.pivot_table(index="sample", columns="variant",
                                     values="delta_f1", aggfunc="first")
    fig, ax = plt.subplots(figsize=(13, 5.5))
    pivot.plot(kind='bar', ax=ax, width=0.75)
    ax.axhline(0, color='black', linewidth=0.6)
    ax.set_ylabel("ΔF1 (after − before)")
    ax.set_title("Step 4 — Z3 ∩ CN-filter with different CN definitions")
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1), fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "step4_oracle_vs_covm_deltaf1.png"),
                dpi=150, bbox_inches='tight')
    plt.close()

    # ── Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY — Variant A (CovM_v2) vs original (CovM_old)")
    print("=" * 70)
    cmp = res2.pivot_table(index="sample", columns="variant",
                            values="delta_f1", aggfunc="first").round(5)
    print(cmp.to_string())

    print("\n--- Variant B oracle (HCC1395 only) ---")
    print(res3.to_string(index=False))
    print("=" * 70)


if __name__ == "__main__":
    main()
