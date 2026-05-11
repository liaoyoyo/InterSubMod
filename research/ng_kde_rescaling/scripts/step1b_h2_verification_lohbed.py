#!/usr/bin/env python3
"""
Step 1b — H2 verification: replicate 20260414 Test A using external LOH.bed
(instead of Potential_LOH fallback), TP-only, to isolate methodology effect
from true KDE-baseline effect on the observed Δ direction reversal.

20260414 canonical conditions:
  - LOH source: external tumor_phased_LOH.bed (TO pipeline)
  - Filter: truth_label == "TP" only
  - CN1: Coverage_Multiple < 0.75
  - Compare: Intermediate AF vs Extreme AF within CN1 LOH TP

This script applies the SAME filters to the NEW KDE master (step0 output)
and compares to the Potential_LOH fallback in step1.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import DATA_DIR, OBS_DIR, SAMPLE_ORDER

# ── 20260414 LOH.bed paths (from loh_subclone_af_paired utils) ──
LOH_BEDS = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "COLO829": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/duplicates/20260317_colo829_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed",
    "H1437": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "H2009": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1937": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1954": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
}


def annotate_loh_bed(df: pd.DataFrame, loh_beds: dict) -> pd.Series:
    """Vectorized LOH.bed hit using per-chr searchsorted."""
    hit = pd.Series(False, index=df.index)
    for sample, bed in loh_beds.items():
        mask = df["sample"] == sample
        if mask.sum() == 0:
            continue
        sub = df.loc[mask]
        for chrom in sub["Chr"].unique():
            chr_mask = mask & (df["Chr"] == chrom)
            if chr_mask.sum() == 0:
                continue
            bed_chr = bed[bed["chr"] == chrom].sort_values("start")
            if bed_chr.empty:
                continue
            starts = bed_chr["start"].values
            ends = bed_chr["end"].values
            pos = df.loc[chr_mask, "Pos"].values
            idx = np.searchsorted(starts, pos, side="right") - 1
            valid = (idx >= 0) & (idx < len(starts))
            in_seg = valid & (pos < ends[np.clip(idx, 0, len(ends) - 1)])
            hit.loc[chr_mask] = in_seg
    return hit


def load_beds() -> dict:
    beds = {}
    for s, path in LOH_BEDS.items():
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["chr", "start", "end", "label", "het_rate", "dot", "ts", "te", "rgb"])
        beds[s] = df
        print(f"  {s}: {len(df)} LOH segments (median size {df['end'].sub(df['start']).median():,.0f} bp)")
    return beds


def cohen_h_ng2(inter_vals: np.ndarray, ext_vals: np.ndarray) -> float:
    """Cohen's h for 'NG >= 2' proportions."""
    p_i = float((inter_vals >= 2).mean()) if len(inter_vals) else np.nan
    p_e = float((ext_vals >= 2).mean()) if len(ext_vals) else np.nan
    if np.isnan(p_i) or np.isnan(p_e):
        return np.nan
    phi_i = 2 * np.arcsin(np.sqrt(np.clip(p_i, 1e-9, 1 - 1e-9)))
    phi_e = 2 * np.arcsin(np.sqrt(np.clip(p_e, 1e-9, 1 - 1e-9)))
    return phi_i - phi_e


def main() -> None:
    master_path = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
    print(f"[step1b] loading {master_path.name}")
    master = pd.read_csv(master_path, sep="\t", low_memory=False)
    paired = master[master["mode"] == "paired_full"].copy()
    paired["Chr"] = paired["Chr"].astype(str)
    print(f"[step1b] paired_full n={len(paired):,}")

    print("[step1b] loading LOH.bed for 7 samples")
    beds = load_beds()

    print("[step1b] annotating loh_bed_hit (vectorized)...")
    paired["loh_bed_hit"] = annotate_loh_bed(paired, beds)
    print(f"  coverage: {paired['loh_bed_hit'].mean()*100:.1f}% regions hit LOH.bed")

    # Per-sample LOH.bed hit fraction
    print("[step1b] per-sample LOH.bed hit rate:")
    for s in SAMPLE_ORDER:
        sub = paired[paired["sample"] == s]
        if len(sub) == 0:
            continue
        hit_rate = sub["loh_bed_hit"].mean() * 100
        pot_rate = (sub["loh_side"] == "Inner").mean() * 100
        print(f"  {s}: LOH.bed={hit_rate:.1f}%  Potential_LOH={pot_rate:.1f}%  n={len(sub):,}")

    # ── Replication Test A (20260414 EXACT): LOH.bed + TP only + CN1 (CovM<0.75) ──
    print("\n[step1b] 20260414 EXACT replication (LOH.bed + TP only + CN1<0.75):")
    tp_label = paired["tp_label"] == 1
    cn1 = paired["CovM_used"] < 0.75
    loh_bed = paired["loh_bed_hit"] == True
    subset = paired[tp_label & cn1 & loh_bed].copy()
    print(f"  subset n={len(subset):,}")

    rows_bed = []
    print("\n  sample | ΔNG | Inter mean | Ext mean | n_i/n_e | MW p (Inter>Ext) | Cohen_h_NG2")
    print("  " + "-" * 95)
    for sample in SAMPLE_ORDER:
        grp = subset[subset["sample"] == sample]
        inter = grp[grp["AF_class"] == "Intermediate"]["HPFineNGroups"].dropna().values
        extr = grp[grp["AF_class"] == "Extreme"]["HPFineNGroups"].dropna().values
        if len(inter) < 5 or len(extr) < 5:
            print(f"  {sample:14s} SKIP  n_i={len(inter)} n_e={len(extr)}")
            rows_bed.append({
                "sample": sample, "method": "LOH_bed_TP_only",
                "n_inter": len(inter), "n_extreme": len(extr),
                "ng_mean_intermediate": float(inter.mean()) if len(inter) else np.nan,
                "ng_mean_extreme": float(extr.mean()) if len(extr) else np.nan,
                "delta_ng": np.nan, "cohen_h_ng2": np.nan, "mw_pvalue": np.nan, "status": "SKIP",
            })
            continue
        d = float(inter.mean() - extr.mean())
        h = cohen_h_ng2(inter, extr)
        try:
            stat, pv = stats.mannwhitneyu(inter, extr, alternative="greater")
        except ValueError:
            pv = np.nan
        sig = "POS" if (d > 0 and pv < 0.05) else ("NEG" if (d < 0) else "NS")
        print(f"  {sample:14s} {d:+.3f} | {inter.mean():.3f} | {extr.mean():.3f} | {len(inter):5d}/{len(extr):5d} | {pv:.2e} | {h:+.3f} [{sig}]")
        rows_bed.append({
            "sample": sample, "method": "LOH_bed_TP_only",
            "n_inter": len(inter), "n_extreme": len(extr),
            "ng_mean_intermediate": float(inter.mean()), "ng_mean_extreme": float(extr.mean()),
            "delta_ng": d, "cohen_h_ng2": float(h), "mw_pvalue": float(pv), "status": sig,
        })

    # ── Cross-check: Potential_LOH + TP only + CN1 (isolate TP filter effect) ──
    print("\n[step1b] Potential_LOH + TP only + CN1<0.75 (isolate TP filter effect):")
    pot = paired["loh_side"] == "Inner"
    subset_pot = paired[tp_label & cn1 & pot].copy()
    print(f"  subset n={len(subset_pot):,}")

    rows_pot = []
    for sample in SAMPLE_ORDER:
        grp = subset_pot[subset_pot["sample"] == sample]
        inter = grp[grp["AF_class"] == "Intermediate"]["HPFineNGroups"].dropna().values
        extr = grp[grp["AF_class"] == "Extreme"]["HPFineNGroups"].dropna().values
        if len(inter) < 5 or len(extr) < 5:
            rows_pot.append({
                "sample": sample, "method": "PotentialLOH_TP_only",
                "n_inter": len(inter), "n_extreme": len(extr),
                "delta_ng": np.nan, "cohen_h_ng2": np.nan, "mw_pvalue": np.nan, "status": "SKIP",
            })
            continue
        d = float(inter.mean() - extr.mean())
        h = cohen_h_ng2(inter, extr)
        try:
            stat, pv = stats.mannwhitneyu(inter, extr, alternative="greater")
        except ValueError:
            pv = np.nan
        sig = "POS" if (d > 0 and pv < 0.05) else ("NEG" if (d < 0) else "NS")
        rows_pot.append({
            "sample": sample, "method": "PotentialLOH_TP_only",
            "n_inter": len(inter), "n_extreme": len(extr),
            "ng_mean_intermediate": float(inter.mean()), "ng_mean_extreme": float(extr.mean()),
            "delta_ng": d, "cohen_h_ng2": float(h), "mw_pvalue": float(pv), "status": sig,
        })

    # ── Jaccard between LOH.bed and Potential_LOH ──
    print("\n[step1b] LOH.bed vs Potential_LOH Jaccard per sample (paired_full, all rows):")
    print("  sample | LOH.bed only | Potential only | Both | Jaccard")
    for sample in SAMPLE_ORDER:
        sub = paired[paired["sample"] == sample]
        if sub.empty:
            continue
        bed_only = ((sub["loh_bed_hit"] == True) & (sub["loh_side"] != "Inner")).sum()
        pot_only = ((sub["loh_bed_hit"] == False) & (sub["loh_side"] == "Inner")).sum()
        both = ((sub["loh_bed_hit"] == True) & (sub["loh_side"] == "Inner")).sum()
        union = bed_only + pot_only + both
        jac = both / union if union > 0 else np.nan
        print(f"  {sample:14s} | {bed_only:6d} | {pot_only:6d} | {both:6d} | {jac:.3f}")

    # Merge and write comparison
    df_cmp = pd.DataFrame(rows_bed + rows_pot)
    out = DATA_DIR / "ng_h2_verification_lohbed_vs_potential.tsv"
    df_cmp.to_csv(out, sep="\t", index=False)
    print(f"\n[step1b] wrote {out.name} ({len(df_cmp)} rows)")

    # Summary verdict
    bed_df = df_cmp[df_cmp["method"] == "LOH_bed_TP_only"]
    n_pos_bed = ((bed_df["status"] == "POS") | (bed_df["delta_ng"] > 0.1)).sum()
    n_neg_bed = ((bed_df["status"] == "NEG") | (bed_df["delta_ng"] < -0.1)).sum()
    print(f"\n[step1b] VERDICT (LOH.bed + TP-only, new KDE):")
    print(f"  ΔNG > +0.1 (POSITIVE, matches 20260414 direction): {n_pos_bed}/7")
    print(f"  ΔNG < -0.1 (NEGATIVE, reversal): {n_neg_bed}/7")
    if n_pos_bed >= 5:
        print("  → 20260414 POSITIVE direction PRESERVED under new KDE")
    elif n_neg_bed >= 5:
        print("  → 20260414 POSITIVE direction REVERSED under new KDE (CONFIRMED)")
    else:
        print("  → MIXED; no clear verdict")


if __name__ == "__main__":
    main()
