#!/usr/bin/env python3
"""
Shared utilities for LOH Subclone AF Paired Mode analysis.
Constants, data loading, classification functions, and LOH.bed annotation.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Paths ──
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUTDIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_subclone_af_paired")
FIGDIR = OUTDIR / "figures"
DATADIR = OUTDIR / "data"

# ── Sample order (fixed, matches TO analysis) ──
SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_SHORT = {
    "HCC1395": "HCC1395", "HCC1395_DORADO": "HCC1395D",
    "COLO829": "COLO829", "H1437": "H1437", "H2009": "H2009",
    "HCC1937": "HCC1937", "HCC1954": "HCC1954",
}
SAMPLE_LABELS = {
    "HCC1395": "HCC1395\n(breast)",
    "HCC1395_DORADO": "HCC1395\n(Dorado)",
    "COLO829": "COLO829\n(melanoma)",
    "H1437": "H1437\n(lung)",
    "H2009": "H2009\n(lung)",
    "HCC1937": "HCC1937\n(BRCA1)",
    "HCC1954": "HCC1954\n(breast)",
}

# ── LOH.bed paths (all 7 samples from TO pipeline) ──
LOH_BEDS = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "COLO829": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/duplicates/20260317_colo829_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed",
    "H1437": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "H2009": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1937": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1954": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
}

# ── AF classification thresholds ──
# Same as TO analysis for direct comparability
AF_EXTREME_LO, AF_EXTREME_HI = 0.1, 0.9
AF_NEARHALF_LO, AF_NEARHALF_HI = 0.4, 0.6

# ── CN tier thresholds (Coverage_Multiple) ──
CN1_MAX = 0.75
CN2_MAX = 1.25
CN3_MAX = 1.75

# ── NR bins for confound control ──
NR_BINS = [(10, 30), (30, 50), (50, 80), (80, 150), (150, 500)]
NR_LABELS = ["10-30", "30-50", "50-80", "80-150", "150+"]

# ── Colors ──
AF_COLORS = {"Extreme": "#1976D2", "Near-half": "#4CAF50", "Intermediate": "#FF9800"}
MODE_COLORS = {"paired": "#E91E63", "to": "#2196F3"}


def classify_af(af):
    """Classify AF into Extreme / Near-half / Intermediate."""
    if af < AF_EXTREME_LO or af > AF_EXTREME_HI:
        return "Extreme"
    elif AF_NEARHALF_LO <= af <= AF_NEARHALF_HI:
        return "Near-half"
    else:
        return "Intermediate"


def cn_tier(cm):
    """Classify Coverage_Multiple into CN tier."""
    if pd.isna(cm):
        return "Unknown"
    if cm < CN1_MAX:
        return "CN1"
    elif cm < CN2_MAX:
        return "CN2"
    elif cm < CN3_MAX:
        return "CN3"
    else:
        return "CN4+"


def load_loh_beds():
    """Load all LOH.bed files into dict of DataFrames."""
    beds = {}
    for sample, path in LOH_BEDS.items():
        try:
            df = pd.read_csv(path, sep="\t", header=None,
                             names=["chr", "start", "end", "label", "het_rate",
                                    "dot", "thick_start", "thick_end", "rgb"])
            df["sample"] = sample
            df["segment_id"] = [f"{sample}_{i}" for i in range(len(df))]
            df["segment_size"] = df["end"] - df["start"]
            beds[sample] = df
            print(f"  {sample}: {len(df)} segments, median size={df['segment_size'].median():,.0f}bp")
        except Exception as e:
            print(f"  {sample}: ERROR loading LOH.bed: {e}")
    return beds


def annotate_loh_bed(df_variants, loh_beds):
    """Annotate variants with LOH.bed hit status using vectorized searchsorted.

    Args:
        df_variants: DataFrame with 'sample', 'Chr', 'Pos' columns
        loh_beds: dict from load_loh_beds()

    Returns:
        Boolean Series indicating LOH.bed hit
    """
    loh_hit = pd.Series(False, index=df_variants.index)

    for sample, bed in loh_beds.items():
        mask = df_variants["sample"] == sample
        if mask.sum() == 0:
            continue

        sample_vars = df_variants.loc[mask]

        for chrom in sample_vars["Chr"].unique():
            chr_mask = mask & (df_variants["Chr"] == chrom)
            if chr_mask.sum() == 0:
                continue

            # LOH segments on this chromosome
            bed_chr = bed[bed["chr"] == chrom].sort_values("start")
            if len(bed_chr) == 0:
                continue

            starts = bed_chr["start"].values
            ends = bed_chr["end"].values
            positions = df_variants.loc[chr_mask, "Pos"].values

            # Find which segment each position falls into
            idx = np.searchsorted(starts, positions, side="right") - 1
            valid = (idx >= 0) & (idx < len(starts))
            in_segment = valid & (positions < ends[np.clip(idx, 0, len(ends) - 1)])

            loh_hit.loc[chr_mask] = in_segment

    return loh_hit


def load_paired_data(extra_cols=None):
    """Load master dataset, filter paired mode, add AF/CN classifications."""
    print("Loading master dataset...")
    base_cols = [
        "RegionID", "Chr", "Pos", "variant_key",
        "sample", "mode", "truth_label", "caller_af",
        "to_loh_bed_hit", "core_loh_like",
        "HPFineNGroups", "HPFineF", "HPFineP",
        "Coverage_Multiple", "NumReads", "HP_Ratio",
        "AlleleDelta", "Potential_LOH", "LOH_Subtype",
        "Quality_Score", "PairwiseMeanDist", "PairwiseMedianDist",
        "CramersV", "ClusterPermanovaF",
        "HP1FamilyN", "HP2FamilyN",
    ]
    if extra_cols:
        base_cols = list(set(base_cols + extra_cols))

    df = pd.read_csv(MASTER, sep="\t", usecols=base_cols, low_memory=False)
    print(f"  Total rows: {len(df):,}")

    df_paired = df[df["mode"] == "paired"].copy()
    print(f"  Paired mode rows: {len(df_paired):,}")

    # Classifications
    df_paired["af_class"] = df_paired["caller_af"].apply(classify_af)
    df_paired["cn_tier"] = df_paired["Coverage_Multiple"].apply(cn_tier)

    # Derived features
    df_paired["hp_imbalance"] = (df_paired["HP_Ratio"] - 0.5).abs()
    total_hp = df_paired["HP1FamilyN"].fillna(0) + df_paired["HP2FamilyN"].fillna(0)
    df_paired["effective_hp_reads"] = total_hp

    return df, df_paired
