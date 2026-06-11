#!/usr/bin/env python3
"""
Phase 2 Cycle 2 Step C1 (H_C1_6) — Build per-BAM master TSVs.

The v1.0 master `step5_master_augmented.tsv` already joins methylation features
for all three BAMs (V3F / V5 / V6) onto the same 35,332-row chassis (caller_af,
LOH, coverage). For per-BAM cross-binary sanity we simply slice this master
into three sub-tables, renaming each BAM's columns to a canonical schema so
that `apply_cycle1_filter_per_bam.py` can reuse one feature list.

Canonical feature names (matched to cycle1_track_a_filter.json `feature_set`):

    bam_off_NG, bam_off_meth_HPMergedDelta, bam_off_meth_HPFineF,
    bam_off_meth_NME_imbalance, bam_off_meth_Epipoly_Delta,
    bam_off_meth_ClusterPermanovaF,
    caller_af, loh_inner_flag, Coverage_Multiple_imp, chr8_flag

Shared (BAM-invariant) columns retained as-is: region_id, chr, pos, label,
caller_af, Coverage_Multiple (+_imp), loh_side (+ loh_inner_flag), chr8_flag.

Output: cycle2/data/per_bam_master_{V3F,V5,V6}.tsv
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

BAM_TAGS = ["V3F", "V5", "V6"]

CYCLE1_FILTER_FEATURES = [
    "V6_off_NG",
    "caller_af",
    "loh_inner_flag",
    "V6_off_meth_HPMergedDelta",
    "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance",
    "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF",
    "chr8_flag",
    "Coverage_Multiple_imp",
]

# Canonical (BAM-agnostic) feature names exposed to downstream apply script.
CANONICAL_FEATURES = [
    "bam_off_NG",
    "caller_af",
    "loh_inner_flag",
    "bam_off_meth_HPMergedDelta",
    "bam_off_meth_HPFineF",
    "bam_off_meth_NME_imbalance",
    "bam_off_meth_Epipoly_Delta",
    "bam_off_meth_ClusterPermanovaF",
    "chr8_flag",
    "Coverage_Multiple_imp",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def prep_shared_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add shared derived columns (BAM-invariant)."""
    df = df.copy()
    df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    cov_median = df["Coverage_Multiple"].median()
    df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    df["y"] = (df["label"] == "TP").astype(int)
    return df


def slice_per_bam(df: pd.DataFrame, bam: str) -> pd.DataFrame:
    """Slice a per-BAM master from the joined master TSV."""
    rename_map = {
        f"{bam}_off_NG": "bam_off_NG",
        f"{bam}_off_meth_HPMergedDelta": "bam_off_meth_HPMergedDelta",
        f"{bam}_off_meth_HPFineF": "bam_off_meth_HPFineF",
        f"{bam}_off_meth_NME_imbalance": "bam_off_meth_NME_imbalance",
        f"{bam}_off_meth_Epipoly_Delta": "bam_off_meth_Epipoly_Delta",
        f"{bam}_off_meth_ClusterPermanovaF": "bam_off_meth_ClusterPermanovaF",
    }
    keep_cols = [
        "region_id", "chr", "pos", "label", "y",
        "caller_af", "Coverage_Multiple", "Coverage_Multiple_imp",
        "loh_side", "loh_inner_flag", "chr8_flag",
        *rename_map.keys(),
    ]
    missing = [c for c in keep_cols if c not in df.columns]
    if missing:
        raise KeyError(f"Missing columns for {bam}: {missing}")
    sub = df[keep_cols].copy().rename(columns=rename_map)
    return sub


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--master",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/"
                "v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/"
                "step5_master_augmented.tsv",
    )
    parser.add_argument(
        "--output-dir",
        default="/big7_disk/liaoyoyo2001/InterSubMod/research/"
                "methyl_augmented_filter_phase2/cycle2/data",
    )
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log(f"Loading master TSV: {args.master}")
    df = pd.read_csv(args.master, sep="\t", low_memory=False)
    log(f"  shape: {df.shape}")
    log(f"  labels: {df['label'].value_counts().to_dict()}")

    df = prep_shared_columns(df)

    for bam in BAM_TAGS:
        sub = slice_per_bam(df, bam)
        out_path = outdir / f"per_bam_master_{bam}.tsv"
        sub.to_csv(out_path, sep="\t", index=False)
        # Lightweight summary
        meth_na = sub["bam_off_meth_HPFineF"].isna().sum()
        ng_med = sub["bam_off_NG"].median()
        hpfine_med = sub["bam_off_meth_HPFineF"].median()
        log(
            f"  {bam}: rows={len(sub)} NG_median={ng_med:.2f} "
            f"HPFineF_median={hpfine_med:.4f} HPFineF_NA={meth_na}"
        )
        log(f"  -> {out_path}")

    log("Done.")


if __name__ == "__main__":
    main()
