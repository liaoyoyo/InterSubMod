#!/usr/bin/env python3
"""R6 — Coverage_Multiple per-sample z-score normalization.

Question: Does HCC1954 Z3 chr5/8/17 Coverage_Multiple remain elevated after
normalizing against each sample's non-LOH baseline? If yes -> true CNV signal;
if no -> the "high CovM in Z3" is artifact of per-sample coverage scaling.

Strategy:
  1) Load master dataset (748K rows, 116 cols).
  2) Per (sample, mode): compute non-LOH CovM mean/std (baseline).
     - TO non-LOH: ~to_loh_bed_hit
     - paired non-LOH: ~tool_potential_loh & ~core_loh_like
  3) Z-score every row's Coverage_Multiple against its sample baseline.
  4) Define Z3 per Step 3 convention and tag extreme (z > 2).
  5) Validate HCC1954 chr5/8/17 FP enrichment under z-score definition.

Outputs:
  - data/r6_covm_per_sample_stats.tsv
  - data/r6_covm_zscore_regions.tsv.gz  (lightweight: sample,mode,Chr,Pos,CovM,z,in_z3,truth_label)
  - data/r6_hcc1954_z3_extreme_validation.tsv
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

MASTER = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz"
)
OUT_DIR = Path(__file__).resolve().parents[1] / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

USE_COLS = [
    "sample", "mode", "Chr", "Pos", "Coverage_Multiple",
    "Potential_LOH", "tool_potential_loh", "core_loh_like", "to_loh_bed_hit",
    "HPFineNGroups", "caller_af", "truth_label", "NumReads",
]


def load_master() -> pd.DataFrame:
    print(f"[R6] Loading {MASTER} ...")
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", usecols=USE_COLS, low_memory=False)
    # coerce bool-ish cols
    for c in ("tool_potential_loh", "core_loh_like", "to_loh_bed_hit"):
        if df[c].dtype != bool:
            df[c] = df[c].astype(str).str.lower().isin({"true", "1", "t"})
    df["Coverage_Multiple"] = pd.to_numeric(df["Coverage_Multiple"], errors="coerce")
    df["caller_af"] = pd.to_numeric(df["caller_af"], errors="coerce")
    df["HPFineNGroups"] = pd.to_numeric(df["HPFineNGroups"], errors="coerce")
    print(f"[R6]   rows={len(df):,}  samples={sorted(df['sample'].unique())}")
    return df


def non_loh_mask(df: pd.DataFrame) -> pd.Series:
    """Per-mode non-LOH mask."""
    to = (df["mode"] == "to") & (~df["to_loh_bed_hit"])
    pa = (df["mode"] == "paired") & (~df["tool_potential_loh"]) & (~df["core_loh_like"])
    return to | pa


def loh_mask(df: pd.DataFrame) -> pd.Series:
    to = (df["mode"] == "to") & (df["to_loh_bed_hit"])
    pa = (df["mode"] == "paired") & (df["tool_potential_loh"] | df["core_loh_like"])
    return to | pa


def compute_stats(df: pd.DataFrame) -> pd.DataFrame:
    nl = df[non_loh_mask(df)].copy()
    g = nl.groupby(["sample", "mode"])["Coverage_Multiple"]
    out = pd.DataFrame({
        "mean": g.mean(), "std": g.std(ddof=1), "median": g.median(),
        "q25": g.quantile(0.25), "q75": g.quantile(0.75), "count": g.count(),
    }).reset_index()
    return out


def apply_zscore(df: pd.DataFrame, stats: pd.DataFrame) -> pd.DataFrame:
    m = df.merge(stats[["sample", "mode", "mean", "std"]], on=["sample", "mode"], how="left")
    m["covm_zscore"] = (m["Coverage_Multiple"] - m["mean"]) / m["std"].clip(lower=1e-6)
    return m


def tag_z3(df: pd.DataFrame) -> pd.DataFrame:
    """Z3 per Step 3 convention: loh & extreme AF & HPFineNGroups <= 1."""
    loh = loh_mask(df)
    af_extreme = (df["caller_af"] < 0.1) | (df["caller_af"] > 0.9)
    low_ng = df["HPFineNGroups"].fillna(0) <= 1
    df["in_z3"] = loh & af_extreme & low_ng
    return df


def hcc1954_validation(df: pd.DataFrame) -> pd.DataFrame:
    sel = (df["sample"] == "HCC1954") & df["in_z3"]
    sub = df[sel].copy()
    sub["is_fp"] = sub["truth_label"].str.upper().eq("FP")
    sub["is_big3"] = sub["Chr"].isin(["chr5", "chr8", "chr17"])
    sub["z_extreme"] = sub["covm_zscore"] > 2
    rows = []
    for mode, g in sub.groupby("mode"):
        total = len(g)
        rows.append({
            "mode": mode, "metric": "total", "value": total,
            "frac": 1.0 if total else 0.0,
        })
        for chr_tag, m in [("big3", g["is_big3"]), ("fp", g["is_fp"]),
                           ("z_extreme", g["z_extreme"]),
                           ("big3_and_z_extreme", g["is_big3"] & g["z_extreme"]),
                           ("fp_big3", g["is_fp"] & g["is_big3"]),
                           ("fp_z_extreme", g["is_fp"] & g["z_extreme"]),
                           ("fp_big3_z_extreme", g["is_fp"] & g["is_big3"] & g["z_extreme"])]:
            n = int(m.sum())
            rows.append({
                "mode": mode, "metric": chr_tag, "value": n,
                "frac": (n / total) if total else 0.0,
            })
    return pd.DataFrame(rows)


def main():
    df = load_master()
    stats = compute_stats(df)
    stats_path = OUT_DIR / "r6_covm_per_sample_stats.tsv"
    stats.to_csv(stats_path, sep="\t", index=False)
    print(f"[R6] Wrote {stats_path}  rows={len(stats)}")

    df = apply_zscore(df, stats)
    df = tag_z3(df)

    light = df[["sample", "mode", "Chr", "Pos", "Coverage_Multiple", "covm_zscore",
                "in_z3", "truth_label", "HPFineNGroups", "caller_af",
                "to_loh_bed_hit", "tool_potential_loh", "core_loh_like"]].copy()
    light_path = OUT_DIR / "r6_covm_zscore_regions.tsv.gz"
    light.to_csv(light_path, sep="\t", index=False, compression="gzip")
    print(f"[R6] Wrote {light_path}  rows={len(light):,}")

    val = hcc1954_validation(df)
    val_path = OUT_DIR / "r6_hcc1954_z3_extreme_validation.tsv"
    val.to_csv(val_path, sep="\t", index=False)
    print(f"[R6] Wrote {val_path}")
    print("[R6] HCC1954 Z3 summary:")
    print(val.to_string(index=False))


if __name__ == "__main__":
    main()
