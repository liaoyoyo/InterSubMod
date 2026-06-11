"""Shared helpers for Step 1 + Step 2 of step5_methyl_filter_pilot.

Loads the augmented master TSV, slices powered cells per V6_off coords (the
canonical partition from Step 2), fits Models A/B per (BAM × flag × cell),
runs LRT df=5, BH-FDR across all 138 (BAM × flag × cell) combinations.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn")
STEP5_DIR = PROJECT_ROOT / "step5_methyl_filter_pilot"
STEP2_GRID = PROJECT_ROOT / "step2_three_axis_grid" / "step2_grid_3d.tsv"
STEP3_CHR8 = PROJECT_ROOT / "step3_fp_zone_zoom" / "chr8_hotspot" / "zone_grid_real_cov.tsv"
STEP3_AUTO = PROJECT_ROOT / "step3_fp_zone_zoom" / "auto_detected_top5pct" / "zone_grid.tsv"
STEP3_CROSS_SUMMARY = PROJECT_ROOT / "step3_fp_zone_zoom" / "step3_cross_zone_summary.tsv"
MASTER_AUG = STEP5_DIR / "step5_master_augmented.tsv"

LOG_FILE = STEP5_DIR / "intermediate" / "step1_log.txt"

POWERED_THRESHOLD = 50
HP_BUCKETS = ["same_HP1", "same_HP2", "cross_het", "cross_het_inv", "other"]
COV_LEVELS = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]
BAMS = ["V3F", "V5", "V6"]
FLAGS = ["off", "on"]

METHYL_COVARIATES = [
    "HPMergedDelta",
    "HPFineF",
    "NME_imbalance",
    "Epipoly_Delta",
    "ClusterPermanovaF",
]


def log_msg(msg: str, also_print: bool = True) -> None:
    import datetime as _dt

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    ts = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    if also_print:
        print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


def matplotlib_setup():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = [
        "Droid Sans Fallback",
        "DejaVu Sans",
        "Liberation Sans",
        "sans-serif",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    return plt


# ---- HP bucket / cov_bin derivation (mirrored from step2/_common.py) ----

def compute_hp_bucket(df: pd.DataFrame, prefix: str) -> pd.Series:
    hp1 = df[f"{prefix}_1"].fillna(0).astype(int)
    hp1s = df[f"{prefix}_1-1"].fillna(0).astype(int)
    hp2 = df[f"{prefix}_2"].fillna(0).astype(int)
    hp2s = df[f"{prefix}_2-1"].fillna(0).astype(int)
    bucket = pd.Series("other", index=df.index, dtype=object)
    bucket.loc[(hp1 > 0) & (hp1s > 0) & (hp2 == 0) & (hp2s == 0)] = "same_HP1"
    bucket.loc[(hp2 > 0) & (hp2s > 0) & (hp1 == 0) & (hp1s == 0)] = "same_HP2"
    bucket.loc[(hp1 > 0) & (hp2s > 0) & (hp1s == 0) & (hp2 == 0)] = "cross_het"
    bucket.loc[(hp2 > 0) & (hp1s > 0) & (hp2s == 0) & (hp1 == 0)] = "cross_het_inv"
    return bucket


def compute_cov_bin(s: pd.Series) -> pd.Series:
    bins = [-np.inf, 0.7, 1.3, 1.8, 2.5, np.inf]
    labels = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]
    return pd.cut(s, bins=bins, labels=labels, right=False).astype(object).fillna("unknown")


def compute_ng(df: pd.DataFrame, prefix: str) -> pd.Series:
    cols = [f"{prefix}_1", f"{prefix}_2", f"{prefix}_1-1", f"{prefix}_2-1"]
    return (df[cols].fillna(0) > 0).sum(axis=1).astype(int)


def annotate_v6off_cell(df: pd.DataFrame) -> pd.DataFrame:
    """Add hp_bucket_v6off + cov_bin + loh_side_norm + cell_id columns
    (the canonical partition used in Step 2)."""
    out = df.copy()
    out["hp_bucket_v6off"] = compute_hp_bucket(out, prefix="V6_off")
    out["cov_bin"] = compute_cov_bin(out["Coverage_Multiple"])
    out["loh_side_norm"] = out["loh_side"].fillna("UNKNOWN").astype(str)
    out["cell_id"] = (
        out["loh_side_norm"] + "|" + out["hp_bucket_v6off"] + "|" + out["cov_bin"]
    )
    return out


def annotate_bam_cell(df: pd.DataFrame, bam: str, flag: str) -> pd.DataFrame:
    """Add hp_bucket + cov_bin per BAM×flag (used for cell partitioning when
    we want a BAM-specific partition; for now we use V6_off-canonical partition
    for cross-BAM comparison per the plan)."""
    out = df.copy()
    prefix = f"{bam}_{flag}"
    out[f"hp_bucket_{prefix}"] = compute_hp_bucket(out, prefix=prefix)
    return out


# ---- Powered cells loader ----

def load_powered_cells() -> pd.DataFrame:
    g = pd.read_csv(STEP2_GRID, sep="\t")
    return g[g["powered_flag"] == True].reset_index(drop=True)


def load_fp_rich_cells(fp_thresh: float = 0.2055) -> pd.DataFrame:
    """Combine step3 chr8 + auto top5pct powered cells with FP_rate > thresh.

    Returns a DataFrame with columns: source, cell_id, n, n_TP, n_FP, FP_rate,
    loh_side, hp_bucket, cov, cov_axis.
    Cell partitioning differs between chr8 (genome-wide cov_bin on chr8 subset)
    and auto top5pct (FP-density top-5% regions × cov_proxy axis on those regions).
    """
    rows = []
    g_chr8 = pd.read_csv(STEP3_CHR8, sep="\t")
    for _, r in g_chr8[(g_chr8["powered"] == True) & (g_chr8["FP_rate"] > fp_thresh)].iterrows():
        rows.append({
            "source": "chr8_hotspot",
            "cell_id": f"chr8|{r['loh']}|{r['hp_bucket']}|{r['cov']}",
            "n": int(r["n"]),
            "n_TP": int(r["n_TP"]),
            "n_FP": int(r["n_FP"]),
            "FP_rate": float(r["FP_rate"]),
            "loh_side": r["loh"],
            "hp_bucket": r["hp_bucket"],
            "cov": r["cov"],
            "cov_axis": "cov_bin",
        })
    g_auto = pd.read_csv(STEP3_AUTO, sep="\t")
    for _, r in g_auto[(g_auto["powered"] == True) & (g_auto["FP_rate"] > fp_thresh)].iterrows():
        rows.append({
            "source": "auto_top5pct",
            "cell_id": f"auto|{r['loh']}|{r['hp_bucket']}|{r['cov']}",
            "n": int(r["n"]),
            "n_TP": int(r["n_TP"]),
            "n_FP": int(r["n_FP"]),
            "FP_rate": float(r["FP_rate"]),
            "loh_side": r["loh"],
            "hp_bucket": r["hp_bucket"],
            "cov": r["cov"],
            "cov_axis": "cov_proxy",
        })
    return pd.DataFrame(rows)
