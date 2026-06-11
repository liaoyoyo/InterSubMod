"""Shared helpers for Step 4 cross-sample extension (mirrors Step 2 _common.py with V6-only)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn")
STEP4_DIR = PROJECT_ROOT / "step4_cross_sample_extension"
PER_SAMPLE_MASTER = STEP4_DIR / "per_sample_master"
PER_SAMPLE_GRID = STEP4_DIR / "per_sample_grid"
INTERMEDIATE_DIR = STEP4_DIR / "intermediate"
FIGURES_DIR = STEP4_DIR / "figures"

SAMPLES = ["HCC1395", "H1437", "H2009", "HCC1954", "HCC1937"]
SAMPLES_NEW = ["H1437", "H2009", "HCC1954", "HCC1937"]

# Power thresholds (same as Step 2)
POWERED_THRESHOLD = 50
MARGINAL_THRESHOLD = 30

# Canonical orderings (LOH × HP × CN)
LOH_LEVELS = ["Inner", "Outer"]
HP_LEVELS = ["same_HP1", "same_HP2", "cross_het", "cross_het_inv", "other"]
COV_LEVELS = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]


def compute_hp_bucket(df: pd.DataFrame, prefix: str = "V6_off") -> pd.Series:
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


def compute_ng(df: pd.DataFrame, prefix: str = "V6_off") -> pd.Series:
    cols = [f"{prefix}_1", f"{prefix}_2", f"{prefix}_1-1", f"{prefix}_2-1"]
    return (df[cols].fillna(0) > 0).sum(axis=1)


def compute_cov_bin(s: pd.Series) -> pd.Series:
    bins = [-np.inf, 0.7, 1.3, 1.8, 2.5, np.inf]
    labels = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]
    return pd.cut(s, bins=bins, labels=labels, right=False).astype(object).fillna("unknown")


def load_sample_master(sample: str) -> pd.DataFrame:
    """Load per-sample master TSV (Step 4 stage 1 output) for the 4 new samples,
    OR reconstruct the HCC1395 view from Step 1 master_three_way for unified analysis.
    """
    if sample == "HCC1395":
        path = (
            PROJECT_ROOT
            / "step1_v3f_v5_v6_three_way"
            / "step1_master_three_way.tsv"
        )
        df = pd.read_csv(path, sep="\t", low_memory=False)
        # Step 1 already has V6_off_* / V6_on_* columns named exactly the same
        return df
    path = PER_SAMPLE_MASTER / f"{sample}_v6_master.tsv"
    return pd.read_csv(path, sep="\t", low_memory=False)


def annotate_axes(df: pd.DataFrame, prefix: str = "V6_off") -> pd.DataFrame:
    out = df.copy()
    out["hp_bucket"] = compute_hp_bucket(out, prefix=prefix)
    out["cov_bin"] = compute_cov_bin(
        pd.to_numeric(out["Coverage_Multiple"], errors="coerce")
    )
    out["loh_side_norm"] = out["loh_side"].fillna("UNKNOWN").astype(str)
    out["NG_v6off"] = compute_ng(out, prefix=prefix).astype(int)
    out["n_reads_v6off"] = pd.to_numeric(
        out[f"{prefix}_n_reads"], errors="coerce"
    ).fillna(0).astype(int)
    out["caller_af_num"] = pd.to_numeric(out["caller_af"], errors="coerce")
    return out


def joined_subset(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["master_join_ok"] == 1].copy()


def wilson_ci(k: int, n: int, z: float = 1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    margin = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (center - margin, center + margin)


def cell_id(loh: str, hp: str, cov: str) -> str:
    return f"{loh}|{hp}|{cov}"


def all_main_cells():
    """50-cell main grid (Inner/Outer × 5 HP × 5 cov, excludes UNKNOWN)."""
    return [(l, h, c) for l in LOH_LEVELS for h in HP_LEVELS for c in COV_LEVELS]


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
