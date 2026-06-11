"""Shared helpers for Step 2 scripts (loading master TSV, HP bucket, cov_bin, NG, etc.)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

# Project root
PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn")
STEP1_MASTER = PROJECT_ROOT / "step1_v3f_v5_v6_three_way" / "step1_master_three_way.tsv"
STEP2_DIR = PROJECT_ROOT / "step2_three_axis_grid"
INTERMEDIATE_DIR = STEP2_DIR / "intermediate"
FIGURES_DIR = STEP2_DIR / "figures"
BUILD_LOG = INTERMEDIATE_DIR / "build_log.txt"

# HP bucket / cov_bin / NG computation prefixes
HP_PREFIX_PRIMARY = "V6_off"  # primary BAM/flag for main analysis
HP_PREFIX_SENSITIVITY = "V6_on"  # sensitivity check
HP_PREFIXES_TRAJECTORY = ["V3F_off", "V5_off", "V6_off"]

# Power thresholds (per 02_methodology_notes §2)
POWERED_THRESHOLD = 50
MARGINAL_THRESHOLD = 30


def log_msg(msg: str, also_print: bool = True) -> None:
    """Append a timestamped message to the build log."""
    import datetime as _dt

    BUILD_LOG.parent.mkdir(parents=True, exist_ok=True)
    ts = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    if also_print:
        print(line)
    with open(BUILD_LOG, "a") as f:
        f.write(line + "\n")


def load_master() -> pd.DataFrame:
    """Load Step 1 master TSV (full)."""
    df = pd.read_csv(STEP1_MASTER, sep="\t", low_memory=False)
    return df


def joined_subset(df: pd.DataFrame) -> pd.DataFrame:
    """Return only master_join_ok==1 rows (Step 2 analysis universe)."""
    return df[df["master_join_ok"] == 1].copy()


def compute_hp_bucket(df: pd.DataFrame, prefix: str = HP_PREFIX_PRIMARY) -> pd.Series:
    """Compute Thread D 4-bucket + other classification per row."""
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


def compute_ng(df: pd.DataFrame, prefix: str = HP_PREFIX_PRIMARY) -> pd.Series:
    """Compute HPFineNGroups (number of non-zero HP family buckets, 0-4)."""
    cols = [f"{prefix}_1", f"{prefix}_2", f"{prefix}_1-1", f"{prefix}_2-1"]
    return (df[cols].fillna(0) > 0).sum(axis=1)


def compute_cov_bin(s: pd.Series) -> pd.Series:
    """Coverage_Multiple → 5 bins. <0.7 / 0.7-1.3 / 1.3-1.8 / 1.8-2.5 / >=2.5."""
    bins = [-np.inf, 0.7, 1.3, 1.8, 2.5, np.inf]
    labels = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]
    return pd.cut(s, bins=bins, labels=labels, right=False).astype(object).fillna("unknown")


def normalise_loh_side(s: pd.Series) -> pd.Series:
    return s.fillna("UNKNOWN").astype(str)


def annotate_axes(
    df: pd.DataFrame, hp_prefix: str = HP_PREFIX_PRIMARY
) -> pd.DataFrame:
    """Add loh_side / hp_bucket / cov_bin / NG / n_reads columns derived from chosen prefix."""
    out = df.copy()
    out["hp_bucket"] = compute_hp_bucket(out, prefix=hp_prefix)
    out["cov_bin"] = compute_cov_bin(out["Coverage_Multiple"])
    out["loh_side_norm"] = normalise_loh_side(out["loh_side"])
    out["NG_v6off"] = compute_ng(out, prefix=hp_prefix).astype(int)
    out["n_reads_v6off"] = df[f"{hp_prefix}_n_reads"].fillna(0).astype(int)
    return out


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


# Canonical cell axis orderings (consistent across all tables / figures)
LOH_LEVELS = ["Inner", "Outer", "UNKNOWN"]
HP_LEVELS = ["same_HP1", "same_HP2", "cross_het", "cross_het_inv", "other"]
COV_LEVELS = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]


def all_cells():
    return [(l, h, c) for l in LOH_LEVELS for h in HP_LEVELS for c in COV_LEVELS]


def matplotlib_setup():
    """Configure matplotlib with Droid Sans Fallback for CJK support."""
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
