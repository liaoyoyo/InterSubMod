"""Shared helpers for Step 3 FP-zone deep dive.

Reads Step 1 master TSV (full 35,332 rows, NOT filtered by master_join_ok)
to retain all FP. Step 2 caveat #1 explicit: 90% FP lost when filtering by
master_join_ok; Step 3 must use full set for FP enrichment analysis.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

# Paths
PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn")
STEP1_MASTER = PROJECT_ROOT / "step1_v3f_v5_v6_three_way" / "step1_master_three_way.tsv"
STEP3_DIR = PROJECT_ROOT / "step3_fp_zone_zoom"
SEQC2_CNV = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv"
)

# HP prefix selection (Step 2 convention)
HP_PREFIX_PRIMARY = "V6_off"
HP_PREFIXES_TRAJECTORY = ["V3F_off", "V5_off", "V6_off"]

# Power thresholds
POWERED_THRESHOLD = 50
MARGINAL_THRESHOLD = 30


def load_master_full() -> pd.DataFrame:
    """Load full Step 1 master TSV (all 35,332 rows, including unjoined FP)."""
    df = pd.read_csv(STEP1_MASTER, sep="\t", low_memory=False)
    return df


def load_seqc2_cnv() -> pd.DataFrame:
    """Load SEQC2 CNV truth for HCC1395 paired mode."""
    df = pd.read_csv(SEQC2_CNV, sep="\t", low_memory=False)
    return df


def compute_hp_bucket(df: pd.DataFrame, prefix: str = HP_PREFIX_PRIMARY) -> pd.Series:
    """Thread D 4-bucket + 'other' classification per row.

    same_HP1 / same_HP2 / cross_het / cross_het_inv / other
    """
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
    cols = [f"{prefix}_1", f"{prefix}_2", f"{prefix}_1-1", f"{prefix}_2-1"]
    return (df[cols].fillna(0) > 0).sum(axis=1)


def compute_cov_bin(s: pd.Series) -> pd.Series:
    """Coverage_Multiple → 5 bins. NaN → 'unknown'."""
    bins = [-np.inf, 0.7, 1.3, 1.8, 2.5, np.inf]
    labels = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]
    out = pd.cut(s, bins=bins, labels=labels, right=False).astype(object)
    return out.fillna("unknown")


def compute_n_reads_cov_proxy(s: pd.Series) -> pd.Series:
    """Step 3 fallback when Coverage_Multiple missing (master_join_ok==0).

    Uses V6_off_n_reads quantile (q33/q67/q90) as cov proxy.
    """
    q33, q67, q90 = s.quantile([0.33, 0.67, 0.90])
    bins = [-np.inf, q33, q67, q90, np.inf]
    labels = ["cov_proxy_low", "cov_proxy_mid", "cov_proxy_high", "cov_proxy_top"]
    return pd.cut(s, bins=bins, labels=labels, right=False).astype(object).fillna("cov_proxy_unknown")


def annotate_axes(df: pd.DataFrame, hp_prefix: str = HP_PREFIX_PRIMARY) -> pd.DataFrame:
    out = df.copy()
    out["hp_bucket"] = compute_hp_bucket(out, prefix=hp_prefix)
    out["cov_bin"] = compute_cov_bin(out["Coverage_Multiple"])
    out["loh_side_norm"] = out["loh_side"].fillna("UNKNOWN").astype(str)
    out["NG_v6off"] = compute_ng(out, prefix=hp_prefix).astype(int)
    out["n_reads_v6off"] = df[f"{hp_prefix}_n_reads"].fillna(0).astype(int)
    out["cov_proxy"] = compute_n_reads_cov_proxy(out["n_reads_v6off"])
    return out


def wilson_ci(k: int, n: int, z: float = 1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    margin = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - margin), min(1.0, center + margin))


def cell_stats(df_cell: pd.DataFrame, global_tp_rate: float, global_fp_rate: float) -> dict:
    """Compute per-cell TP/FP rate, Wilson CI, enrichment, log-odds."""
    n = len(df_cell)
    n_tp = int((df_cell["label"] == "TP").sum())
    n_fp = int((df_cell["label"] == "FP").sum())
    tp_rate = n_tp / n if n else float("nan")
    fp_rate = n_fp / n if n else float("nan")
    tp_lo, tp_hi = wilson_ci(n_tp, n)
    fp_lo, fp_hi = wilson_ci(n_fp, n)
    return dict(
        n=n,
        n_TP=n_tp,
        n_FP=n_fp,
        TP_rate=tp_rate,
        FP_rate=fp_rate,
        TP_wilson_lo=tp_lo,
        TP_wilson_hi=tp_hi,
        FP_wilson_lo=fp_lo,
        FP_wilson_hi=fp_hi,
        TP_enrichment=(tp_rate / global_tp_rate) if global_tp_rate > 0 else float("nan"),
        FP_enrichment=(fp_rate / global_fp_rate) if global_fp_rate > 0 else float("nan"),
        log_odds=(math.log(n_tp / n_fp) if (n_tp > 0 and n_fp > 0) else float("nan")),
        powered=(n >= POWERED_THRESHOLD),
        marginal=(MARGINAL_THRESHOLD <= n < POWERED_THRESHOLD),
        underpowered=(n < MARGINAL_THRESHOLD),
    )


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. Input: 1D p-values. Output: q-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    valid = ~np.isnan(p)
    q = np.full(n, np.nan)
    if valid.sum() == 0:
        return q
    order = np.argsort(p[valid])
    ranked_p = p[valid][order]
    ranked_q = ranked_p * len(ranked_p) / (np.arange(len(ranked_p)) + 1)
    # Cumulative minimum from the right
    ranked_q = np.minimum.accumulate(ranked_q[::-1])[::-1]
    ranked_q = np.clip(ranked_q, 0, 1)
    # Map back
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    q_valid = ranked_q[inv]
    q[valid] = q_valid
    return q


def deviance_decomposition(y: np.ndarray, X_full: np.ndarray, X_reduced: np.ndarray) -> dict:
    """Logistic regression deviance decomposition (replaces AUC).

    Returns dev_full, dev_reduced, dev_explained_full, lrt_p.
    """
    from scipy import stats
    from statsmodels.api import GLM, families

    try:
        m_full = GLM(y, X_full, family=families.Binomial()).fit()
        m_red = GLM(y, X_reduced, family=families.Binomial()).fit()
        d_full = float(m_full.deviance)
        d_red = float(m_red.deviance)
        d_null = float(m_full.null_deviance)
        dev_explained_full = 1 - d_full / d_null if d_null > 0 else float("nan")
        lrt_stat = d_red - d_full
        ddf = m_full.df_model - m_red.df_model
        if ddf <= 0 or lrt_stat < 0:
            lrt_p = float("nan")
        else:
            lrt_p = float(stats.chi2.sf(lrt_stat, ddf))
        return dict(
            dev_full=d_full,
            dev_reduced=d_red,
            dev_null=d_null,
            dev_explained_full=dev_explained_full,
            lrt_stat=lrt_stat,
            lrt_p=lrt_p,
        )
    except Exception as e:
        return dict(
            dev_full=float("nan"),
            dev_reduced=float("nan"),
            dev_null=float("nan"),
            dev_explained_full=float("nan"),
            lrt_stat=float("nan"),
            lrt_p=float("nan"),
            err=str(e),
        )


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


# Canonical orderings
HP_LEVELS = ["same_HP1", "same_HP2", "cross_het", "cross_het_inv", "other"]
COV_LEVELS = ["cov_loss", "cov_normal", "cov_elevated", "cov_gain", "cov_high_gain"]
COV_PROXY_LEVELS = ["cov_proxy_low", "cov_proxy_mid", "cov_proxy_high", "cov_proxy_top"]
LOH_LEVELS = ["Inner", "Outer", "UNKNOWN"]
