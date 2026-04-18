#!/usr/bin/env python3
"""Shared utilities for O1-O10 systematic observation scripts.

Provides:
  - Data loading (master dataset, read training tables)
  - Plot style & save helpers
  - Statistical tools (AUC, effect size, enrichment, bootstrap CI)
  - Workspace output helpers (context JSON, feature stats, report scaffold)

Usage:
    from observation_common import (
        load_master_dataset, setup_plot_style, save_figure,
        compute_auc, compute_effect_size, compute_enrichment,
        bootstrap_ci, write_round_context, write_feature_statistics,
        ensure_dir, safe_div,
    )
"""

from __future__ import annotations

import json
import math
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MASTER_DATASET_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260327_loh_round1_cross_sample_audit/"
    "all_region_rows.tsv.gz"
)

OUTPUT_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces"
)

SAMPLE_ORDER = [
    "HCC1395", "HCC1395_DORADO", "COLO829",
    "H1437", "H2009", "HCC1937", "HCC1954",
]

MODE_ORDER = ["paired", "to"]

# Colour palette
COLOR_TP = "#2196F3"
COLOR_FP = "#F44336"
COLOR_PAIRED = "#1565C0"
COLOR_TO = "#E65100"
COLOR_NEUTRAL = "#757575"

TRUTH_PALETTE = {"TP": COLOR_TP, "FP": COLOR_FP}
MODE_PALETTE = {"paired": COLOR_PAIRED, "to": COLOR_TO}

# ---------------------------------------------------------------------------
# Filesystem helpers
# ---------------------------------------------------------------------------


def ensure_dir(path: Path) -> Path:
    """Create directory (and parents) if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_div(num: float, den: float, default: float = 0.0) -> float:
    """Safe division avoiding ZeroDivisionError."""
    return num / den if den != 0 else default


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_master_dataset(
    path: Optional[Path] = None,
    usecols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Load the 748K-row master dataset (gzipped TSV).

    Parameters
    ----------
    path : Path, optional
        Override the default master dataset path.
    usecols : list of str, optional
        Only load these columns (faster for large-column subsets).

    Returns
    -------
    pd.DataFrame
    """
    p = path or MASTER_DATASET_PATH
    print(f"[observation_common] Loading master dataset from {p} ...")
    df = pd.read_csv(p, sep="\t", compression="gzip", usecols=usecols, low_memory=False)
    print(f"[observation_common]   -> {len(df):,} rows x {len(df.columns)} cols")
    return df


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------


def setup_plot_style() -> None:
    """Configure matplotlib/seaborn defaults for observation figures."""
    sns.set_theme(style="whitegrid", font_scale=1.1)
    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": 180,
        "savefig.bbox": "tight",
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 14,
    })


def save_figure(
    fig: plt.Figure,
    path: Union[str, Path],
    dpi: int = 180,
    close: bool = True,
) -> Path:
    """Save figure to PNG and optionally close it."""
    p = Path(path)
    ensure_dir(p.parent)
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor="white")
    if close:
        plt.close(fig)
    print(f"[observation_common] Saved figure: {p.name} ({p.stat().st_size / 1024:.0f} KB)")
    return p


def add_caption(
    fig: plt.Figure,
    text: str,
    fontsize: int = 7,
    y: float = -0.02,
) -> None:
    """Add a data-source caption at the bottom of the figure."""
    fig.text(0.5, y, text, ha="center", fontsize=fontsize, color="#999999")


# ---------------------------------------------------------------------------
# Statistical tools
# ---------------------------------------------------------------------------


def compute_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute ROC AUC. Returns NaN if fewer than 2 classes present."""
    from sklearn.metrics import roc_auc_score  # lazy import

    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    mask = np.isfinite(y_score) & np.isfinite(y_true)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 10:
        return float("nan")
    try:
        return roc_auc_score(y_true, y_score)
    except Exception:
        return float("nan")


def compute_effect_size(
    group1: np.ndarray,
    group2: np.ndarray,
    n_boot: int = 1000,
    ci: float = 0.95,
) -> Dict[str, Any]:
    """Compute Cohen's d with bootstrap 95% CI.

    Returns dict with keys: d, ci_lo, ci_hi, n1, n2.
    """
    g1 = np.asarray(group1, dtype=float)
    g2 = np.asarray(group2, dtype=float)
    g1 = g1[np.isfinite(g1)]
    g2 = g2[np.isfinite(g2)]
    n1, n2 = len(g1), len(g2)
    if n1 < 2 or n2 < 2:
        return {"d": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"), "n1": n1, "n2": n2}

    pooled_std = np.sqrt(((n1 - 1) * g1.std(ddof=1) ** 2 + (n2 - 1) * g2.std(ddof=1) ** 2) / (n1 + n2 - 2))
    d = (g1.mean() - g2.mean()) / pooled_std if pooled_std > 0 else 0.0

    # Bootstrap CI
    rng = np.random.RandomState(42)
    boot_ds = []
    for _ in range(n_boot):
        b1 = rng.choice(g1, size=n1, replace=True)
        b2 = rng.choice(g2, size=n2, replace=True)
        ps = np.sqrt(((n1 - 1) * b1.std(ddof=1) ** 2 + (n2 - 1) * b2.std(ddof=1) ** 2) / (n1 + n2 - 2))
        boot_ds.append((b1.mean() - b2.mean()) / ps if ps > 0 else 0.0)
    alpha = (1 - ci) / 2
    ci_lo = float(np.percentile(boot_ds, 100 * alpha))
    ci_hi = float(np.percentile(boot_ds, 100 * (1 - alpha)))

    return {"d": float(d), "ci_lo": ci_lo, "ci_hi": ci_hi, "n1": n1, "n2": n2}


def compute_enrichment(
    tp_count: int,
    fp_count: int,
    tp_total: int,
    fp_total: int,
) -> Dict[str, Any]:
    """Compute enrichment ratio (FP rate / TP rate) with Fisher exact test.

    enrichment > 1 means FP-enriched; < 1 means TP-enriched.
    Returns dict with keys: enrichment, p_value, tp_rate, fp_rate.
    """
    tp_rate = safe_div(tp_count, tp_total)
    fp_rate = safe_div(fp_count, fp_total)
    enrichment = safe_div(fp_rate, tp_rate, default=float("nan"))

    # Fisher 2x2: [[fp_count, fp_total-fp_count], [tp_count, tp_total-tp_count]]
    table = np.array([
        [fp_count, fp_total - fp_count],
        [tp_count, tp_total - tp_count],
    ])
    try:
        _, p_value = sp_stats.fisher_exact(table)
    except Exception:
        p_value = float("nan")

    return {
        "enrichment": float(enrichment),
        "p_value": float(p_value),
        "tp_rate": float(tp_rate),
        "fp_rate": float(fp_rate),
        "tp_count": tp_count,
        "fp_count": fp_count,
    }


def bootstrap_ci(
    statistic_fn: Callable[[np.ndarray], float],
    data: np.ndarray,
    n_boot: int = 1000,
    ci: float = 0.95,
    seed: int = 42,
) -> Tuple[float, float, float]:
    """Bootstrap confidence interval for a statistic.

    Returns (point_estimate, ci_lo, ci_hi).
    """
    data = np.asarray(data)
    data = data[np.isfinite(data)]
    if len(data) < 2:
        return (float("nan"), float("nan"), float("nan"))

    point = statistic_fn(data)
    rng = np.random.RandomState(seed)
    boot_vals = [statistic_fn(rng.choice(data, size=len(data), replace=True)) for _ in range(n_boot)]
    alpha = (1 - ci) / 2
    lo = float(np.percentile(boot_vals, 100 * alpha))
    hi = float(np.percentile(boot_vals, 100 * (1 - alpha)))
    return (float(point), lo, hi)


def mannwhitney_test(
    group1: np.ndarray,
    group2: np.ndarray,
) -> Dict[str, Any]:
    """Mann-Whitney U test with effect size (rank-biserial correlation).

    Returns dict: statistic, p_value, effect_size_r, n1, n2.
    """
    g1 = np.asarray(group1, dtype=float)
    g2 = np.asarray(group2, dtype=float)
    g1 = g1[np.isfinite(g1)]
    g2 = g2[np.isfinite(g2)]
    n1, n2 = len(g1), len(g2)
    if n1 < 2 or n2 < 2:
        return {"statistic": float("nan"), "p_value": float("nan"),
                "effect_size_r": float("nan"), "n1": n1, "n2": n2}
    try:
        u_stat, p_val = sp_stats.mannwhitneyu(g1, g2, alternative="two-sided")
        # rank-biserial correlation r = 1 - 2U/(n1*n2)
        r = 1 - 2 * u_stat / (n1 * n2)
    except Exception:
        u_stat, p_val, r = float("nan"), float("nan"), float("nan")
    return {"statistic": float(u_stat), "p_value": float(p_val),
            "effect_size_r": float(r), "n1": n1, "n2": n2}


def chi2_test(
    contingency_table: np.ndarray,
) -> Dict[str, Any]:
    """Chi-squared test with Cramer's V.

    Parameters
    ----------
    contingency_table : 2D array-like

    Returns dict: chi2, p_value, cramers_v, dof.
    """
    table = np.asarray(contingency_table)
    try:
        chi2, p, dof, _ = sp_stats.chi2_contingency(table)
        n = table.sum()
        k = min(table.shape) - 1
        v = np.sqrt(chi2 / (n * k)) if n * k > 0 else 0.0
    except Exception:
        chi2, p, dof, v = float("nan"), float("nan"), 0, float("nan")
    return {"chi2": float(chi2), "p_value": float(p), "cramers_v": float(v), "dof": int(dof)}


# ---------------------------------------------------------------------------
# F1 / precision / recall helpers
# ---------------------------------------------------------------------------


def compute_metrics(tp: int, fp: int, fn: int) -> Dict[str, float]:
    """Compute precision, recall, F1 from counts."""
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    f1 = safe_div(2 * precision * recall, precision + recall)
    return {"precision": precision, "recall": recall, "f1": f1, "tp": tp, "fp": fp, "fn": fn}


# ---------------------------------------------------------------------------
# Workspace output helpers
# ---------------------------------------------------------------------------


def write_round_context(
    out_dir: Path,
    observation_id: str,
    title: str,
    description: str,
    script_path: str,
    data_sources: List[Dict[str, Any]],
    row_count: int,
    col_count: int,
    extra: Optional[Dict[str, Any]] = None,
) -> Path:
    """Write round_context.json with standard metadata."""
    ctx = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S CST"),
        "workspace_type": "systematic_observation",
        "observation_id": observation_id,
        "title": title,
        "description": description,
        "script": script_path,
        "data_sources": data_sources,
        "row_count": row_count,
        "col_count": col_count,
        "samples": SAMPLE_ORDER,
        "modes": MODE_ORDER,
    }
    if extra:
        ctx.update(extra)

    path = out_dir / "round_context.json"
    ensure_dir(out_dir)
    with path.open("w", encoding="utf-8") as f:
        json.dump(ctx, f, indent=2, ensure_ascii=False)
    return path


def write_feature_statistics(
    df: pd.DataFrame,
    columns: List[str],
    out_path: Path,
) -> pd.DataFrame:
    """Write descriptive statistics TSV for given columns.

    Columns in output: feature, dtype, count, missing, missing_pct,
    mean, std, min, p25, median, p75, max, nunique.
    """
    rows = []
    for col in columns:
        if col not in df.columns:
            continue
        s = df[col]
        missing = int(s.isna().sum())
        missing_pct = safe_div(missing, len(s)) * 100
        row = {
            "feature": col,
            "dtype": str(s.dtype),
            "count": len(s),
            "missing": missing,
            "missing_pct": round(missing_pct, 2),
            "nunique": int(s.nunique()),
        }
        if pd.api.types.is_numeric_dtype(s):
            desc = s.describe()
            row.update({
                "mean": round(float(desc.get("mean", float("nan"))), 4),
                "std": round(float(desc.get("std", float("nan"))), 4),
                "min": round(float(desc.get("min", float("nan"))), 4),
                "p25": round(float(desc.get("25%", float("nan"))), 4),
                "median": round(float(desc.get("50%", float("nan"))), 4),
                "p75": round(float(desc.get("75%", float("nan"))), 4),
                "max": round(float(desc.get("max", float("nan"))), 4),
            })
        else:
            top = s.value_counts().head(3)
            row["top_values"] = "; ".join(f"{v}({c})" for v, c in top.items())
        rows.append(row)

    stats_df = pd.DataFrame(rows)
    ensure_dir(out_path.parent)
    stats_df.to_csv(out_path, sep="\t", index=False)
    print(f"[observation_common] Feature statistics: {out_path.name} ({len(rows)} features)")
    return stats_df


def format_p(p: float) -> str:
    """Format p-value for display."""
    if math.isnan(p):
        return "NaN"
    if p < 1e-300:
        return "< 1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def format_ci(lo: float, hi: float, digits: int = 3) -> str:
    """Format confidence interval for display."""
    return f"[{lo:.{digits}f}, {hi:.{digits}f}]"


# ---------------------------------------------------------------------------
# Truth label encoding
# ---------------------------------------------------------------------------


def encode_truth_binary(df: pd.DataFrame, col: str = "truth_label") -> pd.Series:
    """Encode truth_label to binary: TP=1, FP=0. NaN for others."""
    mapping = {"TP": 1, "FP": 0}
    return df[col].map(mapping)


# ---------------------------------------------------------------------------
# QS component reconstruction (for O2)
# ---------------------------------------------------------------------------


def reconstruct_qs_components(df: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct Quality Score components from master dataset columns.

    Mirrors the C++ logic in RegionProcessor.cpp:109-157.
    Returns df with new columns: qs_penalty_*, qs_bonus_*, qs_reconstructed.
    """
    out = pd.DataFrame(index=df.index)

    nr = df["NumReads"] if "NumReads" in df.columns else df.get("num_reads", pd.Series(dtype=float))
    nc = df["NumCpGs"] if "NumCpGs" in df.columns else df.get("num_cpgs", pd.Series(dtype=float))
    cov = df.get("Coverage_Multiple", pd.Series(dtype=float))
    loh = df.get("Potential_LOH", pd.Series(dtype=float))
    hp_sig = df.get("HPMergedSig", pd.Series(dtype=float))
    allele_sig = df.get("AlleleSig", pd.Series(dtype=float))
    cv = df.get("CramersV", pd.Series(dtype=float))

    # Penalties
    out["qs_penalty_low_reads"] = np.where(nr < 30, -20.0, np.where(nr < 50, -10.0, 0.0))
    out["qs_penalty_low_cpgs"] = np.where(nc < 20, -15.0, np.where(nc < 30, -10.0, 0.0))
    out["qs_penalty_cnv_loss"] = np.where(cov < 0.5, -30.0, np.where(cov < 0.8, -15.0, 0.0))
    out["qs_penalty_high_copy"] = np.where(cov > 2.0, -20.0, 0.0)
    out["qs_penalty_loh"] = np.where(
        loh.astype(str).str.lower().isin(["true", "1", "yes"]), -25.0, 0.0
    )

    # Bonuses
    hp_sig_bool = hp_sig.astype(str).str.lower().isin(["true", "1", "yes"])
    allele_sig_bool = allele_sig.astype(str).str.lower().isin(["true", "1", "yes"])
    out["qs_bonus_dual_sig"] = np.where(hp_sig_bool & allele_sig_bool, 10.0, 0.0)
    out["qs_bonus_strong_effect"] = np.where(cv >= 0.3, 15.0, np.where(cv >= 0.2, 5.0, 0.0))

    # Reconstruct
    out["qs_reconstructed"] = (
        100.0
        + out["qs_penalty_low_reads"]
        + out["qs_penalty_low_cpgs"]
        + out["qs_penalty_cnv_loss"]
        + out["qs_penalty_high_copy"]
        + out["qs_penalty_loh"]
        + out["qs_bonus_dual_sig"]
        + out["qs_bonus_strong_effect"]
    ).clip(0, 100)

    return out


# ---------------------------------------------------------------------------
# Markdown helpers
# ---------------------------------------------------------------------------


def markdown_table(columns: Sequence[str], rows: Sequence[Dict]) -> str:
    """Generate a Markdown pipe table."""
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, sep]
    for row in rows:
        vals = [str(row.get(c, "")) for c in columns]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)
