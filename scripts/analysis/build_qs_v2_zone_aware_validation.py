#!/usr/bin/env python3
"""QS v2 Zone-Aware Validation: Compare current vs LOH-zone-aware Quality Score.

Validates whether a LOH zone-aware QS redesign improves discriminative power
(AUC) for tumor-only (TO) mode, where current QS has AUC=0.497 (random).

Three formulas compared:
  - QS_v1: current production formula (mode-aware, mirrors C++ exactly)
  - QS_v2a: zone-aware with caller_ad_ref/ad_alt (needs C++ VCF change)
  - QS_v2b: zone-aware with ISM-only fields (no C++ VCF change needed)

Usage:
    python3 scripts/analysis/build_qs_v2_zone_aware_validation.py
"""

from __future__ import annotations

import sys
import time
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure, add_caption,
    compute_auc, encode_truth_binary, ensure_dir,
    SAMPLE_ORDER, MODE_ORDER,
    COLOR_TP, COLOR_FP,
)

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

# ── Output Paths ──────────────────────────────────────────────────────────
BASE_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation")
FIG_DIR = ensure_dir(BASE_DIR / "figures")
DATA_DIR = ensure_dir(BASE_DIR / "data")
REPORT_DIR = ensure_dir(BASE_DIR / "reports")

PREFIX = "qs_v2"

CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz (7 samples) | "
    "LOH zones: per-sample LongPhase-TO LOH.bed (binary inside/outside) | "
    "Script: build_qs_v2_zone_aware_validation.py"
)

# ── LOH.bed Paths ──────────────────────────────────────────────────────────
LOH_BED_PATHS = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "COLO829": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260317_colo829_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed",
    "H1437": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "H2009": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1937": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1954": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
}

# ── Zone Colors ──────────────────────────────────────────────────────────
COLOR_INSIDE = "#E53935"   # red
COLOR_OUTSIDE = "#1E88E5"  # blue
ZONE_ORDER = ["inside", "outside"]

# ── QS formula colors ───────────────────────────────────────────────────
COLOR_V1 = "#9E9E9E"   # gray
COLOR_V2A = "#1976D2"  # blue
COLOR_V2B = "#388E3C"  # green
QS_COLORS = {"v1": COLOR_V1, "v2a": COLOR_V2A, "v2b": COLOR_V2B}
QS_LABELS = {"v1": "QS_v1 (current)", "v2a": "QS_v2a (AD-based)", "v2b": "QS_v2b (ISM-only)"}

# ── Columns to load ─────────────────────────────────────────────────────
LOAD_COLS = [
    "RegionID", "Chr", "Pos", "sample", "mode", "truth_label",
    "NumReads", "NumCpGs", "Coverage_Multiple", "Potential_LOH",
    "HPMergedSig", "AlleleSig", "CramersV",
    "HP_Ratio", "Quality_Score",
    "caller_af", "caller_gq", "caller_dp", "caller_ad_ref", "caller_ad_alt",
    "HPFineF", "LabelAllelePermanovaF",
    "to_loh_bed_hit",
]


# ══════════════════════════════════════════════════════════════════════════
# LOH Zone Assignment (Binary — per-sample LOH.bed)
# ══════════════════════════════════════════════════════════════════════════

def load_bed_intervals(bed_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """Load BED intervals sorted by start for binary search."""
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def assign_loh_binary_vectorized(
    df: pd.DataFrame,
    intervals: Dict[str, List[Tuple[int, int]]],
) -> pd.Series:
    """Assign binary LOH zone (inside/outside) using vectorized search."""
    zones = pd.Series("outside", index=df.index, dtype="object")
    for chrom, grp in df.groupby("Chr"):
        if chrom not in intervals:
            continue
        ivs = intervals[chrom]
        if len(ivs) == 0:
            continue
        starts = np.array([s for s, e in ivs])
        ends = np.array([e for s, e in ivs])
        positions = grp["Pos"].values
        insert_idx = np.searchsorted(starts, positions, side="right") - 1
        in_range = (insert_idx >= 0) & (insert_idx < len(ivs))
        mask = np.zeros(len(positions), dtype=bool)
        valid = np.where(in_range)[0]
        if len(valid) > 0:
            valid_ins = insert_idx[valid]
            mask[valid] = (positions[valid] >= starts[valid_ins]) & (positions[valid] < ends[valid_ins])
        # Use .loc with actual index values (not iloc with positional indices)
        zones.loc[grp.index[mask]] = "inside"
    return zones


def assign_loh_zones_all_samples(df: pd.DataFrame) -> pd.DataFrame:
    """Assign LOH zones per sample using each sample's own LOH.bed."""
    print("[QS-v2] Assigning LOH zones per sample (binary)...")
    df["loh_zone"] = "outside"

    for sample in SAMPLE_ORDER:
        t0 = time.time()
        bed_path = LOH_BED_PATHS.get(sample)
        if bed_path is None or not Path(bed_path).exists():
            print(f"  [WARN] No LOH.bed for {sample} -- all marked 'outside'")
            continue
        intervals = load_bed_intervals(bed_path)
        mask_sample = df["sample"] == sample
        sub = df.loc[mask_sample]
        zones = assign_loh_binary_vectorized(sub, intervals)
        df.loc[mask_sample, "loh_zone"] = zones.values
        n_in = (zones == "inside").sum()
        n_out = (zones == "outside").sum()
        elapsed = time.time() - t0
        print(f"  {sample}: inside={n_in:,}, outside={n_out:,} ({elapsed:.1f}s)")

    # Cross-check for TO rows: compare with existing to_loh_bed_hit
    to_rows = df[(df["mode"] == "to") & df["to_loh_bed_hit"].notna()]
    if len(to_rows) > 0:
        to_hit_bool = to_rows["to_loh_bed_hit"].astype(str).str.lower().isin(["true", "1", "yes"])
        our_inside = to_rows["loh_zone"] == "inside"
        concordance = (to_hit_bool == our_inside).mean()
        print(f"[QS-v2] Cross-check to_loh_bed_hit concordance: {concordance:.4f} "
              f"({(to_hit_bool == our_inside).sum()}/{len(to_rows)})")

    # Summary
    zone_summary = df.groupby(["sample", "mode", "loh_zone"]).size().unstack(fill_value=0)
    print("[QS-v2] Zone assignment summary:")
    print(zone_summary.to_string())
    return df


# ══════════════════════════════════════════════════════════════════════════
# QS Formulas
# ══════════════════════════════════════════════════════════════════════════

def _bool_col(df: pd.DataFrame, col: str) -> np.ndarray:
    """Convert column to boolean array, handling mixed dtypes."""
    s = df[col]
    if s.dtype == bool:
        return s.values
    return s.astype(str).str.lower().isin(["true", "1", "yes"]).values


def compute_qs_v1(df: pd.DataFrame) -> np.ndarray:
    """Vectorized QS v1 matching C++ RegionProcessor.cpp:133-183.

    Mode-aware: TO has loh_penalty=0 and verify_bonus=0.
    """
    score = np.full(len(df), 100.0)
    nr = df["NumReads"].values.astype(float)
    nc = df["NumCpGs"].values.astype(float)
    cov = df["Coverage_Multiple"].values.astype(float)

    # Read penalty
    score = np.where(nr < 30, score - 20, np.where(nr < 50, score - 10, score))
    # CpG penalty
    score = np.where(nc < 20, score - 15, np.where(nc < 30, score - 10, score))
    # Coverage penalty
    score = np.where(cov < 0.5, score - 30,
            np.where(cov < 0.8, score - 15,
            np.where(cov > 2.0, score - 20, score)))

    is_paired = (df["mode"] == "paired").values
    loh = _bool_col(df, "Potential_LOH")
    hp_sig = _bool_col(df, "HPMergedSig")
    allele_sig = _bool_col(df, "AlleleSig")
    cv = df["CramersV"].values.astype(float)

    # Paired-only: LOH penalty (TO has loh_penalty=0)
    score = np.where(is_paired & loh, score - 25, score)
    # Paired-only: verify bonus (TO has verify_bonus=0)
    score = np.where(is_paired & hp_sig & allele_sig, score + 10, score)
    # Effect bonus (both modes)
    score = np.where(cv >= 0.3, score + 15, np.where(cv >= 0.2, score + 5, score))

    return np.clip(score, 0, 100)


def compute_qs_v2a(df: pd.DataFrame) -> np.ndarray:
    """Zone-aware QS v2a: uses caller_ad_ref/ad_alt inside LOH zone.

    Inside LOH zone:
      - AD ratio bonus: ad_alt/(ad_ref+ad_alt) >= 0.15 -> +10; < 0.05 -> -15
      - HPFineF bonus: >= 5.0 -> +10; >= 2.0 -> +5
      - No LOH penalty, no verify bonus
    Outside LOH zone: same as v1 for that mode
    """
    # Start with base score (shared penalties)
    score = np.full(len(df), 100.0)
    nr = df["NumReads"].values.astype(float)
    nc = df["NumCpGs"].values.astype(float)
    cov = df["Coverage_Multiple"].values.astype(float)

    # Read penalty
    score = np.where(nr < 30, score - 20, np.where(nr < 50, score - 10, score))
    # CpG penalty
    score = np.where(nc < 20, score - 15, np.where(nc < 30, score - 10, score))
    # Coverage penalty
    score = np.where(cov < 0.5, score - 30,
            np.where(cov < 0.8, score - 15,
            np.where(cov > 2.0, score - 20, score)))

    is_paired = (df["mode"] == "paired").values
    is_inside = (df["loh_zone"] == "inside").values
    loh = _bool_col(df, "Potential_LOH")
    hp_sig = _bool_col(df, "HPMergedSig")
    allele_sig = _bool_col(df, "AlleleSig")
    cv = df["CramersV"].values.astype(float)

    ad_ref = df["caller_ad_ref"].values.astype(float)
    ad_alt = df["caller_ad_alt"].values.astype(float)
    ad_total = ad_ref + ad_alt
    ad_ratio = np.where(ad_total > 0, ad_alt / ad_total, 0.0)
    hpf = df["HPFineF"].values.astype(float)

    # --- Outside LOH zone: same as v1 ---
    outside = ~is_inside
    # Paired-only LOH penalty (outside)
    score = np.where(outside & is_paired & loh, score - 25, score)
    # Paired-only verify bonus (outside)
    score = np.where(outside & is_paired & hp_sig & allele_sig, score + 10, score)
    # TO outside: no loh_penalty, no verify_bonus (same as v1 TO)
    # Effect bonus (outside, both modes)
    score = np.where(outside & (cv >= 0.3), score + 15,
            np.where(outside & (cv >= 0.2), score + 5, score))

    # --- Inside LOH zone: zone-aware adjustments ---
    # No LOH penalty, no verify bonus (for either mode inside LOH)
    # AD ratio bonus
    score = np.where(is_inside & (ad_ratio >= 0.15), score + 10,
            np.where(is_inside & (ad_ratio < 0.05), score - 15, score))
    # HPFineF bonus
    score = np.where(is_inside & (hpf >= 5.0), score + 10,
            np.where(is_inside & (hpf >= 2.0), score + 5, score))

    return np.clip(score, 0, 100)


def compute_qs_v2b(df: pd.DataFrame) -> np.ndarray:
    """Zone-aware QS v2b: ISM-only fields, no C++ VCF change needed.

    Inside LOH zone:
      - AF balance bonus: caller_af in [0.3, 0.7] -> +10; >= 0.15 -> +5
      - HPFineF bonus: >= 5.0 -> +10; >= 2.0 -> +5
      - LabelAllelePermanovaF bonus: >= 3.0 -> +8
      - No LOH penalty, no verify bonus
    Outside LOH zone: same as v1 for that mode
    """
    score = np.full(len(df), 100.0)
    nr = df["NumReads"].values.astype(float)
    nc = df["NumCpGs"].values.astype(float)
    cov = df["Coverage_Multiple"].values.astype(float)

    # Read penalty
    score = np.where(nr < 30, score - 20, np.where(nr < 50, score - 10, score))
    # CpG penalty
    score = np.where(nc < 20, score - 15, np.where(nc < 30, score - 10, score))
    # Coverage penalty
    score = np.where(cov < 0.5, score - 30,
            np.where(cov < 0.8, score - 15,
            np.where(cov > 2.0, score - 20, score)))

    is_paired = (df["mode"] == "paired").values
    is_inside = (df["loh_zone"] == "inside").values
    loh = _bool_col(df, "Potential_LOH")
    hp_sig = _bool_col(df, "HPMergedSig")
    allele_sig = _bool_col(df, "AlleleSig")
    cv = df["CramersV"].values.astype(float)

    caf = df["caller_af"].values.astype(float)
    hpf = df["HPFineF"].values.astype(float)
    lapf = df["LabelAllelePermanovaF"].values.astype(float)

    # --- Outside LOH zone: same as v1 ---
    outside = ~is_inside
    score = np.where(outside & is_paired & loh, score - 25, score)
    score = np.where(outside & is_paired & hp_sig & allele_sig, score + 10, score)
    score = np.where(outside & (cv >= 0.3), score + 15,
            np.where(outside & (cv >= 0.2), score + 5, score))

    # --- Inside LOH zone: zone-aware ISM-only adjustments ---
    # AF balance bonus
    score = np.where(is_inside & (caf >= 0.3) & (caf <= 0.7), score + 10,
            np.where(is_inside & (caf >= 0.15), score + 5, score))
    # HPFineF bonus
    score = np.where(is_inside & (hpf >= 5.0), score + 10,
            np.where(is_inside & (hpf >= 2.0), score + 5, score))
    # LabelAllelePermanovaF bonus
    score = np.where(is_inside & (lapf >= 3.0), score + 8, score)

    return np.clip(score, 0, 100)


# ══════════════════════════════════════════════════════════════════════════
# AUC Helpers
# ══════════════════════════════════════════════════════════════════════════

def _auc_safe(y_true, y_score):
    """ROC AUC with graceful fallback."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_score)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 10:
        return float("nan")
    return compute_auc(y_true, y_score)


# ══════════════════════════════════════════════════════════════════════════
# Figure Generation
# ══════════════════════════════════════════════════════════════════════════

def fig01_auc_by_zone(df: pd.DataFrame, y_bin: np.ndarray) -> Path:
    """Grouped bar chart: AUC for v1/v2a/v2b grouped by {mode x LOH_zone}."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    groups = []
    for mode in ["to", "paired"]:
        for zone in ["inside", "outside"]:
            groups.append((mode, zone))

    formulas = ["v1", "v2a", "v2b"]
    qs_cols = {"v1": "QS_v1", "v2a": "QS_v2a", "v2b": "QS_v2b"}
    x = np.arange(len(groups))
    width = 0.25

    for i, formula in enumerate(formulas):
        aucs = []
        for mode, zone in groups:
            mask = (df["mode"] == mode) & (df["loh_zone"] == zone) & np.isfinite(y_bin)
            auc_val = _auc_safe(y_bin[mask], df.loc[mask, qs_cols[formula]].values)
            aucs.append(auc_val)
        bars = ax.bar(x + i * width, aucs, width, label=QS_LABELS[formula],
                      color=QS_COLORS[formula], edgecolor="white", linewidth=0.5)
        for bar, auc_val in zip(bars, aucs):
            if not np.isnan(auc_val):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                        f"{auc_val:.3f}", ha="center", va="bottom", fontsize=7)

    ax.axhline(0.50, color="red", linestyle="--", linewidth=1, alpha=0.7, label="AUC=0.50 (random)")
    ax.set_xticks(x + width)
    ax.set_xticklabels([f"{m}\n{z}" for m, z in groups], fontsize=9)
    ax.set_ylabel("AUC (ROC)")
    ax.set_title("Fig01: AUC by Mode x LOH Zone — QS v1 vs v2a vs v2b")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(0.35, 0.85)
    add_caption(fig, CAPTION_SRC)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig01_auc_by_zone.png")


def fig02_auc_global(df: pd.DataFrame, y_bin: np.ndarray) -> Path:
    """Bar chart: overall AUC (not zone-stratified) for v1/v2a/v2b x 2 modes."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(7, 5))

    formulas = ["v1", "v2a", "v2b"]
    qs_cols = {"v1": "QS_v1", "v2a": "QS_v2a", "v2b": "QS_v2b"}
    modes = ["to", "paired"]
    x = np.arange(len(modes))
    width = 0.22

    for i, formula in enumerate(formulas):
        aucs = []
        for mode in modes:
            mask = (df["mode"] == mode) & np.isfinite(y_bin)
            auc_val = _auc_safe(y_bin[mask], df.loc[mask, qs_cols[formula]].values)
            aucs.append(auc_val)
        bars = ax.bar(x + i * width, aucs, width, label=QS_LABELS[formula],
                      color=QS_COLORS[formula], edgecolor="white", linewidth=0.5)
        for bar, auc_val in zip(bars, aucs):
            if not np.isnan(auc_val):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                        f"{auc_val:.3f}", ha="center", va="bottom", fontsize=8)

    ax.axhline(0.50, color="red", linestyle="--", linewidth=1, alpha=0.7, label="AUC=0.50 (random)")
    ax.set_xticks(x + width)
    ax.set_xticklabels([m.upper() for m in modes], fontsize=10)
    ax.set_ylabel("AUC (ROC)")
    ax.set_title("Fig02: Global AUC — QS v1 vs v2a vs v2b")
    ax.legend(fontsize=8)
    ax.set_ylim(0.35, 0.85)
    add_caption(fig, CAPTION_SRC)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig02_auc_global.png")


def fig03_qs_distribution(df: pd.DataFrame) -> Path:
    """Violin plots: QS distribution for each formula x TP/FP x LOH zone (TO only)."""
    setup_plot_style()
    to = df[df["mode"] == "to"].copy()
    formulas = ["QS_v1", "QS_v2a", "QS_v2b"]
    formula_labels = ["QS_v1", "QS_v2a", "QS_v2b"]
    zones = ["inside", "outside"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True)

    for row_i, zone in enumerate(zones):
        sub = to[to["loh_zone"] == zone]
        for col_i, (fcol, flabel) in enumerate(zip(formulas, formula_labels)):
            ax = axes[row_i, col_i]
            plot_data = sub[["truth_label", fcol]].dropna(subset=[fcol])
            if len(plot_data) < 10:
                ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            sns.violinplot(
                data=plot_data, x="truth_label", y=fcol,
                order=["TP", "FP"], palette={"TP": COLOR_TP, "FP": COLOR_FP},
                split=False, inner="quartile", ax=ax, cut=0,
            )
            # Compute per-group medians
            for label in ["TP", "FP"]:
                vals = plot_data.loc[plot_data["truth_label"] == label, fcol]
                if len(vals) > 0:
                    med = vals.median()
                    ax.text(0.02, 0.98, f"n={len(vals):,}", transform=ax.transAxes,
                            fontsize=7, va="top") if label == "TP" else None
            ax.set_title(f"{flabel} | {zone}", fontsize=10)
            ax.set_xlabel("")
            if col_i == 0:
                ax.set_ylabel(f"QS ({zone} LOH zone)")
            else:
                ax.set_ylabel("")

    fig.suptitle("Fig03: QS Distribution by Formula x Truth Label x LOH Zone (TO mode only)",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    add_caption(fig, CAPTION_SRC, y=-0.01)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig03_qs_distribution.png")


def fig04_per_sample_heatmap(df: pd.DataFrame, y_bin: np.ndarray) -> Path:
    """Heatmap: rows=7 samples, cols=6 (v1/v2a/v2b x to/paired), values=AUC."""
    setup_plot_style()
    formulas = ["QS_v1", "QS_v2a", "QS_v2b"]
    col_labels = []
    for f_short in ["v1", "v2a", "v2b"]:
        for m in ["to", "paired"]:
            col_labels.append(f"{f_short}_{m}")

    data = []
    for sample in SAMPLE_ORDER:
        row = []
        for fcol, f_short in zip(formulas, ["v1", "v2a", "v2b"]):
            for mode in ["to", "paired"]:
                mask = (df["sample"] == sample) & (df["mode"] == mode) & np.isfinite(y_bin)
                auc_val = _auc_safe(y_bin[mask], df.loc[mask, fcol].values)
                row.append(auc_val)
        data.append(row)

    heat_df = pd.DataFrame(data, index=SAMPLE_ORDER, columns=col_labels)

    fig, ax = plt.subplots(figsize=(12, 6))
    # Use diverging colormap centered at 0.5
    vmin = max(0.3, np.nanmin(heat_df.values) - 0.02)
    vmax = min(1.0, np.nanmax(heat_df.values) + 0.02)
    center = 0.5
    if vmin < center < vmax:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
    else:
        norm = None
    sns.heatmap(heat_df, annot=True, fmt=".3f", cmap="RdYlGn", ax=ax,
                linewidths=0.5, linecolor="white", norm=norm,
                cbar_kws={"label": "AUC"})
    ax.set_title("Fig04: Per-Sample AUC Heatmap — QS Formulas x Modes", fontsize=12)
    ax.set_ylabel("Sample")
    ax.set_xlabel("QS Formula_Mode")
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    add_caption(fig, CAPTION_SRC, y=-0.02)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig04_per_sample_heatmap.png")


def fig05_roc_curves(df: pd.DataFrame, y_bin: np.ndarray) -> Path:
    """ROC curves for v1/v2a/v2b (TO mode, LOH-inside only, all samples pooled)."""
    from sklearn.metrics import roc_curve, auc as sk_auc
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(7, 7))

    mask = (df["mode"] == "to") & (df["loh_zone"] == "inside") & np.isfinite(y_bin)
    y_true = y_bin[mask]

    for formula, col in [("v1", "QS_v1"), ("v2a", "QS_v2a"), ("v2b", "QS_v2b")]:
        y_score = df.loc[mask, col].values.astype(float)
        valid = np.isfinite(y_score)
        if valid.sum() < 10 or len(np.unique(y_true[valid])) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true[valid], y_score[valid])
        roc_auc = sk_auc(fpr, tpr)
        ax.plot(fpr, tpr, color=QS_COLORS[formula], linewidth=2,
                label=f"{QS_LABELS[formula]} (AUC={roc_auc:.3f})")

    ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Fig05: ROC Curves — TO Mode, LOH-Inside Only")
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    add_caption(fig, CAPTION_SRC)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig05_roc_curves.png")


def fig06_pr_curves(df: pd.DataFrame, y_bin: np.ndarray) -> Path:
    """Precision-Recall curves (TO mode, LOH-inside only, all samples pooled)."""
    from sklearn.metrics import precision_recall_curve, average_precision_score
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(7, 7))

    mask = (df["mode"] == "to") & (df["loh_zone"] == "inside") & np.isfinite(y_bin)
    y_true = y_bin[mask]

    for formula, col in [("v1", "QS_v1"), ("v2a", "QS_v2a"), ("v2b", "QS_v2b")]:
        y_score = df.loc[mask, col].values.astype(float)
        valid = np.isfinite(y_score)
        if valid.sum() < 10 or len(np.unique(y_true[valid])) < 2:
            continue
        precision, recall, _ = precision_recall_curve(y_true[valid], y_score[valid])
        ap = average_precision_score(y_true[valid], y_score[valid])
        ax.plot(recall, precision, color=QS_COLORS[formula], linewidth=2,
                label=f"{QS_LABELS[formula]} (AP={ap:.3f})")

    # Baseline: prevalence
    prevalence = y_true.mean() if len(y_true) > 0 else 0.5
    ax.axhline(prevalence, color="gray", linestyle=":", linewidth=1, alpha=0.7,
               label=f"Baseline (prevalence={prevalence:.2f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Fig06: Precision-Recall Curves — TO Mode, LOH-Inside Only")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(0.0, 1.05)
    add_caption(fig, CAPTION_SRC)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig06_pr_curves.png")


def fig07_qs_scatter(df: pd.DataFrame) -> Path:
    """Scatter: QS_v2a (y) vs QS_v1 (x), colored by truth_label. TO, LOH-inside."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(8, 8))

    mask = (df["mode"] == "to") & (df["loh_zone"] == "inside") & df["truth_label"].isin(["TP", "FP"])
    sub = df.loc[mask].copy()

    if len(sub) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return save_figure(fig, FIG_DIR / f"{PREFIX}_fig07_qs_scatter.png")

    # Add jitter for discrete-ish values
    rng = np.random.RandomState(42)
    jitter_x = rng.uniform(-1.5, 1.5, size=len(sub))
    jitter_y = rng.uniform(-1.5, 1.5, size=len(sub))

    for label, color in [("FP", COLOR_FP), ("TP", COLOR_TP)]:
        m = sub["truth_label"] == label
        ax.scatter(sub.loc[m, "QS_v1"] + jitter_x[m.values],
                   sub.loc[m, "QS_v2a"] + jitter_y[m.values],
                   c=color, alpha=0.15, s=8, label=f"{label} (n={m.sum():,})",
                   edgecolors="none", rasterized=True)

    ax.plot([0, 100], [0, 100], "k--", linewidth=1, alpha=0.5, label="Identity")
    ax.set_xlabel("QS_v1 (current)")
    ax.set_ylabel("QS_v2a (AD-based)")
    ax.set_title("Fig07: QS_v2a vs QS_v1 — TO Mode, LOH-Inside")
    ax.legend(fontsize=9, markerscale=3)
    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)
    add_caption(fig, CAPTION_SRC)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig07_qs_scatter.png")


def fig08_summary_table(auc_table: pd.DataFrame) -> Path:
    """Rendered table as figure showing AUC comparison with GO/NO-GO verdicts."""
    setup_plot_style()

    # Pivot for display
    display_cols = ["mode", "zone", "sample", "AUC_v1", "AUC_v2a", "AUC_v2b",
                    "delta_v2a", "delta_v2b"]
    display = auc_table[display_cols].copy()
    display = display.sort_values(["mode", "zone", "sample"])

    # Format numeric columns
    for col in ["AUC_v1", "AUC_v2a", "AUC_v2b"]:
        display[col] = display[col].apply(lambda x: f"{x:.3f}" if not np.isnan(x) else "NaN")
    for col in ["delta_v2a", "delta_v2b"]:
        display[col] = display[col].apply(
            lambda x: f"{x:+.3f}" if not np.isnan(x) else "NaN")

    # Create figure with table
    n_rows = len(display) + 1  # header
    fig_height = max(8, 0.35 * n_rows + 2)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    ax.axis("off")

    cell_text = display.values.tolist()
    col_labels = display.columns.tolist()
    table = ax.table(cellText=cell_text, colLabels=col_labels,
                     cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.3)

    # Color header
    for j in range(len(col_labels)):
        table[0, j].set_facecolor("#E0E0E0")
        table[0, j].set_text_props(fontweight="bold")

    # Color delta cells
    for i in range(len(cell_text)):
        for j_offset, col_name in [(6, "delta_v2a"), (7, "delta_v2b")]:
            try:
                val = float(cell_text[i][j_offset].replace("+", ""))
                if val > 0.05:
                    table[i + 1, j_offset].set_facecolor("#C8E6C9")
                elif val < -0.02:
                    table[i + 1, j_offset].set_facecolor("#FFCDD2")
            except (ValueError, IndexError):
                pass

    ax.set_title("Fig08: Full AUC Summary Table — QS v1 vs v2a vs v2b\n"
                 "(Green: delta > +0.05, Red: delta < -0.02)", fontsize=11, pad=20)
    fig.tight_layout()
    add_caption(fig, CAPTION_SRC, y=0.01)
    return save_figure(fig, FIG_DIR / f"{PREFIX}_fig08_summary_table.png")


# ══════════════════════════════════════════════════════════════════════════
# Report Generation
# ══════════════════════════════════════════════════════════════════════════

def generate_report(
    df: pd.DataFrame,
    y_bin: np.ndarray,
    auc_table: pd.DataFrame,
    zone_dist: pd.DataFrame,
    v1_validation: Dict[str, Any],
) -> Path:
    """Generate Markdown report with GO/NO-GO evaluation."""
    lines = []

    lines.append("<!--")
    lines.append(f"建立時間: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("目標: 驗證 LOH zone-aware QS v2 重設計是否改善 TO mode 的 AUC")
    lines.append("處理範圍: 7 samples x 2 modes x 3 QS formulas")
    lines.append("關聯檔案:")
    lines.append("  - scripts/analysis/build_qs_v2_zone_aware_validation.py")
    lines.append("  - research/loh_investigation/figures/qs_v2_fig*.png")
    lines.append("  - research/loh_investigation/data/qs_v2_*.tsv")
    lines.append("-->")
    lines.append("")
    lines.append("# QS v2 Zone-Aware Validation Report")
    lines.append("")
    lines.append(f"**Date**: {datetime.now().strftime('%Y-%m-%d')}")
    lines.append(f"**Script**: `scripts/analysis/build_qs_v2_zone_aware_validation.py`")
    lines.append(f"**Data**: `all_region_rows.tsv.gz` ({len(df):,} rows)")
    lines.append(f"**LOH zones**: per-sample LongPhase-TO LOH.bed (binary inside/outside)")
    lines.append("")

    # V1 validation
    lines.append("## QS v1 Reconstruction Validation")
    lines.append("")
    lines.append(f"- Correlation with existing Quality_Score: **{v1_validation['correlation']:.6f}**")
    lines.append(f"- Max absolute error: **{v1_validation['max_abs_error']:.2f}**")
    lines.append(f"- Mean absolute error: **{v1_validation['mean_abs_error']:.4f}**")
    lines.append(f"- Exact match rate: **{v1_validation['exact_match_rate']:.4f}**")
    lines.append("")

    # Zone distribution
    lines.append("## Zone Distribution")
    lines.append("")
    lines.append("```")
    lines.append(zone_dist.to_string())
    lines.append("```")
    lines.append("")

    # Key AUC tables
    lines.append("## AUC Comparison: By Zone")
    lines.append("")
    lines.append("### TO Mode")
    lines.append("")
    to_table = auc_table[auc_table["mode"] == "to"]
    lines.append("| Zone | Sample | AUC_v1 | AUC_v2a | AUC_v2b | delta_v2a | delta_v2b |")
    lines.append("|------|--------|--------|---------|---------|-----------|-----------|")
    for _, row in to_table.sort_values(["zone", "sample"]).iterrows():
        lines.append(
            f"| {row['zone']} | {row['sample']} | "
            f"{row['AUC_v1']:.3f} | {row['AUC_v2a']:.3f} | {row['AUC_v2b']:.3f} | "
            f"{row['delta_v2a']:+.3f} | {row['delta_v2b']:+.3f} |"
        )
    lines.append("")

    lines.append("### Paired Mode")
    lines.append("")
    paired_table = auc_table[auc_table["mode"] == "paired"]
    lines.append("| Zone | Sample | AUC_v1 | AUC_v2a | AUC_v2b | delta_v2a | delta_v2b |")
    lines.append("|------|--------|--------|---------|---------|-----------|-----------|")
    for _, row in paired_table.sort_values(["zone", "sample"]).iterrows():
        lines.append(
            f"| {row['zone']} | {row['sample']} | "
            f"{row['AUC_v1']:.3f} | {row['AUC_v2a']:.3f} | {row['AUC_v2b']:.3f} | "
            f"{row['delta_v2a']:+.3f} | {row['delta_v2b']:+.3f} |"
        )
    lines.append("")

    # Pooled AUC
    lines.append("## Pooled AUC (All Samples)")
    lines.append("")
    lines.append("| Mode | Zone | AUC_v1 | AUC_v2a | AUC_v2b | delta_v2a | delta_v2b |")
    lines.append("|------|------|--------|---------|---------|-----------|-----------|")
    for mode in ["to", "paired"]:
        for zone in ["inside", "outside", "all"]:
            if zone == "all":
                mask = (df["mode"] == mode) & np.isfinite(y_bin)
            else:
                mask = (df["mode"] == mode) & (df["loh_zone"] == zone) & np.isfinite(y_bin)
            auc_v1 = _auc_safe(y_bin[mask], df.loc[mask, "QS_v1"].values)
            auc_v2a = _auc_safe(y_bin[mask], df.loc[mask, "QS_v2a"].values)
            auc_v2b = _auc_safe(y_bin[mask], df.loc[mask, "QS_v2b"].values)
            d_v2a = auc_v2a - auc_v1 if not (np.isnan(auc_v2a) or np.isnan(auc_v1)) else float("nan")
            d_v2b = auc_v2b - auc_v1 if not (np.isnan(auc_v2b) or np.isnan(auc_v1)) else float("nan")
            lines.append(
                f"| {mode} | {zone} | "
                f"{auc_v1:.3f} | {auc_v2a:.3f} | {auc_v2b:.3f} | "
                f"{d_v2a:+.3f} | {d_v2b:+.3f} |"
            )
    lines.append("")

    # GO/NO-GO evaluation
    lines.append("## GO/NO-GO Evaluation")
    lines.append("")

    # Compute thresholds
    # 1. AUC improvement (TO overall) >= +0.05
    to_all_mask = (df["mode"] == "to") & np.isfinite(y_bin)
    auc_v1_to_all = _auc_safe(y_bin[to_all_mask], df.loc[to_all_mask, "QS_v1"].values)
    auc_v2a_to_all = _auc_safe(y_bin[to_all_mask], df.loc[to_all_mask, "QS_v2a"].values)
    auc_v2b_to_all = _auc_safe(y_bin[to_all_mask], df.loc[to_all_mask, "QS_v2b"].values)

    # 2. AUC improvement (TO LOH-inside) >= +0.10
    to_in_mask = (df["mode"] == "to") & (df["loh_zone"] == "inside") & np.isfinite(y_bin)
    auc_v1_to_in = _auc_safe(y_bin[to_in_mask], df.loc[to_in_mask, "QS_v1"].values)
    auc_v2a_to_in = _auc_safe(y_bin[to_in_mask], df.loc[to_in_mask, "QS_v2a"].values)
    auc_v2b_to_in = _auc_safe(y_bin[to_in_mask], df.loc[to_in_mask, "QS_v2b"].values)

    # 3. Cross-sample consistency >= 5/7
    to_inside_table = auc_table[(auc_table["mode"] == "to") & (auc_table["zone"] == "inside")]
    n_improved_v2a = (to_inside_table["delta_v2a"] > 0).sum()
    n_improved_v2b = (to_inside_table["delta_v2b"] > 0).sum()
    n_samples_v2a = to_inside_table["delta_v2a"].notna().sum()
    n_samples_v2b = to_inside_table["delta_v2b"].notna().sum()

    # 4. No regression (LOH-outside) >= -0.02
    to_out_mask = (df["mode"] == "to") & (df["loh_zone"] == "outside") & np.isfinite(y_bin)
    auc_v1_to_out = _auc_safe(y_bin[to_out_mask], df.loc[to_out_mask, "QS_v1"].values)
    auc_v2a_to_out = _auc_safe(y_bin[to_out_mask], df.loc[to_out_mask, "QS_v2a"].values)
    auc_v2b_to_out = _auc_safe(y_bin[to_out_mask], df.loc[to_out_mask, "QS_v2b"].values)

    def _verdict(val, threshold, direction=">="):
        if np.isnan(val):
            return "N/A"
        if direction == ">=":
            return "PASS" if val >= threshold else "FAIL"
        else:
            return "PASS" if val > threshold else "FAIL"

    lines.append("| Criterion | Threshold | v2a Value | v2a Verdict | v2b Value | v2b Verdict |")
    lines.append("|-----------|-----------|-----------|-------------|-----------|-------------|")

    # Criterion 1
    d_v2a_to = auc_v2a_to_all - auc_v1_to_all if not (np.isnan(auc_v2a_to_all) or np.isnan(auc_v1_to_all)) else float("nan")
    d_v2b_to = auc_v2b_to_all - auc_v1_to_all if not (np.isnan(auc_v2b_to_all) or np.isnan(auc_v1_to_all)) else float("nan")
    lines.append(
        f"| TO overall AUC improvement | >= +0.05 | "
        f"{d_v2a_to:+.3f} | {_verdict(d_v2a_to, 0.05)} | "
        f"{d_v2b_to:+.3f} | {_verdict(d_v2b_to, 0.05)} |"
    )

    # Criterion 2
    d_v2a_in = auc_v2a_to_in - auc_v1_to_in if not (np.isnan(auc_v2a_to_in) or np.isnan(auc_v1_to_in)) else float("nan")
    d_v2b_in = auc_v2b_to_in - auc_v1_to_in if not (np.isnan(auc_v2b_to_in) or np.isnan(auc_v1_to_in)) else float("nan")
    lines.append(
        f"| TO LOH-inside AUC improvement | >= +0.10 | "
        f"{d_v2a_in:+.3f} | {_verdict(d_v2a_in, 0.10)} | "
        f"{d_v2b_in:+.3f} | {_verdict(d_v2b_in, 0.10)} |"
    )

    # Criterion 3
    lines.append(
        f"| Cross-sample consistency (TO inside) | >= 5/7 improved | "
        f"{n_improved_v2a}/{n_samples_v2a} | {_verdict(n_improved_v2a, 5)} | "
        f"{n_improved_v2b}/{n_samples_v2b} | {_verdict(n_improved_v2b, 5)} |"
    )

    # Criterion 4
    d_v2a_out = auc_v2a_to_out - auc_v1_to_out if not (np.isnan(auc_v2a_to_out) or np.isnan(auc_v1_to_out)) else float("nan")
    d_v2b_out = auc_v2b_to_out - auc_v1_to_out if not (np.isnan(auc_v2b_to_out) or np.isnan(auc_v1_to_out)) else float("nan")
    lines.append(
        f"| No regression (TO LOH-outside) | >= -0.02 | "
        f"{d_v2a_out:+.3f} | {_verdict(d_v2a_out, -0.02)} | "
        f"{d_v2b_out:+.3f} | {_verdict(d_v2b_out, -0.02)} |"
    )
    lines.append("")

    # Overall verdict
    v2a_pass_count = sum([
        d_v2a_to >= 0.05 if not np.isnan(d_v2a_to) else False,
        d_v2a_in >= 0.10 if not np.isnan(d_v2a_in) else False,
        n_improved_v2a >= 5,
        d_v2a_out >= -0.02 if not np.isnan(d_v2a_out) else False,
    ])
    v2b_pass_count = sum([
        d_v2b_to >= 0.05 if not np.isnan(d_v2b_to) else False,
        d_v2b_in >= 0.10 if not np.isnan(d_v2b_in) else False,
        n_improved_v2b >= 5,
        d_v2b_out >= -0.02 if not np.isnan(d_v2b_out) else False,
    ])

    v2a_overall = "GO" if v2a_pass_count >= 3 else "CONDITIONAL GO" if v2a_pass_count >= 2 else "NO-GO"
    v2b_overall = "GO" if v2b_pass_count >= 3 else "CONDITIONAL GO" if v2b_pass_count >= 2 else "NO-GO"

    lines.append(f"### Overall Verdict")
    lines.append("")
    lines.append(f"- **QS_v2a**: {v2a_pass_count}/4 criteria passed -> **{v2a_overall}**")
    lines.append(f"- **QS_v2b**: {v2b_pass_count}/4 criteria passed -> **{v2b_overall}**")
    lines.append("")

    # Conclusions
    lines.append("## Conclusions")
    lines.append("")
    lines.append(f"1. QS_v1 reconstruction validated (corr={v1_validation['correlation']:.6f}, "
                 f"max_err={v1_validation['max_abs_error']:.2f})")
    lines.append(f"2. TO overall: v1 AUC={auc_v1_to_all:.3f}, "
                 f"v2a AUC={auc_v2a_to_all:.3f} (delta={d_v2a_to:+.3f}), "
                 f"v2b AUC={auc_v2b_to_all:.3f} (delta={d_v2b_to:+.3f})")
    lines.append(f"3. TO LOH-inside: v1 AUC={auc_v1_to_in:.3f}, "
                 f"v2a AUC={auc_v2a_to_in:.3f} (delta={d_v2a_in:+.3f}), "
                 f"v2b AUC={auc_v2b_to_in:.3f} (delta={d_v2b_in:+.3f})")
    lines.append(f"4. Cross-sample consistency: v2a={n_improved_v2a}/{n_samples_v2a}, "
                 f"v2b={n_improved_v2b}/{n_samples_v2b}")
    lines.append(f"5. LOH-outside regression check: v2a={d_v2a_out:+.3f}, v2b={d_v2b_out:+.3f}")
    lines.append("")

    report_path = REPORT_DIR / "20260402_qs_v2_zone_aware_validation.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"[QS-v2] Report written: {report_path}")
    return report_path


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    print("=" * 72)
    print("QS v2 Zone-Aware Validation")
    print("=" * 72)

    # ── 1. Load data ──────────────────────────────────────────────────────
    df = load_master_dataset(usecols=LOAD_COLS)
    print(f"[QS-v2] Loaded {len(df):,} rows")

    # ── 2. Assign LOH zones ──────────────────────────────────────────────
    df = assign_loh_zones_all_samples(df)

    # Zone distribution table
    zone_dist = df.groupby(["mode", "loh_zone", "truth_label"]).size().unstack(fill_value=0)
    zone_dist.to_csv(DATA_DIR / f"{PREFIX}_zone_distribution.tsv", sep="\t")
    print(f"[QS-v2] Zone distribution saved")

    # ── 3. Compute QS formulas ───────────────────────────────────────────
    print("[QS-v2] Computing QS_v1 (mode-aware reconstruction)...")
    df["QS_v1"] = compute_qs_v1(df)

    print("[QS-v2] Computing QS_v2a (zone-aware, AD-based)...")
    df["QS_v2a"] = compute_qs_v2a(df)

    print("[QS-v2] Computing QS_v2b (zone-aware, ISM-only)...")
    df["QS_v2b"] = compute_qs_v2b(df)

    # ── 4. Validate QS_v1 reconstruction ─────────────────────────────────
    print("[QS-v2] Validating QS_v1 reconstruction against existing Quality_Score...")
    existing_qs = df["Quality_Score"].values.astype(float)
    reconstructed_qs = df["QS_v1"].values.astype(float)
    valid_mask = np.isfinite(existing_qs) & np.isfinite(reconstructed_qs)

    corr = np.corrcoef(existing_qs[valid_mask], reconstructed_qs[valid_mask])[0, 1]
    abs_err = np.abs(existing_qs[valid_mask] - reconstructed_qs[valid_mask])
    max_abs_err = abs_err.max()
    mean_abs_err = abs_err.mean()
    exact_match_rate = (abs_err < 0.5).mean()  # within rounding tolerance

    v1_validation = {
        "correlation": corr,
        "max_abs_error": max_abs_err,
        "mean_abs_error": mean_abs_err,
        "exact_match_rate": exact_match_rate,
    }

    print(f"  Correlation:       {corr:.6f}")
    print(f"  Max abs error:     {max_abs_err:.2f}")
    print(f"  Mean abs error:    {mean_abs_err:.4f}")
    print(f"  Exact match rate:  {exact_match_rate:.4f}")

    if max_abs_err > 1.0:
        print("  [WARN] Max abs error > 1.0 — investigating mismatches...")
        mismatch_mask = valid_mask & (np.abs(existing_qs - reconstructed_qs) > 0.5)
        if mismatch_mask.sum() > 0:
            mismatch_sample = df.loc[mismatch_mask, ["sample", "mode", "Quality_Score", "QS_v1",
                                                      "Potential_LOH", "HPMergedSig", "AlleleSig"]].head(10)
            print(mismatch_sample.to_string())

    # ── 5. Encode truth labels ───────────────────────────────────────────
    y_bin = encode_truth_binary(df).values.astype(float)
    n_tp = np.nansum(y_bin == 1)
    n_fp = np.nansum(y_bin == 0)
    n_other = np.isnan(y_bin).sum()
    print(f"[QS-v2] Truth labels: TP={int(n_tp):,}, FP={int(n_fp):,}, other={int(n_other):,}")

    # ── 6. Build per-sample AUC table ────────────────────────────────────
    print("[QS-v2] Computing per-sample AUC table...")
    auc_rows = []
    for mode in ["to", "paired"]:
        for zone in ["inside", "outside"]:
            for sample in SAMPLE_ORDER:
                mask = ((df["mode"] == mode) & (df["loh_zone"] == zone) &
                        (df["sample"] == sample) & np.isfinite(y_bin))
                n = mask.sum()
                auc_v1 = _auc_safe(y_bin[mask], df.loc[mask, "QS_v1"].values) if n >= 10 else float("nan")
                auc_v2a = _auc_safe(y_bin[mask], df.loc[mask, "QS_v2a"].values) if n >= 10 else float("nan")
                auc_v2b = _auc_safe(y_bin[mask], df.loc[mask, "QS_v2b"].values) if n >= 10 else float("nan")
                delta_v2a = auc_v2a - auc_v1 if not (np.isnan(auc_v2a) or np.isnan(auc_v1)) else float("nan")
                delta_v2b = auc_v2b - auc_v1 if not (np.isnan(auc_v2b) or np.isnan(auc_v1)) else float("nan")
                auc_rows.append({
                    "mode": mode, "zone": zone, "sample": sample, "n": n,
                    "AUC_v1": auc_v1, "AUC_v2a": auc_v2a, "AUC_v2b": auc_v2b,
                    "delta_v2a": delta_v2a, "delta_v2b": delta_v2b,
                })

    auc_table = pd.DataFrame(auc_rows)
    auc_table.to_csv(DATA_DIR / f"{PREFIX}_auc_table.tsv", sep="\t", index=False)
    print(f"[QS-v2] AUC table saved ({len(auc_table)} rows)")

    # Print summary
    print("\n[QS-v2] === Pooled AUC Summary ===")
    for mode in ["to", "paired"]:
        for zone in ["inside", "outside", "all"]:
            if zone == "all":
                mask = (df["mode"] == mode) & np.isfinite(y_bin)
                label = f"{mode:7s} | all"
            else:
                mask = (df["mode"] == mode) & (df["loh_zone"] == zone) & np.isfinite(y_bin)
                label = f"{mode:7s} | {zone:7s}"
            auc_v1 = _auc_safe(y_bin[mask], df.loc[mask, "QS_v1"].values)
            auc_v2a = _auc_safe(y_bin[mask], df.loc[mask, "QS_v2a"].values)
            auc_v2b = _auc_safe(y_bin[mask], df.loc[mask, "QS_v2b"].values)
            print(f"  {label}: v1={auc_v1:.3f}  v2a={auc_v2a:.3f}  v2b={auc_v2b:.3f}")

    # ── 7. Save QS columns for further analysis ─────────────────────────
    qs_export = df[["RegionID", "sample", "mode", "truth_label", "loh_zone",
                    "Quality_Score", "QS_v1", "QS_v2a", "QS_v2b"]].copy()
    qs_export.to_csv(DATA_DIR / f"{PREFIX}_qs_scores.tsv.gz", sep="\t", index=False,
                     compression="gzip")
    print(f"[QS-v2] QS scores exported ({len(qs_export):,} rows)")

    # ── 8. Generate figures ──────────────────────────────────────────────
    print("\n[QS-v2] Generating figures...")
    fig01_auc_by_zone(df, y_bin)
    fig02_auc_global(df, y_bin)
    fig03_qs_distribution(df)
    fig04_per_sample_heatmap(df, y_bin)
    fig05_roc_curves(df, y_bin)
    fig06_pr_curves(df, y_bin)
    fig07_qs_scatter(df)
    fig08_summary_table(auc_table)

    # ── 9. Generate report ───────────────────────────────────────────────
    print("\n[QS-v2] Generating report...")
    generate_report(df, y_bin, auc_table, zone_dist, v1_validation)

    elapsed = time.time() - t_start
    print(f"\n[QS-v2] Done in {elapsed:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
