"""
Shared I/O utilities for NG × KDE rescaling analysis.

Responsibilities:
  - Locate NEW KDE master TSVs (6 samples post-fix + HCC1954 stale rescale)
  - Join per-sample rescaling ratios (D3: quantile_drift_per_sample.tsv)
  - Consolidate TP + FP rows with sample / mode / flag columns

Run order: step0 uses this module to produce a single long-form dataset
for subsequent stratification (step1) and visualization (step3).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Canonical sample ordering (for figures / reports)
SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
    "H2009",
    "H1437",
    "COLO829",
]

# Color palette (Claude Code PPTX standard)
PALETTE = {
    "dark": "#1E2A44",
    "bg": "#F7F3EC",
    "accent": "#A85540",
    "green": "#009E73",
    "red": "#D55E00",
    "blue": "#5B8DB7",
    "grey": "#8A8A8A",
}

# ----------------------------------------------------------------------------
# NEW KDE paired_full master paths (verified 2026-04-22 to contain
# Diploid_Coverage_Used column, confirming post-fix binary)
NEW_KDE_PAIRED_FULL = {
    "HCC1395": "/big7_disk/liaoyoyo2001/InterSubMod/output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/InterSubMod/output/canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full",
    "H2009": "/big7_disk/liaoyoyo2001/InterSubMod/output/canonical/H2009/paired_full/20260421_H2009_paired_full_full",
    "H1437": "/big7_disk/liaoyoyo2001/InterSubMod/output/canonical/H1437/paired_full/20260421_H1437_paired_full_full",
    "COLO829": "/big7_disk/liaoyoyo2001/InterSubMod/output/canonical/COLO829/paired_full/20260421_COLO829_paired_full_full",
    "HCC1937": "/big7_disk/liaoyoyo2001/InterSubMod/output/canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full",
}

# HCC1954 fallback: stale 20260315 + post-hoc rescale (ratio 1.2295)
HCC1954_STALE_PAIRED_FULL = "/big7_disk/liaoyoyo2001/InterSubMod/output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix"

# HCC1395 TO (post-fix) sourced from Phase 1 HP-only workspace
HCC1395_TO_TP = "/tmp/ism_hp_fix_phase1/tp_off/significance_summary.csv"
HCC1395_TO_FP = "/tmp/ism_hp_fix_phase1/fp_off/significance_summary.csv"

# 2026-04-22: 5 additional TO samples reran from archive pilots
# (see research/tpfp_loh_af_kde_discrimination/scripts/obs15_rerun_archive_to_ism.sh)
# Each points to the pilot's step05_intersubmod directory containing
# intersubmod_tp/significance_summary.csv and intersubmod_fp/significance_summary.csv
ARCHIVE_TO_RERUN_ROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots"
NEW_KDE_TO_ARCHIVE = {
    "HCC1395_DORADO": f"{ARCHIVE_TO_RERUN_ROOT}/20260315_hcc1395_dorado_to_pilot/step05_intersubmod",
    "H1437":          f"{ARCHIVE_TO_RERUN_ROOT}/20260318_h1437_to_pilot_fastresume/step05_intersubmod",
    "H2009":          f"{ARCHIVE_TO_RERUN_ROOT}/20260318_h2009_to_pilot_fastresume/step05_intersubmod",
    "HCC1937":        f"{ARCHIVE_TO_RERUN_ROOT}/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod",
    "HCC1954":        f"{ARCHIVE_TO_RERUN_ROOT}/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod",
}

# Per-sample baselines (from research/kde_fix_validation/outputs/step3_quantile_drift/)
QUANTILE_DRIFT_TSV = "/big7_disk/liaoyoyo2001/InterSubMod/research/kde_fix_validation/outputs/step3_quantile_drift/quantile_drift_per_sample.tsv"

# SEQC2 CN empirical calibration (HCC1395 paired_pileup per_cn_bin_metrics)
SEQC2_PER_CN_BIN_TSV = "/big7_disk/liaoyoyo2001/InterSubMod/research/kde_fix_validation/outputs/step1_hcc1395_seqc2/paired_pileup/per_cn_bin_metrics.tsv"

# Output workspace
WORKSPACE = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling")
DATA_DIR = WORKSPACE / "data"
OBS_DIR = WORKSPACE / "observations"
FIG_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled"
)

# ----------------------------------------------------------------------------
# AF bin definitions (user request 2026-04-22: 10-bin 0.1 granularity)
AF_BIN_EDGES_FINE = np.arange(0.0, 1.01, 0.1)
AF_BIN_LABELS_FINE = [f"[{AF_BIN_EDGES_FINE[i]:.1f},{AF_BIN_EDGES_FINE[i+1]:.1f})" for i in range(len(AF_BIN_EDGES_FINE) - 1)]


def af_class_coarse(af: float) -> str:
    """20260414-compatible coarse class."""
    if pd.isna(af):
        return "Unknown"
    if af < 0.1 or af > 0.9:
        return "Extreme"
    if 0.4 <= af <= 0.6:
        return "Near-half"
    return "Intermediate"


# ----------------------------------------------------------------------------
# CN tier cut strategies
CN_STRATEGIES: Dict[str, List[float]] = {
    # A. Legacy 20260414
    "A_Legacy": [0.75, 1.25, 1.75],
    # B. Clinical Strict (FACETS/Battenberg)
    "B_ClinicalStrict": [0.5, 0.85, 1.15, 1.4],
    # C. 6-tier Fine (TITAN-inspired)
    "C_FineGrained": [0.4, 0.7, 0.95, 1.1, 1.4, 2.0],
    # F. SEQC2-grounded (HCC1395 paired_pileup empirical midpoints)
    # boundaries from per_cn_bin_metrics: CN=1-2 mean 0.47 / CN=2 0.83 / CN=3 1.14 / CN=4 1.52 / CN≥5 2.11
    # midpoints: 0.65 (1-2 vs 2), 0.99 (2 vs 3), 1.33 (3 vs 4), 1.82 (4 vs 5+)
    "F_SEQC2_grounded": [0.65, 0.99, 1.33, 1.82],
}


def assign_cn_tier(series: pd.Series, edges: List[float]) -> pd.Series:
    """Bucket CovM values into ordered tiers given boundary list."""
    labels = [f"T{i}" for i in range(len(edges) + 1)]
    return pd.cut(series, bins=[-np.inf] + edges + [np.inf], labels=labels, include_lowest=True)


def quantile_boundaries(series: pd.Series, qs: Tuple[float, ...] = (0.05, 0.25, 0.5, 0.75, 0.95)) -> List[float]:
    """Strategy D: per-sample empirical quantiles."""
    return [float(series.quantile(q)) for q in qs]


# ----------------------------------------------------------------------------
# Ingestion helpers
REQUIRED_COLS = [
    "RegionID", "Chr", "Pos", "NumReads", "NumCpGs",
    "HPFineNGroups",
    "Coverage_Multiple", "Coverage_Category",
    "AlleleDelta",
]

# Optional (availability varies by run): keep when present
OPTIONAL_COLS = [
    "HPFine_NGroups_CF",
    "Diploid_Coverage_Used",
    "Potential_LOH", "NTumorReads", "NNormalReads",
    "LOH_Bed_Overlap", "LOH_Bed_Annotation", "LOH_Subtype",
    "label",
]


def _read_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    keep = [c for c in REQUIRED_COLS + OPTIONAL_COLS if c in df.columns]
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing required columns {missing}")
    return df[keep].copy()


def load_sample_run(tp_path: Path, fp_path: Path, sample: str, mode: str, kde_status: str) -> pd.DataFrame:
    """Load TP + FP for a single sample/mode run."""
    tp = _read_summary(tp_path)
    tp["tp_label"] = 1
    fp = _read_summary(fp_path)
    fp["tp_label"] = 0
    df = pd.concat([tp, fp], ignore_index=True)
    df["sample"] = sample
    df["mode"] = mode
    df["kde_status"] = kde_status  # "new_direct" | "stale_rescaled" | "phase1_new"
    return df


def load_baselines() -> pd.DataFrame:
    df = pd.read_csv(QUANTILE_DRIFT_TSV, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df


def get_sample_ratio(sample: str, baselines: pd.DataFrame, mode: str = "paired_full") -> float:
    sub = baselines[(baselines["sample"] == sample) & (baselines["fixed_mode"] == mode)]
    if sub.empty:
        # fallback to paired_pileup
        sub = baselines[(baselines["sample"] == sample) & (baselines["fixed_mode"] == "paired_pileup")]
    if sub.empty:
        raise KeyError(f"No baseline for {sample}/{mode}")
    return float(sub["ratio_stale_over_new"].iloc[0])


def get_sample_baseline(sample: str, baselines: pd.DataFrame, mode: str = "paired_full") -> float:
    sub = baselines[(baselines["sample"] == sample) & (baselines["fixed_mode"] == mode)]
    if sub.empty:
        sub = baselines[(baselines["sample"] == sample) & (baselines["fixed_mode"] == "paired_pileup")]
    return float(sub["new_baseline_x"].iloc[0])


# ----------------------------------------------------------------------------
# LOH classification helper
def loh_side(row) -> str:
    """Classify region as LOH Inner / Outer / Unknown.

    Preferred source: LOH_Bed_Overlap (external PON LOH.bed, Jaccard=1.0 vs self-phased per 2026-04-19 audit).
    Fallback: Potential_LOH (self-phased LOH detection from ReadParser).

    In the 2026-04-20/21 new KDE runs, LOH.bed was NOT loaded (LOH_Bed_Overlap all False),
    so we fall back to Potential_LOH which matches what the 2026-04-14 POSITIVE analysis used.
    """
    bed = row.get("LOH_Bed_Overlap")
    if pd.notna(bed):
        try:
            if int(bed) == 1 or str(bed).lower() == "true":
                return "Inner"
        except (TypeError, ValueError):
            pass
        if str(bed).lower() == "true":
            return "Inner"

    pot = row.get("Potential_LOH")
    if pd.notna(pot):
        if str(pot).lower() == "true" or (isinstance(pot, (int, float)) and int(pot) == 1):
            return "Inner"
        if str(pot).lower() == "false" or (isinstance(pot, (int, float)) and int(pot) == 0):
            return "Outer"

    return "Unknown"
