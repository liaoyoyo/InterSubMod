#!/usr/bin/env python3
"""
Phase 2 Cycle 2 Step B2' — Build per-sample augmented master TSVs for 4 new samples
{H1437, H2009, HCC1954, HCC1937}.

Inputs per sample:
  - Step4 per-sample master TSV: research/v6_bam_tpfp_hp_loh_cn/step4_cross_sample_extension/
      per_sample_master/{sample}_v6_master.tsv  (region-level chassis with region_id, chr,
      pos, label, caller_af, Coverage_Multiple, loh_side, V6_off_NG, NumReads_master, ...)
  - phaseC V6 significance_summary.csv (V6_off_tp + V6_off_fp) under
      research/paired_priority_bug_audit/phaseC_v6_4sample_with_significance/{sample}/

Output:
  - cycle2/data/{sample}_master_augmented.tsv  (per-sample, V6_off only)

Schema (cycle1 filter compatible):
  region_id, chr, pos, label, y,
  caller_af, NumReads_master,
  loh_inner_flag, Coverage_Multiple, Coverage_Multiple_imp, chr8_flag,
  V6_off_NG,
  V6_off_meth_HPMergedDelta,
  V6_off_meth_HPFineF,
  V6_off_meth_NME_imbalance        (derived: |NME_HP1 - NME_HP2|)
  V6_off_meth_Epipoly_Delta,
  V6_off_meth_ClusterPermanovaF

Plus the rest of the methylation features for richer context:
  V6_off_meth_{HPMergedP, HPMergedSig, HPFineP, HPFineSig, HPFineNGroups,
               ClusterPermanovaP, AlleleDelta, AlleleP, Entropy_Imbalance,
               NME_HP1, NME_HP2, Epipoly_HP1, Epipoly_HP2, Epipoly_imbalance}

Sanity:
  - Methylation join key (chr, pos) match rate per (sample, label)
  - 13 methylation feature non-null rate per sample
"""
from __future__ import annotations

import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
STEP4_DIR = (
    REPO_ROOT
    / "research/v6_bam_tpfp_hp_loh_cn/step4_cross_sample_extension/per_sample_master"
)
SIG_ROOT = (
    REPO_ROOT
    / "research/paired_priority_bug_audit/phaseC_v6_4sample_with_significance"
)
OUT_DIR = (
    REPO_ROOT
    / "research/methyl_augmented_filter_phase2/cycle2"
)
OUT_DATA_DIR = OUT_DIR / "data"
OUT_LOG = OUT_DIR / "intermediate" / "build_4sample_master.log"
SANITY_PATH = OUT_DIR / "intermediate" / "step_b2_prime_sanity.md"

SAMPLES = ["H1437", "H2009", "HCC1954", "HCC1937"]
LABELS = ["tp", "fp"]

# 13 methylation feature names that we want to import + 2 derived imbalance.
# Match cycle 1 v1.0 build_augmented_master.py METHYL_FEATURES.
METHYL_FEATURES = [
    "HPMergedDelta",
    "HPMergedP",
    "HPMergedSig",
    "HPFineF",
    "HPFineP",
    "HPFineSig",
    "HPFineNGroups",
    "ClusterPermanovaF",
    "ClusterPermanovaP",
    "AlleleDelta",
    "AlleleP",
    "Entropy_Imbalance",
    "NME_HP1",
    "NME_HP2",
    "Epipoly_HP1",
    "Epipoly_HP2",
    "Epipoly_Delta",
]


def setup_logger() -> logging.Logger:
    OUT_LOG.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("build_4sample_master")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    fh = logging.FileHandler(OUT_LOG, mode="w")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.handlers.clear()
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def load_sig(sample: str, label: str, logger: logging.Logger) -> pd.DataFrame:
    """Load V6_off significance_summary.csv for one (sample, label)."""
    path = SIG_ROOT / sample / f"V6_off_{label}" / "significance_summary.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    missing = [c for c in METHYL_FEATURES if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    cols = ["Chr", "Pos"] + METHYL_FEATURES
    sub = df[cols].copy()
    sub.rename(columns={"Chr": "chr", "Pos": "pos"}, inplace=True)
    logger.info(f"    sig V6_off_{label}: rows={len(sub):,}")
    return sub


def build_methyl_combo(sample: str, logger: logging.Logger) -> pd.DataFrame:
    """Concat TP + FP significance rows, derive imbalance features, rename to V6_off_meth_*."""
    parts = [load_sig(sample, label, logger) for label in LABELS]
    combo = pd.concat(parts, ignore_index=True)

    # Derived imbalance features (mirrors cycle 1 build_augmented_master.py).
    if {"NME_HP1", "NME_HP2"}.issubset(combo.columns):
        combo["NME_imbalance"] = (combo["NME_HP1"] - combo["NME_HP2"]).abs()
    if {"Epipoly_HP1", "Epipoly_HP2"}.issubset(combo.columns):
        combo["Epipoly_imbalance"] = (combo["Epipoly_HP1"] - combo["Epipoly_HP2"]).abs()

    feature_cols = METHYL_FEATURES + ["NME_imbalance", "Epipoly_imbalance"]
    rename_map = {f: f"V6_off_meth_{f}" for f in feature_cols}
    combo.rename(columns=rename_map, inplace=True)

    dup_count = combo.duplicated(subset=["chr", "pos"]).sum()
    if dup_count > 0:
        logger.warning(f"    {sample}: {dup_count} duplicate (chr,pos) rows; dropping later.")
        combo = combo.drop_duplicates(subset=["chr", "pos"], keep="first")
    logger.info(
        f"    combo {sample}: rows={len(combo):,}, methyl cols={len(feature_cols)}"
    )
    return combo


def add_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add cycle1-canonical derived columns: loh_inner_flag, Coverage_Multiple_imp, chr8_flag, y."""
    df = df.copy()
    df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    cov_median = df["Coverage_Multiple"].median()
    df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(cov_median)
    df["chr8_flag"] = (df["chr"] == "chr8").astype(int)
    df["y"] = (df["label"] == "TP").astype(int)
    return df


def main() -> None:
    logger = setup_logger()
    t0 = time.time()
    logger.info("=" * 78)
    logger.info("Step B2' — Build 4-sample per-sample augmented master TSVs (V6_off only)")
    logger.info("=" * 78)

    OUT_DATA_DIR.mkdir(parents=True, exist_ok=True)

    sanity_records = []  # list of dicts
    methyl_feature_cols = [
        f"V6_off_meth_{f}" for f in METHYL_FEATURES + ["NME_imbalance", "Epipoly_imbalance"]
    ]

    for sample in SAMPLES:
        logger.info(f"\n=== {sample} ===")

        # 1) Load step4 master.
        master_path = STEP4_DIR / f"{sample}_v6_master.tsv"
        master = pd.read_csv(master_path, sep="\t", low_memory=False)
        logger.info(
            f"  step4 master: rows={len(master):,} cols={len(master.columns)}"
            f" labels={master['label'].value_counts().to_dict()}"
        )

        # 2) Build per-sample methylation combo.
        combo = build_methyl_combo(sample, logger)

        # 3) Merge on (chr, pos).
        before_cols = len(master.columns)
        merged = master.merge(combo, on=["chr", "pos"], how="left")
        added = len(merged.columns) - before_cols
        logger.info(
            f"  merged: rows={len(merged):,} cols={len(merged.columns)} "
            f"(+{added} methyl cols)"
        )

        # 4) Derived columns.
        merged = add_derived_columns(merged)

        # 5) Write per-sample TSV.
        out_path = OUT_DATA_DIR / f"{sample}_master_augmented.tsv"
        merged.to_csv(out_path, sep="\t", index=False)
        size_mb = out_path.stat().st_size / (1024 * 1024)
        logger.info(f"  -> {out_path} ({size_mb:.1f} MB)")

        # 6) Sanity: methyl non-null rate per feature.
        nonnull_summary = {}
        for col in methyl_feature_cols:
            if col in merged.columns:
                nn_pct = merged[col].notna().mean() * 100
                nonnull_summary[col] = nn_pct
        min_nonnull = min(nonnull_summary.values()) if nonnull_summary else 0.0
        med_nonnull = float(np.median(list(nonnull_summary.values()))) if nonnull_summary else 0.0
        logger.info(
            f"  methyl non-null: min={min_nonnull:.2f}% median={med_nonnull:.2f}% "
            f"(across {len(methyl_feature_cols)} cols)"
        )

        # 7) Sanity: cycle1 essential feature non-null check.
        cycle1_essential = [
            "V6_off_NG", "caller_af", "loh_inner_flag", "Coverage_Multiple_imp",
            "chr8_flag", "NumReads_master",
            "V6_off_meth_HPMergedDelta", "V6_off_meth_HPFineF",
            "V6_off_meth_NME_imbalance", "V6_off_meth_Epipoly_Delta",
            "V6_off_meth_ClusterPermanovaF",
        ]
        essential_nn = {
            c: float(merged[c].notna().mean() * 100)
            for c in cycle1_essential
            if c in merged.columns
        }
        logger.info(f"  cycle1 essential non-null %: {essential_nn}")

        labels = merged["label"].value_counts().to_dict()
        sanity_records.append({
            "sample": sample,
            "rows": int(len(merged)),
            "cols": int(len(merged.columns)),
            "n_TP": int(labels.get("TP", 0)),
            "n_FP": int(labels.get("FP", 0)),
            "methyl_nonnull_min_pct": float(min_nonnull),
            "methyl_nonnull_median_pct": float(med_nonnull),
            "essential_min_pct": float(min(essential_nn.values())) if essential_nn else 0.0,
            "methyl_per_feature_nonnull": nonnull_summary,
        })

    # ---------- Sanity Markdown ----------
    sanity_lines = [
        "# Step B2' — 4-sample augmented master TSV sanity\n\n",
        f"> Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}  \n",
        f"> Script: `cycle2/scripts/build_4sample_master.py`  \n\n",
        "## 1. Per-sample row counts and label distribution\n\n",
        "| Sample | Rows | Cols | n_TP | n_FP | TP/FP ratio |\n|---|---:|---:|---:|---:|---:|\n",
    ]
    for r in sanity_records:
        ratio = (r["n_TP"] / r["n_FP"]) if r["n_FP"] else float("nan")
        sanity_lines.append(
            f"| {r['sample']} | {r['rows']:,} | {r['cols']} | {r['n_TP']:,} | "
            f"{r['n_FP']:,} | {ratio:.2f} |\n"
        )

    sanity_lines.append("\n## 2. Methylation non-null rate per sample\n\n")
    sanity_lines.append("| Sample | min non-null % | median non-null % | essential cycle1 min % |\n|---|---:|---:|---:|\n")
    for r in sanity_records:
        sanity_lines.append(
            f"| {r['sample']} | {r['methyl_nonnull_min_pct']:.2f}% | "
            f"{r['methyl_nonnull_median_pct']:.2f}% | "
            f"{r['essential_min_pct']:.2f}% |\n"
        )

    sanity_lines.append("\n## 3. Per-feature non-null rate (13 + 2 derived)\n\n")
    sanity_lines.append(
        "| Feature | "
        + " | ".join(SAMPLES)
        + " |\n|---|"
        + "|".join(["---:"] * len(SAMPLES))
        + "|\n"
    )
    for col in methyl_feature_cols:
        cells = []
        for s in SAMPLES:
            for r in sanity_records:
                if r["sample"] == s:
                    pct = r["methyl_per_feature_nonnull"].get(col, np.nan)
                    cells.append(f"{pct:.2f}%" if not np.isnan(pct) else "n/a")
                    break
        sanity_lines.append(f"| `{col}` | " + " | ".join(cells) + " |\n")

    sanity_lines.append("\n## 4. Verdict\n\n")
    min_essential = min(r["essential_min_pct"] for r in sanity_records)
    min_methyl = min(r["methyl_nonnull_min_pct"] for r in sanity_records)
    sanity_lines.append(
        f"- Cycle1 essential features min non-null across samples: **{min_essential:.2f}%** "
        f"(target ≥ 95% for direct apply)\n"
    )
    sanity_lines.append(
        f"- 13+2 methylation features min non-null across samples: **{min_methyl:.2f}%** "
        f"(NaN → cycle1 median impute downstream)\n"
    )

    SANITY_PATH.write_text("".join(sanity_lines))
    logger.info(f"\nWrote sanity md: {SANITY_PATH}")

    elapsed = time.time() - t0
    logger.info("=" * 78)
    logger.info(f"Step B2' COMPLETE in {elapsed:.1f} s ({elapsed/60:.2f} min)")
    logger.info("=" * 78)


if __name__ == "__main__":
    main()
