#!/usr/bin/env python3
"""
Step 0 — Build augmented master TSV for Methylation-Augmented FP Filter Pilot v1.0.

Join 13 methylation features (plus 2 derived imbalance features) from phaseC_genome_three_way_with_significance/
{V3F,V5,V6}_{off,on}_{tp,fp}/significance_summary.csv into the v0.3 master TSV.

Keys:
  - Master TSV: (chr, pos) where pos = midpoint of region_id chr_<start>_<end> (10kb window).
  - significance_summary.csv: (Chr, Pos).

Output naming:
  {BAM}_{flag}_meth_{feature_name}  e.g. V6_off_meth_HPMergedDelta

Derived features:
  {BAM}_{flag}_meth_NME_imbalance     = abs(NME_HP1 - NME_HP2)
  {BAM}_{flag}_meth_Epipoly_imbalance = abs(Epipoly_HP1 - Epipoly_HP2)
"""
from __future__ import annotations

import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER_PATH = (
    REPO_ROOT
    / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_three_way.tsv"
)
SIG_ROOT = (
    REPO_ROOT
    / "research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance"
)
OUT_DIR = REPO_ROOT / "research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot"
OUT_PATH = OUT_DIR / "step5_master_augmented.tsv"
LOG_PATH = OUT_DIR / "intermediate/build_log.txt"
SANITY_PATH = OUT_DIR / "intermediate/step0_sanity_check.md"
SCHEMA_PATH = OUT_DIR / "step0_schema.md"

BAMS = ["V3F", "V5", "V6"]
FLAGS = ["off", "on"]
LABELS = ["tp", "fp"]

# 13 methylation features named in plan + 2 derived.
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


def safe_corr(a: pd.Series, b: pd.Series) -> float:
    """Compute Pearson r robustly; return nan if dtype/variance/length issues."""
    try:
        x = pd.to_numeric(a, errors="coerce")
        y = pd.to_numeric(b, errors="coerce")
        valid = pd.concat([x, y], axis=1).dropna()
        if len(valid) < 3:
            return float("nan")
        xv = valid.iloc[:, 0].to_numpy(dtype=float)
        yv = valid.iloc[:, 1].to_numpy(dtype=float)
        if xv.std() < 1e-12 or yv.std() < 1e-12:
            return float("nan")
        return float(np.corrcoef(xv, yv)[0, 1])
    except (ValueError, TypeError, AttributeError):
        return float("nan")


def setup_logger() -> logging.Logger:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("build_augmented_master")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(LOG_PATH, mode="w")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.handlers.clear()
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def load_sig(bam: str, flag: str, label: str, logger: logging.Logger) -> pd.DataFrame:
    """Load one significance_summary.csv with required columns."""
    path = SIG_ROOT / f"{bam}_{flag}_{label}" / "significance_summary.csv"
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)
    missing = [c for c in METHYL_FEATURES if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    cols = ["Chr", "Pos"] + METHYL_FEATURES
    sub = df[cols].copy()
    sub.rename(columns={"Chr": "chr", "Pos": "pos"}, inplace=True)
    logger.info(f"  loaded {bam}_{flag}_{label}: rows={len(sub):,}")
    return sub


def build_combo_frame(bam: str, flag: str, logger: logging.Logger) -> pd.DataFrame:
    """Concat TP + FP for a (BAM, flag) combo and rename methylation cols."""
    parts = [load_sig(bam, flag, label, logger) for label in LABELS]
    combo = pd.concat(parts, ignore_index=True)

    # Derived imbalance features.
    if {"NME_HP1", "NME_HP2"}.issubset(combo.columns):
        combo["NME_imbalance"] = (combo["NME_HP1"] - combo["NME_HP2"]).abs()
    if {"Epipoly_HP1", "Epipoly_HP2"}.issubset(combo.columns):
        combo["Epipoly_imbalance"] = (combo["Epipoly_HP1"] - combo["Epipoly_HP2"]).abs()

    feature_cols = METHYL_FEATURES + ["NME_imbalance", "Epipoly_imbalance"]
    rename_map = {f: f"{bam}_{flag}_meth_{f}" for f in feature_cols}
    combo.rename(columns=rename_map, inplace=True)

    # Detect duplicate (chr, pos) across TP+FP partitions (should not happen but guard).
    dup_count = combo.duplicated(subset=["chr", "pos"]).sum()
    if dup_count > 0:
        logger.warning(
            f"  {bam}_{flag}: {dup_count} duplicate (chr,pos) rows; dropping later."
        )
        combo = combo.drop_duplicates(subset=["chr", "pos"], keep="first")
    logger.info(f"  combo {bam}_{flag}: rows={len(combo):,}, methyl cols={len(feature_cols)}")
    return combo


def main() -> None:
    logger = setup_logger()
    t0 = time.time()
    logger.info("=" * 70)
    logger.info("Step 0 — Build augmented master TSV (v0.3 + phaseC methylation)")
    logger.info("=" * 70)

    # ------------------------------------------------------------------ Stage 1
    logger.info("[Stage 1] Schema confirmation")
    sample_sig = load_sig("V6", "off", "tp", logger)
    logger.info(f"  V6_off_tp sample sig columns confirmed: {list(sample_sig.columns)[:8]}...")
    expected_tp = 30475
    actual_tp = len(sample_sig)
    logger.info(
        f"  V6_off_tp row count: expected={expected_tp}, actual={actual_tp}, "
        f"diff={actual_tp - expected_tp}"
    )

    # ------------------------------------------------------------------ Stage 2
    logger.info("[Stage 2] Key alignment")
    master = pd.read_csv(MASTER_PATH, sep="\t")
    logger.info(f"  Master TSV: rows={len(master):,}, cols={len(master.columns)}")
    logger.info(f"  Master label dist: {master['label'].value_counts().to_dict()}")

    # Build (chr, pos) -> midpoint check from region_id.
    parts = master["region_id"].str.split("_", expand=True)
    derived_mid = parts[1].astype(int) + 5000
    pos_match_rate = (master["pos"] == derived_mid).mean()
    logger.info(f"  region_id midpoint == pos: {pos_match_rate*100:.4f}%")
    assert pos_match_rate == 1.0, "region_id midpoint must equal pos column"

    # ------------------------------------------------------------------ Stage 3
    logger.info("[Stage 3] Build augmented master TSV")
    augmented = master.copy()
    base_cols = len(augmented.columns)
    logger.info(f"  Base cols: {base_cols}")

    pearson_records = []  # for sanity comparisons later
    feature_cols_added = []

    for bam in BAMS:
        for flag in FLAGS:
            logger.info(f"  -> Processing {bam}_{flag}")
            combo = build_combo_frame(bam, flag, logger)
            # Left-join on (chr, pos).
            before_cols = len(augmented.columns)
            augmented = augmented.merge(combo, on=["chr", "pos"], how="left")
            added = len(augmented.columns) - before_cols
            new_cols = [c for c in augmented.columns[-added:]]
            feature_cols_added.extend(new_cols)
            logger.info(f"     merged: +{added} cols (total={len(augmented.columns)})")

    logger.info(f"  Final augmented shape: rows={len(augmented):,}, cols={len(augmented.columns)}")
    logger.info(f"  Methylation cols added: {len(feature_cols_added)}")

    # ------------------------------------------------------------------ Save TSV
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(OUT_PATH, sep="\t", index=False)
    out_size_mb = OUT_PATH.stat().st_size / (1024 * 1024)
    logger.info(f"  Wrote {OUT_PATH} ({out_size_mb:.1f} MB)")

    # ------------------------------------------------------------------ Stage 4 sanity
    logger.info("[Stage 4] Sanity check")

    sanity_lines: list[str] = []
    sanity_lines.append("# Step 0 Sanity Check\n")
    sanity_lines.append("> Generated by `scripts/build_augmented_master.py`\n")
    sanity_lines.append(
        f"> Augmented master TSV: rows={len(augmented):,}, cols={len(augmented.columns)}, "
        f"size={out_size_mb:.1f} MB\n\n"
    )

    # 4.1 Non-null rate per methylation col.
    sanity_lines.append("## 4.1 Non-null rate per methylation column\n\n")
    sanity_lines.append("| Column | Non-null % | Non-null count |\n|---|---:|---:|\n")
    nonnull_table = []
    for col in feature_cols_added:
        nn = augmented[col].notna().sum()
        pct = nn / len(augmented) * 100
        nonnull_table.append((col, pct, nn))
        sanity_lines.append(f"| {col} | {pct:.2f}% | {nn:,} |\n")
    min_nonnull = min(t[1] for t in nonnull_table)
    logger.info(f"  Min non-null rate across {len(feature_cols_added)} methyl cols: {min_nonnull:.2f}%")

    # 4.2 V5 vs V6 pearson (random 100 rows, then full).
    sanity_lines.append("\n## 4.2 V5 vs V6 Pearson r (same flag, same feature)\n\n")
    sanity_lines.append(
        "Expected r > 0.95 — V6 reuses V5 phased VCF; methylation features should be near-identical.\n\n"
    )
    sanity_lines.append("| Flag | Feature | n_paired | Pearson r (random 100) | Pearson r (full) |\n")
    sanity_lines.append("|---|---|---:|---:|---:|\n")
    np.random.seed(42)
    v5_v6_pairs = []
    for flag in FLAGS:
        for feat in METHYL_FEATURES + ["NME_imbalance", "Epipoly_imbalance"]:
            c5 = f"V5_{flag}_meth_{feat}"
            c6 = f"V6_{flag}_meth_{feat}"
            if c5 not in augmented.columns or c6 not in augmented.columns:
                continue
            valid = augmented[[c5, c6]].dropna()
            if len(valid) < 10:
                sanity_lines.append(f"| {flag} | {feat} | {len(valid)} | n/a | n/a |\n")
                continue
            full_r = safe_corr(valid[c5], valid[c6])
            n_sample = min(100, len(valid))
            sample = valid.sample(n_sample, random_state=42)
            rand_r = safe_corr(sample[c5], sample[c6])
            v5_v6_pairs.append((flag, feat, len(valid), rand_r, full_r))
            sanity_lines.append(
                f"| {flag} | {feat} | {len(valid):,} | {rand_r:.4f} | {full_r:.4f} |\n"
            )
    valid_full_r = [p[4] for p in v5_v6_pairs if not np.isnan(p[4])]
    if valid_full_r:
        median_v5v6_r = float(np.median(valid_full_r))
        logger.info(f"  V5 vs V6 median full Pearson r: {median_v5v6_r:.4f} (expect > 0.95)")
    else:
        median_v5v6_r = float("nan")

    # 4.3 V3F vs V5 (HP-conditional may differ).
    sanity_lines.append("\n## 4.3 V3F vs V5 Pearson r (Layer 1.5 effect on HP bucket — may differ)\n\n")
    sanity_lines.append("| Flag | Feature | n_paired | Pearson r (full) |\n|---|---|---:|---:|\n")
    v3f_v5_pairs = []
    for flag in FLAGS:
        for feat in METHYL_FEATURES + ["NME_imbalance", "Epipoly_imbalance"]:
            c3 = f"V3F_{flag}_meth_{feat}"
            c5 = f"V5_{flag}_meth_{feat}"
            if c3 not in augmented.columns or c5 not in augmented.columns:
                continue
            valid = augmented[[c3, c5]].dropna()
            if len(valid) < 10:
                sanity_lines.append(f"| {flag} | {feat} | {len(valid)} | n/a |\n")
                continue
            full_r = safe_corr(valid[c3], valid[c5])
            v3f_v5_pairs.append((flag, feat, len(valid), full_r))
            sanity_lines.append(f"| {flag} | {feat} | {len(valid):,} | {full_r:.4f} |\n")
    valid_v3v5_r = [p[3] for p in v3f_v5_pairs if not np.isnan(p[3])]
    if valid_v3v5_r:
        median_v3v5_r = float(np.median(valid_v3v5_r))
        logger.info(f"  V3F vs V5 median full Pearson r: {median_v3v5_r:.4f}")
    else:
        median_v3v5_r = float("nan")

    # 4.4 Verdicts.
    sanity_lines.append("\n## 4.4 Verdicts\n\n")
    verdict_h5 = "PASS" if (median_v5v6_r >= 0.95) else "FAIL"
    sanity_lines.append(f"- H5 sanity (V5 vs V6 r > 0.95): **{verdict_h5}** (median r = {median_v5v6_r:.4f})\n")
    sanity_lines.append(
        f"- V3F vs V5 median r = {median_v3v5_r:.4f} (Layer 1.5 HP bucket reshuffling expected to lower r)\n"
    )
    sanity_lines.append(f"- Min non-null across {len(feature_cols_added)} cols: {min_nonnull:.2f}%\n")
    sanity_lines.append(f"- 95% threshold met: {'PASS' if min_nonnull >= 95 else 'INFO (see col details)'}\n")

    SANITY_PATH.write_text("".join(sanity_lines))
    logger.info(f"  Wrote {SANITY_PATH}")

    # ------------------------------------------------------------------ Schema doc
    logger.info("[Stage 5] Write schema doc")
    schema_lines = ["# Step 0 — Augmented master TSV schema\n\n"]
    schema_lines.append(
        f"> Source: `scripts/build_augmented_master.py`  \n"
        f"> Rows: {len(augmented):,}  \n"
        f"> Cols: {len(augmented.columns)} (base {base_cols} + methyl {len(feature_cols_added)})  \n"
        f"> Size: {out_size_mb:.1f} MB  \n"
        f"> Output: `step5_master_augmented.tsv`\n\n"
    )
    schema_lines.append("## All columns (dtype + non-null %)\n\n")
    schema_lines.append("| # | Column | dtype | Non-null % |\n|---:|---|---|---:|\n")
    for i, col in enumerate(augmented.columns, 1):
        nn_pct = augmented[col].notna().mean() * 100
        schema_lines.append(f"| {i} | `{col}` | {augmented[col].dtype} | {nn_pct:.2f}% |\n")
    schema_lines.append("\n## Methylation column naming convention\n\n")
    schema_lines.append("- Pattern: `{BAM}_{flag}_meth_{feature}` where BAM ∈ {V3F, V5, V6}, flag ∈ {off, on}.\n")
    schema_lines.append(
        "- Derived imbalance cols: `*_meth_NME_imbalance = |NME_HP1 - NME_HP2|`, "
        "`*_meth_Epipoly_imbalance = |Epipoly_HP1 - Epipoly_HP2|`.\n"
    )
    SCHEMA_PATH.write_text("".join(schema_lines))
    logger.info(f"  Wrote {SCHEMA_PATH}")

    elapsed = time.time() - t0
    logger.info("=" * 70)
    logger.info(f"Step 0 COMPLETE in {elapsed:.1f} s ({elapsed/60:.2f} min)")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
