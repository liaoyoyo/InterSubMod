#!/usr/bin/env python3
"""
Step 0: Build consolidated 7-sample dataset under NEW KDE baseline.

Strategy:
  - 6 samples (HCC1395, HCC1395_DORADO, H2009, H1437, COLO829, HCC1937):
    use 2026-04-20/21 paired_full runs which contain Diploid_Coverage_Used
    (confirmed post-KDE-fix binary output)
  - HCC1954: fallback to stale 20260315 master with post-hoc CovM rescaling
    by ratio 1.2295 (stale 75 / new 61)
  - HCC1395 TO sidebar: Phase 1 fix dataset `/tmp/ism_hp_fix_phase1/{tp,fp}_off/`
    (post-fix Diploid_Coverage_Used available, flag=off only)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import (
    DATA_DIR,
    HCC1395_TO_FP,
    HCC1395_TO_TP,
    HCC1954_STALE_PAIRED_FULL,
    NEW_KDE_PAIRED_FULL,
    NEW_KDE_TO_ARCHIVE,
    OBS_DIR,
    SAMPLE_ORDER,
    get_sample_baseline,
    get_sample_ratio,
    load_baselines,
    load_sample_run,
    loh_side,
)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OBS_DIR.mkdir(parents=True, exist_ok=True)

    baselines = load_baselines()
    print(f"[step0] baselines loaded: {baselines.shape}")

    frames = []

    # 6 direct new-KDE paired_full samples
    for sample, run_dir in NEW_KDE_PAIRED_FULL.items():
        tp = Path(run_dir) / "intersubmod_tp" / "significance_summary.csv"
        fp = Path(run_dir) / "intersubmod_fp" / "significance_summary.csv"
        print(f"[step0] loading {sample} paired_full NEW KDE from {tp.parent.name}/..")
        df = load_sample_run(tp, fp, sample, "paired_full", "new_direct")
        df["CovM_used"] = df["Coverage_Multiple"].astype(float)
        df["baseline_x"] = float(df["Diploid_Coverage_Used"].median()) if "Diploid_Coverage_Used" in df.columns else np.nan
        frames.append(df)
        print(f"  -> n={len(df)}, NG dist: {df['HPFineNGroups'].value_counts().sort_index().to_dict()}")

    # HCC1954 fallback: stale + rescale by per-sample ratio
    stale_tp = Path(HCC1954_STALE_PAIRED_FULL) / "intersubmod_tp" / "significance_summary.csv"
    stale_fp = Path(HCC1954_STALE_PAIRED_FULL) / "intersubmod_fp" / "significance_summary.csv"
    print(f"[step0] loading HCC1954 paired_full STALE + rescale")
    df1954 = load_sample_run(stale_tp, stale_fp, "HCC1954", "paired_full", "stale_rescaled")
    ratio = get_sample_ratio("HCC1954", baselines, mode="paired_full")
    baseline_new = get_sample_baseline("HCC1954", baselines, mode="paired_full")
    df1954["CovM_used"] = df1954["Coverage_Multiple"].astype(float) * ratio
    df1954["baseline_x"] = baseline_new
    print(f"  -> n={len(df1954)}, ratio={ratio:.4f}, new_baseline={baseline_new}x")
    frames.append(df1954)

    # HCC1395 TO sidebar (Phase 1 fix data, flag=off)
    try:
        df_to = load_sample_run(Path(HCC1395_TO_TP), Path(HCC1395_TO_FP), "HCC1395", "to_pileup", "phase1_new")
        df_to["CovM_used"] = df_to["Coverage_Multiple"].astype(float)
        df_to["baseline_x"] = float(df_to["Diploid_Coverage_Used"].median()) if "Diploid_Coverage_Used" in df_to.columns else np.nan
        print(f"[step0] HCC1395 TO sidebar n={len(df_to)}, baseline={df_to['baseline_x'].iloc[0]}x")
        frames.append(df_to)
    except FileNotFoundError as e:
        print(f"[step0] HCC1395 TO sidebar skipped: {e}")

    # 2026-04-22: 5 archive TO pilots reran with post-KDE-fix binary
    for sample, step05_dir in NEW_KDE_TO_ARCHIVE.items():
        tp = Path(step05_dir) / "intersubmod_tp" / "significance_summary.csv"
        fp = Path(step05_dir) / "intersubmod_fp" / "significance_summary.csv"
        try:
            df_to_arch = load_sample_run(tp, fp, sample, "to_pileup", "archive_rerun_20260422")
            df_to_arch["CovM_used"] = df_to_arch["Coverage_Multiple"].astype(float)
            if "Diploid_Coverage_Used" in df_to_arch.columns:
                df_to_arch["baseline_x"] = float(df_to_arch["Diploid_Coverage_Used"].median())
            else:
                df_to_arch["baseline_x"] = np.nan
            print(f"[step0] {sample} TO archive-rerun n={len(df_to_arch)}, baseline={df_to_arch['baseline_x'].iloc[0]}x")
            frames.append(df_to_arch)
        except (FileNotFoundError, ValueError) as e:
            print(f"[step0] {sample} TO archive-rerun skipped: {e}")

    master = pd.concat(frames, ignore_index=True)

    # Derive analysis columns
    master["sample_order"] = master["sample"].apply(
        lambda s: SAMPLE_ORDER.index(s) if s in SAMPLE_ORDER else 99
    )
    master["loh_side"] = master.apply(loh_side, axis=1)

    # AlleleDelta equals variant AF (caller AF), used for AF class
    master["AF"] = master["AlleleDelta"].clip(lower=0.0, upper=1.0)
    master["AF_bin10"] = pd.cut(master["AF"], bins=np.arange(0.0, 1.01, 0.1), include_lowest=True)
    master["AF_class"] = master["AF"].apply(
        lambda x: "Extreme" if (pd.notna(x) and (x < 0.1 or x > 0.9))
        else ("Near-half" if (pd.notna(x) and 0.4 <= x <= 0.6) else ("Intermediate" if pd.notna(x) else "Unknown"))
    )

    # Summary stats
    print("\n[step0] per-sample summary:")
    summary = master.groupby(["sample", "mode", "kde_status"], observed=True).agg(
        n=("RegionID", "size"),
        n_tp=("tp_label", "sum"),
        covm_median=("CovM_used", "median"),
        baseline=("baseline_x", "median"),
        ng_mean=("HPFineNGroups", "mean"),
    ).reset_index()
    print(summary.to_string(index=False))

    # Verification: per-sample median CovM_used should match D3 fixed_p50 within 5%
    print("\n[step0] verification vs quantile_drift fixed_p50:")
    for _, row in summary.iterrows():
        s = row["sample"]; m = row["mode"]
        if m != "paired_full":
            continue
        sub = baselines[(baselines["sample"] == s) & (baselines["fixed_mode"] == "paired_full")]
        if sub.empty:
            continue
        d3_p50 = float(sub["fixed_p50"].iloc[0])
        ratio_error = abs(row["covm_median"] - d3_p50) / d3_p50
        flag = "OK" if ratio_error < 0.05 else "CHECK" if ratio_error < 0.15 else "FAIL"
        print(f"  {s}: our median={row['covm_median']:.3f} vs D3 p50={d3_p50:.3f} (err={ratio_error*100:.1f}% {flag})")

    # Write outputs
    out_tsv = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
    master.to_csv(out_tsv, sep="\t", index=False, compression="gzip")
    print(f"\n[step0] wrote {out_tsv} ({out_tsv.stat().st_size / 1e6:.1f} MB)")
    summary.to_csv(DATA_DIR / "step0_per_sample_summary.tsv", sep="\t", index=False)
    print(f"[step0] wrote per-sample summary table")


if __name__ == "__main__":
    main()
