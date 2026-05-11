#!/usr/bin/env python3
"""B5 Root cause investigation: COLO829 paired_full S1 fold=0.59x reversal.

Thread B S1 = LOH_Strong + Extreme AF.
Other samples: fold > 1 (FP removal > TP loss).
COLO829: fold = 0.59x (TP loss exceeds FP removal).

Hypotheses to distinguish:
  (a) cross-sample cell-definition inconsistency
  (b) COLO829-specific tumor biology
  (c) denominator artifact (baseline TP rate saturates → no headroom)

Input: research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz
Output:
  - docs/experiments/in_progress/2026/04/20260423_B5_COLO829_S1_Fold_Anomaly_01.md
  - research/ng_kde_rescaling/data/B5_colo829_s1_fold_detail.tsv
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MERGED = ROOT / "research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
OUT_TSV = ROOT / "research/ng_kde_rescaling/data/B5_colo829_s1_fold_detail.tsv"
OUT_DIST_TSV = ROOT / "research/ng_kde_rescaling/data/B5_colo829_vs_cohort_distributions.tsv"

# S1 definition: LOH_Strong ∩ Extreme AF
# Samples for comparison (exclude HCC1954 which is post-hoc flagged by user)
SAMPLES = ["COLO829", "H1437", "H2009", "HCC1395", "HCC1395_DORADO", "HCC1937"]


def tp_rate(df_sub: pd.DataFrame) -> tuple[float, int, int]:
    n = len(df_sub)
    if n == 0:
        return (float("nan"), 0, 0)
    tp = int((df_sub["tp_label"] == 1).sum())
    return (tp / n, tp, n)


def wilson_ci(tp: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = tp / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    margin = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (centre - margin, centre + margin)


def s1_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """For each sample in paired_full mode, compute two fold metrics.

    S1 cells = LOH_Strong ∩ Extreme AF. Interpretation in Thread B:
    S1 is a *retention* scheme (keep S1 cells as high-confidence).

    fold_precision = (TP/FP ratio in S1 cells) / (TP/FP ratio baseline)
        → > 1 means S1 enriches TPs relative to baseline
        → < 1 means S1 cells are LESS pure than baseline (reversal)
        This matches tpfp_per_sample_scheme_full.tsv "fold_improvement"

    fold_fp_over_tp_removal = (FP fraction in S1) / (TP fraction in S1)
        → diagnostic: does removing S1 strip more FPs than TPs?
    """
    rows = []
    for sample in SAMPLES:
        sub = df[(df["sample"] == sample) & (df["mode"] == "paired_full")]
        if sub.empty:
            continue

        base_tp_rate, base_tp, base_n = tp_rate(sub)
        base_fp = base_n - base_tp
        base_ratio = (base_tp / base_fp) if base_fp > 0 else float("nan")

        s1_mask = (sub["LOH_Subtype"] == "LOH_Strong") & (sub["AF_class"] == "Extreme")
        s1_cells = sub[s1_mask]
        s1_tp_rate, s1_tp, s1_n = tp_rate(s1_cells)
        s1_fp = s1_n - s1_tp
        s1_ratio = (s1_tp / s1_fp) if s1_fp > 0 else float("nan")

        # fold_precision: S1 cell precision / baseline precision (matches Thread B)
        fold_precision = (s1_ratio / base_ratio) if base_ratio and not np.isnan(base_ratio) else float("nan")

        # fold_fp_over_tp_removal: diagnostic
        fp_removed_frac = (s1_fp / base_fp) if base_fp > 0 else float("nan")
        tp_removed_frac = (s1_tp / base_tp) if base_tp > 0 else float("nan")
        fold_fp_over_tp = (fp_removed_frac / tp_removed_frac) if tp_removed_frac > 0 else float("nan")

        # Wilson CI on S1 tp rate
        s1_ci_lo, s1_ci_hi = wilson_ci(s1_tp, s1_n)
        base_ci_lo, base_ci_hi = wilson_ci(base_tp, base_n)

        # Headroom diagnostic: how much FP headroom exists at baseline?
        baseline_fp_rate = (base_fp / base_n) if base_n > 0 else float("nan")

        rows.append({
            "sample": sample,
            "baseline_n": base_n,
            "baseline_tp": base_tp,
            "baseline_fp": base_fp,
            "baseline_tp_rate": base_tp_rate,
            "baseline_tp_rate_ci_lo": base_ci_lo,
            "baseline_tp_rate_ci_hi": base_ci_hi,
            "baseline_fp_rate": baseline_fp_rate,
            "baseline_tp_fp_ratio": base_ratio,
            "s1_n": s1_n,
            "s1_tp": s1_tp,
            "s1_fp": s1_fp,
            "s1_tp_rate": s1_tp_rate,
            "s1_tp_rate_ci_lo": s1_ci_lo,
            "s1_tp_rate_ci_hi": s1_ci_hi,
            "s1_tp_fp_ratio": s1_ratio,
            "fold_precision_thread_b": fold_precision,
            "fold_fp_over_tp_removal": fold_fp_over_tp,
            "s1_frac_of_baseline": s1_n / base_n if base_n > 0 else float("nan"),
            "tp_rate_delta_s1_vs_base": s1_tp_rate - base_tp_rate,
        })
    return pd.DataFrame(rows)


def distribution_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Per sample distributions of LOH_Subtype, AF_class, Coverage_Category."""
    rows = []
    for sample in SAMPLES:
        sub = df[(df["sample"] == sample) & (df["mode"] == "paired_full")]
        if sub.empty:
            continue
        n = len(sub)

        # LOH_Subtype
        loh = sub["LOH_Subtype"].fillna("None").value_counts()
        for subtype, count in loh.items():
            rows.append({
                "sample": sample, "dimension": "LOH_Subtype",
                "value": str(subtype), "count": int(count),
                "frac": count / n, "total_n": n,
            })
        # AF_class
        afc = sub["AF_class"].fillna("NA").value_counts()
        for cls, count in afc.items():
            rows.append({
                "sample": sample, "dimension": "AF_class",
                "value": str(cls), "count": int(count),
                "frac": count / n, "total_n": n,
            })
        # Coverage_Category
        cov = sub["Coverage_Category"].fillna("NA").value_counts()
        for cls, count in cov.items():
            rows.append({
                "sample": sample, "dimension": "Coverage_Category",
                "value": str(cls), "count": int(count),
                "frac": count / n, "total_n": n,
            })
    return pd.DataFrame(rows)


def main() -> int:
    print(f"[B5] Reading {MERGED}", file=sys.stderr)
    df = pd.read_csv(MERGED, sep="\t", low_memory=False)
    print(f"[B5] Total rows: {len(df):,}", file=sys.stderr)

    # S1 core table
    s1_tbl = s1_analysis(df)
    s1_tbl.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6f")
    print(f"[B5] Wrote {OUT_TSV}", file=sys.stderr)
    print(s1_tbl.to_string(index=False), file=sys.stderr)

    # Distributions
    dist_tbl = distribution_analysis(df)
    dist_tbl.to_csv(OUT_DIST_TSV, sep="\t", index=False, float_format="%.6f")
    print(f"[B5] Wrote {OUT_DIST_TSV}", file=sys.stderr)

    # Summary JSON (for report embedding)
    summary = {
        "s1_fold_by_sample": {
            row["sample"]: {
                "fold_precision_thread_b": row["fold_precision_thread_b"],
                "baseline_tp_rate": row["baseline_tp_rate"],
                "s1_tp_rate": row["s1_tp_rate"],
                "s1_n": int(row["s1_n"]),
                "s1_frac_of_baseline": row["s1_frac_of_baseline"],
                "tp_rate_delta_s1_vs_base": row["tp_rate_delta_s1_vs_base"],
                "baseline_fp_rate": row["baseline_fp_rate"],
            }
            for _, row in s1_tbl.iterrows()
        },
    }
    out_json = OUT_TSV.with_suffix(".summary.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[B5] Wrote {out_json}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
