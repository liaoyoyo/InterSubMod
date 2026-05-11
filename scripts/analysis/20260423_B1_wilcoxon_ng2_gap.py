#!/usr/bin/env python3
"""B1 Wilcoxon signed-rank test on Inner same_HP1 vs Outer cross_het TP gap.

Tests Thread D observation: NG=2 LOH-constrained phasing produces
same-haplotype inner bias across 6 TO samples.

Input: research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv
Output: obs18_wilcoxon_B1.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TSV_PATH = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv"
)
OUT_JSON = TSV_PATH.parent / "obs18_wilcoxon_B1.json"

SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437"]
INNER_KEY = "same_HP1 (HP1 + HP1-1)"
OUTER_KEY = "cross_het (HP1 + HP2-1)"


def main() -> int:
    df = pd.read_csv(TSV_PATH, sep="\t")

    rows = []
    for sample in SAMPLES:
        inner_mask = (df["sample"] == sample) & (df["loh_side"] == "Inner") & (df["combo"] == INNER_KEY)
        outer_mask = (df["sample"] == sample) & (df["loh_side"] == "Outer") & (df["combo"] == OUTER_KEY)
        inner_rows = df[inner_mask]
        outer_rows = df[outer_mask]
        if len(inner_rows) != 1 or len(outer_rows) != 1:
            print(f"ERROR: sample {sample} inner={len(inner_rows)} outer={len(outer_rows)}", file=sys.stderr)
            return 1
        inner_tp = float(inner_rows.iloc[0]["tp_rate"])
        outer_tp = float(outer_rows.iloc[0]["tp_rate"])
        inner_n = int(inner_rows.iloc[0]["n"])
        outer_n = int(outer_rows.iloc[0]["n"])
        rows.append({
            "sample": sample,
            "inner_same_HP1_tp_rate": inner_tp,
            "inner_same_HP1_n": inner_n,
            "outer_cross_het_tp_rate": outer_tp,
            "outer_cross_het_n": outer_n,
            "gap": inner_tp - outer_tp,
        })

    per_sample = pd.DataFrame(rows)
    gaps = per_sample["gap"].values.astype(float)

    # Wilcoxon signed-rank, alternative='greater' (inner > outer)
    wres = stats.wilcoxon(gaps, alternative="greater", zero_method="wilcox", method="exact")

    # Bootstrap 95% CI on median gap (1000 resamples)
    rng = np.random.default_rng(20260423)
    B = 1000
    medians = np.empty(B)
    n = len(gaps)
    for b in range(B):
        idx = rng.integers(0, n, n)
        medians[b] = np.median(gaps[idx])
    ci_lo = float(np.percentile(medians, 2.5))
    ci_hi = float(np.percentile(medians, 97.5))
    median_gap = float(np.median(gaps))
    mean_gap = float(np.mean(gaps))

    W = float(wres.statistic)
    pval = float(wres.pvalue)
    positive = bool(pval < 0.05)

    payload = {
        "analysis": "B1_Wilcoxon_NG2_Inner_same_HP1_vs_Outer_cross_het",
        "date": "2026-04-23",
        "source_tsv": str(TSV_PATH),
        "n_samples": n,
        "samples": SAMPLES,
        "per_sample": rows,
        "gaps": gaps.tolist(),
        "wilcoxon": {
            "statistic_W": W,
            "p_value": pval,
            "alternative": "greater",
            "zero_method": "wilcox",
            "method": "exact",
            "significant_at_0.05": positive,
        },
        "descriptives": {
            "median_gap": median_gap,
            "mean_gap": mean_gap,
            "min_gap": float(np.min(gaps)),
            "max_gap": float(np.max(gaps)),
        },
        "bootstrap_median_gap_95ci": {
            "ci_low": ci_lo,
            "ci_high": ci_hi,
            "n_resamples": B,
            "seed": 20260423,
        },
        "verdict": "POSITIVE" if positive else "NEGATIVE",
    }

    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
