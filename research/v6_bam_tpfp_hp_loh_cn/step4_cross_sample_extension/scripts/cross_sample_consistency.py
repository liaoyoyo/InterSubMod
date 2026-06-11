#!/usr/bin/env python3
"""Step 4 stage 3 — cross-sample consistency analysis (n=5).

For each cell (LOH × HP × cov), compute across HCC1395 + 4 phaseD samples:
  - direction concordance: # samples with TP rate > global baseline
  - Wilcoxon signed-rank test (TP rate vs each sample's own global baseline)
  - powered_count: # samples where cell is powered
  - signature_candidate flag: ≥4/5 same direction AND Wilcoxon p < 0.0625
    (n=5 exact min p = 0.0625 — only achievable when all 5 same-sign)
  - effect_size_mean: mean(TP_rate - global_TP_rate) across samples
  - n_eff: mean sample size of the cell across samples

Outputs:
  step4_cross_sample_extension/step4_consistency.tsv
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _common import (
    COV_LEVELS,
    HP_LEVELS,
    INTERMEDIATE_DIR,
    LOH_LEVELS,
    PER_SAMPLE_GRID,
    SAMPLES,
    STEP4_DIR,
    cell_id,
)


def main():
    # Combined grid
    combined_path = STEP4_DIR / "step4_per_sample_grid.tsv"
    df = pd.read_csv(combined_path, sep="\t", low_memory=False)
    # Load per-sample global TP rates from summary json
    summary_path = INTERMEDIATE_DIR / "step4_grid_summary.json"
    with open(summary_path) as f:
        summaries = json.load(f)
    global_rates = {s["sample"]: s["global_TP_rate"] for s in summaries}

    rows = []
    for loh in LOH_LEVELS:
        for hp in HP_LEVELS:
            for cov in COV_LEVELS:
                cid = cell_id(loh, hp, cov)
                cell_rows = df[df["cell_id"] == cid]
                sample_data = {
                    row["sample"]: {
                        "TP_rate": row["TP_rate"],
                        "n": row["n"],
                        "powered_flag": bool(row["powered_flag"]),
                        "marginal_flag": bool(row["marginal_flag"]),
                    }
                    for _, row in cell_rows.iterrows()
                }
                # Only consider samples where cell has n >= 10 (non-trivial)
                valid_samples = [s for s in SAMPLES if s in sample_data and sample_data[s]["n"] >= 10]
                rates = np.array(
                    [sample_data[s]["TP_rate"] for s in valid_samples], dtype=float
                )
                baselines = np.array(
                    [global_rates.get(s, np.nan) for s in valid_samples], dtype=float
                )
                deltas = rates - baselines
                n_samples = len(valid_samples)
                n_above = int(np.sum(deltas > 0))
                n_below = int(np.sum(deltas < 0))
                # Direction concordance: max(above, below) / n
                direction_match = max(n_above, n_below)

                # Wilcoxon signed-rank vs zero (n=5 exact)
                wilcoxon_p = float("nan")
                wilcoxon_stat = float("nan")
                if n_samples >= 3 and np.any(deltas != 0):
                    try:
                        w = stats.wilcoxon(
                            deltas[~np.isnan(deltas) & (deltas != 0)],
                            zero_method="wilcox",
                            alternative="two-sided",
                            method="exact",
                        )
                        wilcoxon_stat = float(w.statistic)
                        wilcoxon_p = float(w.pvalue)
                    except Exception:
                        pass

                # Sign-direction (majority sign across non-zero deltas)
                if n_above > n_below:
                    sign = "above_global"
                elif n_below > n_above:
                    sign = "below_global"
                else:
                    sign = "tie"

                # Signature candidate: n=5 with all 5 valid + ≥4/5 same direction + p<0.0625
                # (Wilcoxon n=5 exact two-sided min p = 0.0625)
                powered_count = sum(
                    sample_data[s]["powered_flag"] for s in valid_samples
                )
                signature_flag = (
                    n_samples == 5
                    and direction_match >= 4
                    and not np.isnan(wilcoxon_p)
                    and wilcoxon_p <= 0.0625
                )
                rec = {
                    "cell_id": cid,
                    "loh_side": loh,
                    "hp_bucket": hp,
                    "cov_bin": cov,
                    "n_samples_valid": n_samples,
                    "n_above": n_above,
                    "n_below": n_below,
                    "direction_match": direction_match,
                    "majority_sign": sign,
                    "wilcoxon_stat": wilcoxon_stat,
                    "wilcoxon_p": wilcoxon_p,
                    "mean_delta_vs_global": float(np.nanmean(deltas)) if len(deltas) else float("nan"),
                    "median_delta_vs_global": float(np.nanmedian(deltas)) if len(deltas) else float("nan"),
                    "mean_n": float(np.mean([sample_data[s]["n"] for s in valid_samples])) if valid_samples else 0.0,
                    "powered_count": powered_count,
                    "signature_candidate_flag": signature_flag,
                }
                # Add per-sample rates and ns for transparency
                for s in SAMPLES:
                    rec[f"TP_rate_{s}"] = sample_data.get(s, {}).get("TP_rate", float("nan"))
                    rec[f"n_{s}"] = sample_data.get(s, {}).get("n", 0)
                rows.append(rec)
    cons = pd.DataFrame(rows)
    out = STEP4_DIR / "step4_consistency.tsv"
    cons.to_csv(out, sep="\t", index=False, float_format="%.6f")
    print(f"wrote {out} (cells={len(cons)})")

    # Summary stats
    n_signature = int(cons["signature_candidate_flag"].sum())
    print(f"Signature candidates: {n_signature}")
    print(
        cons[cons["signature_candidate_flag"]][
            [
                "cell_id",
                "n_samples_valid",
                "direction_match",
                "majority_sign",
                "wilcoxon_p",
                "mean_delta_vs_global",
                "mean_n",
                "powered_count",
            ]
        ].to_string(index=False)
        if n_signature
        else "  (none)"
    )

    # Save summary
    summary = {
        "n_cells_total": len(cons),
        "n_cells_with_5_samples": int((cons["n_samples_valid"] == 5).sum()),
        "n_cells_4plus_concordant": int((cons["direction_match"] >= 4).sum()),
        "n_cells_4of5_plus_wilcoxon": int(
            ((cons["n_samples_valid"] == 5) & (cons["direction_match"] >= 4) & (cons["wilcoxon_p"] <= 0.0625)).sum()
        ),
        "n_signature_candidates": n_signature,
    }
    with open(INTERMEDIATE_DIR / "step4_consistency_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
