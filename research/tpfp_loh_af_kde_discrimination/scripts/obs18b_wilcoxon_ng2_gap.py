#!/usr/bin/env python3
"""obs18b — committed Wilcoxon for NG=2 Inner same_HP1 vs Outer cross_het TP-rate gap.

PROVENANCE FIX (2026-05-30): the original B1 Wilcoxon (W=21, p=0.0156, dated 2026-04-23)
was computed ad-hoc and NEVER committed — only its JSON output existed
(data/obs18_wilcoxon_B1.json). This script re-implements that computation from the
committed composition TSV so the Grade-A statistic is reproducible, and is n-agnostic
so COLO829 (or any 7th sample) drops straight in for the n=7 cross-sample Wilcoxon.

Gap per sample = (Inner same_HP1 TP rate) - (Outer cross_het TP rate).
Test: scipy.stats.wilcoxon(gaps, alternative='greater', zero_method='wilcox', method='exact').

Usage:
  python3 obs18b_wilcoxon_ng2_gap.py                 # n=6 (reproduce + verify B1)
  python3 obs18b_wilcoxon_ng2_gap.py --extra COLO829:/path/to/colo829_composition_rows.tsv
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination")
COMP_TSV = ROOT / "data/obs18_NG2_composition_by_sample.tsv"
OUT_JSON = ROOT / "data/obs18b_wilcoxon_ng2_gap.json"

INNER_SAME_HP1 = ("Inner", "same_HP1 (HP1 + HP1-1)")
OUTER_CROSS_HET = ("Outer", "cross_het (HP1 + HP2-1)")

# B1 (2026-04-23) reference for reproducibility verification
B1_REF = {
    "W": 21.0, "p": 0.015625, "median_gap": 0.36503531799845024,
    "samples": ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437"],
}


def load_composition(tsv_path):
    """Return {sample: {(loh_side, combo): {'tp_rate': float, 'n': int, 'n_tp': int}}}."""
    by_sample = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            s = row["sample"]
            key = (row["loh_side"], row["combo"])
            by_sample.setdefault(s, {})[key] = {
                "tp_rate": float(row["tp_rate"]) if row["tp_rate"] else float("nan"),
                "n": int(row["n"]), "n_tp": int(row["n_tp"]),
            }
    return by_sample


def per_sample_gap(comp):
    """Compute Inner-same_HP1 minus Outer-cross_het TP-rate gap per sample (in TSV order)."""
    out = []
    for s, buckets in comp.items():
        inner = buckets.get(INNER_SAME_HP1)
        outer = buckets.get(OUTER_CROSS_HET)
        if not inner or not outer:
            out.append({"sample": s, "gap": None, "reason": "missing bucket"})
            continue
        gap = inner["tp_rate"] - outer["tp_rate"]
        out.append({
            "sample": s, "gap": gap,
            "inner_same_HP1_tp_rate": inner["tp_rate"], "inner_same_HP1_n": inner["n"],
            "outer_cross_het_tp_rate": outer["tp_rate"], "outer_cross_het_n": outer["n"],
        })
    return out


def bootstrap_median_ci(gaps, n_resamples=1000, seed=20260423):
    rng = np.random.default_rng(seed)
    g = np.asarray(gaps, dtype=float)
    meds = [np.median(rng.choice(g, size=len(g), replace=True)) for _ in range(n_resamples)]
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comp", default=str(COMP_TSV))
    ap.add_argument("--extra", action="append", default=[],
                    help="SAMPLE:composition_rows.tsv to append a 7th+ sample (same TSV schema)")
    ap.add_argument("--out", default=str(OUT_JSON))
    args = ap.parse_args()

    comp = load_composition(args.comp)
    # append extra samples (e.g. COLO829) if provided, same schema
    for spec in args.extra:
        name, path = spec.split(":", 1)
        extra = load_composition(path)
        for s, buckets in extra.items():
            comp[s] = buckets

    rows = per_sample_gap(comp)
    valid = [r for r in rows if r.get("gap") is not None]
    gaps = [r["gap"] for r in valid]
    n = len(gaps)

    wil = stats.wilcoxon(gaps, alternative="greater", zero_method="wilcox", method="exact")
    blo, bhi = bootstrap_median_ci(gaps)
    result = {
        "analysis": "obs18b_committed_Wilcoxon_NG2_Inner_same_HP1_vs_Outer_cross_het",
        "provenance_note": "committed reimplementation of the never-committed B1 (2026-04-23)",
        "date": "2026-05-30",
        "source_tsv": args.comp,
        "extra_samples": args.extra,
        "n_samples": n,
        "samples": [r["sample"] for r in valid],
        "per_sample": valid,
        "gaps": gaps,
        "wilcoxon": {
            "statistic_W": float(wil.statistic), "p_value": float(wil.pvalue),
            "alternative": "greater", "zero_method": "wilcox", "method": "exact",
            "significant_at_0.05": bool(wil.pvalue < 0.05),
        },
        "descriptives": {
            "median_gap": float(np.median(gaps)), "mean_gap": float(np.mean(gaps)),
            "min_gap": float(np.min(gaps)), "max_gap": float(np.max(gaps)),
            "n_positive": int(sum(g > 0 for g in gaps)),
        },
        "bootstrap_median_gap_95ci": {"ci_low": blo, "ci_high": bhi, "n_resamples": 1000, "seed": 20260423},
        "verdict": "POSITIVE" if (wil.pvalue < 0.05 and all(g > 0 for g in gaps)) else "CHECK",
    }

    # reproducibility check vs B1 (only meaningful when n==6 default set)
    if n == 6 and set(result["samples"]) == set(B1_REF["samples"]):
        reproduces = (abs(result["wilcoxon"]["statistic_W"] - B1_REF["W"]) < 1e-9
                      and abs(result["wilcoxon"]["p_value"] - B1_REF["p"]) < 1e-9
                      and abs(result["descriptives"]["median_gap"] - B1_REF["median_gap"]) < 1e-6)
        result["reproduces_B1_2026_04_23"] = reproduces

    Path(args.out).write_text(json.dumps(result, ensure_ascii=False, indent=2))
    print(f"n={n} samples: {result['samples']}")
    print(f"gaps: {[round(g,4) for g in gaps]}")
    print(f"Wilcoxon: W={result['wilcoxon']['statistic_W']} p={result['wilcoxon']['p_value']:.6f} "
          f"({'sig' if result['wilcoxon']['significant_at_0.05'] else 'ns'})")
    print(f"median_gap={result['descriptives']['median_gap']:.5f} "
          f"bootstrap CI=[{blo:.5f},{bhi:.5f}] n_positive={result['descriptives']['n_positive']}/{n}")
    print(f"verdict: {result['verdict']}")
    if "reproduces_B1_2026_04_23" in result:
        print(f"reproduces B1 (2026-04-23 W=21 p=0.0156): {result['reproduces_B1_2026_04_23']}")
    print(f"written: {args.out}")


if __name__ == "__main__":
    main()
