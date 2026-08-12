#!/usr/bin/env python3
"""Is SAVANA's CN error a bad purity/ploidy fit, or a bad signal?

The adjudication step found SAVANA disagreeing with SEQC2 on ~70% of the
autosome, with a near-constant per-chromosome bias of about -1 copy and a
non-neutral regression R^2 of 0.87.  That pattern says the *segmentation* and
the *relative* signal are fine and the *absolute calibration* is off.

This script tests that directly.  It takes SAVANA's own segmented log2 ratios
-- the raw signal, before any purity/ploidy fitting -- and re-derives absolute
copy number over a grid of (purity, ploidy) values:

    depth_ratio r = [rho * n + 2(1-rho)] / [rho * psi + 2(1-rho)]
    =>          n = { r * [rho * psi + 2(1-rho)] - 2(1-rho) } / rho

For each grid point the re-derived CN is scored against SEQC2 truth.  If some
(purity, ploidy) recovers high agreement, the signal is good and the published
fit is simply wrong -- i.e. recoverable.  If no grid point does, the signal
itself is inadequate and the other samples' CN cannot be rescued by re-fitting.

Scoring is bin-count weighted over chr1-22 using SAVANA's own 10 kb bins.
"""

from __future__ import annotations

import json
import os
from collections import Counter, defaultdict

SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_gain_cn.bed"
LOSS_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_loss_cn.bed"

SAVANA_DIR = "/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna"
LOG2R_TSV = f"{SAVANA_DIR}/HCC1395_read_counts_mnorm_log2r_segmented.tsv"

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(OUT_DIR, "savana_refit_grid.json")

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
NEUTRAL_CN = 2.0

PUBLISHED_PURITY = 0.76
PUBLISHED_PLOIDY = 1.83


def load_seg_bed(path):
    segs = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            segs[p[0]].append((int(p[1]), int(p[2]), float(p[3])))
    for c in segs:
        segs[c].sort()
    return segs


def truth_lookup(gain, loss, chrom, start, end):
    """SEQC2 total CN at the midpoint of [start, end)."""
    mid = (start + end) // 2
    for s, e, cn in gain.get(chrom, []):
        if s <= mid < e:
            return cn
    for s, e, cn in loss.get(chrom, []):
        if s <= mid < e:
            return cn
    return NEUTRAL_CN


def state(cn):
    return "gain" if cn > 2.5 else "loss" if cn < 1.5 else "cn2"


def cn_from_ratio(r, rho, psi):
    denom_term = 2.0 * (1.0 - rho)
    return (r * (rho * psi + denom_term) - denom_term) / rho


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    gain = load_seg_bed(GAIN_CN_BED)
    loss = load_seg_bed(LOSS_CN_BED)

    # Collapse bins into (seg_log2r, truth_cn) -> count so the grid search is cheap.
    combos = Counter()
    n_bins = 0
    skipped = 0
    with open(LOG2R_TSV) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= col["seg_log2r_copynumber"]:
                continue
            chrom = p[col["chromosome"]]
            if chrom not in AUTOSOMES:
                continue
            if p[col["use_bin"]].strip().lower() != "true":
                skipped += 1
                continue
            raw = p[col["seg_log2r_copynumber"]].strip()
            if not raw or raw.upper() in ("NA", "NAN", ""):
                skipped += 1
                continue
            try:
                seg_l2r = float(raw)
            except ValueError:
                skipped += 1
                continue
            start, end = int(p[col["start"]]), int(p[col["end"]])
            t = truth_lookup(gain, loss, chrom, start, end)
            combos[(round(seg_l2r, 6), t)] += 1
            n_bins += 1

    print(f"usable bins: {n_bins} (skipped {skipped}); distinct combos: {len(combos)}")

    combo_list = [(2.0**l2r, t, c) for (l2r, t), c in combos.items()]
    total_bins = sum(c for _, _, c in combo_list)

    def score(rho, psi):
        agree = 0
        abs_err = 0.0
        int_match = 0
        for r, truth, cnt in combo_list:
            n = cn_from_ratio(r, rho, psi)
            if state(n) == state(truth):
                agree += cnt
            abs_err += abs(n - truth) * cnt
            if round(n) == round(truth):
                int_match += cnt
        return (
            agree / total_bins,
            abs_err / total_bins,
            int_match / total_bins,
        )

    results = []
    rho = 0.40
    while rho <= 1.0001:
        psi = 1.5
        while psi <= 4.5001:
            a, e, im = score(rho, psi)
            results.append(
                {
                    "purity": round(rho, 3),
                    "ploidy": round(psi, 2),
                    "state_agreement": round(a, 5),
                    "mean_abs_cn_error": round(e, 5),
                    "integer_match_rate": round(im, 5),
                }
            )
            psi += 0.05
        rho += 0.02

    best_agree = max(results, key=lambda r: r["state_agreement"])
    best_err = min(results, key=lambda r: r["mean_abs_cn_error"])
    best_int = max(results, key=lambda r: r["integer_match_rate"])

    pa, pe, pi_ = score(PUBLISHED_PURITY, PUBLISHED_PLOIDY)
    published = {
        "purity": PUBLISHED_PURITY,
        "ploidy": PUBLISHED_PLOIDY,
        "state_agreement": round(pa, 5),
        "mean_abs_cn_error": round(pe, 5),
        "integer_match_rate": round(pi_, 5),
    }

    # A pure-tumour cell line reading, for reference
    pure = {}
    for psi_try in (2.9, 3.0):
        a, e, im = score(1.0, psi_try)
        pure[f"purity=1.0,ploidy={psi_try}"] = {
            "state_agreement": round(a, 5),
            "mean_abs_cn_error": round(e, 5),
            "integer_match_rate": round(im, 5),
        }

    top = sorted(results, key=lambda r: -r["state_agreement"])[:15]

    out = {
        "generated_by": os.path.basename(__file__),
        "inputs": {
            "savana_log2r_tsv": LOG2R_TSV,
            "seqc2_gain_cn_bed": GAIN_CN_BED,
            "seqc2_loss_cn_bed": LOSS_CN_BED,
        },
        "model": "n = (r*(rho*psi + 2(1-rho)) - 2(1-rho))/rho ; r = 2^seg_log2r",
        "scope": "chr1-22, SAVANA 10kb bins with use_bin=True",
        "bins_scored": total_bins,
        "bins_skipped": skipped,
        "published_fit": published,
        "best_by_state_agreement": best_agree,
        "best_by_mean_abs_error": best_err,
        "best_by_integer_match": best_int,
        "pure_tumour_reference": pure,
        "top15_by_state_agreement": top,
        "grid": {
            "purity_range": [0.40, 1.00, 0.02],
            "ploidy_range": [1.50, 4.50, 0.05],
            "n_points": len(results),
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"published fit      : {published}")
    print(f"best by agreement  : {best_agree}")
    print(f"best by abs error  : {best_err}")
    print(f"best by int match  : {best_int}")
    print(f"pure-tumour ref    : {json.dumps(pure)}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
