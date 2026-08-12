#!/usr/bin/env python3
"""Adjudicate the SAVANA CN pipeline against SEQC2 benchmark truth on HCC1395.

HCC1395 is the only cohort sample with an external CN truth set, so it is the
only place where the SAVANA pipeline itself can be validated.  Whatever this
script finds propagates to H1437 / H2009 / HCC1954, whose CN annotation comes
from SAVANA with no external check.

Three questions, in order:

1. How much of the genome does SAVANA call correctly (bp-weighted)?
2. Is any disagreement a *scale* error (purity/ploidy mis-fit, recoverable by
   re-calibration) or a *segmentation* error (breakpoints in the wrong place,
   not recoverable)?
3. What ploidy does SEQC2 truth itself imply, and does SAVANA's fitted 1.83
   agree?

The scale-vs-segmentation split matters because it decides whether the other
samples' CN can be rescued by re-fitting or must be discarded.
"""

from __future__ import annotations

import json
import os
from collections import Counter, defaultdict

SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_gain_cn.bed"
LOSS_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_loss_cn.bed"
CATEGORICAL_BED = f"{SEQC2_DIR}/ngs_benchmark_cnvs_gain_loss_loh.bed"

SAVANA_DIR = "/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna"
SAVANA_CN_TSV = f"{SAVANA_DIR}/HCC1395_segmented_absolute_copy_number.tsv"
SAVANA_PURITY_TSV = f"{SAVANA_DIR}/HCC1395_fitted_purity_ploidy.tsv"

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(OUT_DIR, "savana_vs_seqc2_adjudication.json")

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
BIN = 100_000  # 100 kb comparison bins
NEUTRAL_CN = 2.0


def load_seg_bed(path, cast=float):
    segs = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            segs[p[0]].append((int(p[1]), int(p[2]), cast(p[3])))
    return segs


def load_savana(path):
    segs = defaultdict(list)
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                cn = float(p[col["copyNumber"]])
            except (ValueError, IndexError, KeyError):
                cn = None
            try:
                minor = float(p[col["minorAlleleCopyNumber"]])
            except (ValueError, IndexError, KeyError):
                minor = None
            segs[p[col["chromosome"]]].append(
                (int(p[col["start"]]), int(p[col["end"]]), cn, minor)
            )
    return segs


def build_bin_map(segs, value_index=2, default=None):
    """chrom -> {bin_index: bp-weighted value} using 100 kb bins."""
    acc = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0]))  # [sum(v*bp), bp]
    for chrom, items in segs.items():
        if chrom not in AUTOSOMES:
            continue
        for it in items:
            s, e = it[0], it[1]
            v = it[value_index] if len(it) > value_index else None
            if v is None:
                continue
            b0, b1 = s // BIN, (e - 1) // BIN
            for b in range(b0, b1 + 1):
                lo = max(s, b * BIN)
                hi = min(e, (b + 1) * BIN)
                bp = hi - lo
                if bp > 0:
                    cell = acc[chrom][b]
                    cell[0] += v * bp
                    cell[1] += bp
    out = {}
    for chrom, bins in acc.items():
        out[chrom] = {
            b: (tot / bp) for b, (tot, bp) in bins.items() if bp > 0
        }
    return out


def linreg(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0:
        return None
    slope = sxy / sxx
    intercept = my - slope * mx
    syy = sum((y - my) ** 2 for y in ys)
    r2 = (sxy**2 / (sxx * syy)) if syy > 0 else None
    return {
        "n": n,
        "slope": round(slope, 4),
        "intercept": round(intercept, 4),
        "r_squared": round(r2, 4) if r2 is not None else None,
    }


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    gain = load_seg_bed(GAIN_CN_BED)
    loss = load_seg_bed(LOSS_CN_BED)
    cat = load_seg_bed(CATEGORICAL_BED, cast=str)
    sav = load_savana(SAVANA_CN_TSV)

    with open(SAVANA_PURITY_TSV) as fh:
        hdr = fh.readline().split()
        vals = fh.readline().split()
        savana_fit = dict(zip(hdr, vals))

    # ---- SEQC2 truth CN per bin (default CN=2 where no gain/loss record) ----
    gain_bins = build_bin_map(gain)
    loss_bins = build_bin_map(loss)
    sav_bins = build_bin_map(sav, value_index=2)
    sav_minor_bins = build_bin_map(sav, value_index=3)

    # LOH bins
    loh_bins = defaultdict(set)
    for chrom, items in cat.items():
        if chrom not in AUTOSOMES:
            continue
        for s, e, label in items:
            if label != "loh":
                continue
            for b in range(s // BIN, (e - 1) // BIN + 1):
                loh_bins[chrom].add(b)

    xs, ys = [], []
    xs_nonneutral, ys_nonneutral = [], []
    state_confusion = Counter()
    bp_agree = 0
    bp_total = 0
    per_chrom_bias = defaultdict(lambda: [0.0, 0])
    seqc2_cn_bp = 0.0
    seqc2_bp = 0

    # union of bins that SAVANA covers (SAVANA defines the comparable genome)
    for chrom in AUTOSOMES:
        sbins = sav_bins.get(chrom, {})
        gb = gain_bins.get(chrom, {})
        lb = loss_bins.get(chrom, {})
        for b, sav_cn in sbins.items():
            if b in gb:
                truth = gb[b]
            elif b in lb:
                truth = lb[b]
            else:
                truth = NEUTRAL_CN
            xs.append(truth)
            ys.append(sav_cn)
            if truth != NEUTRAL_CN:
                xs_nonneutral.append(truth)
                ys_nonneutral.append(sav_cn)

            t_state = "gain" if truth > 2.5 else "loss" if truth < 1.5 else "cn2"
            s_state = "gain" if sav_cn > 2.5 else "loss" if sav_cn < 1.5 else "cn2"
            state_confusion[(t_state, s_state)] += 1
            bp_total += 1
            if t_state == s_state:
                bp_agree += 1
            per_chrom_bias[chrom][0] += sav_cn - truth
            per_chrom_bias[chrom][1] += 1

        # SEQC2-implied ploidy over the whole autosome, independent of SAVANA
        chrom_bins = set(gb) | set(lb)
        # count every bin SAVANA saw, so the denominator is comparable
        for b in sbins:
            if b in gb:
                seqc2_cn_bp += gb[b]
            elif b in lb:
                seqc2_cn_bp += lb[b]
            else:
                seqc2_cn_bp += NEUTRAL_CN
            seqc2_bp += 1

    fit_all = linreg(xs, ys)
    fit_nonneutral = linreg(xs_nonneutral, ys_nonneutral)

    # What scale factor would best map SAVANA onto SEQC2 (no intercept)?
    num = sum(x * y for x, y in zip(xs, ys))
    den = sum(y * y for y in ys)
    best_scale = num / den if den else None

    # Re-score agreement after applying the best pure-scale correction
    rescaled_agree = 0
    rescaled_conf = Counter()
    if best_scale:
        for truth, sav_cn in zip(xs, ys):
            r = sav_cn * best_scale
            t_state = "gain" if truth > 2.5 else "loss" if truth < 1.5 else "cn2"
            s_state = "gain" if r > 2.5 else "loss" if r < 1.5 else "cn2"
            rescaled_conf[(t_state, s_state)] += 1
            if t_state == s_state:
                rescaled_agree += 1

    seqc2_implied_ploidy = seqc2_cn_bp / seqc2_bp if seqc2_bp else None

    out = {
        "generated_by": os.path.basename(__file__),
        "inputs": {
            "seqc2_gain_cn_bed": GAIN_CN_BED,
            "seqc2_loss_cn_bed": LOSS_CN_BED,
            "savana_cn_tsv": SAVANA_CN_TSV,
            "savana_purity_tsv": SAVANA_PURITY_TSV,
        },
        "bin_size_bp": BIN,
        "scope": "chr1-22, bins where SAVANA reports a copy number",
        "savana_fitted": savana_fit,
        "seqc2_implied_ploidy_autosomal_bp_weighted": (
            round(seqc2_implied_ploidy, 4) if seqc2_implied_ploidy else None
        ),
        "bins_compared": bp_total,
        "state_agreement": {
            "bins_agreeing": bp_agree,
            "denominator": bp_total,
            "percent": round(100.0 * bp_agree / bp_total, 2) if bp_total else None,
        },
        "state_confusion_bins": {
            f"seqc2={a}|savana={b}": c for (a, b), c in sorted(state_confusion.items())
        },
        "regression_savana_on_seqc2": fit_all,
        "regression_savana_on_seqc2_nonneutral_only": fit_nonneutral,
        "best_pure_scale_factor_savana_to_seqc2": (
            round(best_scale, 4) if best_scale else None
        ),
        "state_agreement_after_rescale": {
            "bins_agreeing": rescaled_agree,
            "denominator": bp_total,
            "percent": round(100.0 * rescaled_agree / bp_total, 2) if bp_total else None,
        },
        "state_confusion_bins_after_rescale": {
            f"seqc2={a}|savana={b}": c for (a, b), c in sorted(rescaled_conf.items())
        },
        "per_chrom_mean_bias_savana_minus_seqc2": {
            c: round(v[0] / v[1], 4) for c, v in sorted(per_chrom_bias.items()) if v[1]
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"bins compared: {bp_total}")
    print(f"SAVANA fitted: {savana_fit}")
    print(f"SEQC2-implied autosomal ploidy: {seqc2_implied_ploidy:.4f}")
    print(
        f"state agreement: {bp_agree}/{bp_total} = "
        f"{100.0 * bp_agree / bp_total:.2f}%"
    )
    print(f"regression (all bins): {fit_all}")
    print(f"best pure scale: {best_scale:.4f}" if best_scale else "no scale")
    print(
        f"agreement after rescale: {rescaled_agree}/{bp_total} = "
        f"{100.0 * rescaled_agree / bp_total:.2f}%"
    )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
