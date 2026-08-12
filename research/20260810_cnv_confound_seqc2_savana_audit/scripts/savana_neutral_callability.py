#!/usr/bin/env python3
"""Can SAVANA be trusted to identify copy-neutral ground on the samples with no truth?

Extending the clean-ground analysis to H1437 / H2009 / HCC1954 requires SAVANA
to answer one specific question correctly: "is this segment 1+1?"  That is a
much narrower demand than getting every copy number right, so it deserves its
own measurement rather than inheriting the general verdict.

HCC1395 is the only place this can be measured, because only it has SEQC2
truth.  Two SAVANA variants are scored against that truth:

  as-published   the fit SAVANA actually emitted (purity 0.76 / ploidy 1.83),
                 which the refit grid showed to be wrong
  recalibrated   copy number re-derived from SAVANA's own segmented log2 ratios
                 at the grid optimum (purity 1.0 / ploidy 2.95)

"Neutral" is defined as total CN 2 with both alleles present, since that is the
condition under which a clonal mutation on one haplotype reads out at 1/1 and
read-AF carries no copy-number component.  LOH is detected from the segment's
mean B-allele fraction rather than from SAVANA's minor-allele call, because the
latter is itself a function of the (possibly wrong) purity.

If recalibration does not rescue neutral-calling on HCC1395, the clean-ground
extension to the other samples cannot be trusted either, and that is the
finding rather than a reason to proceed.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict

SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
GAIN_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_gain_cn.bed"
LOSS_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_loss_cn.bed"
CATEGORICAL_BED = f"{SEQC2_DIR}/ngs_benchmark_cnvs_gain_loss_loh.bed"

WORK = "/big7_disk/liaoyoyo2001/cnv_sv_work"
HCC1395_CNA = f"{WORK}/HCC1395/savana_wgs/cna"
LOG2R = f"{HCC1395_CNA}/HCC1395_read_counts_mnorm_log2r_segmented.tsv"

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(OUT_DIR, "savana_neutral_callability.json")

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
BIN = 100_000
BAF_LOH = 0.75  # mean BAF at or above this is treated as loss of heterozygosity
CN_TOL = 0.5  # |CN - 2| <= tol counts as total CN 2

REFIT_PURITY, REFIT_PLOIDY = 1.0, 2.95
PUB_PURITY, PUB_PLOIDY = 0.76, 1.83


def load_bed(path, cast=float):
    segs = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            segs[p[0]].append((int(p[1]), int(p[2]), cast(p[3])))
    for c in segs:
        segs[c].sort()
    return segs


def bin_of_bed(segs, want=None):
    """chrom -> {bin: value} for the bins each segment covers."""
    out = defaultdict(dict)
    for chrom, items in segs.items():
        if chrom not in AUTOSOMES:
            continue
        for s, e, v in items:
            if want is not None and v != want:
                continue
            for b in range(s // BIN, (e - 1) // BIN + 1):
                out[chrom][b] = v
    return out


def cn_from_ratio(r, rho, psi):
    d = 2.0 * (1.0 - rho)
    return (r * (rho * psi + d) - d) / rho


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    gain = bin_of_bed(load_bed(GAIN_CN_BED))
    loss = bin_of_bed(load_bed(LOSS_CN_BED))
    loh = bin_of_bed(load_bed(CATEGORICAL_BED, cast=str), want="loh")

    # ---- SAVANA segment table: BAF per bin ----
    seg_path = f"{HCC1395_CNA}/HCC1395_segmented_absolute_copy_number.tsv"
    baf_bin = defaultdict(dict)
    pub_cn_bin = defaultdict(dict)
    with open(seg_path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            chrom = p[col["chromosome"]]
            if chrom not in AUTOSOMES:
                continue
            try:
                s, e = int(p[col["start"]]), int(p[col["end"]])
            except ValueError:
                continue
            try:
                baf = float(p[col["meanBAF"]])
            except (ValueError, IndexError):
                baf = None
            try:
                cnv = float(p[col["copyNumber"]])
            except (ValueError, IndexError):
                cnv = None
            for b in range(s // BIN, (e - 1) // BIN + 1):
                if baf is not None:
                    baf_bin[chrom][b] = baf
                if cnv is not None:
                    pub_cn_bin[chrom][b] = cnv

    # ---- recalibrated CN per bin from SAVANA's own log2r ----
    refit_cn_bin = defaultdict(lambda: defaultdict(list))
    with open(LOG2R) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            chrom = p[col["chromosome"]]
            if chrom not in AUTOSOMES:
                continue
            if p[col["use_bin"]].strip().lower() != "true":
                continue
            raw = p[col["seg_log2r_copynumber"]].strip()
            if not raw or raw.upper() in ("NA", "NAN"):
                continue
            try:
                l2r = float(raw)
                s = int(p[col["start"]])
            except ValueError:
                continue
            refit_cn_bin[chrom][s // BIN].append(2.0**l2r)
    refit_cn = {
        c: {b: cn_from_ratio(sum(v) / len(v), REFIT_PURITY, REFIT_PLOIDY)
            for b, v in bins.items()}
        for c, bins in refit_cn_bin.items()
    }

    # ---- truth neutral per bin ----
    def truth_neutral(chrom, b):
        if b in gain.get(chrom, {}) or b in loss.get(chrom, {}):
            return False
        if b in loh.get(chrom, {}):
            return False
        return True

    def score(cn_bins, label, use_baf=True):
        tp = fp = fn = tn = 0
        no_baf = 0
        for chrom in AUTOSOMES:
            for b, cnv in cn_bins.get(chrom, {}).items():
                truth = truth_neutral(chrom, b)
                is_cn2 = abs(cnv - 2.0) <= CN_TOL
                if use_baf:
                    baf = baf_bin.get(chrom, {}).get(b)
                    if baf is None:
                        no_baf += 1
                        called = is_cn2
                    else:
                        called = is_cn2 and baf < BAF_LOH
                else:
                    called = is_cn2
                if called and truth:
                    tp += 1
                elif called and not truth:
                    fp += 1
                elif not called and truth:
                    fn += 1
                else:
                    tn += 1
        prec = tp / (tp + fp) if (tp + fp) else None
        rec = tp / (tp + fn) if (tp + fn) else None
        f1 = (2 * prec * rec / (prec + rec)) if (prec and rec) else None
        return {
            "label": label,
            "bins_scored": tp + fp + fn + tn,
            "called_neutral": tp + fp,
            "truth_neutral": tp + fn,
            "true_positive": tp,
            "false_positive": fp,
            "false_negative": fn,
            "precision_percent": round(100 * prec, 2) if prec is not None else None,
            "recall_percent": round(100 * rec, 2) if rec is not None else None,
            "f1_percent": round(100 * f1, 2) if f1 is not None else None,
            "bins_without_baf": no_baf,
        }

    results = {
        "as_published_cn_only": score(pub_cn_bin, "published CN, no BAF gate", use_baf=False),
        "as_published_with_baf": score(pub_cn_bin, "published CN + BAF<%.2f" % BAF_LOH),
        "recalibrated_cn_only": score(refit_cn, "recalibrated CN, no BAF gate", use_baf=False),
        "recalibrated_with_baf": score(refit_cn, "recalibrated CN + BAF<%.2f" % BAF_LOH),
    }

    # ---- how much clean ground each variant would hand downstream ----
    truth_n = sum(
        1
        for chrom in AUTOSOMES
        for b in refit_cn.get(chrom, {})
        if truth_neutral(chrom, b)
    )

    out = {
        "generated_by": os.path.basename(__file__),
        "question": "on HCC1395, how accurately does SAVANA identify copy-neutral (1+1) ground, before and after recalibration?",
        "neutral_definition": {
            "total_cn": "abs(CN - 2) <= %.1f" % CN_TOL,
            "heterozygosity": "segment mean BAF < %.2f" % BAF_LOH,
            "why": "a clonal mutation on a single-copy haplotype reads 1/1 only when the haplotype carries exactly one copy",
        },
        "recalibration": {
            "published": {"purity": PUB_PURITY, "ploidy": PUB_PLOIDY},
            "refit": {"purity": REFIT_PURITY, "ploidy": REFIT_PLOIDY},
            "external_ploidy_references": {
                "seqc2_implied_bp_weighted": 2.9104,
                "refit_grid_optimum": 2.95,
                "knowledge_base_hcc1395": 2.85,
            },
        },
        "bin_size_bp": BIN,
        "truth_neutral_bins": truth_n,
        "scores": results,
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"truth-neutral bins (SEQC2): {truth_n}")
    for k, v in results.items():
        print(
            f"  {v['label']:38s} called={v['called_neutral']:6d} "
            f"P={v['precision_percent']}% R={v['recall_percent']}% F1={v['f1_percent']}%"
        )
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
