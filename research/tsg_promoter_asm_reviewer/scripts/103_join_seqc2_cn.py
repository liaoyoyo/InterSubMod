#!/usr/bin/env python3
"""
103 - Join every ISM-analysed locus to the SEQC2 copy-number GROUND TRUTH, so the
workstation can show the real SEQC2 CN answer (not just ISM's read-depth proxy).

SEQC2 truth files (/big8_disk/data/HCC1395/SEQC2/CNV/):
  ngs_benchmark_cnvs_gain_loss_loh.bed : chrom start end {gain|loss|loh} (660 regions)
  ngs_benchmark_cnv_gain_cn.bed        : chrom start end median_CN (gain)
  ngs_benchmark_cnv_loss_cn.bed        : chrom start end median_CN (loss)

Per locus -> SEQC2 class (neutral/gain/loss/loh) + SEQC2 median CN (2 if neutral;
gain/loss CN from the cn beds; loh = copy-neutral so CN≈2).

Output: display_v2/seqc2_cn.json  {key: [class_idx, median_cn]}  class 0neutral/1gain/2loss/3loh
"""
import csv, json, os
from bisect import bisect_right

CNV = "/big8_disk/data/HCC1395/SEQC2/CNV"
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
CLASS = {"neutral": 0, "gain": 1, "loss": 2, "loh": 3}


def load_bed(path, valcol=3, asfloat=False):
    by_chr = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            ch, s, e = p[0], int(p[1]), int(p[2])
            v = float(p[valcol]) if asfloat else p[valcol]
            by_chr.setdefault(ch, []).append((s, e, v))
    for ch in by_chr:
        by_chr[ch].sort()
    return by_chr


def lookup(by_chr, ch, pos):
    """return value of the region containing pos, else None (regions may overlap; first hit)."""
    regs = by_chr.get(ch)
    if not regs:
        return None
    # linear scan (≤~50 regions/chr) — find any region covering pos
    for s, e, v in regs:
        if s <= pos < e:
            return v
        if s > pos:
            break
    return None


def main():
    cls_bed = load_bed(f"{CNV}/ngs_benchmark_cnvs_gain_loss_loh.bed")          # class
    gain_cn = load_bed(f"{CNV}/ngs_benchmark_cnv_gain_cn.bed", asfloat=True)   # median CN
    loss_cn = load_bed(f"{CNV}/ngs_benchmark_cnv_loss_cn.bed", asfloat=True)
    out = {}
    counts = {0: 0, 1: 0, 2: 0, 3: 0}
    for cls in ("tp", "fp"):
        for r in csv.DictReader(open(f"{EX}/HCC1395_{cls}/significance_summary.csv")):
            ch, pos = r["Chr"], int(r["Pos"]); key = f"{ch}_{pos}"
            klass = lookup(cls_bed, ch, pos) or "neutral"
            ci = CLASS.get(klass, 0)
            if klass == "gain":
                cn = lookup(gain_cn, ch, pos)
            elif klass == "loss":
                cn = lookup(loss_cn, ch, pos)
            else:
                cn = 2.0  # neutral / loh = copy-neutral
            out[key] = [ci, round(cn if cn is not None else 2.0, 2)]
            counts[ci] += 1
    with open(f"{DV}/seqc2_cn.json", "w") as f:
        json.dump(out, f)
    print(f"[103] joined {len(out)} loci to SEQC2 CN truth")
    print(f"   SEQC2 class: neutral={counts[0]} gain={counts[1]} loss={counts[2]} loh={counts[3]}")
    print(f"[103] wrote {DV}/seqc2_cn.json")


if __name__ == "__main__":
    main()
