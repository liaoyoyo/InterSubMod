#!/usr/bin/env python3
"""Pilot: can the blocking input (HP -> allele mapping) actually be derived?

The spec-readiness audit named one blocking input: knowing how many copies the
*observed haplotype* carries in a segment, which needs the link between the
germline HP families and the segment's major/minor allele.

The phased germline VCF turns out to carry everything needed:

    GT   0|1 or 1|0   which haplotype holds ref and which holds alt
    AD   ref,alt      allelic depths
    PS                phase set

so for a phased het site the per-haplotype depth is read off directly --
GT 0|1 means HP1 is supported by the ref count and HP2 by the alt count, and
1|0 is the reverse.  Aggregating those over a copy-number segment gives the
haplotype depth ratio, which maps onto major/minor copy number.

This pilot does not build the mapping for production.  It tests the one thing
that would sink the idea: whether the haplotype depth ratio actually separates
segments that SEQC2 calls LOH from those it does not.  If LOH segments do not
show extreme haplotype imbalance, the derivation does not work and the blocking
input stays blocked.
"""

from __future__ import annotations

import gzip
import json
import math
import os
from collections import defaultdict

SEQC2_DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
CATEGORICAL_BED = f"{SEQC2_DIR}/ngs_benchmark_cnvs_gain_loss_loh.bed"
GAIN_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_gain_cn.bed"
LOSS_CN_BED = f"{SEQC2_DIR}/ngs_benchmark_cnv_loss_cn.bed"

VCF = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
    "20260420_HCC1395_paired_full_full_2/longphase_s/germline_phased_merged.vcf.gz"
)
SAVANA = (
    "/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna/"
    "HCC1395_segmented_absolute_copy_number.tsv"
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(OUT_DIR, "hp_to_allele_pilot.json")

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
MIN_DP = 10
MIN_SITES_PER_SEG = 20


def load_bed(path, cast=float):
    segs = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4 or p[0] not in AUTOSOMES:
                continue
            segs[p[0]].append((int(p[1]), int(p[2]), cast(p[3])))
    for c in segs:
        segs[c].sort()
    return segs


def build_lookup(segs, want=None):
    out = defaultdict(list)
    for c, items in segs.items():
        for s, e, v in items:
            if want is None or v == want:
                out[c].append((s, e, v))
    return out


def in_any(lookup, chrom, pos):
    for s, e, v in lookup.get(chrom, []):
        if s <= pos < e:
            return v
    return None


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    cat = load_bed(CATEGORICAL_BED, cast=str)
    loh = build_lookup(cat, want="loh")
    gain = build_lookup(load_bed(GAIN_CN_BED))
    loss = build_lookup(load_bed(LOSS_CN_BED))

    # SAVANA segments as the aggregation frame
    segs = []
    with open(SAVANA) as fh:
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

            def num(key):
                try:
                    return float(p[col[key]])
                except (ValueError, IndexError, KeyError):
                    return None

            segs.append(
                {
                    "chrom": chrom,
                    "start": s,
                    "end": e,
                    "cn": num("copyNumber"),
                    "minor": num("minorAlleleCopyNumber"),
                    "baf": num("meanBAF"),
                    "hp1": 0,
                    "hp2": 0,
                    "sites": 0,
                }
            )
    by_chrom = defaultdict(list)
    for i, sg in enumerate(segs):
        by_chrom[sg["chrom"]].append(i)

    def find_seg(chrom, pos):
        for i in by_chrom.get(chrom, []):
            if segs[i]["start"] <= pos < segs[i]["end"]:
                return i
        return None

    # ---- walk the phased germline VCF ----
    n_var = n_phased_het = n_used = 0
    with gzip.open(VCF, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom = f[0]
            if chrom not in AUTOSOMES:
                continue
            n_var += 1
            fmt = f[8].split(":")
            val = f[9].split(":")
            if len(val) < len(fmt):
                continue
            d = dict(zip(fmt, val))
            gt = d.get("GT", "")
            if gt not in ("0|1", "1|0"):
                continue
            n_phased_het += 1
            ad = d.get("AD", "")
            if "," not in ad:
                continue
            try:
                ref_d, alt_d = (int(x) for x in ad.split(",")[:2])
            except ValueError:
                continue
            if ref_d + alt_d < MIN_DP:
                continue
            pos = int(f[1])
            i = find_seg(chrom, pos)
            if i is None:
                continue
            # GT 0|1 -> HP1 carries ref, HP2 carries alt
            if gt == "0|1":
                segs[i]["hp1"] += ref_d
                segs[i]["hp2"] += alt_d
            else:
                segs[i]["hp1"] += alt_d
                segs[i]["hp2"] += ref_d
            segs[i]["sites"] += 1
            n_used += 1

    # ---- score: does haplotype imbalance separate LOH from non-LOH? ----
    rows = []
    for sg in segs:
        if sg["sites"] < MIN_SITES_PER_SEG:
            continue
        tot = sg["hp1"] + sg["hp2"]
        if tot == 0:
            continue
        maj = max(sg["hp1"], sg["hp2"]) / tot
        mid = (sg["start"] + sg["end"]) // 2
        is_loh = in_any(loh, sg["chrom"], mid) is not None
        g = in_any(gain, sg["chrom"], mid)
        l = in_any(loss, sg["chrom"], mid)
        truth_cn = g if g is not None else (l if l is not None else 2.0)
        rows.append(
            {
                "chrom": sg["chrom"],
                "start": sg["start"],
                "end": sg["end"],
                "sites": sg["sites"],
                "hp_major_fraction": round(maj, 4),
                "which_hp_is_major": "HP1" if sg["hp1"] >= sg["hp2"] else "HP2",
                "seqc2_loh": is_loh,
                "seqc2_total_cn": truth_cn,
                "savana_cn": sg["cn"],
                "savana_minor": sg["minor"],
                "savana_baf": sg["baf"],
            }
        )

    lohs = [r["hp_major_fraction"] for r in rows if r["seqc2_loh"]]
    nonl = [r["hp_major_fraction"] for r in rows if not r["seqc2_loh"]]

    def desc(v):
        if not v:
            return None
        s = sorted(v)
        n = len(s)
        return {
            "n": n,
            "median": round(s[n // 2], 4),
            "q1": round(s[n // 4], 4),
            "q3": round(s[(3 * n) // 4], 4),
        }

    # simple separation quality: AUC of hp_major_fraction for LOH vs non-LOH
    auc = None
    if lohs and nonl:
        wins = ties = 0
        for a in lohs:
            for b in nonl:
                if a > b:
                    wins += 1
                elif a == b:
                    ties += 1
        auc = (wins + 0.5 * ties) / (len(lohs) * len(nonl))

    # how well does the derived haplotype ratio agree with SAVANA's own BAF?
    agree = [
        (r["hp_major_fraction"], r["savana_baf"])
        for r in rows
        if r["savana_baf"] is not None
    ]
    corr = None
    if len(agree) > 3:
        xs = [a for a, _ in agree]
        ys = [b for _, b in agree]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        sxy = sum((x - mx) * (y - my) for x, y in agree)
        sxx = sum((x - mx) ** 2 for x in xs)
        syy = sum((y - my) ** 2 for y in ys)
        if sxx > 0 and syy > 0:
            corr = sxy / math.sqrt(sxx * syy)

    out = {
        "generated_by": os.path.basename(__file__),
        "question": "can HP -> allele mapping be derived from the phased germline VCF alone, without returning to the BAM?",
        "method": {
            "source": VCF,
            "rule": "GT 0|1 assigns the ref allelic depth to HP1 and the alt depth to HP2; 1|0 is the reverse; depths are summed over each copy-number segment",
            "min_site_depth": MIN_DP,
            "min_phased_het_sites_per_segment": MIN_SITES_PER_SEG,
        },
        "vcf_scan": {
            "autosomal_variants": n_var,
            "phased_het": n_phased_het,
            "used_after_depth_filter_and_segment_hit": n_used,
        },
        "segments_scored": len(rows),
        "haplotype_major_fraction": {
            "seqc2_loh_segments": desc(lohs),
            "seqc2_non_loh_segments": desc(nonl),
            "separation_auc": round(auc, 4) if auc is not None else None,
        },
        "agreement_with_savana_baf": {
            "n_segments": len(agree),
            "pearson_r": round(corr, 4) if corr is not None else None,
            "note": "SAVANA's meanBAF is computed from the same het sites but without using the phase; a high correlation shows the derivation reproduces it while additionally attributing the major allele to a named haplotype",
        },
        "verdict": None,
        "segments": rows[:200],
    }

    if auc is not None:
        out["verdict"] = (
            "DERIVABLE - haplotype imbalance separates SEQC2 LOH from non-LOH with AUC %.4f; "
            "the blocking input can be produced from the phased germline VCF without touching the BAM"
            % auc
            if auc >= 0.8
            else "WEAK - separation AUC %.4f is not strong enough to trust the derived mapping; investigate before relying on it"
            % auc
        )

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"VCF: autosomal variants={n_var:,} phased-het={n_phased_het:,} used={n_used:,}")
    print(f"segments scored (>= {MIN_SITES_PER_SEG} het sites): {len(rows)}")
    print(f"  LOH segments      : {desc(lohs)}")
    print(f"  non-LOH segments  : {desc(nonl)}")
    print(f"  separation AUC    : {auc}")
    print(f"  corr with SAVANA BAF: {corr}")
    print(f"\nverdict: {out['verdict']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
