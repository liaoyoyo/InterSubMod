#!/usr/bin/env python3
"""PRE-REGISTERED pilot: does HP-resolved depth on the TAGGED TUMOUR bam
recover allelic imbalance?

The previous attempt failed because it used the phased germline VCF's AD field,
which comes from the matched normal and is ~50:50 by construction (LOH
separation AUC 0.4341).  That run is retained as the NEGATIVE CONTROL: if the
tumour bam shows separation where the normal did not, the signal is tumour-
derived rather than an artefact of the site set or the segment frame.

DECISION RULES ARE FIXED BEFORE EXECUTION.  Each hypothesis below carries an
explicit support threshold and an explicit refute threshold, with the band
between them declared INCONCLUSIVE rather than being argued either way after
the fact.

  H1  LOH vs non-LOH separation by haplotype major fraction
        support  AUC >= 0.80      refute  AUC < 0.65
  H2  agreement with SAVANA meanBAF
        support  r   >= 0.70      refute  r   < 0.40
  H3  SEQC2 LOH segments show strong imbalance
        support  median >= 0.85   refute  median < 0.70
  H4  SEQC2 copy-neutral non-LOH segments stay near balance  (negative control)
        support  0.45 <= median <= 0.62   refute  median > 0.70
  H5  haplotype-tagged depth tracks total copy number
        support  r   >= 0.40      refute  r   < 0.20

Scope is a deliberate pilot: three chromosomes chosen for contrast --
chr8 (previously characterised as LOH-dominated), chr21 (the cleanest
chromosome, 91.23% CN-neutral units) and chr2 (mixed, 19.75% neutral).
That gives a positive, a negative and an intermediate arm in one run.
"""

from __future__ import annotations

import json
import math
import os
import random
from collections import Counter, defaultdict

import pysam

# ---------------------------------------------------------------- pre-registration
PREREG = {
    "H1_loh_separation": {"metric": "AUC", "support": 0.80, "refute": 0.65, "direction": "higher"},
    "H2_agreement_with_savana_baf": {"metric": "pearson_r", "support": 0.70, "refute": 0.40, "direction": "higher"},
    "H3_loh_imbalance": {"metric": "median_major_fraction_LOH", "support": 0.85, "refute": 0.70, "direction": "higher"},
    "H4_neutral_balance": {"metric": "median_major_fraction_neutral", "support": [0.45, 0.62], "refute": 0.70, "direction": "band"},
    "H5_depth_tracks_cn": {"metric": "pearson_r", "support": 0.40, "refute": 0.20, "direction": "higher"},
}
NEGATIVE_CONTROL = {
    "source": "hp_to_allele_pilot.json",
    "what": "same site set and segment frame scored on the matched-normal AD field",
    "result_auc": 0.4341,
    "role": "if the tumour bam separates and the normal does not, the signal is tumour-derived",
}

BAM = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
    "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
)
VCF = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
    "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/"
    "germline_phased_merged.vcf.gz"
)
SEQC2 = "/big8_disk/data/HCC1395/SEQC2/CNV"
SAVANA = (
    "/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna/"
    "HCC1395_segmented_absolute_copy_number.tsv"
)
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(OUT_DIR, "hp_pileup_pilot_prereg.json")

CHROMS = ["chr8", "chr21", "chr2"]
MAX_SEG_PER_CHROM = 15
MAX_SITES_PER_SEG = 100
MIN_HP_READS = 10
SEED = 20260813


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


def hit(segs, chrom, pos):
    for s, e, v in segs.get(chrom, []):
        if s <= pos < e:
            return v
    return None


def pearson(pairs):
    if len(pairs) < 3:
        return None
    xs = [a for a, _ in pairs]
    ys = [b for _, b in pairs]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxy = sum((x - mx) * (y - my) for x, y in pairs)
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def auc_of(pos_vals, neg_vals):
    if not pos_vals or not neg_vals:
        return None
    w = t = 0
    for a in pos_vals:
        for b in neg_vals:
            if a > b:
                w += 1
            elif a == b:
                t += 1
    return (w + 0.5 * t) / (len(pos_vals) * len(neg_vals))


def desc(v):
    if not v:
        return None
    s = sorted(v)
    n = len(s)
    return {"n": n, "median": round(s[n // 2], 4), "q1": round(s[n // 4], 4),
            "q3": round(s[(3 * n) // 4], 4)}


def verdict(name, value):
    spec = PREREG[name]
    if value is None:
        return "NO_DATA"
    if spec["direction"] == "band":
        lo, hi = spec["support"]
        if lo <= value <= hi:
            return "SUPPORTED"
        if value > spec["refute"]:
            return "REFUTED"
        return "INCONCLUSIVE"
    if value >= spec["support"]:
        return "SUPPORTED"
    if value < spec["refute"]:
        return "REFUTED"
    return "INCONCLUSIVE"


def main():
    random.seed(SEED)
    os.makedirs(OUT_DIR, exist_ok=True)

    cat = load_bed(f"{SEQC2}/ngs_benchmark_cnvs_gain_loss_loh.bed", cast=str)
    loh = defaultdict(list)
    for c, items in cat.items():
        for s, e, v in items:
            if v == "loh":
                loh[c].append((s, e, v))
    gain = load_bed(f"{SEQC2}/ngs_benchmark_cnv_gain_cn.bed")
    loss = load_bed(f"{SEQC2}/ngs_benchmark_cnv_loss_cn.bed")

    # SAVANA segments restricted to the pilot chromosomes
    segs = []
    with open(SAVANA) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if p[col["chromosome"]] not in CHROMS:
                continue
            try:
                s, e = int(p[col["start"]]), int(p[col["end"]])
                baf = float(p[col["meanBAF"]])
                nhet = int(p[col["no_hetSNPs"]])
            except (ValueError, IndexError):
                continue
            try:
                cnv = float(p[col["copyNumber"]])
            except (ValueError, IndexError):
                cnv = None
            if nhet < 50:
                continue
            segs.append({"chrom": p[col["chromosome"]], "start": s, "end": e,
                         "savana_baf": baf, "savana_cn": cnv})

    by_chrom = defaultdict(list)
    for sg in segs:
        by_chrom[sg["chrom"]].append(sg)
    chosen = []
    for c in CHROMS:
        pool = by_chrom.get(c, [])
        random.shuffle(pool)
        chosen.extend(pool[:MAX_SEG_PER_CHROM])

    vcf = pysam.VariantFile(VCF)
    bam = pysam.AlignmentFile(BAM, "rb")

    rows = []
    sites_done = 0
    for sg in chosen:
        sites = []
        try:
            for rec in vcf.fetch(sg["chrom"], sg["start"], sg["end"]):
                s = rec.samples[0]
                gt = s.get("GT")
                if gt is None or not s.phased:
                    continue
                if set(gt) != {0, 1}:
                    continue
                sites.append((rec.pos, gt))
        except ValueError:
            continue
        if len(sites) < 20:
            continue
        random.shuffle(sites)
        sites = sites[:MAX_SITES_PER_SEG]

        hp1 = hp2 = 0
        used = 0
        for pos, gt in sites:
            c1 = c2 = 0
            try:
                for col_ in bam.pileup(sg["chrom"], pos - 1, pos,
                                       truncate=True, min_base_quality=0,
                                       max_depth=500):
                    for pr in col_.pileups:
                        r = pr.alignment
                        if not r.has_tag("HP"):
                            continue
                        fam = str(r.get_tag("HP"))[0]
                        if fam == "1":
                            c1 += 1
                        elif fam == "2":
                            c2 += 1
            except (ValueError, OSError):
                continue
            if c1 + c2 >= MIN_HP_READS:
                hp1 += c1
                hp2 += c2
                used += 1
            sites_done += 1

        tot = hp1 + hp2
        if used < 10 or tot == 0:
            continue
        mid = (sg["start"] + sg["end"]) // 2
        is_loh = hit(loh, sg["chrom"], mid) is not None
        g = hit(gain, sg["chrom"], mid)
        l = hit(loss, sg["chrom"], mid)
        truth_cn = g if g is not None else (l if l is not None else 2.0)
        rows.append({
            "chrom": sg["chrom"], "start": sg["start"], "end": sg["end"],
            "sites_used": used, "hp1_reads": hp1, "hp2_reads": hp2,
            "hp_major_fraction": round(max(hp1, hp2) / tot, 4),
            "major_is": "HP1" if hp1 >= hp2 else "HP2",
            "mean_hp_depth_per_site": round(tot / used, 2),
            "seqc2_loh": is_loh, "seqc2_total_cn": truth_cn,
            "savana_baf": sg["savana_baf"], "savana_cn": sg["savana_cn"],
        })

    bam.close()
    vcf.close()

    lohv = [r["hp_major_fraction"] for r in rows if r["seqc2_loh"]]
    nonv = [r["hp_major_fraction"] for r in rows if not r["seqc2_loh"]]
    neutral = [r["hp_major_fraction"] for r in rows
               if not r["seqc2_loh"] and abs(r["seqc2_total_cn"] - 2.0) < 0.01]

    a1 = auc_of(lohv, nonv)
    r2 = pearson([(r["hp_major_fraction"], r["savana_baf"]) for r in rows])
    m3 = desc(lohv)["median"] if lohv else None
    m4 = desc(neutral)["median"] if neutral else None
    r5 = pearson([(r["mean_hp_depth_per_site"], r["seqc2_total_cn"]) for r in rows])

    results = {
        "H1_loh_separation": {"value": round(a1, 4) if a1 is not None else None,
                              "verdict": verdict("H1_loh_separation", a1)},
        "H2_agreement_with_savana_baf": {"value": round(r2, 4) if r2 is not None else None,
                                         "verdict": verdict("H2_agreement_with_savana_baf", r2)},
        "H3_loh_imbalance": {"value": m3, "verdict": verdict("H3_loh_imbalance", m3)},
        "H4_neutral_balance": {"value": m4, "verdict": verdict("H4_neutral_balance", m4)},
        "H5_depth_tracks_cn": {"value": round(r5, 4) if r5 is not None else None,
                               "verdict": verdict("H5_depth_tracks_cn", r5)},
    }

    supported = sum(1 for v in results.values() if v["verdict"] == "SUPPORTED")
    refuted = sum(1 for v in results.values() if v["verdict"] == "REFUTED")
    inconcl = sum(1 for v in results.values() if v["verdict"] == "INCONCLUSIVE")

    out = {
        "generated_by": os.path.basename(__file__),
        "design": "PRE-REGISTERED pilot; thresholds fixed before execution",
        "task_type": "A_exploratory_pilot",
        "partial_scope_flag": "PILOT - 3 of 22 autosomes (%s), max %d segments each, max %d sites per segment"
        % (", ".join(CHROMS), MAX_SEG_PER_CHROM, MAX_SITES_PER_SEG),
        "chromosome_choice_rationale": "chr8 LOH-dominated (positive arm), chr21 cleanest at 91.23% CN-neutral (negative arm), chr2 mixed at 19.75% neutral (intermediate arm)",
        "prereg_thresholds": PREREG,
        "negative_control": NEGATIVE_CONTROL,
        "inputs": {"tagged_tumour_bam": BAM, "phased_germline_vcf": VCF,
                   "seed": SEED, "min_hp_reads_per_site": MIN_HP_READS},
        "counts": {"segments_scored": len(rows), "sites_pileup_attempted": sites_done,
                   "loh_segments": len(lohv), "non_loh_segments": len(nonv),
                   "cn_neutral_non_loh_segments": len(neutral)},
        "distributions": {"loh": desc(lohv), "non_loh": desc(nonv),
                          "cn_neutral_non_loh": desc(neutral)},
        "results": results,
        "scoreboard": {"supported": supported, "refuted": refuted,
                       "inconclusive": inconcl, "total": len(results)},
        "segments": rows,
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"segments scored: {len(rows)}  (LOH {len(lohv)} / non-LOH {len(nonv)} / CN-neutral {len(neutral)})")
    print(f"sites pileup attempted: {sites_done}")
    print(f"  LOH        : {desc(lohv)}")
    print(f"  non-LOH    : {desc(nonv)}")
    print(f"  CN-neutral : {desc(neutral)}")
    print("\n--- pre-registered verdicts ---")
    for k, v in results.items():
        sp = PREREG[k]
        thr = (f"support {sp['support']} / refute {sp['refute']}")
        print(f"  {k:32s} value={str(v['value']):>8s}  [{thr}]  -> {v['verdict']}")
    print(f"\nscoreboard: SUPPORTED {supported} / REFUTED {refuted} / INCONCLUSIVE {inconcl}")
    print(f"negative control (matched normal): AUC {NEGATIVE_CONTROL['result_auc']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
