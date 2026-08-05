#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
join_cn_tpfp_markers.py — join CN-count + TP/FP markers onto the WG dashboard cases.

For each dashboard case (chr,pos) it joins (all by genomic-coordinate overlap, real
sources only — §13 no fabrication):
  - tp_fp        : "TP"/"FP" from the kism region_table_{tp,fp}.tsv membership; "TP" default
                   (this dashboard is the somatic-TP confirmation set; 0 cases match FP table)
  - k_cn         : SAVANA-derived CN state count  (minorAlleleCopyNumber>=0.5 -> 2 het / <0.5 -> 1 LOH)
  - savana_cn    : SAVANA total copyNumber (rounded 1 dp) of the overlapping segment
  - savana_minor : SAVANA minorAlleleCopyNumber (rounded 2 dp)
  - seqc2_cn     : SEQC2 truth CN value  (gain_cn.bed / loss_cn.bed integer; 2 if neutral; "2(LOH)" if loh)
  - seqc2_dir    : SEQC2 direction        (gain / loss / loh / neutral)

Writes enriched_cn_tpfp.json = {id: {tp_fp,k_cn,savana_cn,savana_minor,seqc2_cn,seqc2_dir}}
and prints a verification summary (coverage + distributions).
"""
import json, csv, bisect, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DASH = os.path.join(HERE, "20260618_HCC1395_observation_dashboard_FULL.standalone.html")
KISM = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/_assets"
SAVANA = "/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv"
SEQC2DIR = "/big8_disk/data/HCC1395/SEQC2/CNV"
OUT = os.path.join(HERE, "enriched_cn_tpfp.json")


def load_cases():
    s = open(DASH, encoding="utf-8").read()
    tag = '<script id="d" type="application/json">'
    a = s.find(tag) + len(tag)
    b = s.find("</script>", a)
    return json.loads(s[a:b])["cases"]


def load_tpfp():
    mem = {}
    for lab, f in (("TP", "region_table_tp.tsv"), ("FP", "region_table_fp.tsv")):
        with open(os.path.join(KISM, f)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                mem[(r["chr"], int(r["pos"]))] = lab
    return mem


def build_intervals(rows):
    """rows: list of (chr,start,end,payload). -> {chr: (sorted_starts, [(end,payload)])}"""
    by = {}
    for c, s, e, p in rows:
        by.setdefault(c, []).append((s, e, p))
    idx = {}
    for c, lst in by.items():
        lst.sort()
        idx[c] = ([x[0] for x in lst], lst)
    return idx


def lookup(idx, chrom, pos):
    if chrom not in idx:
        return None
    starts, lst = idx[chrom]
    i = bisect.bisect_right(starts, pos) - 1   # last interval with start <= pos
    steps = 0
    while i >= 0 and steps < 500:               # scan back over possibly-overlapping intervals
        s, e, p = lst[i]
        if s <= pos <= e:
            return p
        if pos - s > 50_000_000:                # no SEQC2/SAVANA segment is this long -> stop
            break
        i -= 1; steps += 1
    return None


def load_savana():
    rows = []
    with open(SAVANA) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                cn = float(r["copyNumber"]); minor = float(r["minorAlleleCopyNumber"])
            except (ValueError, KeyError):
                continue
            rows.append((r["chromosome"], int(r["start"]), int(r["end"]), (cn, minor)))
    return build_intervals(rows)


def load_seqc2():
    def bed(path, conv):
        out = []
        with open(path) as fh:
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                out.append((p[0], int(p[1]), int(p[2]), conv(p[3])))
        return out
    def cn_num(v):
        f = float(v)
        return int(f) if f == int(f) else round(f, 1)
    gain = bed(os.path.join(SEQC2DIR, "ngs_benchmark_cnv_gain_cn.bed"), lambda v: ("gain", cn_num(v)))
    loss = bed(os.path.join(SEQC2DIR, "ngs_benchmark_cnv_loss_cn.bed"), lambda v: ("loss", cn_num(v)))
    loh = bed(os.path.join(SEQC2DIR, "ngs_benchmark_cnvs_gain_loss_loh.bed"), lambda v: (v, None))
    # priority: explicit gain/loss CN beds first, then loh from the direction bed
    return build_intervals(gain + loss), build_intervals([r for r in loh if r[3][0] == "loh"])


def main():
    cases = load_cases()
    tpfp = load_tpfp()
    sav = load_savana()
    seqc2_gl, seqc2_loh = load_seqc2()

    enr = {}
    cov = Counter(); kcn = Counter(); s2dir = Counter(); tpc = Counter()
    for c in cases:
        chrom, pos = c["chr"], int(c["pos"])
        rec = {}
        # tp/fp
        rec["tp_fp"] = tpfp.get((chrom, pos), "TP")
        tpc[rec["tp_fp"] + ("" if (chrom, pos) in tpfp else "(default)")] += 1
        # savana
        sv = lookup(sav, chrom, pos)
        if sv:
            cn, minor = sv
            rec["savana_cn"] = round(cn, 1)
            rec["savana_minor"] = round(minor, 2)
            rec["k_cn"] = 2 if minor >= 0.5 else 1
            cov["savana"] += 1
            kcn[rec["k_cn"]] += 1
        else:
            rec["savana_cn"] = None; rec["savana_minor"] = None; rec["k_cn"] = None
        # seqc2
        gl = lookup(seqc2_gl, chrom, pos)
        if gl:
            direction, cnv = gl
            rec["seqc2_dir"] = direction
            rec["seqc2_cn"] = cnv
            cov["seqc2_glcn"] += 1
        elif lookup(seqc2_loh, chrom, pos):
            rec["seqc2_dir"] = "loh"; rec["seqc2_cn"] = "2(LOH)"
            cov["seqc2_loh"] += 1
        else:
            rec["seqc2_dir"] = "neutral"; rec["seqc2_cn"] = 2
            cov["seqc2_neutral"] += 1
        s2dir[rec["seqc2_dir"]] += 1
        enr[c["id"]] = rec

    json.dump(enr, open(OUT, "w"), ensure_ascii=False)
    n = len(cases)
    print(f"cases: {n}  -> wrote {OUT} ({os.path.getsize(OUT)//1024} KB)")
    print(f"tp_fp: {dict(tpc)}")
    print(f"SAVANA coverage: {cov['savana']}/{n} ({100*cov['savana']/n:.1f}%)  k_cn dist: {dict(kcn)}")
    print(f"SEQC2 dir dist: {dict(s2dir)}  (gl-cn {cov['seqc2_glcn']} / loh {cov['seqc2_loh']} / neutral {cov['seqc2_neutral']})")
    # spot-check 3 known loci
    for cid in ("chr20_39492762", "chr1_207491367", "chr8_81434745"):
        print(f"  {cid}: {enr.get(cid,'NOT IN DASHBOARD')}")


if __name__ == "__main__":
    main()
