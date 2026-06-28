#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q2: 設下限(min-read)+ ε(2% 雜訊樓地板)是否減少可用/可定位資訊?
對每區重算 genotype 群 at MINREAD∈{1,2,3},數降門檻多救回幾個 ALT 群,
並判其為訊號 vs 雜訊:singleton(n=1) 且與某大群差 1-Hamming = 很可能定序錯誤;
coherent(≥2 read 或 genotype 非單一鄰居) = 可能真稀有群。
ε 效果:同時對 co-occurrence「cell-real」門檻(count>coread×ε)做敏感度。
分塊平行 argv[1]=total argv[2]=idx。輸出 eps_minread_sensitivity[_partN].json。compute batch(§13.0)。
"""
import json, os, sys
from collections import defaultdict, Counter
import pysam

CHUNK_TOTAL = int(sys.argv[1]) if len(sys.argv) > 1 else 1
CHUNK_IDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0
SUF = "" if CHUNK_TOTAL == 1 else f"_part{CHUNK_IDX}"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; READ_CAP = 600
det = json.load(open(f"{DATA}/topology_per_region.json"))["detail"]
items = det[CHUNK_IDX::CHUNK_TOTAL]

_TBX = {}
def ref_alt_map(chrom, s, e):
    m = {}
    for src in ("tp", "fp"):
        p = f"{VCFD}/filtered_snv_{src}_{chrom}.vcf.gz"
        if p not in _TBX: _TBX[p] = pysam.TabixFile(p) if os.path.exists(p) else None
        if _TBX[p] is None: continue
        try:
            for ln in _TBX[p].fetch(chrom, max(0, s - 1), e + 1):
                f = ln.split("\t"); pos = int(f[1])
                if pos not in m: m[pos] = (f[3].upper(), f[4].strip().upper())
        except Exception: pass
    return m

def geno_str(a, som):
    rq = {rr: q for q, rr in a.get_aligned_pairs(matches_only=True)}
    seq = a.query_sequence; g = []
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None or seq is None: g.append("-"); continue
        b = seq[q].upper()
        g.append("A" if b == alt else ("R" if b == ref else "-"))
    return "".join(g)

def hamming1_to_major(g, majors):
    return any(sum(c1 != c2 for c1, c2 in zip(g, mg)) == 1 for mg in majors)

tb = pysam.AlignmentFile(TBAM, "rb")
agg = {"extra_alt_groups_mr1_vs_mr3": 0, "extra_alt_groups_mr2_vs_mr3": 0,
       "extra_singleton": 0, "extra_singleton_1hamming(likely_error)": 0,
       "extra_coherent(possible_rare)": 0, "n_regions": 0}
per_region = []
for r in items:
    chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2: continue
    cnt = Counter()
    nr = 0
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        g = geno_str(a, som)
        if "-" in g: continue
        cnt[g] += 1; nr += 1
        if nr >= READ_CAP: break
    if not cnt: continue
    agg["n_regions"] += 1
    altg = {g: c for g, c in cnt.items() if "A" in g}
    def nalt(mr): return sum(1 for g, c in altg.items() if c >= mr)
    n1, n2, n3 = nalt(1), nalt(2), nalt(3)
    majors = [g for g, c in altg.items() if c >= 3]
    extra = [(g, c) for g, c in altg.items() if c < 3]  # MINREAD=3 會丟的 ALT 群
    sing = [(g, c) for g, c in extra if c == 1]
    sing_err = [(g, c) for g, c in sing if hamming1_to_major(g, majors)]
    coherent = [(g, c) for g, c in extra if c >= 2 or (c == 1 and not hamming1_to_major(g, majors))]
    agg["extra_alt_groups_mr1_vs_mr3"] += (n1 - n3)
    agg["extra_alt_groups_mr2_vs_mr3"] += (n2 - n3)
    agg["extra_singleton"] += len(sing)
    agg["extra_singleton_1hamming(likely_error)"] += len(sing_err)
    agg["extra_coherent(possible_rare)"] += len(coherent)
    if n1 > n3:  # 只記有差的區
        per_region.append({"region": r["region"], "cn": r["cn"], "determinacy": r.get("determinacy", "?"),
                           "nalt_mr1": n1, "nalt_mr3": n3, "n_extra": n1 - n3,
                           "n_singleton_1hamming": len(sing_err), "n_coherent": len(coherent)})

json.dump({"agg": agg, "per_region": per_region}, open(f"{DATA}/eps_minread_sensitivity{SUF}.json", "w"), ensure_ascii=False, indent=1)
print(f"EPS-SENS DONE part={CHUNK_IDX}/{CHUNK_TOTAL}", json.dumps(agg, ensure_ascii=False))
