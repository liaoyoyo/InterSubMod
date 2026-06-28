#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
全資料集 HP-alignment 分類:對所有 ≥2 ALT genotype 群的區,判斷 subclone 群是否
SAME-HP(同 germline 單倍型,normal-HP cis-control 正交) vs CROSS-HP(跨 HP,cis-control 有效)
vs MIXED(unphased/H3 主導,無法判定)。交叉 CN。決定 matched-normal cis-control 的可用 scope。
輸出:hp_alignment_fullscan.json。compute batch(§13.0,背景跑→Read 驗→另批寫報告)。
"""
import json, os
from collections import defaultdict, Counter
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3

det = {r["region"]: r for r in json.load(open(f"{DATA}/topology_per_region.json"))["detail"]}
targets = [r for r in det.values() if len([k for k in r["populations"] if "A" in k]) >= 2]

vcf_cache = {}
def ref_alt_map(chrom, s, e):
    m = {}
    for src in ("tp", "fp"):
        p = f"{VCFD}/filtered_snv_{src}_{chrom}.vcf.gz"
        if not os.path.exists(p): continue
        try:
            if p not in vcf_cache: vcf_cache[p] = pysam.TabixFile(p)
            for ln in vcf_cache[p].fetch(chrom, max(0, s - 1), e + 1):
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

def hp_norm(a):
    if not a.has_tag("HP"): return 0
    try: return int(str(a.get_tag("HP")).split("-")[0])
    except Exception: return 0

tb = pysam.AlignmentFile(TBAM, "rb")
results = []; n_skip = 0
for r in targets:
    chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2: n_skip += 1; continue
    pop_hp = defaultdict(Counter); pop_n = Counter()
    try:
        for a in tb.fetch(chrom, s, e + 1):
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
            g = geno_str(a, som)
            if "-" in g: continue
            pop_hp[g][hp_norm(a)] += 1; pop_n[g] += 1
    except Exception:
        n_skip += 1; continue
    altpops = [g for g, n in pop_n.most_common() if n >= MINREAD and "A" in g]
    if len(altpops) < 2: n_skip += 1; continue
    def dom(c):
        tot = sum(c.values()) or 1
        if c.get(1,0)/tot >= 0.6: return "H1"
        if c.get(2,0)/tot >= 0.6: return "H2"
        if c.get(3,0)/tot >= 0.6: return "H3"
        return "mixed"
    doms = [dom(pop_hp[g]) for g in altpops]
    pure = [d for d in doms if d in ("H1", "H2")]
    if len(pure) == len(altpops) and len(set(pure)) == 1:
        v = "SAME-HP"
    elif len(set(pure)) >= 2:
        v = "CROSS-HP"
    else:
        v = "MIXED"
    results.append({"region": r["region"], "cn": r["cn"], "n_altpop": len(altpops),
                    "topology_type": r["topology_type"], "hp_alignment": v})

# 交叉表 cn × alignment
ct = defaultdict(Counter)
for x in results: ct[x["cn"]][x["hp_alignment"]] += 1
overall = Counter(x["hp_alignment"] for x in results)
out = {"n_targets": len(targets), "n_classified": len(results), "n_skipped": n_skip,
       "overall": dict(overall),
       "by_cn": {cn: dict(c) for cn, c in ct.items()},
       "note": "CROSS-HP=normal-HP cis-control 有效(subclone 跨 germline HP);SAME-HP=正交(germline-ASM 因共用單倍型抵消,殘差=somatic-cis 不可用 normal 解);MIXED=unphased/H3 無法判定。CN-gain 另受 multiplicity 混淆。",
       "regions": results}
json.dump(out, open(f"{DATA}/hp_alignment_fullscan.json", "w"), ensure_ascii=False)
print("FULLSCAN DONE", json.dumps({"classified": len(results), "skipped": n_skip, "overall": dict(overall)}, ensure_ascii=False))
print("by_cn:")
for cn, c in ct.items(): print(f"  {cn}: {dict(c)}")
