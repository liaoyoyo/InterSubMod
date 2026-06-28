#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
軸對齊 robustness:在 CN-NEUTRAL 區(兩條 germline HP 都在)測 tumor genotype cluster 是否仍 SAME-HP。
若 neutral 區也 SAME-HP 主導 → 「subclone 在同一 germline HP 內」是普遍性質,非 LOH artifact。
只算 genotype cluster 的 HP 組成(不需 Δβ)。輸出 axis_alignment_neutral.json。compute batch(§13.0)。
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
# 選 CN-neutral + ≥2 ALT 群 + span<8kb + branched/linear(有結構),跨 chr 取 ~12
cand = []
for r in det.values():
    alt = [k for k in r["populations"] if "A" in k]
    if r["cn"] == "neutral" and len(alt) >= 2 and r["span"] < 8000 and ("branched" in r["topology_type"] or "linear" in r["topology_type"]):
        cand.append(r)
cand.sort(key=lambda r: -sum(r["populations"].values()))  # 覆蓋高優先
pilot = cand[:12]

def ref_alt_map(chrom, s, e):
    m = {}
    for src in ("tp", "fp"):
        p = f"{VCFD}/filtered_snv_{src}_{chrom}.vcf.gz"
        if not os.path.exists(p): continue
        try:
            tb = pysam.TabixFile(p)
            for ln in tb.fetch(chrom, max(0, s - 1), e + 1):
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
    v = a.get_tag("HP")
    try: return int(str(v).split("-")[0])
    except Exception: return 0

tb = pysam.AlignmentFile(TBAM, "rb")
align = []
for r in pilot:
    chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2: continue
    pop_hp = defaultdict(Counter); pop_n = Counter()
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        g = geno_str(a, som)
        if "-" in g: continue
        pop_hp[g][hp_norm(a)] += 1; pop_n[g] += 1
    # 只看含 ALT 的群(subclone),REF-only 群是背景
    altpops = [g for g, n in pop_n.most_common() if n >= MINREAD and "A" in g]
    if len(altpops) < 2: continue
    def hp_frac(c):
        tot = sum(c.values()) or 1
        return {"H1": round(c.get(1,0)/tot,2), "H2": round(c.get(2,0)/tot,2), "H3": round(c.get(3,0)/tot,2), "unph": round(c.get(0,0)/tot,2), "n": tot}
    def dom(f): return "H1" if f["H1"]>=0.6 else ("H2" if f["H2"]>=0.6 else ("H3" if f["H3"]>=0.6 else "mixed"))
    doms = {g: dom(hp_frac(pop_hp[g])) for g in altpops}
    pure = [d for d in doms.values() if d in ("H1","H2")]
    if len(set(pure)) == 1 and len(pure) == len(altpops):
        verdict = "SAME-HP"
    elif len(set(pure)) >= 2:
        verdict = "CROSS-HP(subclone 跨 germline HP→normal-HP 部分可 control)"
    else:
        verdict = "MIXED(含 unphased/H3,無法判定)"
    align.append({"region": r["region"], "cn": r["cn"], "n_altpop": len(altpops),
                  "altpop_doms": doms, "altpop_hp": {g: hp_frac(pop_hp[g]) for g in altpops}, "verdict": verdict})

vc = Counter(a["verdict"].split("(")[0] for a in align)
out = {"n_neutral_checked": len(align), "verdict_dist": dict(vc), "regions": align,
       "interpretation": "neutral 區若仍 SAME-HP 主導 → subclone-intra-haplotype 是普遍性質(somatic 突變 cis 於單條親代染色體),非 LOH artifact;normal-HP cis-control 與 subclone 軸正交為結構性。"}
json.dump(out, open(f"{DATA}/axis_alignment_neutral.json", "w"), ensure_ascii=False, indent=1)
print("NEUTRAL AXIS-ALIGN DONE", json.dumps({"n": len(align), "dist": dict(vc)}, ensure_ascii=False))
for a in align:
    print(f"  {a['region']} cn={a['cn']} altpops={a['n_altpop']} doms={a['altpop_doms']} → {a['verdict']}")
