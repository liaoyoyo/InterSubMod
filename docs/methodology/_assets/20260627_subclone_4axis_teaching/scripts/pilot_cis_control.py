#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
T-GATE-GB cis-control PILOT(5-10 區,驗證設計 + 觀察分布,先不設門檻)。
每區:tumor 依 genotype 分群 → 取 top-2 群的 per-CpG Δβ(subclone 候選訊號);
      normal 依 HP 分群(H1 vs H2)→ 同 CpG 的 per-CpG Δβ(germline-ASM 基線)。
觀察每 CpG (tumor_dbeta, normal_dbeta):tumor 高+normal 高=germline-ASM(非subclone);
      tumor 高+normal≈0=候選 subclone-specific。先輸出分布,不裁決。
輸出:pilot_cis_control.json。BAM compute → 落 JSON(§13.0,跑完才另批讀+寫報告)。
"""
import json, os
from collections import defaultdict, Counter
import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3

# 選 pilot 區:undefined/needs_methyl + ≥2 ALT 群 + 乾淨 CN + span<8kb,跨不同 chr 取前 8
cs = json.load(open(f"{DATA}/candidate_scoring.json"))
det = {r["region"]: r for r in json.load(open(f"{DATA}/topology_per_region.json"))["detail"]}
cand = []
seen_chr = set()
for q in cs["queue"]:
    r = det.get(q["region"])
    if not r: continue
    alt = [k for k in r["populations"] if "A" in k]
    if (r.get("undefined") or q["needs_methyl"]) and len(alt) >= 2 and r["cn"] in ("loh", "neutral") and r["span"] < 8000:
        cand.append(r)
# 取跨不同 chr 前 8(多樣)
pilot = []
for r in cand:
    if r["chrom"] not in seen_chr or len(pilot) < 8:
        pilot.append(r); seen_chr.add(r["chrom"])
    if len(pilot) >= 8: break

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
        except Exception:
            pass
    return m

def read_meth(a):
    try: mb = a.modified_bases
    except Exception: return None
    qr = {q: r for q, r in a.get_aligned_pairs(matches_only=True)}
    meth = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != "m": continue
            for qpos, mlq in lst:
                r = qr.get(qpos)
                if r is not None: meth[r] = mlq / 255.0
    return meth

def geno_str(a, som):
    rq = {r: q for q, r in a.get_aligned_pairs(matches_only=True)}
    seq = a.query_sequence; g = []
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None or seq is None: g.append("-"); continue
        b = seq[q].upper()
        g.append("A" if b == alt else ("R" if b == ref else "-"))
    return "".join(g)

def per_cpg(read_meth_list, idset):
    cpg = defaultdict(list)
    for rn in idset:
        for c, b in read_meth_list[rn].items(): cpg[c].append(b)
    return {c: np.mean(v) for c, v in cpg.items() if len(v) >= MINREAD}

tb = pysam.AlignmentFile(TBAM, "rb"); nb = pysam.AlignmentFile(NBAM, "rb")
out_regions = []; all_pairs = []
for r in pilot:
    chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2: continue
    # tumor: meth + geno + cluster
    tmeth = {}; tgeno = {}
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        m = read_meth(a);
        if m is None: continue
        g = geno_str(a, som)
        if "-" in g: continue
        tmeth[a.query_name] = m; tgeno[a.query_name] = g
    pop = defaultdict(list)
    for rn, g in tgeno.items(): pop[g].append(rn)
    pops = sorted([(g, ids) for g, ids in pop.items() if len(ids) >= MINREAD], key=lambda x: -len(x[1]))
    # normal: meth + HP
    nmeth = {}; nhp = {}
    for a in nb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        if not a.has_tag("HP"): continue
        m = read_meth(a)
        if m is None: continue
        nmeth[a.query_name] = m; nhp[a.query_name] = a.get_tag("HP")
    h1 = [rn for rn, h in nhp.items() if h == 1]; h2 = [rn for rn, h in nhp.items() if h == 2]
    rec = {"region": r["region"], "n_som": len(som), "cn": r["cn"], "haplotypes": r["haplotypes"],
           "tumor_n_pop": len(pops), "tumor_popA": pops[0][0] if pops else None, "tumor_popA_n": len(pops[0][1]) if pops else 0,
           "tumor_popB": pops[1][0] if len(pops) >= 2 else None, "tumor_popB_n": len(pops[1][1]) if len(pops) >= 2 else 0,
           "normal_h1_n": len(h1), "normal_h2_n": len(h2), "cpg_pairs": [], "note": ""}
    if len(pops) < 2:
        rec["note"] = "tumor <2 群(覆蓋全 som 的 read 不足)"; out_regions.append(rec); continue
    if len(h1) < MINREAD or len(h2) < MINREAD:
        rec["note"] = "normal HP1/HP2 reads 不足→無 germline 基線"; out_regions.append(rec); continue
    tA = per_cpg(tmeth, set(pops[0][1])); tB = per_cpg(tmeth, set(pops[1][1]))
    nH1 = per_cpg(nmeth, set(h1)); nH2 = per_cpg(nmeth, set(h2))
    for c in set(tA) & set(tB) & set(nH1) & set(nH2):
        td = round(tA[c] - tB[c], 3); nd = round(nH1[c] - nH2[c], 3)
        rec["cpg_pairs"].append({"cpg": c, "tumor_dbeta": td, "normal_dbeta": nd})
        all_pairs.append({"region": r["region"], "cpg": c, "tumor_dbeta": td, "normal_dbeta": nd})
    rec["n_cpg_both"] = len(rec["cpg_pairs"])
    out_regions.append(rec)

# 觀察分布(先不設門檻)
tg = [abs(p["tumor_dbeta"]) for p in all_pairs]; ng = [abs(p["normal_dbeta"]) for p in all_pairs]
summary = {"n_pilot_regions": len(pilot), "n_regions_with_cpg": sum(1 for r in out_regions if r.get("n_cpg_both")),
           "n_cpg_pairs_total": len(all_pairs),
           "tumor_dbeta_abs": {"median": round(float(np.median(tg)), 3) if tg else None, "p90": round(float(np.percentile(tg, 90)), 3) if tg else None, "max": round(max(tg), 3) if tg else None},
           "normal_dbeta_abs": {"median": round(float(np.median(ng)), 3) if ng else None, "p90": round(float(np.percentile(ng, 90)), 3) if ng else None, "max": round(max(ng), 3) if ng else None},
           "observation": "tumor 高+normal 高=germline-ASM(非subclone);tumor 高+normal≈0=候選 subclone-specific。門檻待看分布散點後定。"}
json.dump({"summary": summary, "regions": out_regions, "all_cpg_pairs": all_pairs},
          open(f"{DATA}/pilot_cis_control.json", "w"), ensure_ascii=False, indent=1)
print("PILOT DONE", json.dumps(summary, ensure_ascii=False))
for r in out_regions:
    print(f"  {r['region']} som={r['n_som']} tumor群={r['tumor_n_pop']}(A{r['tumor_popA_n']}/B{r['tumor_popB_n']}) normalHP1/2={r['normal_h1_n']}/{r['normal_h2_n']} cpg_both={r.get('n_cpg_both',0)} {r['note']}")
