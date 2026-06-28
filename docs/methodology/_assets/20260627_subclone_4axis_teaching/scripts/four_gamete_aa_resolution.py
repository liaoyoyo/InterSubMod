#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4-gamete(RR/RA/AR/AA 全有=成環)時:AA 屬 AR 還是 RA 譜系?
用戶提議:① same-HP 前提(AR/RA/AA 須同 germline HP,否則跨HP=獨立非譜系問題)
         ② VAF/count:AA 的母系應 CCF≥AA;count_AR vs count_RA 大者=較可能母系
         ③ 甲基 read 距離:AA mean-β 距 AR vs RA,近者=母系
測 VAF 與甲基是否一致。從 BAM 重 derive 全基因型(非截斷)。
targets = incompatible + drop_noise 區(4-gamete 候選)。輸出 four_gamete_aa_resolution.json。compute batch(§13.0)。
"""
import json, os
from itertools import combinations
from collections import defaultdict, Counter
import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3

det = json.load(open(f"{DATA}/topology_per_region.json"))["detail"]
targets = [r for r in det if r["determinacy"] == "incompatible" or r.get("drop_noise_frac", 0) > 0]

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

def read_meth(a):
    try: mb = a.modified_bases
    except Exception: return None
    qr = {q: rr for q, rr in a.get_aligned_pairs(matches_only=True)}
    m = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != "m": continue
            for qp, mlq in lst:
                rr = qr.get(qp)
                if rr is not None: m[rr] = mlq / 255.0
    return m

def hp_norm(a):
    if not a.has_tag("HP"): return 0
    try: return int(str(a.get_tag("HP")).split("-")[0])
    except Exception: return 0

READ_CAP = 400; SOM_CAP = 30

tb = pysam.AlignmentFile(TBAM, "rb")
cases = []
for r in targets:
    chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som_all = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som_all) < 2: continue
    # 每 read 一次 aligned_pairs → 各 som 的 base(O(reads×som),非 ×aligned_pairs)
    reads = []; nr = 0
    cov = Counter()
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        rq = {rr: q for q, rr in a.get_aligned_pairs(matches_only=True)}
        seq = a.query_sequence
        if seq is None: continue
        sb = {}
        for pos, ref, alt in som_all:
            q = rq.get(pos - 1)
            if q is None: sb[pos] = "-"; continue
            b = seq[q].upper(); sb[pos] = "A" if b == alt else ("R" if b == ref else "-")
            if sb[pos] in "RA": cov[pos] += 1
        m = read_meth(a)
        reads.append({"hp": hp_norm(a), "meth": (float(np.mean(list(m.values()))) if m and len(m) >= 3 else None), "sb": sb})
        nr += 1
        if nr >= READ_CAP: break
    # 只取覆蓋最高的 SOM_CAP 位點(bound pair 數)
    som = [t for t in som_all if t[0] in dict(cov.most_common(SOM_CAP))]
    for (pi, ri, ai), (pj, rj, aj) in combinations(som, 2):
        groups = defaultdict(list)
        for rd in reads:
            gi, gj = rd["sb"][pi], rd["sb"][pj]
            if gi in "RA" and gj in "RA": groups[gi + gj].append(rd)
        if not all(len(groups[g]) >= MINREAD for g in ("RR", "RA", "AR", "AA")): continue
        # same-HP 前提
        def domhp(rds):
            c = Counter(rd["hp"] for rd in rds); tot = sum(c.values()) or 1
            for k in (1, 2, 3):
                if c.get(k, 0) / tot >= 0.6: return k
            return 0
        hp_ar, hp_ra, hp_aa = domhp(groups["AR"]), domhp(groups["RA"]), domhp(groups["AA"])
        same_hp = (hp_ar in (1, 2) and hp_ar == hp_ra == hp_aa)
        # VAF/count 判母系(母系 count 大)
        n_ar, n_ra, n_aa = len(groups["AR"]), len(groups["RA"]), len(groups["AA"])
        vaf_parent = "AR" if n_ar > n_ra else ("RA" if n_ra > n_ar else "tie")
        # 甲基距離判母系
        def mbeta(rds):
            v = [rd["meth"] for rd in rds if rd["meth"] is not None]
            return float(np.mean(v)) if len(v) >= MINREAD else None
        b_ar, b_ra, b_aa = mbeta(groups["AR"]), mbeta(groups["RA"]), mbeta(groups["AA"])
        meth_parent = None
        if None not in (b_ar, b_ra, b_aa):
            meth_parent = "AR" if abs(b_aa - b_ar) < abs(b_aa - b_ra) else "RA"
        cases.append({"region": r["region"], "cn": r["cn"], "pair": [pi, pj],
                      "n_RR": len(groups["RR"]), "n_AR": n_ar, "n_RA": n_ra, "n_AA": n_aa,
                      "hp_AR": hp_ar, "hp_RA": hp_ra, "hp_AA": hp_aa, "same_hp": same_hp,
                      "vaf_parent": vaf_parent, "beta_AR": round(b_ar, 3) if b_ar else None,
                      "beta_RA": round(b_ra, 3) if b_ra else None, "beta_AA": round(b_aa, 3) if b_aa else None,
                      "meth_parent": meth_parent,
                      "vaf_meth_agree": (meth_parent is not None and vaf_parent == meth_parent)})

# 摘要
same = [c for c in cases if c["same_hp"]]
clean = [c for c in same if c["cn"] in ("neutral", "loh")]
both = [c for c in clean if c["meth_parent"] is not None and c["vaf_parent"] != "tie"]
agree = [c for c in both if c["vaf_meth_agree"]]
summary = {
    "n_targets": len(targets), "n_4gamete_pairs": len(cases),
    "same_hp(前提成立)": len(same), "cross_or_mixed_hp(前提不成立→跨HP獨立)": len(cases) - len(same),
    "same_hp_AND_cn_clean": len(clean),
    "resolvable(same-HP+clean+VAF非tie+甲基可測)": len(both),
    "vaf_meth_agree": len(agree), "vaf_meth_disagree": len(both) - len(agree),
    "verdict": "same-HP 前提必查(跨HP=兩突變在不同染色體=AA 可能 doublet/CN artifact 非譜系後代)。"
               "VAF(count 大=母系)可判但 CN-gain 區 multiplicity 混淆;甲基距離=L3 弱輔助。"
               "VAF 與甲基一致才較可信;不一致或 CN-gain→不可靠。"}
json.dump({"summary": summary, "cases": cases}, open(f"{DATA}/four_gamete_aa_resolution.json", "w"), ensure_ascii=False, indent=1)
print("4-GAMETE AA RESOLUTION DONE")
print(json.dumps(summary, ensure_ascii=False, indent=1))
for c in cases[:8]:
    print(f"  {c['region']} pair={c['pair']} RR/AR/RA/AA={c['n_RR']}/{c['n_AR']}/{c['n_RA']}/{c['n_AA']} same_hp={c['same_hp']} VAF→{c['vaf_parent']} 甲基→{c['meth_parent']} agree={c['vaf_meth_agree']}")
