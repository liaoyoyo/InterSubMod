#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q6 延伸:對 H3-phasing 無法定相的區(HP1≈HP2 無 germline-ASM 45 + 缺一HP 38),
用甲基 read 距離確認「是否就是一群」(負向篩選)。
測:① ALT(somatic) 是否集中單一 HP(若是→分群基因型已定,定相 moot)
    ② 甲基 per-read β 是否雙峰(GMM,扣不掉就用負向:unimodal=無次結構=支持同一群)
③ ALT reads 本身是否與 REF reads 甲基可分(mutation 的 read 群聚是否合理 corroborate)
輸出 h3_unresolved_grouping.json。compute batch(§13.0)。
"""
import json, os
from collections import defaultdict, Counter
import numpy as np
import pysam
from sklearn.mixture import GaussianMixture

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINR = 4; READ_CAP = 600

h3 = json.load(open(f"{DATA}/h3_methyl_phasing.json"))["regions"]
# 取「無法定相」的區:HP不可分 + 資料不足(缺HP)
unresolved = [r["region"] for r in h3 if ("不可分" in r["verdict"]) or ("不足" in r["verdict"]) or ("只有 HP" in r["verdict"])]
det = {r["region"]: r for r in json.load(open(f"{DATA}/topology_per_region.json"))["detail"]}

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

def geno_has_alt(a, som):
    rq = {rr: q for q, rr in a.get_aligned_pairs(matches_only=True)}
    seq = a.query_sequence
    if seq is None: return None
    has = False; cov = 0
    for pos, ref, alt in som:
        q = rq.get(pos - 1)
        if q is None: continue
        cov += 1
        if seq[q].upper() == alt: has = True
    return has if cov > 0 else None

def hp_norm(a):
    if not a.has_tag("HP"): return 0
    try: return int(str(a.get_tag("HP")).split("-")[0])
    except Exception: return 0

def bimodal(vals):
    x = np.array(vals).reshape(-1, 1)
    if len(x) < 12: return None
    try:
        g1 = GaussianMixture(1, random_state=0).fit(x); g2 = GaussianMixture(2, n_init=1, random_state=0).fit(x)
    except Exception: return None
    if g2.bic(x) >= g1.bic(x): return False
    m = g2.means_.flatten(); lab = g2.predict(x)
    return bool(abs(m[0] - m[1]) >= 0.2 and min((lab == 0).sum(), (lab == 1).sum()) >= 4)

tb = pysam.AlignmentFile(TBAM, "rb")
out = []
for region in unresolved:
    r = det.get(region)
    if not r: continue
    chrom, s = r["chrom"], r["start"]; e = int(region.split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if not som: continue
    alt_beta = []; ref_beta = []; all_beta = []; alt_hp = Counter()
    nr = 0
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        m = read_meth(a)
        if not m or len(m) < 3: continue
        mb = float(np.mean(list(m.values())))
        ha = geno_has_alt(a, som)
        if ha is None: continue
        all_beta.append(mb)
        if ha: alt_beta.append(mb); alt_hp[hp_norm(a)] += 1
        else: ref_beta.append(mb)
        nr += 1
        if nr >= READ_CAP: break
    rec = {"region": region, "cn": r["cn"], "n_alt_reads": len(alt_beta), "n_ref_reads": len(ref_beta),
           "alt_hp_dist": dict(alt_hp)}
    # ① ALT 是否集中單一 HP
    phased = {k: v for k, v in alt_hp.items() if k in (1, 2)}
    if phased:
        tot = sum(phased.values())
        domfrac = max(phased.values()) / tot if tot else 0
        rec["alt_single_hp"] = bool(domfrac >= 0.9)
        rec["alt_dom_hp_frac"] = round(float(domfrac), 2)
    # ② 甲基次分群(ALT reads)
    bi_alt = bimodal(alt_beta) if len(alt_beta) >= 12 else None
    rec["alt_methyl_bimodal"] = bi_alt
    # ③ ALT vs REF 甲基可分?(mutation read 群聚 corroborate)
    if len(alt_beta) >= MINR and len(ref_beta) >= MINR:
        da = float(abs(np.mean(alt_beta) - np.mean(ref_beta)))
        rec["alt_vs_ref_dbeta"] = round(da, 3)
        rec["mut_methyl_coherent"] = bool(da >= 0.2)  # ALT 群與 REF 群甲基明顯不同
    out.append(rec)

# 摘要
def cnt(pred): return sum(1 for x in out if pred(x))
summary = {
    "n_unresolved_input": len(unresolved), "n_analyzed": len(out),
    "alt_single_hp(ALT集中單一HP→分群基因型已定)": cnt(lambda x: x.get("alt_single_hp") is True),
    "alt_cross_hp(ALT跨HP)": cnt(lambda x: x.get("alt_single_hp") is False),
    "alt_methyl_unimodal(無次結構→支持同一群)": cnt(lambda x: x.get("alt_methyl_bimodal") is False),
    "alt_methyl_bimodal(疑次結構)": cnt(lambda x: x.get("alt_methyl_bimodal") is True),
    "alt_methyl_untestable(<12 ALT reads)": cnt(lambda x: x.get("alt_methyl_bimodal") is None),
    "mut_methyl_coherent(ALT vs REF 甲基可分≥0.2)": cnt(lambda x: x.get("mut_methyl_coherent") is True),
    "mut_methyl_incoherent(<0.2)": cnt(lambda x: x.get("mut_methyl_coherent") is False),
    "note": "alt_single_hp=ALT 集中單一 HP→『誰帶突變』分群已由基因型定,甲基只需確認無次結構(負向篩選);"
            "alt_methyl_unimodal=ALT reads 甲基無雙峰=無偵測到次群(absence of evidence≠proof);"
            "mut_methyl_coherent=ALT 群甲基與 REF 明顯不同=mutation read 群聚 corroborate(但 cis-confounded,非獨立)。"}
json.dump({"summary": summary, "regions": out}, open(f"{DATA}/h3_unresolved_grouping.json", "w"), ensure_ascii=False, indent=1)
print("H3 UNRESOLVED GROUPING DONE", json.dumps(summary, ensure_ascii=False))
