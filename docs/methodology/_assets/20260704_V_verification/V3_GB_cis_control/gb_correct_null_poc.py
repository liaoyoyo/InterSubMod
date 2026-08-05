#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V3 G-B correct-null PoC — HP-ORIENTED matched-normal cis-control on the 43 CLEAN
CROSS-HP (CN-neutral/loss) regions that the scope pilot only *classified*, never ran.

For a CROSS-HP region the two ALT tumor populations sit on DIFFERENT germline HPs,
so tumor_dbeta = tA - tB genuinely contains a germline-ASM(Ha vs Hb) component, and
the HP-ORIENTED normal baseline nHa - nHb is the CORRECT null:
    residual = tumor_dbeta - normal_dbeta   (germline-ASM removed -> somatic-specific)

Contrast with the SAME-HP pilot: there tA,tB share one HP so nH1-nH2 is ORTHOGONAL
(corr~=0 tautology). Here we EXPECT corr(tumor,normal)>0 if germline-ASM is a real
shared component, and the subtraction should actually change |dbeta|.

Reuses read-parsing verbatim from pilot_cis_control.py + classify_hp_alignment_all.py.
Output: gb_correct_null_poc.json (every number grep-able).
"""
import json, os
from collections import defaultdict, Counter
import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
# data dir holds the frozen inputs; region list written by the repl filter step
DATA = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data"
OUTDIR = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260704_V_verification/V3_GB_cis_control"
os.makedirs(OUTDIR, exist_ok=True)
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3
SIG = 0.2      # tumor "signal" threshold (existing ASM caliber)
CLEAN = 0.1    # residual "collapsed" threshold

regions = json.load(open(f"{DATA}/_v3_clean_cross_hp_regions.json"))

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

def read_meth(a):
    try: mb = a.modified_bases
    except Exception: return None
    qr = {q: r for q, r in a.get_aligned_pairs(matches_only=True)}
    meth = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != "m": continue      # 5mC only
            for qpos, mlq in lst:
                r = qr.get(qpos)
                if r is not None: meth[r] = mlq / 255.0
    return meth

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

def per_cpg(read_meth_list, idset):
    cpg = defaultdict(list)
    for rn in idset:
        for cc, b in read_meth_list[rn].items(): cpg[cc].append(b)
    return {cc: float(np.mean(v)) for cc, v in cpg.items() if len(v) >= MINREAD}

def dom_hp(cnt):
    tot = sum(cnt.values()) or 1
    for h in (1, 2, 3):
        if cnt.get(h, 0) / tot >= 0.6: return h
    return 0  # mixed

tb = pysam.AlignmentFile(TBAM, "rb"); nb = pysam.AlignmentFile(NBAM, "rb")
out_regions = []; all_pairs = []
n_no_som = n_lt2 = n_notcross = n_nonormal = n_ok = 0

for reg in regions:
    chrom = reg.split(":")[0]
    s = int(reg.split(":")[1].split("-")[0]); e = int(reg.split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2:
        n_no_som += 1; out_regions.append({"region": reg, "note": "som<2"}); continue
    # tumor reads: meth + genotype + HP
    tmeth = {}; tgeno = {}; thp = {}
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        m = read_meth(a)
        if m is None: continue
        g = geno_str(a, som)
        if "-" in g: continue
        tmeth[a.query_name] = m; tgeno[a.query_name] = g; thp[a.query_name] = hp_norm(a)
    pop = defaultdict(list)
    for rn, g in tgeno.items(): pop[g].append(rn)
    # ALT pops (contain 'A'), >=MINREAD, largest first
    altpops = sorted([(g, ids) for g, ids in pop.items() if len(ids) >= MINREAD and "A" in g],
                     key=lambda x: -len(x[1]))
    if len(altpops) < 2:
        n_lt2 += 1; out_regions.append({"region": reg, "note": "tumor ALT pops<2"}); continue
    gA, idsA = altpops[0]; gB, idsB = altpops[1]
    hpA = dom_hp(Counter(thp[rn] for rn in idsA))
    hpB = dom_hp(Counter(thp[rn] for rn in idsB))
    if hpA == 0 or hpB == 0 or hpA == hpB:
        n_notcross += 1
        out_regions.append({"region": reg, "note": f"not read-level CROSS (hpA={hpA},hpB={hpB})",
                            "popA": gA, "popB": gB, "hpA": hpA, "hpB": hpB}); continue
    # normal reads on the SAME two HPs (oriented: A->hpA, B->hpB)
    nmeth = {}; nhp = {}
    for a in nb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        if not a.has_tag("HP"): continue
        m = read_meth(a)
        if m is None: continue
        nmeth[a.query_name] = m; nhp[a.query_name] = hp_norm(a)
    nA = [rn for rn, h in nhp.items() if h == hpA]
    nB = [rn for rn, h in nhp.items() if h == hpB]
    if len(nA) < MINREAD or len(nB) < MINREAD:
        n_nonormal += 1
        out_regions.append({"region": reg, "note": f"normal HP{hpA}/HP{hpB} reads<{MINREAD} ({len(nA)}/{len(nB)})",
                            "popA": gA, "popB": gB, "hpA": hpA, "hpB": hpB}); continue
    tA = per_cpg(tmeth, set(idsA)); tB = per_cpg(tmeth, set(idsB))
    nHa = per_cpg(nmeth, set(nA)); nHb = per_cpg(nmeth, set(nB))
    shared = set(tA) & set(tB) & set(nHa) & set(nHb)
    rec = {"region": reg, "cn": None, "popA": gA, "popB": gB, "hpA": hpA, "hpB": hpB,
           "tumor_nA": len(idsA), "tumor_nB": len(idsB), "normal_nA": len(nA), "normal_nB": len(nB),
           "n_cpg_shared": len(shared), "cpg": []}
    for cc in shared:
        td = round(tA[cc] - tB[cc], 3)          # oriented A(hpA) - B(hpB)
        nd = round(nHa[cc] - nHb[cc], 3)         # germline-ASM(hpA - hpB) = CORRECT null
        res = round(td - nd, 3)                   # residual = somatic-specific
        rec["cpg"].append({"cpg": cc, "tumor_dbeta": td, "normal_dbeta": nd, "residual": res})
        all_pairs.append({"region": reg, "cpg": cc, "tumor_dbeta": td, "normal_dbeta": nd, "residual": res})
    if rec["cpg"]:
        n_ok += 1
    out_regions.append(rec)

# ---- distribution / metrics ----
def absmed(xs): return round(float(np.median([abs(x) for x in xs])), 3) if xs else None
def absp90(xs): return round(float(np.percentile([abs(x) for x in xs], 90)), 3) if xs else None
def absmax(xs): return round(float(max(abs(x) for x in xs)), 3) if xs else None

td_all = [p["tumor_dbeta"] for p in all_pairs]
nd_all = [p["normal_dbeta"] for p in all_pairs]
res_all = [p["residual"] for p in all_pairs]

# tumor-significant CpGs and what happens to their residual
sig = [p for p in all_pairs if abs(p["tumor_dbeta"]) >= SIG]
sig_res_survive = [p for p in sig if abs(p["residual"]) >= SIG]     # germline-clean somatic candidate
sig_res_collapse = [p for p in sig if abs(p["residual"]) < CLEAN]   # explained by germline-ASM

corr = None
if len(td_all) >= 3 and np.std(td_all) > 0 and np.std(nd_all) > 0:
    corr = round(float(np.corrcoef(td_all, nd_all)[0, 1]), 3)

summary = {
    "design": "HP-oriented matched-normal cis-control on clean CROSS-HP regions (CORRECT null)",
    "n_regions_input": len(regions),
    "funnel": {"som<2": n_no_som, "tumor_ALTpops<2": n_lt2, "not_read_level_CROSS": n_notcross,
               "normal_HP_reads<3": n_nonormal, "regions_with_cpg": n_ok},
    "n_cpg_pairs_total": len(all_pairs),
    "tumor_dbeta_abs": {"median": absmed(td_all), "p90": absp90(td_all), "max": absmax(td_all)},
    "normal_dbeta_abs": {"median": absmed(nd_all), "p90": absp90(nd_all), "max": absmax(nd_all)},
    "residual_abs": {"median": absmed(res_all), "p90": absp90(res_all), "max": absmax(res_all)},
    "corr_tumor_vs_normal_dbeta": corr,
    "corr_interpretation": "CROSS-HP expects POSITIVE (shared germline-ASM); ~0 would mean subtraction is inert",
    "tumor_sig_cpg": {"n": len(sig), "pct_of_total": round(100*len(sig)/max(1,len(all_pairs)),1)},
    "of_sig_residual_survives_ge0.2": {"n": len(sig_res_survive),
        "pct": round(100*len(sig_res_survive)/max(1,len(sig)),1)},
    "of_sig_residual_collapses_lt0.1": {"n": len(sig_res_collapse),
        "pct": round(100*len(sig_res_collapse)/max(1,len(sig)),1)},
    "thresholds": {"tumor_sig": SIG, "residual_clean": CLEAN, "minread": MINREAD, "mapq": MAPQ},
}
json.dump({"summary": summary, "regions": out_regions, "all_cpg_pairs": all_pairs},
          open(f"{OUTDIR}/gb_correct_null_poc.json", "w"), ensure_ascii=False, indent=1)
print("PoC DONE", json.dumps(summary, ensure_ascii=False, indent=1))
