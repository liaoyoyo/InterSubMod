#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
甲基排序啟發式驗證:在「基因型已能定先後(A_determined)」的 nested 區當 GROUND TRUTH,
測「甲基離 root(RR 群)距離」是否單調隨 基因型譜系深度 增加(Spearman)。
🔑 關鍵:把 CpG 分 near(<1kb 任一 somatic)/distal(>1kb),分開測 —
   若只有 near 隨深度漲=somatic-cis 痕跡累積(假時鐘);distal 也漲=真譜系時鐘。
限 SAME-HP(germline-free) + CN-clean(neutral/loh,排 CN-gain multiplicity)。
compute batch(§13.0):跑完落 JSON,另批讀+寫報告。用法:python3 methyl_ordering_pilot.py [N_CAP]
"""
import json, os, sys
from collections import defaultdict, Counter
import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3; NEAR_BP = 1000
N_CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 0  # 0=全跑

det = {r["region"]: r for r in json.load(open(f"{DATA}/topology_per_region.json"))["detail"]}
hpa = {r["region"]: r["hp_alignment"] for r in json.load(open(f"{DATA}/hp_alignment_fullscan.json"))["regions"]}

def node_depth(edges):
    """從 edges [(parent,child)] BFS 算每個 genotype 節點離 ROOT 的深度。"""
    children = defaultdict(list)
    for p, c in edges: children[p].append(c)
    depth = {}
    stack = [("ROOT", -1)]
    while stack:
        n, d = stack.pop()
        if n != "ROOT": depth[n] = d
        for c in children.get(n, []): stack.append((c, d + 1))
    return depth  # genotype-string -> depth (ROOT 的直接子=0)

# 選區:A_determined + SAME-HP + CN-clean + 有 ≥2 不同深度的 ALT 節點
cand = []
for region, r in det.items():
    if "A_determined" not in r.get("determinacy", ""): continue
    if hpa.get(region) != "SAME-HP": continue
    if r["cn"] not in ("neutral", "loh"): continue
    dep = node_depth(r.get("edges", []))
    if len(set(dep.values())) >= 2:  # 至少兩種深度才能測單調
        cand.append((region, r, dep))
if N_CAP: cand = cand[:N_CAP]

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

def read_meth(a):
    try: mb = a.modified_bases
    except Exception: return None
    qr = {q: rr for q, rr in a.get_aligned_pairs(matches_only=True)}
    meth = {}
    if mb:
        for k, lst in mb.items():
            if k[2] != "m": continue
            for qpos, mlq in lst:
                rr = qr.get(qpos)
                if rr is not None: meth[rr] = mlq / 255.0
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

def spearman(x, y):
    if len(x) < 3: return None
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    if np.std(rx) == 0 or np.std(ry) == 0: return None
    return float(np.corrcoef(rx, ry)[0, 1])

tb = pysam.AlignmentFile(TBAM, "rb")
per_region = []
pool = {"all": [], "near": [], "distal": []}  # (depth, methdist)
nested_pairs = {"all": [0, 0], "near": [0, 0], "distal": [0, 0]}  # [deeper>shallower, total]
n_skip = 0
for region, r, dep in cand:
    chrom, s = r["chrom"], r["start"]; e = int(region.split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2: n_skip += 1; continue
    sompos = [p for p, _, _ in som]
    meth = {}; geno = {}
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        m = read_meth(a)
        if m is None: continue
        g = geno_str(a, som)
        if "-" in g: continue
        meth[a.query_name] = m; geno[a.query_name] = g
    grp = defaultdict(list)
    for rn, g in geno.items(): grp[g].append(rn)
    grp = {g: ids for g, ids in grp.items() if len(ids) >= MINREAD}
    rootg = "R" * len(som)
    if rootg not in grp: n_skip += 1; continue  # 需 RR 祖先群當參考
    # per-group per-CpG β,分 near/distal
    def cpg_profile(ids, which):
        cpg = defaultdict(list)
        for rn in ids:
            for c, b in meth[rn].items():
                near = any(abs(c - p) <= NEAR_BP for p in sompos)
                if which == "all" or (which == "near" and near) or (which == "distal" and not near):
                    cpg[c].append(b)
        return {c: np.mean(v) for c, v in cpg.items() if len(v) >= MINREAD}
    rec = {"region": region, "cn": r["cn"], "n_som": len(som), "n_groups": len(grp), "rho": {}}
    for which in ("all", "near", "distal"):
        root_prof = cpg_profile(grp[rootg], which)
        depths = []; dists = []
        for g, ids in grp.items():
            if "A" not in g: continue  # 只測 ALT 群(subclone)
            if g not in dep: continue
            gp = cpg_profile(ids, which)
            shared = set(gp) & set(root_prof)
            if len(shared) < 3: continue
            d = np.mean([abs(gp[c] - root_prof[c]) for c in shared])
            depths.append(dep[g]); dists.append(d)
            pool[which].append((dep[g], d))
        rho = spearman(depths, dists)
        rec["rho"][which] = rho
        # nested pair: 任兩 ALT 群深度不同→深的距離應較大
        for i in range(len(depths)):
            for j in range(len(depths)):
                if depths[i] > depths[j]:
                    nested_pairs[which][1] += 1
                    if dists[i] > dists[j]: nested_pairs[which][0] += 1
    per_region.append(rec)

# pooled Spearman + permutation null
def pooled_stats(which):
    pairs = pool[which]
    if len(pairs) < 5: return {"n": len(pairs), "rho": None}
    x = np.array([p[0] for p in pairs]); y = np.array([p[1] for p in pairs])
    rho = spearman(x, y)
    # permutation null(打亂 depth)
    rng = list(range(len(y)))
    perm = []
    for k in range(200):
        idx = [(k * 7919 + i * 104729) % len(y) for i in range(len(y))]  # 確定性打亂(無 random)
        seen = {}; order = []
        for v in idx:
            while v in seen: v = (v + 1) % len(y)
            seen[v] = 1; order.append(v)
        pr = spearman(x, y[order])
        if pr is not None: perm.append(pr)
    p = (sum(1 for pr in perm if abs(pr) >= abs(rho)) + 1) / (len(perm) + 1) if rho is not None and perm else None
    return {"n": len(pairs), "rho": rho, "perm_p": p, "perm_mean": float(np.mean(perm)) if perm else None}

out = {"n_candidate_regions": len(cand), "n_used": len(per_region), "n_skipped": n_skip,
       "near_bp": NEAR_BP,
       "pooled": {w: pooled_stats(w) for w in ("all", "near", "distal")},
       "nested_pair_frac": {w: (round(nested_pairs[w][0] / nested_pairs[w][1], 3) if nested_pairs[w][1] else None,
                                f"{nested_pairs[w][0]}/{nested_pairs[w][1]}") for w in ("all", "near", "distal")},
       "per_region_rho_summary": {w: {"median": round(float(np.median([r["rho"][w] for r in per_region if r["rho"][w] is not None])), 3) if [r for r in per_region if r["rho"][w] is not None] else None,
                                      "n": sum(1 for r in per_region if r["rho"][w] is not None)} for w in ("all", "near", "distal")},
       "interpretation": "distal rho>>0 且 perm_p<0.05 → 真譜系時鐘(遠離突變仍隨深度漂移);只有 near rho>0 → somatic-cis 痕跡累積(假時鐘);全 rho≈0 → 甲基不能排序。"}
json.dump({"summary": out, "per_region": per_region}, open(f"{DATA}/methyl_ordering_pilot.json", "w"), ensure_ascii=False, indent=1)
print("ORDERING PILOT DONE", json.dumps(out, ensure_ascii=False, indent=1))
