#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
甲基排序 — ISM read-level BERNOULLI 距離 vs mean-Δβ 對照(2026-07-04)。
完全復用 methyl_ordering_pilot.py 的選區/讀取/near-distal/permutation,唯一差別:
group-to-root 距離同時算兩種 —
  (A) mean-Δβ(pilot 原法):|group per-CpG 平均β − root 平均β| 的平均(group 層 L1)。
  (B) read-level Bernoulli(ISM C++ calculate_bernoulli 同公式,between-group 平均):
      每對(group read i, root read j)在共同 CpG 上 δ=p(1-q)+(1-p)q,權重 w=2|p-.5|·2|q-.5|,
      dist_ij=Σwδ/Σw;group-to-root = 所有 read-pair dist_ij 的平均(= compute_group_distances between)。
測「距離是否單調隨基因型譜系深度(A_determined ground truth)」的 pooled Spearman + permutation p。
若 (B) 的 ρ/p 顯著優於 (A) → Bernoulli 值得取代;若相近/更差 → 確認「Bernoulli 救不了、天花板在訊號本身」。
§13.0 compute batch:跑完落 JSON,另批讀+寫報告。用法:python3 methyl_ordering_bernoulli_compare.py [N_CAP] [READCAP]
"""
import json, os, sys
from collections import defaultdict
import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3; NEAR_BP = 1000; MIN_COMMON = 3  # C_min(Bernoulli 最少共同 CpG)
N_CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 0        # 0=全跑
READCAP = int(sys.argv[2]) if len(sys.argv) > 2 else 50     # 每群取樣上限(bound 成對成本)

det = {r["region"]: r for r in json.load(open(f"{DATA}/topology_per_region.json"))["detail"]}
hpa = {r["region"]: r["hp_alignment"] for r in json.load(open(f"{DATA}/hp_alignment_fullscan.json"))["regions"]}

def node_depth(edges):
    children = defaultdict(list)
    for p, c in edges: children[p].append(c)
    depth = {}; stack = [("ROOT", -1)]
    while stack:
        n, d = stack.pop()
        if n != "ROOT": depth[n] = d
        for c in children.get(n, []): stack.append((c, d + 1))
    return depth

cand = []
for region, r in det.items():
    if "A_determined" not in r.get("determinacy", ""): continue
    if hpa.get(region) != "SAME-HP": continue
    if r["cn"] not in ("neutral", "loh"): continue
    dep = node_depth(r.get("edges", []))
    if len(set(dep.values())) >= 2: cand.append((region, r, dep))
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
pool_beta = {"all": [], "near": [], "distal": []}
pool_bern = {"all": [], "near": [], "distal": []}
n_skip = 0; n_used = 0

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
    if rootg not in grp: n_skip += 1; continue

    # 每 read 預算 near/distal CpG 快取:readcache[rn][which] = {c: beta}
    def is_near(c): return any(abs(c - p) <= NEAR_BP for p in sompos)
    readcache = {}
    for rn in meth:
        allm = meth[rn]; nearm = {}; distm = {}
        for c, b in allm.items():
            (nearm if is_near(c) else distm)[c] = b
        readcache[rn] = {"all": allm, "near": nearm, "distal": distm}

    def cpg_profile(ids, which):
        cpg = defaultdict(list)
        for rn in ids:
            for c, b in readcache[rn][which].items(): cpg[c].append(b)
        return {c: float(np.mean(v)) for c, v in cpg.items() if len(v) >= MINREAD}

    def bern_between(ids_g, ids_root, which):
        gi = ids_g[:READCAP]; gr = ids_root[:READCAP]; dists = []
        for ri in gi:
            vi = readcache[ri][which]
            if not vi: continue
            for rj in gr:
                vj = readcache[rj][which]
                if not vj: continue
                common = vi.keys() & vj.keys()
                if len(common) < MIN_COMMON: continue
                w = 0.0; wd = 0.0
                for c in common:
                    p = vi[c]; q = vj[c]
                    wk = (2.0 * abs(p - 0.5)) * (2.0 * abs(q - 0.5))
                    w += wk; wd += wk * (p * (1.0 - q) + (1.0 - p) * q)
                if w > 1e-9: dists.append(wd / w)
        return float(np.mean(dists)) if dists else None

    for which in ("all", "near", "distal"):
        root_prof = cpg_profile(grp[rootg], which)
        for g, ids in grp.items():
            if "A" not in g or g not in dep: continue
            # (A) mean-Δβ
            gp = cpg_profile(ids, which); shared = set(gp) & set(root_prof)
            db = float(np.mean([abs(gp[c] - root_prof[c]) for c in shared])) if len(shared) >= 3 else None
            # (B) read-level Bernoulli
            bd = bern_between(ids, grp[rootg], which)
            if db is not None: pool_beta[which].append((dep[g], db))
            if bd is not None: pool_bern[which].append((dep[g], bd))
    n_used += 1

def pooled(pairs):
    if len(pairs) < 5: return {"n": len(pairs), "rho": None}
    x = np.array([p[0] for p in pairs], float); y = np.array([p[1] for p in pairs], float)
    rho = spearman(x, y)
    perm = []
    for k in range(200):  # 確定性打亂(無 random,對齊 pilot)
        idx = [(k * 7919 + i * 104729) % len(y) for i in range(len(y))]
        seen = {}; order = []
        for v in idx:
            while v in seen: v = (v + 1) % len(y)
            seen[v] = 1; order.append(v)
        pr = spearman(x, y[np.array(order)])
        if pr is not None: perm.append(pr)
    p = (sum(1 for pr in perm if abs(pr) >= abs(rho)) + 1) / (len(perm) + 1) if rho is not None and perm else None
    return {"n": len(pairs), "rho": rho, "perm_p": p}

out = {
    "n_candidate_regions": len(cand), "n_used": n_used, "n_skipped": n_skip,
    "read_cap_per_group": READCAP, "min_common_cpg": MIN_COMMON, "near_bp": NEAR_BP,
    "mean_delta_beta": {w: pooled(pool_beta[w]) for w in ("all", "near", "distal")},
    "read_level_bernoulli": {w: pooled(pool_bern[w]) for w in ("all", "near", "distal")},
    "note": "(A) mean-Δβ 復現 pilot(sanity: 應≈ρ0.18/p0.06-0.08);(B) read-level Bernoulli(ISM C++ 公式,between-group 平均)。比 ρ/perm_p 判 Bernoulli 是否救得了 ordering。",
}
outp = os.path.join(DATA, "methyl_ordering_bernoulli_compare.json")
with open(outp, "w", encoding="utf-8") as f: json.dump(out, f, ensure_ascii=False, indent=1)
print("OK wrote", outp)
print(json.dumps(out, ensure_ascii=False, indent=1))
