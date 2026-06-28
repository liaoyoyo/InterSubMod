#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
甲基有界軟標記 — 整合兩用途(L3 auxiliary,非定論):
(b) 隱藏次結構旗標:對每個 genotype 群測「群內甲基是否雙峰(≥2 次群)」,扣 HP(germline)+CN(gain) 混淆。
    保守判準:GMM 2-comp BIC<1-comp + 兩峰均值差≥0.2 + 小組件≥4 reads。
    分類:unimodal / hp-explained(雙峰沿 HP=germline) / cn-flagged(CN-gain=multiplicity) / residual-candidate(扣完仍雙峰)。
(a) A_ambiguous 軟提示:對缺中間群(nested vs sibling 待定)的區,用甲基離 root 距離排序給「較可能」提示(標 L3)。
一次 BAM 掃描。compute batch(§13.0)。輸出 methyl_auxiliary_annotation.json。
"""
import json, os, sys
from collections import defaultdict, Counter
import numpy as np
import pysam
from sklearn.mixture import GaussianMixture
# 分塊平行: argv[1]=chunk_total(預設1=全跑), argv[2]=chunk_idx。strided 切分。
CHUNK_TOTAL = int(sys.argv[1]) if len(sys.argv) > 1 else 1
CHUNK_IDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0
OUTSUF = "" if CHUNK_TOTAL == 1 else f"_part{CHUNK_IDX}"

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TBAM = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
VCFD = "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
MAPQ = 20; MINREAD = 3; MINGRP = 12; SEP = 0.2; MINSUB = 4

det = {r["region"]: r for r in json.load(open(f"{DATA}/topology_per_region.json"))["detail"]}

_TBX = {}
def _tabix(p):
    if p not in _TBX:
        _TBX[p] = pysam.TabixFile(p) if os.path.exists(p) else None
    return _TBX[p]

def ref_alt_map(chrom, s, e):
    m = {}
    for src in ("tp", "fp"):
        p = f"{VCFD}/filtered_snv_{src}_{chrom}.vcf.gz"
        tbx = _tabix(p)
        if tbx is None: continue
        try:
            for ln in tbx.fetch(chrom, max(0, s - 1), e + 1):
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

def hp_norm(a):
    if not a.has_tag("HP"): return 0
    try: return int(str(a.get_tag("HP")).split("-")[0])
    except Exception: return 0

def bimodal_test(vals):
    """保守雙峰: GMM 2-comp BIC<1-comp + 峰均值差≥SEP + 小組件≥MINSUB。回 (is_bimodal, info)。"""
    x = np.array(vals).reshape(-1, 1)
    if len(x) < MINGRP: return False, None
    try:
        g1 = GaussianMixture(1, covariance_type="full", random_state=0).fit(x)
        g2 = GaussianMixture(2, covariance_type="full", random_state=0, n_init=1).fit(x)
    except Exception:
        return False, None
    if g2.bic(x) >= g1.bic(x): return False, None
    m = g2.means_.flatten()
    if abs(m[0] - m[1]) < SEP: return False, None
    lab = g2.predict(x)
    c0, c1 = int((lab == 0).sum()), int((lab == 1).sum())
    if min(c0, c1) < MINSUB: return False, None
    return True, {"means": [round(float(m[0]), 3), round(float(m[1]), 3)], "sizes": [c0, c1], "labels": lab.tolist()}

def cpg_meanbeta_per_read(meth_map, ids):
    """每 read 在區內 CpG 的平均 β(read 須覆蓋 ≥MINREAD CpG)。回 {read:meanβ}。"""
    out = {}
    for rn in ids:
        v = list(meth_map[rn].values())
        if len(v) >= MINREAD: out[rn] = float(np.mean(v))
    return out

tb = pysam.AlignmentFile(TBAM, "rb")
hidden = []      # (b) per-group bimodality flags
amb_hints = []   # (a) A_ambiguous soft ordering hints
flag_counter = Counter()
n_groups_tested = 0
items = list(det.items())[CHUNK_IDX::CHUNK_TOTAL]
READ_CAP = 500
for region, r in items:
    chrom, s = r["chrom"], r["start"]; e = int(region.split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som) < 2: continue
    is_amb = "ambiguous" in r.get("determinacy", "")
    meth = {}; geno = {}; hp = {}; n_read = 0
    for a in tb.fetch(chrom, s, e + 1):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < MAPQ: continue
        m = read_meth(a)
        if m is None: continue
        g = geno_str(a, som)
        if "-" in g: continue
        meth[a.query_name] = m; geno[a.query_name] = g; hp[a.query_name] = hp_norm(a)
        n_read += 1
        if n_read >= READ_CAP: break  # 高覆蓋區封頂(雙峰/排序用不到全部)
    grp = defaultdict(list)
    for rn, g in geno.items(): grp[g].append(rn)
    # ---- (b) 每個群測雙峰 ----
    for g, ids in grp.items():
        if len(ids) < MINGRP: continue
        mb = cpg_meanbeta_per_read(meth, ids)
        if len(mb) < MINGRP: continue
        rns = list(mb.keys()); vals = [mb[rn] for rn in rns]
        is_bi, info = bimodal_test(vals)
        n_groups_tested += 1
        if not is_bi:
            flag_counter["unimodal"] += 1; continue
        # 去混淆: 雙峰是否沿 HP?
        lab = info["labels"]
        hp0 = Counter(hp[rns[i]] for i in range(len(rns)) if lab[i] == 0)
        hp1 = Counter(hp[rns[i]] for i in range(len(rns)) if lab[i] == 1)
        def domhp(c):
            tot = sum(c.values()) or 1
            for k in (1, 2, 3):
                if c.get(k, 0) / tot >= 0.7: return k
            return 0
        d0, d1 = domhp(hp0), domhp(hp1)
        hp_explained = (d0 in (1, 2) and d1 in (1, 2) and d0 != d1)
        cn_flag = (r["cn"] == "gain")
        if hp_explained: flag = "hp-explained"
        elif cn_flag: flag = "cn-flagged"
        else: flag = "residual-candidate"
        flag_counter[flag] += 1
        hidden.append({"region": region, "genotype": g, "n_reads": len(ids), "cn": r["cn"],
                       "means": info["means"], "sizes": info["sizes"], "flag": flag,
                       "hp_modes": [d0, d1]})
    # ---- (a) A_ambiguous 軟提示 ----
    if "ambiguous" in r.get("determinacy", ""):
        rootg = "R" * len(som)
        hint = {"region": region, "n_sSNV": len(som), "cn": r["cn"], "ambig_nodes": r.get("ambig_nodes", 0),
                "status": "", "ordering": None}
        if rootg not in grp or len([1 for ids in grp.values() if len(ids) >= MINREAD]) < 2:
            hint["status"] = "無法產生提示(缺 RR root 或 <2 群覆蓋)"
        else:
            root_prof = defaultdict(list)
            for rn in grp[rootg]:
                for c, b in meth[rn].items(): root_prof[c].append(b)
            rootmean = {c: np.mean(v) for c, v in root_prof.items() if len(v) >= MINREAD}
            dists = []
            for g, ids in grp.items():
                if "A" not in g or len(ids) < MINREAD: continue
                gp = defaultdict(list)
                for rn in ids:
                    for c, b in meth[rn].items(): gp[c].append(b)
                gm = {c: np.mean(v) for c, v in gp.items() if len(v) >= MINREAD}
                shared = set(gm) & set(rootmean)
                if len(shared) < 3: continue
                d = float(np.mean([abs(gm[c] - rootmean[c]) for c in shared]))
                dists.append({"genotype": g, "n_mut": g.count("A"), "meth_dist_to_root": round(d, 3)})
            if len(dists) >= 2:
                dists.sort(key=lambda x: x["meth_dist_to_root"])
                hint["ordering"] = dists  # 距離小=較可能近祖先(L3 軟提示)
                hint["status"] = "已產生 L3 軟提示(甲基距離小→較可能近祖先;全域 ρ≈0.18 故低信心)"
            else:
                hint["status"] = "資料不足(<2 ALT 群有足夠共享 CpG)"
        amb_hints.append(hint)

# 統計
resid = [h for h in hidden if h["flag"] == "residual-candidate"]
resid_clean = [h for h in resid if h["cn"] in ("neutral", "loh")]
summary = {
    "b_hidden_substructure": {
        "n_groups_tested": n_groups_tested,
        "flag_dist": dict(flag_counter),
        "n_bimodal": len(hidden),
        "n_residual_candidate": len(resid),
        "n_residual_candidate_cnclean": len(resid_clean),
        "pct_residual_of_tested": round(100 * len(resid) / n_groups_tested, 2) if n_groups_tested else None,
    },
    "a_ambiguous_hint": {
        "n_ambiguous_regions": len(amb_hints),
        "n_hint_generated": sum(1 for h in amb_hints if h["ordering"]),
        "n_no_data": sum(1 for h in amb_hints if not h["ordering"]),
    },
    "params": {"MINGRP": MINGRP, "SEP": SEP, "MINSUB": MINSUB, "MAPQ": MAPQ},
    "interpretation": "residual-candidate = 扣 HP+CN 後群內甲基仍雙峰 = L3『可能隱藏次結構』候選(非確認,需 single-cell/更多 sSNV);A_ambiguous 提示為 L3 低信心(全域 ρ≈0.18)。",
}
json.dump({"summary": summary, "hidden_candidates": hidden, "ambiguous_hints": amb_hints},
          open(f"{DATA}/methyl_auxiliary_annotation{OUTSUF}.json", "w"), ensure_ascii=False, indent=1)
print(f"AUX ANNOTATION DONE part={CHUNK_IDX}/{CHUNK_TOTAL}")
print(json.dumps(summary, ensure_ascii=False, indent=1))
