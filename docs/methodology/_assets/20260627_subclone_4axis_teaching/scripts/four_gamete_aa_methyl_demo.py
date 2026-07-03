#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q1 實證 demo(2026-07-03):四配子(RR/RA/AR/AA 全有=成環)區,甲基能否判「AA 屬 AR 還是 RA 譜系」?
增強 four_gamete_aa_resolution.py:把 whole-read mean-β → 拆 near-i / near-j / distal 三帶 + 對稱統計 + 置換 null。

標籤/方法意義(先驗證再下結論):
  pair(i,j) 四配子:AA=(A,A) 與 AR=(A,R) 共享 locus-i ALT;AA 與 RA=(R,A) 共享 locus-j ALT。
  → cis-ASM 預測:near-i 帶 β(AA)≈β(AR)(共享 i 的 ALT cis 足跡)、near-j 帶 β(AA)≈β(RA)(共享 j)。
  → 「AA 離誰近」在 near 帶天生對稱(=重讀 genotype);唯一非循環 lineage 訊號在 distal 帶(離兩位點>1kb、無突變足跡)。
  同時閘:same-HP(AR/RA/AA 同 germline HP;跨 HP=兩突變不同染色體=AA doublet/CN artifact 非譜系)+ CN 分層。

輸出 four_gamete_aa_methyl_demo.json。§13.0 compute batch(先落檔→Read 驗→才寫報告)。seed 固定。
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
MAPQ = 20; MINREAD = 4; READ_CAP = 600; SOM_CAP = 30
NEAR_BP = 1000; MIN_BAND = 4; N_PERM = 2000
RNG = np.random.default_rng(20260703)

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

def read_betas(a):
    """回 {ref_pos(1-based): beta}"""
    try: mb = a.modified_bases
    except Exception: return {}
    if not mb: return {}
    qr = {q: rr for q, rr in a.get_aligned_pairs(matches_only=True)}
    out = {}
    for k, lst in mb.items():
        if k[2] != "m": continue
        for qp, ml in lst:
            rr = qr.get(qp)
            if rr is not None: out[rr + 1] = ml / 255.0
    return out

def hp_norm(a):
    if not a.has_tag("HP"): return 0
    try: return int(str(a.get_tag("HP")).split("-")[0])
    except Exception: return 0

def band_beta(betas, centers, near):
    """betas={pos:beta}; near=True 取 <=NEAR_BP any center; near=False 取 distal(>NEAR_BP all centers)。回 mean 或 None"""
    vals = []
    for p, b in betas.items():
        d = min(abs(p - c) for c in centers)
        if (near and d <= NEAR_BP) or ((not near) and d > NEAR_BP): vals.append(b)
    return float(np.mean(vals)) if vals else None

def perm_asym(betas_ar, betas_ra, betas_aa):
    """distal 帶:S=|Δ(AA,AR)|-|Δ(AA,RA)|;置換三群 read label 檢定 |S| 是否超越隨機。回 (S_obs, p, muAR,muRA,muAA,nAR,nRA,nAA) 或 None"""
    a = np.array([x for x in betas_ar if x is not None]); r = np.array([x for x in betas_ra if x is not None]); z = np.array([x for x in betas_aa if x is not None])
    if len(a) < MIN_BAND or len(r) < MIN_BAND or len(z) < MIN_BAND: return None
    muAR, muRA, muAA = a.mean(), r.mean(), z.mean()
    S_obs = abs(muAA - muAR) - abs(muAA - muRA)
    pool = np.concatenate([a, r, z]); n1, n2, n3 = len(a), len(r), len(z)
    idx = np.arange(len(pool)); cnt = 0
    for _ in range(N_PERM):
        RNG.shuffle(idx)
        pa = pool[idx[:n1]]; pr = pool[idx[n1:n1+n2]]; pz = pool[idx[n1+n2:]]
        S = abs(pz.mean() - pa.mean()) - abs(pz.mean() - pr.mean())
        if abs(S) >= abs(S_obs): cnt += 1
    return (round(float(S_obs), 4), (cnt + 1) / (N_PERM + 1),
            round(float(muAR), 3), round(float(muRA), 3), round(float(muAA), 3), n1, n2, n3)

def domhp(rds):
    c = Counter(rd["hp"] for rd in rds); tot = sum(c.values()) or 1
    for k in (1, 2, 3):
        if c.get(k, 0) / tot >= 0.6: return k
    return 0

tb = pysam.AlignmentFile(TBAM, "rb")
cases = []
for r in targets:
    chrom, s = r["chrom"], r["start"]; e = int(r["region"].split("-")[-1])
    ra = ref_alt_map(chrom, s, e)
    som_all = sorted([(p, ra[p][0], ra[p][1]) for p in ra if s <= p <= e])
    if len(som_all) < 2: continue
    reads = []; nr = 0; cov = Counter()
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
        reads.append({"hp": hp_norm(a), "betas": read_betas(a), "sb": sb})
        nr += 1
        if nr >= READ_CAP: break
    som = [t for t in som_all if t[0] in dict(cov.most_common(SOM_CAP))]
    for (pi, ri, ai), (pj, rj, aj) in combinations(som, 2):
        groups = defaultdict(list)
        for rd in reads:
            gi, gj = rd["sb"][pi], rd["sb"][pj]
            if gi in "RA" and gj in "RA": groups[gi + gj].append(rd)
        if not all(len(groups[g]) >= MINREAD for g in ("RR", "RA", "AR", "AA")): continue
        hp_ar, hp_ra, hp_aa = domhp(groups["AR"]), domhp(groups["RA"]), domhp(groups["AA"])
        same_hp = (hp_ar in (1, 2) and hp_ar == hp_ra == hp_aa)
        n_ar, n_ra, n_aa = len(groups["AR"]), len(groups["RA"]), len(groups["AA"])
        vaf_parent = "AR" if n_ar > n_ra else ("RA" if n_ra > n_ar else "tie")
        # near-i 帶(共享 locus i):cis 預測 AA≈AR;near-j 帶(共享 j):AA≈RA
        def gb(gname, centers, near): return [band_beta(rd["betas"], centers, near) for rd in groups[gname]]
        def mean_notnone(xs):
            v = [x for x in xs if x is not None]; return float(np.mean(v)) if len(v) >= MIN_BAND else None
        bi_ar, bi_ra, bi_aa = mean_notnone(gb("AR",[pi],True)), mean_notnone(gb("RA",[pi],True)), mean_notnone(gb("AA",[pi],True))
        bj_ar, bj_ra, bj_aa = mean_notnone(gb("AR",[pj],True)), mean_notnone(gb("RA",[pj],True)), mean_notnone(gb("AA",[pj],True))
        # distal 帶 lineage 測試 + 置換
        pa = perm_asym(gb("AR",[pi,pj],False), gb("RA",[pi,pj],False), gb("AA",[pi,pj],False))
        # near-band 綜合「AA 離誰近」(用 near-i,near-j 平均;=whole-read 近似)
        near_parent = None
        if None not in (bi_ar,bi_ra,bi_aa,bj_ar,bj_ra,bj_aa):
            dAR = (abs(bi_aa-bi_ar)+abs(bj_aa-bj_ar))/2; dRA = (abs(bi_aa-bi_ra)+abs(bj_aa-bj_ra))/2
            near_parent = "AR" if dAR < dRA else ("RA" if dRA < dAR else "tie")
        distal_parent = None; distal_p = None; S_obs = None
        if pa is not None:
            S_obs, distal_p = pa[0], pa[1]
            distal_parent = "AR" if S_obs < 0 else ("RA" if S_obs > 0 else "tie")
        cases.append({"region": r["region"], "cn": r["cn"], "pair": [pi, pj],
                      "n_RR": len(groups["RR"]), "n_AR": n_ar, "n_RA": n_ra, "n_AA": n_aa,
                      "same_hp": same_hp, "hp": [hp_ar, hp_ra, hp_aa],
                      "vaf_parent": vaf_parent,
                      # cis 驗證:near-i 帶 AA vs AR(共享 i,應近) vs RA(差 i,應遠);near-j 反之
                      "near_i_AR_RA_AA": [round(bi_ar,3) if bi_ar else None, round(bi_ra,3) if bi_ra else None, round(bi_aa,3) if bi_aa else None],
                      "near_j_AR_RA_AA": [round(bj_ar,3) if bj_ar else None, round(bj_ra,3) if bj_ra else None, round(bj_aa,3) if bj_aa else None],
                      "near_parent": near_parent,
                      "distal_S(neg=AR,pos=RA)": S_obs, "distal_perm_p": distal_p, "distal_parent": distal_parent,
                      "distal_measurable": pa is not None})

# ── 聚合 ──
same = [c for c in cases if c["same_hp"]]
clean = [c for c in same if c["cn"] in ("neutral", "loh")]
distal_ok = [c for c in cases if c["distal_measurable"]]
# cis 驗證:near-i 帶 AA 應更貼 AR(共享 i ALT)而非 RA
cis_i = [(c["near_i_AR_RA_AA"]) for c in cases if None not in c["near_i_AR_RA_AA"]]
cis_j = [(c["near_j_AR_RA_AA"]) for c in cases if None not in c["near_j_AR_RA_AA"]]
def cis_track(rows, shared_idx):  # shared: AA vs (AR if i else RA);  非共享對照
    # rows=[bAR,bRA,bAA]; near-i: AA 共享 AR(idx0);near-j: AA 共享 RA(idx1)
    d_share = np.mean([abs(x[2]-x[shared_idx]) for x in rows]) if rows else None
    d_nonsh = np.mean([abs(x[2]-x[1-shared_idx]) for x in rows]) if rows else None
    return (round(float(d_share),4) if d_share is not None else None,
            round(float(d_nonsh),4) if d_nonsh is not None else None, len(rows))
cis_i_share, cis_i_nonsh, n_ci = cis_track(cis_i, 0)   # near-i 共享=AR
cis_j_share, cis_j_nonsh, n_cj = cis_track(cis_j, 1)   # near-j 共享=RA
# distal 方向可重現性:置換顯著(p<0.05)比例 = distal 真能判別的比例
distal_sig = [c for c in distal_ok if c["distal_perm_p"] is not None and c["distal_perm_p"] < 0.05]
# near-based parent vs distal-based parent 一致性(兩者皆有)
both_pp = [c for c in cases if c["near_parent"] in ("AR","RA") and c["distal_parent"] in ("AR","RA")]
np_dp_agree = [c for c in both_pp if c["near_parent"] == c["distal_parent"]]

summary = {
    "sample": "HCC1395", "partial_flag": "single-sample HCC1395 / 四配子(incompatible)子集",
    "n_targets(incompatible+dropnoise)": len(targets),
    "n_4gamete_pairs": len(cases),
    "funnel": {
        "same_hp(前提成立)": len(same),
        "cross_or_mixed_hp(前提不成立→AA=doublet/CN artifact 非譜系)": len(cases) - len(same),
        "same_hp_AND_cn_clean(neutral/loh)": len(clean),
        "resolvable_clean(same-HP+CN-clean+distal 可測)": len([c for c in clean if c["distal_measurable"]]),
    },
    "cis_label_check(方法/標籤意義驗證)": {
        "near_i_dist_AA_to_shared(AR)": cis_i_share, "near_i_dist_AA_to_nonshared(RA)": cis_i_nonsh, "n": n_ci,
        "near_j_dist_AA_to_shared(RA)": cis_j_share, "near_j_dist_AA_to_nonshared(AR)": cis_j_nonsh, "n": n_cj,
        "解讀": "若 dist_AA_to_shared << dist_AA_to_nonshared → cis 結構成立(AA 在共享位點 track 對應 parent = 重讀 genotype),near 帶對稱不可判譜系",
    },
    "distal_lineage_test": {
        "n_distal_measurable": len(distal_ok),
        "n_distal_perm_sig(p<0.05=distal 真能判別)": len(distal_sig),
        "frac_distal_sig": round(len(distal_sig)/len(distal_ok), 3) if distal_ok else None,
        "near_vs_distal_parent_agree": f"{len(np_dp_agree)}/{len(both_pp)}",
    },
    "verdict": "見報告;預測:cis 結構成立(near 對稱)+ distal 顯著比例≈隨機水準 → 甲基無法判四配子拓撲方向。",
}
json.dump({"summary": summary, "cases": cases}, open(f"{DATA}/four_gamete_aa_methyl_demo.json", "w"), ensure_ascii=False, indent=1)
print("=== FOUR-GAMETE AA METHYL DEMO DONE ===")
print(json.dumps(summary, ensure_ascii=False, indent=1))
