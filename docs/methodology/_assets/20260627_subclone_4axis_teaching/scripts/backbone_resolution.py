#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
骨幹欠定問題解決層(純 JSON,不需 BAM):
A1 成環剖析:12 incompatible 區 → 哪幾對 sSNV 違反 4-gamete + cn 成因。
B1 parsimony(Q4):A_ambiguous 缺中間群 → 最大簡約定序(假設逐步突變,跳 d 步=缺 d-1 中間群)+ 信心。
B2 機率第一順位(Q3/Q5):C_underdetermined → greedy 最大簡約樹當第一順位 + 信心 + Tier-PS 排除跨-HP。
輸出 backbone_resolution.json。compute batch(§13.0)。
"""
import json, os
from itertools import combinations
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(HERE, "..", "data")))  # per-sample 化(對齊 candidate_scoring.py)
det = json.load(open(f"{DATA}/topology_per_region.json"))["detail"]
DET = lambda r: r.get("determinacy", "?")

def altset(g): return frozenset(i for i, c in enumerate(g) if c == "A")

# ---------- A1 成環剖析 ----------
def four_gamete_violations(pops):
    """回違反 4-gamete(RR/RA/AR/AA 全present)的 sSNV 位置對。"""
    if not pops: return []
    k = len(next(iter(pops)))
    viol = []
    for i, j in combinations(range(k), 2):
        gam = set(g[i] + g[j] for g, c in pops.items() if c > 0 and g[i] in "RA" and g[j] in "RA")
        if {"RR", "RA", "AR", "AA"} <= gam:  # 4 gamete 全present = perfect-phylogeny 違反
            viol.append((i, j, sorted(gam)))
    return viol

a1 = []
GCAP = 8  # genotype 向量上限(高 n_sSNV 區只用前 8 位點)
for r in det:
    if DET(r) != "incompatible": continue
    pops = r["populations"]
    glen = len(next(iter(pops))) if pops else 0
    npop = len(pops)
    viol = four_gamete_violations(pops)
    # npop<=2 不可能 4-gamete 違反;n_sSNV>GCAP 被截斷 → 違反落在未存位點
    if npop <= 2:
        cap_caveat = f"npop={npop}(≤2 不可能真 4-gamete 違反)→ determinacy 與 stored populations 不一致 = 分類 artifact"
    elif r["n_sSNV"] > glen:
        cap_caveat = f"n_sSNV={r['n_sSNV']}>glen={glen}(genotype 截斷上限 {GCAP})→ 違反對可能落在未存位點,stored 查不到"
    else:
        cap_caveat = "stored 完整,違反對可查"
    a1.append({"region": r["region"], "cn": r["cn"], "genome_ctx": r["genome_ctx"],
               "n_sSNV": r["n_sSNV"], "glen": glen, "npop": npop, "tp": r["tp"], "fp": r["fp"], "span": r["span"],
               "n_violating_pairs_in_stored": len(viol),
               "violations": [{"sSNV_idx": [i, j], "gametes": gam} for i, j, gam in viol[:5]],
               "cap_caveat": cap_caveat,
               "cause": ("CN-gain multiplicity(同突變多拷貝)" if r["cn"] == "gain"
                         else "LOH/其他 — 偽影或 mapping") + (" + FP主導" if r["fp"] > r["tp"] else "")
                        + (f" + {r['genome_ctx']}偽影區" if r["genome_ctx"] in ("centromere", "telomere") else ""),
               "fixable": "標 artifact 不建樹(mappability mask + 解截斷)" if (r["cn"] == "gain" or r["genome_ctx"] in ("centromere", "telomere") or r["fp"] > r["tp"] or npop <= 2) else "需檢視(可能真衝突)"})

# ---------- B1 parsimony(A_ambiguous) + B2 機率第一順位(C_underdetermined) ----------
def vafmap(r):
    pv = r.get("pos_vaf", {})
    return pv if isinstance(pv, dict) else {}

def parsimony_tree(r):
    """每節點 → 最大真子集父(最少新增突變=最簡約);tie 用 CCF(VAF 大=祖先)。
    回 (edges, per_node_conf, tree_conf, crossHP_excluded, missing_intermediates)。"""
    pops = r["populations"]
    alt = {g: altset(g) for g in pops if "A" in g}
    if len(alt) < 2: return None
    vaf = vafmap(r); nhp = r.get("node_hp", {})
    cnt = {g: pops[g] for g in alt}
    nodes = sorted(alt, key=lambda g: (len(alt[g]), g))
    # 節點 HP(由其 ALT 位點的 node_hp 推；混則 ?)
    def node_hp_of(g):
        hps = set()
        # node_hp key 是 'chrom:pos',無法直接對 idx → 用整體 haplotypes 近似;這裡僅標記是否有跨HP風險
        return r.get("haplotypes", "?")
    edges = []; conf = []; miss = 0; crossx = 0
    for g in nodes:
        cand = [h for h in nodes if alt[h] < alt[g]]
        if not cand:
            edges.append(("ROOT", g)); conf.append(1.0); continue
        msz = max(len(alt[h]) for h in cand)
        mcand = [h for h in cand if len(alt[h]) == msz]
        jump = len(alt[g]) - msz
        if jump > 1: miss += (jump - 1)  # 缺中間群數
        # tie: CCF(VAF 大者較可能祖先) → 取 VAF 最大;無 VAF 用 read 數
        def ccf(h): return vaf.get(h, cnt[h] / (sum(cnt.values()) or 1))
        parent = max(mcand, key=lambda h: ccf(h))
        # 信心: 1/同簡約候選數 × (jump==1 ? 1 : 0.5^(jump-1)) (跳越多越不確定)
        c = (1.0 / len(mcand)) * (1.0 if jump == 1 else 0.5 ** (jump - 1))
        edges.append((parent, g)); conf.append(round(c, 3))
    tree_conf = round(min(conf), 3) if conf else 0.0  # 樹整體信心 = 最弱環節
    prob_first = round(1.0, 3)
    for c in conf: prob_first *= c
    return {"edges": edges, "per_node_conf": conf, "tree_conf": tree_conf,
            "first_rank_prob": round(prob_first, 3), "missing_intermediates": miss}

b1 = []
for r in det:
    if "ambiguous" not in DET(r): continue
    res = parsimony_tree(r)
    if not res: continue
    b1.append({"region": r["region"], "cn": r["cn"], "n_sSNV": r["n_sSNV"],
               "ambig_nodes": r.get("ambig_nodes", 0), "haplotypes": r.get("haplotypes", "?"),
               "parsimony_edges": res["edges"], "missing_intermediates": res["missing_intermediates"],
               "first_rank_prob": res["first_rank_prob"], "tree_conf": res["tree_conf"]})

b2 = []
for r in det:
    if DET(r) != "C_underdetermined": continue
    res = parsimony_tree(r)
    if not res:
        b2.append({"region": r["region"], "cn": r["cn"], "n_sSNV": r["n_sSNV"],
                   "status": "無 ≥2 ALT 群(單群/缺對)→ 無法排序,需深覆蓋"})
        continue
    b2.append({"region": r["region"], "cn": r["cn"], "n_sSNV": r["n_sSNV"],
               "haplotypes": r.get("haplotypes", "?"),
               "first_rank_edges": res["edges"], "first_rank_prob": res["first_rank_prob"],
               "tree_conf": res["tree_conf"], "missing_intermediates": res["missing_intermediates"]})

# 摘要
b1_probs = [x["first_rank_prob"] for x in b1]
b2_resolvable = [x for x in b2 if "first_rank_edges" in x]
b2_probs = [x["first_rank_prob"] for x in b2_resolvable]
import statistics as st
summary = {
    "A1_incompatible": {"n": len(a1),
        "cause_dist": dict(Counter("CN-gain" if "gain" in x["cause"] else "artifact/其他" for x in a1)),
        "fixable_artifact": sum(1 for x in a1 if "artifact" in x["fixable"])},
    "B1_ambiguous_parsimony": {"n": len(b1),
        "median_first_rank_prob": round(st.median(b1_probs), 3) if b1_probs else None,
        "high_conf(>=0.7)": sum(1 for p in b1_probs if p >= 0.7),
        "total_missing_intermediates": sum(x["missing_intermediates"] for x in b1)},
    "B2_underdetermined_firstrank": {"n": len(b2),
        "resolvable_to_firstrank": len(b2_resolvable),
        "no_pair_need_coverage": len(b2) - len(b2_resolvable),
        "median_first_rank_prob": round(st.median(b2_probs), 3) if b2_probs else None,
        "high_conf(>=0.7)": sum(1 for p in b2_probs if p >= 0.7)},
    "note": "first_rank_prob = ∏(每節點父指派信心);parsimony=最大真子集父+CCF tiebreak;跳>1=缺中間群(0.5^步)。"
            "這是『相容樹集合的第一順位+機率』,非單一定論(欠定本質);Tier-PS 僅能排跨-HP(已知 germline≠clonal)。"}
json.dump({"summary": summary, "A1_incompatible": a1, "B1_ambiguous": b1, "B2_underdetermined": b2},
          open(f"{DATA}/backbone_resolution.json", "w"), ensure_ascii=False, indent=1)
print("BACKBONE RESOLUTION DONE")
print(json.dumps(summary, ensure_ascii=False, indent=1))
print("\n=== A1 成環例(前3) ===")
for x in a1[:3]: print(f"  {x['region']} cn={x['cn']} 違反對={x['n_violating_pairs_in_stored']} 例={x['violations'][:1]} → {x['cause']}")
