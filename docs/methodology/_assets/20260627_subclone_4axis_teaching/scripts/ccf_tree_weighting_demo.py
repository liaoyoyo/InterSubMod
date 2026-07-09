#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CCF read-count 加權 demo(2026-07-09,使用者提案):把 enumerate_min_trees 的「N 棵等機率樹」
用每群 read 數(→per-mutation CCF)加權 → 加權 posterior + 較可能者。
約束不變(RRR 根 + 有向 unit-flip);只用 read count(遺傳、非甲基、非循環)。
原理:pigeonhole/crossing rule — 祖先突變 CCF ≥ 後代突變 CCF(早發生=頻率高)。
  tree score = Σ_edge Σ_(anc∈parent)(CCF[anc] − CCF[acquired]);softmax→posterior。
  violation = 後代 CCF > 祖先 CCF + margin(明確違反 pigeonhole)。
🔴 CN-gain 區 read≠CCF(multiplicity)→ 標記、不當有效訊號。只讀不寫源碼。partial flag: HCC1395。
用法: python3 ccf_tree_weighting_demo.py
"""
import json, os, math
from itertools import combinations
from collections import Counter
import tree_enumeration_solver as S

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.normpath(os.path.join(HERE, "..", "data"))
MARGIN = 0.05   # CCF 差 < margin 視為等頻率(不可破 tie)
TEMP = 0.05     # softmax 溫度(CCF 差的尺度)

def mutation_ccf(pops, k):
    tot = sum(pops.values()) or 1
    return {j: sum(n for g, n in pops.items() if g[j] == "A") / tot for j in range(k)}

def _unlab(s):
    if s == "ROOT": return frozenset()
    s2 = s[2:] if s.startswith("H_") else s
    return frozenset(i for i, c in enumerate(s2) if c == "A")

def tree_score_and_viol(edges_lab, ccf):
    """回 (score, n_violation)。score 高=CCF-consistent(祖先突變頻率高);violation=後代>祖先+margin。"""
    score = 0.0; viol = 0
    for (p, c) in edges_lab:
        P, C = _unlab(p), _unlab(c)
        acq = C - P
        if len(acq) != 1: continue
        j = next(iter(acq))
        for anc in P:
            score += (ccf[anc] - ccf[j])
            if ccf[j] > ccf[anc] + MARGIN: viol += 1
    return score, viol

def main():
    R = json.load(open(os.path.join(D, "sm_region_integration.json")))["regions"]
    reg = [r for r in R if r.get("n_sSNV", 0) >= 2 and (r.get("populations") and len(r["populations"]) >= 1)]
    n_amb = 0; broke = 0; stay = 0; winner_clean = 0; cn_gain_amb = 0
    examples = []
    dist_top = Counter()
    for r in reg:
        pops = r["populations"]; k = len(next(iter(pops)))
        if k < 2 or k > 8: continue
        res = S.enumerate_min_trees(pops, [], k)
        if res["capped"] or res["n_trees"] < 2:
            continue
        # 只看有 stored trees(枚舉出的候選)
        trees = res["trees"]
        if len(trees) < 2: continue
        n_amb += 1
        cn = r.get("cn")
        ccf = mutation_ccf(pops, k)
        scored = []
        for t in trees:
            sc, vi = tree_score_and_viol(t["edges"], ccf)
            scored.append((sc, vi, t))
        mx = max(s for s, _, _ in scored)
        exps = [math.exp((s - mx) / TEMP) for s, _, _ in scored]
        tot = sum(exps) or 1
        post = [e / tot for e in exps]
        top = max(post)
        dist_top[("≥0.9" if top >= 0.9 else "0.6-0.9" if top >= 0.6 else "0.4-0.6" if top >= 0.4 else "<0.4")] += 1
        win_idx = post.index(top)
        win_viol = scored[win_idx][1]
        if cn == "gain": cn_gain_amb += 1
        if top >= 0.6:                       # 明確較可能者(tie 被破)
            broke += 1
            if win_viol == 0: winner_clean += 1
        else:
            stay += 1
        if len(examples) < 6 and top >= 0.6 and cn in ("neutral", "loh"):
            examples.append((r["region"], cn, res["n_trees"], round(top, 3), win_viol,
                             {j: round(v, 2) for j, v in ccf.items()}))
    out = {
        "partial_flag": "HCC1395 / ambiguous(n_trees≥2,非capped) 子集",
        "n_ambiguous_regions": n_amb,
        "tie_broken(top posterior≥0.6)": broke, "broke_frac": round(broke / n_amb, 3) if n_amb else None,
        "stay_uniform(<0.6)": stay,
        "winner_pigeonhole_clean(0 violation)": winner_clean,
        "winner_clean_frac_of_broke": round(winner_clean / broke, 3) if broke else None,
        "cn_gain_ambiguous(read≠CCF 需標記)": cn_gain_amb,
        "top_posterior_dist": dict(dist_top),
        "examples(CN-clean,tie 破)": examples,
    }
    json.dump({"summary": out}, open(os.path.join(D, "ccf_tree_weighting_demo.json"), "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))

if __name__ == "__main__":
    main()
