#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CCF read-count 加權 — 全面觀察 v2(2026-07-09,含 partial reads)。
demo v1 用全覆蓋 populations(reach 下界 2.3%);本版用 mlhp col_coverage_by_hp(每位點 nALT,
含 partial reads=更多 read)算 per-mutation CCF,對 layered solver 的等機率樹(full+partial)加權。
全面分類觀察:ambiguity 型(order/structure)× CN(clean/gain)× tie-broken × winner-clean。
約束不變(RRR 根+unit-flip);純 read count(非甲基非循環)。只讀不寫源碼。partial flag: HCC1395。
用法: SM_ML_GLOB="<dir>/mlhp_part_*.json" python3 ccf_tree_weighting_full_observe.py
"""
import json, os, math, glob
from collections import Counter
import tree_enumeration_solver as S

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.normpath(os.path.join(HERE, "..", "data"))
GLOB = os.environ.get("SM_ML_GLOB", os.path.join(HERE, "..", "..", "20260618_subcluster_pilot", "mlhp_part_*.json"))
MARGIN = 0.05; TEMP = 0.05
GERM = ("1", "2", "3")

def ccf_from_colcov(colcov, positions, k):
    """col_coverage {pos:[nREF,nALT]} + positions → 每 bit CCF(用含 partial 的所有 read)。位點對不齊回 None。"""
    if not colcov or len(positions) < k: return None
    ccf = {}
    for j in range(k):
        pos = positions[j]
        v = colcov.get(str(pos)) or colcov.get(pos)
        if not v: return None
        nr, na = v[0], v[1]; tot = nr + na
        ccf[j] = na / tot if tot else 0.0
    return ccf

def _unlab(s):
    if s == "ROOT": return frozenset()
    s2 = s[2:] if s.startswith("H_") else s
    return frozenset(i for i, c in enumerate(s2) if c == "A")

def score_viol(edges_lab, ccf):
    sc = 0.0; vi = 0
    for (p, c) in edges_lab:
        P, C = _unlab(p), _unlab(c); acq = C - P
        if len(acq) != 1: continue
        j = next(iter(acq))
        for anc in P:
            sc += (ccf[anc] - ccf[j])
            if ccf[j] > ccf[anc] + MARGIN: vi += 1
    return sc, vi

def main():
    groups = []
    for f in sorted(glob.glob(GLOB)):
        groups += json.load(open(f)).get("groups", [])
    regs = [g for g in groups if g.get("n_sSNV", 0) >= 2]
    # 分類累計:key=(amb_type, cn_class, broke) → count
    n_amb = 0; broke = 0; stay = 0; win_clean = 0
    have_ccf = 0; no_ccf = 0
    by = Counter()          # (amb_type, cn_class) → total
    by_broke = Counter()    # (amb_type, cn_class) → broke
    top_dist = Counter()
    examples = []
    for r in regs:
        pbh = r.get("populations_by_hp", {}) or {}; sbh = r.get("subread_groups_by_hp", {}) or {}
        cbh = r.get("col_coverage_by_hp", {}) or {}; positions = r.get("positions", []); cn = r.get("cn")
        for fam in sorted(set(pbh) | set(sbh)):
            if fam not in GERM: continue
            full = pbh.get(fam, {}) or {}; part = list((sbh.get(fam, {}) or {}).keys())
            if not full and not part: continue
            k = len(next(iter(full))) if full else len(part[0])
            if k < 2 or k > 8: continue
            res = S.enumerate_min_trees(full, part, k)
            cls = S.classify(res)
            if res["capped"] or "ambiguous" not in cls or res["n_trees"] < 2 or len(res["trees"]) < 2:
                continue
            n_amb += 1
            amb_type = "order(同節點集)" if len(res.get("_feasible_N", [])) == 1 else "structure(多完成)"
            cn_class = "gain(read≠CCF)" if cn == "gain" else ("clean" if cn in ("neutral", "loh", "loss") else "unknown")
            ccf = ccf_from_colcov(cbh.get(fam, {}), positions, k)
            if ccf is None:
                # fallback: populations CCF(全覆蓋)
                tot = sum(full.values()) or 1
                ccf = {j: sum(n for g, n in full.items() if g[j] == "A") / tot for j in range(k)}
                no_ccf += 1
            else:
                have_ccf += 1
            scored = [score_viol(t["edges"], ccf) + (t,) for t in res["trees"]]
            mx = max(s for s, _, _ in scored)
            exps = [math.exp((s - mx) / TEMP) for s, _, _ in scored]; tot = sum(exps) or 1
            post = [e / tot for e in exps]; top = max(post)
            top_dist[("≥0.9" if top >= 0.9 else "0.6-0.9" if top >= 0.6 else "0.4-0.6" if top >= 0.4 else "<0.4")] += 1
            by[(amb_type, cn_class)] += 1
            broke_this = top >= 0.6
            if broke_this:
                broke += 1; by_broke[(amb_type, cn_class)] += 1
                if scored[post.index(top)][1] == 0: win_clean += 1
            else:
                stay += 1
            if broke_this and cn_class == "clean" and len(examples) < 6:
                examples.append({"region": f"{r['chrom']}:{r['start']}-{r['end']}", "fam": fam, "cn": cn,
                                 "amb_type": amb_type, "n_trees": res["n_trees"], "top_post": round(top, 3),
                                 "ccf": {j: round(v, 2) for j, v in ccf.items()}})
    out = {
        "partial_flag": "HCC1395 / layered(full+partial) ambiguous 子集 / col_coverage CCF(含partial reads)",
        "n_ambiguous_units": n_amb, "ccf_from_colcov": have_ccf, "ccf_fallback_populations": no_ccf,
        "tie_broken(≥0.6)": broke, "broke_frac": round(broke / n_amb, 3) if n_amb else None,
        "stay_uniform": stay, "winner_pigeonhole_clean": win_clean,
        "winner_clean_frac_of_broke": round(win_clean / broke, 3) if broke else None,
        "top_posterior_dist": dict(top_dist),
        "by_type_cn_total": {f"{a}|{c}": v for (a, c), v in sorted(by.items())},
        "by_type_cn_broke": {f"{a}|{c}": v for (a, c), v in sorted(by_broke.items())},
        "examples(clean,broke)": examples,
    }
    json.dump({"summary": out}, open(os.path.join(D, "ccf_tree_weighting_full_observe.json"), "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))

if __name__ == "__main__":
    main()
