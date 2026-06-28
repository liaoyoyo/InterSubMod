#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
骨幹結構穩定性審查(落實完整性 workflow 的修正,數字 grep-able §13-A):
- 雙 denominator: detail 3885(有 genotype 向量) vs stats 7143(全區)。
- robustness ladder: L1 A_determined → +TP-backed → +n_sSNV≥3 → +non-CN-gain(真正穩固核心)。
- verified incompatible: 要求 npop≥3 AND stored 4-gamete 違反對(現 12 多為 artifact)。
- core composition: %n_sSNV==2 / %CN-gain / %無 truth label。
- determinacy rate by CN(查 CN-gain inflation)。
- 截斷: n_sSNV>8(=42,非43)。
輸出 backbone_stability_audit.json。compute batch(§13.0,純 JSON)。
"""
import json, os
from itertools import combinations
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
T = json.load(open(f"{DATA}/topology_per_region.json"))
det = T["detail"]; stats = T.get("stats", {})
DET = lambda r: r.get("determinacy", "?")
isdet = lambda r: "A_determined" in DET(r)

# 雙 denominator(G1 修後:stats.determinacy 已 canonical=3885;7143 全區看 region_coverage)
n_detail = len(det)
cov = stats.get("region_coverage", {})
n_full = sum(cov.values()) if cov else None
stats_det_sum = sum(stats.get("determinacy", {}).values())

# robustness ladder
core = [r for r in det if isdet(r)]
L1 = len(core)
L2 = [r for r in core if r["fp"] == 0 and r["tp"] > 0]
L3 = [r for r in L2 if r["n_sSNV"] >= 3]
L4 = [r for r in L3 if r["cn"] != "gain"]

# incompatible 真相(G2 修):has_cycle 是上游 sm_region_integration 在「全 sSNV pairwise 關係圖」算的真 cycle,
# 非 stored 截斷 populations 可重現。stored 查 4-gamete = 截斷後查不到(非「0 真 cycle」)。
def four_gamete(pops):
    if not pops: return 0
    k = len(next(iter(pops))); v = 0
    for i, j in combinations(range(k), 2):
        gam = set(g[i] + g[j] for g, c in pops.items() if c > 0 and g[i] in "RA" and g[j] in "RA")
        if {"RR", "RA", "AR", "AA"} <= gam: v += 1
    return v
incomp = [r for r in det if DET(r) == "incompatible"]
verified_stored = [r for r in incomp if len(r["populations"]) >= 3 and four_gamete(r["populations"]) > 0]
RI = json.load(open(f"{DATA}/sm_region_integration.json"))
has_cycle_all = [r for r in RI["regions"] if r.get("has_cycle")]
hc_cn = Counter(r.get("cn") for r in has_cycle_all)
hc_gain = sum(1 for r in has_cycle_all if r.get("cn") == "gain")

# core composition
def pct(sub, tot): return round(100 * len(sub) / tot, 1) if tot else 0
core_n2 = [r for r in core if r["n_sSNV"] == 2]
core_gain = [r for r in core if r["cn"] == "gain"]
core_notruth = [r for r in core if r["tp"] == 0 and r["fp"] == 0]

# determinacy rate by CN
bycn = {}
for cn in ("gain", "loh", "neutral", "loss"):
    sub = [r for r in det if r["cn"] == cn]
    d = [r for r in sub if isdet(r)]
    bycn[cn] = {"n": len(sub), "determined": len(d), "rate": round(100 * len(d) / len(sub), 1) if sub else 0}

# 截斷
trunc = [r for r in det if r["n_sSNV"] > 8]
trunc_tp = [r for r in trunc if r["tp"] > 0]
trunc_incomp = [r for r in trunc if DET(r) == "incompatible"]

out = {
    "denominators": {"detail_somatic_vector": n_detail, "all_regions(region_coverage)": n_full,
        "region_coverage": dict(cov), "stats_determinacy_sum(now canonical)": stats_det_sum,
        "consistent(stats==detail)": stats_det_sum == n_detail,
        "note": "G1 修(2026-06-29):stats.determinacy 現 canonical = 有向量 3885(== detail,sum 一致);"
                "全 7143 區的覆蓋看 region_coverage(with_vector 3885 / germline_only 371 / no_vector 2887)。"
                "determinacy % 一律以 3885(有向量)為分母;對外提全區時用 region_coverage,勿再混引 7143-determinacy。"
                "47% determined = 1812/3885(canonical) = 25.4%/7143(全區覆蓋脈絡)。"},
    "robustness_ladder": {
        "L1_A_determined": {"n": L1, "pct_of_3885": pct(core, n_detail)},
        "L2_+TP_backed(fp=0,tp>0)": {"n": len(L2), "pct_of_3885": pct(L2, n_detail)},
        "L3_+multisite(n_sSNV>=3)": {"n": len(L3), "pct_of_3885": pct(L3, n_detail)},
        "L4_+non_CN_gain(truly_robust)": {"n": len(L4), "pct_of_3885": pct(L4, n_detail)},
        "interpretation": "L1 47% 是軟上界;同時非循環+TP-backed+多位點+非CN-gain 的真正穩固核心 = L4。"},
    "determined_core_fragility": {"n_core": L1,
        "pct_n_sSNV==2(僅2位點非多節點樹)": pct(core_n2, L1),
        "pct_CN_gain(multiplicity artifact 區)": pct(core_gain, L1),
        "pct_no_truth_label": pct(core_notruth, L1)},
    "incompatible_recompute": {
        "has_cycle_all_regions(upstream pairwise graph)": len(has_cycle_all),
        "with_genotype_vector(in detail)": len(incomp),
        "no_vector(not in detail)": len(has_cycle_all) - len(incomp),
        "cn_dist": dict(hc_cn), "cn_gain_frac": round(hc_gain / len(has_cycle_all), 2) if has_cycle_all else 0,
        "verifiable_from_stored_capped": len(verified_stored),
        "note": "has_cycle = 上游在全 sSNV pairwise 關係圖(nested_edges+sibling_pairs)算的真 cycle(如雙向 nesting),"
                "非 stored 截斷(cap=8)populations 可重現 → verifiable_from_stored=0 是『截斷查不到』非『0 真 cycle』。"
                "22 真 cycle 中 77% CN-gain → 最可能 = CN multiplicity(同突變多拷貝)製造 pairwise 矛盾,非真演化樹衝突。"
                "處置:標 likely-CN-multiplicity-artifact、不建樹;若要真驗證須提高 genotype cap + CN-aware multiplicity 模型。"},
    "determinacy_rate_by_cn": bycn,
    "cn_gain_inflation_note": "determinacy RATE 在 CN-gain 最高 → multiplicity 可能製造假共現→假 determined;tree-level claim 前應 mask/控制 CN-gain。",
    "genotype_truncation": {"n_truncated(n_sSNV>8)": len(trunc), "with_tp": len(trunc_tp),
        "incompatible_among_trunc": len(trunc_incomp),
        "note": "正解=42(n_sSNV>8),非43(43=n_sSNV>=8,含 1 個 n_sSNV==8 剛好不截斷)。core 中僅 3 區受影響→主結論穩。"},
}
json.dump(out, open(f"{DATA}/backbone_stability_audit.json", "w"), ensure_ascii=False, indent=1)
print("STABILITY AUDIT DONE")
print(json.dumps(out, ensure_ascii=False, indent=1))
