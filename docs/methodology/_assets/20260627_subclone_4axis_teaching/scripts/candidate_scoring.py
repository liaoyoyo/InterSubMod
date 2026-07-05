#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
候選拓樸評分 + 確認佇列：找出有衝突/資訊不確定/可能有 2-3 順位組合的區，
用可驗證資訊(determinacy/coverage/CN/VAF-orderable/HP/genome_ctx/FP)算 confidence 分數(0-100)，
並標 resolution_path(需甲基/需覆蓋/需Tier-PS/genetic-足夠)。供左右選項確認。
來源：topology_per_region.json。輸出：candidate_scoring.json。
用法：python3 candidate_scoring.py
"""
import json, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(HERE, "..", "data")))  # per-sample 化
d = json.load(open(os.path.join(DATA, "topology_per_region.json"), encoding="utf-8"))

BASE = {"A_determined(單分子向量)": 80, "A_ambiguous_order(缺中間群)": 55,
        "B_pairwise_structure": 50, "C_underdetermined": 30,
        # 2026-07-04 gap#3: Model A m-通道拆分後的 4 個 recurrence 子類
        "recurrence_candidate(m=1)": 45,              # diploid/single-copy → 真 recurrence 候選(留)
        "recurrence_LOH_unresolved": 30,              # copy-neutral LOH → m 未定(VAF L3 軟旗標)
        "recurrence_required(m通道不可用)": 30,        # 6 樣本 cn=unknown → 無 CN m 通道
        "recurrence_artifact(m>1;CN-amp)": 15,        # 整數 CN≥3 → 多重度 artifact(棄,同 incompatible)
        # 2026-07-04 gap#1: partial-read subcube 救回(無 full-cov single-molecule → 不僭稱 A_determined)
        "E_subcube_recovered(gap#1;pairwise樹)": 42, "E_subcube_recovered(gap#1;欠定)": 28,
        "E_subcube_recovered(gap#1;含衝突)": 18,
        "incompatible": 15, "other": 40}

def score_region(r):
    det = r["determinacy"]; base = BASE.get(det, 40)
    reads = sum(r["populations"].values()) if r["populations"] else 0
    s = base
    s += min(15, reads / 2)                              # 覆蓋加分
    s += 10 if r["cn"] in ("loh", "neutral") else 0      # CN-clean(VAF 可信)
    s += 5 if r["haplotypes"] in ("H1", "H2") else 0     # 單根 HP 乾淨
    s -= 10 if r["fp"] > r["tp"] else 0                  # FP 主導扣分
    s -= min(20, 10 * r["ambig_nodes"])                  # 順序未定扣分
    s -= {"centromere": 15, "telomere": 5}.get(r["genome_ctx"], 0)  # 偽影區扣分
    return max(0, min(100, round(s)))

def resolution(r):
    if r["determinacy"] == "incompatible" or (r["genome_ctx"] == "centromere" and r["fp"] > 0):
        return "likely-artifact(成環/著絲點+FP)→補 mappability mask、不強建樹"
    _d = r["determinacy"]
    if "recurrence_artifact" in _d:
        return "整數 CN≥3(m>1)→ CN-multiplicity artifact,棄(粗 categorical gate 漏掉的 GAINLOH)"
    if "recurrence_candidate" in _d:
        return "diploid/single-copy(m=1)→ Model A 真 recurrence 候選,留(需 orthogonal 佐證確認,非甲基)"
    if "recurrence_LOH" in _d:
        return "copy-neutral LOH total=2、allele-specific CN 無 → m∈{1,2} 未解;VAF L3 軟旗標非硬判"
    if "recurrence" in _d:
        return "cn=unknown(非HCC1395)→ CN m 通道不可用,recurrence 未拆"
    if r["determinacy"] == "A_ambiguous_order(缺中間群)":
        # D4 修(2026-07-01):移除「VAF tie→甲基輔助」overclaim。甲基排序為 L3-weak(ρ≈0.18 未顯著,
        # ordering pilot)、非可信 resolver,且 backbone_resolution compute 中甲基完全缺席 → 不宣稱甲基定序。
        return "順序未定→VAF/CCF(CN-clean)定先後;VAF tie→保持未定(甲基僅 L3-weak 軟旗標、非排序 resolver)"
    if r["determinacy"] == "C_underdetermined":
        return "資訊不足→加深覆蓋觀察缺的對 / Tier-PS 連遠 partner"
    if r["determinacy"] == "B_pairwise_structure":
        return "pairwise 拼接非單分子→需更長 read/單分子整跨確認"
    return "genetic 足夠(單分子向量確定)"

def situation(r):
    d = r["determinacy"]
    if "E_subcube_recovered" in d: return "subcube救回(partial建樹)"   # gap#1: 無 full-cov,partial+pairwise 建
    if d == "incompatible": return "衝突(成環)"
    if "recurrence_artifact" in d: return "recurrence→artifact(m>1棄)"      # gap#3: 整數 CN≥3 多重度
    if "recurrence_candidate" in d: return "recurrence候選(m=1留)"          # diploid → 真 recurrence 候選
    if "recurrence_LOH" in d: return "recurrence·LOH(m未定)"                # copy-neutral LOH,VAF L3
    if "recurrence" in d: return "recurrence(m通道不可用)"                   # 6 樣本 cn=unknown
    if r["ambig_nodes"] > 0: return "順序 2-3 順位待定"
    if r["determinacy"] == "C_underdetermined": return "多樹相容(欠定)"
    if r["haplotypes"] not in ("H1", "H2", "?"): return "跨HP(兩棵樹)"
    if r["determinacy"] == "B_pairwise_structure": return "pairwise 拼接"
    return "已確定"

# 加 3 欄資料源(2026-06-29):backbone_resolution(parsimony 第一順位) + detail cycle_cause(why-conflict)
import os as _os
_br_path = _os.path.join(DATA, "backbone_resolution.json")
prob1 = {}
if _os.path.exists(_br_path):
    _br = json.load(open(_br_path, encoding="utf-8"))
    for x in _br.get("B1_ambiguous", []): prob1[x["region"]] = x.get("first_rank_prob")
    for x in _br.get("B2_underdetermined", []):
        if x.get("first_rank_prob") is not None: prob1[x["region"]] = x["first_rank_prob"]

def why_conflict(r):
    """這區為何無法乾淨定樹(逐 situation 解釋)。"""
    det = r["determinacy"]
    if det == "incompatible":
        return r.get("cycle_cause") or "pairwise 成環(上游 has_cycle)"
    if "recurrence" in det:
        mv = (r.get("m_channel") or {}).get("verdict", "?")
        return f"非-gain 四型違反(某位元須翻≥2次)→Model A recurrence;m-通道判={mv}"
    if r["ambig_nodes"] > 0:
        return "缺中間群(跳>1突變→順序未定)" + ("+截斷" if r.get("truncated") else "")
    if det == "C_underdetermined":
        return "只觀測到單一 ALT 群(缺連接 read)→需深覆蓋"
    if r["haplotypes"] not in ("H1", "H2", "?"):
        return "跨 germline HP(已拆兩棵獨立樹)"
    if det == "B_pairwise_structure":
        return "pairwise 拼接(無單分子整跨)"
    return "已確定(單分子向量唯一)"

def methyl_applic(r):
    """甲基對這區能否幫(基於 bounded-auxiliary 全測結論)。"""
    if "3" in r.get("haplotypes", ""):  # 有 H3-unphased
        return "HP定相: germline-ASM 在→可定相 H3(否則不可);排序 L3 弱;specificity 不可"
    return "負篩可用(確認無次結構);排序 L3 弱;specificity/定群 不可"

queue = []
for r in d["detail"]:
    sc = score_region(r); sit = situation(r)
    # 候選數估計: 順序未定→n_clusters! 級(capped 顯示);欠定→「≥2」;確定→1
    if "recurrence" in r["determinacy"]: ncand = "recurrence(m-通道:%s)" % ((r.get("m_channel") or {}).get("verdict", "?"))
    elif r["ambig_nodes"] > 0: ncand = f"≥{r['ambig_nodes']+1}(順序排列)"
    elif r["determinacy"] == "C_underdetermined": ncand = "≥2(缺對)"
    elif r["determinacy"] == "incompatible": ncand = "0(不支持單樹)"
    else: ncand = "1(唯一)"
    rec = {"region": r["region"], "chrom": r["chrom"], "start": r["start"], "n_sSNV": r["n_sSNV"],
           "n_clusters": r["n_clusters"], "topology_type": r["topology_type"], "haplotypes": r["haplotypes"],
           "cn": r["cn"], "genome_ctx": r["genome_ctx"], "tp": r["tp"], "fp": r["fp"],
           "determinacy": r["determinacy"], "ambig_nodes": r["ambig_nodes"],
           "confidence_score": sc, "situation": sit, "n_candidates": ncand,
           "resolution_path": resolution(r), "needs_methyl": "VAF tie" in resolution(r) or r["determinacy"] == "C_underdetermined",
           "why_conflict": why_conflict(r),                       # 新欄1: 為何無法乾淨定樹
           "parsimony_first_rank_prob": prob1.get(r["region"]),    # 新欄2: parsimony 第一順位機率(ambiguous/underdetermined)
           "methyl_applicability": methyl_applic(r),               # 新欄3: 甲基能否幫(bounded)
           "truncated": r.get("truncated", False), "m_channel": r.get("m_channel"),  # gap#3 m-通道拆分
           "edges": r["edges"], "populations": r["populations"]}
    queue.append(rec)

# 需確認佇列 = 非「已確定 + 高分」的
need = [q for q in queue if q["situation"] != "已確定" or q["confidence_score"] < 70]
need.sort(key=lambda q: q["confidence_score"])  # 最需關注(低分)在前

out = {"n_total": len(queue), "n_need_confirm": len(need),
       "score_formula": "base(determinacy) + 覆蓋 + CN-clean + 單HP − FP主導 − 順序未定 − 偽影區; clamp 0-100",
       "situation_dist": dict(Counter(q["situation"] for q in queue)),
       "resolution_dist": dict(Counter(q["resolution_path"].split("→")[0] for q in queue)),
       "score_buckets": dict(Counter(("≥80高" if q["confidence_score"]>=80 else "60-79中" if q["confidence_score"]>=60 else "40-59低" if q["confidence_score"]>=40 else "<40極低") for q in queue)),
       "needs_methyl_n": sum(1 for q in queue if q["needs_methyl"]),
       "queue": need}
with open(os.path.join(DATA, "candidate_scoring.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False)

print(f"總 {len(queue)} 區;需確認佇列 {len(need)} 區")
print("situation 分布:", out["situation_dist"])
print("score 分布:", out["score_buckets"])
print("resolution 分布:", out["resolution_dist"])
print("需甲基輔助(VAF tie/欠定):", out["needs_methyl_n"])
print("佇列前 5(最低分):")
for q in need[:5]: print(f"  {q['region']} score={q['confidence_score']} [{q['situation']}] cand={q['n_candidates']} → {q['resolution_path'][:40]}")
print("OK wrote candidate_scoring.json")
