#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
匯總 T-GATE-GB cis-control pilot + 全掃描的所有報告數字到單一 JSON(§13-A 注入,確保 grep-able)。
不做新 compute,只讀已落檔 JSON 重算衍生數字(聯合分布/scope/needs_methyl 對接)。
輸出:cis_control_scope_summary.json。報告只引此 JSON 的數字。
"""
import json, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

pilot = json.load(open(f"{DATA}/pilot_cis_control.json"))
axis = json.load(open(f"{DATA}/pilot_axis_alignment.json"))
neu = json.load(open(f"{DATA}/axis_alignment_neutral.json"))
fs = json.load(open(f"{DATA}/hp_alignment_fullscan.json"))
cs = json.load(open(f"{DATA}/candidate_scoring.json"))

# (1) pilot 聯合分布
P = pilot["all_cpg_pairs"]
sig = [p for p in P if abs(p["tumor_dbeta"]) >= 0.2]
germ = [p for p in sig if abs(p["normal_dbeta"]) >= 0.2]
mid = [p for p in sig if 0.1 <= abs(p["normal_dbeta"]) < 0.2]
cand = [p for p in sig if abs(p["normal_dbeta"]) < 0.1]
joint = {"n_cpg_total": len(P), "n_tumor_sig": len(sig), "tumor_sig_pct": round(100*len(sig)/len(P), 1),
         "of_sig_germline_asm": len(germ), "of_sig_germline_asm_pct": round(100*len(germ)/len(sig)),
         "of_sig_graybard": len(mid), "of_sig_graybard_pct": round(100*len(mid)/len(sig)),
         "of_sig_subclone_candidate_if_axes_aligned": len(cand), "of_sig_candidate_pct": round(100*len(cand)/len(sig)),
         "corr_tumor_vs_normal_dbeta": axis["corr_tumor_vs_normal_dbeta"]}

# (2) 全掃描 scope
N = fs["n_classified"]
ov = fs["overall"]
scope = {"n_classified": N, "n_targets": fs["n_targets"], "n_skipped": fs["n_skipped"],
         "CROSS_HP": ov.get("CROSS-HP", 0), "CROSS_HP_pct": round(100*ov.get("CROSS-HP", 0)/N, 1),
         "SAME_HP": ov.get("SAME-HP", 0), "SAME_HP_pct": round(100*ov.get("SAME-HP", 0)/N, 1),
         "MIXED": ov.get("MIXED", 0), "MIXED_pct": round(100*ov.get("MIXED", 0)/N, 1),
         "by_cn": fs["by_cn"]}

# (3) 最乾淨適用集(neutral/loss + CROSS-HP)
clean = [r["region"] for r in fs["regions"] if r["cn"] in ("neutral", "loss") and r["hp_alignment"] == "CROSS-HP"]

# (4) needs_methyl 對接
align = {r["region"]: r["hp_alignment"] for r in fs["regions"]}
acn = {r["region"]: r["cn"] for r in fs["regions"]}
nm = [q["region"] for q in cs["queue"] if q.get("needs_methyl")]
nm_cl = [r for r in nm if r in align]
nm_dist = dict(Counter(align[r] for r in nm_cl))
nm_cross = [r for r in nm_cl if align[r] == "CROSS-HP"]
nm_cross_neu = [r for r in nm_cross if acn[r] == "neutral"]
needs_methyl = {"n_needs_methyl_total": len(nm), "n_with_ge2_altpop_classified": len(nm_cl),
                "hp_alignment_dist": nm_dist, "cn_dist": dict(Counter(acn[r] for r in nm_cl)),
                "n_cross_hp_cis_applicable": len(nm_cross), "n_cross_hp_AND_neutral_clean": len(nm_cross_neu)}

# (5) pilot 區清單(全 LOH 證據)
pilot_regions = [{"region": r["region"], "cn": r["cn"], "n_cpg_both": r.get("n_cpg_both", 0)}
                 for r in pilot["regions"] if r.get("n_cpg_both")]

out = {"pilot_summary": pilot["summary"], "pilot_regions_with_cpg": pilot_regions,
       "joint_distribution": joint, "fullscan_scope": scope,
       "cleanest_applicable_set": {"n": len(clean), "definition": "cn in (neutral,loss) AND CROSS-HP", "examples": clean[:8]},
       "needs_methyl_intersection": needs_methyl,
       "neutral_pilot_axis": {"n": neu["n_neutral_checked"], "dist": neu["verdict_dist"]},
       "verdict": "matched-normal HP cis-control 僅對 CROSS-HP 區有效(全資料 35.4%);needs_methyl ∩ 乾淨可用 ≈ 0。甲基維持 bounded-auxiliary;subclone-specificity 對 SAME-HP 多數區為 structural UNDETERMINED。"}
json.dump(out, open(f"{DATA}/cis_control_scope_summary.json", "w"), ensure_ascii=False, indent=1)
print("SUMMARY WRITTEN")
print(json.dumps({"joint": joint, "scope": {k: scope[k] for k in ("CROSS_HP_pct","SAME_HP_pct","MIXED_pct")},
                  "clean_set": len(clean), "needs_methyl": needs_methyl}, ensure_ascii=False, indent=1))
