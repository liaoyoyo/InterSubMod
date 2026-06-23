#!/usr/bin/env python3
"""[v4 統計聚合] records_v4 → 所有 SVG/數字+%/cross-tab/funnel/原始盤點/甲基 by 狀況 的注入資料。
數字全落 dashboard_stats_v4.json(generator 不手打)。純讀, 零 binary。"""
import json, statistics
from collections import Counter, defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v4.json"))
TP = [r for r in rec if r["set"] == "TP"]; FP = [r for r in rec if r["set"] == "FP"]
N = len(rec)
WSTATES = ["S1", "S2", "S3", "S4", "S6"]
CNSTATES = ["gain", "loss", "loh", "neutral"]
CHRS = [f"chr{c}" for c in range(1, 23)]


def cnt(items, key): return dict(Counter(r[key] for r in items))
def pct(p, w): return round(100 * p / w, 2) if w else 0
def med(vals):
    v = [x for x in vals if x is not None]
    return round(statistics.median(v), 4) if v else None
def structured(items): return sum(1 for r in items if r["coarse_ng"] >= 2)


S = {
    "n": N, "tp": len(TP), "fp": len(FP),
    "tp_pct": pct(len(TP), N), "fp_pct": pct(len(FP), N),
    # ---- ⓪ 原始數據盤點 ----
    "raw": {
        "tp_structured": structured(TP), "tp_structured_pct": pct(structured(TP), len(TP)),
        "fp_structured": structured(FP), "fp_structured_pct": pct(structured(FP), len(FP)),
        "tp_aligned": sum(1 for r in TP if r["align_status"] == "aligned"), "fp_aligned": sum(1 for r in FP if r["align_status"] == "aligned"),
        "median_n": int(statistics.median([r["n"] for r in rec])),
        "median_n_tp": int(statistics.median([r["n"] for r in TP])), "median_n_fp": int(statistics.median([r["n"] for r in FP])),
    },
    # per-chrom TP/FP/structure
    "per_chrom": {c: {"tp": sum(1 for r in TP if r["chrom"] == c), "fp": sum(1 for r in FP if r["chrom"] == c),
                      "struct": sum(1 for r in rec if r["chrom"] == c and r["coarse_ng"] >= 2)} for c in CHRS},
    # ---- 比例條 ----
    "cn_state": {s: cnt(rec, "cn_state").get(s, 0) for s in CNSTATES},
    "wstate": {s: cnt(rec, "wstate").get(s, 0) for s in WSTATES},
    "tier": {str(t): cnt(rec, "structure_tier").get(t, 0) for t in (1, 2, 3)},
    "align_status": cnt(rec, "align_status"),
    "tn_total": {"tumor": sum(r.get("n_tumor", 0) or 0 for r in rec), "normal": sum(r.get("n_normal", 0) or 0 for r in rec)},
    "median_nt": int(statistics.median([r["n_tumor"] for r in rec if r.get("n_tumor") is not None])),
    "median_nn": int(statistics.median([r["n_normal"] for r in rec if r.get("n_normal") is not None])),
    "coarse_dist": {str(k): v for k, v in sorted(cnt(rec, "coarse_ng").items())},
    "fine_dist": {str(k): v for k, v in sorted(cnt(rec, "fine_ng").items())},
    "geom_dist": {str(k): v for k, v in sorted((k, v) for k, v in cnt(rec, "geometry_ng").items() if k is not None)},
    "funnel": [
        {"stage": "全位點", "n": N},
        {"stage": "有結構(coarse≥2)", "n": sum(1 for r in rec if r["coarse_ng"] >= 2)},
        {"stage": "對齊(cis-ASM 跡象)", "n": sum(1 for r in rec if r["align_status"] == "aligned")},
        {"stage": "未對齊候選", "n": sum(1 for r in rec if r["align_status"] == "unaligned")},
        {"stage": "可疑漏切(幾何≥2嚴閘1)", "n": sum(r["geom_divergence"] for r in rec)},
    ],
    "confident_multi_tp_pct": pct(sum(1 for r in TP if r["wstate"] in ("S1", "S2")), len(TP)),
    "confident_multi_fp_pct": pct(sum(1 for r in FP if r["wstate"] in ("S1", "S2")), len(FP)),
    "geom_divergence": sum(r["geom_divergence"] for r in rec),
    # ---- cross-tabs(更多參數) ----
    "wstate_x_set": {w: {"tp": sum(1 for r in TP if r["wstate"] == w), "fp": sum(1 for r in FP if r["wstate"] == w)} for w in WSTATES},
    "cn_x_structure": {s: {"multi": sum(1 for r in rec if r["cn_state"] == s and r["coarse_ng"] >= 2),
                           "single": sum(1 for r in rec if r["cn_state"] == s and r["coarse_ng"] < 2)} for s in CNSTATES},
    "align_x_cn": {s: {a: sum(1 for r in rec if r["cn_state"] == s and r["align_status"] == a) for a in ("aligned", "unaligned", "NA")} for s in CNSTATES},
    # coarse × fine (over-split diagnostic) — 限制到 1..5+ 格
    "coarse_x_fine": {str(c): {str(f): sum(1 for r in rec if min(r["coarse_ng"], 5) == c and min(r["fine_ng"], 5) == f) for f in range(1, 6)} for c in range(1, 6)},
    # 各 flag rate × set
    "flag_x_set": {flag: {"tp": pct(sum(1 for r in TP if r[col]), len(TP)), "fp": pct(sum(1 for r in FP if r[col]), len(FP))}
                   for flag, col in [("unstable", "unstable"), ("hidden_het", "hidden_het"), ("geom_divergence", "geom_divergence")]},
    "other_rate": {"tp": pct(sum(1 for r in TP if r["n_other"] > 0), len(TP)), "fp": pct(sum(1 for r in FP if r["n_other"] > 0), len(FP))},
    # V 主導軸 × set (在有結構位點: hp 主導 / allele 主導 / 皆弱)
    "vaxis_x_set": {st: {"hp": sum(1 for r in items if r["coarse_ng"] >= 2 and r["V_hp"] >= r["V_allele"] and max(r["V_hp"], r["V_allele"]) >= 0.3),
                          "allele": sum(1 for r in items if r["coarse_ng"] >= 2 and r["V_allele"] > r["V_hp"] and max(r["V_hp"], r["V_allele"]) >= 0.3),
                          "weak": sum(1 for r in items if r["coarse_ng"] >= 2 and max(r["V_hp"], r["V_allele"]) < 0.3)}
                    for st, items in [("tp", TP), ("fp", FP)]},
    "reads_by_wstate": {w: int(statistics.median([r["n"] for r in rec if r["wstate"] == w] or [0])) for w in WSTATES},
    # ---- 甲基 aggregate ----
    "meth_overall": {
        "mean_beta_median": med([r["m_mean_beta"] for r in rec]),
        "dbeta_tn_median": med([r["m_dbeta_tn"] for r in rec]),
        "dbeta_group_median": med([r["m_dbeta_group"] for r in rec]),
        "n_with_group": sum(1 for r in rec if r.get("m_dbeta_group") is not None),
    },
    # 甲基覆蓋 / CpG(甲基位點)數量
    "meth_coverage": {
        "n_with_methylation": sum(1 for r in rec if r.get("m_mean_beta") is not None),
        "n_with_dbeta_group": sum(1 for r in rec if r.get("m_dbeta_group") is not None),
        "n_with_dbeta_tn": sum(1 for r in rec if r.get("m_dbeta_tn") is not None),
        "n_with_dbeta_hp": sum(1 for r in rec if r.get("m_dbeta_hp") is not None),
        "ncpg_total": sum(r["m_n_cpg"] for r in rec if r.get("m_n_cpg")),
        "ncpg_median": int(statistics.median([r["m_n_cpg"] for r in rec if r.get("m_n_cpg")])) if any(r.get("m_n_cpg") for r in rec) else 0,
        "ncpg_min": min((r["m_n_cpg"] for r in rec if r.get("m_n_cpg")), default=0),
        "ncpg_max": max((r["m_n_cpg"] for r in rec if r.get("m_n_cpg")), default=0),
    },
    # CpG(甲基位點)數/位點 直方
    "ncpg_hist": {lab: sum(1 for r in rec if r.get("m_n_cpg") and lo <= r["m_n_cpg"] < hi)
                  for lab, lo, hi in [("<20", 0, 20), ("20-50", 20, 50), ("50-100", 50, 100), ("100-200", 100, 200), ("200+", 200, 10**9)]},
    # mean β 直方(0-1, 10 bins)
    "mean_beta_hist": {f"{b/10:.1f}": sum(1 for r in rec if r.get("m_mean_beta") is not None and b/10 <= r["m_mean_beta"] < (b+1)/10 + (0.001 if b == 9 else 0)) for b in range(10)},
    # |dbeta_group| 直方(multi-group 位點; 0..0.3+, 6 bins)
    "dbeta_group_hist": {lab: sum(1 for r in rec if r.get("m_dbeta_group") is not None and lo <= r["m_dbeta_group"] < hi)
                         for lab, lo, hi in [("0-.05", 0, 0.05), (".05-.1", 0.05, 0.1), (".1-.15", 0.1, 0.15), (".15-.2", 0.15, 0.2), (".2-.3", 0.2, 0.3), (".3+", 0.3, 9)]},
    # ---- 甲基 by 狀況(median; 觀察「某些狀況的甲基」) ----
    "meth_by_wstate": {w: {"mean_beta": med([r["m_mean_beta"] for r in rec if r["wstate"] == w]),
                            "dbeta_group": med([r["m_dbeta_group"] for r in rec if r["wstate"] == w]),
                            "dbeta_tn": med([r["m_dbeta_tn"] for r in rec if r["wstate"] == w])} for w in WSTATES},
    "meth_by_cn": {s: {"mean_beta": med([r["m_mean_beta"] for r in rec if r["cn_state"] == s]),
                        "dbeta_group": med([r["m_dbeta_group"] for r in rec if r["cn_state"] == s]),
                        "dbeta_tn": med([r["m_dbeta_tn"] for r in rec if r["cn_state"] == s])} for s in CNSTATES},
    "meth_by_set": {st: {"mean_beta": med([r["m_mean_beta"] for r in items]),
                          "dbeta_group": med([r["m_dbeta_group"] for r in items]),
                          "dbeta_tn": med([r["m_dbeta_tn"] for r in items])} for st, items in [("TP", TP), ("FP", FP)]},
}
S["sumcheck_wstate"] = sum(cnt(rec, "wstate").values()) == N
S["sumcheck_perchrom"] = sum(v["tp"] + v["fp"] for v in S["per_chrom"].values()) == N
json.dump(S, open(f"{A}/dashboard_stats_v4.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps({k: S[k] for k in ("n", "tp", "fp", "raw", "meth_overall", "confident_multi_tp_pct", "confident_multi_fp_pct")}, ensure_ascii=False, indent=1))
print("meth_by_set:", json.dumps(S["meth_by_set"], ensure_ascii=False))
print("wstate_x_set:", json.dumps(S["wstate_x_set"], ensure_ascii=False))
print("sumcheck wstate/perchrom:", S["sumcheck_wstate"], S["sumcheck_perchrom"])
