#!/usr/bin/env python3
"""[v3 dashboard 統計聚合] records_v3 → 所有 SVG 比例圖/cross-tab/funnel 的注入資料(數字落 JSON, generator 不手打)。
輸出 dashboard_stats_v3.json。純讀, 零 compute。"""
import json, statistics
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v3.json"))
TP = [r for r in rec if r["set"] == "TP"]; FP = [r for r in rec if r["set"] == "FP"]
N = len(rec)


def cnt(items, key):
    return dict(Counter(r[key] for r in items))


def pct(part, whole):
    return round(100 * part / whole, 2) if whole else 0


WSTATES = ["S1", "S2", "S3", "S4", "S6"]
CNSTATES = ["gain", "loss", "loh", "neutral"]

S = {
    "n": N, "tp": len(TP), "fp": len(FP),
    # 比例條: CN 狀態 / structure 5 態 / structure tier
    "cn_state": {s: cnt(rec, "cn_state").get(s, 0) for s in CNSTATES},
    "wstate": {s: cnt(rec, "wstate").get(s, 0) for s in WSTATES},
    "tier": {str(t): cnt(rec, "structure_tier").get(t, 0) for t in (1, 2, 3)},
    "align_status": cnt(rec, "align_status"),
    # T/N 組成(總) — 補充圖
    "tn_total": {"tumor": sum(r.get("n_tumor", 0) or 0 for r in rec), "normal": sum(r.get("n_normal", 0) or 0 for r in rec)},
    "median_nt": int(statistics.median([r["n_tumor"] for r in rec if r.get("n_tumor") is not None])),
    "median_nn": int(statistics.median([r["n_normal"] for r in rec if r.get("n_normal") is not None])),
    # 切群數分佈(直方): coarse(null95) / fine(null90) / geometry
    "coarse_dist": {str(k): v for k, v in sorted(cnt(rec, "coarse_ng").items())},
    "fine_dist": {str(k): v for k, v in sorted(cnt(rec, "fine_ng").items())},
    "geom_dist": {str(k): v for k, v in sorted((k, v) for k, v in cnt(rec, "geometry_ng").items() if k is not None)},
    # 判別 funnel
    "funnel": [
        {"stage": "全位點", "n": N},
        {"stage": "有結構(coarse≥2)", "n": sum(1 for r in rec if r["coarse_ng"] >= 2)},
        {"stage": "對齊(cis-ASM 跡象)", "n": sum(1 for r in rec if r["align_status"] == "aligned")},
        {"stage": "未對齊候選", "n": sum(1 for r in rec if r["align_status"] == "unaligned")},
        {"stage": "可疑漏切(幾何≥2嚴閘1)", "n": sum(r["geom_divergence"] for r in rec)},
    ],
    # cross-tab: CN×structure(multi/single)
    "cn_x_structure": {s: {"multi": sum(1 for r in rec if r["cn_state"] == s and r["coarse_ng"] >= 2),
                           "single": sum(1 for r in rec if r["cn_state"] == s and r["coarse_ng"] < 2)} for s in CNSTATES},
    # cross-tab: align×CN
    "align_x_cn": {s: {a: sum(1 for r in rec if r["cn_state"] == s and r["align_status"] == a) for a in ("aligned", "unaligned", "NA")} for s in CNSTATES},
    # reads × wstate (中位 n)
    "reads_by_wstate": {w: int(statistics.median([r["n"] for r in rec if r["wstate"] == w] or [0])) for w in WSTATES},
    # TP vs FP confident-multi (反判別)
    "confident_multi_tp_pct": pct(sum(1 for r in TP if r["wstate"] in ("S1", "S2")), len(TP)),
    "confident_multi_fp_pct": pct(sum(1 for r in FP if r["wstate"] in ("S1", "S2")), len(FP)),
    "geom_divergence": sum(r["geom_divergence"] for r in rec),
    "sumcheck_wstate": sum(S0 for S0 in cnt(rec, "wstate").values()) == N if rec else False,
}
json.dump(S, open(f"{A}/dashboard_stats_v3.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps({k: S[k] for k in ("n", "tp", "fp", "tn_total", "median_nt", "median_nn",
                                    "confident_multi_tp_pct", "confident_multi_fp_pct", "geom_divergence")}, ensure_ascii=False, indent=1))
print("funnel:", S["funnel"])
print("sumcheck:", S["sumcheck_wstate"])
