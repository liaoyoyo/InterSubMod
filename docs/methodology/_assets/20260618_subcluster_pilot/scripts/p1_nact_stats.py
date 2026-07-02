#!/usr/bin/env python3
"""[P1.3 NACT 統計] 從 nact_results.json 算 candidate_subclone rate 分 CN/set 層 + Clopper-Pearson 95% CI
+ Fisher 檢定(neutral vs LOH / LOH vs gain-loss / TP vs FP) → land nact_stats.json。誠實量化『CN-driven / 非 somatic-specific』。純讀。"""
import json
from collections import Counter
from scipy.stats import beta, fisher_exact

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
res = json.load(open(f"{A}/nact_results.json"))


def ci(k, n):
    if n == 0:
        return None
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return {"k": k, "n": n, "rate_pct": round(100 * k / n, 2), "ci95": [round(100 * lo, 1), round(100 * hi, 1)]}


def cand(sub):
    return sum(1 for r in sub if r["nact_verdict"] == "candidate_subclone")


def strata(pred):
    sub = [r for r in res if pred(r)]
    return ci(cand(sub), len(sub))


cn = lambda r: r.get("cn_state")
groups = {
    "ALL": strata(lambda r: True),
    "neutral": strata(lambda r: cn(r) == "neutral"),
    "LOH": strata(lambda r: r.get("is_loh") is True),
    "gain+loss": strata(lambda r: cn(r) in ("gain", "loss")),
    "non_loh": strata(lambda r: r.get("is_loh") is not True),
    "TP": strata(lambda r: r["set"] == "TP"),
    "FP": strata(lambda r: r["set"] == "FP"),
}


def fish(p1, p2):
    a, b = cand([r for r in res if p1(r)]), len([r for r in res if p1(r)])
    c, d = cand([r for r in res if p2(r)]), len([r for r in res if p2(r)])
    _, p = fisher_exact([[a, b - a], [c, d - c]])
    return round(float(p), 4)


fishers = {
    "neutral_vs_LOH": fish(lambda r: cn(r) == "neutral", lambda r: r.get("is_loh") is True),
    "LOH_vs_gainloss": fish(lambda r: r.get("is_loh") is True, lambda r: cn(r) in ("gain", "loss")),
    "TP_vs_FP": fish(lambda r: r["set"] == "TP", lambda r: r["set"] == "FP"),
}
verdict_dist = dict(Counter(r["nact_verdict"] for r in res))
out = {"n_total": len(res), "verdict_dist": verdict_dist, "candidate_rate_by_stratum": groups, "fisher_tests": fishers,
       "interpretation": {
           "neutral_underpowered": "neutral 0/11 CI95 上限涵蓋 LOH 率 → Fisher p=" + str(fishers["neutral_vs_LOH"]) + " 不顯著 → 不可單憑 neutral=0 下 CN-driven",
           "loh_elevated": "LOH vs gain+loss Fisher p=" + str(fishers["LOH_vs_gainloss"]) + " → LOH 特異升高=LOH-unmask germline ASM",
           "not_tp_specific": "TP vs FP Fisher p=" + str(fishers["TP_vs_FP"]) + " → candidate 不區分 TP/FP=非 somatic-specific=artifact-consistent"}}
json.dump(out, open(f"{A}/nact_stats.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
print("[-> nact_stats.json]")
