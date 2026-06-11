#!/usr/bin/env python3
"""obs22 — 統計確認: FP-strong vs TP-strong 的 4 特殊現象是否真顯著 (非 median 巧合)?

之前 obs21 只給 median。本腳本對每維跑正式檢定:
  連續 (n_cpg / |Δβ| / baseline_extremity / germ_beta): Mann-Whitney U (two-sided)
    + rank-biserial effect size r = 1 - 2U/(n1*n2)
  類別 (LOH yes/no, hypo yes/no): Fisher exact + odds ratio
小 n (FP-strong n=26) → 報精確檢定 + 標 caveat (LOH/coverage 互相關非獨立)。
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
TP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
FP = WS / "genome_survey_v2/asm_dualaxis_fp.tsv"
OUT = WS / "genome_survey_v2/obs22_tp_vs_fp_stats.json"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}


def prep_strong(path):
    d = pd.read_csv(path, sep="\t").dropna(subset=["wilcoxon_p", "mean_delta"])
    d = d[d.axis.isin(HP_AXES)].copy()
    bonf = 0.05 / len(d)
    s = d[(d.wilcoxon_p < bonf) & (d.mean_delta.abs() >= 0.1)].copy()
    s["abs_delta"] = s.mean_delta.abs()
    s["extremity"] = (s.mean_germ_beta - 0.5).abs()
    s["is_loh"] = (s.loh_status == "LOH").astype(int)
    s["is_hypo"] = (s.mean_delta < 0).astype(int)
    return s


def mwu(tp_vals, fp_vals):
    u, p = stats.mannwhitneyu(tp_vals, fp_vals, alternative="two-sided")
    n1, n2 = len(tp_vals), len(fp_vals)
    r = 1 - 2 * u / (n1 * n2)  # rank-biserial (sign: +ve => TP>FP)
    return {"U": float(u), "p": float(p), "rank_biserial_r": round(float(r), 3),
            "tp_median": round(float(np.median(tp_vals)), 3), "fp_median": round(float(np.median(fp_vals)), 3)}


def fisher(tp_bin, fp_bin):
    a = int(tp_bin.sum()); b = int(len(tp_bin) - a)
    c = int(fp_bin.sum()); d = int(len(fp_bin) - c)
    odd, p = stats.fisher_exact([[a, b], [c, d]])
    return {"tp_yes": a, "tp_no": b, "fp_yes": c, "fp_no": d,
            "odds_ratio": round(float(odd), 3), "p": float(p),
            "tp_pct": round(100*a/(a+b), 1), "fp_pct": round(100*c/(c+d), 1)}


def main():
    tp, fp = prep_strong(TP), prep_strong(FP)
    res = {
        "analysis": "obs22_TP_strong_vs_FP_strong_significance",
        "date": "2026-05-31", "axis": "HP-axis (somatic-controlled)",
        "n": {"TP_strong": int(len(tp)), "FP_strong": int(len(fp))},
        "continuous_MannWhitneyU": {
            "n_cpg (coverage)": mwu(tp.n_paired_cpg, fp.n_paired_cpg),
            "baseline_extremity": mwu(tp.extremity, fp.extremity),
            "abs_dbeta": mwu(tp.abs_delta, fp.abs_delta),
            "germ_beta": mwu(tp.mean_germ_beta, fp.mean_germ_beta),
        },
        "categorical_Fisher": {
            "LOH": fisher(tp.is_loh, fp.is_loh),
            "hypo": fisher(tp.is_hypo, fp.is_hypo),
        },
    }
    # bonferroni over the 6 tests
    pvals = [res["continuous_MannWhitneyU"][k]["p"] for k in res["continuous_MannWhitneyU"]] + \
            [res["categorical_Fisher"][k]["p"] for k in res["categorical_Fisher"]]
    bonf6 = 0.05 / 6
    res["multiple_testing"] = {"n_tests": 6, "bonferroni_alpha": bonf6,
                               "n_significant_bonf": int(sum(p < bonf6 for p in pvals))}
    # verdict per dimension
    verdict = {}
    for k, v in res["continuous_MannWhitneyU"].items():
        verdict[k] = "CONFIRMED_diff" if v["p"] < bonf6 else ("nominal" if v["p"] < 0.05 else "ns")
    for k, v in res["categorical_Fisher"].items():
        verdict[k] = "CONFIRMED_diff" if v["p"] < bonf6 else ("nominal" if v["p"] < 0.05 else "ns")
    res["verdict_per_dim"] = verdict
    res["caveat"] = ("Small FP-strong n=%d. LOH/coverage/baseline are mechanistically correlated "
                     "(not independent axes) — they co-describe one regression-to-extreme artifact, "
                     "not 4 separate findings. Single-sample HCC1395." % len(fp))
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=2))

    print(f"TP-strong n={len(tp)} vs FP-strong n={len(fp)}  (Bonferroni alpha over 6 tests = {bonf6:.4f})")
    print(f"\n{'dimension':<22}{'TP med':>9}{'FP med':>9}{'p':>11}{'effect':>9}  verdict")
    for k, v in res["continuous_MannWhitneyU"].items():
        print(f"{k:<22}{v['tp_median']:>9}{v['fp_median']:>9}{v['p']:>11.2e}{v['rank_biserial_r']:>9}  {verdict[k]}")
    for k, v in res["categorical_Fisher"].items():
        print(f"{k+' %':<22}{v['tp_pct']:>9}{v['fp_pct']:>9}{v['p']:>11.2e}{('OR='+str(v['odds_ratio'])):>9}  {verdict[k]}")
    print(f"\n{res['multiple_testing']['n_significant_bonf']}/6 survive Bonferroni")
    print(f"caveat: {res['caveat']}")
    print(f"written: {OUT}")


if __name__ == "__main__":
    main()
