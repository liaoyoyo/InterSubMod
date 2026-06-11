#!/usr/bin/env python3
"""obs21 — TP vs FP ASM 差異剖析: 為什麼 strong-ASM 是 anti-discriminative?

問題: strong-ASM 在 FP 富集 ~10×。FP 的 strong-ASM 是真生物訊號還是 artifact?
比較 TP-strong vs FP-strong (HP-axis) 在: n_cpg(coverage) / germline baseline beta /
|Δβ| / 方向 / LOH。若 FP-strong = 低 n_cpg + 極端 baseline + 大 |Δβ| → 是雜訊膨脹
(regression-to-extreme) 而非真 ASM, 即可解釋 anti-discriminative。
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
TP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
FP = WS / "genome_survey_v2/asm_dualaxis_fp.tsv"
OUT = WS / "genome_survey_v2/obs21_tp_vs_fp_asm.json"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}


def prep(path):
    d = pd.read_csv(path, sep="\t").dropna(subset=["wilcoxon_p", "mean_delta"])
    d = d[d.axis.isin(HP_AXES)].copy()
    return d


def describe(d, label):
    if len(d) == 0:
        return {"label": label, "n": 0}
    gb = d["mean_germ_beta"]
    # baseline extremity = distance of germline beta from 0.5 (max=0.5 = fully extreme)
    extremity = (gb - 0.5).abs()
    return {
        "label": label, "n": int(len(d)),
        "median_n_cpg": float(d["n_paired_cpg"].median()),
        "median_abs_dbeta": float(d["mean_delta"].abs().median()),
        "median_germ_beta": float(gb.median()),
        "median_baseline_extremity": float(extremity.median()),
        "pct_germ_extreme(>0.4 from .5)": float((extremity > 0.4).mean() * 100),
        "hypo_pct": float((d["mean_delta"] < 0).mean() * 100),
        "loh_pct": float((d["loh_status"] == "LOH").mean() * 100),
    }


def main():
    tp, fp = prep(TP), prep(FP)
    tp_bonf = 0.05 / len(tp)
    fp_bonf = 0.05 / len(fp)
    tp["strong"] = (tp.wilcoxon_p < tp_bonf) & (tp.mean_delta.abs() >= 0.1)
    fp["strong"] = (fp.wilcoxon_p < fp_bonf) & (fp.mean_delta.abs() >= 0.1)

    out = {
        "analysis": "obs21_TP_vs_FP_strong_ASM_characterization",
        "date": "2026-05-31",
        "axis": "HP-axis (somatic-controlled)",
        "strong_rate": {
            "TP": round(float(tp.strong.mean()), 5),
            "FP": round(float(fp.strong.mean()), 5),
            "FP_over_TP": round(float(fp.strong.mean() / tp.strong.mean()), 2),
        },
        "TP_strong": describe(tp[tp.strong], "TP-strong"),
        "FP_strong": describe(fp[fp.strong], "FP-strong"),
        "TP_all": describe(tp, "TP-all"),
        "FP_all": describe(fp, "FP-all"),
    }
    # interpretation flags
    ts, fs = out["TP_strong"], out["FP_strong"]
    out["interpretation"] = {
        "FP_strong_lower_coverage": fs["median_n_cpg"] < ts["median_n_cpg"],
        "FP_strong_more_extreme_baseline": fs["median_baseline_extremity"] > ts["median_baseline_extremity"],
        "FP_strong_larger_dbeta": fs["median_abs_dbeta"] > ts["median_abs_dbeta"],
        "n_cpg_ratio_TP_over_FP": round(ts["median_n_cpg"] / fs["median_n_cpg"], 2) if fs["median_n_cpg"] else None,
    }
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))

    print("=== TP vs FP strong-ASM (HP-axis) ===")
    print(f"strong rate: TP {out['strong_rate']['TP']} / FP {out['strong_rate']['FP']} = {out['strong_rate']['FP_over_TP']}x FP-enriched")
    print(f"{'metric':<32}{'TP-strong':>12}{'FP-strong':>12}")
    for k in ["n", "median_n_cpg", "median_abs_dbeta", "median_germ_beta",
              "median_baseline_extremity", "pct_germ_extreme(>0.4 from .5)", "hypo_pct", "loh_pct"]:
        tv, fv = ts.get(k), fs.get(k)
        tvs = f"{tv:.3f}" if isinstance(tv, float) else str(tv)
        fvs = f"{fv:.3f}" if isinstance(fv, float) else str(fv)
        print(f"{k:<32}{tvs:>12}{fvs:>12}")
    print(f"\ninterpretation: {out['interpretation']}")
    print(f"written: {OUT}")


if __name__ == "__main__":
    main()
