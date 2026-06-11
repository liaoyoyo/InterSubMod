#!/usr/bin/env python3
"""obs20 — 見樹也見林: per-site ASM 統一輸出 (forest aggregate + ranked trees).

軸 2（甲基 × HP × somatic）的「見樹也見林」統一統計。Forest（全域漏斗/方向/TP-vs-FP）
既有 synthesis 已做透 — 本腳本重現 forest + 補上「見樹」= 排序個案表（4 類：
extreme / canonical / well-explained / FP-enriched），讓個案與大規模統計同一視圖。

口徑：HP-axis (HP1_vs_HP1-1, HP2_vs_HP2-1) = somatic-controlled 乾淨軸（主結論）；
ALLELE-axis (ALT_vs_REF) = baseline allelic 甲基化 confounded（只列不下結論，per
feedback_asm_allele_axis_baseline_confound）。strong-ASM = Bonferroni AND |Δβ|>=0.1。
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
TP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
FP = WS / "genome_survey_v2/asm_dualaxis_fp.tsv"
OUT_TREES = WS / "genome_survey_v2/obs20_ranked_trees.tsv"
OUT_FOREST = WS / "genome_survey_v2/obs20_forest_stats.json"

HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
BRCA2 = ("chr13", 32315128)


def load(path, label):
    df = pd.read_csv(path, sep="\t")
    df["set"] = label
    df = df.dropna(subset=["wilcoxon_p", "mean_delta"])
    return df


def main():
    tp = load(TP, "TP")
    fp = load(FP, "FP")
    alld = pd.concat([tp, fp], ignore_index=True)

    n_valid = len(tp)  # forest funnel is over TP universe (the discovery set)
    bonf = 0.05 / n_valid

    # strong-ASM on TP, HP-axis (clean) vs ALLELE (confounded, reported separately)
    tp["is_bonf"] = tp["wilcoxon_p"] < bonf
    tp["is_strong"] = tp["is_bonf"] & (tp["mean_delta"].abs() >= 0.1)
    hp = tp[tp["axis"].isin(HP_AXES)].copy()
    al = tp[tp["axis"] == "ALT_vs_REF"].copy()

    # ---------- FOREST (aggregate) ----------
    def funnel(d):
        return {
            "n_records": int(len(d)),
            "p<0.05": int((d["wilcoxon_p"] < 0.05).sum()),
            "bonferroni": int(d["is_bonf"].sum()),
            "strong_ASM_bonf_and_dbeta>=0.1": int(d["is_strong"].sum()),
        }

    def direction(d):
        s = d[d["is_strong"]]
        hypo = int((s["mean_delta"] < 0).sum()); hyper = int((s["mean_delta"] > 0).sum())
        return {"n_strong": len(s), "hypo": hypo, "hyper": hyper,
                "hypo_pct": round(100*hypo/len(s), 1) if len(s) else None}

    # TP-vs-FP enrichment of strong-ASM (the anti-discriminative finding)
    fp["is_bonf"] = fp["wilcoxon_p"] < (0.05 / max(len(fp), 1))
    fp["is_strong"] = fp["is_bonf"] & (fp["mean_delta"].abs() >= 0.1)
    tp_strong_rate = tp["is_strong"].mean()
    fp_strong_rate = fp["is_strong"].mean()

    forest = {
        "analysis": "obs20_see_forest_axis2_methylation_ASM",
        "date": "2026-05-30",
        "bonferroni_threshold": bonf,
        "n_valid_TP": n_valid,
        "axis_breakdown_TP": {
            "HP_axis (clean, somatic-controlled)": funnel(hp),
            "ALLELE_axis (CONFOUNDED by baseline allelic meth)": funnel(al),
            "ALL_TP": funnel(tp),
        },
        "direction_strong_ALL": direction(tp),
        "direction_strong_HP": direction(hp),
        "TP_vs_FP_strong_rate": {
            "TP_strong_rate": round(float(tp_strong_rate), 5),
            "FP_strong_rate": round(float(fp_strong_rate), 5),
            "note": "strong-ASM anti-discriminative: enriched in FP (per synthesis OR=0.194)",
        },
        "effect_size_HP_strong": {
            "median_abs_dbeta": round(float(hp[hp.is_strong]["mean_delta"].abs().median()), 4) if hp.is_strong.any() else None,
            "median_n_cpg": int(hp[hp.is_strong]["n_paired_cpg"].median()) if hp.is_strong.any() else None,
        },
    }
    OUT_FOREST.write_text(json.dumps(forest, ensure_ascii=False, indent=2))

    # ---------- TREES (ranked individual cases) ----------
    cols = ["chrom", "somatic_pos", "axis", "loh_status", "n_paired_cpg",
            "mean_germ_beta", "mean_som_beta", "mean_delta", "wilcoxon_p"]
    hp_strong = hp[hp.is_strong].copy()
    hp_strong["abs_delta"] = hp_strong["mean_delta"].abs()

    trees = []

    def tag(rows, category):
        for _, r in rows.iterrows():
            trees.append({**{c: r[c] for c in cols}, "category": category,
                          "abs_delta": round(abs(r["mean_delta"]), 4)})

    # 樹1 extreme outlier: top 8 by |Δβ| (HP, strong)
    tag(hp_strong.sort_values("abs_delta", ascending=False).head(8), "extreme_|dbeta|")
    # 樹2 top significance: top 8 by p (HP, strong)
    tag(hp_strong.sort_values("wilcoxon_p").head(8), "top_significance")
    # 樹3 canonical: 8 nearest the median |Δβ| of HP strong
    if len(hp_strong):
        med = hp_strong["abs_delta"].median()
        hp_strong["dist_med"] = (hp_strong["abs_delta"] - med).abs()
        tag(hp_strong.sort_values("dist_med").head(8), "canonical_median")
    # 樹4 well-explained: BRCA2 (any axis)
    brca = tp[(tp.chrom == BRCA2[0]) & (tp.somatic_pos == BRCA2[1])]
    tag(brca, "well_explained_BRCA2")
    # 樹5 FP-enriched cautionary: top 6 FP strong-ASM by |Δβ| (HP)
    fp_hp_strong = fp[(fp.axis.isin(HP_AXES)) & (fp.is_strong)].copy()
    fp_hp_strong["abs_delta"] = fp_hp_strong["mean_delta"].abs()
    tag(fp_hp_strong.sort_values("abs_delta", ascending=False).head(6), "FP_enriched_cautionary")

    tdf = pd.DataFrame(trees)
    tdf.to_csv(OUT_TREES, sep="\t", index=False)

    # ---------- print summary ----------
    print("=== 見林 (FOREST) — HP-axis (clean) ===")
    f = forest["axis_breakdown_TP"]["HP_axis (clean, somatic-controlled)"]
    print(f"  HP-axis funnel: {f['n_records']} -> p<.05 {f['p<0.05']} -> Bonf {f['bonferroni']} -> strong {f['strong_ASM_bonf_and_dbeta>=0.1']}")
    print(f"  direction(HP strong): {forest['direction_strong_HP']}")
    print(f"  effect(HP strong): {forest['effect_size_HP_strong']}")
    print(f"  TP-vs-FP strong rate: TP {forest['TP_vs_FP_strong_rate']['TP_strong_rate']} / FP {forest['TP_vs_FP_strong_rate']['FP_strong_rate']} (anti-discriminative)")
    print(f"  ALLELE-axis (CONFOUNDED): {forest['axis_breakdown_TP']['ALLELE_axis (CONFOUNDED by baseline allelic meth)']}")
    print()
    print("=== 見樹 (TREES) — ranked individual cases (HP-axis unless noted) ===")
    for cat in ["well_explained_BRCA2", "extreme_|dbeta|", "top_significance", "canonical_median", "FP_enriched_cautionary"]:
        sub = tdf[tdf.category == cat]
        print(f"\n[{cat}] ({len(sub)})")
        for _, r in sub.head(5).iterrows():
            print(f"  {r['chrom']}:{int(r['somatic_pos'])} {r['axis']} {r['loh_status']:>6} "
                  f"Δβ={r['mean_delta']:+.3f} (germ {r['mean_germ_beta']:.3f}->som {r['mean_som_beta']:.3f}) "
                  f"n_cpg={int(r['n_paired_cpg'])} p={r['wilcoxon_p']:.1e}")
    print(f"\nwritten: {OUT_TREES} ({len(tdf)} ranked rows) + {OUT_FOREST}")


if __name__ == "__main__":
    main()
