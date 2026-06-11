#!/usr/bin/env python3
"""obs23 — ASM 是不是「混合多種狀況」? regime 分解 + 整體/單獨/差分 + 觀察圖.

假設: 我們量到的甲基差異 (ASM Δβ) 不是單一現象, 而是多 regime 混合:
  (A) 可信 ASM     = 高 coverage + 中段 baseline + nonLOH  → TP-enriched, 真生物
  (B) regression   = 低 coverage + 極端 baseline           → FP-enriched, 雜訊膨脹
  (C) LOH-confound = LOH 區單 haplotype → HP1-vs-HP1-1 混 CNV/allelic-loss
驗證法: 對全 HP-axis 顯著位點 (p<0.05) 分 regime, 看每 regime 的 TP-rate / Δβ 分布是否
迥異 (迥異 = 證實混合). 同時畫 4-panel 觀察圖.

Output: figures/obs23_asm_mixture_4panel.png + genome_survey_v2/obs23_regime_summary.tsv
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

CJK = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
fp_font = FontProperties(fname=CJK) if Path(CJK).exists() else None
plt.rcParams["axes.unicode_minus"] = False

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
TPP = WS / "genome_survey_v2/asm_dualaxis_tp.tsv"
FPP = WS / "genome_survey_v2/asm_dualaxis_fp.tsv"
FIG = WS / "figures/obs23_asm_mixture_4panel.png"
OUTTSV = WS / "genome_survey_v2/obs23_regime_summary.tsv"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}


def load():
    frames = []
    for path, lab in [(TPP, 1), (FPP, 0)]:
        d = pd.read_csv(path, sep="\t").dropna(subset=["wilcoxon_p", "mean_delta"])
        d = d[d.axis.isin(HP_AXES)].copy()
        d["is_tp"] = lab
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    d["abs_delta"] = d.mean_delta.abs()
    d["extremity"] = (d.mean_germ_beta - 0.5).abs()
    d["sig"] = d.wilcoxon_p < 0.05
    return d


def regime(row):
    hi_cov = row.n_paired_cpg >= 100
    extreme = row.extremity > 0.3
    loh = row.loh_status == "LOH"
    if hi_cov and not extreme and not loh:
        return "A_credible (hiCov+midBase+nonLOH)"
    if (not hi_cov) and extreme:
        return "B_regression (loCov+extremeBase)"
    if loh:
        return "C_LOH_confound"
    return "D_intermediate"


def main():
    d = load()
    sig = d[d.sig].copy()
    sig["regime"] = sig.apply(regime, axis=1)

    # ---- regime summary (整體 vs 單獨) ----
    rows = []
    rows.append(("OVERALL_sig", len(sig), sig.is_tp.mean(), sig.abs_delta.median(),
                 sig.n_paired_cpg.median(), sig.extremity.median(), (sig.loh_status == "LOH").mean()))
    for rg in ["A_credible (hiCov+midBase+nonLOH)", "D_intermediate", "C_LOH_confound", "B_regression (loCov+extremeBase)"]:
        s = sig[sig.regime == rg]
        if len(s):
            rows.append((rg, len(s), s.is_tp.mean(), s.abs_delta.median(),
                         s.n_paired_cpg.median(), s.extremity.median(), (s.loh_status == "LOH").mean()))
    summ = pd.DataFrame(rows, columns=["regime", "n", "tp_rate", "median_abs_dbeta",
                                       "median_n_cpg", "median_extremity", "loh_frac"])
    summ.to_csv(OUTTSV, sep="\t", index=False)

    # ---- 4-panel figure ----
    fig, ax = plt.subplots(2, 2, figsize=(15, 12), dpi=140)
    tp = sig[sig.is_tp == 1]; fp = sig[sig.is_tp == 0]
    BLUE, RED = "#1E3A8A", "#C2410C"

    # Panel A: Δβ vs n_cpg (funnel: low cov -> wide spread)
    a = ax[0, 0]
    a.scatter(tp.n_paired_cpg, tp.mean_delta, s=10, c=BLUE, alpha=0.35, label="TP")
    a.scatter(fp.n_paired_cpg, fp.mean_delta, s=14, c=RED, alpha=0.5, label="FP")
    a.axhline(0, color="grey", lw=0.6); a.axvline(100, color="green", ls="--", lw=1, alpha=0.6)
    a.set_xscale("log"); a.set_xlim(5, 600); a.set_ylim(-0.8, 0.8)
    a.set_xlabel("n_cpg (coverage, log)"); a.set_ylabel("Δβ (som − germ)")
    a.set_title("A. Δβ vs coverage — 低 coverage 區 Δβ 散開 (漏斗=雜訊膨脹)", fontproperties=fp_font, fontsize=12)
    a.legend(fontsize=10)

    # Panel B: Δβ vs germline baseline (regression-to-extreme)
    b = ax[0, 1]
    b.scatter(tp.mean_germ_beta, tp.mean_delta, s=10, c=BLUE, alpha=0.35, label="TP")
    b.scatter(fp.mean_germ_beta, fp.mean_delta, s=14, c=RED, alpha=0.5, label="FP")
    b.plot([0, 1], [0, -0.5], color="grey", ls=":", lw=1)  # if germ high, hypo bias
    b.plot([0, 1], [0.5, 0], color="grey", ls=":", lw=1)
    b.axvspan(0, 0.2, alpha=0.08, color=RED); b.axvspan(0.8, 1.0, alpha=0.08, color=RED)
    b.axhline(0, color="grey", lw=0.6); b.set_xlim(0, 1); b.set_ylim(-0.8, 0.8)
    b.set_xlabel("germline mean β (baseline)"); b.set_ylabel("Δβ")
    b.set_title("B. Δβ vs baseline — 極端 baseline (紅區) 強制 regression", fontproperties=fp_font, fontsize=12)
    b.legend(fontsize=10)

    # Panel C: TP-rate heatmap by (coverage bin × baseline-extremity bin)
    c = ax[1, 0]
    sig["cov_bin"] = pd.cut(sig.n_paired_cpg, [0, 50, 100, 200, 1e9], labels=["<50", "50-100", "100-200", ">200"])
    sig["ext_bin"] = pd.cut(sig.extremity, [0, 0.15, 0.3, 0.5], labels=["mid(<.15)", "(.15-.3)", "extreme(>.3)"])
    piv = sig.pivot_table(index="ext_bin", columns="cov_bin", values="is_tp", aggfunc="mean", observed=False)
    pivn = sig.pivot_table(index="ext_bin", columns="cov_bin", values="is_tp", aggfunc="size", observed=False)
    im = c.imshow(piv.values, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    c.set_xticks(range(len(piv.columns))); c.set_xticklabels(piv.columns, fontsize=10)
    c.set_yticks(range(len(piv.index))); c.set_yticklabels(piv.index, fontproperties=fp_font, fontsize=10)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]; n = pivn.values[i, j]
            if not np.isnan(v):
                c.text(j, i, f"{v:.2f}\nn={int(n)}", ha="center", va="center", fontsize=9,
                       color="black" if 0.3 < v < 0.85 else "white")
    c.set_xlabel("coverage (n_cpg)"); c.set_ylabel("baseline extremity")
    c.set_title("C. TP-rate by coverage × baseline — 左下綠(可信) → 右上紅(雜訊)", fontproperties=fp_font, fontsize=12)
    plt.colorbar(im, ax=c, fraction=0.046, label="TP rate")

    # Panel D: regime composition (n + TP-rate)
    dd = ax[1, 1]
    rs = summ[summ.regime != "OVERALL_sig"].copy()
    short = [r.split(" ")[0] for r in rs.regime]
    yp = np.arange(len(rs))
    cols = ["#15803D" if t > 0.6 else ("#A16207" if t > 0.4 else "#C2410C") for t in rs.tp_rate]
    dd.barh(yp, rs.n, color=cols)
    dd.set_yticks(yp); dd.set_yticklabels(short, fontsize=10)
    for i, (n, t) in enumerate(zip(rs.n, rs.tp_rate)):
        dd.text(n + max(rs.n)*0.01, i, f"n={n}, TP={t:.2f}", va="center", fontsize=10)
    dd.set_xlabel("n sites (顯著 ASM)"); dd.invert_yaxis()
    dd.set_title("D. regime 組成 — TP-rate 迥異 ⇒ ASM 是混合非單一", fontproperties=fp_font, fontsize=12)

    fig.suptitle("ASM 甲基差異 = 多 regime 混合 (見樹也見林: 整體 + 單獨 + 差分)\n"
                 "HCC1395 paired_full · HP-axis · 顯著 ASM (p<0.05) · 單樣本",
                 fontproperties=fp_font, fontsize=15, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(FIG, dpi=140, bbox_inches="tight", facecolor="white")
    print(f"figure: {FIG} ({FIG.stat().st_size//1024} KB)")

    # ---- print 整體/單獨/差分 ----
    print("\n=== regime summary (整體 vs 單獨) ===")
    print(f"{'regime':<38}{'n':>6}{'TP-rate':>9}{'med|Δβ|':>9}{'med n_cpg':>10}{'LOH frac':>9}")
    for _, r in summ.iterrows():
        print(f"{r.regime:<38}{int(r.n):>6}{r.tp_rate:>9.3f}{r.median_abs_dbeta:>9.3f}{r.median_n_cpg:>10.0f}{r.loh_frac:>9.2f}")
    # differential: credible vs regression
    A = summ[summ.regime.str.startswith("A_")]
    B = summ[summ.regime.str.startswith("B_")]
    if len(A) and len(B):
        print(f"\n差分 A_credible vs B_regression: TP-rate {float(A.tp_rate.iloc[0]):.2f} vs {float(B.tp_rate.iloc[0]):.2f} "
              f"(Δ={float(A.tp_rate.iloc[0])-float(B.tp_rate.iloc[0]):+.2f})")
    print(f"written: {OUTTSV}")


if __name__ == "__main__":
    main()
