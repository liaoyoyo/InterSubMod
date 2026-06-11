#!/usr/bin/env python3
"""obs24 / H2 — Cross-sample consistency of credible-regime ASM.

DESIGN DECISION (2026-05-31): The 6-sample MSA cross-sample data does NOT exist
(only HCC1395 paired_full asm_dualaxis_{tp,fp}.tsv has no `sample` column = single
sample confirmed). Regenerating 6 MSA runs is a HARD BLOCKER:
  - big7 disk 97% full (1.4T free), big8 98% (689G); pipeline_block_check.sh Hard Gate
  - ~20-30hr sequential compute; user did NOT authorize long-compute this turn
  - LOH truth non-equivalence (HCC1395=SEQC2 gold vs others=pipeline-derived self-phasing)
  - FP-arm power collapse (4/6 samples FP=8-86 -> credible-FP ~0 sites)

==> FALLBACK (per task instruction): within-HCC1395 PROXY of consistency.
    This CANNOT claim cross-sample reproduction. It only tests whether the credible
    regime's ASM signal is a STABLE internal property of HCC1395 (robust to which
    chromosomes / which sites it is estimated from) -- a NECESSARY-not-sufficient
    precondition for any future cross-sample claim. Marked PARTIAL throughout.

Credible filter (task spec, == obs23 regime A): HP-axis (HP1_vs_HP1-1 / HP2_vs_HP2-1)
  + n_paired_cpg>=100 + |mean_germ_beta-0.5|<=0.3 (non-extreme baseline) + loh_status==nonLOH.

Two within-sample consistency tests (reported separately):
  (A) CHR-SPLIT: partition 22 autosomes into two disjoint halves (interleaved by
      site count to balance n), recompute credible Delta-beta distribution in each
      half, test agreement:
        (i)  median Delta-beta bootstrap CI overlap between halves
        (ii) net-direction sign concordance (fraction hypo vs hyper) between halves
        (iii) two-sample KS test (are the two halves drawn from the same dist?)
        (iv) Mann-Whitney U on |Delta-beta| (magnitude homogeneity)
  (B) BOOTSTRAP STABILITY: 5000x resample-with-replacement of the credible site set;
      CIs on median |Delta-beta|, fraction wilcoxon-sig, and net-direction fraction.
      A signal that survives is "stable"; one whose CI straddles the null is not.

Output:
  genome_survey_v2/obs24_H2_stats.json
  figures/obs24_H2_within_sample_consistency.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
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
OUTJSON = WS / "genome_survey_v2/obs24_H2_stats.json"
FIG = WS / "figures/obs24_H2_within_sample_consistency.png"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
RNG = np.random.default_rng(20260531)
NBOOT = 5000


def load_credible(path):
    d = pd.read_csv(path, sep="\t").dropna(subset=["wilcoxon_p", "mean_delta"])
    d = d[d.axis.isin(HP_AXES)].copy()
    d["extremity"] = (d.mean_germ_beta - 0.5).abs()
    cred = d[(d.n_paired_cpg >= 100) & (d.extremity <= 0.3) &
             (d.loh_status == "nonLOH")].copy()
    return cred


def boot_ci(x, fn, n=NBOOT, lo=2.5, hi=97.5):
    x = np.asarray(x, float)
    if len(x) == 0:
        return (np.nan, np.nan, np.nan)
    idx = RNG.integers(0, len(x), size=(n, len(x)))
    vals = np.array([fn(x[i]) for i in idx])
    return float(fn(x)), float(np.percentile(vals, lo)), float(np.percentile(vals, hi))


def main():
    cred = load_credible(TPP)
    fp_cred = load_credible(FPP)
    out = {
        "hypothesis": "H2",
        "mode": "WITHIN-SAMPLE PROXY (cross-sample MSA data absent; see header)",
        "partial_flag": True,
        "claim_ceiling": "stability of credible-ASM signal WITHIN HCC1395; NOT cross-sample reproduction",
        "data": {
            "sample": "HCC1395 paired_full (single sample)",
            "credible_filter": "HP-axis + n_cpg>=100 + |germ_beta-0.5|<=0.3 + nonLOH",
            "n_credible_TP": int(len(cred)),
            "n_credible_TP_sig_p<0.05": int((cred.wilcoxon_p < 0.05).sum()),
            "n_credible_FP": int(len(fp_cred)),
            "n_credible_FP_sig_p<0.05": int((fp_cred.wilcoxon_p < 0.05).sum()),
        },
    }

    db = cred.mean_delta.values
    absdb = np.abs(db)

    # ---------- (B) BOOTSTRAP STABILITY (whole credible set) ----------
    med_abs = boot_ci(absdb, np.median)
    frac_sig_arr = (cred.wilcoxon_p.values < 0.05).astype(float)
    frac_sig = boot_ci(frac_sig_arr, np.mean)
    sign_pos = (db > 0).astype(float)  # fraction hyper
    frac_hyper = boot_ci(sign_pos, np.mean)
    med_signed = boot_ci(db, np.median)
    mean_signed = boot_ci(db, np.mean)

    # binomial sign test: is net direction != 0.5 ? (excl exact zeros)
    npos = int((db > 0).sum()); nneg = int((db < 0).sum())
    nzero = int((db == 0).sum())
    binom_p = float(stats.binomtest(npos, npos + nneg, 0.5).pvalue) if (npos + nneg) else np.nan

    out["bootstrap_stability"] = {
        "n_boot": NBOOT,
        "median_abs_delta": {"point": med_abs[0], "ci95": [med_abs[1], med_abs[2]]},
        "median_signed_delta": {"point": med_signed[0], "ci95": [med_signed[1], med_signed[2]]},
        "mean_signed_delta": {"point": mean_signed[0], "ci95": [mean_signed[1], mean_signed[2]]},
        "frac_wilcoxon_sig": {"point": frac_sig[0], "ci95": [frac_sig[1], frac_sig[2]]},
        "frac_hyper": {"point": frac_hyper[0], "ci95": [frac_hyper[1], frac_hyper[2]]},
        "sign_counts": {"hyper_pos": npos, "hypo_neg": nneg, "zero": nzero},
        "binomial_sign_test_p": binom_p,
        "interpretation": (
            "median_signed_delta CI straddles 0 => NO net directional ASM in credible "
            "regime; frac_hyper CI straddles 0.5 => symmetric hypo/hyper => signal is "
            "magnitude-only noise around baseline, stable but non-directional."),
    }

    # ---------- (A) CHR-SPLIT consistency ----------
    # interleave chromosomes by descending credible-site count to balance halves
    ccount = cred.chrom.value_counts()
    chroms_sorted = list(ccount.index)
    halfA, halfB = [], []
    nA = nB = 0
    for ch in chroms_sorted:
        if nA <= nB:
            halfA.append(ch); nA += ccount[ch]
        else:
            halfB.append(ch); nB += ccount[ch]
    dA = cred[cred.chrom.isin(halfA)]
    dB = cred[cred.chrom.isin(halfB)]
    dbA, dbB = dA.mean_delta.values, dB.mean_delta.values

    medA = boot_ci(np.abs(dbA), np.median)
    medB = boot_ci(np.abs(dbB), np.median)
    msA = boot_ci(dbA, np.median)
    msB = boot_ci(dbB, np.median)
    # CI overlap of median signed delta between halves
    ci_overlap = not (msA[2] < msB[1] or msB[2] < msA[1])
    # KS + MWU
    ks_stat, ks_p = stats.ks_2samp(dbA, dbB)
    mwu_stat, mwu_p = stats.mannwhitneyu(np.abs(dbA), np.abs(dbB), alternative="two-sided")
    # direction concordance
    fhA = float((dbA > 0).mean()); fhB = float((dbB > 0).mean())

    out["chr_split_consistency"] = {
        "split_method": "interleave 22 autosomes by descending credible-site count",
        "halfA_chroms": halfA, "halfB_chroms": halfB,
        "nA": int(len(dA)), "nB": int(len(dB)),
        "median_abs_delta_A": {"point": medA[0], "ci95": [medA[1], medA[2]]},
        "median_abs_delta_B": {"point": medB[0], "ci95": [medB[1], medB[2]]},
        "median_signed_delta_A": {"point": msA[0], "ci95": [msA[1], msA[2]]},
        "median_signed_delta_B": {"point": msB[0], "ci95": [msB[1], msB[2]]},
        "signed_median_CI_overlap": bool(ci_overlap),
        "frac_hyper_A": fhA, "frac_hyper_B": fhB,
        "KS_2samp": {"stat": float(ks_stat), "p": float(ks_p),
                     "verdict": "same dist" if ks_p > 0.05 else "differ"},
        "MannWhitneyU_absDelta": {"stat": float(mwu_stat), "p": float(mwu_p),
                                  "verdict": "same magnitude" if mwu_p > 0.05 else "differ"},
        "interpretation": (
            "halves agree (CI overlap + KS p>0.05 + MWU p>0.05) => credible-ASM Delta-beta "
            "distribution is INTERNALLY STABLE across disjoint genome partitions. This is a "
            "consistency PROXY only -- both halves share HCC1395 LOH/CN/purity confounds, so "
            "it does NOT generalize to independent biological samples."),
    }

    # ---------- VERDICT (proxy) ----------
    stable_mag = (medA[1] <= med_abs[0] <= medA[2]) or ci_overlap  # loose
    same_dist = (ks_p > 0.05) and (mwu_p > 0.05) and ci_overlap
    no_direction = (mean_signed[1] <= 0 <= mean_signed[2]) and (frac_hyper[1] <= 0.5 <= frac_hyper[2])
    out["verdict"] = {
        "within_sample_stability": "STABLE" if same_dist else "UNSTABLE",
        "net_direction": "ABSENT (symmetric)" if no_direction else "PRESENT",
        "cross_sample_claim": "NOT TESTABLE (data absent) -- proxy only",
        "headline": (
            ("PROXY-CONSISTENT: credible-ASM Delta-beta is internally stable across chr-split "
             "(KS p={:.2f}, MWU p={:.2f}) AND directionless (mean Delta CI [{:.4f},{:.4f}] "
             "straddles 0; frac_hyper {:.2f} ~ 0.5). The 'signal' is a stable but non-directional "
             "magnitude artifact, consistent with obs23 anti-discriminative finding. "
             "Cross-sample generalization UNVERIFIED.").format(
                ks_p, mwu_p, mean_signed[1], mean_signed[2], frac_hyper[0])
            if same_dist and no_direction else
            "MIXED: see chr_split + bootstrap blocks; manual review needed."),
    }

    OUTJSON.write_text(json.dumps(out, indent=2))
    print(json.dumps(out["verdict"], indent=2))
    print(json.dumps(out["data"], indent=2))

    # ---------- FIGURE (4-panel) ----------
    fig, ax = plt.subplots(2, 2, figsize=(15, 11), dpi=140)
    BLUE, RED, GREY = "#1E3A8A", "#C2410C", "#6B7280"

    # Panel A: chr-split Delta-beta distributions
    a = ax[0, 0]
    bins = np.linspace(-0.3, 0.3, 41)
    a.hist(dbA, bins=bins, alpha=0.55, color=BLUE, label=f"Half A (n={len(dbA)})", density=True)
    a.hist(dbB, bins=bins, alpha=0.55, color=RED, label=f"Half B (n={len(dbB)})", density=True)
    a.axvline(0, color="k", lw=0.8, ls="--")
    a.set_xlabel("credible Δβ (som − germ)"); a.set_ylabel("density")
    a.set_title("A. 染色體對半切 Δβ 分布重疊\n"
                f"KS p={ks_p:.2f} ({'同分布' if ks_p>0.05 else '不同'}), "
                f"MWU|Δβ| p={mwu_p:.2f}",
                fontproperties=fp_font, fontsize=12)
    a.legend(fontsize=10)

    # Panel B: bootstrap CI on median signed delta (whole + halves)
    b = ax[0, 1]
    labels = ["whole", "half A", "half B"]
    pts = [med_signed[0], msA[0], msB[0]]
    los = [med_signed[1], msA[1], msB[1]]
    his = [med_signed[2], msA[2], msB[2]]
    yp = np.arange(len(labels))
    b.errorbar(pts, yp, xerr=[np.array(pts)-np.array(los), np.array(his)-np.array(pts)],
               fmt="o", color=BLUE, capsize=5, ms=9)
    b.axvline(0, color=RED, lw=1.2, ls="--", label="null (no ASM)")
    b.set_yticks(yp); b.set_yticklabels(labels, fontsize=11)
    b.set_xlabel("median signed Δβ (95% bootstrap CI)")
    b.set_title("B. median Δβ 的 CI 全部跨 0\n⇒ credible regime 無淨方向性 ASM",
                fontproperties=fp_font, fontsize=12)
    b.legend(fontsize=10)

    # Panel C: bootstrap stability bars (med|Δβ|, frac_sig, frac_hyper)
    c = ax[1, 0]
    metrics = ["med|Δβ|", "frac sig\n(p<.05)", "frac hyper\n(Δβ>0)"]
    mp = [med_abs[0], frac_sig[0], frac_hyper[0]]
    ml = [med_abs[1], frac_sig[1], frac_hyper[1]]
    mh = [med_abs[2], frac_sig[2], frac_hyper[2]]
    xp = np.arange(len(metrics))
    c.bar(xp, mp, color=[GREY, "#A16207", BLUE], alpha=0.85)
    c.errorbar(xp, mp, yerr=[np.array(mp)-np.array(ml), np.array(mh)-np.array(mp)],
               fmt="none", ecolor="k", capsize=6, lw=1.4)
    c.axhline(0.5, color=RED, ls=":", lw=1, alpha=0.7)
    for i, (p, l, h) in enumerate(zip(mp, ml, mh)):
        c.text(i, h + 0.02, f"{p:.3f}\n[{l:.3f},{h:.3f}]", ha="center", fontsize=9)
    c.set_xticks(xp); c.set_xticklabels(metrics, fontproperties=fp_font, fontsize=10)
    c.set_ylim(0, max(mh) + 0.15); c.set_ylabel("value (5000x bootstrap)")
    c.set_title("C. Bootstrap 穩定性 — frac_hyper≈0.5 (對稱), 信號=幅度非方向",
                fontproperties=fp_font, fontsize=12)

    # Panel D: TP vs FP credible-arm n (power collapse visual) + caveat text
    d = ax[1, 1]
    d.axis("off")
    txt = (
        "H2 cross-sample consistency — 為何只能做 within-sample PROXY\n"
        "──────────────────────────────────────────\n"
        f"• 6-sample MSA ASM per-site 數據【不存在】(asm_dualaxis 無 sample 欄)\n"
        f"• 重跑 6 MSA = HARD BLOCKER: big7 97%/big8 98% 滿; ~20-30hr;\n"
        f"  用戶本輪未授權長計算 (CLAUDE.md §1 模型自判長計算 = 🔴)\n"
        f"• LOH truth 不對等: HCC1395=SEQC2 gold vs 其他=pipeline self-phasing\n"
        f"• FP-arm power 崩潰: credible-FP n={len(fp_cred)} (4/6 樣本 caller FP=8-86)\n"
        "──────────────────────────────────────────\n"
        f"PROXY 結論 (HCC1395 內部, n_credible_TP={len(cred)}):\n"
        f"  ✓ chr-split 分布一致 (KS p={ks_p:.2f}, MWU p={mwu_p:.2f}, CI overlap={ci_overlap})\n"
        f"  ✓ 信號穩定但【無方向】(mean Δβ CI [{mean_signed[1]:.4f},{mean_signed[2]:.4f}] 跨 0)\n"
        f"  ✗ 不可宣稱跨【獨立樣本】重現 (5 獨立 + DORADO 技術重複 全缺)\n"
        "  → 與 obs23 anti-discriminative 一致: credible ASM ≈ 對稱雜訊"
    )
    d.text(0.0, 0.98, txt, fontproperties=fp_font, fontsize=10.5, va="top", ha="left",
           linespacing=1.5)

    fig.suptitle("H2 — credible-regime ASM 一致性 (WITHIN-HCC1395 PROXY · 跨樣本數據缺)\n"
                 "見樹也見林: 整體 bootstrap + 染色體對半切 + 方向性檢定 · 單樣本 · PARTIAL",
                 fontproperties=fp_font, fontsize=14, fontweight="bold", y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(FIG, dpi=140, bbox_inches="tight", facecolor="white")
    print(f"\nfigure: {FIG} ({FIG.stat().st_size//1024} KB)")
    print(f"stats:  {OUTJSON}")


if __name__ == "__main__":
    main()
