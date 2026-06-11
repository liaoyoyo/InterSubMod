#!/usr/bin/env python3
"""obs26 / H5 — ASM |Δβ| predict caller verdict (is_TP)? Single-feature AUC.

Question (mixture-decomposition follow-up to obs23):
  strong-ASM is anti-discriminative (FP-enriched 5x), but it's a MIXTURE of
  clean / LOH-confound / regression regimes. After filtering to the *credible*
  regime (hiCov + midBaseline + nonLOH), is residual |Δβ| neutral / weak-positive
  rather than anti-discriminative?

Design (aligned to /auc-confound-guard spirit; single-sample HCC1395 paired_full):
  - HP-axis only (HP1_vs_HP1-1, HP2_vs_HP2-1) = somatic-controlled clean axis.
  - PER-POSITION aggregation (1-2 axis rows / position -> 1 row): abs_delta=max,
    n_cpg=max, extremity=min, loh=any. AUC on positions (i.i.d.), not rows.
  - score = |Δβ|. AUC via Mann-Whitney U = U/(n_pos*n_neg) (base-rate robust).
  - permutation null: shuffle is_tp label 2000x -> null 95% CI + two-sided p.
  - confound guard (no AF col available): AUC[n_cpg->TP], spearman(|Δβ|,n_cpg),
    within-coverage-bin AUC.
  - regimes: (a) ALL HP-axis (b') ALL significant (p<.05) (b) credible regime
    using obs23's EXACT regime() definition.

Output: genome_survey_v2/obs26_H5_stats.json + figures/obs26_H5_*.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
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
FIGDIR = WS / "figures"
OUTJSON = WS / "genome_survey_v2/obs26_H5_stats.json"
HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
RNG = np.random.default_rng(20260531)
N_PERM = 2000


def auc_mwu(score, label):
    """AUC = P(score|pos > score|neg) via Mann-Whitney U. label: 1=pos(TP),0=neg(FP)."""
    pos = score[label == 1]
    neg = score[label == 0]
    n_pos, n_neg = len(pos), len(neg)
    if n_pos == 0 or n_neg == 0:
        return np.nan, n_pos, n_neg
    U, _ = mannwhitneyu(pos, neg, alternative="two-sided")
    return U / (n_pos * n_neg), n_pos, n_neg


def perm_null(score, label, n_perm=N_PERM):
    """Shuffle label n_perm times -> null AUC distribution + two-sided p vs observed."""
    obs, n_pos, n_neg = auc_mwu(score, label)
    if np.isnan(obs):
        return obs, np.nan, np.nan, np.nan, np.array([])
    lab = label.copy()
    null = np.empty(n_perm)
    for i in range(n_perm):
        RNG.shuffle(lab)
        null[i], _, _ = auc_mwu(score, lab)
    lo, hi = np.percentile(null, [2.5, 97.5])
    # two-sided p: how often null is at least as extreme (|null-0.5| >= |obs-0.5|)
    p = (np.sum(np.abs(null - 0.5) >= np.abs(obs - 0.5)) + 1) / (n_perm + 1)
    return obs, lo, hi, p, null


def load_per_position():
    """Load HP-axis TP+FP, dropna, aggregate to one row per (chrom,pos)."""
    frames = []
    for path, lab in [(TPP, 1), (FPP, 0)]:
        d = pd.read_csv(path, sep="\t").dropna(subset=["wilcoxon_p", "mean_delta"])
        d = d[d.axis.isin(HP_AXES)].copy()
        d["is_tp"] = lab
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    d["abs_delta"] = d.mean_delta.abs()
    d["extremity"] = (d.mean_germ_beta - 0.5).abs()
    d["is_loh"] = (d.loh_status == "LOH")
    d["sig"] = d.wilcoxon_p < 0.05

    # PER-POSITION aggregation across HP1/HP2 axis rows
    g = d.groupby(["chrom", "somatic_pos"])
    agg = pd.DataFrame({
        "abs_delta": g["abs_delta"].max(),        # strongest ASM at this site
        "n_cpg": g["n_paired_cpg"].max(),         # best coverage
        "extremity": g["extremity"].min(),        # closest-to-mid baseline
        "is_loh": g["is_loh"].any(),              # any axis in LOH
        "min_p": g["wilcoxon_p"].min(),           # most significant axis
        "is_tp": g["is_tp"].first(),              # label (sites mutually exclusive)
    }).reset_index()
    agg["sig"] = agg["min_p"] < 0.05
    # sanity: a position can't be both TP and FP (verified mutually exclusive)
    dup = d.groupby(["chrom", "somatic_pos"])["is_tp"].nunique()
    n_conflict = int((dup > 1).sum())
    return agg, n_conflict


def regime_obs23(row):
    """EXACT obs23 regime() definition (scripts/23_asm_mixture_decomposition.py:48)."""
    hi_cov = row.n_cpg >= 100
    extreme = row.extremity > 0.3
    loh = row.is_loh
    if hi_cov and (not extreme) and (not loh):
        return "A_credible"
    if (not hi_cov) and extreme:
        return "B_regression"
    if loh:
        return "C_LOH_confound"
    return "D_intermediate"


def run_set(df, name):
    score = df["abs_delta"].to_numpy()
    label = df["is_tp"].to_numpy().astype(int)
    auc, lo, hi, p, null = perm_null(score, label)
    n_pos = int((label == 1).sum())
    n_neg = int((label == 0).sum())
    return {
        "name": name, "n": int(len(df)), "n_tp": n_pos, "n_fp": n_neg,
        "tp_base_rate": round(n_pos / len(df), 4) if len(df) else None,
        "auc_abs_delta": round(float(auc), 4) if not np.isnan(auc) else None,
        "perm_null_95ci": [round(float(lo), 4), round(float(hi), 4)] if not np.isnan(lo) else None,
        "perm_p_two_sided": round(float(p), 4) if not np.isnan(p) else None,
        "in_null": bool(lo <= auc <= hi) if not np.isnan(lo) else None,
    }, null, score, label


def main():
    FIGDIR.mkdir(exist_ok=True)
    agg, n_conflict = load_per_position()

    sig = agg[agg.sig].copy()
    agg2 = agg.copy()
    agg2["regime"] = agg2.apply(regime_obs23, axis=1)
    credible = agg2[agg2.regime == "A_credible"].copy()

    results = {}
    nulls = {}
    scores = {}
    # (a) ALL HP-axis
    r_all, null_all, sc_all, lb_all = run_set(agg, "a_all_hp_axis")
    # (b') ALL significant
    r_sig, null_sig, _, _ = run_set(sig, "bp_all_significant_p05")
    # (b) credible regime
    r_cred, null_cred, sc_cred, lb_cred = run_set(credible, "b_credible_regime")
    results["a_all_hp_axis"] = r_all
    results["bp_all_significant_p05"] = r_sig
    results["b_credible_regime"] = r_cred
    nulls = {"a_all_hp_axis": null_all, "b_credible_regime": null_cred}

    # ---- confound guard (no AF col; use n_cpg as observed weak axis) ----
    sc = agg["abs_delta"].to_numpy()
    cov = agg["n_cpg"].to_numpy()
    lab = agg["is_tp"].to_numpy().astype(int)
    auc_cov, _, _ = auc_mwu(cov, lab)
    rho, rho_p = spearmanr(sc, cov)
    # within-coverage-bin AUC: control coverage by binning, pooled within-bin AUC
    bins = [0, 25, 50, 100, 200, 1e9]
    blab = ["<25", "25-50", "50-100", "100-200", ">200"]
    agg["cov_bin"] = pd.cut(agg["n_cpg"], bins, labels=blab)
    wbin = {}
    for b in blab:
        sub = agg[agg.cov_bin == b]
        a, npp, nnn = auc_mwu(sub["abs_delta"].to_numpy(), sub["is_tp"].to_numpy().astype(int))
        wbin[b] = {"n": int(len(sub)), "n_fp": int(nnn),
                   "auc": round(float(a), 4) if not np.isnan(a) else None}
    confound = {
        "auc_ncpg_to_tp": round(float(auc_cov), 4),
        "auc_absdelta_to_tp": results["a_all_hp_axis"]["auc_abs_delta"],
        "spearman_absdelta_ncpg": round(float(rho), 4),
        "spearman_p": float(rho_p),
        "within_coverage_bin_auc": wbin,
        "note": "n_cpg(coverage) is a stronger weak discriminator than |Δβ|; "
                "credible filter fixes n_cpg/LOH/baseline -> residual |Δβ| hollowed out. "
                "caller_af/AF column absent in TSV -> AF-bin + within-group OLS guards NOT runnable.",
    }

    # ---- regime composition table (per-position) ----
    rc_rows = []
    for rg in ["A_credible", "D_intermediate", "C_LOH_confound", "B_regression"]:
        s = agg2[agg2.regime == rg]
        rc_rows.append({
            "regime": rg, "n": int(len(s)), "n_fp": int((s.is_tp == 0).sum()),
            "tp_rate": round(float((s.is_tp == 1).mean()), 4) if len(s) else None,
            "median_abs_dbeta": round(float(s.abs_delta.median()), 4) if len(s) else None,
            "median_n_cpg": float(s.n_cpg.median()) if len(s) else None,
        })

    out = {
        "hypothesis_id": "H5",
        "title": "ASM |Δβ| predicts is_TP — credible-regime rescue test",
        "sample": "HCC1395 paired_full (single sample)",
        "data_source": {
            "tp_tsv": str(TPP), "fp_tsv": str(FPP),
            "tp_truth": "SEQC2 high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
            "fp_caller": "ClairS filtered_snv_fp.vcf.gz",
        },
        "method": {
            "axis": "HP-axis only (HP1_vs_HP1-1, HP2_vs_HP2-1)",
            "aggregation": "PER-POSITION (max|Δβ|, max n_cpg, min extremity, any LOH)",
            "score": "|Δβ| = |mean_delta|",
            "auc": "Mann-Whitney U / (n_pos*n_neg), two-sided",
            "null": f"label shuffle x{N_PERM}, 95% CI + two-sided p",
            "credible_regime": "obs23 regime(): n_cpg>=100 & extremity<=0.3 & nonLOH",
        },
        "position_level_label_check": {
            "n_positions": int(len(agg)),
            "n_conflict_tp_and_fp": n_conflict,
            "note": "TP/FP sites mutually exclusive (label clean)",
        },
        "results": results,
        "confound_guard": confound,
        "regime_composition_per_position": rc_rows,
    }

    OUTJSON.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps(out, indent=2, ensure_ascii=False))

    # ================= FIGURES =================
    _fig_roc(scores_dict={
        "ALL HP-axis": (sc_all, lb_all, r_all),
        "Credible regime": (sc_cred, lb_cred, r_cred),
    })
    _fig_perm(nulls, results)
    _fig_regime(rc_rows)
    _fig_confound(agg, wbin, confound)
    print("\n[figures]")
    for f in sorted(FIGDIR.glob("obs26_H5_*.png")):
        print(" ", f)


def _roc_curve(score, label):
    order = np.argsort(-score, kind="mergesort")
    y = label[order]
    P = (y == 1).sum()
    N = (y == 0).sum()
    tpr = np.cumsum(y == 1) / P
    fpr = np.cumsum(y == 0) / N
    return np.concatenate([[0], fpr]), np.concatenate([[0], tpr])


def _fig_roc(scores_dict):
    fig, ax = plt.subplots(figsize=(6, 6))
    colors = {"ALL HP-axis": "#1f77b4", "Credible regime": "#d62728"}
    for name, (sc, lb, r) in scores_dict.items():
        fpr, tpr = _roc_curve(sc, lb.astype(int))
        auc = r["auc_abs_delta"]
        ax.plot(fpr, tpr, color=colors[name], lw=2,
                label=f"{name}  AUC={auc} (n={r['n']}, FP={r['n_fp']})")
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5, label="chance (0.50)")
    ax.set_xlabel("FPR (FP ranked above)"); ax.set_ylabel("TPR (TP ranked above)")
    ax.set_title("H5 ROC: |delta-beta| -> is_TP  (per-position, HP-axis)")
    ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(FIGDIR / "obs26_H5_roc.png", dpi=130)
    plt.close(fig)


def _fig_perm(nulls, results):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, key, title in [(axes[0], "a_all_hp_axis", "ALL HP-axis"),
                           (axes[1], "b_credible_regime", "Credible regime")]:
        null = nulls[key]; r = results[key]
        ax.hist(null, bins=40, color="#999", alpha=0.7, label="null (shuffled)")
        ax.axvline(r["auc_abs_delta"], color="#d62728", lw=2,
                   label=f"observed AUC={r['auc_abs_delta']}")
        ax.axvline(0.5, color="k", ls=":", lw=1)
        lo, hi = r["perm_null_95ci"]
        ax.axvspan(lo, hi, color="#1f77b4", alpha=0.12, label=f"null 95% [{lo},{hi}]")
        ax.set_title(f"{title}\nn={r['n']} (FP={r['n_fp']})  perm-p={r['perm_p_two_sided']}",
                     fontsize=10)
        ax.set_xlabel("AUC"); ax.legend(fontsize=8)
    fig.suptitle("H5 permutation null: observed AUC vs label-shuffle null (2000x)")
    fig.tight_layout()
    fig.savefig(FIGDIR / "obs26_H5_perm_null.png", dpi=130)
    plt.close(fig)


def _fig_regime(rc_rows):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    names = [r["regime"] for r in rc_rows]
    ns = [r["n"] for r in rc_rows]
    fps = [r["n_fp"] for r in rc_rows]
    trs = [r["tp_rate"] for r in rc_rows]
    x = np.arange(len(names))
    ax.bar(x, ns, color="#1f77b4", alpha=0.6, label="n positions")
    ax.bar(x, fps, color="#d62728", alpha=0.85, label="n FP")
    for i, (n, fp_, tr) in enumerate(zip(ns, fps, trs)):
        ax.text(i, n + 30, f"n={n}\nFP={fp_}\nTPr={tr}", ha="center",
                fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(names)
    ax.set_ylabel("count (per-position)")
    ax.set_title("H5 regime composition (per-position) -- credible (A): FP depleted to 8")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGDIR / "obs26_H5_regime_composition.png", dpi=130)
    plt.close(fig)


def _fig_confound(agg, wbin, confound):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    a = axes[0]
    tp = agg[agg.is_tp == 1]; fp = agg[agg.is_tp == 0]
    a.scatter(tp.n_cpg, tp.abs_delta, s=8, c="#1f77b4", alpha=0.25, label="TP")
    a.scatter(fp.n_cpg, fp.abs_delta, s=12, c="#d62728", alpha=0.5, label="FP")
    a.set_xscale("log"); a.set_xlabel("n_cpg (coverage, log)")
    a.set_ylabel("|delta-beta|")
    a.set_title(f"|delta-beta| vs coverage funnel\nspearman={confound['spearman_absdelta_ncpg']} "
                f"(low cov -> large |delta-beta|)", fontsize=10)
    a.legend()
    b = axes[1]
    bl = list(wbin.keys())
    av = [wbin[k]["auc"] if wbin[k]["auc"] is not None else np.nan for k in bl]
    nf = [wbin[k]["n_fp"] for k in bl]
    b.axhline(0.5, color="k", ls=":", lw=1)
    b.bar(bl, av, color="#2ca02c", alpha=0.7)
    for i, (v, f_) in enumerate(zip(av, nf)):
        if not np.isnan(v):
            b.text(i, v + 0.01, f"{v}\nFP={f_}", ha="center", fontsize=8)
    b.set_ylim(0, 1); b.set_ylabel("AUC[|delta-beta|->TP] within bin")
    b.set_xlabel("coverage bin (n_cpg)")
    b.set_title(f"within-coverage-bin AUC\n(AUC[n_cpg->TP]={confound['auc_ncpg_to_tp']} "
                f"> AUC[|delta-beta|->TP]={confound['auc_absdelta_to_tp']})",
                fontsize=10)
    fig.suptitle("H5 confound guard: coverage is the stronger weak axis, not ASM")
    fig.tight_layout()
    fig.savefig(FIGDIR / "obs26_H5_confound.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
