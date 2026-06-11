#!/usr/bin/env python3
"""
43_q2_dose_response.py — Q2 CN dose-response (pre-reg H-A).

H-A: HP-axis ASM |Δβ| is INDEPENDENT of ground-truth median CN.
Falsifier: coverage-controlled Spearman rho(|Δβ|, median_cn) > 0.2 AND same-direction p<0.05.
Threshold: partial-rho <= 0.1 OR opposite sign -> NOT CN-driven (H-A SUPPORTED).
Verdict driver = PARTIAL (coverage-controlled) Spearman on HP-axis among sig sites.

Among sig sites (wilcoxon_p<0.05), separately for axis_type HP and ALLELE:
  (a) raw Spearman |Δβ| vs median_cn + signed mean_delta vs median_cn (direction)
  (b) PARTIAL Spearman given coverage: rank-OLS residualize |Δβ| on n_paired_cpg,
      then Spearman(residual, median_cn)  [mirrors 30_H3:332-343]
  (c) within-coverage-stratum Spearman: cov bins <20, 20-40, 40-80, >=80
  (d) box of |Δβ| by integer median_cn

Single sample: HCC1395 paired_full. A-pilot. CN ground truth = SEQC2 (Masood 2024).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
MASTER = ROOT / "genome_survey_v2/cn_confound/master_o1_cn.tsv"
OUT_JSON = ROOT / "genome_survey_v2/cn_confound/q2_dose_response.json"
FIG = ROOT / "figures/cn_confound/q2_dose_response.png"
OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
FIG.parent.mkdir(parents=True, exist_ok=True)

# pre-reg thresholds
RHO_FALSIFIER = 0.2     # partial-rho > 0.2 AND p<0.05 same dir => CN-driven (H-A REFUTED)
RHO_SUPPORT = 0.1       # partial-rho <= 0.1 OR opposite sign => H-A SUPPORTED
ALPHA = 0.05

COV_BINS = [("<20", 0, 20), ("20-40", 20, 40), ("40-80", 40, 80), (">=80", 80, np.inf)]


def fmt(x):
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return None
    return float(x)


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 4:
        return None, None, len(x)
    rho, p = stats.spearmanr(x, y)
    return fmt(rho), fmt(p), int(len(x))


def partial_spearman_cov(sub):
    """rank-OLS residualize |Δβ| on coverage (n_paired_cpg), then Spearman(resid, median_cn).
    Mirrors 30_H3.py:332-343 (rank-transform all, lstsq, residual, spearman vs CN)."""
    cc = sub.dropna(subset=["abs_delta", "n_paired_cpg", "median_cn"]).copy()
    if len(cc) <= 20:
        return {"rho": None, "p": None, "n": int(len(cc)),
                "note": "n<=20, partial Spearman skipped"}
    rk = cc[["abs_delta", "n_paired_cpg", "median_cn"]].rank()
    X = np.column_stack([np.ones(len(rk)), rk.n_paired_cpg.values])
    beta, *_ = np.linalg.lstsq(X, rk.abs_delta.values, rcond=None)
    resid = rk.abs_delta.values - X @ beta
    rho, p = stats.spearmanr(rk.median_cn.values, resid)
    return {"rho": fmt(rho), "p": fmt(p), "n": int(len(cc)),
            "note": "rank-OLS residualize |Δβ| on n_paired_cpg, then Spearman vs median_cn"}


def within_cov_strata(sub):
    out = {}
    for label, lo, hi in COV_BINS:
        s = sub[(sub.n_paired_cpg >= lo) & (sub.n_paired_cpg < hi)]
        rho, p, n = spearman(s.abs_delta.values, s.median_cn.values)
        out[label] = {"rho": rho, "p": p, "n": n}
    return out


def axis_block(sub, label):
    # (a) raw
    rho_abs, p_abs, n = spearman(sub.abs_delta.values, sub.median_cn.values)
    rho_sgn, p_sgn, _ = spearman(sub.mean_delta.values, sub.median_cn.values)
    # (b) partial
    partial = partial_spearman_cov(sub)
    # (c) strata
    strata = within_cov_strata(sub)
    # |db| by integer CN summary
    by_cn = {}
    for cn, g in sub.groupby("median_cn"):
        by_cn[int(cn)] = {
            "n": int(len(g)),
            "median_abs_delta": fmt(g.abs_delta.median()),
            "mean_abs_delta": fmt(g.abs_delta.mean()),
            "median_signed_delta": fmt(g.mean_delta.median()),
        }
    return {
        "axis_type": label,
        "n_sig": int(len(sub)),
        "raw_spearman_absdelta_vs_cn": {"rho": rho_abs, "p": p_abs, "n": n},
        "raw_spearman_signeddelta_vs_cn": {"rho": rho_sgn, "p": p_sgn, "n": n,
                                           "note": "direction: negative=more hypomethylation at higher CN"},
        "partial_spearman_absdelta_vs_cn_given_cov": partial,
        "within_coverage_strata_spearman_absdelta_vs_cn": strata,
        "abs_delta_by_integer_cn": by_cn,
    }


def verdict_for(partial_hp):
    rho = partial_hp.get("rho")
    p = partial_hp.get("p")
    if rho is None or p is None:
        return "INCONCLUSIVE", "partial-rho unavailable (n<=20)"
    # Falsifier: rho > 0.2 (positive, dose-response direction) AND p<0.05
    if rho > RHO_FALSIFIER and p < ALPHA:
        return "REFUTE", (f"partial-rho={rho:.4f} > {RHO_FALSIFIER} AND p={p:.3g} < {ALPHA} "
                          f"=> CN-driven dose-response confirmed; H-A REFUTED")
    # Support: rho <= 0.1 OR opposite (negative) sign
    if rho <= RHO_SUPPORT or rho < 0:
        return "CONFIRM", (f"partial-rho={rho:.4f} <= {RHO_SUPPORT} or opposite sign "
                           f"=> NOT CN-driven; H-A SUPPORTED (HP-axis |Δβ| independent of CN)")
    # gray zone 0.1 < rho <= 0.2
    return "INCONCLUSIVE", (f"partial-rho={rho:.4f} in gray zone ({RHO_SUPPORT}, {RHO_FALSIFIER}] "
                            f"(p={p:.3g}); below falsifier but above support threshold")


def make_figure(hp, allele, df):
    fp = matplotlib.font_manager.FontProperties(
        fname="/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf")
    fig, axes = plt.subplots(2, 4, figsize=(22, 11))
    fig.suptitle("Q2 CN Dose-Response (pre-reg H-A): HP-axis ASM |Delta-beta| vs ground-truth median CN\n"
                 "HCC1395 paired_full (single sample) | sig sites wilcoxon_p<0.05 | SEQC2 CN ground truth",
                 fontsize=14, fontweight="bold")

    for row, (label, sub, block) in enumerate([("HP", df[df.axis_type == "HP"], hp),
                                               ("ALLELE", df[df.axis_type == "ALLELE"], allele)]):
        # col 0: |db| vs CN box by integer CN
        ax = axes[row, 0]
        cns = sorted(sub.median_cn.unique())
        data = [sub[sub.median_cn == c].abs_delta.values for c in cns]
        ax.boxplot(data, positions=cns, widths=0.6, showfliers=False)
        ax.set_xlabel("median CN")
        ax.set_ylabel("|Delta-beta|")
        rr = block["raw_spearman_absdelta_vs_cn"]
        ax.set_title(f"{label}: |db| by integer CN\nraw rho={rr['rho']:.3f} p={rr['p']:.2g} n={rr['n']}")
        ax.grid(alpha=0.3)

        # col 1: partial residual vs CN (residualized on coverage)
        ax = axes[row, 1]
        cc = sub.dropna(subset=["abs_delta", "n_paired_cpg", "median_cn"]).copy()
        if len(cc) > 20:
            rk = cc[["abs_delta", "n_paired_cpg", "median_cn"]].rank()
            X = np.column_stack([np.ones(len(rk)), rk.n_paired_cpg.values])
            beta, *_ = np.linalg.lstsq(X, rk.abs_delta.values, rcond=None)
            resid = rk.abs_delta.values - X @ beta
            jit = cc.median_cn.values + np.random.RandomState(0).uniform(-0.15, 0.15, len(cc))
            ax.scatter(jit, resid, s=6, alpha=0.25, color="tab:blue")
            # trend line: mean residual per CN
            for c in sorted(cc.median_cn.unique()):
                m = resid[cc.median_cn.values == c].mean()
                ax.plot([c - 0.3, c + 0.3], [m, m], color="red", lw=2)
        pb = block["partial_spearman_absdelta_vs_cn_given_cov"]
        prho = pb["rho"] if pb["rho"] is not None else float("nan")
        pp = pb["p"] if pb["p"] is not None else float("nan")
        ax.set_xlabel("median CN")
        ax.set_ylabel("|db| rank residual (cov-controlled)")
        ax.set_title(f"{label}: PARTIAL (cov-controlled)\npartial rho={prho:.3f} p={pp:.2g} n={pb['n']}")
        ax.grid(alpha=0.3)

        # col 2: coverage-colored scatter |db| vs CN
        ax = axes[row, 2]
        jit = sub.median_cn.values + np.random.RandomState(1).uniform(-0.2, 0.2, len(sub))
        cov = np.clip(sub.n_paired_cpg.values, 0, 120)
        sc = ax.scatter(jit, sub.abs_delta.values, c=cov, s=7, alpha=0.4, cmap="viridis")
        plt.colorbar(sc, ax=ax, label="n_paired_cpg (cov, clip120)")
        ax.set_xlabel("median CN")
        ax.set_ylabel("|Delta-beta|")
        ax.set_title(f"{label}: |db| vs CN (coverage-colored)")
        ax.grid(alpha=0.3)

        # col 3: signed db vs CN (direction)
        ax = axes[row, 3]
        cns = sorted(sub.median_cn.unique())
        med_sgn = [sub[sub.median_cn == c].mean_delta.median() for c in cns]
        ax.bar(cns, med_sgn, color="tab:orange", alpha=0.7)
        ax.axhline(0, color="k", lw=0.8)
        sg = block["raw_spearman_signeddelta_vs_cn"]
        ax.set_xlabel("median CN")
        ax.set_ylabel("median signed Delta-beta")
        ax.set_title(f"{label}: signed db direction\nraw rho={sg['rho']:.3f} p={sg['p']:.2g}")
        ax.grid(alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG, dpi=130)
    plt.close(fig)


def main():
    df = pd.read_csv(MASTER, sep="\t")
    # sig sites only
    sig = df[df.wilcoxon_p < ALPHA].copy()
    hp_sub = sig[sig.axis_type == "HP"].copy()
    allele_sub = sig[sig.axis_type == "ALLELE"].copy()

    hp = axis_block(hp_sub, "HP")
    allele = axis_block(allele_sub, "ALLELE")

    verdict, rule_result = verdict_for(hp["partial_spearman_absdelta_vs_cn_given_cov"])

    make_figure(hp, allele, sig)

    out = {
        "question": "Q2: Is HP-axis ASM |Delta-beta| independent of ground-truth median CN (CN dose-response)?",
        "hypothesis": "H-A: HP-axis ASM |Delta-beta| is INDEPENDENT of ground-truth median CN",
        "pre_registration": {
            "falsifier": f"coverage-controlled (partial) Spearman rho(|Δβ|, median_cn) > {RHO_FALSIFIER} AND same-direction p<{ALPHA}",
            "support_threshold": f"partial-rho <= {RHO_SUPPORT} OR opposite (negative) sign => NOT CN-driven (H-A supported)",
            "verdict_driver": "PARTIAL (coverage-controlled) Spearman on HP-axis among sig sites",
        },
        "data": {
            "sample": "HCC1395 paired_full (single sample)",
            "cn_ground_truth": "SEQC2 Masood 2024 (short-read NGS, median CN; default 2 diploid/cnLOH)",
            "master": str(MASTER),
            "n_total": int(len(df)),
            "n_sig": int(len(sig)),
            "n_sig_HP": int(len(hp_sub)),
            "n_sig_ALLELE": int(len(allele_sub)),
            "sig_definition": "wilcoxon_p < 0.05",
        },
        "HP": hp,
        "ALLELE": allele,
        "verdict": verdict,
        "decision_rule_result": rule_result,
        "outputs": {"json": str(OUT_JSON), "figure": str(FIG)},
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    # console summary
    print("=== Q2 CN dose-response (H-A) ===")
    print(f"n_sig HP={len(hp_sub)} ALLELE={len(allele_sub)}")
    print("HP raw    |db|~CN :", hp["raw_spearman_absdelta_vs_cn"])
    print("HP signed db ~CN  :", hp["raw_spearman_signeddelta_vs_cn"])
    print("HP PARTIAL|db|~CN :", hp["partial_spearman_absdelta_vs_cn_given_cov"])
    print("HP strata         :", hp["within_coverage_strata_spearman_absdelta_vs_cn"])
    print("ALLELE PARTIAL    :", allele["partial_spearman_absdelta_vs_cn_given_cov"])
    print("VERDICT:", verdict)
    print("RULE  :", rule_result)
    print("JSON  :", OUT_JSON)
    print("FIG   :", FIG)


if __name__ == "__main__":
    main()
