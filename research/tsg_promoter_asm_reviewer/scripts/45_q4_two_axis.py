#!/usr/bin/env python3
"""Q4 — Two-axis diagnostic (pre-reg H-C): is ALLELE-specific ASM signal an
allele-dosage (copy-number) artifact, or does it agree with the HP-axis?

Pre-reg H-C: HP-axis and ALLELE-axis agree at matched loci.
Falsifier: excess := (ALLELE |Δβ| - HP |Δβ|) INCREASES with median_cn,
  interaction p < 0.05  ->  ALLELE-specific signal = allele-dosage artifact
  (and HP-axis is the clean estimate).

For each somatic_pos that has BOTH a HP-axis row (collapse HP1/HP2 by max |Δβ|
and by mean |Δβ|) AND an ALLELE row (ALT_vs_REF):
  excess = ALLELE_abs_delta - HP_abs_delta   (HP collapsed by MAX, primary)
Then:
  - Spearman(excess, median_cn)
  - OLS excess ~ median_cn  (slope, p)  [statsmodels if available, else scipy]
  - stratify by cn_class (neutral vs gain) -> Mann-Whitney U on excess
  - HP |db| vs ALLELE |db| correlation overall and within cnLOH (cnloh_flag==1)

Verdict per H-C uses the excess~median_cn OLS slope sign + p:
  slope > 0 AND p < 0.05  -> H-C REFUTED (ALLELE = allele-dosage artifact;
                              HP-axis is clean) [falsifier met]
  else                    -> H-C CONFIRMED / INCONCLUSIVE (axes agree;
                              no dosage-driven inflation of ALLELE excess)

Single sample (HCC1395) — A pilot. Detail -> json+png; summary returned by caller.
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from scipy import stats

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
MASTER = WS / "genome_survey_v2/cn_confound/master_o1_cn.tsv"
OUT_JSON = WS / "genome_survey_v2/cn_confound/q4_two_axis.json"
OUT_PNG = WS / "figures/cn_confound/q4_two_axis.png"

HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
ALLELE_AXIS = "ALT_vs_REF"
CJK = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
ALPHA = 0.05

# Labels are ENGLISH (context guidance). Use DejaVu Sans (full Latin glyphs) as
# primary; register CJK as a fallback so any incidental CJK still renders.
try:
    font_manager.fontManager.addfont(CJK)
    _cjk_name = font_manager.FontProperties(fname=CJK).get_name()
except Exception:
    _cjk_name = None
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = (
    ["DejaVu Sans"] + ([_cjk_name] if _cjk_name else [])
)
plt.rcParams["axes.unicode_minus"] = False


def fp(size=11):
    # DejaVu Sans has full Latin glyphs (labels are English); CJK is the
    # rcParams fallback for any incidental non-Latin text.
    return font_manager.FontProperties(family="DejaVu Sans", size=size)


def ols_slope_p(x, y):
    """OLS y ~ x. Try statsmodels (gives proper t-test p); fallback to scipy linregress."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    try:
        import statsmodels.api as sm

        X = sm.add_constant(x)
        m = sm.OLS(y, X).fit()
        return {
            "slope": float(m.params[1]),
            "intercept": float(m.params[0]),
            "p": float(m.pvalues[1]),
            "se": float(m.bse[1]),
            "r2": float(m.rsquared),
            "n": int(len(y)),
            "engine": "statsmodels.OLS",
        }
    except Exception:
        lr = stats.linregress(x, y)
        return {
            "slope": float(lr.slope),
            "intercept": float(lr.intercept),
            "p": float(lr.pvalue),
            "se": float(lr.stderr),
            "r2": float(lr.rvalue ** 2),
            "n": int(len(y)),
            "engine": "scipy.linregress",
        }


def main():
    df = pd.read_csv(MASTER, sep="\t")
    df = df.dropna(subset=["abs_delta"]).copy()

    # ---- Build per-(chrom,somatic_pos) HP-collapsed and ALLELE values ----------
    hp = df[df.axis.isin(HP_AXES)].copy()
    al = df[df.axis == ALLELE_AXIS].copy()

    # HP collapse: per locus take MAX |db| (primary) and MEAN |db| (sensitivity).
    hp_max = (
        hp.groupby(["chrom", "somatic_pos"])["abs_delta"].max().rename("hp_abs_max")
    )
    hp_mean = (
        hp.groupby(["chrom", "somatic_pos"])["abs_delta"].mean().rename("hp_abs_mean")
    )
    hp_n = hp.groupby(["chrom", "somatic_pos"]).size().rename("hp_n_axes")

    # ALLELE: one row per locus expected; if duplicates, take max |db| + carry meta.
    al_g = (
        al.groupby(["chrom", "somatic_pos"])
        .agg(
            al_abs=("abs_delta", "max"),
            median_cn=("median_cn", "first"),
            cn_class=("cn_class", "first"),
            cnloh_flag=("cnloh_flag", "first"),
            is_tp=("is_tp", "first"),
            loh_status=("loh_status", "first"),
        )
        .reset_index()
    )

    m = al_g.merge(hp_max, on=["chrom", "somatic_pos"], how="inner")
    m = m.merge(hp_mean, on=["chrom", "somatic_pos"], how="inner")
    m = m.merge(hp_n, on=["chrom", "somatic_pos"], how="inner")

    # primary excess uses HP MAX collapse
    m["excess"] = m["al_abs"] - m["hp_abs_max"]
    m["excess_meancollapse"] = m["al_abs"] - m["hp_abs_mean"]

    n_both = int(len(m))
    n_al_loci = int(al_g.shape[0])
    n_hp_loci = int(hp_max.shape[0])

    # ---- (1) Spearman excess vs median_cn --------------------------------------
    sp_rho, sp_p = stats.spearmanr(m["excess"], m["median_cn"])

    # ---- (2) OLS excess ~ median_cn --------------------------------------------
    ols = ols_slope_p(m["median_cn"], m["excess"])
    # sensitivity: mean-collapse excess
    ols_meancollapse = ols_slope_p(m["median_cn"], m["excess_meancollapse"])

    # ---- (3) cn_class stratify: neutral vs gain Mann-Whitney on excess ---------
    neu = m.loc[m.cn_class == "neutral", "excess"].values
    gain = m.loc[m.cn_class == "gain", "excess"].values
    if len(neu) > 0 and len(gain) > 0:
        mw_u, mw_p = stats.mannwhitneyu(gain, neu, alternative="two-sided")
        mw = {
            "n_neutral": int(len(neu)),
            "n_gain": int(len(gain)),
            "median_excess_neutral": float(np.median(neu)),
            "median_excess_gain": float(np.median(gain)),
            "mean_excess_neutral": float(np.mean(neu)),
            "mean_excess_gain": float(np.mean(gain)),
            "U": float(mw_u),
            "p": float(mw_p),
            "direction": "gain>neutral" if np.median(gain) > np.median(neu) else "gain<=neutral",
        }
    else:
        mw = {"error": "missing neutral or gain stratum", "n_neutral": int(len(neu)), "n_gain": int(len(gain))}

    # ---- (4) HP vs ALLELE |db| correlation overall and within cnLOH ------------
    pear_all = stats.pearsonr(m["hp_abs_max"], m["al_abs"])
    spear_all = stats.spearmanr(m["hp_abs_max"], m["al_abs"])
    cnloh = m[m.cnloh_flag == 1]
    if len(cnloh) >= 3:
        pear_cnloh = stats.pearsonr(cnloh["hp_abs_max"], cnloh["al_abs"])
        spear_cnloh = stats.spearmanr(cnloh["hp_abs_max"], cnloh["al_abs"])
        cnloh_corr = {
            "n": int(len(cnloh)),
            "pearson_r": float(pear_cnloh[0]),
            "pearson_p": float(pear_cnloh[1]),
            "spearman_rho": float(spear_cnloh[0]),
            "spearman_p": float(spear_cnloh[1]),
            "median_excess": float(cnloh["excess"].median()),
        }
    else:
        cnloh_corr = {"n": int(len(cnloh)), "error": "n<3 cnLOH loci with both axes"}

    # ---- Verdict per H-C (uses OLS slope sign + p) ------------------------------
    slope = ols["slope"]
    slope_p = ols["p"]
    falsifier_met = (slope > 0) and (slope_p < ALPHA)
    if falsifier_met:
        verdict = "REFUTE"
        decision_rule = (
            f"FALSIFIER MET: excess~median_cn OLS slope={slope:+.5f} > 0 with p={slope_p:.3g} < {ALPHA}. "
            "ALLELE excess over HP INCREASES with copy number -> ALLELE-specific signal is an "
            "allele-dosage (copy-number) artifact; HP-axis is the clean estimate."
        )
    else:
        if slope_p < ALPHA and slope < 0:
            verdict = "CONFIRM"
            decision_rule = (
                f"Falsifier NOT met (slope={slope:+.5f}, p={slope_p:.3g}<{ALPHA} but slope<0): "
                "ALLELE excess does not grow with CN (if anything shrinks); axes agree under dosage. H-C supported."
            )
        elif slope_p >= ALPHA:
            verdict = "CONFIRM"
            decision_rule = (
                f"Falsifier NOT met: excess~median_cn OLS slope={slope:+.5f}, p={slope_p:.3g} >= {ALPHA} "
                "(no significant CN dependence of ALLELE excess). Axes agree; no dosage-driven inflation. H-C supported."
            )
        else:
            verdict = "INCONCLUSIVE"
            decision_rule = f"slope={slope:+.5f}, p={slope_p:.3g}; ambiguous."

    result = {
        "question": "Q4 H-C two-axis: does ALLELE-specific ASM excess over HP scale with copy number (allele-dosage artifact)?",
        "sample": "HCC1395",
        "task_type": "A pilot (single sample)",
        "n_loci_both_axes": n_both,
        "n_loci_allele_axis": n_al_loci,
        "n_loci_hp_axis": n_hp_loci,
        "hp_collapse": "MAX of HP1/HP2 |db| (primary); MEAN as sensitivity",
        "alpha": ALPHA,
        "excess_definition": "ALLELE |db| (ALT_vs_REF) - HP |db| (max-collapsed)",
        "excess_summary": {
            "mean": float(m["excess"].mean()),
            "median": float(m["excess"].median()),
            "sd": float(m["excess"].std()),
            "pct_allele_gt_hp": float((m["excess"] > 0).mean() * 100),
        },
        "spearman_excess_vs_median_cn": {"rho": float(sp_rho), "p": float(sp_p)},
        "ols_excess_vs_median_cn": ols,
        "ols_excess_vs_median_cn_meancollapse_sensitivity": ols_meancollapse,
        "cn_class_mannwhitney_excess": mw,
        "hp_vs_allele_absdelta_corr_overall": {
            "n": n_both,
            "pearson_r": float(pear_all[0]),
            "pearson_p": float(pear_all[1]),
            "spearman_rho": float(spear_all[0]),
            "spearman_p": float(spear_all[1]),
        },
        "hp_vs_allele_absdelta_corr_within_cnLOH": cnloh_corr,
        "falsifier_met": bool(falsifier_met),
        "verdict": verdict,
        "decision_rule_result": decision_rule,
        "outputs": {"json": str(OUT_JSON), "figure": str(OUT_PNG)},
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)

    # ---- Figure: 3 panels ------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
    color_map = {"neutral": "#4C72B0", "gain": "#C44E52", "loss": "#55A868"}
    cn_colors = m["cn_class"].map(color_map).fillna("#888888")

    # Panel A: excess vs median_cn (jitter + per-CN mean + OLS line)
    axA = axes[0]
    jit = (np.random.RandomState(0).rand(len(m)) - 0.5) * 0.6
    axA.scatter(m["median_cn"] + jit, m["excess"], s=6, alpha=0.18, c=cn_colors, linewidths=0)
    cn_means = m.groupby("median_cn")["excess"].mean()
    axA.plot(cn_means.index, cn_means.values, "o-", color="black", ms=6, lw=1.8, label="per-CN mean excess", zorder=5)
    xs = np.array([m["median_cn"].min(), m["median_cn"].max()], float)
    axA.plot(xs, ols["intercept"] + ols["slope"] * xs, "--", color="darkorange", lw=2,
             label=f"OLS slope={slope:+.4f}\np={slope_p:.2e}", zorder=6)
    axA.axhline(0, color="gray", lw=0.8, ls=":")
    axA.set_xlabel("median copy number (SEQC2)", fontproperties=fp(12))
    axA.set_ylabel("excess = |db|(ALLELE) - |db|(HP)", fontproperties=fp(12))
    axA.set_title(f"A. ALLELE-over-HP excess vs CN  (n={n_both})\nSpearman rho={sp_rho:+.3f}, p={sp_p:.1e}",
                  fontproperties=fp(12))
    axA.legend(prop=fp(9), loc="best")

    # Panel B: HP |db| vs ALLELE |db| scatter colored by cn_class, y=x line
    axB = axes[1]
    for cls in ["neutral", "gain", "loss"]:
        sub = m[m.cn_class == cls]
        if len(sub):
            axB.scatter(sub["hp_abs_max"], sub["al_abs"], s=7, alpha=0.25,
                        c=color_map[cls], label=f"{cls} (n={len(sub)})", linewidths=0)
    lim = [0, max(m["hp_abs_max"].max(), m["al_abs"].max()) * 1.02]
    axB.plot(lim, lim, "--", color="black", lw=1.2, label="y=x (agreement)")
    axB.set_xlim(lim)
    axB.set_ylim(lim)
    axB.set_xlabel("HP-axis |db| (max of HP1/HP2)", fontproperties=fp(12))
    axB.set_ylabel("ALLELE-axis |db| (ALT_vs_REF)", fontproperties=fp(12))
    axB.set_title(
        f"B. HP vs ALLELE |db|  Pearson r={pear_all[0]:.3f}\n"
        f"within cnLOH r={cnloh_corr.get('pearson_r', float('nan')):.3f} (n={cnloh_corr.get('n')})",
        fontproperties=fp(12),
    )
    axB.legend(prop=fp(9), loc="upper left")

    # Panel C: excess box by cn_class (neutral / gain / loss)
    axC = axes[2]
    order = [c for c in ["neutral", "gain", "loss"] if (m.cn_class == c).any()]
    box_data = [m.loc[m.cn_class == c, "excess"].values for c in order]
    bp = axC.boxplot(box_data, labels=order, showfliers=False, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(color_map[c])
        patch.set_alpha(0.55)
    axC.axhline(0, color="gray", lw=0.8, ls=":")
    mw_p_txt = mw.get("p", float("nan"))
    axC.set_ylabel("excess = |db|(ALLELE) - |db|(HP)", fontproperties=fp(12))
    axC.set_xlabel("cn_class", fontproperties=fp(12))
    axC.set_title(
        f"C. excess by cn_class\nMann-Whitney gain vs neutral p={mw_p_txt:.2e}",
        fontproperties=fp(12),
    )
    for lbl in axC.get_xticklabels():
        lbl.set_fontproperties(fp(11))

    fig.suptitle(
        f"Q4 H-C two-axis diagnostic (HCC1395, A pilot) — verdict: {verdict}",
        fontproperties=fp(14), y=1.02,
    )
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # compact stdout for caller
    print(json.dumps({
        "verdict": verdict,
        "n_both": n_both,
        "ols_slope": slope, "ols_p": slope_p,
        "spearman_rho": float(sp_rho), "spearman_p": float(sp_p),
        "mw_gain_vs_neutral_p": mw.get("p"),
        "pearson_overall": float(pear_all[0]),
        "pearson_cnLOH": cnloh_corr.get("pearson_r"),
        "json": str(OUT_JSON), "png": str(OUT_PNG),
    }, indent=2))


if __name__ == "__main__":
    main()
