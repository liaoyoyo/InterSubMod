#!/usr/bin/env python3
"""
C1 -- Context TP/FP on FULL master (genome-wide, per-position).

'See the forest': comprehensive contingency of TP-rate (is_tp mean) by structural
context {LOH x CN-class x coverage}, plus logistic regression and a targeted Fisher
exact test (LOH-abnormalCN vs LOH-cnLOH FP-enrichment).

PER-POSITION input (independence safe; avoids HP1/HP2 pseudoreplication):
  master_perpos.tsv  (34,154 positions; TP 30,511 + FP 3,643)
  fallback: dedup master_o1_cn.tsv HP-axis per somatic_pos.

Coverage per position = max(hp_n_cpg_max, allele_n_cpg) -> n_cpg (all positions have >=1 axis).

NOTE on cn_class/cnLOH coupling (verified in data, not assumed):
  Within LOH, cn_class==neutral (CN=2) == cnloh_flag==1 exactly (6,486 positions).
  So 'LOH & cnLOH' == 'LOH & neutral CN'; 'LOH & abnormal-CN' == gain|loss within LOH.
  cn_class==neutral on the nonLOH side is ordinary diploid (1,227 positions, cnloh_flag=0).

Outputs:
  genome_survey_v2/cn_confound/c1_context_tpfp_full.json
  figures/cn_confound/c1_context_tpfp_full.png

A pilot; single sample (HCC1395); FP small (3,643 positions) -> per-cell n reported,
exact tests where cells small, underpowered cells flagged.
"""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# ---- paths ----
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CN_DIR = os.path.join(ROOT, "genome_survey_v2", "cn_confound")
FIG_DIR = os.path.join(ROOT, "figures", "cn_confound")
os.makedirs(CN_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
PERPOS = os.path.join(CN_DIR, "master_perpos.tsv")
O1 = os.path.join(CN_DIR, "master_o1_cn.tsv")
OUT_JSON = os.path.join(CN_DIR, "c1_context_tpfp_full.json")
OUT_FIG = os.path.join(FIG_DIR, "c1_context_tpfp_full.png")

CJK = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
fp_cjk = FontProperties(fname=CJK) if os.path.exists(CJK) else None

CN_ORDER = ["loss", "neutral", "gain"]
LOH_ORDER = ["LOH", "nonLOH"]
COV_BINS = [0, 20, 40, 80, np.inf]
COV_LABELS = ["<20", "20-40", "40-80", ">=80"]


def load_perpos():
    """Load per-position master; fallback to dedup HP-axis of o1 if perpos missing."""
    if os.path.exists(PERPOS):
        df = pd.read_csv(PERPOS, sep="\t")
        df["n_cpg"] = df[["hp_n_cpg_max", "allele_n_cpg"]].max(axis=1)
        src = "master_perpos.tsv"
    else:
        o = pd.read_csv(O1, sep="\t")
        hp = o[o["axis_type"] == "HP"].copy()
        # dedup: one row per somatic_pos (max n_paired_cpg = most-covered HP axis)
        hp = hp.sort_values("n_paired_cpg", ascending=False).drop_duplicates("somatic_pos")
        hp = hp.rename(columns={"n_paired_cpg": "n_cpg"})
        df = hp[["chrom", "somatic_pos", "is_tp", "loh_status",
                 "median_cn", "cn_class", "cnloh_flag", "n_cpg"]].copy()
        src = "master_o1_cn.tsv (HP-axis dedup fallback)"
    # normalize types
    df["is_tp"] = df["is_tp"].astype(int)
    df["cn_class"] = pd.Categorical(df["cn_class"], categories=CN_ORDER, ordered=True)
    df["loh_status"] = pd.Categorical(df["loh_status"], categories=LOH_ORDER, ordered=True)
    df["cov_bin"] = pd.cut(df["n_cpg"], bins=COV_BINS, labels=COV_LABELS, right=False)
    return df, src


def cell_stats(sub):
    n = int(len(sub))
    if n == 0:
        return {"n": 0, "n_tp": 0, "n_fp": 0, "tp_rate": None}
    n_tp = int(sub["is_tp"].sum())
    return {"n": n, "n_tp": n_tp, "n_fp": n - n_tp, "tp_rate": round(n_tp / n, 4)}


def main():
    df, src = load_perpos()
    overall = cell_stats(df)
    res = {
        "question": ("C1: genome-wide per-position TP-rate (is_tp mean) contingency by "
                     "structural context {LOH x CN-class x coverage}; which contexts shift "
                     "TP-rate vs baseline; within LOH does abnormal-CN show more FP-enrichment "
                     "than cnLOH(CN=2)?"),
        "data_source": src,
        "sample": "HCC1395 (single sample)",
        "n_positions": int(len(df)),
        "overall_tp_rate": overall,
        "design_note": ("per-position dedup (HP1/HP2 pseudoreplication removed). "
                        "Within LOH, cn_class==neutral(CN=2) == cnLOH (verified 1:1)."),
    }

    # ---------- (a) contingency tables ----------
    # 2-way: LOH x CN
    tab_loh_cn = {}
    for loh in LOH_ORDER:
        tab_loh_cn[loh] = {}
        for cn in CN_ORDER:
            sub = df[(df.loh_status == loh) & (df.cn_class == cn)]
            tab_loh_cn[loh][cn] = cell_stats(sub)

    # 2-way: LOH x coverage
    tab_loh_cov = {}
    for loh in LOH_ORDER:
        tab_loh_cov[loh] = {}
        for cv in COV_LABELS:
            sub = df[(df.loh_status == loh) & (df.cov_bin == cv)]
            tab_loh_cov[loh][cv] = cell_stats(sub)

    # marginal: CN
    marg_cn = {cn: cell_stats(df[df.cn_class == cn]) for cn in CN_ORDER}
    # marginal: LOH
    marg_loh = {loh: cell_stats(df[df.loh_status == loh]) for loh in LOH_ORDER}
    # marginal: coverage
    marg_cov = {cv: cell_stats(df[df.cov_bin == cv]) for cv in COV_LABELS}

    # 3-way: LOH x CN x coverage
    tab_3way = {}
    for loh in LOH_ORDER:
        tab_3way[loh] = {}
        for cn in CN_ORDER:
            tab_3way[loh][cn] = {}
            for cv in COV_LABELS:
                sub = df[(df.loh_status == loh) & (df.cn_class == cn) & (df.cov_bin == cv)]
                tab_3way[loh][cn][cv] = cell_stats(sub)

    res["a_tables"] = {
        "marginal_loh": marg_loh,
        "marginal_cn": marg_cn,
        "marginal_coverage": marg_cov,
        "twoway_loh_x_cn": tab_loh_cn,
        "twoway_loh_x_coverage": tab_loh_cov,
        "threeway_loh_x_cn_x_coverage": tab_3way,
    }

    # ---------- (b) logistic regression ----------
    # is_tp ~ C(loh)*C(cn) + coverage. Baseline = nonLOH (ref), neutral CN (ref).
    # Drop loss-cell sparsity worries by keeping all 3 cn classes; model reports each coef.
    dfm = df.copy()
    dfm["loh_status"] = dfm["loh_status"].astype(str)
    dfm["cn_class"] = dfm["cn_class"].astype(str)
    # set explicit references via Treatment
    formula = ("is_tp ~ C(loh_status, Treatment(reference='nonLOH')) * "
               "C(cn_class, Treatment(reference='neutral')) + n_cpg")
    logit_res = {}
    try:
        model = smf.glm(formula=formula, data=dfm,
                        family=sm.families.Binomial()).fit()
        coefs = []
        for name in model.params.index:
            coefs.append({
                "term": name,
                "coef": float(model.params[name]),
                "std_err": float(model.bse[name]),
                "z": float(model.tvalues[name]),
                "p": float(model.pvalues[name]),
                "or": float(np.exp(model.params[name])),
                "ci_low": float(np.exp(model.conf_int().loc[name, 0])),
                "ci_high": float(np.exp(model.conf_int().loc[name, 1])),
                "sig_0.05": bool(model.pvalues[name] < 0.05),
            })
        logit_res = {
            "formula": formula,
            "baseline": "nonLOH x neutral-CN",
            "n": int(model.nobs),
            "llf": float(model.llf),
            "coefficients": coefs,
            "significant_contexts": [c["term"] for c in coefs
                                     if c["sig_0.05"] and c["term"] != "Intercept"],
        }
    except Exception as e:
        logit_res = {"error": str(e), "formula": formula}
    res["b_logistic"] = logit_res

    # ---------- (c) Fisher: LOH-abnormalCN vs LOH-cnLOH ----------
    # Within LOH: abnormal-CN = gain|loss (cn_class != neutral); cnLOH = neutral (CN=2).
    loh = df[df.loh_status == "LOH"]
    abn = loh[loh.cn_class != "neutral"]    # gainLOH + lossLOH
    cnloh = loh[loh.cn_class == "neutral"]  # == cnloh_flag==1
    a_tp, a_fp = int(abn.is_tp.sum()), int((abn.is_tp == 0).sum())
    c_tp, c_fp = int(cnloh.is_tp.sum()), int((cnloh.is_tp == 0).sum())
    # 2x2: rows = [abnormalCN, cnLOH], cols = [FP, TP]
    # OR for FP-enrichment of abnormalCN relative to cnLOH:
    #   OR = (abn_FP * cnloh_TP) / (abn_TP * cnloh_FP)
    table = [[a_fp, a_tp], [c_fp, c_tp]]
    orr, p = fisher_exact(table, alternative="two-sided")
    res["c_fisher_loh_abnormalCN_vs_cnLOH"] = {
        "question": "Within LOH, is abnormal-CN (gainLOH/lossLOH) MORE FP-enriched than cnLOH(CN=2)?",
        "abnormalCN": {"n": a_tp + a_fp, "tp": a_tp, "fp": a_fp,
                       "fp_rate": round(a_fp / (a_tp + a_fp), 4) if (a_tp + a_fp) else None,
                       "tp_rate": round(a_tp / (a_tp + a_fp), 4) if (a_tp + a_fp) else None},
        "cnLOH": {"n": c_tp + c_fp, "tp": c_tp, "fp": c_fp,
                  "fp_rate": round(c_fp / (c_tp + c_fp), 4) if (c_tp + c_fp) else None,
                  "tp_rate": round(c_tp / (c_tp + c_fp), 4) if (c_tp + c_fp) else None},
        "table_rows_abn_cnloh_cols_FP_TP": table,
        "odds_ratio_FP_enrichment_abn_vs_cnloh": float(orr),
        "p_value": float(p),
        "interpretation": ("OR>1 & p<0.05 => abnormal-CN MORE FP-enriched than cnLOH; "
                           "OR<1 => cnLOH more FP-enriched; n.s. => no difference"),
    }

    # ---------- verdict synthesis ----------
    # Which contexts shift TP-rate: use logistic sig terms + flag notable cells.
    notable = []
    for loh in LOH_ORDER:
        for cn in CN_ORDER:
            c = tab_loh_cn[loh][cn]
            if c["n"] >= 30:  # only call out cells with some power
                notable.append((f"{loh}.{cn}", c["tp_rate"], c["n"], c["n_fp"]))
    notable.sort(key=lambda x: (x[1] if x[1] is not None else 1))
    res["verdict"] = build_verdict(res, notable)
    res["caveats"] = [
        "Single sample (HCC1395) -> no cross-sample generalization.",
        "FP small: 3,643 FP positions genome-wide; loss-CN cells especially sparse "
        "(LOH.loss n=%d, nonLOH.loss n=%d) -> Fisher/exact used, underpowered cells flagged."
        % (tab_loh_cn["LOH"]["loss"]["n"], tab_loh_cn["nonLOH"]["loss"]["n"]),
        "cn_class is tumor median_cn (gain dominates, n=26,029); 'neutral'=CN2. Within LOH, "
        "neutral==cnLOH (1:1); on nonLOH side neutral=ordinary diploid.",
        "TP-rate here is ClairS-TO precision proxy on the ASM-eligible position set (positions "
        "with >=1 paired-CpG axis), NOT genome-wide caller precision.",
        "is_tp definition: SEQC2 HighConf sSNV truth membership (TP=30,511 analyzed; FP=4,842 "
        "caller FP of which 3,643 land on ASM-eligible positions).",
    ]

    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=2)
    print("wrote", OUT_JSON)

    make_figure(df, tab_loh_cn, tab_loh_cov, logit_res, res, OUT_FIG)
    print("wrote", OUT_FIG)

    # ---- console compact summary ----
    print("\n=== C1 SUMMARY ===")
    print("overall TP-rate %.4f (n=%d, FP=%d)" %
          (overall["tp_rate"], overall["n"], overall["n_fp"]))
    print("LOH x CN TP-rate table:")
    for loh in LOH_ORDER:
        row = "  %-7s " % loh
        for cn in CN_ORDER:
            c = tab_loh_cn[loh][cn]
            row += "%s=%.3f(n=%d,fp=%d)  " % (cn, c["tp_rate"] if c["tp_rate"] else 0,
                                              c["n"], c["n_fp"])
        print(row)
    fr = res["c_fisher_loh_abnormalCN_vs_cnLOH"]
    print("Fisher LOH-abnormalCN vs cnLOH FP-enrich: OR=%.3f p=%.3g  abnFP%%=%.3f cnlohFP%%=%.3f"
          % (fr["odds_ratio_FP_enrichment_abn_vs_cnloh"], fr["p_value"],
             fr["abnormalCN"]["fp_rate"], fr["cnLOH"]["fp_rate"]))
    print("logistic sig contexts:", logit_res.get("significant_contexts"))


def build_verdict(res, notable):
    sig = res["b_logistic"].get("significant_contexts", [])
    fr = res["c_fisher_loh_abnormalCN_vs_cnLOH"]
    parts = []
    base = res["overall_tp_rate"]["tp_rate"]
    parts.append("Overall TP-rate=%.3f (n=%d). " % (base, res["n_positions"]))
    if sig:
        parts.append("Logistic: contexts significantly shifting TP-rate vs (nonLOH x neutral) = "
                     + "; ".join(sig) + ". ")
    else:
        parts.append("Logistic: no context significantly shifts TP-rate. ")
    # describe direction of LOH main effect & coverage
    for c in res["b_logistic"].get("coefficients", []):
        t = c["term"]
        if "loh_status" in t and "cn_class" not in t and c["sig_0.05"]:
            dirn = "raises" if c["coef"] > 0 else "lowers"
            parts.append("LOH main effect %s TP-rate (OR=%.2f, p=%.2g). " % (dirn, c["or"], c["p"]))
        if t == "n_cpg" and c["sig_0.05"]:
            dirn = "raises" if c["coef"] > 0 else "lowers"
            parts.append("Higher coverage %s TP-rate (OR/CpG=%.3f, p=%.2g). " % (dirn, c["or"], c["p"]))
    # Fisher c
    if fr["p_value"] < 0.05:
        if fr["odds_ratio_FP_enrichment_abn_vs_cnloh"] > 1:
            parts.append("Within LOH, abnormal-CN is MORE FP-enriched than cnLOH "
                         "(OR=%.2f, p=%.2g). " % (fr["odds_ratio_FP_enrichment_abn_vs_cnloh"], fr["p_value"]))
        else:
            parts.append("Within LOH, cnLOH is MORE FP-enriched than abnormal-CN "
                         "(OR=%.2f, p=%.2g). " % (fr["odds_ratio_FP_enrichment_abn_vs_cnloh"], fr["p_value"]))
    else:
        parts.append("Within LOH, abnormal-CN vs cnLOH FP-enrichment NOT significant "
                     "(OR=%.2f, p=%.2g). " % (fr["odds_ratio_FP_enrichment_abn_vs_cnloh"], fr["p_value"]))
    return "".join(parts)


def make_figure(df, tab_loh_cn, tab_loh_cov, logit_res, res, out):
    fig = plt.figure(figsize=(16, 6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.3], wspace=0.32)

    # --- panel 1: heatmap TP-rate LOH x CN ---
    ax1 = fig.add_subplot(gs[0, 0])
    M = np.full((len(LOH_ORDER), len(CN_ORDER)), np.nan)
    for i, loh in enumerate(LOH_ORDER):
        for j, cn in enumerate(CN_ORDER):
            c = tab_loh_cn[loh][cn]
            if c["tp_rate"] is not None:
                M[i, j] = c["tp_rate"]
    im = ax1.imshow(M, cmap="RdYlGn", vmin=0.5, vmax=1.0, aspect="auto")
    ax1.set_xticks(range(len(CN_ORDER)))
    ax1.set_xticklabels(CN_ORDER)
    ax1.set_yticks(range(len(LOH_ORDER)))
    ax1.set_yticklabels(LOH_ORDER)
    ax1.set_xlabel("CN class (tumor median_cn)")
    ax1.set_ylabel("LOH status")
    ax1.set_title("TP-rate by LOH x CN-class")
    for i, loh in enumerate(LOH_ORDER):
        for j, cn in enumerate(CN_ORDER):
            c = tab_loh_cn[loh][cn]
            if c["tp_rate"] is not None:
                txt = "%.3f\nn=%d\nfp=%d" % (c["tp_rate"], c["n"], c["n_fp"])
                col = "black" if M[i, j] > 0.7 else "white"
                ax1.text(j, i, txt, ha="center", va="center", fontsize=8, color=col)
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label="TP-rate")

    # --- panel 2: coverage panel (TP-rate by coverage bin, split by LOH) ---
    ax2 = fig.add_subplot(gs[0, 1])
    x = np.arange(len(COV_LABELS))
    w = 0.38
    for k, loh in enumerate(LOH_ORDER):
        rates = [tab_loh_cov[loh][cv]["tp_rate"] if tab_loh_cov[loh][cv]["tp_rate"] is not None else 0
                 for cv in COV_LABELS]
        ns = [tab_loh_cov[loh][cv]["n"] for cv in COV_LABELS]
        bars = ax2.bar(x + (k - 0.5) * w, rates, w, label=loh,
                       color=["#3b78c2", "#c2563b"][k], alpha=0.85)
        for b, n in zip(bars, ns):
            ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.005,
                     "n=%d" % n, ha="center", va="bottom", fontsize=6, rotation=90)
    ax2.set_xticks(x)
    ax2.set_xticklabels(COV_LABELS)
    ax2.set_xlabel("Coverage bin (n paired-CpG)")
    ax2.set_ylabel("TP-rate")
    ax2.set_ylim(0, 1.05)
    ax2.axhline(res["overall_tp_rate"]["tp_rate"], ls="--", color="gray", lw=1,
                label="overall %.3f" % res["overall_tp_rate"]["tp_rate"])
    ax2.set_title("TP-rate by coverage x LOH")
    ax2.legend(fontsize=8, loc="lower right")

    # --- panel 3: forest of logistic coefficients (OR) ---
    ax3 = fig.add_subplot(gs[0, 2])
    coefs = [c for c in logit_res.get("coefficients", []) if c["term"] != "Intercept"]
    # nicer labels
    def pretty(t):
        t = t.replace("C(loh_status, Treatment(reference='nonLOH'))", "LOH")
        t = t.replace("C(cn_class, Treatment(reference='neutral'))", "CN")
        t = t.replace("[T.", "=").replace("]", "")
        t = t.replace(":", " x ")
        return t
    labels = [pretty(c["term"]) for c in coefs]
    ors = [c["or"] for c in coefs]
    los = [c["ci_low"] for c in coefs]
    his = [c["ci_high"] for c in coefs]
    yy = np.arange(len(coefs))[::-1]
    for y, c, lo, hi, oo in zip(yy, coefs, los, his, ors):
        col = "#c2563b" if c["sig_0.05"] else "#888888"
        # clip extreme CI for display
        hi_d = min(hi, 50)
        ax3.plot([lo, hi_d], [y, y], color=col, lw=2)
        ax3.plot(oo, y, "o", color=col, ms=6)
        star = "*" if c["sig_0.05"] else ""
        ax3.text(min(hi_d, 50) * 1.05, y, "OR=%.2f%s p=%.1e" % (oo, star, c["p"]),
                 va="center", fontsize=7)
    ax3.axvline(1.0, ls="--", color="gray", lw=1)
    ax3.set_yticks(yy)
    ax3.set_yticklabels(labels, fontsize=8)
    ax3.set_xscale("log")
    ax3.set_xlabel("Odds ratio (is_tp), log scale; ref = nonLOH x neutral-CN")
    ax3.set_title("Logistic coefficients (forest)\nis_tp ~ LOH*CN + coverage")
    ax3.set_xlim(0.1, 120)

    fig.suptitle("C1 -- Context TP/FP on FULL master (HCC1395, per-position, n=%d; FP=%d)"
                 % (res["n_positions"], res["overall_tp_rate"]["n_fp"]),
                 fontsize=13, y=1.02)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
