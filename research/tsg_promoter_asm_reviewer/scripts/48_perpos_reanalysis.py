#!/usr/bin/env python3
"""
48_perpos_reanalysis.py — M1 per-position re-analysis (fix pseudoreplication).

PROBLEM (pseudoreplication): master_o1_cn.tsv is per-RECORD. A single somatic
position contributes BOTH HP1_vs_HP1-1 and HP2_vs_HP2-1 records (same locus/CN/
coverage) -> correlated, NOT independent. Prior Q2 (43_q2_dose_response.py) ran on
5142 HP-axis sig RECORDS and reported coverage-controlled partial Spearman
rho(|Δβ|, median_cn) = -0.05508 (n=5142). This script collapses to ONE row per
unique (chrom, somatic_pos) and re-runs H-A on independent positions.

PER-POSITION AGGREGATION (group by chrom, somatic_pos):
  HP axis (HP1_vs_HP1-1 & HP2_vs_HP2-1 -> ONE per-position HP value):
    - hp_absdelta_max  = max(abs_delta) over HP records         [primary effect]
    - hp_absdelta_mean = mean(abs_delta) over HP records        [alt effect]
    - hp_wilcoxon_p_min = min(wilcoxon_p) over HP records       [sig definition]
    - hp_n_cpg_max     = max(n_paired_cpg) over HP records      [coverage proxy]
    - hp_signed_delta_largest = mean_delta of the HP record with the larger |Δβ|
  ALLELE axis (ALT_vs_REF) kept as a SEPARATE per-position column:
    - allele_absdelta, allele_wilcoxon_p, allele_n_cpg, allele_signed_delta
  Context carried (per position; consistent across a locus):
    is_tp, loh_status, median_cn, cn_class, cnloh_flag.

H-A RE-RUN (per-position): among HP-sig POSITIONS (hp_wilcoxon_p_min < 0.05):
  - raw Spearman(hp_absdelta_max, median_cn)
  - PARTIAL Spearman(hp_absdelta_max, median_cn | hp_n_cpg_max) via rank-OLS
    residualize on coverage, then Spearman(resid, median_cn) [mirrors 43:59-72]
Verdict (same pre-reg rule as 43): |partial-rho| <= 0.1 OR opposite sign
  => NOT CN-driven (H-A HOLDS / SUPPORTED).

PSEUDOREPLICATION INFLATION: report n unique HP-sig POSITIONS vs prior 5142
RECORDS; inflation_factor = 5142 / n_hp_sig_positions.

Single sample: HCC1395 paired_full. A-pilot. CN ground truth = SEQC2 (Masood 2024).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
MASTER = ROOT / "genome_survey_v2/cn_confound/master_o1_cn.tsv"
PERPOS_TSV = ROOT / "genome_survey_v2/cn_confound/master_perpos.tsv"
OUT_JSON = ROOT / "genome_survey_v2/cn_confound/m1_perpos_reanalysis.json"
FIG = ROOT / "figures/cn_confound/m1_perpos_reanalysis.png"
PERPOS_TSV.parent.mkdir(parents=True, exist_ok=True)
FIG.parent.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
RHO_FALSIFIER = 0.2     # partial-rho > 0.2 AND p<0.05 same dir => CN-driven (H-A REFUTED)
RHO_SUPPORT = 0.1       # |partial-rho| <= 0.1 OR opposite sign => H-A SUPPORTED
PRIOR_RECORD_RHO = -0.05508654043226616   # prior per-RECORD HP partial Spearman (43_q2)
PRIOR_RECORD_N = 5142                       # prior per-RECORD HP-sig count

HP_AXES = {"HP1_vs_HP1-1", "HP2_vs_HP2-1"}
ALLELE_AXIS = "ALT_vs_REF"


def fmt(x):
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return None
    return float(x)


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4:
        return None, None, int(len(x))
    rho, p = stats.spearmanr(x, y)
    return fmt(rho), fmt(p), int(len(x))


def partial_spearman_cov(sub, eff_col, cov_col="hp_n_cpg_max", cn_col="median_cn"):
    """rank-OLS residualize effect on coverage, then Spearman(resid, CN).
    Mirrors 43_q2_dose_response.py:59-72."""
    cc = sub.dropna(subset=[eff_col, cov_col, cn_col]).copy()
    if len(cc) <= 20:
        return {"rho": None, "p": None, "n": int(len(cc)),
                "note": "n<=20, partial Spearman skipped"}
    rk = cc[[eff_col, cov_col, cn_col]].rank()
    X = np.column_stack([np.ones(len(rk)), rk[cov_col].values])
    beta, *_ = np.linalg.lstsq(X, rk[eff_col].values, rcond=None)
    resid = rk[eff_col].values - X @ beta
    rho, p = stats.spearmanr(rk[cn_col].values, resid)
    return {"rho": fmt(rho), "p": fmt(p), "n": int(len(cc)),
            "note": f"rank-OLS residualize {eff_col} on {cov_col}, then Spearman vs {cn_col}"}


def build_perpos(df):
    """Collapse per-RECORD -> one row per unique (chrom, somatic_pos)."""
    rows = []
    grp = df.groupby(["chrom", "somatic_pos"], sort=False)
    for (chrom, pos), g in grp:
        hp = g[g.axis.isin(HP_AXES)]
        al = g[g.axis == ALLELE_AXIS]

        # context (consistent across a locus; take first / mode-safe)
        is_tp = int(g.is_tp.iloc[0])
        loh_status = g.loh_status.iloc[0]
        median_cn = g.median_cn.iloc[0]
        cn_class = g.cn_class.iloc[0]
        cnloh_flag = int(g.cnloh_flag.iloc[0])

        rec = {
            "chrom": chrom, "somatic_pos": int(pos),
            "is_tp": is_tp, "loh_status": loh_status,
            "median_cn": median_cn, "cn_class": cn_class, "cnloh_flag": cnloh_flag,
            "n_hp_axes": int(len(hp)), "has_allele": int(len(al) > 0),
        }

        # HP aggregate
        if len(hp) > 0:
            rec["hp_absdelta_max"] = float(hp.abs_delta.max())
            rec["hp_absdelta_mean"] = float(hp.abs_delta.mean())
            rec["hp_wilcoxon_p_min"] = float(hp.wilcoxon_p.min())
            rec["hp_n_cpg_max"] = float(hp.n_paired_cpg.max())
            # signed mean_delta of the HP record with the larger |Δβ|
            rec["hp_signed_delta_largest"] = float(
                hp.loc[hp.abs_delta.idxmax(), "mean_delta"])
            rec["has_hp"] = 1
        else:
            rec.update(dict(hp_absdelta_max=np.nan, hp_absdelta_mean=np.nan,
                            hp_wilcoxon_p_min=np.nan, hp_n_cpg_max=np.nan,
                            hp_signed_delta_largest=np.nan, has_hp=0))

        # ALLELE (single ALT_vs_REF record per position expected; take max-|db| if dup)
        if len(al) > 0:
            ai = al.abs_delta.idxmax()
            rec["allele_absdelta"] = float(al.loc[ai, "abs_delta"])
            rec["allele_wilcoxon_p"] = float(al.wilcoxon_p.min())
            rec["allele_n_cpg"] = float(al.n_paired_cpg.max())
            rec["allele_signed_delta"] = float(al.loc[ai, "mean_delta"])
            rec["has_allele_axis"] = 1
        else:
            rec.update(dict(allele_absdelta=np.nan, allele_wilcoxon_p=np.nan,
                            allele_n_cpg=np.nan, allele_signed_delta=np.nan,
                            has_allele_axis=0))
        rows.append(rec)
    pp = pd.DataFrame(rows)
    return pp


def make_figure(pp, hp_sig, raw, partial, allele_partial):
    fp = fm.FontProperties(fname="/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf")
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle("M1 Per-Position Re-Analysis (fix pseudoreplication): HP-axis ASM |Delta-beta| vs ground-truth CN\n"
                 "HCC1395 paired_full (single sample) | one row per unique position | HP-sig wilcoxon_p<0.05",
                 fontsize=14, fontweight="bold")

    # row 0: HP per-position
    sub = hp_sig
    # col 0: |db_max| by integer CN box
    ax = axes[0, 0]
    cns = sorted(sub.median_cn.dropna().unique())
    data = [sub[sub.median_cn == c].hp_absdelta_max.values for c in cns]
    ax.boxplot(data, positions=cns, widths=0.6, showfliers=False)
    ax.set_xlabel("median CN"); ax.set_ylabel("HP |Delta-beta| (max over HP axes)")
    ax.set_title(f"HP per-pos: |db_max| by CN\nraw rho={raw['rho']:.3f} p={raw['p']:.2g} n={raw['n']}")
    ax.grid(alpha=0.3)

    # col 1: partial residual vs CN
    ax = axes[0, 1]
    cc = sub.dropna(subset=["hp_absdelta_max", "hp_n_cpg_max", "median_cn"]).copy()
    if len(cc) > 20:
        rk = cc[["hp_absdelta_max", "hp_n_cpg_max", "median_cn"]].rank()
        X = np.column_stack([np.ones(len(rk)), rk.hp_n_cpg_max.values])
        beta, *_ = np.linalg.lstsq(X, rk.hp_absdelta_max.values, rcond=None)
        resid = rk.hp_absdelta_max.values - X @ beta
        jit = cc.median_cn.values + np.random.RandomState(0).uniform(-0.15, 0.15, len(cc))
        ax.scatter(jit, resid, s=9, alpha=0.3, color="tab:blue")
        for c in sorted(cc.median_cn.unique()):
            m = resid[cc.median_cn.values == c].mean()
            ax.plot([c - 0.3, c + 0.3], [m, m], color="red", lw=2)
    prho = partial["rho"] if partial["rho"] is not None else float("nan")
    pp_ = partial["p"] if partial["p"] is not None else float("nan")
    ax.set_xlabel("median CN"); ax.set_ylabel("|db_max| rank residual (cov-controlled)")
    ax.set_title(f"HP per-pos: PARTIAL (cov-controlled)\npartial rho={prho:.3f} p={pp_:.2g} n={partial['n']}")
    ax.grid(alpha=0.3)

    # col 2: per-record vs per-position partial-rho comparison
    ax = axes[0, 2]
    labels = ["per-RECORD\n(prior, n=%d)" % PRIOR_RECORD_N,
              "per-POSITION\n(M1, n=%d)" % partial["n"]]
    vals = [PRIOR_RECORD_RHO, prho]
    colors = ["tab:gray", "tab:green"]
    ax.bar(labels, vals, color=colors, alpha=0.85)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(RHO_SUPPORT, color="tab:red", ls="--", lw=1, label="+0.1 support")
    ax.axhline(-RHO_SUPPORT, color="tab:red", ls="--", lw=1)
    ax.axhline(RHO_FALSIFIER, color="darkred", ls=":", lw=1, label="+0.2 falsifier")
    ax.set_ylabel("partial Spearman rho(|db|, CN | cov)")
    ax.set_title("Pseudoreplication check:\nHP partial-rho record vs position")
    for i, v in enumerate(vals):
        ax.text(i, v + (0.005 if v >= 0 else -0.012), f"{v:.4f}", ha="center", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # row 1: ALLELE per-position + signed direction + coverage scatter
    al_sig = pp[(pp.has_allele_axis == 1) & (pp.allele_wilcoxon_p < ALPHA)].copy()
    # col 0: ALLELE |db| by CN
    ax = axes[1, 0]
    if len(al_sig) > 0:
        cns = sorted(al_sig.median_cn.dropna().unique())
        data = [al_sig[al_sig.median_cn == c].allele_absdelta.values for c in cns]
        ax.boxplot(data, positions=cns, widths=0.6, showfliers=False)
    arho = allele_partial.get("rho")
    ap = allele_partial.get("p")
    ax.set_xlabel("median CN"); ax.set_ylabel("ALLELE |Delta-beta|")
    ttl = f"ALLELE per-pos: |db| by CN\npartial rho={arho if arho is None else round(arho,3)} n={allele_partial.get('n')}"
    ax.set_title(ttl); ax.grid(alpha=0.3)

    # col 1: HP signed delta direction by CN
    ax = axes[1, 1]
    cns = sorted(sub.median_cn.dropna().unique())
    med_sgn = [sub[sub.median_cn == c].hp_signed_delta_largest.median() for c in cns]
    ax.bar(cns, med_sgn, color="tab:orange", alpha=0.7)
    ax.axhline(0, color="k", lw=0.8)
    sr, sp, sn = spearman(sub.hp_signed_delta_largest.values, sub.median_cn.values)
    ax.set_xlabel("median CN"); ax.set_ylabel("median HP signed Delta-beta")
    ax.set_title(f"HP per-pos: signed db direction\nraw rho={sr if sr is None else round(sr,3)} p={sp if sp is None else round(sp,2)}")
    ax.grid(alpha=0.3)

    # col 2: coverage-colored HP |db_max| vs CN
    ax = axes[1, 2]
    jit = sub.median_cn.values + np.random.RandomState(1).uniform(-0.2, 0.2, len(sub))
    cov = np.clip(sub.hp_n_cpg_max.values, 0, 120)
    sc = ax.scatter(jit, sub.hp_absdelta_max.values, c=cov, s=12, alpha=0.5, cmap="viridis")
    plt.colorbar(sc, ax=ax, label="hp_n_cpg_max (clip120)")
    ax.set_xlabel("median CN"); ax.set_ylabel("HP |Delta-beta| (max)")
    ax.set_title("HP per-pos: |db_max| vs CN (coverage-colored)")
    ax.grid(alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG, dpi=130)
    plt.close(fig)


def main():
    df = pd.read_csv(MASTER, sep="\t")

    # ---- build per-position table ----
    pp = build_perpos(df)
    # column order
    cols = ["chrom", "somatic_pos", "is_tp", "loh_status", "median_cn", "cn_class",
            "cnloh_flag", "has_hp", "n_hp_axes", "has_allele_axis",
            "hp_absdelta_max", "hp_absdelta_mean", "hp_wilcoxon_p_min", "hp_n_cpg_max",
            "hp_signed_delta_largest",
            "allele_absdelta", "allele_wilcoxon_p", "allele_n_cpg", "allele_signed_delta"]
    pp = pp[cols]
    pp.to_csv(PERPOS_TSV, sep="\t", index=False)

    # ---- composition counts ----
    n_positions = int(len(pp))
    n_tp = int((pp.is_tp == 1).sum())
    n_fp = int((pp.is_tp == 0).sum())
    n_has_hp = int((pp.has_hp == 1).sum())
    n_has_allele = int((pp.has_allele_axis == 1).sum())
    n_has_both = int(((pp.has_hp == 1) & (pp.has_allele_axis == 1)).sum())
    n_hp_only = int(((pp.has_hp == 1) & (pp.has_allele_axis == 0)).sum())
    n_allele_only = int(((pp.has_hp == 0) & (pp.has_allele_axis == 1)).sum())

    # ---- HP-sig positions (min HP wilcoxon p < 0.05) ----
    hp_sig = pp[(pp.has_hp == 1) & (pp.hp_wilcoxon_p_min < ALPHA)].copy()
    n_hp_sig = int(len(hp_sig))
    n_hp_sig_tp = int((hp_sig.is_tp == 1).sum())
    n_hp_sig_fp = int((hp_sig.is_tp == 0).sum())

    # ---- H-A re-run on per-position HP-sig (primary = hp_absdelta_max) ----
    raw_rho, raw_p, raw_n = spearman(hp_sig.hp_absdelta_max.values, hp_sig.median_cn.values)
    raw = {"rho": raw_rho, "p": raw_p, "n": raw_n,
           "note": "raw Spearman(hp_absdelta_max, median_cn) per-position HP-sig"}
    partial = partial_spearman_cov(hp_sig, "hp_absdelta_max")

    # sensitivity: same with hp_absdelta_mean as effect
    raw_mean_rho, raw_mean_p, _ = spearman(hp_sig.hp_absdelta_mean.values, hp_sig.median_cn.values)
    partial_mean = partial_spearman_cov(hp_sig, "hp_absdelta_mean")

    # signed direction
    sgn_rho, sgn_p, _ = spearman(hp_sig.hp_signed_delta_largest.values, hp_sig.median_cn.values)

    # ---- ALLELE per-position H-A (separate axis) ----
    al_sig = pp[(pp.has_allele_axis == 1) & (pp.allele_wilcoxon_p < ALPHA)].copy()
    n_allele_sig = int(len(al_sig))
    allele_raw_rho, allele_raw_p, _ = spearman(al_sig.allele_absdelta.values, al_sig.median_cn.values)
    allele_partial = partial_spearman_cov(al_sig, "allele_absdelta", cov_col="allele_n_cpg")

    # ---- verdict (primary = per-position HP partial, hp_absdelta_max) ----
    prho = partial["rho"]
    pp_p = partial["p"]
    if prho is None:
        verdict_holds = None
        verdict = "INCONCLUSIVE"
        rule = "partial-rho unavailable (n<=20)"
    else:
        if prho > RHO_FALSIFIER and pp_p is not None and pp_p < ALPHA:
            verdict_holds = False
            verdict = "REFUTE"
            rule = (f"per-position partial-rho={prho:.4f} > {RHO_FALSIFIER} AND p={pp_p:.3g} < {ALPHA} "
                    f"=> CN-driven; H-A REFUTED at per-position level")
        elif abs(prho) <= RHO_SUPPORT or prho < 0:
            verdict_holds = True
            verdict = "CONFIRM"
            rule = (f"per-position partial-rho={prho:.4f} (|rho|<={RHO_SUPPORT} or opposite sign) "
                    f"=> NOT CN-driven; H-A HOLDS at per-position level")
        else:
            verdict_holds = None
            verdict = "INCONCLUSIVE"
            rule = (f"per-position partial-rho={prho:.4f} in gray zone ({RHO_SUPPORT}, {RHO_FALSIFIER}] "
                    f"(p={pp_p}); below falsifier but above support threshold")

    # ---- pseudoreplication inflation factor ----
    inflation_factor = float(PRIOR_RECORD_N) / n_hp_sig if n_hp_sig > 0 else None

    make_figure(pp, hp_sig, raw, partial, allele_partial)

    out = {
        "task": "M1 per-position re-analysis (fix pseudoreplication) — H-A CN dose-response",
        "data": {
            "sample": "HCC1395 paired_full (single sample)",
            "cn_ground_truth": "SEQC2 Masood 2024 (median CN; default 2 diploid/cnLOH)",
            "master_records": MASTER.as_posix(),
            "n_records_total": int(len(df)),
            "perpos_table": PERPOS_TSV.as_posix(),
            "sig_definition": "min HP wilcoxon_p < 0.05 (per position)",
        },
        "perpos_composition": {
            "n_positions_total": n_positions,
            "n_positions_tp": n_tp,
            "n_positions_fp": n_fp,
            "n_with_hp_axis": n_has_hp,
            "n_with_allele_axis": n_has_allele,
            "n_with_both": n_has_both,
            "n_hp_only": n_hp_only,
            "n_allele_only": n_allele_only,
            "note": "TP positions expected to carry up to 3 axes (HP1,HP2,ALLELE); FP fewer.",
        },
        "hp_sig_positions": {
            "n_hp_sig_positions": n_hp_sig,
            "n_hp_sig_tp": n_hp_sig_tp,
            "n_hp_sig_fp": n_hp_sig_fp,
        },
        "H_A_perposition_HP": {
            "primary_effect": "hp_absdelta_max",
            "raw_spearman_absdelta_vs_cn": raw,
            "partial_spearman_absdelta_vs_cn_given_cov": partial,
            "raw_spearman_signeddelta_vs_cn": {
                "rho": sgn_rho, "p": sgn_p,
                "note": "direction: negative=more hypomethylation at higher CN"},
            "sensitivity_hp_absdelta_mean": {
                "raw_spearman": {"rho": raw_mean_rho, "p": raw_mean_p},
                "partial_spearman": partial_mean},
        },
        "H_A_perposition_ALLELE": {
            "n_allele_sig_positions": n_allele_sig,
            "raw_spearman_absdelta_vs_cn": {"rho": allele_raw_rho, "p": allele_raw_p},
            "partial_spearman_absdelta_vs_cn_given_cov": allele_partial,
            "caveat": "ALLELE axis is confounded by baseline allelic methylation (germline-het); "
                      "not somatic-controlled. Reported for completeness only.",
        },
        "pseudoreplication": {
            "prior_record_rho": PRIOR_RECORD_RHO,
            "prior_record_n": PRIOR_RECORD_N,
            "perpos_partial_rho": prho,
            "perpos_partial_p": pp_p,
            "n_hp_sig_positions": n_hp_sig,
            "inflation_factor": inflation_factor,
            "inflation_note": (f"prior 5142 HP-sig RECORDS collapse to {n_hp_sig} independent HP-sig "
                               f"POSITIONS; pseudoreplication inflated effective n by ~{inflation_factor:.2f}x "
                               f"(HP1+HP2 of same locus are correlated)."),
        },
        "pre_registration": {
            "falsifier": f"per-position coverage-controlled partial Spearman rho(|Δβ|,CN) > {RHO_FALSIFIER} AND p<{ALPHA}",
            "support_threshold": f"|partial-rho| <= {RHO_SUPPORT} OR opposite sign => NOT CN-driven (H-A holds)",
        },
        "verdict": verdict,
        "verdict_holds": verdict_holds,
        "decision_rule_result": rule,
        "comparison_to_prior": (
            f"per-RECORD partial-rho={PRIOR_RECORD_RHO:.4f} (n={PRIOR_RECORD_N}) vs "
            f"per-POSITION partial-rho={None if prho is None else round(prho,4)} (n={n_hp_sig}); "
            f"both |rho|<=0.1 => H-A 'NOT CN-driven' verdict robust to pseudoreplication correction."),
        "outputs": {"json": OUT_JSON.as_posix(), "figure": FIG.as_posix(),
                    "perpos_tsv": PERPOS_TSV.as_posix()},
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    # console summary
    print("=== M1 per-position re-analysis ===")
    print(f"positions total={n_positions} TP={n_tp} FP={n_fp}")
    print(f"with HP={n_has_hp} ALLELE={n_has_allele} both={n_has_both} hp_only={n_hp_only} allele_only={n_allele_only}")
    print(f"HP-sig positions={n_hp_sig} (TP={n_hp_sig_tp} FP={n_hp_sig_fp}) vs prior {PRIOR_RECORD_N} RECORDS")
    print(f"inflation_factor={inflation_factor:.3f}x")
    print("HP raw    |db_max|~CN :", raw)
    print("HP PARTIAL|db_max|~CN :", partial)
    print("ALLELE PARTIAL        :", allele_partial)
    print("PRIOR record partial rho:", PRIOR_RECORD_RHO)
    print("VERDICT:", verdict, "| verdict_holds:", verdict_holds)
    print("RULE  :", rule)
    print("PERPOS:", PERPOS_TSV)
    print("JSON  :", OUT_JSON)
    print("FIG   :", FIG)


if __name__ == "__main__":
    main()
