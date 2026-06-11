#!/usr/bin/env python3
"""
101 - Does methylation CLUSTERABILITY (whether a locus's methylation reads cluster, and
how many groups) correlate with READ COUNT (coverage) and COPY NUMBER (CN)?

User question (2026-06-10): 甲基可以分群的數量 vs read 數量 + CN 數量 是否相關 (confound check).

CONFOUND CARE: coverage and CN are themselves correlated (amplified region -> more reads),
so we DISENTANGLE with partial Spearman + within-CN-state coverage corr + within-coverage-bin
CN rate + multivariable logistic (is_clustered ~ coverage + n_CpG + CN). auc-confound-guard logic.

Metrics for "甲基可以分群":
  is_clustered (binary) = reliable (Significant) OR latent (script-88 PERMANOVA recovery)
  NGroups   = HPFineNGroups (⚠ HP-tag sub-family count, capped 4; per HD-4 it is phasing/HP-tag
              driven not pure methylation — reported with that caveat)
  cramersV  = clustering strength (continuous)

Predictors:
  coverage  = NumReads ; n_CpG = NumCpGs (clustering POWER control) ; CN = SEQC2 truth (HCC1395)

CN truth = /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed (gain/loss/loh;
not-in-bed = neutral CN=2). CN-only valid for HCC1395 (SEQC2 = HCC1395 truth).

Output: genome_survey_v2/catalog/clusterability_vs_cov_cn.json  (deterministic; §13 layer-A)
"""
import csv, json, glob, os
import numpy as np
from scipy import stats
import statsmodels.api as sm

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
GS = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2"
CNBED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
LATENT_MINN = 5


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def is_latent(r, minn=LATENT_MINN):
    if tb(r, "Significant") or not tb(r, "PassedGating") or (num(r, "NumReads") or 0) < 20:
        return False
    if not (tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05):
        return False
    if not (tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05):
        return False
    return min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0) >= minn


def is_clustered(r):
    return 1 if (tb(r, "Significant") or is_latent(r)) else 0


# ---- CN bed interval lookup ----
def load_cn():
    cn = {}
    for line in open(CNBED):
        p = line.rstrip().split("\t")
        if len(p) < 4:
            continue
        cn.setdefault(p[0], []).append((int(p[1]), int(p[2]), p[3]))
    for c in cn:
        cn[c].sort()
    return cn


def cn_state(cn, chrom, pos):
    for s, e, lab in cn.get(chrom, []):
        if s <= pos <= e:
            return lab  # gain / loss / loh
    return "neutral"


def load_sample(s, classes=("tp",)):
    rows = []
    for cl in classes:
        f = f"{EX}/{s}_{cl}/significance_summary.csv"
        if os.path.exists(f):
            for r in csv.DictReader(open(f)):
                rows.append(r)
    return rows


def spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return None
    rho, p = stats.spearmanr(x[m], y[m])
    return {"rho": round(float(rho), 4), "p": float(f"{p:.2e}"), "n": int(m.sum())}


def partial_spear(x, y, z):
    """partial Spearman of x,y controlling z: corr of rank-residuals."""
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if m.sum() < 20:
        return None
    rx, ry, rz = (stats.rankdata(v[m]) for v in (x, y, z))
    def resid(a, b):
        b1 = np.c_[np.ones_like(b), b]
        return a - b1 @ np.linalg.lstsq(b1, a, rcond=None)[0]
    ex, ey = resid(rx, rz), resid(ry, rz)
    rho, p = stats.pearsonr(ex, ey)
    return {"partial_rho": round(float(rho), 4), "p": float(f"{p:.2e}"), "n": int(m.sum()), "control": "rank"}


def rate_by_bin(vals, flags, nbin=10):
    vals, flags = np.asarray(vals, float), np.asarray(flags, float)
    m = np.isfinite(vals) & np.isfinite(flags)
    vals, flags = vals[m], flags[m]
    qs = np.quantile(vals, np.linspace(0, 1, nbin + 1))
    out = []
    for i in range(nbin):
        lo, hi = qs[i], qs[i + 1]
        sel = (vals >= lo) & (vals <= hi) if i == nbin - 1 else (vals >= lo) & (vals < hi)
        if sel.sum():
            out.append({"bin": i, "cov_lo": round(float(lo), 1), "cov_hi": round(float(hi), 1),
                        "n": int(sel.sum()), "clustered_rate": round(float(flags[sel].mean()), 4)})
    return out


def main():
    cn = load_cn()
    out = {"question": "methylation clusterability (is_clustered / NGroups / CramersV) vs read-count(coverage) + CN",
           "definitions": {"is_clustered": "reliable(Significant) OR latent(PERMANOVA minN5)",
                           "NGroups": "HPFineNGroups (HP-tag sub-family count, HD-4 phasing-driven caveat)",
                           "coverage": "NumReads", "n_CpG": "NumCpGs", "CN": "SEQC2 gain/loss/loh/neutral (HCC1395)"}}

    # ---------- A. genome-wide coverage correlation (6 samples, TP) ----------
    gw = {}
    for s in SAMPLES:
        rows = load_sample(s, ("tp",))
        if not rows:
            continue
        clu = [is_clustered(r) for r in rows]
        cov = [num(r, "NumReads") for r in rows]
        ncpg = [num(r, "NumCpGs") for r in rows]
        ng = [num(r, "HPFineNGroups") for r in rows]
        gw[s] = {
            "n": len(rows),
            "clustered_frac": round(float(np.mean(clu)), 4),
            "rho_clustered_vs_coverage": spear(clu, cov),
            "rho_clustered_vs_nCpG": spear(clu, ncpg),
            "rho_NGroups_vs_coverage": spear(ng, cov),
            "partial_clustered_vs_coverage_given_nCpG": partial_spear(clu, cov, ncpg),
        }
    out["A_coverage_genomewide_TP"] = gw

    # ---------- B. HCC1395 full (TP) coverage + CN disentangle ----------
    rows = load_sample("HCC1395", ("tp",))
    clu = np.array([is_clustered(r) for r in rows], float)
    cov = np.array([num(r, "NumReads") or np.nan for r in rows], float)
    ncpg = np.array([num(r, "NumCpGs") or np.nan for r in rows], float)
    ng = np.array([num(r, "HPFineNGroups") or np.nan for r in rows], float)
    cvv = np.array([num(r, "CramersV") or np.nan for r in rows], float)
    cns = np.array([cn_state(cn, r.get("Chr"), int(r.get("Pos"))) for r in rows])
    logcov = np.log1p(cov)

    out["B_HCC1395_TP"] = {
        "n": len(rows),
        "clustered_frac": round(float(np.nanmean(clu)), 4),
        "coverage": {
            "rho_clustered_vs_coverage": spear(clu, cov),
            "rho_NGroups_vs_coverage": spear(ng, cov),
            "rho_cramersV_vs_coverage": spear(cvv, cov),
            "clustered_rate_by_coverage_decile": rate_by_bin(cov, clu),
        },
        "nCpG": {
            "rho_clustered_vs_nCpG": spear(clu, ncpg),
            "rho_cramersV_vs_nCpG": spear(cvv, ncpg),
        },
        "coverage_vs_nCpG_collinearity": spear(cov, ncpg),
        "partial": {
            "clustered_vs_coverage_GIVEN_nCpG": partial_spear(clu, cov, ncpg),
            "clustered_vs_nCpG_GIVEN_coverage": partial_spear(clu, ncpg, cov),
        },
    }

    # ---------- C. CN (HCC1395) ----------
    cn_levels = ["neutral", "gain", "loh", "loss"]
    by_cn = {}
    for lv in cn_levels:
        sel = cns == lv
        if sel.sum() == 0:
            continue
        by_cn[lv] = {
            "n": int(sel.sum()),
            "clustered_rate": round(float(np.nanmean(clu[sel])), 4),
            "median_coverage": round(float(np.nanmedian(cov[sel])), 1),
            "median_NGroups": round(float(np.nanmedian(ng[sel])), 2),
            "median_nCpG": round(float(np.nanmedian(ncpg[sel])), 1),
        }
    # chi2: is_clustered x CN-state
    cats = [lv for lv in cn_levels if (cns == lv).sum() > 0]
    table = np.array([[int(((cns == lv) & (clu == 1)).sum()), int(((cns == lv) & (clu == 0)).sum())]
                      for lv in cats])
    chi2, pchi, dof, _ = stats.chi2_contingency(table)
    # gain vs neutral clustered-rate (the dosage-up direction)
    g, n0 = cns == "gain", cns == "neutral"
    out["C_CN_HCC1395_TP"] = {
        "clustered_rate_by_CN": by_cn,
        "chi2_clustered_x_CN": {"chi2": round(float(chi2), 2), "p": float(f"{pchi:.2e}"), "dof": int(dof),
                                "levels": cats},
        "gain_vs_neutral": {"gain_rate": round(float(np.nanmean(clu[g])), 4),
                            "neutral_rate": round(float(np.nanmean(clu[n0])), 4),
                            "gain_median_cov": round(float(np.nanmedian(cov[g])), 1),
                            "neutral_median_cov": round(float(np.nanmedian(cov[n0])), 1),
                            "note": "if gain rate>neutral but gain also has higher coverage -> need within-coverage check"},
    }

    # within-coverage-bin: does CN still raise clustered-rate at FIXED coverage?
    covq = np.nanquantile(cov, [0, .33, .66, 1.0])
    strat = []
    for i in range(3):
        lo, hi = covq[i], covq[i + 1]
        inb = (cov >= lo) & (cov <= hi if i == 2 else cov < hi)
        gg, nn = inb & (cns == "gain"), inb & (cns == "neutral")
        if gg.sum() >= 10 and nn.sum() >= 10:
            strat.append({"cov_bin": f"{lo:.0f}-{hi:.0f}", "n_gain": int(gg.sum()), "n_neutral": int(nn.sum()),
                          "gain_clustered_rate": round(float(np.nanmean(clu[gg])), 4),
                          "neutral_clustered_rate": round(float(np.nanmean(clu[nn])), 4)})
    out["C_CN_HCC1395_TP"]["within_coverage_bin_gain_vs_neutral"] = strat

    # ---------- D. multivariable logistic: is_clustered ~ z(logcov)+z(nCpG)+CN dummies ----------
    def z(a):
        a = a.astype(float); mu, sd = np.nanmean(a), np.nanstd(a)
        return (a - mu) / (sd if sd else 1)
    gain_d = (cns == "gain").astype(float)
    loh_d = (cns == "loh").astype(float)
    loss_d = (cns == "loss").astype(float)
    X = np.c_[z(logcov), z(ncpg), gain_d, loh_d, loss_d]
    names = ["z_logCoverage", "z_nCpG", "CN_gain", "CN_loh", "CN_loss"]
    mfin = np.all(np.isfinite(X), axis=1) & np.isfinite(clu)
    Xf, yf = sm.add_constant(X[mfin]), clu[mfin]
    res = sm.Logit(yf, Xf).fit(disp=0)
    coefs = {}
    for i, nm in enumerate(["const"] + names):
        coefs[nm] = {"coef": round(float(res.params[i]), 4), "p": float(f"{res.pvalues[i]:.2e}"),
                     "OR": round(float(np.exp(res.params[i])), 3)}
    out["D_logistic_is_clustered"] = {
        "n": int(mfin.sum()), "pseudo_R2": round(float(res.prsquared), 4), "coefficients": coefs,
        "interpretation": "standardized coef magnitude + p rank the drivers; CN dummies test CN effect AFTER controlling coverage+nCpG",
    }

    json.dump(out, open(f"{GS}/catalog/clusterability_vs_cov_cn.json", "w"), indent=1, ensure_ascii=False)

    # console verdict
    b = out["B_HCC1395_TP"]
    print("=== HCC1395 TP clusterability vs coverage / CN ===")
    print(f"clustered_frac={b['clustered_frac']}")
    print(f"rho(clustered, coverage)={b['coverage']['rho_clustered_vs_coverage']}")
    print(f"rho(clustered, nCpG)    ={b['nCpG']['rho_clustered_vs_nCpG']}")
    print(f"coverage~nCpG collinearity={b['coverage_vs_nCpG_collinearity']}")
    print(f"PARTIAL clustered~coverage | nCpG = {b['partial']['clustered_vs_coverage_GIVEN_nCpG']}")
    print(f"PARTIAL clustered~nCpG | coverage = {b['partial']['clustered_vs_nCpG_GIVEN_coverage']}")
    print("CN clustered-rate:", {k: v["clustered_rate"] for k, v in out["C_CN_HCC1395_TP"]["clustered_rate_by_CN"].items()})
    print("chi2 clustered x CN:", out["C_CN_HCC1395_TP"]["chi2_clustered_x_CN"])
    print("within-cov-bin gain vs neutral:", out["C_CN_HCC1395_TP"]["within_coverage_bin_gain_vs_neutral"])
    print("logistic coefs:", {k: (v["coef"], v["p"]) for k, v in coefs.items()})
    print("genome-wide rho(clustered,coverage):", {s: gw[s]["rho_clustered_vs_coverage"]["rho"] for s in gw})
    print(f"\n[101] wrote {GS}/catalog/clusterability_vs_cov_cn.json")
    print("DONE")


if __name__ == "__main__":
    main()
