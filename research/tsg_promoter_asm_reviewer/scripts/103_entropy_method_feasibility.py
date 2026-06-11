#!/usr/bin/env python3
"""
103 - Is methylation ENTROPY (epipolymorphism/NME) an effective method? Re-test O11 on
current 6-sample data + test the user's DIFFERENTIAL framing (ALT-vs-REF / HP-imbalance).

O11 (2026-03-31, concluded NEGATIVE): ABSOLUTE epipolymorphism AUC 0.845 was 100% n_reads
confound (Epipoly~n_reads r=0.79; residualized -> 0.51-0.58). Externally confirmed (Scherer20:
entropy coverage-dependent >10x).

User idea = DIFFERENTIAL entropy (ALT vs REF / HP1 vs HP2 imbalance). A within-locus DIFFERENCE
may partially cancel the shared-coverage confound -> worth checking if it escapes O11. ISM already
emits Epipoly_HP1/HP2/Delta + NME_HP1/HP2 + Entropy_Imbalance (HP-axis; NOT allele ALT/REF split).

Tests (HCC1395 TP vs FP):
  (1) coverage confound: spearman(feature, NumReads) for ABSOLUTE vs DIFFERENTIAL
  (2) discrimination: AUC(feature, TP-vs-FP) RAW + coverage/nCpG-residualized
  (3) redundancy vs clustering: AUC of is_clustered for reference + corr(entropy, is_clustered)

Output: genome_survey_v2/catalog/entropy_method_feasibility.json
"""
import csv, json, os
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LinearRegression

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
GS = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2"
LATENT_MINN = 5


def num(r, k):
    try:
        v = float(r.get(k, "")); return None if v != v else v
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


def load(s, cl):
    f = f"{EX}/{s}_{cl}/significance_summary.csv"
    return list(csv.DictReader(open(f))) if os.path.exists(f) else []


def auc_safe(y, x):
    m = np.isfinite(x)
    if m.sum() < 30 or len(set(y[m])) < 2:
        return None
    a = roc_auc_score(y[m], x[m])
    return round(float(max(a, 1 - a)), 4)  # report |deviation from 0.5| direction-agnostic


def resid(x, Z):
    """residualize x on columns Z (handle nan by mean-impute within finite mask)."""
    m = np.isfinite(x) & np.all(np.isfinite(Z), axis=1)
    out = np.full_like(x, np.nan, dtype=float)
    if m.sum() < 30:
        return out, m
    lr = LinearRegression().fit(Z[m], x[m])
    out[m] = x[m] - lr.predict(Z[m])
    return out, m


def spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 30:
        return None
    rho, p = stats.spearmanr(x[m], y[m])
    return {"rho": round(float(rho), 4), "p": float(f"{p:.2e}"), "n": int(m.sum())}


def main():
    tp = load("HCC1395", "tp"); fp = load("HCC1395", "fp")
    rows = tp + fp
    y = np.array([1] * len(tp) + [0] * len(fp))
    cov = np.array([num(r, "NumReads") or np.nan for r in rows], float)
    ncpg = np.array([num(r, "NumCpGs") or np.nan for r in rows], float)
    clu = np.array([is_clustered(r) for r in rows], float)
    logcov = np.log1p(cov)
    Z = np.c_[logcov, ncpg]

    ABS = {"Epipoly_HP1": "absolute epipoly (O11-style level)",
           "NME_HP1": "absolute NME level"}
    DIFF = {"Epipoly_Delta": "HP1-HP2 epipoly imbalance (user differential)",
            "Entropy_Imbalance": "NME HP imbalance (user differential)"}

    res = {"note": "re-test O11 + user differential-entropy framing; HCC1395 TP+FP",
           "O11_recall": "absolute epipoly AUC 0.845 -> residualized 0.53 (n_reads confound, concluded NEGATIVE 2026-03-31)",
           "n_tp": len(tp), "n_fp": len(fp)}

    def block(feats):
        out = {}
        for k, desc in feats.items():
            x = np.array([abs(num(r, k)) if num(r, k) is not None else np.nan for r in rows], float)
            rcov = spear(x, cov)
            raw = auc_safe(y, x)
            xr, m = resid(x, Z)
            resid_auc = auc_safe(y[m], xr[m]) if m.sum() > 30 else None
            out[k] = {"desc": desc,
                      "rho_vs_coverage": rcov,
                      "AUC_TPvsFP_raw": raw,
                      "AUC_TPvsFP_resid_cov_nCpG": resid_auc,
                      "rho_vs_is_clustered": spear(x, clu)}
        return out

    res["ABSOLUTE_entropy_O11_style"] = block(ABS)
    res["DIFFERENTIAL_entropy_user_idea"] = block(DIFF)
    # reference: clustering itself
    res["reference_is_clustered_AUC"] = auc_safe(y, clu)
    # ALT vs REF note
    res["ALT_vs_REF_entropy_status"] = ("NOT computed by ISM — entropy columns are HP-axis (HP1 vs HP2), "
                                        "not allele ALT/REF split. Epipoly_Delta/Entropy_Imbalance (HP1 vs HP2) "
                                        "is the closest existing proxy. True ALT-vs-REF entropy would need a new "
                                        "allele-split extraction AND a germline-het null control (ALLELE-axis baseline confound).")

    json.dump(res, open(f"{GS}/catalog/entropy_method_feasibility.json", "w"), indent=1, ensure_ascii=False)

    print("=== O11 recall: absolute epipoly = n_reads confound (AUC 0.845->0.53) ===")
    print(f"TP={len(tp)} FP={len(fp)}  reference is_clustered AUC(TPvsFP)={res['reference_is_clustered_AUC']}\n")
    for grp in ("ABSOLUTE_entropy_O11_style", "DIFFERENTIAL_entropy_user_idea"):
        print(f"--- {grp} ---")
        for k, v in res[grp].items():
            print(f"  {k:18} rho_cov={v['rho_vs_coverage']['rho'] if v['rho_vs_coverage'] else 'NA'} "
                  f"AUC_raw={v['AUC_TPvsFP_raw']} AUC_resid={v['AUC_TPvsFP_resid_cov_nCpG']} "
                  f"rho_clustered={v['rho_vs_is_clustered']['rho'] if v['rho_vs_is_clustered'] else 'NA'}")
    print(f"\n[103] wrote {GS}/catalog/entropy_method_feasibility.json\nDONE")


if __name__ == "__main__":
    main()
