#!/usr/bin/env python3
"""
102 - EXTEND 101 to FP/FN classes + all 6 samples.

101 found (HCC1395 TP): clusterability is NOT a coverage/CN-dosage artifact — weak coverage
corr (controlled -> negative), LOH suppresses clustering (need >=2 haplotypes), n_CpG is the
power driver. This extends:
  A. coverage corr per (sample x class)              -> 6x3 = 18 cells
  B. CN (SEQC2 truth) x {TP,FP,FN}                   -> HCC1395 only (SEQC2 = HCC1395 truth)
  C. ISM Potential_LOH proxy: clustered-rate LOH-vs-nonLOH per sample (TP)  -> tests if the
     "LOH suppresses clustering" mechanism REPRODUCES across all 6 samples (no SEQC2 for the
     other 5, so use ISM's own LOH call). + HCC1395 validation: ISM-LOH vs SEQC2-LOH agreement.

CAVEAT: FP loci sparse in most samples (H1437=8, HCC1954=29) -> low power flagged. ISM
Potential_LOH is a data-derived proxy, not orthogonal CN truth.

Output: genome_survey_v2/catalog/clusterability_cov_cn_extended.json (deterministic)
"""
import csv, json, os
import numpy as np
from scipy import stats

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
GS = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2"
CNBED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CLASSES = ["tp", "fp", "fn"]
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


def load_cn():
    cn = {}
    for line in open(CNBED):
        p = line.rstrip().split("\t")
        if len(p) >= 4:
            cn.setdefault(p[0], []).append((int(p[1]), int(p[2]), p[3]))
    return cn


def cn_state(cn, chrom, pos):
    for s, e, lab in cn.get(chrom, []):
        if s <= pos <= e:
            return lab
    return "neutral"


def load(s, cl):
    f = f"{EX}/{s}_{cl}/significance_summary.csv"
    return list(csv.DictReader(open(f))) if os.path.exists(f) else []


def spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 30:
        return None
    rho, p = stats.spearmanr(x[m], y[m])
    return {"rho": round(float(rho), 4), "p": float(f"{p:.2e}"), "n": int(m.sum())}


def main():
    cn = load_cn()
    out = {"note": "extend 101 to FP/FN + 6 samples; CN truth (SEQC2) HCC1395-only; ISM Potential_LOH proxy for cross-sample LOH test"}

    # ---------- A. coverage corr per (sample x class) ----------
    A = {}
    for s in SAMPLES:
        for cl in CLASSES:
            rows = load(s, cl)
            if len(rows) < 30:
                A[f"{s}_{cl}"] = {"n": len(rows), "note": "too few loci (low power)"}
                continue
            clu = [is_clustered(r) for r in rows]
            cov = [num(r, "NumReads") for r in rows]
            ncpg = [num(r, "NumCpGs") for r in rows]
            A[f"{s}_{cl}"] = {
                "n": len(rows), "clustered_frac": round(float(np.mean(clu)), 4),
                "rho_clustered_vs_coverage": spear(clu, cov),
                "rho_clustered_vs_nCpG": spear(clu, ncpg),
            }
    out["A_coverage_per_sample_class"] = A

    # ---------- B. CN (SEQC2) x class (HCC1395) ----------
    B = {}
    for cl in CLASSES:
        rows = load("HCC1395", cl)
        if not rows:
            continue
        clu = np.array([is_clustered(r) for r in rows], float)
        cov = np.array([num(r, "NumReads") or np.nan for r in rows], float)
        cns = np.array([cn_state(cn, r.get("Chr"), int(r.get("Pos"))) for r in rows])
        levels = {}
        for lv in ["neutral", "gain", "loh", "loss"]:
            sel = cns == lv
            if sel.sum() >= 5:
                levels[lv] = {"n": int(sel.sum()), "clustered_rate": round(float(np.nanmean(clu[sel])), 4),
                              "median_cov": round(float(np.nanmedian(cov[sel])), 1)}
        # LOH vs non-LOH OR (Fisher) within this class
        loh = cns == "loh"; nonloh = ~loh
        a, b = int(clu[loh].sum()), int((1 - clu[loh]).sum())
        c, d = int(clu[nonloh].sum()), int((1 - clu[nonloh]).sum())
        orr, pf = stats.fisher_exact([[a, b], [c, d]])
        B[cl.upper()] = {"by_CN": levels, "n": len(rows),
                         "LOH_vs_nonLOH": {"OR": round(float(orr), 3), "p": float(f"{pf:.2e}"),
                                           "loh_rate": round(a / max(a + b, 1), 4),
                                           "nonloh_rate": round(c / max(c + d, 1), 4)}}
    out["B_CN_HCC1395_by_class"] = B

    # ---------- C. ISM Potential_LOH proxy across 6 samples (TP) ----------
    C = {}
    for s in SAMPLES:
        rows = load(s, "tp")
        if len(rows) < 30:
            continue
        clu = np.array([is_clustered(r) for r in rows], float)
        loh = np.array([tb(r, "Potential_LOH") for r in rows])
        a, b = int(clu[loh].sum()), int((1 - clu[loh]).sum())
        c, d = int(clu[~loh].sum()), int((1 - clu[~loh]).sum())
        if (a + b) < 10 or (c + d) < 10:
            C[s] = {"note": "insufficient LOH/non-LOH split"}; continue
        orr, pf = stats.fisher_exact([[a, b], [c, d]])
        C[s] = {"n": len(rows), "loh_n": a + b, "nonloh_n": c + d,
                "loh_clustered_rate": round(a / (a + b), 4), "nonloh_clustered_rate": round(c / (c + d), 4),
                "LOH_vs_nonLOH_OR": round(float(orr), 3), "p": float(f"{pf:.2e}")}
    out["C_ISM_LOH_proxy_per_sample_TP"] = C
    ors = [C[s]["LOH_vs_nonLOH_OR"] for s in C if "LOH_vs_nonLOH_OR" in C[s]]
    out["C_summary"] = {"n_samples": len(ors), "n_with_OR_lt_1_LOH_suppresses": sum(1 for o in ors if o < 1),
                        "OR_range": [round(min(ors), 3), round(max(ors), 3)] if ors else None,
                        "verdict": "LOH suppresses clusterability reproduces across samples if all OR<1"}

    # ---------- D. HCC1395 ISM-LOH vs SEQC2-LOH agreement (proxy validation) ----------
    rows = load("HCC1395", "tp")
    ism = np.array([tb(r, "Potential_LOH") for r in rows])
    seq = np.array([cn_state(cn, r.get("Chr"), int(r.get("Pos"))) == "loh" for r in rows])
    tp_ = int((ism & seq).sum()); fp_ = int((ism & ~seq).sum())
    fn_ = int((~ism & seq).sum()); tn_ = int((~ism & ~seq).sum())
    out["D_ISM_LOH_vs_SEQC2_LOH_agreement_HCC1395"] = {
        "ism_loh_n": int(ism.sum()), "seqc2_loh_n": int(seq.sum()),
        "agree_both": tp_, "ism_only": fp_, "seqc2_only": fn_, "agree_neither": tn_,
        "precision_ism_vs_seqc2": round(tp_ / max(tp_ + fp_, 1), 3),
        "recall_ism_vs_seqc2": round(tp_ / max(tp_ + fn_, 1), 3),
        "note": "ISM Potential_LOH is a coverage/HP-ratio proxy, partial overlap with SEQC2 truth expected",
    }

    json.dump(out, open(f"{GS}/catalog/clusterability_cov_cn_extended.json", "w"), indent=1, ensure_ascii=False)

    print("=== A. coverage corr per sample x class (clustered_frac | rho_cov) ===")
    for k, v in A.items():
        if "clustered_frac" in v:
            r = v["rho_clustered_vs_coverage"]
            print(f"  {k:14} n={v['n']:6} frac={v['clustered_frac']:.3f} rho_cov={r['rho'] if r else 'NA'}")
        else:
            print(f"  {k:14} n={v['n']:6} ({v['note']})")
    print("\n=== B. CN x class (HCC1395) LOH-vs-nonLOH OR ===")
    for cl, v in B.items():
        print(f"  {cl}: rates {{lv: by_CN[lv]['clustered_rate']}} ->",
              {lv: v["by_CN"][lv]["clustered_rate"] for lv in v["by_CN"]},
              f"| LOH OR={v['LOH_vs_nonLOH']['OR']} p={v['LOH_vs_nonLOH']['p']}")
    print("\n=== C. ISM-LOH proxy across 6 samples (TP): LOH suppresses? ===")
    for s, v in C.items():
        if "LOH_vs_nonLOH_OR" in v:
            print(f"  {s:9} loh_rate={v['loh_clustered_rate']:.3f} nonloh_rate={v['nonloh_clustered_rate']:.3f} OR={v['LOH_vs_nonLOH_OR']} p={v['p']}")
    print("  summary:", out["C_summary"])
    print("\n=== D. ISM-LOH vs SEQC2-LOH agreement:", out["D_ISM_LOH_vs_SEQC2_LOH_agreement_HCC1395"]["precision_ism_vs_seqc2"],
          "prec /", out["D_ISM_LOH_vs_SEQC2_LOH_agreement_HCC1395"]["recall_ism_vs_seqc2"], "recall")
    print(f"\n[102] wrote {GS}/catalog/clusterability_cov_cn_extended.json\nDONE")


if __name__ == "__main__":
    main()
