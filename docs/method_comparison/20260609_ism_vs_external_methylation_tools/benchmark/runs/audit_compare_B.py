#!/usr/bin/env python3
"""Audit B compare: ISM-style Fisher p (per-CpG) vs DSS beta-binomial p (same counts).
   Over-dispersion signature = DSS p systematically >= Fisher p (Fisher anti-conservative)."""
import csv, json
import numpy as np
B = "/big7_disk/liaoyoyo2001/InterSubMod/docs/method_comparison/20260609_ism_vs_external_methylation_tools/benchmark"
res = {}
for locus in ["BRCA2", "TBC1D16"]:
    fish = {}
    for row in csv.DictReader(open(f"{B}/runs/counts_{locus}.tsv"), delimiter="\t"):
        if row["fisher_p"] not in ("", "None", "nan"):
            fish[int(row["pos"])] = float(row["fisher_p"])
    dss = {}
    try:
        for row in csv.DictReader(open(f"{B}/runs/dss_{locus}.tsv"), delimiter="\t"):
            try:
                pv = float(row["pval"])
                if not np.isnan(pv):
                    dss[int(round(float(row["pos"])))] = pv
            except Exception:
                pass
    except FileNotFoundError:
        res[locus] = {"error": "dss output missing"}; continue
    common = [k for k in fish if k in dss]
    if len(common) < 3:
        res[locus] = {"n_cpg": len(common), "note": "too few matched CpG"}; continue
    fp = np.array([fish[k] for k in common]); dp = np.array([dss[k] for k in common])
    res[locus] = {
        "n_cpg": len(common),
        "median_fisher_p": round(float(np.median(fp)), 5),
        "median_dss_p": round(float(np.median(dp)), 5),
        "frac_fisher_sig_p<0.05": round(float(np.mean(fp < 0.05)), 4),
        "frac_dss_sig_p<0.05": round(float(np.mean(dp < 0.05)), 4),
        "fisher_sig_but_dss_NOT": int(np.sum((fp < 0.05) & (dp >= 0.05))),
        "dss_sig_but_fisher_NOT": int(np.sum((dp < 0.05) & (fp < 0.05) * 0 + (dp < 0.05) & (fp >= 0.05))),
        "pct_DSS_p_>=_Fisher_p": round(float(np.mean(dp >= fp)), 4),
        "median_log10(dss_p/fisher_p)": round(float(np.median(np.log10((dp + 1e-300) / (fp + 1e-300)))), 3),
    }
json.dump(res, open(f"{B}/runs/audit_B.json", "w"), indent=2)
print(json.dumps(res, indent=2))
