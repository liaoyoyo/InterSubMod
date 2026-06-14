#!/usr/bin/env python3
"""P1a whole-genome HP-AUC distribution — does distance recover germline-HP across the genome?

Reads a SKIP-default summary that now carries HP_AUC_Normal/Tumor/All columns (cpp 185db26).
Env U1_SKIP (default = whole-genome SKIP run). Arg1=OUTDIR.

Interpretation thresholds: >=0.7 strong, 0.6-0.7 moderate, 0.5-0.6 weak, <0.5 anti.
HP_AUC=-1 means undefined (no same+diff HP pairs, e.g. tumor single-HP at somatic site).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

SKIP = os.environ.get("U1_SKIP", "/tmp/wg_hpauc/significance_summary.csv")
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUTDIR, exist_ok=True)


def num(s):
    return pd.to_numeric(s, errors="coerce")


def truthy(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


ds = pd.read_csv(SKIP)
N = len(ds)


def auc_dist(col, mask=None):
    v = num(ds[col])
    if mask is not None:
        v = v[mask]
    valid = v[v >= 0]
    if len(valid) == 0:
        return {"n_valid": 0, "n_undefined": int((v < 0).sum())}
    return {
        "n_valid": int(len(valid)), "n_undefined": int((v < 0).sum()),
        "median": round(float(valid.median()), 3), "mean": round(float(valid.mean()), 3),
        "frac_ge0.7_strong": round(float((valid >= 0.7).mean()), 3),
        "frac_0.6_0.7_mod": round(float(((valid >= 0.6) & (valid < 0.7)).mean()), 3),
        "frac_0.5_0.6_weak": round(float(((valid >= 0.5) & (valid < 0.6)).mean()), 3),
        "frac_lt0.5_anti": round(float((valid < 0.5).mean()), 3),
    }


report = {"input": SKIP, "n_regions": N}
report["HP_AUC_Normal"] = auc_dist("HP_AUC_Normal")
report["HP_AUC_Tumor"] = auc_dist("HP_AUC_Tumor")
report["HP_AUC_All"] = auc_dist("HP_AUC_All")

# stratify: tumor both-HP regions -> can tumor's own HP be separated by distance?
t1, t2 = num(ds["Tumor_HP1"]), num(ds["Tumor_HP2"])
both = (t1 > 0) & (t2 > 0)
report["tumor_both_hp"] = {
    "n_both_hp": int(both.sum()),
    "tumor_auc_in_both": auc_dist("HP_AUC_Tumor", mask=both),
    "normal_auc_in_both": auc_dist("HP_AUC_Normal", mask=both),
}

# sig vs all
sig = truthy(ds["Significant"])
report["sig_subset"] = {
    "n_sig": int(sig.sum()),
    "normal_auc_in_sig": auc_dist("HP_AUC_Normal", mask=sig),
}

# is high normal HP-AUC associated with being called Significant? (sanity: distance works -> sig)
nv = num(ds["HP_AUC_Normal"])
report["normal_auc_vs_significant"] = {
    "median_normal_auc_when_sig": round(float(nv[sig & (nv >= 0)].median()), 3) if (sig & (nv >= 0)).sum() else None,
    "median_normal_auc_when_notsig": round(float(nv[~sig & (nv >= 0)].median()), 3) if (~sig & (nv >= 0)).sum() else None,
}

with open(f"{OUTDIR}/wg_hpauc.json", "w") as f:
    json.dump(report, f, indent=2, ensure_ascii=False)
print(json.dumps(report, indent=2, ensure_ascii=False))
