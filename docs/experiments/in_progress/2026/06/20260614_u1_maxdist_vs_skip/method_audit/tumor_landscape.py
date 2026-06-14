#!/usr/bin/env python3
"""Whole-genome tumor landscape: single-HP prevalence, LOH, CNV — from a SKIP summary.

Reproducible (env U1_SKIP, default = whole-genome wg run). Behind the tumor-status numbers
in plans/20260614_ism_method_validation_plan_01.md.

Key correction this script documents: tumor single-HP is NOT universal (whole-genome ~32%,
vs 100% in the 4 CramersV=1.0 cherry-picked case regions = selection bias). 67.7% of regions
have tumor both-HP -> room for somatic-subclone signal.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

SKIP = os.environ.get("U1_SKIP", "/tmp/u1_wg_skip/significance_summary.csv")
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUTDIR, exist_ok=True)


def truthy(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def num(s):
    return pd.to_numeric(s, errors="coerce")


ds = pd.read_csv(SKIP)
N = len(ds)
t1, t2 = num(ds["Tumor_HP1"]), num(ds["Tumor_HP2"])
has_t = (t1 + t2) > 0
both = (t1 > 0) & (t2 > 0)
single = has_t & ~both
cm = num(ds["Coverage_Multiple"])

r = {
    "input": SKIP, "n_regions": N,
    "tumor_hp": {
        "has_tumor_hp": int(has_t.sum()),
        "single_hp": int(single.sum()), "single_hp_frac_of_has": round(float(single.sum() / max(1, has_t.sum())), 3),
        "both_hp": int(both.sum()), "both_hp_frac_of_has": round(float(both.sum() / max(1, has_t.sum())), 3),
        "hp1_only": int((has_t & (t1 > 0) & (t2 == 0)).sum()), "hp2_only": int((has_t & (t2 > 0) & (t1 == 0)).sum()),
    },
    "loh": {
        "potential_loh": int(truthy(ds["Potential_LOH"]).sum()),
        "potential_loh_frac": round(float(truthy(ds["Potential_LOH"]).mean()), 3),
        "loh_subtype": {str(k): int(v) for k, v in ds["LOH_Subtype"].value_counts(dropna=False).items()},
        "hp_ratio_median": round(float(num(ds["HP_Ratio"]).median()), 3),
    },
    "cnv": {
        "coverage_category": {str(k): int(v) for k, v in ds["Coverage_Category"].value_counts(dropna=False).items()},
        "coverage_multiple_median": round(float(cm.median()), 3),
        "amplification_gt1.5_frac": round(float((cm > 1.5).mean()), 3),
        "deletion_lt0.7_frac": round(float((cm < 0.7).mean()), 3),
    },
}
with open(f"{OUTDIR}/tumor_landscape.json", "w") as f:
    json.dump(r, f, indent=2, ensure_ascii=False)
print(json.dumps(r, indent=2, ensure_ascii=False))
