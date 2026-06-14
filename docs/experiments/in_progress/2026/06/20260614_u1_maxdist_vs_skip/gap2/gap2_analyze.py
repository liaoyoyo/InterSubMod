#!/usr/bin/env python3
"""gap2: somatic-subclone vs germline-HP discrimination from a paired SKIP summary.

Reproducible — re-run on ANY significance_summary.csv to observe/verify the same metrics.
This is the script behind gap2/README.md numbers; keeps the analysis versioned (was inline).

Env:
  U1_SKIP   path to SKIP-strategy significance_summary.csv (default = whole-genome wg run)
Arg1 = OUTDIR (default current dir).

Method caveats this script documents (see gap2/README.md §方法質疑):
  - HP test is distance-permutation based (NOT methylation-mean); signed-delta below is a
    SEPARATE mean-axis proxy and does NOT align with ISM's distance-based HP significance.
  - ISM's hp_residual_sig = tumor_hp.significant (RegionProcessor.cpp:978), i.e. it does NOT
    test the residual; so HP_Residual_Sig is reported but flagged unreliable here.
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
sig = truthy(ds["Significant"])
sub = ds[sig].copy()
n = int(sig.sum())

r = {"input": SKIP, "n_sig": n}

# (1) normal-control signal prevalence among significant regions
r["sig_signal_prevalence"] = {
    col: {"count": int(truthy(sub[col]).sum()), "frac": round(float(truthy(sub[col]).mean()), 3)}
    for col in ["HPMergedSig", "HP_Residual_Sig", "Normal_HP_Valid", "Tumor_HP_Valid",
                "SampleASM_Sig", "AlleleSig", "Potential_LOH"]
    if col in sub.columns
}

# (2) tumor x normal HP validity cross-tab
tv = truthy(sub["Tumor_HP_Valid"])
nv = truthy(sub["Normal_HP_Valid"])
r["tumor_x_normal_valid"] = {
    "both": int((tv & nv).sum()), "tumor_only": int((tv & ~nv).sum()),
    "normal_only": int((~tv & nv).sum()), "neither": int((~tv & ~nv).sum()),
}

# (3) signed-delta residual (mean-axis proxy; see caveat — not the distance HP test)
res = num(sub["HP_Signed_Residual"])
td = num(sub["Tumor_HP_Signed_Delta"])
nd = num(sub["Normal_HP_Signed_Delta"])
ar = res.abs().dropna()
r["signed_residual"] = {
    "available": int(res.notna().sum()), "nan_single_hp": int(res.isna().sum()),
    "abs_median": round(float(ar.median()), 4), "abs_mean": round(float(ar.mean()), 4),
    "gt0.05": int((ar > 0.05).sum()), "gt0.1": int((ar > 0.1).sum()), "gt0.2": int((ar > 0.2).sum()),
}
both = td.notna() & nd.notna()
r["tumor_vs_normal_signed_corr"] = round(float(np.corrcoef(td[both], nd[both])[0, 1]), 3) if both.sum() > 2 else None
r["abs_tumor_signed_median"] = round(float(td.abs().median()), 4)
r["abs_normal_signed_median"] = round(float(nd.abs().median()), 4)

# (4) discrimination estimates (caveated — mean-axis proxy, not the corrected distance test)
r["somatic_candidate_upper_bound"] = {
    "rule": "tumorValid & |signed_residual|>0.1 (mean-axis proxy)",
    "count": int((tv & (res.abs() > 0.1)).sum()),
    "frac": round(float((tv & (res.abs() > 0.1)).sum()) / n, 3),
}
r["germline_dominant"] = {
    "rule": "normalValid & |signed_residual|<=0.05",
    "count": int((nv & (res.abs() <= 0.05)).sum()),
}

with open(f"{OUTDIR}/gap2_subclone_discrimination.json", "w") as f:
    json.dump(r, f, indent=2, ensure_ascii=False)
print(json.dumps(r, indent=2, ensure_ascii=False))
