#!/usr/bin/env python3
"""U1: MAX_DIST vs SKIP nan-distance-strategy — subclone clustering quality comparison.

Apples-to-apples: same input, only --nan-distance-strategy differs.
MAX_DIST = invalid pair (common CpG < C_min) set to 1.0 (pollutes distance).
SKIP     = invalid pair set to NaN -> read filtered out of clustering (fix 517ed90).

Input paths via env (default = chr1 pilot):
  U1_MAXD  significance_summary.csv for MAX_DIST run
  U1_SKIP  significance_summary.csv for SKIP run
Arg1 = OUTDIR.

Outputs to OUTDIR:
  comparison_summary.json   — aggregate dists per strategy + axis breakdown + change counts + confound metrics
  vc_crosstab.tsv           — VerificationClass MAX_DIST(row) x SKIP(col)
  region_full_compare.tsv   — per-region side-by-side
  region_diff_changed.tsv   — only regions whose VC or Significant changed
"""
import json
import os
import sys

import pandas as pd

MAXD = os.environ.get("U1_MAXD", "/tmp/u1_maxdist/significance_summary.csv")
SKIP = os.environ.get("U1_SKIP", "/tmp/u1_skip/significance_summary.csv")
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUTDIR, exist_ok=True)


def truthy(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def dist(series):
    return {str(k): int(v) for k, v in series.value_counts(dropna=False).items()}


def num(s):
    return pd.to_numeric(s, errors="coerce")


def strat_block(d):
    sig = truthy(d["Significant"])
    sub = d[sig]
    n_sig = int(sig.sum())
    return {
        "VerificationClass": dist(d["VerificationClass"]),
        "Significant": n_sig,
        "PassedGating": int(truthy(d["PassedGating"]).sum()),
        "ClusterPermanovaValid": int(truthy(d["ClusterPermanovaValid"]).sum()),
        "CramersV_mean": round(float(num(d["CramersV"]).mean()), 5),
        "CramersV_median": round(float(num(d["CramersV"]).median()), 5),
        "HeuristicScore_mean": round(float(num(d["HeuristicScore"]).mean()), 5),
        "NumReads_median": round(float(num(d["NumReads"]).median()), 1),
        "NumCpGs_median": round(float(num(d["NumCpGs"]).median()), 1),
        # axis breakdown among Significant regions (germline-HP vs subclone discriminant)
        "sig_DominantLabel": dist(sub["DominantLabel"]),
        "sig_hp_fraction": round((sub["DominantLabel"] == "hp").sum() / max(1, n_sig), 3),
        "sig_HPMergedSig": int(truthy(sub["HPMergedSig"]).sum()),
        "sig_AlleleSig": int(truthy(sub["AlleleSig"]).sum()),
        "sig_SampleASM_Sig": int(truthy(sub["SampleASM_Sig"]).sum()),
    }


dm = pd.read_csv(MAXD)
ds = pd.read_csv(SKIP)
key = "RegionID"
common = sorted(set(dm[key]) & set(ds[key]))
dm = dm[dm[key].isin(common)].sort_values(key).reset_index(drop=True)
ds = ds[ds[key].isin(common)].sort_values(key).reset_index(drop=True)
assert (dm[key].values == ds[key].values).all(), "RegionID misalignment"

report = {
    "inputs": {"MAX_DIST": MAXD, "SKIP": SKIP},
    "n_regions_compared": len(common),
    "MAX_DIST": strat_block(dm),
    "SKIP": strat_block(ds),
}

vc_change = dm["VerificationClass"].astype(str) != ds["VerificationClass"].astype(str)
sig_change = dm["Significant"].astype(str) != ds["Significant"].astype(str)
report["changes"] = {
    "VerificationClass_changed": int(vc_change.sum()),
    "Significant_changed": int(sig_change.sum()),
    "any_changed": int((vc_change | sig_change).sum()),
}

# SKIP-new-significant confound metrics (was MAX_DIST not-sig -> SKIP sig)
sig_m = truthy(dm["Significant"])
sig_s = truthy(ds["Significant"])
new = (~sig_m) & sig_s
cvs = num(ds.loc[new, "CramersV"])
pps = num(ds.loc[new, "ClusterPermanovaP"])
report["skip_new_significant"] = {
    "count": int(new.sum()),
    "lost_by_skip": int((sig_m & ~sig_s).sum()),
    "NumReads_median": round(float(num(dm.loc[new, "NumReads"]).median()), 1),
    "CramersV_median": round(float(cvs.median()), 3),
    "CramersV_gt0.1_frac": round(float((cvs > 0.1).mean()), 3),
    "PermanovaP_median": round(float(pps.median()), 4),
    "PermanovaP_lt0.01_count": int((pps < 0.01).sum()),
}

pd.crosstab(dm["VerificationClass"], ds["VerificationClass"], dropna=False).to_csv(f"{OUTDIR}/vc_crosstab.tsv", sep="\t")
cmp = pd.DataFrame({
    "RegionID": dm["RegionID"], "Chr": dm.get("Chr"), "Pos": dm.get("Pos"),
    "NumReads": dm["NumReads"], "NumCpGs": dm["NumCpGs"],
    "VC_maxdist": dm["VerificationClass"], "VC_skip": ds["VerificationClass"],
    "Sig_maxdist": dm["Significant"], "Sig_skip": ds["Significant"],
    "DomLabel_maxdist": dm["DominantLabel"], "DomLabel_skip": ds["DominantLabel"],
    "CramersV_maxdist": dm["CramersV"], "CramersV_skip": ds["CramersV"],
    "PermanovaP_maxdist": dm["ClusterPermanovaP"], "PermanovaP_skip": ds["ClusterPermanovaP"],
})
cmp.to_csv(f"{OUTDIR}/region_full_compare.tsv", sep="\t", index=False, float_format="%.6g")
cmp[vc_change | sig_change].to_csv(f"{OUTDIR}/region_diff_changed.tsv", sep="\t", index=False, float_format="%.6g")

with open(f"{OUTDIR}/comparison_summary.json", "w") as f:
    json.dump(report, f, indent=2, ensure_ascii=False)
print(json.dumps(report, indent=2, ensure_ascii=False))
