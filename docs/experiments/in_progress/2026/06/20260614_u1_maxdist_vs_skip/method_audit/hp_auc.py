#!/usr/bin/env python3
"""HP-AUC sanity check: does the read-read distance recover the germline HP label?

Ground truth = HP tag (germline haplotype from longphase). A good distance makes
"different-HP read pairs" farther than "same-HP read pairs".
  AUC = P(dist(diff-HP pair) > dist(same-HP pair))  -- 1.0 perfect, 0.5 unrelated.

Variants compared (same case regions):
  BERNOULLI_MAXDIST, NHD_MAXDIST  -> metric comparison (no NaN, apples-to-apples)
  BERNOULLI_SKIP                  -> does dropping NaN-pairs make distance recover HP better?
Subsets:
  normal_only -> cleanest germline-HP truth (normal has both HP1+HP2)
  tumor_only  -> at somatic SNV, tumor reads are haplotagged to ~one HP => often single-HP
                 (n_diff=0 => AUC undefined). This asymmetry IS the finding.
  all

CAVEAT: HP = GERMLINE truth. High AUC proves germline-HP capture, NOT somatic-subclone.
SKIP distance has NaN pairs (dropped); those are skipped in the AUC (uses valid pairs only).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

CASES = [
    ("chr1_73791759", "chr1_73786759_73796759"),
    ("chr1_56791672", "chr1_56786672_56796672"),
    ("chr1_58306441", "chr1_58301441_58311441"),
    ("chr1_67579144", "chr1_67574144_67584144"),
]
# (label, ISM output root, metric subdir)
VARIANTS = [
    ("BERNOULLI_MAXDIST", "/tmp/u1_case_max_dist", "BERNOULLI"),
    ("NHD_MAXDIST", "/tmp/u1_case_nhd", "NHD"),
    ("BERNOULLI_SKIP", "/tmp/u1_case_skip", "BERNOULLI"),
]
SUBSETS = {
    "normal_only": lambda hp, t: hp >= 0 and t == 0,
    "tumor_only": lambda hp, t: hp >= 0 and t == 1,
    "all": lambda hp, t: hp >= 0,
}
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUTDIR, exist_ok=True)


def hp_family(h):
    h = str(h)
    if h.startswith("1"):
        return 0
    if h.startswith("2"):
        return 1
    return -1


def auc_diff_gt_same(D, idx, hp):
    same, diff = [], []
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            d = D[i, j]
            if np.isnan(d):
                continue
            (same if hp[i] == hp[j] else diff).append(d)
    if not same or not diff:
        return None, len(same), len(diff)
    same = np.array(same)
    diff = np.array(diff)
    wins = sum(np.sum(d > same) + 0.5 * np.sum(d == same) for d in diff)
    return round(wins / (len(diff) * len(same)), 3), len(same), len(diff)


report = {"cases": {}, "variant_subset_means": {}}
acc = {v[0]: {s: [] for s in SUBSETS} for v in VARIANTS}

for rid, rsub in CASES:
    report["cases"][rid] = {}
    reads = pd.read_csv(f"{VARIANTS[0][1]}/u1_cases/chr1/{rid}/{rsub}/reads/reads.tsv", sep="\t")
    hp = np.array([hp_family(h) for h in reads["hp"]])
    is_tumor = reads["is_tumor"].values.astype(int)
    n_tumor_hp = int(((hp >= 0) & (is_tumor == 1)).sum())
    n_tumor_hp1 = int(((hp == 0) & (is_tumor == 1)).sum())
    n_tumor_hp2 = int(((hp == 1) & (is_tumor == 1)).sum())
    report["cases"][rid]["_tumor_hp_counts"] = {"hp1": n_tumor_hp1, "hp2": n_tumor_hp2}
    for vlabel, root, metric in VARIANTS:
        mpath = f"{root}/u1_cases/chr1/{rid}/{rsub}/distance/{metric}/matrix.csv"
        if not os.path.exists(mpath):
            report["cases"][rid][vlabel] = "MISSING"
            continue
        D = pd.read_csv(mpath, index_col=0).values.astype(float)
        res = {}
        for sname, sfilter in SUBSETS.items():
            idx = [i for i in range(len(hp)) if sfilter(hp[i], is_tumor[i])]
            auc, ns, nd = auc_diff_gt_same(D, idx, hp)
            res[sname] = {"auc": auc, "n_same": ns, "n_diff": nd}
            if auc is not None:
                acc[vlabel][sname].append(auc)
        report["cases"][rid][vlabel] = res

for vlabel in acc:
    report["variant_subset_means"][vlabel] = {
        s: round(float(np.mean(vals)), 3) if vals else None for s, vals in acc[vlabel].items()
    }

with open(f"{OUTDIR}/hp_auc.json", "w") as f:
    json.dump(report, f, indent=2, ensure_ascii=False)
print(json.dumps(report, indent=2, ensure_ascii=False))
