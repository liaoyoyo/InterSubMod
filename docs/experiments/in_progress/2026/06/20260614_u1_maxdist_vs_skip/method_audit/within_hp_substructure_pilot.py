#!/usr/bin/env python3
"""Pilot: within-HP substructure — does a single germline-HP split into >=2 methylation groups? (S4/S6).

For each region, within each HP-family (HP1=1/1-1, HP2=2/2-1) with >=6 reads, test if the per-read
mean β is bimodal (1-D best 2-partition: both groups >=3, mean gap >0.15, within-var reduction large).
Quantify how many regions have within-HP substructure + cross with CN (Coverage_Category) and the
VerificationClass false-negatives. Output JSON only.
"""
import csv
import glob
import math
import os
import sys

import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_wg_abc/filtered_snv_tp/chr1"  # placeholder; auto-detect below
SUMMARY = sys.argv[2] if len(sys.argv) > 2 else "/tmp/dbeta_wg_abc/significance_summary.csv"
MIN_GROUP = 3
GAP = 0.15

# auto-detect the per-region root (whole-genome run writes <vcfname>/<chr>/<chr>_<pos>/)
base = os.path.dirname(SUMMARY)
roots = glob.glob(os.path.join(base, "*", "chr*"))
roots = [r for r in roots if os.path.isdir(r)]


def best_bimodal(vals):
    """1-D best 2-partition. Return True if bimodal (both>=MIN_GROUP, gap>GAP, var reduction>0.4)."""
    v = np.sort(np.asarray(vals, float))
    m = len(v)
    if m < 2 * MIN_GROUP:
        return False
    tot_var = v.var() * m
    best = None
    for s in range(MIN_GROUP, m - MIN_GROUP + 1):
        a, b = v[:s], v[s:]
        wv = a.var() * len(a) + b.var() * len(b)
        if best is None or wv < best[0]:
            best = (wv, b.mean() - a.mean())
    if best is None or tot_var <= 0:
        return False
    wv, gap = best
    return (gap > GAP) and (1 - wv / tot_var > 0.4)


def region_means(inner, hp_set):
    rt = os.path.join(inner, "reads", "reads.tsv")
    mc = os.path.join(inner, "methylation", "methylation.csv")
    if not (os.path.exists(rt) and os.path.exists(mc)):
        return None
    keep = {}
    with open(rt) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["is_tumor"] == "1" and row["hp"] in hp_set:
                keep[row["read_id"]] = True
    means = []
    with open(mc) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            if p[0] not in keep:
                continue
            b = [float(x) for x in p[1:] if x not in ("NA", "")]
            if b:
                means.append(sum(b) / len(b))
    return means


summ = {r["Pos"]: r for r in csv.DictReader(open(SUMMARY))}
n_proc = 0
hp1_bi = hp2_bi = any_bi = 0
bi_cn = bi_normal = 0
bi_in_fn = 0  # within-HP-substructure region that is a classification false-negative


def f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def asig(r):
    return any(r.get(c) == "true" for c in ["GermlineAsmDbeta_Sig", "SubcloneDbeta_HP1_Sig", "SubcloneDbeta_HP2_Sig"])


region_dirs = []
for rt in roots:
    region_dirs += glob.glob(os.path.join(rt, "chr*"))

for d in region_dirs:
    pos = os.path.basename(d).split("_")[-1]
    if pos not in summ:
        continue
    inner = glob.glob(os.path.join(d, "chr*"))
    if not inner:
        continue
    m1 = region_means(inner[0], {"1", "HP1", "1-1"})
    m2 = region_means(inner[0], {"2", "HP2", "2-1"})
    if m1 is None:
        continue
    n_proc += 1
    b1 = best_bimodal(m1) if m1 else False
    b2 = best_bimodal(m2) if m2 else False
    if b1:
        hp1_bi += 1
    if b2:
        hp2_bi += 1
    if b1 or b2:
        any_bi += 1
        r = summ[pos]
        if r["Coverage_Category"] != "Normal":
            bi_cn += 1
        else:
            bi_normal += 1
        if asig(r) and r["VerificationClass"] in ("Weak", "Noise"):
            bi_in_fn += 1

out = {
    "region_root_detected": roots[:2],
    "n_processed": n_proc,
    "hp1_within_substructure": hp1_bi,
    "hp2_within_substructure": hp2_bi,
    "any_hp_within_substructure": any_bi,
    "any_hp_substructure_frac": round(any_bi / n_proc, 4) if n_proc else None,
    "within_substructure_CN_altered": bi_cn,
    "within_substructure_Normal_cov": bi_normal,
    "within_substructure_CN_frac": round(bi_cn / any_bi, 4) if any_bi else None,
    "within_substructure_is_classification_FN": bi_in_fn,
}
import json  # noqa: E402

print(json.dumps(out, indent=2, ensure_ascii=False))
