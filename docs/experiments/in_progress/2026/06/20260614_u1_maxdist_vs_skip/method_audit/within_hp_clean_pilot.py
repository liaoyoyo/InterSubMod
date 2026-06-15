#!/usr/bin/env python3
"""Refined within-HP substructure: separate CLEAN multi-group (valid) from tiny-outlier/chaotic (noise).

For each region's HP-family (>=8 reads), 2-means on per-read mean β, then require for a CLEAN multi-group:
  - silhouette >= 0.5 (clean separation, not chaotic)
  - min group >= max(3, 20% of n)  (balanced, not a tiny outlier group)
  - mean gap > 0.15
Reports the loose count (any bimodal) vs clean count, and a silhouette-bin breakdown, cross with NME
(summary) to characterize chaotic. Converges the 40.7% upper bound to a credible "clear multi-group" value.
"""
import csv
import glob
import math
import os
import sys

import numpy as np

SUMMARY = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_wg_abc/significance_summary.csv"
base = os.path.dirname(SUMMARY)
roots = [r for r in glob.glob(os.path.join(base, "*", "chr*")) if os.path.isdir(r)]


def two_means_1d(v):
    """1-D 2-means; return (labels, centers, silhouette, sizes, gap)."""
    v = np.asarray(v, float)
    n = len(v)
    # init by median split
    c = np.array([v.min(), v.max()])
    for _ in range(20):
        d = np.abs(v[:, None] - c[None, :])
        lab = d.argmin(1)
        if lab.min() == lab.max():
            return None
        nc = np.array([v[lab == 0].mean(), v[lab == 1].mean()])
        if np.allclose(nc, c):
            c = nc
            break
        c = nc
    # silhouette (1-D)
    a = np.abs(v - c[lab])
    b = np.abs(v - c[1 - lab])
    s = np.where(np.maximum(a, b) > 0, (b - a) / np.maximum(a, b), 0).mean()
    sizes = (int((lab == 0).sum()), int((lab == 1).sum()))
    return lab, c, float(s), sizes, abs(c[1] - c[0])


def hp_means(inner, hp_set):
    rt, mc = os.path.join(inner, "reads", "reads.tsv"), os.path.join(inner, "methylation", "methylation.csv")
    if not (os.path.exists(rt) and os.path.exists(mc)):
        return None
    keep = {row["read_id"] for row in csv.DictReader(open(rt), delimiter="\t")
            if row["is_tumor"] == "1" and row["hp"] in hp_set}
    means = []
    with open(mc) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            if p[0] in keep:
                b = [float(x) for x in p[1:] if x not in ("NA", "")]
                if b:
                    means.append(sum(b) / len(b))
    return means


def clean_multigroup(means):
    if not means or len(means) < 8:
        return False, 0.0
    r = two_means_1d(means)
    if r is None:
        return False, 0.0
    _, _, sil, sizes, gap = r
    n = len(means)
    balanced = min(sizes) >= max(3, int(0.2 * n))
    return (sil >= 0.5 and balanced and gap > 0.15), sil


summ = {r["Pos"]: r for r in csv.DictReader(open(SUMMARY))}
region_dirs = []
for rt in roots:
    region_dirs += glob.glob(os.path.join(rt, "chr*"))

n_proc = clean = loose = 0
sil_bins = {f"{x/10:.1f}-{(x+2)/10:.1f}": 0 for x in range(0, 10, 2)}
clean_cn = clean_fn = 0


def f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def asig(r):
    return any(r.get(c) == "true" for c in ["GermlineAsmDbeta_Sig", "SubcloneDbeta_HP1_Sig", "SubcloneDbeta_HP2_Sig"])


for d in region_dirs:
    pos = os.path.basename(d).split("_")[-1]
    if pos not in summ:
        continue
    inner = glob.glob(os.path.join(d, "chr*"))
    if not inner:
        continue
    m1 = hp_means(inner[0], {"1", "HP1", "1-1"})
    if m1 is None:
        continue
    n_proc += 1
    m2 = hp_means(inner[0], {"2", "HP2", "2-1"})
    best_sil = 0.0
    is_clean = False
    for m in (m1, m2):
        if m and len(m) >= 2 * 3:
            c, sil = clean_multigroup(m)
            best_sil = max(best_sil, sil)
            if c:
                is_clean = True
    # loose = any bimodal-ish (sil>0.3)
    if best_sil > 0.3:
        loose += 1
    if is_clean:
        clean += 1
        r = summ[pos]
        if r["Coverage_Category"] != "Normal":
            clean_cn += 1
        if asig(r) and r["VerificationClass"] in ("Weak", "Noise"):
            clean_fn += 1
    for lo in range(0, 10, 2):
        if lo / 10 <= best_sil < (lo + 2) / 10:
            sil_bins[f"{lo/10:.1f}-{(lo+2)/10:.1f}"] += 1
            break

out = {
    "n_processed": n_proc,
    "clean_multigroup (sil>=0.5 + balanced + gap>0.15)": clean,
    "clean_frac": round(clean / n_proc, 4) if n_proc else None,
    "loose_bimodal (sil>0.3)": loose,
    "loose_frac": round(loose / n_proc, 4) if n_proc else None,
    "best_silhouette_distribution": sil_bins,
    "clean_multigroup_CN_altered": clean_cn,
    "clean_multigroup_is_classification_FN": clean_fn,
}
import json  # noqa: E402

print(json.dumps(out, indent=2, ensure_ascii=False))
