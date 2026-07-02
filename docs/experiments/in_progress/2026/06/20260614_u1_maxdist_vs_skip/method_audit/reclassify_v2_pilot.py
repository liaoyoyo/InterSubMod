#!/usr/bin/env python3
"""Reclassify v2 pilot — full evidence-based VerificationClass (user philosophy: retain if ANY structure).

Per region combines: within-HP clean multi-group (per-region methylation) + summary evidence
(Δβ / HP-AUC / orig cluster-label match) + SEQC2 CN/LOH overlap + noise sub-typing (Epipoly / PairwiseMeanDist).

VALID (any → keep) -> subtyped: Strong / LabelShift / CN-Substructure / MultiGroupNoLabel / LOH-Structure / StructureNoLabel
NOISE (none) -> Noise_Uniform / Noise_Chaotic / Noise_Uncorrelated
Reports new distribution, rescue of orig Weak/Noise, orig-Strong change (FP check), noise sub-type sizes.
"""
import bisect
import csv
import glob
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

SUMMARY = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_wg_abc/significance_summary.csv"
SEQC2_BED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
base = os.path.dirname(SUMMARY)
roots = [r for r in glob.glob(os.path.join(base, "*", "chr*")) if os.path.isdir(r)]

# SEQC2 intervals
seqc2 = defaultdict(list)
for line in open(SEQC2_BED):
    c, s, e, t = line.split("\t")[:4]
    seqc2[c].append((int(s), int(e), t.strip()))
for c in seqc2:
    seqc2[c].sort()
seqc2_starts = {c: [x[0] for x in iv] for c, iv in seqc2.items()}


def seqc2_type(chrom, pos):
    iv = seqc2.get(chrom, [])
    if not iv:
        return "none"
    i = bisect.bisect_right(seqc2_starts[chrom], pos) - 1
    for j in range(max(0, i - 1), min(len(iv), i + 2)):
        if iv[j][0] <= pos <= iv[j][1]:
            return iv[j][2]
    return "none"


def f(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


def two_means_clean(v):
    v = np.sort(np.asarray(v, float))
    n = len(v)
    if n < 8:
        return False
    c = np.array([v[0], v[-1]])
    for _ in range(20):
        lab = (np.abs(v[:, None] - c[None, :])).argmin(1)
        if lab.min() == lab.max():
            return False
        nc = np.array([v[lab == 0].mean(), v[lab == 1].mean()])
        if np.allclose(nc, c):
            break
        c = nc
    sizes = (int((lab == 0).sum()), int((lab == 1).sum()))
    gap = abs(c[1] - c[0])
    return min(sizes) >= max(3, int(0.2 * n)) and gap > 0.15


def hp_means(inner, hp_set):
    rt, mc = os.path.join(inner, "reads", "reads.tsv"), os.path.join(inner, "methylation", "methylation.csv")
    if not (os.path.exists(rt) and os.path.exists(mc)):
        return None
    keep = {row["read_id"] for row in csv.DictReader(open(rt), delimiter="\t")
            if row["is_tumor"] == "1" and row["hp"] in hp_set}
    means = []
    with open(mc) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split(",")
            if p[0] in keep:
                b = [float(x) for x in p[1:] if x not in ("NA", "")]
                if b:
                    means.append(sum(b) / len(b))
    return means


summ = {r["Pos"]: r for r in csv.DictReader(open(SUMMARY))}
region_dirs = []
for rt in roots:
    region_dirs += glob.glob(os.path.join(rt, "chr*"))

orig_vc, new_vc = [], []
n_proc = 0
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
    r = summ[pos]
    m2 = hp_means(inner[0], {"2", "HP2", "2-1"})
    within_clean = (m1 and two_means_clean(m1)) or (m2 and two_means_clean(m2))
    dbeta = any(r.get(c) == "true" for c in ["GermlineAsmDbeta_Sig", "SubcloneDbeta_HP1_Sig", "SubcloneDbeta_HP2_Sig",
                                             "AltSubcloneDbeta_HP1_Sig", "AltSubcloneDbeta_HP2_Sig", "SomaticResidualDbeta_Sig"])
    aucv = f(r.get("HP_AUC_All"))
    hp_auc = aucv is not None and aucv >= 0.7
    ovc = r["VerificationClass"]
    orig_match = ovc in ("Strong", "Subclone")  # cluster matched a label via GlobalTest
    s2 = seqc2_type(r["Chr"], int(f(r["Pos"])))
    loh_struct = s2 == "loh" and (within_clean or hp_auc)

    valid = dbeta or hp_auc or orig_match or within_clean or loh_struct
    if valid:
        if orig_match:
            v = "Strong"
        elif within_clean and s2 == "gain":
            v = "CN-Substructure"
        elif within_clean:
            v = "MultiGroupNoLabel"
        elif loh_struct:
            v = "LOH-Structure"
        elif dbeta:
            v = "LabelShift"
        else:
            v = "StructureNoLabel"
    else:
        ep = f(r.get("Epipoly_HP1"))
        pw = f(r.get("PairwiseMeanDist"))
        if (pw is not None and pw < 0.15) or (ep is not None and ep < 0.2):
            v = "Noise_Uniform"
        elif (ep is not None and ep > 0.5) or (pw is not None and pw > 0.35):
            v = "Noise_Chaotic"
        else:
            v = "Noise_Uncorrelated"
    orig_vc.append(ovc)
    new_vc.append(v)

n = len(orig_vc)
valid_classes = ("Strong", "LabelShift", "CN-Substructure", "MultiGroupNoLabel", "LOH-Structure", "StructureNoLabel")
rescued = sum(1 for i in range(n) if orig_vc[i] in ("Weak", "Noise") and new_vc[i] in valid_classes)
weaknoise = sum(1 for x in orig_vc if x in ("Weak", "Noise"))
strong_changed = sum(1 for i in range(n) if orig_vc[i] == "Strong" and new_vc[i] != "Strong")
out = {
    "n_processed": n,
    "orig_dist": dict(Counter(orig_vc)),
    "new_dist": dict(Counter(new_vc)),
    "rescued_weaknoise": rescued,
    "weaknoise_total": weaknoise,
    "rescue_rate": round(rescued / weaknoise, 4) if weaknoise else None,
    "orig_Strong_changed (FP check)": strong_changed,
    "noise_subtypes": {k: v for k, v in Counter(new_vc).items() if k.startswith("Noise")},
}
import json  # noqa: E402

print(json.dumps(out, indent=2, ensure_ascii=False))
