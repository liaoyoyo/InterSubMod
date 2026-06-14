#!/usr/bin/env python3
"""Discrimination (allele read-grouping) vs window + RELIABILITY test of the correlation plateau.

Goal: show whether the methylation read-clustering DISCRIMINATION (between-allele vs within-allele
methylation difference) is preserved / improved at +-5000 vs +-1000, AND test whether the
distance-independent correlation plateau (from comethyl_window_analysis) is GROUP(allele)-driven.

Reliability built-in:
  - per-region (not pooled) discrimination: mean within/between abs-diff per region, then aggregate
    => large regions do not dominate.
  - require BOTH allele groups present (>=2 reads each) and >=10 valid pairs per side.
  - read-pair distance = mean |methyl_prob_i - methyl_prob_j| over >=3 common CpG (ISM C_min=3),
    interpretable, no std requirement.
  - KEY validity test: CpG-CpG correlation vs distance computed 3 ways:
      pooled (all reads) | within-REF only | within-ALT only.
    If the 1-5kb plateau is allele/haplotype-driven, within-group plateau should DROP toward 0.
  - scope caveat: regions are somatic-SNV-anchored (V6_on_fp FP-candidates + V3F_on_tp TP);
    ALLELE axis (ALT/REF) can be germline-allelic-confounded (orthogonal to window question).

Output: research/flagship_chr2_18086020_20260612/discrimination_window_analysis.json
"""
import glob, random, csv, json, os, re, itertools
import numpy as np

BASE = "research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance"
PATTERNS = [
    f"{BASE}/V6_on_fp/filtered_snv_fp/chr*/chr*/chr*_*",
    f"{BASE}/V3F_on_tp/filtered_snv_tp/chr*/chr*/chr*_*",
]
OUT = "research/flagship_chr2_18086020_20260612/discrimination_window_analysis.json"
N_REGIONS = 600
random.seed(7)
ANCHOR_RE = re.compile(r"/chr[0-9XY]+_(\d+)/chr[0-9XY]+_\d+_\d+$")
WINDOWS = (1000, 5000)
DIST_BINS = [(0, 200), (200, 600), (600, 1000), (1000, 2000), (2000, 3500), (3500, 5000)]


def load_region(d):
    mp, rp = f"{d}/methylation/methylation.csv", f"{d}/reads/reads.tsv"
    m = ANCHOR_RE.search(d)
    if not (m and os.path.exists(mp) and os.path.exists(rp)):
        return None
    anchor = int(m.group(1))
    with open(mp) as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        positions = np.array([int(x) for x in header[1:]])
        rows = [[np.nan if x == "NA" else float(x) for x in r[1:]] for r in rdr]
    if not rows or positions.size == 0:
        return None
    M = np.array(rows, dtype=float)
    alleles = []
    with open(rp) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            alleles.append(r.get("alt_support", "").strip())
    alleles = np.array(alleles)
    if M.ndim != 2 or alleles.shape[0] != M.shape[0] or M.shape[1] != positions.size:
        return None
    return anchor, positions, M, alleles


disc = {W: [] for W in WINDOWS}           # per-region (between_mean - within_mean)
disc_pos = {W: 0 for W in WINDOWS}        # regions with between > within
within_all = {W: [] for W in WINDOWS}     # pooled within-allele dists (for global mean)
between_all = {W: [] for W in WINDOWS}
corr = {g: {b: [] for b in DIST_BINS} for g in ["pooled", "within_ref", "within_alt"]}
nreg = nreads = 0

dirs = []
for pat in PATTERNS:
    dirs += glob.glob(pat)
random.shuffle(dirs)
dirs = dirs[:N_REGIONS]

for d in dirs:
    r = load_region(d)
    if r is None:
        continue
    anchor, pos, M, alleles = r
    isref, isalt = alleles == "REF", alleles == "ALT"
    if isref.sum() < 2 or isalt.sum() < 2:
        continue
    nreg += 1
    nreads += M.shape[0]
    notna = ~np.isnan(M)
    nr = M.shape[0]
    allpairs = list(itertools.combinations(range(nr), 2))
    random.shuffle(allpairs)
    allpairs = allpairs[:400]
    for W in WINDOWS:
        inw = np.abs(pos - anchor) <= W
        Mi, nn = M[:, inw], notna[:, inw]
        wi, be = [], []
        for i, j in allpairs:
            if alleles[i] not in ("REF", "ALT") or alleles[j] not in ("REF", "ALT"):
                continue
            mask = nn[i] & nn[j]
            if int(mask.sum()) < 3:
                continue
            dd = float(np.mean(np.abs(Mi[i][mask] - Mi[j][mask])))
            (wi if alleles[i] == alleles[j] else be).append(dd)
        if len(wi) >= 10 and len(be) >= 10:
            wm, bm = float(np.mean(wi)), float(np.mean(be))
            disc[W].append(bm - wm)
            disc_pos[W] += (bm > wm)
            within_all[W] += wi
            between_all[W] += be
    # correlation reliability test (pooled vs within-group), full +-5000
    ncol = M.shape[1]
    cps = list(itertools.combinations(range(ncol), 2))
    random.shuffle(cps)
    for a, b in cps[:700]:
        dd = abs(int(pos[a]) - int(pos[b]))
        binb = next(((lo, hi) for lo, hi in DIST_BINS if lo <= dd < hi), None)
        if binb is None:
            continue
        va, vb = M[:, a], M[:, b]
        for g, sub in [("pooled", None), ("within_ref", isref), ("within_alt", isalt)]:
            m = ~np.isnan(va) & ~np.isnan(vb)
            if sub is not None:
                m = m & sub
            if int(m.sum()) >= 10 and va[m].std() > 1e-9 and vb[m].std() > 1e-9:
                corr[g][binb].append(float(np.corrcoef(va[m], vb[m])[0, 1]))


def summ(l):
    a = np.array(l, float)
    return {"n": int(a.size), "mean": round(float(np.mean(a)), 4) if a.size else None,
            "median": round(float(np.median(a)), 4) if a.size else None}


out = {
    "n_regions": nreg, "n_reads": nreads,
    "scope_caveat": "regions = somatic-SNV-anchored (V6_on_fp FP-candidates + V3F_on_tp TP); ALLELE axis can be germline-allelic-confounded (orthogonal to the window question).",
    "metric_note": "read-pair distance = mean |methyl prob_i - prob_j| over >=3 common CpG (ISM C_min=3). discrimination = between-allele - within-allele (>0 => methylation separates alleles).",
    "discrimination_per_region": {
        str(W): {"mean_between_minus_within": round(float(np.mean(disc[W])), 4) if disc[W] else None,
                 "median": round(float(np.median(disc[W])), 4) if disc[W] else None,
                 "frac_regions_between_gt_within": round(disc_pos[W] / len(disc[W]), 4) if disc[W] else None,
                 "n_regions_evaluable": len(disc[W])} for W in WINDOWS},
    "pooled_dist": {str(W): {"within_allele": summ(within_all[W]), "between_allele": summ(between_all[W])} for W in WINDOWS},
    "RELIABILITY_corr_by_group": {g: {f"{lo}-{hi}bp": summ(corr[g][(lo, hi)]) for (lo, hi) in DIST_BINS} for g in corr},
}
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2, ensure_ascii=False))
