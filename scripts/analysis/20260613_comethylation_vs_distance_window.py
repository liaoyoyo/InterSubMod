#!/usr/bin/env python3
"""Co-methylation correlation vs genomic distance + per-read / per-pair CpG-call density
in +-1000bp vs +-5000bp windows, from existing phaseC region methylation matrices.

Rationale (window-size decision basis):
  - co-methylation correlation is local (literature: MHB ~r2>=0.5, ~1-2kb); a wide window
    dilutes read-read distance with uncorrelated distant CpGs.
  - BUT a narrow window risks < min_common_coverage(=3) CpGs => distance uncomputable.
  Existing phaseC regions span +-5000 (10kb); we sub-window to +-1000 to measure BOTH.

Output: research/flagship_chr2_18086020_20260612/comethyl_window_analysis.json
  (verified numbers for the paper window-size narrative; cite literature for mechanism.)
"""
import glob, random, csv, json, os, re, itertools
import numpy as np

BASE = "research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance"
PATTERNS = [
    f"{BASE}/V6_on_fp/filtered_snv_fp/chr*/chr*/chr*_*/methylation/methylation.csv",
    f"{BASE}/V3F_on_tp/filtered_snv_tp/chr*/chr*/chr*_*/methylation/methylation.csv",
]
OUT = "research/flagship_chr2_18086020_20260612/comethyl_window_analysis.json"
N_REGIONS = 600
random.seed(42)

ANCHOR_RE = re.compile(r"/chr[0-9XY]+_(\d+)/chr[0-9XY]+_\d+_\d+/methylation")

def anchor_from_path(p):
    m = ANCHOR_RE.search(p)
    return int(m.group(1)) if m else None

files = []
for pat in PATTERNS:
    files += glob.glob(pat)
random.shuffle(files)
files = files[:N_REGIONS]

DIST_BINS = [(0, 100), (100, 300), (300, 600), (600, 1000), (1000, 2000), (2000, 3500), (3500, 5000)]
corr = {b: [] for b in DIST_BINS}
read_tot = 0
read_lt3_1000 = read_lt3_5000 = 0
calls_1000 = []
calls_5000 = []
pair_tot = 0
pair_lt3_1000 = pair_lt3_5000 = 0
ref_lt3_1000 = ref_lt3_5000 = 0
ref_n_1000 = []
ref_n_5000 = []
nreg = 0

for f in files:
    anchor = anchor_from_path(f)
    if anchor is None:
        continue
    try:
        with open(f) as fh:
            rdr = csv.reader(fh)
            header = next(rdr)
            positions = np.array([int(x) for x in header[1:]])
            rows = [[np.nan if x == "NA" else float(x) for x in r[1:]] for r in rdr]
    except Exception:
        continue
    if not rows or positions.size == 0:
        continue
    M = np.array(rows, dtype=float)  # reads x cpg
    if M.ndim != 2 or M.shape[1] != positions.size:
        continue
    nreg += 1
    in1000 = np.abs(positions - anchor) <= 1000
    in5000 = np.abs(positions - anchor) <= 5000

    # reference (covered) CpG count per window = number of CpG columns in window
    n1, n5 = int(in1000.sum()), int(in5000.sum())
    ref_n_1000.append(n1); ref_n_5000.append(n5)
    ref_lt3_1000 += (n1 < 3); ref_lt3_5000 += (n5 < 3)

    notna = ~np.isnan(M)
    c1 = notna[:, in1000].sum(axis=1)
    c5 = notna[:, in5000].sum(axis=1)
    read_tot += M.shape[0]
    read_lt3_1000 += int((c1 < 3).sum())
    read_lt3_5000 += int((c5 < 3).sum())
    calls_1000 += c1.tolist(); calls_5000 += c5.tolist()

    # per-read-pair common coverage (the min_common_coverage=3 metric)
    nr = M.shape[0]
    if nr >= 2:
        pairs = list(itertools.combinations(range(nr), 2))
        random.shuffle(pairs)
        for i, j in pairs[:300]:
            both = notna[i] & notna[j]
            common1 = int((both & in1000).sum())
            common5 = int((both & in5000).sum())
            pair_tot += 1
            pair_lt3_1000 += (common1 < 3); pair_lt3_5000 += (common5 < 3)

    # CpG-pair correlation vs genomic distance
    ncol = M.shape[1]
    colpairs = list(itertools.combinations(range(ncol), 2))
    random.shuffle(colpairs)
    for a, b in colpairs[:600]:
        d = abs(int(positions[a]) - int(positions[b]))
        for lo, hi in DIST_BINS:
            if lo <= d < hi:
                va, vb = M[:, a], M[:, b]
                mask = ~np.isnan(va) & ~np.isnan(vb)
                if mask.sum() >= 10:
                    x, y = va[mask], vb[mask]
                    if x.std() > 1e-9 and y.std() > 1e-9:
                        corr[(lo, hi)].append(float(np.corrcoef(x, y)[0, 1]))
                break


def summ(lst):
    a = np.array(lst, dtype=float)
    if a.size == 0:
        return {"n": 0, "mean": None, "median": None}
    return {"n": int(a.size), "mean": round(float(np.mean(a)), 4), "median": round(float(np.median(a)), 4)}


out = {
    "n_regions": nreg, "n_reads": read_tot, "n_pairs_sampled": pair_tot,
    "note": "phaseC regions span +-5000 (10kb); +-1000 = inner sub-window. CpG cols = covered/called CpGs (union over reads), not all reference CpGs.",
    "P_read_lt3_calls": {"pm1000": round(read_lt3_1000 / read_tot, 4) if read_tot else None,
                          "pm5000": round(read_lt3_5000 / read_tot, 4) if read_tot else None},
    "P_pair_lt3_common_cov": {"pm1000": round(pair_lt3_1000 / pair_tot, 4) if pair_tot else None,
                               "pm5000": round(pair_lt3_5000 / pair_tot, 4) if pair_tot else None,
                               "comment": "min_common_coverage=3: pair distance UNCOMPUTABLE if common CpG < 3"},
    "P_region_lt3_cpg": {"pm1000": round(ref_lt3_1000 / nreg, 4) if nreg else None,
                          "pm5000": round(ref_lt3_5000 / nreg, 4) if nreg else None},
    "calls_per_read": {"pm1000": summ(calls_1000), "pm5000": summ(calls_5000)},
    "covered_cpg_per_region": {"pm1000": summ(ref_n_1000), "pm5000": summ(ref_n_5000)},
    "corr_vs_distance": {f"{lo}-{hi}bp": summ(corr[(lo, hi)]) for (lo, hi) in DIST_BINS},
}
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2, ensure_ascii=False))
