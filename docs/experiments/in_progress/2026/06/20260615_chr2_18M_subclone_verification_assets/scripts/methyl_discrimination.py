#!/usr/bin/env python3
"""
Quantify methylation discrimination between the two subclone lineages.
Lineage is defined by the variant pattern, NOT by HP tag, so this works for
both HKU_MOD (haplotagged) and DORADO (raw, no HP).

Lineage assignment per read (using covered SNV bases):
  alpha (group 1/2):  pos3 == A   (ref at 1,2 ; variant A at 3)
  beta  (group 3/4/5): pos1==G or pos2==C or pos3==G-with-1/2-variant
We use a robust rule:
  - if base(3)=='A'                       -> alpha
  - elif base(1)=='G' or base(2)=='C'     -> beta
  - elif base(3)=='G' and (base(1)=='C' and base(2)=='G') -> ambiguous(skip: could be HP2 ref)
For methylation we report, per CpG site and per lineage:
  n reads with a call, mean 5mC prob, %H (prob>=0.5)
Usage: methyl_discrimination.py <matrix_json> [matrix_json2 ...]
"""
import sys, json
from collections import defaultdict

CPG_ORDER = ["2.1", "2.2", "3.1", "3.2", "3.3", "3.4", "3.5", "4.1", "5.1", "5.2"]


def lineage(r):
    # 3.x CpG sites are physically adjacent to pos3 (4-10kb); pos3 allele (A vs G)
    # is the natural local marker linking them. pos3=A -> alpha(grp1/2); pos3=G -> beta(grp3/4/5).
    b3 = r["bases"].get("3")
    if b3 == "A":
        return "alpha"
    if b3 == "G":
        return "beta"
    return None


def summarize(path):
    j = json.load(open(path))
    label = j["label"]
    # dedup tumor reads by name (tagged BAM may duplicate records)
    seen = set(); reads = []
    for r in j["tumor_reads"]:
        if r["name"] in seen:
            continue
        seen.add(r["name"]); reads.append(r)
    by = defaultdict(lambda: defaultdict(list))   # site -> lineage -> [m_prob]
    lin_count = defaultdict(int)
    for r in reads:
        lin = lineage(r)
        if lin is None:
            continue
        lin_count[lin] += 1
        for s in CPG_ORDER:
            m = r["meth"].get(s)
            if m and m.get("m_prob") is not None:
                by[s][lin].append(m["m_prob"])
    print(f"\n################ {label}  (tumor unique reads={len(reads)}) ################")
    print(f"lineage read counts: " + ", ".join(f"{k}={v}" for k, v in sorted(lin_count.items())))
    print(f"{'CpG':>5} | {'alpha(1/2) meanM  %H  n':>26} | {'beta(3/4/5) meanM  %H  n':>26} | flip?")
    print("-" * 78)
    for s in CPG_ORDER:
        cells = {}
        for lin in ("alpha", "beta"):
            v = by[s].get(lin, [])
            if v:
                mean = sum(v) / len(v)
                ph = 100 * sum(1 for x in v if x >= 0.5) / len(v)
                cells[lin] = (mean, ph, len(v))
            else:
                cells[lin] = None
        def fmt(c):
            return f"{c[0]:.2f}   {c[1]:3.0f}%  {c[2]:>3}" if c else f"{'--':>16}"
        a = cells["alpha"]; b = cells["beta"]
        flip = ""
        if a and b:
            ah = a[0] >= 0.5; bh = b[0] >= 0.5
            flip = "FLIP*" if ah != bh else ("same" if abs(a[0]-b[0]) < 0.2 else "diff")
        print(f"{s:>5} | {fmt(a):>26} | {fmt(b):>26} | {flip}")


for p in sys.argv[1:]:
    summarize(p)
