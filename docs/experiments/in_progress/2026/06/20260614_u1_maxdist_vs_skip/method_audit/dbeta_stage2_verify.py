#!/usr/bin/env python3
"""Δβ module stage-2 verification.

Quantifies the three stage-2 columns:
  - GermlineAsmDbeta  (normal HP1f − HP2f, goal 1) + cross-check vs Normal_HP_Signed_Delta
  - SubcloneDbeta_HP1 (tumor HP1 germline − HP1-1 carrier, goal 3)
  - SubcloneDbeta_HP2 (tumor HP2 germline − HP2-1 carrier, goal 3)
The fine-subclone "valid" count = how many regions have BOTH germline + somatic-carrier
reads on the same haplotype in tumor (= the testable same-hap subclone subset).
Output JSON only (no report writing).
"""
import csv
import json
import math
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_chr1_s2/significance_summary.csv"


def fnum(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


def col(rows, name):
    return [fnum(r.get(name, "")) for r in rows]


def sigcol(rows, name):
    return [r.get(name, "") == "true" for r in rows]


def pct(a, b):
    return round(a / b, 4) if b else None


def describe(vals, sigs, n):
    valid = [i for i in range(n) if vals[i] is not None]
    nv = len(valid)
    nsig = sum(1 for i in valid if sigs[i])
    vv = sorted(vals[i] for i in valid)
    med = vv[len(vv) // 2] if vv else None
    p10 = vv[int(len(vv) * 0.1)] if vv else None
    p90 = vv[int(len(vv) * 0.9)] if vv else None
    return {
        "n_valid": nv,
        "valid_frac_of_all": pct(nv, n),
        "n_sig": nsig,
        "sig_frac_of_valid": pct(nsig, nv),
        "median_signed": round(med, 4) if med is not None else None,
        "p10": round(p10, 4) if p10 is not None else None,
        "p90": round(p90, 4) if p90 is not None else None,
    }


def corr(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return None, 0
    xv = [p[0] for p in pairs]
    yv = [p[1] for p in pairs]
    mx, my = sum(xv) / len(xv), sum(yv) / len(yv)
    cov = sum((x - mx) * (y - my) for x, y in pairs)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xv))
    sy = math.sqrt(sum((y - my) ** 2 for y in yv))
    return (round(cov / (sx * sy), 4) if sx > 0 and sy > 0 else None), len(pairs)


rows = list(csv.DictReader(open(path)))
n = len(rows)

germ = col(rows, "GermlineAsmDbeta")
germ_sig = sigcol(rows, "GermlineAsmDbeta_Sig")
normal_signed = col(rows, "Normal_HP_Signed_Delta")  # cell-weighted legacy; read-weighted germ should correlate high
sub1 = col(rows, "SubcloneDbeta_HP1")
sub1_sig = sigcol(rows, "SubcloneDbeta_HP1_Sig")
sub2 = col(rows, "SubcloneDbeta_HP2")
sub2_sig = sigcol(rows, "SubcloneDbeta_HP2_Sig")

germ_corr, germ_npair = corr(germ, normal_signed)
# union of subclone testable regions (either HP)
sub_any = sum(1 for i in range(n) if sub1[i] is not None or sub2[i] is not None)
sub_any_sig = sum(1 for i in range(n) if (sub1[i] is not None and sub1_sig[i]) or (sub2[i] is not None and sub2_sig[i]))

out = {
    "source": path,
    "n_regions": n,
    "germline_asm": describe(germ, germ_sig, n),
    "germline_asm_corr_vs_normal_signed_delta": germ_corr,
    "germline_asm_corr_npair": germ_npair,
    "subclone_hp1": describe(sub1, sub1_sig, n),
    "subclone_hp2": describe(sub2, sub2_sig, n),
    "subclone_testable_any_hp": sub_any,
    "subclone_testable_any_hp_frac": pct(sub_any, n),
    "subclone_sig_any_hp": sub_any_sig,
}
print(json.dumps(out, indent=2, ensure_ascii=False))
