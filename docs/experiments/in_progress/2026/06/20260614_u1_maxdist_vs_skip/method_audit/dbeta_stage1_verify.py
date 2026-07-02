#!/usr/bin/env python3
"""Δβ module stage-1 completion verification.

Reads significance_summary.csv (must contain SomaticResidualDbeta columns) and
quantifies the somatic residual Δβ distribution vs the legacy (flawed) hp_residual_sig,
to confirm #3 修法 lands correctly. Output JSON only (no report writing).
"""
import csv
import json
import math
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_chr1/significance_summary.csv"


def fnum(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


rows = list(csv.DictReader(open(path)))
n = len(rows)

# Stage 1 columns
db = [fnum(r.get("SomaticResidualDbeta", "")) for r in rows]
dbp = [fnum(r.get("SomaticResidualDbeta_P", "")) for r in rows]
dbsig = [r.get("SomaticResidualDbeta_Sig", "") == "true" for r in rows]
# Legacy value + flawed sig for cross-check
sres = [fnum(r.get("HP_Signed_Residual", "")) for r in rows]
old_sig = [r.get("HP_Residual_Sig", "") == "true" for r in rows]

valid = [i for i in range(n) if db[i] is not None]
nv = len(valid)

# sig fraction (over valid)
new_sig_valid = sum(1 for i in valid if dbsig[i])
old_sig_all = sum(1 for s in old_sig if s)
old_sig_valid = sum(1 for i in valid if old_sig[i])

# correlation Δβ vs legacy HP_Signed_Residual (both present)
both = [i for i in valid if sres[i] is not None]
if len(both) >= 3:
    xs = [db[i] for i in both]
    ys = [sres[i] for i in both]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    corr = cov / (sx * sy) if sx > 0 and sy > 0 else None
    max_abs_diff = max(abs(x - y) for x, y in zip(xs, ys))
else:
    corr, max_abs_diff = None, None

# p-value distribution
pv = [dbp[i] for i in valid if dbp[i] is not None]
p_min = min(pv) if pv else None
p_med = sorted(pv)[len(pv) // 2] if pv else None
# floor check: 999-perm -> min p = 1/1000 = 0.001
p_at_floor = sum(1 for p in pv if p is not None and abs(p - 0.001) < 1e-9)

# Δβ magnitude distribution (signed)
dbvals = sorted(db[i] for i in valid)
db_med = dbvals[len(dbvals) // 2] if dbvals else None
db_p10 = dbvals[int(len(dbvals) * 0.1)] if dbvals else None
db_p90 = dbvals[int(len(dbvals) * 0.9)] if dbvals else None

out = {
    "source": path,
    "n_regions": n,
    "n_valid_dbeta": nv,
    "valid_frac": round(nv / n, 4) if n else None,
    "new_sig_count_valid": new_sig_valid,
    "new_sig_frac_valid": round(new_sig_valid / nv, 4) if nv else None,
    "old_hp_residual_sig_count_all": old_sig_all,
    "old_hp_residual_sig_frac_all": round(old_sig_all / n, 4) if n else None,
    "old_hp_residual_sig_count_valid": old_sig_valid,
    "corr_dbeta_vs_legacy_signed_residual": round(corr, 4) if corr is not None else None,
    "max_abs_diff_dbeta_vs_legacy": round(max_abs_diff, 6) if max_abs_diff is not None else None,
    "n_both_for_corr": len(both),
    "p_min": p_min,
    "p_median": p_med,
    "p_at_floor_0.001_count": p_at_floor,
    "dbeta_median_signed": round(db_med, 4) if db_med is not None else None,
    "dbeta_p10": round(db_p10, 4) if db_p10 is not None else None,
    "dbeta_p90": round(db_p90, 4) if db_p90 is not None else None,
}
print(json.dumps(out, indent=2, ensure_ascii=False))
