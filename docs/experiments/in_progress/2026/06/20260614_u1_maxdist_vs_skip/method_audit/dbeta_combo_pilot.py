#!/usr/bin/env python3
"""Pilot: Δβ label-combination (A alt-axis subclone / B full-combo / C per-CpG localization) + timing.

Reuses per-region reads.tsv (hp, alt_support, is_tumor) + methylation.csv (read x CpG β).
Vectorized permutation (numpy) for speed. Times each component to extrapolate cpp-change cost.

A: alt-axis subclone — tumor HP1(1/1-1) split by alt=ALT vs REF; compare sig to HP-tag subclone.
B: full-combo — groups = sample×HP-family×alt with >=min; all pairwise Δβ + perm + BH-FDR; count combos/region.
C: per-CpG localization — for HP-tag subclone-sig regions, per-CpG germline-vs-carrier mean diff; how localized.
"""
import csv
import glob
import json
import os
import sys
import time

import warnings

import numpy as np

warnings.filterwarnings("ignore")

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_chr1_v2/filtered_snv_tp_chr1/chr1"
SUMMARY = sys.argv[2] if len(sys.argv) > 2 else "/tmp/dbeta_chr1_v2/significance_summary.csv"
NPERM = int(sys.argv[3]) if len(sys.argv) > 3 else 999
MIN_GROUP = 3
SEED = 42
rng = np.random.default_rng(SEED)


def dbeta_perm(mean_beta, grp):
    """Δβ = mean(grp0) − mean(grp1) + vectorized two-sided perm p. grp = bool array (True=grp0)."""
    n0 = int(grp.sum())
    n1 = len(grp) - n0
    if min(n0, n1) < 1:
        return np.nan, 1.0, n0, n1
    obs = mean_beta[grp].mean() - mean_beta[~grp].mean()
    # vectorized: NPERM shuffles of the group assignment
    idx = np.argsort(rng.random((NPERM, len(grp))), axis=1)  # random permutations
    shuffled = grp[idx]  # (NPERM, n)
    m0 = (mean_beta[None, :] * shuffled).sum(1) / shuffled.sum(1)
    m1 = (mean_beta[None, :] * ~shuffled).sum(1) / (~shuffled).sum(1)
    perm = m0 - m1
    n_extreme = 1 + int((np.abs(perm) >= abs(obs)).sum())
    p = n_extreme / (NPERM + 1)
    return obs, p, n0, n1


def load_region(inner):
    rt = os.path.join(inner, "reads", "reads.tsv")
    mc = os.path.join(inner, "methylation", "methylation.csv")
    if not (os.path.exists(rt) and os.path.exists(mc)):
        return None
    meta = {}
    with open(rt) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            meta[row["read_id"]] = (row["hp"], row["alt_support"], row["is_tumor"] == "1")
    ids, betarows = [], []
    with open(mc) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            if p[0] not in meta:
                continue
            row = [np.nan if v in ("NA", "") else float(v) for v in p[1:]]
            ids.append(p[0])
            betarows.append(row)
    if not ids:
        return None
    B = np.array(betarows, float)  # reads x CpG
    mean_beta = np.nanmean(np.where(np.isnan(B), np.nan, B), axis=1)
    keep = ~np.isnan(mean_beta)
    ids = [ids[i] for i in range(len(ids)) if keep[i]]
    hp = np.array([meta[i][0] for i in ids])
    alt = np.array([meta[i][1] for i in ids])
    tum = np.array([meta[i][2] for i in ids])
    return {"mean_beta": mean_beta[keep], "B": B[keep], "hp": hp, "alt": alt, "tum": tum}


summ = {r["Pos"]: r for r in csv.DictReader(open(SUMMARY))}
region_dirs = sorted(glob.glob(os.path.join(ROOT, "chr1_*")))

# ---- accumulators ----
tA = tB = tC = 0.0
# A: alt-axis subclone vs HP-tag subclone concordance
A_cmp = []  # (hp, hptag_sig, altaxis_sig)
# B: combos per region + sig-combo discoveries beyond targeted-4
B_pairs_per_region = []
B_total_sig = 0  # all BH-sig combo pairs
B_altaxis_sig = 0  # NEW: same-sample same-family ALT-vs-REF (the alt axis HP-tag subclone cannot see)
B_altaxis_tested = 0
# C: per-CpG localization of HP-tag subclone-sig regions
C_localize = []  # (n_cpg_total, n_cpg_driving) for sig subclone regions

HP1F = lambda h: h in ("1", "1-1", "HP1")  # noqa: E731
HP2F = lambda h: h in ("2", "2-1", "HP2")  # noqa: E731

for d in region_dirs:
    pos = os.path.basename(d).split("_")[-1]
    inner = glob.glob(d + "/chr*")
    if not inner or pos not in summ:
        continue
    R = load_region(inner[0])
    if R is None or len(R["mean_beta"]) < 6:
        continue
    mb, hp, alt, tum = R["mean_beta"], R["hp"], R["alt"], R["tum"]

    # ===== A: alt-axis subclone (tumor HP-family, ALT vs REF) =====
    t0 = time.time()
    for fam, famfn, tagcol in [("HP1", HP1F, "SubcloneDbeta_HP1_Sig"), ("HP2", HP2F, "SubcloneDbeta_HP2_Sig")]:
        m = tum & np.array([famfn(h) for h in hp]) & np.isin(alt, ["ALT", "REF"])
        if m.sum() >= 2 * MIN_GROUP:
            sub_mb = mb[m]
            g = alt[m] == "ALT"  # True=ALT (carrier-by-allele), False=REF
            if MIN_GROUP <= int(g.sum()) and MIN_GROUP <= int((~g).sum()):
                _, p, n0, n1 = dbeta_perm(sub_mb, g)
                altaxis_sig = p <= 0.05
                hptag_sig = summ[pos].get(tagcol) == "true"
                A_cmp.append((fam, hptag_sig, altaxis_sig))
    tA += time.time() - t0

    # ===== B: full-combo pairwise Δβ (sample × HP-family × alt) =====
    t0 = time.time()
    groups = {}
    for s, sname in [(True, "T"), (False, "N")]:
        for famn, famfn in [("HP1", HP1F), ("HP2", HP2F)]:
            for a in ["ALT", "REF"]:
                m = (tum == s) & np.array([famfn(h) for h in hp]) & (alt == a)
                if m.sum() >= MIN_GROUP:
                    groups[f"{sname}.{famn}.{a}"] = mb[m]
    gkeys = list(groups)
    region_pairs = 0
    raw_p = []
    for i in range(len(gkeys)):
        for j in range(i + 1, len(gkeys)):
            a, b = groups[gkeys[i]], groups[gkeys[j]]
            obs = a.mean() - b.mean()
            allv = np.concatenate([a, b])
            lab = np.array([True] * len(a) + [False] * len(b))
            idx = np.argsort(rng.random((NPERM, len(allv))), axis=1)
            sh = lab[idx]
            perm = (allv[None, :] * sh).sum(1) / sh.sum(1) - (allv[None, :] * ~sh).sum(1) / (~sh).sum(1)
            p = (1 + int((np.abs(perm) >= abs(obs)).sum())) / (NPERM + 1)
            raw_p.append((gkeys[i], gkeys[j], obs, p))
            region_pairs += 1
    B_pairs_per_region.append(region_pairs)
    # count alt-axis tested pairs (same sample, same family, ALT vs REF)
    def is_altaxis(g1, g2):
        a, b = g1.split("."), g2.split(".")
        return a[0] == b[0] and a[1] == b[1] and {a[2], b[2]} == {"ALT", "REF"}

    for g1, g2, _, _ in raw_p:
        if is_altaxis(g1, g2):
            B_altaxis_tested += 1
    # BH-FDR within region
    if raw_p:
        ps = sorted(range(len(raw_p)), key=lambda k: raw_p[k][3])
        mtot = len(raw_p)
        for rank, k in enumerate(ps, 1):
            q = raw_p[k][3] * mtot / rank
            if q <= 0.05:
                B_total_sig += 1
                if is_altaxis(raw_p[k][0], raw_p[k][1]):
                    B_altaxis_sig += 1
    tB += time.time() - t0

    # ===== C: per-CpG localization for HP-tag subclone-sig (HP1) =====
    if summ[pos].get("SubcloneDbeta_HP1_Sig") == "true":
        t0 = time.time()
        Bm = R["B"]
        gmask = tum & (hp == "1")
        cmask = tum & (hp == "1-1")
        if gmask.sum() >= MIN_GROUP and cmask.sum() >= MIN_GROUP:
            ncpg = Bm.shape[1]
            driving = 0
            for c in range(ncpg):
                gv = Bm[gmask, c]
                cv = Bm[cmask, c]
                gv = gv[~np.isnan(gv)]
                cv = cv[~np.isnan(cv)]
                if len(gv) >= MIN_GROUP and len(cv) >= MIN_GROUP:
                    if abs(gv.mean() - cv.mean()) > 0.2:  # localization proxy: substantial per-CpG diff
                        driving += 1
            C_localize.append((ncpg, driving))
        tC += time.time() - t0

# ---- summarize ----
def concord_A():
    if not A_cmp:
        return {}
    n = len(A_cmp)
    both = sum(1 for _, h, a in A_cmp if h and a)
    hptag_only = sum(1 for _, h, a in A_cmp if h and not a)
    alt_only = sum(1 for _, h, a in A_cmp if not h and a)
    neither = sum(1 for _, h, a in A_cmp if not h and not a)
    agree = (both + neither) / n
    return {"n_tests": n, "agreement": round(agree, 4), "both_sig": both,
            "hptag_only_sig": hptag_only, "altaxis_only_sig": alt_only, "neither": neither}


out = {
    "scope": f"{len(region_dirs)} region dirs, NPERM={NPERM}, MIN_GROUP={MIN_GROUP}",
    "A_alt_axis_subclone_vs_hptag": concord_A(),
    "B_full_combo": {
        "regions_with_combos": len(B_pairs_per_region),
        "median_pairs_per_region": int(np.median(B_pairs_per_region)) if B_pairs_per_region else 0,
        "max_pairs_per_region": int(np.max(B_pairs_per_region)) if B_pairs_per_region else 0,
        "total_BH_sig_combos": B_total_sig,
        "altaxis_pairs_tested": B_altaxis_tested,
        "altaxis_BH_sig": B_altaxis_sig,
    },
    "C_per_cpg_localization": {
        "n_subclone_sig_regions_localized": len(C_localize),
        "median_cpg_total": int(np.median([c[0] for c in C_localize])) if C_localize else 0,
        "median_cpg_driving(|Δβ|>0.2)": float(np.median([c[1] for c in C_localize])) if C_localize else 0,
        "frac_regions_localizable(>=1 driving CpG)": round(
            sum(1 for c in C_localize if c[1] >= 1) / len(C_localize), 4) if C_localize else 0,
    },
    "timing_python_seconds": {"A": round(tA, 2), "B": round(tB, 2), "C": round(tC, 2),
                              "total": round(tA + tB + tC, 2)},
}
print(json.dumps(out, indent=2, ensure_ascii=False))
