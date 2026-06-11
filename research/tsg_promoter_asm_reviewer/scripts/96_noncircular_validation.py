#!/usr/bin/env python3
"""
96 - NON-CIRCULAR, TP/FP-free validation that a locus's methylation is genuinely
associated with the (SNV-derived) HP tag — covering BOTH ASM modes (mean-shift AND
clustering) with a single multivariate statistic, validated by INDEPENDENT replication
(not the statistic that defined the tiers).

Statistic (mode-agnostic): on the read×CpG methylation matrix, pairwise read distance
  d(i,j) = mean |β_i - β_j| over shared covered CpGs (>=3 shared). HP-separation
  F = mean(between-HP d) / mean(within-HP d).  F>1 <=> reads differ by HP (shift OR
  cluster both raise it). Non-circular vs LabelHPPermanova: we recompute our OWN
  distance and, crucially, REPLICATE on independent CpG halves & strands.

Per locus (base HP = SNV-phased HP1/HP2 family):
  perm   : shuffle HP labels on the FIXED full-CpG distance -> empirical p on F (BH-FDR)
  splitA : recompute distance on EVEN CpGs only -> F_even + its own perm p
  splitB : recompute distance on ODD  CpGs only -> F_odd  + its own perm p
           (independent CpG subsets; both significant & F>1 = replicated, non-circular)
  strand : F on forward-only reads vs reverse-only reads -> same direction (F>1 both)
  boot   : resample reads, fraction with F>1 (stability)
  VALID  = perm_q<=0.05 & splitA_p<=0.05 & splitB_p<=0.05 & (strand ok or not measurable)
           & boot>=0.90
Output: display_v2/validation3.json
"""
import os, csv, json, glob
os.environ.setdefault("TMPDIR", "/big7_disk/liaoyoyo2001/tmp")
from multiprocessing import Pool
import numpy as np

DV = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample/display_v2")
SCAN = "/big7_disk/liaoyoyo2001/ism_display_scan"
N_PERM, N_BOOT = 1000, 500
MIN_SHARED, MIN_PER_GROUP = 3, 4


def region_dir(cls, chrom, pos):
    for d in glob.glob(f"{SCAN}/HCC1395_{cls}/curated_{cls}/{chrom}/{chrom}_{pos}/*"):
        if os.path.exists(f"{d}/methylation/methylation.csv"):
            return d
    return None


def pdist(M):
    """pairwise mean|Δβ| over shared covered CpGs; NaN if <MIN_SHARED shared."""
    n = M.shape[0]
    D = np.full((n, n), np.nan)
    valid = ~np.isnan(M)
    for i in range(n):
        di = M[i]; vi = valid[i]
        diff = np.abs(di - M)                     # (n, ncpg)
        shared = vi & valid                        # (n, ncpg)
        cnt = shared.sum(1)
        s = np.where(shared, diff, 0.0).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            row = np.where(cnt >= MIN_SHARED, s / cnt, np.nan)
        D[i] = row
    return D


def Fstat(D, lab, mask=None):
    n = D.shape[0]
    idx = np.arange(n) if mask is None else np.where(mask)[0]
    if len(idx) < 2:
        return np.nan
    L = lab[idx]
    sub = D[np.ix_(idx, idx)]
    iu = np.triu_indices(len(idx), 1)
    dvals = sub[iu]
    same = (L[iu[0]] == L[iu[1]])
    btw = dvals[~same]; win = dvals[same]
    btw = btw[np.isfinite(btw)]; win = win[np.isfinite(win)]
    if len(btw) < 1 or len(win) < 1 or np.mean(win) == 0:
        return np.nan
    return np.mean(btw) / np.mean(win)


def perm_p(D, lab, obs, rng, mask=None, nperm=N_PERM):
    if not np.isfinite(obs):
        return np.nan
    n = D.shape[0]
    sub_lab = lab if mask is None else lab[mask]
    cnt = 0
    for _ in range(nperm):
        pl = lab.copy()
        if mask is None:
            pl = rng.permutation(lab)
        else:
            idx = np.where(mask)[0]
            pl = lab.copy(); pl[idx] = rng.permutation(lab[idx])
        if (Fstat(D, pl, mask) or -1) >= obs:
            cnt += 1
    return (cnt + 1) / (nperm + 1)


def work(rec):
    cls, chrom, pos = rec["cls"], rec["chr"], rec["pos"]
    key = f"{chrom}_{pos}"
    rd = region_dir(cls, chrom, pos)
    if rd is None:
        return key, None
    try:
        rows = list(csv.reader(open(f"{rd}/methylation/methylation.csv")))
        rids = [int(r[0]) for r in rows[1:]]
        M = np.array([[np.nan if x == "NA" else float(x) for x in r[1:]] for r in rows[1:]])
        base, strand = {}, {}
        for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t"):
            h = r.get("hp", "")
            base[int(r["read_id"])] = 1 if h in ("1", "1-1") else (2 if h in ("2", "2-1") else 0)
            strand[int(r["read_id"])] = r.get("strand", "")
        lab = np.array([base.get(i, 0) for i in rids])
        st = np.array([strand.get(i, "") for i in rids])
        keep = lab > 0
        M, lab, st = M[keep], lab[keep], st[keep]
        if (lab == 1).sum() < MIN_PER_GROUP or (lab == 2).sum() < MIN_PER_GROUP:
            return key, dict(ok=False, reason="too few per HP")
        rng = np.random.RandomState((int(pos) * 131 + len(chrom)) % (2**31))
        ncol = M.shape[1]
        ev = np.arange(ncol) % 2 == 0
        Dfull, Dev, Dod = pdist(M), pdist(M[:, ev]), pdist(M[:, ~ev])
        F = Fstat(Dfull, lab); pF = perm_p(Dfull, lab, F, rng)
        Fe = Fstat(Dev, lab); pe = perm_p(Dev, lab, Fe, rng)
        Fo = Fstat(Dod, lab); po = perm_p(Dod, lab, Fo, rng)
        # strand
        def strand_F(s):
            m = st == s
            if (lab[m] == 1).sum() >= 3 and (lab[m] == 2).sum() >= 3:
                return Fstat(Dfull, lab, m)
            return np.nan
        Ff, Fr = strand_F("+"), strand_F("-")
        strand_meas = np.isfinite(Ff) and np.isfinite(Fr)
        strand_ok = bool(strand_meas and Ff > 1 and Fr > 1)
        # bootstrap
        idx = np.arange(len(lab)); bc = 0
        for _ in range(N_BOOT):
            s = rng.choice(idx, len(idx), replace=True)
            if (Fstat(Dfull[np.ix_(s, s)], lab[s]) or 0) > 1:
                bc += 1
        boot = bc / N_BOOT
        return key, dict(ok=True, F=round(float(F), 3), perm_p=round(float(pF), 4),
                         F_even=round(float(Fe), 3), pe=round(float(pe), 4),
                         F_odd=round(float(Fo), 3), po=round(float(po), 4),
                         strand_meas=bool(strand_meas), strand_ok=strand_ok,
                         F_fwd=None if not np.isfinite(Ff) else round(float(Ff), 3),
                         F_rev=None if not np.isfinite(Fr) else round(float(Fr), 3),
                         boot=round(boot, 3))
    except Exception as e:
        return key, dict(ok=False, reason=f"{type(e).__name__}: {e}")


def bh(p):
    p = np.array(p); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank, i in enumerate(reversed(o)):
        prev = min(prev, p[i] * n / (n - rank)); q[i] = prev
    return q


def main():
    curated = json.load(open(f"{DV}/curated_loci.json"))
    tier = json.load(open(f"{DV}/tier_assignment.json"))
    man = {m["chr"] + "_" + str(m["pos"]): m for m in json.load(open(f"{DV}/manifest.json"))}
    print(f"[96] non-circular validation of {len(curated)} loci ...")
    res = {}
    with Pool(12) as p:
        for i, (k, r) in enumerate(p.imap_unordered(work, curated, chunksize=6)):
            if r:
                res[k] = r
            if (i + 1) % 300 == 0:
                print(f"   ...{i+1}/{len(curated)}")
    ok = {k: v for k, v in res.items() if v.get("ok")}
    ks = list(ok); q = bh([ok[k]["perm_p"] for k in ks])
    for k, qq in zip(ks, q):
        v = ok[k]; v["perm_q"] = round(float(qq), 4)
        v["validated"] = bool(v["perm_q"] <= 0.05 and v["pe"] <= 0.05 and v["po"] <= 0.05
                              and (v["strand_ok"] or not v["strand_meas"]) and v["boot"] >= 0.90)
    def agg(keys):
        return dict(n=len(keys), validated=sum(1 for k in keys if ok.get(k, {}).get("validated")))
    by_tier = {t or "none": agg([k for k in ok if tier.get(k, "") == t]) for t in ("A", "B", "")}
    dbonly = [k for k in ok if abs(man.get(k, {}).get("db", 0)) >= 0.2 and tier.get(k, "") not in ("A", "B")]
    summary = dict(n_testable=len(ok), by_tier=by_tier,
                   dbeta_only=agg(dbonly),
                   criteria="perm_q<=0.05 & split-half(even&odd) both p<=0.05 & strand(if meas) & boot>=0.90 — all INDEPENDENT replication, non-circular, TP/FP-free")
    json.dump(dict(summary=summary, loci=ok), open(f"{DV}/validation3.json", "w"))
    print(f"\n[96] testable={len(ok)}")
    for t in ("A", "B", "none"):
        a = by_tier[t]; print(f"   Tier {t:4s}: n={a['n']:4d} INDEP-validated={a['validated']:4d} ({100*a['validated']/max(a['n'],1):.0f}%)")
    d = summary["dbeta_only"]; print(f"   Δβ-only: n={d['n']} indep-validated={d['validated']} ({100*d['validated']/max(d['n'],1):.0f}%)")
    print(f"[96] wrote {DV}/validation3.json")


if __name__ == "__main__":
    main()
