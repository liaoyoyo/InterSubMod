#!/usr/bin/env python3
"""Genome-wide per-locus x germline-HP methylation observation scan (HCC1395,
2026-01 w5000 archive). One row per (locus, germline HP family).

All statistics are computed on the COMMON CpG set (CpGs valid in >=90% of the
unit's reads) so that read-span differences cannot manufacture structure.
Emits raw statistics only; flags are derived downstream against the null.
"""
import os, csv, sys, math, random, statistics as st
from multiprocessing import Pool

A = "/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp"
OUT = sys.argv[1]
NPERM = 2
MINR = 6

COLS = ["chrom", "pos", "hp", "n_alt", "n_ref", "n_cpg", "n_common",
        "mean_alt", "mean_ref", "d_altref", "p_altref",
        "mabs_obs", "mabs_null", "exc_altref",
        "sep_alt", "dsplit_alt", "nlo_alt", "sep_ref", "sep_null",
        "strand_p", "cpg_lo", "cpg_hi"]


def opt_split(v, minc=3):
    v = sorted(v); n = len(v)
    if n < 2 * minc: return None
    best = None
    for k in range(minc, n - minc + 1):
        lo, hi = v[:k], v[k:]
        w = (len(lo) * st.pstdev(lo) ** 2 + len(hi) * st.pstdev(hi) ** 2) / n
        if best is None or w < best[0]: best = (w, k, lo, hi)
    w, k, lo, hi = best
    sd = max(math.sqrt(w), 1e-6)
    return ((st.mean(hi) - st.mean(lo)) / sd, st.mean(hi) - st.mean(lo), len(lo))


def germ(h):
    return h.split("-")[0] if h and h not in ("NA", ".") else None


def mwu_p(x, y):
    """Two-sided Mann-Whitney U normal approximation with tie correction."""
    n1, n2 = len(x), len(y)
    if n1 < 3 or n2 < 3: return 1.0
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks = [0.0] * len(allv); i = 0; ties = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]: j += 1
        r = (i + j) / 2.0 + 1
        t = j - i + 1
        ties += t ** 3 - t
        for q in range(i, j + 1): ranks[q] = r
        i = j + 1
    r1 = sum(ranks[q] for q in range(len(allv)) if allv[q][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 / 12.0 * ((n + 1) - ties / float(n * (n - 1))))
    if sd <= 0: return 1.0
    z = (abs(u1 - mu) - 0.5) / sd
    return math.erfc(z / math.sqrt(2))


def work(item):
    ch, l = item
    rng = random.Random(hash(l) & 0xFFFF)
    d = os.path.join(A, ch, l)
    try:
        sub = os.listdir(d)[0]
    except Exception:
        return []
    b = os.path.join(d, sub)
    fr, fm = os.path.join(b, "reads", "reads.tsv"), os.path.join(b, "methylation", "methylation.csv")
    if not (os.path.isfile(fr) and os.path.isfile(fm)): return []
    try:
        pos = int(l.rsplit("_", 1)[1])
    except Exception:
        return []
    info = {}
    try:
        with open(fr) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("is_tumor") != "1": continue
                g = germ(r.get("hp"))
                if g in ("1", "2") and r.get("alt_support") in ("ALT", "REF"):
                    info[r["read_id"]] = (g, r["alt_support"], r.get("strand", ""))
    except Exception:
        return []
    if len(info) < 2 * MINR: return []
    try:
        with open(fm) as fh:
            rd = csv.reader(fh); hdr = next(rd)
            cpg = [int(x) for x in hdr[1:]]
            mat = {r[0]: r[1:] for r in rd if r[0] in info}
    except Exception:
        return []
    ncpg = len(cpg)
    out = []
    for g in ("1", "2"):
        Aid = [k for k, v in info.items() if v[0] == g and v[1] == "ALT" and k in mat]
        Rid = [k for k, v in info.items() if v[0] == g and v[1] == "REF" and k in mat]
        if len(Aid) < MINR or len(Rid) < MINR: continue
        allid = Aid + Rid
        cov = [sum(1 for k in allid if mat[k][i] not in ("NA", "")) for i in range(ncpg)]
        common = [i for i in range(ncpg) if cov[i] >= 0.9 * len(allid)]
        if len(common) < 8: continue

        def pm(ids):
            o = {}
            for k in ids:
                v = [float(mat[k][i]) for i in common if mat[k][i] not in ("NA", "")]
                if len(v) >= max(4, 0.5 * len(common)): o[k] = st.mean(v)
            return o

        ma, mr = pm(Aid), pm(Rid)
        if len(ma) < MINR or len(mr) < MINR: continue
        va, vr = list(ma.values()), list(mr.values())
        # per-CpG |delta| observed vs permutation null
        obs, nul = [], []
        pool = list(ma) + list(mr); na = len(ma)
        perms = []
        for _ in range(NPERM):
            p = pool[:]; rng.shuffle(p); perms.append((set(p[:na]), set(p[na:])))
        for i in common:
            xa = [float(mat[k][i]) for k in ma if mat[k][i] not in ("NA", "")]
            xr = [float(mat[k][i]) for k in mr if mat[k][i] not in ("NA", "")]
            if len(xa) < 3 or len(xr) < 3: continue
            obs.append(abs(st.mean(xa) - st.mean(xr)))
            for pa, pr in perms:
                ya = [float(mat[k][i]) for k in pool if k in pa and mat[k][i] not in ("NA", "")]
                yr = [float(mat[k][i]) for k in pool if k in pr and mat[k][i] not in ("NA", "")]
                if len(ya) >= 3 and len(yr) >= 3:
                    nul.append(abs(st.mean(ya) - st.mean(yr)))
        if not obs or not nul: continue
        sa, sr = opt_split(va), opt_split(vr)
        if not (sa and sr): continue
        mu, sd = st.mean(va), max(st.pstdev(va), 1e-6)
        ns = [opt_split([rng.gauss(mu, sd) for _ in va]) for _ in range(3)]
        ns = [x[0] for x in ns if x]
        # strand association with the ALT split
        thr = sorted(va)[sa[2] - 1]
        hi = [k for k, v in ma.items() if v > thr]; lo = [k for k, v in ma.items() if v <= thr]
        fw = lambda ids: sum(1 for k in ids if info[k][2] in ("+", "forward", "F", "1"))
        a11, a12, a21, a22 = fw(hi), len(hi) - fw(hi), fw(lo), len(lo) - fw(lo)
        try:
            from scipy.stats import fisher_exact
            sp = fisher_exact([[a11, a12], [a21, a22]])[1]
        except Exception:
            sp = 1.0
        out.append([ch, pos, g, len(ma), len(mr), ncpg, len(common),
                    round(st.mean(va), 4), round(st.mean(vr), 4),
                    round(st.mean(va) - st.mean(vr), 4), "%.3e" % mwu_p(va, vr),
                    round(st.median(obs), 4), round(st.median(nul), 4),
                    round(st.median(obs) - st.median(nul), 4),
                    round(sa[0], 3), round(sa[1], 4), sa[2],
                    round(sr[0], 3), round(st.median(ns) if ns else 0, 3),
                    "%.3e" % sp, min(cpg[i] for i in common), max(cpg[i] for i in common)])
    return out


if __name__ == "__main__":
    locs = []
    for ch in sorted(os.listdir(A)):
        cd = os.path.join(A, ch)
        if os.path.isdir(cd):
            for l in os.listdir(cd): locs.append((ch, l))
    sys.stderr.write("loci=%d\n" % len(locs)); sys.stderr.flush()
    n = 0
    with open(OUT, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        with Pool(12) as p:
            for i, rows in enumerate(p.imap_unordered(work, locs, chunksize=40)):
                for r in rows:
                    fh.write("\t".join(str(x) for x in r) + "\n"); n += 1
                if i % 3000 == 0:
                    sys.stderr.write("scanned %d/%d units=%d\n" % (i, len(locs), n)); sys.stderr.flush()
    sys.stderr.write("DONE units=%d\n" % n)
