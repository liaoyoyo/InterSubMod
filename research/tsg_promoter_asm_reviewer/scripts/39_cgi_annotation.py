#!/usr/bin/env python3
"""
39_cgi_annotation.py — does the methylation difference fall in CpG islands?

For key loci, annotate EACH differential CpG with (a) per-CpG d_HP (HP1-1 - HP1),
(b) in-CGI membership, (c) distance to the SNV. Answers: are the high-|Δβ| CpGs
in a CGI (regulatory cis-plausible) or outside (copy/structural)?  Also gives the
spatial profile (|Δβ| vs distance from SNV) as a by-product.
"""
import sys, glob, gzip
import numpy as np
sys.path.insert(0, "pipeline")
from lib.cis_asm_core import load_level1_plus, collapse_modtype, _beta_cpg, MIN_N
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"

LOCI = [("chr17", 79991120, "clean-survivor"),
        ("chr13", 32315128, "BRCA2"),
        ("chr4", 133868344, "strongest d_cis=0.706"),
        ("chr20", 59439285, "strongest ratio=10.9"),
        ("chr7", 79941963, "strong dodged")]

# load CGI intervals per chrom (UCSC cpgIslandExt: col2=chrom col3=0-based start col4=end col5=name)
CGI = {}
with gzip.open(f"{ROOT}/data/cpgIslandExt_hg38.txt.gz", "rt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        CGI.setdefault(p[1], []).append((int(p[2]) + 1, int(p[3]), p[4]))  # 1-based inclusive
for c in CGI:
    CGI[c].sort()


def cgi_hit(chrom, pos):
    for s, e, name in CGI.get(chrom, []):
        if s <= pos <= e:
            return name, 0
    # nearest distance
    best = None
    for s, e, name in CGI.get(chrom, []):
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if best is None or d < best[1]:
            best = (name, d)
    return (None, best[1]) if best else (None, -1)


def per_cpg_dHP(D, spos):
    G = _beta_cpg(D, spos, "tumor", {"1"})
    S = _beta_cpg(D, spos, "tumor", {"1-1"})
    out = []
    for c in set(G) & set(S):
        if G[c][1] >= MIN_N and S[c][1] >= MIN_N:
            out.append((int(c), S[c][0] - G[c][0]))   # cpg_pos(1-based), dHP
    return sorted(out)


print(f"{'locus':22s}{'tag':20s}{'SNV in CGI?':>14}{'nearest CGI dist':>17}")
for chrom, pos, tag in LOCI:
    cache = glob.glob(f"{ROOT}/pipeline/cache/level1/{chrom}_{pos}_*.tsv.gz")
    if not cache:
        print(f"{chrom}:{pos:14d}  NO CACHE"); continue
    name, dist = cgi_hit(chrom, pos)
    snv_str = f"YES ({name})" if name else f"NO"
    print(f"\n{chrom}:{pos:<11d}{tag:20s}{snv_str:>14}{dist:>17}")
    Dp = load_level1_plus(cache[0]); D = collapse_modtype(Dp)
    cpgs = per_cpg_dHP(D, str(pos))
    if not cpgs:
        print("   (no per-CpG d_HP computable)"); continue
    # annotate each CpG
    rows = []
    for cp, dhp in cpgs:
        nm, dd = cgi_hit(chrom, cp)
        rows.append((cp, dhp, nm is not None, dd, abs(cp - pos)))
    n_cgi = sum(1 for r in rows if r[2])
    # weight: does |dHP| concentrate in CGI vs non-CGI?
    in_cgi = [abs(r[1]) for r in rows if r[2]]
    out_cgi = [abs(r[1]) for r in rows if not r[2]]
    print(f"   CpGs total={len(rows)}  in-CGI={n_cgi}  out-CGI={len(rows)-n_cgi}")
    if in_cgi:
        print(f"   mean|dHP| in-CGI  = {np.mean(in_cgi):.3f}  (n={len(in_cgi)})")
    if out_cgi:
        print(f"   mean|dHP| out-CGI = {np.mean(out_cgi):.3f}  (n={len(out_cgi)})")
    # top-5 |dHP| CpGs: where are they?
    top = sorted(rows, key=lambda r: -abs(r[1]))[:5]
    print("   top-5 |dHP| CpGs (pos, dHP, in_CGI, dist_to_CGI, dist_to_SNV):")
    for cp, dhp, inc, dd, dsnv in top:
        print(f"      {cp}  dHP={dhp:+.2f}  CGI={'IN' if inc else f'{dd}bp away'}  SNV_dist={dsnv}bp")
