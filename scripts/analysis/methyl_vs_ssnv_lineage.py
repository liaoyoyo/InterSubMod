#!/usr/bin/env python3
"""Does methylation read-clustering align with / differ by the sSNV multi-locus lineage?

The sSNV co-occurrence lineage (from sm_linkage) is the non-circular subclone backbone.
This asks the direct question the CN / single-locus-HP axes could not: within a multi-sSNV
tree region, do reads of different sSNV-lineages (multi-locus genotype populations) carry
different methylation?

READ-ONLY pysam pileup of the tagged BAM (sSNV alleles + 5mC both live in the BAM -> no
external data, no BAM edit). Per read spanning a region:
  - multi-sSNV genotype across the region's sSNVs (R/A) -> population/lineage
  - mean 5mC methylation over the read's CpG calls (MM/ML)

Per region metrics:
  - delta_beta        : max pairwise |mean_methyl(pop_i)-mean_methyl(pop_j)|, pops >=MIN_POP
  - perm_p            : PERMUTATION null p (shuffle pop labels K times) -> proper significance
                        (replaces n-driven MWU false positives)
  - single_sSNV_max_db: max over each sSNV of allele-Δβ (ALT-reads vs REF-reads mean methyl)
                        = per-locus cis-ASM proxy; lineage_db > single_sSNV_max_db => signal
                        BEYOND single-locus cis (multi-lineage / subclone-level)
  - ari_2grp          : ARI(median-β 2-split, top-2 populations) = alignment proxy
FDR (BH) applied to perm_p across regions.

Usage: methyl_vs_ssnv_lineage.py <SAMPLE> <MSROOT> <TAGGED_BAM> [MAX_REGIONS] [SEED]
"""
import glob
import json
import math
import os
import sys
from collections import defaultdict

import pysam

QUAL_METH = 128
MIN_POP = 4
MIN_CPG = 3
N_PERM = 1000
MAX_SPAN = 100000


def load_census(msroot):
    d = json.load(open(f"{msroot}/sm_linkage_genomewide.json"))
    return {k: (v["ref"].upper(), v["alt"].upper()) for k, v in d["census"].items()}


def load_regions(msroot, max_regions):
    regs = []
    for f in sorted(glob.glob(f"{msroot}/ml_part_*.json")):
        for g in json.load(open(f))["groups"]:
            if (g["n_sSNV"] >= 2 and g.get("n_populations", 0) >= 2
                    and g.get("n_populations_with_ALT", 0) >= 2
                    and g["n_full_cov_reads"] >= 12
                    and (g["end"] - g["start"]) <= MAX_SPAN):
                regs.append(g)
    dedup = {(g["chrom"], g["start"], g["end"]): g for g in regs}
    regs = sorted(dedup.values(), key=lambda g: (-g["n_populations"], -g["n_full_cov_reads"]))
    return regs[:max_regions] if max_regions else regs


def read_methyl(read):
    try:
        mb = read.modified_bases
    except Exception:
        return None
    if not mb:
        return None
    quals = []
    for (canon, strand, mod), calls in mb.items():
        if mod == "m" or str(mod) in ("m", "C+m", "27551"):
            quals.extend(q for _, q in calls if q >= 0)
    if len(quals) < MIN_CPG:
        return None
    return sum(1 for q in quals if q >= QUAL_METH) / len(quals)


def max_pair_delta(groups):
    means = [sum(v) / len(v) for v in groups.values() if len(v) >= MIN_POP]
    return (max(means) - min(means)) if len(means) >= 2 else None


def perm_p(read_meth, read_pop, observed, rng, k=N_PERM):
    """empirical p: shuffle pop labels, recompute delta_beta."""
    vals = list(read_meth)
    labs = list(read_pop)
    n = len(vals)
    ge = 0
    for _ in range(k):
        rng.shuffle(labs)
        g = defaultdict(list)
        for i in range(n):
            g[labs[i]].append(vals[i])
        d = max_pair_delta(g)
        if d is not None and d >= observed - 1e-12:
            ge += 1
    return round((ge + 1) / (k + 1), 4)


def ari_2grp(vals, labs):
    """ARI between median-beta 2-split and the binary (top-2-pop) labels."""
    med = sorted(vals)[len(vals) // 2]
    pred = [1 if v > med else 0 for v in vals]
    # contingency
    from collections import Counter
    n = len(vals)
    a = Counter(zip(pred, labs))
    ai = Counter(pred)
    bi = Counter(labs)
    idx = sum(math.comb(v, 2) for v in a.values())
    ea = sum(math.comb(v, 2) for v in ai.values())
    eb = sum(math.comb(v, 2) for v in bi.values())
    exp = ea * eb / math.comb(n, 2) if n >= 2 else 0
    mx = (ea + eb) / 2.0
    return round((idx - exp) / (mx - exp), 3) if (mx - exp) else 0.0


def analyze_region(bam, cen, g, rng):
    chrom = g["chrom"]
    sites = [(p, cen.get(f"{chrom}:{p}")) for p in g["positions"]]
    sites = [(p, ra) for p, ra in sites if ra]
    if len(sites) < 2:
        return None
    lo, hi = g["start"], g["end"]
    site0 = {p - 1 for p, _ in sites}
    reads = {}
    for read in bam.fetch(chrom, lo - 1, hi + 1):
        if read.is_secondary or read.is_supplementary or read.is_unmapped or not read.query_sequence:
            continue
        pos2b = {}
        for qp, rp in read.get_aligned_pairs(matches_only=True):
            if rp in site0:
                pos2b[rp] = read.query_sequence[qp].upper()
                if len(pos2b) == len(site0):
                    break
        geno, per_site = [], {}
        ok = True
        for p, (ref, alt) in sites:
            b = pos2b.get(p - 1)
            if b is None:
                ok = False
                break
            call = "A" if b == alt else ("R" if b == ref else "?")
            geno.append(call)
            per_site[p] = call
        if not ok or "?" in geno:
            continue
        m = read_methyl(read)
        if m is None:
            continue
        reads[read.query_name] = {"geno": "".join(geno), "methyl": m, "sites": per_site}
    if len(reads) < 8:
        return None
    pops = defaultdict(list)
    for v in reads.values():
        pops[v["geno"]].append(v["methyl"])
    big = {k: v for k, v in pops.items() if len(v) >= MIN_POP}
    if len(big) < 2:
        return {"region": f"{chrom}:{lo}-{hi}", "n_sSNV": len(sites), "n_reads": len(reads),
                "n_pops_ge4": len(big), "delta_beta": None, "cn": g.get("cn"),
                "note": "only 1 population >=4 reads"}
    obs = max_pair_delta(big)
    # permutation null over reads in big populations only
    rm, rp = [], []
    for k, v in big.items():
        rm.extend(v)
        rp.extend([k] * len(v))
    pp = perm_p(rm, rp, obs, rng)
    # single-sSNV allele-Δβ (cis-ASM proxy)
    single_max = 0.0
    for p, (ref, alt) in sites:
        alt_m = [v["methyl"] for v in reads.values() if v["sites"].get(p) == "A"]
        ref_m = [v["methyl"] for v in reads.values() if v["sites"].get(p) == "R"]
        if len(alt_m) >= MIN_POP and len(ref_m) >= MIN_POP:
            single_max = max(single_max, abs(sum(alt_m) / len(alt_m) - sum(ref_m) / len(ref_m)))
    # ARI vs top-2 populations
    order = sorted(big, key=lambda k: -len(big[k]))
    top2 = set(order[:2])
    v2 = [(v["methyl"], v["geno"]) for v in reads.values() if v["geno"] in top2]
    ari = ari_2grp([x[0] for x in v2], [x[1] for x in v2]) if len(v2) >= 8 else None
    means = {k: round(sum(v) / len(v), 3) for k, v in big.items()}
    return {"region": f"{chrom}:{lo}-{hi}", "n_sSNV": len(sites), "n_reads": len(reads),
            "n_pops_ge4": len(big), "pop_sizes": {k: len(v) for k, v in sorted(big.items(), key=lambda x: -len(x[1]))},
            "pop_mean_methyl": {k: means[k] for k in order},
            "delta_beta": round(obs, 3), "perm_p": pp,
            "single_sSNV_max_db": round(single_max, 3),
            "lineage_beyond_cis": round(obs - single_max, 3),
            "ari_2grp": ari, "cn": g.get("cn")}


def bh_fdr(pvals):
    idx = sorted(range(len(pvals)), key=lambda i: pvals[i])
    n = len(pvals)
    q = [1.0] * n
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = idx[rank]
        prev = min(prev, pvals[i] * n / (rank + 1))
        q[i] = round(prev, 4)
    return q


def main():
    sample, msroot, bam_path = sys.argv[1], sys.argv[2], sys.argv[3]
    max_regions = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    import random
    rng = random.Random(seed)
    cen = load_census(msroot)
    regs = load_regions(msroot, max_regions)
    bam = pysam.AlignmentFile(bam_path, "rb")
    out = []
    for g in regs:
        r = analyze_region(bam, cen, g, rng)
        if r:
            out.append(r)
    bam.close()
    comp = [r for r in out if r.get("delta_beta") is not None]
    if comp:
        qs = bh_fdr([r["perm_p"] for r in comp])
        for r, q in zip(comp, qs):
            r["perm_q_fdr"] = q
    outp = os.path.join("/big7_disk/liaoyoyo2001/cnv_sv_work/methyl_lineage",
                        f"methyl_vs_ssnv_lineage_{sample}.json")
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    n_sig = sum(1 for r in comp if r.get("perm_q_fdr", 1) < 0.05)
    n_beyond = sum(1 for r in comp if r["lineage_beyond_cis"] > 0.05)
    med_db = sorted([r["delta_beta"] for r in comp])[len(comp) // 2] if comp else None
    json.dump({"sample": sample, "n_regions": len(comp),
               "n_sig_permFDR": n_sig, "n_lineage_beyond_cis": n_beyond,
               "median_delta_beta": med_db, "regions": comp}, open(outp, "w"),
              ensure_ascii=False, indent=2)
    print(f"{sample}: regions={len(comp)} sig(FDR<.05)={n_sig} "
          f"lineage>cis={n_beyond} median_db={med_db} -> {outp}")


if __name__ == "__main__":
    main()
