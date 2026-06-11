#!/usr/bin/env python3
"""
38_survivor_permutation.py — per-locus residual-copy permutation for the 2 survivors
(chr17:79991120, chr18:11741161) + BRCA2 re-confirm. Same test as 37(3): within somatic
tag, observed ALT-vs-REF dbeta vs random-relabel null. A genuine allele-cis sits OUTSIDE
the null; a copy/within-tag-heterogeneity effect sits INSIDE.
"""
import sys, glob, json, random
import numpy as np
sys.path.insert(0, "pipeline")
from lib.cis_asm_core import load_level1_plus, collapse_modtype, MIN_N
random.seed(20260603); np.random.seed(20260603)
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
LOCI = [("chr17", "79991120", "1-1"), ("chr18", "11741161", "1-1"),
        ("chr20", "7482356", "1-1"), ("chr13", "32315128", "1-1")]  # survivors + a collapser + BRCA2


def perm_locus(chrom, pos, somtag, nperm=2000):
    cache = glob.glob(f"{ROOT}/pipeline/cache/level1/{chrom}_{pos}_*.tsv.gz")
    if not cache:
        # try alt somatic tag 2-1
        return None
    Dp = load_level1_plus(cache[0]); D = collapse_modtype(Dp)
    mat = {}; allele = {}
    for (bam, hp, al), cpgd in D[pos].items():
        if bam != "tumor" or hp != somtag:
            continue
        for cpg, rd in cpgd.items():
            for rid, v in rd.items():
                mat.setdefault(rid, {})[cpg] = v
                allele[rid] = al
    reads = [r for r in mat if len(mat[r]) >= 5]
    if not reads:
        return None

    def diff(labels):
        bycpg = {}
        for rid in reads:
            for cpg, v in mat[rid].items():
                bycpg.setdefault(cpg, {"alt": [], "ref": []})
                if labels[rid] in ("alt", "ref"):
                    bycpg[cpg][labels[rid]].append(1 if v >= 0.5 else 0)
        ds = [np.mean(d["alt"]) - np.mean(d["ref"]) for d in bycpg.values()
              if len(d["alt"]) >= MIN_N and len(d["ref"]) >= MIN_N]
        return np.mean(ds) if len(ds) >= 4 else None

    obs = diff(allele)
    if obs is None:
        return None
    alt_n = sum(1 for r in reads if allele[r] == "alt")
    ref_n = sum(1 for r in reads if allele[r] == "ref")
    labs = [allele[r] for r in reads]
    null = []
    for _ in range(nperm):
        sh = labs[:]; random.shuffle(sh)
        d = diff({reads[i]: sh[i] for i in range(len(reads))})
        if d is not None:
            null.append(d)
    null = np.array(null)
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)
    return {"locus": f"{chrom}:{pos}", "obs": round(float(obs), 4), "alt_n": alt_n, "ref_n": ref_n,
            "null_mean": round(float(null.mean()), 4), "null_sd": round(float(null.std()), 4),
            "ci_lo": round(float(np.percentile(null, 2.5)), 3), "ci_hi": round(float(np.percentile(null, 97.5)), 3),
            "p": round(float(p), 4), "genuine": bool(p < 0.05)}


out = []
print(f"{'locus':18s}{'obs':>9}{'alt/ref':>10}{'null_95%CI':>20}{'p':>8}  verdict")
for chrom, pos, tag in LOCI:
    r = perm_locus(chrom, pos, tag)
    if r is None:
        print(f"{chrom}:{pos:14s}  (uncomputable)")
        continue
    out.append(r)
    v = "GENUINE allele-cis" if r["genuine"] else "inside null (copy/noise)"
    print(f"{r['locus']:18s}{r['obs']:>9}{f'{r[chr(97)+chr(108)+chr(116)+chr(95)+chr(110)]}/{r[chr(114)+chr(101)+chr(102)+chr(95)+chr(110)]}':>10}"
          f"{f'[{r[chr(99)+chr(105)+chr(95)+chr(108)+chr(111)]},{r[chr(99)+chr(105)+chr(95)+chr(104)+chr(105)]}]':>20}{r['p']:>8}  {v}")
json.dump(out, open(f"{ROOT}/genome_survey_v2/survivor_permutation.json", "w"), indent=1)
print(f"\n-> {sum(r['genuine'] for r in out)}/{len(out)} genuine allele-cis (permutation p<0.05)")
