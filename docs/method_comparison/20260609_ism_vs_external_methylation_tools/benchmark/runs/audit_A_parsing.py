#!/usr/bin/env python3
"""Audit A: ISM/MSA Level-1 per-CpG 5mC beta (tumor HP1) vs modkit pileup per-CpG 5mC%
   -> validates ISM MM/ML parsing (the zero-unit-test high-risk area, audit T1).
   Also builds per-CpG counts (HP1 vs HP1-1 from modkit Nmod/Ncanon) + Fisher p
   for Audit B (DSS comparison)."""
import gzip, json, csv, glob, sys
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact, pearsonr, spearmanr

THR = 0.5
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod"
B = f"{ROOT}/docs/method_comparison/20260609_ism_vs_external_methylation_tools/benchmark"
LOCI = {
    "BRCA2":   ("chr13", glob.glob(f"{ROOT}/research/tsg_promoter_asm_reviewer/pipeline/cache/level1/chr13_32315128_*.tsv.gz")[0]),
    "TBC1D16": ("chr17", glob.glob(f"{ROOT}/research/tsg_promoter_asm_reviewer/pipeline/cache/level1/chr17_79991120_*.tsv.gz")[0]),
}

def ism_beta(level1, hap):
    pc = defaultdict(list)
    with gzip.open(level1, "rt") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row["bam_source"] != "tumor" or row["haplotype_tag"] != hap or row["mod_code"] != "m":
                continue
            pc[int(row["methyl_pos"])].append(float(row["meth_call"]))
    out = {}
    for p, calls in pc.items():
        n = len(calls); m = sum(1 for c in calls if c >= THR)
        out[p] = {"beta": m / n, "nmeth": m, "nunmeth": n - m, "n": n}
    return out

def modkit_cpg(bed_gz, chrom):
    out = {}
    with gzip.open(bed_gz, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] != chrom or p[3] != "m":
                continue
            out[int(p[1])] = {"pct": float(p[10]) / 100.0, "nmod": int(p[11]), "ncan": int(p[12])}
    return out

resultA = {}
for locus, (chrom, level1) in LOCI.items():
    ism_hp1 = ism_beta(level1, "1")
    ism_hp1m1 = ism_beta(level1, "1-1")
    mk_hp1 = modkit_cpg(f"{B}/runs/pileup_somatic/hp1.bedmethyl.gz", chrom)
    mk_hp1m1 = modkit_cpg(f"{B}/runs/pileup_somatic/hp2.bedmethyl.gz", chrom)
    # find best coord offset (ISM methyl_pos vs modkit 0-based start)
    best = None
    for off in (-1, 0, 1):
        common = [p for p in ism_hp1 if (p + off) in mk_hp1]
        if best is None or len(common) > best[1]:
            best = (off, len(common))
    off = best[0]
    pairs = [(ism_hp1[p]["beta"], mk_hp1[p + off]["pct"]) for p in ism_hp1 if (p + off) in mk_hp1]
    ism_b = np.array([a for a, _ in pairs]); mk_p = np.array([b for _, b in pairs])
    pr = float(pearsonr(ism_b, mk_p)[0]) if len(pairs) > 2 else None
    sr = float(spearmanr(ism_b, mk_p)[0]) if len(pairs) > 2 else None
    resultA[locus] = {
        "chrom": chrom, "coord_offset": off, "n_matched_cpg": len(pairs),
        "n_ism_cpg": len(ism_hp1), "n_modkit_cpg": len(mk_hp1),
        "pearson_r": pr, "spearman_r": sr,
        "mean_abs_diff": float(np.mean(np.abs(ism_b - mk_p))) if len(pairs) else None,
        "mean_beta_ism(thr0.5)": float(ism_b.mean()) if len(pairs) else None,
        "mean_pct_modkit(thr0.74)": float(mk_p.mean()) if len(pairs) else None,
        "note": "threshold differs (ISM 0.5 vs modkit auto ~0.74) -> systematic offset expected; r tests pattern/position fidelity",
    }
    # build counts for Audit B (HP1 vs HP1-1 from modkit, matched positions)
    with open(f"{B}/runs/counts_{locus}.tsv", "w") as out:
        out.write("chr\tpos\thp1_N\thp1_X\thp2_N\thp2_X\tfisher_p\n")
        for p in sorted(mk_hp1):
            if p not in mk_hp1m1:
                continue
            a = mk_hp1[p]; b = mk_hp1m1[p]
            n1 = a["nmod"] + a["ncan"]; n2 = b["nmod"] + b["ncan"]
            if n1 < 5 or n2 < 5:
                continue
            try:
                fp = float(fisher_exact([[a["nmod"], a["ncan"]], [b["nmod"], b["ncan"]]])[1])
            except Exception:
                fp = None
            out.write(f"{chrom}\t{p+1}\t{n1}\t{a['nmod']}\t{n2}\t{b['nmod']}\t{fp}\n")

json.dump(resultA, open(f"{B}/runs/audit_A.json", "w"), indent=2)
print(json.dumps(resultA, indent=2))
