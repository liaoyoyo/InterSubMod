#!/usr/bin/env python3
"""
Per-variant ASM pivot from MSA Level 1 (the gap MSA C++ doesn't compute).
Reads level1_raw_methylation_details.tsv.gz, computes per-variant:
  HP1 (germline) vs HP1-1 (somatic) per-CpG Δβ + paired Wilcoxon (tumor only).
Also HP2 vs HP2-1 where present.
Output: per-variant ASM summary TSV (appendable to genome-wide master).

Usage: python3 15_genome_asm_pivot.py <level1.tsv.gz> <out_summary.tsv> [meth_thresh]
"""
import sys, gzip
from collections import defaultdict
import numpy as np
from scipy import stats

level1 = sys.argv[1]
out_tsv = sys.argv[2]
# meth_call is already a probability 0-1 in MSA Level1; threshold to binary
METH_THR = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5

# per (somatic_pos) -> per (methyl_pos) -> per HP group -> [meth_count, total]
data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, 0])))

opener = gzip.open if level1.endswith(".gz") else open
with opener(level1, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    try:
        i_spos = idx["somatic_pos"]; i_mpos = idx["methyl_pos"]
        i_hp = idx["haplotype_tag"]; i_meth = idx["meth_call"]
        i_bam = idx["bam_source_id"]; i_chrom = idx["chrom"]
    except KeyError as e:
        print(f"[ERROR] missing column {e} in {level1}", file=sys.stderr)
        sys.exit(1)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= i_meth:
            continue
        if p[i_bam] != "tumor":   # ASM on tumor haplotypes
            continue
        spos = p[i_spos]; mpos = p[i_mpos]; hp = p[i_hp]
        try:
            mc = float(p[i_meth])
        except ValueError:
            continue
        chrom = p[i_chrom]
        key = (chrom, spos)
        cell = data[key][mpos][hp]
        cell[1] += 1
        if mc >= METH_THR:
            cell[0] += 1

# Per-variant ASM stats
rows = []
for (chrom, spos), cpgs in data.items():
    for (g_germ, g_som) in [("1", "1-1"), ("2", "2-1")]:
        paired_g, paired_s = [], []
        for mpos, hpd in cpgs.items():
            cg = hpd.get(g_germ); cs = hpd.get(g_som)
            if cg and cs and cg[1] >= 3 and cs[1] >= 3:
                paired_g.append(cg[0] / cg[1])
                paired_s.append(cs[0] / cs[1])
        n = len(paired_g)
        if n >= 5:
            deltas = [s - gv for gv, s in zip(paired_g, paired_s)]
            mean_delta = float(np.mean(deltas))
            max_abs = float(np.max(np.abs(deltas)))
            try:
                w_p = float(stats.wilcoxon(paired_g, paired_s).pvalue)
            except Exception:
                w_p = float("nan")
            rows.append({
                "chrom": chrom, "somatic_pos": spos,
                "axis": f"HP{g_germ}_vs_HP{g_som}",
                "n_paired_cpg": n,
                "mean_germ_beta": round(float(np.mean(paired_g)), 4),
                "mean_som_beta": round(float(np.mean(paired_s)), 4),
                "mean_delta": round(mean_delta, 4),
                "max_abs_delta": round(max_abs, 4),
                "wilcoxon_p": w_p,
            })

# Append/write
import os
write_header = not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0
cols = ["chrom", "somatic_pos", "axis", "n_paired_cpg", "mean_germ_beta",
        "mean_som_beta", "mean_delta", "max_abs_delta", "wilcoxon_p"]
with open(out_tsv, "a") as f:
    if write_header:
        f.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda x: x["mean_delta"]):
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

print(f"[pivot] {level1}: {len(rows)} variant-axis ASM records appended to {out_tsv}")
