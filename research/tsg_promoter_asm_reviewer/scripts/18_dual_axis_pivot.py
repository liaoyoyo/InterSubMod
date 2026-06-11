#!/usr/bin/env python3
"""
Dual-axis per-variant ASM pivot from MSA Level 1 (extends 15_genome_asm_pivot.py).

Computes TWO comparison axes per somatic variant (tumor reads only):
  1. HP-axis  : HP1 vs HP1-1, HP2 vs HP2-1  (haplotype germline vs somatic-reconstructed)
                — valid only where both haplotypes present (non-LOH).
  2. ALLELE-axis: ALT-allele reads vs REF-allele reads (somatic_allele_type)
                — valid even in LOH (one haplotype lost) since both ALT+REF reads
                  exist at VAF<1. Sidesteps the lost-haplotype problem.

Each somatic_pos is annotated with loh_status (passed in) so downstream can split.

Output columns (appended to master TSV):
  chrom somatic_pos axis axis_type loh_status n_paired_cpg mean_germ_beta
  mean_som_beta mean_delta max_abs_delta wilcoxon_p

  axis_type = HP | ALLELE
  axis      = HP1_vs_HP1-1 | HP2_vs_HP2-1 | ALT_vs_REF
  For ALLELE axis: "germ"=REF-allele reads, "som"=ALT-allele reads (delta=ALT-REF).

Usage: python3 18_dual_axis_pivot.py <level1.tsv.gz> <out_summary.tsv> <loh_bed> [meth_thresh]
"""
import sys, gzip, os, subprocess
from collections import defaultdict
import numpy as np
from scipy import stats

level1 = sys.argv[1]
out_tsv = sys.argv[2]
loh_bed = sys.argv[3]
METH_THR = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
MIN_N = 3      # per-CpG per-group min reads
MIN_PAIRED = 5 # min paired CpGs for a test

# Load LOH regions for annotation (chrom -> list of (start,end))
loh = defaultdict(list)
if os.path.exists(loh_bed):
    with open(loh_bed) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                try:
                    loh[p[0]].append((int(p[1]), int(p[2])))
                except ValueError:
                    continue
def in_loh(chrom, pos):
    pos = int(pos)
    for s, e in loh.get(chrom, []):
        if s <= pos <= e:
            return "LOH"
    return "nonLOH"

# per (chrom,spos) -> per methyl_pos -> per group -> {read_id: max_meth_call}
#   HP groups: '1','2','1-1','2-1'   ALLELE groups: 'ALT','REF'
# CRITICAL (2026-05-30 fix): BAM is 5mCG_5hmCG -> Level1 emits TWO rows per (read,CpG):
#   P(5mC) and P(5hmC). We collapse per-read with MAX across mod rows -> single
#   "any-modification" call per read (P_mod>=THR). This avoids the prior bug where each
#   mod row was counted independently (5hmC low-prob row inflated denominator -> beta halved).
data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

opener = gzip.open if level1.endswith(".gz") else open
with opener(level1, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    try:
        i_spos = idx["somatic_pos"]; i_mpos = idx["methyl_pos"]
        i_hp = idx["haplotype_tag"]; i_meth = idx["meth_call"]
        i_bam = idx["bam_source_id"]; i_chrom = idx["chrom"]
        i_allele = idx["somatic_allele_type"]; i_read = idx["read_id"]
    except KeyError as e:
        print(f"[ERROR] missing column {e} in {level1}", file=sys.stderr)
        sys.exit(1)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(i_meth, i_allele, i_read):
            continue
        if p[i_bam] != "tumor":
            continue
        try:
            mc = float(p[i_meth])
        except ValueError:
            continue
        key = (p[i_chrom], p[i_spos]); mpos = p[i_mpos]; rid = p[i_read]
        # HP group: collapse 5mC/5hmC rows per read via MAX
        cell = data[key][mpos][p[i_hp]]
        if rid not in cell or mc > cell[rid]:
            cell[rid] = mc
        # ALLELE group (alt/ref)
        at = p[i_allele].strip().lower()
        if at in ("alt", "ref"):
            c2 = data[key][mpos]["ALT" if at == "alt" else "REF"]
            if rid not in c2 or mc > c2[rid]:
                c2[rid] = mc

def beta_of(cell):
    """cell = {read_id: max_meth_call}; beta = frac reads with P_mod>=THR."""
    n = len(cell)
    if n == 0:
        return None, 0
    m = sum(1 for v in cell.values() if v >= METH_THR)
    return m / n, n

def paired_test(cpgs, g_germ, g_som):
    pg, ps = [], []
    for mpos, hpd in cpgs.items():
        bg, ng = beta_of(hpd.get(g_germ, {}))
        bs, ns = beta_of(hpd.get(g_som, {}))
        if bg is not None and bs is not None and ng >= MIN_N and ns >= MIN_N:
            pg.append(bg); ps.append(bs)
    n = len(pg)
    if n < MIN_PAIRED:
        return None
    deltas = [s - gv for gv, s in zip(pg, ps)]
    try:
        w_p = float(stats.wilcoxon(pg, ps).pvalue)
    except Exception:
        w_p = float("nan")
    return {
        "n_paired_cpg": n,
        "mean_germ_beta": round(float(np.mean(pg)), 4),
        "mean_som_beta": round(float(np.mean(ps)), 4),
        "mean_delta": round(float(np.mean(deltas)), 4),
        "max_abs_delta": round(float(np.max(np.abs(deltas))), 4),
        "wilcoxon_p": w_p,
    }

rows = []
for (chrom, spos), cpgs in data.items():
    loh_status = in_loh(chrom, spos)
    # HP-axis
    for (g_germ, g_som, axis) in [("1", "1-1", "HP1_vs_HP1-1"), ("2", "2-1", "HP2_vs_HP2-1")]:
        r = paired_test(cpgs, g_germ, g_som)
        if r:
            r.update({"chrom": chrom, "somatic_pos": spos, "axis": axis,
                      "axis_type": "HP", "loh_status": loh_status})
            rows.append(r)
    # ALLELE-axis (germ=REF, som=ALT)
    r = paired_test(cpgs, "REF", "ALT")
    if r:
        r.update({"chrom": chrom, "somatic_pos": spos, "axis": "ALT_vs_REF",
                  "axis_type": "ALLELE", "loh_status": loh_status})
        rows.append(r)

cols = ["chrom", "somatic_pos", "axis", "axis_type", "loh_status",
        "n_paired_cpg", "mean_germ_beta", "mean_som_beta", "mean_delta",
        "max_abs_delta", "wilcoxon_p"]
write_header = not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0
with open(out_tsv, "a") as f:
    if write_header:
        f.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda x: x["mean_delta"]):
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

n_hp = sum(1 for r in rows if r["axis_type"] == "HP")
n_al = sum(1 for r in rows if r["axis_type"] == "ALLELE")
print(f"[dual-pivot] {os.path.basename(level1)}: {len(rows)} rows ({n_hp} HP-axis, {n_al} ALLELE-axis) -> {out_tsv}")
