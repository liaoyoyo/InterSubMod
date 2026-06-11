#!/usr/bin/env python3
"""Genome-wide ASM survey analysis: BRCA2 ranking + spatial clustering + TSG enrichment."""
import subprocess
from collections import defaultdict
import numpy as np

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
MASTER = f"{BASE}/genome_survey/genome_asm_summary.tsv"
TSG_BED = f"{BASE}/output/tsg_promoter_merged.bed"

rows = []
with open(MASTER) as f:
    header = f.readline().rstrip("\n").split("\t")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        try:
            rows.append({
                "chrom": p[0], "pos": int(p[1]), "axis": p[2],
                "n_cpg": int(p[3]), "germ": float(p[4]), "som": float(p[5]),
                "delta": float(p[6]), "max_abs": float(p[7]),
                "p": float(p[8]) if p[8] != "nan" else 1.0,
            })
        except (ValueError, IndexError):
            continue

print(f"Total variant-axis records: {len(rows)}")
n_var = len(set((r['chrom'], r['pos']) for r in rows))
print(f"Unique variants: {n_var}")

# === 1. BRCA2 ranking ===
by_p = sorted(rows, key=lambda r: r["p"])
brca2 = [(i+1, r) for i, r in enumerate(by_p) if r["chrom"]=="chr13" and r["pos"]==32315128]
print(f"\n=== BRCA2 chr13:32315128 全基因組排名 ===")
for rank, r in brca2:
    print(f"  by p-value: rank #{rank}/{len(rows)} (top {100*rank/len(rows):.3f}%) | {r['axis']} p={r['p']:.2e} delta={r['delta']:+.3f} n_cpg={r['n_cpg']}")

# === 2. Significant ASM ===
sig = [r for r in rows if r["p"] < 0.05 and abs(r["delta"]) >= 0.1]
sig_hypo = [r for r in sig if r["delta"] <= -0.1]
print(f"\n=== 顯著 ASM (p<0.05 & |delta|>=0.1) ===")
print(f"  total: {len(sig)} | hypomethylation: {len(sig_hypo)} | hypermethylation: {len(sig)-len(sig_hypo)}")

# === 3. Top hypomethylation (reviewer-relevant direction) ===
print(f"\n=== Top 15 hypomethylation (delta 最負, p<0.05) ===")
top_hypo = sorted([r for r in rows if r["p"]<0.05], key=lambda r: r["delta"])[:15]
print(f"  {'chrom':<6}{'pos':<11}{'axis':<14}{'nCpG':>5}{'germ':>7}{'som':>7}{'delta':>8}{'p':>11}")
for r in top_hypo:
    print(f"  {r['chrom']:<6}{r['pos']:<11}{r['axis']:<14}{r['n_cpg']:>5}{r['germ']:>7.3f}{r['som']:>7.3f}{r['delta']:>+8.3f}{r['p']:>11.1e}")

# === 4. Spatial clustering (1-D gap-based on significant ASM positions) ===
print(f"\n=== Spatial clustering of significant ASM (gap < 10kb) ===")
GAP = 10000
sig_pos = defaultdict(list)
for r in sig:
    sig_pos[r["chrom"]].append(r)
clusters = []
for chrom, items in sig_pos.items():
    items = sorted(items, key=lambda r: r["pos"])
    cur = [items[0]]
    for r in items[1:]:
        if r["pos"] - cur[-1]["pos"] <= GAP:
            cur.append(r)
        else:
            if len(cur) >= 2:
                clusters.append((chrom, cur))
            cur = [r]
    if len(cur) >= 2:
        clusters.append((chrom, cur))
clusters.sort(key=lambda c: -len(c[1]))
print(f"  Total spatial clusters (>=2 sig ASM within 10kb): {len(clusters)}")
print(f"  Top 10 clusters by size:")
for chrom, cl in clusters[:10]:
    pmin = min(c["pos"] for c in cl); pmax = max(c["pos"] for c in cl)
    md = np.mean([c["delta"] for c in cl])
    best_p = min(c["p"] for c in cl)
    direction = "HYPO" if md < 0 else "HYPER"
    print(f"    {chrom}:{pmin}-{pmax} ({pmax-pmin}bp) | {len(cl)} sites | mean_delta={md:+.3f} {direction} | best_p={best_p:.1e}")

# === 5. TSG promoter enrichment ===
print(f"\n=== TSG promoter 富集 ===")
# write sig positions to bed
sig_bed = f"{BASE}/genome_survey/sig_asm.bed"
all_bed = f"{BASE}/genome_survey/all_asm.bed"
with open(sig_bed, "w") as f:
    for r in sorted(set((r["chrom"], r["pos"]) for r in sig)):
        f.write(f"{r[0]}\t{r[1]-1}\t{r[1]}\n")
with open(all_bed, "w") as f:
    for r in sorted(set((r["chrom"], r["pos"]) for r in rows)):
        f.write(f"{r[0]}\t{r[1]-1}\t{r[1]}\n")
def count_in_tsg(bed):
    out = subprocess.run(f"bedtools intersect -a {bed} -b {TSG_BED} -u 2>/dev/null | wc -l",
                         shell=True, capture_output=True, text=True)
    return int(out.stdout.strip() or 0)
n_sig_tsg = count_in_tsg(sig_bed)
n_all_tsg = count_in_tsg(all_bed)
n_sig_uniq = sum(1 for _ in open(sig_bed))
n_all_uniq = sum(1 for _ in open(all_bed))
sig_rate = n_sig_tsg / n_sig_uniq if n_sig_uniq else 0
all_rate = n_all_tsg / n_all_uniq if n_all_uniq else 0
print(f"  sig ASM in TSG promoter: {n_sig_tsg}/{n_sig_uniq} ({100*sig_rate:.2f}%)")
print(f"  all ASM in TSG promoter: {n_all_tsg}/{n_all_uniq} ({100*all_rate:.2f}%)")
print(f"  enrichment: {sig_rate/all_rate:.2f}x" if all_rate>0 else "  enrichment: NA")

# Save top hits + clusters for report
with open(f"{BASE}/genome_survey/genome_asm_top_hits.tsv", "w") as f:
    f.write("rank\tchrom\tpos\taxis\tn_cpg\tgerm\tsom\tdelta\tp\n")
    for i, r in enumerate(by_p[:100], 1):
        f.write(f"{i}\t{r['chrom']}\t{r['pos']}\t{r['axis']}\t{r['n_cpg']}\t{r['germ']}\t{r['som']}\t{r['delta']}\t{r['p']:.3e}\n")
print(f"\nSaved top-100 -> genome_survey/genome_asm_top_hits.tsv")
