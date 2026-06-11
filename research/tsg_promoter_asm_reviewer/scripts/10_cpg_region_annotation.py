#!/usr/bin/env python3
"""
Annotate each differential CpG with functional region + assess BRCA2 promoter impact.

Functional regions (hg38):
  BRCA2 MANE TSS (ENST00000380152):     32,315,508 (+)
  BRCA2 alt-upstream TSS (ENST544455):  32,315,086 (+)
  ZAR1L canonical TSS (ENST533490):     32,315,363 (- strand)
  CpG island "CpG: 43":                 32,315,396 - 32,315,763
  Variant:                              32,315,128
  BRCA2 CDS start (ATG):                32,316,461
  Liu 2010 minimal bidirectional prom:  BRCA2-TSS -187 to +310 (use alt TSS as Liu's reference)
  ZAR2 exon 1 (Liu 2010):               BRCA2-TSS -134 to +111
"""
import json
import numpy as np
from scipy import stats

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"

# Functional landmarks (hg38)
BRCA2_MANE_TSS = 32315508
BRCA2_ALT_TSS = 32315086
ZAR1L_TSS = 32315363       # - strand
CPG_ISLAND = (32315396, 32315763)
VARIANT = 32315128
BRCA2_CDS = 32316461
# Liu 2010 uses a BRCA2 TSS reference; their +1 aligns near alt-upstream isoform
LIU_REF_TSS = BRCA2_ALT_TSS  # 32315086
LIU_PROMOTER = (LIU_REF_TSS - 187, LIU_REF_TSS + 310)   # 32314899 - 32315396
ZAR2_EXON1 = (LIU_REF_TSS - 134, LIU_REF_TSS + 111)      # 32314952 - 32315197

with open(f"{BASE}/output/step4_ism_results.json") as f:
    results = json.load(f)
brca2 = next(r for r in results if r["candidate"]["gene"] == "BRCA2")
cpgs = brca2["cpg_records"]

def annotate_region(pos):
    """Return list of functional regions this CpG falls into."""
    regions = []
    # BRCA2 relative position (MANE)
    d_mane = pos - BRCA2_MANE_TSS
    if d_mane < 0:
        regions.append(f"BRCA2_upstream_promoter(MANE{d_mane:+d})")
    elif pos < BRCA2_CDS:
        regions.append(f"BRCA2_5UTR(MANE{d_mane:+d})")
    # ZAR1L (- strand): TSS at 32315363, gene body extends to smaller coords
    if pos <= ZAR1L_TSS and pos >= 32303699:
        d_zar = ZAR1L_TSS - pos  # downstream into ZAR1L
        regions.append(f"ZAR1L_genebody(+{d_zar}from_TSS)")
    # CpG island
    if CPG_ISLAND[0] <= pos <= CPG_ISLAND[1]:
        regions.append("CpG_island")
    # Liu 2010 bidirectional promoter
    if LIU_PROMOTER[0] <= pos <= LIU_PROMOTER[1]:
        regions.append("Liu2010_bidirectional_promoter")
    # ZAR2 exon 1
    if ZAR2_EXON1[0] <= pos <= ZAR2_EXON1[1]:
        regions.append("ZAR2_exon1")
    return regions

# Build annotated table
print(f"{'CpG_pos':<12} {'dist_var':>8} {'HP1_b':>7} {'HP1-1_b':>8} {'Delta':>7} {'n_HP1':>5} {'n_HP1-1':>7}  regions")
print("-" * 130)

annotated = []
for r in cpgs:
    if r["HP1_n"] >= 3 and r["HP1-1_n"] >= 3 and r["HP1_beta"] is not None and r["HP1-1_beta"] is not None:
        pos = VARIANT + r["dist_to_var"]
        delta = r["HP1-1_beta"] - r["HP1_beta"]
        regions = annotate_region(pos)
        annotated.append({
            "pos": pos, "dist_var": r["dist_to_var"],
            "hp1_beta": r["HP1_beta"], "hp1_1_beta": r["HP1-1_beta"],
            "delta": delta, "n_hp1": r["HP1_n"], "n_hp1_1": r["HP1-1_n"],
            "regions": regions,
        })

# Sort by absolute delta (most differential first)
annotated.sort(key=lambda x: x["delta"])  # most negative (hypomethylated) first

print("\n### TOP 15 most HYPOMETHYLATED CpG (somatic < germline) ###")
for a in annotated[:15]:
    reg = ",".join(a["regions"]) if a["regions"] else "intergenic/downstream"
    print(f"chr13:{a['pos']:<8} {a['dist_var']:>+6}bp {a['hp1_beta']:>6.3f} {a['hp1_1_beta']:>7.3f} {a['delta']:>+7.3f} {a['n_hp1']:>5} {a['n_hp1_1']:>7}  {reg}")

# Region-level aggregation
print("\n### REGION-LEVEL hypomethylation summary ###")
region_stats = {}
for a in annotated:
    for reg in (a["regions"] or ["intergenic/downstream"]):
        # collapse the parametrized region names
        key = reg.split("(")[0]
        region_stats.setdefault(key, []).append(a["delta"])

print(f"{'Region':<40} {'n_CpG':>6} {'mean_Delta':>11} {'min_Delta':>10} {'n_hypo(<-0.1)':>14}")
print("-" * 90)
for key in sorted(region_stats, key=lambda k: np.mean(region_stats[k])):
    d = region_stats[key]
    n_hypo = sum(1 for x in d if x < -0.1)
    print(f"{key:<40} {len(d):>6} {np.mean(d):>+11.4f} {min(d):>+10.3f} {n_hypo:>14}")

# Statistical test per key region: is mean delta significantly < 0?
print("\n### Per-region one-sample Wilcoxon (H0: median Delta = 0) ###")
for key in sorted(region_stats, key=lambda k: np.mean(region_stats[k])):
    d = region_stats[key]
    if len(d) >= 6:
        try:
            stat, p = stats.wilcoxon(d)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"  {key:<40} n={len(d):>3}  mean_Delta={np.mean(d):>+.4f}  Wilcoxon_p={p:.2e} {sig}")
        except Exception as e:
            print(f"  {key:<40} n={len(d):>3}  (test failed: {e})")
    else:
        print(f"  {key:<40} n={len(d):>3}  (n<6, skip test)")

# Key question: do the differential CpGs concentrate near BRCA2 promoter / CpG island?
print("\n### BRCA2 promoter impact assessment ###")
# CpGs within BRCA2 proximal promoter (MANE TSS -250 to +250)
prox = [a for a in annotated if -250 <= (a["pos"] - BRCA2_MANE_TSS) <= 250]
cpgi = [a for a in annotated if CPG_ISLAND[0] <= a["pos"] <= CPG_ISLAND[1]]
liu = [a for a in annotated if LIU_PROMOTER[0] <= a["pos"] <= LIU_PROMOTER[1]]
print(f"  CpGs in BRCA2 proximal promoter (MANE TSS +/-250bp): {len(prox)}, mean Delta = {np.mean([a['delta'] for a in prox]):+.4f}" if prox else "  (none in proximal)")
print(f"  CpGs in CpG island ({CPG_ISLAND[0]}-{CPG_ISLAND[1]}): {len(cpgi)}, mean Delta = {np.mean([a['delta'] for a in cpgi]):+.4f}" if cpgi else "  (none in CpG island)")
print(f"  CpGs in Liu2010 bidirectional promoter: {len(liu)}, mean Delta = {np.mean([a['delta'] for a in liu]):+.4f}" if liu else "  (none in Liu promoter)")

# Save annotated table
import csv
with open(f"{BASE}/output/brca2_cpg_region_annotated.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["chrom", "pos", "dist_to_variant", "HP1_germline_beta", "HP1-1_somatic_beta", "delta_beta", "n_HP1", "n_HP1-1", "functional_regions"])
    for a in sorted(annotated, key=lambda x: x["pos"]):
        w.writerow(["chr13", a["pos"], a["dist_var"], f"{a['hp1_beta']:.4f}", f"{a['hp1_1_beta']:.4f}",
                    f"{a['delta']:+.4f}", a["n_hp1"], a["n_hp1_1"], ";".join(a["regions"]) if a["regions"] else "intergenic"])
print(f"\nSaved annotated TSV -> output/brca2_cpg_region_annotated.tsv ({len(annotated)} CpGs)")
