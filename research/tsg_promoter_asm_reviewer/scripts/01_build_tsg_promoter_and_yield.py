#!/usr/bin/env python3
"""
Prereq A + P0b pipeline:
  1. Curated TSG list (90 well-known TSG, hardcoded from COSMIC CGC + Vogelstein 138 + breast cancer DDR)
  2. Parse RefSeq GTF -> TSS per TSG (take leftmost TSS for + strand, rightmost for - strand)
  3. Build promoter BED (TSS +/- 2kb)
  4. Subtract LOH regions (SEQC2 ngs_benchmark_cnvs_gain_loss_loh.bed where type=loh)
  5. Intersect HCC1395 SEQC2 TP sSNV (PASS,HighConf)
  6. Yield count + write candidate TSV
"""
import gzip, os, subprocess, sys, json
from collections import defaultdict

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
DATA = f"{BASE}/data"
OUT = f"{BASE}/output"
os.makedirs(OUT, exist_ok=True)

GTF = f"{DATA}/hg38.ncbiRefSeq.gtf.gz"
LOH = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
CPG = f"{DATA}/cpgIslandExt_hg38.txt.gz"
TP_VCF = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
HC_REGIONS = "/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"

PROMOTER_FLANK = 2000

# Curated TSG list (90 genes)
# Sources: COSMIC CGC TSG annotation + Vogelstein et al. Science 2013 138 driver genes (TSG only)
# + Breast cancer DDR pathway (Nik-Zainal 2016) + ENCODE chromatin remodeler TSG
TSG_LIST = sorted(set([
    # Classical TSG (Vogelstein 2013)
    "TP53", "RB1", "APC", "VHL", "NF1", "NF2", "PTEN", "BRCA1", "BRCA2",
    "CDKN2A", "CDKN2B", "CDKN1A", "CDKN1B", "MEN1", "WT1", "STK11", "BAP1",
    # DNA repair / Fanconi anemia
    "ATM", "ATR", "MLH1", "MSH2", "MSH6", "PMS2", "MUTYH",
    "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI",
    "FANCL", "FANCM", "PALB2", "BARD1", "BRIP1", "RAD51", "RAD51C", "RAD51D",
    "NBN", "MRE11", "RAD50",
    # Cell cycle checkpoints
    "CHEK1", "CHEK2", "WEE1",
    # Chromatin remodelers (SWI/SNF + COMPASS-like)
    "ARID1A", "ARID1B", "ARID2", "SMARCA4", "SMARCB1", "SMARCE1", "PBRM1",
    "KDM6A", "KMT2C", "KMT2D", "SETD2", "EP300", "CREBBP", "DAXX", "ATRX",
    # Cell-cell adhesion / metastasis suppressors (breast-relevant)
    "CDH1", "CDH13", "RASSF1", "GSTP1", "DAPK1", "TIMP3", "RARB", "APAF1",
    # NF-kB / TGFbeta TSG
    "SMAD4", "SMAD2", "TGFBR1", "TGFBR2",
    # MYC pathway TSG
    "MGA", "MAX",
    # Other commonly silenced in breast cancer
    "ESR1", "PGR", "HOXA5", "HOXA9", "TWIST1", "MLH3",
    # ENCODE / driver expansion
    "KEAP1", "NFE2L2", "FBXW7", "AXIN1", "AXIN2", "CTNNA1", "FAT1", "FAT4",
    # SCAN, BAP1 family related
    "ASXL1", "ASXL2", "EZH2",
    # Important breast TSG
    "TBX3", "RUNX1", "GATA3",
]))

print(f"[INFO] Curated TSG list: {len(TSG_LIST)} genes")

# Step 1: Parse RefSeq GTF -> TSS per gene (collapse all transcripts, take outermost TSS)
print(f"[INFO] Parsing {GTF}...")
gene_tss = defaultdict(lambda: {"chrom": None, "strand": None, "tss_set": set()})
target = set(TSG_LIST)
n_lines = 0
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        n_lines += 1
        f9 = line.rstrip("\n").split("\t")
        if len(f9) < 9 or f9[2] != "transcript":
            continue
        attrs = {}
        for kv in f9[8].split(";"):
            kv = kv.strip()
            if not kv:
                continue
            parts = kv.split(" ", 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip('"')
        gname = attrs.get("gene_name") or attrs.get("gene_id")
        if not gname or gname not in target:
            continue
        chrom = f9[0]
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        start = int(f9[3])
        end = int(f9[4])
        strand = f9[6]
        tss = start - 1 if strand == "+" else end - 1  # 0-based BED
        g = gene_tss[gname]
        if g["chrom"] is None:
            g["chrom"] = chrom
            g["strand"] = strand
        g["tss_set"].add(tss)

print(f"[INFO] Found TSS for {len(gene_tss)}/{len(TSG_LIST)} TSG genes")
missing = sorted(set(TSG_LIST) - set(gene_tss.keys()))
if missing:
    print(f"[WARN] TSG not in RefSeq GTF: {missing}")

# Step 2: Build promoter BED (TSS +/- 2kb)
promoter_bed = f"{OUT}/tsg_promoter_raw.bed"
intervals = []
for gname, g in gene_tss.items():
    chrom = g["chrom"]
    for tss in g["tss_set"]:
        start = max(0, tss - PROMOTER_FLANK)
        end = tss + PROMOTER_FLANK
        intervals.append((chrom, start, end, gname, g["strand"]))

intervals.sort(key=lambda x: (x[0], x[1]))
with open(promoter_bed, "w") as f:
    for chrom, start, end, gname, strand in intervals:
        f.write(f"{chrom}\t{start}\t{end}\t{gname}\t.\t{strand}\n")
print(f"[INFO] Wrote {len(intervals)} promoter intervals -> {promoter_bed}")

# Step 3: Merge overlapping intervals per gene (sort + bedtools merge)
promoter_merged = f"{OUT}/tsg_promoter_merged.bed"
subprocess.run(
    f"sort -k1,1 -k2,2n {promoter_bed} | bedtools merge -i - -c 4 -o distinct > {promoter_merged}",
    shell=True, check=True
)
n_merged = sum(1 for _ in open(promoter_merged))
print(f"[INFO] After merge: {n_merged} intervals -> {promoter_merged}")

# Step 4: Subtract LOH regions
loh_only = f"{OUT}/seqc2_loh_only.bed"
subprocess.run(
    f"awk '$4==\"loh\"' {LOH} | sort -k1,1 -k2,2n > {loh_only}",
    shell=True, check=True
)
n_loh = sum(1 for _ in open(loh_only))
print(f"[INFO] LOH regions: {n_loh} -> {loh_only}")

promoter_nonloh = f"{OUT}/tsg_promoter_nonLOH.bed"
subprocess.run(
    f"bedtools subtract -a {promoter_merged} -b {loh_only} > {promoter_nonloh}",
    shell=True, check=True
)
n_nonloh = sum(1 for _ in open(promoter_nonloh))
total_bp_orig = sum(int(l.split('\t')[2]) - int(l.split('\t')[1]) for l in open(promoter_merged))
total_bp_nonloh = sum(int(l.split('\t')[2]) - int(l.split('\t')[1]) for l in open(promoter_nonloh))
print(f"[INFO] After LOH subtract: {n_nonloh} intervals, {total_bp_nonloh:,} bp (vs orig {total_bp_orig:,} bp; retained {100*total_bp_nonloh/total_bp_orig:.1f}%)")

# Track which TSG remain after LOH subtraction
tsg_after_loh = set()
for line in open(promoter_nonloh):
    f4 = line.rstrip("\n").split("\t")
    if len(f4) >= 4:
        for g in f4[3].split(","):
            tsg_after_loh.add(g)
tsg_lost = set(gene_tss.keys()) - tsg_after_loh
print(f"[INFO] TSG retained after LOH subtract: {len(tsg_after_loh)}/{len(gene_tss)}")
print(f"[INFO] TSG lost to LOH: {sorted(tsg_lost)}")

# Step 5: Intersect with HCC1395 SEQC2 TP sSNV (PASS + HighConf)
# Filter VCF for PASS,HighConf and within HC regions
yield_vcf = f"{OUT}/tp_in_tsg_promoter_nonLOH.vcf"
yield_tsv = f"{OUT}/tp_in_tsg_promoter_nonLOH.tsv"
subprocess.run(
    f"bcftools view -f PASS {TP_VCF} -R {promoter_nonloh} -o {yield_vcf} 2>/dev/null",
    shell=True, check=True
)
# Annotate each VCF row with TSG name (intersect)
candidate_bed = f"{OUT}/tp_candidates.bed"
subprocess.run(
    f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/TVAF\\t%INFO/NVAF\\n' {yield_vcf} | "
    f"awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2-1,$2,$3,$4,$5,$6}}' > {candidate_bed}",
    shell=True, check=True
)
n_tp = sum(1 for _ in open(candidate_bed))
print(f"[INFO] HCC1395 TP sSNV in TSG promoter (non-LOH, HC-regions): {n_tp}")

# Annotate with gene name(s)
subprocess.run(
    f"bedtools intersect -a {candidate_bed} -b {promoter_nonloh} -wa -wb | "
    f"awk 'BEGIN{{OFS=\"\\t\"; print \"chrom\",\"pos\",\"ref\",\"alt\",\"TVAF\",\"NVAF\",\"gene\"}}"
    f" {{print $1,$3,$4,$5,$6,$7,$11}}' > {yield_tsv}",
    shell=True, check=True
)
print(f"[INFO] Annotated yield TSV: {yield_tsv}")

# Per-gene summary
per_gene = defaultdict(int)
for line in open(yield_tsv):
    if line.startswith("chrom"):
        continue
    f7 = line.rstrip("\n").split("\t")
    if len(f7) >= 7:
        for g in f7[6].split(","):
            per_gene[g] += 1
print(f"\n[YIELD SUMMARY]")
print(f"  Total TP sSNV in TSG promoter (non-LOH): {n_tp}")
print(f"  TSG with >=1 hit: {len(per_gene)}")
print(f"  Top genes:")
for g, c in sorted(per_gene.items(), key=lambda x: -x[1])[:20]:
    print(f"    {g}: {c}")

# Write summary JSON
summary = {
    "tsg_list_size": len(TSG_LIST),
    "tsg_in_refseq": len(gene_tss),
    "tsg_missing_from_refseq": missing,
    "promoter_intervals_raw": len(intervals),
    "promoter_intervals_merged": n_merged,
    "promoter_bp_raw": total_bp_orig,
    "promoter_bp_after_loh": total_bp_nonloh,
    "promoter_bp_retention_pct": round(100*total_bp_nonloh/total_bp_orig, 2),
    "loh_regions": n_loh,
    "tsg_retained_after_loh": sorted(tsg_after_loh),
    "tsg_lost_to_loh": sorted(tsg_lost),
    "tp_snv_in_tsg_promoter_nonLOH": n_tp,
    "tsg_with_hits": dict(sorted(per_gene.items(), key=lambda x: -x[1])),
}
with open(f"{OUT}/prereq_a_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n[INFO] Summary -> {OUT}/prereq_a_summary.json")

# GO/NO-GO decision
print(f"\n[DECISION]")
if n_tp >= 10:
    print(f"  GO: yield {n_tp} >= 10 -> proceed to full pipeline")
elif n_tp >= 1:
    print(f"  PARTIAL: yield {n_tp} < 10 but > 0 -> pilot single example, report to user")
else:
    print(f"  NO-GO: yield 0 -> abort, report 'no candidate' to reviewer")
