#!/bin/bash
# R1-Global: Build global VCF + truth table for HCC1395 TO arm
# Merge TP (28,396) + FP (11,843) = 40,239 global variants from step04 benchmark

set -euo pipefail

OUT_DIR="/big7_disk/liaoyoyo2001/InterSubMod/output/hcc1395_normal_pilot_global"
DATA_DIR="/big7_disk/liaoyoyo2001/InterSubMod/research/F_hpfinengroups_deepening/data"

TO_TP="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step04_benchmark_longphase_to/filtered_snv_tp.vcf.gz"
TO_FP="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step04_benchmark_longphase_to/filtered_snv_fp.vcf.gz"

mkdir -p "$OUT_DIR"

# 1. Build truth TSV (Chr, Pos, truth_label) covering all 40,239 variants
TRUTH_TSV="$DATA_DIR/r1_hcc1395_filter_snvs_to_global.tsv"
echo -e "Chr\tPos\ttruth_label" > "$TRUTH_TSV"
zcat "$TO_TP" | grep -v "^#" | awk -v OFS='\t' '{print $1, $2, "TP"}' >> "$TRUTH_TSV"
zcat "$TO_FP" | grep -v "^#" | awk -v OFS='\t' '{print $1, $2, "FP"}' >> "$TRUTH_TSV"
echo "[+] Truth TSV rows: $(( $(wc -l < "$TRUTH_TSV") - 1 ))  ($TRUTH_TSV)"

# 2. Merge TP + FP → HCC1395_TO_global.vcf.gz
OUT_VCF="$OUT_DIR/HCC1395_TO_global.vcf.gz"
echo "[+] Merging TP + FP..."
bcftools concat -a -d all "$TO_TP" "$TO_FP" 2>/dev/null \
    | bcftools sort -O z -o "$OUT_VCF" - 2>/dev/null
bcftools index -t "$OUT_VCF"
echo "[+] Merged VCF: $(bcftools view -H "$OUT_VCF" | wc -l) variants at $OUT_VCF"

ls -la "$OUT_DIR"/*.vcf.gz*
