#!/bin/bash
# Build subset VCFs for R1 HCC1395 Normal BAM pilot
# Subset VCF to the 860 TO / 983 paired filter SNV loci using bcftools

set -euo pipefail

DATA_DIR="/big7_disk/liaoyoyo2001/InterSubMod/research/F_hpfinengroups_deepening/data"
OUT_DIR="/big7_disk/liaoyoyo2001/InterSubMod/output/hcc1395_normal_pilot"
mkdir -p "$OUT_DIR"

# Source VCFs (TP + FP) — to be merged
TO_TP="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step04_benchmark_longphase_to/filtered_snv_tp.vcf.gz"
TO_FP="/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step04_benchmark_longphase_to/filtered_snv_fp.vcf.gz"

PAIRED_TP="/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395/20260211/longphase_s/filtered_snv_tp.vcf.gz"
PAIRED_FP="/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395/20260211/longphase_s/filtered_snv_fp.vcf.gz"

# Build region TSVs (chr\tpos) from filter SNV lists
build_region_tsv() {
    local src=$1 dst=$2
    awk 'NR>1 {print $1"\t"$2}' "$src" | sort -k1,1 -k2,2n | uniq > "$dst"
}

build_region_tsv "$DATA_DIR/r1_hcc1395_filter_snvs_to.tsv" "$OUT_DIR/to_regions.tsv"
build_region_tsv "$DATA_DIR/r1_hcc1395_filter_snvs_paired.tsv" "$OUT_DIR/paired_regions.tsv"

echo "[+] TO regions: $(wc -l < $OUT_DIR/to_regions.tsv)"
echo "[+] Paired regions: $(wc -l < $OUT_DIR/paired_regions.tsv)"

# Merge TP + FP per mode, subset to filter regions
merge_and_subset() {
    local mode=$1 tp=$2 fp=$3 regions=$4 out_vcf=$5
    echo "[+] Merging $mode TP+FP and subsetting to $regions ..."
    bcftools concat -a -d all "$tp" "$fp" 2>/dev/null | \
        bcftools sort -O z -o "$out_vcf.merged.vcf.gz" - 2>/dev/null
    bcftools index -t "$out_vcf.merged.vcf.gz"
    bcftools view -R "$regions" "$out_vcf.merged.vcf.gz" -O z -o "$out_vcf" 2>/dev/null
    bcftools index -t "$out_vcf"
    echo "    $mode subset: $(bcftools view -H $out_vcf | wc -l) variants"
}

merge_and_subset to "$TO_TP" "$TO_FP" "$OUT_DIR/to_regions.tsv" "$OUT_DIR/HCC1395_TO_filter_subset.vcf.gz"
merge_and_subset paired "$PAIRED_TP" "$PAIRED_FP" "$OUT_DIR/paired_regions.tsv" "$OUT_DIR/HCC1395_paired_filter_subset.vcf.gz"

echo "[+] Subset VCFs ready at $OUT_DIR/"
ls -la "$OUT_DIR/"*.vcf.gz
