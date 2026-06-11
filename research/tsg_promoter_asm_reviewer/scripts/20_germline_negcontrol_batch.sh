#!/bin/bash
# Negative control: germline het SNPs through the SAME dual-axis ASM pipeline.
# Biological null = baseline allelic methylation asymmetry at het sites (NOT somatic).
# If somatic strong-ASM rate >> germline het rate -> somatic events add ASM beyond baseline.
# Sample ~N het SNPs on a subset of chromosomes, exclude positions near TP/FP, run MSA dual-axis.
set -uo pipefail
export TMPDIR=/big7_disk/liaoyoyo2001/tmp
MSA=/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/build/bin/msa
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
TUMOR=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam
NORMAL=/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam
GERMLINE=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz
TP=/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
PIVOT=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/scripts/18_dual_axis_pivot.py
LOH=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/output/seqc2_loh_only.bed
CLAIRS=/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/filtered_snv_tp.vcf

WORK=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2
mkdir -p "$WORK/negctrl"
MASTER="$WORK/asm_dualaxis_germline_negctrl.tsv"
rm -f "$MASTER"
CLAIRS_HDR="$WORK/clairs_header.txt"

CHROMS="${1:-chr1 chr2 chr3 chr5 chr11 chr17}"   # subset for speed; ~ sample het across these
PER_CHR="${2:-60}"                                # het SNPs per chr (random)

echo "[$(date +%H:%M)] GERMLINE NEG-CONTROL START. chroms=$CHROMS per_chr=$PER_CHR"
for CHR in $CHROMS; do
  VCF="$WORK/negctrl/${CHR}_germ.vcf"
  cp "$CLAIRS_HDR" "$VCF"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n' >> "$VCF"
  # het SNPs (GT 0/1 or 0|1 or 1|0), SNV only, NOT within 2kb of any TP, random sample PER_CHR
  bcftools view -r "$CHR" "$GERMLINE" 2>/dev/null \
    | bcftools view -i 'GT="het"' 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} length($3)==1 && length($4)==1 {print $1,$2-1,$2,$3,$4}' \
    | bedtools window -w 2000 -a - -b <(bcftools view -r "$CHR" "$TP" 2>/dev/null | bcftools query -f '%CHROM\t%POS\t%POS\n' 2>/dev/null | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}') -v 2>/dev/null \
    | shuf | head -n "$PER_CHR" \
    | awk 'BEGIN{OFS="\t"}{print $1,$3,".",$4,$5,"30","PASS","GERMLINE","GT:GQ:DP:AF:AD","0/1:30:50:0.5:25,25"}' \
    | sort -k2,2n >> "$VCF"
  NV=$(grep -vc "^#" "$VCF")
  echo "[$(date +%H:%M)] $CHR germline het control variants: $NV"
  [ "$NV" -eq 0 ] && continue
  bgzip -f "$VCF" && tabix -f -p vcf "$VCF.gz"
  OUT="$WORK/negctrl/msa_$CHR"; rm -rf "$OUT"; mkdir -p "$OUT"
  "$MSA" -v "$VCF.gz" -r "$REF" -t "$TUMOR" -n "$NORMAL" \
    --window 1000 --meth-high 0.8 --meth-low 0.2 --min-strand-reads 1 \
    --max-read-depth 500 --threads 12 -o "$OUT" \
    > "$WORK/negctrl/${CHR}_msa.log" 2>&1
  L1=$(find "$OUT" -name "level1_raw_methylation_details.tsv.gz" | head -1)
  if [ -n "$L1" ] && [ -f "$L1" ]; then
    python3 "$PIVOT" "$L1" "$MASTER" "$LOH" 0.5 2>> "$WORK/negctrl/${CHR}_pivot.log"
    rm -f "$L1"
  fi
done
echo "[$(date +%H:%M)] NEG-CONTROL COMPLETE. Master: $MASTER ($(($(wc -l < "$MASTER" 2>/dev/null||echo 1)-1)) records)"
