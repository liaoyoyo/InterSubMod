#!/bin/bash
# Full HCC1395 TP (non-LOH + LOH) dual-axis ASM survey + optional FP set.
# Extends 16_genome_asm_batch.sh:
#   - NO LOH subtract (run ALL TP; loh_status annotated by 18_dual_axis_pivot.py)
#   - dual-axis pivot (HP-axis + ALLELE-axis)
#   - VARSET param: "tp" (SEQC2 truth) or "fp" (caller FP)
# Disk-safe: per-chr Level1 deleted after pivot.
set -uo pipefail

export TMPDIR=/big7_disk/liaoyoyo2001/tmp
MSA=/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/build/bin/msa
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
TUMOR=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam
NORMAL=/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam
SEQC2_TP=/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
FP_VCF=/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/filtered_snv_fp.vcf.gz
LOH=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/output/seqc2_loh_only.bed
PIVOT=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/scripts/18_dual_axis_pivot.py
CLAIRS=/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/filtered_snv_tp.vcf

VARSET="${1:-tp}"          # tp | fp
CHROMS="${2:-chr13}"       # default smoke chr13; pass full list for genome-wide

WORK=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2
mkdir -p "$WORK/vcf" "$WORK/msa_tmp" "$WORK/logs"
MASTER="$WORK/asm_dualaxis_${VARSET}.tsv"

if [ "$VARSET" = "tp" ]; then SRC="$SEQC2_TP"; else SRC="$FP_VCF"; fi

CLAIRS_HDR="$WORK/clairs_header.txt"
grep "^##" "$CLAIRS" > "$CLAIRS_HDR"

echo "[$(date +%H:%M)] DUAL-AXIS ASM survey START. VARSET=$VARSET Chroms=$CHROMS master=$MASTER"

for CHR in $CHROMS; do
  echo "[$(date +%H:%M)] === $CHR ($VARSET) ==="
  VCF="$WORK/vcf/${CHR}_${VARSET}_msa.vcf"
  cp "$CLAIRS_HDR" "$VCF"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n' >> "$VCF"
  # ALL TP/FP for this chr (NO LOH subtract — loh_status annotated downstream)
  bcftools view -f PASS -r "$CHR" "$SRC" 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} length($3)==1 && length($4)==1 {print $1,$2,".",$3,$4,"30","PASS","SOMATIC","GT:GQ:DP:AF:AD","0/1:30:50:0.2:25,10"}' \
    >> "$VCF"
  NV=$(grep -vc "^#" "$VCF")
  echo "[$(date +%H:%M)] $CHR $VARSET variants (all, SNV): $NV"
  [ "$NV" -eq 0 ] && { echo "  skip"; continue; }
  bgzip -f "$VCF" && tabix -f -p vcf "$VCF.gz"

  OUT="$WORK/msa_tmp/${CHR}_${VARSET}"
  rm -rf "$OUT"; mkdir -p "$OUT"
  "$MSA" -v "$VCF.gz" -r "$REF" -t "$TUMOR" -n "$NORMAL" \
    --window 1000 --meth-high 0.8 --meth-low 0.2 --min-strand-reads 1 \
    --max-read-depth 500 --threads 16 -o "$OUT" \
    > "$WORK/logs/${CHR}_${VARSET}_msa.log" 2>&1
  RC=$?
  [ $RC -ne 0 ] && { echo "  MSA exit $RC (see logs/${CHR}_${VARSET}_msa.log)"; continue; }

  L1=$(find "$OUT" -name "level1_raw_methylation_details.tsv.gz" | head -1)
  if [ -n "$L1" ] && [ -f "$L1" ]; then
    echo "  Level1 size: $(du -h "$L1" | cut -f1)"
    python3 "$PIVOT" "$L1" "$MASTER" "$LOH" 0.5 2>> "$WORK/logs/${CHR}_${VARSET}_pivot.log"
    rm -f "$L1"
  else
    echo "  no Level1 output"
  fi
  echo "[$(date +%H:%M)] $CHR $VARSET done. master rows: $(wc -l < "$MASTER" 2>/dev/null || echo 0)"
  df -h /big7_disk | tail -1 | awk '{print "  big7 used:",$5}'
done

echo "[$(date +%H:%M)] COMPLETE $VARSET. Master: $MASTER ($(($(wc -l < "$MASTER" 2>/dev/null || echo 1)-1)) records)"
