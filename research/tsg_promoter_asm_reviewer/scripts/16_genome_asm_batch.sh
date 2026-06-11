#!/bin/bash
# Genome-wide non-LOH TP ASM survey — per-chromosome batch.
# Per chr: extract non-LOH PASS TP -> cloned-header reformat -> MSA -> Python pivot -> delete Level1 raw.
# Disk-safe: peak = single chromosome's Level 1 (deleted after pivot; only summary kept).
set -uo pipefail

export TMPDIR=/big7_disk/liaoyoyo2001/tmp
MSA=/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/build/bin/msa
REF=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
TUMOR=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam
NORMAL=/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam
SEQC2_TP=/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
LOH=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/output/seqc2_loh_only.bed
CLAIRS=/big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/filtered_snv_tp.vcf
PIVOT=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/scripts/15_genome_asm_pivot.py

WORK=/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey
mkdir -p "$WORK/vcf" "$WORK/msa_tmp" "$WORK/logs"
MASTER="$WORK/genome_asm_summary.tsv"

CLAIRS_HDR="$WORK/clairs_header.txt"
grep "^##" "$CLAIRS" > "$CLAIRS_HDR"

CHROMS="${1:-chr13}"   # default test chr13; pass full list to run genome-wide

echo "[$(date +%H:%M)] ASM survey START. Chroms: $CHROMS  master=$MASTER"

for CHR in $CHROMS; do
  echo "[$(date +%H:%M)] === $CHR ==="
  VCF="$WORK/vcf/${CHR}_msa.vcf"
  # build cloned-header VCF: non-LOH PASS TP for this chr + legal ClairS-style FORMAT
  cp "$CLAIRS_HDR" "$VCF"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n' >> "$VCF"
  bcftools view -f PASS -r "$CHR" "$SEQC2_TP" 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' 2>/dev/null \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4}' \
    | bedtools subtract -a - -b "$LOH" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"}{print $1,$3,".",$4,$5,"30","PASS","SOMATIC","GT:GQ:DP:AF:AD","0/1:30:50:0.2:25,10"}' \
    >> "$VCF"
  NV=$(grep -vc "^#" "$VCF")
  echo "[$(date +%H:%M)] $CHR non-LOH TP variants: $NV"
  [ "$NV" -eq 0 ] && { echo "  skip"; continue; }
  bgzip -f "$VCF" && tabix -f -p vcf "$VCF.gz"

  OUT="$WORK/msa_tmp/$CHR"
  rm -rf "$OUT"; mkdir -p "$OUT"
  "$MSA" -v "$VCF.gz" -r "$REF" -t "$TUMOR" -n "$NORMAL" \
    --window 1000 --meth-high 0.8 --meth-low 0.2 --min-strand-reads 1 \
    --max-read-depth 500 --threads 16 -o "$OUT" \
    > "$WORK/logs/${CHR}_msa.log" 2>&1
  RC=$?
  [ $RC -ne 0 ] && { echo "  MSA exit $RC (see logs/${CHR}_msa.log)"; continue; }
  grep -E "保留|過濾" "$WORK/logs/${CHR}_msa.log" | tail -1

  L1=$(find "$OUT" -name "level1_raw_methylation_details.tsv.gz" | head -1)
  if [ -n "$L1" ] && [ -f "$L1" ]; then
    L1SZ=$(du -h "$L1" | cut -f1)
    echo "  Level1 size: $L1SZ"
    python3 "$PIVOT" "$L1" "$MASTER" 0.5 2>> "$WORK/logs/${CHR}_pivot.log"
    rm -f "$L1"   # disk safety: drop big raw, keep Level2/3 summaries
  else
    echo "  no Level1 output"
  fi
  echo "[$(date +%H:%M)] $CHR done. master rows: $(wc -l < "$MASTER" 2>/dev/null || echo 0)"
  df -h /big7_disk | tail -1 | awk '{print "  big7 used:",$5}'
done

echo "[$(date +%H:%M)] COMPLETE. Master: $MASTER ($(($(wc -l < "$MASTER")-1)) records)"
