#!/bin/bash
# A0a — 全基因組 HP-tag 分布（真實 paired BAM，全染色體平行）
# 唯讀；輸出 per-chr TSV + 彙總。unphase = 無 HP:Z tag 的 primary read。
set -uo pipefail

BAM=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam
OUT=/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets
mkdir -p "$OUT/per_chr"

count_chr() {
  local chr=$1
  local bam=$2
  local out=$3
  samtools view -F 0x900 "$bam" "$chr" 2>/dev/null | awk -v CHR="$chr" '
    {
      n++
      if (match($0,/HP:Z:[^\t]+/)) { v=substr($0,RSTART+5,RLENGTH-5); hp[v]++ }
      else { notag++ }
    }
    END{
      # 固定欄位順序輸出
      printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
        CHR, n, notag,
        hp["1"], hp["2"], hp["3"], hp["1-1"], hp["2-1"],
        (n - notag - hp["1"] - hp["2"] - hp["3"] - hp["1-1"] - hp["2-1"])
    }' > "$out/per_chr/${chr}.tsv"
  echo "  done $chr"
}
export -f count_chr

CHRS=$(seq 1 22 | sed 's/^/chr/'); CHRS="$CHRS chrX chrY"
echo "[A0a] counting HP-tag distribution across $(echo $CHRS | wc -w) chromosomes (parallel 8)..."
echo "$CHRS" | tr ' ' '\n' | xargs -P 8 -I {} bash -c 'count_chr "$@"' _ {} "$BAM" "$OUT"

# 彙總
HDR="chrom\ttotal\tunphase\tHP1\tHP2\tHP3\tHP1_1\tHP2_1\tother"
echo -e "$HDR" > "$OUT/hp_distribution_by_chr.tsv"
cat "$OUT"/per_chr/chr*.tsv 2>/dev/null | sort -V >> "$OUT/hp_distribution_by_chr.tsv"

# genome total + 百分比
awk -F'\t' 'NR>1{for(i=2;i<=9;i++)s[i]+=$i}
  END{
    t=s[2]
    printf "\n=== GENOME-WIDE TOTAL ===\n"
    printf "total primary reads = %d\n", t
    printf "unphase (no HP) = %d (%.2f%%)\n", s[3], 100*s[3]/t
    printf "HP1   = %d (%.2f%%)\n", s[4], 100*s[4]/t
    printf "HP2   = %d (%.2f%%)\n", s[5], 100*s[5]/t
    printf "HP3   = %d (%.2f%%)\n", s[6], 100*s[6]/t
    printf "HP1-1 = %d (%.2f%%)\n", s[7], 100*s[7]/t
    printf "HP2-1 = %d (%.2f%%)\n", s[8], 100*s[8]/t
    printf "other = %d (%.2f%%)\n", s[9], 100*s[9]/t
    # 守恆檢查
    sum=s[3]+s[4]+s[5]+s[6]+s[7]+s[8]+s[9]
    printf "守恆檢查: 各類和=%d vs total=%d → %s\n", sum, t, (sum==t?"PASS":"FAIL")
  }' "$OUT/hp_distribution_by_chr.tsv" | tee "$OUT/hp_distribution_genome_total.txt"

echo "[A0a] DONE → $OUT/hp_distribution_by_chr.tsv"
