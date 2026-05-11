#!/bin/bash
# Per-chr HP tag aggregation across 4 BAM versions.
# baseline/V5/V6 use HP:i:integer (1, 11, 2, 21, 33)
# paired_T uses HP:Z:string (1, 1-1, 2, 2-1, 3)
#
# Uses samtools view -c with awk-style filter to count HP tags per chr.
# Parallel across (version, chr) = 32 jobs total (8 chr x 4 versions).

set -euo pipefail

SAMTOOLS=/usr/local/bin/samtools
OUT_DIR=/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence
mkdir -p "$OUT_DIR/per_chr_hp"

declare -A BAMS=(
    [baseline]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam
    [V5]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam
    [V6]=/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam
    [paired_T]=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam
)

# Representative chrs: cancer-relevant (1, 3, 5, 8, 11, 17, 19) + X
CHRS=(chr1 chr3 chr5 chr8 chr11 chr17 chr19 chrX)

scan_one() {
    local version=$1
    local chr=$2
    local bam=$3
    local out_file="$OUT_DIR/per_chr_hp/${version}_${chr}.tsv"
    # Skip if already done
    if [[ -f "$out_file" ]]; then
        echo "  SKIP ${version} ${chr} (exists)"
        return
    fi

    # Filter: primary only (-F 256 secondary, -F 2048 supplementary, -F 4 unmapped)
    # Extract HP tag value, count per tag
    # HP:i: for integer, HP:Z: for string (paired_T)
    $SAMTOOLS view -F 2304 -F 4 "$bam" "$chr" \
        | awk 'BEGIN{OFS="\t"} {
            hp = "NA"
            for(i=12;i<=NF;i++){
                if($i ~ /^HP:[iZ]:/){
                    split($i, a, ":")
                    hp = a[3]
                    break
                }
            }
            print hp
        }' \
        | sort | uniq -c | sort -rn \
        | awk -v v="$version" -v c="$chr" 'BEGIN{OFS="\t"; print "version","chr","hp_tag","count"} {print v,c,$2,$1}' \
        > "$out_file"
    echo "  DONE ${version} ${chr} ($(wc -l < $out_file) tags)"
}

export -f scan_one
export SAMTOOLS OUT_DIR

echo "[scan-hp] $(date) — scanning per-chr HP tags"
echo "  BAMs: 4 versions x 8 chrs = 32 jobs, parallel"

# Build job list
JOB_ARGS=()
for version in baseline V5 V6 paired_T; do
    for chr in "${CHRS[@]}"; do
        bam=${BAMS[$version]}
        JOB_ARGS+=("$version|$chr|$bam")
    done
done

# Run with limited parallelism (8 concurrent to avoid IO contention)
printf '%s\n' "${JOB_ARGS[@]}" | xargs -I{} -P 8 bash -c '
    IFS="|" read -r v c b <<< "$1"
    scan_one "$v" "$c" "$b"
' _ {}

echo "[scan-hp] $(date) — DONE"

# Aggregate to single TSV
echo
echo "=== Aggregated per-chr HP distribution ==="
{
    echo -e "version\tchr\thp_tag\tcount"
    for f in "$OUT_DIR/per_chr_hp"/*.tsv; do
        tail -n +2 "$f"
    done
} > "$OUT_DIR/per_chr_hp_all.tsv"

wc -l "$OUT_DIR/per_chr_hp_all.tsv"
echo "  Output: $OUT_DIR/per_chr_hp_all.tsv"
