#!/bin/bash
# region_samtools_snapshot.sh - Collect lightweight samtools evidence for a target region.

set -euo pipefail

SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
BAM=""
REGION=""
OUTPUT_DIR=""
MAX_READS=80

show_help() {
    cat <<'EOF'
Usage: region_samtools_snapshot.sh --bam FILE --region chr:start-end --output-dir DIR [--max-reads N]

Outputs:
  - summary.txt
  - depth.tsv
  - mpileup.txt
  - reads.sam
  - hp_tag_counts.tsv
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam) BAM="$2"; shift 2 ;;
        --region) REGION="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --max-reads) MAX_READS="$2"; shift 2 ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${SAMTOOLS_BIN}" ]]; then
    echo "[ERROR] samtools not found" >&2
    exit 1
fi
if [[ -z "${BAM}" || -z "${REGION}" || -z "${OUTPUT_DIR}" ]]; then
    show_help >&2
    exit 1
fi
if [[ ! -f "${BAM}" ]]; then
    echo "[ERROR] BAM not found: ${BAM}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

READ_COUNT="$("${SAMTOOLS_BIN}" view -c "${BAM}" "${REGION}")"
{
    echo "bam=${BAM}"
    echo "region=${REGION}"
    echo "read_count=${READ_COUNT}"
    echo "max_reads=${MAX_READS}"
} > "${OUTPUT_DIR}/summary.txt"

"${SAMTOOLS_BIN}" depth -aa -r "${REGION}" "${BAM}" > "${OUTPUT_DIR}/depth.tsv"
"${SAMTOOLS_BIN}" mpileup -r "${REGION}" "${BAM}" > "${OUTPUT_DIR}/mpileup.txt"
"${SAMTOOLS_BIN}" view "${BAM}" "${REGION}" | awk -v limit="${MAX_READS}" 'NR <= limit { print }' > "${OUTPUT_DIR}/reads.sam"

awk '
{
    hp="NA";
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^HP:Z:/) {
            hp = substr($i, 6);
            break;
        } else if ($i ~ /^HP:i:/) {
            hp = substr($i, 6);
            break;
        }
    }
    counts[hp]++;
}
END {
    print "hp_tag\tcount";
    for (k in counts) {
        print k "\t" counts[k];
    }
}
' "${OUTPUT_DIR}/reads.sam" | {
    IFS= read -r header || true
    printf "%s\n" "${header}"
    sort
} > "${OUTPUT_DIR}/hp_tag_counts.tsv"

echo "[region_samtools_snapshot] Wrote ${OUTPUT_DIR}" >&2
