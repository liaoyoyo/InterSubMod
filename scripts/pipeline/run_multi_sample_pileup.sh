#!/bin/bash
# run_multi_sample_pileup.sh - Run pileup benchmark for remaining samples
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}/../.."

# HCC1395 is running manually
SAMPLES=(COLO829 H1437 H2009 HCC1937 HCC1954 HCC1395_DORADO)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "================================================================="
    echo "  Processing ${SAMPLE} in PILEUP mode"
    echo "================================================================="
    bash "${ROOT_DIR}/scripts/pipeline/run_benchmark.sh" \
        --mode s-pure-pileup \
        --sample "${SAMPLE}" \
        --vcf-source pileup \
        --threads 112
    
    echo "Finished ${SAMPLE} pileup analysis."
    echo ""
done

echo "All pileup samples completed."
