#!/bin/bash
# =============================================================================
# Run batch VCF analysis inside Docker container
# =============================================================================
# Usage: ./docker/run-analysis.sh [options]
# =============================================================================

set -e

echo "=== InterSubMod Batch Analysis ==="
echo "=================================="

# Default options
MODE="${MODE:-all-with-w5000}"
THREADS="${THREADS:-8}"
PLOT_TYPE="${PLOT_TYPE:-all}"

echo "Mode: ${MODE}"
echo "Threads: ${THREADS}"
echo "Plot Type: ${PLOT_TYPE}"

# Ensure build exists
if [ ! -f "/workspace/build/bin/inter_sub_mod" ]; then
    echo "Error: Binary not found. Run build.sh first."
    exit 1
fi

# Run the batch analysis
cd /workspace
./scripts/run_batch_vcf_analysis.sh \
    --mode "${MODE}" \
    --threads "${THREADS}" \
    --plot-type "${PLOT_TYPE}"

echo ""
echo "=== Analysis Complete ==="
