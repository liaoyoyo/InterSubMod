#!/bin/bash
# Script to run ALL SNV tests from a real VCF file
#
# This script supports multiple modes and automatic cluster heatmap generation:
#   1. baseline: Normal filtering with debug logging
#   2. all-with-w1000: All filters enabled, ±1000bp window
#   3. all-with-w2000: All filters enabled, ±2000bp window
#
# Usage:
#   ./run_full_vcf_test.sh                          # Default baseline mode
#   ./run_full_vcf_test.sh --mode baseline          # Baseline mode explicitly
#   ./run_full_vcf_test.sh --mode all-with-w1000    # ±1000bp window mode
#   ./run_full_vcf_test.sh -o /path/to/output       # Custom output directory
#   ./run_full_vcf_test.sh --no-plots               # Skip cluster heatmap generation

set -e

# Defaults
VCF_PATH="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
THREADS=64
TUMOR_BAM="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL_BAM="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF_FASTA="/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
MODE="baseline"
OUTPUT_DIR=""
METRICS="NHD"
GENERATE_PLOTS=true
PLOT_THREADS=64
date_str=$(date +%Y%m%d)

# Script directory (for finding Python tools)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
PLOT_SCRIPT="${PROJECT_ROOT}/tools/plot_cluster_heatmap.py"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v|--vcf)
            VCF_PATH="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -o|--out)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        --metrics)
            METRICS="$2"
            shift 2
            ;;
        --no-plots)
            GENERATE_PLOTS=false
            shift
            ;;
        --plot-threads)
            PLOT_THREADS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -v, --vcf PATH       Path to VCF file"
            echo "  -t, --threads N      Number of threads for C++ (default: 64)"
            echo "  -o, --out DIR        Output directory"
            echo "  -m, --mode MODE      Test mode: baseline, all-with-w1000, all-with-w2000"
            echo "  --metrics LIST       Space-separated list of metrics (e.g. 'NHD L1')"
            echo "  --no-plots           Skip cluster heatmap generation"
            echo "  --plot-threads N     Number of threads for Python plotting (default: 16)"
            echo "  -h, --help           Show this help message"
            echo ""
            echo "Modes:"
            echo "  baseline:        Normal filtering with debug logging"
            echo "                   Outputs to: output/\${date}_vcf_baseline"
            echo ""
            echo "  all-with-w1000:  All filters enabled, ±1000bp window, no-filter output"
            echo "                   Outputs to: output/\${date}_vcf_all_w1000"
            echo ""
            echo "  all-with-w2000:  All filters enabled, ±2000bp window, no-filter output"
            echo "                   Outputs to: output/\${date}_vcf_all_w2000"
            echo ""
            echo "Cluster Heatmap:"
            echo "  By default, cluster heatmaps are generated after C++ processing."
            echo "  Use --no-plots to skip this step."
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Set output directory based on mode if not specified
if [[ -z "$OUTPUT_DIR" ]]; then
    case $MODE in
        baseline)
            OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output/${date_str}_vcf_baseline"
            ;;
        all-with-w1000)
            OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output/${date_str}_vcf_all_w1000"
            ;;
        all-with-w2000)
            OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output/${date_str}_vcf_all_w2000"
            ;;
        *)
            echo "Unknown mode: $MODE"
            echo "Valid modes: baseline, all-with-w1000, all-with-w2000"
            exit 1
            ;;
    esac
fi

LOG_FILE="${OUTPUT_DIR}/full_execution_analysis.log"

echo "=== InterSubMod FULL VCF Test ==="
echo "Mode: ${MODE}"
echo "VCF: ${VCF_PATH}"
echo "Tumor BAM: ${TUMOR_BAM}"
echo "Normal BAM: ${NORMAL_BAM}"
echo "Reference: ${REF_FASTA}"
echo "Output Dir: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Metrics: ${METRICS}"
echo "Log File: ${LOG_FILE}"
echo "---------------------------------"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Step 1: Run Execution
echo ""
echo "[1] Running InterSubMod in ${MODE} mode..."
EXECUTABLE="/big8_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod"

if [ ! -f "${EXECUTABLE}" ]; then
    echo "Error: Executable not found at ${EXECUTABLE}"
    echo "Please build the project first with:"
    echo "  cd /big8_disk/liaoyoyo2001/InterSubMod/build && cmake .. && make -j"
    exit 1
fi

# Construct Metric Flags
METRIC_FLAGS=""
for m in ${METRICS}; do
    METRIC_FLAGS="${METRIC_FLAGS} --distance-metric ${m}"
done

# Construct Command based on mode
case $MODE in
    baseline)
        # Baseline mode: Normal filtering with debug logging
        # - Default window size (±1000bp)
        # - Default filtering parameters
        # - Debug mode enabled for logging filtered reads
        CMD="${EXECUTABLE} \
            --tumor-bam ${TUMOR_BAM} \
            --normal-bam ${NORMAL_BAM} \
            --reference ${REF_FASTA} \
            --vcf ${VCF_PATH} \
            --output-dir ${OUTPUT_DIR} \
            --threads ${THREADS} \
            --log-level debug \
            --output-filtered-reads \
            ${METRIC_FLAGS}"
        ;;
    
    all-with-w1000)
        # All-with-window mode: All filters + larger window + no-filter output
        # - Window size: ±1000bp
        # - All filters enabled (default)
        # - No-filter flag: output all reads without filtering
        # - This is useful for verification and comparison
        CMD="${EXECUTABLE} \
            --tumor-bam ${TUMOR_BAM} \
            --normal-bam ${NORMAL_BAM} \
            --reference ${REF_FASTA} \
            --vcf ${VCF_PATH} \
            --output-dir ${OUTPUT_DIR} \
            --threads ${THREADS} \
            --window-size 1000 \
            --log-level debug \
            --output-filtered-reads \
            --no-filter \
            ${METRIC_FLAGS}"
        ;;

    all-with-w2000)
        # All-with-window mode: All filters + larger window + no-filter output
        # - Window size: ±2000bp
        # - All filters enabled (default)
        # - No-filter flag: output all reads without filtering
        # - This is useful for verification and comparison
        CMD="${EXECUTABLE} \
            --tumor-bam ${TUMOR_BAM} \
            --normal-bam ${NORMAL_BAM} \
            --reference ${REF_FASTA} \
            --vcf ${VCF_PATH} \
            --output-dir ${OUTPUT_DIR} \
            --threads ${THREADS} \
            --window-size 2000 \
            --log-level debug \
            --output-filtered-reads \
            --no-filter \
            ${METRIC_FLAGS}"
        ;;
esac

echo "Command to run (for IDE Debugging):"
echo "${CMD}"
echo ""

# Use /usr/bin/time for resource monitoring
echo "    Starting execution with ${THREADS} threads..."
echo "    (This may take a few minutes...)"
{
    /usr/bin/time -v ${CMD}
} 2>&1 | tee "${LOG_FILE}"

echo ""
echo "=== Analysis Log Saved to: ${LOG_FILE} ==="

# Step 2: Output summary
echo ""
echo "[2] Output Summary:"
echo "    Output directory: ${OUTPUT_DIR}"

# Count output files
if [ -d "${OUTPUT_DIR}" ]; then
    NUM_REGIONS=$(find "${OUTPUT_DIR}" -name "metadata.txt" 2>/dev/null | wc -l)
    NUM_READS_FILES=$(find "${OUTPUT_DIR}" -name "reads.tsv" 2>/dev/null | wc -l)
    NUM_METHYL_FILES=$(find "${OUTPUT_DIR}" -name "methylation.csv" 2>/dev/null | wc -l)
    NUM_FILTERED=$(find "${OUTPUT_DIR}" -name "filtered_reads.tsv" 2>/dev/null | wc -l)
    NUM_FORWARD=$(find "${OUTPUT_DIR}" -name "methylation_forward.csv" 2>/dev/null | wc -l)
    NUM_REVERSE=$(find "${OUTPUT_DIR}" -name "methylation_reverse.csv" 2>/dev/null | wc -l)
    
    echo "    Regions processed: ${NUM_REGIONS}"
    echo "    Reads files: ${NUM_READS_FILES}"
    echo "    Methylation matrices: ${NUM_METHYL_FILES}"
    echo "    Strand-specific matrices:"
    echo "      - Forward (+): ${NUM_FORWARD}"
    echo "      - Reverse (-): ${NUM_REVERSE}"
    if [[ "$MODE" == "baseline" ]] || [[ "$MODE" == "all-with-window" ]]; then
        echo "    Filtered reads logs: ${NUM_FILTERED}"
    fi
fi

# Step 3: Generate cluster heatmaps (if enabled)
if [[ "${GENERATE_PLOTS}" == "true" ]]; then
    echo ""
    echo "[3] Generating Distance-based Cluster Heatmaps..."
    
    # Use the new distance heatmap script
    DISTANCE_PLOT_SCRIPT="${SCRIPT_DIR}/../tools/plot_distance_heatmap.py"
    
    # Check if Python script exists
    if [ ! -f "${DISTANCE_PLOT_SCRIPT}" ]; then
        echo "    Warning: Distance plot script not found at ${DISTANCE_PLOT_SCRIPT}"
        echo "    Skipping heatmap generation."
    else
        # Check Python dependencies
        if ! python3 -c "import matplotlib, seaborn, scipy, pandas, numpy" 2>/dev/null; then
            echo "    Warning: Python dependencies not available."
            echo "    Please install: pip install matplotlib seaborn scipy pandas numpy"
            echo "    Skipping heatmap generation."
        else
            # Determine distance metric (use first metric)
            FIRST_METRIC=$(echo "${METRICS}" | awk '{print $1}')
            
            echo "    Using ${PLOT_THREADS} threads for plotting..."
            echo "    Distance metric: ${FIRST_METRIC}"
            echo "    Generating Read × Read distance heatmaps with dendrograms..."
            
            PLOT_START=$(date +%s)
            
            python3 "${DISTANCE_PLOT_SCRIPT}" \
                --output-dir "${OUTPUT_DIR}" \
                --threads "${PLOT_THREADS}" \
                --metric "${FIRST_METRIC}" \
                --linkage average \
                --min-reads 10 \
                --format png \
                --dpi 150
            
            PLOT_END=$(date +%s)
            PLOT_ELAPSED=$((PLOT_END - PLOT_START))
            
            # Count generated plots
            NUM_PLOTS=$(find "${OUTPUT_DIR}" -name "distance_heatmap*.png" 2>/dev/null | wc -l)
            
            echo ""
            echo "    Distance heatmaps generated: ${NUM_PLOTS}"
            echo "    Plotting time: ${PLOT_ELAPSED} seconds"
        fi
    fi
else
    echo ""
    echo "[3] Cluster heatmap generation skipped (use without --no-plots to enable)"
fi

echo ""
echo "=== Test Complete (Mode: ${MODE}) ==="
echo "Done."
