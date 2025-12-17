#!/bin/bash
# Script to run ALL SNV tests from a real VCF file
#
# This script supports multiple modes and automatic cluster heatmap generation:
#   1. baseline: Normal filtering with debug logging
#   2. all-with-w1000: All filters enabled, ±1000bp window
#   3. all-with-w2000: All filters enabled, ±2000bp window
#
# Heatmap generation:
#   By default, both distance_heatmap and cluster_heatmap are generated.
#   Use --plot-type to specify which plot(s) to generate.
#
# Usage:
#   ./run_full_vcf_test.sh                          # Default baseline mode
#   ./run_full_vcf_test.sh --mode baseline          # Baseline mode explicitly
#   ./run_full_vcf_test.sh --mode all-with-w1000    # ±1000bp window mode
#   ./run_full_vcf_test.sh -o /path/to/output       # Custom output directory
#   ./run_full_vcf_test.sh --no-plots               # Skip all heatmap generation
#   ./run_full_vcf_test.sh --plot-type distance     # Only distance heatmap
#   ./run_full_vcf_test.sh --plot-type cluster      # Only cluster heatmap
#   ./run_full_vcf_test.sh --plot-type all          # Both heatmaps (default)

set -e

# ============================================================================
# 配置變數
# ============================================================================

# 預設值
VCF_PATH="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
THREADS=120
TUMOR_BAM="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL_BAM="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF_FASTA="/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
MODE="baseline"
OUTPUT_DIR=""
# METRICS="NHD"
METRICS="L1"
GENERATE_PLOTS=true
PLOT_TYPE="all"
PLOT_THREADS=120

# 日期字串
date_str=$(date +%Y%m%d)

# 腳本目錄（用於尋找 Python 工具）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
CLUSTER_HEATMAP_SCRIPT="${PROJECT_ROOT}/tools/plot_cluster_heatmap.py"
DISTANCE_HEATMAP_SCRIPT="${PROJECT_ROOT}/tools/plot_distance_heatmap.py"
BASE_OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output"

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
        --plot-type)
            PLOT_TYPE="$2"
            shift 2
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
            echo "  --no-plots           Skip all heatmap generation"
            echo "  --plot-type TYPE     Type of heatmap to generate:"
            echo "                         all      - Both distance and cluster heatmaps (default)"
            echo "                         distance - Distance-based Read×Read heatmap only"
            echo "                         cluster  - Methylation cluster heatmap only"
            echo "  --plot-threads N     Number of threads for Python plotting (default: 64)"
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
            echo "Heatmap Types:"
            echo "  distance_heatmap: Read×Read distance matrix with dendrograms on both axes"
            echo "  cluster_heatmap:  Read×CpG methylation matrix with Y-axis dendrogram"
            echo ""
            echo "By default, both heatmaps are generated after C++ processing."
            echo "Use --no-plots to skip or --plot-type to select specific types."
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ============================================================================
# 函數定義
# ============================================================================

# 生成輸出目錄名稱（根據模式和參數）
generate_output_dir_name() {
    local mode=$1
    local threads=$2
    local base_name=""
    
    case $mode in
        baseline)
            base_name="${date_str}_vcf_baseline_t${threads}"
            ;;
        all-with-w1000)
            base_name="${date_str}_vcf_all_w1000_t${threads}"
            ;;
        all-with-w2000)
            base_name="${date_str}_vcf_all_w2000_t${threads}"
            ;;
        *)
            echo "Unknown mode: $mode" >&2
            echo "Valid modes: baseline, all-with-w1000, all-with-w2000" >&2
            exit 1
            ;;
    esac
    
    echo "${BASE_OUTPUT_DIR}/${base_name}"
}

# 為輸出目錄添加編號後綴（避免覆蓋）
add_suffix_if_exists() {
    local base_dir=$1
    local final_dir="${base_dir}"
    local counter=1
    
    # 如果目錄已存在，添加編號後綴
    while [[ -d "${final_dir}" ]]; do
        final_dir="${base_dir}_${counter}"
        counter=$((counter + 1))
    done
    
    echo "${final_dir}"
}

# ============================================================================
# 設定輸出目錄
# ============================================================================

# 如果未指定輸出目錄，根據模式和參數自動生成
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR=$(generate_output_dir_name "$MODE" "$THREADS")
    # 檢查並添加編號後綴（避免覆蓋當天已存在的輸出）
    OUTPUT_DIR=$(add_suffix_if_exists "$OUTPUT_DIR")
else
    # 如果指定了輸出目錄，直接使用（不添加編號）
    OUTPUT_DIR="${OUTPUT_DIR}"
fi


# ============================================================================
# 主程序執行
# ============================================================================

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
echo "Plot Type: ${PLOT_TYPE}"
echo "Log File: ${LOG_FILE}"
echo "---------------------------------"

# 創建輸出目錄
mkdir -p "${OUTPUT_DIR}"

# 檢查可執行文件
EXECUTABLE="/big8_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod"
if [ ! -f "${EXECUTABLE}" ]; then
    echo "Error: Executable not found at ${EXECUTABLE}" >&2
    echo "Please build the project first with:" >&2
    echo "  cd /big8_disk/liaoyoyo2001/InterSubMod/build && cmake .. && make -j" >&2
    exit 1
fi

# 構建 Metric Flags
METRIC_FLAGS=""
for m in ${METRICS}; do
    METRIC_FLAGS="${METRIC_FLAGS} --distance-metric ${m}"
done

# Step 1: Run Execution
echo ""
echo "[1] Running InterSubMod in ${MODE} mode..."

# ============================================================================
# 構建執行命令
# ============================================================================

build_command() {
    local mode=$1
    local executable=$2
    local output_dir=$3
    local threads=$4
    local metric_flags=$5
    
    case $mode in
        baseline)
            # Baseline mode: Normal filtering with debug logging
            # - Default window size (±1000bp)
            # - Default filtering parameters
            # - Debug mode enabled for logging filtered reads
            echo "${executable} \
                --tumor-bam ${TUMOR_BAM} \
                --normal-bam ${NORMAL_BAM} \
                --reference ${REF_FASTA} \
                --vcf ${VCF_PATH} \
                --output-dir ${output_dir} \
                --threads ${threads} \
                --log-level debug \
                --output-filtered-reads \
                ${metric_flags}"
            ;;
        
        all-with-w1000)
            # All-with-window mode: All filters + larger window + no-filter output
            # - Window size: ±1000bp
            # - All filters enabled (default)
            # - No-filter flag: output all reads without filtering
            # - This is useful for verification and comparison
            echo "${executable} \
                --tumor-bam ${TUMOR_BAM} \
                --normal-bam ${NORMAL_BAM} \
                --reference ${REF_FASTA} \
                --vcf ${VCF_PATH} \
                --output-dir ${output_dir} \
                --threads ${threads} \
                --window-size 1000 \
                --log-level debug \
                --output-filtered-reads \
                --no-filter \
                ${metric_flags}"
            ;;

        all-with-w2000)
            # All-with-window mode: All filters + larger window + no-filter output
            # - Window size: ±2000bp
            # - All filters enabled (default)
            # - No-filter flag: output all reads without filtering
            # - This is useful for verification and comparison
            echo "${executable} \
                --tumor-bam ${TUMOR_BAM} \
                --normal-bam ${NORMAL_BAM} \
                --reference ${REF_FASTA} \
                --vcf ${VCF_PATH} \
                --output-dir ${output_dir} \
                --threads ${threads} \
                --window-size 2000 \
                --log-level debug \
                --output-filtered-reads \
                --no-filter \
                ${metric_flags}"
            ;;
        *)
            echo "Unknown mode: $mode" >&2
            exit 1
            ;;
    esac
}

CMD=$(build_command "$MODE" "$EXECUTABLE" "$OUTPUT_DIR" "$THREADS" "$METRIC_FLAGS")

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
print_output_summary() {
    local output_dir=$1
    local mode=$2
    
    echo ""
    echo "[2] Output Summary:"
    echo "    Output directory: ${output_dir}"

    # 統計輸出文件
    if [ -d "${output_dir}" ]; then
        NUM_REGIONS=$(find "${output_dir}" -name "metadata.txt" 2>/dev/null | wc -l)
        NUM_READS_FILES=$(find "${output_dir}" -name "reads.tsv" 2>/dev/null | wc -l)
        NUM_METHYL_FILES=$(find "${output_dir}" -name "methylation.csv" 2>/dev/null | wc -l)
        NUM_FILTERED=$(find "${output_dir}" -name "filtered_reads.tsv" 2>/dev/null | wc -l)
        NUM_FORWARD=$(find "${output_dir}" -name "methylation_forward.csv" 2>/dev/null | wc -l)
        NUM_REVERSE=$(find "${output_dir}" -name "methylation_reverse.csv" 2>/dev/null | wc -l)
        
        echo "    Regions processed: ${NUM_REGIONS}"
        echo "    Reads files: ${NUM_READS_FILES}"
        echo "    Methylation matrices: ${NUM_METHYL_FILES}"
        echo "    Strand-specific matrices:"
        echo "      - Forward (+): ${NUM_FORWARD}"
        echo "      - Reverse (-): ${NUM_REVERSE}"
        if [[ "$mode" == "baseline" ]] || [[ "$mode" == "all-with-w1000" ]] || [[ "$mode" == "all-with-w2000" ]]; then
            echo "    Filtered reads logs: ${NUM_FILTERED}"
        fi
    fi
}

print_output_summary "${OUTPUT_DIR}" "${MODE}"

# Step 3: Generate heatmaps (if enabled)
generate_heatmaps() {
    local output_dir=$1
    local metrics=$2
    local plot_threads=$3
    local generate_plots=$4
    local plot_type=$5
    
    if [[ "${generate_plots}" != "true" ]]; then
        echo ""
        echo "[3] Heatmap generation skipped (use without --no-plots to enable)"
        return
    fi
    
    echo ""
    echo "[3] Generating Heatmaps..."
    
    # Check Python dependencies once
    if ! python3 -c "import matplotlib, seaborn, scipy, pandas, numpy" 2>/dev/null; then
        echo "    Warning: Python dependencies not available."
        echo "    Please install: pip install matplotlib seaborn scipy pandas numpy"
        echo "    Skipping heatmap generation."
        return
    fi

    echo "    Using ${plot_threads} threads for plotting..."
    
    PLOT_START=$(date +%s)
    
    # Loop over each metric
    for metric in ${metrics}; do
        echo ""
        echo "    >>> Processing Metric: ${metric} <<<"
        
        # Generate distance heatmap (if requested)
        if [[ "${plot_type}" == "all" ]] || [[ "${plot_type}" == "distance" ]]; then
            echo "    [3.1] Generating Distance Heatmaps (Read x Read)..."
            
            if [ -f "${DISTANCE_HEATMAP_SCRIPT}" ]; then
                python3 "${DISTANCE_HEATMAP_SCRIPT}" \
                    --output-dir "${output_dir}" \
                    --threads "${plot_threads}" \
                    --metric "${metric}" \
                    --linkage average \
                    --min-reads 10 \
                    --format png \
                    --dpi 150
                
                # Count plots in the metric-specific subdirectory if possible, 
                # but since the script handles the directory structure, we just check vaguely or skip counting precision here.
                # Just counting generic success message from script output is usually enough.
            else
                echo "         Warning: Distance heatmap script not found at ${DISTANCE_HEATMAP_SCRIPT}"
            fi
        fi
        
        # Generate cluster heatmap (if requested)
        if [[ "${plot_type}" == "all" ]] || [[ "${plot_type}" == "cluster" ]]; then
            echo "    [3.2] Generating Cluster Heatmaps (Read x CpG with dendrogram)..."
            
            if [ -f "${CLUSTER_HEATMAP_SCRIPT}" ]; then
                python3 "${CLUSTER_HEATMAP_SCRIPT}" \
                    --output-dir "${output_dir}" \
                    --threads "${plot_threads}" \
                    --metric "${metric}" \
                    --linkage average \
                    --min-reads 10 \
                    --min-cpgs 3 \
                    --format png \
                    --dpi 150
            else
                echo "         Warning: Cluster heatmap script not found at ${CLUSTER_HEATMAP_SCRIPT}"
            fi
        fi
    done
    
    PLOT_END=$(date +%s)
    PLOT_ELAPSED=$((PLOT_END - PLOT_START))
    
    echo ""
    echo "    Total plotting time: ${PLOT_ELAPSED} seconds"
}

generate_heatmaps "${OUTPUT_DIR}" "${METRICS}" "${PLOT_THREADS}" "${GENERATE_PLOTS}" "${PLOT_TYPE}"

echo ""
echo "=== Test Complete (Mode: ${MODE}) ==="
echo "Done."
