#!/bin/bash
# 驗證 InterSubMod 輸出的正確性

OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output"

echo "=== InterSubMod Output Verification ==="
echo ""

# 1. 檢查 region 數量
echo "[1] Checking number of regions..."
NUM_REGIONS=$(ls -d ${OUTPUT_DIR}/region_* 2>/dev/null | wc -l)
echo "    Found ${NUM_REGIONS} regions"

if [ ${NUM_REGIONS} -eq 0 ]; then
    echo "    ✗ No regions found!"
    exit 1
fi

# 2. 檢查每個 region 的檔案完整性
echo ""
echo "[2] Checking file completeness..."
MISSING_FILES=0
for region_dir in ${OUTPUT_DIR}/region_*; do
    region_name=$(basename ${region_dir})
    
    if [ ! -f "${region_dir}/metadata.txt" ]; then
        echo "    ✗ ${region_name}: missing metadata.txt"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
    
    if [ ! -f "${region_dir}/reads.tsv" ]; then
        echo "    ✗ ${region_name}: missing reads.tsv"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
    
    if [ ! -f "${region_dir}/cpg_sites.tsv" ]; then
        echo "    ✗ ${region_name}: missing cpg_sites.tsv"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
    
    if [ ! -f "${region_dir}/methylation.csv" ]; then
        echo "    ✗ ${region_name}: missing methylation.csv"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
done

if [ ${MISSING_FILES} -eq 0 ]; then
    echo "    ✓ All files present (${NUM_REGIONS} × 4 = $((NUM_REGIONS * 4)) files)"
else
    echo "    ✗ Found ${MISSING_FILES} missing files"
fi

# 3. 驗證矩陣維度一致性
echo ""
echo "[3] Verifying matrix dimensions..."
INCONSISTENT=0
for region_dir in ${OUTPUT_DIR}/region_*; do
    region_name=$(basename ${region_dir})
    
    # 從 metadata 提取維度
    if grep -q "Num Reads:" "${region_dir}/metadata.txt" && \
       grep -q "Num CpG Sites:" "${region_dir}/metadata.txt"; then
        
        NUM_READS=$(grep "Num Reads:" "${region_dir}/metadata.txt" | awk '{print $3}')
        NUM_CPGS=$(grep "Num CpG Sites:" "${region_dir}/metadata.txt" | awk '{print $4}')
        
        # 驗證 reads.tsv
        ACTUAL_READS=$(($(wc -l < "${region_dir}/reads.tsv") - 1))
        if [ ${ACTUAL_READS} -ne ${NUM_READS} ]; then
            echo "    ✗ ${region_name}: reads.tsv has ${ACTUAL_READS} rows, expected ${NUM_READS}"
            INCONSISTENT=$((INCONSISTENT + 1))
        fi
        
        # 驗證 cpg_sites.tsv
        ACTUAL_CPGS=$(($(wc -l < "${region_dir}/cpg_sites.tsv") - 1))
        if [ ${ACTUAL_CPGS} -ne ${NUM_CPGS} ]; then
            echo "    ✗ ${region_name}: cpg_sites.tsv has ${ACTUAL_CPGS} rows, expected ${NUM_CPGS}"
            INCONSISTENT=$((INCONSISTENT + 1))
        fi
        
        # 驗證 methylation.csv
        ACTUAL_MATRIX_ROWS=$(($(wc -l < "${region_dir}/methylation.csv") - 1))
        if [ ${ACTUAL_MATRIX_ROWS} -ne ${NUM_READS} ]; then
            echo "    ✗ ${region_name}: methylation.csv has ${ACTUAL_MATRIX_ROWS} rows, expected ${NUM_READS}"
            INCONSISTENT=$((INCONSISTENT + 1))
        fi
    fi
done

if [ ${INCONSISTENT} -eq 0 ]; then
    echo "    ✓ All dimensions consistent"
else
    echo "    ✗ Found ${INCONSISTENT} inconsistencies"
fi

# 4. 抽樣檢查矩陣內容
echo ""
echo "[4] Sampling matrix content..."
SAMPLE_REGION="${OUTPUT_DIR}/region_0000"
if [ -f "${SAMPLE_REGION}/methylation.csv" ]; then
    # 計算 NA 數量
    NA_COUNT=$(grep -o "NA" "${SAMPLE_REGION}/methylation.csv" | wc -l)
    TOTAL_CELLS=$(($(head -1 "${SAMPLE_REGION}/methylation.csv" | tr ',' '\n' | wc -l) - 1))
    TOTAL_ROWS=$(($(wc -l < "${SAMPLE_REGION}/methylation.csv") - 1))
    TOTAL_VALUES=$((TOTAL_CELLS * TOTAL_ROWS))
    
    echo "    Region 0000 matrix:"
    echo "      - Dimensions: ${TOTAL_ROWS} × ${TOTAL_CELLS}"
    echo "      - Total cells: ${TOTAL_VALUES}"
    echo "      - NA values: ${NA_COUNT}"
    echo "      - Data values: $((TOTAL_VALUES - NA_COUNT))"
    echo "      - Sparsity: $(echo "scale=2; ${NA_COUNT} * 100 / ${TOTAL_VALUES}" | bc)%"
fi

# 5. 總結
echo ""
echo "=== Summary ==="
echo "Regions processed: ${NUM_REGIONS}"
echo "Missing files: ${MISSING_FILES}"
echo "Dimension inconsistencies: ${INCONSISTENT}"

if [ ${MISSING_FILES} -eq 0 ] && [ ${INCONSISTENT} -eq 0 ]; then
    echo ""
    echo "✓ ALL CHECKS PASSED!"
    exit 0
else
    echo ""
    echo "✗ Some checks failed"
    exit 1
fi

