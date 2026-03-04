#!/bin/bash
# knowledge_check.sh - 偵測使用者提問是否涉及知識庫相關主題
# 用於 Claude Code UserPromptSubmit hook
# 僅在偵測到相關關鍵字時輸出提醒，避免無關提問受干擾

PROMPT=$(cat)
KB="/big8_disk/liaoyoyo2001/Knowledge"
MATCHED=()

# 檔案格式相關
if echo "$PROMPT" | grep -qiE '(VCF|\.vcf|BAM|\.bam|MM.?ML|FILTER|CIGAR|modcall|phased.?vcf|haplotag|HP.?tag|PS.?tag|QUAL|GQ|AF|VAF|allele.?freq)'; then
    MATCHED+=("檔案格式 → ${KB}/03_file_formats/")
fi

# 樣本相關
if echo "$PROMPT" | grep -qiE '(HCC1395|COLO829|H1437|H2009|HCC1937|HCC1954|HG002|tumor.?purity|subsample|truth.?set|SEQC2|NYGC|t[0-9]+_n[0-9]+)'; then
    MATCHED+=("樣本資訊 → ${KB}/02_samples/")
fi

# 工具相關
if echo "$PROMPT" | grep -qiE '(LongPhase|longphase|ClairS|clairs|DeepSomatic|deepsomatic|Clair3|clair3|MethPhaser|WhatsHap|minimap2|samtools|bcftools|dorado)'; then
    MATCHED+=("工具文件 → ${KB}/05_tools/")
fi

# 資料庫相關
if echo "$PROMPT" | grep -qiE '(PON|gnomAD|dbSNP|CoLoRSdb|reference.?genome|GRCh38|panel.?of.?normals)'; then
    MATCHED+=("資料庫 → ${KB}/04_databases/")
fi

# 流程相關
if echo "$PROMPT" | grep -qiE '(somatic.?call|variant.?call|phasing|haplotagg|methylation.?analy|benchmark|pipeline|流程|分析流程)'; then
    MATCHED+=("分析流程 → ${KB}/06_workflows/")
fi

# 資料路徑相關
if echo "$PROMPT" | grep -qiE '(/big8_disk/data|資料路徑|資料位置|目錄結構|檔案在哪|哪個目錄|storage)'; then
    MATCHED+=("資料路徑 → ${KB}/01_data_overview/")
fi

# 腳本相關
if echo "$PROMPT" | grep -qiE '(auto_run|run_benchmark|自動化腳本|batch.?script|run_vcf|run_batch)'; then
    MATCHED+=("腳本說明 → ${KB}/07_scripts/")
fi

# 僅在有匹配時輸出
if [ ${#MATCHED[@]} -gt 0 ]; then
    echo "[Knowledge Base] 偵測到相關主題，建議查閱知識庫確認細節："
    for item in "${MATCHED[@]}"; do
        echo "  - ${item}"
    done
    echo "  快速導航：先讀 ${KB}/README.md"
fi
