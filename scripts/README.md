# scripts 目錄說明

## 主要入口

1. `pipeline/run_benchmark.sh`：建議主入口（longphase-s -> InterSubMod -> filter analysis）
2. `pipeline/run_multi_sample.sh`：多樣本批次入口
3. `analysis/run_purity_and_standard_verification.sh`：純度流程入口
4. `run_vcf_all_snv.sh`：舊版單次分析入口（legacy）
5. `run_batch_vcf_analysis.sh`：舊版 TP/FP 批次入口（legacy）

## 子目錄

1. `analysis/`：分析與報告輔助腳本
2. `pipeline/`：多階段 pipeline 腳本與步驟
3. `hooks/`：流程守門與檢查腳本
4. `mcp/`：MCP 服務腳本（pubmed/biorxiv/ensembl/github）

## Output 搬移工具

1. `analysis/migrate_output_to_bip8.sh`：將 `/big8_disk/.../InterSubMod_runs/output` 舊輸出安全搬移到 `/bip8_disk/.../InterSubMod_out/output`，並輸出 manifest（不刪檔）。

## 注意事項

1. 腳本路徑應以 repo 根目錄為相對基準。
2. 新增腳本請附最小 `--help` 或註解說明。
3. 牽涉固定資料路徑時，請於腳本頂部標註可覆寫參數。
