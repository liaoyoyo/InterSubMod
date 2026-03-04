<!--
建立時間: 2026-03-04 01:32
更新時間: 2026-03-04 01:32
狀態: finalized
目標: 盤點 scripts/tools、整理可重構流程，並完成 output 空間遷移與差異報告
範圍:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/
  - /big8_disk/liaoyoyo2001/InterSubMod/tools/
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/
  - /bip8_disk/liaoyoyo2001/InterSubMod_out/output/
關聯文件:
  - docs/reports/validated/2026/03/assets/20260304_scripts_inventory_01.tsv
  - docs/reports/validated/2026/03/assets/20260304_tools_inventory_01.tsv
  - docs/reports/validated/2026/03/assets/20260303_output_migration_big8_to_bip8_manifest_01.tsv
-->

# scripts/tools/output 流程重整與空間遷移報告

## 1. 依據與查核來源

本次判讀腳本流程與路徑策略時，先查核 Knowledge Base：

1. `Knowledge/07_script_docs/data_processing_scripts.md`
2. `Knowledge/07_script_docs/run_benchmark_sh.md`
3. `Knowledge/05_tools/InterSubMod.md`
4. `Knowledge/01_data_overview/directory_structure.md`

結論：現況流程可分為「舊版單次腳本」與「新 benchmark pipeline」兩條線，且 output 實際上已朝 `bip8` 路徑遷移。

## 2. scripts 盤點結果

盤點檔：`20260304_scripts_inventory_01.tsv`

1. scripts 檔案總數：31（TSV 含表頭共 32 行）
2. 類型分佈：
   - `sh`: 16
   - `py`: 12
   - `pyc`: 2
   - `md`: 1
3. pipeline 關鍵步驟：
   - `00_prepare_germline.sh`
   - `01_longphase_s.sh`
   - `02_intersubmod.sh`
   - `03_filter_analysis.py`
   - `04_cleanup.sh`

### 2.1 執行差異重點（入口腳本）

1. `scripts/run_vcf_all_snv.sh`
   - 直接呼叫 C++ binary，偏單次/手動模式，硬編碼路徑多。
2. `scripts/run_batch_vcf_analysis.sh`
   - 同時跑 TP/FP，後接 `tools/compare_vcf_results.py`，屬舊版批次入口。
3. `scripts/pipeline/run_benchmark.sh`
   - 新版主流程，支援 `--sample/--mode/--skip-*`，含 MM/ML guard 與 steps 拆分，較適合標準化。
4. `scripts/analysis/run_purity_and_standard_verification.sh`
   - 針對 purity/subsample 的批次驗證包裝器，屬進階流程。

## 3. tools 盤點結果

盤點檔：`20260304_tools_inventory_01.tsv`

1. tools 檔案總數：8（TSV 含表頭共 9 行）
2. 類型分佈：
   - `py`: 7
   - `pyc`: 1
3. 功能分群：
   - 視覺化：`plot_cluster_heatmap.py`、`plot_distance_heatmap.py`
   - 比較分析：`compare_vcf_results.py`、`analyze_depth_effect.py`
   - 候選篩選：`find_verification_candidates.py`

## 4. 重構可行性與問題熱點

熱點檔：`20260304_scripts_tools_path_and_cleanup_hotspots_01.txt`

### 4.1 可行性結論

可重構成單一標準流程，建議以 `pipeline/run_benchmark.sh` 為主入口，舊入口標記 legacy。

### 4.2 主要問題

1. 硬編碼路徑引用仍多（124 處命中）。
2. `04_cleanup.sh` 與部分 analysis 腳本仍使用 `rm` 清理。
3. output 路徑有歷史層疊（`output -> .../InterSubMod_runs/output -> bip8_disk_output -> /bip8_disk/...`），對新使用者不夠直觀。

## 5. output 空間遷移（已執行）

### 5.1 遷移前

1. `/big8_disk/liaoyoyo2001/InterSubMod_runs/output`：8.2G
2. 最大目錄：`20260118_vcf_all_w5000_t120`（8.1G）

基線檔：

1. `20260303_big8_output_topdir_size_before_move_01.tsv`
2. `20260303_bip8_output_topdir_name_before_move_01.txt`

### 5.2 遷移動作

1. 目標路徑：`/bip8_disk/liaoyoyo2001/InterSubMod_out/output/`
2. 搬移目錄數：12
3. 搬移總量：8.14GB
4. `bip8_disk_output` symlink 已保留（未搬移）
5. 另外建立 12 個相容性 symlink 回 `InterSubMod_runs/output`，避免舊絕對路徑失效

遷移紀錄：

1. `20260303_output_migration_big8_to_bip8_manifest_01.tsv`
2. `20260303_output_migration_postcheck_01.tsv`
3. `20260303_output_migration_compat_symlink_manifest_01.tsv`

### 5.3 遷移後差異

1. `big8` 來源目錄現況：僅保留 symlink（含相容性連結），目錄大小 4.0K
2. `bip8` top-level 目錄數：24 -> 36（新增 12）

差異檔：

1. `20260303_big8_output_topdir_after_move_01.tsv`
2. `20260303_bip8_output_topdir_name_after_move_01.txt`
3. `20260303_bip8_output_name_diff_after_move_01.md`

## 6. 本次實作調整

1. 新增 `scripts/analysis/migrate_output_to_bip8.sh`
   - 支援 `--dry-run`
   - 預設跳過 symlink
   - 自動輸出 manifest
2. 更新 `scripts/README.md`
   - 明確標示建議入口（pipeline）與 legacy 入口
   - 新增 output 搬移工具說明
