<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# InterSubMod 研究現況與數據盤點報告

- 文件日期：2026-02-28
- 盤點範圍：`/big8_disk/liaoyoyo2001/InterSubMod` 現有程式、輸出、測試、文檔
- 目的：整理「已完成」「待驗證」「現有結論」「可能原因」並建立下一步決策基線

## 1. 已完成研究（可重現）

### 1.1 已完成分析里程碑

| 時間 | 里程碑 | 主要產出證據 |
|---|---|---|
| 2026-01-18 | 全基因組窗口 `w5000` 顯著性分析 | `output/20260118_vcf_all_w5000_t120/significance_summary.csv`、`significance_statistics.txt` |
| 2026-01-19 | chr19 驗證 run + 多層分析 | `output/20260119_vcf_chr19_verif_t120/*`、`output/advanced_analysis_20260119/*` |
| 2026-01-19 | F1 準則比較與閾值探索 | `output/f1_evaluation_20260119/criteria_comparison.csv`、`v_sig_ratio_table.csv` |
| 2026-02-11~02-17 | 多篇內部策略報告（過濾、純度 aware、跨樣本） | `docs/reports/validated/2026/02/*`（含 `20260211_...`、`20260216_...`、`20260217_...`） |

### 1.2 本輪重跑讀取與測試驗證（2026-02-28）

1. `./build/bin/test_phase1_2`：**成功**
- BAM / FASTA 可讀
- 測試區域（chr17:7577000-7579000）抓到 151 reads
- 前 10 reads 全部可解析 MM/ML，平均 24.1 CpGs/read

2. `./build/bin/run_tests`：**105 測試中 105 通過、0 跳過**
- `tests/test_bam_reader.cpp` 已加入 fixture fallback（`data/bam/test.bam` 不存在時改用 `data/bam/HCC1395/tumor.bam`）

3. MCP（Knowledge）連線驗證：**可用**
- 透過 MCP 列表已可讀到 `knowledge://stats`、`knowledge://doc/...`
- 專案 `.mcp.json` 已包含：
  - `"command": "python3"`
  - `"args": ["/big8_disk/liaoyoyo2001/Knowledge/scripts/mcp/knowledge_server.py"]`

4. 知識庫自動提醒 hook：**可用**
- `.claude/settings.local.json` 已啟用 `UserPromptSubmit -> scripts/hooks/knowledge_check.sh`

## 2. 現有結果與結論（以數據為準）

資料來源：`output/20260118_vcf_all_w5000_t120/significance_summary.csv`（共 30,476 筆）

### 2.1 全體統計

| 指標 | 數值 |
|---|---:|
| 總區域數 | 30,476 |
| Significant=True | 1,860 (6.10%) |
| PassedGating=True | 9,636 (31.62%) |
| 平均讀段數（NumReads） | 71.4 |
| 平均 CpG 數（NumCpGs） | 97.4 |

### 2.2 VerificationClass 與顯著性的關係

| VerificationClass | 數量 | Significant=True | 顯著率 |
|---|---:|---:|---:|
| Strong | 7,271 | 1,839 | 25.29% |
| Subclone | 1,416 | 21 | 1.48% |
| Weak | 11,343 | 0 | 0.00% |
| Noise | 10,446 | 0 | 0.00% |

**目前可支持的結論**：
1. 目前顯著性邏輯幾乎只在 `Strong` 類別產生正訊號。
2. `Subclone` 的召回非常弱（僅 1.48%），與「亞克隆解析」的核心目標有落差。

### 2.3 F1 準則比較（既有評估檔）

資料來源：`output/f1_evaluation_20260119/criteria_comparison.csv`

| 準則 | TP_kept | FP_kept | F1 |
|---|---:|---:|---:|
| Current (Significant=True) | 1,860 | 101 | 0.0904 |
| HPMergedDelta<=0.1 | 27,959 | 4,194 | 0.7810 |
| HPMergedDelta<=0.05 | 24,510 | 3,723 | 0.7243 |
| HPMergedSig OR AlleleSig | 18,614 | 2,883 | 0.6109 |

**目前可支持的結論**：
1. 現行 `Significant=True` 準則 precision 高、但 recall 極低，導致 F1 極差（0.0904）。
2. 單用或優先使用 `HPMergedDelta` 在此資料集上顯著優於現行顯著性準則。

## 3. 還需要驗證的事情與數據

| 主題 | 目前狀態 | 缺口 | 優先級 |
|---|---|---|---|
| 跨樣本泛化（HCC1395 以外） | 部分文件有提及，缺統一重跑表 | 未形成同口徑 TP/FP/FN + methylation 的橫向比較 | P0 |
| 純度（purity）對甲基化判定影響 | 有 purity-aware 分析方向 | 目前常用 subsample 多為無 MM/ML，不足以做甲基化 purity 實證 | P0 |
| Subclone 偵測能力 | 已知顯著率偏低 | 缺針對 Subclone 的特徵工程與閾值校正流程 | P0 |
| FP 強訊號來源（例如 chr9 熱區） | 已有個案報告 | 缺系統性黑名單/重複區/比對偏差交叉檢驗 | P1 |
| 測試完整性 | 單元測試已全數通過 | 缺少獨立最小測試 BAM fixture（目前以 fallback 解決） | P2 |
| 文件一致性 | 已進行重整，仍需持續清理歷史連結 | 舊連結與部分 `file://` 引用仍待清理 | P1 |

## 4. 可能造成目前現象的原因（假設）

| 假設編號 | 假設內容 | 目前依據 | 需要的驗證 |
|---|---|---|---|
| H1 | 現行顯著性門檻對真陽性太嚴，造成 recall 崩落 | Current F1=0.0904，TP_kept 1,860/30,476 | 重新校正多目標閾值（Precision/Recall/F1） |
| H2 | `Strong` 分類規則混入大量非目標訊號 | FP 中仍有大量 `Strong` | 分層檢查 `Strong` 的特徵分布與區域偏倚 |
| H3 | 部分高訊號 FP 來自區域結構偏差（重複區/CNV/LOH） | chr9/特定區域反覆出現 | 引入 repeat masker、CNV、mappability 註記交叉 |
| H4 | Subclone 指標受 coverage/HP 不平衡干擾 | Subclone 顯著率僅 1.48% | 增加 HP balance + coverage-aware 特徵 |
| H5 | 標籤一致性（HPMerged vs Allele）規則未充分利用 | 現有檔案顯示組合規則可提升表現 | 將一致性特徵納入正式過濾策略 |

## 5. 目前結論可信度分級

### 高可信（已有重跑或直接數據支持）
1. 工具可讀取實際 BAM/FASTA 並解析 MM/ML（`test_phase1_2`）。
2. 單元測試主體穩定（105 pass，0 skip）。
3. 現行 `Significant=True` 在現有資料集上 F1 很低。

### 中可信（有分析跡象，待跨資料重複）
1. `HPMergedDelta` 類準則可大幅改善整體 F1。
2. LOH/CNV/coverage 異常可能是 FP 強訊號來源。

### 低可信（目前僅為合理推論）
1. 若導入新特徵後即可穩定提升跨樣本泛化。
2. 現有結論可直接外推至臨床樣本情境。

## 6. 文件與流程一致性問題（需修正）

1. 入口文件已改為 `scripts/run_vcf_all_snv.sh`，但歷史文件仍有舊路徑引用需續清理。
2. 早期文檔中的輸出結構（如 `region_*`）與現有 `output/20260118...` 風格不一致，已改為新分層策略。
3. 後續應維持單一主流程命令文件於 `docs/references/manual/`。

## 7. 本報告引用來源

### 專案內部來源
- `output/20260118_vcf_all_w5000_t120/significance_summary.csv`
- `output/20260118_vcf_all_w5000_t120/significance_statistics.txt`
- `output/f1_evaluation_20260119/criteria_comparison.csv`
- `output/advanced_analysis_20260119/multilayer_analysis_report.md`
- `output/deep_analysis_20260119/analysis_report.md`
- `docs/reports/validated/2026/02/20260211_甲基化過濾策略綜合分析報告_01.md`
- `docs/reports/validated/2026/02/20260216_跨樣本跨純度TP_FP_F1綜合分析報告_01.md`
- `docs/reports/validated/2026/02/20260217_purity_aware_analysis_report_01.md`

### Knowledge Base（根據知識庫）
- 根據 `Knowledge/05_tools/InterSubMod.md`
- 根據 `Knowledge/06_workflows/methylation_analysis.md`
- 根據 `Knowledge/02_samples/HCC1395.md`
- 根據 `Knowledge/03_file_formats/vcf_overview.md`
- 根據 `Knowledge/04_databases/seqc2_truth_set.md`
- 根據 `Knowledge/06_workflows/benchmark_workflow.md`
