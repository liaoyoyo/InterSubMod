<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# HCC1395 Subsample Purity 完整驗證報告（LongPhase-S / InterSubMod / Methylation Filter）

日期：2026-02-13  
專案路徑：`/big8_disk/liaoyoyo2001/InterSubMod`

---

## 1. 本次目標與結論摘要

### 目標

1. 確認 HCC1395 每個 subsample purity 版本是否都有執行與分析。  
2. 若有缺漏，修正腳本達到可覆蓋所有 purity，且輸出不覆蓋。  
3. 驗證三階段指標（ClairS 原始、LongPhase-S 後、InterSubMod+Filter 後）與 F1 變化。  
4. 驗證 Methylation Filter 標準是否有效，以及是否需要依 F1 調整。

### 結論（先講重點）

1. **Subsample purity 版本覆蓋狀態**：
   - 已覆蓋：`t10_n40, t20_n30, t30_n20, t40_n10, t50_n00`（20~100%）
   - 阻塞：`t00_n25`（0%）缺少 `ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf`，無法進完整流程。

2. **F1 提升主因**：本批 purity 資料的 F1 提升主要來自 **LongPhase-S**（相較 ClairS 原始）。
   - 20%: `+0.0029`
   - 40%: `+0.0062`
   - 60%: `+0.0149`
   - 80%: `+0.0386`
   - 100%: `+0.0733`

3. **InterSubMod + Methylation Filter 在這批 subsample 無法評估有效性**：
   - `tp_regions_analyzed=0`、`fp_regions_analyzed=0`（全部 purity）
   - `significance_summary.csv` 只有表頭（TP/FP 檔案各 1 行）
   - 原因是 subsample BAM 本身缺 `MM/ML` 標籤，無法建立甲基化統計。

4. **腳本已修正**：`scripts/analysis/run_purity_and_standard_verification.sh` 已改為
   - 自動掃描所有 `t*_n*`
   - 以 `run_tag` 產生不覆蓋輸出
   - 顯示每個 purity 的 `skip/fail` 原因
   - 加入 `MM/ML` 前置檢查
   - 可在流程後自動刪除 `tagged.bam`（可用 `--keep-tagged-bam` 關閉）

---

## 2. 依據 Knowledge 的關鍵前提（已查證）

### 2.1 subsample 結構

根據 `Knowledge/02_samples/subsample_purity.md`：subsample 是「單一混合 BAM」，不是 tumor/normal 雙 BAM。

### 2.2 InterSubMod 甲基化分析前提

根據 `Knowledge/06_workflows/methylation_analysis.md`：InterSubMod 需要 BAM 內有 `MM/ML` 標籤。

### 2.3 本次補充知識

已更新：`Knowledge/02_samples/subsample_purity.md`，加入 2026-02-13 驗證結果：
- `t10_n40~t50_n00` 抽樣讀段下 `MM/ML=0`
- `t00_n25` 缺少 ClairS pileup VCF

另新增 issue：
- `Knowledge/issues/004-2026-02-13-subsample-mmml-and-t00-clairs-gap.md`

---

## 3. Subsample purity 覆蓋檢查

### 3.1 目錄層級

`/big8_disk/data/HCC1395/ONT/subsample/` 共有：
- `t00_n25`
- `t10_n40`
- `t20_n30`
- `t30_n20`
- `t40_n10`
- `t50_n00`

### 3.2 可執行性與實際狀態

| Subsample | 估計 Purity | ClairS pileup VCF | 已有流程輸出 | 狀態 |
|---|---:|---|---|---|
| t00_n25 | 0% | 缺少 | 無 | 阻塞（缺 Somatic VCF） |
| t10_n40 | 20% | 有 | `purity_20` | 已完成（但 InterSubMod 無有效區域） |
| t20_n30 | 40% | 有 | `purity_40` | 已完成（但 InterSubMod 無有效區域） |
| t30_n20 | 60% | 有 | `purity_60` | 已完成（但 InterSubMod 無有效區域） |
| t40_n10 | 80% | 有 | `purity_80` | 已完成（但 InterSubMod 無有效區域） |
| t50_n00 | 100% | 有 | `purity_100` | 已完成（但 InterSubMod 無有效區域） |

---

## 4. 三階段 F1 驗證（ClairS → LongPhase-S → InterSubMod+Filter）

> truth total 皆使用 `39447`（SEQC2 HCC1395 sSNV truth）

| Subsample | Purity | ClairS TP | ClairS FP | ClairS FN | ClairS F1 | LongPhase TP | LongPhase FP | LongPhase FN | LongPhase F1 | Filter TP | Filter FP | Filter FN | Filter F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| t10_n40 | 20% | 12741 | 505 | 26706 | 0.4836 | 12698 | 58 | 26749 | 0.4865 | 12698 | 58 | 26749 | 0.4865 |
| t20_n30 | 40% | 23834 | 952 | 15613 | 0.7421 | 23763 | 298 | 15684 | 0.7483 | 23763 | 298 | 15684 | 0.7483 |
| t30_n20 | 60% | 28160 | 2108 | 11287 | 0.8079 | 28051 | 689 | 11396 | 0.8228 | 28051 | 689 | 11396 | 0.8228 |
| t40_n10 | 80% | 29653 | 4517 | 9794 | 0.8056 | 29540 | 1000 | 9907 | 0.8442 | 29540 | 1000 | 9907 | 0.8442 |
| t50_n00 | 100% | 30233 | 8450 | 9214 | 0.7739 | 30061 | 1457 | 9386 | 0.8472 | 30061 | 1457 | 9386 | 0.8472 |

### 4.1 階段差異（重點）

- `ClairS -> LongPhase-S`：**全部 purity 都提升**（F1 +0.0029 ~ +0.0733）
- `LongPhase-S -> InterSubMod+Filter`：**全部 purity 為 +0.0000**（無提升）

---

## 5. Methylation Filter 有效性驗證

### 5.1 實測證據

在 20/40/60/80/100 全部 purity：

- `intersubmod_tp/significance_summary.csv` 行數 = 1（只有表頭）
- `intersubmod_fp/significance_summary.csv` 行數 = 1（只有表頭）
- `metrics.json`：`tp_regions_analyzed=0`、`fp_regions_analyzed=0`

### 5.2 根因

抽樣檢查（前 500 reads）顯示：

| Subsample | Source BAM MM | Source BAM ML | Tagged BAM MM | Tagged BAM ML |
|---|---:|---:|---:|---:|
| t10_n40 | 0 | 0 | 0 | 0 |
| t20_n30 | 0 | 0 | 0 | 0 |
| t30_n20 | 0 | 0 | 0 | 0 |
| t40_n10 | 0 | 0 | 0 | 0 |
| t50_n00 | 0 | 0 | 0 | 0 |

因此在此資料條件下：
- `(AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24)` 無可用統計輸入
- `(QUAL < 0.75)` 也無法與甲基化統計共同驗證
- 結論不是「Filter 無效」，而是「**沒有可分析的甲基化訊號**」

---

## 6. 特殊現象與解讀

1. **低 purity 到高 purity，LongPhase-S 對 FP 壓制效果越明顯**：
   - 100% purity 下 ClairS FP=8450，LongPhase-S FP=1457，F1 大幅提升。

2. **InterSubMod 看似成功執行，但統計為空**：
   - 這是資料標籤缺失（MM/ML=0）造成，不是程式崩潰。

3. **0% purity（t00_n25）是流程缺口**：
   - 缺 ClairS pileup 輸出，故無法進入 TP/FP 評估鏈。

---

## 7. 腳本修正內容（已完成）

檔案：`scripts/analysis/run_purity_and_standard_verification.sh`

已新增：

1. 自動掃描所有 `t*_n*`，不再硬編碼只跑部分 purity。  
2. 輸出路徑加入 `run_tag`（例如 `purity_t10_n40_20260213_011603`），避免覆蓋舊結果。  
3. `MM/ML` 前置檢查與狀態輸出（`purity_status.tsv`）。  
4. 缺資料（例如 `missing_somatic_vcf`）會明確標註 skip 原因。  
5. 完成後可自動清除 `tagged.bam`（預設啟用，可用 `--keep-tagged-bam` 保留）。

### 7.1 實測驗證（新腳本）

- 指令：`bash scripts/analysis/run_purity_and_standard_verification.sh --only-subdir t00_n25 --run-tag 20260213_t00_check`
- 結果：正確偵測 `missing_somatic_vcf`，並輸出狀態表  
  `output/bip8_disk_output/s-pure-pileup/HCC1395/purity_runs/20260213_t00_check/purity_status.tsv`
- 意義：確認「每個 subsample purity 版本」都可被腳本掃到，且缺資料會有明確可追蹤紀錄。

---

## 8. 任務完成度評估

### 已完成

- [x] 確認 subsample purity 版本清單與可執行性
- [x] 補齊腳本以覆蓋所有 `t*_n*`
- [x] 修正輸出防覆蓋機制
- [x] 驗證 ClairS 與 LongPhase-S 的 F1 差異
- [x] 驗證 InterSubMod 階段是否有實際分析資料
- [x] 產出繁中完整報告

### 仍受資料限制

- [ ] 在目前 subsample 資料上，無法完成「甲基化條件閾值優化」的有效比較（因 MM/ML 缺失）
- [ ] `t00_n25` 缺 ClairS pileup VCF，無法納入完整三階段

---

## 9. 建議下一步（若要真正驗證 Methylation Filter）

1. 先建立「保留 MM/ML」的 purity 混合 BAM（由 `ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam` 與 `HCC1395BL.bam` 重新 downsample/merge）。
2. 對每個 purity 重新跑 ClairS（至少要有 pileup VCF）。
3. 再執行本報告同樣流程，才可公平評估：
   - `(AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24)`
   - `(QUAL < 0.75)`
   - 與優化後閾值（例如 AD/CV/VAF 新組合）之 F1 差異。
