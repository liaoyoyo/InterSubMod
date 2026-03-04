<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# HCC1395 Subsample Purity（run-tag: 20260213_repurity_full）完整分析報告

> 備註：本檔是 `HCC1395/ONT`（非 Dorado）資料的結果，該來源 BAM 缺少 MM/ML，無法有效驗證甲基化過濾。
> 已改用 `HCC1395_DORADO` 完整重跑 purity（含 MM/ML）並產出新報告：
> `docs/ai_sessions/2026/02/2026-02-13_HCC1395_DORADO_subsample_purity_20260213_dorado_purity_full_完整分析報告.md`
- 產生時間：2026-02-13
- 執行腳本：`scripts/analysis/run_purity_and_standard_verification.sh`
- Run-tag：`20260213_repurity_full`（不覆蓋舊結果）
- 資料來源：`/big8_disk/data/HCC1395/ONT/subsample/`
- 真值：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- 高可信區域：`/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`
- Truth SNV（BED 內）總數：**39447**

## 0. 知識庫依據
- 根據 `Knowledge/02_samples/subsample_purity.md`：`t00_n25` 為 0% tumor 的 normal baseline，不納入腫瘤 purity 驗證主流程。
- 根據 `Knowledge/06_workflows/methylation_analysis.md`：InterSubMod 依賴 BAM 中 MM/ML 修飾標籤；若無 MM/ML，無法進行有效甲基化區域統計。

## 1. 執行覆蓋確認
- `t00_n25` 已排除（未執行）。
- 本次執行 purity：`t10_n40`, `t20_n30`, `t30_n20`, `t40_n10`, `t50_n00`（共 5 組）。
- 每組均完成 Step01(LongPhase-S) / Step02(InterSubMod) / Step03(Filter Analysis)；`cleanup=REMOVED`。
- 每組結束後 `tagged.bam` 均已自動刪除（節省空間）。

## 2. 各 purity 三階段指標（TP/FP/FN/P/R/F1）
| Purity | MM/ML | ClairS TP | ClairS FP | ClairS FN | ClairS F1 | LongPhase TP | LongPhase FP | LongPhase FN | LongPhase F1 | InterSubMod 可分析 TP/FP regions | 最終 TP | 最終 FP | 最終 FN | 最終 F1 | ΔF1 (ClairS→LongPhase) | ΔF1 (LongPhase→最終) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20% | 0/0 | 12741 | 505 | 26706 | 0.4836 | 12698 | 58 | 26749 | 0.4865 | 0/0 | 12698 | 58 | 26749 | 0.4865 | +0.0029 | +0.0000 |
| 40% | 0/0 | 23834 | 952 | 15613 | 0.7421 | 23763 | 298 | 15684 | 0.7483 | 0/0 | 23763 | 298 | 15684 | 0.7483 | +0.0062 | +0.0000 |
| 60% | 0/0 | 28160 | 2108 | 11287 | 0.8079 | 28051 | 689 | 11396 | 0.8228 | 0/0 | 28051 | 689 | 11396 | 0.8228 | +0.0149 | +0.0000 |
| 80% | 0/0 | 29653 | 4517 | 9794 | 0.8056 | 29540 | 1000 | 9907 | 0.8442 | 0/0 | 29540 | 1000 | 9907 | 0.8442 | +0.0386 | +0.0000 |
| 100% | 0/0 | 30233 | 8450 | 9214 | 0.7739 | 30061 | 1457 | 9386 | 0.8472 | 0/0 | 30061 | 1457 | 9386 | 0.8472 | +0.0733 | +0.0000 |

## 3. 關鍵結論
1. **所有 subsample purity（排除 t00_n25）皆已完成執行，且不覆蓋舊結果。**
2. **LongPhase-S 在所有 purity 對 ClairS(PASS) 皆有 F1 提升**（+0.0029 到 +0.0733）。
3. **InterSubMod + Methylation Filter 在本批資料無額外 F1 提升**（全部 0.0000）。
4. 原因：5 組 subsample 的來源 BAM 皆無 MM/ML（`MM=0, ML=0`），導致 `tp_regions_analyzed=0, fp_regions_analyzed=0`，Methylation 訊號無法進入判定。
5. 這代表此批 `ONT/subsample` 目前只能驗證「LongPhase-S 改善」，**不能驗證甲基化過濾是否能再提升 F1**。

## 4. 特殊現象與特徵
- 隨 purity 提升（20%→100%），LongPhase 的 TP 與 F1 單調上升：
  - 20%：LongPhase TP=12698, FP=58, F1=0.4865
  - 40%：LongPhase TP=23763, FP=298, F1=0.7483
  - 60%：LongPhase TP=28051, FP=689, F1=0.8228
  - 80%：LongPhase TP=29540, FP=1000, F1=0.8442
  - 100%：LongPhase TP=30061, FP=1457, F1=0.8472
- 100% purity（t50_n00）的 LongPhase F1 為 0.8472，已非常接近純樣本基準。
- `tp_regions_analyzed=0` 與 `fp_regions_analyzed=0` 在 5 組全部一致，顯示不是單一 purity 的偶發問題。

## 5. 對「Methylation Filter 是否有效」的判定
- 在 **本次 subsample 資料條件** 下：**無法判定有效**（不是無效，而是無甲基化可用訊號）。
- 在已有 MM/ML 的資料（先前 pure-bam 實驗）才可有效測試過濾能力。
- 若要完成 purity 條件下的甲基化驗證，需先確保 subsample 生成流程保留 MM/ML 標籤。

## 6. 主要輸出檔
- 執行狀態表：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395/purity_runs/20260213_repurity_full/purity_status.tsv`
- 階段指標表：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395/purity_runs/20260213_repurity_full/purity_stage_metrics.tsv`
- 各 purity 輸出目錄：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395/purity_t*_20260213_repurity_full/`

## 7. 甲基化標籤補充驗證（純樣本 vs subsample）
- 純樣本 tumor BAM（HCC1395/HCC1395_DORADO/COLO829/H1437/H2009/HCC1937/HCC1954）抽樣 500 reads 皆為 `MM=500, ML=500`。
- HCC1395 subsample（`t10_n40`, `t50_n00`）抽樣 500 reads 為 `MM=0, ML=0`。
- 結論：本次無法由 subsample purity 驗證 Methylation Filter，不是因為 InterSubMod 參數失效，而是上游 subsample BAM 缺少甲基化標籤。
