<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# InterSubMod 知識庫與最新文獻對照分析報告

- 文件日期：2026-02-28
- 目的：把本專案現況與 2025-2026 最新研究結論對齊，找出「可直接採納」與「必須再驗證」項目
- 分析原則：先以本地 Knowledge Base 確認流程/資料事實，再對照外部一手文獻

## 1. 知識庫基線（我們現在的已知事實）

1. 根據 `Knowledge/05_tools/InterSubMod.md`
- InterSubMod 目標是以 long-read read-level methylation 做 subclonal resolution。
- 核心依賴包含 MM/ML、HP tag 與距離矩陣分析。

2. 根據 `Knowledge/06_workflows/methylation_analysis.md`
- 甲基化分析流程高度依賴「可用的 haplotagged + MM/ML BAM」。
- 分析輸出應形成可重現的 feature matrix 與報告層。

3. 根據 `Knowledge/02_samples/HCC1395.md`
- HCC1395 有多套 ONT 資料；其中並非所有子集都帶可用甲基化標籤。
- purity subsample 常見資料來源存在「無 MM/ML」限制，這會直接影響甲基化研究可驗證性。

4. 根據 `Knowledge/04_databases/seqc2_truth_set.md` 與 `Knowledge/06_workflows/benchmark_workflow.md`
- 評估需嚴格限制於 High-Confidence Regions。
- benchmark 口徑（PASS、SNV/INDEL 分開、正規化）必須一致，否則結果不可橫向比較。

## 2. 最新外部研究（2025-2026）重點

| 年份 | 文獻 | 主要結論（摘要層） | 對 InterSubMod 的直接意義 |
|---|---|---|---|
| 2025 | ClairS-TO（Nature Communications） | long-read tumor-only somatic small variant calling 在 HCC1395/COLO829 benchmark 上具穩定表現 | 我們的 TP/FP 分流可借鏡其 filter 與 verdict 設計 |
| 2025 | DeepSomatic（Nature Biotechnology） | 多平台（含 ONT）somatic small variant calling 可達高準確度，並支援 tumor-only/tumor-normal 場景 | 支持「多特徵融合優於單一規則」的方向 |
| 2025 | LongPhase-S（bioRxiv） | 以 somatic haplotyping + purity estimation 進行 variant recalibration，在多資料集有 F1 改善 | 支持把 HP/purity 資訊納入 InterSubMod 決策特徵 |
| 2025 | Benchmarking long-read variant calling in diploid and polyploid genomes（BMC Genomics） | 高複雜/重複區域仍是長讀段 variant calling 挑戰，工具偏好存在差異 | 對應我們觀察到的「特定區域 FP 強訊號」問題 |
| 2025 | Toward clinical long-read genome sequencing for rare diseases（Nature Genetics） | 臨床化落地需要可重現流程、品質控管與跨場景穩定度 | 呼應我們「跨樣本同口徑評估」與文件規範化的必要性 |

## 3. 與本專案現有結果的對照

### 3.1 一致處（外部文獻支持我們目前觀察）

1. **單一顯著性規則不足以支撐高 F1**
- 我們：`Current (Significant=True)` F1=0.0904。
- 文獻：DeepSomatic 與多篇 benchmark 都是多訊號整合才達高表現。

2. **區域複雜度會驅動錯誤型態**
- 我們：特定染色體與熱區出現異常行為（如高訊號 FP）。
- 文獻：2025 benchmark 明確指出重複/複雜區域是主要難點。

3. **甲基化是有價值訊號，但不是獨立萬能訊號**
- 我們：`HPMergedDelta` 等規則可顯著改善，但 Subclone 召回仍弱。
- 文獻：主流高效模型都採多特徵融合，不依賴單一 methylation 門檻。

### 3.2 落差處（我們尚未達到文獻建議）

1. 尚未建立統一的跨樣本、跨純度、跨區域評估 protocol（同口徑 metrics + 資料品質標註）。
2. 尚未把區域註解（repeat/CNV/mappability）系統化納入判定。
3. 尚未建立「可重現」的模型比較框架（rule-based vs lightweight ML vs ensemble）。

## 4. 目前可採納的研究方向（由對照推導）

1. **多特徵融合優先於單閾值**
- 把 `Significant` 從「硬過濾」轉成一個 feature。
- 與 `HPMergedDelta`、HP balance、coverage、區域註解一起建模。

2. **區域風險建模（region-aware）**
- 先建立黑名單/灰名單（重複、高背景噪音、已知 problematic hotspots）。
- 在評估報告中分開統計「可比區」與「高風險區」。

3. **建立 purity-aware methylation 測試資料**
- 目前 purity subsample 的 MM/ML 限制，是關鍵瓶頸。
- 應優先補齊「可帶甲基化標籤的 purity 系列資料」。

4. **文獻對齊的實驗設計**
- 目標不是追求單點 F1，而是建立「跨樣本穩定改善 + 可解釋原因」的證據鏈。

## 5. 需要再驗證的假設（外部研究導向）

| 假設 | 驗證方式 | 驗收標準 |
|---|---|---|
| H-A | 多特徵融合可優於 `HPMergedDelta<=0.1` 單規則 | 在至少 3 個樣本中 F1/PR-AUC 同向提升，且非單一樣本特例 |
| H-B | 區域風險註解可降低高訊號 FP | 高風險區 FP 下降明顯，且 TP 損失可控 |
| H-C | purity-aware 特徵可改善低純度場景穩定度 | 20-80% purity 區間內，性能退化曲線斜率降低 |
| H-D | Subclone 特徵補強可提升亞克隆辨識 | `VerificationClass=Subclone` 召回顯著提升且 precision 不崩落 |

## 6. 外部來源連結

1. ClairS-TO（Nature Communications）  
https://www.nature.com/articles/s41467-025-64547-z

2. DeepSomatic（Nature Biotechnology）  
https://www.nature.com/articles/s41587-025-02839-x

3. LongPhase-S（bioRxiv）  
https://www.biorxiv.org/content/10.1101/2025.11.20.689492v1

4. Benchmarking long-read variant calling in diploid and polyploid genomes（BMC Genomics）  
https://link.springer.com/article/10.1186/s12864-025-12259-5

5. Toward clinical long-read genome sequencing for rare diseases（Nature Genetics）  
https://www.nature.com/articles/s41588-025-02160-y

6. Google Research: DeepSomatic（公開摘要）  
https://research.google/pubs/accurate-long-read-somatic-mutation-calling-by-integrating-methylation-signals/

## 7. 限制與說明

1. 外部文獻連結已於 2026-03-01 以 DOI / 期刊官方頁面進行可達性重驗證；若後續期刊調整網址，需再更新。
2. 本文件中的「可採納方向」屬研究規劃推論，非已驗證事實。
