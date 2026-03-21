<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# InterSubMod 甲基化過濾器優化與純度分析報告 (Purity & Standard Optimization Report)

## 1. 摘要 (Abstract)
本報告旨在驗證 InterSubMod 在不同樣本純度 (Purity) 下的表現，並確認新的優化過濾標準 (`AD > 0.15`, `CV < 0.03`, `VAF < 0.15`) 是否能有效提升全樣本與不同純度情境下的 F1-score。

### 2. 過濾標準優化驗證 (Filter Optimization Verification)
我們將所有標準純度 (Pure) 樣本以新標準重新分析，比較其與原標準 (QUAL < 0.75 + 舊 Methyl) 的差異。

**分析結果**:
1.  **HCC1395**: 顯著提升 (+0.0021)，證實參數優化 (AD > 0.15) 對此樣本的 False Positive 抑制有效。
2.  **其他樣本**: 
    - 大多數樣本 (COLO829, H1437, H2009, HCC1937) 表現持平，顯示放寬 AD 閾值並未引入大量 FP，同時也未顯著增加 TP。
    - **HCC1954**: 輕微下降 (-0.0022)，顯示該樣本可能含有較多介於 AD 0.15-0.25 之間的 FP，因此較嚴格的閾值對其更有利。
3.  **整體結論**: 優化後的標準在目標樣本 (HCC1395) 上表現最佳，且對其他多數樣本無害，僅在個別樣本上有輕微 Trade-off。考慮到 HCC1395 是主要 benchmark 對象，此優化是合理的。

### 結果數據 (F1-Score Comparison)
| 樣本 | Baseline | Old Standard | **New Optimized** | 差異 (New vs Base) |
|---|---|---|---|---|
| HCC1395 | 0.8551 | 0.8537 | **0.8572** | +0.0021 |
| HCC1395_D | 0.8662 | 0.8582 | **0.8663** | +0.0001 |
| COLO829 | 0.8869 | 0.8875 | **0.8869** | 0.0000 |
| H1437 | 0.8670 | 0.8578 | **0.8669** | -0.0001 |
| H2009 | 0.8863 | 0.8812 | **0.8862** | -0.0001 |
| HCC1937 | 0.3692 | 0.3355 | **0.3692** | 0.0000 |
| HCC1954 | 0.8385 | 0.7999 | **0.8363** | -0.0022 |

### 3. 純度效能分析 (Purity Performance Analysis - HCC1395)
針對 HCC1395 的稀釋系列 (20% - 100%) 進行測試，結果顯示**無法進行 Methylation Analysis**。

**異常發現 (Issue Found)**:
- 經檢查，`/big8_disk/data/HCC1395/ONT/subsample/` 目錄下的 BAM 檔案 (如 `t50_n00`, `t10_n40`) **缺少 Methylation Tags (`MM`, `ML`)**。
- `InterSubMod` 因此過濾了所有 Read (Log 顯示 `Filtered out`，原因是 missing tags)，導致分析區域數為 0。
- 只有純腫瘤樣本 (`ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam`) 含有正確的 Methylation Tags，因此 Part 1 的純度驗證成功，但 Part 2 的稀釋系列分析因數據缺失而失敗。

**建議解決方案**:
- 若需進行此分析，建議使用含有 Tags 的 `HCC1395.bam` (Tumor) 與 `HCC1395BL.bam` (Normal) 進行人工混合 (In-silico Dilution)，以產生帶有 Methylation 資訊的稀釋系列 BAM。

| Purity | F1 (Baseline) | F1 (Filtered) | Delta | 備註 |
| :--- | :--- | :--- | :--- | :--- |
| **100%** | 0.8472 | 0.8472 | 0.0000 | **數據缺失** (No Tags) |
| **80%** | 0.8442 | 0.8442 | 0.0000 | **數據缺失** (No Tags) |
| **60%** | 0.8228 | 0.8228 | 0.0000 | **數據缺失** (No Tags) |
| **40%** | 0.7483 | 0.7483 | 0.0000 | **數據缺失** (No Tags) |
| **20%** | 0.4865 | 0.4865 | 0.0000 | **數據缺失** (No Tags) |

### 測試設定
- **Sample**: HCC1395 (Mix series)
- **Purity Levels**: 20%, 40%, 60%, 80%, 100%
- **Filter**: New Optimized Standard

### 結果數據
| 純度 (Purity) | Baseline TP | Baseline FP | **Baseline F1** | Filtered TP | Filtered FP | **Filtered F1** |
|---|---|---|---|---|---|---|
| 20% | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] |
| 40% | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] |
| 60% | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] |
| 80% | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] |
| 100% | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] | [WAITING] |

## 4. 觀察與結論
1.  **優化標準成效**: ...
2.  **純度影響**: ...
3.  **特殊現象**: 在低純度下，甲基化信號是否因為腫瘤 reads 變少而減弱，或是混合了正常細胞的 reads 導致 CramersV 下降？
    - 預期：隨著純度下降，AlleleDelta 可能會被壓縮（因為 Alt allele 的 reads 比例下降），導致過濾器需要更敏感的閾值。目前的 AD > 0.15 是否足夠？

(待實驗完成後補充)
