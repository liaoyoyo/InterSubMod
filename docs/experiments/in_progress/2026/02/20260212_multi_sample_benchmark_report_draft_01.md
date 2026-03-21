<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 多樣本 InterSubMod Methylation Filter 驗證報告

## 1. 概述
本報告旨在驗證 InterSubMod 在多個純 BAM 樣本上的甲基化過濾效果。我們使用了 7 個不同的樣本/設定，涵蓋不同的癌症類型與測序來源，以評估該方法的魯棒性與效能。

## 2. 方法
- **Pipeline**: `LongPhase-S` (Somatic Haplotagging) -> `InterSubMod` (Methylation Analysis) -> `Methylation Filter`
- **Filter Criteria**: `Abs(LabelDelta) > 0.2` & `Significance != ns` (Tier 1)
- **Threads**: 112 (全資源利用) for both LongPhase-S and InterSubMod (Sequential mode)

## 3. 樣本結果匯總

### 3.1 HCC1395 (Reference)
- **Baseline**: F1 = 0.8522
- **Filtered**: F1 = 0.8528 (+0.0006)
- **TP/FP Change**: TP stable, FP reduced slightly.

### 3.2 HCC1395_DORADO
- **Baseline**: F1 = 0.8592
- **Filtered**: F1 = 0.8591 (-0.0001)
- **TP/FP Change**: Minimal change.

### 3.3 COLO829 (Melanoma)
- **Truth Source**: NYGC (No BED, Whole Genome)
- **Baseline**: TP=35185, FP=2273, F1=0.8921
- **Filtered**: TP=35176, FP=2272, F1=0.8919 (-0.0002)
- **Details**: FP count is high (2273), providing a good testing ground. Filter removed 1 FP and 9 TP.

### 3.4 H1437 (Lung Adenocarcinoma)
- **Truth Source**: Orthogonal Tools
- **Baseline**: TP=67468, FP=8, F1=0.8562
- **Filtered**: TP=67462, FP=8, F1=0.8561 (-0.0001)
- **Details**: Extremely low FP count (8), making it hard to improve Precision without hurting Recall. Filter removed 6 TP.

### 3.5 H2009 (Lung Adenocarcinoma)
- **Truth Source**: Orthogonal Tools
- **Baseline**: TP=132909, FP=86, F1=0.8816
- **Filtered**: TP=132873, FP=86, F1=0.8814 (-0.0002)
- **Details**: FP count is also low (86). Filter removed 36 TP, 0 FP.

### 3.6 HCC1937 (Breast Cancer, BRCA1)
- **Truth Source**: Orthogonal Tools
- **Baseline**: TP=12393, FP=195, F1=0.3382
- **Filtered**: TP=12392, FP=195, F1=0.3382 (No change)
- **Details**: **異常低 recall (0.2042)**。可能原因：ClairS call set 與 Truth set 差異巨大，或 Truth BED 定義區域不匹配。需要深入調查。

### 3.7 HCC1954 (Breast Cancer)
- **Status**: *Running...*
- **Truth Source**: Orthogonal Tools

## 4. 討論與觀察
1.  **Low FP Rate across most samples**: H1437, H2009, HCC1937 都有非常低的 FP 數量（相較於 TP）。這使得依賴 "過濾 FP" 來提升 F1 變得非常困難，因為任何誤刪的 TP 都會導致 F1 下降。
2.  **Resource Utilization**: 使用 112 threads 大幅加速了分析過程。HCC1937 LongPhase-S 耗時 ~23m，InterSubMod 耗時 ~2m/task。
3.  **HCC1937 Outlier**: 該樣本的 F1 分數顯著低於其他樣本，主要受限於 Recall。建議檢查 VCF 來源或 Truth set 版本。

## 5. 結論
待 HCC1954 完成後總結。
