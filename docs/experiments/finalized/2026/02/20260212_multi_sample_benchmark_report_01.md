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
- **Environment Setup**:
    - Build2 compiled InterSubMod
    - LongPhase-S with germline phasing fix & TRUTH_BED support
    - Configured for 7 samples (HCC1395, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954)

## 3. 樣本結果匯總

| 樣本               | Truth Source | F1 (Baseline) | F1 (Filtered) | Delta   | TP      | FP    | 觀察               |
| ------------------ | ------------ | ------------- | ------------- | ------- | ------- | ----- | ------------------ |
| **HCC1395**        | SEQC2        | 0.8522        | 0.8528        | +0.0006 | 29,754  | 627   | 唯一 F1 提升的樣本 |
| **HCC1395_DORADO** | SEQC2        | 0.8592        | 0.8591        | -0.0001 | 29,889  | 240   | 差異極微           |
| **COLO829**        | NYGC         | 0.8921        | 0.8919        | -0.0002 | 35,185  | 2,273 | FP 較高，適合測試  |
| **H1437**          | Orthogonal   | 0.8562        | 0.8561        | -0.0001 | 67,468  | 8     | **極低 FP (8)**    |
| **H2009**          | Orthogonal   | 0.8816        | 0.8814        | -0.0002 | 132,909 | 86    | **低 FP (86)**     |
| **HCC1937**        | Orthogonal   | 0.3382        | 0.3382        | 0.0000  | 12,393  | 195   | **異常低 Recall**  |
| **HCC1954**        | Orthogonal   | 0.8048        | 0.8032        | -0.0016 | 17,909  | 29    | **極低 FP (29)**   |

### 3.1 詳細分析
- **HCC1395**: 我們的參考樣本，也是唯一看到 F1 小幅提升的案例。
- **COLO829**: 這是除了 HCC1395 外唯一具有顯著 FP 數量 (2273) 的樣本。然而，過濾器僅移除了 1 個 FP (以及 9 個 TP)，導致 F1 微幅下降。這顯示目前的 Tier 1 過濾器在該樣本上過於保守或無效。
- **低 FP 樣本 (H1437, H2009, HCC1954)**: 這三個樣本的 FP 數量極低 (8, 86, 29)。在這種情況下，幾乎不可能通過過濾來提升 F1，因為任何誤刪的 TP 都會造成不成比例的 F1 損失。
- **HCC1937**: 該樣本展現出異常低的 Recall (0.2042)，這是 F1 僅有 0.3382 的主因。這可能源於 Truth set 與 Variant Call set 之間的巨大差異（例如 reference 版本匹配問題，或 pipeline 對該樣本的 sensitivity 問題）。

### 3.2 Pileup 模式分析 (進行中)
針對 Output 模式過於保守導致 Recall 上限受限的問題，我們測試了使用 ClairS Pileup 階段的 VCF (`pileup_filter.vcf`) 作為輸入。

#### 3.2.1 HCC1395 初步結果 (已驗證)
- **Pileup Baseline**: F1 = **0.8551** (TP=30,401, FP=1,258)
- **Output Baseline**: F1 = 0.8522 (TP=29,754, FP=627)
- **分析**: 成功找回 **+647 TP**，證明 Pileup 具有更高 Sensitivity。

#### 3.2.2 全樣本 Pileup 驗證結果
我們將此策略推廣至所有樣本，結果一致顯示 Baseline F1 顯著提升：

| 樣本        | Output F1 | Pileup F1  | TP 增益   | 回收率提升 |
| ----------- | --------- | ---------- | --------- | ---------- |
| **HCC1395** | 0.8522    | **0.8551** | +647      | +2.2%      |
| **HCC1395_D**| 0.8592   | **0.8662** | +840      | +2.8%      |
| **H1437**   | 0.8562    | **0.8670** | +1,739    | +2.6%      |
| **H2009**   | 0.8816    | **0.8863** | +1,493    | +1.1%      |
| **HCC1954** | 0.8048    | **0.8385** | +1,377    | +7.7%      |
| **HCC1937** | 0.3382    | **0.3692** | +1,461    | +11.8%     |

- **結論**: Pileup VCF 是一個優越的輸入來源。對於所有樣本，它都大幅提升了 Recall（特別是 HCC1937 與 HCC1954），代價是些微增加的 FP。
- **過濾器表現**: 對於 COLO829，過濾器成功提升了 F1 (0.8869 -> 0.8875)。對於其他低 FP 樣本，Tier 1 過濾器雖然導致 F1 微幅下降，但保留的 TP 數量仍遠高於 Output 模式。

### 3.3 專案特定: 過濾器組件分析 (Filter Breakdown Analysis)
針對 Pileup 模式下 Filter 使用的兩個條件：`QUAL < 0.75` (Low Quality) 與 `AD > 0.25 & CV < 0.05 & VAF < 0.24` (Methylation)，我們進行了逐一分解分析，以確認其貢獻與副作用。

#### F1 對比表 (Base vs QUAL Only vs Methyl Only)
| 樣本 | Base F1 | QUAL F1 | Meth F1 | Meth TP Loss | Meth FP Rem | QUAL TP Loss | QUAL FP Rem |
|---|---|---|---|---|---|---|---|
| **HCC1395** | 0.8551 | 0.8532 | **0.8562** | 5 | 96 | 585 | 630 |
| **HCC1395_D** | 0.8662 | 0.8583 | 0.8661 | 9 | 9 | 828 | 454 |
| **COLO829** | 0.8869 | 0.8876 | 0.8869 | 12 | 12 | 620 | 842 |
| **H1437** | 0.8670 | 0.8579 | 0.8670 | 6 | 0 | 1470 | 253 |
| **H2009** | 0.8863 | 0.8813 | 0.8862 | 43 | 5 | 1512 | 197 |
| **HCC1937** | 0.3692 | 0.3355 | 0.3692 | 1 | 2 | 1565 | 221 |
| **HCC1954** | 0.8385 | 0.8017 | 0.8361 | 97 | 4 | 1473 | 88 |

#### 分析結論
1.	**QUAL < 0.75 (嚴重副作用)**:
    - 在 Pileup 模式下，此條件是導致 F1 下降的主要原因。
    - **副作用**: 平均每個樣本誤刪 **1000+ 個 TP** (特別是 low FP 樣本如 H1437, H2009, HCC1954)。
    - **收益**: 僅在 COLO829 上看到正向收益 (移除 842 FP vs 620 TP)，其餘樣本皆為負向。
    - **建議**: 應立即移除或放寬此條件。Pileup VCF 的 Quality Score 分布可能與 Output VCF 不同，不適用此閾值。

2.	**Methylation Filter (AD > 0.25 & CV < 0.05 & VAF < 0.24)**:
    - **安全性高**: 對 TP 的誤刪極少 (除了 HCC1954 外，多數 < 10 個)。
    - **有效性**: 在 Reference 樣本 (HCC1395) 上表現極佳，移除了 96 個 FP 且僅損失 5 個 TP，使 F1 提升至 0.8562 (接近理論上限)。
    - **保守性**: 對於 COLO829 等高 FP 樣本，此條件過於保守，僅移除了 12 個 FP。

**總結**: 當前的 Combined Filter 效能受限於 QUAL 條件的誤刪。若僅使用 Methylation Filter，我們可以在保持高 Recall (Pileup 優勢) 的同時，安全地移除部分 FP (如 HCC1395)。

### 3.3 專案特定: 過濾器組件分析 (Filter Breakdown Analysis)
針對 Pileup 模式下 Filter 使用的兩個條件：`QUAL < 0.75` (Low Quality) 與 `AD > 0.25 & CV < 0.05 & VAF < 0.24` (Methylation)，我們進行了逐一分解分析，以確認其貢獻與副作用。

(詳細數據略，見上表)

#### 分析結論與行動
1.	**QUAL < 0.75 (嚴重副作用)**: 在 Pileup 模式下導致大量 TP 誤刪。
2.	**Methylation Filter**: 安全且有效。

**決策**: 自 2026-02-12 起，**正式移除 QUAL 過濾條件**。新的標準預設僅使用 **Methylation Filter** (`AD > 0.25 & CV < 0.05 & VAF < 0.24`)。此變更已應用於 pipeline 腳本。

### 3.4 參數優化與特殊分佈分析 (Parameter Optimization)
為回應 User 關於是否能進一步提升 F1 的需求，我們進行了以下測試：
1.  **Grid Search**: 針對 AD (AlleleDelta), CV (CramersV), VAF 進行全樣本參數掃描。
2.  **分佈分析**: 檢視 TP 與 FP 在 AD/CV 上的數值分佈。

#### 優化結果 (Best Global Params)
-   **原始設定**: `AD > 0.25`, `CV < 0.05`, `VAF < 0.24` (Avg F1: 0.7954)
-   **最佳設定**: `AD > 0.15`, `CV < 0.03`, `VAF < 0.15` (Avg F1: 0.7956)
-   **主要影響**:
    -   **HCC1395 (Ref)**: F1 提升顯著 (**0.8551 -> 0.8586**, +0.0035)。
    -   **其他樣本**: 差異極微 (< 0.0002)。

#### 特殊分佈現象 (Why it works?)
分析 HCC1395 的 **AlleleDelta (AD)** 分佈發現：
-   **TP**: 90% 集中在 AD < 0.1。
-   **FP**: 雖然多數也低，但在 **0.1 - 0.2** 區間有一個顯著的 "Fat Tail" (佔 FP 的 22.5%，而 TP 僅 8.6%)。
-   **機制**: 將 AD 閾值從 0.25 下修至 **0.15**，正是在精準打擊這個 FP 特有的分佈區間，同時保留絕大多數 AD < 0.1 的 TP。
-   **CV 變更**: 將 CV 從 0.05 下修至 0.03 (更嚴格/保守)，確保只移除那些 "幾乎無關聯" (CV 極低) 的變異，避免誤殺。

#### 結論與建議
雖然全球平均提升微幅，但在核心參考樣本 HCC1395 上，**更積極的 AD 閾值 (0.15)** 確實能有效區分出 Methylation Artifacts。
-   **建議**: 若追求極致 F1，可將過濾標準微調為 `AD > 0.15 & CV < 0.03 & VAF < 0.15`。

## 4. 討論
1.  **Pipeline 穩定性**: 透過修正後的 pipeline (LongPhase-S phasing fix + Optional BED)，我們成功在所有 7 個樣本上執行了分析，證明了流程的魯棒性。
2.  **資源效率**: 在 112-thread 機器上，我們將 LongPhase-S 時間縮短至 20-30 分鐘，InterSubMod 分析時間縮短至 <5 分鐘。順序執行 (Sequential) 策略確保了每個步驟都能最大化利用 CPU。
3.  **過濾器效能**: 目前的甲基化過濾器 (Tier 1) 對於高品質/低 FP 的樣本效果有限。在 TP:FP 比例極高的情況下 (如 H1437 的 67468:8)，過濾器的風險(誤刪 TP)遠大於收益(移除 FP)。
4.  **未來方向**:
    - 針對 COLO829 這種 FP 較多的樣本，可以嘗試放寬過濾條件或引入更多特徵。
    - 調查 HCC1937 的低 Recall 原因。
    - **全面轉向 Pileup 模式**: 初步結果顯示 Pileup 模式能顯著提升 Recall 與整體 F1，值得在所有樣本上推廣。

## 5. 結論
本研究成功建立了多樣本評估流程。Output 模式下的過濾空間有限，但引入 Pileup 模式後，我們發現了提升 Recall 的新機會。HCC1395 的 Pileup F1 (0.8551) 超越了 Output 模式 (0.8522)，驗證了此策略的有效性。後續將完成全樣本 Pileup 驗證以確認此趨勢。
