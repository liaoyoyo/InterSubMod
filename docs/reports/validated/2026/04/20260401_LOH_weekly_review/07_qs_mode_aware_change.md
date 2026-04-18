# 07. QS Mode-Aware 程式碼修改

<!--
建立時間: 2026-04-01 00:00
目標: 記錄 QualityScore mode-aware 重構的動機、修改內容、預期效果與佐證資料
處理範圍: src/core/RegionProcessor.cpp commit b9eaba7 — QualityScoreWeights struct 設計
關聯檔案:
  - src/core/RegionProcessor.cpp
  - docs/reports/validated/2026/03/20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md
  - docs/experiments/INDEX.md (O2 QS component decomposition)
-->

---

## 1. 修改動機

### 1.1 TO QS 等同隨機猜測

| 指標 | Paired | Tumor-Only |
|------|--------|------------|
| QS AUC | 0.754 | **0.497** |
| 隨機基準 | 0.500 | 0.500 |
| 解讀 | 有效 | **完全無效** |

TO 模式的 QS AUC = 0.497，甚至略低於隨機猜測（0.500），代表 QS 對 TO 場景不僅無幫助，反而有微幅反向效果。

### 1.2 根因分析：LOH penalty 打擊 TP > FP

O2 觀察（QS component waterfall decomposition）揭示了失效根因：

| QS Component | TP trigger rate | FP trigger rate | 方向 |
|-------------|-----------------|-----------------|------|
| **LOH penalty (-25)** | **44.5%** | **35.8%** | **反向**（懲罰 TP 多於 FP） |
| Verify bonus (+10) | 較低 | 較高 | **反向**（獎勵 FP 多於 TP） |

- **LOH penalty**：TO 模式下 phasing artifact 導致 LOH-like 比率虛高，而 TP 位點因為真實 allelic imbalance 更容易被判定為 LOH，因此 penalty 打擊 TP 的比例（44.5%）高於 FP（35.8%）。
- **Verify bonus**：VerificationClass 在 TO 下 Cramer's V = 0.023（近乎無效），但 bonus 仍然基於 HP-sig 和 Allele-sig 判定。由於 TO phasing 雜訊使 FP 更容易同時觸發兩項顯著性檢驗，bonus 實際上獎勵 FP 多於 TP。

### 1.3 關鍵數據佐證

- TO 同位點 LOH 判定率為 paired 的 16-52 倍（系統性過判）
- VerificationClass Cramer's V = 0.023（閾值 0.1 以下視為無效）
- 所有 TO 單一特徵 AUC < 0.58

---

## 2. 修改內容

### 2.1 QualityScoreWeights struct 設計

Commit `b9eaba7` 將原本硬編碼的 QS 計算重構為 `QualityScoreWeights` struct，支援 mode-aware 權重配置：

```cpp
struct QualityScoreWeights {
    float read_penalty_severe = 20.0f;   // num_reads < 30
    float read_penalty_moderate = 10.0f; // num_reads < 50
    float cpg_penalty_severe = 15.0f;    // num_cpgs < 20
    float cpg_penalty_moderate = 10.0f;  // num_cpgs < 30
    float cov_penalty_cnv_loss = 30.0f;  // coverage_multiple < 0.5
    float cov_penalty_low = 15.0f;       // coverage_multiple < 0.8
    float cov_penalty_high_copy = 20.0f; // coverage_multiple > 2.0
    float loh_penalty = 25.0f;           // potential LOH
    float verify_bonus = 10.0f;          // HP and Allele tests both significant
    float effect_bonus_strong = 15.0f;   // cramers_v >= 0.3
    float effect_bonus_moderate = 5.0f;  // cramers_v >= 0.2
};
```

### 2.2 Paired vs TO 權重配置

| Component | Paired weights | TO weights | 修改理由 |
|-----------|---------------|------------|---------|
| read_penalty_severe | 20.0 | 20.0 | 不變 — read count 在兩種模式下同等重要 |
| read_penalty_moderate | 10.0 | 10.0 | 不變 |
| cpg_penalty_severe | 15.0 | 15.0 | 不變 — CpG count 在兩種模式下同等重要 |
| cpg_penalty_moderate | 10.0 | 10.0 | 不變 |
| cov_penalty_cnv_loss | 30.0 | 30.0 | 不變 — coverage anomaly 在兩種模式下同等重要 |
| cov_penalty_low | 15.0 | 15.0 | 不變 |
| cov_penalty_high_copy | 20.0 | 20.0 | 不變 |
| **loh_penalty** | **25.0** | **0.0** | **TO LOH-like 由 phasing artifact 主導，penalty 方向反轉** |
| **verify_bonus** | **10.0** | **0.0** | **TO VerificationClass 無效（V=0.023），bonus 獎勵 FP > TP** |
| effect_bonus_strong | 15.0 | 15.0 | 不變 — effect size 仍有參考價值 |
| effect_bonus_moderate | 5.0 | 5.0 | 不變 |

### 2.3 Mode detection 邏輯

```cpp
auto qs_weights = normal_bam_path_.empty() ? get_tumor_only_weights() : get_paired_weights();
```

透過 `normal_bam_path_` 是否為空判斷模式：
- 有 normal BAM → Paired → 使用預設權重（不變）
- 無 normal BAM → Tumor-Only → LOH penalty = 0, verify bonus = 0

### 2.4 改動範圍

- **檔案**：`src/core/RegionProcessor.cpp`（1 file changed, 48 insertions, 21 deletions）
- **影響範圍**：僅 QS 計算邏輯，不影響其他模組
- **向後相容**：Paired 模式行為完全不變

---

## 3. 預期效果

### 3.1 AUC 改善預估

| 指標 | 修改前 | 修改後 | Delta |
|------|-------|--------|-------|
| TO QS AUC | 0.497 | ~0.546 | **+0.049** |
| Paired QS AUC | 0.754 | 0.754 | 0.000（不變） |

### 3.2 效果有限的原因

移除 LOH penalty 和 verify bonus 後，TO QS 仍然很低（~0.546），因為 QS 的其他 components 在 TO 也弱：

- **Coverage penalty**：TO 和 Paired 的 coverage 分佈無顯著差異，penalty 對 TP/FP 區分力弱
- **Effect size bonus**：TO 的 Cramer's V 分佈在 TP/FP 之間重疊度高
- **Read/CpG penalty**：與 TP/FP 無關，僅反映 region 品質

**結論**：此修改是「移除幫倒忙的 component」，而非「增加有效 component」。TO QS 要達到實用水準（AUC > 0.7），需要全新的 feature engineering（Phase 1A ML 方向）。

### 3.3 實際意義

- 消除了 QS 在 TO 模式下的**反向效果**
- 為後續 Paired/TO 分離策略打下基礎
- 建立了 mode-aware 架構，未來新增 TO-specific features 更方便

---

## 4. O2 圖片佐證

### 4.1 QS Component Waterfall — Paired

![Paired QS component waterfall](figures/fig01_qs_component_waterfall_paired.png)

Paired 模式下各 component 對 TP/FP 的貢獻方向一致：LOH penalty 正確懲罰 FP，verify bonus 正確獎勵 TP。

### 4.2 QS Component Waterfall — TO

![TO QS component waterfall](figures/fig02_qs_component_waterfall_to.png)

TO 模式下 LOH penalty 和 verify bonus 的方向反轉，是 QS 失效的直接原因。

### 4.3 LOH Penalty Trigger Rate

![LOH penalty trigger rate comparison](figures/fig03_qs_loh_penalty_trigger_rate.png)

TO 模式下 LOH penalty 在 TP 的觸發率（44.5%）高於 FP（35.8%），驗證了 penalty 方向反轉的結論。

---

## 5. Commit 資訊

- **Hash**: `b9eaba7844faa2ac05fed774c4b0e43c570e002b`
- **Message**: `feat(QS): mode-aware quality score -- disable LOH penalty and verify bonus for tumor-only`
- **Date**: 2026-03-31
- **Files changed**: `src/core/RegionProcessor.cpp` (+48, -21)

---

## 待驗證問題（已驗證 / 已更新）

### 已解決

2. **~~移除 LOH penalty 後，是否有 paired-like TO 樣本受損？~~** ✅ **不會。** O2 分析確認 TO 全域 QS AUC = 0.497（隨機水準），LOH penalty 在 TO 是反向作用（AUC=0.457 < 0.5），移除只能改善不能損害。且 C13 確認方向翻轉是 7/7 樣本一致的系統性屬性，不存在「paired-like TO 樣本」。

### 尚未解決

1. **修改後的 TO QS AUC 是否確實達到 ~0.546？** 需實際跑 benchmark 驗證。
3. **Effect size bonus 在 TO 下是否需要調整？** 目前保持不變，需 benchmark 數據評估。AUC=0.505（幾乎無效），調整的預期增益極小。

## 認知門檻補充建議

- 對於不熟悉 QS 計算的讀者，建議先閱讀 O2 報告（QS component decomposition）了解各 component 的作用機制
- LOH penalty 方向反轉是 TO phasing artifact 的下游效應；根本原因在於 LongPhase-TO 缺乏 normal BAM �� phasing anchor
- Mode-aware 架構是通往 Paired/TO 完全分離模型的第一步
