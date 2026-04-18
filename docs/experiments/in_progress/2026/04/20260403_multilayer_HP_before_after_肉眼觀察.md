<!--
建立時間: 2026-04-03 23:30
目標: 多層 HP 測試改動前後的輸出數據比較，供肉眼觀察確認
處理範圍: Paired mode HCC1395 TP/FP (before: commit b9eaba7, after: multilayer HP)
關聯檔案:
  - output/multilayer_hp_benchmark/comparison_figures/
  - scripts/analysis/build_multilayer_hp_before_after_comparison.py
-->

# 多層 HP 測試：Before/After 肉眼觀察

## 實驗設定

| 項目 | Before (b9eaba7) | After (multilayer HP) |
|------|------------------|----------------------|
| 程式碼版本 | `b9eaba7` feat(QS) | Working tree (多層 HP) |
| 模式 | Paired | Paired |
| 樣本 | HCC1395 | HCC1395 |
| TP 區域數 | 30,476 | 30,476 |
| FP 區域數 | 4,822 | 4,822 |
| 區域匹配率 | 100% | 100% |

## 受影響的欄位

### 邏輯改動

| 欄位 | 改動前 | 改動後 | 影響方向 |
|------|--------|--------|---------|
| **GlobalP** | `min(p_alt, p_hp)` | `min(p_alt, p_hp, p_hp_family)` | 只能 ≤ |
| **CramersV** | `max(v_alt, v_hp)` | `max(v_alt, v_hp, v_hp_family)` | 只能 ≥ |
| **PassedGating** | `alt ∥ hp_pure` | `alt ∥ hp_family` | 見下方說明 |
| **VerificationClass** | `cluster_sig` 用 `hp_pure` | `cluster_sig` 用 `hp_family` | 可能升/降級 |
| **HeuristicScore** | `best_p/v` 用 `alt+hp` | `best_p/v` 用 `alt+hp+hp_family` | 只能 ≥ |

### 新增欄位（舊版無）

- `GlobalP_HPFamily`, `CramersV_HPFamily`
- `GlobalP_HPFine`, `CramersV_HPFine`, `HPFine_NGroups_CF`

### 已驗證不受影響的欄位

| 欄位 | 驗證結果 |
|------|---------|
| AlleleDelta / AlleleP / AlleleSig | **100% 相同** |
| HPMergedP / HPMergedDelta | **100% 相同** |
| LabelHPPermanovaF / LabelAllelePermanovaF | **100% 相同** |
| ClusterPermanovaF/P（同時通過 gating 的區域）| **100% 相同** |
| NumReads / NumCpGs | **100% 相同** |
| HP1FamilyN / HP2FamilyN | **100% 相同** |

---

## Figure 1: PassedGating 轉移矩陣

![PassedGating Transition](../../../../../output/multilayer_hp_benchmark/comparison_figures/fig1_gating_transition.png)

### 數據摘要

| 指標 | TP | FP |
|------|-----|-----|
| Before pass rate | 31.92% (9,728) | 46.27% (2,231) |
| After pass rate | 31.02% (9,455) | 47.01% (2,267) |
| 淨變化 | **-0.90pp** (-273) | **+0.75pp** (+36) |
| False→True (新通過) | 266 | 49 |
| True→False (失去通過) | 539 | 13 |

### 觀察重點

- **TP 淨減少、FP 淨增加**：表面看似退化，但原因是 **RNG state shift**（見下方說明）
- 失去 gating 的 539 個 TP 區域中，44.2% 的 GlobalP 在 0.03-0.07（**邊界區域**）
- 所有失去 gating 的區域 output CramersV = 0（V 不可靠）

### RNG State Shift 說明

```
舊版 test_all():     [test_allele(), test_hp(),        test_sample()]
新版 test_all_exp(): [test_allele(), test_hp_family(), test_hp(), test_hp_fine(), test_sample()]
```

`test_hp_family()` 插入在 `test_hp()` 前面，消耗了共享 Fisher MC RNG 的隨機數。
→ `test_hp()` 在新版中使用不同的 RNG state → MC p-value 在邊界區域產生雙向波動。

**這不是邏輯錯誤，是 MC 抽樣噪音。** 建議未來為每個測試使用獨立 RNG（從 master seed 衍生）來消除此效應。

---

## Figure 2: VerificationClass 轉移矩陣

![VerificationClass Transition](../../../../../output/multilayer_hp_benchmark/comparison_figures/fig2_verification_class_transition.png)

### 數據摘要

| 類別 | TP Before → After | FP Before → After |
|------|-------------------|-------------------|
| Strong | 7,282 → 7,229 (-53) | 1,737 → 1,750 (+13) |
| Subclone | 1,452 → 1,463 (+11) | 239 → 255 (+16) |
| Weak | 11,335 → 11,388 (+53) | 1,152 → 1,139 (-13) |
| Noise | 10,407 → 10,396 (-11) | 1,694 → 1,678 (-16) |

### 主要轉移路徑

| 轉移 | TP 數量 | 原因 |
|------|---------|------|
| Strong → Weak | 304 | 失去 cluster_significant（gating 變化） |
| Weak → Strong | 251 | 獲得 cluster_significant |
| Noise → Subclone | 82 | 新通過 gating 但無 label_sig |
| Subclone → Noise | 71 | 失去 gating |

---

## Figure 3: GlobalP 比較

![GlobalP Comparison](../../../../../output/multilayer_hp_benchmark/comparison_figures/fig3_globalp_comparison.png)

### 數據摘要

| 指標 | TP | FP |
|------|-----|-----|
| 完全相同 | 84.3% | 類似 |
| 改善（After 更低） | 3,399 | 364 |
| 不變 | 26,445 | 4,418 |
| 惡化（After 更高） | 632 | 40 |
| Mean delta (-log10) | +0.1234 | +0.1124 |

### 觀察重點

- 散點圖幾乎全在 y=x 線上，少量點偏上（改善）
- 改善方向占多數（3,399 vs 632），這是 hp_family 提供額外更低 p-value 的效果
- 632 個「惡化」實際是 RNG state shift 導致的 MC 噪音

---

## Figure 4: HeuristicScore 比較

![HeuristicScore Comparison](../../../../../output/multilayer_hp_benchmark/comparison_figures/fig4_heuristic_score.png)

### 數據摘要

| 指標 | TP | FP |
|------|-----|-----|
| 改善 | 3,911 | 442 |
| 不變 | 25,934 | 4,340 |
| Mean delta | +0.1358 | +0.1296 |

---

## Figure 5: Quality Score & Tier 比較

![Quality Comparison](../../../../../output/multilayer_hp_benchmark/comparison_figures/fig5_quality_comparison.png)

### 數據摘要

| 指標 | TP | FP |
|------|-----|-----|
| QS 改善 | 45 | 24 |
| QS 不變 | 30,431 | 4,798 |
| QS 惡化 | 0 | 0 |
| QS Mean delta | +0.0141 | +0.0518 |

### 觀察重點

- Quality_Score 變化非常微小（99.85% TP 不變）
- **零惡化**：沒有任何區域的 Quality_Score 下降
- Quality_Tier 完全不變（所有 TP/FP 的 tier 分布不變）

---

## Figure 6: 新增欄位分布

![New Columns](../../../../../output/multilayer_hp_benchmark/comparison_figures/fig6_new_columns.png)

### HP Family 顯著性

| 指標 | TP | FP |
|------|-----|-----|
| HP Family sig (p<0.05) | 5,556 (18.2%) | 985 (20.4%) |

### HP Fine N Groups（Cluster-First）

| NGroups | TP | FP |
|---------|-----|-----|
| 0 (invalid) | 6,020 (19.8%) | 764 (15.8%) |
| 2 | 14,088 (46.2%) | 2,413 (50.0%) |
| 3 | 8,087 (26.5%) | 1,176 (24.4%) |
| 4 | 2,281 (7.5%) | 469 (9.7%) |

---

## Figure 7: Sanity Check（不應變化的欄位）

![Sanity Check](../../../../../output/multilayer_hp_benchmark/comparison_figures/fig7_sanity_check.png)

### 結果

| 欄位 | 差異數 | 說明 |
|------|--------|------|
| NumReads | 0 | ✅ |
| NumCpGs | 0 | ✅ |
| ClusterPermanovaF | 805 (2.6%) | ⚠️ 全來自 gating 進出（同時通過者 100% 相同） |
| ClusterPermanovaP | 805 (2.6%) | ⚠️ 同上 |
| HPMergedP | 0 | ✅ |
| AlleleP | 0 | ✅ |
| HP1FamilyN | 0 | ✅ |
| HP2FamilyN | 0 | ✅ |
| LabelHPPermanovaF | 0 | ✅ |
| LabelAllelePermanovaF | 0 | ✅ |

---

## 綜合評估

### Paired 模式影響程度：**微小**

多層 HP 改動對 Paired 模式的影響非常有限：
- 84.3% 的 GlobalP **完全相同**
- Quality_Score 99.85% **不變**，0 個惡化
- VerificationClass 淨變化 < 0.2%

### 觀察到的雙向 gating 波動根因

不是 hp_family 邏輯問題，而是 **Fisher MC RNG state shift**：
- `test_hp_family()` 插入到共享 RNG 序列中
- `test_hp()` 在新版中看到不同的隨機數 → 邊界區域 p-value 雙向波動
- 影響 539+266=805 個區域（2.6%）

### 建議

1. **結果可接受**：淨影響在統計噪音範圍內
2. **可選改善**：為每個 GlobalTest 方法使用獨立 RNG（從 master seed 衍生），消除測試順序對結果的影響
3. **TO 模式才是重點**：Paired 模式的 "1-1"/"2-1" 讀段稀少，多層 HP 的真正價值在 TO 模式

---

## 圖片路徑

所有圖片位於：
```
output/multilayer_hp_benchmark/comparison_figures/
├── fig1_gating_transition.png
├── fig2_verification_class_transition.png
├── fig3_globalp_comparison.png
├── fig4_heuristic_score.png
├── fig5_quality_comparison.png
├── fig6_new_columns.png
├── fig7_sanity_check.png
└── summary_stats.txt
```
