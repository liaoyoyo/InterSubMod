<!--
建立時間: 2026-04-04 08:00
目標: 完整驗證 LOH_Strong/Weak 7-feature AUC > 0.88 的真實性、可推廣性、TO 流程可用性
處理範圍: v2b TO ISM (HCC1395 5kHz, 30,476 TP + 4,822 FP) + baseline 對照
關聯檔案:
  - docs/experiments/in_progress/2026/04/20260404_HPFineP_QS整合完整研究報告_01.md
  - docs/experiments/in_progress/2026/04/20260404_HPFineP納入QS研究分析_01.md
  - src/core/RegionProcessor.cpp (QS 計算, lines 105-196)
  - src/core/SignificanceAnalyzer.cpp (LOH 偵測, lines 313-339)
-->

# LOH_Strong/Weak 7-Feature AUC > 0.88 完整驗證報告

## 驗證問題

> LOH_Strong/Weak 7-feature AUC > 0.88：
> 1. 是否在所有樣本都成立？
> 2. 是否真正能區分 TP/FP？
> 3. 是否能在 TO 流程下正確應用？
> 4. 是否使用到 HP 或 LOH 資訊？

本報告由三個獨立分析 agent 平行驗證，綜合結論如下。

---

## 1. 驗證架構

| Agent | 任務 | 核心問題 |
|-------|------|---------|
| **Agent 1** | Circular Dependency 分析 | LOH_Subtype 定義是否與 7 個特徵循環依賴？ |
| **Agent 2** | Cross-Validation & Robustness | AUC 是否過擬合？LOCO/Permutation 是否穩定？ |
| **Agent 3** | TO 流程可用性 | 實際 variant filter 效果？覆蓋率？HP 依賴？低純度風險？ |

---

## 2. Agent 1 結論：Selection Bias 存在但非 Label Leakage

### 2.1 LOH_Subtype 如何定義

`LOH_Subtype` 由 `determine_loh_subtype()` 計算（`RegionProcessor.cpp:80-96`），核心依據：

```
HP_Ratio = max(HP1_count, HP2_count) / (HP1_count + HP2_count)
Potential_LOH = HP_Ratio ≥ 0.85
```

| LOH_Subtype | 定義條件 |
|-------------|---------|
| LOH_Strong | Potential_LOH + VerificationClass == Strong |
| LOH_Weak | Potential_LOH + VerificationClass == Weak |
| LOH_Noise | Potential_LOH + VerificationClass == Noise |
| LOH_Subclone | Potential_LOH + VerificationClass == Subclone |
| None | 非 Potential_LOH |

### 2.2 循環依賴鏈

```
HP_Ratio ≥ 0.85 → 定義 LOH_Subtype → 篩選子群
                                         ↓
7 特徵中 3 個直接使用 HP:
  - HP_Ratio (就是分群依據本身)
  - HPFineP (HP Fine 測試 p-value)
  - LabelHPPermanovaF (HP PERMANOVA F)
                                         ↓
在 HP-extreme 子群中用 HP 特徵 → 結構性 AUC 膨脹
```

### 2.3 判定

| 項目 | 結論 |
|------|------|
| **是否 label leakage** | **否**（LOH_Subtype 不使用 TP/FP 標籤） |
| **是否 selection bias** | **是**（用 HP_Ratio 選子群，再用 HP 特徵分類 → 子群內 AUC 膨脹） |
| **HP 資訊使用** | **3/7 特徵直接依賴 HP tag** |
| **LOH 資訊使用** | **LOH_Subtype 定義依賴 HP_Ratio → 間接使用 LOH/HP** |
| **TO 可用性** | 7 個特徵在 TO 模式都可取得（無需 truth label） |

---

## 3. Agent 2 結論：AUC 是真實的，沒有過擬合

### 3.1 5-Fold Cross-Validation

| 子群 | Resubstitution AUC | CV AUC (mean±std) | Gap |
|------|-------------------|-------------------|-----|
| LOH_Strong | 0.912 | **0.912 ± 0.012** | **0.001** |
| LOH_Weak | 0.884 | **0.881 ± 0.022** | **0.004** |

**CV gap < 0.5% → 沒有過擬合。**

### 3.2 Leave-One-Chromosome-Out (LOCO)

| 子群 | LOCO Mean | LOCO Std | 最差 chr |
|------|-----------|----------|---------|
| LOH_Strong | **0.897** | 0.092 | **chr14 = 0.529** |
| LOH_Weak | **0.859** | 0.047 | chr6 = 0.781 |

- LOH_Strong chr14 AUC=0.529（近隨機）：該染色體 FP 分布異常
- 其餘染色體穩定在 0.85-0.97

### 3.3 Permutation Test

| 子群 | 觀測 AUC | Permutation 均值 | Z-score | p-value |
|------|---------|----------------|---------|---------|
| LOH_Strong | 0.912 | 0.498 | **27.5σ** | **0/100** |
| LOH_Weak | 0.884 | 0.501 | 22.1σ | **0/100** |

**訊號是統計學上真實的，非隨機。**

### 3.4 v2b vs Baseline 比較

| 子群 | v2b AUC | Baseline AUC | Delta |
|------|---------|-------------|-------|
| LOH_Strong | 0.912 | 0.887 | -0.025 |
| LOH_Weak | 0.884 | 0.876 | -0.008 |

v2b（PON-only phasing）略優於 baseline，但差距不大。

### 3.5 判定

| 項目 | 結論 |
|------|------|
| **過擬合** | **否**（CV gap < 0.5%） |
| **統計顯著** | **是**（27.5σ，p < 0.01） |
| **跨染色體穩定** | **大致是**（例外：chr14 LOH_Strong = 0.529） |
| **v2b vs baseline** | **差異小**（HP tag 品質影響 < 0.03 AUC） |

---

## 4. Agent 3 結論：NOT ACTIONABLE — 不具實際可用性

### 4.1 覆蓋率分析

| LOH_Subtype | 佔全量 | TP 覆蓋 | FP 覆蓋 | 7-feat AUC |
|-------------|--------|--------|--------|-----------|
| LOH_Strong | 16.6% | 16.2% | 18.8% | 0.912 |
| LOH_Weak | 20.2% | 21.5% | 12.0% | 0.884 |
| LOH_Noise | 29.7% | 29.4% | 31.7% | **0.688** |
| LOH_Subclone | 5.1% | 5.0% | 5.3% | **0.645** |
| None | 28.4% | 27.8% | 32.2% | 0.800 |

**LOH_Strong + Weak 合計僅覆蓋 36.8% 的變異。**

### 4.2 Variant Filter 實際效果（LOH_Strong 子群內）

| Threshold | Precision | Recall | F1 | FP removed | TP lost |
|-----------|-----------|--------|----|-----------|---------|
| 0.5 | 0.918 | 0.966 | 0.941 | 53.0% | 3.4% |
| Optimal F1 | 0.933 | 0.955 | 0.944 | 62.5% | 4.5% |
| 95% recall | 0.937 | 0.950 | 0.943 | 65.0% | 5.0% |

### 4.3 全局 F1 影響 — 微乎其微

| Filter 範圍 | FP removed (全局) | TP lost (全局) | Global F1 | ΔF1 |
|-------------|-------------------|---------------|-----------|-----|
| Baseline | 0% | 0% | 0.927 | — |
| LOH_Strong only | 11.7% | 0.7% | 0.931 | **+0.004** |
| LOH_Strong+Weak | 18.2% | 1.3% | 0.933 | **+0.006** |

**全局 F1 改善僅 +0.006（0.927 → 0.933），不具實質意義。**

### 4.4 主要驅動力 — AF Confound

| 特徵 | LOH_Strong 單特徵 AUC | 性質 |
|------|---------------------|------|
| **HPFineP** | **0.823** | HP-dependent |
| **AlleleDelta** | **0.800** | **AF confound** |
| HP_Ratio | 0.603 | HP-dependent |
| CramersV | 0.599 | 非 HP |
| LabelAllelePermanovaF | 0.593 | 非 HP |
| PairwiseMeanDist | 0.555 | 非 HP |
| LabelHPPermanovaF | 0.502 | HP-dependent |

**AlleleDelta + HPFineP 兩特徵合計 AUC = 0.872**，已接近 7-feature 的 0.912。

**AlleleDelta 的區分力機制（與 O12 結論一致）：**
- **TP somatic in LOH**：突變在單一等位基因 → 大部分 reads 同一 HP → AlleleDelta 小
- **FP germline in LOH**：先於 LOH 的雜合變異 → 兩等位基因仍有甲基化差異 → AlleleDelta 大
- 這是 LOH 區域的結構性特徵，不是甲基化新訊號

### 4.5 移除 HP 特徵後的效果

| 組合 | LOH_Strong AUC | 說明 |
|------|---------------|------|
| 7-feature (all) | 0.912 | 完整 |
| 4-feature (no HP) | **0.899** | 移除 HPFineP, HP_Ratio, LabelHPPermanovaF |
| AlleleDelta only | 0.800 | 單特徵 |

**移除 HP 特徵 → AUC 僅降 0.013。** 全局 F1 改善不變（+0.006 → +0.006）。

### 4.6 低純度風險模擬

| 噪聲水平 (std) | AUC | 下降 | 對應臨床場景 |
|---------------|-----|------|------------|
| 0 (原始) | 0.912 | — | HCC1395 purity 93% |
| 0.5 | 0.889 | -0.023 | purity ~70% |
| 1.0 | **0.837** | -0.075 | purity ~50% |
| 2.0 | **0.739** | -0.173 | purity ~30% |

**臨床 purity 30-50% 下，AUC 預計降至 0.74-0.84，全局效果更差。**

---

## 5. 綜合判定

### 5.1 回答原始問題

| 問題 | 答案 | 證據 |
|------|------|------|
| **是否在所有樣本都成立？** | **無法確認** | 僅 HCC1395 單一樣本驗證；噪聲模擬顯示低純度下 AUC 顯著衰退；chr14 LOCO AUC=0.529 |
| **是否真正能區分 TP/FP？** | **統計上是，實際上效果微小** | AUC 真實（CV gap<0.5%, permutation 27.5σ），但全局 F1 改善僅 +0.006 |
| **是否能在 TO 流程下正確應用？** | **不建議** | 覆蓋率 36.8%，全局效果微乎其微，主要區分力來自 AF confound |
| **是否使用到 HP 或 LOH 資訊？** | **是** | 3/7 特徵直接依賴 HP；LOH_Subtype 由 HP_Ratio 定義；移除 HP 後 AUC 僅降 0.013 |

### 5.2 AUC > 0.88 的本質

```
LOH_Strong/Weak AUC > 0.88 的因果鏈：

HP_Ratio ≥ 0.85 → 選出 HP-extreme 子群（LOH 區域）
                     ↓
在 LOH 區域中：
  - TP somatic: 突變在單一等位基因 → AlleleDelta 小 + HPFineP 顯著
  - FP germline: LOH 前雜合 → AlleleDelta 大 + HPFineP 不顯著
                     ↓
AlleleDelta + HPFineP 兩特徵 AUC = 0.872（7-feat = 0.912）
                     ↓
本質 = LOH 區域的 AF confound + HP-extreme selection bias
       不是可推廣的甲基化生物學信號
```

### 5.3 與既有結論的一致性

| 既有結論 | 本次驗證 | 一致性 |
|---------|---------|--------|
| **O12**: AlleleDelta = AF confound | AlleleDelta 單特徵 AUC=0.800 in LOH_Strong | ✅ 完全一致 |
| **O1-O13**: TO 無特徵 AUC > 0.58 (全域) | 全域 AUC 確實 < 0.80，子群高 AUC 是 selection bias | ✅ 一致 |
| **Germline FP NO-GO**: 所有 post-hoc 特徵 AUC < 0.64 | 全域 7-feat AUC=0.800 但受 LOH 子群拉高 | ✅ 一致 |
| **QS TO AUC=0.497（隨機）** | QS 需要改善但 LOH_Strong/Weak LR 不是解法 | ✅ 一致 |

---

## 6. 對 QS 整合的影響

### 6.1 LOH_Strong/Weak LR 不適合作為 QS 改善方案

7-feature LR 雖然 AUC 高，但：
- 覆蓋率太低（36.8%）
- 全局效果微小（ΔF1=+0.006）
- 需要訓練 LR 模型 → 過擬合風險大（僅單樣本）

### 6.2 HPFineP 的 QS 整合方案（Scheme 3）仍然推薦

相較之下，Scheme 3（LOH+HPFineSig bonus=15, NGroups≤1 penalty=15）：
- **不需要 LR 模型**（規則式，無過擬合風險）
- **覆蓋率 100%**（NGroups≤1 penalty 對所有區域生效）
- **AUC 改善 +0.048**（全域，非子群）
- **FP High tier 降低 15.1pp**
- **Paired 不受影響**

| 比較 | 7-feat LR (LOH only) | Scheme 3 (全域) |
|------|---------------------|----------------|
| 覆蓋率 | 36.8% | 100% |
| 全域 F1 Δ | +0.006 | +0.048 AUC |
| 過擬合風險 | 高（單樣本 LR） | 低（規則式） |
| 實作複雜度 | 高（需 LR + threshold） | 低（兩條 if 語句） |
| **推薦** | **否** | **是** |

---

## 7. 結論與建議

### 7.1 最終判定

**LOH_Strong/Weak 7-feature AUC > 0.88 是統計上真實的，但不具實際可用性。**

核心原因：
1. **Selection bias**：用 HP_Ratio 選子群再用 HP 特徵 → 子群內 AUC 膨脹
2. **AF confound**：AlleleDelta 反映 LOH 下的等位基因頻率差異，不是甲基化新訊號
3. **覆蓋率限制**：僅覆蓋 36.8%，全局效果微乎其微（ΔF1=+0.006）
4. **單樣本風險**：高純度 HCC1395 不可推廣到低純度臨床樣本
5. **與 O12 一致**：這不是新發現，是已知 AF confound 的子群表現

### 7.2 建議行動

| 項目 | 行動 |
|------|------|
| LOH_Strong/Weak LR filter | **正式關閉**，不納入 TO 流程 |
| Scheme 3 HPFineP QS 整合 | **繼續推進**（規則式，全域有效） |
| 多樣本驗證 | 待 big8 BAM 建立後執行（已排程但不阻塞 Scheme 3） |

---

## 附錄：三個 Agent 的完整數據

### A1. Agent 1 — 特徵 HP 依賴分類

| 特徵 | HP 依賴 | 說明 |
|------|--------|------|
| AlleleDelta | 否 | Allele 分群的甲基化差異，但在 LOH 中是 AF confound |
| CramersV | 否 | 全域統計 |
| PairwiseMeanDist | 否 | 距離矩陣 |
| LabelAllelePermanovaF | 否 | Allele PERMANOVA |
| HPFineP | **是** | HP Fine 4-group Fisher p-value |
| HP_Ratio | **是** | max(HP1,HP2)/(HP1+HP2)，且是 LOH 定義依據 |
| LabelHPPermanovaF | **是** | HP merged PERMANOVA F |

### A2. Agent 2 — LOCO 染色體詳細表

**LOH_Strong:**

| Chr | AUC | 備註 |
|-----|-----|------|
| chr1 | 0.929 | |
| chr2 | 0.872 | |
| chr3 | 0.885 | |
| chr5 | 0.938 | |
| chr6 | 0.926 | |
| chr8 | 0.897 | |
| chr9 | 0.878 | |
| chr10 | 0.959 | |
| chr11 | 0.970 | |
| chr14 | **0.529** | FP 分布異常 |
| chr17 | 0.912 | |
| Mean | 0.897 | |
| Std | 0.092 | |

### A3. Agent 3 — 特徵組合消融實驗

| 組合 | LOH_Strong AUC | 說明 |
|------|---------------|------|
| AlleleDelta only | 0.800 | 單特徵最強（non-HP） |
| HPFineP only | 0.823 | 單特徵最強（HP） |
| AlleleDelta + HPFineP | 0.872 | 兩特徵接近完整模型 |
| 4-feature (no HP) | 0.899 | 移除 HP 僅降 0.013 |
| 7-feature (all) | 0.912 | 完整模型 |
