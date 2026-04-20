# Multiple Testing Correction Audit（全 AUC 掃描 FDR 補算）

> **建立日期**: 2026-04-19

## 背景

**問題來源**（C-STAT-4，Phase 2 Role 1 審查）：
> 全 AUC 掃描（60+ 特徵）無 FDR/Bonferroni 校正 → 補充結論假陽性風險

**用戶決定（P1-B）**：對 60+ 特徵套 BH-FDR，補充列表

**涉及範圍**：
- 14+ 主要結論 + 補充結論 15-22
- 系統性觀察 O1-O13
- Zone-Aware Framework Z1-Z5 × 多指標
- B.2 系列 subclone marker 探索

---

## 需套 FDR 的結論清單

| Card | 結論主題 | 特徵數量 | AUC 範圍 | POSITIVE 門檻 | FDR 優先序 |
|------|---------|---------|---------|-------------|-----------|
| **C03** | TO AUC ceiling | 30+ | <0.58 全 | N/A（全 NEGATIVE） | P2（驗證 ceiling 穩定） |
| **C04** | O11 Heterogeneity | ~10 | 0.6-0.7 | n_reads stratify 後消失 | P1（已 NEGATIVE 但應補 FDR） |
| **C05** | O12 LOH scenarios | ~15 | 0.6-0.75 | L2 AF bin 後消失 | P1（已 NEGATIVE） |
| **C06** | O13 Cross-region | ~8 | 0.6-0.7 | shared read stratify 後消失 | P1（已 NEGATIVE） |
| **C07** | G1-G7 Germline | 7 | <0.64 全 | LOSO AUC=0.721 但 FP=0 | P2 |
| **C11** | Phase 1A F1 | 3 ML models | F1 metric | +0.0112（小） | P1（effect size 小） |
| **C16** | HPFineNGroups NG × NR | ~20 組合 | residualized | NG=4+AF<0.4+NR≥80 | **P0（最需 FDR）** |
| **C17** | LOH Subclone AF×NGroups | r=+0.705 + cross-sample | 7/7 | p<10⁻³⁹ | P1（單一 metric） |
| **C22** | Zone × Multi-metric | 5 zones × 多指標 | Zone TP rate | N/A | P1（bug-fix 後） |

---

## 最高優先：C16 HPFineNGroups FDR

### 問題

- **特徵掃描空間**：NG=[2,3,4,5,≥6] × NR cut-off [40,60,80,100,120,150] × AF filter [none, <0.4, <0.5] ≈ 5×6×3 = **90 組合**
- 最終選定 NG=4+AF<0.4+NR≥80 為 best filter
- 單點 AUC 聲明無 FDR → 「最佳組合」可能是多重測試的 artifact

### FDR 補算方案（BH-FDR）

```python
# 虛擬碼
combos = expand_grid(NG_list, NR_list, AF_list)
for combo in combos:
    auc[combo] = compute_auc(data, combo)
    p_values[combo] = permutation_test(data, combo, n_perm=1000)

adjusted_p = benjamini_hochberg(p_values, alpha=0.05)
surviving_combos = [c for c in combos if adjusted_p[c] < 0.05]
```

### 驗收

- **穩固**：NG=4+AF<0.4+NR≥80 survives BH-FDR（adjusted p<0.05）
- **反轉**：最佳組合不 survives → C16 結論需降為 CONDITIONAL
- **部分穩固**：多個組合 survive → 報告 all surviving combos，NG=4 非唯一

---

## O11/O12/O13 FDR 補算（C04/C05/C06）

**雖然結論已 NEGATIVE（歸因於 confound），但：**
- 原本 AUC 0.6-0.7 的觀察仍散佈在文件中
- 應追溯補 FDR，正式關閉這些「pseudo-POSITIVE」的 observational AUC

### 方案

1. 列出 O11/O12/O13 各自的所有 AUC 數值（歷史報告中）
2. 套 BH-FDR 於每組（分別對 O11/O12/O13 獨立校正）
3. 報告 adjusted p value
4. 在 C04/C05/C06 audit card 註明「FDR 補算後 X 個特徵 survive，其中 Y 個經 stratify 後消失」

### 驗收

- 預期：多數特徵 FDR-adjusted p 仍 <0.05（pre-stratify），但 stratify 後效應消失
- 這強化 NEGATIVE 結論：「統計顯著 ≠ 真實 signal，confound 才是根因」

---

## Phase 1A F1 Effect Size 揭露（C11）

**問題**（C-STRAT-1）：ΔF1=+0.0112 過小，未揭露 per-sample CI overlap

### 補算方案

1. Per-sample bootstrap 1000× F1
2. 計算 per-sample F1 95% CI
3. 畫 cross-sample CI overlap 圖
4. 驗收：若 CI overlap 廣 → 實務顯著性存疑，標 CONDITIONAL

---

## Zone-Aware Multi-metric FDR（C22）

**問題**：5 zones × 多指標（TP rate, F1, ΔF1, QS, methylation）= 15-25 組測試

### 方案

1. CovM bug fix 後重跑 Zone-Aware 分析
2. 每個 metric 套 BH-FDR
3. 報告 per-zone per-metric adjusted p
4. 驗收：Z3 TP rate 0.608 vs Z1 0.85 在 FDR 後仍顯著 → characterization 結論穩固

---

## LOH Subclone Cross-Sample Meta p（C17）

**現狀**：7/7 樣本 p<10⁻³⁹，無 FDR 必要（單一 hypothesis per sample）

**但需做**：
- Bootstrap CI for r=+0.705（P1-C 已列）
- Meta-analysis across 7 samples（Fisher's method）
- 驗收：meta p 仍顯著 → 結論不因多樣本增加 Type I 風險

**無需 BH-FDR**（單 hypothesis × 7 獨立樣本 ≠ 多重測試）

---

## G1-G7 Germline FP FDR（C07）

**問題**：G1-G7 7 個特徵 + LOSO AUC = 8 個測試

### 方案

1. 對 G1-G7 AUC 套 Bonferroni（n=7，嚴格但清晰）
2. LOSO AUC 0.721 報告 95% CI（bootstrap）
3. 驗收：LOSO CI lower bound >0.65 → LOSO 結果穩固（儘管 FP=0）

---

## FDR 補算整體執行順序

### Phase A（與 CovM bug fix 並行）

1. **C16 HPFineNGroups BH-FDR**（90 組合）
   - 預期 2-3 天
   - 結果寫回 C16 audit card + 00_INDEX.md

### Phase B（CovM bug fix 後）

2. **C22 Zone-Aware multi-metric FDR**（15-25 測試）
3. **C17 meta p-value + bootstrap CI**

### Phase C（追溯驗證 NEGATIVE 結論）

4. **O11/O12/O13 FDR 補算**（歷史 observational AUC）
   - 文件化 NEGATIVE 結論不因 FDR 改變
5. **G1-G7 Bonferroni + LOSO CI**

### Phase D（Effect size 揭露）

6. **C11 Phase 1A F1 per-sample bootstrap**

---

## 不做 FDR 的項目（止損）

以下結論**不適用** FDR：
- **C01**（Paired/TO FP rate 分離，描述性統計）
- **C02**（PON coverage metric）
- **C08**（LOSO single AUC，用 bootstrap 不用 FDR）
- **C09 / C10 / C14**（causal chain / mechanism，非 feature screening）
- **C12 ASM**（分佈描述）
- **C13 / C21 LOH.bed**（set operation Jaccard）
- **C15 / C18 / C19 / C20**（feature screening 已內含層級化假設）

---

## 驗收標準

| 項目 | 標準 |
|------|------|
| C16 FDR | 90 組合 AUC + adjusted p 表；surviving combos 列表 |
| C22 FDR | 5 zones × multi-metric adjusted p 表（bug-fix 後） |
| C17 meta | Fisher's method meta p + per-sample bootstrap CI |
| O11/12/13 追溯 | FDR 前後 AUC 對比；NEGATIVE 穩固性確認 |
| C07 Bonferroni | 7 特徵 adjusted p + LOSO CI |
| C11 bootstrap | per-sample F1 95% CI + overlap 圖 |
| 文件化 | FDR 方法論寫入 `/known-pitfalls` skill |

---

## 預期 if-then 影響

- **If** C16 NG=4+AF<0.4+NR≥80 fail FDR **then** 降為 CONDITIONAL，保留 NG=4 單一 filter
- **If** Zone Z3 TP rate difference fail FDR **then** C22 characterization 結論降為「trend-level」
- **If** O11/O12/O13 pre-stratify AUC fail FDR **then** 強化「從未真 POSITIVE」敘事
- **If** C07 G1-G7 Bonferroni 後仍全 <0.64 **then** NO-GO 結論完全穩固
- **If** C11 F1 CI overlap 廣 **then** 論文敘事降為「方法示範」而非「實用 filter」

---

## 與其他 cross-cutting 的耦合

- **Pooled_OLS_Audit**：C16/C17 同時需 within-group + FDR，合併執行
- **CovM_Bug_Impact**：C17/C22 FDR 需在 CovM bug fix 後進行
- **Characterization_Functions**：FDR 結果直接影響 F2/F3 成熟度評級

---

## 整體評分

**🟡 Medium-risk systematic gap — C16 的 90 組合 screening 是最大 FDR 風險；Zone-Aware 與 Phase 1A 次之；O11-13 補算強化 NEGATIVE 結論。修正後所有 POSITIVE 聲明具備 reviewer-defensible 嚴謹度。**
