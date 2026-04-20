# Pooled OLS Residualization Audit（跨卡 pooled OLS trap 審查）

> **建立日期**: 2026-04-19

## 背景

**陷阱定義**（`feedback_pooled_ols_residualization_trap.md`）：
- Pooled OLS residualization 在「confound 變量與分組標籤相關」時**無法移除分組信號**
- 正確作法：**within-group OLS**（每組內分別 residualize），或提供分組作為 co-variate
- O12 AlleleDelta vs AF 是此陷阱的標誌案例（L2 collider bias 同時揭露）

**相關 feedback**：
- `feedback_L2_collider_bias.md`（residualize on AF 產生虛假信號）
- `feedback_pooled_ols_residualization_trap.md`（within-group 替代方案）
- `/known-pitfalls` skill

**用戶決定（R-05 建議）**：回溯所有 residualized 結論，一次修完避免後續再補

---

## 跨卡 pooled OLS 風險矩陣

| Card | 結論主題 | Pooled OLS 使用位置 | 風險等級 | 驗證狀態 |
|------|---------|-------------------|---------|---------|
| **C04** | O11 Heterogeneity | epipolymorphism residualize on n_reads | 🟢 已揭露 | NEGATIVE 歸因 n_reads confound |
| **C05** | O12 LOH Scenarios | AlleleDelta residualize on AF | 🔴 **陷阱標誌案例** | NEGATIVE 歸因 AF collider bias |
| **C06** | O13 Cross-region | cross-region metric residualize on shared reads | 🟢 已揭露 | NEGATIVE 歸因 shared read confound |
| **C15** | LOH Methylation failure | 多特徵 residualize | 🟡 待驗證 | 需 within-group 重算 |
| **C16** | HPFineNGroups marker | NGroups residualized AUC | 🟡 **需 P0-B 驗證** | AUC 單點聲明，無 CI |
| **C17** | LOH Subclone AF×Methylation | Inter AF→NGroups r=+0.705 | 🔴 **需 P0-B within-group OLS 驗證** | 可能保留 AF-bin 分組信號 |

**圖例**：🔴 需立即驗證 / 🟡 待驗證 / 🟢 已揭露（作為 NEGATIVE 歸因）

---

## 風險結論詳述

### 🔴 C17 LOH Subclone (Inter AF→NGroups r=+0.705)

**現狀**：
- 報告 r=+0.705（7/7 p<10⁻³⁹）
- 使用 pooled OLS 可能在 residualize 時保留 AF-bin 分組信號
- Confound 問題：NGroups 可能與 AF 分組相關而非獨立 signal

**必要驗證（P0-B）**：
1. 按 AF bins（0.05, 0.1, 0.2, 0.3, 0.4）分層
2. Within-group OLS residualize on NG
3. 重算 Inter AF vs residualized NGroups 相關係數
4. 驗收：若 within-group r 仍 >0.5 且 7/7 顯著 → 結論穩固
5. 驗收：若 within-group r 降至 <0.3 → 結論反轉，pooled OLS trap 確認

**預期結果**（基於 F pilot 2026-04-18 經驗）：
- Within-group r 可能保留 0.5-0.65 範圍
- 方向 robust（正相關）但 effect size 縮小

---

### 🟡 C16 HPFineNGroups（NG=4+AF<0.4+NR≥80）

**現狀**：
- Filter-based，residualized AUC 聲明跨樣本一致
- 缺 bootstrap CI（P1-C 待補）
- Pooled OLS 在 NG residualize on AF 的步驟中可能保留 NG↔AF 耦合

**必要驗證**：
1. 按 AF bins 分層，within-group 驗證 NG=4 的 TP rate 是否仍 >85%
2. 1000× bootstrap CI for residualized AUC
3. 驗收：CI lower bound >0.58 且 within-group 7/7 一致 → POSITIVE 穩固

---

### 🟡 C15 LOH Methylation Failure

**現狀**：
- LOH methylation 系列特徵全 NEGATIVE
- 但文件描述「residualized on X」步驟缺具體方法論說明
- 需確認 NEGATIVE 結論是否因 pooled OLS 把真實 signal 也 residualize 掉

**必要驗證**：
1. 列出 C15 使用的所有 residualization 步驟
2. 對「最接近 POSITIVE 的邊緣特徵」做 within-group 重算
3. 驗收：若 within-group 仍 NEGATIVE → 結論穩固
4. 驗收：若 within-group 出現 POSITIVE 邊緣特徵 → 需升級為 CONDITIONAL

---

## 已揭露結論（Green List，pooled OLS trap 已識別為 NEGATIVE 根因）

### C04 O11 Heterogeneity

- 陷阱已識別：n_reads confound
- AUC 全部因 n_reads 解釋
- NEGATIVE 結論穩固
- **不需重做**

### C05 O12 LOH Scenarios

- 陷阱已識別：AlleleDelta=AF 本身（L2 collider bias）
- 經典 pooled OLS trap 案例
- NEGATIVE 結論穩固
- **不需重做**

### C06 O13 Cross-region

- 陷阱已識別：shared read count confound
- AUC 分層後消失
- NEGATIVE 結論穩固
- **不需重做**

---

## Within-Group OLS 標準作業流程（SOP）

### Step 1: 定義分組變量

- 對 AF-dependent 結論：AF bins = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
- 對 NR-dependent 結論：NR bins = [80, 100, 120, 150, 200, ∞]
- 對 n_reads-dependent 結論：log(n_reads) quintiles

### Step 2: Within-group OLS

```python
# 虛擬碼
for group in af_bins:
    subset = data[data['af_bin'] == group]
    if len(subset) < min_n:  # 避免 underpowered
        continue
    residuals[group] = OLS(y=subset['target'], x=subset['confound']).resid
```

### Step 3: 重算 target 指標

- AUC / r / effect size 在 residuals 上重算
- 每組分別報告，再 cross-group 整合（Fisher's method 或 meta-analysis）

### Step 4: 驗收

- **穩固**：cross-group 效應方向一致 + meta p<0.01 + per-group effect 大小相近
- **反轉**：pooled effect 消失或方向翻轉
- **部分穩固**：方向一致但 effect size 縮小 >50% → CONDITIONAL 標籤

---

## 建議執行優先序

### P0（立即，與 CovM bug fix 同時進行）

1. **C17 within-group OLS 驗證**（R-01 CovM 修正後自然觸發）
   - 預期 1-2 天完成
   - 結果寫回 C17 audit card + 00_INDEX.md

### P1（CovM bug fix 完成後）

2. **C16 HPFineNGroups within-group + bootstrap CI**
   - 預期 2-3 天
   - 結果寫回 C16 audit card

### P2（Phase 2 A+D 準備前）

3. **C15 LOH Methylation within-group 驗證**
   - 預期 3-5 天
   - 若揭露邊緣 POSITIVE → 可能重啟 LOH methylation 方向

---

## 不執行的回溯範圍（止損）

以下結論**不做** within-group 回溯，原因：
- C01 / C02 / C13 / C21：純 VCF coordinate / set operation，無 OLS
- C03 / C07 / C08：ISM feature AUC，不涉及 residualization
- C09 / C10 / C11 / C14：無 residualize on confound 步驟
- C12：ASM 度量直接輸出，無 residualization
- C18 / C19 / C20 / C22：B.2 series 的統計設計不同（per-sample z-score）

---

## 驗收標準

| 項目 | 標準 |
|------|------|
| C17 重算 | Within-group OLS r 值 + meta p value 寫回卡片 |
| C16 重算 | Bootstrap 1000× CI + per-bin AUC 表寫回卡片 |
| C15 重算 | 每個 LOH methylation 邊緣特徵 within-group AUC |
| 結論反轉檢查 | 任何結論若 within-group 與 pooled 結論不一致 → 00_INDEX.md 狀態更新 |
| 文件化 | SOP 寫入 `/known-pitfalls` skill 作為 L2 驗證要求 |

---

## 與 CovM bug fix 的耦合

**C17 同時受兩個問題影響**：
1. Pooled OLS 分組信號保留（此文件）
2. CovM bug 在 step3 CN1 分層（CovM_Bug_Impact.md）

**建議合併執行**：
- CovM bug fix 完成後，重算 step3 CN1 分層時，**同時**使用 within-group OLS
- 一次驗證兩個問題：r=+0.705 是否穩固 vs 變動
- 預期 r 會同時受兩個修正影響，但方向（正相關）應 robust

---

## 整體評分

**🟡 Medium-risk systematic issue — 3 個結論（C15, C16, C17）需 within-group 驗證；C17 同時與 CovM bug 耦合，建議合併修正。O12 作為經典案例已成為警告樣本。修正後 Phase 2 論文 residualized 聲明可具備可辯護性。**
