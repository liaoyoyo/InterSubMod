# CovM Bug Impact Matrix（跨卡 CovM baseline bug 影響清單）

> **建立日期**: 2026-04-19

## 背景

**Bug 識別**：`expected_coverage=75.0` hardcoded default fallback，master dataset 全 7 樣本 × 2 mode（Paired/TO）共用同一值，KDE-based per-sample baseline 未啟用。

**來源**：
- `docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md`
- `research/coverage_multiple_validation/`
- C-INFRA-1（Phase 2 Role 3 審查）

**影響機制**：
1. CovM = observed_coverage / expected_coverage
2. expected_coverage 若 hardcoded 75.0 → CovM 成為 observed_coverage × (1/75.0) 線性轉換
3. 跨樣本比較時樣本間 coverage 差異被抹平 → zone 分層失真
4. Per-sample 排序（單樣本內）相對順序保留 → 樣本內結論相對 robust

**用戶決定（R-01）**：立即修 + 重跑（本審查完成後優先啟動 `/cpp-change`）

---

## 影響矩陣

| Card | 結論主題 | Bug 影響等級 | 影響機制 | 修正後預期 |
|------|---------|-------------|---------|-----------|
| **C17** | LOH Subclone AF (雙模式) | 🔴 **直接** | r=+0.705 Inter AF→NGroups 的 step3 CN1 分層依賴 CovM | r 值可能保留但分層 cutpoint 需重定 |
| **C18** | HCC1954 CNV-driven reversal | 🟡 **部分** | Z3 chr5/8/17 CovM=0.733 per-sample 特徵 | 絕對值變，相對排序 robust |
| **C19** | Z3 amplicon blacklist | 🟡 **邊界** | Z3 zone 定義依賴 CovM 分層 | 邊界微變，結論不變 |
| **C20** | Coverage_Multiple 非獨立 CN proxy | 🔴 **Bug 直接相關** | z_extreme 0.15% 的絕對值與 r=0.831 均受 bug | 需重算 r 值與 z_extreme |
| **C21** | LOH.bed 不受 self-phasing 汙染 | 🟢 **無影響** | 純 VCF coordinate，不讀 CovM | 無需重跑 |
| **C22** | Zone-Aware F1 + Characterization | 🔴 **整條 framework** | Z1-Z5 zone 邊界依賴 CovM 分層 | Zone 需重定義；H2 NEGATIVE 穩定 |
| **B.2-2** | HPFineNGroups CovM 關聯 | 🟡 **間接** | CovM 作為 NGroups 解釋變量 | 關聯仍在，強度可能變 |
| **Q3** | LOH_CN_AF verification | 🟡 **間接** | CovM 作為 CN proxy 驗證 | 需重驗證 r 值 |

**圖例**：🔴 必須重算核心數值 / 🟡 絕對值變但相對結論穩 / 🟢 無影響

---

## 受汙染結論的「預期修正後值」追蹤表

| 結論 | 原始數值 | Bug 性質 | 修正後預期範圍 | 結論是否反轉風險 |
|------|---------|---------|--------------|---------------|
| C20 CovM vs CN r | 0.831 | Hardcoded 扁平化跨樣本差異 | **0.5-0.9 範圍**（不確定，per-sample normalize 可能使 r 上升或下降） | 中度（若 r<0.5 則 CN proxy 地位下降） |
| C20 z_extreme | 0.15% | 跨樣本 baseline 錯誤 | **可能 0.1%-1%** | 低（極端值特徵穩） |
| C17 Inter AF→NGroups r | +0.705 | step3 CN1 分層依賴 | **0.65-0.75**（CN1 分層 cutpoint 變） | 低（7/7 顯著不倚賴絕對 r） |
| C18 HCC1954 CovM=0.733 | 0.733 | Per-sample normalization 差異 | **絕對值變，相對 Z3 特徵穩** | 低 |
| C19 Z3 blacklist ΔF1 | HCC1954 +0.0065 / 其他 -0.0044 | Z3 邊界微變 | **邊界 ±5% 變動；ΔF1 不變** | 極低 |
| C22 Zone TP rate | Z3 0.608, Z1 0.85 | Zone 分層變 | **Zone 重定義；TP rate 差異仍在** | 中（Z3 可能併入 Z2） |

---

## 修正後需重跑的 pipeline

### 必跑（R-01 後立即）

1. **CovM KDE baseline 重算**（per-sample）：7 樣本 × 2 mode × 2 region type
2. **C20 CovM vs CN r 重算**：cross_sample_audit 資料重處理
3. **C22 Zone 重定義**：Z1-Z5 邊界以修正後 CovM 重分層

### 應跑（R-01 後兩週內）

4. **C17 step3 CN1 分層重算**：Inter AF→NGroups 關係重驗證
5. **C19 Z3 blacklist pilot 重跑**：確認 ΔF1 結論穩定
6. **C18 HCC1954 z-score 重算**：per-sample 特徵穩定性驗證

### 可延後（Phase 2 A+D 開始前）

7. **B.2-2 HPFineNGroups CovM 關聯**：解釋變量強度重估
8. **Q3 LOH_CN_AF verification**：cross-dimension 整合更新

---

## 不受影響的結論（Green List）

以下結論**完全不依賴 CovM**，可直接繼承現有結論：

- **C01 Paired/TO FP rate 分離**（純 VCF coordinate）
- **C02 PON coverage**（baseline 度量，非 CovM 衍生）
- **C03 TO AUC ceiling**（ISM 特徵，與 coverage 無關）
- **C04 O11 Heterogeneity**（confound 源是 n_reads 非 CovM）
- **C05 O12 LOH Scenarios**（confound 源是 AF 非 CovM）
- **C06 O13 Cross-region**（confound 源是 shared reads 非 CovM）
- **C07 G1-G7 Germline**（methylation 特徵，非 CovM）
- **C08 Read-level FP**（read-level 特徵）
- **C09 Self-Phasing causal chain**（VCF coordinate）
- **C10 PON-only phasing**（VCF coordinate）
- **C11 Phase 1A F1 optimization**（ML 特徵已 stratify）
- **C12 ASM**（methylation 特徵）
- **C13 LOH.bed SEQC2**（VCF coordinate）
- **C14 QS TO failure**（observational，已歸因到 C07）
- **C15 LOH Methylation failure**（methylation 特徵）
- **C16 HPFineNGroups marker**（NGroups + NR filter，無 CovM）
- **C21 LOH.bed clean**（VCF coordinate）

**計數**：17/22 結論無影響 / 5/22 結論受影響（C17, C18, C19, C20, C22）/ C16 間接（via B.2-2）

---

## R-01 修正優先序

### Phase A（立即，本審查後）

1. `/cpp-change` 觸發：修 `src/core/*.cpp` KDE-based expected_coverage 啟用路徑
2. 編譯 + 單元測試
3. 7 樣本 × 2 mode 重跑 ISM baseline（Master dataset v6 候選）

### Phase B（重跑後一週內）

4. C20 重算 + update audit card
5. C22 Zone 重定義 + update audit card
6. C17 step3 重算 + update audit card

### Phase C（整合）

7. C18, C19, B.2-2 post-fix 確認
8. 00_INDEX.md CovM bug 狀態更新：🔴 active → 🟢 resolved

---

## 驗收標準

| 項目 | 標準 |
|------|------|
| CovM KDE 啟用 | `cmake`/`make` 編譯通過，單樣本 KDE baseline 數值產出 |
| Per-sample baseline 差異化 | 7 樣本 expected_coverage 值**非全部 75.0** |
| 跨樣本排序穩定性 | 修正前後 per-sample z-score rank 相關 >0.8 |
| C22 zone 重定義 | Z1-Z5 邊界文件化；Z3 cardinality 變化 <30% |
| 結論反轉檢查 | C17/C20 r 值若變動 >0.2 → 升級為「待重新論證」 |

---

## 風險與 if-then 分析

- **If** C20 r<0.5 **then** 「CovM ≈ CN proxy」聲明進一步 FALSIFIED → 需改以 CNV caller 輸出取代
- **If** C22 Z3 zone cardinality <30% 原值 **then** Z3 blacklist（C19）結論需重新驗證
- **If** C17 r<0.5 **then** LOH Subclone 雙重證據鏈（AF + Methylation）需降為單重
- **If** 修正後 per-sample baseline 全部落在 70-80 範圍 **then** bug 實質影響小，結論大致保留
- **If** 修正後樣本間差異顯著（例 HCC1395 65, HCC1954 85）**then** 所有 cross-sample CovM 結論需重評

---

## 整體評分

**🔴 Active critical bug — 本審查後最優先修正 + 重跑，影響 5-6 個結論的跨樣本比較，per-sample 相對排序 robust。修正後 Phase 2 A+D 才有乾淨 baseline 起點。**
