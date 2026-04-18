<!--
建立時間: 2026-04-05 01:00
目標: 使用正確的 ClairS-TO VCF 重跑 ISM 後，全面重新驗證所有受影響結論
處理範圍: 所有基於 ism_pononly_v2b_tp/fp 的 TO 模式分析
關聯檔案:
  - docs/experiments/in_progress/2026/04/20260404_VCF來源錯誤矯正報告_01.md
  - docs/experiments/in_progress/2026/04/20260404_HPFineP_QS整合完整研究報告_01.md
  - docs/experiments/in_progress/2026/04/20260404_LOH_Strong_Weak_7feature_AUC驗證報告_01.md
-->

# ClairS-TO VCF 矯正全面重分析報告

## 1. 背景

### 1.1 問題回顧

ISM TO 分析使用了 **ClairS paired pileup VCF**（`##source=ClairS`），而非正確的 **ClairS-TO VCF**（`##source=ClairS-TO`）。已用正確 VCF 重跑 ISM 並取得新結果。

### 1.2 數據集

| 數據集 | 路徑 | TP | FP | VCF Source |
|--------|------|----|----|------------|
| **新（正確）** | `ism_pononly_v2b_clairsto_{tp,fp}` | 28,383 | 11,830 | `##source=ClairS-TO` |
| **舊（錯誤）** | `ism_pononly_v2b_{tp,fp}` | 30,476 | 4,822 | `##source=ClairS` (paired) |

### 1.3 VCF 驗證確認

- `grep "##source=" ClairS-TO VCF` → `ClairS-TO`（正確）
- 無 NAF/NDP/NAD 欄位（TO 模式無 Normal BAM）
- TP: 28,396 input → 28,383 analyzed (99.95%)
- FP: 11,843 input → 11,830 analyzed (99.89%)

---

## 2. 核心結論影響

### 2.1 總覽表

| # | 結論 | 舊值 | 新值 | 狀態 | 說明 |
|---|------|------|------|------|------|
| 1 | **QS 區分 TP/FP** | AUC=0.579 | **AUC=0.507** | **🔴 FAILED** | QS 等同隨機，完全無區分力 |
| 2 | **7-Feature LR** | AUC=0.760 | **AUC=0.627** | **🟡 WEAKENED** | 仍 >0.5 但大幅下降 |
| 3 | **NGroups≤1 指標** | FP/TP=2.54 | **FP/TP=1.66** | **🟡 WEAKENED** | 效力減半 |
| 4 | **LOH+HPFineSig** | ratio=3.00 | **ratio=1.49** | **🟡 WEAKENED** | 大幅弱化 |
| 5 | **Scheme 3 改善** | ΔAUC=+0.026 | **ΔAUC=+0.028** | **🔴 MEANINGLESS** | 基線 0.507，改善無實際意義 |
| 6 | **LOH_S AUC** | 0.895 | **0.721** | **🟡 WEAKENED** | 仍 >0.5 但下降顯著 |
| 7 | **LOH_W AUC** | 0.863 | **0.739** | **🟡 WEAKENED** | 仍 >0.5 但下降顯著 |
| 8 | **F1 改善有限** | — | — | **🟢 CONFIRMED** | 更加確認 ISM filter 效果有限 |
| 9 | **TO/Paired 分離** | 架構決策 | 架構決策 | **🟢 CONFIRMED** | FP 性質完全不同，更強化 |
| 10 | **AlleleDelta 區分** | FP>>TP | **FP≈TP** | **🔴 REVERSED** | 真正 TO FP 無 AlleleDelta 差異 |

### 2.2 詳細分析

#### 2.2.1 QS 完全失效（🔴 最嚴重）

```
舊（Paired FP）: QS AUC = 0.579
  Tier High:  TP 64.3%, FP 47.4%
  Tier Low:   TP  8.1%, FP 17.1%

新（ClairS-TO）: QS AUC = 0.507（≈ 隨機）
  Tier High:  TP 91.2%, FP 90.4%  ← 幾乎完全重疊！
  Tier Low:   TP  0.0%, FP  0.0%  ← 沒有任何 Low tier
```

**根因**：ClairS-TO germline FP 具有真正的生物學信號（HP groups、methylation patterns、高 AF），這些都是 QS 用來給高分的特徵。QS 無法區分「真正 somatic 的生物信號」與「germline 的生物信號」。

**影響**：QS 框架的 TO 模式需要根本重新設計。

#### 2.2.2 7-Feature AUC 下降（🟡 中度影響）

```
舊: Global=0.760, LOH_S=0.895, LOH_W=0.863
新: Global=0.627, LOH_S=0.721, LOH_W=0.739
```

- AUC 仍然 >0.5，表示 ISM 特徵組合有一定區分力
- 但遠低於之前估計的 0.76-0.90
- **LR 係數變化**：
  - 舊: AlleleDelta (-0.864) 是最強特徵
  - 新: HPFineP (-0.346) 取代 AlleleDelta 成為最強；AlleleDelta 降至 -0.149
  - 這證實 AlleleDelta 的舊效力主要來自 paired FP 的特性

#### 2.2.3 AlleleDelta 無差異（🔴 反轉）

```
舊: TP mean=0.034, FP mean=0.082 (FP >> TP)
新: TP mean=0.034, FP mean=0.034 (FP ≈ TP)
```

**根因**：
- Paired FP（mapping artifacts）缺乏真正的 allele 結構 → AlleleDelta 大（noise）
- ClairS-TO FP（germline variants）有正常的 allele structure → AlleleDelta 與 TP 相同
- 之前的 FP>>TP 差異完全是 paired mode artifacts 造成的假象

#### 2.2.4 NGroups 弱化（🟡）

```
舊: NGroups≤1 FP/TP ratio=2.54
新: NGroups≤1 FP/TP ratio=1.66
```

Germline FP 有真正的 HP groups（因為它們是真正的 heterozygous variants），所以 NGroups 分布更接近 TP。但仍然 FP>TP，可能因為部分 germline FP 是 homozygous（NGroups=0-1）。

#### 2.2.5 HPFineSig Ratio 弱化（🟡）

```
LOH: 舊 ratio=3.00, 新 ratio=1.49
All: 舊 ratio=1.75, 新 ratio=1.22 (estimated)
```

Germline FP 在 LOH 區域有真正的 allele-specific methylation（因為 germline variants 位於不同 haplotypes），所以 HPFineSig 觸發率上升。

---

## 3. Quality Score 分布細節

### 3.1 TP vs FP Quality Score

```
NEW (ClairS-TO):
  TP: mean=89.52, median=100.00
  FP: mean=89.06, median=100.00
```

TP 和 FP 的 QS 分布幾乎完全相同。這意味著 QS 的所有組件（VerificationClass、LOH penalty、coverage 等）對 germline FP 的行為與 TP 相同。

### 3.2 VerificationClass

```
NEW (ClairS-TO):
  Strong  : TP=23.8%, FP=19.5%
  Weak    : TP=34.5%, FP=28.8%
  Subclone: TP= 5.3%, FP= 7.9%
  Noise   : TP=36.4%, FP=43.9%
```

差異存在但非常小。Noise FP 稍多（+7.5pp），但不足以驅動有效的 QS 區分。

### 3.3 PassedGating

```
NEW: TP=32.0%, FP=30.6% (差異僅 1.4pp)
```

Gating 對 TO FP 幾乎無效。

---

## 4. Scheme 3 評估

### 4.1 組件效力

| 組件 | 舊 | 新 | 評估 |
|------|----|----|------|
| NGroups≤1 Gap | +0.276 | **+0.121** | 弱化但仍正向 |
| LOH+HPFineSig Ratio | 3.00 | **1.49** | 弱化但仍 >1 |
| Scheme 3 ΔAUC | +0.026 | **+0.028** | 絕對改善相當，但基線 0.507 |

### 4.2 Scheme 3 參數掃描

| Penalty/Bonus | 舊 ΔAUC | 新 ΔAUC |
|---------------|---------|---------|
| ±3 | +0.023 | +0.026 |
| ±5 | +0.026 | +0.028 |
| ±10 | +0.040 | +0.036 |

Scheme 3 在新數據上 ΔAUC 略有增加，但：
- 基線 QS AUC = 0.507（隨機）
- Scheme 3 後 QS AUC = 0.535
- 這個 AUC 仍然接近隨機，**實際上無法用於任何過濾決策**

### 4.3 結論

**Scheme 3 在 TO 模式下不可行**。不是因為 ΔAUC 小，而是因為整個 QS 框架在 TO 模式下就是隨機的。

---

## 5. 仍然有效的信號

### 5.1 7-Feature LR（AUC=0.627）

雖然大幅下降，但仍然 >0.5。主要有效特徵：
- **HPFineP** (-0.346)：HP fine p-value，仍有一定區分力
- **HP_Ratio** (+0.287)：LOH 模式中有區分力
- **HPFineNGroups** (+0.152)
- **AlleleDelta** (-0.149)：殘留的微弱效力

### 5.2 LOH 子類型 AUC

- LOH_Strong: AUC=0.721
- LOH_Weak: AUC=0.739

LOH 區域內仍有一定區分力，但遠不如之前估計。

### 5.3 VerificationClass Noise 比例差異

Noise 在 FP 中佔 43.9% vs TP 36.4%（差 7.5pp），這是唯一仍然明顯的信號。

---

## 6. 對研究策略的影響

### 6.1 必須暫停的工作

| 項目 | 原因 |
|------|------|
| **Task #21: Scheme 3 實作** | QS 框架 TO 模式完全無效，Scheme 3 建立在無效基線上 |
| **所有基於 QS 的 TO 過濾研究** | QS AUC=0.507 |

### 6.2 需要更新的結論

| 結論 | 更新 |
|------|------|
| AlleleDelta FP>>TP | ❌ 錯誤 — 基於 paired FP artifacts |
| NGroups≤1 是強 FP 指標 | ⚠️ 弱化 — ratio 從 2.54 降至 1.66 |
| LOH+HPFineSig ratio=3.00 | ⚠️ 弱化 — 降至 1.49 |
| 7-Feature AUC > 0.76 | ⚠️ 下調至 0.627 |
| LOH_S/W AUC > 0.88 | ⚠️ 下調至 0.72-0.74 |

### 6.3 不受影響的結論

| 結論 | 原因 |
|------|------|
| TO/Paired 必須分離 | 更加強化（FP 性質完全不同） |
| ISM filter F1 改善有限 | 確認（現在更確定） |
| Self-phasing circular dependency | 程式碼邏輯，與 VCF 無關 |
| PON-only phasing 驗證 | 與 VCF 來源無關 |
| 所有 Paired mode 分析 | 使用正確 VCF |

### 6.4 新的關鍵洞察

1. **QS 在 TO 模式下需要根本不同的設計**：不能依賴 VerificationClass/HPFineSig/LOH 等「生物信號指標」，因為 germline FP 擁有同樣的信號。
2. **AlleleDelta 在 TO 模式無用**：germline variants 有真正的 allele structure。
3. **ISM 的 TO 價值可能在於組合特徵而非單一指標**：7-feature AUC=0.627 仍然 >0.5，暗示某種非線性組合可能有用。
4. **方向應轉向 FN rescue 而非 FP filtering**：覆蓋率和 TP loss 仍然是主要瓶頸。

---

## 7. 圖表

| 圖 | 檔名 | 內容 |
|----|------|------|
| Fig 1 | `fig1_overview_comparison.png` | 6 子圖：QS/NGroups/HPFineSig/VerificationClass/AlleleDelta/LOH 分布比較 |
| Fig 2 | `fig2_auc_scheme3_comparison.png` | AUC 條形圖 + Scheme 3 ΔAUC 參數掃描 |
| Fig 3 | `fig3_quality_tier_stacked.png` | QS Tier 堆疊圖（最關鍵：顯示 90% FP 獲得 High tier） |

---

## 8. ISM 輸出路徑

| 輸出 | 路徑 |
|------|------|
| ClairS-TO TP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_clairsto_tp/` |
| ClairS-TO FP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_clairsto_fp/` |
| 分析圖表 | `/big7_disk/liaoyoyo2001/InterSubMod/output/synthesis/concluded/clairsto_correction_analysis/` |
| 舊（Paired）TP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_tp/` |
| 舊（Paired）FP ISM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/ism_pononly_v2b_fp/` |

---

## 9. 時間線

| 時間 | 事件 |
|------|------|
| 2026-04-04 17:30 | 發現 VCF 來源錯誤 |
| 2026-04-04 18:30 | VCF 來源錯誤矯正報告完成 |
| 2026-04-04 23:45 | 用 ClairS-TO VCF 啟動 ISM 重跑（TP+FP 雙線並行） |
| 2026-04-05 00:30 | FP ISM 完成（11,830/11,843 成功） |
| 2026-04-05 00:50 | TP ISM 完成（28,383/28,396 成功） |
| **2026-04-05 01:00** | **本報告：全面重分析完成** |
