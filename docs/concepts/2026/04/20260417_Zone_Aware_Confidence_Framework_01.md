<!--
建立時間: 2026-04-17 21:00
目標: 將 LOH/CN/AF 觀察轉化為可行動的 Zone-Aware Confidence Framework，連結五大研究目標
處理範圍: 全部 LOH/CN/AF 驗證結果 × 五大目標的交叉分析
關聯檔案:
  - docs/architecture/20260327_InterSubMod研究願景定錨_01.md
  - docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md
  - research/loh_cn_af_verification/20260417_LOH_CN_AF_結論驗證報告_01.md
  - docs/experiments/validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md
-->

# Zone-Aware Confidence Framework

> **狀態**: 構想階段（待 H1/H3 驗證）
> **動機**: 將 14 個月 LOH/CN/AF 研究觀察從「已關閉的 filter 方向」轉化為「可行動的分區差異化策略」
> **核心轉換**: Zone Exclusion (failed) → Zone-Aware Confidence Adjustment

---

## 一、背景：為什麼「排除」失敗但「分區」可能有效

### 排除策略失敗的根本原因

所有 zone 排除策略（排除 LOH、排除 CN_Loss、排除 CN_Gain 等）均失敗，因為：

1. **TP >> FP**：paired 模式 TP:FP = 24:1，任何排除都大量損失 TP
2. **一刀切**：排除整個 zone 無差別移除所有 regions，不考慮 zone 內的異質性
3. **絕對數量不利**：最佳策略（排除 CN_Loss）仍然每移除 1 FP 損失 12 TP

### 分區策略為什麼不同

分區策略不移除 regions，而是**調整信心度**：

- LOH 區域 TP-enriched → 不排除，而是「這個區域的 somatic call 更可信」
- NGroups >= 4 的 89.1% TP rate → 不排除其他，而是「高 NGroups 是強 TP 指標」
- 不改變 precision-recall trade-off，而是**提供額外 evidence layer**

---

## 二、五大目標 × 觀察映射

### 映射矩陣

| 觀察結果 | 目標 1 (per-CpG) | 目標 2 (Clone) | 目標 3 (二次打擊) | 目標 4 (TO補強) | 目標 5 (F1) |
|---------|:-:|:-:|:-:|:-:|:-:|
| LOH ALL tiers TP-enriched (0.766-0.903) | - | - | - | - | **核心** |
| Subclone AF x NGroups +0.642 (7/7) | 間接 | **核心** | **核心** | - | 間接 |
| HPFineNGroups N>=4 TP=89.1% | - | **核心** | 間接 | - | **核心** |
| CovM CN proxy r=0.997 | - | 間接 | 間接 | **核心** | 間接 |
| Methylation CN-blind (|r|<0.07) | - | - | - | 排除項 | 排除項 |
| NGroups x LOH inverse (NGroups=1 LOH=76.7%) | **核心** | **核心** | **核心** | - | 間接 |
| Phase D 4-group subclone | **核心** | **核心** | 間接 | **核心** | 間接 |
| Sample ASM 97.3% significant | **核心** | 間接 | 間接 | **核心** | - |

### 各目標的具體收益

**目標 1 (per-CpG)**：Zone 提供 CpG association scoring 的解讀 context
- LOH 區域：HP imbalance 是**預期的**（一個 allele 消失），不是 error
- Normal Diploid：HP imbalance 是**真 ASM 信號**，可信度高
- 影響：per-CpG 分數的可信度需標記 zone，不同 zone 的同一分數有不同含義

**目標 2 (Clone)**：AF gradient 是最直接的 subclone 證據
- Intermediate AF LOH → NGroups 高 → subclonal LOH（部分細胞保留雙 allele）
- Extreme AF LOH → NGroups 低 → complete LOH（所有細胞失去一個 allele）
- Phase D 的 4-group 分類可以進一步在 LOH group 內用 AF gradient 細分

**目標 3 (二次打擊)**：三軸組合推論事件順序
- LOH + Extreme AF + low NGroups + high Sample ASM = 古老 LOH → 甲基化已完成重塑
- LOH + Intermediate AF + high NGroups + high Sample ASM = 近期 LOH → ASM 保留
- Non-LOH + high NGroups + high HP ASM = allele-specific methylation（無 LOH 的表觀遺傳事件）

**目標 4 (TO)**：Coverage_Multiple 和 LOH.bed 是 TO 模式的可靠資訊源
- CovM r=0.997 as CN proxy → 不需外部 CNV truth set
- LOH.bed (VCF AF/VAF) 不受 self-phasing 影響 → TO 模式可用
- Normal methylation reference 提供期望基線

**目標 5 (F1)**：Zone-Aware QS 調整的核心應用場景
- 見下方第三節詳細設計

---

## 三、Zone-Aware Confidence Framework 設計

### 3.1 Zone 定義

```
Zone Z1: LOH Subclonal Active [2026-04-17 已放寬]
  原始條件: LOH_Bed_Overlap=true AND caller_af in [0.1, 0.4] u [0.6, 0.9] AND HPFineNGroups >= 3
    → 覆蓋率極低 (TO 0.04%, N=153)，統計不穩定

  修正條件: LOH=true AND caller_af in [0.1, 0.4] u [0.6, 0.9] AND HPFineNGroups >= 2
    → 覆蓋率提升 100× (TO 4.6%, N=19,266)
    → TO TP Rate: 0.965 (7/7 positive, 6/7 significant, delta +0.27)
    → Paired TP Rate: 0.982

  生物推論: 部分 LOH 保留 ASM 多樣性，≥2 groups = subclonal structure 存在
  策略: 信心度大幅上調

Zone Z2: High Somatic Heterogeneity [2026-04-17 H3 驗證完成]
  條件: HPFineNGroups >= 4 AND NumReads >= 80（不論 LOH）
  Paired TP 率: 98.8%（7/7 ≥ 89.1%，確認）
  TO TP 率: 71.6%（3/7 ≥ 89.1%，絕對值不達標，但 6/7 顯著高於 global）
  生物推論: 4+ methylation groups = active somatic subclone diversity
  策略: 信心度上調（TO 模式效果較弱但方向一致）

Zone Z3: Complete LOH [2026-04-17 驗證完成]
  條件: LOH_Bed_Overlap=true AND (caller_af < 0.1 OR caller_af > 0.9) AND HPFineNGroups <= 1
  Paired TP 率: 0.987（高，與預期一致）
  TO TP 率: **0.608**（所有 zone 中最低！self-phasing artifact 主因）
  生物推論: 完全 LOH，一個 allele 消失，甲基化趨均質化
  策略: Paired 維持標準；**TO 應降低信心度**（Z3 在 TO 是 FP 高風險 zone）

Zone Z4: Normal Diploid
  條件: Subclone_ID = 0（Phase D Group 0）
  TP 預期: 標準（HP ratio=0.503，最低 ASM）
  生物推論: 正常二倍體區域，somatic SNV 未顯著影響甲基化
  策略: 標準判定

Zone Z5: CN Gain with Low Diversity
  條件: Coverage_Multiple > 1.5 AND HPFineNGroups <= 2 AND NOT LOH
  TP 預期: 略偏低
  生物推論: CN gain 增加 read noise，低 NGroups 暗示非 somatic diversity
  策略: 略微加嚴（需更強 evidence）
```

### 3.2 Two-Pass Confirmation 機制

```
Pass 1: Primary QS Scoring (現有系統)
  Input: methyl+context features (Phase 1A 已鎖定)
  Output: QS_primary, TP/FP classification

Pass 2: Zone-Aware Adjustment
  Input: QS_primary + Zone annotation (Phase A-D 輸出)
  Logic:
    QS_final = QS_primary + zone_adjustment

    zone_adjustment (TO mode — Paired 差異太小暫不調整):
      Z1 (LOH Subclonal Active):    +delta_1 (TP rate 0.965, delta +0.27 vs global)
      Z2 (High Somatic Hetero):     +delta_2 (TP rate 0.891, delta +0.20 vs global)
      Z3 (Complete LOH):            -delta_3 (TP rate 0.608, 最低 zone！)
      Z4 (Normal Diploid):           0        (TP rate 0.694, 接近 global baseline)
      Z5 (CN Gain Low Diversity):   -delta_5 (TP rate 0.667, 略低於 global)
```

### 3.3 Rescue 候選機制

Two-Pass 的另一個用途：**邊界 case 的 rescue**

```
if QS_primary in [borderline_low, borderline_high]:
    if zone in [Z1, Z2]:
        decision = "RESCUE as TP"  (zone evidence 支持 TP)
    elif zone in [Z3, Z5]:  # Z3 加入！TO Z3 TP rate 僅 0.608
        decision = "CONFIRM as FP" (zone evidence 不利)
    else:
        decision = "KEEP original"
```

### 3.4 與現有 Phase 2 代碼的整合路徑

Phase A-D 代碼已完成，以下 annotation 欄位已存在於 ISM 輸出中：

| 已有欄位 | Zone 定義對應 | 來源 |
|---------|-------------|------|
| LOH_Bed_Overlap | Z1/Z3 LOH 判定 | LohBedAnnotator.hpp |
| LOH_Source | LOH 可靠度分級 | Phase C |
| Subclone_ID (0-3) | Z4 判定 | SubcloneAnalyzer.hpp |
| HPFineNGroups | Z2 判定 | ISM core features |
| Coverage_Multiple | CN proxy, Z5 判定 | ISM core features |
| caller_af | AF gradient, Z1 vs Z3 | VCF 輸入 |
| NumReads | Z2 門檻 | ISM core features |

**不需新增 C++ 代碼**。Zone-Aware 調整可在 QualityScore 計算階段（`QualityScorer.cpp`）或後處理 Python 腳本中實現。

---

## 四、假設驗證結果（2026-04-17 完成）

> 驗證報告：`research/zone_aware_validation/20260417_H1_H3_Zone_Validation_Report_01.md`
> 腳本：`scripts/analysis/validate_zone_hypotheses_h1_h3.py`

### H1: Zone Z1 TP rate > Global — CONDITIONAL

| 模式 | 正方向 | 顯著 | Mean Z1 TP Rate | Mean Delta | Total Z1 N |
|------|:-:|:-:|:-:|:-:|:-:|
| Paired | 6/7 | 0/7 | 0.971 | -0.014 | 1,106 |
| TO | 5/7 | 4/7 | 0.901 | +0.255 | 153 |

**判定**：方向確認（Z1 在 TO 模式顯著 TP-enriched），但 **Z1 覆蓋率過低**（Paired 0.3%, TO 0.04%）。Z1 定義需放寬（降低 NGroups 門檻或拆開條件）。

### H3: Zone Z2 89.1% TP rate in TO — PARTIAL

| 模式 | ≥89.1% | 顯著 | Mean Z2 TP Rate | Total Z2 N |
|------|:-:|:-:|:-:|:-:|
| Paired | **7/7** | 3/7 | 0.988 | 11,801 |
| TO | 3/7 | **6/7** | 0.716 | 25,752 |

**判定**：Paired 確認（7/7 ≥ 89.1%）。TO **絕對值不成立**（mean 0.716），但**相對提升效應確認**（6/7 significant, delta +0.02 ~ +0.24）。

**修正宣稱**：「NGroups≥4 + NR≥80 → Paired TP rate ≥ 93%; TO 顯著高於 global 但絕對值 ~72%」

### 意外發現：TO Zone 排序

```
TO: Z1 (0.941) > Z2 (0.891) > Z4 (0.694) ≈ Unassigned (0.695) > Z5 (0.667) > Z3 (0.608)
```

- **Z3（Complete LOH）在 TO 是最低 TP rate zone（0.608）**：self-phasing artifact 在 extreme AF LOH 區域大量產生虛假 HP imbalance
- **Z1 vs Z3 在 TO 差距 +0.383**：Subclonal LOH 遠比 Complete LOH 可信
- **Zone-Aware 的主要價值場景是 TO 模式**（Paired 所有 zone > 0.92，差異太小）

### 待驗證（需額外條件）

| # | 假設 | 依賴 | 狀態 |
|---|------|------|------|
| H2 | Zone-aware QS 調整能提升 TO F1 | 現有數據模擬 | ❌ **NEGATIVE**（max delta +0.001） |
| H4b | Z5 是否有更高 FP rate（SEQC2 CNV） | SEQC2 數據 | 初步確認（TO Z5=0.667） |
| H5 | Zone-aware rescue 能在低純度樣本改善 recall | 低純度樣本數據 | 待執行 |

**H2 失敗根因**：TO QS AUC=0.497（隨機），在隨機分數上做線性 zone 調整結果仍然隨機。Zone TP rate 差異真實但無法透過 QS adjustment 傳遞。[模擬報告](../../research/zone_aware_validation/20260417_QS_Simulation_Report_01.md)

---

## 五、影響評估（2026-04-17 模擬完成後更新）

### TO mode F1 — ❌ NEGATIVE（模擬確認）

- 5 種 zone delta 配置 × 21 QS 閾值 × 7 樣本 = 735 組合全面測試
- **Max delta F1 = +0.001**（HCC1954 only），6/7 樣本改善為 0
- **根因**：TO QS AUC=0.497，zone 調整無法修正隨機分數
- **QS 閾值過濾在 TO 模式全面有害**：6/7 樣本最佳 F1 在 T=0（不過濾）
- [模擬報告](../../research/zone_aware_validation/20260417_QS_Simulation_Report_01.md)

### Paired mode（未測試，預期小）

- Paired QS 有區分力（AUC > 0.5），zone 調整可能有正向效果
- 但 TP:FP=24:1，所有 zone TP rate > 0.92，delta 空間 < 0.005
- **可測試但優先序低**

### Characterization（目標 1-3）— ✅ Zone-Aware 的核心價值

- Zone annotation 提供 **region-level biological context**
- **不改 F1 也有獨立價值**：每個 somatic call 的 zone label 為生物學解讀提供框架
- 論文定位：「ISM 不只 call TP/FP，更提供 epigenetic zone context」
- Zone TP rate 差異可作為 confidence annotation（非 QS 調整，而是獨立標註）

### 重新定位

| 原始預期 | 驗證結果 | 調整方向 |
|---------|---------|---------|
| Zone-Aware QS 調整改善 TO F1 | ❌ QS 隨機，無法傳遞 | 放棄 QS 調整路線 |
| Zone 用於 post-filter rescue | ❌ 過濾不如不過濾 | 放棄 post-filter 路線 |
| Zone 作為 characterization annotation | ✅ Zone TP rate 差異真實 | **確認核心價值** |
| Zone 作為未來 ML 分類特徵 | 待測試 | 需重建 QS 模型才有意義 |

---

## 六、實施順序

```
Step 1: H1/H3 驗證 ✅（2026-04-17 完成）
  → H1 CONDITIONAL：Z1 方向確認但覆蓋率不足，需放寬定義
  → H3 PARTIAL：Paired 確認，TO 絕對值不成立但相對提升確認
  → 報告：research/zone_aware_validation/20260417_H1_H3_Zone_Validation_Report_01.md

Step 2: Zone annotation 統計 ✅（2026-04-17 完成）
  → 全 zone Z1-Z5 assignment + TP/FP/size 分佈已計算
  → TO 模式 zone TP rate 範圍 0.61-0.94（差異足以支持策略）
  → 數據：research/zone_aware_validation/zone_*.tsv

Step 3: Z1 定義放寬 ✅（2026-04-17 完成）
  → 7 個變體 Pareto 分析完成
  → Z1b (LOH+IntAF+NG≥2) 為最佳甜蜜點：覆蓋率 100× (4.6%), TP rate 0.965
  → Z3 在 TO 是最低 zone (0.608)，需加入負向調整
  → 腳本：scripts/analysis/zone_z1_relaxation_analysis.py
  → 數據：research/zone_aware_validation/Z1_relaxation_variants.tsv

Step 4: QS 調整模擬 ❌（2026-04-17 完成，NEGATIVE）
  → 5 configs × 21 thresholds × 7 samples = 735 組合
  → Max delta F1 = +0.001（HCC1954 only）
  → 根因：TO QS AUC=0.497 隨機，zone delta 無法修正
  → 報告：research/zone_aware_validation/20260417_QS_Simulation_Report_01.md

Step 5: C++ 整合 — 暫停
  → QS 調整路線已證實無效，C++ 整合無意義
  → 如果未來重建 TO QS 模型（用 zone 作為 primary feature），再重啟

Step 6: Zone Characterization Annotation（待規劃）
  → 在 ISM 輸出中加入 zone label 欄位
  → 用途：論文 evidence layer、生物學解讀框架
  → 不改變 F1，但增加輸出的解讀維度
```

---

## 七、與已關閉方向的區別

| 已關閉方向 | Zone-Aware 為什麼不同 |
|-----------|---------------------|
| LOH binary filter | 不排除 LOH regions，而是在 LOH regions 內調整信心度 |
| CNV zone exclusion | 不排除 CNV zones，而是用 CovM+NGroups 組合判斷 |
| TO germline FP post-hoc | 不依賴甲基化區分 germline FP，而是用 zone context 調整門檻 |
| Option C HP-free features | 不試圖找 HP-free discriminator，而是利用 HP-dependent features 在可信 zone 中發揮作用 |

**核心差異**：已關閉的方向都是「用甲基化特徵直接判別 TP/FP」（AUC < 0.64），Zone-Aware 是「用 region-level context 調整判別門檻」——不要求甲基化有 discriminative power，只要求 zone 有 TP rate 差異。
