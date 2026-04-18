<!--
建立時間: 2026-04-13 18:30
目標: 分析 Germline ASM × Tumor ASM × LOH 交叉關係，判定 Normal BAM 加入是否提供更清楚的 characterization
處理範圍: HCC1395 全基因體 31,659 regions (30,401 TP + 1,258 FP)
關聯檔案:
  - docs/experiments/validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md
  - include/core/RegionProcessor.hpp (新增 4 個 signed delta 欄位)
  - src/core/RegionProcessor.cpp (計算邏輯 + CSV 輸出)
-->

# Germline ASM × Tumor ASM × LOH 交叉分析

## 1. 研究動機

Phase A-D Dual-BAM 架構完成後，ISM 首次同時擁有以下維度：
- **Normal HP Signed Delta**: 正常組織的 HP1-HP2 甲基化方向差異（germline ASM）
- **Tumor HP Signed Delta**: 腫瘤的 HP1-HP2 甲基化方向差異
- **Sample ASM Delta**: 腫瘤 vs 正常整體甲基化差異
- **LOH 標注**: 外部 LOH BED 區域重疊

核心問題：加入 Normal BAM 後，能否更清楚地分類位點？哪些位點的 ASM 是 germline 就存在的？腫瘤是否改變了 ASM 方向？LOH 如何影響 ASM？

## 2. 方法

### 2.1 C++ 新增 Signed HP Delta（4 個欄位）

在 `RegionProcessor` 新增 signed HP delta 計算：

```
Signed Delta = mean(HP1 methylation across all CpGs) - mean(HP2 methylation across all CpGs)
```

**與距離式 HP Delta 的區別**：
- 距離式 delta = `d_between - d_within`（PERMANOVA）：只衡量兩 HP 是否不同，**無方向**
- Signed delta = `mean(HP1_meth) - mean(HP2_meth)`：衡量 HP1 甲基化是否高於 HP2，**有方向**

4 個新欄位：`Tumor_HP_Signed_Delta`, `Normal_HP_Signed_Delta`, `HP_Signed_Residual`, `Combined_HP_Signed_Delta`

### 2.2 驗證

- C++ chr19 輸出 vs Python 全基因體獨立計算：max diff < 0.00002（浮點精度）
- 全基因體 signed delta 使用 Python 從既有 per-region methylation files 計算（16 workers）

### 2.3 數據

| 指標 | 值 |
|:---|:---|
| 總 regions | 31,659 |
| TP / FP | 30,401 / 1,258 |
| LOH regions | 1,713 (5.4%) |
| Tumor HP signed valid | 20,866 (65.9%) |
| Normal HP signed valid | 30,635 (96.8%) |
| Both valid | 20,517 (64.8%) |

## 3. 結果

### 3.1 HP Valid 率 × LOH（穩定性 4/5）

| 子集 | N | Tumor HP valid | Normal HP valid | Both valid |
|:---|:---:|:---:|:---:|:---:|
| 全部 | 31,659 | 65.9% | 96.8% | 64.8% |
| **LOH** | 1,713 | **23.8%** | 72.6% | 18.2% |
| Non-LOH | 29,946 | 68.3% | 98.2% | 67.5% |

**結論**：LOH 使 tumor HP valid 率從 68% 降至 24%。原因：LOH 物理消除一個 allele，tumor reads 嚴重偏向一個 HP（HP1/(HP1+HP2) = 0.592），導致另一 HP reads 不足。

### 3.2 Germline ASM 盤點（穩定性 4/5）

| 門檻 |delta| | 有 ASM 位點 | 比例 |
|:---:|:---:|:---:|
| > 0.01 | 27,311 | 89.1% |
| > 0.02 | 24,190 | 79.0% |
| > 0.05 | 16,028 | 52.3% |
| > 0.10 | 7,371 | 24.1% |
| > 0.15 | 3,260 | **10.6%** |
| > 0.20 | 1,419 | 4.6% |

- 方向平衡：HP1 > HP2 = 39.7%, HP2 > HP1 = 39.2%（完美對稱）
- **文獻比較**：Rosenski 2025 報告 11% CpG 為 bimodal（ASM）。ISM 在 |delta| > 0.15 得到 10.6%，數字接近但度量不同（ISM 是窗口平均甲基化差異，非單 CpG bimodal 判定）

### 3.3 方向一致性：50% = 隨機（穩定性 5/5，確定性結論）

| 比較 | Sign Agreement |
|:---|:---:|
| Tumor vs Normal | **49.6%** |
| Tumor vs Combined | 71.0% |
| Normal vs Combined | 60.8% |

**即使限制在高 delta 位點（|delta| > 0.2 in both），方向一致性仍為 51.0%。**

**根因分析**：
1. HP labels 在 ISM 中是一致的（同一 `hp_merged_labels` 向量）— T-Combined 71% 和 N-Combined 61% 確認此點
2. 非 HP label 問題，而是**腫瘤確實改變了 ASM 的方向**
3. 即使 tumor HP 很平衡（0.4-0.6），方向一致性仍為 50.3%
4. Pearson r(tumor signed, normal signed) = 0.002 — 零相關

**生物學解釋**：腫瘤表觀遺傳重編程影響不同的 CpGs。在 ±2kb 窗口內，germline ASM 的 CpGs 和 tumor 新增 ASM 的 CpGs 是不同位點，導致窗口平均方向獨立。與 Do et al. 2020 報告的 allele-specific hypomethylation 一致。

**對 ISM 的意義**：跨 tumor-normal 的 signed delta 方向比較無意義。HP_Signed_Residual（tumor - normal）反映的是隨機方向差異而非 somatic ASM 改變。

### 3.4 Tumor/Normal ASM 比例（穩定性 3/5）

| 子集 | |Tumor|/|Normal| | 文獻預測 |
|:---|:---:|:---:|
| 全部 | **1.12×** | 5-9× (Do 2020) |
| LOH | 1.70× | - |
| Non-LOH | 1.11× | - |

**與文獻不一致的原因**：ISM 看 somatic SNV 周圍 ±2kb 窗口（非隨機 CpG），且使用連續甲基化值（0-1）的窗口平均，不是二元判定。Do 2020 的 5-9× 是全基因體 CpG level 的比較。

### 3.5 LOH 對 ASM 的影響（穩定性 3/5）

LOH 區域 both valid 僅 311 regions（18.2%），分析受限：
- LOH 區域 |Tumor signed delta| = 0.1265 > |Normal signed delta| = 0.0744
- 預期相反（LOH 應消除 ASM），但 both-valid 樣本是有偏的（LOH 較弱才能兩 HP 都有足夠 reads）
- |Tumor| < |Normal| 的比例只有 39.5%

**結論**：LOH 主要的效果是讓 tumor HP 測試失效（valid 率 24%），而非可觀察的 ASM 消除。要觀察 EVOFLUx 預測的三態→兩態崩塌，需要 per-CpG level 分析（非窗口平均）。

### 3.6 位點分類系統（穩定性 4/5）

門檻 |signed delta| > 0.05：

| 類別 | N | 比例 | FP rate | FP enrichment | 特徵 |
|:---|:---:|:---:|:---:|:---:|:---|
| **No_ASM** | 6,209 | 19.6% | **2.05%** | **0.51×** | 最乾淨；|tumor|=0.020, |normal|=0.023 |
| **Tumor_Homogenized** | 5,808 | 18.3% | 2.69% | 0.68× | Normal 有 ASM 但 tumor 消失 |
| Normal_ASM_Tumor_Unknown | 4,836 | 15.3% | 4.20% | 1.06× | Normal 有 ASM, tumor HP 不足 |
| **Shared_ASM** | 4,718 | 14.9% | **6.53%** | **1.64×** | 兩者都有 ASM，FP 最集中 |
| Low_ASM_Tumor_Unknown | 4,350 | 13.7% | 3.40% | 0.86× | 兩者低 ASM, tumor HP 不足 |
| **Tumor_Acquired** | 3,471 | 11.0% | 6.08% | 1.53× | Tumor 新獲得 ASM |
| LOH_Eliminated_ASM | 666 | 2.1% | 4.20% | 1.06× | LOH + normal 有 ASM |
| LOH_No_Preexisting_ASM | 577 | 1.8% | 5.03% | 1.26× | LOH + normal 無 ASM |

**核心發現**：
- **No_ASM 類別 FP 率最低（2.05%，0.51× enrichment）** — 無 ASM = 最可靠的位點
- **Shared_ASM 類別 FP 率最高（6.53%，1.64× enrichment）** — 兩者都有 ASM = 最不可靠
- FP enrichment 差異（1.64× vs 0.51×）有統計意義但效果量太小，不足以作為 filter

### 3.7 TP vs FP 信號探索（穩定性 4/5，預期 negative 確認）

| Feature | Raw AUC | |abs| AUC |
|:---|:---:|:---:|
| Normal_HP_Signed_Delta | 0.489 | 0.470 |
| **Tumor_HP_Signed_Delta** | **0.600** | 0.332 |
| **HP_Signed_Residual** | **0.614** | 0.371 |
| Combined_HP_Signed_Delta | 0.605 | 0.352 |
| HP_Residual_Delta (unsigned) | 0.436 | 0.477 |
| SampleASM_Delta | 0.430 | 0.429 |
| HPFineF | 0.610 | 0.612 |

- **HP_Signed_Residual raw AUC = 0.614**：最高但低於 HPFineF（0.612）
- Signed delta raw AUC > unsigned delta AUC：方向信息有增量
- 但全部 AUC < 0.62，與先前結論一致：ISM 特徵無法有效區分 TP/FP

**FP tumor ASM 來源確認**：FP 的 tumor HP balance = 0.381（偏斜），TP = 0.500（平衡）。FP 低 AF → reads 偏向一個 HP → 表面更高 ASM。這是 AF confound，非新信號。

## 4. 文獻對照總表

| ISM 觀察 | 文獻預測 | 一致？ | 說明 |
|:---|:---|:---:|:---|
| Germline ASM |delta|>0.15: 10.6% | 11% bimodal (Rosenski 2025) | ≈ | 度量不同但數量級一致 |
| Tumor/Normal ASM = 1.12× | 5-9× (Do 2020) | ✗ | ISM 窗口平均 ≠ 全基因體 per-CpG |
| 方向一致性 50% | 預期同方向 if maintained | ✗→ | 腫瘤 ASM 影響不同 CpG，符合 Do 2020 |
| FP 更高 tumor ASM | AF confound | ✓ | 低 AF → HP imbalance → 表面 ASM |
| LOH tumor HP valid 24% | LOH 消除一個 allele | ✓ | EVOFLUx 三態崩塌需 per-CpG 分析 |

## 5. 結論

### 5.1 加入 Normal BAM 是否使結果更清楚？

**是的，提供了 characterization 增量，但未產生新的 TP/FP 區分信號。**

1. **新增維度有效**：6 類位點分類系統能區分不同表觀遺傳狀態
2. **Germline ASM 基線**：10.6-52% 位點有 germline ASM（依門檻），提供正常參考
3. **LOH 影響量化**：LOH 使 tumor HP 分析失效率從 32% 升至 76%
4. **方向獨立性**：腫瘤 ASM 方向與 germline 完全獨立，這是重要的生物學發現

### 5.2 TP/FP 區分

- 所有新 signed delta 特徵 AUC < 0.62
- FP enrichment 最大 1.64×（Shared_ASM）→ 不足以做 filter（需 > 3× 才有實用價值）
- **與先前所有研究結論一致：ISM 無法有效過濾 FP**

### 5.3 生物學發現

1. **方向 50%**：腫瘤不維持 germline ASM 方向，表觀重編程影響不同 CpGs
2. **No_ASM 最乾淨**：無 ASM 位點 FP rate 最低（0.51×），ASM 本身是噪聲環境指標
3. **Tumor_Homogenized**：Normal 有 ASM 但 Tumor 消失（均質化），FP rate 低（0.68×）— tumor 消除 ASM 反而是好信號

### 5.4 穩定性等級

| 結論 | 穩定性 | 理由 |
|:---|:---:|:---|
| 方向一致性 = 50% | 5/5 | 20,517 regions, 各門檻一致 |
| LOH 降低 tumor HP valid | 4/5 | 明確的物理機制 |
| Germline ASM ~11% at |delta|>0.15 | 4/5 | 與文獻一致 |
| 6 類分類系統有效 | 4/5 | 各類特徵分離清晰 |
| Signed delta 無 filter 價值 | 4/5 | 所有 AUC < 0.62 |
| Tumor/Normal ratio 1.12× | 3/5 | 窗口大小和度量方式影響大 |
| LOH ASM 消除效果 | 3/5 | Both-valid 樣本偏差，需更多數據 |

## 6. 後續建議

1. **7 samples 驗證**：用其他 6 個癌症樣本驗證方向 50% 結論
2. **Per-CpG 分析**：如需驗證 EVOFLUx 三態崩塌，需 per-CpG level（非窗口平均）
3. **Phase 2A baseline**：此分析確認 Normal BAM 的 characterization 價值，支持 Phase 2A 正式分析

## 7. 數據位置

| 數據 | 路徑 |
|:---|:---|
| 全基因體 CSV + signed delta | `/tmp/ism_phase_b_test/output_full_phase_cd/significance_summary_with_signed_delta.csv` |
| chr19 C++ 驗證 | `/tmp/ism_phase_b_test/output_chr19_signed_delta/significance_summary.csv` |
| TP/FP labels | `/tmp/ism_phase_b_test/tp_fp_labels.tsv` |
| LOH BED | `/tmp/ism_phase_b_test/ism_loh_regions.bed` |
| 計算腳本 | `/tmp/ism_phase_b_test/compute_signed_delta.py` |

## 8. C++ 修改清單

| 檔案 | 修改 |
|:---|:---|
| `include/core/RegionProcessor.hpp` | 新增 4 欄位 + constructor 初始化 |
| `src/core/RegionProcessor.cpp` | 計算邏輯（lambda `compute_signed_hp_delta`）+ CSV header + data output |

173/173 tests passed. chr19 C++ vs Python 交叉驗證 max diff < 0.00002。
