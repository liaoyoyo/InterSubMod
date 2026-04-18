<!--
建立時間: 2026-04-17 22:00
目標: Zone-Aware Confidence Framework H1/H3 假設驗證
處理範圍: all_region_rows.tsv.gz 748,391 rows × 7 samples × 2 modes
關聯檔案:
  - docs/concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md
  - scripts/analysis/validate_zone_hypotheses_h1_h3.py
-->

# H1/H3 Zone-Aware Hypothesis Validation Report

> **驗證日期**: 2026-04-17
> **腳本**: `scripts/analysis/validate_zone_hypotheses_h1_h3.py`
> **數據**: `all_region_rows.tsv.gz` (748,391 rows, post-HP-fix)
> **圖表目錄**: `research/zone_aware_validation/figures/`

---

## 總覽

| 假設 | Paired 結論 | TO 結論 | 整體判定 |
|------|-----------|--------|---------|
| **H1**: Z1 TP rate > Global | 方向正確 6/7，但 N 太少（1,106），統計不顯著 | 方向正確 5/7，4/7 significant，但 N 極少（153） | **CONDITIONAL** — 方向確認但覆蓋率不足 |
| **H3**: Z2 TP rate ≥ 89.1% in TO | 7/7 confirmed（mean=0.988） | **3/7 only**（mean=0.716）| **PARTIAL** — Paired 確認，TO 不成立 |

---

## H1: Zone Z1（LOH Subclonal Active）TP Rate Validation

### Z1 定義

```
LOH = True AND caller_af ∈ [0.1, 0.4] ∪ [0.6, 0.9] AND HPFineNGroups ≥ 3
```

LOH 來源：Paired 使用 `core_loh_like`（ISM HP ratio），TO 使用 `to_loh_bed_hit`（LOH.bed，不受 self-phasing 影響）。

### Paired Mode 結果

| Sample | Global TP Rate | Z1 TP Rate | Delta | N(Z1) | p-value | Sig |
|--------|:-:|:-:|:-:|:-:|:-:|:-:|
| H2009 | 0.999 | 1.000 | +0.001 | 355 | >0.05 | ns |
| H1437 | 1.000 | 1.000 | +0.000 | 54 | >0.05 | ns |
| HCC1395 | 0.979 | **0.805** | **-0.174** | 411 | >0.05 | ns |
| HCC1395_DORADO | 0.992 | 0.992 | +0.000 | 132 | >0.05 | ns |
| HCC1937 | 0.985 | 1.000 | +0.015 | 52 | >0.05 | ns |
| HCC1954 | 0.998 | 1.000 | +0.002 | 100 | >0.05 | ns |
| COLO829 | 0.940 | 1.000 | +0.060 | 2 | >0.05 | ns |

**Paired 小結**：6/7 正方向，0/7 顯著。Paired baseline TP rate 太高（0.94-1.00），Z1 無法再提升。HCC1395 是唯一反例（-0.174），原因待查。

### TO Mode 結果

| Sample | Global TP Rate | Z1 TP Rate | Delta | N(Z1) | p-value | Sig |
|--------|:-:|:-:|:-:|:-:|:-:|:-:|
| H2009 | 0.913 | 0.900 | -0.013 | 10 | >0.05 | ns |
| H1437 | 0.772 | 1.000 | +0.228 | 4 | >0.05 | ns |
| HCC1395 | 0.711 | **0.989** | **+0.279** | 93 | <0.001 | *** |
| HCC1395_DORADO | 0.714 | **1.000** | **+0.286** | 20 | <0.01 | ** |
| HCC1937 | 0.512 | **0.714** | **+0.202** | 21 | <0.05 | * |
| HCC1954 | 0.254 | **0.800** | **+0.546** | 5 | <0.05 | * |
| COLO829 | 0.654 | — | — | 0 | — | — |

**TO 小結**：5/7 正方向，4/7 顯著。Z1 在 TO 模式中確實 TP-enriched，delta 範圍 +0.20 ~ +0.55。但覆蓋率極低（僅 153 regions 佔全 TO 的 0.04%）。

### H1 展示圖

![H1 Validation](figures/H1_zone_z1_tp_rate_validation.png)

### H1 關鍵觀察

1. **Z1 定義太嚴格**：要求同時滿足 LOH + Intermediate AF + NGroups≥3，導致覆蓋率極低（Paired 0.3%, TO 0.04%）
2. **TO 的 to_loh_bed_hit 覆蓋率低**：TO 模式 LOH.bed 標記的 region 本身就少（因 self-phasing 不影響 LOH.bed，但 LOH.bed 來自 VCF phasing，TO 模式 VCF 品質不同）
3. **方向一致但實用價值存疑**：即使 Z1 TP rate 高，153 個 region 的 rescue 對 F1 影響微乎其微
4. **HCC1395 Paired 異常**：Z1 TP rate 僅 0.805，是唯一低於 global 的樣本。可能與 HCC1395 chr8 LOH+ASM 熱點有關

### H1 判定

**CONDITIONAL**：方向確認（Z1 確實 TP-enriched），但 Z1 定義需放寬以提升覆蓋率。建議：
- 放寬 NGroups 門檻到 ≥2（而非 ≥3）
- 或：拆開 Z1 為 LOH+Intermediate_AF（不限 NGroups）+ NGroups≥3（不限 LOH），分別評估

---

## H3: Zone Z2（High Somatic Heterogeneity）TP Rate Validation

### Z2 定義

```
HPFineNGroups ≥ 4 AND NumReads ≥ 80
```

### Paired Mode 結果

| Sample | Global TP Rate | Z2 TP Rate | Delta | N(Z2) | p-value | Sig |
|--------|:-:|:-:|:-:|:-:|:-:|:-:|
| H2009 | 0.999 | **1.000** | +0.001 | 7,582 | <0.05 | * |
| H1437 | 1.000 | 1.000 | +0.000 | 1,315 | >0.05 | ns |
| HCC1395 | 0.979 | **0.989** | +0.010 | 1,129 | <0.01 | ** |
| HCC1395_DORADO | 0.992 | 0.993 | +0.001 | 546 | >0.05 | ns |
| HCC1937 | 0.985 | **1.000** | +0.015 | 593 | <0.001 | *** |
| HCC1954 | 0.998 | 1.000 | +0.002 | 621 | >0.05 | ns |
| COLO829 | 0.940 | 0.933 | -0.006 | 15 | >0.05 | ns |

**Paired 小結**：**7/7 above 89.1%**（mean=0.988）。89.1% claim 在 Paired 模式完全確認。

### TO Mode 結果

| Sample | Global TP Rate | Z2 TP Rate | Delta | N(Z2) | p-value | Sig |
|--------|:-:|:-:|:-:|:-:|:-:|:-:|
| H2009 | 0.913 | **0.934** | +0.022 | 19,983 | <0.001 | *** |
| H1437 | 0.772 | **0.921** | +0.149 | 1,410 | <0.001 | *** |
| HCC1395 | 0.711 | 0.810 | +0.099 | 1,173 | <0.001 | *** |
| HCC1395_DORADO | 0.714 | **0.903** | +0.189 | 639 | <0.001 | *** |
| HCC1937 | 0.512 | 0.713 | +0.201 | 891 | <0.001 | *** |
| HCC1954 | 0.254 | 0.497 | +0.243 | 1,622 | <0.001 | *** |
| COLO829 | 0.654 | 0.235 | -0.418 | 34 | >0.05 | ns |

**TO 小結**：**3/7 ≥ 89.1%**，mean=0.716。89.1% claim 在 TO 模式**不成立**。但 6/7 顯著優於 global（delta +0.02 到 +0.24）。

### H3 展示圖

![H3 Validation](figures/H3_zone_z2_tp_rate_validation.png)

### H3 關鍵觀察

1. **Paired 完全確認**：89.1% 的 TP rate 在 Paired 模式中 7/7 樣本均成立（最低 0.933 for COLO829，但 N=15）
2. **TO 絕對 TP rate 不達標**：TO 模式 Z2 mean TP rate = 0.716，遠低於 89.1%。原因：TO FP 基線高（TO 整體 TP rate 僅 0.647），NGroups≥4 的「somatic heterogeneity」判斷在 TO 中被 self-phasing artifact 干擾
3. **TO 相對提升一致**：6/7 樣本 Z2 顯著優於 global，表示 NGroups≥4 仍是有效的 TP 正向指標，只是無法達到 89.1% 的絕對門檻
4. **COLO829 TO 反向**：Z2 TP rate = 0.235 < Global 0.654。COLO829 TO 模式 FP 比例極高，Z2 N=34 過少
5. **HCC1954 TO 值得注意**：Z2 TP rate = 0.497，高於 global 的 0.254（2× 提升），但仍遠低於 89.1%

### H3 判定

**PARTIAL**：
- Paired 模式：89.1% claim **確認**（7/7）
- TO 模式：89.1% 絕對值 **不成立**（3/7），但**相對提升效應確認**（6/7 significant）
- **修正宣稱**：「NGroups≥4 + NR≥80 is a significant positive TP indicator across both modes, with Paired TP rate ≥ 93% and TO delta +2~24pp over global」

---

## Full Zone Assignment Statistics

### Zone 分布

| Zone | Paired N (%) | Paired TP Rate | TO N (%) | TO TP Rate |
|------|:-:|:-:|:-:|:-:|
| **Z1** (LOH Subclonal Active) | 1,106 (0.3%) | 0.927 | 153 (0.04%) | 0.941 |
| **Z2** (High Somatic Hetero) | 11,788 (3.6%) | **0.999** | 25,745 (6.1%) | **0.891** |
| **Z3** (Complete LOH) | 48,125 (14.6%) | 0.987 | 52,594 (12.5%) | **0.608** |
| **Z4** (Normal Diploid) | 147,174 (44.8%) | 0.996 | 187,311 (44.6%) | 0.694 |
| **Z5** (CN Gain Low Diversity) | 11,898 (3.6%) | 0.998 | 18,213 (4.3%) | **0.667** |
| Unassigned | 108,608 (33.0%) | 0.981 | 135,676 (32.3%) | 0.695 |

### Zone Overview

![Zone Overview](figures/zone_overview_statistics.png)

### Zone 排序（按 TP Rate）

**Paired**（所有 zone TP rate > 0.92，差異小）:
```
Z2 (0.999) > Z5 (0.998) > Z4 (0.996) > Z3 (0.987) > Unassigned (0.981) > Z1 (0.927)
```

**TO**（差異顯著，可操作）:
```
Z1 (0.941) > Z2 (0.891) > Z4 (0.694) ≈ Unassigned (0.695) > Z5 (0.667) > Z3 (0.608)
```

### 關鍵發現

1. **TO 模式 Zone 差異遠大於 Paired**：TO Z2 (0.891) vs Z3 (0.608) 差距 0.283，Zone-Aware 在 TO 模式有實質意義
2. **Z3（Complete LOH）在 TO 是最低 TP rate zone**：0.608。這是 self-phasing 在 extreme AF LOH 區域大量產生 artifact HP imbalance 的結果
3. **Z2 在 TO 的 0.891 與 Paired 的 89.1% claim 巧合吻合**：TO Z2 mean 恰好是 0.8912 ≈ 89.1%（全域聚合值，但 per-sample 異質性大）
4. **Z5（CN Gain + Low Diversity）在 TO 確實偏低**：0.667，支持 Zone-Aware Framework 中的「略微加嚴」策略
5. **Unassigned 佔 33%**：需要進一步細化 zone 定義以減少未分配比例

---

## Z1 vs Z3 對比（Subclonal Active vs Complete LOH）

| 模式 | Z1 Mean TP Rate | Z3 Mean TP Rate | Mean Delta (Z1-Z3) |
|------|:-:|:-:|:-:|
| Paired | 0.971 | 0.974 | -0.003 |
| TO | 0.901 | 0.518 | **+0.383** |

**TO 模式 Z1 vs Z3 差距極大（+0.38）**，且 per-sample 方向全部一致（6/6 有數據的樣本 Z1 > Z3）。這強烈支持 Zone-Aware 在 TO 模式的核心假設：**Subclonal LOH（Z1）遠比 Complete LOH（Z3）更可信**。

---

## 對 Zone-Aware Framework 的修正建議

### Z1 定義需放寬

- 當前 Z1 覆蓋率（Paired 0.3%, TO 0.04%）過低，無實用價值
- 建議選項：
  - **A**: 降 NGroups 門檻到 ≥2（仍保留 LOH + Intermediate AF 的核心邏輯）
  - **B**: 分離為 Z1a（LOH + Intermediate AF）和 Z1b（任意 + NGroups≥3），分別評估
  - **C**: 直接以 LOH + Intermediate AF 為 Z1（移除 NGroups 條件），NGroups 獨立作為 Z2

### H3 宣稱需修正

- 原始：「NGroups≥4 + NR≥80 → 89.1% TP rate」
- 修正：「NGroups≥4 + NR≥80 → **Paired**: TP rate ≥ 93% (7/7); **TO**: 顯著高於 global (+2~24pp, 6/7 p<0.001) 但絕對值 ~72%」

### TO 模式 Zone-Aware 是真正的價值場景

- Paired 所有 zone TP rate > 0.92，zone 差異太小無法產生 F1 改進
- TO 模式 zone TP rate 範圍 0.61-0.94，**差異足以支持差異化策略**
- Zone-Aware QS 調整應聚焦 TO 模式

---

## 輸出檔案

| 檔案 | 說明 |
|------|------|
| `H1_zone_z1_tp_rate.tsv` | H1 驗證數據（per sample × mode） |
| `H3_zone_z2_tp_rate.tsv` | H3 驗證數據（per sample × mode） |
| `zone_assignment_statistics.tsv` | 全 zone 統計（per mode） |
| `zone_per_sample_statistics.tsv` | Zone × Sample × Mode 交叉統計 |
| `Z1_vs_Z3_contrast.tsv` | Z1 vs Z3 對比 |
| `figures/H1_zone_z1_tp_rate_validation.png` | H1 驗證圖 |
| `figures/H3_zone_z2_tp_rate_validation.png` | H3 驗證圖 |
| `figures/zone_overview_statistics.png` | Zone 全覽統計圖 |
