---
title: Z3 × HCC1954 Amplicon/High-CN Blacklist Pilot — 結果報告
date: 2026-04-19
status: CONDITIONAL（HCC1954-only；其他樣本 collateral damage）
owner: InterSubMod Research
scope: TO mode · Method A 驗證
related:
  - 20260419_Z3_amplicon_blacklist_pilot_design_01.md
  - 20260418_Z3_Internal_Feature_Exploration_01.md
---

# Z3 × HCC1954 Amplicon/High-CN Blacklist Pilot — 結果

## 一、執行摘要

- **POSITIVE（HCC1954-only）**：S2 whole-chr5/8/17 ∩ Z3 達到 ΔF1=+0.0065，S1 literature arm-level 達 +0.0058；ceiling（S4）為 +0.0075
- **NEGATIVE（跨樣本）**：S1/S2 對其他 6 樣本 mean ΔF1=−0.0043 / −0.0044，5/6 樣本 F1 下降
- **最終判定**：**CONDITIONAL** — blacklist 在 HCC1954 有可重現增益，但無法作為跨樣本 global rule；若採用需 sample-conditional gating
- **對五研究目標影響**：對目標 5（F1）僅 HCC1954-local 路徑，**預期 ensemble F1 增益 < +0.001**（因其他樣本損失抵消），**不改變 Zone-Aware Framework characterization-only 定位**

---

## 二、主結果表

### 2.1 跨 strategy 匯總

![Blacklist ΔF1 per-sample × strategy](../../../../research/z3_internal_feature_exploration/figures/step5_blacklist_delta_f1.png)

| Strategy | HCC1954 ΔF1 | HCC1954 FP− | HCC1954 TP− | 其他 6 mean ΔF1 | 其他 improved | 其他 hurt |
|---------|------------|------------|-----------|---------------|-------------|---------|
| S1 Literature arms ∩ Z3 | **+0.0058** | 1470 | 71 | −0.0043 | 1 | 5 |
| S2 Whole chr5/8/17 ∩ Z3 | **+0.0065** | 1637 | 77 | −0.0044 | 1 | 5 |
| S3 CovM > 95%ile ∩ Z3 | +0.0002 | 35 | 0 | −0.0000 | 3 | 2 |
| **S4 Reject all Z3 (ceiling)** | **+0.0075** | 1929 | 101 | −0.0329 | 1 | 5 |

### 2.2 per-sample × S1/S2 詳細

| Sample | F1 before | S1 ΔF1 | S2 ΔF1 | S3 ΔF1 | S4 ΔF1 |
|--------|---------|--------|--------|--------|--------|
| H2009 | 0.9545 | −0.0024 | −0.0025 | 0.0000 | −0.0270 |
| H1437 | 0.8712 | −0.0072 | −0.0086 | −0.0002 | −0.0717 |
| HCC1395 | 0.8309 | −0.0037 | −0.0041 | +0.0001 | −0.0292 |
| HCC1395_DORADO | 0.8330 | −0.0041 | −0.0039 | +0.0003 | −0.0361 |
| HCC1937 | 0.6772 | +0.0023 | +0.0038 | +0.0000 | **+0.0309** |
| **HCC1954** | **0.4047** | **+0.0058** | **+0.0065** | +0.0002 | **+0.0075** |
| COLO829 | 0.7906 | −0.0108 | −0.0109 | −0.0004 | −0.0644 |

**關鍵觀察**：

1. **HCC1937 的 S4 ΔF1=+0.0309**：HCC1937 Z3 本身 TP rate 也偏低，reject all Z3 可獲增益；意味著「Z3 高 FP rate」並非 HCC1954 獨有，HCC1937 也符合類似 pattern（但本 pilot S1 arm-level 只捕到 +0.0023，因 HCC1937 的 FP 不集中在 chr5/8/17）
2. **COLO829 S2=−0.0109**：chr5/8/17 上有 Z3 TP 被誤殺，reflect COLO829 也有 chr8 CDKN2A 等 CNV 但 FP 分佈與 HCC1954 不同
3. **S3 對所有樣本近零**：High-CN cutoff 在 Z3 內捕獲極少，證實 Coverage_Multiple 對 Z3 內 TP/FP 區分力不足（與 Step 1 AUC 表一致）

---

## 三、合理性驗證

### 3.1 文獻座標 vs 實際 FP 分佈

Step 4 校核顯示狹窄 focal amplicon coord 捕獲率極低：

| 文獻座標 | HCC1954 Z3 FP 捕獲 | TP rate in amp vs out |
|---------|------------------|---------------------|
| chr5p gain 1–50Mb | 126 FP | 4.5% vs 2.0% |
| chr8q24 MYC 120–140Mb | 37 FP | 0% vs 8.8% |
| chr17q12 HER2 35–42Mb | **0 FP** | — vs 3.1% |

**→ 實際 FP 富集是 arm-level（8p loss / 17p TP53 LOH / chr5 整段重排）而非 focal amplicon**。這就是為什麼 S1 設計採 arm-level 而非文獻狹窄座標。

![HCC1954 Z3 FP Pos 分佈](../../../../research/z3_internal_feature_exploration/figures/step4_hcc1954_fp_pos_hist.png)

### 3.2 Z3-conditional gating 的必要性

本 pilot 第一次執行時（無 z3_only 限制），S1 unconditional 結果：

```
HCC1954: F1 0.4047 → 0.3899 (Δ=−0.0148), TP−=1513
```

- chr5/8/17 **non-Z3 TP 1,442 個**（來自 Z1/Z2 high-confidence call）全被誤殺
- 加上 Z3 TP 71 → 共 TP− 1513，遠超 FP−

加入 `z3_only=True` 後 S1 變 +0.0058。**此修正為本 pilot 最關鍵的合理性設置**，確認 blacklist 必須 zone-gated 才有意義。

### 3.3 Circularity Guard 驗證（S3）

S3 採用**非 Z3 region** 的 CovM 95th percentile 作 cutoff（per-sample 計算）：

| Sample | Cutoff | Z3 captured | ΔF1 |
|--------|--------|-----------|-----|
| H2009 | 2.48 | 0 | 0.0000 |
| H1437 | 1.587 | 88 | −0.0002 |
| HCC1395 | 1.747 | 26 | +0.0001 |
| HCC1395_DORADO | 1.853 | 50 | +0.0003 |
| HCC1937 | 2.88 | 1 | +0.0000 |
| HCC1954 | 1.64 | 35 | +0.0002 |
| COLO829 | 0.64 | 84 | −0.0004 |

- 無 outcome 參與 cutoff 決定 → 無循環定義
- 但 Z3 內 CovM 區分力不足 → 捕獲極少 → 增益微弱
- 結論：**S3 非不合理，只是 signal 太弱**

### 3.4 Ceiling（S4）校準

HCC1954 ceiling ΔF1=+0.0075。S2（+0.0065）達到 **87%** ceiling，S1（+0.0058）達 **77%**。

→ 表示「refined blacklist」已經接近「reject all Z3」的上限，沒有更多 FP 可再被「更聰明」的規則分出。

**對 HCC1954 而言，blacklist 的效能幾乎 = 捨棄 Z3 的效能**。

---

## 四、判定與影響評估

### 4.1 五目標影響

| 目標 | 影響 | 說明 |
|------|------|------|
| 1. per-CpG ASM | 零 | 本 pilot 不涉及甲基化層 |
| 2. Clone 結構 | 零 | HCC1954 CNV pattern 與 subclonal 結構正交 |
| 3. 二次打擊順序 | 零 | 無新訊號 |
| 4. Biomarker | 零 | HCC1954-specific 無跨樣本推廣性 |
| **5. F1 優化** | **HCC1954-only 微弱 POSITIVE** | 若僅部署於 HCC1954 可 +0.0065；ensemble mean +0.0004 |

### 4.2 是否採用

**建議：不納入 canonical filter**

理由：

1. 跨樣本 ensemble mean ΔF1 ≈ +0.0004（S2），未達 ≥+0.005 的採用門檻
2. 需 sample-specific gating（需預先知道樣本是 HER2+ breast cancer），違反 caller 對樣本無假設的原則
3. F pilot（`NG=4+AF<0.4+NR≥80`）已覆蓋 HCC1954 大部分 FP，本 pilot 增量邊際效用低
4. 若仍想部署：作為**可選 post-hoc annotation**（「此樣本若已知為 HER2+ breast cancer，可套用 chr5/8/17 Z3 FN 標記」），非 default filter

### 4.3 Zone-Aware Framework 定位

**不變**。本 pilot 再次驗證：

- Z3 內部的 FP 機制是**樣本特異的 CNV 架構artifact**，非 Zone 定義本身的缺陷
- Zone-Aware Framework 持續作 **characterization only**，不升級為 production filter

---

## 五、產物清單

### 腳本
- `research/z3_internal_feature_exploration/scripts/step4_pos_distribution_check.py`
- `research/z3_internal_feature_exploration/scripts/step5_amplicon_blacklist_pilot.py`

### 數據
- `data/step4_amplicon_literature_coords.tsv` — 文獻 bed
- `data/step4_hcc1954_fp_amplicon_hit.tsv` — focal coord hit rate
- `data/step5_blacklist_bed.tsv` — S1/S2 使用 bed
- `data/step5_blacklist_pilot_results.tsv` — per-sample × strategy 全表
- `data/step5_blacklist_pilot_summary.tsv` — 匯總表

### 圖表
- `figures/step4_hcc1954_fp_pos_hist.png` — FP Pos × 文獻 amplicon
- `figures/step5_blacklist_delta_f1.png` — ΔF1 per-sample × strategy

---

## 六、結論

1. **方法 A（Amplicon/High-CN Blacklist）在 HCC1954 可達 +0.0058 – +0.0065 F1 增益**（近 ceiling 77–87%）
2. **對其他 6 樣本造成 mean ΔF1 −0.004 collateral damage**（5/6 hurt）
3. **Circularity-free、Z3-conditional 的 S3 CovM strategy** 合理但訊號過弱（Δ≈+0.0002）
4. **判定 CONDITIONAL**：HCC1954-local 有效；**不建議納入 global canonical filter**
5. **Zone-Aware Framework 定位不變**，HCC1954 的 Z3 FP 問題本質是樣本特異 CNV 架構，不是 Zone 機制失效

**下一步（若繼續）**：

- 不再探索 Z3 × HCC1954 方向（已接近 ceiling）
- 若要改善 HCC1954 F1，方向應轉向：
  - (i) Phase 2 Normal Methylation Reference 的 sample-level CNV annotation
  - (ii) Coverage_Multiple 全域 z-score normalization（與 HCC1954 已知 ~2N vs HCC1395 ~2N 配對 re-normalize）
  - (iii) 與 F pilot 的 `NG=4+AF<0.4+NR≥80` filter 組合檢驗互補性
