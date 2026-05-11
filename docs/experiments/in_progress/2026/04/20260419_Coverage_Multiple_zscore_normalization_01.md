---
title: R6 — Coverage_Multiple 全域 z-score normalization
date: 2026-04-19
status: CONDITIONAL NEGATIVE（對 HCC1954 Z3 high-CovM 訊號不敏感，但揭示 Per-sample baseline 結構）
owner: InterSubMod Research
scope: Batch 2 軌 3 · CovM 是否為精確 CN proxy 的獨立訊號源
related:
  - 20260419_LOH_bed_generation_audit_01.md
  - 20260419_Z3_amplicon_blacklist_pilot_result_01.md
  - 20260418_B2_LOH_Subclone_Deep_Skeptical_Check_01.md
r_id: R6
priority: P1
---

# R6 — Coverage_Multiple 全域 z-score normalization

## 一、執行摘要

- **結論**：Per-sample non-LOH baseline 的 z-score normalization **吸收了**原 Z3 pilot Step 3 觀察到的「HCC1954 Z3 FP Coverage_Multiple=0.73 偏高」訊號。Raw CovM ≠ 精確 CN proxy，主要反映 **per-sample overall coverage scaling + 局部 CN 複合訊號**
- **HCC1954 Z3 FP z-extreme rate（z>2）**：TO 0.15%（3/1929），paired 0%（0/3）— 遠低於 chr5/8/17 FP 富集率（TO 84.4%）
- **真正的 Z3 FP 富集訊號是染色體位置分佈**（big3=chr5/8/17 collectively 80%），非 Coverage_Multiple 極端值
- **對 Z3 pilot 的影響**：強化既有結論「S3 CovM cutoff 對 Z3 FP 訊號微弱（ΔF1≈+0.0002）」。不改變 S1/S2 染色體 arm-level blacklist 的 HCC1954-only CONDITIONAL 定位
- **對 Opus 4.7 plan 假設 4（Coverage_Multiple ≈ CN 代理）**：❌ 否證 — Coverage_Multiple 作為 CN 代理在 per-sample normalization 後失去獨立訊號
- **對後續 R1 HCC1395 Normal BAM pilot 的影響**：**強化必要性** — 因 CovM 單獨訊號弱，Normal BAM annotation 的 Δ_ASM 成為唯一可能的獨立 CN-aware 訊號

---

## 二、審查問題定義

Opus 4.7 plan E.9 R6 原文：

> Coverage_Multiple 受 GC、mappability、purity 影響；若與 NGroups 共同受「local complexity」驅動，產生虛假 LOH×methylation 關聯

以及 Z3 pilot §4.2.4 下一步 (ii)：

> Coverage_Multiple 全域 z-score normalization（與 HCC1954 已知 ~2N vs HCC1395 ~2N 配對 re-normalize）

**驗證標準**：normalize 後 HCC1954 Z3 chr5/8/17 CovM 偏高仍顯著（z>2）→ 原訊號真實；消失 → 原訊號為 per-sample scaling artifact

---

## 三、方法

### 3.1 Z-score 定義

```
per-sample baseline = CovM 在 (sample, mode) × non-LOH regions 的 mean, std
CovM_zscore = (CovM - mean) / max(std, 1e-6)
```

- **TO mode non-LOH**：`~to_loh_bed_hit`
- **Paired mode non-LOH**：`~tool_potential_loh AND ~core_loh_like`

### 3.2 Z3 定義（沿用 Step 3）

```
Z3 = loh_flag AND (caller_af < 0.1 OR caller_af > 0.9) AND (HPFineNGroups <= 1)
  其中 loh_flag TO=to_loh_bed_hit, paired=(tool_potential_loh OR core_loh_like)
```

### 3.3 判定

- **z_extreme** = `covm_zscore > 2`（對稱尾部 ~2.3%）
- HCC1954 Z3 × chr5/8/17 × z_extreme × FP 的 overlap 作為「R6 驗證訊號」

---

## 四、主結果

### 4.1 Per-sample non-LOH CovM baseline

| Sample | Mode | Mean | Std | Median | Q25 | Q75 | N |
|--------|------|------|-----|--------|-----|-----|---|
| COLO829 | paired | 0.417 | 0.131 | 0.413 | 0.333 | 0.493 | 29,647 |
| COLO829 | to | 0.427 | 0.138 | 0.413 | 0.333 | 0.507 | 40,035 |
| H1437 | paired | 1.012 | 0.258 | 0.973 | 0.853 | 1.120 | 53,374 |
| H1437 | to | 1.065 | 0.307 | 1.027 | 0.880 | 1.200 | 43,384 |
| H2009 | paired | 1.451 | 0.588 | 1.320 | 1.067 | 1.640 | 94,938 |
| H2009 | to | 1.426 | 0.593 | 1.293 | 1.053 | 1.627 | 103,554 |
| HCC1395 | paired | 1.080 | 0.375 | 1.013 | 0.827 | 1.280 | 15,962 |
| HCC1395 | to | 1.095 | 0.402 | 1.027 | 0.827 | 1.293 | 22,541 |
| HCC1395_DORADO | paired | 1.114 | 0.396 | 1.053 | 0.840 | 1.333 | 16,643 |
| HCC1395_DORADO | to | 1.161 | 0.423 | 1.080 | 0.880 | 1.360 | 22,643 |
| HCC1937 | paired | 1.789 | 0.666 | 1.653 | 1.293 | 2.160 | 5,541 |
| HCC1937 | to | 1.760 | 0.689 | 1.627 | 1.267 | 2.120 | 12,236 |
| **HCC1954** | **paired** | **0.952** | **0.433** | **0.867** | **0.720** | **1.067** | **15,989** |
| **HCC1954** | **to** | **0.978** | **0.692** | **0.880** | **0.720** | **1.080** | **63,081** |

**觀察**：
1. **HCC1954 TO std=0.692 為所有 (sample, mode) 最高**（僅次 HCC1937 的 0.689）→ HCC1954 non-LOH baseline 本身極度不均勻
2. **HCC1954 non-LOH mean=0.978**（低於 H2009 1.426、HCC1937 1.76、HCC1395 1.095）→ HCC1954 在 non-LOH region 整體 coverage 偏低（這與 ~2N pseudo-tetraploid 的分佈假設一致：很多 "normal copy" 的 2N 區域）
3. **COLO829 baseline 特別低（0.427）**→ ONT_R10 平台差異（已知 out-of-scope）

### 4.2 HCC1954 Z3 Validation

| Mode | Metric | Count | Frac |
|------|--------|-------|------|
| TO | total | 2,030 | 100% |
| TO | big3 (chr5/8/17) | **1,714** | **84.4%** |
| TO | FP | 1,929 | 95.0% |
| TO | z_extreme (z>2) | 3 | **0.15%** |
| TO | big3 ∩ z_extreme | 1 | 0.05% |
| TO | FP ∩ big3 | 1,637 | 80.6% |
| TO | FP ∩ z_extreme | 3 | 0.15% |
| TO | **FP ∩ big3 ∩ z_extreme** | **1** | **0.05%** |
| Paired | total | 140 | 100% |
| Paired | big3 | 110 | 78.6% |
| Paired | FP | 3 | 2.1% |
| Paired | z_extreme | 1 | 0.71% |
| Paired | FP ∩ z_extreme | 0 | 0% |

### 4.3 關鍵解讀

**原 Z3 pilot Step 3 觀察**：HCC1954 Z3 FP CovM=0.733 vs TP=0.493（Mann-Whitney p=4.7e-9）

**R6 z-score 重算**：
- FP mean z = (0.733 − 0.978) / 0.692 = **−0.35**（低於 baseline mean）
- TP mean z = (0.493 − 0.978) / 0.692 = **−0.70**（更低於 baseline mean）
- **兩者皆在 z<2 範圍內** — 以絕對 z-score 尺度看 HCC1954 Z3 沒有 CN 極端值現象

**意義**：
- HCC1954 Z3 內部 FP vs TP 的 CovM 差異（0.733 vs 0.493）是**相對差異**，不是**絕對極端**
- 把訊號放入全局 per-sample baseline 框架後，被 HCC1954 本身寬 std（0.692）吸收
- 這與 Z3 pilot S3 strategy 發現「CovM 95th percentile cutoff ΔF1≈+0.0002 微弱」完全一致

---

## 五、合理性驗證

### 5.1 Coverage_Multiple 訊號的兩個成分分離

```
raw CovM = local_coverage / sample_median_coverage
          ≈ (CN_effect × read_depth_variation) / sample_median
```

Per-sample baseline normalize 後：
```
CovM_zscore = (raw_CovM - per_sample_non_LOH_mean) / per_sample_non_LOH_std
            → 去除「sample overall coverage level」
            → 保留「local CN vs 該樣本典型 region」
```

HCC1954 的特殊性（pseudo-tetraploid + HER2+ 複雜 CNV）表現為：
- **per-sample non-LOH mean 偏低**（0.978）— 因 ~2N 區域多
- **per-sample non-LOH std 偏高**（0.692）— 因局部大範圍 CN 變化劇烈
- LOH region（Z3）CovM 落在此分佈中下段，看不出 "extreme high"

### 5.2 為何 Step 3 Mann-Whitney 仍顯著（p=4.7e-9）

- Mann-Whitney 是 rank-based，對**分布中心偏移**敏感，不需要絕對極端值
- FP 0.733 vs TP 0.493 在 HCC1954 內部是 middle-shift，仍具 rank 差異
- R6 z-score 的 threshold-based 判定（z>2）與 rank-based 統計差異源自不同尺度

### 5.3 與 Step 5 ceiling 校準 consistent

Z3 pilot Step 5 顯示：S2 blacklist（whole chr5/8/17 ∩ Z3）達 ceiling 87%，意味 refined CovM strategy 幾乎無剩餘訊號可利用。R6 的 NEGATIVE 結果正是這一 ceiling observation 的**per-sample normalization 解釋**。

---

## 六、對相關結論的影響

| 結論 / 聲稱 | R6 前 | R6 後 |
|-----------|-------|-------|
| Opus 4.7 plan 假設 4「Coverage_Multiple ≈ CN」 | ❌ 未驗證 | **❌ 否證**（per-sample normalize 後訊號 dissolve） |
| Z3 pilot S3 strategy（CovM 95th cutoff）訊號微弱 | ⭐4（pilot 已定） | ⭐5（R6 獨立機制解釋） |
| Z3 pilot S1/S2 染色體 arm blacklist HCC1954-only CONDITIONAL | ⭐4 | ⭐4（不變，R6 正交） |
| B.2 批次 1 Coverage_Multiple 用作 CN tier 代理 | caveat noted | caveat 升級為「明確受 per-sample baseline 吸收」 |
| B.2 R4 (LOH 類型分層：cnLOH vs deletion-LOH 用 CovM) | ⭐3 | ⭐3（限制保留，但 CN≥3 tier 判定需 R17 CNV caller 整合而非純 CovM） |

---

## 七、對 R1 HCC1395 Normal BAM pilot 的啟示

R6 否證了 Coverage_Multiple 作為獨立 CN 訊號源 → **強化 Normal BAM annotation 的必要性**：

- Normal BAM coverage 可作為真正的 **sample-matched baseline**（不需 per-sample self-reference）
- Δ_coverage = tumor_CovM / normal_CovM 理論上比 raw CovM 更接近精確 CN
- **R1 pilot 的一個 side-benefit**：同時產出 Normal-based CN ratio 作為 R17（未來 HCC1954 Normal BAM pilot）的方法學 template

---

## 八、產物清單

### 腳本
- `research/z3_internal_feature_exploration/scripts/r6_covm_zscore_normalization.py`

### 數據
- `research/z3_internal_feature_exploration/data/r6_covm_per_sample_stats.tsv` — 14 rows (7 samples × 2 modes)
- `research/z3_internal_feature_exploration/data/r6_covm_zscore_regions.tsv.gz` — 748,391 rows (全量 z-score)
- `research/z3_internal_feature_exploration/data/r6_hcc1954_z3_extreme_validation.tsv` — HCC1954 Z3 驗證表

### 報告
- 本文件 `docs/experiments/in_progress/2026/04/20260419_Coverage_Multiple_zscore_normalization_01.md`

---

## 九、結論與後續

### 9.1 結論

- **R6 判定 CONDITIONAL NEGATIVE**：Coverage_Multiple 經 per-sample non-LOH baseline z-score normalize 後，原 Step 3 觀察到的 HCC1954 Z3 FP CovM 偏高訊號**消失**（z_extreme 僅 0.15% TO）
- 真正的 HCC1954 Z3 FP 富集訊號是**染色體位置**（chr5/8/17，84.4%），而非 CovM 極端值
- Opus 4.7 plan 假設 4 ❌ 否證；結論 11（LOH×AF×Methylation POSITIVE）受此影響，Coverage_Multiple 的 CN tier 使用需加註 caveat

### 9.2 後續行動

- **R1 軌不受阻**：R6 的 NEGATIVE 反而**強化** Normal BAM annotation 的必要性（提供 sample-matched CN baseline）
- **穩定度表更新**：`06_結論穩定性審查.md` 新增結論 19「Coverage_Multiple 非獨立 CN proxy」⭐4（per-sample normalize 驗證）
- **B.2 R4 後續**：若要深入 cnLOH vs deletion-LOH 分層，需 **R17（未來）CNV caller 整合**（Delly/Manta/sequenza），不能只靠 CovM
- **Z3 pilot 不需重做**：S3 strategy 微弱訊號在 R6 得到機制解釋，原 CONDITIONAL 判定維持
