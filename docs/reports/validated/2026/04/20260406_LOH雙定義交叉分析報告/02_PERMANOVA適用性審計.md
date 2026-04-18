<!--
建立時間: 2026-04-06
目標: Step 3.1 PERMANOVA 在 LOH 區域的適用性審計
處理範圍: 7 samples × 2 modes, HP/Allele PERMANOVA
關聯檔案:
  - scripts/analysis/build_permanova_loh_audit.py
-->

# Step 3.1: PERMANOVA 適用性審計

## 背景

ISM 使用 PERMANOVA 檢驗甲基化模式在分組（HP/Allele）之間的差異。Wave 1 需確認 PERMANOVA 在 LOH 區域是否仍然有效 — LOH 造成 HP 嚴重偏斜，可能使 HP 維度的 PERMANOVA 前提（每組 ≥ 3 reads）不成立。

---

## 核心結果

| 模式 | LOH 狀態 | N | HP Valid (%) | HP Sig (%) | min(HP)<3 (%) | Allele Valid (%) | Allele Sig (%) |
|------|----------|---|-------------|-----------|--------------|-----------------|---------------|
| PAIRED | LOH | 96,605 | **6.5** | 1.5 | **93.5** | **50.7** | 28.0 |
| PAIRED | Non-LOH | 232,094 | 86.6 | 46.6 | 13.4 | 98.5 | 52.2 |
| TO | LOH | 175,828 | **5.4** | 1.1 | **94.6** | **58.6** | 31.0 |
| TO | Non-LOH | 243,864 | 98.9 | 53.0 | 1.1 | 99.8 | 53.3 |

![HP PERMANOVA valid rate](figures/w1p_01_hp_permanova_valid_rate.png)

![min(HP) group distribution](figures/w1p_02_min_hp_group_distribution.png)

---

## HP PERMANOVA 區分力 (AUC)

| 條件 | HP P → TP/FP AUC | 解讀 |
|------|------------------|------|
| Paired LOH | **0.237** | 反向（比隨機差） |
| Paired Non-LOH | 0.573 | 微弱 |
| TO LOH | **0.461** | 低於隨機 |
| TO Non-LOH | 0.479 | 低於隨機 |

![HP merged sig rate](figures/w1p_03_hp_merged_sig_rate.png)

![Pseudo-F distribution](figures/w1p_04_pseudo_f_distribution.png)

![PERMANOVA p-value CDF](figures/w1p_05_permanova_p_cdf.png)

![TP/FP pseudo-F in LOH](figures/w1p_06_tp_fp_pseudo_f_loh.png)

![Per-sample valid rate heatmap](figures/w1p_07_per_sample_valid_heatmap.png)

![Summary table](figures/w1p_08_summary_table.png)

---

## 判定

### J2: HP PERMANOVA 在 LOH 內不可用（確定）

| 判定標準 | 閾值 | 實際值 | 結果 |
|---------|------|--------|------|
| HP valid rate < 50% | < 50% | **5.4-6.5%** | 遠低於閾值 |
| min(HP) < 3 比例 > 60% | > 60% | **93.5-94.6%** | 遠超閾值 |
| HP AUC < 0.55 | < 0.55 | **0.237-0.479** | 全面低於閾值 |

### J4: Allele 維度在 LOH 內遠優於 HP

| 指標 | HP | Allele | 倍率 |
|------|-----|--------|------|
| Valid rate (TO LOH) | 5.4% | **58.6%** | **~10x** |
| Sig rate (TO LOH) | 1.1% | **31.0%** | **~28x** |

---

## 操作條件與限制

- LOH 定義: Potential_LOH (HP_Ratio based)
- HP PERMANOVA valid = min(HP1N, HP2N) ≥ 3
- Allele valid rate 高但 AUC 尚未計算（→ Step 3.3 驗證）
- Pseudo-F median 在 LOH 內為 1.66-1.93（vs Non-LOH 4.59-4.84）
