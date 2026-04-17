---
title: P4 Second-Hit Order Pilot — CONDITIONAL NEGATIVE
date: 2026-04-17
status: in_progress
verdict: CONDITIONAL_NEGATIVE
ism_goal: 目標 3 (二次打擊順序) — 單 region 特徵無法區分
related:
  - research/P4_second_hit_order_pilot/
  - project_loh_subclone_af_methylation_positive.md (B.2 不同比較)
---

# P4 Second-Hit Order Pilot — CONDITIONAL NEGATIVE

## 1. 目的

ISM 目標 3 = 推斷二次打擊順序。P4 小規模 pilot 驗證：**區分「methylation silencing 先於 LOH」vs「LOH 先於 methylation silencing」是否可在 region-level AF × methylation 聯合分佈上做到**。

## 2. 方法

- 資料：master dataset 748K regions, 6 樣本 TO mode
- 過濾：NR≥40, TP+FP, 有 AF/AlleleDelta/HPFineNGroups
- 分層：AF bin × LOH flag (`to_loh_bed_hit`)
- 核心對比：同 AF bin 下 LOH vs NonLOH 的 AlleleDelta Bimodality Coefficient + HPFineNGroups
- Bimodality Coefficient (SAS): BC = (skew² + 1) / (kurt + 3(n-1)²/((n-2)(n-3)))
  - BC > 0.556 代表 bimodal
- 假設：
  - 若 methylation-first：intermediate AF 的 LOH region 應 AlleleDelta **分散**（biallelic）
  - 若 LOH-first：intermediate AF 的 LOH region 應 AlleleDelta **集中**（monoallelic）
  - 兩者差異應反映在 BC 與 NGroups

## 3. 結果

### 3.1 Bimodality Coefficient: 差異微小

| AF bin | sample | LOH BC | NonLOH BC | Δ |
|---|---|---|---|---|
| 0.3-0.5 | COLO829 | 0.739 | 0.668 | +0.071 |
| 0.3-0.5 | H1437 | 0.392 | 0.447 | -0.054 |
| 0.3-0.5 | H2009 | 0.749 | 0.716 | +0.033 |
| 0.3-0.5 | HCC1395 | 0.548 | 0.534 | +0.013 |
| 0.3-0.5 | HCC1937 | 0.464 | 0.407 | +0.058 |
| 0.3-0.5 | HCC1954 | 0.538 | 0.564 | -0.027 |

Mean |Δ BC| = **0.043**（遠低於 POSITIVE 門檻 0.15）

→ AlleleDelta 分佈在 LOH/NonLOH 都偏 unimodal（BC 多介於 0.4-0.75，邊界 0.556），**不顯現預期的 biallelic vs monoallelic 分離**

### 3.2 HPFineNGroups: LOH 系統性偏低（非 order 訊號）

| sample | AF 0.3-0.5 ΔNGroups (LOH - NonLOH) |
|---|---|
| COLO829 | -0.448 |
| H1437 | -0.688 |
| H2009 | -0.952 |
| HCC1395 | -0.628 |
| HCC1937 | -0.517 |
| HCC1954 | -0.588 |

Mean |Δ NGroups| = **0.637**（遠超門檻，7/6 樣本方向一致）

**但這不是 order 訊號**：
- LOH region 單 haplotype → possible NGroups values 本身較小
- 這是 **LOH 機制副作用**（haplotype 減少），不是 "哪個先發生"

### 3.3 與 B.2 LOH POSITIVE 的關係

- B.2 結論：**Intermediate AF** within LOH vs **Extreme AF** within LOH → ΔNG=+0.705
- P4 比較：**LOH** vs **NonLOH** within **same AF bin** → ΔNG=-0.637
- **不矛盾**：B.2 回答「LOH 內 subclonal vs clonal」；P4 回答「LOH vs Non-LOH 同 AF」
- 兩者揭示不同機制面向

## 4. Verdict: CONDITIONAL NEGATIVE

- 目標訊號（Bimodality）**未達門檻**
- 次要訊號（NGroups）**是 LOH 副作用**，非 order marker
- **單 region 特徵 summary 無法區分二次打擊順序**

## 5. 對 ISM 目標 3 的意涵

### 方法論結論
推斷 two-hit order **需要 per-read epigenotype**：
- 每條 read 的 methylation pattern 跟 variant allele 的 cis/trans 關係
- phased haplotype × methylation sub-pattern
- **對應 ISM 目標 1（per-CpG 多標籤）**

### 目標依賴關係
- **目標 3 ← 依賴 ← 目標 1**
- Phase 2 路線圖調整：目標 3 不應在目標 1 之前做 pilot
- 現有特徵集（HPFineNGroups + AlleleDelta + ASM summaries）對 order 無力

### 下一步建議
- 目標 3 工作**暫緩**，等目標 1（per-CpG 多標籤 + per-read epigenotype）有進度再做
- 現階段優先：HPFineNGroups × LOH × AF 更細緻的特徵化（已確認有訊號）、P0 修復

## 6. 記錄與數據

- 腳本：`research/P4_second_hit_order_pilot/scripts/01_af_methylation_joint_distribution.py`
- 原始數據：`research/P4_second_hit_order_pilot/data/{af_loh_methylation_joint_stats,loh_vs_nonloh_delta}.tsv`
- 視覺化：`research/P4_second_hit_order_pilot/figures/01_af_loh_joint_distribution.png`
- Manifest: `research/P4_second_hit_order_pilot/manifest.yaml`
