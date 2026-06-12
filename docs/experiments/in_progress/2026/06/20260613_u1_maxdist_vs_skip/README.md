---
title: U1 MAX_DIST vs SKIP 全量對比（BERNOULLI ±5000）
date: 2026-06-13
type: experiment
status: DONE（SKIP segfault，對比未完成；揭露新 bug）
data_sources:
  - HCC1395 paired tumor/normal, BERNOULLI ±5000, 全量 TP+FP
---

# U1: 距離不足 pair 的 MAX_DIST vs SKIP（你最關心的分群誤判驗證）

## L0
- 🔴 **SKIP segfault（新 bug U1-bug）** → SKIP 目前不可用，對比無法完成
- ✅ **MAX_DIST BERNOULLI ±5000 運作良好**：Noise 僅 9.4%（vs 我前面 NHD ±1000 的 26.4%）→ 你擔心的「MAX_DIST 大量分群誤判」證據不支持

## 數據（HCC1395 paired，全量）
| 配置 | n | 保留率 | Strong | **Noise** | PermanovaValid |
|---|---|---|---|---|---|
| NHD ±1000（前面觀察）TP | 30490 | 2.42% | 6.5% | **26.4%** | — |
| **BERNOULLI ±5000 MAX_DIST TP** | 30490 | **4.18%**(1276) | **16.8%**(5133) | **9.4%**(2856) | 6971 |
| **BERNOULLI ±5000 MAX_DIST FP** | 4842 | 1.67%(81) | 13.4%(651) | 13.8%(666) | 947 |
| BERNOULLI ±5000 SKIP TP/FP | — | **segfault(exit139)** | — | — | — |

FP 保留 1.67% < TP 4.18% → 方向正確。

## 🔴 U1-bug: SKIP segfault（待修）
- 現象：`--nan-distance-strategy SKIP` → 第一個 region 就 exit 139（SIGSEGV）。
- 根因：SKIP 把不重疊 pair 距離設 **NaN**（DistanceMatrix.cpp:413）→ 下游 **clustering(UPGMA) / PERMANOVA 沒處理 NaN 距離矩陣 → segfault**。
- 修法（走 methodology-audit→cpp-change）：clustering/PERMANOVA 前剔除 NaN-heavy read（或 NaN-aware linkage）。`StructureTest` 已有 `filter_reads_for_complete_matrix`，但 clustering 路徑可能未過此 filter。

## 結論
- **MAX_DIST + BERNOULLI ±5000 = 目前唯一可用且表現好的配置**（Noise 低、Strong 多、TP>FP）。
- **±5000 BERNOULLI 量化證實 > ±1000 NHD**（Noise 26%→9%、Strong 翻倍、保留率更高）。
- SKIP 是否更優 → **必須先修 U1-bug 才能驗證**。
