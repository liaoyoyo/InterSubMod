---
title: 切群重設計 全基因組 TP+FP — 五類比例 + TP/FP 特異性對比
date: 2026-06-22
status: in_progress
tier: 4
sample: HCC1395 tumor-only（單樣本 ⭐2-3 characterization）
scope: 全基因組 22 autosome TP 30077 + FP 4659 = 34736 位點（big7 本機 binary 重跑）
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/cluster_redesign_wg_summary.json
build_commit: f429313
observation_standard: true
---

# 切群重設計 全基因組 TP+FP — 五類 + TP/FP 特異性

> **🔴 核心結論**：germline 對齊的 cis-ASM 結構（CONFIRMED）**somatic-相關**（TP 富集 3.3×）；但「不對齊 germline 的真實結構」（REAL_NOVEL = subclone 候選）**非 somatic 特異**（TP 0.65×，FP 反而更多）→ **單樣本無法宣稱為 subclone 訊號**。
> **圖**：`figs_dashboard/wg_tpfp_contrast.png`

## 0. 規模 + 方法
- big7 本機 binary（TP 30490→30077 splittable + FP 4842→4659）→ 五類切群（四閘 + null Rnull=15 + Pool24）→ 50 min、0 error。
- 五類：CONFIRMED / NEAR_CONFIRMED / REAL_NOVEL(subclone候選) / REAL_DIFFUSE / NO_CLEAR。

## 1. 五類 TP vs FP（全 L1 from summary.json）
| 類 | TP n | TP % | FP n | FP % | **富集 TP/FP** |
|---|---:|---:|---:|---:|---:|
| **CONFIRMED**（對齊 germline）| 6198 | 20.6% | 292 | 6.3% | **3.29×** ✅ |
| NEAR_CONFIRMED | 118 | 0.4% | 9 | 0.2% | 2.03× |
| **REAL_NOVEL**（subclone 候選）| 3206 | 10.7% | 761 | 16.3% | **0.65×** 🔴 |
| REAL_DIFFUSE | 12479 | 41.5% | 1998 | 42.9% | 0.97× |
| NO_CLEAR（真無結構）| 8076 | 26.9% | 1599 | 34.3% | 0.78× |

- coarse 有 germline 骨幹：**TP 22.9% vs FP 8.0% = 2.86×**。
- subclone 候選合計（NOVEL+DIFFUSE）：TP 52.1% vs FP 59.2% = **0.88×（不富集）**。

## 2. 解讀（誠實，論文口徑）
1. **CONFIRMED（germline-cis）TP 富集 3.3×** = 真 somatic 位點旁的甲基結構，凡能對齊 germline 軸者，確實 **somatic-位點相關**（cis-ASM at somatic loci）。這是真實正向訊號。
2. **🔴 REAL_NOVEL（subclone 候選）TP 貧化（0.65×），FP 反更多** = **「不跟 germline 走的真實結構」不是 somatic 特異** → 在 FP（artifact）位點更常見，**很可能是 artifact/比對假象驅動的 read clustering，非 subclone**。
3. **REAL_DIFFUSE ≈ 相等（0.97×）** = diffuse 結構 TP/FP 無差 = 技術/背景 cis-ASM。
4. **NO_CLEAR TP 略低（0.78×）** = FP 更多無結構（合理，FP 品質低）。

→ **單樣本 HCC1395 tumor-only：唯一 somatic-相關的是 germline 對齊 cis-ASM（CONFIRMED）；subclone 候選結構（REAL_NOVEL）非 somatic 特異、無法當 subclone 訊號**。

## 3. 對論文的意涵
- 呼應既有 **tumor-only 軸 NEGATIVE** + 「ASM 真實但不判別 TP/FP」+ baseline-dependence。新方法的 read-level 結構偵測**沒有逃脱此 confound**：unaligned 結構在 FP 一樣多（甚至更多）。
- **強化**：單樣本談 subclone 需 **normal cis-control 排 cis-ASM + 多樣本/COLO829**；REAL_NOVEL 標「候選」是正確的（全基因組證實非 somatic 特異）。
- **正向可用**：CONFIRMED（germline-cis）的 3.3× TP 富集是真實 characterization 訊號（cis-ASM 在 somatic 位點），可作論文 characterization 面材料（非 subclone reconstruction）。

## 4. 驗證表
| 數字 | 值 | 來源 key | L |
|---|---|---|---|
| 總位點 | 34736（TP 30077 / FP 4659）| by_set | L1 |
| CONFIRMED 富集 | 3.29×（TP 20.6/FP 6.3%）| fine_TP/fine_FP.CONFIRMED | L1 |
| REAL_NOVEL 富集 | 0.65×（TP 10.7/FP 16.3%）| fine_TP/fine_FP.REAL_NOVEL | L1 |
| coarse 骨幹 | TP 22.9/FP 8.0% = 2.86× | records coarse_conf | L1 |
| 耗時 / error | 3031s / 0 | summary.elapsed_s | L1 |

## 5. 限制
- **單樣本 ⭐2-3**：富集是 TP/FP 對比，非 normal cis-control（cis-ASM vs subclone 仍未分離）。
- FP 集合本身有限（4659）+ FP 定義依賴 caller；FN 未跑。
- null Rnull=15（提速）vs 小規模 25；門檻同小規模（excess 0.10 / near 0.08 / gap big_gap / e≥5）。
- old(maxclust) vs new delta 未逐位點存（records 只存 new 五類）；如需另跑。

## 6. 重生
```bash
cd InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot
python3 scripts/cluster_redesign_wg.py        # 全 22 chr TP+FP（big7 本機）→ wg_records/summary.json
python3 scripts/plot_wg_tpfp.py               # → wg_tpfp_contrast.png
```

關聯 [[project_tumor_only_axis_negative_subclone_classification]] [[project_ism_complete_tpfpfn_existence_cis]] [[feedback_baseline_dependence_not_result]] [[project_cluster_redesign_three_gate]]。
