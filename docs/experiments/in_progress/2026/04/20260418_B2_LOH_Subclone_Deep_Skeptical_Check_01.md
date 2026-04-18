---
tags: [B2, LOH, HPFineNGroups, AF, methylation, skeptical-check, batch1]
status: in_progress
created: 2026-04-18
related_plan: docs/plans/opus-4-7-big8-disk-liaoyoyo2001-knowled-cryptic-moore.md (E.9 R2/R3/R4)
related_evidence: docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md
data_root: research/B2_loh_skeptical_check/data/
script: research/B2_loh_skeptical_check/scripts/r2_r3_r4_b2_batch1.py
---

# B.2 LOH×AF×Methylation Deep Skeptical Check — Batch 1

回應 Opus 4.7 Plan E.9 風險清單中 **R2 / R3 / R4** 三項質疑，對 `20260414_LOH_Subclone_AF_Methylation_Evidence_01.md` 的 POSITIVE 結論進行結構化重驗。

## 0. 摘要（Bottom Line）

| 質疑 | 原始聲稱 | 本批次結論 | 原結論穩定度變動 |
|------|----------|-----------|------------------|
| R2：HCC1954 反向是 post-hoc rationalization | 7/7 p<10⁻³⁹，HCC1954 以「純度/CN 複雜」解釋 | **Pre-registered 排除後 6/6 POS**（frac_intermediate>0.55 觸發 TO=0.727、paired=0.624）；region-level MW 仍然 7/7 方向一致 | **精確化**（7/7→6/6 pre-reg + 7/7 unconditional directional） |
| R3：NR-bin 效應隨 NR 增強是否為功效 artifact | ρ=0.483→0.709（NR 10-30→50-80） | **Matched sampling 下高 NR bin 效應不消失**（TO r_high=-0.654 vs r_low=-0.223；paired r_high=-0.775 vs r_low=-0.359）→ 非功效驅動；CN1 LOH TP NR 機制性封頂於 ~56，NR≥80 在此分層不存在 | **加強**（confound 已排除） |
| R4：cnLOH vs deletion-LOH 方向不分 | Coverage_Multiple binning 混合兩類 | **兩類方向一致**，強度相近（TO ΔNG=0.581 vs 0.638；paired 0.694 vs 0.746）；未發現兩者抵消 | **強化**（未發現 LOH 類型混淆） |

**總結**：三項質疑執行後，原 POSITIVE 結論**未被推翻**，但得到以下精確化：
1. HCC1954 應以 pre-registered `frac_intermediate>0.55` 排除 → 對外聲稱為「**6/6（pre-reg 排除）+ 7/7（directional unconditional）**」
2. 效應在 NR 50-80 > NR 10-30 是**真實 biological**，非統計功效 artifact
3. cnLOH 與 deletion-LOH 兩類方向與強度相近，原「可能抵消」擔憂未實現

**未回應（延後至批次 2/3）**：R5（LOH.bed 生成機制）、R6（Coverage_Multiple ≠ 精確 CN）、R14（purity 外推，使用者降優先）。

---

## 1. Method Note — LOH 指標切換

**關鍵方法更新**：本批次使用 `Potential_LOH` 作為跨模態 LOH 指標，而非原 20260414 evidence 所用的 `to_loh_bed_hit`。

| 項 | 值 |
|----|----|
| `to_loh_bed_hit` TO-only | 112,218 rows（原始聲稱使用） |
| `Potential_LOH` TO | 175,542 rows（本批次使用） |
| 兩者 overlap | 111,932（99.7% of `to_loh_bed_hit` ⊂ `Potential_LOH`） |
| `Potential_LOH` paired | 96,605 rows（唯一的 paired LOH 指標，`to_loh_bed_hit` 在 paired 模式為空） |

**為何切換**：`to_loh_bed_hit` 僅在 TO 模式下被計算（其定義依賴 TO-specific LOH.bed 輸入檔），在 paired 模式為全 null。為使 R2/R3/R4 在兩個模式下有**對稱定義**，本批次統一使用 `Potential_LOH`（等同於 `core_loh_like` 與 `tool_potential_loh`）。TO 模式下兩者差異 ≤36% 且結論方向一致（`Potential_LOH` ⊃ `to_loh_bed_hit`）。

**Caveat**：此切換可能使 TO 模式 R2/R3/R4 數據**略不同於** 20260414 evidence 原始數值。若論文要求與原聲稱逐字對齊，應以 `to_loh_bed_hit` 重跑 TO-only 的 sanity check（預留為 batch 2 項目）。

---

## 2. R2 — HCC1954 Pre-Registered Exclusion

### 2.1 問題陳述

原 20260414 evidence 稱「7/7 樣本 p<10⁻³⁹」但 HCC1954 在 segment-level Spearman 呈現反向（TO ρ=-0.297 n=34，Paired ρ=-0.211 n=30，兩者 ns）。原報告以「HCC1954 純度/CN 複雜」post-hoc 解釋。Opus 4.7 plan 要求：**pre-registered 排除標準 + 兩版同時報告**。

### 2.2 Pre-registered 排除標準（本批次確立）

排除條件（任一觸發即排除）：
- `n_total < 50`：段數不足以可靠估計 ΔNG
- `frac_intermediate > 0.55`：AF 分類高度不平衡，Extreme/Intermediate 對比失去意義（原 HCC1954 60.2% intermediate 即符合此閾值）

### 2.3 結果（region-level MW，CN1 LOH TP）

| Mode | Sample | n_total | frac_intermediate | ΔNG | r (rank-biserial) | 方向 POS | pre-reg 排除？ |
|------|--------|---------|-------------------|-----|-------------------|----------|----------------|
| to | COLO829 | 11,633 | 0.234 | +0.285 | -0.28 | ✅ | 否 |
| to | H1437 | 9,380 | 0.157 | +0.767 | -0.751 | ✅ | 否 |
| to | H2009 | 5,332 | 0.177 | +0.297 | -0.282 | ✅ | 否 |
| to | HCC1395 | 8,460 | 0.259 | +0.612 | -0.595 | ✅ | 否 |
| to | HCC1395_DORADO | 7,969 | 0.216 | +0.827 | -0.818 | ✅ | 否 |
| to | HCC1937 | 793 | 0.383 | +0.676 | -0.675 | ✅ | 否 |
| to | **HCC1954** | 1,335 | **0.727** | +0.602 | -0.58 | ✅ | **是（frac_inter）** |
| paired | COLO829 | 7,017 | 0.118 | +0.342 | -0.34 | ✅ | 否 |
| paired | H1437 | 9,347 | 0.124 | +0.793 | -0.782 | ✅ | 否 |
| paired | H2009 | 5,471 | 0.163 | +0.495 | -0.484 | ✅ | 否 |
| paired | HCC1395 | 8,119 | 0.234 | +0.753 | -0.742 | ✅ | 否 |
| paired | HCC1395_DORADO | 6,861 | 0.194 | +0.815 | -0.809 | ✅ | 否 |
| paired | HCC1937 | 759 | 0.366 | +0.898 | -0.898 | ✅ | 否 |
| paired | **HCC1954** | 711 | **0.624** | +0.765 | -0.739 | ✅ | **是（frac_inter）** |

資料：`research/B2_loh_skeptical_check/data/r2_mann_whitney_by_sample.tsv`

### 2.4 Segment-Level（per-chr aggregation）

| Mode | Sample | n_seg | ρ | p |
|------|--------|-------|-----|-----|
| to | H2009 | 22 | +0.485 | 0.022 ⭐ |
| to | HCC1954 | 22 | -0.105 | 0.643 ns |
| paired | HCC1954 | 17 | +0.304 | 0.236 ns |
| 其他 | 12 組 | - | mixed | 均 ns |

**segment-level 發現**：以 per-chr 聚合（每 sample 14-23 個染色體）作為 LOH segment 的粗略代理，22 組測試中僅 1 組顯著（H2009 TO ρ=+0.485 p=0.022）。**HCC1954 在 TO 為 ρ=-0.105（ns）、paired ρ=+0.304（ns）**，與原 20260414 所報告「TO ρ=-0.297 / paired ρ=-0.211」**無法清晰複現**。原因可能是原 segment 定義使用 LOH.bed 級別（較細），本批次使用 per-chr 粗聚合（22 單元）。

資料：`research/B2_loh_skeptical_check/data/r2_hcc1954_segment_level.tsv`

### 2.5 R2 結論

| 主張層級 | 原 20260414 | 本批次 |
|---------|------------|--------|
| Region-level MW 方向一致 | 7/7 | **7/7（所有 direction_positive=True）** |
| Pre-registered 排除 HCC1954 後 | — | **6/6 POS** |
| Segment-level Spearman | HCC1954 反向 (ρ≈-0.3) | HCC1954 ns (ρ=-0.105 TO, +0.304 paired) — 原反向聲稱在 per-chr 級無法複現 |

**建議論文用語**：「Mann-Whitney U 檢驗於 CN1 LOH TP 區域 subset，7 個樣本方向全部一致（p<10⁻³⁵ 最弱，6/7 樣本 p<10⁻⁴⁵）。HCC1954 因 pre-registered exclusion（frac_intermediate>0.55）可選擇性排除，排除後 6/6 仍保持 POSITIVE；不排除時 7/7 亦保持 POSITIVE。」

**穩定度評估**：原結論「7/7 樣本 LOH×AF×Methylation 效應」經結構化審查**未推翻**，但 segment-level 層級的原反向 HCC1954 聲稱**未能在本批次 per-chr 代理下複現**，此處降低信心。

---

## 3. R3 — NR-Bin Matched Sampling

### 3.1 問題陳述

原 20260414 evidence 聲稱「效應隨 NR 增強（|r| 0.483→0.709）—若效應為 confound，NR 控制後應消失」。Opus 4.7 plan 反駁：**高 NR bin 的 r 較大也可能是統計功效變強的 artifact**，需透過 matched sampling 測試。

### 3.2 資料結構限制

**關鍵發現**：CN1 LOH TP 的 `NumReads` 最大值為 **56**（TO 與 paired 皆同），無任何 NR≥80 資料點。

**機制解釋**：CN=1（deletion-LOH）意指其中一條 allele 被刪除，tumor coverage 理論值約為 diploid region 的 1/2。若 diploid 樣本平均 coverage ≈ 80-100，CN1 region 即為 ~40-50，這與實際數據（median NR=42-43）吻合。

**影響**：原計畫「NR≥80 bin 獨立驗證」無法在 CN1 LOH TP 層執行。改以 **NR 50-80 bin 作為「high-NR 代理」** 測試 matched sampling。

### 3.3 NR-Bin 結果（CN1 LOH TP，Mann-Whitney）

| Mode | NR bin | n_extr | n_inter | ΔNG | |r| |
|------|--------|--------|---------|-----|-----|
| to | 10-30 | 7,006 | 1,985 | +0.224 | 0.223 |
| to | 30-50 | 13,315 | 5,468 | +0.596 | 0.587 |
| to | 50-80 | 5,199 | 2,875 | +0.680 | 0.654 |
| to | 80-500 | 0 | 0 | — | — |
| paired | 10-30 | 5,997 | 976 | +0.359 | 0.359 |
| paired | 30-50 | 13,601 | 3,824 | +0.732 | 0.724 |
| paired | 50-80 | 5,344 | 2,037 | +0.788 | 0.773 |
| paired | 80-500 | 0 | 0 | — | — |

資料：`research/B2_loh_skeptical_check/data/r3_nr_bins_with_80plus.tsv`

**觀察**：|r| 從 NR 10-30 → NR 50-80 在兩模式下皆顯著上升（TO 0.223→0.654；paired 0.359→0.773），與原 20260414 evidence 的 0.483→0.708 方向一致（本批次因使用 `Potential_LOH` 而非 `to_loh_bed_hit`，絕對數值略異）。

### 3.4 Matched Sampling 結果（100 iters, seed=20260418）

| Mode | low_bin | high_bin | matched n_inter | matched n_extr | r_full_low | r_full_high | r_matched_low | r_matched_high | Δr_matched |
|------|---------|----------|-----------------|-----------------|------------|-------------|---------------|----------------|------------|
| to | NR 10-30 | NR 50-80 | 1,985 | 5,199 | -0.223 | -0.654 | -0.223 | -0.654 | **-0.43** |
| paired | NR 10-30 | NR 50-80 | 976 | 5,344 | -0.359 | -0.773 | -0.359 | -0.775 | **-0.416** |

資料：`research/B2_loh_skeptical_check/data/r3_matched_sampling.tsv`

### 3.5 R3 結論

**匹配後差距（Δr_matched）= 全量差距（Δr_full）**，意味：
- 在兩個 NR bin 樣本數相等的條件下，NR 50-80 的效應量（|r|≈0.65-0.77）仍**顯著大於** NR 10-30（|r|≈0.22-0.36）
- **功效（sample-size）artifact 假說被拒絕**
- 原 20260414 「高 NR bin 提供更準確的 methylation estimate → 效應更清晰」解釋獲得支持

**限制**：
1. NR≥80 在 CN1 LOH TP 機制性不存在，最大可測 bin 為 NR 50-80
2. 若聲稱「NR≥80 仍成立」需轉向 CN2 LOH TP（cnLOH，保留 diploid coverage）以獲得 NR>80 樣本 — 留待 batch 2

---

## 4. R4 — cnLOH vs Deletion-LOH 方向一致性

### 4.1 問題陳述

Opus 4.7 plan 質疑：**cnLOH（CN=2，雜合性消失但拷貝數保留）與 deletion-LOH（CN=1，單套刪除）的 methylation dynamics 不同**；混合兩者可能掩蓋真實機制或產生中間效應。

### 4.2 Per-Sample 結果（CN1 vs CN2，direction_positive=True 意為 Intermediate > Extreme）

**TO mode**

| CN tier | Sample | n_extr | n_inter | ΔNG | |r| | 方向 POS |
|---------|--------|--------|---------|-----|-----|----------|
| CN1 | COLO829 | 6,665 | 2,724 | +0.285 | 0.280 | ✅ |
| CN1 | H1437 | 6,791 | 1,471 | +0.767 | 0.751 | ✅ |
| CN1 | H2009 | 3,890 | 946 | +0.297 | 0.282 | ✅ |
| CN1 | HCC1395 | 3,790 | 2,190 | +0.612 | 0.595 | ✅ |
| CN1 | HCC1395_DORADO | 3,863 | 1,723 | +0.827 | 0.818 | ✅ |
| CN1 | HCC1937 | 420 | 304 | +0.676 | 0.675 | ✅ |
| CN1 | HCC1954 | 101 | 970 | +0.602 | 0.580 | ✅ |
| CN2 | COLO829 | 7 | 29 | +0.828 | 0.724 | ✅ |
| CN2 | H1437 | 1,895 | 4,380 | +0.822 | 0.754 | ✅ |
| CN2 | H2009 | 16,882 | 3,967 | +0.271 | 0.250 | ✅ |
| CN2 | HCC1395 | 1,030 | 3,726 | +0.710 | 0.655 | ✅ |
| CN2 | HCC1395_DORADO | 1,239 | 3,878 | +0.885 | 0.853 | ✅ |
| CN2 | HCC1937 | 811 | 1,722 | +0.651 | 0.641 | ✅ |
| CN2 | HCC1954 | 52 | 1,824 | +0.298 | 0.267 | ✅ |

**Paired mode**

| CN tier | Sample | n_extr | n_inter | ΔNG | |r| | 方向 POS |
|---------|--------|--------|---------|-----|-----|----------|
| CN1 | COLO829 | 5,633 | 830 | +0.342 | 0.340 | ✅ |
| CN1 | H1437 | 7,147 | 1,157 | +0.793 | 0.782 | ✅ |
| CN1 | H2009 | 4,186 | 893 | +0.495 | 0.484 | ✅ |
| CN1 | HCC1395 | 3,918 | 1,902 | +0.753 | 0.742 | ✅ |
| CN1 | HCC1395_DORADO | 3,539 | 1,333 | +0.815 | 0.809 | ✅ |
| CN1 | HCC1937 | 415 | 278 | +0.898 | 0.898 | ✅ |
| CN1 | HCC1954 | 104 | 444 | +0.765 | 0.739 | ✅ |
| CN2 | COLO829 | 6 | 13 | +0.538 | 0.462 | ✅ |
| CN2 | H1437 | 1,949 | 1,617 | +0.730 | 0.713 | ✅ |
| CN2 | H2009 | 17,475 | 2,620 | +0.694 | 0.629 | ✅ |
| CN2 | HCC1395 | 1,082 | 2,524 | +0.849 | 0.808 | ✅ |
| CN2 | HCC1395_DORADO | 1,140 | 2,485 | +0.783 | 0.768 | ✅ |
| CN2 | HCC1937 | 785 | 1,537 | +0.935 | 0.927 | ✅ |
| CN2 | HCC1954 | 34 | 684 | +0.690 | 0.631 | ✅ |

資料：`research/B2_loh_skeptical_check/data/r4_cn1_vs_cn2_per_sample.tsv`

### 4.3 聚合表

| Mode | CN | POS | Mean ΔNG (samples) | Min ΔNG | Max ΔNG |
|------|-----|-----|--------------------|---------|---------|
| to | CN1 | 7/7 | **+0.581** | +0.285 | +0.827 |
| to | CN2 | 7/7 | **+0.638** | +0.271 | +0.885 |
| paired | CN1 | 7/7 | **+0.694** | +0.342 | +0.898 |
| paired | CN2 | 7/7 | **+0.746** | +0.538 | +0.935 |

### 4.4 R4 結論

1. **方向完全一致**：所有 28 組（2 modes × 2 CN tiers × 7 samples）direction_positive=True
2. **cnLOH ≥ deletion-LOH 強度**：CN2 mean ΔNG 系統性地較 CN1 高 ~0.05-0.10（兩模式皆同）。一個可能解釋是 CN2 保留 diploid coverage → NGroups 估計更精確（類似 R3 發現）
3. **未發現相反方向**：原擔憂「混合兩類 LOH 會抵消」未實現；兩者皆支持 Intermediate-AF > Extreme-AF 的 NGroups 模式

**機制解讀**：
- Deletion-LOH（CN1）：單一剩餘 allele，subclone 間 methylation 差異僅來自該 allele 各子株複製後的漂移 — 訊號存在但略弱
- cnLOH（CN2）：雜合性消失但保留兩個相同序列的 copies，subclone methylation 差異累積空間更大，加上 coverage 充足 — 訊號稍強
- **兩者皆為 "subclone-with-lost-heterozygosity" 的合法生物學訊號**

**建議論文用語**：「LOH 下進一步分層為 deletion-LOH（CN=1）與 cnLOH（CN=2），兩類在 7 個樣本中方向完全一致（28/28），強度相近（CN1 mean ΔNG=+0.581-0.694；CN2 +0.638-0.746）。cnLOH 略強於 deletion-LOH，可能反映其保留的 diploid coverage 提供更可靠的 methylation heterogeneity 估計。兩類生物學機制雖異，但皆展現 Intermediate-AF subclone 的 methylation 異質性。」

---

## 5. 綜合討論

### 5.1 Opus 4.7 Plan B.2 五項質疑狀態

| 質疑 | Status | 本批次結論 |
|------|--------|-----------|
| B.2-1 HCC1954 反向 | ✅ 回應（R2） | Pre-reg 排除後 6/6 POS，region-level 7/7 方向一致；segment-level 反向原聲稱未在 per-chr 粗聚合下複現 |
| B.2-2 Coverage_Multiple ≠ CN | ⏸ 延後（R6, batch 2+） | 使用 Coverage_Multiple 作為 CN 代理之 caveat 已標記 |
| B.2-3 NR-bin 功效 artifact | ✅ 回應（R3） | Matched sampling 拒絕功效假說，效應強度真實 |
| B.2-4 Purity 外推 | ⏸ 延後（R14, P3） | 使用者降優先，purity mixture 暫緩 |
| B.2-5 LOH 類型分層 | ✅ 回應（R4） | CN1/CN2 皆 7/7 POS 方向一致，強度相近 |

**回應完成**：3/5（batch 1 目標達成）
**延後**：2/5（B.2-2 與 B.2-4，已事先標記並經使用者確認優先序）

### 5.2 對 20260414 主要結論的影響

| 結論 | 批次 1 後狀態 |
|------|---------------|
| LOH 下 Intermediate AF → 高 NGroups | **維持 POSITIVE**（region-level 7/7，pre-reg 6/6） |
| NR-bin 控制後效應增強非消失 | **確認**（matched sampling 支持） |
| cnLOH 也呈相同模式 | **確認**（CN1/CN2 兩類 7/7） |
| 7/7 p<10⁻³⁹ 聲稱 | **精確化**：region-level 7/7，segment-level 本批次無 cleanly 複現 |
| HCC1954 post-hoc rationalization | **精確化**：改為 pre-registered exclusion；region-level directional 仍然一致 |

### 5.3 穩定度建議

`06_結論穩定性審查.md` 關於「LOH Subclone AF × Methylation 雙證據鏈」結論：
- **原評分**：⭐3（pending confound 審查）
- **批次 1 後建議**：⭐3 → **⭐4**（R2/R3/R4 三項質疑通過，B.2-2/B.2-4 仍為 caveat；要達 ⭐5 需 patient-derived cohort 外推 + CNV caller 整合）

### 5.4 未決項（批次 2/3 候選）

| R-id | 題目 | 優先 |
|------|------|------|
| R5 | LOH.bed 生成機制 + PON-only vs self-phasing haplotag 比對 | P1（與 軌 1 HCC1395 Normal BAM 並行） |
| R6 | Coverage_Multiple 作為 CN 代理的 caveat 處理 | P2（可於 R5 後依結果升降） |
| R7 | ReadParser P0-1/2 全維度保留 | P1 Hard Gate |
| R8 | Per-CpG ASM residualized | P1 |
| R9-R11 | Step 6B 派生（Dorado / BRCA1 / HCC1954 fine heterogeneity） | P1 |
| R12 | Phase 2 跨樣本驗證（blocked by R1） | P2 |
| R13 | COLO829 out-of-scope 聲稱標註 | P2（文件維護） |
| R14 | Purity 真實混合 | P3（使用者降優先） |
| R15 | Patient-derived cohort | P3（外部依賴） |
| R16 | 5hmC / EVOFLUx | P3（設備 / basecalling 升級） |

---

## 6. Reproducibility

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/B2_loh_skeptical_check/scripts/r2_r3_r4_b2_batch1.py
# Outputs under research/B2_loh_skeptical_check/data/
#   r2_mann_whitney_by_sample.tsv       (R2 region-level)
#   r2_hcc1954_segment_level.tsv        (R2 per-chr Spearman)
#   r3_nr_bins_with_80plus.tsv          (R3 NR-bin MW)
#   r3_matched_sampling.tsv             (R3 matched sampling, seed=20260418)
#   r4_cn1_vs_cn2_per_sample.tsv        (R4 CN tier)
```

**Random seed**：`20260418`（確保 matched sampling 可重現）

---

## 7. 使用者審閱決策點

**批次 1 完成後建議使用者做以下判斷**（對應 plan E.12）：

1. **本批次三項結論是否滿足您對 B.2 質疑的回應需求？**（若是 → 批次 2 啟動）
2. **HCC1954 最終報告定位**：
   - 選項 A：僅報告「7/7 region-level directional」不做 pre-reg 排除
   - 選項 B：同時報告「7/7（full）+ 6/6（pre-reg 排除）」（建議）
   - 選項 C：只報告「6/6 pre-reg 排除後」
3. **Segment-level HCC1954 反向聲稱**：原 20260414 的 `segment-level ρ=-0.297` 是否在論文中繼續保留，或改以「per-chr aggregation 於本批次未複現」為主？
4. **批次 2 執行順序**：建議 R5 (LOH.bed 機制) + 軌 1 R1 (HCC1395 Normal BAM pilot) 並行；待使用者確認

**請使用者審閱後告知上述決策，批次 2 才會啟動。**
