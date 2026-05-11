---
title: Batch 2 三軌同步檢核與對齊分析
date: 2026-04-19
status: PARTIAL COMPLETE（R5 + R6 完成；R1 in progress）
owner: InterSubMod Research
scope: Opus 4.7 plan G.3 批次 2 · 三軌交叉驗證
related:
  - 20260419_LOH_bed_generation_audit_01.md
  - 20260419_Coverage_Multiple_zscore_normalization_01.md
  - ../../../research/F_hpfinengroups_deepening/observations/step7_hcc1395_normal_TO_pilot.md
  - ../../../research/F_hpfinengroups_deepening/observations/step7_hcc1395_normal_paired_pilot.md
---

# Batch 2 三軌同步檢核與對齊分析

> 使用者指示：「三軌並行全啟動,同步檢核與較驗對齊與分析」

## 一、三軌狀態表

| 軌 | R-id | Scope | 狀態 | 結論 | 對其他軌影響 |
|----|------|-------|------|------|-------------|
| 1 | R5 | LOH.bed 可信度審查 | ✅ **COMPLETED** | PASS（Jaccard=1.0000 >> 5% 標準）| 解鎖 R1/R6 LOH filter 使用 |
| 2 | R1 | HCC1395 Normal BAM pilot | 🔄 **IN PROGRESS** | BED+subset 就緒；Phase 2 pipeline 待執行 | 依賴 R5（LOH filter）+ R6（CovM caveat）|
| 3 | R6 | CovM z-score normalization | ✅ **COMPLETED** | CONDITIONAL NEGATIVE（CovM 非獨立 CN proxy） | **強化 R1 必要性**（Normal BAM Δ_coverage 成為唯一 sample-matched CN signal） |

## 二、Cross-track 對齊驗證

### 2.1 R5 → R1 & R6：LOH filter 合法性

**R5 發現**：LOH.bed 在 PON-only vs self-phasing 兩版本間 Jaccard=1.0000（1,094 regions 完全一致），F1 vs SEQC2=96.2%。LongPhase-TO 使用 phased VCF 的 region-level GT ratio 判定 LOH（不讀 BAM haplotag），與 ISM HP_Ratio self-phasing 循環分岔。

**對 R1 的影響**：
- R1 filter BED 的 NonLOH 條件（`~to_loh_bed_hit` for TO；`~tool_potential_loh & ~core_loh_like` for paired）**使用合法**
- HCC1395 1,843 filtered regions 的 LOH 判定免除循環 confound

**對 R6 的影響**：
- R6 的 LOH region 定義（Z3 inside LOH）直接使用 `to_loh_bed_hit`，**使用合法**
- HCC1954 Z3 2,030 regions 的 LOH 判定受 R5 保護

### 2.2 R6 → R1：CovM 失效強化 Normal BAM 必要性

**R6 發現**：Per-sample non-LOH baseline z-score normalize 後：
- HCC1954 Z3 內 FP z_extreme（z>2）rate 僅 0.15%
- 原 Step 3「HCC1954 Z3 FP CovM=0.733 偏高」在全局尺度上**不是**極端值
- 意味 raw CovM ≈ (per-sample overall coverage) × (local CN)；per-sample baseline 吸收大部分訊號

**對 R1 的影響（強化）**：
- R1 Normal BAM 提供 **sample-matched coverage baseline** — Δ_coverage(tumor/normal) 可直接分離 sample overall 與 local CN
- 若 R1 成功 → Δ_coverage 可作為 CN-aware 特徵（R6 否證的 CovM 單獨使用 → R1 的 Normal-anchored 替代）
- 若 R1 失敗 → 意味 ASM 本身（not coverage）也無獨立 CN signal，需 R17 CNV caller 整合

### 2.3 R5 ↔ R6：正交驗證

- R5 確認 LOH.bed region-level 可信（non-haplotag 循環）
- R6 確認 CovM 在 per-sample baseline 下 dissolve → 不是 LOH.bed 的 haplotag 問題，而是 coverage scaling 本質
- 兩 findings **不衝突**：R5 保證 region 選擇正確，R6 指出選定 region 的 CovM 訊號無 CN 特異性

### 2.4 三軌整合對 Opus 4.7 plan 假設的影響

| Opus 4.7 plan 假設 | R5/R6/R1 驗證 | 狀態 |
|-------------------|---------------|------|
| 假設 4「Coverage_Multiple ≈ CN 代理」 | R6 ❌ 否證（per-sample normalize 後 dissolve） | **FALSIFIED** |
| 假設 5「LOH.bed 不受 self-phasing 汙染」 | R5 ✅ 驗證（Jaccard=1.0000） | **VALIDATED** |
| 假設 7「HPFineNGroups 非 NR-binned artifact」| F pilot Step 3 + 5 | 已驗證（independent of R5/R6/R1） |
| 假設 8「LOH 單一類型（cnLOH vs deletion）」 | R6 揭示 CovM 不足以分型 → 需 R17 | **PENDING（升級 R17）** |

## 三、對既有結論的穩定度修訂建議

| 結論 | 原穩定度 | 批次 2 後建議 |
|------|---------|-------------|
| 9. 62% ISM HP_Ratio LOH artifact（**非 LOH.bed**） | ⭐3 | ⭐**4** — R5 獨立確認 LOH.bed 不受 self-phasing，原精確化聲稱正確 |
| 11. LOH Subclone AF×Methylation 雙證據鏈 | ⭐3 | ⭐**4** — B.2 批次 1 已升級；R5/R6 不翻轉 |
| 新增 17. HCC1954 CNV-driven reversal 機制 | — | ⭐**4** — B.2 + Z3 pilot 雙獨立驗證（F.6 表）|
| 新增 18. Z3 amplicon blacklist CONDITIONAL-NEGATIVE-for-canonical | — | ⭐**5** |
| 新增 19. Coverage_Multiple 非獨立 CN proxy | — | ⭐**4**（R6 per-sample normalize 驗證）|

## 四、R1 待完成項目 & 對批次 2 收尾影響

### 4.1 R1 剩餘工作

- [ ] samtools extraction 完成（背景執行中；~164MB / 預期 1-3GB）
- [ ] 確認 HCC1395 TO + paired tumor BAM + VCF canonical paths（可能需查 `output/canonical/`）
- [ ] 執行 Phase 2 pipeline TO arm + paired arm（~2h × 2，可並行）
- [ ] 撰寫 Δ_ASM 計算腳本
- [ ] 更新 `step7_hcc1395_normal_TO_pilot.md` 與 `step7_hcc1395_normal_paired_pilot.md` 至 VALIDATED/CONDITIONAL/NEGATIVE

### 4.2 若 R1 需延至下次會話

**不 block** 批次 2 其他成果彙報，因：
- R5（完成）+ R6（完成）已獨立回答 Opus 4.7 plan 假設 4/5
- R1 本身即使未完成 BED+BAM subset 仍可支援下一階段（Phase 2 pipeline 可在下會話啟動）
- 穩定度表更新項目（9, 11, 17, 18, 19）均不依賴 R1 結果

## 五、下批次（Batch 3）候選優先序（per Opus 4.7 plan G.6）

| 候選 R-id | Scope | 優先序 | 依賴 |
|----------|-------|-------|------|
| R1 完成後延伸 | Phase 2 TO + paired pipeline execution | P0 | samtools extraction |
| R8 | Per-CpG ASM fisher_frac_sig residualized（task #19）| P1 | 無 |
| R12 | Phase 2 跨樣本 Sample ASM 擴展至其他 5 in-scope | P1 | R1 成功 |
| R17（新）| HCC1954 Normal BAM pilot 延伸 | P1 | R1 成功 + 使用者授權 HCC1954 normal BAM 複製 |
| R9/R10/R11 | F pilot Step 6B 派生（Dorado/BRCA1/HCC1954 二階）| P2 | 無 |
| Research Chain Registry 新交付物 | Q7 meta-methodological（task #35）| P1 | 批次 2 完成後統整 |

## 六、彙報建議（給使用者）

1. **R5 PASS**：LOH.bed 不受 self-phasing 污染，F pilot / B.2 / Z3 pilot 的 LOH filter 合法性獲得獨立證據支持
2. **R6 CONDITIONAL NEGATIVE**：Coverage_Multiple 非獨立 CN proxy；強化 Normal BAM pilot 必要性；Opus 4.7 plan 假設 4 ❌ 否證
3. **R1 IN PROGRESS**：Filter BED 1,177 regions (12.6 Mb) ready；HCC1395 Normal BAM subset extraction 背景執行中；Phase 2 pipeline 需下一步啟動
4. **批次 2 整體**：兩軌完成，一軌進行中；對既有結論穩定度建議 3 項升級（9, 11）+ 3 項新增（17, 18, 19）
5. **下一步決策點**：
   - (a) 是否先進入 Batch 3 候選項目（不 blocking R1）
   - (b) 或等 R1 Phase 2 pipeline 執行完成（~4h）後再一次彙報
   - (c) Research Chain Registry（task #35）是否現在啟動（可與 R1 Phase 2 並行）
