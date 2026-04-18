<!--
建立時間: 2026-03-31 23:59
目標: 系統性觀察 O1-O10 Level A 結論整合報告
處理範圍: 748,391 rows × 116 columns, 7 samples × 2 modes, 9 個觀察主題
關聯檔案:
  - big7_disk_output/synthesis/observation_workspaces/OBSERVATION_INDEX.md
  - big7_disk_output/synthesis/observation_workspaces/20260401_cross_validation_report.md
-->

# InterSubMod 系統性觀察 Level A 結論整合報告

## 報告摘要

**日期**: 2026-03-31
**數據基準**: all_region_rows.tsv.gz (748,391 rows × 116 cols, post HP-fix)
**觀察範圍**: 9/10 主題完成 (O1-O8, O10), 82 張圖表, 53 份數據表
**交叉驗證**: 通過（零矛盾，數字一致）

---

## Top 10 Level A 發現

### 1. TO 模式下無任何單一特徵有效區分 TP/FP
**來源**: O1, O4, O5, O8（4 個觀察獨立確認）
**數字**: Caller 特徵 TO AUC 0.42-0.60; 甲基化特徵 TO AUC 0.43-0.54; per-sample TO AUC < 0.58
**影響**: TO 模式必須依賴多特徵組合或全新特徵工程

### 2. LOH Penalty 是 TO QS 失效的根本原因
**來源**: O2, O3
**數字**: LOH trigger TP=44.5% vs FP=35.8%（反向）; 移除後 AUC +0.045; TO QS ceiling ~0.55
**影響**: LOH penalty 必須在 TO 模式下移除或反轉

### 3. Paired 與 TO 是根本不同的問題空間
**來源**: O1, O5, O7
**數字**: Paired FP 1.04% vs TO 30.6%; 5/9 甲基化方向反轉; HP_Ratio cross-mode r=0.001
**影響**: 必須分離 paired/TO 模型

### 4. GQ 是 Paired 最強特徵（AUC=0.811），超越 QS Composite
**來源**: O4, O8
**數字**: Per-sample 穩定 0.755-0.947; Cohen's d=1.314; QS paired AUC=0.754
**影響**: Paired ML baseline 以 GQ 為首要特徵

### 5. 樣本間異質性極大
**來源**: O8
**數字**: TO FP rate 8.7% (H2009) to 74.6% (HCC1954); H2009 佔 36.2%
**影響**: 需要 sample-aware 策略或 LOSO-CV 評估

### 6. ISM 甲基化特徵 Region-Level 鑑別力極弱
**來源**: O5, O10
**數字**: 最佳 AUC=0.543; read-level TP/FP AUC=0.737（但受 region clustering 膨脹）
**影響**: 需要新特徵工程方向

### 7. VerificationClass 無法作為 TP/FP 過濾器
**來源**: O6
**數字**: Cramér's V=0.023-0.024; paired-TO kappa=0.854
**影響**: 應以連續的 AlleleDelta/HPMergedDelta 取代

### 8. AF 在 TO 反向（高 AF = 更多 FP）
**來源**: O4, O8
**數字**: TO AF AUC=0.418; FP rate @ AF 0.8-0.9 = 55.2%
**影響**: AF 硬閾值在 TO 有害

### 9. TO 過度判定 LOH（85.5% concordance）
**來源**: O3, O7
**數字**: 95.5% discordant = TO=LOH where paired=nonLOH; partial genotype 驅動
**影響**: TO LOH 判定不可信，effective_hp_reads 是更好的連續替代

### 10. Read-Level 甲基化對 ALT/REF 分類無用
**來源**: O10
**數字**: 12 個特徵 AUC 0.500-0.547; 跨樣本 delta 方向不一致
**影響**: 不應投入 read-level 甲基化 ALT/REF 分類器

---

## 行動建議

| 優先級 | 行動 | 依據 | 預期效果 |
|--------|------|------|---------|
| **P0** | 移除 TO LOH penalty | O2, O3 | QS AUC +0.045 |
| **P0** | 建立 Paired/TO 分離模型策略 | O1, O5, O7 | 避免方向反轉特徵汙染 |
| **P1** | Phase 1A ML 特徵集: GQ + DP + 5 甲基化 + effective_hp_reads | O4, O5, O3 | Paired baseline AUC ~0.85+ |
| **P1** | 移除 VerificationClass 從 QS 決策 | O6 | 降低噪音 |
| **P2** | Sample-aware calibration 或 LOSO-CV | O8 | 處理 8.6× FP 率差異 |
| **P2** | 執行 O9 FN 觀察 | 待 FN ISM | 量化 LOH rescue 潛力 |

---

## 數據來源

**Workspace 路徑**: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/`
**索引**: `OBSERVATION_INDEX.md`
**交叉驗證**: `20260401_cross_validation_report.md`

所有結論均經交叉驗證通過（零矛盾，AUC 數字精確吻合）。
