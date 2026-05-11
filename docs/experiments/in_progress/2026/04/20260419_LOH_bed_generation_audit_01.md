---
title: R5 — LOH.bed 生成機制審查與 self-phasing 污染檢驗
date: 2026-04-19
status: VALIDATED
owner: InterSubMod Research
scope: Batch 2 軌 1 · LOH.bed 可信度確認
related:
  - 20260418_B2_LOH_Subclone_Deep_Skeptical_Check_01.md
  - 20260419_Z3_amplicon_blacklist_pilot_result_01.md
  - research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md
  - research/loh_investigation/scripts/observe_to_ge_vs_seqc2.py
r_id: R5
priority: P0
---

# R5 — LOH.bed 生成機制審查

## 一、執行摘要

- **結論**：LOH.bed **不受** self-phasing 循環依賴污染，可作為下游（B.2 批次 1、Z3 pilot、F pilot NonLOH filter）的**可信前置輸入**
- **關鍵證據**：PON-only haplotag vs self-phasing haplotag 兩個版本產生的 LOH.bed **Jaccard=1.0000**（1,094 regions / 1,632.2 Mb 完全一致）
- **外部驗證**：LOH.bed vs SEQC2 truth set F1=96.2%、Jaccard=0.8470（修改前後無變化）
- **機制根因**：LongPhase-TO 使用 **phased VCF 的 region-level genotype ratio** 判定 LOH（`--loh` flag），**不**依賴 tumor BAM 的 haplotag 狀態
- **對 batch 1 / F pilot / Z3 pilot 的影響**：全部依賴 LOH.bed 的結論免除此條 confound；Opus 4.7 plan 結論 9 隱含假設 5（LOH.bed 不受 self-phasing 汙染）✅ **驗證通過**

---

## 二、審查問題定義

### 2.1 疑慮來源

Opus 4.7 plan E.9 R5 列出：

> 若 LOH.bed 生成**也**依賴 haplotag（self-phasing 下游），LOH.bed 本身可能汙染；F pilot 全部 NonLOH filter 依賴 LOH.bed 正確性。

本專案 2025 年底已確認 ISM HP_Ratio 有 62% self-phasing artifact（結論 9）。若 LOH.bed 的生成流程也「看」了同一組 tumor haplotag → 則 B.2、Z3、F pilot 全部的 NonLOH 過濾條件可能內含相同循環。

### 2.2 驗證標準（E.9 原文）

> bed diff <5% → LOH.bed 可信；>5% → 擴展至所有 LOH-dependent 結論重審

---

## 三、機制追溯（原始資料 → 結論）

### 3.1 LOH.bed 生成流程

依據既有 `research/loh_investigation/` 報告（2026-04-03）與 C++ 程式碼審查：

```
tumor BAM（含 MM/ML + 可選 HP tag）
        │
        ↓
ClairS-TO 呼叫變異 → tumor.vcf（含 caller_af, DP, AD）
        │
        ↓
LongPhase-TO 相位化（phased VCF）
        │
        ↓  (--loh flag)
LongPhase-TO 以 phased VCF 的 region-level GT ratio 判斷 LOH
        │  （不讀 BAM haplotag；純 VCF-level logic）
        ↓
tumor_phased_GE.bed — 輸出 LGE / SGE / LOH 區段
        │
        ↓  (ISM 下游抽取)
to_loh_bed_hit (per-region boolean)  →  F pilot / Z3 pilot / B.2 使用
```

### 3.2 C++ 程式碼佐證

| 檔案 | 行號 | 關鍵邏輯 |
|------|------|---------|
| `include/core/LohBedAnnotator.hpp` | L10-71 `load(bed_path)` | 讀取 BED 檔、索引區間 |
| `include/core/LohBedAnnotator.hpp` | L73-77 `overlaps(chrom, pos)` | 純 coordinate lookup |
| `include/core/LohBedAnnotator.hpp` | L88-95 `classify(chrom, pos, hp_ratio_loh)` | BED-based，不依賴 haplotag 判斷 |
| `src/core/RegionProcessor.cpp` | L456-461 | LOH BED loading 路徑 |
| `src/core/RegionProcessor.cpp` | L1664-1672 | LOH annotation 輸出 |

**關鍵觀察**：`LohBedAnnotator::classify()` 雖接受 `hp_ratio_loh` 參數，但這只用於 **subtype 分類**（cnLOH vs deletion-LOH），**不參與** "是否為 LOH" 本身的判定。"是否為 LOH" 完全由 BED 區間 overlap 決定，而 BED 由 LongPhase-TO 上游產出。

---

## 四、實證檢驗（Jaccard=1.0000）

### 4.1 實驗設計

參見 `research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md`：

- **Arm A**：self-phasing haplotag（tumor BAM → LongPhase → phased VCF → LOH.bed）
- **Arm B**：PON-only haplotag（tumor BAM + PON reference → LongPhase → phased VCF → LOH.bed）

兩者差別：Arm A 讓 tumor 自己參與 phasing（循環風險）；Arm B 只用外部 PON reference（無循環）。

### 4.2 結果

| 指標 | Self-phasing | PON-only | 差異 |
|------|------------|----------|------|
| LOH region 數 | 1,094 | 1,094 | 0 |
| LOH 總長度 | 1,632.2 Mb | 1,632.2 Mb | 0 |
| Jaccard（兩 arm 間） | — | — | **1.0000** |
| Jaccard vs SEQC2 | 0.8470 | 0.8470 | 0 |
| F1 vs SEQC2 | 96.2% | 96.2% | 0 |

**結論**：LOH.bed 在 haplotag 策略切換下**完全一致**，遠優於 E.9 設定的「<5% diff 可信」標準。

### 4.3 為何 LOH.bed 不受 self-phasing 影響？

LongPhase-TO 的 `--loh` 判定邏輯是 **region-level phased VCF GT ratio**：

- 對 phased block 內所有變異計算 `het_ratio = n_het / n_total`
- 若 `het_ratio` 顯著偏離 0.5（例如 <0.2 或 >0.8）→ 該 block 標為 LOH
- 此邏輯**只讀 VCF 的 GT 欄位**（ClairS-TO 輸出），**不**讀 BAM 的 HP tag

換言之，self-phasing 汙染的是 **read-level HP assignment**（B.2 / Z3 FP 的 HP_Ratio 偏移），而非 **region-level LOH 判定**。兩者共用同一 phased VCF 作為 parent，但下游分岔後互不干擾。

---

## 五、對依賴 LOH.bed 結論的影響

| 結論 / 聲稱 | 依賴 LOH.bed 的哪個欄位 | R5 審查後狀態 |
|------------|----------------------|-------------|
| F pilot `NG=4+AF<0.4+NR≥80+NonLOH` filter | `to_loh_bed_hit=False` | ✅ 可信 |
| B.2 批次 1 R2/R3/R4 的 LOH 定義 | `tool_potential_loh`, `core_loh_like`（paired）+ `to_loh_bed_hit`（TO） | ✅ 可信 |
| Z3 pilot Z3 Zone 定義 | `to_loh_bed_hit AND AF extreme AND HPFineNGroups<=1` | ✅ 可信 |
| Z3 amplicon blacklist pilot S1/S2 | Z3 ∩ chr5/8/17 arm | ✅ 可信 |
| Potential_LOH column（結論 9）| HP_Ratio<0.1 OR >0.9（ISM 內部，**非 LOH.bed**）| ⚠ 仍是 self-phasing downstream（既有結論 9 早已鎖定這部分受汙染，本 R5 不改變） |

**關鍵區分**：
- `to_loh_bed_hit` / `tool_potential_loh` / `core_loh_like` → 從 LongPhase-TO LOH.bed 派生 → **R5 驗證可信**
- `Potential_LOH`（ISM internal HP_Ratio binary）→ 從 tumor BAM HP tag ratio 派生 → **既有結論 9 已鎖定為 self-phasing artifact**（本 R5 不改變此狀態）

---

## 六、殘留風險與限制

1. **單樣本驗證**：2026-04-03 的 Jaccard=1.0000 實驗在 **HCC1395**（SEQC2）做。其他 6 樣本（COLO829, H2009, H1437, HCC1937, HCC1954, HCC1395_DORADO）未逐一重做
   - 緩解：LongPhase-TO 的 LOH 邏輯不依賴樣本特性（純 VCF GT ratio），跨樣本一致性有機制保證
   - 升級條件：若未來任何樣本發現 LOH.bed 內部不一致，觸發 sample-specific re-audit
2. **LOH.bed 錯誤非 self-phasing 類**：LOH.bed 本身可能因 phasing block 斷裂、低 AF 變異漏偵、或 PON 品質問題而有**絕對錯誤**（vs SEQC2 F1=96.2% 表示 ~4% 錯誤存在）
   - 定位：此為 **LOH caller accuracy**，非 self-phasing 循環
   - 影響：4% 錯誤率進入下游時屬 irreducible noise，不是 R5 要解決的問題
3. **paired mode 的 `core_loh_like`**：此欄位來自 ClairS Paired 模式 + SEQC2 matched normal，與 TO mode 路徑正交，R5 聚焦 TO 路徑
   - paired mode LOH 判定由 ClairS Paired 自帶（非 LongPhase-TO），本 R5 不涵蓋
   - 若有疑慮需另開 R5b（目前無證據指出 paired LOH 有循環依賴）

---

## 七、產物清單

### 新增
- 本報告 `docs/experiments/in_progress/2026/04/20260419_LOH_bed_generation_audit_01.md`

### 引用既有
- `research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md`（Jaccard=1.0000 原始證據）
- `research/loh_investigation/reports/20260403_pon_only_haplotag_ism_verification_report.md`（ISM 下游對 LOH.bed 一致性的影響）
- `research/loh_investigation/scripts/observe_to_ge_vs_seqc2.py`（LongPhase-TO GE.bed 分析工具）
- `research/loh_investigation/reports/20260402_seqc2_vs_longphase_to_loh_validation.md`（LOH.bed vs SEQC2 F1=96.2%）

---

## 八、結論與後續行動

### 8.1 結論

- R5 **PASS**（Jaccard=1.0000 遠超 E.9 的 <5% diff 標準）
- LOH.bed 可作為**下游所有依賴 LOH 定義結論的可信前置輸入**
- 結論 9 隱含假設 5（LOH.bed 不受 self-phasing 汙染）✅ 驗證通過

### 8.2 後續行動

- **不需新增實驗**：既有 2026-04-03 證據已充分
- **穩定度表更新**：`06_結論穩定性審查.md` 結論 9（62% ISM HP_Ratio LOH artifact）保持精確化聲稱「**非 LOH.bed region-level LOH**」正確，無需修改
- **批次 2 並行軌不受影響**：R6（CovM z-score）與 R1（HCC1395 Normal BAM pilot）的 LOH filter 條件 continue as planned
- **若有樣本-特異 LOH.bed 異常浮現**（例如 HCC1954 Z3 pilot 之外的 FP pattern）→ 新開 R5-sample (sample-specific re-audit)
