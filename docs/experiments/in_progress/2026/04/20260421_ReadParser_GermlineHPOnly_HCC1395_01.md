---
status: validated_conditional_negative
date: 2026-04-21
type: phase1_verification
phase: 1
plan: /bip7_disk/liaoyoyo2001/.claude/plans/streamed-spinning-wilkinson.md
hypothesis_id: H-HP-ONLY-01
priority: P0
pipeline_track: TO
outcome: filter_negative_characterization_positive
---

# ReadParser Germline-HP-Only — Phase 1 HCC1395 全 TO 驗證

## 結論

**修正機制正確**（Phase 0 已確認），但 **Phase 1 Gate 未通過**：無任何 ISM TSV 特徵因修正而獲得 AUC ≥0.02 的改善，部分 HP 衍生特徵 AUC 反而輕微下降。因此：

- **不進 Phase 2**（7 樣本 × 2 模式全量重跑）
- **修正保留**（default=off，不影響既有流程；flag=on 提供研究者獲得乾淨 HP 分群的工具）
- **價值定位轉向**：filter 方向 NEGATIVE；characterization 方向保留（subclone 結構解析不再被 self-phasing 污染）
- **plan R3 條件未觸發**：HP_Ratio 未歸零（停留 0.529），bug 不在上游 assignment 層，不需升級為 LongPhase-TO C++ 修復

---

## 測試設定

| 參數 | 值 |
|------|----|
| BAM | `pononly_v3_fixed/tumor_tagged.bam`（V3-Fixed haplotag TO BAM，287 GB）|
| Normal BAM | HCC1395 normal.bam |
| VCF (TP) | `clairsto_v3fixed_work/clairsto_tp.vcf.gz` — 28,509 sites |
| VCF (FP) | `clairsto_v3fixed_work/clairsto_fp.vcf.gz` — 11,606 sites |
| Reference | GRCh38_no_alt_analysis_set.fasta |
| Window | 5,000 bp |
| Metric | BERNOULLI |
| Threads | 24 × 2 parallel |
| 4 runs | `{tp,fp} × {off,on}` = 4 組 |

執行時間：TP off 21 min、TP on 19 min、FP off 11 min、FP on ~11 min。
Memory 穩定 <28 MB per run。

---

## 驗證 1：Audit 機制一致性（Phase 0 延伸）

| 欄位 | TP sum | FP sum | off/on identical |
|------|--------|--------|-----------------|
| NHP_Somatic11 | 347,213 | 122,307 | ✅ True (TP+FP) |
| NHP_Somatic21 | 308,762 | 149,308 | ✅ True (TP+FP) |
| NHP_Somatic33 | 124,096 | 50,081 | ✅ True (TP+FP) |

全基因 scale audit 獨立性驗證通過。

**Somatic tag 分佈密度**（per site 平均）：
- TP: 11 → 12.2/site，21 → 10.8/site，33 → 4.4/site（合計 27.4/site）
- FP: 11 → 10.5/site，21 → 12.9/site，33 → 4.3/site（合計 27.7/site）

TP/FP somatic tag 密度幾乎相同 → 單憑 somatic tag 本身對 TP/FP 區分無訊號。

---

## 驗證 2：HP_Ratio 全基因 median 比對

### Plan 預期 vs 實測

| 配對 | Plan 預期 baseline | Plan 預期 on | 實測 off | 實測 on | delta |
|------|-------------------|-------------|---------|--------|-------|
| TP HP_Ratio median | **0.836** | 0.55-0.65 | **0.549** | 0.529 | -0.020 |
| FP HP_Ratio median | — | — | 0.516 | 0.500 | -0.016 |

**關鍵發現**：Plan 的 baseline 0.836 與實測 0.549 不符。差異來源為：
- Plan 引用來源（landscape doc 02_Self_Phasing根因.md）可能使用不同 VCF 子集或舊版 haplotag BAM
- V3-Fixed TO BAM 的 HP_Ratio 本就較低 → plan 「修正後 0.55-0.65」的假設已包含 baseline

**分層檢視（TP only）**：

| 分層 | n | off median | on median |
|------|---|-----------|-----------|
| Potential_LOH=True | 2,070 | 0.0909 | 0.0927 |
| Potential_LOH=False | 26,439 | 0.5542 | 0.5312 |

- LOH regions：幾乎不動（已極度不平衡，somatic tags 本就少參與）
- Non-LOH regions：-0.023 shift（方向正確，但幅度遠小於 plan 預期）

**plan R3 條件**：「修正後 HP_Ratio 仍偏 0 → bug 不在 hp_tag 解析」— 實測 HP_Ratio 停留 0.53，**未觸發此條件**；修正定位為「解析正確」，無需升級為 LongPhase-TO 上游 C++ 修復。

---

## 驗證 3：AUC Gate（TP vs FP 分類能力）

**每個特徵取「最佳方向 AUC」（max(auc, 1-auc)）**：

| 特徵 | AUC_off | AUC_on | delta | 判定 |
|------|---------|--------|-------|------|
| HP_Ratio | 0.5257 | 0.5137 | -0.0119 | ⬇ |
| HPFineNGroups | 0.5359 | 0.5099 | -0.0260 | ⬇⬇ |
| HeuristicScore | 0.5066 | 0.5088 | +0.0022 | — |
| **Quality_Score** | 0.5251 | 0.5148 | -0.0103 | ⬇ plan 目標 ≥0.55 **FAIL** |
| Coverage_Multiple | 0.5126 | 0.5126 | ±0.0000 | — |
| CramersV | 0.5004 | 0.5010 | +0.0007 | — |
| CramersV_HPFine | 0.5019 | 0.5012 | -0.0007 | — |
| HPFineF | 0.5117 | 0.5045 | -0.0071 | — |
| HP1FamilyN | 0.5086 | 0.5171 | +0.0084 | — |
| HP2FamilyN | 0.5420 | 0.5295 | -0.0125 | ⬇ |
| NHP3 | 0.5350 | 0.5000 | -0.0350 | ⬇⬇ (constant 0 on on-run) |
| NHP0 | 0.5563 | 0.5528 | -0.0036 | — |
| NumReads | 0.5085 | 0.5085 | ±0.0000 | — |
| Fisher_Frac_Sig | 0.5047 | 0.5127 | +0.0080 | — |
| Stability | 0.5000 | 0.5000 | ±0.0000 | — |
| HPMergedDelta | 0.5406 | 0.5154 | -0.0252 | ⬇⬇ |
| AlleleDelta | 0.6294 | 0.6294 | ±0.0000 | — (HP-independent) |
| HPFine_NGroups_CF | 0.5359 | 0.5099 | -0.0260 | ⬇⬇ |

**關鍵觀察**：
- 無任何特徵 delta AUC ≥ +0.02
- 4 個 HP 衍生特徵（HPFineNGroups、NHP3、HPMergedDelta、HPFine_NGroups_CF）delta AUC ≤ -0.025
- AlleleDelta（唯一 HP-independent 中等訊號特徵，AUC=0.629）完全不動 — 與預期一致

---

## 驗證 4：`within_dom_alt_frac` 特殊說明（plan 首要目標）

**無法在本次 Phase 1 直接驗證**。理由：
- 該特徵（LOSO AUC 0.721，plan landscape doc 03_ISM分析價值界定.md line 109）為 **downstream Python 派生特徵**，非 ISM TSV 直接欄位
- 本專案 scripts/ 與 research/ 目錄 grep 無此特徵計算腳本（可能原分析腳本未提交或已搬移）
- 重建該特徵需從 40,115 個 per-region `reads.tsv` 重新聚合，成本高於 Phase 1 minimum scope

**間接推論**：
- 該特徵 = Dominant HP reads 中 ALT 比例
- 修正後 Dominant HP 定義改變（somatic tags 不再參與 HP label）
- 若 plan 預期「AUC 0.679 → 0.70」成立，應在 TSV 層級看到 HP 衍生特徵同向改善
- 實測 HP 衍生特徵 **全體不改善或下降** → within_dom_alt_frac 改善至 ≥0.70 的後驗機率低
- 嚴格驗證須 Phase 1.5（重建 downstream pipeline + 重算特徵），視研究必要性另行決策

---

## 驗證 5：HPFineNGroups 分佈收斂

| NGroups | TP off | TP on | FP off | FP on |
|---------|--------|-------|--------|-------|
| 0 | 11 | 38 | 1 | 8 |
| 1 | 140 | 742 | 34 | 80 |
| 2 | 5,626 | **27,729** | 3,423 | **11,518** |
| 3 | 18,729 | 0 | 6,385 | 0 |
| 4 | 4,003 | 0 | 1,763 | 0 |

- 97.3% TP regions on-run 歸在 NGroups=2
- 82.0% TP regions 的 HPFineNGroups 值因 flag 改變
- 修正前：NGroups ≥3 被當作 subclone 候選訊號（占 80% TP）
- 修正後：僅 germline HP1/HP2 +unphased → 三類 labels，NGroups 極難超過 2

**意涵**：
- 「HPFineNGroups≥4 標記 subclone marker」的結論（見 memory `project_hpfinengroups_subclone_marker.md`）**實質依賴 somatic HP tags** 的 self-phasing 分類
- 修正後該 marker 失去辨識力 → 需重新評估此 marker 的生物學意義

---

## Phase 1 Gate 判定表

| 檢核項 | Plan 預期 | 實測 | 判定 |
|-------|----------|------|------|
| HP_Ratio TP 相關 (r) | 0.001 → >0.2 | 未計算（需 paired/TO 兩套 haplotag 對比）| — |
| within_dom_alt_frac AUC | 0.679 → ≥0.70 | 無法本地驗證（downstream 派生） | INCONCLUSIVE |
| TO Quality_Score AUC | 0.497 → ≥0.55 | 0.525 → 0.515 (delta -0.010) | ❌ FAIL |
| 任一特徵 AUC delta ≥ 0.02 | 至少 1 項 | 0 項 | ❌ FAIL |
| R3: HP_Ratio 歸零? | 若發生 → 升級 upstream 修復 | median 0.529 ≠ 0 | ✅ 無需升級 |

**Gate 總判定**：**FAIL — 不進 Phase 2**

---

## 結論與後續策略

### 本次修正的定位

1. ✅ **機制正確**：Config / CLI / Parser / TSV audit 欄位全通；unit tests 12/12 pass；demotion 數學守恆
2. ✅ **研究誠實性提升**：移除 LongPhase-TO self-phasing 循環依賴；HP 分群現在僅反映 germline 結構
3. ❌ **Filter 價值 NEGATIVE**：無 TSV 特徵 AUC 顯著改善；甚至多個 HP 衍生特徵 AUC 下降
4. ⚠ **Characterization 價值條件性**：subclone marker (HPFineNGroups≥4) 結論需重新審視 — 可能原先的「2.1% 的 TP 在 NGroups≥4」其實是 somatic tag 污染產生的人工分組，**不代表真實 subclone**

### 本次釐清的研究洞見（比 AUC 更重要）

**HPFineNGroups≥4 的 subclone marker 結論需重新評估**：

- Memory 記錄 `project_hpfinengroups_subclone_marker.md`：「N≥4+NR≥80 TP rate 89.1%」
- 但本 Phase 1 顯示：HPFineNGroups≥3 在 flag=on 時**完全消失**（0 regions）
- 解讀：原先 NGroups≥4 的訊號 **可能來自 somatic HP tags 的細分**（11/21/33 三個額外 labels），而非真實 subclone 結構
- 建議：重新跑「NGroups≥4 = TP marker」在 flag=on 下的行為；若該 marker 完全失效，需撤回 POSITIVE 結論

### 不進 Phase 2 的理由

1. 4 runs × 7 samples × 2 modes = 56 runs 的 bench cost ~3-5 day 機器時
2. Phase 1 已顯示無 AUC 改善潛力 → Phase 2 擴大 sample pool 不會翻盤
3. CovM baseline blocker 與本項非互相支援關係（本修正不依賴正確 CovM）

### 建議的後續研究方向（優先序）

| 優先 | 項目 | 成本 | 預期產出 |
|------|------|------|---------|
| P1 | 重新審視 HPFineNGroups subclone marker（flag=on 下重算 TP rate）| 1-2 hr | 可能撤回 POSITIVE 結論 |
| P2 | within_dom_alt_frac 重建 downstream pipeline（若研究方向必要）| 1 day | 特徵 AUC 最終判定 |
| P3 | 將此修正推廣為 ISM 預設行為（若確認無 filter 傷害）| 思考 | 未來版本 default=on 議題 |

---

## 檔案清單

- **驗證資料**：`/tmp/ism_hp_fix_phase1/{tp,fp}_{off,on}/significance_summary.csv`
- **合併資料**：`/tmp/ism_hp_fix_phase1/merged_{off,on}.csv`
- **Phase 0 報告**：`docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_smoke_01.md`
- **本報告**：`docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md`
- **Commit**：Phase 0 changes at `775027d` on `refactor/phase1-safety`

---

## 對 memory 的影響

**需更新**：
- `project_hpfinengroups_subclone_marker.md` — 加入「flag=on 下 NGroups≥4 訊號消失」警告
- `project_self_phasing_causal_chain_confirmed.md` — 實測 HP_Ratio median 0.549（而非 landscape 引的 0.836）

**需新增**：
- 新 memory：`project_readparser_germline_hp_only_phase1_negative.md`（filter 負面結論）
