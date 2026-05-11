<!--
build_date: 2026-04-29 00:55
revised: 2026-04-29 03:15 (fact-check 全面校正 — 詳見附錄 B)
agent: structured-tech-report skill (首個示範案例 + 反例教材)
status: validated (post fact-check)
report_class: methodology-improvement
audience: PI + engineer
inputs:
  - InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/04_imbalance_ratio_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md
  - InterSubMod/docs/reports/pi_reports/2026/04/20260424_V5_vs_Baseline_complete_comparison_01.md  # PI 報告 4 全基因組 V5=90.5%
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp  # V5 patch 實際位置（非 InterSubMod repo）
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.h
  - /big7_disk/liaoyoyo2001/longphase-to-mod/Util.h
  - InterSubMod/scripts/analysis/v5_imbalance_improvement.py
  - InterSubMod/scripts/analysis/v5_sanity_paired_check.py
outputs:
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
related_memory:
  - memory/project_v5_somatic_fallback_verification.md
  - memory/project_v3_fixed_haplotag_verification.md
  - memory/project_self_phasing_causal_chain_confirmed.md
  - memory/project_pon_only_phasing_verification.md
verdict: VALIDATED — V5 為 production；Layer 1.5 增強 + SNP alt guard 仍在 working tree (建議切 2 commit)；7 樣本擴展 / cnLOH 區仍為 open issue
last_verified: 2026-04-29 (post fact-check vs v5_audit_suite 6 份子報告)
report_template: structured-tech-report v1.0 (13-section)
-->

# longphase-to-mod 4-Commit 演進與 V5 Layer 1.5 Somatic Fallback：Self-Phasing 因果鏈整理

> **TL;DR**
> - **問題**：原始 LongPhase-TO 在 ClairS-TO TO 模式下出現 self-phasing — `--pon-only-phasing=false` 把 somatic variants 當 phasing anchor（PhasingProcess.cpp:154-157 的 `convertNonGermlineToSomatic()` 不被觸發），下游 `getVote()` 又有 priority bug（somatic vote 搶先 germline）+ enum vs HP tag integer literal mismatch（caller 端 `if(hpResult != HAPLOTYPE1_1)` 永不匹配，HP:i:33 永不出現）。HCC1395 5kHz 觀察到 99.9% reads 拿 HP21 的極端 imbalance，下游 ISM HP-依賴特徵全面被汙染。
> - **修復**：在 longphase 的本地 fork **`/big7_disk/liaoyoyo2001/longphase-to-mod/`** 內以**4 個 commit 漸進**修復（**InterSubMod 本 repo 無 C++ 改動**）：(A) `8b8c1fd` 加 `--pon-only-phasing` flag 解 phasing scaffold；(B) `41ff147` 重寫 `getVote()` 為兩層 + 修 enum→int literal bug；(C) `380e8d2` 補 `countINDELHaplotype()` UNDEFINED guard；(D) **V5 working tree（uncommitted）**為 `getVote()` 加 Layer 1.5 somatic fallback、為 `countSNPHaplotype()` 補對稱 alt guard。
> - **結果（HCC1395 5kHz, 15 cherry-picked sites）**：
>   - **SEQC2 calling F1 持平在統計噪音範圍**：ClairS-TO raw（無 ISM）F1=0.7166 對**所有版本完全相同**（V5 不改 caller）；ISM SuggestFilter 後 Baseline=0.7157 / V3F=0.7154 / **V5=0.7154（vs Baseline = -0.0003 噪音）**。**F1 不是衡量 V5 真實價值的指標**（詳 §4.2.2）。
>   - **真實價值 = read-level tag 品質**：clean PS blocks（11/15 sites）paired GT concordance V5 88.2% vs Baseline 74.9%（+13.3 pp）；aggregate 全 15 sites pooled V5 78.85% vs BL 72.20%（+6.65 pp）。**全基因組（PI 報告 4 § 3.7）**：clean PS V5=90.5% / BL=82.2%（+8.3 pp）。
>   - **內部健全性**：AMB% 17.5→8.0%；HP:i:33 reads 239,679→110,197（−54%）；4 項硬性 sanity check 15/15 PASS、0 violation；程式碼 +68/-36 行集中於 3 函式、`HaplotagProcess.h:66-68` 介面契約零變動。
> - **狀態**：`validated`（V5 為 production）；**V5 working tree 未 commit**、Confidence threshold 0.6 未直接驗證、7 樣本擴展未做、cnLOH 雙親同源區仍為 open issue。
> - **InterSubMod 與 V5 的關係**：**ISM 是下游消費者**，讀 V5 產出的 `tumor_tagged.bam`；ISM 自身**不**修補 phasing/haplotag 問題（ISM F1=0.0124 vs Paired 0.0909，對 TO germline FP 仍無效）。

---

## 1. 報告目的

整合 2026-04 內 v5_audit_suite 6 份核心子報告 + 4 個 MEMORY 條目 + PI 報告 4 全基因組數字，回答四件事：

1. `LongPhase-TO`（外部 binary）為何在 TO 模式出現 self-phasing；
2. **`longphase-to-mod` fork**（本地維護的 LongPhase 修改版）以哪 4 個 commit 漸進修復；
3. **InterSubMod 與這條修復鏈是什麼關係**（澄清誤解：ISM 是下游消費者，非實作者）；
4. 哪些已知 open issue（cnLOH 區、7 樣本擴展、V5 working tree commit、Confidence threshold 直接驗證）。

**讀者預設**：具備 NGS / long-read / variant calling 基本背景，但**不熟 longphase-to-mod 與 InterSubMod 的職責切分**的 PI 與工程同事。

## 2. 系統背景

| 名詞 | 定義 |
|------|------|
| **ClairS-TO** | Tumor-Only somatic variant caller（外部 binary），輸出 VCF（FILTER 含 PASS / NonSomatic / LowQual）。 |
| **LongPhase-TO**（上游） | TO 模式 phasing/haplotag 工具的上游 GitHub binary（Lin et al., *Bioinformatics* 2022）。 |
| **longphase-to-mod**（本地 fork） | InterSubMod 同層的 longphase 本地 fork @ `/big7_disk/liaoyoyo2001/longphase-to-mod/`，**獨立 git repo**。本報告討論的 4 個 commit 全部在此 fork 內，**不在 InterSubMod**。 |
| **HaplotagProcess.{h,cpp}** | longphase-to-mod 內處理 read × HP tag 對應的核心檔；位於 `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp`，**InterSubMod `src/core/` 內無此檔**。 |
| **getVote() / judgeHaplotype() / countSNPHaplotype() / countINDELHaplotype()** | longphase-to-mod 的 4 個關鍵函數；本報告 §7.2 列其修改位置。 |
| **InterSubMod (ISM)** | 本 repo 工具，**read-level epigenetic characterization**；**下游消費** longphase-to-mod 產出的 `tumor_tagged.bam`（透過 HP:i tag），自身不修改 haplotag 邏輯。 |

**上下游關係**：

```
ClairS-TO (外部 binary) ────→ VCF ────────────────────────┐
                                                          │
tumor BAM ───────────────────────────────────────────────→├─→ longphase-to-mod (本地 fork)
                                                          │     ├─ HaplotagProcess.cpp (4 commits)
                                                          │     └─ outputs:
                                                          │         ├─ phased VCF
                                                          │         ├─ tumor_tagged.bam (HP:i)
                                                          │         └─ LOH.bed
                                                          │             │
                                                          │             ▼
                                                          │     InterSubMod (本 repo)
                                                          │     ├─ 下游消費 HP:i tag
                                                          │     └─ 計算 ISM features (HP-依賴)
```

ISM 對 HP tag **強依賴**（distance metric、HPSig、HPFineNGroups 全部需要 HP），因此上游 longphase-to-mod 任何 bias 都直接傳遞至下游分析；但 ISM **不修補**上游問題，只是被動受影響。

## 3. 原本流程（Before：原始 LongPhase-TO 行為）

對齊 `01_code_diff_analysis.md` §1.3「PON / Phase / Tag 三階段」：

```
[1] PON 階段（原始 LongPhase-TO 與 V5 完全相同）
    └─ 讀 PON 資料庫（1000g-pon、CoLoRSdb、dbsnp、gnomad）做 caller-level germline / somatic 分類

[2] Phase 階段（原始 LongPhase-TO ❌ vs V5 ✅）
    │
    │ 原始 LongPhase-TO（params.ponOnlyPhasing == false）：
    │     phasing graph 用 germline + somatic + unknown 混合 anchor
    │     somatic variants 互相 phasing → self-phasing 17.3:1 bias
    │
    └─→ tumor BAM 進 [3]

[3] Tag 階段（原始 LongPhase-TO 有兩個 bug）
    │
    │ Bug 1（getVote priority）：variantKeys 順序為 {HP1_1,HP2_1} → {HP3,HP2_1} → {HP1,HP2}
    │     → 一旦有任何 somatic vote，迴圈第一輪 break，germline 完全被忽略
    │
    │ Bug 2（enum vs integer literal）：caller 端比對 hpResult != HAPLOTYPE1_1（enum=3）
    │     但 hpResult 已是 HP tag integer (11/21/33) → 永遠不匹配
    │     → fallback 死分支永不執行；HP:i:33 永不出現
    │
    └─→ tumor_tagged.bam（99.9% reads → HP21 的極端 imbalance）

[結果] InterSubMod 讀此 BAM 計算 ISM 特徵 → HP-依賴特徵全面被汙染
```

**關鍵假設（後被推翻）**：原始 LongPhase 開發時假設 phasing graph 內所有 variants 都是 germline 等價的；TO 模式 somatic-rich region 違反此假設。

## 4. 問題描述

### 4.1 觸發情境

- **樣本**：HCC1395 5kHz（V3F → V5 audit 主驗證樣本，n=1）；其他 6 樣本為後續行動建議（v5_audit_suite/00_INDEX 後續行動 #3）。
- **參數**：ClairS-TO 預設 + LongPhase-TO `--pon-only-phasing=false`（pre-V2b）。

### 4.2 錯誤現象（量化指標）

| 指標 | Baseline（pre-V2b） | V3-Fixed | V5 | 來源 |
|-----|--------------------|----------|----|----|
| 觀察到的 HP imbalance | 99.9% reads → HP21 | HP_Ratio 接近 paired baseline | 同 V3F + Layer 1.5 重分配 | 01_code_diff §2.1, V5 memory |
| HP:i:33 reads (HCC1395 全基因組) | enum bug 導致永不出現（0）| 239,679 | 110,197 (−54%) | V5 memory |
| TP Balanced% | (V2b: 13.0%) | 22.5% (+9.5pp) | (V5 同 V3F 基準) | V3F memory |
| AMB% | N/A（HP33 不存在）| 17.5% | **8.0%**（−9.5 pp） | V5 memory, 00_INDEX |
| **SEQC2 calling F1**（詳 §4.2.2）| 0.7117 / 0.7157 ⁽¹⁾ | 0.7154 | **0.7154**（vs Baseline+ISM = **-0.0003 噪音**）| PI 報告 4 §3.4 第 252-262 行 |
| ↳ ClairS-TO raw F1（不經 ISM） | **0.7166** | **0.7166** | **0.7166** ← 所有版本完全相同（V5 不改 caller） | PI 報告 4 §3.4 第 256 行 |
| ISM F1 | 0.0151 | 0.0124 | 0.0125 | V3F memory, V5 memory |

⁽¹⁾ 「Baseline F1=0.7117」（V5 memory 版本演進表）與「Baseline F1=0.7157」（PI 報告 4 §3.4 完整對比）為**同實驗不同 cut**：0.7117 為早期 ISM 版本對 baseline BAM 的結果；0.7157 為 PI 報告 4 規範化 ISM SuggestFilter 後的結果。**以 PI 報告 4 為準**（更新）。
| 4 項 sanity check | — | — | **15/15 PASS, 0 violation** | 06_sanity §7 |
| Clean PS paired GT concordance（11/15 sites）| — | — | **V5 88.2% / BL 74.9%（+13.3 pp）** | 07_paired §3 |
| Aggregate paired GT（all 15 sites pooled）| — | — | **V5 78.85% / BL 72.20%（+6.65 pp）** | 07_paired §2 |
| 全基因組 clean PS（PI 報告 4 § 3.7）| — | — | **V5 90.5% / BL 82.2%（+8.3 pp）** | 20260424_V5_vs_Baseline_complete_comparison_01.md |

> **解讀**：99.9% HP21 不是 phasing 失敗本身，是 **getVote priority bug + enum mismatch 雙重 bug 在 PON-only phasing 啟用後集中暴露**（因 PON-only 改變了 countMap 分佈，原本只在 paired 模式下不會被觀察到的 priority bug 立刻顯形 — 01_code_diff §2.1）。

### 4.2.1 self-phasing 對 LOH 的兩層影響（重要精確化）

`02_Self_Phasing根因.md`（第 108-114 行 + 第 264 行精確化）區分**兩個 LOH 層次**，本報告之前未明確區分，此處補強：

| LOH 層次 | 路徑 | self-phasing 影響 | 數字 |
|---------|-----|-----------------|------|
| **ISM HP_Ratio LOH** | BAM HP tag → ISM HP_Ratio<0.1 or >0.9 | **嚴重**：62% 是 artifact，AF 0.1-0.8 近 100%；TO TP 中 86.5% 在 paired 下完全平衡（HP_Ratio 0.4-0.6） | **62%** ISM-level LOH artifact |
| **LOH.bed region-level LOH** | VCF allele depth (AD) → LongPhase region detection | **零**：PON-only phasing 實驗 LOH.bed Jaccard=1.0，完全不受 self-phasing 影響 | LOH.bed 不變 |

- **同位點 HP_Ratio 差異 Cohen's d = -1.20**（巨大效應量；02_Self_Phasing根因.md 第 114 行）
- **94.6% somatic reads → HP1**（第 258 行）
- **6,485 個非 LOH 平衡位點 100% 出現此問題**（第 197 行）
- 兩套 LOH 系統使用不同定義（**kappa = 0.670** 的不一致性由此解釋；第 234 行）

> **正確說法**：self-phasing 的因果影響位於「haplotag → ISM HP_Ratio」這條路徑上，**非** LongPhase 的 LOH region detection；本報告 §0 TL;DR 與 §4.2 表的「62%」皆指 ISM-level LOH artifact，與 LOH.bed 不同概念。

### 4.2.2 為何 V5 「沒改 phasing 與 caller，SEQC2 F1 卻有微幅變動」？（用戶質疑的概念釐清）

**核心邏輯（PI 報告 4 §3.4 第 264-268 行）**：

```
ClairS-TO (caller) ──→ VCF (FILTER: PASS / NonSomatic / LowQual)
                           │
                           │ 此時 raw calling F1 = 0.7166（caller 層的真實 F1）
                           │ ←── V5 不改此處
                           ▼
                       longphase-to-mod (4 commits) ──→ tumor_tagged.bam (HP:i tag)
                           │
                           │ V5 改的是 BAM HP tag，不是 VCF
                           ▼
                       InterSubMod ISM ──→ 計算 region-level features (HP_Ratio, HPSig, ...)
                           │
                           │ HP tag 變了 → region 特徵變了
                           ▼
                       ISM SuggestFilter (`src/core/RegionProcessor.cpp:1120, 1269`)
                           │
                           │ 對每個 variant 計算「是否該標 LOW_QUAL」的 SuggestFilter 欄位
                           │ ←── V5 對 F1 的影響進入這裡（透過改變 BAM 而非 VCF）
                           ▼
                       套用 SuggestFilter 過濾後 → 最終 SEQC2 F1
                           │
                           │ Raw 0.7166 → Baseline+ISM 0.7157 → V3F+ISM 0.7154 → V5+ISM 0.7154
                           │ V5 vs Baseline = -0.0003（**噪音**，非真實改善）
```

**為何 F1 變動如此小（-0.0003）**：

| 證據 | 數字 | 含義 |
|------|------|------|
| ClairS-TO raw F1（所有版本） | **0.7166** 完全相同 | V5 不改 caller |
| ISM SuggestFilter 對 F1 的整體影響 | -0.0009 ~ -0.0013（負面，所有版本）| ISM 過濾本身略微負面 |
| V5 vs Baseline 的 F1 差距 | -0.0003（-0.04%）| 在統計噪音範圍 |
| V5 vs V3F 的 F1 差距 | +0.0001（+0.01%）| 同樣噪音 |
| ISM SF Precision（V5）| 113 TP / 74 FP；FP catch rate 0.63% | ISM 過濾無效於 ClairS-TO 的 germline FP |
| 所有版本 F1 區間 | 0.7154 - 0.7166（寬度 0.0012）| 任何版本之間都是噪音級差距 |

**為何 F1 不能衡量 V5 品質（PI 報告 4 直接結論）**：
- ClairS-TO 的 FP **主要是 germline variants**，非 somatic；ISM 甲基化分析設計用來區分 subclone 結構，不是區分 germline vs somatic；**V5 修正不改變此根本限制**（PI 報告 4 第 279 行）
- F1 變動 100% 來自 ISM SuggestFilter；SuggestFilter 對 ClairS-TO 的 FP 過濾力**先天不足**（所有版本都是 -0.0009 ~ -0.0013 級的微影響）
- **F1 持平不代表 tags「沒有更好」，F1 變化也不代表 tags「更好或更差」**（PI 報告 4 第 268 行原文）

**V5 真實價值的正確衡量**：
1. **read-level paired GT concordance**：clean PS +13.3 pp（15-site）/ +8.3 pp（全基因組）
2. **AMB% 修正**：17.5% → 8.0%（從過度模糊化恢復至 paired ground truth 5.4% 量級）
3. **HP:i:33 reads 重分配**：239,679 → 110,197（-54%；保留 phased 方向資訊）
4. **4 項 sanity check**：15/15 PASS、0 violation
5. **介面契約零變動**：`HaplotagProcess.h:66-68` 三函數簽章 baseline→V5 一字未變

> **用戶質疑的本質正確**：V5 確實不改 phasing graph（V2b 解的）也不改 caller。**任何 0.7117 → 0.7154 的「F1 提升」敘事都需要謹慎** — 真正提升從 baseline → V3F = +0.0040 的部分主要來自 V3F 修了 enum vs integer literal bug（讓 ISM 第一次能正確讀到 HP:i:33），不是 V5 的 Layer 1.5 帶來的。**V5 對 F1 的貢獻 = +0.0001（V5 vs V3F），噪音級**。

### 4.3 影響範圍 — ISM 特徵 3-tier 分類（self-phasing 直接 / 間接 / 不受影響）

對齊 `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md`（含 mermaid 影響分級圖）+ `InterSubMod/docs/reports/research_landscape/03_ISM分析價值界定.md` 的特徵分類。**self-phasing 並非汙染所有 ISM 特徵 — 55% 特徵完全不受影響，可作為 V5 之前 archive 的安全結論**：

![Self-Phasing 對 ISM 特徵 3-tier 影響分級總覽（嚴重 38% / 中度 7% / 不受影響 55%）](../../../research_landscape/figures/03_self_phasing_impact.png)

*Figure 0 — Self-Phasing 影響量化總覽（`InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md`）。整合 Somatic HP bias 17.3:1、ISM HP_Ratio LOH 62% artifact、跨樣本 7/7 一致性等指標於單一視圖。*

| 影響等級 | 特徵數（佔比） | 代表特徵 | TO 結果處理 |
|---------|------------|---------|-----------|
| 🔴 **嚴重影響**（直接依賴 HP）| **29 個（38%）** | HP_Ratio → 假 LOH；Potential_LOH → 62% artifact；HPMergedDelta/Sig → 方向反轉；hp_assign_rate → 偏高；effective_hp_reads → 偏離；HPFineNGroups（含 NG=2 LOH-constrained phasing） | **必重跑** V5 BAM；舊結論需加註版本 |
| 🟡 **中度影響**（間接污染）| **14 個（7%）** | QualityScore → AUC 0.497（已移除 LOH penalty）；GlobalP → 取 HP/Allele 最小值（HP 噪音可能偶然壓低）；CramersV → 取 HP/Allele 最大值；VerificationClass → label_sig 含 HP 成分 | 重評；多數影響微弱或已程式碼移除 |
| 🟢 **不受影響**（無 HP 依賴）| **42 個（55%）** | PairwiseMean/MedianDist（全 reads 不分 HP）；AlleleDelta / AlleleP（只用 ALT/REF）；Caller 特徵（AF / GQ / DP / SB 來自 VCF）；甲基化矩陣（BAM MM/ML tag）；CpG 座標、region_methyl_mean | **結論穩固，不需重測** |

**對既有結論的具體影響**：

- 🔴 `HPFineNGroups subclone marker ⭐4 → ⭐3 降級`（2026-04-23 週報）— 機制部分歸於 phasing signature 而非 methylation；待 V5 + flag=on 重驗（後續 F4）
- 🔴 `LOH-constrained phasing discovery (NG=2, 6/6 樣本 Inner 93-99% same-hap)`（2026-04-22 新發現）— 已重新詮釋為 phasing artifact 主導
- 🔴 歷史 archive 數據：所有 V5 之前產出的 ISM **HP-依賴**結果需重跑或加註版本（manifest.yaml `haplotag_version` 欄位 — 後續 F5）
- 🟢 paired-pure delta F1=+0.0112 等與 HP tag 無關的結論不受影響

> 詳細特徵分類表見 `InterSubMod/docs/reports/research_landscape/03_ISM分析價值界定.md`；mermaid 影響分級圖見 `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md` Section「Self-Phasing 對 ISM 各輸出欄位的影響」。

### 4.4 跨樣本一致性與 Simpson's Paradox 釐清（7/7 樣本驗證）

對齊 `InterSubMod/docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md`（CV-1 ~ CV-5 跨樣本驗證）+ `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md` 量化證據總表：

| 跨樣本指標 | 數值 | 解讀 |
|----------|------|------|
| **CV-2 方向一致性** | **7/7 樣本全部觀察到相同方向的 self-phasing 效應** | 排除樣本特異性可能 |
| 同位點 HP_Ratio 跨模式 r | **r = 0.001**（n=288K pairs）| TO 與 paired 在同一變異上的 HP 判定**完全不相關**（不是噪音，是無訊號）|
| TO-only LOH 在 paired 下完全平衡 | **86.5%**（HP_Ratio 0.4-0.6）| TO 看到的 LOH 在 paired 下大多是平衡的 → TO LOH 主要是 artifact |
| Cohen's d (HP_Ratio 差異) | **-1.20** | 巨大效應量（一般 \|d\|≥0.8 即視為大）|
| CV-1 Simpson's paradox | **r = -0.964** 跨樣本 imbalance vs self-phasing fraction | **強負相關，非預期的正相關** |

**Simpson's Paradox 機制釐清**（causal chain report 第 311 行）：
- 結構性 LOH 比例高的樣本（HCC1937、HCC1395）有高整體 imbalance
- 但 self-phasing **fraction** 反而較低（被結構性 LOH 稀釋）
- **不否定 self-phasing 機制** — 是「加法關係非替代關係」：self-phasing 仍是 TO ISM HP_Ratio LOH 的主因（62% artifact），但在樣本層級結構性 LOH 是 imbalance 的更強預測因子
- CV-5：4/7 樣本 self-phasing fraction > 0.30（**PARTIAL** — HCC1395 / HCC1395_DORADO / HCC1937 結構性 LOH 主導，self-phasing 雖存在但不是該樣本 imbalance 的主要驅動）

## 5. 根本原因（多鏈 5 Whys）

self-phasing 不是單因果鏈，是**三層獨立 bug 在 PON-only 模式啟用後集中暴露**。對齊 01_code_diff 拆三鏈：

### 5.A Phase 層 — self-phasing scaffold（V2b 解）

| 層 | 提問 | 答覆 |
|---|------|------|
| 1 | 為何 phasing graph 把 somatic 當 anchor？ | 因 `params.ponOnlyPhasing == false`（預設）→ 不呼叫 `convertNonGermlineToSomatic()` |
| 2 | 為何預設 false？ | 原始 LongPhase 對 paired 模式設計，假設 phasing graph 內所有 variants germline 等價 |
| 3 | 為何此假設在 TO 失敗？（**根因**）| TO 無 normal control，somatic-rich region 違反等價假設；somatic ALT reads 集中於子族系一條 hap → 自我定相 |

#### 5.A.1 「somatic ALT reads 集中於子族系一條 hap」的 3 層獨立證據鏈

對應 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` Section 5（含 6-panel 概念圖與 IGV 真截圖）：

> **⚠ 圖片渲染說明**：以下表格的「視覺證據」欄為**檔案路徑**（用於追蹤來源）；**實際 inline 渲染**請見表格下方 §5.A.1.1 ~ §5.A.1.3 三個獨立圖片區塊。
> Markdown 表格儲存格內若塞 `![alt](path)` 會因儲存格寬度被壓縮（GitHub / VSCode / Obsidian 渲染不一致），無法清楚閱讀；本報告慣例：**表格列路徑作 source citation，獨立區塊用 `![](...)` 完整渲染**。

| 層 | 量化指標 | 數值 | 視覺證據（路徑）|
|---|--------|------|--------|
| **理論層** | Phasing graph edge weight：`weight(A,B) = Σ_reads I(read 帶 A.alt) × I(read 帶 B.alt)` | 同 clone somatic 在同 reads **共現 ≈ 100%**（共享 sub-population）；germline het 隨機分散 ≈ 50% → somatic-somatic edges 比 germline 權重更高 → somatic 反客為主定義 scaffold | `InterSubMod/docs/reports/pi_reports/2026/04/figures/fig2_self_phasing_concept.png`、`fig3_af03_walkthrough.png`（AF=0.3 具體走例）|
| **全基因組層** | HCC1395 baseline TO HP1:HP2 reads | **614,000 : 35,500 = 17.3 : 1**（**94.6%** 集中於 HP1，跨 23 染色體）vs 隨機預期 ~1:1 → 生物學上不可能，必為 artifact | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01d_somatic_bias_explanation.png` Panel D；`fig4_evidence_summary.png` |
| **個別位點層**（IGV 6-BAM 並列） | 3 個 SP-extreme 位點 baseline HP2:HP1 ratio | SP1 chr19:17565944 = **113:0**；SP2 chr19:12452332 = **109:1**；SP3 chr19:12467180 = **108:0** — baseline 與 paired 方向**完全相反**；V5 修正後與 paired 一致（3/3） | `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP{1,2,3}_*.png`（baseline / V2b / V3F / V5 / Paired tumor / Paired normal 6-BAM 並列）|

**3 層獨立但結論一致** — self-phasing artifact 真實存在；17.3:1 為 average，個別位點可達完全失衡（113:0），average 已被「不那麼極端的位點」稀釋。

##### 5.A.1.1 理論層視覺：phasing graph self-phasing 機制

![Self-Phasing 概念流程：同 clone somatic 共享 reads → phasing graph 強連結 → 全部塞同一 phase block](../../../pi_reports/2026/04/figures/fig2_self_phasing_concept.png)

*Figure 1 — Self-phasing 概念流程（PI 報告 1 Section 2）。同一 sub-clone 的 somatic variants 共享 sub-population reads，long read 跨多個 somatic variants 都帶 ALT，phasing graph 看到強連結 → 全部塞同一 phase block (HP1)。*

![AF=0.3 走例：30% ALT reads 在 paired 模式隨機分到 HP1/HP2，但 TO 模式因 self-phasing 全部指向 HP1](../../../pi_reports/2026/04/figures/fig3_af03_walkthrough.png)

*Figure 2 — AF=0.3 具體走例（PI 報告 1 Section 2.2）。同一個 somatic variant 在 paired 模式 HP_Ratio ≈ 0.5（隨機），在 TO 模式 HP_Ratio → 0.94（偏 HP1）。*

##### 5.A.1.2 全基因組層視覺：HCC1395 17.3:1 量化證據（6-panel 概念圖）

![Somatic Bias 17.3:1 — 6-panel 概念與實證視覺化（A: 預期分布 / B: 實測 self-phasing / C: 機制推導 / D: 全基因組 614K vs 35.5K / E: 個別極端位點 / F: V5 修復後）](../../../pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01d_somatic_bias_explanation.png)

*Figure 3 — `v5_audit_suite/10_somatic_bias_explanation.md` Section 2 主圖。Panel D 顯示 HP1=614K vs HP2=35.5K bar chart；Panel E 列 SP1/SP2/SP3 的 113:0 / 109:1 / 108:0 個別失衡；Panel F 顯示 V5 修復後 ~1:1。*

![量化證據總表：Self-Phasing 影響 7 項指標（17.3:1、94.6%、62% LOH artifact、Cohen's d=-1.20、7/7 一致）](../../../pi_reports/2026/04/figures/fig4_evidence_summary.png)

*Figure 4 — 量化證據總表（PI 報告 1 Section 3）。整合 Somatic HP bias、ISM HP_Ratio LOH artifact 比例、跨樣本一致性等 7 項證據。*

##### 5.A.1.3 個別位點層視覺：3 個 SP-extreme 位點 IGV 6-BAM 並列截圖

> 每張截圖由上至下 6 個 panel：**baseline / V2b / V3-Fixed / V5 / Paired tumor / Paired normal**。觀察重點 — baseline 的 HP 主導方向與 paired tumor **完全相反**（self-phasing artifact 在單一位點的實證），V5 修正後與 paired 一致。

![SP1 chr19:17565944 — baseline 113:0 全在 HP1，paired 與 V5 都是 HP2 主導](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png)

*Figure 5 — SP1（chr19:17565944）。baseline panel reads 全部集中於 HP1+HP1-1（粉紅+淡綠群）；V5 / V2b / V3-Fixed reads 整體翻轉至 HP2+HP2-1；Paired tumor 確認 HP2 為真實方向。*

![SP2 chr19:12452332 — baseline 109:1，V5 翻轉至 HP2 主導對齊 paired](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP2_chr19_12452332.png)

*Figure 6 — SP2（chr19:12452332）。baseline HP1+HP1-1 集中 109 reads，HP2 stack 僅 1 read；V5 方向翻轉至 HP2+HP2-1，與 paired tumor 一致。*

![SP3 chr19:12467180 — baseline 108:0 與 SP1 相同模式](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP3_chr19_12467180.png)

*Figure 7 — SP3（chr19:12467180）。與 SP1 同模式 — baseline HP1 主導 → V5 HP2 主導，HP orientation 整體翻轉，V5 與 paired 一致。*

> **完整 6-panel 概念圖細節 + 3 IGV 截圖視覺解讀**詳見 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` Section 2-3（Panel A-F 逐項 + 4-BAM 視覺觀察彙整表）。

### 5.B Tag 層 — getVote priority bug（V3F 解）

| 層 | 提問 | 答覆 |
|---|------|------|
| 1 | 為何 99.9% reads 拿 HP21？ | `getVote()` 迴圈順序 `{HP1_1,HP2_1}` 在 `{HP1,HP2}` 之前；任一 somatic vote 觸發 break |
| 2 | 為何順序這樣設計？ | baseline 用 `std::map<int,int>` + `std::vector<std::pair<int,int>>` 字面量初始化；順序是 enum 數值升序的副作用，非有意 |
| 3 | 為何 paired 模式沒暴露？（**根因**）| paired 模式 germline 訊號充足、somatic 比例小，迴圈通常在 germline pair 命中；TO PON-only 後 germline votes 急減，bug 立刻顯形 |

### 5.C Tag 層 — enum vs integer literal mismatch（V3F 同 commit 解）

| 層 | 提問 | 答覆 |
|---|------|------|
| 1 | 為何 HP:i:33 永不出現？ | caller 端 `if(hpResult != HAPLOTYPE1_1 && hpResult != HAPLOTYPE2_1)` fallback 永不觸發 |
| 2 | 為何不觸發？ | `HAPLOTYPE1_1=3, HAPLOTYPE2_1=4`（enum，`Util.h:21-25`），但 `getVote()` 已寫入 HP tag integer (11/21/33)；型別語意失配 |
| 3 | 為何寫成這樣？（**根因**）| 函數內部用 enum 值，函數輸出當 BAM HP tag integer，**未統一單一語意層** — typed enum 缺失與 integer literal 並用的歷史包袱 |

### 5.D Tag 層 — V3F 過於保守（V5 working tree 解）

| 層 | 提問 | 答覆 |
|---|------|------|
| 1 | 為何 V3F 後仍有 17.5% AMB? | V3F 在 `germlineHP1=germlineHP2=0` 時一律 `hpResult=33`，丟失 `HAPLOTYPE1_1/HAPLOTYPE2_1` 的 phased 方向資訊 |
| 2 | 為何丟失？ | V3F 採嚴格「germline first」原則，未利用 somatic enum 後綴 `_1` 本身已含的 hap 標記 |
| 3 | 為何不利用？（**根因**）| 設計取捨：V3F 優先**確保不誤標**，V5 以 Layer 1.5 + Confidence threshold 0.6 在「不誤標」前提下榨出剩餘訊息 |

### 5.E 因果鏈獨立驗證（5 條獨立路徑）

> **⚠ 表格內路徑為 source citation；實際視覺證據已 inline 渲染於 §5.A.1.1 ~ §5.A.1.3（Figures 1-7）。**

| # | 證據 | 數值／結論 | 對應視覺 | 來源 |
|---|------|-----------|---------|------|
| 1 | **PON-only VCF 層修正** | LOH.bed Jaccard=1.0、Somatic bias 17.3:1 → 消除、Phase block N50 +99.7%、Phased rate +23.6 pp、執行時間 1.36× 快 | Figure 3 Panel F | `memory/project_pon_only_phasing_verification.md`、`InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md` PON-only 修正測試 |
| 2 | **全基因組 17.3:1 量化** | HCC1395 HP1=614K / HP2=35.5K（94.6% to HP1）；跨 23 染色體一致；不可能為真實生物學現象 | Figure 3 Panel D, Figure 4 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` §2 Panel D、PI 報告 1 § 3.1 |
| 3 | **個別位點 IGV 6-BAM 真截圖** | SP1=113:0 / SP2=109:1 / SP3=108:0 baseline 與 paired 方向相反；V5 修正後與 paired 一致 3/3 | **Figures 5, 6, 7** | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` §3 |
| 4 | **V5 working tree sanity check** | 4 項硬性檢查 15/15 PASS、0 violation；守恆律 A/B + Layer 1.5 期望 1/2 + germline → HP33 違規=0 | （fig06a/06b in audit suite） | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md` §7 |
| 5 | **程式碼層最小必要** | +68/-36 行集中 3 函式、`HaplotagProcess.h:66-68` 介面契約零變動、無 over-engineering（無新 enum / 新 member / 新 std 容器 / 新 logging）| （fig01a/01b/01c in audit suite） | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` §4-§5 |

**獨立性說明**：5 條路徑各從不同角度（VCF 層、全基因組統計、單一位點視覺、軟體 sanity、程式碼結構）獨立驗證 self-phasing 機制與 V5 修補的有效性 — 任一條失敗都不會推翻其他四條，但 5 條全部通過構成穩固證據（穩定度 4-5/5）。

## 6. 修改方向（4-commit 漸進，已決定的歷史路徑）

不同於一般技術報告「列候選方案 → 選一條」，本案是 **4 commit 已序列實施**，每 commit 對應 §5 的一條鏈：

| Commit | 解的鏈 | 關鍵修改 | 驗證 |
|--------|--------|---------|-----|
| **`8b8c1fd` (V2b)** | 5.A Phase scaffold | 加 `--pon-only-phasing` flag；新 helper `convertNonGermlineToSomatic` / `resetNonPonOrigin` / `syncPhasingResultOrigins` | 後處理 GT 乾淨；HaplotagProcess.cpp 未動 |
| **`41ff147` (V3-Fixed)** | 5.B getVote priority + 5.C enum mismatch | `getVote()` 重寫為兩層 germline-first；caller 端比對改 integer literal 11/21/33 | HP_Ratio 從 99.9% 偏 HP21 恢復；HP:i:33 開始正確出現（6,793 in HCC1395） |
| **`380e8d2` (INDEL guard)** | UB 補洞 | `countINDELHaplotype()` 兩分支補 `HAPLOTYPE_UNDEFINED` guard | 消除 PON-only somatic-only INDEL 站點的 `countMap[-1]++` UB |
| **V5 working tree（uncommitted）** | 5.D V3F 保守 | (1) `getVote()` 加 Layer 1.5 somatic fallback；(2) `countSNPHaplotype()` 加對稱 alt guard | AMB% 17.5→8.0%、HP33 reads −54%、4 項 sanity 15/15 PASS |

### 6.1 替代候選（已被 audit suite 排除）

| 候選 | 為何不採 | 來源 |
|-----|---------|-----|
| 替換 LongPhase 為 WhatsHap / HapCUT2 | 風險高、外部依賴、整合工作量大；4 commit 漸進已 +13.3 pp clean PS gain，cost-benefit 不對等 | (本報告判讀) |
| 在 InterSubMod 自己加 haplotag 邏輯 | 介面割裂；ISM 設計為下游消費者，不應接 phasing 責任；ISM F1=0.0124 證明 ISM 對 TO germline FP 無修復力 | V5 memory, V3F memory |

## 7. 修改內容（雙語）

### 7.1 非工程版（給 PI / 外部讀者）

> 修補的對象**不是 InterSubMod**，而是另一個獨立的工具 `longphase-to-mod`（longphase 的本地修改版）。我們在那邊以**4 次小修補**漸進解決 self-phasing：第 1 次（V2b）讓上游不再把體細胞變異當定相錨點；第 2 次（V3F）修兩個程式碼 bug — 一個讓「無法判定」標籤永遠不出現、一個讓投票優先順序錯亂；第 3 次（INDEL guard）補一個記憶體越界的洞；第 4 次（V5）讓「投票打平」的 reads 不要全部一刀切標「不確定」，而是還能再榨一點方向訊號。
>
> **InterSubMod 自己沒改任何程式碼**。它讀新版 longphase-to-mod 產出的 BAM 後，下游分析（甲基化距離、HP 群組）就因為輸入更乾淨而連帶受惠。

### 7.2 工程版（給工程同事）

| 維度 | 變更 |
|-----|-----|
| Repo | `/big7_disk/liaoyoyo2001/longphase-to-mod/`（**獨立 git repo，非 InterSubMod**） |
| **介面契約** | `HaplotagProcess.h:66-68` 三個 method signature 從 baseline 到 V5 **一字未變**（01_code_diff §5.1） |
| **總修改量** | **+68 / -36 行**，集中於 3 函式（`getVote()` / `countSNPHaplotype()` / `countINDELHaplotype()`） |
| **Commit 1（`8b8c1fd` V2b）** | 加 `--pon-only-phasing` flag；新 helpers `convertNonGermlineToSomatic` / `resetNonPonOrigin` / `syncPhasingResultOrigins`；改在 `Phasing.cpp` (+9/-2)、`PhasingGraph.cpp` (+34/-0)、`PhasingProcess.cpp` (+25/-3)；**`HaplotagProcess.cpp` 未動** |
| **Commit 2（`41ff147` V3-Fixed）** | `getVote()` 重寫為兩層（`HaplotagProcess.cpp:506-541`，+36/-25 行）；caller 端 `judgeHaplotype()` 內 `if(hpResult != HAPLOTYPE1_1)` 改用 integer literal 11/21（`HaplotagProcess.cpp:697`） |
| **Commit 3（`380e8d2` INDEL guard）** | `countINDELHaplotype()`（`HaplotagProcess.cpp:497-510`）兩分支各加 `if(... != HAPLOTYPE_UNDEFINED)`，+8/-4 行；模式 mirror SNP 端既有 guard |
| **V5 working tree (uncommitted)** | (1) `getVote()` Layer 1.5（`HaplotagProcess.cpp:512-563`，+15/-1 行）— `else if (somaticHP1 > 0 \|\| somaticHP2 > 0)` 純分支插入；(2) `countSNPHaplotype()`（`HaplotagProcess.cpp:489-494`，+9/-6 行）對稱 alt guard `if(haplotypeBase.altHaplotype != HAPLOTYPE_UNDEFINED)` |
| **HP tag 對應表（V5 final）** | `(germlineResult, somaticTotal>0) → hpResult`：(0, F)→0 / (1, F)→1 / (2, F)→2 / (1, T)→11 / (2, T)→21 / (0, T)→33（01_code_diff §3.3） |
| **未引入** | 新 enum、新 class member、新 method signature、新 std 容器、新 logging、新 test path（01_code_diff §4.1） |
| **InterSubMod 改動** | **無**（本 repo `src/core/` 不包含 HaplotagProcess；`/big7_disk/liaoyoyo2001/InterSubMod/src/core/` 列出 20 檔，0 個 phasing/haplotag 相關） |
| **未 commit 修補** | V5 = `380e8d2` + 兩塊 working-tree 修改；audit suite 後續行動 #1 建議切 2 獨立 commits：`feat(haplotag): Layer 1.5 somatic fallback in getVote()` + `fix(haplotag): guard countSNPHaplotype against UNDEFINED on alt path` |

## 8. 新舊流程比較

| 維度 | Baseline（pre-V2b） | V3-Fixed（`41ff147`） | V5（V3F + working tree）| 來源 |
|------|-------|---------|-----|------|
| `--pon-only-phasing` | false（預設） | true（V2b 起）| true | 01_code_diff §1.3 |
| Phasing graph anchor | germline + somatic + unknown 混合 | 僅 PON-confirmed germline | 同 V3F | 同上 |
| `getVote()` 邏輯 | for-loop + `variantKeys` 順序依賴（priority bug）| 兩層顯式（germline first → somatic annotate） | 三層（Layer 1 → Layer 1.5 fallback → Layer 2 encode）| 01_code_diff §3.2 |
| HP:i:33 寫入機制 | enum bug（永不寫）| integer literal（恢復寫入，HCC1395=6,793）| 同 V3F + Layer 1.5 重分配（HCC1395 全基因組: 239,679 → 110,197, −54%）| V3F memory, V5 memory |
| AMB% | N/A | 17.5% | 8.0% | V5 memory |
| **SEQC2 calling F1（ClairS-TO raw, 不經 ISM）** | **0.7166** | **0.7166** | **0.7166** ← V5 不改 caller | PI 報告 4 §3.4 |
| SEQC2 F1（套 ISM SuggestFilter 後）| 0.7157 | 0.7154 | 0.7154（vs Baseline = **-0.0003 噪音**）| PI 報告 4 §3.4 |
| Clean PS paired GT（15-site, 11 clean） | — | — | V5 88.2% vs BL 74.9%（+13.3 pp） | 07_paired §3 |
| 全基因組 clean PS（PI 報告 4） | — | — | V5 90.5% vs BL 82.2%（+8.3 pp）| 20260424_V5_vs_Baseline |
| 4 項 sanity check（15 sites） | — | — | 15/15 PASS, 0 violation | 06_sanity §7 |

## 9. 驗證方式

對齊 v5_audit_suite Agent D（`06_v5_sanity_bug_check.md`）實際採用的 4 項硬性檢查 + 後續行動建議：

```
1. 守恆律 A · Δ-consistency
   → 驗證: 對 15 sites 計算 ΔHP33 + (ΔHP11 + ΔHP21) 是否 = 0
   → 結果: 15/15 PASS（06_sanity §2）

2. 守恆律 B · Germline 不變
   → 驗證: V3F 與 V5 的 HP1 / HP2 read 數逐 site 比對
   → 結果: 15/15 PASS, 0 reads 漂移（06_sanity §3）

3. Layer 1.5 期望 1 · 33→directional 精確守恆
   → 驗證: 對每 site 計算 ΔHP11 == n(V3F=33→V5=11)、ΔHP21 == n(V3F=33→V5=21)
   → 結果: 15/15 PASS（06_sanity §4）；C_V5max1 = 34 reads / C_V5max2 = 22 reads / C_V5max3 = 12 reads 全部精確守恆

4. Layer 1.5 期望 2 · 無 germline → HP33
   → 驗證: 跨 15 sites pool transition table；germline (HP1/HP2) → HP33 的 read 數
   → 結果: 0 reads（06_sanity §5）

5. SEQC2 F1 不退步
   → 驗證: ./scripts/run_vcf_all_snv.sh --mode all-with-w5000；F1 ≥ 0.7150
   → 結果: 0.7154 ≥ 0.7153（V5 memory）

6. Paired ground-truth concordance
   → 驗證: 15 sites × per-read HP 對 paired tumor BAM HP:Z 字串首段（含 PS-block orientation correction）
   → 結果: clean PS（11 sites）V5=88.2% / BL=74.9%；aggregate V5=78.85% / BL=72.20%（07_paired §2-§3）

7. ⚠ Confidence threshold 0.6 直接驗證
   → 狀態: 未做 — 需 V5 binary 加 vote log 才能直接驗證；間接證據（PS-block 內方向一致、無鋸齒）齊全（06_sanity §6）
```

> **統計顯著性**：本報告未直接呼叫 `/auc-confound-guard`（V5 改動屬 haplotag 邏輯而非特徵 AUC 宣告）。下游 HPFineNGroups marker 在 V5 後的重驗**必須**走 `/auc-confound-guard` 三關。

## 10. 影響範圍

| 受影響對象 | 影響 |
|-----------|-----|
| **InterSubMod 本 repo** | **無 C++ 改動**；`src/core/` 無 phasing/haplotag 檔；ISM 為被動受惠者 |
| longphase-to-mod fork | `HaplotagProcess.cpp` 集中改 +68/-36 行；介面契約零變動 |
| 使用者（PI / 合作）| 新版 ISM 報告需註明 V5；舊報告需加版本標籤 |
| BAM HP tag | HP:i 仍為 1/2/11/21/33 五值；HP:i:33 比例顯著下降（HCC1395 全基因組 −54%） |
| ISM 下游 HP-依賴特徵 | 必須在 V5 BAM 上重跑（HPSig、HPFineNGroups、HP_Ratio_norm、SampleASM_Delta） |
| 舊資料 archive | 重跑或標 `haplotag_version: V3F / V5`；建議在 manifest.yaml 加此欄位 |
| ClairS-TO 上游 | 無影響（V5 修補在 longphase 階段，VCF 已產出後） |
| 部署 | 需 rebuild longphase-to-mod；無新外部依賴；InterSubMod 不需重 build |

## 11. 風險與限制

| ID | 風險／限制 | 來源 |
|----|----------|-----|
| **R1** cnLOH 雙親同源 region | 兩 hap 序列接近相同，V5 fallback 也無法區分；需 cnLOH-aware filter（CN 層 + germline-only methylation reference） | (定性論述) |
| **R2** AF > 0.9 真結構性 LOH 與 self-phasing 假性 LOH 的分離邊界 | V5 在 0.1–0.8 區間表現好；AF > 0.9 的真實 LOH 與 phasing artifact 仍可能共存 | (定性論述) |
| **R3** 15 sites cherry-picked | audit suite 為 15 自選位點（5 TP / 4 FP / 3 V5-reassign / 3 SP-extreme），不是隨機抽樣；全基因組數字參照 PI 報告 4（V5=90.5% / BL=82.2%, +8.3 pp） | 00_INDEX caveat #1 |
| **R4** Confidence threshold 0.6 未直接驗證 | 06_sanity §6 標 "需 V5 binary 加 vote log"；目前間接證據（PS-block 內方向一致）齊全但非 hard evidence | 06_sanity §6, 00_INDEX caveat #3 |
| **R5** V5 working tree uncommitted | V5 = `380e8d2` + 兩塊未 commit 修改；可追溯性風險 | 01_code_diff §5.2, 00_INDEX caveat #4 |
| **R6** 7 樣本擴展未做 | 僅 HCC1395 5kHz 一樣本驗證；其他 6 樣本 V5 全量重跑為「中優先後續」 | 00_INDEX 後續行動 #3 |
| **R7** 部分 metric V5 略遜 | L4 orientation-corrected per-site count V5 wins=2 / BL wins=9（受 problem PS 影響）；imbalance ratio mean +0.014（受 B_FPA2 outlier 主導，移除後 ≈持平）| 02_concordance, 04_imbalance §3 |
| **R8** Problem PS / low_germ_n 不適用 read-level 對齊 | A_TP02、D_SP1、D_SP3 在 problem PS 上 BL 看似較好為 PS orientation pick + germline_n<5 的 artifact；非 V5 退步 | 07_paired §4 |
| **R9** Paired 自身仍用 LongPhase | paired-mode ground truth 的 phasing 也由 LongPhase 演算法產生；理想需 trio-phased（normal+tumor+normal）作 second ground truth | 07_paired §5d |

## 12. 後續工作

對齊 audit suite `00_INDEX.md` 後續行動 + 本 fact-check 新增項：

| ID | 動作 | 優先 | 來源 |
|----|------|-----|-----|
| F1 | **commit V5 working-tree 修改**（切 2 獨立 commits）| **高** | 00_INDEX 後續 #1, 01_code_diff §5.2 |
| F2 | 追加 Confidence threshold 0.6 投票 log 驗證（V5 binary 加 log 或 IGV session 看 PS block ALT/REF 投票）| 中 | 00_INDEX 後續 #2, 06_sanity §6 |
| F3 | 7 樣本 V5 BAM 全量重跑（HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829）| 中 | 00_INDEX 後續 #3 |
| F4 | P4 master dataset × 兩 flag 重跑驗證 HPFineNGroups marker | 低（依發表計畫）| 00_INDEX 後續 #4, `memory/project_hpfinengroups_subclone_marker.md` |
| F5 | manifest.yaml 加 `haplotag_version: V3F / V5` 欄位（archive 治理）| 中 | (本報告新增) |
| F6 | cnLOH 區獨立評估方案（CN 層 + germline-only methylation reference）| 中 | `memory/project_research_strategy_2026Q2.md` |
| F7 | trio-phased（normal+tumor+normal）second ground truth | 低 | 07_paired §5d |
| F8 | 跨 50-100 隨機抽樣位點 cross-validate 15-site cherry-picked 結論 | 中 | 07_paired §5d |

## 13. 結論

原始 LongPhase-TO 在 ClairS-TO TO 模式下的 self-phasing 不是單一 bug，是 **Phase 層 scaffold + Tag 層 priority bug + enum-int mismatch + V3F 過度保守** 四鏈獨立故障在 PON-only 啟用後集中暴露；修補在 **`longphase-to-mod` fork（非 InterSubMod）**以 4 個 commit 漸進完成（`8b8c1fd` → `41ff147` → `380e8d2` → V5 working tree），總 +68/-36 行集中於 3 函式、介面契約零變動。**V5 不改 ClairS-TO caller（raw F1=0.7166 所有版本相同）也不改 phasing graph（V2b 解的）；V5 改 BAM HP tag → 透過 ISM SuggestFilter 影響最終 F1，但 V5 vs Baseline F1 = -0.0003（噪音級，非真實改善）— SEQC2 F1 不是衡量 V5 品質的指標**（詳 §4.2.2）。V5 真實價值在 read-level：HCC1395 5kHz 15 sites 4 項硬性 sanity check 15/15 PASS、clean PS paired GT concordance V5 88.2% vs BL 74.9%（+13.3 pp）、aggregate +6.65 pp、AMB% 17.5→8.0%；全基因組（PI 報告 4）clean PS V5=90.5% / BL=82.2%（+8.3 pp）。**InterSubMod 在這條修補鏈是下游消費者而非實作者**，本 repo 無 C++ 改動，ISM F1=0.0124 對 TO germline FP 仍無修復力（FP 主要是 germline variants，ISM 甲基化分析設計用於 subclone 結構而非 germline/somatic 區分）。當前 V5 為 production，但 V5 working tree 未 commit、Confidence threshold 未直接驗證、7 樣本擴展未做、cnLOH 雙親同源區仍為 open issue 為**短期優先後續**。

---

## 附錄 A：引用清單

### A.1 文件（皆以 `InterSubMod/...` 前綴；longphase-to-mod 為獨立 repo 用絕對路徑）

- `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md`
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md`（19 子報告母索引）
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md`（V3F + V5 commits 逐 diff）
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/04_imbalance_ratio_analysis.md`
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md`（4 項硬性檢查 15/15 PASS）
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md`（範本來源）
- `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md`（5 PI 問題答覆）
- `InterSubMod/docs/reports/pi_reports/2026/04/20260424_V5_vs_Baseline_complete_comparison_01.md`（PI 報告 4 全基因組數字）
- `InterSubMod/docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md`
- `InterSubMod/docs/references/manual/20260402_longphase_to_phasing_quality_literature.md`

### A.2 程式碼（皆於獨立 repo `/big7_disk/liaoyoyo2001/longphase-to-mod/`）

- `HaplotagProcess.cpp:484-495`（V5 `countSNPHaplotype()` + alt guard）
- `HaplotagProcess.cpp:497-510`（V5 `countINDELHaplotype()` + UNDEFINED guard）
- `HaplotagProcess.cpp:512-563`（V5 `getVote()` Layer 1 / 1.5 / 2）
- `HaplotagProcess.cpp:697`（V3F 修正後 caller 端 integer literal 比對）
- `HaplotagProcess.cpp:725`（`getVote()` caller in `judgeHaplotype()`）
- `HaplotagProcess.h:66-68`（三函數簽章 — baseline 到 V5 一字未變）
- `Util.h:20-25`（`Haplotype` enum 定義 - HAPLOTYPE_UNDEFINED=-1, HAPLOTYPE1=1, HAPLOTYPE2=2, HAPLOTYPE1_1=3, HAPLOTYPE2_1=4, HAPLOTYPE3=5）
- `PhasingProcess.cpp:55, 154-157`（PON 設定 + `convertNonGermlineToSomatic` 觸發點）

**InterSubMod 內**（下游消費，無修改）：
- `InterSubMod/scripts/analysis/v5_imbalance_improvement.py`（Agent C 產出）
- `InterSubMod/scripts/analysis/v5_sanity_paired_check.py`（Agent D 產出）
- `InterSubMod/scripts/analysis/build_s0_self_phasing_*.py`（多腳本）

### A.3 MEMORY

- `memory/project_v5_somatic_fallback_verification.md`（當前最佳 haplotag）
- `memory/project_v3_fixed_haplotag_verification.md`（V3F 含 commit `41ff147`、`380e8d2`）
- `memory/project_self_phasing_causal_chain_confirmed.md`
- `memory/project_pon_only_phasing_verification.md`（V2b LOH.bed Jaccard=1.0）
- `memory/project_hpfinengroups_subclone_marker.md`（⭐4→⭐3 降級）
- `memory/project_loh_constrained_phasing_discovery.md`（NG=2, 6/6 樣本）
- `memory/feedback_feature_name_vs_definition_rule.md`（**本次校正觸發的規則** — 必查 C++ 原始碼）

### A.4 外部論文

- Lin J-H et al. (2022). LongPhase: an ultra-fast chromosome-scale phasing algorithm. *Bioinformatics* 38(7): 1816–1822. https://github.com/twolinin/longphase

### A.5 框架方法學（本報告結構來源）

- Toyota A3 Thinking — `references/frameworks_cheatsheet.md` §A
- ADR (Nygard, 2011) — 同上 §B
- Google SRE Postmortem — 同上 §A
- Diátaxis（§7 雙語層次）— 同上 §D

---

## 附錄 B：變更歷史

| 日期 | 變更 | 觸發 | 作者 |
|-----|------|-----|-----|
| 2026-04-29 00:55 | 初稿（structured-tech-report skill 首個示範案例）| 用戶請求建立報告 skill | structured-tech-report v1.0 |
| 2026-04-29 02:30 | 用戶質疑 V5 歸屬 | 「V5 應該是 longphase-to 的 HaplotagProcess 吧」 | 用戶 |
| 2026-04-29 02:50 | RCA 認知確認；列 30 項待修清單 | 自我審查 | structured-tech-report skill |
| **2026-04-29 03:15** | **Fact-check 全面校正** — 整體重寫對齊 v5_audit_suite 6 子報告 ground truth；修正：(1) V5 歸屬至 longphase-to-mod fork；(2) 4-commit 演進敘事取代「V5 一次解所有」；(3) 真實函數名／行號／commit hash；(4) 刪除捏造的「Phase block N50 11.9→1.2 Mbp」「test_haplotag_v5.cpp」「7 樣本 ≥0.78」；(5) 引入 R3-R9 真實 caveat（15 sites cherry-picked、Confidence threshold 未驗證、V5 working tree uncommitted、7 樣本未做、L4/imbalance 部分 V5 略遜）；(6) §10 明確「InterSubMod 本 repo 無 C++ 改動」；(7) §13 結論重寫；附錄 A 修正檔名（01_code_diff_analysis.md 等）| 用戶選 A 嚴格 fact-check 全修 | structured-tech-report skill (post fact-check) |
| 2026-04-29 03:35 | ✅ §4.2 「62% / Cohen's d=-1.20」核對完成 — 數字均正確；新增 §4.2.1 LOH 兩層精確化（ISM HP_Ratio LOH vs LOH.bed region-level LOH）+ 補 94.6% somatic→HP1、6,485 非 LOH 平衡位點 100% 出問題、kappa=0.670 等數字 | 待補項 #1 收斂 | structured-tech-report skill |
| **2026-04-29 04:00** | **概念性修正：F1 變動因果鏈澄清** — 用戶質疑「V5 沒改 phasing 與 caller，F1 為何變？」；對齊 PI 報告 4 §3.4 全 F1 對比表發現：(1) ClairS-TO raw F1=0.7166 對所有版本完全相同（V5 不改 caller）；(2) F1 變動 100% 來自 ISM SuggestFilter；(3) V5 vs Baseline F1 = -0.0003（**噪音級，非真實改善**）；(4) 之前報告把「F1 0.7117→0.7154」呈現為改進是誤導敘事 — 0.7117 與 0.7157 為同實驗不同 cut；(5) F1 不能衡量 V5 品質（PI 報告 4 第 268 行原文）。新增 §4.2.2「為何 V5『沒改 phasing 與 caller，F1 卻有微幅變動』」+ 流程圖 + 5 列證據表；§0 TL;DR 重寫；§8 表加 Raw F1 列；§13 結論加 F1 釐清句 | 用戶概念質疑（Hard Gate level）| structured-tech-report skill |
| **2026-04-29 04:30** | **證據鏈補強：3 層獨立證據 + 跨樣本一致性 + ISM 影響 3-tier 分類** — 用戶質疑「somatic ALT reads 集中於子族系一條 hap 是否有量化分析與圖像驗證」；補 4 處：(A) §5.A.1 加「3 層獨立證據鏈」表（理論 phasing graph edge weight + 全基因組 17.3:1 / 614K:35.5K + 個別位點 SP1/2/3 113:0/109:1/108:0）+ 連結至 v5_audit_suite/10_somatic_bias_explanation.md（含 fig01d 6-panel + 3 IGV 真截圖）；(B) §4.3 改寫為 ISM 特徵 3-tier 分類（🔴 38% / 🟡 7% / 🟢 55%）；(C) §4.4 新增跨樣本 7/7 一致性 + Simpson's paradox r=-0.964 + 同位點 r=0.001 + Cohen's d=-1.20；(D) §5.E 從 3 條擴展為 5 條獨立路徑（PON-only VCF / 全基因組統計 / IGV 真截圖 / sanity check / 程式碼最小必要）。**意義**：之前 §5.A L3 為定性論述，補強後變為「3 層獨立證據鏈」可被外部驗證；ISM 3-tier 分類釐清「self-phasing 並非汙染所有 ISM 特徵 — 55% 完全不受影響」 | 用戶因果鏈質疑 + 圖像驗證需求 | structured-tech-report skill |
| **2026-04-29 04:50** | **圖片 inline 渲染 + 表格限制說明** — 用戶要求「圖片用 `![](相對路徑)` 顯示，包含表格限制說明」；補：(A) §5.A.1 表格上方加 ⚠ 圖片渲染說明（表格儲存格寬度限制 + GitHub/VSCode/Obsidian 渲染不一致 + 慣例：表格列路徑作 source citation、獨立區塊用 `![](...)` 完整渲染）；(B) §5.A.1 表格之後新增 §5.A.1.1（理論層 fig2/3）+ §5.A.1.2（全基因組層 fig01d 6-panel + fig4_evidence_summary）+ §5.A.1.3（個別位點 SP1/2/3 IGV 6-BAM 並列）共 7 張 inline 圖（Figures 1-7）；(C) §4.3 表格之上加 ISM 影響 3-tier 總覽圖（Figure 0）；(D) §5.E 表格加「對應視覺」欄連結至已 inline 的 Figures。8 張圖路徑均使用報告位置的相對路徑 `../../../...`，已 ls 驗證實體檔案存在（fig01d 275KB、3 IGV 各 150-170KB、fig2/3/4 130-180KB、fig 03 136KB） | 用戶圖片顯示要求 | structured-tech-report skill |
| ⚠ 待補 | F1 commit V5 working tree 後回填 commit hash | audit suite 後續行動 #1 | TBD |
| ⚠ 待補 | F3 7 樣本擴展結果回填 | audit suite 後續行動 #3 | TBD |
