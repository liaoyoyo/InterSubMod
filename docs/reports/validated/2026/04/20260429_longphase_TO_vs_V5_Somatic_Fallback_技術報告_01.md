<!--
build_date: 2026-04-29 00:55
revised:
  - 2026-04-29 03:15 (fact-check vs v5_audit_suite 6 子報告)
  - 2026-04-29 03:35 (核對 02_Self_Phasing根因 + 補 LOH 兩層精確化)
  - 2026-04-29 04:00 (F1 變動因果鏈澄清 vs PI 報告 4 §3.4)
  - 2026-04-29 04:30 (3 層獨立證據鏈 + ISM 3-tier + 跨樣本 CV)
  - 2026-04-29 04:50 (8 inline 圖 + 表格限制說明)
  - 2026-04-29 05:00 (重寫對齊 20260428 PI 審查風格 + 13 段 → 11 節 + 主讀者切 PI + 範圍縮 self-phasing+V5)
agent: structured-tech-report skill (post 6 輪 fact-check + 對齊 20260428 主軸重寫)
status: validated
report_class: PI-audit (重寫前為 methodology-improvement)
audience: PI / 決策層為主（工程細節保留為 §5 工程審核段）
parent_audit: InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md
inputs:
  - InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md  # 主軸範式
  - InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md
  - InterSubMod/docs/reports/pi_reports/2026/04/20260424_V5_vs_Baseline_complete_comparison_01.md  # PI 報告 4 全基因組 V5=90.5%
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp  # V5 patch 實際位置（非 InterSubMod repo）
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.h
  - /big7_disk/liaoyoyo2001/longphase-to-mod/Util.h
  - InterSubMod/src/core/ReadParser.cpp  # ISM 下游 HP tag demotion
outputs:
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
related_memory:
  - memory/project_v5_somatic_fallback_verification.md
  - memory/project_v3_fixed_haplotag_verification.md
  - memory/project_self_phasing_causal_chain_confirmed.md
  - memory/project_pon_only_phasing_verification.md
verdict: VALIDATED — V5 為 production；6 大判決齊備；caveat 必須保留；HPFineNGroups marker / Thread D 屬另案
last_verified: 2026-04-29
report_template: structured-tech-report v1.0 → 對齊 20260428 PI audit 11 節
-->

# longphase-to-mod V5 Somatic Fallback Haplotag — Self-Phasing 審核報告（PI）

> ⚠️ **ERRATUM 2026-05-09 / 修訂 2026-05-10**：本報告 5 處表述已修訂，**主結論不撤回**。完整 5 條 patch + before/after 對照詳見 [`InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md`](../../05/20260509_PI_Report_4_29_Errata_01.md)：
> - **E1** §3.3.3 chr19 SP1/2/3 為「可重現案例」（占 priority bug 全基因組 **2.16%**，rank 19，非主要 hotspot 分佈）
> - **E2** §5.2 V5 working tree **已 commit**（`d0bcd8c` + `938f0df` 4-30 完成；V5 chain 為 **5 commits** 不是 4）
> - **E3** §5.2 priority bug 機制證據升級至「個案 + 統計 + 機制」三重佐證（read-level 全基因組 **34,855 victims** 鐵證，V3F+V5 修正率 **100%**）
> - **E4** §6.4/§6.5 V5 數值為 **Pass 1 only** 結果（4-12 BAM ploidy bug → purity=0 → Pass 2 從未觸發），主要功勞 **V3F + Layer 1.5**；Pass 2 second round 二次效益尚未獨立量化
> - **E5（5/10 加）** §5.2 V5 Layer 1.5 在 **germline-absent 區域與 baseline 4.19:1 偏 HP1 完全相同**（priority bug 的 feature 化而非修補）；V3F 標 hp=33 保守處理反而更穩健（5/9 paired cross-ref Step D 量化證明）

## 一句結論

**Self-phasing 在 baseline LongPhase-TO 是真實 tag-level artifact（94.6% somatic ALT reads 集中於 HP1，HP1:HP2 = 17.3:1）；修補在 longphase-to-mod fork 內以 4 commit 漸進完成（`8b8c1fd` PON-only flag → `41ff147` getVote 兩層 → `380e8d2` INDEL guard → V5 working tree Layer 1.5），總 +68/-36 行集中於 3 函式、`HaplotagProcess.h:66-68` 介面契約零變動；V5 通過 4 項硬性 sanity check 15/15 PASS、clean PS paired GT concordance +13.3 pp（全基因組 PI 報告 4 +8.3 pp），可作為新 haplotag baseline。但 V5 working tree 未 commit、Confidence threshold 0.6 未直接驗證、7 樣本擴展未做、cnLOH 雙親同源區仍 open、SEQC2 calling F1 完全持平（Raw 0.7166 對所有版本相同；V5 vs Baseline+ISM = -0.0003 噪音，F1 不衡量 tag 品質）。InterSubMod 在這條修補鏈是下游消費者而非實作者，本 repo 無 C++ 改動。**

本次審核以 v5_audit_suite 6 份子報告 + PI 報告 4 全基因組對比 + research_landscape/02 為主要證據，僅做證據審查不做大型重跑。

---

## 1. 審核範圍與判決

| 問題 | 判決 | 重點 |
|------|------|------|
| Self-phasing 問題是否真實 | **已證明成立** | baseline 17.3:1 somatic HP1:HP2 bias；TO HP_Ratio 與 paired r=0.001（n=288K）；3 層獨立證據鏈（理論／全基因組／個別位點 IGV）|
| LOH.bed 是否被 self-phasing 污染 | **目前判斷：否** | LOH.bed baseline vs PON-only Jaccard=1.0000；62% 指 ISM HP_Ratio LOH（BAM HP tag 路徑）非 LOH.bed |
| HP 偏斜的 tag-level 修補是否充分 | **充分（HCC1395 5kHz）** | clean PS paired GT V5=88.2% / BL=74.9%（+13.3pp）；4 項硬性 sanity 15/15 PASS、0 violation；HP:i:33 −54%；AMB% 17.5→8.0% |
| V5 是否改善 SEQC2 calling F1 | **否（噪音範圍）** | ClairS-TO Raw F1=0.7166 對所有版本完全相同（V5 不改 caller）；V5+ISM=0.7154 vs Baseline+ISM=0.7157 = **-0.0003 噪音**；F1 不衡量 tag 品質 |
| InterSubMod 本 repo 是否改 C++ | **否** | 修補在 longphase-to-mod 獨立 fork；ISM 為下游消費者；介面契約零變動；ISM F1=0.0124 對 TO germline FP 仍無修復力 |
| V5 working tree 是否可直接採用 | **可採用，但 caveat 必須保留** | sanity 全 pass + concordance 顯著改善；但 V5 = `380e8d2` + 兩塊未 commit working-tree 修改；Confidence threshold 0.6 未直接驗證；7 樣本擴展未做 |

> **範圍嚴格排除**：Thread D LOH-constrained phasing、HPFineNGroups marker grade 評級、cross-sample F1 ablation 屬另案（見 §8 caveat）。

---

## 2. 必要背景與術語邊界

### 2.1 三層資料不可混用

self-phasing 議題涉及三層獨立資料路徑，**不可在 metric 之間互相外推**：

| 層級 | 檔案 / 欄位 | 正確角色 | 常見誤用 |
|------|------------|---------|---------|
| **Caller** | ClairS-TO VCF `FILTER / AF / GQ / DP` | 決定 PASS candidate 與 caller-level F1 | 把 LongPhase phase 當作 re-caller |
| **Phasing** | phased VCF `GT / PS / GT2 / GT3`、`LOH.bed` | 決定 phase block、sub-genotype、region-level LOH | 把 BAM HP tag skew 寫成 LOH.bed |
| **Haplotag** | BAM `HP:i:1/2/11/21/33`、`PS` | read-level haplotype assignment | 把 `HP:i:21` 直接等同當前位點 ALT |

> 來源：`InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md` §2.1。

### 2.2 HP tag 編碼表（V5 final）

| BAM HP:i | ISM ReadParser 內部值 | 語意 | 觸發條件 |
|---------|---------------------|------|---------|
| `HP:i:0` 或無 tag | `"0"`（unphased）| 無方向證據 | germlineResult=0 且 somaticTotal=0 |
| `HP:i:1` | `"1"` | 純 germline HP1 | germlineHP1>0；無 somatic |
| `HP:i:2` | `"2"` | 純 germline HP2 | germlineHP2>0；無 somatic |
| `HP:i:11` | `"1-1"` | germline=HP1 + 至少一 somatic vote | germlineResult=1 且 somaticTotal>0 |
| `HP:i:21` | `"2-1"` | germline=HP2 + 至少一 somatic vote | germlineResult=2 且 somaticTotal>0 |
| `HP:i:33` | `"3"` | 只有 HAPLOTYPE3 somatic、無 HP1_1/HP2_1 方向（V5 後降 54%） | germlineResult=0 且 somaticTotal>0 |

V5 Layer 1.5 新增：當 germline=0 但 somatic HP1_1/HP2_1 有方向，從 `33` 升級為 `11`/`21`（保留 phased 方向資訊）。

### 2.3 術語表

| 術語 | 定義 | 邊界 |
|------|------|------|
| **Self-phasing** | TO 模式下 somatic / germline-like variants 進入 phasing graph 互相連結，把 read 拉向同一 HP | 屬 BAM HP tag layer；不影響 LOH.bed |
| **V5 Somatic Fallback** | longphase-to-mod fork 在 `HaplotagProcess.cpp::getVote()` 加 Layer 1.5 + `countSNPHaplotype()` alt guard | 屬 haplotag layer；不改 phasing graph 也不改 caller |
| **longphase-to-mod** | longphase 的本地 fork @ `/big7_disk/liaoyoyo2001/longphase-to-mod/`（獨立 git repo）| 4 commit 演進在此；非 InterSubMod 內 |
| **InterSubMod (ISM)** | 本 repo C++ 工具，read-level epigenetic characterization | 下游消費 longphase-to-mod 產出的 BAM；本身不修補 phasing/haplotag |
| **ISM HP_Ratio LOH** | ISM 計算：HP_Ratio<0.1 或 >0.9 | 受 self-phasing 嚴重污染（62% artifact）|
| **LOH.bed region-level LOH** | LongPhase phase 階段 VCF allele depth 推算 | **不**受 self-phasing 影響（PON-only Jaccard=1.0）|

---

## 3. Self-phasing 完整問題鏈與 3 層證據

### 3.1 問題來源

baseline LongPhase-TO 在 tumor-only 條件下沒有 matched normal；PON 可降低 germline leak，但不能完整排除所有 germline-like het。當 somatic 或 germline-like variants 進入 phasing / haplotag 決策時，read-level HP assignment 會被該 variant 自己與附近 variants 的共現 pattern 牽引，產生**自參考**（self-referential）的 phasing scaffold。

關鍵程式碼路徑（PhasingProcess.cpp:154-157）：

```cpp
if(params.ponOnlyPhasing){            // baseline: false → 跳過; V2b 起: true
    vGraph->convertNonGermlineToSomatic();   // 把非 PON-germline 標為 somatic
}
```

baseline 預設 `ponOnlyPhasing=false`，phasing graph 用 germline + somatic + unknown 混合 anchor → somatic 反客為主。

### 3.2 LOH 兩層精確化（**重要釐清**）

`02_Self_Phasing根因.md` 第 111+264 行明確區分**兩個 LOH 層次**，此區分是後續所有結論的基礎：

| LOH 層次 | 計算路徑 | self-phasing 影響 | 量化 |
|---------|---------|------------------|------|
| **ISM HP_Ratio LOH** | BAM HP tag → ISM HP_Ratio<0.1 or >0.9 | **嚴重**：62% 是 artifact，AF 0.1-0.8 近 100%；TO TP 86.5% 在 paired 下完全平衡（HP_Ratio 0.4-0.6） | **62%** ISM-level LOH artifact |
| **LOH.bed region-level LOH** | VCF allele depth → LongPhase region detection | **零**：PON-only Jaccard=1.0 完全不受影響 | LOH.bed 不變 |

兩套 LOH 系統使用不同定義（**kappa=0.670** 不完全一致）；self-phasing 的因果影響位於「haplotag → ISM HP_Ratio」這條路徑上，**非** LongPhase region detection。

### 3.3 3 層獨立證據鏈

對應 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` Section 5：

> **⚠ 圖片渲染說明**：以下表格的「視覺證據」欄為**檔案路徑**（用於追蹤來源）；**實際 inline 渲染**請見表格下方 §3.3.1 ~ §3.3.3 三個獨立圖片區塊。
> Markdown 表格儲存格內若塞 `![alt](path)` 會因儲存格寬度被壓縮（GitHub / VSCode / Obsidian 渲染不一致），無法清楚閱讀；本報告慣例：**表格列路徑作 source citation，獨立區塊用 `![](...)` 完整渲染**。

| 層 | 量化指標 | 數值 | 視覺證據（路徑）|
|---|--------|------|--------|
| **理論層** | Phasing graph edge weight：`weight(A,B) = Σ_reads I(read 帶 A.alt) × I(read 帶 B.alt)` | 同 clone somatic 在同 reads **共現 ≈ 100%**（共享 sub-population）；germline het 隨機分散 ≈ 50% → somatic-somatic edges 比 germline 權重更高 → somatic 反客為主定義 scaffold | `InterSubMod/docs/reports/pi_reports/2026/04/figures/fig2_self_phasing_concept.png`、`fig3_af03_walkthrough.png` |
| **全基因組層** | HCC1395 baseline TO HP1:HP2 reads | **614,000 : 35,500 = 17.3 : 1**（**94.6%** 集中於 HP1，跨 23 染色體）vs 隨機預期 ~1:1 → 生物學上不可能，必為 artifact | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01d_somatic_bias_explanation.png` Panel D；`fig4_evidence_summary.png` |
| **個別位點層**（IGV 6-BAM 並列） | 3 個 SP-extreme 位點 baseline HP2:HP1 ratio | SP1 chr19:17565944 = **113:0**；SP2 chr19:12452332 = **109:1**；SP3 chr19:12467180 = **108:0** — baseline 與 paired 方向**完全相反**；V5 修正後與 paired 一致（3/3） | `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP{1,2,3}_*.png` |

**3 層獨立但結論一致** — self-phasing artifact 真實存在；17.3:1 為 average，個別位點可達完全失衡（113:0），average 已被「不那麼極端的位點」稀釋。

#### 3.3.1 理論層視覺：phasing graph self-phasing 機制

![Self-Phasing 概念流程：同 clone somatic 共享 reads → phasing graph 強連結 → 全部塞同一 phase block](../../../pi_reports/2026/04/figures/fig2_self_phasing_concept.png)

*Figure 1 — Self-phasing 概念流程（PI 報告 1 Section 2）。同一 sub-clone 的 somatic variants 共享 sub-population reads，long read 跨多個 somatic variants 都帶 ALT，phasing graph 看到強連結 → 全部塞同一 phase block (HP1)。*

![AF=0.3 走例：30% ALT reads 在 paired 模式隨機分到 HP1/HP2，但 TO 模式因 self-phasing 全部指向 HP1](../../../pi_reports/2026/04/figures/fig3_af03_walkthrough.png)

*Figure 2 — AF=0.3 具體走例（PI 報告 1 Section 2.2）。同一個 somatic variant 在 paired 模式 HP_Ratio ≈ 0.5（隨機），在 TO 模式 HP_Ratio → 0.94（偏 HP1）。*

#### 3.3.2 全基因組層視覺：HCC1395 17.3:1 量化證據（6-panel 概念圖）

![Somatic Bias 17.3:1 — 6-panel 概念與實證視覺化（A: 預期分布 / B: 實測 self-phasing / C: 機制推導 / D: 全基因組 614K vs 35.5K / E: 個別極端位點 / F: V5 修復後）](../../../pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01d_somatic_bias_explanation.png)

*Figure 3 — `v5_audit_suite/10_somatic_bias_explanation.md` Section 2 主圖。Panel D 顯示 HP1=614K vs HP2=35.5K bar chart；Panel E 列 SP1/SP2/SP3 的 113:0 / 109:1 / 108:0 個別失衡；Panel F 顯示 V5 修復後 ~1:1。*

![量化證據總表：Self-Phasing 影響 7 項指標（17.3:1、94.6%、62% LOH artifact、Cohen's d=-1.20、7/7 一致）](../../../pi_reports/2026/04/figures/fig4_evidence_summary.png)

*Figure 4 — 量化證據總表（PI 報告 1 Section 3）。整合 Somatic HP bias、ISM HP_Ratio LOH artifact 比例、跨樣本一致性等 7 項證據。*

#### 3.3.3 個別位點層視覺：3 個 SP-extreme 位點 IGV 6-BAM 並列截圖

> 每張截圖由上至下 6 個 panel：**baseline / V2b / V3-Fixed / V5 / Paired tumor / Paired normal**。觀察重點 — baseline 的 HP 主導方向與 paired tumor **完全相反**（self-phasing artifact 在單一位點的實證），V5 修正後與 paired 一致。

![SP1 chr19:17565944 — baseline 113:0 全在 HP1，paired 與 V5 都是 HP2 主導](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png)

*Figure 5 — SP1（chr19:17565944）。baseline panel reads 全部集中於 HP1+HP1-1（粉紅+淡綠群）；V5 / V2b / V3-Fixed reads 整體翻轉至 HP2+HP2-1；Paired tumor 確認 HP2 為真實方向。*

![SP2 chr19:12452332 — baseline 109:1，V5 翻轉至 HP2 主導對齊 paired](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP2_chr19_12452332.png)

*Figure 6 — SP2（chr19:12452332）。baseline HP1+HP1-1 集中 109 reads，HP2 stack 僅 1 read；V5 方向翻轉至 HP2+HP2-1，與 paired tumor 一致。*

![SP3 chr19:12467180 — baseline 108:0 與 SP1 相同模式](../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP3_chr19_12467180.png)

*Figure 7 — SP3（chr19:12467180）。與 SP1 同模式 — baseline HP1 主導 → V5 HP2 主導，HP orientation 整體翻轉，V5 與 paired 一致。*

> **完整 6-panel 概念圖細節 + 3 IGV 截圖視覺解讀**詳見 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` Section 2-3。

### 3.4 跨樣本一致性與 Simpson's Paradox 釐清（7/7 樣本）

對應 `InterSubMod/docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md` CV-1～CV-5：

| 跨樣本指標 | 數值 | 解讀 |
|----------|------|------|
| **CV-2 方向一致性** | **7/7 樣本全部觀察到相同方向的 self-phasing 效應** | 排除樣本特異性可能 |
| 同位點 HP_Ratio 跨模式 r | **r = 0.001**（n=288K pairs）| TO 與 paired 在同一變異上的 HP 判定**完全不相關** |
| TO-only LOH 在 paired 下完全平衡 | **86.5%**（HP_Ratio 0.4-0.6）| TO 看到的 LOH 在 paired 下大多是平衡的 |
| Cohen's d (HP_Ratio 差異) | **-1.20** | 巨大效應量（一般 \|d\|≥0.8 即視為大）|
| CV-1 Simpson's paradox | **r = -0.964** 跨樣本 imbalance vs self-phasing fraction | 強負相關，**非**預期的正相關 |

**Simpson's Paradox 機制釐清**：結構性 LOH 比例高的樣本（HCC1937、HCC1395）有高整體 imbalance，但 self-phasing fraction 反而較低（被結構性 LOH 稀釋）；**不否定 self-phasing 機制** — 是「加法關係非替代關係」。

---

## 4. 錯誤分析與處理狀態

### 4.1 舊說法 / 現狀 / 處理充分度 / 建議寫法

對應 `20260428_Self_Phasing_Baseline_V5_Audit_01.md` §5：

| 舊說法 / 風險 | 現在狀態 | 處理是否充分 | 建議寫法 |
|---------------|----------|--------------|----------|
| Self-phasing 造成 TO LOH（含 LOH.bed）| **過寬，已修正** | 充分 | 拆兩句：BAM HP_Ratio LOH 受影響（62% artifact）；LOH.bed Jaccard=1.0 不變 |
| V5 全面優於 baseline | **過寬** | 已在 V5 audit suite 修正 | V5 在 sanity 與 clean PS paired concordance 勝出（+13.3pp）；problem PS / weak directional 位點仍有限 |
| LongPhase phase/tag 改善 caller F1 | **不支持** | 已由 V5 suite §3.4 釐清 | tag 改變 HP 品質；caller F1 由 ClairS-TO PASS set 決定（Raw 0.7166 不變）|
| `HP:i:21` 必然是當前位點 ALT | **錯誤，已修正** | 已列入 code issue inventory | HP tag 是 read-level phasing state；需用 per-site ALT/REF 檢查 |
| InterSubMod 內部修了 V5 | **錯誤，2026-04-29 fact-check 修正** | 充分 | V5 patch 在 longphase-to-mod fork 內；ISM 為下游消費者，本 repo 無 C++ 改動 |
| F1 0.7117 → 0.7154 = V5 改善 | **誤導敘事** | 充分（PI 報告 4 §3.4）| Raw F1=0.7166 對所有版本相同；V5 vs Baseline+ISM = -0.0003 噪音；F1 不衡量 tag 品質 |
| HPFineNGroups 是 methylation bimodality | **錯誤** | 已由 Thread D 修正（屬另案）| phasing bucket occupancy / LOH-constrained phasing signature；需 Phase 2B 重驗 |

### 4.2 ISM 特徵 3-tier 影響分類（self-phasing 直接 / 間接 / 不受影響）

對齊 `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md` + `03_ISM分析價值界定.md` 的特徵分類。**self-phasing 並非汙染所有 ISM 特徵 — 55% 特徵完全不受影響**：

![Self-Phasing 對 ISM 特徵 3-tier 影響分級總覽（嚴重 38% / 中度 7% / 不受影響 55%）](../../../research_landscape/figures/03_self_phasing_impact.png)

*Figure 0 — Self-Phasing 影響量化總覽。整合 Somatic HP bias 17.3:1、ISM HP_Ratio LOH 62% artifact、跨樣本 7/7 一致性等指標於單一視圖。*

| 影響等級 | 特徵數（佔比）| 代表特徵 | TO 結果處理 |
|---------|------------|---------|-----------|
| 🔴 **嚴重影響**（直接依賴 HP）| **29 個（38%）** | HP_Ratio → 假 LOH；Potential_LOH → 62% artifact；HPMergedDelta/Sig → 方向反轉；hp_assign_rate；effective_hp_reads；HPFineNGroups | **必重跑** V5 BAM；舊結論需加註版本 |
| 🟡 **中度影響**（間接污染）| **14 個（7%）** | QualityScore → AUC 0.497（已移除 LOH penalty）；GlobalP；CramersV；VerificationClass | 重評；多數影響微弱或已程式碼移除 |
| 🟢 **不受影響**（無 HP 依賴）| **42 個（55%）** | PairwiseMean/MedianDist；AlleleDelta / AlleleP；Caller 特徵（AF/GQ/DP/SB）；甲基化矩陣；CpG 座標、region_methyl_mean | **結論穩固，不需重測** |

**對既有結論的具體影響**：
- 🟢 paired-pure delta F1=+0.0112 等與 HP tag 無關的結論不受影響
- 🔴 歷史 archive 數據 V5 之前的 HP-依賴 ISM 結果需重跑或加註版本（`haplotag_version` 欄位 — 後續 §9 F5）
- 🔴 HPFineNGroups subclone marker 需在 V5 + flag=on 重驗（屬另案，見 §8 caveat）

---

## 5. GT / Phasing / Tagging 細節審核

### 5.1 VCF GT 與 BAM HP 角色不同

對應 `v5_audit_suite/12_gt_distribution_audit.md`：baseline vs V2b/V5 的 PASS somatic GT class 幾乎一致：

| GT class | Baseline % | V2b % | Δpp |
|----------|-----------:|------:|----:|
| Germline_Het | 13.09 | 12.91 | -0.18 |
| Germline_Hom_or_LOH | 7.48 | 7.48 | 0.00 |
| Somatic_NoLOH | 42.99 | 42.99 | 0.00 |
| Somatic_in_LOH | 24.19 | 24.19 | 0.00 |
| Unphased | 12.26 | 12.44 | +0.18 |

→ 17.3:1 HP skew **不能寫成「PASS somatic GT 被大量改判」**；它主要發生在 **haplotag / read assignment 層**，不是 caller 層。

### 5.2 V5 tag 修補機制 — 4 commit 演進

對應 `v5_audit_suite/01_code_diff_analysis.md`：

| 版本 | Git ref | 修改內容 | 程式碼變動 | HP tag 影響 |
|------|---------|---------|----------|-----------|
| **baseline** | parent of `8b8c1fd` | 原始 LongPhase getVote（priority bug + enum mismatch）| n/a | 99.9% reads 偏 HP21（PON-only 啟用後暴露） |
| **V2b** | `8b8c1fd` | 加 `--pon-only-phasing` flag；`HaplotagProcess.cpp` 不變 | `Phasing.cpp` +9/-2、`PhasingGraph.cpp` +34/-0、`PhasingProcess.cpp` +25/-3 | 暴露 getVote priority bug |
| **V3-Fixed** | `41ff147` | `getVote()` 重寫為兩層（germline first / somatic second）+ 修 enum→int literal bug | `HaplotagProcess.cpp:506-541` +36/-25 | HP21 偏移修正；HP:i:33 開始正確出現（HCC1395=6,793）|
| **INDEL guard** | `380e8d2` | `countINDELHaplotype()` 加 `HAPLOTYPE_UNDEFINED` guard | `HaplotagProcess.cpp:497-510` +8/-4 | 消除 PON-only somatic INDEL 站點 UB |
| **V5（uncommitted working tree）** | working tree | `getVote()` Layer 1.5 somatic fallback + `countSNPHaplotype()` 對稱 alt guard | Layer 1.5：`HaplotagProcess.cpp:512-563` +15/-1；alt guard：`HaplotagProcess.cpp:489-494` +9/-6 | AMB% 17.5→8.0%；HP:i:33 −54%（239,679 → 110,197）|

**總計 +68 / -36 行**，集中於 3 函式（`getVote()` / `countSNPHaplotype()` / `countINDELHaplotype()`）；`HaplotagProcess.h:66-68` 三函數簽章從 baseline 到 V5 **一字未變**；無新 enum / 新 member / 新 std 容器 / 新 logging。

### 5.3 SEQC2 calling F1 變動因果鏈釐清（**關鍵概念**）

對應 PI 報告 4 §3.4 第 252-268 行完整對比表：

```
ClairS-TO (caller) ──→ VCF（FILTER: PASS/NonSomatic/LowQual）
                           │  Raw F1 = 0.7166（所有版本完全相同 ← V5 不改此處）
                           ▼
                       longphase-to-mod (V5 改這層)
                           │  V5 改 BAM 的 HP:i tag，不改 VCF
                           ▼
                       InterSubMod ISM ──→ region 特徵（HP_Ratio…）
                           │  HP tag 變了 → region 特徵變了
                           ▼
                       ISM SuggestFilter（InterSubMod/src/core/RegionProcessor.cpp:1120, 1269）
                           │  對每個 variant 計算「該不該標 LOW_QUAL」
                           │  ←── F1 變動的唯一來源
                           ▼
                       套濾後最終 F1：Baseline+ISM 0.7157 / V3F+ISM 0.7154 / V5+ISM 0.7154
                                       V5 vs Baseline = -0.0003（噪音）
```

| 證據 | 數字 | 含義 |
|------|------|------|
| ClairS-TO Raw F1（所有版本）| **0.7166** | V5 不改 caller |
| Baseline+ISM SuggestFilter | 0.7157 | ISM 過濾本身略微負面 -0.0009 |
| V3F+ISM | 0.7154 | -0.0012 |
| V5+ISM | **0.7154**（vs Baseline+ISM = **-0.0003 噪音**）| V5 vs V3F = +0.0001 噪音 |
| ISM SF Precision（V5）| 113 TP / 74 FP；FP catch rate **0.63%** | ISM 過濾無效於 ClairS-TO 的 germline FP |

**為何 F1 不能衡量 V5 品質**（PI 報告 4 第 268 行原文）：
- ClairS-TO 的 FP **主要是 germline variants**，非 somatic；ISM 甲基化分析設計用於分 subclone 結構，不是分 germline vs somatic
- V5 修正不改變此根本限制；**F1 持平不代表 tags「沒有更好」，F1 變化也不代表 tags「更好或更差」**

V5 真實價值的正確衡量在 **read-level tag 品質**：clean PS paired GT concordance +13.3 pp（HCC1395）/ +8.3 pp（PI 報告 4 全基因組）、AMB% 17.5→8.0%、HP:i:33 −54%、4 項 sanity 15/15 PASS。

---

## 6. baseline vs V5 量化對照

### 6.1 ISM aggregate（HCC1395 5kHz）

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/version_summary.tsv`

| Version | n_regions | TP | FP | TP_rate | Potential_LOH% | HP_Ratio median | HPFineNGroups=2% | HPFineNGroups>=3% |
|---------|----------:|---:|---:|--------:|---------------:|----------------:|------------------:|-------------------:|
| baseline_old | 40,213 | 28,383 | 11,830 | 0.706 | 58.7 | 0.788 | 51.2 | 27.9 |
| v5_new | 40,096 | 28,495 | 11,601 | 0.711 | 62.2 | **0.574** | 49.7 | 29.8 |

→ TP_rate 差異很小（不能宣稱 caller F1 改善）；HP_Ratio median 從 0.788 拉回 **0.574**，符合 tag bias 被修正方向。

### 6.2 HP_Ratio AUC

| Version | Side | n | TP_rate | HP_Ratio_AUC |
|---------|------|--:|--------:|-------------:|
| baseline_old | All | 40,213 | 0.706 | 0.531 |
| v5_new | All | 40,096 | 0.711 | 0.526 |
| baseline_old | Inner | 23,620 | 0.729 | 0.523 |
| v5_new | Inner | 24,949 | 0.730 | 0.525 |

→ V5 修正 tag distribution，但 HP_Ratio 本身仍**不是強 filter feature**；AUC 接近隨機。

### 6.3 methylation feature

| Feature | baseline AUC_All | v5 AUC_All | 判斷 |
|---------|-----------------:|-----------:|------|
| HPMergedDelta | 0.517 | 0.514 | 不變 |
| AlleleDelta | 0.541 | 0.543 | 不變 |
| HPFineF | 0.576 | 0.584 | 小幅變動，仍未成為可靠 filter |
| PairwiseMeanDist | 0.525 | 0.528 | 不變 |
| CramersV | 0.507 | 0.506 | 不變 |
| GlobalP | 0.533 | 0.530 | 不變 |

→ V5 沒有推翻 pure methylation feature 空間耗盡的結論。

### 6.4 全基因組 paired ground-truth concordance（PI 報告 4 §3.7）

| Cohort | Baseline | V5 | Δ |
|--------|:--------:|:---:|:----:|
| **全基因組 clean PS blocks** | **82.2%** | **90.5%** | **+8.3 pp** |
| 15-site Aggregate (pooled) | 72.20% | **78.85%** | **+6.65 pp** |
| 15-site Clean PS（11 sites）| 74.9% | **88.2%** | **+13.3 pp** |
| 15-site Problem PS（2 sites, SP1/SP2/SP3）| 48.5% | 52.0% | +3.5 pp |

### 6.5 全基因組 HP 結構改善

| 指標 | Baseline | V5 | 改善 |
|------|:--------:|:---:|:----:|
| HP1:HP2 ratio | 17.3:1 | ~1:1 | 消除 |
| 94.6% somatic→HP1 | 是 | ~50% 平衡 | 修正 |
| Phase block N50 | 4,061 | 8,109 | **+99.7%** |
| Phased rate | 54.9% | 78.5% | **+23.6 pp** |
| 執行時間 | 2,693s | 1,976s | 1.36× 快 |

---

## 7. InterSubMod ReadParser 與 longphase-to-mod 的關係

### 7.1 上下游架構（**澄清誤解：ISM 為下游消費者**）

```
ClairS-TO (外部 binary) ────→ VCF ────────────────────────┐
                                                          │
tumor BAM ───────────────────────────────────────────────→├─→ longphase-to-mod (本地 fork)
                                                          │     ├─ HaplotagProcess.cpp (4 commits, V5 修補在此)
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

**InterSubMod 在這條修補鏈是下游消費者而非實作者**：本 repo `src/core/` 不包含任何 `HaplotagProcess` / `getVote` / `Phasing*` 檔案；本 repo 在此修補週期**無 C++ 改動**。

### 7.2 ReadParser HP tag demotion 邏輯

`InterSubMod/src/core/ReadParser.cpp:120` 邏輯：

- raw `HP:i:1` → `"1"`
- raw `HP:i:2` → `"2"`
- raw `HP:i:11` → `"1-1"`
- raw `HP:i:21` → `"2-1"`
- raw `HP:i:33` → `"3"`
- `--germline-hp-only` 開啟時：`"1-1" / "2-1" / "3"` 被 demote 為 `"0"`，但 `hp_tag_raw` 保留供 audit

**重要邊界**：`--germline-hp-only` 是 **InterSubMod 下游 demotion**（屬另案，CONDITIONAL NEGATIVE on filter）；**V5 是 longphase-to-mod 上游 phase/tag 修補**；兩者不是同一個 fix。

### 7.3 測試覆蓋

- `InterSubMod/tests/test_read_parser.cpp:83`：HP tag parsing 與 demotion
- `InterSubMod/tests/test_global_local.cpp:404`：HPFine label mapping
- 已驗證：`./build/bin/run_tests --gtest_filter=ReadParserHPTagTest.*:GlobalTestTest.HPFine*:FullLabelTest.HPFineLabel_Mapping:SignificanceAnalyzerTest.HPFamily_GatingPassesWhenPureFails`，**17/17 pass**

### 7.4 ISM 對 TO germline FP 的根本限制

ISM F1 = 0.0124（vs Paired 0.0909），對 TO germline FP **仍無修復力**：
- ClairS-TO 的 FP 主要是 germline variants
- ISM 甲基化分析設計用於 subclone 結構區分（read-level epigenetic context），非 germline/somatic 區分
- V5 BAM HP tag 改善後，ISM SuggestFilter 對 F1 的影響仍只 -0.0003 ~ +0.0001 噪音級

---

## 8. 仍需保留的 caveat

| Caveat | 影響 | 建議 |
|--------|------|------|
| **R1 V5 working tree 未 commit** | V5 = `380e8d2` + Layer 1.5 + SNP alt guard 兩塊未 commit；可追溯性不足 | 切 2 獨立 commits：`feat(haplotag): Layer 1.5 somatic fallback in getVote()` + `fix(haplotag): guard countSNPHaplotype against UNDEFINED on alt path` |
| **R2 Confidence threshold 0.6 未直接驗證** | V5 設計證據仍有一格間接 | 加 V5 binary vote log 或用 IGV session 直接驗證 PS block 內 ALT/REF 投票 |
| **R3 7 樣本擴展未做** | 僅 HCC1395 5kHz 一樣本驗證 | 7 樣本 V5 BAM 全量重跑（HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829）|
| **R4 cnLOH 雙親同源 region** | 兩 hap 序列接近相同，V5 fallback 也無法區分 | 需 cnLOH-aware filter（CN 層 + germline-only methylation reference）|
| **R5 `HP:i:21` 不必然是當前位點 ALT** | HPFineNGroups bucket 可能混合 ALT/REF | 加 8-bucket derived metric 或 per-site ALT/REF audit |
| **R6 paired ground truth 仍由 phasing tool 產生** | paired concordance 不是絕對真值 | 用 clean/problem PS 分層；加 trio-phased 作 second ground truth |
| **R7 0.6 purity simulation 尚未執行** | V5 低 purity generalization 仍是推論 | 按 `20260427_purity06_simulation_plan_01.md` 另案執行 |
| **R8 範圍排除：Thread D LOH-constrained phasing** | 屬另案；TO 層論文主軸 pivot；與本報告 self-phasing+V5 主題分離 | 詳見 `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md` |
| **R9 範圍排除：HPFineNGroups marker grade B** | 屬另案；需 Phase 2B master + flag=on/off 7 樣本重驗 | 詳見 `memory/project_hpfinengroups_subclone_marker.md` |
| **R10 範圍排除：cross-sample F1 ablation** | 屬另案；V5 fallback 閾值未做 7 樣本 ablation；fallback 閾值取自 SEQC2 + 4 樣本經驗值 | 7 樣本擴展（R3）後做閾值 grid scan |

### 證據獨立性說明（本報告 5 條獨立路徑）

| # | 證據 | 數值／結論 | 對應視覺 |
|---|------|-----------|---------|
| 1 | PON-only VCF 層修正 | LOH.bed Jaccard=1.0、Somatic bias 17.3:1 → 消除、Phase block N50 +99.7% | Figure 3 Panel F |
| 2 | 全基因組 17.3:1 量化 | HCC1395 HP1=614K / HP2=35.5K（94.6% to HP1）；跨 23 染色體一致 | Figure 3 Panel D, Figure 4 |
| 3 | 個別位點 IGV 6-BAM 真截圖 | SP1=113:0 / SP2=109:1 / SP3=108:0 baseline 與 paired 方向相反；V5 一致 3/3 | **Figures 5, 6, 7** |
| 4 | V5 working tree sanity check | 4 項硬性檢查 15/15 PASS、0 violation | （fig06a/06b in audit suite）|
| 5 | 程式碼層最小必要 | +68/-36 行集中 3 函式、`HaplotagProcess.h:66-68` 介面契約零變動 | （fig01a/01b in audit suite）|

5 條獨立路徑：任一條失敗都不會推翻其他四條，5 條全通過構成穩定度 4-5/5 結論。

---

## 9. 最終建議口徑

### 9.1 可直接使用

1. **Self-phasing 是 baseline TO tag-level artifact**，會污染 HP-dependent ISM features（38% 嚴重 / 7% 中度 / 55% 不受影響）。
2. **LOH.bed 不因 self-phasing 修正而改變**；LOH.bed 與 BAM HP_Ratio 是不同路徑（kappa=0.670）。
3. **V5 是 longphase-to-mod fork 內的 4 commit 漸進修補**，總 +68/-36 行集中於 3 函式、介面契約零變動；**InterSubMod 本 repo 無 C++ 改動**，ISM 為下游消費者。
4. **V5 真實價值在 read-level tag 品質**：clean PS paired GT +13.3 pp（15-site）/ +8.3 pp（全基因組）；4 項 sanity 15/15 PASS；HP:i:33 −54%；AMB% 17.5→8.0%。
5. **V5 可作為新 haplotag baseline**，但 caveat 必須保留（V5 working tree uncommitted、Confidence threshold 0.6 未直接驗證、7 樣本擴展未做）。

### 9.2 不建議再使用

1. 不要把 SEQC2 calling F1 變動寫成「V5 改善 caller」（V5 不改 caller，Raw F1=0.7166 對所有版本相同；V5 vs Baseline+ISM = -0.0003 噪音）。
2. 不要把「62% LOH 消失」寫成 LOH.bed 消失（限指 ISM HP_Ratio LOH，BAM HP tag 路徑）。
3. 不要把 V5 修補寫成「InterSubMod 內部改動」（V5 在 longphase-to-mod 獨立 fork）。
4. 不要把 `HP:i:11/21/33` 直接等於當前位點 ALT（HP tag 是 read-level phasing state，需 per-site ALT/REF 檢查）。
5. 不要把單一 HCC1395 結論寫成跨樣本 production filter（7 樣本擴展未做；保留 HCC1395 5kHz 主驗證限定）。

### 9.3 後續行動（依優先序）

| ID | 動作 | 優先 |
|----|------|-----|
| F1 | commit V5 working tree 修改（切 2 獨立 commits）| **高** |
| F2 | 追加 Confidence threshold 0.6 投票 log 驗證 | 中 |
| F3 | 7 樣本 V5 BAM 全量重跑 | 中 |
| F4 | manifest.yaml 加 `haplotag_version: V3F / V5` 欄位 | 中 |
| F5 | cnLOH 區獨立評估方案（CN 層 + germline-only methylation reference）| 中 |
| F6 | trio-phased second ground truth | 低 |

---

## 10. 參考文件索引

| 類別 | 文件 |
|------|------|
| **主軸範式** | `InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md` |
| 根因分析 | `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md` |
| ISM 特徵分類 | `InterSubMod/docs/reports/research_landscape/03_ISM分析價值界定.md` |
| V5 audit 母索引 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md`（19 子報告）|
| V5 程式碼 diff | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` |
| V5 sanity check | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md` |
| V5 paired concordance | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md` |
| V5 synthesis | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md` |
| Somatic bias 17.3:1 整合 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` |
| GT distribution audit | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md` |
| Phase vs Tag 演算法 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md` |
| PI 報告 4 全基因組 | `InterSubMod/docs/reports/pi_reports/2026/04/20260424_V5_vs_Baseline_complete_comparison_01.md` |
| Causal chain 因果鏈 | `InterSubMod/docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md` |
| LongPhase 文獻 | `InterSubMod/docs/references/manual/20260402_longphase_to_phasing_quality_literature.md` |
| 程式碼（longphase-to-mod fork）| `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:484-563`、`HaplotagProcess.h:66-68`、`Util.h:20-25`、`PhasingProcess.cpp:154-157` |
| 程式碼（InterSubMod 下游）| `InterSubMod/src/core/ReadParser.cpp:120`（HP tag demotion）|
| MEMORY | `memory/project_v5_somatic_fallback_verification.md`、`project_v3_fixed_haplotag_verification.md`、`project_self_phasing_causal_chain_confirmed.md`、`project_pon_only_phasing_verification.md` |
| 範圍排除（另案）| Thread D：`20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`；HPFineNGroups marker：`memory/project_hpfinengroups_subclone_marker.md` |
| Knowledge Base | `/big8_disk/liaoyoyo2001/knowledge/03_file_formats/bam-format.md`、`vcf-longphase.md`、`05_tools/longphase-to.md`、`06_workflows/phasing-workflow.md` |

---

## 附錄 A：圖表清單

8 張 inline 圖編號表（按章節順序）：

| Figure | 章節 | 主題 | 路徑（相對於本報告）|
|--------|------|------|--------|
| **Figure 0** | §4.2 | Self-Phasing 影響 3-tier 總覽 | `../../../research_landscape/figures/03_self_phasing_impact.png` |
| **Figure 1** | §3.3.1 理論層 | Self-phasing 概念流程 | `../../../pi_reports/2026/04/figures/fig2_self_phasing_concept.png` |
| **Figure 2** | §3.3.1 理論層 | AF=0.3 具體走例 | `../../../pi_reports/2026/04/figures/fig3_af03_walkthrough.png` |
| **Figure 3** | §3.3.2 全基因組層 | 17.3:1 6-panel 概念圖（A-F）| `../../../pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01d_somatic_bias_explanation.png` |
| **Figure 4** | §3.3.2 全基因組層 | 量化證據總表 | `../../../pi_reports/2026/04/figures/fig4_evidence_summary.png` |
| **Figure 5** | §3.3.3 個別位點 | SP1 IGV 6-BAM 113:0 | `../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png` |
| **Figure 6** | §3.3.3 個別位點 | SP2 IGV 6-BAM 109:1 | `../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP2_chr19_12452332.png` |
| **Figure 7** | §3.3.3 個別位點 | SP3 IGV 6-BAM 108:0 | `../../../pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP3_chr19_12467180.png` |

---

## 附錄 B：變更歷史

| 日期 | 變更 | 觸發 | 作者 |
|-----|------|-----|-----|
| 2026-04-29 00:55 | 初稿（structured-tech-report skill 首個示範案例 13 段技術報告）| 用戶請求建立報告 skill | structured-tech-report v1.0 |
| 2026-04-29 02:30 | 用戶質疑 V5 歸屬（「應該是只 longphase-to 的 HaplotagProcess 吧」）| 用戶事實校正質疑 | 用戶 |
| 2026-04-29 03:15 | Fact-check 全面校正：V5 歸屬至 longphase-to-mod fork、4-commit 演進敘事、刪除捏造數字（Phase block N50、test_haplotag_v5.cpp、7 樣本 ≥0.78）| 用戶選 A 嚴格 fact-check 全修 | structured-tech-report skill |
| 2026-04-29 03:35 | §4.2 「62% / Cohen's d=-1.20」核對；新增 LOH 兩層精確化（ISM HP_Ratio LOH vs LOH.bed）+ 補 94.6% / 6,485 / kappa=0.670 | 用戶選 A 收斂待補項 #1 | structured-tech-report skill |
| 2026-04-29 04:00 | F1 變動因果鏈澄清（V5 不改 caller；Raw F1=0.7166 對所有版本相同；F1 100% 來自 ISM SuggestFilter；V5 vs Baseline = -0.0003 噪音）| 用戶概念質疑 | structured-tech-report skill |
| 2026-04-29 04:30 | 證據鏈補強：3 層獨立證據（理論／全基因組／個別位點 SP1/2/3）+ 跨樣本 7/7 + ISM 影響 3-tier 分類（38%/7%/55%）+ §5.E 5 條獨立路徑 | 用戶因果鏈與圖像驗證質疑 | structured-tech-report skill |
| 2026-04-29 04:50 | 8 inline 圖完整渲染（Figures 0-7）+ 表格圖片限制說明 | 用戶圖片顯示要求 | structured-tech-report skill |
| **2026-04-29 05:00** | **重寫對齊 20260428 PI 審查風格**：13 段教材 → 11 節審查；主讀者切 PI / 決策層；範圍縮 self-phasing+V5（Thread D / HPFineNGroups marker 移至 §8 caveat 另案）；保留全部 fact-check 數字 + 8 inline 圖 + 3 層證據鏈；§7.1 非工程版白話刪除（PI 不需要）；補強：4 commit 演進表 / 三層資料不可混用表 / HP tag 編碼表 / GT 分布表 / 完整 F1 對比 / 全基因組 vs 15-site concordance 對照 | 用戶 4 問答決策（plan mode 批准）| structured-tech-report skill (post 主軸對齊) |
| ⚠ 待補 | F1 commit V5 working tree 後回填 commit hash | audit suite 後續行動 #1 | TBD |
| ⚠ 待補 | F3 7 樣本擴展結果回填 | audit suite 後續行動 #3 | TBD |
