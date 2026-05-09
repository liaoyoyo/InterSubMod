<!--
build_date: 2026-05-09
agent: errata patch (從 5/8 整合報告 §9.2 + 5/9 V5 commit 狀態更新)
status: validated
report_class: errata-companion (修訂 4-29 PI 報告 4 處表述)
audience: PI / lab member / 自己未來
parent_report: InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
parent_synthesis: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
inputs:
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md (PI 報告原文)
  - InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md §9.2 (4 條 errata 清單)
  - InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md (34,855 全基因組鐵證)
  - longphase-to-mod commit log (d0bcd8c / 938f0df 已 commit)
outputs:
  - 本檔（獨立 erratum companion）
  - PI 報告頂部加一行 erratum banner 引用本檔
verdict: 4 條 errata，主要 PI 結論不撤回；補強為「機制 + 案例 + 統計」三重佐證；歸因從「V5 整體」精確化為「V3F + Layer 1.5 為主，Pass 2 二次效益尚未獨立量化」
last_verified: 2026-05-09
report_template: errata-companion v1.0
-->

# Errata for PI Report 2026-04-29 — Self-Phasing 整合報告（5/8）後續修訂

## 0. TL;DR

PI 報告 [`InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`](../2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md)（2026-04-29）**主要結論不撤回**（self-phasing 真實存在、V3F + V5 修正確立、V5 可作 production tag baseline）。本 erratum 修訂 4 處表述：

1. **E1** §3.3.3 chr19 SP1/2/3 解讀：從「主要 hotspot」→「可重現案例（chr19 占全基因組 priority bug 僅 2.16%）」
2. **E2** §5.2 V5 commit 狀態：「V5 working tree uncommitted」→ 已 commit (`d0bcd8c` 4-30 含 ploidy fix + bundled Layer 1.5 + countSNP guard；`938f0df` 4-30 含 threshold 0.95→0.9)
3. **E3** §5.2 priority bug 機制證據強度：升級「commit message + 3 IGV 截圖」→「+ 34,855 read-level victims 全基因組鐵證、V3F+V5 修正率 100%」
4. **E4** §6.4/§6.5 V5 數值歸因：精確化「V5 four-commit chain 整體效益」→「**主要 V3F + Layer 1.5（tagging layer）；Pass 2 second round 二次效益尚未獨立量化**」（PI 報告引用之 4-12 V5 BAM 實為 Pass 1 only，因 ploidy bug 未修正前 purity=0 → highPurity=false → Pass 2 從未觸發）

證據來源：5/8 整合報告 §7（V5 Provenance bonus）、§8（顛覆性發現）、§5（修補設計演進）。

---

## 1. E1 — §3.3.3 chr19 SP1/2/3 解讀降級

### 1.1 PI 報告 4-29 原文（line 146 + 170-184）

> **個別位點層**（IGV 6-BAM 並列）｜ 3 個 SP-extreme 位點 baseline HP2:HP1 ratio：SP1 chr19:17565944 = **113:0**；SP2 chr19:12452332 = **109:1**；SP3 chr19:12467180 = **108:0** — baseline 與 paired 方向**完全相反**；V5 修正後與 paired 一致（3/3）

§3.3.3 透過三 IGV 截圖呈現 chr19 三個 SP-extreme 位點為 self-phasing artifact 證據。**隱含敘事是 chr19 為主要 hotspot**。

### 1.2 5/8 整合報告 §8.1 顛覆性發現

T1.2-F1 全基因組 vote audit 顯示 priority bug 在 chr19 占比僅 **2.16%**（752 / 34,855 全基因組 victims，rank 19）：

| chr | victims | 占全基因組 | rank |
|---|---:|---:|---:|
| chr7 | 3,508 | 10.1% | 1 |
| chr2 | 2,792 | 8.0% | 2 |
| chr1 | 2,674 | 7.7% | 3 |
| chr16 | 2,584 | 7.4% | 4 |
| chr20 | 2,101 | 6.0% | 7 |
| **chr19** | **752** | **2.16%** | **19** |

### 1.3 修訂建議

**保留**：chr19 SP1/2/3 IGV 截圖、113:0 / 109:1 / 108:0 數字、baseline vs paired 方向相反、V5 修正 3/3 對齊 paired — 全部仍正確、仍是極佳教學案例。

**修訂**：
- 「個別位點層」段加註：「**chr19 SP1/2/3 是經 IGV 篩選之可重現極端案例**，全基因組 priority bug 主要分佈於 chr7/chr2/chr1/chr16/chr20（chr19 占比僅 2.16%，rank 19）。SP1/2/3 用於 demonstrate 機制，不代表分佈位置。」
- 圖題（Figure 5/6/7 caption）保留現狀。

---

## 2. E2 — §5.2 V5 working tree commit 狀態更新

### 2.1 PI 報告 4-29 原文（line 「V5（uncommitted working tree）」+ 「V5 working tree 未 commit」caveat）

> **V5（uncommitted working tree）** ｜ working tree ｜ getVote() Layer 1.5 somatic fallback + countSNPHaplotype() 對稱 alt guard ｜ Layer 1.5：`HaplotagProcess.cpp:512-563` +15/-1；alt guard：`HaplotagProcess.cpp:489-494` +9/-6 ｜ AMB% 17.5→8.0%；HP:i:33 −54%（239,679 → 110,197）

> 但 V5 working tree 未 commit、Confidence threshold 0.6 未直接驗證、7 樣本擴展未做、cnLOH 雙親同源區仍 open

當時（2026-04-29）V5 改動仍在 working tree，未 commit。

### 2.2 5/9 V5 commit 狀態更新

V5 chain 已完成 commit：

| Commit | 日期 | 內容 |
|---|---|---|
| `8b8c1fd` | 2026-04-09 | feat: --pon-only-phasing flag |
| `41ff147` | 2026-04-10 | fix(haplotag): two-layer getVote |
| `380e8d2` | 2026-04-25 | fix(haplotag): countINDELHaplotype UNDEFINED guard |
| `d0bcd8c` | 2026-04-30 | **fix(purity): collect ploidyRatio + bundled Layer 1.5 + countSNPHaplotype guard** |
| `938f0df` | 2026-04-30 | **Update purity threshold 0.95→0.9** |

V5 為 5 commits（不是 4），全部 commit 完成；HEAD = `938f0df` = 最新有效演算法版本（origin/main 後續 10 commits 全 doc / dead-code）。

### 2.3 修訂建議

- §5.2 表「V5（uncommitted working tree）」→ **「V5（commit `d0bcd8c` + `938f0df` 完成於 2026-04-30）」**
- §1 caveat「V5 working tree 未 commit」→ **「✅ 2026-04-30 已 commit（d0bcd8c + 938f0df）」**
- §5.2 表 4 commit → **5 commit**（補 d0bcd8c, 938f0df）

---

## 3. E3 — §5.2 priority bug 機制證據強度升級

### 3.1 PI 報告 4-29 原文（§5.2 V5 tag 修補機制 + §3.3 3 層證據鏈）

§3.3 「3 層獨立證據鏈」涵蓋理論層 / 全基因組層（17.3:1）/ 個別位點層（SP1/2/3 IGV）。**read-level 個案層尚未涵蓋**。

§5.2 「getVote priority bug」機制證據主要為 commit message 描述 + 3 個 IGV 截圖。

### 3.2 5/8 整合報告 §6 / T1.2 + T1.2-F1 補強

對 baseline / V3F / V5 三版 testing-only binary patch 加 `--debug-vote-dump` flag，dump 每條 read 經 `getVote()` 後的 5-vote countMap + hpResult：

| 規模 | chr19 pilot | 全基因組 |
|---|---:|---:|
| Dump rows | 549,206 | 29,973,253 |
| **Priority bug victims** | **752** | **34,855** |
| V3F 修正比例 | 100.00% | 100.00% |
| V5 修正比例 | 100.00% | 100.00% |
| 4-path 驗證 | 3.5/4 PASS | — |

### 3.3 修訂建議

- §3.3 證據鏈加第 4 層「read-level 個案層」：34,855 全基因組 priority bug confirmed victims，V3F + V5 修正率 100%
- §5.2 「priority bug 機制」段加註：「**read-level audit (T1.2 + T1.2-F1) 進一步驗證**：chr19 752 victims + 全基因組 34,855 victims，全部單向 baseline=11→V3F=21→V5=21，V3F + V5 修正率 100% — 從『理論 + IGV 3 截圖』升級為『個案 + 統計 + 機制』三重佐證。詳見 [T1.2 mechanism report](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md) + [T1.2-F1 全基因組擴展](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md)」

---

## 4. E4 — §6.4/§6.5 V5 數值歸因精確化（最重要 errata）

### 4.1 PI 報告 4-29 原文（§6.4 + §6.5）

> §6.4 全基因組 paired ground-truth concordance（PI 報告 4 §3.7）：全基因組 clean PS blocks Baseline=82.2% / V5=90.5% / Δ=+8.3 pp；15-site Aggregate 72.20%→78.85%（+6.65 pp）；15-site Clean PS（11 sites）74.9%→88.2%（+13.3 pp）

> §6.5 全基因組 HP 結構改善：HP1:HP2 17.3:1→~1:1；Phase block N50 4,061→8,109（+99.7%）；Phased rate 54.9%→78.5%（+23.6 pp）；執行時間 2,693s→1,976s（1.36× 快）

數值正確，但隱含歸因「V5 four-commit chain 的整體效益」不精確。

### 4.2 5/5 V5 audit 揭露 + 5/8 整合報告 §7 補強

PI 報告引用之 V5 BAM 是 `output/pononly_v5_somatic_fallback/tumor_tagged.bam`（4-12，吃 4-3 V2b phased VCF）：

```
4-12 V5 BAM 對應的 phasing.log:
   purity: 0       ← ploidy bug 讓 ploidyRatioMap 留空 → q1=q3=0 → polynomial → clamp 到 0
   highPurity = (0 > 0.9) = false
   Pass 2 second round phasing 從未觸發
   = Pass 1 only 結果
```

→ PI 報告 §6 全部 V5 數值均為 **Pass 1 only** 結果，主要功勞來自：
- **V3F two-layer `getVote()`**（commit `41ff147`）— tagging layer 修對 priority bug
- **V5 Layer 1.5 somatic fallback**（commit `d0bcd8c` bundled）— germline 缺席區補 fallback

**Pass 2 second round 二次效益**（commit `d0bcd8c` ploidy fix + `938f0df` threshold）**尚未在 PI 報告數值中體現**。

5/8 整合報告 §7.3 量化各功勞來源：

| PI 報告數值 | 主要功勞來源 | Pass 2 觸發後預期 |
|---|---|---|
| HP1:HP2 17.3:1 → ~1:1 | V3F tagging fix | 應持平 |
| 94.6% somatic→HP1 → ~50% balanced | V3F | 應持平 |
| sanity 15/15 PASS | V3F + Layer 1.5 | 應持平 |
| clean PS paired GT concordance +13.3 pp | V3F + Layer 1.5 | 持平 or 微升 |
| HP:i:33 −54% (110,197) | V3F + Layer 1.5 | 應持平 |
| Phase block N50 +99.7% | 主要 PON-only Pass 1 | Pass 2 後可能更升 |

獨立量化（4-30 重跑後 Pass 1 vs Pass 1+2 同 sample 對比，5/8 整合報告 §8.5.3）：

| 指標 | Pass 1 only | Pass 1+2 | Δ |
|---|---:|---:|---|
| Phase blocks | 1,808 | 1,631 | −177 (−9.79%) |
| Phased variants | 1,848,538 (58.00%) | 1,756,339 (55.10%) | −92,199 (−2.90 pp) |
| **N50 (bp)** | **11,388,114** | **11,788,053** | **+399,939 (+3.51%)** |

→ Pass 2 incremental N50 = +3.51%（從 4-12 到 4-30 重跑），相比於 baseline → V5 整體 N50 +99.7%（PI §6.5）只是小幅。

### 4.3 修訂建議

- §6.4/§6.5 表後加註：「**歸因精確化**：PI 報告引用之 V5 BAM 為 4-12 `output/pononly_v5_somatic_fallback/tumor_tagged.bam`，實為 Pass 1 only（ploidy bug 讓 purity=0，Pass 2 second round 從未觸發）。表中 V5 數值主要功勞來源為 **V3F two-layer `getVote()`（commit `41ff147`）+ Layer 1.5 somatic fallback（`d0bcd8c` bundled）**，Pass 2 second round（`d0bcd8c` ploidy fix + `938f0df` threshold）二次效益尚未在此數值中體現。Pass 1 vs Pass 1+2 incremental 量化見 [Self-Phasing 完整觀察整合報告 §8.5.3](../2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md)。」
- §1 一句結論加註：「PI 報告 V5 數值為 Pass 1 only；4-30 重跑後完整 V5 (Pass 1 + Pass 2) BAM 已產出但未跑 ISM benchmark（用戶 5/7 決策 cancel — ISM 是下游消費者，longphase-to 端 V3F 修對後 ISM 自動受惠）」
- §1 caveat 「Confidence threshold 0.6 未直接驗證」→ 補「但 6 cell × 4 metric caller F1 vs SEQC2 truth 三版完全相同（HCC1395 5kHz 0.93 = 0.7166 / t30_n20 0.6 = 0.6273），caller F1 與 V5 tag 改動無關（V5 不改 caller）」

---

## 5. 修訂後 §1 一句結論建議全文

（**修訂後完整版**，可整段替換 PI 報告 line 47）：

> **Self-phasing 在 baseline LongPhase-TO 是真實 tag-level artifact（94.6% somatic ALT reads 集中於 HP1，HP1:HP2 = 17.3:1）；修補在 longphase-to-mod fork 內以 5 commits 漸進完成（`8b8c1fd` PON-only flag → `41ff147` getVote 兩層 → `380e8d2` INDEL guard → `d0bcd8c` ploidy fix + bundled Layer 1.5 → `938f0df` threshold 0.95→0.9，全部 commit 完成於 2026-04-30），總修補集中於 3 函式、`HaplotagProcess.h:66-68` 介面契約零變動；V5 通過 4 項硬性 sanity check 15/15 PASS、clean PS paired GT concordance +13.3 pp（全基因組 PI 報告 4 +8.3 pp），可作為新 haplotag baseline。Read-level audit (T1.2 + T1.2-F1) 補強 priority bug 機制因果至 chr19 752 + 全基因組 34,855 victims 鐵證、V3F+V5 修正率 100%。但 PI 報告 §6 V5 數值為 4-12 BAM = Pass 1 only（ploidy bug 讓 purity=0 從未觸發 Pass 2），主要功勞來自 V3F + Layer 1.5（tagging layer），Pass 2 second round 二次效益尚未獨立量化（用戶 5/7 決策 cancel ISM benchmark — ISM 為下游消費者）；caller F1 vs SEQC2 truth 三版完全相同（HCC1395 0.93 = 0.7166 / 0.6 = 0.6273，V5 不改 caller）；7 樣本擴展未做、cnLOH 雙親同源區仍 open。InterSubMod 在這條修補鏈是下游消費者而非實作者，本 repo 無 C++ 改動。**

---

## 6. 修訂歷程

| 時間 | 修訂事件 | 來源 |
|---|---|---|
| 2026-04-29 | PI 報告原文 commit | [`InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`](../2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md) |
| 2026-04-30 | V5 chain 完成 commit (d0bcd8c, 938f0df) | longphase-to-mod git log |
| 2026-05-05 | V5 audit 揭露 PI 報告 V5 數據 = Pass 1 only | [`InterSubMod/docs/reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md`](20260505_self_phasing_V5_data_provenance_audit_01.md) |
| 2026-05-07 | T1.2 chr19 752 victims read-level 鐵證 | [`InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md`](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md) |
| 2026-05-08 | T1.2-F1 全基因組 34,855 + 顛覆三項 chr19 結論 | [`InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md`](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md) |
| 2026-05-08 | 5/8 Self-Phasing 整合報告（4 條 errata 列出於 §9.2） | [`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`](20260508_Self_Phasing_完整觀察整合報告_01.md) |
| **2026-05-09** | **本 erratum companion + PI 報告頂部 banner** | 本檔 |

---

## 7. 引用文件

- [PI 報告原文 (4-29)](../2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md)
- [V5 Audit (5-05)](20260505_self_phasing_V5_data_provenance_audit_01.md)
- [Self-Phasing 完整觀察整合報告 (5-8)](20260508_Self_Phasing_完整觀察整合報告_01.md)
- [T1.2 chr19 priority bug mechanism (5-7)](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md)
- [T1.2-F1 全基因組 audit (5-8)](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md)
