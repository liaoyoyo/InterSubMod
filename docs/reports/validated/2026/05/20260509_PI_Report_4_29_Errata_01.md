<!--
build_date: 2026-05-09
revised: 2026-05-10 (加 E5 V5 Layer 1.5 設計缺陷 — paired germline-absent xref 5/9 Step D 新發現)
agent: errata patch (從 5/8 整合報告 §9.2 + 5/9 V5 commit 狀態更新 + 5/10 Step D 補強)
status: validated
report_class: errata-companion (修訂 4-29 PI 報告 5 處表述)
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
verdict: 5 條 errata，主要 PI 結論不撤回；補強為「機制 + 案例 + 統計」三重佐證；歸因從「V5 整體」精確化為「V3F + Layer 1.5 為主，Pass 2 二次效益尚未獨立量化」；E5 (5/10 加) — V5 Layer 1.5 在 germline-absent 區域與 baseline 4.19:1 偏 HP1 完全相同，是 priority bug 的 feature 化非修補，V3F 標 hp=33 反而更穩健
last_verified: 2026-05-10
report_template: errata-companion v1.0
-->

# Errata for PI Report 2026-04-29 — Self-Phasing 整合報告（5/8）後續修訂

## 0. TL;DR

PI 報告 [`InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`](../2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md)（2026-04-29）**主要結論不撤回**（self-phasing 真實存在、V3F + V5 修正確立、V5 可作 production tag baseline）。本 erratum 修訂 5 處表述：

1. **E1** §3.3.3 chr19 SP1/2/3 解讀：從「主要 hotspot」→「可重現案例（chr19 占全基因組 priority bug 僅 2.16%）」
2. **E2** §5.2 V5 commit 狀態：「V5 working tree uncommitted」→ 已 commit (`d0bcd8c` 4-30 含 ploidy fix + bundled Layer 1.5 + countSNP guard；`938f0df` 4-30 含 threshold 0.95→0.9)
3. **E3** §5.2 priority bug 機制證據強度：升級「commit message + 3 IGV 截圖」→「+ 34,855 read-level victims 全基因組鐵證、V3F+V5 修正率 100%」
4. **E4** §6.4/§6.5 V5 數值歸因：精確化「V5 four-commit chain 整體效益」→「**主要 V3F + Layer 1.5（tagging layer）；Pass 2 second round 二次效益尚未獨立量化**」（PI 報告引用之 4-12 V5 BAM 實為 Pass 1 only，因 ploidy bug 未修正前 purity=0 → highPurity=false → Pass 2 從未觸發）
5. **E5（5/10 加）** §5.2 V5 Layer 1.5 設計描述：「Layer 1.5 = germline 缺席區域的 fallback（隱含「修補」）」→ **在 germline-absent 區域 V5 Layer 1.5 與 baseline 4.19:1 偏 HP1 完全相同**，是 priority bug 的 feature 化非修補；V3F 標 hp=33（純 somatic ambiguous）保守處理反而更穩健（5/9 Step D paired cross-ref 5,789 chr19 events 量化證明）

證據來源：5/8 整合報告 §7（V5 Provenance bonus）、§8（顛覆性發現）、§5（修補設計演進）；5/10 加 §8.6 + paired audit Step D。

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

## 5. E5（5/10 新加） — §5.2 V5 Layer 1.5 設計描述精確化

### 5.1 PI 報告 4-29 原文（§5.2 V5 working tree row + §3.3 三層證據鏈描述）

> **V5（uncommitted working tree）**：getVote() Layer 1.5 somatic fallback + countSNPHaplotype() 對稱 alt guard ｜ AMB% 17.5→8.0%；HP:i:33 −54%（239,679 → 110,197）

§5.2 對 Layer 1.5 的描述隱含「germline 缺席區域 Layer 1.5 提供修補」的敘事，但**未量化「修補」是否真存在於該區域**。

### 5.2 5/9 paired germline-absent xref 揭露（5/8 §8.6 + paired audit Step D）

對 paired chr19 read_name × T1.2 baseline / V3F / V5 vote dump JOIN，篩 `cnt_HP1+cnt_HP2=0` 且 `somatic>0` 的 events（5,789 events）：

| 版本 | hp=11 (HP1 系列) | hp=21 (HP2 系列) | hp=33 (somatic ambiguous) | 比例 / 評語 |
|---|---:|---:|---:|---|
| **baseline** | **3,312** | 791 | 80 | **4.19:1 偏 HP1**（priority bug 次峰）|
| **V3F** | 0 | 0 | **5,789** | 全標 hp=33（保守不選邊）✅ |
| **V5** | **3,313** | 790 | 53 | **4.19:1 偏 HP1（與 baseline 完全相同！）** ⚠️ |

→ **V5 Layer 1.5 在 germline-absent 區域行為與 baseline 完全相同**：

- baseline 4.19:1 偏 HP1 是 priority bug 在 somatic-only 投票場景的次峰偏移
- V5 Layer 1.5 設計（germline 缺席時用 `somaticHP1` vs `somaticHP2` 票數決方向）在 self-phasing 機制下 — sub-clone somatic 100% 共現 → graph 偏向同一 haplotype → 投票偏向同邊
- → **V5 Layer 1.5 = priority bug 的 feature 化**：把「baseline 用 somatic vote 蓋過 germline」buggy 行為，改成「germline 缺席時才用 somatic vote」designed 行為，但**該區域偏移本質沒變**

V3F 在該區域標 hp=33（純 somatic ambiguous，方向不選邊）— **保守正確**避免錯標方向。

### 5.3 修訂建議

- §5.2 表「**V5（uncommitted working tree）**」row 加註：「**caveat (5/10 補)**：在 germline-absent 區域 Layer 1.5 與 baseline 4.19:1 偏 HP1 完全相同（5/9 Step D paired cross-ref 量化），是 priority bug 的 feature 化而非修補；V3F 在該區域標 hp=33 反而更穩健。完整量化見 [Self-Phasing 完整觀察整合報告 §8.6](20260508_Self_Phasing_完整觀察整合報告_01.md)。」
- §5.2 V5 修補敘事改寫：把「V5 含 Layer 1.5 補 germline 缺席區域 fallback」精確化為「V5 Layer 1.5 嘗試在 germline 缺席時用 somatic phased votes 補方向，但該設計在 self-phasing 機制下繼承 priority bug 偏移；germline-absent 區域真正穩健的選擇是 V3F 標 hp=33」
- §1 caveat 加：「V5 Layer 1.5 設計選擇待 ISM 影響量化（F-paired-D3）— 改回 V3F「germline 缺席標 hp=33」可能是更安全 default」

### 5.4 不影響的 PI 結論

- 17.3:1 全基因組偏移確立 ✅（priority bug 主要由 germline_vote>0 區域貢獻，整體 ratio 由該區域主導）
- V3F + V5 修正 chr19 752 / 全基因組 34,855 victims 100% 修正率 ✅（這是 germline_vote>0 區域的修正）
- sanity 15/15 PASS ✅
- caller F1 三版相同 ✅
- E5 僅影響 V5 vs V3F 在 **germline-absent 區域** 的設計選擇，**不撤回任何主結論**

### 5.5 對 region-level downstream 的影響評估（V6-C Phase B 補強，5/10）

V5 Layer 1.5 的 read-level 4.19:1 偏移**在 region-level 後續特徵化（如 HPFineNGroups marker）影響有限**：

| 量化 | 結果 | 說明 |
|---|---|---|
| chr19 marker filter (NG≥3) flag=off TP rate | 94.7% (463/489) | germline-existent 區為主，priority bug 已被 V3F/V5 修正 |
| 對應 flag=on (NG_on=2) cell TP rate | 91.5% (367/401) | bucket schema collapse 後仍保留 ≥0.85 |
| 最強 cell (NG_off=5 → NG_on=2) | 99.2% (122/123) | schema 訊號塌陷後物理屬性仍區辨 |

→ **V5 Layer 1.5 設計缺陷在 read-level 確實存在（germline-absent 區），但 region-level marker 在 flag=on/off 下都通過 0.85 gate**：
- read-level audit（如全基因組 priority bug 統計）仍受 4.19:1 偏移污染
- region-level downstream feature（如 HPFineNGroups subclone marker）對該區的 reads 共現不敏感（germline-absent 區 events 占比小，且 marker filter 多落在 germline-existent 區）

完整量化 → [`InterSubMod/research/paired_priority_bug_audit/05_V6C_phaseB_findings.md`](../../../research/paired_priority_bug_audit/05_V6C_phaseB_findings.md)。

→ E5 對 V5 production usage 的結論：**Layer 1.5 設計缺陷在 read-level 真實，但 region-level 影響輕微**；不阻擋 V5 作為 production tag baseline；germline-absent 區改回 V3F 的 ISM 影響待 Phase C 7 樣本量化（F-paired-D3 follow-up）。

---

## 6. 修訂後 §1 一句結論建議全文

（**修訂後完整版**，可整段替換 PI 報告 line 47）：

> **Self-phasing 在 baseline LongPhase-TO 是真實 tag-level artifact（94.6% somatic ALT reads 集中於 HP1，HP1:HP2 = 17.3:1）；修補在 longphase-to-mod fork 內以 5 commits 漸進完成（`8b8c1fd` PON-only flag → `41ff147` getVote 兩層 → `380e8d2` INDEL guard → `d0bcd8c` ploidy fix + bundled Layer 1.5 → `938f0df` threshold 0.95→0.9，全部 commit 完成於 2026-04-30），總修補集中於 3 函式、`HaplotagProcess.h:66-68` 介面契約零變動；V5 通過 4 項硬性 sanity check 15/15 PASS、clean PS paired GT concordance +13.3 pp（全基因組 PI 報告 4 +8.3 pp），可作為新 haplotag baseline。Read-level audit (T1.2 + T1.2-F1) 補強 priority bug 機制因果至 chr19 752 + 全基因組 34,855 victims 鐵證、V3F+V5 修正率 100%。但 PI 報告 §6 V5 數值為 4-12 BAM = Pass 1 only（ploidy bug 讓 purity=0 從未觸發 Pass 2），主要功勞來自 V3F + Layer 1.5（tagging layer），Pass 2 second round 二次效益尚未獨立量化（用戶 5/7 決策 cancel ISM benchmark — ISM 為下游消費者）；caller F1 vs SEQC2 truth 三版完全相同（HCC1395 0.93 = 0.7166 / 0.6 = 0.6273，V5 不改 caller）。**5/9 paired cross-ref 揭露 V5 Layer 1.5 在 germline-absent 區域與 baseline 4.19:1 偏 HP1 完全相同（priority bug 的 feature 化而非修補），V3F 標 hp=33 反而更穩健 — 該區域 V5 設計選擇待 ISM 影響量化（F-paired-D3）；但 germline-absent 區域佔比小，不阻擋 V5 作為整體 production baseline**。7 樣本擴展未做、cnLOH 雙親同源區仍 open。InterSubMod 在這條修補鏈是下游消費者而非實作者，本 repo 無 C++ 改動。**

---

## 7. 修訂歷程

| 時間 | 修訂事件 | 來源 |
|---|---|---|
| 2026-04-29 | PI 報告原文 commit | [`InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`](../2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md) |
| 2026-04-30 | V5 chain 完成 commit (d0bcd8c, 938f0df) | longphase-to-mod git log |
| 2026-05-05 | V5 audit 揭露 PI 報告 V5 數據 = Pass 1 only | [`InterSubMod/docs/reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md`](20260505_self_phasing_V5_data_provenance_audit_01.md) |
| 2026-05-07 | T1.2 chr19 752 victims read-level 鐵證 | [`InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md`](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md) |
| 2026-05-08 | T1.2-F1 全基因組 34,855 + 顛覆三項 chr19 結論 | [`InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md`](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md) |
| 2026-05-08 | 5/8 Self-Phasing 整合報告（4 條 errata 列出於 §9.2） | [`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`](20260508_Self_Phasing_完整觀察整合報告_01.md) |
| **2026-05-09** | **本 erratum companion + PI 報告頂部 banner（E1-E4）** | 本檔 |
| 2026-05-09 | Paired audit Step A+C cycle 36 — paired 沒 priority bug | [`InterSubMod/research/paired_priority_bug_audit/00_audit_report.md`](../../../research/paired_priority_bug_audit/00_audit_report.md) (commit 6ed8a0d) |
| 2026-05-09 | Paired audit Step D cycle 37 — V5 Layer 1.5 設計缺陷 | [`InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md`](../../../research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md) (commit 766ec5f) |
| 2026-05-10 | 5/8 整合報告補 §8.6 Paired Mode Cross-Reference Audit | [§8.6](20260508_Self_Phasing_完整觀察整合報告_01.md) (commit df5137e) |
| **2026-05-10** | **本檔加 E5 + renumber §5/§6/§7 → §6/§7/§8** | 本檔 |
| **2026-05-10** | **本檔 §5.5 V6-C Phase B chr19 region-level marker robustness 補強** | [V6-C Phase B](../../../research/paired_priority_bug_audit/05_V6C_phaseB_findings.md) |

---

## 8. 引用文件

- [PI 報告原文 (4-29)](../2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md)
- [V5 Audit (5-05)](20260505_self_phasing_V5_data_provenance_audit_01.md)
- [Self-Phasing 完整觀察整合報告 (5-8)](20260508_Self_Phasing_完整觀察整合報告_01.md) — 含 §8.6 Paired Cross-Ref Audit
- [T1.2 chr19 priority bug mechanism (5-7)](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md)
- [T1.2-F1 全基因組 audit (5-8)](../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md)
- [Paired audit Step A+C (5-9)](../../../research/paired_priority_bug_audit/00_audit_report.md)
- [Paired audit Step D germline-absent (5-9)](../../../research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md)
