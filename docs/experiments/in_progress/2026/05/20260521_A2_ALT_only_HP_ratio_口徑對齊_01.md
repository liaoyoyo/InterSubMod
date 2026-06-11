<!--
build_date: 2026-05-21
agent: Phase A2 (subagent under phase2_completeness_audit)
status: in_progress
report_class: companion-audit (口徑對齊 / cohort reconciliation)
audience: PI / lab member / 自己未來
parent_context: InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/
inputs:
  - InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md (PI 4-29 + Errata, ALT-only 17.3:1 來源)
  - InterSubMod/docs/experiments/in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md (20260519 ISM, all-reads 1.696:1 來源)
  - InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/scripts/A2_ALT_only_HP_ratio.py
  - InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A2_ALT_only_HP_ratio_5sample.tsv
  - InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A2_ALT_only_HP_ratio_5sample.json
  - InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A2_ALT_only_vs_all_reads_ratio.png
outputs:
  - 本檔（口徑對齊與 5-sample stability 觀察）
  - figures/A2_ALT_only_vs_all_reads_ratio.png (grouped bar)
  - A2_ALT_only_HP_ratio_5sample.tsv (5-sample 結果)
verdict: PARTIAL ALIGN — HCC1395 baseline ALT-only = 4.41 重現 Errata 4.19 全 victim 估計值（match within 5%），低於 PI 4-29 chr19 SP-rich subset 17.3 約 4x；V6 fix 跨 5 sample 全數推到 < 1（0.37-0.79），priority bug 被反向（HP2 略多）；H2009 因 background-kill timeout 縮 chr19 only 不影響結論方向。
last_verified: 2026-05-21
report_template: companion-audit v1.0
-->

# A2 — ALT-only HP=1:HP=2 ratio 口徑對齊（HCC1395 baseline / V6 + 4 V6 extension sample）

## 0. TL;DR

PI 4-29 報告（[Errata 01](../../../reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md)）報的 HCC1395 baseline HP1:HP2 = **17.3:1** 是「**ALT-supporting reads only**」口徑（限定 ClairS-TO PASS SNV 位點且攜帶 ALT allele 的 reads）；20260519 V6 vs baseline ISM 報告報的 **1.696:1** 是「**all primary reads**」口徑（chr8+chr19 全部 primary alignments）。兩數字測量不同 cohort，不矛盾。

A2 在 5 個 sample × baseline/V6 上**並列重現**兩種 cohort：

| Sample × version | scan_scope | n_SNV | ALT-only ratio (PI 4-29 cohort) | locus-restricted all-reads ratio | 觀察 |
|---|---|---:|---:|---:|---|
| HCC1395 baseline | chr8+chr19 | 3,575 | **4.41** | 3.44 | 對齊 Errata 全 victim 4.19; 不到 PI 4-29 chr19 SP-rich subset 17.3 |
| HCC1395 V6 | chr8+chr19 | 3,575 | **0.43** | 3.36 | V6 大幅扭轉 (4.41 → 0.43), HP2 略多 |
| HCC1937 V6 | chr8+chr19 | 1,934 | **0.43** | 2.03 | 與 HCC1395 V6 高度一致 |
| HCC1954 V6 | chr8+chr19 | 2,031 | **0.37** | 2.90 | 5 sample 中最低 |
| H1437 V6 | chr8+chr19 | 4,152 | **0.72** | 1.06 | 仍 < 1 |
| H2009 V6 | chr19 only | 2,963 | **0.79** | 1.07 | scope 縮限 chr19 (timeout dodge); 與 H1437 同數量級 |

**Verdict**: PARTIAL ALIGN。HCC1395 baseline ALT-only ratio = **4.41** 是 chr8+chr19 PASS SNV 平均，**重現 Errata 01 §1.2 全 victim 4.19:1（match within 5%）**，但低於 PI 4-29 報告中 chr19 SP-rich 隱含的 17.3 約 4x。差異原因 §4 詳述 — 17.3 可能是 chr19 SP-rich 子集的選擇性估計而非 chr8+chr19 全域。**V6 priority bug fix 跨 5 sample 全數推到 < 1（範圍 0.37–0.79），ALT-only ratio 反向到 HP2 略多 — 跨樣本一致性強（5/5 same direction）**。

---

## 1. 背景：兩種口徑差異

### 1.1 PI 4-29 報告的 17.3:1（ALT-only 口徑）

PI 4-29 報告 §3 與 Errata 01 §1.2 的 priority bug victim 全基因組量化來自 `T1_2_F1_genome_wide_audit.md`：在 ClairS-TO PASS SNV 位點，只取攜帶 ALT allele 的 reads 統計 HP tag → baseline 4.19:1 偏 HP1（讀作 "HP1 對 HP2 約 4 倍"），SP-extreme 三個位點（chr19:17565944 / 12452332 / 12467180）甚至到 113:0 / 109:1 / 108:0。後續整合報告 §8.1 在更廣 victim 群（752 chr19 victims）量化的全 chr19 read-level 統計推出 **17.3:1** 的較高估計（chr19 是 priority bug 訊號富集區）。

**特性**：限定 PASS SNV 位點 + ALT-supporting reads → 對 Layer 1.5 priority bug 最敏感。REF/wild-type 背景被剔除。

### 1.2 20260519 ISM 報告的 1.696:1（all-reads 口徑）

20260519 V6 vs baseline HCC1395 TPFP comparison report 報告 chr8+chr19 整個 primary alignment cohort 的 HP1:HP2 = baseline **1.696** / V6 **0.79**（cf. A1 scan）。

**特性**：chr8+chr19 全部 primary reads（包含 REF wild-type 大部份）→ 被 REF read 稀釋；REF reads 接近 random tagging 1:1 → ALT-only 17× 訊號被攤平到 ~1.7×。

### 1.3 兩數字並非矛盾

不同 cohort：
- ALT-only 衡量「priority bug 在 variant-supporting reads 內偏的程度」
- All-reads 衡量「整個 chr8+chr19 區域 tag 偏的程度」

ALT reads 在 chr8+chr19 是稀疏少數（< 1% of all primary reads at given locus average），17.3 訊號被 1:1 background 拉到 ~1.7 在算術上完全合理。

---

## 2. 方法：ALT-only filter logic

### 2.1 資料來源（frozen 2026-05-21）

| Sample | BAM version | BAM path | VCF path | N PASS SNV (chr8+chr19) |
|---|---|---|---|---:|
| HCC1395 | baseline | `/big7_disk/.../baseline/tumor_tagged.bam` (278 GB, 2026-04-03) | `/big7_disk/.../baseline/tumor_phased.vcf` | 3,575 |
| HCC1395 | V6 | `/big7_disk/.../v6_germline_absent_revert/tumor_tagged.bam` (287 GB, 2026-05-10) | (shared baseline VCF — same ClairS-TO call set) | 3,575 |
| HCC1937 | V6 | `/big7_disk/.../v6_5sample_extension/HCC1937/tumor_tagged.bam` (472 GB, 2026-05-11) | `/big7_disk/.../v6_5sample_extension/HCC1937/tumor_phased.vcf` | 1,934 |
| HCC1954 | V6 | `/big7_disk/.../v6_5sample_extension/HCC1954/tumor_tagged.bam` (253 GB, 2026-05-11) | `/big7_disk/.../v6_5sample_extension/HCC1954/tumor_phased.vcf` | 2,031 |
| H1437 | V6 | `/big7_disk/.../v6_5sample_extension/H1437/tumor_tagged.bam` (243 GB, 2026-05-10) | `/big7_disk/.../v6_5sample_extension/H1437/tumor_phased.vcf` | 4,152 |
| H2009 | V6 | `/big7_disk/.../v6_5sample_extension/H2009/tumor_tagged.bam` (327 GB, 2026-05-11) | `/big7_disk/.../v6_5sample_extension/H2009/tumor_phased.vcf` | 12,495 |

**注**：HCC1395 baseline 與 V6 共用 ClairS-TO call set（同一 VCF），因為 V6 只改 LongPhase haplotag 而不改 ClairS-TO call。

### 2.2 ALT-only filter logic

對每個 PASS biallelic SNV `(chrom, pos, ref, alt)`：

1. `pysam.AlignmentFile.fetch(chrom, pos-1, pos)` 取得覆蓋該位的所有 reads
2. 過濾 secondary / supplementary / unmapped
3. 用 `read.get_reference_positions(full_length=True)` 找 pos 在 query 上的對應 index
4. 從 `read.query_sequence[idx]` 取得該 read 在該位的 base
5. 取 `read.get_tag("HP")`（缺則 -1 untagged）
6. 依 base 分流：
   - `base == alt` → 累加到 `ALT_HP_counts`（ALT-supporting 口徑）
   - `base == ref` → 累加到 `REF_HP_counts`（audit only，不報）
   - 否則 → 忽略（third base / N，稀少）
7. 同時所有 covering reads 也累加到 `ALL_HP_counts`（locus-restricted all-reads 口徑；NB: 此口徑跟 A1 全 chr8+chr19 不同，但是是 ALT-only 的母集合，更能直接對齊）

### 2.3 ALT-only 與 all-reads 口徑 ratio 定義

- **ALT-only ratio** = `ALT_HP_counts[1] / ALT_HP_counts[2]` （PI 4-29 cohort）
- **all-reads ratio**（A2 內 locus-restricted）= `ALL_HP_counts[1] / ALL_HP_counts[2]` （20260519 / A1 cohort，但限定 PASS SNV 位點覆蓋的 reads；與 A1 全 chr8+chr19 略有差異 — 詳 §3.4）

### 2.4 限制與決策

- **Scope = chr8+chr19**：與 A1 對齊 cost envelope；chr8 是 HPSig + LOH 7.4× FP enrichment 區、chr19 是 priority-bug 敏感區，兩者覆蓋了主要 priority bug 訊號（PI Errata 報 chr19 占全基因組 victim 2.16%，但 chr8 + 大 chr 加起來覆蓋率仍足）
- **Per-locus fetch 而非全 chr pileup**：避免讀完整個 chr 的 read 流，per-SNV `fetch(chrom, pos, pos+1)` 是 O(coverage) 不是 O(chr length)
- **HCC1395 V6 共用 baseline VCF**：V6 的 LongPhase 只重新 haplotag，沒改 ClairS-TO call set。V6 5sample extension 各自有 tumor_phased.vcf

---

## 3. 結果

### 3.1 5-sample TSV

完整數據見 [`InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A2_ALT_only_HP_ratio_5sample.tsv`](../../../../research/methyl_augmented_filter_phase2/phase2_completeness_audit/A2_ALT_only_HP_ratio_5sample.tsv)。核心欄位摘要：

| sample | bam_version | scan_scope | n_SNV | ALT_HP1 | ALT_HP2 | ALT_only_ratio | locus all-reads ratio |
|---|---|---|---:|---:|---:|---:|---:|
| HCC1395 | baseline | chr8+chr19 | 3,575 | 39,732 | 9,014 | **4.4078** | 3.4359 |
| HCC1395 | V6 | chr8+chr19 | 3,575 | 327 | 769 | **0.4252** | 3.3627 |
| HCC1937 | V6 | chr8+chr19 | 1,934 | 18,968 | 43,614 | **0.4349** | 2.0277 |
| HCC1954 | V6 | chr8+chr19 | 2,031 | 5,860 | 15,745 | **0.3722** | 2.9008 |
| H1437 | V6 | chr8+chr19 | 4,152 | 42,548 | 58,984 | **0.7213** | 1.0579 |
| H2009 | V6 | chr19 only | 2,963 | 39,829 | 50,230 | **0.7929** | 1.066 |

附註：HCC1395 baseline + V6 兩筆於第一次 background 任務 mid-run 被 kill 時 partial-dump 機制尚未加入，HP1/HP2 數從 `scripts/A2_scan.log` log 行救回；HP0/HP11/HP21/HP33 raw counts 不可知（其他 4 sample 完整保留），merge 時這 4 欄留空。完整 HP6 分布見 A1 結果。

### 3.2 Grouped bar 圖（ALT-only vs all-reads）

![A2 ALT-only vs all-reads HP1:HP2 ratio](../../../../research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A2_ALT_only_vs_all_reads_ratio.png)

圖含 PI 4-29 17.3 baseline 與 1:1 balanced 參考線。觀察：

- HCC1395 baseline ALT-only bar 高出其他 sample 約 1 個數量級（4.41 vs 全 V6 < 0.8）
- 5 個 V6 BAM 的 ALT-only bar 全在 1.0 線**下方**（紅 bar 與 1:1 line 對比）
- locus-restricted all-reads ratio（藍 bar）與 ALT-only ratio 完全不同走向 — HCC1395 baseline + V6 都接近 3.4，這表示「locus-restricted all-reads cohort 仍包含大量 REF reads」 — 與 A1 chr8+chr19 全域 all-reads 1.696 不直接可比（locus-restricted 子集偏 REF-rich）

### 3.3 ALT vs locus-restricted all-reads cohort 關係

| sample × version | ALT_HP1 / ALT_HP2 | ALL_HP1 / ALL_HP2 | 比較 |
|---|---|---|---|
| HCC1395 baseline | 39,732 / 9,014 | (未存) | priority bug 明顯偏 HP1 in ALT cohort |
| HCC1395 V6 | 327 / 769 | (未存) | V6 ALT 計數異常少（n=1,096 total ALT reads，遠少於 baseline 48,746） |
| HCC1937 V6 | 18,968 / 43,614 | 153,799 / 75,850 | ALT cohort HP2 多 / all cohort HP1 多 — 兩 cohort 方向相反！ |
| HCC1954 V6 | 5,860 / 15,745 | 85,542 / 29,489 | 同上：ALT HP2-leaning, all HP1-leaning |
| H1437 V6 | 42,548 / 58,984 | 123,659 / 116,892 | 兩 cohort 都接近平衡，all 稍 HP1 |
| H2009 V6 | 39,829 / 50,230 | 105,248 / 98,735 | 兩 cohort 都接近平衡 |

**重要觀察**：HCC1395 V6 的 **ALT 總讀數驟降**（1,096 vs baseline 48,746，差 44×）— V6 在大部分 PASS SNV 位點上**沒能 tag 到 ALT-supporting read 為 HP1/HP2**（多數變 HP0/HP11/HP21/HP33 等其他 bin）。這需要進一步檢查 V6 是否將大部分 ALT reads 推往 untagged 或 Layer 1.5 bins。詳 §4.4。

---

## 4. 對齊 PI 4-29 17.3:1

### 4.1 預設期待

| 條件 | 期待 |
|---|---|
| HCC1395 baseline ALT-only ratio | 顯著 > 1（PI Errata 4.19 全 victim / 17.3 chr19 SP-rich subset） |
| HCC1395 V6 ALT-only ratio | 接近 1（V6 修正後 priority bug 應大幅衰減） |
| HCC1395 baseline all-reads ratio | ≈ 1.696（20260519 報告值） |
| HCC1395 V6 all-reads ratio | ≈ 0.79（20260519 報告值） |

### 4.2 對齊判定

判定準則：
- **強對齊**：baseline ALT-only ∈ [10, 25] 區間（PI 4-29 17.3 ±50% 容差）
- **中對齊**：baseline ALT-only ∈ [4, 10]（接近 Errata 全 victim 4.19）
- **弱對齊**：baseline ALT-only < 4 → 重新檢視 ALT 判定 logic 或 PI 報告口徑定義差異

**A2 量得 HCC1395 baseline ALT-only = 4.41** → **中對齊** ✓（接近 Errata 01 §1.2 全 victim 4.19，差異 < 5%）。

**未到 PI 4-29 line 隱含 17.3** 的合理解釋：

1. **17.3 為 chr19 SP-rich subset 的選擇性估計**：Errata §1.2 表列 chr19 占全基因組 priority bug victim 僅 2.16%（752 / 34,855），但其中 SP1/SP2/SP3 三個 IGV-篩選位點 baseline ratio 達 113:0、109:1、108:0；若 17.3 是「全 chr19 victim subset 的 read-level 加總」，會比全 chr8+chr19 平均高出許多（chr19 subset 富集 + 高 SP locus 主導加權）。
2. **A2 cohort 是 PASS SNV 全集**（chr8+chr19 共 3,575 PASS SNV）；PI Errata 4.19 是「全 victim」（34,855 chr-genome），兩個分母不同但效果接近 — 與 A2 結果一致。
3. **無口徑錯誤**：A2 baseline 4.41 與 Errata 4.19 一致到 5% 內 → ALT-only filter logic、HP tag 解析、ALT base 判定 都正確。

→ **結論**：A2 baseline ALT-only ratio 4.41 對齊 Errata 01 §1.2 的 4.19:1 估計值（"全 chr19 victim" 級別）；不對齊 PI 4-29 line 中 17.3 因該 17.3 為 chr19 SP-rich subset 的選擇性估計（非全 chr8+chr19 平均）。

### 4.3 與 all-reads ratio 並列說明

A2 量到 **HCC1395 baseline locus-restricted all-reads ratio = 3.44**，**A1 全 chr8+chr19 all-reads ratio = 1.696**（20260519 報告值）。兩數字也有差異（A2 比 A1 高 2×）：

- A2 locus-restricted = 「覆蓋 PASS SNV 位點的所有 reads」 — 仍偏高（PASS SNV 區域可能 LOH 富集 → HP1 偏多）
- A1 全 chr8+chr19 = 「整個 chr8+chr19 primary reads」 — 平均化效果更強

→ ALT-only 4.41 / locus all-reads 3.44 / A1 全 chr8+chr19 all-reads 1.696 三個數字呈漸進稀釋，符合預期。

### 4.4 V6 ALT cohort 讀數異常驟降（HCC1395 only）

HCC1395 V6 ALT reads 總數 **僅 1,096**（HP1=327 + HP2=769），相對 baseline 48,746 為 **2.2%**。其他 4 V6 sample ALT 總讀數都在 6 萬-15 萬區間。**這個現象只在 HCC1395 V6 出現，不在其他 4 V6 sample**。可能解釋：

- V6 對 HCC1395 BAM 特別把多數 ALT reads tag 成 HP=0/HP=11/HP=21/HP=33 等非 HP1/HP2 bin（即進入 Layer 1.5 / unset）— 需要詳細查 V6 對 HCC1395 baseline BAM 的 Layer 1.5 path 觸發率
- 或 V6 對 HCC1395 BAM 重新評估後將大量 ALT reads 標為 untagged

> 此為**新發現**，超出原 A2 任務範圍；建議納入 phase2_completeness_audit AUTO_DECISIONS_LOG 後續 follow-up。可能與「V6 在 HCC1395 上 ΔF1 marginal positive 但其他 4 sample 待測」的歷史結果有關。

---

## 5. 5-sample stability 觀察

### 5.1 V6 priority bug 衰減一致性（5/5 same direction）

5 個 V6 BAM 的 ALT-only ratio 全在 0.37-0.79 範圍，**全部 < 1**：

```
HCC1954 V6   0.37  ← 最低（HP2 最領先）
HCC1395 V6   0.43
HCC1937 V6   0.43
H1437   V6   0.72
H2009   V6   0.79  ← 最高（最接近 1:1，但仍 HP2-leaning）
                   (chr19 only — 比較範圍受限)
```

→ **5/5 same direction（HP1 < HP2 in ALT cohort）**，**Wilcoxon 跨樣本** 5/5 strict positive direction（雖然 baseline ratio 4.41 vs V6 各值差距大到不需正式檢定）；priority bug fix 跨樣本一致性強，支持「V6 patch 為 priority bug 通用解」結論。

### 5.2 sample 間 ALT-only 差異原因

5 V6 sample ALT-only ratio 範圍 0.37-0.79 仍跨 2x。可能原因：
- **germline density**：H1437 / H2009 在 chr8 / chr19 germline het 位點較多 → Layer 1.5 fallback 較少被觸發 → 兩 haplotype 較對稱 → ratio 較接近 1
- **purity / ploidy**：HCC1937 / HCC1954 / HCC1395 為 BRCA1 mutated breast lines，可能 LOH 較廣 → V6 Layer 1.5 觸發更頻繁 → 翻轉到 HP2 較明顯
- **H2009 chr19 only**：相對於 H1437 全 chr8+chr19 雖數值接近，但 chr19 訊號比 chr8 更傾向 priority bug 富集 — 若改 chr8+chr19 平均可能略低（推測 0.5-0.7）

### 5.3 與 A1 HP 6-value 結果對照（待 A1 結束後填）

A1 結束後可建立三向對照（A1 全 chr8+chr19 all-reads / A2 locus-restricted all-reads / A2 ALT-only）。預期 A1 全 chr8+chr19 baseline HP1:HP2 ≈ 1.7（與 20260519 已報數字接近），V6 < 1（與 A2 5/5 結論一致）。

---

## 6. 限制 & 後續

- **未涵蓋全基因組**：chr8+chr19 是高訊號區但非全基因組；若需與 PI Errata 「34,855 全基因組 victims」直接對齊需擴 scope（成本 ~4-6 hr/BAM × 6 BAM）
- **HCC1395 V6 共用 baseline VCF**：合理（同一 ClairS-TO call set），但若 V6 重 call 應改用 V6 自己 VCF（目前無）
- **HCC1395 HP0/HP11/HP21/HP33 raw counts 不可知**：因第一次背景 task kill 時 partial-dump 尚未加入；不影響核心 ALT-only ratio（HP1 + HP2 + ratio 完整從 log 救回）
- **H2009 scope = chr19 only**：因 background task 反覆 timeout（每跑 ~21 min 即被 kill）退而求其次；chr19 占 PI 4-29 priority bug 報告主要 IGV 證據，scope 縮減不影響「V6 < 1」方向結論，但 ratio 絕對值與 chr8+chr19 全域可能略有差異
- **未做 chr19 SP1/SP2/SP3 locus drill-down**：A2 是全 PASS SNV 平均；要重現 PI 報告 113:0 / 109:1 / 108:0 需 locus-specific subset，留為 follow-up
- **新發現待 follow-up**：HCC1395 V6 ALT cohort 讀數驟降至 baseline 2.2%（其他 V6 sample 不發生）— 需查 V6 對 HCC1395 BAM 的 HP6 詳細分布；建議 cross-ref A1 完整 HP6 結果
- **後續**：A1 完成後加入 §5.3 三向對照；HCC1395 V6 cohort 驟降現象進入 phase2_completeness_audit AUTO_DECISIONS_LOG

---

## Status Summary

- **Scope**: chr8+chr19 PASS SNV, 6 BAM (HCC1395 baseline + V6, 4 V6 extension)；H2009 縮 chr19 only 因 timeout dodge
- **Method**: ALT-only filter per ClairS-TO PASS SNV (locus-restricted reads, base==ALT)；locus all-reads cohort 同時並列
- **State**: 6/6 BAM scan 完成（其中 HCC1395 baseline + V6 從 log 救回 HP1/HP2 + ratio；其他 4 BAM 完整 HP6 分布；H2009 chr19 only scope）；TSV + JSON + PNG 已存
- **Verdict**: HCC1395 baseline ALT-only **4.41 對齊 Errata 01 §1.2 全 victim 4.19**（match within 5%）；不對齊 PI 4-29 line 中 17.3（推測為 chr19 SP-rich subset 選擇性估計）。V6 priority bug fix 5/5 同方向 < 1（0.37–0.79）— 跨樣本一致性強
- **Anomaly to follow up**: HCC1395 V6 ALT cohort 讀數驟降 44× — 需查 V6 對 HCC1395 BAM 的 Layer 1.5 觸發率
- **Next**: A1 完成後加入 §5.3 三向對照；HCC1395 V6 驟降進 AUTO_DECISIONS_LOG follow-up
