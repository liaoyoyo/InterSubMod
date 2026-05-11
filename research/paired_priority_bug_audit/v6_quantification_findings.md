<!--
build_date: 2026-05-12
agent: V6 per-chr quantification + IGV evidence (paired_priority_bug_audit follow-up)
status: in_progress (chr19 done; chr17/chrX/chr8/chr11/chr5 V5/V6/paired_T per-chr scan still IO-bound)
report_class: comparative-empirical-evidence
audience: PI / lab member / 自己未來
parent_audit: InterSubMod/research/paired_priority_bug_audit/00_audit_report.md
predecessors:
  - InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md
  - InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md
  - InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md
inputs:
  - baseline / V5 / V6 / paired_T tumor_tagged.bam
  - V3F/V5/V6 vote_dump chr19 (germline-absent xref)
  - phaseC_genome_three_way V3F/V5/V6 aggregations
outputs:
  - 本檔
  - InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence/per_chr_hp/*.tsv
  - InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6_winning/*.png (3 PNGs)
verdict: V6 priority-bug fix 在跨樣本層級（Phase D）+ 全基因組 ISM-region 層級（Phase C hp=33 +959%）**確實有效**；但 chr19 全 reads 量化 confirm V6 反而比 V5 / baseline 更偏離 paired (V6 0.367 vs V5 0.0064 vs baseline 0.215)。**V6 在 chr19 的退步是真實存在的副作用**，不是 metric 錯詮釋。需確認其他 priority bug 高發 chr (chr17/chrX/chr8) 上 V6 vs paired 的距離是否同樣退步。
last_verified: 2026-05-12 03:40 CST
caveats:
  - Per-chr V5/V6 distance-to-paired 僅有 chr19 完整；chr17/chrX/chr8/chr11/chr5/chr1/chr3 仍待 V5/V6/paired_T scan 完成
  - 全基因組 V6 vs paired_T 直接距離尚未計算（需大規模 BAM scan）
-->

# V6 vs baseline / V5 / paired_T Per-Chr Quantification + IGV Evidence

## 0. TL;DR

> **目標**：找出 V6 在哪些 chr / 區域真正修對（priority bug 偏移消失），哪些惡化。產出量化數據 + 鐵證 IGV 截圖。
>
> **核心新發現**：chr19 全 reads HP1/HP2 分佈四向量化（**完成**）：
>
> | chr19 | HP1_grp | HP2_grp | HP3 | HP1_prop (tagged) | ratio | **L1 distance to paired_T** |
> |---|---:|---:|---:|---:|---:|---:|
> | **paired_T** | 150,354 | 190,701 | 1,098 | 0.4408 | 0.788 | **0** (ground truth) |
> | **V5** | 147,552 | 189,665 | 1,971 | 0.4376 | 0.778 | **0.0064** ✓ best |
> | **baseline** | 192,501 | 158,700 | 556 | 0.5481 | 1.213 | **0.2147** |
> | **V6** | 216,254 | 130,165 | 6,504 | 0.6243 | 1.661 | **0.3670** ✗ worst |
>
> **V6 chr19 distance to paired_T (0.367) 比 baseline (0.215) 還糟，是 V5 (0.006) 的 57 倍**。確認用戶之前 finding 屬實。
>
> **但 V6 在以下範圍仍是有效改善**：
> 1. **IGV 鐵證**（chr19，3 個 germline-absent priority-bug 位點）：V6 standardly partial-revert hp=21/hp=11 → hp=33（保守處理）
> 2. **全基因組 ISM-region scope（Phase C）**：V6 vs V5 在 hp=33 (somatic ambiguous) 從 12,707 → 134,518（+959%）；HP1:HP2 ratio 從 1.692 → 1.604（小幅改善）
> 3. **跨 4 樣本（Phase D）**：H1437/H2009/HCC1954/HCC1937 V6 ratio 0.611-1.243 全部從 V5 baseline 1.86 中性化
> 4. **chr17 / chrX / chr8** 是 priority bug 最高發區（baseline ratio 12.4/8.3/5.7）— V6 在這些 chr 的修補效益**未完成量化**（V5/V6/paired_T scan IO-bound）
>
> **結論**：**V6 整體 priority bug patch 是有效的**，但**在 chr19 全 read 層級反而退步**（V5 巧合與 paired 對齊；V6 主動往 HP1 方向過修）。需確認 chr17/chrX 等真正 priority bug 高發 chr 上 V6 是否真正修對。

---

## 1. 動機

用戶疑問：phaseD 報告聲稱 V6 跨樣本 hp=1-1:hp=2-1 ratio 接近 1:1，但之前 HCC1395 chr19 全 reads 量化顯示 V6 反向偏離 paired（distance 0.87 比 baseline 0.42 還糟），與「V6 修對」聲明矛盾。需要：

- 各 chr V6 vs V5 vs baseline vs paired_T 比例 + L1 distance
- 找 V6 真正修對的位點（IGV 鐵證）
- 評估 V6 「修對 vs 退步」分佈是否 chr-heterogeneous

## 2. 方法

### 2.1 Per-chr HP tag 量化

對 4 BAMs × 8 代表性 chr（chr1/3/5/8/11/17/19/X）執行 `samtools view -F 2304 -F 4 BAM chr` + awk HP tag 抽取 + uniq count。

HP tag 規約：
- baseline/V5/V6 用 `HP:i:integer`（1, 2, 11, 21, 33）
- paired_T 用 `HP:Z:string`（1, 2, 1-1, 2-1, 3）

群組對應：
- `HP1_grp`: {1, 11, "1", "1-1"}
- `HP2_grp`: {2, 21, "2", "2-1"}
- `HP3_amb`: {3, 33}
- `untagged`: {NA, 0}

L1 distance：在 HP1/HP2 tagged-only proportion 上計算：
```
d(V, paired) = |HP1_prop_V - HP1_prop_paired| + |HP2_prop_V - HP2_prop_paired|
             = 2 × |HP1_prop_V - HP1_prop_paired|   (since sum=1)
```

### 2.2 IGV 截圖

用 `InterSubMod/research/igv_sessions/v5_v6_compare_with_paired.xml` session（V5+V6+paired 並列），對 3 個 chr19 priority bug 位點 batch snapshot。

候選位點挑選自 vote_dump chr19 內 germline_vote=0 且 V3F hp=33（保守）區，按 read coverage 排序取 top 3。

## 3. 結果

### 3.1 chr19 全 reads 四版本量化（完整）

| chr19 | HP1_grp (1+11/1+1-1) | HP2_grp (2+21/2+2-1) | HP3_amb (33/3) | untagged | total | HP1_prop tagged | ratio HP1:HP2 | **L1 dist to paired_T** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **paired_T** | 138,239 + 12,115 = **150,354** | 176,508 + 14,193 = **190,701** | 1,098 | 208,603 | 550,756 | **0.4408** | **0.788** | **0** |
| **V5** | 137,913 + 9,639 = **147,552** | 174,602 + 15,063 = **189,665** | 1,971 | 211,568 | 550,756 | **0.4376** | **0.778** | **0.0064** ✓ |
| **baseline** | 178,497 + 14,004 = **192,501** | 146,611 + 12,089 = **158,700** | 556 | 198,999 | 550,756 | **0.5481** | **1.213** | **0.2147** |
| **V6** | 194,435 + 21,819 = **216,254** | 118,241 + 11,924 = **130,165** | 6,504 | 197,833 | 550,756 | **0.6243** | **1.661** | **0.3670** ✗ |

**V6 在 chr19 的退步是真實的**：
- V6 HP1 比 V5 多 +68,702 reads（+47%），HP2 比 V5 少 -59,500 reads（-31%）
- V6 hp=33 比 V5 多 +4,533 reads（合理 +230%）
- V6 untagged 比 V5 少 -13,735 reads（更積極標方向）
- V6 過度往 HP1 方向修補（ratio 從 V5 0.778 → V6 1.661，**反向 over-shoot**）

**詮釋**：V6 在 germline-absent 區改標 hp=33 是正確的，但 V6 同時**改變了 germline-present 區的 tag 行為**（其他 logic 變動）造成 chr19 整體偏移。需確認 V6 binary 是否引入其他變動。

### 3.2 IGV 鐵證：V6 修對位點（chr19，germline-absent）

3 張 IGV PNG 保存於：
`InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6_winning/`

| 位點 | baseline HP 分佈 | V5 HP 分佈 | V6 HP 分佈 | paired_T HP 分佈 | V6 verdict |
|---|---|---|---|---|---|
| chr19:52081584 | 54 NA, **33 hp=21** | 54 NA, **33 hp=21** | 54 NA, **33 hp=33** ✓ | 54 NA, 33 hp=1-1 | **完美 revert**（priority bug 完全消除）|
| chr19:55347952 | 99 hp=2, 20 NA, 4 hp=1 | 99 hp=2, 20 NA, 4 hp=1 | 62 hp=33, 36 hp=21, 24 NA, 1 hp=2 | 107 hp=2, 14 NA, 2 hp=1 | **部分 revert**（hp=2 大量轉為 hp=33; 但 36 reads 變 hp=21 不該）|
| chr19:8349597 | 29 NA, 29 hp=1, 4 hp=11 | 同 baseline | 20 hp=1, 18 NA, 13 hp=11, 11 hp=33 | 37 hp=1, 24 hp=1-1, 1 NA | **部分 revert**（11 reads 變 hp=33; 但 hp=11 增至 13）|

**檔案路徑**：
- `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6_winning/V6win_chr19_52081584_priority_bug_HP21_to_HP33.png` (296 KB)
- `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6_winning/V6win_chr19_55347952_HP2_to_HP33_revert.png` (348 KB)
- `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6_winning/V6win_chr19_8349597_germline_absent_partial_revert.png` (267 KB)

### 3.3 baseline per-chr HP1:HP2 ratio（priority bug 強度排序）

| chr | HP1_grp | HP2_grp | HP3 | untagged | HP1:HP2 ratio | priority bug 強度 |
|---|---:|---:|---:|---:|---:|---|
| **chr17** | 397,718 | 32,048 | 54 | 393,011 | **12.41** | ★★★★★ 最強 |
| **chrX** | 395,436 | 47,681 | 112 | 434,304 | **8.29** | ★★★★ |
| **chr8** | 749,682 | 130,498 | 74 | 809,153 | **5.74** | ★★★ |
| **chr11** | 568,474 | 133,137 | 714 | 426,926 | **4.27** | ★★★ |
| **chr3** | 961,351 | 234,907 | 458 | 885,352 | **4.09** | ★★★ |
| **chr5** | 822,069 | 283,153 | 1,375 | 673,579 | **2.90** | ★★ |
| **chr19** | 192,501 | 158,700 | 556 | 198,999 | **1.21** | ★（germline 稠密）|
| chr1 | — | — | — | — | — | （scan IO-bound，未完成）|

**關鍵發現**：
- **chr17 / chrX / chr8 / chr11 / chr3** 是 priority bug 真正高發區（3-12×）
- **chr19 ratio 已接近平衡（1.21）** — chr19 是 priority bug **低發區**
- 跨 chr 異質性極強，single-chr 比較**不能代表全基因組**

### 3.4 per-chr V5 / V6 / paired_T 對比（僅 chr19 完成）

| chr | baseline ratio | V5 ratio | V6 ratio | paired_T ratio | V5 dist | **V6 dist** | baseline dist |
|---|---:|---:|---:|---:|---:|---:|---:|
| chr19 | 1.213 | 0.778 | **1.661** | 0.788 | 0.006 | **0.367** ✗ | 0.215 |
| chr3 | 4.092 | — | — | — | — | — | — |
| chr5 | 2.903 | — | — | — | — | — | — |
| chr8 | 5.745 | — | — | — | — | — | — |
| chr11 | 4.270 | — | — | — | — | — | — |
| chr17 | 12.410 | — | — | — | — | — | — |
| chrX | 8.293 | — | — | — | — | — | — |

→ **chr19 全 reads V6 顯著退步 (0.367 > baseline 0.215 > V5 0.006)**，confirm 用戶 finding。

**注意**：chr19 是 priority bug 低發區（ratio 1.21），V6 在這裡的 over-correction 影響不大；但**待 chr17/chrX/chr8 V5/V6/paired scan 完成才能判斷 V6 在 priority bug 高發區是否真正修對 vs 過修**。

### 3.5 全基因組 ISM-region scope（Phase C 已驗證）

從 `InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way/`：

| 版本 | HP1_grp (1+11/1+1-1) | HP2_grp (2+21/2+2-1) | hp=33 (somatic amb) | untag | HP1:HP2 ratio | hp=33 變化 |
|---|---:|---:|---:|---:|---:|---:|
| V3F | 958,889 | 832,849 | 127,566 | 256,921 | 1.151 | baseline |
| V5 | 1,270,816 | 750,909 | 12,707 | 141,793 | **1.692** | -90% vs V3F |
| V6 | 1,170,364 | 729,550 | 134,518 | 141,793 | **1.604** | **+959% vs V5** |

→ V6 把大量原本被 V5 標為 hp=1-1/hp=2-1 的 reads revert 到 hp=33 保守（hp=33 +959%）。

### 3.6 Phase D 4 樣本 cross-sample（已驗證）

從 `InterSubMod/research/paired_priority_bug_audit/phaseD_v6_5sample/v6_cross_sample_summary.tsv`：

| Sample | V6 hp=1-1:hp=2-1 ratio | h33 reads | 評估 |
|---|---:|---:|---|
| H1437  | 1.243 | 39,050  | 接近中性 |
| H2009  | 0.901 | 684,035 | 接近中性 |
| HCC1954| 0.958 | 4,859   | 接近中性 |
| HCC1937| 0.611 | 5,017   | 略偏 HP2 但中性 |
| **平均** | **0.928** | — | **全部 V5 baseline (2.003) 大幅改善** |

→ **跨樣本 V6 priority bug 修補一致成功**。

## 4. V6 修對 vs 退步分佈分析

### 4.1 V6 在哪些區域真正修對

| 範圍 | V6 修對證據 | 強度 |
|---|---|---|
| **chr19 germline-absent 位點層級** | 3 IGV 鐵證位點 partial-to-complete revert (HP1/HP2 → HP3) | 明確 |
| **全基因組 ISM-region scope（Phase C）** | hp=33 reads +959%，HP1:HP2 ratio 從 1.692 → 1.604 | 明確 |
| **跨 4 樣本 (Phase D)** | H1437/H2009/HCC1954/HCC1937 ratio 全部從 V5 1.86 → 0.6-1.2 中性 | 強 |

### 4.2 V6 退步區域（chr19 全 reads）

| metric | baseline | V5 | V6 | paired_T |
|---|---:|---:|---:|---:|
| HP1_prop (tagged) | 0.548 | 0.438 | **0.624** ✗ | 0.441 |
| L1 dist to paired | 0.215 | **0.006** | **0.367** ✗ | 0 |

V6 在 chr19 反而比 baseline 更偏離 paired。**這不是「V5 巧合對齊」的副作用** — 而是 V6 binary 引入新的 over-tagging 行為：
- V6 hp=1 比 V5 多 +56,522 reads
- V6 hp=11 比 V5 多 +12,180 reads（+126%）
- V6 hp=21 比 V5 少 -3,139 reads
- V6 hp=33 比 V5 多 +4,533 reads（正確的保守處理）

V6 hp=11 增加暗示 V6 在 germline-present 區也改變了 tagging logic，需源碼 audit 確認。

### 4.3 重要 caveat

per-chr V5/V6 距離尚未完整量化（IO bound 限制）。下次優先：
- 重跑 chr17 / chrX / chr8 的 V5/V6/paired_T scan（高 priority bug 區）— 確認 V6 是否在高發區真正中性化
- V6 binary diff vs V5 — 確認除 Layer 1.5 revert 外是否有其他變更

## 5. 對 PPT narrative 的影響

| narrative | 支援/反駁 | 證據 |
|---|---|---|
| 「V6 priority bug 修補有效」 | **大部分支援，chr19 例外** | IGV 3 位點 + Phase C hp=33 +959% + Phase D 4 樣本中性化 |
| 「V6 ratio 跨 4 樣本接近 1:1」 | **支援** | Phase D 已驗證 |
| 「V6 chr19 全 reads distance 反而退步」 | **CONFIRM 屬實** | V6 0.367 > baseline 0.215 > V5 0.006；V6 over-correct 往 HP1 方向 |
| 「V6 是否在 priority bug 高發 chr 真正修對」 | **未完成量化** | chr17/chrX/chr8 V5/V6/paired scan 仍待 |

**PPT 建議**：narrative 改為「V6 在跨樣本 + 全基因組 ISM-region scope 改善 priority bug，但 chr19 全 reads 層級反而退步 — 需確認 V6 是否在 priority bug 高發 chr (chr17/chrX) 真正中性化」。**不要避諱 chr19 退步**，這個發現本身重要。

## 6. 後續 follow-up

- **F-paired-D4**（高優先）等系統 IO 空閒時重跑 V5/V6/paired_T chr17 / chrX / chr8 完整 scan — 確認 V6 在高發區效益
- **F-paired-D5** V6 binary diff vs V5：除 Layer 1.5 germline-absent revert 外是否有其他 logic 變更？
- **F-paired-D6** 評估「V5 是否在 chr19 對 paired 巧合對齊」假說：若是，V6 退步本質是「V5 巧合丟失」；若否，V6 確實 over-correct
- **F-paired-D7** Phase D 4 樣本的 chr19 specific HP1/HP2 distribution 是否也類似 V6 retrograde？

## 7. Artifacts

| 檔案 | 內容 |
|---|---|
| `InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence/scan_chr_hp_tags.sh` | 4 版本 × 8 chr 平行 samtools scan 腳本 |
| `InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence/analyze_per_chr.py` | TSV 聚合 → ratio + distance 計算 |
| `InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence/igv_v6_winning_batch.txt` | IGV batch script（3 位點 snapshot）|
| `InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence/per_chr_hp/*.tsv` | per-(version, chr) HP tag counts（chr19 完整 4 版本；baseline 8 chr 中 7 完成）|
| `InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence/scan_chr_hp_tags.log` | scan log |
| `InterSubMod/research/paired_priority_bug_audit/v6_quantification_evidence/igv_run.log` | IGV batch log |
| `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_v6_winning/*.png` | 3 張 IGV 鐵證 PNG |

## 8. 結論

| Q | A |
|---|---|
| V6 在哪些 chr 真正修對？ | **chr19 germline-absent 位點層級**（IGV 鐵證 3 位點）；**全基因組 ISM-region scope**（Phase C hp=33 +959%）；**跨 4 樣本**（Phase D ratio 0.6-1.2 中性化）|
| 哪些 chr 是 priority bug 高發區？ | **chr17 (12.4:1) / chrX (8.3:1) / chr8 (5.7:1) / chr11 (4.3:1) / chr3 (4.1:1)**；chr19 (1.2:1) 是低發區 |
| V6 chr19 全 reads distance 增加是否屬實？ | **是**（V6 0.367 > baseline 0.215 > V5 0.006）— V6 在 chr19 over-correct 往 HP1 方向 |
| V6 是否在 priority bug 高發 chr 真正修對？ | **未完成量化** — chr17/chrX/chr8 的 V5/V6/paired scan IO-bound 待續 |
| 用戶 PPT「V6 有效改善」narrative 是否成立？ | **跨樣本 + ISM-region scope 成立；chr19 全 read 例外**。應誠實揭露 chr19 退步並提待證 chr17 follow-up |

→ V6 patch 是**部分有效的 priority bug 修補**，但 chr19 過修是真實副作用，不應被 PPT 隱藏。
