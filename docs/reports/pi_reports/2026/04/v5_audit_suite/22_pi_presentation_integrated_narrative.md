---
title: V5 / PON-only 修正 HP tag bias — PI 簡報用整合敘事報告
date: 2026-04-30
author: liaoyoyo2001
tags: [pi-presentation, narrative, v5, self-phasing, integrated, ppt-source]
status: presentation_source_validated
audience: PI presentation source
purpose: PPT 簡報素材 — 精準 metric map + 6 步驟敘事 + 重點位點 + 數據彙整
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/14_user_report_integrated.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/20_whole_genome_paired_audit.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/21_v5_change_boldness_audit.md
---

# V5 / PON-only 修正 HP tag bias — PI 簡報用整合敘事

## 📋 報告結構（6 大 Section 對應 PPT 章節）

| Section | PPT 章節 | 主題 |
|:-:|:-:|------|
| 1 | 開場 | 高層次問題 + 精準結論 |
| 2 | 問題定義 | 定義 17.3:1 + 驗證理論 + 分析根因 |
| 3 | 影響評估 | 為何重要：TO pure benchmark + ISM 共用 tag |
| 4 | 解釋與機制 | baseline vs V5 程式碼對比 + 圖示 |
| 5 | 改動驗證 | 多角度證據 + 矛盾與盲點 |
| 6 | 未來規劃 | 接續工作 + 與五大目標關係 |

---

# § Section 1：高層次問題與精準結論

## 1.1 一句話問題定義

> **LongPhase-TO baseline 在 somatic ALT-supporting reads metric 中出現 HP1 family : HP2 family = 17.3 : 1，顯示嚴重偏向 HP1 family。作為 truth anchor 的全 PASS sites whole-genome metric 也顯示 baseline = 2.08 : 1、paired truth = 1 : 1.08、V5 = 1.42 : 1。這代表上游 HP tag bias 會污染下游 ISM 的 HP-dependent 特徵、TO TP/FP characterization 與 subclonal interpretation，因此必須先分層修正並明確標示適用邊界。**

## 1.2 最重要結論（4 行，避免 over-claim）

```
✅ 0.93 高純度場景：全 PASS sites HP1:HP2 從 2.08:1 → 1.42:1，距 paired truth 改善 47%
✅ 不傷 caller：F1 = 0.7157 → 0.7154（噪音範圍 Δ = -0.0003）
✅ 修正邊界清楚：V2b 處理 phase graph anchor；V3F/V5 處理 tag voting / HP encoding
⚠ 在 0.6 低 purity 下 V5 conservative tagging 副作用顯著
```

## 1.3 做了哪些事

| 階段 | 動作 | 產出 |
|------|------|------|
| 1 | 確認問題 | IGV + 全基因組統計，確認 17.3:1 真實存在 |
| 2 | 找根因 | 雙層問題：phase graph anchor policy + tag 階段 somatic-first voting / HP encoding |
| 3 | 修復 | V2b/V3F/V5 core + V5+ purity helper，需分層敘述不可混成單一 V5 |
| 4 | 驗證 | 22 份 audit suite 文件、47 張 figure、20 個 TSV |
| 5 | 後續 | 0.6 simulation + WG audit + nuance + boldness audit |

## 1.4 Metric map（簡報前必先鎖定 denominator）

| Metric | Denominator | Baseline | V5 / V5+ | Paired truth | 可宣稱 |
|--------|-------------|:-:|:-:|:-:|------|
| Somatic ALT-supporting reads bias | whole-genome somatic ALT-supporting read subset | **17.3 : 1** | 不用 1.42 直接對比 | 接近 1:1 | 用來定義 baseline 問題嚴重度 |
| Whole-genome PASS sites aggregate | all PASS sites / total ALT reads | **2.08 : 1** | **1.42 : 1** | **1 : 1.08** | V5 在 0.93 場景讓 aggregate HP balance 更接近 paired truth |
| Site-by-site dominant HP concordance | PASS sites with PA_93 dominant HP1/HP2 | 約 random level | 與 baseline 接近，部分場景略差 | PA_93 | 不能宣稱 V5 全面提升 per-site correctness |
| Caller F1 | SEQC2 truth / ClairS-TO calls | 0.7157 | 0.7154 | N/A | V5 不改 caller，F1 不變是預期 |
| 0.6 purity aggregate | t30_n20 低 purity BAM | baseline 自然減弱 | V5 HP33 增加，可能更偏離 PA_93 | 缺 0.6 paired truth | 只能說 conservative tagging 明顯，不能說 0.6 下全面修復 |

---

# § Section 2：觀察到的問題

## 2.1 變數定義與名詞解釋

### 2.1.1 HP tag 5 種編碼

LongPhase-TO 的 haplotag 在 BAM 寫入 `HP:i:N`：

| HP:i: 值 | 含義 | 規範 (phase.md) |
|:--:|------|------|
| **0** | untagged（無 phasing 證據）| HAPLOTYPE_UNDEFINED |
| **1** | germline → HP1 | HP1 family |
| **2** | germline → HP2 | HP2 family |
| **11** | hp1-1：somatic from HP1 | HP1 family |
| **21** | hp2-1：somatic from HP2 | HP2 family |
| **33** | hp3：ambiguous somatic | 不歸 family |

→ **HP1 family** = {HP:i:1, HP:i:11}；**HP2 family** = {HP:i:2, HP:i:21}

### 2.1.2 Phasing 與 Tag 兩階段

```
[Stage 1: phase] (PhasingProcess.cpp + PhasingGraph.cpp)
  Input: ClairS-TO snv.vcf + tumor BAM
  Output: phased.vcf (含 GT/PS/GT2/PS2/GT3/PS3 標記)
                      ↓
[Stage 2: haplotag] (HaplotagProcess.cpp)
  Input: phased.vcf + tumor BAM
  Output: tumor_tagged.bam (含 HP:i: 標記)
```

→ Phase 在 VCF 上分配 GT；Tag 在 BAM 上分配 HP。**兩階段獨立**。

## 2.2 問題點：HP1 family : HP2 family = 17.3 : 1

### 2.2.1 全基因組 ALT-supporting reads 統計（17.3:1 metric 出處）

從 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` §3：

| 指標 | Baseline @ 0.93 |
|------|:-:|
| Somatic ALT reads → HP1 | **614,000** |
| Somatic ALT reads → HP2 | **35,500** |
| **HP1 : HP2 ratio** | **17.3 : 1** ❌ |

→ **這是 specific subset metric**（whole-genome somatic ALT reads aggregate）。

### 2.2.2 全 PASS sites WG 統計（從 20 號 §2.1）

| Sample | total_alt | HP1 | HP2 | HP1:HP2 |
|--------|:-:|:-:|:-:|:-:|
| Baseline @ 0.93 | 2,344,492 | 1,497,129 | 719,367 | **2.08 : 1** |
| **Paired truth (PA_93)** | 2,344,492 | 325,505 | 350,723 | **1 : 1.08** ✓ |
| V5 @ 0.93 (修復後) | 2,344,492 | 1,218,134 | 856,769 | **1.42 : 1** ✅ |

→ **PA truth ≈ 1:1 證實理論預期**；baseline 偏離 PA 0.352；V5 偏離 0.186（**改善 -47%**）。

![WG 統計圖](figures/20_whole_genome/wg_summary.png)

### 2.2.3 Per-site 極端例子（chr19 self-phasing cluster）

從 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/05_per_site_improvement.md`：

| 位點 | Baseline ratio | 註解 |
|------|:-:|------|
| **SP1** chr19:17565944 | **113:0** (∞) | 完全極端 self-phasing |
| **SP2** chr19:12452332 | **109:1** | 109× imbalance |
| **SP3** chr19:12467180 | **108:0** (∞) | 完全極端 |

→ 個別位點可達 **109× 以上**，遠超平均 17.3×。

## 2.3 驗證理論：HP1:HP2 應為 1:1（paired 證實）

### 2.3.1 理論依據

**Tumor diploid 假設**：每條 read 來自 maternal (HP1) 或 paternal (HP2) 機率 50:50 → ratio ≈ 1:1。

### 2.3.2 Paired 實證（Section 2.2.2 數據）

PA_93（paired-tagged tumor BAM, normal-anchored phasing）：
- HP1 = 325,505 / HP2 = 350,723
- ratio = **1 : 1.08** ≈ 1:1 ✓

→ **paired 數據實證理論預期成立**，baseline 的 17.3:1 是**真實 bug**。

## 2.4 分析根因：不是 caller bug；需區分 phase graph policy 與 tag writing

### 2.4.1 證據 ① — V2b/V5 backbone 的 VCF GT 分布一致

從 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md` §3.1：

| GT_class | Baseline % | V2b/V5 % | Δpp |
|----------|:--:|:--:|:--:|
| Germline_Het | 13.09 | 12.91 | −0.18 |
| **Somatic_NoLOH (`0\|0`)** | **42.99** | **42.99** | **0.00** |
| **Somatic_in_LOH** | **24.19** | **24.19** | **0.00** |

→ **PASS somatic GT 100% 一致**（21,304 / 11,983 個 site 1:1 對應）表示在 V2b/V5 backbone 下，主要差異不在 final GT distribution，而在 BAM HP tag 寫入與讀段分配。

### 2.4.2 結論：問題鏈是雙層，不是單純「phasing 沒問題」

```
Caller 階段（ClairS-TO calls）：✓ 不變
                            ↓
Phase graph policy：⚠ baseline 允許 ClairS-TO PASS somatic 參與 phasing anchor
                            ↓ V2b: --pon-only-phasing / convertNonGermlineToSomatic()
VCF GT distribution：✓ V2b/V5 backbone 下 PASS somatic GT 分布一致
                            ↓
Tag 階段（BAM HP 寫入）：❌ baseline getVote() somatic-first + HP encoding bug
                            ↓ V3F/V5: germline-first + Layer 1.5 fallback + explicit HP33
BAM HP tag：✅ aggregate balance 在 0.93 場景更接近 paired truth
```

→ 17.3:1 bias 是 **phase graph anchor policy + tag voting / encoding** 疊加的結果。對 PI 報告時不要簡化成「phasing 完全沒問題」或「V5 一次修好所有 self-phasing」。

---

# § Section 3：為何這個問題重要

## 3.1 影響 #1：TO pure 數據是核心 benchmark

InterSubMod 的研究主要在 **tumor-only (TO) 模式**：
- 沒有 paired normal 可用
- 完全依賴 longphase-to 的 phasing + tag 結果
- TP / FP 觀察、subclonal analysis 都基於 tumor-only tag

→ 如果 tag 結果偏差 17.3:1，**直接使用 HP tag 的 TO 觀察都可能帶 bias**；不依賴 HP 的 caller-level 指標則需分開解讀。

## 3.2 影響 #2：ISM 等下游分析共用 tag

```
ISM 分析流程：
  Tagged BAM (V5 / Baseline)
    ↓
  按 HP1/HP2 分組 reads
    ↓
  比較甲基化 pattern (allele-specific methylation)
    ↓
  推論 subclonal structure / epigenetic heterogeneity
```

→ ISM 第一步就用 HP tag。**Tag 錯 = ISM 錯**。

### 3.2.1 量化影響範例

如果 baseline 把 94.6% somatic reads 標 HP1（17.3:1 = 94.6:5.4）：
- 「HP1 上的 somatic 甲基化率」會被 614K reads 強烈影響
- 「HP2 上的 somatic 甲基化率」只有 35K reads 支持，統計力極弱
- **HP allele-specific methylation analysis 失準**

## 3.3 影響 #3：與五大研究目標的關聯

| 研究目標 | HP tag 依賴層級 | 影響說法 |
|---------|:----:|------|
| 目標 1: Subclonal structure detection | ✅ 直接依賴 | HP 分組錯會直接污染 clone / subclone 推論 |
| 目標 2: Allele-specific methylation | ✅ 直接依賴 | ASM 需要可信 HP1/HP2 分組 |
| 目標 3: TP/FP characterization (TO mode) | ✅ 直接依賴 | TO TP/FP 的 HP 分布與 HP-dependent feature 需重看 |
| 目標 4: Cross-platform consistency | ⚠ 部分依賴 | 若比較的是 HP-dependent 指標才直接受影響 |
| 目標 5: Clinical interpretation | ⚠ 間接依賴 | 需先經過穩定 feature / model 才影響臨床敘事 |

→ **3/5 目標直接依賴 HP tag，2/5 間接或部分依賴**。修復 tag 是上游基礎建設工作，但不能寫成所有研究結論都會等幅翻轉。

## 3.4 結論：修復 tag 的優先級 = P0

```
若不修：
  ❌ HP-dependent TO 研究結論帶 17.3:1 或 2.08:1 類型偏差
  ❌ ISM 的 HP_Ratio / HPFineNGroups / HPMergedDelta 等欄位失準
  ❌ 直接依賴 HP tag 的研究目標基線不穩
若修復：
  ✅ 0.93 高純度場景下 HP aggregate balance 更接近 paired truth
  ✅ 下游可用 V5 / V5+ tag 作為新 baseline，但需保留 0.6 conservative caveat
  ✅ V5 audit suite 22 份文件提供 traceability
```

---

# § Section 4：解釋問題與發生原因

## 4.1 Baseline 的設計邏輯

### 4.1.1 Phase 階段（baseline）

baseline 接受所有 ClairS-TO PASS variants 進入 phasing graph：
- PON-confirmed germline 進 graph
- ClairS-TO PASS somatic **也**進 graph
- → graph 用所有共現 reads 算 edge weight

**問題**：同一 tumor clone 的 reads 共享相同 ALT pattern → graph 把 somatic sites 推進同一 phase block。

### 4.1.2 Tag 階段（baseline `getVote()`）

baseline 的 priority for-loop：
```cpp
// Baseline (錯誤順序)：
variantKeys = {
    {HAPLOTYPE1_1, HAPLOTYPE2_1},   // ← 第 1 優先（somatic）❌
    {HAPLOTYPE3, HAPLOTYPE2_1},      // ← 第 2 優先
    {HAPLOTYPE1, HAPLOTYPE2}         // ← 第 3 優先（germline）❌
};
for (const auto& pair : variantKeys) {
    if (countMap[key1] > 0 || countMap[key2] > 0) {
        hpResult = haplotypeBase[winner];  // 直接給 11/21/33
        break;
    }
}
```

**問題**：read 涵蓋任何 somatic site → 第 1 優先立刻決定 → germline 投票被無條件覆蓋。

## 4.2 V5 的修法對比

### 4.2.1 V2b Phase graph policy：`--pon-only-phasing` flag

```cpp
// V2b / V5 backbone (PhasingProcess.cpp:154-157)：
if (params.ponOnlyPhasing) {
    vGraph->convertNonGermlineToSomatic();
    // → 把所有非 PON variants 標為 SOMATIC origin
    // → phasing graph 只用 PON germline 當 anchor
}
```

### 4.2.2 V5 Tag 階段：`getVote()` 三層決策

```cpp
// V5 (HaplotagProcess.cpp:512-563)：

// Layer 1: Germline first（優先級反轉）
int germlineResult = 0;
if (germlineHP1 > 0 || germlineHP2 > 0) {
    germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
}
// Layer 1.5: Somatic fallback（只有 germline 為 0 才用）
else if (somaticHP1 > 0 || somaticHP2 > 0) {
    germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;
}

// Layer 2: Encoding（與 voting 拆分）
if (somaticTotal > 0) {
    if      (germlineResult == 1) hpResult = 11;  // HP1 + somatic
    else if (germlineResult == 2) hpResult = 21;  // HP2 + somatic
    else                          hpResult = 33;  // ambiguous
} else {
    hpResult = germlineResult;  // 0/1/2
}
```

### 4.2.3 V5 修法重點

> 用戶定義的修法核心：
> **「依據 HP1 與 HP2 的量然後看是否經過 somatic 這樣去定義」**
> 而不是 baseline 的「依據 HP1-1 與 HP2-1 的量直接定義」

| 對比 | Baseline | V5 |
|------|------|------|
| **投票主體** | somatic-traceable (HP1_1/HP2_1) 為 primary | germline (HP1/HP2) 為 primary |
| **Encoding 觸發** | 直接給 11/21 | 先決 germlineResult，再加 somatic encoding |
| **設計理念** | somatic 投票主導 | germline 主導，somatic 作 annotation |

## 4.3 圖示化解釋

### 4.3.1 17.3:1 self-phasing 機制圖

![Self-phasing 機制](figures/01_code_diff/fig01d_somatic_bias_explanation.png)

→ 6-panel：Expected (1:1) vs Observed (17.3:1) + 機制 + 全基因組量化 + Per-site 極端例子 + V2b/V3F/V5 修正效果

### 4.3.2 getVote 三層決策樹

![getVote 三層決策](figures/13_phase_vs_tag_algo/figB_tag_getVote_decision.png)

→ 對比 baseline 單體 if-else (Layer 1 only) vs V5 三層拆分

### 4.3.3 具體例子（3 reads × 5 variants in LOH region）

![具體例子](figures/13_phase_vs_tag_algo/figC_concrete_example.png)

→ Baseline：Read B/C lost (HP:i:0)；V5：Read B/C 用 Layer 1.5 fallback 救起為 HP:i:11

### 4.3.4 SE 軟工視角

![SE 視角](figures/15_se_perspective/fig15a_se_dimensions.png)

→ 5 SE 維度：SoC / SRP / OCP / Defensive / Fail-Safe

## 4.4 修法清單（分 V2b / V3F / V5 core / V5+）

| 層級 | 修法 | 程式碼位置 | 性質 | PI 報告口徑 |
|------|------|:-:|------|------|
| V2b | `--pon-only-phasing` flag + non-PON origin 轉 somatic | PhasingProcess.cpp:154-157 + Graph.cpp:1139 | 大膽改 | 解 phase graph anchor policy；不是 V5 單獨完成 |
| V3F | `getVote()` 改為 germline-first | HaplotagProcess.cpp:512-563 | 必須改 | 修 tag voting priority |
| V3F | HP:i:33 enum / integer literal 比對 | HaplotagProcess.cpp:556-559 | 必須改 | 讓 ambiguous HP33 可被正確寫出 |
| V3F / V5 | `countSNP/INDEL` UNDEFINED guard | HaplotagProcess.cpp:484-510 | 必須改 | 防禦性修補，避免 undefined haplotype 誤進投票 |
| V5 core | Layer 1.5 somatic fallback + conservative HP33 | HaplotagProcess.cpp:512-563 | 條件性設計 | 0.93 場景改善 aggregate balance；0.6 下需保留 HP33 副作用 caveat |
| V5+ | `collectPloidyRatio()` purity helper | PhasingGraph.cpp:1147-1175 | completeness fix | 本 session 後補；建議與 V5 core 分開命名 |

詳見 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/21_v5_change_boldness_audit.md`。簡報時建議使用 **V5 core** 指 haplotag 修補，用 **V5+** 指包含 `collectPloidyRatio()` 的後補版本。

---

# § Section 5：改動驗證與結論

## 5.1 多角度證據

### 5.1.1 證據 #1 — 先定義 17.3:1 問題，再用 WG audit 量化改善

```
Somatic ALT-supporting reads subset:
  Baseline: HP1=614K, HP2=35.5K → 17.3:1 ❌

Whole-genome PASS sites / total ALT reads:
  Baseline: HP1:HP2 = 2.08:1
  V5_93:    HP1:HP2 = 1.42:1
  PA_93:    HP1:HP2 = 1:1.08
  distance to PA: 0.352 → 0.186（改善 47%）✅
```

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` + `20_whole_genome_paired_audit.md`

### 5.1.2 證據 #2 — 全 PASS sites WG audit

| Sample | HP1:HP2 | distance to PA |
|--------|:-:|:-:|
| PA_93 | 1:1.08 | 0 (truth) |
| BL_93 | 2.08:1 | 0.352 |
| **V5_93** | **1.42:1** | **0.186** ← 改善 −47% |

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/20_whole_genome_paired_audit.md`

### 5.1.3 證據 #3 — VCF GT 分布一致（差異主要落在 tag 階段）

| | PASS Somatic_NoLOH | PASS Somatic_in_LOH |
|--|:--:|:--:|
| Baseline | 21,304 | 11,983 |
| V5 (V2b backbone) | 21,304 | 11,983 |
| Δ | **0** | **0** |

→ 證實在 V2b/V5 backbone 下 PASS somatic GT distribution 一致；這支持「後續可把主要差異聚焦在 tag writing」，但不應倒推成 baseline phase graph policy 完全沒問題。

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md`

### 5.1.4 證據 #4 — F1 不受影響（不傷 caller）

| 版本 | F1 (vs SEQC2 truth) | Δ |
|------|:-:|:-:|
| Baseline | 0.7157 | (基準) |
| V3-Fixed | 0.7153 | −0.0004 |
| **V5** | **0.7154** | **−0.0003 (噪音)** |

→ V5 修法不影響 ClairS-TO calling，純粹 phasing/tag 修復。

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md`

### 5.1.5 證據 #5 — Sanity check 5/5 PASS

| 檢查項 | 結果 |
|--------|:-:|
| 守恆律 A：Δ_HP33 + Δ_(HP11+HP21) = 0 | 15/15 PASS |
| 守恆律 B：germline (HP1, HP2) 不變 | 15/15 PASS |
| Layer 1.5 預期 1：新 HP11/HP21 在 V3F 中是 HP33 | 15/15 PASS |
| Layer 1.5 預期 2：無 germline → HP33 違規 | 0 violations |
| Untagged → directional：V5 不強 tag 原 HP=0 | 0 violations |

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md`

### 5.1.6 證據 #6 — IGV 視覺案例（4 版本對比）

關鍵位點 IGV 截圖（紅=HP1, 藍=HP2, 紫=HP33, 灰=untagged）：

#### **D_SP1** chr19:17565944 — Self-phasing extreme

![SP1 IGV](figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png)

| 版本 | reads HP 分布 | ratio |
|------|------|:-:|
| Baseline @ 0.93 | 全紅 (HP1) | **113:0** ❌ |
| V5 @ 0.93 | 全藍 (HP2) | 0:113 (翻轉) |
| Baseline @ 0.6 | 大部分藍 | 1:70 |
| V5 @ 0.6 | 大部分藍 + 紫 | 0:71 |

#### **A_TP04** chr16:35118902 — V5 唯一強改善

![TP04 IGV](figures/igv_v5_audit/by_HP_4ver/A_TP04_chr16_35118902.png)

| 版本 | exact-site | paired conditional | L4 |
|------|:-:|:-:|:-:|
| Baseline | 2.7:1 | (基準) | (基準) |
| V5 | 1:6.9 | **+0.9737** | 略輸 |

→ **V5 唯一 strong improve 案例**。

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/19_per_site_nuance_audit.md`

## 5.2 為何 V5 修法有效

### 5.2.1 邏輯閉環檢查

```
[Phase 階段]
  --pon-only-phasing
  → convertNonGermlineToSomatic()
  → graph anchor only PON germline
  → 非 PON variants 不再作為等價 germline anchor
  → 降低 somatic-as-anchor 造成的 phase graph bias
  ↓
[Tag 階段]
  Layer 1: germline first
  Layer 1.5: somatic fallback (LOH 救起)
  Layer 2: encoding 拆分
  → read 投票公平 (germline primary, somatic secondary)
  → HP33 顯式標 ambiguous (conservative)
  ↓
[結果]
  HP1:HP2 在 0.93 aggregate metric 更接近 paired truth ✓
  Sanity check 全 PASS ✓
  F1 不變 ✓
```

### 5.2.2 矛盾與盲點

從 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/19_per_site_nuance_audit.md` + `20_whole_genome_paired_audit.md`：

| 矛盾 | 解讀 |
|------|------|
| Aggregate 改善 −47% vs Site-level match% −0.5pp | aggregate 信度高（large-N），site-level 是 conditional metric |
| 13 sites 中 1 strong / 2 regression / 5 tie / 4 NA | V5 不是全面改善，是「特定場景改善」 |
| V5_06 在 0.6 反而偏離 PA +192% | low purity 下 conservative 副作用 |
| HP33% V5_06 11% (BL_06 1.7%) | conservative-aware metric 缺失 |

### 5.2.3 是否符合論文設計

從 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/17_design_consistency_check.md`：

| 設計理念 | V5 是否違反 |
|---------|:--:|
| 6 步驟流程 | ❌ 不違反 |
| Two-Pass 高純度策略 | ❌ 不違反（V5 重用 `convertNonGermlineToSomatic`）|
| GT/GT2/GT3 規範 | ❌ 不違反 |
| HP:i: 5 種編碼 | ❌ 不違反 |
| Phase Set 機制 | ❌ 不違反 |

→ **V5 完全對齊原始 LongPhase-TO 設計理念**，不是發明新邏輯。

## 5.3 業界相關方法支持

| 方法 | 來源 | 對應 V5 哪部分 |
|------|------|------|
| LongPhase v1.7 paired mode | LongPhase paper | germline-anchored phasing |
| ClairS-TO PoN-aware filtering | ClairS-TO docs | --pon-only-phasing flag |
| WhatsHap haplotag (germline-first) | WhatsHap | getVote() Layer 1 = germline first |
| GATK PhaseByTransmission | GATK best practices | 主→輔 priority pattern |

→ V5 的 germline-first 設計**符合常見 phasing / haplotagging 的主輔 priority 慣例**。簡報時建議避免寫成「業界標準已證明本實作正確」，改寫成「LongPhase/WhatsHap/GATK 等相鄰設計都支持先使用可信 germline anchor，再處理 somatic annotation」。

## 5.4 結論

**V5 / V5+ 修法是有效、可信、但有適用邊界的工程改進**：
- ✅ 0.93 高純度場景下，全 PASS sites aggregate 從 2.08:1 → 1.42:1（距 paired truth 改善 47%）
- ✅ Sanity check 全 PASS
- ✅ F1 不變
- ✅ 對齊論文 + 知識庫 + README
- ⚠ 0.6 下 conservative 副作用需後續優化

## 5.5 可宣稱 / 需保留 / 不可宣稱

| 類別 | 內容 |
|------|------|
| 可宣稱 | Baseline 有嚴重 HP tag bias；V2b/V3F/V5 core 在 0.93 場景讓 aggregate HP balance 更接近 paired truth；caller F1 不變符合預期 |
| 需保留 caveat | 0.6 低 purity 下 HP33 conservative tagging 明顯；site-by-site paired concordance 不全面提升；V5+ purity helper 是後補 completeness fix |
| 不可宣稱 | V5 單獨修復所有 self-phasing；17.3:1 直接變成 1.42:1 是同一 denominator；V5 全面提升 per-site correctness；所有五大研究目標等幅受益 |

---

# § Section 6：未來目標與規劃

## 6.1 修復後立即可接續的工作

### 6.1.1 高優先（基於 V5 / V5+ 修正後的 BAM）

| 任務 | 描述 | 預估時間 |
|------|------|:--:|
| **重跑 ISM 分析** | 用 V5 tagged BAM 重做 allele-specific methylation | 1-2 day |
| **TP/FP characterization** | 用 V5 BAM 看 TP / FP 在 HP1/HP2 上的真實分布 | 0.5-1 day |
| **更新 v5_audit_suite §02-08** | 加入本次新發現的 nuance（19/20/21 號）| 0.5 day |

### 6.1.2 中優先（V5+ 改進）

| 任務 | 描述 | 預估時間 |
|------|------|:--:|
| **logging syncOrigins erase** | 補強 `syncPhasingResultOrigins` 的 logging | 1 hr |
| **0.6 paired truth 生成** | 跑 t30_n20 paired-mode longphase | 1-2 hr |
| **HP33-aware metric 設計** | 設計 conservative-aware concordance | 1-2 day |
| **重訓練 purity polynomial** | 加 PON-only mode 的 q1/q3 分布 | 1-2 day |

### 6.1.3 低優先（長期）

| 任務 | 描述 |
|------|------|
| Trio sequencing | 取得 read-level hap truth |
| Cross-sample 驗證 | COLO829 / DORADO / H2009 等樣本重做 V5 audit |
| V5 vs SAVANA / Wakhan | 與其他 caller 整合 |

## 6.2 與五大研究目標關係（研究發展樹）

```
                    [V2b/V3F/V5 修正 HP tag bias]
                            │
         ┌──────────────────┼──────────────────┐
         ▼                  ▼                  ▼
   [目標 1: Subclonal]  [目標 2: ASM]   [目標 3: TP/FP]
   structure detection   methylation     characterization
         │                  │                  │
         ├─ 重跑 ISM       ├─ HP-based ASM   ├─ TO mode FP root cause
         ├─ Sample Het     ├─ Cross-platform  ├─ Low-AF rescue
         └─ subclonal mutation pattern ←──────┴─────────────┘
                            │
                ┌───────────┴───────────┐
                ▼                       ▼
         [目標 4: Cross-platform]  [目標 5: Clinical]
         consistency (DORADO/5kHz)  interpretation
```

→ 這條修補鏈是 **HP-dependent 研究目標的基石**；對非 HP-dependent 結論則主要提供版本標註與風險隔離。

## 6.3 已有初步發現的後續研究方向

從 audit suite 既有 + 本 session 新發現：

| 發現 | 報告 | 後續研究 |
|------|------|---------|
| **0.6 下 V5 conservative 副作用** | 09 + 20 號 | 設計 purity-aware HP33 threshold |
| **13 sites 改善等級分類**（1 strong / 2 reg / 5 tie / 4 NA） | 19 號 | 找改善的真實場景（按 AF/cov/LOH 分層） |
| **V5 雙重修法**（PON-only + getVote 順序） | 16 號 | 是否有更精準的 single-fix 設計 |
| **Purity calculator 在 PON-only mode 失效** | 18 號 | 已修，待重訓練 polynomial |
| **paired truth 71% reads untagged** | 20 號 | trio sequencing 補強 truth |

## 6.4 整體進展圖

```
[2026-04-10] commit 8b8c1fd: --pon-only-phasing flag
                     ↓
[2026-04-15] commit 41ff147: getVote two-layer fix
                     ↓
[2026-04-22] commit 380e8d2: countINDELHaplotype guard
                     ↓
[2026-04-27] V5 working tree: countSNPHaplotype guard
                     ↓
[2026-04-27 ~ 04-30] Audit suite 22 份報告
                     ↓
[2026-04-29] V5+ collectPloidyRatio helper（本 session）
                     ↓
[次步] 重跑 HP-dependent ISM + TP/FP characterization with V5 / V5+ BAM
                     ↓
[長期] 直接依賴 HP tag 的研究目標基於 V5 / V5+ baseline 重新展開
```

---

# § PPT Slide 對應素材清單

## PPT 章節 1：開場（3 slides）

- Slide 1: 一句話結論（Section 1.1）
- Slide 2: 17.3:1 問題定義 + 2.08:1 → 1.42:1 三大成果（Section 1.2）
- Slide 3: 整體流程（Section 1.3）

## PPT 章節 2：問題定義（4 slides）

- Slide 4: HP tag 5 種編碼定義（Section 2.1）
- Slide 5: Phase vs Tag 兩階段（Section 2.1.2）
  - 圖：fig01c_pon_phase_tag_comparison.png
- Slide 6: 17.3:1 全基因組數據（Section 2.2）
  - 圖：fig01d_somatic_bias_explanation.png
- Slide 7: SP1/SP2/SP3 極端 per-site 例子（Section 2.2.3）

## PPT 章節 3：影響評估（2 slides）

- Slide 8: TO pure 與 ISM 影響（Section 3.1-3.2）
- Slide 9: 五大目標關聯（Section 3.3）

## PPT 章節 4：機制與修法（5 slides）

- Slide 10: Baseline getVote 順序錯誤（Section 4.1.2）
- Slide 11: V5 getVote 三層決策（Section 4.2.2）
  - 圖：figB_tag_getVote_decision.png
- Slide 12: 具體例子 3 reads × 5 variants（Section 4.3.3）
  - 圖：figC_concrete_example.png
- Slide 13: V5 PON-only flag（Section 4.2.1）
- Slide 14: 分層修法清單（Section 4.4）

## PPT 章節 5：驗證證據（5 slides）

- Slide 15: 17.3:1 問題 subset + 2.08:1 → 1.42:1 WG 數據（Section 5.1.1-5.1.2）
  - 圖：wg_summary.png
- Slide 16: VCF GT 100% 一致 + F1 不變（Section 5.1.3-5.1.4）
- Slide 17: Sanity check 5/5 PASS（Section 5.1.5）
- Slide 18: IGV 視覺案例（SP1 + TP04）（Section 5.1.6）
- Slide 19: 矛盾與業界對齊（Section 5.2-5.3）

## PPT 章節 6：未來規劃（3 slides）

- Slide 20: 立即接續工作（Section 6.1）
- Slide 21: 五大目標關係樹（Section 6.2）
- Slide 22: 整體進展時間軸（Section 6.4）

---

# § PPT 重點位點與數據速查表

## 必引用的位點

| 用途 | 位點 | 為何引用 |
|------|------|---------|
| **17.3:1 self-phasing 例子** | SP1 chr19:17565944 (113:0) | 最極端，視覺強 |
| **V5 強改善案例** | A_TP04 chr16:35118902 | 唯一 strong improve |
| **0.6 self-phasing 衰減** | SP2 chr19:12452332 (109:1 → 1:41) | 衰減 60% 顯著 |
| **V5 反向 regression** | C_V5max3 chr19:7405500 | 揭露 V5 不是全面改善 |

## 必引用的數據

| 用途 | 數據 | 出處 |
|------|------|------|
| 整體修復成果 | 問題 subset 17.3:1；WG aggregate 2.08:1 → 1.42:1 | 10 號 + 20 號 |
| 修復幅度 | distance −47% | 20 號 §2.2 |
| 不傷 caller | F1 = 0.7157 → 0.7154 (Δ −0.0003) | 07 號 |
| Sanity 5/5 | 守恆律 + Layer 1.5 全 PASS | 06 號 |
| V5+ purity helper 修正 | PON-only purity estimate 0 → 0.634（仍需重訓練/驗證） | 18 號 |

## 必引用的圖

| 圖 | 路徑 |
|----|------|
| 17.3:1 機制 6-panel | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01d_somatic_bias_explanation.png` |
| getVote 三層決策 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/13_phase_vs_tag_algo/figB_tag_getVote_decision.png` |
| 具體例子 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/13_phase_vs_tag_algo/figC_concrete_example.png` |
| WG 統計 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/20_whole_genome/wg_summary.png` |
| SP1 IGV (4ver) | `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png` |
| TP04 IGV (4ver) | `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/A_TP04_chr16_35118902.png` |
| SP1 IGV (paired+LOH/GE) | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/09_purity06/igv_with_paired_loh/SP1_chr19_17565944.png` |
| 影響分類矩陣 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/14_impact_classification/fig14a_impact_matrix.png` |
| SE 5 維度 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/15_se_perspective/fig15a_se_dimensions.png` |

---

# § 跨檔索引（按 PPT 章節）

| PPT 章節 | 對應 audit suite |
|---------|---------|
| Section 1 開場 | 14 整合報告、08 synthesis、20 WG audit |
| Section 2 問題定義 | 10 somatic bias、12 GT distribution、13 phase vs tag、20 WG |
| Section 3 影響評估 | 14 整合報告（自有）、目標關聯需另寫 |
| Section 4 機制與修法 | 01 code diff、11 issues inventory、13 algorithm detail、15 SE、16 baseline subgenotype |
| Section 5 驗證證據 | 02-08 (4 agents)、19 nuance、20 WG、17 design consistency、18 purity fix |
| Section 6 未來規劃 | 21 boldness audit、CURRENT_FOCUS、五大目標 manifest |

完整 22 份 audit suite 索引：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md`
