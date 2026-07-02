---
title: ISM 判定軸的 read-set provenance 拆解（tumor-only vs 需 tumor/normal 對照）
date: 2026-06-16
type: observation
scope: single-sample HCC1395, paired, ±5000bp window, whole-genome TP SNV
status: in_progress / partial（單樣本，需 COLO829 等跨樣本復現）
tier: L2（第一手讀回全基因組輸出 + 程式碼 trace；軸門檻為 first-pass heuristic）
build_branch: feat/summary-nreadsvalid @ 990175d
data_sources: docs/methodology/_assets/20260616_readset_provenance/readset_stats.json, output/_wg_d1_unified/significance_summary.csv
related:
  - docs/methodology/20260616_structure_label_verdict_false_negative_audit_01.md
  - docs/methodology/20260616_significance_verdict_methodology_explainer_for_ai_01.md
  - feedback_asm_allele_axis_baseline_confound (memory)
---

> ⚠️ **PARTIAL / 單樣本觀察** — HCC1395 paired，±5000bp window，全基因組 30,490 個 TP SNV。
> 軸觸發門檻（AUC≥0.7、|Δβ|≥0.1）為 first-pass heuristic，非驗證過的最佳邊界。
> 跨樣本（COLO829 等）復現待補。**勿當判別性 filter 結論**。

# ISM 判定軸的 read-set provenance 拆解

## §0 TL;DR（Data-Showcase 框架）

本文回答 ask (1)：**目前「Significant」判定裡，哪些訊號用 tumor reads 自己就能偵測，哪些需要 tumor/normal 對照才驗證得出？** 三個發現：

1. **判定骨幹是「混池（all-pool）」軸，不是純 tumor 軸。** 主距離矩陣 / clustering / GlobalTest / PERMANOVA 全在 **tumor+normal 混池** 上跑（`read_list = get_reads()`，含全部 reads）。只有 3 個軸（within-HP、HP_AUC_Tumor、Tumor_HP_Delta）是純 tumor。
2. **只有 48.0% 的 Significant 有純 tumor 軸佐證；55.5% 有 tumor-vs-normal 對照軸；28.6% 兩者皆無**（純靠混池軸，最不可信）。
3. **混池「有結構」分不出三件事**：真 tumor 子結構 / germline HP phasing（兩邊都有）/ 就只是把 tumor 和 normal reads 分開（sampleOrphan confound）。→ 判定本身**無法自證 somatic**，需 C 類對照軸補。

---

## §1 前置確認：Δβ 計算範圍 — 有界（±5000bp），非無限延伸 ✅

所有標籤 Δβ（germline_asm / somatic_residual / tumor_hp / normal_hp / allele / hp_merged / subclone / within-HP）+ normal_baseline 都吃**同一個 in-window 甲基矩陣** `meth_mat.raw_matrix`。

該矩陣的 CpG 在解析時就被位置過濾（`src/core/MethylationParser.cpp` 行 152-155 / 190-193）：
```cpp
int ref_offset = ref_pos_0based - ref_start_pos;          // ref_start_pos = region_start
if (ref_offset >= 0 && (ref_offset+1) < ref_seq.size() && ref_seq[ref_offset]=='C' && ...)
    calls.emplace_back(...)                               // 只有 window 內才收
```
`ref_seq` 長度 = `region_end − region_start` = SNV ± window_size。**一條 50kb ONT read 只貢獻落在 window 內的 CpG，window 外的 CpG `ref_offset` 超界 → 不進矩陣。**

- 本次 run window = **±5000bp**（NumCpGs median=76、min=16、`<10 CpG=0%`，與既有 ±5000 特性吻合）。
- 結論：Δβ 由建構保證有界，**「無限延伸」風險不存在**。

---

## §2 方法：判定軸依 read-set 分三類

| 類別 | 軸 | 來源欄位 | read-set | 程式位置 |
|---|---|---|---|---|
| **A 純 TUMOR-only** | within-HP multigroup | `WithinHP_CleanMultigroup` | tumor reads（`if(!is_tumor)continue`）| RegionProcessor.cpp:1169 |
| | HP_AUC_Tumor≥0.7 | `HP_AUC_Tumor` | tumor reads（idx_t）| :2237 |
| | Tumor_HP_Delta≥0.1 | `Tumor_HP_Valid`+`Tumor_HP_Delta` | tumor reads | :970 |
| **B ALL-pool** | cluster_match (Strong) | `VerificationClass_Legacy∈{Strong,Subclone}` | tumor+normal 混池 | :848/:866 |
| | label-allele / label-HP PERMANOVA | `LabelAllele/HPPermanovaValid`+p | 混池 | StructureTest |
| | AlleleSig | `AlleleSig` | 混池（ALT vs REF）| :2450 |
| | HP_AUC_All≥0.7 | `HP_AUC_All` | 混池（idx_a）| :2238 |
| **C TUMOR-vs-NORMAL 對照** | HP_Residual_Sig | `HP_Residual_Sig`（tumor HP − normal HP）| 需 normal | :1006 區 |
| | SomaticResidualDbeta_Sig | `SomaticResidualDbeta_Sig`（扣 normal baseline）| 需 normal | :1040 |
| | GermlineAsmDbeta_Sig | `GermlineAsmDbeta_Sig`（germline 等位特異）| 需 normal | :1066 |

**判據**：`tumor_only = A 任一觸發`；`contrast = C 任一觸發`；`only_pool = 既無 A 也無 C（純靠 B）`。

---

## §3 結果：各軸在 28,327 Significant 上的觸發數

### A. 純 TUMOR-only 軸（不需 normal）
| 軸 | 觸發數 | 佔 Significant |
|---|---:|---:|
| within-HP multigroup | 8,077 | 28.5% |
| HP_AUC_Tumor≥0.7 | 7,871 | 27.8% |
| Tumor_HP_Delta≥0.1 | 2,513 | 8.9% |

### B. ALL-pool 軸（混池；偵測「有結構」但無法自證 somatic）
| 軸 | 觸發數 | 佔 Significant |
|---|---:|---:|
| **cluster_match (Strong)** | **22,310** | **78.8%** |
| label-allele PERMANOVA | 27,009 | 95.3% |
| label-HP PERMANOVA | 24,875 | 87.8% |
| AlleleSig | 16,339 | 57.7% |
| HP_AUC_All≥0.7 | 7,345 | 25.9% |

### C. TUMOR-vs-NORMAL 對照軸（需 normal）
| 軸 | 觸發數 | 佔 Significant |
|---|---:|---:|
| HP_Residual_Sig（somatic）| 11,106 | 39.2% |
| GermlineAsmDbeta_Sig（germline）| 7,733 | 27.3% |
| SomaticResidualDbeta_Sig（somatic）| 3,031 | 10.7% |

---

## §4 結果：決策性拆分（互斥四分桶）

| 桶 | 意義 | 數量 | 佔 Significant |
|---|---|---:|---:|
| **both** | 有純 tumor 軸 + 有對照軸（證據最厚）| 9,103 | 32.1% |
| **tumor-only excl** | 只有純 tumor 軸，無對照 | 4,502 | 15.9% |
| **contrast excl** | 只有對照軸，無純 tumor 軸 | 6,619 | 23.4% |
| **neither（only all-pool）** | 純靠混池軸，無 tumor 也無對照（最弱）| 8,103 | 28.6% |
| 合計 | | 28,327 | 100% |

**兩個 headline cut（含重疊）：**
- **可用 tumor-only 偵測到結構（≥1 純 tumor 軸）= 13,605 = 48.0%**
- **有 tumor-vs-normal 對照軸（需 normal 才驗證）= 15,722 = 55.5%**
- **只靠混池軸、最不可信 = 8,103 = 28.6%**

---

## §5 結果：per-VC 拆解

| VC | n | tumor-only% | contrast% | (somatic% / germline-only%) | only-pool% |
|---|---:|---:|---:|---|---:|
| **Strong** | 22,310 | 48.7% | 56.2% | (40.7% / 15.4%) | 28.6% |
| **MultiGroupNoLabel** | 1,899 | **100.0%** | 57.3% | (50.9% / 6.4%) | **0.0%** |
| **LabelShift** | 1,809 | 22.5% | **81.6%** | (54.1% / 27.5%) | 15.6% |
| **PermanovaLocation** | 1,756 | 14.4% | 27.3% | (27.3% / 0.0%) | **65.1%** |
| **LOH-Structure** | 376 | 33.5% | 33.0% | (6.9% / 26.1%) | 44.9% |
| **StructureNoLabel** | 177 | 29.9% | 14.7% | (14.7% / 0.0%) | **68.9%** |

**模式：**
- **MultiGroupNoLabel = 最純 tumor 軸**（100% tumor-only、0% only-pool；定義上就是 within-HP tumor 子結構）。
- **PermanovaLocation / StructureNoLabel = 最弱**（65–69% only-pool；純混池 PERMANOVA 觸發，無 tumor 也無對照佐證）。
- **LabelShift = 最 contrast-heavy**（81.6%，Δβ 驅動，需 normal 區分 somatic/germline）。
- **Strong（骨幹）= 混合**：48.7% 有 tumor、56.2% 有對照、28.6% 純混池。

---

## §6 解讀與含意

1. **判定的「有結構」≠「tumor 特異結構」**。骨幹 Strong（佔 78.8% Significant）來自混池 clustering，其「有結構」可能是：(a) 真 tumor 子結構、(b) germline HP phasing（兩邊都有）、(c) 就只是把 tumor 和 normal reads 分開（sampleOrphan / tumor-vs-normal confound）。判定本身分不出。
2. **能自證 somatic 的只有 C 類對照軸**（HP_Residual / SomaticResidual = 扣掉 normal 後仍在）。這些軸觸發於 55.5% 的 Significant，但其中真 somatic（非 germline ASM）更少（見 §5 somatic% 欄）。
3. **28.6% 純混池、無任何 tumor 或對照佐證** → 這 8,103 個是「最該被降級審視」的集合。
4. 對 subclonal-reconstruction 論文目標：**應把 read-set provenance 做成判定的一個維度**，把「tumor 可自證 + 對照確認」與「純混池」分層，而非全算成 Significant。

---

## §7 Caveats（誠實邊界）

- **單樣本 HCC1395**，±5000bp，normal median 58×。低 normal 覆蓋處 baseline 扣得不乾淨。需 COLO829 等跨樣本復現。
- **軸門檻（AUC≥0.7、|Δβ|≥0.1）是 first-pass heuristic**，非驗證過的最佳邊界；換門檻會移動比例。
- 「has contrast axis」**不等於 somatic**：C 類含 germline-ASM（GermlineAsmDbeta = germline 訊號）。真 somatic 拆分見 §5 somatic% 欄 + ask (2)。
- 「somatic-specific」仍含 **cis-ASM**（等位特異但非次克隆）；分 cis-ASM vs subclone 需 normal-anchored cis-control（既有 memory 警告）。
- within-HP 軸**無統計 null**（heuristic 幾何切）— 是目前唯一未檢定的軸（另案待補 permutation null）。

---

## §8 下一步（接 ask 2）

C 類對照軸**本身已帶 germline vs somatic 區分**，可現在就做 ask (2) 的 germline/somatic 拆分；**CN 歸因需先做 coverage-peak 自估**（SEQC2 不可當輸入）。建議序：先 germline/somatic（用既有欄位）→ 確認後補 CN coverage-peak → 三方歸因完整。

---

**Provenance**：所有數字注入自 `docs/methodology/_assets/20260616_readset_provenance/readset_stats.json`（本 session 由 `output/_wg_d1_unified/significance_summary.csv` 30,490 列計算並讀回）。程式碼 trace 對應 commit 990175d（feat/summary-nreadsvalid）。
<!-- provenance-verified: 所有 metric 來自同層 _assets/readset_stats.json，本 session 計算讀回；程式行號為 trace 結果 -->
