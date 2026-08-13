---
title: LongLineage 充分性稽核 — 能否承擔「純遺傳共現 → 拓撲 → read 標籤」主角色
date: 2026-08-06
task_type: C production（前置可行性判定）
scope: LongLineage 全 repo + 與 InterSubMod 交界
verdict: 前提錯誤已澄清；三條路線待裁決
data_sources: scratchpad/ll_audit.json, scratchpad/methyl_gate_chain_verified.md, LongLineage/docs/reports/20260720_HCC1395完整科學運算與parity報告_01.json
build_branch: research/subclonal-reconstruction-202606
---

# LongLineage 充分性稽核與路線裁決

## 0. 必須先講的前提錯誤

你的裁決基於「LongLineage 完整可用，只要確認並串接」。**這個前提與實測不符。**

`LongLineage/docs/reports/20260720_HCC1395完整科學運算與parity報告_01.json`
（HCC1395 全 autosome 唯一一次完整跑，已凍結）：

```
79,687 sites
   → m1_stable          12,851
   → m2_eligible             5      ← 甲基 gate
   → joint_signature_pass    0
   → topology_units          0      ← 拓撲產出為零
```
134,278 個 pair row 中 **134,276 被標 `INELIGIBLE_M2_SCREEN`**。

**LongLineage 的 formal topology 端點目前輸出 0 個 unit。**
這不是「先用著再改」，是「現在就是空的」。repo 自己的 `ADR-0007` 明文承認：
> "The current LongLineage formal topology endpoint reports zero units after M2,
> joint-signature and directionality gates."

7 個稽核問題的 verdict：**0 個 SUFFICIENT**，6 個 `INSUFFICIENT_ARCHITECTURE_OK`，
1 個 `INSUFFICIENT_ARCHITECTURE_BLOCKS`（寫 BAM）。

---

## 1. 好消息：你的架構判斷是對的，而且 repo 自己也這樣規定

你要的分工「LongLineage 出 artifact、InterSubMod 做甲基分析、Python 出 HTML」
**是 LongLineage 唯一允許的分工** —— 這由 CI 機械強制：

- `scripts/ci/check_source_boundaries.sh:17`：Python 只能存在於 `presentation/`
- presentation 的 Python **禁止 import pysam/cyvcf2、禁止呼叫 samtools/bcftools**

⇒ 想把 ISM 的分析腳本搬進 LongLineage 會被 BOUNDARY FAIL 擋下。你的分工是唯一可行解。

### 而且最難的部分完全乾淨

`src/solver/` **3,021 行**（obligation_bnb / small_q_oracle / terminal_subset_dp /
topology_record / topology_router / vertex_set_ranker）
grep `methyl|m1_|cpg|label` **零命中**。solver 輸入 `ReadPatternObservation` 只有
`pattern_raox + base_qualities + multiplicity`，且由 **整個 block 的所有 read** 投影而成（非 M1 core）。

⇒ **拓撲求解核心不用動**。

---

## 2. 甲基 gate 的完整鏈條（逐行實證）

| 步 | 位置 | 內容 |
|---|---|---|
| ① | `dataset_producer.cpp:519` | `row.label = assignment.coarse_label`（M1 甲基群） |
| ② | `site_science.cpp:150` | `m1_flagged = analysis->stable_null_multigroup`（coarse_groups≥2 && !unstable && ARI_min≥0.8） |
| ③ | `m2/eligibility.cpp:593-596` | `if (!input.m1_flagged) return {kNotRun, kM1NotFlagged, false, false}` ← 硬短路 |
| ④ | `eligibility.cpp:600-602` | `group_count < 2` 直接 **throw** |
| ⑤ | `joint_signature.cpp:510/516/521` | `categorical_association(labels, signatures)` ← **甲基群是自變數** |
| ⑥ | `topology_membership.cpp:386-397` | `if (!candidate.m2_eligible) continue;` ← 不建 unit |
| ⑦ | `topology_membership.cpp:512-518` | read-pattern adapter 入口**再擋一次**（fail-closed） |

**M2 不是獨立的遺傳篩選器，而是「M1 甲基群的品質關卡」**。

### 純遺傳訊號其實已經算好了，只是被擋住

`site_cooccurrence.cpp:323-341` 的 **endpoint B**（4×4 `JointAlleleCounts` +
`summarize_four_state_relations` 出 focal-ancestor / partner-ancestor / branching 相容性）
**全程零甲基**，只用 `focal.covering_read_indices`（純位置條件）× 兩端 allele call。

但 `site_cooccurrence.cpp:251-258` 在函式開頭就對 `!stable_null_multigroup` **early-return 空 pairs**。

⇒ **「純遺傳訊號被甲基 gate 誤殺」的單一最關鍵位置**。

### 第二層耦合：read universe 本身

`block_reader.cpp:476-487` 硬性丟棄沒有 MM/ML tag 的 read，
`contracts/v1/science_parameters.json` 凍結 `require_mm=true / require_ml=true`，
綁 governance gate `IO_MM_ML`。

⇒ 就算拆掉 M1/M2 gate，production 路徑的 read 母體仍是「有甲基標記的 read」。

### 沒有任何 flag 可繞過

`DatasetProductionOptions` / `SiteScienceOptions` / `ProductionManifest` / CLI
全部查過，無旁路選項。

---

## 3. 🔴 寫 BAM：架構層禁止，不是可調參數

這是唯一 `INSUFFICIENT_ARCHITECTURE_BLOCKS` 的項目。

| 障礙 | 證據 |
|---|---|
| JSON Schema **const false** | `schema/core/production_input_authority.schema.json:69`（`tagged_bam_output_allowed`）、`:108`（`persisted_tagged_bam`） |
| 4 處 runtime 硬拒 | `apps/longlineage_main.cpp:188`；`src/compat/regional_io.cpp:391,465`；`src/validation/regional_compat_validator.cpp:1121` |
| FILE_CENSUS 拒絕多餘檔案 | `artifact_validator.cpp:3091-3101`；16-artifact catalog 無 BAM 格式 |
| 禁止 in-place 加 tag | `artifact_validator.cpp:1255,1614` 要求 `input_snapshot_before_sha256 == after` |
| **自我不可重現** | aux tag 屬 `typed_aux` 摘要契約（`SAM_CORE_AND_ALL_TYPED_AUX_EXCEPT_RG_EXACT_V1`, FAIL_CLOSED）→ **LongLineage 讀回自己寫的 tagged BAM 會改變 read identity 與 semantic SHA** |
| 讀取路徑不適合 | `block_reader.cpp` 用 `sam_itr_queryi` + 5000bp halo → 同一 read 跨 block 重複解碼、遠離 focal 的 read 完全未讀 → **結構上無法產出完整且無重複的 BAM** |

**可行的替代**：tagged BAM 只能是 **post-freeze、run-root 外的獨立 export 工具**
（讀 frozen artifact + raw BAM，重算 `SHA-256(QNAME)` 做 join）。
這不違反上述任何一條，因為它不是 LongLineage run 的產物。

---

## 4. read 標籤輸出：交界目前是空的

| 路徑 | per-read 資訊 | 有無 lineage 標籤 |
|---|---|---|
| production `site_reads` | ✅ read_id + HP + PS + allele_call(R/A/O/X)，**且不受 M1 gate 影響**（每個 focal site 無條件寫） | ❌ 無 clone/lineage 欄位 |
| production `m1_assignments` | ✅ read_ids[] + labels[] | ⚠ 是**甲基**分群標籤 |
| solver 內部 | ✅ `canonical_read_ids`（`topology_membership.cpp:698-701`） | ❌ **從未序列化，只活在記憶體** |
| compat 路徑 | ❌ 只有 pattern→count 聚合 | ❌ |

⇒ 「LongLineage 輸出分群標籤 → InterSubMod 吃」這個交界**需要新寫 writer**。
好消息：`site_reads` 的 `read_id` 不受 gate 影響，可直接當 join key。

⚠ **語意警告**：read → vertex 在三種情況下**無唯一解**（含 X 的 partial-subcube read、
`best_vertex_set_tie_class` 長度 >1、family incomplete/abstain）。
這是模型固有性質，artifact 必須用 status enum 表達，**不可強制唯一標籤**，否則就是 overclaim。

---

## 5. 三條路線

### 路線 A — 拆 production 路徑的 M1 gate

| | |
|---|---|
| 工作量 | **large** |
| 必做 | ① 刪 `site_cooccurrence.cpp:251-258` early-return ② **定義純遺傳的 partner 選擇規則**（唯一必須寫新演算法處）③ 改 `topology_membership.cpp:386` 謂詞 + funnel 守恆 + canonical digest ④ 解除 MM/ML filter（換 contract 版本）⑤ **同步獨立 validator**（`artifact_validator.cpp` 5,367 行，CMake 禁止 link producer core → 改動成本 **2×**）⑥ 更新 parity 測試 |
| 風險 | pair row 從 134,278 暴增到全 partner universe（79,687 × ±5000bp），記憶體需重估；contract SHA 換版 → 既有 VALIDATED_FROZEN run 不再可比；與 P3 M1 legacy parity 目標直接衝突（需判定 P3 是否退場） |

### 路線 B — 扶正已凍結的 compat 純遺傳路徑

`src/compat/` 的 `PYTHON_V2_DESCRIPTIVE_REGIONAL`：

- 8 個 compat 檔 grep `methyl|mm_ml|ml_raw|cpg` **全 0**，read filter 只有 FLAG + MAPQ≥20，**不要求 MM/ML**
- **全 chr1-22 已 VALIDATED_FROZEN：8,222 regions / 20,119 units / 106,559 patterns**
- 已有獨立 binary `longlineage-regional-compat` + 獨立 validator

| | |
|---|---|
| 工作量 | **medium** |
| 阻礙 | 🔴 `ADR-0007 §1` **明文禁止**該端點宣稱 clone / ancestor / temporal order / formal topology discovery → 扶正是**研究誠信層決策**，須你確認，工程不可單方面改 |
| 🔴 科學缺陷 | 它複製的是 **ISM 舊版幾何 50kb 分群**，而 ISM 自己**已用 read-linked 方法取代**。稽核原話：「只能當工程模板不能當科學規格」 |

### 路線 C — 新建 read-linked component builder（稽核推薦的科學正解）

照 ISM **現行**的 `build_strict_ps_hp_regions.py` 合約在 LongLineage 實作：

- container = chrom × exact non-missing PS × HP family
- node = 有 fixed R/A 的 sSNV
- edge = ≥T 條 distinct `read_id` 兩端皆 fixed R/A
- 取無向連通分量；距離只報告不建邊

| | |
|---|---|
| 工作量 | **large**（新模組 `src/linkage/`）但**輸入沿用既有 `site_reads` artifact**（不受 gate 影響） |
| 加分 | `exact_ps_k12_partition.cpp` 已是 C++17 且有 Python parity 對照，移植風險低，可接上 `topology_router` 目前直接 abstain 的路徑 |
| 優點 | 這才是 ISM 的**現行**科學規格，非舊版幾何分群；避開 M1/M2 整條甲基鏈 |

---

## 6. 🔴 一個必須先回答的科學問題

稽核在 Q3 提出警告：

> **ISM 已把「跨 unit 整樹 ranking」判為 CLOSED NEGATIVE（等機率、不可枚舉）。**
> 若目標是單一整樹，這一項可能在科學上就不成立，應先做可行性判定再寫程式。

所以「整體拓撲」要先定義清楚是：
- (a) **forest + 跨 unit 一致性標記**（相容/衝突），或
- (b) **單一整樹**（← 可能撞既有 NEGATIVE 結論）

---

## 7. 其他實測限制

| 項目 | 狀態 |
|---|---|
| `longlineage run` | 未實作 stub（`longlineage_main.cpp:967-970`），7 樣本目標現況無法啟動 |
| `dataset-gate` | 硬編碼 `dataset_id=='HCC1395'`（`:235-238, :258-259`） |
| 7 樣本總時間 | **NOT_FOUND** — `phase_ledger.json:272-276` 明列 `SEVEN_DATASET_W24_W40_RUNS_NOT_EXECUTED`，且稽核文件**明文禁止單樣本外推** |
| 記憶體上界 | 最大樣本 HCC1937（**472 GB** BAM，HCC1395 的 1.62 倍）**完全無實測證據**；`PHYSICAL_GLOBAL_MEMORY_BOUND_NOT_PROVEN` |
| 1.94 TB raw BAM 的 SHA-256 preflight 秒數 | **NOT_FOUND**，但很可能是 cohort 主要 wall time |
| phase 狀態 | P3/P4/P5/P7/P8 全 BLOCKED，`release_attestation = NOT_READY` |

---

## 8. 待裁決

1. **路線 A / B / C 三選一**（阻塞所有後續工作）
2. **「整體拓撲」= forest+一致性 還是 單一整樹？**（後者可能撞 CLOSED NEGATIVE）
3. **tagged BAM**：接受「post-freeze 外部 export 工具」的形式嗎？（run 內產出被架構禁止）
