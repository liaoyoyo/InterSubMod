---
title: KEEP / REMOVE 判定 + 保留分類 + 低信心原因 — 統一綜整（帶絕對數字 + 方法邏輯）
date: 2026-06-16
type: methodology / consolidation
scope: single-sample HCC1395, paired, ±5000bp window, 全基因組 30,490 TP SNV
status: in_progress / partial（單樣本；分類用既有 a-priori 欄，verdict 為 D1 現狀；tumor-only 軸待移除）
tier: L1（絕對數字第一手讀回 consolidation.json）/ L2（生物學詮釋）
build_branch: feat/summary-nreadsvalid @ 990175d
data_sources: docs/methodology/_assets/20260616_readset_provenance/consolidation.json, docs/methodology/_assets/20260616_readset_provenance/residual_axes.json, docs/methodology/_assets/20260616_readset_provenance/subclone_morphology.json, docs/methodology/_assets/20260616_readset_provenance/venn_structural.json
related:
  - docs/methodology/20260616_readset_provenance_observation_01.md
  - docs/methodology/20260616_tumor_only_structure_axis_cpp_change_audit_01.md
---

> ⚠️ **PARTIAL / 單樣本** — HCC1395，±5000bp，30,490 TP SNV。分類用既有 a-priori 軸（合法、會判別）；KEEP/REMOVE 為 D1 現狀 verdict；tumor-only unsupervised 軸（無效）**不納入**本分類、待 C++ 移除。

# KEEP / REMOVE 判定 + 保留分類 + 低信心原因 — 統一綜整

## §0 TL;DR

- **總位點 30,490**。KEEP **28,327 (92.9%)** / REMOVE **2,163 (7.1%)**。
- 保留的 28,327 裡：**有 a-priori 標記（高信心）19,440 (68.6% of keep)** / **無標記（低信心 epiallele）8,887 (31.4% of keep)**。
- 判別性已證：a-priori 軸在 Noise **0%**、結構類 10-42% → 真判別；unsupervised 軸 Strong≈Noise 83% → 不判別、已棄。

---

## §1 KEEP vs REMOVE — 標準 + 數字 + 方法

**標準（D1 統一判定）**：`KEEP = Significant`（VerificationClass ∈ 6 valid 類）；`REMOVE = DispersionStructure（離散度假象）+ Noise（無結構）`。

| 判定 | 數量 | % 全體 | 經過的方法 / 邏輯 |
|---|---:|---:|---|
| **KEEP**（Significant）| **28,327** | **92.9%** | Stage④ override：任一結構軸觸發（cluster/Δβ/PERMANOVA-location/within-HP）|
| REMOVE — Dispersion | 1,700 | 5.6% | PERMANOVA 顯著但 **PERMDISP（analytic-F）判定純離散度差異**（群心無位移）|
| REMOVE — Noise | 463 | 1.5% | distance + level + label + PERMANOVA **全部無訊號** |
| 合計 | 30,490 | 100% | 加總檢核 ✓ |

---

## §2 保留的分類（5 類，依 a-priori 方法優先序）+ 強度 + 用途

> 互斥優先序：subclone > somatic-residual > germline > CN > structure-no-marker。每類標「用哪個方法、null 是什麼、合不合法、對題目用途」。

| 類別 | 數量 | %全體 | %保留 | 方法（標籤如何定義）| null / 合法性 | 對題目用途 |
|---|---:|---:|---:|---|---|---|
| **KEEP_Subclone** | **5,473** | 18.0% | **19.3%** | germline-HP 內 germline vs **1-1 somatic-carrier**（haplotag）/ ALT vs REF | permutation（a-priori 標籤，**無 double-dip**）| ✅ **題目核心**：clone 有獨特甲基 |
| **KEEP_SomaticResidual** | **9,074** | 29.8% | **32.0%** | tumor HP **扣 normal baseline** 後殘差顯著 | permutation（HP a-priori）| ✅ somatic 但未解析到 clone |
| **KEEP_GermlineASM** | 3,407 | 11.2% | 12.0% | HP1 vs HP2 甲基差（normal 也有）| permutation（HP a-priori）| ❌ germline confound（題目要排除）|
| **KEEP_CN_LOH** | 1,486 | 4.9% | 5.2% | HP1/HP2 read 比失衡（**read 高峰/HP-tag 推估**）| 比例門檻 | ⚠ CN 軸（非甲基）|
| **KEEP_StructureNoMarker** | **8,887** | 29.1% | **31.4%** | Significant 但**無任何 a-priori 標記**（靠 unsupervised cluster / pool-PERMANOVA / within-HP）| ❌ **無合法 null**（epiallele 混淆）| ❌ **低信心、不可驗證** |

**強度（StrengthScore，舊公式）per class**：Subclone 0.345 / SomaticResidual 0.346 / GermlineASM **0.444** / CN 0.354 / StructureNoMarker 0.350。
→ 舊強度**分不出** subclone vs germline（都 ~0.34-0.44，germline 反而最高）→ 這正是改 somatic-weighted + 元件輸出的理由（但敏感度已證權重排序 ρ=0.998 幾乎無差）。

**高信心保留合計**（subclone+somatic+germline+CN，有 a-priori 標記）= **19,440 (63.8% 全體 / 68.6% of keep)**。

---

## §3 低信心位點 — 原因 + 數字

| 低信心類別 | 數量 | %全體 | 原因（經過哪些方法判定）|
|---|---:|---:|---|
| **KEEP_StructureNoMarker** | 8,887 | 29.1% | 被 verdict 保留，但**只靠 unsupervised 證據**（pool-clustering / within-HP / pool-PERMANOVA）；這些軸已證**對 Noise 也觸發（Strong≈Noise 83%）= epiallele 混淆、無合法 null** → 不可驗證為 subclone |
| REMOVE_Dispersion | 1,700 | 5.6% | PERMANOVA 顯著但 PERMDISP 判定**只是群內離散度不同、群心無位移**（Anderson 2006 假象）|
| REMOVE_Noise | 463 | 1.5% | 所有軸（距離/level/label/PERMANOVA）**皆無訊號** |
| **低信心合計** | **11,050** | **36.2%** | — |

> 核心邏輯：低信心 ≠ 全部丟。`StructureNoMarker`（8,887）是「有結構跡象但無 a-priori 標記可驗證」——對 reconstruction **無用**（不能 anchor 到 clone/ordering），但保留為「待 cross-sample / cis-control 釐清」。`Dispersion + Noise`（2,163）是「真無可處理訊號」直接 REMOVE。

---

## §4 方法邏輯總表 — 每個判定經過哪些方法

| 軸（方法）| raw 觸發 | %全體 | 標籤來源 | null | 合法？|
|---|---:|---:|---|---|:---:|
| SubcloneDbeta（同單倍型）| 4,596 | 15.1% | haplotag 1-1/2-1（somatic sSNV）| permutation | ✅ |
| AltSubcloneDbeta（allele）| 4,454 | 14.6% | ALT/REF allele | permutation | ✅ |
| SomaticResidual（扣 normal）| 3,031 | 9.9% | germline HP + normal baseline | permutation | ✅ |
| HP_Residual（tumorHP−normalHP）| 11,632 | 38.2% | germline HP | permutation | ✅ |
| GermlineASM | 7,733 | 25.4% | germline HP（normal）| permutation | ✅ |
| Potential_LOH（CN）| 2,605 | 8.5% | HP-tag 比 / coverage | 門檻 | ✅（描述）|
| ~~tumor-only unsupervised~~ | ~~97.6%~~ | — | **clustering 自造** | ~~PERMANOVA~~ | ❌ **double-dip，已棄** |

**判別性鐵證**：a-priori 軸在 Noise_Uniform/Uncorrelated **0%**；unsupervised 軸 Strong 83% ≈ Noise 83%。→ 保留判定**只信 a-priori 軸**。

---

## §4.5 結構偵測軸交集 Venn（距離/level/label/PERMANOVA）

> data_source: venn_structural.json。marginal：距離cluster 23,007 (75.5%) · level 7,910 (25.9%) · label 29,615 (97.1%) · PERMANOVA 29,726 (97.5%)。

| 層級 | 主要組合 | 數量 | %全體 |
|---|---|---:|---:|
| **4 軸全交** | 距離∩level∩label∩PERMANOVA | 6,090 | 20.0% |
| **3 軸**（18,249 / 59.9%）| 距離∩label∩PERMANOVA | 16,565 | 54.3% |
| | level∩label∩PERMANOVA | 1,619 | 5.3% |
| **2 軸**（5,327 / 17.5%）| label∩PERMANOVA | 4,992 | 16.4% |
| **1 軸**（497 / 1.6%）| label 223 / PERMANOVA 193 / level 59 / 距離 22 | 497 | 1.6% |
| **0 軸** | 無任何結構軸 | 327 | 1.1% |

**讀法**：label+PERMANOVA 近 97% 普遍觸發（germline-HP/allele 關聯無所不在）→ 「有結構」≈ 普遍；判別力來自 conditioned 軸非這四原始軸。

### 四軸判定與意義
| 軸 | 怎麼算 | 意義 | 合法性 |
|---|---|---|---|
| 距離 cluster | read×read BERNOULLI 距離 → UPGMA → silhouette | reads 是否依甲基 pattern 自然分群 | ⚠ unsupervised double-dip（Noise 也 83%）|
| level | HP 內依 mean-β 高低切群 | 整體甲基高/低子群（距離漏的 level 位移）| ⚠ heuristic 無 null |
| label 關聯 | 已知標籤(HP/allele) chi-square + Δβ | 已知單倍型/等位與甲基相關 | ✅ a-priori，多為 germline |
| PERMANOVA | labeled 群全-CpG centroid 是否不同（permutation）| labeled 群甲基整體分布分離 | ✅ a-priori（+PERMDISP）|

## §4.6 低信心位點為何低（StructureNoMarker 8,887）
Significant 但**無任何 Δβ 標記**（subclone/somatic/germline/CN 全不顯著）；結構只來自 (a) 距離cluster（unsupervised，Noise 也觸發=epiallele）、(b) level（無 null）、(c) PERMANOVA pattern 差但無 effect-size 錨。→ **「有跡象但不可靠/不可錨定驗證」，非真無訊號** → 對 reconstruction 無用。

## §4.7 Dispersion 判定（REMOVE_Dispersion 1,700）
PERMANOVA 測群心**位置**差異，但 Anderson 2006：群心相同、只是**離散度(spread)不同**時 PERMANOVA 也顯著（假象）。→ 對顯著群跑 **PERMDISP（analytic-F p，對 scipy 1.4e-14）**：PERMANOVA 顯著且 PERMDISP 不警告=真位移（保留）；PERMANOVA 顯著但 PERMDISP 警告=離散度假象 → **REMOVE**。剔除「一群散一群密」的偽結構。

## §5 現狀 vs 計劃 + caveat

- **現狀**：D1 verdict（KEEP 28,327）已落地；本分類用既有 a-priori 欄計算（未進 verdict）；tumor-only 軸在未提交 C++ 中（待移除）。
- **計劃**：(1) 移除無效 tumor-only unsupervised 軸；(2) verdict 細分 KEEP 為上 5 類（Strong 不再鐵板一塊）；(3) StructureNoMarker 降為低信心旗標。
- **caveat**：單樣本 HCC1395；Subclone-specific 仍含 cis-ASM（需 normal-anchored cis-control 分 cis vs subclone，故 5,473 為**上界**）；CN 軸 coverage-peak 自估尚未實作（現用 HP-tag 比）。

---

**Provenance**：絕對數字注入自 `consolidation.json`（由 `output/_wg_d1_unified/significance_summary.csv` 30,490 列計算讀回，本 session）。判別性數字見 readset/residual/tumor_cluster_pilot JSON。
<!-- provenance-verified: 所有數字來自同層 _assets/consolidation.json + residual_axes.json + subclone_morphology.json，本 session 計算讀回 -->
