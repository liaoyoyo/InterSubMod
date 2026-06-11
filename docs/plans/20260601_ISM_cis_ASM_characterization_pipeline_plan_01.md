---
title: ISM cis-ASM characterization pipeline — 改進計劃書（gap 分析 + 路線 + 開發 workflow 架構）
date: 2026-06-01
status: plan / 待用戶確認（未執行）
scope: 單樣本 HCC1395 paired_full；characterization（非 filter）
source: 7-dimension gap-analysis workflow (wf_a5b4649f, 7 agents) + 本 session BRCA2 cis-test 發現
commit: 859c55a
guardrail: cis-ASM characterization；不重開 ASM-as-discriminator / TO-germline-FP filter（concluded dead）
---

# ISM cis-ASM characterization pipeline — 改進計劃書

> **狀態：計劃，未執行。** 任何 code 撰寫 / ISM 改動等用戶確認後才動（Python-first；C++ 走 `/cpp-change` Hard Gate）。

## §0 TL;DR + 定位

**目標**：把現在的 ad-hoc ISM 分析，重新定義成一個能**系統性「找出 + 嚴格分析 + 分層輸出 + 持續驗證」cis-candidate ASM 位點**（像 BRCA2：somatic 子克隆有真實、germline-haplotype-controlled 的甲基差異）的 pipeline。

**定位（與已關閉方向的明確區別，回應 Research Direction Guard）**：
- ✅ 這是 **cis-ASM characterization**（描述「哪些 somatic 位點有真實甲基差異 + 證據強度分層」）。
- ❌ **不是** concluded-dead 的「ASM 當 TP/FP discriminator / filter」（pure-methylation 判別、TO germline-FP G1-G7）。每個 scan record schema **嚴禁含任何 TP/FP classifier score**，不餵 F1 filter。
- 全域 NEGATIVE（甲基不能系統性判別）**仍成立**；本系統找的是「稀少但真實的個別 cis 案例」。

## §1 診斷 — 為何現在的 ISM 無法乾淨分析 BRCA2 這種位點（3 個根因，全 L1 源碼/實測佐證）

| # | 根因 | 證據 | 後果 |
|---|------|------|------|
| **R1 抽取層 type-erased + ephemeral** | MSA C++ `MethylationSiteDetail`（`include/msa/Types.h:82-95`）**無 mod-type 欄**；`MethylHaploExtractor.cpp:201-218` 收 'm'+'h' 但丟棄 mod_code → 5mC/5hmC 變兩列未標記。`19_*.sh` 跑完 `rm -f $L1` 刪 Level-1。Level-3 ASM 統計比 VCF 來源非 HP group。 | **每個下游腳本被迫 max-collapse**（−0.054 砍半 artifact 的根源）；**BRCA2 是唯一能重分析的位點**（只因手動 clone 保住 Level-1）；MSA 內建統計不可用。 |
| **R2 cis-test 是 demo 不是 detector** | `script 18` genome-wide pivot **只讀 tumor**（line 78），無 normal anchor；normal-anchored cis-test（`script 34`）只跑 10 個手挑位點。 | genome-wide 層級**無法區分 somatic-association vs germline-drift**；association→causation 無分層；small-n 離散下限（36% 位點 n_cpg≤20 測不到）未處理。 |
| **R3 LOH/CNV 是標籤不是模型** | LOH 只是 point-in-region 註解（whole-region NGS）；**CNV 完全沒整合**（整數 CN 在 `gain_cn.bed` 但沒人接）；**43% HP-axis 跑在 LOH 內**（10,220/23,840，違反自身 docstring）；β 在多拷貝區是跨拷貝混合估計。 | LOH 區 HP-axis 噪音掩蓋真 nonLOH cis-candidate；CN-gain（如 BRCA2 CN=5）放大未排除。 |
| **R0 可重現性斷裂** | 34 個一次性腳本、無依賴拓樸、34/34 hardcode 絕對路徑、`glob[0]` 抓「碰巧第一個」、Level-1 被刪、無 regression test。 | 重現 BRCA2 d_cis=−0.1418 須人腦記憶腳本順序；換樣本即斷。 |

**一句話**：現在的系統能「測」BRCA2，但**測不乾淨、存不住、找不到更多、也驗不了改動**——因為抽取層丟資訊、cis-test 沒 genome-wide 化、LOH/CNV 沒進模型、整套是探索期 throwaway 腳本。

## §2 目標 — 重新定義的 ISM cis-ASM characterization 系統

```
[版本化抽取]  BAM-direct extractor (pysam, mod_code 保留 5mC/5hmC) + tumor+normal
   │  → Level-1-plus（持久化 cache，鍵入 BAM/REF/params hash）
   ▼
[cis-ASM core]  per-locus 一級輸出：
   │  · HP-axis(HP1 vs HP1-1, nonLOH-only) + ALLELE-axis(LOH 合法)  ← axis-gate
   │  · normal-anchored cis 三點(A=normalHP1/B=tumorHP1/C=tumorHP1-1) → d_cis/d_drift
   │  · per-tag cohesion(silhouette) + k-selection(multi-subclone) + cross-tag membership
   │  · CN-aware null(整數 CN 加權) + power-class(離散下限)
   │  · causation 層: mechanical-cis(CpG 破壞) + location profile + subclone-partition
   ▼
[分層 verdict]  每位點: power_class · cis_tier(T0 no-ASM→T3 cis-candidate) · causation_tier · LOH/CN-class
   ▼
[系統性 scan]  two-stage: 便宜 dual-axis shortlist(~32) → 昂貴 cis-confirm → ranking + 證據卡
   ▼
[持續驗證]  BRCA2 regression test + git-diff hook + provenance manifest
```

## §3 分階段路線（Phase 0-4；P0 先；幾乎全 Python，1 個 conditional C++ 延後）

| Phase | 內容（對應維度 improvement）| cost | C++? | 驗證（exit criteria）|
|-------|------|------|------|------|
| **P0 可重現地基** | 建 `pipeline/`（loci.yaml SoT + stages + Makefile）；**BRCA2 canonical regression test**（凍結 d_cis=−0.1418/d_drift=−0.0217/p_cis≈0）；舊新並存逐位對齊 | high* | ❌ | `pytest test_brca2_regression` 對未改 code PASS；把 max-collapse 改回 buggy 必 FAIL |
| **P1 抽取 substrate** | `35_bam_extract_modtype.py`（mod_code 欄，5mC/5hmC 分離）；Level-1 持久化 cache（停止 rm，寫 `/` 702G）| med | ❌ | BRCA2 5mC-only Δβ 重現 **−0.122 非 −0.054**；5hmC ~0.03；cache HIT 不重抽 |
| **P2 cis-core + LOH/CNV** | `cis_asm_core.py`（genome-wide cis 欄 + cis-ladder T0-T3 + power-class）；`cn_annotate.py`（整數 CN）；**axis-gate**（HP nonLOH-only）| med | ❌ | BRCA2 重現 cis-test；gate 後 HP-axis 100% nonLOH（0 個 LOH 違規）；BRCA2 CN=5 |
| **P3 scan + 子克隆 + causation** | `scan_cis_candidates.py`（two-stage + testability funnel + ranking + 證據卡）；k-selection + cross-tag membership；mechanical-cis CpG-disruption + location profile + subclone-partition | med-high | ❌ | scan grep 到 BRCA2 cis-cand；het-null 當負控落 T0/T1；輸出 cis-candidate 排名清單 |
| **P4 持續驗證閉環** | git-diff hook（`asm_pipeline_regression.sh`，**exit 0 advisory**）；provenance manifest；`/pipeline-manifest` 對接 | low | ❌ | 編輯 stage 觸發 advisory；manifest DAG 無 orphan |
| **(future) conditional C++** | 若 Python cn_annotate/axis-gate 連 ≥2 cycle 穩定 + 需跨樣本上線 → 評估內建進 MSA/inter_sub_mod | high | ✅ | **本輪不動**；屆時走 `/cpp-change` PDD + 編譯 Hard Gate + ctest |

\* P0 cost high 是因為要把 34 腳本逐位對齊重構；但**風險低**（golden test 守門，不改數值）。

**關鍵發現：幾乎全部改進是 Python（needs_cpp=false）**——印證 Python-first。唯一 C++ 項目（CN 內建 MSA）明確標 future-only，不在本輪。

## §4 開發 workflow 架構（git-diff 持續驗證）

開發階段（**確認後**才跑）用第二個 workflow，逐 Phase 推進，每個 stage 改動都對 BRCA2 regression test 驗證：

1. **每 Phase = 一個 workflow 段**：fan-out 寫各 stage 的 Python module → 各自對 golden BRCA2 數值 self-verify → 主 agent 合併 + 跑 `pytest`。
2. **git-diff 研究分析**：每個 stage commit 前，diff 舊輸出 vs 新輸出（BRCA2 三數逐位對齊）；數值若變必 commit message 標 `[golden-update]` + 理由。
3. **regression hook**：`asm_pipeline_regression.sh`（PostToolUse Edit|Write on `pipeline/stages/*.py`）advisory 提醒跑測試（**exit 0**，不 mask）。
4. **Level-1 cache**：避免每次重跑 MSA（BRCA2 從 multi-hour → 秒級 re-read）。

## §5 Guardrails（必遵守）

- **characterization only**：每個 scan record / verdict 嚴禁含 TP/FP discriminator score；不餵 F1 filter。與 concluded-dead「ASM-as-filter / TO-germline-FP G1-G7」明確區隔。
- **Python-first**：本輪 0 個 C++ 改動；唯一 C++ 項 future-only + `/cpp-change` PDD + 編譯 Hard Gate。
- **association ≠ causation**：cis_tier T0-T2 嚴禁用 causation 字眼；只有 T3 用「候選-需-causation-證據」且必帶 missing_evidence 清單。
- **單樣本**：HCC1395；cross-sample 是 future（untestable 位點標 `untestable_in_HCC1395_single_sample`，不丟不標 NEGATIVE）。
- **磁碟紀律**：cache/persist 寫 `/`（702G）非 big7/big8（97/98% 滿）；`TMPDIR=/big7_disk`（/tmp 800GB 前科）。
- **不捏造**：golden 值來自實測 `control_cohesion_cistest.json`（可 grep −0.1418），非捏造。

## §6 待你確認的決策點

1. **診斷（§1 三根因）**是否符合你對「為何現在分析不乾淨」的理解？
2. **目標系統（§2）**範圍對嗎？特別是「characterization 不是 filter」的定位。
3. **路線優先序**：P0 可重現地基**先做**（先有 BRCA2 regression test 當守門，後續改動才安全）——同意嗎？還是你想先看 P1 抽取（5mC/5hmC 分離全位點）或 P3 scan（找更多位點）的結果？
4. **C++ 完全延後**（本輪純 Python）——同意嗎？
5. **重心**：廣度（P3 scan 找更多位點）vs 深度（P2/P3 causation 把 BRCA2 釘死）——你偏哪邊？

確認後我用開發 workflow 從你選的 Phase 開始，每步 git-diff 驗證。
