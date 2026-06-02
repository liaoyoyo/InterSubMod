---
title: 甲基分群拯救 unphase/HP3 與 distant read 的相位 tag — 探索與規劃 scoping
date: 2026-05-31
status: in_progress
task_type: A (pilot / scoping) — partial flag: 探索階段，尚未鎖定 scope，未跑全樣本
phase: brainstorming (design 前置探索)
author: InterSubMod Research (AI session)
partial_flag: true
provenance:
  - workflow_run_id: wf_41624212-f28 (task w6eys2r09)
  - agents: 6 平行 Explore + 1 綜整 (7 agents, ~842k subagent tokens, ~19 min)
  - method: 唯讀探索，未改任何程式碼
related:
  - InterSubMod/docs/CURRENT_FOCUS.md (5/30 audit)
  - memory: project_phase2_cycle1_global_fp_filter (甲基 FP filter DEAD)
  - memory: project_loh_constrained_phasing_discovery (LOH-phasing 活軸)
  - memory: feedback_asm_allele_axis_baseline_confound (ALLELE-axis confound)
  - memory: reference_longphase_getvote_source (getVote 行號)
---

# 甲基分群拯救 unphase/HP3 與 distant read 的相位 tag — 探索與規劃

> **狀態**：brainstorming/scoping。此文件記錄探索發現與待決策清單，**設計尚未鎖定、未改任何程式碼**。
> **partial flag**：探索以既有 master dataset + 單樣本 HCC1395 既有產出為主，未跑全樣本驗證。

## 0. 任務動機（用戶原始需求濃縮）

paired phasing+tag（normal+tumor）結果中，很大一部分 read 被標成 **unphase (HP0)** 與 **HP3** → 可能存在大量 read 被 phasing 認為「難處理 / 資訊不足」。目標：

1. **先觀察 ISM 程式**，確認能否修改達成「清楚改變並紀錄數據 + 整體與個案觀察」，評估後讓用戶確認。
2. **評估 longphase-S 碼級可行性** — 用現有 phase+tag 結果，把低信心 read（unphase/HP3、以及「沒交集到 somatic 但交集到『經過 somatic 的 read』的更遠 read」）用甲基分群，配合有交集 read 與原 read 的 HP tag，判定該歸 **HP n 還是 HP n-1**。
3. **約束**：經過的 read 同時注意是否在 LOH / 多拷貝區；HP3 可能含「無 germline SNP 但經過 somatic alt」或「經過多 germline SNP 但 HP 判斷不一致」的 read，unphase 同理；要有可量化信心數值；注意經過兩 somatic 點 read 的 **Ref-Alt vs Alt-Ref**（HP1-1 多點 read 差異）+ 甲基分群；參考既有 RR/RA/AR/AA subclone 與甲基分群方法。

## 1. 卡關根因（最關鍵結論）

**之前卡住的是「甲基化當 TP/FP filter」這條【具體死路】，不是整個甲基化空間都死。**

三層疊加證據：
- **Phase 2 Cycle 1**：HCC1395 in-distribution ΔF1=+0.02236 → LOSO held-out 反轉 −0.00012、5 樣本平均 −0.00004、Wilcoxon p=0.125 direction-negative = **100% sample-level circularity**。
- **Cycle 3 ablation**：甲基化是 5th-rank vestigial（drop 全 5 個 methyl 只掉 +0.00065 = 3% uplift）。
- **系統根因**：ISM 特徵空間對 TP/FP 標籤無內生判別軸（pure methylation AUC ≤0.58 ceiling，O9/O12/O13/TO-pure 多輪確認）；TO 模式 caller_af 單一特徵 AUC 0.654 即超越全 ISM，FP 集中 low-AF，甲基化無法解耦此 confound。

**本次新任務是【不同方向】，不是那條已死的路**，三點正交性：
- **(a) 目標正交**：死路判別 variant 的 TP vs FP（caller 對不對）；新任務修復 read 的 haplotype 歸屬（read 該歸 HP n 還 n-1）。輸出是 **read tag 不是 variant verdict**。
- **(b) 決策粒度正交**：死路需全基因組 universal threshold（正是 sample-level circularity 來源）；新任務是 **per-region / per-read 局部重指派**，理論上可規避樣本層循環。
- **(c) 訊號用法正交**：死路要甲基化「跨樣本泛化判別 TP/FP」；新任務只要甲基化「在同一 region 內把無 germline-SNP 證據的 read 聚到既有 HP cluster」= **region-local 相對距離**而非 absolute classifier。

**但不是免死金牌**：死路根因瓶頸仍可能咬到新任務 — unphase/HP3 人群極可能同樣高度 low-AF / caller_af-dependent（Phase 2 發現 95.2% lost TP 是 low-AF proxy），且 pure-methylation AUC ≤0.58 ceiling 沒理由因侷限子集就突破。→ **必須先在 HP3/unphase-only 子集做「甲基化內生訊號強度」前置驗證**（見 Phase A gate）。

## 2. 兩條紅線（必須主動避免「滑入」死路）

1. **絕不可把「甲基重指派後 HP tag 變化」回頭包裝成「改善 variant TP/FP F1」** — 一旦這樣 framing 就立即落入已 concluded 死路，需 Productive-Failure C1/C2/C3 reopen 條件。**success metric 必須鎖死在 read tag / phasing accuracy 層級**。
2. **ALLELE-axis confound**：用戶提的 Ref-Alt vs Alt-Ref（HP1-1 多點 read）+ 甲基分群，若用 ALT-vs-REF read 甲基差做重指派證據，會被 germline-het 基線 ASM confound（memory: TP 11.1% < het null 15.2%）。**只有 HP-axis（HP1 vs HP1-1）somatic-controlled**（OR=1.79 modest）。→ 重指派證據軸建議以 HP-axis 為主，ALLELE-axis 僅作描述且強制 germline-het negative control。

**Hard Gate 提醒**（流程約束）：本任務含 C++ 改碼可能性（ISM 與/或 longphase-to-mod）→ 任何實際改碼是 Hard Gate（C++ commit 必編譯、不可包進 workflow、改碼前走 /methodology-audit → /cpp-change）。目前僅探索與規劃，未觸發。

## 3. longphase-S 身分查證（已確認 — 2026-05-31 更正前一輪 agent 誤判）

> ⚠ **更正**：第一輪探索 workflow（wf_41624212）誤判 longphase-S = longphase-to-mod，**此判斷錯誤**。第二輪親自 `ls` + 讀 `01_longphase_s.sh` + 讀 longphase-S CLAUDE.md 後確認真相如下。用戶明確要求「對應的 longphase-S 軟體不要錯誤」，故此節為硬事實，已坐實。

**真正的 longphase-S = `/big8_disk/liaoyoyo2001/longphase-s`**（owner: chenhan112；binary `longphase-s` 已 build，32MB）。其 CLAUDE.md 自述：「LongPhase-S is a bioinformatics tool for somatic haplotagging and tumor purity estimation from tumor-normal paired long-read sequencing data. It extends LongPhase v2.0.0」。

**ISM paired pipeline 確實吃它**：`InterSubMod/scripts/pipeline/steps/01_longphase_s.sh` 呼叫 `${LONGPHASE_S_BIN} somatic_haplotag -s germline -b normal --tumor-snv-file --tumor-bam-file ...` → 產 `${SAMPLE}_tagged.bam`，這就是 ISM 消費的 paired tagged BAM。

**longphase-S 真實結構**（`src/` 模組化，非 longphase-to-mod 的扁平單檔）：
- `src/main.cpp` dispatch 5 命令：`phase` / `haplotag` / **`somatic_haplotag`**（最大模組，本任務核心）/ `estimate_purity` / `modcall`（placeholder）
- `src/somatic_haplotag/`：`SomaticVarCaller.cpp`（2433 行）+ `SomaticHaplotagProcess.cpp`（654）+ `TumorPurityEstimator.cpp`（1190）+ `SomaticBenchmark.cpp`（991）
- `src/haplotag/HaplotagType.h`：定義 `enum ReadHP { unTag=0, H1, H2, H3, H4, H1_1=5, H1_2, H2_1=7, H2_2 }` + `SomaticData` / `ReadHpResult` / `PosBase` struct

**最重大發現 — longphase-S 已內建 inheritance + 量化訊號**（改變設計基礎）：
- **已有 H3 → H1-1/H2-1 inheritance 邏輯**（pipeline 報告提到 `readDistri_beforeInheritance.out` / `afterInheritance.out`）= 正是本任務核心「HP3 read 該歸 HP n 還 n-1」，**不是從零建**。
- `SomaticData` struct 已有大量現成量化訊號：`germlineHaplotypeImbalanceRatio` / `allelicImbalanceRatio` / `somaticHaplotypeImbalanceRatio` / `shannonEntropy` / `pure_H1_1_readRatio` 等 = 用戶要的「可量化 tag 信心」候選。
- `PosSomaticOffsetBase: array<vector<pair<int,char>>,2>`（0:ref, 1:alt）+ `alleleCount` = 正是用戶要的 **Ref-Alt vs Alt-Ref（RR/RA/AR/AA）** 記錄結構。

**其他 longphase 變體**（非本任務目標，僅留存盤點）：`/big8/longphase-to`（baseline）、`/big7/longphase-to-mod`（V5/V6，git 8a90532 2026-05-20）、`/big8/longphase`、`/big8/longphase_methylation_dev`。

**待第二輪 workflow（wf_8a3c1d92）確認**：inheritance 確切決策閾值、甲基是否進 voting（預期否，modcall 獨立）、最小注入點、ISM 端能取得哪些欄位。

## 4. 可行性矩陣（9 能力）

| # | 能力 | 可行性 | 實作位置 | 主要風險 |
|---|------|:---:|---------|---------|
| 1 | ISM 篩 unphase/HP3 + 紀錄數據（整體+個案） | 🟢 | ISM 不改碼：reads.tsv 已有 hp_tag(0/3)+alt_support+mapq；master `all_region_rows.tsv.gz` 已有 hp0/hp3_ratio + A8 per-chr | 低；若要新增 per-read 信心欄位則升 🟡 |
| 2 | ISM C++ 改碼做甲基重指派 | 🟡 | `RegionProcessor::process_single_region()` LabelTest 後、`write_region()` 前；複用 DistanceMatrix + HierarchicalClustering + **LabelTest Stage 4 (test_unassigned_affinity)**；hp_tag 為 mutable string；~50-100 行 | Hard Gate；circular（Stage1 已排 HP3→-1）；建議新增 `reassigned_hp_tag`+confidence 欄位不原地覆寫 |
| 3 | **甲基在 HP3/unphase 子集是否有內生判別訊號（前置 gate）** | 🟡 | 分析腳本：複用 `analyze_methyl_cluster_allele_cooccurrence.py`；HP3/unphase-only 子集算甲基距離能否聚回已知 HP（有 germline-SNP anchor 算 silhouette/affinity AUC） | **紅旗核心**：AUC ≤0.58 ceiling 是系統根因；若此驗證仍 ≤0.58 → 動 C++ 前 NO-GO。**最該先做、最便宜的 gate** |
| 4 | longphase-S 碼級 + 甲基注入 | 🟡 | `judgeHaplotype()` post-vote（V5 行 720-743）對 hpResult==0/33 接 `hts_base_mod_state` 解 MM/ML | longphase 不 parse MM/ML（從零接 htslib）；改 BAM HP:i tag = 高 blast radius（影響全下游）；V5/V6 Layer 1.5 region bias 繼承風險 |
| 5 | per-read tag 信心量化 | 🟢 | 投票邊際 `max/(max+min)` 已算（行 728）；germline 票強度 (HP1+HP2)；MAPQ/BQ 已在 ReadInfo | 語意：原 phasing 信心 ≠ 重指派信心，需分兩個獨立數值 |
| 6 | 兩 somatic 點 read Ref-Alt vs Alt-Ref（HP1-1 多點）+ 甲基 | 🟡 | ISM：ReadParser 已有 per-read alt_support+CIGAR；需擴展記錄 read 跨多 somatic 位點有序 ALT/REF（目前只 region 聚合）；複用 `compare_subclone_validation.py` Fisher | 無現成 read-level 多點 co-occurrence 結構；ALLELE-axis confound；跨兩點 read 稀疏（相鄰 TP pair median common reads 僅 17）統計力受限 |
| 7 | LOH/多拷貝區 constraint | 🟢 | master 已有 loh_coverage_pct/subtype/Coverage_Multiple；HCC1395 A8 per-chr | COLO829/H2009/H1437/HCC1954 缺 A8 級統計需補；LOH-inner NG=2 強制 same-hap 會扭曲 ALT/REF 比例（勿誤當甲基訊號） |
| 8 | **distant read 重指派**（沒交集 somatic、交集到「經過 somatic 的 read」的更遠 read） | 🔴 | 需新建：無 read-overlap graph 顯式結構（DistanceMatrix 是甲基 read×read 非基因組 overlap） | **最高紅旗**：「distant」定義未定；transitive inference 每跳累積誤差；跨 PS block 在 longphase 本就強制 unphase（對抗其保守設計）。建議切獨立後期階段或先純可視化 |

## 5. 建議任務拆解（4 phase, ~7-12 工作天）

- **Phase A — ISM 觀察與 unphase/HP3 數據刻畫（不改碼）**：master+A8 整體刻畫 HP0/HP3（見樹也見林四層）；拆 HP3/unphase 子類（無 germline-SNP+somatic-alt vs 多 germline-SNP 不一致）；**前置可行性驗證（GO/NO-GO gate）**：HP3/unphase 子集甲基 affinity/silhouette AUC；量化子集 caller_af 分布（confound 前檢）；標 region LOH/CN。→ 主回合 / 單 analysis agent；gate 結果回主 agent 等用戶確認。
- **Phase B — 碼級可行性評估（methodology-audit，不改碼）**：確證 longphase-S 身分；ISM vs longphase 注入點碼級評估；信心量化雙數值設計；產 methodology-audit 決策文件供選方案。→ 單 methodology-reviewer agent（read-only）。
- **Phase C — 重指派原型（Hard Gate：含 C++ 改碼）**：用戶選方案後走 /cpp-change；優先 ISM in-memory+TSV（reassigned_hp_tag+confidence，不碰 BAM）；unphase/HP3 局部重指派用 HP-axis 甲基 affinity + LOH/CN gate；（可選後期）RR/RA/AR/AA 多點 read；distant read 列獨立後期先可視化。→ 主 agent 編排（Hard Gate）。
- **Phase D — 驗證（防 overfit）**：重指派正確性（vs anchor/ground-truth 一致率 + phaseable read 提升）；**LOSO 樣本層交叉驗證**（若用 sample-level 權重）；LOH/CN 分層；HP-axis vs ALLELE-axis + germline-het negative control；run-evaluator tier；明文確認 metric 全程是 phasing accuracy。→ 跨 7 樣本 fan-out（無 Hard Gate）；tier/NO-GO 判定回主 agent。

**關鍵 gate 順序**：Phase A 的「HP3/unphase 子集甲基判別力前置驗證」是最便宜、最該先做的 GO/NO-GO 節點 — 若 AUC ≤0.58，應在動任何 C++ 前停下與用戶確認，避免重蹈 pure-methylation ceiling。**首版強烈建議鎖 ISM in-memory 原型 + 僅 unphase/HP3 局部重指派**，longphase BAM 回灌與 distant-read 跨中介傳遞留到證實概念後的獨立後期階段。

## 6. 用戶已拍板決策（2026-05-31 AskUserQuestion）

1. **success metric** = **phasing rescue + 觀察對 LOH-constrained phasing 活軸的增益**（嚴格分離，絕不回頭包裝成 variant F1）
2. **改哪邊** = **先 ISM 原型**（低風險可逆）；用戶補充：聚焦 **paired 模式 + 正確的 longphase-S 軟體**（非 TO），這是基本方法學改動而非 F1 提升問題；paired 有結果後再評估是否套用到 TO
3. **首版範圍** = **unphase/HP3 局部重指派（核心）+ per-read tag 信心量化 + RR/RA/AR/AA 兩點 read 區分**；distant read 留後期；每步拆最小單元 + commit + 驗證 + 可回溯（**TDD 風格**）
4. **驗證 scope** = **HCC1395 pilot**（partial flag）
5. **甲基證據軸** = HP-axis 為主 + ALLELE-axis 描述用（避 confound）
6. **longphase-S** = `/big8_disk/liaoyoyo2001/longphase-s`（已查證更正，見 §3）
7. **新增要求**（用戶）：**抉擇方法時要紀錄觀察數據 + 可量化的抉擇依據與決定**，供檢核理解驗證 → 落實為 method decision ledger（見後續 design spec）

## 6b. ⭐ Phase A 現成 ground truth + 已驗證量化基線（2026-05-31 親自重算）

> 用戶要求「方法抉擇要紀錄觀察數據 + 可量化依據供檢核」。以下數字均親自從 longphase-S 輸出檔重算、欄位逐一驗證、守恆檢查通過。

**資料源**：`/big8_disk/liaoyoyo2001/longphase_somatic_output_ssrs/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter_*`
源 BAM = `ONT_5khz_simplex_5mCG_5hmCG`（**甲基資料已在 BAM 內，無需重跑 basecalling**）；percentageThreshold=0.6（HP3.out header）。

**現成可直接觀察的 per-read / per-position 輸出檔（Phase A 零改碼即可分析）**：
| 檔案 | 內容 | 對本任務用途 |
|------|------|------|
| `_totalRead.out`（1.07 GB）| per-read：germlineVarSimilarity / deriveByHpSimilarity / germlineSnpCount / tumorSnpCount / Haplotype / somaticVariant,HP | **per-read 信心值現成**（germlineVarSimilarity = 投票相似度）|
| `_readDistri_beforeInheritance.out` / `_afterInheritance.out`（各 9.7 MB, 120,834 位點）| per-position：DeriveHP / DeriveHPsimilarity / somaticBase_{H1-1,H2-1,H3} / {HP1,HP2,HP1-1,HP2-1,HP3,untag}read / 各 ratio | **inheritance 前後對照現成** |
| `_readDistri_Scaller_derive_by_H1_H2.out`（25,653 位點）| DeriveHP=H0 的難繼承位點子集 | **位點級救援目標清單** |
| `_alleleCount.out`（39,447 變異）| per-variant REF/ALT/DEL count | RR/RA/AR/AA 基礎 |
| `_HP3.out`（21 MB）| somatic calling log + filter 參數 | HP3 個案 |
| `_crossHighConSnpRead.out`（247 MB）| 跨 high-con SNP 的 read | 跨多點 read 線索 |

**已驗證量化基線 — longphase-S 既有 inheritance 效果（HCC1395, 120,834 somatic SNP 位點，somaticBase read-instances）**：

| somaticBase read-instances | BEFORE inheritance | AFTER inheritance | 變化 |
|---|---:|---:|---:|
| H1-1（繼承到 HP1 family）| 1,296,335 | 2,595,184 | **+1,298,849** |
| H2-1（繼承到 HP2 family）| 1,007,638 | 1,531,667 | **+524,029** |
| H3（卡住、未繼承）| 4,343,202 | 2,520,324 | **−1,822,878** |

- **守恆驗證 ✓**：H1-1+H2-1 增加量 (1,298,849+524,029=1,822,878) = H3 減少量 (1,822,878)，完全相等。
- **longphase-S 自己已救 42% 的 H3** somatic read-instance 進 HP n-1（4,343,202 → 2,520,324）。
- **位點級 DeriveHP 分布**（120,834 somatic 位點）：H0（無法定相，最難）=**65,033 (53.8%)** / H1=33,947 / H2=21,854；其中 25,653 個 H0 位點屬 `derive_by_H1_H2` 子集。
- **read 級分布（per-position read tally，inheritance 前→後）**：HP3read 4,953,587 → 2,920,520；**untag(unphase) read 1,777,511 → 1,777,511（完全不變）**。

**⭐ 關鍵發現 1 — longphase-S 的 inheritance 完全不碰 unphase read**：untag read 在 inheritance 前後一模一樣（1,777,511）。longphase-S 只把 H3（有 somatic 證據但未定 germline family）往 H1-1/H2-1 救，**對 unphase read 零處理**。→ **unphase 是 longphase-S 主動放棄的族群，是甲基分群最大、最乾淨的機會空間**（不與既有 inheritance 邏輯衝突）。

**⭐ 關鍵發現 2 — 本任務不是「從零建 inheritance」**，而是**「在 longphase-S 既有 inheritance 已盡力後，對 (a) 它放棄的 unphase read、(b) 它救不動仍卡 H3 的 read、(c) 65,033 個 H0 難定相位點，用甲基分群再救一層」**。Phase A 第一步零改碼就能量化這三類待救援族群的規模、caller_af 分布、LOH/CN 位置、甲基覆蓋率。

- ⚠ 註：上述多為 read-instance（一條 read 跨 N 位點計 N 次），非 unique read；Phase A 需另從 `_totalRead.out` 算 unique-read 層級。

## 7. 策略分岔已定案（2026-05-31，第二輪深度分析 wf_7cbe71fc + 用戶拍板）

第二輪深度分析（8 agent，~186 萬 token，逐位元組 + 雙重對抗驗證 confirmed）確認 longphase-S 既有 inheritance 細節後，用戶拍板：
- **MD1 介入層 = C（先 ISM post-hoc 旁路驗證，過 GO 門檻才開 C++）**
- **Phase A = 純 Python 旁路分析**（連 ISM C++ 都不改，最快驗證判別力）

→ 完整設計、method decision ledger（10 項）、Phase A 詳細規劃見獨立 design spec：
**`InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_unphase_HP3_phasing_rescue_design_02.md`**

本 scoping doc 自此凍結為「探索證據紀錄」，後續設計與執行以 design_02 為準。

---
*探索原始輸出：workflow wf_41624212-f28 transcript dir（6 Explore + 1 綜整 agent 結構化 JSON）。本文件為 brainstorming 階段 scoping，設計鎖定後另寫 design spec。*
