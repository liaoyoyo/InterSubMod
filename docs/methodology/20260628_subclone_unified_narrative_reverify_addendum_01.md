<!--
建立時間: 2026-06-28
類型: 獨立大規模重驗 addendum — 對 20260627_subclone_unified_verified_narrative_01.md 的數字 / 框架 / deck 第三方查核
狀態: in_progress（addendum；不取代 20260627 unified narrative，作為其 verified 補丁 + deck reframe 依據）
build_branch: research/subclonal-reconstruction-202606
provenance: 本檔由本 session 獨立 workflow wk0n0tjxk（run wf_398e23b3-cf1；9 agent = 5 數字叢集 + F1 框架 + F2 跨doc + F3 deck關係 + 1 綜合）對「實際 data 檔」親自重算產出；非引用前述報告自述。
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/{sm_summary.json,sm_completeness_ledger.json,regions.tsv,sm_region_integration.json,region_shape_distribution.json,clean_subset.json,sm_hp_contribution.json,sm_ccf_tiers.json,sm_locus_master_summary.json,sm_phaseset_extension.json,sm_methyl_corroboration.json,sm_methyl_reextract_ALL.json,sm_methyl_genetic_concordance.json,chr17_subclone_data.json,chr17_tree_data.json,sm_configuration_census.json}, /big7_disk/liaoyoyo2001/wt-merge-summary/docs/methodology/20260627_clone_subclone_integrated_report/data/{sm_methyl_sufficiency_audit.json,sm_single_locus_methyl.json}
-->

<!-- provenance-verified: 本檔每個 metric 數字皆由 workflow wk0n0tjxk 之 agent 從上列 data 檔親自 grep/python 重算命中（39 項檢查 35 MATCH）；2 個 MISMATCH 與 2 個 UNVERIFIABLE 已逐一標記為「不可寫死」。 -->

# Subclone 統一敘述 — 獨立大規模重驗 addendum（2026-06-28）

> **本檔職責**：對 `InterSubMod/docs/methodology/20260627_subclone_unified_verified_narrative_01.md`（06-27 統一敘述）與 `20260628_reconstruction_model_verification_01.md`（06-28 模型驗證，另一 session 產）做**第三方獨立重驗**——親自重算 headline 數字、抓捏造、收斂跨-doc 衝突、定 deck reframe。
> **不取代**上述兩檔；作為其 verified 補丁。框架 SoT 仍指向 06-27 unified narrative。
> **敘述框架**：Verdict-Pyramid + scientific-rigor §2（claim 標 L1-L5）。

---

## §1 框架裁決 — CONFIRMED（高信心）

**verdict = `SUPPORTED_WITH_CAVEATS`**（框架信心 = **高**；個別 headline 數字信心 = 中）。

**39 項檢查 → MATCH 35 / MISMATCH 2 / UNVERIFIABLE 2**（90% 為 agent 從 `regions.tsv` / per-locus TSV / `sm_summary.json` / SEQC2 bed **獨立重算命中**，非僅引用 JSON）。分叢集：A 骨幹 7/7 · B HP 3+1U · C VAF/CCF/PS/CN 11+1U · D 甲基 9/9 · E chr17/碼層 5+2M。

**你的主軸框架「somatic sSNV 共現 = 非循環骨幹 / HP = 鑑別器（非確認器）/ 甲基 = off-ladder bounded-auxiliary」有碼級支撐：**

- **[L1] 骨幹非循環**：`classify()`（`sm_linkage_genomewide.py:103-112`）關係**純由 2×2 零格拓樸指派、完全不碰 VAF / read 數**；HP tag 直接讀 BAM longphase-S 既有 tag（非從受測 sSNV 衍生）→ 非循環成立。
- **[L1] HP = 鑑別器**：互斥 `mutual_excl` DEPLETED **0.86×**（診斷力所在）vs nested / co_linked / independent **1.7–1.87×**（85–94% same-HP = 區域背景非克隆特異）；HP ablation 移除 **57%** allelic 假 subclone。
- **[L1] 甲基 = corroborate-not-detect**：`new_partitions_detected_by_methylation = 0`、read×read PERMANOVA 754 區僅 **1** 可算、**cis-control `hp_control_evaluable = 0/740`（逐 1610 列 TSV 確認 structural zero）→ subclone-specificity UNDETERMINED（非已驗陰性）**。
- **[L1] 七紅線全 GUARDED**（甲基偵測器 / LOH 當乾淨 subclone / genome-wide tree / 對手缺檢定 皆未 committed）。
- **[L1] 誠實上界**：骨幹 **53–69% 落 CN-gain / LOH**（multiplicity confound）→ 該子集非循環性 **UNDETERMINED**；宣稱須限 **CN-clean ∩ non-segdup**、**regional（≤ read-span）**、**⭐3 single-pipeline 封頂**。

---

## §2 數字安全台帳（對外敘述前必查）

| 類別 | 數字 | 敘述指引 |
|---|---|---|
| ✅ **VERIFIED（可講）** | sSNV 宇宙 35,332 = TP 30,490 + FP 4,842（sum✓）· 7,143 區 · 677 full_tree · ~205 CN-clean full_tree（可由 `regions.tsv` filter `cn∈{loh,neutral}&full_tree` 復現＝真衍生非硬編）· HP 0.86× / 1.7–1.87× / 移除 57% · 甲基 49/740 = 6.6% / new_part = 0 / cis-control 0/740 · CN-gain 52.8%（12,569）· PS reliable 92.7% · VAF 69.8 / 5.7 / 24.5 · CCF GMM best_n=3 · Fisher≈permutation 97.7% | 直接用，grep-able |
| 🔴 **MISMATCH（捏造，必修不可用）** | chr17 VAF 第三值「**0.18**」→ 實 **0.47/0.48**（`chr17_subclone_data.json` snv_stat 重算）；「0.18」+ 不存在的 γ∥α 第 4 sibling 僅硬編於 `sm_summary.py:156`。 · 成環「**19/23**」CN-gain → 實 **18/23**（對 SEQC2 bed first/majority/any 三規則一律 18 gain+5 loh；off-by-one 硬編）。 | 講「VAF 0.82 > 0.47/0.48」、成環「18/23」；`0 violations` 部分為真（3 對 off-diagonal 皆乾淨零格） |
| ⚠ **UNVERIFIABLE（無來源，不可寫死）** | 「HP3 污染 **6.5%**」（`sm_hp_contribution.json` 無 `sibling_sameHP` key 也無 6.5 字面；HP3=1317 但任何乾淨分母湊不到：1317/19497=6.75%、1317/20815=6.33%）· 「**77.8%** segment-level CN-gain」（無 data 檔，SAVANA per-segment 不在資料集） | 對外不寫；77.8% **永不與 52.8% 並列** |

---

## §3 跨-doc 衝突收斂（含用戶 2026-06-28 拍板）

| # | 衝突 | 收斂決定 |
|---|---|---|
| 1 | chr17 canonical = **3**（data）vs **4**（headline 硬編）sSNV | 🟡 **先擱置，標「待 somatic 定義統一」**——因 `is_somatic` 碼級分裂（全基因組 `normal VAF<5%` vs chr17 per-read 嚴格 `==0`，差 **2,230 sSNV / 9.05%**）牽動「到底幾 sSNV」的根；統一前**不寫死 2/3/4** |
| 2 | structure 確認共現結構區 = **3,820**（53%，排除 858 單 lineage）vs **4,678**（65%）| ✅ **對外口徑採 4,678（65%）**（含單 lineage）；honest 註：含 858 個僅單一 lineage 的 co_linked 區（06-27 防禦口徑曾排除為 3,820）|
| 3 | CN-gain **52.8%**（somatic-weighted）vs **68.9%**（linked 子集，實分母 linked_somatic 20815；linked 全體 21554 → 66.6%）vs 77.8%（segment-level，UNVERIFIABLE）| ✅ **錨點 = 52.8% somatic-weighted**；68.9%/66.6% 用時**必標分母 = linked 子集**；77.8% 不可用、永不與 52.8% 並列（§13 跨分母禁並列）|
| 4 | 06-28 揭碼層問題影響 06-27 headline 但 06-27 無 | 記入 §5 待修：somatic 定義 ==0 vs <5%（2,230 sSNV）· `classify()` `aa<2` vs off-diagonal `==0` 不對稱（50.61% nested 有 off-diag==1 可翻 co_linked → 系統性低估 lineage 桶）· HP3('3') 未排除 |
| 5 | provenance 脆弱性 | `regions.tsv` / `sm_region_integration.json` 現已 materialise 至主樹 `_assets/20260627_subclone_4axis_teaching/data/`（06-28 已可主樹 drill），但**數字 SoT 仍依賴 pending-merge branch 5308d9e**；FF-merge 進 trunk 前須留意清 worktree 風險 |

---

## §4 06-23 週報 deck 處置 = RE-SCOPE（非 REFUTE）

verdict = `PARTIALLY`（部分過時，主訊息方向反轉須改，多數方法論頁仍成立）。

**保留（新框架仍認同）**：
- **S4** 方法 SOUND（read×read 距離 → PERMANOVA 成熟非 ad-hoc；僅證「存在性」非「鑑別性」）
- **S6** 判別機制（a-priori 標籤 CramérV≥0.7 交叉 = 唯一合法非循環鑑別）
- **S7** chr2:18M（locus-subclone 未證、18,086,020 落 SEQC2 HC 空隙 unevaluable、天花板 L2）
- **S8**（對齊 germline = cis-ASM 非 subclone、分類 ≠ 確認、normal cis-control = 未測的決定性 gate）
- **S9** 收束（已含 genetic-first 種子）

**必改（甲基-主秀方向被反轉）**：
- **S1 封面 + S2 BLUF**：「甲基為主秀的位點分類 pipeline」→ 改 **genetic 起手 / 甲基收尾**（整頁 reframe）
- **S3**：甲基距離→UPGMA→切群當「the pipeline」**過時**（= 循環 double-dip、recover 0 新 partition、754 區僅 1 可算、6.6% 弱佐證）→ 降為輔助節；主流程改 **genetic 骨幹**（2×2 共現 census → four-gamete → nesting → 6,146 component → 7,143 區）
- **S5**：「S1 對齊 18.41% 當主要發現」→ 降為 characterization 副結果，明標 circular / double-dip 限制

**新增（缺席的 genetic 骨幹頁群）**：35,332 sSNV → 7,143 區 → 677 full_tree → ~205 CN-clean 子集；four-gamete / perfect-phylogeny；CN-gain 52.8% multiplicity confound；HP 只在互斥 0.86× DEPLETED 有診斷力。

**對齊**：scope 數字漂移（deck `34,736 / TP 30,077` = **舊 build**，vs 新框架 `35,332 / TP 30,490 + FP 4,842`；§13 跨 build 禁並列 → 統一錨點 + 標 provenance）；S4 reclassify「C++ 53.9% 待落檔」可更新為**已落**（commit b6dde2b）。

---

## §5 open gaps（依 load-bearing 排序）

1. 🔴 **對外不可宣稱**（七紅線守住但前提未測）：甲基為 subclone 偵測器（new_part=0）· LOH 當乾淨 subclone（corroborated 41/49 正落 LOH = 假陽最高發處）· genome-wide tree（僅 regional ≤ read-span）· 對手缺顯著性檢定（ASMS 有 permutation 檢定，口徑必改）。甲基 subclone-specificity = **UNDETERMINED**（不可寫成陰性、也不可寫成佐證）。
2. 🔴 **最 load-bearing 下一步 = `T-GATE-GB` normal cis-control**（cis-control 0/740 structural zero = 判別力的決定性 gate；未測前甲基鎖在 bounded-auxiliary）。
3. **碼層待修**（影響桶計數但結論未翻）：`classify()` `aa<2` vs off-diagonal `==0` 不對稱（系統性低估 lineage 桶）· Fisher over-dispersion → beta-binomial · somatic 定義 `==0` vs `<5%` 統一（差 2,230 sSNV 影響宇宙計數 + chr17 sSNV 數）· 清 `sm_summary.py:156` 硬編「0.18」+ 不存在的 γ sibling。
4. **骨幹非循環性**在 ≥53% 骨幹（CN-gain / LOH）為 UNDETERMINED → 宣稱限 CN-clean ∩ non-segdup；clean 桶把 8,697 LOH 計為 clean = 寬鬆定義選擇（LOH-unmask cis-ASM confound）非已驗 confound-free。
5. **確認層天花板 ⭐3 single-pipeline**；subclone 確認黃金標準 = single-cell / multi-region，單 bulk 只能 characterize；需 COLO829（normal 缺）補單樣本。
6. **Provenance**：FF-merge `5308d9e` 進 trunk 前，trunk 內無原始 drill-down；須先回流避免清 worktree 即證據全失。

---

## §6 知識追溯（§8.4）

- **重驗工具**：workflow `wk0n0tjxk`（run `wf_398e23b3-cf1`），9 agent，唯讀對抗式重算，1,043s，784K subagent tokens。
- **本檔關係**：addendum；框架 SoT = `20260627_subclone_unified_verified_narrative_01.md`；模型/碼層細節 = `20260628_reconstruction_model_verification_01.md`；任務節點 = `T-SUBCLONE-VERIFIED`。
- **confirm 範例**：`python3 -c "import json; d=json.load(open('docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_summary.json')); print(d)"` 查骨幹計數；chr17 VAF：`python3 -c "import json; d=json.load(open('docs/methodology/_assets/20260627_subclone_4axis_teaching/data/chr17_subclone_data.json')); print(d.get('snv_stat'))"`。

---

## §7 「表面矛盾」狀況的驗證與定義處理原則（co_linked × diff-HP case study）

> 起因：2026-06-28 教學討論質疑 co_linked（同分子）卻標 diff-HP（不同單倍型）是否真矛盾。逐對重驗後確認**非矛盾**，並據此訂出可複用的處理原則。本節 = 推論鏈紀錄 + 驗證 recipe + 定義處理規則。

### §7.1 案例推論鏈與結論（verified）

**觀察**：co_linked 對中 **12.7%（1,496/11,750，canonical coread≥6∩both-somatic 層）** 標 diff-HP（coread≥2 寬層為 15.8%）。表面矛盾＝同一條分子怎會跨兩單倍型。

**推論鏈**：
1. **證據性質不同**：co-linkage = **直接觀測**（2×2 表 RA=AR=0，直接從 read 讀到兩 ALT 同分子）；HP tag = **衍生推斷**（longphase 用鄰近 germline 雜合 SNP 推每條 read 的單倍型）。
2. **HP tag 已知不可靠條件**：CN-gain multiplicity（≥3 拷貝使二元 HP 指派歧義）/ HP3 未定相 / phase-switch。
3. **驗證**：cross-tab `lists/co_linked_diffHP.tsv`（n=1,496）的 `cn_a` × `hp_a`/`hp_b`。

**結論（L1，grep-able）**：

| 指標 | co_linked diffHP (1,496) |
|---|---|
| CN-gain | **1,314 = 87.8%** |
| LOH | 166 = 11.1% |
| loss | 13 = 0.9% |
| **CN-neutral** | **3 = 0.2%（= 0.026% of 全 11,750 co_linked）** |
| 含 HP3('3') | hp_a=3:303 / hp_b=3:371（~25–45% 涉及未定相）|

→ **非矛盾**：diff-HP co_linked 幾乎全落在「HP tag 本就不可靠」之處（CN-gain 87.8% + HP3 未定相）；**CN-clean 乾淨區 co-linkage 與 HP 幾乎從不衝突（3/11,750）**。

**裁決（L1/L2）**：`co-linkage（直接）> HP tag（衍生）`；HP 鑑別力宣稱限 **CN-clean ∩ phased** 子集。此發現**強化**框架（骨幹 > 輔助、CN-gain 為 precondition confound）而非削弱——先前「🔴 16% 內部矛盾須查」的框定已撤回。

### §7.2 通用「表面矛盾」驗證 recipe（可複用於其他軸）

當「直接觀測」與「衍生標註」表面衝突時，按此 5 步：
1. **分清證據性質**：哪個是直接 read-level 觀測、哪個是衍生推斷。
2. **列衍生標註的已知不可靠條件**（HP→CN-gain/HP3/phase-switch；VAF→CN；甲基→cis/double-dip）。
3. **cross-tab 衝突案例 × 不可靠條件**（用逐對 TSV / 逐位點資料，非報告自述）。
4. **判定**：衝突是否集中在不可靠條件？乾淨子集殘留多少？
5. **結論**：殘留 ≈ 0 → 非真矛盾（衍生標註不可靠，信直接觀測）；殘留顯著 → 真異常須深查。

### §7.3 定義處理規則（落地，對外/SoT 一律遵循）

- **證據階層公理**：`直接分子觀測 > 衍生推斷標註`；衝突時信前者、標後者「不可靠」。
- **same_hp / diff_hp 條件式可靠**：`diff_hp` 在 CN-gain/HP3 **不代表** trans/不同 lineage，只代表 **HP tag 在該處不可靠**。
- **scope 紀律**：HP 鑑別、VAF magnitude、甲基 corroboration 的宣稱各限其「可靠子集」（HP→CN-clean∩phased；VAF→CN-clean；甲基→需 cis-control）。
- **輸出補強建議**：逐對資料應帶 `cn_a` + HP-reliability flag。現 `lists/co_linked_diffHP.tsv` 已有 `cn_a/hp_a/hp_b`（可重現本節），建議補 `hp3_involved` + `cn_b` + `phase_set` 一致性欄，讓「表面矛盾」可被機械 audit。
- **confirm 指令**：`python3 -c "import csv; from collections import Counter; r=list(csv.DictReader(open('docs/methodology/_assets/20260627_subclone_4axis_teaching/data/lists/co_linked_diffHP.tsv'),delimiter=chr(9))); print(len(r), Counter(x['cn_a'] for x in r))"`
