<!--
建立時間: 2026-07-11 16:48 +0800
目標: read-group C、candidate tree T、topology shape 與 VAF ordering 報告實作 living document
處理範圍: HCC1395 教學案例 + chr1-22 × 7 datasets historical layered-v2 census
cycle_id: 20260710-2244-layered-reconstruction-v2 (derivative report correction)
spec_id: read_group_C_tree_T_topology_report
status: report_complete_partial_evidence
advisory: on
build_branch: research/subclonal-reconstruction-202606
build_commit: 6067568637088838a9f518955e41d222f057e4f1
關聯檔案:
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/pre-decision-audit.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/implementation-notes.md
-->

# Implementation Notes：read-group C、tree T 與 topology report

> 本檔記錄使用者在對話中鎖定的定義、舊口徑偏離與仍未解的證據邊界。

## 🔵 設計決定

### 2026-07-11 21:04 — strict >0.50 × top-tie × Topo census（完成）

- **Status**：Completed；HCC1395 probe PASS → 全 7 rows PASS → 兩位 post-fix reviewers PASS / 0 blocker。
- **研究問題**：重算原 `top weight<0.60` 的 14,497 regions，區分 unique first、tied first/same Topo、tied first/different Topo，並比較 0.50 與 0.60。
- **Scope**：chr1-22 × 7 dataset rows；historical layered-v2 engineering snapshot；不替換 clean layered-v3。
- **主 tie contract**：依 reviewer 建議，改用整數 REF/ALT counts 推導 `Fraction` exact raw score，完全相等才算 exact co-top；另以 absolute score gap `1e-12`、`1e-9`、`1e-6` 報 near-tie 敏感度。
- **必要邊界**：`>=0.50` 單獨不保證唯一；每個 region 必須同報 `n_top_joint`、`n_top_shapes_joint`、top/second weight 與 direction flag。
- **使用者門檻定案**：主表使用嚴格 `top weight >0.50`；`>=0.50` 僅作邊界對照，避免 `.50/.50` 並列被標成唯一。
- **驗證門檻**：舊 0.60 結果精確重現；14,497/29,053/39,885 守恆；candidate count 與 canonical shape count逐 unit吻合。
- **全部可評估區域**：28,183 = unique first 19,092（67.74%）+ tied/same Topo 6,886（24.43%）+ tied/different Topo 2,205（7.82%）。
- **原 14,497 below-0.60 拆分**：unique 5,406（37.29%）+ tied/same 6,886（47.50%）+ tied/different 2,205（15.21%）。
- **原始 Topo 交叉**：`T>1,Topo=1` 的 5,498 = 1,570 unique + 3,928 tied/same + 0 tied/different；`T>1,Topo>1` 的 8,999 = 3,836 + 2,958 + 2,205。
- **0.50 邊界**：原 14,497 中 strict `>0.50`=2,261（全部 unique）；`=0.50`=2,638（全部 tied/same Topo）；`<0.50`=9,598。
- **全 ambiguous 比較**：legacy `>=0.60`=13,686/29,053（47.11%）；strict `>0.50`=15,947/29,053（54.89%），增加 2,261／7.78 pp。
- **Temperature sensitivity**：strict `>0.50` 在 T=0.025/0.05/0.10 為 17,458/15,947/13,884；因此 direct Fraction-score rank/tie 作主描述，weight threshold 只作 concentration 輔助。
- **數值穩定性**：exact class 在 absolute score gap `1e-12/1e-9/1e-6` 相同；只代表 numerical tolerance 穩定，不代表 read-sampling uncertainty。
- **方向修正**：1,738 exact-top candidates 無 ancestry comparisons，已輸出 `direction_evaluable=False`、`direction_consistent=NA`；legacy vacuous-true 欄只用於重現舊結果。
- **方法缺口**：21,646/35,457 ranked units 的候選 comparison count 不同；normalized-score／likelihood 與 read/molecule bootstrap 尚未完成。
- **輸出規模**：47,377 region rows、68,544 unit rows、385,790 candidate ranking rows；56/56 checks PASS。
- **獨立 QA**：逐 region TSV group-by 重算三類、strict >0.50、逐 dataset、原 Topo、cardinality bounds、score_rank 與 direction-NA 全 PASS；科學 reviewer 判定 direct ranking 主描述與 strict >0.50 次要定位正確。

**正式執行命令**：

```bash
python InterSubMod/research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_top_tie_census.py \
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 \
  --input-manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/input_manifest.json \
  --corrected-report InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json \
  --legacy-vaf-summary InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_region_resolution_summary.json \
  --method-script-dir InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts \
  --output-json InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_census.json \
  --output-region-tsv InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_regions.tsv \
  --output-unit-tsv InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_units.tsv \
  --output-candidate-tsv-gz InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_candidates.tsv.gz \
  --output-summary-tsv InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_summary.tsv \
  --output-checks-tsv InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_checks.tsv \
  --output-report-md InterSubMod/research/20260711_read_group_C_tree_T_topology_report/20260711_VAF_top050並列樹與拓撲完整重算_01.md
```

**最終 SHA-256**：script=`6e9786be070e361601cdb1f71d76c2d1060f0d041506f9f6265383a4c9f47005`；JSON=`c0879b878899c0abfd635d0ae77c9153f86d1991ed3e04dc8b04ec66dfc5d3a7`；region TSV=`a6649bf12dba686fad3d74d33243e16f4ec2b74361cac5a6b6d09eb09b794996`；candidate TSV.GZ=`6c70f4a1114e163c415736f1542e8e550e7c9b7240c36df666b9d9e211e33493`；report=`9f6f7b4301ac86ca9bb720d04a9269bc13b788859ba60d2f3ee0ee4c6f873841`。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: 依使用者後續修正，主表使用 strict `top weight >0.50`，並把原 14,497 regions 拆成唯一第一、並列第一且 Topo 一致、並列第一且 Topo 不一致；direct exact-score ranking 作主描述。
**DO NOT change**: 不把 0.50 升格為 calibrated probability、biological threshold 或真樹確認；`>=0.50` 只作 `.50/.50` 邊界對照。
**Rationale**: 使用者要求用「最高順位候選樹有幾棵」取代只有 PASS/FAIL 的粗分類。
<!-- END USER-SPECIFIED -->

### 2026-07-11 20:03 — Slide 8 改為「判定依據 → 類別內分流」圖

- **Status**：Accepted；artifact/browser/no-JS/print QA PASS。
- **舊圖問題**：原本上方只顯示 complete regions 的 `T=1,Topo=1 / T>1,Topo=1 / T>1,Topo>1` 組成，下方再另列 VAF 數字；沒有用同一圖說明哪些 region 進 gate、依據什麼通過、以及真正改變的欄位。
- **新圖起點**：`39,885 complete = 10,832 T=1/Topo=1 已唯一 + 29,053 T>1 多解區域`；只有 29,053 進入 VAF gate。
- **Unit gate**：區域內 `T>1` primary HP unit 需 top relative weight ≥0.60；`T=1` HP unit 無排序歧義，自動通過。
- **Region gate**：所有 primary HP units 都通過；雙 HP 為 AND，不用 weight product。
- **圖表分流**：以原始 `T>1,Topo=1` 與 `T>1,Topo>1` 為兩列，各自分成 `Region gate 通過（各 T>1 HP unit top≥0.60） / Top<0.60 / 缺 VAF or mismatch / recurrence 未評估`；每列加總回 11,144 與 17,909。
- **不變項**：VAF 只新增 candidate-priority status，不刪除候選 T、不改寫原始 Topo class、不稱拓撲已矯正或確認。
- **守恆**：`11,144=5,417+5,498+135+94`；`17,909=8,269+8,999+607+34`；合計 `29,053=13,686+14,497+742+128`。
- **Post-fix 精修**：speaker notes 將 10,832 明確標為 `T=1/Topo=1 regions（不進 VAF 排序）`；完整 gate 圖例改成 2×2 排列，避免投影時重疊。
- **獨立複核**：visual/data reviewer 與 scientific wording reviewer 均 `PASS / 0 blocker`；1920×1080 人工檢視及 browser QA 均無裁切或溢位。
- **最終 HTML SHA-256**：`84718ebbdc3b99459c85b80630c14a0ffad120478e3f71458e7d8122d3b52c72`。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: 第 8 頁圖示必須直接呈現「判定依據」與「真正改變的狀態」，不能只把結構組成和 VAF 數字並排。
**DO NOT change**: 不得用圖形暗示 VAF 將 `Topo>1` 矯正成已確認的 `Topo=1`。
**Rationale**: 使用者於 2026-07-11 指出原 Slide 8 圖示沒有正確解釋改變依據。
<!-- END USER-SPECIFIED -->

### 2026-07-11 18:34 — VAF-prioritized region-level 收斂規則

- **Status**：Accepted for descriptive PPT comparison；不是 biological confirmation。
- **決定**：只在 region 的所有 primary HP units 都 `analysis_candidate_set_complete=true` 時評估。`T_unit=1` 自動不需排序；每個 `T_unit>1` ambiguous HP unit 必須以全候選 read-AF scoring 得到 top relative weight ≥0.60。雙 HP region 兩邊皆通過才算 `VAF-prioritized`，不以 weight product 取代 AND gate。
- **排除／未評估**：缺 VAF、candidate mismatch、top weight <0.60、或含 `recurrence_required,T>1` unit 的 region 不算收斂。
- **結果**：結構上 ambiguous 的 29,053 complete regions 中，13,686（47.11%）可取得 region-level VAF-prioritized top；14,497 低於門檻、742 缺 VAF、128 為 recurrence scope-blocked。
- **依原始結構類別**：`T>1,Topo=1` 為 5,417/11,144（48.61%）；`T>1,Topo>1` 為 8,269/17,909（46.17%）。後者只表示從既有候選 shapes 中選一個 prioritized shape，不表示真實 topology 已確認。
- **同源一致性**：13,686 prioritized regions 中 13,674 的所有 winner 通過同一 read-VAF direction check；此檢查不是獨立驗證。
- **Evidence tier**：L2 ⭐⭐⭐⭐（決定性工程重算）；biological true topology 仍為 L5 ⭐。
- **來源**：`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_region_resolution_summary.json`。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: Slide 8 要比較原始 T/Topo 與 VAF 後可優先選擇的結果；Slide 10 要有所有樣本整體圖表。
**DO NOT change**: 不得把 VAF-prioritized 寫成 VAF-corrected/confirmed true topology。
**Rationale**: 2026-07-11 使用者要求補充比較；現有證據只支援 candidate prioritization。
<!-- END USER-SPECIFIED -->

### 2026-07-11 16:48 — C 是 read-supported mutation-state group 數

- **Status**：Accepted（user-specified）
- **決定**：We will define `C_{region,HP}` as the number of full-span, MINREAD-supported, ALT-containing R/A genotype groups in `populations_by_hp[HP]`。
- **排除**：virtual `ROOT`、all-reference group、hidden Steiner node、partial-only subcube 均不算入 C。
- **雙 HP**：保留 `(C_HP1,C_HP2)`，並可報 `C_region=C_HP1+C_HP2`；不相乘。
- **影響**：舊報告把 candidate-tree joint combinations 稱為 C 的欄位全部改名為 `T_joint`。
- **Revisit if**：使用者明確要求把 all-reference observed group 納入 C；屆時另列 `C_all_full`，不可覆蓋主 C。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐（工程定義）；生物 clone 解讀不在本決定範圍。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: C 表示 reads 認為有多少群；它是布爾超立方體上被觀測支持的 mutation-state vertices。
**DO NOT change**: 不得再把候選樹組合數命名為 C。
**Rationale**: 2026-07-11 使用者修正。
<!-- END USER-SPECIFIED -->

### 2026-07-11 16:48 — T 表示 exact candidate tree structures

- **Status**：Accepted（user-specified）
- **決定**：We will use `T_unit=n_trees` for one primary HP and `T_joint=product(T_unit)` for a region with one or two primary HP trees。
- **Topology**：用 rooted unordered AHU canonical form 忽略 mutation label 與 sibling order；每個 shape 配穩定 `Topo_###` ID。
- **數量口徑**：主 shape quantity 為「允許此 shape 的 primary HP units 數」；另列「此 shape 為唯一可行形狀的 units 數」，避免多候選 unit 被重複加權造成誤讀。
- **Revisit if**：clean layered-v3 schema 直接輸出完整 canonical-shape histogram。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐（決定性同構計算）。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: 樹結構以 T 表示；報告要清楚統計有哪些 topology shapes 與數量。
**DO NOT change**: T 不得再與 C 混名。
**Rationale**: 2026-07-11 使用者修正。
<!-- END USER-SPECIFIED -->

### 2026-07-11 16:48 — 多 T 用逐位點 VAF 作相對排序

- **Status**：Accepted with scientific boundary
- **決定**：We will compute family-specific read-derived `VAF=ALT/(REF+ALT)` per locus and rank every complete non-capped candidate T by ancestor–descendant VAF ordering。
- **Winner rule**：default temperature 0.05、softmax relative weight ≥0.60、violation margin 0.05；未達門檻歸多種/擴散，不硬選。
- **命名限制**：只稱 `VAF-prioritized most-likely T` 或 relative ordering；不稱 calibrated posterior、CCF、independent confirmation 或 true clone tree。
- **Revisit if**：完成 purity/CN/multiplicity correction、depth uncertainty model 與 orthogonal truth validation。
- **Evidence tier**：L2 ⭐⭐⭐⭐（heuristic implementation）；真實演化確認為 L4 ⭐⭐。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: 報告必須說明多顆 T 的最可能者是利用每個位點 VAF 判斷。
**DO NOT change**: 不得省略 VAF ranking 的邏輯。
**Rationale**: 2026-07-11 使用者補充。
<!-- END USER-SPECIFIED -->

## 🟠 偏離之處

### 2026-07-11 16:48 — 完整樣本列使用 historical engineering snapshot

- **Status**：Accepted for this report only
- **規範要求**：Task B 應使用 clean raw-all → layered-v3 7/7 canonical output。
- **實作偏離**：active clean producer 尚未完成；本報告使用 historical layered-v2 7/7 raw JSON，並顯著標 `PARTIAL / scientific NO-GO`。
- **理由**：使用者要立即取得 HCC1395 範例與完整樣本列；現有 historical root 是唯一同時具有 7/7 C/T/Topo raw artifacts 的資料鏈。
- **回退路徑**：clean layered-v3 7/7 `_SUCCESS` 後以相同分析腳本重生，不沿用舊數字。
- **Evidence tier**：L2 ⭐⭐⭐⭐（工程描述）；正式 scientific claim 為 U。

### 2026-07-11 17:24 — HTML 只完成 structural-only QA

- **Status**：Accepted with visible limitation
- **實際結果**：portable artifact validation 與 package PASS；38 blocks、6 charts、20 tables、2 metric strips 的結構驗證 PASS。
- **偏離**：環境內 Chromium 在 static-chart extraction 階段提前退出，因此 source dialog、interaction 與 viewport QA 未執行。
- **處置**：使用 builder 內建 structural-only fallback；`qa_receipt.json` 顯式保存未驗證項目，不宣稱 browser QA PASS。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐（結構）；browser rendering 為 U。

## 🟡 折衷考量

### 2026-07-11 16:48 — shape quantity 採 unit-support，不採 candidate occurrence

- **Status**：Accepted
- **方案 A**：每 primary HP unit 對 shape 最多計一次，並另列 shape-only units。
- **方案 B**：每棵 candidate tree 都計一次；會讓高 `T` ambiguity unit 過度加權。
- **方案 C**：只列 coarse single/linear/branch；會丟失使用者要求的實際形狀 catalog。
- **採用 A 理由**：同時保存 exact canonical shape 與可解釋母數。
- **Revisit if**：未來需要 candidate-space entropy 報告，再加 occurrence-weighted appendix。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐（統計 grain 決定）。

## 🔴 未決問題

### 2026-07-11 16:48 — VAF missing units 的根因

- **Status**：open
- **Question**：各樣本 3–603 個 ambiguous primary units 為何缺完整 family-specific site VAF？
- **Default**：主報告列 `not_ranked_missing_VAF`，不插值、不歸類成多 winner。
- **Revisit if**：coordinate join、zero-depth 或 family coverage 稽核完成。
- **Priority**：major
- **Evidence tier**：L5 ⭐

## 📚 Lore

### 2026-07-11 — C=0 不等於沒有樹

- **Constraint**：primary unit 可只有 partial R/A/X groups，因而 `C=0`，但 solver 仍能由 subcube constraints 建 T。
- **Why it matters**：C=0 必標 partial-only evidence，不能解讀為 true negative。
- **Evidence**：historical mlhp `subread_groups_by_hp` + layered primary units。

### 2026-07-11 — VAF winner 與獨立確認不同

- **Constraint**：ranking score與 winner consistency均使用同一組 read-derived VAF，非 orthogonal validation。
- **Why it matters**：高一致率不能當第二條獨立證據。
- **Evidence**：`InterSubMod/docs/reports/validated/2026/07/20260710_分層重建端到端流程程式碼判斷稽核_01/20260710_分層重建端到端流程程式碼判斷稽核_01.md`。

### 2026-07-11 — `Topo_1/Topo_n` 與 `TopoShape-ID` 是不同欄位

- **Constraint**：前者描述一個候選集合內有 1 或 >1 個 shape；後者識別全資料中的具體 rooted shape。
- **Why it matters**：HCC1395 有 32 種 unit shapes，但單一 region 的候選集合仍可能是 Topo_1；雙 HP region 另有 311 種 ordered forest signatures。
- **Evidence**：`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json` 與 `regional_topology_composition.json`。

## ✅ 交付結果（2026-07-11 17:24）

- 7/7 dataset data validation：PASS；failures=0。
- 原始巢狀需求：HCC1395 71-row 與全樣本 457-row `HP scope × ordered C state × H × T/Topo` 長表，逐樣本 regions 皆守恆回 `W_primary`。
- HCC1395：10,355 primary HP units、53,467 exact candidate T、32 unit topology shapes。
- 全 7 rows：68,544 primary HP units、421,687 exact unit T、46 global rooted shapes。
- VAF 教學例：`chr4:40979546-40998095` HP2，T=2、Topo=1；site read-AF=0.9844/0.5098，top relative weight=0.999999994。
- HTML：`InterSubMod/docs/reports/in_progress/2026/07/20260711_HCC1395_read群C與Steiner樹T拓撲形狀全樣本數據報告_01/20260711_HCC1395_read群C與Steiner樹T拓撲形狀全樣本數據報告_01.html`。
- QA：artifact validation PASS、package PASS、structural-only PASS；詳見同目錄 `qa_receipt.json`。

## ✅ 教授版 Slide 8 / 10 更新（2026-07-11 18:52）

- Slide 8 保留 complete regions 的原始 `T/Topo` 三類圖，並於下方加入類別對應：
  - `T>1,Topo=1`：5,417/11,144（48.61%）可得 VAF-prioritized exact T；Topo 原本已唯一。
  - `T>1,Topo>1`：8,269/17,909（46.17%）可得 VAF-prioritized T 與其 Topo；不稱真實拓撲確認。
  - 合計 13,686/29,053（47.11%）。
- Slide 10 改為 7 dataset rows 的 region-level 100% 堆疊圖；每列四類狀態守恆回自身 ambiguous complete region 分母。
- 簡報在 `S` 來源抽屜新增第四張逐樣本 VAF-region 表，`D` 定義抽屜新增 `VAF-prioritized region`，Slide 8/10 speaker notes 保留完整規則與科學上限。
- 資料驗證：C/T `182/182 PASS`；VAF region `81/81 PASS`；semantic SHA-256=`fe02fcf689f25973e2d56a20b2ac1abb6725b6ec8202ab605bbc6f1639469509`。
- HTML 驗證：10 slides、6 charts/inline SVG、official chart contract PASS、0 remote/runtime dependencies。
- Browser QA：1920×1080、1366×768、390×844、320×568、no-JS desktop/mobile、10-page print PDF 全部 PASS。
- Fail-closed 資料 gate：C/T 必須 182 個 unique checks + frozen hash；VAF 必須 81 個 exact unique keys + frozen hash；summary 必須 7 rows + frozen hash + 逐列 JSON 對帳。
- Fail-closed 反事實：1-row checks、6-row summary、HTTPS image、protocol-relative image、CSS `url(//)`、JS `fetch("//")` 全部被拒絕。Browser QA 開始前另綁定 HTML/version/validation path、SHA-256 與 clean artifact PASS。
- 獨立終審：數字/分母、generator/browser、scientific narrative 皆 PASS，3/3 reviewers 無 blocker。
- 最終 HTML：`/big7_disk/liaoyoyo2001/Meeting/interSubMod_reports_workspace/20260711_分層重建教授重點簡報_01/slides/20260711_分層重建教授重點簡報_01.standalone.html`。
- 驗證收據：`/big7_disk/liaoyoyo2001/Meeting/interSubMod_reports_workspace/20260711_分層重建教授重點簡報_01/slides/20260711_分層重建教授重點簡報_01.validation.json` 與 `.browser_qa.json`。

## Provenance

- Commit：`6067568637088838a9f518955e41d222f057e4f1`
- Build time：2026-07-11 16:48 +0800
- Pre-decision：`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/pre-decision-audit.md`
- Target report：`InterSubMod/docs/reports/in_progress/2026/07/20260711_HCC1395_read群C與Steiner樹T拓撲形狀全樣本數據報告_01/`
- Region-level VAF 重算命令：

```bash
python InterSubMod/research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_region_resolution.py \
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 \
  --input-manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/input_manifest.json \
  --read-af-summary InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/read_af_tree_ordering_historical.json \
  --corrected-report InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json \
  --method-script-dir InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts \
  --output-json InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_region_resolution_summary.json \
  --output-summary-tsv InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_region_resolution_summary.tsv \
  --output-checks-tsv InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_region_resolution_checks.tsv
```
