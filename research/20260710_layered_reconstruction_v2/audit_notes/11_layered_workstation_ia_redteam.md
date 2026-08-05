<!--
建立時間: 2026-07-14T20:30:00+08:00
目標: 對 layered_workstation 入口、7 個 dataset 頁與 generators 做全頁 IA／network-topology display／responsive／accessibility 紅隊稽核，提出 generator-first redesign spec
處理範圍: Task Type B comprehensive validation；7 datasets / 6 biological samples / chr1-22；桌機 1440×1000 + 手機 390×844；服務 G3/G4/G5
關聯檔案:
  - InterSubMod/docs/methodology/_assets/layered_workstation/index.html
  - InterSubMod/docs/methodology/_assets/layered_workstation/README.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation.py
  - InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
限制: 唯讀稽核；未修改任何既有 HTML、generator、README 或 data
-->

# Layered workstation 全頁 IA／generator／network-topology display 紅隊稽核

用 SCQA + Assertion–Evidence：先判斷是否可當 current canonical 入口，再把每個問題綁到頁面 selector／generator line／最新 canonical SoT，最後給可直接落地與驗收的 index + per-sample redesign spec。

> **TL;DR：目前視覺殼與離線互動健康，但不能繼續以「CURRENT / CANONICAL」發布；共有 4 個 P0、9 個 P1、2 個 P2。** 最大風險不是排版，而是工作站仍由 7/10 historical sources 生成、把 distinct full-span genotype 的 `c` 寫成 clone 數、用 Read/VAF softmax 排序等機率候選，以及從 stored prefix 推導「forced 邊」。在這 4 個 P0 關閉前，建議標為 **STALE / NOT FOR CANONICAL CITATION**。（影響：高；信心：高）

## 1. 稽核範圍與判定

### 1.1 Task classification 與研究邊界

- **Task Type B — comprehensive validation**：入口頁、README、7 個 dataset HTML、兩個 generator 全納入，不用 subset。
- **G3**：檢查 read-level sSNV／HP／甲基輔助的證據邊界是否清楚。
- **G4**：檢查 7 datasets 的一致口徑、全 chr1–22 scope、技術 replicate 與 reproducibility provenance。
- **G5**：提供 selector/line、當次 Chromium 證據、可量測 acceptance gates。
- Thread D：**有關**，但 L3 methylation 只可顯示 pending／auxiliary，不能承重或排序。
- Thread B 撤回範圍：**無**；本稽核不重啟 methylation filter 假說。
- KDE-corrected：**不適用**，本任務不重算 methylation metric。
- VCF caller AF：**不使用**；但 UI 現行自行以 read/VAF 排樹是本次 P0。
- 長計算／C++／搬移／NO-GO gate：未執行 C++ 或研究重跑；只讀現有 canonical v5 與 HTML。

### 1.2 最終 verdict

| 面向 | 判定 | 說明 |
|---|---|---|
| Artifact/browser health | **PASS** | 8 頁 × 2 viewport 均載入；0 console error、0 page error、0 remote request；無 body 級水平溢出 |
| Static integrity | **PASS** | HTML parse、embedded JSON、inline JS `node --check`、local link、duplicate ID、兩個 generator `py_compile` 均通過 |
| Current canonical custody | **FAIL / P0** | HTML 合計仍是 historical `568,080` tree-input records / `W_tree=48,959`，不是 canonical v5 的 `582,820` / `51,815` |
| Network/topology semantics | **FAIL / P0** | `c=單 clone`、Read/VAF posterior、stored-prefix forced edge 都超過現行 claim ceiling或違反 formal spec |
| Genome-wide scope legibility | **FAIL / P1** | 資料確實含 chr1–22，但頁面沒有 scope funnel 或 22-chromosome overview，使用者只能由染色體下拉選單反推 |
| Responsive/accessibility | **PARTIAL** | reflow 與焦點基本可用；高-k 位點表在手機會被裁切，動態 detail／tree carousel 的 ARIA 不完整 |

**發布 gate：NO-SHIP AS CURRENT CANONICAL。** 可以保留作 2026-07-10 historical engineering artifact，但 badge、入口文案與 README 必須先降級；重生 current 版時必須 generator-first 關閉 P0。

## 2. 當次視覺證據與閱讀路線

以下截圖全部由本次 2026-07-14 本機 Chromium／Playwright 重新擷取；不是沿用舊 QA 圖。截圖左上角可見的 skip link 是刻意 focus，用來驗證鍵盤可達性。

### Step 1 — cohort entry（桌機）：視覺層級清楚，authority 錯置

![桌機入口與七個 dataset 全頁證據](11_layered_workstation_ia_redteam_assets/contact_desktop_all_pages_v2.png)

- **健康度：視覺 PASS／內容 authority FAIL。** Hero → interpretation boundary → L0–L3 → cohort table 的節奏清楚。
- `index.html:105-114` 同時顯示 `Current · canonical` 與 2026-07-10 historical backbone/snapshot，形成自相矛盾的首屏信任訊號。
- 桌機 cohort table 可一次比較 7 datasets；但欄位仍是舊 `all-determined / ambiguous / capped`，沒有 canonical `complete / C_region / Topo_region`。

### Step 2 — per-dataset investigation（桌機）：兩欄 route 可用，候選集合語意錯位

![HCC1395 桌機 region 選取與 topology detail](11_layered_workstation_ia_redteam_assets/desktop_02_hcc1395.png)

- **健康度：互動 PASS／topology claim FAIL。** 篩選 → exact region search → Enter 選取 → detail 更新 → carousel 切換成功。
- 左側 region list、右側 evidence detail 在 1440px 易理解；Observed / Inferred、forced / variable 的視覺冗餘值得保留。
- 同一 detail 同時宣稱「等機率最小樹」與「Read/VAF posterior 98% 最一致」，直接破壞候選集等價性。

### Step 3 — cohort entry（手機）：reflow 正常，但跨樣本比較退化成橫向尋欄

![手機入口頁](11_layered_workstation_ia_redteam_assets/mobile390_01_index.png)

- **健康度：可操作 PARTIAL。** Hero、claim cards、L0–L3 均單欄 reflow，root 390/390 無 overflow。
- cohort table 為 1,160px 寬的水平 scroll region；初始只看得到 sample / sSNV / regions，最重要的 complete/C/Topo 狀態要多次橫滑。
- scroll region 有 `role=region`、`aria-label`、`tabindex=0`，這個 accessibility pattern應保留。

### Step 4 — per-dataset investigation（手機）：選取後 focus 正確，長 detail 與高-k table 失控

![HCC1954 手機 capped region detail](11_layered_workstation_ia_redteam_assets/mobile390_08_hcc1954.png)

- **健康度：PARTIAL。** 7/7 頁在手機點選 region 後都把焦點移到 `#detail`；無 body 級 overflow。
- HCC1954 selected detail 高 4,542px；family/tree/raw trace 連續展開，缺少 sticky「返回清單／上一區／下一區」與 detail section index。
- `8 sSNV` 位點證據表超出 `.note` 約 241–252px，因 `.detail{overflow:hidden}` 被裁掉，沒有可操作的水平 scroll container。

完整 screenshot 與機械 receipt：

- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/11_layered_workstation_ia_redteam_assets/playwright_full_page_receipt.json`
- `InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/11_layered_workstation_ia_redteam_assets/static_redteam_receipt.json`

## 3. 已確認的優點（redesign 應保留）

1. **離線可攜與零外部依賴**：16/16 browser runs 的 remote request 都是 0。
2. **整體 responsive shell 健康**：1440 與 390 的 `root.scrollWidth == clientWidth`；入口與 sample header 都能 reflow。
3. **基本調查路線成立**：7/7 dataset 的 determinacy filter、chromosome filter、sort、exact region search、keyboard Enter 與 hash link 都可運作。
4. **全基因資料並未被裁成 subset**：每個 sample 的 chromosome selector 都是「全部 + chr1–22」共 23 個 options；要改的是 scope 呈現，不是刪資料。
5. **證據狀態有雙重編碼**：Observed / Inferred / Confounded / Unavailable 不只靠顏色；節點也用填色／空心虛線／`ᴴ` 區分。
6. **不確定性沒有被完全隱藏**：ambiguous、capped、recurrence、tree count、shape count、stored prefix 都有局部文案。
7. **L3 boundary 很醒目**：入口與 7/7 sample 都顯示 `PENDING / 不 rank 樹`；這個 contract 必須保留。
8. **鍵盤基礎可用**：入口與 sample skip link 可 focus；mobile selection 會把 focus 送到 detail。
9. **靜態完整性良好**：8/8 HTML 無 duplicate ID 或 broken local links；7/7 embedded JSON 與 inline JS 可解析。

## 4. P0 findings — 發布前必關

### P0-01 — 「CURRENT / CANONICAL」指向 historical sources，數字與 backbone 已失效

**問題。** 入口與 7 個 sample page 都宣稱 current/canonical，但 builder 固定讀 7/10 舊輸入；最新 canonical v5 明確要求 LongPhase-S recalibrated `FILTER=PASS`。

**證據。**

- Generator ownership：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py:28-48` 固定 `MSROOT=/big7.../multisample_subclone`，HCC1395 甚至指向另一個 20260618 pilot path。
- Artifact selector：`index.html .eyebrow`、`.provenance`（`index.html:105-114`）；sample `.status.current`、`.statusline`（各 sample `:77`）。
- Latest SoT：`InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md:12-14,18-36,112-114`。
- 機械總和：artifact=`568,080 / W_tree 48,959`；canonical=`582,820 / 51,815`。逐 dataset 全部有 drift：

| Dataset | Artifact tree input | Canonical v5 | Δ | Artifact W_tree | Canonical v5 | Δ |
|---|---:|---:|---:|---:|---:|---:|
| COLO829 | 38,196 | 39,458 | +1,262 | 7,774 | 8,007 | +233 |
| H1437 | 75,578 | 79,655 | +4,077 | 8,865 | 9,238 | +373 |
| H2009 | 157,405 | 161,595 | +4,190 | 9,717 | 9,674 | -43 |
| HCC1395 | 113,997 | 113,061 | -936 | 7,927 | 8,222 | +295 |
| HCC1395_DORADO | 112,387 | 113,145 | +758 | 7,958 | 8,385 | +427 |
| HCC1937 | 49,548 | 52,115 | +2,567 | 2,695 | 3,612 | +917 |
| HCC1954 | 20,969 | 23,791 | +2,822 | 4,023 | 4,677 | +654 |

**根因。** Builder 以手寫 path + mtime 當 authority；`CURRENT / CANONICAL` 是 template literal，不是通過 `_SUCCESS`、verifier、selected role 與 hash gate 後才生成的狀態。

**Generator-first fix。**

- `build_layered_per_sample.py` 改收 `--run-root` 或 signed manifest，從 `<run>/samples/<dataset>/layered_region_view_<dataset>.json` 發現 7 datasets。
- build 前 fail-closed：root `_SUCCESS`、`verification_summary.all_pass=true`、7/7 `selected_tree_role=longphase_recalibrated_pass_vcf`、summary/input hashes 都必須吻合。
- badge 只能由 release manifest 產生；historical root 自動輸出 `HISTORICAL / NOT CURRENT`，不得靠 mtime 推測。
- README 與首頁 SoT 指向 2026-07-14 current summary，不再把已被 override 的 20260706 數字表當唯一 current data source。

**Acceptance。** 入口精確顯示 `7 datasets / 6 biological samples / chr1–22 / canonical v5`，且 totals 必須為 `582,820 → 469,849 → 194,149 → 51,815 → 50,215 → 42,240 complete + 7,975 incomplete`；任一 gate 失敗時 build exit non-zero 且不產生 CURRENT badge。

### P0-02 — UI 以 Read/VAF softmax 排序「等機率」候選，違反 formal red line

**問題。** Header、README 與 region detail 都稱候選為「等機率最小樹」，但 `ccfBlock()` 又計算 softmax posterior、標 `≥60% 最一致`；這不是純 annotation，而是對不可辨識集合做 total-order。

**證據。**

- Generator：`build_layered_workstation.py:350-370`；DOM selector `#detail .lin .note` 內文字 `Read/VAF 排序` / `posterior(softmax·TEMP0.05)`。
- 視覺：`desktop_02_hcc1395.png` 中候選 #1 被標 98%。
- Formal SoT：`InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md:115-121` 明定「枚舉，非最佳化，非排序」。
- CN 邊界也不一致：程式只在 `gain/amp` 停用排序（`build_layered_workstation.py:363-370`），`loss/loh/unknown` 仍可排名；但 dashboard 又說只有 neutral 可把 read/VAF 當未校正 CCF 參考（`:196-206`）。

**根因。** 7/9 的 exploratory CCF heuristic 被直接混進 canonical renderer，沒有受 formal-method claim contract 管理。

**Generator-first fix。** 從 canonical UI 完全移除 `posterior`、winner badge 與 tree reorder。若需保留研究史，只能放進明確的 `EXPERIMENTAL / DOES NOT ALTER CANDIDATE SET` drawer，且預設關閉、不可改 carousel order、不可寫「最一致」。

**Acceptance。** 8 個 HTML source 與 rendered body 都不得出現 `Read/VAF 排序`、`posterior(`、`softmax`、`最一致`；candidate tabs 永遠以 deterministic digest/order 呈現。

### P0-03 — 「forced 邊」從展示用 stored prefix 計算，對 3,311 個現行頁面 non-capped lineages 可能是假穩定

**問題。** Generator 只對 `L.trees`（最多 stored prefix）取交集；即使 `L.n_trees > L.trees.length`，仍把交集邊稱為「在全部等機率樹一致」。這把「所有已顯示的樹」誤寫成「完整 candidate set」。

**證據。**

- Generator：`build_layered_workstation.py:421-432` 先以 `L.trees` 求 `stable`，再於 `n_trees>ns` 只補一句「顯示前 N」。
- Caveat 只警告 capped：`build_layered_workstation.py:385`，沒有處理 non-capped 但 display prefix 不完整。
- 現行 HTML payload 有 **3,311** 個 non-capped lineage `n_trees > stored trees`；若直接改指 canonical v5，為 **3,468** 個。詳 `static_redteam_receipt.json`。
- Canonical v5 已提供 `analysis_candidate_set_complete`、`display_trees_complete`、`n_trees_stored`、`n_distinct_shapes_exact`、`verification_status`，但兩個 builders 的 term occurrence 都是 0。

**根因。** Analysis completeness 與 display completeness 沒有分層；renderer 自行從 display payload 重算科學狀態。

**Generator-first fix。**

- `forced_edges_exact` 必須在 upstream 對完整 candidate set 計算並附 digest；renderer 只顯示，不重算。
- 若 `analysis_candidate_set_complete=false` 或 `display_trees_complete=false`，forced badge 固定顯示 `NOT EVALUATED FROM PARTIAL DISPLAY`。
- 顯示 `stored / exact total`，例如 `32 shown / 125 exact candidates`；capped 顯示 `prefix only / lower bound`，不得顯示唯一或 forced。

**Acceptance。** 對所有 `n_trees_stored < n_trees` 的 unit，DOM 中 0 個 `forced` claim；candidate-complete + display-complete 的 unit 才可顯示 upstream verified forced edge，並以 digest regression test 對帳。

### P0-04 — `c` 被寫成 clone 數，與 canonical `C_region` 衝突且對 partial-only evidence 做錯推論

**問題。** `regionC()` 計數「有 full-span observed population 的 distinct ALT genotype」，UI 卻寫 `c=1·單 clone`、`c≥2 subclonal 候選`。現行 canonical 的大寫 `C_region` 是各 primary units `product(n_trees)`，兩者不是同一物件；更不能由 `0 full-span` 推成沒有 clone。

**證據。**

- Generator：`build_layered_workstation.py:183-206,347-348,377`；DOM selector `#detail .kv .b[title*="distinct"]`。
- Canonical definition：`InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md:55-62`。
- Canonical evidence：72,994 primary units 中 **44,672 是 partial-only**，但仍由 overlapping reads 約束候選集（同檔 `:38-53`）；`static_redteam_receipt.json` 獨立重算同為 44,672。
- Node 文案也把 mutation state 說成 clone：`build_layered_workstation.py:159,381` 的「未觀測 clone／中間 clone」。

**根因。** 舊 topology「clone 數 c」概念殘留，沒有隨 7/14 candidate-combination ontology 升級。

**Generator-first fix。**

- 全面停止產生 `clone 數 c`、`單 clone`、`中間 clone`。
- 使用兩個明確欄位：`C_region = exact candidate combinations`、`Topo_region = distinct topology combinations`。
- observed full-span genotype 另名 `n_full_span_observed_states`，只做 evidence coverage；partial-only 明示 `overlap-constrained`。
- node 一律叫 **mutation state**；hidden node 叫 **inferred intermediate state**。Edge 是單一 sSNV acquisition constraint，不是真細胞祖先證明。

**Acceptance。** UI 全文 0 個 `單 clone` / `clone 數 c` / `中間 clone`；每個 complete region 都滿足並顯示 `C_region ≥ Topo_region ≥ 1`，`C=1/Topo>1` 永遠為 0；partial-only 仍可正常呈現 candidate set。

## 5. P1 findings — 下一版核心 IA 與可及性

### P1-01 — 全 chr1–22 資料存在，但沒有可辨識的 scope funnel 或 genome overview

- **Evidence**：sample `#fchr` 有 23 options；但 `index.html:110-114` provenance 沒有 `chr1–22`，sample `.statusline` 也沒有；首頁直接從 sSNV 跳到 eligible regions。
- **Risk**：使用者無法分辨 all contigs、chr1–22、retained sSNV、W_tree、W_primary 的 scope，容易把 region browser 當全位點 census。
- **Fix**：入口與 sample 都增加固定 scope ribbon + disposition funnel；另做 22-row chromosome navigator，所有 chr1–22 永遠保留，filter 只改 active view、不改 scope truth。
- **Acceptance**：任一頁在首屏可讀到 `chr1–22 biallelic sSNV`；22 個 chromosome row/chip 都存在並有 count，總和回到 canonical aggregate。

### P1-02 — `region_determinacy` 單一 precedence label 混合了 candidate completeness、recurrence 與 topology identifiability

- **Evidence**：filter 只做 `r.region_determinacy===fd`（`build_layered_workstation.py:264-277`）。最新 canonical H2009 有 4 regions 同時 recurrence + capped，legacy label 顯示 recurrence，但 candidate-level 必須歸 incomplete（current report `:83-95`）。
- **Risk**：這 4 區在 capped filter 中消失；使用者可能把 recurrence badge誤讀成完整 candidate set。
- **Fix**：拆成正交 facets：`candidate completeness`、`exact cardinality C`、`topology cardinality Topo`、`recurrence present`、`verification`、`hidden-state`。
- **Acceptance**：H2009 四個 overlap regions 同時可由 `incomplete` 與 `recurrence` filter 找到，且永遠不進 complete topology denominator。

### P1-03 — HP primary／auxiliary／reference control 混在同一「germline-HP family」語言

- **Evidence**：`build_layered_workstation.py:377` 用 `r.lineages.length` 顯示「N germline-HP 家族」；實際 lineages 可含 HP3、none、reference-only。HCC1954 screenshot 中 HP1+HP2+none 被寫成 3 個 germline-HP 家族。
- **Risk**：HP multiplicity、primary topology denominator 與 auxiliary evidence 被混為一談。
- **Fix**：主畫布只放 primary HP1/HP2 lanes；H3/H4、unphased none、reference-only control 放到 auxiliary drawer，並顯示 `unit_role` / `is_primary_lineage`。
- **Acceptance**：family count 只等於 `n_primary_lineages`；H3/none 不增加 HP multiplicity；7 dataset totals 對回 72,994 primary units。

### P1-04 — PS 完全缺席；CN 有 badge 但缺 source/availability contract；L3 boundary則應保留

- **Evidence**：兩個 generator 對 `phase_set` / `mixed_PS` 都是 0 occurrences；最新 canonical 有 5,623 mixed-PS regions，PS 明定只作 QC、不是 topology edge（current report `:40-53`）。CN 只顯示 `neutral/gain/...`，沒有 source、misfit/unavailable；README `:18-20,38` 已要求 unavailable 不可當 neutral。
- **Risk**：使用者看不到跨 PS 的 QC 風險，且可能把 CN badge 或 HP family 當成 tree edge authority。
- **Fix**：region header 增加三個獨立 sidecars：`PS QC`、`CN post-hoc`、`L3 methyl auxiliary`；三者都放在 tree candidate set 外框，視覺上不得連成 edge。
- **Acceptance**：mixed-PS aggregate=5,623；PS 只顯示 counts/warning，0 個 PS-derived edge；CN 每區都顯示 source + available/unavailable/misfit；L3 繼續 `PENDING / DOES NOT RANK`。

### P1-05 — Region browser 能找 exact ID，但不足以支援 genome-scale investigation 與可分享查詢

- **Evidence**：`#fq` 只做 `r.region.toLowerCase().indexOf(fq)`（`build_layered_workstation.py:264-275`）；初始只 render 前 400；URL 只保存 `#region=`（`:440-456`），不保存 filters/sort。深連結呼叫 `show(i,null,false)`，不會同步 list selection。
- **Risk**：無法用 interval overlap、gene、candidate completeness、C/Topo、PS、CN、partial-only 找區；分享後 reviewer 無法重現查詢狀態。
- **Fix**：支援 `chr:start-end` overlap parser、region ID、可選 gene annotation；filter facets 全寫入 URL query；結果表可 pagination/virtualize；deep link 要同步 active row、breadcrumb 與 genome navigator。
- **Acceptance**：refresh/back/forward 後 dataset、run、chrom、interval、all filters、sort、region 全復原；查 `chr8:79992384-80149786` 可重現 H2009 recurrence+capped case。

### P1-06 — 手機 detail 太長且高-k tables 被裁切；入口 cohort table 不適合手機比較

- **Evidence**：`.detail{overflow:hidden}`（sample `:62`）；位點 evidence table 直接注入（generator `:397-407`），只有 `#dash table` 在 mobile 有 overflow rule（CSS `:67`）。HCC1954 8-site detail 截圖高 4,542px，內部 overflow 241–252px。
- **Risk**：S6–S8 等證據不可見／不可操作；使用者選完 region 後很難回清單。
- **Fix**：所有 wide table 包進具名 scroll region或改成 stacked site cards；mobile cohort 改 dataset cards；detail 增 sticky region bar（Back to results / Prev / Next）與 section accordions。
- **Acceptance**：390/320 下每個 overflowing element都有可操作 container、`role=region`、accessible name、keyboard focus；任何 site/column 都能到達；selected detail 首屏可見 region verdict + candidate summary。

### P1-07 — Region list、detail live region 與 tree carousel 的 ARIA model 不成立

- **Evidence**：`#list[role=list]` 的 child 是無 `role=listitem` 的 buttons，卻設 `aria-selected`（generator `:121,272-277`）；`#detail[aria-live=polite]` 每次替換 1,186–3,379 characters；tree SVG 無 `role/title/desc`，thumb 無 `aria-pressed/current`，counter 無 live state（`:215-257,419-435`）。Playwright receipt 7/7 都是 `svg role=null`, `activeThumbPressed=null`, `treeCountersLive=null`。
- **Risk**：assistive tech 可能一次朗讀整份 detail，卻無法知道 active candidate tree 或圖的整體含意。
- **Fix**：region results 用 native buttons + `aria-current`，或正式 listbox/option；detail 只以短 status announce「已選 X、complete/incomplete、C/Topo」，正文移除 aria-live；SVG 加 `role=img`、unique title/desc；tabs/thumbnail 使用 tablist 或 `aria-pressed`。
- **Acceptance**：Accessibility tree 中 active region、active topology shape、candidate index 都有唯一名稱與 state；切樹只 announce 短 counter；無 invalid ARIA。

### P1-08 — 7 datasets 被寫成 7 samples，技術 replicate 未分組

- **Evidence**：index instrument strip 寫 `Sample workstations 7/7`（`index.html:125-129`）；current SoT scope 是 `7 datasets / 6 biological samples`；HCC1395 與 HCC1395_DORADO 只靠小字說 same cell line，且表中未相鄰分組。
- **Risk**：cross-dataset inventory 被誤讀成 7 個獨立 biological samples，形成 pseudoreplication 語意。
- **Fix**：主詞改 Dataset；增加 biological sample group；HCC1395 兩個 dataset 放在同一 group，comparison 標 `technical/basecaller sensitivity`。
- **Acceptance**：首頁明示 7/6；任何 aggregate legend 都分 dataset count 與 biological sample count。

### P1-09 — 單檔可攜性達成，但 138 MB monolith 與 400-row nested scroll 讓 genome-scale UI 成本過高

- **Evidence**：8 HTML 合計 143,834,193 bytes；最大 H2009 38,274,656 bytes；本機 load 0.8–3.0s。每頁初始 render 400 buttons，list scrollHeight 約 21.7–22.0k px。
- **Risk**：低階手機解析 20–38 MB inline JSON、重建 400 buttons 與長 nested scroll，會放大 memory/interaction latency；目前測試未做 throttling，真實設備風險更高。
- **Fix**：保留單 HTML offline contract，但拆 embedded payload：lightweight region index + 22 個 chromosome detail chunks；先載 scope/genome/results index，選 chr/region 才解碼 detail。結果列表用 pagination/virtualization。
- **Acceptance**：repo machine H2009 initial interactive ≤2s；initial rendered region rows ≤100；切換 filter ≤100ms p95；仍可在斷網狀態檢索全部 chr1–22。

## 6. P2 findings — 語意與 polish

### P2-01 — Table semantics 不完整

- `index.html table` 無 `<caption>`，11/11 `<th>` 都無 `scope=col`；動態 dashboard tables 也沒有 caption。
- Fix：每表有 visible/visually-hidden caption、header scope；mobile scroll cue 由 `aria-describedby` 綁定。
- Acceptance：screen-reader table navigation 可讀表名、row/column headers 與分母。

### P2-02 — Share/copy control 與 touch targets 可再收斂

- `#copy-link` 在未選 region 前仍寫「複製目前 Region 連結」；使用 deprecated `document.execCommand('copy')`（generator `:454`）。sample select/filter 高 31–36px，雖高於 WCAG 2.5.8 的 24px AA 門檻，但不利手機觸控。
- Fix：未選時 disable；選取後使用 Clipboard API + fallback；主要 mobile controls target ≥44px，並 announce copy success。
- Acceptance：無 region 時不可 copy 假連結；copy URL 含完整 query/hash；觸控主控制達 44×44px。

## 7. 建議的 information architecture 與 section order

### 7.1 Index page — cohort / authority first

1. **Authority ribbon**：`CANONICAL v5`、run ID/hash、7 datasets / 6 biological samples、chr1–22、generated time、claim ceiling。
2. **Question / answer / limit**：保留現有三張 claim cards，但縮成首屏一列；用語改成 regional mutation-state candidate tree。
3. **Canonical scope funnel**：582,820 tree-input → 469,849 autosomal biallelic → 194,149 retained → 51,815 W_tree → 50,215 W_primary → 42,240 complete + 7,975 incomplete。
4. **C/Topo state matrix**：三個合法狀態 11,582 / 10,737 / 19,921 + incomplete 7,975；impossible=0 只作 invariant badge。
5. **Whole-genome navigator**：固定顯示 chr1–22，每 chr 顯示 W_primary、complete/incomplete、mixed-PS；可點進 dataset/chr，但不隱藏其他 chromosome 的 scope。
6. **Dataset matrix**：以 biological sample 分 group；HCC1395 內含 5kHz 與 DORADO datasets。桌機 table，手機 cards。
7. **Evidence boundary**：Primary HP1/2 → read-state constraints → candidate set；CN/PS/L3 是 sidecar，不畫進主 edge。
8. **Canonical vs sensitivity**：LongPhase-S PASS 主結果 vs ClairS PASS sensitivity，顯示 backbone-sensitive verdict與 Jaccard/digest；預設不混在主表。
9. **Methods / provenance drawer**：verification receipts、hash、README、method specs、download/portable notes。

### 7.2 Per-sample page — genome → region → candidate set → evidence

1. **Sticky dataset context**：dataset、biological group、canonical run、scope、selected region、Back/Prev/Next/Share。
2. **Sample funnel + status**：tree input、autosomal、retained、W_tree、W_primary、complete/incomplete、V1–V7。
3. **22-chromosome genome overview**：全 chr1–22 miniature rows；active chr 有清楚 highlight；all-genome totals 永遠可見。
4. **Search/filter bar**：interval/gene/region、chrom、candidate complete、C/Topo class、recurrence、primary HP multiplicity、full+partial/partial-only、PS、CN、hidden state、verification。
5. **Region results**：virtualized/paged table；每列顯示 region、k/span、primary HP、evidence coverage、complete、C、Topo、hidden、recurrence、PS、CN、verification。
6. **Region verdict header**：先用一句 assertion 回答「candidate set complete 嗎？Exact 幾個？Topology 幾種？主要限制是什麼？」
7. **Observed evidence**：位點表、full-span states、partial overlaps、read-span gaps；PS 只在此作 QC。
8. **Candidate-set explorer**：先按 topology shape 分組，再展開 exact state trees；顯示 `shown/total` 與 digest。
9. **Mutation-state tree**：nodes/edges/legend/edge acquisition；只顯示 upstream verified forced/variable status。
10. **Independent sidecars**：CN post-hoc、PS QC、L3 methylation pending，各自 source/status，不改 candidate order。
11. **Raw evidence drawers**：observed read groups、trace、V1–V7 receipt、source paths；預設折疊。

## 8. Network/topology display component contract

| Component | 必須回答 | 必備欄位／狀態 | 禁止語意 |
|---|---|---|---|
| Scope badge | 這是哪個 run、哪些 chromosome、哪個 denominator？ | run ID/hash、selected role、chr1–22、dataset/bio-sample counts | 僅靠 `CURRENT` literal |
| Evidence coverage | 候選由 full-span 還是 overlapping partial reads約束？ | full+partial / partial-only、n_full_pops、n_partial、read-span gap | `0 full-span = 無 clone` |
| HP lanes | 哪些是 primary reconstruction units？ | HP1、HP2 primary；H3/H4/none/reference auxiliary | 用 `r.lineages.length` 當 germline count |
| Candidate-set header | 集合完整嗎？exact / shape 各多少？ | complete/incomplete、C、Topo、stored/total、verification | capped prefix 宣稱唯一 |
| State node | 這個 binary mutation state是 observed 或 inferred？ | state vector、Observed/Hidden、support count | 稱 cellular clone |
| Edge | 哪個 sSNV acquisition 連接兩個 states？ | `+S_i`、forced/variable/not-evaluated、support scope | 稱真實祖先關係 |
| Shape tabs | exact trees 如何歸成 topology shapes？ | topology digest、exact count per shape、active state | 平鋪 125 trees 無分組 |
| PS QC | 是否混 phase blocks？ | PS counts、mixed flag、QC warning | PS 作 topology edge |
| CN sidecar | recurrence/multiplicity 判讀是否可用？ | source、available/misfit/unavailable、class | missing 當 neutral；用 CN 排候選 |
| L3 sidecar | methylation 是否參與？ | PENDING/AVAILABLE、`DOES NOT RANK` | 以 methylation 選樹 |

### Tree visual grammar

- Root：`0^k / all-reference mutation state`，灰色方形。
- Observed state：實心節點 + readable state label + read support。
- Hidden state：空心虛線節點 + `H` + 明文「inferred state；not directly observed」。
- Edge：標 `+S1` / `+chr:pos ALT`；solid=verified across complete exact set，dashed=varies，dotted gray=not evaluated because incomplete display。
- Exact candidate 與 topology shape 分開：`C` 不能拿來代替 `Topo`；同 shape 的 exact candidates 收在同一 tab。
- 不用顏色單獨承載狀態；line style、icon/label、文字三者至少兩種同時存在。

## 9. 保留／改寫／停止生成

### 保留

- 單檔離線、零外部依賴與 per-dataset portability。
- index → sample、sample switch、region deep link。
- 桌機兩欄、手機單欄基本 layout。
- skip link、focus outline、filter count `aria-live`。
- Observed / Inferred / Confounded / Unavailable vocabulary。
- HP1/HP2 color identity、hidden node dashed outline、forced/variable edge line-style redundancy。
- V1–V7、exact tree count、distinct shape count、raw read drawers。
- 「ambiguous/capped 是結果」與 L3 `PENDING / 不 rank`。

### 改寫

- `Sample` → `Dataset`，另加 biological sample group。
- `somatic sSNV` → `operational tree-input records` 或明確 LongPhase-S selected records。
- `region_determinacy` → independent candidate completeness + C/Topo + recurrence facets。
- `germline-HP 家族數` → `primary HP lineages`，auxiliary 分 drawer。
- 入口 wide table → desktop table + mobile cards；sample nested 400-row list → virtualized result set。
- source filename → run-root-relative path + SHA-256 + verifier status。

### 停止生成（不是刪除歷史檔；是 current renderer 不再輸出）

- Hard-coded `CURRENT / CANONICAL`。
- `clone 數 c`、`c=1·單 clone`、`中間 clone`。
- `Read/VAF posterior/softmax/winner`。
- 從 display prefix 計算的 forced edges。
- 用單一 precedence label 取代 completeness/recurrence 多標籤。
- 把 H3/none/reference controls 算入 germline family count。
- 在整個 `#detail` 上使用 `aria-live=polite`。

## 10. 全頁 acceptance criteria（Task B，不可 subset）

### A. Scientific/data custody

- [ ] 7/7 pages 來源皆是 canonical v5 root，且 `_SUCCESS`、7/7 verifier、selected role、hash 全 PASS。
- [ ] Index totals 精確等於 `582,820 / 469,849 / 194,149 / 51,815 / 50,215 / 42,240 / 7,975`。
- [ ] Topology classes 精確等於 `11,582 / 10,737 / 19,921 / incomplete 7,975 / impossible 0`。
- [ ] 7 datasets / 6 biological samples / chr1–22 在 8/8 首屏可見。
- [ ] 44,672 partial-only primary units 保留並明示 overlap-constrained；不能被歸為 no-clone/no-data。
- [ ] H2009 4 recurrence+capped regions 同時可由兩個 facets 找到，且計入 incomplete。
- [ ] 5,623 mixed-PS regions 可查；PS 不產生任何 edge。

### B. Candidate/topology semantics

- [ ] 0 個 `Read/VAF 排序`、`posterior`、`softmax`、`單 clone`、`中間 clone`。
- [ ] Candidate header 同時顯示 completeness、C、Topo、stored/total、verification。
- [ ] `analysis_candidate_set_complete=false` 或 `display_trees_complete=false` 時，0 個 forced-edge claim。
- [ ] Complete regions 全滿足 `C≥Topo≥1`；impossible state=0。
- [ ] Nodes、edges、HP lanes、PS/CN/L3 sidecars 的角色符合 §8 contract。

### C. Browser and responsive — 8 pages × 1440 / 390，另加 320 smoke

- [ ] `documentElement.scrollWidth == clientWidth`。
- [ ] 無 clipped site/table/tree；每個 overflow region 都可 keyboard focus、有 role/name/cue。
- [ ] Mobile selected region 首屏先看 verdict + C/Topo，而非直接落入 raw evidence。
- [ ] Back to results / Prev / Next 在長 detail 中保持可達。
- [ ] 所有 chr1–22 可由 genome navigator 到達；filter 不靜默改 denominator。
- [ ] Console error=0、page error=0、remote request=0。

### D. Accessibility

- [ ] 1 個 H1；heading order 可導覽 region、HP lane、candidate shape。
- [ ] Table 有 caption、`th scope`；scroll container 有 `aria-describedby`。
- [ ] Region selection 使用合法 state model；detail 不整段 live announce。
- [ ] Tree SVG 有 unique accessible name/description；active candidate/shape state 可由 AT 讀出。
- [ ] Carousel 可用 Tab + Enter/Space + arrow keys；active tab/thumbnail 有 programmatic state。
- [ ] Mobile primary controls target ≥44×44px；visible focus 不被 sticky header 遮住。

### E. Performance / portability

- [ ] 斷網仍可查全部 22 chromosomes 與每個 region detail。
- [ ] H2009 initial interactive ≤2s（同一 repo machine、warm server、3 runs median）。
- [ ] 初始 rendered region rows ≤100；filter interaction p95 ≤100ms。
- [ ] 每個 embedded chromosome chunk 有 hash；解碼失敗顯示明確 error state，不顯示空白 detail。

## 11. 建議實作順序

1. **先關 P0 custody**：manifest/run-root/fail-closed + canonical v5 rebuild；在此之前先把現行 badge 降為 historical。
2. **移除 P0 semantics**：刪除 renderer 內 CCF ranking、clone-c、prefix forced edge；改讀 upstream completeness/digest fields。
3. **建立 scope/C/Topo data contract**：index 和 sample 共用同一個 machine summary adapter。
4. **重建 IA**：scope funnel → genome navigator → region facets → candidate-set explorer。
5. **mobile/accessibility pass**：wide evidence containers、sticky return route、ARIA tree/tab model。
6. **全頁 Chromium gate**：8 pages × 1440/390/320 + keyboard + no-network + hash receipts。

## 12. 驗證命令、輸入、輸出與實際結果

### 12.1 Browser capture

**輸入**：`InterSubMod/docs/methodology/_assets/layered_workstation/{index.html,7 sample HTML}`。

**命令**：

```bash
python3 /bip7_disk/liaoyoyo2001/.codex/skills/webapp-testing/scripts/with_server.py \
  --server "python3 -m http.server 8765 --bind 127.0.0.1" --port 8765 -- \
  python3 /tmp/capture_layered_workstation_redteam.py
```

**輸出**：`InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/11_layered_workstation_ia_redteam_assets/`。

**實際輸出摘要**：16/16 screenshot；load 0.522–3.029s；root overflow 0/16；console/page/remote errors 0/0/0；7/7 exact search 回 `1 區`；7/7 mobile selection focus 到 detail。

### 12.2 Canonical comparison / generator field audit

**輸入**：現行 7 HTML embedded payload、canonical v5 7 個 region views、`current_layered_topology_v3_raw_all_v1.json`、兩個 generators。

**命令**：

```bash
PYTHONDONTWRITEBYTECODE=1 \
  python3 /tmp/collect_layered_workstation_redteam_static.py
```

**輸出**：`InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/11_layered_workstation_ia_redteam_assets/static_redteam_receipt.json`。

**實際輸出摘要**：artifact/canonical tree input=`568080/582820`；W_tree=`48959/51815`；artifact/canonical non-capped truncated lineages=`3311/3468`；canonical partial-only primary units=`44672`；generator occurrences for PS/candidate-completeness/exact-shape/verification fields均 0。

### 12.3 Syntax / local dependency checks

**輸入**：8 HTML + 2 Python generators。

**命令類型**：Python `HTMLParser`、embedded JSON `json.loads`、inline JS `node --check -`、local href existence、`python3 -m py_compile`。

**實際輸出摘要**：HTML parse 8/8 PASS；JSON 7/7 PASS；JS 7/7 PASS；broken links 0；duplicate IDs 0；generator compile 2/2 PASS。Index semantic follow-up：table caption 0、`th scope` 0/11。

## 13. Evidence limits

- 沒有用 NVDA／VoiceOver／TalkBack 實機朗讀；accessibility 結論分成「DOM/keyboard 已確認」與「需 screen-reader regression」兩級，未宣稱 WCAG 全面合規。
- 未做 CPU/memory throttling；load time 只代表本機 Chromium，本報告不把它外推成真實手機效能。
- 因任務是唯讀，沒有把 canonical v5 真正重生為新 HTML；redesign acceptance 的 current numbers 來自已驗證 canonical JSON/報告，不是假裝新 UI 已完成。
- 截圖用 localhost 服務靜態檔；另以 DOM 確認 0 外部 runtime dependency，與 offline artifact contract一致。
- 本次沒有修改 `index.html`、7 sample HTML、README、兩個 generator 或任何 data；新增內容只有本 audit 文件、screenshots 與 receipts。
