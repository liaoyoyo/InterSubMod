<!--
建立時間: 2026-07-15 13:00
目標: current-v5 layered workstation 全基因組拓撲／形態與多選互動 living document
處理範圍: chr1-22、7 datasets、read-AF 描述性排序、desktop/mobile QA
cycle_id: cycle_20260715-1300-layered-workstation-genome-topology
spec_id: layered_workstation_genome_topology_multiselect
status: validated
advisory: on
關聯檔案:
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/pre-decision-audit.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
-->

# Implementation Notes：全基因組拓撲圖層與多選互動

> **Purpose**：即時記錄本次 current-v5 數據接合、語意界線、互動決定與 QA。
> **Task type**：B Comprehensive validation；全 chr1–22 × 7 datasets。

## 🔵 設計決定（Design Decisions）

### 2026-07-16 — 視覺誠信與 answer-first remediation addendum

- **Status**：Accepted / validated
- **背景**：全量 Chromium 稽核雖通過既有 219/219 checks，但新的人眼與 DOM 稽核發現 overview 長條以 panel 最大值當 100%，與文字所用 canonical denominator 不一致；morphology 的單支與無 primary 也共用同一 slate 實線。
- **決定**：We will derive every overview bar and printed percentage from the same panel denominator, and give N/A a dashed/hatched non-colour encoding distinct from the solid single-branch category.
- **理由**：前者是數據視覺誤述，後者會讓類別映射失去一對一關係；兩者皆列 P0。
- **影響範圍**：7 個 sample pages、renderer UI contract、Playwright assertions。
- **Revisit if**：無；這是固定資料誠信契約。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-16 — Index 與 sample page 改成 answer-first 閱讀順序

- **Status**：Accepted / validated
- **決定**：Index 先顯示 HCC1395 technical verdict 與三維跨樣本比較，再放方法教學；sample page 先顯示 canonical 摘要、GRCh38 ideogram 與 region browser，再放三維解釋與七組完整分布。兩層都加入可見 sticky local nav。
- **理由**：核心研究問題原本在 desktop 第 3–5 個 viewport、mobile 第 7–11 個 viewport，無法在第一輪掃讀完成判斷。
- **影響範圍**：資訊架構、responsive CSS、鍵盤導覽；既有 section ID 保留以避免 deep-link 失效。
- **Revisit if**：後續正式使用者測試顯示另一個任務序列更常見。
- **Evidence tier**：L2 ⭐⭐⭐⭐

### 2026-07-16 — 比較視圖的三組分母／聚合切換

- **Status**：Accepted / validated
- **決定**：We will expose 7-dataset operational versus 6-biological-ID macro profiles, full W_primary versus conditional-evaluable TVD, and HCC confusion count versus row percentage.
- **資料契約**：Full TVD 讀既有 7×7 matrix；conditional TVD 僅由 21 筆 pair records 的 `conditional_evaluable_tvd` 組裝；row-% 僅由 confusion row count / `left_counts[row]` 計算。
- **限制**：Conditional 是排除 incomplete／unavailable／unresolved 後的 renormalized subset，不能取代 full-profile verdict；HCC pair 是 technical comparison，不是 biological replication。
- **影響範圍**：index comparison payload、matrix/confusion controls、ARIA 與 tests。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-16 — Claude Code 雙階段唯讀共識 gate

- **Status**：Accepted / validated
- **初審 verdict**：`AGREE_WITH_CHANGES`。
- **共同 blocker**：sample renderer 必須改 `build_layered_workstation_v5.py` 並全量重建；overview 分母 bug 必修；USER-SPECIFIED Set 聯集／再點取消／零選全部不得改動。
- **共同驗收**：bar 分母、非色彩類別編碼、contrast、answer-first、sticky nav、matrix scale/header、macro 實圖、full/conditional、confusion count/%、caption/ARIA/claim/hidden JSON 與全量 Playwright 全數通過後，才可宣布 `AGREE`。
- **最終 verdict**：Claude Code 初審 `AGREE_WITH_CHANGES`；修正後審查為 `AGREE_WITH_NONBLOCKING_NOTES`；補上 index UI contract 並以位元層級排除 tab 誤判後，收斂複審為 `AGREE`，`BLOCKERS / MAJOR / MINOR = none`。
- **Evidence tier**：L2 ⭐⭐⭐⭐

### 2026-07-15 13:00 — current canonical v5 是唯一數據來源

- **Status**：Accepted
- **背景**：historical layered-v2 region keys 與 current-v5 只有部分重疊。
- **決定**：We will recompute read-AF ranking from the current canonical v5 run and never join the historical VAF TSV into current pages.
- **理由**：避免版本錯接；代價是需一次 durable exhaustive build。
- **影響範圍**：current ranking script、artifact、7 sample HTML。
- **Revisit if**：current run root 或 manifest 改版。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-15 13:00 — 把 structural、read-AF selection 與 morphology 分開

- **Status**：Accepted
- **背景**：`Topo=1`、exact tree unique、read-AF unique-first 與直系／旁系是不同維度。
- **決定**：We will expose separate ideogram layers for structural determinacy, descriptive read-AF selection status, and morphology/clone-compatibility.
- **理由**：避免把可辨識度、順位與生物 clone 混成一個 label。
- **影響範圍**：ideogram mode、legend、region detail badges、overview cards。
- **Revisit if**：有 calibrated CCF/purity/CN model 可支持更高 claim。
- **Evidence tier**：L2 ⭐⭐⭐⭐

### 2026-07-15 13:00 — 類別圖例使用聯集多選

- **Status**：Accepted
- **背景**：用戶要求可一次顯示多類；同類第二次點擊取消；沒有選取就是全部。
- **決定**：We will store selected category keys in a Set, render their union, toggle one key per click, and interpret an empty Set as all categories.
- **理由**：語意直接、可逆，也適合 keyboard/aria-pressed。
- **影響範圍**：所有 GRCh38 ideogram modes 與可點擊 summary bins。
- **Revisit if**：未來需要 AND/intersection filter（目前類別互斥，聯集才合理）。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

<!-- BEGIN USER-SPECIFIED -->
**Decision**：圖例類別可多選；單個已選類別再點一次取消；零選取表示全部。
**DO NOT change**：不得改回單選 radio 語意。
**Rationale**：2026-07-15 user request。
<!-- END USER-SPECIFIED -->

### 2026-07-15 13:00 — 無 read denominator 顯示 N/A

- **Status**：Accepted
- **背景**：目前 0/0 coverage 會被格式化成 `0 ALT/read ALT 0.0%`。
- **決定**：We will display `N/A · 此 evidence scope 無有效 read coverage` whenever REF+ALT=0 and exclude it from read-AF ordering. Aggregate site evidence 與 family-specific ranking input 會明確分開。
- **理由**：0/0 不是 0%，也不能支持候選排序。
- **影響範圍**：site evidence、ranking eligibility、tooltip。
- **Revisit if**：無。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

## 🟠 偏離之處（Deviations）

### 2026-07-15 13:00 — 「最可能」改成「read-AF 唯一第一順位」

- **Status**：Accepted
- **規範要求**：用戶希望顯示「有最可能性」。
- **實作偏離**：介面使用 `read-AF 唯一第一順位（描述性）`；說明中可對照「最可能性提示」，但不稱 probability／confirmed。
- **理由**：ordering score 未經 CN/purity/CCF 校準，不能解讀為機率。
- **風險評估**：降低誤讀；不影響排序資訊完整度。
- **回退路徑**：若未來完成 calibration，可新增獨立 calibrated likelihood layer。
- **Revisit if**：有 pre-registered calibrated model 與外部驗證。
- **Evidence tier**：L2 ⭐⭐⭐⭐

## 🟡 折衷考量（Trade-offs）

### 2026-07-15 13:00 — durable ranking artifact 必須包含 top edges

- **Status**：Accepted
- **方案 A**：We will rerun exhaustive candidates once and save exact-top candidate edges plus ranks for HTML consumption。
- **方案 B**：只重排 current JSON stored first 32 trees；rejected，HCC probe 有 241 units 的 top candidate index >32。
- **方案 C**：在每張 HTML build 時即時計算；rejected，重複 solver 7 次且 provenance 不清。
- **採用 A 理由**：正確、可重用、可稽核。
- **Tier check**：HCC current probe L1 engineering observation。
- **Revisit if**：upstream layered JSON 改為保存所有 candidates 或原生 ranking。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-15 13:00 — clone/subclone 只顯示相容形態

- **Status**：Accepted
- **方案 A**：We will label `單支／直系／旁系／直系+旁系相容` and retain an unresolved category。
- **方案 B**：標 confirmed clone/subclone；rejected，mutation-state trees 不是 clonal truth。
- **方案 C**：完全不提 clone；rejected，會失去使用者需要的解讀橋梁。
- **採用 A 理由**：保留可讀性且守住 scientific claim ceiling。
- **Revisit if**：有 purity/CN/multiplicity/CCF 與 orthogonal validation。
- **Evidence tier**：L2 ⭐⭐⭐⭐

## 🔴 未決問題（Open Questions）

### 2026-07-15 14:00 — read-AF first tree 的詳細視圖排序

- **Status**：resolved
- **Question**：top tree edges 是否可無損接回現有 per-family network viewer，並同時保留其他候選瀏覽？
- **Context**：現有 JSON display cap 32，current top candidate 可能在 cap 外。
- **Options**：A 重建 read-AF 排序後 top32；B 只插入 exact-top；C 顯示獨立 selected-tree card。
- **Decision**：採 C。每個 primary HP lane 先顯示獨立 read-AF ranking card（exact Fraction、preview、co-top representative edges 與原 candidate index），下方仍保留原 structural stored candidate browser。
- **Reason**：不重排或截斷 structural viewer；第一順位落在舊 top-32 以外仍能無損顯示，也能讓使用者分清「排序結果」與「完整／stored 候選瀏覽」。
- **Revisit if**：HTML contract 或 bundle size 超標。
- **Priority**：closed
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

## ✅ 實作與驗證紀錄

### 2026-07-15 13:49 — current-v5 exhaustive sidecar 完成

- **輸入**：
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/`
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json`
  - `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`
- **命令**：`python3 research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py --run-root … --input-manifest … --current-summary … --method-script-dir docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts --output-dir research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology`
- **輸出**：7 個 compact sample sidecar + `current_v5_read_af_topology.index.json`，共約 77 MB。
- **實際結果**：7/7 `all_checks_pass=true`；W_tree=51,815、W_primary=50,215；candidate／shape mismatch=0；所有 ranked unit 均有 exact-top edges。完整執行約 8:03，最大 RSS 約 665 MB。

| 互斥分布（W_primary） | 數量 | 比例 |
|---|---:|---:|
| 結構已 exact 唯一 | 11,582 | 23.1% |
| read-AF 唯一第一順位 | 20,732 | 41.3% |
| 並列第一、同一 Topo | 6,768 | 13.5% |
| 並列第一、多 Topo | 2,152 | 4.3% |
| read-AF 不可排序 | 1,006 | 2.0% |
| 候選未完整 | 7,975 | 15.9% |

| Mutation-state morphology（W_primary） | 數量 | 比例 |
|---|---:|---:|
| 單支／無 HP 內關係 | 6,369 | 12.7% |
| 直系鏈 | 28,855 | 57.5% |
| 旁系／分支 | 1,715 | 3.4% |
| 直系＋旁系 | 2,428 | 4.8% |
| 形態未解 | 10,848 | 21.6% |

兩組分布各自加總 50,215；它們是 7-dataset aggregate（HCC1395 與 DORADO 是同一 biological sample 的兩個 dataset），不是 biological-sample 加權。

### 2026-07-15 14:28 — 7 頁與 index 最終重建

- **輸入**：current summary、7 份 current-v5 region-view、7 份 read-AF/morphology sidecar、最新版 renderer。
- **命令**：`python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`
- **輸出**：`InterSubMod/docs/methodology/_assets/layered_workstation/{COLO829,H1437,H2009,HCC1395,HCC1395_DORADO,HCC1937,HCC1954,index}.html`
- **實際片段**：`[sample] BUILT hash-bound page` × 7；`OK canonical index + 7/7 hash-bound sample pages`；exit 0；elapsed 2:08；max RSS 808,740 KB。
- **Fail-closed readback**：`python3 …/build_layered_per_sample.py --index-only` 顯示 `[sample] VERIFIED hash-bound page` × 7，exit 0。
- **Provenance**：index 另綁 read-AF index SHA；driver 同時核對 ranking generator、exhaustive solver 與 input manifest SHA。

### 2026-07-15 14:30 — 手機觸控 preflight 修正

- **觀察**：390 px viewport 的 read-AF overview bins 為 40.6 px，高於文字可讀門檻但低於 44 px touch target。
- **修正**：共用 `.overview-bin-button` 加上 `min-height:44px`；mode、legend 與 chromosome controls 原已 ≥44 px。
- **驗證狀態**：正式 7 datasets × desktop/mobile Playwright receipt 進行中。

### 2026-07-15 15:00 — 第一順位 network 初始可見性修正

- **人眼與座標證據**：HCC1395 mobile fixture 的 network scroller `clientWidth=288`、`scrollWidth=736`、`scrollLeft=0`；5 個節點全在初始 viewport 右側，`visible_nodes=0`。雖可橫捲抵達，但不符合「打開即看見樹」。
- **根因**：network SVG 為保留標籤可讀性，最小寬度 720 px；單鏈節點置於 SVG 中央，overflow scroller 初始停在最左端。
- **修正**：保留 720 px 可讀畫布與局部橫捲，DOM 插入後把每個 `.network-scroll` 的初始 `scrollLeft` 設為可捲範圍中點；使用者後續仍可自由捲動，candidate 切換時只初始化新 scroller。
- **新 regression check**：每個 sample × desktop/mobile 的 read-AF top tree 都量測 node bounding box 與 scroller viewport 交集，要求每棵初始 `visible_nodes>=1`。
- **重建**：7/7 hash-bound sample pages，exit 0；elapsed 1:17；max RSS 808,960 KB。

### 2026-07-15 15:30 — auto-fit 與最終 Playwright 全量通過

- **版面收斂**：read-AF 第一順位只有一個 shape 時卡片自動撐滿；兩個以上才分欄。HCC1395 desktop 單樹由半寬改為全寬，`clientWidth=scrollWidth=750`、`visible_nodes=5`；mobile 維持局部捲動並自動置中，`clientWidth=288`、`scrollWidth=736`、`scrollLeft=224`、`visible_nodes=5`。
- **輸入**：最終 7 個 hash-bound standalone HTML + 7 個 current-v5 sidecar。
- **命令**：`python3 research/20260715_layered_workstation_genome_topology_multiselect/scripts/audit_layered_workstation_playwright.py --workstation-dir docs/methodology/_assets/layered_workstation --sidecar-dir research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology --output-dir research/20260715_layered_workstation_genome_topology_multiselect/qa/full --timeout-ms 180000`
- **輸出**：`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/qa/full/validation_receipt.json` + 18 張 desktop/mobile screenshots。
- **實際結果**：Chromium 147.0.7727.15；index 2 runs + samples 14 runs = 16/16，219/219 checks、console/page errors=0、exit 0；elapsed 2:08.30；max RSS 143,720 KB。
- **Receipt schema／SHA-256**：`1.1.0`／`7e915dec5224d062c89f64633cbbef8de58b9ac2efc5c944eef66e6b1e5fdfd3`。
- **Index 專項**：read-AF index meta SHA exact match；兩張 aggregate cards 各守恆到 50,215；desktop/mobile body overflow=0；17/17 `.json` links 全位於 `<details>`。
- **Quick cross-workstation regression**：`validate_workstation_ui.py --quick` 另通過 8/8 runs、108/108 checks，含 index desktop/mobile 與 archived topology route。
- **Claude Code 唯讀第二審**：縮小至關鍵函式後完成；多選空集合、claim ceiling、SHA fail-closed、初始置中一次性 guard 四項皆 PASS，無 blocker/high/medium。
- **演算法 contract tests**：`python3 research/20260715_layered_workstation_genome_topology_multiselect/scripts/test_current_v5_read_af_topology_contract.py -v`；5/5 PASS、exit 0，覆蓋 sibling-order invariant、四種 morphology、Fraction score、0/0 N/A 與跨 HP 禁止建 edge。

### 2026-07-17 00:14 — answer-first remediation、三視窗全量稽核與 Claude 共識完成

- **輸入**：
  - `InterSubMod/docs/methodology/_assets/layered_workstation/index.html`
  - `InterSubMod/docs/methodology/_assets/layered_workstation/{COLO829,H1437,H2009,HCC1395,HCC1395_DORADO,HCC1937,HCC1954}.html`
  - `InterSubMod/research/20260715_sample_topology_comparison/artifacts/sample_topology_comparison.json`
  - `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology/`
- **全站命令**：`python3 research/20260715_layered_workstation_genome_topology_multiselect/scripts/audit_layered_workstation_playwright.py --workstation-dir docs/methodology/_assets/layered_workstation --sidecar-dir research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology --output-dir research/20260715_layered_workstation_genome_topology_multiselect/qa/full --timeout-ms 180000`
- **全站輸出**：`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/qa/full/validation_receipt.json` + 27 張 desktop 1440／mobile 390／narrow 320 截圖。
- **全站結果**：Chromium 147.0.7727.15；8 documents、24/24 page runs、464/464 checks、0 console/page errors、exit 0；receipt SHA-256 `4b94ef59446ab20e1cbb8ec533e0c736582d6e1223c2b281cb1000bd25e42156`。
- **比較頁命令**：`python3 research/20260715_sample_topology_comparison/scripts/audit_sample_topology_comparison_playwright.py --index docs/methodology/_assets/layered_workstation/index.html --comparison research/20260715_sample_topology_comparison/artifacts/sample_topology_comparison.json --output-dir research/20260715_sample_topology_comparison/qa/playwright --timeout-ms 180000`
- **比較頁輸出**：`InterSubMod/research/20260715_sample_topology_comparison/qa/playwright/validation_receipt.json` + 15 張 desktop／mobile／narrow 截圖。
- **比較頁結果**：3/3 runs、66/66 checks、exit 0；receipt SHA-256 `1e39f500d5758b7778b71c6f4ebd52ed57f76f39e53772d8de357a4d83ad00b1`；index SHA-256 `fbd9ff22b0556ac37e1d207b7a8fb4a7df5c7b4eb13ba54610b1bf379811bca5` 與 receipt input binding 相同。
- **新固定契約**：overview bar 與文字百分比共用 canonical denominator；N/A 以斜紋／虛線加文字非色彩編碼；index 與 sample pages 共用 `layered-workstation-v5-grch38-topology-multiselect-3`；sticky nav 不遮標題；七種 ideogram view 保留 Set 聯集多選與零選全部。
- **資訊層次**：index 先回答 HCC1395 technical reproducibility 與跨樣本比較，再展開方法；sample page 先顯示摘要、GRCh38、region，再放維度教學與完整分布。HCC1395 × DORADO 僅支撐 `PARTIAL TECHNICAL REPRODUCIBILITY`，不支撐相同 ancestry／真樹／clone。
- **Claude Code 最終唯讀複審**：`VERDICT: AGREE`；`BLOCKERS: none`；`MAJOR: none`；`MINOR: none`。複審重新核對 index UI-contract、兩份 receipt、index hash、U+0009=0 與六檔 `git diff --check`，確認前一審 13 項結論仍成立。

### 2026-07-23 07:15 — 完整候選邊聯集、selected-tree overlay 與 7 頁重建

- **Task type／研究目標**：B（Comprehensive validation）；全基因組、7 datasets，不採 subset；服務 G4（多樣本一致性／可重現性）與 G5（可外部驗證的工程證據）。
- **根因**：HCC1395 `chr10:87818272-87928739` HP1 有 74 個完整 exact candidates，但舊 candidate browser 只儲存／展示前 32 個；read-AF 第一順位是原 candidate `#63`，因此 `H_RRAA → H_ARAA` 雖存在於最終選擇，舊 stored union 無法顯示。
- **Python sidecar 修正**：schema 升為 `1.1.0`；每個成功完整重枚舉的 T>1 unit 新增 compact `full_edge_census`：`candidate_total`、`top_candidate_total`、`node_total`、`edge_total` 與 `[parent, child, candidate_count, top_candidate_count]`。計數是 candidate membership，不是 posterior、機率或 read-support weight；不完整 unit 不宣稱完整 census。
- **全量命令**：`python3 research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --input-manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json --current-summary research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json --method-script-dir docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts --output-dir research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology`。
- **Sidecar 輸出／結果**：`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology/`；7/7 schema `1.1.0`、`all_checks_pass=true`；38,748 個完整 multi-candidate units、372,315 個 edge rows；index SHA-256 `b8b8f0e0de4035f58fba623bd9dee08ee0ea9c1ba2e7e03b58169d45cdffc6dc`，7 個 sample hash readback 全 PASS。
- **目標 fixture**：candidate `74/74`、top `1/1`、16 nodes、32 union edges；`H_RRAA → H_ARAA = 2/74` 且 `top = 1/1`。這證明該邊不是無來源輸出，而是舊 32/74 display cap 遺漏。
- **Renderer 修正**：UI contract 升為 `layered-workstation-v5-grch38-topology-multiselect-4`；完整 census 畫 exhaustive edge union，forced／variable 由 `n/N` 決定，read-AF co-top 以藍色 halo 疊加；另提供可摺疊逐邊表與明示 `Stored candidate preview 32/74`。union 文案明示「不是單棵樹」，selected 也不是 truth／probability。
- **演算法 regression**：`python3 -m unittest discover -s research/20260715_layered_workstation_genome_topology_multiselect/scripts -p 'test_current_v5_read_af_topology_contract.py' -v`；11/11 PASS，涵蓋完整／top incidence、denominator drift、selected-edge drift、missing read-AF、recurrence 與 incomplete absence。
- **跨面板 regression**：sidecar schema／row closure → full union 32 unique edges → target base edge／selected halo → read-AF top tree → stored preview 32/74；桌面、390 px、320 px 三視窗均驗證同一 fixture。
- **重建命令／輸出**：`python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`；輸出 `InterSubMod/docs/methodology/_assets/layered_workstation/{COLO829,H1437,H2009,HCC1395,HCC1395_DORADO,HCC1937,HCC1954,index}.html`。`--index-only` readback 顯示 7/7 `VERIFIED hash-bound page`；HCC1395 HTML SHA-256 `6f40e3ec850953c94e319775031560e8ff32fe40c819c8a3f6ae876419f9686c`。
- **Chromium 最終結果**：工作站 8 documents × 3 viewports = 24/24 runs、491/491 checks、27 screenshots、0 failures；receipt `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/qa/full/validation_receipt.json`，SHA-256 `1326b4d90062220a9cb674ea61cf255820e2be6715c49f321c8fed2ca6bfd6bd`。比較頁另為 3/3 runs、66/66 checks、15 screenshots。
- **視覺 remediation**：實際瀏覽器先抓到 locus `R`／`·` 對比 3.39／1.63 與 HCC1395 mobile SHA digest 撐寬；修正後 WCAG AA failure=0，1440／390／320 的 body overflow 都為 0。完整聯集人工截圖：`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/qa/full/HCC1395__desktop__edge_union.png`、`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/qa/full/HCC1395__mobile__edge_union.png`。

## 📚 Lore

### 2026-07-15 — Exact tie 與 softmax weight 不能混用

- **Constraint**：主分類以 Fraction score 完全相等判 exact tie；softmax weight 只作描述性 concentration。
- **Why it matters**：避免 arbitrary epsilon 或 posterior 語意混入 determinacy。
- **Evidence**：`research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_top_tie_census.py:269-294`。

### 2026-07-15 — Topo=1 不等於 clone 已確認

- **Constraint**：shape 是 rooted unordered mutation-state topology，忽略 mutation label 與 sibling order。
- **Why it matters**：只能說形態唯一／相容，不能說 biological lineage 已確認。
- **Evidence**：`research/20260712_vaf_selected_shape_four_class_census/20260712_VAF最終單一Shape四類全樣本重算_01.md`。

## Provenance Footer

- **Base commit before remediation**：`b7826318a1813878a9156bf63401e2b45fbf74eb`
- **Final implementation commit**：`68a2e3545e6bf5b82ce5e3413173a0f36dc00a27`（`fix(workstation): clarify topology comparison and responsive audit`）。
- **Build time**：2026-07-15 13:00:30 +0800
- **Skill**：`/implementation-notes` v0.1
- **Pre-decision**：`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/pre-decision-audit.md`
- **Cycle**：`InterSubMod/state/cycles/cycle_20260715-1300-layered-workstation-genome-topology/`
