<!--
建立時間: 2026-07-23
目標: 稽核三條可能混用的方法版本、k>12/區域/區塊語意，以及 A-B、B-C read-link chain 的可主張範圍
處理範圍: 只讀方法與實作證據；不指定單一 SoT；不修改程式
服務研究目標: G1、G5
關聯檔案:
  - InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/20260715_GRCh38拓撲形態多選工作站完整檢視_01.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/00_INDEX.md
  - InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/pre-decision-audit.md
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md
-->

# 方法版本、k>12 與 A–B–C 串聯有效性稽核

用「版本矩陣＋claim conflict audit」：

> **TL;DR：A–B 與 B–C 各有達門檻的 direct endpoint evidence 時，把 A、B、C 放進同一個 read-linked connected component，在圖論與目前 strict-region 定義下是合理的；但它只建立「同一分析區域」的傳遞連通性，不建立 A–C 直接共讀、同一 molecule、同一 cell、同一 clone 或祖先關係。`k>12` 也不是「密集區域」的正式定義：目前 strict route 保留原始 `W` 為一個 source region，再切成 `k≤12` 的 computational blocks；跨 block 只能報 local topology，不能宣稱 global tree。（影響：高；信心：高）**

## 1. 背景摘要

目前 repository 至少同時存在三種「方法」與兩個過渡版本。它們回答的問題不同：

1. 2026-07-04 formal spec 定義的是 application-agnostic Boolean-lattice 最小相容樹集合；
2. 2026-07-15 current-v5 在舊 geometric region 上，對 formal solver 候選另加描述性的 read-AF 順位；
3. 2026-07-16 M2 嘗試用 exact molecule likelihood 比較不同最小 vertex sets，但正式 full run 尚未完成；
4. 2026-07-22 exact-PS k≤12 是 HCC1395 partial pilot；
5. 2026-07-23 strict endpoint route 正在把「區域怎麼形成」改成 exact PS×HP 內的 read-link graph connected component，目前仍是 in-progress production route。

所以「region、edge、k、cap、selection」不能只說一套而省略版本。若論文把 current-v5 的 read-AF ranking、M2 likelihood、strict endpoint region 和 formal spec 的「禁止 rank」寫成同一個 finalized 方法，會形成可被口委直接指出的內部矛盾。

本稽核不替使用者指定哪條是論文 SoT；以下把需要決策的選項與安全措辭明列。

| 路線 | 成熟度判定 | 現在可扮演的角色 |
|---|---|---|
| A formal spec＋old solver | formal contract 候選；`k≤8` tractability bound | 定義 structural candidate set；不能定義 genomic region |
| B current-v5 | 已完成的 legacy descriptive baseline | 舊 unit semantics 的比較基準，不是 calibrated truth |
| C1 M2 | v2 historical NO-GO；v4 current PROBE/WAIT | in-progress exact likelihood research |
| C2 exact-PS k≤12 | HCC1395-only partial pilot，技術完成 | side-route 技術證據，不是 cohort/paper final |
| C3 strict endpoint | production route in progress | 最新 region contract 候選；all-7 final receipt 尚待完成 |

## 2. 核心術語表

| 術語 | 本稽核採用的精確意思 | 不可直接等同 |
|---|---|---|
| fixed `R/A` call | molecule 在指定 sSNV endpoint 通過品質門檻，且觀測為 REF 或 ALT | clone identity |
| direct endpoint edge | 同一 canonical molecule 在兩 endpoint 皆為 fixed R/A；support 達門檻 | 祖先—後代 edge |
| read-linked region `W` | 同一 exact PS×HP 容器內，合格 endpoint edges 的 maximal connected component，且 `k≥2` | 一個真實 clone、單一 molecule 橫跨整區 |
| computational block `B` | `k>12` 的 `W` 為配合 solver 上限切出的 `k≤12` 局部運算單位 | 新的生物 region |
| `k_total` | `W` 或 `B` 中列入分析的 locus 數 | observed-ALT effective k |
| `k_active` | M2 observed-ALT model 中至少有 pattern 明示 ALT 的 active bit 數 | source region 的全部 locus 數 |
| linkage graph | node 是 sSNV locus；edge 是 molecule endpoint co-observation | Boolean-state tree |
| mutation-state graph | node 是 `k` 位點的 R/A state vector；unit-flip edge 是一次 0→1 state change | read-link edge、已確認 clone tree |
| structurally co-optimal candidates | 同一 formal objective 下的全部最小相容候選 | 「等機率」候選；若沒有 probability model，不應說等機率 |

Strict route 已正式區分 linkage graph 與 mutation-state graph；前者決定哪些 loci 可共同分析，後者才枚舉 state trees（`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md:28-37`）。

## 3. 版本矩陣

### 3.1 A — Formal spec v5＋`tree_enumeration_solver`

| 維度 | 目前文件／實作 | 判讀 |
|---|---|---|
| 狀態 | 文件自標 `formal spec v5`、Model A 定案（`InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md:1-14`） | 是 formal contract 候選，不是 region-building production spec |
| region | formal spec 只接收 `O⊆{0,1,?}^k`，沒有 genomic region membership 規則（同檔:20-35） | 不能用它回答 50 kb、PS、HP 或 A–B–C 如何成區 |
| tree edge | Boolean lattice 的合法 edge 只翻一個 bit `0→1`（同檔:26-30） | 不是 read endpoint edge |
| `k` | 文件硬上界 `k≤8`、狀態空間 `2^k≤256`（同檔:96-113） | 與後來 `k≤12` route 尚未正式對齊 |
| 最小性 | 先覆蓋全部合格觀測，再依 hidden nodes、edge 等字典序最小化，列出全部同分最佳（同檔:50-55、84-94） | 輸出 candidate set，不強迫唯一 |
| selection | 明文「enumerate，非 optimize，非 rank」；禁止把不可辨識集強行全序（同檔:115-121） | 不允許把 read-AF／likelihood winner 偷寫成 formal solver 本身 |
| solver cap | 函式預設 `extra_cap=4`、`per_level_budget=150000`、`tree_cap=32`（`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py:99-108`） | `tree_cap=0` 只解除儲存樹數上限，不解除 combinatorial budget |
| capped 結果 | layer 超 budget 時轉 greedy fallback 單樹並標 `capped`；不能稱完整（同檔:138-165、188-207、242-257） | fallback tree 不是全部最小樹 |
| effective universe | solver 只把曾被完整／部分 pattern 明示為 `A` 的 bit 放進 `universe`（同檔:120-121） | 和 nominal `k` 的關係需在論文明示 |

**狀態結論**：A 是仍具價值的 formal model，但其可計算性主張綁定 `k≤8`。後來拿同一 solver 跑 `k≤12`，若仍引用「`2^k≤256` 因此可完整窮舉」，就是不成立的版本混用。`enumerate_min_trees(..., k, ...)` 的入口到 universe 建立前沒有 `k≤8` 拒絕條件（`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py:99-121`），因此「程式能跑」不等於「formal tractability claim 已擴充」。

### 3.2 B — 2026-07-15 current-v5 read-AF topology

| 維度 | 定義 | 證據 |
|---|---|---|
| 狀態 | 7 datasets comprehensive descriptive workstation 已完成 | `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/20260715_GRCh38拓撲形態多選工作站完整檢視_01.md:14-21` |
| region | 舊 route 是相鄰 sSNV gap≤50 kb 的 connected component，總 span 可 >50 kb；每 region×HP family 處理 | `InterSubMod/docs/methodology/20260709_quantified_topology_clone_subclone_overview_01.md:65-71` |
| `k` | 舊每區送 solver `≤8 sSNV` | 同檔:70 |
| candidate enumeration | `tree_cap=0` 重枚舉 candidate universe | `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py:382-413` |
| selection | `score=Σ(AF_ancestor−AF_newly-acquired)`，同 HP unit 內 exact tie | `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/20260715_GRCh38拓撲形態多選工作站完整檢視_01.md:43-51`；實作同 topic script:123-138、438-466 |
| output ceiling | regional mutation-state candidate topology；不是 posterior、CCF、confirmed clone 或 true ancestry | 同報告:21、25-30 |

**狀態結論**：B 不是「未完成方法」，而是已完成的**舊 unit semantics 描述性 baseline**。但 strict route 已明確指出舊 production grouping 先用 50 kb、又可能跨 exact PS 聚合，不符合新的 region 問題（`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md:18-24`）。因此它可作 baseline，不宜直接稱為目前 production region SoT；是否保留為論文主要方法，待使用者決定。

### 3.3 C1 — 2026-07-16 M2 exact likelihood

| 維度 | 定義／狀態 | 證據 |
|---|---|---|
| 正式狀態 | frozen v2 pilot `NO-GO`；第一個 ranking 8 小時 timeout，full 154 tasks 未啟動 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/00_INDEX.md:21-29` |
| current route | v4 為 `PROBE`，formal actual-data `WAIT`，雙 pilot 通過前不能進 full | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260719_M2_exact_preserving資源閘門修補正式執行Runbook_04.md:14-18` |
| region boundary | primary unit 是 HP×exact PS×read-linked component；M2 文件用「相鄰 site cut 兩側皆有 fixed call」的 cut-support 定義 bridge | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/00_INDEX.md:61-68、109-117` |
| `k` | gate 用 `k_observed_alt_active`；`≤8` exhaustive oracle、9–12 MILP、`>12` local-only／abstain | 同檔:76-82 |
| candidate object | 先依 vertex set 合併；同 vertex set、不同 parent edges 不可由 state likelihood 區分 | 同檔:67-72 |
| selection | error-aware molecule-pattern likelihood，只比較不同 vertex sets | 同檔:69-74 |
| output ceiling | solver-certified regional mutation-state candidate sets；不是真實 clone 數、唯一 parent edge 或完整演化史 | 同檔:157-165 |

**狀態結論**：M2 的 likelihood 不是 B 的 read-AF heuristic；兩者不能都簡稱「VAF 排序」。M2 是設計較完整但尚未完成正式真實資料 full validation 的研究路線，不能當成目前已交付結果。

### 3.4 C2 — 2026-07-22 exact-PS k≤12 HCC1395 pilot

| 維度 | 定義／狀態 | 證據 |
|---|---|---|
| 狀態 | `PARTIAL / exploratory pilot / HCC1395 chr1-22 only` | `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/pre-decision-audit.md:1-14` |
| unit | exact PS×HP×upstream read-linked component | 同檔:22-30 |
| block | `k>12` 在同一 unit 內切成 contiguous、non-overlap `k≤12` blocks；block 不能宣稱無損拼回 global tree | 同檔:28-30、57-59 |
| partition objective | ordered-hypergraph DP 最大化完整保留於同一 block 的 constraint weight | `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/read_support_partition.py:1-21、224-235` |
| tree selection | HCC1395 local enumeration/T/Topo 已跑；VAF、count-weighted likelihood、CN、methylation 未完成 | `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/implementation-notes.md:191-210、296-303` |
| claim ceiling | 技術 checks PASS，但上游 receipt 仍 `validation_evidence_eligible=false`，不可作 cohort/paper final | 同檔:296-303 |

**狀態結論**：C2 是已完成 HCC1395 技術驗證的 side-route pilot，不是 finalized cohort method。它的 upstream component 是 cut-span semantics，並非 C3 的 strict endpoint component；兩者 unit 數不可直接相減。

### 3.5 C3 — 2026-07-23 exact-PS strict endpoint production route

| 維度 | 定義／狀態 | 證據 |
|---|---|---|
| 狀態 | `in_progress`；HCC1395 完成，七資料集未全完成 | `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md:1-14、169-175` |
| container | dataset×chromosome×exact nonmissing PS×primary HP | 同檔:18-24、39-48 |
| direct edge | 同 canonical molecule 對兩 endpoint 都 fixed R/A，distinct support≥3 | 同檔:20-24、50-56 |
| source region `W` | 合格 edges 的 maximal connected component，`k≥2`；singleton abstain | 同檔:39-46 |
| block `B` | `k>12` 的 `W` 切成 `k≤12` computational blocks；不是新 region | 同檔:43-46、114-126 |
| selection | 目前這一層只完成 grouping／segmentation；尚不能宣稱 topology、VAF、CCF、CN 或 methylation 支持 | 同檔:131、177-181 |
| threshold | primary=3 是工程設定，已做 1/2/3/5 sensitivity，但不是完成統計校準 | 同檔:103-112 |

**狀態結論**：C3 是最接近使用者最新口語描述的 region contract，但目前仍在 production promotion 中，不應稱「七樣本終版」。同日進度文件也有時間差：主報告仍寫 HCC1395＋132 tasks 待補（同報告:14、169-175），living notes 已寫四個 dataset summary validated、三個待處理（`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/implementation-notes.md:9-15`）。論文引用進度前應以最終 immutable receipt 重查，而不是自行挑其中一個數字。

## 4. 互相矛盾或容易誤用的 claims

### 4.1 Formal「不排序」對上 current-v5／M2「排序」

- formal spec 禁止對不可辨識樹集建立 total order（`InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md:115-121`）。
- current-v5 確實依 read-AF score 排 rank（`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py:438-466`）。
- M2 則以 error-aware likelihood 比較不同 vertex sets（`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/00_INDEX.md:67-72`）。

這三句若都放在「本研究方法」而沒有分層，就是矛盾。可選的安全整合方式是：formal structural solver 只產生 co-optimal set；read-AF 只列為 historical/descriptive post-hoc sensitivity；M2 likelihood 列為 in-progress future validation。是否採此整合，需使用者定案。

formal spec 自身也有措辭衝突：§7.5.3 說不做 total-order（同檔:115-121），§8 卻說通道 M 可做「相容樹集內弱排序」（同檔:133-140）。若「弱排序」只指 annotation／partial preference，應改名並說不移除候選；若會決定第一順位，就已跨過前述紅線。

### 4.2 「等機率」沒有 probability model

formal spec 把 ambiguous 集合寫成「等機率相容集」（同檔:123-129），但該層沒有定義 prior 或 calibrated likelihood。論文宜改為：

> 「全部 structurally co-optimal compatible candidates；本結構層不賦予候選機率。」

### 4.3 `k≤8` 與 `k≤12` 不是只改一個常數

- formal tractability 論證明確綁定 `k≤8`、`2^k≤256`（formal spec:96-113）。
- exact-PS pilot 把 segmentation 上限改為 12，並明示是工程界線（`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/pre-decision-audit.md:22-30、72-77`）。
- adapter 雖用 `ANALYSIS_TREE_CAP=0`，實際仍有 420 units 被 solver 標 capped（`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/implementation-notes.md:203-210`）。

因此，`k≤12` 表示「允許送入 solver 的 block 大小」，不表示「所有 `k≤12` 一定完整枚舉」。需同時報 `capped`、budget 與 complete flag。

### 4.4 M2 的 `k_active` 與 strict route 的 `k_total` 不同

M2 用 observed-ALT effective k（`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/00_INDEX.md:76-82`）；strict report 的 k 是 `W` 中 loci 數（strict report:39-46）。總 k>12 但 active-k≤12 在 M2 model 內可能仍 exact；strict production route 則會依總 k 切 block。兩者都合理，但回答的是不同 model，論文不可只寫「k>12 就切」而不標 k 的定義。

### 4.5 C2 cut-span component 與 C3 strict endpoint component 不同

production pre-decision audit 已量化兩者不等價：cut-span threshold=3 不等於 pairwise endpoint threshold=3，HCC1395 有 component／block 會裂開（`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/pre-decision-audit.md:27-36`）。所以 C2 的 39,544 components／11,542 route-blocks與 C3 的 11,462 `W`／11,712 `B` 不可直接作 before-after improvement；C2 自己也說新舊分母不能直接相減（`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/20260723_HCC1395_exactPS拓撲重建觀察_01.md:71-87`）。

### 4.6 `RR` edge 的研究語意尚需決策

strict graph 把 RR、RA、AR、AA 都算作 fixed endpoint support；報告也明示 RR-only 只證明共同 callable 的 reference state，不證明 somatic ALT 共現（strict report:95-99、177-181）。

- 若 region 的問題是「哪些 loci 有足夠同分子可觀測關係，可放入共同分析容器」，RR edge 可保留。
- 若論文主張是「mutation-bearing molecules 的 somatic co-occurrence」，RR-only edge 太寬，至少要做 `ALT-informative edge` sensitivity（RA/AR/AA，或明定至少一 endpoint 為 ALT）。

這是高優先級、尚不能默認的定義決策。

## 5. `k>12`：切開，還是當特別密集區域？

### 5.1 目前 strict route 的精確答案

兩個動作同時成立，但在不同層：

1. **科學分析區域層**：不切掉 source identity。原始 `W` 仍是一個 read-linked candidate region。
2. **運算層**：若 `k_total>12`，把 `W` 分成若干 `k≤12` computational blocks `B_ij`。
3. **解讀層**：每個 block 的 tree 是 local result；保留 `source_component_id` 回指 `W`，但不能把多個 local trees 自動拼成一棵 global tree。

此合約已明寫在 strict report：`B` 不是新 region（`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md:39-46`），HCC1395 的 90 個 `k>12 W` 經切割後仍是 90 個 source regions，而不是把 block 數當 region 數（同檔:114-131）。

### 5.2 為何不應直接叫「密集區域」

repository 目前沒有一個跨版本一致、具有量綱與門檻的「密集區域」正式定義。現有「密集／太密」至少指三件不同事情：

1. historical `densest-8`：在連續 8 個 locus windows 中，選 coordinate span 最小者（`InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/build_k_gt8_partitions.py:222-229`）；
2. old solver `太密`：某 hidden-state layer 的組合數超過 `per_level_budget`，是計算爆量，不是 bp 密度（`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py:138-165`）；
3. strict `large region`：只表示 `k>12`，而 strict membership 不由 genomic distance 決定（strict report:80、95-99、114-126）。

因此安全用詞是：

> 「large-k read-linked region（`k_total>12`）」  
> 或「超出目前 exact/local solver block 上限的 read-linked region」。

若要正式使用「密集」，需另定 metric，例如 `loci per kb`、median adjacent gap、max adjacent gap 或 graph edge density，並預先給 threshold；`k>12` 本身只代表 loci 數超過工程上限，不代表空間密度、clone 數或演化複雜度。

### 5.3 future exact route 的另一選項

M2 efficiency audit 指出 `k_total>12` 不應永遠切：若 observed-ALT active-k≤12，仍可能 exact；若要對真正 large active-k 做 global exact decomposition，必須保存 separator boundary assignment 與 certificate，任意分塊後各解一棵再拼不具 global certificate（`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_k大於12_exact效率與certified策略_01.md:93-110、114-125、155-180、239-245`）。

這是 future certified router 選項，不是目前 strict production 已完成能力。

## 6. A–B、B–C 串聯是否合理？

### 6.1 合理的層次：region membership

令 direct edge relation 為：

```text
E(i,j) = 同一 exact PS×HP 容器內，
         至少 threshold 條 distinct canonical molecules
         對 i、j 皆有 fixed R/A call
```

如果 `E(A,B)=true` 且 `E(B,C)=true`，connected-component 演算法會把 `{A,B,C}` 放在同一 `W`。這是標準的無向圖 connectedness，不需要 `E(A,C)=true`。實作以 DSU 對所有 support≥threshold edges union（`InterSubMod/tools/strict_endpoint_graph.py:274-295`）；合成測試也明確驗證 A–B 與 B–C 形成一個 component、同時 A–C edge 不存在（`InterSubMod/tests/test_strict_endpoint_graph.py:26-33`）。

所以以下說法正確：

> 「A–B 與 B–C 各自有足夠 direct read evidence，因此三個位點屬於同一個 read-linked connected component，可作為同一局部候選區域共同分析。」

### 6.2 不合理的跳躍：把 graph transitivity 當 molecule／clone transitivity

A–B edge 與 B–C edge 可以由完全不同的 molecule 集合支持；因此不可推出：

- 有一條 molecule 同時 fixed-observe A、B、C；
- A 與 C 有 direct read evidence；
- A、B、C 同屬一個 cell；
- 三個突變一定同屬一個 clone；
- A→B→C 是祖先順序；
- 三個 mutation state 一定共同存在於同一 biological lineage。

strict 方法文件已逐字寫出「兩條 edge 可以由不同 molecule 集合支持」以及 `W` 不代表單一 read／cell／clone 橫跨所有 k（`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md:58-64`）。

精確講法是：

> 「傳遞的是『同一 connected component 的分析資格』，不是『分子共現、細胞共現或 clone 身分』。」

### 6.3 到 computational block 之後會發生什麼？

- 若 source `W` 的 `k≤12`，pilot partition 會保留為單一 block（`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py:30-32、771-776`），A–B–C 可一起送入 local solver。
- 若 `W` 的 `k>12`，contiguous DP 可能把 A–B 或 B–C constraint 切在不同 blocks。DP 只保證在 `k≤12` 非重疊 blocks 下最大化完整保留的 constraint weight，並逐條標 retained/cut/unavoidable（`InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/read_support_partition.py:1-21、224-235、345-374`）。
- 一旦 B 位於 block boundary，不能因 source `W` 仍連通，就宣稱各 block 的 local trees 已形成 global tree。要做 global claim，需 exact separator/boundary-state 合併 certificate（M2 k>12 audit:114-125）。

### 6.4 mutation-state solver 能說到哪裡？

把 A、B、C 放入同一 `W` 只決定「共同分析範圍」。之後 Boolean-state solver可根據所有完整與部分 R/A patterns，列出滿足 coverage axiom 的最小相容 state-tree candidates（formal spec:20-35、50-55、84-94）。

因此可說：

> 「在目前觀測與模型下，A、B、C 被納入同一局部候選問題；solver 會列出能共同解釋這些 observations 的最小相容 mutation-state candidates。」

不可說：

> 「因為 A–B、B–C 能串，所以三者已被證明是同一 clone 的最小組合。」

「最小」是模型 objective 下的 state-set/tree minimality，不是 biological clone identity 的證明。

## 7. 論文與口試可直接使用的精確措辭

### 7.1 Region formation

> 本研究先在同一染色體、exact phase set 與 primary haplotype family 內建立 sSNV read-linkage graph。每個節點代表一個 sSNV 位點；只有當至少三條不同 canonical molecules 在一對位點上都具有通過品質門檻的 REF/ALT 判讀時，才建立 direct endpoint edge。由這些 edge 形成、且至少包含兩個位點的 maximal connected component，定義為 read-linked candidate region。

### 7.2 A–B–C chain

> 若 A–B 與 B–C 各自有足夠 direct endpoint support，即使沒有 A–C direct support，A、B、C 仍可透過圖的連通性歸入同一候選區域。這個傳遞性只用於界定共同分析範圍；兩條 edge 可由不同 molecules 支持，因此不代表同一 read、同一 cell、同一 clone 或 A→B→C 的演化順序。

### 7.3 `k>12`

> 當一個 source read-linked region 含超過 12 個位點時，我們保留其 source-region 身分，但因目前 solver 的工程上限，將它切成 `k≤12` 的 computational blocks 進行局部重建。block 不是新的 biological region；若沒有跨 block 的 exact decomposition certificate，結果只能解讀為 local candidate topology，不能宣稱完整 global tree。

### 7.4 Candidate tree

> 對每個可計算 block，我們列出在既定 coverage 與 unit-flip model 下的最小相容 mutation-state tree candidates。若候選不唯一，我們保留整個 structurally co-optimal set，不把結構上的不可辨識性誤寫成唯一真實 clone tree。

### 7.5 最短口試版

> 「A–B 和 B–C 的 read evidence 可以把三個位點連成同一分析區域，但不能證明三個突變在同一個 clone。區域超過 12 個位點時，我們保留原始區域 ID，只為了運算切成小 block；目前只對 block 內做局部候選樹重建。」

## 8. 已知前提、未知與風險

### 已知前提

1. exact PS 與 HP 是目前 strict primary container boundary；跨 PS 不合併（strict implementation notes:18-28）。
2. coordinate span 或 `X/O/D/S/L` 不建立 endpoint edge（同檔:37-42）。
3. singleton 沒有 co-read relation，不進 tree denominator（同檔:44-54）。
4. 50 kb 只作 QC，不作 primary membership boundary（同檔:56-63）。
5. `W` 和 `B` 是不同 grain（同檔:30-35）。

### 尚未解決

1. **論文 structural SoT**：採 A 的 no-rank contract、B 的 descriptive read-AF rank，或等待 M2 exact likelihood？
2. **`k` 定義**：論文主 gate 用 `k_total` 還是 M2 `k_observed_alt_active`？
3. **`RR` edge**：作 general callability linkage 保留，或對 mutation-bearing claim 改用 ALT-informative sensitivity？
4. **large-k globality**：接受 local blocks＋明示 local-only，或開發 exact separator/router？
5. **formal spec 內部文字**：「弱排序」和「不排序」、「等機率」和「未定 probability」需統一。
6. **current progress**：C3 的 all-7 狀態須由最終 receipt 定案；目前報告與 living notes 是不同時間快照。
7. **threshold=3 校準**：目前是工程門檻；仍需 held-out/statistical edge confidence 或至少完整 sensitivity。
8. **graph robustness**：A–B–C chain 可能由單一 articulation edge 維繫；應輸出 articulation／bridge 與 leave-one-molecule sensitivity。此風險已在 pre-decision audit 明列（`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/pre-decision-audit.md:73-78`）。

## 9. 待使用者決策的五個 gate

| Gate | 選項 | 本稽核傾向 | 原因 |
|---|---|---|---|
| D1：論文主要 region | legacy 50-kb vs strict endpoint `W` | strict endpoint `W` | 定義直接綁 molecule evidence；但需等 all-7 final receipt |
| D2：候選 selection | no rank vs read-AF rank vs M2 likelihood | 主文 no-rank；read-AF 作 historical sensitivity；M2 作 future/in-progress | 最符合現有驗證成熟度與 claim ceiling |
| D3：k>12 | source region 全丟／改名密集／保留 W 並切 B | 保留 `W`、切 `B`、local-only | 與現行 strict contract一致且不混 grain |
| D4：k 的門檻 | `k_total` vs `k_active` | 論文兩者都報，主 route 明示採哪個 | 兩者代表不同 model，不能混用 |
| D5：edge allele semantics | RR+RA+AR+AA vs ALT-informative | primary可保留全 fixed-call；生物 mutation claim 必做 ALT-informative sensitivity | RR-only 不支持 ALT co-occurrence |

這些是建議，不是本稽核自行指定的 SoT；最終需由使用者確認。

## 10. 建議閱讀順序

1. `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md`  
   最新 region、edge、W/B 與 A–B–C 語意。
2. `InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md`  
   Boolean-state tree、coverage、minimal candidate set 與 no-rank contract。
3. `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/00_INDEX.md`  
   M2 exact likelihood、effective k、current NO-GO/PROBE 狀態。
4. `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_k大於12_exact效率與certified策略_01.md`  
   k>12 exact/local certificate 邊界。
5. `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/20260723_HCC1395_exactPS拓撲重建觀察_01.md`  
   HCC1395 partial topology output 與未完成項目。
6. `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/20260715_GRCh38拓撲形態多選工作站完整檢視_01.md`  
   current-v5 descriptive baseline 與 read-AF claim ceiling。

## 11. 本次稽核的輸入、命令與輸出

### 輸入

- `InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py`
- `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/`
- `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/`
- `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/`
- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/`
- `InterSubMod/tools/strict_endpoint_graph.py`
- `InterSubMod/tests/test_strict_endpoint_graph.py`

### 執行命令

```bash
rg -n "k≤8|k>12|tree_cap|densest|密集|read-AF|likelihood|connected component|weak ranking" \
  docs/methodology research/20260715_layered_workstation_genome_topology_multiselect \
  research/20260716_read_linked_hypercube_exact_likelihood_validation \
  research/20260718_k_gt8_read_supported_segmentation \
  research/20260722_exact_ps_k12_hcc1395_pilot \
  research/20260723_production_exact_ps_strict_read_linkage
nl -ba <上述命中文件> | sed -n '<證據行範圍>'
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  -m pytest -q tests/test_strict_endpoint_graph.py
```

### 輸出

- `InterSubMod/research/20260723_k_gt12_chain_validity_audit/agents/version_claim_audit.md`

### 實際可觀察結果

- 找到 A、B、M2、exact-PS pilot、strict production 五個不同成熟度／語意層。
- 找到至少六組不可無標示混用的 claim：no-rank vs rank、k≤8 vs k≤12、`k_total` vs `k_active`、cut-span vs endpoint component、region vs block、fixed-call linkage vs ALT co-occurrence。
- 找到 A–B／B–C chain 的直接程式與測試證據；它只保證 connected component membership。
- 未找到跨版本一致的「密集區域」正式定義。
- `tests/test_strict_endpoint_graph.py` 實際執行結果：`13 passed in 0.06s`；其中 A–B＋B–C→單一 component 且 A–C edge 不存在的 fixture 通過。
