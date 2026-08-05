<!--
建立時間: 2026-07-16 06:18 +08:00
目標: 記錄 read-linked region、symbolic partial-read Hypercube Steiner、k>8 exact solver 與 lossless read-pattern likelihood 的設計、偏離、折衷及驗證結果
處理範圍: chr1-22、7 datasets / 6 biological samples、LongPhase-S recalibrated FILTER=PASS canonical v5
cycle_id: cycle_20260716-0618-read-linked-hypercube-exact-likelihood
spec_id: read_linked_hypercube_exact_likelihood_validation
status: in_progress
advisory: on
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md
  - InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md
-->

# Implementation Notes：read-linked Hypercube Steiner 全量驗證

> **Purpose**：把每項資料口徑、solver objective、不可辨識性與全量驗證收據保留下來，避免把工程最優誤寫成真實clone ancestry。
> **Task type**：B Comprehensive validation；全 chr1–22 × 7 datasets。服務 G3／G4／G5。

## 🔵 設計決定（Design Decisions）

### 2026-07-16 06:18 — partial read以symbolic subcube處理，不做2^u實體展開

- **Status**：Accepted
- **決定**：每個pattern保存 `(fixed_mask, alt_mask, free_mask)`；candidate vertex `v`只需檢查 `((v XOR alt_mask) AND fixed_mask)==0`。每個partial group必須被候選vertex set至少一個相容state覆蓋。
- **理由**：語意與列舉所有X補值後取聯集相同，但每個pattern只需一條group constraint，不建立 `2^u` 個terminal，也不把一個read複製成多筆獨立證據。
- **實作細節**：現行MILP為該單一group的sparse hit row列出至多`2^u_eff`個相容active-vertex indices；這是同一條限制式的nonzeros，不是`2^u_eff`次獨立求樹，也不會跨reads做Cartesian product。
- **驗證**：k=1-8共9,840 patterns、2,015,538個state checks，explicit與symbolic相容集合0 mismatch。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-16 06:18 — 先求minimum-hidden可行state set，再以read likelihood比較不同vertex sets

- **Status**：Accepted
- **決定**：結構層先固定root、full observed states、partial group coverage與monotone unit-step reachability，lexicographic第一目標為hidden-node數最少；資料層只在同一結構複雜度或明示penalty後比較error-aware read-pattern likelihood。
- **理由**：partial read只能說「真實state屬於相容集合」，通常不能直接支持某一條parent edge；candidate generation與data discrimination需分層。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-16 06:18 — 相同vertex set、不同parent edges一律標記不可由bulk snapshot辨識

- **Status**：Accepted
- **決定**：likelihood以mutation-state mixture為單位；先以`vertex_set_id`合併候選。若只差parent assignment，分數必須相同並輸出`EDGE_NONIDENTIFIABLE`。
- **驗證**：synthetic same-vertex/different-edge score差 `1.42e-14`，數值容差內為0。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-16 06:18 — 結構evidence與likelihood evidence保留兩層

- **Status**：Accepted
- **決定**：`structural_patterns`可沿用MINREAD>=3作硬結構門檻；`scoring_patterns`必須保留門檻前所有canonical molecule patterns，低頻pattern交給error model處理。
- **理由**：用MINREAD後aggregate做likelihood會把低頻資料審查掉；只能稱thresholded alignment-exposure profile，不是lossless molecule likelihood。
- **Evidence tier**：L2 ⭐⭐⭐⭐

<!-- BEGIN USER-SPECIFIED -->
**Decision**：最終驗證必須包含全chr1-22、所有7 datasets、清楚的最後數字與各分母；partial-read做法必須確認是否等同多種補值、如何輸出候選與如何後續排序。
**DO NOT change**：subset只可作pilot，不能代替最終數據；不能把可行tree直接稱為真實clone tree。
**Rationale**：2026-07-16 user request。
<!-- END USER-SPECIFIED -->

## 🟠 偏離之處（Deviations）

### 2026-07-16 06:18 — 不把read數均分到流入邊，也不將VAF重複加入同read likelihood

- **Status**：Accepted
- **規範／原始想法**：把partial-read數量分配給相容state或流入邊，再用加權Steiner直接選樹。
- **實作偏離**：每個read只進一次marginal likelihood；相容state由mixture權重求和，不做人為均分。VAF若由同一批reads計算，只作posterior predictive check，不再加分。
- **理由**：均分會假設未知state等機率且把snapshot state evidence誤當edge evidence；read likelihood與VAF同源時相加會double count。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

## 🟡 折衷考量（Trade-offs）

### 2026-07-16 06:18 — exact optimality與biological stability分開報告

- **Status**：Accepted
- **方案 A**：solver保存incumbent／bound／gap／status，回答「固定模型下是否最優」。
- **方案 B**：read bootstrap、HP／bridge threshold、error rate、CN-clean sensitivity，回答「資料擾動後是否穩定」。
- **理由**：gap=0不能證明資料足以唯一辨識，也不能證明生物演化為真。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

### 2026-07-16 06:18 — 全量BAM掃描延後到既有Task-B釋放資源

- **Status**：Accepted
- **決定**：先完成低資源solver、schema、synthetic與canonical aggregate cross-check；不與現有48-core全量任務競爭。之後才依7-dataset manifest排程lossless extraction。
- **理由**：避免兩個全量任務互相拖慢或產生不完整輸出；不是縮小最終scope。
- **Revisit if**：現有Task-B完成且preflight確認CPU／disk安全。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐

## 🔴 未決問題（Open Questions）

### 2026-07-16 06:18 — read-linked boundary與k>12 decomposition全量穩定性

- **Status**：open
- **Question**：bridge support門檻1／2／3／5下區域數與fragmentation如何；k>12能否用articulation／separator或overlapping exact subproblems保存全域相容性與certificate？
- **Resolve by**：7-dataset full molecule table、threshold sensitivity、synthetic adversarial tests與full terminal receipts。
- **Priority**：P0
- **Evidence tier**：L5 ⭐

### 2026-07-16 06:18 — CN／LOH／purity後的likelihood claim ceiling

- **Status**：open
- **Question**：哪些regions可在copy-number-clean條件下比較mixture likelihood；哪些必須`ABSTAIN_CN`？
- **Resolve by**：接入具provenance的CN/LOH track並做分層sensitivity；不得用raw VAF單調性代替。
- **Priority**：P1
- **Evidence tier**：L5 ⭐

## ✅ 實作與驗證紀錄

### 2026-07-16 06:18 — symbolic／MILP／likelihood先導驗證

- **Input**：canonical v5 frozen solver與7-dataset MLHP JSON；synthetic k=1-12 fixtures。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_exact_symbolic_pilot.py --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot --canonical-per-dataset 5`
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot/pilot_receipt.json`，SHA256 `9939a9aeb89758d17657876089faac20f0e4fb248dfdb2091b61fea54149d8b0`。
- **Actual**：exit=0、`all_pass=true`；symbolic 2,015,538 checks零差異；legacy/MILP 66 vertex-set checks零差異；9 legacy-capped案例契約性SKIP；k=9-12為12/12 PASS、max variables=4,096、max runtime=0.225s；likelihood controls全部PASS。

### 2026-07-16 06:18 — 單元測試

- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/{scripts,tests}`。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`
- **Output**：terminal receipt。
- **Actual**：exit=0；9/9 tests PASS；包含partial masks、all-X常數likelihood、same-V edge invariance與EM monotonicity。

### 2026-07-16 07:02 — canonical M0 分母與候選集合稽核

- **Input**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/samples/*/layered_region_view_*.json` 與 frozen solver source bundle。
- **Command**：唯讀 schema／join／count 稽核；display 截斷案例以 frozen solver 重建完整 analytical vertex sets。
- **Output**：本文件與進行中的 `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/`。
- **Actual**：72,994 個 mutation-bearing HP1/HP2 units；其中 capped 8,021（10.99%），eligible non-capped complete 64,973（89.01%），來自 46,385 regions。eligible 中原始 `T=1` 26,225、`T>1` 38,748；61,702 個未受 display cap 影響的 units 中，26,326 個只有一種 vertex set、35,376 個有多種 vertex sets，另 3,271 個需由 frozen solver重建。M0 thresholded patterns 保存 7,086,491／7,341,571 alignment exposures（96.5255%），漏掉 MINREAD<3 的 255,080（3.4745%）；因此只能稱 M0，不能稱 lossless M2。

### 2026-07-16 07:02 — COLO829 M0 全樣本 performance gate

- **Input**：canonical v5 COLO829，11,401 eligible HP units。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time -v python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_m0_likelihood_census.py --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_colo829_v3 --sample COLO829`。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_colo829_v3/`；rows SHA256 `b3e964...`。
- **Actual**：exit=0、elapsed=4m50.71s、peak RSS=274,616 KB。分類為 `T1/V1` 4,290（37.6283%）、多 vertex sets likelihood 並列 6,572（57.6441%）、likelihood 唯一第一 438（3.8418%）、optimizer fail-closed 101（0.8859%）；沒有把 non-convergence 誤列為唯一第一。

### 2026-07-16 07:02 — lossless M2 extractor 與資源 gate

- **Input**：canonical v5 manifest、tree VCF、約 1.806 TiB 的 7 dataset BAM、LongPhase-S sidecars。
- **Command**：先以 HCC1954 chr22、1 samtools thread、`nice -n 10`／`ionice` 執行 performance pilot；全量 runner 僅作 preflight。
- **Output**：進行中的 `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/`；正式全量將落在 `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_chr1_22_7datasets/`，目前尚未建立，避免把大型正式輸出放入repo。
- **Actual**：extractor 已保存 raw alignment funnel、逐 CIGAR 的 R/A/O/D/S/L/X、BQ、exact sidecar join、unique-molecule bridge threshold 1/2/3/5、component membership、hash receipt。資源 gate 已修正為匹配 canonical 與 `.pinned_<sha>.py` 命名；目前偵測到使用者 44-worker all-sSNV 工作，故 full M2 必須等待，不能以 pilot 代替全量結果。

### 2026-07-16 07:02 — 單元測試擴充

- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/{scripts,tests}`。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Output**：terminal receipt。
- **Actual**：exit=0；29/29 tests PASS。

### 2026-07-16 07:10 — partial-read 方法獨立反例稽核

- **Input**：新 symbolic MILP、canonical-v5 frozen solver、`RRA + AAA + RAX` synthetic case。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py -v`，並用兩個 solver及 k=3 exhaustive brute force交叉驗證。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md`，SHA256 `e9ac09bea2ad06b3a0ca3daf9b6d0c0af3a00e20ed3e330249852df4672fed2a`。
- **Actual**：10/10 tests PASS。`RAX={RAR,RAA}`；全域唯一 minimum-hidden vertex set為 `{RRR,RRA,RAA,AAA}`、hidden=1。若強迫 `RAR` 雖可行但 hidden=2，直接證明「逐 completion、任一成功就停止」錯誤；正確做法是整個 group constraint共同求全域 optimum。

### 2026-07-16 07:15 — HCC1954 chr22 M2 extraction pilot v1

- **Input**：HCC1954 tree VCF 220 chr22 sSNV、253.2 GB BAM（storage identity）、635.3 MB LongPhase-S read-tag sidecar。
- **Command**：`nice -n 10 ionice -c2 -n7 /usr/bin/time -v python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/extract_lossless_read_linkage.py --manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json --sample HCC1954 --chrom chr22 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22 --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5 --samtools-threads 1`。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/`。
- **Actual**：exit=0、all_pass=true、elapsed=16m00.77s、peak RSS=93,444 KB。raw overlaps 21,540 = flag-excluded 5,506 + MAPQ-rejected 363 + eligible 15,671；15,671/15,671 exact sidecar join、missing=0。18,938 sparse calls = R 11,586 + A 4,951 + O 8 + D 212 + S 0 + L 2,181 + X 0。156 個相鄰 gap>50 kb cuts中有3個仍有read bridge，threshold≥3仍有1個（gap 52,478 bp、3 molecules），直接否定「50 kb 是ONT無法連接上限」。此v1只有pooled components；QA要求補per-HP basis後另跑v2，故不可作final全量比例。

### 2026-07-16 07:20 — M2 extractor schema 1.1.0 correctness／resource／receipt QA

- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{lossless_read_contract.py,extract_lossless_read_linkage.py,run_full_m2_extraction.py}`、對應 tests、tiny synthetic BAM／VCF／tabix sidecar。
- **Command**：`python3 -m unittest discover -v -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py'`；另以 full runner `--preflight-only` 對 canonical v5 manifest 執行資源 gate。
- **Output**：terminal test receipt；`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2_extractor程式與合成資料QA_01.md`；extractor SHA256 `3d62c68d8e97818224c0ecb66ced2ab4d5c9a5936a3882316cc3f3df8b49f446`。
- **Actual**：exit=0、50/50 tests PASS。schema 1.1.0 將每個 threshold 分別輸出 `pooled/HP1/HP2/HP3/HP4/unphased` component partition；raw alignment class、flag/MAPQ funnel、written rows、R/A/O/D/S/L/X call mass、fixed R/A 與 ALT mass 全守恆；unknown HP 與同 molecule 多 primary alignment 均 fail-closed。`receipt.json` 改由 exact-byte `.sha256` sidecar驗證，resume reuse另驗 scope、參數、manifest/extractor hash與每個 output hash。資源 preflight exit=2，實際偵測50個衝突process並聚合為2個root；未建立正式output root。

### 2026-07-16 07:25 — M2 symbolic pattern／BQ likelihood／full-ranking契約

- **Input**：schema 1.1.0 chromosome extraction的sparse molecule calls、site catalog、`linkage_basis × threshold` component membership與exact-byte receipt。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py`、`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_ranking.py`、`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_m2_patterns_and_rank.py`。
- **Actual**：exit=0；53/53 tests PASS，其中新增12項涵蓋`RAX` symbolic group、all-X守恆排除、MINREAD只作用於structure、same-V edge不可辨識、k>12 abstain、per-call BQ emission、O/D/S/L/X missing policy、per-HP component basis、seeded molecule bootstrap、topology coarse class與154-task分層aggregation。
- **Ranking contract**：primary score為每筆R/A fixed call的Phred-BQ conditional-biallelic emission；O/D/S/L/X保留reason count但以missing marginalize。fixed-error grid是sensitivity，與optional deterministic molecule bootstrap分開；同reads衍生VAF不重複加分。family 1合併`1/1-1/1-2`、family 2合併`2/2-1/2-2`，somatic fine tag不可作獨立evidence或clone label。
- **Topology contract**：輸出四類coarse class，但`coarse_topology_class_unique`不等於canonical exact shape unique；`exact_topology_unique`只在唯一winning parent-edge assignment時證明true，多edge情況標NA等待shape-isomorphism，禁止以V或coarse class冒充Topo唯一。
- **Full schema**：`intersubmod.m2_full_ranking_receipt` v1.0.0固定scope為7 datasets × chr1–22＝154 tasks，依`linkage_basis × threshold`及dataset分層聚合call funnel、k route、solver completeness、T/V、likelihood、BQ/error-grid、bootstrap與topology；每個child receipt綁定full extraction child receipt、參數、producer與output hashes。
- **Fail-closed gate**：現有`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/receipt.json`是舊schema 1.0.0 pooled-only pilot，不能進primary ranking；需用schema 1.1.0 extractor重跑HCC1954 chr22 v2後才可產生pilot ranking數字。

### 2026-07-16 07:38 — 獨立 red-team 找到並修正 V>1 runtime blocker

- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py` 的真實 `distinct_vertex_sets_V>1` 分支。
- **Finding**：原實作在 likelihood fit 前，以尚未賦值的 local `all_converged/all_monotone` 作條件，V>1 component 會觸發 `UnboundLocalError`；原 53 tests 只走到 V=1，未覆蓋這條分支。
- **Fix**：V>1 時先無條件 fit 全部 distinct vertex sets，之後才由 `_rank_fits()` 回傳 convergence/monotonicity並 fail-closed。新增 `AA`-only 反例：兩條 minimum paths `RR→AR→AA` 與 `RR→RA→AA` 必須形成 V=2 likelihood tie。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_m2_patterns_and_rank.py -v`；另跑 `git diff --check`。
- **Output**：terminal receipt；ranker SHA256 `358f3985135e6dbeb660e59145313e2b8647fc900adeaafd37c5755650e7092c`，test SHA256 `ebe95f7cc5e0c93892e7cd44091eff435d87b5d339f9c04616d8acf50e4de490`。
- **Actual**：exit=0，13/13 tests PASS；新增 V>1 regression PASS；`git diff --check` exit=0。正式 pilot 仍必須再走真實 schema 1.1 data，不能只靠單元測試放行。

### 2026-07-16 07:45 — schema 1.1 M2 因跨 phase-set 混合暫停 final，改為 PS-aware primary

- **Status**：P0 accepted；M2 full verdict由GO暫回PROBE，M0／symbolic solver不受影響。
- **Finding**：extractor已保存每條 molecule 的`phase_set`，但 schema 1.1 的 HP1／HP2 cut support、component與ranking只依HP family聚合。不同PS block的HP1/HP2 orientation可任意翻轉，故不能把不同PS的同號HP當成同一條haplotype lineage。
- **Read-only audit input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/HCC1954.chr22.molecule_sparse_calls.tsv.gz`及component membership。
- **Actual diagnostic**：15,671 molecules；HP1有46個PS類別（含missing）、HP2有44個；HP1固定R/A sites中33/167涉及>1 PS，HP2為26/207。舊pooled threshold=3的164 components中，42個被>1個HP1/2 PS觸及；進一步限制k>1 component×HP unit時仍有實際跨PS混合（threshold 1/2/3/5皆非全0）。
- **Decision**：primary unit升級為`HP family × known PS × read-linked component × threshold`；不同PS不得共享bridge support或進同一candidate tree。PS missing另列diagnostic/sensitivity，不與known PS合併。正在執行的schema 1.1 chr22 v2允許完成以量化差異，但標記`DIAGNOSTIC_NOT_FINAL`；正式154-task full不得以它放行。
- **Additional red-team fixes**：V>1 primary nonconvergence時不得執行bootstrap；唯一vertex-set摘要必包含`T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED`並做partition守恆。targeted tests 15/15 PASS。
- **Uncertainty contract**：現有molecule bootstrap固定candidate sets，只能稱`conditional-on-fixed-candidate-set ranking bootstrap`，不可稱整條建樹流程穩定性。BQ primary需採明示的conditional biallelic substitution假設；由同reads衍生VAF仍禁止重複加分。

### 2026-07-16 08:01 — 終止已被取代的 schema 1.1 chr22 v2

- **Input**：進行中的`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22_v2/`。
- **Command**：對本研究自己的 unified exec session送出`SIGINT`；未碰觸使用者的all-sSNV processes。
- **Output**：同目錄`ABORTED.md`；保留不完整artifact，不刪除。
- **Actual**：exit=130、wall=38m21.64s、peak RSS=93,436 KB；中途gzip明示缺EOF marker、無receipt，故fail-closed。停止理由不是計算失敗，而是schema 1.1已被PS-block P0新證據淘汰，繼續掃BAM只會競爭I/O且產生不可用結果。

### 2026-07-16 08:19 — cycle state 與 PS amendment 對齊

- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md` v0.3 與 `InterSubMod/state/cycles/cycle_20260716-0618-read-linked-hypercube-exact-likelihood/audit.json`。
- **Command**：以 JSON parser 驗證 state audit，將已被 amendment 取代的 `verdict_GO` 改為 `verdict_PROBE_M2_PS_contract`；保留原始 90/100 分數與 initial method probe 作歷史證據。
- **Output**：`InterSubMod/state/cycles/cycle_20260716-0618-read-linked-hypercube-exact-likelihood/audit.json` v0.3。
- **Actual**：目前 overall probe 為 `PENDING_PS_AWARE_M2_PILOT`；symbolic／MILP initial probe 仍為 PASS，但不得再讓舊 GO 狀態誤放行 schema 1.1 full M2。只有 schema 1.2 extraction＋2.0 ranking pilot 完成 PS 隔離、hash 與守恆 QA 後才能重新評分。

### 2026-07-16 08:35 — PS-aware M2 schema 1.2／ranking 2.0 contract 完成

- **服務目標**：G3（read-level mutation evidence）、G4（7 technical datasets／6 biological samples 一致口徑）、G5（可重算與fail-closed證據鏈）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{extract_lossless_read_linkage.py,run_full_m2_extraction.py,build_m2_patterns_and_rank.py,run_full_m2_ranking.py}`、tiny synthetic BAM／VCF／PS sidecar與合成 sparse-call/component fixtures；未讀寫 canonical 結果、未重跑 BAM pilot/full。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Output**：上述4支 scripts、`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/{test_extract_lossless_read_linkage.py,test_full_m2_extraction.py,test_build_m2_patterns_and_rank.py,test_full_m2_ranking.py}`與本 living note。
- **Actual**：exit=0；84/84 tests PASS。extractor／full-extraction schema固定`1.2.0`；child／full-ranking固定`2.0.0`。producer SHA256依序為`d30e15cd175f7b812a062a45962b3ebf18ceccee0b52657c697899522e471a3a`、`fc5af8383880013207845487d7b1c5d91a4759a59a29953db326681a8224a30b`、`fa651efb6e3b0e38942d2611fbd422713a1b774a5516fa69c1ea722756649a9c`、`ff03ef3087ff829ecc5d665377e26be6601bc6b2b00e707ced8826714a24eedf`。

#### Schema 1.2 extraction contract

- primary unit固定為`HP family × exact known PS × read-linked component × bridge threshold`；不同PS永不共享support或membership。相同site可出現在多個PS units，key不衝突。
- primary component只含該`HP×PS`至少有一筆fixed R/A的active sites；missing PS輸出到`MISSING_PS_HP1/2` sensitivity-abstain，舊`pooled/HP1/HP2/...`明標diagnostic non-primary。
- receipt量化known/missing PS molecule/site memberships、舊跨PS或missing-PS聚合才成立的boundaries、legacy component跨多PS數與PS-aware component分裂量；每threshold 1/2/3/5分開。

#### Schema 2.0 ranking contract

- primary BQ model不把`e=10^(-Q/10)`直接當R↔A flip；在明示symmetric substitution後，condition on observed∈{R,A}：`match=(1-e)/(1-2e/3)`、`flip=(e/3)/(1-2e/3)`。O/D/S/L/X仍作missing，保留claim ceiling；fixed grid明示為conditional binary-flip sensitivity。
- structural gate是「同一個exact R/A/X pattern count ≥ minread」，不把`AXR/ARX/AXX`相容pattern默默pool；full contract固定minread grid 1/2/3/5，primary=3。
- exact gate使用minread-specific `k_observed_alt_active`，另列`k_component_sites`、`k_scoring_alt_observed`與`n_not_structural_alt_active_sites`；raw component k>12但effective k≤12仍可exact求解。
- partial read不展開：每個pattern列conceptual `2^u`與structural observed-alt universe下的effective completions；receipt分開unique RAX groups、BQ-quality groups、molecule projections三個分母，並列full/partial、u分布、coverage denominator／covered／unsatisfied（要求0）。
- bootstrap正式改名`conditional_candidate_ranking_bootstrap`，只對固定candidate set重抽molecule ranking，不宣稱整條建樹穩定性。
- topology正式aggregate只納入likelihood reliable units；unique class互斥計數、ambiguous class-set計數與inclusion counts分開，nonconvergence不滲入正式topology分母。
- 新增compressed candidate table與primary-winning pattern→state BQ-weighted posterior responsibility；responsibility明標latent expectation，不能解讀為read/cell/clone硬指派。full runner再產生report-gated canonical candidate table。
- resume reuse除receipt sidecar／semantic hash外，重新hash ranker實際消費的calls/sites/components/membership，並驗`hp_families`與`component_basis_mode`；upstream TSV被改但receipt未改時仍拒絕reuse。
- full aggregate逐dataset×PS-aware basis×threshold與global驗selection、solver、projection、structural/scoring、BQ、T≥V、topology、partial coverage等守恆；任何cell失敗則`all_pass=false`。scope明示7 technical datasets／6 biological samples，HCC1395與HCC1395_DORADO為同biological sample的technical datasets。

#### 尚未完成／不可宣稱

- 本次只有contract與synthetic regression驗證，沒有產生154個dataset×chr的實際數字；full extraction／ranking仍須等resource gate放行後由主流程執行。
- 因此目前不可把schema 1.1 pilot數字升級成PS-aware final，也不可宣稱真實clone數或唯一parent edge。

### 2026-07-16 08:37 — k>12 exact efficiency 靜態稽核；正式 benchmark 資源 gate

- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`（SHA256 `db85c94f...`）、既有 pilot receipt（SHA256 `9939a9ae...`）、Mahapatra et al. directed-hypercube Steiner 原始論文。
- **Command**：`python3 -m py_compile .../scripts/benchmark_k_gt12_exact.py .../tests/test_benchmark_k_gt12_exact.py`；`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m unittest .../tests/test_benchmark_k_gt12_exact.py -v`。另跑 k=3 三路 subprocess/hard-timeout smoke；正式 k13–16 命令已寫入專題文件但**未執行**。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_k大於12_exact效率與certified策略_01.md`；`.../scripts/benchmark_k_gt12_exact.py`（SHA256 `bdfcfda2...`）；`.../tests/test_benchmark_k_gt12_exact.py`（SHA256 `b55db6d9...`）。正式 `results/k_gt12_efficiency_pilot/` 未建立。
- **Actual**：4/4 unit tests PASS；k=3 smoke 3/3 solve gap=0、3/3 enumeration complete，wall 2.779 s（`/tmp`，非研究證據）。靜態稽核確認現行 variables=`2^active-k`、partial row有`2^u` coefficients、每個no-good row有`2^active-k` coefficients；提出 exact active/model reduction、group dominance、sparse row generation、separator boundary enumeration、ordinary-terminal DP、exact BDD/ZDD/branch-price與 h* 雙 certificate。
- **Resource gate**：08:37 有 48 logical CPUs、46個相關使用者 workers、約3,162% CPU、load=61.11/58.55/151.33；依不競爭政策以 `RESOURCE-GATED` 結案，沒有啟動k13–16，也沒有捏造 objective/gap/runtime。
- **Exactness amendment**：observed-alt compression只在該 structural contract內稱exact；相對all-loci只主張 objective／投影最優解存在，不主張保留所有full-universe V-set identities。正式script已預留同input的k13 all-loci sensitivity。

### 2026-07-16 08:55 — 主 agent 獨立 code review 與 86-test 全套回歸

- **服務目標**：G3（read likelihood 不重複計算）、G4（全資料集使用相同 fail-closed 規則）、G5（報告數字可由 receipt 追溯）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/{scripts,tests}` 與 `InterSubMod/state/cycles/cycle_20260716-0618-read-linked-hypercube-exact-likelihood/audit.json`。
- **Command 1**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Command 2**：`rg --files InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests | rg '\\.py$' | xargs -r python3 -m py_compile`。
- **Command 3**：`python3 -m json.tool InterSubMod/state/cycles/cycle_20260716-0618-read-linked-hypercube-exact-likelihood/audit.json`；另以 `rg -n '[ \\t]+$' ... -g '*.py'` 檢查 Python 行尾空白。
- **Output**：terminal receipts；更新後的 ranker、full-ranking runner、HTML builder 與兩個 regression tests。
- **Actual**：全套 86/86 tests PASS、exit=0；所有 Python `py_compile` PASS；cycle JSON parser PASS；Python trailing-whitespace 0 matches。主 agent 額外修正三點：
  1. ranker文件明確以conditional BQ likelihood為primary、effective active k為solver gate、fixed-error grid僅為sensitivity；
  2. 若同一unit任一candidate fit nonconverged／nonmonotone，candidate table的整個unit一律標`ABSTAIN_UNIT_OPTIMIZER`，禁止只排除單一壞candidate後宣稱winner；
  3. HTML漏斗將舊欄名`hp_included_molecules`的畫面語意修正為「進入至少一個known-PS primary unit」，不再誤稱所有HP-tagged molecules都已進入primary。
- **Producer SHA256**：`build_m2_patterns_and_rank.py=f20bbd6583f579816c2b4b7d7ab00b91d5b5997c7698faa7774ba3e0da947424`；`run_full_m2_ranking.py=a00ca1f3b65df3a49c6f90bef800c91ac4375ca96a0efb663cfb820ea2f47788`；`build_validated_html_report.py=1581e8aee9a41091b9dc34f8a70df887cb9ff809c182fd4349156866b2378515`。
- **Remaining gate**：M0正在跑最後資料集；M2 BAM pilot、154-task full與k=13–16 benchmark仍因既有all-sSNV workload超過48 logical CPUs而resource-gated，未以synthetic PASS冒充全量數字。

### 2026-07-16 09:01 — M0 七資料集 census 完成並由獨立 verifier 深度核對

- **服務目標**：G4（7 technical datasets／6 biological samples 全範圍）、G5（逐 byte、逐 unit、候選集合可追溯）。
- **Input**：canonical v5 的7份`layered_region_view_<dataset>.json`，root=`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`。
- **Command 1**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time -v python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_m0_likelihood_census.py --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4`。
- **Command 2**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 nice -n 10 ionice -c2 -n7 /usr/bin/time -v python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_m0_likelihood_census.py --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4 --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --candidate-mode deep --json-output InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json`。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/{m0_receipt.json,m0_hp_lineage_likelihood_census.tsv.gz}`與`.../results/m0_audit_full_v4/independent_verification.json`。
- **Actual M0**：exit=0、wall=1:50:10、peak RSS=640,028 KB。64,973個eligible HP-lineage units、46,385個至少含一個eligible unit的regions；T=1/V=1 26,225，T>1/V=1 edge-unresolved 101，T>1/V>1 38,647。V>1中likelihood unique vertex set 8,751、likelihood tied 27,270、optimizer abstain 2,626；partition sum=64,973。raw T=444,007、distinct V=443,745。
- **Actual independent verification**：verdict=`PASS`、exit=0、wall=4:54.97、peak RSS=1,259,236 KB；64,973/64,973 canonical eligible units深查，stored 61,702、由frozen solver重建3,271；missing=0、extra=0；eligible 64,973 + capped 8,021 = primary mutation HP units 72,994；所有scope、hash、aggregate、T/V及selection partitions通過。
- **SHA256**：M0 receipt=`eba081a70f16c008f70cd97c85ec3bcbce41d3982eb6a55c5915e89149197699`；gzip rows=`9df74db30bc930fee4e0a6941f371bfe12b069808423bb53be5eaf3fc77c1a6c`；independent verification=`912760104cb3b7dca4e18d56ac429f0cbc901c81800719425ebaf228c2ae735c`。
- **Claim ceiling**：M0使用canonical v5 thresholded alignment exposure與固定error rate，只能作既有候選歧義的engineering baseline；它不是PS-aware lossless-molecule final、不是校準後model selection、不可解讀為真clone數或唯一parent edge。正式生物結論仍以M2 schema1.2／ranking2.0全量結果為gate。

### 2026-07-16 09:08 — M0 optimizer 假失敗診斷與 global-KKT robust contract pilot

- **服務目標**：G3（read-pattern likelihood 數值可信）、G4（全7 datasets同一停止準則）、G5（全域誤差上界可驗證）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/{m0_receipt.json,m0_hp_lineage_likelihood_census.tsv.gz}`、canonical v5 frozen candidate solver與前10個`RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE` units；不讀BAM、不加VAF。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 nice -n 15 ionice -c3 python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/diagnose_m0_optimizer_nonconvergence.py --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --m0-output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_optimizer_audit_pilot10_v1 --max-units 10`。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_optimizer_audit_pilot10_v1/`；unit TSV SHA256=`957d33c9e513b037856051d8dcb4f7a7f52a0c478b8668c78923364284d21842`，candidate TSV SHA256=`db21aafa911efcd0dfa212e603e22f5554b9b4350a736a9411086ab4c5c40bd6`。
- **Actual**：10/10舊abstain units、690/690 candidate fits全數取得global LL gap≤`1e-8`證明；31 fits的舊SLSQP status=8，但新contract全救回。warm-start gap max=`5.616e-4`→final max=`9.962e-9`；LL調整絕對值max=`7.924e-8`；pairwise refinement iterations p50=1、p99=14、max=19。舊/新top vertex-set IDs 0/10改變，10個counterfactual皆為likelihood tied。
- **Method decision**：負log-likelihood是simplex上convex（等價正log-likelihood concave）；SLSQP只當uniform warm start，後續用deterministic pairwise mass-transfer + exact concave line search單調精煉，並以Frank–Wolfe/KKT gap作`LL(optimum)-LL(current)`全域上界。`scipy.success` 不再是收旂權威。
- **Claim ceiling**：只證明已宣告read-pattern mixture的數值全域最佳；若`rank([Q;1])<V`，則π不可唯一識別。即使可識別，π仍是latent mutation-state exposure proportion，不是cellular clone fraction。
- **Tests**：`test_hypercube_exact.py + test_build_m2_patterns_and_rank.py` 40/40 PASS，包含boundary optimum、rank-deficient π、V=125 deterministic stress、SLSQP status=8 rescue、EM LL cross-check、simplex/KKT證明與BQ likelihood。
- **Remaining gate**：M0 2,626個舊abstain units的full re-fit預估單核約40分鐘；09:08 load=`109.01/87.53/86.08`且仍有47個all-sSNV相關process，依資源政策等使用者jobs清空後才執行，不以10-unit pilot充當full。

### 2026-07-16 09:12 — PS-aware 稀疏化、串流整合與主 agent 99-test 回歸

- **服務目標**：G3（exact-PS read evidence不混塊）、G4（154 tasks相同聚合契約）、G5（效能修正不改科學語意且可回歸）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{extract_lossless_read_linkage.py,build_m2_patterns_and_rank.py,run_full_m2_extraction.py,run_full_m2_ranking.py,hypercube_exact.py}`與同目錄全部tests；未掃描BAM、未改canonical資料。
- **Command 1**：`/usr/bin/time -v env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Command 2**：以`find ... -name '*.py' -print0 | xargs -0 python3 -m py_compile`、`python3 -m json.tool`及`rg -n '[ \\t]+$' ... -g '*.py'`檢查編譯、cycle JSON與行尾空白。
- **Output**：更新後的4個M2 scripts、robust optimizer、99項regression tests與terminal receipts。
- **Actual**：主agent獨立重跑99/99 PASS，exit=0、wall=3.68s、peak RSS=84,388 KB；py_compile PASS、cycle JSON PASS、Python trailing whitespace=0。
- **效能修正**：extractor的PS bridge／legacy-union改為稀疏two-pointer support，不建立`#PS × chromosome-sites`矩陣；ranker建立`HP family × exact PS`route index，不再讓每個molecule掃描全染色體所有PS/component；full extraction的`max_k`改取task-wise最大值而非相加；candidate／responsibility detail改逐unit gzip串流及增量semantic hash。
- **語意回歸**：新增old-vs-sparse equivalence、PS isolation、route lookup complexity counter、max-vs-sum、streaming digest與sink API tests。穩定版100,000列microbenchmark中，僅就明細list/hash層，materialized peak 35,175,566 bytes、streamed peak 5,392 bytes、digest一致；此數字不代表MILP整體RSS。（開發中的另一個不同row fixture曾量到47,575,566／5,502 bytes，因輸入物件不同，不混入穩定版receipt。）
- **Producer SHA256**：extractor=`84fb041a965acd667e08a10e6aaaf5951d6de929657a5b05c81d531c2701f222`；ranker=`b7a19439b037a36f25d74c5bce337fbdeec1e303a95d8bbcc460a4874c25f03c`；full extraction=`465668a14b89e7b5f01acb84d7986ac84435057eab840670871981cbfb899005`；full ranking=`1a760045123aecdb63adff71ac7ae55fbb03054a0cc5993efca07e7ef4d59505`；optimizer core=`6576910b52cd8a3463b9bbe3d22fde5d735f5af78ae2c512f0bf4b67031d5fb8`。
- **Remaining gate**：串流解除全域detail list的RAM風險，但154-task磁碟量仍須由PS-aware單染色體pilot量測四個gzip的rows/bytes、wall/RSS後才能放行full。

### 2026-07-16 09:13 — 全量資源 preflight 仍 fail-closed

- **Input**：canonical v5 manifest=`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json`（SHA256 `16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45`）；預定output=`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_full_schema_1_2`。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_extraction.py --manifest <manifest> --output-root <output> --workers 2 --samtools-threads 1 --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5 --preflight-only`。
- **Output**：stdout JSON preflight receipt；因gate失敗，正式output root未啟動。
- **Actual**：exit=2、`resource_gate_pass=false`；scope正確列出7 technical datasets／6 biological samples、chr1–22、154 tasks。偵測46個conflict processes、2個root：HCC1954 all-sSNV group 4 processes／約187.1% CPU，以及first-six all-sSNV group 42 processes／約3,192.6% CPU。未使用`--ignore-resource-gate`，也未終止使用者工作。

### 2026-07-16 09:14 — PARTIAL 教授版 HTML 重新生成與真實瀏覽器 QA

- **Input**：current canonical topology JSON、M0 full receipt、symbolic pilot receipt、partial-method／literature audits與舊HCC1954 chr22 diagnostic pilot；full M2 extraction/ranking receipts刻意缺省。
- **Command 1**：`python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py --canonical-json InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json --m0-receipt .../results/m0_full_v4/m0_receipt.json --pilot-receipt .../results/pilot/pilot_receipt.json --method-audit .../partial-read-method-audit.md --literature-audit .../20260716_方法學原始文獻與主張邊界_01.md --m2-pilot-extraction-receipt .../results/m2_pilot_HCC1954_chr22/receipt.json --output .../20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.html --allow-partial --overwrite`。
- **Command 2**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/qa_validated_html_report.py --html .../20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.partial-preview.html --output-dir .../results/html_qa_partial_full_m0_v5 --expect-status partial`。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.partial-preview.html`（62,541 bytes；SHA256 `7275573cd1a2e23baf69dfa9aa8c5cc143437a4e20f21954c0bae18546710e85`）與`.../results/html_qa_partial_full_m0_v5/{browser_qa_receipt.json,desktop_full.png,mobile_full.png,print_qa.pdf}`。
- **Actual**：generator exit=0但`all_pass=false/final_ready=false`，明列缺M2 full extraction與ranking receipts；不會建立final-looking HTML。headless Chromium QA `all_pass=true`：11 sections、5 SVG、9 details；桌機1440px與手機390px均無overflow；anchors/local links全通過；remote requests、console/page errors均為0。主agent另實際檢視desktop/mobile full-page screenshots；版面與PARTIAL標示清楚。

### 2026-07-16 09:58 — Exact solver reductions、source-backed funnel 與 FINAL fail-closed gates

- **服務目標**：G3（partial-read結構語意）、G4（7 datasets守恆與獨立重算）、G5（教授版數字可追溯且不能由fixture假通過）。
- **Solver decision**：保留「每個達 structural exact-pattern minread 的 distinct partial pattern＝一個joint group constraint」語意；低於門檻的informative molecules只進likelihood。加入duplicate／dominance、mandatory-hit、singleton forcing、active-bit predecessor、fixed-hidden sparse no-good與`h=0` early complete。每次MILP construction的sparse row可列出同一reduced group的`2^u_eff`個相容active vertices，但不把它們當成分開tree worlds，也不採first-success。
- **Solver verification**：`test_hypercube_exact.py` 21/21 PASS；k=1..3共18 seeded cases與獨立exhaustive oracle的objective與完整optimal vertex-set集合全同；k=4仍列出24 optima且digest前後相同。Solver SHA256=`891c0e469c8b24f37a4e98668a564f46c2cd505ed3da633f9bdd63b6f7294aa7`；audit=`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_Hypercube_exact_solver縮減與列舉效能稽核_01.md`。
- **Current funnel input**：`InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`與canonical v5的7份`ssnv_site_ledger_*.summary.json`。
- **Current funnel command**：`python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_current_funnel_receipt.py --canonical-json InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json --ledger-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --output InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/current_funnel_receipt.json`。
- **Current funnel actual**：exit=0、`all_pass=true`；638,259 raw ClairS → 582,820 LongPhase-S PASS → 469,849 chr1–22 biallelic → 194,149 retained。分支為55,439 LongPhase-filter excluded、112,971 non-autosomal/non-biallelic、225,268 MAX_SNV-excluded sSNV records、50,432 positional singleton、194,149 retained；raw與autosomal兩層皆守恆。Receipt SHA256=`a6282ecb73df314782d39a6ae4df410cebca05a42425bd7eebe28f12b1d35d75`。
- **Professor-report red-team**：獨立audit指出6個P0：M0 verifier未設gate、HP1/HP2任選、T/V誤印比例、漏斗不完整、per-dataset M2空cell、extraction可接受154 PASS+額外SKIP。已逐一加入builder與15項針對性tests；另明示legacy `C_region`、7,975 candidate-incomplete與MAX_SNV sSNV branch不可混稱、legacy schema1.0 S8只作diagnostic。
- **Independent evidence gates**：FINAL現在額外要求source-backed funnel receipt、M0 row TSV＋deep independent verifier，以及M2 full independent verifier。M2 verifier不匯入production aggregator/ranker、不讀BAM，會從154+154 child receipts重算aggregates並逐列重建candidate table；implementation tests 8/8 PASS，真實full154 receipt仍PENDING。
- **Full regression command**：`/usr/bin/time -v env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py'`。
- **Full regression actual**：121/121 PASS，exit=0、wall=5.42s、peak RSS=86,216 KB；全scripts/tests `py_compile` PASS。此結果只證明程式契約；真實M2 154-task數字仍未完成。
- **Resource pilot plan**：先跑`HCC1395_DORADO×chr6`（最大VCF 27,127）壓extraction/PS routing，再跑`H2009×chr2`（舊代理2,133 family units、919 T>1、309 capped）壓solver tail；兩者未過之前不啟動full154。文件：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2_schema1_2_2_0安全資源pilot規劃_01.md`。

### 2026-07-16 10:16 — M2 per-unit solver/runtime 計時與獨立重算契約

- **服務目標**：G4（全154 tasks使用同一計時口徑）、G5（可從per-unit表重算p50/p95/p99）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{build_m2_patterns_and_rank.py,run_full_m2_ranking.py,verify_full_m2_receipts.py}`與對應tests；未讀BAM、未啟動pilot/full M2。
- **計時契約**：`time.monotonic_ns()`在每個unit分開量渪candidate generation、likelihood fitting、unit total；另存candidate/likelihood是否實際invoked。數值換算為seconds，只作當下機器/load的performance diagnostic，不作跨機可重現數據。
- **Output contract**：child新增`m2_unit_runtime_diagnostics.tsv.gz`；receipt對primary/all-minread與primary-invoked segments保存`n,sum,max,p50,p95,p99`。quantile是exact empirical nearest-rank：一起算`ceil(p*n)`，`p=0.50/0.95/0.99`。
- **Full aggregate**：154個child TSV逐列讀取，不materialize rows；只保留5個packed float64 vectors，上限40 bytes/primary unit。同時輸出all-primary與invoked-only tail，避免not-run units的0秒稀釋solver p95/p99。
- **Independent verification**：`verify_full_m2_receipts.py`不匯入production ranker/aggregator，使用full sort從per-unit TSV獨立重算child primary/all-minread/invoked summaries與full aggregate；TSV/header/hash/summary/invocation flag篡改均fail-closed。`SUM_FIELDS`無需修改，因runtime是獨立diagnostic contract，不是可加的科學count。
- **Semantic isolation**：runtime欄位只存獨立TSV與underscore internal payload，不進`unit_semantic_sha256`；同一科學row在不同wall time下semantic hash仍完全相同。
- **Partial-read wording correction**：修正成「不materialize completion-wise tree worlds/cross-read Cartesian products；但每個reduced group的sparse MILP hit row會materialize active compatible vertex indices（最多`2^u_eff`）」，不再誤寫成「任何completion都未materialize」。
- **Command**：`/usr/bin/time -v env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Actual**：exit=0；130/130 tests PASS；wall=`7.12s`；maximum RSS=`86,848 KB`。全scripts/tests `py_compile` PASS，本次touched Python files行尾空白=0。
- **Source SHA256**：ranker=`476539baac01e4b82a73c414288e99be557cef52baef32c9dee73be6f2dbdab3`；full ranking=`da1896e50ab8b0aa65f38eaab763b71a36e19f26c2ebc400148fbdabe4095d94`；independent verifier=`7259d3f32ceb27a9bfe151153127a249c0d76f40671d0f9cbdafef9c88d964ca`。
- **Remaining gate**：這只解除「無正式p95/p99口徑」的planner blocker；真實數值仍必須從PS-aware single-chromosome pilot產生，且通過資源/磁碟gate後才能啟動full154。

### 2026-07-16 10:20 — Current-solver identity-bound pilot、漏斗獨立驗證與資源再閘門

- **服務目標**：G3（partial-read solver語意精確）、G4（7 datasets漏斗獨立守恆）、G5（程式／receipt身分可重現）。
- **Solver wording-only update**：`hypercube_exact.py`的module contract修正為「每個reduced group的sparse MILP row會列出active compatible vertex indices，但不建立completion-wise tree worlds或cross-read Cartesian product」；演算法邏輯未改。Current solver SHA256=`9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`。
- **Pilot input**：canonical v5 root=`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`；manifest SHA256=`16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45`。
- **Pilot command**：`/usr/bin/time -v env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_exact_symbolic_pilot.py --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot_v3_identity_bound`。
- **Pilot output／actual**：`results/pilot_v3_identity_bound/`；exit=0、wall=8.60s、peak RSS=107,948KB、`all_pass=true`。Symbolic/explicit為9,840 patterns、2,015,538 checks、0 mismatch；legacy/MILP為66個vertex-set checks、0 mismatch；k=9–12為12/12 PASS；likelihood controls PASS。Receipt SHA256=`5cab8c37f31908d13b0c921619264598a4e6e3ef29cbcf723d1a6c421b6297c8`，內綁runner SHA256=`816086246c4102f419216685ee76996aa2e859dedf16fe2d5e55583f29c5f865`與current solver SHA256。
- **Independent funnel actual**：獨立standard-library verifier未import producer，從canonical＋7份site ledgers重算；322/322 checks PASS、0 failures，確認`638,259 → 582,820 → 469,849 → 194,149`與`469,849 = 50,432 + 225,268 + 194,149`。Receipt SHA256=`a0a098f103980204269d92a1d75ac148b408be3c0ea408349a1a84f7487eb796`。
- **HTML QA v6**：輸入`20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.partial-preview.html`，輸出`results/html_qa_partial_full_m0_v6/`；browser receipt `all_pass=true`，11 sections、5 SVG、11 details，desktop 1440px／mobile 390px均無overflow，PARTIAL狀態吻合；HTML SHA256=`eae8a4285870dfcdf5cb23587a0760221dce6e860bda206642e34ea51df47187`。此版仍待report red-team修補後重建，不作final。
- **Resource preflight再驗證**：同schema1.2 preflight於10:19 exit=2，`resource_gate_pass=false`；scope仍為154 tasks，但既有all-sSNV作業為2 roots／57 conflict processes。未使用`--ignore-resource-gate`，未啟動DORADO chr6、H2009 chr2或full154。

### 2026-07-16 10:30 — 教授版語意gate收斂、135-test回歸與HTML v7

- **Builder fixes**：FINAL新增S14 current-funnel independent verification hard gate；candidate rows／units／mean rows-per-unit與parent choices改為count/ratio，禁止誤印百分比；legacy S8分母明示raw HP value含missing `.`；canonical逐dataset新增`complete+incomplete=W_primary`、`W_primary+no_primary=W_tree`、mutually-exclusive topology sum=`W_primary`，可抓兩dataset `+1/-1` aggregate cancellation。
- **Pilot identity gate**：builder會重算`pilot_v3_identity_bound`內runner與current research solver的live SHA，舊receipt或事後改檔皆fail closed。
- **Command 1**：`/usr/bin/time -v env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Actual 1**：135/135 PASS、exit=0、wall=8.48s、peak RSS=87,188KB；全部Python `py_compile` PASS、Python trailing whitespace=0、cycle JSON可解析。
- **Command 2**：`python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py --canonical-json .../current_layered_topology_v3_raw_all_v1.json --funnel-receipt .../current_funnel_receipt.json --funnel-verification-receipt .../independent_verification.json --m0-receipt .../m0_receipt.json --m0-verification-receipt .../independent_verification.json --pilot-receipt .../pilot_v3_identity_bound/pilot_receipt.json --method-audit .../partial-read-method-audit.md --literature-audit .../20260716_方法學原始文獻與主張邊界_01.md --m2-pilot-extraction-receipt .../m2_pilot_HCC1954_chr22/receipt.json --output .../20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.html --allow-partial --overwrite`。
- **Actual 2**：`generation_pass=true`、`final_ready=false`；唯一final issues為缺M2 full extraction、full ranking、full independent verification。輸出`20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.partial-preview.html`，75,803 bytes，SHA256=`e22c6eeb05dbfe7afa054413c535097fff06ee898618616a9fe3bc7a2682ff9a`。
- **Browser QA**：輸出`results/html_qa_partial_full_m0_v7/`；receipt `all_pass=true`，11 sections、5 SVG、11 details，desktop 1440px／mobile 390px無overflow，local links/anchors全通過，無remote requests、console/page errors。Receipt SHA256=`f619a87154f5b00535b412d634ca0754a88ae1f8653b8a16fe3ca8326c621278`；另由主agent目視desktop/mobile full screenshots通過。
- **Claim status**：此版可作方法與既有漏斗/M0的內部教授預覽，但明示`PARTIAL PREVIEW · NOT VALIDATION EVIDENCE`；不得在三份M2 terminal receipts前改名或聲稱final。

### 2026-07-16 10:50 — [決策][偏離][驗證] Partial-read 方法文件與 current implementation 對齊

- **[決策] Structural／scoring 分層**：只有`count >= structural_exact_pattern_minread`的distinct exact R/A/X patterns形成structural full/group constraints；低於門檻的informative molecule projections仍與達門檻資料一起進quality-aware likelihood。門檻是結構限制門檻，不是read-discard filter。
- **[決策] Partial materialization 精確語意**：`u`個X仍代表full cube中`2^u`個conceptual completions；每次MILP build／rebuild、每個reduced active group會materialize至多`2^u_eff`個compatible active indices作為該次sparse row係數，但不建立completion-wise tree worlds或cross-read Cartesian products。source-bound method contract以`active_compatible_vertex_indices_materialized_for_sparse_rows=true`與`completion_wise_tree_worlds_materialized=false`同時表達這兩層。
- **[偏離] 修正舊文件過強／過期主張**：不再寫成「completion完全不materialize」或「整個workflow只materialize一次」；current pilot改綁`results/pilot_v3_identity_bound/pilot_receipt.json`（SHA256=`5cab8c37...`）、current solver=`9dbaf8ec...`、k9–12最大runtime=`0.191109856 s`。`891c0e...`與108/108明標為exact-reduction audit-time historical snapshot。
- **[偏離] 複雜度與VAF證據層**：縮減前dense no-good的`2^m` nonzeros明標baseline；current fixed-hidden sparse row約為`h*` nonzeros，但每個optimum約再build／solve一次的tail仍在。同read VAF不重複加分是method declaration＋score construction／regression契約，不能把hard-coded receipt boolean本身當獨立驗證證據。
- **[驗證] Input／command／output**：輸入為current solver與identity-bound receipt；執行`sha256sum .../scripts/hypercube_exact.py .../results/pilot_v3_identity_bound/pilot_receipt.json`，實際得到`9dbaf8ec...`與`5cab8c37...`。再以`rg`掃描5份current方法文件，未見舊`results/pilot/pilot_receipt.json`、`9939a9...`或`0.225 s`被當current claim；`git diff --check` exit=0。輸出為本cycle的`partial-read-method-audit.md`、k>12稽核、exact-solver稽核、M2 extractor QA、`00_INDEX.md`與本living note；未修改任何數據值或BAM/full-M2結果。
- **[驗證] Source-bound method contract migration**：current ranker=`1b54fe07...`、full runner=`576e8c6f...`、independent verifier=`de3af412...`；method-contract migration audit-time 149/149 focused tests PASS。HTML builder仍由主流程進行教授報告修正，故本entry不釘暫時hash。四處exact-compare同一method contract，並把`active compatible indices materialized`與`completion-wise tree worlds not materialized`分成兩個不矛盾的欄位；未讀BAM、未產生154-task M2結果。

### 2026-07-16 11:30 — [驗證] Profile likelihood 獨立重算與教授版 hard gate

- **服務目標**：G3（read-pattern likelihood不重複使用VAF）、G4（154 dataset×chrom child均可重算）、G5（數值與winner/tie分類有獨立實作證據）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py`、M2 child的`m2_symbolic_pattern_counts.tsv.gz`與`m2_compressed_vertex_set_candidates.tsv.gz`契約，以及教授版builder/tests；此階段仍未讀BAM或產生真實M2全量數字。
- **獨立方法**：verifier不import production ranker／solver／aggregator；按child與primary unit lockstep streaming，直接由R/A/X、fixed BQ、count、candidate states與persisted π重算conditional symmetric-substitution profile LL、simplex residual、Frank–Wolfe/KKT gap、relative weights與winner/tie partition。`n_scoring_molecules`明定為跨unit/threshold的molecule projections加總，不是全域unique molecules。
- **負向測試**：額外LL／類VAF score、BQ、π、state、KKT gap、winner/tie任一篡改都會fail；另含真正不同state sets並列第一的positive case。Full receipt新增4個profile checks，全部共同決定`all_pass`。
- **Professor HTML gate**：FINAL現在要求`profile_likelihood_recomputed_from_patterns_states_pi=true`、154 child、profile units/candidates與canonical candidate table完全相等、3個match flags、154份child summaries及count/numeric/streaming contracts；並展示重算coverage、最大差值與streaming peaks。Builder SHA256=`f1ebb17f66d7856bfd3161ff20562ed98fdd8706874c280784d6ddc113054352`；builder test SHA256=`b84f38131b061bcff510af04b2dba63b1e910dfa7c8c4e109a7fce188651d075`。
- **Verifier identity**：`verify_full_m2_receipts.py` SHA256=`4859598d74486f4eba6e4af6fa2dec2b4c0eb5c4e8ed86feac82483e5a7f32d8`；test SHA256=`ae35510a44235ebeed4ffe6e84b385850cb2a29912146f9b283112820bfafb59`。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /usr/bin/time -v python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py'`。
- **Actual**：加入固定completion-world/joint-group反例前為159/159 PASS；其後新增`RRA+AAA+RAX` global-world comparison、`AX+XA` joint coverage/dedup及deterministic exhaustive-reduction audit，current full suite為162/162 PASS、exit=0、wall=12.96s、peak RSS=90,940KB；`test_hypercube_exact.py` 23/23 PASS，builder targeted 23/23 PASS；`py_compile`與`git diff --check` exit=0。這些只驗證方法與fail-closed contract；真實雙pilot與154-task receipt仍等待resource gate。
- **Presolve/no-good receipt**：執行`audit_hypercube_reductions_exhaustive.py --json-output .../results/hypercube_reduction_exhaustive_v1/receipt.json`；61,340 presolve cases、1,979,356 selected-set predicates、23,909 sparse與21,844 dense no-good pairs均0 mismatch，receipt `all_pass=true`、SHA256=`a60037e7cf527d82abfc1b676ba7405a7b8a66625df73491c57ad8e7bd7b88b8`，checksum PASS。原始membership使用獨立R/A/X字元oracle；scope只到k≤3與指定抽象變數上限，不外推成全樣本生物結果。
- **數值分母交叉稽核**：current funnel verifier 322/322 PASS，另用stdlib Python跨funnel/canonical/M0重算22/22 PASS；正式文件為`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_最新數字分母與守恆獨立稽核_01.md`。確認225,268是sSNV records、7,975是incomplete regions、8,021是capped HP units，M0 unique/tied/abstain為HP-lineage units。

### 2026-07-16 11:31 — PARTIAL HTML v9 與資源閘門再驗證

- **HTML input/command**：以current canonical、current funnel＋獨立verification、M0＋deep verifier、identity-bound symbolic pilot與方法／文獻audits執行`build_validated_html_report.py ... --allow-partial --overwrite`；full M2三份terminal receipts刻意缺省。
- **HTML output/actual**：在method-audit加入explicit-world與presolve/no-good可重跑證據後再次生成；`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.partial-preview.html`，76,500 bytes，SHA256=`61fe52cef92a9e5080db707be04d3e731ea7859a95762f77f5f7b0fdc7d3463c`；`generation_pass=true`、`final_ready=false`，唯一final issues為缺M2 extraction/ranking/independent verification。
- **Browser QA command/output**：`qa_validated_html_report.py --html <v9 HTML> --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/html_qa_partial_full_m0_v9 --expect-status partial`；receipt SHA256=`43adb59d61c17aa0de4b22fa96bf885c6a6ed30c9b6218186a68c10902d91a5c`。19/19 checks PASS：11 sections、5 SVG、11 details，desktop 1440px／mobile 390px無overflow，local links/anchors全通過，remote requests、console/page errors為0，HTML在QA期間identity不變；主agent另目視desktop full screenshot通過（mobile layout已由相同QA量測）。
- **Resource preflight input/command**：manifest SHA256=`16f2ef...`，執行schema1.2 full runner `--workers 1 --samtools-threads 1 --bridge-thresholds 1,2,3,5 --preflight-only`，未掃BAM。
- **Resource actual**：2026-07-16 11:29 CST exit=2、`resource_gate_pass=false`；57 conflict processes／2 roots，CPU group sums約1152.5%與2776.2%。未使用`--ignore-resource-gate`，雙pilot/full154均未啟動。

### 2026-07-16 11:50 — 雙 pilot 啟動前的 source-bound preflight

- **Input**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json`；預定 output root=`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_resource_pilot_HCC1395_DORADO_chr6_schema1_2_rank2_0_v1`。
- **Command**：`python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_extraction.py --manifest <manifest> --output-root <pilot>/preflight_unused --workers 1 --samtools-threads 1 --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5 --preflight-only`；其後執行`df -Pk`、`test ! -e <pilot>`與`sha256sum`。
- **Expected output**：stdout preflight JSON；若resource gate通過才可新建pilot root。此命令本身不得掃BAM或留下pilot目錄。
- **Actual**：exit=2、`resource_gate_pass=false`；scope精確為7 technical datasets／6 biological samples、chr1–22、154 tasks。manifest SHA256=`16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45`、extractor SHA256=`84fb041a965acd667e08a10e6aaaf5951d6de929657a5b05c81d531c2701f222`、ranker SHA256=`1b54fe07f5bc49dd041bd9d6fdba056e59a2b20fb9092b260712bc66e88aa59a`；available=`751,482,108 KiB`，pilot root absent PASS。
- **阻塞證據**：57 conflict processes／2 roots；HCC1954 all-sSNV group為15 processes、約1102.8% CPU，first-six group為42 processes、約2522.5% CPU。未使用`--ignore-resource-gate`、未終止使用者作業、未啟動任何BAM讀取。背景監視只等待兩個root自然退出，退出後必須再跑同一preflight取得exit=0才放行。

### 2026-07-16 13:02 — [驗證] 正式輸入身分與154-task runner fail-closed紅隊收斂

- **服務目標**：G4（154 tasks固定同一資料／程式契約）、G5（PRE／POST與失敗收據可稽核）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{verify_m2_input_identity.py,run_full_m2_extraction.py,run_full_m2_ranking.py}`及三份對應tests；本階段沒有讀取真實BAM payload或建立正式M2 output。
- **輸入身分契約**：production只能接受frozen canonical-v5 manifest與完整canonical-v3 schema；42個直接輸入roles必須為42個不同storage identities，其中7個BAM只驗metadata＋manifest指定的21個固定chunks，35個非BAM roles由同一verified FD重算full SHA-256。每次full identity都重新雜湊並完成logical-path／boundary reopen check；PRE的`validation_evidence_eligible`必須是JSON boolean且authority mode精確相同。
- **Runner契約**：full extraction與ranking都只在154個唯一且完整的`7 datasets × chr1–22` child pairs後產生terminal receipt；partial batches只能產生checkpoint。ranking preflight現在會直接驗embedded 154 child receipts、PASS狀態、schema 1.2.0與scope，`results=[]`或duplicate pair即exit 1。一般FAIL／TIMEOUT均留下exclusive durable marker；TIMEOUT marker保存return code、elapsed與stdout/stderr tail。
- **Command 1**：`/usr/bin/time -v python3 -m unittest -v InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_m2_input_identity.py`。
- **Actual 1**：19/19 PASS；獨立唯讀red-team判定P0=0、P1=0；驗證器SHA256=`ac4cb4107d3b3a012ee222ab2bcfb57879197e17948ab3ebc1f0c6516c56a0d7`。真實canonical PRE尚未執行，故actual-data gate仍為WAIT。
- **Command 2**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time -v python3 -m unittest -v InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_extraction.py InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_ranking.py`。
- **Actual 2**：36/36 PASS、exit=0、wall=2.61s、peak RSS=68,316KB。先前malformed ranking preflight反例已由獨立red-team實跑確認：empty results與duplicate pair皆exit 1且不再宣稱upstream pass。current SHA256：extraction runner=`e2824e2d191b3db9810d40011e445c617b0ef7fc4e9eb883afa155b72858d5af`；ranking runner=`90abe30e690f949d01adbf02c8b273be3b8bdda228b03547e80e00add2776ecc`。
- **已揭露限制**：PRE／POST是逐檔snapshot，不能證明分析期間沒有「改後復原」；SHA sidecar不是可信時間戳或數位簽章。process-group timeout也不能殺死主動`setsid` detach的惡意descendant，但current extractor/ranker無detach路徑。正式run manifest必須綁PRE receipt、source snapshot、命令與啟動時間，分析後再做POST exact comparison。
- **Resource status**：13:02仍有一個使用者first-six all-sSNV root與40個active workers；未使用ignore gate、未終止該作業。M2 formal PRE、雙pilot及154-task full仍等待其自然退出。

### 2026-07-16 13:10 — [驗證] Release contract、immutable code snapshot與227-test全套回歸

- **Input**：正式9支execution/verification scripts、canonical manifest schema、`freeze_m2_release_contract.py`與本研究全部tests；真實PRE尚不存在，故只用明確標為non-evidentiary的synthetic fixtures測試。
- **Release contract**：production CLI只接受canonical-v5 authority與正式PRE；PRE verifier path/SHA必須等於當次freeze所讀到的input verifier。9支程式＋schema以原repo相對路徑複製為10個不同inode regular files，file mode=`0444`、directory mode=`0555`，以`renameat2(RENAME_NOREPLACE)`一次publish；run manifest綁7 technical datasets／6 biological samples／22 autosomes／154 tasks、B20/seed、全部source/copy SHA、runtime與snapshot entrypoints。verify-only另記錄自身path/SHA。
- **Focused command**：`PYTHONDONTWRITEBYTECODE=1 python3 -m unittest -v InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_freeze_m2_release_contract.py`。
- **Focused actual**：17/17 PASS；freezer SHA256=`bfd467306dcf726de9b1cecb01afde63dfe4874bba9b0983e0be9e61062e496e`，test SHA256=`9843d3267e4e082d402fa7d042617494dd0a70a8c70126511469881834c423ce`。負測包含PRE／sidecar／verifier path+SHA、source/snapshot symlink或duplicate inode、permission、byte tamper、post-freeze manifest/PRE drift與overwrite。
- **Full-suite command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /usr/bin/time -v python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`。
- **Full-suite actual**：227/227 PASS、exit=0、wall=17.48s、peak RSS=92,052KB；全部Python `py_compile` PASS、cycle audit JSON parser PASS、Python trailing whitespace=0。獨立release-contract紅隊仍在執行，通過前不建立formal contract root。
- **Output**：僅repo內兩支release-contract source/test與本living note；沒有建立`/big7_disk/.../20260716...`正式PRE、contract、pilot或full結果目錄。

### 2026-07-16 15:02 — [驗證][未決] 154-task session證據鏈主代理重播；checkpoint canonical-path P1 尚未收旂

- **服務目標**：G4（7 technical datasets × chr1–22 共154 tasks的中斷續跑一致性）、G5（release、PRE/POST、session、batch、child receipt可追溯）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{run_full_m2_extraction.py,run_full_m2_ranking.py,verify_full_m2_receipts.py}`與三份對應tests；本次只用synthetic 154-task chains，沒有讀取真實BAM或產生M2生物數字。
- **Command**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m unittest -v InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_extraction.py InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_ranking.py InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_verify_full_m2_receipts.py`。
- **Actual**：主代理獨立重播87/87 PASS，exit=0，wall=15.83s。當時source SHA256為extraction runner=`35c214463ae53eac4d6ae2067fbe65b405633692a423202eb5ed01bc541d7b67`、ranking runner=`8cf224d3e4ea6dc7a52cb241dd56775501db59518cef4b579fc9919db30785e6`、final verifier=`c96751c5cda94fcb7375e0b78a893ebe657eea0d51375c3336f312157291c7c3`；tests為`d7b0cdce...`、`7a786792...`、`91754e62...`。
- **已驗證範圍**：production helper可產生extraction與ranking各154 tasks／11 batches的session→batch→grant→completion→checkpoint/terminal chain；final verifier可重新雜湊child receipts/outputs、驗POST精確等於immutable PRE、驗兩stage parent binding，並在發布前二次force-fresh revalidation。Final receipt與sidecar使用exclusive writer且讀取權限凍結，已存在時不覆寫。
- **[未決] P1／為何當時仍NO-GO**：獨立QA發現runner resume雖驗checkpoint內容與sidecar，但尚未強制「`cumulative_attested_tasks` ↔ 唯一canonical checkpoint filename/path」；改名後重建一致sidecar可能被runner接受，卻會被terminal verifier拒絕，造成後續batch自毀。這是orchestration可重現性P1，不影響partial-read數學契約，但修補前不可啟動formal run。兩runner正加入exact canonical path與renamed/resigned negative regression；上述87 PASS在source改變後必須作廢重跑。
- **Resource gate actual**：15:02 CST尚有一個使用者all-sSNV root (`PID 1654819`) 與40-worker child (`PID 1654821`)。Child處於不可中斷disk wait，`read_bytes=3,626,762,240`，相較14:55仍持續增加，無死鎖證據；available disk=`769,472,671,744` bytes。未使用ignore gate、未終止使用者作業、formal PRE/pilot/full仍WAIT。
- **威脅模型邊界**：這些attestations證明的是持久化protocol與內容一致性，前提是無hostile same-UID actor；不是作業系統process ancestry的密碼學證明，也不是可信時戳。

### 2026-07-16 15:08 — [驗證][決策] Release/session hard gate收旂；程式GO、actual-data仍WAIT

- **P1修補**：extraction與ranking runners現均強制`cumulative_attested_tasks`對應唯一canonical checkpoint filename/path，並重驗receipt integrity、session/batch/previous-chain-head與child output identity。新增「checkpoint改名後重建sidecar」negative regression，兩stage均fail closed。
- **Final fixed source SHA256**：`run_full_m2_extraction.py=fffaf70a4a967d489d2b6864b574ccfea02442fd77e0129faa6bffe04b5cd399`；`run_full_m2_ranking.py=4355a7963c7d772a9d860d81ccb81f648b4771d6a35f1e41a2bf6aef79dc8d92`；`verify_full_m2_receipts.py=43a3988f64cdd81e85b0b8027098debceaf408117e970cb4e06dd46f73eab3a9`。Final test SHA256為`008969df...`、`5711f4fb...`、`91754e62...`；freezer=`9690c23e...`，input verifier=`ac4cb410...`。所有負責agents已停止修改。
- **Focused replays**：runner tests 57/57 PASS；verifier 34/34 PASS；主代理三檔重播91/91 PASS，exit=0，wall=16.154s，peak RSS=80,752KB。Production helpers的extraction與ranking各154-task／11-batch (`8→24→…→152→154`) golden chains均由final verifier接受；中斷、orphan、path/inode不一致等negative inputs均被拒絕。
- **Full-suite command/actual**：以具Playwright的`/bip7_disk/liaoyoyo2001/miniconda3/bin/python3` 3.9.12執行`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /usr/bin/time -v python3 -m unittest discover -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py' -v`；261/261 PASS，exit=0，wall=35.49s，peak RSS=100,000KB。先用`cnvtools` Python 3.11的一次重播因該environment未安裝Playwright，258 tests通過但1個HTML QA module import error，該輪明確作廢，沒有排除測試後宣稱PASS。
- **Static QA**：全部research Python `py_compile` PASS；cycle audit JSON parser PASS；Python trailing whitespace=0；`git diff --check -- <research cycle>` exit=0。Markdown檢查命中的是CommonMark強制換行所需的行尾兩空格，不是程式格式錯誤。
- **Independent QA**：固定8個source/test SHA重播109/109 PASS，exit=0；未解P0=0、P1=0，結論`GO`，限本release-contract／receipt／schema與orchestration範圍。
- **當前gate**：這個`GO`不等於154個真實dataset×chromosome tasks已完成。Formal PRE、immutable contract、DORADO chr6／H2009 chr2雙pilot、full extraction/ranking與POST尚等使用者all-sSNV作業自然結束後再執行；不使用resource-gate override。

### 2026-07-16 16:08 — [驗證][未決] 正式執行協定、教授版P1收旂與PARTIAL瀏覽器QA

- **服務目標**：G3（partial/read-pattern語意不誤讀）、G4（7 datasets×chr1-22數字可重算）、G5（正式receipt→summary→HTML證據鏈）。
- **正式協定**：新增`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式凍結全量執行與驗證協定_01.md`，固定PRE→freeze→frozen verify-only→DORADO chr6/H2009 chr2雙pilot（B0+B20）→extraction/ranking各11 batches（`8→24→…→152→154`）→POST→independent final verifier→numeric summary→presentation snapshot/HTML/QA。所有full/pilot批次開始前保留至少300 GiB；雙pilot實測外推若跨過保留線即停止。
- **Pilot provenance P1 closure**：pilot前後均重跑frozen release verify-only；11個snapshot roles的path/SHA/mode/inode全重驗，child extractor/ranker receipt主producer path/SHA與immutable entrypoint精確相等。在已宣告的non-hostile same-UID threat model下關閉P1；不宣稱可防惡意same-UID執行中swap/restore。
- **Final numeric summary**：新增`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_final_numeric_summary.py`。獨立red-team修正「candidate rows只對V、未把`parent_choice_count`回綁raw T卻宣稱T/V conserved」P1；現已逐dataset×HP×threshold對T/V、winner topology、parent-edge uniqueness、h*、runtime與分母exact對帳。Targeted 11/11 PASS；整合102/102 PASS（wall 18.099s、peak RSS 70,804KB）；四類byte tamper全fail closed。Final script SHA256=`8ecb62a6e42e359a4f7d2209fe8d1225032f89d6327fd96be6e101673a9174cc`。
- **HTML builder P1 closure**：完整驗證/顯示partial conceptual/effective completion四種grain、structural/scoring funnel、solver/effective-k、HP1/HP2、exact topology/parent-choice、runtime、zero-conflict session attestation與比例；bridge≥3明標「professor display convention」。主代理獨立重跑28/28 PASS（wall 12.144s）；`py_compile`、fixture JSON、Python trailing whitespace均PASS。Builder SHA256=`6d03a39f2b057ac288b70eb32cff0d31935e29126930515d650ca75a31535b02`。
- **PARTIAL preview input**：current canonical JSON、current funnel+independent verifier、M0+deep verifier、identity-bound symbolic pilot、partial-method/literature audits；刻意不提供尚未存在的M2 full三件receipts。
- **PARTIAL preview command/output**：執行`build_validated_html_report.py ... --output InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_read_linked_hypercube_exact_likelihood全量驗證報告_02.html --allow-partial`；實際`generation_pass=true`、`final_ready=false`，只列缺M2 extraction/ranking/final verifier。產物`.../全量驗證報告_02.partial-preview.html`，73,369 bytes，SHA256=`91ff7c2c7c39e43029592a49837bcc0eaf5ed3a60664f4accab81bcac54a0e07`。
- **Browser QA command/output**：`qa_validated_html_report.py --html <preview> --output-dir InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/html_qa_partial_full_m0_v10 --expect-status partial`；19/19 checks PASS，11 sections、5 SVG、10 details，desktop 1440px/mobile 390px均無overflow，local links/anchors全通過，remote requests、console/page errors全為0，列印PDF成功。QA receipt SHA256=`b42695a1d417be635fdae72051b3de3186e543d2b7bc6426139570f7a18f00f5`；主代理實際目視desktop/mobile full screenshots通過。
- **Presentation provenance snapshot**：新增`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py`，與exact-11科學release分離。固定3個presentation programs、12個evidence roles的absolute path/size/SHA-256、Python/Playwright/Chromium runtime，staged verify後用`renameat2(RENAME_NOREPLACE)`發布；files `0444`、dirs `0555`。主代理獨立重跑11/11 PASS；`py_compile` PASS；destructive delete call=0。初版失敗清理曾使用`rmtree`，主代理code review依治理規範要求修改；現在失敗staging只會原子封存為`.failed-staging.<UTC>.<pid>`，封存也失敗則留原處，不刪除。Final tool SHA256=`639898a882b5d72f9c58ad95701a5afdfa22c0ad37f1aebc8748320d2efbd625`。
- **Current full regression**：主代理以具Playwright的`/bip7_disk/liaoyoyo2001/miniconda3/bin/python3`重跑整個research test suite；287/287 PASS，exit=0，wall=43.41s，peak RSS=97,884KB。全部scripts/tests Python `py_compile` PASS、cycle audit JSON PASS、Python trailing whitespace=0、`git diff --check` PASS；exact-11十支程式SHA與15:08固定值完全相同。
- **[未決] Actual-data gate**：16:05 CST使用者all-sSNV root/child與40 workers仍在，15:54→16:05 physical `read_bytes`+275,656,704 B、iowait約33%，判定持續運算非死鎖。磁碟可用769,466,286,080 B；resource gate維持NO-GO，沒有ignore/中止使用者作業，formal PRE/pilot/full仍未啟動。

### 2026-07-16 16:28 — [驗證][未決] Partial-group／likelihood核心重播與資源等待

- **服務目標**：G3（partial molecule不展開重複計數）、G4（symbolic solver與likelihood可重跑）、G5（正式數字前先鎖定方法契約）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/{scripts/hypercube_exact.py,scripts/build_m2_patterns_and_rank.py,tests/test_hypercube_exact.py,tests/test_build_m2_patterns_and_rank.py}`。
- **Command**：`/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest -v tests/test_hypercube_exact.py tests/test_build_m2_patterns_and_rank.py`，工作目錄為本research cycle；stdout為唯一輸出，未建立衍生數據檔。
- **Actual**：49/49 PASS、exit=0、wall約1.59s。包含`RRA+AAA+RAX` global completion-world minimum、`AX+XA` joint coverage/dedup、k=1..6 symbolic membership、small-k independent exhaustive oracle、sparse no-good、per-call BQ missing marginalization、simplex/KKT certificate、同vertex set的parent-edge score invariance、k>12 abstain與T/V/topology contract。
- **判讀邊界**：此重播只證明程式與小型oracle契約，不能替代7 technical datasets×chr1-22的154個真實extraction/ranking task。16:27 CST user-owned all-sSNV root PID 1654819與child PID 1654821仍存在；child physical `read_bytes=5,719,195,648`且持續上升，因此resource gate仍為NO-GO，formal root尚未建立。

### 2026-07-16 19:31 — [驗證][決策] Resource/provenance reader-test GO；actual-data 仍 WAIT

- **服務目標**：G4（7 technical datasets × chr1–22 的一致 resume 口徑）、G5（每層 receipt／checkpoint／FINAL HTML 可追溯）。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/{run_full_m2_extraction.py,run_full_m2_ranking.py,verify_m2_single_task_pilot.py,verify_full_m2_receipts.py}`、對應 tests 與 `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式執行Runbook_02.md`。本段只用 synthetic／`/tmp` fixtures，沒有讀 BAM／VCF，也沒有建立 formal RUN_ROOT。
- **修補的 5 組 P1**：(1) extraction terminal `REUSED_FINAL` 深驗；(2) ranking terminal `REUSED_FINAL` 深驗；(3) PRE／POST／direct child／pilot GO receipt 的 exclusive immutable sealing；(4) extraction／ranking completed checkpoint 由已驗 child receipts 重建 aggregate／checks／results 後 deep equality；(5) ranking consolidated candidate gzip 必須由 154 child candidate artifacts 串流逐列 lockstep 重建，不只驗 self-consistent metadata。
- **負向反例**：checkpoint aggregate／checks 改寫後重簽、terminal duplicate-for-missing、candidate 實體列／raw SHA／semantic SHA／terminal／sidecar 同時改寫後重簽，現均被 runner 本身拒絕；main 不會誤印 `REUSED_FINAL`。
- **Command／Actual**：獨立 resource targeted suites `117/117 PASS`（elapsed `31.942 s`，wall `32.88 s`，RSS `76,812 KiB`，exit `0`）；全 research suite `327/327 PASS`（elapsed `161.383 s`，wall `162.52 s`，RSS `111,492 KiB`，exit `0`）；exact-11 `11/11`；Runbook 13 Bash fences `bash -n` exit `0`；CLI `4/4`。測試前後兩 runner／兩 test SHA `4/4` 一致。
- **Final source SHA256**：extraction runner=`cf016b9a046c214bbefb6a4b2509955910710fce73d3186dce27b666d5c40fc4`；ranking runner=`66bb175404c207ef320f213c650bb10c6d5fcf3c84cbc40b8ca25e68604da767`；pilot verifier=`9d15ce2bf15af5cc2c4c690cd7718b131108fd8e3946f6a72da40487b06f1578`。
- **Documents**：Runbook SHA256=`576caae20be596f2ba879a8c5b6f9ecce49d401da41335a705c691b35da9cba1`；independent audit SHA256=`e07066345ff6554e9e9137f8b184c44854417ad96ade7978a7f6d06d8933aaf1`。Verdict：reader-test `GO`，P0=`0`、open P1=`0`、P2=`2`；P2 為後續 hardening，不阻擋 live gate。
- **Actual-data gate**：19:30 CST 使用者 all-sSNV root PID `1654819`、child PID `1654821` 與 40 workers 仍存在；child=`wait_on_buffer`，aggregate worker CPU 約 `2096.1%`，root `read_bytes=5,998,473,216`，判定為持續運算，不是只以 gzip 檔案大小判定停滯。Formal actual-data仍 `WAIT`，未繞過 resource gate，未終止使用者作業。

### 2026-07-17 04:43 — [驗證][決策] Frozen v2 正式 pilot B0 達 8 小時上限；terminal verdict＝NO-GO

- **服務目標**：G3（真實 read-pattern likelihood）、G4（154 tasks 同一停止準則）、G5（不以未封存中間檔冒充全量結果）。
- **Task type／scope**：B Comprehensive validation；正式範圍仍為 chr1–22 × 7 technical datasets，但 full 154 extraction＋154 ranking 必須先通過兩個 frozen single-chromosome pilots。HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 technical datasets。
- **Input**：frozen release root `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/`；第一個真實 pilot 為 `HCC1395_DORADO × chr6`。Extraction 已通過 receipt gate；B0 ranking 讀取其 `extraction/`，使用 frozen `build_m2_patterns_and_rank.py`、thresholds `1,2,3,5`、primary structural minread `3`、`exact-k-max=12`、`max-vertex-sets=256`、solver per-call limit `30 s`、bootstrap `0`。
- **Command**：依 `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式執行Runbook_02.md` 的 Bash fences 串接執行；B0 child command 為 release snapshot ranker加 `--require-existing-empty-output-dir`，完整參數保存於 process command line與後續 receipt。Runbook session=`65173`；ranker PID=`2103507`。
- **Output（terminal failed evidence）**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v2/pilots/HCC1395_DORADO_chr6/ranking_bootstrap0/`。沒有完成 `receipt.json`；四個 `.tsv.gz` 均為 timeout 留下的截斷串流，`gzip -t` 各回傳 `1`／`unexpected end of file`，只可作效能診斷，不具 final 科學數據資格。
- **Terminal time evidence**：B0 於 `2026-07-16 20:32:33 +08:00` 啟動；正式 wall gate 為 `≤4 h=GO`、`>4–8 h=PROBE`、`8 h timeout=停止`。`ranking_bootstrap0.time.txt` 保存 `elapsed=8:00:00`、user=`29,267.25 s`、system=`899.07 s`、CPU=`104%`、max RSS=`510,208 KiB`、exit=`124`。程序已停止，沒有進入B20。
- **Fail-closed verifier**：凍結 single-task verifier沒有讀BAM／VCF、也沒有重啟pilot；因 `ranking_bootstrap0/receipt.json` 缺失而輸出 `release_gate.verdict=NO_GO`、`all_pass=false`。Receipt=`.../HCC1395_DORADO_chr6/pilot_gate_verification_receipt.json`，SHA-256=`dc1381a517ba70b0a3f4a0a90b8dbf143ff4b515758a631e3a4ede7fbe456926`，其外部sidecar由正確工作目錄重驗為`OK`。
- **Progress diagnostic（截斷、非 final 數字）**：此 pilot 每個structural-minread設定有 `7,058` 個 `component basis × threshold × component` units，四個設定合計規劃 `28,232` unit-evaluations。截斷runtime串流可讀 `19,270` rows（`68.255880%`）：minread 1=`7,058/7,058`、minread 2=`7,058/7,058`、primary minread 3=`5,154/7,058`（`73.023519%`）、minread 5=`0/7,058`。另可讀pattern=`433,302`、candidate=`20,233`、responsibility=`73,426` rows；這些分子不是完整輸出且無receipt，禁止用於候選樹、拓撲或生物比例。
- **Independent scaling audit＋主 agent 重算（仍為 intermediate performance evidence）**：minread=1 的 `7,058` units 已完整刷新；其中 `4,852` 個 candidate-generation＋likelihood皆 invoked、`33` 個 candidate-generation invoked但likelihood未 invoked、`2,173` 個未啟動 solver。前述33個長尾units的candidate-generation合計=`14,715.729568 s`，全部candidate-generation合計=`15,304.384699 s`，占比=`96.153683%`；其餘4,852個完整units candidate-generation=`588.655131 s`、likelihood=`338.985172 s`。所有33個長尾都有partial structural patterns；判讀是 partial ambiguity 間接產生大量同分最優 vertex sets，而不是程式明列 completion-wise tree worlds。
- **Code-path cause**：`solver_time_limit_seconds=30` 是每次MILP的上限，不是每unit總上限。現行 exact enumeration先求 `h*`，再逐一加入no-good cut，但每個新vertex set都重新build sparse model／重新呼叫solver，最多保留256 sets；少數unit因此可累積遠高於30秒。四個minread都會執行完整candidate generation、BQ likelihood與fixed-error ranking，只有bootstrap限於primary；timeout時primary minread 3只完成可讀runtime rows的73.0%，minread 5尚未開始。
- **Exact-first remediation order（下一個release候選，不在本次run修改）**：(P0) canonical structural-problem cache，完全相同full/group signatures跨threshold/minread重用已證明結果；(P0) persistent solver＋incremental no-good cuts；(P1) shared emission matrix與`unit-evidence hash × vertex-set × model` fit cache；(P1) B20綁定已驗B0，只補conditional bootstrap；(P1) deterministic unit sharding；另加入真正per-unit deadline＋zero-candidate abstain及atomic unit checkpoint/resume。primary-minread-only完整ranking可保留primary estimand，但會改變sensitivity deliverable，必須另立新release contract。
- **效能判讀**：RAM／I/O不是 terminal gate failure 主因；minread 1 中33個 candidate-generation invoked但likelihood未 invoked的長尾units占該minread candidate-generation時間 `96.153683%`。核心瓶頸是少數partial-heavy問題產生大量同分最優vertex sets後，現行程式逐解重建MILP；不是把每條partial read的每個completion各建一棵樹。Likelihood screening可作P2，但無法先解決96%發生在candidate generation的主瓶頸。
- **Terminal decision**：本 frozen v2 以 **NO-GO** 結案並保留所有失敗證據；沒有啟動B20、H2009 chr2、154-task full extraction/ranking、final numeric summary、presentation snapshot或FINAL教授版HTML。下一步必須在新release實作exact-first修復、重新pre-decision audit／freeze／pilot；不得resume或清空本失敗root。

### 2026-07-17 05:17 — [驗證][決策] 教授版NO-GO HTML交付；official portable QA全通過

- **服務目標**：G3（partial-read語意不誤讀）、G4（技術／生物樣本範圍不混淆）、G5（每個教授版數字可回到receipt）。此HTML是第一pilot的blocked／NO-GO驗證摘要，不是FINAL topology報告。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/html_no_go_report_v1/artifact.json`；有效數據只來自frozen extraction receipt與fail-closed gate；截斷ranking只標`diagnostic_only`。
- **Command**：`node /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/deliver_portable_artifact.mjs --input InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/html_no_go_report_v1/artifact.json --output InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2_pilot_NO_GO驗證報告_01.html`。
- **Output**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2_pilot_NO_GO驗證報告_01.html`；成功receipt保存於`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/html_no_go_report_v1/delivery_receipt.json`。
- **Actual**：exit=0；`validation=passed`、`package=passed`、`verification=passed`；source dialog與`keyboard_menu_semantic_click`通過；1440與390兩個viewport通過；rendered counts=`2 blocks／1 chart`。HTML明示Extraction PASS、Ranking TIMEOUT／NO-GO、T／Topo不可用、partial read joint-group constraint與exact-first下一步。
- **版面修正紀錄**：官方reader在有垂直scrollbar時，其`100vw` top bar會造成8 px overflow；早期failure screenshots保留為QA證據。最後以answer-first結論＋一張有效extraction native funnel收斂單頁內容，未更動任何數值、證據等級或科學主張。

### 2026-07-17 17:56 — [驗證][文件] Joint groups → minimum-extra candidates → likelihood 教學HTML

- **服務目標**：G3（正確說明read-level partial evidence）、G4（固定V／T／complete／abstain口徑）、G5（教授可讀且可回到source/test/receipt）。
- **Task type／scope**：F Demo / Illustration。這是方法教學，不是chr1–22 × 7 technical datasets的正式數據驗證；頁首、證據卡與頁尾均明示`DEMO`、正式ranking pilot=`NO-GO`、full-data=`UNAVAILABLE`。
- **Input**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`、`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py`、`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py`、`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md`與frozen v2正式pilot NO-GO紀錄。
- **方法結論**：partial pattern不逐completion建樹；每個`G(p)`形成一條group-hit限制，全部groups在同一MILP中聯合滿足。先求全域`h*`，再固定objective並用no-good cuts列舉distinct optimal vertex sets；只有solver證明再無其他解時`complete=true`。Likelihood以全部informative molecules、BQ與error model比較不同`V`，`X`作missing marginalization；同一`V`的不同parent-edge assignments `T`目前同分。
- **HTML command／output**：以`apply_patch`建立`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_聯合GroupConstraints到Likelihood排序方法教學_01.html`；單檔、無外部resource。Final HTML SHA256=`c7afc44a831f668f24ce362309e90d870976c4a2af7df6619aa1564d1b1271dd`。
- **Method-test command**：工作目錄`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/`；`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest tests/test_hypercube_exact.py -v`。
- **Method-test actual**：exit=`0`；23/23 PASS，1.156 s。涵蓋`RRA+AAA+RAX`全域minimum-extra反例、`AX+XA` joint coverage、k≤6 symbolic/explicit等價、independent exhaustive oracle、likelihood edge invariance與global KKT bound。Current test SHA256=`a86c23e981cf92b2a41c1d86e778d5bba6e6a56d0ba1743180de5534fb3313de`；舊method audit表內test SHA已過期，引用測試版本時應使用current SHA。
- **Browser/print QA command**：工作目錄`InterSubMod/`；`/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/qa_joint_group_method_html.py`。
- **Browser/print QA output／actual**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/20260717_joint_group_method_html_qa/`；exit=`0`。1440／1024／390／320 px皆無document overflow；7/7 details展開／收合PASS；print details visible；SVG title/desc、table caption/scope、單一h1皆PASS；external resources、console errors、page errors皆0。QA receipt SHA256=`2830252a631ef16aba4a6a09595d795eeee2bd9cf01bcdf3a468b8752d298585`。
- **三視角獨立審閱**：算法契約審閱確認不可稱`2^u` tree worlds或first-success；效能審閱將canonical cache／persistent solver列P0 exact-preserving，per-unit deadline與primary-only ranking列為coverage／contract change；教授版審閱要求方法AUDITED與formal pilot NO-GO並列呈現。三者均已反映於HTML。
- **改進優先序**：(P0 exact-preserving) canonical structural cache、persistent solver＋incremental no-good cuts；(P1) atomic checkpoint/resume＋deterministic sharding、shared emission/fit cache、B20只補bootstrap；(P2 research) implicit subcube/lazy constraints但須先證明等價。First-success、硬補X、未證明beam/top-N與read/VAF edge cost不得包裝成目前exact方法的透明加速。
- **證據邊界**：`33 units／96.153683%`只解釋截斷pilot的candidate-generation長尾，不是生物比例。現有formal ranking仍NO-GO，故本HTML不提供最終T／Topo或跨樣本比例。

### 2026-07-19 11:04 — [偏離][決策][驗證] all-sSNV audit 資源閘門假陰性修補與 v4 source identity

- **服務目標**：G4（全量工作不與既有高負載研究作業互相污染）、G5（preflight decision 可由程式與receipt獨立重現）。
- **偏離／根因**：v3 runbook文字要求任何既有all-sSNV作業自然結束後才啟動M2，但兩個runner的`conflict_kind()`只辨識`analyze_all_ssnv_focal_alt_multigroup*.py`；現場新版`audit_cooccurrence_task_contract_preflight.py --workers 40`未被涵蓋，造成preflight回報`process_count=0`的假陰性。
- **Hotfix**：`run_full_m2_extraction.py`與`run_full_m2_ranking.py`同步辨識`audit_cooccurrence_task_contract_preflight*.py`（含source-locked basename），歸類為`all_ssnv_cooccurrence_audit`；extraction既有synthetic process-family測試加入parent＋worker，ranking新增對稱回歸測試。
- **Focused command／actual**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_extraction.py research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_full_m2_ranking.py`；68/68 PASS，16.641 s，exit 0。
- **真實preflight input／command**：origin manifest=`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json`；以live extraction runner、`--preflight-only --workers 2 --samtools-threads 1 --mapq-min 20 --baseq-min 20 --bridge-thresholds 1,2,3,5`測試v4 bootstrap probe。
- **真實preflight actual**：exit 2；`active_conflict_process_count=41`、`active_conflict_root_count=1`、kind=`all_ssnv_cooccurrence_audit`、group CPU約`3900.2%`；`zero_conflict_pass=false`、`disk_pass=true`、可用`767,689,805,824 B`、300 GiB reserve gate不變。沒有建立formal run root、沒有使用ignore gate、沒有終止既有作業。
- **Source identity decision**：兩個frozen sources SHA已改為extraction=`1cbdb3ba603de060a6ddb5143c89c796787378d305ee1a4ae9cc102c0765e0b5`、ranking=`d2b060a6798e9ae966bf2a46c8755ff73ac71c5216ad7507362c518e9ac793a0`。保留歷史v3文件，新增`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260719_M2_exact_preserving資源閘門修補正式執行Runbook_04.md`，formal roots改為`m2_frozen_release_v4{,_attempts}`，不沿用v3 remediation binding。
- **Full regression command／actual**：同Python環境執行`python3 -m unittest discover -s research/20260716_read_linked_hypercube_exact_likelihood_validation/tests -p 'test_*.py'`；342/342 PASS，174.228 s，exit 0；v4 exact-11 allowlist test PASS。
- **未決**：formal Task B仍等待現有41-process audit自然完成；其後才可執行v4 PRE→freeze→DORADO chr6／H2009 chr2雙pilot。另需先完成44 hard-tail／25 structural-key optimized backend stress，不能把三個bounded fixtures的加速外推成full 154。

## 📚 Lore

### 2026-07-16 — `RRA + AAA + RAX`的正確解讀

- **Constraint**：`RAX`代表`{RAR, RAA}`這一個group terminal；candidate vertex set只要包含其中至少一個state即覆蓋該read pattern。
- **Consequence**：它不等於同時觀察兩種state，也不直接支持`RAR -> RAA`。在目前minimum-hidden objective下，若root/full constraints固定，可能由`RAA`得到最小可行vertex set；其他最小可行sets若存在也應全部保留。
- **Why it matters**：避免把一筆partial read展開後重複計數，也避免製造不存在的edge support。

## Provenance Footer

- **Commit hash**：pending；working-tree base `fd4c7c2cfbe8ffc0081cb7addba0313507379676`
- **Build time**：2026-07-16 06:18 +08:00
- **Skill**：`/implementation-notes` v0.1
- **Pre-decision**：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md`
- **Cycle**：`InterSubMod/state/cycles/cycle_20260716-0618-read-linked-hypercube-exact-likelihood/`
