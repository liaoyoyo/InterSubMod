<!--
建立時間: 2026-07-16 07:20 +08:00
目標: 驗證 M2 lossless read extractor 的漏斗、CIGAR、sidecar join、per-HP linkage、receipt 與 resource gate
處理範圍: 純程式審查、tiny synthetic BAM/VCF/sidecar、canonical v5 manifest preflight；未執行新的大型 BAM 掃描
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/lossless_read_contract.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/extract_lossless_read_linkage.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_extraction.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/implementation-notes.md
-->

# M2 extractor 程式與合成資料 QA

## 結論

extractor schema `1.2.0` 與 M2 ranking schema `2.0.0` 已完成第二輪 static red-team與P0/P1修正；本文件該輪audit-time synthetic regression為130/130。它已能把 pooled genomic scope 與 HP×exact-PS linkage 分開，並以 fail-closed receipt/reuse 契約保存 raw alignment、filter、call-reason、component partition、T/V/topology、optimizer certificate 與 per-unit monotonic wall-time 診斷。舊 HCC1954 chr22 pilot v1 是 schema 1.0／pooled-only，只能保留作歷史 performance pilot；新版結果必須另用新目錄重跑。**本文件尚未產生 154-task 全量數據，因此程式/語義判定 PASS，但全量資源與最終研究數字仍為 PROBE。**

## 修正與判定

| 檢查面向 | 原始風險 | 修正後判定 | 驗證證據 |
|---|---|---|---|
| pooled vs HP linkage | pooled read bridge 可能把 HP1、HP2 或不同PS各自不連通的位點合成一個可建樹區 | PASS；每個 threshold 輸出diagnostic `pooled/HP1/HP2/HP3/HP4/unphased`及primary `PS_HP1/PS_HP2` partition；primary component ID綁定exact PS | 反例：pooled `(1,1)`形成 k=3；HP1 `(1,0)`與HP2 `(0,1)`分別切開；同site可屬不同PS但read支持不混合 |
| raw funnel | 只檢查 raw=flag+MAPQ+eligible，未檢查 alignment class、writer row與call reason mass | PASS；新增 class、written row、R/A/O/D/S/L/X、fixed R/A、ALT mass守恆 | tiny fixture：raw 8 = flag 4 + MAPQ 1 + eligible 3；eligible/written/sidecar match皆3 |
| CIGAR／座標 | insertion、deletion、refskip、低BQ、other base與span boundary可能混淆 | PASS | synthetic tests分別驗 `I/D/N`、O、L、X與0-based half-open alignment span |
| HP vocabulary | `startswith("1")`可能將`1foo`靜默歸HP1 | PASS；僅接受 `.,1,2,3,4,1-1,2-1,1-2,2-2` | `0/1foo/2-3/5`全部raise `ValueError` |
| multi-primary molecule | 同 dataset＋RG＋QNAME有兩個eligible primary alignment可能重複計數 | PASS；第二筆fail-closed | tiny BAM實際經`samtools -M -L`觸發預期`RuntimeError` |
| sidecar exact join | 同 exact identity conflicting HP/PS可能被覆寫 | PASS；衝突直接停止 | duplicate exact key、不同HP synthetic test |
| receipt 自身hash | JSON不可能穩定內含其最終bytes的SHA；舊版只在記憶體補receipt identity | PASS；`receipt.json.sha256`保存exact-byte外部identity | receipt tamper後reuse必定失敗 |
| resume reuse | 舊版只看`all_pass=true`，可能重用錯scope／舊參數／壞output | PASS；驗schema、scope、參數、manifest hash、extractor hash、每個output size+SHA | output、receipt、參數、scope與producer tamper tests全部fail-closed |
| resource gate | `.pinned_<sha>.py`漏抓；`/usr/bin/time` parent可能被誤判self-conflict；44 children灌滿JSON | PASS；以實際argv token辨識、排除current ancestor/descendant、依root聚合 | canonical preflight抓50 processes、聚合2 roots、exit=2、未建立output root |

## 執行命令與實際輸出

### 全測試

```bash
python3 -m unittest discover -v \
  -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests \
  -p 'test_*.py'
```

- 輸入：上述 scripts/tests；整合測試另於系統暫存目錄建立8–9 alignments的tiny BAM、5-site VCF與tabix sidecar。
- 輸出：terminal test receipt；暫存fixture自動移除。
- 實際結果：exit `0`；`Ran 50 tests in 3.285s`；`OK`。

### Canonical resource preflight

```bash
python3 InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_extraction.py \
  --manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json \
  --output-root /tmp/m2_preflight_qa_not_created \
  --workers 1 --samtools-threads 1 --mapq-min 20 --baseq-min 20 \
  --bridge-thresholds 1,2,3,5 --preflight-only
```

- 輸入：canonical v5 manifest與當下process table；不開啟BAM。
- 輸出：stdout JSON preflight；正式output root未建立。
- 實際結果：exit `2`、`resource_gate_pass=false`、50 conflict processes聚合成2 roots（4與46 members）。

## 來源identity（第一輪歷史快照）

| 檔案 | SHA256 |
|---|---|
| `scripts/lossless_read_contract.py` | `e378f7186b262c3629c45eb388204e65a302acf3b9f4f95ffd0be835c47862de` |
| `scripts/extract_lossless_read_linkage.py` | `3d62c68d8e97818224c0ecb66ced2ab4d5c9a5936a3882316cc3f3df8b49f446` |
| `scripts/run_full_m2_extraction.py` | `27fe3f6b74e069300932019195278cecfe1ffbd67f74f8bb94be11d2cc9bac16` |

## 限制與下一個強制步驟

1. 本QA證明工程契約與合成反例正確，不是 chr1–22 × 7 datasets的最終數據。
2. 新版 HCC1954 chr22 pilot必須寫入新目錄，且receipt需為schema 1.2.0；不可覆寫v1。
3. 正式樹推論應以`PS_HP1/PS_HP2 × exact known PS` component為主要單位；pooled/legacy HP component只可作genomic scope或sensitivity，不可直接宣稱within-PS linkage。
4. `O/D/S/L/X`在extractor保留原因；下游likelihood目前將其marginalize為missing。若未來導入校準過的indel/other-base emission model，需另升schema與重跑驗證。

---

## 第二輪 schema 1.2／2.0 static red-team 與修正（2026-07-16）

### 1. 任務範圍與判定

- 任務類型：**B — Comprehensive validation**；服務 G3（read-level證據）、G4（跨樣本可重現性）、G5（外部可驗證）。
- 處理範圍：逐行 static audit、tiny synthetic BAM/VCF/TSV、symbolic/MILP/likelihood tests、記憶體 microbenchmark。
- 未處理：未讀 canonical 大型 BAM、未啟動 154 tasks、未產生全樣本最終數字。
- 最終判定：
  - 方法語義、PS isolation、partial-read 契約、T/V/topology守恆：**PASS**。
  - 已發現的 P0/P1 程式問題：**已修正且 regression PASS**。
  - 154-task 的 RAM materialization：**PASS（已串流化）**。
  - 154-task 的實際磁碟量、wall time、solver tail：**PROBE；需先做最大型單一 chromosome resource pilot**。

### 2. Static audit 問題—修正—驗證矩陣

| 嚴重度 | 原始問題 | 修正 | 驗證結果 |
|---|---|---|---|
| P0 | 每個 `(HP family, PS)` 都配置並掃描整條 chromosome site array，形成 `O(#PS × #sites)` RAM/時間 | `sparse_support_at_active_boundaries`只排序非零 bridge events，two-pointer 查詢 active cuts；legacy cross-PS audit改為 threshold interval union，不再保存 dense PS×site matrix | sparse 與舊 dense reference完全等價；`n_sites=1,000,000,000`、2 events可立即完成；tiny extraction與PS isolation通過 |
| P0 | `build_evidence` 每個 molecule 掃描全 chromosome 的所有 PS/component mappings，形成 `O(#molecules × #PS × thresholds)` | 一次建立 `(raw family, exact PS)` 與 family wildcard route index；每個 molecule只查2個index bucket | 100個PS反例：naive reference 100 mapping visits，實際 route visits=1；僅PS=42 unit收到1筆，其他99 units為0 |
| P1 | full extraction aggregate用`Counter.update`把`max_k_component_sites/max_k`跨154 tasks相加 | count欄位才sum；兩個max欄位獨立用`max()` | synthetic child max=5與7：aggregate max=7、`n_components`仍正確相加為4 |
| P0/P1 | child ranker同時materialize全部 pattern、candidate、responsibility rows，可能在單一task先耗盡RAM | 每個unit直接寫入staging gzip；semantic hash、row count、optimizer diagnostics在線累積；成功後原子改名 | sink輸出與舊retained rows逐列相同；100,000列probe峰值由35,175,566 bytes降至5,392 bytes，digest相同 |
| P1 | full candidate table先讀入整個child candidate TSV再group | 改為逐列讀取、每次只buffer一個unit，並檢查child/global unit排序 | 4/4 full-ranking tests通過；receipt新增`max_buffered_candidate_rows` |

對應程式證據：

- chromosome site catalog直接來自整條指定染色體PASS sSNV：`scripts/extract_lossless_read_linkage.py:208-226`。
- unique molecule bridge event：`scripts/extract_lossless_read_linkage.py:529-560`。
- sparse PS支援與threshold union：`scripts/extract_lossless_read_linkage.py:312-397`、`:811-882`。
- 50 kb僅診斷欄位，不是component boundary：`scripts/extract_lossless_read_linkage.py:622-658`；真正切割只依bridge threshold（`:707-790`）。
- exact-PS route index：`scripts/build_m2_patterns_and_rank.py:1493-1536`；每molecule只取matching exact/wildcard buckets：`:1606-1624`。
- max-vs-sum聚合：`scripts/run_full_m2_extraction.py:291-366`。
- child detail streaming：`scripts/build_m2_patterns_and_rank.py:2049-2179`。
- full candidate per-unit streaming：`scripts/run_full_m2_ranking.py:902-1021`。

### 3. Partial read 的實際方法：不是建立 `2^u` 份完整資料

使用者提出的直覺「partial read有幾種可能，就建立幾種完整情況；只要其中一種能被最小樹涵蓋便成功」在**邏輯上的存在條件**是對的，但目前實作更有效率，並不真的建立所有 completion 組合。

1. 一條molecule投影為一個固定長度R/A/X字串；例如`RAX`保留為一個獨立pattern，`AXR`、`ARX`、`AXX`不互相合併。只有相同exact R/A/X pattern的count達`structural_exact_pattern_minread`才形成structural full/group constraint；所有低於門檻與達門檻的informative molecule projections都仍進後續likelihood（`scripts/build_m2_patterns_and_rank.py:921-979`）。
2. `RAX`被編成`fixed_mask + alt_mask + free_mask`。若有`u`個X，概念上有`2^u` completions；production保留一個symbolic subcube group，不建立`2^u`個獨立tree worlds。每次MILP build／rebuild（`_build_problem()`）時，會對每個reduced active group列出active universe內至多`2^u_eff`個相容vertex indices作為該次sparse hit row係數；不是整個workflow只materialize一次。source-bound `METHOD_CONTRACT`因此同時宣告`active_compatible_vertex_indices_materialized_for_sparse_rows=true`與`completion_wise_tree_worlds_materialized=false`，兩者不矛盾。程式參考：`scripts/hypercube_exact.py:36-91,208-282,285-408`。
3. MILP先固定root `RR…R`與所有達structural門檻的完整observed states；每個被選取的non-root state必須至少有一個Hamming distance=1且少一個ALT的已選predecessor，因此局部限制即保證有向連到root。程式參考：`scripts/hypercube_exact.py:302-345`。
4. 對每個達門檻的partial pattern，加入「至少選到一個compatible vertex」的group-terminal限制：`sum(x_v for compatible v) >= 1`。所以不是先為`RAX`另跑`RAR`版與`RAA`版；而是在同一個MILP中用一條存在性限制表達兩者。程式參考：`scripts/hypercube_exact.py:347-354`。
5. 目標函數最小化非mandatory selected vertices，也就是最少hidden states。先求最小hidden數，再固定此最佳值、逐次加入no-good cut，枚舉所有不同但同樣最小的vertex sets，直到證明無更多解，或達`max_vertex_sets=256`／time limit。程式參考：`scripts/hypercube_exact.py:325-388,412-529`。
6. 只有`complete=true`時才可說「所有最小vertex sets已列完」。若達256上限或time limit，必須abstain於全域完整性，不可把已找到的部分誤稱為全部合理解。

因此，更精確的說法是：

> 每個partial read pattern形成一個布爾子立方體的group-terminal限制；每棵候選最小樹必須在該子立方體中選到至少一個相容狀態。系統以單一joint MILP求解並枚舉所有已證明完整的最小vertex sets；每條group row可列出相容active vertices，但不建立completion-wise tree worlds，也不建立跨partial reads的Cartesian completion組合。

### 4. Read count、BQ、VAF 與拓撲各自只負責一層

1. 結構層：exact R/A/X pattern的read count只用來判斷該pattern是否達structural minread；partial compatible patterns不合併湊門檻。
2. 排序層：所有informative molecule projections依exact pattern+BQ signature分組，likelihood中每個固定R/A call的BQ emission乘一次，再乘該group count；X/O/D/S/L視為missing，emission factor=1。production emission／fit串接見`scripts/build_m2_patterns_and_rank.py:595-681`；共用certified optimizer見`scripts/hypercube_exact.py:705-839`，固定錯誤率sensitivity見`:926-973`。
3. VAF層：同一批read衍生的VAF**不再加成第二個分數**，避免把同一證據重複計算。source-bound method declaration為`parameters.method_contract.same_molecule_vaf_added_as_separate_term=false`，receipt另保留`optimizer_contract.same_read_vaf_added=false`（`scripts/build_m2_patterns_and_rank.py:70-87,2446-2466`）；它們需與score construction／regression及producer SHA共同維護。方法宣告不是可獨立量測「一定沒double count」的證據，故不能只引用hard-coded boolean作驗證。
4. 候選層：`V`是不同minimum vertex sets數；同一vertex set可能有多個parent-edge assignments，因此`T = Σ parent_choice_count(vertex set)`，必須有`T >= V`。
5. 拓撲層：snapshot read-pattern likelihood可排序vertex sets，但同一vertex set內的parent edges不應再用相同snapshot reads假裝可辨識。故edge ambiguity需保留，不能把`V=1`自動寫成exact topology唯一。

主要守恆檢查包括：

- molecule projection = informative scoring + all-X excluded；
- projected cells = fixed R/A + marginalized non-R/A；
- `T >= V`；
- winner vertex counts、winner edge counts與selection status一致；
- coarse topology class partition、optimizer convergence/KKT certificate、candidate/responsibility semantic hashes一致；
- exact known PS不互相交換read支持，missing PS不進primary。

### 5. 驗證命令與實際輸出

#### 5.1 完整 synthetic regression

```bash
/usr/bin/time -v env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 -m unittest discover \
  -s InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests \
  -p 'test_*.py' -v
```

- 輸入：上述scripts/tests與測試即時建立的tiny fixtures；未讀canonical大型BAM。
- 輸出：terminal test receipt；fixture離開temporary directory後自動移除。
- 該輪實際：exit `0`；`Ran 130 tests in 6.101s`；`OK`；wall `7.12s`；maximum RSS `86,848 KB`。

新增的直接反例：

- sparse-vs-dense bridge old/new等價；
- billion-site complexity invariant；
- same site跨兩個PS但primary read support不混合；
- 100-PS route-index complexity counter；
- count sum／max max aggregate語義；
- streaming hash與retained-row完全等價；
- full candidate per-unit streaming與optimizer-abstain語義。
- `time.monotonic_ns()`三段計時、exact nearest-rank p50/p95/p99、runtime TSV tamper fail-closed，以及production/full-independent雙重重算。

#### 5.2 100,000列RAM microbenchmark

```bash
python3 - <<'PY'
import gc, pathlib, sys, tracemalloc
scripts = pathlib.Path(
    'research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts'
).resolve()
sys.path.insert(0, str(scripts))
from build_m2_patterns_and_rank import _SemanticListHasher, semantic_digest
n = 100_000
tracemalloc.start()
rows = [{"index": i, "pattern": "RAX", "count": 1} for i in range(n)]
materialized_digest = semantic_digest(rows)
materialized_peak = tracemalloc.get_traced_memory()[1]
tracemalloc.stop(); del rows; gc.collect()
tracemalloc.start()
hasher = _SemanticListHasher()
for i in range(n):
    hasher.update({"index": i, "pattern": "RAX", "count": 1})
streamed_peak = tracemalloc.get_traced_memory()[1]
tracemalloc.stop()
print(f"n_rows={n}")
print(f"materialized_peak_bytes={materialized_peak}")
print(f"streamed_peak_bytes={streamed_peak}")
print(f"peak_ratio={materialized_peak / streamed_peak:.1f}x")
print(f"digest_equal={materialized_digest == hasher.hexdigest()}")
PY
```

實際輸出：

```text
n_rows=100000
materialized_peak_bytes=35175566
streamed_peak_bytes=5392
peak_ratio=6523.7x
digest_equal=True
```

此probe只證明detail-row list materialization已移除，不代表MILP本身的峰值記憶體，也不會減少必須寫出的磁碟列數。

### 6. 來源 identity（第二輪 audit-time snapshot）

| 檔案 | SHA256 |
|---|---|
| `scripts/extract_lossless_read_linkage.py` | `84fb041a965acd667e08a10e6aaaf5951d6de929657a5b05c81d531c2701f222` |
| `scripts/run_full_m2_extraction.py` | `465668a14b89e7b5f01acb84d7986ac84435057eab840670871981cbfb899005` |
| `scripts/build_m2_patterns_and_rank.py` | `b7a19439b037a36f25d74c5bce337fbdeec1e303a95d8bbcc460a4874c25f03c` |
| `scripts/run_full_m2_ranking.py` | `1a760045123aecdb63adff71ac7ae55fbb03054a0cc5993efca07e7ef4d59505` |
| `scripts/hypercube_exact.py` | `6576910b52cd8a3463b9bbe3d22fde5d735f5af78ae2c512f0bf4b67031d5fb8` |

### 6.1 2026-07-16 per-unit solver/runtime instrumentation 增補

- child ranker新增`m2_unit_runtime_diagnostics.tsv.gz`，每個`component×HP×PS×threshold×minread` evaluation保存：
  - `candidate_generation_elapsed_seconds`：只包exact candidate vertex-set enumeration；
  - `likelihood_fit_elapsed_seconds`：包primary BQ fits、fixed-error sensitivity fits，以及conditional bootstrap的resampling+fits；
  - `unit_total_elapsed_seconds`：完整`rank_unit` wall time；
  - 另存`candidate_generation_invoked`與`likelihood_fit_invoked`，避免未執行單位的0秒稀釋solver tail。
- clock固定為`time.monotonic_ns()`，換算單位為seconds。數值只是同一process/機器/當下load的performance diagnostic，**不是跨機可重現的科學數據**。
- child receipt對primary units與all-minread evaluations各存`n,sum,max,p50,p95,p99`；full runner從154個child TSV逐列重算primary aggregate，同時報告all-primary與invoked-only tail。
- quantile定義：exact empirical nearest-rank，一起算`rank=ceil(p×n)`，`p∈{0.50,0.95,0.99}`。full aggregate只保留5個packed float64 vectors，上限40 bytes/primary unit，不materialize TSV rows。
- independent verifier不匯入production ranker/aggregator，改用full sort獨立重算child與full runtime summaries；任一TSV、summary、invocation flag或hash被篡改即fail-closed。
- wall-time值不放入`unit_semantic_sha256`，因此統計診斷不會破壞科學row的deterministic semantic identity。
- 本次未讀BAM；真實p50/p95/p99必須由PS-aware single-chromosome pilot與後續154-task receipts產生，不從synthetic數值外推。

本次instrumentation source identities：

| 檔案 | SHA256 |
|---|---|
| `scripts/build_m2_patterns_and_rank.py` | `476539baac01e4b82a73c414288e99be557cef52baef32c9dee73be6f2dbdab3` |
| `scripts/run_full_m2_ranking.py` | `da1896e50ab8b0aa65f38eaab763b71a36e19f26c2ebc400148fbdabe4095d94` |
| `scripts/verify_full_m2_receipts.py` | `7259d3f32ceb27a9bfe151153127a249c0d76f40671d0f9cbdafef9c88d964ca` |

### 6.2 2026-07-16 source-bound partial-read method contract 增補

為避免把hard-coded receipt boolean誤當成可量測的驗證，child ranker、full runner、independent verifier與HTML builder現在各自固定並exact-compare同一份method contract。關鍵欄位是：

- `structural_group_scope=DISTINCT_EXACT_RAX_COUNT_GE_THRESHOLD`；
- `active_compatible_vertex_indices_materialized_for_sparse_rows=true`；
- `completion_wise_tree_worlds_materialized=false`；
- `cross_read_cartesian_products_materialized=false`；
- `same_molecule_vaf_added_as_separate_term=false`；
- `parent_edges_or_transitions_scored=false`。

contract本身是source-bound方法宣告，不是生物真值或獨立量測；可信度來自四個獨立位置的exact comparison、producer SHA、negative regression與full independent verifier共同fail-closed。此次migration未讀BAM、未產生新M2數據。

| current source | SHA256 |
|---|---|
| `scripts/build_m2_patterns_and_rank.py` | `1b54fe07f5bc49dd041bd9d6fdba056e59a2b20fb9092b260712bc66e88aa59a` |
| `scripts/run_full_m2_ranking.py` | `576e8c6fc1d29c9953e5cb149ef035aec141fc3df7cb5ceb7ca34d39ed693ec1` |
| `scripts/verify_full_m2_receipts.py` | `de3af4128c93496e35ac48ee3632f065b4f5a00801767c4557a3aeec02e34350` |

Method-contract migration audit-time focused regression：149/149 PASS。此數字只證明該次code／receipt contract，不替代154-task真實M2 receipts；HTML builder後續仍有active report修正，其final identity與完整回歸由主流程完成後另行綁定。

### 6.3 2026-07-16 independent profile-likelihood recomputation 增補

- Verifier現在不只比對method boolean與candidate rows；它逐child、逐primary unit串流`m2_symbolic_pattern_counts.tsv.gz`與candidate states／π，獨立重算conditional R/A BQ emission、profile LL、simplex residual、Frank–Wolfe/KKT gap、relative weights與winner/tie partition。
- 不import production ranker／solver／aggregator，也不重跑SLSQP；是在持久化π上重算concave objective與global-gap certificate。
- Full receipt要求154 child全覆蓋，profile `n_units/n_candidates`分別等於canonical candidate table `n_units/n_rows`。`n_scoring_molecules`是跨unit/threshold的molecule projections加總，不是全域unique molecules。
- 額外LL／類VAF score、BQ、π、state、KKT gap或winner/tie被篡改均由負向測試拒絕。Profile-verifier targeted 23/23；加入completion-world/joint-group與deterministic exhaustive-reduction audit後，full research suite 162/162 PASS；仍不替代真實full154 receipts。

| current source | SHA256 |
|---|---|
| `scripts/verify_full_m2_receipts.py` | `4859598d74486f4eba6e4af6fa2dec2b4c0eb5c4e8ed86feac82483e5a7f32d8` |
| `scripts/build_validated_html_report.py` | `f1ebb17f66d7856bfd3161ff20562ed98fdd8706874c280784d6ddc113054352` |

### 7. 154-task前的唯一剩餘工程 gate

串流化解決RAM list爆量，卻不會消除科學契約要求的輸出列：

- pattern rows：約為`minread grid數 × Σ quality-pattern groups`；
- candidate rows：上限受`grid × units × max_vertex_sets(256)`控制；
- responsibility rows：約為`Σ primary winners × informative patterns × candidate states`。

正式154 tasks前應先跑最大型單一chromosome task，記錄五個gzip（pattern、unit、candidate、responsibility、runtime）的實際bytes/rows、wall time、maximum RSS與solver complete/cap/time-limit比例，再用task-size與新增per-unit invoked-only p95/p99外推。只有在預估總輸出低於可用磁碟安全水位、並行RSS低於可用RAM安全水位後才可全量啟動。

**目前可以宣稱：schema 1.2/2.0方法與程式契約已通過static/synthetic驗證。現在不能宣稱：7 datasets × chr1–22的最終區域、T/V、Topo數字已完成；那些數字必須來自後續154-task passing receipts。**
