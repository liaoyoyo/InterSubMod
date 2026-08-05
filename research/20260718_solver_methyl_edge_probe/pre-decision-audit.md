<!--
建立時間: 2026-07-18 02:24 +08:00
目標: 小步驗證 persistent MILP、uncovered-group exact B&B、vertex-set-level read-AF cache 與 methylation edge tie-breaker 是否值得進入正式實作
處理範圍: synthetic exact cases + 至多兩個既有 real capped units + HCC1395_DORADO chr6 methylation feasibility scan；PARTIAL exploratory pilot
cycle_id: cycle_20260718-solver-methyl-edge-probe
topic: 20260718_solver_methyl_edge_probe
status: verdict_PROBE
audit_version: 1.0
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py
  - InterSubMod/docs/CURRENT_FOCUS.md
-->

# Pre-Decision Audit：solver 與甲基 edge 小步驗證

> **Verdict：PROBE（50/100）**。允許在全新隔離目錄做 synthetic 與極小 real pilot；不允許替換 canonical solver、全樣本執行，或把甲基距離寫成祖先／真實 edge 證明。
>
> **服務目標**：G3（read-level epigenetic evidence）與 G5（可重現、可外部檢查的 exact solver）。

## §0 任務分類與範圍

- **Task Type**：A — Exploratory pilot。
- **PARTIAL scope**：synthetic k≤5 exact oracle、`AAAA,k=4` 24-optima stress、最多兩個既有 capped real units，以及 HCC1395_DORADO chr6 methylation feasibility。
- **Domain**：Complex；採 safe-to-fail probe。
- **禁止推廣**：本輪結果不得代表 chr1-22、全樣本 runtime、edge accuracy 或 clone/subclone accuracy。

### 研究啟動五問

1. **Thread D？** 是；甲基只作同一 vertex set 內的 exploratory edge regularizer。
2. **Thread B 撤回範圍？** 不直接涉及 FP filter；但既有「methylation cluster ≠ biological subclone」結論仍是硬邊界。
3. **KDE-corrected？** 不適用於結構 solver；甲基資料直接使用 MM/ML，不用 KDE 結論。
4. **VCF caller AF？** 不需要；現行 primary score 是 read-pattern mixture likelihood。不得用相同 reads 再重複加入 VAF。
5. **長計算／C++／搬移／NO-GO？** 不改 C++、不搬檔；每個 pilot 設短時限，只寫本研究目錄與 `/tmp` 隔離依賴。

## §1 已知觀察

| 觀察 | 狀態 | 證據 |
|---|---:|---|
| 現行 SciPy `milp` 每找到一個 candidate 會重建整個模型 | ✓ | `hypercube_exact.py:285,412,466` |
| `AAAA,k=4` 有 24 個最佳 vertex sets，因此現行路徑 build/solve 25 次 | ✓ | 獨立 read-only audit |
| SciPy 1.13.1 wrapper 沒有可增量 `addRow` 的 persistent handle | ✓ | 環境與 API audit |
| 本環境目前沒有 `highspy` | ✓ | `ModuleNotFoundError` |
| persistent model 只能減少 structural rebuild，不能消除 V+1 次求解與 V 個輸出 | ✓ | no-good enumeration contract |
| 現行 `rank_unit` 已先按 `vertex_sets` 做一次 mixture fit，再解析計數 parent assignments | ✓ | `build_m2_patterns_and_rank.py:886+` |
| 同 vertex set、不同 parent edges 的 read likelihood 必須同分 | ✓ | 現行模型註解、既有 regression |
| HCC1395_DORADO chr6 在 MINREAD 1/2/3/5 都沒有同一 component×HP 同時出現 10、01、11 exact states | ✓ | 本輪 lossless sidecar feasibility scan |
| HCC1395 是 hyper-diploid（ploidy 2.85），大量 gain/LOH；甲基 edge 正例需 CN gate | ✓ | KB `02_samples/hcc1395.md`，last verified 2026-07-11 |
| read×CpG cluster 是 methylation similarity，不是已驗證 biological subclone | ✓ | KB `05_tools/intersubmod.md`，last verified 2026-07-11 |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | MILP incremental cuts、exact branching completeness、vertex-set likelihood invariance皆可形式化 |
| 觀察支撐 | 10 | 已知 build bottleneck、real MILP optimum pilot與現行 cache路徑，但尚無新 backend實測 |
| 機制清晰度 | 10 | solver三層可分；甲基 edge score仍依賴可逆性較弱的相似性假設 |
| 反例風險 | 0 | 舊 full pilot曾8小時 timeout；tie數量仍可爆炸；CNA/LOH與甲基收斂可使 edge誤導 |
| 所需資源 | 10 | synthetic短；新 backend與real enumeration仍有版本、尾端runtime及輸出量風險 |
| **TOTAL** | **50 / 100** | **PROBE** |

## §3 Falsifiers 與停止條件

1. persistent 與現行 SciPy/exhaustive oracle 的完整 vertex-set digest 任一不一致 → **FAIL**。
2. 未先固定 `sum(extra)=h*` 就加入 sparse no-good，或把高一層合法解錯刪 → **FAIL**。
3. B&B 在 seeded k≤5 任一 case 的 optimum／完整最佳 set 不等於 MILP/exhaustive oracle → **FAIL**。
4. B&B node count/runtime在強 partial 小案例沒有改善，或記憶體因 memo/output迅速放大 → 僅保留替代 backend，不升 canonical。
5. primary read-AF fit呼叫次數不等於 distinct vertex sets，或不同 edge得到不同 read-AF → **FAIL**。
6. synthetic methylation「11接近10／接近01／等距」不能分別回復 T10／T01／TIE → **FAIL**。
7. 真實資料缺 MM/ML、三個完整 state、共同 CpG、HP×PS一致或 CN-clean gate → **ABSTAIN**，不得補值。

## §4 Assumption Map

| 假設 | 重要性 | 已知？ | Action |
|---|---:|---:|---|
| 新 `highspy` 與 SciPy內嵌 HiGHS 給相同整數解集合 | High | No | pin version；完整 set digest雙 backend對照 |
| sparse no-good在固定h*後只排除一個 set | High | Yes | regression含高一層反例 |
| uncovered-group branching涵蓋每個可行延伸 | High | Yes | formal proof + k≤5 exhaustive oracle |
| hidden node少、partial constraints強是常見實例 | High | Partial | small real capped units只作描述 |
| snapshot read likelihood可分parent edge | High | Yes=False | 同V只算一次，edge維持 unresolved |
| methylation沿真edge保持相似且不收斂／反轉 | High | No | synthetic sensitivity；real資料只稱 regularizer |
| HP/PS已排除copy-number混淆 | High | No | CNA/LOH/allele-specific amplification gate |

## §5 Quick Pilot

1. **Persistent solver**  
   一次建模求 `h*`，固定 cardinality，再增量加入 sparse no-good rows。  
   → 驗證：`AAAA,k=4` 完整24 sets；model builds=1；sets digest等於現行25-build結果。
2. **Exact uncovered-group B&B**  
   對 orphan predecessor 或最小 uncovered group 做完備分枝，使用安全 lower bound與 memo。  
   → 驗證：seeded k≤5 optimum與完整sets 0 mismatch；回報 visited/pruned/runtime。
3. **Vertex-set-level read-AF**  
   對相同V的兩個 parent assignments插入instrumentation。  
   → 驗證：fit calls=distinct V=1、parent choices=2、score完全相同。
4. **Methyl edge tie-breaker**  
   `Δ = D_M(01,11) − D_M(10,11)`；`Δ>0`只表示模型偏好 `10→11`。  
   → 驗證：正向、鏡像、等距、uninformative與missing-CpG synthetic controls。
5. **Real feasibility gate**  
   先查同HP×PS的10/01/11、MM/ML、共同CpG與CN/LOH。  
   → 驗證：不合格即明確 `ABSTAIN_NO_EVALUABLE_REAL_UNIT`。

## §6 Conflict Scan 與 Red Team

| 既有結論／反例 | 關係 | 處理 |
|---|---|---|
| HCC1395 chr6 formal M2排名曾8小時timeout | 直接風險 | 本輪只跑短 synthetic／單元，不升full |
| methylation-only subclone usability為negative | claim conflict | 只測edge-similarity scorer，不稱clone |
| same-V/different-edge對snapshot sSNV不可辨識 | 方法硬限制 | read-AF只算一次；甲基另列exploratory score |
| allele-specific amplification可在同cell/同HP產生10與11 | 生物反例 | CN/LOH confounded即ABSTAIN |
| methylation可逆、收斂或受strand/CpG coverage驅動 | 模型反例 | coverage matching、bootstrap/permutation與TIE |

獨立 red team 的最強反方：

- **Persistent**：仍需 V+1 次 `run()`；候選全集很大時，重用模型不能解決 output explosion。
- **B&B**：若 group subcubes高度重疊且下界弱，仍會接近 `O(2^M)`。
- **甲基**：即使 `Δ` 顯著，也可能是 CNA、LOH、甲基狀態收斂或 technical geometry，不是真parent。

## §7 決策

- **允許**：全新 research scripts/tests、`/tmp` pinned dependency、synthetic與至多兩個small real units。
- **不允許**：修改 canonical C++/Python、全資料執行、publication claim、真實edge／clone／subclone宣稱。
- **升級條件**：正確性全PASS，且至少一項 solver在預註冊case有實質改善；甲基需另找到滿足三state+共同CpG+CN-clean的real unit才可升下一個PROBE。
- **Claim ceiling**：solver-certified regional mutation-state candidate sets；甲基最多是 methyl-similarity regularized candidate-edge preference。

## 2026-07-19 Amendment — frozen hard-tail stress（仍為 PROBE）

> **Verdict：PROBE_STRESS。** 核准從 frozen HCC1395_DORADO chr6 diagnostic pattern snapshot 建立不可變 stress panel與generic dual-backend harness；不核准修改canonical solver、啟動full 154，或把pattern replay稱為raw-BAM全量驗證。

- **Primary scope**：structural minread=1的44個hard unit instances／25個distinct structural keys；選取規則固定為33個candidate-incomplete units加11個complete但`V≥100` controls。
- **Sensitivity**：minread=2的36個hard instances／21個distinct keys；不得與primary合併成獨立樣本。
- **可恢復邊界**：pattern snapshot可無損恢復hard units的k、full/partial structural patterns、observed-ALT universe與BQ/count groups；33 tails尚未求出的true V/family digest與30個非factorial cases的歷史停止原因不可恢復，必須由新solver重算。
- **預先固定成功門檻**：11個complete controls的objective/V/family digest零mismatch；所有complete output均通過coverage、rooted-predecessor、objective與duplicate validator；任何deadline/cap均`family_complete=false`、`ranking_allowed=false`、winner=0；primary total candidate wall至少5×改善、p95不劣於current、peak RSS不超過current 1.5×。
- **大family oracle**：三個single-terminal keys以`V=w!`獨立核對`720／5,040／720`；允許`max_sets=None`的bounded materialization，但不得宣稱對任意V具有限記憶體保證，也不得materialize parent-tree Cartesian product。
- **H2009／H1437**：目前沒有M2 extraction artifacts或預註冊unit selection；只能列為後續panel。必須先建立100% route census與固定抽樣／全tail規則，不能事後挑有利cases。
- **資源條件**：現場另有41-process all-sSNV audit時只允許小型tests／fixture construction；正式stress timing與M2 pilot必須等zero-conflict resource gate通過，否則wall speedup無可比性。
- **Promotion gate**：hard 25/25 keys與H2009/H1437後續panel皆通過correctness/fail-closed/resource gates後，才可提出dual pilot adapter；在此之前optimized backend維持隔離research source。

## 2026-07-19 Amendment — compressed VAF/read-AF ranking feasibility

> **Verdict：PROBE_COMPRESSED_RANK（Complex / Task A）。** 只核准在隔離
> `compressed_vaf_rank_probe_v1` 中測試 DP traceback lazy enumeration 與
> certified upper-bound pruning；不修改 canonical ranker，不跑 production，
> 也不把有界未完成搜尋稱為 exact winner。

- **研究問題**：當 exact perfect-family counter 已證明一個 hard key 有
  `122,281,152` 個最簡 vertex sets 時，能否不列完整 family，仍以現行
  BQ-aware profile mixture likelihood 找到全域最佳 vertex set與完整
  `tie_tolerance` tie class？
- **同 endpoint 範圍**：只比較
  `build_m2_patterns_and_rank.py::fit_quality_aware_mixture` 的 primary
  `log_likelihood` 與 `_rank_fits` tie 規則；fixed-error grid、bootstrap、
  responsibility TSV與所有候選逐列輸出不在本 probe 的同 endpoint 宣稱內。
- **安全上界**：對一個未展開 traceback branch，取該 branch 所有可能 state
  的 union \(U\)，在 \(U\) 上擬合相同 mixture likelihood。因每個 completion
  的 convex hull 都包含於 `conv(U)`，`LL(U)+KKT_gap` 是可證的 branch upper
  bound。只有當 `UB < incumbent_lower_bound - tie_tolerance` 才可剪枝。
- **成功條件**：`m<=5` exhaustive oracle 中，最佳 score、best vertex set IDs、
  完整 tie IDs皆零 mismatch；所有分數皆 KKT-certified；hard key只跑短
  deadline／candidate budget並產生可機械檢查的 complete 或 fail-closed receipt。
- **失敗／abstain 條件**：任一 optimizer 非 certified、branch upper bound
  非 finite、deadline/candidate/tie-output cap、traceback count不守恆，均不得回
  exact；改回 `INCOMPLETE_*`、`ranking_complete=false`、winner/tie不發布。
- **最強反例**：若 reads不具資訊或所有候選 emission hull近乎相同，全部
  `122,281,152` candidates都可能落在 tie tolerance 內；一般 current likelihood
  非 edge-additive，因此 DP count 本身不能壓縮完整 tie membership。這個 probe
  可以證明「某一 instance 被安全剪枝後 exact」，不能預先保證所有 hard cases
  多項式時間完成。
- **資源界線**：synthetic `m<=5` + 單一 frozen hard `09b` 的秒級 bounded probe；
  不啟動全25-key stress、M2雙pilot或154 tasks。
- **Review後實作界線**：v1只完成small oracle；現有UB只切在top-level traceback
  branches，不是recursive DP-state B&B。後續hard probe須寫全新v2 receipt，
  且先驗R3 pointer/sidecar、schema/status、path containment、manifest
  byte/semantic integrity、frozen checks及case/unit records；不得把v1 receipt
  外推成current-source hard evidence。

### 2026-07-19 P1 amendment — numerical certificate與TOCTOU fail-closed

> **Verdict：PROBE_P1_HARDENING。** 服務G3/G5；只修改隔離research
> prototype/runner/tests與v2文件，不修改canonical、不執行09b hard。

- **研究問題**：普通IEEE-754 float relaxation即使數學上是upper bound，未使用
  directed rounding時，是否可能因roundoff而錯剪tie/winner並冒稱machine-exact？
- **假設**：結構DP/count可維持exact integer certificate；數值UB只能標
  `numerical_bound_certified=false`。若任何候選靠此UB剪除，core不得發布
  authoritative best/ties；只有`evaluated==logical_family_count`的全實評結果可
  machine-publish。小型獨立exhaustive oracle可另標`oracle_confirmed=true`，
  但不能改寫一般core證書。
- **成功條件**：
  1. 有float-bound pruning時，`ranking_complete=false`且authoritative fields隱藏。
  2. 無prune且全實評時，best/tie可發布。
  3. deadline/candidate/tie controls拒絕bool、非原生int、非finite與越界值；
     candidate LL/gap/status全部finite且符合current optimizer contract。
  4. branching、mandatory-full、mixed R/A/X、gapped-active m≤5 fixtures逐一驗證
     branch互斥、count守恆、每completion vertices包含於其branch
     `possible_vertices`，並與獨立parent-vector oracle一致。
  5. runner執行前後重算所有source/input hashes；任一漂移即不寫acceptance
     receipt。v2綁`hypercube_exact.py`與NumPy版本。
- **失敗／abstain條件**：任何hash drift、非finite分數、控制型別錯誤、
  branch overlap／containment錯誤、partial-search或float-pruned search都不得
  標machine-exact。
- **複雜度口徑**：目前會materialize全部top-level branches與每個branch的
  `possible_vertices`；最壞約`O(m·3^m)`時間／materialized state work，
  並非recursive B&B。
- **Deadline界線**：core deadline目前在DP與all-branch UB建構後才檢查，故不是
  hard wall bound。09b維持NOT_RUN，直到有外層process watchdog；禁止把
  `--hard-deadline 5`宣稱成真正5秒上限。

### Cross-sample historical-v5 amendment（PROBE_CORE_ONLY）

- **目的**：在真正M2 H2009／H1437 extraction尚未存在前，先用canonical v5 `mlhp_part_1..5.json` 的`subread_groups_by_hp`驗solver core；不可替代M2 read-linked component／PS block／threshold route。
- **母集合**：H1437 historical primary incomplete `1,773/1,773`、H2009 `3,759/3,759`，pattern join missing必須為0；distinct structural keys分別`1,763/3,746`，union=`5,506`。
- **Deterministic panel**：100%納入48個reduced `q>8` keys；其餘依`sample×k×effective_m×reduced_q×historical_cap_class`每層取structural SHA字典序最小者，預期union=`208` keys（H1437=86、H2009=122），key-list digest=`7ae4f85ba9f05f146303d3f7921e5f333f395c690cdcd1dbd84ea39037b75692`。
- **允許**：建立source-bound exporter／manifest／tests與短dry-run；zero-conflict前不跑正式208-key timing。
- **Gate**：12個source files與兩份output manifest SHA全吻合、母集合守恆、selection重算digest一致、任何缺pattern／重複unit／事後抽樣均fail-closed。
- **Claim ceiling**：historical-v5-derived solver-core stress；不得用其complete率或runtime宣稱M2全樣本end-to-end。

### Cross-sample perfect-family count census amendment（PROBE_CORE_ONLY）

- **允許範圍**：只對已凍結的208-key historical-v5 panel執行`perfect_family_count.py`解析計數；每個key只能回`exact`或明確`abstain`，不得呼叫current／optimized顯式family solver。
- **用途限制**：本census只回答recurrence-free最小層的family count與objective certificate；所有rows固定`ranking_allowed=false`，不得進行VAF排序或宣稱候選winner。
- **效能限制**：可記錄本次bounded counter的finite elapsed作可重現診斷，但現場仍有外部作業，固定`formal_speed_claim=false`；不得與既有solver runtime換算正式speedup。
- **完整性Gate**：輸入manifest byte/content hash、counter／runner source hash皆須綁定；208/208 rows必須total且finite，不相容case或未知status一律fail-closed。
- **Promotion限制**：結果維持`HISTORICAL_V5_SOLVER_CORE_ONLY`；不得取代M2 end-to-end驗證、不得修改canonical、不得作production promotion。

## 2026-07-19 Final Gate — NO-GO_PRODUCTION

> **Verdict：NO-GO_PRODUCTION。** Hard25 controller雖有50/50 worker rows，
> optimized family completion只有16/25；9個incomplete中6個perfect family
> 最多達122,281,152個，另3個需要recurrence-aware method。故不能確認所有最簡
> 候選列完、所有VAF/read-AF排名完成或25-key overall completion speedup。

- **通過證據**：雙方皆objective-certified的20/20完全一致；雙方皆family
  complete的3/3 objective/count/digest一致；factorial 4/4 PASS；
  incomplete-ranked=0；compressed m≤5 full-enumeration oracle PASS。
- **失敗證據**：current complete=3/25、optimized complete=16/25；
  controls未全部paired-complete；兩backend皆有deadline-censored runs，
  `reported_ratio_kind=NOT_COMPARABLE`、ratio=null。
- **數值剪枝限制**：普通float UB沒有outward-rounded interval certificate；
  若`pruned>0`必須`ranking_complete=false`且不得發布winner/ties。
- **停止動作**：不得啟動production adapter、雙actual-data pilots或full 154。
  這些未執行項目是gate要求，不是遺漏。
- **重新開啟條件**：recurrence-aware traceback、recursive DP-state B&B、
  interval-certified UB、outer process watchdog完成後，interleaved repeated
  hard25須達25/25 family/ranking complete與零correctness mismatch。

## Provenance

- **Current-source warning（2026-07-19）**：以下commit與solver/ranker SHA是本
  research cycle早期歷史快照，不能代表main compressed-VAF revalidation的
  current source。Current綁定必須讀
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v2/main_reverify_r3/receipt.json`
  的`source_bindings`；其中hypercube=`1def0de1952d127d8d013820ac7b0eabe33d6f1a66fd8c6d47a1985b14a32f77`、
  ranker=`b28d494563bea70cbfea7b13853d99c0b3aeafed7ed4f0c0a6fff3b92c84f0a9`。
- **Git preservation warning**：本cycle目錄目前未被Git追蹤；
  `git ls-files research/20260718_solver_methyl_edge_probe`為0。下列commit不包含
  本輪frozen receipts/evidence，不能單靠commit還原。
- Git commit：`0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Current exact solver SHA256：`9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`
- Current ranker SHA256：`c82210f25870d1405cc070c18096fb7d1c2b908fb4d8a7858aece7ac4b151d28`
- Runtime baseline：Python / SciPy 1.13.1；`highspy` absent at audit。
- Resource：`/big7_disk` 99% used、716G available；pilot outputs必須保持小型。
