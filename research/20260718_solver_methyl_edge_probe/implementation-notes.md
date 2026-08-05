<!--
建立時間: 2026-07-18 02:26 +08:00
目標: 持續記錄 solver 與 methyl edge exploratory pilot 的設計決定、偏離、折衷與未決問題
處理範圍: synthetic + bounded real pilot；PARTIAL；不修改 canonical solver
關聯檔案:
  - InterSubMod/research/20260718_solver_methyl_edge_probe/pre-decision-audit.md
  - InterSubMod/state/cycles/cycle_20260718-solver-methyl-edge-probe/audit.json
status: bounded_optimized_backend_validated_not_canonical
-->

# Implementation Notes：solver 與 methyl edge probe

> 狀態：`bounded_optimized_backend_validated_not_canonical`／Task Type A／PARTIAL。

## 設計決定

### 2026-07-18 02:26 — 隔離 prototype

- **決定**：所有新程式、測試與收據只寫入
  `InterSubMod/research/20260718_solver_methyl_edge_probe/`；新 solver dependency 只裝到 `/tmp`。
- **理由**：現行 research source 未追蹤且 worktree 有大量使用者變更；不能把探索性 backend 混入 canonical。
- **驗證**：新目錄在啟動時為空；現行 solver/ranker SHA256 已寫入 audit。

### 2026-07-18 02:26 — solver correctness 優先於速度

- **決定**：persistent與B&B都必須對完整 sorted optimal vertex-set digests，而不只對 objective。
- **理由**：找到一個最低解不等於列舉全部同分解；後續 ambiguity分析依賴完整集合。
- **驗證**：`AAAA,k=4` 預期 `h*=3`、24 sets；seeded k≤5要求0 mismatch。

### 2026-07-18 02:26 — 甲基只作獨立 tie-break endpoint

- **決定**：本輪只回報
  `Δ_M = D_M(01,11) − D_M(10,11)`，不先與sSNV likelihood以任意lambda相加。
- **理由**：同一V的sSNV likelihood理論上完全同分；lambda若用同資料調整會製造結果。
- **驗證**：positive-10、positive-01、tie、uninformative、missing-CpG synthetic controls。

## 偏離

### 2026-07-18 02:26 — 真實甲基正例不可評估

- **原計畫**：用HCC1395_DORADO chr6找一個10／01／11 real edge unit。
- **實際**：MINREAD 1、2、3、5皆為0個同component×HP三state單元。
- **處理**：真實部分降為feasibility negative；不從其他不相容輸出拼接正例。

## 折衷

### 2026-07-18 02:26 — 新 backend版本

- **折衷**：`highspy`目前不存在；允許在`/tmp`安裝pin版本做API prototype。
- **代價**：它與SciPy 1.13.1內嵌HiGHS版本不同，runtime不能作單純同版本比較。
- **保護**：只以solution-set一致性作primary gate；速度只標exploratory。

## 未決

- ~~Persistent backend是否能在`AAAA,k=4`把structural builds從25降為1，且保持24-set digest一致？~~  
  **已答**：可以，25→1，24-set digest一致；但wall time未穩定改善。
- ~~B&B在partial constraints強的小型real unit是否比rebuild MILP少wall time／少states？~~  
  **已答（pilot限定）**：H2009 M31為0.232s vs 16.809s；COLO829 M31為0.143s vs 8.949s，完整sets一致。
- ~~現行ranker是否所有primary與sensitivity fits都遵守「per V一次」，抑或fixed-error/bootstrap有意重算？~~  
  **已答**：primary fit每V一次；fixed-error grid與bootstrap因模型／資料不同而有意重算。
- 何處能取得clean CN、同HP×PS、每state足量且有MM/ML的real edge unit？

## Lore / Why not

- 不只比較`h*`：不同optimal families會改變後續候選與ambiguity比例。
- 不把parent assignments逐一重算read-AF：snapshot read states只依V，不依E。
- 不把甲基相似度稱為ancestry：甲基可逆、可收斂，也會受CNA/LOH與technical geometry影響。

## Probe 結果

### 2026-07-18 02:57 — exact solver

- 9個新tests與23個既有hypercube tests全PASS。
- `AAAA,k=4`：`h*=3`、24個最佳sets；current SciPy 25 builds，persistent HiGHS 1 build／25 solve calls。
- H2009 M31：242個最佳sets，三backend完整集合一致。
- COLO829 M31：216個最佳sets，三backend完整集合一致。
- B&B單次pilot分別比current SciPy快72.5×與62.5×；這是兩個偏好案例，不是全資料runtime估計。
- Persistent單次wall time一快一慢；只證明API與build reuse可行，尚未證明值得替換。

### 2026-07-18 02:57 — read-AF cache

- 同V、2個parent assignments：`T=2,V=1,primary_fit_calls=1`。
- 兩個不同V：`V=2,primary_fit_calls=2`。
- 結論：第3項已存在，不需重新實作；只需在未來refactor保留回歸測試。

### 2026-07-18 02:57 — methyl edge

- 7/7 synthetic controls PASS：正向、鏡像、兩種tie、低信心、無共同CpG、混HP×PS。
- HCC1395_DORADO chr6：threshold 1/2/3/5的primary k=2 components分別229/204/191/168；三個exact states各≥1者均為0。
- 結論：scorer程式邏輯可行；real edge仍`ABSTAIN_NO_EVALUABLE_REAL_EDGE_UNIT`。

## 2026-07-18 09:11 — Optimized backend bounded implementation

### [決策] objective 與 family 證書拆開

- **決定**：reduced terminal 數 `q<=8` 時，先以 group-terminal subset DP 一次證明
  `h*`；再把 certified `h*` 傳給 bitset obligation B&B 列完整 family。
- **理由**：舊 one-pass B&B 在 `max_sets/deadline` 前可能只取得非最優 incumbent；
  `objective` 欄名不足以表示是否已證明。新路由分成
  `objective_certified`、`incumbent_objective`、`family_complete`。
- **反例**：`full=RAA, partial=ARX/AXA` 的舊 cap-first traversal可先看到 cost 3，
  真值為 `h*=2`；新路由固定回 `h*=2`，且唯一 optimal family可完整結束。
- **Fail-closed**：`AAAA,k=4,max_sets=1` 在已有 certified `h*` 時回
  `CANDIDATE_SET_INCOMPLETE_MAX_SETS`；total deadline 在尚無 feasible
  incumbent 時回 `NO_FEASIBLE_CERTIFICATE_DEADLINE`。兩者皆
  `ranking_allowed=false`，不輸出 winner。

### [決策] exact-preserving bitset／antichain／downward closure

- **決定**：selected、excluded、groups、predecessors皆改為Python integer bitsets；
  raw mutation state與dense index保留雙向 mapping。
- **決定**：每個search state重新計算動態 obligation antichain；只移除包含較小
  obligation domain的superset，不移除subset。
- **決定**：universe只保留mandatory／active group states的Boolean ancestor closure。
- **驗證**：raw `k=12`、active bits `{0,11}` 保留state
  `[0,1,2048,2049]`，objective與2組完整family digest皆與brute force一致。

### [決策] fixed-N edge 與biology邊界

- fixed node set 的edge-local additive score改以每個child獨立parent mapping，
  一次計算best score、best-tree count、tie與softmax posterior；不materialize
  parent Cartesian product。
- 無read evidence的unary connector只產生
  `MULTI_MUTATION_EDGE_EQUIVALENCE` projection；mutation order標
  `UNRESOLVED_NO_READ_EVIDENCE`，原始candidate與objective仍保留。
- `M1_STRICT_INFINITE_SITES`只有 allele-specific `total/major/minor CN=2/1/1`，
  clonal/subclonal LOH、deletion、amplification、CNA、WGD、mutation loss 與
  duplicated-copy evidence 全為 absent，且 phasing／mapping／BQ／strand／
  read-independence／somatic-variant QC 全通過時才執行；缺欄位或矛盾 adverse
  evidence 一律 `ABSTAIN_CN_LOH_GATE`。
- `M2_LOSS_SUPPORTED_DOLLO`尚未建立copy/loss state model，固定回
  `UNRESOLVED_LOSS_SUPPORTED_MODEL_NOT_IMPLEMENTED`，不可發表biology claim。

### [偏離] 未修改 canonical／production router

- **原要求可解讀為**直接替換current solver；但既有pre-decision audit與最新研究計畫
  都是`PROBE_ADAPTIVE_RELAXATION`，明確要求先做P1–P3a bounded dual validation。
- **實際處理**：新增
  `InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/optimized_hypercube_backend.py`，
  保持
  `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`
  byte-identical。
- **Promotion blocker**：尚未跑33-tail、H2009/H1437 stress panel與正式route census；
  因此不能啟動full run或宣稱production improvement。

### [偏離] stale build-count assertion

- current solver已完成prepared-base remediation，`AAAA`的base build實際為1，
  不是歷史報告的25。
- 測試已改為驗證current/persistent皆build一次；科學gate仍是`h*`與完整family digest，
  build count只作diagnostic。

### 實際輸入、命令與輸出

**輸入**：

- `InterSubMod/research/20260718_solver_methyl_edge_probe/tests/fixtures/real_units.json`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`
- `InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/optimized_hypercube_backend.py`

**有限域oracle命令**：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/verify_optimized_backend_oracles.py \
  --output InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/oracle_receipt.json \
  --seeded-k4-cases 300
```

**實際輸出**：exit 0；`27,054`個exhaustive `k<=3` cases＋`300`個seeded
`k=4`＋1個noncontiguous `k=12` case，合計`27,355`，mismatch=`0`；
cap/deadline ranked count=`0`。輸出：
`InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/oracle_receipt.json`。

**同輸入benchmark命令**：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/run_optimized_backend_benchmark.py \
  --output-dir InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/runs \
  --receipt InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/benchmark_receipt.json \
  --import-current-dir /tmp/intersubmod_solver_baseline_20260718 \
  --repeats 5 --warmups 1 --time-limit 30 --q-max 8
```

**實際輸出**：三案皆5/5 repeat digest穩定且current/new完整family完全相同：

| Case | `h*` | `V` | current median solver s | optimized median solver s | bounded speedup | optimized/current RSS |
|---|---:|---:|---:|---:|---:|---:|
| AAAA | 3 | 24 | 0.291476 | 0.010585 | 27.54× | 0.285 |
| H2009_M31 | 5 | 242 | 22.686913 | 0.276785 | 81.97× | 0.251 |
| COLO829_M31 | 5 | 216 | 12.280389 | 0.157009 | 78.21× | 0.254 |

以上是同機bounded cases，不是33-tail或production speedup。Current為per-MILP-solve
30秒上限，新 backend 為 prepare＋subset DP＋bitset B&B 共用的 per-unit 30秒上限；
所有案例均在限制前完整，receipt明示 `time_limit_semantics_same=false`，不可解讀成
相同production SLA。

**完整probe tests**：

```bash
INTERSUBMOD_HIGHSPY_PATH=/tmp/intersubmod_highspy_1_15_1 \
python3 -m unittest discover \
  -s InterSubMod/research/20260718_solver_methyl_edge_probe/tests \
  -p 'test_*.py' -v
```

實際輸出：33/33 PASS，exit 0，wall 12.348s。Current hypercube regression另為
27/27 PASS，exit 0，wall 1.643s。

### [驗證] red-team 與獨立重算

- 先前 red-team 指出的 M1 fail-open／selected-state不連root、DP allocation/resource
  gate、final-mask／final-leaf／NaN total deadline、無-incumbent狀態、edge
  `tie_tolerance`／單edge與跨child累加overflow，以及 unary multiplicity 均已補上
  regression。反例重播為 `BLOCKER=NONE`。
- 獨立 AI 在 `/tmp` 重新跑完整 `27,355` oracle cases：mismatch=`0`、
  incomplete-ranked=`0`；另逐欄重算30次raw runs、126個summary scalars與9個效能比，
  mismatch=`0`。
- 最新 bounded receipts：
  - oracle SHA-256：
    `52fb731797e6c0ab5fe327d7316f3c7352050110b35beb5713fb6120a8688a25`
  - benchmark SHA-256：
    `09b80734619c1dc02e9c9b08b8b167979845ba09ab49b2e419e4eb5b8d0bc7cb`
- **判定**：`PASS_FOR_BOUNDED_DUAL_PILOT`；`canonical_or_production_promotion_allowed=false`。

### [驗證] standalone HTML 與瀏覽器 QA

- **輸入**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/{benchmark_receipt,oracle_receipt}.json`
- **輸出**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/20260718_HypercubeOptimizedBackend完整流程與效能驗證_01.html`
- **瀏覽器 QA**：Chromium headless；desktop 1440×1200＋mobile 390×844；18/18 PASS。
  已檢查 title、單一 H1、10 sections、scope ribbon、embedded certificates、TOC targets、
  theme toggle、details controls、desktop/mobile body overflow、console/page errors與實際截圖。
- **QA receipt**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/html_qa/browser_qa_receipt.json`
- **Frozen policy**：HTML、QA、current/optimized raw runs、receipts、sources與tests將由
  `freeze_optimized_backend_evidence.py`以no-clobber＋`0444`＋manifest sidecars一次封存；
  bounded bundle不等於production promotion。

### [未決] promotion 前必做

1. 以不可變source/input snapshot重播33-tail與H2009/H1437 stress panel。
2. 量測總wall、p95/p99、peak RSS與complete率；完整controls digest mismatch必須為0。
3. single-terminal factorial compact-family與complete cache/checkpoint是共享workspace
   最新plan新增的優先工作，本次未冒充已完成。
4. 只有上述gate通過，才可提出feature-flag router；本輪禁止替換canonical或啟動full。

## 2026-07-19 — [決策][驗證][未決] Frozen hard-tail harness 與資源互斥

- **服務目標**：G3（保留read/subcube exact estimand）、G4（hard-tail panel可重播）、G5（加速主張有不可變input與fail-closed receipt）。
- **Pre-decision amendment**：仍為`PROBE_STRESS`；只核准frozen-pattern replay與隔離dual harness，不核准canonical替換／full 154。Primary固定minread=1的44 units／25 structural keys；sensitivity固定minread=2的36／21。
- **新增程式**：`InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/export_solver_stress_panel.py`與`run_solver_stress_panel.py`；測試為`InterSubMod/research/20260718_solver_methyl_edge_probe/tests/test_solver_stress_panel.py`。Runner使用單一per-unit總deadline、可接受`max_sets=none`，任何deadline/cap均`family_complete=false`、`ranking_allowed=false`、`incomplete_ranked=false`；不materialize parent trees。
- **Frozen input/output**：權威pointer為`InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_stress_panel_v1/AUTHORITATIVE_MANIFEST.json`；manifest為`.../authoritative_r2/manifest.json`，content SHA-256=`013b784bcd9e0f29b6f478a7ae92aef42357255836b6efb44657180cf8f4b7d1`，byte sidecar重驗PASS。三個source gzip均截斷，但selected units缺失=0；只稱selected complete-record pattern replay，不稱全chr6完整。
- **Objective authority hotfix**：optimized final objective以`family_certificate.objective/objective_certified`為權威；subset-DP若因q/resource gate未證明，僅作diagnostic。Regression確認DP未認證時，complete B&B仍正確回objective=3。
- **測試命令／actual**：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest discover -v -s research/20260718_solver_methyl_edge_probe/tests -p 'test_*.py'`；主代理重跑40/40 PASS，12.615 s，exit 0。
- **Dry-run command／actual**：`.../python3 export_solver_stress_panel.py --census research/20260718_solver_methyl_edge_probe/results/m2_solver_route_census_receipt.json --dry-run`；`all_checks_pass=true`、primary=`44/25`、sensitivity=`36/21`、union keys=`46`、exit 0。
- **Bounded smoke**：同一`k8_m5` control current／optimized皆`h*=5,V=120`、family digest=`1dbb67a5575a50ebb96d3f6238a9bf5f181a1880c4173ad7c645b1fcbaae62e5`；solver wall=`3.818275/0.050356 s`（75.82×）。外層process記錄約`4.78/0.29 s`；此時系統有41-process衝突，只作correctness smoke，不作正式speed gate。
- **Factorial oracle**：optimized unbounded materialization對三案回`720/5,040/720`，solver wall=`0.4978/5.6765/0.5138 s`，全部objective certified、family complete、oracle PASS；`V=5,040`在`max_sets=256`時正確回`CANDIDATE_SET_INCOMPLETE_MAX_SETS`且禁止ranking。
- **未決／禁止外推**：現場all-sSNV cooccurrence audit仍為1 root／41 processes、約3900% aggregate CPU；完整25-key dual timing未執行，早期optimized-only partial run已標`PARTIAL_RUN.json`，不可作promotion evidence。Zero-conflict後才執行`runs/primary_dual_unbounded_r1`，要求control mismatch=0、total wall≥5×、optimized p95不劣、control RSS ratio≤1.5、optimized全部complete、incomplete-ranked=0。
- **仍缺**：H2009／H1437尚無M2 frozen stress units；single-terminal compact descriptor只能exact回h*/V/shape，不能替代現行逐vertex-set read-pattern likelihood winner。故尚不可宣稱「所有區域排名完成」或「production整體加速」。

### 2026-07-19 — [驗證][決策] \(h^*+d=m\) exact perfect-family counter

- **根因定位**：R2 partial diagnostics中四個optimized hard keys皆已`objective_certified=true`，但30秒未列完family；三案符合`h^*=m`、一案為`h^*=m+1`。問題是family輸出量，不是最低解未證明。
- **一般式**：令`d`為distinct mandatory non-root full states數；selected tree的非根節點／邊數為`h*+d`。只有`h*+d=effective m`時，所有active mutation剛好各取得一次，才可使用exact perfect-event counter；不是無條件把recurrence-allowed模型改成strict model。
- **隔離實作**：新增`InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/perfect_family_count.py`、`run_perfect_family_count_validation.py`與`tests/test_perfect_family_count.py`；使用mutation-event subset DP，時間`O(3^m poly(m))`、空間`O(2^m+m^2)`，只count／certificate，不materialize family、不做VAF ranking。
- **Hard results**：a094=`104,640`（約0.0265s）、09b=`122,281,152`（約0.0613s）、cec=`27,360`（約0.0257s）；e5b之perfect count=`0`且objective evidence=`m+1`，正確回`ABSTAIN_UNSUPPORTED_RECURRENCE_COUNT`。
- **外部controls**：8個source-complete controls中，7個符合perfect gate者exact count逐一吻合；1個objective=`m+1`的136-family control正確abstain。46/46 frozen keys皆有total status，failure=0。
- **主代理tests**：完整隔離probe 57/57 PASS，14.332s，exit0；receipt與sidecar mode均`0444`、hardlink count=1、`sha256sum -c` PASS。
- **權威receipt修正**：v1 receipt綁R2 manifest，保留不可覆寫；主代理以相同source重新執行並建立`InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_stress_panel_v1/perfect_family_count_v2/receipt.json`，綁定R3 manifest SHA=`6c59345c7657f19136f5f8b855ec0ffa2aacd525422a697a0ea36b416c1f74cc`／content SHA=`5daac21fac8b5e4a617e700d61fd2fb87afe708f3a0d039e31fd6a3d31bb1fa1`；v1/v2的summary、hard counts、controls、method contract完全相同。
- **方法結論**：09b有1.22億個同樣最簡vertex sets，顯式TSV＋逐候選read-AF fit不是可行的「全部完成」定義。下一步必須用DP traceback DAG/ZDD表示logical family，並把`logical_family_complete`與`ranking_complete`分開；e5b另需`excess=1` recurrence-aware DP，現階段不可排名。

### 2026-07-19 — [驗證][決策] H1437／H2009 historical-v5跨樣本壓力面板

- **適用範圍**：`PROBE_CORE_ONLY`；輸入是historical-v5
  `subread_groups_by_hp`與既有layered solver狀態，不是M2 component／phase-set
  extraction，也不是end-to-end或生物clone驗證。
- **Frozen panel**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/cross_sample_solver_stress_panel_v1/manifest.json`。
  H1437有1,773個historical incomplete units／1,763 structural keys；H2009為
  3,759／3,746；union=5,506。選取全部48個`q>8` keys與每個非tail
  sample×k×m×q×cap stratum的字典序最小key，得到208個distinct keys
  （H1437=86、H2009=122）。
- **資料綁定**：12/12 data source SHA、兩份sample output manifests與origin
  manifest皆重驗PASS；missing pattern joins=0。Key-list SHA-256=
  `7ae4f85ba9f05f146303d3f7921e5f333f395c690cdcd1dbd84ea39037b75692`，
  payload-list SHA-256=
  `30a1a767f2bdf6c83b30c383d941748bd933207bc93cd80a4b8a49ca0acbd682`。
- **主代理驗證**：完整isolated probe 65/65 PASS，12.846s；manifest byte
  sidecar、semantic digest與12/12 source hashes獨立重算PASS。

### 2026-07-19 — [驗證][限制] 跨樣本perfect-family count census

- **輸入／輸出**：以上208-key frozen manifest →
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/cross_sample_solver_stress_panel_v1/perfect_count_census_v1/receipt.json`。
- **結果**：208/208均得到finite total status；186個recurrence-free exact
  （root-only 167、mandatory-state 19），22個明確
  `ABSTAIN_UNSUPPORTED_RECURRENCE_COUNT`。最大exact family count為82,396。
  H1437 exact/abstain=82/4；H2009=104/18。
- **時間口徑**：runner wall=1.087669s、counter sum=0.447240s、
  per-key max=0.006808s；當時另有正式hard25 timing，因此只作diagnostic，
  `formal_speed_claim=false`。
- **不可外推**：此端點是logical perfect-family count certificate，不是family
  materialization、current read-AF/VAF ranking或recurrence-family completion；
  receipt固定`ranking_allowed=false`且未執行current／optimized explicit solver。

### 2026-07-19 — [驗證][決策][未決] Compressed VAF/read-AF exact lazy probe

- **Endpoint**：只維持current M2 primary BQ-aware vertex-set profile mixture
  likelihood與`tie_tolerance` membership；不包含fixed-error、bootstrap、
  responsibilities或all-candidate TSV。歷史ancestor-AF edge heuristic不是同一分數。
- **隔離實作**：新增
  `InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/compressed_vaf_rank_probe.py`、
  `run_compressed_vaf_rank_probe.py`與
  `tests/test_compressed_vaf_rank_probe.py`。DP traceback保留logical family，
  branch只在convex-hull/KKT-certified UB低於incumbent減tie tolerance時剪除。
- **Exact oracle**：m=1…5的family、best score、winner IDs與完整tie count對獨立
  parent-vector exhaustive oracle全部零mismatch；對稱2-chain正確保留2-way tie。
  Focused tests 8/8 PASS（3.862s test、4.204s wall）。m=5 synthetic family
  `120`中實評`24`、安全剪除`96`，只是5× candidate-evaluation reduction，
  不作formal wall-speed claim。
- **Fail closed**：candidate cap、deadline、tie-output cap、optimizer非certified、
  UB非finite或traceback count不守恆皆不得發布winner；diagnostic incumbent固定
  `authoritative=false`。
- **不可變receipt**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v1/receipt.json`；
  runner wall=`3.331930s`、small oracle PASS，receipt/sidecar均`0444`、
  nlink=`1`、`sha256sum -c` PASS。
- **未決／限制**：formal hard25 timing持有資源鎖，09b bounded probe未執行並
  明確記為`NOT_RUN_RESOURCE_MUTEX`。最壞全部1.22億候選同分時，完整顯式tie
  membership仍有`\Omega(T)`輸出下界；沒有額外exact compression theorem前，
  不得宣稱任意hard family可多項式時間完整排名。
- **Review amendment（未執行）**：prototype目前只對top-level traceback
  branches計算一次UB，不是recursive DP-state B&B；正確性條件不變，但剪枝力
  可能不足。原v1 receipt保留為immutable historical small-oracle evidence，
  不得拿來代表current source或hard結果。後續runner預設改寫全新
  `compressed_vaf_rank_probe_v2/receipt.json`，執行時重新綁current prototype、
  counter與R3 authority validator hashes。
- **R3 hard-loader修補（靜態完成、未測）**：改為fail-closed驗證pointer
  byte sidecar、schema/status、relative path containment、manifest byte SHA與
  sidecar、semantic content digest、`AUTHORITATIVE_R3`、全部frozen checks、
  case/unit structural records及frozen pattern SHA。依資源鎖要求，此修補尚未
  執行任何test或hard command，必須等formal timing結束後才驗證。

### 2026-07-19 — [驗證][決策][未決] Compressed ranking P1 fail-closed

- **Claim修正**：普通NumPy/SciPy float UB沒有directed rounding，固定
  `numerical_bound_certified=false`。若`pruned>0`，core回
  `INCOMPLETE_UNCERTIFIED_FLOAT_BOUND_PRUNING`且隱藏best/ties；只有全family
  實評且tie輸出完整才`exact_publishable=true`。v1的machine-exact pruning claim
  已撤回。
- **Focused tests**：14/14 PASS，5.953s，exit0；新增branching、mandatory full、
  mixed R/A/X、gapped active fixtures，逐branch驗互斥、count與
  `possible_vertices` containment；另驗strict controls、LL/gap finite/nonnegative、
  cap/tie fail-closed。
- **Small oracle**：m=5 full evaluation=`120/120`、prune=0，best score與24-way
  tie IDs對exhaustive current endpoint完全相同；float diagnostic=`24 evaluated /
  96 pruned`，但正確`ranking_complete=false`，不得當加速後exact結果。
- **v2 hardened receipt**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v2/p1_hardened_r2/receipt.json`；
  wall=`5.764620s`、hard=`NOT_RUN_P1_REVIEW_GATE`，receipt/sidecar均0444、
  nlink=1、SHA sidecar PASS；綁`hypercube_exact.py`與NumPy 1.23.5，JSON
  `allow_nan=false`。
- **複雜度**：目前是top-level branch materialization，含所有
  `possible_vertices`最壞約`O(m*3^m)`，不是recursive B&B；無machine-certified
  pruning時完整排名最壞仍`Theta(F*Fit)`，完整tie輸出仍`Omega(T)`。
- **P2／未決**：hard unit唯一ID與unit↔case structural逐欄一致尚未補完，更強
  same-FD/inode TOCTOU證明亦未完成；因此09b hard繼續NOT_RUN。Core deadline
  不涵蓋DP／branch建構／單次optimizer，不得稱5秒hard bound。

### 2026-07-19 — [驗證][修補] Perfect counter authority／evidence fail-closed

- 獨立red-team以2,331個full／partial／gapped-active cases、358個canonical
  brute-force cases與全46-key重算確認計數核心P0=0；hard counts
  `104,640／122,281,152／27,360`與recurrence abstention皆一致。
- P1修補後，validation runner優先解析R3 authority pointer；direct R2、pointer
  escape、sidecar／byte／semantic tamper與stale evidence皆拒絕。歷史optimized
  evidence只有在source manifest與R3 case的structural key、完整payload、k與m
  全相同時，才標`VERIFIED_STRUCTURAL_EQUIVALENCE`。
- Focused 19/19與當時完整probe 79/79均PASS；舊v2 receipt未覆寫。因source
  已再修正`h*=m-d`命名與嚴格整數型別，新權威receipt必須另建v3，不能沿用
  舊source hash。

### 2026-07-19 — [驗證][NO-GO] Resource-gated hard25 dual run

- **正式資源閘門**：all-sSNV 40-worker audit自然結束後，v4 bootstrap preflight
  exit 0；active conflict process/root=`0/0`，available bytes=
  `767,608,098,816`，required reserve=`322,122,547,200`，且未建立formal
  M2 root。Preflight attempt：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260716_read_linked_hypercube_exact_likelihood_validation/m2_frozen_release_v4_attempts/20260719T035234.620261042Z-765565-6313/`。
- **輸入／命令**：R3 authority pointer的primary 25 structural keys；current與
  optimized各1次，per-unit total deadline=30s、`max_sets=none`、`q_max=9`。
  完整命令與GNU time保存在
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_stress_panel_v1/runs/authoritative_r3_primary_all_both_r1.time.txt`。
- **完整執行狀態**：50/50 worker rows落盤；current complete=3/25、
  optimized complete=16/25；current／optimized incomplete=22/9；任何
  incomplete-ranked=`0`；controls objective／family digest mismatch=`0`；
  factorial checks=4/4 PASS。
- **同端點可比子集**：只有3 keys雙方都完整且objective/digest相同。
  Solver speedup min/median/max=`63.17×/90.55×/263.04×`；
  isolated subprocess speedup=`16.39×/42.11×/69.49×`。
  這只是paired-complete子集，不能當25-key或production整體speedup。
- **整體不可比較**：兩backend皆有deadline-censored runs，因此receipt正確回
  `reported_ratio_kind=NOT_COMPARABLE`、ratio=`null`。Observed censored wall
  current/optimized=`703.77/327.54s`不能轉成completion speedup；p95 wall
  `30.579/30.221s`，RSS control ratio=`0.2491`。
- **9個optimized incomplete根因分層**：perfect counter為6個exact perfect
  families（counts=`104,640；122,281,152；27,360；72,576；2,680；9,408`）
  及3個recurrence-required abstentions。顯式family輸出與recurrence-aware
  family尚未解決。
- **Final gate**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_stress_panel_v1/runs/authoritative_r3_primary_all_both_r1/receipt.json`
  為`FAIL_CLOSED`；50 rows存在、p95與RSS通過，但
  `optimized_all_selected_families_complete=false`、controls complete=false、
  overall performance ratio不可比較。因此禁止production adapter、雙actual-data
  pilots與full 154啟動。Receipt／time與sidecars已重驗並封成`0444`。
- **獨立audit與artifact hygiene**：50/50 raw JSON與receipt fields逐欄吻合，
  stdout JSON 50/50一致、stderr全空、46/46 R3 structural hashes一致。第一次
  錯誤建立的0-byte sidecar已保留改名為
  `time.txt.sha256.INVALID_EMPTY_DO_NOT_USE`並設`0444`；正確time sidecar仍在
  runs父目錄且PASS。R3 pointer／manifest／sidecars及本次150個raw files合計
  154檔均已封為`0444`、nlink=1，未改任何內容。

### 2026-07-19 — [驗證] Current-source exact count receipts

- HCC hard-tail v3：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_stress_panel_v1/perfect_family_count_v3/receipt.json`；
  focused 22/22 PASS，46/46 total、7 exact control matches＋1 expected
  recurrence abstention，R3 pointer與4個hard evidence皆通過
  `VERIFIED_STRUCTURAL_EQUIVALENCE`。Receipt SHA-256=
  `7f19d7940a4299a4a2d9b5429515132d610d1bcf68a8be90ce6fb85846bbebea`。
- Cross-sample v2：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/cross_sample_solver_stress_panel_v1/perfect_count_census_v2/receipt.json`；
  208/208 total，186 exact＋22 abstain，current live counter SHA=
  `721f24bb2aee270d2e242fe65c315b56590fd6570518a7d9214240191c090520`。
  Receipt與sidecar已重驗並封成`0444`。

### 2026-07-19 — [驗證][最終決策] 全流程完成度與效能閘門

- **主代理獨立重跑**：compressed focused 14/14 PASS（5.796s）；
  整個isolated solver probe suite 96/96 PASS（17.832s）。
- **Main-agent receipt**：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v2/main_reverify_r3/receipt.json`；
  runner elapsed=`5.609974s`、small oracle PASS、hard=
  `NOT_RUN_P1_REVIEW_GATE`、receipt/sidecar `0444`、nlink=1、
  `sha256sum -c` PASS。
- **m=5 machine-safe結果**：full evaluation=`120/120`且對exhaustive oracle；
  float diagnostic=`24 evaluated/96 pruned`，但
  `numerical_bound_certified=false`、`ranking_complete=false`並隱藏winner。
- **完成度結論**：hard25 controller 50/50 rows均完成回報，但family completion
  僅current=`3/25`、optimized=`16/25`；因此不能稱「所有候選與VAF ranking
  都完成」。
- **效能結論**：3個paired-complete cases的solver speedup為
  `63.17×–263.04×`（median `90.55×`），只作診斷。含censoring的表面
  `703.769/327.540=2.149×`不是completion speedup；正式整體ratio維持null。
- **Final gate**：`NO-GO_PRODUCTION`。依預註冊條件禁止production adapter、
  雙actual-data pilots與full 154。完整說明：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/20260719_全流程完成度與效能閘門最終驗證_01.md`。

### 2026-07-19 — [交接] AI 單一入口

- 已建立：
  `InterSubMod/research/20260718_solver_methyl_edge_probe/20260719_最簡演化樹與VAF排序AI交接重點_01.md`。
- 文件包含背景摘要、術語、數學邊界、已證明／未證明、程式碼地圖、
  權威receipts、重驗命令、下一步P0–P5、禁止事項、可直接貼用的AI prompt與
  接手驗收清單。
- 交接狀態保持`NO-GO_PRODUCTION`；本文件不改變既有promotion gate，也不把
  historical-v5 solver-core證據升級為M2 end-to-end證據。
- **保存風險**：本cycle目錄目前整體為Git untracked，`git ls-files`計數為0；
  已在交接文件標示commit不可還原此bundle，並要求以receipt-specific
  `source_bindings`與SHA sidecars為準。`docs/CURRENT_FOCUS.md`及
  `docs/experiments/INDEX.md`亦尚未收錄本輪Final Gate，冷啟動AI必須直接讀
  交接文件。
