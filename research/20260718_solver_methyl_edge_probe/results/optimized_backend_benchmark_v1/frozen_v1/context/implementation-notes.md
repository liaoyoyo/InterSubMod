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
