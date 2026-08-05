<!--
建立時間: 2026-07-24
目標: 逐步記錄 exact-PS C++ topology/read-AF 實作的設計、偏離、折衷與未決
處理範圍: Task Type B / comprehensive validation / HCC1395 gate then all seven datasets
關聯檔案: InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/pre-decision-audit.md
-->

# Exact-PS C++ topology + read-AF implementation notes

## 設計決定

- [2026-07-24] 使用 Stage-Gate：HCC1395 exact parity與加速通過前，不啟動七資料集。
- [2026-07-24] 模型凍結為 recurrence-allowed、unit-edge-cost、single-root、
  upward Hamming-1 directed Boolean hypercube；不暗中改成 strict infinite-sites。
- [2026-07-24] Python `Fraction` AF score是語意 oracle，不是效能基準的唯一後端。
- [2026-07-24] C++ 可對固定 vertex set 直接計算最佳 score、最佳樹並列數與唯一性；
  預設不 materialize所有 parent trees。只有使用者要求完整樹或 shape census需要時，
  才展開 tied-best parent assignments，並受獨立 output cap約束。
- [2026-07-24] objective certificate、optimal vertex-set family certificate、
  AF ranking certificate分離；任何後層不得補救前層 incomplete。
- [2026-07-24] 七資料集輸出只保存壓縮 unit records與必要代表樹，不複製大型輸入。

## 偏離

- 尚未將 AF 分數稱為 likelihood；舊公式沒有 error model、CN、LOH、purity或
  clone prevalence normalization，只能稱 exact ordering score。
- HCC歷史 7/23 topology的420個 capped units不是 oracle winner；只能作
  `historically_incomplete` 對照。

## 折衷

- 為達 Python exact tie parity，C++ 使用 arbitrary-precision rational，而非 double。
- shape 若只需判斷唯一性，可在 compressed family上計算；若需完整 shape census，
  只展開全域 tied-best subset，不展開所有低分 parent trees。
- k<=12限制 vertex universe最多4096，但不保證 optimal vertex-set數或時間有固定上限；
  仍保留 per-unit deadline、memory/output budget與 fail-closed。

## 未決

- 已決：legacy parity 與本輪 production 均使用 MLHP 內
  exact-PS supported-pattern HP-specific `R/(R+A)`、`A/(R+A)` 分母；
  不稱 caller VAF 或 CCF。
- 已決：research CLI 連結 LongLineage exact obligation B&B 與 parent-mapping
  backend；不修改 production target。
- 已決：154/154 extraction receipts與154/154 strict receipts通過；七樣本
  input readiness成立。
- 未決：按需輸出全部 tied trees 的 enumerator尚未實作；目前只輸出 exact
  tie count與一棵 deterministic representative。
- 未決：guard=1,000下10,717個 mutation-bearing units保守棄權；k>=6長尾
  需要更強的 exact pruning或明示成本的adaptive second pass。

## 執行紀錄

### 2026-07-24｜環境與前置狀態

- Repo commit：`0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- HCC input authority：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2/`
- HCC historical Python topology：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/`
- 可用 `/big7_disk` 空間：約 675 GB；採 compact-output與preflight estimate。
- 一個由本輪誤啟動的遞迴 `rg` 高 CPU程序已終止；未觸碰其他工作。

### 2026-07-24｜舊 Python AF 初步公式稽核

- Source：
  `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py`
- Site AF：`ALT/(REF+ALT)`，Python `Fraction`。
- Edge score：child新取得 mutation `j` 時，對 parent已有每個 mutation `i`
  加 `AF_i-AF_j`。
- 初判：對固定 vertex set且recurrence-allowed時為edge-additive，可依 child分解。
- 待完成：逐行 source contract、exhaustive property oracle與HCC paired parity。

### 2026-07-24｜舊 Python AF 全量 factorization oracle

- Oracle scope：7/23 HCC exact-PS `display_all`＋相同 MLHP；只讀。
- Complete ambiguous T>1：4,601 units。
- AF 分母可用：4,505；zero denominator：96。
- Explicit trees：26,942；逐樹與 factorized 的 candidate count、top score、
  top tie count均為0 mismatch。
- Provisional exact-PS supported-pattern read-AF：unique first 2,473、
  tied same topology 1,804、tied different topology 228。
- 重要限制：26,942 trees同時也是26,942 distinct vertex sets；本批每個
  recurrence-free vertex set只有一個 parent mapping。因此 HCC實際加速主要仍應
  來自 C++ optimal vertex-set solver與streaming，不能把加速單獨歸因於
  parent-factorization。
- 分母固定為 MINREAD後 supported patterns的HP-specific R/A；X與低支持pattern
  不進分母。正式名稱使用
  `exact-PS supported-pattern read-AF`，不得稱 caller VAF。

### 2026-07-24｜七資料集 L2 input readiness

- 154/154 extraction receipts與154/154 strict receipts：schema、all_pass、
  checks及receipt-sidecar SHA全PASS。
- L2必要檔完整：每樣本44/44 extraction files與66/66 strict files存在且非空。
- 全七資料集 strict W=85,621；k>12 W=5,323；bounded blocks數學下限97,797。
- 最低壓縮讀取量約2.8 GiB；H2009最重（W=24,193、k>12=4,206、
  原始max k=567）。
- HCC 7/22 pilot與全七 strict production分母不得混用：
  - 7/22 raw/exact-PS pilot供 Python oracle parity；
  - 全七採7/23 strict W，HCC正確production partition為11,462 units→
    11,712 blocks，不納入28,384個unlinked singleton memberships。
- HCC production partition逐chr 22/22 Python/C++已PASS，但root缺
  `run_receipt.json`；本輪需重建/驗證root receipt。其餘六樣本尚未partition。
- Strict provenance有兩個builder SHA：前六樣本多數為`7260a7…`，
  HCC1954使用missing-PS hotfix `912721…`；既有獨立D010證明未觸發分支的
  132 chromosomes數值等價。本輪揭露mixed provenance，不冒充單一builder bytes。

### 2026-07-24｜HCC1395 Stage-Gate

- C++ binary：
  `/bip7_disk/liaoyoyo2001/build-exact-af-20260724/bin/exact_ps_topology_af`
  （SHA256 `ba13ccf23d091854c191f81dd97fa891368d11179df9a69e915f0340b7233b2e`）。
- 11,122/11,122 historical complete units：active k、最小 hidden-node 數、
  minimum family count/digest與parent-tree count全部0 mismatch。
- 4,601 read-AF oracle units：exact score與top tie count全部0 mismatch。
- Historical capped 420 units：C++完整補解146、fail-closed abstain274、
  invalid 0。
- Controlled wall time：Python 101.54 s；C++ 6.73 s；工程吞吐約15.09x。
- Gate receipt：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/hcc1395_7_22_pilot_gate_guard1000_v3/HCC1395.cpp_gate.receipt.json`。

### 2026-07-24｜七樣本 chr1–22 完整執行

- Output root：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/`。
- Technical PASS：7/7 samples；cohort receipt SHA與八項守恆全部PASS。
- 98,955 groups；85,941 mutation-bearing；75,224 family complete；
  10,717 `ABSTAIN_RESOURCE_LIMIT`。
- 71,955 AF-ranked：39,648 unique、32,307 tied；3,224 zero denominator；
  45 recurrence-required。
- C++ topology stage 161.76 s；adapter 696.23 s；成功resume 19m51.49s；
  含首次receipt欄位不相容的fail-closed診斷，共20m36.78s。
- H2009占8,457/10,717 abstains；guard=1,000下k>=6長尾大多未完整列完，
  因此`validation_evidence_eligible=false`，不得聲稱全cohort topology-complete。
- Regression：28 passed in 20.55s。
- 完整報告：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/20260724_C++_exact_PS拓撲與read_AF全七樣本驗證_01.md`。
