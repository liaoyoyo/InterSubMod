<!--
建立時間: 2026-07-18 01:00 +08:00
目標: 確認 read-linked hypercube 候選是否需要資源上限、strict perfect phylogeny 是否能安全縮小候選，以及提出可驗證的分階段 exact 流程
處理範圍: Task Type A exploratory pilot；2 個樣本、8 個 deterministic capped HP units、3 個 synthetic controls；不是全樣本驗證
關聯檔案: InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/scripts/run_pilot.py；InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/scripts/verify_receipt.py；InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/results/pilot_receipt.json
-->

# Perfect phylogeny 限制與候選上限：小型 exact pilot

用 SCQA（Situation → Complication → Question → Answer）：先固定現行模型，再分開「資源上限」、「perfect-tree 生物假設」與「exact 完整性」。

> **TL;DR：候選數／總時間必須設工程上限，但達上限後只能標記 incomplete 並停止 read-AF 排名，不能把前 N 個當成全部。Strict perfect phylogeny 可作條件式限制或敏感度分析，不宜直接取代現行 recurrence-allowed 模型。本 pilot 最有效的新做法是 cross-model certificate：6/8 真實案例可由一般模型的最優解直接證明 strict \(h^*\)，不必等待較慢的 strict MILP；但兩個全候選列舉案例在 15 秒內仍只找到至少 202／172 組，皆未完成。**（影響：中；算法判讀信心：高；跨樣本與生物結論信心：低）

> **PARTIAL／EXPLORATORY**：只測 H2009、COLO829 各 4 個 chr1 capped HP units。不能外推全基因組比例、真實 clone/subclone 數或最終拓撲。

本任務服務 G4（可重現、fail-closed 的候選處理）與 G5（可由 receipt 驗證的方法決策）。

## 0. 研究 round 定義與必要背景

### 背景摘要

目前流程先以 HP 與同 read 的 sSNV 共現狀態建立硬限制，再尋找最少額外 mutation states。真正造成長尾的不是單次找一個可行解，而是必須證明 \(h^*\) 並列完所有同分 node sets；先前 HCC1395_DORADO chr6 正式 pilot 已因逐解重建 MILP、缺少 per-unit 總 deadline，在 8 小時上限 fail-closed。

本輪只回答算法決策：資源 cap 如何不破壞 exact 語意、strict perfect-tree 限制是否真的加速，以及能否用較便宜的 certificate 先得到同等數學結論。它不重算全樣本拓撲，也不把 molecule states 翻譯成 cell clones。

### 固定研究模板

| 欄位 | 本輪定義 |
|---|---|
| 本次研究問題 | 是否必須限制候選；strict perfect constraint 能否安全、有效地縮小搜尋；何時仍可證明 exact |
| 本次假設 | H1：global cap 必要但須 fail-closed；H2：strict 可排除 recurrence-compatible sets，但不保證更快；H3：Model A optimum 可提供 strict lower bound |
| 成功條件 | 合成正反例符合預期；real incumbents 通過獨立 checker；任何 timeout 都不標 complete；能形成可重播 receipt |
| 失敗條件 | 把 limit 當 optimal、把前 N 組當全部、strict checker 與 small-k parent oracle 不一致，或修改 canonical output |
| 評估指標 | solver status、objective、dual bound、gap、wall time、optimal-set count、complete flag、checker pass |
| 主要資料與限制 | 2 samples × 4 capped HP units；全部 chr1；非隨機、非全基因組 |
| 本次改動 | 新增隔離的 pilot、receipt verifier、small-k oracle；未改 canonical／production solver |
| 預計輸出 | JSON receipt、可重播 scripts、本 Markdown 方法決策紀錄 |

研究啟動 5 問：Thread D 僅間接相關；不在 Thread B 撤回範圍；KDE 與 VCF caller AF 不適用；本輪讀 read-state JSON，不讀 BAM／VCF；屬短時間 Python/MILP pilot，不觸發 full release 或 NO-GO gate。

### 重要術語

| 術語 | 中文詞義 | 本輪角色 |
|---|---|---|
| node set \(V\) | 被選取的 binary mutation states 集合 | Phase A／B 真正列舉的候選單位 |
| parent tree \(T\) | 對固定 \(V\) 指派每個非 root 節點的 parent | 同一 \(V\) 可能有多個 edge assignments |
| \(h^*\) | 最少額外／hidden node 數 | 第一層 exact objective |
| complete | 已證明沒有其他同分 node set | read-AF 全候選排名的必要條件 |
| incumbent／bound | 目前可行上界／已證明下界 | time limit 時不能合稱 optimum |
| strict perfect-compatible | 每個 mutation 在全樹最多取得一次 | 條件式 graph property，不是 clone 證明 |

### 已知前提與閱讀建議

- 已成立：現行 partial groups 是聯合 group-hit constraints；不是逐 X completion 選第一個成功樹。
- 已成立：canonical v5 是 recurrence-allowed Model A，rooted three-gamete 目前主要作 post-hoc classification。
- 尚未成立：full chr1–22 × 7 technical datasets 的 exact candidate ranking；先前正式 ranking pilot 為 NO-GO。
- 尚未成立：HP molecule-state topology 等於 cell clone/subclone genealogy。

建議補讀：

- [CURRENT_FOCUS](/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
- [Partial-read 方法稽核](/big7_disk/liaoyoyo2001/InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md)
- [2026-07-17 frozen-v2 8 小時 NO-GO 紀錄](/big7_disk/liaoyoyo2001/InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md)

## 1. 先釐清現行模型到底保證什麼

令每個 HP 內有 \(u\) 個有效 sSNV；binary vertex \(v\in\{0,1\}^u\) 表示一種 read/molecule mutation state。`R=0`、`A=1`、`X=未知`。

現行 Model A 使用 \(z_v\in\{0,1\}\) 選節點，最小化額外節點：

\[
h_A^*=\min \sum_{v\notin M} z_v
\]

其中 \(M\) 是 ROOT 與完整 read states。硬限制為：

\[
z_{\mathbf 0}=1,\qquad z_m=1\ \forall m\in M
\]

\[
\sum_{v\in G_g}z_v\ge 1
\quad\text{（每個 partial R/A/X group 至少被一個相容 state 命中）}
\]

\[
z_v\le \sum_{p\in Pred(v)}z_p
\quad\forall v\ne\mathbf 0
\quad\text{（每個節點有一個只少一個 mutation 的已選前驅）}
\]

因此每條 root-to-leaf path 都是 \(0\rightarrow1\)，但**同一 mutation 仍可在不同分枝取得兩次**。這是「單一路徑不回復」的 monotone tree，不等於全樹 infinite-sites／perfect phylogeny。

Canonical v5 solver 已對每個候選 node set 做 rooted three-gamete **事後標記**，但不把它當成最小化時的硬限制；本 pilot 才另外建立 strict Model B。

## 2. Strict perfect 限制

在 root 為全 0 時，對任一 mutation pair \(i,j\)，候選節點中不能同時存在：

\[
10,\quad 01,\quad 11
\]

即

\[
\neg(C_{ij}^{10}\land C_{ij}^{01}\land C_{ij}^{11})
\]

本 pilot 對每組 pair 加入三個 category indicators，並限制三類最多啟用兩類。若通過，代表在目前 binary-state contract 下，每個 mutation 可全樹只取得一次。

對固定 node set \(V\)，現行 parent-tree 數可直接計算：

\[
T(V)=\prod_{v\in V\setminus\{\mathbf 0\}}
\left|Pred(v)\cap V\right|
\]

若 \(V\) 通過 strict three-gamete，任一節點不會同時有兩個已選 unit predecessors，因此在本模型內 \(T(V)=1\)。這可以省去 parent-edge Cartesian-product 展開；但不同 strict node sets 仍可能很多。

## 3. 為什麼候選必須設上限，但不能把上限當答案

有效 state universe 為 \(2^u\)。若可選 extra pool 大小為 \(M\)，逐層 exhaustive search 在 hidden count \(e\) 需面對：

\[
\binom{M}{e}
\]

而「列出全部同分最優解」至少需要 \(\Omega(S)\) 時間與輸出空間，\(S\) 是實際最優 node-set 數；若 \(S\) 很大，任何顯式列舉法都無法繞過 output-size 下限。

因此要分成兩種 cap：

| Cap | 可以做什麼 | 達上限後的正確語意 |
|---|---|---|
| 分析 cap：總秒數、max sets、RAM／bytes | 防止少數 unit 拖垮全 run | `CANDIDATE_SET_INCOMPLETE_*`；不得宣稱第一名 |
| 顯示 cap：只存前 20／32 個例圖 | 控制 HTML／JSON 大小 | 若完整總數、digest 與分析結果另有保存，不影響 exact claim |

`max_sets=256` 或本 pilot 的 `512` 都只是資源政策，不是「最多只可能有 256／512 棵生物樹」。

建議至少輸出下列互斥狀態：

1. `OPTIMAL_VALUE_CERTIFIED`：只證明 \(h^*\)。
2. `CANDIDATE_SET_COMPLETE`：固定 \(h^*\) 後，solver 再證明沒有其他解。
3. `CANDIDATE_SET_INCOMPLETE_MAX_SETS`：只知道至少找到 N 組。
4. `CANDIDATE_SET_INCOMPLETE_DEADLINE`：時間到，總數未知。
5. `FEASIBLE_NOT_PROVEN`：只有 incumbent 與 lower bound，連 \(h^*\) 都未證明。

只有第 2 類可以把 read-AF/VAF 第一名寫成「在全部 minimum candidates 中的第一名」。第 3–5 類必須 abstain。

## 4. 比「直接重跑 strict MILP」更有效的 staged certificate

因為 strict 可行集合 \(F_B\subseteq F_A\)，一定有：

\[
h_A^*\le h_B^*
\]

若 Model A 已證明 \(h_A^*\)，而找到任一 strict-compatible 解 \(V\) 滿足 \(h(V)=h_A^*\)，則：

\[
h_B^*\le h(V)=h_A^*
\]

上下界合併可得：

\[
\boxed{h_A^*=h_B^*}
\]

這個證明不要求 strict solver 自己跑到 gap=0。

建議流程：

1. Model A 先求 \(h_A^*\)。  
   → 驗證：status=`OPTIMAL`、objective=dual bound、gap=0。
2. 對 Model A 回傳的最優 node set 做 \(O(|V|u^2)\) three-gamete 檢查。  
   → 驗證：若 PASS，直接寫入 cross-model strict certificate。
3. 若第一個最優解不 strict，固定 objective=\(h_A^*\) 只解「是否存在 strict 解」。  
   → 驗證：找到一個即完成 strict \(h^*\) 證明；若要說不存在，必須有 `INFEASIBLE` certificate。
4. 只有需要 read-AF 比較所有 strict candidates 時，才啟動 Phase B 全列舉。  
   → 驗證：global deadline／max sets 任一觸發即 incomplete、ranking abstain。
5. Parent tree 只用分析式計數；strict node set 不展開 parent assignments。  
   → 驗證：small-k oracle 與分析式一致。

這個順序不改變 Model A 的 estimand，也不會把 strict 假設偷偷寫入所有區域。

## 5. 小規模真實資料測試

### 5.1 輸入

- Canonical read-only root：  
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`
- Current exact solver：  
  `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`
- Solver SHA-256：  
  `9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`
- Deterministic scope：H2009、COLO829；每個樣本各取 M=31／63／127／255 的第一個 capped primary-HP unit。
- Phase A：每次 solve 8 秒；Phase B：每 case 總計 15 秒、每次 solve 最多 3 秒、max sets=512。

### 5.2 執行命令

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/scripts/run_pilot.py \
  --phase-a-time-limit 8 \
  --phase-b-total-time-limit 15 \
  --phase-b-per-solve-time-limit 3 \
  --max-sets 512 \
  --output \
/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/results/pilot_receipt.json
```

### 5.3 真實案例結果

`inc/LB` 是目前可行解 objective／已證明 lower bound；`OPT` 才是 gap=0。

| 樣本／區域／HP | \(u\)／M | Model A | Strict solver | Cross-cert strict \(h^*\) |
|---|---:|---:|---:|---:|
| H2009 chr1:3624808-3734966 HP2 | 5／31 | OPT 5；0.072 s | OPT 5；0.138 s | 5 |
| H2009 chr1:8025843-8127514 HP1 | 6／63 | OPT 6；0.631 s | OPT 6；0.357 s | 6 |
| H2009 chr1:6704079-6781584 HP1 | 7／127 | OPT 7；1.604 s | limit 7／6；8.006 s | **7** |
| H2009 chr1:2751181-2810948 HP1 | 8／255 | limit 8／4；8.035 s | limit 8／4；8.016 s | 未證明 |
| COLO829 chr1:13182304-13213106 HP2 | 5／31 | OPT 5；0.031 s | OPT 5；0.072 s | 5 |
| COLO829 chr1:38314729-38487374 HP1 | 6／63 | OPT 6；0.547 s | OPT 6；2.432 s | 6 |
| COLO829 chr1:18185468-18318039 HP2 | 7／127 | OPT 7；2.804 s | limit 7／3；8.003 s | **7** |
| COLO829 chr1:96400782-96560535 HP1 | 8／255 | limit 8／3；8.010 s | limit 8／3；8.004 s | 未證明 |

重點：

- Model A：6/8 自身證明最優；2/8 只有 incumbent／bound。
- Strict solver：4/8 自身 gap=0；4/8 達 8 秒上限。
- 但 6 個 Model A 已證明案例所回傳 node set 全部 strict-compatible，因此用上下界夾證後，**6/8 strict \(h^*\) 已證明**。
- 這 6 個案例均沒有證據顯示 strict 改變 minimum hidden objective。
- 兩個 M=255 案例只知有 \(h=8\) 可行解，下界分別為 4 與 3；不可稱最優。

### 5.4 Synthetic 反例

| Case | Model A | Strict | 意義 |
|---|---:|---:|---|
| full `10,01` | \(h^*=0\) | \(h^*=0\) | 正常 sisters，可為 perfect tree |
| full `10,01,11` | \(h^*=0\) | infeasible | 任選 11 的 parent 都會使某 mutation 在兩枝取得 |
| partial `AX,XA` | \(h^*=2\) | \(h^*=2\) | 未知位置要聯合滿足；不等於兩個 cell clones |

### 5.5 Phase B 全候選列舉

| Case | 已找到 strict optimal node sets | 15 秒後 | 可報告的樹數 |
|---|---:|---|---:|
| H2009 M=31 | 202 | incomplete | \(\ge 202\) |
| COLO829 M=31 | 172 | incomplete | \(\ge 172\) |

兩個案例的每個 strict node set 都只有一個 parent tree，但 node-set ties 本身仍未列完。Strict 限制在這個 bounded test 中**沒有解決 all-ties explosion**。

完整 pilot wall time=95.106 秒。Receipt SHA-256：  
`5c5ac26aa8d11eddfdfd0dbeb5293c09e081fabfada5818c5278437cbdd4695c`。

## 6. 獨立驗證

### Receipt invariants

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/scripts/verify_receipt.py \
/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/results/pilot_receipt.json
```

實際輸出：

```text
PASS checks=65
SUMMARY synthetic=3 real=8 both_certified=4 strict_cross_certified=6
real_proven_model_changes=0 phase_b_complete=0/2
```

### Independent small-k math oracle

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/scripts/verify_small_k_oracle.py
```

實際輸出：

```text
PASS subsets=138 root_connected=91 mismatches=0 k_max=3
```

它對 \(k\le3\) 的全部 node subsets 做枚舉；91 個 root-connected sets 中，rooted three-gamete 與「存在全樹每 mutation 只取得一次的 parent assignment」0 mismatch。

### Current exact solver regression

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
  -m unittest tests/test_hypercube_exact.py -v
```

實際輸出：23/23 PASS，1.109 秒。

## 7. 時間與空間複雜度

令 \(g\) 為 partial groups、\(u\) 為有效 mutation 數、\(S\) 為同分最優 node-set 數。

| 階段 | 大致大小／下限 |
|---|---|
| Hypercube state variables | \(O(2^u)\) |
| Model A sparse constraints | closure \(O(2^u)\) rows，加上 \(O(g)\) group rows |
| Strict indicators | \(3\binom{u}{2}=O(u^2)\) variables、\(4\binom{u}{2}=O(u^2)\) rows；matrix nonzeros最壞 \(O(u^2 2^u)\) |
| MILP worst case | 指數時間；實務依 instance／presolve／branching 而異 |
| All-optimal enumeration | 至少 \(\Omega(S)\) 時間與輸出 |
| Three-gamete post-check | \(O(|V|u^2)\) |
| Parent count | \(O(|V|u)\)，不需 materialize trees |

標準 binary perfect-phylogeny existence 與 incomplete directed perfect-phylogeny completion 有線性或近線性演算法；但本問題還包含「每個 partial group 至少命中一個 state」、「最少額外 distinct states」與「列出全部最優 sets」，不能直接把標準演算法的複雜度套過來。[Gusfield 的 classical perfect-phylogeny 結果](https://faculty.engineering.ucdavis.edu/gusfield/publications/)；[Pe’er et al. 的 incomplete directed perfect phylogeny](https://epubs.siam.org/doi/10.1137/S0097539702406510)。

## 8. 生物假設與 CN/LOH gate

Strict perfect compatibility 只是一個 mutation-state graph property；它不自動證明：

- 每個 state 是不同 cell clone；
- HP1／HP2 各自只有一個實體 chromosome copy；
- molecule fraction 等於 cell fraction；
- 真實歷史沒有 deletion、LOH、subclonal CNA、allele-specific amplification、mutation loss 或 convergent mutation。

尤其 HP 是 germline-haplotype read grouping。CNA 後同一 HP-derived allele 可有多份 copies；同一細胞即可在同一 HP 產生不同 molecule states。LOH／deletion 也可讓既有 mutation 消失，違反 strict infinite-sites。腫瘤演化工具本身也常需處理 ISA violations；PhISCS 明確把 LOH、deletion 與 convergent evolution列為 perfect-tree 假設可能失效的原因，並以 ILP/CSP 允許 subperfect 解。[PhISCS](https://pubmed.ncbi.nlm.nih.gov/31628256/)。

因此 strict 模式建議只在下列條件成立時作較強解釋：

1. allele-specific CN 與 LOH 已知且相容；
2. 無 mutation-loss／recurrence 證據；
3. mapping、phasing、read-error QC 通過；
4. 結論仍寫成「perfect-compatible molecule-state topology」，不直接寫成 clone truth。

Branch-and-bound 可作效能方向；例如 PhISCS-BnB 在 single-cell genotype matrix 上比其比較的通用方法快 10–100 倍且保有最優性。但它的資料單位與本專案的 HP-read group Steiner objective 不同，只能作設計參考，不能直接宣稱可替換。[PhISCS-BnB](https://pubmed.ncbi.nlm.nih.gov/32657358/)。

## 9. 建議是否採用

| 項目 | 建議 | 理由 |
|---|---|---|
| 真正 per-unit global deadline＋max sets／bytes | **採用** | 防止長尾拖垮全 run；必須 fail-closed |
| \(h^*\)、candidate completeness、parent topology 三層狀態分開 | **採用** | 避免 timeout／display cap 被誤寫成唯一解 |
| Model A optimum 後做 cheap three-gamete check＋cross certificate | **下一輪 pilot 採用** | 本次 6/8 真實案例直接完成 strict \(h^*\) 證明 |
| Strict perfect 當預設硬模型 | **暫不採用** | 改變生物假設；CN/LOH/recurrence 可使其錯誤；本次也未改善 all-ties |
| Strict perfect 作 CN/LOH-gated sensitivity mode | **可試用** | 能回答「在 infinite-sites 成立時結果是否穩定」 |
| Persistent solver／incremental no-good＋structural cache | **優先工程 pilot** | exact-preserving，對既有 8 h 長尾原因最直接 |
| 現在立刻替換 full pipeline | **NO-GO** | 只有 8 units；M=255 仍未證明，Phase B 0/2 complete |

下一輪最小風險的工程順序是：

1. 加入 per-unit **總** deadline、互斥狀態與 zero-winner abstain。  
   → 驗證：故障注入時 incomplete unit 不產生 read-AF winner。
2. 實作 cross-model strict certificate。  
   → 驗證：small-k exhaustive oracle 0 mismatch；本 8-case panel重播一致。
3. 換成 persistent model＋incremental no-good，不逐 candidate 重建 MILP。  
   → 驗證：所有 complete controls 的 \(h^*\)、完整 candidate digests、complete flags 逐案相同。
4. 再測 33 個既有長尾 units；只在 coverage／wall gate 達標後才評估全量採用。  
   → 驗證：達 cap 全部 abstain，沒有 incomplete ranking；完整率與 wall time另立正式 receipt。

本 pilot 沒有修改 canonical solver、production pipeline 或既有輸出。
