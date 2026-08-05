<!--
建立時間: 2026-07-18 18:00 +08:00
目標: 定義 read-linked Boolean-hypercube solver 的 edge、subcube、一次求解與條件式演化限制改良路線
處理範圍: Task Type A exploratory design；方法契約、small-k oracle、2 個既有 real fixtures 與後續 33-tail pilot 設計；非 production 實作、非全樣本驗證
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py
  - InterSubMod/research/20260718_solver_methyl_edge_probe/20260718_solver與甲基edge小步驗證報告_01.md
  - InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/verify_parent_mapping_and_subcube_reductions.py
  - InterSubMod/research/20260718_solver_methyl_edge_probe/results/parent_mapping_subcube_plan_oracle_receipt.json
  - InterSubMod/research/20260718_solver_methyl_edge_probe/20260718_Hypercube理論縮減與自適應放寬評估_01.md
  - InterSubMod/research/20260718_solver_methyl_edge_probe/results/adaptive_solver_limit_assessment_receipt.json
  - InterSubMod/research/20260718_solver_methyl_edge_probe/results/m2_solver_route_census_receipt.json
  - InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/20260718_PerfectPhylogeny限制與候選上限小型測試_01.md
  - InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/20260718_Hypercube_exact_preserving修補與驗證結果_01.md
-->

# Hypercube edge、subcube 與條件式演化模型改良研究計畫

用 SCQA＋假說驗證三層樓：

> **TL;DR：下一步不應直接放寬全域 \(m\) cap，也不應用「達爾文天擇」或strict infinite-sites取代現行模型。最新route census顯示0.6234%的 \(V\)-hard-tail吃掉99.8376% candidate time；最高優先應改為single-terminal \(V=w!\) 解析route、current complete cache full-run驗證／incomplete checkpoint，再用obligation B&B處理其餘長尾。固定node set後，局部可加edge score仍可用一次 \(O(E_N)\) 掃描；\(m>12\)只對resource gate通過的少數cases作research probe。**（影響：高；數學與單一pilot診斷信心：高；跨樣本效能與生物效度信心：待驗證）

> **PARTIAL／EXPLORATORY**：本文件是方法決策與 bounded pilot 計畫，不是 canonical solver 已修改，也不是 clone/subclone truth 驗證。

本計畫服務 G3（read-level epigenetic 的可解釋角色）、G4（可重現與 fail-closed）及 G5（可由 oracle、receipt 與外部方法比較）。

---

## 0. 先把改善目標說清楚

### 0.1 建議對教授或方法章節使用的主敘述

> 我們把同一 HP×PS 區域中的 read-linked sSNV 狀態表示為 observed-ALT Boolean subcube。完整 read states 是 mandatory states；含 X 的 partial states 則形成「至少選中一個相容 state」的 group constraints。第一階段在 rooted monotone mutation-state graph 上求最少額外 states，並保留候選完整性證明；第二階段以 read-pattern likelihood 排序不同 node sets。改良工作的目標，是在不改變上述 estimand 的前提下，減少重複建模與無效分枝，再把 parent-edge inference 與條件式演化假設分開加入。

### 0.2 不應使用的敘述

- 不說「達爾文天擇保證祖先 VAF 一定較高」。
- 不說「每個新 hidden node 只選最近 parent 就能得到全域最優」。
- 不說「HP 已區分，所以 CNA／LOH 不再影響 molecule-state 解讀」。
- 不說「strict perfect tree 一定比 recurrence-allowed solver 快」。
- 不說「找到前 256 組候選，就已經找到最可能的真實樹」。
- 不把 `minimum_hidden_nodes` 直接翻成真實 hidden clones；較精確名稱是 **minimum-extra-state count**。

### 0.3 本輪研究啟動五問

| 問題 | 判定 |
|---|---|
| Thread D read-level epigenetic 相關？ | 是；但甲基只進 edge soft score／ranking，不作未驗證硬限制 |
| Thread B 撤回範圍？ | 否 |
| KDE-corrected？ | 不適用於結構 solver；若後續使用既有甲基 feature，需沿用各 feature 的既有校正契約 |
| 需要 VCF caller AF？ | edge／cell-fraction constraint 不能直接使用 raw caller VAF；需 CN／purity／multiplicity 校正後的 CCF |
| 長計算／C++／搬移／NO-GO？ | 本計畫只允許 bounded Python pilot；production replacement 與 full run 另設 Gate |

---

## 1. 現行 estimand 與三個不能混寫的證明

令 active sSNV 維度為 \(m\)，候選 state universe：

\[
Q_m=\{0,1\}^m,\qquad |Q_m|=2^m.
\]

`M` 為 root 與達 structural minread 的完整 states；每個 retained partial pattern \(p_j\) 形成相容 subcube \(G_j\)。現行 objective：

\[
h^*=\min_{N\subseteq Q_m}|N\setminus M|
\]

subject to：

\[
M\subseteq N,
\]

\[
N\cap G_j\neq\varnothing\quad\forall j,
\]

\[
\forall v\in N\setminus\{0\},\quad Pred(v)\cap N\neq\varnothing.
\]

方法報告必須分開三層 certificate：

| Certificate | 真正證明的事 | 沒有證明的事 |
|---|---|---|
| Objective certificate | 已證明 minimum-extra \(h^*\) | 未必列完所有同分 node sets |
| Candidate-family certificate | 固定 \(h^*\) 後已列完或壓縮表示全部 optimal \(N\) | 尚未決定 parent edges |
| Edge／biology certificate | 固定 \(N\) 後存在、唯一或最優的 parent assignment；可能另通過演化模型 | 不自動證明每個 state 是 cell clone |

若任一層達 deadline、memory 或 set cap，該層必須標 `incomplete`；不得用候選 prefix 宣稱唯一 winner。

---

## 2. Edge 問題：什麼可以一次算完

### 2.1 Current model 中 edge 為何可後處理

對固定 root-connected node set \(N\)，每個非 root state 的合法父節點集合：

\[
P_N(v)=Pred(v)\cap N.
\]

由於每條 edge 都讓 Hamming weight 增加 1，graph 自動為 DAG。只要 \(P_N(v)\neq\varnothing\)，每個 child 任選一個 parent 即形成 rooted arborescence；不需額外 cycle-elimination constraints。

parent-tree 數：

\[
T(N)=\prod_{v\in N\setminus\{0\}}|P_N(v)|.
\]

因此 current snapshot read likelihood 不依賴 edges 時，正確做法是對每個 distinct \(N\) fit 一次，不展開 \(T(N)\) 棵 parent trees。

### 2.2 Parent mapping 已足以建立整棵樹

rooted tree 建議只保存一種 canonical edge representation：

\[
parent(v)=p,\qquad v\neq0.
\]

因為每個非 root 節點恰有一個 parent，所以一旦全部 `parent(v)` 決定，children可直接反向整理：

\[
Children(p)=\{v:parent(v)=p\}.
\]

因此不應同時獨立搜尋「child的parent」與「parent的children」；兩套獨立決策可能重複計算或產生不一致。

精確規則是：

- root：沒有 parent，可以有零到多個 children；
- 非 root：必須恰有一個 parent；
- leaf：可以沒有 child；
- internal node：可以有一個或多個 children。

使用 current Hamming-1 edge model時：

\[
parent(v)\in Pred(v)\cap N.
\]

每個 child只需從合法 predecessors中選一個 parent；children lists不需另外最佳化。

#### 新增 hidden node 時的兩種方向

1. **Root-outward：parent先存在**

   若 connected parent \(p\) 已存在，新增 child \(h\) 且

   \[
   p\in Pred(h),
   \]

   則 \(h\) 的root連線義務已完成。除非 \(h\) 必須繼續連接某個pending terminal／group，否則不要求它再有 child。

2. **Leaf-inward：child先存在**

   若已有尚未連回root的 child \(c\)，新增其 parent候選 \(h\) 且

   \[
   h\in Pred(c),
   \]

   則 \(c\) 的parent義務完成；但只要 \(h\neq0\)，\(h\) 自己又產生「需要parent」的新義務。不能因為已接到 child 就停止。

若要把 hidden node \(h\) 插在兩個相距兩個 mutations的 endpoints \(p,c\) 中間，必須同時驗證：

\[
p\in Pred(h),\qquad h\in Pred(c).
\]

這正是 obligation B&B適合維護的狀態：每次新增點只更新受影響的 group-hit與orphan-parent obligations，不必重建整棵樹。

Root-outward搜尋可以強制「新增節點當下已有parent」而保持exact，但必須：

- 把一開始尚未連通的mandatory terminals視為pending endpoints；
- 對所有合法frontier children完整分枝；
- 不可只挑排序第一個child。

否則會因插入順序而漏掉合法解。

### 2.3 局部可加 edge score：可用一次掃描求最佳樹

若每條 edge 有不依賴其他 edges 的分數 \(s(p,v)\)，整棵樹分數為：

\[
S(T)=\sum_{(p,v)\in T}s(p,v),
\]

則：

\[
\max_T S(T)
=
\sum_{v\in N\setminus\{0\}}
\max_{p\in P_N(v)}s(p,v).
\]

也就是每個 child 各自選最高分 parent，時間：

\[
O\left(\sum_v|P_N(v)|\right)=O(E_N),
\]

不需要計算 parent Cartesian product。

若需要保留同分 ambiguity，對每個 child 計算最高分 parents 的個數 \(t_v\)，最高分 parent trees 數為：

\[
T_{\mathrm{best}}(N)=\prod_v t_v.
\]

若把edge score視為log weight，還可以一次計算全部parent trees的partition function：

\[
\log Z(N)=
\sum_{v\neq0}
\log\sum_{p\in P_N(v)}e^{s(p,v)},
\]

以及每個child的parent posterior：

\[
P(p\rightarrow v\mid N)
=
\operatorname{softmax}_{p\in P_N(v)}s(p,v).
\]

因此固定 \(N\) 後，可以一次輸出best parent、edge ambiguity與soft posterior；不必materialize所有trees。

本輪 deterministic small-k oracle 對 \(k\le3\) 的 91 個 root-connected node sets，比較「展開全部 parent assignments」與上述解析式，0 mismatch。

### 2.4 何時不能獨立選 parent

下列條件會讓不同 edge choices 互相耦合，不能再逐 child 獨立求最大：

- 每個 mutation 全樹最多 gain 一次；
- parent 的 child cellular prevalence sum 不可超過 parent；
- parent degree、branch 數或深度有 global penalty／cap；
- edge score依賴整條 path、sibling 或跨 edge methylation consistency；
- 同一 molecule／cell 只能被配置到一條 edge；
- CN event 與 sSNV event 必須共同解釋。

此時引入 edge variable：

\[
y_{pv}\in\{0,1\},
\]

\[
y_{pv}\le x_p,\qquad y_{pv}\le x_v,
\]

\[
\sum_{p\in Pred(v)}y_{pv}=x_v\quad(v\neq0).
\]

在 \(m=12\) full active cube，最多有：

\[
m2^{m-1}=12\times2048=24,576
\]

個 unit-flip edge variables。數量仍可 bounded pilot，但不保證 MILP 容易。

若current 4,096個node variables全部保留，node＋edge formulation會上升到28,672個binary variables。當edge沒有獨立score或global coupling時，消去 \(y\) 後仍只得到原本的predecessor condition；它不會強化整數vertex-set family，卻會增加model大小。因此不應為「只是輸出tree」而預設加入edge variables。

### 2.5 一次 joint node＋edge 最佳化的精確邊界

建議分兩次 lexicographic solve：

1. 先求並固定 \(h=h^*\)。
2. 再在所有 \(h^*\) node sets 與 edges 中最大化可信 edge score。

這可以用一次第二階段 MILP 找到 **top-1 node＋edge solution**，但不能宣稱：

- 已列完所有同分 node sets；
- edge winner 有獨立生物證據；
- likelihood winner 與 edge-score winner相同；
- top-1 之外沒有跨 topology ambiguity。

若目標是報告完整不確定性，仍需 solution family 或其壓縮表示。

### 2.6 無證據 unary chain 可考慮壓成「未定序 mutation-set edge」

Strict infinite-sites仍不能消除單一路徑的mutation-order排列。例如只有 terminal `1111` 時，current Hamming-1展開仍有：

\[
4!=24
\]

條合法 chains。

若一段 path 的中間 states：

- 沒有 full-read直接支持；
- 沒有 partial group需要命中；
- 沒有低支持read在likelihood層提供可區分訊號；
- 在所有候選中都只是 in-degree=1、out-degree=1 connector；

可另建新模型，把：

```text
0000 → hidden → hidden → hidden → 1111
```

壓成：

```text
0000 ──{mutation 1,2,3,4；order unresolved}──▶ 1111
```

輸出應是「1個edge-equivalence class、24個未定順序」，不是選其中一條chain當真實歷史。

此做法可能最直接減少無資料支持的 \(k!\) edge爆炸，但它會改變：

- Hamming-1 edge contract；
- minimum-extra-state objective；
- 可供mixture likelihood配置權重的states。

因此只能作獨立 `MULTI_MUTATION_EDGE_EQUIVALENCE` pilot，不能稱為current solver的純加速；只要任一中間state有read／likelihood evidence就不得collapse。

---

## 3. Subcube 還能如何 exact-preserving 改善

### 3.1 已有且應保留

1. observed-ALT active dimensions：

   \[
   2^k\rightarrow2^m.
   \]

2. 相同 symbolic group 合併。
3. mandatory-hit、singleton forcing。
4. group subset dominance：若 \(G_a\subseteq G_b\)，structural stage 可移除 \(G_b\)。
5. sparse group-hit row，而非跨 read completion Cartesian product。
6. complete-only structural cache。

### 3.2 新增 downward-closure universe reduction

令：

\[
D=
\downarrow
\left(
F\cup\bigcup_jG_j
\right),
\]

其中 \(\downarrow v\) 是 \(v\) 的所有 Boolean ancestors。

在 current rooted monotone、minimum-extra objective 下，minimum solution 中任何非 mandatory、非 group representative 的 vertex 都必須服務於某個 terminal／representative 的 root path；因此 optimum 不需要使用 \(D\) 外的 vertex。

可先把候選 universe 從 \(Q_m\) 限制成 \(D\)：

\[
Q_m\rightarrow D.
\]

本輪獨立 exhaustive check：

- \(k\le3\)；
- 6,243 個 mandatory＋最多兩個 partial groups cases；
- 11,432 個 minimum node sets；
- 所有 minimum sets 均位於 \(D\)，0 mismatch。

但兩個既有 M31 real fixtures 的 \(D\) 都仍為完整 32-state active cube，縮減 0%。因此這是安全 reduction probe，不可預先宣稱 real data 會顯著加速。

### 3.3 Bitset＋antichain

在 \(m\le12\) 時，一個完整 group domain 僅需 4,096 bits，即 512 bytes。建議：

- 每個 \(G_j\)、selected、excluded、predecessor domain 都用 fixed-size bitset；
- duplicate 用 hash；
- subset dominance 用 bit operation；
- obligation hit、MRV domain size、singleton propagation 使用 bitwise AND/popcount；
- scoring counts仍保留，不因 structural dominance遺失。

此改動只改表示與運算，不改 estimand，適合最先做 regression oracle。

### 3.4 Obligation-driven exact B&B

search state：

\[
(I,E)
\]

分別是 selected 與 excluded vertices。未滿足義務：

\[
Q(I)=
\{G_j:I\cap G_j=\varnothing\}
\cup
\{Pred(v):v\in I\setminus\{0\},\
I\cap Pred(v)=\varnothing\}.
\]

必要元件：

- singleton propagation；
- MRV：優先分枝最小 domain；
- root reachability prune；
- safe lower bounds 取 `max`，不可將重疊 bounds直接相加；
- `(selected, excluded)` memo；
- 對 chosen obligation 的全部 candidates完整分枝；
- 只有 `LB > best` 才剪枝，`LB == best` 必須保留以列同分解。

這個 traversal 可在一次 search 中同時找 \(h^*\) 與全部 minimum sets，避免 current MILP 的約 \(V+1\) 次 solves。現有兩個 M31 probe 的完整 candidate digest 已與 MILP 一致，但尚未重播 33 個正式長尾 units。

### 3.4B 若保留MILP，可加入lazy rank-layer cuts

current predecessor row對整數解已充分，但LP relaxation可能把fractional support分散在多個levels。對任一selected \(v\) 與中間rank \(r\)，任何root-to-\(v\) path都必須通過一個rank-\(r\) ancestor，因此可加入：

\[
x_v
\le
\sum_{\substack{u\subset v\\|u|=r}}x_u,
\qquad 1\le r<|v|.
\]

這些cuts不改整數可行集合，只可能改善lower bound／branching。但若全部預建，總NNZ可接近 \(3^m\)；建議只在root LP或incumbent顯示violation時lazy加入，並與node-only baseline做A/B。

### 3.5 ZDD／BDD：壓縮「答案家族」，不是保證變快

若全部 minimum sets 數 \(V\) 很大，但具有大量共同 selected/excluded suffix，可用 ZDD 表示 family：

- 不先 materialize 每個 \(N\)；
- 可計數、抽樣、套額外 filter；
- 可查各 vertex inclusion probability／topology class；
- 需要時才逐解輸出。

限制：

- decision diagram 大小高度依賴 variable order；
- 某些 family 本身不可壓縮，ZDD 仍會爆；
- current group-Steiner-like problem不是一般 Steiner tree，需自建 frontier state；
- 若最後仍要逐一 fit 非線性 mixture likelihood，輸出量瓶頸仍存在。

因此 ZDD 應在 B&B dual validation 後作第二線 pilot。

### 3.6 Group-to-terminal subset DP：q 小時可一次求 \(h^*\)

令 \(q\) 為 reduction 後的 mandatory terminals與partial groups數。可將每個 group \(G_j\) 轉成一個人工 terminal \(\tau_j\)，並從每個相容 state加入零成本edge：

\[
v\rightarrow\tau_j,\qquad v\in G_j.
\]

再把vertex cost設為：

- root／mandatory full state：0；
- 其他 mutation states：1；
- artificial group terminals：0。

原問題的 objective即可轉成rooted node-weighted directed Steiner objective。以 terminal subset作DP state，可在一次dynamic program中求一個 \(h^*\)：

\[
DP[S,v]=\text{覆蓋terminal subset }S\text{ 並在 }v\text{ 合併的最小成本}.
\]

典型時間與空間量級：

\[
O(3^q2^m+2^q m2^{m-1}),
\]

\[
O(2^q2^m),
\]

另加group-terminal edges \(\sum_j|G_j|\) 的處理。

這對「\(m\)中等、reduced groups \(q\)很小」可能比對 \(2^m\) 個vertices作通用MILP更好。既有兩個M31 fixtures經current exact reduction後：

```text
H2009_M31：12 input groups → 4 reduced groups
COLO829_M31：4 input groups → 2 reduced groups
```

因此值得small pilot。

限制：

- DP先解一個optimum／objective certificate，不自動輸出全部optimal node sets；
- 若要exact count／all sets，需保留backpointer DAG或接ZDD，仍受output family大小限制；
- \(q\)大時 \(3^q\) 會反向爆炸；
- node cost在subset merge時不可重複計算，必須用small-k oracle驗證。

### 3.7 Incremental SAT／PB 與 column generation

目前 constraints 多為 Boolean clauses＋cardinality objective，適合評估 incremental SAT／pseudo-Boolean backend：

- group-hit：一個 positive clause；
- predecessor：\(\neg x_v\lor\bigvee_{p\in Pred(v)}x_p\)；
- mandatory：unit clauses；
- minimum extras：cardinality／PB objective；
- all-optimal enumeration：固定 \(h^*\) 後 incremental blocking。

這可能比通用 MILP 更貼近模型，但仍有 all-solutions output lower bound。

若 \(m>12\) 且只需要少數 states，可研究 column generation／branch-and-price，只生成能改善 reduced cost 的 vertices；但 predecessor closure與group coverage使 pricing problem本身不簡單，列為後期研究，不作第一實作。

---

## 4. Darwinian selection 與 infinite-sites 能限制什麼

### 4.1 達爾文式腫瘤演化不是可直接寫入 solver 的 hard constraint

Darwinian／clonal evolution提供的生物背景是：

- 腫瘤內存在可遺傳變異；
- 不同細胞族群在特定環境下有差異性增殖／存活；
- clone frequency 隨時間、治療、微環境與隨機漂變改變。

它沒有保證：

- 後出現的 mutation一定使 fitness上升；
- ancestor當下頻率一定高於所有 descendants；
- 同一 mutation不會平行出現；
- mutation不會因 deletion／LOH消失；
- 單一時間點 bulk／read fractions能唯一決定時間方向。

因此「達爾文天擇」只能作研究動機或需要 longitudinal data 的 fitness prior，不能直接作目前 cross-sectional tree 的剪枝定理。

### 4.2 Infinite-sites 是額外模型假設，不是天擇推論

Strict infinite-sites／perfect phylogeny 要求每個 mutation：

- 全樹只 gain 一次；
- 不 loss；
- 不 recurrence／back mutation。

若使用 edge variables，可寫成：

\[
\sum_{(p,v):p_j=0,\ v_j=1}y_{pv}\le1
\quad\forall j.
\]

對有觀測 ALT 的 locus，可視 contract 改成恰好 1。

更便宜的current-cube實作可先使用node-only rooted three-gamete condition，或只對實際violation加入lazy witness cut：

\[
x_{10}+x_{01}+x_{11}\le2.
\]

若root-connected \(N\) 通過three-gamete，任一非root node不可能同時有兩個selected unit predecessors，因此：

\[
T(N)=1.
\]

也就是strict-compatible node set的parent edge可直接解析得到，不需要為此全面加入edge variables。

它能排除 `00、10、01、11` diamond：不管 `11` 接 `10` 或 `01`，都會讓另一個 mutation gain 第二次。本輪 small-k oracle 對該 diamond 的兩個 parent assignments均判 strict=false。

但 strict constraint 不保證較快。既有 bounded pilot：

- 8 個 real capped units中，cross-certificate 可證明 6 個 strict \(h^*\)；
- 沒有證明 strict 改變這 6 個 units的 minimum objective；
- 兩個全 strict candidate enumeration 在 15 秒內皆 incomplete。

### 4.3 建議使用四層演化模型

| 模型層 | 假設 | 建議角色 |
|---|---|---|
| M0 recurrence-allowed monotone | 每條 path只 0→1，但全樹可多次 gain | current baseline；所有區域保留 |
| M1 strict infinite-sites | gain一次、不可 loss | 只在 CN/LOH eligibility 通過時作 sensitivity／certificate |
| M2 loss-supported Dollo | gain最多一次；CN deletion／LOH可支持 loss | CNA／LOH 區域的下一個生物模型 |
| M3 copy-aware／finite-sites | 可多次 gain/loss；state含 allele-specific CN | 長期模型；最接近複雜癌症，但計算與可識別性最困難 |

### 4.4 M1 eligibility gate

只有同時滿足下列條件，才允許把 strict infinite-sites 提升為較強的 hard sensitivity model：

1. allele-specific CN顯示 copy-neutral、兩個 homolog均保留；`total CN=2` 本身不足，因 2+0 copy-neutral LOH 仍不合格；
2. 無 clonal／subclonal LOH、deletion、amplification、WGD boundary；
3. 同一 HP×exact PS內 phasing、mapping、BQ、strand與read independence QC通過；
4. 無 allele-specific duplicated copies 的證據；
5. sSNV本身為高信心 somatic，不是 mapping artifact／germline leakage；
6. 報告名稱仍是 `strict-perfect-compatible molecule-state topology`，不是 clone truth。

HP tag只能區分 germline-derived haplotypes，不能區分同一 HP在 CNA後衍生的多份 copies，因此 HP存在不等於 M1 eligibility自動通過。

### 4.5 VAF／CCF constraint 的正確層級

若 parent mutation存在於所有 descendant cells，理論上的 **cellular prevalence** 滿足：

\[
\phi_{\mathrm{ancestor}}\ge\phi_{\mathrm{descendant}}.
\]

對互斥 child lineages，還有 sum rule：

\[
\phi_p\ge\sum_{c\in Child(p)}\phi_c.
\]

但這裡的 \(\phi\) 必須是經 tumor purity、allele-specific CN、mutation multiplicity校正的 CCF／subclonal frequency，不是 raw VAF，也不是同一批 reads的 molecule fraction。

因此：

- raw read-AF／VAF：ranking／diagnostic；
- CN-corrected CCF＋uncertainty interval：可作 conditional edge feasibility；
- uncertainty intervals overlap：TIE／ABSTAIN；
- CNA／LOH不明：不得作 hard ancestor constraint。

---

## 5. 是否能把 structure＋likelihood 一次求完

### 5.1 可以，但會改變輸出契約

現行流程先列完 minimum-extra \(N\)，再對每個 \(N\) fit mixture likelihood。若只想要 top-1，可考慮：

\[
\text{lexicographic minimize }h(N),
\quad
\text{then maximize }\ell(N,\pi)+\lambda S_{\mathrm{edge}}.
\]

其中：

\[
\pi_v\ge0,\qquad \sum_v\pi_v=1,\qquad \pi_v\le x_v.
\]

read log-likelihood含：

\[
\sum_r n_r\log\left(\sum_v\pi_vq_{rv}\right).
\]

它不是 current SciPy MILP 的線性 objective；需要 mixed-integer convex optimization、outer approximation或具有可證 upper bound 的 custom B&B。

### 5.2 建議定位

| 目標 | 建議 |
|---|---|
| 找一個 \(h^*\) | current MILP或exact B&B均可 |
| 列完所有 minimum \(N\) | exact B&B／incremental PB；必要時ZDD壓縮 |
| 固定 \(N\) 找最佳local additive edges | 解析式一次 \(O(E_N)\) |
| 在全部 \(h^*\) solutions 中找 top-1 edge-linear score | 固定 \(h^*\) 後一個 edge-variable MILP |
| 在全部 \(h^*\) solutions 中找 top-1 mixture likelihood | 後期 MICP/custom B&B pilot |
| 保留完整 likelihood ambiguity | 仍需枚舉或能對 likelihood作動態規劃的壓縮family；目前尚無 |

若研究目標是向教授呈現「還有多少未解 ambiguity」，不可只輸出 joint top-1；否則會把 solver選擇誤寫成資料已確認。

---

## 6. Pre-decision 判定

| 改良項目 | 判定 | 原因 |
|---|---|---|
| Obligation exact B&B dual backend | **PROBE，最高優先** | 同 estimand；兩個 M31完整digest一致且顯著少訪問 states |
| Bitset＋downward closure | **PROBE，高優先** | exact-preserving；small-k oracle通過；real縮減幅度未知 |
| q-parameterized group-terminal subset DP | **PROBE，高優先** | 可一次求objective；兩個M31 reduced q僅4／2，但尚無實作oracle |
| 固定 \(N\) 的 additive edge解析式 | **GO for isolated implementation** | 可由數學分解與 exhaustive oracle直接驗證；不改current structure |
| Edge-variable top-1 MILP | **PROBE，中優先** | 可避免parent tree展開，但會新增edge estimand |
| ZDD壓縮all-optimal family | **PROBE，中後期** | 可能解決儲存／重複suffix；不保證diagram不爆 |
| Incremental SAT／PB backend | **PROBE，中後期** | constraints吻合；需與MILP全family oracle比較 |
| Joint structure＋mixture likelihood | **DEFER** | 非線性、改輸出契約，先解結構長尾 |
| Strict infinite-sites作全域預設 | **NO-GO** | CN/LOH/recurrence可違反；既有pilot未解all-ties |
| CN/LOH-gated strict sensitivity | **PROBE** | 可提供穩健性與cross-model certificate |
| Darwinian selection作hard prune | **NO-GO** | 非可觀察的拓撲定理；單時間點資料不可識別fitness |

---

## 7. 分階段實驗計畫

### Phase 0：凍結 contract 與 oracles

**研究問題**：新 backend／reduction 是否完全保留 current M0 estimand？

1. 固定 solver、ranker、fixture SHA。  
   → 驗證：receipt含source bytes、Python／SciPy／HiGHS版本。
2. 建立 \(k\le4\) exhaustive oracle。  
   → 驗證：objective、全部 minimum set digest、parent count逐案一致。
3. 固定五種狀態：objective-only、family-complete、incomplete-cap、incomplete-deadline、feasible-unproven。  
   → 驗證：故障注入不產生winner。

**失敗即停止**：任何 complete case digest mismatch或 incomplete case仍進ranking。

### Phase 1：Exact B&B dual validation

**假設 H1**：obligation B&B可在不改答案下減少33個正式長尾的candidate-generation wall。

1. B&B只作research backend，不替換production。  
   → 驗證：同一 unit 同時跑 current MILP與B&B。
2. 先跑所有既有 MILP-complete small units。  
   → 驗證：\(h^*\)、ordered／unordered family digest、complete flag 100%一致。
3. 再跑33個既有長尾 units。  
   → 驗證：記錄visited、propagation、prune、wall、peak RSS、candidate count。

**成功門檻**：

- complete controls：0 digest mismatch；
- incomplete不產生ranking winner；
- 33-tail candidate-generation總wall至少下降 5×，且p95不劣於current；
- peak RSS不超過current 1.5×；
- 若答案本身 \(V>cap\)，正確標family incomplete。

### Phase 2：Subcube exact reductions

**假設 H2**：bitset＋downward closure能降低部分units的active vertices／row nonzeros，且不改candidate family。

1. small-k exhaustive驗證downward closure。  
   → 驗證：擴至 \(k\le4\) seeded＋edge cases，0 mismatch。
2. 對33 tails與隨機complete controls計算 reduction-only census。  
   → 驗證：輸出 raw \(2^m\)、closure \(|D|\)、groups before/after、NNZ before/after。
3. 只有有縮減的units才跑 on/off solver。  
   → 驗證：digest、objective、complete flags一致。

**判定**：

- correctness是必要 Gate；
- real median reduction即使為0也算有效負結果；
- 不以兩個M31 0% reduction否定所有資料，也不預設一定加速。

### Phase 2B：q-parameterized subset-DP objective backend

**假設 H2b**：當reduced terminal/group數 \(q\)小時，group-to-terminal subset DP可一次證明 \(h^*\)，並與current MILP／B&B一致。

1. 實作node-weighted artificial-terminal轉換。  
   → 驗證：每個group terminal恰由至少一個selected compatible state到達。
2. \(k\le4\) exhaustive比較。  
   → 驗證：objective與MILP／brute-force 0 mismatch。
3. 路由只先測 \(q\le8\) units。  
   → 驗證：記錄DP states、wall、RSS、objective certificate；不冒充family complete。
4. 若要重建全部sets，另測backpointer DAG／ZDD。  
   → 驗證：可枚舉small cases的family digest一致。

**停止條件**：DP objective mismatch、\(q\le8\)仍常比B&B慢／耗記憶體，或為取得family completeness必須展開同樣大的V。

### Phase 3：Edge layer

**假設 H3a**：固定 \(N\) 的 additive edge scorer可用解析式取代parent Cartesian enumeration。

1. \(k\le5\) exhaustive parent oracle。  
   → 驗證：best score、best-tree count、unique/tie狀態0 mismatch。
2. current read likelihood保持vertex-set only。  
   → 驗證：開關edge scorer不改primary likelihood。
3. 甲基／CCF edge evidence各自校正。  
   → 驗證：held-out或cross-fit；CI跨0即TIE；eligibility不足即ABSTAIN。

**假設 H3b**：需要跨edge constraint時，edge-variable MILP可在固定 \(h^*\) 下找top-1且與small-k exhaustive一致。

1. 加 \(y_{pv}\) 與parent equality。  
   → 驗證：每個selected child恰一parent、未選edge為0。
2. 加 strict gain-once或CCF sum rule作獨立mode。  
   → 驗證：diamond strict infeasible；正常chain feasible。
3. 不列完family時輸出 `TOP1_ONLY`。  
   → 驗證：UI／TSV不可顯示candidate-family complete。

### Phase 4：演化模型 sensitivity

**假設 H4**：M1只在CN/LOH合格區域能安全縮小ambiguity；M2在loss-supported cases比M1少產生假性infeasible。

1. 建立CN/LOH eligibility表。  
   → 驗證：每unit有CN state、minor CN、LOH、subclonal-CNA、source與uncertainty。
2. M0→M1 cross-certificate。  
   → 驗證：若M0 optimum本身strict，直接證明 \(h^*_{M0}=h^*_{M1}\)。
3. 建立synthetic CN/LOH反例：2+0 cnLOH、deletion loss、allele-specific amplification。  
   → 驗證：M1正確fail／abstain，M2能以明示loss event解釋。
4. 比較M0/M1/M2的objective、family size、topology、abstain。  
   → 驗證：分層報告，不把模型間差異當truth。

### Phase 5：Joint top-1 research prototype

**前提**：Phase 1–4均通過才啟動。

**假設 H5**：在只求top-1時，joint optimizer可避免顯式列完所有 \(h^*\) sets，且在small complete units與enumerate-then-rank winner一致。

1. 固定 \(h^*\)，先只整合edge-linear score。  
   → 驗證：small exhaustive winner／tie一致。
2. 再評估mixture likelihood upper bound。  
   → 驗證：global gap certificate，不接受只有local optimizer success。
3. 對高ambiguity units同時保留family count lower bound與top-1-only標記。  
   → 驗證：不把效率模式混入完整拓撲比例。

---

## 8. 固定輸出與評估指標

每個 unit 至少輸出：

| 類別 | 欄位 |
|---|---|
| Scope | sample、region、HP family、PS、raw k、effective m |
| Subcube | raw vertices、closure vertices、groups before/after、group NNZ |
| Objective | incumbent、lower bound、gap、\(h^*\) certificate |
| Search | backend、wall、CPU、peak RSS、visited、propagated、pruned、solve calls |
| Family | optimal sets found、complete、family digest或ZDD size |
| Edges | parent count、best edge score、best-edge count、unique/tie/abstain |
| Biology | CN、minor CN、LOH、CNA eligibility、M0/M1/M2 mode |
| Ranking | primary log-likelihood、delta、relative weight、convergence／global gap |
| Failure | deadline、cap、memory、numerical、model-infeasible、QC-abstain |

Aggregate需同時報：

- completed／objective-only／incomplete／abstain數與比例；
- median／p95／p99 wall；
- peak RSS；
- candidate-generation與likelihood分開時間；
- strict eligibility分母；
- M0→M1 candidate縮減與infeasible率；
- edge unique／tie／abstain；
- 所有oracle mismatch數，必須為0。

---

## 9. 主要風險與停損條件

| 風險 | 停損／處理 |
|---|---|
| B&B對weak constraints退化 | 保留MILP router與per-unit deadline；依instance特徵選backend |
| 答案family本身巨大 | ZDD／count-only／top-1-only明示模式；不偽稱complete |
| Downward closure real縮減為0 | 記錄負結果，不繼續複雜化 |
| Edge score重用同一批reads | 不把same-read VAF當獨立term；甲基用cross-fit／held-out |
| CN/LOH錯誤分類 | strict只作sensitivity；eligibility帶uncertainty與abstain |
| HP誤解為physical copy | 報告明示HP-derived copies在CNA後仍不可區分 |
| Darwinian敘事過度 | 不寫fitness／selection hard constraint；需longitudinal或functional資料才估fitness |
| Joint optimizer只找到local optimum | 沒有global bound即不升級 |
| Source／receipt drift | 每輪綁SHA；historical wall不可當current production benchmark |

---

## 10. 建議的實作順序

```text
P0 contract／oracle
  ↓
P0a single-terminal V=w! analytic／compact-family route
  ↓
P0b 驗證current complete cache＋新增incomplete checkpoint
  ↓
P1 exact B&B dual backend＋其餘33-tail重播
  ↓
P2 bitset＋downward-closure census
  ↓
P2B q-small subset-DP objective pilot
  ↓
P3a fixed-N additive edge解析式＋parent mapping
  ↓
P3b edge-variable top-1 mode
  ↓
P3c unary-chain equivalence pilot
  ↓
P4 CN/LOH-gated M1／loss-supported M2 sensitivity
  ↓
P5 joint structure＋likelihood research prototype
```

最近一步只建議批准：

1. P0a single-terminal analytic／compact-family route；
2. P0b current complete cache full-run驗證＋incomplete checkpoint；
3. P1 exact B&B dual validation；
4. P2B q-small subset-DP只對adaptive gate通過cases作objective pilot；
5. P3a fixed-\(N\) edge解析式oracle。

不建議現在批准：

1. strict infinite-sites取代M0；
2. Darwinian hard pruning；
3. joint nonlinear optimizer直接進production；
4. 未完成candidate family仍輸出VAF／甲基winner。

---

## 11. 本輪已執行的兩個獨立數學檢查

輸入：

- `InterSubMod/research/20260718_solver_methyl_edge_probe/tests/fixtures/real_units.json`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`
- 所有 \(k\le3\) synthetic Boolean node sets／R-A-X group cases。

執行命令：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/verify_parent_mapping_and_subcube_reductions.py \
  --output \
/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_solver_methyl_edge_probe/results/parent_mapping_subcube_plan_oracle_receipt.json
```

輸出：

- [可重播oracle腳本](/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/verify_parent_mapping_and_subcube_reductions.py)
- [JSON receipt](/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_solver_methyl_edge_probe/results/parent_mapping_subcube_plan_oracle_receipt.json)

實際stdout：

```text
PASS parent_sets=91 closure_cases=6243 closure_optimal_sets=11432 real_cases=2
```

### 11.1 Edge-local additive score

輸入：

- 所有 \(k\le3\) root-connected Boolean node sets；
- deterministic local edge scores。

執行方法：

- 明列所有 parent Cartesian assignments；
- 與 \(\sum_v\max_p s(p,v)\) 比較。

實際輸出：

```text
PASS edge_local_additive analytic_equals_exhaustive
root_connected_sets=91 k_max=3
```

Diamond：

```text
00,10,01,11
11→10：mutation 2 gain兩次，strict=false
11→01：mutation 1 gain兩次，strict=false
```

### 11.2 Downward-closure reduction

輸入：

- \(k\le3\)；
- mandatory states；
- 0–2個R/A/X partial groups。

實際輸出：

```text
PASS downward_closure_contains_all_minimum_sets
cases=6243 optimal_sets=11432 k_max=3
```

既有兩個M31 fixture：

```text
H2009_M31   active cube=32  downward closure=32  reduction=0.0%
COLO829_M31 active cube=32  downward closure=32  reduction=0.0%
```

這證明目前只適合把downward closure列為safe probe，不能宣稱已在real fixtures縮小。

---

## 12. 外部方法學對照

- Nowell 的 clonal evolution提出腫瘤細胞族群中的遺傳變異與選擇，但並未給出「mutation只發生一次」或「單一時間點頻率唯一決定祖先」的演算法約束：  
  [Nowell 1976, Science](https://pubmed.ncbi.nlm.nih.gov/959840/)
- LICHeE在perfect-phylogeny框架使用mutation presence、VAF ordering與sum rule；原文同時指出CNV會破壞VAF ordering，並建議使用已校正CN／LOH／purity的cellular prevalence：  
  [Popic et al. 2015, Genome Biology](https://link.springer.com/article/10.1186/s13059-015-0647-8)
- PhyloWGS把CNV影響納入VAF／population-frequency模型，顯示copy-number correction與subclonal reconstruction不可分離：  
  [Deshwar et al. 2015, Genome Biology](https://link.springer.com/article/10.1186/s13059-015-0602-8)
- Pairtree使用pairwise ancestry與tree-constrained subclonal frequencies，並明示subclone frequency是其所有descendant populations之和；這需要cellular frequency，不是raw read-AF：  
  [Wintersinger et al., Pairtree](https://pmc.ncbi.nlm.nih.gov/articles/PMC9780082/)
- PhISCS與SCARLET說明LOH、deletion、convergent evolution及copy-number error會違反infinite-sites；因此以subperfect／loss-supported model取代全域strict假設較合理：  
  [PhISCS](https://pmc.ncbi.nlm.nih.gov/articles/PMC6836735/)；[SCARLET](https://pmc.ncbi.nlm.nih.gov/articles/PMC7451135/)
- 有實證研究直接比較infinite-sites與finite-sites，觀察到parallel hits與mutation loss；strict ISA應被視為可檢驗模型而非永真公理：  
  [Kuipers et al. 2017, Genome Research](https://pmc.ncbi.nlm.nih.gov/articles/PMC5668945/)
- ZDD／BDD能壓縮並延後展開Steiner solution families，但diagram大小與變數順序有關，不能保證消除所有組合爆炸：  
  [Sasaki 2021, cost-constrained Steiner-tree BDD](https://arxiv.org/abs/2104.06696)

---

## 13. 本輪閱讀建議

1. [Partial-read方法獨立稽核](/big7_disk/liaoyoyo2001/InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md)
2. [Exact B&B與甲基edge bounded probe](/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_solver_methyl_edge_probe/20260718_solver與甲基edge小步驗證報告_01.md)
3. [Perfect-phylogeny限制pilot](/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_perfect_phylogeny_constraint_pilot/20260718_PerfectPhylogeny限制與候選上限小型測試_01.md)
4. [Prepared-base與complete-only cache修補](/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/20260718_Hypercube_exact_preserving修補與驗證結果_01.md)
5. [互動式Boolean-hypercube教學](/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_solver_methyl_edge_probe/docs/20260718_布爾超立方體k與組合數量公式教學_01.standalone.html)
6. [理論縮減與自適應放寬評估](/big7_disk/liaoyoyo2001/InterSubMod/research/20260718_solver_methyl_edge_probe/20260718_Hypercube理論縮減與自適應放寬評估_01.md)

---

## 14. 量化後的放寬決策

完整計算與分母稽核見第13節第6份文件；核心數字如下：

- active-bit real fixture：raw \(k=7\to m=5\)，vertices `128→32`（-75.00%）、edges `448→80`（-82.14%）。
- exact group reduction：H2009 `12→4`（-66.67%）；COLO829 `4→2`（-50.00%）。
- \(m=5,h^*=5\) naive subset layers `206,368`；B&B只visited `814/531`（-99.606%／-99.743%），完整 optimal-set digests仍一致。
- proposed small-\(q\) DP在50M proxy operations、512 MiB、16 B/cell假設下：\(q\le6\)可到 \(m=15\)，\(q\le4\)可到 \(m=17\)，\(q\le2\)理論可到 \(m=19\)。這些是pilot gate，不是runtime保證。
- M2 frozen-v2 diagnostic中，只有31/7,058（0.439%）因 \(m>12\) 未啟動solver，resource gate只接受其中2個；既有33個tails全部在 \(m=6\)–11。
- 11個已完成 \(V\ge100\) units加33 tails僅44/7,058（0.6234%），卻消耗99.8376% candidate-generation time。
- 4,380個completed single-terminal cases均滿足 \(V=w!\)，0 mismatch；另3個tails的 \(V=720/5,040/720\)在現有cap下結構上不可能complete。
- 4,852個completed evaluations只有165個exact structural keys，solver-call重用機會96.599%；tails 33→19 keys，但incomplete只可checkpoint，不可當complete cache。
- current-v5的44,672/72,994（61.20%）primary units為partial-only；H2009＋H1437占7,975個incomplete的68.89%，應列入stress panel。

因此限制政策定為：

1. raw \(k\) 不單獨作拒絕；先做exact active-bit compression。
2. production effective \(m\le12\)暫不變。
3. 先實作single-terminal exact count／compact family，再驗證current complete cache並新增incomplete checkpoint。
4. 其餘hard-tail優先走obligation B&B；small-\(q\) DP先作objective oracle。
5. \(m>12\)只對預估operations、memory與 \(V\) gate通過的cases放行；目前先測2／31，不設固定新cap。
6. 不直接提高`max_sets=256`；改成objective、uniqueness、count/compressed-family與full-enumeration模式。
7. Historical layered-v2 VAF的25,978/28,183=92.18%只作ranking endpoint設計；不可與current-v5的42,240 complete denominator混用。

---

**目前決策**：`PROBE_ADAPTIVE_RELAXATION`，但實作優先序改為P0a factorial route → P0b cache/checkpoint → P1 obligation B&B → adaptive \(m>12\) probe。HCC1395 chr6、H2009／H1437 stress panel及small-\(m\) oracle通過前，不修改canonical／production solver。Strict infinite-sites維持CN/LOH-gated sensitivity；Darwinian selection不作hard constraint。
