<!--
建立時間: 2026-07-12 +0800
目標: 在實作前質疑、形式化並驗證「以 ONT read 支持度加權 Boolean-hypercube rooted group-Steiner tree」的可行性
處理範圍: Task Type A exploratory method probe；現行 solver/read-AF 契約、k=3 反例、primary literature、bounded pilot 設計；未修改 solver、未執行全樣本長計算
服務目標: G3 read-level evidence；G4 reproducibility；G5 externally verifiable contribution
證據狀態: PARTIAL / METHOD PROBE；不是 biological validation，不是 confirmed clone tree
cycle_id: none_precycle
status: verdict_PROBE_reframed_model; verdict_NO_GO_raw_edge_count_and_uniform_inflow
audit_version: 0.2
關聯檔案:
  - InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_linkage_genomewide.py
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/pre-decision-audit.md
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_top_tie_census.py
-->

# Pre-Decision Audit：read-pattern likelihood 與 Boolean-hypercube 候選樹

> **TL;DR — 原始寫法 `read count → Hamming-1 edge weight → 最佳路徑` 不成立；改寫為「保留全部最小候選樹，再用 partial/full read-pattern whole-tree likelihood 排序並保留不確定性」值得做 bounded pilot（影響：高；信心：中）。**

> **決策**：原始 raw edge-count 版本 **NO-GO as stated**；重新形式化的 read-pattern likelihood 版本 **PROBE（40/100）**。本輪只完成前驗與計畫，不改 solver、不把任何候選稱為真實演化樹。

## 0. 任務分類、範圍與研究前置 5 問

- **Task Type**：A — exploratory pilot。
- **Scope**：先做 `k=3` synthetic／hand-solvable cases，再做 HCC1395 bounded subset；不是全基因組、全樣本 validation。
- **Partial flag**：本文件不可作論文效能證據；後續若進 Type B，必須補 chr1–22 × 全 dataset rows 的 verification chain。
- **G3**：把 read-level joint genetic evidence 放進候選樹評分。
- **G4**：用 HCC1395／HCC1395_DORADO 檢查技術重現，之後才可擴其他 biological samples。
- **G5**：以可否證 simulation、held-out prediction 與 calibration 形成外部可驗證方法。

研究任務 5 問：

1. **Thread D？** 間接相關；是 read-level genetic evidence，不是 methylation 主題。
2. **Thread B 撤回範圍？** 否；不重啟已撤回的 methylation／variant-filter 假說。
3. **KDE-corrected？** 本問題不使用 KDE，N/A。
4. **需要 caller AF？** read-pattern likelihood 應以原始 read observations 為主；caller AF 只可作外部校準欄位。若 read-AF/VAF 來自同一 BAM，不可再乘一次造成 double counting。
5. **長計算／C++／檔案搬移／NO-GO gate？** 本輪無；若要改 solver objective，屬新 spec 與高影響方法變更，須先通過本 audit 的 pilot gate。

Cynefin：**Complex**。相同的 raw count weighting 尚未有可預測、可重現的 truth 結果，因此必須 probe-first。

## 1. 先把三個不同概念拆開

### 1.1 現行 Hamming=1 是「合法性」，不是數值 cost

現行狀態為 `x ∈ {0,1}^k`，root 為 `0^k`。合法有向邊只允許一次 `0→1`：

\[
x\to y \iff y=x+e_j,\qquad d_H(x,y)=1.
\]

因此每條合法邊的 Hamming 距離都是 1，但程式不是在最小化 `Σ edge weight`。實際 objective 是：

1. 滿足所有 full／partial observations 的 coverage；
2. 找到最少 hidden／extra nodes 的所有 feasible node sets；
3. 對每個 node 枚舉全部合法 unit predecessor；
4. 保留所有同層最小候選樹。

程式證據：

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py:55`：unit predecessor 定義。
- 同檔 `:66`：full 與 partial-subcube coverage。
- 同檔 `:113`：`full_counts` 被讀入，但後續 objective 完全未使用 count。
- 同檔 `:142`：依 hidden-node 數逐層搜尋。
- 同檔 `:211`：對固定 node set 枚舉所有 parent choices。

### 1.2 read 觀察「現存 genotype state」，沒有觀察演化 transition

ONT read 是某個時間點抽到的一條 DNA molecule。它可提供：

- 哪些 mutation 在同一 molecule 上共現；
- 哪些 mutation 在該 molecule 上沒有被觀察到；
- 某 genotype pattern 在樣本中的觀察頻率。

它沒有直接提供：

- 某個 clone 曾沿 `parent→child` 邊發生轉換；
- 兩個突變的時間先後；
- 某個相鄰 Hamming-1 transition 被穿越多少次。

精確映射應為：

```text
partial/full read  → genotype node 或相容 node group
read multiplicity  → observed pattern abundance / uncertainty
tree edge          → 由跨 read、結構假設與演化模型推定的 transition
```

### 1.3 partial read 是 subcube constraint，不是 subcube-edge support

形式化規格只要求 partial observation 的相容子立方體與最終 node set 相交：

\[
V(T)\cap G_p\ne\varnothing.
\]

它不要求整個 face 入樹，也不表示 face 內的邊都得到支持。規格來源：
`InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md:84`。

## 2. `k=3` 的最小反例：為何「二點 read 支持所有相應邊」會錯

用 `A/R/X = 1/0/?`，三個 sSNV 為 A、B、C。

### 2.1 `AAX = (1,1,?)`

若 C 位於 solver universe，這條 read 的相容完整 states 是：

\[
G_{AAX}=\{AAR,AAA\}=\{110,111\}.
\]

兩 states 之間的唯一 cube edge 是：

```text
AAR → AAA
```

但這條 edge 翻轉的正是 read **未觀察的 C**。所以 `AAX` 對這條 edge 是「沒有資訊」，不是「支持一次」。把 `n(AAX)` 全加到該邊，等於把 missingness 錯當演化證據。

同時，`AAX` 也無法區分：

```text
ROOT → ARR → AAR    （A 先、B 後）
ROOT → RAR → AAR    （B 先、A 後）
```

把 `AAX` 重複 10、100 或 1,000 次，只會增加「A 與 B 可共現」的信心，不會產生 A/B 先後資訊。

### 2.2 `AAA × 100` 仍不能決定三個突變的順序

實跑現行 solver：

```text
AAA_x100 {'universe': [0, 1, 2], 'n_hidden': 2, 'n_trees': 6}
```

100 條 read 都支持 terminal node `AAA`，但三個 mutation 有 `3! = 6` 種 unit-flip order。`counts_used_in_objective=False` 是由 code audit 得到的判讀，不是 solver return schema 的欄位。這是後續任何新模型都必須保留的 negative control：若只靠這 100 條 reads 選出唯一順序，新方法就是在製造資訊。

### 2.3 更強的不變式

若兩棵候選樹 `T1,T2` 有完全相同的 genotype node set，只差 parent edges：

\[
V(T_1)=V(T_2),\quad E(T_1)\ne E(T_2),
\]

則任何只依 cross-sectional genotype-state mixture 的 likelihood 都應滿足：

\[
P(D\mid T_1)=P(D\mid T_2).
\]

若程式對這兩棵樹給不同 read likelihood，來源只能是：

1. edge attribution／double-count leakage；或
2. 額外 evolutionary branching prior。

第 2 項可以研究，但必須明寫是 prior，不可稱為 read 直接證據。

## 3. 對原始想法的九個必要質疑

| 質疑 | 為何會影響結論 | 必須怎麼驗證 |
|---|---|---|
| count 的方向 | Steiner 通常最小化 cost；直接令 `cost=count` 會讓支持越多越不易選 | 只能用明確 probability model 形成 `−log likelihood`，不可直接用 raw count |
| node vs edge | read 觀察 state，不觀察 transition | `AAX`、`AAA×100` 與 same-V/different-E tests 必須保持 tie |
| coverage opportunity | read 多可能只是區域 coverage 高／span 短 | 對 observed mask 分層，以 callable opportunity 正規化 |
| path-length bias | 同一 child-bearing read 若灌給所有 ancestor edges，深樹被重複計分 | 每 molecule 只進一次 generative likelihood |
| hidden／extinct ancestor | 沒觀察到 intermediate node 不等於不存在 | 允許 hidden node abundance `π=0`，另報 structural support vs observed support |
| low-frequency truth | 真 subclone 本來可能只有少 reads | 用 binomial/beta-binomial uncertainty，不以低 count 直接判假 |
| VAF double counting | 現行 read-AF 來自同一批 BAM reads | 同一 read 不可同時進 joint likelihood 與額外 VAF likelihood |
| CN／LOH／purity／multiplicity | read AF 不等於 CCF，祖先 AF 也不一定大於後代 AF | discovery 先限 CN-neutral/non-recurrence；其他 strata 分開標 confounded |
| HP／alignment grain | HP 可能錯分；現行 counts 是 alignment key，不保證 unique molecule | primary-only/QNAME collapse sensitivity；HP error/missing strata |

## 4. 現有資料能用，但尚不能直接拿來相加

### 4.1 `subread_groups_by_hp` 目前包含 full-cover reads

在 upstream 程式中，每個有至少一個觀察位點的 alignment 都先進 `subread_fam`；若恰好 full-cover，又再進 `combo_fam`：

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py:149`
- 同檔 `:164`：所有 pattern 進 subread count。
- 同檔 `:171`：full-cover pattern 再進 population count。

Frozen HCC1954 真實例：

```json
{"chrom":"chr1","start":1625526,"end":1659152,"family":"1",
 "pattern":"ARR","population_count":4,"subread_count":4}
```

因此 `populations_by_hp + subread_groups_by_hp` 不能直接相加，否則同 4 個 full-cover observations 被算兩次。第一版 likelihood 應使用一份 mutually exclusive pattern table，或回到 MINREAD 前的 per-molecule observations。

### 4.2 現行 pattern table 有 selection truncation

`MINREAD=3` 後才保留 pattern；低頻 pattern 被丟掉。這對 presence-only coverage 是既有設計，但對 probability likelihood 會造成截尾偏誤。正式模型需要：

- MINREAD 前的所有 R/A/X patterns；
- observed mask；
- base quality／MAPQ；
- HP tag 或 HP posterior；
- stable molecule identity；
- alignment class。

### 4.3 count 的 grain 是 alignment exposure，不一定是 molecule

`per_read_calls` 以 `(QNAME, chrom, start, end, FLAG, CIGAR digest)` 為 key，用來防止 primary／supplementary segments 被拼成一個虛構 vector；但同一 QNAME 的不同 alignments 仍可成為多筆 exposure：

`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_linkage_genomewide.py:91`、`:128`。

因此 likelihood 的獨立單位需先定案為 molecule/QNAME 或 alignment，並做 primary-only／QNAME-collapse sensitivity。

## 5. 建議模型：不是 weighted edge，而是 candidate-tree read-pattern likelihood

### 5.1 保留現行 exact candidate generator

第一版不改 `enumerate_min_trees`：

1. 明確呼叫 `enumerate_min_trees(..., tree_cap=0)`，以現行 coverage + min-hidden 規則枚舉完整候選集；
2. 候選 `T`、`Topo`、candidate digest 全部不變；
3. 新增獨立 `L1b_read_pattern_score` annotation；
4. 報 unique top、tied same Topo、tied different Topo、near-tie 與 abstain。

Evaluator 必須在評分前 hard-assert：

```python
assert not result["capped"]
assert result["trees_complete"] is True
assert result["n_trees"] == len(result["trees"])
```

任一失敗就不評分；不可把 solver 預設的前 32 棵 stored prefix 誤當完整候選集。

這可把「結構可行性」和「資料相對支持度」分開，避免一開始就破壞現行 identifiability contract。

### 5.2 Whole-tree likelihood

對候選樹 `T_h`，node `v` 的 genotype 為 `g_v∈{0,1}^k`，其 mixture proportion 為 `π_v`。read `r` 覆蓋 mask `M_r`，觀察 alleles `y_r`。第一版在每個 observed HP family `h` 內條件化建模：

\[
P(r\mid T_h,\pi,h_r=h)=
\sum_{v\in V(T_h)}\pi_{v\mid h}
\prod_{j\in M_r}P(y_{rj}\mid g_{vj},Q_{rj},MAPQ_r,M_r,h_r=h).
\]

HP uncertainty 不塞進每個 locus emission。第一版以 HP-missing／reassignment sensitivity處理；後續若要顯式建模，才改成在 true family `h` 上求和，加入 read-level `P(h_r|h,η_r)`，並對整個 HP forest 聯合評分。

對固定候選樹：

\[
\ell(T)=
\max_{\pi\in\Delta}
\sum_r\log P(r\mid T,\pi).
\]

相同 mask／pattern／quality-bin 的 reads 可先聚合：

\[
\ell(T)=\max_{\pi\in\Delta}
\sum_g n_g\log\left[
\sum_{v\in V(T)}\pi_vP(g\mid v)
\right].
\]

partial read 的 missing bits 由 emission 中自然 marginalize，不必把它指派給某條 edge。

上式定義的是 **profile likelihood**，尚未定義 `P(T|D)`。第一版因此使用經 molecule／parametric bootstrap 校準的 likelihood support set：

\[
\mathcal S_{1-\alpha}=
\left\{T:2\left[\ell(T_{\max})-\ell(T)\right]\le c_{1-\alpha}^{\mathrm{boot}}\right\}.
\]

`c_boot` 必須由 synthetic truth／bootstrap校準，不能直接假設一般卡方近似。只有未來另行明定 tree prior 與 `π` prior，並計算

\[
P(T\mid D)\propto P(T)\int P(D\mid T,\pi)P(\pi)\,d\pi,
\]

才可使用 posterior、credible set 或 Brier score 等 Bayesian 用語。

### 5.3 VAF 是否已包含？

若所有 full／single／partial masks 都放入同一份 mutually exclusive per-molecule table，則 HP family `h` 內、在位點 `j` callable 且通過相同 QC 條件下：

\[
E[\mathbf 1(Y_{rj}=ALT)\mid j\in M_r,h,QC]
=\sum_v\pi_{v\mid h}P(Y_{rj}=ALT\mid g_{vj},j\in M_r,h,QC).
\]

只有在 CN-neutral、error-free或已校正、mask/callability 與 genotype state條件獨立時，才簡化成 `Σ_v π_v g_vj`。在相同 HP、QC、read universe與 denominator 下，read-AF 才是 joint pattern distribution 的 marginal summary；joint data 理論上至少不比該 marginal 少資訊。這是新方法可能比現行 VAF ordering 更有效的理由。

但有三個界線：

1. 若只使用 spanning reads，並未包含 non-spanning reads 的完整 marginal AF；single-locus reads 應直接放進同一 observed-mask likelihood，不另建重疊 VAF term。
2. 若 joint likelihood 已含同一批 reads，再加入同 BAM 計算的 VAF score就是 double counting。
3. HP-specific read-AF 不等於 global caller AF；若要比較，必須先統一 filter、callability與 denominator。

### 5.4 為何它不等於「找一條最高權重路徑」

`log Σ_v π_v P(r|v)` 會同時耦合多個 nodes、所有 reads 與 `π`，通常不能分解成互相獨立的 edge costs。因此一般情況下：

```text
argmax whole-tree likelihood ≠ shortest / maximum-weight path on independently weighted edges
```

若未來要導出 additive edge approximation，必須另外證明 factorization，且以 full likelihood 作校準基準。

## 6. 與現行 VAF ordering 的精確比較

現行程式的 exact score 是：對每條 `parent→child`，找 child 新增 mutation，將 parent 已有的每個 ancestor mutation 的 read-AF 減去 child mutation read-AF後加總：

`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_top_tie_census.py:100`、`:116`。

它有三個已知限制：

1. 使用 per-site marginals，沒有直接使用 read 內 joint pattern frequency；
2. 不同 candidates 的 ancestry comparisons 數量可能不同，raw score-sum 可比性未證；
3. raw read-AF 未校正 purity／CN／multiplicity，selection 不等於 confirmation。

因此 read-pattern likelihood **有合理機會**優於 VAF heuristic，但不能以「Topo>1 變少」當成功。真正成功是：

- truth 排名更準；
- held-out reads 預測更好；
- 無法識別的案例仍保持 tie／abstain；
- 不製造過度自信。

## 7. 四個方法的預註冊比較

| ID | 方法 | 角色 | 預期 |
|---|---|---|---|
| M0 | 現行 exact read-AF/VAF ordering | baseline | 可排序部分 candidates，但受 marginal/CN/score-length 限制 |
| M1 | raw／normalized edge-count heuristic | negative-control ablation | 預期在 missingness、coverage、same-V/different-E 反例失敗 |
| M2 | co-covered pair `00/10/01/11` relation composite likelihood | 中間 baseline | 可利用 joint pair evidence；重疊 reads 需 cluster bootstrap，不能稱 direct edge likelihood |
| M3 | full partial-read candidate-tree mixture likelihood | 主提案 | 可解析 node-set/completion ambiguity；對 pure parent-edge ambiguity應保持 tie |

M3 第一版只估固定或簡化的 error model；後續完整版才加入 CN／purity／multiplicity、context-specific base-error、HP switching與branching prior。不能在第一版一次混入，否則無法知道增益來自哪裡。

## 8. Step → Verify：bounded pilot 計畫

### Phase 0 — 已完成的 contract／反例檢查（<30 min）

1. 執行 solver golden cases  
   → **驗證**：8/8 cases PASS；V1–V7 全通過。
2. 執行 `AAA×100`  
   → **驗證**：`n_trees=6`，count 不得創造順序。
3. 執行 `{AAX,XXA}`  
   → **驗證**：`universe=[0,1,2]`、`n_trees=10`；partial face 只作 node-group constraint。
4. 查 frozen JSON 的 full/subread overlap  
   → **驗證**：HCC1954 `chr1:1625526-1659152 fam1 ARR` 在兩表皆為 4，證明直接相加會 double count。

### Phase 1 — synthetic falsification suite

預計輸入／輸出：

- Input：程式產生的 `k=3–6` truth tree、clone proportions、read masks、ONT-like error／HP switching／coverage。
- Proposed script：`InterSubMod/research/20260712_read_weighted_hypercube_tree/scripts/evaluate_read_pattern_likelihood.py`。
- Proposed output：
  - `InterSubMod/research/20260712_read_weighted_hypercube_tree/data/synthetic_cases.tsv`
  - `InterSubMod/research/20260712_read_weighted_hypercube_tree/data/model_comparison.tsv`
  - `InterSubMod/research/20260712_read_weighted_hypercube_tree/data/verification.json`

1. 建 chain／sibling／hidden／partial／recurrence truth cases  
   → **驗證**：truth candidate 必在現行 candidate set；candidate digest before/after 100% 相同。
2. `AAX` only、`AAA×100`、same-V/different-E  
   → **驗證**：M3 false unique-call ≤5%；same-V/different-E base likelihood exact tie。
3. 加 `ARR` intermediate reads  
   → **驗證**：`ROOT→ARR→AAR` 的 node-set likelihood 應上升。
4. 改加 `RAR` reads  
   → **驗證**：相反 node-set 的 likelihood 應上升。
5. Uniform depth ×2/×5 與 locus-specific coverage distortion  
   → **驗證**：uniform scaling 不改 winner；M1 應在 distortion 下劣化。
6. 加 error／HP switching／CN/LOH  
   → **驗證**：各 confound 的性能下降與 calibration 都被量化，不准只報 aggregate。

### Phase 2 — strict held-out HCC1395 bounded probe

第一批只納入：

- CN-neutral；
- HP1/HP2 primary family；
- candidate-complete、non-capped；
- non-recurrence；
- `k=3/4`；
- `T>1`，並分 `Topo=1/Topo>1`、`C=0/C>0`；
- 最多 200 regions，明示 partial scope。

資料分割以 stable QNAME hash，不能以 alignment row 隨機切：

```text
50%：只建 candidate set
25%：fit π/error、選 tree／likelihood support set
25%：凍結 candidate、tree 與所有參數，只計 held-out predictive NLL
```

做 5 次 cross-fitting／rotation。若測試 reads 參與 candidate generation，就不是嚴格 held-out，必須另標 transductive sensitivity。

Step → Verify：

1. 建包含 full／single／partial masks 的 mutually exclusive MINREAD-before pattern table  
   → **驗證**：每 molecule 在每 region×HP 只計一次；pattern counts 回到 callable molecule denominator。
2. 對相同 candidate set 算 M0–M3  
   → **驗證**：所有方法都呼叫 `tree_cap=0`，且 hard-assert `not capped`、`trees_complete`、`n_trees==len(trees)`；`T/Topo/candidate digests` 完全不變，任何 mismatch 直接 FAIL。
3. 以 fitting partition 選 tree／bootstrap-calibrated likelihood support set；held-out partition 禁止重估 `π`、error 或門檻  
   → **驗證**：M0–M3 選出的結果一律用同一個 frozen observation model算逐 region predictive NLL，並同報 winner、abstain、read depth、mask、CN、HP。
4. molecule bootstrap／downsample  
   → **驗證**：identifiable unique cases 報 top stability；tied cases 報 support-set Jaccard／relation stability；support-set coverage 有 95% CI，不只看 exact tie。

### Phase 3 — 技術重現與後續全量門檻

1. HCC1395 vs HCC1395_DORADO matched `region×HP`  
   → **驗證**：只比較共同 candidate universe/digest；固定相同 call rate或畫 selective agreement–coverage curve，並同報 both-called `n`、agreement、abstention、unmatched與excluded。這是同 biological sample 的 platform reproducibility，不當成兩個生物複本。
2. 通過後才升 Type B，擴 chr1–22 × 全 dataset rows  
   → **驗證**：每樣本同報 aggregate、canonical case、extreme outlier、well-explained case；所有母數守恆。

## 9. 指標與預註冊判定門檻

### 9.1 主要指標

- exact labeled tree accuracy；
- canonical Topo accuracy；
- node-set/completion accuracy 與 parent-edge accuracy分開；
- true tree／Topo 是否位於 bootstrap-calibrated 90%／95% likelihood support set；
- mean reciprocal rank；
- held-out per-read negative log-likelihood；
- support-set empirical coverage／width與 selective accuracy–coverage curve；
- identifiable unique cases 的 top stability；tied cases的 support-set Jaccard／relation stability；
- unique／tied-same-Topo／tied-different-Topo／abstain；
- false uniqueness in known non-identifiable cases；
- runtime／non-convergence。

### 9.2 暫定 GO gate（須在跑結果前鎖定）

M3 才能升為 GO，必須同時滿足：

1. identifiable synthetic cases 在相同 call/abstention rate 下，Topo accuracy 比 M0 至少 +10 percentage points；paired-bootstrap 95% CI lower bound >0；
2. bootstrap-calibrated nominal 95% likelihood support-set coverage 落在 90–98%；
3. non-identifiable controls false unique-call ≤5%；
4. identifiable／declared-unique cases 的 cross-fit top stability ≥0.80；不可識別／tied cases改看 support-set Jaccard與 coverage；
5. HCC1395/DORADO 只在相同 matched `region×HP`、共同 candidate digest與相同 call rate下比較；agreement 不低於 M0且改善至少 +5 percentage points，並完整列 both-called／abstention／excluded分母；
6. 同方向結果在 coverage／read length／MAPQ／HP strata 不翻轉。

這些是 pre-registered engineering thresholds，不是領域既有公認標準；pilot 後需做 calibration review。

### 9.3 任一出現即 NO-GO／退回模型

- `AAX` 或 `AAA×100` 被判唯一順序；
- same-V/different-E 被 base read likelihood強行分開；
- in-sample 有增益，但 held-out 增益消失；
- 增益只由 coverage／read length／HP／CN confound 驅動；
- bootstrap likelihood support set 嚴重 under-cover；
- candidate set 被 scorer 靜默刪減；
- 同一 molecule 被 joint pattern、pair terms與 VAF 重複計分。

## 10. Assumption Map

| Assumption | Importance | Known? | 處置 |
|---|---|---|---|
| read pattern 比 marginal VAF 多保留 joint information | High | Known mathematically | 仍需 truth simulation證明能轉為 topology accuracy |
| partial read 可直接支持 face edge | High | **False** | 原始 edge-count formulation NO-GO |
| pattern counts 是 independent molecules | High | Unknown／現行為 alignment-key | P0：QNAME／primary-only sensitivity |
| missing mask 可視為 ignorable | High | Unknown | 以 mask conditional likelihood + span/read-length strata 驗證 |
| raw read-AF 可當 CCF | High | **False** | CN/purity/multiplicity correction；先限 neutral stratum |
| M3 能辨識 same-V/different-E | High | **False under basic model** | 必須保持 tie；若加 branching prior，分開報 prior-driven result |
| HCC/DORADO 是獨立 biological replication | High | **False** | 只作技術重現 |
| current solver candidate set complete | High | Known only for non-capped verified units | pilot 僅納 candidate-complete units |

## 11. Credibility score

本 score 評的是「重新形式化的 M3 pilot」，不是已被反例否決的 raw edge-count 版本。

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 10/20 | mixture/read likelihood合理；直接 edge interpretation 不成立 |
| 觀察支撐 | 10/20 | 真實 joint pattern counts存在，但未有 orthogonal truth／held-out validation |
| 機制清晰度 | 10/20 | node-emission DAG清楚；CN/HP/missingness仍待建模 |
| 反例風險 | 0/20 | 已知 same-V/different-E、coverage、CN、alignment duplication 等高風險；conflict scan扣 10 |
| 所需資源 | 10/20 | hand tests <30 min；合格 evaluator＋bounded probe約 1–6 h，不含 raw data重產 |
| **TOTAL** | **40/100** | **PROBE** |

Falsifier observable：若 M3 在 `AAX`／`AAA×100`／same-V-different-E 資料上製造唯一拓撲，或其 synthetic／held-out增益消失，核心假說即被否證。

## 12. Evidence conflict scan

| Prior evidence／decision | Relation | 判讀 |
|---|---|---|
| Formal spec §7.5.3 明定 enumerate/non-rank | Conflict | 新 scorer只能先作外掛 annotation；若要升 canonical，必須先修 spec，不可靜默改 objective |
| Formal spec 舊 §5 曾寫 max observed support | Internal inconsistency | 文件與 code不同步；現行 code未實作 support tie-break |
| `20260709_ccf_readcount_tree_weighting_demo_01.md` 已更正為 historical family read-AF heuristic | Dependent | 舊 66.8%「tie 被破」不能當真值 accuracy，selection ≠ confirmation |
| `analyze_read_likelihood.py` 的 `Σlog(parent reads+1)` | Negative-control candidate | 名稱雖為 likelihood，實際不是 normalized generative likelihood；不可直接採用 |
| 20260711 VAF audit | Support／ceiling | exact ranking可作 descriptive ordering，但 raw AF非 CCF、未做 molecule bootstrap、不可稱真樹確認 |
| repo root `MEMORY.md` 不存在 | Gap | 以 validated reports／method docs／recent git log替代 conflict scan；後續應補 SoT 路由 |

## 13. 原始文獻對齊

| Primary source | 支持什麼 | 不支持什麼 |
|---|---|---|
| [Pe'er et al., Incomplete Directed Perfect Phylogeny](https://doi.org/10.1137/S0097539702406510) | 在 species×character matrix 情境中，missing binary states可作 incomplete directed perfect phylogeny | 與 bulk ONT partial reads僅為形式類比；不提供 read count 作 transition edge count |
| [TreeClone, Annals of Applied Statistics](https://doi.org/10.1214/18-AOAS1224)／[PairClone](https://doi.org/10.1111/rssc.12328) | mutation-pair／same-read joint information可比單點 marginals提供更多 subclone資訊 | 沒有把一條 partial read灌給所有相容 hypercube edges |
| [GCtree](https://doi.org/10.1093/molbev/msy020) | 可用明確 branching-process abundance likelihood排序 equally parsimonious trees | genotype abundance與bulk read depth不是同一量；需要生成模型，不是 raw count |
| [SCITE](https://doi.org/10.1186/s13059-016-0936-x) | noisy/missing observations應對整棵 tree計算 ML/MAP likelihood | 不支持獨立 raw edge加權 |
| [SCIΦ](https://doi.org/10.1038/s41467-018-07627-7) | read counts與coverage可用 beta-binomial聯合建模 | 其 single-cell資料模型不可原封不動套到bulk ONT |
| [PhyloSub](https://doi.org/10.1186/1471-2105-15-35) | SNV frequencies、clone proportions與tree uncertainty應聯合推論並保留多解 | 不保證單一 bulk sample能唯一辨識 topology |
| [BAMSE](https://doi.org/10.1186/s12859-019-2824-3) | 癌症 read counts可在含error的 Bayesian model中排序候選 subclone trees | 沒有利用 long-read partial multi-locus subcube作直接 edge traversal |

截至 2026-07-12 的 primary-source 查核中，未找到「把 partial ONT read 所定義的 Boolean subcube multiplicity，直接當 Hamming-1 evolutionary edge cost」的正式先例。較一致的文獻路線是 observation likelihood、node/genotype abundance model、whole-tree profile likelihood 或 posterior。

## 14. Decision path

- **NO-GO**：直接以 raw read count 當 edge cost，或聲稱 read 直接支持 evolution edge。
- **PROBE**：保留完整 min-hidden candidate set，新增 M2/M3 read-pattern scorer，執行 synthetic falsification + strict held-out bounded probe。
- **GO only if**：§9.2 全部通過，且獨立 red-team／data QA再次確認。
- **Decision lock**：在 §9.2 通過前，不改現行 solver objective、不刪候選、不把 unique top寫成 confirmed true topology。

教授版一句話：

> **read pattern 可以提高候選樹的相對支持度，但 read 沒有直接觀察演化邊；因此應以節點／整樹 likelihood 排候選並保留並列，而不能把 read count 當成 edge traversal count。**

## 15. 本輪實際執行與輸出

### Inputs

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py`
- `InterSubMod/research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_top_tie_census.py`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1954/mlhp_part_1.json`

### Commands and observed fragments

```bash
cd /big7_disk/liaoyoyo2001
python3 InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py
# observed: GOLDEN 總結: ALL PASS ✓

cd /big7_disk/liaoyoyo2001/InterSubMod
python3 -c 'import sys; sys.path.insert(0,"docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"); import tree_enumeration_solver as S; a=S.enumerate_min_trees({"AAA":100},[],3,tree_cap=0); b=S.enumerate_min_trees({},["AAX","XXA"],3,tree_cap=0); print("AAA_x100", {"universe":a["universe"],"n_hidden":a["n_hidden"],"n_trees":a["n_trees"]}); print("AAX_XXA", {"universe":b["universe"],"partial_subcubes":[sorted(S.label(x,set(),3) for x in q) for q in b["_partial_subcubes"]],"n_hidden":b["n_hidden"],"n_trees":b["n_trees"]})'
# observed: AAA×100 n_trees=6; AAX/XXA n_trees=10

jq -c 'first(.groups[] as $g | (($g.populations_by_hp // {}) | to_entries[]) as $fam | ($fam.value | to_entries[]) as $pat | select((($g.subread_groups_by_hp // {})[$fam.key][$pat.key] // null) != null) | {chrom:$g.chrom,start:$g.start,end:$g.end,family:$fam.key,pattern:$pat.key,population_count:$pat.value,subread_count:(($g.subread_groups_by_hp // {})[$fam.key][$pat.key])})' \
  /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1954/mlhp_part_1.json
# observed: ARR population_count=4, subread_count=4
```

### Output

- `InterSubMod/research/20260712_read_weighted_hypercube_tree/pre-decision-audit.md`
- 未產生 solver patch、未產生新 biological result、未啟動長計算。

## 16. 獨立分視角審查

| Reviewer 視角 | Final readback | 結論 |
|---|---|---|
| Solver contract／k=3 red-team | **PASS** | raw edge support語義不成立；`tree_cap=0`＋完整候選 hard asserts；外掛 whole-tree scorer |
| Statistical likelihood／identifiability | **PASS** | M3可解析 node-set ambiguity，但 same-V/different-E 必須保持 tie；profile likelihood、cross-fit與support-set定義一致 |
| Primary literature／novelty | **PASS** | joint read evidence有先例；直接 subcube-count→edge cost未找到先例；文獻與claim ceiling正確 |

三視角 final readback 為 **3/3 PASS**；因此文件可作 bounded pilot 的前驗依據，但研究 verdict 仍是 **PROBE**，不是 GO。

## 17. 2026-07-12 addendum：compatible-node 流入邊均分法

> **新增提案**：將 pattern count 均分至所有相容 nodes 的流入 Hamming-1 edges，讓全 cube 的流入質量等於 mutation-bearing read 數，再以 tree 捕捉的 edge weight 判斷較可能結構。

> **新增裁決**：質量守恆成立，但不足以賦予 edge evidence 語意。`all-inflow` 版本 **NO-GO as main score**；排除 X-flip 後的 `observed-ALT boundary` 版本只保留為 M1 negative-control／relation heuristic；candidate-conditioned node responsibility仍是首選 PROBE。

### 17.1 均分規則與「守恆」究竟證明了什麼

令 partial pattern `p` 的 count 為 `n_p`，相容 nodes 為：

\[
G_p=\{x:x\text{ 與 }p\text{ 的所有 observed bits 相容}\}.
\]

所有相容 node 的流入邊集合：

\[
I_p=\bigcup_{x\in G_p}\delta^-(x).
\]

若採「所有 edges 平分」：

\[
a_{p,e}=\frac{\mathbf 1[e\in I_p]}{|I_p|},\qquad
w_e=\sum_p n_pa_{p,e},
\]

此式只定義於 mutation-bearing pattern 且 `|I_p|>0`；`q=0` 的 REF-only patterns 必須保留在 root/node likelihood，不可硬塞入 mutation-edge weighting。

則每個 pattern 確實滿足：

\[
\sum_e n_pa_{p,e}=n_p.
\]

若採「先均分 compatible nodes，再均分該 node 的 inflow」，也同樣守恆。兩者差別只是每條 edge 分到多少，核心問題不變：**守恆只證明一條 read 在整個 cube 的總帳仍算一次，沒有證明它被分到資料支持的 transition。**

### 17.2 `AXX`／`AAX` 的精確 edge census

bit order 為 `(A,B,C)`。

| Pattern | Compatible nodes | 所有 inflow | 翻轉 X 的錯誤 edges | flat-edge 錯誤質量 | 排除 X 後 boundary edges |
|---|---|---:|---:|---:|---:|
| `AXX` | `100,101,110,111` | 8 | 4 | 4/8 = **50%** | 4 |
| `AAX` | `110,111` | 5 | 1 | 1/5 = **20%** | 4 |

`AXX` 的 8 條所有 inflow：

```text
000→100  flip A   ✓ observed ALT
001→101  flip A   ✓ observed ALT
100→101  flip C   ✗ X
010→110  flip A   ✓ observed ALT
100→110  flip B   ✗ X
011→111  flip A   ✓ observed ALT
101→111  flip B   ✗ X
110→111  flip C   ✗ X
```

`AAX` 的 5 條所有 inflow：

```text
010→110  flip A   ✓ observed ALT
100→110  flip B   ✓ observed ALT
011→111  flip A   ✓ observed ALT
101→111  flip B   ✓ observed ALT
110→111  flip C   ✗ X
```

所以如果「排除沒有資訊支持的 edge」包含排除所有 X-flips，方向是正確的；但它只修掉最明顯的錯誤，尚未使 edge weight 成為 likelihood。

若採「先均分 nodes、再均分該 node inflow」而不是 flat-edge 均分，結論也不變：`AAX` 先給 `110/111` 各 `n/2`，再分 inflow後，未觀察 C 的 `110→111` 仍取得 `n/6` 假支持。

### 17.3 實算：守恆仍可產生假順位

令 `n(AXX)=80`、`n(AAX)=50`，比較：

```text
T_A：000→100→110→111  （A→B→C）
T_B：000→010→110→111  （B→A→C）
```

flat all-inflow 均分結果：

| Candidate | AXX 貢獻 | AAX 貢獻 | 總分 |
|---|---:|---:|---:|
| `T_A` | 3/8×80 = 30 | 2/5×50 = 20 | **50** |
| `T_B` | 2/8×80 = 20 | 2/5×50 = 20 | **40** |

`AXX` 並沒有觀察 B 或 C，理論上不能由它判 A/B order；`50 vs 40` 是 candidate 含幾條被均分 edges 的 cube-geometry artifact。

只保留 observed-ALT flips 後：

| Candidate | AXX 貢獻 | AAX 貢獻 | 總分 |
|---|---:|---:|---:|
| `T_A` | 1/4×80 = 20 | 1/4×50 = 12.5 | **32.5** |
| `T_B` | 1/4×80 = 20 | 1/4×50 = 12.5 | **32.5** |

這次保持 tie 是正確行為，但也顯示均分本身沒有提供 order 的新資訊。

另一個 invariant：`AA` 在 `k=2` 時只有 2 條 inflow，每邊 `1/2`；加入完全未觀察的 C 變成 `AAX`，所有 inflow 變 5 條，每邊 `1/5`。資料對 A/B 沒有改變，權重卻因多一個 X 維度而稀釋，稱為 missing-dimension bias。

### 17.4 排除 X-flip 後，仍不能同時保存 read mass 與 VAF mass

令 observed ALT 數為 `q`、missing 數為 `m`。排除 X-flip 後的 boundary edges 數為：

\[
|B_p|=q\,2^m.
\]

若每邊分 `n_p/(q2^m)`：

- 全 edge read mass 合計為 `n_p`；
- 但每個 observed ALT mutation 的總量只有 `n_p/q`，不再等於該 mutation 的 ALT numerator；
- `q=0` 的 `RRX/RRR` 沒有 boundary edge，reference denominator全部丟失。

若改成每邊 `n_p/2^m`，讓每個 observed ALT 都保存 `n_p`：

- edge mass 合計變成 `q·n_p`；
- multi-ALT read被依 mutation load 重複計算。

因此「mutation-bearing read mass守恆」與「逐位點 VAF numerator守恆」不是同一件事。edge table 本身也不是原始 pattern table 的充分統計量，不能宣稱已完整包含 VAF。

### 17.5 為何 REF／all-reference reads 不能從模型排除

`RRR`／`RRX` 不會直接支持 mutation acquisition edge，但它們提供：

- 每位點 callable REF denominator；
- root／ancestral／normal-background proportion；
- 判斷 low-frequency ALT 是真 subclone 或 error 的基準。

尤其 `RRX` 仍可能相容 `000/001`；只因 observed 部分皆 REF 就排除，會把 unknown state 誤當確定 root。可在「mutation-bearing edge-width視覺化」中不畫 root mass，但不可從 likelihood／VAF denominator刪除。

### 17.6 與現行 VAF 的差別

| 問題 | 現行 VAF ordering | Uniform inflow edge | Whole read-pattern likelihood |
|---|---|---|---|
| 使用的基本量 | 每位點 `ALT/callable` marginal | pattern count投影後的人工 edge mass | full/single/partial molecule patterns |
| 是否保留同 read 共現 | 否 | 部分，但被 uniform completion壓縮 | **是** |
| X/missing 處理 | 位點 denominator不含 X | 容易產生 X-flip或 `2^m` dilution | 對相容 nodes marginalize |
| REF denominator | **有** | mutation-only規則會丟失 | **有** |
| 是否為 calibrated probability | 否，現行為 ordering heuristic | **否** | 可建立 profile likelihood＋bootstrap support set |
| 可否直接證明 edge | 否 | 否 | 否；只解析資料可識別的 node-set/completion |
| CN／purity敏感 | 高 | 高 | 仍需明確建模／分層 |

VAF 只看 marginal，因此可能漏掉 joint structure。兩組資料可以有相同 VAF：

```text
D_linked：50 AAX + 50 RRX  → AF_A=AF_B=0.5，A/B 共現
D_sister：50 ARX + 50 RAX  → AF_A=AF_B=0.5，A/B 互斥
```

VAF 無法區分兩者；read-pattern likelihood可以。可是 uniform edge projection不是保留 joint table的必要方法，且可能丟掉 denominator與製造 geometry bias。

### 17.7 若要保留「read 流入邊」的直觀，建議改成 candidate-conditioned responsibility

對每棵候選樹 `T`，先把 read 分配到相容 latent nodes：

\[
\gamma_{rv}^{(T)}=
\frac{\pi_vP(y_r\mid v,M_r,Q_r)}
{\sum_{u\in V(T)}\pi_uP(y_r\mid u,M_r,Q_r)},
\qquad \sum_v\gamma_{rv}^{(T)}=1.
\]

每棵候選樹都必須以完全相同規則重新最佳化 `π`（或事先明確固定同一組 `π`）；emission model須有 nonzero error floor，確保所有 read 的 responsibility denominator大於零。

node 的 expected read mass：

\[
N_v^{(T)}=\sum_r\gamma_{rv}^{(T)}.
\]

此時 `Σ_v N_v=N_reads`。root mass須另列；真正畫在非 root incoming edges 上的總量是 `N_reads−N_ROOT`。若 node `v` 在 candidate `T` 中有唯一 parent，可把 `N_v` 畫在 `parent_T(v)→v` 的 edge width。這符合使用者想要的「流入質量守恆」直觀，但其正確名稱是：

> **expected reads assigned to the child genotype state**，不是 observed transition count。

比較候選樹時，不能比較 `Σedge N_v`，因為總質量受守恆約束；必須比較：

\[
\ell(T)=\sum_r\log\sum_{v\in V(T)}\pi_vP(y_r\mid v,M_r,Q_r).
\]

若兩棵樹 node set完全相同、只差 parent edges，基本 likelihood仍應同分。任何進一步排序必須來自明示的 branching／transition prior或外部時間／多區域／single-cell evidence。

### 17.8 方法優先序與新 variant 決策

| 方法 | 角色 | 裁決 |
|---|---|---|
| 所有 compatible-node inflow 均分 | M1a negative control | **NO-GO 主方法**；會支持 X-flip、產生假順位 |
| 只對 observed-ALT boundary 均分 | M1b ablation／教學圖 | 可測，但仍是 uniform-completion prior，不是 edge likelihood |
| 現行 VAF ordering | M0 baseline | 暫時保留；簡單可解釋，但只用 marginal且受 CN 影響 |
| Candidate-conditioned EM／read-pattern likelihood | M3 主提案 | **首選 PROBE**；理論上同時保留 joint patterns與VAF marginal |

在「現行 VAF vs 均分 edge」二選一時，**VAF 較適合當主 baseline**；均分 edge 只能作 negative-control。真正值得投入驗證的升級不是均分 edge，而是 M3 whole-tree likelihood。

### 17.9 此 variant 的額外 falsifiers

1. **X-edge test**：`AXX` 對 B/C flip、`AAX` 對 C flip必須為 0。
2. **Missing-dimension invariance**：新增完全未觀察的 X 維度，不得單獨改變 A/B conclusion。
3. **Same-V/different-E**：基本 read likelihood必須同分。
4. **Mass dual-ledger**：分開報 read-mass conservation與per-mutation ALT/REF conservation，禁止混稱。
5. **Linked-vs-sister**：在兩種 truth candidates 都被完整枚舉、informative read patterns足夠且 emission已校準時，相同 VAF、不同 joint patterns應可由 M3區分；VAF/M1不得製造錯誤順序。
6. **Reference denominator**：移除 RRR/RRX 後的 winner若改變，證明 mutation-only edge score有 selection bias。

### 17.10 新 variant 分視角 red-team

| Reviewer 視角 | Final readback | 結論 |
|---|---|---|
| Cube/solver contract | **PASS** | all-inflow與observed-ALT boundary均分都不是 transition observation；前者直接支持 X-flip，後者仍有 uniform-completion／recurrence bias |
| Statistical likelihood | **PASS** | 質量守恆只是 bookkeeping；edge table不保存 callable REF／完整 VAF；candidate-conditioned responsibility才有機率語意 |

兩視角 final readback 為 **2/2 PASS**：**均分規則可當對照實驗，不應取代現行 VAF，也不應直接進正式 solver objective。**
