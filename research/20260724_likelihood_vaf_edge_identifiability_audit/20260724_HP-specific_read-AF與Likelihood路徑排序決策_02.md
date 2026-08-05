<!--
建立時間: 2026-07-24
目標: 定義 HP-specific read-AF、joint read-pattern likelihood 與 directed path/topology 的可辨識邊界，形成可供教授與後續實作使用的方法決策
處理範圍: Task Type B comprehensive method validation；exact PS × HP × read-linked component；本文件不宣稱 current strict 7 datasets 已完成 likelihood/topology
關聯檔案: InterSubMod/research/20260724_likelihood_vaf_edge_identifiability_audit/20260724_Likelihood與VAF樹路徑可辨識性驗證_01.md
-->

# HP-specific read-AF 與 likelihood 的路徑排序決策

> **TL;DR**：HP-specific read-AF 比 global VAF 更適合解釋「同一 HP 內的相對盛行率」，但不比 joint read-pattern likelihood 保留更多資料。正式流程應以 likelihood 選 mutation-state vertex set，再以 HP-specific read-AF／校正後 CCF 作「祖先順序相容性檢查」。它可以排除不相容路徑或給出條件式偏好，不能單獨證明 direct parent、hidden clone 或唯一真實 topology。

> **Scope ribbon**：這是方法定義與驗證規格，不是 7 datasets 新版結果。current strict exact-PS 全樣本 likelihood／directed-topology receipt 尚未完成。

## 1. Task type 與服務目標

- Task Type：B — comprehensive validation。
- 服務目標：G1（longphase 家族整合）、G3（read-level 資訊）、G4（跨樣本可重現）、G5（可外部驗證）。
- Thread D：間接相關；本輪處理的是 genetic read-state，不使用 methylation 作主要 tree score。
- Thread B 撤回範圍：否。
- KDE：不適用。
- caller AF：不作主要輸入；主要量為 exact-PS×HP molecule calls。caller VAF 只能作外部比較。
- 長計算／C++：本輪不啟動；先凍結方法語意與驗證 gate。

## 2. 三個量必須分開命名

### 2.1 HP-specific read-AF

分析單位 \(U\) 定義為：

\[
U=\text{dataset}\times\text{chromosome}\times\text{exact PS}
\times\text{HP}\times\text{read-linked component}.
\]

對位點 \(j\)、HP family \(h\)：

\[
q_{j,h}=
\frac{N_{\mathrm{ALT}}(j,h)}
{N_{\mathrm{ALT}}(j,h)+N_{\mathrm{REF}}(j,h)}.
\]

- 分子單位應為去重後的 primary informative molecule。
- ALT／REF 必須使用相同 QC、mapping、base-quality 與 HP/PS 契約。
- X、O、D、S、L 等無法判定的 call 不可硬補成 REF 或 ALT，也不進該位點分母。
- 比較不同位點時，必須處理 denominator/callability 不同與共同 molecule 造成的相關性。

它的正確名稱是：

> **HP-conditioned ALT molecule fraction**  
> 或 **HP-specific read-AF**

它不是：

- 全樣本 global VAF；
- 經 purity／copy number／multiplicity 校正的 CCF；
- 真實 clone fraction；
- 突變發生時間。

### 2.2 Joint read-pattern likelihood

每條 molecule 以 \(r\in\{R,A,X\}^k\) 表示；候選 mutation-state vertex set 為 \(V\)。現行 M2 的核心是：

\[
\ell(V)=
\max_{\pi}
\sum_r n_r
\log\left(
\sum_{v\in V}\pi_vP(r\mid v,\mathrm{BQ})
\right).
\]

它保留：

- 同一 molecule 上的 R/A 共現；
- partial read 的 observed mask；
- 每個 call 的 base quality／error probability；
- 各 pattern 的 molecule count；
- candidate state mixture \(\pi\) 與 optimizer/KKT 證書。

理想無誤差條件下：

\[
q_{j,h}=\sum_{v\in V}\pi_v\,\mathbf{1}[v_j=A].
\]

因此 HP-specific read-AF 是 joint pattern distribution 的 **marginal summary**。若二者來自同一批 molecules，不能把 AF score 再加到 likelihood 當第二份獨立證據。

### 2.3 CN/purity/multiplicity-aware CCF

只有把 normal contamination／tumor DNA fraction、allele-specific CN、LOH、mutation multiplicity 與 sampling uncertainty 納入後，才可把 allele fraction轉成較接近 mutation cellular prevalence 的量。

正式高 claim 應優先使用：

> **帶不確定區間的 HP-compatible cellular prevalence／CCF**

而不是 raw HP-specific read-AF。

## 3. 「AF 越大代表越早」的正確邏輯

### 3.1 可以成立的單向命題

若突變 \(A\) 真的是突變 \(B\) 的祖先，且符合：

1. infinite-sites：無 loss、back mutation、recurrence；
2. 局部 copy-neutral；
3. 無 LOH、WGD 或 allele-specific gain/loss；
4. mutation multiplicity 相同；
5. HP/PS、mapping、callability 與錯誤率可信；
6. 深度足夠，且比較的是同一局部 lineage；

則後代帶有 \(B\) 的細胞集合應包含於帶有 \(A\) 的集合：

\[
\mathcal{C}_B\subseteq\mathcal{C}_A
\quad\Rightarrow\quad
f_A\ge f_B.
\]

在上述條件下，HP-specific read-AF 才可近似保留相同順序：

\[
q_{A,h}\gtrsim q_{B,h}.
\]

因此正確用語是：

> 「較高的校正後盛行率 **支持** A 早於 B，或使 B→A 不相容。」

### 3.2 不可反推的命題

\[
q_A>q_B
\quad\not\Rightarrow\quad
A\text{ 是 }B\text{ 的祖先}.
\]

原因：

- A、B 可能是兩條 sister branches，只是 A 所在 clone 較大；
- 較晚突變可能因正向選擇而快速擴張；
- CN、LOH、multiplicity、purity、depth 或 mapping bias 可改變 AF；
- 不同位點由不同 molecules 覆蓋；
- 即使 A 早於 B，也不能只靠 AF 證明 A 是 B 的 **direct parent**，中間仍可能有 hidden state。

所以「早」只可解釋為：

> 在同一候選 lineage 內，較早取得且未丟失的 mutation 應存在於不小於後代的 carrier population。

它不是絕對發生時間，也不能跨 sister branches 或不相干區域做全域排序。

## 4. likelihood 與 HP-specific read-AF 各自回答什麼

| 問題 | Joint likelihood | HP-specific read-AF／CCF |
|---|---|---|
| 哪些 R/A/X states 能解釋 molecules？ | 主要方法 | 只保留 marginal，較弱 |
| partial reads 如何處理？ | 對 missing bits marginalize | 每位點分母不同，較難處理 |
| BQ／錯誤率如何處理？ | 直接進 emission | raw AF 通常未完整處理 |
| 不同 vertex sets 哪個較符合 reads？ | 可以比較 | 可作診斷，不應重複加分 |
| 同 vertex set 的不同 parent edges？ | edge 不在 emission，通常同分 | 可檢查條件式 prevalence order，但不是獨立識別 |
| 生物敘述是否直觀？ | 較抽象 | 較直觀 |
| 能否證明 direct parent？ | 通常不能 | 不能 |
| 能否證明 hidden clone？ | connector-only state 不能 | 不能 |

結論：

> **read-AF 在解釋上較簡單；likelihood 在資料使用上較完整。兩者不是互相替代，而是不同層級。**

## 5. 為何同一 vertex set 的 edge tie 不能被同批 VAF 真正破解

考慮：

\[
V=\{00,10,01,11\},\quad
\pi=(0.45,0.30,0.15,0.10).
\]

則：

\[
AF_1=0.30+0.10=0.40,\qquad
AF_2=0.15+0.10=0.25.
\]

兩個 parent assignments：

- \(E_1\)：`11` 從 `10` 流入；
- \(E_2\)：`11` 從 `01` 流入。

兩者具有完全相同的 \(V\)、\(\pi\)、read distribution 與 marginal AF，所以 joint likelihood 必須同分，AF 本身也沒有新增觀測。

舊程式的 ordering score 會給：

```text
E1 score = +3/20
E2 score = -3/20
```

它確實能「選」E1，但選擇來自以下先驗：

> 較高 AF 的 mutation 應先取得、較低 AF 的 mutation 應後取得。

因此應標示為：

`HP_AF_CONDITIONAL_PREFERENCE`

不可標示為：

`LIKELIHOOD_VAF_CONFIRMED_EXACT_TOPOLOGY`

## 6. 正式建議：lexicographic pipeline，不直接加權相加

### Step 1 — 候選完整性

枚舉所有 read-compatible、minimum-hidden candidates。

→ 驗證：`candidate_generation_complete=true`；未列完、timeout 或 cap-hit 一律 `ABSTAIN_ENUMERATION`。

### Step 2 — likelihood 選 vertex set

以全部 informative R/A/X molecules、BQ 與 missing marginalization 比較 vertex sets。

→ 驗證：KKT/optimizer certificate、winner tolerance、held-out／bootstrap stability 均有 receipt。

### Step 3 — 保留 winning vertex sets 的全部 parent assignments

不得因同一 snapshot likelihood 無法分 edge 就偷偷任選一棵。

→ 驗證：same-\(V\)/different-\(E\) fixture 必須保持 likelihood tie。

### Step 4 — AF/CCF eligibility gate

只有下列條件通過才做順序檢查：

- exact PS × HP；
- read-linked component；
- CN-neutral、無 LOH/WGD/allele-specific CN confound；
- recurrence 不必要；
- mutation multiplicity 可比；
- coverage/callability 足夠；
- molecule bootstrap 或 posterior uncertainty 可估。

→ 驗證：每個 unit 產出 `AF_EVALUABLE` 或具體不可評估原因。

### Step 5 — prevalence-order compatibility

對每條候選 ancestry edge 檢查：

\[
P(f_{\mathrm{parent}}\ge f_{\mathrm{child}}\mid D)
\]

以及互斥 children 的 sum condition：

\[
f_{\mathrm{parent}}\ge
\sum_c f_{\mathrm{child},c}.
\]

應用 read-level bootstrap 保留共同 molecules 的 covariance；若用 beta-binomial，必須說明它是否忽略跨位點相關性。

→ 驗證：報告每條 edge 的相容機率／穩定率與 uncertainty；不可因 1-read 微差直接稱唯一。

### Step 6 — 結論分級

只有一條 path 在完整候選、likelihood winner set、AF/CCF constraints 與 sensitivity analysis 後仍可行時，才可稱：

`AF_CONDITIONAL_UNIQUE_PREFERENCE`

若還有 single-cell、multi-region/timepoint 或其他獨立 lineage evidence，才可再升高 claim。

## 7. 建議輸出分類

### Candidate generation

- `ENUM_COMPLETE`
- `ENUM_INCOMPLETE_ABSTAIN`

### Vertex-set inference

- `V_UNIQUE`
- `V_TIED`
- `OPTIMIZER_ABSTAIN`

### Edge/path inference

- `E_STRUCTURALLY_UNIQUE`
- `E_MULTIPLE_AF_CONDITIONAL_ONE_PREFERRED`
- `E_MULTIPLE_AF_TIED`
- `E_MULTIPLE_AF_DISCORDANT`
- `E_AF_NOT_EVALUABLE_CNV_LOH`
- `E_AF_NOT_EVALUABLE_RECURRENCE`
- `E_AF_NOT_EVALUABLE_LOW_POWER`

### Topology claim

- `EXACT_TOPOLOGY_UNIQUE_WITHIN_MODEL`
- `COARSE_TOPOLOGY_UNIQUE_EXACT_EDGE_AMBIGUOUS`
- `AF_PREFERRED_TOPOLOGY_ASSUMPTION_CONDITIONAL`
- `TOPOLOGY_AMBIGUOUS`
- `ABSTAIN`

connector-only hidden node必須另標：

- `HIDDEN_FULL_OBSERVED`
- `HIDDEN_PARTIAL_SUPPORTED`
- `HIDDEN_CONNECTOR_ONLY`

`HIDDEN_CONNECTOR_ONLY` 不可直接解釋為已觀察到的 clone。

## 8. 舊 VAF ordering 應如何定位

舊分數：

\[
S(T)=
\sum_{p\rightarrow p\cup\{j\}}
\sum_{i\in p}
(q_i-q_j).
\]

可保留作：

- historical baseline；
- 教授可理解的方向示意；
- sensitivity sidecar；
- candidate path 的 assumption-conditional ordering。

不可直接作：

- calibrated probability；
- posterior tree probability；
- CCF；
- 真實唯一 topology；
- likelihood 的額外加分項。

目前程式另有三個限制：

1. exact Fraction 的任何非零差都可能變成 unique-first；
2. 不同 hidden vertex sets 的 ancestry comparison 數可能不同；
3. raw HP-AF 未校正 CN/purity/multiplicity，且 recurrence-required candidates 不可可靠評估。

若要升級，不應只把 score threshold 改成 0.5；應改成 uncertainty-aware compatibility probability，並用 simulation／known mixture 校準 decision threshold。

## 9. 教授版正式敘述

> 我們先在 exact phase block 與同一 HP 內，利用 long-read 上的 sSNV 共現、partial-read mask、base quality 與 molecule counts，找出哪些 mutation-state combinations 最能解釋資料。對於仍有多條合法 parent path 的候選，再以 HP-specific read-AF 作方向相容性檢查：在 copy-neutral、無 LOH 且 mutation multiplicity 相同時，祖先突變的 carrier fraction 理論上不應低於後代。因此這個步驟可以排除明顯不合理的方向，或標出一條「條件式較偏好」的 path；但不能單獨證明 direct parent、hidden clone 或唯一真實演化樹。無法由資料區分時，我們保留多解或 abstain，而不強迫選一棵。

一句話版本：

> **Likelihood 決定哪些 mutation states 能解釋 reads；HP-specific read-AF 檢查候選路徑的祖先盛行率是否合理。AF 是方向約束，不是獨立真值。**

## 10. 本輪決策

| 選項 | 判定 | 理由 |
|---|---|---|
| 以 HP-specific read-AF 取代 likelihood | 不採用 | 丟失 joint pattern、BQ 與 missingness 資訊 |
| likelihood 後再加同 BAM 的 VAF score | 不採用 | 重複使用同一 allele-count evidence |
| likelihood 選 \(V\)，AF/CCF 檢查 \(E\) | 採用 | 分工清楚、可解釋、保留 non-identifiability |
| raw AF 微差直接決定唯一 topology | 不採用 | depth/CN/LOH/coverage confound，且無 uncertainty |
| 校正後 CCF＋bootstrap 作條件式 edge 偏好 | 建議驗證 | 可排除不相容方向，但仍需外部 truth 校準 |
| 無法區分時保留多樹／abstain | 採用 | 符合資料辨識上限 |

## 11. 驗證計畫

1. 建立 synthetic same-\(V\)/different-\(E\) fixtures  
   → 驗證：likelihood 保持 tie，AF layer 只輸出 conditional preference。
2. 建立 CN-neutral known-mixture simulation  
   → 驗證：ancestry constraint 的方向錯誤率、coverage sensitivity 與 calibration curve。
3. 加入 CN gain、LOH、multiplicity、mapping bias 與 partial coverage  
   → 驗證：eligibility gate 能 abstain，不製造假唯一。
4. 在 HCC1395 可對應的 CN-neutral strict units 做 bounded pilot  
   → 驗證：報告 `evaluable / preferred / tied / discordant / abstain` 完整分母。
5. 與 single-cell、multi-region/timepoint、synthetic truth 或 known mixture 比較  
   → 驗證：被偏好的 path 在外部 truth 的 accuracy 高於 likelihood-only tie 的任選基準。
6. 凍結 Python oracle 後才移植 C++  
   → 驗證：candidate digest、LL、KKT gap、AF eligibility、edge status 與 topology claim 全欄位 parity。

### 11.1 建議預註冊的成功 gate

- synthetic truth 必須有至少 99% 被完整 candidate family 包含；
- 對 matched call-rate units，exact-path accuracy 相對 likelihood-only／raw-AF baseline 至少提升 5 percentage points，且 block-bootstrap 95% CI 下限大於 0；
- false-unique rate 不高於 5%；
- 理論不可辨 same-\(V\)/different-\(E\) cases 被強迫唯一的比例不高於 5%；
- winner 對 CN source、purity 與 multiplicity plausible scenarios 的穩定率至少 90%；
- HCC1395 與 HCC1395_DORADO 的 split-half／depth-matched concordance 改善，且 CCF conflict 不增加。

這些是預註冊建議，不是已通過的結果；threshold 必須在 production 前以 simulation／known mixture 校準。

## 12. HCC1395 現況：可做 pilot，尚不可升格 production

### 可用資產

- SAVANA allele-specific CN：
  `/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv`
  - 737 行（含表頭）；
  - SHA256：`5b1d525eb398d897eee0d935b5ab7a7b7feb4ab13f207e8b414ce33855ddb4d7`；
  - primary purity/ploidy solution：0.96／2.79，第二解為 0.76／1.84。
- Wakhan HP-specific CN 可作 sensitivity，但 Wakhan HP1/HP2 不可直接假定等於 LongPhase-S HP1/HP2；必須在每個 exact PS 內以 germline SNP 建立 orientation authority。
- 既有 PyClone-VI 可提供 bulk cellular-prevalence 比較，但不是 HP-specific，且使用 legacy regions／共享 CN context，不是 strict path truth。
- current strict exact-PS read-linked regions 可作新 pilot 的 unit universe。

### 現有反證與缺口

- 歷史 HCC1395 paired-source VAF heuristic 在 3,860 個可評估 legacy regions 中，有 1,234 個 candidate-space conflicts：

\[
\frac{1,234}{3,860}=31.97\%.
\]

這證明 VAF 增加「決斷性」時也可能增加來源敏感性，不能直接視為 accuracy。

- 尚缺：
  - HCC1395_DORADO source-specific CN／purity；
  - mutation multiplicity posterior；
  - Wakhan 與 LongPhase-S 的 exact-PS HP orientation authority；
  - single-cell、colony、multi-region/timepoint 或 synthetic exact-edge truth；
  - current strict 7-dataset likelihood／directed-topology receipts。

因此目前結論是：

> **GO for bounded HCC1395 CN-neutral pilot；NO-GO for raw-AF production winner rule。**

## 13. 證據與來源

### Repository evidence

- `InterSubMod/research/20260712_read_weighted_hypercube_tree/pre-decision-audit.md`
- `InterSubMod/research/20260713_hp_readstate_clone_identifiability/20260713_HP分群Read共現與區域Clone可辨識性結論_01.md`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/read_af_tree_ordering_multisample.py`
- `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py`

### Primary literature

- Roth et al., PyClone, *Nature Methods* 2014: https://doi.org/10.1038/nmeth.2883
- Deshwar et al., PhyloWGS, *Genome Biology* 2015: https://doi.org/10.1186/s13059-015-0602-8
- Popic et al., LICHeE, *Genome Biology* 2015: https://doi.org/10.1186/s13059-015-0647-8
- Niknafs et al., SCHISM, *PLOS Computational Biology* 2015: https://doi.org/10.1371/journal.pcbi.1004416
- Qi et al., bulk phylogeny non-uniqueness, *Algorithms for Molecular Biology* 2019: https://doi.org/10.1186/s13015-019-0155-6
- Wintersinger et al., Pairtree, *Blood Cancer Discovery* 2022: https://doi.org/10.1158/2643-3230.BCD-21-0092

## 14. Claim ceiling

目前可說：

> 在完整候選 family 內，joint read-pattern likelihood 可評估不同 mutation-state vertex sets；HP-specific read-AF／校正後 CCF 可在嚴格前提下檢查 ancestry order 的相容性，並標出條件式 preferred path。

目前不可說：

> 單憑 bulk HP-specific read-AF 已證明真實直接 parent-child、hidden clone、唯一 exact topology、真實 clone 數或完整腫瘤演化史。
