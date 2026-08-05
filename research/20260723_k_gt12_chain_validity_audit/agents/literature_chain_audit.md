<!--
建立時間: 2026-07-23 18:20
目標: 以外部原始文獻與目前 InterSubMod 實作，審查 k>12 分區、A-B/B-C 串聯、partial observation 與局部系統發生重建的合理性
處理範圍: 方法定義、文獻支持邊界、目前 HCC1395 exact PS×HP pilot；不包含新樣本實驗、臨床 clone 驗證或全 cohort 結論
關聯檔案:
  - InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/pre-decision-audit.md
  - InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/implementation-notes.md
  - InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/20260723_HCC1395_exactPS拓撲重建觀察_01.md
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/read_support_partition.py
服務目標: G1、G5
任務類型: B（限定於方法定義與文獻的完整驗證；不是全樣本生物驗證）
狀態: literature-audited / local-method-audited / biological-validation-pending
-->

# `k > 12` 與 A–B／B–C 串聯合理性：文獻與方法邊界審查

> **TL;DR：保留完整 source component，同時只把 exact solver 的輸入切成 `k≤12` blocks；不要把 `k>12` 直接叫「高密度區域」。A–B 與 B–C 足以建立同一個 read-linked analysis component，但不等於 A、B、C 位於同一條分子、同一個細胞或同一個 clone。沒有跨 block bridge reconstruction 時，結果只能稱為 local candidate topology。**

敘述框架：SCQA（問題情境 → 方法衝突 → 判定問題 → 建議答案）加上 Claim–Evidence–Boundary（主張 → 證據 → 可宣稱邊界）。

## 1. 審查結論

### 1.1 三個直接答案

| 問題 | 判定 | 精確說法 |
|---|---|---|
| `k > 12` 要切開，還是視為特殊密集區域？ | **兩者不是互斥選項。** | 原始 read-linked component 應完整保留並標記 `k>12`；為了 bounded exact computation，再產生 `k≤12` 的 solver blocks。只有另算位點/距離密度後，才能再標記 `high-density`。 |
| A–B 與 B–C 是否可串成一個 block／region？ | **可作共同分析範圍，但不可直接作共同 clone。** | graph connectedness 是合理的 analysis-scope 定義；A–B 與 B–C 是兩筆直接分子證據，A–C 只是間接連通，並非同一條 read 的直接觀察。 |
| 分塊後是否仍可宣稱整個 component 的演化樹？ | **目前不可以。** | 現行非重疊分塊會把跨切點 constraint 分成 `cut` 或 `unavoidable_span_gt_max_block_size`；若未建構跨 block bridge likelihood／separator state，只能報每個 block 的 local candidate topology。 |

### 1.2 最推薦的雙層資料模型

```text
source component（永遠保留）
  exact PS × HP × read-linked connected component，k 可大於 12
  ├─ source_component_id
  ├─ 完整 sites / molecules / pair or hyperedge evidence
  ├─ density 指標（另算，非由 k>12 代替）
  └─ solver partition（計算用）
       ├─ block_1, k≤12 → local candidate topology
       ├─ block_2, k≤12 → local candidate topology
       └─ cut / unavoidable constraints → 跨 block 未解證據
```

因此，最安全的政策不是「保留或切開」二選一，而是：

1. **資料層保留完整 component**：避免失去原始關係與可追溯性。
2. **計算層切成 bounded blocks**：讓 exact solver 可完成。
3. **解讀層區分 local 與 global**：未完成跨 block 重建前，不宣稱 global topology。

## 2. 必須先固定的定義

### 2.1 建議術語

| 術語 | 本研究建議定義 | 不應混寫成 |
|---|---|---|
| site／character | 一個納入 binary mutation-state model 的 sSNV 位點 | 一個 clone |
| `k` | 某一 source component 或 solver block 內的 site／character 數 | read 數、clone 數、密度 |
| direct molecule link | 同一條通過 QC 的 read 在至少兩個 site 上具有 fixed `R/A` call | 同一個細胞的永久 barcode |
| read-linked component | 以 site 為 vertex、同一 read 的多位點 fixed call 為 edge／hyperedge 所形成的 connected component | 一條 DNA 分子、單一 haplotype、單一 clone |
| source component | 分區前完整的 exact `PS×HP×read-linked component` | 可直接求解的 block |
| bounded block | 為 exact solver 產生、連續且不重疊、`k≤12` 的計算單位 | 無損的完整 region |
| local candidate topology | 單一 bounded block 內，在指定模型與證據下相容或最小的候選拓撲 | 已確認的病人 clone tree |
| global topology | 同一 source component 所有 blocks 經跨邊界 evidence 共同約束後的拓撲 | 各 local tree 的任意拼接 |
| partial observation | read 未對某位點提供 fixed `REF/ALT` call；應保留為 missing／unknown | `REF=0` |
| density | site 數相對於 genomic span 的密集程度 | `k>12` 本身 |

### 2.2 `k > 12` 不等於「高密度」

至少要分開三個量：

\[
k = \text{component 內納入的 sSNV 數}
\]

\[
\text{span}_{kb} = \frac{\max(pos)-\min(pos)+1}{1000}
\]

\[
\text{site density} =
\begin{cases}
k/\text{span}_{kb}, & \text{若以每 kb 位點數表示}\\
(k-1)/\text{span}_{kb}, & \text{若以相鄰間隔事件率表示}
\end{cases}
\]

例如，13 個位點分布於 5 kb 與 13 個位點分布於 1 Mb，都有 `k=13`，但生物與技術上的密度完全不同。建議使用：

- `oversized_for_exact_solver`：`k>12`；
- `high_site_density`：另有事先定義且經資料分布校準的 density threshold；
- 兩個 flag 可以同時存在，但不可互相代替。

### 2.3 目前 `R/A/O/D/S/L/X` 的正確語意

目前 extractor 保留 `R/A/O/D/S/L/X` 原因碼；segmentation linkage 只使用 fixed `R/A`：

| code | 本地程式語意 | 是否可建立 linkage |
|---|---|---|
| `R` | base 與 REF 相同且通過品質門檻 | 是 |
| `A` | base 與 ALT 相同且通過品質門檻 | 是 |
| `O` | 觀察到 REF/ALT 以外 base | 否 |
| `D` | deletion CIGAR operation | 否 |
| `S` | reference skip（CIGAR `N`）原因碼；在其他報表中 `S` 又可能代表候選 sSNV 總數 | 否 |
| `L` | base quality 太低或缺品質 | 否 |
| `X` | alignment span 內仍無法取得可判定 call | 否 |

**重要校正**：未知狀態不應說成 `S`。在本流程的 partial pattern 敘述中建議使用 `X` 或數學記號 `?`；`S` 已有其他語意，容易造成口試時混淆。`X/?` 也不可直接補成 `R=0`，而應在後續 likelihood 或 completion model 中 marginalize／枚舉相容狀態。

## 3. A–B 與 B–C 串聯：哪一層成立？

### 3.1 成立：connected component 作為共同分析範圍

定義 site graph：

\[
G=(V,E),\quad
V=\{\text{sSNV sites}\},\quad
(u,v)\in E \iff \exists\text{ 一條通過 QC 的 read 對 }u,v\text{ 都有 fixed R/A call}
\]

若有：

- read 1 支持 A–B；
- read 2 支持 B–C；

則 A、B、C 在 graph theory 上位於同一個 connected component。把三者放進同一個「待共同推論的 source component」是合理的，而且與 haplotype assembly 以 read overlap 延伸 phase block 的思想一致。

### 3.2 不成立：把 transitive connectivity 當成同一分子

`same-read(A,B)` 與 `same-read(B,C)` **不推出** `same-read(A,C)`。前兩個 relation 由不同 reads 支持時，A 與 C 從未被同一條 molecule 直接觀察。

因此，報告中必須區分：

- **直接證據**：A–B、B–C；
- **間接 graph relation**：A 與 C 同 component；
- **尚未觀察**：A–C 是否位於同一 molecule。

### 3.3 不成立：把 connected component 當成同一 clone

在 infinite-sites／perfect-phylogeny 直觀下，同一條 molecule 同時看到兩個 ALT mutation，可約束兩 mutation 位於某條共同 root-to-leaf path；但兩個 pair constraint 的連通性仍不保證三者同一路徑。

反例：

```text
      B
     / \
    A   C
```

- 某個 B→A 後代細胞可提供 `A=ALT, B=ALT` read；
- 另一個 B→C 後代細胞可提供 `B=ALT, C=ALT` read；
- A 與 C 是 sister branches，未必存在同時帶 A、B、C 的 clone。

所以：

\[
AB\text{ comparable} \land BC\text{ comparable}
\not\Rightarrow AC\text{ comparable}
\]

這也是為什麼 A–B／B–C 可以建立共同推論問題，卻不能直接宣稱「三者是同一個 clone 的最小組合」。

### 3.4 在 diploid germline phasing 與 bulk tumor 中，推論強度不同

| 場景 | A–B＋B–C 可推出什麼 | 額外前提 |
|---|---|---|
| diploid germline haplotype phasing | 若 allele orientation 一致，常可間接推定 A–C phase，形成較長 phase block | 恰有兩條主要 haplotypes、genotype 與 read-error model、無衝突／chimera |
| bulk tumor subclone inference | 只足以形成 joint constraint component；可能是 trunk mutation B 連到兩個 sister subclones | clone mixture、purity、CN state、sequencing error、VAF／cell fraction、phylogeny model |
| single cell | 同一細胞 barcode 內的多位點可提供較直接的 cellular genotype | dropout、doublet、allelic amplification bias 仍需建模 |

本研究是 bulk tumor long-read 資料，因此不可直接借用「兩條 germline haplotype」的 transitivity，當成 clone identity。

## 4. 現行 `k≤12` 分區實際做了什麼

### 4.1 已確認的工程行為

目前 `read_support_partition.py`：

1. 將排序後的 genomic sites 視為 vertices；
2. 將一條 aggregated read pattern constraint 視為 hyperedge；
3. 產生連續、非重疊、完整覆蓋 sites、每塊 `k≤max_block_size` 的 blocks；
4. 依序最佳化：
   - retained hyperedge weight 最大；
   - retained pattern count 最大；
   - block 數最少；
   - cut genomic gap 總和最大；
   - cut index 字典序固定，以保持 deterministic；
5. 一條 constraint 只有在所有 sites 都落在同一 block 時才算 retained；
6. 跨 cut constraint 被標成 `cut`；其 site index span 本來就大於 block 上限時，標成 `unavoidable_span_gt_max_block_size`。

這個演算法是**證據保留最佳化的 bounded partition**，不是跨 block 的 global phylogeny solver。

### 4.2 現行做法合理的條件

若研究問題明確限定為：

> 「在每個 bounded block 內，列舉與 retained molecule patterns 相容的 local mutation-state candidate topologies。」

則目前做法理論上合理，且工程上可稽核。它服務的是局部推論、計算可完成性與 evidence accounting。

### 4.3 現行做法不足的情況

若要回答：

> 「整個原始 `k>12` source component 的單一 global topology 是什麼？」

目前非重疊分區仍不足，因為：

- `cut`／`unavoidable` constraints 不再同時進入任何 local solver；
- 各 block 的 hidden state 命名沒有自動的 shared-cell identity；
- local minimum 的組合不一定等於 global minimum；
- 每塊各自存在多個 equally optimal topology 時，拼接會產生組合爆炸；
- CNV、LOH、purity 與 bulk mixture 可讓 binary perfect-phylogeny 前提失效。

### 4.4 `k=12` 的正確定位

`k=12` 目前是**工程與 exact enumeration 的上限**，不是領域公認的生物 cutoff：

| `k` | binary state vertices `2^k` |
|---:|---:|
| 10 | 1,024 |
| 11 | 2,048 |
| 12 | 4,096 |
| 13 | 8,192 |
| 14 | 16,384 |
| 20 | 1,048,576 |

這只顯示 state space 的指數增長；候選 Steiner subgraph／tree search 還可能比單純列舉 vertices 更昂貴。外部文獻支持「需 bounded exact computation」，但**沒有替本研究證明 12 是最佳值**。`12` 必須由本地 benchmark 支持，例如 `k=10..14` 的：

- wall time；
- peak RSS；
- complete-enumeration rate；
- capped／timeout rate；
- topology ambiguity；
- retained evidence 與 cut evidence；
- 結果對 cutoff 的 sensitivity。

## 5. 外部原始文獻逐篇查核

### 5.1 HapCUT2：read-fragment graph 與 phase block

來源：Edge、Bafna、Bansal，*Genome Research* 2017，Methods「Haplotype likelihood」「Likelihood-based HapCUT2 algorithm」「Complexity of HapCUT2」。  
[原始論文（PMC）](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411775/)；[DOI](https://doi.org/10.1101/gr.213462.116)

**支持本方法的部分**

- 一條 read／fragment 覆蓋多個 heterozygous variants 時，可提供 partial haplotype information。
- variant 為 graph nodes，fragment 支持的 variant pair 可形成 edges。
- 可按 read-haplotype graph 的 connected component 分開求解，因此 A–B、B–C 用於形成同一 analysis component 是標準思想。
- 未覆蓋位置以 missing 表示，likelihood 只使用 read 實際覆蓋的 sites。

**不支持的部分**

- connected component 並不是一條實體 molecule。
- 連通性本身不證明 A–C 被同一 read 直接觀察。
- phase block 不等於 tumor clone／subclone lineage。

**與本研究的差異**

- HapCUT2 的核心場景是 diploid haplotype assembly；本研究是 bulk tumor、多 latent cellular states。
- HapCUT2 使用全域 likelihood 與 sequencing-error model，不是只用 thresholded connectivity 後直接指定 clone。

### 5.2 LongPhase：相鄰 variant graph 與有證據的 block stitching

來源：Lin 等，*Bioinformatics* 2022，Methods §2.4.1、§2.4.2、Figure 2。  
[原始論文](https://academic.oup.com/bioinformatics/article/38/7/1816/6519151)；[DOI](https://doi.org/10.1093/bioinformatics/btac058)

**支持本方法的部分**

- heterozygous variants 可形成依 genomic order 排列的 graph；
- read-supported allele connections 可把局部關係串成較長 haplotypes；
- initial blocks 可以再由跨 block long-read evidence 重新 phasing 與 concatenation。

**不支持的部分**

- LongPhase 不是僅因兩 blocks 相鄰就拼接；它需要實際 spanning-read evidence。
- graph component 不是 clone；
- 跨區拼接仍有 chimera、alignment artifact 與距離上限等防護。

**與本研究的差異**

- LongPhase 解的是兩條 germline haplotypes 的 path optimization。
- 本研究若要從 local mutation-state trees 形成 global tumor topology，仍需要另外定義跨 block hidden-state correspondence、bulk mixture 與 error model。

### 5.3 PairClone：same-fragment mutation pairs 可增加 subclone 訊息

來源：Zhou 等，*Journal of the Royal Statistical Society Series C* 2019，§1.2、§2.1、Figure 1。  
[原始論文](https://academic.oup.com/jrsssc/article/68/3/705/7058355)；[DOI](https://doi.org/10.1111/rssc.12328)

**支持本方法的部分**

- proximal mutation pairs 若被同一 paired-end fragment phase，可比兩個獨立 marginal VAF 提供更多 subclone 訊息。
- 觀察型態明確包含 `00,01,10,11,-0,-1,0-,1-`；`-` 是 missing，不應當成 REF。
- partial reads 可納入 multinomial sampling model，而不是直接丟掉。

**不支持的部分**

- 不支持把由不同 molecules 提供的 A–B 與 B–C，自動視為一筆三位點 observed genotype。
- 不支持「connected component = clone」。

**與本研究的差異**

- PairClone 共同估計 subclone 數、genotype 與 proportions，並含 background error component。
- 本研究目前的 local exact/minimal topology 並沒有等價的完整 posterior mixture model。

### 5.4 TreeClone：mutation-pair data 必須放入 mixture 與 tree model

來源：Zhou 等，*Annals of Applied Statistics* 2019，Introduction、§2.2、§2.3。  
[原始論文 DOI](https://doi.org/10.1214/18-AOAS1224)；[作者專案頁](https://ccte.uchicago.edu/treeclone/)

**支持本方法的部分**

- mutation-pair data 與 partial observations 可用於 subclone phylogeny inference。
- 以 Bayesian latent-feature model 共同推論 subclone 數、genotype matrix、sample proportions 與 tree，可比獨立 VAF 更直接使用 phase information。

**不支持的部分**

- pair graph 的 connectedness 本身不是一棵 clone tree；
- 不支持忽略 copy number、sampling error、cell fraction 後直接將最小 state graph 當作真實 clones。

**與本研究的差異**

- TreeClone 使用明確 tree prior 與 mixture likelihood，且其主要設定限制於 copy-number-neutral region。
- 本研究的 candidate Steiner topology 比較接近「相容解空間與 minimal hidden states」，而非 posterior-confirmed cellular lineage。

### 5.5 Bulk phylogenetic deconvolution 的非唯一性：long-read pair constraints 仍不足以唯一決定樹

來源：Pradhan、El-Kebir，*Algorithms for Molecular Biology* 2019，「Additional constraints on the solution space」、Figures 2/5、Conclusions。  
[原始論文（PMC）](https://pmc.ncbi.nlm.nih.gov/articles/PMC6719395/)；[DOI](https://doi.org/10.1186/s13015-019-0156-x)

**支持本方法的部分**

- 同一條 long read 上的 mutation pair 可提供共同 root-to-leaf path 的額外 constraint。
- single-cell 與 long-read evidence 都可縮小 bulk deconvolution 的解空間。

**不支持的部分**

- long-read pair constraints 只能減少、不能保證消除 non-uniqueness；
- 真實 read error 必須另建模；
- mutation／解越多時，可行 phylogenies 仍可能很多。

**對 A–B／B–C 的關鍵意義**

- A 與 B comparable、B 與 C comparable，不強迫 A 與 C comparable；
- B 可以是 A、C 的共同祖先，A、C 位於分支兩側；
- 因此 connected component 是合理 problem scope，但不是 unique root-to-leaf clone path。

### 5.6 Incomplete Directed Perfect Phylogeny：unknown 是待完成狀態

來源：Pe’er、Pupko、Shamir、Sharan，*SIAM Journal on Computing* 2004。  
[原始論文／DOI](https://doi.org/10.1137/S0097539702406510)

**支持本方法的部分**

- binary characters 中可有 unknown entries；
- 問題是尋找一個 missing-value completion，使完成後資料承認 directed perfect phylogeny；
- 一般可能有多個 completions／solutions。

**不支持的部分**

- 不處理 bulk sample proportions、ONT read error、copy number 或 methylation；
- 不支持把 `X` 靜默補成 REF；
- 不支持「某一個最小 completion 就是唯一真實 clone」。

**與本研究的差異**

- 這篇文獻提供 partial state 的形式化基礎；實際 bulk ONT 應再加 likelihood、error 與 CN-aware boundary。

### 5.7 Directed hypercube／Steiner tree：minimal hidden states 是數學候選，不是已觀察 clone

來源：

- Mahapatra、Narayanan、Narayanaswamy，*Acta Informatica* 2025。  
  [原始論文／DOI](https://doi.org/10.1007/s00236-024-00474-8)
- Blelloch 等，ICALP 2006。  
  [原始論文／DOI](https://doi.org/10.1007/11786986_58)

**支持本方法的部分**

- `k` 個 binary characters 對應 `2^k` 個 hypercube vertices；
- observed binary states 可作 terminals；
- Steiner arborescence／tree 可加入未觀察的 internal states；
- binary maximum-parsimony phylogeny 與 hypercube Steiner formulation 有直接數學關係；
- exact computation 具指數性，因此設定 bounded `k` 有工程合理性。

**不支持的部分**

- 文獻未支持特定 cutoff `k=12`；
- 把一個 global instance 分成互不重疊 blocks，不保證保留 global optimum；
- Steiner node 是最佳化所需的 latent state，不自動等於真實存在且被取樣的 clone。

## 6. 「最小圖」可以怎麼說才精確

原說法：

> 「三者不能說一定是同一個 clone，但可以用最小圖將它們判定為同一種最小組合。」

建議修正為：

> 「A–B 與 B–C 的 molecule-supported edges 使三個位點落在同一個 read-linked analysis component。我們接著在指定的 binary-state、partial-observation 與 parsimony 假設下，尋找能解釋 retained read patterns 的最小候選拓撲。這個結果代表一個最簡相容模型，不等於已證實三個突變存在於同一 clone，也不保證候選拓撲唯一。」

再更短的口試版本：

> 「串聯先定義共同分析範圍，最小圖再提供最簡相容解；兩者都不是 clone 身分的直接觀察。」

## 7. 可宣稱與不可宣稱

### 7.1 目前可宣稱

1. 同一條通過 QC 的 long read 對兩個以上位點提供 fixed `R/A` calls，可形成 direct molecule-supported constraint。
2. A–B、B–C 可以透過 shared vertex B 形成同一 read-linked component。
3. component connectedness 可用來決定哪些 sites 進入同一個 joint inference problem。
4. `k>12` source component 可用 read-support objective 分為 `k≤12` blocks，以進行 bounded exact local reconstruction。
5. 每一 block 可輸出 local candidate mutation-state topology 與 ambiguity／completeness 狀態。
6. cut 與 unavoidable constraints 可量化分區造成的 evidence loss。

### 7.2 目前不可宣稱

1. A、B、C 被同一條 molecule 直接跨越。
2. A、B、C 一定屬於同一個 clone。
3. 同一 read-linked component 就是一個 haplotype 或 clone。
4. 每個 hidden Steiner state 都是真實存在的 cellular clone。
5. 各 `k≤12` local trees 可無證據直接拼成原 component 的 global tree。
6. `k>12` 等於生物學高密度／hotspot。
7. `k=12` 是領域標準或最佳生物 cutoff。
8. 在未處理 CNV／LOH、purity、read error 與 non-uniqueness 前，把 local candidate topology 稱為臨床 clone evolution tree。

## 8. `k>12` 決策表

| 目標 | 建議處理 | 輸出命名 | 可否宣稱 global |
|---|---|---|---|
| 只需 exact local topology | 保留 source component；依 retained-read objective 切 `k≤12` | `local_candidate_topology` | 否 |
| 研究高密度／hotspot | 另算 site density、span、mappability、CN／repeat context | `high_site_density`（需門檻） | 與 topology 無直接等價 |
| 要整個 component 的 global topology | 使用 overlapping blocks、boundary separator states、cross-block bridge likelihood 或 approximate whole-component solver | `global_candidate_topology` | 只有在 bridge validation 通過後 |
| `k` 極大且 evidence 稀疏 | 先做 graph decomposition／articulation separator，再保留所有跨 separator constraints | `decomposed_global_inference` | 需證明 decomposition 的 losslessness 或量化 approximation |
| 高衝突／高 ambiguity | 不強迫輸出 winner，回報 `ABSTAIN/UNRESOLVED` | `unresolved_component` | 否 |

## 9. 下一步驗證計畫

### P0：先修正資料與簡報定義

1. 將 `k>12` 改稱 `oversized_for_exact_solver`。  
   → 驗證：所有 report/schema 不再把 `k>12` 等同 `dense`。
2. 固定 `source_component_id` 與 `solver_block_id` 兩層 identity。  
   → 驗證：每個 block 可回溯原 component；所有 `cut/unavoidable` constraints 可回溯原 molecule patterns。
3. 將 unknown 統一寫成 `X/?`，避免以 `S` 表示。  
   → 驗證：簡報、schema、speaker notes 對 `R/A/O/D/S/L/X` 無衝突定義。

### P0：建立 A–B／B–C 反例測試

至少建立四個 synthetic fixtures：

1. 同一 read 直接跨 A–B–C；
2. 不同 reads 支持 A–B 與 B–C，但 A–C 無 direct read；
3. trunk B、sister branches A/C；
4. A–B 與 B–C 的 allele orientation 互相衝突。

→ 驗證：系統只對情境 1 標 direct three-site molecule evidence；情境 2/3 只能標 connected；情境 4 應 conflict／ABSTAIN，而非強制拼接。

### P1：證明 `k=12` 是可辯護的工程 cutoff

在相同輸入與 budget 下測 `k=10,11,12,13,14`：

1. wall time；
2. peak RSS；
3. exact complete rate；
4. capped／timeout rate；
5. candidate tree count；
6. local topology stability；
7. retained／cut evidence ratio。

→ 驗證：能指出 `k=12` 相對 11/13 的明確 Pareto trade-off，而不是只寫「因為 2^12」。

### P1：量化分區對結果的影響

對可在更高 budget 下完整求解的小型 truth fixtures：

1. full-component exact result；
2. non-overlap `k≤12` local result；
3. overlap／bridge-aware reconstruction；
4. 比較 terminal coverage、hidden-state count、edge set 與 topology ambiguity。

→ 驗證：能量化 local partition 對 global answer 的偏差；若無法恢復，明確維持 local-only claim。

### P2：若要升級成 clone claim

至少補：

1. VAF／cellular prevalence 與 sample mixture constraint；
2. purity 與 CNV／LOH-aware modeling；
3. read error、chimera、mapping artifact model；
4. 多區域或多時間點一致性；
5. single-cell、orthogonal assay 或已知 mixture truth 的外部驗證。

→ 驗證：hidden states 不只滿足 parsimony，也能與 abundance、copy number 及外部 cellular evidence 相容。

## 10. 建議放進論文的方法敘述

> 本研究先在同一 dataset、chromosome、exact phase set 與 haplotype tag 範圍內，將 sSNV 視為 vertices，並以同一條通過品質控制之 long read 上至少兩個 fixed REF/ALT observations 建立 molecule-supported hyperedges。由這些 hyperedges 所形成的 connected component 定義共同分析範圍；connectedness 僅表示位點可由 read evidence chain 間接連結，不代表所有位點由同一分子、同一細胞或同一 clone 共同觀察。當 source component 的 character 數 `k>12` 時，本研究保留完整 component 與所有 constraint identity，另以 retained-read-weight objective 分割為連續、非重疊且 `k≤12` 的 bounded blocks，供 exact local topology reconstruction 使用。跨切點 evidence 會被明確標記為 cut 或 unavoidable；在尚未完成 cross-block bridge reconstruction 前，輸出僅解讀為 local candidate mutation-state topology，不宣稱原 component 的 global clone lineage。

## 11. 建議口試回答

教授若問「A–B 和 B–C 可以串起來嗎？」：

> 「可以串成同一個 read-linked analysis component，因為 B 是共同節點；但要分清楚，A–B 和 B–C 是直接 read 證據，A–C 只有間接連通。對 diploid haplotype，額外假設成立時可以做 transitive phasing；但對 bulk tumor，B 也可能是 trunk mutation，而 A、C 是不同分支，所以不能僅由連通性說三者是同一 clone。」

教授若問「為什麼 `k>12` 要切？」：

> 「12 是目前 exact solver 的工程上限，不是生物 cutoff，也不是密度定義。我們保留原始完整 component，只把計算輸入切成 `k≤12` blocks，並記錄被切斷的 read constraints。因此現在可宣稱的是局部候選拓撲；若要整區 global tree，還需要跨 block bridge reconstruction。」

教授若問「最小圖就是 clone tree 嗎？」：

> 「不是。最小圖是在既定 binary-state、partial observation 和 parsimony 假設下，解釋觀察資料所需的最簡候選模型。hidden state 可能是尚未取樣的祖先狀態，也可能只是模型補點；必須再結合 VAF、copy number、purity 與外部細胞證據，才能往 clone lineage 解讀。」

## 12. 證據與再現紀錄

### 12.1 本地輸入

- `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/pre-decision-audit.md`
- `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/implementation-notes.md`
- `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/20260723_HCC1395_exactPS拓撲重建觀察_01.md`
- `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py`
- `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/read_support_partition.py`
- `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/extract_lossless_read_linkage_collapsing.py`

### 12.2 查核命令

```bash
rg -n "k.?12|exact|PS|component|local|global|clone" \
  InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/pre-decision-audit.md

rg -n "max_block_size|unavoidable_span|retained|hyperedge|objective|contiguous|cut|blocks|span" \
  InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/read_support_partition.py

rg -n "R/A/O/D/S/L/X|NON_LINKING_CALLS|FIXED_CALLS" \
  InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot \
  InterSubMod/research/20260718_k_gt8_read_supported_segmentation
```

### 12.3 證據層級

| 層級 | 本報告用途 |
|---|---|
| Primary literature | 判定領域方法「支持什麼／不支持什麼」 |
| Local code and audit | 判定目前實作實際如何分區、哪些 constraint 會失去 block-local eligibility |
| Logical inference | A–B／B–C 的 transitivity 反例、local/global claim boundary |
| 尚缺 empirical validation | `k=12` 最佳性、跨 block global recovery、hidden state 的 clone identity |

## 13. 最終 verdict

**方法核心可保留，但宣稱必須收斂：**

- A–B／B–C chaining：**GO**，只用於建立 joint analysis component。
- `k>12` bounded partition：**GO**，只用於 exact local computation，且完整 source component 與 cut evidence 必須保留。
- 把 `k>12` 稱高密度：**NO-GO**，除非另有 density 定義與門檻。
- 把 connected component 稱同一 clone：**NO-GO**。
- 把各 block local minimal trees 直接拼成 global clone tree：**NO-GO**，直到 cross-block bridge reconstruction 與 sensitivity validation 完成。
- 將結果稱為「local candidate mutation-state topology」：**GO**，最符合目前證據。

