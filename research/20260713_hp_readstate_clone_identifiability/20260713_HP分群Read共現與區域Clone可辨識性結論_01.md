<!--
建立時間: 2026-07-13T05:45:02+08:00
目標: 凍結目前對「germline HP 分群後，以同-read sSNV 聯合狀態判讀區域 clone/subclone」的共識、公式、適用條件與 claim ceiling
處理範圍: 兩位點至多位點 sSNV；bulk ONT；LongPhase-S germline HP family；full/partial read states；allele-specific CN；regional mutation-state tree
任務類型: A — 聚焦式方法學整理（PARTIAL；不是全癌症演化驗證）
服務目標: G4（可重現口徑）/ G5（外部可驗證）
狀態: current decision record / scientific claim PARTIAL
關聯檔案:
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md
  - InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py
-->

# HP 分群 Read 共現與區域 Clone 可辨識性：目前結論

> **PARTIAL METHOD INTERPRETATION / NOT BIOLOGICAL CLONE TRUTH**
>
> 用 Assertion–Evidence：**germline HP 分群後的同-read sSNV 聯合狀態，能直接觀察／支持 homolog-specific regional mutation states；若 same PS 或跨 PS orientation 已獨立驗證、每個 state 均由多條獨立 full-span molecules 支持、錯誤與取樣偏差受控、採用§5.1的標準穩定cell-genotype模型，且完整分析窗中所有 contributing cells 的相關 HP 均為 allele-specific CN=1、無 subclonal CNA，互異的同-HP states 才可進一步推論來自不同來源細胞。Mutation-set 的包含或互不包含關係分別提供 descendant/subclone-compatible 或 sister/branch-compatible evidence，但不能單獨確認 genome-wide biological clone、extant clone 數或真實祖先樹。**

## 0. 本輪決議（先讀這一節）

1. **方法方向合理**：sSNV same-read joint states 比 VAF-only 多出位點間聯合資訊；先依 germline HP 分開，可避免把兩條 homolog 的差異誤畫成細胞分枝。
2. **可直接確認的對象是 molecule/haplotype state**：例如同一 HP-family 的 full reads 直接支持 `10`、`11` 兩種 regional mutation states。
3. **「不同細胞」是有條件的強推論**：只有在 same PS／可靠跨-block orientation、每個state均有多條獨立full-span molecules、錯誤與取樣偏差受控、採用§5.1的標準穩定cell-genotype模型，且完整分析窗中所有contributing cells的相關HP均為allele-specific CN=1、無subclonal CNA時，同一HP的互異exact states才不可能出自同一細胞的單一homolog copy。
4. **Nested 不等於祖先已確認**：`10` 與 `11` 在 infinite-sites、無 loss、無 recurrence 下與 `10→11` 相容；這是必要模型條件下的方向性，不是 lineage truth。
5. **Incomparable 不等於姐妹 clone 已確認**：`10` 與 `01` 在相同假設下與分枝／互斥細胞群相容；CN gain、recurrence、loss 或 phase error 仍可產生替代解釋。
6. **現行 InterSubMod 輸出仍是 regional mutation-state forest**：每個 `region × HP-family` 獨立建樹；跨 HP 不做 cell pairing，也不建立 ancestry edge。
7. **正式用途**：找出 regional clone/subclone candidates、計算可辨識性與 cell-overlap bounds、排序後續 allele-specific CN、多樣本、single-cell 或 colony validation；不可直接報成 clone/subclone 發生率。

## 1. 為何這個研究問題重要

VAF-only 只提供各 mutation 的邊際頻率：

\[
P(A),\qquad P(B)
\]

但同一 molecule/read 可提供兩位點聯合狀態：

\[
P(00),\quad P(10),\quad P(01),\quad P(11)
\]

因此 same-read evidence 能回答「兩個 sSNV 是否在同一 homolog molecule 上共現、互斥或呈 nested pattern」，這是 VAF-only 無法保留的資訊。再將 reads 依 germline HP family 分開，可把 **allelic phase 軸**與**細胞演化軸**分離，避免把 parental homolog 差異誤稱 clone branching。

但 bulk ONT 對細胞裂解後的 DNA molecules 取樣，標準 HP tag 沒有 cell barcode。它通常知道一條 molecule 屬 HP1 或 HP2，卻不知道哪一條 HP1 molecule 與哪一條 HP2 molecule 原本來自同一細胞。因此 homolog phase、cell genotype 與 clone lineage 必須分層表述。

## 2. 兩個互相正交的判讀軸

| Allelic phase | 同一群細胞 | 不同細胞群 |
|---|---|---|
| **同一 HP** | cis local genotype | 不同 cells/clones 恰在同一 germline HP 累積不同 states |
| **不同 HP** | trans local genotype，例如 `10｜01` | HP1-only 與 HP2-only 的不同細胞群 |

Germline phasing 主要決定表格的「同一 HP／不同 HP」；cell-level linkage 才能直接決定「同一群細胞／不同細胞群」。只有§5.1全部條件成立時，互異同-HP exact states才可由物理限制間接推出不同來源細胞。

## 3. 重要術語

| 術語 | 中文詞義 | 本輪精確角色 |
|---|---|---|
| homologous chromosome | 同源染色體 | 二倍體 CN-neutral 時，一顆細胞通常各有一份 germline HP1／HP2 |
| haplotype | 單倍型／同一 homolog 上的 allele 組合 | HP family 所代表的 allelic background；不是 clone ID |
| PS / phase set | 相位集合 | 只有同 PS 或已驗證跨 PS orientation 的 HP1/2 才可直接比較 |
| HP tag | read-level 單倍型標籤 | 標記 molecule 的 haplotype background；不是 cell barcode |
| LongPhase-S `1-1/1-2` | HP1 背景上的 somatic haplotype read tag | 後綴是 read-level somatic haplotype 索引，不是跨 HP 通用 clone ID |
| read state | read 多位點狀態 | 依位置編碼為 R/A/X 或 0/1/?；例如 `AR=10` |
| full observation | 完整 read 狀態 | 同一 alignment 覆蓋該分析窗全部 sSNV；可直接指定 exact state |
| partial observation | 部分 read 狀態 | 含 X/?；只限制相容 state subcube，不能單獨確認 exact genotype |
| regional mutation state | 區域突變狀態 | 一條 homolog molecule 在本區域的 ALT-set；目前最安全的基本單位 |
| local cell-genotype class | 局部細胞基因型類別 | 只在本區域相同；不保證 genome-wide 是同一 clone |
| clone | 克隆細胞群 | 共同祖先衍生且共享一組可遺傳事件的細胞群；不是單一 read 或單一 HP |
| subclone | 亞克隆 | 相對祖先 clone 多出事件的後代細胞群；需細胞群存在與祖先關係兩項證據 |
| node | Boolean mutation state 節點 | plain node 可由 full state 支持；`H_*` 可為 hidden／Steiner／partial-supported，不等於 extant clone |
| topology / shape | 圖形幾何 | 表示 mutation-state tree 的包含／分枝關係；不能直接翻譯成 biological clone 數 |

## 4. Formal model：一個 clone 的兩條 homolog

### 4.1 命名決議

本文件中的 `10｜01` 情境不得再簡稱「模型 A」。Repo 既有 **Model A** 已專指「允許同一 bit 在不同分支 recurrence 的 Boolean-lattice state-tree solver」。為避免把 diploid cell-genotype hypothesis 與演算法模型混為一談，本情境統一稱：

> **H-trans-1clone：CN-neutral diploid 假設下的 single-clone trans candidate genotype。**

令區域有 \(k\) 個 sSNV。Clone \(c\) 的局部 diploid genotype 可寫為：

\[
G_c=(z_{c,1},z_{c,2}),\qquad z_{c,h}\in\{0,1\}^k
\]

其中 \(h=1,2\) 分別代表 germline HP1、HP2。兩位點案例：

\[
G_c=(10\mid01)
\]

表示 mutation A 在 HP1、mutation B 在 HP2，兩者為 trans。這是合法且已有文獻正式建模的 diploid clone genotype；PairClone 直接以「兩 homolog × 兩 loci」的 \(2\times2\) binary genotype 表示 subclone。

但 final genotype 不決定事件先後。以下歷史得到同一結果：

\[
00|00\rightarrow10|00\rightarrow10|01
\]

或：

\[
00|00\rightarrow00|01\rightarrow10|01
\]

因此「trans」是 allelic configuration；「A 早於 B」是 cellular lineage assertion，兩者不可混用。

## 5. 單一 HP 內，何時可推論不同來源細胞

### 5.1 條件式命題

假設對所有 contributing cells：

1. 在完整分析窗與所有 contributing cells 中，相關 germline HP 的 allele-specific copy number 恆為 1，且無 subclonal CNA；
2. variants位於same PS，或跨PS orientation已經獨立驗證，且read HP/PS assignment正確；
3. 每個exact state均由多條獨立full-span molecules支持；
4. base、mapping、alignment、sample contamination、phase-switch error、mapping／allelic sampling bias與callability差異已受控；
5. 採用標準穩定 cell-genotype 模型；不把 post-replicative lesion、同一細胞內 chromatid discordance 等未固定事件當 clone genotype。

則每顆細胞在該 HP 上只有一個固定 local state \(z_{c,h}\)。若觀察到兩個互異 exact states \(x\neq y\)，它們不可能同時來自同一細胞的唯一 HP copy，因此至少存在兩種來源 cell genomes／regional cell-genotype classes。

這可確認的是 **regional cellular heterogeneity**；仍不能只靠單區域判定完整 genome-wide clone 身分。

### 5.2 兩位點判讀矩陣

| 同一 HP 的 full-read states | Read 直接觀察／支持 | §5.1全部條件成立時 | 合理的演化口徑 |
|---|---|---|---|
| 只有 `11` | A、B 在同一 homolog molecule 共現 | 一種觀察到的 local state | cis co-occurrence；不能證明 subclone |
| `10`＋`11` | 兩種 homolog-specific exact states | 不同來源細胞 | `10→11` descendant/subclone-compatible |
| `01`＋`11` | 兩種 homolog-specific exact states | 不同來源細胞 | `01→11` descendant/subclone-compatible |
| `10`＋`01` | 兩種互不包含的 exact states | 不同來源細胞 | sister/independent-branch-compatible |
| `00`＋`10` | reference與mutation state並存 | 不同來源 cells；其中 `00` 可能代表正常背景 | `00` 可能是正常細胞，不能直接稱兩個 tumor clones |
| 只有 `H_11` | solver 需要／允許一個 completion | 未直接觀察 | 不可當 extant subclone |

令 state \(x\) 的 ALT-set 為 \(M(x)\)：

- \(M(x)\subset M(y)\)：在 infinite-sites、無 mutation loss、無 recurrence 下，\(x\rightarrow y\) 為 descendant-compatible。
- \(M(x)\nsubseteq M(y)\) 且 \(M(y)\nsubseteq M(x)\)：同一模型下不相容於單一路徑，為 branch/sister-compatible。
- 上述皆為模型相容性；CN、loss、recurrence 或錯誤未排除時，不可升格成真實 lineage。

此表只判讀指定 state pair 的 **pairwise compatibility**。若同一 HP 還存在其他 error-controlled states，必須把完整 state set 一起送入 solver；pairwise relation 本身不足以決定整體 topology。

## 6. 跨 HP 的 cell-overlap bounds

對 `HP1=10、HP2=01`，令：

\[
a=P(HP1=10),\qquad b=P(HP2=01)
\]

真正的 same-cell trans proportion 為：

\[
q=P(HP1=10,HP2=01\text{ in the same cell})
\]

在 balanced diploid 1+1、共同 cell denominator 與無 sampling bias 下，bulk HP marginals 雖不能唯一決定 \(q\)，但給出 Fréchet bounds：

\[
\boxed{\max(0,a+b-1)\le q\le\min(a,b)}
\]

四種 local cell-genotype proportions 可寫成：

| Local genotype | 比例 |
|---|---:|
| `10｜01` | \(q\) |
| `10｜00` | \(a-q\) |
| `00｜01` | \(b-q\) |
| `00｜00` | \(1-a-b+q\) |

因此：

- \(a=b=0.30\)：\(q\in[0,0.30]\)，one-trans-clone 與兩個完全互斥 clones 都相容。
- \(a=b=0.80\)：\(q\in[0.60,0.80]\)，完全互斥模型被排除，至少 60% cells 必須重疊。
- \(a=b=1\)：模型內 \(q=1\)，代表此區域 genotype 在 contributing cells 中固定；仍不保證只有一個 genome-wide clone。

實算 \(a,b\) 時必須使用兩位點皆 callable 的 common full-read denominator，不可直接混用兩個位點各自不同 denominator 的 site-specific AF。

## 7. 現行 InterSubMod 的實際能力

現行流程為：

```text
L0 germline HP family
  → L1 per-(region × HP-family) mutation-state tree enumeration
  → L2 post-tree CN annotation（目前主要處理 recurrence）
  → L3 bounded auxiliary annotation
```

關鍵實作語義：

1. `1/1-1/1-2 → family 1`，`2/2-1/2-2 → family 2`；不同 germline families 分開建樹。
2. Read 依區域內 sSNV 順序編碼成 R/A/X；只有完整覆蓋全部位點的 alignment 才進 full-state collection。
3. 每個 `region × HP-family` 的 solver node 是 Boolean mutation state；tree edge 是一個 bit 的 0→1 unit flip。
4. Region 結果是保留 HP component 的 ordered forest；跨 HP 不建立 mutation ancestry edge，也不把兩條 homolog 配回同一細胞。
5. `model-determined` 只表示 retained full/partial constraints 下候選唯一；不保證整棵樹有 full-span molecule，也不保證 hidden node 存在於取樣時的細胞群。
6. 現行 driver 的正式 claim scope 是 `regional mutation-state trees, not confirmed cell clones`。

因此最簡單的 `HP1=10、HP2=01` 在現行輸出是：

```text
HP-family1: ROOT → 10
HP-family2: ROOT → 01
region forest: (T_HP1, T_HP2)
```

這不是 `10` 與 `01` 的 sister tree，也不是已完成 cell-genotype pairing。

### 7.1 現有欄位名稱的語義校正

- `populations_by_hp`：實際是達 MINREAD 的 homolog-specific full-read patterns；不是已確認的 cell populations。
- `lineage-unit`：是 `region × HP-family` 的重建單位；不是已確認 biological lineage。
- `Single/no-within-HP-relation`：沒有在任一 HP component 看到分枝或多層路徑；不等於單一 clone。
- `Direct`／`Sister`：mutation-state graph geometry；只能稱 subclone-compatible／branch-compatible。
- `H_*`：hidden／Steiner／partial-supported state；不等於 hidden clone。

## 8. CNV 為何仍是硬邊界

只有在完整分析窗與所有 contributing cells 中 allele-specific CN=1、無 subclonal CNA，並使用共同 callable denominator、沒有顯著 mapping／allelic sampling bias時，HP-specific molecule fraction才可近似cell fraction（另見§5.1與§6）。一般CNV下：

\[
P(s\mid h)=
\frac{\sum_c\pi_c C_{c,h,s}}
     {\sum_c\pi_c\sum_x C_{c,h,x}}
\]

其中 \(C_{c,h,s}\) 是 clone \(c\) 中，源自 HP \(h\)、帶 state \(s\) 的 physical copy 數。此時 HP 是 parental haplotype family，不再是一對一的 physical copy ID。

| CN 情境 | 對同-HP多 state 的影響 |
|---|---|
| balanced allele-specific CN=1+1 | 只是§5.1的必要條件之一；§5.1全部條件成立後才可使用「互異同-HP state → 不同來源細胞」推論 |
| HP-specific amplification | 同一細胞可有多份同源自 HP1 的 copies，且 copies 可帶不同 states |
| copy-neutral LOH 2+0 | total CN仍是2，但同一 parental HP可能已有兩份 copies |
| deletion／subclonal loss | mutation state 可在後代消失，破壞單調 ALT-set 模型 |
| subclonal CNA | molecule abundance不再等於 cell prevalence |
| male chrX／chrY | 通常不具兩條 homolog，`HP1｜HP2` diploid model 不適用 |

所以需要的是整段、allele-specific、最好 clone-aware 的 CN；不能只用 total CN=2 或 region midpoint 的 `neutral` 標籤當充分證明。

## 9. Claim ladder：可說與不可說

| Evidence level | 可使用的正式說法 | 不可升格 |
|---|---|---|
| E1：HP assignment | mutation/read 被分到 germline HP1 或 HP2 | 父源／母源，除非有 pedigree orientation |
| E2：full-read joint state | 同一 HP molecule 直接支持 `10/01/11` | cell clone |
| E3：同-HP多 exact states＋§5.1全條件 | same/verified PS、每state多條獨立full-span molecules、錯誤與偏差受控、穩定genotype模型、完整窗逐細胞HP-copy=1且無subclonal CNA時，支持regional cellular heterogeneity | genome-wide clone 數 |
| E4：nested/incomparable ALT-sets | descendant/subclone-compatible 或 branch-compatible | 真實祖先／姐妹關係已確認 |
| E5：跨 HP bounds | same-cell overlap 的上下界；必要時排除完全互斥模型 | 唯一 diploid clone pairing |
| E6：多樣本 CN-aware mixture | model-supported clone genotypes/frequencies | calibrated biological truth，除非有外部 truth |
| E7：single-cell／colony | 相同細胞中兩條 homolog 的 joint genotype | 仍需足夠細胞與技術錯誤控制 |

### 9.1 本階段核准句

> 利用 germline HP 分群後的同-read sSNV 聯合狀態，我們可在單一 homolog family 內辨識多個 regional mutation states。當 variants 位於same PS或跨PS orientation已獨立驗證、每個state均由多條獨立full-span molecules支持、錯誤與取樣偏差受控、符合§5.1的標準穩定cell-genotype模型，且完整分析窗中所有contributing cells的相關HP均為allele-specific CN=1、無subclonal CNA時，互異states可推論來自不同來源細胞；mutation-set的包含或互不包含關係分別提供descendant/subclone-compatible或sister/branch-compatible evidence。此方法適合建立區域clone/subclone candidates並支援後續CN-aware、多樣本與single-cell驗證，但目前不等同確認完整biological clone或真實演化祖先樹。

### 9.2 禁止句

- 「Single-only = 一個 clone。」
- 「Direct-only = 已確認 subclone。」
- 「Sister-only = 已確認兩個姐妹 clones。」
- 「每個 read state 就是一個細胞群。」
- 「HP1-1／HP1-2 是已確認 clone 1／clone 2。」
- 「Hidden node 是取樣時仍存在的 hidden clone。」
- 「VAF 第一順位就是最可能真實 biological tree。」
- 「total CN=2 已排除 copy-number confounding。」

## 10. 後續可計算的分析流程（尚未執行）

1. **Phase gate**：限定 same PS；跨 PS 必須有獨立 orientation/stitching 驗證。
2. **Full-state census**：逐 `region × HP-family` 輸出未被 MINREAD 刪除的 raw `00/10/01/11/...` counts與兩位點皆 callable denominator。
3. **Error gate**：去除低 BQ/MQ、secondary/supplementary重複、phase conflict與alignment identity conflict；以獨立 molecules計數。
4. **Allele-specific CN gate**：要求整個分析窗 CN=1+1；gain、LOH、missing/misfit 分開標記，不能預設 neutral。
5. **State heterogeneity分層**：只有§5.1全部條件成立時，同HP至少兩個互異exact full states（可含`00`）才標`regional_source_cell_heterogeneity_supported`；若至少兩個皆為mutation-bearing states，再標`regional_mutation_bearing_cellular_heterogeneity_supported`。前者可能只是normal/tumor mixture，不等於兩個tumor clones。
6. **Geometry**：ALT-set nested → `descendant_compatible`；incomparable → `branch_compatible`；partial/hidden-only → `model_only`。
7. **跨 HP bounds**：計算 \(q_{min},q_{max}\)，並以 molecule bootstrap／beta-binomial interval傳播不確定性。
8. **多樣本一致性**：檢查 candidate states是否在獨立樣本／時間點同步變動；不可把同一 cell line 的不同 basecaller當 biological replicate。
9. **Orthogonal validation**：優先挑選高支持、CN乾淨、可辨識性高的 region做 single-cell DNA或single-cell-derived colony驗證。

## 11. 已知前提與未決問題

### 已驗證／已確立

- 現行流程確實先按 germline HP family 分開 reads，再各自建 regional mutation-state trees。
- Full observation 與 partial observation 在 solver 中有不同證據地位。
- 跨 HP component 不建立 ancestry edge，也不提供 cell pairing。
- PairClone 類方法允許以兩 homolog × 多 loci 的 latent clone genotype建模。
- CN/mutation multiplicity與mutation loss會混淆 bulk VAF、CCF與phylogeny。

### 尚未驗證／本輪不可假設

- 每個現行 region 均有 same-PS、低 switch-error 的完整 phase continuity。
- 每個 region 的 allele-specific CN均為1+1，且沒有 subclonal CNA。
- 現行 `populations_by_hp` 已保留計算 prevalence 所需的所有低頻 raw states。
- `Direct`／`Sister` geometry對應的 cell lineage已由single-cell或synthetic truth確認。
- VAF ranking分數是 calibrated probability或獨立validation。
- 現行全樣本 topology比例可改報為clone/subclone發生率。

## 12. 外部研究背景

1. [PairClone: A Bayesian Subclone Caller Based on Mutation Pairs](https://academic.oup.com/jrsssc/article/68/3/705/7058355) — 直接以兩 homolog × 兩 loci 的 \(2\times2\) genotype表示 latent subclone；主分析限制copy-number-neutral區域。
2. [Somatic mutation phasing and haplotype extension using linked reads](https://pmc.ncbi.nlm.nih.gov/articles/PMC11326269/) — joint mutation states比marginal VAF多資訊；haplotype resolution仍不等於cell identity或clone truth。
3. [Long-read phasing of lung cancer genomes](https://pmc.ncbi.nlm.nih.gov/articles/PMC9203510/) — germline anchors可將somatic mutations定相到HP，但研究本身未提供bulk cross-homolog cell pairing。
4. [DeCiFering the elusive cancer cell fraction](https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/) — allele-specific CN、mutation multiplicity與loss會混淆CCF及下游phylogeny。
5. [VCF v4.4 specification](https://samtools.github.io/hts-specs/VCFv4.4.pdf) — PS是可彼此定相的genotype集合；不同phase sets不可默認共享相同HP orientation。

## 13. 本輪閱讀建議

- [目前研究焦點與正式 claim scope](/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md:21) — clean-v3 lifecycle與 `region tree ≠ biological clone`。
- [HP／LOH 術語校正手冊](/big7_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260316_HP_LOH術語校正手冊_01.md:10) — read-level HP、variant-level線索與region-level LOH的分層。
- [Formal problem statement](/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md:32) — full/partial observation、既有 Model A與minimal-compatible tree定義。
- [Layered data model](/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md:98) — region、HP-family、lineage-unit、tree與node的操作型單位。
- [HP-family與R/A/X read-state實作](/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py:85) — HP-family mapping與full-read state建立。
- [Layered tree reconstruction driver](/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py:376) — per-family solver、CN post-hoc角色與正式claim scope。
- [LongPhase-S本地知識庫](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md:183) — 完整9態HP詞彙、somatic haplotype與限制。

## 14. 執行與驗證紀錄

### 輸入

- Repo root：`/big7_disk/liaoyoyo2001/InterSubMod`
- 本輪對話中確認的HP、same-read、clone/subclone、CN與identifiability結論
- §13列出的本地SoT與§12原始研究

### 建立命令

```bash
mkdir -p InterSubMod/research/20260713_hp_readstate_clone_identifiability
# 文件內容以 apply_patch 新增，未覆寫既有檔案
```

### 輸出

- `InterSubMod/research/20260713_hp_readstate_clone_identifiability/20260713_HP分群Read共現與區域Clone可辨識性結論_01.md`

### 驗證狀態

- 完成時間：`2026-07-13T05:57:05+08:00`。
- Markdown結構／metadata／路徑：**PASS**；本地引用檔案 `7/7` 存在，code fences `8`（偶數、成對）。
- claim ceiling關鍵詞檢查：**PASS**；CN gate縮寫過度宣稱掃描為0，所有「不同來源細胞」推論均回指§5.1完整gate。
- 無上下文reader測試：前3輪依序找出摘要條件壓縮、stable-genotype遺漏與gate不一致，均已修正；第4位全新reader final gate：**PASS**。
- Worktree範圍：只新增本目錄與本文件，未修改既有檔案；repo原有dirty changes均未觸碰。
- 未執行生物數據重算；本文件是方法學current decision record，不是new experimental result。

### 實際QA輸出片段

```text
LOCAL_REFERENCES PASS
CODE_FENCES 8 PASS
CN1_OVERCLAIM_SCAN PASS
FRESH_READER_FINAL PASS
```
