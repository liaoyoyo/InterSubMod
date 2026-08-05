<!--
建立時間: 2026-07-25
目標: 固定 exact-PS×HP subclonal reconstruction 全流程的代號、每步數學定義、計算範圍與 claim ceiling
處理範圍: 7 technical datasets / 6 biological IDs / GRCh38 chr1-22；彙整既有 validated 定義，不產生新數字
類型: 定義規格（definitions & notation SoT）
data_sources:
  - InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/data/20260724_exactPS_k_HP與分母重算紀錄_01.md
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/20260724_exactPS最佳樹拓撲Signature精確統計_01.md
  - InterSubMod/research/20260725_crosssample_topology_matrix_report/results_v2/measure_definitions.json
  - InterSubMod/research/20260724_all_ssnv_representative_locus_explanation/20260724_單位點sSNV甲基多群代表位點定義與停止判定_01.md
-->

# exact-PS×HP 流程：代號、每步數學、計算範圍

> 本檔不產生任何新數字。所有數值均自上列 data_sources 逐一讀回抄錄，供跨 session 統一口徑用。

---

## §0 計算範圍（scope）— 所有數字的共同前提

| 維度 | 範圍 | 排除了什麼 |
|---|---|---|
| 參考基因體 | GRCh38 | — |
| 染色體 | **chr1–22（autosomes）** | chrX / chrY / chrM |
| 資料集 | **7 technical datasets / 6 biological IDs** | HCC1395 與 HCC1395_DORADO 是同一生物株的兩個技術來源 |
| 變異母體 | 同一次 **LongPhase-S recalibrated `FILTER=PASS` biallelic sSNV** | raw ClairS 全部候選；非同批 recalibration |
| Haplotype | **HP1 = raw {1, 1-1, 1-2}；HP2 = raw {2, 2-1, 2-2}** | HP3、HP4、unphased、missing PS |
| 比對品質 | mapped **primary** alignment、**MAPQ ≥ 20** | unmapped / secondary / supplementary / QC-fail / duplicate |
| 鹼基品質 | fixed R/A 需 **BQ ≥ 20** 且鹼基等於 REF 或 ALT | 其餘記 X/L/O/D/S，不建邊 |
| Edge 門檻 | **distinct canonical molecule support ≥ 3**（primary） | sensitivity 另跑 1 / 2 / 5 |
| 樹規模上限 | 計算切割 **k ≤ 12** | k>12 的 source region 仍保留，只是切成 blocks |
| 未做的修正 | **CN、LOH、purity、multiplicity、CCF 校正皆未納入** | 因此不可由形狀或 HP 數量宣稱 clone prevalence |

**顯示別名**：互動工作站對外顯示 `HCC1395_HKU` / `HCC1395_NYGC`；內部 authority key 與檔名仍為 `HCC1395` / `HCC1395_DORADO`。別名只在 presentation layer，不參與 join。

---

## §1 代號表（notation）

### 1.1 資料層

| 代號 | 全名 | 精確定義 |
|---|---|---|
| `S` | candidate sites | LongPhase-S PASS autosomal biallelic sSNV loci 總數（all7 = 469,849） |
| `P` | phase set | `(dataset, chromosome, exact nonmissing PS)`；PS 數字只在同一 dataset×chromosome 內有意義 |
| `H` | container | `(dataset, chromosome, exact nonmissing PS, primary HP family)` — **建 linkage graph 的容器** |
| canonical molecule | — | `SHA256(dataset, RG-or-".", QNAME)`；calls/quality/MAPQ/HP/PS 全一致才可 collapse，衝突則 fail-closed |
| fixed `R`/`A` | — | 該 endpoint BQ≥20 且鹼基 == REF（R）或 ALT（A） |
| `X/L/O/D/S` | — | 未覆蓋 / 低品質 / 其他鹼基 / deletion / ref-skip — **均不建邊** |
| `RR/RA/AR/AA` | — | 一條 molecule 對一組 endpoint pair 的觀測型態，**依 genomic 座標順序編碼**（不是祖先→後代方向） |
| allele key | — | `chrom:pos1:ref:alt` |

### 1.2 圖層（🔴 兩種「節點」不可混用）

| 層 | 圖 | 節點是什麼 | 邊是什麼 |
|---|---|---|---|
| **L-linkage** | sSNV read-linkage graph | **一個 sSNV 位點**（membership） | 同一 molecule 在兩端都 fixed R/A |
| **L-tree** | Boolean hypercube / mutation-state graph | **k 位元的 R/A 狀態向量**（如 `000`、`100`、`110`） | Hamming distance = 1 的一次突變 |

| 代號 | 定義 |
|---|---|
| membership | `(locus, exact PS, HP)`；**同一 physical locus 可重複屬於不同容器** → memberships ≠ unique sSNV |
| component | 容器內由合格 edges 形成的 maximal connected component |
| `k` | component 內的 sSNV 位點數 |
| **`W_i`** | `k_i ≥ 2` 的 component = **唯一被稱為 read-linked candidate region 的單位** |
| k=1 abstain | `ABSTAIN_SINGLETON_UNLINKED` — 此 membership 在該容器、threshold=3 下無合格 edge；保留在 audit registry，**不進下游分母** |
| `start/end/span` | `min/max(position in W)` 與其差 — **只是顯示用座標包絡**，W 是位點集合不是連續區間 |
| **`B_ij`** | `W_i` 經 k≤12 限制切出的 **computational block**。🔴 `B ≠ W`，不可把 block 數當 region 數 |
| final group | 通過 post-cut 與 pattern support 後，實際進入建樹的單位（all7 = 98,955） |

### 1.3 樹層

| 代號 | 定義 |
|---|---|
| `X = {0,1}^k` | 狀態空間；根 `r = 0^k` |
| unit-flip 邊 | `x → y` 合法 ⟺ `y = x + e_j`（恰多一個 1，方向 0→1） |
| `O ⊆ {0,1,?}^k` | 觀測集合；`?` = 未觀測 |
| 相容 `x ~ p` | 對所有 `p[j] ≠ ?` 有 `x[j] = p[j]` |
| Model A | **本專案採用**：同一位元允許在不同分支重複 0→1（recurrence） |
| Model B | 全樹每位元 ≤1 次翻轉；**現行 C++ `classify()` 實作的是 Model B** |
| `m_f(x)` | 獨立多重度通道 — Model A 的**硬前置需求**，用來拆「多重度 artifact」vs「真 recurrence」 |
| 通道 A（`φ`） | 外部控制集錨定的 2-著色，只做互斥/sibling 過濾 |
| 通道 M（`d`） | 實值輪廓距離（**含甲基化**），只做無結構負篩 / 相容集內弱標註 / 佐證已定邊 |
| ranked_complete | 已完成 ranking 的 unit（all7 = 71,955）— **拓撲 determinacy 的分母** |
| global-best tree | 在 exact rational AF score 下並列最佳的完整樹（all7 = 680,527 棵） |

---

## §2 每一步的意義與數學

### Step 1 — 建容器 `H`

把 sSNV 依 `dataset × chromosome × exact nonmissing PS × primary HP family` 分箱。

**為什麼**：PS 只提供 LongPhase-S 定義的**局部 HP orientation 座標系**。不同 PS 之間 HP1/HP2 的相對方向在本層**未定義**，所以跨 PS 聚合同一數字 HP 是錯的（舊版 production 的缺陷）。

### Step 2 — 建 read-linkage graph（L-linkage）

- 節點：容器中至少被一條 canonical molecule 以 fixed R/A 觀察到的 sSNV membership
- 邊：同一 canonical molecule 在**兩個 endpoint 都是 fixed R/A**
- 投票：每條 molecule 對同一 endpoint pair **最多投一票**，分類為 RR/RA/AR/AA
- 保留條件：`total support = RR+RA+AR+AA ≥ 3`

**🔴 距離不參與連邊。** 三條 molecule 若都 fixed-observe 相距 >50 kb 的兩端，該邊仍合法。50 kb 只作 QC。

**Partial read 規則**：對排序位點 `s1<s2<s3`，pattern `R-X-A` **只**對 `(s1,s3)` 投一票 `RA`；不會因 read 的座標 span 被動補 `(s1,s2)` 或 `(s2,s3)`。

**Transitivity**：`(s1,s2)≥3` 與 `(s2,s3)≥3` 即使 `(s1,s3)=0`，三點仍併成同一 component，且兩條邊可由**不同 molecule 集合**支持。
→ **W 不代表有一條 read／一個細胞／一個 clone 橫跨全部 k 個位點。**

### Step 3 — component 切分與 k=1 排除

```
all7:  255,752 source components = 170,131 (k=1 abstain) + 85,621 (W, k≥2)
       613,480 active memberships = 170,131 + 443,349 (W memberships)
```

**排除必須在 membership 層執行**：HCC1395 有 5,717 個 physical loci 在某個 PS×HP 是 singleton、在另一個 PS×HP 卻有合格連接。不可把該 genomic locus 全域刪除。

### Step 4 — bounded partition（k ≤ 12）

```
all7:  85,621 W → 99,966 bounded blocks
       98,955 final groups = 99,966 − 419 (post-cut k=1) − 592 (pattern unsupported)
       439,685 final memberships = 443,349 − 419 − 3,245 (unsupported memberships)
       保留率 99.17%
```

每個 block 必須保留 `source_component_id` 以回聚到 `W_i`。final groups 的 original k 全部落在 **k = 2–12**（`k=1 = 0`、`k>12 = 0`）。

### Step 5 — mutation-state tree 枚舉（L-tree）

在 Boolean lattice 上找 rooted directed tree：根 `0^k`、每邊一位元 `0→1`、每筆觀測至少相容於樹上一節點。

**最小化為字典序**（不是單一 cost function）：

```
(1) 解釋所有 O  ≻  (2) min 隱藏(Steiner)節點  ≻  (3) min 邊數
                  ≻  (4) max 觀測支持  ≻  (5) 列出所有同分
```

**覆蓋公理（硬約束）**：完整觀測 `p ∈ {0,1}^k` ⇒ `p ∈ V`。solver **不得**為求更小的樹而丟棄任一通過雜訊門檻的觀測；覆蓋失敗 → 判 `conflict`，不是默默省略。

**🔴 三條紅線**：**枚舉（enumerate），非最佳化（optimize），非排序（rank）**。
`∂(tree set)/∂φ = ∂(tree set)/∂d = 0` —— 甲基化與 HP **無條件不改變樹集**，這是非循環性的核心保證。

**複雜度正當化**：`k ≤ 12` 有界 ⇒ `|X| = 2^k` 有限 ⇒ 窮舉 tractable。NP-hardness 是 `k → ∞` 的漸近性質，在 k 有界時不生效。solver 踩的是「小 k 窮舉島」，不是「近似大 k 最佳化」。

**「定不出來」即答案**：determined ⇒ 唯一最小樹；否則 ⇒ 輸出全相容集 + 不確定性標註。

### Step 6 — exact AF 最佳樹展開

1. 以 **exact rational AF** 計算每個 child 的所有 Hamming-1 parents
2. 保存該 child 的**所有**最佳 parents
3. 展開各 child 最佳 parent choices 的 **Cartesian product**
4. 收集每棵 global-best tree 的 canonical signature 與 coarse class

本批每 unit 最多 **276 棵** global-best trees → 可直接精確展開，不需近似或抽樣。

### Step 7 — topology signature 與 determinacy 分類

**Signature 建構**：從 root 遞迴轉成括號字串；同一 parent 的 child signatures **排序後**串接。

因此 signature：**保留** root 與 parent-child geometry；**忽略** mutation/site label；**忽略** sibling 順序。
→ `00→10→11` 與 `00→01→11` 都是同一條三節點 chain `((()))`，雖然 acquisition label 不同。

| 類別 | 條件 | all7 |
|---|---|---:|
| `UNIQUE_TREE` | 最佳完整樹只有一棵 | 39,648 / 71,955 = **55.1011%** |
| `TIED_SAME_TOPOLOGY` | 多棵但 rooted-unlabeled signature 只有一種 | 23,858 / 71,955 = **33.1568%** |
| `TIED_CROSS_TOPOLOGY` | 多棵且 signature 超過一種 | 8,449 / 71,955 = **11.7421%** |
| **單一 exact topology** | 前兩者合計 | 63,506 / 71,955 = **88.2579%** |

🔴 舊 legacy grouping 的 `92.18%` **不可**移植到新版；同口徑值是 `88.26%`。

**四類 coarse graph geometry**：

| 類別 | 定義 |
|---|---|
| `Single-only` | 沒有 branching，也沒有 root 之下兩層以上路徑 |
| `Sister-only` | 有 branching，但沒有兩層以上路徑 |
| `Direct-only` | 有兩層以上路徑，但沒有 branching |
| `Sister+direct` | 兩者同時有 |

可唯一指派 coarse class：`63,511 / 71,955 = 88.2649%`。其中只有 **5 units** 是「exact shape 超過一種但 coarse class 相同」。

**這四類是 mathematical graph geometry，不等同 cellular single-clone / multi-clone / confirmed subclone。**

### Step 8 — 跨樣本比較

| 度量 | 數學定義 |
|---|---|
| coordinate signature | `(chrom, sorted active_positions)` — **只用座標，不含 REF/ALT** |
| fail-closed 1:1 | 共享 signature 只有在兩邊各恰好一個 ranked unit 時才可評估；重複 key 排除 |
| candidate sites | 所有 `site_catalog` rows |
| active sites | threshold=3 membership 中的 unique allele keys |
| W sites | threshold=3 tree_eligible components 中的 unique allele keys |
| primary edges | `passes_primary_threshold=true` 的**無向** allele endpoint pairs |
| exact components | threshold=3 tree_eligible components 的 exact allele-member 集合 |
| phase invariance | **PS 與 HP labels 跨資料集刻意排除** |
| TVD similarity | `1 − 0.5·Σ|Δp|` |
| JS distance | `sqrt(JSD_nats / ln 2)`；**JS similarity = 1 − JS distance** |
| `js_divergence_complement` | `1 − JSD_nats/ln2` — 保留但**明確不等於** 1 減 JS distance |
| directional containment | `M[row, col] = \|row ∩ col\| / \|row\|` — 有方向，不必對稱 |
| exact label permutation | 6 biological IDs 枚舉 `6!/(3!·2!·1!) = 60` 個 3/2/1 分組；統計量 = 4 個 within-group pairs 平均 − 11 個 between-group pairs 平均；one-sided p = 統計量不小於觀察值的 assignment 比例 |

**技術來源處理**：涉及 HCC1395 的 biological pair，先把 HCC1395 與 HCC1395_DORADO 的 technical-pair 值 macro-average，再進 permutation。

### Step 9 — 甲基化軸（平行、有界、非骨幹）

| 名詞 | 精確定義 | 不能回答 |
|---|---|---|
| focal ALT reads | 在指定 sSNV 位點支持 ALT allele 的 reads | 這些 reads 必然來自同一 clone |
| positional singleton | 同 dataset/chrom、相鄰 gap ≤ 50 kb 的 transitive component 中 size = 1 | 全基因組 read-sharing degree = 0 |
| **M1** stable multigroup | 甲基距離分群**跨 seeds 穩定** + 通過 **column-shuffle null** + 最小群門檻 | 技術/生物 confound 已排除 |
| **M2 PASS** | strand / start / end / length / MAPQ / CpG-called / latest HP 等已測軸**皆可判定**，且**無軸與群組對齊** | 未量測 confound、CN/purity、細胞身分已排除 |
| latent molecular substructure | 同一 focal ALT read 集內，穩定但尚未綁定細胞身分的甲基分層 | 等同已確認 subclone |
| Bernoulli read distance | 兩條 read 在**共同可比較 CpG** 上的甲基狀態差異 | 突變時間或演化方向 |
| UPGMA | 對 read×read 甲基距離做 average-linkage hierarchy | cellular phylogeny / mutation tree |
| HP/PS | 最新 same-run LongPhase-S sidecar 的 haplotype/phase-set tag | HP1/HP2 等於 clone 1/clone 2 |

**分母（皆為 operational screen yield，非 clone prevalence）**：

| 問題 | 分子/分母 | 比例 |
|---|---:|---:|
| 全 sSNV M1 穩定多群 | 102,842 / 469,849 | 21.8883% |
| singleton M1 | 5,961 / 50,432 | 11.8199% |
| singleton M2 PASS | 30 / 50,432 | 0.05949% |
| 選入且可完成 M2 者的 PASS 率 | 30 / 48 | 62.5%（高度選擇後條件比例，**不可外推**） |

---

## §3 分母守恆式（all7）

```
255,752 source components  = 170,131 singleton + 85,621 W
613,480 active memberships = 170,131 singleton memberships + 443,349 W memberships
 98,955 final groups       =  99,966 bounded blocks − 419 post-cut k1 − 592 pattern unsupported
439,685 final memberships  = 443,349 − 419 − 3,245 unsupported memberships
 71,955 ranked_complete units  →  680,527 global-best trees
```

cohort checks 全部 `true`（`all_pass: true`）。

---

## §4 三個 50 kb 指標必須分開（不可互相代稱）

以 HCC1395 chr1–22 為例：

| 指標 | 精確定義 | 數量 | 分母 |
|---|---|---:|---:|
| direct edge > 50 kb | qualifying endpoint pair 兩端座標距離 > 50 kb | 47 | 76,202 retained direct edges |
| W span > 50 kb | `max(pos) − min(pos) > 50 kb` | 1,064 | 11,462 W |
| W adjacent gap > 50 kb | W 成員排序後至少一個相鄰 gap > 50 kb | 4 | 11,462 W |

1,064 個長 span W 中只有 4 個有相鄰 gap > 50 kb → **多數長 W 是由多段合格 read chain 連成，不是由 50 kb 幾何區間定義**。

（all7 尺度的對應數字見 `InterSubMod/docs/CURRENT_FOCUS.md` §2026-07-23：16,537 條 direct edges > 50 kb；移除後 1,172 個 W partition 改變、228 個 W 完全失格、961 memberships 回到 k=1；W 加總因 component splitting 由 85,621 變 85,872。兩者 scope 不同，不可混用。）

---

## §5 Claim ceiling（統一口徑）

| 已成立 | 不成立 |
|---|---|
| L1 exact-PS×HP strict read-linkage **7/7** | production strict **directed** topology **0/7** |
| 88.26% ranked units 有單一 exact rooted-unlabeled topology | 一個 tree node = 一個 cellular clone |
| coarse geometry 可唯一指派 88.26% | Sister = 多 clone；Direct = 已確認 subclone ancestry |
| biological-ID specificity 極強 | 每條精確邊都穩定（primary-edge Jaccard 僅 0.198050） |
| M1/M2 位點可定位、可回查、可重現 | confirmed cellular subclone（= 0）、clone 數、linear ancestry（= 0） |

**W 是 read-linkage region；endpoint edge 是無向 linkage，不是 evolutionary parent→child。**
**AF 最佳樹不是 calibrated biological posterior；CNA / LOH / allele-specific amplification 未被排除。**

---

## §6 待確認的定義風險（open）

1. **`MAX_SNV` 上界不一致**：formal spec §7.5.2 寫 `k ≤ 8 (MAX_SNV)`，production bounded partition 用 `k ≤ 12`。需確認何者為現行權威。
2. **Model A 尚未落地**：現行 C++ `classify()` 是 Model B。改採 Model A 需要獨立多重度通道 `m_f(x)`（非由 B 導出），該通道目前不存在。
3. **memberships ≠ unique sSNV**：active/HP memberships 可重複計入同一 genomic locus，不可與 unique sSNV 混用分母。
4. **「合計」列含技術重複**：7 technical datasets 的合計含 HCC1395 + HCC1395_DORADO，**不可解讀成 7 個獨立生物樣本**（真 n = 6 biological IDs）。
5. **上游 S→active 漏斗未稽核**：chr6 與 chr16 合計占 upstream S 的 57.6229%，卻只占 active loci 的 9.2073%。論文引用此漏斗前必須另稽核 upstream VCF、exact PS 可用性與 fixed R/A callable 率，**不可直接解釋為生物差異**。
