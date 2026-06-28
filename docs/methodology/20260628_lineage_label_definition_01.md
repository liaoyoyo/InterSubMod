<!--
建立時間: 2026-06-28
類型: methodology 定義 — clone/subclone lineage 標籤通用定義 (HP{h}-path) + situation 分類 + 甲基裁決
狀態: in_progress(定義定稿 + chr17 驗證通過;genome-wide v1 標籤產出)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data
provenance: 凍結資料 @ branch feat/summary-nreadsvalid@5308d9e;甲基裁決據 chr17 per-read 實測
-->

# Clone/Subclone Lineage 標籤通用定義（HP{h}-path）+ ≥3 位點處理

> 框架：Diátaxis(reference/定義) + scientific-rigor §2(每節標 L1-L5)。每數字 grep-able。

## §0 一句話定義

**每條 read / 每個克隆節點標 `HP{h}` + Dewey 路徑 `-{b1}-{b2}-…`**：h = germline haplotype 根（H1/H2/H3-unphased）；路徑 = 該 read 在 haplotype h 的 somatic 樹上從根往下的 lineage（nested=往下加層、sibling=同層新分支、co_linked=同節點；分支編號 = VAF 遞減）。

## §1 理論依據（為何 pairwise 合法）[L1 理論]

sSNV = 二元字元（REF=祖先/ALT）。**Four-Gamete / Perfect-Phylogeny 定理**（Estabrook 1975；Gusfield 1991）：一組二元字元存在完美系統發生樹 ⟺ **每一對都相容**（無四配子違反）。→ **成對相容對二元字元即足以保證整棵樹存在；不需單分子整跨。** `classify()` 的 `independent`（四格全有）= 四配子違反偵測。
- 實證：≥3 區 **99.4% 相容**（3,414/3,436；cycle 不相容僅 0.6%）。

## §2 2-locus 基礎（HP 鑑別）[L1]

每對 sSNV 看共讀 read 的 2×2（RR/RA/AR/AA），cell 為真 ⟺ `count > coread×2%`（ε=2% 定案，見 `20260628_sSNV_linkage_threshold_decision_eps2_01.md`）：

| 2×2 模式 | 關係 | + HP | 標籤 |
|---|---|---|---|
| AA only | co_linked | — | 同節點（同 lineage）|
| 單側零格(AR=0 或 RA=0) | nested | same-HP | 祖先→後代（往下加層）|
| 互斥(AA=0, RA+AR>0) | sibling/allelic | **diff-HP → allelic** | `H1 vs H2`（**非 subclone**）|
| 互斥 | sibling/allelic | **same-HP → sibling subclone** | `H{h}-1 / H{h}-2`（VAF 高=1）|

- diff-HP 互斥佔 mutual_excl **59.1%**（其中 88.2% 乾淨 H1-vs-H2、11.8% 涉 HP3）；same-HP 互斥 H1 1322 / H2 1318 / H3 185。
- VAF 排序可分 **72%**（CN-clean 才可信）；28% tie 標 `?`。

## §3 通用定義（遞迴，推廣到任意位點數）[L1]

**Label(node) = `H{h}` + `-{b1}-{b2}-…-{bk}`**
- **h**：germline haplotype 根。正常人 = 2 根 H1、H2（HP tag："1-1"→H1、"2-1"→H2、"3"→H3-unphased）。
- **路徑**：somatic 樹從根往下：
  - level-1 somatic 事件（nested DAG 的根、無入邊）= `H{h}-1, H{h}-2, …`（VAF 遞減編號）
  - nested 後代 = 在父路徑後接新數字 `…-j`
  - sibling = 同層另一分支（同位新數字）
  - co_linked = 同一節點（同標籤）
- **read 標籤** = 它攜帶突變所達的**最深節點**；germline（全 REF）= `H{h}` 根。
- **深度 k** = read 攜帶的 nested somatic 事件數。

**遞迴 → 天然推廣**：1 突變→`H{h}-1`；2 突變分叉→`H{h}-1`/`H{h}-2`，巢狀→`H{h}-1-1`；k 突變→最深 k 層（`H1-1-1`/`H1-2`/`H2-1-1`/`H2-2`… 皆合法）。

### chr17 worked example 驗證（per-read，通過）[L1]
3 sSNV 全 H1：α(48365089,VAF 0.82)=祖先、β1≡β2(co_linked,VAF~0.48)=後代。
- 全 REF → `H1`（germline 根）：**11 reads**（L0 9 + other 2）
- α-only → `H1-1`：**20 reads**（=L1 ✓）
- α+β → `H1-1-1`：**19 reads**（=L2 ✓）
產生器 `scripts/build_lineage_labels.py` 輸出與 lineage_counts 一致。

## §4 ≥3 位點：先分 situation 再處理（必須）[L1]

≥3 不可一視同仁；先依可驗證判準分類：

| Situation | 判準 | 處理 | 證據層級 | v1 計數 |
|---|---|---|---|--:|
| **A 單分子整跨** | span≤34kb + 有跨全位點 ALT read | 直接讀基因型向量定樹 | 最強 | 1,549 |
| **B 可整跨/pairwise** | span≤76kb | pairwise 建樹 + 相容檢定 | 強 | 1,603 |
| **C 必鏈接** | span>76kb | pairwise transitivity；**沿 read-gap 切子區** | 中（統計相位）| 668 |
| **HP-mixed** | 區內混 H1/H2/H3 | 各 haplotype 分樹；HP3 待 cis-control | 分層 | （標於 haplotypes 欄）|
| **不相容** | cycle/independent | flag，不強建樹 | 排除 | 22 |

> ⚠ v1 tier 門檻為啟發式（A 的「multi」採 ≥2 distinct pop + ≥1 ALT），與「multi_supported 嚴格 ≥2 ALT vectors（20%）」定義略異，待精煉統一。

## §5 甲基是否有用 — chr17 實測裁決 [L1]

chr17 L1(α-only) vs L2(α+β)：16/67 CpG 顯著差異（看似能分 lineage），**但**：
- 甲基 k=2 分群實際分的是 **germline(L0) vs mutant**（cluster2={L0:7,L1:1}、cluster1={L1:19,L2:19 混}）。
- best axis = **ALT@α 對齊 23 CpG ≫ lineage_L1_vs_L2 只 6** → 甲基結構**對齊 α-genotype 軸（cis-ASM），非獨立 subclone 軸**。

**裁決**：
- 🔴 **甲基不能當 HP3 定相 / 分叉偵測的獨立驗證器**——被突變的 cis 效應主導（double-dip）。
- Part 1（HP3 甲基定相）：**只在有 germline-ASM（H1≠H2 穩定甲基差、與突變無關）的區可行，且須 normal cis-control 確認**（T-GATE-GB）；非通用。
- Part 2（甲基確認分叉）：genetics 已判處冗餘、genetics 不能判處被 cis 汙染 → **不提供可靠獨立確認**。
- 甲基定位維持 **bounded-auxiliary / off-ladder**（corroborate 非 detect）。

## §6 通用定義成立的前提（誠實邊界）[L1/L3]

1. **haplotype 可靠**：HP3（somatic-unphased，linked-somatic 1317/20815=6.3%）標 `H3?`，待 cis-control 定相。
2. **樹被觀察決定**：pairwise 觀察完整度中位僅 **33%**（observed same-HP pairs / C(n,2)）→ 稀疏處拓樸欠定，報「一棵相容樹」非唯一。
3. **VAF 排序**：僅 CN-clean、72% 可分；CN-altered/tie 標 `?`。
4. **深度 ≥3 路徑**：需 nesting 鏈被觀察；C 類（span>76kb，668 區）為統計鏈接非單分子。
5. **⭐3 單樣本** HCC1395 single-pipeline；regional 非 genome-wide tree；分子共現≠single-cell。

## §7 可驗證產物

- `scripts/build_lineage_labels.py`（產生器；chr17 per-read 驗證 + genome-wide per-region）
- `data/lineage_labels_regions.tsv`（3,842 區：region/n_sSNV/span/tree_shape/haplotypes/situation_tier/vaf_order_ok/labels）
- 理論/門檻/驗證互引：`20260628_sSNV_linkage_threshold_decision_eps2_01.md`、`20260628_reconstruction_model_verification_01.md`、`_assets/.../PROVENANCE.md`。
