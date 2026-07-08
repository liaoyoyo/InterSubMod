---
title: subclone 重建演算法 — 形式化變數定義 + 流程 + 計算量（naive→收斂→實際）
date: 2026-07-06
sample: HCC1395 (frozen)
tier: L1（數字 grep-verified from sm_linkage_genomewide.json + topology_per_region.json，2026-07-06 canonical）
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json
build_branch: research/subclonal-reconstruction-202606
---

# subclone 重建演算法：形式化定義 + 流程 + 計算量

![L0→L5 資料流 + 計算量流程圖](_assets/20260706_algorithm_flow.svg)

> 流程圖：L0→L5 每階段的定義 + 計算量（樸素上界 → 三重收斂 → HCC1395 實際）+ 右側收斂機制 + 甲基不進計算的守則。

## 0. 變數定義（notation）

| 符號 | 定義 | 備註 |
|---|---|---|
| **G** | germline 雜合 SNP 集；G₁,G₂,… | **只用於 HP 定相**（外部 LongPhase-S），不進 somatic 組合 |
| **S** | somatic sSNV 宇宙；Sᵢ 第 i 個 | N = \|S\| = 全體 sSNV 數 |
| **W** | region（linkage 區）；Wⱼ 第 j 個 | = 相鄰 gap ≤ τ_R 的極大 sSNV 鏈 |
| **k** | 一個區內 sSNV 數 = \|S ∩ W\| | 區的「複雜度」 |
| **R** | reads；R = 區內 read 數；r = 單一 read | ⚠ 等位向量中的 R=REF、A=ALT（與 read 數 R 同字母，依上下文） |
| **ρᵢⱼ** | coread(Sᵢ,Sⱼ) = 同時覆蓋兩位點的 read 數 | powered ⟺ ρᵢⱼ ≥ ρ_min |
| **v** | genotype 向量 ∈ {R,A}^k（單 read 在 k 個位點的等位）| 例 RRA、AAA |
| **c** | cluster/population 數 = 區內 distinct v（≥m_min reads）| 「幾個細胞狀態」 |
| **H** | HP（germline 單倍型）標籤 H1/H2/H3 | 切 read 成 HP family |
| **T** | clone 樹；T₁,…,T_m 候選樹 | root = 全 R（germline） |

**參數（門檻）**：τ_R = 50,000 bp（same-read 窗）· ρ_min = 6（powered）· m_min = 3（population 支持）· κ = 8（genotype-向量 cap）· ε = 0.02（ONT 噪聲底線）。

---

## 1. 流程 L0→L5（每階段：定義 + 計算量 + HCC1395 實際）

### L0 — 宇宙 S
S = per-chrom (TP ∪ FP)。**計算量** O(N log N) 排序。**實際**：N = **35,332**（TP 30,490 + FP 4,842）。

### L1 — 定區 W（linkage 群）
以相鄰 sSNV gap > τ_R 切 S → {Wⱼ}。k_j = \|Wⱼ\|。**計算量** O(N) 掃描。
**實際**：k median **3** / p90 **6** / **max 150**；k>8 的區 = **171**。census 三桶（每 sSNV 一態）：linked **21,554** / underpowered **5,458** / isolated **8,320**（= N ✓）。

### L2 — pairwise 共現
對每對 (Sᵢ,Sⱼ)（dist ≤ τ_R）建 2×2 N(RR/RA/AR/AA) → classify（co_linked/nested/mutual_excl/independent）。
**計算量**：對數 = Σⱼ C(k_j,2)（理論）；但**每區單次 region-pileup** → 有效 O(R·k) 非 O(k²·pileup)。
**實際**：Σ C(k,2) = **79,858** 對；**n_pairs_recorded = 53,094**（powered/記錄）；coread ρ median **32** / p90 **91** / max **428**。

### L3 — multi-locus genotype 向量
k' = min(k, κ)（k>κ 取最密 κ 個）；對「全覆蓋 k' 位點」的 read 建 v；population = distinct v 且 ≥ m_min reads。
**計算量**：v ∈ {R,A}^{k'} → 上限 2^κ = **256** 種，但實際 ≤ R（read 數）。**這是 κ=8 cap 的來源**（k>8 時全覆蓋 read→0）。

### L4 — cluster 數上界 c（核心收斂）
**c ≤ min(R+1, k'+1)**：
- 樸素基因型空間 = 2^{k'}（指數）；
- **perfect-phylogeny**（每位點只突變一次，laminar）→ **c ≤ k'+1**；
- read 數限制 → c ≤ R+1；
- HP 分家 → 各 family 內再降。
**實際**：c median **1** / max **6**；**c > k+1 違反 = 0 / 6,288**（perfect-phylogeny 上界實證成立）。

### L5 — 建樹 T
root = 全 R（germline，如 RRR…）。從 pairwise 建：co_linked→併節點、nested→祖裔邊、mutual_excl→姊妹、independent→四配子違反。
- **候選樹數**（缺中間群順序歧義）：樸素 ≤ (k+1)!；但 **gained-pair pairwise 定序**把排列 prune 成 partial-order 的 linear extensions（CAP=24）→ block/ordered 收斂單一解。
- **k>κ（去 8-cap）**：`gt8` 用**全 sSNV pairwise** 建 position-level 樹（union-find 併 co_linked + 祖裔 DAG + 環偵測 + transitive reduction），O(k²)/區；節點數 ≤ k。

---

## 2. 計算量收斂全景（樸素 → 限制 → 實際）

| 階段 | 樸素上界 | 收斂機制 | HCC1395 實際 |
|---|---|---|---|
| genotype 空間 | 2^k | perfect-phylogeny + read | c ≤ k+1；median c=**1**、max **6**（0 違反）|
| 樹拓樸 | (k+1)! 排列 | pairwise 定序 → linear extensions | 缺中間群 69 區 → 多數 block/ordered、真等機率≈0 |
| pairwise 比較 | Σ C(k,2)=79,858 | region 單次 pileup | 有效 O(R·k)、53,094 對 |
| >8 區 | 全覆蓋 read→0（不可建 vector）| 全 pairwise position 樹 | 171 個 >8 區可 pairwise 建（乾淨子集）|

**核心洞見**：樸素搜尋是指數/階乘，但**真實資料被三重壓縮到極小**——perfect-phylogeny（c≤k+1，0 違反）、read 數（median 32 coread）、HP 分家。∴ 每區搜尋空間中位數只有 **c=1**（單一群；max 也僅 6）。

---

## 2b. 上限推導（約束下）+ 演算法處理

### (a) 樸素上限（無約束）
- 基因型空間　|{v}| = **2^k**（每 sSNV 兩態）
- 拓樸/順序　缺 g 個中間群 → **g!** 排列；整區樹排序 ≤ **(k+1)!**
- pairwise　全對 = **C(N,2)**（若不限窗）

### (b) 約束下的上限（算式）
令一區有 k 個 sSNV、R 條 read、覆蓋兩位點之 coread ρ、germline HP family {H}。

- **B1 — perfect-phylogeny（infinite-sites）**：每 sSNV 只突變一次 ⇒ 觀測基因型（樹節點）的突變集合構成 laminar family，每 sSNV 落唯一樹邊。⇒
  $$c \le k+1$$
  （推導：k 個突變事件 → 根 + 至多 k 條帶突變邊 → ≤ k+1 個相異節點；co_linked 併邊則更少。實測 **0 違反 / 6,288 區**。）

- **B2 — 覆蓋/read 支持**：每 population 需 ≥ m_min reads，總 R reads，且相異觀測向量 ≤ R。⇒
  $$c \le \min\!\big(R,\ \lfloor R/m_{\min}\rfloor\big)$$

- **B3 — HP 分家**：somatic read 依 germline HP 切；各 family 內獨立建樹。⇒
  $$c \le \sum_{H} (k_H + 1)$$

- **合併群數上限**：
  $$\boxed{\,c \le \min\!\Big(k+1,\ \lfloor R/m_{\min}\rfloor,\ \textstyle\sum_H (k_H{+}1)\Big)\,}$$

- **B4 — same-read 窗 τ_R**：L2 只配 dist ≤ τ_R 的對 ⇒ 每 sSNV 的 partner 數 ≤ w := |τ_R 窗內 sSNV|，pairwise 對數 ≤ **N·w / 2**（非 C(N,2)）。物理上 τ_R ≤ read span。

- **B5 — 拓樸候選數**：缺中間群 jump 邊帶 g 個 gained 突變，先以 gained-pair pairwise 建 partial order（co_linked 併群、nested 定序、4-gamete flag）；候選 = 該 poset 的 **linear extensions** e(poset) ≤ g!；整區 T ≤ Π_edges e(poset)，演算法截於 **CAP=24**。多數區收斂到 **1**（真 order-ambiguous ≈ 0）。

### (c) 演算法處理：以「觀察 c + 從 pairwise 建樹」取代枚舉指數空間
核心：**不列舉 2^k 基因型、不列舉 (k+1)! 樹**；c 由資料「數」出、樹由 pairwise「建」出。

```text
ReconstructRegion(W):                       # W: sSNV S_1..S_k, reads R
  # L2 pairwise（B4：只掃窗內對，單次 region pileup）
  for (S_i,S_j) with |pos_i-pos_j| ≤ τ_R:
      2×2 ← co-reads;  rel_ij ← classify(RR,RA,AR,AA)   # O(R·k) 總計

  # L3–L4 觀察 c（B1+B2：c 是量出來的、非枚舉）
  P ← multiset of genotype vectors over k'=min(k,κ) full-cover reads
  c ← |{ v : P[v] ≥ m_min }|                # 直接 ≤ min(k'+1, ⌊R/m_min⌋)

  # L5 從 pairwise 建樹（perfect-phylogeny，O(k²)，非 2^k）
  UF ← ∅;  anc ← ∅
  for each powered (i,j):
      co_linked → UF.union(i,j)             # 併節點
      nested    → anc.add(祖→裔)  (ε 去噪)   # 定祖裔
      4-gamete  → conflict++
  T ← transitive_reduction( DAG(UF, anc) )  # O(k²) 邊 + 約簡
  if cycle(T): return CONFLICT              # 不成單樹

  # 缺中間群：只列 linear extensions（B5），非 g!
  cand ← linear_extensions(partial_order)[:CAP]
  return (T, cand, c)
```

**演算法如何「吃掉」每個上限**：
| 上限 | 樸素 | 演算法手段 | 有效複雜度 |
|---|---|---|---|
| 基因型空間 | 2^k | 掃 read 數 distinct v（不列舉）| O(R·k′) |
| 群數 c | — | 直接量測，受 B1/B2 夾 | c ≤ k+1（實測 median 1）|
| 建樹 | (k+1)! | union-find + 祖裔 DAG + transitive reduction | O(k²) |
| 順序候選 | g! | poset linear extensions + CAP | 多數 = 1 |
| pairwise | C(N,2) | τ_R 窗 + 單次 pileup | O(R·k)/區 |
| k>κ（去 8-cap）| vector 不可建 | 全 pairwise position 樹（gt8）| O(k²)/區 |

⇒ **每區實際計算 = O(R·k + k²)**，k median 3（max 150）、c median 1 ⇒ 全基因組線性於 (Σread × sSNV)，非指數。

---

## 3. HCC1395 實際資料量（grep-verified，2026-07-06 canonical）

| 量 | 值 |
|---|---|
| N = \|S\| | 35,332（TP 30,490 / FP 4,842）|
| census 三桶 | linked 21,554 / underpowered 5,458 / isolated 8,320 |
| n_pairs 記錄 | 53,094 |
| Σ C(k,2) | 79,858 |
| N_W（有向量區）| 6,288 |
| k | median 3 / p90 6 / max 150 / >8: 171 |
| c | median 1 / max 6 / **c>k+1: 0** |
| coread ρ | median 32 / p90 91 / max 428 |
| determinacy | A_determined 1,764 · B_pairwise 943 · C_underdetermined 544 · A_ambiguous 69 · E_subcube_recovered 2,403 · incompatible 18 · recurrence(cand 9/artifact 11/LOH 50) · other 477 |

> ⚠ determinacy 分類含其他 session 2026-07-04/05 的 incompatible 重分類（E_subcube_recovered / recurrence_*）；N_W 由 3,885→6,288。此表為當前 canonical。

---

## 4. 限制（physical / algorithmic bounds）
- **τ_R = 50kb**：物理上限 = 單 ONT read aligned-span（p95~34kb、max~76kb）；跨此需 Tier-PS（germline-only 統計相位，非克隆連鎖，未做）。
- **κ = 8 vector cap**：k>8 時全覆蓋 read→0（工程限制，非物理）→ 由 L5 全 pairwise 建樹繞過（但那是 B_pairwise 組合、非 A_determined 單分子）。
- **c ≤ k+1**：perfect-phylogeny 理論上界（實證 0 違反）。
- **甲基不進 L0-L5 任何計算**（🔴 循環；只事後 bounded-auxiliary 註記）。

---

## 5. provenance
數字來源：`sm_linkage_genomewide.json`（universe/pairs/census）+ `topology_per_region.json`（regions/k/c/determinacy）。演算法：`sm_linkage_genomewide.py`（L0-L2）·`sm_multilocus_combinations.py`（L3）·`sm_region_integration.py`+`topology_analysis.py`（L4-L5）·`enumerate_candidate_trees.py`（定序）·`gt8_pairwise_tree.py`（>8）。
