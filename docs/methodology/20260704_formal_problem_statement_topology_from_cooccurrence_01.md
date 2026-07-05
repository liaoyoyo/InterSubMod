<!--
建立時間: 2026-07-04
更新: 2026-07-04 v3 (使用者定案 Model A) — unit-flip-tree 核心 + Boolean lattice + Model A(recurrence 允許) + 多重度擔子 + 四解法 + minimal-compatible-set; v2(共同設計) v1(初形式化) 併入
類型: 形式化問題敘述（純數學/演算法,application-agnostic）
狀態: formal spec v3 (application terms stripped)
用途: formal methods framing 候選 / 自足抽象問題定義
-->

# Rooted Unit-Flip Tree Reconstruction from Incomplete Binary Observations
### 從不完整二元觀測重建根式單位位元變化樹（Model A：允許 recurrence）
> 副標：*Enumerate all minimal Boolean-lattice trees compatible with partial 0/1 patterns; report identifiability rather than force a unique tree.*
> 不強迫唯一樹；輸出所有**最小相容樹**與可辨識性。**採 Model A（狀態樹，允許同位元在不同分支重複 0→1）**。

> 純數學/演算法；不使用任何應用/實體詞彙。

---

## §1 主問題（一句話）

> 給定一組含缺失值的二元字串 $O\subseteq\{0,1,?\}^k$，尋找嵌入 Boolean lattice 的 rooted directed tree：根 $=0^k$、每邊只允許一個位元 $0\to1$、每筆觀測至少相容於樹上一節點；輸出所有**最小相容樹**及可辨識性——**唯一 → determined；多個 → 輸出整個相容拓撲集合與不確定性**。

---

## §2 狀態空間與 unit-flip DAG

- 狀態空間 $X=\{0,1\}^k$；根 $r=0^k$。
- **合法有向邊** $x\to y$ ⟺ $y=x+e_j$（$y$ 比 $x$ 多恰一個 $1$，即 $\mathrm{Ham}(x,y)=1$ 且方向 $0\to1$）。
- 全體候選邊 = **Boolean lattice**（$k$-維超立方體只允許往 $1$ 增加的有向圖）。分層依 $\mathrm{popcount}(x)$（$1$ 的個數）：Level $\ell$ = $\{x:\mathrm{popcount}=\ell\}$。

## §3 觀測與相容性

- $O=\{p_1,\dots,p_n\}$，$p_i\in\{0,1,?\}^k$（$?$=未觀測）。
- **相容** $x\sim p$ ⟺ 對所有 $p[j]\neq?$ 有 $x[j]=p[j]$。（完整觀測強：$011$ 只相容 $011$；部分弱：$01?$ 相容 $\{010,011\}$。）

## §4 模型選擇（定案 Model A）

| | 每位元翻轉次數 | support set | 四型違反 |
|---|---|---|---|
| Model B（單次翻轉 / laminar-realizable）| 全樹 $\le1$ 次 | 必 laminar | → **conflict** |
| **Model A（本問題採用）** | 可在不同分支**重複** $0\to1$ | 無 laminar 約束 | → **recurrence-required（非 conflict）** |

**Model A 定義**：樹 $T=(V,E)$，$V\subseteq X$，$\mathrm{root}=0^k$，每邊為合法 unit-flip；**同一位元 $j$ 允許在多條邊上翻轉**（recurrence）。

🔴 **Model A 的擔子（consequence，必記）**：四型違反在 Model A 下是合法 recurrence，**不被 flag**。但**多重度假象**（座標 $x$ 有效 $m_f(x)>1$）與真 recurrence **在觀測上不可分**（皆表現為某位元翻兩次）。故 **Model A 必須有一個獨立多重度通道 $m_f(x)$（不得由 B 導出）** 來區分：
$$\text{recurrence-required 邊}\ \begin{cases}m>1\ \Rightarrow\ \text{多重度 artifact（棄）}\\ m\le1\ \Rightarrow\ \text{候選真 recurrence（留）}\end{cases}$$
無獨立 $m$ → Model A **over-build**（默默接受 artifact 樹）。此即 v2 §14「$m$ 由何獨立通道估」在 Model A 下**升為必需**。

## §5 相容樹與最小性

樹 $T$ **解釋** $O$ ⟺ $\forall p\in O,\ \exists x\in V,\ x\sim p$。
相容樹通常多（且可含多餘隱藏節點）→ 加最小化（**字典序**）：
$$\text{(1) 解釋所有 }O\ \succ\ \text{(2) min 隱藏(Steiner)節點}\ \succ\ \text{(3) min 邊數}\ \succ\ \text{(4) max 觀測支持}\ \succ\ \text{(5) 列出所有同分最佳}.$$
> 目標**不是**任一相容樹，是**所有最小相容樹**。唯一 → 拓撲可辨識；多個 → 問題欠定，輸出整個集合。

## §6 可辨識性輸出分類

| 類別 | 條件 |
|---|---|
| **determined** | 最小相容樹唯一（$0$ 隱藏節點且無同分）|
| **ambiguous** | 多棵同分最小樹 |
| **partial-order** | 資料只支持部分偏序（非完整樹）|
| **underpowered** | 覆蓋/長度不足以判斷 |
| **recurrence-required** | 存在位元須翻 $\ge2$ 次 → 送 $m$-通道拆 artifact vs 真 recurrence（§4）|
| **conflict** | （Model A 下罕見：即使允許 recurrence 仍無相容樹）|

## §7 工程解法（四路，皆納入）

**解法 1 — 小 $k$ 直接枚舉**：① 各部分觀測列 completion；② 建 Boolean lattice DAG；③ 選節點集 $V$ 使每觀測 $\ge1$ 節點 cover（set-cover）；④ 在 DAG 求 root→選定節點的**最小 arborescence**；⑤ 枚舉所有同分。本質 ≈ **directed Steiner tree / arborescence**（root $=0^k$；terminals=必解釋 completion；Steiner=隱藏狀態）。

**解法 2 — ILP**：變數 $z_x\in\{0,1\}$（狀態選入）、$e_{xy}$（邊選入）、$a_{px}$（觀測 $p$ 指派狀態 $x$）。約束：每 $p$ 至少一相容 $x$（$a_{px}\Rightarrow z_x$）；$e_{xy}\Rightarrow z_x\wedge z_y$；非根最多一 parent；每選節點 root-可達；每邊 unit-flip。目標 $\min$（隱藏節點 $+\lambda_1$ 邊數 $+\lambda_2$ 未支持節點 $+\lambda_3$ recurrence-with-$m{>}1$ 罰）。多最佳 → **no-good cut** 重解列全部。

**解法 3 — 分層 DAG-DP**：依 $\mathrm{popcount}$ 分層，每子節點只從上一層一 parent 來 → 「在分層 DAG 選 root-connected subtree，cover 所有 partial obs」的 DP。

**解法 4 — 保守 deterministic core**（主推穩固核心）：只接受 $0$-隱藏節點、laminar-或-recurrence-乾淨（$m\le1$）的樹；pairwise-only 串接標 weakly-determined；conflict/underpowered/multiplicity 分開標。抗 overfitting、可解釋。

---

## §8 輔助通道（沿用 v2；有界、非循環）

主結構由**通道 B**（上述 unit-flip-tree）定。另兩通道**只精修、不重定義**：
- **通道 A**（外部控制集 $\mathcal C$ 錨定的 2-著色 $\varphi$）：對互斥/sibling 候選過濾（異色→棄）。$\varphi$ 須外錨（$\varphi\perp$ B）否則循環。
- **通道 M**（實值輪廓的距離 $d$）：**只**做 (1) 無結構負篩、(2) 相容樹集**內**弱排序、(3) 佐證 B 已定邊。**禁**用 $d$ 把不可辨識強行全序化。
- **獨立多重度通道 $m_f(x)$**（§4 必需）：拆 recurrence-required 的 artifact vs 真 recurrence。

**循環性判準**：$O=g(\psi(D),D)$ 且對 $\psi$ 敏感 → 循環。核心保證：$\partial(\text{tree set})/\partial\varphi=\partial/\partial d=0$（樹集只由 B + $m$ 定）→ $\varphi,d$ 無條件不改樹集 → 非循環。

---

## §9 研究含意（Model A 對現行 pipeline 的翻案）

現行程式碼 = **Model B**（`classify()` 四型 = conflict；118 incompatible 一律歸多重度 artifact）。改採 **Model A** →
- 這 118 區從「conflict/artifact」重判為 **recurrence-required**；
- 送**獨立 $m$-通道**拆：$m>1$ → 確為多重度 artifact（與舊結論一致）；$m\le1$ → **候選真 recurrence 樹**（新結構，舊 Model B 漏掉）；
- **前置硬需求**：$m_f(x)$ 的獨立來源（如外部拷貝數/多重度估，非自 B）。無此則 Model A 無法安全落地。

---

## §10 輸出格式（每 region 一列）+ evidence levels

`region_id / n_events / n_fragments / n_multiway_fragments(≥3) / relation_class / minimal_tree_set / n_minimal_topologies / hidden_nodes / recurrence_edges / m_check_status / q_A / M_structure_p / final_status(robust|weak|unresolved|rejected)`

**evidence_level**：$L1$=pairwise only；$L2$=single-fragment $\ge2$ events；$L3$=$\ge3$ events；$L4$=$L3$+truth-backed+$m$-clean(-recurrence-verified)。

---

## §11 一句話版（Model A）

> 給定含缺失的二元字串，找嵌入 Boolean lattice 的 rooted tree（根 $0^k$、每邊一位元 $0\to1$、允許位元 recurrence、每觀測相容於某節點）；最小相容樹唯一 → determined；多個 → 輸出全相容集 + 不確定性；含 recurrence-required 邊 → 送獨立 $m$-通道拆 artifact vs 真 recurrence。

---

## §12 開放問題（優先序）

1. **獨立多重度通道 $m_f(x)$**（Model A 硬需求）——由何非-B 來源估。
2. **$q_A$**（B-事件片段未覆蓋 A 比例）——A-精修可用範圍。
3. **multiway fragment 比例**（覆蓋 $\ge2/3/4$ events）——證 long-fragment 有用。
4. **determined vs ambiguous 全域比例** + **recurrence-required 區的 $m$-拆分結果**（118 中有多少翻案成真 recurrence）。
5. $\mathbb P_{\mathcal T}$ 先驗；$d$ 排序互資訊 $I(d;\text{true-order})$；$\tau,W,C_{\min}$ 敏感度。

---

## §13 演算法定位（citation-verification 待驗；使用者提供）

> ⚠ 演算法層定位，非核心形式化；論文採用前跑 `/citation-verification`。
- laminar / 全三型相容性核心（Gusfield, perfect-phylogeny）——對應 Model B 的 conflict 判準（本問題採 Model A 放寬）。
- 同片段近端事件對加強推論（TreeClone, arXiv:1703.03853）——本框架擴到多位點片段標籤 + 可辨識性。
- 相容樹不確定性呈現（PhyloSub, arXiv:1210.3384）——對應 §5「輸出即測度」。
- 組合最佳化 / ILP 形式（PhISCS PMC6836735 / TUSV-ext）——對應 §7 解法 2。
- directed Steiner tree / arborescence——對應 §7 解法 1。
