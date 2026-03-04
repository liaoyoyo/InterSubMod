<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# Bernoulli 距離計算方法詳解

**建立日期**: 2025-12-24
**適用模組**: InterSubMod Distance Calculation

---

## 1. 核心概念 (Core Concept)

傳統的甲基化距離計算（如 Hamming Distance）通常先將甲基化機率二元化（>0.5 為 1，<0.5 為 0）後再計算差異。然而，Nanopore 測序產生的甲基化判讀是機率性的 (Probabilistic)，許多位點的機率值落在 0.5 附近的模糊地帶（Low Confidence / Ambiguous）。

**Bernoulli 距離方法**的設計初衷是：
1.  **直接使用原始機率**：不進行強制二元化，保留數據的不確定性資訊。
2.  **抑制雜訊**：降低低信心（模糊）位點對整體距離的影響。
3.  **統計意義**：將每個位點視為一個 Bernoulli 分布，計算兩條 Read 在該位點「取樣結果不一致」的期望機率。

---

## 2. 計算流程 (Calculation Process)

假設有兩條 Read $i$ 和 $j$，我們只針對它們**共同覆蓋 (Common Coverage)** 且**數據有效 (Non-NaN)** 的 CpG 位點集合 $K_{ij}$ 進行計算。

對於任意一個共同位點 $k$，令：
*   $p_{ik}$ 為 Read $i$ 在位點 $k$ 的甲基化機率 ($0 \le p \le 1$)
*   $p_{jk}$ 為 Read $j$ 在位點 $k$ 的甲基化機率

計算步驟分為三部分：

### 步驟一：計算期望不一致率 (Expected Disagreement)

我們定義差異量 $\delta_k$ 為「假設兩條 Read 的狀態分別服從參數為 $p_{ik}$ 與 $p_{jk}$ 的 Bernoulli 分布，其結果**不同**的機率」：

$$ \delta(p_{ik}, p_{jk}) = P(\text{Read } i \ne \text{Read } j) $$
$$ = P(i=1, j=0) + P(i=0, j=1) $$
$$ = p_{ik}(1 - p_{jk}) + (1 - p_{ik})p_{jk} $$

**特點**：
*   當 $p_i=1, p_j=0$ 時，$\delta = 1$ (完全不同)
*   當 $p_i=1, p_j=1$ 時，$\delta = 0$ (完全相同)
*   當 $p_i=0.5, p_j=0.5$ 時，$\delta = 0.5$ (最大不確定性)

### 步驟二：計算信心權重 (Confidence Weighting)

為了避免模糊位點（$p \approx 0.5$）主導距離計算，我們引入權重函數。單一位點的信心權重 $c(p)$ 定義為該機率偏離 0.5 的程度：

$$ c(p) = 2 \cdot |p - 0.5| $$

*   $p=0$ 或 $1$ (高信心) $\rightarrow$ $c(p) = 1$
*   $p=0.5$ (無資訊) $\rightarrow$ $c(p) = 0$

對於一對 Read 在位點 $k$ 的**聯合權重 $w_k$**，我們取兩者信心的乘積：

$$ w_k = c(p_{ik}) \cdot c(p_{jk}) $$

這意味著：**只有當兩條 Read 在該位點都具有一定信心時，該位點的差異才會有顯著貢獻。**

### 步驟三：加權平均 (Weighted Average)

最終的 Bernoulli 距離 $D_{ij}$ 是所有共同位點不一致率的加權平均：

$$ D_{Bernoulli}(i, j) = \frac{\sum_{k \in K_{ij}} w_k \cdot \delta_k}{\sum_{k \in K_{ij}} w_k} $$

---

## 3. 邊界情況處理 (Edge Cases)

在實作中，我們必須處理以下特殊情況：

1.  **覆蓋度不足 (Insufficient Coverage)**：
    若共同有效位點數 $|K_{ij}|$ 小於設定閾值（預設 `min_common_coverage` = 3），則視為資訊不足，回傳無效距離（-1 或 NaN）。

2.  **總權重為零 (Zero Total Weight)**：
    即使有共同覆蓋，若所有共同位點的機率都非常接近 0.5（導致 $\sum w_k \approx 0$），這表示雖然有重疊，但**沒有任何有效資訊**可以用來判斷相似度。
    *   **處理方式**：視同覆蓋度不足，回傳無效距離。

---

## 4. 範例演示 (Examples)

### 案例 A：明確的差異
*   Read 1: $p=0.9$ (高甲基化) $\rightarrow c(0.9) = 0.8$
*   Read 2: $p=0.1$ (低甲基化) $\rightarrow c(0.1) = 0.8$
*   聯合權重 $w = 0.64$
*   差異 $\delta = 0.9(0.9) + 0.1(0.1) = 0.82$
*   **結果**：高權重，大差異 $\rightarrow$ 距離顯著增加

### 案例 B：模糊位點 vs 明確位點
*   Read 1: $p=0.5$ (模糊) $\rightarrow c(0.5) = 0$
*   Read 2: $p=0.0$ (低甲基化) $\rightarrow c(0.0) = 1$
*   聯合權重 $w = 0 \times 1 = 0$
*   **結果**：權重為 0，此位點完全不影響最終距離計算（被忽略）。

### 案例 C：微小波動
*   Read 1: $p=0.01$ $\rightarrow c=0.98$
*   Read 2: $p=0.02$ $\rightarrow c=0.96$
*   聯合權重 $w \approx 0.94$
*   差異 $\delta \approx 0.03$
*   **結果**：高權重，極小差異 $\rightarrow$ 距離接近 0（正確識別為相似）。

---

## 5. 總結

Bernoulli 距離方法透過引入**機率模型**與**信心加權**，將 Nanopore 數據中的「模糊性」從雜訊轉化為可控的權重因子。相比於 L1 或 Hamming 距離，它能提供更穩健 (Robust) 的子克隆分群結果，特別是在測序深度較低或甲基化訊號較弱的區域。
