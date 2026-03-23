<!--
建立時間: 2026-03-23
問題類型: 統計方法 | 語意設計
影響 track: TO | 兩者
狀態: pending_decision
-->

# BERNOULLI 距離與 Weak/Subclone 語意定義 審查文件

## 問題 1：BERNOULLI 距離 p≈0.5 問題

**描述**：BERNOULLI 距離 = Σ |p_i(1-p_i)| × |x_i - y_i|，當某位點的甲基化率 p≈0.5（高雜訊位點）時，權重 p(1-p) 最大（≈0.25）。這使雜訊位點對距離貢獻過大。

**程式碼位置**：`src/core/DistanceMatrix.cpp`（BERNOULLI 計算函式）

**實際影響**：
- 目前距離矩陣無 NaN（0% NaN for Cases A/B/D）
- 但高雜訊位點的距離被過度放大，可能導致 UPGMA 聚類不穩定

**量化**：需從甲基化矩陣分析各樣本中 p≈0.5 位點的比例（待執行）。

---

## 問題 2：Weak vs Subclone 語意定義

**現況**（`src/core/SignificanceAnalyzer.cpp:331-334`）：

```
Subclone = label_sig=False AND cluster_significant=True
Weak = label_sig=True AND cluster_significant=False
```

**語意問題**：
- `Subclone` 的定義是「cluster 有結構，但 label（HP/Allele）無法解釋」
  - 這代表可能有第三維度的 subclone（非 HP、非 Allele 分組）
  - 或者是 PERMANOVA 停用導致的錯誤分類

- `Weak` 的定義是「label 有信號，但 cluster 結構不顯著」
  - 可能是信號很弱但確實存在
  - 也可能是 global test 閾值設太緊

**數據觀察**（HCC1395-5kHz TP）：
- Subclone TP: 1,373（4.8%）→ 全部 ClusterPassedGating=True，但 label_sig=False
- Weak TP: 10,141（35.6%）→ 其中 4.9% ClusterPassedGating=True，95.1% ClusterPassedGating=False

**問題**：Weak TP 中有 95.1% 是 ClusterPassedGating=False，代表沒有通過 global gating。但被標為 Weak（而非 Noise）是因為 `label_sig=True`（label-first 找到了信號）。

---

## 量化影響

### Subclone 分析（HCC1395-5kHz）

| 指標 | TP Subclone | FP Subclone |
|------|-------------|-------------|
| ClusterPassedGating | 100% True | 100% True |
| label_sig(LabelSig) | False | False |
| 比例 | 4.8% | ? |

Subclone 的存在暗示：cluster 有結構（passed_gate=True），但 HP/Allele 標籤無法解釋這個結構。

在 TO mode 中，這可能代表：
1. 真實的表觀遺傳亞克隆（真 TP 亞克隆）
2. 隨機聚類（global test 假陽性）

### Weak TP 的 label_sig 機制

Weak TP 有 label_sig=True。這來自哪裡？
- `label_sig = result.label_significant` = LabelTest 判斷有信號
- LabelTest 結合 HP PERMANOVA、LocalTest（HPMergedDelta, AlleleP）
- Weak TP 可能是 LocalTest（HPMergedDelta 或 AlleleSig）= True，但 ClusterPassedGating=False

---

## 修改選項

### 方案 A：不修改（接受現有定義）
- Weak/Subclone 定義有其邏輯合理性
- 目前 VerificationClass 未連接最終過濾，改動影響有限
- **F1 影響**：0

### 方案 B：重新審查 BERNOULLI 距離的 p≈0.5 問題
- 從實際甲基化矩陣量化 p≈0.5 位點比例
- 若比例 > 20%，考慮改用 NHD 或 truncated Bernoulli
- **修改位置**：`src/core/DistanceMatrix.cpp`
- **成本**：中
- **F1 影響**：不確定

### 方案 C：區分「Weak(label)」和「Weak(cluster)」
- 目前 Weak = label_sig AND NOT cluster_sig
- 可細分為：LabelAllele-Weak（allele 信號）vs LabelHP-Weak（HP 信號）
- 有助於理解哪種 Weak 更可靠
- **成本**：低（分析，無需改 C++）

---

## 建議

**方案 C**（純分析）→ 了解 Weak TP 的構成，再決定是否需要改設計。

---

## 用戶決策

**選擇**：[ ] A / [ ] B / [ ] C（分析 Weak 構成）
**日期**：
**理由**：
