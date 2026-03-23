<!--
建立時間: 2026-03-23
問題類型: C++ 方法 | 統計設計
影響 track: TO | 兩者
狀態: pending_decision
-->

# PERMANOVA 設計梳理與 Label-HP 失效率 審查文件

## 問題描述

ISM 的顯著分析包含兩個 PERMANOVA 路徑：

| 路徑 | 說明 | 狀態 |
|------|------|------|
| **cluster_structure** PERMANOVA | 用 cluster 分配標籤（k=2）檢驗距離矩陣結構 | **刻意停用**（`enable_permanova=false`） |
| **label_structure** HP PERMANOVA | 用 HP tag 分組（HP1/HP2）檢驗距離結構 | 70–97% 失效 |
| **label_structure** Allele PERMANOVA | 用 ALT/REF 分組檢驗距離結構 | 16–34% 失效 |

### 問題 1：Cluster PERMANOVA 刻意停用

**程式碼位置**：`src/core/RegionProcessor.cpp:1146`
```cpp
sig_config.enable_permanova = false;  // 刻意停用
```

`significance.json` 中的 `cluster_structure` 區塊永遠是：
```json
"permanova_valid": false, "permanova_f": 0.0, "permanova_p": 1.0
```

### 問題 2：Label HP PERMANOVA 失效率

失效原因（`src/core/SignificanceAnalyzer.cpp:187-216`）：
- 沒有 HP tag 的 reads 被排除（`raw_group_labels[i] < 0`）
- 排除後每組讀數 < `min_reads_per_group=3`
- 主要出現在：TO mode 中未相分的 region（缺乏 HP1/HP2 標籤）

---

## 量化影響

### Label HP PERMANOVA 失效率（3樣本統計）

| 樣本 | TP LabelHP invalid | FP LabelHP invalid | TP LabelAllele invalid | FP LabelAllele invalid |
|------|-------------------|-------------------|----------------------|----------------------|
| HCC1395-5kHz | 69.6% | 73.3% | 16.2% | 34.3% |
| COLO829 | 97.3% | 97.1% | 21.3% | 21.5% |
| H1437 | 90.4% | 95.2% | 19.0% | 24.1% |

**結論**：HP PERMANOVA 在 COLO829/H1437 中幾乎完全失效（>90%），主因是 TO mode 的 haplotagging 覆蓋不足。

### 各特徵 AUROC（3樣本合計，TP vs FP 二元分類）

| 特徵 | AUROC | 備註 |
|------|-------|------|
| LabelAllelePermanovaF | **0.569** | 最佳，但仍低 |
| LabelHPPermanovaF | 0.515 | 接近隨機 |
| ClusterCramersV | 0.511 | 接近隨機 |
| ClusterPassedGating | 0.529 | 微弱 |
| **Significant（最終判決）** | **0.511** | 幾乎隨機 |

### VerificationClass 分佈

| Class | TP 比例 | FP 比例 |
|-------|---------|---------|
| Noise | **46.6%** | **55.8%** |
| Weak | 28.5% | 24.7% |
| Strong | 22.1% | 16.2% |
| Subclone | 2.8% | 3.2% |

**關鍵觀察**：TP 中有 46.6% 被分為 Noise，說明 ISM 目前對這些位點無法產生有效判斷。TP 和 FP 的 class 分佈差異相當小（Strong TP=22.1% vs Strong FP=16.2%）。

---

## 根本問題分析

### 為何 Cluster PERMANOVA 停用？

可能原因（需確認）：
1. cluster-first 路徑的 cluster label 品質不穩定（k=1 仍可能通過，但 PERMANOVA 需 k≥2）
2. 計算成本（permanova 需 n_permutations 次 shuffle）
3. 設計決策：先觀察 label-first 效果

**影響**：`cluster_structure` 欄位完全沒有資訊貢獻，但還是在 significance.json 輸出（佔 I/O 資源但無用）。

### 為何 HP PERMANOVA 失效率這麼高？

TO mode 的讀數中，只有已通過 longphase haplotagging 的讀數才有 HP tag。在低 VAF 或高 noise 的 region，HP1/HP2 讀數可能：
- 數量不足（每組 < 3 reads）
- 分組後 dist matrix 過小

### 最終顯著判決（Significant）AUROC=0.511 的含義

ISM 的 `Significant=True` 標籤目前的 TP/FP 辨別力只略高於隨機猜測。這與跨樣本 F1 delta 接近 0 的觀察一致。

---

## 修改選項

### 方案 A：不修改（接受現狀，補文件）
- **理由**：cluster PERMANOVA 刻意停用可能有技術原因（待確認）；HP PERMANOVA 失效是資料問題而非邏輯錯誤
- **後續**：補充文件說明 `cluster_structure` 欄位的語意（目前 reserved/disabled）
- **F1 影響**：0（不改程式）

### 方案 B：重新評估 cluster PERMANOVA 的用途
- **修改位置**：`src/core/RegionProcessor.cpp:1146`
- **實作摘要**：
  1. 先確認停用的原因（git blame / 提問用戶）
  2. 若無技術障礙，嘗試開啟 `enable_permanova=true`
  3. 在 `validation` 測試中比較開啟前後的 F1 和計算時間
- **預估成本**：低（改一行）但需要完整驗證（時間成本中）
- **預估 F1 影響**：不確定，可能 ±0.001（cluster_permanova 目前沒有連接到最終分類決策）
- **風險**：cluster PERMANOVA 結果如何連接到 `VerificationClass`？需先確認邏輯

### 方案 C：改進 HP PERMANOVA 對低 HP 覆蓋的處理
- **修改位置**：`src/core/SignificanceAnalyzer.cpp:194` 或 `LabelTest.cpp`
- **實作摘要**：
  - 當 HP 覆蓋不足時，自動降階到 HP-merged（HPMergedDelta）作為替代
  - 或放寬 `min_reads_per_group` 到 2（從 3 降到 2）
- **預估成本**：中（需要測試各種 edge case）
- **預估 F1 影響**：可能 +0.001~+0.003（更多 label HP 有效，但精確度影響不確定）

### 方案 D：重新設計 Significant 判決邏輯（優先）
- **問題**：當前 `Significant=True` 的 AUROC 只有 0.511，代表判決邏輯本身有問題
- **修改位置**：`src/core/SignificanceAnalyzer.cpp:306-420`（VerificationClass 決策樹）
- **實作摘要**：
  - 梳理當前 VerificationClass 決策樹的每個分支條件
  - 找出 Strong FP 的主要進入路徑
  - 找出 Noise TP（被錯誤分為 Noise）的原因
- **預估成本**：高（需要完整方法學分析）
- **預估 F1 影響**：可能 +0.002~+0.010

---

## 建議優先順序

1. **立即確認**：cluster PERMANOVA 停用原因（低成本，釐清設計意圖）
2. **優先分析**：VerificationClass 決策樹（方案 D）— 這是 Significant AUROC=0.511 的根本原因
3. **中期**：HP PERMANOVA 低覆蓋處理（方案 C）

---

## 驗收標準

若選擇方案 B 或 C：
- [ ] test-quick 通過
- [ ] test-data F1 delta ≥ 0（至少不回退）
- [ ] 跨 3 個樣本（HCC1395-5kHz, COLO829, H1437）確認無負向影響

---

## 用戶決策

**選擇**：[ ] A / [ ] B / [ ] C / [ ] D / [ ] 其他
**日期**：
**理由**：
**下一步**：若選 B/C/D → 呼叫 `cpp-change` skill
