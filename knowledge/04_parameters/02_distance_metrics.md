---
id: ism-kb-04-parameters-distance-metrics
name: "Distance Metrics 規格"
description: "6 種距離度量數學定義：NHD / L1 / L2 / CORR / JACCARD / BERNOULLI；含公式、實作行號、應用場景、預設值。BERNOULLI 為 canonical benchmark 預設。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "distance metric formulas against src/core/DistanceMatrix.cpp:1-262"
related_ids:
  - ism-kb-04-parameters-index
  - ism-kb-04-parameters-cli-arguments
  - ism-kb-04-parameters-statistical-methods
  - ism-kb-05-data-formats-per-region-files
tags: [parameters, distance, metric, NHD, L1, L2, CORR, JACCARD, BERNOULLI]
canonical_paths: [04_parameters/02_distance_metrics.md]
alias_paths: []
---

# Distance Metrics 規格

- 一句結論：6 種距離度量，canonical 為 BERNOULLI（confidence-weighted），min_common_coverage 預設 3
- 適用對象：選擇距離度量時、客製新度量時
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ./build/bin/inter_sub_mod --distance-metric NHD L1 BERNOULLI -t ... -r ... -v ... -o out
  ```

---

## 六大度量總表

| 度量 | 公式 | 資料要求 | 實作 | 預設？ |
|------|------|---------|------|--------|
| **NHD** | `diff_count / common_valid_sites` | 二元化甲基化 | DistanceMatrix.cpp:31-55 | ✗ |
| **L1** | `mean(|p_i - p_j|)` | 原始機率 | DistanceMatrix.cpp:62-83 | ✗ |
| **L2** | `sqrt(mean((p_i-p_j)²))` | 原始機率 | DistanceMatrix.cpp:90-112 | ✗ |
| **CORR** | `(1 - pearson_corr) / 2` | 原始機率，≥3 重疊 | DistanceMatrix.cpp:121-175 | ✗ |
| **JACCARD** | `1 - |A∩B|/|A∪B|` | 甲基化/未甲基化集合 | DistanceMatrix.cpp:185-234 | ✗ |
| **BERNOULLI** | 見下方 | 原始機率 + confidence 加權 | DistanceMatrix.cpp:254-302 | **✓ canonical** |

---

## 詳細定義

### 1. NHD (Normalized Hamming Distance)
**公式**：
```
NHD(read_i, read_j) = sum_k[I(b_i_k ≠ b_j_k)] / n_common
```
其中 `b_i_k` 是 read_i 在 CpG 位點 k 的**二元**甲基化狀態（由 `--methyl-high` / `--methyl-low` 閾值決定）；`n_common` 是兩 reads 都有覆蓋的位點數。

**應用**：標準二元甲基化比較，最直觀
**實作**：`src/core/DistanceMatrix.cpp:31-55`

### 2. L1 (Manhattan)
**公式**：
```
L1(read_i, read_j) = mean_k |p_i_k - p_j_k|
```
`p` 為原始甲基化機率 [0, 1]，不二元化。

**應用**：全機率空間距離，對漸進變化敏感
**實作**：`src/core/DistanceMatrix.cpp:62-83`

### 3. L2 (Euclidean)
**公式**：
```
L2(read_i, read_j) = sqrt(mean_k (p_i_k - p_j_k)²)
```

**應用**：幾何距離，對極端差異敏感（平方放大）
**實作**：`src/core/DistanceMatrix.cpp:90-112`

**⚠ 警告**：O12 研究顯示，若用於 residualization on AF 會產生 collider bias（見 MEMORY feedback_L2_collider_bias）

### 4. CORR (Pearson Correlation Distance)
**公式**：
```
CORR(read_i, read_j) = (1 - pearson(p_i, p_j)) / 2
```
限 `n_common ≥ 3`（Pearson 相關性需至少 3 點）。

**應用**：模式相似度，scale-independent
**實作**：`src/core/DistanceMatrix.cpp:121-175`

### 5. JACCARD
**公式**（以甲基化/未甲基化集合）：
```
Jaccard = 1 - |A∩B| / |A∪B|
A = {k : b_i_k = methylated}
B = {k : b_j_k = methylated}
```

**應用**：位點集合相似度
**實作**：`src/core/DistanceMatrix.cpp:185-234`

### 6. BERNOULLI（canonical 預設）
**公式**：
```
BERNOULLI = sum_k (w_k × δ_k) / sum_k w_k

其中：
  w_k     = 2|p_i_k - 0.5| + 2|p_j_k - 0.5|    (confidence 權重)
  δ_k     = p_i_k × (1 - p_j_k) + (1 - p_i_k) × p_j_k   (disagreement 機率)
```

**直觀解讀**：
- `w_k`：兩個 reads 在位點 k 的信心度總和（越接近 0 或 1 越有信心）
- `δ_k`：兩個獨立 Bernoulli 分布不一致的機率
- **低信心位點 (p≈0.5) 被降權**，avoid 放大雜訊

**應用**：能量度量；降低低信心位點影響；benchmark canonical
**實作**：`src/core/DistanceMatrix.cpp:254-302`

---

## min_common_coverage (C_min)

**參數**：`--min-common-coverage`（預設 3）

**意義**：兩個 reads 需在至少 C_min 個 CpG 位點都有覆蓋，才計算距離。

**若不足**：由 `--nan-distance-strategy` 決定：
- `MAX_DIST`（預設）：填 `--max-distance-value`（預設 1.0）
- `SKIP`：跳過不納入距離矩陣

---

## 輸出檔案

每個度量產生獨立子目錄：
```
region_<N>/
├── distance_BERNOULLI/
│   ├── matrix.csv            # 完整
│   ├── matrix_forward.csv    # 僅前鏈 reads
│   └── matrix_reverse.csv    # 僅反鏈 reads
├── distance_NHD/             # 若同時指定 NHD
│   └── ...
└── ...
```

**格式**：對稱 N×N matrix；首行首列為 read_id。

---

## 選擇建議

| 目的 | 推薦度量 |
|------|---------|
| Canonical benchmark | **BERNOULLI** |
| 二元狀態比較 | NHD |
| 機率空間比較 | L1 > L2（避免 collider bias） |
| 模式相似度 | CORR |
| 集合式比較 | JACCARD |
| 多度量並行 | `--distance-metric NHD L1 BERNOULLI`（單次執行產生多 matrix） |

---

## 常見陷阱

- ❌ **誤把 L2 residualize 當獨立信號**：產生 collider bias（O12 已驗證）
- ❌ **NHD 不設 `--methyl-high`/`--methyl-low` 就用**：二元化閾值影響結果
- ❌ **CORR 在稀疏 CpG region 用**：小樣本 Pearson 不穩

---

## 相關

- CLI 參數：[01_cli_arguments.md](01_cli_arguments.md)
- 統計方法（用距離矩陣）：[03_statistical_methods.md](03_statistical_methods.md)
- Per-region 檔案：[../05_data_formats/03_per_region_files.md](../05_data_formats/03_per_region_files.md)
- 原始碼：`src/core/DistanceMatrix.cpp`（全檔 ~631 行，六個度量分布如上表）
