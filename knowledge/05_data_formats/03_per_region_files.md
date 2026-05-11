---
id: ism-kb-05-data-formats-per-region-files
name: "Per-Region 檔案結構"
description: "每個 SNV region 子目錄內的檔案：reads.tsv、methylation.csv、distance_<METRIC>/matrix.csv 的格式規範。"
status: active
last_verified: 2026-04-22
content_nature: reference
doc_type: reference
verified_scope: "per-region file structure against src/io/RegionWriter.cpp"
related_ids:
  - ism-kb-05-data-formats-index
  - ism-kb-05-data-formats-significance-summary-schema
  - ism-kb-04-parameters-distance-metrics
  - ism-kb-05-data-formats-merged-dataset-pitfalls
tags: [data-format, per-region, reads, methylation, distance-matrix]
canonical_paths: [05_data_formats/03_per_region_files.md]
alias_paths: []
---

# Per-Region 檔案結構

- 一句結論：每 region 子目錄含 reads.tsv、methylation/（3 個 CSV）、distance_<METRIC>/（3 個 matrix CSV）
- 適用對象：分析單一 region、驗證 ISM 計算、客製視覺化
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/*/TP/region_0/ | head
  ```

---

## 目錄樹

```
{output_dir}/{TP|FP}/region_{id}/
├── reads.tsv                            # read-level metadata
├── methylation/
│   ├── methylation.csv                  # 全 reads 的甲基化矩陣
│   ├── methylation_forward.csv          # 僅正股
│   └── methylation_reverse.csv          # 僅反股
└── distance_{METRIC}/                    # 每個指定度量一個子目錄
    ├── matrix.csv                       # 全 reads 的距離矩陣
    ├── matrix_forward.csv               # 僅正股
    └── matrix_reverse.csv               # 僅反股
```

**產出管理**：`src/io/RegionWriter.cpp:14-400`

---

## 1. reads.tsv

**內容**：該 region 每個通過過濾的 read 的 metadata

**欄位**：
| 欄位 | 類型 | 說明 |
|------|------|------|
| `read_id` | int | 從 0 開始的內部編號 |
| `read_name` | string | BAM 檔案內的 QNAME |
| `hp_tag` | int | 標準化 HP tag（0/1/2/3） |
| `hp_tag_raw` | string | 原始 HP tag（包含 HP:i:11/21/33） |
| `alt_support` | bool | 該 read 是否支持 ALT |
| `strand` | char | `+` / `-` |
| `align_start` | int | 對齊起點（0-based） |
| `align_end` | int | 對齊終點 |
| `mapq` | int | 映射品質 |
| `depth_at_snv` | int | SNV 位置深度 |
| `base_quality_at_snv` | int | SNV 位置的鹼基品質 |

**範例**：
```
read_id	read_name	hp_tag	hp_tag_raw	alt_support	strand	align_start	align_end	mapq	depth_at_snv	base_quality_at_snv
0	read_abc123	1	HP:i:1	0	+	12345	50000	60	45	35
1	read_def456	2	HP:i:2	1	-	12300	49800	58	45	33
```

---

## 2. methylation.csv

**內容**：reads × CpG 的甲基化矩陣

**格式**：
```
,cpg_12345,cpg_12378,...,cpg_49900
read_0,0.85,NaN,0.15,...
read_1,0.92,0.88,NaN,...
```

**欄位**：
- 第 1 欄：`read_id`（與 reads.tsv 對應）
- 後續欄：每個 CpG 位點（欄名 `cpg_<position>`）

**數值**：
- `[0.0, 1.0]`：原始 ML/(ML+MM) 機率
- `NaN`：該 read 在該 CpG 位點無覆蓋

**正反股分開檔**：
- `methylation_forward.csv`：只含正股 reads
- `methylation_reverse.csv`：只含反股 reads

**用途**：分析前後股甲基化不對稱性

---

## 3. distance_<METRIC>/matrix.csv

**內容**：read × read 成對距離矩陣（對稱）

**格式**：
```
,0,1,2,...,N-1
0,0.0,0.23,0.45,...
1,0.23,0.0,0.38,...
2,0.45,0.38,0.0,...
```

**欄位**：
- 第 1 行首列：read_id（0 ~ N-1）
- `(i, j)` 元素：`distance(read_i, read_j)`
- 對稱矩陣，對角線為 0

**METRIC** 取值：`NHD`, `L1`, `L2`, `CORR`, `JACCARD`, `BERNOULLI`（大寫）

**正反股分開檔**：
- `matrix_forward.csv`、`matrix_reverse.csv`

**數值範圍**：
| 度量 | 範圍 |
|------|------|
| NHD, JACCARD | [0, 1] |
| L1, L2, CORR, BERNOULLI | [0, 1]（canonical 設定下） |

詳見 [../04_parameters/02_distance_metrics.md](../04_parameters/02_distance_metrics.md)

---

## 載入範例

```python
import pandas as pd
import numpy as np

# reads.tsv
reads = pd.read_csv('TP/region_0/reads.tsv', sep='\t')

# methylation.csv
methyl = pd.read_csv('TP/region_0/methylation/methylation.csv', index_col=0)

# distance matrix
dist = pd.read_csv('TP/region_0/distance_BERNOULLI/matrix.csv', index_col=0)
assert dist.shape[0] == dist.shape[1] == len(reads)  # 三者 N 一致
```

---

## 檔案大小估算

| 檔案 | 每 region 大小 | 備註 |
|------|---------------|------|
| reads.tsv | 10-50 KB | 依 read 數 |
| methylation.csv | 100 KB - 5 MB | N_reads × N_CpGs |
| distance_matrix | O(N²) | 100 reads ≈ 40 KB；1000 reads ≈ 4 MB |

**全 7 樣本 × TP+FP**：約 **50-200 GB** 的 per-region 檔案。

---

## 常用對應關係

```
significance_summary.csv 的 RegionID 欄
  ↓ 對應
TP|FP/region_<RegionID>/ 目錄
  ↓ 進入
reads.tsv (NumReads 欄 = 行數)
methylation.csv (NumCpGs 欄 = 欄數 - 1)
distance_<METRIC>/matrix.csv (N × N)
```

---

## 相關

- significance_summary：[01_significance_summary_schema.md](01_significance_summary_schema.md)
- 距離度量：[../04_parameters/02_distance_metrics.md](../04_parameters/02_distance_metrics.md)
- 原始碼：`src/io/RegionWriter.cpp:14-400`
