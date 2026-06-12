---
title: ISM 取樣設定依據 — window_size / C_min / read 過濾（數據支持的設定理由）
date: 2026-06-12
type: methodology / settings-rationale
status: in-progress
branch: chore/ism-review-param-governance-202606
data_sources:
  - 本次 ±1000bp 全量: 20260612_2201_sampling_threshold_audit/raw_summaries/tp_significance_summary.csv
  - canonical ±5000bp: output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_tp/significance_summary.csv
  - ±2000bp chr1 真實: 本輪 binary -w 2000 chr1 重跑（臨時已刪）
  - 既有: docs/methodology/20260324_H006_窗口大小TO調整分析_01.md
---

# ISM 取樣設定依據（為什麼這樣設）

> 目的：把 window_size / C_min / read 過濾的設定**用數據佐證理由**記錄下來，當作設定依據。

## 1. window_size = ±5000bp（為什麼不是 ±1000/±2000）

### 1.1 三 window 真實 CpG 數量（HCC1395 TP）

| window（±bp）| total bp | NumCpGs median | NumCpGs<30 | NumCpGs<10（距離隨機/無法分析）|
|---|---|---|---|---|
| ±1000 | 2kb | 15 | 83.6% | **24.9%** |
| ±2000 | 4kb | 31 | 47.5% | 1.3% |
| **±5000（canonical 預設）** | 10kb | 76 | 0.6% | **0.0%** |

線性關係 ≈ **7.5 CpG/kb**（人類基因組 CpG 密度 ~9.3/kb，扣除覆蓋/共同因素）。
- ±1000bp 真實 median 15 ↔ H006 理論估算 ~15 **完美驗證**。

### 1.2 設定依據
- 🔴 **±1000bp 有 25% 位點 CpG<10** → 距離矩陣/PERMANOVA 在這些位點不可靠或無法算 → 強制移除/被當噪訊。±5000bp 幾乎消除（0%）。
- **辨別力不因大 window 變差**（H006 實測：±5000 的 TP/FP PassedGating 比值 1.07 > ±1000 的 0.98）——「大 window 稀釋 SNV-local 訊號」的疑慮**在數據上未成立**，因 ONT long read（15-50kb）覆蓋 SNV 的 read 數不隨 window 變，只是 CpG 分析範圍變大。
- **判決（承襲 H006 + 本次真實驗證）**：**保持 ±5000bp**。±1000 大幅降覆蓋率（FN 7%→1.4%）且無辨別力提升。

### 1.3 外部做法 + 升級方向
- 多數外部工具（modkit/DSS/methylKit）**無 SNV-anchored window 概念**（per-CpG 率差 / DMR region），ISM 的 SNV±window 是軸 C 路線特有。
- 💡 **潛在升級**：借 MHB 的 **LD-r² greedy block** → 資料驅動定 window（用 haplotype block 邊界取代固定 ±5000bp）。見 method_comparison 03:90。

## 2. C_min（min_common_coverage）= 3

### 2.1 設定依據（避免距離太隨機）
兩條 read 的 NHD 距離基於「共同覆蓋的 CpG」。共同 CpG 太少 → 距離解析度低、隨機性大：

| 共同 CpG | NHD 可能值 | 隨機性 |
|---|---|---|
| 1 | {0, 1} | 最高（二元，1 個 CpG 差 = 完全不像）|
| 2 | {0, 0.5, 1} | 高 |
| **3** | {0, 0.33, 0.67, 1} | 開始穩定（≥4 級解析）|

- **C_min=3 = 「距離至少 4 個解析等級」的合理下限**，平衡距離穩定性 vs pair 覆蓋。
- 數據：±1000bp median 共同 CpG/pair=10，C_min=3 排除 ~10-12% pair（±5000bp 因 CpG 更多，排除更少）。
- ⚠ 注意 C_min 預設文檔不一致（Config.hpp=3 實跑 / 01_spec 寫 5）→ **統一改 3**。

### 2.2 距離不足 pair 的處理策略 — 方案抉擇（未定案）

被排除的 pair（共同 CpG<C_min）目前設 **MAX_DIST=1.0**。所有方案：

| 方案 | 做法 | 優點 | 缺點 | 工程 |
|---|---|---|---|---|
| MAX_DIST=1.0（現狀）| 沒交集當最遠 | 保留所有 read | 假距離污染分群 | 無 |
| SKIP | 沒交集設 NaN→剔除 read | 距離真實 | **剔除邊緣 read（漏觀察）** | flag |
| 加權距離 | 附「基於幾 CpG」信賴權重 | 保留+不污染 | 需改 clustering 算法 | 中 |
| soft/機率距離 | 機率非硬二值（cvlr 式）| 原則性最佳 | 大改 | 大 |
| imputation | 補缺失 CpG | 完整矩陣 | 引入假設偏差 | 中 |
| 動態 window(LD-block)| haplotype block 定窗 | 源頭減少不重疊 | 需 LD-r² | 中-大 |

**SKIP 的代價（實測前提）**：SKIP 不是漏 pair，而是 PERMANOVA `filter_reads_for_complete_matrix` 會**貪婪剔除 NaN 最多的 read**（覆蓋窗口邊緣、與多數 read 不重疊者）→ ±5000bp 16.3% pair 沒交集 → 剔除多少邊緣 read 需實測。

**抉擇維度**：保留 read vs 距離真實性 / 是否改 clustering 算法 / 工程量 / 對 subclone 分群的實際影響（必須實測）。

## 4. 🔴 不確定點（需用戶判斷 + 需實測，不假裝確定）

| # | 不確定點 | 現狀 | 定案需要 |
|---|---|---|---|
| U1 | MAX_DIST vs SKIP 哪個對 subclone 分群更好 | 只有方向推論 | **跑對比實驗**看實際 clustering 差異 |
| U2 | C_min=3 是否「最佳」（非僅合理）| 理論 ≥4 級解析合理，未優化 | C_min vs 分群品質曲線 |
| U3 | window=5000 vs 動態 LD-block | 5000 有據，LD-block 未測 | LD-block 實作 + 對比 |
| U4 | 上述參數對 subclone clustering 的敏感度 | 未直接測 | 參數 × 分群品質 grid |

> 寫論文/軟體框架時：window=5000 / C_min=3 理由 / 共同 CpG 影響分群 = **可直接寫（數據+理論支持）**；U1-U4 = **標為待驗證或需實測**，不可寫成定論。

## 3. read 過濾（沒達標準的 read）

### 3.1 過濾規則（ReadParser.cpp:26-69 + ReadAggregator）
| 原因 | 條件 |
|---|---|
| 低品質 | MAPQ < 20 |
| 太短 | read_len < 1000 |
| 無甲基資訊 | 缺 MM 或 ML tag |
| 不覆蓋 SNV | tumor read 在 SNV 非 ALT/REF（UNKNOWN）|
| 技術重複 | 同 qname dedup |
| secondary/supp/dup/unmapped | flag 過濾 |

### 3.2 記錄機制（已存在）
- binary 有 `--output-filtered-reads`（debug 模式輸出被過濾 read + 原因）。
- canonical 20260420 run **已用此 flag** → 可能已有 filtered-read 記錄（待 check）。
- **建議**：跑一次統計「各原因過濾掉多少 read」，對應 §參數實驗框架 A 的 read 層排除記錄（位點層 = RejectReason，read 層 = FilterReason）。

## L0 設定依據總結
| 設定 | 值 | 依據 |
|---|---|---|
| window_size | **±5000bp** | CpG 充足（median 76, <10=0%）消除強制移除；辨別力不降（H006 1.07>0.98）|
| C_min | **3** | 距離至少 4 級解析、避免隨機；排除 ~10% pair |
| read 過濾 | MAPQ20/len1000/MM+ML/SNV | 標準 ONT somatic-anchored；`--output-filtered-reads` 可記錄 |
