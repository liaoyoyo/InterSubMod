<!--
建立時間: 2026-02-09 14:40
目標: InterSubMod 專案全面審視 — 目標文件、執行架構、程式碼、輸出結果、結果分析的完整評估
處理範圍: 專案所有層面的改進建議與研究空間分析
關聯檔案:
  - docs/architecture/system_overview.md
  - docs/research/methylation_f1_optimization/2026/01/20260115_final_conclusion_report_01.md
  - CLAUDE.md
  - /big8_disk/liaoyoyo2001/Knowledge/README.md
-->

# InterSubMod 專案完整分析報告

**報告日期**: 2026-02-09
**分析範圍**: 目標文件、程式執行架構、程式碼品質、程式輸出結果、結果分析程式

---

## 目錄

1. [專案目標與定位](#1-專案目標與定位)
2. [程式執行架構分析](#2-程式執行架構分析)
3. [程式碼品質審視](#3-程式碼品質審視)
4. [程式輸出結果評估](#4-程式輸出結果評估)
5. [結果分析程式審視](#5-結果分析程式審視)
6. [可改進的方向](#6-可改進的方向)
7. [具研究空間的方向](#7-具研究空間的方向)
8. [整體評估與優先建議](#8-整體評估與優先建議)

---

> 校正註記（2026-02-09）：本版已對齊目前程式碼現況，修正輸出結構、PERMANOVA 預設行為、測試案例數與 Python 腳本規模等不一致敘述。

## 1. 專案目標與定位

### 1.1 核心目標

InterSubMod 旨在利用 ONT 長讀定序的**單分子甲基化模式 (Read-level Methylation Pattern)** 結合**體細胞變異 (Somatic SNV)** 與**單倍型 (Haplotype)** 資訊，達成三項分析目標：

| 目標 | 說明 | 目前狀態 |
|:-----|:-----|:---------|
| **Somatic SNV 驗證** | 透過 TP/FP 位點附近的甲基化分群模式差異，輔助判定變異真實性 | 已完成，但效果有限 (CramersV AUC = 0.52) |
| **等位特異性甲基化 (ASM)** | 比較 HP1 vs HP2 在局部 CpG 的穩定差異 | 已實作 LabelTest 多階段驗證 |
| **亞克隆結構分析** | 利用 per-SNV window 的 read 甲基化距離矩陣 → 聚類 → 演化樹 | 架構已建立，尚未有系統性的臨床驗證 |

### 1.2 目標文件的完整度評估

**優勢**：
- `docs/architecture/system_overview.md` 提供了詳盡的系統設計文件，包含完整 Mermaid 流程圖、9 個模組步驟 (S0-S8) 的規格說明
- 明確的全域參數預設值 (window_size_bp=1000, min_mapq=20, binary thresholds 等)
- 系統假設與限制有清楚標註 (僅支援 biallelic SNV、僅處理 ONT 長讀)

**不足之處**：
- **缺少明確的科學假說陳述**：文件描述了「做什麼」但缺乏「為什麼預期甲基化能區分 TP/FP」的理論依據
- **缺少量化的成功指標**：除了 F1 research 中的 baseline (0.8155)，沒有定義什麼數值代表「成功」
- **PMD/ChromHMM gating** 在架構文件中有描述，但實際程式碼中尚未看到完整實現 (CpGSite 中有 `is_pmd` 和 `is_repressive` 欄位，但 gating 邏輯未在 MatrixBuilder 中明確啟用)

### 1.3 參考知識庫的對應性

根據 Knowledge Base (`/big8_disk/liaoyoyo2001/Knowledge/`)：
- 樣本資訊完整覆蓋 HCC1395 (SEQC2 truth set) 和 COLO829 (NYGC truth set)
- 工具鏈描述 (ClairS, LongPhase, DeepSomatic) 與程式碼中的 VCF 處理邏輯一致
- **缺口**：Knowledge Base 中提到的 subsample purity series (t50_n00 ~ t10_n40) 和多樣本驗證，在目前的分析中僅使用了 HCC1395 單一樣本

---

## 2. 程式執行架構分析

### 2.1 整體架構

```
main.cpp → ArgParser → Config validation
    ↓
RegionProcessor (OpenMP 平行處理)
    ├── S0: SomaticSnv.load_from_vcf()     → SomaticSnvTable
    ├── S1: RegionOfInterest (POS ± window) → per-SNV 視窗
    ├── S2: BamReader.fetch_reads()         → 原始 BAM reads
    ├── S3: ReadParser + MethylationParser  → ReadInfo + MethylCall[]
    ├── S4: MatrixBuilder.build()           → Read × CpG 甲基化矩陣
    ├── S5: DistanceMatrix.compute()        → Read-Read 距離矩陣
    ├── S6: HierarchicalClustering          → 聚類樹
    ├── S7: SignificanceAnalyzer            → 統計顯著性
    └── S8: RegionWriter                    → 輸出所有結果
```

### 2.2 架構優勢

| 面向 | 評價 | 說明 |
|:-----|:-----|:-----|
| **模組化** | 優秀 | 每個步驟 (S0-S8) 對應獨立的 class，職責清晰 |
| **平行化** | 良好 | OpenMP per-region 平行，thread-local BamReader/FastaReader 避免鎖競爭 |
| **記憶體管理** | 良好 | RAII 模式、jemalloc 靜態連結、move semantics |
| **可測試性** | 良好 | 核心 library + GoogleTest（`tests/` 目前約 105 個 TEST 宏）與 `src/test/` phase 驅動測試 |
| **可擴展性** | 良好 | 距離度量、聚類方法皆為策略模式 (strategy pattern) |

### 2.3 架構風險點

| 風險 | 嚴重度 | 說明 |
|:-----|:------:|:-----|
| **每個 SNV 獨立分析** | 中 | 相鄰 SNV 的 window 可能高度重疊，但目前不共享 reads/CpG 資訊，重複計算造成效能浪費 |
| **O(n³) 聚類** | 中 | HierarchicalClustering 每次迭代做 O(n²) 最小值搜尋；n=500 reads 時約 1.25×10⁸ 比較 |
| **單一執行模式** | 低 | 整個流程從 VCF 到圖表是一條固定 pipeline，無法只重跑統計或只重繪圖 |
| **缺少 checkpoint** | 低 | 長時間全基因組分析無中斷恢復機制，失敗需從頭開始 |

### 2.4 資料流瓶頸分析

```
BAM I/O (HTSlib) → 甲基化解析 (O(seq_len)) → 矩陣建構 (O(reads × CpGs))
                                                    ↓
                                            距離計算 (O(n²k)) ← 主要瓶頸
                                                    ↓
                                            聚類 (O(n³)) ← 次要瓶頸
                                                    ↓
                                            統計檢驗（主要為 Fisher/LabelTest；若啟用 PERMANOVA 則近似 O(permutations × n²)）
```

對典型 region (100 reads, 50 CpGs)：
- 距離計算: 100² × 50 = 500K 運算 → 快速
- 聚類: 100³ = 10⁶ → 可接受
- 統計: Global/Local Fisher + LabelTest permutation 為主，預設流程可接受

對極端 region (1000 reads, 200 CpGs)：
- 距離計算: 10⁶ × 200 = 2×10⁸ → 需秒級
- 聚類: 10⁹ → 需數十秒，**可能成為瓶頸**
- 統計: 若額外啟用 PERMANOVA/Bootstrap，將成為**明確瓶頸**

---

## 3. 程式碼品質審視

### 3.1 程式碼規模

| 類別 | 檔案數 | 行數 (估算) |
|:-----|------:|----------:|
| Core C++ (src/core/) | 18 | ~6,438 |
| Headers (include/core/) | 24 | ~3,518 |
| I/O 模組 (src/io + include/io) | 4 | ~820 |
| Utils (src/utils + include/utils) | 6 | ~651 |
| main.cpp | 1 | 82 |
| **C++ 小計** | **53** | **~11,509** |
| Tests (`tests/` + `src/test/test_phase*`) | 13 | ~3,085 |
| Python tools (`tools/*.py`) | 7 | ~3,246 |
| Shell scripts (`scripts/*.sh`) | 4 | ~985 |

### 3.2 程式碼品質矩陣

| 指標 | 評分 | 說明 |
|:-----|:----:|:-----|
| **一致性** | 8/10 | 遵循 `.clang-format` (Google style, 120 col, 4-space indent) |
| **RAII 與資源安全** | 9/10 | BamReader 正確 destroy、move semantics、thread-local 設計 |
| **錯誤處理** | 7/10 | Filter reason bitmask 完善，但部分函式回傳空 vector 而非 std::optional |
| **文件化** | 7/10 | 架構文件完善，但程式碼內註解不均勻 |
| **測試覆蓋** | 6/10 | 有約 105 個 GoogleTest + phase 驅動測試，但缺少 MethylationParser、MatrixBuilder、LabelTest 的獨立單元測試 |

### 3.3 關鍵模組深入分析

#### 3.3.1 MethylationParser — 最複雜的模組

**關鍵風險**：
1. **反向股 MM/ML 解碼** (`MethylationParser.cpp:114-151`)
   - 反向股需從 SEQ 末端往回迭代，匹配 5'→3' 的 MM tag delta encoding 順序
   - 反向股的 CpG 驗證邏輯 (`ref[offset-1]='C' && ref[offset]='G'`) 正確但**非顯而易見**
   - 建議：加入更詳細的 edge case 單元測試（尤其是讀段起止恰好在 CpG 位點上的情境）

2. **MM tag 多修飾類型處理** (`MethylationParser.cpp:42-77`)
   - 以逗號計數來定位 `C+m` 修飾在 ML 陣列中的偏移量
   - **脆弱性**：假設 MM 字串格式良好；若含有 `C+h` (5hmC) 修飾，comma counting 可能出錯
   - 建議：考慮使用正則或更穩健的 parser

3. **缺少獨立單元測試**
   - 現有 `test_phase1_2.cpp` 是整合測試 (需實際 BAM)，無法獨立驗證解碼邏輯
   - 建議：建立 mock BAM record 測試 forward/reverse 各種 CIGAR pattern

#### 3.3.2 DistanceMatrix — 6 種度量的實現

| 度量 | 實現品質 | Edge Case 處理 | 備註 |
|:-----|:--------:|:-------------:|:-----|
| NHD | 良好 | common_sites < min_cov → -1.0 | Binary matrix 預設度量 |
| L1 | 良好 | 標準 | 連續概率 |
| L2 | 良好 | 標準 | 歐氏距離 |
| CORR | 良好 | zero-variance → 1.0 (max dist) | 需 ≥3 sites |
| BERNOULLI | 良好 | sum_weights < 1e-9 → -1.0 | 最複雜，confidence-weighted |
| JACCARD | 良好 | 可含/不含 unmethylated | Set-based |

**BERNOULLI 度量的特殊設計值得關注**：
```
w(p) = 2×|p - 0.5|        // confidence weight：概率越極端 (接近 0 或 1) 權重越大
delta = p(1-q) + (1-p)q   // expected disagreement：二項分布的不一致概率
dist = Σ(w_i × w_j × delta) / Σ(w_i × w_j)
```
這個設計理念是：只有高信心的甲基化位點才應對距離有大貢獻。這比簡單的 NHD 更適合 ONT 數據 (因 ML 值常有中間不確定值)。

#### 3.3.3 SignificanceAnalyzer — 多層統計框架

**5 層測試設計**：

```
SignificanceAnalyzer
├── GlobalTest:   Fisher-Freeman-Halton (RxC 列聯表 Monte Carlo)
├── LocalTest:    One-vs-Rest Fisher (每個 cluster 的個別顯著性)
├── StructureTest: PERMANOVA (基於距離的群組差異檢驗)
├── LabelTest:    多階段 HP/Allele 驗證
│   ├── Merged HP test (HP1 vs HP2 vs OTHER)
│   ├── Fine-grained HP test (HP1, HP2, HP1-1, HP2-1, HP3)
│   ├── Allele test (ALT vs REF)
│   └── Unassigned affinity test
└── Bootstrap:    叢集穩定性檢驗 (subsampling)
```

**品質觀察**：
- Fisher FFH 使用 Patefield 算法生成固定邊際的隨機表格，統計學上正確
- Monte Carlo 的 early stopping (99% CI) 在極端 p-value 時有效節省計算
- PERMANOVA 的實現正確 (pseudo-F = SS_Between/SS_Within 的 ratio)
- **現況補充**：主流程 `RegionProcessor` 目前預設 `enable_permanova = false`、`enable_dispersion = false`、`enable_bootstrap = false`，實際以 Global/Local + LabelTest 為主
- **潛在問題（啟用 PERMANOVA 時）**：PERMANOVA 要求完整距離矩陣 (無 NaN)，需確認 NaN strategy 與 filtering 在不同設定下行為一致

#### 3.3.4 HierarchicalClustering — O(n³) 的實現

實現了 4 種 linkage：UPGMA (預設)、WARD、SINGLE、COMPLETE。

- 使用 naive algorithm（每步 O(n²) 尋找最小距離）
- 可替換為 nearest-neighbor chain 或 MST-based 加速至 O(n² log n)
- 對本專案而言，n 通常 < 200 reads，**O(n³) 尚可接受**
- 若未來擴展到全基因組同時分析大量 reads，則需優化

### 3.4 測試覆蓋缺口

| 模組 | 有單元測試 | 缺口說明 |
|:-----|:---------:|:---------|
| Config | ✅ | 89 行，驗證配置合法性 |
| BamReader | ✅ | 139 行，測試區域擷取 |
| DistanceMatrix | ✅ | 464 行，各度量測試 |
| HierarchicalClustering | ✅ | 536 行，4 種 linkage |
| GlobalTest / LocalTest | ✅ | 366 行 |
| StructureTest / Bootstrap | ✅ | 348 行 |
| SignificanceAnalyzer | ✅ | 226 行 |
| SNV Loading | ✅ | 58 行 (較薄) |
| **MethylationParser** | ❌ | **最複雜模組無獨立測試** |
| **ReadParser** | ⚠️ | 有基本測試（配置與整合情境），但 **CIGAR 解析、ALT 判定邊界案例** 仍缺 |
| **MatrixBuilder** | ❌ | **矩陣建構邏輯無獨立測試** |
| **LabelTest** | ❌ | **多階段 HP 驗證無獨立測試** |
| **RegionProcessor** | ❌ | 整合測試依賴實際 BAM |

---

## 4. 程式輸出結果評估

### 4.1 輸出結構

每個 SNV region 產出完整的分析目錄：

```
output/{vcf_name}/{chr}/{chr_POS}/{chr_START_END}/
├── metadata.txt                    # 區域統計與品質資訊
├── reads/reads.tsv                 # Read 詳細資訊 (ID, QNAME, pos, MAPQ, HP, strand, allele)
├── methylation/
│   ├── cpg_sites.tsv              # CpG 座標
│   ├── methylation.csv            # Read × CpG 矩陣 (主矩陣)
│   ├── methylation_forward.csv    # 正向股
│   └── methylation_reverse.csv    # 反向股
├── clustering/
│   ├── tree.nwk                   # Newick 格式演化樹
│   └── linkage_matrix.csv         # 合併記錄
├── distance/
│   └── {METRIC}/
│       ├── matrix.csv
│       ├── matrix_forward.csv
│       └── matrix_reverse.csv
└── clustering/significance.json   # 區域顯著性結果（JSON）
```

Run-level（output root）另外產出：

```
output/
├── significance_summary.csv       # 全 run 彙總（非每個 region 各一份）
└── significance_statistics.txt
```

### 4.2 輸出品質評估

**優勢**：
- 輸出結構組織有系統，每個 region 自成完整分析單元
- run-level `significance_summary.csv` 包含豐富統計量 (GlobalP, CramersV, VerificationClass, DominantLabel, HP ratio 等)
- `metadata.txt` 包含品質評分系統 (100 分制)，分 High/Medium/Low 三級
- 支援 strand-aware 分析，forward/reverse 分別輸出

**不足之處**：
- **檔案量龐大**：30,000+ SNV × 每個 region 5-10 個檔案 = 150K+ 檔案，I/O 壓力大
- **跨批次彙總不足**：雖有單次 run 的全域 summary，但缺少跨樣本/跨 run 的自動合併報表

### 4.3 甲基化 F1 研究結果

根據 `docs/research/methylation_f1_optimization/2026/01/20260115_final_conclusion_report_01.md`：

**Baseline Performance** (ClairS 0.4.1 + longphase-s):

| 指標 | 值 |
|:-----|---:|
| TP | 30,490 |
| FP | 4,842 |
| FN | 8,957 |
| Precision | 0.8630 |
| Recall | 0.7729 |
| **F1** | **0.8155** |

**特徵區分能力**：

| 特徵 | AUC | 判定 |
|:-----|:---:|:-----|
| QUAL (VCF) | 0.9668 | 極佳 |
| AF (VCF) | 0.9235 | 極佳 |
| NumReads (甲基化) | 0.6303 | 有效 |
| GlobalP (甲基化) | 0.5614 | 弱 |
| CramersV (甲基化) | 0.5194 | 無效 |
| HeuristicScore (甲基化) | 0.4437 | 無效 |

**核心發現**：
1. VCF 原始特徵 (QUAL, AF) 的區分能力是甲基化特徵的 **2-3 倍**
2. 最佳策略 `QUAL < 0.8` 可將 F1 提升 **+2.70%** (0.8155 → 0.8424)
3. 甲基化顯著性的最大價值在**低品質邊界區域** (低 AF 時 NumReads AUC = 0.84)
4. 甲基化顯著的位點 94.8% 為 TP，可作為「高信心保留」依據

---

## 5. 結果分析程式審視

### 5.1 Python 分析工具概覽

| 腳本 | 用途 | 品質 |
|:-----|:-----|:-----|
| `compare_vcf_results.py` | TP/FP 比較分析 (ROC, AUC, distribution) | 良好，功能完整 |
| `plot_distance_heatmap.py` | 距離矩陣熱圖 + 樹狀圖 | 良好，支援多種標註 |
| `plot_cluster_heatmap.py` | 甲基化模式熱圖 + 聚類 | 良好 |
| `visualization_common.py` | 共用繪圖工具 | 良好，一致性高 |
| `analyze_depth_effect.py` | 覆蓋深度影響分析 | 基本功能 |
| `find_verification_candidates.py` | 候選位點篩選 | 實用 |
| `font_config.py` | CJK 字型配置 | 功能性工具 |

### 5.2 分析程式的優勢

- **`compare_vcf_results.py`** 整合了 sklearn 的 ROC/AUC 計算，支援多組比較
- 視覺化工具支援中文標註 (CJK font support)
- 熱圖同時顯示 HP tag、Allele、Tumor/Normal 標籤，資訊密度高
- `run_batch_vcf_analysis.sh` 支援批次比較多個 VCF

### 5.3 分析程式的不足

1. **缺少自動化彙總流程**
   - C++ 已輸出單次 run 全域 `significance_summary.csv`
   - 但跨 run 比較仍需手動指定多個路徑給 `compare_vcf_results.py`
   - 建議：在 batch script 中新增跨 run 自動 merge 與統一報表輸出

2. **缺少統計嚴謹的多重檢驗校正**
   - 30,000+ regions 各做獨立 Fisher 測試（PERMANOVA 目前預設關閉，啟用時同樣需要校正）
   - 沒有 Bonferroni 或 FDR (Benjamini-Hochberg) 校正
   - 可能大幅增加 false discovery

3. **缺少跨 region 的模式識別**
   - 每個 region 獨立分析，但亞克隆結構應跨越多個 SNV
   - 需要「多 region 一致性分析」來識別真正的亞克隆 vs 隨機噪聲

4. **缺少自動化 benchmark 流程**
   - F1 計算需要手動操作
   - 建議：自動化 TP/FP/FN 分類 + F1 計算 pipeline

---

## 6. 可改進的方向

### 6.1 演算法與效能改進

| 項目 | 優先級 | 預期收益 | 難度 |
|:-----|:------:|:---------|:----:|
| **A1. 重疊 window 的 read 共享** | 高 | 減少 BAM I/O 30-50%；相鄰 SNV 不必重複 fetch 同一批 reads | 中 |
| **A2. 跨批次 significance_summary 彙總** | 高 | 自動合併多個 run/樣本結果到單一 CSV | 低 |
| **A3. 多重檢驗校正 (FDR)** | 高 | 避免 30K+ independent tests 的 false discovery | 低 |
| **A4. MethylationParser 單元測試** | 高 | 確保 reverse strand 和 multi-modification 解碼正確 | 中 |
| **A5. HierarchicalClustering 加速** | 低 | O(n³)→O(n² log n)；目前 n 通常 <200，暫不急迫 | 中 |
| **A6. Checkpoint 機制** | 低 | 全基因組分析中斷恢復 | 中 |

### 6.2 程式碼品質改進

| 項目 | 說明 |
|:-----|:-----|
| **B1. MM tag parser 穩健化** | 用結構化 parser 取代 comma counting；支援 `C+m;C+h` 混合修飾 |
| **B2. ReadParser / MatrixBuilder 單元測試** | 建立 mock BAM record 測試各種 CIGAR pattern |
| **B3. NaN distance 策略統一** | 確保 -1.0 在所有下游模組 (聚類、PERMANOVA) 中被正確處理 |
| **B4. 錯誤處理改用 std::optional** | 區分「無資料」與「空結果」，避免空 vector 的歧義 |
| **B5. Python 分析工具模組化** | `compare_vcf_results.py` 雖僅 ~624 行，但職責偏多；可按載入/統計/繪圖拆分 |

### 6.3 輸出與分析改進

| 項目 | 說明 |
|:-----|:-----|
| **C1. 跨批次彙總 CSV** | 將多個 run 的 `significance_summary.csv` 自動合併為跨樣本總表 |
| **C2. 自動化 benchmark pipeline** | 整合 SEQC2/NYGC truth set 的 F1 計算到 CI |
| **C3. 機器學習分類器** | 在 `compare_vcf_results.py` 中加入 sklearn 分類器，利用組合特徵 |
| **C4. 互動式視覺化** | 考慮 Plotly/Dash 取代靜態 PNG，支援 zoom/filter |

---

## 7. 具研究空間的方向

### 7.1 甲基化特徵的「條件價值」— 最具前景的近期方向

F1 研究的核心發現揭示了一個重要模式：

> **甲基化特徵在整體上 AUC 低 (0.52)，但在低品質邊界區域 AUC 高達 0.84**

這意味著甲基化資訊不應作為獨立過濾器，而應作為 **VCF 品質不足時的輔助判斷依據**。

**研究方向**：
1. **分層過濾策略 (Tiered Filtering)**
   - Tier 1: QUAL ≥ 0.8 → 直接接受
   - Tier 2: QUAL < 0.8 但甲基化顯著 → 保留 (「救援」機制)
   - Tier 3: QUAL < 0.8 且甲基化不顯著 → 過濾
   - 預期效果：在不犧牲 recall 的前提下提升 precision

2. **甲基化信心分數 (Methylation Confidence Score)**
   - 結合 NumReads、GlobalP、CramersV 為單一分數
   - 用邏輯回歸或 gradient boosting 訓練
   - 在低 QUAL/AF 區域作為補充證據

### 7.2 跨 Region 亞克隆一致性分析 — 核心研究目標

目前每個 SNV window 獨立分析，但**真正的亞克隆應在多個 SNV 上展現一致的甲基化模式**。

**研究方向**：
1. **Multi-region 叢集比對**
   - 對同一 read 跨越多個 SNV window 的甲基化 pattern 做一致性檢驗
   - 若 Read A 在 SNV1 和 SNV2 都被分到 cluster1，強化亞克隆證據

2. **全基因組甲基化向量**
   - 不以 per-SNV window 為單位，而是建立 read 的全長甲基化 profile
   - 計算 read-level 的全基因組距離矩陣
   - 需要大幅改變架構，但科學上更有意義

3. **Phase block 層級的分析**
   - 利用 LongPhase 的 phase block (PS tag) 定義分析區域
   - 在同一 phase block 內的多個 SNV 共享 haplotype 資訊

### 7.3 5hmC 的整合分析

目前程式碼僅處理 5mC (`C+m`)，但 ONT 5kHz data (HCC1395) 同時提供 5hmC (`C+h`)。

**研究方向**：
- 5hmC 是主動去甲基化的中間產物 (TET 酶介導)
- 腫瘤中 5hmC 的分布可能更有亞克隆特異性
- 需修改 MethylationParser 支援同時解碼 `C+m` 和 `C+h`
- 建立 5mC + 5hmC 的雙通道距離矩陣

### 7.4 LOH 與 Copy Number 的整合

根據 Knowledge Base，LongPhase-TO 支援 LOH 偵測 (`--loh` flag)。

**研究方向**：
- LOH 區域的甲基化模式應更單一 (只剩一個 allele)
- HP ratio < 0.1 或 > 0.9 可能反映 LOH，目前程式碼已計算 HP ratio
- 可將 LOH 資訊整合到 significance scoring 中
- LOH 區域的 MultiHap variant 可能是真正的亞克隆標記

### 7.5 多樣本交叉驗證

目前所有分析僅在 HCC1395 (乳腺癌) 上進行。

**研究方向**：
1. **COLO829** (黑色素瘤): TMB 高，可測試高突變負荷下的甲基化分析
2. **H1437/H2009** (肺腺癌): Google 數據，僅有 5mCG
3. **Purity series**: t50_n00 ~ t10_n40，測試不同腫瘤純度下的表現
4. **跨平台**: ONT R10 vs PAO chemistry 的一致性

### 7.6 PMD/ChromHMM Gating 的實際啟用

架構文件中描述了 PMD (Partially Methylated Domain) 和 ChromHMM gating，程式碼 `CpGSite` 結構已有 `is_pmd` 和 `is_repressive` 欄位，但 **實際 gating 邏輯是否生效需要驗證**。

**研究方向**：
- PMD 區域的甲基化本底噪聲高，可能導致 false clustering
- 實作 BED annotation lookup，在 MatrixBuilder 中排除 PMD CpGs
- 比較 gating on/off 對 TP/FP 區分的影響

### 7.7 機器學習方法取代統計檢驗

目前以 Fisher 為主、PERMANOVA 為可選分析，有明確的統計學意義，但：

**研究方向**：
- 用 Random Forest / XGBoost 整合所有 15 個特徵 (F01-F15)
- 訓練 TP/FP 二元分類器
- 測試甲基化特徵是否在 ensemble 中貢獻 information gain
- **風險**：需防止過擬合，建議用 COLO829 做 held-out validation

---

## 8. 整體評估與優先建議

### 8.1 總體評分

| 面向 | 評分 | 說明 |
|:-----|:----:|:-----|
| 程式架構 | **8/10** | 模組化良好、平行化正確、職責清晰 |
| 程式碼品質 | **7/10** | RAII 優秀、部分模組缺測試、MM parser 有脆弱點 |
| 科學方法 | **7/10** | 統計框架完善、但缺多重校正和跨 region 分析 |
| 實驗驗證 | **6/10** | 僅 HCC1395 單樣本、缺跨樣本和跨平台驗證 |
| 分析工具 | **7/10** | 功能完善，但缺跨批次自動化彙總 |
| 文件品質 | **8/10** | 架構文件和研究報告詳盡 |

### 8.2 優先行動建議 (依投入產出比排序)

**第一優先 — 低成本高收益**：
1. ✅ 跨批次 significance_summary 自動彙總 (shell 腳本或 batch pipeline)
2. ✅ 多重檢驗校正 (FDR) 加入分析流程
3. ✅ MethylationParser 單元測試建立

**第二優先 — 中等投入中等收益**：
4. 🔧 分層過濾策略實驗 (Tiered Filtering with methylation rescue)
5. 🔧 COLO829 + purity series 交叉驗證
6. 🔧 MM tag parser 穩健化 (支援 C+m;C+h)

**第三優先 — 高投入高收益 (研究方向)**：
7. 🔬 跨 Region 亞克隆一致性分析
8. 🔬 5hmC 雙通道距離矩陣
9. 🔬 機器學習組合特徵分類器

### 8.3 關鍵結論

InterSubMod 已建立了一個**技術上成熟、架構合理**的甲基化分析平台。F1 研究清楚揭示了當前甲基化特徵的瓶頸 (整體 AUC 低)，但也發現了**甲基化在邊界案例中的獨特價值** (低 AF 時 AUC = 0.84)。

專案的下一步突破點不在於程式碼優化 (已足夠好)，而在於**科學方法論的升級**：
- 從 per-SNV 獨立分析 → 跨 Region 一致性分析
- 從單一樣本驗證 → 多樣本交叉驗證
- 從傳統統計 → 分層決策 + 機器學習輔助
- 從 5mC only → 5mC + 5hmC 雙通道

這些方向不僅能提升工具的實用性，更可能帶來值得發表的科學發現。

---

*本報告基於對專案完整程式碼、知識庫、架構文件、研究報告的全面審視所撰寫。*
