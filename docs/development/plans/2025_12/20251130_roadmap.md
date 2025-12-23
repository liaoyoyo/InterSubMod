* [x] **P0-4**: 設定基本 Unit Test 架構 (GoogleTest) 與 CI 整合 (VS Code Tasks)。

## Phase 1: 體細胞變異與區域定義 (S0, S1)

**目標**：處理 VCF 輸入，建立後續分析的錨點 (Anchor)。

* [ ] **P1-1**: 定義 `SomaticSnv` 與 `SomaticSnvTable` 資料結構。
* [ ] **P1-2**: 實作 VCF Parser (使用 htslib)，支援 Filter Config (`PASS` only, `min_qual`)。
* [ ] **P1-3**: 實作 `ChromIndex` 類別，處理 `chr_name` 到 `int ID` 的映射。
* [ ] **P1-4**: 實作 `RegionOfInterest` 生成邏輯 (Window Size 擴展)。
* [ ] **P1-5**: 單元測試：載入測試 VCF，驗證過濾邏輯與座標正確性。

## Phase 2: 讀段擷取與甲基化解析 (S2, S3)

**目標**：核心資料解析，這是最困難且最易出錯的部分。

* [ ] **P2-1**: 實作 `BAMFetcher`，使用 `sam_itr_multi` 批次抓取 Region Reads。
* [ ] **P2-2**: 實作 CIGAR Parser，正確處理 `Soft-clip`, `Insertion`, `Deletion`，映射 Read Base 到 Genome SNV。
* [ ] **P2-3**: 實作 `AltSupport` 判定邏輯 (含 `min_base_quality` 檢查)。
* [ ] **P2-4**: 實作 `MethylationParser`，解析 `MM`/`ML` tags，輸出 `(chr, pos, prob)` 列表。
* [ ] **P2-5**: **關鍵驗證**：人工構造含 Indel/Softclip 的 SAM record，驗證座標轉換是否精確無誤。

## Phase 3: 矩陣建構與距離計算 (S4, S5)

**目標**：將非結構化 Reads 轉為結構化矩陣並計算距離。

* [ ] **P3-1**: 實作 `MethylationMatrix` 類別，包含 `CpGSite` 列表維護與 Gating 邏輯 (PMD/Coverage)。
* [ ] **P3-2**: 實作矩陣填充邏輯 (Sparse to Dense/Bitset)。
* [ ] **P3-3**: 實作距離計算器 `DistanceCalculator`，支援 NHD (Normalized Hamming Distance)。
* [ ] **P3-4**: 優化：導入 Bitset (SIMD) 加速 Hamming Distance 計算。
* [ ] **P3-5**: 處理 Edge Cases：`NaN` distance (低重疊) 的處理策略 (`MAX_DIST`)。

## Phase 4: 聚類與統計檢定 (S6, S7)

**目標**：執行分群並進行生物學意義關聯。

* [ ] **P4-1**: 實作/串接 Hierarchical Clustering (UPGMA/Ward)。
* [ ] **P4-2**: 實作 `TreeCut` 邏輯，產生 Cluster Labels。
* [ ] **P4-3**: 實作 `ClusterStats`，統計每個 Cluster 的 HP/Tumor 比例。
* [ ] **P4-4**: 實作 Fisher Exact Test (需處理 log-factorial 避免溢位)。
* [ ] **P4-5**: 整合測試：輸入模擬的 "Perfect Tumor/Normal" 混合資料，驗證能否完美分群。

## Phase 5: 輸出與整合 (Output & Integration)

**目標**：產出使用者可讀報告與視覺化檔案。

* [ ] **P5-1**: 實作 JSONL Writer (Per-read detail)。
* [ ] **P5-2**: 實作 CSV Summary Writer (Per-region stats)。
* [ ] **P5-3**: 實作 Newick Tree String 生成。
* [ ] **P5-4**: 全流程整合測試 (End-to-End Test)。
* [ ] **P5-5**: 效能 Profiling 與優化 (Memory/CPU usage)。
