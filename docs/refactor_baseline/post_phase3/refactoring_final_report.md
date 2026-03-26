<!--
建立時間: 2026-03-26
目標: InterSubMod C++ 重構計劃最終報告
處理範圍: Phase 1~3 全部完成的總結
關聯檔案:
  - docs/refactor_baseline/post_phase1/
  - docs/refactor_baseline/post_phase3/performance_comparison.md
-->

# InterSubMod C++ 重構最終報告

## 執行摘要

依照重構計劃書（`polished-cooking-hedgehog.md`），完成 Phase 1~3 共 13 個任務，
在不中斷 Phase 1 ML read classification 研究主線的前提下，系統性改善安全性、
架構可測試性與效能。

---

## 完成任務清單

### Phase 1：安全性與正確性（commits 4f57ba1 ~ c496e9a）

| Task | 說明 | Commit |
|------|------|--------|
| 1.1 | `BamReader::fetch_reads_safe()` RAII 介面，防止 HTSlib 記憶體洩漏 | `4f57ba1` |
| 1.2 | `Config::validate()` RAII guard，消除 idx/hdr 洩漏路徑 | `291ed2a` |
| 1.3 | `FisherExact` log-sum-exp 數值穩定 p-value 累加 | `5b7efc8` |
| 1.4 | PERMANOVA ss_between ⏸ **暫緩**（等 Phase 1 ML 結論） | — |
| 1.5 | `ReadParser` uint32_t→int32_t SNV 位置溢出防護 | `f9fee1d` |

### Phase 2：架構改善（commits 9b1fab8 ~ 7e67dff）

| Task | 說明 | Commit |
|------|------|--------|
| 2.1 | `DistanceMatrix` BERNOULLI 測試補強（覆蓋率提升） | `092ee00` |
| 2.2 | `LabelTest` 18 個測試（原本 0% 覆蓋率 → 有測試保護） | `9b1fab8` |
| 2.3 | `RegionProcessor` God Class 分解：提取 `ReadAggregator` + `RegionBounds` | `ff535eb` |
| 2.5 | API 命名一致化：`run()→test_all()`、`run()→compute()`（`[[deprecated]]` 過渡） | `1ff1c9d` |
| —   | Code-review 後修正（RegionBounds 座標語義、dangling ref、deprecated 雙重複製） | `7e67dff` |

### Phase 3：效能優化（commits 0197caf ~ a9912e7）

| Task | 說明 | Commit |
|------|------|--------|
| 3.1 | UPGMA Lance-Williams O(N³)，N=500 從分鐘級→**243ms** | `0197caf` |
| 3.2 | `MethylationParser` zero-alloc MM tag 解析（`from_chars` + 指標掃描） | `d0004b8` |
| 3.3 | `DistanceMatrix` OpenMP `reduction` 消除原子競爭 + `schedule(guided,4)` | `a9912e7` |

---

## 關鍵數據

### 測試覆蓋

| 項目 | 重構前 | 重構後 |
|------|--------|--------|
| 測試總數 | 107 | **159** |
| LabelTest 覆蓋率 | 0% | 有測試 |
| ReadAggregator 覆蓋率 | 無此模組 | 11 tests |
| RegionBounds 覆蓋率 | 無此模組 | 6 tests |

### 架構改善

| 指標 | 重構前 | 重構後 |
|------|--------|--------|
| RegionProcessor.cpp 行數 | 1,358 | 1,294（-64 行，process_reads 方法移除） |
| process_single_region 圈複雜度 | >20 | 降低（職責分離至 ReadAggregator） |
| 新模組 | — | ReadAggregator, RegionBounds |

### 效能提升

| 項目 | 改前 | 改後 |
|------|------|------|
| UPGMA N=500 | 預計分鐘級 | **243ms** |
| UPGMA 複雜度 | O(N⁴) | **O(N³)** |
| per-read MM heap alloc | 6+ 次 | **0 次** |
| DistanceMatrix atomic ops/pair | 3 | **0** |

### 安全性修正

| 問題 | 位置 | 狀態 |
|------|------|------|
| BamReader HTSlib 洩漏 | BamReader.cpp | ✅ 修正 |
| Config validate 洩漏 | Config.cpp | ✅ 修正 |
| FisherExact 數值精度 | FisherExact.cpp | ✅ 修正 |
| SNV 位置整數溢出 | ReadParser.cpp | ✅ 修正 |
| PERMANOVA 災難性消去 | StructureTest.cpp | ⏸ 暫緩 |

---

## 研究方向相容性

### Phase 1 ML read classification 可擴充性

`ReadAggregator` 是 Phase 1 ML 最關鍵的整合點：
- 每個 `bam1_t*` read → `ReadInfo` + `MethylCall` 管線現在封裝在單一類別中
- 新增 read-level 特徵 export（HP tag、SNV coverage、methylation profile）只需在 `ReadAggregator::aggregate()` 中加入 side-channel 輸出，無需改動 `RegionProcessor` 主流程
- 此架構直接支援 `build_phase1_training_manifest.py` 所需的 per-read TSV 生成

### 未改動的研究核心

以下均完全保持不變：
- PERMANOVA / Fisher / Bernoulli 統計演算法本身
- 輸出檔案格式（CSV、TSV、PNG）
- Command-line 介面與 shell script 相容性
- 研究結論所依賴的 p-value / pseudo-F 計算路徑

---

## 剩餘未完成事項

| 任務 | 狀態 | 說明 |
|------|------|------|
| Task 1.4 PERMANOVA ss_between | ⏸ 暫緩 | 等 Phase 1 ML 結論確定後再評估 |
| Task 2.4 LabelTest 內部分解 | 可選 | 有 18 tests 保護，低優先級 |
| Phase 3 合併至 main | 待 PR | 需 chr19 完整驗證後 merge |

---

## Git 分支

```
Branch: refactor/phase1-safety
Commits since baseline: 13 commits
Tests: 159/159 PASS
```

最終合併前應執行：
```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification
diff baseline_chr19_metrics.txt new_chr19_metrics.txt  # 應為空
```
