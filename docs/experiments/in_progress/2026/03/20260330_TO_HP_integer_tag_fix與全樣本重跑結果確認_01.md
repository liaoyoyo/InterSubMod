<!--
建立時間: 2026-03-30 06:30
目標: 記錄 LongPhase-TO HP integer tag 解析修正的完整驗證過程與全樣本影響評估
處理範圍: 全7個 TO 樣本，ReadParser.cpp HP tag 修正，before/after HP統計比較
關聯檔案:
  - src/core/ReadParser.cpp (修正位置: lines 121-141)
  - scripts/analysis/resume_to_pilot_from_tagged_outputs.sh
  - /tmp/run_hp_fix_rerun.sh (重跑腳本)
  - /tmp/full_hp_fix_analysis.py (分析腳本)
-->

# TO HP Integer Tag 修正與全樣本重跑結果確認

**日期**：2026-03-30
**修正位置**：`src/core/ReadParser.cpp:121-141`
**Binary 重建時間**：2026-03-30 03:19（163/163 tests passed）

---

## 1. 問題根源

### 1.1 Bug 描述

LongPhase-TO 輸出的 BAM 檔案使用 **HP integer format**（`HP:i:N`），而 LongPhase-S 使用 **HP string format**（`HP:Z:N`）。

舊的 `ReadParser.cpp` 在解析 HP 整數格式時，直接使用 `std::to_string(bam_aux2i(hp_aux))`，導致：

| LongPhase-TO 寫入 | 舊解析結果 | 正確結果 | 問題 |
|-------------------|-----------|---------|------|
| `HP:i:1`  | `"1"` | `"1"` | ✓ 正確 |
| `HP:i:2`  | `"2"` | `"2"` | ✓ 正確 |
| `HP:i:11` | `"11"` | `"1-1"` | ✗ **LabelTest 無法識別** |
| `HP:i:21` | `"21"` | `"2-1"` | ✗ **LabelTest 無法識別** |
| `HP:i:33` | `"33"` | `"3"` | ✗ **LabelTest 無法識別** |

被 LabelTest 無法識別的 reads（HP:i:11, HP:i:21, HP:i:33）全部被歸類為 `untracked`（既不算入 HP1Family 也不算入 HP0），在所有 HP 統計中完全消失。

### 1.2 修正方案

在 `src/core/ReadParser.cpp` 加入 switch 映射：

```cpp
} else if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I') {
    int hp_int = bam_aux2i(hp_aux);
    switch (hp_int) {
        case 1:  info.hp_tag = "1";   break;
        case 2:  info.hp_tag = "2";   break;
        case 11: info.hp_tag = "1-1"; break;
        case 21: info.hp_tag = "2-1"; break;
        case 33: info.hp_tag = "3";   break;
        default: info.hp_tag = std::to_string(hp_int); break;
    }
}
```

同時加入 `HP:Z:` string format 的處理，使 binary 對兩種格式均相容（longphase-s 和 longphase-to）。

---

## 2. 全樣本重跑結果

### 2.1 HP Read 回收（Table 1）

| 樣本 | Subset | Untrk%（修正前） | HP:i:11→HP1 恢復 | HP:i:21→HP2 恢復 | HP:i:33→HP3 恢復 | EffHP Δ% |
|-----|--------|----------------|-----------------|-----------------|-----------------|---------|
| H2009 | TP | **51.4%** | 3,294,573 | 3,000,121 | 112,063 | **+113%** |
| H2009 | FP | 44.5% | 272,392 | 250,661 | 10,040 | +84% |
| H1437 | TP | **68.8%** | 1,267,030 | 981,131 | 36,363 | **+271%** |
| H1437 | FP | 69.3% | 376,874 | 293,653 | 12,406 | +281% |
| HCC1937 | TP | 30.8% | 289,393 | 142,067 | 6,229 | +47% |
| HCC1937 | FP | 29.1% | 196,852 | 184,523 | 9,354 | +42% |
| HCC1954 | TP | 32.4% | 182,206 | 210,206 | 8,679 | +52% |
| HCC1954 | FP | 30.0% | 440,622 | 631,182 | 27,358 | +46% |
| COLO829 | TP | **61.3%** | 345,520 | 250,188 | 11,155 | **+251%** |
| COLO829 | FP | 60.5% | 184,505 | 136,652 | 5,731 | +229% |
| HCC1395D | TP | **49.0%** | 690,739 | 354,889 | 16,230 | **+114%** |
| HCC1395D | FP | **66.1%** | 365,118 | 197,632 | 8,192 | +227% |
| HCC1395 | TP | 38.4% | 494,705 | 275,268 | 9,572 | +70% |
| HCC1395 | FP | 39.4% | 186,542 | 129,166 | 5,753 | +70% |

> **關鍵發現**：修正後所有樣本 untracked reads = 0（完全解決）。修正前 29-69% 的 reads 是不可見的。

### 2.2 LOH Eligible Region 變化（Table 2，eff_hp ≥ 30）

| 樣本 | Subset | 修正前 LOH 數 | 修正前 % | 修正後 LOH 數 | 修正後 % | Δ% |
|-----|--------|-------------|---------|-------------|---------|-----|
| H2009 | TP | 81,023 | 64.5% | 122,858 | **97.7%** | +52% |
| H2009 | FP | 8,579 | 71.6% | 11,734 | **97.9%** | +37% |
| H1437 | TP | 11,144 | 24.5% | 43,440 | **95.5%** | +290% |
| H1437 | FP | 3,319 | 24.7% | 12,740 | **94.8%** | +284% |
| HCC1937 | TP | 10,948 | 86.7% | 12,290 | **97.4%** | +12% |
| HCC1937 | FP | 10,027 | 83.3% | 11,639 | **96.7%** | +16% |
| HCC1954 | TP | 11,944 | 70.0% | 15,876 | **93.0%** | +33% |
| HCC1954 | FP | 35,074 | 69.8% | 47,084 | **93.8%** | +34% |
| COLO829 | TP | 232 | **0.7%** | 11,475 | **34.7%** | **+4846%** |
| COLO829 | FP | 227 | 1.3% | 6,576 | 37.5% | +2797% |
| HCC1395D | TP | 13,521 | 46.9% | 26,857 | **93.1%** | +99% |
| HCC1395D | FP | 3,797 | 32.8% | 11,051 | **95.5%** | +191% |
| HCC1395 | TP | 16,665 | 58.5% | 26,276 | **92.2%** | +58% |
| HCC1395 | FP | 6,926 | 59.7% | 10,687 | **92.1%** | +54% |

> **COLO829 特別嚴重**：修正前只有 0.7% 的 regions LOH-eligible，LOH 分析幾乎完全失效。修正後 34.7%（雖仍低於其他樣本，可能是 COLO829 本身的生物特性）。

---

## 3. 各樣本嚴重程度評估

| 樣本 | EffHP Δ% (TP) | LOH Elig Δ% (TP) | 舊分析可信度 | 備注 |
|-----|--------------|----------------|------------|------|
| H1437 | +271% | +290% | **極低** | 69% reads 不可見，LOH 分析需全部重做 |
| COLO829 | +251% | +4846% | **極低** | LOH 分析完全失效（0.7% → 34.7%） |
| H2009 | +113% | +52% | **低** | 51% reads 不可見 |
| HCC1395D | +114% | +99% | **低** | 49% reads 不可見 |
| HCC1954 | +52% | +34% | **中** | 32% reads 不可見 |
| HCC1937 | +47% | +12% | **中** | 31% reads 不可見，LOH eligibility 接近天花板 |
| HCC1395 | +70% | +58% | **中-低** | 38% reads 不可見 |

---

## 4. 對過往研究的影響評估

### 4.1 Paired 資料（longphase-s, HP:Z: format）→ **不受影響**

LOH Round 1-4、Phase 1A Round 2-3 使用的 paired 樣本資料，HP tags 為 `HP:Z:` string format，舊 ReadParser 正確處理。這些分析結論仍然有效。

### 4.2 TO 資料（longphase-to, HP:i: format）→ **嚴重受影響**

所有以 TO significance_summary.csv 為基礎的 HP 統計（HP1FamilyN, HP2FamilyN, effective_hp, LOH eligibility）均是基於錯誤資料。

**受影響的分析**：
- TO LOH enrichment 分析（如果有）
- TO HP 分層分析
- TO LOH eligible region counts
- 任何使用 TO effective_hp 的閾值決策

**不受影響的分析**：
- F1 / precision / recall / TP / FP / FN counts（來自 VCF benchmark，不依賴 HP tags）
- ALT/REF read support（基於 CIGAR + base quality，不依賴 HP tags）
- 甲基化統計（基於 MM/ML tags，不依賴 HP tags）

### 4.3 LOH Round 1-4 結論有效性

LOH Rounds 1-4 主要使用 **paired** 資料（rescue_joined_features.tsv 來自 paired 樣本），HP tags 格式為 longphase-s（HP:Z:），**不受此 bug 影響**。

TO 側的 LOH 分析若有，需重新評估。

---

## 5. 下一步

1. **用新 TO 數據重新評估 TO LOH enrichment**：
   - 特別是 COLO829（從幾乎無 LOH-eligible 到有 34.7%）
   - H1437（LOH eligible 從 24.5% 到 95.5%）
2. **重新跑 TO LOH analysis rounds**（若有依賴 TO HP 統計的分析）
3. **更新 Phase 1A TO track**（如果 TO track 是進行中）
4. **備份確認**：所有樣本的舊數據已保存在 `intersubmod_tp_before_hp_fix/` 目錄

---

## 6. 執行紀錄

| 樣本 | 開始時間 | 完成時間 | 備注 |
|-----|---------|---------|------|
| H2009 | 03:25 | 04:04 | 第一個完成，用來驗證 fix |
| H1437 | 04:04 | 04:25 | - |
| HCC1937 | 04:25 | 04:44 | - |
| HCC1954 | 04:44 | 05:36 | - |
| COLO829 | 05:36 | 05:46 | - |
| HCC1395_DORADO | 05:46 | 06:01 | - |
| HCC1395 | 06:01 | 06:21 | 最後完成 |

總耗時：~3 小時（包含 fix 確認、binary 切換問題排除）
