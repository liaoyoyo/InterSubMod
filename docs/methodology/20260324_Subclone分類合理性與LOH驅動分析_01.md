<!--
建立時間: 2026-03-24
目標: 確認 Subclone class 分類合理性、LOH 驅動比例、與是否可更明確定義
處理範圍: 6 個 canonical paired_full 樣本的 significance_summary.csv 分析
關聯檔案:
  - src/core/SignificanceAnalyzer.cpp (lines 304-337, VerificationClass 賦值)
  - src/core/GlobalTest.cpp (lines 127-141, gate 邏輯)
  - docs/methodology/20260323_VerificationClass決策樹跨樣本量化_01.md
-->

# Subclone 分類合理性與 LOH 驅動機制分析

**狀態：分析完成 — Subclone 機制已釐清，建議分層使用 LOH_Subtype**

---

## 核心問題

1. Subclone class 分類合理嗎？
2. 是否主要是 LOH 驅動？各樣本差異為何？
3. 分類後是否可更明確定義與觀察？

---

## 機制釐清

### Subclone 的產生條件

```
cluster_significant = True   → passed_gate AND (global_alt_p ≤ 0.05 OR global_hp_p ≤ 0.05)
label_significant  = False   → NOT (HPMergedSig OR AlleleSig)
→ VerificationClass = Subclone
```

**gate 使用 raw（unreliable allowed）CramersV**，而 CSV 輸出只顯示 reliable CramersV：

| 樣本 | passed_gate TP | 其中 reliable CramersV ≈ 0 |
|------|---------------|---------------------------|
| HCC1395 | 31.7% (9421/29740) | **80.5%**（7587/9421） |
| COLO829 | 14.9% (5209/34971) | **91.0%**（4740/5209） |
| H1437 | 35.3% (23786/67467) | **66.3%**（15761/23786） |

→ 大多數 passed_gate 案例的「統計顯著性」來自 unreliable CramersV（期望次數不足的列聯表）。

---

## 兩類 Subclone 機制

### 類型 1：LOH_Subclone（生物學意義明確）

- **特徵**：`LOH_Subtype = 'LOH_Subclone'`，HP_Ratio ≈ 0 或 ≈ 1（極端 HP 偏移）
- **機制**：LOH 區域全部 reads 位於同一單倍體 → HP 分佈完全偏移 → raw CramersV 高 → gate passes
- **為何 label_significant=False**：同一單倍體內無跨 HP 對比可能 → HPMergedSig=False；純 LOH 甲基化差異未達 allele 顯著
- **意義**：真實的 LOH 驅動聚類結構，語意為「**LOH 區域中所有 reads 集中於一個 haplotype**」

### 類型 2：Non-LOH Subclone（可能為雜訊）

- **特徵**：`LOH_Subtype = 'None'`，HP_Ratio ≈ 0.5，AlleleDelta ≈ 0
- **機制**：k=2 聚類找到某種結構，但 HP 和 allele 均衡 → 真正原因未知（可能是 read length / strand / 其他混雜因素）
- **為何 passed_gate=True**：unreliable CramersV ≥ 0.1（low expected cell counts → CramersV inflated by chance）
- **意義**：大多數案例可能是**隨機聚類結構**，非生物學訊號

---

## 跨樣本 Subclone 統計

| 樣本 | SC-TP | SC-FP | LOH_Subclone% | Non-LOH% | 類型 |
|------|-------|-------|---------------|---------|------|
| HCC1395 | 1399 (4.7%) | 8 | **89.1%** | 10.9% | LOH 驅動（合理） |
| COLO829 | 526 (1.5%) | 49 | 10.1% | **89.9%** | 非 LOH（可能雜訊） |
| H1437 | 1303 (1.9%) | 0 | 19.6% | **80.4%** | 非 LOH（可能雜訊） |
| H2009 | 1954 (1.5%) | 0 | 30.7% | 69.3% | 混合 |
| HCC1937 | 163 (1.3%) | 1 | **87.1%** | 12.9% | LOH 驅動（合理） |
| HCC1954 | 475 (2.7%) | 1 | 16.0% | **84.0%** | 非 LOH（可能雜訊） |

**觀察**：
- HCC1395/HCC1937 高純度 LOH_Subclone → 這兩個樣本的 Subclone 分類語意清晰
- COLO829/H1437/HCC1954 以 Non-LOH 為主 → Subclone 語意模糊，可能是 false detection

---

## Non-LOH Subclone 特徵

| 樣本 | avg HP 偏移 | AlleleDelta 中位數 | GlobalP 中位數 |
|------|------------|------------------|----------------|
| COLO829 | 0.143 | 0.008 | 0.0195 |
| H1437 | 0.119 | 0.005 | 0.0035 |
| H2009 | 0.159 | 0.006 | 0.0190 |

→ HP 幾乎均衡，AlleleDelta ≈ 0，GlobalP 接近顯著邊界（0.02-0.05）。
→ 「顯著性」主要來自 unreliable CramersV 的統計雜訊，非真實甲基化訊號。

---

## 結論與建議

### 結論 1：Subclone 語意因樣本而異

- **HCC1395/HCC1937**：Subclone = LOH 區域的 HP-confined 聚類，語意合理
- **COLO829/H1437/HCC1954**：Subclone = unreliable CramersV 觸發的假陽性聚類，語意模糊

### 結論 2：LOH_Subtype 欄位已可分層

`LOH_Subtype = 'LOH_Subclone'` 已在 CSV 輸出中存在，可直接用於：
- `VC=Subclone AND LOH_Subtype=LOH_Subclone` → 明確的 LOH 驅動信號（建議保留）
- `VC=Subclone AND LOH_Subtype=None` → 疑似雜訊（下游分析時降權或排除）

### 結論 3：不建議修改 C++ VerificationClass 邏輯

原因：
1. LOH_Subtype 欄位已提供分層能力（無需 C++ 改動）
2. Non-LOH Subclone 在 TP/FP 辨別中影響有限（大多數樣本 FP Subclone 很少）
3. 修改 gate 邏輯可能影響 Strong/Weak 分類穩定性

### 建議：下游 Python 腳本中使用分層條件

```python
# 可靠的 Subclone 信號
df['Subclone_LOH'] = (df['VerificationClass'] == 'Subclone') & (df['LOH_Subtype'] == 'LOH_Subclone')

# 疑似雜訊 Subclone
df['Subclone_Noise'] = (df['VerificationClass'] == 'Subclone') & (df['LOH_Subtype'] == 'None')
```

---

## 與整體 VerificationClass 合理性的關係

### 通過顯著性判斷的邏輯合理性確認

| 類別 | 機制 | 是否合理 |
|------|------|---------|
| Strong | global_sig + label_sig → 強雙重確認 | ✓ 合理（但受 germline ASM 干擾） |
| Subclone-LOH | HP 極端偏移觸發 gate，無 label → LOH 結構 | ✓ 合理（已有 LOH_Subtype 標記） |
| Subclone-nonLOH | unreliable CramersV 觸發，無真實訊號 | △ 邊界（可能雜訊） |
| Weak | label_sig 但未過 gate → 弱整體結構 | △ 合理但需更多驗證 |
| Noise | 未過 gate → 無結構性訊號 | ✓ 合理（但包含弱信號 TP） |

→ **整體顯著性判斷框架合理**，問題在於 unreliable CramersV 造成 Subclone 過度寬鬆。
→ 如需改進：提高 gate 可靠性（例如要求 reliable CramersV ≥ 0.1，或設最低期望次數），但影響需評估。
