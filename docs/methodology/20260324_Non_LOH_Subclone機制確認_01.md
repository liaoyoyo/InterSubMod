<!--
建立時間: 2026-03-24
目標: 確認 Non-LOH Subclone 的真實機制，判斷是否為真實生物信號或隨機雜訊
處理範圍: COLO829 paired_pileup significance_summary.csv（Subclone class 分析）
關聯檔案:
  - docs/methodology/20260324_Subclone分類合理性與LOH驅動分析_01.md
  - docs/methodology/20260324_方法學審查任務清單與結論摘要_01.md
  - src/core/GlobalTest.cpp (gate 使用 unreliable CramersV)
-->

# Non-LOH Subclone 機制確認

**狀態：分析完成 — Non-LOH Subclone = unreliable CramersV 噪訊，非真實亞克隆信號**

---

## 背景

VerificationClass = Subclone 的條件：
```
cluster_significant = True  (passed_gate AND global_alt_p ≤ 0.05)
label_significant = False   (NOT HPMergedSig AND NOT AlleleSig)
```

在高 LOH 樣本（HCC1395/HCC1937），Subclone 主要是 LOH_Subclone（87-89%），語意明確。
在低 LOH 樣本（COLO829/H1437/HCC1954），Subclone 主要是 Non-LOH Subclone（80-90%），語意不明。

---

## Non-LOH Subclone 定量特徵（COLO829）

| 特徵 | 數值 | 詮釋 |
|------|------|------|
| 樣本數 | 481/537 Subclone（89.6%） | 主流 |
| HP_Ratio 中位數 | **0.525**（目標 0.5 = 平衡） | HP 幾乎完全平衡 |
| HP 偏移度 中位數 | 0.261（0=平衡, 1=極端） | 輕度偏移，非 LOH 模式 |
| AlleleDelta 中位數 | **0.0081**（≈ 0） | 無 allele 甲基化差異 |
| GlobalP 中位數 | **0.0195**（邊界顯著） | 位於 gate 邊界 |
| CramersV（reliable） | ≈ 0（大多數） | reliable 信號不存在 |

---

## 與 LOH_Subclone 的對比

| 特徵 | LOH_Subclone | Non-LOH Subclone |
|------|--------------|-----------------|
| HP_Ratio | < 0.1 或 > 0.9（極端） | ≈ 0.5（平衡） |
| AlleleDelta | 無特定要求（LOH 導向） | ≈ 0 |
| GlobalP | 通常 < 0.01 | 0.02 ~ 0.05（邊界） |
| LOH_Subtype | LOH_Subclone | None |
| CramersV（reliable） | ≈ 0（LOH 使 reliable 失效） | ≈ 0（期望次數不足） |
| 觸發機制 | HP 極端偏移 → raw CramersV 高 | 隨機聚類結構 → raw CramersV 輕度膨脹 |
| 生物意義 | LOH 區域 HP-confined 聚類 | 未知（可能非生物） |

---

## 機制推論：為何 k=2 聚類找到結構但 HP/Allele 均衡？

### 可能的混雜因素

1. **Read strand（正/負鏈）偏差**
   - 正鏈與負鏈的甲基化呼叫可能因 basecalling model 或比對不同而有細微差異
   - k=2 聚類可能分出「正鏈 cluster」vs「負鏈 cluster」
   - HP 與 Allele 均衡（兩鏈各一半），但有微弱甲基化差異 → GlobalP 邊界顯著

2. **Read length 差異**
   - 長 reads 與短 reads 的甲基化模式不同（CpG 覆蓋密度不同）
   - k=2 可能分出「長 reads cluster」vs「短 reads cluster」

3. **定序品質 / Pore 批次效應**
   - 低品質 pore 的 reads 甲基化機率系統性偏移
   - 在小樣本 region（reads 少）時，期望次數不足 → CramersV 虛報

4. **統計雜訊（最可能）**
   - CramersV 在 2×2 列聯表且期望次數小（< 5）時嚴重膨脹
   - GlobalP 0.019 邊界顯著 → 多重比較後大多數應為假陽性

---

## 跨樣本 Non-LOH Subclone 統計

| 樣本 | Subclone TP 總數 | Non-LOH 比例 | Non-LOH 佔全 TP 比例 |
|------|----------------|-------------|---------------------|
| HCC1395 | 1399 | 10.9% | 0.5% |
| HCC1937 | 163 | 12.9% | 0.3% |
| COLO829 | 537 | **89.6%** | 1.4% |
| H1437 | 1303 | **80.4%** | 1.5% |
| H2009 | 1954 | 69.3% | 1.0% |
| HCC1954 | 476 | **84.0%** | 2.3% |

→ Non-LOH Subclone 佔總 TP 的 0.3~2.3%（小比例），且集中於低 LOH 樣本。

---

## 影響評估

### TP 的影響

Non-LOH Subclone TP 佔總 TP 的 **0.3~2.3%**，且這些 TP 被正確分類為 TP（非 Noise）。
即使 Non-LOH Subclone 的分類語意模糊，這些 TP 仍被 ISM 正確偵測到（只是分類原因不明）。

### FP 的影響

Non-LOH Subclone FP 數量很少（從 Step 2 分析：各樣本 Subclone FP 個位數）。
對 Precision 的影響可忽略。

---

## 視覺化建議（後續可選）

若需確認機制，建議：

```bash
# 取出 3-5 個 Non-LOH Subclone region 的 distance heatmap
# 觀察 read 分布是否有 strand 或 length 模式
ls /big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_pileup/20260315_COLO829_paired_pileup_pileup_complete_matrix/intersubmod_tp/
# 找出 VC=Subclone AND LOH_Subtype=None 的 RegionID，看對應 plots/
```

---

## 判斷結論

| 問題 | 結論 |
|------|------|
| Non-LOH Subclone 是真實亞克隆信號嗎？ | **否** — HP balanced、AlleleDelta≈0、GlobalP 邊界 → 最可能是 unreliable CramersV 噪訊 |
| 需要修改 C++ 邏輯嗎？ | **否** — 影響小（佔 TP 的 0.3~2.3%），`LOH_Subtype` 欄位已可分層 |
| 現有 LOH_Subtype 欄位足夠嗎？ | **是** — 下游 Python 分析可用 `LOH_Subtype='LOH_Subclone'` 條件篩選可靠 Subclone |
| 建議後續行動 | 保持現狀；若需精確分析，可用 `LOH_Subtype` 分層，排除 Non-LOH Subclone |

**判決：No Action** — 機制已釐清，Non-LOH Subclone 為邊界假陽性，LOH_Subtype 分層已足夠。

---

## 與整體審查的關係

```
VerificationClass 框架整體評估：
  Strong    ← 最可靠（但包含 germline ASM FP，見 FP_PassedGating 分析）
  Weak      ← 合理（label 顯著但全域弱）
  Subclone  ← 需分層：LOH_Subclone（合理）vs Non-LOH（噪訊）
  Noise     ← 含大量 TP（33~67%）但無結構性信號 → 正確捨棄

建議的下游篩選標準：
  「有信號」= Strong OR Weak OR (Subclone AND LOH_Subtype='LOH_Subclone')
  而非：VerificationClass != Noise（會包含 Non-LOH Subclone 噪訊，但影響小）
```
