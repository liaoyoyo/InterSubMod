<!--
建立時間: 2026-03-22（補建文件：2026-03-24）
目標: 記錄 TO track 中 annotation_only 特徵的觀察結果（PairwiseMedianDist、CramersV、hp_assign_rate）
處理範圍: HCC1395 5kHz TO pilot（evidence_ledger H008, H009, H010）
關聯檔案:
  - research/autoresearch/evidence_ledger.jsonl (H008, H009, H010)
  - docs/methodology/20260322_TO_QS_GQ_Rescue規則採納分析_01.md (H011/H012 主規則)
-->

# TO Track：Annotation 特徵觀察表（H008/H009/H010）

**狀態：分析完成 — 三個特徵均為 annotation_only，不作為獨立 rescue 規則**

---

## 背景

在 H011（QS>=50）和 H012（GQ>=3）採納後，評估三個輔助特徵是否能提供獨立辨別力或與主規則正交組合：

1. **H008**：PairwiseMedianDist（PMD）作為附加 rescue 條件
2. **H009**：CramersV>0.3 作為高信心 annotation
3. **H010**：hp_assign_rate<0.5 作為低 haplotagging 品質標記

---

## H008：PairwiseMedianDist>=0.20

### 測試結果

| 指標 | 數值 |
|------|------|
| delta F1（PMD>=0.10 alone） | +0.008787（略優於 H011） |
| delta F1（GQ gate + PMD>=0.20） | +0.006662（遠遜於 H011） |
| TP rescued | 665 |
| FP introduced | 208 |
| H011 的正交 TP 貢獻 | 16/408（**3.9%**） |

### 關鍵觀察

- **PMD 與 QS 高度相關**（96.1% 重疊）：PMD>=0.20 的候選中 96.1% 也通過 QS>=50
- **PMD 獨立信號極弱**：只有 3.9% TP 是 PMD 能找到但 QS>=50 找不到的
- TP median PMD=0.210 vs FP median PMD=0.172 — 差距存在但不如 QS 精確
- 組合（QS>=50 OR PMD>=0.20）只額外增加 delta F1 +0.000209

### 判斷

`decision = annotation_only`

**可作為**：QS>=50 candidate 的「高甲基化聚類距離」支持 annotation（PMD 高 → 更可信的 cluster 結構）
**不作為**：獨立 rescue 規則（與 QS 高度冗餘）

---

## H009：CramersV>0.3

### 測試結果

| 指標 | 數值 |
|------|------|
| delta F1 | 0.0000（無直接 F1 貢獻） |
| TO CramersV>0.3 候選數 | 48（占全 TO candidates 4.5%） |
| TO CramersV>0.3 Precision | **85.4%** |
| Paired CramersV>0.3 Precision | 99.72% |
| H011（QS>=50）的子集？ | **是（48/48 = 100%）** |

### 關鍵觀察

- **CramersV 在 TO 中接近無效**：95.5% TO candidates 的 CramersV=0（中位數=0）
- **CramersV>0.3 全部包含於 QS>=50**：無正交信號
- **Precision 跨 track 不穩定**：paired=99.72% vs TO=85.4%（差 14.3%）
  → CramersV>0.3 在 TO 無法達到 paired 的高精確性，原因是 TO 中 haplotagging 把 reads 聚集到一個 HP，人工造成 CramersV 信號膨脹
- **TOtrack 中 CramersV 信號解讀**：TO 中 high CramersV 主要來自 HP 不均（非甲基化聚類差異），與 paired 的含義不同

### 判斷

`decision = annotation_only`

**可作為**：`CramersV>0.3 in TO` = "strong ISM methylation cluster signal"，但需注意其 Precision 僅 85.4%（vs paired 99.7%）
**不作為**：TO rescue 的精確度指標（跨 track 不一致）

---

## H010：hp_assign_rate<0.5（Low Phase Quality）

### 測試結果

| 指標 | 數值 |
|------|------|
| delta F1（直接 rescue） | 0.0000 |
| H011+hp>=0.5 filter F1 delta | H011 基礎上 **−0.000403**（降低） |
| hp<0.2 group FP rate | **46.1%**（高於整體 27.8%） |
| hp<0.5 group FP rate | 略高於整體 |
| TP median hp_assign_rate | 0.963 |
| FP median hp_assign_rate | 0.920 |

### 關鍵觀察

- **hp_assign_rate<0.2 是更有意義的閾值**（非 <0.5）：46.1% FP rate 相比整體 27.8%，顯著升高
- **hp<0.5 filter 降低 F1**：TP/FP 受 hp>=0.5 影響比例相似（4.8% vs 5.2%），無差別效果
- **hp<0.2 group 的解讀**：這些 candidates 的 haplotagging 幾乎失敗（大多數 reads 無 HP tag），甲基化信號不可靠
- **候選數小**：hp_assign_rate<0.2 group 只有 193 個，占 TO candidates 18%

### 判斷

`decision = annotation_only`

**可作為**：
```python
# 標記低 haplotagging 品質的 TO rescue candidate
df['low_phase_quality'] = (df['hp_assign_rate'] < 0.2)
# 這些 candidates 的 rescue FP rate 是正常 candidates 的 1.7 倍（46.1% vs 27.8%）
```

**不作為**：直接 rescue filter（TP 損失 ≈ FP 損失，F1 不升反降）

---

## 三特徵對比總表

| 特徵 | 決策 | 獨立 delta F1 | H011 子集？ | 主要用途 |
|------|------|-------------|------------|---------|
| PairwiseMedianDist>=0.20 | annotation_only | +0.008787（但 96% 與 H011 重疊） | 是（96.1% 重疊） | 聚類距離支持標記 |
| CramersV>0.3 | annotation_only | 0.0000 | 是（100% 子集） | 高甲基化信號標記（注意 TO 精度只有 85%） |
| hp_assign_rate<0.2 | annotation_only | 0.0000 | — | 低 haplotagging 品質 QC 標記 |

---

## Annotation 欄位建議使用方式

```python
# 加入 TO rescue candidate 的輔助 annotation 欄位
df['to_rescue_supporting'] = {
    'PMD_high': df['PairwiseMeanDist'] >= 0.20,     # 聚類距離支持
    'CV_high': df['CramersV'] > 0.3,                 # 甲基化信號強（注意 TO 精度 85%）
    'low_phase_quality': df['hp_assign_rate'] < 0.2, # haplotagging 失效（FP 風險高）
}
# 使用：在報告中標記 low_phase_quality=True 的 rescue candidate，提示下游慎用
```

---

## 判斷結論

所有三個特徵均已確認為 **annotation_only**：
- 無法提供正交於 H011/H012 的獨立 rescue 信號
- 作為 candidate 品質 annotation 有一定價值（尤其 hp_assign_rate<0.2 作為 FP 風險標記）
- 不建議作為 hard rescue filter（降低 F1 或無效果）
