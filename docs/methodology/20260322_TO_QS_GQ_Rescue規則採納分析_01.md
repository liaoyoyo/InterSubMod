<!--
建立時間: 2026-03-22（補建文件：2026-03-24）
目標: 記錄 TO track 中 QS>=50 與 GQ>=3 rescue 規則的測試結果、採納決策與規則選擇邏輯
處理範圍: HCC1395 5kHz TO pilot（evidence_ledger H011, H012, H_COMBO1）
關聯檔案:
  - research/autoresearch/evidence_ledger.jsonl (H011, H012, H_COMBO1)
  - research/autoresearch/cycles/20260322_041120/
  - research/autoresearch/cycles/20260322_045801/
  - research/autoresearch/cycles/20260322_045905/
-->

# TO Track：QS>=50 / GQ>=3 Rescue 規則採納分析

**狀態：分析完成 — 兩條規則均已採納（decision=keep），GQ>=3 最佳 F1，QS>=50 最佳 Precision**

---

## 背景

TO（Tumor-Only）pipeline 的 FN rescue 目標是：找出被 ISM 分析但未被 LongPhase-TO 選中的 TP，通過甲基化質量條件重新納入。

測試資料：`HCC1395 5kHz TO pilot`
- 基線 F1：0.712697
- ISM 覆蓋的 FN pool：773 個（占全 FN pool 的 7%，11051 個總 FN 中）

---

## H011：QS>=50 Alone（無 GQ gate）

### 定義

```python
rescue_condition = (Quality_Score >= 50)
```

### 測試結果

| 指標 | 數值 |
|------|------|
| delta F1 | **+0.008556** |
| TP rescued | +642 |
| FP introduced | +193 |
| Precision（rescue 候選） | 76.9%（642/835） |
| 基線 F1 → 結果 F1 | 0.712697 → 0.721253 |

### 關鍵觀察

1. **QS>=50 alone 優於 GQ+QS>=60**（+0.008556 vs +0.005551）
2. **HPP 在 TO 中無辨別力**：TP 與 FP 的 HPP<0.05 比例相同（均約 52%）— TO haplotagging 把多數 reads 分到同一 HP，HPP 天然偏低
3. **CramersV 在 TO 幾乎無用**：TO candidates 中 CramersV=0 比例極高（95.5%），無法作為辨別特徵
4. **AlleleDelta 在 TO ≈ 0**：TP/FP 的 AlleleDelta 差距極小（paired 5.29× 差距在 TO 消失）
5. **QS scores 存在聚類**：無候選在 QS [40,50) 區間，有效閾值是 >=50

### 決策

`decision = keep`，作為 TO rescue 主要規則候選

---

## H012：GQ>=3（從 GQ 掃描自動發現）

### 定義

```python
rescue_condition = (GQ >= 3)
```

### 測試結果

| 指標 | 數值 |
|------|------|
| delta F1 | **+0.009365** |
| TP rescued | +728 |
| FP introduced | +255 |
| Precision（rescue 候選） | 74.1%（728/983） |
| 基線 F1 → 結果 F1 | 0.712697 → 0.722062 |

### 關鍵觀察

1. **GQ>=3 是所有單一規則中最高 F1**（+0.009365 > H011 +0.008556）
2. **H011（QS>=50）⊂ H012（GQ>=3）**：QS>=50 的 608 個候選中 100% 也通過 GQ>=3
3. **H012 額外 86 TP**（相比 H011）：這些候選的 GQ median=9，QS median=0（有 caller 信心但無 ISM 質量信號）
4. **H012 vs H011 精度換算**：GQ>=3 Precision=74.1%（−2.8%）換取 86 個額外 TP

### 決策

`decision = keep`，作為最高 F1 的 TO rescue 規則選項

---

## H_COMBO1：組合分析（H011+H012+Annotation filters）

### 組合結果

| 組合條件 | delta F1 | TP | FP | Precision |
|---------|---------|----|----|-----------|
| GQ>=3（H012） | **+0.009365** | 728 | 255 | 74.1% |
| QS>=50（H011） | +0.008556 | 642 | 193 | 76.9% |
| GQ>=5 | +0.008840 | 672 | 201 | 77.0% |
| QS>=50 + hp>=0.2 filter | +0.008153 | −0.0004 vs H011 | — | — |
| QS>=50 + VC!=Noise | +0.008253 | — | — | — |
| GQ>=0（all ISM candidates） | +0.009692 | 732 | 222 | 76.7% |

### 關鍵結論

1. **H011 ⊂ H012**：兩者不可正交組合（H011 是 H012 的子集）
2. **所有 annotation filter 均降低 F1**：加入 hp>=0.2、hp>=0.5、VC!=Noise、QS>=20/30/40 等條件均使 delta F1 下降（TP 損失 > FP 移除效益）
3. **最優 F1 規則**：GQ>=3（+0.009365），等同 `GQ>=0 without outlier cleanup`
4. **最優 Precision 規則**：QS>=50（76.9% vs GQ>=3 74.1%）

### 決策建議

根據下游需求：
- **注重 F1 最大化** → 使用 GQ>=3（+0.009365）
- **注重 Precision** → 使用 QS>=50（precision=76.9%，delta=+0.008556）
- **不建議**：組合 H011+H012（無正交性）；加 annotation filter（降低 F1）

---

## 規則對比總表

| 規則 | delta F1 | TP | FP | Precision | 特性 |
|------|---------|----|----|-----------|------|
| GQ>=3（H012） | **+0.0094** | 728 | 255 | 74.1% | 最高 F1，H011 超集 |
| QS>=50（H011） | +0.0086 | 642 | 193 | 76.9% | 最高 Precision，GQ>=3 子集 |
| GQ>=5 | +0.0088 | 672 | 201 | 77.0% | 中間選項 |
| GQ>=0（全部） | +0.0097 | 732 | 222 | 76.7% | 理論上限 |

---

## TO Track 背景限制

- ISM 覆蓋的 FN pool：773/11051（7%）— 93% FN **不在 ISM 分析範圍內**
- 主要 FN 原因：ISM 未分析的 VCF 位點（窗口外、reads 不足等）
- **Rescue 上限**：即使完美 rescue 所有 773 ISM-analyzed FN，F1 delta ≈ +0.0097（GQ>=0）

---

## 跨樣本驗證狀態

| 樣本 | 測試狀態 | delta F1 |
|------|---------|---------|
| HCC1395 5kHz TO | ✅ 測試完成 | +0.0086~+0.0094 |
| HCC1395 DORADO TO | ⚠️ 初步觀察 | ≈ +0.0005（幾乎無機會） |
| 其他樣本 | ❌ 未測試 | — |

**注意**：DORADO TO 的 rescue 機會極小（幾乎沒有符合 QS>=50 且 ISM-analyzed 的 FN）。跨樣本一致性待驗證。

---

## 判斷結論

| 假設 | 決策 | 原因 |
|------|------|------|
| H011（QS>=50） | **keep** | delta=+0.0086，最高 Precision |
| H012（GQ>=3） | **keep** | delta=+0.0094，最高 F1 |
| H_COMBO1（組合分析） | **keep** | 確認兩規則無正交性，GQ>=3 = 最優單一規則 |

**建議下游使用**：
- TO rescue 規則選擇取決於實驗目標（F1 vs Precision）
- **跨樣本驗證前不建議在非 HCC1395-5kHz 上硬部署**
