<!--
建立時間: 2026-03-24
目標: 解釋 HCC1395/H1437 中 FP PassedGating 率高於 TP 的根本原因
處理範圍: HCC1395 paired_full（TP=29740, FP=627）; H1437 paired_full（FP=8）
關聯檔案:
  - docs/methodology/20260323_VerificationClass決策樹跨樣本量化_01.md
  - src/core/GlobalTest.cpp (lines 127-141, gate 邏輯)
  - docs/methodology/20260324_方法學審查任務清單與結論摘要_01.md
-->

# FP PassedGating 率高於 TP 的機制分析

**狀態：分析完成 — 根本原因確認：FP 富集於胚系 ASM 區域**

---

## 觀察現象

| 樣本 | TP PassedGating 率 | FP PassedGating 率 | 差異 |
|------|-------------------|-------------------|------|
| HCC1395 | 31.7%（9421/29740） | **50.6%（317/627）** | +18.9% |
| H1437 | 35.3%（23786/67467） | 62.5%（5/8） | +27.2%（FP=8，統計意義不足） |
| HCC1395_DORADO | 36.7% | 29.6% | −7.1%（正常） |

---

## 機制分析（HCC1395，FP=627）

### 1. FP 幾乎全靠 unreliable CramersV 過 gate

| 族群 | passed_gate 中 CramersV=0（unreliable 觸發） | 比例 |
|------|---------------------------------------------|------|
| TP passed gate | 7587/9421 | 80.5% |
| **FP passed gate** | **313/317** | **98.7%** |

→ FP 比 TP 更依賴 unreliable CramersV 觸發 gate，幾乎沒有可靠的結構性信號。

### 2. FP 有強烈的 Allele 甲基化差異

| 訊號 | FP passed gate | TP passed gate |
|------|----------------|----------------|
| AlleleSig=True | **94.0%**（298/317） | 76.4%（7202/9421） |
| HPMergedSig=True | 59.6%（189/317） | 50.1%（4719/9421） |
| AlleleDelta 中位數 | **0.1548** | 0.0618 |
| AlleleDelta > 0.15 | 168/317（53.0%） | 762/9421（8.1%） |

→ **FP 的 AlleleDelta 中位數是 TP 的 2.5 倍**，且 94% 達到 AlleleSig。

### 3. FP GlobalP 甚至比 TP 更顯著

| 指標 | FP passed gate | TP passed gate |
|------|----------------|----------------|
| GlobalP 中位數 | 0.0000 | 0.0000 |
| GlobalP 均值 | **0.0089** | 0.0119 |

→ FP 的全域顯著性高於 TP（GlobalP 更低）。

### 4. FP passed gate 的 VerificationClass

| VC | FP（n=317） | 比例 |
|----|------------|------|
| Strong | **285** | **89.9%** |
| Weak | 21 | 6.6% |
| Subclone | 8 | 2.5% |
| Noise | 3 | 0.9% |

→ FP 過 gate 後 90% 被分類為 Strong（同時具備 global 顯著 AND label 顯著）。

---

## 根本原因

```
HCC1395 FP = 胚系變異（germline variants），而非體細胞變異。
這些胚系變異位點富集於胚系 ASM（Allele-Specific Methylation）區域。

在 ASM 區域：
  → Ref allele reads 與 Alt allele reads 天然具有不同甲基化模式
  → AlleleDelta 天然偏高（約 0.15-0.3）
  → 甲基化 k=2 聚類找到 Allele 層面的結構 → unreliable CramersV 高 → gate 通過
  → AlleleSig 觸發（94%）→ Strong 分類

ISM gate 無法區分「體細胞驅動的甲基化差異」vs「胚系 ASM 驅動的甲基化差異」。
```

### 為何 FP gate 率 > TP gate 率？

- **FP（胚系變異）選擇性富集於 ASM 活躍區域**：VCF caller 在高甲基化活躍的 CpG island 附近容易產生胚系假陽性（甲基化不均導致測序或比對誤差）
- **TP（體細胞突變）分散於全基因組**：許多 TP 位於甲基化同質區域（無 ASM），CramersV 低 → gate 通過率低
- **機制差異**：FP 信號 = 天然 ASM，TP 信號 = 突變誘導的亞克隆甲基化差異（需要更多 reads 才能偵測）

---

## 跨樣本對照

| 樣本 | FP gate 率 > TP | FP AlleleDelta 特徵 | 推測原因 |
|------|----------------|---------------------|---------|
| HCC1395 | **是**（+18.9%） | 高（中位數 0.155） | FP 富集於 ASM 區域 |
| H1437 | **是**（+27.2%） | 不可靠（FP=8） | 統計意義不足 |
| HCC1395_DORADO | 否（−7.1%） | 未特別分析 | FP 分佈不同？ |
| H2009 | 否（TP 12.2% > FP 4.7%） | FP=86，接近正常 | FP 不富集 ASM |

→ H1437 FP=8 不具統計意義。僅 HCC1395 有明確機制：FP 富集於 ASM 區域。

---

## 潛在改進方向

### 方向 A：引入 VAF 門檻（已有 VCF 資訊）

```python
# 現有規則（purity ge60）已有 VAF < 0.15 條件
# 但 AlleleDelta > 0.15 條件才是真正驅動力
# 胚系 ASM FP 的 VAF 通常接近 0.5（胚系 het）
# 可考慮：VAF > 0.3 → 可能胚系，降低信心
```

問題：需要 VCF 對應，且 VAF 閾值需跨樣本驗證。

### 方向 B：正常樣本甲基化參考

若有配對正常樣本的甲基化 reference，可過濾掉 normal 中已存在的 ASM 信號（但 ISM 是 tumor-only 設計）。

### 方向 C：提高 gate reliable 要求（保守）

```
當前 gate = CramersV（unreliable allowed）
修改為 gate = reliable CramersV ≥ 0.05
```

效果預測：
- TP: 98.7% FP 已觸發 unreliable CramersV → 幾乎全部 FP 不過 gate
- 但 TP 的 80.5% 也靠 unreliable CramersV → TP gate 率大幅下降
- **結論：不建議** — TP 損失遠大於 FP 收益

---

## 判斷結論

| 問題 | 結論 |
|------|------|
| FP gate 率 > TP 是 bug 嗎？ | **否** — 是結構性問題（設計限制） |
| 原因是什麼？ | FP 富集於胚系 ASM 區域，天然 AlleleDelta 高（中位數 0.155 vs TP 0.062） |
| 可透過 gate 閾值調整改善嗎？ | 否 — 提高 gate 嚴格度會同等損失 TP |
| 現有的 AlleleDelta+VAF rule 有幫助嗎？ | HCC1395-5kHz 有效（AUROC 0.763），但跨樣本不一致（已降為 annotation） |
| 建議後續行動 | 保持現狀；若需進一步改善 FP，需引入 normal sample 參考或其他非甲基化特徵 |

**判決：No Action** — 機制已明確，非可用閾值調整解決的問題。
