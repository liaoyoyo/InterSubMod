<!--
建立時間: 2026-04-06 23:30
目標: PPTX 完整敘事稿（偵探故事線）
處理範圍: 2026-03-31 ~ 2026-04-06 LOH 雙定義與特徵探索全面關閉
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/04/20260406_研究週報_20260331_20260406_LOH雙定義與特徵探索全面關閉_01.md
-->

# 完整敘事稿：LOH 雙定義與特徵探索全面關閉

## 敘事主軸：偵探故事 — 從 LOH 定義存疑到全面關閉特徵探索

---

### 開場：問題設定（~2 分鐘）

上週我們修正了 HP integer tag bug、完成 LOH Evidence Panel 四輪研究、鎖定 Phase 1A methyl+context 方向。但一個關鍵問題懸而未決：**LOH 定義本身可信嗎？ISM 在 LOH 區域的分析還有意義嗎？**

本週，我們從 LOH 的根本定義出發，系統性驗證兩種 LOH 定義的可信度與一致性，量化 LOH 對 ISM 的影響，並窮盡所有剩餘特徵方向——最終得到一個決定性結論：**TO 特徵探索全面關閉，研究方向正式轉向 read-level epigenetic characterization。**

---

### 第一幕：LOH 有兩種定義，它們一致嗎？（~3 分鐘）

ISM 分析中存在兩種 LOH 定義：
- **LOH.bed**（區域層級）：LongPhase-TO 輸出，基於 phased genotype ratio
- **HP_Ratio LOH**（位點層級）：ISM 計算的 HP1/(HP1+HP2) > 0.9

分析 7 樣本 419,692 筆 TO mode 數據，建立**四象限分佈**：
- Q1 Both LOH：26.7%（兩者一致）
- Q2 HP-only：15.2%（ISM 判 LOH，LOH.bed 不判）
- Q3 LOH.bed-only：0.07%（僅 286 筆，近乎為零）
- Q4 Neither：58.1%

**核心發現**：HP Imbalance 是 LOH.bed 的超集（Sensitivity 99.7-100%）。Q2 的 AF median=0.47，與真 LOH 的 AF 分佈完全不同（d=-1.04）——Q2 是 phasing 偏差，不是弱 LOH。

---

### 第二幕：SEQC2 — 找到外部金標準（~3 分鐘）

但 LOH.bed 本身可信嗎？我們找到了 SEQC2（FDA 多平台共識）CNV benchmark：
- 320 個 autosomal LOH entries，覆蓋 1,490.5 Mb
- 基於 WGS short-read + SNP array + 多中心共識——與 LongPhase-TO（ONT long-read phasing）**方法學零重疊**

**驗證結果**：Jaccard=0.928（HCC1395 5kHz: 0.9277, Dorado: 0.9290）。Sensitivity=0.961, Precision=0.964, F1=0.963。

**這是本週最重要的發現之一**：兩種完全獨立的技術對同一細胞株的 LOH 判定幾乎完全吻合。LOH.bed 可以作為後續所有分層分析的可靠基礎。

---

### 第三幕：LOH 對 ISM 的影響 — PERMANOVA 全面失效（~2 分鐘）

LOH.bed 可信後，接下來的問題是：LOH 區域內，ISM 還能正常工作嗎？

答案是**不能**：
- LOH 區域 PERMANOVA valid rate 僅 **5-6%**（需兩群各 ≥3 reads）
- 94% 位點 min(HP) < 3（只有一條 haplotype 有 reads）
- LOH 內甲基化特徵 **7/7 樣本 AUC~0.50**（與隨機無異）

根因很清楚：LOH 區域只有一條 haplotype 上有 reads，ISM 的 HP1 vs HP2 比較失去統計意義。

**已實作修正**：QS mode-aware（commit b9eaba7）——TO 模式下停用 LOH penalty 和 verify bonus。

---

### 第四幕：LOH 不可作 Filter — 10/10 全面 FAIL（~2 分鐘）

既然 LOH 內甲基化無效，能否直接「移除 LOH 區域的 variants」來提升 FP 過濾？

**不可以**：
- LOH 區域 FP rate=0.239 vs Non-LOH 0.338 — LOH 內反而 TP 比例更高
- 10 種篩選策略在 7 樣本上全面 FAIL（TP loss > 2% 安全線）

生物學解釋：LOH 是 tumor suppressor 失活的常見機制，LOH 內大多數 variants 是真正的 somatic TP。移除 LOH 會大量丟失 TP。

---

### 第五幕：Non-LOH 也無突破（~2 分鐘）

LOH 內無解，那 Non-LOH（58.1%）呢？

最高 AUC=0.643（HPFineNGroups），但驗證後是 read count proxy。排除 confound 後**無特徵 > 0.58 門檻**。

多特徵組合：AlleleDelta + CramersV + PairwiseMedianDist 三特徵 Voting AUC=**0.577 < 0.58**。

cnLOH 的 PairwiseMeanDist 看似 0.5865 > 0.58，但這是 **Simpson's Paradox**：per-sample mean=0.4987，被 H2009（48.6%）主導。

**結論**：問題是全域性的——不是 LOH 特異的，LOH 與 Non-LOH 的 AUC 差距 < 0.06。

---

### 第六幕：平行線索 — 甲基化三維度 NEGATIVE + TO 60+ 特徵 NO-GO（~2 分鐘）

同時進行的其他探索方向：

**O11-O13 甲基化三維度全 NEGATIVE**：
- O11 Heterogeneity：AUC 0.845→0.530 after n_reads correction
- O12 LOH 甲基化場景：corrected AUC < 0.58 + L2 collider bias
- O13 跨區域 correlation：shared read count confound

**G1-G7 TO VCF 特徵 NO-GO**：48 圖 60+ 特徵全 AUC<0.64，安全約束下 FP removal=0%。

**唯一亮點**：ASM（allele-specific methylation）32-66% SNV 位點驗證，ISM PERMANOVA 是唯一量化 haplotype 間 entropy imbalance 的工具。ISM 不是無用——用途在 characterization 而非 filtering。

---

### 收尾：戰略結論與下一步（~2 分鐘）

**本週 14 結論中 5 個新 NEGATIVE + 2 個新 CONFIRMED**，研究圖景大幅收斂：

- SEQC2 驗證 LOH.bed 可信 → 分層分析基礎穩固
- LOH 內甲基化失效 + filter FAIL → LOH 用途在 stratification 非 filtering
- Non-LOH + 多特徵組合都無突破 → TO post-hoc 特徵探索正式全面關閉
- ASM 是 ISM 唯一亮點 → 論文定位為 read-level epigenetic characterization

**下一步**：
1. Phase 2A 規劃：PON-only phasing 全量重跑 + ISM 重跑
2. Normal methylation reference baseline 建立
3. 研究全景文件 v2.0 更新
