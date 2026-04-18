<!--
建立時間: 2026-04-02 02:00
目標: LOH 研究週報完整敘事稿（v2 大修版），敘事重心：偵探故事線 — 發現→懷疑→驗證→追查→修正
處理範圍: 2026-03-25 ~ 2026-03-31 全部研究內容
關聯檔案:
  - 02_ppt_slide_outline.md
  - 03_slide_layout_and_script.md
  - 05_loh_read_depth_and_copy_number_analysis.md
-->

# TO 模式 LOH 系統性偏差調查：從發現到機制剖析

**報告期間**：2026-03-25 ~ 2026-03-31
**報告人**：[報告人]
**資料規模**：748,391 regions × 116 features × 7 個癌症細胞系 × 2 種分析模式（Paired / TO）
**敘事主線**：偵探故事 — 發現異常 → 懷疑偏差 → 驗證確認 → 追查根因 → 提出修正

**一句話摘要**：發現 TO 模式的 LOH 判定存在系統性偏差，完整驗證機制並確認是 LongPhase-TO self-phasing 造成，提出修正方向。

---

## Section 1 — 為什麼要研究 LOH：從 Paired 延伸到 Tumor-Only

*（對應 Slide 2）*

### 動機：Paired FP 太少，需要 TO 的高 FP Rate 來深度驗證 ISM

InterSubMod (ISM) 是我們開發的 read-level 表觀遺傳分析工具，透過整合甲基化模式、somatic SNV 和 haplotype 資訊，嘗試區分 somatic variant calling 的 TP（真陽性）和 FP（假陽性）。

上一階段在 **Paired 模式**下建立基線後，遇到一個根本問題：**FP 數量太少，無法可靠評估 ISM 的特徵區分能力**。

| 模式 | FP 數量 | FP Rate | TP:FP 比例 | 統計力 |
|------|---------|---------|-----------|--------|
| **Paired** | 3,429 / 328,699 | **~1%** | 95:1 | 不足 |
| **TO** | 128,382 / 419,692 | **~31%** | 2.3:1 | 充足 |

**FP 來源根本不同**：
- **Paired FP** = calling 殘留錯誤（~1%）——caller 在有 normal 對照的情況下仍然誤判的少數位點
- **TO FP** = germline variants 被誤判為 somatic（~31%）——因為**沒有 normal 樣本對照**，caller 無法有效過濾 germline variants

**本週核心問題**：ISM 在 TO 模式下能否區分 TP/FP？LOH 在其中扮演什麼角色？

![TP/FP 比例對比](figures/O01_fig04_truth_label_composition.png)

---

## Section 2 — 技術背景：ISM 分析框架與 HP Tag 差異

*（對應 Slide 3）*

### ISM 如何分析甲基化模式

ISM 分析每個 somatic variant 位點周邊 **±5000bp** 區域內的 CpG 甲基化模式。對每個 region：

1. **讀取分群**：依據 **LongPhase haplotagged BAM** 中的 HP tag 將 reads 分為 HP1（父源）和 HP2（母源）。HP tag 由 LongPhase 預先計算，**不是 ISM 自行判定**。
2. **甲基化矩陣**：每條 read × 每個 CpG 位點的甲基化狀態（0/1/NA）
3. **距離計算**：reads 間的甲基化模式相似度（L1, L2, Bernoulli 等）
4. **統計檢驗**：PERMANOVA 檢驗不同 haplotype 之間的甲基化模式是否顯著不同
5. **品質評分 (QS)**：先前嘗試使用 **Quality Score 加分表**整合多組件評分（Penalties + Bonuses），Paired 模式下 AUC=0.754

### Paired 與 TO 的 Haplotype 差異

這是這個偵探故事的起點——Paired 和 TO 使用**不同的 phasing 策略和 HP tag 格式**：

| 特性 | Paired 模式 | TO 模式 |
|------|-----------|---------|
| **Phasing 工具** | LongPhase-s | LongPhase-TO |
| **Phasing anchor** | **Germline SNV** | **Tumor SNV（含 somatic）** |
| **HP tag 格式** | `HP:Z:1-1`（字串，haplotype-block） | `HP:i:11`（整數編碼） |
| **潛在問題** | 無 | Somatic variant 參與自身的 phasing → circular dependency |

**LOH-like 操作定義**：
- `HP_Ratio = min(HP1, HP2) / (HP1 + HP2)`
- `HP_Ratio < 0.1 或 > 0.9` 且 `effective_hp_reads ≥ 30` → 判定為 **LOH-like**
- Enrichment = FP LOH% / TP LOH%（>1 = FP 標記，<1 = TP 標記）

---

## Section 3 — 數據品質檢查：HP Tag Bug 發現與修正

*（對應 Slide 4）*

### AI 執行無報錯，但看結果圖片發現數據異常

這是一段典型的「表面正常、內裡有毒」的除錯經歷：

1. **程式執行無報錯**：讓 AI 運行 ISM pipeline 處理 TO 數據，程式完整執行完畢，沒有任何 error 或 warning
2. **檢視結果圖片發現異常**：但在檢視結果圖片時，發現 TO 的 HP 數據異常——大量 reads 缺少 haplotype 歸屬，TO 的 Tier A/A+ 佔比只有 ~50%（Paired 是 ~90%）
3. **追查根因**：追查發現 `ReadParser.cpp` 只處理字串格式（`HP:Z:1-1`），TO BAM 的整數格式（`HP:i:11` / `HP:i:21` / `HP:i:33`）被**靜默忽略**
4. **修正程式碼**：修改 ISM 程式碼支援兩種輸入格式：

| BAM Tag | 修正前 | 修正後 | 語義 |
|---------|--------|--------|------|
| `HP:i:11` | 忽略 | `"1-1"` | Haplotype 1, somatic-supported |
| `HP:i:21` | 忽略 | `"2-1"` | Haplotype 2, somatic-supported |
| `HP:i:33` | 忽略 | `"3"` | Ambiguous haplotype |

5. **驗證修正效果**：修正後 TO 的 Tier A/A+ 從 ~50% 升至 ~88%（接近 Paired 的 ~90%）

**處置**：修正前所有 TO 數據作廢，7 個樣本全量重跑 ISM pipeline → 重建 master dataset

![HP Fix 後 TO Tier 分佈](figures/O07_fig06_to_tier_distribution_post_fix.png)

---

## Section 4 — 發現異常：LOH Enrichment 方向翻轉

*（對應 Slide 5）*

### 同一個指標，Paired 和 TO 的含義完全相反

修正 HP tag bug 並全量重跑後，第一個驚人發現浮現——**同一個 LOH enrichment 指標，在 Paired 和 TO 之間方向完全翻轉**：

| 模式 | LOH Enrichment | 含義 | 機制 |
|------|:--------------:|------|------|
| **Paired** | **1.194x** | FP 更多 LOH-like | Paired FP 的 HP support 低 → 隨機波動產生假 LOH |
| **TO** | **0.805x** | **TP 更多 LOH-like** | TO phasing 中 somatic allele 造成系統性 HP 偏移 |

**7/7 樣本全部一致**：
- Paired 全為 FP-enriched（1.02-3.18x）
- TO 全為 TP-enriched（0.852-0.956x）
- COLO829 因低深度（~30x）enrichment 接近中性，但方向仍為 TP-enriched

**這直接導致 TO QS 失效**：QS 的 LOH penalty 假設「LOH = FP」，在 TO 下此假設不成立，反而懲罰 TP。

![LOH Enrichment 熱圖 — 7 樣本 × 2 模式](figures/O03_fig04_loh_enrichment_heatmap.png)

> 此熱圖清楚展示各樣本的 LOH 比例差異與 Paired/TO 方向翻轉的對應關係

---

## Section 5 — 量化異常：288K 同位點配對分析

*（對應 Slide 6）*

### 同一基因組位點，TO 單獨判 LOH 是 Paired 的 16-52 倍

為了量化翻轉的規模，我們對 **288,609 個同位點 TP 配對**進行分析——同一基因組位點，分別觀察 Paired 和 TO 的 LOH 判定是否一致：

| 類別 | 數量 | 比例 |
|------|-----:|:----:|
| 共識：非 LOH | 157,774 | 55.0% |
| 共識：LOH | 81,738 | 30.6% |
| **TO-only LOH** | **39,724** | **13.8%** |
| Paired-only LOH | 1,856 | 0.6% |

關鍵數字：
- **一致率 85.5%**（共識 LOH + 共識非 LOH）
- 不一致 41,580 位點中 **95.5% 是 TO 單獨判 LOH**
- Per-sample TO excess ratio = **16-52 倍**
- 共識 LOH 30.6% 代表在雙方都認定 LOH 的情況下，LOH 確實在所有數據中佔約三成

![LOH Concordance 綜合分析 — 5 面板](figures/concordance_summary_5panel.png)

> (A) 全域 2×2 concordance matrix；(B) 每個樣本四分類比例；(C) 每個樣本 TO excess ratio（16-52x）；(D) Extreme LOH rate 方向翻轉；(E) 7 樣本 Enrichment 全部一致翻轉

---

## Section 6 — 懷疑 TO 系統性偏差：所有 TO 特徵都無效

*（對應 Slide 7）*

### 9 項系統性觀察 → TO 模式全面 Negative（AUC < 0.58）

LOH 翻轉不是孤立事件。我們對 **748,391 regions × 116 features** 進行了九項系統性觀察（O1-O8 + O10），加上四項假說驗證（O11-O13 + N4），共產出 **82 張圖表**。結論是：**TO 模式下所有特徵鑑別力不足**，不只是 LOH。

| 特徵 | Paired AUC | TO AUC | 備註 |
|------|:----------:|:------:|------|
| GQ | **0.811** | — | TO 無 GQ |
| AF | 0.665 | 0.418 | 反轉 |
| ISM 甲基化 | 0.543 | 0.530 | 微弱 |
| LOH | 0.579 | 0.537 | 微弱 |
| **QS 綜合** | **0.754** | **0.497** | 隨機 |

- **TO 全面 Negative**：全部 AUC < 0.58，5/9 甲基化特徵 Paired→TO 方向反轉
- **QS 完全失效**：Paired AUC=0.754 → TO AUC=0.497（等於隨機猜測）
- **LOH penalty 觸發率反轉**：TO TP **44.5%** vs FP **35.8%** → 反向懲罰 TP

![TO QS Waterfall](figures/O02_fig02_qs_component_waterfall_to.png)

---

## Section 7 — 排除假說：不是 Read Depth 不足造成的

*（對應 Slide 8）*

### TO-only LOH 的讀取深度反而更高

看到這些異常，一個自然的假設是：「也許 TO 的 read depth 較低，導致 HP ratio 隨機偏移到極端值？」數據直接否定了這個假說。

**四象限 Read Depth 比較**（288,609 matched loci）：

| 象限 | n | TO eff_hp median | Paired eff_hp median |
|------|-------|-----------------|---------------------|
| both_LOH | 88,109 | **60** | **60** |
| TO_only_LOH | 39,978 | **68** | **61** |
| Paired_only_LOH | 1,874 | 65 | 43 |
| neither_LOH | 158,648 | **79** | **79** |

**關鍵發現**：
- 如果是 depth 不足造成：TO-only LOH 的 reads 應顯著**低於** both_LOH → **不成立**
- TO-only median effective_hp_reads = **68**，反而**高於** both_LOH 的 60，Cohen's d = +0.29（方向相反）
- **決定性數據**——同一基因組位點：
  - Paired HP_Ratio 中位數 = **0.509**（完全平衡）
  - TO HP_Ratio 中位數 = **0.026**（極端 LOH）
  - Cohen's d = **-1.20**（巨大效應量）
- **86.5%** 的 TO-only LOH 位點在 Paired 下完全平衡（HP_Ratio 0.2-0.8）

![四象限 Effective HP Reads 分佈](figures/concordance_quadrant_depth.png)

> Box plot 清楚顯示 TO-only LOH（橘色）的讀取深度與 both_LOH（綠色）相當甚至更高

---

## Section 8 — 確認偏差：HP Ratio 散點圖的直接可視化

*（對應 Slide 9）*

### 同一基因組位點，Paired 看到平衡，TO 被推到極端

HP Ratio 散點圖是 self-phasing circular dependency 的**最直觀視覺證據**：

- **X 軸** = Paired HP_Ratio，**Y 軸** = TO HP_Ratio
- **綠色（both_LOH）**：聚集在四角（雙方都判 LOH）
- **藍色（neither）**：聚集在中央（雙方都不判 LOH）
- **橘色（TO-only LOH）**：分佈在**上下兩條水平帶** — X 軸（Paired）在 0.3-0.7（平衡），Y 軸（TO）在 <0.1 或 >0.9（極端）

這是 LongPhase-TO self-phasing 的最直觀證據：在同一基因組位點，Paired phasing 看到對稱的 haplotype 分配，但 TO phasing 因為將 somatic variant 自身納入 phasing anchor，將 reads **系統性推向一側**。

![HP Ratio Paired vs TO — 四象限上色](figures/concordance_hp_ratio_scatter.png)

![Read Depth vs HP Ratio（2D hexbin）](figures/concordance_depth_vs_hpratio.png)

> Hexbin 圖顯示：即使在高 depth（>100）區域，TO HP ratio 仍被推到極端值 — 進一步排除「低 depth 隨機偏移」假說

---

## Section 9 — 補充驗證：LOH Read Depth 與 Copy Number 分析

*（對應 Slide 10）*

### LOH 大部分不是 Allelic Deletion

為了完整排除 copy number 的干擾，我們進一步分析了 LOH 區域的 read depth 和 copy number 組成。

#### Read Depth：LOH ≈ 0.73x（不是理論的 0.5x）

如果 LOH 代表一個 allele 真正被刪除（copy-loss LOH），region 的總 read 數應該約為 non-LOH 的 50%。實際觀察到的 ratio 遠高於此：

**Per-Sample TP Read Depth Ratio（LOH / non-LOH）**：

| Sample | Paired TP | TO TP | 備註 |
|--------|-----------|-------|------|
| HCC1395 | 0.68 | 0.72 | 最低，接近 copy-loss |
| H1437 | 0.69 | 0.74 | |
| HCC1395_DORADO | 0.71 | 0.74 | |
| HCC1937 | 0.72 | 0.75 | |
| COLO829 | 0.74 | 0.77 | |
| H2009 | 0.77 | 0.82 | |
| HCC1954 | **0.97** | **1.02** | LOH depth ≈ non-LOH，大量 copy gain |

大多數 sample 落在 0.68-0.82 之間，遠高於純 copy-loss 的理論值 0.5。

![Per-Sample LOH Read Depth Ratio](figures/loh_read_depth_ratio_per_sample.png)

> 紅線 0.5 = 純 copy-loss 理論值，紅線 1.0 = copy-neutral。大多數 sample 遠高於 0.5

#### Copy Number 組成：60% Copy-Neutral

| Mode | Copy-neutral | Copy-loss + Deep del | Copy-gain + High amp |
|------|-------------|---------------------|---------------------|
| **Paired** | **61.8%** | 29.8% | 8.4% |
| **TO** | **60.1%** | 24.2% | **15.8%** |

- **60% copy-neutral**（可能是 UPD 或 phasing artifact），24-30% copy-loss，8-16% copy-gain
- HCC1954 最極端：ratio ≈ 1.0（LOH depth ≈ non-LOH），已知大量基因組重排
- **TO copy-gain LOH 比例 = Paired 的 2 倍**（15.8% vs 8.4%）→ TO 在 copy-gain 區域有額外的系統性偏差

![CN LOH Composition Per Sample](figures/cn_loh_composition_per_sample.png)

---

## Section 10 — 根因確認：Self-Phasing Circular Dependency

*（對應 Slide 11）*

### LongPhase-TO 的 Somatic Variant 參與自身 Phasing — 造成循環依賴

所有線索指向同一個根因：**TO 使用 somatic variant 做 phasing anchor，TP 的 ALT reads 驅動 HP 偏斜，形成循環依賴**。

#### 四步邏輯鏈

1. **TO joint phasing**：LongPhase-TO 將 germline + somatic variant 共同建構 phasing graph（`PhasingGraph.cpp` addEdge() 確認），但 somatic variant 只有 ALT allele 貢獻 haplotype voting
2. **Self-phasing circular dependency**：TP 的 somatic ALT reads 天然偏向一個 haplotype → 形成 haplotype block → HP tagging 繼承偏斜 → HP_Ratio 偏離 0.5
3. **低 read 區域被串入長 block**：TO 在 read 數量低的區域也嘗試串接，造成 phasing block 特別長，但 haplotag 結果品質很差，不如高可信度的 Paired phasing
4. **TP 更容易被判 LOH-like**：TO TP LOH 觸發率 44.5% vs FP 35.8% → LOH 成為 TP 弱標記 → Enrichment 0.805（7/7 一致）

#### 量化驗證

- 39,724 位點從 Paired non-LOH → TO LOH-like
- 其中 **71.6%** 在 TO 的 min(HP) = 0（完全單 haplotype），同位點 Paired HP1:HP2 = 8:8
- **Shapley 分解**：TP rate 變動貢獻 105.6%，FP rate 僅 -5.4%

![LOH Penalty 觸發率對比](figures/O02_fig03_qs_loh_penalty_trigger_rate.png)

![HP Ratio 分佈（per sample × mode）](figures/O03_fig01_hp_ratio_distribution_by_sample_mode.png)

---

## Section 11 — 否決的假說：甲基化方向全部關閉

*（對應 Slide 12）*

### O11/O12/O13 — 三項甲基化假說全部否決，問題在 Phasing 不在 Methylation

在確認 LOH 翻轉後，我們嘗試用甲基化特徵解決 TO 問題，但**三項假說全部因 confound 被否決**。真正的問題在 phasing 層面，不是 methylation。

| 假說 | 結果 | 否決原因 |
|------|------|---------|
| **O11**：Within-group heterogeneity | NEGATIVE | epipolymorphism AUC 0.845 → n_reads 校正後 **0.530** → read count confound |
| **O12**：LOH 甲基化場景區分 | NEGATIVE | AlleleDelta = AF confound + L2 collider bias → 全 corrected AUC < 0.58 |
| **O13**：跨區域甲基化 correlation | NEGATIVE | shared read count confound → 控制後差異消失 |
| **N4**：HP0 filter | NEGATIVE | TP loss ≤2% 下 FP removal = 0% |

**教訓**：L1→L2→L3 三層驗證架構，任何 L2 有效的特徵（AUC > 0.6）必須通過 L3 AF-bin 分層交叉驗證，否則可能是 confound 造成的假象。

---

## Section 12 — 已完成修正：QS Mode-Aware

*（對應 Slide 13）*

### 第一步修正 — 消除 TO LOH Penalty 的反效果

基於完整的機制證據，我們在 TO 模式下停用 LOH penalty 和 verify bonus。

| QS 組件 | Paired | TO（修改後） | 理由 |
|---------|:------:|:----------:|------|
| LOH penalty | -25 | **0** | 方向反轉，懲罰 TP |
| Verify bonus | +15 | **0** | 依賴 LOH 判定 |
| 其他 | 不變 | 不變 | 方向一致 |

- **先前做法**：使用 Quality Score 加分表整合多組件，在 Paired 有效（AUC=0.754），TO 失效（AUC=0.497）
- **修正效果**：TO QS AUC 0.497 → ~0.546（理論推算），from harmful to neutral
- **已提交**：commit `b9eaba7`
- **注意**：0.546 仍接近隨機，此修正意義是「消除反效果」，不是讓 QS 變有效

---

## Section 13 — 修正方向：TO Phasing 程式碼改善可能性

*（對應 Slide 14）*

### 根因在 LongPhase-TO，需分析程式碼層面的改善空間

既然問題的根因已經定位在 LongPhase-TO 的 phasing 策略，下一步需要分析是否有程式碼層面的修正空間。

#### 分析方向

1. **LongPhase-TO phasing anchor 策略**：目前將 somatic variant 納入 phasing graph，可否排除或降權？
2. **Phasing block 品質過濾**：TO 在低 read 區域串出的長 block 品質差，可否加 block quality threshold？
3. **HP tagging 後處理**：在 ISM 層面，可否根據 phasing block 大小 / read support 動態調整 LOH 判定門檻？
4. **Paired phasing 作為 benchmark**：利用 288K 同位點數據，量化 TO phasing 在哪些條件下可以信任

#### 已有的改善基礎

- 288,609 matched loci 數據可以做 TO phasing quality calibration
- 已知 86.5% 的 TO-only LOH 在 Paired 下平衡 → 這些是可以被校正的 false LOH
- Copy-gain 區域額外偏差（2x Paired）→ 需要 depth-aware 校正

#### 待執行

- 分析 LongPhase-TO 的 phasing block 長度分佈 vs Paired
- 比較 block 內 read count 和 haplotag 一致性
- 評估排除 somatic variant from phasing anchor 的可行性

---

## Section 14 — 核心結論

*（對應 Slide 15）*

### Paired 和 TO 是根本不同的問題空間

| 維度 | Paired | TO |
|------|:------:|:--:|
| FP rate | 1% | 31% |
| LOH enrichment | FP-enriched (1.194x) | TP-enriched (0.805x) |
| 特徵方向 | — | 5/9 反轉 |
| QS AUC | 0.754 | 0.497 (隨機) |
| TO-only LOH excess | — | 16-52x over Paired-only |
| 自我驗證 | — | Cohen's d = -1.20 |

### 結論層次

1. **Paired/TO 分離**不是選項，是**必要條件**
2. TO LOH 偏差源自 **self-phasing circular dependency**，非 read depth 或 copy number
3. 60% LOH 是 copy-neutral → 部分可能是 phasing artifact
4. 甲基化特徵無法解決 TO 問題（O11-O13 全否決）→ 問題在 phasing 層
5. 下一步：TO phasing 程式碼分析與修正方案設計

### Limitations

- 所有結論基於 7 個高 purity 細胞系（purity ≈ 100%），臨床低 purity 需額外驗證
- LOH-like 閾值（0.1/0.9）為操作定義，需 sensitivity analysis

---

## Section 15 — 討論與下週方向

*（對應 Slide 16）*

### 下週行動

| 優先級 | 行動 | 依據 |
|--------|------|------|
| **P0** | TO phasing block 品質分析 | 根因在 phasing |
| **P0** | Paired/TO 分離策略框架 | 方向完全相反 |
| **P1** | Paired ML 特徵集（GQ + AF + LOH subtype） | GQ AUC=0.811 |
| **P2** | TO phasing anchor 修正方案評估 | 排除 somatic variant |

### 開放討論問題

1. LongPhase-TO 排除 somatic variant 後，phasing 品質是否足夠？
2. TO 獨立 QS 框架應以什麼為 anchor？
3. Paired phasing 可以做 TO 的 calibration reference 嗎？
