<!--
建立時間: 2026-04-02 01:00
目標: Phase 2 — PPT 投影片段落主軸（v2 大修版），敘事重心：LOH 偏差調查的偵探故事線
處理範圍: 2026-03-25 ~ 2026-03-31 LOH 研究週報
關聯檔案:
  - 01_full_narrative_report.md
  - 03_slide_layout_and_script.md
  - 05_loh_read_depth_and_copy_number_analysis.md
-->

# LOH 研究週報 PPT — v2 大修版

**報告期間**：2026-03-25 ~ 2026-03-31
**預估頁數**：16 張（含封面與結尾）
**建議時長**：20-25 分鐘
**敘事主線**：偵探故事 — 發現異常 → 懷疑偏差 → 驗證確認 → 追查根因 → 提出修正

---

## Slide 1 — 封面

**標題**：TO 模式 LOH 系統性偏差調查：從發現到機制剖析

- 報告人 / 日期 / 資料規模（748,391 regions × 116 features × 7 cancer cell lines × 2 modes）
- 一句話摘要：「發現 TO 模式的 LOH 判定存在系統性偏差，完整驗證機制並確認是 LongPhase-TO self-phasing 造成，提出修正方向」

---

## Slide 2 — 為什麼要研究 LOH

**標題**：從 Paired 延伸到 Tumor-Only — 為了統計檢定力

**核心訊息**：Paired FP 太少無法有效評估 ISM，借 TO 的高 FP rate 深度驗證。

**要點**：
- Paired FP = 3,429 / 328,699（~1%）→ TP:FP ≈ 95:1，統計力不足
- TO FP = 128,382 / 419,692（~31%）→ TP:FP ≈ 2.3:1，充足
- **FP 來源差異**：Paired FP 主要是 calling 殘留錯誤（~1%）；TO FP 主要是 **germline variants 被誤判為 somatic**（~31%，因無 normal 對照過濾）
- **本週核心問題**：ISM 在 TO 模式下能否區分 TP/FP？LOH 扮演什麼角色？

![TP/FP 比例對比](figures/O01_fig04_truth_label_composition.png)

---

## Slide 3 — 技術背景：ISM 分析框架與 HP Tag 差異

**標題**：InterSubMod 分析架構與 Paired/TO Haplotype 差異

**核心訊息**：ISM 以 ±5000bp 區域分析甲基化模式，HP tag 來自 LongPhase 的 haplotagged BAM，Paired 與 TO 的 tag 格式不同。

**要點**：
- ISM 分析每個 somatic variant 位點周邊 **±5000bp** 區域內的 CpG 甲基化模式
- Reads 分群依據 **LongPhase haplotagged BAM** 中的 HP tag（非 ISM 自行判定）
  - **Paired**: LongPhase-s 使用 germline SNV 做 phasing anchor → HP tag 格式 `HP:Z:1-1`（字串，haplotype-block）
  - **TO**: LongPhase-TO 使用 tumor SNV（含 somatic）做 joint phasing → HP tag 格式 `HP:i:11`（整數編碼）
- **關鍵差異**：TO 的 phasing anchor 包含 somatic variant 本身 → 潛在 circular dependency
- LOH-like 定義：`HP_Ratio < 0.1 or > 0.9`（effective_hp_reads ≥ 30）
- Enrichment = FP LOH% / TP LOH%（>1 = FP 標記，<1 = TP 標記）
- 先前嘗試：使用 **Quality Score 加分表** 整合多組件評分（Penalties + Bonuses），Paired AUC=0.754

---

## Slide 4 — 數據品質檢查：HP Tag Bug 發現與修正

**標題**：結果看起來不對 — 追查到 HP 輸入格式問題

**核心訊息**：AI 自動執行完成無錯誤，但檢視結果圖片發現數據異常，追查到 ISM 無法解析 TO 的 HP:i 整數格式，修正後重跑所有數據。

**要點**（簡化故事線）：
1. **現象**：讓 AI 運行 ISM pipeline，程式執行無報錯，但**檢視結果圖片時**發現 TO 的 HP 數據異常 — 大量 reads 缺少 haplotype 歸屬
2. **診斷**：追查發現 `ReadParser.cpp` 只處理字串格式（`HP:Z:1-1`），TO BAM 的整數格式（`HP:i:11`/`HP:i:21`/`HP:i:33`）被靜默忽略
3. **修正**：修改 ISM 程式碼支援兩種輸入格式（`HP:i:11→"1-1"`, `HP:i:21→"2-1"`, `HP:i:33→"3"`）
4. **驗證**：修正後 TO 的 Tier A/A+ 從 ~50% 升至 ~88%（接近 Paired 的 ~90%），修正前所有 TO 數據作廢，7 個樣本全量重跑

![HP Fix 後 TO Tier 分佈](figures/O07_fig06_to_tier_distribution_post_fix.png)

---

## Slide 5 — 發現異常：LOH Enrichment 方向翻轉

**標題**：Paired 和 TO 的 LOH 含義完全相反 — 7/7 樣本一致

**核心訊息**：同一個 LOH 指標，Paired 說「FP 更 LOH」，TO 說「TP 更 LOH」——方向完全翻轉。

**要點**：

| 模式 | LOH Enrichment | 含義 | 機制 |
|------|:--------------:|------|------|
| **Paired** | **1.194x** | FP 更多 LOH-like | Paired FP 的 HP support 低 → 隨機波動產生假 LOH |
| **TO** | **0.805x** | **TP 更多 LOH-like** | TO phasing 中 somatic allele 造成系統性 HP 偏移 |

- 7/7 樣本全部一致：Paired = FP-enriched (1.02-3.18x)，TO = TP-enriched (0.852-0.956x)
- COLO829 因低深度（~30x）enrichment ≈ 1（接近中性），但方向仍為 TP-enriched
- **這直接導致 TO QS 失效**：QS 的 LOH penalty 假設「LOH = FP」，TO 下此假設不成立

![LOH Enrichment 熱圖 — 7 樣本 × 2 模式](figures/O03_fig04_loh_enrichment_heatmap.png)

> 此熱圖清楚展示各樣本的 LOH 比例差異與 Paired/TO 方向翻轉的對應關係

---

## Slide 6 — 量化異常：288K 同位點配對分析

**標題**：同一基因組位點，TO 單獨判 LOH 是 Paired 的 16-52 倍

**核心訊息**：288,609 個同位點 TP 配對分析 → TO-only LOH 是不一致的主體（95.5%）。

**要點**：

| 類別 | 數量 | 比例 |
|------|-----:|:----:|
| 共識：非 LOH | 157,774 | 55.0% |
| 共識：LOH | 81,738 | 30.6% |
| **TO-only LOH** | **39,724** | **13.8%** |
| Paired-only LOH | 1,856 | 0.6% |

- 一致率 85.5%（共識 LOH + 共識非 LOH）
- 不一致 41,580 位點中 **95.5% 是 TO 單獨判 LOH**
- Per-sample TO excess ratio = **16-52 倍**
- 共識 LOH 30.6% 代表在雙方都認定 LOH 的情況下，LOH 在所有數據中確實佔三成

![LOH Concordance 綜合分析 — 5 面板](figures/concordance_summary_5panel.png)

> (A) 全域 2×2 concordance matrix；(B) 每個樣本四分類比例；(C) 每個樣本 TO excess ratio（16-52×）；(D) Extreme LOH rate 方向翻轉；(E) 7 樣本 Enrichment 全部一致翻轉

---

## Slide 7 — 懷疑 TO 系統性偏差：所有 TO 特徵都無效

**標題**：9 項系統性觀察 → TO 模式全面 Negative（AUC < 0.58）

**核心訊息**：對 748,391 regions × 116 features 進行九項系統性觀察，TO 模式下所有特徵鑑別力不足，不只是 LOH。

**要點**：
- 9 項觀察（O1-O8 + O10）+ 4 項假說驗證（O11-O13 + N4），產出 82 張圖表
- **TO 全面 Negative**：全部 AUC < 0.58，5/9 甲基化特徵 Paired→TO 方向反轉
- **QS 完全失效**：Paired AUC=0.754 → TO AUC=0.497（等於隨機）
- LOH penalty 觸發率：TO TP **44.5%** vs FP **35.8%** → 反向懲罰 TP

| 特徵 | Paired AUC | TO AUC | 備註 |
|------|:----------:|:------:|------|
| GQ | **0.811** | — | TO 無 GQ |
| AF | 0.665 | 0.418 ⚠️ | 反轉 |
| ISM 甲基化 | 0.543 | 0.530 | 微弱 |
| LOH | 0.579 | 0.537 | 微弱 |
| **QS 綜合** | **0.754** | **0.497** ⚠️ | 隨機 |

![TO QS Waterfall](figures/O02_fig02_qs_component_waterfall_to.png)

---

## Slide 8 — 排除假說：不是 Read Depth 不足造成的

**標題**：TO-only LOH 的讀取深度反而更高 — 排除 Low-Depth 假說

**核心訊息**：TO-only LOH 位點的 effective_hp_reads 中位數=68，高於 both_LOH=60。同位點 Paired HP_Ratio=0.509（平衡），TO=0.026（極端）→ read depth 不是原因。

**要點**：
- **如果是 depth 不足造成**：TO-only LOH 的 reads 應顯著低於 both_LOH → **不成立**
- TO-only median eff_hp_reads = **68**（高於 both_LOH 的 60），Cohen's d = +0.29（反向）
- 決定性數據：同一位點
  - Paired HP_Ratio 中位數 = **0.509**（完全平衡）
  - TO HP_Ratio 中位數 = **0.026**（極端 LOH）
  - Cohen's d = **-1.20**（巨大效應量）
- 86.5% 的 TO-only LOH 位點在 Paired 下完全平衡（HP_Ratio 0.2-0.8）

![四象限 Effective HP Reads 分佈](figures/concordance_quadrant_depth.png)

> Box plot 清楚顯示 TO-only LOH（橘色）的讀取深度與 both_LOH（綠色）相當甚至更高

---

## Slide 9 — 確認偏差：HP Ratio 散點圖的直接可視化

**標題**：同一基因組位點，Paired 看到平衡，TO 被推到極端

**核心訊息**：HP Ratio 散點圖是 self-phasing circular dependency 的直接可視化證據。

**要點**：
- X 軸=Paired HP_Ratio，Y 軸=TO HP_Ratio
- **綠色（both_LOH）**：四角（雙方都判 LOH）
- **藍色（neither）**：中央（雙方都不判 LOH）
- **橘色（TO-only LOH）**：**水平帶** — Paired 在 0.3-0.7（平衡），TO 被推到 <0.1 或 >0.9
- 這是 LongPhase-TO self-phasing 的**最直觀證據**：同一位點的 reads，Paired phasing 看到對稱分配，TO phasing 將它們系統性推向一側

![HP Ratio Paired vs TO — 四象限上色](figures/concordance_hp_ratio_scatter.png)

![Read Depth vs HP Ratio（2D hexbin）](figures/concordance_depth_vs_hpratio.png)

> Hexbin 圖顯示：即使在高 depth（>100）區域，TO HP ratio 仍被推到極端值 — 進一步排除「低 depth 隨機偏移」假說

---

## Slide 10 — 補充驗證：LOH Read Depth 與 Copy Number 分析

**標題**：LOH ≈ 0.73x Read Depth，60% Copy-Neutral — 不是 Allelic Deletion

**核心訊息**：LOH 區域 reads 不是少一半（0.73x 而非 0.5x），大部分是 copy-neutral LOH，TO 在 copy-gain 區域有額外偏差。

**要點**：
- LOH 區域 read depth ≈ non-LOH 的 **0.73x**（非理論 0.5x）
- **60% copy-neutral**（可能是 UPD 或 phasing artifact），24-30% copy-loss，8-16% copy-gain
- HCC1954 最極端：ratio ≈ 1.0（LOH depth = non-LOH），已知大量基因組重排
- **TO copy-gain LOH 比例 = Paired 的 2 倍**（15.8% vs 8.4%）→ TO 在 copy-gain 區域有額外的系統性偏差

![Per-Sample LOH Read Depth Ratio](figures/loh_read_depth_ratio_per_sample.png)

> 此圖清楚展示各 sample 的 LOH depth 分佈：紅線 0.5 = 純 copy-loss 理論值，紅線 1.0 = copy-neutral。大多數 sample 落在 0.68-0.82 之間，遠高於 0.5

![CN LOH Composition Per Sample](figures/cn_loh_composition_per_sample.png)

---

## Slide 11 — 根因確認：Self-Phasing Circular Dependency

**標題**：LongPhase-TO 的 Somatic Variant 參與自身 Phasing — 造成循環依賴

**核心訊息**：TO 使用 somatic variant 做 phasing anchor，TP 的 ALT reads 驅動 HP 偏斜，形成循環依賴。低 read 區域更容易被串入長 phasing block，但 haplotag 品質遠不如 Paired。

**機制推導（四步邏輯鏈）**：
1. **TO joint phasing**：LongPhase-TO 將 germline + somatic variant 共同建構 phasing graph（`PhasingGraph.cpp` addEdge() 確認），但 somatic variant 只有 ALT allele 貢獻 haplotype voting
2. **Self-phasing circular dependency**：TP 的 somatic ALT reads 天然偏向一個 haplotype → 形成 haplotype block → HP tagging 繼承偏斜 → HP_Ratio 偏離 0.5
3. **Low-read 區域被串入長 block**：TO 在 read 數量低的區域也嘗試串接，造成 phasing block 特別長，但 haplotag 結果品質很差，不如高可信度的 Paired phasing
4. **TP 更容易被判 LOH-like**：TO TP LOH 觸發率 44.5% vs FP 35.8% → LOH 成為 TP 弱標記 → Enrichment 0.805（7/7 一致）

**量化驗證**：
- 39,724 位點從 Paired non-LOH → TO LOH-like
- 其中 71.6% 在 TO 的 min(HP) = 0（完全單 haplotype），同位點 Paired HP1:HP2 = 8:8
- Shapley 分解：TP rate 變動貢獻 105.6%，FP rate 僅 -5.4%

![LOH Penalty 觸發率對比](figures/O02_fig03_qs_loh_penalty_trigger_rate.png)

![HP Ratio 分佈（per sample × mode）](figures/O03_fig01_hp_ratio_distribution_by_sample_mode.png)

---

## Slide 12 — 否決的假說：甲基化方向全部關閉

**標題**：O11/O12/O13 — 三項甲基化假說全部否決，問題在 Phasing 不在 Methylation

**核心訊息**：嘗試用甲基化特徵解決 TO 問題，但三項假說全部因 confound 被否決。真正的問題在 phasing 層面。

**要點**：
- **O11**（Within-group heterogeneity）：epipolymorphism AUC 0.845 → n_reads 校正後 0.530 → read count confound
- **O12**（LOH 甲基化場景）：AlleleDelta = AF confound + L2 collider bias → 全 corrected AUC < 0.58
- **O13**（跨區域 correlation）：shared read count confound → 控制後差異消失
- **N4**（HP0 filter）：TP loss ≤2% 下 FP removal = 0%
- **教訓**：L1→L2→L3 三層驗證架構，任何 L2 有效的特徵必須通過 L3 AF-bin 交叉驗證

---

## Slide 13 — 已完成修正：QS Mode-Aware

**標題**：第一步修正 — 消除 TO LOH Penalty 的反效果

**核心訊息**：基於機制證據，TO 模式下停用 LOH penalty 和 verify bonus。AUC 從有害恢復到中性。

**要點**：

| QS 組件 | Paired | TO（修改後） | 理由 |
|---------|:------:|:----------:|------|
| LOH penalty | -25 | **0** | 方向反轉，懲罰 TP |
| Verify bonus | +15 | **0** | 依賴 LOH 判定 |
| 其他 | 不變 | 不變 | 方向一致 |

- **先前做法**：使用 Quality Score 加分表整合多組件，在 Paired 有效（AUC=0.754），TO 失效（AUC=0.497）
- **修正效果**：TO QS AUC 0.497 → ~0.546（理論推算），from harmful to neutral
- 已提交：commit `b9eaba7`
- **注意**：0.546 仍接近隨機，此修正意義是「消除反效果」，不是讓 QS 變有效

---

## Slide 14 — 修正方向：TO Phasing 程式碼改善可能性

**標題**：如何改正 TO 的 Phasing 品質 — 分析與修正方向

**核心訊息**：根因在 LongPhase-TO 的 phasing 策略，需分析是否有程式碼層面的改善空間。

**分析方向**：
1. **LongPhase-TO phasing anchor 策略**：目前將 somatic variant 納入 phasing graph，可否排除或降權？
2. **Phasing block 品質過濾**：TO 在低 read 區域串出的長 block 品質差，可否加 block quality threshold？
3. **HP tagging 後處理**：在 ISM 層面，可否根據 phasing block 大小 / read support 動態調整 LOH 判定門檻？
4. **Paired phasing 作為 benchmark**：利用 288K 同位點數據，量化 TO phasing 在哪些條件下可以信任

**已有的改善基礎**：
- 288,609 matched loci 數據可以做 TO phasing quality calibration
- 已知 86.5% 的 TO-only LOH 在 Paired 下平衡 → 這些是可以被校正的 false LOH
- Copy-gain 區域額外偏差（2x Paired）→ 需要 depth-aware 校正

**待執行**：
- 分析 LongPhase-TO 的 phasing block 長度分佈 vs Paired
- 比較 block 內 read count 和 haplotag 一致性
- 評估排除 somatic variant from phasing anchor 的可行性

---

## Slide 15 — 核心結論

**標題**：本週結論 — Paired/TO 是根本不同的問題空間

**核心訊息**：LOH 在 Paired 和 TO 的行為方向完全相反，必須分離建模，且 TO 需要 phasing 層面的根本修正。

**核心數字**：

| 維度 | Paired | TO |
|------|:------:|:--:|
| FP rate | 1% | 31% |
| LOH enrichment | FP-enriched (1.194x) | TP-enriched (0.805x) |
| 特徵方向 | — | 5/9 反轉 |
| QS AUC | 0.754 | 0.497 (隨機) |
| TO-only LOH excess | — | 16-52x over Paired-only |
| 自我驗證 | — | Cohen's d = -1.20 |

**結論層次**：
1. ✅ Paired/TO 分離不是選項，是**必要條件**
2. ✅ TO LOH 偏差源自 **self-phasing circular dependency**，非 read depth 或 copy number
3. ✅ 60% LOH 是 copy-neutral → 部分可能是 phasing artifact
4. ⚠️ 甲基化特徵無法解決 TO 問題（O11-O13 全否決）→ 問題在 phasing 層
5. 📋 下一步：TO phasing 程式碼分析與修正方案設計

**Limitations**：
- 所有結論基於 7 個高 purity 細胞系（purity ≈ 100%），臨床低 purity 需額外驗證
- LOH-like 閾值 (0.1/0.9) 為操作定義，需 sensitivity analysis

---

## Slide 16 — 討論與下週方向

**標題**：Discussion & Next Steps

**下週行動**：

| 優先級 | 行動 | 依據 |
|--------|------|------|
| **P0** | TO phasing block 品質分析 | 根因在 phasing |
| **P0** | Paired/TO 分離策略框架 | 方向完全相反 |
| **P1** | Paired ML 特徵集（GQ + AF + LOH subtype） | GQ AUC=0.811 |
| **P2** | TO phasing anchor 修正方案評估 | 排除 somatic variant |

**開放討論問題**：
1. LongPhase-TO 排除 somatic variant 後，phasing 品質是否足夠？
2. TO 獨立 QS 框架應以什麼為 anchor？
3. Paired phasing 可以做 TO 的 calibration reference 嗎？

---

## 圖片使用總覽

| Slide | 圖片 | 說明 |
|:-----:|------|------|
| 2 | O01_fig04_truth_label_composition.png | TP/FP 比例對比 |
| 4 | O07_fig06_to_tier_distribution_post_fix.png | HP fix 後 tier 分佈 |
| 5 | O03_fig04_loh_enrichment_heatmap.png | 7×2 enrichment 熱圖 |
| 6 | concordance_summary_5panel.png | 同位點 concordance 5 面板 |
| 7 | O02_fig02_qs_component_waterfall_to.png | TO QS waterfall |
| 8 | concordance_quadrant_depth.png | 四象限 read depth |
| 9 | concordance_hp_ratio_scatter.png | HP ratio 散點圖 |
| 9 | concordance_depth_vs_hpratio.png | Depth vs HP ratio hexbin |
| 10 | loh_read_depth_ratio_per_sample.png | Per-sample LOH depth ratio |
| 10 | cn_loh_composition_per_sample.png | CN composition per sample |
| 11 | O02_fig03_qs_loh_penalty_trigger_rate.png | LOH penalty 觸發率 |
| 11 | O03_fig01_hp_ratio_distribution_by_sample_mode.png | HP ratio 分佈 |
