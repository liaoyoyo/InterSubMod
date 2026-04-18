<!--
建立時間: 2026-04-02 04:00
目標: Phase 3 — 逐張投影片視覺佈局描述 + 口述講稿（v2 大修版，對齊偵探故事線）
處理範圍: 16 張投影片的版面設計與口頭講解稿
關聯檔案:
  - 02_ppt_slide_outline.md (v2 大修版段落主軸)
  - 01_full_narrative_report.md (完整敘事)
  - 05_loh_read_depth_and_copy_number_analysis.md (新數據)
-->

# PPT 投影片佈局說明與講解腳本（v2 大修版）

**版面規範**：16:9 寬螢幕，深色背景（#1a1a2e 或類似），白色主文字，強調色 #00d4ff（藍）/ #ff6b6b（紅）/ #4ecdc4（綠）
**字型**：標題 28-32pt 粗體，內文 18-22pt，數字強調 36-48pt
**通則**：每張投影片一個核心訊息，不超過 6 行文字，圖片佔版面 40-60%
**預估總時長**：20-25 分鐘（16 張投影片）
**敘事主線**：偵探故事 — 發現異常 → 懷疑偏差 → 驗證確認 → 追查根因 → 提出修正

---

## Slide 1 — 封面

**版面配置**：
- 上半（60%）：標題區塊置中，主標題 + 副標題
- 下半（40%）：左側放報告資訊（報告人/日期/資料規模），右側放一句話摘要（框線強調）

```
┌────────────────────────────────────────────────────┐
│                                                    │
│   TO 模式 LOH 系統性偏差調查：                       │
│       從發現到機制剖析                               │
│                                                    │
│   ────────────────────────────────────              │
│                                                    │
│   報告期間：2026-03-25 ~ 2026-03-31                 │
│   資料規模：748,391 regions × 116 features           │
│             × 7 cancer cell lines × 2 modes         │
│   [報告人]  [日期]                                   │
│                                                    │
│   ┌── 一句話摘要 ──────────────────────────────┐     │
│   │ 發現 TO 模式的 LOH 判定存在系統性偏差，完整    │     │
│   │ 驗證機制並確認是 LongPhase-TO self-phasing     │     │
│   │ 造成，提出修正方向                             │     │
│   └──────────────────────────────────────────┘     │
│                                                    │
└────────────────────────────────────────────────────┘
```

**講稿**（預估 1 分鐘）：
> 各位好，這是本週 InterSubMod 的研究週報。本週的主題是一個偵探故事——我們在擴展 Tumor-Only 模式的分析時，發現了一個看起來不對勁的現象，花了整整一週追查它的來源。最後確認這是 LongPhase-TO 在做 phasing 時引入的系統性偏差。我會用 15 張投影片帶大家走過這個調查過程。

**過渡句**：
> 首先，讓我說明為什麼我們要把分析擴展到 Tumor-Only 模式。

---

## Slide 2 — 為什麼要研究 LOH

**版面配置**：
- 上方（15%）：標題「從 Paired 延伸到 Tumor-Only — 為了統計檢定力」
- 中間主體（55%）：左右各一個對比卡片（Paired vs TO），數字加大強調
- 下方（30%）：圖片（TP/FP 比例堆疊條形圖）+ 右側三個核心問題

```
┌────────────────────────────────────────────────────┐
│  從 Paired 延伸到 Tumor-Only — 為了統計檢定力        │
│                                                    │
│  ┌─── Paired ──────────┐  ┌─── TO ──────────────┐  │
│  │  FP = 3,429 (1%)    │  │  FP = 128,382 (31%) │  │
│  │  TP:FP = 95:1       │  │  TP:FP = 2.3:1      │  │
│  │  [紅] 統計力不足     │  │  [綠] 統計力充足     │  │
│  │                     │  │                      │  │
│  │ FP 來源：           │  │ FP 來源：             │  │
│  │ calling 殘留（~1%） │  │ germline 誤判（~31%）│  │
│  └─────────────────────┘  └──────────────────────┘  │
│                                                    │
│  ┌── 圖 40% ──────────┐  本週核心問題：             │
│  │ O01_fig04_truth_    │  1. ISM 能否在 TO 區分     │
│  │ label_composition   │     TP/FP？                │
│  │      .png           │  2. LOH 扮演什麼角色？     │
│  └─────────────────────┘  3. Paired/TO 一致嗎？     │
│                                                    │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![TP/FP 比例對比](figures/O01_fig04_truth_label_composition.png)

**講稿**（預估 1.5 分鐘）：
> 先說背景。上週我們建立了 Paired 模式的基線，但遇到一個根本的問題——Paired 的 False Positive 只有百分之一，大約三千四百個，TP 對 FP 是 95 比 1。這麼懸殊的比例讓我們沒辦法可靠地評估特徵的區分能力。
>
> 相比之下，Tumor-Only 模式因為沒有正常樣本對照，大量 germline variants 被錯誤判定為 somatic，False Positive 率高達百分之三十一。[停頓] 這提供了十二萬八千多個 FP，統計力完全充足。
>
> 所以本週的策略很清楚：借 TO 的大量 FP 來深度驗證 ISM 的區分能力，特別是 LOH 這個組件在裡面扮演什麼角色。

**過渡句**：
> 在進入結果之前，先讓我快速說明 ISM 的分析框架和一個重要的技術細節——HP tag 的差異。

---

## Slide 3 — 技術背景：ISM 分析框架與 HP Tag 差異

**版面配置**：
- 上方（12%）：標題
- 左側（55%）：ISM 流程示意 + HP tag 差異對比表
- 右側（45%）：先前 QS 嘗試的加分表概念圖

```
┌────────────────────────────────────────────────────┐
│  InterSubMod 分析架構與 Paired/TO Haplotype 差異      │
│                                                    │
│  ┌─── 流程圖 ──────────────────────────┐           │
│  │ [BAM] → [Region ±5000bp] → [HP 分群] │           │
│  │  → [CpG 甲基化矩陣] → [距離計算]      │           │
│  │  → [PERMANOVA] → [QS 評分]           │           │
│  └──────────────────────────────────────┘           │
│                                                    │
│  HP tag 來自 LongPhase haplotagged BAM              │
│  ┌──────────────┬──────────────────────┐           │
│  │ Paired       │ HP:Z:1-1（字串）      │           │
│  │ (LongPhase-s)│ haplotype-block 格式  │           │
│  ├──────────────┼──────────────────────┤           │
│  │ TO           │ HP:i:11（整數編碼）    │           │
│  │(LongPhase-TO)│ joint phasing        │           │
│  └──────────────┴──────────────────────┘           │
│                                                    │
│  關鍵差異：TO phasing anchor 包含 somatic variant     │
│  LOH-like: HP_Ratio < 0.1 or > 0.9 (eff_hp ≥ 30)  │
│  先前嘗試：Quality Score 加分表（Penalties + Bonuses） │
│  → Paired AUC=0.754                                │
└────────────────────────────────────────────────────┘
```

**講稿**（預估 1.5 分鐘）：
> 快速過一下技術框架。ISM 針對每個 somatic variant，擷取周邊正負五千個 base pair 範圍的 reads，依照 haplotype 分群後比較甲基化模式。
>
> 這裡有一個重要的細節：reads 的 haplotype 歸屬不是 ISM 自己判定的，而是來自 LongPhase 做 haplotagging 之後 BAM 檔案裡的 HP tag。[停頓] 但 Paired 和 TO 的 HP tag 格式完全不同。Paired 用的是 LongPhase-s，HP tag 是字串格式，像是 HP:Z:1-1；TO 用的是 LongPhase-TO，HP tag 是整數格式，像是 HP:i:11。
>
> 另一個關鍵差異是：TO 的 phasing anchor 不只用 germline variants，還包含 somatic variant 本身——這造成一個潛在的循環依賴，我們待會會深入討論。
>
> 最後，我們先前已經嘗試過用 Quality Score 加分表做整合評分，在 Paired 模式拿到 AUC 0.754。

**過渡句**：
> 有了這個背景，接下來就是偵探故事的開始——我們拿到了第一批 TO 結果，但發現事情不太對。

---

## Slide 4 — 數據品質檢查：HP Tag Bug 發現與修正

**版面配置**：
- 上方（12%）：標題（紅色警示風格）
- 左側（50%）：四步故事線（帶圖示的時間線，由上到下）
- 右側（50%）：修正前後對比圖

```
┌────────────────────────────────────────────────────┐
│  ⚠ 結果看起來不對 — 追查到 HP 輸入格式問題           │
│                                                    │
│  ┌─── 故事線 ──────────┐  ┌─── 圖 ──────────────┐  │
│  │                     │  │                      │  │
│  │ 1. AI 執行 pipeline │  │  O07_fig06_to_tier_  │  │
│  │    → 程式無報錯     │  │  distribution_post_  │  │
│  │                     │  │  fix.png             │  │
│  │ 2. 檢視結果圖片     │  │                      │  │
│  │    → 發現 TO HP     │  │  [圖佔右半]           │  │
│  │      數據異常 ⚠     │  │                      │  │
│  │                     │  │                      │  │
│  │ 3. 追查根因：       │  │                      │  │
│  │    ReadParser.cpp   │  │                      │  │
│  │    只認字串格式     │  │                      │  │
│  │    HP:i 被靜默忽略  │  │                      │  │
│  │                     │  │                      │  │
│  │ 4. 修正 + 全量重跑  │  │                      │  │
│  │    TO Tier A/A+:    │  │                      │  │
│  │    ~50% → ~88%      │  │                      │  │
│  └─────────────────────┘  └──────────────────────┘  │
│                                                    │
│  修正：HP:i:11→"1-1", HP:i:21→"2-1", HP:i:33→"3"   │
│  7 個樣本全量重跑，修正前 TO 數據全部作廢             │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![HP Fix 後 TO Tier 分佈](figures/O07_fig06_to_tier_distribution_post_fix.png)

**講稿**（預估 1.5 分鐘）：
> 接下來是第一個轉折。[語氣轉為敘事] 我們讓 AI 自動執行整個 ISM pipeline，程式跑完了，沒有任何錯誤訊息。但當我們打開結果圖片來看的時候，發現 TO 的 HP 數據明顯不對——大量 reads 缺少 haplotype 歸屬。
>
> [停頓] 追查之後發現，ISM 的 ReadParser 只處理字串格式的 HP tag，也就是 Paired 用的 HP:Z。TO 的整數格式 HP:i 被靜默忽略了，完全沒有報錯。
>
> 修正之後效果很明顯：TO 的 Tier A 加 A+ 比例從大約百分之五十升到百分之八十八，接近 Paired 的百分之九十。修正前的 TO 數據全部作廢，七個樣本全量重跑。
>
> [重點強調] 這個 bug 的教訓是——程式不報錯不代表結果是對的。是看圖片的時候才發現的。

**過渡句**：
> 修正之後我們拿到了乾淨的數據，但真正的異常才剛剛開始浮現。

---

## Slide 5 — 發現異常：LOH Enrichment 方向翻轉

**版面配置**：
- 上方（12%）：標題（強調「完全相反」）
- 中間（45%）：左側放對比表格（Paired vs TO），右側放圖片
- 下方（43%）：三點說明 + 紅色警示框（「直接導致 TO QS 失效」）

```
┌────────────────────────────────────────────────────┐
│  Paired 和 TO 的 LOH 含義完全相反 — 7/7 樣本一致     │
│                                                    │
│  ┌─── 對比表 ────────────────┐  ┌── 圖 ────────┐   │
│  │                           │  │              │   │
│  │  模式    │ Enrichment     │  │ O03_fig04_   │   │
│  │  ────────┼──────────────  │  │ loh_enrich-  │   │
│  │  Paired  │ 1.194x → FP   │  │ ment_heat-   │   │
│  │          │ 更多 LOH-like  │  │ map.png      │   │
│  │  ────────┼──────────────  │  │              │   │
│  │  TO      │ 0.805x → TP   │  │ [圖佔右 50%] │   │
│  │          │ 更多 LOH-like  │  │              │   │
│  │  ────────┴──────────────  │  │              │   │
│  └───────────────────────────┘  └──────────────┘   │
│                                                    │
│  • 7/7 樣本全部一致：Paired FP-enriched (1.02-3.18x)│
│    TO TP-enriched (0.852-0.956x)                    │
│  • COLO829 低深度 (~30x) → enrichment ≈ 1           │
│                                                    │
│  ┌─ ⚠ ──────────────────────────────────────────┐  │
│  │ 直接後果：QS 的 LOH penalty 假設「LOH = FP」  │  │
│  │ → TO 下此假設不成立 → QS 失效                  │  │
│  └──────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![LOH Enrichment 熱圖 — 7 樣本 x 2 模式](figures/O03_fig04_loh_enrichment_heatmap.png)

**講稿**（預估 2 分鐘）：
> [語氣提升] 這一頁是整個調查的起點。當我們用修正後的數據計算 LOH enrichment——也就是 FP 的 LOH 比例除以 TP 的 LOH 比例——我們發現了一個令人意外的結果。
>
> 在 Paired 模式下，enrichment 是 1.194 倍，意思是 FP 比 TP 更容易被判為 LOH-like。這符合直覺，因為 FP 的 haplotype 支持度低，隨機波動容易產生假的 LOH 信號。
>
> [停頓] 但在 TO 模式下，enrichment 是 0.805 倍——方向完全翻轉了。這代表 TP 反而比 FP 更容易被判為 LOH-like。[重點強調] 而且這不是某個樣本的特例，七個樣本全部一致。
>
> 這直接導致一個嚴重的後果：我們的 Quality Score 裡有一個 LOH penalty，假設「LOH 代表可疑」。在 TO 下，這個懲罰反而在懲罰 TP。

**過渡句**：
> 方向翻轉聽起來嚇人，但到底有多嚴重？讓我們用同位點配對分析來量化。

---

## Slide 6 — 量化異常：288K 同位點配對分析

**版面配置**：
- 上方（12%）：標題
- 中間（50%）：左側放四分類表格，右側放重點數字（大字 + 強調色）
- 下方（38%）：五面板圖（全版寬）

```
┌────────────────────────────────────────────────────┐
│  同一基因組位點，TO 單獨判 LOH 是 Paired 的 16-52 倍  │
│                                                    │
│  ┌─── 表格 ──────────────┐  ┌── 重點數字 ────────┐ │
│  │                       │  │                    │ │
│  │ 共識：非 LOH │ 55.0%  │  │  一致率 85.5%      │ │
│  │ 共識：LOH    │ 30.6%  │  │                    │ │
│  │ TO-only LOH  │ 13.8%  │  │  不一致中           │ │
│  │ Paired-only  │  0.6%  │  │  95.5%             │ │
│  │              │        │  │  是 TO 單獨判       │ │
│  │ n = 288,609  │        │  │                    │ │
│  └───────────────────────┘  └────────────────────┘ │
│                                                    │
│  ┌── concordance_summary_5panel.png ─────────────┐  │
│  │                                               │  │
│  │  (A) 2×2 concordance  (B) per-sample 比例     │  │
│  │  (C) TO excess ratio  (D) LOH rate 翻轉      │  │
│  │  (E) Enrichment 一致                          │  │
│  │                                               │  │
│  └───────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![LOH Concordance 綜合分析 — 5 面板](figures/concordance_summary_5panel.png)

**講稿**（預估 2 分鐘）：
> 我們有二十八萬八千多個同時存在於 Paired 和 TO 的 TP 位點。用這些同位點配對來看 LOH 判定的一致性。
>
> 結果分成四類：百分之五十五雙方都說不是 LOH，百分之三十點六雙方都說是 LOH——這兩類是共識，加起來一致率百分之八十五點五。
>
> [停頓] 重點在不一致的部分。四萬一千多個不一致位點中，[重點強調] 百分之九十五點五是 TO 單獨判為 LOH，而 Paired 認為不是。反過來，Paired 單獨判 LOH 的只有百分之零點六。
>
> 看子圖 C，每個樣本的 TO excess ratio 在 16 到 52 倍之間。這不是小數目——三萬九千多個位點在 Paired 看起來完全正常，但 TO 把它們判成了 LOH。

**過渡句**：
> 那現在問題來了：這三萬九千多個 TO-only LOH 是怎麼回事？是因為 reads 太少，還是有更深層的原因？

---

## Slide 7 — 懷疑 TO 系統性偏差：所有 TO 特徵都無效

**版面配置**：
- 上方（12%）：標題
- 左側（55%）：AUC 對比表格（5 列），紅色標記反轉/失效
- 右側（45%）：TO QS Waterfall 圖

```
┌────────────────────────────────────────────────────┐
│  9 項系統性觀察 → TO 模式全面 Negative（AUC < 0.58） │
│                                                    │
│  ┌─── AUC 對比表 ───────┐  ┌── 圖 ──────────────┐ │
│  │                      │  │                     │ │
│  │ 特徵     │ Pa  │ TO  │  │  O02_fig02_qs_     │ │
│  │ ─────────┼─────┼──── │  │  component_water-  │ │
│  │ GQ       │.811 │ (無)│  │  fall_to.png       │ │
│  │ AF       │.665 │.418⚠│  │                     │ │
│  │ ISM 甲基化│.543 │.530 │  │  [圖佔右半]         │ │
│  │ LOH      │.579 │.537 │  │                     │ │
│  │ QS 綜合  │.754 │.497⚠│  │                     │ │
│  │          │     │     │  │                     │ │
│  └──────────────────────┘  └─────────────────────┘ │
│                                                    │
│  • 9 項觀察 + 4 項假說驗證 → 82 張圖表              │
│  • TO 全面 AUC < 0.58，5/9 方向反轉                 │
│  • LOH penalty 觸發：TO TP 44.5% vs FP 35.8%       │
│    → 反向懲罰 TP                                    │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![TO QS Waterfall](figures/O02_fig02_qs_component_waterfall_to.png)

**講稿**（預估 1.5 分鐘）：
> 問題不只是 LOH。我們做了九項系統性觀察加上四項假說驗證，總共產出八十二張圖表。結論是——TO 模式下所有特徵的 AUC 都低於 0.58，接近隨機。
>
> 看這張表：Paired 的 GQ 有 0.811，是最強的特徵，但 TO 根本沒有 GQ。AF 在 Paired 是 0.665，到了 TO 變成 0.418，方向翻轉了。最嚴重的是 QS 綜合評分：Paired 0.754，[重點強調] TO 只有 0.497，比丟銅板還不如。
>
> 看右邊的 waterfall 圖，TO 的 QS 各組件貢獻幾乎都是反向的。特別是 LOH penalty 的觸發率——TO 的 TP 觸發百分之四十四點五，FP 只有百分之三十五點八，等於在懲罰好的位點。

**過渡句**：
> 我們開始懷疑 TO 有系統性偏差。但第一個要排除的假說是：會不會只是 reads 太少造成的？

---

## Slide 8 — 排除假說：不是 Read Depth 不足造成的

**版面配置**：
- 上方（12%）：標題（綠色「排除」標記）
- 左側（45%）：三個關鍵數字（大字卡片風格）
- 右側（55%）：四象限 read depth 箱型圖

```
┌────────────────────────────────────────────────────┐
│  ✓ 排除 — TO-only LOH 的讀取深度反而更高              │
│                                                    │
│  ┌─── 關鍵數字 ──────────┐ ┌── 圖 ──────────────┐  │
│  │                       │ │                     │  │
│  │  TO-only LOH          │ │ concordance_        │  │
│  │  eff_hp_reads         │ │ quadrant_depth.png  │  │
│  │  中位數 = 68          │ │                     │  │
│  │  （高於 both_LOH=60） │ │ [圖佔右半]           │  │
│  │                       │ │                     │  │
│  │  同位點：              │ │                     │  │
│  │  Paired HP_Ratio      │ │                     │  │
│  │  = 0.509（平衡）      │ │                     │  │
│  │  TO HP_Ratio          │ │                     │  │
│  │  = 0.026（極端）      │ │                     │  │
│  │                       │ │                     │  │
│  │  Cohen's d = -1.20    │ │                     │  │
│  │  （巨大效應量）        │ │                     │  │
│  └───────────────────────┘ └─────────────────────┘  │
│                                                    │
│  86.5% 的 TO-only LOH 位點在 Paired 下完全平衡       │
│  （HP_Ratio 0.2-0.8）                               │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![四象限 Effective HP Reads 分佈](figures/concordance_quadrant_depth.png)

**講稿**（預估 1.5 分鐘）：
> 第一個要排除的假說：TO-only LOH 是不是因為 reads 太少，HP ratio 因為隨機波動跑到極端值？
>
> [停頓] 答案是否定的。TO-only LOH 的 effective HP reads 中位數是 68，反而高於雙方都判 LOH 的 60。看右邊的箱型圖，橘色（TO-only）和綠色（both LOH）的讀取深度相當，甚至更高。
>
> [重點強調] 決定性的證據是同位點比較：同一個基因組位置，Paired 看到的 HP ratio 中位數是 0.509——幾乎完全平衡；但 TO 看到的是 0.026——極端的 LOH。Cohen's d 是負 1.20，這是巨大的效應量。
>
> 更具體地說，百分之八十六點五的 TO-only LOH 位點，在 Paired 模式下 HP ratio 都在 0.2 到 0.8 之間，完全平衡。Reads 的數量夠，問題在於 reads 被分配到哪邊。

**過渡句**：
> 如果不是 depth 的問題，那讓我們直接看 HP ratio 的分佈——一張散點圖就能看出答案。

---

## Slide 9 — 確認偏差：HP Ratio 散點圖的直接可視化

**版面配置**：
- 上方（12%）：標題
- 左側主圖（55%）：HP ratio 散點圖（大）
- 右側輔圖（45%）：上方放散點圖解讀圖例，下方放 hexbin 圖

```
┌────────────────────────────────────────────────────┐
│  同一基因組位點，Paired 看到平衡，TO 被推到極端        │
│                                                    │
│  ┌─── 主圖 ──────────────┐  ┌── 解讀 ───────────┐  │
│  │                       │  │                    │  │
│  │ concordance_hp_ratio_ │  │ 綠 = both LOH     │  │
│  │ scatter.png           │  │   （四角）          │  │
│  │                       │  │ 藍 = neither       │  │
│  │ X = Paired HP_Ratio   │  │   （中央）          │  │
│  │ Y = TO HP_Ratio       │  │ 橘 = TO-only LOH  │  │
│  │                       │  │   （水平帶）⚠      │  │
│  │ [圖佔左 55%]          │  │                    │  │
│  │                       │  ├────────────────────┤  │
│  │                       │  │ concordance_depth_ │  │
│  │                       │  │ vs_hpratio.png     │  │
│  │                       │  │                    │  │
│  │                       │  │ 即使 depth>100     │  │
│  │                       │  │ TO 仍被推到極端    │  │
│  └───────────────────────┘  └────────────────────┘  │
│                                                    │
│  → Self-phasing 的最直觀證據                         │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![HP Ratio Paired vs TO — 四象限上色](figures/concordance_hp_ratio_scatter.png)
- ![Read Depth vs HP Ratio（2D hexbin）](figures/concordance_depth_vs_hpratio.png)

**講稿**（預估 2 分鐘）：
> [語氣強調] 這張圖是整個調查裡最有說服力的一張。X 軸是 Paired 的 HP ratio，Y 軸是 TO 的 HP ratio。每個點是一個同位點配對。
>
> 先看顏色。綠色是雙方都判 LOH，集中在四個角落——雙方一致認為不平衡。藍色是雙方都不判 LOH，聚集在中央——雙方一致認為平衡。
>
> [停頓] 現在看橘色——TO-only LOH。它們形成兩條水平帶。[重點強調] X 軸在 0.3 到 0.7，表示 Paired 看到的是平衡的；但 Y 軸在小於 0.1 或大於 0.9，表示 TO 把它們推到了極端值。同一批 reads，同一個基因組位置，因為 phasing 方式不同，haplotype 的分配完全不一樣。
>
> 右下方的 hexbin 圖進一步確認：即使在 depth 超過一百的區域，TO 的 HP ratio 仍然被推到極端。這徹底排除了低深度假說。

**過渡句**：
> 散點圖是視覺化證據，接下來我們從 copy number 的角度做補充驗證。

---

## Slide 10 — 補充驗證：LOH Read Depth 與 Copy Number 分析

**版面配置**：
- 上方（12%）：標題
- 左圖（50%）：per-sample LOH read depth ratio
- 右圖（50%）：CN composition per sample
- 下方（20%）：四個重點 bullet

```
┌────────────────────────────────────────────────────┐
│  LOH ≈ 0.73x Read Depth，60% Copy-Neutral           │
│  — 不是 Allelic Deletion                             │
│                                                    │
│  ┌── 左圖 ──────────────┐  ┌── 右圖 ──────────┐    │
│  │                      │  │                   │    │
│  │ loh_read_depth_      │  │ cn_loh_compo-     │    │
│  │ ratio_per_sample.png │  │ sition_per_       │    │
│  │                      │  │ sample.png        │    │
│  │ 紅線 0.5 = copy-loss │  │                   │    │
│  │ 紅線 1.0 = neutral   │  │ 灰=neutral 60%   │    │
│  │ 大多 0.68-0.82       │  │ 紅=loss 24-30%   │    │
│  │ HCC1954 ≈ 1.0       │  │ 藍=gain 8-16%    │    │
│  │                      │  │                   │    │
│  └──────────────────────┘  └───────────────────┘    │
│                                                    │
│  • LOH 區域 read depth ≈ non-LOH 的 0.73x          │
│  • 60% copy-neutral → UPD 或 phasing artifact      │
│  • HCC1954 ratio ≈ 1.0 → 已知大量基因組重排          │
│  • TO copy-gain LOH = Paired 的 2 倍（15.8% vs 8.4%）│
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![Per-Sample LOH Read Depth Ratio](figures/loh_read_depth_ratio_per_sample.png)
- ![CN LOH Composition Per Sample](figures/cn_loh_composition_per_sample.png)

**講稿**（預估 1.5 分鐘）：
> 這一頁從 copy number 的角度做補充。如果 LOH 是真正的一個 allele 被刪除——也就是 copy-loss LOH——read depth 應該只有正常的一半，也就是 0.5 倍。但實際上，LOH 區域的 depth 是 non-LOH 的 0.73 倍，遠高於 0.5。
>
> 右邊的圖更清楚：大約百分之六十的 LOH 是 copy-neutral，也就是讀取深度正常但 haplotype 不平衡。這可能是 uniparental disomy，也可能是 phasing artifact。只有百分之二十四到三十是真正的 copy-loss。
>
> [停頓] 值得注意的是，TO 的 copy-gain LOH 比例是 Paired 的兩倍——百分之十五點八對百分之八點四——表示 TO 在 copy-gain 區域有額外的系統性偏差。

**過渡句**：
> 所有的線索都指向同一個方向。現在讓我們把機制完整拼出來。

---

## Slide 11 — 根因確認：Self-Phasing Circular Dependency

**版面配置**：
- 上方（12%）：標題
- 左側（55%）：四步邏輯鏈（垂直流程圖，帶箭頭）
- 右上（22%）：LOH Penalty 觸發率圖
- 右下（22%）：HP Ratio 分佈圖

```
┌────────────────────────────────────────────────────┐
│  LongPhase-TO 的 Somatic Variant 參與自身 Phasing    │
│  — 造成循環依賴                                      │
│                                                    │
│  ┌─── 四步邏輯鏈 ──────────┐  ┌── 圖 1 ─────────┐  │
│  │                         │  │                  │  │
│  │ ❶ TO joint phasing      │  │ O02_fig03_qs_   │  │
│  │   germline + somatic    │  │ loh_penalty_    │  │
│  │   共建 phasing graph    │  │ trigger_rate    │  │
│  │           ↓             │  │    .png         │  │
│  │ ❷ Self-phasing          │  │                  │  │
│  │   somatic ALT reads     │  ├──────────────────┤  │
│  │   偏向一個 haplotype    │  │                  │  │
│  │           ↓             │  │ O03_fig01_hp_   │  │
│  │ ❸ Low-read 區域         │  │ ratio_distri-   │  │
│  │   被串入長 block        │  │ bution_by_      │  │
│  │   品質差               │  │ sample_mode     │  │
│  │           ↓             │  │    .png         │  │
│  │ ❹ TP 更容易判 LOH-like  │  │                  │  │
│  │   44.5% vs 35.8%       │  │                  │  │
│  └─────────────────────────┘  └──────────────────┘  │
│                                                    │
│  量化：39,724 位點中 71.6% TO min(HP)=0（完全單側）   │
│  同位點 Paired HP1:HP2 = 8:8（完全平衡）              │
│  Shapley 分解：TP rate 貢獻 105.6%，FP rate 僅 -5.4% │
└────────────────────────────────────────────────────┘
```

**圖片**：
- ![LOH Penalty 觸發率對比](figures/O02_fig03_qs_loh_penalty_trigger_rate.png)
- ![HP Ratio 分佈（per sample x mode）](figures/O03_fig01_hp_ratio_distribution_by_sample_mode.png)

**講稿**（預估 2 分鐘）：
> [語氣嚴肅] 現在把機制完整串起來。這是一個四步的邏輯鏈。
>
> 第一步：LongPhase-TO 在做 phasing 的時候，把 germline 和 somatic variant 一起丟進 phasing graph。我們從 LongPhase 的原始碼 PhasingGraph.cpp 確認了這件事。
>
> 第二步：因為 TP 的 somatic variant 有真實的 ALT allele reads，這些 reads 天然偏向一個 haplotype。[停頓] 也就是說，variant 自己的 reads 在影響自己被歸到哪個 haplotype——這是一個循環依賴。
>
> 第三步：在 read 數量少的區域，TO 仍然嘗試把它串接到 phasing block 裡面，但品質遠不如 Paired。
>
> 第四步：結果就是 TP 更容易被判成 LOH-like——觸發率百分之四十四點五，比 FP 的百分之三十五點八還高。
>
> [重點強調] 定量上：三萬九千多個 TO-only LOH 位點中，百分之七十一點六在 TO 下某一側 HP 的 read 數是零——完全單側。但同一個位點，Paired 看到的 HP1 比 HP2 是 8 比 8，完美平衡。

**過渡句**：
> 在確認機制之前，我們也嘗試過用甲基化特徵來解決這個問題——但全部被否決了。

---

## Slide 12 — 否決的假說：甲基化方向全部關閉

**版面配置**：
- 上方（12%）：標題（灰色「否決」風格）
- 主體（70%）：四個假說卡片（2x2 排列），每張標紅色叉
- 下方（18%）：教訓框

```
┌────────────────────────────────────────────────────┐
│  O11/O12/O13 — 三項甲基化假說全部否決                 │
│  問題在 Phasing 不在 Methylation                     │
│                                                    │
│  ┌─── O11 ──────────────┐  ┌─── O12 ──────────┐    │
│  │ ✗ Within-group       │  │ ✗ LOH 甲基化場景  │    │
│  │   heterogeneity      │  │                   │    │
│  │                      │  │ AlleleDelta =     │    │
│  │ epipolymorphism      │  │ AF confound       │    │
│  │ AUC 0.845 → 0.530   │  │ + L2 collider     │    │
│  │ (n_reads confound)   │  │ bias              │    │
│  │                      │  │ 全 corrected      │    │
│  │                      │  │ AUC < 0.58        │    │
│  └──────────────────────┘  └───────────────────┘    │
│                                                    │
│  ┌─── O13 ──────────────┐  ┌─── N4 ───────────┐    │
│  │ ✗ 跨區域 correlation │  │ ✗ HP0 filter     │    │
│  │                      │  │                   │    │
│  │ shared read count    │  │ TP loss ≤2% 下    │    │
│  │ confound             │  │ FP removal = 0%   │    │
│  │ 控制後差異消失        │  │                   │    │
│  └──────────────────────┘  └───────────────────┘    │
│                                                    │
│  教訓：L1→L2→L3 三層驗證，L2 有效必須通過 L3 AF-bin   │
│  交叉驗證。表面顯著的特徵可能只是 confound。           │
└────────────────────────────────────────────────────┘
```

**講稿**（預估 1.5 分鐘）：
> 在確認機制之前，我們其實花了不少時間試圖用甲基化特徵來解決 TO 的問題。結果三項假說全部被否決。
>
> O11 嘗試用 within-group methylation heterogeneity 來區分，表面上 epipolymorphism 的 AUC 有 0.845，看起來很強。但控制 read count 之後掉到 0.530——原來是 read 數量的 confound。
>
> O12 嘗試用 LOH 區域的甲基化場景分類，結果 AlleleDelta 本身就等同於 AF，加上我們發現了 L2 collider bias，校正後全部 AUC 低於 0.58。
>
> O13 嘗試跨區域的甲基化 correlation，看起來 FP 比 TP 的 correlation 高，但控制 shared read count 之後差異完全消失。
>
> [停頓] 還有 N4 — HP0 filter，在容許 TP loss 百分之二以下的條件，FP removal 是零。
>
> 結論很明確：問題在 phasing 層面，不在 methylation。

**過渡句**：
> 既然原因確認了，我們已經做了第一步修正。

---

## Slide 13 — 已完成修正：QS Mode-Aware

**版面配置**：
- 上方（12%）：標題
- 中間（50%）：左側放修改前後對比表，右側放 QS AUC 變化圖示（大數字）
- 下方（38%）：三點說明 + commit 資訊

```
┌────────────────────────────────────────────────────┐
│  第一步修正 — 消除 TO LOH Penalty 的反效果            │
│                                                    │
│  ┌─── 修改對比 ──────────────────────────────────┐  │
│  │                                               │  │
│  │  QS 組件        │ Paired  │ TO（修改後）       │  │
│  │  ────────────────┼─────────┼──────────────────  │  │
│  │  LOH penalty    │  -25    │  0  ← 方向反轉     │  │
│  │  Verify bonus   │  +15    │  0  ← 依賴 LOH     │  │
│  │  其他           │  不變    │  不變               │  │
│  │                                               │  │
│  └───────────────────────────────────────────────┘  │
│                                                    │
│  ┌── AUC 變化 ──────────────────────────────────┐   │
│  │                                              │   │
│  │  先前做法：Quality Score 加分表               │   │
│  │  Paired AUC = 0.754 ✓                        │   │
│  │  TO AUC = 0.497（隨機）→ ~0.546（中性）       │   │
│  │                                              │   │
│  │  意義：消除反效果，不是讓 QS 變有效            │   │
│  │                                              │   │
│  └──────────────────────────────────────────────┘   │
│                                                    │
│  已提交：commit b9eaba7                              │
│  注意：0.546 仍接近隨機，後續需要全新的 TO 策略       │
└────────────────────────────────────────────────────┘
```

**講稿**（預估 1 分鐘）：
> 基於以上所有證據，我們做了第一步修正：在 TO 模式下停用 LOH penalty 和 verify bonus。
>
> 原本的做法是用一張 Quality Score 加分表，整合多個組件的懲罰和加分。在 Paired 模式效果不錯，AUC 0.754。但在 TO，LOH penalty 的方向完全反了，所以 TO QS 降到 0.497——比隨機還差。
>
> 停用之後，TO QS 理論上回到 0.546 左右。[停頓] 但我要強調，這個修正的意義是「不再傷害」，而不是「變好」。0.546 還是接近隨機。TO 需要的是全新的策略，不是修修補補。

**過渡句**：
> 那接下來怎麼辦？讓我分享一下修正方向的思考。

---

## Slide 14 — 修正方向：TO Phasing 程式碼改善可能性

**版面配置**：
- 上方（12%）：標題
- 主體（65%）：四個分析方向（編號卡片，垂直排列）
- 下方（23%）：已有的改善基礎（三個 bullet）

```
┌────────────────────────────────────────────────────┐
│  如何改正 TO 的 Phasing 品質 — 分析與修正方向         │
│                                                    │
│  ┌─── 四個方向 ─────────────────────────────────┐   │
│  │                                              │   │
│  │ ❶ Phasing Anchor 策略                        │   │
│  │   排除或降權 somatic variant 在 phasing graph │   │
│  │                                              │   │
│  │ ❷ Phasing Block 品質過濾                      │   │
│  │   低 read 區域的長 block 品質差 → 加門檻       │   │
│  │                                              │   │
│  │ ❸ HP Tagging 後處理                           │   │
│  │   根據 block 大小/read support 動態調整門檻    │   │
│  │                                              │   │
│  │ ❹ Paired Phasing 作為 Benchmark               │   │
│  │   288K 同位點做 TO phasing quality calibration │   │
│  │                                              │   │
│  └──────────────────────────────────────────────┘   │
│                                                    │
│  已有改善基礎：                                      │
│  • 288,609 matched loci → TO phasing calibration    │
│  • 86.5% TO-only LOH 在 Paired 下平衡 → 可校正      │
│  • Copy-gain 區域 2x 偏差 → 需 depth-aware 校正     │
│                                                    │
│  待執行：block 長度分佈 / haplotag 一致性 / 排除      │
│  somatic anchor 的可行性評估                         │
└────────────────────────────────────────────────────┘
```

**講稿**（預估 1.5 分鐘）：
> 根因在 LongPhase-TO 的 phasing 策略，所以修正方向也要從 phasing 著手。我們列出四個分析方向。
>
> 第一，分析能不能把 somatic variant 從 phasing anchor 中排除或降權。第二，TO 在低 read 區域串出的長 phasing block 品質很差，能不能加一個品質門檻。第三，在 ISM 這一層，能不能根據 phasing block 的資訊動態調整 LOH 判定的門檻。第四，利用二十八萬八千個同位點數據，用 Paired phasing 作為 benchmark 來校正 TO。
>
> [停頓] 好消息是，我們已經有很好的基礎。百分之八十六點五的 TO-only LOH 在 Paired 下是平衡的，這代表這些是理論上可以被校正的 false LOH。Copy-gain 區域有額外的偏差，需要結合 depth 資訊做校正。

**過渡句**：
> 最後讓我把本週的結論做一個完整的整理。

---

## Slide 15 — 核心結論

**版面配置**：
- 上方（12%）：標題
- 左側（50%）：核心數字對比表（Paired vs TO，6 列）
- 右側（50%）：五層結論（帶 icon 的列表，逐層遞進）
- 下方（15%）：Limitations 小字

```
┌────────────────────────────────────────────────────┐
│  本週結論 — Paired/TO 是根本不同的問題空間            │
│                                                    │
│  ┌─── 核心數字 ─────────┐  ┌── 結論層次 ─────────┐ │
│  │                      │  │                     │ │
│  │ 維度        Pa   TO  │  │ 1. Paired/TO 分離   │ │
│  │ ──────── ───── ───── │  │    是必要條件        │ │
│  │ FP rate   1%   31%  │  │                     │ │
│  │ LOH enr. 1.19x 0.81x│  │ 2. Self-phasing     │ │
│  │ 方向反轉   —  5/9反轉│  │    circular dep.    │ │
│  │ QS AUC  0.754 0.497 │  │    是根因            │ │
│  │ TO-only   —  16-52x │  │                     │ │
│  │ Cohen's d  — -1.20  │  │ 3. 60% copy-neutral │ │
│  │                      │  │    → phasing artifact│ │
│  │                      │  │                     │ │
│  │                      │  │ 4. 甲基化無法解決    │ │
│  │                      │  │    O11-O13 全否決   │ │
│  │                      │  │                     │ │
│  │                      │  │ 5. 下一步：TO       │ │
│  │                      │  │    phasing 分析修正  │ │
│  └──────────────────────┘  └─────────────────────┘ │
│                                                    │
│  Limitations: 7 個高 purity 細胞系（~100%），臨床低   │
│  purity 需額外驗證；LOH 閾值 0.1/0.9 需 sensitivity  │
│  analysis                                           │
└────────────────────────────────────────────────────┘
```

**講稿**（預估 2 分鐘）：
> [語氣總結] 最後做一個完整的整理。這張表列出 Paired 和 TO 在六個維度的對比——每一項都指向同一個結論：這兩個模式是根本不同的問題空間。
>
> 五個層次的結論。第一，[重點強調] Paired 和 TO 分離建模不是選項，是必要條件。第二，TO LOH 偏差的根因是 self-phasing circular dependency——LongPhase-TO 讓 somatic variant 參與自己的 phasing，而不是 read depth 或 copy number 的問題。第三，百分之六十的 LOH 是 copy-neutral，其中一部分可能是 phasing artifact，這對 Paired 模式也值得關注。
>
> [停頓] 第四，甲基化特徵無法解決 TO 的問題——O11、O12、O13 三項假說全部被否決。問題在 phasing 層面。第五，下一步是分析 TO 的 phasing 程式碼，設計修正方案。
>
> 底下列了兩個 limitations：所有結論基於七個高 purity 細胞系，臨床低 purity 需要額外驗證；LOH 的閾值 0.1 和 0.9 是操作定義，需要做 sensitivity analysis。

**過渡句**：
> 最後一頁，下週的行動計畫和想跟大家討論的問題。

---

## Slide 16 — 討論與下週方向

**版面配置**：
- 上方（12%）：標題
- 左側（55%）：下週行動表格（4 列，帶優先級顏色）
- 右側（45%）：三個開放討論問題（帶編號框）

```
┌────────────────────────────────────────────────────┐
│  Discussion & Next Steps                            │
│                                                    │
│  ┌─── 下週行動 ────────────────────────────────┐    │
│  │                                             │    │
│  │ P0 │ TO phasing block 品質分析              │    │
│  │    │ → 根因在 phasing                      │    │
│  │ ───┼───────────────────────────────────────  │    │
│  │ P0 │ Paired/TO 分離策略框架                 │    │
│  │    │ → 方向完全相反                         │    │
│  │ ───┼───────────────────────────────────────  │    │
│  │ P1 │ Paired ML 特徵集                       │    │
│  │    │ GQ + AF + LOH subtype                  │    │
│  │ ───┼───────────────────────────────────────  │    │
│  │ P2 │ TO phasing anchor 修正方案評估          │    │
│  │    │ 排除 somatic variant                   │    │
│  └─────────────────────────────────────────────┘    │
│                                                    │
│  ┌─── 開放討論 ─────────────────────────────────┐   │
│  │                                              │   │
│  │ 1. LongPhase-TO 排除 somatic variant 後，    │   │
│  │    phasing 品質是否足夠？                     │   │
│  │                                              │   │
│  │ 2. TO 獨立 QS 框架應以什麼為 anchor？        │   │
│  │                                              │   │
│  │ 3. Paired phasing 可以做 TO 的 calibration   │   │
│  │    reference 嗎？                            │   │
│  │                                              │   │
│  └──────────────────────────────────────────────┘   │
│                                                    │
│  謝謝！歡迎提問。                                    │
└────────────────────────────────────────────────────┘
```

**講稿**（預估 1.5 分鐘）：
> 下週的行動按優先級排。P0 有兩件：第一是分析 TO 的 phasing block 品質，比較 block 長度和 haplotag 一致性，因為根因在 phasing。第二是開始設計 Paired 和 TO 的分離策略框架，因為方向完全相反，不能用同一套規則。
>
> P1 是 Paired 的 ML 特徵集，GQ 的 AUC 有 0.811 是最強特徵，加上 AF 和 LOH subtype 應該能建出不錯的 Paired 模型。P2 是評估能否從 LongPhase-TO 排除 somatic variant 做 phasing anchor，這需要深入分析原始碼。
>
> 最後三個問題想跟大家討論。第一，排除 somatic variant 之後 TO 的 phasing 品質還夠不夠？第二，TO 獨立的評分框架應該以什麼為基準？第三，Paired 的 phasing 結果能不能作為 TO 的校正參考？
>
> 以上就是本週的報告，謝謝。歡迎提問。

---

## 圖片使用總覽

| Slide | 圖片 | 說明 |
|:-----:|------|------|
| 2 | O01_fig04_truth_label_composition.png | TP/FP 比例對比 |
| 4 | O07_fig06_to_tier_distribution_post_fix.png | HP fix 後 tier 分佈 |
| 5 | O03_fig04_loh_enrichment_heatmap.png | 7x2 enrichment 熱圖 |
| 6 | concordance_summary_5panel.png | 同位點 concordance 5 面板 |
| 7 | O02_fig02_qs_component_waterfall_to.png | TO QS waterfall |
| 8 | concordance_quadrant_depth.png | 四象限 read depth |
| 9 | concordance_hp_ratio_scatter.png | HP ratio 散點圖 |
| 9 | concordance_depth_vs_hpratio.png | Depth vs HP ratio hexbin |
| 10 | loh_read_depth_ratio_per_sample.png | Per-sample LOH depth ratio |
| 10 | cn_loh_composition_per_sample.png | CN composition per sample |
| 11 | O02_fig03_qs_loh_penalty_trigger_rate.png | LOH penalty 觸發率 |
| 11 | O03_fig01_hp_ratio_distribution_by_sample_mode.png | HP ratio 分佈 |

**未使用圖片（備用）**：

| 圖片 | 可用於 |
|------|--------|
| cn_chr_depth_heatmap.png | Slide 10 補充染色體層級 CN |
| cn_hp_balance_by_category.png | Slide 10 補充 HP balance per CN |
| cn_loh_rate_vs_depth.png | Slide 10 補充 LOH rate vs depth 相關 |
| concordance_depth_diff.png | Slide 8 補充 depth 差異分佈 |
| concordance_per_sample_depth.png | Slide 8 補充 per-sample 分佈 |
| hp_ratio_core_distribution_post_fix.png | Slide 4 補充修正後分佈 |
| loh_chr_depth_heatmap.png | Slide 10 補充染色體 depth |
| loh_hp_balance_scatter.png | Slide 9 補充 HP1 vs HP2 |
| loh_read_depth_distribution.png | Slide 10 補充 violin plot |
| O02_fig01_qs_component_waterfall_paired.png | Slide 7 補充 Paired waterfall |
| O03_fig03_core_loh_like_rate_by_tier.png | Slide 5 補充 tier-wise LOH rate |
| O03_fig06_to_vs_paired_loh_concordance.png | Slide 6 補充 concordance |
| O03_fig10_loh_feature_importance_for_tp_fp.png | Slide 7 補充 feature importance |
| O04_fig04_caller_feature_auc_by_mode.png | Slide 7 補充 caller feature AUC |
| O08_fig09_hcc1395_chr8_hotspot.png | Appendix 用 |
| O10_fig01_read_methyl_mean_by_alt_support.png | Appendix 用 |
| O10_fig05_read_methyl_by_truth_label.png | Appendix 用 |

---

## 時間分配總覽

| 段落 | Slides | 時間 | 敘事角色 |
|------|:------:|:----:|---------|
| 開場 | 1-2 | 2.5 min | 設定動機 |
| 技術與 Bug | 3-4 | 3 min | 建立信任 |
| 發現與量化 | 5-6 | 4 min | 呈現異常 |
| 診斷 | 7-9 | 5 min | 追查偏差 |
| 驗證 | 10-11 | 3.5 min | 確認根因 |
| 否決與修正 | 12-13 | 2.5 min | 關閉錯誤方向 |
| 結論與展望 | 14-16 | 5 min | 收束方向 |
| **總計** | **16** | **~25 min** | |
