# LOH_weekly_report_v5.pptx

- Slides: 25
- Date: N/A

## Slide 1: TO 模式 LOH 系統性偏差調查：
從發現到機制剖析

TO 模式 LOH 系統性偏差調查：
從發現到機制剖析

Systematic Investigation of LOH Bias in TO Mode: From Discovery to Mechanism

InterSubMod LOH 研究週報 — 偵探故事線

報告人：廖子游 | 日期：2026-04-01
報告期間：2026-03-25 ~ 2026-03-31
資料規模：748,391 regions × 116 features × 7 cancer cell lines × 2 modes

發現 TO 模式的 LOH 判定存在系統性偏差，完整驗證機制並確認是 LongPhase-TO self-phasing 造成，提出修正方向

1 / 24

使用AI後

> **Notes:** 各位好，這是本週 InterSubMod 的研究週報。本週的主題是一個偵探故事——我們在擴展 Tumor-Only 模式的分析時，發現了一個看起來不對勁的現象，花了整整一週追查它的來源。最後確認這是 LongPhase-TO 在做 phasing 時引入的系統性偏差。我會用二十三張投影片帶大家走過這個調查過程，分成三幕：發現異常、深入追查、根因與修正。

[過渡] 首先，讓我說明為什麼我們要把分析擴展到 Tumor-Only 模式。

---

## Slide 2: 從 Paired 延伸到 TO — 為了統計檢定力

從 Paired 延伸到 TO — 為了統計檢定力

From Paired to TO — Boosting Statistical Power

本週分析對象：純腫瘤樣本 (Tumor-Only, TO)
This week's focus: Tumor-Only (TO) samples
TO: 僅有腫瘤 BAM，無配對正常樣本
TO: Only tumor BAM, no matched normal sample
Paired FP = 3,429 (~1%) → 統計力不足
Paired FP = 3,429 (~1%) → Insufficient statistical power
TO FP = 128,382 (~31%) → 統計力充足
TO FP = 128,382 (~31%) → Sufficient statistical power
核心問題：ISM 在 TO 能否區分 TP/FP？
Key question: Can ISM distinguish TP/FP in TO mode?

2 / 24

> **Notes:** 先說背景。上週我們建立了 Paired 模式的基線，但遇到一個根本的問題——Paired 的 False Positive 只有百分之一，大約三千四百個，TP 對 FP 是 95 比 1。這麼懸殊的比例讓我們沒辦法可靠地評估特徵的區分能力。

相比之下，Tumor-Only 模式因為沒有正常樣本對照，大量 germline variants 被錯誤判定為 somatic，False Positive 率高達百分之三十一。這提供了十二萬八千多個 FP，統計力完全充足。

所以本週的策略很清楚：借 TO 的大量 FP 來深度驗證 ISM 的區分能力，特別是 LOH 這個組件在裡面扮演什麼角色。

[過渡] 在進入結果之前，先讓我快速說明 ISM 的分析流程。

---

## Slide 3: InterSubMod 五步分析流程

InterSubMod 五步分析流程

InterSubMod Five-Step Analysis Pipeline

BAM 讀取
(HP tag)

→

Region 擷取
(±5000bp)

→

HP 分群
(Haplotype)

→

CpG 甲基化
距離計算

→

統計檢驗
(PERMANOVA + QS)

HP tag 來自 LongPhase haplotagged BAM，非 ISM 自行判定 → phasing 品質直接決定下游結果

3 / 24

> **Notes:** 這是 ISM 的五步分析流程。第一步從 BAM 檔案讀取 reads 和它們的 HP tag——也就是 haplotype 歸屬。第二步擷取每個 somatic variant 周邊正負五千 base pair 的區域。第三步依照 HP tag 把 reads 分成兩組。第四步計算每組內的 CpG 甲基化距離矩陣。第五步做統計檢驗——PERMANOVA 和 Quality Score。

重點在這裡：HP tag 不是 ISM 自己判定的，而是來自 LongPhase 做 haplotagging 後的結果。所以 phasing 的品質直接決定了下游所有分析的可靠性。

[過渡] 接下來看 Paired 和 TO 在 phasing 上的關鍵差異。

---

## Slide 4: Paired vs TO Phasing 差異 — 關鍵技術細節

Paired vs TO Phasing 差異 — 關鍵技術細節

Paired vs TO Phasing — Key Technical Details

Paired: HP:Z:"1-1" 字串格式（germline only）
Paired: HP:Z:"1-1" string format (germline only)
TO: HP:i:11 整數格式（germline + somatic）
TO: HP:i:11 integer format (germline + somatic)
關鍵：TO anchor 含 somatic → 循環依賴
Key: TO anchor includes somatic → circular dependency
QS: Paired AUC=0.754 vs TO AUC=0.497
QS: Paired AUC=0.754 vs TO AUC=0.497

Paired Mode

Tumor-Only Mode

Normal + Tumor BAM
LongPhase joint phasing

vs

Tumor BAM only
LongPhase-TO self phasing

HP tag: 字串格式
HP:Z:"1-1" / HP:Z:"2-1"

vs

HP tag: 整數格式
HP:i:11 / HP:i:21

Phasing anchor:
germline het SNPs only

vs

Phasing anchor:
germline + somatic variants

FP rate ~1%
Calling residual errors

vs

FP rate ~31%
Germline misclassified as somatic

QS AUC = 0.754
LOH penalty works correctly

vs

QS AUC = 0.497
LOH penalty direction reversed

4 / 24

> **Notes:** 快速過一下技術框架。ISM 針對每個 somatic variant，擷取周邊正負五千個 base pair 範圍的 reads，依照 haplotype 分群後比較甲基化模式。

這裡有一個重要的細節：reads 的 haplotype 歸屬不是 ISM 自己判定的，而是來自 LongPhase 做 haplotagging 之後 BAM 檔案裡的 HP tag。但 Paired 和 TO 的 HP tag 格式完全不同。Paired 用的是 LongPhase-s，HP tag 是字串格式，像是 HP:Z:1-1；TO 用的是 LongPhase-TO，HP tag 是整數格式，像是 HP:i:11。

另一個關鍵差異是：TO 的 phasing anchor 不只用 germline variants，還包含 somatic variant 本身——這造成一個潛在的循環依賴，我們待會會深入討論。

最後，我們先前已經嘗試過用 Quality Score 加分表做整合評分，在 Paired 模式拿到 AUC 0.754。

[過渡] 有了這個背景，接下來就是偵探故事的開始——我們拿到了第一批 TO 結果，但發現事情不太對。

---

## Slide 5: Act I: 發現異常

Act I: 發現異常

Act I: Anomaly Discovery

HP Tag Bug 修正與 LOH 方向翻轉的發現

5 / 24

> **Notes:** [短暫停頓] 接下來進入第一幕——發現異常。我們在擴展 TO 分析的過程中，先後碰到了兩個意外：一個是程式 bug，另一個是數據本身的系統性異常。

---

## Slide 6: HP Tag Bug 發現與修正

HP Tag Bug 發現與修正

HP Tag Bug Discovery and Fix

修正前（Bug）

• ReadParser 只處理 HP:Z（字串）

• HP:i:11/21/33 被靜默忽略

• TO Tier A/A+ ≈ 50%

• 大量 reads 無 haplotype 歸屬

→

InterSubMod
ReadParser.cpp 修正

修正後

• HP:i:11→"1-1", HP:i:21→"2-1"

• 整數格式完整支援

• TO Tier A/A+ ≈ 88%（≈ Paired 90%）

• 7 樣本全量重跑，舊數據作廢

教訓：程式不報錯不代表結果正確 — 是看圖片時才發現異常

6 / 24

> **Notes:** 接下來是第一個轉折。我們讓 AI 自動執行整個 ISM pipeline，程式跑完了，沒有任何錯誤訊息。但當我們打開結果圖片來看的時候，發現 TO 的 HP 數據明顯不對——大量 reads 缺少 haplotype 歸屬。

追查之後發現，ISM 的 ReadParser 只處理字串格式的 HP tag，也就是 Paired 用的 HP:Z。TO 的整數格式 HP:i 被靜默忽略了，完全沒有報錯。

修正之後效果很明顯：TO 的 Tier A 加 A+ 比例從大約百分之五十升到百分之八十八，接近 Paired 的百分之九十。修正前的 TO 數據全部作廢，七個樣本全量重跑。

這個 bug 的教訓是——程式不報錯不代表結果是對的。是看圖片的時候才發現的。

[過渡] 修正之後我們拿到了乾淨的數據，但真正的異常才剛剛開始浮現。

---

## Slide 7: LOH 定義與 LOH Enrichment 計算

LOH 定義與 LOH Enrichment 計算

LOH Definition and LOH Enrichment Calculation

LOH (Loss of Heterozygosity)

在某個基因組區域，原本應有的兩條 haplotype 中，一條的 reads 數量極端偏低或完全缺失。
In a genomic region, reads from one haplotype are extremely depleted or absent.

HP_Ratio = min(HP1, HP2) / (HP1 + HP2)

LOH-like 判定

當 HP_Ratio < 0.1 或 HP_Ratio > 0.9 時，判定為 LOH-like。
Classified as LOH-like when HP_Ratio < 0.1 or > 0.9.

LOH-like ⟺ HP_Ratio < 0.1 or HP_Ratio > 0.9

LOH Enrichment

FP 群中 LOH-like 的比例 除以 TP 群中 LOH-like 的比例。
Ratio of LOH-like proportion in FP group to TP group.

Enrichment = (FP LOH-like %) / (TP LOH-like %)

判定標準 Classification Criteria

• HP_Ratio 0.1/0.9 閾值為操作定義（需 sensitivity analysis）

• Enrichment > 1: FP 更容易 LOH → penalty 有效

• Enrichment < 1: TP 更容易 LOH → penalty 反向！

右圖：HP Ratio 散點圖 — X=Paired, Y=TO | 橘=TO-only LOH（水平帶 = self-phasing）| 綠=both LOH | 藍=neither

7 / 24

> **Notes:** 在看 LOH enrichment 結果之前，先定義幾個關鍵概念。

第一，LOH——Loss of Heterozygosity。在一個基因組區域裡，理論上兩條 haplotype 應該都有 reads 覆蓋。如果其中一條的 reads 極端偏低甚至完全沒有，就是 LOH。我們用 HP_Ratio 來量化：取 HP1 和 HP2 讀取數中較小的除以總和。

第二，LOH-like 判定。當 HP_Ratio 低於 0.1 或高於 0.9 時，我們判定這個區域為 LOH-like。注意這是一個操作定義，0.1 和 0.9 的閾值需要做 sensitivity analysis。

第三，LOH Enrichment。FP 群中 LOH-like 的比例，除以 TP 群中 LOH-like 的比例。如果大於 1，代表 FP 更容易表現出 LOH 特徵——這是我們期望的方向，因為代表 LOH penalty 可以有效懲罰 FP。如果小於 1，代表 TP 反而更容易出現 LOH——這就是問題了。

右邊這張散點圖是整個研究的核心視覺化之一，稍後會詳細討論。

[過渡] 有了這些定義，讓我們看看實際的 LOH enrichment 數據。

---

## Slide 8: LOH Enrichment 方向翻轉 — 7/7 樣本一致

LOH Enrichment 方向翻轉 — 7/7 樣本一致

LOH Enrichment Direction Reversal — Consistent Across 7/7 Samples

Paired: FP enriched (1.02-3.18x) | TO: TP enriched (0.852-0.956x) — 方向完全相反，QS LOH penalty 在 TO 下懲罰 TP

8 / 24

> **Notes:** 這一頁是整個調查的起點。當我們用修正後的數據計算 LOH enrichment——也就是 FP 的 LOH 比例除以 TP 的 LOH 比例——我們發現了一個令人意外的結果。

在 Paired 模式下，enrichment 是 1.194 倍，意思是 FP 比 TP 更容易被判為 LOH-like。這符合直覺，因為 FP 的 haplotype 支持度低，隨機波動容易產生假的 LOH 信號。

但在 TO 模式下，enrichment 是 0.805 倍——方向完全翻轉了。這代表 TP 反而比 FP 更容易被判為 LOH-like。而且這不是某個樣本的特例，七個樣本全部一致。

這直接導致一個嚴重的後果：我們的 Quality Score 裡有一個 LOH penalty，假設「LOH 代表可疑」。在 TO 下，這個懲罰反而在懲罰 TP。

[過渡] 方向翻轉聽起來嚇人，但到底有多嚴重？讓我們用同位點配對分析來量化。

---

## Slide 9: Act II: 深入追查

Act II: 深入追查

Act II: Deep Investigation

以 Paired LOH 與 TO LOH 交集觀察 × 288K 同位點配對 × 系統性假說排除

9 / 24

> **Notes:** [短暫停頓] 進入第二幕——深入追查。我們動用了二十八萬八千個同位點配對、九項系統性觀察、以及多個假說排除實驗，來釐清 TO LOH 偏差的規模和根因。

---

## Slide 10: 同位點配對分析規模

同位點配對分析規模

Matched Loci Analysis Scale

288,609

Paired-TO matched loci

85.5%

一致率

39,724

TO-only LOH

16-52x

Per-sample TO excess

10 / 24

> **Notes:** 我們有二十八萬八千多個同時存在於 Paired 和 TO 的 TP 位點——這是本週分析的核心數據集。用這些同位點配對來看 LOH 判定的一致性。

一致率百分之八十五點五，聽起來還不錯。但重點在不一致的部分——三萬九千七百多個位點，TO 單獨判為 LOH，而 Paired 認為不是。每個樣本的 TO excess ratio 在 16 到 52 倍之間。

[過渡] 讓我們看看完整的 concordance 分析圖。

---

## Slide 11: 五面板 Concordance 綜合分析

五面板 Concordance 綜合分析

Five-Panel Concordance Comprehensive Analysis

(A) 2x2 concordance matrix | (B) per-sample 四分類比例 | (C) TO excess ratio | (D) LOH rate 翻轉 | (E) 7 樣本 Enrichment 全部一致

11 / 24

> **Notes:** 這張五面板圖是 concordance 分析的全景。子圖 A 是二乘二的 concordance matrix——百分之五十五雙方都不判 LOH，百分之三十點六雙方都判 LOH，這是共識部分。

重點在不一致：百分之十三點八是 TO-only LOH，只有百分之零點六是 Paired-only。看子圖 C，每個樣本的 TO excess ratio 在 16 到 52 倍之間。這不是小數目——三萬九千多個位點在 Paired 看起來完全正常，但 TO 把它們判成了 LOH。

子圖 E 確認七個樣本全部一致，沒有例外。

[過渡] 那現在問題來了：這三萬九千多個 TO-only LOH 是怎麼回事？先看系統性觀察的全貌。

---

## Slide 12: 9 項系統性觀察 → TO 全面 Negative

9 項系統性觀察 → TO 全面 Negative

9 Systematic Observations → TO Entirely Negative

Paired vs TO QS Component Waterfall 對照

9 項觀察 + 4 項假說驗證 → 82 張圖表
9 observations + 4 hypothesis tests → 82 figures
TO 全部 AUC < 0.58，5/9 方向反轉
All TO AUC < 0.58, 5/9 direction reversed
LOH penalty 觸發：TP 44.5% vs FP 35.8%
LOH penalty trigger: TP 44.5% vs FP 35.8%
GQ: 0.811→(無) | AF: 0.665→0.418
GQ: 0.811→(none) | AF: 0.665→0.418

Paired QS Waterfall (AUC = 0.754)

TO QS Waterfall (AUC = 0.497)

左圖 Paired 各組件方向正確，右圖 TO 幾乎全部反轉 — LOH penalty 是最大反效果來源

12 / 24

> **Notes:** 問題不只是 LOH。我們做了九項系統性觀察加上四項假說驗證，總共產出八十二張圖表。結論是——TO 模式下所有特徵的 AUC 都低於 0.58，接近隨機。

看這張表：Paired 的 GQ 有 0.811，是最強的特徵，但 TO 根本沒有 GQ。AF 在 Paired 是 0.665，到了 TO 變成 0.418，方向翻轉了。最嚴重的是 QS 綜合評分：Paired 0.754，TO 只有 0.497，比丟銅板還不如。

看右邊的 waterfall 圖，TO 的 QS 各組件貢獻幾乎都是反向的。特別是 LOH penalty 的觸發率——TO 的 TP 觸發百分之四十四點五，FP 只有百分之三十五點八，等於在懲罰好的位點。

[過渡] 我們開始懷疑 TO 有系統性偏差。但第一個要排除的假說是：會不會只是 reads 太少造成的？

---

## Slide 13: QS Waterfall — TO AUC = 0.497（隨機）

QS Waterfall — TO AUC = 0.497（隨機）

QS Waterfall — TO AUC = 0.497 (Random)

Paired QS AUC = 0.754（有效）
Paired QS AUC = 0.754 (effective)
TO QS AUC = 0.497（等於丟銅板）
TO QS AUC = 0.497 (coin flip)
LOH penalty 方向完全反轉
LOH penalty direction completely reversed
QS 整合各組件全部失效
All QS components failed in TO mode

13 / 24

> **Notes:** 這張圖更直接地展示 QS 失效的原因。LOH penalty 的觸發率在 TO 模式下，TP 比 FP 更高——百分之四十四點五對百分之三十五點八。這意味著 Quality Score 裡的 LOH 懲罰機制，在 TO 下反而在懲罰 True Positive。

結合上一頁的 waterfall 圖，TO 的 QS 各組件的貢獻方向幾乎全部反轉。這不是某個組件的問題，是整體框架在 TO 模式下完全不適用。

[過渡] 那 TO-only LOH 是不是因為 reads 太少造成的？讓我們用數字來排除這個假說。

---

## Slide 14: 排除 Read Depth 假說 — TO-only LOH 深度反而更高

排除 Read Depth 假說 — TO-only LOH 深度反而更高

Ruling Out the Read Depth Hypothesis — TO-only LOH Has Higher Depth

TO-only LOH eff_hp_reads 中位數

68

高於 both_LOH 的 60，Cohen's d = +0.29

同位點 Paired HP_Ratio 中位數

0.509

完全平衡（TO = 0.026 極端 LOH）

HP_Ratio Cohen's d

-1.20

巨大效應量 | 86.5% 在 Paired 下完全平衡

Reads 數量夠，問題在 reads 被分配到哪邊，不是數量不夠

14 / 24

> **Notes:** 第一個要排除的假說：TO-only LOH 是不是因為 reads 太少，HP ratio 因為隨機波動跑到極端值？

答案是否定的。TO-only LOH 的 effective HP reads 中位數是 68，反而高於雙方都判 LOH 的 60。看右邊的箱型圖，橘色（TO-only）和綠色（both LOH）的讀取深度相當，甚至更高。

決定性的證據是同位點比較：同一個基因組位置，Paired 看到的 HP ratio 中位數是 0.509——幾乎完全平衡；但 TO 看到的是 0.026——極端的 LOH。Cohen's d 是負 1.20，這是巨大的效應量。

更具體地說，百分之八十六點五的 TO-only LOH 位點，在 Paired 模式下 HP ratio 都在 0.2 到 0.8 之間，完全平衡。Reads 的數量夠，問題在於 reads 被分配到哪邊。

[過渡] 如果不是 depth 的問題，那讓我們直接看 HP ratio 的分佈——一張散點圖就能看出答案。

---

## Slide 15: 排除低 Read Depth 影響

排除低 Read Depth 影響

14 / 24

> **Notes:** 第一個要排除的假說：TO-only LOH 是不是因為 reads 太少，HP ratio 因為隨機波動跑到極端值？

答案是否定的。TO-only LOH 的 effective HP reads 中位數是 68，反而高於雙方都判 LOH 的 60。看右邊的箱型圖，橘色（TO-only）和綠色（both LOH）的讀取深度相當，甚至更高。

決定性的證據是同位點比較：同一個基因組位置，Paired 看到的 HP ratio 中位數是 0.509——幾乎完全平衡；但 TO 看到的是 0.026——極端的 LOH。Cohen's d 是負 1.20，這是巨大的效應量。

更具體地說，百分之八十六點五的 TO-only LOH 位點，在 Paired 模式下 HP ratio 都在 0.2 到 0.8 之間，完全平衡。Reads 的數量夠，問題在於 reads 被分配到哪邊。

[過渡] 如果不是 depth 的問題，那讓我們直接看 HP ratio 的分佈——一張散點圖就能看出答案。

---

## Slide 16: HP Ratio 散點圖 — Self-Phasing 最直觀證據

HP Ratio 散點圖 — Self-Phasing 最直觀證據

HP Ratio Scatter — Most Direct Evidence of Self-Phasing

X=Paired HP_Ratio, Y=TO HP_Ratio | 綠=both LOH（四角）| 藍=neither（中央）| 橘=TO-only LOH（水平帶 = self-phasing 證據）

15 / 24

> **Notes:** 這張圖是整個調查裡最有說服力的一張。X 軸是 Paired 的 HP ratio，Y 軸是 TO 的 HP ratio。每個點是一個同位點配對。

先看顏色。綠色是雙方都判 LOH，集中在四個角落——雙方一致認為不平衡。藍色是雙方都不判 LOH，聚集在中央——雙方一致認為平衡。

現在看橘色——TO-only LOH。它們形成兩條水平帶。X 軸在 0.3 到 0.7，表示 Paired 看到的是平衡的；但 Y 軸在小於 0.1 或大於 0.9，表示 TO 把它們推到了極端值。同一批 reads，同一個基因組位置，因為 phasing 方式不同，haplotype 的分配完全不一樣。

[過渡] 散點圖是視覺化證據，接下來從 copy number 角度做補充。

---

## Slide 17: Read Depth Ratio + Copy Number 組成

Read Depth Ratio + Copy Number 組成

Read Depth Ratio + Copy Number Composition

Per-Sample LOH Read Depth Ratio（0.73x，非理論 0.5x）

CN LOH Composition（灰=neutral 60%, 紅=loss, 藍=gain）

60% copy-neutral → UPD 或 phasing artifact
60% copy-neutral → UPD or phasing artifact
TO copy-gain LOH = Paired 的 2 倍（15.8% vs 8.4%）
TO copy-gain LOH = 2x of Paired (15.8% vs 8.4%)

16 / 24

> **Notes:** 這一頁從 copy number 的角度做補充。如果 LOH 是真正的一個 allele 被刪除——也就是 copy-loss LOH——read depth 應該只有正常的一半，也就是 0.5 倍。但實際上，LOH 區域的 depth 是 non-LOH 的 0.73 倍，遠高於 0.5。

右邊的圖更清楚：大約百分之六十的 LOH 是 copy-neutral，也就是讀取深度正常但 haplotype 不平衡。這可能是 uniparental disomy，也可能是 phasing artifact。只有百分之二十四到三十是真正的 copy-loss。

值得注意的是，TO 的 copy-gain LOH 比例是 Paired 的兩倍——百分之十五點八對百分之八點四——表示 TO 在 copy-gain 區域有額外的系統性偏差。

[過渡] 所有的線索都指向同一個方向。進入第三幕——根因與修正。

---

## Slide 18: Act III: 根因與修正

Act III: 根因與修正

Act III: Root Cause and Correction

Self-Phasing 因果鏈確認 × 假說否決 × QS 修正

17 / 24

> **Notes:** [短暫停頓] 進入最後一幕——根因與修正。我們確認了 self-phasing circular dependency 的因果鏈，排除了所有甲基化替代方案，並完成了第一步 QS 修正。

---

## Slide 19: Self-Phasing 四步因果鏈 — 根因確認

Self-Phasing 四步因果鏈 — 根因確認

Self-Phasing Four-Step Causal Chain — Root Cause Confirmed

TO joint phasing
germline + somatic
共建 phasing graph

→

Self-phasing
somatic ALT reads
偏向一個 haplotype

→

Low-read 區域
被串入長 block
品質差

→

TP 更容易
判 LOH-like
44.5% vs 35.8%

量化：39,724 位點中 71.6% TO min(HP)=0（完全單側）| 同位點 Paired HP1:HP2 = 8:8 | Shapley: TP rate 貢獻 105.6%

18 / 24

> **Notes:** 現在把機制完整串起來。這是一個四步的邏輯鏈。

第一步：LongPhase-TO 在做 phasing 的時候，把 germline 和 somatic variant 一起丟進 phasing graph。我們從 LongPhase 的原始碼 PhasingGraph.cpp 確認了這件事。

第二步：因為 TP 的 somatic variant 有真實的 ALT allele reads，這些 reads 天然偏向一個 haplotype。也就是說，variant 自己的 reads 在影響自己被歸到哪個 haplotype——這是一個循環依賴。

第三步：在 read 數量少的區域，TO 仍然嘗試把它串接到 phasing block 裡面，但品質遠不如 Paired。

第四步：結果就是 TP 更容易被判成 LOH-like——觸發率百分之四十四點五，比 FP 的百分之三十五點八還高。

定量上：三萬九千多個 TO-only LOH 位點中，百分之七十一點六在 TO 下某一側 HP 的 read 數是零——完全單側。但同一個位點，Paired 看到的 HP1 比 HP2 是 8 比 8，完美平衡。

[過渡] 在確認機制之後，我們也嘗試過用甲基化特徵來解決——但全部被否決了。

---

## Slide 20: 量化驗證 — 觸發率與 HP Ratio 分佈

量化驗證 — 觸發率與 HP Ratio 分佈

Quantitative Verification — Trigger Rate and HP Ratio Distribution

Self-Phasing Circular Dependency 的定量證據

1

TO TP LOH 觸發率 44.5% vs FP 35.8%（反向）

2

71.6% TO-only LOH 位點 min(HP) = 0

3

同位點 Paired HP1:HP2 = 8:8（完全平衡）

4

Shapley 分解：TP rate 貢獻 105.6%

Self-phasing 不是假說，是已確認的機制 — 62% LOH 消失、d=-1.20、31.2% self-phasing LOH

19 / 24

> **Notes:** 這一頁提供定量驗證。左圖是 LOH penalty 的觸發率——TO 的 TP 觸發百分之四十四點五，FP 只有百分之三十五點八，方向完全反轉。右圖是 HP ratio 的分佈，可以看到 TO 模式下 reads 被極端地推向一側。

四個關鍵數字：第一，TO TP 的 LOH 觸發率高於 FP——反向。第二，百分之七十一點六的 TO-only LOH 位點，某一側的 HP read 數是零——完全單側。第三，同一個位點在 Paired 看到的是 8 比 8 的完美平衡。第四，Shapley 分解顯示 TP rate 貢獻了百分之一百零五點六。

底線：self-phasing 不是假說，是已確認的機制。

[過渡] 在確認根因之前，我們也嘗試過用甲基化特徵來解決——但全部被否決了。

---

## Slide 21: 四項替代方案全部否決

四項替代方案全部否決

Four Alternative Hypotheses All Rejected

問題在 Phasing 不在 Methylation

O11: Within-group Heterogeneity

AUC 0.845 AUC 0.530

n_reads confound

✗

O12: LOH 甲基化場景分類

AlleleDelta 顯著 AUC < 0.58

AF confound + L2 collider bias

✗

O13: 跨區域 Correlation

FP > TP 顯著 差異消失

shared read count confound

✗

N4: HP0 Filter

TP loss ≤ 2% FP removal = 0%

無法在安全約束下移除 FP

✗

教訓：L1→L2→L3 三層驗證 | 表面顯著的特徵可能只是 confound | 問題在 phasing 層，非 methylation

20 / 24

> **Notes:** 在確認機制之前，我們其實花了不少時間試圖用甲基化特徵來解決 TO 的問題。結果四項替代方案全部被否決。

O11 嘗試用 within-group methylation heterogeneity 來區分，表面上 epipolymorphism 的 AUC 有 0.845，看起來很強。但控制 read count 之後掉到 0.530——原來是 read 數量的 confound。

O12 嘗試用 LOH 區域的甲基化場景分類，結果 AlleleDelta 本身就等同於 AF，加上我們發現了 L2 collider bias。校正後全部 AUC 低於 0.58。

O13 嘗試跨區域的甲基化 correlation，看起來 FP 比 TP 的 correlation 高，但控制 shared read count 之後差異完全消失。

還有 N4 — HP0 filter，在容許 TP loss 百分之二以下的條件，FP removal 是零。

結論很明確：問題在 phasing 層面，不在 methylation。

[過渡] 既然原因確認了，我們已經做了第一步修正。

---

## Slide 22: QS Mode-Aware 修正 — 消除反效果

QS Mode-Aware 修正 — 消除反效果

QS Mode-Aware Correction — Eliminating Adverse Effects

修正前（TO 失效）

• LOH penalty = -25（反向懲罰 TP）

• Verify bonus = +15（依賴錯誤 LOH）

• TO QS AUC = 0.497（比隨機還差）

→

commit b9eaba7

修正後（TO 中性）

• LOH penalty = 0（停用）

• Verify bonus = 0（停用）

• TO QS AUC ≈ 0.546（中性）

意義：「不再傷害」而非「變好」| 0.546 仍接近隨機 | TO 需要全新策略，不是修修補補

21 / 24

> **Notes:** 基於以上所有證據，我們做了第一步修正：在 TO 模式下停用 LOH penalty 和 verify bonus。

原本的做法是用一張 Quality Score 加分表，整合多個組件的懲罰和加分。在 Paired 模式效果不錯，AUC 0.754。但在 TO，LOH penalty 的方向完全反了，所以 TO QS 降到 0.497——比隨機還差。

停用之後，TO QS 理論上回到 0.546 左右。但我要強調，這個修正的意義是「不再傷害」，而不是「變好」。0.546 還是接近隨機。TO 需要的是全新的策略，不是修修補補。

[過渡] 那接下來怎麼辦？讓我分享修正方向的路線圖。

---

## Slide 23: 修正方向路線圖 — 四項行動

修正方向路線圖 — 四項行動

Correction Roadmap — Four Action Items

1

Phasing Anchor 策略

排除或降權 somatic variant 在 LongPhase-TO phasing graph 中的角色

2

Phasing Block 品質過濾

TO 在低 read 區域串出的長 block 品質差 → 加 block quality threshold

3

HP Tagging 後處理

根據 phasing block 大小/read support 動態調整 LOH 判定門檻

4

Paired Phasing 作為 Benchmark

288,609 matched loci 數據做 TO phasing quality calibration

已有的改善基礎

288,609 matched loci → TO phasing quality calibration 數據
86.5% TO-only LOH 在 Paired 下平衡 → 可校正的 false LOH
Copy-gain 區域 2x 偏差 → 需要 depth-aware 校正

22 / 24

> **Notes:** 根因在 LongPhase-TO 的 phasing 策略，所以修正方向也要從 phasing 著手。我們列出四個分析方向。

第一，分析能不能把 somatic variant 從 phasing anchor 中排除或降權。第二，TO 在低 read 區域串出的長 phasing block 品質很差，能不能加一個品質門檻。第三，在 ISM 這一層，能不能根據 phasing block 的資訊動態調整 LOH 判定的門檻。第四，利用二十八萬八千個同位點數據，用 Paired phasing 作為 benchmark 來校正 TO。

好消息是，我們已經有很好的基礎。百分之八十六點五的 TO-only LOH 在 Paired 下是平衡的，這代表這些是理論上可以被校正的 false LOH。Copy-gain 區域有額外的偏差，需要結合 depth 資訊做校正。

[過渡] 最後讓我把本週的結論做一個完整的整理。

---

## Slide 24: 核心結論 — Paired/TO 是根本不同的問題空間

核心結論 — Paired/TO 是根本不同的問題空間

Core Conclusions — Paired/TO Are Fundamentally Different Problem Spaces

維度 | Paired | TO
--- | --- | ---
FP rate | 1% | 31%
LOH enrichment | 1.194x (FP) | 0.805x (TP)
特徵方向反轉 | — | 5/9 反轉
QS AUC | 0.754 | 0.497 (隨機)
TO-only LOH excess | — | 16-52x
HP_Ratio Cohen's d | — | -1.20

1 Paired/TO 分離建模是必要條件

2 TO LOH 偏差源自 self-phasing circular dependency

3 60% LOH 是 copy-neutral → 部分是 phasing artifact

4 甲基化特徵無法解決（O11-O13 全否決）

5 下一步：TO phasing 程式碼分析與修正方案

Limitations

所有結論基於 7 個高 purity 細胞系（purity ≈ 100%），臨床低 purity 需額外驗證。LOH-like 閾值 (0.1/0.9) 為操作定義，需 sensitivity analysis。

23 / 24

> **Notes:** 最後做一個完整的整理。這張表列出 Paired 和 TO 在六個維度的對比——每一項都指向同一個結論：這兩個模式是根本不同的問題空間。

五個層次的結論。第一，Paired 和 TO 分離建模不是選項，是必要條件。第二，TO LOH 偏差的根因是 self-phasing circular dependency——LongPhase-TO 讓 somatic variant 參與自己的 phasing，而不是 read depth 或 copy number 的問題。第三，百分之六十的 LOH 是 copy-neutral，其中一部分可能是 phasing artifact，這對 Paired 模式也值得關注。

第四，甲基化特徵無法解決 TO 的問題——O11、O12、O13 三項假說全部被否決。問題在 phasing 層面。第五，下一步是分析 TO 的 phasing 程式碼，設計修正方案。

底下列了 limitations：所有結論基於七個高 purity 細胞系，臨床低 purity 需要額外驗證；LOH 的閾值 0.1 和 0.9 是操作定義，需要做 sensitivity analysis。

[過渡] 最後一頁，下週的行動計畫和想跟大家討論的問題。

---

## Slide 25: 討論與下週方向

討論與下週方向

Discussion and Next Week's Direction

下週行動計畫

優先級 | 行動 | 依據
--- | --- | ---
P0 | TO phasing block 品質分析 | 根因在 phasing
P0 | Paired/TO 分離策略框架 | 方向完全相反
P1 | Paired ML 特徵集 (GQ+AF+LOH) | GQ AUC=0.811
P2 | TO phasing anchor 修正方案 | 排除 somatic variant

開放討論問題

Q1

LongPhase-TO 排除 somatic variant 後，phasing 品質是否足夠？

Q2

TO 獨立 QS 框架應以什麼為 anchor？

Q3

Paired phasing 可以做 TO 的 calibration reference 嗎？

謝謝！歡迎提問。

24 / 24

> **Notes:** 下週的行動按優先級排。P0 有兩件：第一是分析 TO 的 phasing block 品質，比較 block 長度和 haplotag 一致性，因為根因在 phasing。第二是開始設計 Paired 和 TO 的分離策略框架，因為方向完全相反，不能用同一套規則。

P1 是 Paired 的 ML 特徵集，GQ 的 AUC 有 0.811 是最強特徵，加上 AF 和 LOH subtype 應該能建出不錯的 Paired 模型。P2 是評估能否從 LongPhase-TO 排除 somatic variant 做 phasing anchor，這需要深入分析原始碼。

最後三個問題想跟大家討論。第一，排除 somatic variant 之後 TO 的 phasing 品質還夠不夠？第二，TO 獨立的評分框架應該以什麼為基準？第三，Paired 的 phasing 結果能不能作為 TO 的校正參考？

本週最重要的結論：Paired 和 TO 是根本不同的問題空間。Self-phasing circular dependency 是已確認的機制，不是假說。下週我們會從根源——phasing 演算法層——著手修正。以上是本週報告，謝謝，歡迎提問。

---
