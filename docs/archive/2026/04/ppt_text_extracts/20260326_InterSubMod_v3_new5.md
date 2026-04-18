# 20260326_InterSubMod_v3_new5.pptx

- Slides: 31
- Date: 2026-03-26

## Slide 1: InterSubMod 整體概念、現有觀察
、未來方向 報告

InterSubMod 整體概念、現有觀察
、未來方向 報告

InterSubMod: Concept Overview, Current Observations, and Future Directions

報告人

廖子游

日期

2026-03-26

單位

中正大學資工系 黃耀廷 教授 Lab405 實驗室

封面與簡介

技術方法

管線與品質

實驗結果

文獻對照

未來方向

【請插入圖片】

Insert your figure here

> **Notes:** S01 封面。右側為空白區域，請自行插入示意圖。報告人：廖子游，單位：中正大學資工系Lab405。

---

## Slide 2: S02 ｜ 封面與簡介

S02 ｜ 封面與簡介

目錄

Table of Contents

1

封面與簡介
Cover & Introduction

ISM 定義

研究動機

五個目標

2

技術方法
Technical Methods

Read×CpG 矩陣

距離度量

階層聚類

3

管線與品質
Pipeline & Quality

LOH 偵測

品質分層

ISM 管線

4

實驗結果
Experimental Results

AUROC+LOH

FP 分析

已確認發現

5

文獻對照
Literature Review

Long-Read POG

CAMDAC

EVOFLUx

6

未來方向
Future Directions

四階段策略

研究路線圖

總結

> **Notes:** S02 目錄。六章節各以品牌色區分，底部關鍵詞條方便快速導覽。

---

## Slide 3: S03 ｜ 封面與簡介

S03 ｜ 封面與簡介

為何需要 ISM？

Why InterSubMod? — The Problem Cascade

如何從 Long-read 數據中辨別腫瘤亞克隆結構？

How to resolve tumor subclonal structure from long-read sequencing?

TP/FP 辨識困難

TP/FP Discrimination

短讀長測序工具將
Germline ASM 誤判為
體細胞甲基化變化

Short-read tools misidentify
Germline ASM as somatic
methylation changes

Germline ASM

VerificationClass

等位基因特異性甲基化

Allele-Specific Methylation

HP1/HP2 讀序甲基化
差異是否源於體細胞
突變而非正常 ASM？

Is HP1/HP2 methylation diff
from somatic mutation or
normal ASM?

ASM

Second Hit

亞克隆異質性

Subclonal Heterogeneity

不同亞克隆的甲基化
特徵如何在 reads 中
被精確識別與分離？

How are distinct subclone
methylation signatures
identified in reads?

subclone

LOH

ISM：整合 Read × CpG 矩陣 + 距離聚類 + 統計顯著性，系統性解決上述三大問題
ISM: Read×CpG matrix + distance clustering + significance testing → systematic solution

> **Notes:** S03 問題層級展開。中央大問題框展開三個子問題：TP/FP、ASM、亞克隆。底部 ISM 答案框。

---

## Slide 4: S04 ｜ 封面與簡介

S04 ｜ 封面與簡介

ISM 五個研究目標

ISM Five Research Objectives — A Comprehensive Framework

1

🔬

CpG 關聯分析

Methylation-SNV Association

分析單個 SNV 視窗內每個 CpG
位點與 ALT/REF/HP 標籤的關聯

輸入

BAM + VCF

輸出

per-CpG label + strength

2

🧬

癌症克隆狀態

Cancer Clone Architecture

找出腫瘤中的亞克隆位置、
排列與克隆階層結構

輸入

multi-SNV + methylation

輸出

subclone architecture map

3

⏱

Second Hit 時序推斷

Second Hit & Event Ordering

推斷甲基化演化事件的先後
順序（LOH/CNV + SNV/INDEL）

輸入

HP + LOH + CNV + SNV

輸出

methylation timeline

4

🏥

TO / Purity 場景

Tumor-Only Scenarios

無 normal sample 下，從甲基
化模式推斷正常細胞背景

輸入

TO BAM only

輸出

purity/normal estimate

5

📈

整合提升 F1

Somatic F1 Boost

整合甲基、HP、LOH 等特徵
提升 TP retention + FP reduction

輸入

all above features

輸出

F1 score improvement

> **Notes:** S04 五個目標。上排2塊：CpG關聯分析、癌症克隆狀態。下排3塊：Second Hit、TO/Purity、整合F1。

---

## Slide 5: S05 ｜ 封面與簡介

S05 ｜ 封面與簡介

ONT 長讀序三合一優勢

ONT Long-Read Three-in-One Advantages

ONT Long Read (~10–50 kb)

① 鹼基序列

Base Sequence
(SNV/INDEL calling)

② MM/ML 甲基化

Methylation
(CpG probability)

③ HP 單倍型

Haplotype
(HP1 / HP2)

短讀長 (Illumina)

長讀長 (ONT)

鹼基序列

✓ 高精度

✓ 含 SNV

甲基化

✗ 需額外 BSseq

✓ MM/ML 標籤

單倍型

△ 需額外定相

✓ HP1/HP2 標籤

讀序長度

150 bp

~15 kb 平均

MM/ML：ModKit 甲基化機率標籤格式

HP tag：LongPhase 單倍型分相標籤

> **Notes:** S05 ONT 三合一。一條長讀序覆蓋三種資訊，底部比較表說明短讀長的侷限。

---

## Slide 6: S06 ｜ 技術方法

S06 ｜ 技術方法

核心假設鏈

Core Hypothesis Chain — ISM Theoretical Foundation

核心假設：體細胞突變（SNV）會改變鄰近 CpG 位點的甲基化模式

Core: somatic SNVs alter methylation patterns at neighboring CpG sites

假設 A：TP/FP 辨識

TP/FP Discrimination

ALT reads 的甲基化
模式與 REF reads 顯著
不同 → TP 標誌

ALT read methylation differs
significantly from REF reads
→ somatic TP marker

✓ AUROC > 0.75
✓ 甲基化方向一致

假設 B：ASM & Second Hit

ASM & Event Ordering

HP1-1（ALT+HP1）
讀序的甲基化異常
源於 somatic 事件

HP1-1 (ALT+HP1) read
methylation abnormality
from somatic events

→ Second Hit 推斷
→ LOH 甲基化保留

假設 C：亞克隆結構

Subclonal Architecture

不同甲基化模式的
Read 群集代表不同
腫瘤亞克隆

Read clusters with distinct
methylation represent
different tumor subclones

→ VerificationClass
→ Subclone 分離

補充假設 D：順式調控與特殊演化路徑（Cis-acting factors & Epigenetic Evolutionary Paths）

腫瘤克隆的甲基化格局可能受 SNV 附近的 cis 調控元件影響，形成獨立的演化路徑 → 可作為克隆溯源依據
Tumor clone methylation may be influenced by cis-regulatory elements near SNVs, forming independent evolutionary paths → clonal tracing basis

Second Hit：第二次打擊事件：LOH 後甲基化異常即為 Second Hit 的表觀遺傳證據

ASM：Allele-Specific Methylation — 等位基因特異甲基化

> **Notes:** S06 核心假設鏈。中央大框說明核心假設，三個分支分別說明三大應用場景。底部補充假設D（順式調控）。

---

## Slide 7: S07 ｜ 技術方法

S07 ｜ 技術方法

SNV 視窗與 Read 標籤

SNV Window & Read Labeling — Cross-section View

Reference Genome (chr17:...)

SNV

SNV Window (±2500 bp)

HP1-1

ALT

HP1-1

ALT

HP1

REF

HP2

REF

HP2

REF

HP2

REF

★ CpG位點3、4（橘色高亮）= HP1-1 特有高甲基位點（HP1 reads 此處低甲基）

HP1讀序(REF)

HP2讀序(REF)

ALT讀序(HP1-1)

HP1-1：HP1 背景 + ALT 等位基因 read（體細胞標記）

SNV Window：SNV 位點前後 ±2500 bp 分析範圍

> **Notes:** S07 SNV視窗。基因組橫條 + 更大CpG三角標記 + 6條讀序分HP1-1/HP1/HP2 + 甲基化dots。HP1-1術語首次出現。

---

## Slide 8: S08 ｜ 技術方法

S08 ｜ 技術方法

Read × CpG 甲基化矩陣

Read × CpG Methylation Matrix — Core ISM Data Structure

CpG1

HP1-meth

CpG2

HP2-meth

CpG3

HP1-meth

CpG4

HP1-meth

CpG5

HP2-meth

CpG6(SNV)

ALT-meth

HP1-1

R1

ALT

0.90

0.10

0.88

0.85

0.12

0.90

HP1-1

R2

ALT

0.88

0.12

0.85

0.82

0.10

0.87

HP1

R3

REF

0.85

0.15

0.80

0.85

0.15

0.10

HP2

R4

REF

0.10

0.88

0.12

0.10

0.85

0.08

HP2

R5

REF

0.12

0.85

0.10

0.12

0.82

0.10

HP2

R6

REF

0.08

0.90

0.15

0.08

0.88

0.12

> 0.80 高甲基化

0.50–0.80

0.20–0.50

< 0.20 未甲基化

HP1-1：HP1 背景
攜帶 ALT（體細胞克隆）

HP1：HP1 背景
REF 讀序（野生型）

HP2：HP2 背景
REF 讀序

甲基化概率：0.00=未甲基化, 1.00=完全甲基化（MM/ML標籤）

> **Notes:** S08 Read×CpG矩陣。值格式已修正為0.90等小數格式。HP1-1術語已在S07介紹。

---

## Slide 9: S09 ｜ 技術方法

S09 ｜ 技術方法

距離計算：六種度量

Distance Metrics — Six Measures for Read Comparison

L1

曼哈頓距離
Manhattan

Σ|xᵢ - yᵢ|

L2

歐氏距離
Euclidean

√Σ(xᵢ-yᵢ)²

NHD

正規化海明
Norm. Hamming

Σ|round(xᵢ)-round(yᵢ)|/n

Bernoulli

伯努利 KL
Bernoulli KL

Σ[xᵢlog(xᵢ/yᵢ)+...]

★ 推薦

Jaccard

Jaccard 距離
Jaccard

1 - |A∩B|/|A∪B|

Cosine

餘弦距離
Cosine

1 - (x·y)/(|x||y|)

Bernoulli 計算範例

HP1-1 (ALT)

0.90

CpG1

0.10

CpG2

0.88

CpG3

0.85

CpG4

HP2  (REF)

0.10

0.88

0.12

0.10

D_Bernoulli(HP1-1, HP2) = Σ KL(p₁‖p₂) ≈ 2.14

正規化至 [0,1]：d_norm = tanh(d/n_CpG) ≈ 0.85

為何選用 Bernoulli？

• 甲基化概率為 [0,1] 連續值，Bernoulli 天然匹配

• 非對稱敏感（高甲基 vs 低甲基差異更明顯）

• 比 L1/L2 更能區分 HP1-1 vs HP2 模式

Bernoulli KL：甲基化概率向量的 KL 散度，非對稱距離

> **Notes:** S09 六種距離度量。左側表格列出6種，右側 Bernoulli 計算範例。Bernoulli術語首次出現。

---

## Slide 10: S10 ｜ 技術方法

S10 ｜ 技術方法

階層聚類：讀序分群概念圖

Hierarchical Clustering — R1-R6 Reads, 6×6 Distance Matrix & UPGMA Tree

讀序甲基化模式

R1-R6 read methylation patterns

R1

HP1-1

R2

HP1-1

R3

HP1

R4

HP2

R5

HP2

R6

HP2

HP1
群

HP2
群

高

中

低

距離矩陣（標準化 0-1）

Normalized Bernoulli distance [0,1]

R1

R2

R3

R4

R5

R6

R1

0.00

0.05

0.22

0.89

0.91

0.87

R2

0.05

0.00

0.20

0.88

0.90

0.86

R3

0.22

0.20

0.00

0.86

0.88

0.84

R4

0.89

0.88

0.86

0.00

0.06

0.09

R5

0.91

0.90

0.88

0.06

0.00

0.07

R6

0.87

0.86

0.84

0.09

0.07

0.00

0.0 同一讀序

≤0.30 相似（同群）

≥0.75 相異（跨群）

UPGMA 聚類樹 (6 reads)

HP1群(R1+R2+R3) vs HP2群(R4+R5+R6)

R1

R2

R3

R4

R5

R6

0.05

0.21

0.06

0.08

0.87

HP1 群

HP2 群

切割

UPGMA：Unweighted Pair Group Method — 分層聚類算法

> **Notes:** S10 三面板：左=6條讀序(R1-R6)甲基化模式，中=標準化6×6距離矩陣，右=UPGMA樹(HP1群R1+R2+R3 / HP2群R4+R5+R6)。

---

## Slide 11: S11 ｜ 技術方法

S11 ｜ 技術方法

雙向顯著性框架

Bidirectional Significance Framework — Two Complementary Approaches

方法 A：依標籤驗證群內外顯著

Label-First Approach

① 依 HP/Allele 標籤分組讀序

② 比較群內 vs 群外甲基化分布

③ 卡方檢驗 + Cramér's V 效應量

結果示例

標籤

顯著比例

判定

Allele

100%

✓ 高顯著

HP

98%

✓ 高顯著

Strand

0%

✗ 不顯著

→ 驗證 HP/Allele 是有效甲基化標籤
   Validates HP/Allele as effective methylation label

方法 B：依分群來驗證標籤顯著

Cluster-First Approach

① 計算讀序距離矩陣（Bernoulli）

② UPGMA + Silhouette Score 分群

③ Fisher 檢驗標籤純度

特性分析

優點

依數據分群，不受先驗標籤影響；能偵測 subclone

缺點

受參數影響，可能切出離群值；計算量較大

適用

HP1-1 獨立成群時 → 偵測 Subclone 結構

→ Fisher 純度顯著 → 分群有意義
   Fisher purity significance confirms cluster validity

VerificationClass

Strong

Weak

Subclone

Noise

PERMANOVA：置換型多變量方差分析，用於讀序群落比較

Cramér's V：卡方效應量指標 [0,1]，≥0.3 為中度關聯

> **Notes:** S11 雙向顯著性。方法A(Label-First)含結果表，方法B(Cluster-First)含特性分析，匯入VC。

---

## Slide 12: S12 ｜ 技術方法

S12 ｜ 技術方法

VerificationClass：四類判斷

VerificationClass — Four Outcome Categories

判斷流程

兩方法結果

p值 + 效應量

方法A顯著?

Cramér's V≥0.3

方法B顯著?

Fisher p<0.05

Subclone?

Silhouette k=3

Strong: 雙方皆顯著

Weak: 僅一方顯著

Subclone: HP1-1獨立成群

Noise: 雙方皆不顯著

✓

Strong

Strong Signal

兩種方法均顯著
Cramér's V ≥ 0.3
Fisher p < 0.05

Both methods significant
Cramér's V ≥ 0.3 & Fisher p < 0.05

Precision ~88%

△

Weak

Weak Signal

單一方法顯著
效應量較小
需更多樣本確認

Single method significant
Small effect size; needs more samples

Precision ~64%

⑂

Subclone

Subclonal Event

HP1-1 獨立成群
Silhouette k=3 顯著
→ 可能有亞克隆

HP1-1 clusters independently
Silhouette k=3 significant → subclone

Subclone flag ✓

✗

Noise

Noise / Artifact

無顯著差異
Germline ASM 干擾
排除此 SNV

No significant difference
Germline ASM interference; exclude SNV

FP filter ✓

Strong：AUROC>0.8 + Cramér's V≥0.3，最高信心體細胞變異

Subclone：HP1-1 讀序聚類與 HP1 分離 → 亞克隆結構證據

> **Notes:** S12 VerificationClass。上方判斷流程（決策樹），下方2×2四象限卡片含閾值標準。

---

## Slide 13: S13 ｜ 管線與品質

S13 ｜ 管線與品質

LOH 偵測與品質分層

LOH Detection & Quality Score Stratification

LOH 示意（Read 層級）

Loss of Heterozygosity — read-level evidence

正常細胞（平衡）

HP1

HP1

HP1

HP2

HP2

HP2

HP1: 3 reads (50%)

HP2: 3 reads (50%)

腫瘤化（LOH 發生）

LOH 腫瘤細胞（HP2 缺失）

HP1

HP1

HP1

HP1

HP1

HP2?

HP1: 5 reads (83%)

HP2: 1 read (17%)

LOH Ratio = (HP1 + HP1-1) reads / (HP1 + HP1-1 + HP2) ≥ 0.8 → LOH 判定

※ HP1-1（HP1背景+ALT讀序）在LOH偵測中計入HP1群——LOH反映HP2缺失，而非ALT allele

甲基化保留效應：

HP2 缺失 → HP2 甲基化模式在剩餘讀序中消失
→ HP1 甲基化異常更突出（Second Hit 偵測窗口）

品質分層（Quality Score）

QualityScore — composite metric for read group reliability

1

讀序數量
Read Count

HP1/HP1-1/HP2 各自讀序數

≥ 5 reads 為門檻

2

ALT 比例
ALT Fraction

HP1-1 reads / HP1 total

SNV 體細胞比例

3

甲基化一致性
Meth Consensus

群內 CpG 甲基化方差

< 0.15 為高品質

4

統計顯著性
P-value

PERMANOVA / 卡方 p-value

< 0.05 顯著

QualityScore = f(ReadCount, ALTfrac, MethConsensus, Pvalue) → [0, 1]

LOH：Loss of Heterozygosity — 雜合性缺失，HP2 reads 比例降低

QualityScore：綜合讀序數量、ALT比例、甲基化一致性、統計顯著性

> **Notes:** S13 LOH偵測。正常3+3均衡，LOH腫瘤5+1，右側QualityScore四組成元件。LOH術語首次出現。

---

## Slide 14: S14 ｜ 管線與品質

S14 ｜ 管線與品質

ISM 完整分析管線

ISM Full Analysis Pipeline — Four-Layer Architecture

輸入層
Input Layer

BAM (HP-tagged)

VCF (phased SNV)

Reference genome

→ SNV window extraction

Read 解析層
Read Parsing Layer

HP/Allele 標籤讀取
HP/Allele tag reading

MM/ML 甲基化解碼
MM/ML methylation decode

CpG 位點定位
CpG site localization

→ Read×CpG matrix

分析層
Analysis Layer

距離矩陣計算（6種）
Distance matrix (6 metrics)

UPGMA 階層聚類
Hierarchical clustering

PERMANOVA + 卡方檢驗
PERMANOVA + Chi-sq test

→ VerificationClass

輸出層
Output Layer

VerificationClass 判斷
VC decision

QualityScore 分層
QS stratification

TSV + 視覺化報告
TSV + visualization

→ Downstream use

管線統計

7

樣本

<300ms

per region

6

距離度量

4

VC 類別

16

文件記錄

MM/ML：ONT 甲基化概率標籤，由 ModKit 解碼

HAP tag：LongPhase 分相後的 HP1/HP2 標籤

> **Notes:** S14 完整管線。四層架構：輸入→Read解析→分析→輸出，右側管線統計數字。

---

## Slide 15: S15 ｜ 管線與品質

S15 ｜ 管線與品質

方法學審查總覽

Methodology Review — Implementation Audit

1,200+

git commits

程式碼提交

< 300ms

per region

單 Region 處理時間

7

樣本

腫瘤/正常樣本對

16

份技術文件

方法論記錄文件

6

分析模組

核心 C++ 模組

開發時間軸

2025-12

ISM 概念
設計

2026-01

C++ 核心
實作

2026-02

7樣本
全量分析

2026-03

FP分析+
論文定位

2026-Q2

Phase 1-4
策略啟動

已確認方法與現況

✓

Bernoulli 距離

最佳讀序甲基化距離度量  (best methylation distance)

✓

PERMANOVA 框架

讀序群落顯著性分析  (read group significance test)

✓

LOH 偵測邏輯

HP ratio ≥ 0.8 門檻  (LOH threshold)

✓

AUROC 評估

HP1-1 vs REF 分類效能 > 0.75  (AUC > 0.75)

當前瓶頸

• Germline ASM 干擾（98.7% FP 根因）
  Germline ASM interference (98.7% FP root cause)

• TO 場景無 normal sample 背景
  Tumor-only lacks normal sample background

• Phase 2 多 SNV 合併分析尚未實作
  Phase 2 multi-SNV joint analysis not yet implemented

> **Notes:** S15 方法學審查。頂部5個大數字，時間軸，已確認方法，當前瓶頸。

---

## Slide 16: S16 ｜ 實驗結果

S16 ｜ 實驗結果

本週實驗概覽（關鍵頁）

This Week's Experimental Summary — Key Results

★ 最重要結果頁

任務 1：SNV 分類

廖子游

Task 1: SNV Classification

做了什麼 / What

7樣本全量 AUROC 評估
HP1-1 vs HP1/HP2 讀序
甲基化分類效能

得到什麼 / Result

AUROC 0.75–0.88
（依 LOH 程度分層）

結論 / Conclusion

LOH ≥ 0.8 區間效果最佳
甲基化特徵可有效分離HP1-1

✓ 達成

任務 2：FP 根因分析

廖子游

Task 2: FP Root Cause

做了什麼 / What

分析 98.6% FP 屬於
raw_absent 類型
追蹤 Germline ASM 干擾

得到什麼 / Result

HCC1395_HKU germline
ASM 影響 94% 位點
0.155 CpG 甲基化差異

結論 / Conclusion

FP 無法透過 ISM 過濾
FN rescue 才是正確方向

△ 部分達成

任務 5：F1 提升

秉庭

Task 5: F1 Score Boost

做了什麼 / What

Paired F1 ΔF1 計算
TO 甲基化特徵分析
HCC1395_HKU 場景驗證

得到什麼 / Result

ΔF1 = +0.0094
（TO 策略 HCC1395 樣本）
統計顯著 p < 0.05

結論 / Conclusion

TO 策略在特定樣本有效
需更多樣本驗證普遍性

✓ 達成

> **Notes:** S16 最重要實驗頁。三欄：任務1/任務2/任務5，各有做了什麼→得到什麼→結論。

---

## Slide 17: S17 ｜ 實驗結果

S17 ｜ 實驗結果

AUROC 與 LOH 數據詳覽

AUROC & LOH Analysis — Tumor Samples (HG002 excluded as germline reference)

AUROC（依 LOH 程度分層）

分析目的：量化 HP1-1 讀序甲基化特徵能否有效區分 ALT（HP1-1）vs REF（HP1/HP2）讀序

LOH ≥ 0.8

0.88

LOH 0.6–0.8

0.80

LOH 0.4–0.6

0.74

LOH < 0.4

0.65

全量平均

0.76

Random line
(AUC=0.75)

LOH 為何重要：
高 LOH → HP1 背景純化 → HP1-1 甲基化信號更清晰
低 LOH → HP1/HP2 混雜 → AUROC 下降
Why LOH matters: high LOH purifies HP1 background → clearer HP1-1 signal; low LOH → AUROC↓

腫瘤樣本 AUROC × VC 類別

腫瘤樣本

Strong

Weak

Subcl.

Noise

AUROC

HCC1395_HKU

42%

28%

8%

22%

0.88

HCC1395_5kHz

35%

31%

5%

29%

0.82

COLO829

38%

25%

12%

25%

0.85

H1437

29%

33%

6%

32%

0.74

H2009

31%

30%

7%

32%

0.76

AUROC：ROC曲線下面積，隨機分類=0.5，越高越好

> **Notes:** S17 AUROC+LOH。移除HG002（非腫瘤），只保留5腫瘤樣本，加入LOH重要性說明。AUROC術語首次出現。

---

## Slide 18: S18 ｜ 實驗結果

S18 ｜ 實驗結果

FP 分類流量分析

False Positive Classification Flow — Root Cause Analysis

全部 FP

100%

真值集合 vs ISM輸出，標記為FP的variants

All FP variants

raw_absent

98.6%

ISM覆蓋讀序 < 5條，無法觀測

ISM cannot see

Germline ASM

94%

matched normal中也存在相同甲基化差異

Germline ASM origin

HKU 特有

88%

僅在HCC1395_HKU中出現的ASM

HKU-specific ASM

無法過濾

~88%

不符合以上任何標準，ISM無法識別

Cannot filter

TO 策略的效果（秉庭分析）

HCC1395 樣本：ΔF1 = +0.0094

甲基化特徵補充 normal background
(methylation features supplement normal BG)

Germline ASM 仍為主要雜訊來源
(Germline ASM remains primary noise)

→ FN rescue 比 FP filter 更有效
→ FN rescue > FP filter strategy

TO 可介入：Layer 3（Germline ASM 建模）

關鍵洞察

98.6% raw_absent：由 ISM 覆蓋分析統計
（非 PON 過濾）——ISM 嘗試分析 FP 位點
時，98.6% 的位點在 BAM 中讀序 < 5 條。

94% Germline ASM：另一來源，由 matched
normal 比對確認（非同一統計）。

百分比 = 佔全部FP的比例 | 98.6% 為ISM覆蓋統計；94% 為matched normal比對統計（兩者獨立）

raw_absent：SNV 位點在 BAM 中無 reads 覆蓋，ISM 無法觀測

FN rescue：透過額外特徵恢復被漏掉的真陽性（提高 Recall）

> **Notes:** S18 FP流量分析。五層+標準欄位。右側TO策略+關鍵洞察。raw_absent/FN rescue術語首次出現。

---

## Slide 19: S19 ｜ 實驗結果

S19 ｜ 實驗結果

已確認正面發現

Confirmed Positive Findings — Validated Results

1

✓ HP1-1 甲基化可區分 TP vs FP

HP1-1 methylation discriminates TP from FP

為何測試

Why test

核心假設：體細胞 SNV 會改變鄰近 CpG 甲基化格局

數據

Evidence

AUROC = 0.85（LOH≥0.8 樣本），一致性方向在 6/7 樣本

推斷

Inference

TP reads 的甲基化有真實可區分信號 → 確認假設 A 成立。支持 Phase 1 繼續優化 Bernoulli 距離參數

2

✓ LOH 程度與 AUROC 正相關

LOH level positively correlates with AUROC

為何測試

Why test

預期：LOH 純化 HP1 背景，降低 germline 干擾，信號更清晰

數據

Evidence

LOH≥0.8 → AUROC 0.88；LOH<0.4 → AUROC 0.65

推斷

Inference

LOH 是最重要的分層變量 → Phase 2 purity-aware clustering 有直接理論支撐；高LOH場景為優先應用

3

✓ VerificationClass 精確度層次分明

VerificationClass shows clear precision stratification

為何測試

Why test

驗證：雙向顯著性框架是否確實能區分事件可信度

數據

Evidence

Strong~88%, Weak~64%, Noise<10%（HCC1395 全量統計）

推斷

Inference

Strong VC 可作為高信心 somatic methylation 指標直接輸入下游整合；Noise<10% 證明分類穩健

4

✓ TO 策略在 HCC1395 有顯著 ΔF1

Tumor-only strategy shows significant F1 gain on HCC1395

為何測試

Why test

秉庭分析：缺乏 normal sample 時，甲基化特徵能否補足缺失信號

數據

Evidence

ΔF1 = +0.0094，統計顯著（p < 0.05），HCC1395_HKU

推斷

Inference

TO 場景有真實應用潛力，但 HCC1395_HKU 的 germline ASM 偏差（94%）是主要干擾 → 需更多樣本確認 p 值穩定性

> **Notes:** S19 四項正面發現。每項有：勾標記+標題，及為何測試/數據/推斷三欄說明。

---

## Slide 20: S20 ｜ 實驗結果

S20 ｜ 實驗結果

FP 根因：Germline ASM 干擾

FP Root Cause — Germline ASM Interference

機制說明

正常細胞

HP1

HP2

腫瘤化（Germline ASM 保留）

腫瘤細胞（Germline ASM 模式保留）

HP1-1

HP1

HP2

Germline ASM 甲基化模式在腫瘤中保留
→ HP1/HP1-1 甲基化差異來自 Germline，非 Somatic

量化數據（HCC1395_HKU）

98.7%

FP 受 Germline ASM 影響
FPs affected by Germline ASM

0.155

HP1 vs HP2 甲基化差異中位數
Median methylation diff (HP1 vs HP2)

94%

HKU 樣本中 ASM 影響比例
ASM-affected fraction in HKU

18/19

Chr 標記 ASM 位點
ASM-marked chromosomal loci

98.7% 量化方法：比對7樣本全量FP中，matched normal也顯示相同CpG位點甲基化差異（|Δmeth|≥0.15）
的位點比例。HCC1395_HKU樣本的數據，排除了VAF/AD過低和raw_absent後的有效FP分析結果。

解決方向

使用 normal ASM 資料庫（CoLoRSdb）排除 Germline ASM
Use CoLoRSdb to exclude Germline ASM loci

FN rescue：從 ISM 輸出 strong signal 協助恢復真陽性
FN rescue: strong ISM signals boost True Positive recall

論文定位：ISM 不是 filter，是 epigenetic context provider
Paper: ISM as epigenetic context provider, not filter

Germline ASM：胚系等位基因特異性甲基化，在所有細胞中保留

CoLoRSdb：正常長讀長甲基化變異資料庫（Phase 1 參考）

> **Notes:** S20 FP根因Germline ASM。左側機制圖，右側4個關鍵數字與解決方向。

---

## Slide 21: S21 ｜ 實驗結果

S21 ｜ 實驗結果

已排除方向與當前最佳策略

Ruled-Out Directions & Current Best Strategy

已排除方向（非無效嘗試）

✗

PassedGating 作為 FP 過濾器

[為何測試]

測試 ISM gating 是否能直接排除假陽性

[排除原因]

反效果：HCC1395_5kHz 候選池偏差，全量統計無效果

[重要發現]

揭示了候選池選擇偏差問題，確認 FP 無法透過 ISM 過濾

✗

Jaccard 距離作為主要度量

[為何測試]

Jaccard 適合二值化甲基化，理論上匹配

[排除原因]

連續甲基化概率下，Bernoulli 顯著更優（AUROC 差 0.08）

[重要發現]

排除了次優方案，確立 Bernoulli 為最佳選擇

✗

Per-window 整體甲基化差異

[為何測試]

窗口均值可能比 per-read 更穩健

[排除原因]

失去 Read-level 異質性資訊，掩蓋亞克隆結構

[重要發現]

確認 Read-level 分析不可降維至窗口均值

當前最佳策略

E → A

高信心 TP 擴充

Strong signal → FN rescue
AUROC 0.85 閾值篩選

Strong signal → FN rescue
Filter by AUROC 0.85 threshold

A → D

Germline ASM 排除

CoLoRSdb + 正常樣本 ASM DB
建立 Germline 背景模型

CoLoRSdb + normal sample ASM DB
Build Germline background model

D → B

多 SNV 合併

柏宇任務：2 SNV within 10kb
PS block 合併分析

Task: 2 SNVs within 10kb
Merge PS block analysis

B → C

論文定位確立

ISM = epigenetic context
不是 variant filter

ISM = epigenetic context provider
not a variant filter

> **Notes:** S21 排除方向與策略。左側3個排除項，右側E→A→D→B→C最佳策略路徑。

---

## Slide 22: S22 ｜ 文獻對照

S22 ｜ 文獻對照

文獻 1：Long-Read POG（Cell Genomics 2024）

Literature: Long-Read Pan-Cancer Analysis & IMPALA Integration

Long-Read POG

期刊：

Cell Genomics, 2024

重點：

IMPALA：整合 CNV + aDMR + ASE

方法：

Long-read ONT 全基因組 pan-cancer

發現：

甲基化+CNV+基因表達三維整合
可提升腫瘤 clonal 解析度

一行重點：IMPALA 整合 CNV+aDMR+ASE，長讀長分析可達基因層級 clonal 解析

DOI: 10.1016/j.xgen.2024.100XXXX

ISM 啟示

Phase 3 目標

ISM 可作為 gene-level evidence integrator
整合甲基化+HP+CNV 資訊

技術啟發

aDMR 概念延伸到 ISM 的 Read-level ASM
CpG cluster 比 single-site 更穩健

論文定位

ISM 補充 IMPALA 缺少的 Read-level
亞克隆分辨能力

aDMR：Allele-specific Differentially Methylated Region

> **Notes:** S22 Long-Read POG。左側論文資訊+一行重點，右側ISM啟示：Phase 3整合器定位。

---

## Slide 23: S23 ｜ 文獻對照

S23 ｜ 文獻對照

文獻 2：TRACERx CAMDAC（Nature Genetics 2025）

Literature: Purity-Aware CN & LOH-SNV Integration

TRACERx CAMDAC

期刊：

Nature Genetics, 2025

重點：

purity/CN 校正 + LOH-SNV 聯合驗證

方法：

甲基化陣列 + WGS，TRACERx 縱向隊列

發現：

腫瘤純度校正後，LOH 位點甲基化
異常可作為克隆狀態指標

一行重點：purity/CN 校正 + LOH-SNV 驗證，提供 ISM Phase 2 的方法學模板

DOI: 10.1038/s41588-025-XXXX

ISM 啟示

Phase 2 方法學

採用 purity-aware clustering：
依腫瘤純度調整 HP1-1 threshold
避免低純度樣本誤判

LOH 驗證啟發

CAMDAC LOH-SNV 聯合驗證
→ ISM 加入 LOH ratio 作為
甲基化分析的前置過濾

資料整合

TRACERx 縱向設計 → ISM 未來
可延伸到時序分析（Second Hit 時間點）

purity-aware：依腫瘤純度校正，避免正常細胞污染稀釋信號

> **Notes:** S23 CAMDAC。purity/CN校正+LOH-SNV驗證 → ISM Phase 2方法學模板。

---

## Slide 24: S24 ｜ 文獻對照

S24 ｜ 文獻對照

文獻 3：EVOFLUx + MethylBERT

Literature: fCpG Evolutionary Barcodes & Transformer Read Classification

EVOFLUx — fCpG 演化條碼

重點：fCpGs（功能性 CpG）作為腫瘤演化條碼
可追蹤克隆演化路徑，區分 somatic vs germline 甲基化

ISM 啟示：Phase 4 CpG 篩選——以 fCpG 優先，濾除 Germline ASM 干擾

EVOFLUx DOI: 10.1038/XXXX

MethylBERT — Transformer 讀序分類

重點：Transformer 直接對 methylation string 編碼
無需手工設計特徵，直接學習讀序甲基化模式

ISM 啟示：Phase 1 ML 替代方案——以 MethylBERT 取代距離+UPGMA

MethylBERT DOI: 10.1038/XXXX

整合啟示：ISM 四階段文獻對應

Phase 1

MethylBERT 替代距離+聚類
→ 端對端讀序分類

Phase 2

CAMDAC purity-aware 框架
→ 依純度校正分群

Phase 3

IMPALA 基因層級整合
→ ISM + CNV + ASE

Phase 4

EVOFLUx fCpG 篩選
→ Germline ASM 排除

fCpG：功能性 CpG：受轉錄因子調控，演化上保守

MethylBERT：甲基化字串的 BERT Transformer 分類器

> **Notes:** S24 EVOFLUx+MethylBERT。左側兩篇論文各有重點+ISM啟示，右側四階段文獻對應。

---

## Slide 25: S25 ｜ 未來方向

S25 ｜ 未來方向

ISM 定位 vs 競爭研究工具

ISM Positioning — Competitive Analysis

工具比較矩陣

ISM（本研究）

IMPALA

CAMDAC

MethylBERT

DeepSomatic

Read-level 甲基化

✓

✓

✗

✓

✗

HP 單倍型整合

✓

✗

✗

✗

✗

LOH 偵測

✓

✓

✓

✗

✗

亞克隆分解

✓

✗

✓

✗

✗

Somatic SNV 驗證

✓

✗

✓

✗

✓

長讀長 (ONT)

✓

✓

✗

✓

✗

ISM 組合獨特性  |  Why ISM is Unique

唯一整合三維度

Read-level 甲基化＋HP 單倍型＋LOH 偵測 → 三維度同步分析
Unique 3D integration: read-level methylation + haplotype + LOH

SNV 視窗聚焦

不是全基因組，而是 SNV 周邊 ±2500 bp 精準窗口的甲基化格局
Focused on ±2500 bp SNV window methylation, not whole genome

VerificationClass

系統化 4 類輸出（Strong / Weak / Subclone / Noise），提供信心層級
Systematic 4-tier output with confidence stratification

論文定位差異化

非 variant filter，而是 epigenetic context provider（讀序層級語境）
Not a variant filter — epigenetic context provider

> **Notes:** S25 ISM定位。左側工具比較矩陣5工具×6特徵，底部橫向四卡片ISM組合獨特性說明。

---

## Slide 26: S26 ｜ 未來方向

S26 ｜ 未來方向

四階段研究策略概覽

Four-Phase Research Strategy — Overview

Phase 1

讀序分類優化

Read Classification

2026 Q2

文獻：MethylBERT
CoLoRSdb

• Bernoulli 距離確立
• MethylBERT ML 替代方案
• VerificationClass 精準化
• CoLoRSdb Germline ASM 排除

Phase 2

克隆結構解析

Clone Architecture

2026 Q3

文獻：CAMDAC
TRACERx

• Purity-aware clustering
• 多 SNV 合併（柏宇任務）
• LOH-SNV 聯合分析
• Second Hit 時序推斷

Phase 3

基因層級整合

Gene-level Integration

2026 Q4

文獻：Long-Read POG
IMPALA

• ISM + CNV + ASE 整合
• Gene-level evidence 彙總
• 腫瘤 clonal 解析報告
• API 設計與輸出格式

Phase 4

論文完成

Paper Submission

2027 Q1

文獻：EVOFLUx
SEQC2

• fCpG 篩選（EVOFLUx）
• 完整 benchmark 7樣本
• ISM 定位確立
• 期刊投稿

← 你在這裡

> **Notes:** S26 四階段策略概覽。水平時間軸Phase1→4，上方卡片顯示目標/文獻，下方列出具體任務。

---

## Slide 27: S27 ｜ 未來方向

S27 ｜ 未來方向

Phase 1 & 2 詳細規劃

Phase 1 & 2 — Detailed Implementation Plan

Phase 1：讀序分類優化

Read Classification Enhancement

時程：

2026 Q2（3個月）

文獻：

MethylBERT + CoLoRSdb

負責：

廖子游

ISM 實作方式

1

Bernoulli 距離 → 確認最佳化參數

2

CoLoRSdb：Germline ASM 資料庫整合

3

VerificationClass 門檻重新校準（7樣本）

4

MethylBERT：作為 Bernoulli 的 ML 對照組

預期收益

FP率 降低 >50%（Germline ASM 排除後）
AUROC 平均提升至 0.82+
FP rate reduced >50% (after Germline ASM exclusion); avg AUROC → 0.82+

Phase 2：克隆結構解析

Clonal Architecture Resolution

時程：

2026 Q3（3個月）

文獻：

CAMDAC + TRACERx

負責：

柏宇（多SNV任務）+ 廖子游

ISM 實作方式

1

多 SNV 合併：2 SNV within 10kb + PS block

2

Purity-aware clustering（依腫瘤純度調整）

3

LOH-SNV 聯合驗證框架（參考 CAMDAC）

4

Second Hit 時序推斷 → 甲基化事件排序

預期收益

亞克隆結構識別率 >60%
Subclone VC 類別精確度提升
Subclone identification rate >60%; Subclone VC precision improved

> **Notes:** S27 Phase 1&2詳細。兩欄：Phase1（廖子游，讀序分類，Q2）Phase2（柏宇+廖子游，克隆結構，Q3）。

---

## Slide 28: S28 ｜ 未來方向

S28 ｜ 未來方向

Phase 3 & 4 + 論文定位

Phase 3 & 4 Detail + Paper Positioning

Phase 3：基因層級整合

Gene-level Evidence Integration  |  2026 Q4

1. ISM + CNV + ASE 三維整合框架
3D integration framework

2. Gene-level clonal evidence 彙總
Gene-level evidence aggregation

3. 整合 IMPALA 輸出格式
Integrate IMPALA output format

4. API 設計：下游工具接口
API design for downstream tools

Phase 4：論文完成與投稿

Manuscript & Submission  |  2027 Q1

1. fCpG 篩選（參考 EVOFLUx）
fCpG filtering (inspired by EVOFLUx)

2. 完整 7 樣本 benchmark
Full 7-sample benchmark

3. ISM 工具文件化與 GitHub 發布
ISM documentation & GitHub release

4. 期刊投稿（目標：Nature Methods）
Submission (target: Nature Methods)

★ 論文定位宣言

✗ ISM 不是什麼

• Variant filter（直接排除假陽性）
  (not a direct FP eliminator)
• 精確度過高的 SNV caller
  (not an ultra-precise SNV caller)
• 替代 ClairS 或 DeepSomatic
  (not a replacement for ClairS/DeepSomatic)

→

✓ ISM 是什麼

• Read-level epigenetic context provider
• 為每個 SNV 提供甲基化語境（HP + CpG + LOH）
  Provides per-SNV methylation context (HP + CpG + LOH)
• 協助下游工具提升 somatic calling 信心
  Boosts somatic calling confidence for downstream tools
• 揭示腫瘤亞克隆結構的甲基化特徵
  Reveals methylation signatures of tumor subclonal structure

> **Notes:** S28 Phase3&4+論文定位。上半兩欄Phase3/4，下半橘色邊框論文定位宣言：NOT filter，IS context provider。

---

## Slide 29: S29 ｜ 未來方向

S29 ｜ 未來方向

研究路線圖（三軌）

Research Roadmap — Three Parallel Tracks

2025-12

2026-01

2026-02

2026-03

2026-Q2

2026-Q3

2026-Q4

2027-Q1

← 你現在在這裡

① 個人實驗進度
Personal Experiments

ISM
概念設計

✓

C++ 核心
implementation

✓

7樣本
全量分析

✓

★ FP根因
+AUROC

Phase 1
Germline ASM

Phase 2
柏宇多SNV

Phase 3
整合框架

口試與
論文撰寫

② 外部文獻對照
External Literature

MethylBERT
(preprint)

✓

EVOFLUx
(bioRxiv)

✓

★ 開始
文獻整合

CAMDAC
(Nat Genet)

Long-Read
POG (CellGen)

系統性
benchmark

論文
撰寫

③ 程式開發進度
Code Development

Python
pptx腳本

✓

C++ HTSlib
parser

✓

OpenMP
平行化

✓

★ VC class
+QualScore

CoLoRSdb
integration

Multi-SNV
merge

API
設計

GitHub
public

→ 下一步：Phase 1 Germline ASM 排除（Q2 2026）｜ 柏宇啟動多 SNV 合併任務

> **Notes:** S29 三軌路線圖。橫軸時間2025-12到2027-Q1，縱軸三條軌道：個人實驗/文獻/程式開發。2026-03標記你在這裡，已完成事項有✓標記，下一步橘色條。

---

## Slide 30: S30 ｜ 未來方向

S30 ｜ 未來方向

總結：五個目標進度

Summary — Progress on Five ISM Objectives

對應 S04 五個目標 | Current status of each objective

1

CpG 甲基化關聯分析
Methylation-SNV Association

70%

部分完成

✓ AUROC > 0.75（LOH分層樣本）  ✓ 6/7 樣本甲基化方向一致

→ Germline ASM 排除（Phase 1）

2

癌症克隆狀態
Cancer Clone Architecture

25%

進行中

✓ VerificationClass Subclone 類別建立

→ 柏宇任務：多 SNV 合併（Phase 2）

3

Second Hit 時序推斷
Second Hit & Event Ordering

30%

初步框架

✓ LOH 偵測邏輯確立  ✓ 甲基化保留效應觀察

→ HP1-LOH 聯合時序模型（Phase 2）

4

TO / Purity 場景
Tumor-Only Scenarios

20%

秉庭分析中

✓ ΔF1 = +0.0094（HCC1395 樣本）

→ 更多樣本驗證 TO 策略普遍性

5

整合提升 F1
Somatic F1 Boost

15%

方向確立

✓ FN rescue 策略確認（取代FP filter）

→ FN rescue + Germline ASM DB（Phase 1+2）

> **Notes:** S30 總結。五個目標各有進度條+已完成+下一步，對應S04目標清單，形成首尾呼應。

---

## Slide 31: 謝謝聆聽

謝謝聆聽

Thank You for Your Attention

Q & A

廖子游

中正大學資工系 黃耀廷 教授 Lab405 實驗室

InterSubMod v3 | 2026-03-26

今日報告重點

✓  ISM 五個目標框架確立

✓  AUROC > 0.75（LOH分層）

✓  FP根因：Germline ASM（98.7%）

✓  ΔF1 +0.0094（秉庭 TO 分析）

→  Phase 1–4 路線圖規劃完成

InterSubMod · Inter-Subclonal Methylation Analysis · ONT Long-Read

> **Notes:** S31 Q&A感謝頁。上半深藍底謝謝，下半HP1藍底聯絡資訊+今日重點回顧。

---
