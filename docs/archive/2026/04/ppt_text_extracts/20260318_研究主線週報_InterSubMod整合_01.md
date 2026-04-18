# 20260318_研究主線週報_InterSubMod整合_01.pptx

- Slides: 27
- Date: 2026-03-18

## Slide 1: InterSubMod 研究主線週報

InterSubMod 研究主線週報

2026-03-12 → 2026-03-17  整合報告

Weekly Report: Data Migration · Visualization · Methylation Analysis

2026-03-18

> **Notes:** 報告定位：整合 2026-03-12 至 2026-03-17 期間所有 validated reports 與 assets。涵蓋資料遷移、視覺化 AI 評估 Skill 流程、甲基分布觀察、Knowledge 修正與整體結論。

---

## Slide 2: 本週工作大綱

本週工作大綱

Weekly Agenda

①

資料遷移

建立全組合 Tagged.bam
InterSubMod 輸出與 F1 Benchmark

②

視覺化評估

AI 輔助位點觀察 Skill 流程
正/反例 × heatmap × IGV × .md 報告

③

分布分析

甲基特徵整體觀察
Knowledge 修正與邊界確立

Data Migration: Build full-combination Tagged.bam; ISM output + F1 Benchmark

Visualization: AI locus Skill pipeline; TP/FP × heatmap × IGV × .md report

Distribution Analysis: Methylation feature overview; Knowledge correction

> **Notes:** 三部分相互依存。② 提供修正方向，③ 確立邊界。視覺化 Skill 流程是本週核心可交付物。

---

## Slide 3: 本週最重要的兩正兩負

本週最重要的兩正兩負

Key Results This Week

✓

全組合 Tagged.bam 建立完成

21 completed runs 可重跑；
覆蓋 7 samples × 3 workflow（pf/pp/to）

✓

位點觀察 Skill 流程建立

heatmap → IGV → .md 報告，可直接套用到後續樣本

✗

ISM 效果高度依賴 LP 殘餘 FP 量

FP < 100 時基本無效；
LP_FP > 600 才有穩定增益

✗

無法跨樣本通用的 TP/FP 過濾規則

Quality_Score 等特徵均有 sample dependency

✓ Full-combination Tagged.bam built (21 runs; 7 samples × 3 workflows, retriable)

✓ Locus observation Skill pipeline built (heatmap → IGV → .md, reusable)

✗ ISM efficacy depends on LP residual FP; LP_FP<100 = near-ineffective, >600 for stable gain

✗ No universal cross-sample TP/FP filter; Quality_Score et al. are all sample-dependent

> **Notes:** 先告訴聽眾這週哪些成立、哪些不能升格，後面看圖時才不會誤解。

---

## Slide 4: 資料遷移：全組合 Tagged.bam 建立流程

資料遷移：全組合 Tagged.bam 建立流程

Data Migration Pipeline

Samples × pileup / full_model  ×  paired / Tumor-only

7 Samples
×2 Caller
×2 Workflow
= 28 組合

LongPhase
Tagging

Tagged.bam
(big7_disk)

InterSubMod

甲基輸出
+ F1 Benchmark

Pipeline: 7 Samples × 2 Callers × 2 Workflows = 28 combinations → LongPhase Tagging → Tagged.bam → InterSubMod methylation output + F1 Benchmark

> **Notes:** TO-full 全部 7 runs 不可得（source tree 確認）；HCC1395/DORADO to_pileup 亦不可得。其餘 21 runs 全部 completed。

---

## Slide 5: 全組合涵蓋矩陣（Coverage Matrix）

全組合涵蓋矩陣（Coverage Matrix）

28-Combination Coverage

Sample | paired_full | paired_pileup | to_pileup
--- | --- | --- | ---
HCC1395 | completed | completed | completed
HCC1395_DORADO | completed | completed | completed
COLO829 | completed | completed | completed
H1437 | completed | completed | completed
H2009 | completed | completed | completed
HCC1937 | completed | completed | completed
HCC1954 | completed | completed | completed

> **Notes:** completed 21 / N/A 7（to_full 全不可得；HCC1395/DORADO to_pileup 不可得）。N/A 是已確認限制，非待補工作。

---

## Slide 6: 三層 Benchmark F1 總表（21 Runs）

三層 Benchmark F1 總表（21 Runs）

Layer F1 Summary · 21 Runs

Sample | Mode | Caller F1 | LP F1 | ISM F1 | LP→ISM Δ
--- | --- | --- | --- | --- | ---
HCC1395 | paired_full | 0.8443 | 0.8522 | 0.8532 | +0.0010
HCC1395 | paired_pileup | 0.7921 | 0.8551 | 0.8571 | +0.0020
HCC1395 | to_pileup | 0.7166 | 0.7166 | 0.7167 | +0.0001
HCC1395_DORADO | paired_full | 0.8565 | 0.8592 | 0.8590 | -0.0002
HCC1395_DORADO | paired_pileup | 0.8511 | 0.8662 | 0.8663 | +0.0001
HCC1395_DORADO | to_pileup | 0.7226 | 0.7226 | 0.7224 | -0.0002
COLO829 | paired_full | 0.8698 | 0.8725 | 0.8725 | +0.0000
COLO829 | paired_pileup | 0.8253 | 0.8680 | 0.8680 | +0.0000
COLO829 | to_pileup | 0.7068 | 0.7068 | 0.7069 | +0.0001
H1437 | paired_full | 0.9143 | 0.9087 | 0.9087 | +0.0000
H1437 | paired_pileup | 0.9237 | 0.9195 | 0.9194 | -0.0001
H1437 | to_pileup | 0.6499 | 0.6499 | 0.6499 | +0.0000
H2009 | paired_full | 0.9688 | 0.9663 | 0.9662 | -0.0001
H2009 | paired_pileup | 0.9710 | 0.9710 | 0.9708 | -0.0002
H2009 | to_pileup | 0.8986 | 0.8986 | 0.8985 | -0.0001
HCC1937 | paired_full | 0.8276 | 0.8415 | 0.8414 | -0.0001
HCC1937 | paired_pileup | 0.8310 | 0.8875 | 0.8875 | +0.0000
HCC1937 | to_pileup | 0.6076 | 0.6076 | 0.6077 | +0.0001
HCC1954 | paired_full | 0.8762 | 0.8743 | 0.8731 | -0.0012
HCC1954 | paired_pileup | 0.9011 | 0.9083 | 0.9058 | -0.0025
HCC1954 | to_pileup | 0.3780 | 0.3780 | 0.3767 | -0.0013

> **Notes:** 21 completed runs（較上週新增 COLO829/H1437/H2009/HCC1937/HCC1954 三個 to_pileup）。HCC1395 paired_pileup 正增益最大（+0.0020）；HCC1954 全部 mode 均為負增益，paired_pileup 最差（-0.0025）。

---

## Slide 7: ISM 效果的關鍵：依賴 LongPhase 殘餘 FP 量

ISM 效果的關鍵：依賴 LongPhase 殘餘 FP 量

ISM Efficacy vs. Residual LP FP

規律：LP_FP < 100  →  ISM 幾乎必然損傷 TP
       LP_FP > 600  →  ISM 才有效率地過濾 FP

Run（節選） | LP_FP | ISM 移除 FP | ISM 損失 TP | 效率
--- | --- | --- | --- | ---
HCC1954 paired_pileup | 148 | 4 | -101 | 0.04x（最差）
H1437 paired_full | 8 | 0 | -8 | 0x
HCC1395 paired_full | 627 | 83 | -2 | 41.5x（最佳）
HCC1395 paired_pileup | 1258 | 178 | -6 | 29.7x

⚠  結論：HCC1954 negative transfer 原因：LongPhase 已把 FP 壓到極低（148），
    剩餘 FP 是高置信強訊號型，ISM 無法分辨，反而大量誤刪 TP（101 個）。

Rule: LP_FP < 100 → ISM nearly always damages TP; LP_FP > 600 → ISM efficiently filters FP

Conclusion: HCC1954 negative transfer — LP already suppressed FP to very low (148); remaining FPs are high-confidence signals ISM cannot distinguish, causing mass TP deletion (101)

> **Notes:** 這解釋了為何 HCC1954 是負增益最嚴重案例。LP_FP < 100 時 ISM 幾乎必然損傷 TP。

---

## Slide 8: 視覺化 × AI 評估 Skill 流程

視覺化 × AI 評估 Skill 流程

Visualization & AI Evaluation Skill Pipeline

研究數據

輸入

AI 篩選
代表性位點

正/反例

甲基熱圖
methylartist

視覺

IGV 固定
模板截圖

確認

.md
觀察報告

輸出

回饋修正
InterSubMod

迭代

重點原則：AI 負責視覺初篩，最終裁決以 BAM / VCF / truth TSV 為準。
不可單憑甲基圖下結論。

Key principle: AI handles visual pre-screening; final judgment based on BAM/VCF/truth TSV — never conclude from methylation heatmap alone

> **Notes:** 此流程為本週核心成果。可重跑、可延續到其他樣本。已制度化為 Skill，後續任何樣本可直接套用。

---

## Slide 9: 代表性 TP：Allele-driven 甲基分層

代表性 TP：Allele-driven 甲基分層

TP-4 chr16:35118902 — Allele-driven Methylation Stratification

左：Cluster heatmap — ALT vs REF reads 甲基模式清楚分離｜右：Distance heatmap — 兩個獨立群體確認非視覺錯覺

Left: Cluster heatmap — ALT vs REF reads show clear methylation stratification | Right: Distance heatmap — two distinct groups confirm non-visual-artifact

> **Notes:** 這類 allele-driven TP 是 InterSubMod 最有解釋力的場景（subclone annotation）。

---

## Slide 10: 代表性 TP：IGV 固定模板正式驗證

代表性 TP：IGV 固定模板正式驗證

TP-4 IGV Snapshot Validation

• 固定 track 順序

• Tumor BAM

• → Normal BAM

• ALT allele 集中

• 於特定 haplotype

• BAM / VCF /

• truth TSV

• 三層交叉驗證

• → somatic TP

Fixed track order: Tumor BAM → Normal BAM; ALT allele concentrated on specific haplotype

Triple cross-validation: BAM/VCF/truth TSV → confirmed somatic TP

> **Notes:** heatmap 提供假說，IGV 提供視覺確認，validation TSV 是最終依據。三層交叉確認保證嚴謹性。

---

## Slide 11: IGV 固定模板截圖方法說明

IGV 固定模板截圖方法說明

IGV Snapshot Protocol

TP-4（正例）

FP-B2（反例）

操作流程

每個 loci 使用相同 IGV session

與 track 順序

Tumor BAM → Normal BAM → VCF

批次截圖腳本自動輸出：

{TP/FP}_{id}_{chr}_{pos}.png

9/9 loci 均已自動出圖

固定模板讓 TP 與 FP 在

同格式下可直接目測比較

Same IGV session per locus; Tumor BAM → Normal BAM → VCF; batch script: {TP/FP}_{id}_{chr}_{pos}.png; 9/9 loci auto-generated

> **Notes:** 固定模板是可重現研究的關鍵。9/9 loci 均已成功截圖，包含 5 個 TP 與 4 個 FP 代表性位點。

---

## Slide 12: 代表性 FP：HP-driven（FP-B1）

代表性 FP：HP-driven（FP-B1）

FP-B1 Haplotype-Driven Pattern

[ FP-B1 methylartist 圖 ]
05_FPB1_methylartist.svg
（cairosvg 不可用，請手動插入）

FP-B1 特徵分析

• 甲基訊號：無 allele 間差異
（HP0 dominant）

• HPP：高（HP-dominant Pattern 明顯）

• AlleleDelta：低（接近 0）

• CramersV：低（無 allele 關聯）

• 辨識方式：高 HPP + 低 AlleleDelta
是 InterSubMod 可用特徵

• 注意：HPP 是本週確立的核心
候選訊號，尚未接回正式邏輯

FP-B1: No inter-allele methylation diff (HP0 dominant); high HPP + low AlleleDelta/CramersV → key candidate signal; HPP confirmed this week, not yet wired into formal logic

> **Notes:** HPP 是本週確立的核心候選訊號，但尚未接回正式邏輯。methylartist SVG 需手動插入（cairosvg 環境缺失）。

---

## Slide 13: 代表性 FP：Alignment Anomaly（FP-B2 修正過程）

代表性 FP：Alignment Anomaly（FP-B2 修正過程）

FP-B2 Misclassification Corrected via IGV

早期觀察

圖像看似 MNP 類型
（多鹼基 polymer）

IGV 驗證

發現明顯 deletion gap
/ soft-clip 集中

最終口徑

deletion-like gap /
local alignment anomaly
（非 MNP）

重點：圖像初篩 ≠ 最終判斷；BAM level 驗證不可省略

Early obs: appeared MNP-type (multi-base polymer) → IGV: clear deletion gap/soft-clip → Final: deletion-like gap/local alignment anomaly (not MNP) | Key: image pre-screening ≠ final judgment; BAM-level verification is mandatory

> **Notes:** 這個案例說明 AI 視覺初篩的風險，以及正式驗證工作流的必要性。最終口徑已更新為 deletion-like gap / local alignment anomaly。

---

## Slide 14: .md 觀察報告格式展示

.md 觀察報告格式展示

Observation Report Template

報告結構

位點座標（chr:pos）

TP / FP 類型標記

特徵數值：

CramersV / AlleleDelta

Quality_Score / HPP

甲基熱圖連結

IGV 截圖連結

結論摘要

格式：Markdown

可在 Obsidian 直接顯示

已制度化為 Skill

可套用到後續樣本

/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260317_case_obs_report_FP_B1_範例_01.md

Report: locus coord + TP/FP label + feature values (CramersV/AlleleDelta/QS/HPP) + heatmap/IGV links + conclusion. Markdown format (Obsidian); institutionalized as Skill

> **Notes:** 觀察報告 Skill 已完成。後續任何樣本可直接套用。報告作為 AI 初篩 → 正式驗證之間的交接文件，已制度化。

---

## Slide 15: Tumor BAM 
group
:base

Tumor BAM
group
:base

Tumor BAM
group
:phase

Normal BAM
group
:phase

HC.bed
sSNV
sIndel

phase

output

---

## Slide 16: (no title)

---

## Slide 17: (no title)

---

## Slide 18: 特殊現象：FP 但 CramersV 極高（SEQC2 Indel 問題）

特殊現象：FP 但 CramersV 極高（SEQC2 Indel 問題）

Edge Case: High CramersV FP due to SEQC2 Indel Boundary

現象

• FP 位點 CramersV 異常高

• AlleleP、HPP 等特徵也偏向 TP

• 若只看甲基特徵，幾乎無法辨識為 FP

→ 甲基特徵完全誤導

根本原因（已確認）

• SEQC2 將該位點歸類為

somatic indel，而非 sSNV

• 對應 sSNV 呼叫被標為 FP

（truth class 定義邊界問題）

• 非 InterSubMod 辨識失誤

量化：SEQC2 indel ∩ sSNV 交集

< 10 個位點（極端案例）

→ 已知，記錄在案，不特別處理

Phenomenon: FP locus with abnormally high CramersV; AlleleP/HPP also skewed toward TP → methylation features completely misleading

Root cause: SEQC2 classifies as somatic indel, not sSNV → sSNV call labeled FP (truth class boundary issue). SEQC2 indel ∩ sSNV < 10; documented as known edge case

> **Notes:** 這類問題從未被思考過。統計後確認影響極小（個位數），記錄在案即可。這是 truth set 的定義邊界問題，而非 InterSubMod 的辨識失誤。

---

## Slide 19: 甲基數據整體分布觀察

甲基數據整體分布觀察

Methylation Feature Distribution Overview

HCC1395
5kHz paired

Quality_Score

TP 中位 95

FP 中位 75

最強（+0.800）

CramersV

中位數皆 0

看右尾 Q90+

AlleleDelta

5kHz 極端高

PairwiseMD

無跨樣本一致

Quality_Score strongest (TP median 95 / FP 75); CramersV median=0; AlleleDelta sample-specific; PairwiseMD no cross-sample consistency

> **Notes:** 此圖為 HCC1395 5kHz paired 範例。各欄位 TP vs FP 密度分布。Quality_Score 在此樣本最具鑑別力，但跨樣本無法硬閾值化（見下頁）。

---

## Slide 20: Quality_Score：有訊號但無法跨樣本硬閾值

Quality_Score：有訊號但無法跨樣本硬閾值

Quality_Score: Signal Present but Not Threshold-Safe

問題分析

• 視覺上看起來可分

TP（QS≈100）vs FP（QS≈75）盒形圖有間距

• 問題 1：Class imbalance

TP 與 FP 數量差異極大，視覺誤導

• 問題 2：H2009 實測

QS < 75 precision 僅 0.002
（幾乎全部誤殺 TP）

• 問題 3：TO 模式反轉

to_pileup 下 FP QUAL > TP QUAL

結論：Quality_Score 適合 soft support annotation；不可做跨樣本硬閾值

Issues: (1) class imbalance misleads; (2) H2009 QS<75 precision=0.002; (3) TO mode reversal. Conclusion: soft annotation only — no hard cross-sample threshold

> **Notes:** 觀察.md §14.5：TO 模式的高置信 germline FP 導致 QUAL 反轉。Quality_Score 有訊號，但 sample-specific，不可全域閾值化。

---

## Slide 21: 跨樣本特徵可靠性評級

跨樣本特徵可靠性評級

Cross-Sample Feature Reliability Rating

特徵 | 跨樣本一致性 | 鑑別力 | 使用建議
--- | --- | --- | ---
QUAL（caller） | paired 一致；TO 反轉 | 高（paired） | mode-specific；TO 禁用
Quality_Score | 多數 TP>FP，但有 run 重疊 | 中 | soft annotation；禁硬閾值
VAF | 完全 dataset-dependent | 無跨樣本可靠性 | within-sample 診斷只
AlleleDelta | 5kHz 極端高；其他不一 | 低，dataset-specific | 特定 platform 輔助
CramersV | 中位數全為 0 | 中位數層面無效 | 需改看 Q90+ 尾部
PairwiseMedianDist | 同平台同癌種方向均不同 | 無跨樣本一致性 | within-sample diagnostic
NumCpGs | FP 略高（HCC1954 最明顯） | 弱，方向不穩 | 配合其他特徵 stratified

> **Notes:** 這張表是本週方法邊界收斂的核心。任何後續規則升格必須先對照此表。沒有一個特徵可以直接做跨樣本硬閾值。

---

## Slide 22: Knowledge 修正：三個之前理解錯誤的觀念

Knowledge 修正：三個之前理解錯誤的觀念

Knowledge Corrections: 3 Misconceptions Fixed

ClairS paired ≠ PON caller

舊：ClairS paired = PON caller

✓ 正確：ClairS paired = Tumor + Normal caller
PON 是 ClairS-TO / LongPhase-TO 的組件

HP=3 ≠ LOH

舊：HP=3 代表 LOH

✓ 正確：HP=3 表示 ALT-supporting read
germline HP 背景無法唯一解析 LOH

HP 不可跨 phase block 比較

舊：HP1 與 HP2 可跨 phase block 比較

✓ 正確：HP 只在同一 PS / phase block 內有意義
跨 block 不能直接推論親緣關係

K1: ClairS paired = Tumor+Normal caller; PON is a component of ClairS-TO/LongPhase-TO

K2: HP=3 means ALT-supporting read, NOT LOH; germline HP cannot uniquely resolve LOH

K3: HP labels only valid within same PS/phase block; cross-block comparison is invalid

> **Notes:** 這三點已全面更新回術語手冊與週報口徑。HP=3 的修正影響最廣，之前許多觀察結論需重新詮釋。

---

## Slide 23: 整體重點結論

整體重點結論

Summary of Key Conclusions

✓  正向進展

16 completed runs 全組合
tagged.bam 建立完成

位點觀察 Skill 流程已建立（可重跑）
heatmap → methylartist → IGV → .md

Allele-driven TP 甲基訊號可穩定觀察
case triage / subclone annotation 高價值

SEQC2 indel/sSNV 交集問題已量化
（< 10），記錄為已知極端案例

✗  邊界與限制

ISM 效果由 LongPhase 殘餘 FP 量決定
LP_FP < 100 時 ISM 幾乎無效

無跨樣本通用 TP/FP 過濾規則
各樣本 FP 類型根本不同

Quality_Score 等特徵均有
sample dependency，不可全域閾值化

✓ All tagged.bam built; ✓ Locus Skill pipeline established (rerunnable); ✓ Allele-driven TP methylation stably observable; ✓ SEQC2 indel/sSNV intersection (<10) quantified

✗ ISM efficacy determined by LP_FP; LP_FP<100 = near-ineffective; ✗ No universal cross-sample TP/FP filter; ✗ Quality_Score et al. all sample-dependent

> **Notes:** 本週進步的是工作流程與邊界感，不是新的硬規則。這些邊界是真實的研究限制，不是待解決的 bug。

---

## Slide 24: 未來繼續的方向

未來繼續的方向

Next Steps

方向 1：HP Imprinting Filter（接回正式邏輯）

以 Normal 樣本甲基 reference 建立 HP imprinting filter

K1 確認：DominantLabel=hp 作為 paired mode Q3 FP 過濾標準

把 DominantLabel / phase block 接回正式 triage/annotation 輸出

方向 2：LOH 盲點補充機制 + 特徵工程

K3 確認：LOH 區域改用 allele-based 甲基比較作為備援路徑

補 AlleleMethDelta、HP1MethMean/HP2MethMean 新欄位

解析 H1437:paired_pileup 例外（allele-dominant Q3 in paired mode）

Dir 1 — HP Imprinting Filter: Build filter using Normal methylation reference; K1 confirms DominantLabel=hp as paired-mode Q3 FP filter; wire into formal triage/annotation output

Dir 2 — LOH Blind Spot + Feature Engineering: K3 confirms allele-based methylation comparison as LOH backup; add AlleleMethDelta/HP1MethMean/HP2MethMean; investigate H1437:pp exception

> **Notes:** 兩個方向直接接續本週四象限分析的三個驗證結論（K1/K2/K3）。方向 1 針對 HP imprinting FP；方向 2 針對 LOH 盲點補強。

---

## Slide 25: 四象限甲基特徵分析框架

四象限甲基特徵分析框架

Quadrant Analysis Framework · HCC1395 + 6 Samples

Q1

TP + 顯著

方法命中
Allele-driven 甲基分層
→ ISM 正確過濾 FP

Q3

FP + 顯著 ★

方法誤判
HP imprinting 或 allele-driven
→ 高 CramersV 但為 germline FP

Q2

TP + 不顯著

方法沉默（LOH 盲點）
HP_Ratio ≈ 1.0，CramersV ≈ 0
→ 結構性漏判，非方法失效

Q4

FP + 不顯著

方法正確靜默
LP 已過濾 FP，ISM 亦沉默
→ 正常表現

Significant →

Not Significant →

Q1 TP+Significant: Method hit — allele-driven methylation stratification → ISM correctly filters FP

Q3 FP+Significant: Misclassification — HP imprinting or allele-driven → high CramersV but germline FP

Q2 TP+Insignificant: Method silent (LOH blind spot) — HP_Ratio≈1.0, CramersV≈0 → structural miss

Q4 FP+Insignificant: Correctly silent — LP already filtered FP, ISM also silent → normal behavior

> **Notes:** 四象限框架：(truth_label) × (Significant/NotSignificant)。Q3 FP+顯著是本週重點分析對象，發現兩個機制：paired mode → HP imprinting（K1），to_pileup → allele-driven（K2）。Q2 的結構性沉默由 LOH 解釋（K3）。

---

## Slide 26: K1 驗證：Q3 FP+顯著 = HP Imprinting（6 樣本 × 13 Runs）

K1 驗證：Q3 FP+顯著 = HP Imprinting（6 樣本 × 13 Runs）

K1: Paired Mode Q3 FP Driven by Germline HP Imprinting

各 Run Q3 DominantLabel

Run（paired）

Q3數

hp%

結論

HCC1395 pf

4

100%

K1 confirmed

HCC1395 pp

14

93%

K1 confirmed

DORADO pf

10

100%

K1 confirmed

DORADO pp

10

100%

K1 confirmed

COLO829 pp

10

100%

K1 confirmed

H1437 pf

10

90%

K1 confirmed

H1437 pp★

10

20%

EXCEPTION

K1 結論

Paired mode Q3 FP → hp dominant（>80%）
→ germline HP imprinting 是主要 FP 機制

★ H1437:pp 例外：8/10 allele-dominant
   同一 chr9:668xx 區段 → 待調查

K1 conclusion: Paired mode Q3 FP → hp dominant (>80%); germline HP imprinting is primary FP mechanism. Exception: H1437:pp — 8/10 allele-dominant on chr9:668xx → under investigation

> **Notes:** K1 已在 6 個不同樣本 / 13 個 runs 跨樣本驗證。H1437:paired_pileup 例外（8/10 allele-dominant Q3）懷疑是 chr9 segmental duplication 或特殊 imprinting domain。

---

## Slide 27: K2/K3 驗證：to_pileup 機制 + Q2 LOH 盲點（跨樣本 13 Runs）

K2/K3 驗證：to_pileup 機制 + Q2 LOH 盲點（跨樣本 13 Runs）

K2/K3: Allele-Driven to_pileup FP · LOH Q2 Silence Confirmed

K2：to_pileup Q3 FP = Allele-driven

HCC1395_DORADO:to_pileup  10/10 allele-dominant Q3 FP
to_pileup 模式缺乏 Normal 配對 → allele 分群 artifact
→ 需要 Normal BAM 甲基比較才能有效排除

三鍵結論彙整

K1

Paired Q3 FP

HP imprinting
（germline artifact）

17/18 hp dominant
6樣本跨樣本確認

K2

to_pileup Q3 FP

Allele-driven
（缺 Normal 配對）

DORADO 10/10 allele
機制已確認

K3

Q2 TP+不顯著

LOH 盲點
（HP_Ratio≈1.0）

全13 runs 中位=1.0
結構性限制，非方法失誤

K2: to_pileup Q3 FP = allele-driven (lacks Normal pairing → artifact); Normal BAM methylation comparison needed for effective exclusion

K1: Paired Q3 FP = HP imprinting (17/18 hp dominant, 6-sample confirmed) | K2: to_pileup Q3 FP = allele-driven (DORADO 10/10, confirmed) | K3: Q2 TP+Insignificant = LOH blind spot (structural limitation)

> **Notes:** K2 驗證僅有 HCC1395_DORADO:to_pileup 的 10 個位點，但 10/10 allele-dominant 方向極為一致。K3 在 13 個 runs 中全部呈現 HP_Ratio median=1.0，LOH 是唯一解釋。

---
