<!--
建立時間: 2026-04-14 23:30
更新時間: 2026-04-15 (V2 重構)
目標: V2 完整敘事稿 — 每張 slide 講者筆記與口述腳本
處理範圍: 2026-04-14 ~ 2026-04-20
關聯檔案:
  - 00_director_storyboard.md (V2)
  - research/loh_subclone_af/reports/20260414_LOH_Subclone_AF_完整分析報告_01.md
  - docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md
  - docs/presentations/validated/2026/04/20260408_研究週報v3_完整推論鏈與獨立建模/pptx_config.json
  - docs/presentations/validated/2026/04/20260401_LOH_weekly_report_draft/02_ppt_slide_outline.md
數據出處:
  - v5 LOH report (為什麼研究 TO, Slide 2)
  - v3 完整推論鏈 (SEQC2 驗證, 雙定義, filter FAIL)
  - 0413 LOH×CNV 報告 (CN Tier 分層, Coverage_Multiple)
  - validated report §4.1-§4.3 (Step 1-3 結果)
  - validated report §5 (綜合結論)
  - research/loh_subclone_af/manifest.yaml (假說、參數)
-->

# V2 完整敘事稿：LOH Subclone AF × Methylation 雙重證據鏈

## 報告規格

| 項目 | 說明 |
|------|------|
| 主題 | LOH 區域 Subclonal 結構偵測 — Intermediate AF × 甲基化多樣性雙重證據鏈 |
| 時間範圍 | 2026-04-14 ~ 2026-04-20 |
| 受眾 | PI（全面教學深度） |
| 總 slides | 39（35 content + 4 dividers） |
| 數據圖 | 12 張 |
| 教學概念圖 | 12 張（8 沿用 + 4 新增） |
| 預估報告時間 | ~55 min |
| 版本 | V2（從 V1 27 slides 擴充，新增 Part B 背景、Part C CN 發現、4 張新增教學、Paired 驗證、Normal ASM 願景、驗證需求） |

---

## Part A：封面與導覽（~1 min）

### Slide 1: 封面

> **講者筆記**: 開場直接點出本週核心結論 — POSITIVE。這是 InterSubMod 研究中第二個 POSITIVE 結論（第一個是 HPFineNGroups subclone marker），且是第一個建立完整六層證據鏈的結論。本週的核心發現是：在 LOH 區域中，intermediate allele frequency 的 variants 有顯著更高的甲基化多樣性（NGroups），這在 7 個獨立 cell line 中一致出現，構成 subclonal LOH 的遺傳+表觀遺傳雙重證據鏈。

---

### Slide 2: 本週重點與報告路線圖

> **講者筆記**: 先讓 PI 看到 4 個核心數字，建立預期：
>
> 1. **24.6% intermediate AF** — 在 purity=1.0 的 cell line 中，LOH 區域的 TP 有四分之一的 AF 不在期望值（0 或 1.0），這在完全 clonal 的 LOH 中不應出現
> 2. **ΔNGroups = +0.705** — 這些 intermediate AF variants 的甲基化群組數比 extreme AF 高 0.705，且 7 個樣本全部 p < 10⁻³⁹
> 3. **NR-bin 控制效應增強** — 控制 NumReads confound 後，效應量從 |r|=0.48 增加到 0.71，排除 confound 為主因
> 4. **Segment ρ=0.270** — 在 LOH segment 層級，AF 變異性與 NGroups 正相關，6/7 樣本方向一致
>
> **V2 更新**：路線圖擴充為「背景→LOH 驗證→CN 發現→教學→數據→整合→競爭者→未來」，標注報告預估 ~55 min。讓 PI 了解本次報告先花約 6 分鐘回顧為什麼研究 TO×LOH 的完整背景。

---

## Section Divider D1: 背景與 LOH 驗證

> 過渡頁，停頓 2 秒讓觀眾準備。

---

## Part B：為什麼研究 TO × LOH — 完整背景（~6 min）

### Slide 3: 從 Paired 到 TO — 為了統計檢定力

> **講者筆記（口述稿）**:
>
> 在進入本週的 subclone 分析之前，先讓我回顧一下為什麼我們在 Tumor-Only (TO) 模式下做研究。
>
> ISM 的原始設計是用 Paired 模式（有 normal BAM 對照），Paired 模式下 FP 非常少 — 在 HCC1395 中只有 3,429 個 FP，佔總 variant 的約 1%。TP:FP 比大約是 95:1。這個比例太懸殊了，我們沒有足夠的 FP 來做有統計力的比較分析。
>
> 而 TO 模式因為沒有 normal BAM 對照，germline variant 會被誤判為 somatic，FP 暴增到 128,382 個（約 31%）。TP:FP 比大約 2.3:1，這讓我們有足夠的 FP 數量來驗證 ISM 特徵是否能區分 TP 和 FP。
>
> 簡單說：**我們借用 TO 的高 FP rate 作為驗證工具**。TO 的 FP 主要是 germline variant，它們在分子層面有不同的特性（例如 AF 接近 1.0、甲基化 pattern 穩定），這讓我們可以研究 ISM 特徵在不同生物學背景下的表現。
>
> 出處：v5 LOH report Slide 2; Paired=3,429/328,699, TO=128,382/419,692

### Slide 4: LOH.bed 金標準驗證

> **講者筆記（口述稿）**:
>
> 接下來確認 LOH.bed 這個定義本身是否可信。我們用 FDA 的 SEQC2 金標準做了驗證。
>
> SEQC2 (Sequencing Quality Control Phase 2) 是一個多中心共識的 benchmark — 結合了 WGS、SNP array、多個 bioinformatics pipeline 的結果，定義了 HCC1395 的 LOH 真值區域。
>
> 我們的 LongPhase-TO 產出的 LOH.bed 與 SEQC2 truth 的比較結果：
> - **Jaccard index = 0.928** — 表示兩個 LOH 集合有 92.8% 的重疊
> - **Sensitivity = 0.961** — 96.1% 的真實 LOH 被偵測到
> - **Precision = 0.964** — 96.4% 的偵測結果是正確的
> - **F1 = 0.963** — 綜合指標接近完美
>
> 這個結果確認 LOH.bed 是高度可信的。右邊的染色體 ideogram 圖可以看到，LOH 的位置幾乎完全吻合 SEQC2 truth，只有少數 FP（紅色區塊）和 FN（橙色區塊）。
>
> 出處：v3 seqc2_validation section; J=0.928, Sens=0.961, Prec=0.964, F1=0.963

### Slide 5: LOH 雙定義 — HP Imbalance 包含 LOH.bed

> **講者筆記（口述稿）**:
>
> 在 ISM 中，我們有兩種方式定義「LOH 相關」的 variants：
>
> 1. **HP Imbalance**：在位點層級（per-variant），HP_Ratio < 0.1 或 > 0.9。HP_Ratio 是一個 variant 的 reads 在兩個 haplotype 間的分配比例。如果幾乎所有 reads 都在同一邊，代表這個位點附近可能有 LOH。
>
> 2. **LOH.bed**：在區域層級（genomic region），LongPhase-TO 輸出的連續 LOH 區塊。
>
> 兩者的關係用四象限表來看：
> - Q1（兩者皆有）= 26.7% — LOH.bed 內且 HP Imbalance
> - Q2（只有 HP）= 15.2% — HP Imbalance 但不在 LOH.bed 內
> - Q3（只有 LOH.bed）= 0.07% — 在 LOH.bed 內但 HP 不 imbalance → 只有 286 筆
> - Q4（都沒有）= 58.1%
>
> Sensitivity = 99.7%。也就是說，HP Imbalance 幾乎完全覆蓋 LOH.bed。只有極少數（0.07%）的情況是在 LOH.bed 區域內但 HP 沒有 imbalance。這兩個定義是層級互補的 — LOH.bed 定義更嚴格的連續區域，HP Imbalance 捕捉所有零星位點。
>
> 出處：v3 loh_dual_definition section; Q1=26.7%, Q2=15.2%, Q3=0.07%, Q4=58.1%

### Slide 6: 前期結論橋樑 — LOH 是分層工具，非過濾器

> **講者筆記（口述稿）**:
>
> 上面確認了 LOH.bed 可信且 HP Imbalance 覆蓋完整。但上週的結論是：**LOH 不能用來 filter TP/FP**。讓我快速回顧原因：
>
> 1. **Filter 10/10 全部失敗** — 我們嘗試了 10 種基於 LOH 的 filter 策略（zone filter, LOH-aware QS, 閾值組合等），全部 FP/TP removal ratio < 2.0×，不符合安全約束（TP loss ≤2%）
> 2. **LOH 內甲基化 AUC ≈ 0.50** — 在 LOH 區域內，所有 ISM 甲基化特徵對 TP/FP 的判別力為零
> 3. **轉向 characterization** — ISM 的價值不在 filter，在於 read-level epigenetic characterization
>
> 但有一個重要觀察沒有被否定：**LOH FP rate = 0.239 < Non-LOH FP rate = 0.338**。這表示 LOH 區域確實有 TP enrichment。雖然不能用來 filter，但作為分層工具，它揭示了一個問題：LOH 區域的 variants 是否全部 clonal？如果不是，為什麼？
>
> 另外，caller_af (AUC=0.654) 是唯一超越全部 ISM 特徵組合的單一特徵。它來自 variant caller (ClairS-TO)，與 ISM 無關。這啟發了我們從 AF 的角度重新審視 LOH 區域。
>
> 這就是本週的起點：**LOH 區域有 TP enrichment + AF 有異常值 → 可能有 subclonal LOH → 甲基化能否偵測？**
>
> 出處：v3 C14 結論; 0413 LOH×CNV 總結; LOH FP rate 0.239 vs 0.338

---

## Part C：LOH+CNV 發現 — TP 分離（~5 min）

### Slide 7: Coverage_Multiple 與 CN Proxy

> **講者筆記（口述稿）**:
>
> 要把 LOH 區域分得更細，我們需要知道每個區域的 copy number。ISM 沒有直接做 CNV calling，但我們有一個替代指標 — Coverage_Multiple。
>
> Coverage_Multiple 的計算很簡單：每個 region 的 coverage 除以全基因組的平均 coverage。如果一個區域的 coverage 是平均值的一半，Coverage_Multiple ≈ 0.5，對應 CN=1 (deletion)。如果接近 1.0，對應 CN=2 (diploid or cnLOH)。
>
> 我們用 HCC1395 的 SEQC2 truth CN 做了校準，Pearson r = 0.831，確認 Coverage_Multiple 是合理的 CN proxy。不是精確的（不做 GC correction、不考慮 ploidy），但作為分層工具足夠可靠。
>
> 基於此，定義了 4 個 CN Tier：
> - CN1 (CM < 0.75): 接近 1 copy → Deletion LOH
> - CN2 (0.75 - 1.25): 接近 2 copies → cnLOH 或 diploid
> - CN3 (1.25 - 1.75): 接近 3 copies → single-allele gain
> - CN4+ (> 1.75): 4 copies 以上 → high gain
>
> 閾值選擇理由：整數 CN boundary ±25%，足夠寬鬆以容納 coverage 波動。
>
> 出處：validated report §3.4; HCC1395 SEQC2 CN correlation r=0.831

### Slide 8: CN Tier 分層揭示 AF 分離

> **講者筆記（口述稿）**:
>
> 這張 slide 是本週研究的關鍵發現之一。把 LOH 的 TP 按 CN Tier 分層後，intermediate AF 的比例呈現戲劇性的差異：
>
> - **CN1: 16.9%** intermediate — 最低比例
> - **CN2: 24.8%**
> - **CN3: 45.2%**
> - **CN4+: 73.1%** — 最高比例
>
> CN4+ 的 73.1% 高比例很大程度是 **allele dosage effect**。當一個區域有 4+ copies，不同 allele 的拷貝數自然不均等（例如 3 copies of allele A + 1 copy of allele B），AF 自然不是 0 或 1.0。這是正常現象，不需要 subclone 來解釋。
>
> 但 **CN1 的 16.9% 沒有 allele dosage 解釋**。Deletion LOH 只留 1 個 copy，在 clonal LOH + purity=1.0 下，AF 只能是 0 或 1.0。那 16.9% 的 intermediate AF 從哪來？唯一的解釋是 subclonal LOH — 只有一部分細胞經歷了 deletion。
>
> 這就是為什麼 CN1 是「最乾淨的 subclone 信號」— 沒有 confound，物理限制只允許一個解釋。
>
> 出處：validated report §4.1.2; CN1=16.9%, CN2=24.8%, CN3=45.2%, CN4+=73.1%

### Slide 9: 生物學解釋 — 為什麼 CN1 可以分離

> **講者筆記（口述稿）**:
>
> 讓我用數學推導來說明 CN1 的特殊性。
>
> 在 Deletion LOH (CN=1) + purity=1.0 的情況下：
> - 正常細胞有兩條染色體 (A和B)，AF = 0.5
> - Clonal LOH 後只剩一條 (假設 A)，AF = 0 或 1.0
>
> 如果是 Subclonal LOH，fraction s 的細胞有 deletion (CN=1)，(1-s) 的細胞正常 (CN=2)：
>
> AF = s × 1.0 + (1-s) × 0.5 = 0.5 + 0.5s
>
> 當 0 < s < 1，AF 落在 0.5-1.0 的 intermediate 區間。
>
> 反過來，CN4+ 的情況：假設 3 copies of allele A + 1 copy of allele B，AF = 3/4 = 0.75。這是 allele dosage 自然產生的 intermediate AF，完全不需要 subclone。
>
> 所以 **CN1 的 intermediate AF 是物理約束下的唯一解釋 = subclonal LOH**。CN4+ 有替代解釋 (dosage)，CN1 沒有。這是本研究聚焦 CN1 的核心理由。
>
> 接下來我們需要教學幾個關鍵概念，然後用甲基化數據作為獨立的第二條證據鏈來驗證。
>
> 出處：validated report §2.1 推導; §4.1.2 CN1 分析

---

## Section Divider D2: 教學 — LOH 到 Subclonal LOH 完整推導

> 過渡頁。提醒觀眾接下來約 18-22 分鐘的教學，涵蓋從 LOH 基礎到 ISM 驗證機制的完整推導。

---

## Part D：教學（~18-22 min）

### Slide 10: 什麼是 LOH — 三種機制

> **講者筆記（口述稿）**:
>
> LOH，全名 Loss of Heterozygosity，是腫瘤中非常常見的基因組事件。正常細胞有兩套染色體（父源和母源），每個基因座有兩個等位基因。LOH 就是其中一個等位基因丟失，只剩下一個。
>
> LOH 有三種主要機制：
>
> 1. **Deletion LOH (CN=1)**：物理性刪除一條染色體的一段。細胞從 2n 變成 1n，定序的 coverage 下降到平均值的一半左右。這是最「乾淨」的 LOH — 一個 allele 完全不見了。
>
> 2. **cnLOH (Copy-Neutral LOH, CN=2)**：先刪除一條，然後把剩下的複製一份。最終 CN 還是 2，coverage 不變，但兩條染色體序列完全相同（都是同一個 allele）。
>
> 3. **Gain + LOH (CN>2)**：刪除加上額外的增益。Coverage 可能上升。
>
> ISM 使用 LongPhase-TO 的 `tumor_phased_LOH.bed` 偵測 LOH 區域。LOH.bed 的可信度已在 Slide 4 確認（SEQC2 J=0.928）。
>
> 出處：validated report §2.1 生物學原理

### Slide 11: Clonal vs Subclonal LOH

> **講者筆記（口述稿）**:
>
> 重要的一點：LOH 不一定發生在所有腫瘤細胞中。
>
> **Clonal LOH** 指的是 100% 的腫瘤細胞都經歷了同一個 LOH 事件。如果是 deletion (CN=1)，所有 reads 只來自保留的那個 allele，AF = 0 或 1.0。
>
> **Subclonal LOH** 指的是只有一部分比例 s 的腫瘤細胞有 LOH，其餘 (1-s) 保留兩個 allele。讀出來的 reads 是兩群細胞的混合，AF 落在中間某個值。
>
> 關鍵前提：我們用的是 cell line，purity = 1.0，沒有 normal cell 汙染。如果看到 intermediate AF，不可能是 normal cell 造成的稀釋 — 只能是 subclonal LOH。
>
> 但 cell line 不代表所有腫瘤細胞基因組都一樣。例如 HCC1954 已知是高度 subclonal 的乳癌 cell line，後面會看到它有 60.2% 的 CN1 TP 是 intermediate AF。而 COLO829 是相對 clonal 的黑色素瘤 cell line，只有 7.0%。
>
> 出處：validated report §2.1, §4.1.3

### Slide 12: AF 的物理意義與數學推導

> **講者筆記（口述稿）**:
>
> 現在推導為什麼 intermediate AF = subclonal LOH。
>
> 考慮 Deletion LOH (CN=1)，purity = 1.0：
>
> - **Clonal LOH** (s=1)：AF 只能是 0 或 1.0
> - **Subclonal LOH** (fraction s, 0<s<1)：
>
>   AF_obs = s × 1.0 + (1-s) × 0.5 = 0.5 + 0.5s
>
>   當 s 在 0 到 1 之間，AF 落在 0.5 到 1.0 — intermediate 區間。
>
> AF 三分類定義：
> - **Extreme**: AF < 0.1 or > 0.9 → clonal LOH 期望值
> - **Near-half**: 0.4-0.6 → heterozygous 或 cnLOH 複製後
> - **Intermediate**: 0.1-0.4 or 0.6-0.9 → 異常值，可能 subclonal
>
> 閾值選擇：0.1/0.9 與 TITAN BAF 閾值一致；0.4/0.6 涵蓋 heterozygous ±0.1 波動。
>
> 佐證數字：LOH TP 24.6% intermediate vs FP 4.1% → 6× enrichment。FP 的 AF 接近 1.0 是正常的（cell line 中 germline 幾乎都是 homozygous）。
>
> 出處：validated report §3.3, §4.1.1

### Slide 13: Self-Phasing Circular Dependency ← 新增

> **講者筆記（口述稿）**:
>
> 在繼續之前，需要回答一個重要的質疑：LongPhase-TO 的 self-phasing 問題是否影響本研究的 subclone 發現？
>
> Self-phasing 是 TO 模式特有的問題。LongPhase-TO 在做 haplotype phasing 時，會使用 somatic variant 作為 phasing anchor。這造成了一個循環依賴：
>
> 1. LongPhase-TO 用 somatic variant 做 phasing anchor
> 2. TP 的 ALT reads 天然偏向一個 haplotype（因為 somatic mutation 是單 allele 的）
> 3. HP_Ratio 被人為推向極端值（d = -1.20，大效應）
> 4. 62% 的 LOH 因此消失（PON-only 修正後）
>
> 這是一個已確認的 artifact，已經在 v5 和 v3 報告中詳細討論過。
>
> **關鍵問題**：這是否影響本研究的 AF 分析？
>
> **答案是否**。caller_af 來自 variant caller（ClairS-TO），是在 variant calling 階段計算的，完全獨立於 haplotype phasing。intermediate AF 是 caller 層面的觀察，不是 HP 層面的。Self-phasing 影響的是 HP_Ratio 和 LOH.bed 的精確邊界，不影響 caller 報告的 AF 值。
>
> 所以本研究的 subclone 發現（基於 caller_af）在邏輯上獨立於 self-phasing 問題。
>
> 出處：v5 Slide 11; v3 self-phasing 分析; d=-1.20; 62% LOH 消失

### Slide 14: Coverage_Multiple 校準 ← 新增

> **講者筆記（口述稿）**:
>
> Coverage_Multiple 是本研究的關鍵分層工具，需要確認它的可靠性。
>
> 計算方法很簡單：ISM 計算每個 region 的 coverage，除以 genome-wide average coverage，得到 Coverage_Multiple (CM)。CM ≈ 1.0 表示 diploid region，CM ≈ 0.5 表示 deletion (CN=1)，CM ≈ 1.5 表示 3 copies。
>
> 校準用的是 HCC1395 的 SEQC2 truth CN — 這是 FDA 金標準，結合了多個 pipeline 的 consensus CN call。散佈圖顯示 CM vs truth CN 的 Pearson r = 0.831，表示 CM 可以合理反映 CN。
>
> CN Tier 的閾值設在 0.75, 1.25, 1.75 — 接近整數 CN boundary ±25%。這個 margin 足夠寬鬆以容納 coverage 的隨機波動。
>
> 但要明確的是，CM 不是精確的 CN。它不做 GC content correction，不考慮 ploidy，也不做 segmentation。所以 CN Tier 是一個 proxy，不是 ground truth。這就是為什麼我們聚焦 CN1 — 它離 CN=2 (diploid) 最遠，被錯誤分類的風險最低。
>
> 出處：validated report §3.4; HCC1395 SEQC2 校準 r=0.831

### Slide 15: CAMDAC 原理 — 甲基化如何反映 LOH

> **講者筆記（口述稿）**:
>
> 這是本報告最重要的教學概念 — 為什麼甲基化可以獨立驗證 LOH 的 clonal/subclonal 狀態。
>
> CAMDAC（Cancer Allele-specific Methylation in Deconvolution After Correction）是 2020 年 TRACERx 計劃發表的方法。核心原則：
>
> - **Clonal LOH → ASM 消失**：所有細胞只有一個 allele，所有 reads 來自同一個等位基因。同一 allele 內甲基化 pattern 相對均一。結果：HPFineNGroups = 1，AlleleDelta ≈ 0。
>
> - **Subclonal LOH → ASM 部分保留**：只有一部分細胞有 LOH，其餘保留兩個 allele。Reads 是混合的。如果 allele A 和 allele B 有不同甲基化（ASM），混合後出現多群結構。結果：NGroups > 1，AlleleDelta > 0。
>
> ISM 的數據完美驗證了 CAMDAC 原理：
> - Extreme AF: AlleleDelta = 0.003 → 幾乎為零，ASM 消失
> - Intermediate AF: AlleleDelta = 0.031 → **增加 12 倍**
> - Cohen's d = +0.724（大效應）
>
> 這就是「雙重證據鏈」的核心：遺傳證據（AF 偏離期望值）和表觀遺傳證據（NGroups/AlleleDelta 升高）**同時指向同一個解釋**。
>
> 出處：validated report §2.1 CAMDAC 原理; §4.2.4 AlleleDelta 比較

### Slide 16: AlleleDelta / ASM 原理 ← 新增

> **講者筆記（口述稿）**:
>
> 進一步解釋 AlleleDelta 的物理意義。
>
> ASM (Allele-Specific Methylation) 是指同一個基因座的兩個 allele 有不同的甲基化水平。這由 imprinting、CTCF 結合差異、chromatin accessibility 等分子機制驅動。
>
> ISM 如何量測 ASM？每個 variant 的 reads 已經被 LongPhase-TO 分配到 HP1 或 HP2。AlleleDelta 計算兩個 haplotype 的平均甲基化率差異：
>
> AlleleDelta = |mean(HP1 methylation) - mean(HP2 methylation)|
>
> 在 clonal LOH 下，所有 reads 來自同一 allele（即使 HP1/HP2 標籤不同，它們的甲基化 pattern 應該一樣），所以 AlleleDelta ≈ 0。
>
> 在 subclonal LOH 下，部分 reads 來自只有 allele A 的細胞，部分來自有 A+B 的細胞。如果 A 和 B 甲基化不同（ASM），HP1 和 HP2 的平均甲基化率就會有差異 → AlleleDelta > 0。
>
> 在本研究中：
> - Extreme AF: AlleleDelta = 0.003（幾乎零）
> - Intermediate AF: AlleleDelta = 0.031（+12× 增加）
> - Cohen's d = 0.724（大效應）
>
> 這個 +12× 的劇變直接支持 CAMDAC 機制 — subclonal LOH 保留了 ASM，ISM 能偵測到。
>
> 出處：validated report §4.2.4; Extreme=0.003, Inter=0.031, d=0.724

### Slide 17: Read-Level 分析原理 — ISM 如何工作 ← 新增

> **講者筆記（口述稿）**:
>
> 在看數據之前，先理解 ISM 是如何工作的 — 它跟傳統的甲基化分析有根本性的不同。
>
> ISM 的 pipeline 流程：
> 1. 從 BAM 檔案提取 MM/ML 標籤（ONT basecaller 輸出的甲基化修飾資訊）
> 2. 透過 CIGAR alignment 進行 CpG 座標校正（精確定位每個 read 上的 CpG 位置）
> 3. 建構 Read × CpG 甲基化矩陣（每個元素是該 read 在該 CpG 的甲基化機率，0-255 scale）
> 4. 計算 reads 間的距離（Bernoulli, L1, NHD 等多種距離度量）
> 5. Agglomerative clustering 分群
> 6. 輸出 NGroups（最佳群組數，1-4）
>
> **ISM vs 傳統方法**的關鍵差異：
> - 傳統方法：對一個 region 算平均甲基化率 → 丟失了 read 間的差異資訊
> - ISM：保留每一條 read 的個別甲基化模式 → 能偵測 read 群之間的異質性
>
> 這就是為什麼 ISM 能偵測 subclone — 不同細胞群的 reads 有不同的甲基化 pattern，ISM 透過 clustering 把它們分開。
>
> 出處：ISM README; validated report §3.1

### Slide 18: HPFineNGroups — ISM 的甲基化分群指標

> **講者筆記（口述稿）**:
>
> HPFineNGroups 是 ISM 核心產出的指標。計算過程是：在每個 somatic variant 的 1kb 視窗內，建構 Read × CpG 甲基化矩陣，執行 agglomerative clustering，用 F-statistic 判斷最佳群組數（1-4）。
>
> NGroups 的生物學意義：
> - **NGroups=1**: 所有 reads 甲基化 pattern 一致 → 均質（clonal LOH 或單 allele）
> - **NGroups=2**: 有兩群不同 pattern → 可能是雙 allele ASM 或 subclone
> - **NGroups=3-4**: 更複雜結構 → 多個 subclone 或複雜的甲基化 heterogeneity
>
> Phase BCD 的先前驗證：NGroups ≥ 4 + NR ≥ 80 → TP rate 89.1%。NGroups 已被確認是 somatic heterogeneity marker。
>
> 在本研究的 CN1 LOH TP 中：
> - Extreme AF: 90.7% 是 NGroups=1（符合 clonal LOH 預期 — 單 allele，甲基化均一）
> - Intermediate AF: 79.6% 是 NGroups≥2（符合 subclonal LOH 預期 — 混合群，甲基化異質）
>
> 這個 90.7% vs 79.6% 的對比直接支持 CAMDAC 機制的預測。
>
> 出處：validated report §3.2 NGroups 定義; §4.2.2 NGroups 分布; Phase BCD NGroups≥4+NR≥80 TP rate

### Slide 19: NumReads Confound — 為什麼必須控制

> **講者筆記（口述稿）**:
>
> 在聲稱 NGroups 差異是「真信號」之前，必須排除一個重要的 confound — NumReads (NR)。
>
> 問題是：讀數少的 region（例如 NR=15），即使真的有 2 群甲基化 pattern，也可能因為統計力不足而偵測不到 → NGroups 被低估。反過來，讀數多的 region（NR=100），群組結構更容易被偵測 → NGroups 偏高。
>
> 如果 Intermediate AF 的 variants 碰巧有更多 reads，那 NGroups 差異可能只是 NR confound。
>
> 控制方法：NR-bin 分層。把 NR 分成 5 個 bins（10-30, 30-50, 50-80, 80-150, 150-500），在 NR 完全匹配的條件下重新比較。
>
> 判斷邏輯很簡單：
> - 如果 NR 是主因 → 控制後效應應該消失
> - 如果真信號 → 控制後效應不消失，甚至增強
>
> 事先確認：Intermediate 和 Extreme 的 NR 比值只有 1.01-1.34（差異很小），NR confound 本身不太可能解釋 ΔNG=+0.705 的大差異。控制後的結果在 Part E 展示。
>
> 出處：validated report §3.5 NR-bin 方法; NR ratio 1.01-1.34

### Slide 20: 證據鏈架構預告 — 六層金字塔

> **講者筆記（口述稿）**:
>
> 在進入數據之前，先預告我們的證據鏈架構。本研究建立了六層交叉驗證：
>
> - **L1 Genetic**: Intermediate AF 的存在本身（purity=1.0 下不應有 → subclonal）
> - **L2 Epigenetic**: Intermediate AF → NGroups 升高（甲基化獨立驗證）
> - **L3 Confound exclusion**: NR-bin 控制後效應不消失反增強
> - **L4 Mechanistic**: AlleleDelta +12× 支持 CAMDAC ASM 機制
> - **L5 Spatial**: LOH segment 內 AF-SD 與 NGroups 正相關
> - **L6 Technical**: HCC1395 ONT vs DORADO 技術重複一致
>
> 為什麼需要多層？因為單一證據可能有 alternative explanation。例如 NGroups 升高可能是 NR confound (L3 排除)；AF 偏離可能是 random noise (L5 的 segment 一致性排除)。多層交叉排除所有已知 confound 後，結論才夠強。
>
> 接下來的 Part E 逐層展示具體數據。
>
> 出處：validated report §5 六層證據鏈

---

## Section Divider D3: 核心結果 — 三步驟驗證

> 過渡頁。提醒觀眾：教學結束，接下來是數據。

---

## Part E：核心結果（~10-12 min）

### Slide 21: Step 1 — AF 分布 Baseline

> **講者筆記（口述稿）**:
>
> 第一步先看 AF 的基本分布。左邊是圖 01，LOH 的 TP 和 FP 的 AF 分布密度圖。可以明顯看到：
> - FP 的 AF 極度集中在接近 1.0（germline variant 在 cell line 中幾乎是 homozygous）
> - TP 的 AF 有明顯的 intermediate 尾巴，尤其在 0.5-0.9 之間
>
> 右邊圖 02 是各樣本的 intermediate AF 比例。Overall 24.6%，但樣本間差異很大：
> - COLO829: 7.0%（相對 clonal）
> - HCC1954: 60.2%（高度 subclonal）
> - 技術重複：HCC1395 ONT = 23.3%, DORADO = 20.0%（一致）
>
> 對比 FP 的 4.1% → 6× enrichment。這意味著 intermediate AF 在 TP 中是顯著高於 FP 的。
>
> **此圖證明**：Intermediate AF 在 LOH TP 中普遍存在，在 purity=1.0 下不應出現於 clonal LOH → 提示 subclonal LOH（L1 Genetic）。
>
> 出處：validated report §4.1.1; Overall=24.6%, FP=4.1%; 圖 01+02

### Slide 22: Step 1b — CN Tier 分層聚焦 CN1

> **講者筆記（口述稿）**:
>
> 第二步把 LOH TP 按 CN Tier 分層。這是 Slide 8 數據的詳細呈現。
>
> 左邊圖 05 的 4 panel 密度圖清楚顯示：CN1 的 AF 集中在 0 和 1.0 附近（clonal LOH 的期望），但有一個小但明確的 intermediate bump；CN4+ 則幾乎整個分布都在 intermediate 區間。
>
> 右邊圖 04 的 scatter 把 AF vs Coverage_Multiple 畫在一起，可以看到 CN1 區域（CM < 0.75）的 intermediate AF 點清楚可辨。
>
> 關鍵數字已在 Slide 8 討論：CN1=16.9%, CN4+=73.1%。CN1 是最乾淨的 subclone 信號。
>
> **此圖證明**：CN1 的 intermediate AF 無 allele dosage 解釋 → 最可靠的 subclone 信號 → 後續聚焦 CN1 分析（L1 refined）。
>
> 出處：validated report §4.1.2; 圖 04+05

### Slide 23: Step 2 — NGroups × AF 交叉分析（核心結果）

> **講者筆記（口述稿）**:
>
> **這是全報告最重要的一張 slide。**
>
> 左邊圖 06 是 NGroups 按 AF class 和 CN Tier 的交叉分析。聚焦 CN1：
> - Extreme AF: 平均 NGroups = 1.091（幾乎都是 1，符合 clonal LOH 預期）
> - Intermediate AF: 平均 NGroups = 1.796
> - **ΔNGroups = +0.705** — 這是巨大的差異
>
> 更直觀的對比：
> - Extreme AF 中 90.7% 是 NGroups=1（甲基化均一 → clonal）
> - Intermediate AF 中 79.6% 是 NGroups≥2（甲基化異質 → subclonal）
>
> 右邊圖 07 是 7 個樣本各自的結果。**全部 7 個樣本都顯著**：
> - 最低 p 值：7.8×10⁻⁴⁰（H2009，樣本最少的）
> - 最高效應量：|r| = 0.822（HCC1395_DORADO）
>
> NR ratio 只有 1.01-1.34，不足以解釋 +0.705 的差異（NR 差 30% 最多造成 ~0.1 的 NGroups 差異）。
>
> **此圖證明**：Intermediate AF variants 有顯著更高的甲基化多樣性 → H1 (AF 分離) + H4 (NGroups 差異) 都 supported（L2 Epigenetic）。
>
> 出處：validated report §4.2.1-§4.2.2; ΔNG=+0.705; 7/7 p<10⁻³⁹; 圖 06+07

### Slide 24: Step 2b — NR Confound 排除

> **講者筆記（口述稿）**:
>
> 控制 NR confound 後的結果。
>
> 左邊圖 08 是 NR-bin 分層的 NGroups 分布：
> - NR 10-30: ΔNG = +0.484, |r| = 0.483
> - NR 30-50: ΔNG = +0.715, |r| = 0.709
> - NR 50-80: ΔNG = +0.711, |r| = 0.708
>
> 關鍵觀察：**效應隨 NR 增加而增強，而非消失**。這正好和 confound 假說的預測相反（如果是 confound，控制後效應應該消失）。低 NR bin 的效應較弱是 floor effect — reads 太少偵測不到群組結構，不是沒有群組。
>
> 右邊圖 09 是 6 個甲基化特徵的比較：
> - **AlleleDelta**: d = +0.724（大效應，ASM +12×）— 直接支持 CAMDAC 機制
> - **HPFineF**: d = +0.639（大效應，F-statistic +23×）
> - **CramersV**: d = +0.318（中效應）
> - **PairwiseMeanDist**: d = +0.031（無差異）— 這是 negative control。如果所有甲基化特徵都升高，可能是某種系統性 bias；但 PairwiseMeanDist 不變，說明效應是特異性的。
>
> **此圖證明**：NGroups 差異獨立於 NumReads（L3 Confound exclusion）；AlleleDelta 的 +12× 增加直接支持 CAMDAC ASM 機制（L4 Mechanistic）。
>
> 出處：validated report §4.2.3 NR-bin; §4.2.4 甲基化特徵; r 0.48→0.71; AlleleDelta d=0.724; 圖 08+09

### Slide 25: Step 3 — Segment 空間一致性

> **講者筆記（口述稿）**:
>
> 最後一步：如果 intermediate AF 真的反映 subclonal LOH 事件，那同一個 LOH segment 內應該有一致的 AF pattern。AF 變異性高的 segment（混合了 clonal 和 subclonal variant）應該有更高的 NGroups。
>
> 左邊圖 11 是 AF-SD vs mean NGroups 的 segment-level scatter。CN1 overall：
> - Spearman ρ = 0.270, p = 5.6×10⁻²² — 高度顯著的正相關
> - 6/7 樣本方向一致
> - 最強：H1437 ρ=0.809, COLO829 ρ=0.763
> - 唯一反向：HCC1954 ρ=-0.297（但只有 n=34 segments，underpowered + 可能 subclone saturation）
>
> 右邊圖 12 是 Uniform vs Mixed segment 的比較：
> - Uniform segments（AF 一致）: 平均 NGroups = 1.292
> - Mixed segments（AF 變異大）: 平均 NGroups = 1.717 → +0.425
>
> **此圖證明**：AF 變異性高的 LOH segment 有更高的甲基化多樣性 → 是 segmental 事件而非 random noise（L5 Spatial）。
>
> 出處：validated report §4.3; ρ=0.270, p=5.6e-22; 6/7 positive; 圖 11+12

### Slide 26: Per-Sample 一致性總覽

> **講者筆記（口述稿）**:
>
> 彙整所有 7 個樣本的結果。這張表格呈現每個樣本的 ΔNG、Mann-Whitney p、rank-biserial |r|、和 segment ρ。
>
> 全部 7 個樣本的 ΔNG 都是正值（+0.329 到 +0.823），全部 p 值都極顯著（最高 7.8×10⁻⁴⁰ 在 H2009）。技術重複 HCC1395 ONT 和 DORADO 的結果一致（ΔNG 0.581 vs 0.823，方向相同）。
>
> Segment 一致性 6/7 positive。唯一反向的 HCC1954 有明確的解釋：只有 34 個 segments（underpowered），且已知是極度 subclonal（saturation effect — 幾乎所有 segment 都是 mixed）。
>
> **此圖證明**：這是跨 7 cancer cell lines × 2 技術平台的 universal effect（L6 Technical replication）。
>
> 出處：validated report §4.1.3, §4.2.5, §4.3.2; 圖 10+13

---

## Section Divider D4: 總結與展望

> 過渡頁。

---

## Part F：整合（~5 min）

### Slide 27: 六層證據鏈匯總

> **講者筆記（口述稿）**:
>
> 現在把六層證據鏈匯總。與 Slide 20 的空金字塔相呼應，這次每層都填入具體數據：
>
> | 層級 | 證據 | 核心數字 |
> |------|------|---------|
> | L1 Genetic | 16.9-60.2% intermediate AF in CN1 | purity=1.0 排除 dilution |
> | L2 Epigenetic | ΔNGroups = +0.705 | 7/7 p<10⁻³⁹ |
> | L3 Confound | NR-bin 控制後效應增強 | r: 0.48→0.71 |
> | L4 Mechanistic | AlleleDelta +12× (d=0.724) | CAMDAC ASM 機制 |
> | L5 Spatial | AF-SD ∝ NGroups (ρ=0.270) | 6/7 positive |
> | L6 Technical | HCC1395 ONT vs DORADO | 23.3% vs 20.0% |
>
> 結論：**POSITIVE** — Intermediate AF 在 LOH 區域對應 subclonal LOH 事件，且 ISM 的甲基化指標（NGroups, AlleleDelta）能獨立偵測這些事件。這是 ISM 作為 read-level epigenetic characterization 工具的直接佐證。
>
> 出處：validated report §5 綜合結論

### Slide 28: Paired 模式驗證可行性分析 ← 新增

> **講者筆記（口述稿）**:
>
> PI 可能會問：「如果在 Paired 模式下重做，結果會一樣嗎？」
>
> 三個層面的分析：
>
> 1. **caller_af 不受 self-phasing 影響** ✅
>    caller_af 由 ClairS-TO 的 variant calling 模組計算，完全獨立於 haplotype phasing。即使 HP_Ratio 被 self-phasing bias 推偏，caller_af 不受影響。所以我們在 TO 下觀察到的 intermediate AF 和 NGroups 差異，在邏輯上是 self-phasing independent 的。
>
> 2. **Paired LOH.bed 更可靠** ✅
>    Paired 模式有 normal BAM 作為 phasing reference，不存在 self-phasing 問題。用 Paired 的 LOH.bed 可以獨立確認 LOH 區域定義的準確性。
>
> 3. **但 Paired FP 太少** ⚠
>    Paired 模式只有 ~3,429 FP（~1%），做大規模統計分析的 power 不足。Subclone 分析需要足夠多的 LOH 內 variant，而 TO 的 ~128,382 FP 提供了這個統計力。
>
> 結論：TO 的 subclone 發現在邏輯上獨立於 self-phasing 問題。Paired 可以做 LOH.bed 定義的獨立驗證，但不適合做大規模的 AF × NGroups 統計分析。
>
> 出處：v5 Paired/TO 比較; self-phasing 獨立性分析

### Slide 29: Normal ASM + Tumor ASM 共分析願景 ← 新增

> **講者筆記（口述稿）**:
>
> 當前研究有一個未解決的 confound：NGroups 升高可能是 germline ASM 而非 subclone 造成的。也就是說，有些位點天生就有兩個 allele 不同的甲基化（這是正常的 germline ASM），不需要 subclone 就會有 NGroups > 1。
>
> 解決方案是引入 Normal BAM 作為 reference。概念很簡單：
> - 如果一個位點在 Normal 和 Tumor 都有高 NGroups → 這是 Germline ASM → 排除
> - 如果一個位點只在 Tumor 有高 NGroups → 這是 Tumor-specific ASM → subclone 信號
>
> Phase 2 方向 A+D 正在做這件事：
> - Normal BAM integration 已開始
> - LOH BED annotation 已完成
> - 下一步是 per-variant Normal NGroups 計算
>
> Normal reference 是提升 subclone NGroups 可信度的關鍵下一步。它能把「可能是 germline ASM」這個 confound 從 L4 提升到真正的機制驗證。
>
> 出處：Phase 2 A+D 設計文件; germline ASM confound 分析

---

## Part G：TumorLens 競爭者比較（~3 min）

### Slide 30: 功能比較表

> **講者筆記（口述稿）**:
>
> TumorLens 是 2026 年發表在 medRxiv 的 preprint，是第一個 long-read unified tumor analysis pipeline。它同時偵測 SNV, indel, SV, CNV, LOH, 和 CpG methylation。
>
> 功能比較表顯示兩個工具的互補性：
> - TumorLens 是全基因組多模態偵測 → **macro-level**。它做 variant calling, CNV, LOH, 和 region-level 甲基化分析。
> - ISM 是 per-variant read-level 甲基化 clustering → **micro-level**。它不做 variant calling，但給每個 variant 附加精細的甲基化結構 context。
>
> ISM 的獨特優勢：**per-variant read-level methylation clustering**。沒有其他工具在做這件事。TumorLens 的甲基化分析是 region-level 的平均值，ISM 保留每條 read 的個別模式。
>
> 出處：TumorLens preprint (medRxiv 2026), doi: 10.64898/2026.03.18.26348569

### Slide 31: 定位圖 — Macro vs Micro 互補

> **講者筆記（口述稿）**:
>
> 這張定位圖把 ISM 和 TumorLens 放在同一個框架中。TumorLens 做 macro（全基因組總覽），ISM 做 micro（個別 variant 的表觀遺傳 context）。兩者是互補的，不是競爭的。
>
> 2026 年的學界趨勢是 LOH + CNV + methylation 的 unified approach。TumorLens 和 ISM 都在朝這個方向走，但切入角度不同：
> - TumorLens：從 genome-wide 偵測出發，整合多個分析模組
> - ISM：從 per-variant read-level 出發，提供其他工具無法產出的精細 context
>
> ISM 的差異化定位：**唯一做 read-level per-variant epigenetic context 的工具**。本週的 subclone 偵測結果正是這個定位的具體佐證。
>
> 出處：TumorLens preprint; ISM 定位分析

---

## Part H：未來方向（~4 min）

### Slide 32: ISM 工具定位

> **講者筆記（口述稿）**:
>
> 總結 ISM 的工具定位：Read-level Subclone Characterization。
>
> 它不是 variant filter — Phase 1A 已經確認 F1 增益有限（paired-pure ΔF1 = +0.0112，TO 模式甲基化增益為負）。
>
> 它是 epigenetic characterization 工具：
> - 給每個 somatic variant 附加甲基化結構 context（NGroups, AlleleDelta, CramersV 等）
> - 偵測 subclonal LOH/heterogeneity
> - 輸出 structural annotations 供下游工具使用
>
> 本週的結果直接佐證了這個定位：ISM 的 NGroups 能區分 clonal vs subclonal LOH（ΔNG = +0.705, 7/7 p<10⁻³⁹）。這正是 ISM 作為 characterization 工具的核心能力。
>
> 出處：ISM 工具定位; Phase 1A 結論; 本週 ΔNG=+0.705

### Slide 33: 區塊分割策略 + Timeline

> **講者筆記（口述稿）**:
>
> 基於本週的發現，規劃了短期和中期計劃：
>
> **短期計劃**：
> 1. 在 ISM output 加入 LOH segment annotation（LOH.bed + CN tier tag）
> 2. Per-segment summary statistics（AF mean/SD, NGroups mean, AlleleDelta mean）
> 3. Segment-level subclone classification：Uniform（clonal）vs Mixed（subclonal）
>
> **中期計劃**：
> 1. 跨 segment 的 subclone fraction estimation（利用 AF = 0.5 + 0.5s 反推 s）
> 2. Normal BAM reference 整合（Phase 2 A+D）→ 區分 germline ASM vs tumor ASM
> 3. 連結甲基化 block structure：同一 subclone 的 segments 甲基化 pattern 是否一致
>
> Timeline：Phase 1A (完成) → Phase 2 A+D (進行中) → Subclone annotation (下一步) → Block segmentation (遠程)
>
> 出處：研究路線圖; Phase 2 設計

### Slide 34: 驗證需求與風險 ← 新增

> **講者筆記（口述稿）**:
>
> 最後列出驗證需求的優先序和已知風險。
>
> **驗證優先序**：
> 1. Phase 2 A+D Normal BAM integration（進行中）— 解決 germline ASM confound
> 2. 臨床低純度樣本擴展驗證（待規劃）— cell line purity=1.0 的結論是否在 purity<1.0 時成立
> 3. Haplotag 重跑（blocker: 計算資源）— 用不同 phasing 策略驗證 LOH.bed 定義
>
> **已知風險/限制**：
> 1. Coverage_Multiple 不是精確 CN（無 GC correction、不考慮 ploidy）
> 2. Cell line ≠ clinical sample（purity=1.0 簡化了分析，臨床樣本更複雜）
> 3. LOH.bed 不區分 LOH 類型（Deletion vs cnLOH vs Gain — 需要 CN tier 輔助判斷）
> 4. HCC1954 的 segment 分析反向（n=34 underpowered + subclone saturation）
> 5. NGroups 上限 = 4（ISM clustering 目前只分 1-4 群，可能低估高 subclonality）
>
> 出處：研究限制清單; Phase 2 計劃

---

## Part I：結尾

### Slide 35: 討論與 Q&A

> **講者筆記**: 回顧核心結論 POSITIVE，列出 5 項限制，預告下週工作（Phase 2 A+D 推進、Subclone annotation 初步整合、LOH segment summary 功能），開放 Q&A。
>
> 重點答辯準備：
> - 如果被問「self-phasing 影響」→ 指向 Slide 13，caller_af 獨立於 phasing
> - 如果被問「CN1 以外呢」→ CN2/3/4+ 有 allele dosage confound，需進一步方法排除
> - 如果被問「HCC1954 反向」→ n=34 underpowered，且 HCC1954 已知極度 subclonal（可能所有 segment 都是 mixed）
> - 如果被問「NGroups 上限 4」→ 是 ISM clustering 的設計選擇，可考慮未來放寬
> - 如果被問「臨床價值」→ 需低純度樣本驗證，但機制上 subclone 在低純度下應更明顯（normal cell 混合增加 heterogeneity）
