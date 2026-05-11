# Self-Phasing 完整教授報告 · Storyboard（v1, 2026-04-29）

## 0. 報告定位

| 項目 | 設定 |
|------|------|
| 受眾 | 教授（具 NGS / long-read / variant calling 基本背景，但不熟 LongPhase-TO 與 InterSubMod 的職責切分） |
| 風格 | 演講者模式：slide 文字最少（標題 + 視覺為主），細節放 speaker notes |
| 目標 | 讓教授從基本定義 → 問題 → 改動敘述 → 差異與結論數據都能聽懂 |
| 結構 | 工程報告（problem → root cause → fix → verify）+ 研究報告（claim → evidence → caveat）混合 |
| 投影片數 | 30–34 slides（含封面、章節分隔、結論） |
| 比例 | 16:9 |
| 字型 | Latin: Arial / Helvetica；CJK: Droid Sans Fallback；圖內中文：matplotlib + DroidSansFallbackFull.ttf |
| 圖片規則 | 等比 fit-within，禁止強制 width+height 擠壓 |

---

## 1. PPT 大架構（7 段、34 slides）

```
Section 0 · 開場（1 slide · S1）
Section 1 · 基本定義（5 slides · S2-S6）          ← 認知同步、減少負荷
Section 2 · 問題敘述（4 slides · S7-S10）
Section 3 · 量化證據（4 slides · S11-S14）
Section 4 · 修補方向（5 slides · S15-S19）
Section 5 · 改動細節（5 slides · S20-S24）
Section 6 · 驗證與差異（5 slides · S25-S29）
Section 7 · 結論與下一步（5 slides · S30-S34）
```

---

## 2. 逐 slide storyboard

### Section 0 · 開場（S1）

| # | title | 視覺 | 文字（slide 上）| speaker note 重點 |
|---|------|-----|-----|-----|
| S1 | Self-Phasing：原因、修補、驗證 | Pipeline 視覺（caller→phasing→tagging→ISM 上下游） | 標題 + 副標 + 日期 | 開場：今天從基本定義開始，介紹 self-phasing 是什麼、它如何在 TO 模式爆發、我們怎麼修、修完之後驗證了哪些東西、哪些還沒解決 |

### Section 1 · 基本定義（S2-S6）

| # | title | 視覺 | speaker note |
|---|------|-----|-----|
| S2 | 三條 pipeline 上下游關係 | `fig1_pipeline_comparison.png` | ClairS-TO 是 caller，輸出 VCF；LongPhase-TO 是 phasing/haplotag 工具，輸出 phased VCF + tumor_tagged.bam + LOH.bed；InterSubMod 是下游消費者，讀 BAM 的 HP tag 算 region-level 特徵 |
| S3 | Phasing 是什麼？ | 自繪 phasing graph schematic（在 build_pptx 內畫） | 把 reads 上的 variants 串成同一條 haplotype，每個 phase block 內部 GT 用 `\|` 表示，PS 是 phase-set ID |
| S4 | 三層資料的角色（caller / phasing / haplotag） | 三層彩色框圖 | Caller 看 VCF FILTER/AF/GQ；Phasing 看 phased VCF GT/PS + LOH.bed；Haplotag 看 BAM HP:i tag。**強調：三者不可混用**，這也是過去把 ISM HP_Ratio LOH 與 LOH.bed 混為一談的根因 |
| S5 | HP tag 的 5 個整數值 | `fig17_hp_tag_5versions.png` | 1 = germline HP1；2 = germline HP2；11 = somatic on HP1；21 = somatic on HP2；33 = somatic ambiguous（兩 hap 都不確定）；0 = unphased |
| S6 | InterSubMod 的下游角色 | 簡化 pipeline 圖 + ISM 強依賴標記 | ISM 強依賴 HP tag，因此上游任何 bias 直接傳遞下來；ISM 自己**不修補**，只是受影響 |

### Section 2 · 問題敘述（S7-S10）

| # | title | 視覺 | speaker note |
|---|------|-----|-----|
| S7 | Self-phasing 是什麼？ | `fig2_self_phasing_concept.png` | 同一 sub-clone 的 somatic variants 共享 read population → long-read 跨多個 somatic 都帶 ALT → phasing graph 看到強連結 → 全部塞同一 phase block；本應隨機 50:50 變成偏 |
| S8 | AF=0.3 走例 | `fig3_af03_walkthrough.png` | 這張是教授必看的「為什麼 TO 會出現 self-phasing」核心圖：paired 模式 HP_Ratio≈0.5（有 normal 對照），TO 模式 HP_Ratio→0.94（無 normal，phasing 反客為主） |
| S9 | TO 模式為什麼特別嚴重？ | 文字框：兩條件對比 | 1) TO 無 matched normal；2) PON 不能完全排除 germline-like het；3) 預設 `--pon-only-phasing=false`（V2b 之前）讓 somatic 進 phasing graph |
| S10 | 三層 bug 同時暴露 | 三色框圖（Phase / Tag-priority / Tag-enum） | 不是單一 bug，是三層獨立故障在 PON-only 啟用後集中暴露：(A) Phase scaffold 用 somatic anchor、(B) getVote priority bug、(C) enum vs integer literal mismatch |

### Section 3 · 量化證據（S11-S14）

| # | title | 視覺 | speaker note |
|---|------|-----|-----|
| S11 | 全基因組層證據 17.3:1 | `fig4_evidence_summary.png` 或 fig01d Panel D | HCC1395 全基因組 baseline：HP1=614,000 / HP2=35,500（94.6% 集中於 HP1）vs 隨機預期 ~1:1；跨 23 染色體一致 |
| S12 | 個別位點層證據 | `D_SP1_chr19_17565944.png` | 個別位點可達 113:0 完全失衡；baseline 與 paired 方向相反 → V5 修正後與 paired 一致 |
| S13 | 跨樣本一致性 7/7 | 跨樣本 bar chart（build_pptx 自繪）+ Cohen's d=-1.20 | 7/7 樣本同方向；同位點 HP_Ratio 跨模式 r=0.001（無訊號）；TO-only LOH 在 paired 下 86.5% 是平衡的 |
| S14 | LOH 的兩個層次（精確化）| 兩欄表（ISM HP_Ratio LOH vs LOH.bed region-level） | 重點區分：ISM HP_Ratio LOH 受 self-phasing 影響 62% 是 artifact；LOH.bed Jaccard=1.0 完全不變 |

### Section 4 · 修補方向（S15-S19）

| # | title | 視覺 | speaker note |
|---|------|-----|-----|
| S15 | 修補位置在哪裡？ | 上下游圖 + 標紅 longphase-to-mod | 修補的對象**不是 InterSubMod**（本 repo 無 C++ 改動），而是另一個獨立的工具 `longphase-to-mod` |
| S16 | 4-commit 漸進演進 | `fig01a_commit_evolution.png` | V2b（8b8c1fd）→ V3-Fixed（41ff147）→ INDEL guard（380e8d2）→ V5（working tree） |
| S17 | 4 commit 對應的 4 條 bug | 4-row 表格 | V2b 解 phase scaffold；V3F 解 priority bug + enum mismatch；INDEL guard 補 UB；V5 解 V3F 過於保守 |
| S18 | V5 三層投票邏輯 | `fig01b_three_layer_logic.png` 或 `fig16_v5_threelayer_logic.png` | Layer 1：germline first；Layer 1.5：somatic fallback（V5 新增）；Layer 2：encode hpResult（11/21/33） |
| S19 | 候選方案排除 | 兩欄表（替代方案 vs 為何不採） | 替換 LongPhase 為 WhatsHap/HapCUT2：風險高、外部依賴；ISM 自己加 haplotag 邏輯：介面割裂、ISM F1=0.0124 對 TO germline FP 無修復力 |

### Section 5 · 改動細節（S20-S24）

| # | title | 視覺 | speaker note |
|---|------|-----|-----|
| S20 | 介面契約零變動 | 程式碼 diff 簡化視覺 | HaplotagProcess.h:66-68 三個 method signature 從 baseline 到 V5 一字未變；總修改量 +68/-36 行集中於 3 函式 |
| S21 | getVote 兩層 → 三層 | `fig01b_three_layer_logic.png` | V3F：兩層 germline-first → V5：插入 Layer 1.5 「else if (somaticHP1>0 \|\| somaticHP2>0)」純分支 |
| S22 | HP tag 對應表（V5 final） | 表格 (germlineResult, somaticTotal>0) → hpResult | (0,F)→0；(1,F)→1；(2,F)→2；(1,T)→11；(2,T)→21；(0,T)→33 |
| S23 | InterSubMod 的位置 | InterSubMod 與 longphase-to-mod 切分圖 | InterSubMod `src/core/` 不包含 HaplotagProcess；ISM 只是讀新版 BAM 後下游分析自動受惠 |
| S24 | V5 working tree 未 commit caveat | 警告框 | V5 = 380e8d2 + 兩塊 working-tree 修改；建議切 2 獨立 commits（後續 F1）|

### Section 6 · 驗證與差異（S25-S29）

| # | title | 視覺 | speaker note |
|---|------|-----|-----|
| S25 | 4 項硬性 sanity check | 4-row 表格 + PASS/FAIL 標記 | (1) Δ-consistency 守恆律 PASS；(2) Germline 不變 PASS；(3) Layer 1.5 期望 1 PASS；(4) 無 germline→HP33 PASS；15/15 PASS、0 violation |
| S26 | Baseline / V3F / V5 量化對比 | `fig18_concordance_amb_f1.png` | AMB% 17.5→8.0%、HP:i:33 reads −54%、clean PS paired GT concordance V5 90.5% vs Baseline 82.2%（+8.3 pp） |
| S27 | IGV 4-BAM 並列：V5max1 戲劇化 | `C_V5max1_chr19_4639528.png` | chr19:4639528：V3-F panel 紫色 HP33 群（39 reads） → V5 全部正確重分配為 HP11 → 守恆律 PASS |
| S28 | IGV 4-BAM 並列：SP1 orientation 翻轉 | `D_SP1_chr19_17565944.png` | chr19:17565944：baseline HP1 主導 → V2b/V3F/V5 全部翻轉為 HP2 主導；paired tumor 確認 HP2 為真實方向；3/3 SP-extreme 一致翻轉 |
| S29 | F1 為什麼幾乎不變？（重要釐清） | 流程圖：raw F1=0.7166 三版本相同 | ClairS-TO raw F1=0.7166 對所有版本完全相同（V5 不改 caller）；F1 變動來自 ISM SuggestFilter；V5 vs Baseline F1 = -0.0003（噪音級）→ **F1 不能衡量 V5 真實價值**，真實價值在 read-level tag 品質 |

### Section 7 · 結論與下一步（S30-S34）

| # | title | 視覺 | speaker note |
|---|------|-----|-----|
| S30 | ISM 影響 3-tier 分類 | `03_self_phasing_impact.png` 或 `fig5_impact_matrix.png` | 🔴 嚴重 38%（HP-依賴特徵，必重跑）／🟡 中度 7%（QualityScore 等）／🟢 不受影響 55%（PairwiseDist、Allele、Caller、甲基化矩陣）|
| S31 | 修正了哪些舊說法 | 兩欄表（舊說法 → 修正後）| HPFineNGroups 不是 methylation bimodality（是 phasing bucket）；HP:i:21 不必然是當前位點 ALT；Self-phasing 不污染 LOH.bed；V5 修補不等於 caller F1 改善 |
| S32 | 已知限制 R1-R9 | 9-row 表 | cnLOH 區、AF>0.9 邊界、15 sites cherry-picked、Confidence threshold 未驗證、V5 working tree 未 commit、7 樣本擴展未做、L4 部分指標 V5 略遜、problem PS 不適用 read-level、paired ground truth 自身用 LongPhase |
| S33 | 後續工作 F1-F8 | 表格 | F1 commit V5 working-tree（高）；F2 vote log 直接驗證；F3 7 樣本擴展；F4 master×flag 重跑 HPFineNGroups；F5 manifest haplotag_version；F6 cnLOH 獨立評估；F7 trio-phased 第二 ground truth；F8 100 隨機位點 cross-validate |
| S34 | TL;DR 一句話結論 | 標語式大字 | Self-phasing 在 LongPhase-TO 階段以 4-commit 漸進修補；InterSubMod 是下游消費者；V5 真實價值在 read-level tag 品質（+8.3 pp clean PS concordance），SEQC2 F1 不變是預期行為 |

---

## 3. 截圖驗證檢查項

每張 slide 必須通過：

1. **文字無溢出**（Latin + CJK 雙語不互相覆蓋）
2. **圖片等比 fit-within**（無強制縱橫比破壞）
3. **CJK 字型正確**（無方塊或亂碼）
4. **speaker note 存在且非空**
5. **圖片連結有效**（無 broken image placeholder）

驗證方式：`scripts/screenshot_all.py` 必須輸出 `[ISSUES] none detected`。

---

## 4. PPT 主軸設計原則

| 原則 | 落實方式 |
|------|---------|
| 投影片文字越少越好 | slide 標題 ≤ 15 字、bullet ≤ 5 條、每條 ≤ 12 字；長文移至 speaker notes |
| 認知同步 / 減少負荷 | Section 1 用 5 slides 補完所有定義（HP tag、phasing、PS、LOH 兩層）後才進入問題 |
| 工程報告骨架 | Problem (Section 2-3) → Root cause (Section 4) → Fix (Section 5) → Verify (Section 6) → Limit/Next (Section 7) |
| 研究報告骨架 | Claim（修補有效）→ Evidence（5 條獨立路徑）→ Caveat（R1-R9）三段 |
| 清楚定義 | 每個首次出現的術語在 speaker note 內必有 1-2 句解釋 |
