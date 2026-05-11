# 重點數據速查表（v2 PPT 配套）

> PPT 內所有提到的數字 + 來源行號交叉驗證表。Slide 編號對應 `00_storyboard_v2.md`。

## 0. ★ Purity threshold 0.95 觸發鎖（v3 新增 · S11 核心）

| 指標 | 數值 | 含義 | 來源 |
|------|------|------|------|
| 純樣本標準閾值 | **0.95** | `PhasingProcess.cpp:197` 硬編碼 `purity > 0.95` | 04_purity_calculator §X.2 |
| Baseline 估 HCC1395 purity | **0.927**（四捨五入 0.93）| **正確估計**，非誤判；但 ≤ 0.95 → **未觸發 Two-Pass** | 04 §1 |
| Baseline 估 0.6 sample purity | 0.607 | 與真實值 0.6 接近 | 04 §1 |
| V5 PON-only 修復前 purity（兩 sample）| **0** | ploidyRatioMap=nullptr → q1=q3=0 → 強制 0 | 04 §1, §3.3 |
| V5 PON-only 修復後 purity（0.93 sample）| 0.871 | vs baseline 0.927 誤差 -0.06 | 04 §X.2 |
| V5 PON-only 修復後 purity（0.6 sample）| 0.634 | vs baseline 0.607 誤差 +0.03 | 04 §X.2 |
| 修復函式 | `collectPloidyRatio()` `PhasingGraph.cpp:1147-1175` | 於 syncPhasingResultOrigins 後呼叫 | 04 |
| 修復是否需重跑 tag | **否** | 修復 side-effect free，既有 BAM 完全有效 | 04 §X.2 |

**因果鏈（v3 PPT S11 新增主者）**：純樣本 → baseline 估 0.927 ≤ 0.95 → **未觸發 Two-Pass** → 走「三條路」baseline 標準流程 → 流程假設有 normal contamination 但實為純 → 暴露 tag 層 somatic-first 投票 bias → self-phasing 17.3:1。

---

## 1. 觀察到的問題（baseline 量化）

| 指標 | 數值 | 樣本 | 來源 .md（行號）|
|------|------|-----|--------------|
| HP1 family : HP2 family ratio | **17.3 : 1** | HCC1395 全基因組 | source 03 §4.2、§5.A.1 第 282 行；source 02 §3.1 |
| HP1 reads（baseline）| **614,000** | HCC1395 全基因組 | source 03 §5.A.1 第 282 行 |
| HP2 reads（baseline）| **35,500** | HCC1395 全基因組 | source 03 §5.A.1 第 282 行 |
| 94.6% 集中於 HP1 | **94.6%** | HCC1395 跨 23 染色體一致 | source 03 §4.2.1 第 161 行 |
| HP1:HP2 randomized prediction | **~1 : 1** | 生物學期望 | source 03 §5.A，理論層 |
| Paired 模式 HP_Ratio | **~0.5** | 同樣本 paired 對照 | source 03 §4.2.1 |
| TO 模式 HP_Ratio | **0.94** | 同樣本 TO | source 03 §4.2.1 |
| 同位點跨模式 r | **r = 0.001** | n = 288K pairs | source 03 §4.4 第 251 行 |
| TO-only LOH 在 paired 下平衡 | **86.5%** | HP_Ratio 0.4-0.6 | source 03 §4.4 第 251 行 |
| Cohen's d (HP_Ratio 跨模式) | **−1.20** | 巨大效應量 | source 03 §4.2.1 第 158 行；§4.4 第 251 行 |
| 跨樣本一致性 (CV-2) | **7/7 樣本** | 同方向 | source 03 §4.4 第 250 行 |
| Simpson's paradox r | **−0.964** | imbalance vs self-phasing fraction | source 03 §4.4 第 251 行 |

## 2. SP1/SP2/SP3 個別位點失衡

| 位點 | Baseline HP2:HP1 | V5 HP2:HP1 | Paired 方向 | 來源 |
|------|---------------|-----------|----------|------|
| **SP1** chr19:17565944 | **113:0** | 113:0（V5 不修） | HP2 | source 01 §1, §4.1；source 03 §5.A.1 |
| **SP2** chr19:12452332 | **109:1** | 109:1 | HP2 | source 01 §1, §4.2 |
| **SP3** chr19:12467180 | **108:0** | 108:0 | HP2 | source 01 §1, §4.3 |

**注意**：V5 不修 self-phasing 本身（Phase 層 V2b 已解）；SP1-3 視覺顯示 Baseline → V2b/V3F/V5 整體 orientation 翻轉至 Paired 方向。

## 3. V5max1/V5max2/V5max3 重分配熱點

| 位點 | V3-F HP33 reads | V5 HP33 reads | Δ HP33 | Δ direction |
|------|---------------|--------------|:------:|:-----------:|
| **V5max1** chr19:4639528 | 39 | 0 | **−39** | **+39（→ HP11）** |
| V5max2 chr19:2235521 | 26 | 0 | −26 | +26（→ HP11）|
| V5max3 chr19:7405500 | 16 | 0 | −16 | +16（→ HP21）|

來源：source 01 §3.1、§5.1 表；source 03 §5.E table。

## 4. ISM 影響 3-tier 分類

| Tier | 特徵數 | 占比 | 代表特徵 | 來源 |
|------|:----:|:----:|---------|------|
| 🔴 嚴重影響（HP-依賴） | **29** | 38% | HP_Ratio、Potential_LOH、HPMergedDelta/Sig、HPFineNGroups | source 03 §4.3 |
| 🟡 中度影響（間接污染） | **14** | 7% | QualityScore、GlobalP、CramersV、VerificationClass | source 03 §4.3 |
| 🟢 不受影響（無 HP 依賴） | **42** | 55% | PairwiseMean/MedianDist、AlleleDelta、Caller、甲基化矩陣、CpG 座標 | source 03 §4.3 |
| **總計** | **85** | 100% | - | - |

**注意**：百分比與 count 比例不一致（14/85 = 16.5% 而非 7%），來源報告本身呈現問題；以 count 為準。

## 5. ISM HP_Ratio LOH 量化

| 指標 | 數值 | 含義 |
|------|------|------|
| ISM HP_Ratio LOH 中 self-phasing artifact | **62%** | TO mode 的 ISM-level LOH 大部分是假象 |
| AF 0.1-0.8 中 self-phasing artifact | **近 100%** | 中等 AF 範圍最嚴重 |
| TO TP 中 paired 平衡比例 | **86.5%** | TO 看到的 LOH 在 paired 下其實平衡 |
| 6,485 個非 LOH 平衡位點 | **100%** | 全部出現此問題 |
| LOH.bed Jaccard (baseline vs PON-only) | **1.0000** | LOH.bed 完全不變 |
| ISM HP_Ratio LOH vs LOH.bed kappa | **0.670** | 兩套系統不一致性 |

來源：source 02 §3.1 表；source 03 §4.2.1。

## 6. V5 修補後驗證

### 6.1 4 項硬性 sanity check（HCC1395 5kHz 15 sites）

| # | 檢查 | 結果 | 來源 |
|---|------|------|------|
| 1 | 守恆律 A · ΔHP33 + (ΔHP11 + ΔHP21) = 0 | **15/15 PASS** | v5_audit_suite/06 §2 |
| 2 | 守恆律 B · germline HP1/HP2 不變 | **15/15 PASS, 0 reads 漂移** | 06 §3 |
| 3 | Layer 1.5 期望 1 · 33→directional 精確守恆 | **15/15 PASS** | 06 §4 |
| 4 | Layer 1.5 期望 2 · 無 germline → HP33 = 0 | **0 violation** | 06 §5 |

### 6.2 Paired ground-truth concordance

| 對照 | Baseline | V5 | Δ | 來源 |
|------|---------|-----|---|------|
| 15-site clean PS（11 sites）| 74.9% | **88.2%** | **+13.3 pp** | v5_audit_suite/07 §3；source 03 §0 |
| 15-site aggregate（all 15 pooled）| 72.20% | 78.85% | +6.65 pp | 07 §2；source 03 §0 |
| **全基因組 clean PS** | **82.2%** | **90.5%** | **+8.3 pp** | PI 報告 4 § 3.7（source 03 §0）|

### 6.3 量化指標變化

| 指標 | Baseline | V5 | Δ | 來源 |
|------|---------|-----|---|------|
| AMB% | 17.5% | **8.0%** | **−9.5 pp** | source 03 §4.2 表 |
| HP:i:33 reads（HCC1395 全基因組）| 239,679 | 110,197 | **−54%** | source 03 §4.2 表 |
| TP Balanced% | 13.0% (V2b) | (V5 同 V3F 基準) | - | source 03 §4.2 表 |

## 7. F1 因果鏈關鍵數字

| 指標 | 數值 | 含義 |
|------|------|------|
| **ClairS-TO raw F1**（不經 ISM）| **0.7166** | **三版本完全相同 — V5 不改 caller** |
| Baseline + ISM SuggestFilter | 0.7157 | 早期實驗值 |
| V3F + ISM | 0.7154 | - |
| V5 + ISM | 0.7154 | - |
| V5 vs Baseline ΔF1 | **−0.0003** | 噪音範圍 |
| V5 vs V3F ΔF1 | +0.0001 | 同樣噪音 |
| 所有版本 F1 區間 | 0.7154 - 0.7166 | 寬度 0.0012 |
| ISM SuggestFilter 對 F1 整體影響 | -0.0009 ~ -0.0013 | 負面（所有版本）|
| ISM F1（單獨對 TO germline FP）| 0.0124 | 對 caller FP 無修復力 |
| Paired F1 | 0.0909 | paired 模式 ISM F1 |

來源：source 03 §4.2.2 + §8 表；PI 報告 4 §3.4。

## 8. 程式碼修改量

| Commit | Hash | 修改檔案 | 行數 | 備註 |
|--------|------|---------|-----|------|
| V2b | `8b8c1fd` | Phasing.cpp / PhasingGraph.cpp / PhasingProcess.cpp | +9/-2、+34/-0、+25/-3 | HaplotagProcess.cpp **未動** |
| V3-Fixed | `41ff147` | HaplotagProcess.cpp（getVote 重寫 + caller 端 enum 修）| +36/-25 | line 506-541, 697 |
| INDEL guard | `380e8d2` | HaplotagProcess.cpp（countINDELHaplotype）| +8/-4 | line 497-510 |
| V5 working tree | (未 commit) | HaplotagProcess.cpp（Layer 1.5 + countSNPHaplotype alt guard）| +24/-7 | line 489-494, 512-563 |
| **V5 累計修改量** | - | 3 函式集中 | **+68/-36 行** | 介面契約 H:66-68 一字未變 |

## 9. 跨樣本與業界對照

### 9.1 跨樣本數據（CV 系列）

| CV-X | 指標 | 結論 |
|------|-----|------|
| CV-1 | imbalance vs self-phasing fraction Spearman | r = -0.964（強負相關，Simpson's paradox）|
| CV-2 | 跨樣本方向一致性 | 7/7 樣本同方向 |
| CV-5 | self-phasing fraction > 0.30 | 4/7 樣本（PARTIAL）|

### 9.2 業界對照（LongPhase-S 2025 paired）

| 指標 | LongPhase-S（paired）| 本實作（TO + PON）|
|------|--------------------|-----------------|
| F1 改善（ClairS SNV）| **+4.5%** | -0.0003（caller 不變，read-level concordance +8.3 pp 是真實價值）|
| F1 改善（ClairS indel）| +7.1% | - |
| F1 改善（DeepSomatic SNV）| +1.2% | - |
| F1 改善（DeepSomatic indel）| +0.5% | - |
| 設計哲學 | somatic 錨在 germline scaffold | 同（PON 替代 normal 作 anchor）|

來源：industry_references.md → bioRxiv 2025.11.20.689492v1。

## 10. 五大研究目標銜接

| 目標 | self-phasing 影響 | 緩解後可解鎖 |
|------|-----------------|-------------|
| 目標 1 — per-CpG 多標籤關聯 | 直接依賴 HP tag | ✓ |
| 目標 2 — clone 結構分析 | 依賴 4-bucket | ✓ |
| 目標 3 — 二次打擊順序 | 間接依賴 | 部分 |
| 目標 4 — TO normal 補強 | 直接依賴（self-phasing 修補是前提）| ✓ |
| 目標 5 — F1 panel | 部分依賴 | 部分 |

來源：`InterSubMod/docs/architecture/20260327_InterSubMod研究願景定錨_01.md` + memory `project_research_vision_five_goals.md`。

## 11. Thread D 主軸量化（已有初步發現）

| 指標 | 數值 | 來源 |
|------|------|------|
| Wilcoxon signed-rank p-value（NG=2 cross-sample 6/6）| **p = 0.0156** | B1 報告 2026-04-23 |
| W 統計量 | 21 | B1 |
| Median gap | 0.365 | B1 |
| Cross-sample direction | 6/6 正向 | B1 |
| Paired control gap | median 0.00003, p = 0.578 | B3 |
| Inner same-hap %（NG=2）| 93-99% | obs18, 6 樣本 |

## 12. fact-check 警示（reviewer agent 發現）

| 項目 | storyboard 寫的 | 來源確認 |
|------|---------------|---------|
| enum 行號 | ~~Util.h:21-25~~ | **Util.h:20-25**（HAPLOTYPE_UNDEFINED=-1 起於 line 20）|
| +8.3 pp | 不夠完整 | **+8.3 pp（clean PS blocks，全基因組）** — 不是全基因組 raw |
| ISM 3-tier 百分比 | 14 個 (7%) | 14 個 / 共 85 features（百分比與 count 不一致來自來源報告）|
| LongPhase-S F1 改善 | "SNV +4.5%" | **"在 ClairS 上 SNV +4.5%、indel +7.1%"** 必須註明哪個 caller |
| "業界共識" | over-claim | 改為 **"同實驗室相鄰工作"** 更準確 |

詳見 reviewer 報告 — 已套用入 storyboard_v2.md。
