# Self-Phasing 背景概念與前後文（v2 PPT 配套）

> 本文件為自包含背景說明，讓只取得 `v2/` 資料夾的人也能理解整體研究脈絡。包含 **v5_audit_suite 6 件主要子報告濃縮**（21 子報告中最 PPT 直接相關的）。

## 1. 一句話總結

> 原始 LongPhase-TO 在 ClairS-TO TO 模式下出現 self-phasing：`--pon-only-phasing=false` 把 somatic variants 當 phasing anchor，下游 `getVote()` 又有 priority bug + enum vs HP tag integer literal mismatch，導致 HCC1395 5kHz 觀察到 99.9% reads 拿 HP21 的極端 imbalance（17.3:1）。本工作在 longphase-to-mod fork（**非 InterSubMod 本 repo**）以 4 個 commit 漸進修補（V2b → V3F → INDEL guard → V5），總 +68/-36 行集中於 3 函式、介面契約零變動。**V5 真實價值在 read-level tag 品質**（clean PS paired concordance +8.3 pp 全基因組），SEQC2 F1 不變是預期（V5 不改 caller、不改 phasing graph，只改 BAM HP tag）。

## 2. 三個系統與其切分（必讀）

```
ClairS-TO (外部 binary) ──→ VCF (FILTER: PASS / NonSomatic / LowQual)
                                        │
tumor BAM ────────────────────────────→ longphase-to-mod (本地 fork, 4-commit 修補)
                                        │     ├─ HaplotagProcess.cpp（V5 Layer 1.5 fallback）
                                        │     └─ outputs:
                                        │         ├─ phased VCF (GT/PS/GT2/GT3)
                                        │         ├─ tumor_tagged.bam (HP:i)  ← V5 改這裡
                                        │         └─ LOH.bed (region-level)
                                        │
                                        ▼
                                        InterSubMod (本 repo, 無 C++ 改動)
                                        ├─ 下游消費 HP:i tag
                                        └─ 計算 ISM 特徵 (HP-依賴)
```

| 系統 | 角色 | 與本工作關係 |
|------|------|--------------|
| **ClairS-TO** | Tumor-Only somatic variant caller（外部 binary）| V5 不修改、raw F1 對所有版本完全相同 |
| **longphase-to-mod** | LongPhase 本地 fork @ `/big7_disk/liaoyoyo2001/longphase-to-mod/`（**獨立 git repo**）| 4-commit 修補在此 |
| **InterSubMod (ISM)** | read-level epigenetic characterization（本 repo）| 下游消費者，**無 C++ 改動** |

**強調**：教授容易誤以為「ISM 改了 phasing」，但實際上 ISM `src/core/` 不包含 HaplotagProcess.cpp；ISM 是被動受惠者。

## 3. HP tag 的兩套編碼（後續 enum bug 的伏筆）

| 系統 | 編碼 | 範例 |
|------|------|------|
| **C++ enum**（`Util.h:20-25`）| `HAPLOTYPE_UNDEFINED=-1, HAPLOTYPE1=1, HAPLOTYPE2=2, HAPLOTYPE1_1=3, HAPLOTYPE2_1=4, HAPLOTYPE3=5` | enum 值是 1/2/3/4/5 |
| **BAM HP:i tag**（外部標籤）| `0=unphased, 1=germline HP1, 2=germline HP2, 11=somatic on HP1, 21=somatic on HP2, 33=somatic ambiguous` | tag 值是 0/1/2/11/21/33 |

**enum vs integer literal mismatch bug**：caller 端 `if(hpResult != HAPLOTYPE1_1)` 用 enum=3 比較，但 hpResult 已是 HP tag integer 11/21/33，型別語意失配 → fallback 死分支永不執行 → HP:i:33 永不出現（baseline）。V3F 修為 integer literal `if(hpResult != 11)`。

## 4. PON-only phasing 是什麼？為何啟用後 bug 才暴露？

| 概念 | 說明 |
|------|------|
| **PON** (Panel of Normals) | 群體 germline 變異資料庫（1000g、CoLoRSdb、dbSNP、gnomAD）|
| **`--pon-only-phasing=false`**（baseline 預設）| phasing graph 用 germline + somatic + unknown 混合 anchor → somatic 反客為主 |
| **`--pon-only-phasing=true`**（V2b 起）| phasing graph 只用 PON-confirmed germline 作 anchor → 解 phase scaffold self-phasing |

**為何 V2b 啟用後 V3F+V5 bug 才暴露**：PON-only 啟用後 germline votes 急減，原本在 paired 模式下被 germline 訊號充足掩蓋的 `getVote()` priority bug + enum bug 立刻顯形。**這是「修一個 bug 暴露另兩個 bug」的因果鏈**。

## 5. 三層獨立 bug 詳解（v5_audit_suite/01_code_diff_analysis.md 濃縮）

### 5.A Phase 層 — self-phasing scaffold（V2b 解）

- **位置**：`PhasingProcess.cpp:55, 154-157` `convertNonGermlineToSomatic()` 不被觸發
- **症狀**：somatic-rich region 內 somatic ALT reads 集中於 sub-clone 一條 hap，phasing graph 邊權重 ≈ 100% 共現 → somatic 互相 phasing 形成 self-phasing block
- **修法**：V2b commit `8b8c1fd` 加 `--pon-only-phasing` flag

### 5.B Tag 層 — getVote priority bug（V3F 解）

- **位置**：`HaplotagProcess.cpp:506-541`（baseline `getVote()`）
- **症狀**：`variantKeys` 順序為 `{HP1_1,HP2_1} → {HP3,HP2_1} → {HP1,HP2}`，任一 somatic vote 觸發 break，germline 完全被忽略
- **為何 paired 沒事**：paired 模式 germline 訊號充足，迴圈通常在 germline pair 命中；TO PON-only 後 germline votes 急減立刻顯形
- **修法**：V3F commit `41ff147` 重寫為兩層：Layer 1 germline first → Layer 2 somatic annotate

### 5.C Tag 層 — enum vs integer literal mismatch（V3F 同 commit 解）

- **位置**：`HaplotagProcess.cpp:697`（caller 端）
- **症狀**：`if(hpResult != HAPLOTYPE1_1)` 用 enum=3 比較，hpResult 已是 HP tag integer 11/21/33；fallback 死分支永不執行；HP:i:33 永不出現（baseline 0/15 vs V3F 8/15 site 有 HP33≥1）
- **修法**：V3F 同 commit 改為 integer literal 11/21/33 比較

### 5.D Tag 層 — V3F 過於保守（V5 解）

- **位置**：`HaplotagProcess.cpp:512-563`（V5 `getVote()`）
- **症狀**：V3F 在 `germlineHP1=germlineHP2=0` 時一律 `hpResult=33`，丟失 `HAPLOTYPE1_1/HAPLOTYPE2_1` 的 phased 方向資訊；HCC1395 全基因組 17.5% reads 拿 HP33（過度保守）
- **修法**：V5 加 Layer 1.5 somatic fallback — 當 germline 平手時參考 somaticHP1/somaticHP2 vote，confidence ≥ 0.6 才採用

## 6. V5 三層投票邏輯（v5_audit_suite/13_phase_vs_tag_algorithm_detail.md 濃縮）

```
Layer 1 (germline first):
    if (germlineHP1 > germlineHP2)  → hpResult = 1   (germline HP1)
    elif (germlineHP2 > germlineHP1) → hpResult = 2   (germline HP2)
    else                            → 進入 Layer 1.5

Layer 1.5 (somatic fallback, V5 新增):
    if (somaticHP1 + somaticHP2 == 0) → hpResult = 0   (unphased)
    elif (somaticHP1 > somaticHP2 + 0.6 * total) → 對齊 germline = HP1 方向
    elif (somaticHP2 > somaticHP1 + 0.6 * total) → 對齊 germline = HP2 方向
    else                                          → hpResult = 33 (ambiguous, 真正的)

Layer 2 (encode hpResult):
    germline=1 + somatic_total>0 → HP:i:11
    germline=2 + somatic_total>0 → HP:i:21
    germline=0 + somatic_total>0 → HP:i:33
    germline only                → HP:i:1 or 2
    nothing                      → HP:i:0
```

## 7. 4-commit 漸進演進

| Commit | 名稱 | 解的 bug | 行數 |
|--------|------|---------|-----|
| `8b8c1fd` | V2b | 5.A Phase scaffold | 改 Phasing.cpp +9/-2、PhasingGraph.cpp +34/-0、PhasingProcess.cpp +25/-3；HaplotagProcess.cpp **未動** |
| `41ff147` | V3-Fixed | 5.B priority + 5.C enum | HaplotagProcess.cpp +36/-25 |
| `380e8d2` | INDEL guard | 補 UB | HaplotagProcess.cpp +8/-4 |
| **working tree（未 commit）** | **V5** | 5.D 過度保守 | HaplotagProcess.cpp +24/-7（Layer 1.5 + countSNPHaplotype alt guard）|

**總計**：+68/-36 行集中於 3 函式（getVote / countSNPHaplotype / countINDELHaplotype）；介面契約 `HaplotagProcess.h:66-68` 三函數簽章一字未變。

## 8. 量化驗證（v5_audit_suite/06_v5_sanity_bug_check.md + 07_paired_ground_truth_concordance.md 濃縮）

### 8.1 4 項硬性 sanity check（15/15 PASS, 0 violation）

| # | 檢查 | 結果 |
|---|------|------|
| 1 | **守恆律 A · Δ-consistency**：ΔHP33 + (ΔHP11 + ΔHP21) = 0 | 15/15 PASS |
| 2 | **守恆律 B · Germline 不變**：V3F 與 V5 的 HP1/HP2 read 數逐 site 比對 | 15/15 PASS, 0 reads 漂移 |
| 3 | **Layer 1.5 期望 1 · 33→directional 精確守恆**：ΔHP11 == n(V3F=33→V5=11)、ΔHP21 == n(V3F=33→V5=21) | 15/15 PASS（V5max1=34 reads 精確守恆）|
| 4 | **Layer 1.5 期望 2 · 無 germline → HP33**：跨 15 sites pool transition table | 0 reads violation |

### 8.2 paired ground-truth concordance

| 對照 | Baseline | V5 | Δ |
|------|---------|-----|---|
| 15-site clean PS（11 sites）| 74.9% | **88.2%** | +13.3 pp |
| 15-site aggregate（all 15 sites pooled）| 72.20% | 78.85% | +6.65 pp |
| **全基因組 clean PS（PI 報告 4）** | **82.2%** | **90.5%** | **+8.3 pp** |

## 9. F1 變動的因果鏈釐清（v5_audit_suite/08_synthesis_conclusions.md + 12_gt_distribution_audit.md 濃縮）

**核心邏輯（PI 報告 4 §3.4）**：

```
ClairS-TO (caller) ──→ VCF
                       │
                       │ raw calling F1 = 0.7166（caller 層真實 F1）← V5 不改此處
                       ▼
                   longphase-to-mod ──→ tumor_tagged.bam (HP:i tag)
                       │
                       │ V5 改的是 BAM HP tag，不是 VCF
                       ▼
                   InterSubMod ISM ──→ region-level features (HP_Ratio, HPSig, ...)
                       │
                       │ HP tag 變了 → region 特徵變了
                       ▼
                   ISM SuggestFilter（src/core/RegionProcessor.cpp:1120, 1269）
                       │
                       │ 對每個 variant 計算「是否該標 LOW_QUAL」的 SuggestFilter
                       │ ←── V5 對 F1 的影響進入這裡
                       ▼
                   套用 SuggestFilter → 最終 SEQC2 F1 = 0.7154
                       │
                       │ Raw 0.7166 → Baseline+ISM 0.7157 → V3F+ISM 0.7154 → V5+ISM 0.7154
                       │ V5 vs Baseline = -0.0003（噪音）
```

| 證據 | 數字 | 含義 |
|------|------|------|
| ClairS-TO raw F1（所有版本）| **0.7166** 完全相同 | V5 不改 caller |
| ISM SuggestFilter 對 F1 的整體影響 | -0.0009 ~ -0.0013（負面，所有版本）| ISM 過濾本身略微負面 |
| V5 vs Baseline 的 F1 差距 | -0.0003 | 統計噪音範圍 |
| V5 vs V3F 的 F1 差距 | +0.0001 | 同樣噪音 |
| 所有版本 F1 區間 | 0.7154 - 0.7166（寬度 0.0012）| 任何版本之間都是噪音級差距 |

**為何 F1 不能衡量 V5 品質**：ClairS-TO 的 FP **主要是 germline variants**，非 somatic；ISM 甲基化分析設計用來區分 subclone 結構，不是區分 germline vs somatic；**V5 修正不改變此根本限制**。F1 變動 100% 來自 ISM SuggestFilter；SuggestFilter 對 ClairS-TO 的 FP 過濾力先天不足。

## 10. ISM 特徵 3-tier 影響分類（v5_audit_suite/04_imbalance_ratio_analysis.md 濃縮）

| 影響等級 | 特徵數（佔比）| 代表特徵 | TO 結果處理 |
|---------|--------------|---------|-----------|
| 🔴 **嚴重影響**（直接依賴 HP）| **29 個（38%）** | HP_Ratio → 假 LOH；Potential_LOH → 62% artifact；HPMergedDelta/Sig → 方向反轉；HPFineNGroups（含 NG=2 LOH-constrained phasing） | **必重跑** V5 BAM；舊結論需加註版本 |
| 🟡 **中度影響**（間接污染）| **14 個（7%）** | QualityScore → AUC 0.497（已移除 LOH penalty）；GlobalP → 取 HP/Allele 最小值；CramersV → 取 HP/Allele 最大值；VerificationClass → label_sig 含 HP 成分 | 重評；多數影響微弱 |
| 🟢 **不受影響**（無 HP 依賴）| **42 個（55%）** | PairwiseMean/MedianDist；AlleleDelta / AlleleP；Caller 特徵（AF / GQ / DP / SB）；甲基化矩陣（BAM MM/ML tag）；CpG 座標、region_methyl_mean | **結論穩固，不需重測** |

總數 **85 features**。

## 11. LOH 兩層的精確區分（必讀，最易混淆）

| LOH 層次 | 路徑 | self-phasing 影響 | 數字 |
|---------|-----|-----------------|------|
| **ISM HP_Ratio LOH** | BAM HP tag → ISM HP_Ratio<0.1 or >0.9 | **嚴重**：62% 是 artifact；AF 0.1-0.8 近 100% artifact；TO TP 中 86.5% 在 paired 下完全平衡 | 62% ISM-level LOH artifact |
| **LOH.bed region-level LOH** | VCF allele depth (AD) → LongPhase region detection | **零**：PON-only phasing LOH.bed Jaccard=1.0，完全不受 self-phasing 影響 | LOH.bed 不變 |

**關鍵**：兩套系統使用**不同數據源**（BAM HP tag vs VCF AD），kappa = 0.670 的不一致性由此解釋。本工作 self-phasing 修補只影響 BAM HP tag 路徑（ISM HP_Ratio LOH），**不改變 LOH.bed**。

## 12. 跨樣本一致性（v5_audit_suite + 03 報告 §4.4 濃縮）

| 指標 | 數值 |
|------|------|
| **CV-2 方向一致性** | **7/7 樣本同方向** |
| 同位點 HP_Ratio 跨模式 r | r = 0.001（n=288K pairs，TO 與 paired 完全不相關）|
| TO-only LOH 在 paired 下完全平衡 | 86.5%（HP_Ratio 0.4-0.6）|
| **Cohen's d** | **−1.20**（巨大效應量） |
| Simpson's paradox（跨樣本） | r = −0.964（imbalance vs self-phasing fraction） |
| CV-5 self-phasing fraction > 0.30 | 4/7 樣本（PARTIAL — HCC1395 / DORADO / HCC1937 結構性 LOH 主導）|

## 13. 與五大研究目標的關聯（`InterSubMod/docs/architecture/20260327_InterSubMod研究願景定錨_01.md`）

1. **目標 1 — per-CpG 甲基位點多標籤關聯性評分**：每個 CpG 與 ALT/REF/HP1/HP2 的關聯 → **依賴 HP tag 正確（受 self-phasing 影響）**
2. **目標 2 — clone 結構分析**：sub-clone + 共演化 → **依賴 4-bucket 分群（受影響）**
3. **目標 3 — 二次打擊與事件順序推論**：依賴目標 1/2 + LOH/CNV/HP → **間接依賴**
4. **目標 4 — TO normal 資訊補強**：無配對 normal 下估計 normal 背景 → **直接依賴（self-phasing 修補是前提）**
5. **目標 5 — 整合 evidence panel 提升 F1**：補充 caller、保留 TP、過濾 FP → **部分依賴**

**Self-phasing 修補是 5 大目標解鎖的前提**（特別是目標 1/2/4）。

## 14. 已知開放議題（R1-R9，v5_audit_suite/00_INDEX.md 濃縮）

| ID | 議題 | 緩解 |
|----|-----|------|
| **R1** | cnLOH 雙親同源 region 仍未解 | 需 cnLOH-aware filter（CN 層 + germline-only methylation reference）|
| **R2** | AF > 0.9 真結構性 LOH 與 phasing artifact 邊界 | V5 在 0.1–0.8 表現好；AF > 0.9 仍可能共存 |
| **R3** | 15 sites cherry-picked | 全基因組數字參照 PI 報告 4（V5=90.5% / BL=82.2%）|
| **R4** | Confidence threshold 0.6 未直接驗證 | 需 V5 binary 加 vote log 或 IGV session 看 PS block ALT/REF 投票 |
| **R5** | V5 working tree uncommitted | 切 2 獨立 commits 為短期 P0 行動 |
| **R6** | 7 樣本擴展未做 | 中優先後續（COLO829、DORADO、HCC1937、HCC1954、H1437、H2009 + HCC1395）|
| **R7** | 部分 metric V5 略遜 | L4 orientation-corrected 受 problem PS 影響；imbalance ratio 受 outlier 主導 |
| **R8** | Problem PS / low_germ_n 不適用 read-level | A_TP02、D_SP1、D_SP3 為 PS orientation pick artifact |
| **R9** | Paired 自身仍用 LongPhase | 理想需 trio-phased 第二 ground truth |

## 15. 後續行動 F1-F8（v5_audit_suite/00_INDEX.md 濃縮）

| ID | 動作 | 優先 |
|----|------|-----|
| **F1** | commit V5 working-tree 修改（切 2 獨立 commits）| **高** |
| **F2** | Confidence threshold 0.6 投票 log 驗證 | 中 |
| **F3** | 7 樣本 V5 BAM 全量重跑 | 中 |
| **F4** | P4 master dataset × 兩 flag 重跑驗證 HPFineNGroups marker | 低 |
| **F5** | manifest.yaml 加 `haplotag_version` 欄位 | 中 |
| **F6** | cnLOH 區獨立評估方案 | 中 |
| **F7** | trio-phased 第二 ground truth | 低 |
| **F8** | 跨 50-100 隨機抽樣位點 cross-validate | 中 |

## 17. ★ Purity threshold 0.95 觸發鎖機制（v3 新增 · PPT S11 核心）

**位置**：`/big7_disk/liaoyoyo2001/longphase-to-mod/PhasingProcess.cpp:197`

**邏輯**：

```cpp
if (purity > 0.95) {
    // Two-Pass 路徑（兩條路）
    vGraph->convertNonGermlineToSomatic();  // Pass 1: 純 somatic
    vGraph->phasingProcess(..., nullptr);   // Pass 2: 重跑 phasing
}
else {
    // Baseline 標準流程（三條路）
    vGraph->somaticCalling(...);
    vGraph->phasingProcess(..., &chrInfo.ploidyRatioMap);  // 用 somatic+germline 混合
}
```

**HCC1395 5kHz 觸發狀態**：
- 真實 purity > 0.95（純 tumor）
- Baseline 估計 purity = **0.927**（四捨五入 0.93；非誤判，是邊緣 case）
- 0.927 ≤ 0.95 → **未觸發 Two-Pass** → 走「三條路」標準流程
- 標準流程假設樣本含 normal contamination → 但實為純樣本 → somatic 進 phasing graph 反客為主 → 暴露 tag 層 somatic-first 投票 bias
- 結果：**self-phasing 17.3:1 artifact**

**V5 修正**：
- V5 用 `--pon-only-phasing` flag 強制啟用 Two-Pass（不管 purity）
- V5 PON-only 自身有 purity calculator failure（ploidyRatioMap=nullptr → q1=q3=0 → purity=0），2026-04-29 已修復（新增 `collectPloidyRatio()` 於 `PhasingGraph.cpp:1147-1175`）
- 修復後雙 sample 驗證：0.93 sample → 0.871（vs baseline 0.927 誤差 -0.06）；0.6 sample → 0.634（vs baseline 0.607 誤差 +0.03）；**不需重跑 tag**（修復 side-effect free）

**PPT 整合位置**：S11 root-cause tree 把這個觸發條件設為「主者 / 根」，三 bug 為「葉」。

來源：`v2/source_materials/04_purity_calculator_failure_root_cause.md`（v5_audit_suite/18 複製）

---

## 18. v5_audit_suite 21 子報告對應（哪份是 PPT 主要來源）

| 子報告 | 主題 | 本 PPT 用到 |
|-------|------|-----------|
| 00_INDEX.md | 母索引 | 整體 |
| **01_code_diff_analysis.md** | V3F+V5 commits 逐 diff | **§5, §6, S15 程式碼 diff** |
| 02_read_intersection_concordance.md | per-read 比對 | S20 caveat |
| 03_hp_family_vs_exact.md | HP family vs exact 編碼 | S5 HP tag 定義 |
| **04_imbalance_ratio_analysis.md** | 17.3:1 量化 | **S5, §10 ISM 3-tier** |
| 05_per_site_improvement.md | 逐 site 改善 | S18 V5max1 |
| **06_v5_sanity_bug_check.md** | 4 項硬性檢查 15/15 PASS | **S16, §8.1** |
| **07_paired_ground_truth_concordance.md** | clean PS +13.3 pp | **S17, §8.2** |
| **08_synthesis_conclusions.md** | 5 PI 問題答覆綜合 | **§9 F1 釐清** |
| 09_purity06_simulation.md | 0.6 純度模擬 | R6 caveat |
| **10_somatic_bias_explanation.md** | 6-panel + 3 IGV 真截圖 | **S5 fig01d, S18-19 IGV** |
| 11_code_issues_inventory.md | 已知 code issue | R5 caveat |
| 12_gt_distribution_audit.md | GT 分布審計 | §9 F1 釐清補強 |
| 13_phase_vs_tag_algorithm_detail.md | phase vs tag 演算法細節 | §6 V5 三層 |
| 14_user_report_integrated.md | 用戶整合報告 | 主敘事 |
| 15_software_engineering_perspective.md | 工程觀點 | 介面契約 |
| 16-21 | 設計一致性、boldness 等 | caveat |

**6 件主要子報告**（粗體標 ★）：01、04、06、07、08、10。本背景文件已濃縮這 6 份。

---

## 文件結束

如需更深細節，閱讀順序：
1. `v2/source_materials/03_longphase_TO_vs_V5_技術報告.md` 全文（13 段結構）
2. `v2/source_materials/01_V5_IGV_session_visual_audit.md`（IGV 視覺證據）
3. `v2/source_materials/02_Self_Phasing_Baseline_V5_Audit.md`（一句結論審核版）

如需 fact-check：見 `v2/notes/key_metrics_table.md`、`v2/notes/code_references.md`。
如需術語補充：見 `v2/notes/glossary.md`。
如需教授提問預備：見 `v2/notes/qa_11_questions.md`。
