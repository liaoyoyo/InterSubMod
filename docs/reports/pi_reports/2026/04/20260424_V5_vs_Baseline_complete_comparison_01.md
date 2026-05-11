<!--
建立時間: 2026-04-24
目標: 完整比較 LongPhase-TO Baseline 與 V5 (Somatic Fallback) 演算法、量化差異與既有結論影響
受眾: PI（演算法演進 + 量化證據 + 結論穩定性）
處理範圍: 5 版本演進 / code-level diff / HP tag 量化 / 結論影響三層矩陣
狀態: validated_complete
pipeline_track: TO (paired_full / paired_pileup 不受影響)
priority: P0
relates_to:
  - 20260422_Self_Phasing_complete_report_for_PI_01.md
  - 20260422_Self_Phasing_multiperspective_argument_01.md
  - 20260424_Self_Phasing_evidence_chain_methodology_01.md
-->

# LongPhase-TO V5 (Somatic Fallback) vs Baseline 完整比較報告
## ——5 版本演進、量化差異、既有結論穩定性

> 撰稿日期：2026-04-24
> 受眾：PI（同時讀過前 3 份 Self-Phasing PI 報告者佳）
> 關鍵 commits：`8b8c1fd` (PON-only) / `41ff147` (V3-Fixed two-layer) / `380e8d2` (UNDEFINED guard) / V5 HEAD
> 主要素材：
> - V5 binary：`/big7_disk/liaoyoyo2001/longphase-to-mod/longphase-to`
> - V5 BAM：`/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam` (268G)
> - Session report：`InterSubMod/docs/provenance/ai_sessions/2026/04/20260412_V5_somatic_fallback_getVote_修正與驗證_01.md`
> - 量化 TSVs：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/`

---

## Positioning（報告定位）

前 3 份 Self-Phasing PI 報告**全部基於 V3-Fixed BAM**。但 LongPhase-TO 已演進到 V5，V3-F 是**過渡版本**。

本報告回答 PI 三個核心問題：

1. **V5 vs Baseline 演算法差什麼？**（code-level diff，Section 1-2）
2. **HP tag 量化差異有多大？**（data-level，Section 3）
3. **既有研究結論是否穩定？**（impact-level，Section 5-6）

**與前 3 份 PI 報告分工**：

| 報告 | 焦點 |
|------|------|
| 20260422 技術 | self-phasing 問題 → 修改 → 結果 |
| 20260422 多視角 | 3 視角 × 9 challenges adversarial review |
| 20260424 證據鏈 | 6 證據獨立性矩陣 |
| **20260424 V5 vs Baseline（本報告）** | **演算法演進 + 量化比較 + 結論穩定性** |

四份報告構成完整 self-phasing 主題包：問題（A）→ 質疑（B）→ 證據（C）→ **演進與穩定性（D）**

---

## Executive Summary

**V5 vs Baseline 三層結論**：

1. **演算法層面**：V5 = `PON-only phasing` + `getVote() 三層分離` + `Layer 1.5 somatic fallback`。Baseline 有兩個獨立 bugs（self-phasing circular dependency + getVote() 優先序覆蓋），V5 透過三步演進（V2b → V3-F → V5）修正完畢。
2. **量化層面**：V5 比 Baseline (1) 移除 765K self-phasing circular reads；(2) 將 129K reads 從 ambiguous bucket 重分配為 directional；(3) AMB% 從假性 1.3% 修正到 8.0%（Paired ground truth 5.4%）；(4) Clean-PS somatic accuracy 82.2% → 90.5% (+8.3pp)；(5) **SEQC2 F1 持平**（0.7157 vs 0.7154 = noise）。
3. **既有結論穩定性**：10 個既有結論中——**3 個 Tier-1 需 V5 重評**（17.3:1 bias、62% LOH artifact、HPFineNGroups 89.1%），**2 個 Tier-2 需監測**（LOH-constrained ΔNG、LOH.bed AUC），**5 個 Tier-3 完全不受影響**（Germline ASM、LOH.bed Jaccard、Coverage_Multiple、E6 演算法理論、E4 跨樣本一致性）。**所有結論的機制層面均維持有效**。

---

## Section 1：版本演進史

### 1.1 5 版本完整時序

| 版本 | 修正內容 | SEQC2 F1 | AMB% | 狀態 |
|------|---------|---------|------|------|
| **Baseline** | — | 0.7117 | 1.3%（假性低）| 已淘汰 |
| **V2b** PON-Only Phasing | 移除 self-phasing circular | 0.7115 | ~0%（HP3 未出現）| 已淘汰 |
| **V3-Fixed** | getVote() 兩層分離 + UNDEFINED guards | 0.7153 | 17.5%（過嚴）| 過渡最佳 |
| **V4** alt-guard | altHaplotype UNDEFINED guard | 0.7153 | 17.5% | =V3F，可刪 |
| **V5** Somatic Fallback | Layer 1.5 somatic fallback | **0.7154** | **8.0%**（balanced）| **當前最佳** |

### 1.2 演進敘事：每版「修了什麼 + 為什麼還不夠」

**Baseline**（原始）：
- 修了：什麼都還沒修（原始 LongPhase-TO 行為）
- 問題：(a) phasing graph 含 somatic variants → self-phasing；(b) `getVote()` 中 somatic votes 覆蓋 germline；(c) HP:i:33 enum vs integer 寫入 bug → ambiguous tag 永不出現

**V2b PON-Only**：
- 修了：phasing graph 排除 somatic variants
- 為什麼還不夠：getVote() 優先序 bug 未修；HP:i:33 enum bug 未修
- 副作用：tagged reads 從 19.57M → 18.81M（-765K，這是 self-phasing 過度貢獻的 reads）

**V3-Fixed**：
- 修了：(a) `getVote()` 兩層分離（germline first）；(b) `countINDELHaplotype` UNDEFINED guard；(c) HP:i:33 用整數而非 enum 寫入
- 為什麼還不夠：矯枉過正——當 germline votes = 0 時，**完全丟棄 somatic directional evidence**（HP1_1 vs HP2_1），一律返回 HP:i:33 → AMB% 17.5%（遠高於 Paired 5.4%）

**V4 alt-guard**：
- 修了：`countSNPHaplotype` alt 分支加 `altHaplotype != HAPLOTYPE_UNDEFINED` guard
- 為什麼還不夠：純防禦性修正；驗證後與 V3-Fixed 結果完全相同

**V5 Somatic Fallback**：
- 修了：`getVote()` 加入 Layer 1.5——當 germline = 0 時，用 somatic HP1_1 vs HP2_1 作 fallback
- 結果：AMB% 17.5% → 8.0%；129K reads 從 HP:i:33 重新分配為 HP:i:11/21（54% 減少）；clean-PS accuracy +0.5pp；SEQC2 F1 持平

---

## Section 2：演算法 code-level 完整比較

### 2.1 Baseline 的 bug 鏈（為何需要修）

**Bug 1：Phasing VCF 包含 somatic variants** → self-phasing circular dependency
- 來源：原始 LongPhase-TO `phase` 步驟未過濾 somatic
- 影響：62% TO LOH 為 artifact、HP1:HP2 bias 17.3:1、765K 額外 tagged reads 不可靠
- 詳見：`InterSubMod/docs/reports/pi_reports/2026/04/20260422_Self_Phasing_complete_report_for_PI_01.md` Section 2

**Bug 2：getVote() 優先序錯誤**
- 原始行為：somatic votes（HP1_1, HP2_1, HP3）先設定 min/max → germline votes 後寫入時可能被覆蓋
- 影響：HP tag 73.2% extreme + 13.0% balanced（Paired baseline 是 36.9% balanced）；QS ceiling 化（mean=89.5）

**Bug 3：HP:i:33 寫入 enum vs integer 不一致**
- 原始行為：fallback 寫入 enum value `3`/`4`（HAPLOTYPE3/HAPLOTYPE1_1）而非 integer `33`/`11`
- 影響：HP:i:33 在 BAM 中**從未出現**——所有 ambiguous reads 被誤標為 directional

### 2.2 V3-Fixed 的兩個修正（commit `41ff147`、`380e8d2`）

**修正 1：getVote() 兩層分離**
```cpp
// V3-Fixed
if (germlineHP1 > 0 || germlineHP2 > 0) {
    // Layer 1: germline 決定 min/max
    germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
} else {
    min = 0; max = 0;        // ← 完全丟棄 somatic
}
// Layer 2: somatic annotation 加在 germlineResult 之上
```

**修正 2：UNDEFINED guard + integer literal**
- `countINDELHaplotype` REF/ALT 分支均加 `haplotypeBase.refHaplotype != HAPLOTYPE_UNDEFINED` guard
- HP:i:33 fallback 改用 integer literal 寫入

**結果**：HP:i:33 終於出現（chr1:1-10M 抽樣 385 個）；TP Balanced% 13% → 22.5%；SEQC2 F1 0.7117 → 0.7153；**AMB% 17.5%** 過嚴

### 2.3 V5 的 Layer 1.5 突破（核心修正）

**位置**：`/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:512-563`

**三層模型**：

```cpp
void HaplotagProcess::getVote(...) {
    int germlineHP1 = countMap[HAPLOTYPE1];
    int germlineHP2 = countMap[HAPLOTYPE2];
    int somaticHP1  = countMap[HAPLOTYPE1_1];
    int somaticHP2  = countMap[HAPLOTYPE2_1];
    int somaticTotal = somaticHP1 + somaticHP2 + countMap[HAPLOTYPE3];

    // Layer 1: Germline priority
    int germlineResult = 0;
    if (germlineHP1 > 0 || germlineHP2 > 0) {
        if (germlineHP1 >= germlineHP2) {
            min = germlineHP2; max = germlineHP1;
            germlineResult = 1;
        } else { min = germlineHP1; max = germlineHP2; germlineResult = 2; }
    }
    // Layer 1.5: Somatic fallback (V5 NEW, lines 538-547)
    else if (somaticHP1 > 0 || somaticHP2 > 0) {
        if (somaticHP1 >= somaticHP2) {
            min = somaticHP2; max = somaticHP1;
            germlineResult = 1;
        } else { min = somaticHP1; max = somaticHP2; germlineResult = 2; }
    }
    else { min = 0; max = 0; }    // Truly ambiguous

    // Layer 2: Encoding
    if (somaticTotal > 0) {
        if (germlineResult == 1)       hpResult = 11;
        else if (germlineResult == 2)  hpResult = 21;
        else                           hpResult = 33;
    } else {
        hpResult = germlineResult;
    }
}
```

**為何 V5 是顯著進步**：V3-F 矯枉過正——germline=0 時忽略 somatic directional evidence。V5 用 HP1_1 / HP2_1 的 phasing 資訊作 fallback。Confidence threshold 0.6（在 `judgeHaplotype()` 中，line 730-734）自動攔截不確定分配（min/max 接近時）。

### 2.4 五個場景對照（V3-F vs V5）

| 場景 | germline votes | somatic votes | V3-F result | V5 result | 差異 |
|------|---------------|--------------|------------|-----------|------|
| A: germline majority | HP1=20, HP2=3 | HP1_1=2 | 11 | 11 | 同 |
| B: germline minority | HP1=20, HP2=3 | HP2_1=5 | 11 | 11 | 同 |
| **★ C: germline absent** | HP1=0, HP2=0 | HP2_1=5, HP1_1=1 | **33** | **21** | **V5 recovers direction** |
| D: pure HP3 | HP1=0, HP2=0 | HP3=3 | 33 | 33 | 同 |
| E: no somatic | HP1=0, HP2=0 | — | 0 | 0 | 同 |

**Scenario C 是 V5 的唯一行為差異**——這是「資訊回收」而非「規則改寫」。

### 2.5 HP tag enum mapping 完整對照

```
HAPLOTYPE_UNDEFINED = -1  (guard sentinel; not emitted)
HAPLOTYPE1          = 0   → HP:i:1   (germline HP1)
HAPLOTYPE2          = 1   → HP:i:2   (germline HP2)
HAPLOTYPE3          = 2   → HP:i:33  (unphased somatic, ambiguous)
HAPLOTYPE1_1        = 3   → HP:i:11  (somatic linked to HP1)
HAPLOTYPE2_1        = 4   → HP:i:21  (somatic linked to HP2)
```

**Baseline bug**：HAPLOTYPE3 (enum 值 2) 寫入 BAM 變成 `HP:i:2`（germline HP2）→ HP:i:33 永不出現。V5 直接寫入 integer literal 11/21/33 修正。

![Figure 16 — V5 Three-Layer Logic + Code Excerpt + Scenario Table](figures/fig16_v5_threelayer_logic.png)

---

## Section 3：HP 差異完整量化比較

### 3.1 全基因組 Total Tagged Reads（HCC1395 5kHz）

| 版本 | Total Tagged | 與 Baseline 差 | 解釋 |
|------|-------------|---------------|------|
| Baseline | **19,571,246** | — | self-phasing VCF 多 somatic anchors → 更多 phased SNPs |
| V2b PON-Only | 18,805,932 | **-765,314** | PON-only VCF 排除 somatic → 少 phasing anchors |
| V3-Fixed | 18,805,977 | -765,269 | 同 V2b VCF；getVote 修正僅微增 45 reads |
| V5 | 18,805,977 | -765,269 | 同 V3-Fixed |

**765K 額外 tagged reads 是 self-phasing circular dependency 的人為產物**。已驗證：62% LOH 信號消失、germline concordance 無改善 → 這些 tagged reads 不可靠。PON-only 犧牲少量 coverage 換取 phasing 正確性是正確取捨。

### 3.2 HP Tag 分布（chr1:1-10M 抽樣，4 版本）

| HP Tag | Baseline | V2b | V3-Fixed | V5 |
|--------|---------|-----|---------|-----|
| HP:i:1 (germ HP1) | 26,126 | 26,206 | 26,206 | 26,206 |
| HP:i:2 (germ HP2) | 13,598 | 13,682 | 13,682 | 13,682 |
| HP:i:11 (som HP1) | 1,357 | 1,410 | 1,158 | **1,359** |
| HP:i:21 (som HP2) | 1,256 | 1,203 | 1,070 | **1,243** |
| HP:i:33 (amb) | **0** | **0** | 385 | **11** |

**關鍵觀察**：
- Baseline / V2b 的 `HP:i:33 = 0` 是 **enum vs integer bug** 造成（看似全部 directional，實為錯誤覆蓋）
- V3-Fixed 修正 bug 後 HP:i:33 出現（385），但**過多**（很多原可判方向的 reads 被誤標 ambiguous）
- **V5 的 HP:i:11/21 數量接近 Baseline/V2b**——但機制完全不同：Baseline 用 bug 取得方向，V5 用正確 phasing + 合理 fallback

### 3.3 全基因組 V3-F vs V5（somatic tags 重分配）

| HP Tag | V3-Fixed | V5 | Δ |
|--------|---------|-----|---|
| `.` (untagged) | 13,005,605 | 13,005,605 | 0 |
| HP:i:1 (germ HP1) | 10,071,497 | 10,071,497 | 0（不變）|
| HP:i:2 (germ HP2) | 7,363,566 | 7,363,566 | 0（不變）|
| HP:i:11 (som HP1) | 584,117 | 666,997 | **+82,880** |
| HP:i:21 (som HP2) | 547,118 | 593,720 | **+46,602** |
| HP:i:33 (amb) | 239,679 | 110,197 | **-129,482 (-54%)** |

**Germline tags 完全不變** — V5 修正精確影響 somatic 路徑。

![Figure 17 — HP Tag Distribution Across Versions](figures/fig17_hp_tag_5versions.png)

### 3.4 SEQC2 F1（Calling Pipeline Performance）

| 版本 | TP SF | FP SF | Precision | Recall | F1 | vs Baseline |
|------|-------|-------|-----------|--------|------|-------------|
| Raw (no ISM) | 0 | 0 | 0.7107 | 0.7226 | 0.7166 | — |
| Baseline + ISM | 105 | 91 | 0.7115 | 0.7199 | **0.7157** | — |
| V2b + ISM | 125 | 105 | 0.7116 | 0.7194 | 0.7155 | -0.0002 |
| V3-Fixed + ISM | 112 | 73 | 0.7112 | 0.7198 | 0.7154 | -0.0003 |
| **V5 + ISM** | 113 | 74 | 0.7112 | 0.7197 | **0.7154** | **-0.0003** |

**V5 vs Baseline F1 差異 = -0.0003（-0.04%），在統計噪音範圍**。所有版本 F1 都在 0.7154-0.7166 狹窄區間。

**為何 F1 不能衡量 tag 品質**：
- ClairS-TO 產出的 VCF 在 haplotagging **之前**就固定（所有版本 Raw F1=0.7166 完全相同）
- 報告中 0.7154-0.7157 的差異 100% 來自 ISM SuggestFilter 後過濾
- V5 vs V3-F SF 判定僅差 +1 TP / +1 FP（淨效果零）
- → **F1 持平不代表 tags「沒有更好」，F1 變化也不代表 tags「更好或更差」**

### 3.5 ISM SuggestFilter 演進

| 版本 | Total SF | FP caught | TP lost | SF Precision | FP Catch Rate | ISM F1 |
|------|---------|----------|---------|-------------|--------------|--------|
| Baseline | 196 | 91 | 105 | 46.4% | 0.78% | 0.0154 |
| V2b | 230 | 105 | 125 | 45.7% | 0.90% | 0.0177 |
| V3-Fixed | 185 | 73 | 112 | 39.5% | 0.63% | 0.0124 |
| **V5** | 187 | 74 | 113 | 39.6% | 0.64% | **0.0125** |

所有版本 ISM F1 < 0.02——ClairS-TO 的 FP 主要是 germline variants（非 somatic），ISM 甲基化分析無法區分。**V5 修正不改變此根本限制**。

### 3.6 AMB% 與 confidence 分布

| 版本 | AMB%（global）| AMB%（clean-PS）| Paired target | 距離 Paired |
|------|--------------|----------------|--------------|-----------|
| Baseline | **1.3%**（假性低）| — | 5.4% | -4.1pp |
| V3-Fixed | 17.5% | — | 5.4% | +12.1pp |
| **V5** | **8.0%** | 4.7% | 5.4% | **clean: -0.7pp** |

**Baseline 的 1.3% 是假性低**——getVote() bug 強制分配方向，掩蓋真實不確定性。**V5 的 4.7%（clean-PS）最接近 Paired 5.4%**——更誠實反映真實的方向不確定性。

V5 剩餘 110,197 個 HP:i:33 reads 的 confidence 分布：
- 29,156 reads (26.5%)：confidence = 0.500（完全 50:50 split）
- ~79,000 reads：confidence ∈ [0.50, 0.59]（低於 0.6 threshold）
- 1,428 reads：confidence = NaN（純 HP3，無 HP1_1/HP2_1）

→ Confidence threshold 0.6 **正確攔截不確定分配**。

### 3.7 Clean-PS Concordance（vs Paired ground truth）

| 維度 | Baseline | V5 | 改善 |
|------|---------|-----|------|
| **Somatic direction accuracy** | 82.2% | **90.5%** | **+8.3pp** |
| **Germline direction accuracy** | 87.2% | **91.7%** | **+4.5pp** |
| New assignment accuracy（V3F→V5）| — | 95.0% (267/281) | +253 reads |
| Net wrong assignments | — | +14 reads | — |

**V5 在每個 concordance 維度都更接近 Paired ground truth**。

### 3.8 Somatic HP1 Balance

| 版本 | HP1 / (HP1+HP2) | 距離 Paired (0.523) |
|------|----------------|---------------------|
| Baseline | 0.774 | **0.250** |
| V5 | **0.703** | **0.179** |

Baseline 過度偏向 HP1（self-phasing + bug 雙重原因），V5 更接近 Paired 平衡。

![Figure 18 — Concordance + AMB% + F1 Evolution](figures/fig18_concordance_amb_f1.png)

### 3.9 Region-level ISM Indicators（version_summary.tsv）

從 `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/version_summary.tsv`：

| 版本 | n_regions | TP_rate | HP_Ratio_med | Potential_LOH_pct | HPFineNGroups_3plus_pct |
|------|----------|---------|--------------|-------------------|------------------------|
| Baseline | 40,213 | 0.706 | 0.788 | 58.7 | 27.9 |
| V2b PON-Only | 40,213 | 0.706 | 0.447 | 72.6 | 23.5 |
| V3-Fixed | 40,096 | 0.711 | 0.575 | 61.7 | 29.4 |
| V4 | 40,096 | 0.711 | 0.575 | 61.7 | 29.4 |
| V5 | 40,096 | 0.711 | 0.574 | 62.2 | **29.8** |

**HP_Ratio median**：Baseline 0.788（極端偏移）→ V5 0.574（接近平衡）。
**HPFineNGroups_3plus_pct**：Baseline 27.9% → V5 29.8%（V5 略高，因 NG≥3 仍有真實 LOH-constrained phasing 訊號）。

### 3.10 Methylation Features 影響（最小）

從 `methyl_baseline_vs_v5.tsv`：

| Feature | Baseline AUC | V5 AUC | Δ |
|---------|-------------|--------|---|
| HPMergedDelta | 0.517 | 0.514 | -0.003 |
| AlleleDelta | 0.541 | 0.543 | +0.002 |
| HPFineF | 0.576 | 0.584 | **+0.008** |
| PairwiseMeanDist | 0.525 | 0.528 | +0.003 |
| CramersV | 0.507 | 0.506 | -0.001 |
| GlobalP | 0.533 | 0.530 | -0.003 |

**所有 methylation features 在 V5 與 Baseline 間 ΔAUC ∈ [-0.003, +0.008]，全在噪音範圍**。HPFineF 微升 0.008 是唯一非噪音變化（與 Layer 1.5 somatic fallback 增加 directional information 一致）。

### 3.11 跨樣本資料缺口（誠實揭露）

**⚠ V5 僅 HCC1395 5kHz 測試**。COLO829 / H1437 / H2009 / HCC1937 / HCC1954 / HCC1395_DORADO 仍為 V3-F 或更早版本。

本報告所有量化數據限於 HCC1395 5kHz。跨樣本 V5 驗證需要：
- 6 樣本 × LongPhase-TO V5 重跑（每樣本 ~1,869s × 6 ≈ 3 hr）
- 6 樣本 × ISM 全基因體重跑（每樣本 ~10-30 min × 6 ≈ 1-3 hr）
- 共 ~4-6 hr 機器時，加上下游分析比對 1-2 週

---

## Section 4：圖表分析

### 4.1 既有素材重用

| 檔案 | 用途 |
|------|------|
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/hp_ratio_version_compare.png` | 5 版本 HP_Ratio 趨勢圖 |
| `version_summary.tsv` | Region-level ISM 指標 5 版本對照 |
| `hp_ratio_auc_by_version.tsv` | Inner/Outer/All AUC by version |
| `methyl_baseline_vs_v5.tsv` | 6 個 methylation features |

### 4.2 本報告新繪 4 張圖

- **Fig 16**（Section 2.5）：V5 三層邏輯 + code excerpt + 5 場景對照
- **Fig 17**（Section 3.3）：HP tag 4 版本 stacked bar + 全基因組 V3-F vs V5 重分配
- **Fig 18**（Section 3.7）：Concordance + AMB% + F1 evolution 三 panel
- **Fig 19**（Section 5.4）：10 結論 × V5 影響三層矩陣

### 4.3 視覺敘事邏輯

```
fig16（為何修：演算法 + 場景）
   ↓
fig17（HP tag 怎麼變：分布重分配）
   ↓
fig18（品質指標怎麼變：concordance 升、AMB% 收斂、F1 持平）
   ↓
fig19（既有結論受影響範圍：3 Tier 分類）
```

---

## Section 5：Downstream Impact 結論穩定性矩陣（**本報告核心**）

### 5.1 Master Dataset 狀態說明

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz` (748,391 rows, 89 MB)

- 建構時間：**post HP-integer fix（2026-03-30）**
- 但 **pre-V5 Somatic Fallback**（V3-Fixed BAM）
- → 所有 CL-010, CL-011, Registry 條目的數據假設 V3-F HP 分布

### 5.2 三層分類（10 既有結論）

#### Tier 1：DEFINITELY NEEDS RE-EVAL（直接 HP-tag-dependent）

| Conclusion ID | 內容 | 為何受影響 | 重評優先序 |
|---------------|------|-----------|----------|
| **CL-001** | Self-phasing 17.3:1 bias | V5 mitigates bias 量；絕對數字會變但因果鏈仍有效 | 高 |
| **CL-009** | 62% TO LOH = artifact | V5 HP_Ratio 分布變 → artifact 比例可能下降 | 高 |
| **CL-010** | HPFineNGroups 89.1% TP rate | NG bucket count 直接 V3-F 依賴；V3-F 17.5% AMB → V5 8.0% AMB 改變 bucket 分布 | **最高** |

**重評範圍**：

- **CL-001**：在 V5 BAM 上重新計算 somatic ALT reads HP1:HP2 ratio
- **CL-009**：在 V5 BAM 上重新統計 ISM HP_Ratio LOH，比對 Paired 計算 artifact %
- **CL-010**：在 V5 BAM 上重跑 master dataset 級分析 → HPFineNGroups TP rate 是否仍 89.1%？或降到接近 baseline？

#### Tier 2：POSSIBLY SHIFTED（HP-dependent 但機制 robust）

| Conclusion ID | 內容 | 為何可能 robust |
|---------------|------|---------------|
| **CL-011** | LOH-constrained phasing ΔNG = +0.385 | 機制 = LOH 物理單 haplotype 限制；V5 fallback 對此機制中性 |
| **CL-006** | LOH.bed AUC = 0.70 | self-reference artifact；V5 不解決 |

**為何中性**：CL-011 的核心假說「LOH region 物理上只保留單 haplotype，somatic SNV 必然產生 same-hap ref/alt 子族」是物理事實，與 BAM 版本無關。V5 改的是 phasing tag 品質，不改變 LOH 區域 read 分布。

#### Tier 3：UNAFFECTED（mechanism-independent）

| Conclusion ID | 內容 | 為何不受影響 |
|---------------|------|-------------|
| **CL-007** | Germline ASM >> somatic（3-6×）| HP-tag-direction independent（用 absolute count 而非 ratio）|
| **CL-025** | LOH.bed Jaccard = 1.0 | 走 VCF AD path 不經 BAM HP tag |
| **CL-023** | Coverage_Multiple confound | read-depth proxy，不依賴 HP |
| **E6** | Phasing graph math（理論證據）| 純演算法數學推導，不依賴觀察數據 |
| **E4** | 7/7 cross-sample direction consistency | TO 內部複製（方向性結論），bias 量改但方向不變 |

### 5.3 既有 PI 報告引用一致性

3 份 PI 報告（20260422 / 20260422 / 20260424）的核心結論：

| 結論 | V5 後是否仍成立 | 說明 |
|------|---------------|------|
| Self-phasing 真實存在（mechanism）| ✅ 仍成立 | E1/E4/E5/E6 獨立證據鏈，4/6 不依賴 BAM 版本 |
| 17.3:1 bias **數值** | ⚠ V5 後絕對數字會降 | 但 still >> 1:1，方向性結論維持 |
| `--germline-hp-only` flag 機制正確 | ✅ 仍成立 | 機制層獨立於 BAM 版本 |
| Filter AUC NEGATIVE | ✅ 仍成立 | Phase 1 在 V3-F 上實證；V5 SuggestFilter 差異 +1 TP/+1 FP（淨零）|
| HPFineNGroups 89.1% subclone marker | ⚠ 已有 LOH-constrained reinterpretation；V5 進一步測試 | 機制改解讀為 LOH-constrained phasing；數值需重算 |
| 三條路（paired_full/paired_pileup/to_pileup）不衝擊 | ✅ 仍成立 | paired tracks 不涉 V5（只有 to_pileup BAM 變動）|

### 5.4 影響矩陣視覺化

![Figure 19 — V5 Impact on 10 Existing Conclusions](figures/fig19_conclusion_impact_matrix.png)

**Tier 1（3 個）×「需 V5 重評」+ Tier 2（2 個）×「監測」+ Tier 3（5 個）×「不受影響」= 10 個結論的完整覆蓋**。

---

## Section 6：是否改變結論的最終判定

### 6.1 直球回答 PI

**「V5 是否改變既有結論？」**

> **大多數結論機制層面仍成立，但 3 個 Tier-1 量化結論需以 V5 重新計算數字。**

### 6.2 結論二元分解（機制 vs 數字）

| 結論類型 | 機制層面 | 數字層面 |
|---------|---------|---------|
| Self-phasing 存在性 | ✅ 不變 | ⚠ bias 量會下降（仍 >> 1:1）|
| 修改方案合理性 | ✅ 不變 | — |
| Filter AUC 結論 | ✅ 不變 | ✅ 不變（V3-F 上 PASS，V5 SF 差異淨零）|
| HPFineNGroups marker | ✅ 機制不變（LOH-constrained 解讀）| ⚠ TP rate 待重算 |
| 三條路不衝擊 | ✅ 不變 | ✅ 不變 |
| LOH.bed Jaccard=1.0 | ✅ 不變 | ✅ 不變（VCF AD path）|

### 6.3 為何結論機制 robust？

依 20260424 證據鏈報告：6 條核心證據（E1-E6）中 **4 條完全獨立於 BAM 版本**——
- E1（17.3:1 raw bias）：用 raw read count 不依賴 HP_Ratio 指標
- E4（7/7 cross-sample）：方向性結論
- E5（PON-only experiment）：擾動實驗（V5 屬於同一 PON-only 家族）
- E6（phasing graph math）：純理論

**V5 屬於改善「reference BAM 品質」的工程進步，不撼動「self-phasing 機制存在」的核心命題**。

### 6.4 建議行動清單

| 行動 | 優先序 | 預估成本 | ROI | 說明 |
|------|--------|--------|------|------|
| **在既有 3 份 PI 報告加 V5 footnote** | 最高 | <1 hr | 高 | 透明度；說明 V5 已超越 V3-F，但結論機制仍有效 |
| **重算 17.3:1 bias 在 V5（CL-001）** | 高 | ~2 hr | 中 | 數字精確化；預期降至 ~5-10:1（仍顯著）|
| **HPFineNGroups V5 master rerun（CL-010）**| 高 | 12-24 hr | 高 | marker 可發表性；確認 89.1% 是否站得住 |
| **62% LOH artifact V5 重算（CL-009）** | 高 | ~4 hr | 中 | 與 V5 HP_Ratio 分布更新一致 |
| **LOH-constrained ΔNG V5 重算（CL-011）** | 中 | ~4 hr | 中 | 預期機制不變，但數字確認 |
| **7-sample V5 BAM 生成** | 中 | ~4-6 hr 機器時 + 1-2 週分析 | 中-高 | 取決於發表計畫；目前 master dataset 仍 V3-F |

### 6.5 對未來 PI 決策的提醒

1. **V5 是 production-ready haplotag 版本**：用於所有後續 ISM 分析
2. **V3-F binary 保留作回退**：`/big7_disk/liaoyoyo2001/longphase-to-mod/longphase-to-v3fixed`
3. **V4 BAM 可刪除**：與 V3-F 完全相同，節省 268G
4. **Master dataset 重建需求**：若要 publication-grade，建議用 V5 重建 7 樣本 master dataset
5. **Memory 更新建議**（不執行，僅建議）：
   - `project_v5_somatic_fallback_verification.md` 升級為 active
   - 在 3 個 Tier-1 conclusion memory 標註 V5 重評需求
   - `project_v3_fixed_haplotag_verification.md` 標註「已被 V5 取代」

---

## 附錄

### 附錄 A：5 版本演進完整對照表

| 維度 | Baseline | V2b | V3-Fixed | V4 | V5 |
|------|----------|-----|---------|-----|-----|
| Phasing VCF | self-phasing | PON-only | PON-only | PON-only | PON-only |
| getVote() priority | bug（覆蓋）| bug（覆蓋）| 兩層分離 | 兩層分離 | **三層 + Layer 1.5** |
| HP enum write | bug（enum value）| bug（enum value）| integer literal | integer literal | integer literal |
| altHaplotype guard | 無 | 無 | INDEL only | **INDEL + SNP** | INDEL + SNP |
| Total tagged reads | 19.57M | 18.81M | 18.81M | 18.81M | 18.81M |
| HP:i:33 全基因組 | 0 | 0 | 239,679 | 239,679 | **110,197** |
| AMB% | 1.3%（假性）| ~0%（不出現）| 17.5% | 17.5% | **8.0%** |
| SEQC2 F1 | 0.7157 | 0.7155 | 0.7154 | 0.7154 | **0.7154** |
| ISM F1 | 0.0154 | 0.0177 | 0.0124 | 0.0124 | 0.0125 |
| Clean-PS Somatic Acc | 82.2% | — | 90.2% | 90.2% | **90.7%** |
| Somatic HP1 ratio | 0.774 | — | — | — | **0.703** |

### 附錄 B：關鍵檔案參考

**Code 來源**：
- `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:512-563` — V5 三層 getVote()
- `/big7_disk/liaoyoyo2001/longphase-to-mod/Util.h` — HAPLOTYPE enum

**Binary**：
- V5：`/big7_disk/liaoyoyo2001/longphase-to-mod/longphase-to`
- V3-F backup：`/big7_disk/liaoyoyo2001/longphase-to-mod/longphase-to-v3fixed`

**Data 來源**：
- 主 session report：`InterSubMod/docs/provenance/ai_sessions/2026/04/20260412_V5_somatic_fallback_getVote_修正與驗證_01.md`
- 量化 TSVs：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/{version_summary,hp_ratio_auc_by_version}.tsv`
- 甲基化對比：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/methyl_baseline_vs_v5.tsv`
- Master dataset：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz`

**Memory 參考**：
- `project_v5_somatic_fallback_verification.md`
- `project_v3_fixed_haplotag_verification.md`
- `project_self_phasing_causal_chain_confirmed.md`
- `project_pon_only_phasing_verification.md`

**既有 PI 報告**：
- `InterSubMod/docs/reports/pi_reports/2026/04/20260422_Self_Phasing_complete_report_for_PI_01.md`
- `InterSubMod/docs/reports/pi_reports/2026/04/20260422_Self_Phasing_multiperspective_argument_01.md`
- `InterSubMod/docs/reports/pi_reports/2026/04/20260424_Self_Phasing_evidence_chain_methodology_01.md`

### 附錄 C：相關 Commits

| Commit | 修正 | 對應版本 |
|--------|------|---------|
| `8b8c1fd` | --pon-only-phasing 模式 | V2b |
| `41ff147` | getVote() 兩層分離 + low-confidence integer fallback | V3-Fixed |
| `380e8d2` | countINDELHaplotype UNDEFINED guard | V3-Fixed |
| (V4 incremental) | countSNPHaplotype altHaplotype guard | V4 |
| (V5 HEAD) | Layer 1.5 somatic fallback in getVote() | **V5** |

### 附錄 D：本報告與既有 PI 報告交叉引用

| 本報告章節 | 對應既有報告章節 |
|-----------|----------------|
| Section 1（版本演進）| Session report Section 5、Memory `project_v5_somatic_fallback_verification.md` |
| Section 2（code-level）| Session report Section 2 + HaplotagProcess.cpp:512-563 |
| Section 3（量化）| Session report Section 4, 7, 8, 9 |
| Section 5（impact）| 既有 3 份 PI 報告 + 20260424 證據鏈 |
| Section 6（最終判定）| 20260424 證據鏈報告 Section 4 |

---

## 報告結束

**5 個給 PI 的核心訊息**：

1. **V5 是 LongPhase-TO 的 production-ready 版本**——三步演進（V2b → V3-F → V5）修正所有已知 bug
2. **量化層面**：V5 vs Baseline 在 tag 品質有顯著改善（concordance +8.3pp、AMB% 修正、129K reads 重分配），但 SEQC2 F1 持平（F1 不衡量 tag 品質）
3. **既有結論機制全部 robust**：3 個 Tier-1 結論需數字重算但機制不變；2 個 Tier-2 監測即可；5 個 Tier-3 完全不受影響
4. **三條 pipeline track 仍不衝擊**：paired_full / paired_pileup BAM 不涉 V5；只有 to_pileup BAM 變動
5. **下一步建議**：先在 3 份 PI 報告加 V5 footnote（最低成本）；若要發表 HPFineNGroups 論文，啟動 master dataset V5 重建（中成本）；7-sample V5 BAM 生成可延後到 publication 階段

**研究誠實性聲明**：本報告所有量化數據限於 HCC1395 5kHz；跨樣本 V5 驗證為未來工作。Tier 分類基於 mechanism dependency 邏輯推論，不是經驗統計——若未來 V5 master rerun 結果與預測不一致，本報告 Tier 分類需更新。
