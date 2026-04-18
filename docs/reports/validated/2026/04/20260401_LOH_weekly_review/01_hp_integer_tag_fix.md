<!--
建立時間: 2026-04-01 12:00
目標: 完整記錄 HP Integer Tag Bug 的發現過程、修正方式、影響量化與文件修正分級，作為 LOH 週報審閱的第一章
處理範圍:
  - ReadParser.cpp HP tag 解析邏輯
  - LongPhase-TO 的 HP:i:11/21/33 整數格式
  - 7 個 TO 樣本修正前後的 LOH eligibility、Tier 分佈、hp_assign_rate 變化
  - 受影響欄位清單與文件修正分級
關聯檔案:
  - src/core/ReadParser.cpp (lines 121-142)
  - docs/reports/validated/2026/03/20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md
  - docs/reports/validated/2026/03/20260330_TO_HP_fix重跑驗證與validated文件修正分級_01.md
  - docs/reports/validated/2026/04/20260401_LOH_weekly_review/00_background.md
-->

# HP Integer Tag Bug 發現與修正

**文件定位**：本文件為 LOH 週報審閱系列的第一章，記錄本週最重要的數據品質事件——TO 模式下 HP tag 解析的系統性錯誤，以及修正後對所有 TO downstream 分析的影響。

**前置閱讀**：`00_background.md`（理解 HP tag、LOH-like、Tier 分層等術語）。

---

## 1. 問題定義

### 1.1 根本原因

ISM 的 HP tag 解析邏輯位於 `src/core/ReadParser.cpp`，負責從 BAM 檔案中讀取每條 read 的 haplotype 標籤。這個解析邏輯需要處理兩種不同的 HP tag 格式：

- **LongPhase-S（Paired 模式）**：寫入**字串格式**的 HP tag，例如 `HP:Z:1`、`HP:Z:2`。ISM 一直都能正確解析。
- **LongPhase-TO（Tumor-Only 模式）**：寫入**整數格式**的 HP tag，例如 `HP:i:11`、`HP:i:21`、`HP:i:33`。

問題在於：**舊版 ReadParser.cpp 對整數格式的 HP tag 只做了簡單的 `std::to_string()` 轉換**，把 `HP:i:11` 轉成字串 `"11"`、`HP:i:21` 轉成 `"21"`、`HP:i:33` 轉成 `"33"`。但是 ISM 的下游模組（LabelTest、HP family 統計等）只認識 `"1"`、`"2"`、`"1-1"`、`"2-1"`、`"3"` 這些 canonical 字串。因此，`"11"`、`"21"`、`"33"` 這些非標準字串在下游被歸類為**未知/未追蹤（untracked）**的 HP 值。

### 1.2 直接後果

由於 TO 模式下大部分 reads 帶有 `HP:i:11` 或 `HP:i:21` 標籤，這個 bug 導致：

- **30% 到 99% 的 TO reads 被歸類為 HP0（未分配）**，即使它們實際上已經被 LongPhase-TO 成功 phasing。
- ISM 認為這些 region 沒有（或幾乎沒有）HP 資訊，因此在計算 `HP_Ratio`、`effective_hp`、`Potential_LOH` 等欄位時產生了系統性偏差。
- 大量本應在 Tier A/A+（eff_hp >= 30）的 regions 被錯誤降級到 Tier C0（eff_hp = 0）或 Tier C（eff_hp 1-9）。
- 所有依賴 HP family 的統計指標（`VerificationClass`、`Quality_Score`、`agreement_type` 等）都受到汙染。

### 1.3 不受影響的範圍

以下分析**不受此 bug 影響**：

- **Paired 模式的所有分析**：LongPhase-S 使用字串格式 `HP:Z:1/2`，ISM 一直都能正確解析。
- **Somatic caller 的結果**：ClairS / ClairS-TO 的 variant calling、AF、GQ、DP、VCF benchmark 指標（precision / recall / F1）與此 bug 完全無關。
- **ISM 的甲基化原始數據**：`MM`/`ML` 標籤的解碼不受 HP tag 解析影響。
- **不依賴 HP 的 ISM 特徵**：`PairwiseMedianDist`、`AlleleDelta`、`CramersV` 等在修正前後完全一致。

---

## 2. 修正方式

### 2.1 修正後的程式碼

修正後的 `ReadParser.cpp`（lines 121-142）加入了 switch-case mapping，將 LongPhase-TO 的整數格式映射到 ISM 使用的 canonical 字串格式：

```cpp
// Extract HP tag (haplotype)
info.hp_tag = "0";  // Default: unknown/unphased
uint8_t* hp_aux = bam_aux_get(b, "HP");
if (hp_aux) {
    char type = hp_aux[0];
    if (type == 'Z' || type == 'H') {
        // String format (longphase-s): HP:Z:1, HP:Z:2, HP:Z:1-1, HP:Z:2-1, HP:Z:3
        info.hp_tag = bam_aux2Z(hp_aux);
    } else if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I') {
        // Integer format (longphase-to): 1=HP1, 2=HP2, 11=somatic-supported HP1, 21=somatic-supported HP2, 33=ambiguous somatic HP
        // Map to canonical string format used throughout LabelTest
        int hp_int = bam_aux2i(hp_aux);
        switch (hp_int) {
            case 1:  info.hp_tag = "1";   break;
            case 2:  info.hp_tag = "2";   break;
            case 11: info.hp_tag = "1-1"; break;
            case 21: info.hp_tag = "2-1"; break;
            case 33: info.hp_tag = "3";   break;
            default: info.hp_tag = std::to_string(hp_int); break;
        }
    }
}
```

### 2.2 每個 mapping 的語義解釋

LongPhase-TO 的整數 HP tag 編碼遵循特定的語義規則：

| LongPhase-TO 寫入 | ISM 映射結果 | 語義解釋 |
|:---|:---|:---|
| `HP:i:1` | `"1"` | Haplotype 1（標準 phasing 結果，與 Paired 模式的 `HP:Z:1` 等價） |
| `HP:i:2` | `"2"` | Haplotype 2（標準 phasing 結果，與 Paired 模式的 `HP:Z:2` 等價） |
| `HP:i:11` | `"1-1"` | Somatic-supported haplotype 1。read 攜帶 somatic ALT allele 且追溯到 germline haplotype 1。`11` 的編碼規則是：十位數 = haplotype（1），個位數 = somatic support indicator（1）。`"1-1"` 中第二個數字代表 somatic support，不是 phase block ID |
| `HP:i:21` | `"2-1"` | Somatic-supported haplotype 2。read 攜帶 somatic ALT allele 且追溯到 germline haplotype 2。`21` = haplotype 2 + somatic support indicator 1。`"2-1"` 中第二個數字同樣代表 somatic support，不是 phase block ID |
| `HP:i:33` | `"3"` | Ambiguous somatic haplotype。LongPhase-TO 無法確定該 read 的 somatic allele 歸屬於 haplotype 1 或 2 時使用此標籤。**這不是第三個 haplotype**，而是歸屬不明確的分類。映射為 `"3"` 後，ISM 將其排除於 HP1/HP2 family 統計之外 |

### 2.3 防護機制

修正中保留了 `default` 分支（`info.hp_tag = std::to_string(hp_int)`），確保未知的整數值不會導致程式崩潰，而是被轉為字串留待下游處理。同時，程式碼在進入 switch 之前先檢查了 BAM auxiliary data 的 type byte（`'c'`/`'C'`/`'s'`/`'S'`/`'i'`/`'I'`），覆蓋了 HTSlib 所有可能的有號/無號整數類型。

---

## 3. 影響量化

### 3.1 LOH Eligibility 變化（per-sample）

`LOH eligible` 定義為 `effective_hp >= 30`（Tier A 及以上），是進行有意義 LOH 分析的最低門檻。修正前後的變化顯示了 bug 的嚴重程度：

| 樣本 | 修正前 LOH eligible (TP) | 修正後 LOH eligible (TP) | 變化倍數 |
|:---|---:|---:|---:|
| COLO829 | 0.7% | 34.7% | 49.6x |
| H1437 | 24.5% | 95.5% | 3.9x |
| H2009 | 64.5% | 97.7% | 1.5x |
| HCC1395 | — | — | — |
| HCC1395_DORADO | 46.9% | 93.1% | 2.0x |
| HCC1937 | — | — | — |
| HCC1954 | — | — | — |

> 注：以 COLO829 為例，修正前僅 0.7% 的 TP regions 有足夠的 HP 資訊做 LOH 分析。修正後提升到 34.7%，代表舊版資料中 **99.3% 的 LOH 判斷都建立在不充分的 HP 資訊之上**。

### 3.2 Tier 分佈變化

修正前後，TO 模式的 Tier 分佈發生了劇烈變化：

| Tier | 修正前佔比 | 修正後佔比 | 變化方向 |
|:---|---:|---:|:---|
| C0 (eff_hp = 0) | ~15% | ~0.2% | 大幅下降 |
| C (1-9) | ~17% | ~2.5% | 大幅下降 |
| B (10-29) | ~17% | ~9.5% | 下降 |
| A (30-49) | ~22% | ~16.3% | 略降 |
| A+ (>= 50) | 29.5% | 71.8% | **大幅上升** |

**核心結論**：修正後 **88% 的 TO regions 進入 Tier A/A+**（修正前僅 51%）。這表示 TO phasing（LongPhase-TO）的實際品質遠優於修正前的認知。先前對 TO phasing 品質的悲觀評估，很大程度上是 HP tag 解析錯誤造成的假象。

### 3.3 圖示：Tier 分佈修正前後對比

以下圖片位於 workspace 路徑，展示修正前後的 Tier 分佈變化：

**Tier 修正前後總覽**：

![Tier 分佈修正前後對比](../../../../../../output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation/figures/fig04_tier_before_after.png)

**Per-sample Tier 位移**：

![各樣本 Tier 位移](../../../../../../output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation/figures/fig06_per_sample_tier_shift.png)

> 圖片絕對路徑：
> - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation/figures/fig04_tier_before_after.png`
> - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_post_hp_fix_to_loh_investigation/figures/fig06_per_sample_tier_shift.png`

### 3.4 hp_assign_rate 與 agreement 重新量化

`hp_assign_rate` 衡量一個 region 中有多少比例的 reads 被成功分配到某個 haplotype。修正前後的變化極為顯著：

| 樣本 | 舊 hp_assign_rate 中位數 | 新中位數 | 舊 <0.05 比例 | 新 <0.05 比例 | 新 agreement_positive 比例 |
|:---|---:|---:|---:|---:|---:|
| HCC1395 | 0.5634 | 0.9863 | 4.27% | 0.24% | 58.63% |
| HCC1395_DORADO | 0.5309 | 1.0000 | 7.09% | 0.05% | 53.33% |
| COLO829 | 0.3659 | 1.0000 | 5.61% | 0.49% | 28.68% |
| H1437 | 0.3506 | 1.0000 | 9.60% | 0.18% | 60.47% |
| H2009 | 0.4737 | 1.0000 | 6.86% | 0.09% | 49.40% |
| HCC1937 | 0.7468 | 1.0000 | 2.84% | 0.12% | 55.65% |
| HCC1954 | 0.7097 | 1.0000 | 2.21% | 0.26% | 67.51% |

**解讀**：

- 修正前 `hp_assign_rate` 的中位數最低僅 0.35（H1437），表示超過一半的 reads 沒有被認出 HP tag。修正後中位數普遍達到 1.0，表示幾乎所有 reads 的 HP tag 都被正確識別。
- 修正前有 2-10% 的 regions `hp_assign_rate < 0.05`（幾乎無 HP），修正後降至 0.05-0.49%。
- `agreement_positive`（ISM 判斷與 HP label 一致的比例）修正後在 28-68% 之間，提供了可靠的 HP-methyl agreement baseline。

### 3.5 欄位翻轉比例（14 個 TP/FP 子集）

以 7 個 TO run 的 TP/FP `significance_summary.csv` before/after 直接比對，統計各欄位有多少 rows 的值發生了改變：

| 欄位 | 最小翻轉比例 | 中位翻轉比例 | 最大翻轉比例 |
|:---|---:|---:|---:|
| `HP_Ratio` | 39.26% | 75.20% | 95.40% |
| `Quality_Score` | 8.35% | 27.53% | 80.93% |
| `LOH_Subtype` | 8.18% | 27.24% | 79.64% |
| `Potential_LOH` | 7.72% | 26.92% | 79.54% |
| `Quality_Tier` | 2.28% | 8.47% | 79.12% |
| `HPFineSig` | 24.67% | 30.55% | 54.64% |
| `DominantLabel` | 9.17% | 27.85% | 47.66% |
| `HPMergedSig` | 8.12% | 19.03% | 41.65% |
| `VerificationClass` | 1.03% | 2.97% | 10.15% |
| `PairwiseMedianDist` | 0.00% | 0.00% | 0.00% |
| `CramersV` | 0.00% | 0.00% | 0.27% |
| `AlleleDelta` | 0.00% | 0.00% | 0.00% |

**解讀**：

- `HP_Ratio` 受影響最劇烈，中位翻轉率 75%——這是預期中的，因為 HP tag 是 `HP_Ratio` 的直接計算基礎。
- `Quality_Score`、`LOH_Subtype`、`Potential_LOH` 的中位翻轉率都在 27% 左右，代表約四分之一的 rows 數值已改變。
- `PairwiseMedianDist` 和 `AlleleDelta` 完全不受影響（0% 翻轉），因為它們不依賴 HP tag。這兩個特徵是目前最穩定、可跨版本比較的欄位。

---

## 4. 受影響欄位清單

### 4.1 直接受影響的欄位

以下欄位的計算過程直接或間接依賴 HP tag 的正確解析，在修正前的 TO 數據中**不可信**：

| 欄位 | 影響機制 |
|:---|:---|
| `HP_Ratio` | 直接由 HP1/HP2 read count 計算 |
| `effective_hp` (eff_hp) | 直接由具有 HP tag 的 reads 數量決定 |
| `Potential_LOH` | 基於 `HP_Ratio` 判斷 |
| `LOH_Subtype` | 基於 `HP_Ratio` 和 `effective_hp` 分類 |
| `VerificationClass` | 綜合 HP family 統計、label agreement 等判定 |
| `DominantLabel` | 依賴 HP family 的 read-label 分配 |
| `Quality_Score` | 包含 LOH penalty 和 HP 一致性加分 |
| `Quality_Tier` | 基於 `Quality_Score` 分層 |
| `hp_assign_rate` | 直接衡量 HP tag 被辨識的比例 |
| `agreement_type` | 依賴 HP label 與 ISM label 的比較 |
| `agreement_positive` | `agreement_type` 的二元化版本 |
| `HPMergedSig` | HP family merged 的顯著性指標 |
| `HPFineSig` | HP family fine-grained 的顯著性指標 |

### 4.2 不受影響的欄位

| 欄位 | 原因 |
|:---|:---|
| `PairwiseMedianDist` | 基於甲基化距離矩陣，不依賴 HP |
| `AlleleDelta` | 基於 ALT/REF allele 的甲基化差異，不依賴 HP |
| `CramersV` | 基於 label contingency table，極少數 rows 有 < 0.3% 翻轉 |
| Caller AF / GQ / DP | 來自 VCF，與 ISM 解析完全無關 |
| `MM`/`ML` 原始甲基化值 | BAM 中的甲基化標籤，不受 HP 解析影響 |
| TP / FP / FN benchmark 指標 | 來自 VCF benchmark，與 ISM 完全無關 |

---

## 5. 文件修正分級

修正完成後，所有 2026-03-30 之前基於 TO `significance_summary.csv` 的 HP/LOH 相關結論都不能直接沿用舊數字。根據影響程度，文件被分為三個等級：

### 5.1 A 級：必須出新版（不建議只加勘誤）

這些文件的核心結論直接依賴已失效的 TO LOH downstream 數字，需要基於修正後數據重新撰寫：

| 文件 | 原因 |
|:---|:---|
| `20260328_LOH_evidence_panel_final_report_01.md` | 整合所有 Round 的 TO feature 定位、Round 4 enrichment 已實質改變 |
| `20260329_LOH_round1_cross_sample_audit_v2_01.md` | 引用舊 Round 1 workspace，6 張 v2 圖全部改變 |
| `20260327_LOH_round2_support_hp0_analysis_01.md` | 混入局部修正，但仍大量連到舊 workspace、舊圖、舊 TO tier 表 |
| `20260327_LOH_round3_methyl_hp0_filter_01.md` | TO 計數、joint table、附錄 key numbers 都已改 |

### 5.2 B 級：主結論可保留，但內容與圖表需修正

| 文件 | 原因 |
|:---|:---|
| `20260330_LOH_round2_ps_export_and_to_block_audit_01.md` | Paired PS 邊界與 TO block-linkage 結論大致不變，但 top-block 表、少數圖與所有引用舊 Round 1 路徑需更新 |
| `20260327_LOH_round1_cross_sample_audit_01.md` | 已有 2026-03-30 勘誤聲明，可作歷史記錄保留，但需補上 post-fix 指向 |

### 5.3 C 級：可直接保留

| 文件 | 原因 |
|:---|:---|
| `20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md` | 本身就是 post-fix 影響評估文件 |
| `20260330_TO_LOH_enrichment_post_hp_fix_01.md` | 本身就是 post-fix TO LOH 正式結論 |

### 5.4 修正順序建議

1. **先出 Evidence Panel Final Report 的 post-fix 版**——這是對外最像「總結定稿」的文件
2. **再補 Round 1 Cross-Sample Audit v2 的 post-fix 版**——圖資已全部重跑完成
3. **接著重寫 Round 2 Support HP0 Analysis 與 Round 3 Methyl HP0 Filter**——整理 inline 勘誤為乾淨版本
4. **最後修 Round 2 PS-Block Audit**——更新 block-level `LOH-like_frac`、圖檔與 Round 1 path

---

## 6. 修正後的關鍵口徑變化

HP Integer Tag Fix 不只是一個技術修正，它改變了幾個重要的研究結論：

### 6.1 TO LOH 方向翻轉

| 口徑 | 修正前 | 修正後 |
|:---|:---|:---|
| TO LOH enrichment | 接近 1.0 或略 > 1.0（個別樣本） | 全面 < 1.0（0.85-0.96，TP 富集） |
| TO LOH 定位 | 可能是 FP filter 的候選特徵 | **不是 FP filter，而是 TP enrichment 的結構訊號** |
| TO Tier A(30-49) enrichment | 接近中性 | 0.706（TP 富集） |
| TO Tier A+(≥50) enrichment | 接近中性（~1.03） | 0.766（TP 富集） |

### 6.2 TO Phasing 品質評估

| 口徑 | 修正前 | 修正後 |
|:---|:---|:---|
| TO regions 在 Tier A/A+ 的比例 | ~51% | ~88% |
| TO hp_assign_rate 中位數 | 0.35-0.75（大量 reads 「未分配」） | 0.99-1.00（幾乎全部正確辨識） |
| 對 TO phasing 品質的評價 | 悲觀（大量 reads 似乎無法被 phase） | 樂觀（LongPhase-TO 的 phasing 品質遠優於先前認知） |

### 6.3 仍然成立的結論

以下結論在修正前後一致，可以直接保留：

1. **Caller benchmark 與 F1 分數**——完全不受影響
2. **`PairwiseMedianDist` 與 `AlleleDelta` 的特徵空間觀察**——0% 翻轉
3. **Paired 主線的所有 LOH 結論**——Paired 使用字串格式 HP tag，不受此 bug 影響
4. **HP0 方向性結論**：TO LOH-like 比 non-LOH-like 有更高 HP0 比例的方向仍成立，但效果幅度縮小，且數值基礎已改變

---

## 7. 已完成的重跑項目

修正後，所有需要重跑的 downstream workspace 已全部完成：

| 項目 | 狀態 | 新輸出路徑 |
|:---|:---|:---|
| Round 1 corrected rerun | 已完成 | `observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix` |
| Round 2 corrected rerun | 已完成 | `observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix` |
| Round 3 corrected rerun | 已完成 | `observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix` |
| Round 4 corrected rerun | 已完成 | `observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix` |
| PS-block audit corrected rerun | 已完成 | `research_rounds/20260330_loh_round2_ps_export_and_to_block_audit_post_to_hp_fix` |
| Round 1 v2 figure set | 已完成 | `docs/reports/validated/2026/03/assets/20260330_loh_round1_cross_sample_audit_v2_figures_post_to_hp_fix` |
| TO hp_fix reanalysis | 已完成 | `research_rounds/20260330_to_hp_fix_reanalysis/` |
| TO LOH corrections | 已完成 | `research_rounds/20260330_to_loh_corrections/` |
| Master dataset rebuild | 已完成 | `all_region_rows.tsv.gz` 已重建 |

> 以上 `observation_workspaces/` 和 `research_rounds/` 路徑的完整前綴為：
> `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/`

---

## 8. 最終判決

本次修正**沒有推翻 TO track 的全部研究**，但**明確推翻了所有依賴 TO HP family / LOH / label agreement 舊數字的結論**。

最合理的影響邊界是：

- **保留**：caller benchmark、AF、`PairwiseMedianDist`、`AlleleDelta`、Paired 主線
- **重算**：`VerificationClass` / `agreement` / `hp_assign` / `Quality_Score` / `LOH-like` 相關 TO 分析
- **改寫口徑**：TO LOH 是 TP enrichment，不是 FP filter

這不是「全部重做」，而是一次**邊界非常清楚的 TO downstream correction**。所有新基準數據已產出，可立即用於後續分析。

---

## 待驗證問題（已驗證 / 已更新）

### 已解決

3. **~~修正後 TO LOH TP enrichment 的機制~~** ✅ C13 確認：TO phasing 缺乏 normal reference → somatic allele 造成系統性 HP imbalance → TP LOH rate +51.6%（29.3%→44.5%），FP 幾乎不變。7/7 樣本一致。

4. **~~Tier enrichment 穩定性~~** ✅ C12 確認：6/7 樣本 A+ 實際優於 A（Simpson's Paradox）。Tier A 表觀優勢來自 hp_ratio 離散化放大。Per-sample 驗證完成，非少數樣本驅動。

### 尚未解決

1. **HP:i:33 的真實語義**：`HP:i:33` 映射為 `"3"` 是否最佳？需確認 LongPhase-TO 精細語義。定位 P2。
2. **Default 分支的觸發頻率**：是否有非預期整數值？定位 P2。
5. **HP0 filter 效果弱化**：修正後 +0.014 vs 修正前 +0.030。已知 HP0 filter 在 TO 方向相反（N3 否決），此問題優先級降為 P2。

---

## 認知門檻補充建議

1. **BAM auxiliary data 格式**：理解此 bug 需要知道 BAM 檔案中 auxiliary tags 有不同的 type encoding（`Z` = 字串、`i`/`I`/`c`/`C`/`s`/`S` = 各種整數）。同一個 tag name（`HP`）可以有完全不同的 data type，取決於上游工具的實作。
2. **LongPhase-S vs LongPhase-TO 的差異**：這兩個是同一個 phasing 工具（LongPhase）的不同模式，但輸出格式不同。TO 模式使用整數格式編碼 somatic support 資訊——`HP:i:11` 表示 somatic-supported haplotype 1（第二位數字是 somatic support indicator，不是 phase block ID）；`HP:i:33` 表示 ambiguous somatic haplotype（不是第三個 haplotype）。
3. **「翻轉比例」的意義**：欄位翻轉比例 75%（HP_Ratio）不代表 75% 的 rows「從 TP 變成 FP」。它代表 75% 的 rows 的 HP_Ratio 數值改變了，可能從 0.01 變成 0.45（從 LOH-like 變成 non-LOH），也可能只是微小的數值變動。影響的嚴重程度取決於這些變動是否跨越了判斷門檻（如 0.1 / 0.9）。
4. **Post-fix 不等於 ground truth**：修正後的數據比修正前更接近真實狀態，但 phasing 本身仍有錯誤（如 switch error）。HP_Ratio 接近 0 或 1 不一定代表真正的 LOH，也可能是 phasing 錯誤造成的假象。
5. **文件版本管理**：本專案的文件修正分級（A/B/C 級）是一個特定的管理框架。A 級不代表「內容全錯」，而是「核心數字基礎已改變，patch 修正的風險高於重寫」。理解這個分級邏輯有助於判斷哪些歷史文件的結論仍可引用。
