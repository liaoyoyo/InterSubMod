<!--
建立時間: 2026-04-28 00:00
目標: 審核 Self-phasing 完整問題、最新結論、GT/HP/PS 細節處理，以及 baseline 與 V5 版本判斷
處理範圍: 現有文件、Knowledge Base、程式碼定義、TSV 證據與小型單元測試；不做 LongPhase/ISM 大型重跑
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md
  - /big8_disk/liaoyoyo2001/knowledge/03_file_formats/bam-format.md
  - /big8_disk/liaoyoyo2001/knowledge/03_file_formats/vcf-longphase.md
  - /big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md
  - /big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing-workflow.md
-->

# Self-phasing / Baseline / V5 審核報告

## 一句結論

**目前 Self-phasing 的最新判斷是成立且已被修正到較精確版本：TO baseline 的主要問題不是 LOH.bed 邊界錯，也不是 methylation bimodality，而是 LongPhase-TO 在 tumor-only 條件下把 somatic / germline-like signals 放進 phasing/tagging 路徑後，於 BAM HP tag 產生自參考偏斜；Thread D 已把 HPFineNGroups 正確重詮釋為 phasing bucket signature。V5 對 haplotag 行為的修補通過 sanity 與 paired concordance 檢查，可作為新的 tag-level baseline，但它不等於 caller-level F1 改善證據，也不能解除 HPFineNGroups marker 需重驗的風險。**

本次審核只做證據審查，不做大型重跑。

---

## 1. 審核範圍與判決

| 問題 | 判決 | 重點 |
|------|------|------|
| Self-phasing 問題是否真實 | **已證明成立** | baseline 17.3:1 somatic HP1:HP2 bias；TO HP_Ratio 與 paired 幾乎不相關；PON-only / V5 可消除大部分 tag bias |
| LOH.bed 是否被 self-phasing 污染 | **目前判斷：否** | LOH.bed baseline vs PON-only Jaccard=1.0000；正式 LOH 需看 phased genotype / LOH.bed / CNV，不可只看 BAM HP tag |
| HPFineNGroups 是否代表 methylation bimodality | **舊說法錯誤，已修正** | C++ 定義是 4 bucket occupancy：`HP1 / HP1-1 / HP2 / HP2-1`；與 methylation 無直接計算關係 |
| Thread D LOH-constrained phasing signature | **成立，grade B** | 6/6 TO gap 同方向；paired negative control gap 消失；flag-on 使 same-hap bucket 物理塌陷 |
| `--germline-hp-only` 是否可當 filter 修復 | **不建議** | 機制正確，但 HCC1395 TO filter gate FAIL；主要價值是機制診斷與 negative control |
| baseline vs V5 | **V5 tag 修補可信，但 caveat 必須保留** | sanity 15/15 PASS、clean PS paired concordance +13.3pp；但 V5 binary 仍含未 commit working-tree 修改，且不代表 F1 改善 |

---

## 2. 必要背景與術語邊界

根據 Knowledge Base：

- [bam-format.md](/big8_disk/liaoyoyo2001/knowledge/03_file_formats/bam-format.md)：LongPhase-TO 的 BAM `HP:i:1/2/11/21/33` 是 read-level haplotag；`11/21/33` 對應 `1-1/2-1/3` 語意，但為整數格式。
- [vcf-longphase.md](/big8_disk/liaoyoyo2001/knowledge/03_file_formats/vcf-longphase.md)：VCF `GT` 使用 `|` 表示 phased genotype，`PS` 表示 phase set；相同 PS 屬於同一 phase block。
- [longphase-to.md](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md)：LongPhase-TO `phase` 輸出 phased VCF 與 `LOH.bed`，`haplotag` 才輸出 BAM HP tag；正式 LOH 不由單一 HP tag 定義。
- [phasing-workflow.md](/big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing-workflow.md)：LongPhase-S 的 `HP=3` 是 read-level ALT tag，不等於 LOH；LongPhase-TO 的 LOH 來自 phased VCF / `LOH.bed`。

### 2.1 三層資料不可混用

| 層級 | 檔案/欄位 | 正確角色 | 常見誤用 |
|------|-----------|----------|----------|
| Caller | ClairS-TO VCF `FILTER/AF/GQ/DP` | 決定 PASS candidate 與 caller-level benchmark | 把 LongPhase phase 寫成重新 caller |
| Phasing | phased VCF `GT/PS/GT2/GT3`、`LOH.bed` | 決定 phase block、sub-genotype、region-level LOH | 把 HP tag skew 寫成 LOH.bed |
| Haplotag | BAM `HP:i:1/2/11/21/33`、`PS` | read-level haplotype assignment | 把 `HP:i:21` 直接等同該位點 ALT |

---

## 3. Self-phasing 完整問題鏈

### 3.1 問題來源

baseline LongPhase-TO 在 tumor-only 條件下沒有 matched normal；PON 可降低 germline leak，但不能完整排除所有 germline-like het。當 somatic 或 germline-like variants 進入 phasing / haplotag 決策時，read-level HP assignment 會被該 variant 自己與附近 variants 的共現 pattern 牽引，產生自參考。

歷史文件 [02_Self_Phasing根因.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md) 的核心數字仍可作為背景：

| 指標 | 數值 | 審核判斷 |
|------|------|----------|
| Somatic HP1:HP2 bias | 17.3:1 | self-phasing artifact 的核心證據 |
| TO TP ISM HP_Ratio LOH 中消失比例 | 62% | 僅限 **ISM HP_Ratio LOH**，不可外推為 LOH.bed |
| 全 TO ISM HP_Ratio LOH 中 artifact | 31.2% | 同上，屬 BAM HP tag path |
| LOH.bed baseline vs PON-only | Jaccard=1.0000 | 表示 region-level LOH.bed 不受此 artifact 改變 |

### 3.2 已修正的歷史誤判

早期文件曾把「self-phasing 造成 TO LOH」寫得過寬。最新正確表述應拆成兩句：

1. **BAM/ISM HP_Ratio LOH 會被 self-phasing 影響**，因此 HP-dependent features 在 baseline TO 下需降級或重驗。
2. **LongPhase-TO `LOH.bed` region boundary 不因 PON-only 修正而改變**，目前證據是 Jaccard=1.0000；LOH.bed 可作 region-level reference，但仍需 CNV / SEQC2 背景確認。

---

## 4. Thread D 最新結論審核

最新主軸文件 [20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md) 已把 4/23 週報中的待驗事項補成四層證據鏈。

### 4.1 HPFineNGroups 的正確定義

程式碼審核：

- [LabelTest.cpp](/big7_disk/liaoyoyo2001/InterSubMod/src/core/LabelTest.cpp) line 265：`hp_to_fine_labels()` 只依 `hp_tag` 分成 4 類。
- [Stats.hpp](/big7_disk/liaoyoyo2001/InterSubMod/include/core/Stats.hpp) line 323：`d_hp1_hp1s` 註解明確是 "Same haplotype, germline vs somatic"。

| Bucket | ISM 內部值 | 意義 |
|--------|-----------:|------|
| HP1 | 0 | germline haplotype 1 |
| HP1-1 | 1 | somatic-on-HP1 |
| HP2 | 2 | germline haplotype 2 |
| HP2-1 | 3 | somatic-on-HP2 |

**審核判決**：HPFineNGroups 不是 methylation bimodality，也不是 methylation cluster count。它是 phasing / variant-presence 的 bucket occupancy。舊的 "Haplotype-loss-dependent methylation bimodality" TO 主軸應撤回，改為 "LOH-constrained phasing signatures"。

### 4.2 四層證據鏈

| Layer | 證據 | 關鍵結果 | 判斷 |
|------|------|----------|------|
| X5 KDE-corrected TO | 6 TO 樣本 Inner same_HP1 vs Outer cross_het | 6/6 gap 正向，Wilcoxon exact p=0.0156 | 支持跨樣本 TO signature |
| B1 pre-KDE master | 舊 master 版本重算 | median gap 約 +0.37，與 X5 ordering 一致 | 排除 KDE artifact |
| B3 paired negative control | 7 paired_full 樣本 | paired gap median 約 0，p=0.578 | matched normal 排除 germline het 後 gap 消失 |
| X3 flag-on collapse | HCC1395 TO `--germline-hp-only` | NG>=3 30,880 -> 0；same_HP1 n=219 -> 0 | 直接證明 somatic HP attribution 支撐 bucket |

**審核判決**：Thread D 的主要機制判斷成立，但 evidence grade 應維持 B 而非 A，因為 n=6/7 仍小，且 R-SELFREF 尚未用 Phase 2B master + flag=on 完全解除。

---

## 5. 錯誤分析與處理狀態

| 舊說法 / 風險 | 現在狀態 | 處理是否充分 | 建議寫法 |
|---------------|----------|--------------|----------|
| HPFineNGroups 是 methylation bimodality | **錯誤，已更正** | 充分；已有 C++ 定義與 Thread D 主軸 | 寫成 phasing bucket occupancy / LOH-constrained phasing signature |
| Thread B cross-sample whitelist filter 可用 | **撤回** | 充分；已有撤回宣告 | 只保留 HCC1395 case study / characterization |
| `--germline-hp-only` 可修掉 filter AUC | **否** | 充分；HCC1395 TO ΔAUC gate FAIL | 寫成 negative control / diagnostic flag |
| LOH.bed 可能被 self-phasing 循環污染 | **目前排除** | 充分到 HCC1395 層；跨樣本仍可補強 | LOH.bed 不受 BAM HP self-phasing 影響，但不是 CNV truth 的替代 |
| V5 全面優於 baseline | **過寬** | 已在 V5 suite 修正 | V5 在 sanity 與 clean PS paired concordance 勝出；problem PS / weak directional 位點仍有限 |
| LongPhase phase 改善 caller F1 | **不支持** | 已由 V5 suite 釐清 | phase/tag 改變 HP 品質；caller-level F1 由 ClairS-TO PASS set 決定 |
| `HP:i:21` 必然是當前位點 ALT | **錯誤** | 已列入 code issue inventory | HP tag 是 read-level phasing state；需用 per-site ALT/REF 檢查 |

---

## 6. GT / Phasing / Tagging 細節審核

### 6.1 VCF GT 與 BAM HP 的角色不同

根據 [12_gt_distribution_audit.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md)，baseline vs V2b/V5 的 PASS somatic GT class 幾乎一致：

| GT class | Baseline % | V2b % | Δpp | 審核判斷 |
|----------|-----------:|------:|----:|----------|
| Germline_Het | 13.09 | 12.91 | -0.18 | 小幅差異 |
| Germline_Hom_or_LOH | 7.48 | 7.48 | 0.00 | 一致 |
| Somatic_NoLOH | 42.99 | 42.99 | 0.00 | 一致 |
| Somatic_in_LOH | 24.19 | 24.19 | 0.00 | 一致 |
| Unphased | 12.26 | 12.44 | +0.18 | 小幅差異 |

因此 17.3:1 HP skew 不能寫成「PASS somatic GT 被大量改判」；更準確地說，它主要發生在 haplotag / read assignment 層。

### 6.2 V5 tag 修補機制

根據 [13_phase_vs_tag_algorithm_detail.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md) 與 [01_code_diff_analysis.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md)：

| 版本 | 重點 | 影響 |
|------|------|------|
| baseline | 原始 LongPhase-TO，somatic / unknown 可污染 phasing/tagging | 17.3:1 somatic HP1:HP2 bias |
| V2b (`8b8c1fd`) | 加 `--pon-only-phasing`，但未改 `HaplotagProcess.cpp` | 暴露 getVote priority bug |
| V3-Fixed (`41ff147`) | two-layer `getVote()`：germline first, somatic second | 修正大規模 HP21 偏移 |
| `380e8d2` | INDEL guard | 消除 undefined array access 風險 |
| V5 working tree | Layer 1.5 somatic fallback + SNP guard | germline-poor reads 得到更合理 HP11/HP21 assignment |

**審核判決**：V5 的修補是 tag-quality fix，不是 caller F1 fix。報告或簡報中應避免把 "V5 更可信" 寫成 "V5 讓 somatic calling 更準"。

---

## 7. baseline vs V5 量化對照

### 7.1 ISM aggregate 對照

來源：[version_summary.tsv](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/version_summary.tsv)

| Version | n_regions | TP | FP | TP_rate | Potential_LOH% | HP_Ratio median | HPFineNGroups=2% | HPFineNGroups>=3% |
|---------|----------:|---:|---:|--------:|---------------:|----------------:|------------------:|-------------------:|
| baseline_old | 40,213 | 28,383 | 11,830 | 0.706 | 58.7 | 0.788 | 51.2 | 27.9 |
| v5_new | 40,096 | 28,495 | 11,601 | 0.711 | 62.2 | 0.574 | 49.7 | 29.8 |

解讀：

- TP_rate 差異很小，不能宣稱 caller-level F1 改善。
- HP_Ratio median 從 0.788 拉回 0.574，符合 tag bias 被修正的方向。
- HPFineNGroups 分布沒有被簡單清空；因此 V5 與 `--germline-hp-only` 是不同層次的修補，不能混用結論。

### 7.2 HP_Ratio AUC 對照

來源：[hp_ratio_auc_by_version.tsv](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/hp_ratio_auc_by_version.tsv)

| Version | Side | n | TP_rate | HP_Ratio_AUC |
|---------|------|--:|--------:|-------------:|
| baseline_old | All | 40,213 | 0.706 | 0.531 |
| v5_new | All | 40,096 | 0.711 | 0.526 |
| baseline_old | Inner | 23,620 | 0.729 | 0.523 |
| v5_new | Inner | 24,949 | 0.730 | 0.525 |
| baseline_old | Outer | 16,593 | 0.673 | 0.502 |
| v5_new | Outer | 15,147 | 0.679 | 0.510 |

解讀：V5 修正 tag distribution，但 HP_Ratio 本身仍不是強 filter feature；AUC 仍接近隨機。

### 7.3 methylation feature 對照

來源：[methyl_baseline_vs_v5.tsv](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/methyl_baseline_vs_v5.tsv)

| Feature | baseline AUC_All | v5 AUC_All | 審核判斷 |
|---------|-----------------:|-----------:|----------|
| HPMergedDelta | 0.517 | 0.514 | 不變 |
| AlleleDelta | 0.541 | 0.543 | 不變 |
| HPFineF | 0.576 | 0.584 | 小幅變動，仍未成為可靠 filter |
| PairwiseMeanDist | 0.525 | 0.528 | 不變 |
| CramersV | 0.507 | 0.506 | 不變 |
| GlobalP | 0.533 | 0.530 | 不變 |

**審核判決**：V5 沒有推翻 pure methylation feature 空間耗盡的結論。

---

## 8. InterSubMod ReadParser 與 `--germline-hp-only`

InterSubMod 端的邏輯在 [ReadParser.cpp](/big7_disk/liaoyoyo2001/InterSubMod/src/core/ReadParser.cpp) line 120：

- raw `HP:i:1` -> `"1"`
- raw `HP:i:2` -> `"2"`
- raw `HP:i:11` -> `"1-1"`
- raw `HP:i:21` -> `"2-1"`
- raw `HP:i:33` -> `"3"`
- `--germline-hp-only` 開啟時，`"1-1" / "2-1" / "3"` 會被 demote 為 `"0"`，但 `hp_tag_raw` 保留供 audit。

測試覆蓋：

- [test_read_parser.cpp](/big7_disk/liaoyoyo2001/InterSubMod/tests/test_read_parser.cpp) line 83：HP tag parsing 與 demotion。
- [test_global_local.cpp](/big7_disk/liaoyoyo2001/InterSubMod/tests/test_global_local.cpp) line 404：HPFine label mapping。
- 本次已跑小型測試：`./build/bin/run_tests --gtest_filter=ReadParserHPTagTest.*:GlobalTestTest.HPFine*:FullLabelTest.HPFineLabel_Mapping:SignificanceAnalyzerTest.HPFamily_GatingPassesWhenPureFails`，17/17 pass。

**重要邊界**：`--germline-hp-only` 是下游 ISM demotion，V5 是 upstream LongPhase-TO phase/tag 修補；兩者不是同一個 fix。

---

## 9. 仍需保留的 caveat

| Caveat | 影響 | 建議 |
|--------|------|------|
| V5 binary 來自 `380e8d2` + 未 commit working-tree 修改 | 可追溯性不足 | 將 Layer 1.5 fallback 與 SNP guard 切獨立 commit |
| confidence threshold 0.6 尚未直接以 vote log 驗證 | V5 設計證據仍有一格間接 | 加 vote log 或 debug TSV |
| `HP:i:21` 可包含當前位點 REF read | HPFineNGroups bucket 可能混合 ALT/REF | 未來加 8-bucket derived metric 或 per-site ALT/REF audit |
| 0.6 purity simulation 尚未執行 | V5 低 purity generalization 還是推論 | 按既有 `20260427_purity06_simulation_plan_01.md` 另案執行 |
| HPFineNGroups marker 仍有 R-SELFREF | grade B -> A 阻塞 | Phase 2B master + flag=on/off 7 樣本重驗 |
| paired ground truth 本身也由 phasing tool 產生 | paired concordance 不是絕對真值 | 用 clean/problem PS 分層並加外部 CN/SV 輔助 |

---

## 10. 最終建議口徑

### 可直接使用

1. **Self-phasing 是 baseline TO tag-level artifact**，會污染 HP-dependent ISM features。
2. **LOH.bed 不因 self-phasing 修正而改變**；LOH.bed 與 BAM HP_Ratio 是不同路徑。
3. **HPFineNGroups 應寫成 phasing bucket occupancy / LOH-constrained phasing signature**，不是 methylation bimodality。
4. **Thread D 是目前 TO 層主軸**，evidence grade B，已由 TO、paired negative control 與 flag-on collapse 支持。
5. **V5 可作為新的 haplotag baseline**，因為 sanity 全 pass 且 clean PS paired concordance 顯著改善。

### 不建議再使用

1. 不要寫「HPFineNGroups 是 methylation subclone marker」而不加 phasing caveat。
2. 不要把 `HP:i:11/21/33` 直接等於當前位點 ALT。
3. 不要把 V5 tag 修補寫成 caller F1 提升。
4. 不要把單一 HCC1395 whitelist filter 寫成跨樣本 production filter。
5. 不要把 `--germline-hp-only` 寫成正式修復；它是 diagnostic / negative control。

---

## 11. 參考文件索引

| 類別 | 文件 |
|------|------|
| 主軸結論 | [Thread D main axis](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md) |
| Thread B 撤回 | [Thread B retraction](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md) |
| V5 audit index | [V5 audit suite index](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md) |
| V5 synthesis | [08_synthesis_conclusions.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md) |
| GT audit | [12_gt_distribution_audit.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md) |
| phase/tag 細節 | [13_phase_vs_tag_algorithm_detail.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md) |
| Knowledge BAM | [bam-format.md](/big8_disk/liaoyoyo2001/knowledge/03_file_formats/bam-format.md) |
| Knowledge LongPhase VCF | [vcf-longphase.md](/big8_disk/liaoyoyo2001/knowledge/03_file_formats/vcf-longphase.md) |
| Knowledge LongPhase-TO | [longphase-to.md](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md) |
| Knowledge phasing workflow | [phasing-workflow.md](/big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing-workflow.md) |
