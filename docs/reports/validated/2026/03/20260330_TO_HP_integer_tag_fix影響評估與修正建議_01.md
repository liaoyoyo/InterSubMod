<!--
建立時間: 2026-03-30 10:30
目標: 整理 InterSubMod 修正 TO HP integer tag 後，對既有 TO 甲基分析、LOH 結論、candidate-specific rescue、Phase 1 衍生產物的影響範圍，並給出重跑與修正文案建議
處理範圍:
  - TO longphase-to / InterSubMod step05 輸出
  - 7 個 TO 樣本（HCC1395, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954）
  - 舊版 vs 修正版 significance_summary / method_design_validation / label_cluster_agreement 比較
關聯檔案:
  - src/core/ReadParser.cpp
  - docs/experiments/in_progress/2026/03/20260330_TO_HP_integer_tag_fix與全樣本重跑結果確認_01.md
  - docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_hp_fix_reanalysis/
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_corrections/
  - Knowledge/05_tools/longphase_to.md
  - Knowledge/06_workflows/phasing_workflow.md
-->

# TO HP Integer Tag Fix 影響評估與修正建議

**分析日期**：2026-03-30  
**適用範圍**：所有使用 `longphase-to` 後續 `InterSubMod step05` 輸出的 TO 分析  
**不適用範圍**：paired `longphase-s` 分析、caller benchmark 指標、原始 MM/ML 甲基讀值本身

---

## 1. 開頭結論

1. **需要重算的是 TO 的「HP/LOH/label 衍生結果」，不是原始甲基資料本身。**  
   `MM/ML` 讀值沒有改，`ClairS-TO` 呼叫結果也沒有改；真正變的是 `InterSubMod` 對 `HP:i:11/21/33` 的解析，進而影響 `VerificationClass / Potential_LOH / LOH_Subtype / Quality_Score / hp_assign_rate / agreement_type` 等衍生欄位。

2. **所有 2026-03-30 之前基於 TO `significance_summary.csv` 的 HP/LOH 結論，都不能直接沿用舊數字。**  
   尤其是任何依賴 `VerificationClass`、`agreement_positive`、`hp_assign_rate`、`Potential_LOH`、`Quality_Score` 的 TO 文稿與實驗表格。

3. **不需要重跑 caller 與 benchmark。**  
   `TP / FP / FN / precision / recall / F1` 來自 VCF benchmark，與此次 bug 無關；`PairwiseMedianDist` 與 `AlleleDelta` 也在 before/after 比較中維持不變。

4. **TO LOH 的舊口徑需要改寫。**  
   修正後 TO `LOH-like` 在 `eff_hp >= 30` 的區域呈現 **TP 富集**，不是先前一些文件暗示的 FP marker。  
   根據新版彙整，Tier A `FP/TP enrichment = 0.706x`，Tier A+ `= 0.766x`，方向已明確。

5. **建議的重跑邊界很清楚：只要從 `step05_intersubmod` 之後開始重建即可。**  
   不必重跑 `step01_clairs_to`、`step02_benchmark_clairs_to`、`step03_longphase_to`、`step04_benchmark_longphase_to`；需要重跑的是 `step05_intersubmod` 後所有依賴其 summary / region join / validation TSV 的下游分析。

---

## 2. Bug 根源與實際影響邊界

### 2.1 根因

根據 [Knowledge/05_tools/longphase_to.md] 與 [Knowledge/06_workflows/phasing_workflow.md]，TO phasing 走的是 `longphase-to` 路徑；這次問題出在 `InterSubMod` 對 TO BAM 中 `HP` integer tag 的解析。

`src/core/ReadParser.cpp:121-141` 的修正將：

| TO BAM 寫入 | 舊解析 | 新解析 |
|---|---|---|
| `HP:i:11` | `"11"` | `"1-1"` |
| `HP:i:21` | `"21"` | `"2-1"` |
| `HP:i:33` | `"33"` | `"3"` |

舊版直接把整數轉成字串，導致這些 reads 在 `LabelTest` 與後續 HP family 統計中無法被辨識，最後落入 `untracked`。

### 2.2 受影響與不受影響

**直接受影響**

- `VerificationClass`
- `DominantLabel`
- `Potential_LOH`
- `LOH_Subtype`
- `HP_Ratio`
- `Quality_Score`
- `Quality_Tier`
- `hp_assign_rate`
- `agreement_type`
- `agreement_positive`
- 所有依賴上述欄位的 candidate-specific rescue、support annotation、LOH audit

**不直接受影響**

- `TP / FP / FN / precision / recall / F1`
- `PairwiseMedianDist`
- `AlleleDelta`
- 幾乎全部 `CramersV`
- 原始 `MM/ML` 甲基資訊
- caller 端 AF / GQ / DP / VCF benchmark

---

## 3. 定量影響總覽

### 3.1 HP read 回收與 LOH eligibility

根據 `20260330_TO_HP_integer_tag_fix與全樣本重跑結果確認_01.md`：

- 修正前各 TO run 有 **29% 到 69%** 的 reads 落入 `untracked`
- 修正後 **所有樣本 untracked reads = 0**
- `LOH eligible (eff_hp >= 30)` 大幅提升  
  - `COLO829 TP: 0.7% -> 34.7%`
  - `H1437 TP: 24.5% -> 95.5%`
  - `H2009 TP: 64.5% -> 97.7%`
  - `HCC1395_DORADO TP: 46.9% -> 93.1%`

這代表舊版 TO LOH 統計不是「偏差一點」，而是**很多樣本本來就看不到大半 HP family reads**。

### 3.2 14 個 TP/FP 子集的欄位翻轉比例

以 7 個 TO run 的 TP/FP `significance_summary.csv` before/after 直接比對，得到：

| 欄位 | 最小翻轉比例 | 中位翻轉比例 | 最大翻轉比例 |
|---|---:|---:|---:|
| `VerificationClass` | 1.03% | 2.97% | 10.15% |
| `DominantLabel` | 9.17% | 27.85% | 47.66% |
| `LOH_Subtype` | 8.18% | 27.24% | 79.64% |
| `Potential_LOH` | 7.72% | 26.92% | 79.54% |
| `Quality_Tier` | 2.28% | 8.47% | 79.12% |
| `Quality_Score` | 8.35% | 27.53% | 80.93% |
| `HP_Ratio` | 39.26% | 75.20% | 95.40% |
| `HPMergedSig` | 8.12% | 19.03% | 41.65% |
| `HPFineSig` | 24.67% | 30.55% | 54.64% |
| `PairwiseMedianDist` | 0.00% | 0.00% | 0.00% |
| `CramersV` | 0.00% | 0.00% | 0.27% |
| `AlleleDelta` | 0.00% | 0.00% | 0.00% |

### 3.3 `hp_assign_rate` 與 agreement 重新量化

新版 `label_cluster_agreement.tsv` 已重建於  
`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_hp_fix_reanalysis/`

| 樣本 | 舊 `hp_assign_rate` 中位數 | 新中位數 | 舊 `<0.05` 比例 | 新 `<0.05` 比例 | 新 `agreement_positive` 比例 |
|---|---:|---:|---:|---:|---:|
| HCC1395 | 0.5634 | 0.9863 | 4.27% | 0.24% | 58.63% |
| HCC1395_DORADO | 0.5309 | 1.0000 | 7.09% | 0.05% | 53.33% |
| COLO829 | 0.3659 | 1.0000 | 5.61% | 0.49% | 28.68% |
| H1437 | 0.3506 | 1.0000 | 9.60% | 0.18% | 60.47% |
| H2009 | 0.4737 | 1.0000 | 6.86% | 0.09% | 49.40% |
| HCC1937 | 0.7468 | 1.0000 | 2.84% | 0.12% | 55.65% |
| HCC1954 | 0.7097 | 1.0000 | 2.21% | 0.26% | 67.51% |

結論很直接：**所有舊版 `hp_assign_rate` 與 `agreement_positive` 判讀都需要重做。**

---

## 4. 哪些舊結論仍成立，哪些必須改寫

### 4.1 仍可沿用的結論

以下內容原則上可以保留，只需在必要時補一句「已確認不受 HP integer tag fix 影響」：

1. TO caller benchmark 與 InterSubMod 最終 `F1 / precision / recall`
2. 以 caller AF、GQ、DP 為主的結論
3. 以 `PairwiseMedianDist`、`AlleleDelta` 為主的 feature-space 觀察
4. paired track 的 LOH round、方法學審查與 canonical paired 報告
5. `scripts/pipeline/steps/03_filter_analysis.py` 的 benchmark/filter 結論

### 4.2 必須降級或重寫的結論

以下 TO 結論不能直接沿用舊數字：

1. 任何把 `VerificationClass in Strong/Subclone` 當作 support 的論述
2. 任何使用 `agreement_positive` 的 rescue 效果表
3. 任何使用 `hp_assign_rate >= ...` 的 TO trust/QC 門檻
4. 任何使用 `Quality_Score >= ...` 的 TO rescue / annotation 統計
5. 任何使用 `Potential_LOH`、`LOH_Subtype`、`HP_Ratio`、`effective_hp` 的 LOH 富集分析
6. 任何把 TO `LOH-like` 解讀成 FP enrichment 的文字

### 4.3 修正後可成立的新口徑

1. **TO LOH 不是 FP filter，而是偏 TP enrichment 的結構訊號。**
2. **`hp_assign_rate` 在 TO 中更接近「phase completeness / mapping completeness」，不是可直接沿用舊值的真偽特徵。**
3. **`Quality_Score` 雖然不是全部翻盤，但因其計算已受 HP family 修正牽動，舊閾值表必須重算後才能再引用。**
4. **`PairwiseMedianDist` 與 `AlleleDelta` 是這次修正後最穩定、可跨版本比較的欄位。**
5. **HP0 的方向性結論仍成立，但應改引用修正後數字。**  
   根據 `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_corrections/corrected_numbers_summary.json`：  
   `HP0 越高 TP 比例越高趨勢仍成立（0.748→0.757→0.767）`。

---

## 5. 需要重新執行的項目

### 5.1 不必重跑的 upstream

以下流程不用重做：

- `step01_clairs_to`
- `step02_benchmark_clairs_to`
- `step03_longphase_to`
- `step04_benchmark_longphase_to`
- 各樣本 benchmark comparison 與 stage metrics 的 caller-side F1 摘要

### 5.2 必須重跑的 downstream

| 優先級 | 類別 | 項目 | 建議動作 | 狀態 |
|---|---|---|---|---|
| P0 | TO summary 後處理 | 7 個 TO run 的 `step05_intersubmod` 衍生驗證表 | 以修正版 `significance_summary.csv` 重建 `method_design_validation.tsv` 與 `label_cluster_agreement.tsv` | **已完成** |
| P0 | TO LOH 更正 | `build_to_loh_enrichment_post_hp_fix.py` 與相關 LOH correction workspace | 作為新的 TO LOH 正式口徑 | **已完成** |
| P0 | LOH cross-sample workspace | `build_loh_round1_cross_sample_audit.py` | 重新建立含新版 TO rows 的 round1 workspace | 待執行 |
| P0 | LOH round 2 | `build_loh_round2_support_hp0_analysis.py` | 重新輸出 TO support/HP0 分層 | 待執行 |
| P0 | LOH round 2 | `build_loh_round2_ps_export_and_to_block_audit.py` | 重新檢視 TO-only block audit 與 block-level統計 | 待執行 |
| P0 | LOH round 3 | `build_loh_round3_methyl_hp0_filter.py` | 重新檢查 Tier A LOH + HP0 + methyl 聯合規則 | 待執行 |
| P0 | Phase 1 dataset | `build_phase1_training_manifest.py` | 重建含 TO summary 欄位的 manifest | 待執行 |
| P0 | Phase 1 dataset | `export_phase1_read_training_table.py` | 重建含 `VerificationClass / agreement_type / Quality_Score / hp_assign_rate` 的讀層表 | 待執行 |
| P0 | TO provenance | `build_to_fp_provenance_analysis.py` | 重新輸出 `quality_score / verification_class / potential_loh / agreement_type` 對應的 TO provenance 分析 | 待執行 |
| P1 | TO AF gradient | `analyze_cross_sample_ism_af_gradient.py` | 若仍引用 `Quality_Score`、`VerificationClass`，需重跑；若只保留 `Pairwise/CramersV` 可局部修文 | 待判定 |
| P2 | AF only | `analyze_cross_sample_to_af_loh.py` | 腳本本身只看 caller AF，不受本次修正影響 | **不需重跑** |
| P2 | F1 summary | `collect_baseline_metrics.py` | 只讀 stage metrics/F1，與此次修正無關 | **不需重跑** |

### 5.3 已完成的修正輸出

1. 新版 `method_design_validation.tsv` / `label_cluster_agreement.tsv`  
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_hp_fix_reanalysis/{run_id}/`

2. TO LOH 修正後 workspace  
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_enrichment_post_hp_fix/`

3. TO LOH 校正文案與數字  
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_corrections/`

---

## 6. 需要修正或重寫的文件

### 6.1 必須重寫數字與結論的 TO 文件

- `docs/experiments/in_progress/2026/03/20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md`
- `docs/experiments/in_progress/2026/03/20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md`
- `docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md`
- `docs/experiments/in_progress/2026/03/20260308_TO_cluster_only_StrongWeak_diagnostics與scheme調整分析_01.md`
- `docs/experiments/in_progress/2026/03/20260309_5kHz_TO與DORADO_paired_甲基rescue特徵空間分析_01.md`
- `docs/experiments/in_progress/2026/03/20260311_phase2_finer_interval回接annotation_layer驗證_01.md`
- `docs/reports/validated/2026/03/20260322_TO_rescue_rules_autoresearch_report_01.md`
- `docs/reports/validated/2026/03/20260322_cross_sample_TO_ISM_gradient_analysis_01.md`
- `docs/reports/validated/2026/03/20260322_TO_FP_provenance_methylation_analysis_01.md`
- `docs/reports/validated/2026/03/20260322_TO_TP_FP_germline_LOH_characterization_01.md`
- `docs/reports/validated/2026/03/20260322_TO_FP來源分解與paired對照分析_01.md`
- `docs/reports/validated/2026/03/20260322_TO_FP來源分解摘要_01.md`
- `docs/reports/validated/2026/03/20260323_TO_residual_FP_deep_dive_01.md`
- `docs/reports/validated/2026/03/20260324_研究週報_20260318_20260324_方法學審查與TO收尾整合_01.md`

### 6.2 可保留但需補註修正的文件

- 任何只把 TO `PairwiseMedianDist / AlleleDelta / caller AF` 當作主要證據，且未直接引用 `VerificationClass / agreement_positive / hp_assign_rate / Quality_Score / LOH-like` 舊數字的文件
- paired-only 文件
- 修正後新寫的 `20260330_TO_LOH_enrichment_post_hp_fix_01.md`

---

## 7. 建議的修正順序

1. **先凍結新資料源**  
   以 7 個修正版 TO `significance_summary.csv` 與 `20260330_to_hp_fix_reanalysis/` 為唯一新基準。

2. **先補 TO LOH 主線**  
   重建 `Round1 -> Round2 -> Round3` 的 TO 相關工作區，確保 LOH 與 HP0 的正式口徑一致。

3. **再補 candidate-specific rescue 與 Phase 1 匯出**  
   因為這些表直接使用 `VerificationClass / agreement_type / Quality_Score / hp_assign_rate`。

4. **最後修正文稿與週報**  
   把所有舊版 TO 文字統一改成：「舊版 HP integer tag 解析錯誤，2026-03-30 後數字已更正」。

---

## 8. 最終判決

本次修正**沒有推翻 TO track 的全部研究**，但它**明確推翻了所有依賴 TO HP family / LOH / label agreement 舊數字的結論**。

最合理的邊界是：

- **保留**：caller benchmark、AF、`PairwiseMedianDist`、`AlleleDelta`、paired 主線
- **重算**：`VerificationClass / agreement / hp_assign / Quality_Score / LOH-like` 相關 TO 分析
- **改寫口徑**：TO LOH 是 TP enrichment，不是 FP filter

因此，這不是「全部重做」，而是一次**邊界非常清楚的 TO downstream correction**。

---

## 9. 主要依據

1. `src/core/ReadParser.cpp:121-141`
2. `docs/experiments/in_progress/2026/03/20260330_TO_HP_integer_tag_fix與全樣本重跑結果確認_01.md`
3. `docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md`
4. `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_hp_fix_reanalysis/`
5. `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_corrections/corrected_numbers_summary.json`
6. Knowledge Base:
   - `Knowledge/05_tools/longphase_to.md`
   - `Knowledge/06_workflows/phasing_workflow.md`
