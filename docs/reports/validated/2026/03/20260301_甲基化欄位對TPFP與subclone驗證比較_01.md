<!--
建立時間: 2026-03-01 16:40
狀態: validated
資料來源:
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/advanced_analysis_20260119/advanced_samples.csv
  - data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz
  - data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz
  - /big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/tmp_meth_annot_test/subclone_validation_comparison_20260301.tsv
驗證命令:
  - python3 scripts/analysis/compare_subclone_validation.py
  - bcftools view -H -i 'INFO/Verdict_SubclonalSomatic=1' /big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz | wc -l
關聯文件:
  - scripts/analysis/compare_subclone_validation.py
  - Knowledge/03_file_formats/vcf_clairs_to.md
  - Knowledge/04_databases/seqc2_truth_set.md
-->

# 甲基化欄位對 TP/FP 與 subclone 驗證比較

## 1. 目標

回答問題：
1. TP/FP（有 truth set 時）是否能透過甲基化欄位得到可量化提升？
2. 特殊現象富集（LOH-like、高聚集區、低 AF 子群）是否有穩定結論？
3. 是否更能驗證 subclone？

## 2. 依據與前提

1. 根據 `Knowledge/04_databases/seqc2_truth_set.md`，HCC1395 的 truth set 為 SEQC2。
2. `scripts/pipeline/steps/01_longphase_s.sh` 顯示 TP/FP 由 `bcftools isec` 與 truth VCF 交集切分（0002=TP、0000=FP）。
3. `scripts/pipeline/config.sh` 中 HCC1395 truth total 固定為 39,447（SEQC2 HC sSNV）。
4. 根據 `Knowledge/03_file_formats/vcf_clairs_to.md` 與 VCF header，`Verdict_SubclonalSomatic` 為 `INFO Flag (Number=0, Type=Flag)`，不可直接承載連續數值。

## 3. 資料與方法

1. 主比較集：`advanced_samples.csv`（448 筆；TP 288、FP 160）。
2. 以 `filtered_snv_tp.vcf.gz / filtered_snv_fp.vcf.gz` 回填每筆位點的 `AF` 與 `H` flag。
3. 以各自 TP 或 FP 全集合，計算每筆位點 ±25kb 鄰域變異密度（`ClusterCount50kb`）。
4. 定義特徵：
   - `LOH_like`: `Potential_LOH` 或 `H=1`
   - `HighCluster`: `ClusterCount50kb > 20`
   - `LowAF`: `AF < 0.3`
   - `MethylSupport_SubcloneLike`: `VerificationClass in {Subclone, Strong}` 或 `HPMergedDelta <= 0.10`
5. 統計：Fisher exact test（OR + p-value）與規則式 F1。

## 4. 結果

### 4.1 TP/FP 基線（SEQC2 truth-label）

1. TP variants：30,490
2. FP variants：4,842
3. Truth total：39,447
4. 基線（TP/FP）：
   - Precision = 0.8630
   - Recall = 0.7729
   - F1 = 0.8155

### 4.2 特殊現象富集（448 筆比較集）

1. `LOH_like` 與 TP 顯著正相關  
   OR = 7.3297，p = 3.04e-16  
   TP rate(cond)=0.7415，TP rate(non-cond)=0.2813
2. `HighCluster` 與 TP 顯著負相關（偏向 FP）  
   OR = 0.0185，p = 5.75e-26  
   TP rate(cond)=0.0492，TP rate(non-cond)=0.7364
3. `LowAF` 與 TP 顯著負相關（低 AF 子群整體偏 FP）  
   OR = 0.0458，p = 4.44e-40  
   TP rate(cond)=0.3411，TP rate(non-cond)=0.9188

### 4.3 低 AF 子群下，甲基化/次克隆訊號是否有幫助

1. 在 LowAF 子群（n=214）：
   - 基線 TP rate = 0.3411
2. 使用 `MethylSupport_SubcloneLike`：
   - 支持組 TP rate = 0.3512
   - 非支持組 TP rate = 0.1111
   - OR = 4.3308，p = 0.1709（方向正向，但未達顯著）
3. 使用更嚴格 `VerificationClass=Subclone`：
   - 支持組 TP rate = 0.7778
   - 非支持組 TP rate = 0.3010
   - OR = 8.1271，p = 9.88e-05（顯著）

### 4.4 規則式比較（448 筆）

1. `AF>=0.3`：F1 = 0.8238
2. `AF>=0.3 OR VerificationClass=Subclone`：F1 = 0.8481（提升 +0.0244）
3. `SignificantOnly`：F1 = 0.1910（明顯不足）

## 5. Verdict_SubclonalSomatic 欄位可行性

1. `Verdict_SubclonalSomatic` 是 Flag，不可直接放甲基化數值。
2. `/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz` 中 `Verdict_SubclonalSomatic=1` 計數為 0。
3. 實作上可行方案：保留 `Verdict_*`，另外新增 `INFO` 數值欄位（例如 `ISM_CV`, `ISM_HPMD`, `ISM_VC`, `ISM_SF`）做交叉驗證。

## 6. 結論（對問題的直接回答）

1. 可以得到特殊結論：  
   `LOH_like` 富集 TP、`HighCluster` 與 `LowAF` 富集 FP，三者皆有強統計顯著。
2. 可以提升 subclone 驗證：  
   在低 AF 子群，`VerificationClass=Subclone` 可顯著提高 TP 富集，且規則 `AF>=0.3 OR Subclone` 的 F1 高於單純 AF 門檻。
3. 若問「是否可直接把甲基化值加到 Verdict_SubclonalSomatic」：  
   不可（欄位型別不符），應改為新增平行 ISM 欄位。

## 7. 限制與再實驗建議

1. 本次核心比較集為 `advanced_samples.csv`（分層抽樣，每類 32 筆），不等於全量分佈。
2. 建議再做一次全量實驗：
   - 對全 TP/FP 位點批次補齊 ISM 欄位。
   - 在全量資料重算 OR/F1 與門檻敏感度。
3. 建議新增校準：
   - 以 `VerificationClass=Subclone` + `HPMergedDelta` + `ClusterCount50kb` 建立 logistic / tree 模型，與現行規則比較。

## 8. 產物路徑

1. `scripts/analysis/compare_subclone_validation.py`
2. `/big8_disk/liaoyoyo2001/InterSubMod_runs/output/tmp_meth_annot_test/subclone_validation_comparison_20260301.tsv`
3. `/big8_disk/liaoyoyo2001/InterSubMod_runs/output/tmp_meth_annot_test/subclone_validation_comparison_20260301_rules.tsv`
4. `/big8_disk/liaoyoyo2001/InterSubMod_runs/output/tmp_meth_annot_test/subclone_validation_comparison_20260301_samples.tsv`
