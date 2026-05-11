---
title: Step 7 — HCC1395 Normal BAM pilot 完整推理敘述與視覺化（TO + paired 雙模式獨立包裝）
date: 2026-04-19
status: CONDITIONAL（單樣本；跨樣本驗證前不升級 CL-025a/b）
owner: InterSubMod Research
scope: Batch 2 R1 · G.1.2 雙模式獨立包裝 · Opus 4.7 plan F+G 決策執行
related:
  - step7_hcc1395_normal_to_pilot.md
  - step7_hcc1395_normal_paired_pilot.md
  - ../../../docs/reports/research_landscape/10_Research_Chain_Registry.md
  - ../../../docs/reports/research_landscape/06_結論穩定性審查.md
---

# Step 7 — HCC1395 Normal BAM pilot 完整推理敘述

**Date**: 2026-04-19  
**Status**: CONDITIONAL（HCC1395 single sample；殘差化後部分訊號保留，但跨樣本驗證前不升級 ⭐）

## 0. 動機與 Registry 位置

**研究鏈**：Opus 4.7 plan F.7 鎖定決策 Q6r → 「TO 與 paired 兩模式並行、獨立包裝」→ G.1.2 雙 finding 輸出規範。本文件為 R1 pilot 完整推理敘述與視覺化產出。

**目標**：Phase 2 A+D 程式碼（2026-04-13 完工，173/173 tests）在 HCC1395 上的 **生物學聲稱驗證** 第一步：

- **TO arm（主路）**：Normal BAM 整合後新增的特徵（`NormalBaseline_Coverage`, `SampleASM_Delta`, `HP_Residual_Delta`）是否為 F1-improvement 提供 filter-potential 訊號？
- **Paired arm（驗證路）**：Paired mode 飽和（CL-017）聲稱是否只針對 region-level？Per-CpG ASM 是否在 paired 下有 characterization-beyond-F1 訊號？

**Registry 連結**：CL-013（Phase 2 程式碼完工）/ CL-014（Phase B/C/D 單樣本）→ 本 pilot 新增 CL-025a（TO）/ CL-025b（paired）。

## 1. Pipeline 與子集設計

### 1.1 輸入資料流

```
HCC1395 master dataset (all_region_rows.tsv.gz, 748K rows)
   ↓  F pilot canonical filter: NG=4 + AF<0.4 + NR≥80 + NonLOH
   ├── TO:     860 filter SNV loci (TP=763, FP=97)
   └── paired: 983 filter SNV loci (TP=971, FP=12)
   ↓  bcftools merge TP+FP VCFs, view -R regions
   ├── TO subset VCF:     801 variants (after dedup)
   └── paired subset VCF: 983 variants
   ↓  inter_sub_mod Phase 2 pipeline (dual-BAM)
   ├── TO:     build/bin/inter_sub_mod  -t TO_tumor_tagged.bam  -n normal_subset.bam  --loh-bed tumor_phased_LOH.bed
   └── paired: build/bin/inter_sub_mod  -t paired_Tmode.bam     -n normal_subset.bam  (no --loh-bed; ClairS paired 內部 LOH)
   ↓  significance_summary.csv + Sample ASM / HP Residual / Per-CpG Fisher 輸出
   ↓  Δ_ASM analysis + residualization + bootstrap
```

### 1.2 Normal BAM region-subset

**關鍵決策（使用者授權）**：從 `/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam` (136GB) 僅提取 F pilot region 子集 → `output/hcc1395_normal_pilot/HCC1395_normal_subset.bam` (1.8GB, 76× 壓縮)。

**技術可行性 gate**：TO arm 採 **inference-time annotation** 策略（normal ASM 作為 per-region feature 加入，不改 TO variant calling 本身）— 維持 TO mode 純粹性（G.1.3 檢核）。

### 1.3 執行時間

| Arm | Regions | Wall-time | Processing avg |
|-----|---------|-----------|---------------|
| TO | 801 | 64.0s | 1222 ms/region |
| paired | 983 | 72.8s | 1156 ms/region |

完全落在 G 段「<2h initial」時間 budget 內；region-subset 策略證明有效。

## 2. 結果：全子集 per-feature AUC

![Figure 1 — per-feature AUC (bootstrap 500× 95% CI) for both arms](../figures/fig1_feature_auc_bar.png)

### 2.1 TO arm（F pilot subset n=743, TP=665, FP=78）

**Top 5 特徵**（按 bootstrap mean）：

| Feature | AUC | 95% CI |
|---------|-----|--------|
| `NormalBaseline_Coverage` | **0.686** | [0.609, 0.753] |
| `SampleASM_Delta` | **0.645** | [0.580, 0.704] |
| `Normal_HP_Delta` | 0.587 | [0.513, 0.652] |
| `HP_Residual_Delta` | 0.578 | [0.521, 0.639] |
| `HPFineF` | 0.578 | [0.508, 0.655] |

**觀察**：
- **兩個 Normal BAM 原生特徵（`NormalBaseline_Coverage`, `SampleASM_Delta`）佔據前二名**。
- 傳統 ISM 特徵在 F pilot pre-enriched subset 內飽和（`HPFineNGroups` AUC=0.500 = 完全退化），**符合 CL-008（Beyond-AUC 耗盡）預測**。
- CI 下界皆 >0.5 → 兩特徵訊號非 random。

### 2.2 Paired arm（F pilot subset n=942, TP=932, FP=10）

**Top 5 特徵**：

| Feature | AUC | 95% CI |
|---------|-----|--------|
| `Fisher_Frac_Sig` | **0.736** | [0.554, 0.890] |
| `Epipoly_Delta` | **0.699** | [0.518, 0.905] |
| `Tumor_HP_Delta` | 0.640 | [0.511, 0.783] |
| `Normal_HP_Delta` | 0.623 | [0.506, 0.819] |
| `HPFineF` | 0.617 | [0.505, 0.799] |

**觀察**：
- **`Fisher_Frac_Sig` AUC=0.736 首次突破 0.70** — 此前 pure methylation 特徵從未達此水準（CL-008 ≤0.58 上限）。
- **CI 極寬**（0.554-0.890）— 原因是 FP=10 小樣本（Wilson 分布尾端發散）。
- **Epipoly_Delta 0.699** — 第二獨立 per-CpG epiallele 指標亦 >0.65。

## 3. 反駁推論：是否為 confound artifact？

**核心質疑**：已知 `NormalBaseline_Coverage` 受 total coverage 影響；`Fisher_Frac_Sig` 受 `Fisher_N_Tested`（被測 CpG 數）影響。若純粹為 coverage drift artifact → 殘差化後 AUC 應大幅衰減至 <0.55。

### 3.1 殘差化設計

```
NormalBaseline_Coverage ~ NumReads + Coverage_Multiple (OLS)
SampleASM_Delta         ~ NumReads + Coverage_Multiple
Fisher_Frac_Sig         ~ Fisher_N_Tested + NumReads
Epipoly_Delta           ~ NumReads + Coverage_Multiple
```

### 3.2 殘差化結果

![Figure 2 — raw vs residualized AUC](../figures/fig2_residualized_auc_comparison.png)

| Feature | Mode | Raw AUC | Residualized AUC | Δ |
|---------|------|---------|------------------|---|
| `NormalBaseline_Coverage` | TO | 0.686 | **0.604** | −0.082 |
| `SampleASM_Delta` | TO | 0.645 | **0.610** | −0.035 |
| `Fisher_Frac_Sig` | TO | 0.545 | 0.532 | −0.013 |
| `Epipoly_Delta` | TO | 0.536 | 0.529 | −0.007 |
| `NormalBaseline_Coverage` | paired | 0.572 | 0.569 | −0.003 |
| `SampleASM_Delta` | paired | 0.591 | 0.591 | 0.000 |
| `Fisher_Frac_Sig` | paired | 0.736 | **0.698** | −0.038 |
| `Epipoly_Delta` | paired | 0.699 | **0.693** | −0.006 |

### 3.3 推理鏈

**TO arm**：
- `NormalBaseline_Coverage` 衰減 0.082 但殘差化後仍 **0.604 > 0.60 實用閾值** → **部分** coverage drift 驅動，但有獨立訊號。解讀：normal BAM 在 tumor region 的 coverage 差異，除了反映 sampling artifact，也反映 **normal baseline 內在異質性**（可能源自 germline het site 週邊的 ±Coverage 擺動）。
- `SampleASM_Delta` 衰減僅 0.035，殘差化後 **0.610** → 訊號對 coverage 幾乎獨立 → **Δ_ASM(tumor − normal) 提供 tumor-specific ASM 資訊，非 coverage artifact**。

**Paired arm**：
- `Fisher_Frac_Sig` 衰減 0.038，殘差化後 **0.698 > 0.65** → 對 `Fisher_N_Tested` confound **穩健**。
- 但 95% CI 下界 **0.534** → FP=10 太少，統計上**不能拒絕 AUC=0.5**。

### 3.4 TP vs FP 分布直觀比較

![Figure 3 — TO arm TP vs FP distribution (top-4 features)](../figures/fig3_tp_fp_distribution_TO.png)

![Figure 4 — paired arm TP vs FP distribution (top-4 features)](../figures/fig4_tp_fp_distribution_paired.png)

**觀察**：
- TO arm `NormalBaseline_Coverage` TP 與 FP 的中位數分離明顯（虛線可視化）→ 符合 AUC 0.604 殘差化結果。
- Paired arm `Fisher_Frac_Sig` TP 分布較集中於低值、FP 分布較偏向高值，但 FP 僅 10 點 → 直方圖不穩。

## 4. 整合結論（HCC1395 single sample）

### 4.1 結論升級 / 維持 / 降級

| CL id | 本次狀態變動 | 理由 |
|-------|------------|------|
| CL-008（Beyond-AUC 耗盡 ≤0.58） | **保留但加 caveat** | Phase 1A scope 內耗盡聲稱有效；Phase 2 新增 Per-CpG ASM Fisher + Normal BAM 特徵在 paired 下 AUC ≥0.60（FP 小樣本保留疑慮）；結論改為「Phase 1A 特徵空間耗盡；Phase 2 新特徵待跨樣本確認」 |
| CL-025a（TO arm Normal BAM POSITIVE） | **⭐3 CONDITIONAL**（未升級） | `SampleASM_Delta` 殘差化 AUC 0.610 穩定；但 single sample，FP=78 雖夠但僅 HCC1395 |
| CL-025b（paired Fisher CONDITIONAL） | **⭐3 → ⭐2 INSUFFICIENT**（降級） | CI 下界 0.534 跨入隨機；FP=10 小樣本不足支撐 |
| CL-013（Phase 2 code 完工） | 維持 ⭐4 | 程式碼 173/173 通過本次生物學聲稱部分支持 |
| CL-014（Phase B/C/D 單樣本驗證） | 維持 ⭐3 | 跨樣本 R12 仍未啟動 |

### 4.2 對 Opus 4.7 plan 隱含假設的回應

| # | 假設 | 本 pilot 驗證狀態 |
|---|------|-----------------|
| 1 | ClairS Paired 為 ground truth | 部分（SEQC2 truth set 有限；paired 子集 TP rate 0.988 顯示 pre-filtered truth） |
| 4 | Coverage_Multiple ≈ CN | ❌ FALSIFIED（R6 證實，本 pilot 殘差化進一步確認 coverage 需單獨處理） |
| 10 | Phase B/C/D 單樣本驗證足以支撐 Phase 2 生物學聲稱 | ❌ 本 pilot 顯示 HCC1395 有訊號但 CI 寬；**須跨樣本** |

### 4.3 新質疑

1. **Normal BAM region-subset 的 sampling bias**：只提取 F pilot region 的 normal reads → normal baseline 可能不含 background 對照；完整 genome-wide normal reference 是否會改變 AUC？
2. **Paired mode FP=10 是樣本 artifact 還是 pipeline 成熟度**：若 paired mode FP 如此稀少，filter 本身價值有限；characterization value 才是 paired arm 核心（CL-017 聲稱仍有效）。
3. **TO + paired 訊號獨立性未測**：TO `SampleASM_Delta` 與 paired `Fisher_Frac_Sig` 是否捕捉同一 biological signal？ensemble 是否有 additional gain？

## 5. 5 軸 Registry matrix 背景對照

![Figure 5 — 5-axis Registry matrix (mode × zone × AF band) TP rate](../figures/fig5_5axis_tp_rate_heatmap.png)

**觀察**：
- **Paired mode 幾乎全區 TP rate ≥0.98** → paired 本質飽和，filter 增益空間極小（CL-017 再次確認）。
- **TO mode Z5 NonLOH Half/Intermediate** TP rate 0.56-0.66 → 這是當前 TO FP 主要累積區（`data/10_Registry_5axis_matrix.tsv` 顯示 Z5 NonLOH Half CN1 n=85,725）。
- **Z3 LOH Extreme TO** TP rate 0.62 → CL-024 Z3 amplicon blacklist CONDITIONAL 結論對應區（HCC1954 特例）。

**本 pilot 訊號所在區**：F pilot subset 主要落在 **Z1/Z2 NonLOH × AF 0.05-0.4 × CN2** 帶（validated_main_signal + pilot_subclone）；Normal BAM 特徵的效應主要在這些 cell 內累積。

## 6. 下一步行動優先序

### P0（不需決策，可立即啟動）

- **無**（R1 已完成）

### P1（需決策）

- **R12 跨樣本擴展**：把 Normal BAM pilot 擴展至其他 4 in-scope 樣本（HCC1937, HCC1954, H2009, HCC1395_DORADO；COLO829 R10 排除）
  - **阻塞**：需使用者授權逐一複製 Normal BAM（每個 136-296GB；region-subset 策略降至 1-3GB）
  - **成功標準**：TO `SampleASM_Delta` 至少 3/5 樣本保持 AUC ≥0.60（residualized）
  - **失敗後果**：CL-025a 降級為 HCC1395-specific；CL-008 保留完整聲稱

- **R8 Per-CpG ASM residualize 完整驗證**（批次 3 候選）：`Fisher_Frac_Sig` 在 paired arm 是否在更大 FP 數（其他樣本 FP 較多）下維持 AUC ≥0.65
  - 依賴 R12 執行完成

### P2（非阻塞）

- TO + paired cross-mode concordance 分析（本 pilot 已有數據，只需額外腳本）
- `NormalBaseline_Coverage` 與 `SampleASM_Delta` information gain 重疊分析（Spearman 相關 + 條件 AUC）

## 7. Provenance（完整路徑，可重現）

### 7.1 Scripts（本次產生）

| 路徑 | 功能 |
|-----|------|
| `research/F_hpfinengroups_deepening/scripts/r1_build_subset_vcfs.py` | F pilot filter + mode-aware NonLOH → filter SNV TSV |
| `research/F_hpfinengroups_deepening/scripts/r1_build_subset_vcf_bcftools.sh` | bcftools concat TP+FP + sort + view -R regions |
| `research/F_hpfinengroups_deepening/scripts/r1_run_phase2_both_arms.sh` | 並行啟動 TO + paired Phase 2 pipeline |
| `research/F_hpfinengroups_deepening/scripts/r1_delta_asm_analysis.py` | 合併 truth + 15 特徵 AUC |
| `research/F_hpfinengroups_deepening/scripts/r1_residualize_visualize.py` | OLS residualize + bootstrap 500× + 5 圖表 |
| `research/F_hpfinengroups_deepening/scripts/registry_5axis_matrix.py` | 5 軸 450-cell 組合矩陣（Fig 5 來源） |

### 7.2 Inputs

| 路徑 | 說明 |
|-----|------|
| Master dataset | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz` |
| TO tumor BAM | `/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam` (278GB) |
| Paired tumor BAM | `/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` |
| Normal BAM source | `/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam` (136GB) |
| TO LOH BED | `/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed` |
| Reference | `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta` |
| TO TP VCF | `/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step04_benchmark_longphase_to/filtered_snv_tp.vcf.gz` |
| TO FP VCF | 同上 `_fp` |
| Paired TP VCF | `/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395/20260211/longphase_s/filtered_snv_tp.vcf.gz` |
| Paired FP VCF | 同上 `_fp` |

### 7.3 Outputs

| 路徑 | 說明 |
|-----|------|
| `output/hcc1395_normal_pilot/HCC1395_normal_subset.bam` | 1.8GB region-subset normal BAM |
| `output/hcc1395_normal_pilot/HCC1395_TO_filter_subset.vcf.gz` | 801 variants TO subset |
| `output/hcc1395_normal_pilot/HCC1395_paired_filter_subset.vcf.gz` | 983 variants paired subset |
| `output/hcc1395_normal_pilot/TO/significance_summary.csv` | TO arm 801 rows × 103 columns |
| `output/hcc1395_normal_pilot/paired/significance_summary.csv` | paired arm 983 rows × 103 columns |
| `output/hcc1395_normal_pilot/logs/{to,paired}_arm.log` | inter_sub_mod 執行 log |

### 7.4 Figures + Data（本文件引用）

| 路徑 | 說明 |
|-----|------|
| `research/F_hpfinengroups_deepening/figures/fig1_feature_auc_bar.png` | 12 特徵 × 2 arm AUC + bootstrap 95% CI |
| `research/F_hpfinengroups_deepening/figures/fig1_feature_auc_bar.tsv` | Fig 1 數據 |
| `research/F_hpfinengroups_deepening/figures/fig2_residualized_auc_comparison.png` | raw vs residualized AUC 對比 |
| `research/F_hpfinengroups_deepening/figures/fig2_residualized_auc_comparison.tsv` | Fig 2 數據 |
| `research/F_hpfinengroups_deepening/figures/fig3_tp_fp_distribution_TO.png` | TO top-4 特徵 TP vs FP histogram |
| `research/F_hpfinengroups_deepening/figures/fig4_tp_fp_distribution_paired.png` | paired top-4 特徵 TP vs FP histogram |
| `research/F_hpfinengroups_deepening/figures/fig5_5axis_tp_rate_heatmap.png` | 5 軸 matrix mode × zone × AF band TP rate |
| `research/F_hpfinengroups_deepening/figures/summary.txt` | 所有 AUC 純文字彙總 |
| `docs/reports/research_landscape/data/10_Registry_5axis_matrix.tsv` | 240 cell 完整矩陣 |

## 8. 需使用者決策事項

**1. R12 跨樣本擴展是否啟動？**
   - 選項 A：漸進複製（先 HCC1954 Normal BAM 296GB → region-subset 2-3GB）；與 CL-016/018 HCC1954 CNV architecture 特例對照
   - 選項 B：一次複製 4 in-scope 樣本 Normal BAM region-subsets（總 ~10GB）
   - 選項 C：暫緩 R12，先做 P2 cross-mode concordance（TO vs paired 訊號獨立性）

**2. CL-008 聲稱是否應主動降級？**
   - 選項 A：保留完整聲稱，加 caveat「Phase 1A scope 內」
   - 選項 B：降級為 ⭐3 並標「Phase 2 新特徵待驗證」
   - 選項 C：等 R12 跨樣本後再決

**3. Paired mode CONDITIONAL 結論（CL-025b）處理**：
   - 本 pilot FP=10 統計不足；R12 其他樣本 paired FP 可能仍小（paired 飽和）
   - 選項 A：放棄 paired arm F1-filter 方向，轉為 characterization-only
   - 選項 B：R12 執行後若 2/4 樣本 CI 下界 >0.60 則保留
   - 選項 C：合併 4 樣本 FP 做 meta-analysis（避免個別樣本 FP 不足）

---

## 附錄 — 與 Registry 的連結

本 finding 更新 `docs/reports/research_landscape/10_Research_Chain_Registry.md` 以下 entries：
- CL-025a 狀態 active → 更新 `last_reviewed`，rating 維持 ⭐3
- CL-025b 狀態 active → rating 降級 ⭐3 → ⭐2 INSUFFICIENT，加 caveat「FP=10 CI 下界跨入隨機」
- 新增 R17（HCC1954 Normal BAM pilot 候選，待使用者 Q1 答覆啟動）

同步更新 `docs/reports/research_landscape/06_結論穩定性審查.md`（待 batch commit）。
