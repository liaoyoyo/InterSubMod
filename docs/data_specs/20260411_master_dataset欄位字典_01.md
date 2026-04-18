<!--
建立時間: 2026-04-11 21:30
目標: Master Dataset (all_region_rows.tsv.gz) 全部 116 欄完整欄位字典
處理範圍: 59 欄 C++ 原始欄位 + 57 欄 Python 聚合欄位，含計算邏輯、值域、模式差異
關聯檔案:
  - scripts/analysis/build_loh_round1_cross_sample_audit.py (load_dataset_rows, line 586)
  - scripts/analysis/observation_common.py (MASTER_DATASET_PATH, line 41)
  - docs/data_specs/20260411_significance_summary欄位字典_01.md
-->

# Master Dataset 欄位字典 (all_region_rows.tsv.gz)

> **版本**: v1.0 | **日期**: 2026-04-11 | **欄位數**: 116 (59 原始 + 57 聚合)

**位置**: `big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`
**定義**: `scripts/analysis/observation_common.py:41-44` (`MASTER_DATASET_PATH`)
**規格**: 748,391 rows × 116 columns, 84.9 MB (gzipped TSV), 建立日期 2026-03-30
**建構腳本**: `scripts/analysis/build_loh_round1_cross_sample_audit.py` (`load_dataset_rows()`, line 586)

---

## 欄位總覽

### Part A: C++ 原始欄位 (#1-59)

來自 `significance_summary.csv`，由 C++ 核心處理直接輸出。完整說明見 `20260411_significance_summary欄位字典_01.md`。

| 群組 | 欄位 | 說明 |
|------|------|------|
| 區域識別 | RegionID, Chr, Pos, Ref, Alt | SNV 座標 |
| 基本統計 | NumReads, NumCpGs | 深度與 CpG 數 |
| 全域檢驗 | GlobalP, CramersV, HeuristicScore, PassedGating | Fisher + V + 門檻 |
| 距離度量 | PairwiseMeanDist, PairwiseMedianDist | Read 間距離 |
| 分群 PERMANOVA | ClusterPermanovaF/P/Valid, ClusterDispersionP/Warn | 無監督分群檢驗 |
| Stage 1: HP Merged | HPMergedDelta/P/Sig, HP1FamilyN, HP2FamilyN | Haplotype 合併 |
| Stage 2: HP Fine | HPFineF/P/Sig, HPFineNGroups | 精細 4 群 |
| Stage 3: Allele | AlleleDelta/P/Sig | ALT vs REF |
| Label PERMANOVA | LabelHP/Allele PermanovaF/P/Valid, DispersionP/Warn | 標籤驗證 (10 欄) |
| Stage 4: Unassigned | UnassignedAffinity/P/Dir, NHP3, NHP0 | 未分配 read |
| 判定 | DominantLabel, Stability, VerificationClass | 分類結果 |
| 品質 | HP_Ratio, Potential_LOH, Coverage_Multiple/Category, LOH_Subtype, Quality_Score/Tier | 品質指標 |
| 最終判定 | LocalBestCluster, LocalBestP, Significant, SuggestFilter | 顯著性 |

---

### Part B: Python 聚合欄位 (#60-116)

由 `build_loh_round1_cross_sample_audit.py:load_dataset_rows()` 在聚合階段加入。

---

## B1. 身份與標籤 (#60-68)

| # | 欄位名 | 型態 | 值域 | 說明 | 計算來源 |
|---|--------|------|------|------|---------|
| 60 | **truth_label** | string | TP / FP | Truth set 比對結果。TP = 在 truth VCF + high-confidence BED 內 | `bcftools isec` 比對結果 |
| 61 | **variant_key** | string | `chr:pos:ref:alt` | 唯一變異識別碼，格式 `chr1:877772:G:C`。用於跨表 merge | `Chr:Pos:Ref:Alt` 拼接 |
| 62 | **sample** | string | HCC1395, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954 | 樣本名稱 | `DatasetConfig.sample` |
| 63 | **sample_label** | string | 如 `HCC1395 (5kHz)`, `HCC1395 (Dorado)` | 含平台標籤的樣本顯示名 | `DatasetConfig.sample_label` |
| 64 | **platform** | string | ONT_5kHz, ONT_Dorado | 測序平台/basecaller 版本 | `DatasetConfig.platform` |
| 65 | **mode** | string | paired / to | 分析模式。paired=有配對正常樣本，to=tumor-only | `DatasetConfig.mode` |
| 66 | **mode_label** | string | Paired, TO | 顯示用模式名稱 | `DatasetConfig.mode_label` |
| 67 | **dataset_id** | string | 如 `HCC1395_paired_full` | 唯一 dataset 識別碼 | `DatasetConfig.dataset_id` |
| 68 | **dataset_label** | string | 如 `HCC1395 Paired Full` | 顯示用 dataset 名稱 | `DatasetConfig.dataset_label` |

---

## B2. 來源追蹤 (#69-77)

| # | 欄位名 | 型態 | 值域 | 說明 |
|---|--------|------|------|------|
| 69 | **source_kind** | string | canonical / archive | 資料來源類型：canonical=正式 run，archive=歷史歸檔 |
| 70 | **source_priority** | int | 1-3 | 資料優先級。1=最高（canonical full），3=最低（archive） |
| 71 | **source_base_dir** | string | 絕對路徑 | 該 run 的根目錄路徑 |
| 72 | **source_summary_file** | string | 絕對路徑 | significance_summary.csv 的完整路徑 |
| 73 | **source_region_root** | string | 絕對路徑 | Per-region 輸出的根目錄（TP 或 FP 子目錄） |
| 74 | **source_vcf_file** | string | 絕對路徑 | 輸入 VCF 檔案路徑（TP VCF 或 FP VCF） |
| 75 | **source_tagged_bam** | string | 絕對路徑 | Haplotagged tumor BAM 路徑 |
| 76 | **source_phased_vcf** | string | 絕對路徑或空 | Phased VCF 路徑。Paired 模式為空 |
| 77 | **source_loh_bed** | string | 絕對路徑或空 | LOH BED 檔案路徑（如有） |

---

## B3. Context (#78)

| # | 欄位名 | 型態 | 值域 | 說明 |
|---|--------|------|------|------|
| 78 | **context_truth_total** | string/int | 數字或空 | Truth set 中的總變異數（從 context JSON 讀取） |

---

## B4. 衍生 HP 指標 (#79-86)

由原始 HP 欄位重新計算，解決 C++ 版 Laplace smoothing 與 Python 版直接比例的差異。

| # | 欄位名 | 型態 | 值域 | 說明 | 計算公式 |
|---|--------|------|------|------|---------|
| 79 | **effective_hp_reads** | int | ≥0 | HP1FamilyN + HP2FamilyN。不含 HP0/HP3 | `HP1FamilyN + HP2FamilyN` |
| 80 | **hp_ratio_core** | float | [0, 1] 或 NaN | HP1 比例（不含 Laplace smoothing）。NaN 當 effective_hp_reads=0 | `HP1FamilyN / effective_hp_reads` |
| 81 | **hp_assign_rate** | float | [0, 1] 或 NaN | HP 分配率。低值表示大量 unphased reads | `effective_hp_reads / NumReads` |
| 82 | **hp0_ratio** | float | [0, 1] 或 NaN | Unphased reads 比例 | `NHP0 / NumReads` |
| 83 | **hp3_ratio** | float | [0, 1] 或 NaN | Multi-mapped reads 比例 | `NHP3 / NumReads` |
| 84 | **tool_potential_loh** | bool | True/False | C++ Potential_LOH 欄位轉為 Python bool | `Potential_LOH` 轉型 |
| 85 | **core_loh_like** | bool | True/False | 以 hp_ratio_core 重新判定的 LOH-like。**與 C++ 版差異**: 無 Laplace smoothing，且需 effective_hp_reads > 0 | `(hp_ratio_core < 0.1) \| (hp_ratio_core > 0.9)` |
| 86 | **tool_core_loh_match** | bool | True/False | C++ LOH 判定與 Python 重算是否一致。不一致的原因通常是 Laplace smoothing 差異 | `tool_potential_loh == core_loh_like` |

**hp_ratio_core vs HP_Ratio 差異說明**:
- `HP_Ratio` (C++ #49): 使用 Laplace smoothing `(HP1+0.001)/(HP1+HP2+0.002)`
- `hp_ratio_core` (Python #80): 直接比例 `HP1/(HP1+HP2)`，effective_hp_reads=0 時為 NaN
- 兩者在 read 數少時差異較大；研究中通常使用 `hp_ratio_core`

---

## B5. 分箱欄位 (#87-88)

| # | 欄位名 | 型態 | 值域 | 說明 |
|---|--------|------|------|------|
| 87 | **hp_ratio_bin** | string | `[0.0-0.1)`, `[0.1-0.2)`, ..., `[0.9-1.0]`, `too_few_hp` | hp_ratio_core 十分位分箱。effective_hp_reads < 10 時為 `too_few_hp` |
| 88 | **verification_class** | string | Strong / Subclone / Weak / Noise / Unknown | VerificationClass 的小寫複製。**注意**: 與 #48 大寫版相同內容，小寫欄名方便 Python 存取 |

---

## B6. 小寫複製欄位 (#89-94)

原始 C++ 欄位的小寫別名，便於 Python 分析腳本一致使用 snake_case。

| # | 欄位名 | 型態 | 對應原始欄位 | 說明 |
|---|--------|------|-------------|------|
| 89 | **loh_subtype** | string | LOH_Subtype (#53) | 小寫複製 |
| 90 | **allele_delta** | float | AlleleDelta (#28) | 小寫複製 |
| 91 | **pairwise_median_dist** | float | PairwiseMedianDist (#13) | 小寫複製 |
| 92 | **quality_score** | float | Quality_Score (#54) | 小寫複製 |
| 93 | **cramers_v** | float | CramersV (#9) | 小寫複製 |
| 94 | **num_reads** | int | NumReads (#6) | 小寫複製 |
| 95 | **num_cpgs** | int | NumCpGs (#7) | 小寫複製 |

---

## B7. Caller 欄位 (#96-105)

從原始 somatic VCF（ClairS / ClairS-TO）提取的 variant caller 特徵。由 `load_split_vcf_features()` 解析 VCF INFO/FORMAT 欄位。

| # | 欄位名 | 型態 | 值域 | 說明 |
|---|--------|------|------|------|
| 96 | **caller_filter** | string | PASS / LowQual / ... | VCF FILTER 欄位值 |
| 97 | **caller_af** | float | [0, 1] | Tumor allele frequency (VAF)。**研究發現**: TO 模式中 AUC=0.654，超越所有 ISM 特徵 |
| 98 | **caller_dp** | int | ≥0 | Tumor total depth (DP) |
| 99 | **caller_gq** | float | ≥0 | Genotype quality (GQ) |
| 100 | **caller_ad_ref** | int | ≥0 | Tumor REF allele depth (AD[0]) |
| 101 | **caller_ad_alt** | int | ≥0 | Tumor ALT allele depth (AD[1]) |
| 102 | **caller_naf** | float | [0, 1] 或 NaN | Normal allele frequency。TO 模式為 NaN（無 normal sample） |
| 103 | **caller_ndp** | int 或 NaN | ≥0 或 NaN | Normal total depth。TO 模式為 NaN |
| 104 | **caller_nad_ref** | int 或 NaN | ≥0 或 NaN | Normal REF depth。TO 模式為 NaN |
| 105 | **caller_nad_alt** | int 或 NaN | ≥0 或 NaN | Normal ALT depth。TO 模式為 NaN |

---

## B8. 分箱欄位 (#106-107)

| # | 欄位名 | 型態 | 值域 | 說明 |
|---|--------|------|------|------|
| 106 | **af_bin** | string | `[0.0-0.1)`, `[0.1-0.2)`, ..., `[0.9-1.0]`, `NA` | caller_af 十分位分箱 |
| 107 | **effective_hp_bin** | string | `[0-10)`, `[10-30)`, `[30-50)`, `[50-80)`, `[80+)` | effective_hp_reads 分箱 |

---

## B9. TO Phasing 欄位 (#108-116)

Tumor-Only 模式的 phase block 資訊。從 phased VCF 提取。Paired 模式這些欄位皆為空字串。

| # | 欄位名 | 型態 | 值域 | 說明 |
|---|--------|------|------|------|
| 108 | **to_phase_filter** | string | PASS / ... 或空 | TO phased VCF 的 FILTER |
| 109 | **to_phase_ps** | string | phase set ID 或空 | Phase Set ID (PS tag)。同一 PS 內的變異被認為在同一 phase block |
| 110 | **to_phase_gt** | string | 0\|1 / 1\|0 / ... 或空 | Phased genotype (GT) |
| 111 | **to_phase_gt2** | string | 備用 GT 或空 | 第二組 phased genotype（multi-sample 場景） |
| 112 | **to_phase_ps2** | string | 備用 PS 或空 | 第二組 PS |
| 113 | **to_phase_gt3** | string | 備用 GT 或空 | 第三組 GT |
| 114 | **to_phase_ps3** | string | 備用 PS 或空 | 第三組 PS |
| 115 | **to_loh_bed_hit** | bool | True/False | 該位點是否落在 LOH BED 區域內。Paired 模式固定為 False | `BedIndex.contains(chr, pos)` |
| 116 | **phase_block_status** | string | to_ps_available / to_ps_missing / paired_vcf_ps_missing | Phase block 狀態。TO 有 PS=available，TO 無 PS=missing，Paired=paired_vcf_ps_missing |

---

## 模式差異速查

| 欄位群組 | Paired 模式 | TO 模式 |
|---------|------------|---------|
| Normal caller 欄位 (#102-105) | 有值 | NaN |
| TO phasing 欄位 (#108-114) | 空字串 | 有值（如 phased VCF 存在） |
| to_loh_bed_hit (#115) | False | True/False |
| phase_block_status (#116) | `paired_vcf_ps_missing` | `to_ps_available` 或 `to_ps_missing` |
| HP_Ratio 可靠度 | 高（germline phasing） | 低（self-phasing artifact，94.6% HP1 bias） |
| Quality_Score 有效性 | 有效 | 無效 (AUC=0.497) |

---

## 大寫/小寫重複欄位對照

以下欄位存在大寫（C++ 原始）和小寫（Python 複製）兩個版本。**內容相同**，使用小寫版即可。

| C++ 原始 (PascalCase) | Python 複製 (snake_case) |
|-----------------------|-------------------------|
| VerificationClass (#48) | verification_class (#88) |
| LOH_Subtype (#53) | loh_subtype (#89) |
| AlleleDelta (#28) | allele_delta (#90) |
| PairwiseMedianDist (#13) | pairwise_median_dist (#91) |
| Quality_Score (#54) | quality_score (#92) |
| CramersV (#9) | cramers_v (#93) |
| NumReads (#6) | num_reads (#94) |
| NumCpGs (#7) | num_cpgs (#95) |
| Potential_LOH (#50) | tool_potential_loh (#84) — **注意**: bool 轉型後可能有差異 |
