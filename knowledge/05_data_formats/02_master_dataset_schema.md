---
id: ism-kb-05-data-formats-master-dataset-schema
name: "master_dataset 116 欄 Schema"
description: "Python 聚合 7 樣本 × 2 模式 × 748K regions 產生的 master dataset (all_region_rows.tsv.gz) 欄位總覽：59 原始 + 57 Python 衍生。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "schema against docs/data_specs/20260411_master_dataset欄位字典_01.md"
related_ids:
  - ism-kb-05-data-formats-index
  - ism-kb-05-data-formats-significance-summary-schema
  - ism-kb-07-derived-features-coverage-multiple
  - ism-kb-06-workflows-analysis-scripts-index
  - ism-kb-05-data-formats-merged-dataset-pitfalls
tags: [data-format, master-dataset, tsv, schema, aggregated]
canonical_paths: [05_data_formats/02_master_dataset_schema.md]
alias_paths: []
---

# master_dataset 116 欄 Schema

- 一句結論：Python 聚合產出，含 59 C++ 原始欄 + 57 Python 衍生欄；共 748,391 rows × 7 samples × 2 modes
- 適用對象：跨樣本分析、衍生特徵研究
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz
  ```

---

## 生成腳本

- **腳本**：`scripts/analysis/build_loh_round1_cross_sample_audit.py`
- **函式**：`load_dataset_rows()`, `aggregate_master_dataset()`
- **輸入**：多個 `significance_summary.csv` + VCF metadata + LongPhase outputs
- **輸出**：`output/synthesis/observation_workspaces/.../all_region_rows.tsv.gz`

---

## 欄位分兩部分

### Part A: C++ 原始欄（#1-59）
對應 `significance_summary.csv` 59 欄，詳見 [01_significance_summary_schema.md](01_significance_summary_schema.md)。

### Part B: Python 衍生欄（#60-116）57 欄

#### 身份與標籤（9 欄）
| 欄位 | 內容 |
|------|------|
| `truth_label` | TP / FP |
| `variant_key` | `{chr}_{pos}_{ref}_{alt}` |
| `sample` | HCC1395, COLO829, H1437, H2009, HCC1937, HCC1954, HCC1395_DORADO |
| `sample_label` | 人類可讀標籤 |
| `platform` | ONT / ONT_R10 |
| `mode` | paired_full, paired_pileup, TO |
| `mode_label` | 人類可讀 |
| `dataset_id` | `{sample}_{mode}` |
| `dataset_label` | — |

#### 來源追蹤（9 欄）
| 欄位 | 內容 |
|------|------|
| `source_kind` | canonical / archive / legacy |
| `source_priority` | 1-5（canonical = 1） |
| `source_base_dir` | 路徑根 |
| `source_summary_file` | significance_summary.csv 路徑 |
| `source_region_root` | region 子目錄根 |
| `source_vcf_file` | VCF 來源 |
| `source_tagged_bam` | haplotagged BAM |
| `source_phased_vcf` | phased VCF |
| `source_loh_bed` | LOH BED |

#### Context（1 欄）
| 欄位 | 內容 |
|------|------|
| `context_truth_total` | 該樣本 truth set 總變異數 |

#### 衍生 HP 指標（8 欄）
| 欄位 | 定義 |
|------|------|
| `effective_hp_reads` | 實際有 HP tag 的 reads 數 |
| `hp_ratio_core` | HP1-family / effective_hp_reads（無 smoothing） |
| `hp_assign_rate` | effective_hp_reads / NumReads |
| `hp_skew` | `|HP1 - HP2| / (HP1 + HP2)` |
| `hp1_fraction`, `hp2_fraction` | 各半型比例 |
| `hp_somatic_fraction` | HP:i:11/21/33 比例 |
| `hp_unassigned_fraction` | HP:i:0/3 比例 |

#### 其他衍生（30+ 欄）
- LOH 相關：`loh_zone`, `cnLOH_proxy`, `LOH_confidence_class`, ...
- 覆蓋度代理：`coverage_quartile`, `log10_numreads`, ...
- 交叉樣本一致性：`per_sample_fisher_p`, `cross_sample_concordance`, ...
- Fisher Frac Sig 相關：`fisher_frac_sig`, `Fisher_Frac_CI_low/high`, ...
- Zone-Aware：`zone_label`, `zone_confidence`, ...

---

## 載入範例

```python
import pandas as pd

df = pd.read_csv(
    'output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz',
    sep='\t',
    low_memory=False
)
print(df.shape)        # (748391, 116)
print(df.columns[:60]) # Part A (59 欄)
print(df.columns[60:]) # Part B (57 欄)

# 典型篩選：HCC1395 paired_full TP
subset = df.query("sample == 'HCC1395' and mode == 'paired_full' and truth_label == 'TP'")
```

---

## 典型分析

- **跨樣本一致性**：`groupby(['sample', 'mode']).agg(...)` 計算 cross-sample Fisher
- **TP/FP 特徵比較**：`groupby('truth_label').describe()`
- **Zone-Aware 分析**：按 `zone_label` 分層
- **HPFineNGroups subclone marker**：篩選 `HPFineNGroups >= 4 and NumReads >= 80`

---

## 權威全欄字典

👉 **`docs/data_specs/20260411_master_dataset欄位字典_01.md`**

---

## 相關

- C++ 原始：[01_significance_summary_schema.md](01_significance_summary_schema.md)
- 衍生特徵：[../07_derived_features/](../07_derived_features/)
- Cross-sample 分析腳本：[../06_workflows/](../06_workflows/)
