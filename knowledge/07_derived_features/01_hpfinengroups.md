---
id: ism-kb-07-derived-features-hpfinengroups
name: "HPFineNGroups"
description: "Fine-grained PERMANOVA 有效群組數（HP1-G/HP1-S/HP2-G/HP2-S，各 ≥3 reads）；N≥4+NR≥80 TP rate 92.81%；characterization only，不可做 variant filter。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "HPFineNGroups against research/F_hpfinengroups_deepening/ and docs/reports/research_landscape/09_Part_B.md"
related_ids:
  - ism-kb-07-derived-features-index
  - ism-kb-04-parameters-statistical-methods
  - ism-kb-05-data-formats-significance-summary-schema
  - ism-kb-09-conclusions-characterization-only
  - ism-kb-09-conclusions-positive-findings
tags: [features, hpfinengroups, characterization, subclone, marker]
canonical_paths: [07_derived_features/01_hpfinengroups.md]
alias_paths: []
---

# HPFineNGroups

- 一句結論：四分群（HP1-G/HP1-S/HP2-G/HP2-S）有效組數；subclone marker 用途（🟡 characterization only）；**不可做 variant filter**
- 適用對象：Subclone 特徵化研究、read-level 分析
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 從 significance_summary.csv 抓 HPFineNGroups 分布
  awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="HPFineNGroups") c=i} NR>1{print $c}' \
    output/canonical/HCC1395/paired_full/*/significance_summary.csv | sort | uniq -c
  ```

---

## 定義

**HPFineNGroups**：在 Fine-grained PERMANOVA 中，有效的群組數量

**四個可能的 groups**：
- `HP1-G`：germline HP1（HP:i:1）
- `HP1-S`：somatic HP1（HP:i:11）
- `HP2-G`：germline HP2（HP:i:2）
- `HP2-S`：somatic HP2（HP:i:21）

**有效性條件**：每 group 需 ≥3 reads

**取值**：0, 1, 2, 3, 4

---

## Canonical Filter（Part B 升級，2026-04）

**條件**：`HPFineNGroups >= 4 AND NumReads >= 80 AND AF < 0.4`

**效果**（HCC1395 paired_full）：
- TP rate: **92.81%**（vs 舊 89.12%）
- 提升 +3.7 pp

**研究路徑**：`research/F_hpfinengroups_deepening/`

---

## 🟡 Characterization Only 警告

**狀態**：確認為**真實** subclone marker，但**不可**用作 variant filter

**原因**：
- 若用 `HPFineNGroups <= N` 作 FP filter → 產生**負 F1 增益**
- 原因：subclone 結構在 TP 中也有，filter 會傷害 recall

**結論來源**：
- [../09_conclusions/02_characterization_only.md](../09_conclusions/02_characterization_only.md)
- [../../docs/reports/research_landscape/09_Part_B.md](../../docs/reports/research_landscape/09_Part_B.md)

---

## ⚠️ 2026-04-21 警告（flag=on 下的 artifact）

**現象**：當 `--germline-hp-only=on` 時，HPFineNGroups N≥3 **完全消失**

**可能原因**：
1. Somatic HP tag (HP:i:11/21/33) 本身是 LongPhase 人工分組的 artifact
2. 關閉 somatic HP tag → fine-grained 降為純 germline → N 降為 1-2

**影響**：subclone marker 結論需在 flag=on 下重新驗證

**MEMORY 記錄**：
- `project_readparser_germline_hp_only_phase1_negative`
- `project_hpfinengroups_subclone_marker`（POSITIVE 但需重驗證）

---

## 欄位位置

- **在 significance_summary.csv**：G 群「Stage 2: HP Fine-Grained」，4 欄中的 `HPFineNGroups`
- **計算**：`src/core/RegionProcessor.cpp`（Fine-grained PERMANOVA 階段）

---

## 典型分析

```python
import pandas as pd

df = pd.read_csv('significance_summary.csv')

# Subclone marker 分布
print(df.groupby('HPFineNGroups').size())

# Canonical filter 適用 region
candidates = df.query('HPFineNGroups >= 4 and NumReads >= 80 and AlleleDelta < 0.4')
print(f"Candidates: {len(candidates)}")

# TP rate（需 join truth_label）
master = pd.read_csv('all_region_rows.tsv.gz', sep='\t', low_memory=False)
tp_rate = master.query('HPFineNGroups >= 4 and NumReads >= 80')['truth_label'].eq('TP').mean()
print(f"TP rate: {tp_rate:.4f}")  # 預期 ~0.928 (paired_full HCC1395)
```

---

## 7 樣本一致性（精確化）

2026-04-18 Part B 全量驗證實際結果：
- **5/7 POSITIVE**（3 medium + 2 small effect size）
- **2/7 特殊**：COLO829（out-of-scope，無 methylation）、H2009（ceiling effect，caller 已飽和）
- ⚠️ 先前「7/7 direction consistency」為**過度宣稱**；精確描述見 `docs/reports/research_landscape/09_Part_B.md` §3.3

---

## ⚠️ 機制更正（2026-04-22）

**過去描述**：HPFineNGroups 是 methylation subclone marker
**更新理解**：HPFineNGroups 實為 `{HP1-G, HP1-S, HP2-G, HP2-S}` **四 bucket 的 occupancy count**，反映 LOH 區域內 phasing 可分性，**不直接代表 methylation subclone**
**真實機制**：LOH-constrained phasing discovery（見 MEMORY `project_loh_constrained_phasing_discovery` 若存在）
**影響**：過去基於「subclone marker」語義的推論需重新審視；TP rate 差異仍真實但 interpretation 需修正

---

## 相關

- 索引：[00_index.md](00_index.md)
- 統計方法：[../04_parameters/03_statistical_methods.md](../04_parameters/03_statistical_methods.md)
- Characterization only：[../09_conclusions/02_characterization_only.md](../09_conclusions/02_characterization_only.md)
- 權威：[../../docs/reports/research_landscape/09_Part_B.md](../../docs/reports/research_landscape/09_Part_B.md)
- Research 路徑：`research/F_hpfinengroups_deepening/`
