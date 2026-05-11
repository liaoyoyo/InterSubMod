---
id: ism-kb-07-derived-features-coverage-multiple
name: "Coverage_Multiple"
description: "CN (copy number) proxy 特徵：NumReads / expected_coverage；r≈0.83 vs 真實 CN；expected_coverage 預設 75.0 有 infra bug（未自動估計）。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "Coverage_Multiple correlation against research/coverage_multiple_validation/"
related_ids:
  - ism-kb-07-derived-features-index
  - ism-kb-05-data-formats-significance-summary-schema
  - ism-kb-05-data-formats-master-dataset-schema
  - ism-kb-07-derived-features-loh-af-methylation
  - ism-kb-10-research-status-blockers-and-risks
tags: [features, coverage, cn, proxy, kde]
canonical_paths: [07_derived_features/02_coverage_multiple.md]
alias_paths: []
---

# Coverage_Multiple

- 一句結論：CN proxy，定義為 `NumReads / expected_coverage`；r=0.831 (Paired) / 0.827 (TO) vs 真實 CN
- 適用對象：CN-aware 分析、LOH Zone-Aware 輸入
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # Coverage_Multiple 分布
  awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="Coverage_Multiple") c=i} NR>1{print $c}' \
    output/canonical/HCC1395/paired_full/*/significance_summary.csv | \
    python3 -c "import sys; vals=[float(x) for x in sys.stdin if x.strip()]; print(f'Median: {sorted(vals)[len(vals)//2]:.2f}')"
  ```

---

## 定義

```
Coverage_Multiple = NumReads / expected_coverage
```

**意義**：region 實際讀取數與期望二倍體覆蓋度的比值，近似 CN

**取值**：
- `~1.0`：CN = 2（正常二倍體）
- `~0.5`：CN = 1（LOH / 單倍體）
- `>2.0`：CN > 4（amplification）
- `~0`：CN = 0（deletion）

---

## `expected_coverage` 參數

**CLI**：`--expected-coverage`（預設 `0.0` = 自動估計）
**位置**：`include/utils/ArgParser.hpp:122-124`

**自動估計方法**：KDE mode（peak density 估計）

---

## ⚠️ 已知 Infra Bug（2026-04-20）

**問題**：master dataset 全 7 樣本共用 default `75.0`（而非自動估計）
**影響**：
- HCC1395 (>75x)：不受影響或輕微影響
- COLO829 (9x)：**嚴重影響**，Coverage_Multiple 全部偏低
- 其他中等覆蓋度樣本：部分影響

**狀態**：需 `/cpp-change` 修 KDE 啟用並重跑
**MEMORY**：`project_expected_coverage_baseline_bug`

**如何確認中招**：
```bash
# 看 run 的 expected_coverage log
grep -i "expected_coverage" output/canonical/<sample>/paired_full/*/logs/*.log
```

---

## KDE Fix 驗證（2026-04-20）

- **結論**：KDE auto-estimation 機制正確
- **實測**：精度從 6.2% → 43.8%（HCC1395）
- **commit**：`374fad4`, `12d9b3e`（見 research/kde_fix_validation/）

**下游量化**（2026-04-20）：
- H-CN1 PARTIAL POSITIVE
- Z3 scale-invariant 不變
- COLO829 9 筆衝擊待下 cycle 重跑

---

## Correlation 與 CN

| Pipeline | r vs CN | 來源 |
|----------|---------|------|
| Paired | **0.831** | `research/coverage_multiple_validation/` |
| TO | **0.827** | 同上 |

**GC bias 驗證**：delta r = -0.0002（全 NO-GO；無 GC 衝擊）
**Zone 排除驗證**：全不可行

---

## 欄位位置

- **在 significance_summary.csv**：L 群「品質評估」的 `Coverage_Multiple`
- **Python 衍生**：`coverage_quartile`, `log10_numreads` 等在 master_dataset 116 欄中

---

## 典型分析

```python
import pandas as pd
import numpy as np

df = pd.read_csv('significance_summary.csv')

# LOH 候選（CN ≈ 1）
loh_candidates = df.query('Coverage_Multiple.between(0.4, 0.7)', engine='python')

# Amplification 候選（CN > 4）
amp_candidates = df.query('Coverage_Multiple > 2.0')

# 分布
print(f"Median Coverage_Multiple: {df['Coverage_Multiple'].median():.2f}")
print(f"IQR: {df['Coverage_Multiple'].quantile(0.25):.2f} - {df['Coverage_Multiple'].quantile(0.75):.2f}")
```

---

## 相關

- 索引：[00_index.md](00_index.md)
- LOH × AF：[04_loh_af_methylation.md](04_loh_af_methylation.md)
- CLI 參數：[../04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md)
- KDE Fix research：`research/kde_fix_validation/`
