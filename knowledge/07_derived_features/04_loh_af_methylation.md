---
id: ism-kb-07-derived-features-loh-af-methylation
name: "LOH × AF × Methylation"
description: "LOH 區域的 AF→NGroups 關係：Inter +0.705 (7/7 p<10^-39)；最強 positive finding；雙重證據鏈；文獻空白。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "LOH × AF × Methylation against research/loh_subclone_af_methylation/ and MEMORY"
related_ids:
  - ism-kb-07-derived-features-index
  - ism-kb-09-conclusions-positive-findings
  - ism-kb-07-derived-features-coverage-multiple
tags: [features, loh, af, methylation, positive, subclone]
canonical_paths: [07_derived_features/04_loh_af_methylation.md]
alias_paths: []
---

# LOH × AF × Methylation

- 一句結論：🟢 最強 positive finding：Inter AF→NGroups **+0.705** (7/7 樣本 p<10^-39 to 10^-65)；雙重證據鏈；文獻空白
- 適用對象：主要 positive findings 理解、論文核心數據
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/research/loh_subclone_af/
  ls /big7_disk/liaoyoyo2001/InterSubMod/research/loh_subclone_af_paired/
  ```

---

## 核心發現

### Inter AF → NGroups：+0.705 相關性
- **7 樣本一致**：全 7 樣本正相關
- **p-value**：p<10^-39 到 10^-65
- **方向**：LOH region 內，AF 升高 → HPFineNGroups 數值升高

### 雙重證據鏈
1. **Paired_full 證據**：7/7 positive
2. **TO 證據**：7/7 positive
3. **Direction consistency**：100%
4. **文獻空白**：無先行研究描述此關聯

---

## 生物學意義

**假說**：
- LOH region 中，AF 代表 tumor 純度或 subclone 比例
- 高 AF → 更多 somatic reads → 更多 subclone fine-grained 可分性
- 甲基化模式展現 subclone heterogeneity

**影響**：
- Read-level epigenetic characterization 核心價值的體現
- ISM 定位：non-variant-filter 用途的學術重要性

---

## 🟢 Positive（Subclone Characterization）

**狀態**：
- ✅ 真實信號（7/7 跨樣本驗證）
- ✅ 雙 pipeline 一致（paired + TO）
- ⚠️ 不是 variant filter（ISM F1 主表已 locked +0.0112）

**定位**：做 read-level subclonal characterization 的 gold finding

---

## Research 歷史

- Research 路徑：`research/loh_subclone_af/`（TO）, `research/loh_subclone_af_paired/`（Paired）
- 2 個獨立研究專案（Paired + TO）
- MEMORY：`project_loh_subclone_af_methylation_positive`

---

## 欄位依賴

- `AF` / `AlleleDelta`（caller 輸出）
- `HPFineNGroups`（C++ 計算）
- `Potential_LOH`（C++ 計算）
- Zone-Aware labels（Python 衍生）

---

## 典型分析

```python
import pandas as pd
from scipy.stats import pearsonr

df = pd.read_csv('all_region_rows.tsv.gz', sep='\t', low_memory=False)

# Inter AF → NGroups（LOH region 內）
for sample in df['sample'].unique():
    sub = df.query(f"sample == '{sample}' and mode == 'paired_full' and Potential_LOH == 1")
    if len(sub) < 100: continue
    r, p = pearsonr(sub['AlleleDelta'], sub['HPFineNGroups'])
    print(f"{sample}: r={r:.3f}, p={p:.2e}")
# 預期：全 7 樣本 r > 0.5, p < 10^-39
```

---

## 相關

- 索引：[00_index.md](00_index.md)
- HPFineNGroups：[01_hpfinengroups.md](01_hpfinengroups.md)
- Coverage：[02_coverage_multiple.md](02_coverage_multiple.md)
- Positive findings：[../09_conclusions/01_positive_findings.md](../09_conclusions/01_positive_findings.md)
- LOH 統整：[../../docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md](../../docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md)
