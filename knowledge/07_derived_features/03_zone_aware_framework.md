---
id: ism-kb-07-derived-features-zone-aware-framework
name: "Zone-Aware Framework"
description: "基於 AF/LOH/HP 的 5 分區（Z1-Z5）；Zone TP rate 差異真實但 QS 調整 NEGATIVE；characterization only。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "Zone-Aware against docs/concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md"
related_ids:
  - ism-kb-07-derived-features-index
  - ism-kb-09-conclusions-characterization-only
tags: [features, zone-aware, framework, zone, characterization]
canonical_paths: [07_derived_features/03_zone_aware_framework.md]
alias_paths: []
---

# Zone-Aware Framework

- 一句結論：基於 AF × LOH × HP 將 region 分為 5 zone（Z1-Z5）；Zone TP rate 差異真實（characterization）但 QS 調整 NEGATIVE
- 適用對象：理解 TP/FP 分布、Zone-specific 分析
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/research/zone_aware_validation/
  ```

---

## Zone 定義（5 分區）

| Zone | 條件 | TP rate 趨勢 | 說明 |
|------|------|-------------|------|
| **Z1** | High AF (>0.4) + germline HP tag | 高 | Canonical somatic in germline context |
| **Z2** | Low AF + mono-allelic LOH | 中 | Subclonal / LOH somatic |
| **Z3** | Low AF + AF germline contamination | 低 | Germline FP 風險 |
| **Z4** | High AF + LOH | 高 | Clonal LOH somatic |
| **Z5** | Mixed / 難分類 | 變動 | Residual bucket |

**權威文件**：`docs/concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md`

---

## Zone 計算欄位

- `zone_label` (master_dataset)：Z1 / Z2 / Z3 / Z4 / Z5
- `zone_confidence`：zone 分類信心度
- 衍生自：AF, Potential_LOH, HP_Ratio, Coverage_Multiple

---

## H1 / H3 驗證結果

### H1：Zone TP rate 有統計顯著差異 → ✅ POSITIVE
- 7 樣本一致：Z1/Z4 (high) 明顯高於 Z3 (low)
- 效應量：TP rate 差異 ~0.3-0.5

### H3：QS 調整可改善 F1 → ❌ NEGATIVE
- QS simulation 結果：ΔF1 無顯著改善
- 結論：Zone 是**真實現象**，但**不可**用作 F1-filter

---

## 🟡 Characterization Only

**狀態**：Zone-Aware 為 characterization framework，**不可做 variant filter**

**具體結論**：
- Zone TP rate 差異真實（H1 POSITIVE）
- QS 調整無法提升 F1（H3 NEGATIVE）
- Amplicon blacklist pilot：HCC1954 需額外處理

**MEMORY**：`project_zone_aware_framework`

---

## H2009 / HCC1954 特殊情況

- **H1437 TO**：5/6 zone TP rate 0.61-0.94（較高差異）
- **HCC1954**：amplicon artifact 干擾 Z3；AF germline 分層僅 1/7 有效

---

## 典型分析

```python
import pandas as pd

df = pd.read_csv('all_region_rows.tsv.gz', sep='\t', low_memory=False)

# Zone TP rate
zone_tp = df.groupby(['sample', 'mode', 'zone_label'])['truth_label'].apply(
    lambda x: (x == 'TP').mean()
).reset_index()
print(zone_tp.pivot_table(index=['sample', 'mode'], columns='zone_label', values='truth_label'))
```

---

## 相關

- 索引：[00_index.md](00_index.md)
- 權威概念：[../../docs/concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md](../../docs/concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md)
- Research landscape：[../../docs/reports/research_landscape/08_Zone_Aware.md](../../docs/reports/research_landscape/08_Zone_Aware.md)
- Research 路徑：`research/zone_aware_validation/`
