# Step 8: 跨樣本 Coverage_Category PCA — baseline fix 下游影響

## 背景

原 O1-O10 Fig 10 結論：「COLO829 在 PCA 中孤立 → Depth 是第一驅動因素」
疑慮：stale 75× 讓 COLO829 CovM 壓縮至 0.39，造成 CNV_Loss 泛濫，可能是 PCA 孤立的 artifact。

## 特徵向量（per-sample, 6-D Coverage_Category %）

### COLO829 vector

| Category | Stale % | Fixed % | Δpp |
|----------|--------:|--------:|----:|
| CNV_Loss | 79.43 | 3.93 | -75.50 |
| Low | 19.83 | 23.67 | +3.85 |
| Normal | 0.64 | 42.74 | +42.10 |
| Elevated | 0.08 | 21.54 | +21.46 |
| CNV_Gain | 0.02 | 7.23 | +7.22 |
| High_Copy | 0.00 | 0.88 | +0.88 |

### PCA 孤立度對比

| 指標 | Stale 75× | KDE-fixed | 變化 |
|------|----------:|----------:|------|
| COLO829 NN-dist | 81.08 | 12.38 | -84.7% |
| 其他樣本 NN-dist 中位 | 6.83 | 7.30 | — |
| COLO829 isolation ratio | 11.88× | 1.70× | — |
| PC1/PC2 var explained (stale) | 78.2% + 18.2% | — | — |
| PC1/PC2 var explained (fixed) | — | 81.7% + 16.2% | — |

**Verdict**: COLO829 no longer extreme outlier after KDE fix

## 結論

- 「Depth 是第一驅動」屬於 stale baseline artifact
- KDE fix 後 COLO829 Category 分佈已向其他樣本靠攏
- 原 O1-O10 Fig 10 結論需**降級**：從「depth driver」改為「stale baseline artifact」
