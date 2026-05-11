# Step 7: M2 mask 重算 — KDE baseline 修正下游影響

## 背景

M2 mask: `Coverage_Multiple < 0.5 OR > 2.0`（CNV 異常區域遮罩）
原 `20260408_TO_LOH額外研究` 報告 COLO829 91.7% 被 M2 排除，
歸因於 stale baseline 75× vs true 29×（ratio 2.59）。
本步驟量化 KDE fix 後各樣本 M2 排除率變化。

## M2 排除率對比表（paired 模式）

| Sample | stale_paired M2% | fixed_paired_pileup M2% | fixed_paired_full M2% | Δ (pileup−stale) |
|--------|-----------------:|------------------------:|----------------------:|-----------------:|
| HCC1395 | 10.49% | 12.71% | 13.66% | +2.22pp |
| HCC1395_DORADO | 8.30% | 16.49% | 14.04% | +8.19pp |
| HCC1937 | 18.95% | 12.55% | 12.05% | -6.40pp |
| HCC1954 | 6.32% | 6.98% | 5.54% | +0.66pp |
| H2009 | 8.97% | 8.41% | 8.23% | -0.55pp |
| H1437 | 2.66% | 2.15% | 2.07% | -0.51pp |
| COLO829 | 79.40% | 4.87% | 4.75% | -74.53pp |

## 關鍵觀察

### COLO829 戲劇性下降

- **Stale (baseline 75×)**: 79.40% 被 M2 排除
  - 低界 (CovM<0.5): 79.40%
  - 高界 (CovM>2.0): 0.00%
- **Fixed (baseline 29× per-sample KDE)**: 4.87% 被 M2 排除
  - 低界 (CovM<0.5): 3.94%
  - 高界 (CovM>2.0): 0.93%
- **變化**: -74.53pp

### 結論對既有 M2 觀察的影響

- 原 `20260408` 報告 M2 「偏向 COLO829（91.7% 被排除）」的結論**屬於 baseline artifact**
- KDE fix 後 COLO829 M2 排除率與其他樣本相近
- **原 M2 做為 `唯一合理性價比遮罩`** 的定性結論仍成立（不依賴 sample-specific 偏差）
- 但**「M2 偏向 COLO829」這句話需要撤回**
