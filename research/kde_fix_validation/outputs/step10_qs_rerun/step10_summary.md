# Step 10: Paired-mode QS rerun — P1/P2/P3 部分驗證

## 重要前提

本步驟僅涵蓋 **paired 模式**（14 combos rerun 範圍）。
原 P1/P2/P3 結論基於 **TO 模式** 觀察；paired 結果僅為方向性驗證，
TO 模式需另行 rerun 才能完成最終 8/8 驗收。

## P1 COLO829 QS × Coverage_Category（CNV_Loss penalty 驗證）

### COLO829 QS median by Category (paired)

| Category | Stale % | Stale QS | Fixed_pileup % | Fixed_pileup QS | Δ% | ΔQS |
|----------|--------:|---------:|---------------:|----------------:|----:|-----:|
| CNV_Loss | 79.43% | 50.0 | 3.94% | 50.0 | -75.50pp | 0.0 |
| Low | 19.83% | 75.0 | 23.72% | 65.0 | +3.89pp | -10.0 |
| Normal | 0.64% | 100.0 | 42.73% | 80.0 | +42.09pp | -20.0 |
| Elevated | 0.08% | 100.0 | 21.47% | 90.0 | +21.39pp | -10.0 |
| CNV_Gain | 0.02% | 100.0 | 7.22% | 90.0 | +7.20pp | -10.0 |
| High_Copy | 0.00% | nan | 0.93% | 80.0 | +0.93pp | — |

## P2/P3 跨樣本 QS median

| Sample | Stale_paired median | Fixed_pileup median | Fixed_full median | Δ(pileup-stale) |
|--------|--------------------:|--------------------:|------------------:|-----------------:|
| HCC1395 | 85.0 | 100.0 | 100.0 | +15.0 |
| HCC1395_DORADO | 85.0 | 100.0 | 100.0 | +15.0 |
| HCC1937 | 75.0 | 100.0 | 100.0 | +25.0 |
| HCC1954 | 100.0 | 100.0 | 100.0 | +0.0 |
| H2009 | 100.0 | 100.0 | 100.0 | +0.0 |
| H1437 | 100.0 | 100.0 | 100.0 | +0.0 |
| COLO829 (COLO829=melanoma) | 60.0 | 80.0 | 80.0 | +20.0 |

## 初步結論（paired only；TO 待驗證）

- **P1** COLO829 paired QS：stale median 60.0 → fixed 80.0（Δ=+20.0）
  - ✅ QS **明顯回升**，CNV_Loss penalty 影響確認
- **P2** COLO829 vs 其他樣本 paired QS median gap：80.0 vs 100.0（差 20.0）
  - ⚠️ 仍有顯著 gap（但小於 stale），ISM-for-melanoma 假說 paired 下部分成立
- **P3** 跨樣本 QS ranking paired 模式下 COLO829 仍最低，但 gap 縮小（stale→fixed）

### TO 模式 follow-up 必要性

原 P1/P2/P3 結論主要基於 TO 模式（QS median 35-60 vs 75+）。
需要 TO 模式 rerun 才能做最終判定。Paired 結果顯示方向正確（QS 回升），
但量級差異需 TO 實測確認。
