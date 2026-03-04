# 20241224_Phase1_MathUtils_FisherExact_Test

## 測試日期
2025-12-24

## 測試範圍
Phase 1: 基礎統計函式庫
- `MathUtils`: log-factorial, log-binomial, odds ratio, chi-square, Cramér's V, hypergeometric PMF
- `FisherExact`: 2x2 精確檢定, RxC Monte Carlo (Fisher-Freeman-Halton)
- `Stats`: 資料結構 (ContingencyTable2x2, ContingencyTableRxC, FullLabel)

## 測試結果

**全部通過：25/25 測試**

### MathUtilsTest (13 tests)

| 測試名稱 | 結果 | 驗證內容 |
|---------|------|---------|
| LogFactorial_SmallValues | ✅ PASS | n=0~20 與 R lfactorial() 誤差 < 1e-6 |
| LogFactorial_LargeValues | ✅ PASS | n=100, 1000, 10000 Stirling 近似正確 |
| LogFactorial_EdgeCases | ✅ PASS | n=0 → 0, n=-1 → NaN |
| LogBinomial_BasicCases | ✅ PASS | C(5,2)=10, C(10,5)=252 |
| LogBinomial_EdgeCases | ✅ PASS | k>n → -inf |
| OddsRatio_BasicCases | ✅ PASS | Haldane correction 計算正確 |
| OddsRatio_WithZeros | ✅ PASS | 零值格子不造成除零錯誤 |
| LogOddsRatio_Symmetry | ✅ PASS | log(OR(a,b,c,d)) = -log(OR(b,a,d,c)) |
| ChiSquare_IndependentData | ✅ PASS | 獨立資料 χ² = 0 |
| ChiSquare_DependentData | ✅ PASS | 強關聯資料 χ² = 40 |
| CramersV_Range | ✅ PASS | V ∈ [0, 1] |
| CramersV_Reliability | ✅ PASS | 小樣本標記 unreliable |
| HypergeomPMF_BasicCases | ✅ PASS | 機率和 = 1.0 |

### FisherExactTest (8 tests)

| 測試名稱 | 結果 | 驗證內容 |
|---------|------|---------|
| Fisher2x2_LadyTastingTea | ✅ PASS | 經典案例 p = 0.4857143 與 R 一致 |
| Fisher2x2_SignificantAssociation | ✅ PASS | 強關聯 p < 0.001 |
| Fisher2x2_NoAssociation | ✅ PASS | 無關聯 p ≈ 1.0 |
| Fisher2x2_OddsRatio | ✅ PASS | OR > 1 正確計算 |
| Fisher2x2_EmptyTable | ✅ PASS | 空表格返回 p = 1.0 |
| FisherRxC_BasicTest | ✅ PASS | Monte Carlo p ∈ [0, 1] |
| FisherRxC_StrongAssociation | ✅ PASS | 強關聯 p < 0.01 |
| FisherRxC_EarlyStopping | ✅ PASS | 99% CI 提早停止機制運作正常 |

### StatsTest (4 tests)

| 測試名稱 | 結果 | 驗證內容 |
|---------|------|---------|
| ContingencyTable2x2_Totals | ✅ PASS | row/col/grand totals 正確 |
| ContingencyTableRxC_Operations | ✅ PASS | get/set/totals 正確 |
| FullLabel_ComboCode | ✅ PASS | 組合碼編碼正確 |
| FullLabel_TestLabels | ✅ PASS | HP 合併規則正確 |

## R 驗證參考

2x2 Fisher's Exact Test (Lady Tasting Tea):
```r
> fisher.test(matrix(c(3,1,1,3), nrow=2))$p.value
[1] 0.4857143
```

Log-factorial 參考:
```r
> lfactorial(0:10)
[1]  0.0000000  0.0000000  0.6931472  1.7917595  3.1780538  4.7874917
[7]  6.5792512  8.5251614 10.6046029 12.8018275 15.1044126
```

## 結論

Phase 1 基礎統計函式庫實作完成：
1. **MathUtils**: 所有數學函式與 R 結果一致
2. **FisherExact**: 2x2 精確檢定與 R 誤差 < 1e-5，RxC Monte Carlo 正常運作
3. **Stats**: 資料結構設計完成，支援雙層 Label 編碼

可進入 Phase 2 開發。
