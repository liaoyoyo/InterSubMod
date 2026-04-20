# C06: O13 跨區域 correlation = shared read confound（結論 6）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 6
- **來源文件**: `05_證據鏈總覽.md` / `06_結論穩定性審查.md:152-160`
- **穩定度評分**: ⭐5 堅固
- **狀態**: NEGATIVE（方向正式關閉）
- **所屬證據鏈**: 🔗3（TO AUC ceiling, shared read 子鏈）
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ median 36 vs 21 (FP-FP vs TP-TP shared reads) 一致
- ✅ OLS p=0.464, d=-0.071 / matching n=500, p=0.403 一致

## D2 方法論健康度
- ✅ 三種方法一致（分層 6 strata + OLS + matching）
- ✅ **confound 被識別並修正**（典範）
- ⚠️ OLS 是 pooled 還是 within-stratum？文件未明示 → **R-05**
- ✅ 50kb 距離閾值合理（ISM single-region 架構）
- ❌ 無 bootstrap CI、無 FDR

## D3 證據鏈
- **依賴結論**: C01, C04（方法論沿用）
- **被依賴結論**: 無（已關閉方向）
- **鏈完整度**: ✅ 完整但為 dead-end branch

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: 無

## D5 統計嚴謹度
- **Effect size**: OLS d=-0.071 = 接近 null
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ✅ n=500 matching + 6 strata
- **Multiple testing**: ⚠️ 多 stratum 掃描無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: Sample overlap confound 在 GWAS / meta-analysis 為已知問題
- **挑戰文獻**: 無
- **缺口**: 無

## 修正建議
- **R-05**: 確認 OLS 類型
- **R-04**: 若接受 NO-GO 結論重審方案，納入 shared read confound 統一檢查（C04-C06 群組）
- **對應 R-id**: R-04, R-05
- **對應 Q-id**: Q-06

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: shared read count confound 已於 2026-04-02 透過三種方法一致驗證（分層後信號消失）。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**✅ 方向已關閉，方法論穩固 — 僅需補 OLS 類型確認。**
