# C22: Zone-Aware F1 (H2) NEGATIVE + Zone Characterization CONFIRMED

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 22（相當於 INDEX 結論 20 + 21 雙重）
- **來源文件**: `08_Zone_Aware.md` / `06_結論穩定性審查.md`
- **穩定度評分**: ⭐4（F1 H2）+ ⭐3（Characterization）
- **狀態**: F1 NEGATIVE / Characterization CONFIRMED
- **所屬證據鏈**: 🔗10
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ Z1-Z5 zone 定義一致
- ✅ TO Z3 最低 TP rate 0.608 一致
- ✅ **`20260417_Zone_Aware_Confidence_Framework_01.md:49,73` r=0.997 已修正為 r=0.831（2026-04-20 P1-A）**；原誤標註記保留**（Coverage_Multiple）→ **P1-A 已標**
- ⚠️ 「Characterization」應用場景模糊（C-STRAT-3）

## D2 方法論健康度
- ✅ H1-H5 驗證矩陣完整
- ❌ **Zone 定義受 CovM bug 汙染**（zone 邊界依賴 CovM 分層）
- ⚠️ QS 調整 NEGATIVE 已確立，但 zone 邊界變動會影響 characterization 結論
- ❌ 缺 bootstrap CI / FDR

## D3 證據鏈
- **依賴結論**: C17, C18, C19, C20, C21
- **被依賴結論**: Phase 2 A+D 的 zone-aware Normal BAM 設計
- **鏈完整度**: ⚠️ CovM bug 影響整條

## D4 數據信任度
- **dataset 版本**: Haplotag v5 + Zone-Aware run
- **CovM bug 影響**: 🔴 **Zone 定義依賴 CovM baseline → bug 修正後需重定義**
- **重跑必要性**: **必需**（R-01 後）

## D5 統計嚴謹度
- **Effect size**: Zone TP rate 差異真實（Z3 vs Z1 差 ~0.15）
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ✅
- **Multiple testing**: ❌ 5 zone × 多指標無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: Variant confidence zone / region stratification 在 GATK/DeepVariant 有類似概念
- **挑戰文獻**: 無
- **缺口**: Characterization 的具體 downstream 應用

## 修正建議
- **P0-A**: R-01 CovM bug 修正 + Zone 重定義
- **P1-A**: 修正 r=0.997 誤標
- **P1-B**: Zone TP rate 差異補 bootstrap + FDR
- **P2-B**: Characterization 的具體功能清單化（用戶主軸定位下）→ 交由 `cross_cutting/Characterization_Functions.md`
- **對應 R-id**: R-01, R-02
- **對應 Q-id**: Q-01, Q-02, Q-03

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: Zone-Aware Framework NEGATIVE-for-QS / POSITIVE-for-Characterization 雙軸結論；QS NEGATIVE 穩固（跨 7 樣本驗證）。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**🔴 CovM bug 影響 zone 定義 — 修正後整個 zone-aware framework 需重跑；characterization 定位與用戶主軸一致。**
