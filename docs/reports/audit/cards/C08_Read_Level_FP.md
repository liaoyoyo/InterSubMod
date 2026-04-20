# C08: Read-level LOSO AUC=0.721 但 FP removal=0%（結論 8）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 8
- **來源文件**: `04_暫停判定與重評估.md` / `06_結論穩定性審查.md:175-206`
- **穩定度評分**: ⭐3 需注意
- **狀態**: CONDITIONAL NO-GO（最可能反轉）
- **所屬證據鏈**: 🔗5 / 🔗6（暫停判定子鏈）
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ within_dom_alt_frac AUC=0.679 / LOSO=0.721 / TP loss ≤2% 跨文件一致
- ⚠️ AUC 0.721 vs C07 0.638 門檻不一致 → **P2-B**

## D2 方法論健康度
- ✅ LOSO 跨樣本驗證
- ⚠️ **within_dom_alt_frac 依賴 HP → P0-C 修正後可能突破**
- ⚠️ 高純度 cell line 根因，低純度樣本可能更佳 → 未驗證
- ⚠️ 安全約束 TP loss ≤2% 過嚴？未測 5% 放寬
- ❌ Bootstrap CI 缺
- ❌ 多個閾值/特徵掃描無 FDR

## D3 證據鏈
- **依賴結論**: C01, C03, C07
- **被依賴結論**: 無（暫停中）
- **鏈完整度**: ⚠️ 依賴 P0-C 修正

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: **需**（P0-C 完成 + 低純度樣本）

## D5 統計嚴謹度
- **Effect size**: LOSO AUC 0.721 中等強度
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ✅ LOSO 設計
- **Multiple testing**: ❌ 閾值/特徵掃描無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: Read-level classification 在 cancer variant call 少見，但 epigenetic signature 文獻支持（knowledge 02/03）
- **挑戰文獻**: 無
- **缺口**: 低純度樣本（purity 30-70%）

## 修正建議
- **P0-C**: haplotag 修正後 within_dom_alt_frac 重測
- **P2**: 低純度樣本 + 安全約束放寬至 5%
- **P1-B**: 補 bootstrap CI + FDR
- **對應 R-id**: R-04
- **對應 Q-id**: Q-06

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: LOSO AUC=0.721 但 FP removal=0% 已確認安全約束下無效；P0-C ReadParser 修正後可能翻轉（等 rebuild 後再評估）。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**🟡 暫停判定中 — 最具翻轉潛力，必須等 P0-C 修正後再評估。**
