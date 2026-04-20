# C14: QS TO AUC=0.497（結論 14）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 14
- **來源文件**: `03_ISM分析價值界定.md` / `06_結論穩定性審查.md:287-295`
- **穩定度評分**: ⭐4 穩固
- **狀態**: NEGATIVE
- **所屬證據鏈**: 🔗3（TO AUC ceiling）
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ AUC=0.497 一致
- ✅ 非 LOH QS AUC=0.549（LOH 解釋 35%）一致

## D2 方法論健康度
- ✅ 多樣本驗證
- ⚠️ QS 設計本身的 LOH penalty 反向 → **需獨立重設計**
- ✅ 無 confound 疑慮（AUC 接近 random = null result）
- N/A bootstrap/FDR（null 結果）

## D3 證據鏈
- **依賴結論**: C01, C07（germline leak 根因）, C09（LOH penalty 反向）
- **被依賴結論**: QS 重設計（未啟動）
- **鏈完整度**: ✅ 根因歸因完整

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: 無（QS 已廢棄）

## D5 統計嚴謹度
- **Effect size**: AUC 0.497 = null
- **CI 覆蓋率**: N/A（null）
- **Power 評估**: ✅
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: Variant QS 在 germline leak 場景失效為已知機制
- **挑戰文獻**: 無
- **缺口**: 無

## 修正建議
- **R-04**: 接受 NEGATIVE 現狀
- **對應 R-id**: R-04
- **對應 Q-id**: 無

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: AUC=0.497 隨機已確認 TO QS 失效；根因歸因到 C07（germline leak）閉環完整。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**✅ NEGATIVE 穩固 — 根因（germline leak）已歸因到 C07，結論閉環。**
