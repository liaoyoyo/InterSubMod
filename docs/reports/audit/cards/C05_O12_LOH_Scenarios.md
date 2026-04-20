# C05: O12 TO LOH 三場景不可區分（結論 5）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 5
- **來源文件**: `05_證據鏈總覽.md` / `06_結論穩定性審查.md:140-148`
- **穩定度評分**: ⭐4 穩固
- **狀態**: NEGATIVE
- **所屬證據鏈**: 🔗4（L2 collider bias 範例）
- **所屬樣本**: 7/7 (175K rows × 22 features × 4-level confound control)

## D1 內部一致性
- ✅ 22 特徵、4-level confound control 數字一致
- ✅ L3 AUC < 0.59 閾值一致

## D2 方法論健康度
- ✅ **L2 collider bias 被識別並修正**（典範）
- ✅ AlleleDelta=AF confound 被識別
- ⚠️ **P0-C haplotag 汙染影響**：LOH 區域 HP tag 汙染最嚴重（結論 9 self-phasing 62% artifact）→ 修正後可能微改變
- ⚠️ 22 特徵是否涵蓋所有鑑別維度？未測 spatial（CGI/shore, enhancer）
- ❌ Pooled OLS 疑慮：AlleleDelta residualize on AF 是否 within-group（按 L3 AF bin）？→ **R-05 回溯**
- ❌ 無 bootstrap CI、無 FDR

## D3 證據鏈
- **依賴結論**: C01, C04（confound 控制方法）
- **被依賴結論**: C15（LOH 區域甲基化全面失效）
- **鏈完整度**: ⚠️ haplotag 修正前為「中間狀態」

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無（O12 不用 CovM）
- **重跑必要性**: **需**（P0-C 修正後）

## D5 統計嚴謹度
- **Effect size**: L3 AUC < 0.59 = 無效
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ✅ n=175K
- **Multiple testing**: ❌ 22 特徵 × 4 層 = 88 tests 無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: L2 collider bias 在 epidemiology 為經典陷阱（Hernán & Robins）
- **挑戰文獻**: 無
- **缺口**: 無

## 修正建議
- **P0-C**: haplotag 修正後 L2 AlleleDelta 重測
- **P1-B**: 88 tests 套 BH-FDR + bootstrap
- **R-05**: 確認 AlleleDelta residualize on AF 使用 within-group OLS
- **對應 R-id**: R-04, R-05
- **對應 Q-id**: Q-06

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: AlleleDelta=AF confound + L2 collider bias 已於 2026-04-02 透過分層 + AF-bin 交叉驗證識別（L3 全<0.59）；within-group OLS 可於 R-05 盤點後補做，但不影響 NEGATIVE 結論。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**⚠️ 方法論精神正確但需明確 OLS 類型 + haplotag 修正後重測。**
