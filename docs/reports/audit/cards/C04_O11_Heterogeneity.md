# C04: O11 within-group heterogeneity = n_reads confound（結論 4）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 4
- **來源文件**: `05_證據鏈總覽.md` / `06_結論穩定性審查.md:115-137`
- **穩定度評分**: ⭐5 堅固
- **狀態**: NEGATIVE
- **所屬證據鏈**: 🔗3（TO AUC ceiling, n_reads confound 子鏈）
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ AUC 0.845（原）vs 0.560（matched）vs 0.509-0.578（residualized）數字一致
- ✅ Spearman r=0.79 (n_reads vs epipolymorphism) 一致
- ✅ TP median 157 vs FP 92 (n_reads) 一致

## D2 方法論健康度
- ✅ 無 pooled OLS trap（三種方法：分層 + matching + residualization 交叉驗證）
- ✅ 無 L2 collider（n_reads 是 ancestor 非 collider）
- ✅ **confound 被主動識別並修正**（典範方法論）
- ✅ 無 spatial auto-correlation
- ⚠️ 殘差化使用的 OLS 類型未明示 within-group vs pooled → **R-05 全面回溯時應標註**
- ❌ 無 bootstrap CI（single-point）→ 低優先補

## D3 證據鏈
- **依賴結論**: C01
- **被依賴結論**: C05 (O12 L2 collider)、C06 (O13 shared read)
- **鏈完整度**: ✅ 完整（三種方法一致）

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無（epipolymorphism 不用 CovM）
- **重跑必要性**: 無

## D5 統計嚴謹度
- **Effect size**: 原 AUC 0.845 → matched 0.560 Δ=-0.285（confound 幅度巨大）
- **CI 覆蓋率**: ❌ 無 bootstrap CI
- **Power 評估**: ✅ n=148K regions
- **Multiple testing**: ⚠️ 6 特徵殘差化掃描無 FDR → **P1-B**

## D6 知識庫交叉驗證
- **支持文獻**: 知識庫 05（統計陷阱類）無直接 counterpart；epipolymorphism 在 epigenetics 文獻常見（Fang 2012 scEPI）
- **挑戰文獻**: 無
- **缺口**: 單樣本 coverage matching 對人類腫瘤的普適性

## 修正建議
- **R-05**: 確認殘差化用 within-group OLS（否則 C04 與 C17 C15 一樣有 pooled OLS 疑慮）
- **P1-B**: 6 特徵殘差化套 BH-FDR + bootstrap CI
- **對應 R-id**: R-05
- **對應 Q-id**: 無

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: n_reads confound 已於 2026-04-02 透過分層 + matching + residualization 三種方法交叉驗證識別並修正（AUC 0.845→0.560）。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**✅ 方法論典範 — confound 被成功識別並修正；唯一需確認 OLS 類型。**
