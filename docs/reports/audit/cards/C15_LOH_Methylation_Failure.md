# C15: LOH 區域內甲基化全面失效（補充結論 15）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 15（補充）
- **來源文件**: `07_LOH_CN_AF_研究總整理.md` / `06_結論穩定性審查.md:298-308`
- **穩定度評分**: ⭐3 需注意
- **狀態**: NEGATIVE with caveat
- **所屬證據鏈**: 🔗8 LOH 子鏈
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ CramersV AUC~0.53, PairwiseMedianDist~0.54, HPMergedDelta~0.50 一致
- ❌ **措辭「全面失效」過強**：06_審查 已建議改為「依賴 HP 分群的甲基化特徵在 LOH 區域天然無效」→ **P2-A 跨文件統一**

## D2 方法論健康度
- ⚠️ LOH 區域本質為 HP 分群退化，依賴 HP 特徵 by design 無效
- ❌ **未測試 HP-independent 甲基化特徵**（region-level entropy, mean methylation）在 LOH 區域
- ❌ 無 bootstrap CI、無 FDR
- ⚠️ **Pooled OLS 疑慮**（C-STAT-1 ref）：若有 residualize on LOH 狀態，需 within-LOH OLS → **R-05**

## D3 證據鏈
- **依賴結論**: C05（O12 L2 collider）、C09（Self-phasing）
- **被依賴結論**: C17（LOH Subclone 的對比基礎）
- **鏈完整度**: ⚠️ HP-independent 測試空白

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無（LOH 區域甲基化不用 CovM）
- **重跑必要性**: **需**（HP-independent 特徵）

## D5 統計嚴謹度
- **Effect size**: AUC ~0.50-0.54 = 接近 null
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ✅
- **Multiple testing**: ❌ 3+ 特徵無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: LOH 區域 allelic imbalance 影響 HP-based methylation 分析（knowledge 03, 07）
- **挑戰文獻**: epiTRACERx / EVOFLUx 在 LOH 區域仍能擷取 subclone signal（需查確認是否 HP-based）
- **缺口**: HP-independent 甲基化在 LOH 區域表現的文獻

## 修正建議
- **P2-A**: 措辭修正「依賴 HP 的甲基化特徵在 LOH 區域天然無效（non-failure, by-design）」
- **P2**: 測試 HP-independent 特徵（region entropy, mean meth）在 LOH 區域
- **R-05**: 確認 OLS residualization 類型
- **對應 R-id**: R-05
- **對應 Q-id**: Q-10

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: HP-independent 驗證空白已註記；Pooled OLS 疑慮由 R-05 全面盤點處理；NEGATIVE 結論（AlleleDelta=AF confound）不動搖。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## P0-B 補註（2026-04-21）— O15 Scripts 已為 Stratified 設計

**Script inventory**：
- `scripts/analysis/build_observation_O15_loh_zone_metrics_hcc1395.py`（1245 行）
- `scripts/analysis/build_observation_O15b_loh_zone_metrics_cross_sample.py`（1478 行）

**驗證結果**：
- AUC / Cohen's d 計算：**per mode × loh_group × truth_label stratified**（line 960-1001），效果量與 CI 均 within-stratum 計算
- pearsonr 呼叫：**per-zone（LOH-inside / LOH-outside）**（line 468），非 pooled
- `pooled = np.mean(ds)` in O15b（line 807）：僅為 forest plot 視覺化的「各 metric 平均 Cohen's d」，非 OLS residualization
- **無 LinearRegression / sm.OLS / residualize 呼叫**

**結論**：C15 現有 NEGATIVE 結論基於 by-design HP 退化，分析腳本本身已為 stratified，無 Pooled OLS trap。R-04 補註 2026-04-19 的 within-group 不重跑決定維持有效。

## 整體評分
**⚠️ 措辭需精確化 + HP-independent 驗證空白 — 對應 C17 LOH Subclone 雙證據鏈是正向的鏡像案例；P0-B 已確認分析設計已為 stratified 無需重算。**
