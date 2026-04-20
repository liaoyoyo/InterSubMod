# C20: Coverage_Multiple 非獨立 CN proxy（補充結論 20，原 "CN proxy" FALSIFIED）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 20（補充，2026-04-19 R6 新增）
- **來源文件**: `06_結論穩定性審查.md:426-441` / `20260419_Coverage_Multiple_zscore_normalization_01.md`
- **穩定度評分**: ⭐4 穩固
- **狀態**: FALSIFIED（「≈ CN proxy」原聲明）/ DOWNGRADED（僅 per-sample 排序用）
- **所屬證據鏈**: 🔗10
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ **`20260417_Zone_Aware_Confidence_Framework_01.md:49,73` r=0.997 已修正為 r=0.831（2026-04-20 P1-A）**；原誤標註記保留**（實為 CovM vs CN r=0.831）→ **P1-A 最高優先**
- ✅ R6 per-sample baseline 驗證數字一致

## D2 方法論健康度
- ✅ Per-sample z-score normalization 方法正確
- ❌ **發現 CovM baseline expected_coverage=75.0 hardcoded**（C-INFRA-1）→ **P0-A**（R-01 已決定立即修）
- ⚠️ C17 的 CovM vs CN r=0.831 可能本身受 bug 影響 → 修正後 r 值待測
- ❌ 無 bootstrap CI for r

## D3 證據鏈
- **依賴結論**: C13, C17
- **被依賴結論**: C17 step3 CN1 分層 / C19 Z3 zone 定義 / C22 Zone-Aware
- **鏈完整度**: ⚠️ **CovM bug 汙染下游整條 CN-based 分層**

## D4 數據信任度
- **dataset 版本**: Haplotag v5 cross_sample_audit
- **CovM bug 影響**: 🔴 **Bug 直接產生此結論**
- **重跑必要性**: **必需**（R-01 修正後重跑）

## D5 statistics rigor
- **Effect size**: z_extreme 0.15% << 5%（強否證）
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ✅
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: CNV caller 文獻（Delly, Manta, sequenza, CNVkit）—正確作法
- **挑戰文獻**: 無
- **缺口**: 整合 CNV caller 的具體 pipeline 設計

## 修正建議
- **P1-A**: 修正 20260417_Zone_Aware_Confidence_Framework_01.md:49,73 誤標
- **P0-A**: R-01 CovM KDE baseline 修正 + 全量重跑
- **Phase 2 A+D**: Normal BAM Δ_coverage 作為 sample-matched baseline（用戶主軸定位下的 characterization 功能）
- **P2**: CNV caller 整合設計（Phase 2 或 Part C.1）
- **對應 R-id**: R-01
- **對應 Q-id**: Q-01, Q-02

## 整體評分
**🔴 CovM bug 直接相關的結論 — R-01 立即修為最高優先。**
