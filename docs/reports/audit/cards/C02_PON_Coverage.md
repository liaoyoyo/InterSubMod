# C02: PON 移除 99.48% raw FP（結論 2）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 2
- **來源文件**: `docs/reports/research_landscape/01_TO_FP問題全貌.md`
- **穩定度評分**: ⭐3 需注意
- **狀態**: CONFIRMED（單樣本）
- **所屬證據鏈**: 🔗1（FP Provenance）
- **所屬樣本**: HCC1395-only（1/7）

## D1 內部一致性
- ✅ 數值 2,717,339/2,731,541 在各文件一致
- ✅ 百分比 99.48% 一致

## D2 方法論健康度
- ✅ 無 pooled OLS / L2 collider / n_reads / spatial confound
- ✅ 單純集合運算，不需 bootstrap CI
- N/A FDR（非 AUC 掃描）

## D3 證據鏈
- **依賴結論**: C01（Paired/TO 分離前提）
- **被依賴結論**: C03（PON 過濾後剩 14,202 FP = germline leak）
- **鏈完整度**: ✅ 單樣本下完整

## D4 數據信任度
- **dataset 版本**: Haplotag v5 HCC1395
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: **需**跨樣本驗證（非 bug 驅動）

## D5 統計嚴謹度
- **Effect size**: 99.48%（壓倒性）
- **CI 覆蓋率**: ❌ 單樣本無法計算跨樣本 CI
- **Power 評估**: 絕對數字夠大（n=2.73M），但無法外推
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: PON-based germline filter 在 ClairS-TO 官方設計即為核心（knowledge base L1-L4）
- **挑戰文獻**: 不同族群的 rare germline 組成不同 → PON 覆蓋率可能降低
- **缺口**: HCC1954, COLO829 等其他 6 樣本 PON 移除率未驗證

## 修正建議
- **P2**: 在 6 個非 HCC1395 樣本計算 PON 移除率（即使跌到 95%，結論「PON 是最強 filter」不變）
- **對應 R-id**: 無（結論穩固，僅缺外推驗證）
- **對應 Q-id**: Q-07（單樣本外推）

## 整體評分
**✅ 結論穩固但外推不足 — 建議納入 P2 補驗證。**
