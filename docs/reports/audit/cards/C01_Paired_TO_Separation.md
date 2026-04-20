# C01: Paired/TO 必須分離（結論 1）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 1
- **來源文件**: `docs/reports/research_landscape/01_TO_FP問題全貌.md` / `05_證據鏈總覽.md` 🔗1
- **穩定度評分**: ⭐5 堅固
- **狀態**: CONFIRMED
- **所屬證據鏈**: 🔗1（FP Provenance）、🔗2（Paired/TO structural difference）
- **所屬樣本**: 7/7（結構性差異，不隨樣本變化）

## D1 內部一致性
- ✅ 跨文件數值一致（FP rate 1.04% vs 30.6% 所有引用一致）
- ✅ 術語定義一致（Paired = 含 normal, TO = tumor-only）
- ⚠️ TP/FP 絕對數字：最新矯正版 TO=28,509/11,606（2026-04-04 矯正）；部分 archive 仍用舊 30,490/4,842 → **P3-D 同步**

## D2 方法論健康度
- ✅ 無 pooled OLS trap（結構性差異非 regression）
- ✅ 無 L2 collider bias
- ✅ 無 n_reads confound（FP rate 比例比較不受 coverage 影響）
- ✅ 無 spatial auto-correlation
- ✅ 有 bootstrap 支撐（HP r=0.001, n=288K pairs）
- ⚠️ 無 FDR 校正（但結論不基於 AUC 掃描，4 條獨立證據線）

## D3 證據鏈
- **依賴結論**: 無（最基礎）
- **被依賴結論**: C02, C03, C07, C09, C14, C16, C17（所有 TO/Paired mode 區分）
- **鏈完整度**: ✅ 完整（4 條證據線：FP rate / HP r / 甲基化方向 / QS AUC）

## D4 數據信任度
- **dataset 版本**: Haplotag v5 (canonical)
- **CovM bug 影響**: 🟢 無影響（不依賴 CovM）
- **重跑必要性**: 無

## D5 統計嚴謹度
- **Effect size**: FP rate 差 29.6pp（極大）；HP Spearman r=0.001（無相關）
- **CI 覆蓋率**: ✅ 大樣本 n=288K → 窄 CI
- **Power 評估**: ✅ 極大效應量，power→1
- **Multiple testing**: N/A（非掃描結論）

## D6 知識庫交叉驗證
- **支持文獻**: ClairS/ClairS-TO 原論文（tumor-only 設定 FP 結構性升高）
- **挑戰文獻**: 無
- **缺口**: Paired vs TO 的 HP consistency 在人類其他 cell line panel 是否普遍（knowledge base: ONT basecalling / SEQC2 區塊可查）

## 修正建議
- **P3-D**: 早期 archive 的舊數據 30,490/4,842 → grep 全檔替換為 28,509/11,606
- **對應 R-id**: 無
- **對應 Q-id**: 無

## 整體評分
**✅ Clean conclusion — 唯一需修為文件同步。**
