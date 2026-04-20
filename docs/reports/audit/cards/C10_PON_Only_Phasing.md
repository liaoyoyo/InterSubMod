# C10: PON-only phasing VCF 正確但 haplotag 缺陷（結論 10）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 10
- **來源文件**: `02_Self_Phasing根因.md` / `04_暫停判定與重評估.md` / `06_結論穩定性審查.md:221-230`
- **穩定度評分**: ⭐4 穩固
- **狀態**: PARTIAL SUCCESS
- **所屬證據鏈**: 🔗7（Self-Phasing）
- **所屬樣本**: HCC1395-only（1/7）

## D1 內部一致性
- ✅ LOH.bed Jaccard=1.0 / N50 翻倍 / purity 0 一致
- ✅ haplotag GT 解析缺陷與 C21 LOH.bed audit 路徑相容

## D2 方法論健康度
- ✅ 嚴謹（bias 消除量化）
- ⚠️ **單樣本外推**（僅 HCC1395）
- ⚠️ haplotag GT 解析缺陷 C++ code review 未完成（P0-C 依賴）
- ❌ 無 bootstrap CI

## D3 證據鏈
- **依賴結論**: C02, C09
- **被依賴結論**: C11, C13, C21, Phase 1A F1（P0-C 解鎖條件）
- **鏈完整度**: ⚠️ 單樣本限制

## D4 數據信任度
- **dataset 版本**: PON-only phasing HCC1395 專屬 run
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: **需**多樣本（COLO829 / HCC1954）

## D5 統計嚴謹度
- **Effect size**: Jaccard 1.0 / N50 翻倍（極大）
- **CI 覆蓋率**: ❌ 單樣本無 CI
- **Power 評估**: ✅ n=1094 regions
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: LongPhase PON-only 使用模式在 ClairS 官方文件有（knowledge L1）
- **挑戰文獻**: 不同 PON 覆蓋率樣本的適用性
- **缺口**: 多樣本 PON-only 驗證

## 修正建議
- **P0-C**: haplotag ReadParser 修正
- **P1**: 多樣本 PON-only 驗證（COLO829, HCC1954）
- **P1**: Purity estimation 0→修復排查
- **對應 R-id**: 無
- **對應 Q-id**: Q-07（單樣本外推）

## 整體評分
**🟡 PARTIAL — 機制成立但需多樣本 + 下游 C++ 修正。**
