# C13: LOH.bed Jaccard=0.847 vs SEQC2（結論 13）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 13
- **來源文件**: `06_結論穩定性審查.md:276-283`
- **穩定度評分**: ⭐3 需注意
- **狀態**: 觀察值
- **所屬證據鏈**: 🔗1 / 🔗7
- **所屬樣本**: HCC1395-only（1/7）

## D1 內部一致性
- ✅ Jaccard=0.847 一致
- ✅ **`20260417_Zone_Aware_Confidence_Framework_01.md:49,73` r=0.997 已修正為 r=0.831（2026-04-20 P1-A）**；原誤標註記保留（實為 LOH.bed concordance）** → **P1-A**

## D2 方法論健康度
- ✅ 獨立 truth set 對比
- ⚠️ **單樣本**（HCC1395）
- ⚠️ SEQC2 主要為 SNV benchmark，LOH zone 準確性與 SNV 不同
- ✅ Jaccard 非 AUC → 不適用 bootstrap
- N/A Multiple testing

## D3 證據鏈
- **依賴結論**: C01
- **被依賴結論**: C21（LOH.bed 不受 self-phasing 汙染）、C17（LOH Subclone 前提）
- **鏈完整度**: ⚠️ 單樣本

## D4 數據信任度
- **dataset 版本**: HCC1395 LOH.bed + SEQC2 truth
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: **需**多樣本

## D5 統計嚴謹度
- **Effect size**: Jaccard 0.847（高）
- **CI 覆蓋率**: ❌ 單樣本
- **Power 評估**: N/A
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: SEQC2 truth set 在 LOH 方面的限制（knowledge L1-L4）
- **挑戰文獻**: 15.3% 不重疊可能含重要生物學差異
- **缺口**: 多樣本 LOH 獨立 truth（CNV caller Delly/Manta/sequenza）

## 修正建議
- **P1-A**: 修正 `20260417_Zone_Aware_Confidence_Framework_01.md:49,73` 誤標 r=0.997（應為 Jaccard 0.847 或 C21 Jaccard 1.0）
- **P2**: 多樣本 LOH.bed 驗證（Phase 2 A+D 整合 Normal BAM）
- **對應 R-id**: 無
- **對應 Q-id**: Q-07

## 整體評分
**⚠️ 單樣本觀察值 — 主要需修正 P1-A 數值錯標。**
