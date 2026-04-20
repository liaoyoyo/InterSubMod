# C21: LOH.bed 不受 self-phasing 汙染（補充結論 21）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 21（補充，2026-04-19 R5 新增）
- **來源文件**: `06_結論穩定性審查.md:445-460` / `20260419_LOH_bed_generation_audit_01.md`
- **穩定度評分**: ⭐5 堅固
- **狀態**: CONFIRMED
- **所屬證據鏈**: 🔗7（Self-Phasing） / 🔗10
- **所屬樣本**: 7/7（C++ code 層面 sample-invariant）

## D1 內部一致性
- ✅ Jaccard=1.0000 / F1=96.2% vs SEQC2 一致
- ✅ C09 精確化對接（P2-A 修正後完全對齊）

## D2 方法論健康度
- ✅ Jaccard 1.0 極端穩健
- ✅ C++ code review 機制驗證
- ✅ `LohBedAnnotator::classify` coordinate-based（不讀 BAM HP tag）
- ⚠️ Paired mode `core_loh_like` 由 ClairS Paired 自帶，路徑正交（R5 不涵蓋）→ 範圍限制
- N/A bootstrap（Jaccard 1.0 不需）

## D3 證據鏈
- **依賴結論**: C09, C10, C13
- **被依賴結論**: C17（LOH.bed filter 合法性前提）、C22
- **鏈完整度**: ✅ C++ 機制審查完整

## D4 數據信任度
- **dataset 版本**: PON-only vs self-phasing 對比 HCC1395
- **CovM bug 影響**: 🟢 無（LOH.bed 純 VCF coordinate）
- **重跑必要性**: 無

## D5 統計嚴謹度
- **Effect size**: Jaccard 1.0000（極端）
- **CI 覆蓋率**: N/A（set operation）
- **Power 評估**: ✅ n=1094 regions
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: LongPhase `--loh` flag 官方 document（knowledge 02）
- **挑戰文獻**: 無
- **缺口**: 無

## 修正建議
- **P2-A**: 在 C09 精確化同步引用 C21 作為 evidence（「LOH.bed 證實不受 self-phasing 影響」）
- **對應 R-id**: 無
- **對應 Q-id**: 無

## 整體評分
**✅ ⭐5 結論穩固 — C09 精確化的關鍵佐證，對應 C13 P1-A 修正。**
