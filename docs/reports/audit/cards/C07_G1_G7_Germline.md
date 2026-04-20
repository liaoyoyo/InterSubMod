# C07: G1-G7 TO Germline FP 不可識別（結論 7）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 7
- **來源文件**: `03_ISM分析價值界定.md` / `05_證據鏈總覽.md` / `06_結論穩定性審查.md:164-173`
- **穩定度評分**: ⭐4 穩固
- **狀態**: NO-GO
- **所屬證據鏈**: 🔗5（Germline FP identification）
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ 40+ VCF 特徵 × 7 模組 × 7 樣本一致
- ✅ LOSO AUC=0.638 / bootstrap=0.503 一致
- ✅ 48 圖表引用一致

## D2 方法論健康度
- ✅ 無 pooled OLS（LOSO 是跨樣本 hold-out）
- ✅ 無 L2 collider
- ⚠️ n_reads confound 未全面檢查 40+ VCF 特徵 → **R-04 統一檢查建議**
- ⚠️ LR 為主，未測 GBM/RF（但 LOSO 0.638 已接近 LR 天花板）
- ⚠️ FP removal=0% 的「安全約束 TP loss ≤2%」定義是否過嚴？→ C08 類似議題
- ❌ 40+ 特徵 AUC 掃描無 FDR → **P1-B**

## D3 證據鏈
- **依賴結論**: C01, C03
- **被依賴結論**: C08（read-level FP 的 VCF 對照）、C14（QS 失效根因）
- **鏈完整度**: ✅ 完整

## D4 數據信任度
- **dataset 版本**: Haplotag v5（VCF 特徵不受 HP 影響）
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: 僅 R-04 統一檢查時

## D5 statistics rigor
- **Effect size**: LOSO AUC 0.638（弱）/ bootstrap 0.503（null）
- **CI 覆蓋率**: ⚠️ bootstrap=0.503 可能暗含 CI 但未明報
- **Power 評估**: ✅ LOSO 設計穩健
- **Multiple testing**: ❌ 40+ 特徵無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: Germline vs somatic classification 文獻廣泛，但在 FP rescue 場景罕見專論（knowledge 05）
- **挑戰文獻**: GBM 可能非線性突破（但 LR 已 0.638 空間有限）
- **缺口**: Mutational signature features 未納入

## 修正建議
- **P1-B**: 40+ VCF 特徵套 BH-FDR
- **P2**: 補 GBM/RF 測試（P2 優先序）
- **R-04**: 統一 n_reads confound 跨 13 NO-GO 檢查
- **對應 R-id**: R-04, R-05
- **對應 Q-id**: Q-06

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: 60+ VCF 特徵全 AUC<0.64 已於 48 圖表跨 7 樣本驗證；TP loss ≤2% 下 FP removal=0%；n_reads 非主驅動因素（VCF-level 特徵）。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**✅ NO-GO 穩固 — 唯一風險是 GBM/RF 未測試的理論可能性。**
