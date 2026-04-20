# C09: Self-phasing 是 TO HP_Ratio LOH 主因（結論 9）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 9
- **來源文件**: `02_Self_Phasing根因.md` / `06_結論穩定性審查.md:209-217`
- **穩定度評分**: ⭐4 穩固
- **狀態**: CONFIRMED（精確化後）
- **所屬證據鏈**: 🔗7（Self-Phasing Causal Chain）
- **所屬樣本**: 7/7 (288K pairs)

## D1 內部一致性
- ❌ **「62% TO LOH 消失」措辭層次混用**：LOH.bed（region-level VCF）vs ISM HP_Ratio LOH（site-level BAM HP tag）
  - `docs/CURRENT_FOCUS.md:53` 未補註層次
  - `07_LOH_CN_AF_研究總整理.md:102` 未補註層次
  - `06_結論穩定性審查.md:215` **已精確化**
  - → **P2-A 跨文件統一**
- ✅ 17.3:1 bias / r≈0 / 86.5% 平衡 跨文件一致
- ✅ Cohen's d=-1.20 一致

## D2 方法論健康度
- ✅ 多樣本 (7/7) × 多指標（HP ratio, r, balance）
- ⚠️ **演算法層面根因缺**：為何 LongPhase-TO 94.6% somatic → HP1？無 source code 分析 → **P2-E**
- ✅ 無 pooled OLS（相關係數跨樣本 × 位點，非 residualized）
- ⚠️ Paired phasing 作為 ground truth 有自身限制（switch error），31.2% 比例可能需修 20-25%

## D3 證據鏈
- **依賴結論**: C01
- **被依賴結論**: C03 HP-dependent, C05 LOH, C08 within_dom, C10 PON-only, C15 LOH 甲基化（所有 HP 相關）
- **鏈完整度**: ✅ 完整，但措辭混用影響下游理解

## D4 數據信任度
- **dataset 版本**: Haplotag v5（本身是被質疑的對象）
- **CovM bug 影響**: 🟢 無（self-phasing 驗證不用 CovM）
- **重跑必要性**: 無（但需配合 P0-C 後多樣本 PON-only 驗證）

## D5 統計嚴謹度
- **Effect size**: d=-1.20 超大；bias 17.3:1
- **CI 覆蓋率**: ⚠️ r≈0 未明報 CI
- **Power 評估**: ✅ n=288K pairs
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: LongPhase / WhatsHap 等 somatic phasing 文獻（knowledge 02）
- **挑戰文獻**: 若 LongPhase-TO 官方文獻指出此 bias 是預期行為，則結論為已知（需查）
- **缺口**: LongPhase-TO `--loh` flag 演算法層面文獻

## 修正建議
- **P2-A**: CURRENT_FOCUS.md:53 / 07_LOH_CN_AF.md:102 補註「62% 特指 ISM HP_Ratio LOH，不含 LOH.bed」（對應 C21）
- **P2-E**: 補 LongPhase-TO 機制層分析（source code 審查或 `--loh` flag 文件）
- **對應 R-id**: 無（結論穩固）
- **對應 Q-id**: 無

## 整體評分
**⚠️ 結論正確但措辭混用 — P2-A/P2-E 為關鍵文件精確化。**
