# C12: ASM 32-66% POSITIVE（結論 12）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 12
- **來源文件**: `03_ISM分析價值界定.md` / `06_結論穩定性審查.md:266-273`
- **穩定度評分**: ⭐4 穩固
- **狀態**: POSITIVE
- **所屬證據鏈**: 🔗8（ASM Cross-validation）
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ 32-66% 範圍一致
- ✅ 5 方法 Jaccard 0.78-0.83 一致

## D2 方法論健康度
- ✅ 5 方法交叉驗證
- ✅ 7/7 樣本一致
- ⚠️ **未區分 germline vs somatic ASM**（C-BIO-1）→ **P2-C**
- ⚠️ 32% 可能含統計顯著但效應微弱位點（|delta|<0.2）
- ❌ 效應量門檻分析缺
- ❌ Multiple testing correction（跨 ASM 判定方法）未明示

## D3 證據鏈
- **依賴結論**: C01
- **被依賴結論**: C15, C17（LOH 區域 ASM 行為）
- **鏈完整度**: ✅ 完整但 germline/somatic 區分空白

## D4 數據信任度
- **dataset 版本**: Haplotag v5
- **CovM bug 影響**: 🟢 無
- **重跑必要性**: 需 P0-C 後輕度重測（HP 汙染影響 ASM 分類）

## D5 統計嚴謹度
- **Effect size**: 32-66% 範圍大 → 跨樣本異質
- **CI 覆蓋率**: ⚠️ 跨樣本範圍可作 CI 代理
- **Power 評估**: ✅ 跨樣本穩定
- **Multiple testing**: ⚠️ 跨 5 方法比較未明示 FDR

## D6 知識庫交叉驗證
- **支持文獻**: 文獻 germline ASM ~15-30%（Onuchic 2018）；somatic ASM 在癌症文獻多樣（Do 2017, epiTRACERx）
- **挑戰文獻**: 若 32% 多數為 germline ASM → 與 C-BIO-1 一致（未 overlap 文獻範圍）
- **缺口**: 個別樣本 germline vs somatic ASM 拆解

## 修正建議
- **P2-C**: 引用 germline ASM 文獻（15-30%）+ 註明 ASM 類型
- **P2**: 效應量門檻分析（|delta|>0.2 後比例變化）
- **Phase 2 A+D**: 透過 Normal BAM 區分 germline vs somatic ASM → **用戶新決策的 characterization 核心功能**
- **對應 R-id**: 無
- **對應 Q-id**: Q-08（germline vs somatic 區分依賴 Normal BAM）

## 整體評分
**✅ POSITIVE 穩固 — 文獻對照後需加 germline/somatic 拆解；此為「主軸 characterization」的關鍵功能。**
