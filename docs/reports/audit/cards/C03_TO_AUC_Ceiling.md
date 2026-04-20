# C03: TO 無單一特徵 AUC > 0.58（結論 3）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 3
- **來源文件**: `docs/reports/research_landscape/03_ISM分析價值界定.md` / `05_證據鏈總覽.md`
- **穩定度評分**: ⭐4 穩固
- **狀態**: NEGATIVE
- **所屬證據鏈**: 🔗3（TO AUC ceiling）
- **所屬樣本**: 7/7 (748K regions)

## D1 內部一致性
- ✅ 60+ 特徵數字跨文件一致
- ⚠️ AUC 門檻 0.58：跨結論 cutoff 不一（C08 LOSO 0.721 稱 NO-GO；C11 +0.011 稱 POSITIVE）→ **P2-B 統一定義**

## D2 方法論健康度
- ✅ 大樣本掃描減低 pooled OLS 風險（單特徵 AUC）
- ✅ 無 L2 collider（單特徵非 residualized）
- ⚠️ **HP-dependent 特徵仍在 self-phasing 汙染狀態**（P0-C 未解決）
- ✅ 無 spatial auto-correlation（region-level）
- ❌ 無 bootstrap CI（single-point AUC claim）→ **P1-B 需補**
- ❌ 無 FDR（60+ 特徵掃描）→ **P1-B 需補**

## D3 證據鏈
- **依賴結論**: C01, C02
- **被依賴結論**: C07, C08, C14, C16（作為 ceiling reference）
- **鏈完整度**: ⚠️ 28% HP-dependent 特徵在汙染條件下測量（P0-C 前為「中間狀態」結論）

## D4 數據信任度
- **dataset 版本**: Haplotag v5 (canonical 7 × 3 modes)
- **CovM bug 影響**: 🟡 部分特徵含 Coverage_Multiple → 重算後 AUC 可能微變（<0.02）
- **重跑必要性**: **需**（P0-C 完成後）

## D5 統計嚴謹度
- **Effect size**: AUC ≤ 0.58 = 弱分離
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ✅ 大樣本
- **Multiple testing**: ❌ 60+ 特徵無 BH-FDR 校正

## D6 知識庫交叉驗證
- **支持文獻**: 知識庫無直接 counterpart（ISM 獨特 feature space）
- **挑戰文獻**: 無
- **缺口**: 無

## 修正建議
- **P0-C**: haplotag ReadParser 修正後 29 HP-dependent 特徵重測
- **P1-B**: 對 60+ 特徵套 BH-FDR + bootstrap CI
- **P2-B**: 統一 AUC 0.58 門檻定義
- **對應 R-id**: R-04（NO-GO 結論重審）、R-05（pooled OLS 全面回溯）
- **對應 Q-id**: Q-06

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: AUC<0.58 已於多輪分析跨 7 樣本確認（60+ 特徵、全維度失敗）；n_reads confound 並非主根因（TO AUC ceiling 根因為 methylation feature space 耗盡）。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**⚠️ 結論方向穩固但統計嚴謹度不足 — 需補 bootstrap + FDR + haplotag 修正後重測。**
