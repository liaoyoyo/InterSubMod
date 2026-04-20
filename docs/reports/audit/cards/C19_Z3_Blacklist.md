# C19: Z3 amplicon blacklist CONDITIONAL-NEGATIVE-for-canonical（補充結論 19）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 19（補充，2026-04-19 新增）
- **來源文件**: `06_結論穩定性審查.md:407-422` / `20260419_Z3_amplicon_blacklist_pilot_result_01.md`
- **穩定度評分**: ⭐5 堅固
- **狀態**: CONDITIONAL-NEGATIVE-for-canonical / characterization-only
- **所屬證據鏈**: 🔗10
- **所屬樣本**: 7/7

## D1 內部一致性
- ✅ HCC1954-only +0.0065 / 其他 5/6 mean=-0.0044 一致
- ✅ S2 ceiling 87% / S4 reject-all-Z3 +0.0075 一致

## D2 方法論健康度
- ✅ 4 策略 × 7 樣本完整矩陣
- ✅ Ensemble average + per-sample 並呈
- ⚠️ TO mode NonLOH 外推至 Paired mode 未驗
- ✅ 其他 HER2+ 樣本 blacklist 效應未驗（限制已標）
- N/A bootstrap/FDR（簡單 F1 對比）

## D3 證據鏈
- **依賴結論**: C01, C18, C22
- **被依賴結論**: 無（已關閉 canonical 方向）
- **鏈完整度**: ✅ 完整

## D4 數據信任度
- **dataset 版本**: Z3 pilot 獨立 run
- **CovM bug 影響**: 🟡 Z3 zone 定義（CovM 分層）受影響 → 若 Z3 邊界微變，blacklist 範圍微變但結論不變
- **重跑必要性**: CovM bug 修後 zone 重定義 + Z3 blacklist 重跑

## D5 統計嚴謹度
- **Effect size**: ΔF1 ±0.005（小）
- **CI 覆蓋率**: ❌ 缺
- **Power 評估**: ⚠️ 小效應量但跨樣本一致否定 → 結論穩
- **Multiple testing**: N/A

## D6 知識庫交叉驗證
- **支持文獻**: Region blacklist 在 variant filter 常見策略
- **挑戰文獻**: 無
- **缺口**: 無

## 修正建議
- **P1**: R-01 CovM bug 修後 zone 重定義 + Z3 blacklist 重跑確認（預期結論不變）
- **對應 R-id**: R-01, R-04
- **對應 Q-id**: Q-02

## R-04 補註（2026-04-19）
**Within-group OLS not re-run**: ⭐5 結論（Z3 amplicon 已知 artifact）穩固；CovM bug 修正後僅影響 zone 邊界微調。採用 CHECKLIST R-04 選項 A：接受現有 NEGATIVE 結論；within-group OLS 不再重跑。

## 整體評分
**✅ ⭐5 結論穩固 — CovM bug 修正後預期結論穩定（只微調 zone 邊界）。**
