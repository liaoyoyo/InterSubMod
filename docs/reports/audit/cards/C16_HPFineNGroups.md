# C16: HPFineNGroups somatic heterogeneity marker（補充結論 16）

> **建立日期**: 2026-04-19

## 基本資訊
- **結論編號**: 16（補充）
- **來源文件**: `09_Part_B.md` / `06_結論穩定性審查.md:312-342`
- **穩定度評分**: ⭐4 穩固 ⬆️（2026-04-18 F pilot 升級）
- **狀態**: POSITIVE REFINED
- **所屬證據鏈**: 🔗9（Part B 新證據鏈）
- **所屬樣本**: 7/7（5 ≥0.85, HCC1954 +21pp 挽救）

## D1 內部一致性
- ❌ **新 canonical filter (NG=4+AF<0.4+NR≥80) vs 舊 (NG≥4+NR≥80)**：
  - F pilot 2026-04-18 新 filter
  - `docs/CURRENT_FOCUS.md:98-103` R4 仍用舊 → **P3-A**
  - `docs/experiments/INDEX.md` 部分段落仍引用舊
- ✅ TP rate 0.9281（新）/ 0.8912（舊）數字一致

## D2 方法論健康度
- ✅ chr-shuffle null Z=43.5（非 spatial auto-correlation）
- ✅ Coverage_Multiple 跨 CN tiers 0.90-0.94 穩定
- ❌ **Pooled OLS residualization trap 疑慮**（C-STAT-1）：residualized AUC=0.617 是 pooled 還是 within-group？→ **P0-B**
- ❌ **HPFineNGroups residualized AUC 缺 bootstrap CI**（C-STAT-2）→ **P1-C**
- ⚠️ **NG=3 非單調**：NG=2 TP rate 0.643 < NG=1 0.763，根因 germline AF confound → **P2-D**
- ⚠️ HCC1954 AF<0.4 挽救機制需 Normal BAM 驗證
- ❌ 無 FDR（跨 NG bin × 樣本掃描）

## D3 證據鏈
- **依賴結論**: C01, C09, C12, C15
- **被依賴結論**: C17（LOH Subclone 機制共享）、C22（Zone Characterization）
- **鏈完整度**: ✅ 🔗9 新鏈完整

## D4 數據信任度
- **dataset 版本**: Haplotag v5 + F pilot 新 filter
- **CovM bug 影響**: 🟡 **B.2-2 Coverage_Multiple 跨 CN tiers 0.90-0.94 使用 expected_coverage=75.0 hardcoded baseline** → 修正後需重算
- **重跑必要性**: **需**（CovM bug 修正後）

## D5 統計嚴謹度
- **Effect size**: Δ +3.7pp (0.8912→0.9281)；HCC1954 +21pp；per-sample ΔAUC 1/7 ROBUST
- **CI 覆蓋率**: ❌ 缺 → P1-C
- **Power 評估**: ⚠️ COLO829 n=34 小樣本（out-of-scope 5mC basecall）
- **Multiple testing**: ❌ 多 NG × NR × AF bin 無 FDR

## D6 知識庫交叉驗證
- **支持文獻**: HPFineNGroups 概念與 single-cell heterogeneity 文獻相容（knowledge 03）
- **挑戰文獻**: epiTRACERx 的 methylation subclone marker 與此類似嗎？需查
- **缺口**: 非 SEQC2 patient-derived cohort 外推

## 修正建議
- **P0-B**: residualized AUC within-group OLS 重算
- **P1-C**: 補 1000× bootstrap CI
- **P2-D**: NG=3 非單調機制 / germline AF confound 補充討論
- **P3-A**: 同步舊→新 filter 版本
- **Phase 2 A+D**: Normal BAM 驗證 HCC1954 AF<0.4 機制
- **對應 R-id**: R-02（主軸 characterization 下的功能）、R-05
- **對應 Q-id**: Q-09（NG=3 機制）

## 整體評分
**✅ POSITIVE REFINED — 主軸 characterization 下的核心 Subclone Marker 功能；急需 P0-B within-group OLS 驗證 + P1-C bootstrap CI。**
