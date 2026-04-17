# F Pilot Step 3 — 新 Filter (NG=4+AF<0.4) 驗證

**Date**: 2026-04-18
**Script**: `scripts/step3_af_cutoff_validation.py`

---

## 🎯 核心結論：新 filter `NG=4 + AF<0.4 + NR≥80 NonLOH` 顯著優於舊 filter

| Filter | n | TP rate | CI | Coverage loss |
|---|---|---|---|---|
| Overall NonLOH | 307,474 | 0.6699 | — | — |
| Old (NG≥4+NR≥80) | 25,744 | 0.8912 | [0.887, 0.895] | — |
| **New (NG=4+AF<0.4+NR≥80)** | **14,197** | **0.9281** | [0.924, 0.932] | 45% |

**Overall gain: +3.7pp TP rate, 55% coverage retained**

---

## Per-sample 改善量（最關鍵發現）

| sample | old TP rate | new TP rate | Δ | coverage loss |
|---|---|---|---|---|
| **HCC1954** | 0.497 | **0.707** | **+21.0pp** ⭐ | 36% |
| **HCC1937** | 0.714 | **0.867** | **+15.4pp** ⭐ | 30% |
| HCC1395 | 0.810 | 0.887 | +7.7pp | 27% |
| H1437 | 0.921 | 0.965 | +4.3pp | 20% |
| H2009 | 0.935 | 0.957 | +2.3pp | 50% |
| HCC1395_DORADO | 0.903 | 0.919 | +1.6pp | 11% |
| COLO829 | 0.235 | 0.235 | 0pp | 0% (n=34 已 all AF<0.4) |

**5/7 樣本達 TP rate ≥ 0.85**（舊 filter 僅 4/7）

HCC1954 與 HCC1937 的大幅改善完全符合 Step 2 根因分析：其 FP 集中在 AF≥0.4，移除後 TP rate 大躍進。

---

## AF cutoff sensitivity 分析

| AF cut | 5/7 達 TP≥0.85 | 6/7 達 TP≥0.85 | n total | Sum TP rate mean |
|---|---|---|---|---|
| 0.20 | ✅ (H2009=0.985, H1437=0.958, HCC1395_DORADO=0.895, HCC1395=0.891, HCC1937=0.949, HCC1954=0.874) | **6/7** (只差 COLO829) | 3,027 | 0.923 (excl COLO) |
| 0.30 | ✅ | 6/7 | 8,283 | 0.907 |
| 0.35 | ✅ | 6/7 | 11,674 | 0.899 |
| **0.40** | ✅ | 5/7 (HCC1954 0.707) | **14,197** | 0.881 |
| 0.45 | ✅ | 5/7 | 16,217 | 0.869 |
| 0.50 | ✅ | 5/7 | 18,087 | 0.856 |

**AF<0.2 是高 precision 極端**（6/7 ≥0.87，coverage 極少）
**AF<0.4 是平衡點**（5/7 ≥0.85，coverage 14K regions，可解釋性強）

建議：
- **論文 figure / 論證用**：AF<0.2（6/7 high-confidence，但須註明 coverage 限制）
- **通用 filter / 大規模分析**：AF<0.4（平衡）
- **COLO829 仍無效**（n<50，ONT_R10 無 methylation 背景）

---

## Confound 檢查（全 PASS）

### C1: Coverage_Multiple 分層（CN confound）

| CN bin | n | TP rate |
|---|---|---|
| 0.0-1.1 (CN≤1) | 472 | 0.900 |
| 1.1-1.5 | 5,494 | 0.923 |
| 1.5-2.0 (CN~2) | 5,183 | 0.944 |
| 2.0-3.0 (CN~3) | 1,885 | 0.913 |
| 3.0+ (CN≥4) | 1,163 | 0.915 |

所有 CN bins TP rate 0.90-0.94 → **不是 CN confound**（新 filter 跨 CN tiers 穩定）

### C2: Chr-shuffle null (spatial confound check)

- Observed TP rate: 0.9281
- Null mean (within-chr label shuffle, 20 iter): 0.6898 ± 0.0055
- **Z-score: 43.5**（完全壓倒 null）
- → PASS: 新 filter 的 gain 不是 spatial auto-correlation artifact

避免了 P3 pilot 的教訓。

---

## COLO829 問題診斷

- COLO829 NG≥4+NR≥80 NonLOH 只有 n=34，**全部 AF<0.4**
- 所以 new filter 對 COLO829 無效應（無可刪除 FP）
- COLO829 TP rate 0.235 是基本面：ONT_R10 無 methylation → HPFineNGroups 本身 artifact
- memory 已記：**COLO829 ONT_R10 無 methylation、ONT_PAO 才有 5mCG+5hmCG**
- **結論：COLO829 不在 HPFineNGroups filter 的適用範圍**，應明確排除

---

## 對其他任務的影響

### 需更新的結論

1. **Memory `project_hpfinengroups_subclone_marker.md`**：
   - **Old**: "N≥4+NR≥80 NonLOH TP rate 89.1%"
   - **New**: "NG=4+AF<0.4+NR≥80 NonLOH TP rate 92.8% (5/7 samples >=0.85)"
   - 需加: "AF<0.2 for high-precision 6/7 ≥0.87"
   - 需加: "COLO829 不適用"

2. **結論穩定度補充結論 16**：
   - 可升級 ⭐3 → ⭐4（add Step 3 驗證證據：chr shuffle Z=43.5、CN 穩定、5/7 跨樣本）

3. **B.1-3 Cohen's d per-sample 效應**：
   - 若 Cohen's d 重算 with AF<0.4，HCC1954/COLO829 可能從「特殊」→「POS」
   - 建議後續以 AF<0.4 stratified Cohen's d 作為新 canonical effect size

4. **Phase 2 特徵化主線**：
   - HPFineNGroups 作為 **biology-level subclone marker** 得到強化
   - 具體 operational definition: NG=4 + AF<0.4 (＋NR≥80, NonLOH)
   - HCC1954/HCC1937 的「樣本特異性失效」解釋為 AF confound，非訊號本身不存在

### 不受影響（確認）

- LOH Subclone AF×Methylation (B.2)：正交，AF bin 框架不變
- Zone-Aware Framework：無 NGroups 依賴
- Per-CpG ASM characterization：獨立研究

### 未解決 / 新產生問題

1. **AF<0.4 為什麼有效？** 生物學解釋：
   - subclonal somatic SNV 天然 AF 較低（純度<1 或 subclonal events）
   - germline het 在 AF≈0.5；移除中 AF 剔除 germline-like FP
   - **注意**: 若 LOH 讓 somatic AF 漂到 >0.5, AF<0.4 會誤除真 somatic (但 NonLOH filter 已排除此情境)

2. **COLO829 應否有獨立 filter？** 考慮：
   - ONT_R10 無 methylation → 任何 HPFineNGroups-based filter 失效
   - 應增 sample-level gating：優先 ONT_5mCG 或 5mCG+5hmCG basecall 樣本
