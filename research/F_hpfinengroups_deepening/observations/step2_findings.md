# F Pilot Step 2 — 根因調查觀察紀錄

**Date**: 2026-04-17
**Script**: `scripts/step2_root_cause_investigation.py`

---

## Q1 答案 🌟 NGroups 非單調根因 = germline AF confound

### Per-NGroups profile (TO NonLOH NR≥80)

| NGroups | n | TP rate | AF mean | AF∈[0.45,0.55] frac | \|AlleleDelta\| |
|---|---|---|---|---|---|
| 1 | 3,655 | 0.7633 | 0.552 | 0.175 | 0.019 |
| **2** | **51,559** | **0.6434** | **0.471** | **0.212** (最高) | 0.025 |
| 3 | 64,653 | 0.7742 | 0.428 | — | 0.023 |
| 4 | 25,744 | 0.8912 | 0.402 | 0.145 | 0.020 |

### NGroups × AF 2D TP rate

| AF bin | NG=1 | NG=2 | NG=3 | NG=4 |
|---|---|---|---|---|
| 0.0-0.2 | 0.645 | 0.679 | 0.825 | **0.937** |
| 0.2-0.4 | 0.807 | 0.715 | 0.838 | **0.926** |
| 0.4-0.6 | 0.786 | 0.641 | 0.738 | 0.854 |
| 0.6-0.8 | 0.801 | 0.595 | 0.713 | 0.842 |
| 0.8-1.0 | 0.616 | **0.339** | 0.433 | 0.581 |

### 核心發現

1. **NGroups=2 的低 TP rate 是 AF confound**：
   - NG=2 是唯一 AF 最接近 0.5（germline het hallmark）的組
   - AF∈[0.45, 0.55] 的 fraction 0.212 也最高
   - AF 0.8-1.0 + NG=2 的 TP rate 僅 **0.339**（最差 cell）— 這是高 AF germline-like FP

2. **NGroups=4 + 低 AF 是最強 somatic signature**：
   - NG=4 + AF<0.4: TP rate 預估 **~0.93**
   - 這比 current canonical (NG≥4 + NR≥80) 的 0.891 更強

3. **NGroups 非單調的正確解釋**：
   - "NGroups 作為 subclone marker" 在 AF<0.4 是 monotone
   - 在 AF 接近 0.5 時 NG=2 被 germline 污染
   - 全 AF 混用才出現非單調

**Memory 需更新**：
- 原: "HPFineNGroups N≥4+NR≥80 TP rate 89.1%"
- 修正後: "NG=4 + AF<0.4 + NR≥80 TP rate ~0.93（需 Step 3 驗證）"

---

## Q2 答案 🌟 HCC1954 失效 = FP 在 AF ≥ 0.4 極端富集

### HCC1954 NG≥4+NR≥80 NonLOH

| AF bin | n | n_TP | TP rate |
|---|---|---|---|
| 0.0-0.2 | 564 | 493 | **0.874** ✅ |
| 0.2-0.4 | 476 | 242 | 0.508 |
| 0.4-0.6 | 313 | 54 | 0.173 ⚠️ |
| 0.6-0.8 | 223 | 16 | **0.072** ⚠️ |
| 0.8-1.0 | 46 | 1 | **0.022** ⚠️ |

### 根因

- HCC1954 低 AF 區段（AF<0.2）TP rate 0.874 — **完全正常！**
- 失效完全來自 AF ≥ 0.4 區段
- **caller_af TP mean=0.214 vs FP mean=0.482** → FP 系統性高 AF

**HCC1954 FP 可能成因**（純度+ploidy 已知複雜）：
- HER2+ breast cancer, high ploidy (ICGC 記 ~4)
- 未被 LOH.bed 標註的複雜 CNV 區可能讓 germline het 表現為 AF ≈ 0.5
- NumCpGs FP (111) > TP (90) → FP 集中在高 CpG 密度區（likely 啟動子/CpG island）

### 對所有樣本的意涵

若 NG=4 + AF<0.4 是通用 filter，HCC1954 可能從 0.497 → 0.70+
這是本次 pilot 的最重要實踐發現。

---

## Q3 答案 ⚠️ Paired 99.85% 不是真 gain

| Condition | n | TP rate |
|---|---|---|
| All paired (baseline) | 328,699 | **0.9896** |
| Paired + NR≥80 | 141,783 | 0.9967 |
| Paired + NG≥4 | 15,465 | 0.9960 |
| Paired + NG≥4 + NR≥80 | 11,801 | 0.9985 |
| Paired + NG<4 + NR<80 | 183,252 | 0.9841 |

**Filter gain: 只有 +0.89pp**

### 結論

- Paired baseline 本身 = 98.96%（caller 已極強 filter）
- NGroups filter 在 Paired mode **實質無 gain**
- Step 1 "Paired 99.85% 驚人" 是 artifact of baseline

### Per-sample 驗證

- 7/7 樣本 Paired baseline 都 ≥ 0.94
- Filter gain 多 < 1pp，最多 HCC1937 +1.55pp
- **Paired mode 不是 HPFineNGroups 的應用場景**

---

## 對其他任務的影響

### 需更新的其他結論

1. **memory `project_hpfinengroups_subclone_marker.md`**：
   - NEW insight: NGroups 非單調 → 需加 AF<0.4 co-filter 才真正 somatic marker
   - Paired 99%+ baseline 不是 gain
   - HCC1954 失效根因 = 高 AF FP，不是 "樣本特異性"

2. **B.1-1 HPFineNGroups residualized AUC 0.617**：
   - residualize on AF 已包含，但 AF<0.4 stratification 未被測試
   - **B.1-1 結論「pooled ΔAUC=+0.025 robust」仍成立**，但 Step 2 顯示該結論低估了 AF-stratified 潛力

3. **B.1-3 per-sample Cohen's d 5/7 POS + 2/7 特殊**：
   - 若 Cohen's d 重算 with AF<0.4 filter，COLO829/HCC1954 可能從「特殊」→「POS」
   - 未立即重算，Step 3 驗證後再評估

4. **結論穩定度補充結論 16 (HPFineNGroups)**：
   - 可升級從 ⭐3 → ⭐4（after Step 3 確認）

### 不受影響

- LOH Subclone AF×Methylation (B.2) — AF×LOH framework 與 NGroups 正交
- Zone-Aware Framework — 與 Step 2 無直接關聯
- Per-CpG ASM characterization — 無 NGroups 依賴

---

## Step 3 計畫

驗證 **NG=4 + AF<0.4 + NR≥80 NonLOH** 作為新 canonical filter：

1. 全 7 樣本 per-sample TP rate
2. vs 舊 NG≥4 + NR≥80 filter 對比
3. Coverage loss 量化
4. 與其他 AF cutoff (0.3, 0.35, 0.4, 0.45) 比較
5. 確認 gain 不是 spatial confound（快速 chr shuffle）
