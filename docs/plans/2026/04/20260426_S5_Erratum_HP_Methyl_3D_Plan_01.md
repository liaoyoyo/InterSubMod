---
id: plan_20260426_s5_erratum_hp_methyl_3d
date: 2026-04-26
author: AI session (Opus 4.7)
priority: P0 (Plan A, C) / P1 (Plan B, D)
status: planning
estimated_effort: A=3h, B=1h, C=30min, D=4-6h, total=8-11h
related_artifacts:
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_s5b_hp_ratio_color/
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_v2_1_corrected/
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_obs18_af_split/
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_hp_ratio_3d_auc/
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/
related_pptx:
  - InterSubMod/docs/presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/
related_memory:
  - project_loh_constrained_phasing_discovery.md
  - project_loh_subclone_af_methylation_positive.md
  - project_p4_target3_depends_on_target1.md
  - project_zone_aware_framework.md
---

# Plan · S5 Erratum + HP_Ratio × LOH × AF × CN × Methylation 三層整合

## 0. Executive Summary

2026-04-25 至 04-26 的連續分析揭露了三個結構性發現，需要透過 Plan A-D 落實到報告、MEMORY、PPT、Filter prototype 四個層面：

1. **「Inner × Extreme AF = TP-pure」是錯誤觀察**（對應 S5 簡報）。真正的 high-TP 區是 **Inner × Mid AF (0.1–0.9) × Loss/Sub CN**，AF 兩端反而 FP-rich 或反向。
2. **HP_Ratio 在 Inner × Mid AF × Loss/Sub/Diploid/Gain 是有效二級切分器**（AUC 0.65–0.77，6 樣本 4–5/6 一致）；甲基訊號加成 +1~+6pp，combo 上限 AUC ≈0.79。
3. **LongPhase-TO baseline (舊 HP) 不可信**：v5 揭露 38% LOH 翻轉、60% NGroups 改變、Inner HP_Ratio 從單峰 1.0 變雙峰 0/1（揭穿 self-phasing artifact）；但甲基特徵 AUC 完全不變，證明**甲基訊號獨立於 phasing**。

**主要 deliverables**: 1 整合報告、4 MEMORY 修訂、1 PPT 勘誤、1 Filter prototype 跨樣本 F1 驗證。

---

## 1. 背景與證據彙整（為何需要這個 plan）

### 1.1 S5 觀察錯誤的鐵證

| 證據 | 內容 |
|------|------|
| Q4 30/30 cell 鐵證 | 6 樣本 × 5 CN tier 全部 AF-Mid TP > AF-Extreme TP，HCC1395 在 Gain/High 差距高達 +45–55pp |
| Q6 Inner-Outer gap | NearHalf 6/6 樣本 +5.2 至 +42.1pp（最強），Extr_lo 5/6 反向，Extr_hi 普遍微弱 |
| Q2 HCC1395 高 TP 區 | 真正 TP≥0.88 cell 是 Inner × NearHalf × Loss(0.980) / Sub(0.958) 與 Inner × Inter_lo × Sub-Diploid-Gain(0.898–0.929)，**沒有任何 Inner × Extreme AF cell 達 0.88** |

### 1.2 obs18 反轉發現

| AF bin | NG=2 Inner-Outer gap (median) | Direction |
|--------|------------------------------:|-----------|
| Extr_lo (<0.1) | **+23.35pp** | 強（NG=2 condition 後恢復） |
| Inter_lo (0.1–0.4) | **+29.10pp** | 最強 |
| NearHalf (0.4–0.6) | +21.70pp | 強 |
| Inter_hi (0.6–0.9) | +14.95pp | 中 |
| Extr_hi (≥0.9) | **−3.45pp** | 失效 |

→ obs18 「LOH-constrained phasing」結論成立區間應限定為 **AF<0.9**

### 1.3 HP_Ratio + 甲基的 AUC ceiling

| Cell | HP | Methyl | Combo | uplift |
|------|---:|-------:|------:|-------:|
| Inner × NearHalf × Sub | 0.770 | 0.628 | **0.795** | +0.025 |
| Inner × NearHalf × Loss | 0.655 | 0.588 | 0.692 | **+0.037** |
| Inner × Inter_hi × Gain | 0.680 | 0.593 | 0.673 | -0.007 |
| Inner × Inter_hi × High | 0.650 | 0.602 | 0.679 | +0.029 |
| HCC1395_DORADO Inner × NearHalf × Loss | 0.832 | 0.739 | **0.914** | +0.082 |
| HCC1954 Inner × NearHalf × Loss | 0.652 | **0.767** | 0.804 | +0.152 |

### 1.4 LongPhase 版本對比鐵證（HCC1395）

| 指標 | baseline_old | v5_new | Δ |
|------|------------:|------:|--:|
| LOH% | 58.7% | 62.2% | +3.5pp |
| Inner HP_Ratio_p25 | **0.962** | **0.036** | -0.926（劇變）|
| LOH 翻轉 region | — | — | 24,379 (38.5%) |
| HPFineNGroups 改變 | — | — | 38,082 (60.1%) |
| 甲基 6 特徵 AUC | 0.507–0.576 | 0.506–0.584 | **±0.008（不動）** |

---

## 2. Plan A: Erratum + 整合彙總報告（P0, ~3h）

### 2.1 目標
產生一份結構化的 erratum 報告，作為「2026-04-25 至 04-26 連續發現」的 canonical 紀錄，連結所有 figures、tsv 數據、與後續 plan B/C/D。

### 2.2 輸出位置
`InterSubMod/docs/experiments/in_progress/2026/04/20260426_HP_Methyl_3D_Integration_01.md`

### 2.3 必要章節結構
1. **Frontmatter**（id, date, sample list, mode=to, hypothesis_ids, tier）
2. **§1 Erratum**：S5 觀察錯誤的根因物理 + 三鐵證表
3. **§2 obs18 × AF bin 反轉**：NG=2 Inner-Outer gap × 5 AF bin 完整表 + 物理詮釋
4. **§3 HP_Ratio 二級切分**：14 高一致性 cell 表 + 失效區清單
5. **§4 LongPhase 版本對比**：5 版本 HP_Ratio/LOH%/NGroups 表 + Inner HP_Ratio 雙峰物理
6. **§5 甲基訊號 phasing-independent**：baseline vs v5 甲基 AUC 不動表 + AUC ceiling 0.79
7. **§6 對既有結論的影響**：列出每條受影響 MEMORY 與簡報的修訂建議
8. **§7 下一步**：連結 Plan B/C/D
9. **§8 Artifacts**：所有 .png .tsv 路徑（用 InterSubMod/ 前綴）

### 2.4 執行步驟（給未來執行者的精確指令）
```bash
# Step 1: 收集 figures 路徑
ls InterSubMod/docs/experiments/in_progress/2026/04/figures/2026042{5,6}_*

# Step 2: 收集 tsv 數據
ls /tmp/s5_quant/*.tsv  # 需要先把它移入 InterSubMod/research/.../data/
ls InterSubMod/docs/experiments/in_progress/2026/04/figures/2026042{5,6}_*/

# Step 3: 從本 plan §1 與本輪 session 對話直接複製鐵證表
# Step 4: 用 doc-standards skill 確認檔名/元數據
# Step 5: 寫入 docs/experiments/INDEX.md 新一行
```

### 2.5 驗收標準
- [ ] §1-§8 全部存在，每節有對應數據表/圖
- [ ] 所有 .md/.png/.tsv 路徑用 `InterSubMod/...` 前綴
- [ ] 連結到 Plan B/C/D（本檔案）
- [ ] docs/experiments/INDEX.md 新增 1 行
- [ ] verdict 章節：S5 erratum 為 **CONFIRMED CORRECTION**，combo AUC ceiling 為 **POSITIVE pilot (HCC1395)**

### 2.6 風險
- 對 obs18 「LOH-constrained phasing」主軸的限定（AF<0.9）可能與 paper 主軸衝突 → 需在報告中明確「不撤回主軸，但加註適用範圍」
- 若 erratum 與既有 PPT 的 S22 「TP gap +0.37」表述衝突，需說明這個 gap 是 aggregated（包含 AF<0.9 主導 + AF≥0.9 微弱），統計上仍正確

---

## 3. Plan B: HCC1395_DORADO outlier CV 驗證（P1, ~1h）

### 3.1 目標
檢驗 HCC1395_DORADO Inner × NearHalf × Loss combo AUC = 0.914 是否為真實訊號或 overfit。

### 3.2 動機
- 這是所有 sample × cell 組合中最高 combo AUC
- n=201（小樣本，overfit 風險高）
- 若真實，是 **同一基因組（HCC1395）不同 basecaller（DORADO vs Guppy）** 的對比 — DORADO 為什麼能達 0.914 而 Guppy 只到 0.724？

### 3.3 執行步驟
```python
# 1. Stratified 5-fold CV (現在用的是 LR 內建 5-fold，要改 nested CV 避免 leakage)
# 2. Permutation test (n=1000) — H0: combo AUC = HP-only AUC
# 3. Feature importance (logistic regression coefficients)
# 4. 對比 HCC1395 (Guppy) 同 cell 的 feature 分佈，找 DORADO 特殊特徵
# 5. 限制 cell n>=100 重跑 — 看是否 robust to subsampling
```

### 3.4 輸出
`InterSubMod/docs/experiments/in_progress/2026/04/20260426_B_DORADO_NearHalf_Loss_Outlier_01.md`

### 3.5 驗收標準
- [ ] 5-fold CV AUC 95% CI（如 [0.85, 0.95]）報告
- [ ] Permutation p-value 報告（H0: combo = HP-only）
- [ ] Feature importance ranking
- [ ] HCC1395 vs HCC1395_DORADO 該 cell 的 feature 分佈對比表
- [ ] verdict: TRUE_SIGNAL / OVERFIT / INCONCLUSIVE

### 3.6 風險
- n=201 即使 5-fold CV 也可能 overfit（每 fold 約 40 samples）→ 改用 leave-one-region-out 或 stratified 10-fold
- DORADO basecaller 在低 CN region 可能有未知 systematic effect → 需查 knowledge base `02_samples/HCC1395.md` 是否有 DORADO 特殊紀錄

---

## 4. Plan C: MEMORY 修訂（P0, ~30min）

### 4.1 目標
落實本次發現到 MEMORY 系統，避免下次 session 重蹈覆轍。

### 4.2 必須執行的 4 條

#### C-1 新增 `feedback_inner_extreme_af_not_tp_pure.md`
- **type**: feedback
- **內容要點**：
  - Rule: 不要把 Inner × Extreme AF (<0.1 或 ≥0.9) 視為 TP-pure 區
  - Why: 2026-04-25 數據 6 樣本 × 5 CN tier (30/30 cell) 全部 AF-Mid TP > AF-Extreme TP；HCC1395 Gain/High 差距 +45–55pp
  - How to apply: 觀察 LOH × AF × CN heatmap 時，high-TP cell 在 NearHalf/Inter_lo × Loss/Sub，不是 Extreme AF；Extreme AF 含 germline-homozygous-on-retained-hap 污染

#### C-2 新增 `project_hp_ratio_secondary_discrimination_inner_midaf.md`
- **type**: project
- **內容要點**：
  - Fact: HP_Ratio 在 Inner × Mid AF (NearHalf/Inter_lo/Inter_hi) × Loss/Sub/Diploid/Gain 是有效 secondary discriminator (AUC 0.65–0.77, 6 樣本 4–5/6 一致)
  - Why: 機制是 LOH 區的 HP-imbalance 反映 retained hap 上的 read 分布，與 somatic vs germline-leak 區分相關
  - How to apply: 設計 filter scheme 時，在 high-TP cell (Inner × NearHalf × Loss/Sub) 內可加 HP_Ratio≥0.5 cutoff 進一步提純；Outer 與 Inner × Extr_hi 區無效

#### C-3 新增 `feedback_baseline_hp_ratio_unreliable.md`
- **type**: feedback
- **內容要點**：
  - Rule: 不要直接信任 baseline LongPhase-TO 的 HP_Ratio（self-phasing artifact）
  - Why: 2026-04-26 對比 baseline vs v5_somatic_fb，HCC1395 LOH% 58.7→62.2%（38.5% LOH 翻轉），Inner HP_Ratio_p25 從 0.962 變 0.036 揭露雙峰分佈；baseline 把 phasing 失敗區偽造成 HP_Ratio≈1
  - How to apply: 涉及 HP_Ratio/HP1FamilyN/HP2FamilyN 分析時必須用 v3fixed/v4altguard/v5 數據，不可用 baseline

#### C-4 新增 `project_methylation_phasing_independent.md`
- **type**: project
- **內容要點**：
  - Fact: 甲基特徵（HPMergedDelta, AlleleDelta, HPFineF, CramersV）的 TP/FP AUC 在 baseline vs v5_somatic_fb 完全不變（Δ ≤ 0.008），不受 phasing 修正影響
  - Why: 甲基訊號是 read-level 直接觀察，phasing 修正僅改變 HP-tag assignment，不改變 read 的 5mC 訊號
  - How to apply: 甲基特徵可作為 phasing-failed sample (HCC1954) 的 fallback discriminator；HCC1954 Inner × NearHalf × Loss combo AUC 0.804（甲基主導 0.767）

### 4.3 同時更新（加註而非覆寫）
- `project_loh_constrained_phasing_discovery.md`: 加註「機制成立區間限定 AF<0.9，AF≥0.9 區因 germline-homozygous-on-retained 污染失效」
- `project_loh_subclone_af_methylation_positive.md`: 加註「ΔNG=+0.705 是 NG ratio 而非 TP rate；TP 訊號區為 Mid AF (0.1-0.9)」
- `MEMORY.md`: 新增 4 條 entry（保持 <200 行）

### 4.4 驗收標準
- [ ] 4 個新 .md 檔產生於 `InterSubMod/.claude/projects/.../memory/`
- [ ] MEMORY.md 新增 4 條一行 entry
- [ ] 既有 2 條加註（loh_constrained_phasing_discovery, loh_subclone_af_methylation_positive）

---

## 5. Plan D: Combined Filter F1 驗證（P1, ~4-6h）

### 5.1 目標
用 6 個 TO 樣本驗證以下 filter 對 F1 的影響：

```
Filter X1: Pass = Inner × NearHalf × (Loss OR Sub) × HP_Ratio≥0.5
Filter X2: Pass = (Inner × NearHalf × Loss/Sub × HP_Ratio≥0.5) OR (Inner × Inter_lo × Sub × HP_Ratio≥0.5)
Filter X3: 加上甲基 LR score：Pass if combo_LR_prob >= optimal_threshold
```

### 5.2 動機
- 任務 1 證明 HP_Ratio 在 Inner × NearHalf × Sub 達 AUC 0.770
- 甲基 + HP_Ratio combo 達 AUC 0.795
- 但 AUC 是 region-level 區分能力，**實際 F1（precision × recall trade-off）需要量化**
- 是否能在不犧牲 recall 的情況下提升 precision？

### 5.3 執行步驟
```python
# Step 1: Load 6 sample × ground truth (TP/FP from ClairS-TO)
# Step 2: 套用 Filter X1/X2/X3 到每樣本
# Step 3: 計算 ΔF1 = F1(after filter) - F1(before filter), per sample
# Step 4: 分層分析 by sample purity (HCC1395=58%, HCC1954=22%, etc.)
# Step 5: 計算 cross-sample ΔF1 mean ± std + Wilcoxon signed-rank
# Step 6: 對 HCC1954 (phasing-fail) 改用 Filter X3'（甲基主導，跳過 HP_Ratio）
```

### 5.4 必要的 Filter scheme 比較表（範本）

| Sample | Baseline F1 | X1 ΔF1 | X2 ΔF1 | X3 ΔF1 | X3' ΔF1 (HCC1954-only) |
|--------|------------:|-------:|-------:|-------:|----------------------:|
| HCC1395 | TBD | TBD | TBD | TBD | — |
| HCC1395_DORADO | TBD | TBD | TBD | TBD | — |
| HCC1937 | TBD | TBD | TBD | TBD | — |
| HCC1954 | TBD | TBD | TBD | TBD | TBD |
| H2009 | TBD | TBD | TBD | TBD | — |
| H1437 | TBD | TBD | TBD | TBD | — |

### 5.5 輸出
`InterSubMod/docs/experiments/in_progress/2026/04/20260426_D_Combined_Filter_F1_Pilot_01.md`

### 5.6 驗收標準
- [ ] 6 樣本 × 4 filter scheme F1 量化表
- [ ] cross-sample ΔF1 中位數 + Wilcoxon p-value
- [ ] 跨 purity 層的 trend 分析（高 purity 增益 vs 低 purity）
- [ ] verdict: ACCEPTABLE_FILTER / NO_F1_GAIN / NEEDS_PER_SAMPLE_TUNING

### 5.7 風險
- AUC 0.79 不必然轉換為 F1 提升 — 之前 [InterSubMod/.claude/.../memory/project_paired_f1_filter_abandoned.md](.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/project_paired_f1_filter_abandoned.md) 顯示 paired F1-filter 已放棄，TO 可能類似結果
- HCC1954 phasing-fail 需要特殊處理，不能套用通用 filter
- 若 ΔF1 < 0.005 則為 NO_F1_GAIN，須回到 characterization-only 定位

---

## 6. 執行順序與依賴

```
Plan A (報告) ─┬─→ Plan C (MEMORY)  [P0 並行可]
               │
               ├─→ Plan B (DORADO outlier)  [P1 獨立可]
               │
               └─→ Plan D (Filter F1)  [P1 依賴 A 的數據檔]
```

### 建議排程
- **Day 1 (今天)**: Plan A + Plan C（半天，文件工作）
- **Day 2**: Plan B（單樣本 outlier）+ Plan D pilot 1 樣本驗證
- **Day 3-4**: Plan D 全 6 樣本 + cross-sample 整合

### Hard Gate 檢查點
- 完成 Plan A 後：用戶確認 erratum 是否影響 PPT (要不要更新 PPT outline S8 標籤錯誤)
- 完成 Plan D 後：若 ΔF1 < 0.005 → Hard Gate，回報為 NO-GO 等用戶決策是否轉 characterization-only

---

## 7. 受影響的文件清單

### 7.1 PPT 簡報
- [InterSubMod/docs/presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/02_ppt_slide_outline.md](../../../presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/02_ppt_slide_outline.md)
  - S8 「綠區 TP 88-96% (LOH×Extreme)」→ 「LOH × Mid AF 0.1-0.9 × Loss/Sub」
  - S22 加註「gap 主要來自 AF 中段，AF≥0.9 失效」
- [InterSubMod/docs/presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/CHANGELOG.md](../../../presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/CHANGELOG.md)
  - 新增 v3 erratum entry

### 7.2 Research Landscape Reports
- [InterSubMod/docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md](../../../reports/research_landscape/07_LOH_CN_AF_研究總整理.md) — A3 升級為已完成
- [InterSubMod/docs/reports/research_landscape/08_Zone_Aware.md](../../../reports/research_landscape/08_Zone_Aware.md) — 第 117 行 self-phasing artifact 與本次發現吻合，加引用
- [InterSubMod/docs/reports/research_landscape/10_Research_Chain_Registry.md](../../../reports/research_landscape/10_Research_Chain_Registry.md) — Z1-Z5 zone 5-axis matrix 加 5-bin 細化

### 7.3 MEMORY (Plan C)
- 新增：4 個 .md 檔
- 加註：2 個 .md 檔
- 索引：MEMORY.md 4 行

### 7.4 Experiments Index
- [InterSubMod/docs/experiments/INDEX.md](../../../experiments/INDEX.md) — 新增 4 行（A/B/D 報告 + 本 plan）

### 7.5 CURRENT_FOCUS
- [InterSubMod/docs/CURRENT_FOCUS.md](../../../CURRENT_FOCUS.md) — 更新「進行中」事項：S5 erratum 完成、HP_Ratio 二級切分驗證

---

## 8. Artifacts 索引（已產出，連結到 Plan）

### 圖表
| 編號 | 路徑 | 描述 |
|------|------|------|
| F-25a | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_s5b_hp_ratio_color/fig_s5b_hp_ratio_scatter.png](../../../experiments/in_progress/2026/04/figures/20260425_s5b_hp_ratio_color/fig_s5b_hp_ratio_scatter.png) | S5b: HP_Ratio 著色散點 |
| F-25b | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_s5b_hp_ratio_color/fig_s5c_hp_ratio_hist.png](../../../experiments/in_progress/2026/04/figures/20260425_s5b_hp_ratio_color/fig_s5c_hp_ratio_hist.png) | S5c: HP_Ratio 直方圖 (HCC1954 雙峰) |
| F-25c | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_v2_1_corrected/fig_v2_1_to_tp_heatmap_corrected.png](../../../experiments/in_progress/2026/04/figures/20260425_v2_1_corrected/fig_v2_1_to_tp_heatmap_corrected.png) | fig_v2_1 修正版 (high-TP 綠框正確位置) |
| F-25d | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_obs18_af_split/obs18_NG_AF_gap_lines.png](../../../experiments/in_progress/2026/04/figures/20260425_obs18_af_split/obs18_NG_AF_gap_lines.png) | obs18 NG × AF gap 線圖 |
| F-26a | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_hp_ratio_3d_auc/hp_ratio_auc_heatmap.png](../../../experiments/in_progress/2026/04/figures/20260426_hp_ratio_3d_auc/hp_ratio_auc_heatmap.png) | HP_Ratio AUC × LOH × AF × CN heatmap |
| F-26b | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/hp_ratio_version_compare.png](../../../experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/hp_ratio_version_compare.png) | baseline vs v5 HP_Ratio 分佈 |
| F-26c | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/version_3d_tp_heatmap.png](../../../experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/version_3d_tp_heatmap.png) | baseline vs v5 LOH×AF×CN TP heatmap |
| F-26d | [InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/auc_ceiling_per_cell.png](../../../experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/auc_ceiling_per_cell.png) | HP_Ratio + Methyl AUC ceiling bar chart |

### 數據檔（TSV）
| 路徑 | 描述 |
|------|------|
| `/tmp/s5_quant/t1_inner_outer_AF08.tsv` | Stage1 AF<0.8 Inner-Outer gap |
| `/tmp/s5_quant/t1c_5bin_inner_outer.tsv` | 5 AF bin Inner-Outer per sample |
| `/tmp/s5_quant/t_corr1_inner_extreme_af_by_cn.tsv` | Inner × Extreme AF × CN tier |
| `/tmp/s5_quant/t_corr2_full_3d_HCC1395.tsv` | HCC1395 LOH × AF × CN full grid |
| `/tmp/s5_quant/t_corr4_inner_extreme_vs_mid.tsv` | Inner Mid vs Extreme AF (Q4) |
| `/tmp/s5_quant/t_corr6_inner_outer_gap_by_af.tsv` | Inner-Outer gap × AF (Q6) |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_obs18_af_split/obs18_NG2_gap_by_af.tsv` | obs18 NG=2 × AF gap |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260425_obs18_af_split/obs18_NG_by_af_summary.tsv` | obs18 by NG×AF medians |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_hp_ratio_3d_auc/hp_ratio_auc_3d_per_sample.tsv` | HP_Ratio AUC per sample×cell |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_hp_ratio_3d_auc/hp_ratio_auc_cross_sample.tsv` | HP_Ratio cross-sample consistency |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/version_summary.tsv` | 5 版本全局統計 |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_longphase_old_vs_new/hp_ratio_auc_by_version.tsv` | HP_Ratio AUC per version |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/methyl_high_cells_top.tsv` | 甲基 high-cell consistency |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/hp_methyl_lr_combo.tsv` | HP+Methyl LR uplift |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/auc_ceiling_summary.tsv` | AUC ceiling 跨樣本 |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260426_methylation_3d_addon/methyl_baseline_vs_v5.tsv` | 甲基 baseline vs v5 不變對比 |

### 分析腳本（暫存於 /tmp/，建議移入 scripts/analysis/）
- `/tmp/s5_quant_analysis.py` — Stage 1-3 分層 + secondary feature + HCC1954 outlier
- `/tmp/s5_correction.py` — LOH × AF × CN 三維修正驗證
- `/tmp/fig_s5b_hp_ratio_per_sample.py` — fig_s5b/c 產生
- `/tmp/fig_v2_1_redrawn.py` — fig_v2_1 重畫
- `/tmp/obs18_af_bin_split.py` — obs18 × AF bin 拆分
- `/tmp/hp_ratio_3d_auc.py` — HP_Ratio 3D AUC
- `/tmp/longphase_version_compare.py` — LongPhase 5 版本對比
- `/tmp/methylation_3d_analysis.py` — 甲基訊號延伸分析

**Plan A 第一步必做：將 /tmp/ 腳本與 .tsv 移入 InterSubMod/scripts/analysis/ 與 InterSubMod/research/.../data/，避免下次 reboot 遺失。**

---

## 9. 證據鏈（連到 evidence_ledger / hypothesis_queue）

### 9.1 假設 ID 對應
- H_S5_erratum: Inner × Extreme AF ≠ TP-pure (本輪 NEW，待寫入 hypothesis_queue.json)
- H_HP_Ratio_secondary: HP_Ratio 在 Inner × Mid AF × Loss/Sub 是 secondary discriminator (NEW)
- H_methyl_phasing_indep: 甲基 AUC ≠ phasing version (NEW)
- H_baseline_self_phasing: baseline LongPhase HP_Ratio 是 self-phasing artifact (與既有 project_self_phasing_causal_chain_confirmed 共證)

### 9.2 Tier 評估
- S5 erratum: **Tier 5**（6 樣本 × 5 CN tier 30/30 一致）
- HP_Ratio secondary: **Tier 4**（pilot PASS 6 樣本一致 4-5/6 cell）
- 甲基 phasing-independent: **Tier 4**（5 版本對比，但僅 HCC1395）— Plan B 跨樣本驗證可升 Tier 5
- Filter F1: **Tier 0** 待 Plan D 驗證

### 9.3 evidence_ledger 寫入建議
Plan A 完成後，在 [InterSubMod/research/autoresearch/evidence_ledger.jsonl](../../../research/autoresearch/evidence_ledger.jsonl) 加 4 行 JSON entry（每個假設 1 行）。

---

## 10. 風險與不確定性

### 10.1 Plan-level 風險
| Plan | 主要風險 | 緩解 |
|------|---------|------|
| A | obs18 主軸與 erratum 衝突 | 加註「不撤回主軸，限定 AF<0.9」 |
| B | n=201 overfit | 用 leave-one-region-out CV + permutation |
| C | MEMORY.md 超 200 行 | 同時 consolidate concluded entries |
| D | ΔF1 < 0.005 | Hard Gate，回報 NO-GO |

### 10.2 跨 Plan 系統性風險
- 本輪 archive 數據是 post-HP-fix（即 v3fixed/v5），但**並非全 6 樣本都重跑過 v5**（非 HCC1395 樣本可能是 v3fixed）→ 需在報告中明示
- HP_Ratio AUC 0.65-0.77 受限於 archive 數據 quality，pilot 結論 Tier 4 不能直接升 Tier 5

### 10.3 用戶決策需求
- Plan A 完成後：是否同步更新 PPT outline S8 標籤錯誤？
- Plan D 結果若 NO-GO：研究方向是否轉 characterization-only（如 paired_f1_filter_abandoned 一致）？

---

## 11. 成功標準（Plan-level）

| Plan | 成功 = | 失敗 = |
|------|-------|-------|
| A | 整合報告產出 + INDEX 更新 + 連結 4 條 MEMORY 修訂 | 報告未涵蓋 §1-§8 全部 |
| B | DORADO outlier verdict 確定（TRUE_SIGNAL / OVERFIT / INCONCLUSIVE） | CV 無法收斂 |
| C | 4 新 MEMORY + 2 加註 + MEMORY.md 更新 | MEMORY.md 行數 >200 或缺 entry |
| D | 6 樣本 × 4 scheme ΔF1 量化表 + verdict | 無法 reproduce baseline F1 |

---

## 12. 給未來執行者的快速指南

### 12.1 進入這個 plan 前必讀
1. 本檔（本 plan）
2. [InterSubMod/docs/experiments/in_progress/2026/04/20260423_B4_S4_Secondary_Discrimination_01.md](../../../experiments/in_progress/2026/04/20260423_B4_S4_Secondary_Discrimination_01.md)（先前 S4 secondary 工作）
3. [InterSubMod/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/project_loh_constrained_phasing_discovery.md](../../../../.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/project_loh_constrained_phasing_discovery.md)
4. [InterSubMod/docs/presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/01_full_narrative_report.md](../../../presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/01_full_narrative_report.md)

### 12.2 第一個 5 分鐘要做什麼
```bash
# 移植 /tmp/ 腳本到永久位置
cp /tmp/s5_quant_analysis.py InterSubMod/scripts/analysis/20260425_s5_quant_analysis.py
cp /tmp/methylation_3d_analysis.py InterSubMod/scripts/analysis/20260426_methylation_3d_analysis.py
# ... 同上對 8 個腳本

# 移植 /tmp/s5_quant/*.tsv 到 research/
mkdir -p InterSubMod/research/s5_erratum_2026_04/data
cp /tmp/s5_quant/*.tsv InterSubMod/research/s5_erratum_2026_04/data/

# 確認所有 figures 路徑可解（用 InterSubMod/ 前綴）
find InterSubMod/docs/experiments/in_progress/2026/04/figures/2026042* -name "*.png" | wc -l
# 預期: 8+
```

### 12.3 第一個 30 分鐘要做什麼
- 寫 Plan A 的 §1 Erratum 章節（直接從本 plan §1 複製鐵證表）
- 寫 Plan C 的 4 個 MEMORY 檔（參考既有 memory 檔範本格式）

### 12.4 預期完成日
- Plan A + C: 2026-04-26（今天）
- Plan B: 2026-04-27
- Plan D: 2026-04-28 至 04-30

---

## 13. 附錄 A：HCC1954 PON-only Phasing 資源（2026-04-26 盤點完成）

當 Plan B 或後續延伸需要重做 HCC1954 phasing 時，下列資源已備齊：

### 13.A.1 HCC1954 Tagged BAM（既有 canonical post-HP-fix）
- `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_pileup/20260315_HCC1954_paired_pileup_pileup_complete_matrix/longphase_s/HCC1954_tagged.bam`
- `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/longphase_s/HCC1954_tagged.bam`

### 13.A.2 HCC1954 原始 BAM（Google_somatic_data，可重做 phasing）
- Tumor: `/big8_disk/Google_somatic_data/bams/HCC1954/HCC1954_Tumor_ONT.GRCh38.sorted.bam`
- Normal: `/big8_disk/Google_somatic_data/bams/HCC1954/HCC1954_Normal_ONT.GRCh38.sorted.bam`
- Tumor (alt): `/big8_disk/data/HCC1954/ONT/HCC1954.bam`
- Normal (alt): `/big8_disk/data/HCC1954/ONT/HCC1954BL.bam`

### 13.A.3 HCC1954 subsample 系列（純度梯度可用於 purity-vs-phasing-fail 因果驗證）
- `/big8_disk/data/HCC1954/ONT/subsample/t{10,20,30,40,50}_n{00,10,20,25,30,40}/HCC1954_t*_n*.bam`
- 應用：若 HCC1954 phasing-fail 是低 purity 問題，t50_n00 (高 purity) 應顯著改善 HP_Ratio

### 13.A.4 LongPhase-TO 三版本（均位於 /big7_disk/liaoyoyo2001/longphase-to-mod/）
- `longphase-to-baseline` — 舊 HP（self-phasing 受污染）
- `longphase-to` — 當前主版（含 PON-only 模式）
- `longphase-to-v3fixed` — V3 fixed 版

### 13.A.5 既有 HCC1954 outlier 調查文件（必讀）
- [InterSubMod/docs/experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md](../../../experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md) — 早期 reversal 調查
- [InterSubMod/docs/experiments/in_progress/2026/04/20260423_B2_HCC1954_Outlier_RootCause_01.md](../../../experiments/in_progress/2026/04/20260423_B2_HCC1954_Outlier_RootCause_01.md) — 04-23 root cause 草稿（待本次發現注入）
- [InterSubMod/docs/experiments/in_progress/2026/04/20260411_PON跨樣本移除率驗證_01.md](../../../experiments/in_progress/2026/04/20260411_PON跨樣本移除率驗證_01.md) — PON 跨樣本驗證

### 13.A.6 啟動 HCC1954 PON-only phasing 的最小腳本骨架
參考 `/big7_disk/liaoyoyo2001/longphase-to-mod/run_v4_ism_comparison.sh`（HCC1395 版），改：
```
SAMPLE_NAME="HCC1954"
TUMOR_BAM="/big8_disk/Google_somatic_data/bams/HCC1954/HCC1954_Tumor_ONT.GRCh38.sorted.bam"
NORMAL_BAM="/big8_disk/Google_somatic_data/bams/HCC1954/HCC1954_Normal_ONT.GRCh38.sorted.bam"
# 替換 ClairS-TO TP/FP VCF 路徑為 HCC1954 對應
```

預估 wall-clock：LongPhase-TO PON-only ~1.5h + Haplotag ~30min + ISM TP+FP ~1h = **總 3h**

### 13.A.7 現有 phasing pipeline 腳本（2026-04-27 補強盤點）

InterSubMod/scripts/analysis/ 下既有腳本可重用或改寫：
- `InterSubMod/scripts/analysis/run_longphase_to_intersubmod_pilot.sh` — LongPhase-TO + ISM 完整 pilot 流程
- `InterSubMod/scripts/analysis/resume_to_pilot_from_tagged_outputs.sh` — 從 tagged BAM resume（節省時間）
- `InterSubMod/scripts/analysis/resume_paired_run_from_tagged_outputs.sh` — paired 模式 resume
- `InterSubMod/scripts/analysis/resume_failed_complete_matrix_paired.sh` — 失敗 region 重做
- `InterSubMod/scripts/analysis/run_purity_and_standard_verification.sh` — purity 驗證（與 13.A.3 subsample 系列搭配）
- `InterSubMod/scripts/analysis/run_complete_big7_experiments.sh` — 跨樣本完整實驗
- `InterSubMod/scripts/analysis/run_pure_research_round.sh` — 純研究 round
- `InterSubMod/scripts/analysis/build_purity_stage_metrics_bcftools.sh` — purity stage 指標
- `InterSubMod/scripts/run_vcf_all_snv.sh` + `run_batch_vcf_analysis.sh` — 主流程

**Plan B 推薦最小路徑**：用 `resume_to_pilot_from_tagged_outputs.sh` 從 13.A.1 既有 tagged BAM 重做 ISM TP/FP，避免重做 LongPhase（節省 ~1.5h，總 wall-clock 縮為 1.5h）

---

## 14. 附錄 B：原始對話可追溯到的關鍵決策點

| 決策點 | 對話節點 | 影響 |
|--------|----------|------|
| 用戶質疑「Inner Extreme AF = TP-pure」是否錯誤 | 2026-04-25 user msg #4 | 觸發 Plan A erratum |
| 用戶要求 HP_Ratio × LOH × AF × KDE 分析 | 2026-04-26 user msg #6 | 觸發 Plan B/C/D |
| 用戶要求加上甲基訊號 | 2026-04-26 user msg #7 | 完成甲基層分析（已併入 Plan A） |
| 用戶要求「整理寫成計劃報告.md」 | 2026-04-26 user msg #8 | 本檔產生 |

---

**End of Plan**
