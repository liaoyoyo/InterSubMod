# B.2-1 HCC1954 反向排除標準調查

**日期**：2026-04-17
**來源**：Plan Part B.2-1（對結論 11 LOH Subclone AF × Methylation POSITIVE 的質疑）
**相關**：
- `docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md`
- `research/loh_subclone_af/` + `research/loh_subclone_af_paired/`

---

## 0. 質疑背景

原結論 11：TO ΔNG=+0.705 / Paired ΔNG=+0.787 @ 7/7 p<10⁻³⁹、Spearman ρ=0.270-0.382

**發現的矛盾**：
- Step 2 (AF class 比較)：HCC1954 Δ=+0.69 POS (p=2.2e-39)、direction "Inter > Ext"
- Step 3 (segment-level Spearman)：**HCC1954 TO ρ=-0.297 / Paired ρ=-0.211 方向反向**（且不顯著，p>0.05）
- n_segments HCC1954 = 34 (TO) / 30 (Paired)，遠低於其他樣本（105-435）

**另一交叉矛盾**（B.1-1 新發現）：
- HCC1954 在 HPFineNGroups TO NonLOH 是 **7 樣本中唯一 ROBUST**（ΔAUC=+0.022）
- HCC1954 在 LOH×AF step3 是 **7 樣本中唯一反向**
- → 是否 HCC1954 生物學特性造成 "LOH 外 NGroups 強、LOH 內 AF-NG 關聯反向"？

**原報告的 post-hoc 解釋**：
- HCC1954 "純度/CN 複雜" — 質疑：post-hoc rationalization，未 pre-registered

---

## 1. 質疑點（逐條驗證）

| # | 質疑 | 驗證方法 |
|---|------|----------|
| Q1 | HCC1954 step3 ρ=-0.30 是**真反向**還是**小樣本雜訊**？ | Bootstrap 95% CI for Spearman ρ per sample |
| Q2 | HCC1954 n_segments=34 的來源？（少 LOH / 少符合 NumReads / 少 NGroups≥2） | 全程 filter chain 重建 + 樣本間比較 |
| Q3 | 移除 HCC1954 後「7/7 一致」是否變成「6/7」？ | meta-analysis forest + sensitivity |
| Q4 | HCC1954 step2 POS vs step3 NEG — 生物學能一致解釋嗎？ | Segment-level vs region-level 殘差分析 |
| Q5 | HCC1954 的特殊性是否可事先預測（pre-registered exclusion criterion）？ | 知識庫 + CN/ploidy/LOH 分佈比較 |

---

## 2. 資料與程式溯源（Git 可追溯）

### 資料
- `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`
- LOH.bed paths in `research/loh_subclone_af/scripts/step3_spatial_analysis.py:38-44`
- 既有統計：
  - `research/loh_subclone_af/data/step3_per_sample_segment_correlation.tsv`（TO）
  - `research/loh_subclone_af_paired/data/step3_per_sample_segment_correlation.tsv`（Paired）
  - `research/loh_subclone_af/data/step2_per_sample_consistency.tsv`
  - `research/loh_subclone_af_paired/data/step2_per_sample_consistency.tsv`

### Git baseline
- 起始 HEAD：`ab61ad1`（B.1-1 完成後）
- 本次分析不改 C++，只新增 `research/hcc1954_reversal_investigation/`

---

## 3. 方法（Step → Verify）

| Step | 方法 | 驗證標準 |
|------|------|---------|
| S1 | 讀 step3 segment-level TSV + master rows，重建 HCC1954 的 segments 清單 | segments n 吻合 34 (TO) / 30 (Paired) ✅ |
| S2 | Bootstrap 95% CI for Spearman ρ 每個樣本（B=2000） | 檢查 HCC1954 CI 是否含 0（噪音）或不含 0（真反向）✅ |
| S3 | HCC1954 segments filter chain 拆解：（a）LOH segments（b）有 ≥2 regions（c）有 NGroups≥2 | 找到 n_segments=34 的流失點 ✅ |
| S4 | Sensitivity：移除 HCC1954 後 meta-analysis ρ̄ 與 p；對比原 7/7 | 若 6/6 仍顯著 → HCC1954 可 pre-register exclusion ✅ |
| S5 | HCC1954 step2 vs step3 調和：region-level 內 AF→NG ρ（step3 用 segment 平均，step2 用全 region）| 若 region-level ρ 為正 → step3 平均效應翻轉是聚合 artifact ✅ |
| S6 | 知識庫對照：HCC1954 HER2+ hyper-diploid 特性 vs 其他樣本 | 書面比較 ploidy / LOH Mb / segment 數，看是否能事先定 exclusion ✅ |

**成功標準**：
- HCC1954 ρ 95% CI 含 0 → 反向為統計噪音 → 結論 11 實為 **7/7 方向一致 + 1/7 功效不足**
- HCC1954 ρ 95% CI 不含 0 → **真反向** → 需生物學 pre-registered 排除標準
- 任一情況下：移除 HCC1954 的敏感性分析必須做

---

## 4. 預期產出

- `scripts/01_hcc1954_reversal_check.py`（主分析）
- `data/rho_bootstrap_ci.tsv`（7 樣本 × {TO, Paired} × [ρ, CI_low, CI_high, n_segments]）
- `data/hcc1954_filter_chain.tsv`（HCC1954 segments 流失拆解）
- `data/sensitivity_leave_one_out.tsv`（移除每個樣本後的 meta ρ）
- `data/hcc1954_region_level_rho.tsv`（step2 vs step3 調和）
- `figures/01_bootstrap_ci_forest.png`（7 樣本 ρ + CI）
- `figures/02_hcc1954_segments_distribution.png`（HCC1954 vs 其他樣本 segment 數分佈）
- 更新報告：`docs/experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md`
- 更新結論 11 穩定度（補充結論或修正穩定度評分）
