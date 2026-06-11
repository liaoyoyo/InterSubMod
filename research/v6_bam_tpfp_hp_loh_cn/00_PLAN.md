# V3F/V5/V6 ISM 三向 × LOH.bed × Coverage_Multiple × HP bucket — 後分析計畫

> 計畫書版本：v0.3（2026-05-15）
> 目的：在現有 V3F / V5 / V6 三向 ISM 結果（已產出於 `research/paired_priority_bug_audit/phaseC_genome_three_way/` + `phaseD_v6_5sample/`）之上，系統量化 LOH.bed × HP bucket × Coverage_Multiple 三軸對 TP/FP 的區分能力，並比較三個 longphase haplotag 版本之間的訊號變化。
> 範圍：**characterization-only**（不評估 filter / ΔF1）。HCC1395 pilot 確認現象 → 已有 4 樣本擴展（H1437/H2009/HCC1954/HCC1937）。
> **不需重跑任何 BAM 或 ISM**：所有上游資料已存在，本計畫是 post-hoc 後分析。
> 對應 CURRENT_FOCUS：Phase 1 W1 F-paired-D3（V5→V6 ISM 影響量化）的延伸與 Thread D 主軸的三軸擴充。

## v0.3 重要修正

| 版本 | 認知 | 現實 |
|------|------|------|
| v0.1 | 需重跑 baseline/V5/V6 三方 BAM | V3F/V5/V6 三向 BAM 與 ISM 結果都已存在 |
| v0.2 | 5 軸 1000-cell grid | 降為 3 軸 50-cell grid + NG/AF covariate |
| v0.3 | baseline = upstream longphase | 改 V3F (PON-only V3-Fixed) 為三方對照組（與 phaseC 命名一致） |

## 關鍵資料 inventory（v0.3 新增）

| 資料 | 路徑 | 狀態 |
|------|------|------|
| V3F BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam` | ✅ 已存在 |
| V5 BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_tagged.bam` | ✅ 已存在 |
| V6 BAM (268 GB) | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam` | ✅ 已存在 |
| HCC1395 V3F/V5/V6 三向 ISM (12 runs) | `research/paired_priority_bug_audit/phaseC_genome_three_way/{V3F,V5,V6}_{on,off}_{tp,fp}/` | ✅ 已存在 |
| Three-way summary | `phaseC_genome_three_way/v3f_vs_v5_vs_v6_genome_summary.tsv` | ✅ 已存在 |
| Three-way per-region NG | `phaseC_genome_three_way/v3f_vs_v5_vs_v6_region_ng.tsv` (3.7 MB) | ✅ 已存在 |
| 4 樣本 V6 ISM | `phaseD_v6_5sample/{H1437,H2009,HCC1954,HCC1937}/{on,off}_{tp,fp}/` | ✅ 已存在 |
| Cross-sample summary | `phaseD_v6_5sample/v6_cross_sample_summary.tsv` | ✅ 已存在 |
| COLO829 | deferred (truth set 0600 權限) | 🟡 待解 |
| V6 caller F1 驗證 | `09_V6_caller_F1_verification.md` (5/12) | ✅ 已驗 V6 caller F1 不變（重用 V5 phased VCF） |

---

## Context — 為什麼現在做這個

### 歷史警訊（須避免重蹈覆轍）

| 嘗試 | 結論 | 穩定度 | 失敗根因 |
|------|------|--------|---------|
| LOH binary filter 10/10 策略 | NEGATIVE | ⭐4 | LOH 反而 TP-enriched，filter 方向錯 |
| CN Zone-aware filter | 跨 7 樣本 mean AUC ≤ 0.641 | ⭐4 | HCC1395 chr8 樣本特異性 |
| LOH+HPMergedSig 7.4× | 87.5% 來自 HCC1395 chr8，崩塌 1.3× | ⭐4 | 樣本/染色體特異性 |
| TO 60+ 特徵 AUC<0.64 | NEGATIVE | ⭐4 | germline FP ≈ somatic TP in methylation space |

**這些都是 "filter" 失敗，不是 "characterization" 失敗。** 本計畫嚴守 characterization 邊界。

### 3 個新變量讓 HP × LOH × CN 三軸交叉值得重做

1. **V6 phased BAM**（5/10 commit 鏈完成）— V5 在 germline-absent 區域繼承 4.19:1 priority bug 偏移（`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md` §8.6）；V6 修補後 HP tag 更乾淨。**從未量化過 V5→V6 在 LOH.bed 內外的 TP/FP 區分能力差異**。
2. **Thread D 已有「LOH 內外 × HP 4-bucket」框架**（`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`）— NG=2 Inner same-hap TP rate +0.37 gap（6/6 樣本 Wilcoxon p=0.0156）。但 Thread D 沒加入 CN 軸，也沒做 FP-rich 區域 zoom-in。
3. **Coverage_Multiple KDE-corrected** 已修並驗證（r=0.831 vs SEQC2，commits `374fad4`+`12d9b3e`），**幾倍體軸從未與 LOH 內外 × HP bucket 交叉**。

### 預期 outcome

- 對 V6 BAM 量化 LOH.bed 內/外 × HP 4-bucket × CN 區間（4-6 bins）的 TP/FP rate 立體格網（HCC1395 chr19+全基因體）
- 找出 TP/FP 區分最強的 cell（候選 characterization marker）與最弱的 cell（confounds）
- 比較 V5 vs V6 同樣 grid 的差異 → 量化 self-phasing 修補對 read-level signature 的純增益
- HCC1395 chr8 7.4× FP hotspot + Gain+LOH zone + Outer cross_het + 自動偵測 top 5% FP regions 4 種定義交叉比對 → 看哪一種 zone 定義可泛化

---

## 假設

| ID | 假設 | 驗證指標 | 預期方向 |
|---|------|---------|--------|
| H1a | V5 BAM Inner same-hap TP gap > baseline（Layer 1.5 修補增益） | Inner − Outer TP rate Δ(V5 − baseline) | ≥ 0.03 |
| H1b | V6 BAM Inner same-hap TP gap > V5（priority bug fix 增益） | Inner − Outer TP rate Δ(V6 − V5) | ≥ 0.03 |
| H1c | V6 對 baseline 累積增益最大 | Inner − Outer TP rate Δ(V6 − baseline) | ≥ 0.06 |
| H2 | LOH.bed 內 × HP1+HP1S same-hap cell 對 TP 富集 | TP rate vs 全域 baseline | +0.10 ~ +0.30 |
| H3 | Coverage_Multiple ≥ 1.5（Gain+LOH）× Outer cross_het 對 FP 最富集 | FP rate vs 全域 baseline | ≥ 2.0× enrichment |
| H4 | chr8 hotspot 的 TP/FP 區分主要由 LOH+CN 兩軸貢獻（HP bucket 增量 < 0.05） | 5 軸 LR ablation: full vs HP-removed | HP 增量 ≤ 0.05 |
| H5 | Auto-detect top 5% FP regions 與 4 個 known zones 重疊率 ≥ 60% | Jaccard | ≥ 0.60 |
| H6 | 1000-cell grid 中 powered cells (n≥50) 比例 ≥ 30% | Power curve | 確認 grid 不過 sparse |
| H7 | Top 20 extreme cells 中 ≥ 5 個通過所有 confound guard | within-group OLS + permutation + spatial | 找出真實 multivariate effect |

**全部結論限定 HCC1395 pilot**，跨樣本擴展放 step4（依 KDE-only，無 SEQC2）。

---

## 4-step 研究流程

### Step 1：V3F vs V5 vs V6 三向 ISM 結果整合（HCC1395）

**目的**：在 phaseC_genome_three_way 已存在的 12 個 ISM 結果上做後分析，分離兩段修補貢獻 — V3F→V5（Layer 1.5 加入）vs V5→V6（Layer 1.5 移除回 V3F 保守標 hp=33）vs V3F→V6（兩段累積；V6 = V3F + marker engineering 改善）。

**三方版本定義（與 phaseC 命名一致）**：

| 版本 | longphase 設計 | 預期 haplotag 行為 |
|------|---------------|----------------|
| V3F (PON-only V3-Fixed) | upstream PON-only + 移除 Layer 1.5 + germline-absent → hp=33 保守 | 早期 baseline 設計；hp=1-1:hp=2-1 ratio = 1.138（接近中性）|
| V5 (Layer 1.5) | V3F + 加 Layer 1.5 somatic fallback (HaplotagProcess.cpp:537-548) | hp=1-1:hp=2-1 ratio = 1.86；germline-absent 區繼承 priority bug 偏 HP1（4.19:1） |
| V6 (production candidate) | V5 + 移除 Layer 1.5 回 V3F 保守標；重用 V5 phased VCF（phasing 層不變） | marker coverage 23,980 > V3F 21,997 > V5 18,382；hp=33 reads V6=138,317 > V3F=132,060；caller F1 不變 |

**已驗證事實（07_V6_validation_findings.md, 5/10）**：
- V6 hp=1-1:hp=2-1 ratio = 1.838（介於 V3F 1.138 與 V5 1.86 之間；germline-existent 區因重用 V5 phased VCF 殘留）
- V6 marker coverage 比 V3F +9.0%、比 V5 +30.5%
- caller F1 vs SEQC2 truth：V6 = V5（重用 phased VCF）
- 4 樣本（H1437/H2009/HCC1954/HCC1937）V6 跨樣本驗證已通過

**Step 1 後分析任務（read-only）**：
1. **讀現成 ISM 結果**：phaseC_genome_three_way 的 12 個 master.tsv（V3F/V5/V6 × on/off × tp/fp），分別 reconstruct 全基因組 region-level TSV
2. **三方比對**：對每個 region 計算 (V3F, V5, V6) 三點 HP tag 分布、HPFineNGroups、Coverage_Multiple、LOH.bed 標記
3. **per-region trajectory**：V3F → V5 → V6 是否 monotonic（5 類 region：兩段都改善 / 只 V5 改善 / 只 V6 改善 / 無改善 / 反向）
4. **off / on flag 對照**：`--germline-hp-only=off`（V5 Layer 1.5 active）vs `=on`（V5 Layer 1.5 mask 為 unphased）— 看 ISM 端 demotion 是否能近似 V6 的效果

**輸出**：`step1_v3f_v5_v6_three_way/`
- `step1_master_three_way.tsv`（V3F/V5/V6 × on/off × tp/fp = 12 個 ISM 結果整合至 region-level）
- `step1_delta_v5_vs_v3f.tsv` / `step1_delta_v6_vs_v5.tsv` / `step1_delta_v6_vs_v3f.tsv`
- `step1_trajectory.tsv`（每 region 三點 + monotonic flag + 5-類 region classification）
- `step1_three_panel_heatmap.png`（V3F / V5 / V6 三聯 heatmap）
- `step1_off_vs_on_compare.tsv`（flag mask 是否近似 V6）
- `step1_findings.md`（H1a/H1b/H1c 判定 + 兩段修補貢獻分解）

### Step 2：3 軸 grid + NG/AF covariate（HCC1395 — 降軸版）

**修正歷史**：v0.1 設 5 軸 1000-cell，但 Plan agent 審計指出 5 軸 effective dim ≈ 2.5-3（caller_af 與 LOH.bed 強相關、HPFineNGroups 與 HP bucket 重疊），1000 cells 大半 underpowered。v0.2 降為 3 軸主 grid + 2 軸 covariate。

**3 軸 主 grid 定義**：

| 軸 | Bins | 數量 | 來源欄位 |
|----|------|------|---------|
| LOH (LOH.bed) | inner, outer | 2 | `PhasingGraph.cpp:1817` 輸出 BED |
| HP bucket | same_HP1, same_HP2, cross_het, cross_het_inv, other | 5 | `LabelTest.cpp:265-302` 4-bucket + other |
| Coverage_Multiple (KDE) | cov_loss (<0.7), cov_normal (0.7-1.3), cov_elevated (1.3-1.8), cov_gain (1.8-2.5), cov_high_gain (≥2.5) | 5 | master TSV `Coverage_Multiple` |

**總 cell 數**：2 × 5 × 5 = **50 cells per BAM version per sample**

三方對照（baseline / V5 / V6）→ 每樣本 **3 個 50-cell grid = 150 cells**（V6 BAM 待 step1 路徑確認）

**Covariate 軸（不在 grid 內分 bin，但每 cell 內做 LR 控制）**：

| Covariate | 處理方式 | 用途 |
|-----------|---------|------|
| HPFineNGroups (NG) | 在 cell 內當 numeric covariate 加入 LR | 控制 somatic heterogeneity confound |
| caller_af | 在 cell 內當 numeric covariate 加入 LR | 控制 AF gradient confound |

**Cell-level LR model**（每 powered cell 跑一次）：
```
logit(P(TP)) = β0 + β1 · NG + β2 · caller_af + β3 · NumReads + ε
```
- 拿 `β0`（intercept after adjusting covariates）作為「cell 內固有 TP rate」
- 拿 NG / AF / NumReads 三個 coefficient 看每個 covariate 對 cell 內 TP 預測的貢獻
- Likelihood ratio test (LRT) 比較 full vs 移除每個 covariate

**Cell × covariate 完整 cross-tab（供窮舉觀察）**：
- 主 grid 50 cells 之外，每 cell 內額外輸出 NG/AF 2D contingency（5 × 4 = 20 sub-cells）
- 不在主分析，但放 `step2_cell_internal_breakdown.tsv` 供事後 zoom-in
- 等效於提供「全 1000 組合」的 hierarchical 觀察視角，但統計權威只給 50-cell 主 grid

**Power gate**（依 Plan agent 建議加入）：
- 預先 dry-run：算每 cell 預期 n（從 master TSV marginal 計算）
- **Gate 條件**：powered cells (n ≥ 50) ≥ 15 個（50 cells 的 30%）才進入完整 step2；否則合併 NG/cov 鄰近 bin 或降至 2 軸（LOH × HP bucket = 10 cells）

**Collinearity audit（必做，依 Plan agent 建議）**：
- 計算 5 軸（含 covariate）pairwise Cramér's V + VIF
- 若 caller_af vs LOH inner/outer Cramér's V > 0.6（預期會發生）→ 文件化但不刪 caller_af（它是 covariate 不在 grid 內 → 不會自我重疊）
- 若 HPFineNGroups vs HP bucket Cramér's V > 0.6 → 文件化

**三方 BAM 對照**：
- 每 BAM 版本獨立跑 50-cell grid + LR
- 3-point trajectory（baseline → V5 → V6）：per-cell intercept β0 變化 + monotonic flag
- 兩段 Δ heatmap：Δ1 = V5 β0 − baseline β0；Δ2 = V6 β0 − V5 β0
- 識別 5 類 cell（兩段都改善 / 只 Layer 1.5 / 只 priority fix / 無改善 / 反向）

**Power stratification**（避免 sparse cells 誤導）：
- **Powered**：n ≥ 50 → 主分析（rates + Wilson 95% CI）
- **Marginal**：30 ≤ n < 50 → 輔助分析（標記不確定性）
- **Underpowered**：n < 30 → 只列出 cell 存在，不做 rate inference
- 統計：powered / marginal / underpowered cells 各佔幾 %（power curve）

**每 cell 計算**：
1. n, n_TP, n_FP
2. TP rate, FP rate + Wilson 95% CI（避免 n 小時點估計誤導）
3. TP enrichment (vs 全域 baseline), FP enrichment
4. log-odds (TP / FP) + Fisher exact p（vs 全域）
5. powered_flag, marginal_flag, underpowered_flag

**Cell 排序與 top-list（窮舉後選 top）**：
- **Top TP-extreme**：powered cells 中 TP rate 最高的 top 20（候選 TP-marker cells）
- **Top FP-extreme**：powered cells 中 FP rate 最高的 top 20（候選 FP-marker cells）
- **Top log-odds extreme**：兩端各 top 20（最強 TP/FP discrimination）
- **Most populated**：n 最大的 top 20（最 robust 但可能 baseline-like）
- **Surprise cells**：observed rate vs marginal expectation Δ 最大的 top 20（多軸交互效應）

**Confound guard**（對所有 top-list cells 必做）：
1. within-group OLS residualize on `NumReads`（read count confound）
2. within-group OLS residualize on `caller_af`（AF confound）
3. L3 cross-check：在 AF bin 內 cell 是否仍顯著（內生 in-grid，已部分覆蓋）
4. Permutation test：隨機 shuffle labels 後 cell rate 分布的 99 percentile 為閾值
5. Marginal expectation 計算：用 each-axis marginal rate 乘積與 observed rate 比較（找 non-multiplicative cells）

Reference：[L2 collider bias 警告](InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md) + memory `feedback_L2_collider_bias.md` / `feedback_pooled_ols_residualization_trap.md` / `feedback_spatial_autocorrelation_confound.md`

**附加 spatial autocorr guard**（依 memory `feedback_spatial_autocorrelation_confound.md`）：
- 對 top-list cells 額外做 mid-TP-rate chr+pos window 驗證，避免 chr8 hotspot 純空間 artifact

**輸出**：`step2_three_axis_pilot/`（保留資料夾名但內容擴增）
- `step2_grid_5d.tsv`（1000 cells × {n, n_TP, n_FP, rates, Wilson CI, enrichment, log-odds, fisher_p, power flags}）
- `step2_power_curve.png`（n 分布 + powered/marginal/underpowered 比例）
- `step2_top_lists.tsv`（5 個 top-20 lists 合併，含 confound guard 結果欄）
- `step2_marginal_vs_observed.tsv`（multiplicative null vs 真實 rate 差異）
- `step2_facet_grid.png`（高維 grid 投影：固定 2 軸 facet 看其他 3 軸）
- `step2_findings.md`（H2/H3 判定 + powered/marginal cells 統計）

### Step 3：4 個 FP-rich 區域 deep dive

**目的**：對已知的 FP-rich 區域分別評估三軸組合，看哪個區域內三軸能解釋 FP 集中。

**4 個 zone 定義**：

| Zone | 定義 | 來源 | 預期 |
|------|------|------|------|
| Z-OCH | Outer cross_het (LOH.bed outer × HP1+HP2-1) | Thread D | FP rich（germline het 滲漏）|
| Z-CHR8 | HCC1395 chr8 全染色體 | LOH+HPMergedSig hotspot 87.5% 來源 | FP-rich + 樣本特異 |
| Z-GL | CN=3 (Coverage_Multiple 1.8-2.5) ∩ LOH.bed inner | SEQC2 CNV pilot Phase 1 FP rate 12.9% | FP rich (allele imbalance) |
| Z-AUTO | Top 5% FP-density regions（per-region FP/(TP+FP), 用 KDE smoothing） | V6 BAM rerun 後自動偵測 | 跨樣本可泛化的 FP pocket |

**每 zone 做**：
- TP / FP rate vs 全域 baseline
- 在此 zone 內**子 5 軸 grid**（同 step2 1000-cell 結構）：看哪些 cell 在 zone 內仍 powered + 富集
- 5 軸 LR ablation：full 5-feature LR vs leave-one-axis-out（量化每軸增量 AUC）
- HCC1395 chr8 額外做 chr-by-chr 對照（chr8 vs 其他染色體的同 zone 比對）→ 確認 chr8 特異性 vs 機制泛化
- Z-AUTO 與其他 3 個 known zone 的 Jaccard / overlap（H5）

**輸出**：`step3_fp_zone_zoom/`
- `outer_cross_het/` / `chr8_hotspot/` / `gain_loh_zone/` / `auto_detected_top5pct/`
- 每子目錄含：`{zone}_grid.tsv` + `{zone}_lr_ablation.tsv` + `{zone}_findings.md`
- `step3_synthesis.md`（4 zone 對照表 + H3/H4/H5 判定）

### Step 4：HCC1395 有效訊號擴展 4 樣本（已存在資料，後分析）

**前提**：step1-3 在 HCC1395 找到 ≥1 個 robust 信號 cell（after confound guard）。

**現實**：4 樣本 V6 ISM 結果**已存在**於 `phaseD_v6_5sample/{H1437, H2009, HCC1954, HCC1937}/{on,off}_{tp,fp}/`，不需重跑。COLO829 truth set 0600 權限 deferred（雖加上 HCC1395 是 5 樣本，但 V3F/V5 三方數據只有 HCC1395 有 phaseC，其他 4 樣本只有 V6）。

**分析**：
- 對 H1437/H2009/HCC1954/HCC1937 的 V6 ISM 結果，重做 step2 3 軸 grid（注意：只有 V6，沒有 V3F/V5 三方對照）
- Per-cell Wilcoxon signed-rank（n=5 = HCC1395 + 4 樣本），判斷哪些 cell 跨樣本同方向
- HCC1954 outlier 預期會在 Outer cross_het 仍呈異常（Thread D 已知 caller ceiling）—保留作 case study
- HCC1937 marker rate 0.817 邊緣（BRCA1 mutant + 高 ploidy）— 看 grid 是否能解釋此 outlier

**Decision rule**：
- 若 ≥1 個 cell 在 ≥4/5 樣本同方向 + Wilcoxon p<0.05 → 升級為 cross-sample characterization signature（n=5 exact min p=0.0625）
- 否則仍 sample-specific characterization，HCC1395 case study 結論不變

**輸出**：`step4_cross_sample_extension/`
- `step4_per_sample_grid.tsv`（5 樣本 × 50 cells）
- `step4_consistency.tsv`（per-cell direction + Wilcoxon）
- `step4_HCC1937_outlier_analysis.md`（BRCA1 mutant 邊緣 case）
- `step4_signature_candidates.md`

**Follow-up（不在 step4）**：COLO829 truth set 權限解後再補（與 CURRENT_FOCUS Archive TO rerun 連動）

---

## 新資料夾結構

```
research/v6_bam_tpfp_hp_loh_cn/
├── 00_PLAN.md                     ← 本計畫書（從 /bip7_disk/.../plans 複製進來）
├── 01_data_provenance.md          ← V6 BAM / LOH.bed / master TSV / KDE binary commit 來源
├── 02_methodology_notes.md        ← confound guard 協議 + LR ablation 公式
├── step1_baseline_v5_v6_delta/
│   ├── README.md
│   ├── step1_grid_baseline.tsv
│   ├── step1_grid_v5.tsv
│   ├── step1_grid_v6.tsv
│   ├── step1_delta_v5_vs_baseline.tsv
│   ├── step1_delta_v6_vs_v5.tsv
│   ├── step1_delta_v6_vs_baseline.tsv
│   ├── step1_trajectory.tsv
│   ├── step1_three_panel_heatmap.png
│   ├── step1_delta_heatmap.png
│   └── step1_findings.md
├── step2_three_axis_pilot/
│   ├── README.md
│   ├── step2_grid_3d.tsv
│   ├── step2_facets.png
│   ├── step2_confound_guard.tsv
│   └── step2_top_cells.md
├── step3_fp_zone_zoom/
│   ├── README.md
│   ├── outer_cross_het/
│   ├── chr8_hotspot/
│   ├── gain_loh_zone/
│   ├── auto_detected_top5pct/
│   └── step3_synthesis.md
├── step4_cross_sample_extension/
│   ├── README.md
│   ├── step4_per_sample_grid.tsv
│   ├── step4_consistency.tsv
│   └── step4_signature_candidates.md
├── scripts/
│   ├── build_three_axis_grid.py   ← 從 master TSV + LOH.bed 建 50-cell 表
│   ├── confound_guard.py           ← within-group OLS + AF-bin 交叉驗證
│   ├── lr_ablation.py              ← leave-one-axis-out
│   └── auto_detect_fp_zones.py    ← top 5% FP-density KDE smoothing
└── figures/
    ├── step1/, step2/, step3/, step4/
    └── final_synthesis.png
```

---

## Multi-agent 執行階段分工

執行階段（plan approve 後）啟動 4 個 parallel subagent，每個 agent 獨立子目錄避免衝突：

| Agent | 任務 | 輸入 | 輸出 |
|-------|------|------|------|
| **A: V3F/V5/V6 三向整合** | Step 1（讀 phaseC 12 個 ISM 結果整合到 region-level + trajectory + off/on 對照） | phaseC_genome_three_way/ + v3f_vs_v5_vs_v6_*.tsv | `step1_v3f_v5_v6_three_way/` |
| **B: 3 軸 grid + covariate** | Step 2（50-cell grid + power gate + LR with NG/AF/NumReads covariate + collinearity audit + confound guard） | step1 整合 master + LOH.bed + Coverage_Multiple | `step2_three_axis_grid/` |
| **C: FP zone zoom** | Step 3（4 zone deep dive + LR ablation + chr-by-chr 對照） | step2 grid + SEQC2 CNV truth（HCC1395 only） | `step3_fp_zone_zoom/` |
| **D: Auto-detect zone + prior art** | Step 3 Z-AUTO 子任務（KDE FP-density top 5%）+ 背景讀 TumorLens / ROCIT prior art | V6 ISM region-level TP/FP + medRxiv/bioRxiv papers | `step3_fp_zone_zoom/auto_detected_top5pct/` + `02_prior_art_notes.md` |

**Coordinator agent**（主 session）：彙整 4 agent 結果，做 H1-H5 判定 + step3_synthesis.md + 決定 step4 是否啟動。

**Step 4 視 step1-3 結果決定是否啟動**（plan-and-gate），不在初始 fan-out 內。

---

## 關鍵檔案路徑（執行時 reuse）

| 用途 | 路徑 | 備註 |
|------|------|------|
| HCC1395 baseline BAM (upstream longphase) | `/big8_disk/longphase-to/` 產出（路徑 step1 agent 確認） | 原版 longphase-to，未修補 |
| HCC1395 V5 BAM | `/big7_disk/longphase-to-mod/` Layer 1.5 版本產出 | `HaplotagProcess.cpp:512-560` |
| HCC1395 V6 BAM | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/.../tumor_tagged.bam` | 路徑需 step1 agent 重新確認最新 V6 binary 產出 |
| LOH.bed 邏輯 | `InterSubMod/src/core/PhasingGraph.cpp:1817` | VCF VAF ≥0.8 → HOM |
| BAM/HP 讀取 | `InterSubMod/src/io/ReadParser.cpp:123` | HP_Ratio 計算入口 |
| 4-bucket 分類 | `InterSubMod/src/core/LabelTest.cpp:265-302` | Thread D 用同套 |
| Phase 2 BCD code | `LohBedAnnotator.cpp` / `SubcloneAnalyzer.cpp` | 可直接 reuse 4-group 分層 |
| Master TSV | `InterSubMod/output/synthesis/.../master.tsv.gz` | 含 caller_af / NumReads / Coverage_Multiple / Diploid_Coverage_Used |
| SEQC2 CNV truth (HCC1395) | `research/seqc2_cnv_stratification/figures/` 附近 | 6 callers × 21 replicates 共識 |
| Thread D 主報告 | `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md` | 4-bucket × Inner/Outer 框架已建立 |
| Phase BCD 驗證 | `InterSubMod/docs/experiments/validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md` | LOH concordance 98.4% 已驗 |
| 完整觀察整合報告 | `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md` | V5 Layer 1.5 缺陷 §8.6 |

---

## 驗證方式（end-to-end）

1. **Step 1 通過條件**：HCC1395 V6 BAM 跑出 inner/outer × 4-bucket × TP/FP 表，與 Thread D X5 報告（chr19 同範圍）數值偏差 < 5%。
2. **Step 2 通過條件**：3D grid 50 cells 全跑完，confound guard 對 ≥1 cell 顯示 residualized rate vs raw rate 差異 < 0.05（即非 confound 驅動）。
3. **Step 3 通過條件**：4 zone 各自有 LR ablation 結果，H3/H4/H5 任一獲得明確判定（POSITIVE / NEGATIVE / MIXED）。
4. **整體 deliverable**：HCC1395 pilot 報告 → 寫入 `InterSubMod/docs/experiments/in_progress/2026/05/`（依 doc-standards skill 命名）+ 同步 `docs/experiments/INDEX.md`。

**驗證命令模板**：
```bash
# Step 1 跑完後
python research/v6_bam_tpfp_hp_loh_cn/scripts/build_three_axis_grid.py --bam-version v6 --sample HCC1395 --output step1_grid_v6.tsv

# Step 2 confound guard
python research/v6_bam_tpfp_hp_loh_cn/scripts/confound_guard.py --input step2_grid_3d.tsv --residualize caller_af,NumReads

# Step 3 LR ablation per zone
python research/v6_bam_tpfp_hp_loh_cn/scripts/lr_ablation.py --zone chr8_hotspot --features LOH,HP,CN
```

---

## Out-of-scope（明確排除避免 scope creep）

- ❌ **不評估 filter / ΔF1**（用戶明示 characterization-only；Plan agent §5 警告 LR ablation 不出 AUC 改報 deviance decomposition）
- ❌ 不修改 C++ pipeline（純 post-hoc TSV 分析）
- ❌ 不重跑 V6 BAM 或 ISM（phaseC 三向 + phaseD 4 樣本已存在）
- ❌ 不做 Wakhan / SAVANA external CN（KDE-only 即可，prior art 背景讀做差異化說明）
- ❌ 不做 Phase 2A Normal Methylation Reference（另一條線）
- ❌ 不重新 phasing（V6 重用 V5 phased VCF，caller F1 不變已驗證）
- ❌ 不加 contrast 樣本同步跑 step2（用戶明示保持 HCC1395 single-pilot；接受跨樣本崩風險到 step4 才驗）

---

## 預期時程（v0.3 大幅縮短，因 BAM/ISM 不需重跑）

| Step | 估時 | 依賴 |
|------|------|------|
| Step 1 | 0.5 天 | phaseC_genome_three_way master.tsv 整合 |
| Step 2 | 1 天 | Step 1 完成 + power dry-run gate 通過 |
| Step 3 | 1.5 天 | Step 2 完成 + SEQC2 CNV truth |
| Step 4 | 1 天（HCC1395 + 4 樣本 V6 ISM 已存在） | Step 2 完成 |
| HCC1395 pilot 總計 | **~3 天**（含 step4） | — |

**Gate 1（power dry-run，step1 後）**：marginal 計算 powered cells <15 → 降至 2 軸（LOH × HP）10-cell grid 或合併 cov bin
**Gate 2（step3 後 mid-review）**：根據 H1a/H1b/H1c/H2-H7 判定，決定 step4 跨樣本擴展 / 終止 / 修正

---

## 開放問題（執行階段需確認）

1. ~~V6 HCC1395 BAM 路徑~~ ✅ 已確認 `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam`
2. ~~是否需重跑 BAM/ISM~~ ✅ 不需，phaseC 三向 + phaseD 4 樣本已存在
3. phaseC master.tsv 內 Coverage_Multiple 是否用新 KDE binary（agent A 確認 Diploid_Coverage_Used 欄位 ≠ 75.0）
4. Z-AUTO top 5% 的 region 大小 / smoothing kernel（agent D 與 coordinator 討論）
5. Power threshold（預設 n≥50 powered / 30 marginal / <30 underpowered）是否需依 grid sparsity 調整
6. 是否要加 chr 軸（chr 23 × 50 = 1150 cells，仍 sparse；預設不加，但 zone Z-CHR8 已涵蓋 chr8）
7. NGroups bin 是否合併 NG=0+NG=1（兩者都接近 LOH 區域；預設保留分開作 covariate）
8. COLO829 truth set 0600 權限解後是否補入 step4（與 CURRENT_FOCUS Archive TO rerun 連動）
