<!--
建立時間: 2026-04-02 15:00
目標: LOH 判定機制研究與觀察的集中管理區塊
處理範圍: ISM LOH 判定邏輯、LongPhase-TO LOH.bed、self-phasing 因果鏈、read threshold 問題
關聯檔案:
  - src/core/RegionProcessor.cpp (is_potential_loh, compute_hp_ratio)
  - src/core/SignificanceAnalyzer.cpp (potential_loh in verification_class)
  - docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md
  - docs/reports/validated/2026/04/20260402_loh_read_threshold_visual_argument_01.md
-->

# LOH Investigation — 研究區塊

## 目的

集中管理 ISM 中 LOH 判定相關的研究、觀察、數據與圖表。涵蓋：

1. **ISM LOH 判定邏輯** — `is_potential_loh()` 的 HP_Ratio 閾值、無最低 read 門檻問題
2. **LongPhase-TO LOH.bed** — 區域級 LOH 判定、與 ISM 的差異
3. **Self-phasing 因果鏈** — TO LOH 過量的根本原因
4. **跨模式比較** — TO vs Paired 的 LOH 行為差異
5. **改善方案評估** — read threshold fix、scaffold 改善、re-phasing 策略

---

## 目錄結構

```
research/loh_investigation/
├── README.md              # 本文件
├── data/                  # 分析用中間數據（TSV, CSV）
├── figures/               # 關鍵圖表（從 workspace 複製的精選版）
├── reports/               # 研究筆記與分析報告
└── scripts/               # 專用分析腳本（輕量級，非 pipeline 級）
```

---

## 已確認的關鍵結論

### ⚠️ SEQC2 驗證修正（2026-04-02）

**先前假設已被修正**：LOH.bed 在體染色體上高度準確。

| 假設 | 修正前 | 修正後（SEQC2 驗證） |
|------|--------|---------------------|
| LOH.bed 受 self-phasing 嚴重汙染 | 認為過量 | **體染色體 F1=96.2%, Jaccard=0.928** |
| TO LOH 過量來自 self-phasing | 全歸因 self-phasing | **73% TO-only 來自 chrX 半合子** |
| HCC1395 LOH 覆蓋率異常 | 認為不正常 | **SEQC2 確認 49.2% genome 是 LOH** |

### LOH.bed 準確性（SEQC2 Benchmark）

| 指標 | 含 chrX | 體染色體 only |
|------|---------|-------------|
| Sensitivity | 96.1% | 96.1% |
| Precision | 87.7% | **96.4%** |
| F1 Score | 91.7% | **96.2%** |
| Jaccard | 84.7% | **92.8%** |

**四類分類（bp-level）**：TP=1,432 Mb | FP=200 Mb（73% chrX）| FN=58 Mb | TN=1,340 Mb

### 三套 LOH 判定系統

| # | 系統 | 定義 | TO rate | Paired rate |
|---|------|------|---------|-------------|
| 1 | ISM `Potential_LOH` | HP_Ratio < 0.1 or > 0.9 | 41.8% | 29.4% |
| 2 | 分析 `core_loh_like` | 同上（無 smoothing） | 41.8% | 29.4% |
| 3 | LongPhase-TO `LOH.bed` | Phased genotype ratio | 26.7% | N/A（無輸出） |

- ISM 兩套 100% 一致（epsilon smoothing 不影響 0.1/0.9 閾值判定）
- ISM 比 LOH.bed 多判 63,610 cases（+15.2%）— **site-level HP noise，非 LOH.bed 問題**
- Paired 模式完全沒有 LOH.bed — LOH 判定 100% 來自 ISM

### LOH 無法區分 TP/FP

| 分類器 | TO AUC | Paired AUC |
|--------|--------|------------|
| LOH.bed hit | 0.542 | N/A |
| ISM Potential_LOH | 0.544 | 0.472（反向） |

原因：LOH 區域內的 FP 是 germline variants（avg AF=0.964），是真實 variant 非 artifact。

### Self-Phasing 影響重新定位

- **LOH.bed region-level**：體染色體準確率 >96%，self-phasing 影響極小
- **ISM site-level HP_Ratio**：self-phasing 仍是主因（ehp≥200 時 TO LOH 仍為 Paired 的 14 倍）
- **結論**：self-phasing 主要影響 ISM 的 HP_Ratio 計算，而非 LOH.bed 的 region-level 判定

### PON-Only Phasing 驗證（2026-04-03）

| 指標 | Baseline | PON-only | 判定 |
|------|----------|----------|------|
| LOH.bed Jaccard (兩者間) | — | **1.0000** | 完全一致 |
| Somatic phasing bias (HP1:HP2) | **17.3:1** | 消除 | Self-phasing 確認 |
| Phase block N50 (variants) | 4,061 | **8,109** | +99.7% |
| Phased variant rate | 54.9% | **78.5%** | +23.6 pp |
| 執行時間 | 2,693s | **1,976s** | 1.36x 快 |

**Self-phasing 直接證據**：baseline 的 614,471 somatic variants 中 94.6% 被分配到 HP1（`1|.`），僅 5.4% 到 HP2（`.|1`）。PON-only 修正後 somatic variants 改為 `0|.`/`.|0`（無 HP 分配）。

### 雙問題分解（ISM site-level 限定）

| 問題 | 影響 | 修復 |
|------|------|------|
| A: 無 read 門檻 | ~3% LOH | 簡單（加 ehp≥10） |
| B: Self-phasing scaffold | ~97% 過量 LOH（ISM HP_Ratio level） | **需 haplotag + ISM 重跑驗證** |

### QS v2 Zone-Aware 驗證（2026-04-02）

| 公式 | TO inside AUC | TO overall AUC | 備註 |
|------|--------------|----------------|------|
| QS v1 (current) | 0.578 | 0.546 | baseline |
| QS v2a (caller features) | 0.590 (+0.012) | 0.554 (+0.008) | 需新增 VCF 欄位 |
| QS v2b (ISM only) | 0.609 (+0.031) | 0.557 (+0.012) | 不需改 C++ VCF |

- GO threshold: TO overall ≥ +0.05, TO LOH-inside ≥ +0.10 → **均 FAIL**
- Cross-sample: 7/7 improved → PASS
- 結論：CONDITIONAL GO — 改善真實但不足，rule-based QS 有天花板

### Read-Level methyl_mean LOH 驗證（2026-04-02）

| 模式 | LOH Zone | methyl_mean AUC | Cohen's d | TP median | FP median |
|------|----------|----------------|-----------|-----------|-----------|
| **TO** | **inside** | **0.440** | **-0.12** | 0.771 | 0.845 |
| TO | outside | 0.374 | -0.45 | 0.572 | 0.697 |
| Paired | outside | 0.219 | -1.08 | 0.544 | 0.830 |

- **先前 Paired inside LOH AUC=0.785 不可複現**：僅 9 TP regions，HCC1395-only 有零個 paired TP
- **TO methyl_mean 在 LOH 內完全無法區分 TP/FP**（6 特徵全部 AUC < 0.56）
- Per-sample AUC = 0.61-0.65（Simpson's paradox），但仍低於 G1-G7 門檻 0.64
- **結論**：methyl_mean 不可用於 TO LOH 區域 TP/FP 區分，甲基化方向全面關閉

---

## 相關外部資源

| 資源 | 位置 |
|------|------|
| Master dataset | `big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz` |
| Causal chain workspace | `big7_disk_output/synthesis/observation_workspaces/20260402_phasing_causal_chain/` |
| Read threshold workspace | `big7_disk_output/synthesis/observation_workspaces/20260402_loh_read_threshold_analysis/` |
| LOH.bed files | `big7_disk_output/synthesis/research_rounds/{round_id}/step03_longphase_to/tumor_phased_LOH.bed` |
| ISM 原始碼 | `src/core/RegionProcessor.cpp:35-49` (compute_hp_ratio, is_potential_loh) |
| 因果鏈報告 | `docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md` |
| Read threshold 報告 | `docs/reports/validated/2026/04/20260402_loh_read_threshold_visual_argument_01.md` |
| 文獻調查 | `docs/references/manual/20260402_longphase_to_phasing_quality_literature.md` |
| SEQC2 CNV Benchmark | `research/loh_investigation/data/seqc2_cnv_benchmark_v4/` |
| SEQC2 vs TO 驗證數據 | `research/loh_investigation/data/seqc2_vs_to_validation/` |
| 四類分類 BED | `research/loh_investigation/data/seqc2_vs_to_validation/loh_classified_{TP,FP,FN}.bed` |
| SEQC2 論文整理 | `research/loh_investigation/reports/20260402_seqc2_cnv_benchmark_analysis.md` |
| SEQC2 交集驗證 | `research/loh_investigation/reports/20260402_seqc2_vs_longphase_to_loh_validation.md` |
| 四類分類報告 | `research/loh_investigation/reports/20260402_loh_four_class_classification.md` |
| O14 Haplotag Bias 腳本 | `research/loh_investigation/scripts/observe_loh_haplotag_bias.py` |
| O15 Phase 1 腳本 | `scripts/analysis/build_observation_O15_loh_zone_metrics_hcc1395.py` |
| O15 Phase 2 腳本 | `scripts/analysis/build_observation_O15b_loh_zone_metrics_cross_sample.py` |
| O15 Phase 1 報告 | `research/loh_investigation/reports/20260402_O15_phase1_loh_zone_metrics_hcc1395.md` |
| O15 Phase 2 報告 | `research/loh_investigation/reports/20260402_O15_phase2_loh_zone_metrics_cross_sample.md` |
| O15 圖表 (32 張) | `research/loh_investigation/figures/o15_p{1,2}_fig*.png` |
| O15 數據 (6 TSV) | `research/loh_investigation/data/o15_p{1,2}_*.tsv` |
| GE.bed 分析腳本 | `research/loh_investigation/scripts/observe_to_ge_vs_seqc2.py` |
| QS v2 驗證腳本 | `scripts/analysis/build_qs_v2_zone_aware_validation.py` |
| QS v2 驗證報告 | `research/loh_investigation/reports/20260402_qs_v2_zone_aware_validation.md` |
| QS v2 數據 | `research/loh_investigation/data/qs_v2_*.tsv` |
| TO Read-Level 甲基化腳本 | `research/loh_investigation/scripts/analyze_to_read_methyl_in_loh.py` |
| TO Read-Level 甲基化報告 | `research/loh_investigation/reports/20260402_to_read_methyl_in_loh.md` |
| TO Read-Level 數據 | `research/loh_investigation/data/to_read_methyl_stats.tsv` |
| TO Read-Level 圖表 (12 張) | `research/loh_investigation/figures/to_read_methyl_fig*.png` |
| LOH methyl level 腳本 | `research/loh_investigation/scripts/analyze_loh_methyl_level.py` |
| LOH methyl level 數據 | `research/loh_investigation/data/loh_methyl_level_stats.tsv` |
| 三方比較腳本 | `research/loh_investigation/scripts/analyze_hp_ratio_vs_loh_bed_vs_seqc2.py` |
| 三方比較報告 | `research/loh_investigation/reports/20260403_hp_ratio_vs_loh_bed_vs_seqc2_threeway_comparison.md` |
| 三方比較數據 (7 TSV) | `research/loh_investigation/data/loh3way_*.tsv` |
| 三方比較圖表 (9 張) | `research/loh_investigation/figures/loh3way_fig01-09*.png` |
| PON-only phasing 驗證報告 | `research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md` |
| PON-only phasing 修改碼 | `/big7_disk/liaoyoyo2001/longphase-to-mod/` (PhasingProcess.cpp, Phasing.cpp, PhasingProcess.h) |
| PON-only baseline 輸出 | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/` |
| PON-only pononly 輸出 | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly/` |
| HCC1395 測試數據 (big7) | `/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/` (BAM 272G + VCF) |
| PON/REF 測試數據 (big7) | `/big7_disk/liaoyoyo2001/data/PON/`, `/big7_disk/liaoyoyo2001/data/ref/` |

---

## 研究時間線

| 日期 | 里程碑 |
|------|--------|
| 2026-03-27 | LOH Round 1 Cross-Sample Audit — 發現 TO LOH 過量 |
| 2026-03-28 | LOH Round 2-4 — HP fix 重跑、enrichment 分析 |
| 2026-03-30 | O1-O13 系統性觀察 — 確認單一特徵 AUC < 0.60 |
| 2026-04-01 | G1-G7 Germline FP 鑑別 — NO-GO（60+ 特徵全 < 0.64） |
| 2026-04-02 | Self-phasing 因果鏈確認 — 五步因果鏈 + 23 TSV + 7 PNG |
| 2026-04-02 | Read threshold 視覺論證 — 6 圖確認雙問題分解 |
| 2026-04-02 | 三套 LOH 系統比對 — ISM 100% 自洽、LOH.bed AUC 0.54 |
| 2026-04-02 | **SEQC2 驗證** — LOH.bed 體染色體 F1=96.2%，修正先前 self-phasing 假設 |
| 2026-04-02 | 四類分類 — TP/FP/FN/TN BED 切割 + 4 張視覺化 + chrX 主因確認 |
| 2026-04-02 | O14 Haplotag Bias — 99.5% extreme HP_Ratio in TP LOH, 100% region 一致性 |
| 2026-04-02 | GE.bed/CNV 分析 — LGE 覆蓋 97.8% genome，LongPhase-TO 無 CNV calling 能力 |
| 2026-04-02 | **O15 Phase 1** — HCC1395 + SEQC2 4-class，16 張圖，LOH 內甲基化 AUC ~0.50 |
| 2026-04-02 | **O15 Phase 2** — 7 samples × LOH.bed binary，16 張圖，跨樣本一致失效確認 |
| 2026-04-02 | **QS v2 Zone-Aware** — 3 公式 × 8 圖，TO LOH-inside AUC +0.031，CONDITIONAL GO |
| 2026-04-02 | **methyl_mean LOH 分析** — Paired inside AUC=0.785（9 TP regions，不可靠） |
| 2026-04-02 | **TO Read-Level methyl_mean** — TO inside LOH AUC=0.440（NEGATIVE），先前 paired 結果不可複現 |
| 2026-04-03 | **三方 LOH 比較** — ISM ⊃ LOH.bed（recall 99%+），excess 10-16% 來自低 ehp 噪音；TP/FP rate 無差異；9 圖 7 TSV |
| 2026-04-03 | **PON-only phasing 驗證** — LOH.bed 完全一致（Jaccard=1.0）；somatic phasing bias 17.3:1 消除；N50 翻倍；phased rate +23.6pp；self-phasing 不影響 region-level LOH |
