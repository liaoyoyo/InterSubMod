---
title: LOH × AF × KDE TP/FP Discrimination 研究專案索引
status: in_progress
owner: liaoyoyo2001
created: 2026-04-22
last_updated: 2026-04-22
parent_report: /big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/radiant-moseying-turtle.md
---

# LOH × AF × KDE TP/FP Discrimination 研究專案

本資料夾統合「在新 KDE baseline 下，透過 LOH_Subtype × AF_class × CN-tier × NG × NumReads 五維參數切分 somatic caller 輸出 TP/FP 差異區域」的完整研究過程與證據。

---

## 📌 2026-04-22 重要更新（ISM rerun 完成）

**14:40 完成**：5 個 archive TO pilot（HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437）重跑 post-HP-fix ISM binary，補齊 LOH_Subtype / Diploid_Coverage_Used / Phase 2 features 欄位。

**14:49 完成**：`obs09-12` 全面重跑，涵蓋 **13 sample-mode（6 TO + 7 paired_full）**，原限於 8 sample-mode 的分層分析全面升級。

**Hard Gate 全面觸發**（Plan §6.3）：`n_samples_high ≥ 5` 的 cells 達 **16 個**（閾值 10，**超標 60%**），`n_samples_high ≥ 3` 達 64 個（閾值 50）。**LOSO 正式驗證強烈建議啟動**。

**核心結論逆轉**：原「LOH/NG 梯度為 HCC1395_TO-specific」結論，在新 6 TO 樣本對比下，**LOH 與 NG 梯度跨 5-6/6 TO 樣本一致**（HCC1954_TO AF 方向結構性例外）。H4（跨樣本泛化性）判決從「弱支持」升級為「中-強支持」。

**權威分析**：`07_figure_layers.md`（rerun 後權威版本 V2，**37 張 research-question 導向 figures**）；**08 已 superseded**。

---

## 🗺 依時間選讀路徑

| 時間 | 閱讀內容 | 能回答的問題 |
|:-----|---------|-------------|
| **30 秒** | 本頁 📌 更新區 + §關鍵結論速查 | H4 判決？共識 cells 數？下一步？ |
| **5 分鐘** | 上述 + `07_figure_layers.md §7.0 TL;DR` + 看 `figures/new/L3/L3_loh_af_cn.png` | 上述 + 哪張圖最代表結論？ |
| **30 分鐘** | 上述 + `07 §7.5 consistency (Tier A/B/C)` + `07 §7.7 Hard Gate` + `data/cell_consistency_matrix.tsv` 前 16 row | 上述 + 16 cells 的具體 pattern？為何 LOSO 該做？ |
| **完整** | 全部 01-08（含 37 張圖 + 04 §4.9 H4 升級）| 全部研究脈絡 |

---

## 研究問題

在 ClairS/ClairS-TO 輸出的 VCF 上，依照我們定義的 5 個 biology-inspired 維度切分，**能否找出 TP:FP 比值顯著優於 caller baseline 的 operating points？是否存在超越既定 scheme 的 Pareto-optimal 組合？**

**關鍵測試場**：HCC1395 TO mode（baseline TP=71.1%，11,606 FP 最高統計力）+ COLO829 paired_full（baseline TP=94.1%，2,273 FP）。

---

## 研究階段與文件對應

| Step | 研究問題 | 腳本 | 主輸出 |
|------|---------|------|--------|
| Step 1 | 各 sample × mode baseline TP:FP | `step4_tpfp_discrimination.py` | `tpfp_stratified_filter_schemes.tsv`（paired_full） |
| Step 2 | TO mode TP/FP 區域切分 | `step4b_tpfp_tomode.py` | `tpfp_cube_coarse_TO.tsv`、`tpfp_stratified_filter_schemes_TO.tsv` |
| Step 3 | S1-S6 基礎圖表 | `step5_visualize_tpfp.py` | fig_v2_1 ~ fig_v2_6 |
| Step 4 | baseline vs scheme fold-improvement + 統計力 + 特徵差 | `step6_tpfp_detailed.py` | `tpfp_baseline_ratio.tsv`、`tpfp_per_sample_scheme_full.tsv`、`tpfp_feature_diffs.tsv`、`tpfp_subschemes.tsv` |
| Step 5 | F7-F10 深度圖 | `step7_visualize_detailed.py` | fig_v2_7 ~ fig_v2_10 |
| Step 6 | 5D cube Pareto envelope | `step8_multidim_panorama.py` | `tpfp_5d_cube_*.tsv`、fig_v2_11a/b |
| Step 7 | 多角度觀察圖（本資料夾新增）| `obs01-obs08_*.py` | `figures/new/*.png` |

---

## 文件導覽

| 檔案 | 用途 | 優先閱讀順序 |
|------|------|-------------|
| `00_INDEX.md`（本檔）| 專案索引 + 導覽 | 1 |
| `01_research_question.md` | 研究問題 + 假說定義 | 2 |
| `02_methodology.md` | 五維參數定義 + scheme 設計 + 統計指標 | 3 |
| `03_data_inventory.md` | 資料來源 + 限制 + 每 sample 基本量 | 4 |
| `04_comparison_narrative.md` | **核心敘述比較分析**（主要成果） | 5 ⭐ |
| `05_biology_interpretation.md` | 生物學機制詮釋（為何每個 cell 富集 TP/FP） | 6 |
| `06_limitations_future_work.md` | 限制 + 泛化性質疑 + 未來工作 | 7 |
| `07_figure_layers.md` | **分層 combinatorial 圖集**（L1→L2→L3→Consistency × **13 sample-mode** 跨樣本驗證，**2026-04-22 ISM rerun 後權威版**）| 8 ⭐ |
| `08_archive_to_crosssample.md` | Archive TO pre-rerun 彙整紀錄（**已 superseded by 07**；LOH/CN 已補齊）| 9 |
| `09_TO_sample_af_lohside_ng.md` | **TO 模式 4D 分析**（sample × AF × LOH Inner/Outer × NG）：NG=2 Inner filter 新候選 + HCC1954 anomaly 偵測 | 10 ⭐ |

---

## 圖表索引（精簡版）

### 🌟 快速 Gallery（5 張代表圖）

若只能看 5 張圖：

![fig_v2_1 既有 HCC1395 TO](figures/existing/fig_v2_1_to_tp_heatmap.png)
*`fig_v2_1`（既有）：HCC1395 TO 單樣本 LOH×AF×CN 熱圖 —— 原始研究起點*

![fig_v2_7 post-rerun 6 TO panels](figures/existing/fig_v2_7_baseline_tp_fp_fold.png)
*`fig_v2_7`（2026-04-22 rerun 後）：6 個 TO 樣本 baseline vs scheme fold-improvement bar （完整路徑亦在 `docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled/`）*

![L3_loh_af_cn ⭐ 核心圖](figures/new/L3/L3_loh_af_cn.png)
*`L3_loh_af_cn`（⭐ 核心圖）：fig_v2_1 的 13 sample-mode 擴展版；per-sample CN slice 跨 TO/paired 對比*

![L1_loh_TO](figures/new/L1/L1_loh_TO.png)
*`L1_loh_TO`：LOH_Subtype 在 6 TO 樣本的 TP rate 梯度 —— 首次證實跨樣本一致（非 HCC1395-only）*

![consistency 5D Top-30](figures/new/consistency/cell_consistency_5d.png)
*`cell_consistency_5d`：16 cells n_samples_high ≥ 5/13 的 Top-30 stacked bar（Hard Gate 觸發核心證據）*

**完整 37 張 L1/L2/L3/consistency + 12 張既有 fig_v2_* + 8 張 obs01-08 圖** → 見附錄 A。
**詳細分層說明與每張解讀** → `07_figure_layers.md`。

### Archive TO 跨樣本觀察（2026-04-22 ~03:00 建立，已 superseded）

**詳細說明**：見 `08_archive_to_crosssample.md`（標記 superseded）。當時的 4 張 archive_to 圖保留作 pre-rerun vs post-rerun 對照；核心結論已被 `figures/new/L1/` 新圖取代。詳見 07 §7.6 升級對照表。

---

## 關鍵結論速查（v2 §17 證實）

| Operating Point | Purity | Recall | TP:FP Ratio | Fold vs Baseline (2.47:1) |
|-----------------|-------:|-------:|------------:|--------------------------:|
| Caller baseline | 71.1% | 100% | 2.47 | 1.0× |
| S3 Diploid Het | 95.5% | 1.3% | 21.35 | 8.69× |
| S5 Combo | 91.8% | 2.9% | 11.14 | 4.53× |
| **17-cell Pareto envelope** | **96.1%** | **3.7%** | **24.56** | **10.0×** ⭐ |
| **28-cell Pareto envelope** | **90.4%** | **7.4%** | **11.68** | **4.73×** |

**重要修正**：§14 的「單特徵 refinement 飽和」僅限於既定 scheme 內；聯合 5D cube **並未飽和**。

---

## 相關資源

- **母報告**：`docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md`（v2，682 行）
- **母計劃**：`.claude/plans/radiant-moseying-turtle.md`
- **v1 偏題報告**：`docs/experiments/in_progress/2026/04/20260421_NG_KDE_Rescaled_Multi_CN_Analysis_01.md`（已標 superseded）
- **資料 master**：`data/master.tsv.gz` → `/big7_disk/.../ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz`

---

## 更新歷史

- **2026-04-22 15:31** — V2 圖表重繪 + 敘述重構：37 張新圖（L1/L2 拆 TO/paired + L3 caption + consistency TO-only）；每張圖 research-question suptitle + takeaway caption box + panel metadata badge；07/00_INDEX/08/04 §4.9 全面更新
- **2026-04-22 14:49** — obs09-12 首輪重跑完成；涵蓋 13 sample-mode；Hard Gate 全面觸發
- **2026-04-22 14:40** — obs16 下游重生成（step0→step7）完成；新 master 16.9 MB 698k rows；fig_v2_7_baseline_tp_fp_fold.png 改為 6 TO panel 版
- **2026-04-22 14:30** — obs15/obs15b 完成：5 archive TO pilot 的 post-HP-fix ISM rerun（共 328,979 regions 處理完成）
- **2026-04-22 ~05:00** — 本資料夾建立，連結 v2 既有 16 圖 + 16 TSV + 6 scripts；新增 obs01-08（每角度觀察）+ obs09-12（L1-L3 + consistency）+ obs13-14（pre-rerun archive 彙整）

---

## 附錄 A：完整圖表 Gallery

### A.1 既有圖（ng_kde_rescaling/ 產出）

| 檔名 | 主題 |
|------|------|
| `figures/existing/fig_v2_1_to_tp_heatmap.png` | TO TP rate heatmap（LOH×AF×CN）|
| `figures/existing/fig_v2_2_to_fp_heatmap.png` | TO FP rate heatmap |
| `figures/existing/fig_v2_3_filter_scheme_bar.png` | Scheme S1-S8 bar |
| `figures/existing/fig_v2_4_per_sample_schemes.png` | paired_full × S1-S6 |
| `figures/existing/fig_v2_5_operating_points.png` | Scheme operating points |
| `figures/existing/fig_v2_6_biology_module_interpretation.png` | 生物模組詮釋 heatmap |
| `figures/existing/fig_v2_7_baseline_tp_fp_fold.png` | baseline vs scheme fold-improvement（**rerun 後 6 TO panel**）|
| `figures/existing/fig_v2_8_per_sample_scheme_heatmap.png` | Per-sample × scheme TP rate |
| `figures/existing/fig_v2_9_feature_violin_in_S4.png` | S4 內 TP/FP 特徵 violin |
| `figures/existing/fig_v2_10_subscheme_operating_points.png` | Sub-scheme operating points |
| `figures/existing/fig_v2_11_panorama_HCC1395_TO.png` | 5D cube + Pareto envelope |
| `figures/existing/fig_v2_11b_panorama_COLO829.png` | COLO829 panorama |

### A.2 obs01-08（pre-rerun 觀察）

| 檔名 | 主題 |
|------|------|
| `figures/new/obs01_af_distribution_per_sample.png` | AF 分佈（TP vs FP 疊加）|
| `figures/new/obs02_cn_distribution_per_sample.png` | CovM 分佈（TP vs FP）|
| `figures/new/obs03_feature_pairplot.png` | 多特徵 pairwise 散點 |
| `figures/new/obs04_chr_spatial_scheme.png` | Scheme 在 chr 上的空間分佈 |
| `figures/new/obs05_dim_contribution.png` | 單維度邊際貢獻 |
| `figures/new/obs06_ng_by_sample_stacked.png` | NG 分佈堆疊 |
| `figures/new/obs07_pareto_trajectory.png` | 17-cell → 32-cell 軌跡 |
| `figures/new/obs08_overlap_venn_schemes.png` | S3 / Top-17 / Top-28 Venn |

### A.3 L1 單維度（V2 2026-04-22，10 張：5 dim × 2 mode）

| TO (6 sample 2×3) | paired (7 sample 3×3) |
|------|------|
| `figures/new/L1/L1_loh_TO.png` | `figures/new/L1/L1_loh_paired.png` |
| `figures/new/L1/L1_af_TO.png` | `figures/new/L1/L1_af_paired.png` |
| `figures/new/L1/L1_cn_TO.png` | `figures/new/L1/L1_cn_paired.png` |
| `figures/new/L1/L1_ng_TO.png` | `figures/new/L1/L1_ng_paired.png` |
| `figures/new/L1/L1_nr_TO.png` | `figures/new/L1/L1_nr_paired.png` |

### A.4 L2 雙維度（V2 2026-04-22，20 張：10 pair × 2 mode）

10 pairs: LOH×AF / LOH×CN / LOH×NG / LOH×NR / AF×CN / AF×NG / AF×NR / CN×NG / CN×NR / NG×NR
每 pair 兩張：`figures/new/L2/L2_{x}_x_{y}_TO.png` + `{...}_paired.png`

### A.5 L3 三維度（4 張，13 sample-mode 合併）

| 檔名 | 主題 |
|------|------|
| `figures/new/L3/L3_loh_af_cn.png` | ⭐ **fig_v2_1 × 13 sample 擴展版** |
| `figures/new/L3/L3_loh_af_ng.png` | LOH × AF × NG |
| `figures/new/L3/L3_af_cn_nr.png` | AF × CN × NR |
| `figures/new/L3/L3_loh_cn_ng.png` | LOH × CN × NG |

### A.6 Consistency（3 張）

| 檔名 | 主題 |
|------|------|
| `figures/new/consistency/cell_consistency_5d.png` | Top-30 5D cells (13 sample)|
| `figures/new/consistency/cell_consistency_3d_locn.png` | LOH×AF×CN 空間一致性 scatter |
| `figures/new/consistency/cell_consistency_TO_only.png` | **TO-only 版 Top-30**（6 sample，V2 新增）|

### A.7 Archive TO（pre-rerun 歷史，已 superseded by A.3-A.6）

`figures/new/archive_to/fp_scale_comparison.png` / `L1_af_class_to.png` / `L1_hpfinengroups_to.png` / `L2_af_x_ng_to.png`
