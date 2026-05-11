---
title: 資料清冊 — 來源、欄位、樣本計量
status: in_progress
last_updated: 2026-04-22
---

# 03. 資料清冊

## 3.1 Master 資料集

**主資料檔**：`data/master.tsv.gz`
  - symlink → `/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz`
  - 建立於：2026-04-22（NG × KDE rescaling 研究）
  - 格式：TSV gzip，包含所有 TP (intersubmod_tp/) + FP (intersubmod_fp/) 的 significance_summary.csv rows

**總筆數**：368,997 variants

### 3.1.1 欄位清單（28 欄）

| 欄位 | 型別 | 來源 | 用途 |
|------|-----|------|-----|
| RegionID | int | caller | 唯一識別 |
| Chr, Pos | str/int | VCF | 定位 |
| NumReads | int | ReadParser | 單 variant 覆蓋讀數 |
| NumCpGs | int | ReadParser | ISM 特徵基礎 |
| HPFineNGroups | int | C++ ISM | **D4 分箱** (flag=off 值) |
| Coverage_Multiple | float | KDE | raw diploid coverage ratio |
| Coverage_Category | str | KDE | 分類標籤 |
| AlleleDelta | float | ReadParser | TP/FP 特徵比較重要 |
| HPFine_NGroups_CF | float | C++ ISM | NG 信心度 |
| Diploid_Coverage_Used | float | KDE | 2026-04 新 baseline |
| Potential_LOH | bool | ReadParser | fallback LOH 來源 |
| NTumorReads | float | ReadParser | Tumor 讀數子集 |
| NNormalReads | float | ReadParser | Normal 讀數子集（僅 paired mode 有） |
| LOH_Bed_Overlap | bool | PON-only LOH bed | 主 LOH 來源 |
| LOH_Bed_Annotation | str | PON-only LOH bed | 詳細 LOH 類別 |
| **LOH_Subtype** | str | 衍生 | **D1 分箱** |
| tp_label | int | SEQC2 truth | 0=FP / 1=TP |
| **sample**, **mode** | str | pipeline | 樣本分組 |
| kde_status | str | 標籤 | new_direct / stale_rescaled / phase1_new |
| **CovM_used** | float | 衍生 | **D3 輸入** — post-fix KDE 標準化後的 coverage multiple |
| baseline_x | float | KDE | 樣本 baseline X 值 |
| sample_order | int | 排序鍵 | |
| loh_side | str | 衍生 | Inner / Outer / Unknown |
| **AF** | float | VCF | 原始 allele frequency |
| AF_bin10 | str | 衍生 | 10-bin 細切 |
| **AF_class** | str | 衍生 | **D2 分箱** |

粗體欄位為本研究 5D cube 的直接輸入。

---

## 3.2 樣本 × Mode 矩陣（baseline TP:FP）

取自 `data/tpfp_baseline_ratio.tsv`：

| Sample | Mode | n_total | n_TP | n_FP | TP:FP ratio | baseline TP% | Power |
|--------|------|--------:|-----:|----:|----------:|-------------:|-------|
| **HCC1395** | paired_full | 30,381 | 29,754 | 627 | 47.5:1 | 97.9% | 🟡 中 |
| **HCC1395** | to_pileup | 40,115 | 28,509 | **11,606** | **2.47:1** | **71.1%** | ⭐ 極高（主測試場）|
| HCC1395_DORADO | paired_full | 30,129 | 29,889 | 240 | 124.5:1 | 99.2% | 🟡 低 |
| HCC1937 | paired_full | 12,588 | 12,393 | 195 | 63.5:1 | 98.5% | 🟡 中 |
| HCC1954 | paired_full (stale rescaled) | 17,938 | 17,909 | 29 | 617.6:1 | 99.8% | 🔴 極低 |
| H2009 | paired_full | 132,995 | 132,909 | 86 | 1545.5:1 | 99.9% | 🔴 極低 |
| H1437 | paired_full | 67,476 | 67,468 | 8 | 8433.5:1 | 99.99% | 🔴 不可用 |
| **COLO829** | paired_full | 37,458 | 35,185 | **2,273** | **15.5:1** | 93.9% | ⭐ 高（次要）|

**關鍵觀察**：
- **FP 集中**：HCC1395 TO mode (11,606 FP) + COLO829 paired_full (2,273 FP) = 13,879 FP = **所有樣本 FP 的 90%**
- **其他 5 個 paired_full FP 稀少**：加總 <1,400 FP
- **H1437/H2009 不可用**：baseline TP% ≥ 99.9%，scheme 改進無意義（頂板效應）

---

## 3.3 TO mode 資料來源

HCC1395 TO 的 TP / FP 不來自 master，而是源自：
- `/tmp/ism_hp_fix_phase1/tp_off/significance_summary.csv`
- `/tmp/ism_hp_fix_phase1/fp_off/significance_summary.csv`

這是 HP Fix Phase 1（2026-04-20）的 HP-only flag=off 重跑產物，確認 KDE 修正後 TO mode 的完整 ISM 特徵。

**限制**：只有 HCC1395 有 TO mode 的 ISM 特徵輸出，其他 6 樣本 TO mode 需 C++ pipeline rerun（見 `01_research_question.md` §1.5 非本研究範疇）。

---

## 3.4 Scheme Table 資料位置

| TSV | 樣本 | 內容 |
|-----|-----|------|
| `tpfp_stratified_filter_schemes.tsv` | paired_full 7 樣本 | S1-S7 n, n_tp, n_fp, tp_rate, ci |
| `tpfp_stratified_filter_schemes_TO.tsv` | HCC1395 TO | 同上 |
| `tpfp_per_sample_scheme_full.tsv` | 全部 sample-mode | E2 per-sample fold + CI + power |
| `tpfp_feature_diffs.tsv` | 全部 | E4 Mann-Whitney + Cliff's δ |
| `tpfp_subschemes.tsv` | HCC1395 TO | E5 S4 內 sub-scheme 候選 |
| `tpfp_5d_cube_HCC1395_TO.tsv` | HCC1395 TO | E6 5D cell table |
| `tpfp_5d_pareto_HCC1395_TO.tsv` | HCC1395 TO | E6 排序 + ≥20 filter |
| `tpfp_5d_cumulative_envelope_HCC1395_TO.tsv` | HCC1395 TO | E6 累積封套 |
| 同上 3 個 COLO829_pf 版本 | COLO829 paired_full | E6 次要測試場 |

---

## 3.5 Scheme 重疊統計（HCC1395 TO）

來自 `tables/obs08_scheme_overlap_matrix.tsv`：

| 行 ∩ 欄 | S3 | S5 | Top-17 | Top-28 |
|---------|----|----|--------|--------|
| **S3** | 380 | 380 | 337 | 366 |
| **S5** | 380 | 886 | 540 | 730 |
| **Top-17** | 337 | 540 | 1,099 | 1,099 |
| **Top-28** | 366 | 730 | 1,099 | 2,285 |

**重要解讀**：
- S3 (380) 幾乎完全被 Top-28 吸收（366/380 = 96.3%）
- Top-17 有 1,099 regions，但只 337 (31%) 屬 S3；**多出 762 regions 是 S3 捕不到的**
- Top-28 有 2,285 regions，但只 730 (32%) 屬 S5；**多出 1,555 regions 是 S5 捕不到的**

這證明 5D cube 不只是「切分得更細」，而是**發現 S3/S5 無法覆蓋的新高純度區**。

---

## 3.6 資料血脈（Provenance）

1. ISM C++ pipeline (2026-04-20) 產出 per-region significance_summary.csv（canonical/）
2. VCF-level intersubmod_tp/ + intersubmod_fp/ 作 TP/FP 標註
3. `step0_consolidate_dataset.py`（ng_kde_rescaling）合併 6 paired_full + HCC1954 rescaled + HCC1395 TO → master TSV
4. `step4/4b/6/8` 產 table, figures
5. 本研究 `obs01-08` 使用 master TSV + 現成 5D pareto TSV 作第二層觀察

**未記錄但可追溯**：master TSV 的 HCC1395 TO 行來自 Phase 1 HP-only 重跑；paired_full 6 樣本是 2026-04-20/21 Canonical baseline；HCC1954 是 2026-03-15 stale 版 + post-hoc KDE rescale。

---

## 3.7 版本標籤

- 本分析執行於：2026-04-22
- Master TSV 上次更新：2026-04-22（NG × KDE rescaling 研究）
- 5D cube TSV 產出於：2026-04-22 step8 execution
- KDE baseline：post-fix（2026-04-20 PR，CL-kde-fix-hcc1395）
