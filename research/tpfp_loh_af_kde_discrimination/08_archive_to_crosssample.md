---
title: Archive TO 資料跨樣本觀察（已由 2026-04-22 ISM rerun 取代）
status: superseded
superseded_by: 07_figure_layers.md + data/cell_consistency_matrix.tsv
superseded_at: 2026-04-22 14:40
owner: liaoyoyo2001
created: 2026-04-22
last_updated: 2026-04-22 14:55
parent_project: 00_INDEX.md
---

# 08 · Archive TO 資料跨樣本觀察（SUPERSEDED）

> ⚠ **此文件之結論與限制已於 2026-04-22 14:40 完成的 ISM rerun 後被取代**。
> 當時本檔記錄的「LOH/CN 欄位 NaN」「無法支援完整 L3」等限制已全部解除。
> **請改看 `07_figure_layers.md`** 的 rerun 後分析。
> 本檔保留作為**研究歷程紀錄**，以及說明「為何要做 obs15 ISM rerun」的動機依據。

---

## 取代摘要（rerun 前 ❌ vs 後 ✅）

| 項目 | rerun 前（03:xx）| rerun 後（14:49）|
|------|:----------------:|:---------------:|
| 5 archive TO 樣本資料來源 | ❌ label_first_metrics.tsv (36 欄) | ✅ significance_summary.csv (117 欄, post-HP-fix) |
| LOH_Subtype 欄位 | ❌ 全 NaN | ✅ 五類別完整 None/Noise/Weak/Strong/Subclone |
| Coverage_Multiple / Diploid_Coverage_Used | ❌ NaN | ✅ post-KDE-fix（DiplCov 53-91x）|
| CN tier / AF_class 衍生 | ⚠ 部分支援（AF only）| ✅ 完整五維 cube 可分析 |
| Phase 2 features（NHP_Somatic, Fisher_Frac_Sig, PerCpgASM_Valid, NME_HP1/2）| ❌ 不可得 | ✅ 完整可得 |
| 下游分析覆蓋 | ❌ L1+L2 有限（AF × NG）| ✅ **L1 + L2 + L3 + Consistency 完整**（obs09-12 V2 重跑）|
| master dataset | ⚠ master_extended.tsv.gz (含 NaN 欄) | ✅ merged_7samples_paired_full_plus_hcc1395_to.tsv.gz (16.9 MB, obs16 重生成) |
| 圖表設計 | ⚠ label-first technical titles | ✅ research-question suptitles + takeaway captions + metadata badges |

**新主 master**：`/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz`（由 `step0_build_master.py` 於 2026-04-22 14:40 重建）。

---

## 歷史背景（為後續人員保留 context）

`master.tsv.gz` 的 `to_pileup` 模式之前**只有 HCC1395 一個樣本**（Phase 1 HP-only 工作區的 post-fix 資料），導致 §7 分層分析（L1-L3、consistency）無法回答跨樣本 TO 問題。

在 `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/` 下，實則**已經跑過**下列 5 樣本的 TO pipeline（ClairS-TO → LongPhase-TO → InterSubMod），但**使用的是 pre-HP-fix 的 ISM binary**：

| 樣本 | Pilot 目錄 | 時間 | archive 狀態 |
|------|-----------|-----:|------|
| HCC1395_DORADO | `20260315_hcc1395_dorado_to_pilot` | 2026-03-30 (post-HP-fix binary) | 已彙整 |
| H1437 | `20260318_h1437_to_pilot_fastresume` | 2026-03-30 | 已彙整 |
| H2009 | `20260318_h2009_to_pilot_fastresume` | 2026-03-30 | 已彙整 |
| HCC1937 | `20260318_hcc1937_to_pilot_fastresume` | 2026-03-30 | 已彙整 |
| HCC1954 | `20260318_hcc1954_to_pilot_fastresume` | 2026-03-30 | 已彙整 |
| COLO829 | `20260317_colo829_to_pilot` | step05 空 | SKIPPED |

這 5 個 archive pilot 在 2026-03-30 當時**未保存完整 `significance_summary.csv`**（僅保留 `label_first_metrics.tsv` 的 36 欄），這是因為當時研究重點在 paired mode，TO 僅作為 pilot 參考；此外 archive 時期的 ISM binary 亦未產生 KDE (Diploid_Coverage_Used) 與 Phase 2 feature 欄位（該二進位於 2026-04-21 重編後才支援）。

---

## 原 obs13/14 彙整的結果（僅供歷史比對）

**原本限制條件**：archive 的 `label_first_metrics.tsv`（36 欄）**無法直接產生**當前 master 的 LOH/CN/KDE 相關欄位，因此只能支援 L1-L2 AF × HPFineNGroups × NumReads 的跨樣本觀察。

原跨樣本觀察（基於 label_first_metrics，**以下結論經 2026-04-22 rerun 後，有部分被更新**）：

### 原觀察 1：TO baseline TP% 跨樣本差異極大（25%-91%）

| 樣本 × TO | n_total | n_TP | n_FP | baseline TP% |
|-----------|-------:|-----:|-----:|------:|
| H2009 | 137,695 | 125,706 | 11,989 | **91.3%** |
| H1437 | 58,915 | 45,473 | 13,442 | 77.2% |
| HCC1395_DORADO | 40,428 | 28,856 | 11,572 | 71.4% |
| HCC1395 | 40,115 | 28,509 | 11,606 | 71.1% |
| HCC1954 | 67,286 | 17,068 | 50,218 | **25.4%** |
| HCC1937 | 24,655 | 12,623 | 12,032 | 51.2% |

✅ **仍然有效**，數字在新 master 中完全一致（rerun 並未改變變異結構，僅補齊 LOH/CN 欄位）。

### 原觀察 2：HCC1954 AF 方向反向

舊結論：HCC1954 Extreme AF 86% vs Near-half 14%，與其他樣本相反。

✅ **仍然有效**，在 rerun 後 §07 §L1.2 / §L2.1 確認。新增：其他 4 TO 樣本（HCC1395_DORADO, HCC1937, H2009, H1437）在 AF 方向上與 HCC1395_TO 一致。

### 原觀察 3：HCC1395 vs DORADO Extreme AF 佔比 80.9% vs 0.4%

✅ **仍然有效**。此跨-basecaller 差異對 filter 移植性仍是關鍵警訊。

### 原觀察 4：NG≥4 在 5/6 TO 樣本 ≥69%

✅ **已升級為更強結論**。rerun 後 §07 §L1.4 顯示 NG 梯度不僅在 NG≥4 穩健，**NG=2 vs NG=3 之間的遞增也跨 5/6 TO 樣本一致**，這讓 NG 維度從「部分穩健」升級為「核心穩健 discriminator」。

### 原限制：無法跑完整 L3（LOH/CN 缺）

❌ **已解除**。2026-04-22 14:40 的 ISM rerun 補齊所有 LOH/CN/KDE 欄位，`07_figure_layers.md` 的 4 張 L3 + 2 張 consistency 圖已全數可繪。

---

## 衍生輸出（仍可參考）

本檔原本產出的檔案多數仍保留在 `figures/new/archive_to/` 下，可用作 **pre-rerun vs post-rerun** 的對照：

| 檔名 | 備註 |
|------|------|
| `figures/new/archive_to/fp_scale_comparison.png` | paired vs TO FP 膨脹比（18× to 1731×）— **仍是有效結論** |
| `figures/new/archive_to/L1_af_class_to.png` | 6 TO 樣本 AF × TP-rate — 與新 `figures/new/L1/L1_af_TO.png` + `L1_af_paired.png` 結論一致 |
| `figures/new/archive_to/L1_hpfinengroups_to.png` | 6 TO 樣本 NG × TP-rate — 與新 `figures/new/L1/L1_ng_TO.png` + `L1_ng_paired.png` 結論一致，**後者更新更強** |
| `figures/new/archive_to/L2_af_x_ng_to.png` | 6 TO 樣本 AF × NG heatmap — 已被 `L2_af_x_ng.png` 取代（新版含 13 sample-mode）|
| `data/master_extended.tsv.gz` | pre-rerun extended master（含 NaN 欄）— **保留作 provenance 用**，實際分析改用新 master |

**新分析主 output**（rerun 後）：
- Master: `ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz`
- Consistency: `tpfp_loh_af_kde_discrimination/data/cell_consistency_matrix.tsv`（408 cells × 36 欄）
- Figures: `figures/new/L1|L2|L3|consistency/` 共 21 張（每張涵蓋 13 sample-mode）
- 圖檔更新: `figures/existing/fig_v2_7_baseline_tp_fp_fold.png`（6 TO panels）

---

## 研究追溯鏈

1. **2026-04-22 ~03:00** — 發現 master TSV 僅 HCC1395_TO 一個 TO 樣本，手動彙整 5 archive pilot 的 `label_first_metrics.tsv` 得 `master_extended.tsv.gz`（即本檔原始分析）
2. **2026-04-22 ~04:00** — 撰寫本檔，記錄 LOH/CN 缺失與跨樣本原始觀察
3. **2026-04-22 ~05:00** — 決定執行 ISM rerun，撰寫 `20260422_Archive_TO_Rerun_ISM_Requirement_01.md` 需求文件
4. **2026-04-22 10:52** — 啟動 `obs15_rerun_archive_to_ism.sh`（5 pilot parallel，9 threads each）
5. **2026-04-22 12:06** — obs15 driver 被 session compaction 殺掉；DORADO 已完成
6. **2026-04-22 12:17** — 啟動 `obs15b_resume_failed_pilots.sh`（detached，4 pilot 續跑）
7. **2026-04-22 14:30** — 5 pilots 全數完成（TP + FP 兩輪 ISM）
8. **2026-04-22 14:40** — `obs16_regen_after_ism_rerun.sh` 完成（step0→step7 regen）
9. **2026-04-22 14:49** — obs09-12 parallel rerun 完成（21 張新圖 + consistency matrix）
10. **2026-04-22 14:55** — 本檔標記 superseded，`07_figure_layers.md` 改寫為 post-rerun 權威版本

---

## 後續（參考用）

LOSO 驗證目標應改從原計畫「HCC1395 Top-17」改為**新 16 cells 共識集**（見 `07_figure_layers.md` §7.5.1 與 §7.7）。本檔不再需要維護；所有分析接續在 `07` 與 `00_INDEX.md` 進行。
