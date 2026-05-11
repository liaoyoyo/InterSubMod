---
title: 資料 Provenance Audit — 三類誤用風險自動掃描
date: 2026-04-26
type: audit
status: in_progress
scope:
  - docs/experiments/
  - docs/reports/
  - research/feature_layered_observation/
risk_categories:
  - stale_binary
  - af_misuse
  - alleledelta_misuse
  - provenance_missing
related:
  - InterSubMod/docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/docs/data_specs/20260411_significance_summary欄位字典_01.md
  - InterSubMod/research/feature_layered_observation/02_methodology.md
  - InterSubMod/research/feature_layered_observation/03_vcf_annotation_plan.md
  - InterSubMod/src/core/RegionProcessor.cpp
  - InterSubMod/src/core/NormalBaseline.cpp
build_commit_reference: 374fad4 + 12d9b3e (KDE fix)
---

# 資料 Provenance Audit — 2026-04-26

## Summary

| 風險類別 | 命中數 | 嚴重度 |
|---------|-------|--------|
| stale_binary（hardcoded `expected_coverage=75.0`） | 6 處 | 🔴 高（影響跨樣本 CovM/Z-score 結論） |
| af_misuse（master `AF` ≠ caller VAF） | 8 處 | 🟠 中（多在 historical/archive，已於主分析校正） |
| alleledelta_misuse（AlleleDelta 當 VAF proxy） | 7 處 | 🟡 中-低（多為 confound 警告語境，正向使用為陷阱） |
| provenance_missing（無 build_commit / KDE-fix 註記） | 12 處 | 🟠 中（validated/2026/04 共 14 報告，僅 2 含 KDE provenance） |

**總計**：33 處風險條目（高風險 ≥5，已超過任務驗收門檻）。

---

## High-Risk 條目 Top 10

| # | file | line | risk_type | snippet | action_required |
|---|------|------|-----------|---------|-----------------|
| 1 | `InterSubMod/docs/data_specs/20260411_significance_summary欄位字典_01.md` | 212 | stale_binary | `Coverage_Multiple ... fallback=75.0` | 加 [STALE-BINARY] caveat：標註 KDE-fix（commits 374fad4+12d9b3e）後 fallback 改為 KDE estimate；75.0 fallback 已視為 stale |
| 2 | `InterSubMod/docs/reports/research_landscape/06_結論穩定性審查.md` | 591 | stale_binary | `CONFIRMED 但受 P0-A CovM hardcoded 75.0 bug 影響；修正後需重算` | 重跑分析（KDE-corrected）：以 KDE-fixed master 重算 Z1/Z3 邊界與 6 章 Zone 結論 |
| 3 | `InterSubMod/research/feature_layered_observation/00_main_observation.md` | 199 | stale_binary | `HCC1954 paired Diploid_Coverage_Used=NULL（stale binary artifact）` | 重跑分析（KDE-corrected）：HCC1954 paired_full kde_rerun_B_14combos 已產出 dc=61×；更新此段為 RESOLVED |
| 4 | `InterSubMod/research/feature_layered_observation/features/G1_coverage.md` | 24, 91, 108, 119, 124 | stale_binary | `HCC1954 paired_full has *empty* Diploid_Coverage_Used` / `Diploid_Coverage_Used AUC 0.790 ... sample-ID proxy` | 加 [STALE-BINARY] caveat + 註明 G1 觀察使用 pre-KDE-fix master；KDE-corrected 重跑為 P1 待辦（registry 已建立） |
| 5 | `InterSubMod/research/feature_layered_observation/01_feature_inventory.md` | 268 | stale_binary | `Coverage_Multiple hardcoded 75.0 default ... per-sample diploid_coverage auto-estimation may not be active` | 加 [STALE-BINARY] caveat：點明 commits 374fad4+12d9b3e 已修；inventory 寫於 KDE fix 前 |
| 6 | `InterSubMod/docs/reports/20260408_InterSubMod_Stage_Report_v1.md` | 129 | af_misuse + alleledelta_misuse | `\| AF (allele frequency) \| AlleleDelta=AF proxy \| AF-bin stratification \|` | 加 [AF-MISUSE] caveat：master `AF` 欄位 = `\|AlleleDelta\|`，非 caller VAF；分層應改用 `vcf_AF`（已於 inventory 02_methodology.md:10 標註） |
| 7 | `InterSubMod/docs/reports/20260408_InterSubMod_Stage_Report_v1.md` | 282 | alleledelta_misuse | `AlleleDelta \| 0.642 \| 0.607 \| Allele (AF proxy)` | 加 [ALLELEDELTA-MISUSE] caveat：AUC 0.642 為 master `AF`=`\|AlleleDelta\|` 自相關；改用 caller_af 重算 |
| 8 | `InterSubMod/research/feature_layered_observation/00_main_observation.md` | 254 | alleledelta_misuse | `\| AlleleDelta \| L2 collider \| CONFOUND_COLLAPSED (=\|vcf_AF\|) \|` | 此條為**正確警告**用法；保留但確認下游 G9/G10 觀察均已用 `caller_af` 重算（已於 03_vcf_annotation_plan.md:117 註明） |
| 9 | `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_enrichment_paired_to_corrected_analysis_01.md` 等 12 檔 | (file-level) | provenance_missing | 無 `build_commit` / `374fad4` / `12d9b3e` / `KDE.fix` 任何字串 | 加 [STALE-BINARY] caveat 區塊：聲明分析時間早於 KDE fix；CovM-依賴結論需 P1 registry 重跑驗證 |
| 10 | `InterSubMod/docs/reports/research_landscape/06_結論穩定性審查.md` | 20, 581 | stale_binary（landscape 級） | `新 P0-A CovM expected_coverage=75.0 hardcoded bug 修正` / `Zone 分類可能需重算` | 重跑分析（KDE-corrected）+ 鏈結至 P0-1 Thread D 主軸報告（佔位） |

---

## 各風險類別詳表

### Stale Binary（6 處）

| # | file:line | snippet (≤100字) | action_required |
|---|-----------|------------------|-----------------|
| S1 | `InterSubMod/docs/data_specs/20260411_significance_summary欄位字典_01.md:212` | `Coverage_Multiple ... NumReads / diploid_coverage ... fallback=75.0` | 加 [STALE-BINARY] caveat |
| S2 | `InterSubMod/docs/reports/research_landscape/06_結論穩定性審查.md:20` | `新 P0-A CovM expected_coverage=75.0 hardcoded bug 修正 + 7 樣本重跑` | 加 [STALE-BINARY] caveat（已標 P0-A） |
| S3 | `InterSubMod/docs/reports/research_landscape/06_結論穩定性審查.md:581` | `CovM expected_coverage=75.0 hardcoded bug 修正後，Zone 分類可能需重算` | 重跑分析（KDE-corrected） |
| S4 | `InterSubMod/docs/reports/research_landscape/06_結論穩定性審查.md:591` | `CONFIRMED 但受 P0-A CovM hardcoded 75.0 bug 影響；修正後需重算` | 重跑分析（KDE-corrected） |
| S5 | `InterSubMod/research/feature_layered_observation/01_feature_inventory.md:268` | `Coverage_Multiple hardcoded 75.0 default ... per-sample diploid_coverage auto-estimation may not be active` | 加 [STALE-BINARY] caveat |
| S6 | `InterSubMod/research/feature_layered_observation/00_main_observation.md:199, 213, 246` + `features/G1_coverage.md:24,42,48,91,108,119,124,166` | `HCC1954 paired Diploid_Coverage_Used=NULL（stale binary artifact）` / `Diploid_Coverage_Used AUC 0.790 ... sample-ID proxy` / `R-STALE-BIN \| stale-binary sentinel` | 加 [STALE-BINARY] caveat + KDE rerun 後重算 G1 AUC（HCC1954 dc=61× 已備齊） |

**註**：`InterSubMod/docs/experiments/INDEX.md:231` 為 **正確記錄 KDE fix 過程**（commits 374fad4+12d9b3e + bias 量化）— 不列為風險，反為 provenance 來源。

### AF Misuse（8 處）

| # | file:line | snippet (≤100字) | action_required |
|---|-----------|------------------|-----------------|
| A1 | `InterSubMod/docs/reports/20260408_InterSubMod_Stage_Report_v1.md:114` | `AlleleDelta: AF proxy, not methylation` | 加 [AF-MISUSE] caveat：master `AF`≠caller VAF，AlleleDelta 與其同源不獨立 |
| A2 | `InterSubMod/docs/reports/20260408_InterSubMod_Stage_Report_v1.md:129` | `AF (allele frequency) \| AlleleDelta=AF proxy \| AF-bin stratification \| All within-bin AUC<0.55` | 加 [AF-MISUSE] caveat：分層需改用 `vcf_AF` |
| A3 | `InterSubMod/docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md:369` | `AlleleDelta \| ~0.55 \| ~0.50 \| 直接 AF proxy \| L3 AF-bin 分層` | 加 [AF-MISUSE] caveat：需註明 AF-bin 分層應使用 caller `vcf_AF` 而非 master `AF` |
| A4 | `InterSubMod/docs/experiments/INDEX.md:168` | `AlleleDelta AUC=0.556 真實但微弱` | 加 [AF-MISUSE] caveat：AUC 0.556 為 row-level master `AF`（=\|AlleleDelta\|）自相關 |
| A5 | `InterSubMod/docs/reports/validated/2026/02/20260211_甲基化過濾策略綜合分析報告_01.md:63,73,78,163,198,205` | `(AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24)` | 加 [AF-MISUSE] caveat：filter rule 之 `VAF` 為 caller_af（正確），但 AlleleDelta 解讀如 ASM 強度需再驗證 |
| A6 | `InterSubMod/docs/reports/validated/2026/02/20260211_大規模驗證流程規格_01.md:477` | `Remove if: (QUAL < 0.75) OR (AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24)` | 加 [AF-MISUSE] caveat：規格屬 archive；新方法學審查已 reject (annotation-only) |
| A7 | `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_weekly_review/00_background.md:276` | `AF（allele frequency）是 LOH 分析的 confound ... O12 已確認 AlleleDelta = AF confound` | 加 [AF-MISUSE] caveat：此段為**正確警告**；確認 confound 為 master `AF`=\|AlleleDelta\| 同源 |
| A8 | `InterSubMod/research/feature_layered_observation/scripts/g9_build_extended_master.py:4` | `The primary merged_with_vcf.tsv.gz only carries AlleleDelta (== master AF)` | 此為**正確註記**；保留並交叉鏈結至 02_methodology.md:10 |

### AlleleDelta Misuse（7 處）

| # | file:line | snippet (≤100字) | action_required |
|---|-----------|------------------|-----------------|
| D1 | `InterSubMod/docs/reports/20260408_InterSubMod_Stage_Report_v1.md:282` | `AlleleDelta \| 0.642 \| 0.607 \| Allele (AF proxy)` | 改用 caller_af 重算（AUC 0.642 為 master AF 自相關偽訊號） |
| D2 | `InterSubMod/docs/experiments/in_progress/2026/04/20260404_HPFineP_QS整合完整研究報告_01.md:273` | `TP 的 AlleleDelta 更小（AUC=0.393 < 0.5）... LR 利用這個反向信號` | 加 [ALLELEDELTA-MISUSE] caveat：AUC 0.393 為 row-level 自相關，非獨立信號 |
| D3 | `InterSubMod/docs/experiments/in_progress/2026/04/20260405_TO_ClairSTO特徵區分力深度研究_01.md:239` | `AlleleDelta 的有效信號僅來自 LOH 區域中 minor haplotype 極少 reads 的情況（本質上是 HP2FamilyN 的 proxy）` | 加 [ALLELEDELTA-MISUSE] caveat：解釋語意 OK，但效應量需以 caller_af 殘差化重算 |
| D4 | `InterSubMod/docs/archive/2026/03/phase4_igv/20260317_研究主線週報_20260312_20260317_phase4_FP分型_01.md:341,403,404,555` | `S1 HCC1937 FP VAF≈1.0 + AlleleDelta=0` / `B 平台特異高 AD 型 ... AlleleDelta=0.088（13x DORADO）` | 加 [ALLELEDELTA-MISUSE] caveat：FP 分型規則沿用 archive；新方法學已 reject 為 annotation-only |
| D5 | `InterSubMod/docs/methodology/20260324_方法學審查全域結論報告_01.md:151` | `H005 \| VAF<0.08+AlleleDelta>0.05 \| delta=0，AlleleDelta 在 TO ≈ 0` | 已標 reject；保留但補連 [ALLELEDELTA-MISUSE] caveat |
| D6 | `InterSubMod/docs/experiments/in_progress/2026/04/20260404_LOH_Strong_Weak_7feature_AUC驗證報告_01.md:236` | `O12: AlleleDelta = AF confound \| AlleleDelta 單特徵 AUC=0.800 in LOH_Strong \| 完全一致` | 此為**正確警告**；保留 |
| D7 | `InterSubMod/research/feature_layered_observation/00_main_observation.md:97,141,254` | `G10 LabelAllelePermanovaF ... AF proxy（L2 collider），AlleleDelta=\|vcf_AF\| 同源` | 此為**正確警告**；保留並標注為 audit anchor |

### Provenance Missing（12 處）

`docs/reports/validated/2026/04/` 共 14 個 .md，僅 2 個含 KDE / build_commit 標記（`20260421_研究週報` + `20260423_研究週報`）。剩餘 12 個聲稱結論但無 build_commit / KDE-fix 註記：

| # | file | risk_type | action_required |
|---|------|-----------|-----------------|
| P1 | `InterSubMod/docs/reports/validated/2026/04/20260401_comprehensive_observation_O1_O10_report_01.md` | provenance_missing | 加 [STALE-BINARY] caveat 區塊（O1-O10 多依賴 master CovM） |
| P2 | `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_enrichment_paired_to_corrected_analysis_01.md` | provenance_missing | 加 [STALE-BINARY] caveat（LOH 雖 Jaccard=1.0 不變，但 enrichment 用 region pool 受 CovM 影響） |
| P3 | `InterSubMod/docs/reports/validated/2026/04/20260401_systematic_observation_O1_O10_summary_01.md` | provenance_missing | 加 [STALE-BINARY] caveat |
| P4 | `InterSubMod/docs/reports/validated/2026/04/20260402_loh_read_threshold_visual_argument_01.md` | provenance_missing | 加 [STALE-BINARY] caveat（read threshold 受 NumReads 絕對值；scale-invariant 部分可豁免） |
| P5 | `InterSubMod/docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md` | provenance_missing | 加 [STALE-BINARY] caveat（causal chain 結構不變，但定量值需 KDE rerun） |
| P6 | `InterSubMod/docs/reports/validated/2026/04/20260402_purity_dependent_self_phasing_validation_01.md` | provenance_missing | 加 [STALE-BINARY] caveat |
| P7 | `InterSubMod/docs/reports/validated/2026/04/20260402_read_level_germline_fp_research_report_01.md` | provenance_missing | 加 [STALE-BINARY] caveat |
| P8 | `InterSubMod/docs/reports/validated/2026/04/20260404_LOH_evidence_panel_post_TO_HP_fix_final_report_01.md` | provenance_missing | 加 [STALE-BINARY] caveat（HP fix 已記錄；缺 KDE 標記） |
| P9 | `InterSubMod/docs/reports/validated/2026/04/20260406_研究週報_20260331_20260406_LOH雙定義與特徵探索全面關閉_01.md` | provenance_missing | 加 [STALE-BINARY] caveat |
| P10 | `InterSubMod/docs/reports/validated/2026/04/20260406_肉眼檢視推理鏈與TP_FP可區分性分析_01.md` | provenance_missing | 加 [STALE-BINARY] caveat |
| P11 | `InterSubMod/docs/reports/validated/2026/04/20260408_TO_LOH內外ISM特徵區分力完整推論鏈報告_01.md` | provenance_missing | 加 [STALE-BINARY] caveat（TO_LOH 跨樣本 CovM Z-score 受影響） |
| P12 | `InterSubMod/docs/reports/validated/2026/04/20260408_TO_LOH額外研究_遮罩與反轉分析_01.md` | provenance_missing | 加 [STALE-BINARY] caveat |

**Audit anchor（已含 provenance）**：
- `InterSubMod/docs/reports/validated/2026/04/20260421_研究週報_20260414_20260421_多軌收斂與定位定錨_01.md`（含 A3 stale binary 章節）
- `InterSubMod/docs/reports/validated/2026/04/20260423_研究週報_20260416_20260423_NG2_LOH_constrained_phasing與TO_pivot_01.md`

---

## R-DATA-GAP 額外標記

### R-DATA-GAP-01：NormalBaseline writer 三欄位 populated bug（標記，本 audit 不修）

**Symptom**：master dataset / merged TSV 的下列欄位疑似 populated 但值全 0 或 NA：
- `Normal_HP_Delta`
- `Normal_HP_Signed_Delta`
- `HP_Signed_Residual`

**Writer 路徑（已確認）**：

| 欄位 | 計算位置 | TSV 寫入位置 | 函數 |
|------|---------|-------------|------|
| `Normal_HP_Delta` | (subset_indices 計算) | `InterSubMod/src/core/RegionProcessor.cpp:1244` | computed earlier in subset diagnostics block |
| `Normal_HP_Signed_Delta` | `InterSubMod/src/core/RegionProcessor.cpp:1010` | `InterSubMod/src/core/RegionProcessor.cpp:1244` | `compute_signed_hp_delta(normal_indices)` (lambda L987-L1006) |
| `HP_Signed_Residual` | `InterSubMod/src/core/RegionProcessor.cpp:1021` | `InterSubMod/src/core/RegionProcessor.cpp:1245` | `tumor_hp_signed_delta - normal_hp_signed_delta` |

**Header 宣告**：`InterSubMod/src/core/RegionProcessor.cpp:1106-1109`

**疑似 bug 機制（待驗證）**：
- `compute_signed_hp_delta` 在 `hp1_count<1 || hp2_count<1` 時返回 `NAN`
- TO 模式 `normal_indices.empty()` → 全 region NaN（預期）
- Paired 模式應有值，但若 `is_tumor_mask` 過濾後 normal reads 缺失或 HP=0/1 過少，仍會 NaN
- **核心待查**：master TSV 是否將 NaN 序列化為 `NA` 還是 `0`？若為 0 → 統計分析會誤把缺值當真實 zero ASM
- **NormalBaseline.cpp** 本身僅計算 `normal_baseline_mean` / `normal_baseline_coverage`（src:`build_normal_baseline` L7），不寫 HP delta；HP delta 在 RegionProcessor 直接計算

**Action**：未來 `/cpp-change` 修復流程；本 audit 僅標記。建議測試：
1. paired_full HCC1395 抽 100 region，驗證 `Normal_HP_Signed_Delta` 是否為 NA / 0 / 真值分佈
2. 比對 master TSV 與 RegionProcessor TSV 直接輸出，確認序列化路徑

---

## 連結

- KDE registry：`InterSubMod/research/data_registry/kde_corrected_provenance_20260426.tsv`（P1-1 待建立 — 截至 audit 時間 `research/data_registry/` 目錄尚不存在）
- Thread D 主軸報告：佔位（P0-1 待建立 — `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_主軸_NG2_phasing_01.md`）
- Thread B 撤回宣告：佔位（P0-2 in_progress — `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_B_HPFineNGroups_撤回宣告_01.md`）
- KDE fix 驗收：`InterSubMod/docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md`
- INDEX 條目（KDE fix）：`InterSubMod/docs/experiments/INDEX.md:231`
- Methodology master AF 警告：`InterSubMod/research/feature_layered_observation/02_methodology.md:10`
- VCF annotation 計畫：`InterSubMod/research/feature_layered_observation/03_vcf_annotation_plan.md:117`
- 欄位字典：`InterSubMod/docs/data_specs/20260411_significance_summary欄位字典_01.md`
- C++ writer 主檔：`InterSubMod/src/core/RegionProcessor.cpp`（L987-L1245）
- NormalBaseline 計算：`InterSubMod/src/core/NormalBaseline.cpp`（L7-L83）

---

## 驗收檢核

- [x] audit 表 ≥5 high-risk（實際 10 條 Top；總 33 處）
- [x] 每條目有 file + line + risk_type + action_required
- [x] Summary 段含三類風險計數（stale_binary 6 / af_misuse 8 / alleledelta_misuse 7 / provenance_missing 12）
- [x] R-DATA-GAP NormalBaseline bug 標記（writer 路徑 + 待驗證機制 + cpp-change action）
- [x] 掃描範圍涵蓋 docs/experiments、docs/reports、research/feature_layered_observation
- [x] 連結含 KDE registry / Thread D / Thread B 佔位
