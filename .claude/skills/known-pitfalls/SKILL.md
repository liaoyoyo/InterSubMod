---
name: known-pitfalls
description: InterSubMod 已知 AI 陷阱清單。每條記錄具體錯誤、正確做法、觸發場景。避免重複犯錯。USE WHEN「known pitfalls」「踩雷」「avoid mistake」「common bug」「之前怎麼錯的」「歷史教訓」、涉及 OLS/residualization、VCF 來源、特徵設計、AUC 分析、binary commit / KDE fix / working tree、證據鏈 / single-track / ⭐4-5 升級。SKIP WHEN 純 build / commit / docs 寫作、無對應陷阱類別的場景（先查表，無命中即略）、新 feature 探索初期（先設計再查陷阱）、純 UI / 視覺 / PPTX 製作。
allowed-tools: Read, Write, Edit
user-invocable: true
---

# InterSubMod Known Pitfalls（AI 已知陷阱）

> 每條記錄來自過去對話中 AI 犯的具體錯誤。新陷阱發現時追加到對應分類下。

## 使用規則

| 觸發場景 | 必讀陷阱 |
|----------|---------|
| 涉及 OLS / residualization / confound control | P-01, P-02 |
| 涉及 VCF 來源識別 / 數據溯源 | P-03, P-04, P-12 |
| 涉及特徵設計 / AUC 分析 | P-05, P-06, P-09, P-10 |
| 涉及證據鏈 / single-track / ⭐4-5 升級 | P-07 (anchor #1) |
| 涉及 binary commit / KDE fix / working tree | P-08, P-13 |
| 涉及跨樣本 n_passed / saturation | P-11 |

---

## 統計方法陷阱

### P-01: L2 Collider Bias

**錯誤**：對 near-constant 特徵（如 AlleleDelta in LOH regions）做 OLS residualization on AF，產生虛假 AUC 信號（表面 AUC 從 0.50 跳到 0.59）。

**正確做法**：L2 residualized AUC 必須用 L3 AF-bin 交叉驗證。若 L2 與 L3 差距 > 0.10，即為 collider bias，該特徵應判定 CONFOUND。

**來源**：O12 LOH 甲基化場景分析。Memory: `feedback_L2_collider_bias.md`

### P-02: Pooled OLS Residualization Trap

**錯誤**：Pooled OLS（TP+FP 合併後 fit）殘差仍保留分組信息，因為 TP/FP 在特徵空間中佔據不同位置，殘差 = 組間差 + 組內差。

**正確做法**：必須使用 within-group OLS（TP/FP 分別 fit），殘差才真正移除 confound。

**來源**：Beyond-AUC M2 驗證。Memory: `feedback_pooled_ols_residualization_trap.md`

---

## 數據來源陷阱

### P-03: VCF 來源錯誤歸因

**錯誤**：將 canonical VCF 錯誤歸因為 "chenhan112 pipeline"。實際上 canonical TO pipeline VCF 是 liaoyoyo2001 使用 ONT_5kHz BAM（有 5mCG+5hmCG MM/ML）執行 ClairS-TO 產生的。

**正確做法**：確認 VCF 來源時必須追蹤：(1) 誰執行了 caller，(2) 使用哪個 BAM，(3) BAM 是否有 MM/ML tags。查閱 Knowledge/02_samples/ 和 Knowledge/03_file_formats/ 交叉驗證。

**來源**：2026-04-14 TO pipeline staging v2 修正。

### P-04: pileup Symlink 指向錯誤 Caller

**錯誤**：pileup 模式的 output symlink 實際指向 ClairS paired（非 TO），導致 TO 分析使用了 paired caller 的輸出。

**正確做法**：追蹤 symlink 實際目標（`readlink -f`），確認 caller pipeline 與分析模式匹配（TO 分析必須用 ClairS-TO VCF，paired 分析用 ClairS paired VCF）。

**來源**：Memory: `project_vcf_source_error_correction.md`

---

## 特徵分析陷阱

### P-05: CramersV 93% Zero Artifact

**錯誤**：將 CramersV 視為連續區分特徵使用。實際上 CramersV 在 2×2 contingency table（ISM 的標準框架）中只有 {0, 1} 兩個值，93% 的 regions 為 0。

**正確做法**：CramersV 不適合作為連續特徵使用。使用 HPFineNGroups（已克服此限制，AUC 提升 +0.125）作為替代。

**來源**：R1-R5 特徵設計研究。Memory: `project_feature_design_limitations_r1r5.md`

### P-06: n_reads / NumReads Confound

**錯誤**：忽略 read count 對所有統計量的系統性影響。較多 reads → 較高統計功效 → PERMANOVA p-value 更小、HPFineNGroups 更大，但這反映檢測力而非生物效應。

**正確做法**：所有特徵分析必須控制 n_reads（residualize 或分層）。任何 AUC > 0.58 的特徵都需排除 read count confound 後才能判定。

**來源**：O11 heterogeneity 分析。Memory: `project_O11_heterogeneity_negative.md`

---

## 證據鏈陷阱（v1.6 anchor #1 hardening）

### P-07: Single-Track Validated Cycle (Missing Orthogonal Evidence)

**錯誤**：cycle 標 ⭐4 / ⭐5 但缺第四軌「Orthogonal」證據（archive comparison / replicate run / alternate caller）。L4 mandatory 要求 4-track：(i) Statistical (ii) Cross-sample (iii) Mechanism (iv) Orthogonal — 缺 (iv) 是最常見的撤回原因（如 04-26 thread B whitelist）。

**正確做法**：宣告 ⭐4 / ⭐5 之前，至少加 1 個 orthogonal-track artifact 並寫入 `plan.preconditions.upstream_reports`。否則降為 ⭐3（described, single-track）。

**來源**：plan v1.6 §4.5.4-G batch 5a anchor #1；validation-protocol L4 mandatory；Drill 1 retro hpfinengroups + thread_b 案例驗證

---

## 數據新鮮度與整合陷阱（v1.8 T1-5 從 2026-04 retract events 編寫）

### P-08: KDE-fix Stale Binary Downstream

**錯誤**：binary commit (KDE fix) 之後，下游 master_*.tsv.gz dataset 仍由舊版 binary 重建，bias 持續到下游分析（如 HCC1395 S3 TP 95.5%→58.3% post-fix）。

**正確做法**：使用 master dataset 前必驗 dataset 重建時間 ≥ binary fix commit time。`/check-staleness` 已自動檢查 `precheck.checks.binary.stale_distance`；≥1 時 BLOCK，要求重建。

**來源**：20260420_KDE_Fix_Acceptance_Validation。Memory: `project_kde_fix_downstream_quantification.md`

### P-09: Spatial Autocorrelation Confound (chr+pos aggregation)

**錯誤**：以 chr+pos 視窗聚合的特徵會帶入 spatial autocorrelation（linkage / hotspot density / replication timing），AUC 看似顯著但 mid-TP-rate window 分層後消失。

**正確做法**：所有 chr+pos 聚合特徵必須跑 mid-TP-rate window stratification；若 AUC 下降 >0.05 → 判定 spatial autocorr confound，要求 within-window control。

**來源**：Memory: `feedback_spatial_autocorrelation_confound.md`

### P-10: Feature Name Literal Interpretation

**錯誤**：把含生物學語意的特徵名（例如 HPFineNGroups → "methylation subclone marker"）當成生物學意涵直接論證。實際 C++ 實作可能只是 phasing-derived 統計量，與 methylation 無直接關係（HPFineNGroups 04-22 撤回事件）。

**正確做法**：宣告生物學語意前必讀 `src/include/` 對應 feature 定義；plan.preconditions 應加 `source_code_refs` 欄列出 C++ 路徑。

**來源**：04-22 HPFineNGroups 撤回事件。Memory: `feedback_feature_name_vs_definition_rule.md`

### P-11: Saturation Artifact (1/N Samples Drives All)

**錯誤**：cross-sample n_passed=1/N（單一樣本飽和如 H2009 thread B 04-26）但 pooled / mean aggregation 隱藏這個 1-vs-rest pattern。pooled metric 看似 promising，實則由 1 樣本拉動。

**正確做法**：若 `generalize.consistency.n_samples_passed / n_samples_total ≤ 0.2` → BLOCK 不論 pooled metric 多好。深入分析哪個樣本驅動信號（saturation? batch effect?）。

**來源**：20260426_Thread_B_Whitelist_Retraction

### P-12: Merged Dataset AF Schema Drift

**錯誤**：merged_*.tsv.gz dataset 的 `AF` 欄位語意可能不等於 caller_af（5/7 樣本 archive: AF p75 < 0.06）。下游用 `AF` 當 caller_af 處理→錯誤結論。

**正確做法**：用 merged dataset 前驗 schema：`AF.p75 > 0.10` AND `dataset_id contains "caller_af_separate"`。否則用 canonical/ 拉 caller_af join 替換。

**來源**：20260424_X6_Caller_AF_S3S5_CrossSample。Memory: `feedback_merged_dataset_af_and_loh_pitfalls.md`

### P-13: Working Tree Dirty / Uncommitted Binary

**錯誤**：pipeline 跑時 git working tree dirty（modified/uncommitted）。結果不可重現；effect size 通常 near-noise（如 longphase 04-29 ΔF1=-0.0003）卻被當訊號。

**正確做法**：`/check-staleness` 已加 `precheck.checks.git.working_tree_clean` 檢查。dirty → BLOCK 直到 commit OR `binary_version` 寫成顯式 SHA + diff snapshot 歸檔。

**來源**：20260429_longphase_TO_vs_V5_Somatic_Fallback
