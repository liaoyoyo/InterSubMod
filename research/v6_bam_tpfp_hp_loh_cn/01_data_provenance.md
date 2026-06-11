<!--
建立時間: 2026-05-15
更新時間: 2026-05-15
目標: V6 BAM TPFP HP-LOH-CN 計畫所有資料來源與 commit hash 的單一彙整檔
資料來源:
  - InterSubMod/research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3)
  - InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md
  - InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md
  - InterSubMod/research/paired_priority_bug_audit/09_V6_caller_F1_verification.md
  - Explore agent 檔案存在性核對（2026-05-15）
狀態: in_progress
audience: 後續 step1-4 agent + reviewer
-->

# 01 — Data Provenance（資料來源彙整）

> 本檔記錄計畫所有依賴資料的**確切路徑、commit hash、生成日期、驗證狀態**。Step 1-4 agent 引用時必查本檔，不憑記憶。

## 1. Longphase Binary 三版本

| 版本 | Binary 路徑 | Commit / 日期 | 修補狀態 | 預期 phasing 行為 |
|------|-----------|-------------|---------|----------------|
| **V3F** (PON-only V3-Fixed) | 在 V5 之前的版本（commit 確切 hash agent 待補） | ~2026-04 完成 | upstream PON-only + Layer 1.5 未加 + germline-absent → hp=33 保守 | hp=1-1:hp=2-1 ratio = **1.138**（接近中性）|
| **V5** (Layer 1.5) | `/big7_disk/liaoyoyo2001/longphase-to-mod/longphase-to` | commit `938f0df` (2026-05-10 20:59, V5 HEAD) | V3F + 加 Layer 1.5 somatic fallback (HaplotagProcess.cpp:537-548) | hp=1-1:hp=2-1 ratio = **1.86**；germline-absent 區繼承 4.19:1 priority bug |
| **V6** (production candidate) | `/big7_disk/liaoyoyo2001/longphase-to-mod/longphase-to-v6` (22.55 MB) | V5 上單點 patch（-30 KB due to removed code）；2026-05-10 17:13 完成 | V5 + 移除 Layer 1.5 (HaplotagProcess.cpp:537-548)；重用 V5 phased VCF | hp=1-1:hp=2-1 ratio = **1.838**（接近 V3F；germline-existent 區因 V5 phased VCF 殘留）|

**V6 patch 設計（07_V6_validation_findings.md §1.1）**：
```diff
 void HaplotagProcess::getVote(...){
-    // Three-layer haplotype determination:
-    // Layer 1.5: Somatic fallback (HP1_1 vs HP2_1) — when germline absent
+    // V6 two-layer haplotype determination (Layer 1.5 reverted):
+    // V6: germline absent → conservative ambiguous (V3F behavior)
     ...
-    // Layer 1.5: Somatic fallback
-    else if (somaticHP1 > 0 || somaticHP2 > 0) { ... }
+    else { min = 0; max = 0; }
}
```

## 2. Tagged BAM 三版本（HCC1395 TO）

| BAM | 路徑 | 大小 | 用途 |
|-----|------|------|------|
| V3F BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam` | ~268 GB | V3F baseline 三方對照組 |
| V5 BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_tagged.bam` | ~268 GB | V5 Layer 1.5 |
| V6 BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam` | 268 GB | V6 production candidate |

**Note**：本計畫**不直接讀 BAM**，只讀已產出的 ISM 結果 TSV（節省 IO + 確保再現性）。

## 3. ISM 三向結果（HCC1395 — phaseC_genome_three_way）

**路徑根**：`/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way/`

**12 個 ISM run 目錄**（3 BAM × 2 flag × 2 label）：

| Dir | BAM | flag | label | 預期 region 數 |
|-----|-----|------|-------|--------------|
| V3F_off_tp | V3F | --germline-hp-only=off | TP | ~30,490 |
| V3F_off_fp | V3F | off | FP | ~4,842 |
| V3F_on_tp | V3F | --germline-hp-only=on | TP | ~30,490 |
| V3F_on_fp | V3F | on | FP | ~4,842 |
| V5_off_tp | V5 | off | TP | 同上 |
| V5_off_fp | V5 | off | FP | 同上 |
| V5_on_tp | V5 | on | TP | 同上 |
| V5_on_fp | V5 | on | FP | 同上 |
| V6_off_tp | V6 | off | TP | 同上 |
| V6_off_fp | V6 | off | FP | 同上 |
| V6_on_tp | V6 | on | TP | 同上 |
| V6_on_fp | V6 | on | FP | 同上 |

**每 dir 結構**：
```
{V*}_{on,off}_{tp,fp}/
├── debug/
├── filtered_snv_{tp,fp}/
│   └── {chr1-22,X,Y}/
│       └── {chr}_{pos-5000}_{pos+5000}/   ← region_id format
│           ├── reads.tsv                   ← read-level (read_id, name, chr, start, end, mapq, hp, alt_support, is_tumor, strand)
│           ├── methylation.csv             ← Read × CpG matrix
│           ├── distance_matrix_NHD.csv
│           ├── significance_summary.csv    ← per-region ISM 特徵 (~100 欄位)
│           └── *.png                       ← visualization
├── hp_summary.txt
├── run.log
├── significance_statistics.txt
└── significance_summary.csv                ← **HEADER-ONLY** (0 regions aggregated; 真實 per-region 在子目錄)
```

⚠️ **重要陷阱**：`{dir}/significance_summary.csv` 是 header-only（0 rows），真實 per-region data 在 `filtered_snv_*/chr*/region_id/significance_summary.csv`。Agent A 已在 build_three_way_master.py 處理。

## 4. Aggregate Summary TSV

| 檔案 | 行數 | 內容 |
|------|------|------|
| `phaseC_genome_three_way/v3f_vs_v5_vs_v6_genome_summary.tsv` | 4 (header + 3 BAM) | BAM-level summary：marker_off_n / marker_off_rate / ng_on2_tp / ng_on2_rate |
| `phaseC_genome_three_way/v3f_vs_v5_vs_v6_region_ng.tsv` | 105,997 (= TP 30,490 + FP 4,842 × 3 BAM) | region-level：BAM, label, region_id, ng_off, ng_on |
| `phaseD_v6_5sample/v6_cross_sample_summary.tsv` | 5 (header + 4 樣本) | 4 樣本 V6 metrics（H1437/H2009/HCC1954/HCC1937）|

**genome_summary.tsv 數值**（07_V6_validation_findings.md §4 已記載）：
```
BAM  marker_off_n  marker_off_rate  ng_on2_rate
V3F  21,997        0.9175           0.8579
V5   18,382        0.8937           0.8285
V6   23,980        0.9093           0.8285
```

**phaseD cross_sample_summary.tsv 數值**（08_phaseD_v6_cross_sample_findings.md）：
```
sample    marker_off_rate  ng_on2_rate  hp_1_1:hp_2_1_ratio
H1437     0.9923           0.9914       1.243
H2009     0.9930           0.9917       0.901
HCC1954   0.9535           0.9672       0.958
HCC1937   0.8165           0.9039       0.611
```

## 5. 跨樣本 V6 ISM 結果（phaseD_v6_5sample）

**路徑根**：`/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/phaseD_v6_5sample/`

| Sample | 子目錄 | 4 個 run (on/off × tp/fp) | 狀態 |
|--------|------|--------------------------|------|
| H1437 | `H1437/` | ✅ 完整 | n_tp=70,191 / n_fp=773 |
| H2009 | `H2009/` | ✅ 完整 | n_tp=135,359 / n_fp=1,342 |
| HCC1954 | `HCC1954/` | ✅ 完整 | n_tp=19,449 / n_fp=687 |
| HCC1937 | `HCC1937/` | ✅ 完整 | n_tp=13,910 / n_fp=2,697 |
| COLO829 | — | 🟡 deferred（truth set 0600 權限）| step4 待解後補 |

**Note**：phaseD 只有 V6（無 V3F/V5 三方對照）。Step 4 跨樣本 grid 只能在 V6 跑，trajectory 分析限 HCC1395。

## 6. Caller AF 來源（VCF）

**HCC1395 paired pileup VCF**：
- TP: `/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz`
- FP: `/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz`

**Schema**: 標準 VCF + FORMAT/AF/AD/DP 欄位

**用途**：caller_af covariate（Step 2 LR 控制）

## 7. LOH / CN / Coverage 來源（KDE-corrected master.tsv.gz）

**主檔**：`/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz`

**狀態**：
- HCC1395 paired_full：✅ 96.8% TP positions 覆蓋（KDE-corrected, commit `374fad4` + `12d9b3e`, Diploid_Coverage_Used=61）
- HCC1395 paired_full FP：**partial** — phaseC FP n=4,842 vs master FP n=627（master 沒涵蓋所有 phaseC FP）
- HCC1395 TO Phase 1：✅ KDE-corrected
- 其他 6 樣本 TO：🟡 stale-binary（Diploid_Coverage_Used 可能 = 75.0），等 Archive TO rerun

**重要欄位（schema 全欄見 phaseC significance_summary.csv 開頭）**：
```
region_id, chr, pos, NumReads, NumCpGs,
HPFineNGroups, HPFineN_HP1, HPFineN_HP1S, HPFineN_HP2, HPFineN_HP2S,
HP_Ratio, Potential_LOH, LOH_Bed_Overlap, LOH_Source, LOH_Subtype,
Coverage_Multiple, Diploid_Coverage_Used, Coverage_Category,
caller_af (via VCF lookup), label (TP/FP via truth)
```

**驗證需求**（Step 1.5 power gate 前必做）：
- `Diploid_Coverage_Used median == 61` (HCC1395 KDE baseline)
- 若 == 75.0 → stale-binary artifact，flag in findings.md，繼續但備註

## 8. SEQC2 CNV truth (HCC1395)

**路徑根**：`/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data/`

| 檔案 | 大小 | 用途 |
|------|------|------|
| `annotated_hcc1395_cnv.tsv` | 9.7 MB | per-variant TSV (paired mode)，含 SEQC2_CN + SEQC2_Zone |
| `annotated_hcc1395_to_cnv.tsv` | 15 MB | per-variant TSV (TO mode) |

**Note**：6 callers × 21 replicates × 3 technologies 共識；用於 Step 3 Z-GL (Gain+LOH) zone 校準

## 9. Phase 2 BCD C++ 模組（可重用）

| 模組 | 路徑 | 大小 | 用途 |
|------|------|------|------|
| `LohBedAnnotator` | `InterSubMod/include/core/LohBedAnnotator.hpp` + `src/core/LohBedAnnotator.cpp` | 2.8K / 4.0K | LOH.bed 載入 + 位點重疊判定 |
| `NormalBaseline` | `include/core/NormalBaseline.*` + `src/core/NormalBaseline.cpp` | — | Normal BAM baseline 層 |
| `PerCpgAsm` | `include/core/PerCpgAsm.hpp` + `src/core/PerCpgAsm.cpp` | 3.9K / 12.3K | per-CpG ASM delta 計算 |
| `SubcloneAnalyzer` | `include/core/SubcloneAnalyzer.hpp` + `src/core/SubcloneAnalyzer.cpp` | 5.8K / 11.2K | 4-group 亞克隆分層 |

**Note**：本計畫 **不修改 C++** (Out-of-scope §不修改 C++ pipeline)，但分析 reference 可用。

## 10. Driver 報告（生成本計畫的證據鏈）

| 報告 | 路徑 | 內容 |
|------|------|------|
| Self-Phasing 完整觀察整合報告 | `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md` | V5 Layer 1.5 缺陷揭露 §8.6 |
| Thread D 主軸 | `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md` | 4-bucket × Inner/Outer 框架 |
| V6 Proposal Evaluation | `InterSubMod/research/paired_priority_bug_audit/02_V6_proposal_evaluation.md` | V6 設計理由 |
| V6 Validation Findings | `InterSubMod/research/paired_priority_bug_audit/07_V6_validation_findings.md` | V6 chr19+全基因組三向驗證 |
| V6 Cross-Sample Findings | `InterSubMod/research/paired_priority_bug_audit/08_phaseD_v6_cross_sample_findings.md` | 4 樣本驗證 |
| V6 Caller F1 Verification | `InterSubMod/research/paired_priority_bug_audit/09_V6_caller_F1_verification.md` | V6 caller F1 不變 |
| Phase BCD Dual-BAM Validation | `InterSubMod/docs/experiments/validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md` | Sample ASM 97.3% sig / LOH concordance 98.4% |

## 11. Prior Art Reference（Agent D 已完成）

詳見 `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/02_prior_art_notes.md`：
- TumorLens (medRxiv 2026-03-19, doi:10.64898/2026.03.18.26348569)
- ROCIT (bioRxiv 2026-03-05, PMC12991090)
- SGZ (PLoS Comp Biol 2018, doi:10.1371/journal.pcbi.1005965)
- Wakhan (medRxiv 2025-12-15, doi:10.64898/2025.12.11.25342098)
- SAVANA (Nat Methods 2025, PMC12240814)

## 12. Memory & KB 引用

**Memory（auto-loaded MEMORY.md 內）**：
- `feedback_L2_collider_bias.md` — pooled OLS 殘差化陷阱
- `feedback_pooled_ols_residualization_trap.md` — within-group OLS 必須
- `feedback_spatial_autocorrelation_confound.md` — chr+pos 聚合 confound
- `feedback_outside_claim_must_query_kb.md` — 外部工具 claim 必先查 KB
- `project_loh_constrained_phasing_discovery.md` — Thread D NG=2 6/6 same-hap
- `project_v5_v6_tradeoff_sp123.md` — V5 vs V6 SP1/2/3 trade-off
- `project_v5_layer15_design_caveat.md` — V5 Layer 1.5 germline-absent 缺陷

**KB**（`/big8_disk/liaoyoyo2001/Knowledge/`）：
- `05_tools/longphase-to.md` (2026-04-17)
- `05_tools/longphase-s.md`
- `05_tools/variant-callers.md` (2026-04-01)
- `03_file_formats/vcf-clairs-to.md`
- `08_references/paper-index.md` (2026-04-03)

**待補建 KB**（Agent D §6 建議）：`05_tools/{wakhan,savana,tumorlens,rocit}.md`

## 13. 環境變數 & 路徑前綴

| 變數 | 值 | 用途 |
|------|---|------|
| `TMPDIR` | `/big7_disk/tmp`（必設）| 避免 /tmp 800GB 災情，見 memory `feedback_tmp_disk_full_pipeline_pitfall.md` |
| Python | `/usr/bin/python3` (system) | 跑 step1-4 script |
| Project root | `/big7_disk/liaoyoyo2001/InterSubMod` | 計畫所有相對路徑起點 |

---

## 引用本檔的規則

- Step 1-4 agent 引用資料時必用本檔列出的**確切路徑**，不憑記憶
- 若路徑與本檔不符 → 先更新本檔（記錄發現），再繼續分析
- 所有 `.md` 路徑輸出給用戶時必以 `InterSubMod/` 開頭（依 UserPromptSubmit hook 規則）
