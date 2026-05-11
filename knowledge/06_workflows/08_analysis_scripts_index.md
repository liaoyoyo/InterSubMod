---
id: ism-kb-06-workflows-analysis-scripts-index
name: "Analysis Scripts 索引（151 個）"
description: "`scripts/analysis/` 151 個 Python 腳本按前綴分 6 大類索引；15 個高頻腳本詳細表；發現新腳本的流程。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "scripts/analysis 實際檔案列表於 HEAD"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-05-data-formats-master-dataset-schema
  - ism-kb-03-pipelines-f1-baseline-canonical
  - ism-kb-07-derived-features-index
  - ism-kb-06-workflows-pptx-and-weekly-report
tags: [workflow, scripts, analysis, index, python]
canonical_paths: [06_workflows/08_analysis_scripts_index.md]
alias_paths: []
---

# Analysis Scripts 索引（151 個）

- 一句結論：`scripts/analysis/*.py` 151 個；按前綴分 6 大類；本頁列 15 個高頻腳本詳細表 + 其他按類別
- 適用對象：找可用分析工具、寫新腳本前檢查重複
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/*.py | wc -l
  # 預期 151
  ```

---

## 📊 分類總覽

| 前綴 | 數量 | 用途 | 代表腳本 |
|------|-----|------|---------|
| `build_*` | 79 | **資料聚合/矩陣構建**（從 canonical run 生成 master dataset / 分析表） | `build_loh_round1_cross_sample_audit.py` → master_dataset |
| `run_*` | 7 | **執行類**（pipeline wrapper / 跨樣本 benchmark） | `run_batch_vcf_analysis.sh`（上游） |
| `analyze_*` | 14 | **分析類**（TP/FP 特徵、AF gradient、cross-sample） | `analyze_to_tp_fp_characterization.py` |
| `compare_*` / `summarize_*` / `export_*` | 11 | **對比/輸出**（F1、rescue strategies、candidate pool） | `compare_benchmark_f1.py` |
| `verify_*` / `validate_*` / `check_*` / `evaluate_*` | 9 | **驗證類** | `methylation_cn_within_group_validation.py` |
| 其他專題 | ~31 | zone / purity / seqc2 / pon / methylation / haplotag / PPTX 等 | `collect_baseline_metrics.py`、`build_weekly_report_pptx.py` |

> **Note**：以 `ls scripts/analysis/*.py | wc -l` 為唯一權威（當前 151）；分類統計每次大改動後重算。

---

## ⭐ 高頻腳本詳細表（15 個）

### 1. `build_loh_round1_cross_sample_audit.py`
- **用途**：聚合 7 樣本 × 2 modes × 748K regions 成 master_dataset (`all_region_rows.tsv.gz`, 116 欄)
- **輸入**：`output/canonical/*/{paired_full,paired_pileup,TO}/*/significance_summary.csv`
- **輸出**：`output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`
- **對應 KB**：[05_data_formats/02_master_dataset_schema.md](../05_data_formats/02_master_dataset_schema.md)

### 2. `compare_benchmark_f1.py`
- **用途**：跨 run 聚合 F1（precision/recall/f1/tp/fp/truth_total）
- **輸入**：canonical run 目錄
- **輸出**：`benchmark_comparison.tsv`
- **對應 KB**：[03_pipelines/05_f1_baseline_canonical.md](../03_pipelines/05_f1_baseline_canonical.md)

### 3. `collect_baseline_metrics.py`
- **用途**：PDD Step 1/5 的基線 F1 記錄
- **對應 KB**：[07_cpp_change_pdd.md](07_cpp_change_pdd.md) Step 1 + Step 5c

### 4. `analyze_to_tp_fp_characterization.py`
- **用途**：TO 模式 TP/FP 特徵分析
- **對應 KB**：[03_pipelines/03_tumor_only.md](../03_pipelines/03_tumor_only.md)

### 5. `analyze_cross_sample_ism_af_gradient.py`
- **用途**：跨樣本 AF gradient 分析
- **對應 KB**：[07_derived_features/04_loh_af_methylation.md](../07_derived_features/04_loh_af_methylation.md)

### 6. `audit_provenance_chain.py`
- **用途**：研究證據鏈追溯審計（對應 `/provenance-tier-audit` skill）
- **對應 KB**：[10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)

### 7. `check_md_links.py`
- **用途**：Markdown 連結有效性檢查（docs/ 用）
- **類似**：`knowledge/scripts/check_canonical_paths.py`（但後者僅檢 frontmatter）

### 8. `extract_label_first_metrics.py`
- **用途**：從 canonical run 抽 `label_first_metrics.tsv`（供 methodology-audit 用）
- **對應 KB**：[00_governance/09_confirmation_protocol.md](../00_governance/09_confirmation_protocol.md)（methodology-audit 前置數據）

### 9. `build_germline_fp_G{1-7}_*.py`（7 個腳本群）
- **用途**：G1-G7 germline FP 特徵研究（已 NEGATIVE）
- **對應 KB**：[09_conclusions/03_concluded_negative.md](../09_conclusions/03_concluded_negative.md) N8 TO Germline FP NO-GO

### 10. `build_allele_loh_auc_f1_analysis.py`
- **用途**：LOH × Allele × AUC × F1 綜合分析
- **對應 KB**：[07_derived_features/04_loh_af_methylation.md](../07_derived_features/04_loh_af_methylation.md)

### 11. `methylation_cn_within_group_validation.py`
- **用途**：CN-aware within-group OLS 驗證（避免 collider bias）
- **對應 KB**：[09_conclusions/03_concluded_negative.md](../09_conclusions/03_concluded_negative.md) N6 O12 L2 collider bias

### 12. `build_weekly_report_pptx.py`（**週報 PPTX 主腳本**）
- **用途**：每週 PI 週報 PPTX 生成主流程；對應 weekly-report skill Phase 4
- **相關**：`build_weekly_report_pptx_oral_v2.py`（口試版）、`build_loh_weekly_pptx.py`（LOH 主題）、`build_ai_f1_research_pptx.py`
- **對應 KB**：[09_pptx_and_weekly_report.md](09_pptx_and_weekly_report.md)
- **注意**：各 `docs/presentations/*/` 下有 instance-local `build_pptx.py`（模板複製的一次性檔），**非**全域腳本

### 12b. `generate_pi_report_figures_self_phasing.py`
- **用途**：PI 報告用 Self-Phasing 視覺化（Fig 12 類）
- **對應 KB**：[06_workflows/09_pptx_and_weekly_report.md](09_pptx_and_weekly_report.md)

### 13. `haplotag_qc.py`
- **用途**：HP tag 品質 QC（AMB%、clean blocks）
- **對應 KB**：[09_conclusions/01_positive_findings.md](../09_conclusions/01_positive_findings.md) P4 V5 Somatic Fallback

### 14. `gc_correction_to_validation.py`
- **用途**：TO mode GC 校正驗證
- **對應 KB**：[07_derived_features/02_coverage_multiple.md](../07_derived_features/02_coverage_multiple.md)

### 15. `analyze_longphase_rescue_with_methylation.py`
- **用途**：LongPhase rescue × 甲基化整合分析
- **對應 KB**：[09_conclusions/03_concluded_negative.md](../09_conclusions/03_concluded_negative.md) N17 TO Feature Deep Study

---

## 🔍 發現新腳本的流程

### 1. 按關鍵字搜尋
```bash
# 例：找 benchmark 相關
ls /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/*.py | grep -i benchmark

# 例：找 LOH 相關
ls /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/*.py | grep -i loh
```

### 2. 讀腳本註解（前 50 行通常有用途）
```bash
head -50 /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/<script>.py
```

### 3. 看 argparse 定義
```bash
grep -A 3 "add_argument" /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/<script>.py | head -20
```

### 4. 跑 `--help` 若有
```bash
python3 /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/<script>.py --help
```

---

## 📁 次目錄（非頂層）

`scripts/analysis/legacy/`：舊版腳本（不建議使用，僅保留歷史）

---

## 🆕 貢獻新腳本指引

### 命名慣例
- **動詞_名詞** pattern：`analyze_to_tp_fp_characterization.py`
- **前綴選擇**（照上方分類表）：
  - 聚合/矩陣構建 → `build_*`
  - 分析/描述性統計 → `analyze_*`
  - 執行 pipeline → `run_*`
  - 比較多個輸出 → `compare_*`
  - 驗證/檢核 → `verify_*` / `validate_*` / `check_*`
- **專題性高** → 用主題名（不強制加前綴；例：`haplotag_qc.py`）

### 必備元素
1. **Docstring**：前 50 行說明用途、輸入、輸出、對應 KB 概念
2. **argparse**：所有參數明確定義
3. **對應 evidence_ledger cycle**（若為研究類）：在 docstring 記錄 cycle_id

### 整合到 KB（若為高頻腳本）
- 加入本頁「高頻腳本詳細表」
- 若對應既有 KB 主題 → 在該主題頁加「典型分析」章節引用腳本

---

## ⚠️ 常見陷阱

- ❌ **寫新腳本前不查重**：已有 81 個 `build_*`，可能已有類似功能
- ❌ **不寫 docstring**：半年後沒人記得做什麼用
- ❌ **硬編碼路徑**：應用 argparse 或環境變數
- ❌ **忘記 evidence_ledger 紀錄**：研究類腳本產出必記入 ledger

---

## 相關

- Master dataset（最常用輸入）：[../05_data_formats/02_master_dataset_schema.md](../05_data_formats/02_master_dataset_schema.md)
- F1 SoT（benchmark 相關腳本）：[../03_pipelines/05_f1_baseline_canonical.md](../03_pipelines/05_f1_baseline_canonical.md)
- PDD 6 steps（腳本整合 C++ workflow）：[07_cpp_change_pdd.md](07_cpp_change_pdd.md)
- 衍生特徵 index（feature 計算腳本對應）：[../07_derived_features/00_index.md](../07_derived_features/00_index.md)
- Evidence ledger：[../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)
