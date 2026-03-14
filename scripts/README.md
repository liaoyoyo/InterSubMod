# scripts 目錄說明

## 主要入口

1. `pipeline/run_benchmark.sh`：建議主入口（longphase-s -> InterSubMod -> filter analysis）
2. `pipeline/run_multi_sample.sh`：多樣本批次入口
3. `analysis/run_purity_and_standard_verification.sh`：純度流程入口
4. `run_vcf_all_snv.sh`：舊版單次分析入口（legacy）
5. `run_batch_vcf_analysis.sh`：舊版 TP/FP 批次入口（legacy）

## 子目錄

1. `analysis/`：分析與報告輔助腳本
   - `summarize_three_way_comparison.py`：從三路對照輸出彙整 `three_way_comparison_summary.tsv` 與 `three_way_pairwise_delta.tsv`
   - `extract_label_first_metrics.py`：將 `significance_summary.csv` 轉成 `label_first_metrics.tsv`
   - `haplotag_qc.py`：統計 tagged BAM 的 HP 分派率，輸出 `haplotag_qc.tsv`
   - `build_purity_aware_tables.py`：由 `purity_status.tsv` 與 `metrics.json` 產出 `purity_qc.tsv`、`purity_rule_eval.tsv`
   - `run_pure_research_round.sh`：研究 round 入口，整理 pure-sample run 為標準化 bundle
   - `compare_benchmark_f1.py`：統一輸出 `benchmark_comparison.tsv` / `.md`
   - `validate_method_design.py`：輸出 `method_design_validation.tsv` 與 `label_cluster_agreement.tsv`
   - `methodology_sensitivity_matrix.py`：彙整多個 round sample bundle 成 `methodology_sensitivity.tsv`
   - `export_region_diagnostics.py`：由既有 region 輸出生成 `region_summary.tsv` 與 heatmaps
   - `build_round_dashboard.py`：生成 `round_summary.md`、plots 與 `failure_diagnosis.md`
   - `analyze_candidate_rules.py`：跨樣本比較候選規則的 `TP/FP/F1 delta`，支援 `--skip-longphase` run 的 TP/FP VCF fallback
   - `analyze_shift_patterns.py`：彙整 `Weak->Strong / Noise->Strong / Weak->Subclone / Noise->Subclone` 的 TP/FP 分層摘要
   - `region_samtools_snapshot.sh`：對指定 region 輸出 `depth.tsv`、`mpileup.txt`、`reads.sam`、`hp_tag_counts.tsv`
   - `refine_strong_labels.py`：將 `Strong` 細分為 suspect/trusted 子群，量化 refined rule 對 `TP/FP/F1` 的影響
   - `summarize_samtools_snapshots.py`：彙整多個 snapshot 目錄，輸出 `snapshot_summary.tsv` / `.md`
   - `analyze_longphase_rescue.py`：比較 caller TP/FP 與 LongPhase kept TP/FP，分析 `LongPhase-S/TO` 誤刪 TP 的 caller rescue 規則
   - `analyze_longphase_rescue_with_methylation.py`：在 `LongPhase rescue` 問題上加入 InterSubMod 特徵，評估 caller 與 methylation 組合救回 TP 的效果
2. `pipeline/`：多階段 pipeline 腳本與步驟
3. `hooks/`：流程守門與檢查腳本
4. `mcp/`：MCP 服務腳本（pubmed/biorxiv/ensembl/github）

## Output 搬移工具

1. `analysis/migrate_output_to_bip8.sh`：將 `/big8_disk/.../InterSubMod_runs/output` 舊輸出安全搬移到 `/bip8_disk/.../InterSubMod_out/output`，並輸出 manifest（不刪檔）。

## 注意事項

1. 腳本路徑應以 repo 根目錄為相對基準。
2. 新增腳本請附最小 `--help` 或註解說明。
3. 牽涉固定資料路徑時，請於腳本頂部標註可覆寫參數。
