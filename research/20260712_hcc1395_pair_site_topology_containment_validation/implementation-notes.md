<!--
建立時間: 2026-07-12T04:54:04+08:00
目標: 紀錄 HCC1395 逐位點 topology containment 實作決策、偏離、折衷與未決事項
處理範圍: chr1-22；5,720 exact-coordinate complete-both regions
關聯檔案:
  - InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/pre-decision-audit.md
-->

# Implementation Notes：HCC1395 site-aware topology containment

## 設計決定

- [決策] 2026-07-12：Task B；以 5,720 exact-coordinate complete-both regions 為固定母群。
- [決策] 2026-07-12：read full candidate set 與 VAF exact-score top set分開；兩者都保留 candidate uncertainty。
- [決策] 2026-07-12：子結構只由 genomic shared-site induced projection 定義；禁止任意 graph subtree existential match作 headline。
- [決策] 2026-07-12：fixed region outcome 每個 evidence layer 預先以 identity 或 whole-region single swap 擇一；site-pair ledger 固定 Read-selected mapping 給 Read/VAF；HP-count mismatch fail-closed。
- [決策] 2026-07-12：shared-site recurrence 不壓平；primary endpoint fail-closed。

## 偏離

- [偏離] Canonical portable verifier 在 classic-scrollbar Chromium 偵測到 reader `100vw` header 溢位；沿用本 repo 已驗證的 `apply_portable_header_compat.py`，只將 sticky header 改為 `width/max-width:100%` 並重跑完整 verifier。

## 折衷

- [折衷] 現階段使用 historical layered-v2 engineering snapshot，因 clean-v3 producer雖 7/7 PASS，release receipt仍 `E_SIDECAR_VALIDATION`，canonical/sensitivity clean roots尚未可用。
- [折衷] VAF exact-sum照 frozen official 定義重現，但因 comparison count varies，另做 constant-count／normalized sensitivity，不升格為獨立 confirmation。

## 未決

- [未決] Depth-matched read downsampling與 within-dataset split-half需要重跑 solver／read抽樣；若本輪未完成，必列為後續 biological/technical validation缺口。
- [未決] Cross-cell-line negative control缺相同 genomic sSNV identity；若採 shape-only control，必明示它只測 algorithmic prior。
- [未決] 目前 site-label null 直接檢驗 projected candidate-set relation，不直接檢驗 fixed `strict full exact + true induced` endpoint；後者仍需 locked direct null。

## 執行結果

### Core topology comparison

- 輸入：`InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv`、historical layered-v2 HCC1395/HCC1395_DORADO、VAF units/candidates。
- 命令：`python3 InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/scripts/build_site_topology_containment.py ...`
- 輸出：`InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data/`
- 實際片段：`regions=5720 units=15545 trees=97663 checks=18/18 PASS`。

### Postprocess / site ledger

- 命令：`python3 InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/scripts/summarize_site_topology_results.py ...`
- 實際片段：`regions=5720 evidence_pairs=8096 shared_alleles=15713/15713 checks=23/23 PASS`。
- Read strict exact＋真 induced：1,599/5,720；VAF：1,790/5,720。

### Site-label permutation null

- 命令：`python3 InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/scripts/run_site_label_permutation_null.py --permutations 5000 --seed 20260712 --bootstrap-replicates 5000 ...`
- 實際片段：`status=PASS regions=5720 checks=25/25 PASS`。
- Swap-tolerant exact excess：Read +8.07 pp；VAF +23.86 pp；兩者 empirical p=1/5001、block CI>0、LOCO 22/22。

### HTML artifact

- Builder：`InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/scripts/build_site_topology_report_artifact.py`。
- 輸出：`InterSubMod/docs/reports/in_progress/2026/07/20260712_HCC1395逐位點拓撲子結構跨來源驗證_01/`。
- Builder 片段：`regions=5720 report_checks=53 charts=5 tables=12 null=available`。
- Browser QA：41 blocks、5 charts、12 tables、3 metrics；1440/390 viewports 與 source dialog PASS。

## Provenance Footer

- Commit hash：`6067568637088838a9f518955e41d222f057e4f1`
- Build time：2026-07-12T04:54:04+08:00
- Skill：`/implementation-notes` v0.1
