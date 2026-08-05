<!--
建立時間: 2026-07-13T03:50:29.920937+08:00
目標: 產生 region-level possible clone/subclone conditional annotations
處理範圍: chr1-22；historical layered-v2；5 CN-ready samples + 2 fail-closed samples
關聯檔案: InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/scripts/build_region_possible_clone_annotations.py
-->

# Region possible clone/subclone annotation execution record

## 結論

PASS WITH CLAIM CEILING — 32/32 checks PASS。`C_read_groups` 與 `external_cluster_count` 分離；possible state 不是 clone truth。

## 輸入路徑

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_regions.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/region_cluster_patterns.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/hcc1395_exact_complete_pair_gene_drug_flags.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/raw_vaf_validation_v1/data/raw_vaf_records.tsv.gz`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/VALIDATION_SCOPE.md`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/hcc1395_pair_primary_separate_HCC1395_main/results.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/hcc1395_pair_primary_separate_HCC1395_DORADO_main/results.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/H1437_individual_main/results.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/H2009_individual_main/results.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/HCC1954_individual_main/results.tsv.gz`

## 執行命令

```bash
python3 /big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/scripts/build_region_possible_clone_annotations.py --output-dir /big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun
```

## 輸出路徑

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun/region_possible_clone_annotations.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun/sample_possible_clone_summary.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun/hcc1395_pair_region_possible_clone_comparison.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun/checks.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun/summary.json`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun/provenance.json`

## 實際輸出片段

```text
regions=47377
hcc_fixed_pair_regions=5720
checks=32/32 PASS
pair_state_counts={'DORADO_only_subclone_signal': 65, 'HCC1395_only_subclone_signal': 66, 'both_subclone_signal_different_state': 11, 'not_evaluable_missing_external_cluster': 1282, 'same_possible_state': 4296}
```

## Claim ceiling

PyClone-VI CP／cluster只提供條件式 possible-state annotation；不能確認 clone 數、subclone identity、祖先方向或真實演化樹。High-assignment subset 可能因 subclonal union 退化而產生 selection-induced perfect agreement。
