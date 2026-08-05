<!--
建立時間: 2026-07-13 03:22 +0800
目標: 驗證 region possible clone/subclone annotation 是否可完全重現
處理範圍: chr1-22；47,377 historical layered-v2 primary regions
關聯檔案: InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_reproducibility/receipt.json
-->

# Region possible clone/subclone annotation 重現性 receipt

## 結論

PASS — 兩次獨立執行各自為 32/32 checks PASS；四個 canonical data/QA outputs 皆 byte-identical。

## 輸入與命令

輸入與完整 SHA-256 列於兩次 run 的 `provenance.json`。執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260713_crosssample_structure_bulk_sampling_adjustment/scripts/build_region_possible_clone_annotations.py \
  --output-dir research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1
python3 research/20260713_crosssample_structure_bulk_sampling_adjustment/scripts/build_region_possible_clone_annotations.py \
  --output-dir research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1_rerun
```

## 實際輸出

```text
regions=47377; hcc_fixed_pair_regions=5720; checks=32/32 PASS
region_possible_clone_annotations.tsv.gz: ffcdb7fb5effa9fceb1e18622567de198610b36327cce41fe2144cd3746afe50
sample_possible_clone_summary.tsv: d41f7e246047ac8c7317ca26d81a008387044ebbe46ac0af1343793cf7ce333b
hcc1395_pair_region_possible_clone_comparison.tsv.gz: a5fbb43a9a03ce244abd41823b2703b471e2ac26f664f4a5d2fdaac5fb864943
checks.tsv: 9f91b7479f5537cf57bbfde344c27198018314202a5f86f85c3e51ae1f337921
```

此 receipt 只證明 deterministic annotation pipeline；不把 PyClone cluster、CP 或 VAF 最可能樹升格為真實 clone／演化樹。空 subclonal union 的 Jaccard 會輸出 `NA + both_absent`，single-cluster partition 會輸出 `external_partition_informative=False`，避免 vacuous perfect agreement。
