<!--
建立時間: 2026-07-12T09:33:06+08:00
目標: 重建 HCC1395 兩技術來源固定 5,720-region 的位點、區域與拓撲精確分層
處理範圍: GRCh38 chr1-22; exact-coordinate complete-both regions; historical layered-v2 context
關聯檔案: InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/scripts/build_regional_precision_validation.py
-->

# Regional precision validation execution record

## 結論

**PASS WITH CLAIM CEILING** — 32/32 conservation checks PASS。本輪支持固定 population 的技術再現性，不支持唯一真實演化樹 accuracy。

## 輸入路徑

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv` (sha256 `85dc631e211757d5d5fdf97d2132ace15e23961bc41c46e4bd0113ee3a171f16`)
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_regions.tsv` (sha256 `b185a98518e2d4cdb7e6f3163463c06122c5d6febb1d29895f3c2272a6ee6ce2`)
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data/hcc1395_site_topology_pair_outcomes.tsv` (sha256 `febc072ec51c38681f2d9f29db572be5b46e127ea23a77855aa5d6d36b8a21bb`)
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data/hcc1395_site_topology_regions.tsv` (sha256 `7701ea9343f2f073467d85f4826d482de126007c31b352d5ac3c80d1b404a773`)
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data/hcc1395_site_allele_identity.tsv.gz` (sha256 `00a00788033a91b2e4b76d54fb232cb79b8ae87d72ab473dda9552d2dbea249f`)
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/integrated_topology_context_v1.json` (sha256 `d3bfd43a96e72942ec3798a1a458491cf533f49d0b85f2c8f3fd4519c04d354a`)

## 執行命令

```bash
/bip7_disk/liaoyoyo2001/miniconda3/bin/python /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/scripts/build_regional_precision_validation.py --coarse-dir /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data --site-dir /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data --integrated-context /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/integrated_topology_context_v1.json --output-dir /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1
```

## 輸出路徑

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/regional_precision_outcomes.tsv.gz`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/regional_precision_metrics.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/k_region_counts.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/site_set_relations_by_k.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/allele_identity_by_k.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/coarse_confusion_by_k.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/availability.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/input_schema_inventory.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/validation_checks.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/regional_precision_summary.json`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/regional_precision_validation_v1/execution_record.md`

## 實際輸出片段

```text
k=1: 0
k=2: 3,121
k=3: 1,506
k>=4: 1,093
site-set equal: 5,534/5,720
shared allele identity: 15,713/15,713
coarse 5-state: 3,969/5,720
read compatible: 4,036/4,038
read strict+induced: 1,599/5,720 fixed; /4,038 evaluable
VAF strict+induced: 1,790/5,720 fixed; /3,860 evaluable
structure-first+VAF shape swap: 3,667/5,720 fixed; /5,168 evaluable
VAF-unique exact forest swap: 949/5,720 fixed; /2,543 evaluable
validation: 32/32 PASS
```

## 分母與可用性說明

- `global_fixed_denominator=5,720` 是固定 population；不因 endpoint 可否評估而改變。
- `stratum_fixed_denominator` 是 k 分層內全部 regions；`evaluable_denominator` 才是 endpoint-specific gate 後分母。
- k 以 `caller_shared_k` 定義。k=1 在這 5,720 regions 中是資料支持的 0，不是缺值。
- VAF/read endpoints 復用同批 reads；shape 移除 mutation labels/direction；二者都不是 orthogonal truth。
