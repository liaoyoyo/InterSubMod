<!--
建立時間: 2026-07-17T11:58:00+08:00
目標: 保存 ClairS 到 LongPhase-S、tree input、read tag join 與 BAM/FIFO 的獨立現檔重算
處理範圍: 7 datasets / 638,259 raw-all records / 469,849 frozen autosomal dataset-sites / 38,345,639 site-read observations
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/latest_tree_input_contract_audit.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/frozen_input_immutability.post_screen_recovery_v1.json
-->

# ClairS / LongPhase-S / tag / tree / BAM upstream 獨立稽核

## Verdict

- Agent：`019f6e24-18c6-7233-a15a-79136de803e2`（Poincare）。
- 模式：唯讀、現檔重算；未修改 repo。
- 結論：**GO WITH CLAIM CEILING；無 High upstream blocker。** 整體科學release仍須等待downstream完成。

## 位點守恆

| 階段 | records | 守恆／排除口徑 |
|---|---:|---|
| Normalized ClairS | 638,259 | 568,080 PASS + 70,179 LowQual |
| LongPhase-S all-output | 638,259 | record key missing/extra/mismatch=0 |
| LongPhase-S recalibrated PASS | 582,820 | 32,184 LowQual->PASS；17,444 PASS->LowQual；final LowQual=55,439 |
| Frozen chr1-22 | 469,849 | 排除chrX=112,877、chrY=94；other=0；non-biallelic/non-sSNV=0 |

等式：`568,080 + 32,184 - 17,444 = 582,820`；`582,820 - 112,877 - 94 = 469,849`。
`225,268 max_snv_excluded + 50,432 positional_singleton + 194,149 retained = 469,849`是frozen universe內的
下游tree branch分派，不是再從469,849移除。

本輪使用receipt-bound locally patched LongPhase-S，不可描述成未修改stock binary。Patch只讓無eligible read的
LowQual site保留並繼續；兩個實際warning都在chrX，不進autosomal 469,849 scope。

## Tree input

7/7 tree VCF均為`20260711_longphase_s_raw_all_production_sidecars_v2`同一producer run的recalibrated
`FILTER=PASS` VCF；現檔SHA-256全MATCH、nonPASS=0、non-biallelic-sSNV=0：

- HCC1395 113,061 -> 79,687
- HCC1395_DORADO 113,145 -> 79,739
- COLO829 39,458 -> 37,788
- H1437 79,655 -> 77,080
- H2009 161,595 -> 154,465
- HCC1937 52,115 -> 18,690
- HCC1954 23,791 -> 22,400

`same-run`只指LongPhase-S producer VCF/sidecar；正式v10 screen是source-locked first-six加HCC1954 exact-equivalent
recovery merge，不是單一monolithic screen execution。

## Read tag join

現檔逐列重算：469,849 rows、469,849 unique dataset-site keys、38,345,639 site-read observations、
38,345,639 authoritative sidecar exact joins、projection multimatch=0、zero-read sites=0；PS present=24,694,402，
source HP replaced=33,693,248。

限制：38,345,639不是globally unique read數，也不是每列皆有PS。Sidecars另有12,681,029個HP/PS相同的exact
duplicate rows，依契約collapse；multimatch=0只表示collapse後沒有projection對到多個不同full identities。

## BAM / FIFO

7/7 producer receipts與現檔類型皆顯示`transport=named_fifo`、`persisted_bam=false`、`regular_bam_count=0`；
FIFO節點保留、size=0，但tagged BAM payload未持久化。raw/frozen inputs post-screen為77/77 match，既有canonical
tagged BAM audit為7/7 `all_match=true`，未觀察到原始或既有tagged BAM遭覆寫。

限制：大型BAM未做full-content SHA-256；assurance為path/inode/size/mtime/ctime、首中末sampled chunks及BAI full hash。
因此只能發布「在`storage_identity_v1` assurance下未觀察到覆寫」，不能發布「逐位元證明完全不變」。

## 可發布 upstream 敘述

在frozen 7-dataset / 6-biological-sample / chr1-22契約內，638,259個normalized ClairS records無損進入
receipt-bound locally patched LongPhase-S；582,820為recalibrated PASS，其中469,849個autosomal dataset-sites
形成frozen scope。7/7 tree inputs與HP/PS sidecars來自同一producer run；38,345,639個site-read observations全數
完成authoritative sidecar exact join，distinct-identity projection multimatch=0。Tagged-BAM bytes只經FIFO串流，
未持久化regular BAM；在既定sampled-identity assurance下未觀察到原始或既有tagged BAM遭覆寫。

此敘述只證明upstream provenance與execution integrity，不證明somatic truth、科學效果、globally unique reads、
stock LongPhase-S等價性或完整end-to-end scientific release。
