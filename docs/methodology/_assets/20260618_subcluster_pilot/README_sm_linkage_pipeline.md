<!--
建立時間: 2026-06-26
類型: pipeline manifest — 全基因組 sSNV 單分子連鎖 + 演化重建 reproducible record
狀態: in_progress（HCC1395 ⭐3, Tier-R）
build_branch: feat/summary-nreadsvalid
-->

# 全基因組 sSNV 單分子連鎖 + 演化重建 — pipeline manifest（可重現 + 可查詢）

> **單位**：最大可關聯區域（≤50kb somatic-sSNV chain）。從 2-locus pairwise → multi-locus 組合 → per-region 樹。
> **樣本**：HCC1395 ⭐3 單樣本，Tier-R（same-read ≤50kb；same-PS deferred）。所有數字落檔，每個可 grep 重算。
> **對抗稽核**：`sm_adversarial_audit.md`（NEEDS_WORK → 已套 5 修正：per-sSNV / 偽影量化 / γ-reframe / Fisher 更正 / Tier-R 標）。

## 資料來源（verified 路徑）
| 角色 | 路徑 |
|---|---|
| sSNV 宇宙 | `/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_{tp,fp}_{chr}.vcf.gz` |
| tumor BAM | `/big8_disk/.../data/bam/HCC1395/tumor.bam` (HCC1395_Tmode_tagged_ClairS_pileup) |
| normal BAM | `.../normal.bam` (HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged) |
| methyl_PASS resolver | `.../pileup/HCC1395_methyl_PASS_fixed.vcf.gz` |
| SEQC2 truth | `/big8_disk/data/HCC1395/SEQC2/{high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz, sSNV.MSDUKT.superSet.v1.2.vcf.gz, High-Confidence_Regions_v1.2.bed}` |
| CN/LOH | `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed` |

## Pipeline（script → input → output；依序）
| # | script | input | output | 關鍵數字所在 |
|---|--------|-------|--------|-------------|
| S0-S2 | `sm_linkage_genomewide.py {chroms} {out}` (×5 parallel) | TP∪FP VCF + tumor/normal BAM | `sm_linkage_gw_part_{1-5}.json` | universe / pairs / census |
| merge | `sm_merge_parts.py` | part_{1-5} | **`sm_linkage_genomewide.json`** | universe 35,332 / pairs 53,094 / census(35,332) |
| S3 | `sm_evolution_build.py` | sm_linkage_genomewide.json | `sm_evolution.json` | components / nested / sibling / cycles |
| S4 | `sm_completeness_ledger.py` | 同上 | `sm_completeness_ledger.json` | buckets(linked/underpowered/isolated) sum-check |
| S5 null | `sm_verify_null.py` | 同上 | `sm_verify_null.json` | Fisher per_rel + diffHP negctrl |
| S5 SEQC2 | `sm_seqc2_overlay.py` | 同上 + SEQC2 | `sm_seqc2_overlay.json` | truth_class + gamma_like |
| 修正+summary | `sm_summary.py` | 全部 + CN | **`sm_summary.json`** | **headline SoT（per-sSNV / 偽影 / Fisher-vs-MC / CN 分層）** |
| 配置census | `sm_configuration_census.py` | sm_linkage_genomewide.json + CN | `sm_configuration_census.json` | 2×2 cell-pattern 計數 |
| pair lists | `sm_pair_lists.py` | 同上 + CN | `lists/{config}_{HP}.tsv` (+CLEAN) + `per_sSNV_census.tsv` | 各 config 對列表 |
| multi-locus | `sm_multilocus_combinations.py {chroms} {out}` (×5) | TP∪FP + tumor BAM + census | `ml_part_{1-5}.json` | per-region populations |
| **per-region 整合** | `sm_region_integration.py` | sm_linkage_genomewide.json + ml_part_* + CN | **`sm_region_integration.json`** + `lists/regions.tsv` | **每區域樹形 + 比例（最終 SoT）** |
| 對抗稽核 | `evaluator` agent (read-only) | 全部 | `sm_adversarial_audit.md` | NEEDS_WORK + 5 findings |

## 重現（reproduce）
```bash
cd .../scripts
# 1) genome-wide linkage (5 parallel, ~40min) → merge
for p in "chr1,chr6,chr11,chr16,chr21:1" "chr2,chr7,chr12,chr17,chr22:2" "chr3,chr8,chr13,chr18:3" "chr4,chr9,chr14,chr19:4" "chr5,chr10,chr15,chr20:5"; do
  python3 sm_linkage_genomewide.py ${p%:*} ../sm_linkage_gw_part_${p#*:}.json & done; wait
python3 sm_merge_parts.py
# 2) downstream (秒級)
python3 sm_evolution_build.py; python3 sm_completeness_ledger.py
python3 sm_verify_null.py; python3 sm_seqc2_overlay.py; python3 sm_summary.py
python3 sm_configuration_census.py; python3 sm_pair_lists.py
# 3) multi-locus (5 parallel, ~15min) → region integration
for p in ...; do python3 sm_multilocus_combinations.py ${p%:*} ../ml_part_${p#*:}.json & done; wait
python3 sm_region_integration.py
```

## 驗證任一數字（verify）
- **核心問句**：報告任一數字 → 「在哪個檔案 grep 到？」grep 不到 = 捏造。
- headline 數字 SoT = `sm_summary.json`（小檔，可獨立 recompute）+ `sm_region_integration.json`。
- 逐位點查：`lists/per_sSNV_census.tsv`（每 sSNV：src/normal/somatic/vaf/links/CN）。
- 逐對查：`lists/{config}_{HP}.tsv`。逐區域查：`lists/regions.tsv`。
- sum-check：census 數 == universe；ledger buckets 加總 == universe（`sm_completeness_ledger.json.bucket_sum_check`）。

## 🔴 限制（誠實）
- ⭐3 單樣本；Tier-R only（same-PS deferred）；**~69% linked-somatic 在 CN-gain（混淆）**，乾淨集 = LOH+neutral。
- 偽影：dense-uniform-cluster（chr8/9/14）未清（缺 mappability track）；Fisher-sig 不分 subclone/allelic（HP-gate 才分）。
- regional（≤read-span）非 genome-wide tree；分子證據非 single-cell confirmation。
- 大 JSON（sm_linkage_genomewide.json 等）不入版控（可由本 manifest 重生）；腳本 + 小 summary + list TSV + 文件入版控。
