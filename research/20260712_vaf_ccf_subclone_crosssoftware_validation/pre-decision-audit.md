<!--
建立時間: 2026-07-12
目標: 在執行 7-dataset VAF/CCF 與外部 subclone 工具驗證前，鎖定可否證假設、資料 gate 與 claim ceiling
處理範圍: chr1-22；HCC1395/HCC1395_DORADO 技術來源比較；其餘五個 cell line 個別 characterization
關聯檔案:
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/docs/method_comparison/20260619_subclone_analysis_interpretation_full_framework_01.md
  - InterSubMod/docs/method_comparison/20260630_ism_positioning_vs_prior_work_01.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/input_manifest.json
狀態: verdict_PROBE
-->

# VAF／CCF／外部 subclone 工具跨來源驗證：Pre-decision audit

> **Task type B — comprehensive validation；scope=chr1–22、7 datasets。** 這不是 subset pilot。整體科學 verdict 為 **PROBE**：可完成 raw-VAF 技術再現性與具可靠 CN/purity 樣本的外部工具 characterization；在沒有 orthogonal truth、source-specific allele-CN 與 within-source noise ceiling 前，不得宣稱 clone 演化樹已被證明相同或方法 accuracy 已被驗證。

## §0 Cynefin front-gate

- Domain：**Complex**。VAF extraction 是可預測的 complicated workflow，但「相同 VAF/cluster 是否代表相同 clone evolution」受 CN、purity、multiplicity、caller selection 與 tree non-identifiability 共同影響。
- 行動：safe-to-fail probes，所有結果保留分母、missingness、uncertainty 與 abstention；禁止由 aggregate distribution 直接跳到 tree truth。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---:|---|
| 7/7 canonical `somatic_pass.vcf.gz` 可讀，chr1–22 均有 DP/AF/AD | ✓ | L1 | `InterSubMod/.../input_manifest.json` 所列 7 VCF；本輪 `bcftools` 唯讀盤點 |
| ClairS `AD[0]` 在 455,210 loci 為 0；ref count 必須用 `DP-AD[1]` | ✓ | L1 | 本輪 7-VCF schema/count audit |
| HCC pair 全 caller exact-allele shared=76,848，Jaccard=0.9314 | ✓ | L1 | 本輪 shared-locus raw-VAF probe |
| HCC pair raw VAF：Spearman=0.8997、CCC=0.9291、MAE=0.0627 | ✓ | L1 | 本輪 shared-locus raw-VAF probe |
| HCC pair source-specific allele-CN 均可獨立取得 | ✗ | L1 | DORADO 目前只能共用 HCC1395 CN，屬 sensitivity，不是正式獨立 CN |
| COLO829/HCC1937 有可信 allele-CN | ✗ | L1 | SAVANA mis-fit；repo 明訂 unavailable |
| 有 single-cell／colony／synthetic tree truth 可驗 ancestry | ✗ | L1 | `InterSubMod/docs/CURRENT_FOCUS.md` 與既有 method-comparison documents |
| PyClone-VI/MOBSTER/Pairtree 已安裝並通過 smoke test | □ | L1 | 本輪環境盤點：目前未安裝 |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20/20 | VAF→CN/purity/multiplicity-aware CCF→clustering 是成熟框架 |
| 觀察支撐 | 10/20 | raw VAF 技術一致性訊號強，但尚無獨立 CCF/cluster/tree truth |
| 機制清晰度 | 10/20 | confound DAG 清楚；multiplicity 與 tree identifiability 仍不可觀測 |
| 反例風險 | 0/20 | 同 reads 循環、共同 caller/CN、selection bias、相似 marginal distribution 均可製造假一致 |
| 所需資源 | 10/20 | raw probe <1h；正式 recount、split-half、外部工具與 sensitivity 需 1–6h 以上 |
| **TOTAL** | **50/100** | **PROBE** |

**Falsifier observable**：若 HCC pair 的 cross-source paired VAF/CCF 差異明顯大於各來源 split-half/downsample noise、shared mutations 的獨立 cluster assignment 不穩定，或 coordinate permutation 保留相同 marginal peaks 卻破壞 mutation/edge agreement，則「clone/subclone 結構跨來源再現」假設不成立。

## §3 Assumption map

| Importance | Known | Unknown |
|---|---|---|
| High | 兩列 dataset 的 biological ID 都是 HCC1395；VCF 為 ClairS PASS SNV；shared exact alleles 可固定 | DORADO source-specific allele-CN；mutation multiplicity；within-source split-half ceiling；orthogonal tree truth |
| Low | chr1–22 scope、caller AF/DP schema | 外部工具超參數對非主要 endpoints 的小幅影響 |

High × unknown 必須先驗：不能把共用 HCC1395 CN sensitivity 寫成 source-independent CCF confirmation；不能把外部工具彼此同意寫成 biological truth。

## §4 Quick pilot 與 checkpoint

1. 固定 HCC pair exact shared alleles，比較 paired raw VAF。
   - Checkpoint：Jaccard≥0.8、Spearman與CCC≥0.8、觀察 MAE/二項抽樣預期 MAE≤1.25 才進外部 clustering probe。
   - 實際：0.9314、0.8997、0.9291、ratio=1.085，**PASS**。
2. 只在可靠 allele-CN/purity 範圍建 PyClone-VI 輸入，缺失 fail closed。
   - Checkpoint：ref+alt=DP、CN join≥95%、所有 mutation×sample rows 守恆；否則不 fit。
3. 比較 cross-source 與 within-source noise floor，並以 MOBSTER／Pairtree 作模型敏感度而非真值。
   - Checkpoint：cross-source 差異不得顯著大於 split-half ceiling，且 cluster matching CI 需穩定，才可稱 partial technical reproducibility。

## §5 Gap diagnosis

| Missing | 影響 | 工時 | 優先級 |
|---|---|---:|---:|
| DORADO source-specific allele-CN/purity | 無法正式比較兩來源 CCF | 中 | P0 |
| COLO829/HCC1937 可信 allele-CN | 7/7 PyClone-VI comprehensive claim blocked | 高 | P0 |
| HCC pair union-site BAM recount | 只用 VCF 交集會 selection-bias toward agreement | 中 | P0 |
| within-source split-half/downsample | 無法判斷 cross-source 差異是否超過技術噪音 | 中 | P0 |
| single-cell/colony/synthetic tree truth | 無法驗 biological ancestry accuracy | 高 | P1 |
| PyClone-VI/MOBSTER/Pairtree isolated environment | 尚不能執行外部工具 | 低–中 | P1 |

## §6 Evidence conflict scan

| Prior conclusion | Tier | Relation | Source |
|---|---:|---|---|
| 生物 subclone 可 infer，非確認；每區唯一真樹不可由現有 reads 保證 | L1–L2 | 限縮 claim | `InterSubMod/docs/CURRENT_FOCUS.md` |
| legacy `CCF` 實為未校正 read-AF；winner-clean 是自洽率非驗證 | L1 | 直接衝突於「用 VAF 證明樹」 | `InterSubMod/docs/methodology/20260709_ccf_readcount_tree_weighting_demo_01.md` |
| CN/purity/multiplicity 造成系統性不可識別 | L2 | 依賴關係 | `InterSubMod/docs/method_comparison/20260619_subclone_analysis_interpretation_full_framework_01.md` |
| HCC pair VAF-selected shape 77.74%，但 mutation-labeled exact forest 37.32%，皆非 accuracy | L1 | 支持 partial technical reproducibility、反對 truth claim | `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/vaf_pair_comparison.md` |

Repository 無 `MEMORY.md`；已改查 live SoT `InterSubMod/docs/CURRENT_FOCUS.md` 與 validated/in-progress method reports。

## §7 Red-team gate 與 decision path

獨立 red-team failure modes：

1. **同證據循環**：read topology 與 read-AF/VAF ranking 共用 reads；工具一致不是獨立 truth。
2. **marginal spectrum ≠ coordinate/edge identity**：座標置換可保留相同峰。
3. **CN/purity/multiplicity confounding**：共用 CN 可能人工壓低兩來源差異。
4. **technical selection bias**：兩側可評估 subset 偏向高深度與簡單候選區。
5. **forced winner**：候選多重性越高，VAF conflict 越高；需保留 set/abstain。

**Verdict：PROBE（engineering/characterization）；scientific proof of clone-tree identity = NO-GO。**

- GO：7/7 raw caller-VAF distributions、HCC exact-locus paired metrics、binomial/noise diagnostics。
- PROBE：HCC pair shared-CN PyClone-VI sensitivity；H1437/H2009/HCC1954/HCC1395 CN-aware individual fits；MOBSTER orthogonal clustering。
- BLOCKED：COLO829/HCC1937正式 CCF/PyClone；DORADO source-independent CCF；以 Pairtree top tree 當真樹；宣稱兩來源 clone evolution 已「證明相同」。
- Decision lock：**Y**。只有新增 source-specific CN、within-source ceiling 或 orthogonal truth（C1/C2/C3）才可升級 claim。

## 執行紀錄

### 輸入

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/input_manifest.json`
- manifest 所列 7 個 canonical `somatic_pass.vcf.gz`

### 命令

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AF]\n' <VCF>
python3 <inline shared-locus paired-VAF probe>
```

### 實際輸出片段

```text
N_A=80234 N_B=79120 SHARED=76848 UNION=82506 JACCARD=0.9314
SHARED_PEARSON=0.9295 SPEARMAN=0.8997 CCC=0.9291
MEDIAN_ABS_DELTA=0.0497 MAE=0.0627
OBS_MAE=0.06269 BINOMIAL_EXPECTED_MAE=0.05780 RATIO=1.085
```
