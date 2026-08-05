<!--
建立時間: 2026-07-12 00:50 +08:00
目標: 盤點 HCC1395 pair 粗拓撲整合所需癌症基因／藥物來源，驗證版本、鍵、join 膨脹與可主張邊界
處理範圍: chr1-22；全 47,377 dataset×region 母表；HCC1395 exact-coordinate complete-both 5,720 pairs
狀態: completed / scientific NO-GO
關聯檔案:
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/source_inventory.tsv
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/gene_drug_source_profile.json
-->

# 癌症基因／藥物來源與 HCC1395 pair join 稽核

> **TL;DR — 可用 GENCODE v46＋COSMIC CGC v104＋DGIdb snapshot 做「註記分層的技術一致性」；目前不能用它們證明拓撲為真或方法具臨床有效性（影響：高，信心：高）。** Exact complete pair 的全體五類 agreement 為 **3,969/5,720 = 69.39%，κ=0.497**；CGC gene-body 區 raw agreement 72.39%，DGIdb approved-antineoplastic gene-body 區 72.25%，但都只有約 +3.1 percentage points，尚未做 matched null，且後者 κ 反而未高於背景。

用 SCQA（Situation → Complication → Question → Answer）：現有本地癌症基因與藥物資料足以畫 gene-region view；但版本、many-to-many、gene length、資料庫 degree、共同 COSMIC 來源與同一 biological sample 會製造表面一致。因此本輪只把註記當 biological-plausibility sensitivity，不當 topology truth。

## 1. 結論先行

1. **Canonical region 母體已固定**：`coarse_topology_all_regions.tsv` 有 47,377 個唯一 `sample+region`；HCC1395=7,590、HCC1395_DORADO=7,268。Pair 母體 `hcc1395_pair_matches.tsv` 有 exact-coordinate 6,252 pairs，其中 complete-both 5,720。
2. **Region→gene primary join**：GENCODE v46 gene body，1-based inclusive；TSS±2 kb promoter 僅作 secondary。HGNC ID 優先，symbol 只在沒有 stable ID 時 fallback。
3. **癌症基因**：COSMIC CGC v104 共 768 genes，768 unique HGNC，無 duplicate key；可作 cancer-gene/oncogene/TSG/fusion annotation，**不可作 HCC1395 topology truth**。
4. **藥物必分三層**：
   - DGIdb `interaction`：任一 gene–drug source interaction；不表示藥物有效。
   - DGIdb `approved`：`approved=TRUE`；仍不表示該腫瘤／突變適用。
   - DGIdb `approved ∩ anti_neoplastic`：本報告最強 display subset；仍只是可追溯 interaction annotation，**不是 clinical actionability**。
5. **CLP sample-specific sensitivity 沒有形成驗證**：COSMIC CLP v104 的 1,141 個 HCC1395 chr1–22 unique alleles 中，只有 33 個 Confirmed somatic；落 exact complete pair 的 all-allele 333 regions agreement=69.97%，無命中背景=69.35%，差僅 +0.62 pp。
6. **最重要的 claim ceiling**：兩列是同一 HCC1395 biological sample 的技術處理，不是兩個獨立癌症樣本。註記分層最多支持 basecalling/pipeline robustness 的描述性合理性；不能證明 clone tree、藥物反應或方法 biological validity。

## 2. Canonical source inventory

完整 hash、mtime、絕對路徑與 decision 見 `agents/source_inventory.tsv`。

| Source | 本地版本／build | Grain／鍵 | Profile | 本輪決定 |
|---|---|---|---:|---|
| GENCODE basic GTF | v46, Ensembl 112, GRCh38, annotation 2024-03-26 | gene feature；`ENSG`/HGNC | 63,086 gene rows；61,471 symbols；40,950 HGNC | **USE**：region→gene body；promoter secondary |
| COSMIC Cancer Gene Census | v104, GRCh38 | one cancer gene；COSG→HGNC | 768 rows/768 HGNC；Tier1=592、Tier2=176 | **USE**：cancer role annotation only |
| COSMIC Genes | v104, GRCh38 | COSG mapping bridge | 58,343 rows/58,343 COSG | **USE**：避免 CGC symbol-only join |
| DGIdb interactions | overall release 未附；local mtime 2026-06-29 | gene×drug×source；HGNC | 98,239 rows；5,011 gene IDs；27,428 drug keys | **USE WITH HIGH CAVEAT**：三層 flags |
| COSMIC Resistance Mutations | v104, GRCh38 | sample×variant×transcript×drug | 5,201 rows；323 unique genomic-allele×drug | **BLOCK gene-only join**；須 exact allele |
| COSMIC CLP Genome Screens Mutant | v104, GRCh38 | sample×variant×transcript | 6,055,938 rows；HCC1395/COSS749712 另抽取 | **USE sensitivity only** |
| WAKHAN bundled `COSMIC_cancer_genes.tsv` | repo 0.4.3；檔內無版本/header | interval×symbol | 99 rows/98 symbols；僅 7 symbols 是 CGC v104 | **DO NOT USE** 作 canonical CGC |

### 2.1 版本與授權／coverage caveat

- COSMIC 各 README 只附 schema；本地資料夾沒有可直接證明再散布權限的 terms receipt。本輪只做內部衍生統計，不複製完整原始內容到報告。
- DGIdb TSV 旁沒有 overall release manifest、license 或下載 receipt；只能稱「local snapshot」，不能稱 2026 最新 DGIdb。嵌入來源版本包括 CIViC `10-Apr-24`、OncoKB `23-Jul-20`、FDA `10-Apr-24`、PharmGKB `4/5/24`、CGI `2/1/22`、ChEMBL `33`。
- GENCODE header 可確認 GRCh38/v46/Ensembl112/2024-03-26；本地 annotation 目錄未附 license file，若要對外散布需另做 terms check。
- WAKHAN code 的 MIT license **不能自動延伸**到其 bundled COSMIC-named data。

## 3. Source data quality findings

### 3.1 GENCODE stable-key 與 symbol alias

- 63,086 gene rows 的 `ENSG` 全唯一。
- 61,471 unique symbols，但有 1,615 個 duplicate symbol gene rows、488 symbols 對應多個 gene IDs。
- 21,898 gene rows 無 HGNC ID；因此 join 規則為 `HGNC > symbol fallback`，不得只以 symbol 當 primary key。
- 本輪 CGC region-gene rows 無 symbol fallback；DGIdb 僅各樣本 4 個 region-gene rows 用 symbol fallback。

### 3.2 DGIdb 三層必須分開

| 層 | Raw rows | Unique genes | Unique drugs | Unique gene–drug pairs | 正確解讀 |
|---|---:|---:|---:|---:|---|
| Interaction | 98,239 | 5,011 | 27,428 | — | 任何來源聲稱 interaction；不是療效 |
| Approved | 37,970 | 3,424 | 2,997 | 24,950 | 藥物被標 approved；不是本癌別適應症 |
| Anti-neoplastic | 14,567 | 2,089 | 716 | 9,678 | 被標抗腫瘤；可能未 approved |
| **Approved ∩ anti-neoplastic** | **9,966** | **1,873** | **263** | **6,559** | 本報告最保守的 drug display flag；仍非 actionability |

其他 QA：

- 98,239 rows 無 exact duplicate，但以 `(gene,drug,source,interaction_type)` 看仍有 5,695 duplicate-composite rows；不可直接把 rows 當獨立 evidence count。
- 8,177 rows 同時缺 canonical gene name/concept，10,572 rows 缺 canonical drug name。
- 已觀察 identity normalization anomaly：多個 ERBB2 rows 的 `drug_claim_name=Trastuzumab/Pertuzumab`，但 `drug_concept_id=rxcui:121191`、canonical `drug_name=RITUXIMAB`。因此 HTML 應保留 claim name、source、version；不得只顯 canonical drug_name，也不得據此下臨床結論。

### 3.3 COSMIC Resistance 不能 gene-only join

- 5,201 transcript/sample rows collapse 後只有 323 unique genomic-allele×drug keys，expansion ratio=16.10×。
- 檔內沒有 `SAMPLE_NAME=HCC1395` resistance rows。
- Coarse topology mother table只有 region interval，沒有完整 REF/ALT allele key；若把「region overlap 某 gene」直接接 resistance drug，會錯把別人的特定 resistance mutation套到本樣本。故本輪 **hard block**。

### 3.4 WAKHAN bundled 99-row table 不可當 CGC truth

- 99 rows、98 symbols、1 duplicate symbol；只有 7/98 symbols 出現在本地 CGC v104。
- 0 rows 的 interval 精確等於 GENCODE v46 gene body；median interval=7,842,318 bp，median interval/gene length=204.0×。
- 它可是 WAKHAN plotting subset／segment envelope，但不適合作本輪 canonical region→cancer-gene annotation。

## 4. Join 設計與守恆

1. `coarse_topology_all_regions.tsv` → 驗證：47,377 unique `sample+region`，不從 exact-tree complete-unit 子集另建母體。
2. exact complete pair → 驗證：6,252 exact-coordinate 中 5,720 complete-both；5,720 全部有一個 annotation flag record。
3. GENCODE body primary → 驗證：1-based inclusive interval overlap；每個 region 先 collapse 為 binary/list flags。
4. TSS±2 kb promoter secondary → 驗證：獨立判定 promoter-only overlap，不要求 gene body 先 overlap。
5. CGC → 驗證：COSG→HGNC bridge，768/768 CGC 有 HGNC。
6. DGIdb → 驗證：HGNC join；interaction/approved/approved-antineoplastic 分開；最終統計在 region grain，不在 exploded drug-row grain。

### 4.1 全 HCC dataset rows 的 join 膨脹

| Dataset | Region 母體 | Region→gene rows | Blowup | Exploded region×gene×drug rows | Drug blowup vs region | Max drugs/region |
|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | 7,590 | 12,555 | 1.654× | 26,414 | 3.480× | 898 |
| HCC1395_DORADO | 7,268 | 12,130 | 1.669× | 24,403 | 3.358× | 898 |

因此報告必須先 collapse：`region_has_CGC`、`region_has_DGIdb_interaction`、`region_has_approved`、`region_has_approved_antineoplastic`；不能用 raw drug-row count 比較拓撲。

### 4.2 全 region 的 body-primary coverage

| Dataset | Any gene body | CGC body | DGIdb interaction body | DGIdb approved body | Approved anti-neoplastic body |
|---|---:|---:|---:|---:|---:|
| HCC1395 | 5,964/7,590 (78.58%) | 348 (4.59%) | 1,464 (19.29%) | 1,108 (14.60%) | 607 (8.00%) |
| HCC1395_DORADO | 5,727/7,268 (78.80%) | 340 (4.68%) | 1,413 (19.44%) | 1,064 (14.64%) | 576 (7.93%) |

這些 aggregate coverage 很接近，但只能說兩套 region universes 的註記比例相近；不能取代 matched-region topology agreement。

## 5. Exact complete pair annotation-stratified agreement

母體：同座標、兩邊皆 complete 的 5,720 regions。Baseline：3,969 agree，69.3881%，Cohen's κ=0.497329。

| Primary body stratum | Present n | Present agreement | Present κ | Absent agreement | Present−absent |
|---|---:|---:|---:|---:|---:|
| Any GENCODE gene body | 4,459 | 69.61% | 0.479 | 68.60% | +1.02 pp |
| COSMIC CGC gene body | 268 | 72.39% | 0.530 | 69.24% | +3.15 pp |
| DGIdb interaction gene body | 1,104 | 70.20% | 0.473 | 69.19% | +1.01 pp |
| DGIdb approved gene body | 830 | 70.12% | 0.459 | 69.26% | +0.86 pp |
| DGIdb approved anti-neoplastic gene body | 454 | 72.25% | 0.486 | 69.14% | +3.11 pp |

判讀：

- CGC 的 κ略高於背景，但 n=268，尚無 matched null／multiple-testing correction。
- Approved anti-neoplastic 的 raw agreement 高 3.11 pp，**κ反而略低於 absent 0.498**，顯示 category mix／dominant class 可提高 raw agreement；不能宣稱 druggable regions 更穩定。
- Promoter secondary：2,427 present regions raw agreement=71.32% vs 67.96%（+3.36 pp），但 κ=0.392 vs 0.526；同樣是 prevalence-confounded 訊號，不能升格。
- Exact coordinate 代表 annotation 本身在兩列必然相同；這裡測的是「同一 annotation stratum 內 topology 是否更常一致」，不是兩套 annotation 的複現。

## 6. COSMIC CLP HCC1395 sample-specific sensitivity

輸入：`CellLinesProject_GenomeScreensMutant_v104_GRCh38.tsv.gz`。

| 階層 | 數量 |
|---|---:|
| 全檔 transcript rows | 6,055,938 |
| HCC1395 chr1–22 transcript rows | 3,373 |
| HCC1395 unique genomic alleles | 1,141 |
| Confirmed somatic unique alleles | 33 |
| Variant of unknown origin | 1,108 |
| All alleles 落 exact complete pair | 381 alleles / 333 regions |
| Confirmed somatic 落 exact complete pair | 12 alleles / 12 regions |

| CLP stratum | n regions | Agreement | κ | 無命中背景 | 差異 |
|---|---:|---:|---:|---:|---:|
| Any CLP allele | 333 | 69.97% | 0.449 | 69.35% | +0.62 pp |
| Confirmed somatic | 12 | 75.00% | 0.544 | 69.38% | +5.62 pp |
| Unknown origin | 325 | 70.15% | 0.447 | 69.34% | +0.81 pp |

**裁決**：all/unknown 的差異接近 0；Confirmed somatic 只有 12 regions，無法作穩定推論。CLP 可能來自不同 assay、cell passage、library，且 coding-screen ascertainment 與本 ONT regional reconstruction 不同；最多是 position-overlap sanity check，**不是 clone/topology truth**。CGC 與 CLP 同屬 COSMIC，不能把兩者互相富集稱獨立驗證。

## 7. 使用者指定 notable loci

| Gene | HCC1395 region | DORADO region | Pair 結果 | 可寫／不可寫 |
|---|---|---|---|---|
| BRCA2 | chr13:32315128-32345630 | chr13:32317522-32345630 | reciprocal overlap=0.9215；兩邊 complete；皆 `Topo>1 未定`，coarse agree | 可寫 coarse unresolved 狀態一致；不可寫 exact shape/tree 一致或 BRCA2 藥物反應 |
| TBC1D16 | chr17:79992680-79997005 | chr17:79992680-79995249 | reciprocal overlap=0.5941；兩邊 incomplete | 不可評五類 topology agreement；TBC1D16 也不是 CGC v104、DGIdb 本 snapshot 無 interaction |
| ERBB2 | chr17:39673348-39720521 | 無 DORADO body-overlap region | HCC1395 該區 incomplete；unmatched | 不可作 pair concordance；DGIdb 高 degree 更需避免 interaction-count bias |
| MYC | 無 | 無 | 無 region evidence | 只能列 database context，不可列本 pair finding |

這也說明為何只用 exact complete pair 的 gene-level表時，上述四基因皆是 n=0：BRCA2/TBC1D16 的 boundary 不完全相同，ERBB2 不成對，MYC 不在 region universe。Gene HTML 應把「不存在／incomplete／overlap-matched」清楚分開。

## 8. 舊 annotation 不可重用

| Dataset | 舊 annotation keys | 最新 keys | Exact overlap | 最新 coverage |
|---|---:|---:|---:|---:|
| HCC1395 | 3,885 | 7,590 | 2,341 | 30.84% |
| HCC1395_DORADO | 2,379 | 7,268 | 954 | 13.13% |

另外，舊 `region_gene_annotation.py` 有三個方法風險：

1. `ge < region_start` 時在 promoter test 前 `continue`，會漏掉 promoter-only overlap。
2. `genes[:30]`／`protein_coding[:30]` 與每基因 drugs `[:8]` 會截斷大區域與高-degree gene。
3. 直接以 uppercase symbol 接 DGIdb/COSMIC，沒有 HGNC bridge。

本輪重算已修正這三點；舊 JSON 僅可當歷史 UI asset。

## 9. 必要 matched-background null

目前 +0.9–3.1 pp 都是 unadjusted descriptive differences。若要測 enrichment，建議：

1. **分析單位**：5,720 exact complete pairs；outcome=`category_agree`，不能使用 exploded gene/drug rows。
2. **Primary matching**：chromosome、log(region length) decile、primary HP-unit count、region sSNV count、read support/depth、CN state、callable bp／local mutation density。
3. **Gene opportunity**：匹配 region 的 GENCODE gene-covered bp、gene count；gene-level分析再匹配 gene length與 mutation opportunity。
4. **Drug database degree**：按 gene 的 DGIdb unique-drug degree分層／匹配，或只用 binary flags；ERBB2 等高-degree genes不能讓 192 drugs 當 192 個獨立證據。
5. **空間相依**：用 chromosome-preserving、≥1 Mb block permutation／block bootstrap，避免相鄰區域與重疊 window 被當獨立。
6. **負控制**：等長 shifted regions、非 CGC degree-matched genes、隨機 approved non-antineoplastic gene sets。
7. **多重比較**：五個 primary annotation strata預先註冊；若展開逐 gene，做 FDR；報 effect size與 CI，不只 p-value。
8. **來源獨立性**：COSMIC CGC 與 COSMIC CLP 的 enrichment 不算獨立 replication；需 non-COSMIC orthogonal truth 或 single-cell/multi-region tree。

## 10. 可 claim 與不可 claim

### 可以寫

- 本地來源與 join 工程可追溯；HCC1395 pair exact complete coarse-category agreement=69.39%，κ=0.497。
- CGC-body／approved-antineoplastic-body strata 的 raw agreement 約高 3 pp，是待 matched-null 檢驗的 hypothesis-generating signal。
- BRCA2 兩套 boundary 高度重疊，兩邊皆 complete 且 coarse class 都是 `Topo>1 未定`。
- CLP all-allele命中區的 agreement 與背景近似，未見明顯 validation uplift。

### 不可以寫

- 「兩個獨立樣本重現」：兩列是同一 biological HCC1395。
- 「癌症基因／藥物資料證明樹正確」：annotation 沒有 topology label。
- 「方法已證明合理有效」：只有 moderate technical concordance，且缺 clean-v3＋single-cell/multi-region truth。
- 「DGIdb interaction 表示有效治療」或「approved 表示適應於此樣本」。
- 「COSMIC resistance gene 命中表示抗藥」：必須 exact allele，gene-only join 已阻斷。
- 「CLP 是獨立驗證」：不同 assay/passage與共同 COSMIC ascertainment 未排除，且 97.1% HCC alleles 是 unknown origin（1,108/1,141）。

## 11. 執行紀錄

**主要輸入**：

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv`
- `/big7_disk/liaoyoyo2001/gene_annotation/gencode.v46.basic.annotation.gtf.gz`
- `/big7_disk/liaoyoyo2001/gene_annotation/Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz`
- `/big7_disk/liaoyoyo2001/gene_annotation/cosmic_v104/Cosmic_Genes_v104_GRCh38.tsv.gz`
- `/big7_disk/liaoyoyo2001/gene_annotation/dgidb_interactions.tsv`
- `/big7_disk/liaoyoyo2001/gene_annotation/cosmic_v104/CellLinesProject_GenomeScreensMutant_v104_GRCh38.tsv.gz`

**執行命令**：

```bash
python3 -m py_compile research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/profile_gene_drug_sources.py
python3 research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/profile_gene_drug_sources.py
```

**實際輸出片段**：

```text
coarse_region_rows=47377
exact_coordinate_rows=6252
exact_complete_pair_rows=5720
baseline_agreement=3969/5720=69.3881%; kappa=0.497329
HCC1395_CLP_unique_alleles=1141; confirmed=33; unknown_origin=1108
CLP_overlap_exact_complete=381 alleles / 333 regions
```

**主要輸出**：

- `agents/source_inventory.tsv`
- `agents/gene_drug_source_profile.json`
- `agents/hcc1395_pair_region_gene_join.tsv.gz`
- `agents/hcc1395_exact_complete_pair_gene_drug_flags.tsv`
- `agents/hcc1395_exact_complete_annotation_agreement.tsv`
- `agents/hcc1395_cosmic_clp_v104_chr1_22_alleles.tsv.gz`
- `agents/hcc1395_exact_complete_clp_sensitivity.tsv`
- `agents/hcc1395_exact_complete_gene_summary.tsv`
- `agents/hcc1395_notable_gene_locus_sensitivity.tsv`

## 12. Data-quality verdict

| Finding | Severity | Confidence | Remediation |
|---|---|---|---|
| 舊 annotation coverage 僅 30.84%/13.13% | High | High | 已用 canonical mother tables 重算 |
| DGIdb overall version/license receipt 缺、identity anomaly | High | High | 只作 exploratory flags；對外前重抓有 receipt snapshot並逐 drug核對 |
| Many-to-many drug blowup最高 898/region | High | High | region binary collapse＋degree-matched null |
| CGC/CLP 同源循環與 coding-screen ascertainment | High | High | 不稱獨立驗證；另找 orthogonal truth |
| Symbol 非唯一／alias drift | Medium | High | HGNC/COSG bridge優先；symbol fallback另標 |
| Confirmed CLP overlap只有 12 regions | High | High | 不做 effect claim；需較大 orthogonal truth set |

**Final verdict：資料工程 PASS；生物／臨床驗證 NO-GO。**

## 13. Parent integration update：matched-null 已完成

> 本段更新 supersedes §5／§9 中「尚未做 matched null」的暫時描述；inventory 原始數量不變。

主分析後續以 `agents/hcc1395_exact_complete_pair_gene_drug_flags.tsv` 的 5,720 exact complete pairs，固定 seed=20260712，做 5,000 次 chromosome＋global region-length-decile conditional hypergeometric permutation。Canonical 輸出為：

- `data/hcc1395_annotation_reproducibility.tsv`
- `data/hcc1395_annotation_reproducibility.json`
- `scripts/analyze_annotation_reproducibility.py`

| Stratum | Descriptive delta | Conditional-null p |
|---|---:|---:|
| GENCODE gene body | +1.02 pp | 0.4271 |
| COSMIC CGC body | +3.15 pp | 0.5855 |
| DGIdb interaction body | +1.01 pp | 0.7600 |
| DGIdb approved body | +0.86 pp | 0.8486 |
| DGIdb approved∩anti-neoplastic body | +3.11 pp | 0.4295 |
| COSMIC CLP all-status allele-containing region | +0.62 pp | 0.3195 |
| COSMIC CLP confirmed-somatic region | +5.62 pp（n=12） | 0.7313 |

全部 p>0.05；TSV/JSON 第二次執行 byte-identical。這使 claim ceiling 更明確：目前沒有 evidence 顯示 cancer-gene、drug-annotated 或 CLP-containing regions 的 topology agreement 超越可比 chromosome/length background。CN、read depth、HP count、k、gene opportunity 與 block dependence 仍未完全調整，因此這不是「證明無關」，但足以否定把本輪約 +3 pp 描述差異當成 validation evidence。
