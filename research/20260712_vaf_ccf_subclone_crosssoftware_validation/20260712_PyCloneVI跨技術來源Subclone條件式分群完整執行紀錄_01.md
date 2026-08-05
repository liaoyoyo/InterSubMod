<!--
建立時間: 2026-07-12
目標: 以 PyClone-VI 對 HCC1395 跨技術來源及具備可審核 SAVANA allele-specific CN 的樣本做條件式 subclone clustering 驗證
處理範圍: GRCh38 chr1-22、biallelic sSNV、exact CHROM:POS:REF:ALT；14 個 main/sensitivity fits
關聯檔案: InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/{config,data,runs,report}/
-->

# PyClone-VI 跨技術來源 subclone 條件式分群：完整執行紀錄

> **Task type B — comprehensive validation；服務 G4/G5。** 14/14 預註冊 fit 均 PASS；可支持「HCC1395 主 clonal 骨架與低 CCF 成分訊號跨來源再現」，但獨立 fit 的 minor-clone 位點歸屬只有中度一致，**不可**宣稱唯一、準確且相同的真實演化樹。

## 1. 輸入契約

- Source manifest：`InterSubMod/research/20260710_layered_reconstruction_v2/data/layered_input_manifest_v3_raw_all_lps_pass_v1.json`。
- HCC primary universe：兩側 latest LongPhase-S recalibrated PASS tree VCF exact intersection，再交 SEQC2 high-confidence sSNV exact intersection，共 **30,265** 位點；3 位點無可用 SAVANA measured CN，PyClone input 為 **30,262**。
- HCC PASS-union sensitivity：兩側 PASS union ∩ SEQC2 = 31,429；兩側 recalibrated-all VCF 都有 DP/AF counts = **30,393**；再排除 3 個無 measured CN 位點，input = **30,390**。此 counts 來自 VCF，**不是完整 BAM recount**。
- H1437／H2009／HCC1954 tree PASS ∩ orthogonal-tools benchmark exact universe分別為 **74,328／150,726／20,992**；CN fail-closed 後 main input 為 **73,124／150,714／20,981**。
- HCC1395 SAVANA：`/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv`，purity 0.96。DORADO 共用同細胞株 CN，明確標為 **shared-CN sensitivity**，不是獨立 CN measurement。
- H1437／H2009／HCC1954 purity 分別為 0.95／0.95／0.66，使用各自 SAVANA `cna_normalhet/*_segmented_absolute_copy_number.tsv`。
- COLO829、HCC1937 因缺 reviewed measured allele-specific CN 而 **BLOCKED**；禁止用 tumour CN=2 代填。
- Counts：`alt=floor(DP×AF+0.5)`、`ref=DP-alt`；14 inputs 的守恆、complete sample×mutation matrix、mutation/sample uniqueness 全部 PASS。
- Main CN 離散化：`total_i=floor(total+0.5)`；`minor_i=min(floor(minor+0.5), floor(total_i/2))`；`major_i=total_i-minor_i`。Near-integer sensitivity 另要求 total/minor 距最近整數都 ≤0.25。

## 2. 實際命令

### 2.1 建立 inputs

```bash
python3 InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/scripts/build_pyclone_inputs.py \
  --config InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/config/pyclone_validation_config.json \
  --output-dir InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/data/pyclone_inputs
```

### 2.2 Fit 命令契約

每個 job 的完整 argv、stdout、stderr、HDF5、results 與 SHA-256 都保存在 `InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/<bundle>/status.json`。核心命令為：

```bash
InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/external_tools/pyclone_vi_env/bin/pyclone-vi fit \
  -i <bundle>.pyclone_input.tsv.gz -o <bundle>/model.h5 \
  -c 40 -d beta-binomial -r <20-or-10> -t 4 --seed 20260712
```

- Primary joint main、primary joint near、兩側 separate main/near、三個 individual main、PASS-union main：`r=20`。
- 非 primary sensitivity（H1437/H2009/HCC1954 near、PASS-union near）：明示使用 `r=10`；沒有默默變更。
- PyClone-VI 版本：0.2.0。14/14 `seed=20260712`、`c=40`、beta-binomial。

### 2.3 最終分析與 artifact

```bash
InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/external_tools/pyclone_vi_env/bin/python \
  InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/scripts/analyze_pyclone_results.py \
  --run-dir InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi \
  --output-dir InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/analysis

python3 InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/build_report_artifact.py \
  --raw-dir InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/raw_vaf_validation_v1 \
  --pyclone-analysis-dir InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/runs/pyclone_vi/analysis \
  --topology-context InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/integrated_topology_context_v1.json \
  --output InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/artifact.json
```

Artifact builder：**PASS**，24 datasets、46 blocks、7 charts、15 tables、7 sources；snapshot 為 partial，因完整主張仍受 CN 缺口、conditional clustering 與 historical topology evidence ceiling 限制。

### 2.4 Portable HTML package 與瀏覽器 QA

```bash
node InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/package_portable_report.mjs \
  --input InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/artifact.json \
  --output InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/20260712_HCC1395跨來源VAF_CCF與subclone驗證_01.html \
  --plugin-root /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599 \
  --receipt InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/portable_package_receipt.json

python3 InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/qa_portable_report.py \
  --html InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/20260712_HCC1395跨來源VAF_CCF與subclone驗證_01.html \
  --output-dir InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/qa
```

Portable package 使用 Data Analytics 官方 `buildPortableArtifact`、`extractPortableChartSvgs` 與 `verifyPortableArtifact`。因 reader 的 `100vw` topbar 與 Recharts legend 在 non-overlay scrollbar 環境造成 8–97 px 假性水平溢位，package script 只加入 `.analytics-top-bar` 與 `.chart-legend` 的寬度相容性 override；**沒有修改資料或 chart spec**。官方驗證為 49 rendered blocks、7 charts、4 metric cards、15 tables，1440/390 viewports 與 source dialog PASS；獨立 Playwright QA 為 18/18 PASS、console/page errors=0。

## 3. HCC1395 跨來源結果

### 3.1 Joint fit 是描述，不是獨立一致率

Primary joint main 以 30,262 位點得到 4 群：

| mutation fraction | HCC1395 CP | DORADO CP | 分類 |
|---:|---:|---:|---|
| 95.29% | 1.0000 | 1.0000 | clonal |
| 4.03% | 0.6268 | 0.6461 | subclonal |
| 0.54% | 0.3054 | 0.3049 | subclonal |
| 0.14% | 0.8437 | 0.8844 | subclonal |

mutation-weighted mean absolute cluster-CP delta = **0.000838**；CP Spearman = 1.0。但 joint fit 共享 cluster label，不能把 label agreement 當獨立 reproducibility。

### 3.2 真正關鍵：兩側各自 fit

| 指標 | Main 30,262 | Near-integer 22,051 | 解讀 |
|---|---:|---:|---|
| ARI | 0.539 | 0.573 | 全體僅中度 label agreement |
| NMI | 0.339 | 0.390 | 同上 |
| Hungarian overall | 98.29% | 98.96% | 被約 98% clonal majority 支配，不可單獨引用 |
| subclonal mutation Jaccard | **0.381** | **0.411** | minor-clone 位點歸屬僅中度再現 |
| clonal/subclonal κ | 0.544 | 0.577 | 中度一致 |
| both-subclonal cluster ARI | 0.325 | 0.552 | minor 群內結構仍會變 |

Main separate fit 的 subclonal contingency：兩側都 subclonal 277、僅 HCC 235、僅 DORADO 215、either 727。這是判斷「是否同一 subclone 結構」時應使用的分母，而不是 98.29% overall Hungarian agreement。

### 3.3 Model-mode sensitivity

同一來源在 joint vs separate mode 的 subclonal-set Jaccard：HCC1395 **0.3375**、DORADO **0.3347**。這證明精確 subclone 歸屬對「joint 共享群結構」或「各自獨立 fit」的模型選擇相當敏感；不是生物真值。

### 3.4 Universe 與 CN sensitivity

- Primary joint main vs PASS-union main：共同位點 ARI **0.9932**、NMI 0.9596、CP Spearman 0.9942，對 mutation-universe 小幅擴張很穩。
- Primary joint main vs near-integer：ARI **0.8656**、NMI 0.7720、subclonal Jaccard **0.7810**、κ 0.8727；大致穩定，但不能說完全不受 CN discretization 影響。
- Separate fit 各自 main→near 的 subclonal Jaccard：HCC1395 0.8896、DORADO 0.8410。CN 篩選不是兩側 minor-clone 只中度重現的唯一原因。

## 4. 其他樣本 conditional profiles

| Dataset | Main mutations | 群數 | Clonal fraction | Subclonal fraction | assignment prob mean | <0.8 | Near sensitivity |
|---|---:|---:|---:|---:|---:|---:|---|
| H1437 | 73,124 | 5 | 95.04% | 4.96% | 0.935 | 8.11% | ARI 0.974；subclonal J=0.956；κ=0.976 |
| H2009 | 150,714 | 4 | 96.71% | 3.29% | 0.726 | 98.30% | ARI 0.951；subclonal J=0.911；κ=0.952 |
| HCC1954 | 20,981 | 3 | 83.99% | 16.01% | 0.596 | 97.07% | **ARI 0.313；subclonal J=0.253；κ=0.365** |

- H1437/H2009 的 common near-integer subset labels 整體穩定；H2009 main 的 posterior assignment confidence 卻很低，故仍不能把 4 群當確定 clone 數。
- HCC1954 near-integer 僅保留 4,173/20,981（19.9%），且 main assignment confidence 已低；其 clone profile 對 CN gate 不穩健。
- 不同 cell line 不是 biological replicate；只能比較 profile，不能以群數相近證明相同演化。

## 5. QA 與 claim ceiling

- 14/14 status PASS；14/14 seed=20260712；14/14 `num_clusters=40`。
- 10 jobs r20、4 個明示 non-primary sensitivity jobs r10。
- Input mutation count、result mutation count、result rows、expected complete rows：**0 mismatch**。
- 所有 fit/write-results stderr：**0 個非空檔**。
- `analysis_summary.json`：`all_ready=true`、`pass_bundle_count=14`、`pending_or_failed_bundles=[]`。
- Required analysis tables 與 blocked ledger 全部非空；report artifact builder fail-closed PASS。
- Portable package 與官方 browser verifier PASS：49 rendered blocks、7 charts、4 metric cards、15 tables，viewports 1440/390，source dialog PASS。
- 獨立 Playwright QA 18/18 PASS；desktop/mobile full-page screenshots 已保存；console/page errors=0。
- PyClone-VI 只產 CCF/mutation clustering，**不產演化樹**。這是 conditional external clustering；不能當 InterSubMod 的獨立 gold truth。

主要 SHA-256：

| Artifact | SHA-256 |
|---|---|
| input provenance | `2398a2bcea8487b1d0f1b273fa89b12a3db9d20a68ed9399ed4b9ba8f0f0d343` |
| consolidated status/command manifest | `d65391894b6b6e1b134754eb86cb151656082310baa2a27f12951543aacc9fbe` |
| final analysis summary | `ec33513b200df4c1273bfa97929271f6f1dcfbc86d1bcddf46680334cabde3ad` |
| report artifact | `4e6bd9c95f01d9b311d641441cc7b0148492a246638909ecc1b969ee58d38138` |
| portable HTML | `6ae05e9d9368e9daca7468cc09bfb9185ee2f9be735f07323c79bc11c9ddc67b` |
| portable package receipt | `581c24fe9e3f71627fca66efa0151193dcf464352ac51e2472bf7ac7ec6e5db5` |
| browser QA receipt | `9a34427420623b35698bd493b4b917565b68e329f9075e9d7a97383e2c2bbdd1` |

## 6. 結論

PyClone-VI 佐證兩個 HCC1395 技術來源都有一個占絕大多數的 clonal 主幹與低 CCF 成分；在相同 universe 上的獨立 fit，clonal/subclonal 二分 κ 約 0.54–0.58，但精確 subclonal 位點 Jaccard 只有 0.38–0.41。這支持「**主要結構訊號可跨技術來源再現**」，不支持「**每區或每位點能重建唯一且相同的真實 clone／演化樹**」。方法可作結構訊號偵測與候選生成；精確 subclone/樹仍需獨立 CN、單細胞、多時間點或真正多區域樣本驗證。
