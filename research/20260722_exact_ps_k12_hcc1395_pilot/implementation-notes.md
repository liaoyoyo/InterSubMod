<!--
建立時間: 2026-07-22
目標: 逐步記錄 HCC1395 exact PS × HP k<=12 pilot 的設計、偏離、折衷與未決事項
處理範圍: PARTIAL / exploratory pilot / HCC1395 chr1-22 only
關聯檔案: InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/pre-decision-audit.md
-->

# HCC1395 exact PS × HP、k≤12 implementation notes

> **PARTIAL / EXPLORATORY PILOT**：HCC1395 chr1–22 only；等待使用者確認後才可擴樣本。

## 設計決定

- [2026-07-22] primary key 固定為 exact known `PS × HP1/HP2`，不同 PS 不合併。
- [2026-07-22] k＞12 僅在同一 evidence unit 內分割；上限 12 是 solver 工程界線。
- [2026-07-22] Python 為 producer，C++ 為獨立 parity kernel；兩者讀同一份由 big7 BAM 產生的 normalized sparse evidence。
- [2026-07-22] partial call `X` 不複製成多個 read；segmentation link 只用 fixed R/A，後續 tree likelihood 再直接 marginalize X。
- [2026-07-22] `constraints.tsv.gz` 是 **segmentation-only** aggregate，不是後續 likelihood 的完整輸入；完整 R/A/O/D/S/L/X 與 molecule identity 仍以 extraction 的 `molecule_sparse_calls.tsv.gz` 為 authority。未建立 block→partial-molecule likelihood adapter 前，不可宣稱新版已重跑 tree／VAF ranking。

## 偏離

- 目前不直接修改已 dirty 的 InterSubMod core CMake／C++，也不碰正在進行中的 LongLineage all-dataset compatibility write-set。先在隔離 research topic 建立可編譯 C++ pilot，避免污染既有 authority。
- 本輪不實作 signed cross-PS SAME／FLIP bridge；現有 scalar sidecar PS 沒有足夠資訊支持該判定。

## 折衷

- bounded contiguous DP 可最大化局部保留 read constraints，但跨切點 constraint 仍會被切斷；因此報告 retained、cut、unavoidable，而不宣稱 lossless global topology。
- HCC1395 全常染色體是單樣本完整 pilot，但 cohort 結論仍是 partial。

## 未決

- signed germline-anchor bridge 是否足以安全合併不同 PS，需另立方法與 golden fixtures。
- k＞12 是否改採 separator-aware exact decomposition，而非 bounded local blocks，需獨立效率與 certificate 研究。
- downstream tree／likelihood adapter 須在 segmentation parity 通過後才啟動。
- 新版 whole-flow runner 本輪的終點固定為 extraction → exact-PS segmentation → C++ parity；不把尚未實作的 topology／VAF adapter 偽裝成已完成 stage。

## 執行紀錄

### 2026-07-22｜輸入綁定

- 輸入 BAM：`/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam`
- 輸入 BAI：`/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam.bai`
- Pilot manifest：`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/input_contract/big7_hcc1395_pilot_manifest.json`
- 輸出 receipt：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/input_binding/big7_hcc1395_input_binding_receipt.json`
- 驗證：37/37 checks PASS；BAM fixed first/middle/last chunks、BAI full SHA-256、tree VCF/index 與 HP/PS sidecar/index fresh full SHA-256、`samtools quickcheck` 全通過。

### 2026-07-22｜第一次 direct-big7 chr22 attempt（保留失敗證據）

- 命令使用原 20260716 M2 extractor。
- 輸出：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/smoke_chr22/`
- 結果：exit 1；`samtools view -M -L` 因重疊 target intervals 重複輸出同一 canonical alignment identity，原 extractor fail closed。
- 實際錯誤：`duplicate alignment emitted by multi-region iterator`。
- 判定：不是 PS 或 BAM identity 錯誤；不可刪除此失敗證據，也不可把半成品續作成功輸出。

### 2026-07-22｜identity-safe collapsing extractor 修正路徑

- 稽核 20260716 原版與 20260718 collapsing 版 source diff。
- collapsing 版只在 canonical alignment identity 重複時，比對 target calls、BQ、MAPQ、HP、PS；全部相同才保留一筆，任一衝突仍 fail closed。
- Targeted regression：`10 passed`。

### 2026-07-22｜direct-big7 chr22 E2E PASS

- 輸出：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/smoke_chr22_collapsing_v1/`
- Extraction：exit 0；9:12.23；53,764 unique molecules；9,660 duplicate identities一致折疊；conflict=0；17/17 checks PASS。
- Python exact PS×HP partition：exit 0；2.47 秒；530 units、925 unit-site memberships、535 blocks、940 patterns。
- Constraint molecule-component incidence mass：5,576 = retained 5,522 + cut 5 + unavoidable 49。
- C++ kernel：exit 0；0.02 秒；530 units、535 blocks、940 patterns。
- Python/C++ parity：`mismatch_count=0`；blocks、membership、dispositions、counts、weights逐列一致。
- 與 20260718 cached chr22：五個 extractor TSV 的 decompressed SHA-256 全相同；physical gzip bytes 因 gzip header 不同，不能拿來判定科學差異。新版五個 partition gzip artifacts 則 byte-for-byte 相同，semantic SHA-256 皆為 `08125edb97610cd882e66d4b66a942c5b867a3ae73fadac794b083f4aca67146`。

### 2026-07-22｜測試

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests
```

第一次完整 topic regression：`45 passed in 24.96s`，exit 0；包含 Python producer、C++ arbitrary-precision DP、input binding、Python/C++ comparator 與 tamper negative cases。

### 2026-07-22｜whole-flow runner preflight 與中止紀錄

- 第一次 full-direct runner 目標：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v1/`。
- 程序在 pre-run input verifier 尚未完成時，由主 agent 因獨立 code review 回報 4 個 provenance／denominator P1 而主動中止（exit 130）；**未開始任何 chromosome extraction**。該目錄保留作 aborted evidence，不續跑、不覆寫。
- 修正項目：
  1. input binding receipt 強制驗證 schema、sample、manifest SHA、37 checks、semantic integrity、實際 samtools；並新增 pre/post-run snapshot equality。
  2. `S` 與 active unique sites 改以 catalog `(chrom,pos1)` 驗證；拒絕重複 locus、重複 site index 與 membership/catalog mismatch。
  3. aggregate totals 只納入 `all_pass=true` 且 metrics 完整的 chromosome，另列 requested／with-metrics／included／excluded denominator。
  4. collapsing extractor 新增 explicit `--samtools`；receipt 同時記錄 executable 與實際 command，runner 再驗證兩者等於 input verifier 使用的同一絕對路徑。
  5. historical comparison 指定後，6 個 artifact 任一 decompressed semantic mismatch 即 fail closed；gzip header-only 差異不算科學差異。
- 修正後驗證：topic `52 passed in 27.38s`；20260718 segmentation suite `88 passed in 12.26s`；`git diff --check` 無輸出，皆 exit 0。

### 2026-07-22｜HCC1395 chr1–22 direct-big7 v2 PASS

- 輸入 BAM：`/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam`
- 輸入 candidate VCF：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395.longphase_s.recalibrated.pass.vcf.gz`
- 輸入 HP/PS sidecar：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395.read_tags.tsv.gz`
- 輸出 root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2/`

執行命令：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  research/20260722_exact_ps_k12_hcc1395_pilot/scripts/run_hcc1395_exact_ps_k12_pilot.py \
  --manifest research/20260722_exact_ps_k12_hcc1395_pilot/input_contract/big7_hcc1395_pilot_manifest.json \
  --output-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2 \
  --workers 1 \
  --historical-extraction-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1/chromosomes
```

實際結果：

```text
22/22 chromosomes PASS
cross_ps_violations=0
cross_hp_violations=0
python_cpp_mismatch_count=0
historical semantic mismatches=0/132
```

- Run receipt：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2/run_receipt.json`
- Receipt SHA-256：`7a562fbcb8dfa22abfd81b67f8eb092b6bdef519f47d71512309f2a978ff30c0`；`_SUCCESS` 綁定相同 path／size／SHA。
- Pre/post input-binding snapshot 相同：`5976a574849525e434004c0717402da562b890923e7bd068ec6a6b2ec3a58afb`。
- 22/22 stage sidecars、22/22 extraction sidecars與 902 個 receipt-declared file identities 現場重 hash 全吻合。
- 獨立 final audit：5,342 checks、0 failures、verdict=GO。

### 2026-07-22｜direct gzip 重算與分母

逐染色體直接讀取 `site_catalog`、`units`、`membership`、`blocks`、`constraints`、`dispositions` gzip TSV，最後才與 run receipt 比對：

| 指標 | 數量 | 單位／分母 |
|---|---:|---|
| Manifest-bound candidate loci S | 79,687 | unique `(chrom,pos1)`；不是本輪 BAM de novo calling truth |
| Primary unique loci | 36,384 | 至少進一個 exact PS×HP1/2 unit |
| Unit memberships | 62,651 | locus–unit incidences，不是 unique loci |
| Exact PS×HP read-linked units | 39,544 | evidence units |
| Output blocks | 39,807 | k≤12 blocks |
| Constraint rows | 58,014 | aggregated patterns |
| Molecule weight | 285,897 | molecule–unit incidences，不是全域 unique reads |

k 分布：k=1 28,109；k=2–8 11,253；k=9–12 88；k＞12 94。

全部 units 的 molecule-weight disposition：

- retained 281,967／285,897＝98.6254%
- cut 1,254／285,897＝0.4386%
- unavoidable 2,676／285,897＝0.9360%

只限原始 k＞12 units：94 units、3,277 memberships、2,070 unique loci、357 blocks；retained 9,295／13,225＝70.2836%、cut 1,254＝9.4820%、unavoidable 2,676＝20.2344%。全部 cut／unavoidable 都來自這 94 units。

### 2026-07-22｜exact PS census 與不跨接依據

- `(chrom, exact PS)` containers：3,971；HP1+HP2 3,147（79.2496%）、HP1-only 441、HP2-only 383。
- `PS×HP` routes：7,118。
- 1,668 個 unique `locus×HP` keys 落在兩個 PS；合併 HP 後是 1,244 個 unique genomic loci，佔 primary loci 3.4191%。
- 這 1,668 個案例全部恰好兩個 PS；因此若只依 HP family 跨 PS 聚合，存在重複與未校正 SAME／FLIP orientation 的實際風險。

### 2026-07-22｜chr6／chr16 upstream scope outlier

- chr6＜60 Mb：25,672 candidate loci，primary 653（2.54%），落 SEQC2 HC BED 0；known-PS HP1/2 molecule rows 0.76%。
- chr16≥45 Mb：18,030 candidate loci，primary 585（3.24%），落 SEQC2 HC BED 0；known-PS HP1/2 molecule rows 1.33%。
- 兩臂合計佔 S 43,702／79,687＝54.84%，卻只佔 primary 1,238／36,384＝3.40%；排除兩個診斷臂後 primary coverage 35,146／35,985＝97.67%。
- Catalog 與 manifest 綁定 LongPhase-S PASS VCF 逐列一致、無 duplicate loci、targets 為 1-bp；不是 extractor 製造。
- SEQC2 KB 把 chr6p／chr16q 記為 complete-loss／benchmark-excluded。可判定為 candidate-universe／phase-eligibility／CN-exclusion contract 的 P0 QA 訊號；不可直接宣稱 43,702 loci 全是 germline、caller FP 或應刪除。
- Pilot manifest 目前把 CN BED 未列位置解讀為 neutral 且未綁 exclusion BED；segmentation 未使用 CN，故本輪技術數字不變，但 downstream CN／clone 解讀前必須修成 excluded/unknown。

### 2026-07-22｜HTML technical report 與 QA

- Artifact：`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/report/artifact.json`
- Evidence snapshot：`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/report/evidence_snapshot.json`
- HTML：`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/report/report.html`
- Delivery receipt：`InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/report/report_delivery_receipt.json`
- Browser QA：2 charts、7 tables、6 metric cards；1440px desktop、390px mobile、source dialog keyboard interaction 全 PASS。
- 共用 portable reader 在 Linux non-overlay scrollbar 下因 top bar `width:100vw` 固定產生 8px overflow。重用 20260718 已測試的 scrollbar-safe delivery wrapper（Node 3/3 tests PASS）；未改 artifact payload 或 chart runtime，只在 packaged reader copy 加入具名 CSS correction，然後由 canonical deliver/verify lifecycle 通過完整 QA。失敗 screenshots 保留供 renderer bug 稽核。

### 尚未完成／不得宣稱

- block→完整 partial-molecule likelihood adapter 尚未實作。
- tree enumeration、VAF ranking、T／Topo、clone／subclone before/after 尚未重跑。
- 七樣本擴展尚未啟動；等待使用者確認 HCC1395 policy 與 outlier 處置後才可進行。
- Topic 程式與報告目前尚未 commit；正式 handoff 前需由使用者授權 commit。

### 2026-07-23｜正式 topology adapter 整合啟動

- [決策] 使用者確認不同 exact PS 默認必須分開；本輪由 exploratory segmentation 進入 HCC1395-only production-path adapter。
- [決策] 新輸入單位是 `exact PS × HP1/HP2 × read-linked component × bounded block`；不同 PS 不因為 HP label 相同而合併。
- [折衷] 不直接改 dirty `RegionProcessor` 為 topology engine；它是單位點甲基化分析，並沒有 multi-sSNV R/A/X schema。exact-PS C++ partitioner 作為獨立正式 build/validation gate。
- [未決] 本輪完成 T/Topo 之後，VAF/read-pattern likelihood ranking 仍需以新 unit key 重新 join，不沿用舊 region 的 ranking 數字。

### 2026-07-23｜HCC1395 exact-PS topology 整合與全量驗證 PASS

#### [決策] 新版權威分析單位

- Tree input 固定為 `exact PS × HP × read-linked component × bounded block (k=2..12)`。
- 相同 HP label 只在同一 exact PS 內聚合；不同 PS 不合併。
- k＞12 先由 Python/C++ parity 已通過的 retained-read-weight DP 切成 k≤12 blocks，再由完整 molecule R/A/O/D/S/L/X rows 建立 partial-read patterns。
- `region_id = chrom|PS|HP|unit:block`，不再以 coordinate-only key join，避免同座標不同 PS 覆寫。
- C++ exact partition kernel 已加入獨立 CMake target `exact_ps_partition`；不把單 anchor 甲基化 `inter_sub_mod` 偽裝成 multi-sSNV topology engine。

#### Step → Verify 與實際輸入／輸出

1. Partition → mlhp adapter。  
   → 驗證：adapter receipt `all_pass=true`；groups=11,542；cross-PS=0；cross-HP=0；Python/C++ mismatch=0。
2. Layered tree 全候選重算。  
   → 驗證：`ANALYSIS_TREE_CAP=0`、`DISPLAY_TREE_CAP=0`；11,122 non-capped units V1–V7 全通過，420 capped 為 not-applicable。
3. Region view 與 summary。  
   → 驗證：adapter/layered/region-view 三階段 stable ID set 都是 11,542，逐區 exact PS/HP 一致且 ID 唯一。
4. Build、regression、receipt。  
   → 驗證：CMake target build exit 0；topic 59/59 tests PASS；repo CTest 251/251 PASS；final receipt 9/9 SHA identities PASS。

輸入：

```text
BAM=/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam
VCF=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395.longphase_s.recalibrated.pass.vcf.gz
PARTITION_ROOT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2
```

核心執行命令：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/exact_ps_partition_to_mlhp.py \
  --partition-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2 \
  --output /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_mlhp_HCC1395_chr1_22.json \
  --sample HCC1395 \
  --chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
  --min-read 3

SM_ML=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_mlhp_HCC1395_chr1_22.json \
SM_OUT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/layered_reconstruction_HCC1395_chr1_22_display_all.json \
SM_ANALYSIS_TREE_CAP=0 SM_DISPLAY_TREE_CAP=0 SM_VERIFY_EVERY=1 \
SM_CN_INT_GAIN= SM_CN_INT_LOSS= PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py

SM_LAYERED=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/layered_reconstruction_HCC1395_chr1_22_display_all.json \
SM_OUT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/layered_region_view_HCC1395_chr1_22_display_all.json \
SM_ML_GLOB=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_mlhp_HCC1395_chr1_22.json \
SM_SOMATIC_VCF=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/HCC1395/HCC1395.longphase_s.recalibrated.pass.vcf.gz \
SM_BACKBONE_SOURCE=longphase_s_recalibrated_FILTER_PASS SM_SAMPLE=HCC1395 \
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_region_view.py
```

權威輸出：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_mlhp_HCC1395_chr1_22.json
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/layered_reconstruction_HCC1395_chr1_22_display_all.json
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/layered_region_view_HCC1395_chr1_22_display_all.json
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_topology_summary.json
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_topology_regions.tsv.gz
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_topology_summary.receipt.json
```

#### 實際漏斗與樹數字

| 指標 | 數量 | 單位／分母 |
|---|---:|---|
| Candidate sSNV S | 79,687 | autosomal PASS loci |
| Primary unique loci | 36,384 | unique loci；45.6586% of S |
| Exact PS×HP units | 39,544 | read-linked components |
| Bounded blocks | 39,807 | k≤12 blocks |
| k=1，不進樹 | 28,124 | blocks |
| Pattern support 不足 | 141 | blocks |
| Tree-input route-blocks | 11,542 | 28.9949% of all blocks |
| Primary mutation-bearing | 9,600 | T/Topo 主分母 |
| Reference-only controls | 1,942 | 不進 primary T/Topo 分母 |

Constraint weight 守恆：`285,897 = 281,967 retained + 1,254 cut + 2,676 unavoidable`；retained=98.6254%。

| Primary route-block 結構狀態 | n | 占 9,600 primary | 占 9,180 complete |
|---|---:|---:|---:|
| T=1、Topo=1 | 4,579 | 47.6979% | 49.8802% |
| T＞1、Topo=1 | 2,208 | 23.0000% | 24.0523% |
| T＞1、Topo＞1 | 2,393 | 24.9271% | 26.0675% |
| capped／incomplete | 420 | 4.3750% | — |

嚴格只分類 Topo=1 的 6,787 個：無分枝/無多層 3,239、sister-only 442、direct-only 2,977、sister+direct 129。另有 136 個 Topo＞1 的所有候選粗形態都為 sister+direct；若使用 `candidate-invariant coarse morphology`，可分類 6,923／9,600，但不得說這 136 個 exact Topo 已唯一。

#### QA 結果

```text
CMake exact_ps_partition: built target, exit 0
Topic regression: 59 passed in 28.55s
Repository CTest: 251/251 passed in 8.54s
py_compile: exit 0
git diff --check: exit 0, no output
Final receipt: 9/9 input/output/generator SHA identities PASS
```

#### [狀態更新] 可主張與仍未完成

- 已完成：exact PS segmentation → full molecule partial-pattern adapter → local Steiner candidate enumeration → T/Topo census → region-level machine table。
- 已完成：HCC1395 chr1–22 技術驗證；不是只做 chr22 smoke。
- 尚未完成：count-weighted read-pattern likelihood、caller VAF join/ranking、CN 與 methylation 輔助確認、跨 PS signed SAME/FLIP bridge、legacy W matched candidate-forest transition、其餘 6 樣本。
- `direct` 包含 hidden H_* nodes，不能直接等同 confirmed subclone；local mutation-state tree 不能直接等同 cellular clone tree。
- 上游 segmentation receipt 仍是 `exploratory_pilot / PARTIAL / validation_evidence_eligible=false`。本輪技術 checks PASS，但仍不具 cohort/paper final evidence 資格。
- 正式觀察只採 receipt 列出的 underscore display-all paths；同目錄較早、相似名稱的輸出視為 superseded。
- 本輪尚未 commit，等待使用者檢視結果後決定。
