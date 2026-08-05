<!--
建立時間: 2026-07-12 10:00
目標: HCC1395 兩技術流程粗拓撲與基因／藥物整合驗證的 living document
處理範圍: chr1-22、7 dataset rows、region×HP topology、gene/drug annotations
cycle_id: 20260712-1000-hcc1395-pair-coarse-topology-validation
spec_id: hcc1395_pair_coarse_topology_gene_drug_validation
status: in_progress
advisory: on
關聯檔案:
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/pre-decision-audit.md
-->

# Implementation Notes：HCC1395 兩流程粗拓撲與基因／藥物驗證

> 本檔記錄本輪分析的設計決定、偏離、折衷與未決問題；最終數字以可重跑輸出為準。

## 🔵 設計決定

### [2026-07-12 10:00] 五類粗拓撲為使用者鎖定的主要分類

- **Status**：Accepted。
- **背景**：使用者指定 `無 HP 內關係／姐妹 only／直系 only／姐妹＋直系／Topo>1 未定`。
- **決定**：主輸出在 `complete primary region` grain 五選一；每個 HP component 先保留 ordered class/shape signature，雙 HP region 再以 sister/direct feature OR 折疊為五類之一。跨 HP 永不建立姐妹／直系關係。
- **影響**：所有全 dataset 統計、pair concordance、基因圖與表格使用同一 enum。
- **Revisit if**：observed-state relation 無法互斥守恆或 Topo>1 候選其實全落同一 coarse class。
- **Evidence tier**：L5（user-specified definition，待 classifier 驗證）。

<!-- BEGIN USER-SPECIFIED -->
**Decision**：粗拓撲只使用五類：無 HP 內關係、姐妹 only、直系 only、姐妹＋直系、Topo>1 未定。
**DO NOT change**：未經使用者確認不可另增第六個 reader-facing 主類。
**Rationale**：2026-07-12 使用者明確指定。
<!-- END USER-SPECIFIED -->

### [2026-07-12 10:00] HCC1395 pair 視為技術重現性而非獨立生物重現

- **Status**：Accepted。
- **決定**：We will call HCC1395 and HCC1395_DORADO two technical processing datasets from one biological sample, not two independent patients。
- **理由**：避免 pseudoreplication；一致性只能支持 basecalling／pipeline robustness。
- **Revisit if**：provenance 顯示不同 aliquot、culture passage 或獨立 library。
- **Evidence tier**：L2。

### [2026-07-12 10:00] 同時報 exact 與 overlap-matched region agreement

- **Status**：Accepted。
- **決定**：We will use exact `(chrom,start,end)` matches as strict comparison and reciprocal-overlap／locus mapping as sensitivity analysis。
- **理由**：只用全域比例會漏掉 region-level disagreement；只用 exact boundary 會把 boundary drift 當 topology drift。
- **Revisit if**：source region keys 本身以 shared anchor ID 定義。
- **Evidence tier**：L3。

## 🟠 偏離之處

### [2026-07-12 10:00] 最新 clean-v3 未確認前不覆寫 historical snapshot

- **Status**：Accepted。
- **規範要求**：使用最新完整資料。
- **實作偏離**：若 clean-v3 aggregate verification 尚未 PASS，保留 historical layered-v2 作 engineering baseline，並將 clean-v3 可用子集獨立標示，不混算。
- **風險評估**：無法給 final biological rates，但可回答資料版本差異與方法重現性。
- **回退路徑**：clean-v3 PASS 後用同一 classifier 全量重跑。
- **Revisit if**：找到 terminal `_SUCCESS`、receipt 與 verify PASS。
- **Evidence tier**：L1。

## 🟡 折衷考量

### [2026-07-12 10:00] 癌症基因／藥物資料作外部合理性支持，不作 ground truth

- **Status**：Accepted。
- **方案 A**：We will test annotation-stratified concordance against matched genomic background and label it biological plausibility support。
- **方案 B**：把 cancer/drug overlap 當 topology truth — rejected，因註記不含此樣本的真樹。
- **方案 C**：完全不整合 annotation — rejected，因使用者要求 gene-region view。
- **Revisit if**：取得 HCC1395 single-cell、multi-region 或 validated clonal tree truth。
- **Evidence tier**：L2。

## 🔴 未決問題

### [2026-07-12 10:00] clean layered-v3 是否已在 7/12 完成

- **Status**：open。
- **Question**：producer、canonical/sensitivity tree、comparison、immutability verification 是否已有 terminal PASS？
- **Default if no answer**：historical snapshot + `PARTIAL / SCIENTIFIC NO-GO`。
- **Priority**：critical。
- **Evidence tier**：L5。

### [2026-07-12 10:00] 癌症基因與藥物 annotation 的 canonical source

- **Status**：open。
- **Question**：repo 既有哪個版本／來源可作本輪 SoT，join key 是 gene symbol、Ensembl 還是 interval？
- **Default if no answer**：只用可追溯本地資料，逐來源標版本；不臨時拼接不明資料。
- **Priority**：critical。
- **Evidence tier**：L5。

### [2026-07-12 00:48] clean-v3 point-in-time gate 已確認

- **Status**：resolved。
- **Evidence**：producer 6/7 terminal PASS；HCC1954 active；producer aggregate `_SUCCESS`／verification absent；canonical 與 sensitivity layered-v3 roots absent。
- **Decision**：本輪 quantitative SoT 保留 historical layered-v2；報告凍結時需再讀 marker，但即使 producer 7/7，也要等兩套 layered-v3 `_SUCCESS` 與 verification 才能替換。
- **Source**：`agents/latest_data_claim_audit.md`、`data/latest_pipeline_status.json`。
- **Evidence tier**：L1。

### [2026-07-12 00:50] 基因／藥物 canonical sources 已確認

- **Status**：resolved。
- **Primary gene interval**：GENCODE v46 basic gene body（GRCh38；2024-03-26）。
- **Cancer context**：COSMIC Cancer Gene Census v104，768 unique genes；gene-level membership，不是 variant／topology truth。
- **Drug context**：local DGIdb interaction snapshot；分成 any interaction、approved、approved∩anti_neoplastic。來源版本混合且 indication 不固定，不能稱 actionable therapy。
- **Sample-specific sensitivity**：COSMIC CLP v104 HCC1395（COSS749712）1,141 unique chr1–22 genomic alleles，其中 33 confirmed somatic；仍可能受 assay／passage／catalogue 差異影響。
- **Evidence**：`agents/source_inventory.tsv`、`agents/gene_drug_source_profile.json`。
- **Evidence tier**：L2。

## ✅ 實作結果更新

### [2026-07-12 00:40] 五類與 HCC1395 pair canonical 重算完成

- **Inputs**：historical `exact_topology_catalog.json`、`c_t_topology_report_data.json`、`vaf_top_tie_regions.tsv`。
- **Command**：`python3 research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/scripts/topology_pair_analysis.py [...]`；完整參數與 stdout 在 `agents/topology_pair_analysis.md`。
- **Outputs**：`data/topology_pair_analysis.json` 與 9 個 TSV。
- **Actual**：47,377 primary regions；exact matched=6,252；complete-both=5,720；agreement=69.3881%；κ=0.4973；chr-preserving null=39.5115%，p=1/5001。
- **Sensitivity**：四個 overlap/anchor scenarios agreement=69.22%–69.43%、κ=0.492–0.496。
- **QA**：64/64 PASS；獨立 `/tmp` rerun 的 9 TSV byte-identical。
- **Interpretation**：partial/moderate reproducibility；不是完全一致。

### [2026-07-12 00:50] 嚴格結構重現率揭露解析度落差

- **Resolved/unresolved binary**：agreement=73.15%、κ=0.4587；R/R=1,677、R/U=1,096、U/R=440、U/U=2,507。
- **條件式 coarse class**：兩邊都 resolved 後 agreement=87.18%，不可當全體 accuracy。
- **Ordered shape-set**：30.65%；phase-swap tolerant=47.22%。
- **Exact tree-set digest**：19.74%；phase-swap tolerant=35.21%。
- **Aggregate shift**：HCC vs Dorado 五類 TV distance=20.03%。
- **Decision**：方法抓到可重現 coarse signal，但 exact T 及 determinacy 依賴技術資料，不可互換。

### [2026-07-12 00:55] Annotation matched-null 完成

- **Input**：`agents/hcc1395_exact_complete_pair_gene_drug_flags.tsv`，exact complete-both n=5,720。
- **Command**：`python3 scripts/analyze_annotation_reproducibility.py --input ... --output-tsv data/hcc1395_annotation_reproducibility.tsv --output-json data/hcc1395_annotation_reproducibility.json`。
- **Null**：seed=20260712；5,000 次；chromosome＋global region-length decile conditional hypergeometric permutation。
- **CGC body**：n=268，agreement=72.39%，present−absent=+3.15 pp，p=0.5855。
- **DGIdb approved∩anti-neoplastic body**：n=454，agreement=72.25%，+3.11 pp，p=0.4295。
- **COSMIC CLP all-status region**：n=333，agreement=69.97%，p=0.3195。
- **CLP confirmed somatic**：n=12，agreement=75.0%，95% CI 46.8%–91.1%，p=0.7313。
- **QA**：TSV/JSON 第二次 rerun byte-identical；所有特徵 p>0.05。
- **Decision**：gene/drug context 沒有提供超越可比背景的 topology validation；face validity only。

### [2026-07-12 00:56] Notebook fallback QA

- **Output**：`notebooks/20260712_HCC1395_coarse_topology_gene_drug_validation.ipynb`。
- **Environment gap**：主 Python 與 repo `.venv` 皆無 `jupyter`／`nbformat`／`nbclient`／`ipykernel`，未安裝新依賴。
- **Fallback verification**：JSON nbformat=4、10 cells；依 cell 順序以共享 namespace 執行所有 code cells，exit 0，重現 64/64、3,969/5,720 與所有 annotation p>0.05。
- **Limitation**：notebook 尚未嵌入 executed outputs；可在具 Jupyter 環境用 nbconvert 重跑。

## 📚 Lore

### [2026-07-12] 三個不可互換的 grain

- **Constraint**：candidate-T occurrence、unit support、region support 不可互換；雙 HP 是 ordered forest。
- **Why it matters**：一致性分析必須固定在 region×HP 或 region grain，不能用 candidate count 當 prevalence。
- **Evidence**：`research/20260711_read_group_C_tree_T_topology_report/data/`。

### [2026-07-12] VAF 只能排序候選 T

- **Constraint**：HP-specific read-AF 受 CN、purity、coverage 影響。
- **Why it matters**：即使兩流程選到同一 top T，也只能稱 VAF-supported agreement，不能稱真樹確認。
- **Evidence**：2026-07-11 topology report。

### [2026-07-12 02:18] Portable HTML 與瀏覽器 QA 完成

- **Inputs**：`InterSubMod/docs/reports/in_progress/2026/07/20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01/artifact.json`。
- **Build**：canonical builder 輸出 47 blocks、6 charts、18 tables、4 metrics；topology 64/64 與 confusion/metrics n=5,720 守恆 PASS。
- **Compatibility**：runtime sticky header 的 `100vw` 在 classic scrollbar 環境溢位；report-local CSS 改為 `width:100%`。390 px chart legend 僅縮短顯示標籤，完整五類名稱與數據未變。
- **Final browser verify**：unwrapped chrome-headless-shell；1440/390 px 均 PASS；source dialog 與 keyboard interaction PASS；external requests、console/page errors、global horizontal overflow 皆為 0。
- **Visual receipts**：`qa_desktop.png`、`qa_mobile.png`；舊 `qa_failure.png`／`qa_packaging_failure.png` 保留為 superseded audit evidence。
- **Receipt**：`InterSubMod/docs/reports/in_progress/2026/07/20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01/qa_receipt.json`。
- **Final HTML SHA-256**：`1e4432ccfa1b4022505c35dc995a8ae5b3576b3906f2eeb5877f411a19922dd6`；document language=`zh-TW`。
- **Verdict unchanged**：SHARE WITH CAVEATS / PARTIAL technical reproducibility；SCIENTIFIC NO-GO for proof-of-effectiveness。

### [2026-07-12 02:23] Final-freeze pipeline 狀態更新

- **Producer**：HCC1954 於 01:31 完成，raw-all producer 現為 **7/7 PASS**；aggregate `_SUCCESS` 與 producer verification present。
- **Closeout**：latest retry 仍 **FAILED**，error=`E_SIDECAR_VALIDATION`；continuation execution=`FAILED`（exit 4）。
- **Clean tree outputs**：canonical／sensitivity roots仍 absent，故沒有 clean-v3 topology 可替換 historical-v2。
- **Report action**：只更新 point-in-time pipeline status；所有五分類／HCC pair／gene-drug quantitative rates維持 historical-v2，避免混算。
- **Final QA**：browser 1440/390 PASS；delivery validator 49/49 PASS；獨立唯讀 auditor PASS。

### [2026-07-12 03:01] HCC1395 pair 的 VAF-supported 推測樹比對已補入

- **使用者要求**：確認兩個 HCC1395 dataset 的樹結構比對確實有使用每個位點 VAF/read-AF 的可能性結果，且讀者知道它是推測。
- **Inputs**：`vaf_top_tie_units.tsv`、`vaf_top_tie_candidates.tsv.gz`、`vaf_top_tie_regions.tsv`、historical layered-v2 `mlhp_part_*.json`、`data/hcc1395_pair_matches.tsv`。
- **Command**：`python3 research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/vaf_pair_comparison.py --vaf-units ... --vaf-candidates ... --vaf-regions ... --pair-matches ... --pair-analysis ... --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 --output-dir research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data`。
- **Outputs**：`data/hcc1395_vaf_pair_regions.tsv`（5,720 rows）、`data/hcc1395_vaf_pair_metrics.tsv`、`data/hcc1395_vaf_pair_summary.json`、`data/hcc1395_vaf_pair_checks.tsv`。
- **Primary selection**：same-HP read-AF ordering 使用整數 read counts 轉 `Fraction` 的 exact score argmax；所有 exact ties 保留。Softmax 只作 candidate-set 濃度敏感度，不是 posterior。
- **Shape conservation**：HCC1395 6,798=3,438 structural Topo=1＋3,360 VAF-selected；DORADO 6,082=2,444＋3,638。原 Topo>1 必須所有 ambiguous units 可評估才允許 selected shape，採 fail-closed。
- **VAF shape-rescue pair**：兩側原始都 Topo>1、各自被 VAF 縮至單一 shape，n=2,089；ordered=919/2,089=43.99%，phase-swap=1,624/2,089=77.74%。這是 rooted-unlabeled shape，不能分辨 mutation direction。
- **VAF unique exact forest pair**：兩側都實際使用 VAF 且各自唯一，n=2,543；ordered=524/2,543=20.61%，phase-swap=949/2,543=37.32%。Signature 包含每個 HP 的 genomic positions＋candidate IDs，但仍是同批 reads 的 heuristic inference。
- **QA**：20/20 PASS；fixed-T1 shape components 6,554/6,554、candidate-detail join 10,446/10,446；3 TSV rerun byte-identical，JSON 排除 output path 後 semantic-identical。
- **Report action**：保留未排名 exact candidate-set 19.74%/35.21% 作 baseline，另加 VAF 推測層；HTML 49 blocks、6 charts、19 tables、4 metrics，1440/390 browser PASS；delivery validator 66/66 PASS。
- **Claim ceiling**：VAF-selected shape/exact forest 都不是 accuracy、calibrated posterior、independent validation 或 biological clone truth；整體仍為 PARTIAL / SCIENTIFIC NO-GO。

## Provenance Footer

- **Commit hash**：`6067568637088838a9f518955e41d222f057e4f1`
- **Build time**：2026-07-12 10:00:00 +08:00
- **Skill version**：`/implementation-notes` v0.1
- **Line count**：本檔應維持 ≤400 行。
