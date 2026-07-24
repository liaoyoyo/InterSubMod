<!--
建立時間: 2026-07-24
目標: 記錄 exact-PS 新權威接入 layered_workstation 的設計、偏離、折衷與未決
處理範圍: index + 7 sample HTML / data contract / browser QA
關聯檔案:
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/pre-decision-audit.md
  - InterSubMod/docs/methodology/_assets/layered_workstation/
-->

# Exact-PS layered workstation implementation notes

## 設計決定

- [2026-07-24] 使用者已決：`20260724 exact PS × primary HP × strict
  read-linkage threshold=3` 升為 `layered_workstation/` 預設 authority。
- [2026-07-24] 預設分母分成三層：98,955 final groups、85,941
  mutation-bearing groups、71,955 `ranked_complete` groups；圖表不得跨分母混算。
- [2026-07-24] 精確 determinacy 採 signature census：
  `UNIQUE_TREE`、`TIED_SAME_TOPOLOGY`、`TIED_CROSS_TOPOLOGY`。
  代表樹只是 detail exemplar，不代替完整 tie census。
- [2026-07-24] `technical_all_pass` 與 `topology-complete` 分開：
  10,717 `ABSTAIN_RESOURCE_LIMIT` 必須在主畫面可見。
- [2026-07-24] GRCh38 chr1–22、搜尋、resolution/morphology/HP/status
  多選、再次點擊取消、無選擇即全部，列為互動 contract。
- [2026-07-24] JSON/receipt 路徑置於預設收合的 provenance 區，
  不插入主要數字敘事。
- [2026-07-25] Region detail 的資訊鏈固定為：
  `locus 排列 → read 狀態/連結限制 → minimum candidate family →
  AF 排序 → global-best tie census → selected exemplar`；不得只畫代表樹。
- [2026-07-25] 候選空間以 factorized sidecar 保存 minimum vertex sets 與
  每個 child 的 best-parent choices，瀏覽器需要時再展開單棵樹；
  避免把同一候選樹重複展開成近百萬筆 edge list。
- [2026-07-25] cohort 相似度用同分母內的分布 profile 比較，至少拆開
  status、resolution、morphology、active-k；若顯示 overview，
  必須標明等權彙總與不確定區間，且不得解讀為相同 clone tree。
- [2026-07-25] HCC1395／HCC1395_DORADO 另設 technical-pair evidence：
  locus/component overlap 與共同 exact units 的 resolution、signature、
  selected-tree concordance分層呈現；HP 標籤不可任意配對。
- [2026-07-25] Cohort presentation aliases 固定為
  `HCC1395_HKU`／`HCC1395_NYGC`，比較順序固定為
  `HKU→NYGC→HCC1937→HCC1954→H1437→H2009→COLO829`；authority
  keys、檔名、receipt 與 hash 仍使用 `HCC1395`／`HCC1395_DORADO`。
- [2026-07-25] 圖上 hidden state 採 `H:S1+S2` compact label；
  完整 `H_ARR...` state 僅放 hover、vertex-set 與 edge census，避免
  高 k 樹節點文字重疊。
- [2026-07-25] 癌症基因／藥物只以 actual sSNV locus 落入 GENCODE v46
  gene body 為位置 authority，CGC v104 與 DGIdb 以同一 HGNC gene
  交集；region span 不算 hit。`NO_HIT_EVALUATED` 與 `HIT` 必須分開。
- [2026-07-25] topology legend 與 gene/drug quick filter 採 AND；六個
  genome modes 都維持「多選聯集、第二次取消、未選即全部」。切 filter
  後不存在的 legend key 必須 prune，避免孤兒 selection。

## 偏離

- 舊 v5 `adjacent_gap<=50000` region keys 不可與 exact-PS region keys
  假裝一對一；因此不沿用舊逐區域 topology/annotation cache。
- 舊 92.18% topology-unique 指標不沿用；新 exact rooted-unlabeled
  topology 口徑為 88.2579%。
- 新 cohort 技術 PASS，但不是 topology-complete；HTML 不使用
  「全部已解析」「唯一 clone tree」等語句。
- 第一版 exact-PS adapter 只帶入 representative vertices/edges 與
  signature counts；上游 census 計數後丟棄 minimum candidate family，
  導致新頁面看不到候選集合、locus rail、read-state matrix 與替代樹。
  這是資料契約缺口，不是單純 CSS 收合問題。

## 折衷

- 每個 sample 頁嵌入完整 compact group index；detail 只保留重建展示所需欄位，
  原始完整 JSON 仍由收合 provenance 提供。
- topology signature 完整集合只存在 71,955 ranked rows；其餘群組保留
  `read_af_status`/ABSTAIN 原因，不臆造 resolution。
- HCC1395 與 HCC1395_DORADO 同時展示，但比較標籤固定為 technical
  concordance，不稱 2 個 biological samples。
- 舊頁面的 50 kb baseline 僅可作歷史說明，不得參與新 cohort 圖表。
- cohort overview 採「分布相似度」而非把不同分母的 raw counts 拼成一條
  cosine vector；各維度先正規化，再以 Jensen–Shannon profile similarity
  比較，overview 等權，避免大分母狀態表支配結論。

## 未決

- DGIdb 本地 export 沒有完整 release-level provenance，只能顯示逐列
  source name/version 的 gene-level context；不得升格成 clinical
  actionability。
- 10,717 resource-limit ABSTAIN 仍繼承上游 claim ceiling；本輪不補造
  candidate family。

## 執行紀錄

### 2026-07-24｜前置 authority audit

- Census root：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1/`
- Topology root：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/`
- `receipt.v2.json`：7 sample exit status 0、71,955 rows canonical
  reproduction PASS、680,527 best trees。
- `cohort_receipt.json`：7/7 technical PASS；但
  `all_mutation_bearing_families_complete=false`。

### 2026-07-25｜layered workstation authority promotion

- Default driver：
  `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`
  已改為委派 exact-PS builder；`--legacy-v5` 才重現舊 50-kb profile。
- Primary builder：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py`。
- 實際輸出：
  `InterSubMod/docs/methodology/_assets/layered_workstation/index.html` 與 7 個
  sample standalone HTML。
- Fail-closed 驗證：98,955 MLHP/topology 1:1、71,955 ranked/census 1:1、
  非 ranked 禁止接 census；schema、group index、region/unit/block/PS/HP、
  score/tie/signature/coarse tree counts 全部守恆。
- `chr10:87818272-87928739` regression：
  HCC1395 只允許 positions `[87818272,87840023]` 的
  `UNIQUE_TREE/Direct-only` unit；DORADO final topology overlap 必須為 0。
- Python regression：5 tests PASS。
- Chromium/Playwright：index + 7 sample pages × desktop/mobile，共 16
  page states；console error=0、external request=0、horizontal overflow=0，
  圖例聯集/再次取消、座標搜尋與 detail renderer 全 PASS。
- Browser receipt：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit/receipt.json`。
- Claude Code 唯讀審查：P0=0、可提交；指出 README 後半 legacy 口徑殘留。
  修正後再次複核為 `RESOLVED`、無 BLOCKER。

### 2026-07-25｜完整候選空間 sidecar v2

- Input：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/`
  與 frozen LongLineage C++ source/header。
- Command：
  `python3 InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_candidate_factorization.py
  --output-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260725_exact_ps_candidate_factorization/all7_v2`。
- Output：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260725_exact_ps_candidate_factorization/all7_v2/`。
- Receipt：schema `1.1.0`，SHA-256
  `54a8a0a00c1cdb1a40fe96b3a528e9142d21c76af55dc8eae553a6f1432b8164`；
  71,955 ranked units、972,592 minimum vertex sets／minimum trees、
  680,527 global-best trees，9 checks 全真。
- Claim ceiling：`validation_evidence_eligible=false` 與
  `all_mutation_bearing_families_complete=false` 明示保留。

### 2026-07-25｜cohort similarity sidecar

- Command：
  `python3 InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_cohort_similarity.py
  --output /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260725_exact_ps_cohort_similarity/all7_v1/cohort_similarity.v1.json`。
- Output：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260725_exact_ps_cohort_similarity/all7_v1/cohort_similarity.v1.json`。
- SHA-256：
  `926ff9bf636b1596a4a4123bcac10095c3ab4e20b2c3b0f7d4dbe504d8c40abc`。
- HCC technical pair：五維 profile similarity `0.926381`（21 pairs
  rank 1），chr-block bootstrap 95% CI `0.912684–0.934393`；
  strict selected labeled tree `443/564 = 78.5%`。頁面明示
  profile similar 不等於同一 clone/tree。

### 2026-07-25｜candidate/detail 與 cohort UI v2

- Command：
  `python3 InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py`。
- Output：`InterSubMod/docs/methodology/_assets/layered_workstation/index.html`
  與本機 7 個 standalone sample HTML。
- 實際輸出片段：
  `OK exact-PS authority verified: samples=7 groups=98955 ranked=71955`；
  elapsed `3:23.41`、max RSS `3,748,996 kB`。
- Region detail 已恢復 locus rail、read-state matrix、complete candidate
  edge union、global-best edge union、deterministic selected overlay、
  shape exemplar、vertex-set carousel 與 selected result；selected edges
  必須是兩層 union 的子集，tied global-best union 不得冒充一棵樹。
- HCC1395 目標 group：2 loci、5 solver-visible pattern rows、3 minimum
  sets、4 union edges、2 selected edges；S3/S4 明示 singleton ABSTAIN。
- Python regression：
  `python3 -m unittest ...test_exact_ps_candidate_factorization.py
  ...test_exact_ps_cohort_similarity.py ...test_exact_ps_layered_workstation.py`
  → `Ran 22 tests in 64.066s`、`OK`。

### 2026-07-25｜UI v2 Chromium baseline（後由 UI v3 final 取代）

- Command：
  `python3 InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py
  --output-dir InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v2_all7_05`。
- Output receipt：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v2_all7_05/receipt.json`，
  SHA-256
  `31f391efa33efd6e22e4914bf0bff3ad1db665adce1b338827377271a05cfdd2`。
- 16 page states、44 PNG；`all_pages_audited`、`no_console_errors`、
  `no_external_requests`、`no_horizontal_overflow`、
  `screenshots_recorded` 全真。桌面圖不需內層捲動；手機 overflow
  必須可實際捲到最右端。
- Claude Code 2.1.216／Haiku 唯讀結論：
  `AGREE WITH NOTES`，P0/P1=0。3 個 P2（手機可達性、截圖 hash/尺寸、
  alias 註解）已全部落地。

### 2026-07-25｜exact-locus 癌症基因／藥物 sidecar 與 UI v3

- Input：
  `/big7_disk/liaoyoyo2001/gene_annotation/gencode.v46.basic.annotation.gtf.gz`、
  `Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz`、
  `cosmic_v104/Cosmic_Genes_v104_GRCh38.tsv.gz`、
  `dgidb_interactions.tsv`。
- Command：
  `python3 InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_gene_drug_annotation.py
  --gencode /big7_disk/liaoyoyo2001/gene_annotation/gencode.v46.basic.annotation.gtf.gz
  --cgc /big7_disk/liaoyoyo2001/gene_annotation/Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz
  --cosmic-genes /big7_disk/liaoyoyo2001/gene_annotation/cosmic_v104/Cosmic_Genes_v104_GRCh38.tsv.gz
  --dgidb /big7_disk/liaoyoyo2001/gene_annotation/dgidb_interactions.tsv
  --output-dir /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260725_exact_ps_gene_drug_annotation/all7_v2`。
- Output：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260725_exact_ps_gene_drug_annotation/all7_v2/annotation.v1.json`
  （859,706 bytes；SHA-256
  `148f74e9f958cccffe14a2297e9453e28b2346a424adcb30e350f79e91ae4f50`）
  與 `receipt.json`（3,834 bytes；SHA-256
  `9639b36df7a2412afc2aca97eb1f20819e7a41394394ef40b6f3108e809dd99e`）。
- Producer contract：721 CGC genes、279 drug-linked；224 non-HGNC DGIdb
  rows、GENCODE gate accepted=0、CGC overlap key=0。Receipt 的 10 checks
  全由 payload/count/source identities 推導。
- HTML all7 actual-locus結果：CGC `3,554`、同基因 CGC∩drug `1,252`；
  HCC1395 target 為 `NO_HIT_EVALUATED`，不得以舊 span 誤稱 PTEN hit。
- `--verify-only` 實際輸出：
  `OK exact-PS authority verified: samples=7 groups=98955 ranked=71955`，
  exit 0、filesystem outputs=0。
- 正式 build：elapsed `3:41.93`、max RSS `3,775,812 kB`、exit 0；
  index `80,644` bytes，7 個 sample HTML 約 14–182 MB。
- Python regression：gene/drug 11、candidate 9、similarity 8、
  workstation 6，共 34 tests 全部 PASS。

### 2026-07-25｜Opus 5 複核與 6-mode Chromium final

- Claude Code 2.1.219 實際 `canonicalModel=claude-opus-5`、
  `--effort max`、`--permission-mode plan`、唯讀 `Read/Grep/Glob`。
- 初審：`AGREE_WITH_NOTES`、P0=0；指出 4 個 P1：non-HGNC DGIdb
  contamination、symbol fallback、literal receipt checks、filter-switch
  orphan selection。四項修正後複核全部關閉，P0/P1=0，明示
  「合理，可提交」。
- Opus 剩餘 P2 中，「多選只測 resolution」已立即補成六個 genome
  modes 全測；DGIdb count key 名稱與 AND browser expected 同源則記為
  非阻塞 maintenance note。
- Final command：
  `python3 InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py
  --output-dir InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v3_all7_07`。
- Final output：
  `receipt.json` SHA-256
  `0b76be6b84e5007e9356b422f78157068bc96111192ba0dc001f4d9b84c4385c`；
  16 page states、44 PNG、84 mode records，所有 checks 全真，Chromium
  `148.0.7778.96` binary SHA
  `adc1c21ceed5c2a67184766376fe816ac03e556cc0ca3f782e8212235fe05c6f`。
