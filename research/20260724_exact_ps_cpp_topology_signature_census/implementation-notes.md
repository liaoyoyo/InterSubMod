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

## 偏離

- 舊 v5 `adjacent_gap<=50000` region keys 不可與 exact-PS region keys
  假裝一對一；因此不沿用舊逐區域 topology/annotation cache。
- 舊 92.18% topology-unique 指標不沿用；新 exact rooted-unlabeled
  topology 口徑為 88.2579%。
- 新 cohort 技術 PASS，但不是 topology-complete；HTML 不使用
  「全部已解析」「唯一 clone tree」等語句。

## 折衷

- 每個 sample 頁嵌入完整 compact group index；detail 只保留重建展示所需欄位，
  原始完整 JSON 仍由收合 provenance 提供。
- topology signature 完整集合只存在 71,955 ranked rows；其餘群組保留
  `read_af_status`/ABSTAIN 原因，不臆造 resolution。
- HCC1395 與 HCC1395_DORADO 同時展示，但比較標籤固定為 technical
  concordance，不稱 2 個 biological samples。
- 舊頁面的 50 kb baseline 僅可作歷史說明，不得參與新 cohort 圖表。

## 未決

- 全量頁面在 Chromium 的初次 parse/filter latency，待真實 browser audit。
- `chr10:87818272-87928739` 在新 exact-PS authority 下的拆分與可檢視
  region keys，待逐區域 trace 後寫入說明。
- 舊 cancer-gene annotation 若無可證座標 join，本輪維持 out-of-scope，
  不因展示需求犧牲 authority。

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
