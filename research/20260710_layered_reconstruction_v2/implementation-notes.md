<!--
建立時間: 2026-07-10 22:44
目標: layered reconstruction v2 實作過程 living document
處理範圍: core schema、runner、7 dataset outputs、verification、reports
cycle_id: 20260710-2244-layered-reconstruction-v2
spec_id: layered-reconstruction-v2
status: in_progress
advisory: on
build_commit: 4fb9e742482b63a660de19a1f1bd07d49d713111
-->

# Implementation Notes: Layered Reconstruction v2

> **Purpose**：記錄從旁路報告修正升級為 canonical pipeline 修正的設計決定與限制。

## 🔵 設計決定

### [2026-07-11 16:24] 全景 HTML 格式／數據 hotfix 與隱藏 JSON 來源

- **Task type / goals**: E bugfix；服務 G3/G5；scope 為整份 standalone HTML 的展示、current/historical 語意與可重算數字，不修改 upstream research artifacts。
- **User decision**: 各資料 `.json` 不得干擾數字閱讀；改為預設收合、只在人類可讀「資料與驗證／來源」連結的 `href` 中保留。
- **Root cause**: active normalized raw-all、aborted PASS-only 與 historical 7-row snapshot 缺 evidence-track 分層；同時 5,657 個 inline source tooltip 讓 provenance 標記壓過數字，historical 14-chart atlas 與精確表又全數預設展開。
- **Inputs**: historical root `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/`；active root `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/`；builder `InterSubMod/research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py`。
- **Build command**: `python3 research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 --production-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2 --output-dir docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01 --shell-template research/20260710_layered_reconstruction_v2/report_runtime/html-report-shell.v0.2.6.html --embed-helper research/20260710_layered_reconstruction_v2/report_runtime/embed_html_report_runtime.v0.2.6.py --hcc-hc-bed /big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed --hcc-cn-bed /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`；exit `0`。
- **Renderer custody**: 16:01 plugin cache 從 0.2.6 輪替到 0.2.8，原 ephemeral shell/helper 路徑消失；已依既有 hash 與 final artifact 反向驗證後，把 exact shell、compressed runtime source/style、source-tooltip、embed helper、chart contract 與 `package_utils.py` 固定於 `InterSubMod/research/20260710_layered_reconstruction_v2/report_runtime/`。重生不再依賴 plugin cache。
- **Outputs**: `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html`；companions；`audit/after/metrics.json`；`audit/after/data_quality.json`；desktop/mobile screenshots。
- **Regression evidence**: data contract `70/70 PASS`；Chromium 147 desktop/mobile `37/37 PASS`；console/page/request errors=`0/0/0`；body overflow=`0/0 px`；inline tooltip=`0`；5 個 JSON href target 全存在且 link label 不顯示 `.json`。
- **Readability delta**: HTML 2,914,792→985,978 bytes；DOM 36,741→19,716；desktop 22,034→9,142 px；mobile 30,098→14,434 px；historical charts、datasets、精確表皆改為 progressive disclosure。
- **Data verdict**: frozen aggregate仍為 `U=568,080`、`retained=182,400`、`W_tree=48,959`、`W_primary=47,377`、complete/incomplete=`39,885/7,492`；active producer於 16:23 為 `3/7 PASS`、`H1437 START`、aggregate `_SUCCESS` absent。Artifact PASS；Scientific NO-GO。
- **Evidence hashes**: HTML=`069980a2de9c236aa075c6a17ab1150ce0df16a17856a3f39a411c89b5e9ffe6`；data=`2a22585776e23b66af1e80169540bb50df28e33382255e59132fd33d808cc4d8`；builder=`250a41d22287b6605361b3e677da4ec766e6143fddbc18a62ba4296318fab0e7`；embed helper=`fb76020ad214d34d094695db7584a90e01a9c7e9e32f54d353a3359dc34dcd5d`；browser QA=`b62bde8ae02c66311684eb7284e24fa9af2f273937679039b1e52c295c328aa9`；data QA=`0a9bb8ad46b85608ec5928e30751880826f898ca2d4c68f0359e34f5391cb70b`。

### [2026-07-11 07:13] 教授版 S→W→HP→C→Topo 報告完成三視角驗證

- **Status**: Internal PI-share ready；artifact validation PASS；scientific release 仍為 `NO_GO`，輸出只留在 `InterSubMod/docs/reports/in_progress/`。
- **服務目標**: G3（read/HP/C/Topo/hidden 與 auxiliary role 可解釋）、G4（7 datasets 同口徑守恆與可重算）、G5（來源、分母、圖表 fallback、QA 與 claim ceiling 可外部稽核）。
- **形式定義**:
  - 單一 primary HP：`C_region=n_trees`；HP1+HP2：`C_region=n_trees(HP1)×n_trees(HP2)`。
  - `Topo_region=∏n_distinct_shapes_exact`；完整區域必滿足 `C≥Topo≥1`。
  - H3? presence 只依 canonical `n_H3_auxiliary>0`；reference-only、`mutation_bearing=false` 的 family=3 不得算 H3+。
  - 可行狀態只有 `C=1/Topo=1`、`C>1/Topo=1`、`C>1/Topo>1`；`C=1/Topo>1` 數學上不可能。
  - `n_hidden,region=Σn_hidden(primary HP)` 只在候選全集完整時使用；hidden 是 parsimony-inferred mutation-state node，不是 hidden clone。
  - most-likely 與 candidate-tree read-AF 本輪為 `NOT EVALUATED`；0/7 artifacts、0 ranked regions，不能選 winner/posterior。
- **獨立重算**: 7-dataset row aggregate `W_tree=48,959`、`W_primary=47,377`、complete=`39,885`、incomplete=`7,492`；三個 complete topology states 分別 `10,832 / 11,144 / 17,909`，impossible=`0`；hidden `0 / >0 = 7,325 / 32,560`；HP×H3 六格=`23,167 / 3,043 / 16,775 / 4,392 / 459 / 1,123`。C=1…6/>6/incomplete 與 site/region funnel 逐 dataset 均守恆。
- **抓到的既有 label 例外**: H2009 有 3 個 recurrence+capped 混合區被舊 `region_determinacy` 優先標為 recurrence；本報告直接依 primary-unit flags 計算，所以 incomplete region=`3,699`，不是 legacy `has_capped=3,696`。
- **Production evidence 更新**: 原 `20260711_longphase_s_production_sidecars_v1` 已終止並改名 `20260711_longphase_s_production_sidecars_PASS_ONLY_ABORTED_v1`；`_FAILED` reason=`E_METHOD_SCOPE_PASS_ONLY`。PASS-only input 無法驗證 LongPhase-S native non-PASS→PASS rescue；normalized raw-all bounded probes其後已過，7-dataset raw-all rerun獲 `GO_WITH_FAIL_CLOSED` 授權但尚未啟動完成。
- **輸入路徑**:
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2`
  - `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_PASS_ONLY_ABORTED_v1`
  - `/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`
  - `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`
- **執行命令**: `python3 InterSubMod/research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 --production-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_PASS_ONLY_ABORTED_v1 --output-dir InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01 --shell-template /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.6-d37358633e00/assets/html-report-shell.html --embed-helper /bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.6-d37358633e00/skills/build-report/scripts/embed_html_report_runtime.py --hcc-hc-bed /big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed --hcc-cn-bed /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`。
- **主輸出**: `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html`；SHA-256=`2098e2493e621843b1e9692053c4396d17af0060a508abe4d52a75f18fa842af`；data SHA-256=`78cec57adff20d6d3c21acd29e75499acee3be220fdc27012a0029bcbe9d3042`。
- **報告內容**: 14 個 Recharts 圖 + 14 個 same-data no-JS fallback、2 個教學 SVG、完整雙單位 funnel、11 個摺疊名詞、7 個 dataset standalone details、HP×H3、single/double C、Topo 三類、hidden/topology 交叉與 most-likely N/A panel。
- **驗證結果**:
  - artifact hard gate PASS：官方 `chart_contract.py` PASS（errors/warnings=[]；1 stackedBar + 9 stackedBar100 + 2 grouped bar + 2 horizontalBar）、14 hosts/14 fallbacks、5,671 source marks（5,657 inline + 14 focusable chart controls）、duplicate IDs=0、unnamed figures=0、tables without caption=0、headers without scope=0、runtime marker=0。
  - browser QA PASS：1440 light/dark、768、390、360、320、no-JS desktop/mobile；14/14 live SVG、7/7 HP labels、first-tab skip link、0 source/title overlap、0 console/page errors。
  - 三個 subagent 分別由 funnel/HP×H3、C/Topo/hidden invariants、科學敘述與 HTML accessibility 視角獨立稽核；修正後所有 must checks 通過。
- **程式雜湊**: builder=`94fa61a500d28de94874c35c3cc6e182b454290c88f520a81c317e0aae60bbda`；browser QA=`288e6614c235b06b098853d0d13f469fe1fa9a059b5655993fc0da6205f26948`；official chart contract=`7d2f7cfdf3726bd5930970547c878d3809e455dffcadb57ac36a47637140fb86`。
- **Evidence tier**: L1 deterministic recomputation + L1 browser/artifact QA + independent reviewer agreement；生物主張仍受 historical input mismatch、production abort、read-AF/L3 unavailable 與無 single-cell truth 限制。

### [2026-07-11 06:18] 唯一 reviewed handoff supervisor 已進入等待態

- **Status**: Active wait only；full-run root與scientific workers均尚未建立。
- **協作確認**: 另一個同父 Codex AI於`2026-07-11 06:16:43 +0800`啟動PID `929656`；主代理的真實dry-run因同一handoff flock回`E_HANDOFF_LOCKED`，因此確認不是第二套可並行啟動的supervisor。
- **Reviewed plan**: `InterSubMod/research/20260710_layered_reconstruction_v2/launch_plans/20260711_layered_v3_handoff_launch_plan_reviewed_v1.json`；SHA-256=`7e0ab871ee2fa15e772d17da75ebe0e836cc7dcc6c2726e7b136a3b305e3da6a`，與unreviewed template除authorization三欄外semantic diff為空，permissions已freeze為`0444`。
- **Execution gate**: CLI `--execute-reviewed-plan` + JSON `execute=true` + out-of-band plan SHA三層成立；supervisor自身也在25個executable pins內，SHA=`cc91239368bb79732956f2e08fb608de9dea6dcb70d4d50bf724a1e9cf95e49f`。
- **主代理readback**: 5 inputs + 25 executable/source/schema/tool pins全部重算PASS；controlled child PATH固定`/usr/local/bin:/usr/bin`，samtools/bgzip/tabix/stat/bcftools皆解析至pin。Supervisor tests主代理重跑28/28 PASS（10.548s）。
- **目前行為**: 持有`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/.20260711_longphase_s_production_sidecars_v1.handoff.lock`，每30秒只讀exact authority；7/7與quiescent前不建handoff workspace、不產receipt、不啟動runner。
- **Evidence tier**: L1 code/plan/test + L2 live PID/flock/filesystem observation；production completion仍0/7。

### [2026-07-11 05:49] v3 runner 正式資源門檻與 TOCTOU 收斂

- **Status**: Component-ready；real producer 7/7與300秒baseline未過，full run仍NO-LAUNCH。
- **初次退回原因**: runner第一版預設為RAM 16 GiB、disk 50 GiB、iowait 0.25秒，低於red-team已核准門檻；process list亦未含LongPhase producer，因此不能接受「component無缺」。
- **Frozen production policy**: logical CPU≥8、MemAvailable≥128 GiB、big7 free≥500 GiB、load/CPU≤1.25、300秒iowait≤20%、`/big8_disk` NFS per-op READ bytes_recv <80 decimal MB/s；production inventory禁止CLI/env放寬。
- **排他性**: 固定`<run-parent>/.layered_chr1_22_7dataset_full.lock`跨run ID；start/end process snapshots皆拒絕LongPhase/sidecar producer與layered workers，封閉baseline期間新process TOCTOU。
- **證據**: 主代理獨立重跑 `env LC_ALL=C TZ=UTC PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 python3 InterSubMod/scripts/test_run_layered_v3.py`，12/12 PASS（9.854s）、exit 0；py_compile/diff-check PASS。Runner SHA-256=`39d27369bf1ccb5a61cd06a41217fe7e4d78ca89b370f1f8a7e3aa275d497b90`。
- **現場反證**: active `longphase-s`、capture與producer wrapper仍在；`nfsiostat 1 2 /big8_disk`第二窗為115,588.445 kB/s，因此目前preflight應fail，不可launch。
- **Evidence tier**: L1 adversarial mechanism + L2 active resource/process observation。

### [2026-07-11 05:45] v3 producer receipt 的 symlink inventory 語意修正

- **Status**: Accepted；真實 HCC1395_DORADO regression 已修，receipt schema 不變。
- **問題**: active producer 的 `stat -c` 未使用 `-L`，所以 `input_files.tsv` 保存 requested symlink 的 `lstat` metadata；normalizer 卻用 `Path.stat()` 比 target metadata。真實 normal BAM 為 symlink size `80`、target size `91,359,905,868`；tumor BAM為 `79` vs `250,150,882,482`，因此完成後會被錯誤拒絕。
- **決定**: inventory readback精確比對 logical path `lstat(size,mtime)`；另由 `storage_identity_v1` 同時鎖 requested path、`logical_is_symlink`、realpath、target dev/inode/size/mtime/ctime、三段chunk與BAI full SHA，不能只鎖link或只鎖target。
- **程式**: `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_longphase_production_capture_receipt_v3.py`；新增symlink regression於同目錄 `test_build_longphase_production_capture_receipt_v3.py`。
- **驗證**: normalizer `8/8 PASS`（0.375s）；frozen contract `17/17 PASS`（32.641s）；`py_compile`與`git diff --check`皆exit 0。Builder SHA-256=`9f3b7f11393091fec8dc4b0ab2b1b8ca2c5bfd9376443cedc0ae060a1b8ccd3a`。
- **Evidence tier**: L1 code/regression + L2 active path metadata readback；production receipt 7/7仍pending。

### [2026-07-11 05:08] Active incomplete producer 的真實 fail-closed readback

- **Status**: PASS as negative/adversarial evidence；production仍in progress。
- **執行**: 將active run root交給 `prepare_clean_layered_manifest_v3.py`，要求建立新manifest。
- **實際結果**: exit `3`、`E_INDEX_MISMATCH`；HCC1395 recalibrated-all VCF/index尚未產生。正式manifest path不存在，failure report明示 `valid_lock_written:false`。
- **Evidence**: `InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_active_incomplete_layered_v3_manifest.failure.json`。
- **判斷**: 此PASS只證preflight能阻止未完成producer；不增加production完成度，也不解鎖full run。

### [2026-07-11 05:05] 教授版全層數據觀察 HTML 與雙軌證據邊界

- **Status**: Internal PI-share ready；`PARTIAL SCIENTIFIC VALIDATION / PUBLICATION NO-GO`。
- **服務目標**: G3（read/HP/候選樹與 auxiliary role 可解釋）、G4（7 datasets 同口徑與 QA）、G5（來源、分母、限制與 self-contained artifact 可稽核）。
- **設計決定**: 報告明確拆成兩條 evidence track：
  1. `20260710_232501_layered_reconstruction_v2` 是 **hash-verified historical engineering snapshot**；7/7 output hashes PASS，但 6/7 HP/PS inputs 受 truth-BED conditioning，且 `code.sha256` 現況只有 2/7 entries match。
  2. 最新 canonical flow 為 raw ClairS all → ClairS PASS LongPhase-S input → `_sc.vcf` all → `_sc.vcf` PASS tree input；truth-unrestricted production-sidecar probe 在 2026-07-11T05:04:05+08:00 仍為 2 START、0/7 PASS、無 aggregate。
- **量化範圍**: 9 張圖 + semantic fallback tables；site、W、k、read、HP、L1、exact/shape、solver cap、chrX，以及 candidate/hidden/L2/read-AF/L3 狀態表。HCC1395 額外保留 113,997→80,234→25,639 site funnel、7,931→7,928→7,927→7,590 W funnel與 1 region/4 sSNV gap。
- **主張上限**: regional mutation-state candidate sets；不可建立 cell clone/subclone 的存在、比例或祖先關係。candidate-tree read-AF 0/7 artifacts；L3 methylation 7/7 not_evaluated、0 units evaluated。
- **輸入**: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2`；`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_v1`。
- **執行命令**: `python3 InterSubMod/research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py --run-root <historical_run> --production-root <production_probe> --output-dir InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01 --shell-template <data-analytics-html-shell> --embed-helper <runtime-embed-helper> --hcc-hc-bed <SEQC2-HC-BED> --hcc-cn-bed <SEQC2-CN-BED>`。
- **主輸出**: `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html`。
- **驗證**: artifact hard gate PASS；9 chart hosts/9 live SVG/9 no-JS fallback tables；932 source tooltips；duplicate DOM IDs=0；runtime marker=0；1440/768/390、light/dark、keyboard tooltip、console/page errors 全部 PASS。browser QA：`InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.browser_qa.json`。
- **Evidence tier**: L1 artifact/data recomputation + browser QA；scientific release 仍 NO-GO。

### [2026-07-11 04:57] 最終流程說明書敘述框架預註冊

- **Status**: Accepted for report construction；證據狀態仍依run gate動態更新。
- **5W**: Who=使用者/技術 reviewer/future agent；Why=解釋+重驗+決策；What=production tagging→sidecar→validator→frozen run全流程；When=離線長文；How=Markdown SoT + standalone HTML。
- **主框架**: Hybrid-PI-Report L0-L3 + A3/ADR + Diátaxis；L0結論、L1流程、L2判斷/設定、L3命令/證據/限制。
- **備框架**: WRAP+Pre-Mortem用於launch gate；不以純SCQA作主框架，因sidecar equivalence與producer完成度仍在Complex/PROBE域。
- **預期sections**: TL;DR、範圍/名詞、歷史問題、producer命令、sidecar契約、fail-closed、frozen provenance、resource gate、Step→Verify、7-dataset狀態、復原/限制、evidence index。
- **防HARKing**: HTML只呈現source `.md`已有數字；status改變時先更新Markdown與source mapping，再更新HTML。

### [2026-07-11 04:55] Frozen execution 與 atomic lifecycle

- **Status**: Accepted；tiny canonical success 與三個 adversarial cases 已通過。
- **決定**: preflight 在正式 run root 外執行；成功後 freeze manifest 與 source bundle，workers 只讀 frozen copies；失敗建立 `_FAILED`，全部 verification/hash 完成後才原子建立 `_SUCCESS`。
- **驗證**: 7/7 synthetic datasets、35 splits、tree/site-ledger/verifier/source hashes PASS；empty dataset、missing preflight input、existing run root 均 non-zero 且無 `_SUCCESS`，preflight failure 不建立正式 root。
- **附帶修正**: 無 multi-locus group 的 split 現在明確輸出 exposure/match/missing/conflict=0，避免缺欄位 `null` 破壞固定 schema。
- **實作**: `run_layered_7samples_newbb.sh`、`test_layered_runner_fail_closed.py`、`sm_multilocus_combinations.py`。
- **Evidence tier**: L1 synthetic execution；real 7-dataset launch receipt pending。

### [2026-07-11 04:43] Producer binary 與 runtime snapshot

- **Status**: Observed at active-run cutoff；final readback pending。
- **LongPhase-S**: version `1.0.0`；binary SHA-256 `07cbd0aa0c9f33ed59c5baff45fbe3554ef96d55457635de7348c4501b283f54`；size `32,860,600` bytes。
- **Runtime**: Python `3.9.12`、pysam `0.23.3`、samtools `1.17-25-geb3be52`（htslib 1.18）、bcftools `1.13`（htslib 1.18）、bgzip/tabix `1.18`。
- **Run-start source readback**: active production `code.sha256` 於 04:34 仍 4/4 OK；需在7/7完成後再驗一次並保存source bytes，才構成frozen provenance。
- **限制**: binary hash已記錄但active producer原始 `code.sha256`未包含binary；v3 provenance builder必補，不得只引用版本字串。

### [2026-07-11 04:40] coordinate_join_v1 bounded 等價 gate

- **Status**: Bounded PASS；production/full-scope pending。
- **真實資料**: COLO829 `chr1:4,386,684-4,388,348`，35 alignment exposures；direct tagged BAM與raw BAM+sidecar的 calls/HP/PS/family combined digest皆為 `b466e584ab820448769bd1aede31434021f3c6ce1762788afc7363c0fe2d71d7`，exact=35/35、missing=0、conflict=0、duplicate=0。
- **Synthetic edge cases**: same-QNAME primary/supplementary、missing row、conflicting duplicate 3/3 PASS。
- **證據限制**: historical tags只作join方法對照，不是clean production evidence；真實 probe只有primary alignments，supplementary由synthetic fixture覆蓋；contract名稱固定 `coordinate_join_v1`，不得寫成payload identity v2。
- **Gate影響**: U5的方法層由UNKNOWN升為L2 bounded evidence；仍須production 7/7全域unique/hash/BAM binding與35/35 runtime split conservation，才可解鎖full run。
- **Evidence**: `InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/equivalence_probe.json`。

### [2026-07-11 04:24] Sidecar acceptance scope 與 LongPhase filtering policy 澄清

- **Status**: PROVISIONAL / PROBE；此條覆寫本文件 04:08 條目中容易誤讀的「Accepted for clean rerun」措辭，但不停止已在執行的 attach-only producer。
- **已核准範圍**: 既有 FIFO run 可繼續產生 7-dataset alignment-level HP/PS sidecar，作為容量可行性、validator 與 exact-join 的 probe evidence。
- **尚未核准範圍**: sidecar 尚不可宣稱為 production tagged-BAM replacement；尚不可據此啟動新的 layered full run。解鎖仍須 7/7 producer validation、frozen schema-2.1 manifest、COLO829 bounded direct-vs-sidecar equivalence、adversarial fail-closed tests。
- **目前 producer 語意**: 命令無 `--truth-vcf/--truth-bed`，有 `--output-somatic-vcf`，且沒有 `--disableFilter`；因此是 `production_default_filter` tagging + recalibration，不是 `caller_pass_tag_all`。
- **不可混寫**: `_sc.vcf` 保留所有 ClairS PASS record keys，只證明輸出 key conservation；它不證明每個 ClairS PASS variant 都通過 LongPhase-S 內部 filter並參與 tagging。
- **替代 policy**: 若目標被明確定義為「所有 ClairS PASS 都作 tagging seed」，命令必須另開新 role 並加 `--disableFilter`；不得把目前執行中的 run 事後改名成該 role。
- **Evidence tier**: L1 command/source audit；production equivalence L2 pending。

### [2026-07-11 04:12] ClairS 全位點與 LongPhase-S role 分離

- **Status**: Partially superseded by 2026-07-11 04:50；raw/PASS role 分離保留，primary backbone 句改判。
- **決定**: We will preserve every raw ClairS record in a per-site ledger; only ClairS FILTER=PASS records enter LongPhase-S. LongPhase-S `_sc.vcf` FILTER=PASS, not ClairS PASS directly, enters the canonical primary sSNV universe. ClairS PASS input and LongPhase-S recalibrated-all/pass outputs use distinct manifest roles and filenames。
- **理由**: 7/7 record-multiset audit證明 PASS input與 `_sc.vcf` all keys相等，但 raw ClairS另有 non-PASS records；role相等與key相等不是同一命題。
- **驗證**: `clairs_longphase_ssnv_contract_audit.tsv`；HCC1954 raw-all ledger 26,823/26,823 rows、5 checks PASS。
- **Evidence tier**: L1 artifact/code audit。

### [2026-07-11 04:08] Production HP/PS sidecar exact alignment contract

- **Status**: Accepted for clean rerun（7-dataset production generation in progress）。
- **決定**: We will stream the LongPhase-S tagged BAM through FIFO and persist only `chrom/start/end/QNAME/FLAG/MAPQ/CIGAR-digest/HP/PS`; layered reconstruction consumes the original tumor BAM plus this sidecar。
- **理由**: 避免 1.94 TB duplicate BAM payload，同時保留 primary/supplementary alignment identity；不得用QNAME最後一筆覆寫。
- **Smoke evidence**: HCC1954 chr22 execution/capture alignments=217,023；tagged=124,449；unknown HP=0；exact conflicts=0；sSNV exposures/matches=4,432/4,432，missing/conflict=0/0。
- **Revisit if**: clean full run任一 split exact join不守恆，立即 fail-closed並回退到實體 tagged BAM storage plan。
- **Evidence tier**: L1 smoke + L2 full validation pending。

### [2026-07-11 04:05] 大型 workspace 查詢紀律

- **Status**: Accepted。
<!-- BEGIN USER-SPECIFIED -->
**Decision**: We will navigate large data trees through known indexes, manifests, `run_context.json`, `upstream_dependency.tsv`, audit inventories, and exact referenced paths; we will not use broad recursive searches across `/big7_disk` or `/big8_disk` as the default discovery method.
**DO NOT change**: 使用者明示要求避免因資料夾與內容過多造成大量搜尋。
**Rationale**: 優先用關聯導航可降低 I/O、不必要 metadata scan、與把 archive 誤當 active source 的風險。
<!-- END USER-SPECIFIED -->
- **查詢階層**: L1 `CURRENT_FOCUS/00_INDEX`→L2 manifest/run root metadata→L3 exact producer/consumer files→L4 單點 `rg file:pattern` 或限深度列舉。
- **Revisit if**: 索引或manifest缺失；當時先建 bounded inventory，仍不直接搜整顆磁碟。
- **Evidence tier**: L1 user-specified governance。

### [2026-07-10 22:44] Primary lineage 定義

- **Status**: Accepted。
- **背景**: 舊 driver 將 HP1/HP2/HP3 的 root-only family 都算 lineage。
- **決定**: We will define primary lineage as mutation-bearing HP1/HP2 family only。
- **理由**: root-only 是 local reference control；HP3 是 unresolved somatic-integrated haplotype，兩者均不能支持 primary parental lineage claim。
- **影響範圍**: layered schema、region determinacy、multi-HP、CCF/read-AF、所有 denominator。
- **Revisit if**: 有獨立 parental phasing 或 single-cell truth 可解析 HP3。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

<!-- BEGIN USER-SPECIFIED -->
**Decision**: chr1-22 為 primary scope；HP3 標 `H3?` 且退出 primary denominator；論文主張限 regional mutation-state trees。
**DO NOT change**: 2026-07-10 使用者以「全部完成」核准前次建議預設。
**Rationale**: pre-decision audit + repository KB semantics。
<!-- END USER-SPECIFIED -->

### [2026-07-10 22:44] Missing CN 語意

- **Status**: Accepted。
- **決定**: We will emit `unavailable` when no CN source is loaded; `neutral` is allowed only when a declared segmentation source covers the interpretation contract。
- **理由**: missing measurement 不等於 biological diploidy。
- **影響範圍**: mlhp `cn`、L2 verdict、recurrence census、CN-clean denominator。
- **Revisit if**: COLO829/HCC1937 取得可用 orthogonal CN truth。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

### [2026-07-10 22:44] Verification truthfulness

- **Status**: Accepted。
- **決定**: We will separate `full_pass`, `partial_pass`, and `fail`; only zero-skipped V1-V7 may produce `all_V1V7_pass=true`。
- **理由**: executed-check pass 與 complete verification 是不同命題。
- **影響範圍**: unit schema、aggregate、runner與 verifier wording。
- **Revisit if**: V4/V5 independent oracle 被等價且更快的正式證明取代。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

### [2026-07-10 22:44] Candidate completeness 與 display storage 分離

- **Status**: Accepted。
- **決定**: We will compute shape and read-AF weighting over all candidate trees; any display cap must not affect analysis and must be explicit in schema。
- **理由**: UI 儲存限制不能改變方法學結果。
- **影響範圍**: layered output、shape determinacy、read-AF weighting、workstation。
- **Revisit if**: full candidate output exceeds twice the current storage/runtime budget。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

## 🟠 偏離之處

### [2026-07-11 04:02] 實體 tagged BAM 改為 streaming HP/PS sidecar probe

- **Status**: Consumer contract conditionally accepted；7-dataset production completion 仍 pending。
- **規範要求**: We will regenerate 7 production tagged BAMs without `--truth-vcf/--truth-bed`, then feed them to the fail-closed layered runner。
- **實作偏離**: 外部 AI 已啟動 LongPhase-S FIFO streaming，保留每個 mapped alignment 的 HP/PS sidecar 與 recalibrated VCF，不保留 BAM payload。
- **理由**: 7 個 tumor BAM 合計 1,939,360,284,288 bytes，但 `/big7_disk` 只剩 932 GB；完整副本無法放入正式 output root。
- **風險評估**: layered code 已支援原 BAM＋sidecar exact join，且 supplementary/duplicate identity、BAM drift 與 tabix boundary 都有 fail-closed gate；剩餘風險是 7-dataset producer 尚未全部完成。
- **回退路徑**: exact-equivalence 失敗時，不接 sidecar；另申請≥2 TB 正式儲存空間，或設計可審核的 per-chromosome materialize→consume→archive 流程。
- **Revisit if**: 7/7 producer validation、manifest lock 或 downstream exact join 任一失敗。
- **Evidence tier**: COLO829 real-BAM 35/35 exact-equivalence + HCC1954 chr22 4,017/4,017 exact joins；7/7 production pending。

### [2026-07-10 23:25] Second-caller backbone 只能做單一 biological sample

- **Status**: Accepted。
- **規範要求**: Comprehensive validation 原則上不應使用 subset。
- **實作偏離**: DeepSomatic sensitivity 只在 HCC1395_DORADO；LPS PASS／SEQC2-HC 只在 HCC1395，顯著標 partial scope。
- **理由**: workspace 只有 HCC1395 ONT_Dorado DeepSomatic callset，其他 5 biological samples 無第二 caller VCF；不可捏造不存在的 cross-caller evidence。
- **風險評估**: 只能診斷 backbone dependence，不能宣告跨樣本 caller robustness。
- **回退路徑**: 其他樣本 DeepSomatic／Mutect2 callset 到位後，以同 manifest runner 擴充。
- **Revisit if**: ≥3 biological samples 有 orthogonal caller PASS VCF。
- **Evidence tier**: L1 data-availability audit。

### [2026-07-10 23:00] Capped units 不納入 full V1-V7 分母

- **Status**: Accepted。
- **規範要求**: 初版 pre-registration 寫所有 primary units 的 V4/V5 都不得 skipped。
- **實作偏離**: full verification denominator 改為所有 non-capped eligible units；capped units 標 `not_applicable_capped`，不得進 exact-tree／shape claim。
- **理由**: capped 的定義就是枚舉未完成，V4/V5 completeness 無可驗命題；硬算 pass 或 fail 都會失真。
- **風險評估**: 全體 unit 的單一「ALL PASS」口號將被淘汰，改報 eligible coverage、not-applicable 與 fail 三個數。
- **回退路徑**: 若未來 solver 能完整處理 dense units，capped 歸零後自然進 full denominator。
- **Revisit if**: capped solver 被 exact algorithm 取代。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐。

### [2026-07-10 22:44] 無法完成 single-cell confirmation

- **Status**: Accepted。
- **規範要求**: ⭐3→⭐4 需 orthogonal confirmation。
- **實作偏離**: 現有 repo 只能完成 bulk-read regional mutation-state reconstruction；不宣告 biological clone confirmation。
- **理由**: repo 無 single-cell／multi-region truth input，不能用計算補造證據。
- **風險評估**: 限制 paper claim tier，但不阻擋工程與可重現性修正。
- **回退路徑**: 外部資料到位後依 manifest 新增 orthogonal validation cycle。
- **Revisit if**: 新資料可用。
- **Evidence tier**: L1 absence audit（資料盤點待最終重驗）。

## 🟡 折衷考量

### [2026-07-11 04:50] Canonical tree backbone 改為 LongPhase-S recalibrated PASS

- **Status**: Accepted；使用者明確改判，取代 2026-07-10 22:44 的 canonical-backbone 決策。
- **新契約**: raw ClairS all 作 lossless ledger；ClairS PASS 作 LongPhase-S input；LongPhase-S `_sc.vcf` all 保存 recalibrated FILTER；`_sc.vcf` FILTER=PASS 作 primary sSNV tree input。
- **理由**: 若研究主張是 LongPhase-S 整合後的演化樹，primary variant backbone 必須消費 LongPhase-S output；只用 ClairS input 會讓 LongPhase-S 僅提供 read tags，方法角色不完整。
- **驗證**: `raw PASS=LPS input`、`LPS input=_sc all`、`_sc PASS=tree input` 三個 record-key multiset gates；每個 raw record 另有唯一 disposition。
- **限制**: LongPhase-S filtering 與 tree read evidence 使用同批 reads，屬 operational pipeline，不是獨立 biological validation。
- **回退/比較**: 原 ClairS PASS tree backbone 降為 7-dataset explicit sensitivity，不作 canonical claim；它檢驗 LongPhase-S recalibration 影響，但不等同 independent-caller robustness。
- **Evidence tier**: L1 method/data-flow contract；7-dataset artifact validation pending。
- **Bounded real-data verification**: HCC1954 chr22 raw/LPS-input/`_sc`-all/tree = 426/181/181/167；LongPhase-filter exclusion=14；tree funnel 167=singleton 99+retained 68；sidecar exact 4,017/4,017、missing/conflict=0/0。

### [2026-07-11 05:13] 全 7-dataset backbone sensitivity fail-closed gate

- **Status**: Accepted；取代把 ClairS tree sensitivity 限為 HCC1395 partial 的舊做法。
- **決定**: canonical 一律使用 LongPhase-S `_sc.vcf` FILTER=PASS；同一批 7 datasets 另以 ClairS PASS input 重跑，量化 retained-position Jaccard、determined-primary、multi-HP 與 region-determined delta。
- **資料範圍**: 7 datasets / 6 biological samples；不再用 3 個單樣本 historical probes 滿足正式報告 gate。
- **Fail-closed 條件**: base/sensitivity 都要 `_SUCCESS`、7/7 verifier PASS、frozen manifest sample order/set 相同、tree contract 分別為 `longphase_recalibrated_PASS` / `clairs_PASS_input`。
- **測試**: `run_layered_v2_contract_tests.sh` 已擴為 26/26 PASS；backbone 3 tests 驗證 7/7 成功、缺 dataset 失敗、非 canonical summary 被 report 拒絕。
- **限制**: 此比較是同一 ClairS→LongPhase-S pipeline 的上游/下游 backbone sensitivity，不可改寫為 independent-caller robustness；舊 DeepSomatic/SEQC2-HC 仍是 historical partial artifact。
- **Evidence tier**: L1 code/test contract；真實 7-dataset sensitivity results pending production completion。

### [2026-07-11 05:28] Producer closeout 與下游 handoff 原子 gate

- **Status**: Accepted；補 active v1 producer 原始 runner 未內建 `_SUCCESS`/binary lock 的 provenance 缺口。
- **執行中 binary 證據**: 05:15 對 PID 656465/656470 的 `/proc/<pid>/exe` 與磁碟 LongPhase-S binary 重算 SHA-256，三者皆為 `07cbd0aa0c9f33ed59c5baff45fbe3554ef96d55457635de7348c4501b283f54`，device/inode 皆 `59/67474570`；receipt 已寫入 active production root。
- **Closeout 契約**: 7/7 status/verification、啟動 code hashes、每樣本 input/output hashes、inventory size/mtime、sidecar tabix、VCF CSI、ClairS PASS=`_sc` all、`_sc` PASS=tree VCF、HP vocabulary/PS、無 truth flags 全部重驗；freeze source+binary；最後才原子建立 `_SUCCESS`。
- **Handoff 契約**: `prepare_clean_layered_manifest.py` 必須重新核對 `_SUCCESS`→closeout→artifacts manifest hashes，並逐一 readback 下游會用的 sidecar/VCF/index/validation/input inventory；任何 post-closeout drift 都拒絕產生 clean manifest。
- **真實路徑修正**: HCC1395_DORADO BAM 是 symlink；producer `stat -c` 鎖 logical symlink metadata，因此 closeout inventory 必用 `lstat()`，不能用 target `stat()`。已加入 symlink regression fixture。
- **測試**: producer/handoff 6 tests PASS（success last、input mutation、binary mismatch、closeout accept、unclosed reject、post-closeout drift reject）；全 suite 26/26 PASS。
- **限制**: active producer 本身不是 v3 tracked-child lifecycle；closeout 能證完成結果與 provenance，但不能回溯改變它啟動時的 signal/process-control 實作。7/7 未完成前不得發布 `_SUCCESS`。
- **Evidence tier**: L1 synthetic adversarial tests + active `/proc` executable receipt；7/7 production artifact pending。

### [2026-07-11 05:59] Canonical full run 升級採 v3 lifecycle

- **Status**: Accepted；canonical scientific run 優先使用 `InterSubMod/scripts/run_layered_v3.py`，v2 runner保留 ClairS-backbone/parameter sensitivity 相容路徑。
- **理由**: v3 補齊 parent flock、staging publish、source/environment/input lock、tracked child fail-fast、signal reap、heartbeat、NFS/load resource gate、output-manifest binding、runner-only `_SUCCESS`；比 v2 post-hoc closeout 更能保證 U6/U7。
- **元件驗證**: lifecycle 13/13、launcher 12/12、producer receipt 8/8、input/manifest contract 17/17、verifier 16/16、handoff supervisor 20/20、producer closeout 6/6，合計 92/92 PASS；另有 v2/跨版本 contract 34/34 PASS。
- **Handoff 驗證**: supervisor 預設 dry-run，須另有 reviewed JSON authorization 與 out-of-band plan SHA-256 才能執行；active/partial producer、6/7 verification、artifact/source drift、既存 run ID 或外部 v2 launcher 都會 fail closed。
- **受控啟動**: 使用者已要求「全部完成」；reviewed plan 為 `InterSubMod/research/20260710_layered_reconstruction_v2/launch_plans/20260711_layered_v3_handoff_launch_plan_reviewed_v1.json`，bytes SHA-256=`7e0ab871ee2fa15e772d17da75ebe0e836cc7dcc6c2726e7b136a3b305e3da6a`。supervisor session 62875 目前只持 lock 並等待 active producer，尚未建立 layered run root；producer 7/7、receipt/closeout/resource gate 任一未通過即停止。
- **混合執行契約**: canonical base 可為 v3 frozen lock；ClairS-backbone、MAPQ30、BaseQ10、MAX_SNV6 與 derived MINREAD4 可為明示 v2 sensitivity。Summarizers 已要求兩種 root 都有 `_SUCCESS` 且各自 frozen contract 正確，並直接比較 site/unit/topology digest。
- **Resource gate**: 正式 v3 launch不得放寬 frozen production policy（load/core≤1.25、iowait≤25%、NFS read<80 MB/s、固定 300s baseline）；若未通過，不建立 formal run root。
- **PS 契約**: HP 作 family stratification；PS 逐 alignment 保存並輸出 per-region PS/HP×PS census，只作 phase-block QC。HCC1954 chr22 已觀察 1/31 mixed-PS retained region。
- **Evidence tier**: L1 adversarial component tests；real 7-dataset v3 preflight/full run pending producer completion。

### [2026-07-11 07:00] LongPhase-S input 改為 normalized raw ClairS all

- **Status**: `GO_WITH_FAIL_CLOSED` for 7-dataset validation；取代「ClairS PASS 作唯一 LongPhase-S input」的 production contract，但不改變 canonical tree 使用 `_sc.vcf FILTER=PASS` 的決定。這是 launch authorization，不是 production complete。
- **新證據**: LongPhase-S source `HaplotagVcfParser.cpp:589-598` 原生支援 non-PASS→PASS rescue；7 datasets raw ClairS 共 70,179 non-PASS records。PASS-only input 使 rescue branch 永遠不可達。
- **停止錯誤長跑**: `20260711_longphase_s_production_sidecars_v1` 於 0/7 terminal PASS 中止，保留為 `20260711_longphase_s_production_sidecars_PASS_ONLY_ABORTED_v1`；`_FAILED` 存在、`_SUCCESS` 不存在、無殘留 process。舊 reviewed handoff plan已 superseded，未建立 layered run root。
- **HCC1954 chr22 raw-all**: 426→426 unique biallelic sSNV；LowQual→PASS 53、PASS→LowQual 14、PASS→PASS 167、LowQual→LowQual 192；217,023 sidecar rows，unknown/conflict=0。
- **已知 crash**: 原 binary在 HCC1395 chrX:72880028（QUAL=0 LowQual、無 q20 eligible read）exit 1。copied source只把 low-confidence/no-read `calibrateReadHP` fatal改為 warning+no-op；該位點 patched run 1→1、LowQual不變、64 reads全 HP=.、exit 0。
- **Regression**: patched vs original binary在 HCC1954 chr22 的 sidecar bytes SHA-256皆 `49eca6ff331981ee55710a15a5f3135e1dceb35086eadb1fe845eb4e4c1c041c`，426 VCF record lines digest皆 `cdad75060d2396da7a1bf9b784d95f05c9b12beac135dbf87d95105a040ea522`，purity皆0.94308，transition matrix完全相同。
- **HCC1395 whole chrX patched**: 35,823→35,823；LowQual→PASS 959、PASS→LowQual 1,616；2 no-read warnings；1,090,570 sidecar rows、unknown/conflict=0；artifact hashes全PASS。
- **Current gate**: 7-dataset whole-genome raw-all producer/receipt/closeout、U0–U7與 fresh downstream run尚未啟動完成；任一新 fatal/key loss/payload drift/unknown/conflict/hash failure均不得產生 `_SUCCESS`。
- **Evidence tier**: L1 source/runtime/record/read-tag probes；custom binary尚未獲 full-dataset validation或 upstream review。

### [2026-07-11 07:58] raw-all receipt v2、layered v3 gate 與資源排隊

- **Status**: Accepted / execution queued；不是 production complete。
- **Receipt v2**: 明確凍結 `caller_raw_vcf`、normalized raw-all LongPhase input、ClairS PASS sensitivity baseline、LongPhase `_sc` all、LongPhase `_sc` PASS 五個 VCF 角色；另綁 normalization/FILTER transition/sample verification、patch approval、全域 sidecar identity與 FIFO-only BAM policy。
- **BAM policy**: producer 只以 named FIFO 消費 LongPhase-S tagged BAM stream；每樣本 closeout 強制 `regular *.bam = 0`、FIFO 存在。既有 canonical tagged BAM不在新 `-o` 路徑下，未被寫入或覆蓋。
- **Layered v3**: canonical tree role固定為 normalized raw-all run 的 LongPhase-S recalibrated `FILTER=PASS`；ClairS PASS獨立保留為 sensitivity-only baseline。site ledger 需逐筆保存 raw→recalibrated FILTER transition，並驗證 raw=input=recalibrated-all record count。
- **Canonical tests**: `InterSubMod/scripts/run_layered_v3_raw_all_contract_tests.sh` 實跑 **84/84 PASS**（producer 5、receipt 4、read-tag 3、site-ledger 2、canonical-BAM immutability 3、backbone sensitivity 6、lifecycle 13、launcher 13、manifest 18、independent verifier 17）；另 v2/跨版本 **45/45 PASS**。
- **Canonical BAM baseline**: 7 個既有 paired-full LongPhase-S tagged BAM 已凍結 logical/real path、inode、size、mtime/ctime與首/中/尾各 1 MiB SHA-256；identity set=`ce6c63d42e3f334d6847a1a6d3e46ead165b59a03197acb098319be5c67fcf90`。收據：`InterSubMod/research/20260710_layered_reconstruction_v2/audit_receipts/20260711_canonical_longphase_tagged_bam_baseline_v1.json`；正式 run 後必須產生 exact-match verification receipt。
- **Resource gate**: 07:57 實測 load1=64.80、iowait=72%，不適合啟動 NFS-heavy LongPhase-S；已排隊，僅當 load1<60 且 iowait<25% 才以 1 sample×12 threads 啟動，輸出 root 尚不存在。
- **Production target**: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2`。
- **Evidence tier**: L1 code/adversarial tests + live resource observation；7/7 runtime evidence pending。

### [2026-07-11 08:32] canonical 與 sensitivity 的 v3 雙契約

- **Status**: Accepted；runtime 尚待 production input。
- **Canonical**: selected tree=`LongPhase-S recalibrated FILTER=PASS`；task=`comprehensive_validation`。
- **Sensitivity**: selected tree=`ClairS FILTER=PASS`；task=`backbone_sensitivity`；只替換 tree backbone，不重跑/替換 LongPhase-S read tags。
- **不可交換角色**: frozen lock 同時保存 caller raw、normalized raw-all、caller PASS、LongPhase-S all、LongPhase-S PASS、selected tree 六個角色。sidecar subject 永遠綁 `longphase_recalibrated_pass_vcf_sha256`，不綁可替換的 selected tree。
- **Verifier**: 依明示 tree contract重算 selected-role equality、site-ledger tree contract與 output somatic roles；未明示 sensitivity 時 ClairS PASS 替換仍 fail-closed。
- **Comparison**: `summarize_backbone_sensitivity.py` 現可獨立驗證兩個 v3 `_SUCCESS` run，輸出 retained-position Jaccard、primary-unit Jaccard、shared topology digest concordance與三組比例差。
- **Evidence tier**: L1 synthetic/adversarial contract tests；7-dataset runtime sensitivity pending。

### [2026-07-11 08:45] producer 與 downstream 等待期 TOCTOU pinning

- **Status**: Accepted / queued；尚未建立 production root。
- **Producer gate**: `InterSubMod/scripts/launch_longphase_raw_all_when_idle.sh` pins 20 個 authority（Python、producer/helpers/schema、manifest、patched binary、patch evidence、jq/samtools/bcftools/bgzip/tabix），每輪 resource check與真正 launch 前都跑 `sha256sum -c`。
- **Downstream supervisor**: `InterSubMod/scripts/complete_raw_all_layered_v3_validation.sh` pins 18 個 source與 5 個 input authority；producer terminal後及 canonical manifest、sensitivity manifest、兩個 full run、comparison、BAM post-check每一步前重驗。
- **Superseded evidence**: 第一版 downstream supervisor只記錄hash、等待後不重驗，於科學執行前主動 SIGINT fail-closed；保留於 `InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260711_raw_all_to_layered_v3_completion_SUPERSEDED_UNFROZEN_v1/`，無 production/downstream root。
- **Active logs**: `InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260711_raw_all_producer_resource_gate_v1/execution.log` 與 `InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260711_raw_all_to_layered_v3_completion_v1/execution.log`。
- **Evidence tier**: L1 hash/lifecycle mechanism + L2 live resource observation。

### [2026-07-11 09:00] raw-all full producer正式啟動

- **Status**: Running；0/7 terminal PASS，不是production complete。
- **Launch gate**: `09:00:27+08:00` load1=1.98、iowait=0%；20/20 launch authorities逐一 `OK`，`09:00:28` LAUNCH。
- **HCC1395 input**: normalization 134,122→134,122、rescue=0、remove=0；LongPhase runtime解析Tumor SNP count=134,122、truth VCF/BED皆空、tag region=all、supplementary enabled、native filtering enabled。
- **BAM policy runtime**: `HCC1395_production.bam` 實測 type=`fifo`、size=0、inode=614514135；capture process讀同一FIFO，canonical BAM不是`-o`目標。
- **Phase timing**: matched-normal BAM extraction=1280s；tumor BAM extraction進行中。sidecar writer尚未開始，沒有sample `_SUCCESS`。
- **Production root**: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2`。
- **Evidence tier**: L2 live runtime；terminal filter/read-tag/artifact gates pending。

### [2026-07-11 15:00] 教授重點 PPT-style HTML 簡報

- **Status**: internal PI-share artifact PASS；scientific/publication release NO-GO。
- **Task type / goals**: D handoff；服務 G3/G4/G5。
- **輸入**: `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.data.json`（SHA-256=`78cec57adff20d6d3c21acd29e75499acee3be220fdc27012a0029bcbe9d3042`）；producer progress 另讀 raw-all `run_status.tsv`。
- **輸出**: `/big7_disk/liaoyoyo2001/Meeting/interSubMod_reports_workspace/20260711_分層重建教授重點簡報_01/slides/20260711_分層重建教授重點簡報_01.standalone.html`；10 張 16:9、5 個 packaged Recharts 圖表、same-data no-JS fallbacks、名詞/講稿/來源/總覽抽屜與列印模式。
- **重生命令**: `python /big7_disk/liaoyoyo2001/Meeting/interSubMod_reports_workspace/20260711_分層重建教授重點簡報_01/scripts/build_slide_deck.py` 後接 `qa_slide_deck.py`。
- **強驗證**: official chart contract 0 errors/0 warnings；1920×1080、1366×768、390×844、320×568 皆 5/5 live SVG；no-JS desktop/mobile 皆 10/10 slides + 5/5 fallbacks；print 10/10 pages；console/page errors=0。
- **數字邊界**: 投影片保留 historical 7-row aggregate `U=568,080`、`W_tree=48,959`、`W_primary=47,377`、complete=`39,885`、topology three-way=`10,832/11,144/17,909`、hidden=`7,325/32,560`；partial raw-all producer 進度只改 status card，不替換 tree rates。
- **Claim ceiling**: HP ≠ clone；H3? auxiliary；C/Topo 是模型候選與樹形；hidden state ≠ hidden clone；read-AF/methylation 皆 0 evaluated units。
- **獨立 reviewer cycle**: 第一輪找到 print PDF 第 11 張空白頁與 chrX notes 不存在路徑；修正後 `pdfinfo Pages=10`、`pdftotext` 末頁含 `Scientific NO-GO`、證據路徑存在，數字/generator/scientific 三視角最終 **3/3 PASS**。
- **Evidence tier**: artifact/browser/reproducibility L2；biological confirmation L5 pending。

### [2026-07-11 04:02] 儲存量與 production equivalence

- **Status**: Proposed。
- **方案 A**: We will persist 7 full tagged BAMs；拒絕作為當前預設，因 big7 空間小於 1.94 TB input payload，且 bip7 只剩 2.1 TB，無安全餘裕。
- **方案 B**: We will retain original BAM bases and externalize only clean LongPhase-S HP/PS tags；採為 PROBE，因可減少重複 payload，且理論上可以 alignment identity exact join。
- **方案 C**: 逐 chromosome 實體化 BAM 後刪除；拒絕，直接刪除與 repo governance 衝突，且增加 partial-state recovery 難度。
- **採用 B 前提**: exact join、content identity、schema、hash、fault injection 全部過關；未過前不啟動 layered full run。
- **Revisit if**: 有新正式儲存空間或 compressed CRAM 經實測能在預留餘裕下存放。
- **Evidence tier**: L1 disk/input sizes；L3 consumer equivalence pending。

### [2026-07-10 22:44] Canonical backbone 與 sensitivity（已被 2026-07-11 04:50 取代）

- **Status**: Superseded；保留作決策歷史。
- **方案 A**: We will retain ClairS paired PASS as operational primary backbone and evaluate alternatives as sensitivity only。
- **方案 B**: 以 LongPhase-S recalibrated PASS 取代 primary；拒絕，因會混入 downstream recalibration 定義。
- **方案 C**: 以 SEQC2-HC 當全域 universe；拒絕，它是可評估 truth subset，不是完整 somatic caller universe。
- **採用 A 理由**: caller provenance 清楚且可跨 7 datasets；不把 operational universe 說成 biological truth。
- **Tier check**: L2，敏感度尚待本 cycle 產物。
- **Revisit if**: second caller 在多 biological samples 顯示穩定結構差異。
- **Evidence tier**: L2 ⭐⭐⭐⭐。

### [2026-07-10 22:44] L3 methylation contract

- **Status**: Accepted。
- **方案 A**: We will expose an explicit `not_evaluated/bounded_auxiliary` contract and prohibit tree ranking／confirmation。
- **方案 B**: 保留 `pending`；拒絕，pending 容易被誤讀為未來會升級主證據。
- **方案 C**: 從 pipeline 移除 L3；延後，會失去 bounded residual flag 的研究接口。
- **採用 A 理由**: 與既有 NEGATIVE filter evidence 相容，且保留非循環輔助角色。
- **Revisit if**: 預註冊的 orthogonal methylation validation 通過。
- **Evidence tier**: L2 ⭐⭐⭐⭐。

## 🔴 未決問題

### [2026-07-10 22:44] 外部 truth availability

- **Status**: open。
- **Question**: COLO829/HCC1937 CN truth 與 single-cell／multi-region資料何時可取得？
- **Context**: 這是唯一不能由目前 workspace 自行補齊的證據層。
- **Options**: 外部提供；公開資料映射；維持 unavailable／claim ceiling。
- **Default if no answer**: 維持 unavailable，明示不升級 biological clone claim。
- **Revisit if**: 新資料或 collaborator handoff。
- **Priority**: major。
- **Evidence tier**: L5。

## 📚 Lore

### [2026-07-10] 文件修正不等於 core 修正

- **Constraint**: CURRENT_FOCUS 已寫 P0 已修，但 upstream scripts 仍保留 neutral fallback、legacy lineage 與 sampled verification。
- **Why it matters**: status 必由 fresh code＋artifact verification 支持，不能只由報告文字繼承。
- **Evidence**: `InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md` §1。

### [2026-07-10] Region gap 不是 total span

- **Constraint**: region 是相鄰 sSNV gap <=50 kb 的 connected component；總 span 可超過 50 kb。
- **Why it matters**: 方法與 UI 不得寫「region <=50 kb」。
- **Evidence**: `sm_linkage_genomewide.py` grouping implementation + prior 7-sample audit。

## Provenance

- Commit at init: `4fb9e742482b63a660de19a1f1bd07d49d713111`
- Cycle: `InterSubMod/state/cycles/20260710-2244-layered-reconstruction-v2/`
- Pre-decision audit: `InterSubMod/research/20260710_layered_reconstruction_v2/pre-decision-audit.md`
- Status: in_progress。
