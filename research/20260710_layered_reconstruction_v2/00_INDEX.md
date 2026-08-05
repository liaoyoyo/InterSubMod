<!--
建立時間: 2026-07-10 22:44
目標: 註冊 layered reconstruction v2 全量修正與驗證門檻
處理範圍: chr1-22 x 7 datasets（6 biological samples）；程式、重跑、驗證與文件同步
build_branch: master
build_commit: 4fb9e742482b63a660de19a1f1bd07d49d713111
worktree: /big7_disk/liaoyoyo2001/InterSubMod
data_sources: InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts,InterSubMod/docs/CURRENT_FOCUS.md
驗證方式: 本檔 decision thresholds 對照 canonical v5 與 sensitivity v6 verification_summary.json、current_layered_topology_v3_raw_all_v1.json
證據等級: L2 ⭐⭐⭐⭐（7 datasets 全量工程/模型內驗證；無 single-cell/multi-region orthogonal truth）
-->

# Layered Reconstruction v2

> **TL;DR（2026-07-14）**：canonical LongPhase-S recalibrated PASS v5與ClairS PASS sensitivity v6皆7/7 PASS；最新sSNV共現/C/Topo彙總已完成，主結果固定使用LongPhase-S PASS，clone/ancestry claim仍受orthogonal truth限制。

> **2026-07-11 P0 scope correction**：歷史 tagged BAM 中 6/7 受 LongPhase-S `--truth-bed` 限制；`20260710_232501_layered_reconstruction_v2` 雖內部 7/7 PASS，只能作 engineering baseline。focused audit：`InterSubMod/research/20260710_layered_reconstruction_v2/20260711_ClairS_LongPhaseS_sSNV位點與ReadTag契約稽核_01.md`。

> **2026-07-12 live status（18:16 +08:00）**：raw-all producer/receipt closeout已7/7 PASS（638,259 input records；582,820 `_sc PASS`；32,184 rescue；17,444 remove；164,253,537 sidecar rows；unknown/conflict=0）。fresh canonical layered-v3為4/7完成（COLO829/H1437/H2009/HCC1395）、2 active（HCC1395_DORADO/HCC1937）、HCC1954未啟動，且無root `_SUCCESS`；sensitivity/comparison/final verifier/post-run readback尚未開始。producer gate已關，但完整Results仍NO-GO。

> **2026-07-14 terminal closeout（取代上列live snapshot）**：canonical `20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`與sensitivity `...clairs_pass_sensitivity_v6`均有root `_SUCCESS`且7/7 PASS；BAM post-readback 7/7 match。最新machine summary=`current_layered_topology_v3_raw_all_v1.json`，canonical W=51,815、W_primary=50,215、complete/incomplete=42,240/7,975、C/Topo三類=11,582/10,737/19,921、impossible=0；comparison=`backbone_sensitive`。

## Current Canonical Outputs（2026-07-14）

- Canonical root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`
- Sensitivity root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6`
- 最新完整分析：`InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md`
- Machine JSON/TSV：`InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.{json,tsv}`
- Backbone comparison：`InterSubMod/research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.json`
- BAM receipt：`InterSubMod/research/20260710_layered_reconstruction_v2/audit_receipts/20260713_canonical_longphase_tagged_bam_post_full_validation_v6.json`

## Pre-registration (Confirmatory)

| 預測 H | 否證條件 | Decision threshold |
|---|---|---|
| H1：root-only family unit 是 reference control，不是 lineage | 任一 `reference_only=true` 單位仍有 `is_primary_lineage=true` | 7/7 datasets 必須為 0 |
| H2：HP3 是 `H3?` unresolved auxiliary | 任一 family=3 單位進入 primary lineage／region determinacy 分母 | 7/7 datasets 必須為 0 |
| H3：缺 CN 資料必為 unavailable | 無 CN dataset 出現 `candidate_keep` 或 `neutral` 的量測式裁決 | COLO829、HCC1937 必須為 0 |
| H4：full verification 必須真的執行 V1-V7 | `all_eligible_V1V7_pass=true` 但任一 non-capped eligible unit 的 V4/V5 skipped | 7/7 datasets 必須為 0 eligible-skipped、0 fail；capped 另列 not-applicable |
| H5：六層 funnel 守恆 | `autosomal != singleton + cap_excluded + read_unsupported + retained` | 7/7 datasets exact equality |
| H6：shape 與 read-AF weighting 使用完整候選樹集 | exact shape 只由 stored prefix 計算，或 weighting 候選數小於 `n_trees` | 7/7 datasets 必須為 0 truncated analyses |
| H7：主要 topology 對 backbone 選擇不應劇烈翻轉 | determined-primary、multi-HP、region-determined 任一絕對差 ≥10pp，或 retained-site Jaccard <0.80，或 primary-unit key Jaccard <0.70，或 shared-unit candidate-tree digest concordance <0.80 | robust：<5pp + site J≥0.95 + unit J≥0.90 + topology concordance≥0.95；moderate：<10pp + J/concordance達下限；否則 backbone-sensitive |
| H8：主要 topology 對合理 read/cap 參數應穩定 | MAPQ30、BaseQ10、MINREAD4、MAX_SNV6 任一 dataset 主比例差 ≥10pp、site J<0.80、unit J<0.70 或 shared topology concordance<0.80 | H7 同一 compound gate；不能只用 aggregate rate 宣告 robust |
| H9：read-AF ordering 只可作探索性排序 | 文件或 schema 稱 purity/CN-corrected CCF，或只用 stored prefix | 禁用 CCF claim；全部 non-capped candidates；報完整參數 grid |

## Post-audit P0 Remediation Gates（2026-07-11；不得倒填成 pre-registration）

| Gate | Fail condition | Acceptance |
|---|---|---|
| U0 production tagging scope | 任一 LongPhase-S command含 truth VCF/BED或執行 benchmark BED removal | 7/7 truth flags absent |
| U1 caller record continuity | normalized ClairS raw-all biallelic sSNV multiset≠LPS input，或LPS input≠`_sc.vcf` all multiset | missing=0、extra=0（7/7） |
| U1b tree backbone continuity | `_sc.vcf` PASS multiset≠primary tree input | missing=0、extra=0（7/7）；ClairS PASS 不直接作 canonical tree input |
| U2 exact alignment tags | FIFO capture count≠execution count；未知HP；identity conflict | 全部 exact；unknown/conflict=0 |
| U3 sSNV sidecar join | split中 sidecar exact match≠alignment exposure | missing=0、conflict=0、match=exposure；bounded probe 35/35 alignment-group exposures，正式 7×chr splits pending |
| U3b PS 使用邊界 | PS 未逐 alignment 保存、mixed-PS region 未計數，或把 PS 當 topology edge/lineage label | per-region PS 與 HP×PS census守恆；PS 僅 phase-block QC；mixed-PS 逐樣本報告 |
| U4 all-site ledger | raw record未寫出、ClairS/LPS/output key不守恆，或`_sc` PASS sSNV無唯一 funnel disposition | raw rows守恆；U0 PASS=U1、U1=U2 all、U2 PASS=tree input（7/7） |
| U5 sidecar equivalence | bounded direct tagged-BAM與raw-BAM+sidecar的 eligible observations、R/A/X、family 或 digest任一不等 | synthetic adversarial 3/3 + COLO829 bounded real-data digest全部 exact；否則 sidecar只保留probe role |
| U6 frozen provenance | worker/verifier重讀mutable manifest/source；sidecar未綁BAM/producer；任一鎖定hash漂移 | exact 7/6 strict lock、source/env bundle、artifact identity及launch receipt readback全部PASS |
| U7 atomic lifecycle | preflight失敗仍形成正式run root、empty set可PASS、失敗/中斷建立`_SUCCESS` | adversarial matrix全部non-zero且無正式root/worker/`_SUCCESS`；成功時`_SUCCESS`最後原子建立 |

## Scope

- Task type: **B Comprehensive validation**。
- Primary genome scope: chr1-22；chrX/chrY 只列 out-of-scope census。
- Dataset scope: HCC1395、HCC1395_DORADO、COLO829、H1437、H2009、HCC1937、HCC1954。
- Biological scope: 6 biological samples；HCC1395 與 HCC1395_DORADO 是同一細胞株的兩個資料處理版本。
- Claim scope: regional mutation-state trees；bulk molecule evidence 不等同 single-cell clone confirmation。
- Backbone sensitivity scope: 7 datasets 比較 canonical raw-all run 的 LongPhase-S `_sc.vcf` PASS 與同一 paired ClairS 的 PASS **subset**；這可檢驗 LongPhase-S recalibration 對 tree 結果的影響，但不是 independent-caller robustness。舊 SEQC2-HC/DeepSomatic 單樣本比較保留為 historical partial artifact，不進正式 gate。

## Goals

- **G4**：全樣本一致 schema、manifest、checksum、完整驗證與可重跑命令。
- **G5**：分母、候選樹完整性、限制與證據來源可供外部 reviewer 稽核。
- 次要服務 **G1/G3**：保留 HP 與甲基化的正確輔助角色，不把它們寫成獨立確認。

## 教授版全層數據觀察（2026-07-11）

- **主報告 HTML**：`InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html`
- **可重算數據快照**：`InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.data.json`
- **圖表與 claim-to-source mapping**：`InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.source_notes.json`
- **Artifact hard-gate 結果**：`InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.validation.json`
- **瀏覽器 QA**：`InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.browser_qa.json`
- **狀態**：internal PI-share ready；historical artifact validation PASS，scientific release NO-GO。官方 chart contract、14 個 live charts/14 個 no-JS fallback、2 個流程/拓撲 SVG、7 個 dataset standalone details 與 browser/accessibility QA 全部 PASS。報告內數值仍是 historical ClairS-PASS snapshot；production PASS-only probe 已以 `E_METHOD_SCOPE_PASS_ONLY` fail-closed。兩個 normalized raw-all bounded probes與 patch regression已過，但 7-dataset raw-all producer/closeout 尚未啟動完成，故未升為 canonical layered run。
- **正式分類**：`C_region=∏n_trees(primary HP)`、`Topo_region=∏n_distinct_shapes_exact`；可行狀態為 `C=1/Topo=1`、`C>1/Topo=1`、`C>1/Topo>1`，`C=1/Topo>1` 不可能。7-dataset row aggregate：`W_primary=47,377`、complete=`39,885`、incomplete=`7,492`、三狀態=`10,832/11,144/17,909`、impossible=`0`。

### 教授重點 PPT-style HTML（2026-07-11）

- **互動簡報**：`/big7_disk/liaoyoyo2001/Meeting/interSubMod_reports_workspace/20260711_分層重建教授重點簡報_01/slides/20260711_分層重建教授重點簡報_01.standalone.html`
- **講稿與交付說明**：`/big7_disk/liaoyoyo2001/Meeting/interSubMod_reports_workspace/20260711_分層重建教授重點簡報_01/docs/`
- **用途**：10 張結論先行簡報，用 S→W→HP→C→Topo、雙單位漏斗、HP×H3?、complete/incomplete、三類 topology、hidden 邊界與 current gate 在 8–10 分鐘內說明。
- **QA**：artifact/browser/official chart contract/no-JS/print 全 PASS；仍為 internal PI-share，scientific release NO-GO。

## Reproducibility Checklist

- [x] 固定隨機 stress seeds（既有 5 x 800 cases）。
- [x] 鎖定歷史 7 dataset input manifest 與 SHA-256；已識別 scope mismatch。
- [x] PASS-only production probe 依方法 scope fail-closed：`E_METHOD_SCOPE_PASS_ONLY`，未誤產 `_SUCCESS`。
- [x] normalized raw-all bounded probes通過 parsing、record conservation、rescue/remove、zero-read patch regression與 read-tag contracts；pre-decision=`GO_WITH_FAIL_CLOSED`。
- [x] 7-dataset normalized raw-all production sidecar與receipt/closeout通過（2026-07-12 01:31/03:09 +08:00）。
- [x] fresh canonical layered consumer、U3–U7 final verifier、sensitivity/comparison與immutable closeout全部通過（2026-07-13/14）。
- [x] COLO829 bounded real-data direct-vs-sidecar equivalence與synthetic adversarial contract tests通過。
- [x] HCC1954 chr22 PS policy probe：31 retained regions 中 1 mixed-PS；正式 schema 已加入 per-region PS/HP×PS census 與 verifier conservation gate。
- [x] frozen input lock、source bundle、environment lock、launch receipt與real 7-dataset receipt readback全部通過。
- [x] atomic lifecycle：tiny canonical success + empty/missing/existing-root 3-case adversarial matrix通過。
- [x] 保存 params、環境版本、每 stage log 與退出碼。
- [x] 保存五片 mlhp、中間 layered、region view、site ledger 與 verification summary。
- [x] 7/7 full-genome thresholds 全部通過。
- [x] evidence ledger、topic/experiment/output INDEX、docs README 與 CURRENT_FOCUS 已同步至 terminal closeout。
