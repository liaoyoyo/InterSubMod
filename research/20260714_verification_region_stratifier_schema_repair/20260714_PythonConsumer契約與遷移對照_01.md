<!--
建立時間: 2026-07-15
目標: 建立 Verification/RegionStratifier Python consumer 的契約、遷移 disposition、測試 seam 與來源行號對照，供 production/historical reader 實作與審查。
處理範圍: scripts/lib/verification_schema_contract.py、封閉 74-consumer inventory、targeted Python tests、三支任務前既有 untracked historical consumers。
任務類型: (B) Comprehensive validation — bounded inventory 全量 source map；不以 subset 代替 74-consumer 母體。
服務研究目標: G4 reproducibility、G5 外部可驗證工程品質。
證據標記: L1=repo 內可逐行重算的 code/test/TSV；L2=session audit 或 implementation-notes 所記錄、但缺完整 in-repo primary artifact 的邊界。
關聯檔案:
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv
  - InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md
  - InterSubMod/scripts/lib/verification_schema_contract.py
-->

# Python Consumer 契約與遷移對照

**TL;DR — 74/74 個 bounded Python/executable consumers 已有逐列 disposition；production 預設只接受明示且一致的 v2/R1 provenance，unversioned historical H1 必須顯式授權，任何 class 都不得被當成 truth 或 cellular clone。** `[L1]` 來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:148-154`、`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:1-75`。

**敘述框架：SCQA + contract/source map。** Situation 是 v2 taxonomy、legacy view、evidence、LOH、region status 與 pass-through provenance 同時存在；Complication 是舊 consumer 曾混用 current/legacy/truth/clone 語意；Question 是每個入口應讀什麼、何時拒絕；Answer 是下列 C2/L4/E/R1/P/LOH-L/H1 單一 helper 契約與 74-row migration ledger。`[L1]` 來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md:61-73`、`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:78-86`。

## 1. 證據層與閱讀規則

- `[L1]` 表示可直接由 repo 內程式碼、測試或 TSV 逐行驗證；它是本文件的主要證據層。`[L1]` 來源：`InterSubMod/scripts/lib/verification_schema_contract.py:2-6`、`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:1`。
- `[L2]` 只用於 session audit／dirty-worktree 邊界；若沒有完整 before digest，就不能提升成完整 byte-for-byte 證明。`[L2]` 來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:123-129`。
- 本文件不貼出 74 列全文；權威逐檔狀態在 [consumer_migration_status.tsv](./consumer_migration_status.tsv)，其欄位是 `path/inventory_contract/blocking/disposition/status/test_or_guard/notes`。`[L1]` 來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:1`。
- 74 是封閉 inventory 的唯一 executable path 聯集：62 個 VerificationClass paths 加 23 個 LOH paths，扣除 11 個重疊 path；不是一次新的無界全 repo 掃描。`[L1]` 來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:148-153`。

## 2. 七個契約代碼：輸入、語意與最小判斷

| Code | Authoritative input | 最基本判斷與設定 | Unknown / missing 策略 | 可安全宣稱的語意 |
|---|---|---|---|---|
| `[L1]` `C2` | `VerificationClass` + 單一 `VerificationSchemaVersion=2` | v2 enum 固定 11 類；保留固定排序與零頻類；缺 `VerificationClass` 直接拒絕 | 未知或空值轉成 `UnknownCurrentClass` 並保留 count；不折疊為 Noise | current verification taxonomy；不是 truth/evidence boolean |
| `[L1]` `L4` | `VerificationClass_Legacy` | 只接受 `Strong/Subclone/Weak/Noise`；canonical legacy 欄一旦存在，任何未知或 missing 立即拒絕 | canonical mode 無容錯 bucket；H1 才能另選 `fail/exclude` | historical four-state cohort；不是 current v2 taxonomy |
| `[L1]` `E` | 4 typed booleans + `EvidencePath` + `EvidenceDerivation` | `LabelFirstSupport/ClusterFirstSupport` 不可 NA；`WithinHPSupport/DispersionWarning` 可保留 NA；path/derivation 必須命中 enum | 非 `true/false`、不允許的 NA、未知 path/derivation 全部拒絕 | 明示 evidence axes；禁止從 class 名稱反推 |
| `[L1]` `R1` | `RegionStratificationSchemaVersion=1`、ID/Label/Reason + status row | status 欄完整；status/版本/count/sentinel/label/reason/alias 互相一致；非 `VALID` 必須零 assignment | 不存在「猜測」路徑；任何矛盾 fail closed | region-stratum assignment/status；不是 cellular clone |
| `[L1]` `LOH-L` | `LOH_Subtype_LegacyVC` + deprecated exact alias `LOH_Subtype` | canonical 欄存在時 deprecated alias 也必須存在且逐列完全相同；enum 固定 | unknown/missing/alias mismatch 全部拒絕 | legacy-VerificationClass-derived LOH subtype |
| `[L1]` `P` | 完整 12 欄 `VERIFICATION_PROVENANCE_COLUMNS` | 同時驗 C2、L4、E、LOH-L，以及 current↔deprecated↔legacy↔evidence 的列級關係 | P 比 C2 更嚴：current unknown 即拒絕；欄缺、混版、關係不一致皆拒絕 | derived table 可攜、可回溯的完整 provenance bundle |
| `[L1]` `H1` | 無 `VerificationSchemaVersion` 的 frozen historical input | 只能由 `--allow-unversioned-v1` 或等價顯式參數啟用；有 version 卻缺 canonical legacy 欄時禁止 fallback | L4 可選 `fail` 或 `exclude`；LOH-L 固定 fail；輸出/metadata 必須留 `UNVERSIONED_V1` 與 warning | historical replication compatibility；不可升格為 production v2 證據 |

`[L1]` 上表 authority：C2/L4/E/R1/LOH-L/P/H1 的原始規格對照在 `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md:63-71`；實作 authority 分別在 `InterSubMod/scripts/lib/verification_schema_contract.py:18-126`、`:196-292`、`:321-429`、`:432-573`、`:576-684`。

### 2.1 C2 的 11-class v2 設定

- `[L1]` 固定順序為 `Strong_Bidirectional`、`ClusterFirstOnly`、`LOH-Structure`、`MultiGroupNoLabel`、`LabelShift`、`PermanovaLocation`、`StructureNoLabel`、`DispersionStructure`、`Noise_Uniform`、`Noise_Chaotic`、`Noise_Uncorrelated`；machine value 不得用 display alias 取代。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:21-33`、`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md:75-81`。
- `[L1]` `select_current_view()` 缺 version 時預設拒絕；只有 `allow_unversioned_raw=True` 才回傳 `UNVERSIONED` raw view，而且明示「without v2 interpretation」。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:196-217`。
- `[L1]` v2 mode 強制 version 恰為整數 2；未知及 missing 都正規化到 `UnknownCurrentClass`，並分別計數。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:184-193`、`:219-233`。

### 2.2 L4 與 H1 的分界

- `[L1]` canonical L4 只讀 `VerificationClass_Legacy`，未知或 missing 直接 `SchemaContractError`；不能由 current `VerificationClass` 猜回四態。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:236-265`。
- `[L1]` H1 同時要求「顯式授權」與「輸入沒有 VerificationSchemaVersion」；versioned input 缺 legacy 欄仍必須拒絕。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:261-271`。
- `[L1]` H1 unknown policy 只有 `fail` 與 `exclude`；`exclude` 會保留 unknown count、將不合法值變為 NA，且發出 `UNVERSIONED` warning。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:243-244`、`:272-292`。

### 2.3 E 的 typed evidence 判斷

- `[L1]` boolean parser 只接受 Python bool 或 exact lowercase `true/false`；字串 `NA`／missing 只有在欄位契約允許時才保留為 nullable boolean。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:321-345`。
- `[L1]` `LabelFirstSupport` 與 `ClusterFirstSupport` 不允許 NA；`WithinHPSupport` 與 `DispersionWarning` 允許 NA，避免 historical/offline 未分析被誤記為 false。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:348-366`。
- `[L1]` `EvidencePath` 與 `EvidenceDerivation` 均採封閉 enum；未知或 missing 全部拒絕。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:71-84`、`:368-385`。

### 2.4 LOH-L 的 canonical/alias 判斷

- `[L1]` canonical input 是 `LOH_Subtype_LegacyVC`；`LOH_Subtype` 是 required deprecated exact alias，不是另一個可獨立選擇的 truth。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:388-412`。
- `[L1]` `None/LOH_Noise/LOH_Weak/LOH_Strong/LOH_Subclone` 以外或 missing 直接拒絕；canonical 與 alias 任一列不等也直接拒絕。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:86-92`、`:396-411`。
- `[L1]` deprecated-only fallback 只允許 H1 且輸入不得帶 VerificationSchemaVersion；成功仍標 `UNVERSIONED_V1` 並 warning。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:414-429`。

### 2.5 R1 的 status、sentinel 與 count 判斷

- `[L1]` status 必含 `status/reason/schema_version/eligible_region_count/min_regions_required/assignment_count/n_occupied_region_strata/warning_count/generated_at`，status 只接受 `VALID/INSUFFICIENT_REGIONS/NOT_APPLICABLE_TUMOR_ONLY/FAILED`。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:94-99`、`:438-454`。
- `[L1]` version/count 必須是非負整數；version 必須是 1；reason/generated_at 不可空；`n_occupied_region_strata <= assignment_count`。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:456-486`。
- `[L1]` `VALID` 要求 assignment count 等於 eligible count、eligible 非零且達最低門檻；`INSUFFICIENT_REGIONS` 不可與已達門檻的 eligible count 矛盾。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:487-504`。
- `[L1]` row 必含 schema version、`RegionStratum_ID/Label/Reason`；ID 只接受 `-1,0,1,2,3`，0–3 的 label/reason 固定，`-1/Unassigned` reason 由 status 決定。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:100-111`、`:506-545`。
- `[L1]` deprecated `Subclone_ID` 若存在只能是 exact alias；實際 assignment/occupied count 必須等於 status，非 `VALID` 必須兩者皆為零。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:547-564`。

### 2.6 P 的完整 provenance 關係

- `[L1]` P 固定複製 12 欄：version、current、v1 deprecated、legacy、四個 evidence booleans、EvidencePath、EvidenceDerivation、canonical/deprecated LOH。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:113-126`。
- `[L1]` versioned P 依序執行 C2、L4、E、LOH-L，要求完整欄；若 C2 有 unknown，P 不予輸出。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:599-608`。
- `[L1]` P 逐列檢查 current↔v1 deprecated、current↔EvidencePath、current↔legacy，以及 label/cluster support↔legacy；`LEGACY_CLASS` derivation 還要求 WithinHP/Dispersion 為 NA。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:610-680`。
- `[L1]` 全部通過後才複製固定 12 欄並標 `schema_status=V2`；unversioned H1 只回傳可驗的 legacy/optional deprecated LOH 子集，不冒充完整 v2 P。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:576-597`、`:682-684`。

## 3. Shared helper 入口與拒絕矩陣

| 入口 | Production 呼叫 | Historical 呼叫 | 主要拒絕條件 | Unknown policy | 精確來源 |
|---|---|---|---|---|---|
| `[L1]` `select_current_view` | `select_current_view(df)` | 僅 raw 顯示可用 `allow_unversioned_raw=True` | 缺 class；缺/混/非 2 version | v2 → bucket+count；raw → 不宣稱 v2 | `InterSubMod/scripts/lib/verification_schema_contract.py:196-233` |
| `[L1]` `select_legacy_view` | `select_legacy_view(df)` | `allow_unversioned_v1=True`，並明選 `fail/exclude` | canonical legacy 缺欄或非法；versioned fallback | canonical fail；H1 fail/exclude | `InterSubMod/scripts/lib/verification_schema_contract.py:236-292` |
| `[L1]` `ordered_class_crosstab` | selector 回傳的 `ClassView` 直接餵入 | 同 production；不自行猜 schema | label/view index 不同；normalized table 要求 margins | 不改 unknown；依 view categories 保留零頻類 | `InterSubMod/scripts/lib/verification_schema_contract.py:295-318` |
| `[L1]` `read_evidence` | `read_evidence(df)` | 只有已具 typed evidence 欄才可呼叫 | required 欄缺；boolean/path/derivation 非法 | 不設 unknown bucket，全部 fail | `InterSubMod/scripts/lib/verification_schema_contract.py:321-385` |
| `[L1]` `select_loh_legacy` | `select_loh_legacy(df)` | `allow_unversioned_v1=True` | canonical/deprecated 缺欄、alias 不等、非法 enum、versioned fallback | 全部 fail | `InterSubMod/scripts/lib/verification_schema_contract.py:388-429` |
| `[L1]` `validate_region_strata` | `validate_region_strata(df, status)` | H1 不可假裝有 R1；應明示 unavailable | status/schema/row/sentinel/count/alias 任一矛盾 | 全部 fail | `InterSubMod/scripts/lib/verification_schema_contract.py:432-573` |
| `[L1]` `extract_provenance_frame` | `extract_provenance_frame(df)` | `allow_unversioned_v1=True` 只產有限 historical provenance | index 非唯一、欄缺、混版、unknown、列級關係矛盾 | versioned 全部 fail；H1 走 L4/LOH-L | `InterSubMod/scripts/lib/verification_schema_contract.py:576-684` |

- `[L1]` helper 故意不掌握 InterSubMod data path；caller 先載入 DataFrame/status，再交由 helper 驗證，避免 discovery fallback 把錯檔當輸入。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:2-6`。
- `[L1]` 所有 fail-closed 契約共用 `SchemaContractError`；缺 required columns 會明列缺欄。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:129-130`、`:160-163`。
- `[L1]` C2/L4/LOH-L selector 回傳 `ClassView`，metadata 固定含 selection field、schema status、categories、unknown counts 與 warnings；consumer 不需重造一套 provenance 字典。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:133-149`。

## 4. 74 consumers 的 disposition 分類

以下 count 是以 `disposition` 欄對 rows 2–75 做 deterministic aggregation；16 類互斥且總和為 74。完整逐檔 contract/blocking/test seam 請直接查 [consumer_migration_status.tsv](./consumer_migration_status.tsv)。`[L1]` 來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:1-75`。

| # | Count | Disposition 摘要 | 代表檔（不展開全部 74 列） | Source row |
|---:|---:|---|---|---|
| `[L1]` 1 | 12 | typed E 或明示 schema view；missing/unknown fail 或獨立保留 | `analyze_clairs_borderline_fn.py`、`analyze_gq_methylation_rescue_matrix.py`、`validate_method_design.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:12` |
| `[L1]` 2 | 12 | canonical LOH-L + alias equality + historical gate | `20260423_B4_S4_discrimination.py`、`build_loh_feature_validation.py`、`run_phase1a_round3_loh_feature_test.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:3` |
| `[L1]` 3 | 11 | 明示 L4 historical cohort，必要時 H1 + unknown policy | `analyze_cross_sample_ism_af_gradient.py`、`compare_subclone_validation.py`、`tools/find_verification_candidates.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:13` |
| `[L1]` 4 | 10 | C2 v2 固定 11 類 + unknown bucket | `analyze_candidate_rules.py`、`ari_typology.py`、`tools/compare_vcf_results.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:11` |
| `[L1]` 5 | 10 | P/pass-through 保留 schema identity；拒絕不完整或混版 | `20260424_X6_merge_caller_af_S3S5.py`、`build_phase1_training_manifest.py`、`summarize_phase1a_prediction_failures.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:10` |
| `[L1]` 6 | 4 | historical LOH reader 加 canonical LOH + H1 + provenance | `20260423_B3_paired_obs18.py`、`20260423_B5_colo829_s1_fold.py`、`20260423_B7_loh_noise_signal.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:2` |
| `[L1]` 7 | 4 | remaining pass-through adapter 保留完整 v2 provenance | `export_phase1_manifest_shard.py`、`export_phase1_read_training_table.py`、`run_to_support_feature_diagnostics.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:53` |
| `[L1]` 8 | 3 | 移除 VerificationClass-derived truth；要求 explicit truth 欄 | `build_allele_deep_dive.py`、`build_feature_direction_map.py`、`build_non_loh_discrimination.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:23` |
| `[L1]` 9 | 1 | X5：LOH-L + H1 + fail-closed + per-input JSON provenance | `20260424_X5_crosssample_obs18_S3S5.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:8` |
| `[L1]` 10 | 1 | 純文字 PI 圖說標成 historical v1；不新增 runtime selector | `generate_pi_report_figures_self_phasing.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:58` |
| `[L1]` 11 | 1 | frozen golden schema-aware；stale/non-v2 在 C++ 前拒絕，自動 rebuild 禁止 | `scripts/regression/regression_check.sh` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:72` |
| `[L1]` 12 | 1 | metrics 驗 C2/E、unknown bucket、mixed/implicit schema | `scripts/validation/metrics/compute_metrics.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:73` |
| `[L1]` 13 | 1 | derived figures 要求 upstream taxonomy metadata；H1 必須本地明示 | `tp_fp_structure_label_figures.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:69` |
| `[L1]` 14 | 1 | HTML 要求 upstream schema_contract + region status | `subcluster_split_allsamples/build_hcc1395_html.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:65` |
| `[L1]` 15 | 1 | phase2 annotation 使用 explicit support mode + schema-aware loader | `build_phase2_annotation_layer.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:43` |
| `[L1]` 16 | 1 | HCC accounting 同時驗 verification/LOH/R1 完整 provenance | `subcluster_split_allsamples/hcc1395_full_accounting.py` | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:66` |

`[L1]` 算術閉合：12+12+11+10+10+4+4+3+8×1 = 74；ledger 的 74 列均為 `MIGRATED`、0 unresolved。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:148-154`、`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:2-75`。

## 5. Targeted test → contract/source mapping

| Test source | 覆蓋的 seam | 對應 production risk |
|---|---|---|
| `[L1]` `tests/test_verification_schema_consumers.py` | C2 全 11 類/unknown/零頻/unversioned raw；L4/H1；E typed NA；LOH alias/H1；R1；完整 P；candidate/replication selector | shared helper 的 table-driven contract |
| `[L1]` `tests/test_verification_schema_b0_evidence.py` | phase2 v2/legacy equivalence、missing/NA/unknown fail、非 analyzed row 保持 NA、兩 evidence axes、annotation loader | 防止 class-name evidence inference |
| `[L1]` `tests/test_verification_schema_b0_clean_consumers.py` | C2 unknown、L4 required、LOH alias、H1 unknown、truth/evidence decision tree、mixed version | B0 cohort/decision 錯選 |
| `[L1]` `tests/test_verification_provenance_b2_consumers.py` | adapters fail closed/H1 metadata、12 欄 P、mixed version/status、LOH alias、derived metadata | B2 pass-through 丟 provenance |
| `[L1]` `tests/test_verification_schema_b1_guards2.py` | phase1/metrics/figures taxonomy metadata、direct evidence、H1 local flag、frozen golden refuse | B1 圖表/metrics/golden schema 漂移 |
| `[L1]` `tests/test_verification_schema_historical_consumers.py` | explicit L4；H1 必須 flag 且保留 unknown bucket | historical association/replication 邊界 |
| `[L1]` `tests/test_verification_schema_loh_round3.py` | canonical LOH model feature、alias、H1、missing/unknown 不得變 None | LOH cohort 靜默誤分 |
| `[L1]` `scripts/tests/test_remaining_pass_historical_consumers.py` | 4 個 pass-through consumer 的 complete P/mixed/missing；4 個 historical LOH consumer 的 explicit H1；禁止 raw deprecated LOH | remaining-pass migration closeout |

`[L1]` 上表精確 test seam：`InterSubMod/tests/test_verification_schema_consumers.py:39`、`:55`、`:76`、`:86`、`:101`、`:128`、`:152`、`:192`、`:231`、`:280`、`:303`、`:326`；`InterSubMod/tests/test_verification_schema_b0_evidence.py:58`、`:70`、`:95`、`:108`、`:118`、`:125`、`:149`、`:161`、`:215`；`InterSubMod/tests/test_verification_schema_b0_clean_consumers.py:70`、`:92`、`:113`、`:124`、`:150`、`:172`、`:198`、`:219`。

`[L1]` provenance/guard/historical/LOH 精確 test seam：`InterSubMod/tests/test_verification_provenance_b2_consumers.py:94`、`:115`、`:135`、`:168`、`:191`、`:200`、`:265`；`InterSubMod/tests/test_verification_schema_b1_guards2.py:69`、`:80`、`:94`、`:106`、`:115`、`:168`、`:182`、`:196`、`:236`、`:247`、`:260`、`:273`；`InterSubMod/tests/test_verification_schema_historical_consumers.py:45`、`:56`；`InterSubMod/tests/test_verification_schema_loh_round3.py:26`、`:49`、`:67`；`InterSubMod/scripts/tests/test_remaining_pass_historical_consumers.py:81`、`:111`、`:135`、`:161`、`:186`。

`[L1]` durable execution record 顯示 B0-Clean、shared/evidence/historical 分批 targeted suites 已通過；本 source map 的責任是將測試 seam 指回原始碼，不以文件文字取代重新執行。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:131-137`。

## 6. Production 與 historical 的使用說明

### 6.1 Production/default

1. `[L1]` current 分析先 `select_current_view(df)`；缺 version、混版或非 2 立即停止，未知 class 可保留在 `UnknownCurrentClass` 做 diagnostics，但不可默認 Noise。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:196-233`。
2. `[L1]` 需要 support/rescue 判斷時直接 `read_evidence(df)`，不從 `VerificationClass` 名稱重建 label-first/cluster-first/within-HP/dispersion。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:348-385`。
3. `[L1]` 歷史四態 panel 即使在 v2 檔也只讀 explicit `VerificationClass_Legacy`；LOH 只讀 `LOH_Subtype_LegacyVC` 並驗 deprecated alias。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:246-259`、`:394-412`。
4. `[L1]` derived/export table 用 `extract_provenance_frame(df)`，不能只帶裸 `VerificationClass`；P 會在 unknown、缺欄或關係不一致時拒絕。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:599-684`。
5. `[L1]` region report 必須同時載入 status row 並呼叫 `validate_region_strata(df,status)`；只有 status 與 row/count/sentinel 全一致才可顯示 R1 語意。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:432-573`。
6. `[L1]` production default 不開 `allow_unversioned_*`；migration ledger 明確記錄 production v2 default fail-closed。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:148-153`。

### 6.2 Historical replication/H1

1. `[L1]` 只有確認輸入確實是 unversioned frozen v1 時，才由 CLI `--allow-unversioned-v1` 傳入 helper；任何帶 version 的 input 不得走 H1。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:261-281`。
2. `[L1]` 四態 cohort 使用 `select_legacy_view(..., allow_unversioned_v1=True, unversioned_unknown_policy='fail'|'exclude')`；`exclude` 必須保留 unknown count，不得把 unknown 變 Noise。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:236-292`。
3. `[L1]` deprecated LOH-only historical input 使用 `select_loh_legacy(..., allow_unversioned_v1=True)`；unknown/missing 仍 fail，不得 `fillna('None')`。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:414-429`。
4. `[L1]` H1 current raw view只可做 raw 值展示；其 metadata 是 `UNVERSIONED`，不是 v2 C2。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:202-217`。
5. `[L1]` H1 不具 R1 時應明示 unavailable，不能解讀 deprecated `Subclone_ID`；HCC accounting 的 historical branch即回報 `UNVERSIONED_REGION_SCHEMA/HISTORICAL_SUBCLONE_ID_NOT_INTERPRETED`。來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/hcc1395_full_accounting.py:80-94`。

## 7. 三支既有 untracked historical consumer：before-hash 邊界

`[L2]` 這三檔在任務啟動前已是 user-owned untracked artifacts；接受的改動邊界只有局部 schema/status/legacy/LOH guards，不改資料路徑、門檻或統計方法。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:123-129`。

| 檔案 | Before SHA-256 | 目前 modified SHA-256 | 可驗的局部邊界 |
|---|---|---|---|
| `[L2]` `scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py` | session audit 只留 prefix `4b5d8d3…`；完整 digest unavailable | `[L1]` `43b233ed39cd2c72e599690922d1206a633daff5adb9a240cc8ba423fce8660a` | JSON 必須有 schema contract；class 不重建 truth；current/legacy label 不稱 cellular clone。來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:14-26`、`:44-60`、`:261-265`、`:290-302` |
| `[L2]` `scripts/analysis/subcluster_split_allsamples/hcc1395_full_accounting.py` | session audit 只留 prefix `af386d0…`；完整 digest unavailable | `[L1]` `289add1c0c42319fa8221286f7f8c56f760805f31abe39641c9d27ade9e0683b` | production 驗 C2/L4/LOH-L/P/R1；H1 必須 flag，R1 明示 unavailable。來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/hcc1395_full_accounting.py:19-25`、`:39-41`、`:63-96` |
| `[L2]` `scripts/analysis/tp_fp_structure_label_association.py` | session audit 只留 prefix `fcec9b0…`；完整 digest unavailable | `[L1]` `ce428678c48fcc2985573389ba1a14d562243748b81d12cb0affc75965aa0233` | historical panel 明示 L4/H1；TP/FP schema一致；within-HP 直接讀 canonical row field，不 parse deprecated stub。來源：`InterSubMod/scripts/analysis/tp_fp_structure_label_association.py:11-17`、`:79-123`、`:152-154` |

- `[L2]` before prefixes 來自 session audit，repo 內 implementation-notes 只記錄「原始與修改後 SHA 均留存」，沒有可引用的完整 before digest；因此 prefix 只能作人工識別，不足以作 cryptographic byte-preservation proof。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:126-129`。
- `[L1]` 目前 modified SHA-256 是對上表三個現存檔案完整 bytes 的重算值；它只能鎖定「目前版本」，不能補回缺失的完整 before digest。檔案來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:1`、`InterSubMod/scripts/analysis/subcluster_split_allsamples/hcc1395_full_accounting.py:1`、`InterSubMod/scripts/analysis/tp_fp_structure_label_association.py:1`。

## 8. 語意硬界線：class 不是 truth，也不是 cellular clone

### 8.1 不可把 class 當 truth

- `[L1]` `VerificationClass` 是 taxonomy/annotation，不是 binary truth；三支高風險 consumer 現在都要求 explicit truth field，否則拒絕。來源：`InterSubMod/scripts/analysis/build_allele_deep_dive.py:33-38`、`InterSubMod/scripts/analysis/build_feature_direction_map.py:54-59`、`InterSubMod/scripts/analysis/build_non_loh_discrimination.py:50-55`。
- `[L1]` HTML 的 `Significant` 直接讀 canonical boolean；`DispersionStructure`/Noise/current class 只作 taxonomy 描述，不再重建 truth。來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:261-265`。
- `[L1]` E 的支持判斷必須讀 typed fields；class 名稱不能替代 `LabelFirstSupport` 或 `ClusterFirstSupport`。來源：`InterSubMod/scripts/lib/verification_schema_contract.py:348-385`、`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:131-137`。

### 8.2 不可把 class/region-stratum 當 cellular clone

- `[L1]` `ClusterFirstOnly` 只表示 cluster-first evidence 成立、label-first 不成立；不是 cellular clone。來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:44-46`。
- `[L1]` legacy `Subclone` 是 historical four-state label／within-HP multi-group candidate，名稱不能當 normal-anchored cellular subclone 證據。來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:57-60`、`:295-297`。
- `[L1]` R1 的 `RegionStratum_ID/Label/Reason` 是 region stratification assignment；`-1` 必須依 status/reason 解讀，deprecated `Subclone_ID` 只可驗 exact alias，不能沿用「無 subclone」語意。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_inventory.md:228-234`、`InterSubMod/scripts/lib/verification_schema_contract.py:521-564`。
- `[L1]` 無監督 `cansplit k>=2`、within-HP multigroup、Verification evidence 互不蘊含 cellular clone。來源：`InterSubMod/scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:295-302`。

## 9. 完整性檢核與已知限制

| 檢核項 | 結果 | Evidence |
|---|---|---|
| `[L1]` 七個代碼 C2/L4/E/R1/P/LOH-L/H1 | 7/7 已定義入口、拒絕條件與 unknown policy | `InterSubMod/scripts/lib/verification_schema_contract.py:196-684` |
| `[L1]` bounded consumer disposition | 16 類、count 總和 74；74 `MIGRATED`、0 unresolved | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:1-75` |
| `[L1]` targeted test/source map | 8 個 test modules；shared/B0/B1/B2/historical/LOH/remaining-pass seams 均有 source line | `InterSubMod/tests/test_verification_schema_consumers.py:39`、`InterSubMod/scripts/tests/test_remaining_pass_historical_consumers.py:81` |
| `[L1]` production vs historical | production v2 fail-closed；H1 explicit opt-in且留 metadata/warning | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:148-153` |
| `[L1]` truth/clone hard boundary | explicit truth guard + HTML/RegionStratum wording guard 已定位 | `InterSubMod/scripts/analysis/build_allele_deep_dive.py:33-38`、`InterSubMod/scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:261-302` |
| `[L2]` untracked before-hash | 3/3 只有 before prefix；完整 before digest 不在 repo durable source，不能宣稱完整 byte-preservation proof | `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:123-129` |

- `[L1]` inventory closure 是本任務既定 62 VC/23 LOH readers 的聯集；本文件沒有用大量新搜尋擴張母體，符合「以關聯與既有 inventory 查找」的限制。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:148-153`。
- `[L1]` frozen corpus 能證明 schema/serialization regression 邊界，不能單獨證明 v1 upstream classification 或 biology；consumer source map 同樣不得把 engineering PASS 擴寫為生物 truth。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:232-236`。

---

**Source-map verdict：COMPLETE WITH ONE EXPLICIT L2 LIMITATION。** `[L1]` C2/L4/E/R1/P/LOH-L/H1、74 dispositions、test seams、production/H1 使用方式與 truth/clone 禁線均已逐行映射；`[L2]` 唯一未能提升為 L1 的項目是三支既有 untracked consumer 缺完整 before SHA-256，已明示為不可回推邊界。來源：`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv:1-75`、`InterSubMod/research/20260714_verification_region_stratifier_schema_repair/implementation-notes.md:123-154`。
