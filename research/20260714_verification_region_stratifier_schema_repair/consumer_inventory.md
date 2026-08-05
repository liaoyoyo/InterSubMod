<!--
建立時間: 2026-07-14
目標: 唯讀盤點 Verification schema v2、RegionStratification schema v1 與 LOH provenance 的 executable downstream consumers，提供 P1-B 遷移清單
處理範圍: InterSubMod/tools、InterSubMod/scripts、InterSubMod/tests；排除 docs、archive/deep、資料與既有輸出
關聯檔案:
  - /big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md
  - InterSubMod/tools/compare_vcf_results.py
  - InterSubMod/tools/find_verification_candidates.py
  - InterSubMod/scripts/analysis/analyze_methyl_cluster_allele_cooccurrence.py
  - InterSubMod/scripts/analysis/compare_subclone_validation.py
狀態: read-only source audit；尚未代表 consumers 已修復
服務目標: G4 reproducibility、G5 external verifiability
-->

# Verification / RegionStratifier consumer inventory

## 0. TL;DR 與任務邊界

**結論：P1-B 的實際 executable-code surface 明顯大於規格列出的四支已知 consumer。** 在限定 `tools/`、`scripts/`、`tests/` 且排除獨立文件、資料與 archive 後，共找到下列 token-hit files；其中包含 runtime selector、dynamic reader、pass-through projection與少數只在 executable docstring/plot文字命中的檔案，並不把62支全部冒稱為 active runtime readers：

| Schema token | Executable-code files containing token | 判讀 |
|---|---:|---|
| `VerificationClass` | 62 | 只有 2 支已讀 `VerificationClass_Legacy`；多支仍固定四態或用 class 推導 evidence/truth，少數是pass-through/文字命中 |
| `VerificationClass_Legacy` | 2 | 皆位於 `subcluster_split_allsamples/`，不是共用 schema guard |
| `LOH_Subtype` | 23 | live v2 雖保留 exact deprecated alias，但 consumer provenance 仍未顯式化 |
| `Subclone_ID` / `RegionStratum` | 2 | 都在 HCC1395 accounting/HTML；仍用舊 biological wording；尚無 canonical `RegionStratum_*` reader |
| `n_subclones` / `n_occupied_region_strata` | 0 | 未找到 Python/shell executable consumer |
| `subclone_structure.txt` | 1 | 只有註解引用，沒有 parser；註解仍把 exact state 說成在舊 artifact |
| `tests/` 對上述 schema tokens | 0 | 目前沒有 Python/shell targeted regression 接住 consumer migration |

這是 **Task Type E（schema hotfix）內的 comprehensive consumer audit**，不是大型 BAM、F1 或生物學重驗。沒有讀取 frozen corpus data rows，沒有執行 consumer，也沒有修改 consumer。`truth_total/calls_total/TP/FP/FN/precision/recall/F1` 對本次 source audit 均為 N/A。

## 1. 啟動問題、假設與驗收

### 1.1 本次研究問題

Stage④ v2 與 RegionStratifier 修復落地後，哪些 executable consumers 會靜默：

1. 把 11-class current taxonomy 當成舊四態；
2. 把 `Strong/Subclone` 重新推成 evidence booleans；
3. 把 legacy region stratum 當 cellular subclone；
4. 繼續用無 provenance 的 `LOH_Subtype`；
5. 對未知 class 自動歸入 Noise、丟棄、或產生空 cohort。

### 1.2 關鍵假設

- Canonical contract 以外部修復規格 §4 P0-A/P0-B/P0-C/P1-A/P1-B 為準；字串與 precedence 不自行改名。
- 「固定 historical replication」可讀 unversioned v1，但必須由 explicit flag 開啟、寫 warning/metadata；不得把 fallback 當 canonical v2。
- 對新 canonical output 執行的 consumer 必須驗 `VerificationSchemaVersion=2`；current 未知值進 `UnknownCurrentClass`，legacy 缺欄 fail closed。
- `LOH_Subtype` 在 v2 是 deprecated exact alias；新分析應讀 `LOH_Subtype_LegacyVC`，只有 historical unversioned mode 可顯式讀舊欄。
- 不把一般文字中的 `Strong`/`Noise` 當 consumer；inventory 必須同時命中 schema 欄名或位於該欄位的實際決策路徑。

### 1.3 Step → Verify

1. 精確讀修復規格 P1-B 與四支指定 consumer → 驗證：每支有 file:line、old assumption、canonical source、fallback、test seam。
2. 對 `tools scripts tests` 做 token-level `rg` → 驗證：62/23/2/0/1 的 file counts 可重算，沒有 docs/ 資料列。
3. 將 consumer 分成 current / legacy / evidence / region-stratum / historical-only → 驗證：每列有 mode、blocking 與 unknown policy。
4. 檢查固定 input headers 與 CLI call sites → 驗證：兩個固定 historical inputs 都是 unversioned；`run_batch_vcf_analysis.sh` 是已知 compare caller。
5. 重驗四支 consumer SHA-256 與 targeted git diff → 驗證：四支 consumer 無本次工作樹修改。

## 2. 分類碼與共用 fallback 契約

| Code | Canonical source | Required behavior |
|---|---|---|
| `C2` | `VerificationClass` + `VerificationSchemaVersion=2` | 11 類固定順序；未知值 → `UnknownCurrentClass` 並計數；零頻類仍保留 |
| `L4` | `VerificationClass_Legacy` | 只允許 `Strong/Subclone/Weak/Noise`；缺欄 fail closed；未知值 → `UnknownLegacyClass` 並警告/停止 cohort |
| `E` | `LabelFirstSupport`、`ClusterFirstSupport`、`WithinHPSupport`、`DispersionWarning`、`EvidencePath` | 不再從 class 名稱反推 evidence；offline `NA` 不得當 false |
| `R1` | `RegionStratum_ID/Label/Reason` + `RegionStratificationSchemaVersion=1` | `-1/Unassigned` 依 status 解讀；不得稱 cellular clone |
| `LOH-L` | `LOH_Subtype_LegacyVC` | 新 canonical 必須讀此欄；若 deprecated `LOH_Subtype` 同時存在要 assert exact equality |
| `P` | pass-through provenance | export/manifest 至少帶 schema version、current、legacy、evidence path/booleans；不能只複製裸 current class |
| `H1` | unversioned historical input | 只有 explicit `--allow-unversioned-v1` 或等價 mode 才允許；輸出寫 `UNVERSIONED`、來源欄與 warning |

Blocking 分級：`B0` = 會選錯 cohort、推錯 evidence/truth 或在 v2 靜默變空；`B1` = 圖/統計遺漏或 schema 混用；`B2` = 主要是 provenance/pass-through；`H` = 固定歷史腳本，不阻塞 live pipeline，但必須 guard 才能宣稱 P1-B 完整。

### 2.1 Exact enum authority（不可由 display 名稱反向污染 code）

本 audit 額外發現外部比較文件 `20260714_InterSubMod與subclone重建方法全面比較_01.md:314-315` 使用 display 名稱 `ClusterOnlyAssociation`，但鎖定修復規格 `20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md:185,207-210,236-250,581` 明定 canonical code value 是 **`ClusterFirstOnly`**、evidence value是 **`CLUSTER_FIRST_ONLY`**。因此：

- shared helper、C++ enum、CSV、tests與consumer不得接受或輸出 `ClusterOnlyAssociation` 作 canonical alias；
- 若報表想用 `Cluster-only association` 作 display label，必須保留 machine value `ClusterFirstOnly`，且 display mapping不可回寫資料；
- `rg` 在 scoped executable files目前沒有找到 `ClusterOnlyAssociation`、`ClusterFirstOnly`、`Strong_Bidirectional` 或 `VerificationSchemaVersion`，符合「修復尚未落地」，也代表後續第一次實作必須由單一 helper提供 exact constants。

## 3. 四支規格鎖定 consumer：逐檔最小修法

### 3.1 `tools/compare_vcf_results.py`

**證據位置**：CLI `49-54`；輸入與 recursive fallback `57-91`；raw crosstab `143-166`；固定四態 plot `244-272`；heatmap `340-359`；report `527-580`；batch caller 是 `scripts/run_batch_vcf_analysis.sh:207-237`。

- 舊假設：`class_order=['Strong','Subclone','Weak','Noise']`（`255-260`）。v2 的 11 類會在 stacked plot 被全部丟掉；table/heatmap雖會接受 raw value，卻沒有 schema version、unknown bucket 或固定零頻類。
- Canonical 行為：預設 current `C2`；集中建立 `verification_view` 欄，v2 值先 validate 再以 11 類 + `UnknownCurrentClass` 固定 reindex。legacy plot 只能由 optional `--verification-view legacy` 明確選擇。
- Fallback/unknown：schema version 缺失時只畫 raw values，所有表/標題/metadata 寫 `UNVERSIONED`；不可套四態排序猜版本。v2 未知字串歸 `UnknownCurrentClass`，保留原值於 diagnostics。
- 最小改法：保留既有 `--output-dir/--labels/--paths`，新增 **optional** view flag；將 stats、stacked plot、heatmap、report 全部吃同一 validated view，避免四處各自 crosstab。
- Test seam：對純 helper 餵 11 classes + 1 unknown + 零頻 label，assert table columns/legend 都完整；另測 unversioned raw values 標記。既有 batch caller不可因新增 required arg 破壞。
- 相鄰 I/O 風險：`67-75` 會在 exact file 缺失時 `os.walk` 並拿第一個結果，與本任務「精確關聯查找」原則衝突。最小做法是 fail exact path；若保留 discovery，需 opt-in、排序且多結果時 fail ambiguous。

### 3.2 `tools/find_verification_candidates.py`

**證據位置**：CLI `7-11`；四態 selector `31-37`；report title `39-44`；row label `94-99`；required columns `102-123`。

- 舊假設：直接從 current `VerificationClass` 找 `Weak/Subclone/Strong/Noise`。
- Canonical 行為：這份既有四態報告固定讀 `L4`，標題改為 `Legacy verification candidates`；row metadata 也明示 `VerificationClass_Legacy`。
- Fallback/unknown：缺 `VerificationClass_Legacy` 立即以清楚訊息退出，不可由 current/v1 deprecated 猜；未知 legacy 值列出 count 後 fail，不納入四類。
- 最小改法：只把 selector/required field/report label 集中成一個 `class_field` 常數；不必新增 current mode。若新增，必須用 CLI flag 並走獨立 `C2` path。
- Test seam：四列 fixture（legacy 四類各一）+ unknown + missing-column；capture stdout，四組各命中 1，title/欄名正確，兩個錯誤皆非零退出。
- CLI 風險：未找到 repo 內 caller；既有 unversioned v1 使用者會從「可跑」變成 fail，這是規格要求的 intentional break。

### 3.3 `scripts/analysis/analyze_methyl_cluster_allele_cooccurrence.py`

**證據位置**：固定 input/output `26-45`；全域 warning suppression `24`；result class `170-175`；selector `222-240`；所有圖/文件四態敘述 `278-448`；entry point `456-482`。

- 舊假設：固定 historical HCC1395 summary 的 current `VerificationClass` 直接做 Strong vs Noise。
- Canonical 行為：historical replication 固定讀 `L4` 的 `Strong` vs `Noise`，輸出 metadata 寫 `selection_field=VerificationClass_Legacy` 與 values；不得把 v2 classes折疊。
- Fallback/unknown：只有 `--allow-unversioned-v1` 才可在固定舊檔讀裸 `VerificationClass`；輸出 `UNVERSIONED` warning與 unknown/excluded counts。預設無 flag 必須 fail。
- 實際相容風險：目前硬編碼 input header只有 `VerificationClass`/`LOH_Subtype`，沒有 schema version或 legacy欄；因此現有 no-arg invocation 在修復後會 intentional fail，使用者必須顯式加 flag。`warnings.filterwarnings('ignore')` 會吞掉 `warnings.warn`，所以 warning 必須用 logger/stderr 或移除 blanket suppression。
- 最小改法：新增小型 argparse（保留既有 input/output defaults）與純 `select_historical_classes(df, allow_unversioned_v1)`；將 `analyze_locus` 寫出的欄位改成來源中立 `selection_class`，plots/docs 由 metadata 產生 title。
- Test seam：不碰 BAM/region files；直接測 selector fixture 的 Strong/Noise sample counts、unknown exclusion、warning metadata與未授權 fallback error。

### 3.4 `scripts/analysis/compare_subclone_validation.py`

**證據位置**：CLI/default paths `25-65`；input `148-177`；old support `188-195`；metrics `198-234`；rule definitions `238-258`；outputs `260-270`。

- 舊假設：current `VerificationClass in {Subclone,Strong}` 或 `HPMergedDelta<=0.10` 合併成 `MethylSupport_SubcloneLike`，並把 `VerificationClass==Subclone` 稱為 `SubcloneClass`。
- Canonical 行為：依鎖定 contract，`LegacyVerificationSupport = legacy in {Strong,Subclone}`；`LegacyClusterFirstOnly = legacy==Subclone`。不得加上 HPMergedDelta OR；若仍要分析低 delta，另立獨立欄與 rule。
- Fallback/unknown：缺 legacy欄 fail；不可用 v2 `ClusterFirstOnly` 偷換歷史母體。未知 legacy值列 count 後 fail/exclude，不能被 bool false 吞掉。
- 最小改法：讀檔後先跑 shared `L4` validator；集中 derive 兩個 boolean；將 metrics group、rules與 sample output欄名全部含 `Legacy`，移除 cellular clone wording。
- Test seam：四類 fixture assert support=`T,T,F,F`、cluster-first-only=`F,T,F,F`；missing/unknown tests；再 assert所有 output header/rule name包含 `Legacy`。
- CLI/路徑風險：未找到 repo caller；目前 default input header是 unversioned且沒有 legacy欄，所以 no-arg會 intentional fail。default output仍指向 `/big8_disk/.../InterSubMod_runs/output`，若要執行新 run 應改為 required output或 big7 synthesis；不可新增寫入舊 output 的行為。

## 4. 62 支 `VerificationClass` executable inventory

> `scope`：`LIVE` 代表可對任意/新 summary執行；`DERIVED` 代表吃 manifest/training/joined table，若上游重產就會接觸 v2；`HIST` 代表固定或明顯針對舊 cohort。`?` 表示用途不足以安全決定 cohort，因此建議先加 explicit mode/schema error，不默選。

| # | Consumer evidence | scope | Old assumption / current use | Canonical source | Fallback / unknown | Block |
|---:|---|---|---|---|---|---|
| 1 | `scripts/analysis/20260423_B7_loh_noise_signal.py:4` | HIST | `VerificationClass`只在docstring解釋 LOH_Noise；runtime實際讀 LOH_Subtype | 無 direct VC migration；LOH路徑用 `LOH-L` | `H1`; unknown計數 | H |
| 2 | `scripts/analysis/analyze_candidate_rules.py:216` | DERIVED | raw Counter，未綁 schema | `C2` | unversioned raw + metadata；unknown bucket | B1 |
| 3 | `scripts/analysis/analyze_clairs_borderline_fn.py:249-289` | LIVE | Strong/Subclone作 rescue support | `E.ClusterFirstSupport`; legacy replication用 `L4` | mode必須顯式；不可 current猜legacy | B0 |
| 4 | `scripts/analysis/analyze_cross_sample_ism_af_gradient.py:176-181` | DERIVED | 只算 plain Noise/Strong 百分比 | `L4`（historical metric）或 explicit `C2` mode | 缺欄/unknown不得變 0% | B0 |
| 5 | `scripts/analysis/analyze_gq_methylation_rescue_matrix.py:329-356,439-442` | DERIVED | 缺 class填 NA，再用 Strong/Subclone support | `E.ClusterFirstSupport`; replication `L4` | NA/unknown separate count | B0 |
| 6 | `scripts/analysis/analyze_kism_vs_cn_perread.py:205` | LIVE | 只複製 current字串 | `P`（version/current/legacy/evidence） | 缺 schema標 unversioned，不造空字串 | B2 |
| 7 | `scripts/analysis/analyze_longphase_rescue_with_methylation.py:49-61` | LIVE | 多條 Strong/Subclone rescue rule | `E.ClusterFirstSupport`; old rule用 `L4` | explicit cohort mode | B0 |
| 8 | `scripts/analysis/analyze_methyl_cluster_allele_cooccurrence.py:44,174,225` | HIST | Strong vs Noise讀 current | `L4` | `H1` only；unknown excluded+count | B0 |
| 9 | `scripts/analysis/analyze_methylation_rescue_feature_space.py:245-266,449-452` | DERIVED | Strong/Subclone support | `E.ClusterFirstSupport`; replication `L4` | NA/unknown separate | B0 |
| 10 | `scripts/analysis/analyze_to_feature_recall.py:154-160,335-339` | DERIVED | pass-through current class；另有舊 class_shift vocabulary | `P`; class_shift另立 schema | schema/version一起輸出；unknown不填空 | B2 |
| 11 | `scripts/analysis/analyze_to_tp_fp_characterization.py:437-443` | LIVE | 只認 plain Noise/Strong/Subclone | `L4`或 explicit current distribution | 缺值不可回 0%；mode metadata | B0 |
| 12 | `scripts/analysis/analyze_to_verification_scheme_adjustments.py:145,188-192` | HIST? | 從 current Strong推證據狀態 | `E`；v1 A/B用 `VerificationClass_V1_Deprecated` | 無 evidence時只允許 `H1` | B0 |
| 13 | `scripts/analysis/ari_typology.py:97-100` | LIVE | dynamic crosstab，無 schema | `C2` | unknown bucket / unversioned raw | B1 |
| 14 | `scripts/analysis/build_allele_deep_dive.py:496-499` | LIVE | **由 VerificationClass捏造 truth_label** | 必須真實 truth欄；VC只作 `C2` annotation | 缺 truth直接 fail，絕不 map class | B0 |
| 15 | `scripts/analysis/build_cross_sample_methylation_observation_workspace.py:167,495-503,715-739` | LIVE | 文件/plot固定四態且丟掉其他值 | historical panel=`L4`; current panel=`C2` | explicit panel；unknown保留 | B0 |
| 16 | `scripts/analysis/build_feature_direction_map.py:489-492` | LIVE | **由 VerificationClass捏造 truth_label** | 必須真實 truth欄；VC只作 annotation | 缺 truth fail | B0 |
| 17 | `scripts/analysis/build_loh_expanded_observation.py:59,514-516` | HIST | VC_ORDER固定四態 | `L4` | `H1`；unknown列出 | B1 |
| 18 | `scripts/analysis/build_loh_quadrant_explanation.py:175-200` | DERIVED | dynamic current plot，無版本 | `C2` | unknown bucket；禁止混版本 | B1 |
| 19 | `scripts/analysis/build_loh_round1_cross_sample_audit.py:669-670` | DERIVED | current/LOH alias只改名 pass-through | `C2` + `LOH-L` | version/unknown/provenance欄 | B2 |
| 20 | `scripts/analysis/build_loh_round1_cross_sample_audit_v2_figures.py:315` | DERIVED | derived VC圖，未帶 schema | 上游 `C2` metadata | input metadata缺失標 unversioned | B2 |
| 21 | `scripts/analysis/build_multilayer_hp_before_after_comparison.py:93-103,395-410` | HIST | before/after只容四態 | historical v1明確用 deprecated/raw view；新比較用 `C2` | 禁止 v1/v2直接同矩陣 | H |
| 22 | `scripts/analysis/build_non_loh_discrimination.py:404-406` | LIVE | **VC map成 TP/FP truth** | 必須真實 truth欄 | 缺 truth fail | B0 |
| 23 | `scripts/analysis/build_observation_O01_master_distribution.py:41-42` | DERIVED | pass-through current + LOH alias | `P` + `LOH-L` | provenance欄；unknown保留 | B2 |
| 24 | `scripts/analysis/build_observation_O03_loh_features_post_fix.py:615-649` | DERIVED | dynamic current圖 + LOH alias | `C2` + `LOH-L` | unknown bucket | B1 |
| 25 | `scripts/analysis/build_observation_O06_verification_class.py:43,82-87` | DERIVED | categorical固定四態；v2全變 NaN | historical report=`L4`; current另開 `C2` mode | 缺 legacy fail；unknown不歸 Noise | B0 |
| 26 | `scripts/analysis/build_observation_O10_read_level_methyl.py:331-346` | DERIVED | read-level plot固定四態 | `L4`（保持歷史 cohort） | legacy provenance隨 training table；unknown計數 | B0 |
| 27 | `scripts/analysis/build_observation_O15_loh_zone_metrics_hcc1395.py:360-411` | HIST | vc_order固定四態；另讀 LOH alias | `L4` + `LOH-L` | `H1`; unknown不丟 | B1 |
| 28 | `scripts/analysis/build_observation_O15b_loh_zone_metrics_cross_sample.py:517-567,1180-1193` | HIST | 固定四態且直接算 plain Noise% | `L4` + `LOH-L` | 缺欄/unknown不可顯示 0% | B0 |
| 29 | `scripts/analysis/build_phase1_training_manifest.py:282-329,459-462` | DERIVED | training export只帶裸 current/LOH alias | `P` + `LOH-L` | schema欄必帶；historical標 H1 | B2 |
| 30 | `scripts/analysis/build_phase1a_head_to_head_baseline_table.py:85-100` | DERIVED | dynamic count後只取 Noise/Weak | historical metric=`L4`; current table=`C2` | 明確 view；缺 class不可當 0 | B0 |
| 31 | `scripts/analysis/build_phase1a_split_manifest.py:40,95` | DERIVED | pass-through current | `P` | schema version/current/legacy/evidence一起帶 | B2 |
| 32 | `scripts/analysis/build_phase2_annotation_layer.py:136,324-327` | DERIVED | Strong/Subclone flag | `E.ClusterFirstSupport` | missing/NA separate；不從 current猜 | B0 |
| 33 | `scripts/analysis/build_post_hp_fix_to_loh_investigation.py:247-267` | DERIVED | dynamic current group，無 schema | `C2` | same-version guard + unknown bucket | B1 |
| 34 | `scripts/analysis/build_s5_downstream_impact_and_cv.py:66,133-176` | DERIVED | 比 class concordance但不檢版本 | `C2`; paired sides須同 schema | mixed-version fail；unknown保留 | B0 |
| 35 | `scripts/analysis/build_snv_methylation_association_analysis.py:364` | DERIVED | output欄帶裸 current | `P` | version/provenance一併輸出 | B2 |
| 36 | `scripts/analysis/build_to_feature_study_Q1Q5.py:368,659-661` | HIST | VC/LOH作 strata/plot，LOH用 deprecated alias | current `C2`; LOH=`LOH-L` | historical explicit H1 | B1 |
| 37 | `scripts/analysis/build_to_fp_provenance_analysis.py:422,599-602` | LIVE | rule硬編碼 current Strong | 若語意是雙向則 `C2==Strong_Bidirectional`; 若 cluster support則 `E` | 用途需 explicit rule schema，不能猜 | B0 |
| 38 | `scripts/analysis/compare_phase1a_model_errors.py:36-48,153-168` | DERIVED | VC作 pairing/group key，未檢版本 | `C2` + schema version入 key/metadata | mixed-version fail | B1 |
| 39 | `scripts/analysis/compare_subclone_validation.py:191-194,242` | HIST | Strong/Subclone與Subclone作 methyl/subclone support | `L4` exact contract | 缺/unknown legacy fail | B0 |
| 40 | `scripts/analysis/evaluate_rescue_with_methylation.py:119-145,178-180` | LIVE | Strong/Subclone支援 rule | `E.ClusterFirstSupport`; legacy A/B=`L4` | explicit mode；unknown separate | B0 |
| 41 | `scripts/analysis/export_phase1_manifest_shard.py:37,148,225` | DERIVED | pass-through current | `P` | schema欄不得遺失 | B2 |
| 42 | `scripts/analysis/export_phase1_read_training_table.py:130,388` | DERIVED | pass-through current | `P` | schema/version/current/legacy/evidence帶到 read rows | B2 |
| 43 | `scripts/analysis/extract_label_first_metrics.py:43,112` | LIVE | label-first export只帶 final current class | `P`，尤其 `LabelFirstSupport/EvidencePath` | v2缺 evidence fail；H1另標 | B1 |
| 44 | `scripts/analysis/fn_verdict_readback_audit.py:71,120-123` | HIST | 以 current Noise/Weak作四態 cohort | `L4` | `H1`; unknown excluded+count | B0 |
| 45 | `scripts/analysis/fn_verdict_reclassify_t1.py:13-18,83-107` | HIST | 從 current四態反推 label/cluster booleans | `E.LabelFirstSupport/ClusterFirstSupport` | 無 evidence只可H1；NA不可 false | B0 |
| 46 | `scripts/analysis/generate_pi_report_figures_self_phasing.py:296` | HIST | 只有圖中文字引用 VC | historical-only wording；不作 selector | 不讀新 canonical | H |
| 47 | `scripts/analysis/pon_cross_sample_and_h2009_diagnosis.py:315-322` | DERIVED | dynamic distribution | `C2` | unknown bucket / version metadata | B1 |
| 48 | `scripts/analysis/run_phase1a_read_classifier_benchmark.py:66,165,322` | DERIVED | training/eval pass-through current | `P` | schema欄進 feature provenance；禁止混版本 | B2 |
| 49 | `scripts/analysis/run_phase2_paired_model_feature_analysis.py:205-220,423-427` | DERIVED | Strong/Subclone support | `E.ClusterFirstSupport` | missing/unknown separate | B0 |
| 50 | `scripts/analysis/run_to_support_feature_diagnostics.py:205-260` | DERIVED | pass-through current到 diagnostics | `P` | schema/evidence欄一起輸出 | B2 |
| 51 | `scripts/analysis/seqc2_cnv_cross_sample_and_rootcause.py:515-529` | DERIVED | dynamic VC distribution，未綁 schema | `C2` | unknown bucket / version metadata | B1 |
| 52 | `scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:23,136,205` | HIST | 標成 current 9-state；仍說 Subclone_ID | `C2` + `L4` + `R1` | input JSON須有 schema/status；unknown呈現 | B0 |
| 53 | `scripts/analysis/subcluster_split_allsamples/hcc1395_full_accounting.py:91-105,171-172,229-233` | HIST | current「9state」、手寫 clean set、舊 Subclone_ID | `C2`; Significant直接讀；`R1`; `LOH-L` | schema/status缺失 fail，不能手算 clean | B0 |
| 54 | `scripts/analysis/summarize_phase1a_prediction_failures.py:31,66` | DERIVED | current class作 grouping key | `P/C2` + version進 group | mixed-version fail；unknown bucket | B1 |
| 55 | `scripts/analysis/tp_fp_structure_label_association.py:114,345-356` | LIVE | 固定四態 CMH/matching；舊 artifact註解 | historical panel=`L4`; exact within-HP用 canonical欄，不讀 stub | explicit mode；unknown保留 | B0 |
| 56 | `scripts/analysis/tp_fp_structure_label_figures.py:59-65` | DERIVED | 固定四態畫圖 | 跟上游明確 L4或C2 metadata | missing metadata fail | B1 |
| 57 | `scripts/analysis/validate_method_design.py:286` | LIVE | 缺 current class默認 Noise | 視目的用 `C2`或`E`; unknown不可Noise | missing/unknown fail or explicit bucket | B0 |
| 58 | `scripts/analysis/verify_class_decision_tree_audit.py:124-176` | LIVE | 從四態推cluster_sig、not-noise；missing=Noise | `E` + `C2` 11類 | missing/unknown fail，不能Noise | B0 |
| 59 | `scripts/regression/regression_check.sh:28` | TEST | golden只抽裸 current class | 加 schema/version/v1 alias/evidence；更新 golden需人工審核 | schema drift應 fail，不自動重建 | B1 |
| 60 | `scripts/validation/metrics/compute_metrics.py:45-69` | LIVE | dynamic Counter；缺欄用 Unknown | `C2` + schema version | `UnknownCurrentClass`與 missing分開 | B1 |
| 61 | `tools/compare_vcf_results.py:143-166,244-272` | LIVE | plot固定四態 | `C2`; optional legacy view | unversioned raw +標籤；unknown bucket | B0 |
| 62 | `tools/find_verification_candidates.py:31-37,117-123` | LIVE | current欄找四態 | `L4` | 缺/unknown legacy fail | B0 |

## 5. `LOH_Subtype` provenance inventory（23 支）

這一組在 v2 因 deprecated exact alias 暫時不一定改變數值，但若繼續讀舊欄，報告無法表達「此 subtype 是由 legacy VerificationClass 衍生」。共用策略應是：新 canonical一律 `LOH_Subtype_LegacyVC`；若同時有 `LOH_Subtype` 則逐列 assert equality；只有 explicit H1可讀舊欄。

| Consumer evidence | 主要用途 | 建議 | Blocking |
|---|---|---|---|
| `scripts/analysis/20260423_B3_paired_obs18.py:132` | historical feature select | `LOH-L`; H1 warning | H |
| `scripts/analysis/20260423_B4_S4_discrimination.py:74-76` | `None` cohort filter | `LOH-L`; unknown不得當 None | B1 |
| `scripts/analysis/20260423_B5_colo829_s1_fold.py:79,126-138` | LOH_Strong filter/distribution | `LOH-L`; H1 | H |
| `scripts/analysis/20260423_B7_loh_noise_signal.py:92-139,198-212` | LOH_Noise/Strong comparison | `LOH-L`; H1 | H |
| `scripts/analysis/20260423_phase3_synthesis.py:82` | fillna None後分組 | `LOH-L`; missing與None分開 | B1 |
| `scripts/analysis/20260423_s5_loh_af_cn_scatter.py:73` | feature select/pass-through | `P` + `LOH-L` | B2 |
| `scripts/analysis/20260424_X5_crosssample_obs18_S3S5.py:29` | feature select | `LOH-L`; H1 | H |
| `scripts/analysis/20260424_X5v2_corrected_S3S5.py:56-67,153-156` | S1-S4 rules | `LOH-L`; H1 | B1 |
| `scripts/analysis/20260424_X6_merge_caller_af_S3S5.py:65,95` | merge/pass-through | `P` + `LOH-L` | B2 |
| `scripts/analysis/build_allelesig_loh_study.py:138,464` | subtype cohort filter | `LOH-L`; unknown separate | B1 |
| `scripts/analysis/build_loh_feature_validation.py:42,114,683-688` | fixed subtype validation | `LOH-L`; include UnknownLegacyLOH count | B1 |
| `scripts/analysis/build_loh_round1_cross_sample_audit.py:670` | rename/pass-through | `LOH-L` + provenance | B2 |
| `scripts/analysis/build_observation_O01_master_distribution.py:42` | master pass-through | `P` + `LOH-L` | B2 |
| `scripts/analysis/build_observation_O03_loh_features_post_fix.py:235-237,376-386` | crosstab/paired-TO | `LOH-L`; mixed provenance fail | B1 |
| `scripts/analysis/build_observation_O06_verification_class.py:55,85-87` | fixed LOH categories | `LOH-L`; unknown不變 NaN | B1 |
| `scripts/analysis/build_observation_O15_loh_zone_metrics_hcc1395.py:573-595` | LOH heatmap | `LOH-L`; H1 | H |
| `scripts/analysis/build_observation_O15b_loh_zone_metrics_cross_sample.py:118` | feature select | `LOH-L`; H1 | H |
| `scripts/analysis/build_phase1_training_manifest.py:285,329,462` | training export | `P` + `LOH-L` | B2 |
| `scripts/analysis/build_to_feature_study_Q1Q5.py:118-120,659-661` | subtype strata | `LOH-L`; unknown separate | B1 |
| `scripts/analysis/build_to_feature_study_Q6.py:96-97` | subtype strata | `LOH-L`; unknown separate | B1 |
| `scripts/analysis/run_phase1a_round3_loh_feature_test.py:53,118-119` | ML categorical feature | `LOH-L`; missing不可填 None | B0 |
| `scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:25,270` | derived HTML | 上游顯式 `LOH-L`，標 legacy provenance | B1 |
| `scripts/analysis/subcluster_split_allsamples/hcc1395_full_accounting.py:105` | distribution | `LOH-L`; schema/status guard | B1 |

## 6. Region-stratum / stale artifact consumers

| Evidence | Old behavior | Canonical behavior | Policy |
|---|---|---|---|
| `scripts/analysis/subcluster_split_allsamples/hcc1395_full_accounting.py:229-233` | 把 `Subclone_ID` 放在 normal-anchored欄並用「無 subclone」語境解釋 | 讀 `RegionStratum_ID/Label/Reason` 與 `region_stratification_status.tsv`; `-1`依 reason解讀 | v2缺 canonical欄 fail；deprecated alias只做 equality check |
| `scripts/analysis/subcluster_split_allsamples/build_hcc1395_html.py:205` | HTML顯示 `Subclone_ID=-1` 的舊敘述 | 改寫為 `RegionStratum=Unassigned` + reason/status；不得稱 cellular clone | derived JSON缺 schema/status fail |
| `scripts/analysis/tp_fp_structure_label_association.py:113-116` | 註解稱 exact within-HP state在 `subclone_structure.txt` | 舊檔將成 deprecation stub；只可讀 canonical per-row/evidence欄或 `region_stratification.tsv` | 不要新增 parser讀 stub |

未找到 executable 直接讀 `n_subclones`、`n_occupied_region_strata` 或 parse `subclone_structure.txt`。這代表 class counter/n_occupied修復主要由 C++ unit/integration tests承接，不需要為不存在的 Python caller造 compatibility layer。

同樣地，`region_stratification_status`、`region_stratification.tsv`、`RegionStratificationSchemaVersion`、`RegionStratum` 在 `tools/scripts/tests` 的 executable命中為 **0**：目前沒有任何 status-aware consumer，也沒有把 main CSV與 canonical region-stratification artifact join後再檢 status的 reader。P0-B落地時需要新增這個契約，不能假設現有 consumer已會處理 `VALID/INSUFFICIENT_REGIONS/NOT_APPLICABLE_TUMOR_ONLY/FAILED`。

### 6.1 Positional / exact-header compatibility audit

在 scoped consumers 內，沒有找到以 numeric CSV column index直接定位 `VerificationClass`、`LOH_Subtype`、`Subclone_ID` 或 `n_subclones` 的 reader；現有 `iloc/iat` 命中都是對已建 pivot matrix取 cell，不是讀原始 schema位置。也沒有 downstream 把 `Subclone_ID` 當 vector index或用 `max_id+1` 配置容器：兩個命中都只是 HCC1395 report文字/空值盤點。sparse-ID與 occupied-count的 index safety因此是 producer/C++ test責任，不是 Python compatibility理由。

但下列是 **exact header projection**，若不更新會在 read/export時靜默丟掉 append-only v2欄：

- `scripts/analysis/extract_label_first_metrics.py:9-46,72`：固定 `OUTPUT_COLUMNS`；
- `scripts/analysis/run_to_support_feature_diagnostics.py:170`、`analyze_to_verification_scheme_adjustments.py:44`、`build_to_fp_provenance_analysis.py:669,713`：固定 `DictWriter.fieldnames`；
- `scripts/analysis/build_snv_methylation_association_analysis.py:356`、`analyze_gq_methylation_rescue_matrix.py:269-287`、`build_cross_sample_methylation_observation_workspace.py:395-396`、`build_observation_O06_verification_class.py:80`：name-based `usecols`/projection；
- phase1 manifest/read-table/shard scripts：固定 output field lists，已在 inventory標為 `P/B2`。

這些不是 column-order break，但仍要透過 `schema_provenance_columns()`集中 append欄，否則上游live CSV正確、derived table仍會失去 provenance。

## 7. 可由 shared helper 一次解決的群組

建議先建立單一 Python contract helper（例如 `tools/verification_schema_contract.py`；實際落點由主實作者依 import topology 決定），不要在 62 檔各自複製 enum與 fallback：

1. `validate_current_v2(df)`
   - 驗 `VerificationSchemaVersion` 單一且為 2；回傳 11-class ordered categorical + `UnknownCurrentClass` diagnostics。
   - 提供完整零頻 reindex，直接解 `compare_vcf_results`、O06/O10/O15、cross-sample plots與 metrics counters。
2. `select_legacy4(df, allow_unversioned_v1=False)`
   - canonical只讀 `VerificationClass_Legacy`；H1才讀裸 `VerificationClass`。
   - 回傳 selection field/version/warning/unknown counts，直接解 candidates、methyl replication與所有固定四態報告。
3. `read_evidence(df)`
   - typed parse `LabelFirstSupport/ClusterFirstSupport/WithinHPSupport/DispersionWarning/EvidencePath`；保留 offline NA。
   - 直接解所有 `Strong/Subclone` support/rescue scripts與兩支 verdict audit，避免 class反推 evidence。
4. `select_loh_legacy(df, allow_unversioned_v1=False)`
   - 優先 `LOH_Subtype_LegacyVC`；assert deprecated alias exact equality；unknown/missing分開。
   - 可一次遷移 23 支 LOH consumers。
5. `select_region_stratum(df, status_path)`
   - 驗 schema v1、status、sentinel/reason、deprecated ID equality；只回傳 region-stratum語彙。
6. `schema_provenance_columns()` / `copy_schema_provenance(row)`
   - 給 training/export/manifest scripts；避免 derived tables只留下裸 current class，讓下游再度無法判斷 cohort。

Helper 本身必須有 table-driven tests；consumer只測「是否呼叫正確 view」與 output labels。這能把 62 檔散改縮成約六類 integration，而不是 62 套 enum/fallback。

## 8. 建議實作批次與 acceptance seams

### Batch A — P1-B release blocker（先做）

1. 四支規格鎖定 consumer。
2. 所有 `B0` evidence/truth/cohort selectors。
3. shared helper + unit tests。
4. `scripts/regression/regression_check.sh` 的 schema-aware golden extraction。

Acceptance：

- v2 11 classes + 1 unknown表/legend完整；零頻類不消失。
- legacy四類fixture精確各 1；missing/unknown fail closed。
- evidence support直接取 boolean，不從 current/legacy class推導。
- 三支 `truth_label = VerificationClass.map(...)` 路徑被移除或在缺真 truth時 fail。

### Batch B — reports / historical replication

- O06/O10/O15、cross-sample workspace、multilayer/TP-FP figures等固定四態報告。
- 每支標 `historical replication`、selection field與H1 warning；若要 current view，另開明確 mode而不是自動折疊。

### Batch C — provenance/pass-through

- training manifest、read table、shards、diagnostics、joined tables全部保留 schema version/current/legacy/evidence。
- LOH consumers改用 `LOH_Subtype_LegacyVC`，alias只做 equality check。

## 9. CLI 相容性總結

| Consumer | Existing caller/default | Safe compatibility choice | Intentional break |
|---|---|---|---|
| `compare_vcf_results.py` | `run_batch_vcf_analysis.sh:234-237` 傳既有三組 args | 新 view flag保持 optional，default current/schema-aware | 不應有 |
| `find_verification_candidates.py` | repo未找到 caller | 不必新增 required CLI；default report改 legacy | 無 legacy欄即 fail |
| `analyze_methyl_cluster_allele_cooccurrence.py` | 無 caller；hardcoded unversioned old input | 保留 path defaults，新增 optional `--allow-unversioned-v1` | 現有 no-arg因缺 legacy fail |
| `compare_subclone_validation.py` | 無 caller；default old big8 input/output | input可保留但會被validator擋；新 output不可再落舊 big8 root | no-arg缺 legacy fail；output names改含 Legacy |

## 10. 稽核命令、實際輸出與修改證明

### 10.1 Targeted inventory

輸入：`InterSubMod/tools`、`InterSubMod/scripts`、`InterSubMod/tests`。

```bash
rg -l 'VerificationClass|VerificationClass_Legacy|Subclone_ID|RegionStratum|n_subclones|n_occupied_region_strata|subclone_structure\.txt|LOH_Subtype' \
  tools scripts tests \
  --glob '!**/__pycache__/**' --glob '!**/*.pyc' \
  --glob '!**/*.md' --glob '!**/*.html' --glob '!**/*.json' \
  --glob '!**/*.csv' --glob '!**/*.tsv' --glob '!**/*.txt'
```

實際摘要：`VerificationClass=62 token-hit executable-code files`、`VerificationClass_Legacy=2`、`LOH_Subtype=23`、`Subclone_ID/RegionStratum=2`、`n_subclones/n_occupied=0`、`subclone_structure.txt=1 comment-only`、`tests=0`。62中已在逐列inventory區分runtime selector、pass-through與文字-only，不把docstring引用當科學cohort consumer。

計數口徑對照：若把 producer implementations `src/core/RegionProcessor.cpp` 與 `src/core/SubcloneAnalyzer.cpp` 一併納入，`VerificationClass` token-hit code files是 **64**；本表依 delegated P1-B consumer範圍只列 `tools/scripts/tests` 的62支，避免把producer重複算成downstream consumer。

### 10.2 固定 input header readback

```bash
head -n 1 /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv
head -n 1 /big8_disk/liaoyoyo2001/InterSubMod_runs/output/advanced_analysis_20260119/advanced_samples.csv
```

實際結果：兩者都有裸 `VerificationClass`；前者另有裸 `LOH_Subtype`；兩者都沒有 `VerificationSchemaVersion` 或 `VerificationClass_Legacy`。這支持兩個 historical fallback必須由 explicit H1授權，而不能自動猜。

### 10.3 四支 consumer integrity

```text
706b40ea7a431e1e882b6a4d93e7e76945b7d02789565c58bf3149fa1111c0ac  tools/compare_vcf_results.py
c941bbeac6056d31303f31163103b83d07e37d06bdbb3c35942dc868aef4f443  tools/find_verification_candidates.py
c4c53225a7c809cf40f223f3273f559f71a2db1046503620eed48f30cda7e165  scripts/analysis/analyze_methyl_cluster_allele_cooccurrence.py
736c54fd381a6140a66f0c0355b1724fd43488a68409ce910a70becf415a30d2  scripts/analysis/compare_subclone_validation.py
```

在寫本 audit 前，四支檔案的 targeted `git status --short` 與 `git diff --` 皆無輸出。本文件是唯一授權新增物；沒有修改 consumer。

## 11. 限制與主 agent 決策點

1. 此 inventory是 static source audit；沒有執行62支 scripts，也沒有聲稱其歷史數值正確。
2. `scope=HIST/DERIVED/LIVE` 由 path、CLI與用途判讀；標 `?` 的 scientific cohort必須由實作者用 explicit mode保留，不可默選。
3. 某些 derived tables可能已有自己的 `verification_class` 小寫欄；若其 upstream未同步帶 schema metadata，consumer即使改名仍無法判斷 v1/v2。
4. 規格只明文鎖定四支 consumer 的 exact行為；其他62支中的 scientific cohort若要從 legacy切到 current，屬研究定義變更，應另立決策，不在 hotfix中默做。
5. Deprecated `LOH_Subtype` 在一個 migration release仍與新欄逐列相等，所以 LOH migration主要是 provenance correctness，不代表立即數值改變。

**Release判定：只修四支已知 consumer仍不足以勾選「完整 P1-B」。** 最低可接受做法是：四支 + shared helper + 所有 B0 paths遷移；其餘B1/B2至少加 schema guard/metadata或明確列入 follow-up，不得以搜尋未見而宣稱不存在。

## 12. 2026-07-14 migration closeout

- 本 inventory 的 62 個 VerificationClass paths 與 23 個 LOH paths，合併後為 **74 個唯一 executable consumers**（11 個同時帶兩種契約）。
- 最終狀態：**74 MIGRATED / 0 unresolved**。
- Production v2 consumer 預設要求明確 schema/provenance；歷史 H1 只能透過 explicit flag 授權並保留 warning/receipt。
- 逐 consumer 的 contract、blocking、disposition、test/guard 與 notes 位於 `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/consumer_migration_status.tsv`。
- `generate_pi_report_figures_self_phasing.py` 僅有 historical wording，沒有 selector；已明標 `Historical legacy VerificationClass (v1)`，未虛構新 taxonomy 行為。
- 完成範圍是本文件鎖定的 inventory，沒有用一次新的全 repo 大掃描來宣稱未知檔案不存在。

Release gate：P1-B consumer migration **PASS**；P2 frozen migration 仍須獨立通過 atomic publication 與 golden invariance gate。
