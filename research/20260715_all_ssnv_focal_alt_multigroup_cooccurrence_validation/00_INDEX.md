<!--
建立時間: 2026-07-15T08:01:43+08:00
目標: 提供全 sSNV focal-ALT 甲基多群、sSNV 共現與拓撲驗證的跨 session 單一入口
處理範圍: 7 datasets / 6 biological samples / chr1-22 / 469,849 LongPhase-S recalibrated FILTER=PASS biallelic sSNVs
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/pre-decision-audit.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/implementation-notes.md
  - InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md
-->

# 全 sSNV focal-ALT 甲基多群 x sSNV 共現驗證索引

> **狀態：IN_PROGRESS / Task Type B。** 全量 extraction、469,849-site M1/M2 screen與recovered tumor-REF已PASS；cooccurrence v5因raw-BAM RG-only duplicate identity fail-closed，formal output不存在；partial preflight v4已主動中止並排除，fresh dependency-attested v5必須全102,842 tasks通過後才允許v6。任何沒有`pass=true` final receipt的局部輸出都不可引用為完整end-to-end結果。
> **Release runner gate：** 只使用`run_m2v5_recovered_completion_chain.sh`與source-attestation receipt schema 1.2.0；已凍結的歷史`run_source_attested_release_chain.sh`仍綁v1 verifier，未處理repo-relative command token，不得用於本次release。

最新的全體/positional-singleton分母、M1/M2比例與claim ceiling：
`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/20260719_latent_molecular_substructure證據分層與比例_01.md`。
此文件明確區分positional singleton、local-partner-unavailable與read-sharing degree-zero；
在正式cooccurrence/result/report receipts完成前狀態為PRE-RELEASE。

## 研究問題與目標

- 問題：共同 focal-ALT 所涵蓋的 reads，能否由甲基資訊分成 latent groups，且這些 groups 是否與同一 read 的其他 sSNV R/A 狀態相關？
- 目標：G3（read-level epigenetic value）、G4（跨樣本與跨平台重現）、G5（可被外部重算的證據鏈）。
- 證據上限：甲基群單獨只支持 read-level epigenetic heterogeneity；通過 genetic association 後可稱 local methyl-genetic co-segregation。沒有 normal/CN/CCF/cellular truth，不得稱 confirmed subclone 或 linear evolution。

## 權威輸入

| 層 | 路徑 | 已驗證狀態 |
|---|---|---|
| Frozen manifest | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json` | PASS；469,849 = 335,296 TP + 7,745 FP + 126,808 UNASSESSED |
| LongPhase-S producer | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/` | 7/7 producer PASS；latest HP/PS sidecar |
| Layered PASS context | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/` | 7/7 canonical PASS |
| Full terminal HP/PS join | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_summary.json` | 469,849 sites；38,345,639 observations exact；multimatch=0 |
| Latest-tag regression | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/latest_tag_projected_join_audit.json` | 7,745/7,745 sites；434,759/434,759 reads exact；full identity multiplicity=1 |
| Matched-normal pre-candidate audit | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/matched_normal_methyl_tag_audit.v3_pre_candidate.json` | 7/7 readable；各1,000/1,000 primary reads具MM+ML |
| Frozen input pre-run immutability | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/frozen_input_immutability.pre_v2.json` | 77/77 artifacts PASS；7/7 samples PASS |
| Extraction reference identity | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/extraction_reference_identity_audit.v1.json` | reference+FAI full SHA；7/7 receipt binding PASS |
| Tree input contract | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/latest_tree_input_contract_audit.json` | same-run LongPhase-S PASS VCF；7/7 byte identity PASS |

## 執行版本

| 版本 | 狀態 | 路徑 / receipt |
|---|---|---|
| v1 original binary | **ABORTED / EXCLUDED** | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/v1_aborted_run_incident.json`；不得分析 |
| Verification fix smoke | PASS 2/2 | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/verification_fix_smoke.json` |
| v2 full extraction | **PASS 469,849/469,849** | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/intersubmod_all_ssnv_v2_verification_fix/` |
| v2 batch receipt | PASS 7/7 | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_intersubmod_batch.v2_verification_fix.json` |
| Exact artifact reconciliation | PASS；3類missing/extra/duplicate/empty=0 | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_output_reconciliation.v2_verification_fix.json` |
| Frozen input post-extraction audit | PASS 77/77 | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/frozen_input_immutability.post_v2.json` |
| v3 historical screen | **PARTIAL / SUPERSEDED** | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v3_terminal_tag/`；不可發布 |
| v10 recovered full screen | **PASS 469,849/469,849** | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/` |
| M1 / M2 screen | **PASS** | M1=102,842/469,849；M2=919/1,867；`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v3.json` |
| Cooccurrence attempts 1-5 | **ABORTED / EXCLUDED；formal output=0** | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_execution_incident.v2_through_attempt5.json`；v5 source另存`results/aborted_sources/` |
| Raw-BAM identity preflight v4 | **PRECHECK_ABORTED / EXCLUDED；24,000/102,842；JSON receipt=0** | partial log=`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v4_raw_identity_full_runtime.log`；不得續用 |
| Raw-BAM identity preflight v5 | **REVIEW GATE；fresh 102,842-task全量，不可subset** | 唯一允許輸出=`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v5_dependency_attested_raw_identity_full_runtime.json` |
| Cooccurrence v6 | **BLOCKED ON PREFLIGHT PASS；fresh output尚不存在** | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked/` |
| Tumor-REF recovered control | **COMPLETE；102,842/102,842；pass=true** | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/` |
| Tumor-REF source attestation v2 | **STRICT REPO-RELATIVE PRE-PROBE PASS；official receipt由current completion runner建立** | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json` |

v2 binary SHA-256：`8d082edfea79b1e75600e377d6f00f90dabdf071b739054d65fdfa4806a10eb7`。輸入 BAM 只讀；沒有產生或覆蓋 LongPhase-S tagged BAM。raw BAM 提供 allele/MM/ML，latest HP/PS 由同一次 producer sidecar exact join。Source-attestation reviewer要求的完整`command_binding`與可信v2 verifier身份已由final builder獨立重驗；relative token必須是完整repo-relative exact path，不能只用basename/suffix。

## 分析 Gates

1. v2 batch 7/7 exit=0，且 log 無 `Region ... failed`。
2. reads/methylation/BERNOULLI 的 `(sample,chrom,pos)` 集合與 frozen 469,849 keys exact equal。
3. focal methyl groups truth-blind、cooccurrence-blind建立；latest HP/PS sidecar projection必須只有一個distinct full identity。raw BAM若同projection有多筆record，只能在SAM core與除`RG`外全部typed auxiliary tags完全相同時折疊；任何其他差異hard fail並禁止release。
4. partner universe為 focal ±5 kb 內所有 frozen PASS sSNVs；`O`與`X`不得併入REF。
5. genetic association正式 discovery 使用依賴保守的 global BY，global BH 只作對照；G2另要求同一complete-read set至少2個effect-supported top markers且相距>=20 bp。
6. focal ancestry另用所有焦點REF+ALT reads的RR/AR/RA/AA；只看ALT reads不能識別祖先。
7. tumor-REF、matched-normal、HP/strand、read geometry、strict methyl null、HCC1395/DORADO replication完成後才形成最終候選。

## 預定最終產物

- Machine results：compressed site/pair/read-assignment tables、summary JSON、run manifests與checksums。
- 分析文件：完整中文 Markdown 報告與跨 session handoff。
- 視覺化：aggregate funnel、truth/sample效果、coverage、null、共現四狀態、association、dependency、topology及個案熱圖。
- HTML：standalone研究報告，通過desktop/mobile browser QA。
- 審查：多 agent consensus + 外部 Claude Code read-only重算與問題清單。

完整執行細節、偏離與即時證據見`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/implementation-notes.md`。
