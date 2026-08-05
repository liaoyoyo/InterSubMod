<!--
建立時間: 2026-07-17T14:22:00+08:00
更新時間: 2026-07-22T10:54:15+08:00（加入fresh singleton audit、supplemental report與signed receipt固定路徑；20項稽核題目不變）
目標: 在最終數據產生前預先固定外部 Claude Code 的唯讀結果稽核題目，避免事後挑選驗證項目
處理範圍: 全 sSNV M1/M2/G1/G2/R1/B1/C1/L1/L2、provenance、個案與 portable HTML
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_m2v5_recovered_completion_chain.sh
-->

# External Claude Code final audit prompt

You are an external, read-only scientific-computing reviewer. Audit the completed InterSubMod all-sSNV focal-ALT methyl multigroup and same-read sSNV cooccurrence release. Do not edit, create, delete, rename, chmod, or overwrite any file. Do not run the large producers. Small read-only shell/Python recalculations that write nothing are allowed.

## Scope

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Read these contracts and reports first:

1. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md`
2. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/00_INDEX.md`
3. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/implementation-notes.md`
4. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260717_生物學claim_ceiling與clone_lineage獨立稽核_01.md`
5. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260717_tumor_REF來源attestation與release接縫獨立稽核_01.md`

Audit these completed outputs. Report a release blocker when a required/applicable output is absent or lacks its required `pass=true` receipt. If G2 has zero selected candidates, accept a schema-valid `NOT_APPLICABLE_ZERO_SELECTED_CANDIDATES` receipt and the contractually absent downstream analysis table/directory; independently verify that zero-candidate branch instead of demanding fabricated positive outputs.

- Cooccurrence: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity/`
- Tumor-REF: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/`
- Strict: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/strict_methyl_candidate_confirmation_v3_m2v5_source_authority_v5/`
- Matched normal run/analysis: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/matched_normal_candidate_controls_v3_m2v5_source_authority_v5/` and, when applicable, `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5/`; otherwise verify the native zero-candidate N/A receipt.
- CN/CCF: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5/`
- Final machine dataset: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_final_report_dataset_v5_m2v5_source_attested/`
- Final report: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_final_report_v5_m2v5_source_attested/`
- Official source receipt: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tumor_ref_recovery_source_identity_v1/post_run_source_identity.receipt.json`
- Signed source authority: `docs/provenance/source_authorities/20260722_all_ssnv_focal_alt_release_source_authority.v7.json`
- Signed final dataset release receipt: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/task_b_final_dataset_release_receipt.v1.json`
- Signed final report release receipt: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/task_b_final_report_release_receipt.v1.json`
- Positional-singleton audit: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v2_source_authority_v7/`
- Positional-singleton supplemental report: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_report_v5_source_authority_v7/`
- Signed positional-singleton supplemental release receipt: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/positional_singleton_supplemental_release_receipt.v1.json`

## Required independent checks

Run and report each check separately as PASS / FAIL / NOT VERIFIABLE, with exact evidence paths and recomputed values:

1. Verify the frozen population is exactly 7 dataset rows / 6 biological samples / chr1-22 / 469,849 LongPhase-S recalibrated FILTER=PASS biallelic dataset-sites, with truth split 335,296 TP + 7,745 FP + 126,808 UNASSESSED.
2. Verify the latest same-run HP/PS read-tag join and tree-input receipts, and confirm that no regular LongPhase-S tagged BAM was persisted or used as the current tag source.
3. Recompute M1 FLAGGED count and rate over all 469,849 dataset-sites; do not reinterpret it as biological prevalence.
4. Recompute M2 eligible/evaluable/NOT_EVALUABLE reason conservation, including the unique K=11 case and aligned-below-negative-power exception count.
5. Recompute the complete G1 exact-testable pair family, global BH/BY discoveries, effect/callability/conditional gates, and formal pair count. Confirm all M2 exact-testable pairs enter the global family before FDR.
6. Recompute G2 using only same-complete-read effect-supported top markers, at least two positions separated by >=20 bp, exact-BY formal pairs, and the global-BY joint conditional test. Confirm 20 bp is not described as statistical independence.
7. Recompute R1 strict robustness status and all NOT_EVALUABLE/failure reasons; confirm it is a post-selection robustness audit, not a second FDR confirmation.
8. Recompute tumor-REF and matched-normal B1 controls. Confirm the lenient background predicate lacks the M1 membership-ARI requirement and is only a superset on the same background payload.
9. Verify B1 uses exactly one outcome-blind prespecified G1 pair per site, with deterministic depth/distance/identity tie-break, and never selects a pair based on favorable four-state outcome.
10. Recompute four-state RR/AR/RA/AA totals and relation decisions. Confirm zero focal-REF depth and fixed-error ceiling are NOT_EVALUABLE, not evidence for or against ancestry.
11. Recompute C1 CN/CCF coverage and statuses from authority-locked exact-locus annotations. Confirm any missing/unreviewed model remains NOT_EVALUABLE.
12. Verify L1 and L2 cannot PASS without orthogonal same-cell identity and identifiable >=3-state order evidence. Confirm bulk long reads do not pair HP1 and HP2 back to the same cell.
13. Recompute candidate counts at four grains: dataset-site, unique biological site, unique regional component, and unique ALT readset. Flag any report text that treats these as independent replicates.
14. Recompute per-dataset and TP/FP/UNASSESSED tables and identify Simpson/pseudoreplication risks, especially HCC1395 versus HCC1395_DORADO biological n=1.
15. Trace up to three final candidates (all of them if fewer than three; verify the structured zero-candidate branch if none), one no-candidate extreme M1 case, the K=11 case, and one fully explained negative case from source assignment/pair rows through the final catalog and visible report.
16. For each traced candidate, verify read identity, focal ALT restriction, partner R/A/O/X calls, latest HP/PS, methyl group, pair table, complete-read signature, four-state table, tumor-REF, normal, and CN fields where applicable.
17. Verify every final receipt's declared input/output SHA-256, row counts, key-set reconciliation, source identities, chronology, and the post-downstream frozen-input audit.
18. Verify source-attestation schema 1.2 independently checks the hash-bound snapshot/manifest command, exact repo-relative analyzer token, trusted v2 verifier identity, and official receipt chronology.
19. Compare `final_report_dataset.json`, `candidate_catalog.tsv`, `candidate_witness_pairs.tsv`, `report.md`, `artifact.json`, and rendered HTML for numerical and semantic consistency.
20. Inspect desktop/mobile QA JSON and screenshots. Flag clipped/overlapping text, missing charts/tables, misleading titles/denominators, inaccessible source affordances, or any claim above the evidence ceiling.

## Required output

Use findings-first format:

1. High / Medium / Low findings with exact file and line or record keys.
2. A 20-row verification matrix for the checks above.
3. Independently recomputed headline numbers and at least five traced cases where available; for a zero-candidate branch, include all available negative/extreme/edge cases and explain the reduced count.
4. Claim-ceiling verdict for M1, M2, G1, G2, R1, B1, C1, L1, L2.
5. Final verdict: `GO`, `GO_WITH_CORRECTIONS`, or `NO_GO` for sharing the report.

Do not treat a `pass=true` integration receipt as biological confirmation. The strongest permitted wording without orthogonal cellular identity is a read-level/local multi-marker molecular-haplotype candidate compatible with one or more bulk mutation-order models, not a confirmed subclone, clone count, or linear evolution tree.
