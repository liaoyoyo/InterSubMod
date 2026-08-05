<!--
建立時間: 2026-07-18
目標: 獨立檢查完整 Task-B input/output/command path alignment 與 source authority
處理範圍: P0 pre-run source authorization；只讀與 bounded 重算
關聯檔案: completion runner, final dataset builder, result finalizer, M2 summary tests
-->

# External Reviewer A v11: complete Task-B chain alignment

Perform a fresh independent, read-only P0 review of the exact live source.
Two prior review rounds found fail-closed path drift before signature:

1. v5 shell output directories versus v4 Python output constants.
2. current v3/v4/v7/v8 shell inputs versus stale v2/v6 builder inputs.

The authority has never been created or signed. Both issues are now reported
fixed. Matched-normal has two deliberate terminal branches: analysis-v3 for
nonzero selected candidates, or controls-v3 carrying an explicit
NOT_APPLICABLE receipt for zero candidates. Review the complete chain rather
than only output names.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Inspect:

- `scripts/run_m2v5_recovered_completion_chain.sh`
- `scripts/build_all_ssnv_final_report_dataset.py`
- `scripts/finalize_task_b_result_release.py`
- `scripts/run_matched_normal_candidate_controls.py`
- `scripts/analyze_methyl_ssnv_cooccurrence.py`
- `scripts/m2_screen_gate.py`
- `scripts/release_source_authority.py`
- `scripts/run_cooccurrence_v6_source_locked.sh`
- `tests/test_build_all_ssnv_final_report_dataset.py`
- `tests/test_terminal_validation_chain.py`
- `tests/test_analyze_methyl_ssnv_cooccurrence.py`
- `logs/pytest_full_pre_authority_v14_input_path_alignment_canonical.xml`
- `claim-contract-v5.md`
- `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v5_summary_hotfix/ed25519_public.pem`

Independently verify:

- HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- Exactly 23 protected roles, all mode `0444`, real-newline source digest
  `20232fc5433fdc5a5210c731bfd3b2cbf8ed0c1c04abc2781196ea67b08aad12`.
- JUnit SHA-256
  `511366ce5fe27a8382820b99fd1ddbf6d199c32f1ec58f1865e637ab0bd427fd`,
  with 452 tests and zero failures, errors, or skips.
- Public-key SHA-256
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  private key remains mode `0400` before the one allowed signature.
- Runner and builder agree exactly on manifest, screen, tumor-REF,
  independent-M2 audit, tumor-REF source receipt, cooccurrence v7,
  preflight v8, strict v3, primary-pre v4, primary-post v3, and CN/CCF v3.
- `CANONICAL_MATCHED_NORMAL_DIRS` equals exactly the runner's analysis-v3 and
  controls-v3 directories. Builder path validation accepts only those two.
- Builder produces exactly two canonical direct commands differing only in
  matched-normal directory. Finalizer accepts exactly those commands and still
  independently checks the signed builder receipt.
- Runner, builder, and finalizer agree on dataset v5; runner and finalizer
  agree on report v5. Create, verify, result signature, report signature, and
  report verification bind those paths.
- The builder/finalizer live sources contain no prior canonical v2/v4-output/
  v6-input directory constants. Tests assert fresh fixed inputs, both matched
  branches, two commands, stale-name absence, and v5 outputs.
- M2 empty-map allowlist remains exactly the three true pre-axis gates;
  accepted rows receive a distinct not-evaluated census status. Dedicated
  tests reject partial axes, extra axes, and post-axis empty maps.
- K>10 remains M1 PASS, M2 NOT_EVALUABLE, G1/G2/B1 NOT_RUN.
- Primary audit v4, preflight v8, cooccurrence v7, final dataset/report v5
  refuse overwrite; failed v6 evidence is separate.
- v9/v10 reviews bind superseded digests and cannot authorize this source set.

Any shell/Python path disagreement, unhandled matched branch, stale canonical
input, v4 output redirection, v6 reuse, broad M2 status prefix, denominator
omission, negative recoding, or stale-review acceptance is a blocker.

The v5 authority JSON, approval, and signature must be absent before review.
Do not edit/create/chmod/delete/sign, run pytest or producers, install packages,
or use network. Use bounded direct reads and hash commands only.

Return exactly one valid JSON object and no prose or code fence. Do not repeat
top-level keys. Use exactly:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only with no blocker.
