<!--
建立時間: 2026-07-18
目標: 獨立檢查所有受保護 producer 的 current release path/command alignment
處理範圍: P0 pre-run source authorization；只讀與 bounded 重算
關聯檔案: all 23 EXPECTED_SOURCE_PATHS, v15 canonical JUnit, claim-contract-v5
-->

# External Reviewer A v12: all protected producers

Perform a fresh independent, read-only P0 review of the exact live source.
Prior reviews found path drift at three layers before signature: final
dataset/report outputs, final-builder inputs, and leaf-producer direct-command
constants. No v5 authority artifact has ever been created or signed.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Begin from all 23 live `EXPECTED_SOURCE_PATHS` in
`scripts/release_source_authority.py`. At minimum inspect:

- `run_cooccurrence_v6_source_locked.sh`
- `finalize_cooccurrence_release_receipt.py`
- `run_m2v5_recovered_completion_chain.sh`
- `build_all_ssnv_final_report_dataset.py`
- `finalize_task_b_result_release.py`
- `run_matched_normal_candidate_controls.py`
- `analyze_matched_normal_candidate_controls.py`
- `annotate_candidate_cn_ccf.py`
- `build_all_ssnv_report_artifact.py`
- `analyze_methyl_ssnv_cooccurrence.py`
- `m2_screen_gate.py`
- tests for release authority, terminal chain, cooccurrence finalizer,
  matched normal, CN/CCF, builder, and M2 summary
- `logs/pytest_full_pre_authority_v15_leaf_producer_alignment_canonical.xml`
- `claim-contract-v5.md`
- the external v5 Ed25519 public key

Required anchors:

- HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- 23 protected roles, all `0444`, source digest
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`.
- JUnit digest
  `b671f9b92d081a318e5ffe60beb60ee60ced010ff27d8424ab9ee9ab2c3ff94e`,
  457 tests, zero failures/errors/skips.
- Public-key digest
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  private key `0400` pre-signature.

Verify all of the following:

- Cooccurrence release finalizer binds primary-pre v4, preflight v8, producer
  v7 summary-hotfix, and the same exact files passed by its shell runner.
- Matched-normal runner binds cooccurrence v7 candidate table and controls-v3
  output. Matched analyzer binds controls-v3 input and analysis-v3 output.
  Their canonical direct commands equal completion-runner invocations.
- CN/CCF annotator binds cooccurrence v7 input and CN/CCF-v3 output.
- Final builder binds all v3/v4/v7/v8 inputs, two matched branches, and v5
  dataset. Finalizer accepts exactly the two builder commands and report v5.
- No protected source contains any superseded exact token listed by
  `test_protected_release_sources_contain_no_superseded_path_tokens`.
  Independently repeat that search rather than trusting only the test.
- Shell option/value bindings are tested, not only presence of version names.
- Zero-candidate matched NOT_APPLICABLE and nonzero analysis branches both
  remain executable and distinguishable.
- Direct CLI, source checks, real-newline digest, overwrite refusal, raw-read
  identity preflight, final result signature, and final report signature stay
  fail-closed.
- M2 accepts an empty map only for three true pre-axis gates, conserves every
  row, rejects malformed schemas, and leaves K>10 as M1 PASS / M2
  NOT_EVALUABLE / G1-G2-B1 NOT_RUN.
- Failed v6 evidence is separate; current v7/v8/v5 outputs are fresh.
- Reviews for all earlier digests cannot authorize the current source set.

Any stale producer constant, shell/direct-command mismatch, unhandled branch,
v4 redirect, v6 reuse, denominator omission, negative recoding, or stale-review
acceptance is a blocker.

The v5 authority JSON, approval, and signature must be absent. Do not
edit/create/chmod/delete/sign, run pytest or producers, install packages, or use
network. Use bounded direct reads and hash/stat commands only.

Return exactly one valid JSON object and no prose/code fence, with unique
top-level keys:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only with no blocker.
