<!--
建立時間: 2026-07-18
目標: 對所有受保護 producer 做攻擊面 stale-path/direct-command 稽核
處理範圍: P0 pre-run source authorization；只讀、不得執行 producer
關聯檔案: release authority, both runners, all leaf producers, v15 tests
-->

# External Reviewer B v12: adversarial all-producer review

Perform a fresh independent adversarial review of the exact live source.
Assume that a green builder test can still hide stale direct-command constants
inside leaf producers. Trace every current completion step from shell arguments
through producer canonical constants, receipts, downstream readers, and both
signature finalizers.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Start with all 23 `EXPECTED_SOURCE_PATHS`; explicitly inspect:

- cooccurrence runner + cooccurrence release finalizer
- completion runner + final dataset builder + result/report finalizer
- matched-normal runner + matched-normal analyzer
- CN/CCF annotator, strict producer, primary/frozen auditors
- report builder and portable delivery/QA sources
- M2 analyzer/gate and all new path/stale-token tests

Required anchors:

- HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- 23 roles, all `0444`, source digest
  `db0c33f89c31e6b1a3e4cd68a2943df3c109f0e76ddaff40e9e0d2bf78faa3ed`.
- v15 JUnit digest
  `b671f9b92d081a318e5ffe60beb60ee60ced010ff27d8424ab9ee9ab2c3ff94e`,
  457 tests and zero failures/errors/skips.
- Public key digest
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  private key `0400` before signing.

Attack these contracts:

- Cooccurrence finalizer must use primary-pre v4, preflight v8, and producer v7.
- Matched runner/analyzer must use v7 and controls/analysis v3; both zero and
  nonzero branches must survive formal direct-command enforcement.
- CN/CCF must use v7 input and v3 output.
- Builder/finalizer must use only current inputs, two matched commands, dataset
  v5, report v5, and signed receipts bound to those paths.
- Independently search every protected source for the exact stale tokens in
  `test_protected_release_sources_contain_no_superseded_path_tokens`.
- Check shell options are attached to the correct values, not merely that the
  same versioned name appears somewhere in each file.
- Ensure no alternate CLI, arbitrary matched directory, older authority,
  stale review digest, v6 producer dir, v4 final dir, or overwrite can pass.
- Ensure M2 pre-axis labels cannot become negative evidence and K>10 cannot
  enter a PASS/FAIL denominator or G1/G2/B1.
- Ensure result/report signatures, private-key retirement, source identity,
  raw-read identity, and browser QA remain required.

Treat any stale constant anywhere in the 23-source set, branch-dependent
failure, direct-command mismatch, missing test coverage, denominator omission,
negative recoding, rollback, or stale-review acceptance as a blocker.

The v5 authority artifacts must be absent. Do not reuse prior verdicts. Do not
edit/create/chmod/delete/sign, run pytest/producers, install packages, or use
network. Use bounded direct reads only.

Return exactly one valid JSON object and no prose/code fence, with unique
top-level keys:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only with no blocker.
