<!--
建立時間: 2026-07-18
目標: 對完整 Task-B source-attested chain 做攻擊面 path/receipt/denominator 稽核
處理範圍: P0 pre-run source authorization；只讀、不得執行 producer
關聯檔案: release authority, completion runner, builder/finalizer, M2 gate
-->

# External Reviewer B v11: adversarial complete-chain review

Perform a fresh independent adversarial review of the exact live source.
Previous reviews found two operationally fatal but fail-closed mismatches:
v4/v5 output drift and v2/v6 versus v3/v4/v7/v8 input drift. Determine whether
the complete current chain is now executable without weakening any scientific
or signature gate. Include both matched-normal terminal branches.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Start from live `EXPECTED_SOURCE_PATHS` in `release_source_authority.py`.
Inspect all path-producing and path-consuming protected sources, especially:

- `run_m2v5_recovered_completion_chain.sh`
- `build_all_ssnv_final_report_dataset.py`
- `finalize_task_b_result_release.py`
- `run_matched_normal_candidate_controls.py`
- `analyze_matched_normal_candidate_controls.py`
- `build_all_ssnv_report_artifact.py`
- `analyze_methyl_ssnv_cooccurrence.py`
- `m2_screen_gate.py`
- `run_cooccurrence_v6_source_locked.sh`
- `release_source_authority.py`
- the three relevant test files and v14 JUnit XML
- `claim-contract-v5.md`

Required anchors:

- HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- 23 protected roles, all `0444`, source digest
  `20232fc5433fdc5a5210c731bfd3b2cbf8ed0c1c04abc2781196ea67b08aad12`.
- JUnit digest
  `511366ce5fe27a8382820b99fd1ddbf6d199c32f1ec58f1865e637ab0bd427fd`,
  452 tests, zero failures/errors/skips.
- Public-key digest
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  private key `0400` pre-signature.

Adversarially verify:

- Every completion-runner builder argument has the same canonical resolved path
  in the builder. No fixed input is left at v2/v6.
- Matched-normal analysis-v3 and controls-v3 are the only accepted alternatives.
  The zero-candidate NOT_APPLICABLE branch can build and be signed; an arbitrary
  third directory or a noncanonical direct command cannot.
- Finalizer accepts only the two builder-generated canonical commands and binds
  dataset v5/report v5 through create, verify, both detached signatures, and
  report verification.
- No stale constant can redirect final dataset, report, candidate tables,
  receipts, screenshots, or QA files to v4.
- `require_absent`, exclusive create, source checks before/after, real-newline
  digest, raw-read identity preflight, result receipt, and report receipt all
  remain fail-closed.
- The exact three-status M2 pre-axis allowlist conserves all rows without
  inventing evidence; extra/partial/post-axis-empty/malformed states fail.
- K>10 is not promoted or silently omitted.
- v6 failed log is preserved, v7/v8/v5 paths are new, and reviews for digest
  `4685ff85...` or earlier cannot satisfy the current authority.
- Tests cover the full fixed input map, two matched branches, two commands,
  stale-name absence, v5 outputs, and M2 schema negatives.

Treat any unreviewed alternate command, branch-dependent path failure,
shell/Python mismatch, stale input, v4 redirect, v6 reuse, denominator
omission, negative recoding, authority rollback, or stale-review acceptance as
a blocker.

The v5 authority artifacts are expected absent. Do not reuse prior verdicts.
Do not edit/create/chmod/delete/sign, run pytest or producers, install packages,
or use network. Use bounded direct reads only.

Return exactly one valid JSON object and no prose or code fence. Do not repeat
top-level keys. Use exactly:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only with no blocker.
