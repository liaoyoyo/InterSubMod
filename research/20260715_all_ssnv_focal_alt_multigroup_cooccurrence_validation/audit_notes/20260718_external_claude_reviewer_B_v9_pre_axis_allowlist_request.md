<!--
建立時間: 2026-07-18
目標: 以攻擊者視角檢查 K>10 allowlist、authority rollback 與全量重跑契約
處理範圍: P0 pre-run source authorization；只讀、不得執行 producer
關聯檔案: claim-contract-v5.md, source authority validator, cooccurrence/final release runners
-->

# External Reviewer B v9: adversarial pre-axis allowlist release review

Perform a fresh independent adversarial review of the exact live source. The
previous formal computation reached the summary phase after all workers
completed, then failed closed because one 11-group site was intentionally
rejected by the M2 K<=10 gate before axis classification. Determine whether the
repair closes only that totality defect without converting unmeasured axes into
negative evidence or creating a release-authority bypass.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Start from the live `EXPECTED_SOURCE_PATHS` in
`scripts/release_source_authority.py`, then inspect the relevant protected
sources, especially:

- `scripts/analyze_methyl_ssnv_cooccurrence.py`
- `scripts/m2_screen_gate.py`
- `scripts/audit_cooccurrence_task_contract_preflight.py`
- `scripts/finalize_cooccurrence_release_receipt.py`
- `scripts/build_all_ssnv_final_report_dataset.py`
- `scripts/run_cooccurrence_v6_source_locked.sh`
- `scripts/run_m2v5_recovered_completion_chain.sh`
- `scripts/finalize_task_b_result_release.py`
- `tests/test_analyze_methyl_ssnv_cooccurrence.py`
- `tests/test_terminal_validation_chain.py`
- `logs/pytest_full_pre_authority_v12_pre_axis_allowlist_canonical.xml`
- `claim-contract-v5.md`

Independently verify:

- HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- Exactly 23 protected roles, all mode `0444`, source digest
  `088ccd47ecee4462876c0ce4bb7a4f054f50d28dc52ce995434a0b0cd9221295`.
- JUnit digest
  `2a6489fb4c2386331e7e3db0684c524b148c9e4d2cd9b40f83822b7041a79551`,
  449 tests, zero failures/errors/skips.
- Public-key digest
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  private key is `0400` before the one allowed signature.
- Every axis census has exactly one status per selected site. An empty map is
  allowed only for `NOT_EVALUABLE_M2_SCREEN_FIELDS_MISSING`,
  `NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM`, or
  `NOT_EVALUABLE_M2_CATEGORICAL_LEVEL_COUNTS_MISSING`.
- A site stopped by one of those three pre-axis gates is visibly labelled
  `NOT_EVALUATED_M2_GATE_PRE_AXIS:<gate-status>` and cannot be interpreted as
  effect-below-threshold, p-value failure, aligned, or unaligned evidence.
- Empty maps with an eligible/evaluable or post-axis status, including
  `NOT_EVALUABLE_M2_AXIS_INDETERMINATE`, non-object evidence, missing status,
  partial axes, extra axes, and non-conservation fail closed.
- The only observed K>10 site is not promoted to M2/G1/G2 and downstream
  contracts keep group-count >10 as not evaluable.
- v7/v8/v5 paths are new and reject overwrite; the failed v6 log remains
  separate. Fresh v5 source authority supersedes but cannot masquerade as v4.
- Source identity, real-newline digest, direct canonical commands, raw-read
  identity preflight, result signature, and report signature gates remain
  fail-closed. No rollback to v1-v4 can satisfy v5 validation.

Treat omission of an unmeasured pre-axis category, acceptance by a broad status
prefix, silent exclusion from the denominator, or recoding as negative evidence
as blockers. Treat manual reuse or overwriting of the failed v6 path as a
blocker.

The v5 authority artifacts are expected to be absent at P0. Do not reuse prior
reviews. Do not edit/create/chmod/delete/sign, run pytest, launch formal
producers, install packages, or use network. Use bounded direct reads only.

Return exactly one valid JSON object and no prose or code fence. Do not repeat
any top-level key. Use exactly these top-level fields:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if no blocker remains.
