<!--
建立時間: 2026-07-18
目標: 獨立檢查 v5 producer path alignment、pre-axis allowlist 與 canonical 測試證據
處理範圍: P0 pre-run source authorization；只讀與 bounded 重算
關聯檔案: final dataset builder, result finalizer, completion runner, M2 summary tests
-->

# External Reviewer A v10: v5 output-path alignment and source review

Perform a fresh independent, read-only P0 review of the exact live source.
The v9 review approved the explicit three-status pre-axis allowlist but found
that the completion shell declared v5 final dataset/report paths while the
Python dataset builder and release finalizer still pinned v4 canonical
directories. The authority was not created or signed. The two Python constants
were changed to v5, a cross-producer path-alignment regression was added, and
the previously shared extra-axis fail-closed branch now has a dedicated test.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Inspect these exact files:

- `scripts/build_all_ssnv_final_report_dataset.py`
- `scripts/finalize_task_b_result_release.py`
- `scripts/run_m2v5_recovered_completion_chain.sh`
- `scripts/analyze_methyl_ssnv_cooccurrence.py`
- `scripts/m2_screen_gate.py`
- `scripts/release_source_authority.py`
- `scripts/run_cooccurrence_v6_source_locked.sh`
- `tests/test_analyze_methyl_ssnv_cooccurrence.py`
- `tests/test_terminal_validation_chain.py`
- `logs/pytest_full_pre_authority_v13_output_path_alignment_canonical.xml`
- `claim-contract-v5.md`
- `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v5_summary_hotfix/ed25519_public.pem`

Independently verify:

- Git HEAD is `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID is `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- Exactly 23 protected roles are declared, all live sources are mode `0444`,
  and the real-newline source-set digest is
  `4685ff85bbad8c3174dcfcd3b60dcd7d20a940a6caf11d43558ca079e22029b1`.
- Canonical JUnit SHA-256 is
  `a2d4e94ac9cf2cc411602d28ae4c8734904351340b923ff63a29bf2e1c308dab`;
  it records 451 tests and zero failures, errors, or skips.
- Public-key SHA-256 is
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  the private key remains `0400` before the one permitted signature.
- The completion runner, dataset builder
  `CANONICAL_FINAL_DATASET_DIR`, and result finalizer
  `FINAL_DATASET_DIR` all resolve to
  `all_ssnv_final_report_dataset_v5_m2v5_source_attested`.
- The completion runner and result finalizer `REPORT_DIR` both resolve to
  `all_ssnv_final_report_v5_m2v5_source_attested`. The affected live sources
  contain no v4 canonical dataset/report directory constant.
- `test_m2v5_dataset_and_report_paths_align_across_producers` would fail on the
  prior v4/v5 mismatch and is present in the canonical suite.
- The M2 empty-map allowlist remains exactly
  `SCREEN_FIELDS_MISSING`,
  `GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM`, and
  `CATEGORICAL_LEVEL_COUNTS_MISSING`; accepted rows use
  `NOT_EVALUATED_M2_GATE_PRE_AXIS:<full-status>` without invented statistics.
- Dedicated tests reject partial axes, extra axes, and an empty map paired with
  post-axis `NOT_EVALUABLE_M2_AXIS_INDETERMINATE`.
- Formal output versions remain fresh and non-overwriting: primary audit v4,
  preflight v8, cooccurrence v7 summary-hotfix, final dataset/report v5.
- Both runners bind only v5 authority and retain source checks before, during,
  and after computation. v9 reviews bind the superseded source digest and are
  not eligible for the new authority.

The v5 authority JSON, approval, and signature must be absent before review.
Do not reuse any prior verdict. Do not edit, create, chmod, delete, sign,
launch a signer, run pytest, run formal producers, install packages, or use
network. Use direct reads and bounded hash commands only.

Return exactly one valid JSON object and no prose or code fence. Do not repeat
any top-level key. Use exactly these top-level fields:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if no blocker remains.
