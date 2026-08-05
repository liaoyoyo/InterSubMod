<!--
建立時間: 2026-07-18
目標: 獨立檢查 K=11 pre-axis allowlist、summary totality 與 canonical 測試證據
處理範圍: P0 pre-run source authorization；只讀與 bounded 重算
關聯檔案: scripts/analyze_methyl_ssnv_cooccurrence.py, scripts/run_cooccurrence_v6_source_locked.sh, scripts/release_source_authority.py
-->

# External Reviewer A v9: pre-axis allowlist and source review

Perform a fresh independent, read-only P0 review of the exact live source.
The previous full-scope cooccurrence attempt completed all worker computation
but failed closed during summary construction. Exactly one selected M1 site had
11 methyl groups. The M2 gate correctly stopped that site before per-axis
classification and returned an empty `m2_axis_statuses` object, but the old
summary unconditionally indexed every axis and raised `KeyError: 'hp_exact'`.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Inspect these exact files:

- `scripts/analyze_methyl_ssnv_cooccurrence.py`
- `scripts/m2_screen_gate.py`
- `scripts/release_source_authority.py`
- `scripts/run_cooccurrence_v6_source_locked.sh`
- `scripts/run_m2v5_recovered_completion_chain.sh`
- `tests/test_analyze_methyl_ssnv_cooccurrence.py`
- `tests/test_terminal_validation_chain.py`
- `logs/pytest_full_pre_authority_v12_pre_axis_allowlist_canonical.xml`
- `claim-contract-v5.md`
- `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v5_summary_hotfix/ed25519_public.pem`

Independently verify:

- Git HEAD is `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID is `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- `release_source_authority.py` declares exactly 23 protected roles, all live
  sources are mode `0444`, and the real-newline source-set digest is
  `088ccd47ecee4462876c0ce4bb7a4f054f50d28dc52ce995434a0b0cd9221295`.
- Canonical JUnit XML SHA-256 is
  `2a6489fb4c2386331e7e3db0684c524b148c9e4d2cd9b40f83822b7041a79551`;
  it records 449 tests and zero failures, errors, or skips.
- Public-key SHA-256 is
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  the private key is still mode `0400` before the one permitted signature.
- `summarize_m2_axis_status_counts` derives the exact axis set from
  `M2_GATE.AXIS_SPECS` and accepts an empty axis object only for these three
  explicit pre-axis statuses:
  `NOT_EVALUABLE_M2_SCREEN_FIELDS_MISSING`,
  `NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM`, and
  `NOT_EVALUABLE_M2_CATEGORICAL_LEVEL_COUNTS_MISSING`.
- An accepted pre-axis row is labelled
  `NOT_EVALUATED_M2_GATE_PRE_AXIS:<full-gate-status>` for every axis. It does
  not invent an effect, p-value, alignment, or negative result.
- Empty maps with a post-axis status such as
  `NOT_EVALUABLE_M2_AXIS_INDETERMINATE`, partial or extra axes, malformed
  evidence, missing statuses, and non-conservation all fail closed.
- Positive and negative tests cover the 11-group pre-axis case and the
  post-axis-empty, partial-axis, and extra-axis schema defects.
- New formal paths do not overwrite the failed v6 evidence: cooccurrence output
  is v7 summary-hotfix, preflight is v8, fresh primary audit is v4, and final
  dataset/report paths are v5. The wrapper success message names v7.
- Both runners bind only v5 authority and retain source checks before, during,
  and after computation. Existing v4 artifacts remain historical evidence only.

The v5 authority JSON, approval, and signature are correctly absent before
review and are not blockers. Do not reuse any prior verdict.

Do not edit, create, chmod, delete, sign, launch a signer, run pytest, run formal
producers, install packages, or use network. Use direct reads and bounded hash
commands only.

Return exactly one valid JSON object and no prose or code fence. Do not repeat
any top-level key. Use exactly these top-level fields:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if no blocker remains.
