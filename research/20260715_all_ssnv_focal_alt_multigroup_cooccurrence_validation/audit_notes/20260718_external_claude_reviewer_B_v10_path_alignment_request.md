<!--
建立時間: 2026-07-18
目標: 以攻擊者視角檢查 v4/v5 path rollback、M2 allowlist 與 release authority
處理範圍: P0 pre-run source authorization；只讀、不得執行 producer
關聯檔案: claim-contract-v5.md, source authority validator, completion/final release sources
-->

# External Reviewer B v10: adversarial path-alignment release review

Perform a fresh independent adversarial review of the exact live source. The
v9 review found a fail-closed but operationally fatal mismatch between v5 shell
output paths and v4 Python canonical paths. Determine whether the repair closes
that mismatch across the complete dataset/report/signature chain without
weakening source identity, overwrite refusal, M2 denominator semantics, or the
one-time release authority.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Start from live `EXPECTED_SOURCE_PATHS` in
`scripts/release_source_authority.py`, then inspect the relevant protected
sources, especially:

- `scripts/build_all_ssnv_final_report_dataset.py`
- `scripts/finalize_task_b_result_release.py`
- `scripts/run_m2v5_recovered_completion_chain.sh`
- `scripts/analyze_methyl_ssnv_cooccurrence.py`
- `scripts/m2_screen_gate.py`
- `scripts/audit_cooccurrence_task_contract_preflight.py`
- `scripts/finalize_cooccurrence_release_receipt.py`
- `scripts/release_source_authority.py`
- `tests/test_analyze_methyl_ssnv_cooccurrence.py`
- `tests/test_terminal_validation_chain.py`
- `logs/pytest_full_pre_authority_v13_output_path_alignment_canonical.xml`
- `claim-contract-v5.md`

Independently verify:

- HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID `20260718_all_ssnv_focal_alt_task_b_release_v5`.
- Exactly 23 protected roles, all mode `0444`, source digest
  `4685ff85bbad8c3174dcfcd3b60dcd7d20a940a6caf11d43558ca079e22029b1`.
- JUnit digest
  `a2d4e94ac9cf2cc411602d28ae4c8734904351340b923ff63a29bf2e1c308dab`,
  451 tests, zero failures/errors/skips.
- Public-key digest
  `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`;
  private key is `0400` before the one allowed signature.
- Shell, builder, and finalizer agree on dataset v5; shell and finalizer agree
  on report v5. No live canonical constant can redirect create, verify,
  signature, or report verification to v4.
- CLI equality guards, `require_absent`, output-directory exclusive creation,
  signed result receipt, and signed report receipt all bind those same v5
  directories. A v4 path or v1-v4 authority cannot satisfy v5 validation.
- The new cross-producer regression checks both positive v5 equality and
  absence of prior v4 constants.
- Every M2 axis census conserves every selected row. Empty axis status is
  accepted only for the three explicit pre-axis gates and is visibly labelled
  as not evaluated; partial axes, an extra axis, a post-axis empty map,
  malformed evidence, missing status, or non-conservation fails closed.
- The one K>10 site remains M1 PASS, M2 NOT_EVALUABLE, G1/G2/B1 NOT_RUN, and
  cannot enter a PASS/FAIL denominator through path or summary changes.
- Primary audit v4, preflight v8, cooccurrence v7, and final dataset/report v5
  are fresh paths that refuse overwrite; the failed v6 log remains separate.
- v9 review records bind digest `088ccd47...` and cannot authorize digest
  `4685ff85...`.

Treat any surviving v4 canonical create/verify path, shell/Python output
disagreement, broad status-prefix acceptance, silent denominator omission,
negative-evidence recoding, stale-review acceptance, or v6 output reuse as a
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
