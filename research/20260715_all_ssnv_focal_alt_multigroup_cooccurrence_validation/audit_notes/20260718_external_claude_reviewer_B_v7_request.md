<!--
建立時間: 2026-07-18
目標: 以攻擊者視角獨立檢查 v4 修補後的完整 release chain
處理範圍: P0 pre-run source authorization；只讀、不得執行 producer
關聯檔案: claim-contract-v5.md, source authority, result/report finalizers
-->

# External Reviewer B v7: adversarial v4 release-chain review

Perform a fresh independent adversarial review of the exact live source. A prior
v3 runner failed closed before computation because shell aggregate hashing used
a literal `\n` delimiter while Python used real newlines. Confirm that v4 fixes
the defect without creating an authorization bypass, mutable-source window, or
rollback path.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Start from the exact live authority module's `EXPECTED_SOURCE_PATHS`, then inspect
the relevant paths among its 23 roles, especially:

- `scripts/release_source_authority.py`
- `scripts/run_cooccurrence_v6_source_locked.sh`
- `scripts/run_m2v5_recovered_completion_chain.sh`
- `scripts/audit_cooccurrence_task_contract_preflight.py`
- `scripts/audit_stable_primary_artifacts.py`
- `scripts/analyze_methyl_ssnv_cooccurrence.py`
- `scripts/finalize_cooccurrence_release_receipt.py`
- `scripts/finalize_task_b_result_release.py`
- `tests/test_terminal_validation_chain.py`
- `tests/test_release_source_authority_fd_signature.py`
- `logs/pytest_full_pre_authority_v10_newline_digest_fix_canonical.xml`
- `claim-contract-v5.md`
- `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v4/ed25519_public.pem`

Independently verify:

- HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID `20260718_all_ssnv_focal_alt_task_b_release_v4`.
- 23 roles, all `0444`, source digest
  `9542e7f0f1d12794f7b7736e106f927556d1cc95612d4e99576bb8dd33521a17`.
- Shell and Python source-set calculations both use real newline delimiters and
  yield the same digest.
- Public-key digest
  `b71855f5fd9d0e97df0f6186420b5cec95f85d8b462fde0a890443846271bee4`;
  private key is `0400` before the one permitted signature.
- Canonical XML digest
  `8d2fc24d70fc2b23694cc3fbc8eb5fed65105951cef686533141d644d7bf853e`,
  381 tests, zero failures/errors/skips.
- v1/v2/v3 cannot satisfy v4 validation.
- The new primary-audit v3 and preflight v7 paths avoid overwriting v3-run
  evidence and are consumed by both formal runners.
- Source identity is checked before, during, and after producer execution.
- F1 result release and F2 report release remain detached-signature,
  source-authority, exact-command, and artifact-identity fail-closed gates.
- The already-running result/report signers cannot sign the source approval,
  and the new source signer cannot sign result/report receipts.

Treat any manual bypass of the failed v3 shell check as unacceptable. The v4
authority artifacts are expected to be absent at P0. Do not reuse prior reviews.

Read-only constraints: do not edit/create/chmod/delete/sign, run pytest, launch
formal producers, install packages, or use network. Use direct absolute-path
reads; no recursive `find`. Use at most 20 tool calls.

Return only one valid JSON object with exactly these top-level fields:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if no blocker remains.
