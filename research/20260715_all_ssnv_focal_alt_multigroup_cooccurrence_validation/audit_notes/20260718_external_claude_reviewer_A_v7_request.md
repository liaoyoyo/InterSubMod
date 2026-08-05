<!--
建立時間: 2026-07-18
目標: 獨立檢查 v4 source-set newline digest 修補與 canonical 測試證據
處理範圍: P0 pre-run source authorization；只讀與 bounded 重算
關聯檔案: scripts/release_source_authority.py, scripts/run_cooccurrence_v6_source_locked.sh, scripts/run_m2v5_recovered_completion_chain.sh
-->

# External Reviewer A v7: v4 newline digest and source authority review

Perform a fresh independent, read-only P0 review. A v3 source-authorized runner
failed closed before computation because its shell source-set digest used the
literal two-character delimiter `\n` instead of a real newline. The repair must
not bypass source authorization.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Inspect these exact files:

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_cooccurrence_v6_source_locked.sh`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_m2v5_recovered_completion_chain.sh`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_terminal_validation_chain.py`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_release_source_authority_fd_signature.py`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v10_newline_digest_fix_canonical.xml`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md`
- `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v4/ed25519_public.pem`

Machine claims to independently recompute:

- Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- Authority ID: `20260718_all_ssnv_focal_alt_task_b_release_v4`.
- Exactly 23 protected source roles, all mode `0444`.
- Python source-set digest:
  `9542e7f0f1d12794f7b7736e106f927556d1cc95612d4e99576bb8dd33521a17`.
- The corrected shell real-newline calculation must equal that same digest.
- Public-key SHA-256:
  `b71855f5fd9d0e97df0f6186420b5cec95f85d8b462fde0a890443846271bee4`.
- Private key exists at the sibling path and must still be mode `0400` at P0.
- JUnit XML SHA-256:
  `8d2fc24d70fc2b23694cc3fbc8eb5fed65105951cef686533141d644d7bf853e`.
- JUnit counts: 381 tests, zero failures/errors/skips.
- Both wrappers must reject v1/v2/v3 authorities and bind only v4.
- Primary audit and preflight outputs must use new non-overwriting v3/v7 names.
- The regression test must distinguish a real newline from a literal `\n`.

Check real-FD Ed25519 verification, retained descriptors, `pass_fds`, pinned
OpenSSL/Git identities, source modes, and fail-closed behavior. The v4 authority
JSON, approval, and signature are correctly absent before review and are not
blockers. Do not reuse any prior verdict.

Do not edit, create, chmod, delete, sign, launch a signer, run pytest, or run
formal producers. Do not use network, recursive `find`, or package installation.
Use at most 16 tool calls.

Return only one valid JSON object with exactly these top-level fields:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if no blocker remains.
