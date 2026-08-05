<!--
建立時間: 2026-07-17
目標: 以全絕對路徑修正 v6b wrong-root 誤判並重作 bounded source-integrity review
處理範圍: P0 pre-run source authorization；所有 live artifact 都以 absolute path 指定
關聯檔案: 20260717_external_claude_reviewer_A_v6b_rejected_wrong_relative_root.json
-->

# External Reviewer A v6c: absolute-path v3 real-FD review

Perform a fresh independent review. The previous v6b review resolved topic-relative paths against the repository root and therefore did not inspect the live implementation. Do not reuse its verdict.

Repository root:
`/big7_disk/liaoyoyo2001/InterSubMod`

Topic root:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Read the complete original scope from this exact path:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260717_external_claude_reviewer_A_v6_request.md`

Critical exact paths:

- Source authority module: `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`
- Result/report finalizer: `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/finalize_task_b_result_release.py`
- Real-FD tests: `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_release_source_authority_fd_signature.py`
- Terminal-chain tests: `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_terminal_validation_chain.py`
- Canonical JUnit XML: `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v8_real_fd_signature.xml`
- Claim contract: `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md`
- Public key: `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260717_all_ssnv_v3/ed25519_public.pem`

Do not use `git grep` to decide whether live content exists: this release explicitly permits a dirty/ignored research workspace only through exact live-content hashes. Use direct absolute-path reads/stats and the live module's `EXPECTED_SOURCE_PATHS` list. Do not use `find`, recursive grep/ls, conda discovery, package installation, network, pytest reruns, or formal producers. Do not edit/create/chmod/delete/sign anything. Use at most 16 tool calls.

Machine-verified facts to recompute, not merely trust:

- HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- 23 roles, all mode `0444`, digest `19c373f7895034e303deda11f9b62dd17bf8df5f07d9e33b3953b65560a173ec`.
- XML SHA-256 `259a5e0007b3e4a748acdb0377383f7cf09eba8c98589bec58ce0020d66167c0`, 380 tests, zero failures/errors.
- Public-key SHA-256 `60ebac3ee2ebfbf69a80331b40410b365c0315d41363d59ee0b44f6dbf5040e4`.
- Each of those three hexadecimal digests is exactly 64 characters; verify with code rather than visual counting.
- Authority ID `20260717_all_ssnv_focal_alt_task_b_release_v3`.

Check seekable `/proc/self/fd/*` Ed25519 verification, retained descriptors, `pass_fds`, no `/dev/stdin`, use by source and result/report validators, positive/tamper tests, F1/F2 fail-closed behavior, signer ordering, and archived v1/v2 exclusion. The source authority/signature and formal runtime outputs are correctly absent at P0 and are not blockers; post-run review remains mandatory.

Return only one valid JSON object with exactly these top-level fields: `schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`, `verdict`, `findings_closed`, `f1_status`, `f2_status`, `reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`, `blocking_findings`, `nonblocking_findings`, `evidence`. Use `schema_name=intersubmod.external_claude_source_review`, `schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if no blocker remains.
