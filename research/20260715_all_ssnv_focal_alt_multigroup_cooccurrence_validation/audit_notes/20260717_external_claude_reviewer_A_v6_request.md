<!--
建立時間: 2026-07-17
目標: 對 v3 source authority 的真實 FD Ed25519 驗證與 fail-closed release chain 做唯讀複核
處理範圍: P0 pre-run source authorization；不包含尚未允許產生的 runtime artifacts
關聯檔案: ../logs/pytest_full_pre_authority_v8_real_fd_signature.xml
-->

# External Reviewer A v6: v3 real-FD source-integrity review

Perform a fresh bounded read-only review. Do not edit, create, chmod, delete, or execute a formal producer.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Topic: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Expected Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`

Expected exact 23-role source-set digest: `19c373f7895034e303deda11f9b62dd17bf8df5f07d9e33b3953b65560a173ec`

Canonical test XML: `logs/pytest_full_pre_authority_v8_real_fd_signature.xml`, SHA-256 `259a5e0007b3e4a748acdb0377383f7cf09eba8c98589bec58ce0020d66167c0`, `380 passed`, `0 failed`, mode `0444`.

## Phase definition

This is P0 pre-run source authorization. The v3 authority, approval and detached signature must not exist until both independent APPROVE records have been embedded. Formal result/report artifacts are also absent because their producers must fail closed until this source authority validates. Do not treat these phase-correct absences as blockers.

The signed-but-validator-incompatible v2 authority is historical evidence only. It was archived under `docs/archive/2026/07/20260717_signed_but_validator_incompatible_all_ssnv_source_authority_v2/` and must not be accepted by any protected source.

## Required checks

1. Independently recompute the source digest; confirm exactly 23 roles and all mode `0444`.
2. Confirm all protected source trust anchors require authority ID `20260717_all_ssnv_focal_alt_task_b_release_v3`, v3 authority paths, public-key SHA-256 `60ebac3ee2ebfbf69a80331b40410b365c0315d41363d59ee0b44f6dbf5040e4`, and accept neither v1 nor v2.
3. Confirm `verify_ed25519_signature_fds` invokes pinned `/usr/bin/openssl pkeyutl -verify -rawin` against seekable `/proc/self/fd/*` descriptors retained by `BoundArtifactReader`, passes all descriptors with `pass_fds`, and does not use pipe-backed `/dev/stdin`.
4. Confirm both source-authority and result/report receipt validators use the real-FD helper; confirm real-key positive and tampered-payload negative integration tests exist and passed.
5. Confirm the one-time signer signs and verifies the target path before publication, retains the key on failure, and retires it only after success.
6. Reconfirm F1 producer authority/command/source-lock binding and F2 terminal report binding remain fail closed.

Return only one valid JSON object with exactly these required top-level fields: `schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`, `verdict`, `findings_closed`, `f1_status`, `f2_status`, `reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`, `blocking_findings`, `nonblocking_findings`, `evidence`. Use `schema_name=intersubmod.external_claude_source_review`, `schema_version=1.0.0`, a fresh UUID as `reviewer_id`, and `verdict=APPROVE` only if no blocking finding remains. Approval covers source implementation only, not absent biological results.
