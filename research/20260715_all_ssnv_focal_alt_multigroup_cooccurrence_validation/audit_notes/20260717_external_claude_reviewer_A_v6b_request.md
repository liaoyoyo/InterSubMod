<!--
建立時間: 2026-07-17
目標: 以明確 bounded IO 重作 v3 source-integrity review
處理範圍: P0 pre-run source authorization；禁止廣域檔案系統探索
關聯檔案: 20260717_external_claude_reviewer_A_v6_failed_unbounded_search.json
-->

# External Reviewer A v6b: bounded v3 real-FD review

Perform the same independent review specified in `20260717_external_claude_reviewer_A_v6_request.md`, with these mandatory execution bounds:

- Work only under `/big7_disk/liaoyoyo2001/InterSubMod` and the exact paths named by that request or by `EXPECTED_SOURCE_PATHS` in `scripts/release_source_authority.py`.
- Do not run `find`, recursive `grep`, recursive `ls`, filesystem-wide discovery, `conda`, package installation, network access, plugin discovery, or any full pytest run.
- Do not edit, create, chmod, delete, sign, or execute a formal producer.
- Use `/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I` only when Python is needed.
- Verify the existing JUnit XML by reading/hash/parsing it; do not regenerate it.
- Use no more than 16 tool calls. A missing non-required convenience tool is not a reason to search the filesystem.

The exact v3 facts to independently verify remain:

- Git HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
- 23 protected roles, all `0444`, digest `19c373f7895034e303deda11f9b62dd17bf8df5f07d9e33b3953b65560a173ec`.
- JUnit XML `logs/pytest_full_pre_authority_v8_real_fd_signature.xml`, SHA-256 `259a5e0007b3e4a748acdb0377383f7cf09eba8c98589bec58ce0020d66167c0`, 380 tests and zero failures.
- Authority ID `20260717_all_ssnv_focal_alt_task_b_release_v3`; v3 public-key SHA-256 `60ebac3ee2ebfbf69a80331b40410b365c0315d41363d59ee0b44f6dbf5040e4`.
- Seekable `/proc/self/fd/*` Ed25519 verification with retained descriptors and `pass_fds`; no `/dev/stdin` signature verification.
- Both source and result/report validators use the helper; real-key positive and tampered-payload negative integration tests are present.
- F1/F2 remain fail closed; archived v1/v2 cannot authorize execution.

This is P0 pre-run. Absent v3 authority/signature and formal runtime artifacts are phase-correct and must not be blockers. Post-run review remains required.

Return only one valid JSON object with exactly these top-level fields: `schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`, `verdict`, `findings_closed`, `f1_status`, `f2_status`, `reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`, `blocking_findings`, `nonblocking_findings`, `evidence`. Use `schema_name=intersubmod.external_claude_source_review`, `schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if no blocker remains.
