<!--
建立時間: 2026-07-17
目標: 外部 Claude Code 對正式 Task-B source release 的安全與 provenance 做唯讀複核
處理範圍: 23 protected sources、F1 producer binding、F2 report release signature
關聯檔案: ../logs/pytest_full_pre_authority_v5_producer_report_release.xml
-->

# External Reviewer A v4: release security and provenance

Perform a fresh read-only audit. Do not edit, create, chmod, or delete any file.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Topic: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Expected Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`

Expected 23-role source-set digest: `33fd92e5589ce727a2addf2d10639a3eec00486c69b3ebb95172e8c6156717db`

Canonical test XML: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v5_producer_report_release.xml` (`376 passed`, `0 failed`; 38 deprecation warnings).

Independently recompute the digest from `release_source_authority.EXPECTED_SOURCE_PATHS`, confirm exactly 23 roles and all modes `0o444`, and confirm claim-contract-v5 remains 9,144 bytes, SHA-256 `da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af`, mode `0o444`.

Audit the two previously blocking findings:

1. F1 producer binding. Check `run_matched_normal_candidate_controls.py`, `analyze_matched_normal_candidate_controls.py`, `annotate_candidate_cn_ccf.py`, and `build_all_ssnv_final_report_dataset.py`. A formal producer must validate the same signed source authority before output; reject injected argv; require exact `/proc/self/cmdline` with `python -I`; capture unchanged live source identity/mode; publish a mode-0444 terminal receipt; and the final builder must independently require exact live authority, producer code, canonical command and source lock. Self-consistent arbitrary producer receipts must not pass. Confirm wrong-authority, wrong-code and missing-`-I` regressions.

2. F2 terminal report release. Check `finalize_task_b_result_release.py`, `build_all_ssnv_report_artifact.py`, `deliver_portable_report_scrollbar_safe.mjs`, `qa_portable_report_layout.py`, and `run_m2v5_recovered_completion_chain.sh`. The second one-time Ed25519 report receipt/signature must bind the already signed final dataset release plus report.md, artifact.json, report build receipt, portable HTML, portable delivery receipt, official/desktop/mobile screenshots, and desktop/mobile QA. All must be mode 0444; create/wait/verify must occur after QA; verification must rehash every artifact; HTML and QA tamper tests must fail closed.

Also look for bypasses, circular trust, signing-order errors, mutable unsigned deliverables, missing path/size/SHA bindings, and post-signature artifacts. This is pre-run source approval, not biological confirmation.

Return only the requested structured result. `APPROVE` only when both original blockers are closed and no release-blocking issue remains.

