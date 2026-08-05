<!--
建立時間: 2026-07-17
目標: 關閉 reviewer A v2 portable layout QA isolation advisory 後的最終來源稽核
處理範圍: Task Type B 7 datasets / 469,849 sSNV release source set
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
-->

# External Claude Code Reviewer A v3 Request

You are independent reviewer A v3. Work read-only in
`/big7_disk/liaoyoyo2001/InterSubMod`; do not modify any file.

This is the final source review and supersedes A v1/v2. A v2 approved the
formal source set but reported one non-blocking advisory: the two terminal
desktop/mobile `qa_portable_report_layout.py` invocations used `QA_PYTHON`
without `-I`. Both invocations now include `-I`, and a regression assertion
requires every `${PYTHON}` and `${QA_PYTHON}` invocation in the completion
runner to be isolated. The full canonical topic suite remains 342 passed,
0 failed.

Independently verify:

1. the v2 advisory is closed and no Python invocation in
   `run_m2v5_recovered_completion_chain.sh` lacks `-I`;
2. the prior formal `main(argv)` bypass remains closed by authority-aware argv
   rejection plus exact `/proc/self/cmdline` matching;
3. the exact 23-role source set is mode `0444`, the claim contract is unchanged,
   and the digest independently recomputes to
   `442a741caccc2e6253120e7eb05cc364be007406b4c177cf71c40a79f28666a5`;
4. Git HEAD is `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`;
5. no remaining issue blocks source signing or the formal full-scope run.

Read `AGENTS.md`, `claim-contract-v5.md`, `scripts/release_source_authority.py`,
all 23 `EXPECTED_SOURCE_PATHS`, the v2 A output, and the v4 canonical pytest XML.

Return concrete findings and these exact final fields:

```text
Verdict=APPROVE|REJECT
review_scope=<explicit scope>
findings_closed=true|false
reviewed_source_set_sha256=<64 hex>
reviewed_git_head=<40 hex>
reviewer_id=33c27da9-b56b-4f10-b6e7-25684e1b54e0
remaining_findings=<none or concrete findings>
```

Give `APPROVE` only if the source set is ready to sign and execute. This is a
pre-run execution/release approval, not biological confirmation.
