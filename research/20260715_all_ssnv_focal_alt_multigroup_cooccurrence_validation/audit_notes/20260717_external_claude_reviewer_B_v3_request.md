<!--
建立時間: 2026-07-17
目標: portable layout QA isolation advisory 修補後的最終科學統計來源稽核
處理範圍: Task Type B 7 datasets / 469,849 sSNV scientific claim contract
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
-->

# External Claude Code Reviewer B v3 Request

You are independent reviewer B v3. Work read-only in
`/big7_disk/liaoyoyo2001/InterSubMod`; do not modify any file.

This final review supersedes B v1/v2. Since v2, the only protected-source delta
is adding `-I` to the two terminal desktop/mobile portable-layout QA Python
invocations in `run_m2v5_recovered_completion_chain.sh`. Confirm this closes A
v2's advisory without changing M1/M2/G1/G2/R1/B1/C1 definitions, denominators,
statistics, selection, final-result bytes, or claim ceiling. The canonical
topic suite remains 342 passed, 0 failed.

Read `AGENTS.md`, `claim-contract-v5.md`, `scripts/release_source_authority.py`,
all 23 `EXPECTED_SOURCE_PATHS`, the v2 A/B outputs, and the v4 canonical pytest
XML. Independently recompute the source digest and verify:

- expected digest `442a741caccc2e6253120e7eb05cc364be007406b4c177cf71c40a79f28666a5`;
- expected Git HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`;
- 23/23 protected sources mode `0444` and unchanged claim-contract-v5;
- all science/claim-boundary safeguards from B v2 remain intact;
- no remaining finding blocks signing and executing the full 7-dataset,
  469,849-site Task Type B run.

Return concrete findings and these exact final fields:

```text
Verdict=APPROVE|REJECT
review_scope=<explicit scope>
findings_closed=true|false
reviewed_source_set_sha256=<64 hex>
reviewed_git_head=<40 hex>
reviewer_id=5483591d-842f-4880-9d40-1332100439f1
remaining_findings=<none or concrete findings>
```

`APPROVE` means the source design and claim boundaries are ready for the formal
run; it does not assert that the biological hypothesis is confirmed.
