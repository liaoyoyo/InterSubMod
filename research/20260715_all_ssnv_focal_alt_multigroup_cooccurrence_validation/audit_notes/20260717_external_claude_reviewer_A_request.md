<!--
建立時間: 2026-07-17
目標: 外部 Claude Code reviewer A 正式來源與可重現性稽核
處理範圍: Task Type B 7 datasets / 469,849 sSNV release source set
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
-->

# External Claude Code Reviewer A Request

You are independent reviewer A. Work read-only in
`/big7_disk/liaoyoyo2001/InterSubMod`; do not modify any file.

Review the formal Task Type B source set for the full seven-dataset, 469,849-site
LongPhase-S recalibrated `FILTER=PASS` sSNV analysis. Read `AGENTS.md`, the topic
`claim-contract-v5.md`, `scripts/release_source_authority.py`, and every path in
its `EXPECTED_SOURCE_PATHS` mapping.

Focus on release engineering and reproducibility:

- exact 23-role source identity and mode `0444`;
- one-time Ed25519 source approval and final-result receipt, including private-key retirement;
- same-descriptor / same-byte verification and path/inode binding;
- canonical `python -I` CLI and `/proc/self/cmdline` enforcement;
- fixed Python, Git, OpenSSL, Node, Chromium, runtime-library hashes and thread environment;
- strict-statistic execution replay over every strict output field;
- independent recount of 102,842 stable sites and 308,526 primary artifacts;
- fail-closed full scope of 469,849 sites and immutable formal outputs;
- whether synthetic test seams can be reached by a formal CLI.

Independently recompute the source-set digest with the exact
`release_source_authority.source_set_digest` implementation. The expected live
value at dispatch is
`855e206fdae0fa7549106a721b7c581a69cd2037e2a7050e0121554e5a9c2642` and
the expected Git HEAD is `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
The local topic test result at dispatch is 340 passed, 0 failed.

Return a concise review with these exact final fields:

```text
Verdict=APPROVE|REJECT
review_scope=<explicit scope>
findings_closed=true|false
reviewed_source_set_sha256=<64 hex>
reviewed_git_head=<40 hex>
reviewer_id=f7f56a76-f2ac-4627-8555-802c7ad5ce31
remaining_findings=<none or concrete findings>
```

Give `APPROVE` only when no finding blocks executing or trusting the formal
result. Approval concerns the execution/release contract, not proof of the
biological hypothesis.
