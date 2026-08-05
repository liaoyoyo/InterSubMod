<!--
建立時間: 2026-07-17
目標: 修補 formal CLI provenance bypass 後的外部 Claude Code reviewer A 重新稽核
處理範圍: Task Type B 7 datasets / 469,849 sSNV release source set
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
-->

# External Claude Code Reviewer A v2 Request

You are independent reviewer A v2. Work read-only in
`/big7_disk/liaoyoyo2001/InterSubMod`; do not modify any file.

This review supersedes the earlier A review. The main agent rejected the prior
approval after independently finding that the cooccurrence analyzer and
preflight recorded a synthetic canonical command and could be called through
`main(argv)` from an isolated host process. The source set was changed to:

- include `-I` in both canonical commands;
- reject injected `argv` under the formal source authority;
- read and exactly compare `/proc/self/cmdline` for formal execution;
- run every Python stage in the completion runner with `-I`;
- add regression tests for the import bypass and missing `-I` case.

Review the formal Task Type B source set for the full seven-dataset, 469,849-site
LongPhase-S recalibrated `FILTER=PASS` sSNV analysis. Read `AGENTS.md`, the topic
`claim-contract-v5.md`, `scripts/release_source_authority.py`, and every path in
its `EXPECTED_SOURCE_PATHS` mapping.

Focus on release engineering and reproducibility:

- independently reproduce the former `main(argv)` command-provenance bypass and
  verify the patched formal path now rejects it before output creation;
- verify exact `python -I` plus `/proc/self/cmdline` enforcement for preflight,
  analyzer, primary/strict/final producers and release finalizers;
- exact 23-role source identity and mode `0444`;
- one-time Ed25519 source approval and final-result receipt, including private-key retirement;
- same-descriptor / same-byte verification and path/inode binding;
- fixed Python, Git, OpenSSL, Node, Chromium, runtime-library hashes and thread environment;
- strict-statistic execution replay over every strict output field;
- independent recount of 102,842 stable sites and 308,526 primary artifacts;
- fail-closed full scope of 469,849 sites and immutable formal outputs;
- whether any synthetic test seam remains reachable by a formal CLI.

Independently recompute the source-set digest with the exact
`release_source_authority.source_set_digest` implementation. The expected live
value at dispatch is
`c9017942c9ac526df317f4a21bbefca71dff554143471f7e76212978e399ab9e` and
the expected Git HEAD is `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
The local canonical topic test result at dispatch is 342 passed, 0 failed.

Return a concise review with these exact final fields:

```text
Verdict=APPROVE|REJECT
review_scope=<explicit scope>
findings_closed=true|false
reviewed_source_set_sha256=<64 hex>
reviewed_git_head=<40 hex>
reviewer_id=e1aa2e28-21f4-4a09-83e8-6c85e23a2841
remaining_findings=<none or concrete findings>
```

Give `APPROVE` only when no finding blocks executing or trusting the formal
result. Approval concerns the execution/release contract, not proof of the
biological hypothesis.
