<!--
建立時間: 2026-07-17
目標: CLI provenance 修補後的外部 Claude Code reviewer B 科學統計重新稽核
處理範圍: Task Type B 7 datasets / 469,849 sSNV scientific claim contract
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
-->

# External Claude Code Reviewer B v2 Request

You are independent reviewer B v2. Work read-only in
`/big7_disk/liaoyoyo2001/InterSubMod`; do not modify any file.

This review supersedes the earlier B review because the protected source set
changed after the main agent found and fixed a formal `main(argv)` command-
provenance bypass. Confirm that the fix changes execution provenance only and
does not alter the scientific gates, denominators, statistics, or claim ceiling.

Review the scientific and statistical integrity of the formal Task Type B
source set for all seven datasets and all 469,849 latest LongPhase-S
recalibrated `FILTER=PASS` sSNVs. Read `AGENTS.md`, the topic
`claim-contract-v5.md`, `scripts/release_source_authority.py`, and all 23 paths
in `EXPECTED_SOURCE_PATHS`.

Focus on:

- exact population denominators and TP/FP/UNASSESSED separation;
- M1 and M2 definitions, including axis-indeterminate, K>10, low-power, and aligned-below-negative-power states;
- exact latest HP/PS read-tag use before focal-ALT selection and read identity reconciliation;
- sSNV cooccurrence, complete-read markers, four-state compatibility, and technical-replicate limits;
- strict-stage multiple testing, effect floors, robustness, runtime replay, tumor-REF, matched-normal, and CN/CCF confounds;
- whether methyl groups under a shared ancestral ALT are described only as latent molecular substructure / subclone candidates rather than proof of clones, linear ancestry, or a phylogenetic tree;
- whether G2/B1/C1 claims preserve NOT_EVALUABLE and compatibility-only semantics;
- whether the final builder can independently reject producer-consistent but scientifically invalid outputs;
- whether the CLI-provenance patch leaves all scientific behavior and tests intact.

Independently recompute the source-set digest using
`release_source_authority.source_set_digest`. The expected live value at
dispatch is
`c9017942c9ac526df317f4a21bbefca71dff554143471f7e76212978e399ab9e` and
the expected Git HEAD is `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`.
The local canonical topic test result at dispatch is 342 passed, 0 failed.

Return concrete blocking and non-blocking findings, then these exact fields:

```text
Verdict=APPROVE|REJECT
review_scope=<explicit scope>
findings_closed=true|false
reviewed_source_set_sha256=<64 hex>
reviewed_git_head=<40 hex>
reviewer_id=d7f7a3ed-c08e-419c-a4e2-0431443a11b9
remaining_findings=<none or concrete findings>
```

`APPROVE` means the design and claim boundaries are adequate to begin the
formal run. It does not mean the biological hypothesis is already confirmed.
