<!--
建立時間: 2026-07-17
目標: 外部 Claude Code 對 v2 修正後的 Task-B 科學口徑與報告完整性做唯讀複核
處理範圍: 7 datasets / 469,849 sites、claim ceiling、signer-only delta
關聯檔案: ../logs/pytest_full_pre_authority_v7_fixed_signer.xml
-->

# External Reviewer B v5: scientific-contract non-regression

Perform a fresh independent read-only review. Do not edit, create, chmod, or delete any file.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Topic: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Expected Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`

Expected exact 23-role source-set digest: `b00a1c3605520af5fcc314d9c55b9a00ef90c651244b8c3666a30a67349a6add`

Canonical test XML: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v7_fixed_signer.xml` (SHA-256 `580fa83b75f5167189f4ea50f3939829e7fe022bc218db214ee07a64985f34b0`, mode `0o444`; `377 passed`, `0 failed`).

Independently recompute the source digest and inspect the v2 delta. Confirm it changes only cryptographic trust-anchor paths/hashes, authority ID, fail-safe one-time signer implementation, and matching regression assertions. It must not change M1/M2/G1/G2/R1/B1/C1 selection definitions, denominators, thresholds, multiple-testing logic, four-state compatibility logic, normal/tumor-REF/CN interpretation, or the claim ceiling.

Required scientific boundaries:

- Input remains all 7 datasets and 469,849 latest LongPhase-S recalibrated `FILTER=PASS` chr1-22 biallelic sSNV dataset-sites.
- M1 is descriptive focal-ALT methyl multigroup evidence, not subclone confirmation.
- M2 retains asymmetric measured-axis evaluability, including axis-indeterminate, K>10, low-power and aligned-below-negative-power states.
- G2 remains a multi-marker molecular-haplotype candidate only.
- Four-state and cooccurrence topology remain compatibility/identifiability evidence, not proof of linear ancestry, clone count, a cellular subclone, or a phylogenetic tree.
- Matched-normal non-replication, tumor-REF controls and CN/PyClone annotations retain NOT_EVALUABLE/conditional semantics and cannot independently form cellular C1.
- Cryptographic PASS is execution-integrity evidence only and cannot turn a candidate into biological confirmation.

Check source roles are exactly 23 and mode `0444`; claim-contract-v5 remains SHA-256 `da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af`, size 9,144, mode `0444`; the report release still covers all ten user-facing artifacts; and the v1 unsigned attempt is excluded rather than silently reused.

Return only the requested structured result. `APPROVE` means the v2 source design and scientific claim boundaries are adequate to begin the formal full run; it does not assert a biological result.

