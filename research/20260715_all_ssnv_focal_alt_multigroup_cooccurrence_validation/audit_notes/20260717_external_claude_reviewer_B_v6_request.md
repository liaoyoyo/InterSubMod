<!--
建立時間: 2026-07-17
目標: 對 v3 真實 FD 修正後的科學口徑與報告契約做獨立非回歸複核
處理範圍: 7 datasets / 469,849 sites、claim ceiling、cryptographic-only delta
關聯檔案: ../claim-contract-v5.md
-->

# External Reviewer B v6: v3 scientific-contract non-regression

Perform a fresh independent read-only review. Do not edit, create, chmod, delete, or execute a formal producer.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Topic: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Expected Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`

Expected exact 23-role source-set digest: `19c373f7895034e303deda11f9b62dd17bf8df5f07d9e33b3953b65560a173ec`

Canonical test XML: `logs/pytest_full_pre_authority_v8_real_fd_signature.xml`, SHA-256 `259a5e0007b3e4a748acdb0377383f7cf09eba8c98589bec58ce0020d66167c0`, mode `0444`; `380 passed`, `0 failed`.

This is a P0 pre-run source review. Formal runtime outputs are intentionally absent until signed source authorization validates; do not require them here. Post-run live-artifact review remains mandatory.

Independently recompute the source digest and compare the v3 delta. Confirm it is restricted to the seekable real-FD Ed25519 verification implementation, matching trust-anchor path/hash/ID changes, and regression tests. It must not change M1/M2/G1/G2/R1/B1/C1 selection definitions, denominators, thresholds, multiple-testing logic, four-state compatibility logic, matched-normal/tumor-REF/CN semantics, or the claim ceiling.

Required scientific boundaries:

- Input remains all 7 datasets and 469,849 latest LongPhase-S recalibrated `FILTER=PASS` chr1-22 biallelic sSNV dataset-sites.
- M1 is descriptive focal-ALT methyl multigroup evidence, not subclone confirmation.
- M2 retains asymmetric measured-axis evaluability, including axis-indeterminate, K>10, low-power and aligned-below-negative-power states.
- G2 remains a multi-marker molecular-haplotype candidate only.
- Four-state and cooccurrence topology remain compatibility/identifiability evidence, not proof of linear ancestry, clone count, a cellular subclone, or a phylogenetic tree.
- Matched-normal non-replication, tumor-REF controls and CN/PyClone annotations retain NOT_EVALUABLE/conditional semantics and cannot independently form cellular C1.
- Cryptographic PASS is execution-integrity evidence only and cannot upgrade biological certainty.

Also confirm claim-contract-v5 remains SHA-256 `da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af`, size 9,144, mode `0444`; the report release covers ten user-facing artifacts; and archived v1/v2 attempts cannot be reused.

Return only one valid JSON object with exactly these required top-level fields: `schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`, `verdict`, `findings_closed`, `f1_status`, `f2_status`, `reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`, `blocking_findings`, `nonblocking_findings`, `evidence`. Use `schema_name=intersubmod.external_claude_source_review`, `schema_version=1.0.0`, a fresh UUID as `reviewer_id`, and `verdict=APPROVE` only if no blocking finding remains. Approval means only that the exact v3 source design and scientific boundaries are adequate to begin the formal run.
