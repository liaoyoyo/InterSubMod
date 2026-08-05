<!--
建立時間: 2026-07-17
目標: 外部 Claude Code 對 Task-B 科學口徑、統計 gate 與報告完整性做唯讀複核
處理範圍: 7 datasets / 469,849 sites、claim ceiling、producer/report release delta
關聯檔案: ../logs/pytest_full_pre_authority_v5_producer_report_release.xml
-->

# External Reviewer B v4: scientific contract and report completeness

Perform a fresh independent read-only review. Do not edit, create, chmod, or delete any file.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Topic: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Expected Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`

Expected 23-role source-set digest: `33fd92e5589ce727a2addf2d10639a3eec00486c69b3ebb95172e8c6156717db`

Canonical test XML: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v5_producer_report_release.xml` (`376 passed`, `0 failed`; warnings are SciPy/NumPy deprecations).

Independently recompute the digest and inspect the implementation delta. Confirm that the matched-normal/CN provenance changes and independent report-signature layer alter release evidence only, not M1/M2/G1/G2/R1/B1/C1 selection definitions, denominators, statistical thresholds, multiple-testing logic, four-state compatibility logic, normal/tumor-REF/CN interpretation, or the claim ceiling.

Required scientific boundaries:

- Input remains all 7 datasets and 469,849 latest LongPhase-S recalibrated `FILTER=PASS` chr1-22 biallelic sSNV dataset-sites.
- M1 is descriptive focal-ALT methyl multigroup evidence, not subclone confirmation.
- M2 retains asymmetric measured-axis evaluability, including axis-indeterminate, K>10, low-power and aligned-below-negative-power states.
- G2 remains a multi-marker molecular-haplotype candidate only.
- Four-state and cooccurrence topology outputs remain compatibility/identifiability statements, not proof of linear ancestry, clone count, a cellular subclone, or a phylogenetic tree.
- Matched-normal non-replication, tumor-REF controls and CN/PyClone annotations retain NOT_EVALUABLE/conditional semantics and cannot independently form cellular C1.
- Report release signing must not turn execution-integrity PASS into biological confirmation.

Review whether the report release covers every user-facing final artifact and whether report/portable QA tamper tests are adequate. Check that the full suite really records 376 tests and zero failures, source roles are exactly 23 and mode 0444, and claim-contract-v5 is unchanged at SHA-256 `da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af`, size 9,144, mode 0444.

Return only the requested structured result. `APPROVE` means source design and claim boundaries are adequate to begin the formal full run; it does not assert a biological result.

