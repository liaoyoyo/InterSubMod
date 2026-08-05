<!--
建立時間: 2026-07-18
目標: External Reviewer B adversarial review of v5/v16/v9 ceremony closure and the full Task-B release chain
處理範圍: path/FD/content races, failure truth, stale-key isolation, producer bindings, and scientific semantics
關聯檔案: ../logs/pytest_full_pre_authority_v20_v5_signer_v16_normalizer_v9_assembler_canonical.xml
-->

# External Reviewer B v16: adversarial ceremony and release review

Perform a fresh, independent, adversarial, read-only P0 review. Do not rely on
prior verdicts. Do not edit, chmod, sign, execute producers, run pytest,
install packages, or use web/network tools. Bounded reads, stat, hashing, and
Git inspection are allowed.

Authoritative source enumeration:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py::EXPECTED_SOURCE_PATHS`

Review the v5 signer, v16 normalizer, v9 assembler, production validator,
their tests, and the v4 key archive script/ledger listed in Reviewer A's
request.

Fixed anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Protected roles: 23, all mode `0444`
- Source-set SHA-256:
  `9a7e9ebfca0c3e0daceffdc30def073b10247d80502ce37a56e2b812e5203924`
- Canonical JUnit SHA-256:
  `95b354f77025d53e01a4f294bfd7147cf0cce87fda8401132ddd7d3432204eb9`
- Canonical JUnit counts: 658 tests, 0 failures, 0 errors, 0 skipped
- v5 signer SHA-256:
  `f39e35eedd7ea02fba0632d18fac29928b78831f55e370c5ccfebfb94426cabb`
- v16 normalizer SHA-256:
  `5dfb663be024b8e8f4c6d6c8f265b24a4a73340989cfc8794bc2b6836a2ba251`
- v9 assembler SHA-256:
  `3fc260bcc36f69caef8bf5619d5cd1d150751a848fdd8fa02c83f2ad0cc24a59`
- Source/result/report/supplemental public-key SHA-256 values:
  `677f73c56f8e3fb998ca4f1e360dd58dd73b90c6389f2f3d3d17350c2f5ec698`,
  `1c7abc3a62da417ed6c94ab47ec28e58f3a9c0c4249f413a5997cadc53e366c3`,
  `0dcf55f8a532ac9f506e1aba703fc3caee1a0cde5de602bfd5dd5a3d35e40a88`,
  `7413426516497930be3f73c42bac451a5e0c85021f1cd7660e1c7fc7a5e1cf61`
- All four live key directory/public/private modes: `0700/0444/0400`
- v4 archive ledger SHA-256:
  `b7720d1ab680f9f87b23f238a2f9295df1f78acf9419a1133fb44f492028a9a9`
- Live v5 authority, approval, signature, and pending signature are absent.

Adversarial checks:

1. Trace the exact third `_verify_signature()` path in signer v5. Model a
   rename/unlink plus replacement of the canonical signature path after that
   verification but after private mode becomes `000`. Confirm the new terminal
   pathname/FD and directory checks reject success, with truthful
   post-retirement failure state.
2. Model target replacement, key-directory chmod/replacement, publication
   collision, signature chmod/fsync failure, and private fsync failure. Confirm
   each state is fail closed and no stale record can be reported as `SIGNED`.
3. In normalizer v16, model parent replacement, output entry replacement,
   preexisting output, partial write, and an already-open writable FD doing a
   same-inode/same-size pwrite after the second path check. Confirm terminal
   metadata plus byte reread rejects every split.
4. Replace, delete, chmod, duplicate-key edit, non-finite edit, UUID edit, or
   semantic edit of either v16 review after assembly. Confirm production
   validation reopens and binds both fixed files through final validation.
5. Verify the v4 archive ledger has four never-signed no-replace completions,
   chains to v3, and old public hashes are absent from protected/live
   consumers. Confirm historical v8 cannot accept the rotated release module.
6. Independently recompute the exact 23-role digest using real newlines.
7. Trace every option-to-value/path binding through both formal runners,
   finalizers, matched-normal, CN/CCF, report delivery, and layout QA. Current
   versions must remain primary-pre v4, preflight v8, cooccurrence v7,
   independent M2 v3, matched controls/analysis v3, CN/CCF v3, and final
   dataset/report v5.
8. Confirm the ceremony/security rotation did not alter truth-blind M1,
   M2 pre-axis/K>10 semantics, row conservation, G1/G2/B1 denominators,
   four-state compatibility, negative-control interpretation, raw-read
   identity, same-run LongPhase-S HP/PS handling, 7-dataset/469,849-site
   scope, or the claim ceiling that forbids calling candidates cellular
   subclones or lineage topology.
9. Parse JUnit independently and require attributes/elements both 658/0/0/0.
   Check focused post-fix evidence: signer 17, normalizer 19, validator 11,
   assembler v9 post-anchor 12.

Return exactly one JSON object with these 15 top-level keys and no others:

`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Contract:

- `schema_name`: `intersubmod.external_claude_source_review`
- `schema_version`: `1.0.0`
- `reviewer_label` contains `Reviewer B`
- `reviewer_id`: a fresh canonical lowercase UUIDv4
- `model`: starts with `claude-`
- `verdict`: `APPROVE` or `REQUEST_CHANGES`
- Clean approval requires `findings_closed=true`,
  `f1_status=RESOLVED_VERIFIED`, `f2_status=RESOLVED_VERIFIED`, and no blocking
  findings
- Digest and Git fields equal the fixed anchors
- `evidence` is a JSON object

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
