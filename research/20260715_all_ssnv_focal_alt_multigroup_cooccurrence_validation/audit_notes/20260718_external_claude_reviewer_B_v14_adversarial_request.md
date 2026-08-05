<!--
建立時間: 2026-07-18
目標: 外部 Reviewer B 對 v3/v7/v14 closure 與 Task-B release chain 做 fresh adversarial 唯讀審查
處理範圍: race/failure injection、producer-consumer binding、科學語義、stale-source/key isolation
關聯檔案: ../logs/pytest_full_pre_authority_v18_v3_signer_v14_normalizer_v7_assembler_canonical.xml
-->

# External Reviewer B v14: adversarial ceremony and release review

Perform a fresh, independent, adversarial, read-only P0 review. Do not rely on
v13 findings. Do not edit, chmod, sign, execute producers, run pytest, install
packages, or use the network. Bounded reads, stat, hashing, Git inspection, and
non-mutating shell checks are allowed.

Authoritative source enumeration:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py::EXPECTED_SOURCE_PATHS`

Review signer v3, assembler v7, normalizer v14 and all corresponding tests in
the same paths listed by Reviewer A.

Fixed anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Protected roles: 23, all mode `0444`
- Source-set SHA-256:
  `1d7166f0a192848a9b6ad812e93dac4404b65caaffaa6509fb53156c2ca8eab4`
- Canonical JUnit:
  `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v18_v3_signer_v14_normalizer_v7_assembler_canonical.xml`
- JUnit SHA-256:
  `2bce0db29e92864a740df7b2fa1492d49a624fe5fe5ca5ba3f6f55a0e998dff7`
- JUnit counts: 544 tests, 0 failures, 0 errors, 0 skipped
- v3 signer source:
  `8243ced80242ace2ceae4ae71dad7c685a11f8ebe754dec56e7a1baa509c69d1`
- v7 assembler source:
  `113345e9574725787c67f3cec888aa67328188746b2d2e5daa4ebdf5a4bc6573`
- v14 normalizer source:
  `a4592b19a2cf56c9b50882f413f41cf1bd20d7596e919d7d831d9612aa73f89f`
- Source public-key SHA-256:
  `e4a09d9e43f76efed408330f43cee66ac1fde457f5b4ee6889d11338935e5b6c`
- Source key-directory/public/private modes: `0700/0444/0400`
- Result/report/supplemental public-key SHA-256 values:
  `2d37e58d837c216061ac50c6b123442593a2f87b28db71d162090c98a463d318`,
  `3fb508f62ee2f476fe217f9ade8d3ee7ffcf05d431cbddb7e00e0ec92f79e585`,
  `a67d41ff36af6c3ff92c776ed04c2256cc0547b59d1db50013bab66f5920aca3`
- Result/report/supplemental key-directory/public/private modes:
  `0700/0444/0400`
- v2 archive ledger SHA-256:
  `49860985d5c20c3801f2bb92d7124ea97d079f64f6a500a418244d697b05b94b`
- Live v5 authority, approval, signature and pending signature are absent.

Adversarial checks:

1. Try to disprove signer fail-closed semantics for target replacement,
   preexisting/racing signature output, signature chmod failure, directory
   fsync failure, publication failure, verification failure and private-key
   retirement failure. Confirm tests cover the important observable states.
2. Trace every signer subprocess and prove the binary, key, target and
   signature bytes are the same bound FDs that were hashed or verified.
3. Try to find any assembler identity/payload split read, path reopen,
   `importlib` execution, mutable reviewer payload, JUnit count trust, output
   overwrite, or directory replacement window that can authorize different
   bytes than those recorded.
4. Confirm v14 rejects `NaN`, `Infinity`, `-Infinity`, duplicate keys, extra
   top-level keys, trailing data, wrong digest/HEAD/UUID/reviewer identity and
   inconsistent approval states.
5. Verify the failed first v2 archive attempt ledger is preserved and the v2
   batch-complete ledger proves four no-signature moves; inspect all four
   per-key archive records and all live v3 key modes.
6. Independently recompute the exact 23-role digest using real newlines.
7. Trace every option-to-value/path binding across both formal runners,
   finalizers, matched-normal, CN/CCF, report delivery and QA. Confirm current
   versions remain primary-pre v4, preflight v8, cooccurrence v7, independent
   M2 v3, matched controls/analysis v3, CN/CCF v3, final dataset/report v5.
8. Confirm new key rotation changed only intended trust anchors and did not
   alter M2 pre-axis semantics, K>10 treatment, row conservation, G1/G2/B1
   denominators, four-state compatibility, negative-control interpretation,
   truth-blind selection, raw-read identity, HP/PS handling, full-scope counts,
   or claim ceiling.
9. Confirm v13 reviews bind `8dac...` rather than `1d716...` and therefore
   cannot authorize the new source set.

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
  findings.
- Digest and Git fields equal the fixed anchors.
- `evidence` is a JSON object.

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
