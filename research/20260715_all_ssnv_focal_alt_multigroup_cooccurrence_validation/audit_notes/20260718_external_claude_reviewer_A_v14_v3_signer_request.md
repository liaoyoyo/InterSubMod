<!--
建立時間: 2026-07-18
目標: 外部 Reviewer A 對 Raman blocker closure、v3 signer、v7 assembler 與新 23-role source set 做 fresh 唯讀審查
處理範圍: FD binding、failure propagation、key rotation、canonical tests、authority pre-signature state、release contracts
關聯檔案: ../logs/pytest_full_pre_authority_v18_v3_signer_v14_normalizer_v7_assembler_canonical.xml
-->

# External Reviewer A v14: FD-bound ceremony closure

Perform a fresh, independent, read-only P0 review. Do not reuse any v13 verdict
or evidence. Do not edit, chmod, sign, run producers, run pytest, install
packages, or use the network. Bounded reads, stat, hashing, Git inspection, and
non-mutating shell checks are allowed.

Review all 23 live `EXPECTED_SOURCE_PATHS` from:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`

Also review these ceremony controls and their tests:

- `audit_notes/run_one_time_ed25519_signer_v3.py`
- `tests/test_one_time_ed25519_signer_v3.py`
- `audit_notes/assemble_release_source_authority_v7.py`
- `tests/test_assemble_release_source_authority_v7.py`
- `audit_notes/normalize_external_source_review_v14.py`
- `tests/test_normalize_external_source_review_v14.py`

Required anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Exactly 23 protected roles, all mode `0444`
- Source-set SHA-256:
  `1d7166f0a192848a9b6ad812e93dac4404b65caaffaa6509fb53156c2ca8eab4`
- Canonical JUnit:
  `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v18_v3_signer_v14_normalizer_v7_assembler_canonical.xml`
- JUnit SHA-256:
  `2bce0db29e92864a740df7b2fa1492d49a624fe5fe5ca5ba3f6f55a0e998dff7`
- JUnit counts: 544 tests, 0 failures, 0 errors, 0 skipped
- v3 signer source SHA-256:
  `8243ced80242ace2ceae4ae71dad7c685a11f8ebe754dec56e7a1baa509c69d1`
- v7 assembler source SHA-256:
  `113345e9574725787c67f3cec888aa67328188746b2d2e5daa4ebdf5a4bc6573`
- v14 normalizer source SHA-256:
  `a4592b19a2cf56c9b50882f413f41cf1bd20d7596e919d7d831d9612aa73f89f`
- The three sources above are mode `0444`.
- Source public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v5_summary_hotfix/ed25519_public.pem`
- Source public-key SHA-256:
  `e4a09d9e43f76efed408330f43cee66ac1fde457f5b4ee6889d11338935e5b6c`
- Source key-directory/public/private modes before signing: `0700/0444/0400`
- Result public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/20260717_all_ssnv_result_v2/ed25519_public.pem`
- Result public-key SHA-256:
  `2d37e58d837c216061ac50c6b123442593a2f87b28db71d162090c98a463d318`
- Report public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/20260717_all_ssnv_report_v2/ed25519_public.pem`
- Report public-key SHA-256:
  `3fb508f62ee2f476fe217f9ade8d3ee7ffcf05d431cbddb7e00e0ec92f79e585`
- Supplemental public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_supplemental_report_integrity/20260718_positional_singleton_v1/ed25519_public.pem`
- Supplemental public-key SHA-256:
  `a67d41ff36af6c3ff92c776ed04c2256cc0547b59d1db50013bab66f5920aca3`
- Result/report/supplemental key-directory/public/private modes before signing:
  `0700/0444/0400`.
- v2 key archive ledger:
  `audit_notes/20260718_unused_v2_signer_key_archive.v2.jsonl`
- v2 key archive ledger SHA-256:
  `49860985d5c20c3801f2bb92d7124ea97d079f64f6a500a418244d697b05b94b`
- Current v5 authority JSON, approval JSON, detached signature, and pending
  signature are all absent before assembly/signing.

Explicitly close the three prior high-severity engineering findings:

1. Signer v3 propagates every write/chmod/fsync/link/verify/retirement failure
   as a non-zero error. It never prints `SIGNED` after a failed step, and the
   encrypted private key remains mode `0400` unless final publication and the
   second FD-bound verification have completed.
2. Signer v3 creates key directories at `0700`, uses exclusive key files,
   executes a pre-hashed bound OpenSSL FD, signs/verifies the same target,
   public-key and anonymous signature FDs, and publishes the signature with
   FD-derived `linkat` no-replace semantics.
3. Assembler v7 executes the release module and v14 normalizer from already
   bound bytes, parses each review from the same bytes used for its identity,
   retains every input/output FD through final checks, parses JUnit counts
   independently, and creates outputs through bound parent-dir FDs plus
   `O_EXCL`. It must contain no v6 two-read loader.

Also verify:

4. v14 rejects duplicate keys, trailing prose and non-finite JSON constants at
   any depth.
5. All four unused v2 key directories were never signed, were archived without
   deletion or overwrite, and cannot be reused at the live paths.
6. Independently recompute the 23-role digest with real newline delimiters.
7. The four new v3 public-key hashes are correctly bound by all live
   consumers; old hashes are absent from every protected source and live
   non-historical supplemental consumer.
8. The rotation changed exactly the intended trust-anchor consumers and did
   not alter M1/M2/G1/G2/R1/B1/C1 definitions, denominators, thresholds,
   multiple testing, four-state compatibility, matched-normal/tumor-REF/CN
   interpretation, raw-read identity, full scope, or claim ceiling.
9. Every formal command remains direct-CLI, source-attested, non-overwriting,
   full 7-dataset/469,849-site scope, and bound to the same authority ID.

Return exactly one JSON object with these 15 top-level keys and no others:

`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Contract:

- `schema_name`: `intersubmod.external_claude_source_review`
- `schema_version`: `1.0.0`
- `reviewer_label` contains `Reviewer A`
- `reviewer_id`: a fresh canonical lowercase UUIDv4
- `model`: starts with `claude-`
- `verdict`: `APPROVE` or `REQUEST_CHANGES`
- Clean approval requires `findings_closed=true`,
  `f1_status=RESOLVED_VERIFIED`, `f2_status=RESOLVED_VERIFIED`, and an empty
  `blocking_findings` array.
- Digest and Git fields equal the anchors above.
- `evidence` is a JSON object.

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
