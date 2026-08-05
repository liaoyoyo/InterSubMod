<!--
建立時間: 2026-07-18
目標: External Reviewer A fresh review of v5 signer, v16 normalizer, v9 assembler, live-review reopen, key rotation, and the frozen Task-B source set
處理範圍: ceremony races, failure truth, source/key/JUnit anchors, producer contracts, and scientific claim guards
關聯檔案: ../logs/pytest_full_pre_authority_v20_v5_signer_v16_normalizer_v9_assembler_canonical.xml
-->

# External Reviewer A v16: terminal ceremony closure

Perform a fresh, independent, read-only P0 review. Do not reuse a prior
verdict. Do not edit, chmod, sign, execute producers, run pytest, install
packages, or use web/network tools. Bounded reads, stat, hashing, and Git
inspection are allowed.

Review all 23 live `EXPECTED_SOURCE_PATHS` from:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`

Also review:

- `audit_notes/run_one_time_ed25519_signer_v5.py`
- `tests/test_one_time_ed25519_signer_v5.py`
- `audit_notes/normalize_external_source_review_v16.py`
- `tests/test_normalize_external_source_review_v16.py`
- `audit_notes/assemble_release_source_authority_v9.py`
- `tests/test_assemble_release_source_authority_v9.py`
- `tests/test_release_source_authority_fd_signature.py`
- `audit_notes/archive_unused_v4_signer_keys.py`
- `audit_notes/20260718_unused_v4_signer_key_archive.v1.jsonl`

Required anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Exactly 23 protected roles, all mode `0444`
- Source-set SHA-256:
  `9a7e9ebfca0c3e0daceffdc30def073b10247d80502ce37a56e2b812e5203924`
- Release validator SHA-256:
  `3f56d23c4e23ad90aac6e5d6090c278e78b4cdc2f4ff648c9a57f01803ce076c`
- Canonical JUnit:
  `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v20_v5_signer_v16_normalizer_v9_assembler_canonical.xml`
- JUnit SHA-256:
  `95b354f77025d53e01a4f294bfd7147cf0cce87fda8401132ddd7d3432204eb9`
- JUnit counts: 658 tests, 0 failures, 0 errors, 0 skipped
- v5 signer SHA-256:
  `f39e35eedd7ea02fba0632d18fac29928b78831f55e370c5ccfebfb94426cabb`
- v16 normalizer SHA-256:
  `5dfb663be024b8e8f4c6d6c8f265b24a4a73340989cfc8794bc2b6836a2ba251`
- v9 assembler SHA-256:
  `3fc260bcc36f69caef8bf5619d5cd1d150751a848fdd8fa02c83f2ad0cc24a59`
- The three ceremony sources above are mode `0444`
- Source public-key SHA-256:
  `677f73c56f8e3fb998ca4f1e360dd58dd73b90c6389f2f3d3d17350c2f5ec698`
- Result public-key SHA-256:
  `1c7abc3a62da417ed6c94ab47ec28e58f3a9c0c4249f413a5997cadc53e366c3`
- Report public-key SHA-256:
  `0dcf55f8a532ac9f506e1aba703fc3caee1a0cde5de602bfd5dd5a3d35e40a88`
- Supplemental public-key SHA-256:
  `7413426516497930be3f73c42bac451a5e0c85021f1cd7660e1c7fc7a5e1cf61`
- Each live key directory/public/private mode before signing:
  `0700/0444/0400`
- v4 archive script SHA-256:
  `38b8d2a78087ad6be41f366cb920b34bb812e8a3dbdea0249bec8b2fd1cbdd0d`
- v4 archive ledger SHA-256:
  `b7720d1ab680f9f87b23f238a2f9295df1f78acf9419a1133fb44f492028a9a9`
- The v4 ledger chains to frozen v3 ledger SHA-256:
  `93bceb16d247f021a52fc6b17083dbbf2ba26cbc37da15d730b1de793a9f57e5`
- Current v5 authority, approval, detached signature, and pending signature
  are absent before assembly/signing.

Required closure checks:

1. Assembler v9 emits one canonical `SIGN` instruction containing the exact
   approval path, size, SHA-256, and mode `0o444`. Signer v5 parses and reopens
   the exact approved target; a same-path replacement cannot be signed.
2. Signer v5 retains the key-directory FD and checks path/inode/mode `0700`
   around OpenSSL use. Failure JSON reports actual publication, retirement
   start, and private-key mode.
3. After private-key retirement and the final signature verification, signer
   v5 rebinds the canonical signature pathname to the held signature FD and
   rechecks the parent directory. The regression that replaces the path after
   the third verify must produce `SIGN_FAILED`, never `SIGNED`.
4. Normalizer v16 uses bound parent/output FDs and `O_EXCL`, then checks final
   dev/inode/mode/size/mtime/ctime, rereads all terminal bytes, and hashes
   those bytes. Same-inode, same-size mutation is rejected.
5. The production validator reopens the fixed A/B v16 review files, requires
   exact signed artifact records and clean payloads, rejects duplicate and
   non-finite JSON, and independently requires canonical lowercase UUIDv4
   reviewer IDs.
6. The v4 ledger proves four unused, never-signed, no-replace archive moves;
   old live directories are absent and new v5 directories have the four
   anchored public keys.
7. Independently recompute the 23-role digest using real newline delimiters.
8. Confirm key rotation changed only intended trust anchors and review paths.
   M1/M2/G1/G2/R1/B1/C1 definitions, thresholds, denominators, four-state
   compatibility, matched-normal/tumor-REF/CN interpretation, raw-read
   identity, HP/PS handling, full scope, and claim ceiling remain unchanged.
9. Trace both formal runners and all downstream producers. Commands remain
   direct-CLI, source-attested, non-overwriting, and full
   7-dataset/469,849-site scope.
10. Independently parse JUnit attributes and testcase elements as 658/0/0/0.
    Separately verify the recorded focused runs: signer 17/17, normalizer
    19/19, validator 11/11, and post-anchor assembler 12/12.

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
  `blocking_findings` array
- Digest and Git fields equal the anchors above
- `evidence` is a JSON object

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
