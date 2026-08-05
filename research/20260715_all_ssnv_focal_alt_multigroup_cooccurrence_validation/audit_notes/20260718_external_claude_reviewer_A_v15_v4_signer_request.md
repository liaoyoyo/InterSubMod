<!--
建立時間: 2026-07-18
目標: 外部 Reviewer A 對 Raman 三項 blocker closure、v4 signer、v8 assembler、v15 normalizer 與新 23-role source set 做 fresh 唯讀審查
處理範圍: digest handoff、key-directory lifecycle、failure-state truth、FD output binding、live-review re-open、canonical tests、release contracts
關聯檔案: ../logs/pytest_full_pre_authority_v19_v4_signer_v15_normalizer_v8_assembler_canonical.xml
-->

# External Reviewer A v15: digest-bound ceremony closure

Perform a fresh, independent, read-only P0 review. Do not reuse any v14 verdict
or evidence. Do not edit, chmod, sign, run producers, run pytest, install
packages, or use the network. Bounded reads, stat, hashing, Git inspection, and
non-mutating shell checks are allowed.

Review all 23 live `EXPECTED_SOURCE_PATHS` from:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`

Also review these ceremony controls and their tests:

- `audit_notes/run_one_time_ed25519_signer_v4.py`
- `tests/test_one_time_ed25519_signer_v4.py`
- `audit_notes/assemble_release_source_authority_v8.py`
- `tests/test_assemble_release_source_authority_v8.py`
- `audit_notes/normalize_external_source_review_v15.py`
- `tests/test_normalize_external_source_review_v15.py`
- `tests/test_release_source_authority_fd_signature.py`

Required anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Exactly 23 protected roles, all mode `0444`
- Source-set SHA-256:
  `b5e432fbd98e657ae48d6296fe3050eb01b3ccbc57ecc4aa0e5a25e37600ac9e`
- Canonical JUnit:
  `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v19_v4_signer_v15_normalizer_v8_assembler_canonical.xml`
- JUnit SHA-256:
  `094e64cd7c7c45e0e45c31f2485a8934ca2ebee9e2c6da92cee8f5afcfaf2179`
- JUnit counts: 586 tests, 0 failures, 0 errors, 0 skipped
- v4 signer source SHA-256:
  `6698393db61acf3c8fc242ac6e011eca6f96a128adf23f3d83664c95db654595`
- v8 assembler source SHA-256:
  `e1a43d7e30febd041759901aefd6bd68c5f7d4cbed25673606ee74cbf8ee2904`
- v15 normalizer source SHA-256:
  `840534e11f86fe04d22b2440c1830bdee0b1c774b3a6c1a6af9f4a99f160d0c0`
- The three sources above are mode `0444`.
- Source public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v5_summary_hotfix/ed25519_public.pem`
- Source public-key SHA-256:
  `ae8377ff1753e98e8490e44b58d41503d7484bb19704247be222123e10f4801a`
- Source key-directory/public/private modes before signing: `0700/0444/0400`
- Result public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/20260717_all_ssnv_result_v2/ed25519_public.pem`
- Result public-key SHA-256:
  `1c9cea01da7106ac502ce708831e8fd5b19854fe59956016b836301f100f0404`
- Report public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_report_release_authority/20260717_all_ssnv_report_v2/ed25519_public.pem`
- Report public-key SHA-256:
  `25face031089f413d564c93db4491a4433c81f97f339a616091beba077eca848`
- Supplemental public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_supplemental_report_integrity/20260718_positional_singleton_v1/ed25519_public.pem`
- Supplemental public-key SHA-256:
  `db20e9ca41cb78201578d01dd724a82d3894574a2ae63b583509dbdec9a1a6e8`
- Result/report/supplemental key-directory/public/private modes before signing:
  `0700/0444/0400`.
- v3 key archive ledger:
  `audit_notes/20260718_unused_v3_signer_key_archive.v1.jsonl`
- v3 key archive ledger SHA-256:
  `93bceb16d247f021a52fc6b17083dbbf2ba26cbc37da15d730b1de793a9f57e5`
- Prior v2 key archive ledger:
  `audit_notes/20260718_unused_v2_signer_key_archive.v2.jsonl`
- Prior v2 key archive ledger SHA-256:
  `49860985d5c20c3801f2bb92d7124ea97d079f64f6a500a418244d697b05b94b`
- Current v5 authority JSON, approval JSON, detached signature, and pending
  signature are all absent before assembly/signing.

Explicitly close the three prior high-severity engineering findings:

1. Assembler v8 emits one canonical v4 `SIGN` instruction containing the exact
   approval `path`, `size_bytes`, `sha256`, and mode `0o444`. Signer v4 parses
   that instruction canonically, opens the target once, and refuses a target
   whose live identity differs from the approved record. A same-path
   replacement after assembly therefore cannot be signed.
2. Signer v4 retains the key-directory FD from creation through signing,
   requires its path/inode/mode `0700` at every OpenSSL boundary, and refuses
   a post-READY chmod or replacement. It also reports the actual private-key
   mode, retirement-started state, and signature-publication state after every
   failure; an fsync error after `fchmod(000)` must not be described as a
   retained `0400` key.
3. Normalizer v15 creates output through a bound parent-directory FD plus
   `O_EXCL`, retains the output FD, and verifies FD/directory-entry/canonical
   path identity and read-back bytes. The live release validator re-opens both
   fixed v15 review files and requires their live artifact records and clean
   payload fields to equal both signed approval summaries and authority
   review-evidence records.

Also verify:

4. v15 rejects duplicate keys, trailing prose, non-finite JSON constants,
   output-entry replacement, parent replacement, wrong digest/HEAD/UUID, and
   inconsistent approval states.
5. All four unused v3 key directories were never signed, were archived without
   deletion or overwrite, and cannot be reused at the live paths. Confirm the
   v3 ledger chains to the frozen successful v2 archive ledger.
6. Independently recompute the 23-role digest with real newline delimiters.
7. The four new v4 public-key hashes are correctly bound by all live
   consumers; v3 hashes are absent from every protected source and live
   non-historical supplemental consumer.
8. The rotation changed exactly the intended trust-anchor consumers and did
   not alter M1/M2/G1/G2/R1/B1/C1 definitions, denominators, thresholds,
   multiple testing, four-state compatibility, matched-normal/tumor-REF/CN
   interpretation, raw-read identity, full scope, or claim ceiling.
9. Every formal command remains direct-CLI, source-attested, non-overwriting,
   full 7-dataset/469,849-site scope, and bound to the same authority ID.
10. The 586-test JUnit has matching aggregate attributes and testcase elements,
    and the post-anchor v8 targeted suite is 12/12 passing.

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
