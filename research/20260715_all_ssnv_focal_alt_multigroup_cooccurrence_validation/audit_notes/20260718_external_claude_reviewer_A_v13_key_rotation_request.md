<!--
建立時間: 2026-07-18
目標: 外部 Reviewer A 對 source signer key rotation 後的完整 23-role source set 做全新唯讀審查
處理範圍: source identity、key rotation、canonical tests、authority pre-signature state、release command contracts
關聯檔案: ../logs/pytest_full_pre_authority_v17_key_rotation_v13_normalizer_canonical.xml
-->

# External Reviewer A v13: source key rotation and full source set

Perform a fresh, independent, read-only P0 review. Do not reuse the v12 verdict
or evidence. Do not edit, chmod, sign, run producers, run pytest, install
packages, or use the network. Bounded reads, stat, hashing, Git inspection, and
non-mutating shell checks are allowed.

Review all 23 live `EXPECTED_SOURCE_PATHS` from:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`

At minimum, inspect the source-authority validator, both release runners, all
leaf producers/consumers/finalizers, M2 gate and independent audit, matched
normal and CN/CCF paths, raw-read identity checks, report builder, portable
delivery/QA, and relevant tests.

Required anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Exactly 23 protected roles, all mode `0444`
- Source-set SHA-256:
  `8dac1eaab5f624568b11a36eeb8b8698d371b253a368194f22f9ca73ba5cb99b`
- Canonical JUnit:
  `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v17_key_rotation_v13_normalizer_canonical.xml`
- JUnit SHA-256:
  `f11bd1427719e3c24ae0e6792df2d004736d7b3093ec44345edfcc4352f62067`
- JUnit counts: 510 tests, 0 failures, 0 errors, 0 skipped
- Public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v5_summary_hotfix/ed25519_public.pem`
- Public-key SHA-256:
  `e91076b68495d9d879b378f3af8f16c00b3aa7bd3711943bf1dc4bff66d85cbe`
- Public/private key modes before signing: `0444` / `0400`
- Current v5 authority JSON, approval JSON, detached signature, and pending
  signature must all be absent before assembly/signing.

Explicitly verify:

1. Independently recompute the 23-role digest using real newline delimiters.
2. The old public-key SHA
   `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920`
   is absent from every protected source.
3. The failed unsigned v5 attempt is only in
   `docs/archive/2026/07/20260718_unsigned_source_authority_v5_signer_failure_01/`
   and cannot authorize or enter the live chain.
4. `run_matched_normal_candidate_controls.py` defaults `--normal-audit` to
   `CANONICAL_NORMAL_AUDIT`; the superseded
   `matched_normal_methyl_tag_audit.v2.json` token is absent from every
   protected source and covered by the stale-token test.
5. The key rotation and default cleanup did not change M1/M2/G1/G2/R1/B1/C1
   definitions, denominators, thresholds, multiple testing, four-state
   compatibility, normal/tumor-REF/CN interpretation, path versions, output
   scope, raw-read identity, or claim ceiling.
6. Every formal command remains direct-CLI, source-attested, non-overwriting,
   full 7-dataset/469,849-site scope, and bound to the same authority ID.

Return exactly one JSON object with these 15 top-level keys and no others:

`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Contract:

- `schema_name`: `intersubmod.external_claude_source_review`
- `schema_version`: `1.0.0`
- `reviewer_label` must contain `Reviewer A`
- `reviewer_id`: a fresh canonical lowercase UUIDv4
- `model`: starts with `claude-`
- `verdict`: `APPROVE` or `REQUEST_CHANGES`
- Clean approval requires `findings_closed=true`,
  `f1_status=RESOLVED_VERIFIED`, `f2_status=RESOLVED_VERIFIED`, and an empty
  `blocking_findings` array.
- `reviewed_source_set_sha256` and `reviewed_git_head` must equal the anchors.
- `evidence` must be a JSON object, not an array.

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
