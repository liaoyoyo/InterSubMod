<!--
建立時間: 2026-07-18
目標: 外部 Reviewer B 對 key rotation 後 Task-B release chain 做全新 adversarial 唯讀審查
處理範圍: producer-consumer command/path binding、科學語義、signing trust anchor、stale-source isolation
關聯檔案: ../logs/pytest_full_pre_authority_v17_key_rotation_v13_normalizer_canonical.xml
-->

# External Reviewer B v13: adversarial key-rotation review

Perform a fresh, independent, adversarial, read-only P0 review. Do not rely on
v12 findings. Do not edit, chmod, sign, execute producers, run pytest, install
packages, or use the network. Bounded reads, stat, hashing, Git inspection, and
non-mutating shell checks are allowed.

Authoritative source enumeration:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py::EXPECTED_SOURCE_PATHS`

Fixed anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Protected roles: 23, all mode `0444`
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
- Pre-signing key modes: public `0444`, private `0400`
- Live v5 authority, approval, signature, and pending signature are absent.

Adversarial checks:

1. Independently recompute the exact 23-role digest with real newlines.
2. Trace every option-to-value binding across the cooccurrence runner/finalizer,
   completion runner, matched-normal runner/analyzer, CN/CCF annotator, final
   dataset builder, result/report finalizer, and portable report QA.
3. Confirm current versions remain primary-pre v4, preflight v8,
   cooccurrence v7, independent M2 v3, matched controls/analysis v3, CN/CCF v3,
   final dataset/report v5, and that fresh output guards remain fail-closed.
4. Confirm the source key rotation changed only the trust-anchor hash in the
   validator and two runners. The old key hash must be absent from all protected
   roles; the archived unsigned attempt must be non-live and non-authorizing.
5. Confirm the matched-normal default now uses `CANONICAL_NORMAL_AUDIT`, and
   the stale-token test rejects `matched_normal_methyl_tag_audit.v2.json`.
6. Confirm no source/key change altered M2 pre-axis semantics, K>10 treatment,
   row conservation, formal G1/G2/B1 denominators, four-state compatibility,
   negative-control interpretation, truth-blind selection, raw-read identity,
   HP/PS handling, full-scope counts, or claim ceiling.
7. Confirm old v12 reviews bind a different source digest and therefore cannot
   authorize this source set.

Return exactly one JSON object with these 15 top-level keys and no others:

`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Contract:

- `schema_name`: `intersubmod.external_claude_source_review`
- `schema_version`: `1.0.0`
- `reviewer_label` must contain `Reviewer B`
- `reviewer_id`: a fresh canonical lowercase UUIDv4
- `model`: starts with `claude-`
- `verdict`: `APPROVE` or `REQUEST_CHANGES`
- Clean approval requires `findings_closed=true`,
  `f1_status=RESOLVED_VERIFIED`, `f2_status=RESOLVED_VERIFIED`, and no blocking
  findings.
- `reviewed_source_set_sha256` and `reviewed_git_head` must equal the anchors.
- `evidence` must be a JSON object, not an array.

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
