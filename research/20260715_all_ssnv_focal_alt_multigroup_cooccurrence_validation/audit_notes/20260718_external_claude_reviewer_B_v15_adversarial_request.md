<!--
建立時間: 2026-07-18
目標: 外部 Reviewer B 對 v4/v8/v15 closure 與 Task-B release chain 做 fresh adversarial 唯讀審查
處理範圍: digest/race/failure injection、live-review binding、producer-consumer binding、科學語義、stale-source/key isolation
關聯檔案: ../logs/pytest_full_pre_authority_v19_v4_signer_v15_normalizer_v8_assembler_canonical.xml
-->

# External Reviewer B v15: adversarial ceremony and release review

Perform a fresh, independent, adversarial, read-only P0 review. Do not rely on
v14 findings. Do not edit, chmod, sign, execute producers, run pytest, install
packages, or use the network. Bounded reads, stat, hashing, Git inspection, and
non-mutating shell checks are allowed.

Authoritative source enumeration:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py::EXPECTED_SOURCE_PATHS`

Review signer v4, assembler v8, normalizer v15 and all corresponding tests in
the same paths listed by Reviewer A.

Fixed anchors:

- Git HEAD:
  `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID:
  `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Protected roles: 23, all mode `0444`
- Source-set SHA-256:
  `b5e432fbd98e657ae48d6296fe3050eb01b3ccbc57ecc4aa0e5a25e37600ac9e`
- Canonical JUnit:
  `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v19_v4_signer_v15_normalizer_v8_assembler_canonical.xml`
- JUnit SHA-256:
  `094e64cd7c7c45e0e45c31f2485a8934ca2ebee9e2c6da92cee8f5afcfaf2179`
- JUnit counts: 586 tests, 0 failures, 0 errors, 0 skipped
- v4 signer source:
  `6698393db61acf3c8fc242ac6e011eca6f96a128adf23f3d83664c95db654595`
- v8 assembler source:
  `e1a43d7e30febd041759901aefd6bd68c5f7d4cbed25673606ee74cbf8ee2904`
- v15 normalizer source:
  `840534e11f86fe04d22b2440c1830bdee0b1c774b3a6c1a6af9f4a99f160d0c0`
- Source public-key SHA-256:
  `ae8377ff1753e98e8490e44b58d41503d7484bb19704247be222123e10f4801a`
- Source key-directory/public/private modes: `0700/0444/0400`
- Result/report/supplemental public-key SHA-256 values:
  `1c9cea01da7106ac502ce708831e8fd5b19854fe59956016b836301f100f0404`,
  `25face031089f413d564c93db4491a4433c81f97f339a616091beba077eca848`,
  `db20e9ca41cb78201578d01dd724a82d3894574a2ae63b583509dbdec9a1a6e8`
- Result/report/supplemental key-directory/public/private modes:
  `0700/0444/0400`
- v3 archive ledger SHA-256:
  `93bceb16d247f021a52fc6b17083dbbf2ba26cbc37da15d730b1de793a9f57e5`
- Prior v2 archive ledger SHA-256:
  `49860985d5c20c3801f2bb92d7124ea97d079f64f6a500a418244d697b05b94b`
- Live v5 authority, approval, signature and pending signature are absent.

Adversarial checks:

1. Try to make signer v4 sign a same-path replacement after assembler v8 has
   approved different bytes. Confirm the canonical `SIGN` instruction binds
   exact path/size/SHA-256/mode and noncanonical or mismatched records fail.
2. Chmod or replace the key directory after `SIGNER_READY`; prove the retained
   directory FD and mode checks stop signing at every OpenSSL boundary.
3. Inject failure before and after signature publication and before and after
   private-key `fchmod(000)`. Confirm the failure JSON reports actual
   signature publication, retirement-started and private mode rather than
   always claiming a retained key.
4. Try normalizer v15 output-parent replacement, output-entry replacement,
   preexisting output, partial write, chmod/fsync failure and byte mutation.
   Confirm parent/output FDs, `O_EXCL`, path-entry identity and read-back close
   every path/payload split.
5. Replace, delete, chmod or semantically alter a v15 review after approval
   assembly. Confirm the live release validator re-opens both fixed review
   paths, requires mode `0444`, exact signed artifact records, exact schema and
   clean-approval fields, and retains the review FDs through final validation.
6. Verify the v3 archive ledger proves four never-signed no-replace moves,
   chains to the successful v2 archive ledger, and all live v4 key dirs have
   modes `0700/0444/0400`.
7. Independently recompute the exact 23-role digest using real newlines.
8. Trace every option-to-value/path binding across both formal runners,
   finalizers, matched-normal, CN/CCF, report delivery and QA. Confirm current
   versions remain primary-pre v4, preflight v8, cooccurrence v7, independent
   M2 v3, matched controls/analysis v3, CN/CCF v3, final dataset/report v5.
9. Confirm new key rotation changed only intended trust anchors and did not
   alter M2 pre-axis semantics, K>10 treatment, row conservation, G1/G2/B1
   denominators, four-state compatibility, negative-control interpretation,
   truth-blind selection, raw-read identity, HP/PS handling, full-scope counts,
   or claim ceiling.
10. Confirm all v14 reviews bind an earlier digest and all v3 public-key hashes
   are absent from protected/live consumers, so neither can authorize v15.
11. Independently parse the canonical JUnit and require attributes/elements
    both equal 586 tests, 0 failures, 0 errors and 0 skipped; separately verify
    post-anchor assembler v8 tests are 12/12 passing.

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
