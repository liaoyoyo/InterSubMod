<!--
建立時間: 2026-07-19
目標: External Reviewer B adversarial review of v6/v17/v10 independent authority and the full Task-B release chain
處理範圍: false authority, path/FD/content races, review/JUnit binding, stale-key isolation, and scientific semantics
關聯檔案: ../logs/pytest_full_pre_authority_v21_v6_signer_v17_normalizer_v10_assembler_canonical.xml
-->

# External Reviewer B v17: adversarial independent-authority review

Perform a fresh, independent, adversarial, read-only P0 review. Do not rely on
prior verdicts. Do not edit, chmod, sign, execute producers, run pytest,
install packages, or use web/network tools. Bounded reads, stat, hashing, and
Git inspection are allowed.

Authoritative source enumeration:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py::EXPECTED_SOURCE_PATHS`

Review all 23 sources plus the v6 signer, v17 normalizer, v10 assembler, their
tests, `tests/test_release_source_authority_fd_signature.py`,
`audit_notes/archive_unused_v5_signer_keys.py`, and
`audit_notes/20260719_unused_v5_signer_key_archive.v1.jsonl`.

Fixed anchors:

- Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID: `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Protected roles: 23, all mode `0444`
- Source-set SHA-256:
  `6abaa35ee6045d9cb483623d1a8d08628100a162b3d1b64b91d0ba869963af77`
- Release validator SHA-256:
  `eeb0d048bae3733d7cb7c44cbf0d669e8d471653f3cfc180074f7872c876105d`
- Signer/normalizer/assembler SHA-256 values:
  `e032a5c56dcb40b81b4c963413d7190a97894c7c5258c67f29700bc64503dd49`,
  `14fee1f38683490f8389beb3c66e726433a5eaba8e8a48468114fc27099accd8`,
  `424caff771c6c7fa5695e6e96f775baca5e6f835a684425490216f3093c67336`
- Assembler path/size/mode:
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/assemble_release_source_authority_v10.py`,
  `28911`, `0o444`
- Canonical JUnit path:
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v21_v6_signer_v17_normalizer_v10_assembler_canonical.xml`
- JUnit size/SHA/mode/counts:
  `159743`,
  `4cd9746825f3d50fbb08ae907ced6a3c7b3204a4a2c52812a9a30a81b42f53bc`,
  `0o444`, `722/0/0/0`
- Source/result/report/supplemental public-key SHA-256 values:
  `415f1d4d05e09cbd6168f06c714d2ac825dd5449f532a479927b17f142526b7b`,
  `8f9a1533bd13c123e7e4024ee379cc39873c5e9d6b3344028de6e291ab5d8023`,
  `e511b5b8f5d05f2515ecefb8d8685e08690057e0a4856f022d95dc740b288273`,
  `d3bad977ab5fa737b1d1088ec5ff8400ebfe99ea438b40487c970743d2823b24`
- All four live key directory/public/private modes: `0700/0444/0400`
- v5 archive script/ledger SHA-256:
  `4397c303c40e091bc6e63423bf64e846b97ab26b156b6f106b0ef02c10b14674`,
  `58aabdd0417159b42a01bda2286fd84c8fc21d839cecbd460e94717169b886f8`
- Live authority, approval, signature, pending signature, and formal downstream
  targets are absent

Adversarial checks:

1. Inject a canonical signature or target replacement during signer descriptor
   teardown after private mode becomes `000`. Confirm the signer may only
   report a provisional ceremony event, never `SIGNED` or authority.
2. Search every caller for any path that treats signer exit `0`, signature
   existence, or point-in-time verification as release approval without a
   post-process production consumer reopen.
3. Attempt a self-consistent A/B forgery: both reviews agree on a wrong
   assembler or wrong JUnit. Confirm v10 and production validator compare them
   against live bound artifacts and fail.
4. Attempt split-view races between review artifact identity and parsed bytes,
   assembler self identity, JUnit artifact/counts, approval summary, and signed
   authority evidence. Confirm bound descriptors and late checks close each
   path.
5. Check that assembler v10 does not embed its own hash and does not embed a
   pre-freeze JUnit hash. Confirm assembler bytes were frozen before v21 XML,
   while both reviews must bind exact identities.
6. Replace/chmod/delete either normalized review, v10 assembler, v21 XML,
   source signature, public key, private key mode, any protected source, or Git
   HEAD after assembly/signing. Confirm production validation fails.
7. Verify v5 key archive has four never-signed no-replace completions; no old
   key path/hash remains accepted by the six live protected consumers.
8. Recompute the 23-role digest using real newline delimiters. Trace all
   source/result/report/supplemental key paths through both formal runners,
   finalizers, report builder, and supplemental seal.
9. Trace option-to-value/path bindings through full downstream execution:
   primary-pre v4, preflight v8, cooccurrence v7, independent M2 v3, controls
   v3, CN/CCF v3, and final dataset/report v5. Check no overwrite or subset
   fallback.
10. Confirm trust-chain edits did not alter M1/M2/G1/G2/R1/B1/C1,
    denominators, four-state semantics, negative controls, raw-read identity,
    same-run LongPhase-S HP/PS, 7-dataset/469,849-site scope, or the claim
    ceiling forbidding confirmed subclone/lineage claims.
11. Independently parse JUnit attributes/elements as 722/0/0/0. Check focused
    runs: signer 22, normalizer 26, assembler 15, consumer/key suite 104.

Return exactly one JSON object with these 17 top-level keys and no others:

`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `reviewed_assembler`,
`canonical_junit`, `review_scope`, `blocking_findings`,
`nonblocking_findings`, `evidence`.

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
- Source digest, Git HEAD, `reviewed_assembler`, and `canonical_junit` exactly
  equal the fixed anchors above
- `canonical_junit` has only `artifact` and `counts`
- `evidence` is a JSON object

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
