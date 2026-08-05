<!--
建立時間: 2026-07-19
目標: External Reviewer A fresh review of the v6 two-stage signer, v17 normalizer, v10 assembler, production consumer, and full Task-B source set
處理範圍: source/key/JUnit anchors, independent consumer authority, producer contracts, and scientific claim guards
關聯檔案: ../logs/pytest_full_pre_authority_v21_v6_signer_v17_normalizer_v10_assembler_canonical.xml
-->

# External Reviewer A v17: two-stage authority and full-chain review

Perform a fresh, independent, read-only P0 review. Do not reuse a prior
verdict. Do not edit, chmod, sign, execute producers, run pytest, install
packages, or use web/network tools. Bounded reads, stat, hashing, and Git
inspection are allowed.

Review all 23 live `EXPECTED_SOURCE_PATHS` from:

`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`

Also review:

- `audit_notes/run_one_time_ed25519_signer_v6.py`
- `tests/test_one_time_ed25519_signer_v6.py`
- `audit_notes/normalize_external_source_review_v17.py`
- `tests/test_normalize_external_source_review_v17.py`
- `audit_notes/assemble_release_source_authority_v10.py`
- `tests/test_assemble_release_source_authority_v10.py`
- `tests/test_release_source_authority_fd_signature.py`
- `audit_notes/archive_unused_v5_signer_keys.py`
- `audit_notes/20260719_unused_v5_signer_key_archive.v1.jsonl`

Required anchors:

- Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Authority ID: `20260718_all_ssnv_focal_alt_task_b_release_v5`
- Exactly 23 protected roles, all mode `0444`
- Source-set SHA-256:
  `6abaa35ee6045d9cb483623d1a8d08628100a162b3d1b64b91d0ba869963af77`
- Release validator SHA-256:
  `eeb0d048bae3733d7cb7c44cbf0d669e8d471653f3cfc180074f7872c876105d`
- v6 signer SHA-256:
  `e032a5c56dcb40b81b4c963413d7190a97894c7c5258c67f29700bc64503dd49`
- v17 normalizer SHA-256:
  `14fee1f38683490f8389beb3c66e726433a5eaba8e8a48468114fc27099accd8`
- v10 assembler identity:
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/assemble_release_source_authority_v10.py`,
  size `28911`, SHA-256
  `424caff771c6c7fa5695e6e96f775baca5e6f835a684425490216f3093c67336`,
  mode `0o444`
- Canonical JUnit identity:
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v21_v6_signer_v17_normalizer_v10_assembler_canonical.xml`,
  size `159743`, SHA-256
  `4cd9746825f3d50fbb08ae907ced6a3c7b3204a4a2c52812a9a30a81b42f53bc`,
  mode `0o444`
- JUnit attributes and elements: 722 tests, 0 failures, 0 errors, 0 skipped
- The signer, normalizer, assembler, and their tests are mode `0444`
- Source/result/report/supplemental public-key SHA-256 values:
  `415f1d4d05e09cbd6168f06c714d2ac825dd5449f532a479927b17f142526b7b`,
  `8f9a1533bd13c123e7e4024ee379cc39873c5e9d6b3344028de6e291ab5d8023`,
  `e511b5b8f5d05f2515ecefb8d8685e08690057e0a4856f022d95dc740b288273`,
  `d3bad977ab5fa737b1d1088ec5ff8400ebfe99ea438b40487c970743d2823b24`
- The four live key directory/public/private modes before signing are
  `0700/0444/0400`
- v5 archive script SHA-256:
  `4397c303c40e091bc6e63423bf64e846b97ab26b156b6f106b0ef02c10b14674`
- v5 archive ledger SHA-256:
  `58aabdd0417159b42a01bda2286fd84c8fc21d839cecbd460e94717169b886f8`
- The v5 ledger records four `ARCHIVED_UNUSED_NEVER_SIGNED` completions and
  chains to v4 ledger SHA-256
  `b7720d1ab680f9f87b23f238a2f9295df1f78acf9419a1133fb44f492028a9a9`
- Authority, approval, detached signature, pending signature, and all formal
  downstream targets are absent before assembly/signing

Required closure checks:

1. Signer v6 never emits `SIGNED` or grants authority. Its only successful
   event is `CEREMONY_OUTPUT_AVAILABLE_REQUIRES_INDEPENDENT_VERIFICATION`,
   with `consumer_verification_required=true` and
   `release_authority_granted=false`.
2. Model mutation during descriptor teardown. Confirm signer output remains
   provisional and the independent production consumer must reopen live paths.
3. Production `release_source_authority.py` reopens approval, signature,
   public/private key, all 23 sources, both normalized reviews, v10 assembler,
   and canonical JUnit. Confirm crypto, path, mode, digest, review semantics,
   UUID, and JUnit checks are fail closed.
4. Normalizer v17 requires exactly 17 top-level keys and exact nested
   `reviewed_assembler` and `canonical_junit` records; it rejects duplicate,
   non-finite, trailing, malformed, and semantically inconsistent inputs.
5. Assembler v10 binds its own live bytes without embedding its own SHA,
   parses the post-freeze JUnit, requires A/B records to agree with live
   assembler/JUnit identities, and emits one canonical digest-bound `SIGN`
   instruction. Confirm there is no self-reference or stale-JUnit cycle.
6. Verify v5 key directories were archived no-replace, never signed, and old
   hashes/paths are absent from the six live protected consumers.
7. Independently recompute the 23-role digest with real newline delimiters.
8. Trace both formal runners and all downstream producers. Commands remain
   direct-CLI, source-attested, non-overwriting, and full
   7-dataset/469,849-site scope.
9. Confirm trust-anchor changes did not alter M1/M2/G1/G2/R1/B1/C1
   definitions, denominators, four-state compatibility, matched-normal and
   tumor-REF interpretation, raw-read identity, same-run LongPhase-S HP/PS,
   or the claim ceiling. No candidate may be called a confirmed cellular
   subclone or lineage topology from these data alone.
10. Independently parse JUnit attributes and testcase elements as 722/0/0/0.
    Focused evidence is signer 22/22, normalizer 26/26, assembler 15/15,
    key-rotation/consumer suite 104/104, and full canonical 722/722.

Return exactly one JSON object with these 17 top-level keys and no others:

`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `reviewed_assembler`,
`canonical_junit`, `review_scope`, `blocking_findings`,
`nonblocking_findings`, `evidence`.

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
- `reviewed_source_set_sha256` and `reviewed_git_head` equal the anchors above
- `reviewed_assembler` is the exact four-field identity above
- `canonical_junit` has exactly `artifact` and `counts`, using the exact
  identity/counts above
- `evidence` is a JSON object

Do not add prose, a Markdown fence, or trailing text outside the JSON object.
