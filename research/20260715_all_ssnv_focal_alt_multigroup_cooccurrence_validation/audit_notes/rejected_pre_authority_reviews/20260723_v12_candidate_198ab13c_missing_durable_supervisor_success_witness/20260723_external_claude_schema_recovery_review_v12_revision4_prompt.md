<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v12 revision 4 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: authority/continuation one-time-key lifecycle、provisional signature、v11 full-stat 修復與 append-only failure chain
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v12.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v12 revision 4. This review is READ-ONLY. Never edit, write, create,
move, chmod, or delete files. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v12.py`: size `133355`, SHA-256
   `56f6c9bae7685fb6bb6f479061e3c18088e12c973646661fce62a21ef7bf15a6`
2. `verify_tumor_ref_receipt_promotion_recovery_v12.py`: size `127669`, SHA-256
   `b6ef6752569f53e53ce602edf122ac6e87dba1b38f75b0ef497b9812f6d0fa74`
3. `replay_m2v5_runner_only_gates_recovery_v12.py`: size `132161`, SHA-256
   `5a58730c5fe1d0ccc9e08034ba41c8d8dc6d42245644a04380866e3e2f37bbed`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v12.py`: size `297640`, SHA-256
   `7d71c2702750c1af78a4c2485dfacab4ba6323a6f21ddb066fd2a73be4f2f410`
5. `probe_tumor_ref_schema_recovery_sources_v12.py`: size `45305`, SHA-256
   `3d61cf1b6a0ce361c042a38deb7e44043d1b584f7a31c12bb8796fcd6e567ccb`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v12.py`: size `79361`, SHA-256
   `4d737411f722f64221bbc05b62245cdb4277da67a201085a30b1aed5f3bea0d4`
7. `build_tumor_ref_schema_recovery_authority_v12.py`: size `51664`, SHA-256
   `126f8642998c91cbe8126fe92d365e243b44856523d33a1df60b49a46a94234a`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v12.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`,
`no_output_writes_scope=protected_recovery_and_downstream_namespaces_only`,
`252 passed`, 160 forbidden slots absent, seven reviewed sources,
`review_evidence_state=all_absent`, and cryptographic validation of archived
signed v9, v10, and v11 failures.

Decisive revision-4 review points:

- Revision-3 source-set `16914e84...` was rejected before authority. Mendel
  found that the ceremony builder created signature bytes while the persistent
  key was still mode 0400, that the review prompt incorrectly promised no
  signature file after every post-link failure, and that the lifecycle test was
  only a source-order check. Its stopped external attempt and rejection evidence
  are archived under
  `audit_notes/rejected_pre_authority_reviews/20260723_v12_candidate_16914e84_builder_retirement_and_provisional_signature_semantics/`.
- The authority builder must now copy verified key bytes to a mode-0400,
  nlink-0, four-seal memfd, recheck the persistent key, arm BaseException
  retirement, retire and terminally verify the persistent key at mode 000,
  and only then invoke OpenSSL. It must close the memfd before durable signature
  staging. Confirm that no persistent-mode-0400 execution can create either
  ephemeral or durable signature bytes.
- The continuation follows the same stage-retire-sign-close-publish sequence.
  A signature already linked with no-replace is explicitly provisional until
  independent `--verify-signed-terminal` succeeds. Post-link steps are
  observational and mutate no signed payload. A caught lease-holding supervisor
  failure records a blocking continuation incident; any incident prevents
  release authority. Do not impose the impossible old requirement that every
  failure after an already successful atomic link removes or never leaves the
  file. Instead, test whether independent verification plus the incident/live
  state gates can ever grant false authority.
- `CONTINUATION_POLICY` is historical signed content and must remain byte-exact.
  New recovery semantics belong only to `TERMINAL_GOVERNANCE`,
  `EXIT_ATTESTATION_CHECKS`, pass semantics, and source-bound authority checks.
- Regression tests now generate real Ed25519 keys and execute the production
  sealed-memfd staging, persistent-key mode-000 retirement, OpenSSL 64-byte
  signing, and public-key verification paths for both builder and continuation.
- v11 R emitted four-key records while C required nine-key full-stat records.
  v12 R must emit retained-descriptor full-stat records and C must reject
  missing, extra, wrong-type, and wrong-value records.
- Confirm exact reviewer IDs in the builder: Mendel
  `019f8cf9-af38-7d71-8f02-b1eb09742504`, Nash
  `019f8cf9-b625-7762-a018-b4cb0000bbd7`, and this external session ID.
- Inspect interruption and partial-retirement behavior, retained-FD/no-reopen
  behavior, exact 160-slot cardinality, no-replace bundle publication,
  signatures/commit binding, process/supervisor checks, circular trust, and
  partial-state paths.
- v12 changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated FILTER=PASS sSNVs,
  scientific inputs, BAMs, downstream gates, or claim ceiling.

Canonical bindings:

- `reviewed_source_set_sha256`: `198ab13c4bf3cda52263e5eb168f8ee5a95f1f2b8c9b2c42297a5d83c421dfb4`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `prior_failed_signed_recovery_sha256`: `5cafe2893c0d7d414bd859f10ca29afa317317111f1f90d004c1b2044a05e92a`
- `rejected_generations_sha256`: `50584babd06bc051154585657d1b238573651e8774e6270e879af794aca385c2`
- `scope_sha256`: `5904846f7d4ececa67e4d93cfa5ba1ecc7029ca97df95426091fd088ae812f80`
- trusted v12 public-key SHA-256:
  `b6c0734543c9608f8af830c89bdde071a36b3cc38ab9d5c53399fd63a2d46d3d`

Verdict rule: `APPROVE` only if the probe passes without protected-namespace
writes and high findings `[]`, medium findings `[]`, and unresolved conditions
`[]`; otherwise `REJECT`. Low findings may be recorded only when demonstrably
non-blocking.

Return only one JSON object with these exact keys:
`reviewer`, `reviewer_agent_id`, `verdict`, `reviewed_source_set_sha256`,
`probe_exit_code`, `probe_no_output_writes`, `probe_regression_summary`,
`probe_forbidden_output_slots_checked`, `high_findings`, `medium_findings`,
`low_findings`, `unresolved_conditions`, `summary`, and `pass`.

Use reviewer `External Claude Opus` and reviewer/session ID
`7724f023-2491-47ed-9a1e-7856ae194add`.
