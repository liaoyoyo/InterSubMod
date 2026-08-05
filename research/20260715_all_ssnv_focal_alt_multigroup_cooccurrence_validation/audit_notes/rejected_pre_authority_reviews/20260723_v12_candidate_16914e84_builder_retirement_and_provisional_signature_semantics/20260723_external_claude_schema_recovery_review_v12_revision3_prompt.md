<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v12 revision 3 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v11 R-to-C full-stat 修復、key-lifecycle 修復、fresh reviewer identity binding、append-only failure chain
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v12.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v12 revision 3. This review is READ-ONLY. Never edit, write, create,
move, chmod, or delete files. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v12.py`: size `132555`, SHA-256
   `75c0acada41318ca9f68be8085f49db7561ced6a760a52f749826158966e4cd5`
2. `verify_tumor_ref_receipt_promotion_recovery_v12.py`: size `127669`, SHA-256
   `5a1efae77fb2c2ac591f71d71e12dd3ed6be0d54bb55fcf7a38e63135214dcbe`
3. `replay_m2v5_runner_only_gates_recovery_v12.py`: size `132161`, SHA-256
   `0b18fc0cc41ef80291be6f832d44c8aa78649bee3376ccf9f4cd20c3e6720cc2`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v12.py`: size `296497`, SHA-256
   `6e5e09c0c1a40da75d1b232cd280ccac3ad378fc89d221eab081177db6f0295c`
5. `probe_tumor_ref_schema_recovery_sources_v12.py`: size `45305`, SHA-256
   `e317ca34caf41270b978554e62f17620347205350cb49739cde71711a768f96f`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v12.py`: size `76764`, SHA-256
   `bf5564825d7abdb45b38845f30eac35ea07e3eef96d1e84f73513dbe9af26cd3`
7. `build_tumor_ref_schema_recovery_authority_v12.py`: size `50342`, SHA-256
   `94e3a36bf48fdd3458f36508455b38467c82ff370ff6c6892e034e7d00ac064a`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v12.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`,
`no_output_writes_scope=protected_recovery_and_downstream_namespaces_only`,
`252 passed`, 160 forbidden slots absent, seven reviewed sources,
`review_evidence_state=all_absent`, and cryptographic validation of archived
signed v9, v10, and v11 failures.

Decisive revision-3 review points:

- v11 R emitted four-key records while C required exact nine-key full-stat
  records. v12 R must emit source/canonical full-stat records from retained
  descriptors; C must continue rejecting missing, extra, wrong-type, and
  wrong-value records.
- The unsigned v12 source-set `fee25cdb...` was rejected before authority after
  an internal adversarial review found that signature bytes could be produced
  before persistent-key retirement. Its prior external approval is archived at
  `audit_notes/rejected_pre_authority_reviews/20260723_v12_candidate_fee25cdb_retirement_window/`
  and is not valid for this source set.
- The unsigned source-set `c04ecd89...` was superseded before authority because
  the builder still pinned stale reviewer IDs. Its interrupted external attempt
  and explicit no-verdict internal stops are archived at
  `audit_notes/rejected_pre_authority_reviews/20260723_v12_candidate_c04ecd89_stale_reviewer_ids/`.
- Revision 3 pins the actual fresh agent IDs: Mendel
  `019f8cf9-af38-7d71-8f02-b1eb09742504` and Nash
  `019f8cf9-b625-7762-a018-b4cb0000bbd7`. Confirm these exact IDs in the
  frozen builder and reject any stale or self-reported substitute.
- Revision 2 must stage verified private-key bytes in a mode-0400, nlink-0,
  write-sealed memfd; retire and verify the persistent key at mode 0000 before
  invoking OpenSSL; sign only from the inherited memfd; close the memfd before
  any signature publication; and leave no signature file on any failure.
- Builder `required_module_paths()` must explicitly FD-lease all v11 source,
  authority/review/V/R evidence and the unused continuation key state.
- `no_output_writes` refers only to the enumerated protected recovery and
  downstream namespaces, not the entire filesystem.
- v12 changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated FILTER=PASS sSNVs,
  scientific inputs, BAM files, downstream gates, or claim ceiling.
- Inspect retained-FD/no-reopen behavior, exact 160-slot cardinality,
  no-replace publication, key lifecycle, BaseException cleanup,
  signatures/commit binding, process/supervisor checks, and circular-trust or
  partial-state paths.

Canonical bindings:

- `reviewed_source_set_sha256`: `16914e8414d064d23abc657e65c7afc1c12d7db5d18eba67fba0d1c64e97d555`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `prior_failed_signed_recovery_sha256`: `5cafe2893c0d7d414bd859f10ca29afa317317111f1f90d004c1b2044a05e92a`
- `rejected_generations_sha256`: `50584babd06bc051154585657d1b238573651e8774e6270e879af794aca385c2`
- `scope_sha256`: `1141ff26b6ef6673ba2cac9196b2157e1104f0d0606727124cad738567a49b90`
- trusted public key SHA-256:
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
