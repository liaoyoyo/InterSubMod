<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v12 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v11 R-to-C source full-stat projection mismatch、append-only failure evidence、v12 producer-consumer parity
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v12.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v12. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Use only Read, Grep, Glob, and read-only Bash commands.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v12.py`: size `132112`, SHA-256
   `a2a4b4edb5ab0a1dc9d069099d595e85878c1d98ddc7c8063ecc8876a35210e4`
2. `verify_tumor_ref_receipt_promotion_recovery_v12.py`: size `127669`, SHA-256
   `a5c5e899dc2ee281418d2f6058e8be935fdd963bb0102297ab4ed7920a6735fc`
3. `replay_m2v5_runner_only_gates_recovery_v12.py`: size `132161`, SHA-256
   `8fd71b4806fcb623928fdc4dcedcba6b08045ef967a9262d3d366ad922f6b2db`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v12.py`: size `294082`, SHA-256
   `2f74cd33aa2de7c0a412803637c8b33b41bc034802780cb84429a2c9440828af`
5. `probe_tumor_ref_schema_recovery_sources_v12.py`: size `45170`, SHA-256
   `c84a5d886d0852acba1c82b6339b6b6a498918a53bb710dbde0f78e822879770`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v12.py`: size `74943`, SHA-256
   `c9484c6f8b6597fb6bbe83f66400caace22ce869e883799f3e764512539d9a6b`
7. `build_tumor_ref_schema_recovery_authority_v12.py`: size `49987`, SHA-256
   `5c8cb2cad5839a5a7a38737a565fc5d62d26059ef5dfd519024c6c948c7d31aa`

Run this exact probe and confirm zero output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v12.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `250 passed`,
160 forbidden slots absent, seven reviewed sources,
`review_evidence_state=all_absent`, and cryptographic validation of archived
signed v9, v10, and v11 failures.

Decisive v12 delta:

- Signed v11 authority, V receipt, and R receipt were valid. C failed before its
  first downstream write because R emitted four-key basic records for
  `promotion_trust_chain.historical_source_receipt` and
  `canonical_source_receipt`, while C required nine-key full-stat records.
- v11 authority/reviews/V/R were archived append-only. The authority private key
  is mode 0000. Its unused continuation private key remains mode 0400 and must
  not be represented as retired. No C or scientific downstream output exists.
- v12 must fix the R producer to emit the full-stat records obtained from its
  retained source descriptors. C must continue requiring exact nine-key records;
  lowering C to a basic projection is forbidden.
- Regression tests must read the real archived v11 V/R receipts, reproduce each
  old mismatch independently, prove the repaired pair passes, and reject missing
  or extra fields.
- v12 changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated FILTER=PASS sSNVs,
  scientific inputs, downstream gates, BAM files, or claim ceiling.
- Inspect exact slot cardinality, retained-FD/no-reopen behavior, no-replace
  publication, key lifecycle, signatures/commit binding, process/supervisor
  checks, v9/v10/v11 disclosure, and circular-trust or partial-state paths.

Canonical bindings:

- `reviewed_source_set_sha256`: `fee25cdb575e5b3c82cb45e7e65deb45cc9a936560c14092c12f553878bb30c1`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `prior_failed_signed_recovery_sha256`: `5cafe2893c0d7d414bd859f10ca29afa317317111f1f90d004c1b2044a05e92a`
- `rejected_generations_sha256`: `50584babd06bc051154585657d1b238573651e8774e6270e879af794aca385c2`
- `scope_sha256`: `8e340b13133e83c963f7f636ceb73a1739f3c7723c48a192c12a1a4fe054217f`
- trusted public key SHA-256:
  `b6c0734543c9608f8af830c89bdde071a36b3cc38ab9d5c53399fd63a2d46d3d`

Verdict rule: `APPROVE` only if the probe passes without writes and high findings
`[]`, medium findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`.
Low findings may be recorded only when demonstrably non-blocking.

Return only one JSON object with these exact keys:
`reviewer`, `reviewer_agent_id`, `verdict`, `reviewed_source_set_sha256`,
`probe_exit_code`, `probe_no_output_writes`, `probe_regression_summary`,
`probe_forbidden_output_slots_checked`, `high_findings`, `medium_findings`,
`low_findings`, `unresolved_conditions`, `summary`, and `pass`.

Use reviewer `External Claude Opus` and reviewer/session ID
`7724f023-2491-47ed-9a1e-7856ae194add`.
