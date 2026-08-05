<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v11 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v10 C-stage historical incident projection mismatch 修正、producer-consumer parity、v9/v10 append-only failure evidence
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v11.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v11. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Use only Read, Grep, Glob, and read-only Bash commands.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v11.py`: size `108096`, SHA-256
   `697d028f87a5cbf1e5ef4f14e6add82532807337f4005b0913e9ad3768ef559d`
2. `verify_tumor_ref_receipt_promotion_recovery_v11.py`: size `127669`, SHA-256
   `7a6d75419f9cdfb008e59c3ef5aef47c121c63d7519d8966bb37c8fb5ca5ace0`
3. `replay_m2v5_runner_only_gates_recovery_v11.py`: size `131768`, SHA-256
   `6e0d67ce10c6fa696eed74b8d71fd6cd72b04a434da9c9086a14a5554a01587b`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v11.py`: size `293123`, SHA-256
   `5e61c8b8de84122da367d6e9eb0882ed6b6e6f42e081b6cc97b851304777670d`
5. `probe_tumor_ref_schema_recovery_sources_v11.py`: size `40451`, SHA-256
   `252322bb7427972dd6fc54dd121c609effdba072148e00369fc5181ba9eddcd0`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v11.py`: size `70811`, SHA-256
   `876dc8b758af8ee5ca371cecdf4369ef43f530559a7d8adedab1fa876f799efb`
7. `build_tumor_ref_schema_recovery_authority_v11.py`: size `49987`, SHA-256
   `4eddcc4ec2dec94cf60e266f67dca3a21036c9258eb846dda2c3f48498f8cc03`

Run this exact probe and confirm zero output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v11.py
```

Required probe result: exit `0`, `pass=true`, `no_output_writes=true`, `232
passed`, 147 forbidden output slots absent, seven reviewed sources,
`review_evidence_state=all_absent`, and cryptographic validation of the archived
signed v9 and v10 failures.

Context and decisive v11 delta:

- Signed v10 authority, V receipt, and R receipt were valid. C failed before its
  first downstream write because it compared
  `historical_incident_disclosure.archive_receipt` with the complete archived
  failure object instead of the nested `failed_attempt_archive["receipt"]`
  projection actually emitted by V.
- All v10 formal outputs and reviews were archived append-only. The v10 private
  key is mode 0000. No strict, matched-normal, CN/CCF, final-dataset, final-report,
  C receipt, incident, or terminal receipt was created.
- v11 must accept only the exact nested receipt projection. It must reject the
  former whole-object expectation, a wrong nested receipt, and a receipt-hash
  substitution.
- Regression coverage must use the real archived signed v10 authority, V receipt,
  R receipt, replay log, failure evidence, three reviews, signature, commit, and
  public key. The tests must prove producer-consumer parity across V, R, and C.
- v11 preserves the canonical release key root
  `20260722_all_ssnv_v10_strict_command_parity_bootstrap` and canonical scientific
  screen `all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full`.
  Recovery generation changes must not rename those scientific authorities.
- V, R, and C write only distinct `.recovery.v11` slots. Reviews are pre-existing
  exact mode-0444 APPROVE JSONs; the builder cannot synthesize or normalize them.
- The recovery changes only provenance/control flow. It must not change the seven
  datasets, 469,849 latest same-run LongPhase-S PASS sSNV scope, scientific
  payloads, downstream gates, or the claim ceiling of read-level residual
  epigenetic partition.
- Inspect ordering, retained-FD/no-reopen behavior, exact slot cardinality,
  fail-closed no-replace publication, key retirement, signatures, commit binding,
  process/supervisor checks, v9/v10 failure disclosure, and any circular-trust or
  partial-state path.

Canonical review bindings:

- `reviewed_source_set_sha256`: `b7e4963db70e0b1be3624818f152597e13febd3f0e592d85080fc7470b61487a`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `50584babd06bc051154585657d1b238573651e8774e6270e879af794aca385c2`
- `scope_sha256`: `b61ad4399b89fa7f258f2a6467166fd54b6f1d05134c7fa24b17b299ce556427`
- trusted public key SHA-256:
  `2c1e3b609cb0f58dd23233d257ef66f3950ee5953dcda18474ce209525180b7b`

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
