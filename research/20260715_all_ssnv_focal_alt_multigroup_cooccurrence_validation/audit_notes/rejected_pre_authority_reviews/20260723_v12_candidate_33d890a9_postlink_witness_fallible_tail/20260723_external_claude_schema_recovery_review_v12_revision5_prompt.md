<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v12 revision 5 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: authority/continuation one-time-key lifecycle、provisional signature、durable supervisor-success witness、v11 full-stat 修復與 append-only failure chain
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v12.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v12 revision 5. This review is READ-ONLY. Never edit, write, create,
move, chmod, or delete files. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v12.py`: size `133914`, SHA-256
   `87a6cfce0124105924f18da5b561720108f78ec0cdecad19a5c9bfe44a2d7936`
2. `verify_tumor_ref_receipt_promotion_recovery_v12.py`: size `127760`, SHA-256
   `c0e790c65238cc1da97a6dcbfc9e945dedfbfc8f66817247bb50c2cb6624f8cf`
3. `replay_m2v5_runner_only_gates_recovery_v12.py`: size `132252`, SHA-256
   `28e21f1c833cfd6062e5bcc8841cd195d464221c54710aa56ec75d84e5799d87`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v12.py`: size `314768`, SHA-256
   `5e1717d1cfff380d434a3c0f3970f963041dd8ef74ca66b557c9eb5d75c23b94`
5. `probe_tumor_ref_schema_recovery_sources_v12.py`: size `45396`, SHA-256
   `ab8de6beba122430054262836cee9b2b119fa392c4efafadbbd3943bd7a4df15`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v12.py`: size `80383`, SHA-256
   `81493fa3518a049193123485200da759a8bc44fdc2a343b54ae033c000297c0b`
7. `build_tumor_ref_schema_recovery_authority_v12.py`: size `51664`, SHA-256
   `8fdf8f21290f5d19fa4961a1f789200bbcf2c46e9b3b3a85c12ebf7ad3252736`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v12.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`,
`no_output_writes_scope=protected_recovery_and_downstream_namespaces_only`,
`253 passed`, 161 forbidden slots absent, seven reviewed sources,
`review_evidence_state=all_absent`, and cryptographic validation of archived
signed v9, v10, and v11 failures.

Decisive revision-5 review points:

- Revision-4 source-set `198ab13c...` was rejected before authority even though
  its external review approved it. Both independent internal reviewers found
  the same false-authority path: after signature publication, a killed or
  failed supervisor plus failed best-effort incident publication could leave a
  valid signature with no durable proof of clean supervisor completion. Review
  the archived evidence at
  `audit_notes/rejected_pre_authority_reviews/20260723_v12_candidate_198ab13c_missing_durable_supervisor_success_witness/rejection_evidence.json`.
- The current continuation adds the no-replace artifact
  `results/m2v5_downstream_continuation_supervisor_success.recovery.v12.json`.
  The signature is explicitly provisional. The read-only
  `--verify-signed-terminal-prewitness` path requires the witness to be absent,
  validates the otherwise complete signed terminal state, and returns
  `release_authority=false`. Only after that preverification and all terminal
  rechecks does the lease-holding supervisor atomically commit the success
  witness.
- The final independent `--verify-signed-terminal` path requires the witness,
  validates its exact schema and bindings to source/self-binding, execution
  receipt, exit attestation, full-stat signature record, public/retired key,
  supervisor capability and child wait evidence, complete artifact inventory,
  and preverification stdout/stderr digests. It may return
  `release_authority=true` only after this check.
- Verify the central invariant directly: a signature-linked but witness-absent
  state must never grant release authority, even when incident publication
  itself fails. Conversely, witness publication must be impossible before
  successful preverification and all terminal checks. Check concurrent verifier
  access, interruption points, no-replace behavior, and witness exactness.
- Regression tests must use real Ed25519/OpenSSL paths and include fault
  injection for `signature present + incident commit failure + witness absent`
  (final gate rejects), `witness present` (final gate can proceed), and
  `prewitness invoked after witness exists` (rejects). Do not accept only a
  source-order assertion.
- The authority and continuation signing paths must stage verified key bytes in
  a mode-0400, nlink-0, four-seal memfd; arm BaseException retirement; retire
  and verify the persistent private key at mode 000 before invoking OpenSSL;
  close the memfd before durable signature publication; and fail closed on
  interruption or partial retirement.
- The builder must use one validator FD, no-replace atomic bundle publication,
  retained stage/member FD validation, and the exact reviewer IDs. The
  continuation must retain source/interpreter FDs and must not re-open mutable
  source paths for authority decisions.
- v11 R emitted four-key records while C required nine-key full-stat records.
  v12 R must emit retained-descriptor full-stat records and C must reject
  missing, extra, wrong-type, and wrong-value records.
- `CONTINUATION_POLICY` is historical signed content and must remain byte-exact.
  New recovery semantics belong only to recovery governance/checks. Confirm the
  current code does not alter scientific payload, seven datasets, 469,849
  same-run LongPhase-S recalibrated FILTER=PASS sSNVs, BAMs, downstream gates,
  or biological claim ceiling.

Canonical bindings:

- `reviewed_source_set_sha256`: `33d890a9b4143c70cc70d1c5259305f00fcc99751ddc35b78506590dc7f287a0`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `50584babd06bc051154585657d1b238573651e8774e6270e879af794aca385c2`
- `scope_sha256`: `f297d7b4f0cb9fd0025e7256ec1329bd5b3745ac3e6b115bde8b890a9a409855`
- trusted v12 public-key SHA-256:
  `b6c0734543c9608f8af830c89bdde071a36b3cc38ab9d5c53399fd63a2d46d3d`

Verdict rule: `APPROVE` only if the exact probe passes without protected writes
and high findings `[]`, medium findings `[]`, and unresolved conditions `[]`;
otherwise `REJECT`. Low findings may be recorded only when demonstrably
non-blocking. Do not rely on this prompt's remediation description; inspect the
exact source and tests.

Return only one JSON object with these exact keys:
`reviewer`, `reviewer_agent_id`, `verdict`, `reviewed_source_set_sha256`,
`probe_exit_code`, `probe_no_output_writes`, `probe_regression_summary`,
`probe_forbidden_output_slots_checked`, `durable_success_witness_closes_prior_blocker`,
`high_findings`, `medium_findings`, `low_findings`, `unresolved_conditions`,
`summary`, and `pass`.

Use reviewer `External Claude Opus` and reviewer/session ID
`7724f023-2491-47ed-9a1e-7856ae194add`.
