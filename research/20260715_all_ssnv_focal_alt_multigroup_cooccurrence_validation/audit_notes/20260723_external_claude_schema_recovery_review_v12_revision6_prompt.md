<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v12 revision 6 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: one-time-key lifecycle、provisional signature、no-return success-witness commit、writer-lease release、v11 full-stat 修復
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v12.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v12 revision 6. This review is READ-ONLY. Never edit, write, create,
move, chmod, or delete files. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v12.py`: size `134395`, SHA-256
   `6096b5823144b65af926e5d438f5480236e4bbc1eccec708bcfa6f92e5dddcc9`
2. `verify_tumor_ref_receipt_promotion_recovery_v12.py`: size `127760`, SHA-256
   `42fc979b01f7e8cbe68f2f11321f6a384d3ef769a2e8d4bdf460ba7cf847bde2`
3. `replay_m2v5_runner_only_gates_recovery_v12.py`: size `132252`, SHA-256
   `10ca48e72f76cb445e67298860ab290ab66a2521db6587e34a04ca238b9b6b65`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v12.py`: size `319106`, SHA-256
   `ab423dced07620b667745122f911a6c90f83c0a672f0cf945cefc6ce7a265300`
5. `probe_tumor_ref_schema_recovery_sources_v12.py`: size `45396`, SHA-256
   `7b4853700b27303fb45da1a26dd5bf16b867050393902a536011a61ada440ca7`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v12.py`: size `83583`,
   SHA-256 `508530fe569ce0aa48dcfefd45221f2f9bad81f476f331d78deff5f11f6622e7`
7. `build_tumor_ref_schema_recovery_authority_v12.py`: size `51664`, SHA-256
   `b29bff632e360d88dee5496ce488dda0bc6fa2e39c61b27231094f2504c593d9`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v12.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`,
`no_output_writes_scope=protected_recovery_and_downstream_namespaces_only`,
`255 passed`, 161 forbidden slots absent, seven reviewed sources,
`review_evidence_state=all_absent`, and cryptographic validation of archived
signed v9, v10, and v11 failures.

Decisive revision-6 review points:

- Revision-4 source-set `198ab13c...` was rejected because a provisional
  signature plus failed incident publication lacked a durable positive
  supervisor-success witness. Revision-5 source-set `33d890a9...` added a
  witness but was rejected because validation, directory fsync, and postcommit
  work remained fallible after the witness became visible. Inspect both archived
  rejection records under `audit_notes/rejected_pre_authority_reviews/`.
- The current `commit_terminal_success_witness` must stage the complete witness
  in an unlinked O_TMPFILE, write/chmod/fsync and verify it while `nlink=0`, run
  all input/source/incident checks, fsync the parent before publication, block
  all catchable signals, rerun terminal checks and exact parent/stage tokens,
  then execute `link_fd_no_replace`. After successful link, the only operation
  must be `os._exit(0)`; there must be no validation, fsync, cleanup, callback,
  logging, return, or exception-producing Python work.
- Link success is the final no-return authority commit point. A SIGKILL before
  link leaves no witness. A SIGKILL after link occurs after commit. If a host
  crash loses the directory entry because no post-link directory fsync occurs,
  the later witness-absence check fails closed. If the atomic link survives,
  its file data was already fsynced and must validate exactly.
- Final `verify_signed_terminal(require_success_witness=True)` must first acquire
  a shared nonblocking flock on the exact v12 recovery public key. The supervisor
  holds the exclusive lock for its lifetime, so final verification cannot grant
  authority while the writer is still alive; the retained shared lease excludes
  any new cooperating writer through verifier exit. Prewitness verification
  intentionally does not acquire this lease because it executes while the
  supervisor owns the exclusive lock.
- The final verifier must still validate the witness schema and exact bindings to
  active source/self-binding, execution receipt, signed exit attestation,
  full-stat signature, public/retired private key, supervisor capability and
  child wait, complete child output inventory, and prewitness stdout/stderr.
  Incident absence and all live artifact checks remain mandatory.
- New tests must run the production O_TMPFILE/link/no-return path in a real fork,
  prove writer/verifier flock mutual exclusion, and invoke the final verifier
  entrypoint in the witness-absent state. Existing tests must retain real
  Ed25519/OpenSSL signing, key retirement, incident-publication failure, malformed
  witness/full-stat records, and recursive source self-binding.
- Authority and continuation signing must stage verified private-key bytes in a
  mode-0400, nlink-0, four-seal memfd; arm BaseException retirement; retire and
  verify the persistent key at mode 000 before OpenSSL; and close the memfd before
  durable signature publication.
- v11 R emitted four-key records while C required nine-key full-stat records.
  v12 R must emit retained-descriptor full-stat records and C must reject missing,
  extra, wrong-type, and wrong-value records.
- Confirm exact reviewer IDs in the builder: Mendel
  `019f8cf9-af38-7d71-8f02-b1eb09742504`, Nash
  `019f8cf9-b625-7762-a018-b4cb0000bbd7`, and this external session ID.
- Same-UID malicious mutation is explicitly outside the signed trusted-account
  threat model; do not silently broaden that model. Within its stated boundary,
  find any false-authority state involving crashes, interruptions, lease races,
  missing incident publication, partial key retirement, or mutable path reopen.
- v12 changes provenance/control flow only. It must not alter the seven datasets,
  469,849 same-run LongPhase-S recalibrated FILTER=PASS sSNVs, scientific inputs,
  BAMs, downstream gates, or biological claim ceiling.

Canonical bindings:

- `reviewed_source_set_sha256`: `e738d06ad0730854a0d72ac29ad9f7dd1d45abba1d4c5ae7aeda3465d791ce93`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `50584babd06bc051154585657d1b238573651e8774e6270e879af794aca385c2`
- `scope_sha256`: `56d03544995cd01919e94761816e76b141ff6714e8d51fe2b86e80af71770f5a`
- trusted v12 public-key SHA-256:
  `b6c0734543c9608f8af830c89bdde071a36b3cc38ab9d5c53399fd63a2d46d3d`

Verdict rule: `APPROVE` only if the exact probe passes without protected writes
and high findings `[]`, medium findings `[]`, and unresolved conditions `[]`;
otherwise `REJECT`. Low findings may be recorded only when demonstrably
nonblocking. Inspect the exact source and tests rather than trusting this prompt.

Return only one JSON object with these exact keys:
`reviewer`, `reviewer_agent_id`, `verdict`, `reviewed_source_set_sha256`,
`probe_exit_code`, `probe_no_output_writes`, `probe_regression_summary`,
`probe_forbidden_output_slots_checked`, `no_return_witness_commit_closes_prior_blocker`,
`writer_lease_release_verified`, `high_findings`, `medium_findings`,
`low_findings`, `unresolved_conditions`, `summary`, and `pass`.

Use reviewer `External Claude Opus` and reviewer/session ID
`7724f023-2491-47ed-9a1e-7856ae194add`.
