<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v15 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v14 metadata-only relation 失敗根因、v15 metadata-plus-size 修復、完整 recovery chain 與 scientific claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v15.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v15. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Use only Read, Grep, Glob, and explicitly allowed read-only
Bash commands. Do not trust this prompt as evidence: independently inspect the
frozen sources and run the exact probe below.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v15.py`: size `198269`, SHA-256 `e09220da97cffcdf8edcdd44f447fac88292eb2b989543ea774a8547f85fe090`
2. `verify_tumor_ref_receipt_promotion_recovery_v15.py`: size `127906`, SHA-256 `aa8c5ed29a46ba4a4618e9b40bf7d627df7cf15ffc6daea240ad16aafe3bd043`
3. `replay_m2v5_runner_only_gates_recovery_v15.py`: size `150445`, SHA-256 `0f0811a9171653eea2f59e81b1e5283f281457f40b6c64960c9f32f34d52a4d9`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v15.py`: size `344261`, SHA-256 `df5e93539b7fe33d4b8ccc992656ed5688f265e1287d84176564231c30db989f`
5. `probe_tumor_ref_schema_recovery_sources_v15.py`: size `62434`, SHA-256 `737b1d31fbbe59c589c5cda275d5d219ce5a76d5163f38b7afef97093ce03bdf`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v15.py`: size `113042`, SHA-256 `9e379b4f1ba2f9d27228b9bba5c8ea03ad9f28d6ca94c8a93d8774a0db754bbf`
7. `build_tumor_ref_schema_recovery_authority_v15.py`: size `53376`, SHA-256 `bfc4de402d3f1d6fac836d179ed3a080b8008c3874f1f48b634bc4b621778582`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v15.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `350 passed`,
206 forbidden slots absent, 14 staging patterns absent, seven current frozen
sources, and the append-only v9/v10/v11/v12/v14 signed-failure chain plus v13
unsigned rejection verified. Before formal review publication the review evidence
state is expected to be `all_absent`.

Audit these decisive points rather than merely repeating them:

- v14 formally signed and passed V, then R failed closed with worker exit `70`
  because an exact `{path,mode,size_bytes,state}` private-key state relation was
  incorrectly classified as byte-readable after that path had already been bound
  metadata-only. Confirm its immutable failure archive, signature, reviews,
  verification receipt, stderr, source identities, key lifecycle, and absence of
  replay/continuation/canonical downstream outputs are all bound by v15.
- In v15 R, SHA-bearing records remain full content relations. Only the exact
  four-key `{path,mode,size_bytes,state}` schema with private mode `0o0` or `0o400`
  is metadata-plus-size: it must use `O_PATH` plus `fstat`, validate exact size,
  and never read the private key bytes. Confirm the per-reader path/schema
  registry still rejects incompatible reclassification.
- Confirm regression coverage proves an actual archived v14 V-receipt private-key
  state fails under the v14 canary and passes under v15 with cached bytes `None`.
  Check rejection of missing `state`, extra keys, wrong mode, boolean size, and
  empty state.
- Preserve all previous security properties: exact typed relation schemas,
  bool-as-int rejection, six-key executable alias equality, process-lifetime
  alias/target inotify, supervisor `waitpid` normal-exit-zero witness, complete
  V/R/C state machines, persistent ceremony watch, no-clobber publication,
  descriptor/source binding, signal/crash windows, and private-key retirement.
- Confirm v15 uses fresh authority and terminal keys while v14/v4 failed evidence
  remains immutable and the legacy v2 signed terminal contract is not rewritten.
- Check the exact 25 recovery output slots, 206 forbidden ceremony slots, 14
  staging patterns, review-state transition, and that the readonly probe cannot
  produce protected output.
- Same-UID malicious process mutation is outside the declared trusted-account
  boundary; do not silently broaden the threat model. Identify any false
  authority inside the declared crash/interruption/race boundary.
- The repair changes provenance and control flow only. It must not change the
  seven datasets, 469,849 same-run LongPhase-S recalibrated FILTER=PASS sSNVs,
  BAM/read-tag inputs, scientific results, or claim ceiling. The ceiling remains:
  methylation can nominate reproducible latent molecular substructure beneath an
  ancestral ALT, but one focal locus does not prove a cellular subclone or linear
  ancestry; currently confirmed cellular subclones and linear ancestry calls are
  both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `1c6eff6464d514721d9b8275cea42892311e2b995999ac0841bb764bd24fbef4`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `957cb61147ef9cd0106437b5649b20b5933fedc9a9dfacc8b10de48f37a89966`
- `scope_sha256`: `9d2ff1d384e37d3b1525391f36c3a38b103667190a937be743fd37901bed661b`
- `terminal_key_rotation_sha256`: `aab7d446e2102f7325bf772eafc21c4a78746b6aa108b65a4883ad70ba7e10b5`
- trusted authority public-key SHA-256: `797c9c174b72ca588d6a63b16c4ab1f0b1bf465763a4ee030e367e1e807aaf4d`
- recovery terminal public-key SHA-256: `17a0b67f22b706da120ca46ce37e9104cc53ba2fd524fd37c92b261668e00f84`

Verdict rule: `APPROVE` only if the exact probe passes without protected writes
and `high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`8540c6b8-509c-404c-9575-07662444c433`.
