<!--
建立時間: 2026-07-22
目標: 對 recovery v5 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: BaseException退役、retained-FD publication validation、terminal stage recheck
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260722_recovery_v4_formal_rejection_and_key_retirement.v1.json
-->

You are the external formal reviewer for InterSubMod Task-B transition-relation recovery v5.
This is READ-ONLY. Never edit, write, create, move, chmod, or delete files. Use only Read,
Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v5.py`: size `55850`, SHA-256
   `193a6234e7029c709de62c9ec674c67331a4f41346c18a9dbe002f25b41cfd75`
2. `verify_tumor_ref_receipt_promotion_recovery_v5.py`: size `125930`, SHA-256
   `7db85f5cc0db709e17c2d5df4a81d6b970b50d504417e8a144b84749090e1c17`
3. `replay_m2v5_runner_only_gates_recovery_v5.py`: size `131056`, SHA-256
   `a96d36ca192cccf8cac96f8c250cb1232689e291c9bb61b2969195380b53a438`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v5.py`: size `290849`, SHA-256
   `74b2d5ed95b07f028473b05ceeaef8ab0cb65139650dc04d4a6343a5f26ecd1f`
5. `probe_tumor_ref_schema_recovery_sources_v5.py`: size `23828`, SHA-256
   `159443e47432d4f68fae75a03152043823afcadfd8155e4822354f6ae57c37bc`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v5.py`: size `35921`, SHA-256
   `9a0fb4e4bb8a0b99344e3cef3664d1b6fc3d2a7679951d0115b583476c1043e8`
7. `build_tumor_ref_schema_recovery_authority_v5.py`: size `31932`, SHA-256
   `cf4bc482eeb5a6068ea0d4e4ec5df0da92c5f80593d5ce074df46493b55e245f`

Run this exact probe and confirm zero output writes:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v5.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `45 passed`, 60
forbidden output slots absent, seven reviewed sources, eleven prior inputs, and
`review_evidence_state=all_absent`.

Append-only history: signed v1 produced authority and Vrec; Rrec v1 failed before output.
Unsigned v2, v3, and v4 were rejected before signing and their private keys are mode 0000.
Their rejection-evidence SHA-256 values are respectively
`f0c64e46ab34872149ca0cae830e79f48d818b591f8100160ac0cce44ed5ba84`,
`c31009328e2130449422e01a5fd766446f6673e6f2625be0f4f380e3f41e4ef5`, and
`a1044ae1a0580b9a6587e30ddfeea22afeac84128cae0d8a28ee3a64619b9fb3`.
Confirm all rejected authority/review/runtime/staging outputs remain absent.

The decisive v4 review found one HIGH and two MEDIUM defects. Verify these v5 repairs:

1. `signature_state.retirement_armed` is set before the signing subprocess. The outer
   publication handler catches `BaseException`, not only `Exception`, and unconditionally
   re-runs the idempotent retirement sequence. This must cover KeyboardInterrupt,
   SystemExit, and interruption after `fchmod(0000)` but before fsync/recheck.
2. Builder payload construction uses the bootstrap validator FD record and excludes that
   path from `_records`. After atomic publication, builder does not call
   `validate_recovery_authority()` and therefore cannot reopen the validator pathname after
   committing an irreversible bundle. Instead it verifies final directory identity,
   inventory, exact retained member-FD bytes, signature and commit, and returns records
   computed from those retained bytes. Decide explicitly whether this retained-FD
   publication validation is complete and whether independent consumer runtime validation
   remains correctly fail-closed.
3. `recheck_staged_bundle()` performs signature/commit semantics first. Then
   `terminal_stage_identity_recheck()` revalidates parent, stage path-to-FD inode/mode/nlink,
   exact inventory, every member path-to-FD identity and exact bytes. Within
   `publish_staged_bundle()`, the next operation is `renameat2(RENAME_NOREPLACE)` with no
   subprocess, callback, or validator between terminal recheck and rename.
4. The 45 tests must behaviorally cover two BaseException classes, partial-retirement
   retry, absence of post-publish validator call, two after-external-verification mutation
   cases, and all earlier strict JSON/schema/ctime/source-binding canaries.
5. v5 validator must bind v4 rejection evidence, the retired v4 key, and absence of all v4
   output slots. Probe must cover exactly 60 forbidden slots, eleven prior inputs, v1-v4
   retired keys, active v5 key mode 0400, and all three reviews absent.
6. The original signed promotion authorization/completion key roots intentionally remain
   `20260722_tumor_ref_receipt_promotion_v4`; only schema-recovery generation paths become
   v5. Confirm no scientific artifact or original signed-chain path was accidentally
   version-renamed.
7. Historical transition treatment remains exact eight keys only at the focal
   `during_execution` context, strict historical ctime before current, while current
   identity remains exact nine keys, descriptor bound, mode 0444, nlink 1.
8. Vrec/Rrec/Crec write only distinct `.recovery.v5` outputs. No v1-v4 output is overwritten.
   Scientific payload, seven datasets/469,849 sites, source receipt bytes, downstream gates,
   and biological claim ceiling remain unchanged.
9. Reviews must be pre-existing exact mode-0444 APPROVE JSONs; builder cannot synthesize,
   edit, normalize, or infer them. Confirm no circular trust or partial-as-complete state.

Canonical review bindings:

- `reviewed_source_set_sha256`: `6cdbe386816fdfc22a6605d2a97110fb3cc4bcb2990c35c04188c2b840e02950`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `359c88781c8f41ec8152c7d205ee586d19707fb0c66d132f8e6019241b77ab5a`
- `scope_sha256`: `dd0d1b6d7e4d5ca8027f52da2fb58e4b067a4f59f942e636d4e174bee1fe11bc`
- trusted public key SHA-256:
  `a2e30f9472b41f6eeefa6ddd862d4b640af428960114b55ad2fc3a7cd64f6510`

Verdict rule: `APPROVE` only if probe passes without writes and high findings `[]`, medium
findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`. Low findings may be
recorded only when demonstrably non-blocking.

Return only the caller-supplied structured JSON schema. Use reviewer
`External Claude Opus`, model `claude-opus-4-8`, session ID
`1050e1fa-18e9-4d69-802f-7f2cff654f78`.
