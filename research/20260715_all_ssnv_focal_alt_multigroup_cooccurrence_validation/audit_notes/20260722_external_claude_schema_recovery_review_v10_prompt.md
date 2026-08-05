<!--
建立時間: 2026-07-22
目標: 對 recovery v6 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v5 blocker closure、successful-publish behavioral test、append-only recovery chain
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260722_recovery_v5_formal_rejection_and_key_retirement.v1.json
-->

You are the external formal reviewer for InterSubMod Task-B transition-relation recovery v6.
This is READ-ONLY. Never edit, write, create, move, chmod, or delete files. Use only Read,
Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v6.py`: size `58279`, SHA-256
   `f571d7c0a58242485ac78439277d601502e38236b63db5ded0d38a9d8aa2e2cc`
2. `verify_tumor_ref_receipt_promotion_recovery_v6.py`: size `125930`, SHA-256
   `c5331fa55705b594dafe936a3b51ab9379291c66367b0f81334735246f79150c`
3. `replay_m2v5_runner_only_gates_recovery_v6.py`: size `131056`, SHA-256
   `5d30c53d2a143b7d0192f87fa52d6fdc68ccdb2995704de27f2f95b2fd0c692f`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v6.py`: size `290849`, SHA-256
   `77fe10613c38752fa7770205e04fdf28036e32dc9210b1ac2ec42c50f625dcd0`
5. `probe_tumor_ref_schema_recovery_sources_v6.py`: size `24933`, SHA-256
   `c803b8423dbc9f003de49a053a32af4b37f145dd057800cb143dd2dfddaf6879`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v6.py`: size `38376`, SHA-256
   `735d9612b2e9eb40ab58f2f5b04a4d7b0ecb981160ffa8376bcc70fbf9db0709`
7. `build_tumor_ref_schema_recovery_authority_v6.py`: size `32017`, SHA-256
   `b5d34c22d580da7643ed45ef0ed72002cbb02c9db0b933a0a2088ae2e75416fd`

Run this exact probe and confirm zero output writes:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v6.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `46 passed`, 70
forbidden output slots absent, seven reviewed sources, twelve prior inputs, and
`review_evidence_state=all_absent`.

Append-only history: signed v1 produced authority and Vrec; Rrec v1 failed before output.
Unsigned v2, v3, v4, and v5 were rejected before signing and their private keys are mode
0000. Their rejection-evidence SHA-256 values are respectively
`f0c64e46ab34872149ca0cae830e79f48d818b591f8100160ac0cce44ed5ba84`,
`c31009328e2130449422e01a5fd766446f6673e6f2625be0f4f380e3f41e4ef5`,
`a1044ae1a0580b9a6587e30ddfeea22afeac84128cae0d8a28ee3a64619b9fb3`, and
`0495af623ba463b822f4f823ce28d49745c11c056ae391a36deeda06e0e78047`.
Confirm all rejected authority/review/runtime/staging outputs remain absent.

The decisive v5 review had no HIGH finding and one MEDIUM finding:
`M1_missing_behavioral_post_publish_reopen_test`. Audit these v6 requirements:

1. `test_v6_successful_publish_uses_retained_fds_without_runtime_validator` must execute
   the real successful `publish()` path, atomically expose the exact three-member bundle,
   retire the one-time key, assert retained-FD publication evidence, and sentinel-test
   that `validate_recovery_authority()` is called zero times. A static `co_names` check
   alone is insufficient.
2. The complete 46-test suite must cover the new success path plus KeyboardInterrupt,
   SystemExit, partial-retirement retry, post-external-verification member/directory drift,
   strict JSON/schema/ctime/source-binding canaries, and all prior failure generations.
3. Builder payload construction uses the bootstrap validator FD record and excludes that
   pathname from `_records`. Post-publication validation uses retained stage/member FDs;
   independent consumer runtime validation remains fail-closed.
4. `recheck_staged_bundle()` performs signature/commit semantics first, followed by the
   terminal directory/member identity and exact-byte recheck; the next operation is
   `renameat2(RENAME_NOREPLACE)` with no callback, subprocess, or validator in between.
5. Signing arms retirement before the signing subprocess. The outer handler catches
   `BaseException` and idempotently retries retirement using the same bound FD, including
   interruption after chmod(0000) but before fsync/recheck.
6. v6 validator binds the exact v5 rejection evidence, retired v5 key, and absence of all
   v5 output slots. Probe covers exactly 70 forbidden slots, twelve prior inputs, v1-v5
   retired keys, active v6 key mode 0400, and all three v6 reviews absent.
7. The original signed promotion key roots intentionally remain
   `20260722_tumor_ref_receipt_promotion_v4`; only schema-recovery generation paths become
   v6. Scientific artifact names ending `source_authority_v5` remain unchanged.
8. Historical transition treatment remains exact eight keys only at the focal
   `during_execution` context, with strict historical ctime before current. Current
   identity remains exact nine keys, descriptor bound, mode 0444, nlink 1.
9. Vrec/Rrec/Crec write only distinct `.recovery.v6` outputs. No v1-v5 output is
   overwritten. Scientific payload, seven datasets/469,849 sites, source receipt bytes,
   downstream gates, and biological claim ceiling remain unchanged.
10. Reviews must be pre-existing exact mode-0444 APPROVE JSONs. Builder cannot synthesize,
    edit, normalize, or infer reviews. Confirm there is no circular trust or partial state.

Canonical review bindings:

- `reviewed_source_set_sha256`: `87d68c542d9461d59d96a6d32b006467e4401f9a80e2c65ad878dbeffd632622`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `6cbfb83b214900d281ea675898333f613a42e020798ed22d4eb35ccf3eb0ca45`
- `scope_sha256`: `08e1aa0af13275848aa81fc0983c6bfb6fda87bc83c44a7b044acba7edcca45f`
- trusted public key SHA-256:
  `9fb78aa6bb87e87c40ebdf20451d3503e7e497e5ab6cf9001644054275ed83ee`

Verdict rule: `APPROVE` only if the probe passes without writes and high findings `[]`,
medium findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`. Low findings
may be recorded only when demonstrably non-blocking.

Return only the caller-supplied structured JSON schema. Use reviewer
`External Claude Opus`, model `claude-opus-4-8`, session ID
`fcaed944-0b5f-44d7-8521-843f1a38db29`.
