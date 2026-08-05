<!--
建立時間: 2026-07-22
目標: 對 recovery v3 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: single-FD builder、terminal recheck、strict-forward ctime、atomic authority bundle
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/pre-decision-audit-recovery.md
-->

You are the external formal reviewer for InterSubMod Task-B transition-relation recovery v3.
This is READ-ONLY. Never edit, write, create, move, chmod, or delete files. Never run git
mutation commands. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v3.py`: size `48966`, SHA-256
   `ceedcb2319cb3c6bc1908a0aaa4e67a318d9e4031b7208dce68cd94ff3358b45`
2. `verify_tumor_ref_receipt_promotion_recovery_v3.py`: size `125930`, SHA-256
   `36c312341795030c4f5202f8e028c7671697071a60bb26848c324d50b9d7a157`
3. `replay_m2v5_runner_only_gates_recovery_v3.py`: size `131056`, SHA-256
   `5757476f47ea70bc53b3930a8f144e1354c5351335aff7e6d48795a02907ed5c`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v3.py`: size `290849`, SHA-256
   `bb07d295724d65f80ed6c7c1af95611ba81b4290039f98e19faa402cdab32641`
5. `probe_tumor_ref_schema_recovery_sources_v3.py`: size `21449`, SHA-256
   `c0ec9f86c08ebb4b7895a6cb1d1907d99843ed8990bb02b6a7c416d9541e180f`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v3.py`: size `25442`, SHA-256
   `45e20cc0e2f4e77807b11bc343e24d186e37ce44b7a2e21580105f405637e4e0`
7. `build_tumor_ref_schema_recovery_authority_v3.py`: size `26884`, SHA-256
   `daaedcee09e2c1ff64f85d17fe7722b888d6c02784b8002ff42e0f683ec44e0e`

Run this exact probe and confirm zero output writes:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v3.py
```

Required pre-review result: exit `0`, `pass=true`, `no_output_writes=true`, `34 passed`,
40 forbidden output slots absent, seven reviewed sources, nine prior inputs, and
`review_evidence_state=all_absent`.

The prior signed recovery v1 produced authority and Vrec, then Rrec v1 failed before any
Rrec output because it rebound the signed historical eight-key `during_execution` record
as current live identity. Recovery v2 corrected that relation but was rejected before
signing by Nash because its builder executed a second unhashed validator pathname read,
did not terminally recheck all inputs before signing/publication, and required only ctime
inequality. Confirm v2 authority/reviews/runtime outputs remain absent, its private key is
mode `0000`, and rejection evidence SHA-256 is
`f0c64e46ab34872149ca0cae830e79f48d818b591f8100160ac0cce44ed5ba84`.

Verify all v3 corrections and unchanged contracts:

1. Builder opens validator exactly once with `O_NOFOLLOW`, hashes and compiles the same
   retained FD bytes, and never executes a second pathname read.
2. Recovery, legacy, original chain, v1 predecessor, rejected v2 evidence, reviews, public
   key, and retired keys are held as descriptor leases with terminal stat and SHA rechecks.
   The active v3 private key is one retained FD from mode `0400` through signing and mode
   `0000` retirement.
3. Input barriers run before signing and again before atomic publication. Drift prevents
   the next side effect. OpenSSL uses retained public/private/data/signature FDs.
4. Authority, signature, and commit are first written and verified inside a mode-0700
   hidden staging directory. Final publication is one `renameat2(RENAME_NOREPLACE)` of the
   complete directory; there is no sequential-file fallback. Consumers accept only the
   final bundle with exact three-member inventory, commit hashes, signature, transaction ID,
   retired key, and runtime validator pass.
5. If any failure occurs after signing, the v3 key is retired and staging evidence is kept;
   no cleanup or unsafe retry occurs. Preflight writes nothing and leaves key mode `0400`.
6. Historical transition exemption remains only at exact JSON pointer
   `evidence.focal_source_identity_transition.during_execution`; historical record is exact
   eight keys without inferred link count; current record is exact nine keys, live bound,
   mode `0444`, nlink `1`.
7. Both Vrec and Rrec require `historical ctime_ns < current ctime_ns`. Equal and reversed
   ctime fail even when the difference map is made structurally consistent. Authorization
   and promotion transitions remain type-strictly identical.
8. Tests retain the v1 walker failure canary and add v2 second-read reproduction, v3
   single-open, stat/content drift, private retirement, atomic rename occupied/success, and
   commit-binding tests. Duplicate keys, NaN, `1e400`, floats, nested bool/int, moved context,
   wrong keys/mode/link count/difference map/interpretation all fail closed.
9. Vrec/Rrec/Crec use active v3 paths and distinct `.recovery.v3` outputs. v1 is the only
   signed predecessor; v2 is unsigned rejected evidence only. Builder cannot synthesize or
   transform reviews and accepts only three pre-existing mode-0444 exact APPROVE JSONs.
10. Confirm no circular authority, exploitable TOCTOU, review fabrication, old output write,
    partial-as-complete state, post-sign inevitable failure, or broad mode exemption.
11. Confirm sealed memfd, live-parent binding, exact-child normal exit zero, detached terminal
    signature, continuation key retirement, incident fail-closed, runner lines 1-358 only,
    all downstream output inventory, seven datasets/469,849 sites, and the biological claim
    ceiling are unchanged.

Verdict rule: `APPROVE` only if probe exit `0`, `pass=true`, no writes, high findings `[]`,
medium findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`. Low findings may
be recorded only when demonstrably non-blocking.

Return only the caller-supplied structured JSON schema. Use reviewer `External Claude Opus`,
model `claude-opus-4-8`, session ID `614c2043-d3ee-4eb9-8c8d-a4ef48e53e43`.
