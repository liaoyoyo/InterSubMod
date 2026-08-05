<!--
建立時間: 2026-07-22
目標: 對第四組 schema-recovery frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: recovery validator、Vrec、Rrec、Crec、fault-injection probe
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260722_tumor_REF_schema_recovery預簽審查與修補紀錄_01.md
-->

You are the external formal reviewer for an InterSubMod Task-B verifier schema recovery. This is READ-ONLY. Never edit, write, create, move, chmod, or delete files and never run git mutation commands. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, single-link frozen sources under `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v1.py`: size `28163`, SHA-256 `1dfb3f3efd62068535788a24bf711cabb1ff425a7b09458dcbd182193cb8ac89`
2. `verify_tumor_ref_receipt_promotion_recovery_v1.py`: size `125470`, SHA-256 `f1b29adf1908d58dec4520dba816036bdf27574ceab8dc3c5603c43b63e053c9`
3. `replay_m2v5_runner_only_gates_recovery_v1.py`: size `126273`, SHA-256 `dee756a0ca55884ca88260b3b5f1af71b397b26952a80ad1813855d03733ed3c`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v1.py`: size `290849`, SHA-256 `81175713f1a124d4d60347581723bd8ae8f79f3e4cf9539e14edd1354dbab94d`

Legacy signed sources must remain authoritative only for old promotion assertions: V2 `03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8`, R1 `10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694`, C1 `f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd`.

Run this exact probe and confirm zero output writes:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v1.py
```

The prior formal review rejected the superseded third set because nested `False == 0` could pass, duplicate JSON keys were last-wins, and the retired recovery private key was only pathname-statted. Verify that the fourth set now:

1. Implements recursive type-strict equality for mappings/lists/scalars and rejects bool/int substitutions at every depth.
2. Rejects duplicate JSON keys and NaN/Infinity.
3. Requires review `exit_code` to be exact `int 0`, not `False`.
4. Binds the retired key using `O_PATH|O_NOFOLLOW`, canonical path, single-link regular inode, exact mode `0000`, retains its FD lease, and terminally rechecks descriptor plus live path before returning.
5. Has negative tests that exercise all those failure paths without mutating the real key.

Reverify all earlier gates: full `/usr/bin/env -i` command list; `authority_validator` in Vrec evidence and Rrec recursive relations; historical exact eight-key record without inferred link count; current exact nine-key identity and live nlink one; unchanged old signatures/canonical bytes; Vrec only relaxes the historical schema mismatch; Rrec FD-runs fresh Vrec and replays physical lines 1-358 only; Crec separates legacy `authorized_*` from active recovery `runtime_source`, calls the pure helper from real signed-terminal verification, preserves sealed memfd/waitpid/signature/key retirement/incident guards, and keeps all 7 datasets/469,849 sites plus 32 data contracts. Search for circular authority, TOCTOU, old receipt writes, post-sign inevitable failure, partial-as-complete, and any scientific or claim-ceiling change.

Verdict rule: `APPROVE` only if probe exit `0`, `pass=true`, `no_output_writes=true`, high findings `[]`, medium findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`.

Return only the structured schema supplied by the caller. Use reviewer `External Claude Opus`, model `claude-opus-4-8`, session ID `783f2548-4de3-4c68-a3f5-eec0d0fb0626`.
