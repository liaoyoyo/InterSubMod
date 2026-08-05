<!--
建立時間: 2026-07-22
目標: 對 schema-recovery 修補版執行外部 Claude Opus 唯讀正式審查
處理範圍: recovery authority validator、Vrec、Rrec、Crec 與 regression probe
關聯檔案: probe_tumor_ref_schema_recovery_sources_v1.py
-->

You are the external formal reviewer for an InterSubMod Task-B verifier schema recovery.

This is a READ-ONLY review. Do not edit, write, create, move, chmod, or delete any file. Do not run git mutation commands. You may use Read, Grep, Glob, and read-only Bash commands. Execute the exact read-only probe below and inspect all four frozen recovery sources plus the legacy signed sources they reference.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Frozen recovery source identities:

1. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v1.py`
   SHA-256: `f79c6321eee7efb0748bb73e9689ea0a39bf7ac28345275c2114838f81191198`, size 25651
2. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/verify_tumor_ref_receipt_promotion_recovery_v1.py`
   SHA-256: `ec4197cf58655fdbed0af363d5a6bb0cb8589c984b0943341d6044f234801b9c`, size 125470
3. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/replay_m2v5_runner_only_gates_recovery_v1.py`
   SHA-256: `0fc1721029eac265112cf89875e5cd30680c7cae7500072322ea50fafab0532f`, size 126273
4. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/continue_m2v5_after_tumor_ref_promotion_recovery_v1.py`
   SHA-256: `ab72c308777e16565dadd469c396e27c68bdc410fd9cb9dae7ae1854a7bbd139`, size 290849

Legacy signed source identities that remain authoritative only for old promotion assertions:

- V2 `03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8`
- R1 `10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694`
- C1 `f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd`

Exact probe command, which must create zero output files:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v1.py
```

The previous pre-sign review rejected the superseded source set for three reasons. Verify the fixes semantically, not by trusting comments or names:

1. The review validator previously accepted any non-empty probe command. It must now require the exact command list.
2. Vrec evidence previously omitted the recovery authority validator artifact, while Rrec required that path in recursively bound relations. The runtime evidence must now include it, and the regression test must exercise relation discovery.
3. Crec terminal verification previously compared legacy signed authorization fields directly to active recovery source and output paths. Legacy assertions must use `authorized_*`; active Crec must be independently bound by recovery authority `runtime_source`. The helper must be on the real `verify_signed_terminal()` path, and tests must prove correct separation passes while old-as-active fails.

Review all of these requirements:

1. Historical `during_execution` has exactly eight fields and did not record `link_count`; no code may infer historical `link_count=1`.
2. Current artifacts retain nine-field identity and live `st_nlink == 1`.
3. Recovery authority signature, pinned public key, retired private key, exact source/review maps, old-chain hashes, commands, schemas, and receipt paths are fail-closed.
4. Vrec relaxes only the exact historical schema mismatch; all old signatures, canonical bytes, runtime sources, and gates remain enforced.
5. Rrec and Crec keep legacy signed bindings distinct from active recovery bindings and never write old V/R/C receipt paths.
6. Rrec FD-runs fresh Vrec and replays only completion-runner physical lines 1-358 as non-authoritative evidence.
7. Crec preserves sealed memfd capability, `waitpid` normal exit-zero proof, immutable inputs/sources, detached terminal signature, one-time key retirement, incident guards, and all 7 datasets / 469,849 sites downstream contracts.
8. The recovery changes no scientific payload, claim ceiling, or biological interpretation.
9. Search for old/new path confusion, schema mismatches, circular authority, TOCTOU, irreversible post-sign failure, and partial-as-complete paths.

Verdict rule: `APPROVE` only if high findings = 0, medium findings = 0, unresolved conditions = 0, probe exit = 0, probe pass = true, and `no_output_writes = true`; otherwise `REJECT`.

Return the exact structured result requested by the caller. Use reviewer `External Claude Opus`, model `claude-opus-4-8`, and session ID `17e187cd-2ccb-444a-bbbd-252c25c9c717`.
