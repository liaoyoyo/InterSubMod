You are the external formal reviewer for an InterSubMod Task-B verifier schema recovery.

This is a READ-ONLY review. Do not edit, write, create, move, chmod, or delete any file. Do not run git mutation commands. You may use Read, Grep, Glob, and read-only Bash commands. You must execute the exact read-only probe command below and inspect the four frozen recovery sources plus the legacy signed sources they reference.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Frozen recovery source identities:

1. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v1.py`
   SHA-256: `b93566dd49fb530f3ba48b9d710a9ed6013be83d194f1d16cbda432c743c3abd`
2. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/verify_tumor_ref_receipt_promotion_recovery_v1.py`
   SHA-256: `77701c30f3929a844af28f95b2b8ea60cf0e87e930929d8fbd8afaff7d3afe18`
3. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/replay_m2v5_runner_only_gates_recovery_v1.py`
   SHA-256: `42b6b0995d53001ca72e47c0d8a2b83645ff27d37e1249b7c84cf122e1bc9632`
4. `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/continue_m2v5_after_tumor_ref_promotion_recovery_v1.py`
   SHA-256: `369eb9462fc1f95315dc4732a1354210e2453a2e26a368ce8266a7e3510870ca`

Legacy signed source identities that must remain authoritative only for the old promotion assertions:

- V2 `03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8`
- R1 `10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694`
- C1 `f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd`

Exact probe command (must run; it is designed to make zero output writes):

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v1.py
```

Review requirements:

1. Confirm the historical `during_execution` record has exactly eight fields and did not record `link_count`; reject any design that invents or infers `link_count=1` for that historical observation.
2. Confirm current artifacts still require the nine-field identity and live `st_nlink == 1`.
3. Confirm the recovery authority signature, pinned public key, retired private key, exact source map, review map, old-chain hashes, commands, schemas, and receipt paths are fail-closed.
4. Confirm Vrec relaxes only the exact historical 8-vs-9 schema mismatch. All old promotion signatures, canonical bytes, runtime sources, and existing gates must remain enforced.
5. Confirm Rrec and Crec separate legacy signed bindings from active recovery bindings. New code must not write old V/R/C receipt paths or impersonate old source output.
6. Confirm Rrec still FD-runs a fresh Vrec and only replays completion-runner physical lines 1-358 as non-authoritative evidence.
7. Confirm Crec preserves the sealed-memfd parent capability, wait/exit-zero proof, immutable downstream source/input checks, terminal detached signature, key retirement, incident guards, and full 7-dataset/469,849-site downstream commands.
8. Confirm the recovery does not alter scientific payload, claim ceiling, or biological interpretation.
9. Search for hard-binding mistakes, old/new path confusion, schema mismatches, circular authority assumptions, TOCTOU gaps introduced by the patch, and any path where a partial result could be presented as complete.

Verdict rule:

- `APPROVE` only if there are zero high findings, zero medium findings, no unresolved conditions, the probe exits 0, and it reports `no_output_writes=true`.
- Otherwise `REJECT`.

Your structured result must use reviewer `External Claude Opus`, model `claude-opus-4-8`, and session_id `2cb4daf4-29f4-473f-b2c8-0df3b5ba9ab9`.
