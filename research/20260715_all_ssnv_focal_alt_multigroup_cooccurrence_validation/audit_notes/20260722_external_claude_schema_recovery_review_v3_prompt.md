<!--
建立時間: 2026-07-22
目標: 對最終 schema-recovery frozen source set 執行外部 Claude Opus 唯讀正式審查
處理範圍: recovery validator、Vrec、Rrec、Crec、regression probe
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260722_tumor_REF_schema_recovery預簽審查與修補紀錄_01.md
-->

You are the external formal reviewer for an InterSubMod Task-B verifier schema recovery. This is READ-ONLY: never edit, write, create, move, chmod, or delete files; never run git mutation commands. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these frozen mode-0444, single-link sources:

1. validator: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v1.py`, size `26003`, SHA-256 `208e3e598ee9b2445b8aebd443a174fb5f18fc50ac09a7b6d46cb5e2bdc83351`
2. Vrec: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/verify_tumor_ref_receipt_promotion_recovery_v1.py`, size `125470`, SHA-256 `def82760ce80cc41c3c696397ac2eccfb468b6f3525d95dde82ccdcf62056b42`
3. Rrec: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/replay_m2v5_runner_only_gates_recovery_v1.py`, size `126273`, SHA-256 `6ad5e79e939ded1849e1bf59cf3db2e69c07b90ac5b8d1e928ae6c95b15c8035`
4. Crec: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/continue_m2v5_after_tumor_ref_promotion_recovery_v1.py`, size `290849`, SHA-256 `e5c685aba84673eee0d5ee5682e40800b904e77d5fbe3fe6f2e3629b7e92d31d`

Legacy signed source SHA-256 values must remain authoritative only for old promotion assertions: V2 `03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8`, R1 `10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694`, C1 `f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd`.

Run this exact probe and confirm it creates zero outputs:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v1.py
```

Adversarially verify these corrected defects:

1. Review artifacts must bind the full `/usr/bin/env -i ...` command list exactly, including every clean-room environment token; any shorter/non-empty list must fail.
2. Runtime recovery evidence must include `authority_validator`, so recursive Vrec-to-Rrec relation discovery contains `SCHEMA_RECOVERY_VALIDATOR`.
3. Crec signed-terminal legacy assertions must use only `authorized_*` source/path/commands, while active Crec is independently required as recovery authority `runtime_source`. The pure helper must be called by the real `verify_signed_terminal()` path. Tests must prove valid separation passes and old-as-active fails.

Also verify: historical `during_execution` is exact 8-key with no inferred `link_count`; current artifacts remain exact 9-key with live `st_nlink=1`; authority signature/public key/retired private key/source and review maps/old-chain hashes/commands/schemas/receipt paths are fail-closed; Vrec relaxes no condition except the exact historical schema mismatch; Rrec FD-runs fresh Vrec and replays only completion-runner physical lines 1-358; recovery never writes old V/R/C receipt paths; Crec preserves sealed memfd, waitpid normal exit zero, immutable inputs, detached terminal signature, key retirement, incident guards, all 7 datasets and 469,849 sites, and no partial-as-complete path. Confirm no scientific payload, canonical receipt, claim ceiling, or biological interpretation changes. Search for circular authority, TOCTOU, old/new path confusion, and irreversible post-sign failures.

Verdict rule: `APPROVE` only if probe exit `0`, probe `pass=true`, `no_output_writes=true`, high findings `[]`, medium findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`.

Return only the caller's structured schema. Use reviewer `External Claude Opus`, model `claude-opus-4-8`, session ID `0b3fefd6-3187-4e1a-913d-af5d3ea44d4d`.
