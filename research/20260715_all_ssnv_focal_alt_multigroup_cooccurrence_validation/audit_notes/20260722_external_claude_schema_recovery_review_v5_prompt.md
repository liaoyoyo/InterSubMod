<!--
建立時間: 2026-07-22
目標: 對 transition-relation recovery v2 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: validator、Vrec、Rrec、Crec、probe、tests、ceremony builder
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/pre-decision-audit-recovery.md
-->

You are the external formal reviewer for an InterSubMod Task-B transition-relation
recovery v2. This is READ-ONLY. Never edit, write, create, move, chmod, or delete files;
never run git mutation commands. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, single-link frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v2.py`: size `36899`, SHA-256 `f216d4c1d18fc364b3393308c015ce73d06fc335774bf4d5d5e8a75b8c765e9b`
2. `verify_tumor_ref_receipt_promotion_recovery_v2.py`: size `125930`, SHA-256 `9d3ac1621871a439a8d0ae56fd4b71f986bf496287d0d731a3d17b2c183bfd85`
3. `replay_m2v5_runner_only_gates_recovery_v2.py`: size `131056`, SHA-256 `fe22cf926d398eaff95aeba37f833a8671087eec7404ce1f22665377905c83ad`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v2.py`: size `290849`, SHA-256 `42a6fe01dcf8105dcf4133ec064816d481401d530f34923d67b76312c8def1c7`
5. `probe_tumor_ref_schema_recovery_sources_v2.py`: size `19453`, SHA-256 `b236cb660a0bb0407c3f8001b43dba5826ef2a621a513ceb8d0edca83e2a807d`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v2.py`: size `19245`, SHA-256 `6816c87e158afa9c191385ff05ce8be5e5a0cb3e4f5283552d5b0be7b0beaf19`
7. `build_tumor_ref_schema_recovery_authority_v2.py`: size `15321`, SHA-256 `e5bfc9fa23daeb21bb52bba6dba63103df69ae6abcc338ef149be5a9bb5ad549`

Run this exact probe and confirm zero output writes:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v2.py
```

The required pre-review result is exit `0`, `pass=true`, `no_output_writes=true`,
`25 passed`, 30 forbidden output slots absent, seven reviewed sources, eight prior inputs,
and `review_evidence_state=all_absent`.

The prior signed recovery v1 correctly produced an authority and Vrec receipt, then Rrec
v1 failed before producing log/receipt because its generic walker rebound the signed
historical `during_execution.mode=0o664` record as current live identity `0o444`.
Verify the v2 correction:

1. The exemption exists only at exact JSON pointer
   `evidence.focal_source_identity_transition.during_execution`.
2. Historical identity is exact eight keys with no inferred link count; current identity
   is exact nine keys, descriptor-bound live, mode `0444`, nlink `1`.
3. Path/hash/size/device/inode/mtime remain exact; only signed `0664 -> 0444` mode and
   corresponding differing ctime are allowed. Interpretation/difference mapping is exact.
4. Authorization and promotion transitions must type-strictly match.
5. Moved/duplicate/missing/extra keys, bool-for-int, mode/path/hash/stat drift, ctime equality,
   difference-map drift, and unexpected transition contexts all fail closed before writes.
6. The old v1 walker failure remains a regression canary.

Verify the full source cohort and trust chain:

1. Vrec/Rrec/Crec all use active v2 source paths, v2 authority, and distinct
   `.recovery.v2` receipts. The signed legacy `authorized_*` paths remain unchanged.
2. v1 authority/signature/reviews/Vrec receipt/failure evidence are immutable predecessor
   inputs only; v1 private key remains retired mode `0000`; v2 key is `0400` pre-sign.
3. Validator independently verifies v1 Ed25519 signature and exact hashes, rejects duplicate
   keys, NaN/Infinity, every float literal including `1e400`, and nested bool/int substitutions.
4. Probe/test/builder are included in the seven-source signed review set.
5. Builder cannot synthesize review evidence. It must ingest three already-existing mode-0444
   review JSONs, validate exact source/prior/scope digests and severity arrays, then publish
   only authority/signature with O_EXCL before retiring the v2 key.
6. Search for circular authority, TOCTOU, review fabrication, old receipt writes, post-sign
   inevitable failure, partial-as-complete, and broad mode exemption.
7. Confirm sealed memfd, live-parent binding, exact-child normal exit zero, detached terminal
   signature, continuation key retirement, incident fail-closed, runner lines 1-358 only,
   all downstream output inventory, seven datasets/469,849 sites, and the biological claim
   ceiling are unchanged.

Verdict rule: `APPROVE` only if probe exit `0`, `pass=true`, no writes, high findings `[]`,
medium findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`. Low findings may be
recorded only when demonstrably non-blocking.

Return only the structured schema supplied by the caller. Use reviewer
`External Claude Opus`, model `claude-opus-4-8`, session ID
`c1bc858d-7bfb-4289-a286-91b4e7d05f9a`.
