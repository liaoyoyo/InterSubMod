<!--
建立時間: 2026-07-24
目標: 對 tumor-REF schema recovery v28 的十個 frozen source 執行外部 Claude Opus 唯讀正式重審
處理範圍: failed-v26 projection、rejected-v27 runtime/transport findings、14-role registry、v1/v6 lineage、完整 inventory regressions 與 terminal-v18 rotation
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v28.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v28. This is a READ-ONLY remediation re-review after the v27
pre-authority candidate was rejected because independent reviewers reproduced
stale current-runtime hashes, a review-transport label mismatch, and incomplete
v27 staging coverage. Never edit, write,
create, move, chmod, or delete files. Do not use network tools. Independently
inspect the frozen sources and run the exact probe; do not treat this prompt as
evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v28.py`: size `520943`, SHA-256 `9c1b4b93a60c28a5ab16832f01ef4366bcf994411e2d090a4a1dcf8e941ee946`
2. `verify_tumor_ref_receipt_promotion_recovery_v28.py`: size `129377`, SHA-256 `d6de79445cb6c1be329cf8ebd4b4907a89d7caaf80e99fdceb997636a212ce83`
3. `replay_m2v5_runner_only_gates_recovery_v28.py`: size `153502`, SHA-256 `6a80d8dd5eac8b47bac48453d394e0fcc5d685613418b718c23fdf9d37acdef9`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v28.py`: size `407346`, SHA-256 `6d59f6d30edd01ea0cb220e93cae3cf8d970b96bd4a646c2e8580ea8a97d8cfc`
5. `probe_tumor_ref_schema_recovery_sources_v28.py`: size `144464`, SHA-256 `3526b9677e4f58a81c6258ae10ccca5dd8034eab7d3f178d68ddb76cb5b694c0`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v28.py`: size `196816`, SHA-256 `2153cd24be023d79052144cc723a78bc3408541c332409ac8fbbf63182eed5b2`
7. `build_tumor_ref_schema_recovery_authority_v28.py`: size `61065`, SHA-256 `b5dd46bdd7e771572731b8e7a738c79c4a0a344007688bedd11cf2077b6fb9f4`
8. `build_all_ssnv_final_report_dataset_schema_recovery_v28.py`: size `351388`, SHA-256 `cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f`
9. `finalize_task_b_result_release_schema_recovery_v28.py`: size `33500`, SHA-256 `da84f1f956ce06d6d2640ee56e739a17e98859eae3cfbb15419e86d397aee55c`
10. `build_all_ssnv_report_artifact_schema_recovery_v28.py`: size `238722`, SHA-256 `f0afa66823b63d692ea519a639e65c62ed01acc5c4c6be38f77f1391ed98b9f3`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v28.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `674 passed`,
405 forbidden slots absent, 27 staging patterns absent, ten current frozen
sources, `review_evidence_state=all_absent`, signed failed-v26 quarantined, and
rejected-v27 sources/reviews/transport/two unused key pairs hash-bound.

Audit these decisive points:

- The immutable failed-v26 archive must preserve the signed authority, signature,
  commit, terminal-v16 key state, exact failure evidence, and every protected
  formal/runtime output state. The rejected-v27 archive must preserve all ten
  frozen sources, all three candidate reviews, all six review-transport files,
  and both unused key archive records. No v26/v27 key may authorize v28 output.
- Confirm v28 represents the immutable tumor-REF v1 pre-audit and current v6
  primary artifacts as a chronological lineage rather than pretending they are
  the same path. Their data planes must remain exactly 102,842 sites, 308,526
  artifacts, and artifact-set SHA-256
  `195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca`.
- Legacy v1 and current v6 tumor-REF inputs must be opened with `O_NOFOLLOW` and
  parsed from the same descriptor-bound bytes used for their SHA-256 records.
  The continuation gate and sentinel must bind `LEGACY_PRIMARY_PRE`, validate
  the manifest crosslink, and detect mutation before downstream execution.
- The v26 signed commit must be validated as an exact transaction: member set,
  authority and signature digests, retired key records, and bundle identity must
  all match. Mutation regressions must fail closed.
- Regression tests must execute the descriptor-bound primary Python, QA Python,
  canonical Node v22.22.1, and QA Chromium binaries. Production continuation
  validation must perform the same descriptor-bound checks before downstream
  work and use the canonical v5 runtime output slots.
- Occupied-state regressions must cover every one of the 405 forbidden output
  slots and every one of the 27 staging patterns. Source, probe, and bound test
  execution contracts must remain self-enforcing.
- The exact canonical builder source must remain bound before validator load;
  `__file__`, `sys.argv[0]`, `/proc/self/cmdline`, path, and inode must agree.
  A copied builder must fail before any output.
- Review attribution is deliberately limited: reviewer IDs are
  orchestrator-recorded transport identifiers, not reviewer-held
  cryptographic signatures and not independent proof of reviewer authorship.
- Confirm all fourteen terminal public keys are pairwise distinct: legacy v2;
  failed-v16 v6, failed-v17 v7, failed-v18 v8, failed-v19 v9, rejected-v20 v10,
  failed-v21 v11, rejected-v22 v12, rejected-v23 v13, rejected-v24 v14, and
  failed-v25 v15, failed-v26 v16, and rejected-v27 v17 remain
  historical/quarantined; only fresh v18 may authorize v28 terminal outputs.
- Confirm the historical signed runtime projection is an explicit exact
  11-role allowlist, recovery-only roles are an explicit disjoint 3-role set,
  and the current reviewed runtime set is their exact 14-role union. The three
  recovery roles must never enter historical authorization/completion payloads.
- Confirm `GATE_INPUT_PATHS` binds all 14 current roles under canonical role
  names, including both builders and the result/report finalizer, with no alias,
  missing role, swapped path, or duplicate-path role accepted.
- Confirm V/R/C state machines, descriptor and directory leases, mutation
  watches, waitpid witness, no-clobber publication, signal/crash handling,
  exact relation registries, and one-time private-key retirement remain fail
  closed. Same-UID hostile runtime injection remains outside the declared
  trusted-account boundary.
- The repair changes provenance/runtime validation only. It must not change the
  seven datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS`
  biallelic sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific
  values. The claim ceiling remains latent molecular substructure candidates
  only; confirmed cellular subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `0e5aabe293c7b1ff27597285a2677ac2498e980a7ee244de748776b35604a0a6`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `bbd3ccb4920730ec23f08ddb21550b79b67ea0dd4d48a14908818bd9cd22ffee`
- `scope_sha256`: `a3e2c44caf0db17cc2cef6a4341ce179b58f1b5baa48d251951a651d0c0b905f`
- `terminal_key_rotation_sha256`: `60afefc38bd7ecc275fe9c5fbd3d4ca8c7327a136a1dc0e4c54adfb89164e805`
- trusted authority public-key path: `/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v28/ed25519_public.pem`
- trusted authority public-key SHA-256: `0a4d3a99befa388255c506e5a2a77c2acfd831a8b586d7a5b140585015584e3e`
- recovery terminal public-key SHA-256: `ded65f60581476b08c530c04d795c55bcb393558f393acb0f171629dce47add7`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. The attribution object records transport only and
must explicitly set `cryptographic_reviewer_authorship_proven=false`. Return
only the JSON object required by the supplied schema. Use reviewer
`External Claude Opus` and reviewer/session ID
`0149f9af-ea8a-4ce3-81a4-952a469100c9`.
