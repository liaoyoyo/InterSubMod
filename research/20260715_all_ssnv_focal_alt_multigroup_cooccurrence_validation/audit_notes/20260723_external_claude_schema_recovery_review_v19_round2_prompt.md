<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v19 round-2 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v18 signed failure、round-1 review rejection、v19 alias/inotify 修復、terminal key rotation 與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v19.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v19 round-2. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Independently inspect the frozen
sources and run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v19.py`: size `314623`, SHA-256 `db693ed4e3dc6d66ddab58d5b3c1fa98981abe71e21905c477ef617646e8ad07`
2. `verify_tumor_ref_receipt_promotion_recovery_v19.py`: size `127906`, SHA-256 `026b62c58d549cd53b63c42dc89145778ef97103a6e2437c4314d48c8107f2ad`
3. `replay_m2v5_runner_only_gates_recovery_v19.py`: size `152178`, SHA-256 `10b5301fe52c075733f347b459740a3e316c39e92b69989de7bae440ec1443b4`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v19.py`: size `361552`, SHA-256 `17e0aa4d83181dbca20ea407020621e0032d66c58bbf2e1a21cb9c89b4a7a702`
5. `probe_tumor_ref_schema_recovery_sources_v19.py`: size `83457`, SHA-256 `4cba3d3bbf3bd1fecd5d2b798d84bdda98f2b0b4f85a60be9f20a422e1996957`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v19.py`: size `147278`, SHA-256 `a8e3152d9e77361eb66ba8486e22d90cb094bce6535da37066a96851906047f1`
7. `build_tumor_ref_schema_recovery_authority_v19.py`: size `53672`, SHA-256 `c45cb9cf69dc59496493aea8e5de21260a4bd12c398b72fb01ad48e50cd3c1f4`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v19.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `436 passed`,
266 forbidden slots absent, 18 staging patterns absent, seven current frozen
sources, and `review_evidence_state=all_absent` before review publication.

Audit these decisive points:

- The immutable v19 round-1 rejection archive binds both independent
  `REQUEST_CHANGES` reviews, the old seven-source snapshot, and the aborted
  external envelope. Confirm the active source signs both all-event directories
  and the real transient replace-and-restore inotify regression passes.

- The immutable v18 archive proves a signed authority, successful V/R and fresh
  V, then a C child failure before downstream output creation because the live
  tumor-REF `summary.json` is a relative symlink, not a canonical regular path.
- v19 must bind the canonical `all_ssnv_tumor_ref_control_summary.json` bytes and
  descriptor first, then bind `summary.json` as an exact relative basename
  symlink to that target using `O_PATH|O_NOFOLLOW`. Both alias and target identity
  must remain stable for process lifetime and be covered by mutation monitoring.
- Confirm real regressions accept only the exact relative alias and reject an
  unbound target, wrong target, absolute link text, alias replacement, and
  canonical-target replacement. Confirm `GATE_INPUT_PATHS["tumor_ref_summary"]`
  points to the canonical target rather than the symlink.
- Confirm terminal keys are pairwise distinct: legacy v2 remains preserved,
  failed-v16 v6, failed-v17 v7, and failed-v18 v8 remain mode 0400 but
  quarantined, and only fresh v9 may authorize v19 terminal outputs.
- Confirm V/R/C state machines, descriptor leases, mutation watches, waitpid
  witness, no-clobber publication, signal/crash handling, and one-time private
  key retirement remain fail closed.
- Confirm the builder requires exact equality between probe-reported sources
  and subsequently authority-bound sources. Same-UID hostile replacement
  remains outside the declared trusted-account boundary.
- Confirm v1-through-v18 artifacts are never overwritten and all 266 protected
  slots, 18 staging patterns, and 436 tests are bound.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains latent molecular substructure candidates only;
  confirmed cellular subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `bf34658b2485ae84b08d7ff32db5c588383e68ac99c005a77f6747c37253defa`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `3460173c983daedd1c41ecc13b4881022bbe00f54ee39f82a4d0a3dd3ef80fd3`
- `scope_sha256`: `f27fa134947a72bc2c983dbde931cfdf77055a91cd7966a0b0e7cbda0a460b51`
- `terminal_key_rotation_sha256`: `0030b6315f26b884e6644514c813c98b9e41f5c64e3426ca1cd25bb3f9136d81`
- trusted authority public-key SHA-256: `d494f66e8aea206e37e0e803ccb4e0ceb9cf2b244e71eccba6b370f96ddee2e0`
- recovery terminal public-key SHA-256: `34d11ff6a699b96aaa8624b30f45fbc8845e6958dc04122b19c413a1a2c12c2d`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`db7885d1-6d8c-4bdb-8fcc-7bb3f6f03124`.
