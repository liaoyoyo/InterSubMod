<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v20 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v19 signed failure、mode-000 inotify fallback、terminal-v10 rotation 與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v20.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v20. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Independently inspect the frozen
sources and run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v20.py`: size `341929`, SHA-256 `4e34624e5054b683da9946cf7b843ced01c9fb4f9aab7800df34c2bddd080731`
2. `verify_tumor_ref_receipt_promotion_recovery_v20.py`: size `127906`, SHA-256 `7805b8e14982031d391ddc98d51736b39bb47a9a77f95150afb66bc7a6dd4ce7`
3. `replay_m2v5_runner_only_gates_recovery_v20.py`: size `152178`, SHA-256 `f1de1c5e9a4893a47d401ff1b4c41d469e5c686b3ec80286ad0a7afb1fbc1e51`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v20.py`: size `366859`, SHA-256 `687e92e6d0ed1c4ee8777435b558da03de914b280d7dcbf4483ad66e8d8b1fb4`
5. `probe_tumor_ref_schema_recovery_sources_v20.py`: size `86719`, SHA-256 `432fd8f676f89456a753d5bed86388c1dbfec4aa8d9e9e9c39c94519a0b46a2f`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v20.py`: size `151124`, SHA-256 `f970de837c22c51e068ebd35996d4dbf09ec59eedc5f9676c0c7e5e66a85f3f5`
7. `build_tumor_ref_schema_recovery_authority_v20.py`: size `53672`, SHA-256 `b8f7b47fa0937e39dd071162c2a0b34c62ba33940488364c600772dbdf7ed080`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v20.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `439 passed`,
281 forbidden slots absent, 19 staging patterns absent, seven current frozen
sources, and `review_evidence_state=all_absent` before review publication.

Audit these decisive points:

- The immutable signed v19 failure archive binds the signed authority, three
  approving reviews, successful V/R, exact child stderr and syscall trace. Confirm
  C failed only because a direct inotify watch on an already-retired mode-000 key
  returned EACCES, before any downstream or terminal output was created.
- Confirm v20 permits parent-only fallback only for the exact conjunction
  `lstat mode == 000` and `inotify_add_watch errno == EACCES`; the named parent
  watch remains mandatory. Real transient chmod and replacement regressions must
  detect mutation, while mode 0400 or EPERM must remain ineligible.

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
  failed-v16 v6, failed-v17 v7, failed-v18 v8, and failed-v19 v9 remain mode
  0400 but quarantined, and only fresh v10 may authorize v20 terminal outputs.
- Confirm V/R/C state machines, descriptor leases, mutation watches, waitpid
  witness, no-clobber publication, signal/crash handling, and one-time private
  key retirement remain fail closed.
- Confirm the builder requires exact equality between probe-reported sources
  and subsequently authority-bound sources. Same-UID hostile replacement
  remains outside the declared trusted-account boundary.
- Confirm v1-through-v19 artifacts are never overwritten and all 281 protected
  slots, 19 staging patterns, and 439 tests are bound.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains latent molecular substructure candidates only;
  confirmed cellular subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `595e4169194b87cace651f07d532e0a395be4d5fb2f434cdba519c07bec01467`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `3460173c983daedd1c41ecc13b4881022bbe00f54ee39f82a4d0a3dd3ef80fd3`
- `scope_sha256`: `27aff69c3ecbd6d882b61d7578ffb16e91fd13863d9bbf125def5a5ef87865b6`
- `terminal_key_rotation_sha256`: `0f5508615d371612475e815825ebcfbe2bdf38be44c87eb21976b88746ae3668`
- trusted authority public-key SHA-256: `122fb38c395c37d1a4a2786c385110397db30f4b2db0ae3e4944f55355656fa9`
- recovery terminal public-key SHA-256: `09794c2ce162af3bf2f3117f6d11dea0c4bd626cbe50946267609058c6c0c291`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`84c32c61-b227-4378-b2fb-1169c0821db9`.
