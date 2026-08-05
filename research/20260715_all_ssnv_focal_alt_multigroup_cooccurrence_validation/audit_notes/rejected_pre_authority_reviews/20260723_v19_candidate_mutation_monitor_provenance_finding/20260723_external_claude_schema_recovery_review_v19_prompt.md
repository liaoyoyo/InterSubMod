<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v19 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v18 signed failure、v19 summary alias/canonical target 修復、terminal key rotation 與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v19.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v19. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Independently inspect the frozen
sources and run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v19.py`: size `303580`, SHA-256 `618382f2dd3b6da74a6a965a8d5236acb86cfaff477e35d7e36fd0936936c92d`
2. `verify_tumor_ref_receipt_promotion_recovery_v19.py`: size `127906`, SHA-256 `a24abd4b90b997fc8d5ff16115a195922969260c8cb39aa7c682cde7fea06c4d`
3. `replay_m2v5_runner_only_gates_recovery_v19.py`: size `152178`, SHA-256 `b3486543089a72d6d2a987f054b5261fd307e26ef6301d5959883b487eea9cd9`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v19.py`: size `361428`, SHA-256 `ade2f9b9331b709487e06ca9a72bf46a248f7dc196edc7f01ec655d2472ccc2f`
5. `probe_tumor_ref_schema_recovery_sources_v19.py`: size `82617`, SHA-256 `26ccf2c1f43c83a9fdc1b30234116eafd90939d8b0bd7d486c184b3ff1961829`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v19.py`: size `144833`, SHA-256 `e05ab8bf210a231658bc2240d12427709b8b100fdd1bfc8359e6e05518a6983c`
7. `build_tumor_ref_schema_recovery_authority_v19.py`: size `53317`, SHA-256 `36ac487de3a4e687086eafdcb2170f8c5ae45209ee71f6ea25e250b926cc1129`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v19.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `434 passed`,
266 forbidden slots absent, 18 staging patterns absent, seven current frozen
sources, and `review_evidence_state=all_absent` before review publication.

Audit these decisive points:

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
  slots, 18 staging patterns, and 434 tests are bound.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains latent molecular substructure candidates only;
  confirmed cellular subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `30c6b0c47ecaf59d00dff31ba99bc536a5c75cde32dda15a8f20a7b52e1b99e9`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `2747c8271f46a3cf7583f0e9315d1e6dc8a8fc4f64eefa92c89fac88fda2d59c`
- `scope_sha256`: `641679e2ca6a3dab5369aac978b2cf4cc84757e31ad2958a617a35efa9fe0349`
- `terminal_key_rotation_sha256`: `0030b6315f26b884e6644514c813c98b9e41f5c64e3426ca1cd25bb3f9136d81`
- trusted authority public-key SHA-256: `d494f66e8aea206e37e0e803ccb4e0ceb9cf2b244e71eccba6b370f96ddee2e0`
- recovery terminal public-key SHA-256: `34d11ff6a699b96aaa8624b30f45fbc8845e6958dc04122b19c413a1a2c12c2d`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`006a1050-2d27-4f63-9ee0-80a42fb8f26d`.
