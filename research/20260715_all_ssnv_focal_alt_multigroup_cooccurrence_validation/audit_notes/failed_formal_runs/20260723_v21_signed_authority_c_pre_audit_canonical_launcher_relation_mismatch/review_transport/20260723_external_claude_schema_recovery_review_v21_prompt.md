<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v21 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v20 pre-authority rejection、parent-watch setup TOCTOU、terminal-v11 rotation 與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v21.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v21. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Independently inspect the frozen
sources and run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v21.py`: size `363327`, SHA-256 `b0761a96d1a96a32805e0a14412b5afb9bb403fb06bc300a9f507b73c3afeb75`
2. `verify_tumor_ref_receipt_promotion_recovery_v21.py`: size `127906`, SHA-256 `e95d646da009ad4457c0c0f71d74f57fd207208d4cf541dec90e3476e729cde6`
3. `replay_m2v5_runner_only_gates_recovery_v21.py`: size `152178`, SHA-256 `4f181e07f9ae2175d3cbb92e80e8f17478d820941fa9f3518c793bbd75007a08`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v21.py`: size `373493`, SHA-256 `82b2b49e1fd3675fb962735ceb3c041ef7f16ad29e8eac3d4027674789bc50eb`
5. `probe_tumor_ref_schema_recovery_sources_v21.py`: size `90656`, SHA-256 `8158a737279078fe6192dd6cf599dbfb16b3e8f068f8018689b924cd01eb3d9a`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v21.py`: size `154038`, SHA-256 `43026ee29f2741d7dc0edc0fc54bbdd759713e10a5710422fbd28c34adc08dc6`
7. `build_tumor_ref_schema_recovery_authority_v21.py`: size `53672`, SHA-256 `f6cbbd41f191114f1b8c7baa79156a05f33cf82fbdd83053bfde06ff39003074`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v21.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `506 passed`,
296 forbidden slots absent, 20 staging patterns absent, seven current frozen
sources, and `review_evidence_state=all_absent` before review publication.

Audit these decisive points:

- The immutable signed v19 failure archive binds the signed authority, three
  approving reviews, successful V/R, exact child stderr and syscall trace. Confirm
  C failed only because a direct inotify watch on an already-retired mode-000 key
  returned EACCES, before any downstream or terminal output was created.
- The immutable v20 pre-authority rejection must bind all seven archived source
  identities, two internal `REQUEST_CHANGES` reviews, one external `APPROVE`
  envelope, one high and two medium blocking findings, both unused mode-0400 keys,
  and all 15 formal v20 slots absent. Confirm strictest-review-wins prevented any
  v20 authority, formal review, V/R/C, or terminal output and both v20 keys are
  quarantined and must never be reused.
- Confirm v21 installs every named-parent watch before taking any protected leaf
  snapshot or attempting direct watches. After all watches exist it must compare
  every leaf identity and mode against the initial snapshot and drain events.
  Deterministic setup-window regressions must reject both transient chmod-restore
  and replacement-restore. Parent-only fallback remains limited to the exact
  conjunction `lstat mode == 000` and direct-watch `errno == EACCES`; mode 0400
  and EPERM remain ineligible.

- The immutable v18 archive proves a signed authority, successful V/R and fresh
  V, then a C child failure before downstream output creation because the live
  tumor-REF `summary.json` is a relative symlink, not a canonical regular path.
- v21 must bind the canonical `all_ssnv_tumor_ref_control_summary.json` bytes and
  descriptor first, then bind `summary.json` as an exact relative basename
  symlink to that target using `O_PATH|O_NOFOLLOW`. Both alias and target identity
  must remain stable for process lifetime and be covered by mutation monitoring.
- Confirm real regressions accept only the exact relative alias and reject an
  unbound target, wrong target, absolute link text, alias replacement, and
  canonical-target replacement. Confirm `GATE_INPUT_PATHS["tumor_ref_summary"]`
  points to the canonical target rather than the symlink.
- Confirm seven terminal keys are pairwise distinct: legacy v2 remains preserved;
  failed-v16 v6, failed-v17 v7, failed-v18 v8, failed-v19 v9, and rejected-v20
  v10 remain mode 0400 but quarantined; only fresh v11 may authorize v21 terminal
  outputs.
- Confirm V/R/C state machines, descriptor leases, mutation watches, waitpid
  witness, no-clobber publication, signal/crash handling, and one-time private
  key retirement remain fail closed.
- Confirm the builder requires exact equality between probe-reported sources
  and subsequently authority-bound sources. Same-UID hostile replacement
  remains outside the declared trusted-account boundary.
- Confirm v1-through-v20 artifacts are never overwritten and all 296 protected
  slots, 20 staging patterns, and 506 tests are bound. Every slot and every
  staging pattern must have an occupied-state fail-closed regression.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains latent molecular substructure candidates only;
  confirmed cellular subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `6cd9bcae1276f5d25ae0048fbb0c4c3f487d79276b3df97272d223cf4720dc82`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `1c63d1d8b12b4f3fa65b5840726bd116464b598b8d74f7a8843e209b03874c05`
- `scope_sha256`: `9b6242aaa6642b30ef62ccbd2ad807ffb6b3a79ff35a3feebe6e169616786ea8`
- `terminal_key_rotation_sha256`: `fc0ee34e846e1a04e73f50dd01ef5c42a6ea7caf134b01777b9d221ba70a149e`
- trusted authority public-key SHA-256: `b66aa240810122b7d9be0ff89840339a7a58563163a55b5d2a837178dffebafb`
- recovery terminal public-key SHA-256: `825111aa7dcd25e60c6357243e228422d2cf07855d682aae30d90ebdcd5559d2`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`ab3ab34d-76f2-4da6-ad2a-af5fb88aee2d`.
