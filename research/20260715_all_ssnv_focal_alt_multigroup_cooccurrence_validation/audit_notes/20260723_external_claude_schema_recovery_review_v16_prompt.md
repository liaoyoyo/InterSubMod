<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v16 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v14 signed failure、v15 pre-authority rejection、v16 relation-state/runtime-provenance 修復與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v16.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v16. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Use only Read, Grep, Glob, and the
explicitly allowed read-only Bash commands. Do not trust this prompt as
evidence: independently inspect the frozen sources and run the exact probe.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v16.py`: size `218492`, SHA-256 `54716f6e3da56f9b47d60972ec2ed42469412a62c789026abd698890fe94171a`
2. `verify_tumor_ref_receipt_promotion_recovery_v16.py`: size `127906`, SHA-256 `750738cef8a8a375f64ed1a8cda90db40e11759e06b962772fc88ba821fb7980`
3. `replay_m2v5_runner_only_gates_recovery_v16.py`: size `152178`, SHA-256 `5586bf8a0050d32c1b9417f73b0ece08baa2dbcdf30fb49d8843f8dcaea2e860`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v16.py`: size `345527`, SHA-256 `17b369b7463d59a68259436fd3de925160758011f5bf000bb97ced709387651e`
5. `probe_tumor_ref_schema_recovery_sources_v16.py`: size `71642`, SHA-256 `1c93d11a7923f60e39b3686c16cd4d18237789e471442e2fa297db47e86cb07c`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v16.py`: size `122670`, SHA-256 `75bc7f4141d9b00ff5e14580bff6d65200a877d5e0f2d8d975fb3bc50b92e518`
7. `build_tumor_ref_schema_recovery_authority_v16.py`: size `52849`, SHA-256 `f80ea6d97d3ea12d0745319714fab06315b3633d8e1788352061eb35d120866d`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v16.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `380 passed`,
221 forbidden slots absent, 15 staging patterns absent, seven current frozen
sources, and the complete immutable rejection/failure chain verified. Before
formal review publication the review evidence state is expected to be
`all_absent`.

Audit these decisive points rather than merely repeating them:

- v14 was formally signed and passed V, then R failed closed with worker exit
  `70`: an exact `{path,mode,size_bytes,state}` private-key relation was wrongly
  treated as byte-readable after the same path was bound metadata-only. Confirm
  v16 binds the immutable v14 failure archive, signatures, reviews, verification
  evidence, stderr, source identities, key lifecycle, and absence of downstream
  outputs.
- v15 was rejected before authority creation even though its external review
  approved: the two independent internal reviews found runtime relation-schema
  and provenance gaps. Confirm the rejection evidence and all v15 source/review/
  raw-envelope/key archives are immutable and bound, both original key roots are
  absent, no v15 authority exists, and v16 cannot reinterpret v15 as authorized.
- In v16 R, SHA-bearing relations remain full-content relations. Only the exact
  four-key `{path,mode,size_bytes,state}` private-key schema is metadata-plus-size.
  It must use `O_PATH` plus `fstat`, validate exact size without reading private
  key bytes, and reject missing/extra keys, wrong mode, boolean size, empty state,
  conflicting state, and full-metadata reclassification in either order.
- Confirm relation declarations are isolated per live `BoundArtifactReader`
  instance through a weak-key registry. Independent readers must not leak state,
  while conflicting declarations inside one reader fail closed.
- Confirm C consumes exact metadata-plus-size through descriptor metadata only,
  never reads private bytes, and keeps the complete V/R/C state machines,
  process-lifetime alias/target watches, supervisor `waitpid` witness,
  no-clobber publication, signal/crash handling, and private-key retirement.
- Confirm the regression suite covers the actual archived v14 key state and all
  v16 order-symmetry/reader-independence/runtime-provenance cases, and that the
  readonly probe itself runs all `380` tests without creating protected output.
- Confirm fresh v16 authority and terminal-v6 keys are used, v1-through-v15
  terminal/recovery artifacts are never overwritten, and all 221 protected slots
  plus 15 staging patterns are fail-closed.
- Same-UID hostile mutation remains outside the declared trusted-account
  boundary. Do not silently broaden the threat model; report false authority
  within the declared crash/interruption/race boundary.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains: methylation may nominate reproducible latent
  molecular substructure beneath a shared ancestral ALT, but one focal locus
  does not prove a cellular subclone or linear ancestry. Confirmed cellular
  subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `ccf71f51b3da602f26e865d1ca5a83ae7b24e11e68556e9410d88ea94d7aef79`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `2747c8271f46a3cf7583f0e9315d1e6dc8a8fc4f64eefa92c89fac88fda2d59c`
- `scope_sha256`: `94c2e5ddad1e43186f1ef4356032cbf1ab8856368d157347c71469c730d235a8`
- `terminal_key_rotation_sha256`: `da4d4d65504a6ddf0e56197e77744f27df9f29279f63b9e3a4d4839356623c30`
- trusted authority public-key SHA-256: `540b64ed3615618efed89637069f772787fdd025acadbbd27b6e334423d2345e`
- recovery terminal public-key SHA-256: `066949c0c36be413cd2fb60670e5a2fbc583ab9a3a4264ebf4d3766aba39867f`

Verdict rule: return `APPROVE` only if the exact probe passes without protected
writes and `high_findings=[]`, `medium_findings=[]`, and
`unresolved_conditions=[]`. Otherwise return `REQUEST_CHANGES` and
`pass=false`. Low findings must be demonstrably nonblocking. Return only the
JSON object required by the supplied schema. Use reviewer
`External Claude Opus` and reviewer/session ID
`981fda5f-bebe-47cd-83b1-7a152421d474`.
