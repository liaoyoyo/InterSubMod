<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v17 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v16 signed failure、v17 legacy eight-key stat 修復、terminal key rotation 與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v17.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v17. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Independently inspect the frozen
sources and run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v17.py`: size `245648`, SHA-256 `28c03f8cddf0315bd51ae6c86b712ec5ecb74c74d81d9808da45c5d48052de5c`
2. `verify_tumor_ref_receipt_promotion_recovery_v17.py`: size `127906`, SHA-256 `53dc3e1f4784b97acd1364542f0482cd177d260a2a20c89feb2f8b55196bf2c3`
3. `replay_m2v5_runner_only_gates_recovery_v17.py`: size `152178`, SHA-256 `45a463c3f9a6828ecaeaed46454bd28ea2511950da803b4bdb8064c3ad6a9e32`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v17.py`: size `350476`, SHA-256 `3aac571e4f0d3f80693f01213fd8c72cf400e2beb4839caf3a5e135851cfc4c4`
5. `probe_tumor_ref_schema_recovery_sources_v17.py`: size `75511`, SHA-256 `085c60e18380281cf02bc71c9602615a2a8a5fa77436a92ea2bd4dd113004bd9`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v17.py`: size `131243`, SHA-256 `78908ae59a316db78b46e7c74ca8440fffeb5573dcf3251bf94d86c015f38aff`
7. `build_tumor_ref_schema_recovery_authority_v17.py`: size `53317`, SHA-256 `1839f11a5f7ff8b212f714720c733023fa56b92cf156474da237d70568776671`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v17.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `403 passed`,
236 forbidden slots absent, 16 staging patterns absent, seven current frozen
sources, and `review_evidence_state=all_absent` before review publication.

Audit these decisive points:

- v16 was formally signed and passed V/R, then C failed before child start with
  `Ambiguous identity-like relation schema`. Confirm the immutable archive
  proves that the two historical focal-source records use exactly eight stat
  keys and omit `link_count`, all v16 live slots are absent, and no scientific
  payload changed.
- Confirm v17 recognizes only the exact legacy eight-key stat schema, validates
  every recorded field without inventing `link_count`, and excludes historical
  relations from the live identity registry. Missing/extra keys, booleans,
  conflicting content, and same-path live downgrades must fail closed.
- Confirm full-content, full-stat, executable-alias, and metadata-plus-size
  schemas retain their prior exact behavior and per-reader state isolation.
- Confirm terminal keys are pairwise distinct: legacy v2 remains preserved,
  failed-v16 terminal-v6 remains mode 0400 but is quarantined and forbidden from
  reuse, and only fresh terminal-v7 may authorize v17 terminal outputs.
- Confirm V/R/C state machines, descriptor leases, mutation watches, waitpid
  witness, no-clobber publication, signal/crash handling, and one-time private
  key retirement remain fail closed.
- Confirm the builder requires exact equality between the probe-reported seven
  source records and the subsequently authority-bound source records. Same-UID
  hostile replacement remains outside the declared trusted-account boundary.
- Confirm v1-through-v16 recovery artifacts are never overwritten, all 236
  protected slots and 16 staging patterns are guarded, and all 403 tests cover
  the actual v16 archived payload plus schema near-misses and both declaration
  orders.
- Same-UID hostile mutation remains outside the declared trusted-account
  boundary. Do not silently broaden the threat model.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains: methylation may nominate reproducible latent
  molecular substructure beneath a shared ancestral ALT, but one focal locus
  does not prove a cellular subclone or linear ancestry. Confirmed cellular
  subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `6da641e2845dbb2b46d8e47526abe4ffe381c0b166fd71f519c18a60f0aaad11`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `2747c8271f46a3cf7583f0e9315d1e6dc8a8fc4f64eefa92c89fac88fda2d59c`
- `scope_sha256`: `c52cfd78b7e7c5e536aa75ae2662c0b425afa0d7a965659c8ec792ae9f618bd8`
- `terminal_key_rotation_sha256`: `855c506ae9b689aec3fa0af4820ef0a44890eb70caf6eb6f7abb8cf42faf479f`
- trusted authority public-key SHA-256: `ef3413162e1cb9850fb966a1505d71cf185fe79812415b27a93f4d5887154eac`
- recovery terminal public-key SHA-256: `4b83e655f1a7a778691e27dc3df2257a230001c702afdf6e703e578c706e0b03`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`c0b814e1-1046-49d0-864d-ea37f933c35b`.
