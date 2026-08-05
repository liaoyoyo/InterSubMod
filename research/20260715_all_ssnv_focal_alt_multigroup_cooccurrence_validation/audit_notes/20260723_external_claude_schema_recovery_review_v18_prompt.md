<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v18 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v17 signed failure、v18 metadata enrichment 修復、terminal key rotation 與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v18.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v18. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Independently inspect the frozen
sources and run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v18.py`: size `273464`, SHA-256 `269dedf57cdfffebb4849f276986fc4d72a9ecc268df8063f356c3c3e697e77b`
2. `verify_tumor_ref_receipt_promotion_recovery_v18.py`: size `127906`, SHA-256 `0b443b7757b66990f1418562ac375921cdea29a4e2ca4c70d9440fb2b5282e61`
3. `replay_m2v5_runner_only_gates_recovery_v18.py`: size `152178`, SHA-256 `add0447e765b87677578c355c04db0f0db2bf79ad58d6ca9190cd692ead3adb7`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v18.py`: size `354472`, SHA-256 `523c594d6c9edaa9403a674aa0186ed4880b54646d7f675eda63cf2089cabe99`
5. `probe_tumor_ref_schema_recovery_sources_v18.py`: size `79062`, SHA-256 `94a070d0862880ee8e7cf6d0c39c9876a48d4a0e8aca88ba0b513f10aaf9c663`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v18.py`: size `138786`, SHA-256 `23c3bbf99c6f4d916c86575cdcf033ba5cd8e11efd46e4847f16856e86114cc6`
7. `build_tumor_ref_schema_recovery_authority_v18.py`: size `53317`, SHA-256 `05e3c6a0dfa796b4d09c77c1c8dee8f4c79bd26f91a0cae41b33cba6ab09f30d`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v18.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `427 passed`,
251 forbidden slots absent, 17 staging patterns absent, seven current frozen
sources, and `review_evidence_state=all_absent` before review publication.

Audit these decisive points:

- The immutable v17 archive proves a signed authority, successful V/R and fresh
  V, then a C failure before downstream output creation on an exact cross-kind
  metadata enrichment declaration for the terminal-v2 private key.
- v18 may reconcile only exact `metadata` and `metadata_plus_size` declarations
  when shared `path` and `mode` agree. It must retain the richer size/state
  record regardless of declaration order. Conflicting path, mode, size, state,
  unrelated relation-kind mixing, missing keys, extra keys, and booleans must
  fail closed.
- Confirm tests cover the actual archived v17 replay declarations and the
  terminal-v2/v6/v7/v8 private-key metadata in both declaration orders, while
  preserving all prior exact relation schemas and per-reader state isolation.
- Confirm terminal keys are pairwise distinct: legacy v2 remains preserved,
  failed-v16 v6 and failed-v17 v7 remain mode 0400 but quarantined, and only
  fresh v8 may authorize v18 terminal outputs.
- Confirm V/R/C state machines, descriptor leases, mutation watches, waitpid
  witness, no-clobber publication, signal/crash handling, and one-time private
  key retirement remain fail closed.
- Confirm the builder requires exact equality between probe-reported sources
  and subsequently authority-bound sources. Same-UID hostile replacement
  remains outside the declared trusted-account boundary.
- Confirm v1-through-v17 artifacts are never overwritten and all 251 protected
  slots, 17 staging patterns, and 427 tests are bound.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains latent molecular substructure candidates only;
  confirmed cellular subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `ecb6da093f3d669802b42160bebca10c0304d07ca9b3e2c92e31066bab9f835b`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `2747c8271f46a3cf7583f0e9315d1e6dc8a8fc4f64eefa92c89fac88fda2d59c`
- `scope_sha256`: `fb0524dc9f6653759e003c7af163d3b54c5828d6e6d742a58cfcacc59a7cc65f`
- `terminal_key_rotation_sha256`: `bad0cf4ca9bbddf1da439eeda0c76d8c96528d3c9178bc4096d1a2f9fc6fde82`
- trusted authority public-key SHA-256: `57062811b5b67f514f1454cac3e96fe55f7ab94be2fbbbace43431f68fbb65e4`
- recovery terminal public-key SHA-256: `9fa7667e8076cee90a93fee44cfe08bc45a68e9ce28d01507b5c057463102a93`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`bf87dba8-975f-485b-9954-1f28aad9dd42`.
