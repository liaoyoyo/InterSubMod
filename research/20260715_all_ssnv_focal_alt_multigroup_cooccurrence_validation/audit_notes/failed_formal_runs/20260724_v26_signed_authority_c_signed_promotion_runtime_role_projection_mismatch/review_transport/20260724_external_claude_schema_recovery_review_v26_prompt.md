<!--
建立時間: 2026-07-24
目標: 對 tumor-REF schema recovery v26 的十個 frozen source 執行外部 Claude Opus 唯讀正式重審
處理範圍: failed-v25 signed authority relation、v1/v6 lineage、完整 inventory regressions 與 terminal-v16 rotation
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v26.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v26. This is a READ-ONLY remediation re-review after the signed v25
authority failed its tumor-REF pre-audit path relation. Never edit, write,
create, move, chmod, or delete files. Do not use network tools. Independently
inspect the frozen sources and run the exact probe; do not treat this prompt as
evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v26.py`: size `476922`, SHA-256 `ead51f4b3134268f7ef163ede38827d38d1f07383e2d7c44f949d971b2c315ab`
2. `verify_tumor_ref_receipt_promotion_recovery_v26.py`: size `129377`, SHA-256 `e32a0d57d49af05a123ae0225f8b6a39fd221e9efccfadf96445a62d08e9b885`
3. `replay_m2v5_runner_only_gates_recovery_v26.py`: size `153502`, SHA-256 `2f19d2cdcc96c65ed96c2d0aacbe7fcc24cc303adc01f70d77b9e3150b72355b`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v26.py`: size `398596`, SHA-256 `a79d3630570f523068b1220b20aa78369b4ba7b665672a6b01fcbcdb47eedd54`
5. `probe_tumor_ref_schema_recovery_sources_v26.py`: size `126050`, SHA-256 `4fcf5df25224d2e18a663bc846957cc4c43913162e70d792cef85d7fbe25b187`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v26.py`: size `188158`, SHA-256 `1670c508bb0ba2be3b5e8439b0465b54c8415a240f1d4a5b928a443a76a87f99`
7. `build_tumor_ref_schema_recovery_authority_v26.py`: size `61065`, SHA-256 `4ab964b97b63e2903c9c3e431545942743768b5485fda0fbeaabd686f29cf840`
8. `build_all_ssnv_final_report_dataset_schema_recovery_v26.py`: size `351388`, SHA-256 `cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f`
9. `finalize_task_b_result_release_schema_recovery_v26.py`: size `33500`, SHA-256 `ab758306b80c26c17c169285273f99e235ba2c16496cd1b9a2a8fdce974f73d1`
10. `build_all_ssnv_report_artifact_schema_recovery_v26.py`: size `238722`, SHA-256 `2bf11448db06c7729e94bd28bfee2f2d83269728eedd09dac185a87b10488527`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v26.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `628 passed`,
377 forbidden slots absent, 25 staging patterns absent, ten current frozen
sources, `review_evidence_state=all_absent`, and the signed failed-v25 generation
quarantined before review publication.

Audit these decisive points:

- The immutable failed-v25 archive must preserve the signed authority, signature,
  commit, terminal-v15 key state, exact failure evidence, and every protected
  formal/runtime output state. The v25 private authority and terminal keys must
  remain retired or quarantined and must never authorize v26 output.
- Confirm v26 represents the immutable tumor-REF v1 pre-audit and current v6
  primary artifacts as a chronological lineage rather than pretending they are
  the same path. Their data planes must remain exactly 102,842 sites, 308,526
  artifacts, and artifact-set SHA-256
  `195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca`.
- Legacy v1 and current v6 tumor-REF inputs must be opened with `O_NOFOLLOW` and
  parsed from the same descriptor-bound bytes used for their SHA-256 records.
  The continuation gate and sentinel must bind `LEGACY_PRIMARY_PRE`, validate
  the manifest crosslink, and detect mutation before downstream execution.
- The v25 signed commit must be validated as an exact transaction: member set,
  authority and signature digests, retired key records, and bundle identity must
  all match. Mutation regressions must fail closed.
- Regression tests must execute the descriptor-bound primary Python, QA Python,
  canonical Node v22.22.1, and QA Chromium binaries. Production continuation
  validation must perform the same descriptor-bound checks before downstream
  work and use the canonical v5 runtime output slots.
- Occupied-state regressions must cover every one of the 377 forbidden output
  slots and every one of the 25 staging patterns. Source, probe, and bound test
  execution contracts must remain self-enforcing.
- The exact canonical builder source must remain bound before validator load;
  `__file__`, `sys.argv[0]`, `/proc/self/cmdline`, path, and inode must agree.
  A copied builder must fail before any output.
- Review attribution is deliberately limited: reviewer IDs are
  orchestrator-recorded transport identifiers, not reviewer-held
  cryptographic signatures and not independent proof of reviewer authorship.
- Confirm all twelve terminal public keys are pairwise distinct: legacy v2;
  failed-v16 v6, failed-v17 v7, failed-v18 v8, failed-v19 v9, rejected-v20 v10,
  failed-v21 v11, rejected-v22 v12, rejected-v23 v13, rejected-v24 v14, and
  failed-v25 v15 remain historical/quarantined; only fresh v16 may authorize
  v26 terminal outputs.
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

- `reviewed_source_set_sha256`: `d08b24b49817563e3b6a806dd91125ef17af24bad32085ec1365af87bdcb3e3e`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `1c229e9319b9776f92b078cd97e522335b419a0a277a3259444d1fd3abc61738`
- `scope_sha256`: `e7173fbec7ee5893c4aca7908bade5685c2f67f6abf1988edd2eb51e75135a83`
- `terminal_key_rotation_sha256`: `4c76c30a41ddc6068ad180ef4b7ea5e1fc7185c48ff25d58c7350817e2c91fe7`
- trusted authority public-key path: `/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v26/ed25519_public.pem`
- trusted authority public-key SHA-256: `36fde24fb209c62465d74b4f42a70ca2361bab3c0e326cac25539e0f42049a41`
- recovery terminal public-key SHA-256: `b61e7f75ba3e418098c3c7f9fb19da0b261380d7985833eb6a59b99d6e1aeaee`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. The attribution object records transport only and
must explicitly set `cryptographic_reviewer_authorship_proven=false`. Return
only the JSON object required by the supplied schema. Use reviewer
`External Claude Opus` and reviewer/session ID
`b20db30c-669d-4c6d-b606-6ccc98d761e8`.
