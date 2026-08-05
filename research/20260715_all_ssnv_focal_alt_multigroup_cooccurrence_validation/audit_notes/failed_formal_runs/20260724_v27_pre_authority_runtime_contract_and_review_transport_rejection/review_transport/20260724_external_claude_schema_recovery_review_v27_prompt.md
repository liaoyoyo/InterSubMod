<!--
建立時間: 2026-07-24
目標: 對 tumor-REF schema recovery v27 的十個 frozen source 執行外部 Claude Opus 唯讀正式重審
處理範圍: failed-v26 historical/current runtime projection、14-role gate registry、v1/v6 lineage、完整 inventory regressions 與 terminal-v17 rotation
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v27.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v27. This is a READ-ONLY remediation re-review after the signed v26
authority failed because three recovery-only runtime roles contaminated the
historical 11-role signed projection. Never edit, write,
create, move, chmod, or delete files. Do not use network tools. Independently
inspect the frozen sources and run the exact probe; do not treat this prompt as
evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v27.py`: size `497952`, SHA-256 `12d194b97f6dac302d3e5234c858ceff3ff00170e5d99344bc9af8f29f601801`
2. `verify_tumor_ref_receipt_promotion_recovery_v27.py`: size `129377`, SHA-256 `5e83ed0add71ba832bdc6c5b3729caba198d1c1bc9d831426526d1059016af25`
3. `replay_m2v5_runner_only_gates_recovery_v27.py`: size `153502`, SHA-256 `5302d3307048fd24df3bd82015b4729248c3226a057c1cf873e6487569ff4aa9`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v27.py`: size `403619`, SHA-256 `1005bde75fcbc35fbe53457195bda6e55e01a8df52a3974ab61776c69a170b95`
5. `probe_tumor_ref_schema_recovery_sources_v27.py`: size `134418`, SHA-256 `8be86c609904d4d891db872a1339b34d33f3dfab3a61002b60d47f8360ea2a75`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v27.py`: size `194893`, SHA-256 `a9bb4e654643e82758a6687a6c37f0fb323cc480e0d44c1a9f9c8ee1fbd05891`
7. `build_tumor_ref_schema_recovery_authority_v27.py`: size `61065`, SHA-256 `23137d0f174e435c7cdfcd3c9b4b35efd91063104ff06c17ff0df75c6aa8e5c0`
8. `build_all_ssnv_final_report_dataset_schema_recovery_v27.py`: size `351388`, SHA-256 `cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f`
9. `finalize_task_b_result_release_schema_recovery_v27.py`: size `33500`, SHA-256 `35508a5049e5bba2b1e72e8b08a2e9d0f35665cfb9208a082549923d6a2b3fc7`
10. `build_all_ssnv_report_artifact_schema_recovery_v27.py`: size `238722`, SHA-256 `f16cc51bbec3ee23e7598f8ded3eb50f85e13e1634d5c795c8335024cbfeac32`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v27.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `656 passed`,
390 forbidden slots absent, 26 staging patterns absent, ten current frozen
sources, `review_evidence_state=all_absent`, and the signed failed-v26 generation
quarantined before review publication.

Audit these decisive points:

- The immutable failed-v26 archive must preserve the signed authority, signature,
  commit, terminal-v15 key state, exact failure evidence, and every protected
  formal/runtime output state. The v26 private authority and terminal keys must
  remain retired or quarantined and must never authorize v27 output.
- Confirm v27 represents the immutable tumor-REF v1 pre-audit and current v6
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
- Occupied-state regressions must cover every one of the 390 forbidden output
  slots and every one of the 26 staging patterns. Source, probe, and bound test
  execution contracts must remain self-enforcing.
- The exact canonical builder source must remain bound before validator load;
  `__file__`, `sys.argv[0]`, `/proc/self/cmdline`, path, and inode must agree.
  A copied builder must fail before any output.
- Review attribution is deliberately limited: reviewer IDs are
  orchestrator-recorded transport identifiers, not reviewer-held
  cryptographic signatures and not independent proof of reviewer authorship.
- Confirm all thirteen terminal public keys are pairwise distinct: legacy v2;
  failed-v16 v6, failed-v17 v7, failed-v18 v8, failed-v19 v9, rejected-v20 v10,
  failed-v21 v11, rejected-v22 v12, rejected-v23 v13, rejected-v24 v14, and
  failed-v25 v15 and failed-v26 v16 remain historical/quarantined; only fresh
  v17 may authorize v27 terminal outputs.
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

- `reviewed_source_set_sha256`: `ef396448143c3de86bde18c97f074a8463e24b129e1d2b08d6ebc9332f9b866e`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `1c229e9319b9776f92b078cd97e522335b419a0a277a3259444d1fd3abc61738`
- `scope_sha256`: `bc330cea21146888cd695316792ba0a89784ef73d9e6de0622245c9f29bfae4a`
- `terminal_key_rotation_sha256`: `9052acf5ae91d596f2cdfb9d42448d535eb7b892eab1b701afd37a719a19e0e7`
- trusted authority public-key path: `/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v27/ed25519_public.pem`
- trusted authority public-key SHA-256: `db8c8b5c5ea1a7e4235efb59966c6d41981ed9b6085af6226c07226b56ab3cb3`
- recovery terminal public-key SHA-256: `355979601138f8dca29534db58e3862f30f30b49526c01d10a97b01ba91c26f9`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. The attribution object records transport only and
must explicitly set `cryptographic_reviewer_authorship_proven=false`. Return
only the JSON object required by the supplied schema. Use reviewer
`External Claude Opus` and reviewer/session ID
`949fbc0f-1af7-4186-9ac6-eefbeb7e76be`.
