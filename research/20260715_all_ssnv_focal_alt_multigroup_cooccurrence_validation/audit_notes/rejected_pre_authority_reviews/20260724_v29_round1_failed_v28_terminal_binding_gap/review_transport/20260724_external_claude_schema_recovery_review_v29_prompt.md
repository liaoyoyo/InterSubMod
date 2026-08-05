<!--
建立時間: 2026-07-24
目標: 對 tumor-REF schema recovery v29 與 v28 signed-dataset/report-failure 修復執行外部 Claude Opus 唯讀正式重審
處理範圍: 10-source frozen set、v28 742-item archive、JSON mapping key-set 修復、424-slot/28-pattern regressions、四角色 key separation
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v29.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v29. This is a READ-ONLY remediation review after v28 passed authority,
V, and R, signed the complete final dataset, and then failed during report build
because JSON mapping insertion order was compared even though canonical JSON is
serialized with sorted keys. Never edit, write, create, move, chmod, or delete
files. Do not use network tools. Independently inspect the frozen sources and
run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v29.py`: size `544041`, SHA-256 `02d84a4c766b86276073e23c6dd2cca8a082103b2daf631babba6fa6252568a7`
2. `verify_tumor_ref_receipt_promotion_recovery_v29.py`: size `129377`, SHA-256 `0eb282df7d8fa5a3b2694c9d71a8db3612dd487adf3960d5b892dfd0507a9861`
3. `replay_m2v5_runner_only_gates_recovery_v29.py`: size `153502`, SHA-256 `60074e0c91f2a95e73cc29e92087b21e87d90db3c7d55421bba5d450473b7563`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v29.py`: size `407270`, SHA-256 `e6c4180a8c9fb3e7a3dd87e767c9cc6f12b6e76706ba75edfd4c71bec7094294`
5. `probe_tumor_ref_schema_recovery_sources_v29.py`: size `154320`, SHA-256 `ffb1827b98a3d24d6821a52fb8e568b2f79d25b3738823d3f4cb1b676379452f`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v29.py`: size `201553`, SHA-256 `027b6498dd3ea81f66bdbd55223774ff5947c7b0aca3ab4c3c84e46d304a3d76`
7. `build_tumor_ref_schema_recovery_authority_v29.py`: size `61065`, SHA-256 `9a4e22ed79814c5b0b696dff5504d9feb29f24c459943458d9008ec12d2a5be9`
8. `build_all_ssnv_final_report_dataset_schema_recovery_v29.py`: size `351388`, SHA-256 `cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f`
9. `finalize_task_b_result_release_schema_recovery_v29.py`: size `33424`, SHA-256 `46c093c4d86dcb1ac69693195875f1ba9c75ca6fd251141088f259449b711a9e`
10. `build_all_ssnv_report_artifact_schema_recovery_v29.py`: size `238719`, SHA-256 `fe19be151bfd72978f87f4a003e8fb8732f54c21e49f702ea5c575b92f4ee9ae`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v29.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `699 passed`,
424 forbidden slots absent, 28 staging patterns absent, ten frozen current
sources, `review_evidence_state=all_absent`, and v28 explicitly included in the
verified failed signed-generation chain.

Audit these decisive points:

- The v28 archive at
  `audit_notes/failed_formal_runs/20260724_v28_signed_dataset_c_report_metric_key_order_mismatch/`
  must bind all 742 inventory rows, evidence SHA-256
  `a49fd9339448f75834e18b71f6b43a257ab99aaf26d833f8e9a2644d4e65fbcb`,
  summary SHA-256
  `a4b98ed12d2f7566e8e9eed293d0090a3fb9d269bd986f02e87670053c7453c8`,
  the signed authority, V/R receipts, incident record, complete provisional
  dataset, and both authority and dataset signatures. The failure must be
  classified as report-build only with `scientific_payload_changed=false`.
- v28 authority, terminal-v18, result-v5, and report-v5 keys must all be
  archived and prohibited from reuse. The consumed result-v5 and authority
  private keys must be mode 000; unused terminal-v18 and report-v5 private keys
  must remain quarantined mode 0400. None may authorize v29 output.
- The report repair must compare exact mapping key sets, not insertion order,
  for `per_sample`, `by_truth_label`, and `by_sample_truth_label`. Canonical
  sorted-JSON roundtrip must pass, while missing and extra keys must fail.
  Scientific values and mapping contents must not change.
- Confirm all fifteen terminal public keys are pairwise distinct: legacy v2;
  failed/rejected v16-v28 terminal generations through archived v18; and only
  fresh terminal-v19 may authorize v29 terminal outputs.
- Confirm four active v29 roles use distinct public keys and correct paths:
  authority `819c0ee9d6c6729ad41f3ea1e8071b90ddf506981c68b6927b8abc98e78e9cda`,
  terminal `04d6acab01d56b0bfe25726a242904afd38bbc3ee47d4e3db29e9eb154e23e8b`,
  result `84438d0a91200108ee06ad7600a3c5428804f37567c351ba33036843ae837864`,
  and report `79a684d855ee2d0010691c2a42439389d0e9148d0b84157e04b322c188df6c59`.
- Confirm the historical signed runtime projection remains an exact 11-role
  allowlist, recovery-only roles remain an explicit disjoint 3-role set, and
  the current reviewed runtime set remains their exact 14-role union.
- Confirm V/R/C sequencing, descriptor and directory leases, mutation watches,
  waitpid witness, no-clobber publication, signal/crash handling, exact relation
  registries, and private-key retirement remain fail closed. Same-UID hostile
  runtime injection remains outside the declared trusted-account boundary.
- The immutable tumor-REF v1 pre-audit and current v6 primary artifacts must be
  represented as a chronological lineage. Their data planes remain exactly
  102,842 sites, 308,526 artifacts, and artifact-set SHA-256
  `195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca`.
- The repair changes governance/runtime validation only. It must not change the
  seven datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS`
  biallelic sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific
  values. The claim ceiling remains latent molecular substructure candidates
  only; confirmed cellular subclones and linear ancestry calls are both zero.
- Review attribution is transport attribution only, not reviewer-held
  cryptographic proof. Set `cryptographic_reviewer_authorship_proven=false`.

Canonical bindings:

- `reviewed_source_set_sha256`: `0087b195ce2b7bfb0495f4f3e2c879b11484a645218a318402d14b47879a019e`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `bbd3ccb4920730ec23f08ddb21550b79b67ea0dd4d48a14908818bd9cd22ffee`
- `scope_sha256`: `a0cec5275fc6bd44502303716c225a8fd5c83ff24273203c7d91e6f3629792d2`
- `terminal_key_rotation_sha256`: `bd50b136fa113530d00c33aa878b965772871fd760444495e3c5a49af4c04a0f`
- trusted authority public-key path: `/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v29/ed25519_public.pem`
- trusted authority public-key SHA-256: `819c0ee9d6c6729ad41f3ea1e8071b90ddf506981c68b6927b8abc98e78e9cda`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`2b0e25af-34e4-4ae4-b171-b41ea89a37b8`.
