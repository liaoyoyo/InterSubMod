<!--
建立時間: 2026-07-24
目標: 對 tumor-REF schema recovery v29 與 v28 signed-dataset/report-failure 修復執行外部 Claude Opus 唯讀正式重審
處理範圍: 10-source frozen set、v28 742-item archive、v29-round1 rejection archive、JSON mapping key-set 修復、424-slot/28-pattern regressions、15-key terminal state、四角色 key separation
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

1. `validate_tumor_ref_schema_recovery_authority_v29.py`: size `558733`, SHA-256 `2cad3ab1b5b6cf59f340d5920aa9d76e68dda52c39e56159a81ae9335c6df629`
2. `verify_tumor_ref_receipt_promotion_recovery_v29.py`: size `129377`, SHA-256 `fac6e873039f7c5c30c8695ec75e1da114ce475d75f60431b32486114a891ce5`
3. `replay_m2v5_runner_only_gates_recovery_v29.py`: size `153502`, SHA-256 `1c1b0f2f87bfa406e91237ac7946b71a5130794bb39e383be092aa512a0303cf`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v29.py`: size `412688`, SHA-256 `69fa119b305ad3831b40902569a02562c4e3b08e7b9449f8f99391be7a8be1ab`
5. `probe_tumor_ref_schema_recovery_sources_v29.py`: size `159078`, SHA-256 `20b7bd99609498a386e23ea2aa4b133b0f1ac345ea0f706a8a961c1886bc6bca`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v29.py`: size `202835`, SHA-256 `c016320e9a618219e464284a8b2c1fc230432734fe45a8d91743243a63a166df`
7. `build_tumor_ref_schema_recovery_authority_v29.py`: size `61065`, SHA-256 `fbc6c4a3258f1aa9ca2226fcebce831be7e6d7e58b812acdbebb24775b9f046d`
8. `build_all_ssnv_final_report_dataset_schema_recovery_v29.py`: size `351388`, SHA-256 `cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f`
9. `finalize_task_b_result_release_schema_recovery_v29.py`: size `33424`, SHA-256 `46c093c4d86dcb1ac69693195875f1ba9c75ca6fd251141088f259449b711a9e`
10. `build_all_ssnv_report_artifact_schema_recovery_v29.py`: size `238719`, SHA-256 `fe19be151bfd72978f87f4a003e8fb8732f54c21e49f702ea5c575b92f4ee9ae`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v29.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `700 passed`,
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
- The rejected v29-round1 archive at
  `audit_notes/rejected_pre_authority_reviews/20260724_v29_round1_failed_v28_terminal_binding_gap/`
  must bind exactly 20 payload files plus `rejection_evidence.json` SHA-256
  `6d44a75b32737c0e16a7c05c3801ade734977336982a7f08260b70610442cafa`
  and `SUMMARY.md` SHA-256
  `56122d97ee422d9d80d5508ed0655cd7c7fb9be75887b2b4bc38ac0be8c4fa92`.
  The strictest reproducible finding must remain validator state `30/15`
  versus round1 continuation state `28/14`, with no authority, formal review,
  runtime output, or signature created in that rejected generation.
- v28 authority, terminal-v18, result-v5, and report-v5 keys must all be
  archived and prohibited from reuse. The consumed result-v5 and authority
  private keys must be mode 000; unused terminal-v18 and report-v5 private keys
  must remain quarantined mode 0400. None may authorize v29 output.
- The report repair must compare exact mapping key sets, not insertion order,
  for `per_sample`, `by_truth_label`, and `by_sample_truth_label`. Canonical
  sorted-JSON roundtrip must pass, while missing and extra keys must fail.
  Scientific values and mapping contents must not change.
- Confirm all fifteen terminal public keys are pairwise distinct: legacy v2;
  failed/rejected v16-v28 terminal generations including the archived failed-v28
  terminal-v18 key; and only fresh terminal-v19 may authorize v29 outputs. The
  validator and continuation must expose exact matching 30-record state and
  15-distinct-key projections before promotion.
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

- `reviewed_source_set_sha256`: `eaa21a024bf871e2d6bdbcdb149d728bee91917e1f748036c225f5c5ecc2de69`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `8c929332d4b473c3a7bf6ad45959e131eced81a674b6ce1b806b168d48dee6b5`
- `scope_sha256`: `b614cd41c93c672263d47f6e1503b23f1a29b673b9c93268630dc7638cd17768`
- `terminal_key_rotation_sha256`: `bd50b136fa113530d00c33aa878b965772871fd760444495e3c5a49af4c04a0f`
- trusted authority public-key path: `/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v29/ed25519_public.pem`
- trusted authority public-key SHA-256: `819c0ee9d6c6729ad41f3ea1e8071b90ddf506981c68b6927b8abc98e78e9cda`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`8e6d6be7-d3bd-400d-ae0a-43a7960deec3`.
