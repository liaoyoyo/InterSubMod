<!--
建立時間: 2026-07-24
目標: 對 tumor-REF schema recovery v30 四角色重新 bootstrap、失敗 key 封存與完整 V/R/C 發布契約執行外部 Claude Opus 唯讀正式審查
處理範圍: round-1 rejection remediation、11-source frozen set、v30 bootstrap、pre-authority 與 failed-v8 archives、439-slot/29-pattern regressions、terminal-v21 rotation、result/report signer readiness
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v30.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v30. This is a READ-ONLY comprehensive-validation review after v29
passed authority and V but failed at R because the continuation still referenced
the archived v29 terminal key through a live path. v30 separates prior failed
signed recoveries from a fresh four-role bootstrap and binds new authority,
terminal, result, and report keys. An internal round-1 review then rejected the
candidate because formal reviews did not digest-bind the separately signed
`prior_failed_signed_recovery` aggregate. Round 2 adds that exact binding, three
mutation canaries, corrects the terminal-v21 scope name, and keeps final
post-rename durability/identity checks inside the blocked-signal interval. Never
edit, write, create, move, chmod, or delete files. Do not use network tools.
Independently inspect the frozen sources and run the exact probe; do not treat
this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`.
Confirm each is a regular file with actual link count 1 and exact size/SHA-256:

1. `validate_tumor_ref_schema_recovery_authority_v30.py`: size `593141`, SHA-256 `424541b75e6e8e47c2c7dff0783b5c19e877e4a32d71b2d22238f979ef492fbc`
2. `verify_tumor_ref_receipt_promotion_recovery_v30.py`: size `129377`, SHA-256 `07543693e204f5837b9ef6403ee9b9d2a9b61c7da7422d571ec7c3c176f4b204`
3. `replay_m2v5_runner_only_gates_recovery_v30.py`: size `166914`, SHA-256 `50025d9012d51e0df6a3f4d3c3bfecfca182923f2ff78645e0cbb6c6540614fc`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v30.py`: size `429261`, SHA-256 `4124faf9ba7a7ec030f1e263a44c2194231047f79a161b85f4474c2e2c75cb5d`
5. `probe_tumor_ref_schema_recovery_sources_v30.py`: size `175248`, SHA-256 `a741cc70014328a6f6f040f6257a7176f6970917922d759a4b6fa7e2efd53b27`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v30.py`: size `224466`, SHA-256 `93b1a9a16dc1ec84f00a61a35e9f41e8d41e5397561e11e6e49ce6923d48b3ad`
7. `build_tumor_ref_schema_recovery_authority_v30.py`: size `61403`, SHA-256 `8058f4f650e7a4dc30ae5a9f8adad8cef56803be41dc69fd60f89b832b8fb156`
8. `bootstrap_v30_four_role_keys.py`: size `63736`, SHA-256 `e971d5db502ebeadca3df28c71be42875de957c16dea2fbc2caeb6579bca8f4a`
9. `build_all_ssnv_final_report_dataset_schema_recovery_v29.py`: size `351388`, SHA-256 `cce5eab32cc17b3d1a46e57684dc04b3ce038a44caaa01a4d4910201cf39e91f`
10. `finalize_task_b_result_release_schema_recovery_v29.py`: size `33424`, SHA-256 `46c093c4d86dcb1ac69693195875f1ba9c75ca6fd251141088f259449b711a9e`
11. `build_all_ssnv_report_artifact_schema_recovery_v29.py`: size `238719`, SHA-256 `fe19be151bfd72978f87f4a003e8fb8732f54c21e49f702ea5c575b92f4ee9ae`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v30.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `733 passed`,
439 forbidden slots absent, 29 staging patterns absent, eleven frozen current
sources, `review_evidence_state=all_absent`, static contract pass, all recursive
self-checks pass, archived-state runner replay pass, historical pilot pass, and
prior recorded verifier recheck pass.

Audit these decisive points:

- `fresh_key_bootstrap` must be a separate object from
  `prior_failed_signed_recovery`. It must bind prepare SHA-256
  `3fb0232160a517d1b202462d64a73643e5932ed541845d58cca778240460b91b`,
  progress `328776bf5124c6e11d15236891e8a3e393bc05b4d90a3ac7f6b7f5d0d1c46343`,
  receipt `6c86539fd8c15bb3e98299474fcd27ebe426bf06d99366d4ea36f0d9212fe4be`,
  success `3d4d09eb64f561fa6e2c1d5bd544f8e371aa44eb44a5daf47ac53b69843dd898`,
  exactly seven ordered progress events, and `release_authority_granted=false`.
- Formal reviews must bind the exact separately validated 14-generation
  `prior_failed_signed_recovery` aggregate with SHA-256
  `f47258fc4bad5d15b5f6a1fd0032b8e39c7c05e359aafaa59a8a35d24936cdda`.
  Confirm the validator and authority builder pass the same aggregate and the
  canonical/sidecar tests reject a missing field, a wrong digest, and a mutated
  aggregate while preserving the review digest.
- The pre-authority archive must bind the four unused initial roots
  (authority-v30-initial, terminal-v20-initial, result-v7-initial,
  report-v7-initial), prove all original roots absent, prohibit key reuse, and
  validate private keys by metadata only without reading bytes.
- The failed result-v8 pre-READY signer must be independently archived at
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_result_release_authority/archive/20260724_all_ssnv_result_v8_failed_pre_ready_parent_nlink_contract_01`.
  Its receipt SHA-256 is
  `19aa68f7ed9d43b86d85d027377857423a03e4c6c4c8d8df744c37a91fd8e6c5`,
  failure-record SHA-256 is
  `d6208b736c5eb40b8f13a79c5ace2e20007beaf489d86e46fb14a441b6599bd2`,
  archived public-key SHA-256 is
  `04d1775e1dacb2cb0222816fb25244f8cdd13cf92eab86fa92d0b172898d918c`,
  source root must be absent, no partial witness may exist, and the private key
  must be inspected by metadata only.
- Confirm all four active v30 roles are pairwise distinct and use exact paths:
  authority `a5b0e0b2c2a9f220d988597b47c8eb1d5446de401a932102948d829ffd0611ed`,
  terminal-v21 `3ea7624ed42caba9bd51ade25a4c9a037f0b84689b4e2a3563c8205bbb136fcd`,
  result-v9 `0d985d9afce029c06275b6932d51db050f807db50ed27c3bf66cd6f9e201f267`,
  and report-v9 `8eaf44b95e216320b45ca4109cf34b2c86b8fd4dcf06bb58025c1e27126bf5b5`.
- The result/report signer prepare and READY witnesses must be exact and must
  not imply that a signature already exists. Result prepare/ready SHA-256 are
  `36e42cc2a4230e974c0da2f1a92b26642444b2449b61a389b136738feb1fbad5`
  and `ab939941c7b887f82c28d384879b90ec6e161b1e52de5ad64c768446835f6751`.
  Report prepare/ready SHA-256 are
  `518160b14aff3220f35e1fff7a73bb32af1c7ac44312fb541c524e6e65d4f3b9`
  and `4aabc90a23688271979e9a4ca33866b982f20d5954d1705bb341cb238b50bbd5`.
- The v29 formal failure archive must remain verifiable while its original
  output slots and four live key roots remain absent. Historical v29 validation
  must use archived historical authority/terminal identities, never substitute
  the current v30 keys. Its failed terminal policy must remain the historical
  authorized policy, while current v30 uses the archive-rebased recovery policy.
- Confirm terminal rotation has exact, pairwise-distinct historical/current
  key projections and only fresh terminal-v21 may authorize recovery-v30
  terminal outputs. No failed/rejected generation key may authorize v30.
- Confirm the historical signed runtime projection remains an exact allowlist,
  recovery-only roles remain explicit and disjoint, and current reviewed runtime
  roles are the exact authorized union. V/R/C sequencing, descriptor and
  directory leases, mutation watches, waitpid witness, no-clobber publication,
  signal/crash handling, exact relation registries, and private-key retirement
  must remain fail closed. Same-UID hostile runtime injection remains outside
  the declared trusted-account boundary.
- Confirm authority publication keeps signals blocked through rename, parent
  `fsync`, final directory identity comparison, and terminal recheck via the
  `post_publish` callback, then restores the prior mask. The existing successful
  publish test must observe `SIGTERM` blocked during the final-name recheck.
- The immutable tumor-REF v1 pre-audit and current v6 primary artifacts must
  remain a chronological lineage. Their data planes remain exactly 102,842
  sites, 308,526 artifacts, and artifact-set SHA-256
  `195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca`.
- The repair changes governance/runtime validation only. It must not change the
  seven datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS`
  biallelic sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific
  values. `persisted_bam=false` for all seven datasets: no tagged BAM overwrote
  the original LongPhase-S BAM. The claim ceiling remains reproducible latent
  molecular substructure candidates in shared ancestral ALT reads only;
  confirmed cellular subclones and linear ancestry calls are both zero.
- Review attribution is transport attribution only, not reviewer-held
  cryptographic proof. Set `cryptographic_reviewer_authorship_proven=false`.

Canonical bindings:

- `reviewed_source_set_sha256`: `d6fe0a64c2346fa2fb3300f5e9c58874022fa4918f3becc7c036e89cf7bbad3f`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `8c929332d4b473c3a7bf6ad45959e131eced81a674b6ce1b806b168d48dee6b5`
- `prior_failed_signed_recovery_sha256`: `f47258fc4bad5d15b5f6a1fd0032b8e39c7c05e359aafaa59a8a35d24936cdda`
- `fresh_key_bootstrap_sha256`: `c9955a2e561f5011f703d09bff0a34cad9a4391fb80137dda804d8154f99ac4d`
- `scope_sha256`: `7d8cc38bbfaa080d29357fa87be2a5c899912d026de877c74c4ec3c1fa610223`
- `terminal_key_rotation_sha256`: `d4da8aba4448ac44cd40bdb0a7ea2b83a4a7bedf53b36595513083094a90f124`
- trusted authority public-key path: `/bip7_disk/liaoyoyo2001/.config/intersubmod_verifier_schema_recovery/20260724_v30_four_role_rebootstrap/ed25519_public.pem`
- trusted authority public-key SHA-256: `a5b0e0b2c2a9f220d988597b47c8eb1d5446de401a932102948d829ffd0611ed`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. Return only the JSON object required by the supplied
schema. Use reviewer `External Claude Opus` and reviewer/session ID
`f321c1c6-eb0e-4aa2-b088-68f0227e714e`.
