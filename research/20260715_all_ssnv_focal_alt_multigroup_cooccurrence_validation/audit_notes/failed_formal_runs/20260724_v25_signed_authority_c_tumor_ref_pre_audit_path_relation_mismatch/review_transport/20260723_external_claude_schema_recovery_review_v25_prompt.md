<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v25 的七個 frozen source 執行外部 Claude Opus 唯讀正式重審
處理範圍: v24 pre-authority rejection、四個 runtime 實際執行、完整 inventory regressions 與 terminal-v15 rotation
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v25.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v25. This is a READ-ONLY remediation re-review after the v24 rejection.
Never edit, write, create, move, chmod, or delete files. Do not use network
tools. Independently inspect the frozen sources and run the exact probe; do not
treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v25.py`: size `461398`, SHA-256 `58d19165741dfec488fa0603dcbb01118ed0f9addddf982abbdeb16b609ed2cc`
2. `verify_tumor_ref_receipt_promotion_recovery_v25.py`: size `129383`, SHA-256 `f5e6f00b936b7bd315ce8652c441d121f6ca0d1e5c1d679f8e946efc86f062f0`
3. `replay_m2v5_runner_only_gates_recovery_v25.py`: size `153508`, SHA-256 `871c675a6c5967b0593212ffb6b5db5749879a7ef7b821f94eaa0f09b2e1ecfb`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v25.py`: size `392379`, SHA-256 `151d6810a4f3dc27a50ba5e99aa23ad1f297ce9fb9de2053b816f6b55f2e8e66`
5. `probe_tumor_ref_schema_recovery_sources_v25.py`: size `119153`, SHA-256 `c73f67e2dfae33770165b691702332e2fed73fb497fbfdd381aae056d798e1fb`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v25.py`: size `174699`, SHA-256 `d5f753e85269f55477a8680a3b15581219560b41e8ba1e9a0aacb86feb2284e4`
7. `build_tumor_ref_schema_recovery_authority_v25.py`: size `61065`, SHA-256 `9da7c2a2b221450908685270924a1ab65c14a0dda6df096664c9550dd40fde7c`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v25.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `597 passed`,
364 forbidden slots absent, 24 staging patterns absent, seven current frozen
sources, `review_evidence_state=all_absent`, and rejected pre-authority
generations including v23 and v24 before review publication.

Audit these decisive points:

- The immutable v24 rejection archive must bind both internal
  `REQUEST_CHANGES` reviews, the external approval envelope, all seven v24
  source identities, exact unused authority-v24 and terminal-v14 keys, and all
  absent formal output slots. Strictest-review-wins must have prevented every
  v24 review/authority/V/R/C output. Those unused mode-0400 keys must remain
  quarantined and must never be reused.
- Confirm v25 closes every v24 finding. Regression tests must execute the
  descriptor-bound primary Python, QA Python, canonical Node v22.22.1, and QA
  Chromium binaries. They must fail closed when Node is missing or Chromium is
  substituted. Production continuation validation must perform the same
  descriptor-bound Node and Chromium execution checks before downstream work.
- Occupied-state regressions must cover every one of the 364 forbidden output
  slots and every one of the 24 staging patterns. `RECOVERY_SCOPE` must include
  both rejected v23 and rejected v24, and `AUTHORITY_CHECKS` must make exact,
  supportable 364-slot and 24-pattern assertions.
- The exact canonical builder source must remain bound before validator load;
  `__file__`, `sys.argv[0]`, `/proc/self/cmdline`, path, and inode must agree.
  A copied builder must fail before any output. Probe and bound regression-test
  execution contracts must remain self-enforcing.
- Review attribution is deliberately limited: reviewer IDs are
  orchestrator-recorded transport identifiers, not reviewer-held
  cryptographic signatures and not independent proof of reviewer authorship.
- The immutable signed v21 archive and every rejected/failed successor must
  remain continuously watched. Historical relation failures remain evidence;
  they cannot be silently normalized into success.
- Confirm eleven terminal public keys are pairwise distinct: legacy v2 remains
  preserved; failed-v16 v6, failed-v17 v7, failed-v18 v8, failed-v19 v9,
  rejected-v20 v10, failed-v21 v11, rejected-v22 v12, rejected-v23 v13, and
  rejected-v24 v14 remain quarantined; only fresh v15 may authorize v25
  terminal outputs.
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

- `reviewed_source_set_sha256`: `262c62a2f990a393961e1b86c279011b3a3df9ad7f1db74ac8a31dd10fd0950f`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `1c229e9319b9776f92b078cd97e522335b419a0a277a3259444d1fd3abc61738`
- `scope_sha256`: `aff281b87ab190326eb5c658d1dc6020f1b4df470f9e2059b9cb1cd9d8ca7120`
- `terminal_key_rotation_sha256`: `421e188e7aac6b6aff309a97f6986d8a0f93eed39ff358f31626b7e5fb03f6cc`
- trusted authority public-key SHA-256: `8b5dc9f1715a3b8b0dafee138a66135f1a29a8bcc324ae15f28a0fd7f3ea54a6`
- recovery terminal public-key SHA-256: `b0056c3f60d7a7204d782ac1cea31e1e9411200b371295e38b9afc3d5f67a1d1`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. The attribution object records transport only and
must explicitly set `cryptographic_reviewer_authorship_proven=false`. Return
only the JSON object required by the supplied schema. Use reviewer
`External Claude Opus` and reviewer/session ID
`8ddaf10e-14e4-4e03-ae73-92c0c7051e69`.
