<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v23 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: v22 pre-authority rejection、builder/probe/test FD execution binding、terminal-v13 rotation 與科學 claim ceiling
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v23.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v23. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Do not use network tools. Independently inspect the frozen
sources and run the exact probe; do not treat this prompt as evidence.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, link-count-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v23.py`: size `424582`, SHA-256 `9232148c6dadf015658a14312c44f2c1893aa4ce13c54aacab273c9844ccabe1`
2. `verify_tumor_ref_receipt_promotion_recovery_v23.py`: size `129383`, SHA-256 `7c9138d503eb58e4bcb8eb264e6acd28bb8cd02dc2790267c62c306ce8edfdfd`
3. `replay_m2v5_runner_only_gates_recovery_v23.py`: size `153508`, SHA-256 `8e03bb6fdd0b5cc2cd15e9aaaa20daeaa2ecb6baa81fca89c18efbd2baac3a9b`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v23.py`: size `383965`, SHA-256 `327934cbd6827566a32ffa6e54ea44b584140677c40ec9dd596f41737c2f5b9c`
5. `probe_tumor_ref_schema_recovery_sources_v23.py`: size `105946`, SHA-256 `8ea403f2bc0ce60b02614c0af02e80cd57e78329a4a83410bdedcc069d494426`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v23.py`: size `167122`, SHA-256 `de26903ace2627cecdc5b87938be1fd130e8139025c4bd0746d0aa049cf7a17f`
7. `build_tumor_ref_schema_recovery_authority_v23.py`: size `60231`, SHA-256 `f2f97256f1e6419fe98dc5908181765885c679e5fddd06e830d3876fa528e115`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v23.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `558 passed`,
334 forbidden slots absent, 22 staging patterns absent, seven current frozen
sources, and `review_evidence_state=all_absent` before review publication.

Audit these decisive points:

- The immutable v22 rejection archive must bind both internal
  `REQUEST_CHANGES` reviews, the external approval envelope, all seven v22
  source identities, exact unused authority-v22 and terminal-v12 keys, and all
  absent formal output slots. Strictest-review-wins must have prevented every
  v22 review/authority/V/R/C output. Those unused mode-0400 keys must be
  quarantined and must never be reused.
- Confirm v23 closes v22 H1: the exact canonical builder source must be bound
  before validator load; `__file__`, `sys.argv[0]`, `/proc/self/cmdline`, path,
  and inode must agree. A copied builder must fail before any output. The
  authority must bind the same builder record observed by the probe/source set.
- Confirm v23 closes the probe/test execution gap. The probe executes from a
  bound source FD through a bound Python FD with canonical argv0. Regression
  tests are read from a bound FD, hash/inode/mode checked, compiled into the
  exact importlib-mode dotted module name, preloaded, and then collected by
  pytest. The test-side marker and cmdline assertion must prove the executed
  module came from the bound FD; string-only assertions are insufficient.
- Confirm review attribution is deliberately limited: reviewer IDs are
  orchestrator-recorded transport identifiers, not reviewer-held
  cryptographic signatures and not independent proof of reviewer authorship.
  No report or authority may overclaim stronger provenance.
- The immutable signed v21 archive must bind authority/signature/commit,
  approving reviews, successful V/R, replay witness, C incident, exact source
  and key identities, and every observed downstream artifact. All original v21
  outputs, including the eight observed C outputs, must remain continuously
  watched among the 334 absent slots. The archived C failure remains the
  canonical-launcher relation mismatch; science and claim ceiling are unchanged.
- Confirm nine terminal public keys are pairwise distinct: legacy v2 remains
  preserved; failed-v16 v6, failed-v17 v7, failed-v18 v8, failed-v19 v9,
  rejected-v20 v10, failed-v21 v11, and rejected-v22 v12 remain quarantined;
  only fresh v13 may authorize v23 terminal outputs.
- Confirm V/R/C state machines, descriptor and directory leases, mutation
  watches, waitpid witness, no-clobber publication, signal/crash handling,
  exact relation registries, and one-time private-key retirement remain fail
  closed. Same-UID hostile runtime injection remains outside the declared
  trusted-account boundary.
- The repair changes provenance/control flow only. It must not change the seven
  datasets, 469,849 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic
  sSNVs, latest sidecar HP/PS read tags, BAM identities, or scientific values.
  The claim ceiling remains latent molecular substructure candidates only;
  confirmed cellular subclones and linear ancestry calls are both zero.

Canonical bindings:

- `reviewed_source_set_sha256`: `7304704c23a5cf73c6a4f83f28345a7999a57eddfe507b350f1e894affab25fd`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `a70d62fc189469a0eeece4aacc8c8d8c92aebcb0e4c7e65694e0ddc97c9e710d`
- `scope_sha256`: `1cd023602ab9f794de0a0a69452a3a3f0e55a122474ef53c94f76e84393d512d`
- `terminal_key_rotation_sha256`: `c4c37c5756b5158b363dc584cd5cebffa30edacb89cf7df893b77c8e2acdd5c9`
- trusted authority public-key SHA-256: `468fb1526e08b4e2cee44ee27d3cd20438beec2298769fca2a413e164759fbe1`
- recovery terminal public-key SHA-256: `d050c0dfea29b469e271967ea7759d4cbb36c3cbec493010163be5938e81e54a`

Return `APPROVE` only if the exact probe passes without protected writes and
`high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`.
Otherwise return `REQUEST_CHANGES` and `pass=false`. Low findings must be
demonstrably nonblocking. The attribution object records transport only and
must explicitly set `cryptographic_reviewer_authorship_proven=false`. Return
only the JSON object required by the supplied schema. Use reviewer
`External Claude Opus` and reviewer/session ID
`8ddaf10e-14e4-4e03-ae73-92c0c7051e69`.
