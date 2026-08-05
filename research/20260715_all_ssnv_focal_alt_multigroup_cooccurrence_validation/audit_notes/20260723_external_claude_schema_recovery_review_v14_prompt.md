<!--
建立時間: 2026-07-23
目標: 對 tumor-REF schema recovery v14 的七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: exact typed relation、process-lifetime inotify、waitpid success witness、persistent ceremony watch、v13 rejection chain
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v14.py
-->

You are the external formal reviewer for InterSubMod Task-B tumor-REF schema
recovery v14. This review is READ-ONLY. Never edit, write, create, move, chmod,
or delete files. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v14.py`: size `179003`, SHA-256 `1d9480b0268fe17c1d6bfb5467671c76a1623d197835821d453ae75a99e9410c`
2. `verify_tumor_ref_receipt_promotion_recovery_v14.py`: size `127906`, SHA-256 `13b159e3f1fcd961fd693331cd8f6f3a3a0862a663cceb648d755f22cb91bd3a`
3. `replay_m2v5_runner_only_gates_recovery_v14.py`: size `148758`, SHA-256 `14a012a3407bbe8d8b0855545f7060763107400bc4206571d2612b11a61adcb7`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v14.py`: size `344261`, SHA-256 `e254a21b4dc95d6c25085256080c1d66f95dd3b474f47aa07096bd8367ca90c1`
5. `probe_tumor_ref_schema_recovery_sources_v14.py`: size `56196`, SHA-256 `69742abbdb033f7a004e0a780364cb301b96a7ef95019d64392e853e1e521dd8`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v14.py`: size `109952`, SHA-256 `933ad2162b62cf58da9d924aed2ad081a6b9e1d5b5e896a79eac497b3481e860`
7. `build_tumor_ref_schema_recovery_authority_v14.py`: size `53376`, SHA-256 `5d9bb98d0192ef01cf3e85ab239c519958b9ec5667b9887adeb1949dea50c326`

Run this exact probe and confirm zero protected-namespace output writes:

```bash
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v14.py
```

Required result: exit `0`, `pass=true`, `no_output_writes=true`, `329 passed`,
191 forbidden slots absent, seven reviewed sources, 69 prior inputs,
`review_evidence_state=all_absent`, and v13 rejection verified.

Audit these decisive points rather than trusting this prompt:

- The prior v14 candidate source-set `8a1b8ba97d45c83e0cf385887f306edb76525760dd7f51d673ab862061b85d84`
  was withdrawn before formal review publication because its authority scope
  retained stale `24/175/12` cardinalities. Confirm the current scope is bound
  to the live `25/191/13` inventories and that regression coverage enforces it.

- C uses four mutually exclusive exact typed relation schemas, rejects bool-as-int,
  has a per-reader path/schema registry, and requires the three Python alias JSON
  pointers to be exact six-key, strict-identical alias records.
- R retains process-lifetime alias/target inotify monitoring. Its receipt/log are
  provisional; only a supervisor-observed `waitpid` normal exit zero can publish
  the success witness. The witness binds source, command, receipt, log, wait
  status, and has no fallible work after terminal no-replace publication.
- V accepts only all-three-absent or receipt+log+witness complete state. C refuses
  receipt/log without the supervisor witness.
- The authority builder retains a single persistent inotify watch from the final
  two-pass absence scan through the last pre-`renameat2(RENAME_NOREPLACE)` check,
  so a same-generation late mutation is still rejected.
- v13 remains immutable and rejected pre-authority. Its two unused keypairs are
  archived never-signed, its formal slots remain absent, and v14 uses fresh
  authority and terminal keys.
- Check key retirement, descriptor/source binding, exact JSON schemas, signal and
  crash windows, partial output states, no-clobber publication, and test coverage.
- Same-UID malicious process mutation is outside the declared trusted-account
  boundary; do not silently broaden the threat model. Look for false authority
  within the declared crash/interruption/race boundary.
- The recovery changes provenance/control flow only. They must not change the
  seven datasets, 469,849 same-run LongPhase-S recalibrated FILTER=PASS sSNVs,
  BAM/read-tag inputs, scientific results, or claim ceiling.

Canonical bindings:

- `reviewed_source_set_sha256`: `d3eb22fa1e565beb5e7888d92d9e440439a277e99b93fc4dfd228492c6ae0feb`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `957cb61147ef9cd0106437b5649b20b5933fedc9a9dfacc8b10de48f37a89966`
- `scope_sha256`: `54f6860820a8d75dab79f677e2e09bce9baedb8f2d7642ac56a2772912d18bd1`
- `terminal_key_rotation_sha256`: `5a11f5e24e58af751886fef1253db09bb1fe83e87f2d54e50d4b945564ca1e12`
- trusted public-key SHA-256: `91f3b81a0dfab1911b492269dc40ef150eed76bf3b16c4143ba541d16ffdc8a3`

Verdict rule: `APPROVE` only if the exact probe passes without protected writes
and `high_findings=[]`, `medium_findings=[]`, and `unresolved_conditions=[]`;
otherwise `REQUEST_CHANGES`. Low findings must be demonstrably nonblocking.

Return only one JSON object with these exact keys:
`reviewer`, `reviewer_agent_id`, `verdict`, `probe_exit_code`,
`probe_no_output_writes`, `probe_regression_summary`,
`probe_forbidden_output_slots_checked`, `reviewed_source_set_sha256`,
`high_findings`, `medium_findings`, `low_findings`, `unresolved_conditions`,
`summary`, and `pass`.

Use reviewer `External Claude Opus` and reviewer/session ID
`828d4dcb-1aed-4ade-bbf8-433ee06a895e`.
