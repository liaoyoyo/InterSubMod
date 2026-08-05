<!--
建立時間: 2026-07-18
目標: 補做 v7 遺漏的 live 23-source digest 重算並重新給獨立 verdict
處理範圍: P0 source authority；只讀、最多 10 次工具呼叫
關聯檔案: 20260718_external_claude_reviewer_A_v7_rejected_missing_live_digest_recompute.json
-->

# External Reviewer A v7b: required live source-set recomputation

The preceding v7 response is rejected because it did not recompute the specific
live 23-role source-set digest. Perform a fresh review with a fresh UUID. Do not
reuse the prior verdict.

Repository:
`/big7_disk/liaoyoyo2001/InterSubMod`

Authority module:
`/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py`

Required checks:

1. Load the live module by absolute path with Python `importlib`.
2. For every item in `EXPECTED_SOURCE_PATHS`, independently read the regular
   file, calculate absolute path, byte size, SHA-256, and permission mode.
3. Require exactly 23 unique roles and every mode equal to `0o444`.
4. Reimplement the module's digest independently: sort by role, join
   `role|path|size_bytes|sha256|mode` lines with one real newline byte and no
   terminal newline, then SHA-256 the bytes. Do not call the module's
   `source_set_digest` as the only computation.
5. Confirm the independent value equals
   `9542e7f0f1d12794f7b7736e106f927556d1cc95612d4e99576bb8dd33521a17`.
6. Independently run the exact shell pipeline from both wrappers against the
   same generated records and require the same digest.
7. Confirm the three changed protected sources are mode `0444`, Git HEAD is
   `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`, and neither source file identity
   nor mode changes between the first and last observation.

Also recheck the v7 evidence paths:

- Canonical XML:
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v10_newline_digest_fix_canonical.xml`
- Public key:
  `/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/20260718_all_ssnv_v4/ed25519_public.pem`
- Cooccurrence wrapper:
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_cooccurrence_v6_source_locked.sh`
- Completion wrapper:
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_m2v5_recovered_completion_chain.sh`

The v4 authority JSON is expected to be absent and is not needed for this live
recomputation. Do not edit/create/chmod/delete/sign, run pytest, launch
producers, use network, or use recursive search. Use at most 10 tool calls.

Return only one valid JSON object with exactly these top-level fields:
`schema_name`, `schema_version`, `reviewer_label`, `reviewer_id`, `model`,
`verdict`, `findings_closed`, `f1_status`, `f2_status`,
`reviewed_source_set_sha256`, `reviewed_git_head`, `review_scope`,
`blocking_findings`, `nonblocking_findings`, `evidence`.

Use `schema_name=intersubmod.external_claude_source_review`,
`schema_version=1.0.0`, a fresh UUID, and `APPROVE` only if every required live
digest check passes.
