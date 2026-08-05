<!--
建立時間: 2026-07-23T06:57:00+08:00
目標: 記錄 v9 authority 第一次 publish 在簽章前因 Python memfd 常數缺失而 fail-closed
處理範圍: schema-recovery authority ceremony；不涉及科學資料或 downstream output
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/build_tumor_ref_schema_recovery_authority_v9.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/schema_recovery_tests/test_tumor_ref_schema_recovery_v9.py
-->

# v9 memfd constants pre-signature failure

- Command: exact clean-environment `build_tumor_ref_schema_recovery_authority_v9.py --publish-after-three-approvals`.
- Exit: `1`; `sealed_memfd()` raised `AttributeError` because this Python runtime does not export `os.MFD_CLOEXEC` or `os.MFD_ALLOW_SEALING`.
- Signing state: failure occurred before `signature_state["retirement_armed"] = True` and before the OpenSSL signing subprocess.
- Key state after failure: v9 private key remained mode `0400`, size 119 bytes; public key remained mode `0444`, size 113 bytes.
- Output state after failure: authority bundle, staging directory, verification/replay/continuation receipts, and downstream outputs remained absent.
- Review state: the three approvals for source-set `61941777c6e3b39c75c9875284c2f9244490055b85d12b43203052623d02efa0` were moved here and are superseded after the source patch.
- Follow-up regression: the active Python runtime also lacks `os.memfd_create`; the new sealed-memfd test rejected the constants-only patch with `190 passed, 1 failed` before any second publication attempt.
- A second follow-up regression showed that the runtime also omits Python `fcntl` seal constants; it again rejected the patch with `190 passed, 1 failed` before publication.
- Fix: retain runtime-exported functions/constants when available, otherwise use the Linux ABI for memfd flags, prebound libc `memfd_create`, `F_ADD_SEALS=1033`, `F_GET_SEALS=1034`, and seal bits `0x1/0x2/0x4/0x8`; execute a sealed-read and kernel-returned complete-seal-mask regression against the active Python runtime.
