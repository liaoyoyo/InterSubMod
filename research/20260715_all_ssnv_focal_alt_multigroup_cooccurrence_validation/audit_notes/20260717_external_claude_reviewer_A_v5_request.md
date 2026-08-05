<!--
建立時間: 2026-07-17
目標: 外部 Claude Code 對 v2 一次性簽章修正及完整 Task-B release security 做唯讀複核
處理範圍: 23 protected sources、source/result/report v2 trust anchors、F1/F2
關聯檔案: ../logs/pytest_full_pre_authority_v7_fixed_signer.xml
-->

# External Reviewer A v5: v2 release security and provenance

Perform a fresh read-only audit. Do not edit, create, chmod, or delete any file.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Topic: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Expected Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`

Expected exact 23-role source-set digest: `b00a1c3605520af5fcc314d9c55b9a00ef90c651244b8c3666a30a67349a6add`

Canonical test XML: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v7_fixed_signer.xml` (82,686 bytes, SHA-256 `580fa83b75f5167189f4ea50f3939829e7fe022bc218db214ee07a64985f34b0`, mode `0o444`; `377 passed`, `0 failed`, 38 deprecation warnings).

Independently recompute the digest from `release_source_authority.EXPECTED_SOURCE_PATHS`; confirm exactly 23 roles, all sources mode `0o444`, and unchanged claim-contract-v5 identity: 9,144 bytes, SHA-256 `da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af`, mode `0o444`.

Audit the v2 signing repair. The prior v1 source signer inherited an obsolete Conda OpenSSL and failed before creating any signature; its unsigned manifest/approval are archived and must never validate. Confirm:

1. `release_source_authority.py`, both formal runners and `finalize_task_b_result_release.py` consistently bind only v2 source/result/report public keys, hashes and authority ID; no protected source retains a v1 trust-anchor path.
2. `audit_notes/run_one_time_ed25519_signer_v2.sh` pins `/usr/bin/openssl`, signs to a pending path, verifies before publication, retains the encrypted key after sign/verify failure, and retires it to mode `000` only after successful verification.
3. Source public key SHA-256 is `8773f88579b01a97fe4e39682b42ca5ae55f74872d31b5a5f8585e8942644970`; result is `54598225ab57a52393fbe63a29b24a19f39998bdf5d951fb61f4edab67bfeb24`; report is `98e7ca01a67ce73ac3eea0a18599db78233fd66104f1922db773396ef85f56fb`; all public keys are mode `0444` and the three v2 private keys are still pre-sign mode `0400`.
4. The v1 authority is explicitly `NOT_SIGNED / NOT_AUTHORIZED` under `docs/archive/2026/07/20260717_failed_all_ssnv_source_authority_v1/`; no v1 detached signature exists.

Re-audit F1 producer binding and F2 terminal report release. Formal producers must validate the same authority, reject injected argv, require exact `/proc/self/cmdline` with `python -I`, lock live source identity/mode, and publish mode-0444 receipts. The independent report receipt/signature must bind the signed dataset release plus all ten report/HTML/QA/screenshot outputs, and rehash each artifact on verification. Check wrong-authority, wrong-code, missing-`-I`, HTML-tamper and QA-tamper regressions.

Look for bypasses, circular trust, signing-order errors, mutable unsigned deliverables, missing path/size/SHA bindings, and post-signature artifacts. This is pre-run source approval, not biological confirmation.

Return only the requested structured result. `APPROVE` only when v2 is internally consistent, F1/F2 are closed, and no release-blocking issue remains.

