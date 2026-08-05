<!--
建立時間: 2026-07-17
目標: 修正 phase 誤解後，對 v2 source implementation 做 bounded pre-run re-review
處理範圍: source authorization gate；不包含尚未執行的 post-run artifact closure
關聯檔案: 20260717_external_claude_reviewer_A_v5_rejected_scope_mismatch.json
-->

# External Reviewer A v5b: phase-correct pre-run source approval

Perform a fresh bounded read-only re-review. Do not edit, create, chmod, or delete files.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Topic: `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation`

Expected Git HEAD: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`

Expected 23-role source-set digest: `b00a1c3605520af5fcc314d9c55b9a00ef90c651244b8c3666a30a67349a6add`

Canonical test XML: `logs/pytest_full_pre_authority_v7_fixed_signer.xml`, SHA-256 `580fa83b75f5167189f4ea50f3939829e7fe022bc218db214ee07a64985f34b0`, `377 passed`, `0 failed`, mode `0444`.

## Phase definition

This is **P0 pre-run source authorization**, not P1/P2 runtime release closure.

At this phase the following absence is required and must **not** be treated as a blocker:

- v2 authority manifest/approval/signature do not yet exist;
- v2 source private key is mode `0400`, because it is retired to `000` only after this review plus a second independent approval are embedded and signed;
- result/report receipts, signatures and ten report outputs do not yet exist, because formal execution is forbidden until source authorization validates.

Requiring those artifacts now creates a circular dependency: the authority validator requires two reviewer `APPROVE` records before signing, while formal producers require the signed authority before generating runtime outputs. Post-run live-artifact closure will be a separate review after both detached signatures exist.

For this review, `F1_status` and `F2_status` mean the previously identified **source-implementation blockers**:

- F1: producer authority/command/source-lock binding is adequate and fail closed.
- F2: terminal report release implementation binds and rehashes the signed dataset plus ten report artifacts after QA.

They do not mean that runtime outputs already exist.

## Required checks

1. Reconfirm exact source digest, 23 roles, all mode `0444`, HEAD and claim-contract-v5 identity.
2. Confirm all protected trust anchors consistently use v2 source/result/report paths, hashes and authority ID; no protected source accepts v1.
3. Confirm `run_one_time_ed25519_signer_v2.sh` pins `/usr/bin/openssl`, verifies before publication, retains the key on failure, and retires only after successful verification.
4. Confirm the archived v1 files are now named `UNSIGNED_NOT_AUTHORIZED.*`, remain outside all v2 validator paths, and no v1 signature exists.
5. Reconfirm F1/F2 source implementation and their wrong-authority/wrong-code/missing-`-I`/HTML-tamper/QA-tamper tests.

Return only the requested structured result. `APPROVE` means the exact source implementation is safe to embed in the v2 authority and begin the formal run. It is explicitly not approval of absent biological/runtime results.

