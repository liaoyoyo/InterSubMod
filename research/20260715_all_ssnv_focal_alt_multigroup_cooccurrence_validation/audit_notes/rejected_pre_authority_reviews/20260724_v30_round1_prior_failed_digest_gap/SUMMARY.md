<!--
建立時間: 2026-07-24
目標: 封存 tumor-REF schema recovery v30 round-1 pre-authority 審查拒絕與修正依據
處理範圍: Mendel REQUEST_CHANGES、未完成 external Claude transport、未發布 authority/review/runtime outputs
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v30.py
-->

# v30 round-1 pre-authority rejection

## 結論

Mendel 於 `2026-07-24T04:38:14Z` 回傳 `REQUEST_CHANGES`。round-1 不具發布權限，且沒有建立正式 review、authority、V/R/C runtime output 或 signature。

## 阻斷問題

`prior_failed_signed_recovery` 會由 authority builder 重新計算並簽入 authority，但 round-1 review schema 只綁定 `fresh_key_bootstrap_sha256`，沒有 `prior_failed_signed_recovery_sha256`。因此 reviewer approval 不能證明其檢視的 prior-failure aggregate 與最終 authority 相同。

必要修正：

1. review exact schema 新增 `prior_failed_signed_recovery_sha256`。
2. `_validate_review()` 傳入並 exact 比對 aggregate digest。
3. 三份 review、publisher、prompt/schema 同步新增 binding。
4. canonical 與 sidecar mutation tests 覆蓋缺欄、錯 digest 與 aggregate mutation。
5. `RECOVERY_SCOPE` 的 `fresh_v29_terminal_key_rotation` 更名為 v30/v21 語意。

## 低風險觀察

- authority publish rename 後即解除 signal mask，而 parent `fsync` 與 final identity recheck 位於外層。此項不改變已發布 bundle 的完整性，但需確認 invocation-outcome ambiguity 是否已有可恢復 witness；重審時再判定是否需要擴大修正。

## 外部審查狀態

舊 source set 的 external Claude runner 已主動中止，退出碼 `130`。未建立 envelope/stderr，該 session UUID 無殘留程序。修正與全量測試完成後將使用新 UUID 重新執行。

## 綁定

- round-1 validator SHA-256: `b7eec895b62034ce147eb09a9950436e093e74bb601abb127bf7681b6786ac8a`
- round-1 reviewed source-set SHA-256: `80e9794c97643164b3ebb54a39989161de27658429ada5537e42178e76176174`
- reviewer transport: `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/rejected_pre_authority_reviews/20260724_v30_round1_prior_failed_digest_gap/20260724_tumor_ref_schema_recovery_mendel.v30.round1.multi_agent.transport.json`
