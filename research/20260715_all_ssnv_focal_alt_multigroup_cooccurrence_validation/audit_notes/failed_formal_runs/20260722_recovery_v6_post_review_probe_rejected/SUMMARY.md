<!--
建立時間: 2026-07-22T20:57:17Z
目標: 封存 recovery v6 三方審查後 regression 狀態依賴失敗
處理範圍: v6 frozen sources、三份 review、post-review probe、未使用私鑰退役
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/20260722_recovery_v6_formal_rejection_and_key_retirement.v1.json
-->

# Recovery v6 rejection

v6 在 reviews 尚未建立時通過 read-only probe（46/46 tests、70 output slots、12 prior inputs），
Mendel、Nash 與 External Claude Opus 亦都對 frozen source set 回覆 `APPROVE`。三份 review 以
mode `0444` 建立後，mandatory post-review probe fail-closed：
`test_v6_full_builder_reuses_bootstrap_validator_record` 固定預期缺少 review 所造成的
`FileNotFoundError`，但 reviews 齊全後 `build_payloads()` 正常成功，因此 regression 成為
`45 passed, 1 failed`。

這是 state-dependent test lifecycle defect，不是 reviewer verdict 或科學資料失敗。v6 authority
從未建立、簽署從未嘗試、所有 `.recovery.v6` runtime outputs 與 staging 均不存在。三份 review
已搬到本目錄的 `reviews/` 保留原 bytes；原 review slots 留空。v6 一次性私鑰由 `0400` 退役為
`0000`，不得重用。

v7 必須讓相同 regression 在 `review_evidence_state=all_absent` 與 `all_present` 兩個合法狀態都通過，
再重新執行三方審查與 post-review probe。
