<!--
建立時間: 2026-07-24
目標: 封存 v26 已簽 authority 的 C-stage signed-promotion runtime-role projection 失敗
處理範圍: v26 authority/reviews/V/R、C fail-closed、authority/terminal key quarantine
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260724_v26_signed_authority_c_signed_promotion_runtime_role_projection_mismatch/failure_evidence.json
-->

# v26 Signed C Failure Archive

- v26 的三方 review、signed authority、V 與 R 均通過；C 在任何 downstream producer 執行前 fail-closed。
- 直接原因：歷史 signed promotion authorization 只綁定 11 個 runtime roles；v26 將 3 個新的 recovery builders
  加入同一個 `REVIEWED_RUNTIME_SOURCE_CONTRACTS`，使歷史 authorization 的 expected projection 多出三個未曾簽入的角色。
- 多出的角色為 `recovery_final_dataset_builder`、`recovery_report_builder`、
  `recovery_result_report_finalizer`；這是 provenance projection 分層錯誤，不是資料或科學計算差異。
- 沒有 strict/matched-normal/CN/final dataset/report producer 輸出，沒有 final dataset/report/supplemental receipt 或 signature。
- authority key 已在單次 authority 簽署後退役；未使用的 terminal-v16 key 一併封存並禁止重用。
- v27 必須分開歷史 signed-authorization runtime projection 與 current downstream reviewed-runtime set，並新增回歸測試。

