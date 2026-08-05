<!--
建立時間: 2026-07-24
目標: 封存 v28 已簽 dataset 後的 report metric key-order schema 失敗
處理範圍: v28 authority/reviews/V/R/C incident/provisional dataset/key quarantine
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260724_v28_signed_dataset_c_report_metric_key_order_mismatch/failure_evidence.json
-->

# v28 Signed Dataset C Failure Archive

- v28 authority、V、R 與 final dataset builder 均 PASS；dataset receipt 的 Ed25519 signature 亦獨立驗證成功。
- C 在 report builder `validate_metrics` fail-closed，尚未建立 final report 或 terminal continuation receipt。
- 根因是把 JSON object key insertion order 當作 schema；canonical JSON 會排序 keys，內容未改變。
- provisional dataset、receipt、signature 與所有本輪 downstream outputs 已完整封存，不提供 release authority。
- authority、terminal、result-release、report-release keys 全部封存並禁止重用。
- v29 必須改用 exact key-set 驗證並增加 canonical round-trip regression。
