<!--
建立時間: 2026-07-24
目標: 封存 v29 已簽 authority 與 V 通過後的 historical runner-key live-path 失敗
處理範圍: v29 authority/reviews/V/R failure/key quarantine；科學資料未重算、未改變
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260724_v29_signed_authority_v_pass_r_archived_key_live_path_mismatch/failure_evidence.json
-->

# v29 Signed Authority R Failure Archive

- v29 authority 與 promotion verification（32/32）PASS；正式 R 在第一個 replay output 前 fail-closed。
- exact runner 第 1-358 行仍要求已封存 v28 result/report keys 的舊 live 路徑；第一個缺檔為 result-v5 public key。
- v29 probe/test 驗證了 replayer source 與 helper contract，但未在真實 archived-key live state 執行 runner prefix。
- 未建立 R receipt/log/witness、未啟動 C、未建立 dataset/report 或任何 v29 terminal/result/report signature。
- 不採用舊私鑰 live projection；authority、terminal-v19、result-v6、report-v6 keys 全數封存且禁止重用。
- v30 必須使用 authority-bound archive-path rebase、真實 archived-state regression 與四組全新 keys。
