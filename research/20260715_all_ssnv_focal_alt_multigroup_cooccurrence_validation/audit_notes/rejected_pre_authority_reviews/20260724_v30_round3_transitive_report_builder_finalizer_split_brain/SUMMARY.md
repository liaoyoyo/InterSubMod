<!--
建立時間: 2026-07-24 13:54 +08:00
目標: 記錄 v30 第三輪 pre-authority 審查、停止決策與未簽署狀態
處理範圍: C runner current-key rebase、v30 finalizer、report builder transitive dependency、signer shutdown
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/implementation-notes.md
-->

# v30 round 3：遞迴 report-builder finalizer split-brain

## 結論

本候選在正式 authority 建立前停止。頂層 C runner 的 3 個舊路徑與 2 個舊 digest 已正確替換，
direct v30 finalizer 也已綁定 result/report-v9；但 `build_all_ssnv_report_artifact_schema_recovery_v29.py`
仍遞迴執行 v29 finalizer，該 finalizer 綁定已不存在的 result/report-v6 keys。Mendel 因此判定
`REQUEST_CHANGES`（1 HIGH、1 MEDIUM）。依 bounded-attempt 停止政策，不進第四輪修補。

## 輸入與命令

- 輸入 source set：11 個 mode `0444` sources，digest
  `83bf28b8c9455df8bc74390841ff235de39e3e78b486537a6decc6877ddf8c70`。
- Validator：`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/validate_tumor_ref_schema_recovery_authority_v30.py`，SHA-256
  `74c3e445fec567cf81184611b983d260b1b60344b8def63fd2feac099004daf1`。
- Exact probe：`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v30.py`。
- Sidecar suite：`/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B -m pytest -q -p no:cacheprovider /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/schema_recovery_tests/test_validate_tumor_ref_schema_recovery_authority_v30.py`。

## 驗證結果

- Exact probe：exit `0`、`pass=true`、`no_output_writes=true`、canonical `733 passed`。
- Sidecar：`775 passed`。
- Rewritten C prefix：實際 shell exit `0`、stderr 空白。
- Mendel：`REQUEST_CHANGES`；transitive report-builder binding=`false`。
- Authority bundle、三份 formal reviews、result/report receipts 與 signatures：全部不存在。
- 兩個 signer 已以 `KeyboardInterrupt` 安全停止；沒有 signer process 殘留。
- 兩把未使用 private keys：mode `0400`、link count `1`、size `290`；未簽署、未消耗。

## 科學狀態

本阻塞只影響正式簽章發布，不改變 7 datasets、469,849 個同次 LongPhase-S recalibrated
`FILTER=PASS` sSNV、latest sidecar HP/PS read tags、共現/four-state/control/CN 結果或 claim ceiling。
資料分析可稱已完成並經重算；不可稱正式 signed Task-B release 已完成。
