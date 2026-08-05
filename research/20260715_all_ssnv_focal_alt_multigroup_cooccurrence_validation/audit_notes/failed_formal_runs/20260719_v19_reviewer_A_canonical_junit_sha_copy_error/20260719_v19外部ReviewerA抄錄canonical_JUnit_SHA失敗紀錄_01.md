<!--
建立時間: 2026-07-19 13:20 +08:00
目標: 封存外部 Reviewer A v19 第一次正式嘗試因 canonical JUnit SHA-256 抄錄錯誤而未發布 attestation 的完整證據
處理範圍: Reviewer A command/request/stream/stderr；不包含 Reviewer B、正式 authority 或 downstream producer
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v23_command_parity_attested_canonical.xml
狀態: failed_formal_attempt_archived_not_authorized
-->

# v19 外部 Reviewer A canonical JUnit SHA 抄錄失敗紀錄

## 一句結論

Reviewer A 的 Claude process 回傳 `APPROVE`、`blocking_findings=[]`，但 structured payload 抄錯
canonical JUnit SHA-256；本地 attested validator 因 exact identity 不相等而 fail-closed，因此本次嘗試
**沒有發布 attestation JSON、沒有簽章、沒有授權任何正式執行**。

## 輸入與命令

- runner：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_external_claude_review_v19_attested.py`
- canonical JUnit：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v23_command_parity_attested_canonical.xml`
- 執行型態：以同一個 mode `0444` runner 的 bound FD 啟動 `--reviewer A`；完整 argv 與 request
  保存在本目錄的 `.command.json` 與 `.request.md`。

## 實際結果

| 欄位 | 值 |
|---|---|
| Claude process | success |
| Reviewer verdict | `APPROVE` |
| Blocking findings | `0` |
| 預期 JUnit SHA-256 | `f4842aa9aae54b94122d08c966b8fb199de4f8b3f0c5737459531d92754eb090` |
| Reviewer payload SHA-256 | `f4842aa9aae54b94b94b199de4f8b3f0c5737459531d92754eb090` |
| Local validator | `EvidenceValidationError: External review is not a clean bound approval` |
| Formal attestation JSON | 未建立 |
| Detached signature | 未建立 |
| Release authority | 未授予 |

Reviewer payload 的其他主契約欄位仍為 source-set
`2e5da8df553662c0d77ab6fcabc634079d648a0d70e27fe1b9133b93216df377`、
Git HEAD `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`、738 tests、0
failure/error/skip；但 exact artifact identity 是硬閘門，不能用語意近似取代。

## 封存檔 SHA-256

```text
4ae1ae93275e09e440fe60e13732722126fd8dc4cfa381f939d586ebb69d6e5a  command.json
45e47cf0a0aa2c0ac642949ef39ad03e1747fb96bbb8230d09d95f9d7814b3ca  request.md
806193648350fcd906d99c9cac51b05036195ab8d48653e237b8efa2e826a3e8  stream.jsonl
222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a  stderr.txt
```

## 判定與下一步

這是 reviewer structured output 的 transcription/copy failure，不是 source-level 科學或工程 blocker。
保留本次完整 transcript 供稽核；下一次使用完全相同的 frozen v19 runner、source-set、JUnit 與
Reviewer A key，重新取得 exact-bound review。只有新 attestation 通過本地 validator、完成 detached
signature 並經獨立 consumer 驗證後，才可組裝 v6 source authority。
