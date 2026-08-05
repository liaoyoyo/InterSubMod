<!--
建立時間: 2026-07-19
目標: 保存第一次 formal v22 fail-closed 執行證據
處理範圍: pytest runner / Python FD executable 繼承
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_attested_canonical_pytest_v1.py
-->

# Formal v22 第一次執行失敗紀錄

## 結果

- 狀態：`FAIL_CLOSED_ARCHIVED_UNSIGNED`
- JUnit：735 tests、15 failures、0 errors、0 skipped
- 未產生：canonical JUnit、test-run manifest、Ed25519 signature
- test private key：未使用、未退休

## 根因

pytest process 由 bound Python FD 啟動，因此 child 的
`sys.executable=/proc/self/fd/86`。15 個測試會再以 `sys.executable`
啟動 subprocess，但未傳遞 FD 86，因而得到 `FileNotFoundError`。

## 修正

pytest 仍以 `executable=/proc/self/fd/N` 執行；argv0 改為同一
Python artifact 的 canonical resolved path，並在 signed manifest 的
`runtime.pytest_process` 同時記錄 argv0、FD executable 與 Python identity。

## 封存內容

- `.pytest_full_pre_authority_v22_bound_test_fds/`
- `pytest_full_pre_authority_v22_raw_evidence/`
- `pytest_full_pre_authority_v22_attested_canonical.stdout.txt`
- `pytest_full_pre_authority_v22_attested_canonical.stderr.txt`

以上均為第一次失敗執行的原始內容，未刪除或覆寫。
