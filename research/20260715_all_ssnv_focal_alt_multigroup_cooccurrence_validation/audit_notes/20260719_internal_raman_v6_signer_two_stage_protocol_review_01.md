<!--
建立時間: 2026-07-19 00:24 CST
目標: 保存 v6 one-time signer 兩階段 authority protocol 的 fresh adversarial internal review
處理範圍: signer ceremony semantic boundary；不包含 production source-authority consumer 的最終驗證
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_one_time_ed25519_signer_v6.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_one_time_ed25519_signer_v6.py
-->

# v6 signer two-stage protocol internal review

## Review identity

| 欄位 | 值 |
|---|---|
| Reviewer | Raman subagent `019f72d5-599c-7323-bd51-45f46e58fda0` |
| Review mode | fresh adversarial, read-only |
| Verdict | **APPROVE** |
| Confidence | high |
| Signer SHA-256 | `e032a5c56dcb40b81b4c963413d7190a97894c7c5258c67f29700bc64503dd49` |
| Test SHA-256 | `431f8412a3e86200e3d30fbabe489a514d429886487471b6bc6d187ccbce4897` |

## Finding summary

沒有 HIGH 或 MEDIUM blocking finding。前版在 terminal snapshot 後、Python `finally`
真正返回前仍可替換 canonical path；v6 不宣稱消除任意 filesystem race，而是明確切分責任：

1. signer ceremony 只發布輸出並退休 one-time private key；
2. 成功事件固定為
   `CEREMONY_OUTPUT_AVAILABLE_REQUIRES_INDEPENDENT_VERIFICATION`；
3. 固定欄位為 `consumer_verification_required=true` 與
   `release_authority_granted=false`；
4. signer source 不含成功常值 `SIGNED`；
5. release authority 必須由 signer process 結束後的獨立 consumer 重新開啟 live
   target、signature、public/private key 與 key directory，完成 crypto、identity、mode、
   protected-source 及 late-guard 驗證後才可授予。

## Closure matrix

| 項目 | 判定 | 理由 |
|---|---|---|
| False `SIGNED` 語意 | CLOSED | source 無 `SIGNED` 成功常值，事件明確為 provisional |
| signer 誤授 release authority | CLOSED | `release_authority_granted=false` 為固定 literal |
| teardown close race | semantic finding CLOSED | race 可破壞 live path，但不能再被 signer 宣稱為 authority |
| close 後 signature replacement | CLOSED | regression 注入成功，signer只回 provisional，獨立 OpenSSL consumer拒絕 |
| target replacement / same-inode pwrite / fsync replacement | CLOSED | 原注入點維持 fail-closed |
| failure truth | CLOSED | failure record依 live private-key FD回報狀態 |
| AST/behavior tests | SUFFICIENT FOR SIGNER SEMANTICS | 禁止 `SIGNED`，並驗證 provisional authority booleans |
| production independent authority | REQUIRED OUTSIDE THIS REVIEW | 必須由 fresh canonical與實際 consumer chain另行證明 |

## Residual gate

LOW scope limitation：目前 regression 的 independent consumer 只示範 OpenSSL crypto rejection，
尚未替代 production `release_source_authority.py` 的完整 path/mode/digest/receipt chain。舊 v20
canonical JUnit 早於本次 signer/assembler bytes，不能作為最終證據。兩者均維持正式 release
hard gate，不因本次 `APPROVE` 自動關閉。

