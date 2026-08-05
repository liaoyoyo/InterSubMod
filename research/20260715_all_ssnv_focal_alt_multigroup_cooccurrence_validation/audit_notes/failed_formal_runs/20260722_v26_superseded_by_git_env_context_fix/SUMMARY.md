<!--
建立時間: 2026-07-22
目標: 封存 v26 canonical 與未簽第二輪外審，保留 Git subprocess env context 修正證據
處理範圍: Task-B v7 source-authority pre-authorization chain
關聯檔案: canonical/, reviews/, source_snapshot/
-->

# v26 因 Git env context 修正而封存

## 判定

v26 canonical 為 743 passed、0 failures/errors/skips，A/B 外審皆 APPROVE，且先前 v7 consumer rotation 的 MEDIUM 缺口已關閉。然而 Reviewer A 指出 `release_source_authority.py` 的 live `git rev-parse HEAD` 尚未使用 `GIT_CONFIG_NOSYSTEM=1` 與 `GIT_CONFIG_GLOBAL=/dev/null`。

該項為 LOW defense-in-depth，不形成 false-authority path；但主 agent 驗證發現原修正誤套到 OpenSSL subprocess，且 regression 只檢查字串存在、未鎖定所在 code block。故本輪 review attestation 未簽、v7 authority 未鑄造。

## 後續門檻

- OpenSSL signature verification 恢復使用精確 `CLEAN_SUBPROCESS_ENV`。
- Git HEAD subprocess 使用 clean env 加兩個 Git config 隔離變數。
- regression 必須分別擷取 signature block 與 `current_head` block 驗證，不接受全檔字串命中。
- 使用新 test key 重跑 canonical 與 A/B 外審。
