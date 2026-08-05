<!--
建立時間: 2026-07-17
目標: 保存 source-authority v2 簽章有效但 live validator 不相容的完整證據
處理範圍: Task Type B source authority v2 manifest、approval、signature
關聯檔案: InterSubMod/docs/provenance/source_authorities/ARCHIVED.md
-->

# All-sSNV source authority v2 validator incompatibility

## 狀態

`SIGNATURE_VALID / LIVE_VALIDATOR_REJECTED / NOT_AUTHORIZED_FOR_EXECUTION`

## 已驗證事實

- v2 signer 固定使用 `/usr/bin/openssl`，簽章成功並把 source private key 退休為 mode `000`。
- detached signature SHA-256：`e7926154781dd0ba77f29a61006324226a742196b7e7cb81d3dd3004e477563c`。
- 以 regular-file path 執行 `/usr/bin/openssl pkeyutl -verify -rawin` 時回報 `Signature Verified Successfully`。
- live Python validator 把 payload bytes 傳入 pipe，卻令 OpenSSL 以 `-in /dev/stdin` 重新開啟；OpenSSL 回報 `Could not allocate 0 bytes for oneshot sign/verify buffer`。

## 判定

v2 signature 的密碼學內容有效，但正式 validator fail closed，因此 v2 **不得**啟動任何 producer。修正版必須：

1. 將已開啟的 regular-file payload fd 直接傳給 OpenSSL；
2. 對 source、result、report 三層 verifier 使用同一實作；
3. 加入真實 Ed25519 key/signature integration test；
4. 重新建立 source v3 key、全套測試、兩位 reviewer 核准與簽章。

