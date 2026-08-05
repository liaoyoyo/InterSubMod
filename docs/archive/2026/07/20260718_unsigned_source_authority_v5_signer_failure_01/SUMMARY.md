<!--
建立時間: 2026-07-18
目標: 封存 source-authority v5 第一次未簽署 ceremony 與失去 passphrase 的 key pair
處理範圍: authority/approval JSON、外部 Ed25519 key directory、signer failure provenance
關聯檔案: InterSubMod/docs/provenance/source_authorities/ARCHIVED.md
-->

# Source authority v5 未簽署 ceremony 封存

## 狀態

`UNSIGNED / NOT_AUTHORIZED / NEVER_USED_FOR_FORMAL_EXECUTION`

2026-07-18 第一次 v5 authority 組裝成功，但 one-time source signer 的
OpenSSL 呼叫以 `pkeyutl: Option unknown option -rawin` 結束。沒有產生
`.ed25519.sig` 或 `.pending`；signer process 隨後退出，隨機 passphrase
未落盤且已遺失，因此 encrypted private key 不可再用。

不得把本目錄中的 JSON 或 key pair 當成正式 source authority。正式 v5
只能由新的 key pair、新 source-set digest、新 canonical test evidence、
兩位新的外部 reviewer approval 與新的 detached signature 建立。

## 封存物

| 物件 | SHA-256 | Mode | 結果 |
|---|---|---:|---|
| authority JSON | `99d743da27779721c75bfe015b1e0d094f69345518279fe564928e9cad1df1b9` | `0444` | 未簽署 |
| approval JSON | `60e15e5ce92c36d025de24199e5c71b768abcd3b9d6b421b2c480cf202095aa9` | `0444` | 未簽署 |
| public key | `cd14abe493c146c226ffeea81df571a79ea374e996e59e5d26b06c0fa908b920` | `0444` | 已失去可用 private signer |
| encrypted private key | `bac690113ca9931058fbb037e9992bcccbb715f64c26eb56e03ee841a8dd3265` | `0400` | passphrase 未落盤且已遺失 |

舊 key directory 封存於：

`/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/archive/20260718_all_ssnv_v5_summary_hotfix_failed_signer_01/`

## 驗證

後續以正式 `run_one_time_ed25519_signer_v2.sh`、乾淨 OpenSSL 環境完成獨立
smoke ceremony：簽章驗證成功，private key 退役為 mode `000`。smoke
target 與 detached signature保留在本研究的 `audit_notes/`，不屬於正式
source authority。
