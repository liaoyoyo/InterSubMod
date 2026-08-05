<!--
建立時間: 2026-07-24
目標: 封存 tumor-REF schema recovery v30 round-2 C/finalizer key split-brain 拒絕與修正依據
處理範圍: Mendel approval、External Claude approval、Nash REQUEST_CHANGES、未發布 authority/review/runtime outputs
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/continue_m2v5_after_tumor_ref_promotion_recovery_v30.py
-->

# v30 round-2 pre-authority rejection

## 裁決

Nash 的可重現 high/medium findings 優先於 Mendel 與 External Claude 的 approval。round-2 狀態為 `REJECTED_PRE_AUTHORITY`，沒有建立 formal review、authority、V/R/C output 或 signature。

## High findings

1. C 的 `compose_downstream_script()` 僅重寫 Python launcher 與三個 script path，沒有重寫 key paths/hashes。組成後仍有 3 個不存在的 v28 live paths，而 fresh v30 result/report key paths 出現 0 次；正式 C 會在 runner prefix 失敗。
2. C 指向 frozen v29 finalizer；該 finalizer 仍硬綁已封存的 v29 result/report-v6 keys。active signers 則等待 fresh v30 result/report-v9 targets，形成 split-brain。

## Medium finding

733/775 tests 與 full probe 只驗證 R archive rebase 及一般 composition，沒有驗證完整 C/finalizer key integration，因此可在 formal C 必敗時仍出現假陽性 pass。

## 必要修正

1. C composition 必須 exact rebase 所有 downstream result/report key paths 與 SHA 到 fresh v30 v9 keys，且拒絕任何 stale v28/v29 live binding。
2. 建立 source-bound v30 finalizer，使用 active result/report-v9 public/private paths、SHA 與 exact signer targets；不得修改 frozen historical v29 finalizer。
3. validator/source inventory/V/R/C 全面 pin 新 source identities。
4. canonical 與 sidecar tests 加入完整 composed-script key inventory、raw stale path rejection、fresh key acceptance，以及 finalizer/C signer-target一致性測試。
5. full probe 必須實際覆蓋 C/finalizer integration，再送 round-3 三方重審。

## Superseded approvals

- Mendel round-2: `APPROVE_SUPERSEDED`。
- External Claude Opus round-2: `APPROVE_SUPERSEDED`；envelope SHA-256 `3e30405b37a5fac018107c04b92568b70c210e881d4a6440fadfcf49c956e70b`，CLI stderr 僅有 sandbox dependency warning。

## 綁定

- round-2 validator SHA-256: `424541b75e6e8e47c2c7dff0783b5c19e877e4a32d71b2d22238f979ef492fbc`
- round-2 source-set SHA-256: `d6fe0a64c2346fa2fb3300f5e9c58874022fa4918f3becc7c036e89cf7bbad3f`
- round-2 prior-failure digest: `f47258fc4bad5d15b5f6a1fd0032b8e39c7c05e359aafaa59a8a35d24936cdda`
