<!--
建立時間: 2026-07-24 07:48 Asia/Taipei
目標: 封存 v29 round1 在 authority 建立前發現的 failed-v28 terminal key projection 缺口
處理範圍: 10 份 frozen source、6 份外部 Claude transport、Mendel 初始與修正審查、Nash 審查、orchestrator 重現
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/rejected_pre_authority_reviews/20260724_v29_round1_failed_v28_terminal_binding_gap/rejection_evidence.json
-->

# v29 round1 pre-authority rejection

## 結論

v29 round1 source-set `0087b195ce2b7bfb0495f4f3e2c879b11484a645218a318402d14b47879a019e`
在 authority 建立前遭拒絕。validator 的真實 `terminal_key_state` 有 30 fields／15 public keys，
continuation 的 exact projection 只有 28／14，缺少 `failed_v28_public_key` 與
`failed_v28_private_key`；因此 C 必然拋出 `Promotion signature/key identity drift`。

完整 evidence SHA-256：
`6d44a75b32737c0e16a7c05c3801ade734977336982a7f08260b70610442cafa`。

## 審查鏈

- Mendel 初始判定 `APPROVE`，但未做 validator→continuation 真實跨模組比較；該判定已撤銷。
- Mendel 看過可重現反例後改判 `REQUEST_CHANGES`：1 high、2 medium、1 low。
- Nash 獨立判定 `REQUEST_CHANGES`：1 high、2 medium、1 low。
- External Claude Opus transport 回傳 `APPROVE`，但未發布成 formal review；依 strictest reproducible review wins，不具 authority 效力。
- Orchestrator AST／isolated-import 重現得到 30 vs 28、15 vs 14、exact equality=false。

## 仍有效證據

- v28 archive 742-entry inventory、authority signature 與 provisional dataset signature 仍有效。
- report mapping 的 exact key-set、canonical sorted roundtrip、missing/extra-key 修補仍有效。
- round1 canonical read-only probe 曾為 699 passed、424 forbidden slots、28 staging patterns、無 protected writes。
- 科學 payload 未改：7 datasets、469,849 個 same-run LongPhase-S recalibrated `FILTER=PASS` biallelic sSNV。
- claim ceiling 未改：只支持 latent molecular substructure candidates；confirmed cellular subclone=0、linear ancestry=0。

## Key lifecycle

round1 未建立 authority、未發布 formal reviews、未建立 C/result/report signature。authority、terminal、result、report
四把 active v29 private keys 均為 mode `0400`、link-count 1、未消耗，因此只允許在修正後的同版
`v29 round2 pre-authority candidate` 重用；不得跨版本或在任何簽章後重用。

## Round2 必要條件

1. continuation 完整綁定 failed-v28 rotation/state/return-map/private/protected paths。
2. regression 必須直接使用 validator 真實 30-field terminal state，不再使用手工 14-key fixture。
3. validator 與 probe 必須綁定本 archive 的 20-file inventory、evidence 與 summary。
4. 重新 freeze source hashes，重跑完整 tests、canonical read-only probe與三份 fresh reviews。
5. 只有三方皆 `APPROVE` 才可建立 authority，之後依序執行 V → R → C。
