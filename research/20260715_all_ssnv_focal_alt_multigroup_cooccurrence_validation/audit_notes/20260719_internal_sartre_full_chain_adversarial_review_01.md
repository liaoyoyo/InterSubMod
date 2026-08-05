<!--
建立時間: 2026-07-19
目標: 保存 v17/v10/v21 trust-chain 的獨立 adversarial review
處理範圍: source authority、external review、canonical JUnit、supplemental release
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/
-->

# v17/v10/v21 trust-chain 獨立 adversarial review

## Review identity

- Reviewer agent: Sartre (`019f73fb-a74b-7920-a857-e884b04354d0`)
- Review method: STRIDE + claim-evidence trace
- Reviewed source-set SHA-256: `6abaa35ee6045d9cb483623d1a8d08628100a162b3d1b64b91d0ba869963af77`
- Reviewed assembler SHA-256: `424caff...7336`
- Reviewed canonical JUnit SHA-256: `4cd974...f53bc`
- Canonical JUnit counts: `722 passed / 0 failed / 0 errors / 0 skipped`
- Verdict: **REQUEST_CHANGES**
- Release state at review: **UNSIGNED / NOT_AUTHORIZED**

## Blocking findings

### F1 HIGH: external review 與 JUnit 缺少獨立固定錨點

v17 normalizer 只能確認 review JSON 的格式與自我宣告欄位。v10 assembler
讀取當下的 assembler、JUnit 與 A/B review，再要求彼此一致；production
validator 也只把 live identity 與 review-controlled identity 比對。父目錄可寫時，
攻擊者可以一致置換 assembler、JUnit 與兩份自陳 `APPROVE` review，形成可簽署但
未經真實 reviewer/test-run 證明的 authority。

只讀 predicate probe:

```text
normalizer_accepts_self_asserted_A_B=True
assembler_accepts_self_consistent_wrong_live_identity=True
production_shape_gate_accepts_wrong_assembler_identity=True
```

修補要求:

- protected validator 必須固定 assembler identity；
- JUnit 必須由一次性簽章的 test-run manifest 綁定 source-set、HEAD、測試來源、
  執行命令與 XML identity；
- A/B review 必須有獨立 process/session evidence 與不可在 review JSON 內自陳的
  attestation binding。

### F2 HIGH: supplemental receipt 存在 ABA snapshot mismatch

supplemental finalizer 先計算 input identities，再由 pathname 重新解析 JSON，
最後只確認 pathname 已恢復原 identity。攻擊者可在解析前換成 `S_bad`，在 terminal
identity check 前恢復 `S_good`，使 receipt 的 counts 來自 `S_bad`、宣告的 hash
卻是 `S_good`。

修補要求:

- 所有 input identity 與 payload 必須來自同一組 retained FD snapshot；
- 需有獨立 post-sign verify-only consumer，重新驗簽並重算所有 transitive inputs。

## Medium findings

### F3: private-key retirement 未納入 descriptor binding

source/result/report consumer 只做一次 pathname mode check，無法證明 terminal
state 仍是 mode `000`。

### F4: singleton live consumers 仍接受舊 v4 authority

singleton audit 與 report builder 固定舊 v4 authority/public key。這與 fresh v5
core key rotation 不一致；不能宣稱所有 live consumers 已旋轉。

### F5: canonical JUnit 缺少 authority-positive attack tests

缺少 self-consistent A/B substitution、wrong assembler/JUnit identity、private-key
race、supplemental ABA、post-sign supplemental verification 及 old-key revocation
scan。

## 已確認的正向事實

- 23/23 protected source roles 與 `6abaa35e...af77` 重算一致，全部 mode `0444`。
- signer、normalizer、assembler 與審查指定 SHA 一致。
- v21 JUnit identity 與 `722/0/0/0` counts 一致。
- v10 沒有 literal fixed-point self-reference。
- 一次執行內的 review payload FD binding 可抵抗一般 TOCTOU。
- 四個 signer process 當時只在等待 target，尚未產生任何正式 signature。

## Disposition

v17/v10/v21 保留為可追溯、但失效且不可授權的迭代。只有在 F1-F5 均有攻擊
回歸測試、內部 reviewer 重審及外部 Claude Code review 後，才可建立下一版
authority 與送出 `SIGN`。
