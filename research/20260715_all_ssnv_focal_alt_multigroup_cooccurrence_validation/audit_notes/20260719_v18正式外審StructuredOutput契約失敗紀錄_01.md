<!--
建立時間: 2026-07-19
目標: 留存 External Reviewer A/B v18 各次 fail-closed 原因、封存位置與最終成功證據
處理範圍: Claude CLI transport、signed test-run 淘汰、StructuredOutput command/consumer parity 與 transcript payload binding
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/attested_release_evidence_v1.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_external_claude_review_v18_attested.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_attested_release_evidence_v1.py
-->

# External Reviewer v18 StructuredOutput 契約失敗紀錄

## 狀態

- Task type：B comprehensive validation。
- 歷史：前三輪 reviewer process 都 **fail-closed**，未把未簽文字當 evidence。
- 目前：第四輪 A/B 均正式 `APPROVE`，各自 process attestation 與簽章已由
  consumer 驗證；v5 source authority 已簽署並通過全角色重驗。
- 當前授權範圍：`all_7_datasets_469849_sites_and_release_chain`。

## 輸入與命令

- 輸入 source-set SHA-256：
  `c8bd2f49932fd5ee165c1b7ff46cf0d0b410245d52880af312c75170e874980a`。
- signed test manifest SHA-256：
  `72e10fe9cb4ca04daf620035dfbbc73b42b6bc10b6ac1e11e3ef8d6ee35caa16`。
- 執行方式：釘選 Python 以 `-I /proc/self/fd/<runner-fd>`、clean environment、
  `--model claude-opus-4-8 --effort max --output-format stream-json --json-schema`
  執行 A/B。

## 第一輪失敗：工具未允許 StructuredOutput

兩個 process 均在 Claude 完成後由
`parse_claude_stream_transcript()` 拒絕：

```text
EvidenceValidationError: Claude transcript process identity drift
```

根因是 producer command 只允許 `Read,Glob,Grep`，但 consumer 正確要求
`Read,Glob,Grep,StructuredOutput`，且要求 `result.structured_output`。實際 init event
因此只有三個工具，Claude 以一般文字回傳 JSON，不能形成可簽的 structured attestation。

## 第二輪失敗：transport schema 過於複雜

producer 加入 `StructuredOutput` allowlist 後，仍把完整授權規則直接編碼成深層
JSON Schema。Claude Code 2.1.202 在這個 schema 下把 init tools 降為
`Read,Glob,Grep`，沒有註冊 `StructuredOutput`；A/B 再次由相同 consumer
條件拒絕，沒有 attestation 或簽章。

實際 probe 顯示：

- 最小 schema：四工具存在，`result.structured_output` 正常。
- 完整深層 schema，以及僅移除 `allOf` 的版本：都降為三工具。
- 17 欄、只限制基本型別的淺層 schema：四工具與 structured output 正常。

這與 Claude Code 官方文件對複雜 schema 可能失敗、應縮小 schema 範圍的說明一致：
[Structured outputs](https://code.claude.com/docs/en/agent-sdk/structured-outputs)、
[CLI usage](https://code.claude.com/docs/en/cli-usage)。

## 第三輪失敗：探索性 Read 被誤當必要 Read

A/B 都已成功以全檔、非空、無 error 的 `Read` 讀完七個指定 source，
並各自產生一次有效 StructuredOutput；但舊 parser 對**所有** Read 都要求
全檔成功：

- A 的非必要 `offset/limit` Read 被誤拒為 bounded required Read。
- B 對尚未建立的 v5 authority 做存在性探測，得到預期的 file-not-found，
  被誤拒為 failed required Read。

完整 transcript 與修正說明：

`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v18_optional_exploratory_read_overstrict/`

修正後仍要求七個指定路徑各有至少一次 canonicalized、只含 `file_path`、
成功且非空的 Read；其他探索 Read 可 bounded/error，但不能替代 required evidence。

## 封存

- A/B request、command、stream transcript、stderr：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v18_structured_output_not_allowed/`
- 該 source-set 的 signed test-run：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v22_signed_test_superseded_by_structured_output_fix/`
- StructuredOutput 修正後、prompt parity 修正前的 signed test-run：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v22_signed_test_superseded_by_prompt_parity_fix/`
- 第二輪 schema 過度複雜的 A/B transport：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v18_structured_schema_too_complex/`
- 深層 schema source-set 對應、已淘汰的 signed test-run：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v22_signed_test_superseded_by_transport_schema_fix/`
- 探索性 Read 語意修正前、已淘汰的 signed test-run：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v22_signed_test_superseded_by_optional_read_semantics_fix/`

所有舊 test private keys 均已由一次性 signer 退役為 mode `000`；沒有重用舊簽章。

## 修復與驗證

1. allowlist 同步為 `Read,Glob,Grep,StructuredOutput`。
2. runner 直接使用 bound evidence module 的 `REVIEW_SYSTEM_PROMPT`，移除 prompt 雙份來源。
3. transport schema 只固定 17 個欄位、基本型別與 `additionalProperties=false`；
   精確值、UUID、verdict、F1/F2、空 blocker、source-set、HEAD、JUnit、
   evidence 完整性仍由 `validate_clean_review_payload()` fail-closed 驗證。
4. transcript consumer 要求 StructuredOutput tool use 恰好一次、其 input
   逐欄等於 `result.structured_output`，且要有成功非空的匹配 tool result；
   所有 tool use/result 亦須一一配對。
5. 既有 testcase 新增 tool/prompt producer-consumer parity，以及
   zero/duplicate/mismatched/error StructuredOutput 負向案例。
6. 真實最小 Claude transport probe：
   init tools 含四工具、產生 `StructuredOutput` tool use、結果含
   `structured_output={"ok":true}`、退出碼 0。
7. 最新聚焦回歸：`25 passed`。
8. 凍結檢查：protected source `23/23`、canonical tests `50/50`、
   supporting source `6/6` 均為 mode `0444`。
9. required/optional Read 分流：required 七檔仍需 full/success/nonempty；
   optional bounded/error Read 不造成 false-negative。
10. 補上 missing result 與 unexpected result 負向 regression。

先前 transport-schema 修正前的正式測試曾得到
`735 passed / 0 failed / 0 error / 0 skipped`，但該簽章與
source-set SHA-256
`2f76320735b29d210d1c460e882e5ea4294ab1c424fb96ebf1b48a3447fd7c62`
已因後續 parser/schema 修正淘汰並封存，**不得**作為目前 release evidence。

## 最終有效證據

- v6 canonical test：`735 passed / 0 failed / 0 error / 0 skipped`。
- source-set SHA-256：
  `9e425725162471dfda8df47adf00cae1b779e83931748bd9ab7eebfc11b6b323`。
- JUnit SHA-256：
  `a0dfa9cb1ef91c43d01637d8f6db0041cbf6a1a90b84c48107a168efdf852a39`。
- signed test manifest SHA-256：
  `aff0679eed70b5eba936ef4e5200f159adba4a0c8d230fb97698d734a95693cd`。
- Reviewer A attestation/signature SHA-256：
  `98f9df97833b0b6ccd92582e035de2a3781d719cb942f449fc1827ade6be1f65` /
  `f1e857c271f0701354320f0a80f55f24cf24b6276ab29e6ac9c1ce972526fa1d`。
- Reviewer B attestation/signature SHA-256：
  `09bb7b7b6a8934abb9c02f8f9e03662ba373d45b655bdbd9e47dba5e02c429d5` /
  `70a2281bf8ebd5265782593c93fa1c503cf1a4e11ca82c426a0cc39ba23e6bd9`。
- source authority manifest/approval/signature SHA-256：
  `0b68438b1dfc8d4123a03089e98f79bdd9d3ee0ab77110eb40937e9426e851fc` /
  `034328817b3e480880f9c51088fe4f00e36d1947fae6f8528cc67be853a8b4fb` /
  `3fe83497b4e5a38ac119b5b51265d51ae53197cc030322a4cc6cad5552dd0a89`。

本紀錄不把任何 fail-closed round 的未簽文字 verdict 當作正式 review evidence。
