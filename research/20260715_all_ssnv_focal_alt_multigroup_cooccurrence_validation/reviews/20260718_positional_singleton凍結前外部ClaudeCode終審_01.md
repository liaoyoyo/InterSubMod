<!--
建立時間: 2026-07-18
目標: 留存 positional-singleton 報告與 supplemental release code 凍結前外部終審
處理範圍: report builder、source-authority trust、supplemental finalizer、50-test suite
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/build_positional_singleton_report.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/finalize_positional_singleton_supplemental_release.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_build_positional_singleton_report.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_finalize_positional_singleton_supplemental_release.py
-->

# Positional singleton 凍結前外部 Claude Code 終審

## 審查資訊

- Reviewer: 外部 Claude Code CLI
- 模式: `Read,Grep,Glob` 唯讀
- 原始 verdict: **APPROVE_WITH_RESIDUAL_RISKS**
- Reviewer 明示: **沒有 blocker；無 Critical/High live correctness defect**
- 審查時測試狀態: `40 passed`

## 前輪 findings 判定

| Finding | Reviewer 判定 |
|---|---|
| H1 claim ladder 的 L1/L2 同名反義 | CLOSED |
| M1 M1/M2 funnel 守恆 | CLOSED |
| M2 status/distribution/breakdowns pin | CLOSED |
| M3 summary 連 signed source authority | CLOSED |
| M4 negative tests | PARTIAL |
| L1 恆真 guard | CLOSED |
| L2 amber 對比 | CLOSED |
| L3 顯示硬編碼 | PARTIAL，主要為設計 pin |
| L4 HTML escaping | CLOSED |

Reviewer 確認：

1. source-authority chain 以固定 public-key SHA 為 trust anchor，重新驗
   detached signature、authority ID、source-set SHA、authorized HEAD 與兩份
   external `APPROVE`，邏輯自洽。
2. 三層守恆、dataset/truth 雙邊際、圖表 census 與正文 counts 已對齊。
3. HTML claim ceiling 無過度敘述，atomic publication 與 before/after identity
   皆 fail-closed。

## Reviewer 新 residual risks

1. `validate_source_authority_chain` 缺真簽章 happy/negative tests。
2. supplemental finalizer 除 happy path 與 no-overwrite 外，主要
   `ReleaseError` 分支未測。
3. 必須明示 v4 簽章是 source/commit authority，不是 supplemental report data
   本體的直接封印。
4. supplemental receipt 建立時尚未簽，應明示 detached signature 由
   out-of-band one-time signer 後續產生。
5. 低風險：`m1_flagged <= m1_evaluable`、data-version 重複字面值、圖表零分母、
   claim ladder 非語意 table、missing matrix path 的 KeyError。

## 審查後修正

1. source-authority 新增：
   - 真 Ed25519 happy path
   - tampered signed approval
   - review 數不足
   - review verdict 非 `APPROVE`
   - source-set drift
   - authority-ID drift
2. supplemental finalizer 新增：
   - formal-release link drift
   - audit-summary identity drift
   - detached signature tamper
   - supplemental public-key SHA drift
   - 真雙 formal signatures E2E
   - atomic no-overwrite
3. supplemental receipt 新增：
   - `signature_contract.state_at_receipt_creation =
     UNSIGNED_PENDING_OUT_OF_BAND`
   - 預期 signature path、algorithm、public-key identity、簽後 mode
   - `sealed_artifact_scope`
   - `not_v4_release_reauthorization` 與「pass 需後續簽章驗證」語意
4. 其他修正：
   - `m1_flagged <= m1_evaluable`
   - `data_version` 引用 `SOURCE_AUTHORITY_ID`
   - 圖表使用受控 `ratio()` 零分母錯誤
   - claim ladder 改為語意 `<table>`
   - missing methylation matrix path 改拋 `ReportBuildError`

修正後：

```text
50 passed in 4.21s
```

## 尚待 artifact 終審

- 完成 formal final dataset/report 雙簽章。
- 產生 supplemental Markdown/standalone HTML 與 receipt。
- 以 one-time signer 產生 detached signature，確認私鑰 mode `000`。
- 建立 detached signature verification receipt。
- 執行 desktop/mobile screenshot、overflow/overlap、no-external-asset QA。
- 外部 Claude Code 對最終 receipt、HTML 文字與關鍵數值做唯讀終審。
