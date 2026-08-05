<!--
建立時間: 2026-07-24
目標: 封存 tumor-REF schema recovery v27 的 pre-authority 審查拒絕證據
處理範圍: 十個 frozen sources、三方 review transport、未使用 authority/terminal keys、零正式輸出證明
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260724_v27_pre_authority_runtime_contract_and_review_transport_rejection/rejection_evidence.json
-->
# v27 Pre-Authority Rejection

## 結論

v27 在任何 authority、V/R/C 或下游科學輸出建立前被拒絕。Mendel 與 Nash 均提出可直接重現的阻斷問題；外部 Claude Opus 雖回覆 `APPROVE`，但其唯讀檢核未實際觸發 continuation runtime binding，因此不能推翻 deterministic failure。

## 阻斷問題

1. `recovery_result_report_finalizer` 與 `recovery_report_builder` 的 runtime contract 保留 v26 SHA，會在下游執行前 fail closed。
2. validator 要求 `multi_agent_v1`，而正式 agent transport payload 使用 `multi_agent`，無法發布 review。
3. authority 宣稱 26 個 staging patterns 皆有 occupied-state regression，但參數化測試只涵蓋 index 0-23。

## 封存狀態

- 十個 v27 frozen sources 保持 mode `0444`，位於 `sources/`。
- 三份審查結果與 Claude CLI envelope 位於 `reviews/`、`review_transport/`。
- authority-v27 與 terminal-v17 private keys 保持 mode `0400`，已移至 `.config/.../archive/` 並標記永不重用；兩者皆未簽署。
- exact probe：`656 passed`、390 個輸出槽位未寫入、26 個 staging patterns 不存在。
- 科學資料未變：7 datasets、469,849 same-run LongPhase-S recalibrated `FILTER=PASS` sSNVs；confirmed cellular subclone 與 linear ancestry 仍均為 0。
