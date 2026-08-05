<!--
建立時間: 2026-07-22
目標: 封存 promotion v2 第一次正式審查拒絕時的 frozen source
處理範圍: P/V/R/C 與外部唯讀 probe，不含任何研究輸出
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/
-->

# Promotion v2 正式審查拒絕封存

此目錄保存 2026-07-22 正式審查時的五份 `0444` 原始 bytes。這些 source 不得作為後續 release authority；用途僅為重現拒絕判定。

## 拒絕原因

1. Portable report 動態載入的三個 plugin module 未納入 reviewed source、gate input 與 mutation monitor。
2. Downstream primary/QA Python 由可變 symlink pathname 重新解析，未以 retained FD 執行。
3. 已退休為 mode `0000` 的 private key 被加入 inotify file watch，實際會回傳 `EACCES`。
4. Nested numeric JSON 與 portable QA 使用非 type-strict 比較，可接受 `bool` 或 numeric string。
5. Promotion producer 的 artifact publication 未由 parent wait-normal-exit-zero attestation 約束。

## 審查狀態

- 內部正式 reviewer A：`FORMAL_VERDICT: REJECT`，1 HIGH、2 MEDIUM。
- 內部正式 reviewer B：`FORMAL_VERDICT: REJECT`，2 HIGH。
- 外部 Claude Code 對較窄 scope 回傳 APPROVE，但未涵蓋上述缺口，因此不得用於授權。

後續必須修改 source、重新凍結，並由全新 reviewer 從零審查。
