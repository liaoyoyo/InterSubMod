<!--
建立時間：2026-08-13
目標：記錄 LongLineage private-first gate 被繞過後的 containment 與後續發布限制
處理範圍：GitHub repository visibility；不判定法律授權完成，也不聲稱可追回既有下載
關聯檔案：InterSubMod/research/20260813_complete_research_handoff/receipts/longlineage_visibility_containment_20260813.json
-->

# LongLineage 公開可見性 incident 與 containment

> **Verdict：`FAIL_CLOSED`。LongLineage 已恢復 `PRIVATE`，但曾被觀察為 `PUBLIC`；來源、授權與 public-safety gate 未通過前不得再次公開。**

## Claim–Evidence–Verdict

| Claim | Evidence | Verdict |
|---|---|---|
| private-first 邊界曾被違反 | GitHub CLI 與 REST 同時回傳 `visibility=PUBLIC`、`private=false` | `CONFIRMED` |
| 已完成持續暴露的 containment | 執行 visibility change 後回讀為 `visibility=PRIVATE`、`isPrivate=true` | `CONFIRMED` |
| 可推論沒有任何第三方取得內容 | GitHub 顯示 0 forks、無 Pages，但有 1 stargazer；平台指標不能證明無 clone/download | `UNVERIFIED` |
| 已可重新公開 research preview | source-origin、license、notices、SBOM、history-safety、P8 均未完成 | `BLOCKED` |

## 執行與驗證

- 輸入：GitHub repository `liaoyoyo/LongLineage`。
- 執行命令：`gh repo edit liaoyoyo/LongLineage --visibility private --accept-visibility-change-consequences`。
- 輸出：GitHub repository visibility（外部狀態，無本機產物）。
- 實際回讀：`visibility=PRIVATE`、`isPrivate=true`、`updatedAt=2026-08-13T02:37:56Z`。

結構化 receipt 位於
`InterSubMod/research/20260813_complete_research_handoff/receipts/longlineage_visibility_containment_20260813.json`。

## 後續 Hard Gate

重新公開前必須同時關閉：21 source mappings 的 origin/license disposition、
`THIRD_PARTY_NOTICES`、SBOM、完整 Git history public-safety findings，以及 P8 verification。
即使這些 gate 通過，定位仍只能是 `research preview / non-production`，production `run`
仍應安全回傳 `KernelBlocked` exit 6。
