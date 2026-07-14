<!--
建立時間: 2026-07-15 04:22
目標: 保存 Claude Code 對 layered workstation GRCh38 與樣本全貌改版的唯讀終審結論
處理範圍: renderer / driver / persistent validator / one-off audit / metrics；不讀大型 generated HTML
關聯檔案:
  - InterSubMod/research/20260715_layered_workstation_sample_overview/claude_code_postimplementation_review_prompt.md
  - InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260715_sample_overview_after/metrics.json
status: reviewed_and_followups_closed
-->

# Claude Code 終審結果

> **VERDICT: PASS；P0=0、P1=0、COMMIT READY=YES。** Claude Code 使用 Opus / high effort，允許工具限制為 `Read,Grep/Glob`，全程唯讀。

## Contract checks

| # | 檢核 | 結果 | 核心證據 |
|---:|---|---:|---|
| 1 | GRCh38 真實長度 + midpoint 落點 | PASS | renderer 以 `(start + end) // 2`；track ratio 與 x-error 實測近零 |
| 2 | 一 W_tree region 一 mark + 座標 fail-closed | PASS | mark/unique/index 都閉合 W_tree；越界即 build fail |
| 3 | 五面板 grain + denominator | PASS | Topo/C/determinacy→W_primary；HP×H3/n_sSNV→W_tree；panel failures=0 |
| 4 | n_sSNV 2–8 與 cap boundary | PASS | build 限 2–8；UI 明示 8 無法拆 natural / cap-compressed |
| 5 | 不搬舊 universe / ontology | PASS | 舊數字不進 renderer；舊分類只在 closed boundary drawer |
| 6 | CN 限定為 region sidecar | PASS | 沒有宣稱連續 CN segment |
| 7 | JSON 隱藏但可追溯 | PASS | HCC1395 4 links，初始 visible=0 |
| 8 | 5 modes / legend / keyboard / mobile / overflow | PASS | 4 viewports；mode unmapped=0；body overflow=0 |
| 9 | freshness gate | PASS | summary SHA + region SHA + renderer SHA + UI contract + mtime |
| 10 | P0/P1 科學、數值、a11y、互動 | PASS | 32/32 page-runs；0 runtime/network errors |

## Reviewer P2 與處理

1. **P2：常駐 smoke 原只深入切換 evidence / determinacy。** 已處理：`validate_workstation_ui.py` 現在逐一切換 `determinacy / evidence / primary-hp / n-ssnv / cn-sidecar`，每模式驗證 active、pressed、mark count、unmapped 與 legend；quick regression 106/106、full regression 538/538 通過。
2. **P2：一次性 capture 只檢查 renderer SHA 長度。** 已處理：audit runner 現在直接計算當前 renderer SHA-256，28/28 sample-page runs 逐頁比對相等。

## 保留限制

- current stored payload 無法 exhaustive 恢復 `single/linear/branched/star`；頁面不報假 count，只保留 retired glossary 與未來 producer gate。
- CN 是 region sidecar、MAX_SNV=8 無法拆 natural/cap-compressed。
- candidate tree 是 regional mutation-state candidate，不是 confirmed clone / lineage；read ALT fraction 不是 CCF。

## 實際命令

```text
Input:  InterSubMod/research/20260715_layered_workstation_sample_overview/claude_code_postimplementation_review_prompt.md
Command: claude -p --model opus --effort high --allowedTools 'Read,Grep,Glob' < prompt.md
Output: 本檔
Actual: VERDICT=PASS; P0=0; P1=0; COMMIT READY=YES
```
