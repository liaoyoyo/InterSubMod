<!--
建立時間: 2026-05-18
目標: 標明 .claude/agents/archive/ 目錄用途 — 保留已被 skill 系統取代的歷史 agent
處理範圍: Q9 archive 過時 agent 清理
關聯檔案:
  - InterSubMod/.claude/skills/weekly-report/SKILL.md (取代者)
  - InterSubMod/AGENTS.md §13 (不可直接刪除原則)
-->

# `.claude/agents/archive/` — 已棄用 Agent 歸檔

> 此目錄保留**已被 skill 系統取代的歷史 agent**。
> 依 `InterSubMod/AGENTS.md §13` AI Agent 預設操作政策「不可直接刪除檔案」，保留作歷史紀錄而非清除。

---

## 歸檔清單

| Agent 檔案 | 棄用日期 | 替代 | 原因 |
|----------|--------|----|----|
| `intersubmod-weekly-research-agent.md` | ~2026-04（推估，archive 目錄建於 2026-04-17）| `/weekly-report` skill | skill 系統提供更靈活的 W1-W7 母稿流程，agent 形式不再需要 |

---

## 為何不直接刪除

依 `InterSubMod/AGENTS.md §13`:
- 不可直接刪除檔案（含 `rm`, `find -delete`, 覆寫式清空）
- 需先搬移到 Archive 區 — **本目錄即 archive**

依 `InterSubMod/.claude/CLAUDE.md §1` Hard Gate：
- 刪除檔案是 Hard Gate 必停場景

依 `/scientific-rigor §8.4 知識追溯`:
- 每個結論必能回答「來自哪份數據 / 哪個 cycle / 哪個 commit」
- 保留歷史 agent 提供「為何 skill 系統取代它」的可追溯依據

---

## 重新啟用步驟（若未來需要）

1. 移回 `.claude/agents/` 主目錄
2. 確認 frontmatter 與當前 Claude Code 慣例一致
3. 評估是否與 `/weekly-report` skill 重複；若衝突需先取捨

---

## 參考替代

| 棄用 agent | 替代 skill | 場景對應 |
|----------|----------|--------|
| intersubmod-weekly-research-agent | `/weekly-report` (W1-W7) → `/pptx-build` 或 `/html-report-build` | 週報母稿 + 多格式輸出 |
