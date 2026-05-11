<!--
建立時間: 2026-05-09 03:09
目標: 記錄 .agents/skills/skills/ 巢狀重複副本的封存原因與保留位置
處理範圍: .agents/skills/skills/ 重複 Codex skills 副本封存
關聯檔案:
  - ../../../../../.agents/skills/README.md
  - ../../../../../.agents/skills/skills/ARCHIVED.md
-->

# Codex Nested Skills Duplicate Archive

## 封存原因

`.agents/skills/skills/` 是 `.agents/skills/` 根層 20 個 Codex skills 的巢狀重複副本。Codex skill discovery 會把巢狀 `SKILL.md` 一併列出，造成同名 skill 重複出現。

使用者已確認此副本為重複、可移除。依 repo 規範不直接刪除檔案，因此改採封存搬移。

## 封存位置

原路徑：

```text
.agents/skills/skills/
```

封存後路徑：

```text
docs/archive/2026/05/20260509_codex_nested_skills_duplicate/nested_skills_duplicate/
```

原位置只保留：

```text
.agents/skills/skills/ARCHIVED.md
```

此 redirect notice 不包含 `SKILL.md`，不會再被 Codex 當作 skill 掃描。

## 差異處理

封存前已比較根層與巢狀副本。巢狀副本中 3 處較正確的 Codex wording 已併回根層：

- `.agents/skills/provenance-tier-audit/SKILL.md`
- `.agents/skills/research-loop/references/SCALE_LADDER.md`
- `.agents/skills/validation-protocol/SKILL.md`

`.agents/skills/README.md` 以根層版本為準，並補入 Claude Code harness 到 Codex 的轉譯規則。
