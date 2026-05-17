<!--
建立時間: 2026-05-17
目標: SRE-style Blameless Postmortem 範本 — 對應 /scientific-rigor §9.2 + §11.6 雙環學習
處理範圍: 所有 NEGATIVE / NO-GO 結論的事後分析
關聯檔案:
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md §9.2 (inline template source) + §11.6 (Argyris double-loop)
  - InterSubMod/.claude/skills/known-pitfalls/SKILL.md (P-XX 對應)
  - InterSubMod/.claude/skills/conclude-research/SKILL.md (P6 收尾呼叫)
-->

# Postmortem Template — Blameless SRE-style

> **用法**: 複製本檔到 `InterSubMod/docs/postmortems/<YYYYMMDD>_<topic>.md` 並填入。
> **觸發**: 所有 NEGATIVE / NO-GO / 重大失誤結論強制走此 template。
> **依據**: [Google SRE Postmortem Culture](https://sre.google/sre-book/postmortem-culture/) + [Argyris 1977 Double-loop Learning](https://hbr.org/1977/09/double-loop-learning-in-organizations)

---

# Postmortem: <topic> NO-GO @ cycle_<id>

**Date**: YYYY-MM-DD
**Cycle ID**: cycle_<YYYYMMDD-HHMM-slug>
**Author**: <agent / user>
**Severity**: critical / high / medium / low
**Status**: draft / under-review / resolved

---

## Summary（≤3 行）

（為何此方向無法繼續，核心阻塞 + 學到什麼）

---

## Timeline

- **YYYY-MM-DD**: 啟動方向 + 初始假設 — `cycle_<id>` 註冊
- **YYYY-MM-DD**: Pre-registration commit hash: `<hash>`
- **YYYY-MM-DD**: 中間檢查點 — 觀察到 X
- **YYYY-MM-DD**: 達 NO-GO threshold — 證據 `<artifact_path>`
- **YYYY-MM-DD**: Postmortem 撰寫

---

## Root Cause（**machinery, not blame**）

### Single-loop（執行層）
（哪步驟做錯 / 哪 metric 偏差 / 哪 confound 未控）

### Double-loop（質疑根本假設 — `/scientific-rigor §11.6`）

必須回答 3 問:
1. **我的研究假說背後假設是否有誤？**（非僅「執行哪步錯」）
   - 答:
2. **若假設正確，是否方法論本身需要 pivot？**
   - 答:
3. **若方法正確，是否問題定義需要 reframe？**
   - 答:

---

## What Went Well

- 哪些步驟做得對，可以保留到下次方向

---

## What Went Poorly

- 哪些步驟踩了 pitfall（**必引用 `/known-pitfalls` P-XX**）
  - P-XX: <描述>
- 哪些假設事後證實有誤

---

## Action Items

| # | Action | Owner | Due | Status |
|---|--------|------|-----|-------|
| 1 | (具體可驗證行動) | @user / @agent | YYYY-MM-DD | open / in-progress / done |
| 2 | ... | ... | ... | ... |

---

## Reopen Threshold（Productive Failure 機制 — Kapur 2008-2016）

**若未來再啟動此方向，必須先滿足**:
- [ ] <新證據條件 1>（如：新 dataset / 新方法 / 新 metric 口徑）
- [ ] <新證據條件 2>
- [ ] 重新跑 `/scientific-rigor §7 Pre-registration` 3 欄

**若不滿足上述條件即重啟** → 視為**重複犯錯**，記入 `/known-pitfalls` 新編號 P-XX。

---

## Links

- **Evidence ledger entry**: `state/cycles/<cycle_id>/evaluation.json`
- **對應 memory file**: `memory/project_<slug>.md`
- **/scientific-rigor §8.3 REFLECTION.md** 已建立: yes / no
- **相關 /known-pitfalls 條目**: P-XX (描述)
- **相關 commit / PR**:
- **相關 cycle**: `state/cycles/<cycle_id>/`

---

## Lessons for Future Cycles

- 對 hypothesis_queue 影響: <priority 調整 / 新假設 inject>
- 對 evidence_ledger schema 影響: <若有 schema 演進需求>
- 對 skill 系統影響: <若需新增 /known-pitfalls 條目 / 改 SKILL.md>

---

> **Blameless framing 提醒**: focus on system / mechanism / methodology, **not individual or agent blame**. 個人/agent 錯誤是系統設計的徵兆，不是根因。
