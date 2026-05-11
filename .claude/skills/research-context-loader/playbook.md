# Research Context Loader Playbook

## When to invoke (per SKILL.md USE WHEN)

Triggers:
- 「延續研究」「上次的 LOH 研究」「pivot 到 X」「結論驗證」
- 「過去 X 試過嗎」「我們之前是不是測過」
- inject-hypothesis, pivot-direction, research-loop, conclude-research starting up
- User mentions a Thread (A/B/C/D), AUC threshold, NG2 phasing, etc.

## 3-Tier Decision Tree

### Tier 1 — already loaded by CLAUDE.md
- `docs/CURRENT_FOCUS.md` (live主軸, always loaded)
- No action needed

### Tier 2 — light load (this skill, ~3-5k tokens)

Load these 3 files when user starts ANY research task:
1. `InterSubMod/docs/experiments/INDEX.md` (research history index)
2. `InterSubMod/docs/reports/research_landscape/00_INDEX.md` (landscape navigator)
3. `InterSubMod/docs/concepts/2026/04/20260409_研究構想總索引_01.md` (concept index)

After loading these, evaluate: does the user's question need deeper landscape files?

### Tier 3 — deep load (this skill, on-demand per topic)

Match user's question to landscape files (use the table in SKILL.md):

- TO FP / FP analysis → 01
- Self-phasing / V5 fallback / phasing cycles → 02
- ISM feature value / which feature reliable → 03
- Pause/resume judgment / decision matrix → 04
- Evidence chain / inference logic → 05
- Conclusion stability / tier downgrades → 06
- LOH × CN × AF integrated → 07
- Zone-Aware framework history → 08
- Part B / HPFineNGroups upgrade questions → 09
- Research chain registry / which cycle ID → 10
- Truth set / F1 protocol → docs/data_specs/20260422_truth_sets...

Load only the files matching the question. Don't blanket-load all 10.

## Anti-patterns

- DO NOT load all 10 landscape files preemptively — wastes 23k tokens
- DO NOT skip Tier 2 INDEX files — they're cheap and tell you what Tier 3 to load
- DO NOT load this skill for pure code edits (USE SKIP WHEN)
- DO NOT load this skill twice in same session (Claude Code already cached)
