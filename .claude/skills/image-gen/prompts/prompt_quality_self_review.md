# Prompt Quality Self-Review (run before invoking codex)

For each prompt YAML about to be sent to AI track, the agent should self-check:

1. **Subject is declarative single sentence** — not imperative, not vague. "LOH region triggers ..." not "draw a chromosome".
2. **Constraints include all 4 Anthropic taboos** — `no gradient overuse`, `no glass morphism`, `no multi-indigo stacking`, `colorblind-friendly`.
3. **Output size matches use case** — 1024×1024 default; only larger for hero figures.
4. **No emoji or decorative unicode** in `subject` or `labels`.
5. **Palette ≤ 3 colors** explicit hex codes.

If any check fails, fix the YAML before invoking codex. Cost saving: prompts/ is in git, fixing now is 0 token; fixing after generation is ~30k token waste.
