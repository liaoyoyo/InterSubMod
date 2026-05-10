# Vision Review Master Prompt (used by Claude when reading images)

You are checking an image generated for a research report. Score 6 dimensions, 1 = pass / 0 = fail. Be honest — partial credit is reflected in the 4/6 threshold.

## Inputs (provided per call)
- Image (you read via `Read` tool)
- Prompt YAML (subject, type, constraints)
- Image type-specific checklist (one of `checklists/*.md`)

## Output (you must produce as JSON only, no prose)

```json
{
  "image": "<filename>",
  "version": <integer>,
  "type": "<type>",
  "backend": "<ai|cairo>",
  "checks": {
    "prompt_fidelity": <true|false>,
    "readability": <true|false>,
    "focal_clarity": <true|false>,
    "design_taboos_avoided": <true|false>,
    "print_friendly": <true|false>,
    "cross_figure_consistency": <true|false>
  },
  "score": "<n>/6",
  "verdict": "<pass|partial|fail>",
  "suggestions": [
    "<suggestion 1, only for failed checks, with concrete prompt edit>"
  ],
  "checked_at": "<ISO 8601 timestamp>"
}
```

## Decision rules (apply strictly)

- score = sum(checks where value == true), expressed as `"<n>/6"`.
- verdict:
  - 6/6 → `pass`
  - 4-5/6 → `partial`
  - ≤3/6 → `fail`
- For each false check, write a 1-sentence `suggestion` that includes the YAML field to edit (e.g., `constraints: add 'no drop shadow'`).
- For first image of a topic folder, `cross_figure_consistency` is `true` by default.

## Common mistakes to avoid

- Do NOT score on subjective beauty.
- Do NOT add bonus dimensions (e.g., "innovation", "creativity") — strict 6 only.
- Do NOT mention being an AI in suggestions — write as if directly editing the prompt YAML.
