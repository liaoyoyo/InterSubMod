# Design Taste Check Prompt (used by Claude when auditing rendered HTML)

You are auditing a rendered HTML document for Anthropic / thariqs design taboos.
Read the inline `<style>` block AND the body markup. Score each taboo present (1) or absent (0).

## Output format (JSON ONLY)

```json
{
  "taboos_detected": {
    "gradient_overuse": <true|false>,
    "glass_morphism": <true|false>,
    "multi_indigo_stacking": <true|false>,
    "emoji_decorated_headers": <true|false>,
    "drop_shadow_text": <true|false>,
    "glow_effects": <true|false>
  },
  "n_taboos": <integer>,
  "verdict": <"clean"|"minor_issues"|"redesign">,
  "suggestions": [
    "<concrete edit, e.g., 'Remove gradient on .hero — use --color-bg flat fill'>"
  ]
}
```

## Decision rules

- `n_taboos == 0` → `verdict: clean`
- `n_taboos == 1` → `verdict: minor_issues`
- `n_taboos >= 2` → `verdict: redesign`

For each detected taboo, write a 1-sentence suggestion.
