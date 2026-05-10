# Anthropic / thariqs Design Taboos (audit reference)

Per [thariqs/html-effectiveness](https://thariqs.github.io/html-effectiveness/) FAQ:

| # | Taboo | Detection regex / signal |
|---|-------|--------------------------|
| 1 | gradient overuse (>=2 gradients) | `linear-gradient` or `radial-gradient` count >= 2 |
| 2 | glass morphism | `backdrop-filter:` or `backdrop-blur` |
| 3 | multi-indigo stacking (multiple indigo shades) | `#4F46E5` `#6366F1` `#818CF8` (>=2 of) or `--color-indigo` chains |
| 4 | emoji-decorated headers | `<h[1-6][^>]*>[^<]*[\u{1F300}-\u{1FAFF}]` |
| 5 | drop shadows on text | `text-shadow:` (any non-`none`) |
| 6 | glow effects | `filter: drop-shadow` outside icons |

Also flagged: garish neon colors (#FF00FF, #00FFFF), 3D effects (`transform: rotateY`), comic sans / decorative fonts.
