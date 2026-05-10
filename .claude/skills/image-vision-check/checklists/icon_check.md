# icon Check Checklist

Score each dimension 1 = pass, 0 = fail. Need ≥4 to pass.

1. **Prompt fidelity** — Shape matches `shape` field; glyph matches `glyph` field exactly.
2. **Readability** — Glyph occupies > 40% of inner area, clearly visible against fill.
3. **Focal clarity** — Single dominant element; no unintended texture/border.
4. **Design taboos avoided** — Solid fill (not gradient), no drop shadows, no glow effects.
5. **Print-friendly** — Glyph contrast vs fill ≥ 4.5:1; works in B&W (luminance contrast OK).
6. **Cross-figure consistency** — Same icon set style across multiple icons (size, padding, glyph weight).

Note: cairo-rendered icons should ALWAYS pass 1, 2, 3, 4 (deterministic). If any fails, cairo render bug.
