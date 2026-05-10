# flow_diagram Check Checklist

Score each dimension 1 = pass, 0 = fail. Need ≥4 to pass.

1. **Prompt fidelity** — All `nodes` rendered with their `label`s; all `edges` connect from→to correctly with arrows.
2. **Readability** — Node labels not clipped, font size readable at intended display size (~14pt baseline).
3. **Focal clarity** — Edges visually connect parent → child without crossing important nodes; arrowheads visible.
4. **Design taboos avoided** — Same as concept_diagram (no gradient/glass-morphism/multi-indigo).
5. **Print-friendly** — Black-on-white printable; contrast OK.
6. **Cross-figure consistency** — Same node style / same arrow style across figures.

Note: cairo-rendered diagrams should ALWAYS pass dimension 2 (readability) since text is rasterized exact. If 2 fails, cairo render bug.
