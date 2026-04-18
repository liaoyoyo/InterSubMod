"""Theme constants for Weekly Report v2.

Design principles:
- Assertion-Evidence (Michael Alley, Penn State)
- Okabe-Ito colorblind-safe palette (8 colors)
- Tufte data-ink ratio (minimal chartjunk)
- 16:9 canvas with 15/70/15 vertical split (headline/evidence/footer)
"""

# Okabe-Ito colorblind-safe palette
# Reference: https://jfly.uni-koeln.de/color/
OKABE_ITO = {
    "black":      "#000000",
    "orange":     "#E69F00",
    "sky_blue":   "#56B4E9",
    "blue_green": "#009E73",
    "yellow":     "#F0E442",
    "blue":       "#0072B2",
    "vermillion": "#D55E00",
    "pink":       "#CC79A7",
}

# Semantic role assignment
ROLES = {
    "primary":    OKABE_ITO["blue"],         # main data series / LOH
    "secondary":  OKABE_ITO["vermillion"],   # comparison / Non-LOH / FP
    "accent":     OKABE_ITO["orange"],       # highlights / callouts
    "success":    OKABE_ITO["blue_green"],   # validated positive
    "neutral":    OKABE_ITO["sky_blue"],     # intermediate / context
    "warn":       OKABE_ITO["pink"],         # caveats / limitations
    "text":       "#1A1A1A",
    "text_light": "#666666",
    "text_muted": "#999999",
    "bg":         "#FFFFFF",
    "bg_panel":   "#F5F5F5",
    "rule":       "#D0D0D0",
    "rule_light": "#E8E8E8",
}

# Typography hierarchy (points)
FONTS = {
    "title_zh":   28,   # claim headline (Chinese primary)
    "title_en":   17,   # English subtitle (~60% of zh)
    "section":    22,   # section divider
    "body":       18,
    "body_en":    11,
    "label":      14,
    "stat":       20,   # large stat numbers
    "caption":    11,
    "footer":     10,
}

# 16:9 canvas (inches)
CANVAS = {
    "width_in":  13.333,
    "height_in": 7.5,
}

# Layout regions (inches from top-left)
LAYOUT = {
    "margin_x":        0.60,
    "headline_top":    0.40,
    "headline_height": 1.10,      # ~15% vertical
    "evidence_top":    1.70,
    "evidence_height": 5.05,      # ~70% vertical
    "footer_top":      6.90,
    "footer_height":   0.45,      # ~6% vertical
    "en_offset":       0.40,      # EN subtitle y-offset below ZH title
}

# Section-level colors for the 5 findings (consistent across headlines and figures)
FINDING_COLORS = {
    "F1_af_enrichment":   OKABE_ITO["blue"],
    "F2_delta_ngroups":   OKABE_ITO["vermillion"],
    "F3_allele_delta":    OKABE_ITO["blue_green"],
    "F4_segment_rho":     OKABE_ITO["orange"],
    "F5_4group_strat":    OKABE_ITO["pink"],
}


def apply_matplotlib_rcparams():
    """Apply Okabe-Ito theme to matplotlib globally."""
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family":        "sans-serif",
        "font.sans-serif":    ["DejaVu Sans", "Arial"],
        "axes.edgecolor":     ROLES["text"],
        "axes.labelcolor":    ROLES["text"],
        "axes.titlecolor":    ROLES["text"],
        "xtick.color":        ROLES["text"],
        "ytick.color":        ROLES["text"],
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.grid":          True,
        "grid.color":         ROLES["rule_light"],
        "grid.linewidth":     0.6,
        "axes.axisbelow":     True,
        "figure.facecolor":   ROLES["bg"],
        "axes.facecolor":     ROLES["bg"],
        "axes.prop_cycle":    mpl.cycler(color=[
            OKABE_ITO["blue"],
            OKABE_ITO["vermillion"],
            OKABE_ITO["blue_green"],
            OKABE_ITO["orange"],
            OKABE_ITO["sky_blue"],
            OKABE_ITO["pink"],
            OKABE_ITO["yellow"],
            OKABE_ITO["black"],
        ]),
        "savefig.dpi":        300,
        "savefig.bbox":       "tight",
        "figure.dpi":         100,
    })
