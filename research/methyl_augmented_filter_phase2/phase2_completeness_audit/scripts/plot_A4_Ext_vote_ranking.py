#!/usr/bin/env python3
"""Plot A4-Ext 6 algorithm vote-score ranking (horizontal bar).

Output:
  research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A4_Ext_algorithm_vote_ranking.png
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# CJK font fallback per MEMORY.md feedback_matplotlib_cjk_font_rule
rcParams["font.sans-serif"] = ["DejaVu Sans", "Droid Sans Fallback", "Noto Sans CJK TC", "sans-serif"]
rcParams["axes.unicode_minus"] = False

REPO = Path(__file__).resolve().parents[4]
OUT = REPO / "research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A4_Ext_algorithm_vote_ranking.png"

# (algo, vote, expected_outcome_short, color)
ROWS = [
    ("Threshold-only\n(caller_af + V6_off_NG)", 3, "1 hr robustness check", "#2a9d8f"),
    ("Per-zone LR\n(LOH-in / LOH-out)", 3, "HCC1937 case study", "#2a9d8f"),
    ("Boolean rule /\nzone gate", 2, "cycle 3 partial cover", "#e9c46a"),
    ("Per-AF-bin LR\n(low/mid/high)", 2, "collider bias risk", "#e9c46a"),
    ("Ensemble\n(LR + RF voting)", 1, "math 0+0=0", "#e76f51"),
    ("Calibration\n(Platt / isotonic)", 1, "F1 invariant", "#e76f51"),
]

# Sort by vote desc, stable
ROWS_sorted = sorted(ROWS, key=lambda r: -r[1])

labels = [r[0] for r in ROWS_sorted]
votes = [r[1] for r in ROWS_sorted]
notes = [r[2] for r in ROWS_sorted]
colors = [r[3] for r in ROWS_sorted]

fig, ax = plt.subplots(figsize=(11, 6.2))
y_pos = range(len(labels))[::-1]  # top = highest vote
bars = ax.barh(list(y_pos), votes, color=colors, edgecolor="black", linewidth=0.7, alpha=0.9)

# Vote score annotation
for bar, v, note in zip(bars, votes, notes):
    ax.text(v + 0.06, bar.get_y() + bar.get_height() / 2,
            f"{v}/5  -  {note}",
            va="center", ha="left", fontsize=10)

ax.set_yticks(list(y_pos))
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlim(0, 5.5)
ax.set_xticks([0, 1, 2, 3, 4, 5])
ax.set_xlabel("Vote score (1=not worth, 3=worth conditional, 5=must do)", fontsize=11)
ax.set_title("A4-Ext: 6 algorithm vote-score ranking (PI Goal 4 due-diligence)\n"
             "anchor: LOSO LR mean dF1 = -0.00004; caller_af direction inconsistent across 5 samples; F1 ceiling-bounded 4/5",
             fontsize=11)

# Vertical guides for vote tiers
for x, label, color in [(1.5, "low", "#e76f51"), (2.5, "marginal", "#e9c46a"), (3.5, "worth", "#2a9d8f")]:
    ax.axvline(x, linestyle=":", color=color, alpha=0.5, linewidth=0.8)

# Footer
ax.text(0.0, -1.0,
        "Top 1 (tied): Threshold-only quick check (1 hr) + Per-zone LR HCC1937 (1.5 hr) = <3 hr cycle-5 supplement.\n"
        "No algorithm in this set is expected to flip LOSO sample-level ~0 finding; root cause is distribution shift not model bias.",
        fontsize=9, color="#555", transform=ax.get_yaxis_transform())

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT, dpi=160, bbox_inches="tight")
print(f"Saved: {OUT}")
