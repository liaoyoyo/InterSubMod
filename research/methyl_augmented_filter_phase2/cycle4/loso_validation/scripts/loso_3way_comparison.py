#!/usr/bin/env python3
"""3-way LOSO comparison: baseline (10-feature) vs H_NEW_2 (2-feature) vs H_NEW_4 (9-feature drop caller_af)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
LOSO_DIR = REPO / "research/methyl_augmented_filter_phase2/cycle4/loso_validation"
FIG = LOSO_DIR / "figures/loso_3way_comparison.png"

samples = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]

baseline_df = pd.read_csv(LOSO_DIR / "data/loso_cv_results.tsv", sep="\t").set_index("test_sample")
hnew2_df = pd.read_csv(LOSO_DIR / "data/loso_hnew2_results.tsv", sep="\t").set_index("test_sample")
hnew4_df = pd.read_csv(LOSO_DIR / "data/loso_hnew4_results.tsv", sep="\t").set_index("test_sample")

baseline_vals = [float(baseline_df.loc[s, "delta_F1"]) for s in samples]
hnew2_vals = [float(hnew2_df.loc[s, "delta_F1"]) for s in samples]
hnew4_vals = [float(hnew4_df.loc[s, "delta_F1"]) for s in samples]

fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(len(samples))
width = 0.25

bars1 = ax.bar(x - width, baseline_vals, width, label="Baseline LOSO (10 features)",
                color="#6B7280", edgecolor="black", linewidth=0.4)
bars2 = ax.bar(x, hnew2_vals, width, label="H_NEW_2 LOSO (2 features: LOH + HPFineF)",
                color="#1E3A8A", edgecolor="black", linewidth=0.4)
bars3 = ax.bar(x + width, hnew4_vals, width, label="H_NEW_4 LOSO (9 features, drop caller_af)",
                color="#15803D", edgecolor="black", linewidth=0.4)

# Annotate values
for bars, vals in [(bars1, baseline_vals), (bars2, hnew2_vals), (bars3, hnew4_vals)]:
    for b, v in zip(bars, vals):
        if abs(v) > 0.0001:
            ax.text(b.get_x() + b.get_width()/2,
                    v + (0.0003 if v > 0 else -0.0003),
                    f"{v:+.5f}", ha="center", va="bottom" if v > 0 else "top",
                    fontsize=7, fontweight="bold" if abs(v) > 0.005 else "normal")

ax.axhline(0, color="black", linewidth=0.7)
ax.axhline(0.005, color="#A16207", linestyle="--", linewidth=0.8, label="Cohen +0.005 small effect")
ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("LOSO ΔF1 (held-out sample)", fontsize=11)
ax.set_title("3-way LOSO comparison — baseline vs H_NEW_2 vs H_NEW_4\n"
             "Key finding: H_NEW_4 (drop caller_af) HCC1395 = +0.00699 (sanity violated, unexpected)",
             fontsize=12, fontweight="bold")
ax.legend(loc="upper left", fontsize=9)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(FIG, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"✓ {FIG}")
