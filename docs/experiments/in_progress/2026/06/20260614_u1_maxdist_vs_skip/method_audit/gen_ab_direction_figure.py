#!/usr/bin/env python3
"""A/B two-direction comparison figure — reads whole-genome significance_summary, self-computes.

Panel 1: HP-AUC normal vs tumor distribution (B-clean structure-first foundation + tumor weakening).
Panel 2: method sig-rate comparison — B-clustering (over-sensitive) vs A-Δβ / HP-AUC (calibrated).
Panel 3: A vs B 2x2 agreement (germline ASM Δβ [A] vs GlobalTest sig [B]) — stringency mismatch.
All numbers computed from the CSV (no hand-typed values).
"""
import csv
import math
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CSV = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_wg_abc/significance_summary.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/ab_direction_comparison.png"


def f(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


rows = list(csv.DictReader(open(CSV)))
n = len(rows)
auc_n = [f(r["HP_AUC_Normal"]) for r in rows if f(r["HP_AUC_Normal"]) is not None and f(r["HP_AUC_Normal"]) >= 0]
auc_t = [f(r["HP_AUC_Tumor"]) for r in rows if f(r["HP_AUC_Tumor"]) is not None and f(r["HP_AUC_Tumor"]) >= 0]


def sig_rate(pred):
    return sum(1 for r in rows if pred(r)) / n


def pv(r, c):
    return f(r.get(c)) is not None and f(r.get(c)) <= 0.05


rates = {
    "ClusterPERMANOVA\n(B struct, circular)": sum(
        1 for r in rows if r["ClusterPermanovaValid"] == "true" and pv(r, "ClusterPermanovaP")
    ) / max(1, sum(1 for r in rows if r["ClusterPermanovaValid"] == "true")),
    "HPFine 4-grp\n(B struct)": sig_rate(lambda r: pv(r, "HPFineP")),
    "GlobalTest\n(B cluster×HP)": sig_rate(lambda r: pv(r, "GlobalP")),
    "HP-AUC≥0.7 normal\n(B clean)": sum(1 for x in auc_n if x >= 0.7) / max(1, len(auc_n)),
    "germline ASM Δβ\n(A label, calibrated)": sig_rate(lambda r: r.get("GermlineAsmDbeta_Sig") == "true"),
}

# A vs B 2x2
a = [r.get("GermlineAsmDbeta_Sig") == "true" for r in rows]
b = [pv(r, "GlobalP") for r in rows]
both = sum(1 for i in range(n) if a[i] and b[i])
aonly = sum(1 for i in range(n) if a[i] and not b[i])
bonly = sum(1 for i in range(n) if not a[i] and b[i])
neither = sum(1 for i in range(n) if not a[i] and not b[i])

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

# Panel 1: HP-AUC distributions
ax[0].hist(auc_n, bins=40, alpha=0.6, color="#2ca02c", label=f"normal (med {np.median(auc_n):.3f})")
ax[0].hist(auc_t, bins=40, alpha=0.6, color="#e8820e", label=f"tumor (med {np.median(auc_t):.3f})")
ax[0].axvline(0.5, color="gray", ls=":", lw=1)
ax[0].axvline(0.7, color="red", ls="--", lw=1, label="strong=0.7")
ax[0].set_title("B direction (clean): HP-AUC\n= does distance recover germline-HP", fontsize=10)
ax[0].set_xlabel("HP-AUC (0.5=none, 1.0=perfect)")
ax[0].set_ylabel("# regions")
ax[0].legend(fontsize=8)

# Panel 2: sig-rate stringency
names = list(rates.keys())
vals = [rates[k] * 100 for k in names]
colors = ["#d02020", "#d02020", "#f0a020", "#2ca02c", "#1f5fd0"]
ax[1].barh(range(len(names)), vals, color=colors)
ax[1].set_yticks(range(len(names)))
ax[1].set_yticklabels(names, fontsize=8)
ax[1].invert_yaxis()
for i, v in enumerate(vals):
    ax[1].text(v + 1, i, f"{v:.0f}%", va="center", fontsize=8)
ax[1].set_xlim(0, 110)
ax[1].set_xlabel("% regions significant")
ax[1].set_title("Stringency mismatch:\nB-clustering over-sensitive vs A/HP-AUC calibrated", fontsize=10)

# Panel 3: A vs B 2x2
mat = np.array([[both, aonly], [bonly, neither]])
ax[2].imshow(mat, cmap="Blues")
ax[2].set_xticks([0, 1])
ax[2].set_xticklabels(["B sig", "B not"], fontsize=9)
ax[2].set_yticks([0, 1])
ax[2].set_yticklabels(["A sig", "A not"], fontsize=9)
for (i, j), v in np.ndenumerate(mat):
    ax[2].text(j, i, f"{v}\n({v/n*100:.0f}%)", ha="center", va="center",
               color="white" if v > mat.max() / 2 else "black", fontsize=10)
agree = (both + neither) / n * 100
ax[2].set_title(f"A (germline Δβ) vs B (GlobalTest)\nagreement {agree:.0f}% — B-only {bonly} = B over-calls", fontsize=10)

fig.suptitle(f"ISM A label-first vs B structure-first — whole-genome {n} regions (HCC1395 paired)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT, dpi=120, bbox_inches="tight")
print(f"saved {OUT}; HP-AUC normal med {np.median(auc_n):.3f} tumor {np.median(auc_t):.3f}; "
      f"A/B agree {agree:.0f}% both={both} Aonly={aonly} Bonly={bonly} neither={neither}")
