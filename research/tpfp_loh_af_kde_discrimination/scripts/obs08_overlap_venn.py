#!/usr/bin/env python3
"""
obs08 | Region-level overlap between S3, Top-17 and Top-28 cells (HCC1395 TO).

Question: how much of S3 is already inside Top-17?  Do Top-17 and Top-28
add net-new regions or only slice S3 further?  Implemented as quantitative
matrix + proportional concentric bars rather than classic Venn.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (DATA_DIR, PALETTE, apply_style, build_scheme_masks,
                         ensure_fig_dir, load_master)


def build_top_k_mask(df: pd.DataFrame, pareto_path: Path, k: int) -> pd.Series:
    pareto = pd.read_csv(pareto_path, sep="\t")
    top = pareto.head(k)[["LOH_Subtype", "AF_class", "cn_tier_F", "HPFineNGroups", "nr_band"]]
    edges = [0.65, 0.99, 1.33, 1.82]
    tiers = pd.cut(df["CovM_used"], bins=[-np.inf] + edges + [np.inf],
                    labels=[f"T{i}" for i in range(len(edges) + 1)], include_lowest=True).astype(str)
    loh = df["LOH_Subtype"].fillna("None")
    nr = df["NumReads"].apply(lambda x: "low" if x < 60 else ("mid" if x < 120 else "high"))
    key = pd.DataFrame({
        "LOH_Subtype": loh, "AF_class": df["AF_class"],
        "cn_tier_F": tiers, "HPFineNGroups": df["HPFineNGroups"],
        "nr_band": nr,
    })
    top_set = set(tuple(r) for r in top.astype(object).itertuples(index=False, name=None))
    mask = key.apply(lambda r: tuple(r) in top_set, axis=1)
    return mask


def main() -> None:
    apply_style()
    fig_dir = ensure_fig_dir()
    df = load_master()

    to = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")].copy()

    masks = build_scheme_masks(to)
    s3_mask = masks["S3"].values
    s5_mask = masks["S5"].values
    t17_mask = build_top_k_mask(to, DATA_DIR / "tpfp_5d_pareto_HCC1395_TO.tsv", 17).values
    t28_mask = build_top_k_mask(to, DATA_DIR / "tpfp_5d_pareto_HCC1395_TO.tsv", 28).values

    groups = {"S3": s3_mask, "S5": s5_mask, "Top-17": t17_mask, "Top-28": t28_mask}
    sizes = {k: int(v.sum()) for k, v in groups.items()}
    tp_counts = {k: int((v & (to["tp_label"].values == 1)).sum()) for k, v in groups.items()}

    # Pairwise overlap matrix (row=A, col=B -> |A and B|)
    labs = list(groups.keys())
    M = np.zeros((len(labs), len(labs)))
    for i, a in enumerate(labs):
        for j, b in enumerate(labs):
            M[i, j] = int((groups[a] & groups[b]).sum())

    # Build 2-panel: heatmap + bar with TP/FP breakdown
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: pairwise overlap heatmap (absolute counts)
    ax = axes[0]
    im = ax.imshow(M, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(range(len(labs))); ax.set_yticks(range(len(labs)))
    ax.set_xticklabels(labs); ax.set_yticklabels(labs)
    for i in range(len(labs)):
        for j in range(len(labs)):
            val = int(M[i, j])
            frac_a = val / max(sizes[labs[i]], 1)
            ax.text(j, i, f"{val:,}\n({frac_a:.0%}A)", ha="center", va="center",
                     fontsize=10, color=PALETTE["dark"] if M[i, j] < M.max() * 0.6 else "white")
    ax.set_title("Panel A | |A ∩ B| count (row=A, column=B)\nPercent = |A∩B| / |A|",
                  fontsize=11, fontweight="bold", color=PALETTE["dark"], loc="left")
    fig.colorbar(im, ax=ax, shrink=0.7)

    # Panel B: TP/FP bar per scheme
    ax = axes[1]
    xs = np.arange(len(labs))
    tp = np.array([tp_counts[k] for k in labs])
    fp = np.array([sizes[k] - tp_counts[k] for k in labs])
    ax.bar(xs, tp, color=PALETTE["tp"], edgecolor=PALETTE["dark"], label="TP")
    ax.bar(xs, fp, bottom=tp, color=PALETTE["fp"], edgecolor=PALETTE["dark"], label="FP")
    for i, k in enumerate(labs):
        tot = sizes[k]
        purity = tp_counts[k] / tot if tot > 0 else 0
        ax.text(i, tot + 20, f"{tot:,}\np={purity:.3f}",
                 ha="center", fontsize=9, color=PALETTE["dark"])
    ax.set_xticks(xs); ax.set_xticklabels(labs)
    ax.set_ylabel("region count")
    ax.set_title("Panel B | Set size + TP/FP split + purity",
                  fontsize=11, fontweight="bold", color=PALETTE["dark"], loc="left")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=9, loc="upper right")

    fig.suptitle("obs08 | HCC1395 TO — region-level overlap between S3, S5, Top-17, Top-28",
                  fontsize=14, fontweight="bold", color=PALETTE["dark"], y=1.01)
    fig.tight_layout()
    out = fig_dir / "obs08_overlap_venn_schemes.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs08] wrote {out}")

    # also write a TSV of the overlap matrix for the report
    tsv_path = DATA_DIR.parent / "tables" / "obs08_scheme_overlap_matrix.tsv"
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(M, index=labs, columns=labs).astype(int).to_csv(tsv_path, sep="\t")
    print(f"[obs08] wrote {tsv_path}")


if __name__ == "__main__":
    main()
