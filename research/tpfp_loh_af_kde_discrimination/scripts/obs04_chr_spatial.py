#!/usr/bin/env python3
"""
obs04 | Chromosome spatial distribution — where do S3 and Top-17 cells live?

Question: are the high-TP 5D cells clustered on specific chromosomes or
scattered?  Clustering might hint at biological hotspots (e.g. HLA, centromere)
versus general structural signals.

Layout: HCC1395 TO variants colored by membership in {S3, S5, Top-17, Top-28},
plotted as a strip by chromosome.
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


CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def build_top_k_mask(df: pd.DataFrame, pareto_path: Path, k: int) -> pd.Series:
    pareto = pd.read_csv(pareto_path, sep="\t")
    top = pareto.head(k)[[
        "LOH_Subtype", "AF_class", "cn_tier_F", "HPFineNGroups", "nr_band"
    ]]
    # Recompute the same 5-D categorical columns on df
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
    to = to[to["Chr"].isin(CHROM_ORDER)].copy()

    masks = build_scheme_masks(to)
    s3_mask = masks["S3"]
    s5_mask = masks["S5"]
    top17 = build_top_k_mask(to, DATA_DIR / "tpfp_5d_pareto_HCC1395_TO.tsv", 17)
    top28 = build_top_k_mask(to, DATA_DIR / "tpfp_5d_pareto_HCC1395_TO.tsv", 28)

    to["_chr_idx"] = to["Chr"].map({c: i for i, c in enumerate(CHROM_ORDER)})
    to["_tp"] = to["tp_label"] == 1

    fig, ax = plt.subplots(figsize=(16, 8))

    # Background: all HCC1395 TO variants
    ax.scatter(to["_chr_idx"] + np.random.uniform(-0.3, 0.3, len(to)),
                np.log10(to["Pos"].clip(1)),
                s=1.5, c=PALETTE["grey"], alpha=0.12, label=f"all TO variants n={len(to):,}")

    # S3 cells
    s3 = to[s3_mask]
    ax.scatter(s3["_chr_idx"] + np.random.uniform(-0.3, 0.3, len(s3)),
                np.log10(s3["Pos"].clip(1)),
                s=8, c=PALETTE["blue"], alpha=0.9, label=f"S3 (Diploid Het) n={len(s3):,}")

    # Top-17
    t17 = to[top17]
    ax.scatter(t17["_chr_idx"] + np.random.uniform(-0.3, 0.3, len(t17)),
                np.log10(t17["Pos"].clip(1)),
                s=14, c=PALETTE["green"], alpha=0.9, label=f"Top-17 5D envelope n={len(t17):,}",
                edgecolors=PALETTE["dark"], linewidths=0.3)

    # Top-28 additional
    t28_only = to[top28 & (~top17)]
    ax.scatter(t28_only["_chr_idx"] + np.random.uniform(-0.3, 0.3, len(t28_only)),
                np.log10(t28_only["Pos"].clip(1)),
                s=10, c=PALETTE["accent"], alpha=0.75,
                label=f"Top-28 \\ Top-17 n={len(t28_only):,}")

    ax.set_xticks(range(len(CHROM_ORDER)))
    ax.set_xticklabels(CHROM_ORDER, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("chromosome")
    ax.set_ylabel("log10(Pos)")
    ax.set_title("obs04 | HCC1395 TO — chromosome spatial distribution of S3 / Top-17 / Top-28",
                  fontsize=13, fontweight="bold", color=PALETTE["dark"])
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, axis="x", alpha=0.2, linestyle=":")

    fig.tight_layout()
    out = fig_dir / "obs04_chr_spatial_scheme.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs04] wrote {out}")


if __name__ == "__main__":
    main()
