#!/usr/bin/env python3
"""V6-C Phase B figure: NG_off × NG_on cross-tab heatmap with TP rate annotation.

Output: research/paired_priority_bug_audit/v6c_phaseB_runs/figures/v6c_phaseB_ng_crosstab.png
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

# Inject CJK font (per project rule on Matplotlib CJK)
for cand in [
    "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf",
    "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
]:
    if Path(cand).exists():
        font_manager.fontManager.addfont(cand)

plt.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "WenQuanYi Micro Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
SUMMARY = REPO / "research/paired_priority_bug_audit/v6c_phaseB_runs/cross_tab_summary.tsv"
FIGURE_DIR = REPO / "research/paired_priority_bug_audit/v6c_phaseB_runs/figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def load_summary() -> list[dict]:
    rows: list[dict] = []
    with SUMMARY.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            r["ng_off"] = int(r["ng_off"])
            r["ng_on"] = int(r["ng_on"])
            r["n_reads"] = int(r["n_reads"])
            rows.append(r)
    return rows


def main() -> int:
    rows = load_summary()
    if not rows:
        return 1

    cell_tp: dict[tuple[int, int], int] = defaultdict(int)
    cell_fp: dict[tuple[int, int], int] = defaultdict(int)
    for r in rows:
        key = (r["ng_off"], r["ng_on"])
        if r["label"] == "TP":
            cell_tp[key] += 1
        else:
            cell_fp[key] += 1

    ng_offs = sorted({k[0] for k in set(cell_tp) | set(cell_fp)})
    ng_ons = sorted({k[1] for k in set(cell_tp) | set(cell_fp)})

    # Build matrices: TP count, FP count, TP rate
    tp_mat = np.zeros((len(ng_offs), len(ng_ons)))
    fp_mat = np.zeros((len(ng_offs), len(ng_ons)))
    rate_mat = np.full((len(ng_offs), len(ng_ons)), np.nan)
    for i, noff in enumerate(ng_offs):
        for j, non in enumerate(ng_ons):
            tp = cell_tp.get((noff, non), 0)
            fp = cell_fp.get((noff, non), 0)
            tp_mat[i, j] = tp
            fp_mat[i, j] = fp
            if tp + fp > 0:
                rate_mat[i, j] = tp / (tp + fp)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel 1: TP count heatmap
    ax = axes[0]
    im0 = ax.imshow(tp_mat, cmap="Greens", aspect="auto")
    ax.set_xticks(range(len(ng_ons)))
    ax.set_xticklabels(ng_ons)
    ax.set_yticks(range(len(ng_offs)))
    ax.set_yticklabels(ng_offs)
    ax.set_xlabel("NG_on (flag=on)")
    ax.set_ylabel("NG_off (flag=off)")
    ax.set_title("TP region count\n(chr19 HCC1395 ClairS-TO)")
    for i in range(len(ng_offs)):
        for j in range(len(ng_ons)):
            v = tp_mat[i, j]
            if v > 0:
                color = "white" if v > tp_mat.max() * 0.5 else "black"
                ax.text(j, i, f"{int(v)}", ha="center", va="center", color=color, fontsize=11)
    plt.colorbar(im0, ax=ax, fraction=0.046, pad=0.04)

    # Panel 2: FP count heatmap
    ax = axes[1]
    im1 = ax.imshow(fp_mat, cmap="Reds", aspect="auto")
    ax.set_xticks(range(len(ng_ons)))
    ax.set_xticklabels(ng_ons)
    ax.set_yticks(range(len(ng_offs)))
    ax.set_yticklabels(ng_offs)
    ax.set_xlabel("NG_on (flag=on)")
    ax.set_ylabel("NG_off (flag=off)")
    ax.set_title("FP region count\n(chr19 HCC1395 ClairS-TO)")
    for i in range(len(ng_offs)):
        for j in range(len(ng_ons)):
            v = fp_mat[i, j]
            if v > 0:
                color = "white" if v > fp_mat.max() * 0.5 else "black"
                ax.text(j, i, f"{int(v)}", ha="center", va="center", color=color, fontsize=11)
    plt.colorbar(im1, ax=ax, fraction=0.046, pad=0.04)

    # Panel 3: TP rate heatmap with 0.85 gate annotation
    ax = axes[2]
    im2 = ax.imshow(rate_mat, cmap="RdYlGn", aspect="auto", vmin=0.5, vmax=1.0)
    ax.set_xticks(range(len(ng_ons)))
    ax.set_xticklabels(ng_ons)
    ax.set_yticks(range(len(ng_offs)))
    ax.set_yticklabels(ng_offs)
    ax.set_xlabel("NG_on (flag=on)")
    ax.set_ylabel("NG_off (flag=off)")
    ax.set_title("TP rate per cell\n(green = pass 0.85 gate)")
    for i in range(len(ng_offs)):
        for j in range(len(ng_ons)):
            v = rate_mat[i, j]
            if not np.isnan(v):
                tp = int(tp_mat[i, j])
                fp = int(fp_mat[i, j])
                gate_marker = " ✓" if v >= 0.85 else ""
                color = "white" if v < 0.7 else "black"
                ax.text(
                    j, i,
                    f"{v:.2f}{gate_marker}\n({tp}/{tp + fp})",
                    ha="center", va="center", color=color, fontsize=9,
                )
    plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(
        "V6-C Phase B chr19 — NG_off × NG_on cross-tab (germline-hp-only flag off vs on)\n"
        "Marker 真實性測試: NG≥3 marker filter TP rate=94.7% (off) → 91.5% (on, NG=2 cell), "
        "兩者 ≥0.85 gate ✓ POSITIVE",
        fontsize=13,
    )
    plt.tight_layout(rect=(0, 0, 1, 0.95))

    out_path = FIGURE_DIR / "v6c_phaseB_ng_crosstab.png"
    plt.savefig(out_path, dpi=140, bbox_inches="tight")
    print(f"[V6-C Phase B figure] wrote {out_path}")

    # Second figure: bucket schema collapse (per-bucket read counts before/after)
    fig2, ax = plt.subplots(figsize=(10, 5))
    buckets_off = {"1": 8604, "2": 12289, "1-1": 5040, "2-1": 8924, "3": 2599, "0": 5609}
    buckets_on = {"1": 8604, "2": 12289, "1-1": 0, "2-1": 0, "3": 0, "0": 22172}
    keys = ["1", "2", "1-1", "2-1", "3", "0"]
    x = np.arange(len(keys))
    width = 0.4
    off_vals = [buckets_off[k] for k in keys]
    on_vals = [buckets_on[k] for k in keys]
    ax.bar(x - width / 2, off_vals, width, label="flag=off (original schema)", color="steelblue")
    ax.bar(x + width / 2, on_vals, width, label="flag=on (somatic-tag demoted to 0)", color="lightcoral")
    ax.set_xticks(x)
    ax.set_xticklabels([f"hp={k}" for k in keys])
    ax.set_ylabel("Read count (chr19 全部 reads)")
    ax.set_title(
        "Bucket schema collapse — flag=on 下 hp=1-1/2-1/3 全 demote 為 0\n"
        "16,563 reads shifted (= 8,924 + 5,040 + 2,599)，conservation OK"
    )
    ax.legend()
    for i, (off_v, on_v) in enumerate(zip(off_vals, on_vals)):
        if off_v != on_v:
            ax.annotate(
                f"Δ{on_v - off_v:+d}",
                xy=(i, max(off_v, on_v) + 500),
                ha="center", fontsize=10, color="darkred", fontweight="bold",
            )
    plt.tight_layout()
    out_path2 = FIGURE_DIR / "v6c_phaseB_bucket_collapse.png"
    plt.savefig(out_path2, dpi=140, bbox_inches="tight")
    print(f"[V6-C Phase B figure] wrote {out_path2}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
