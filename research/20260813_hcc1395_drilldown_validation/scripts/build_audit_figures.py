#!/usr/bin/env python3
"""Build compact audit figures from the machine-readable CSV outputs."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


GREEN = "#176B58"
BLUE = "#2D638D"
ORANGE = "#C56A12"
RED = "#A94336"
GRAY = "#7A817C"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def cohort_figure(rows: list[dict], out: Path) -> None:
    rows = sorted(rows, key=lambda r: float(r["tree_pct_all_regions"]), reverse=True)
    labels = [r["sample"].replace("HCC1395_DORADO", "HCC1395\nDORADO†") for r in rows]
    tree = [float(r["tree_pct_all_regions"]) for r in rows]
    unique = [float(r["unique_pct_among_tree"]) for r in rows]
    colors = [GRAY if r["technical_replicate"] == "True" else GREEN for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.3), sharey=False)
    y = range(len(rows))
    for ax, values, title, xlabel in (
        (axes[0], tree, "Representative-tree coverage", "regions with tree / all regions (%)"),
        (axes[1], unique, "Unique solution among tree-bearing regions", "unique best tree / regions with tree (%)"),
    ):
        bars = ax.barh(list(y), values, color=colors, edgecolor="white", linewidth=0.8)
        ax.set_yticks(list(y), labels)
        ax.invert_yaxis()
        ax.set_xlim(0, 100)
        ax.grid(axis="x", color="#D9DDD9", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.set_title(title, loc="left", fontsize=12, weight="bold")
        ax.set_xlabel(xlabel)
        for bar, value in zip(bars, values):
            ax.text(value + 1.2, bar.get_y() + bar.get_height() / 2, f"{value:.1f}%",
                    va="center", fontsize=9, color="#24312C")
    fig.suptitle("Seven comparable topology datasets — descriptive only", x=0.08, ha="left",
                 fontsize=15, weight="bold")
    fig.text(0.08, 0.015,
             "† HCC1395 DORADO is a technical replicate (biological n remains 6). "
             "Schema/model/AF basis match, but no cross-sample ISM/lineage dashboard is available.",
             fontsize=9, color="#4F5954")
    fig.tight_layout(rect=(0.05, 0.06, 1, 0.94))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def coverage_figure(rows: list[dict], out: Path) -> None:
    order = ["1", "2", "3", "4", "5", "6", "7", "8+"]
    by_bundle = {}
    for row in rows:
        by_bundle.setdefault(row["bundle"], {})[row["k_bin"]] = float(row["percent"])

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    styles = {"v1": (BLUE, "o"), "v3": (ORANGE, "s")}
    for bundle in sorted(by_bundle):
        color, marker = styles.get(bundle, (GRAY, "o"))
        values = [by_bundle[bundle][k] for k in order]
        ax.plot(order, values, marker=marker, linewidth=2.3, markersize=6,
                color=color, label=bundle)
        for x, value in zip(order, values):
            ax.annotate(f"{value:.1f}", (x, value), textcoords="offset points",
                        xytext=(0, 7 if bundle == "v3" else -14), ha="center",
                        fontsize=8, color=color)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Topology max active k at distinct sSNV")
    ax.set_ylabel("sSNV with methylation summary / topology sSNV (%)")
    ax.set_title("ISM availability falls sharply at high k", loc="left", fontsize=15, weight="bold")
    ax.grid(axis="y", color="#D9DDD9", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, ncol=2, loc="lower left")
    ax.text(0, -0.23,
            "Selection-bias diagnostic, not biological effect. v1 uses ±1 kb; v3 uses ±5 kb, so versions are not a controlled A/B.",
            transform=ax.transAxes, fontsize=9, color="#4F5954")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    cohort_figure(read_csv(args.results_dir / "cohort_topology_metrics.csv"),
                  args.output_dir / "02_cohort_topology_overview.png")
    coverage_figure(read_csv(args.results_dir / "methylation_coverage_by_k.csv"),
                    args.output_dir / "03_methylation_coverage_by_k.png")
    print("wrote 02_cohort_topology_overview.png and 03_methylation_coverage_by_k.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

