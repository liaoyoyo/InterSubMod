#!/usr/bin/env python3
"""Export region-level diagnostics tables and plots from existing region outputs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-dir", default=None, help="Existing region directory")
    parser.add_argument("--search-root", default=None, help="Search root when using --region-key")
    parser.add_argument("--region-key", default=None, help="Region key in chr:pos:ref:alt format")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: <region-dir>/diagnostics)")
    parser.add_argument("--metric", default="BERNOULLI", help="Distance metric directory name")
    parser.add_argument("--cluster-k", type=int, default=2, help="Number of clusters for diagnostic assignment")
    return parser.parse_args()


def parse_metadata(path: Path) -> Dict[str, str]:
    info: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("Region ID:"):
                info["region_id"] = line.split(":", 1)[1].strip()
            elif line.startswith("Region:"):
                info["region"] = line.split(":", 1)[1].strip()
            elif line.startswith("SNV Position:"):
                info["snv_position"] = line.split(":", 1)[1].strip()
            elif line.startswith("SNV:"):
                info["snv_change"] = line.split(":", 1)[1].strip().replace(" ", "")
            elif line.startswith("Num Reads:"):
                info["num_reads"] = line.split(":", 1)[1].strip()
            elif line.startswith("Num CpG Sites:"):
                info["num_cpgs"] = line.split(":", 1)[1].strip()
    return info


def region_key_from_metadata(info: Dict[str, str]) -> str:
    chrom, pos = info["snv_position"].split(":")
    ref, alt = [part.strip() for part in info["snv_change"].split("->")]
    return f"{chrom}:{pos}:{ref}:{alt}"


def find_region_dir(search_root: Path, region_key: str) -> Optional[Path]:
    for metadata in search_root.rglob("metadata.txt"):
        info = parse_metadata(metadata)
        if "snv_position" not in info or "snv_change" not in info:
            continue
        if region_key_from_metadata(info) == region_key:
            return metadata.parent
    return None


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df = df.replace("NA", np.nan)
    return df.astype(float)


def load_reads(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "read_id" in df.columns:
        df = df.set_index("read_id")
    return df


def derive_clusters(dist: pd.DataFrame, cluster_k: int) -> np.ndarray:
    if len(dist.index) < 2:
        return np.ones(len(dist.index), dtype=int)
    condensed = squareform(dist.values, checks=False)
    linkage_matrix = linkage(condensed, method="average")
    return fcluster(linkage_matrix, t=max(2, cluster_k), criterion="maxclust")


def plot_methylation_heatmap(meth: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(meth, cmap="viridis", ax=ax, cbar_kws={"label": "Methylation"})
    ax.set_title("Methylation Matrix")
    ax.set_xlabel("CpG")
    ax.set_ylabel("Read")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_distance_heatmap(dist: pd.DataFrame, clusters: np.ndarray, out_png: Path) -> None:
    order = np.argsort(clusters)
    ordered = dist.iloc[order, order]
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(ordered, cmap="mako", ax=ax, cbar_kws={"label": "Distance"})
    ax.set_title("Distance Matrix (cluster ordered)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_label_vs_cluster(reads: pd.DataFrame, clusters: np.ndarray, out_png: Path) -> None:
    data = reads.copy()
    data["cluster"] = clusters

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    hp_ct = pd.crosstab(data["cluster"], data["hp"])
    hp_ct.plot(kind="bar", stacked=True, ax=axes[0], colormap="tab20")
    axes[0].set_title("Cluster vs HP")
    axes[0].set_xlabel("Cluster")

    allele_ct = pd.crosstab(data["cluster"], data["alt_support"])
    allele_ct.plot(kind="bar", stacked=True, ax=axes[1], colormap="Set2")
    axes[1].set_title("Cluster vs Allele Support")
    axes[1].set_xlabel("Cluster")

    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def summarize_region(reads: pd.DataFrame, dist: pd.DataFrame) -> Dict[str, str]:
    hp_assigned = reads["hp"].astype(str).isin(["1", "2", "1-1", "2-1", "3"]).sum()
    allele_assigned = (~reads["alt_support"].astype(str).str.upper().eq("UNKNOWN")).sum()
    valid = dist.to_numpy()
    valid = valid[np.triu_indices_from(valid, k=1)]
    valid = valid[~np.isnan(valid)]
    return {
        "reads_total": str(len(reads.index)),
        "hp_assign_rate": f"{(hp_assigned / len(reads.index)) if len(reads.index) else 0.0:.6f}",
        "allele_assign_rate": f"{(allele_assigned / len(reads.index)) if len(reads.index) else 0.0:.6f}",
        "pairwise_mean_dist": f"{np.mean(valid) if len(valid) else math.nan:.6f}" if len(valid) else "",
        "pairwise_median_dist": f"{np.median(valid) if len(valid) else math.nan:.6f}" if len(valid) else "",
    }


def main() -> None:
    args = parse_args()
    region_dir: Optional[Path] = Path(args.region_dir).resolve() if args.region_dir else None

    if region_dir is None:
        if not args.search_root or not args.region_key:
            raise SystemExit("Provide either --region-dir or both --search-root and --region-key")
        region_dir = find_region_dir(Path(args.search_root).resolve(), args.region_key)
        if region_dir is None:
            raise SystemExit(f"Region {args.region_key} not found under {args.search_root}")

    output_dir = Path(args.output_dir).resolve() if args.output_dir else region_dir / "diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = parse_metadata(region_dir / "metadata.txt")
    reads = load_reads(region_dir / "reads" / "reads.tsv")
    meth = load_matrix(region_dir / "methylation" / "methylation.csv")
    dist = load_matrix(region_dir / "distance" / args.metric / "matrix.csv")
    clusters = derive_clusters(dist, args.cluster_k)

    summary = summarize_region(reads, dist)
    summary_path = output_dir / "region_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "region_key",
                "region_id",
                "reads_total",
                "hp_assign_rate",
                "allele_assign_rate",
                "pairwise_mean_dist",
                "pairwise_median_dist",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "region_key": region_key_from_metadata(metadata),
                "region_id": metadata.get("region_id", ""),
                **summary,
            }
        )

    plot_methylation_heatmap(meth, output_dir / "heatmap_methylation.png")
    plot_distance_heatmap(dist, clusters, output_dir / "heatmap_distance.png")
    plot_label_vs_cluster(reads, clusters, output_dir / "label_vs_cluster.png")

    notes = []
    notes.append(f"# Region Diagnostics: {region_key_from_metadata(metadata)}")
    notes.append("")
    notes.append(f"- Region dir: `{region_dir}`")
    notes.append(f"- Region ID: `{metadata.get('region_id', '')}`")
    notes.append(f"- Reads: `{summary['reads_total']}`")
    notes.append(f"- HP assign rate: `{summary['hp_assign_rate']}`")
    notes.append(f"- Allele assign rate: `{summary['allele_assign_rate']}`")
    notes.append("")
    notes.append("## Outputs")
    notes.append("")
    notes.append("- [region_summary.tsv](region_summary.tsv)")
    notes.append("- [heatmap_methylation.png](heatmap_methylation.png)")
    notes.append("- [heatmap_distance.png](heatmap_distance.png)")
    notes.append("- [label_vs_cluster.png](label_vs_cluster.png)")
    (output_dir / "region_notes.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[export_region_diagnostics] Wrote {summary_path}")
    print(f"[export_region_diagnostics] Wrote {output_dir / 'region_notes.md'}")


if __name__ == "__main__":
    main()
