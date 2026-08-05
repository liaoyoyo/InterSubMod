#!/usr/bin/env python3
"""Compare truth-confirmed raw-VAF distributions across all seven datasets.

This post-process intentionally compares marginal distributions rather than exact
loci.  It is a technical profile and cannot establish shared clone identity.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import platform
import shlex
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import jensenshannon, squareform


SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
SOURCE_ORDER = ["latest_lps_pass", "baseline_caller_pass"]
METRIC_ORDER = ["js_divergence_50bin_nats", "wasserstein_vaf", "ks_statistic"]
METRIC_TITLES = {
    "js_divergence_50bin_nats": "Jensen–Shannon divergence (50 bins, nats)",
    "wasserstein_vaf": "Wasserstein distance (raw-VAF units)",
    "ks_statistic": "Kolmogorov–Smirnov effect size D",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260712)
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def pair_metrics(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    bins = np.linspace(0.0, 1.0, 51)
    hist_x, _ = np.histogram(x, bins=bins)
    hist_y, _ = np.histogram(y, bins=bins)
    js_divergence = float(
        jensenshannon(hist_x.astype(float) + 1e-12, hist_y.astype(float) + 1e-12, base=np.e) ** 2
    )
    return {
        "js_divergence_50bin_nats": js_divergence,
        "wasserstein_vaf": float(stats.wasserstein_distance(x, y)),
        "ks_statistic": float(stats.ks_2samp(x, y, alternative="two-sided", method="auto").statistic),
    }


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".png"), dpi=180, bbox_inches="tight", metadata={"Software": "InterSubMod raw VAF pairwise postprocess"})
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def plot_heatmaps(matrices: Mapping[Tuple[str, str], pd.DataFrame], output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for row_index, source in enumerate(SOURCE_ORDER):
        for column_index, metric in enumerate(METRIC_ORDER):
            axis = axes[row_index, column_index]
            matrix = matrices[(source, metric)].loc[SAMPLE_ORDER, SAMPLE_ORDER]
            artist = axis.imshow(matrix.to_numpy(dtype=float), cmap="viridis")
            axis.set_xticks(range(len(SAMPLE_ORDER)), SAMPLE_ORDER, rotation=45, ha="right", fontsize=8)
            axis.set_yticks(range(len(SAMPLE_ORDER)), SAMPLE_ORDER, fontsize=8)
            axis.set_title(METRIC_TITLES[metric], fontsize=10)
            for i in range(len(SAMPLE_ORDER)):
                for j in range(len(SAMPLE_ORDER)):
                    value = float(matrix.iloc[i, j])
                    axis.text(j, i, f"{value:.3f}", ha="center", va="center", fontsize=6.5,
                              color="white" if value > matrix.to_numpy().max() * 0.48 else "black")
            fig.colorbar(artist, ax=axis, fraction=0.046, pad=0.04)
            if column_index == 0:
                axis.set_ylabel("Latest LongPhase-S PASS" if source == "latest_lps_pass" else "Caller-PASS baseline")
    fig.suptitle("Pairwise truth-confirmed marginal raw-VAF distribution distances\nCross-cell-line similarity is descriptive, not clone identity", fontsize=14)
    fig.tight_layout()
    save_figure(fig, output_dir / "all_sample_truth_confirmed_pairwise_distance_heatmaps")


def plot_dendrograms(
    linkage_records: Mapping[Tuple[str, str], np.ndarray], output_dir: Path
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    for row_index, source in enumerate(SOURCE_ORDER):
        for column_index, metric in enumerate(METRIC_ORDER):
            axis = axes[row_index, column_index]
            linkage = linkage_records[(source, metric)]
            hierarchy.dendrogram(linkage, labels=SAMPLE_ORDER, leaf_rotation=40, leaf_font_size=8, ax=axis, color_threshold=None)
            title = METRIC_TITLES[metric]
            if metric == "js_divergence_50bin_nats":
                title += "\n(clustered on sqrt(JSD))"
            axis.set_title(title, fontsize=9)
            axis.set_ylabel("Average-linkage height")
            axis.grid(axis="y", alpha=0.2)
            if column_index == 0:
                axis.text(-0.20, 0.5, "Latest" if source == "latest_lps_pass" else "Baseline", transform=axis.transAxes,
                          rotation=90, va="center", ha="center", fontsize=11, weight="bold")
    fig.suptitle("Hierarchical clustering of truth-confirmed raw-VAF distribution profiles", fontsize=14)
    fig.tight_layout()
    save_figure(fig, output_dir / "all_sample_truth_confirmed_pairwise_dendrograms")


def main() -> None:
    args = parse_args()
    records_path = args.records.resolve()
    output_dir = args.output_dir.resolve()
    data_dir = output_dir / "data"
    figure_dir = output_dir / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.hashsalt": str(args.seed)})

    records = pd.read_csv(
        records_path,
        sep="\t",
        compression="gzip",
        usecols=["sample", "source", "af", "truth_confirmed"],
        dtype={"sample": "category", "source": "category", "af": "float64", "truth_confirmed": "bool"},
    )
    records = records[records["truth_confirmed"]].copy()
    expected_pairs = {(source, sample) for source in SOURCE_ORDER for sample in SAMPLE_ORDER}
    observed_pairs = set(zip(records["source"].astype(str), records["sample"].astype(str)))
    if observed_pairs != expected_pairs:
        raise RuntimeError(f"Dataset/source coverage mismatch: missing={sorted(expected_pairs-observed_pairs)} extra={sorted(observed_pairs-expected_pairs)}")

    arrays: Dict[Tuple[str, str], np.ndarray] = {}
    sample_sizes: List[Dict[str, Any]] = []
    for source in SOURCE_ORDER:
        for sample in SAMPLE_ORDER:
            values = records[(records["source"] == source) & (records["sample"] == sample)]["af"].to_numpy(dtype=float)
            if not len(values):
                raise RuntimeError(f"Empty truth-confirmed distribution: {source}/{sample}")
            arrays[(source, sample)] = values
            sample_sizes.append({"source": source, "sample": sample, "truth_confirmed_n": int(len(values))})

    long_rows: List[Dict[str, Any]] = []
    matrices: Dict[Tuple[str, str], pd.DataFrame] = {}
    linkage_records: Dict[Tuple[str, str], np.ndarray] = {}
    nearest_rows: List[Dict[str, Any]] = []
    clustering_rows: List[Dict[str, Any]] = []
    hcc_summaries: List[Dict[str, Any]] = []
    for source in SOURCE_ORDER:
        matrix_values = {metric: np.zeros((len(SAMPLE_ORDER), len(SAMPLE_ORDER)), dtype=float) for metric in METRIC_ORDER}
        for i, sample_a in enumerate(SAMPLE_ORDER):
            for j in range(i + 1, len(SAMPLE_ORDER)):
                sample_b = SAMPLE_ORDER[j]
                values = pair_metrics(arrays[(source, sample_a)], arrays[(source, sample_b)])
                for metric, value in values.items():
                    matrix_values[metric][i, j] = value
                    matrix_values[metric][j, i] = value
                    long_rows.append(
                        {
                            "source": source,
                            "sample_a": sample_a,
                            "sample_b": sample_b,
                            "metric": metric,
                            "distance": value,
                            "n_a": int(len(arrays[(source, sample_a)])),
                            "n_b": int(len(arrays[(source, sample_b)])),
                            "distribution_scope": "truth_confirmed_marginal_raw_vaf",
                        }
                    )
        for metric in METRIC_ORDER:
            matrix = pd.DataFrame(matrix_values[metric], index=SAMPLE_ORDER, columns=SAMPLE_ORDER)
            matrix.index.name = "sample"
            matrices[(source, metric)] = matrix
            matrix.to_csv(data_dir / f"pairwise_{metric}_{source}_matrix.tsv", sep="\t", float_format="%.10g")

            # sqrt(JSD) is a metric; the other two matrices are used directly.
            clustering_matrix = np.sqrt(matrix_values[metric]) if metric == "js_divergence_50bin_nats" else matrix_values[metric]
            condensed = squareform(clustering_matrix, checks=True)
            linkage = hierarchy.linkage(condensed, method="average", optimal_ordering=True)
            linkage_records[(source, metric)] = linkage
            leaves = hierarchy.leaves_list(linkage)
            leaf_names = [SAMPLE_ORDER[int(index)] for index in leaves]
            for merge_index, row in enumerate(linkage):
                clustering_rows.append(
                    {
                        "source": source,
                        "metric": metric,
                        "linkage_method": "average",
                        "jsd_transform_for_linkage": "sqrt" if metric == "js_divergence_50bin_nats" else "none",
                        "merge_index": merge_index,
                        "left_cluster": int(row[0]),
                        "right_cluster": int(row[1]),
                        "height": float(row[2]),
                        "member_count": int(row[3]),
                        "leaf_order": "|".join(leaf_names),
                    }
                )

            for i, sample in enumerate(SAMPLE_ORDER):
                distances = [(SAMPLE_ORDER[j], float(matrix_values[metric][i, j])) for j in range(len(SAMPLE_ORDER)) if j != i]
                distances.sort(key=lambda item: (item[1], item[0]))
                for rank, (neighbor, distance) in enumerate(distances, start=1):
                    nearest_rows.append(
                        {
                            "source": source,
                            "metric": metric,
                            "sample": sample,
                            "rank": rank,
                            "neighbor": neighbor,
                            "distance": distance,
                            "is_nearest_neighbor": rank == 1,
                        }
                    )

            hcc_index = SAMPLE_ORDER.index("HCC1395")
            dorado_index = SAMPLE_ORDER.index("HCC1395_DORADO")
            dorado_distance = float(matrix_values[metric][hcc_index, dorado_index])
            all_peers = [
                (SAMPLE_ORDER[j], float(matrix_values[metric][hcc_index, j]))
                for j in range(len(SAMPLE_ORDER))
                if j != hcc_index
            ]
            all_peers.sort(key=lambda item: (item[1], item[0]))
            dorado_rank = next(rank for rank, (name, _) in enumerate(all_peers, start=1) if name == "HCC1395_DORADO")
            biological_others = [(name, value) for name, value in all_peers if name != "HCC1395_DORADO"]
            nearest_other, nearest_other_distance = biological_others[0]
            median_other_distance = float(np.median([value for _, value in biological_others]))
            hcc_summaries.append(
                {
                    "source": source,
                    "metric": metric,
                    "hcc1395_nearest_neighbor": all_peers[0][0],
                    "hcc1395_nearest_neighbor_distance": all_peers[0][1],
                    "dorado_rank_among_6_peers": dorado_rank,
                    "dorado_distance": dorado_distance,
                    "nearest_other_cell_line": nearest_other,
                    "nearest_other_cell_line_distance": nearest_other_distance,
                    "nearest_other_over_dorado_distance_ratio": nearest_other_distance / dorado_distance if dorado_distance > 0 else float("inf"),
                    "median_other_cell_line_distance": median_other_distance,
                    "median_other_over_dorado_distance_ratio": median_other_distance / dorado_distance if dorado_distance > 0 else float("inf"),
                }
            )

    long_frame = pd.DataFrame(long_rows).sort_values(["source", "metric", "sample_a", "sample_b"])
    nearest_frame = pd.DataFrame(nearest_rows).sort_values(["source", "metric", "sample", "rank"])
    clustering_frame = pd.DataFrame(clustering_rows).sort_values(["source", "metric", "merge_index"])
    hcc_frame = pd.DataFrame(hcc_summaries).sort_values(["source", "metric"])
    sample_size_frame = pd.DataFrame(sample_sizes).sort_values(["source", "sample"])
    long_frame.to_csv(data_dir / "pairwise_truth_confirmed_distribution_distances.tsv", sep="\t", index=False, float_format="%.10g")
    nearest_frame.to_csv(data_dir / "pairwise_truth_confirmed_nearest_neighbor_rankings.tsv", sep="\t", index=False, float_format="%.10g")
    clustering_frame.to_csv(data_dir / "pairwise_truth_confirmed_hierarchical_clustering.tsv", sep="\t", index=False, float_format="%.10g")
    hcc_frame.to_csv(data_dir / "hcc1395_distribution_neighbor_summary.tsv", sep="\t", index=False, float_format="%.10g")
    sample_size_frame.to_csv(data_dir / "pairwise_distribution_sample_sizes.tsv", sep="\t", index=False)

    plot_heatmaps(matrices, figure_dir)
    plot_dendrograms(linkage_records, figure_dir)

    output_hashes = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "pairwise_distribution_analysis.json":
            output_hashes[str(path.relative_to(output_dir))] = {
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    payload = {
        "schema": "intersubmod.raw_vaf_pairwise_distribution.v1",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "input": {"path": str(records_path), "sha256": sha256_file(records_path)},
        "scope": {
            "sources": SOURCE_ORDER,
            "datasets": SAMPLE_ORDER,
            "records": "truth-confirmed autosomal biallelic SNVs",
            "comparison": "marginal raw-VAF distributions; loci need not be shared",
        },
        "metric_definitions": {
            "js_divergence_50bin_nats": "Jensen-Shannon divergence on 50 equal-width bins over [0,1], natural-log units",
            "wasserstein_vaf": "first Wasserstein distance on raw VAF values",
            "ks_statistic": "two-sample Kolmogorov-Smirnov D statistic (effect size)",
            "hierarchical_clustering": "average linkage; sqrt(JSD) used for JSD clustering because Jensen-Shannon distance is sqrt(divergence)",
        },
        "hcc1395_neighbor_summary": hcc_summaries,
        "claim_ceiling": {
            "level": "marginal_raw_vaf_distribution_similarity_only",
            "allowed": "rank technical distribution similarity and test whether the same-cell-line DORADO profile is unusually close",
            "not_allowed": "claim shared clone/subclone truth or shared phylogeny from marginal raw VAF",
            "major_caveats": [
                "raw VAF is not corrected for purity, copy number, LOH, or mutant multiplicity",
                "truth definitions differ across cell lines",
                "marginal distribution similarity ignores mutation identity",
                "HCC1395 and HCC1395_DORADO are technical datasets from one cell line, not biological replicates",
            ],
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "outputs": output_hashes,
    }
    (output_dir / "pairwise_distribution_analysis.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    )

    print("PASS pairwise truth-confirmed marginal distributions")
    print(hcc_frame.to_string(index=False))
    print(f"OUTPUT {output_dir}")


if __name__ == "__main__":
    main()
