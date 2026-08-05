#!/usr/bin/env python3
"""Build a source-backed H2009 terminal-ALT methyl-distance and UPGMA HTML artifact."""

from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import json
import math
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


REPO = Path(__file__).resolve().parents[3]
WORK = Path(__file__).resolve().parents[1]
FIGURES = WORK / "figures"
DATA = WORK / "data"
RESULTS = WORK / "results"

SOURCE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
M1_ROOT = SOURCE_ROOT / "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
SITE_RESULTS = M1_ROOT / "all_ssnv_site_results.tsv.gz"
ASSIGNMENTS = M1_ROOT / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
REGION = (
    SOURCE_ROOT
    / "intersubmod_all_ssnv_v2_verification_fix/H2009/"
    "H2009.longphase_s.recalibrated.pass.autosomal_biallelic_snv/"
    "chr18/chr18_563687/chr18_558687_568687"
)
READS = REGION / "reads/reads.tsv"
DISTANCE = REGION / "distance/BERNOULLI/matrix.csv"
METHYLATION = REGION / "methylation/methylation.csv"
TOPOLOGY = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/H2009/H2009.topology.jsonl"
)
FACTORIZATION = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260725_exact_ps_candidate_factorization/all7_v1/H2009.candidate_factorization.jsonl"
)

FIGURE_OUT = FIGURES / "01_h2009_chr18_563687_alt_read_distance_upgma.png"
LEAF_ORDER_OUT = DATA / "20260730_H2009_chr18_563687_UPGMA_leaf_order_01.tsv"
GROUP_DISTANCE_OUT = DATA / "20260730_H2009_chr18_563687_group_distance_summary_01.tsv"
SUMMARY_OUT = RESULTS / "20260730_H2009_chr18_563687建圖驗證摘要_01.json"
STAGING_DB = RESULTS / "artifact_staging.sqlite"
ARTIFACT_OUT = WORK / "artifact.json"

EXPECTED_PRIMARY_SHA = {
    READS: "8c79bfb9b1afd682771c11ad47415c11712461c42b6e03b3f7ca157b45c51855",
    DISTANCE: "1af73f6537609ae66390ae3a5ee353f0539232a0927810b2ac98eb3642633acf",
    METHYLATION: "a6e419afbe3f449f4fb116380c5600054ddcee942412b12b9e47a995624da1af",
}
EXPECTED_ASSIGNMENT_SHA = "82af076d8ce70c66c7f75c4ed4bdeacb7d4c5c43d0905859303b121df483a4ba"
EXPECTED_SITE_RESULTS_SHA = "a8871af3a8c3955bf31aec5eeef0c93aca0683f52cf6d6f1e06fbbb713324f74"

GROUP_META = {
    "1-1": {"alias": "MG-A", "prefix": "A", "color": "#167D78", "expected_n": 81},
    "1-2": {"alias": "MG-B", "prefix": "B", "color": "#D0832D", "expected_n": 5},
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def load_tsv_exact_row(path: Path) -> dict[str, str]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        matches = [
            row
            for row in reader
            if row["dataset"] == "H2009" and row["chrom"] == "chr18" and int(row["pos"]) == 563687
        ]
    if len(matches) != 1:
        raise AssertionError(f"Expected one H2009 chr18:563687 site row, found {len(matches)}")
    return matches[0]


def load_jsonl_exact(path: Path, predicate: Callable[[dict[str, Any]], bool], compressed: bool = False) -> dict:
    opener = gzip.open if compressed else Path.open
    mode = "rt"
    matches: list[dict] = []
    if compressed:
        handle_context = opener(path, mode)
    else:
        handle_context = opener(path, mode, encoding="utf-8")
    with handle_context as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if predicate(row):
                matches.append(row)
    if len(matches) != 1:
        raise AssertionError(f"Expected one JSONL match in {path}, found {len(matches)}")
    return matches[0]


def load_inputs() -> tuple[dict, dict, dict, dict, pd.DataFrame, np.ndarray]:
    for path, expected in EXPECTED_PRIMARY_SHA.items():
        observed = sha256_file(path)
        if observed != expected:
            raise AssertionError(f"Primary artifact SHA mismatch: {path}: {observed} != {expected}")

    assignment_sha = sha256_file(ASSIGNMENTS)
    site_sha = sha256_file(SITE_RESULTS)
    if assignment_sha != EXPECTED_ASSIGNMENT_SHA:
        raise AssertionError(f"Assignment source SHA mismatch: {assignment_sha}")
    if site_sha != EXPECTED_SITE_RESULTS_SHA:
        raise AssertionError(f"Site source SHA mismatch: {site_sha}")

    site = load_tsv_exact_row(SITE_RESULTS)
    assignment = load_jsonl_exact(
        ASSIGNMENTS,
        lambda row: row.get("dataset") == "H2009"
        and row.get("chrom") == "chr18"
        and int(row.get("pos", -1)) == 563687,
        compressed=True,
    )
    topology = load_jsonl_exact(
        TOPOLOGY,
        lambda row: row.get("sample") == "H2009"
        and row.get("chrom") == "chr18"
        and row.get("active_positions") == [563687, 567920, 581932, 632412],
    )
    factorization = load_jsonl_exact(
        FACTORIZATION,
        lambda row: row.get("sample") == "H2009"
        and row.get("chrom") == "chr18"
        and row.get("active_positions") == [563687, 567920, 581932, 632412],
    )

    read_ids = [str(read_id) for read_id in assignment["all_after_peel_read_ids"]]
    labels = [str(label) for label in assignment["coarse_labels_all_after_peel"]]
    if len(read_ids) != 86 or len(set(read_ids)) != 86:
        raise AssertionError("The formal M1 after-peel read set must contain 86 unique IDs")
    if Counter(labels) != Counter({"1-1": 81, "1-2": 5}):
        raise AssertionError(f"Unexpected M1 labels: {Counter(labels)}")

    raw_distance = pd.read_csv(DISTANCE, index_col=0)
    raw_distance.index = raw_distance.index.map(str)
    raw_distance.columns = raw_distance.columns.map(str)
    if raw_distance.shape != (124, 124):
        raise AssertionError(f"Raw region distance matrix must be 124x124, found {raw_distance.shape}")
    missing_rows = sorted(set(read_ids) - set(raw_distance.index))
    missing_columns = sorted(set(read_ids) - set(raw_distance.columns))
    if missing_rows or missing_columns:
        raise AssertionError(f"Missing formal M1 IDs: rows={missing_rows}, columns={missing_columns}")
    distance = raw_distance.loc[read_ids, read_ids].to_numpy(dtype=float)
    distance = np.maximum(distance, distance.T)
    np.fill_diagonal(distance, 0.0)
    if distance.shape != (86, 86):
        raise AssertionError(f"Formal distance matrix must be 86x86, found {distance.shape}")
    if not np.isfinite(distance).all():
        raise AssertionError("Formal distance matrix contains non-finite values")
    if not np.allclose(distance, distance.T, atol=1e-12):
        raise AssertionError("Formal distance matrix is not symmetric")
    if not np.allclose(np.diag(distance), 0.0, atol=1e-12):
        raise AssertionError("Formal distance matrix diagonal is not zero")

    read_frame = pd.DataFrame(assignment["core_reads"]).copy()
    read_frame["read_id"] = read_frame["read_id"].astype(str)
    read_frame = read_frame.set_index("read_id").loc[read_ids].reset_index()
    if read_frame.shape[0] != 86:
        raise AssertionError("Core-read metadata join failed")
    read_frame["formal_label"] = labels
    if not (read_frame["label"].astype(str) == read_frame["formal_label"]).all():
        raise AssertionError("Core-read labels differ from coarse_labels_all_after_peel")

    return site, assignment, topology, factorization, read_frame, distance


def group_pair_values(distance: np.ndarray, labels: np.ndarray, group_a: str, group_b: str) -> np.ndarray:
    idx_a = np.flatnonzero(labels == group_a)
    idx_b = np.flatnonzero(labels == group_b)
    if group_a == group_b:
        block = distance[np.ix_(idx_a, idx_a)]
        return block[np.triu_indices(len(idx_a), k=1)]
    return distance[np.ix_(idx_a, idx_b)].ravel()


def build_group_distance_summary(distance: np.ndarray, labels: np.ndarray) -> tuple[list[dict], dict]:
    group_rows: list[dict] = []
    group_order = ["1-1", "1-2"]
    for row_group in group_order:
        for column_group in group_order:
            values = group_pair_values(distance, labels, row_group, column_group)
            group_rows.append(
                {
                    "row_group": f"{GROUP_META[row_group]['alias']} ({row_group})",
                    "column_group": f"{GROUP_META[column_group]['alias']} ({column_group})",
                    "mean_distance": round(float(values.mean()), 9),
                    "pair_count": int(values.size),
                }
            )

    within = np.concatenate(
        [group_pair_values(distance, labels, "1-1", "1-1"), group_pair_values(distance, labels, "1-2", "1-2")]
    )
    between = group_pair_values(distance, labels, "1-1", "1-2")
    stats = {
        "within_mean": float(within.mean()),
        "between_mean": float(between.mean()),
        "between_within_ratio": float(between.mean() / within.mean()),
        "within_pair_count": int(within.size),
        "between_pair_count": int(between.size),
    }
    return group_rows, stats


def build_upgma(distance: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, list[int], dict]:
    z_matrix = linkage(squareform(distance, checks=True), method="average")
    leaf_order = [int(index) for index in dendrogram(z_matrix, no_plot=True)["leaves"]]
    ordered_labels = labels[leaf_order]
    transition_count = int(np.sum(ordered_labels[1:] != ordered_labels[:-1]))
    ordered_counts = Counter(ordered_labels.tolist())
    if z_matrix.shape != (85, 4):
        raise AssertionError(f"Expected 85 UPGMA merges, found {z_matrix.shape}")
    if transition_count != 1 or ordered_counts != Counter({"1-1": 81, "1-2": 5}):
        raise AssertionError(
            f"UPGMA root split does not form one contiguous 81/5 block: transitions={transition_count}, "
            f"counts={ordered_counts}"
        )
    return z_matrix, leaf_order, {
        "root_height": float(z_matrix[-1, 2]),
        "transition_count": transition_count,
        "ordered_group_counts": dict(ordered_counts),
    }


def assign_display_aliases(read_frame: pd.DataFrame, leaf_order: list[int]) -> pd.DataFrame:
    ordered = read_frame.iloc[leaf_order].copy().reset_index(drop=True)
    counters = Counter()
    aliases = []
    for label in ordered["formal_label"].astype(str):
        counters[label] += 1
        aliases.append(f"{GROUP_META[label]['prefix']}{counters[label]:02d}")
    ordered.insert(0, "display_order", np.arange(1, len(ordered) + 1))
    ordered.insert(1, "read_alias", aliases)
    ordered["methyl_group"] = ordered["formal_label"].map(
        lambda value: f"{GROUP_META[value]['alias']} ({value})"
    )
    ordered["read_name_sha8"] = ordered["read_name"].map(
        lambda value: hashlib.sha256(str(value).encode("utf-8")).hexdigest()[:8]
    )
    return ordered


def pure_group_color(members: set[str]) -> str:
    if len(members) == 1:
        return GROUP_META[next(iter(members))]["color"]
    return "#55616F"


def draw_upgma_heatmap(
    z_matrix: np.ndarray,
    distance: np.ndarray,
    ordered_reads: pd.DataFrame,
    leaf_order: list[int],
    output_path: Path,
) -> None:
    ordered_distance = distance[np.ix_(leaf_order, leaf_order)]
    ordered_labels = ordered_reads["formal_label"].astype(str).to_numpy()
    aliases = ordered_reads["read_alias"].tolist()
    n_reads = len(ordered_reads)

    figure = plt.figure(figsize=(15.5, 10.8), facecolor="#F7F8F5")
    grid = figure.add_gridspec(1, 3, width_ratios=[3.4, 0.22, 8.6], wspace=0.035)
    tree_ax = figure.add_subplot(grid[0, 0])
    strip_ax = figure.add_subplot(grid[0, 1])
    heat_ax = figure.add_subplot(grid[0, 2])

    leaf_positions = {leaf_index: float(position) for position, leaf_index in enumerate(leaf_order)}
    node_height: dict[int, float] = {index: 0.0 for index in range(n_reads)}
    node_y: dict[int, float] = dict(leaf_positions)
    node_members: dict[int, set[str]] = {
        index: {str(ordered_reads.iloc[int(leaf_positions[index])]["formal_label"])} for index in range(n_reads)
    }

    for merge_index, merge in enumerate(z_matrix):
        left, right = int(merge[0]), int(merge[1])
        height = float(merge[2])
        y_left, y_right = node_y[left], node_y[right]
        members = node_members[left] | node_members[right]
        tree_ax.plot(
            [node_height[left], height],
            [y_left, y_left],
            color=pure_group_color(node_members[left]),
            linewidth=1.05,
            solid_capstyle="round",
        )
        tree_ax.plot(
            [node_height[right], height],
            [y_right, y_right],
            color=pure_group_color(node_members[right]),
            linewidth=1.05,
            solid_capstyle="round",
        )
        tree_ax.plot(
            [height, height],
            [y_left, y_right],
            color=pure_group_color(members),
            linewidth=1.2 if len(members) > 1 else 1.05,
            solid_capstyle="round",
        )
        node_id = n_reads + merge_index
        node_height[node_id] = height
        node_y[node_id] = (y_left + y_right) / 2.0
        node_members[node_id] = members

    max_height = float(z_matrix[-1, 2])
    tree_ax.set_xlim(max_height * 1.06, 0.0)
    tree_ax.set_ylim(n_reads - 0.5, -0.5)
    tree_ax.set_yticks([])
    tree_ax.set_xlabel("Average-linkage merge height", fontsize=10, color="#3D4854")
    tree_ax.set_title("UPGMA dendrogram", loc="left", fontsize=13, weight="bold", pad=12)
    tree_ax.grid(axis="x", color="#D7DDD9", linewidth=0.65, alpha=0.9)
    tree_ax.spines[["top", "right", "left"]].set_visible(False)
    tree_ax.spines["bottom"].set_color("#AAB4B0")
    tree_ax.tick_params(axis="x", labelsize=8, colors="#55616F")

    strip_values = np.array([0 if label == "1-1" else 1 for label in ordered_labels], dtype=float)[:, None]
    strip_ax.imshow(
        strip_values,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap([GROUP_META["1-1"]["color"], GROUP_META["1-2"]["color"]]),
        vmin=0,
        vmax=1,
    )
    strip_ax.set_ylim(n_reads - 0.5, -0.5)
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])
    strip_ax.set_title("MG", fontsize=10, weight="bold", pad=14)
    for spine in strip_ax.spines.values():
        spine.set_visible(False)

    distance_cmap = LinearSegmentedColormap.from_list(
        "ism_distance", ["#FBFBF7", "#C9DFDA", "#6CA6A1", "#2A6A73", "#16364D"]
    )
    image = heat_ax.imshow(
        ordered_distance,
        aspect="equal",
        interpolation="nearest",
        cmap=distance_cmap,
        vmin=0.0,
        vmax=float(distance.max()),
    )
    heat_ax.set_title("86 × 86 Bernoulli methylation distance", loc="left", fontsize=13, weight="bold", pad=12)
    heat_ax.set_xlabel("ALT reads in the same UPGMA leaf order", fontsize=10, color="#3D4854")
    heat_ax.set_ylabel("ALT reads", fontsize=10, color="#3D4854")

    transition_indices = np.flatnonzero(ordered_labels[1:] != ordered_labels[:-1])
    for transition_index in transition_indices:
        boundary = float(transition_index) + 0.5
        heat_ax.axhline(boundary, color="#F7F8F5", linewidth=2.2)
        heat_ax.axvline(boundary, color="#F7F8F5", linewidth=2.2)

    tick_positions = sorted(
        set(
            [0, n_reads - 1]
            + list(range(0, n_reads, 10))
            + [int(index) for index in transition_indices]
            + [int(index + 1) for index in transition_indices]
        )
    )
    heat_ax.set_xticks(tick_positions)
    heat_ax.set_xticklabels([aliases[index] for index in tick_positions], rotation=90, fontsize=7)
    heat_ax.set_yticks([])
    heat_ax.tick_params(length=0, colors="#55616F")
    for spine in heat_ax.spines.values():
        spine.set_color("#AAB4B0")
        spine.set_linewidth(0.7)

    colorbar = figure.colorbar(image, ax=heat_ax, fraction=0.035, pad=0.035, shrink=0.78)
    colorbar.set_label("Bernoulli methylation distance", fontsize=9)
    colorbar.ax.tick_params(labelsize=8)
    colorbar.outline.set_edgecolor("#AAB4B0")

    first_label = ordered_labels[0]
    first_count = int(np.sum(ordered_labels == first_label))
    second_label = ordered_labels[-1]
    first_alias = GROUP_META[first_label]["alias"]
    second_alias = GROUP_META[second_label]["alias"]
    figure.suptitle(
        "H2009 chr18:563687 A>G · terminal-ALT read methylation structure",
        x=0.055,
        y=0.982,
        ha="left",
        fontsize=19,
        weight="bold",
        color="#152734",
    )
    figure.text(
        0.055,
        0.948,
        f"All 86 after-peel ALT reads · {first_alias} n={first_count} then "
        f"{second_alias} n={n_reads - first_count} · same leaf order on tree and both heatmap axes",
        ha="left",
        fontsize=11,
        color="#55616F",
    )
    figure.text(
        0.055,
        0.02,
        "Average linkage visualizes read-level methylation distance; branch height is not time, "
        "and this dendrogram is not the mutation-ancestry tree.",
        ha="left",
        fontsize=9.5,
        color="#55616F",
    )
    figure.savefig(output_path, dpi=170, bbox_inches="tight", facecolor=figure.get_facecolor())
    plt.close(figure)


def select_representatives(
    distance: np.ndarray, ordered_reads: pd.DataFrame, leaf_order: list[int], max_large_group: int = 5
) -> list[int]:
    selected_original_indices: set[int] = set()
    labels_original = np.empty(len(ordered_reads), dtype=object)
    for display_position, original_index in enumerate(leaf_order):
        labels_original[original_index] = ordered_reads.iloc[display_position]["formal_label"]

    for label in ["1-1", "1-2"]:
        indices = np.flatnonzero(labels_original == label)
        submatrix = distance[np.ix_(indices, indices)]
        medoid_local = int(np.argmin(submatrix.mean(axis=1)))
        medoid_original = int(indices[medoid_local])
        ranked = sorted(
            (int(index) for index in indices),
            key=lambda index: (float(distance[medoid_original, index]), index),
        )
        keep_n = min(max_large_group, len(ranked)) if label == "1-1" else len(ranked)
        selected_original_indices.update(ranked[:keep_n])

    selected_in_leaf_order = [index for index in leaf_order if index in selected_original_indices]
    selected_labels = [str(labels_original[index]) for index in selected_in_leaf_order]
    if Counter(selected_labels) != Counter({"1-1": 5, "1-2": 5}):
        raise AssertionError(f"Unexpected representative set: {Counter(selected_labels)}")
    return selected_in_leaf_order


def representative_distance_rows(
    distance: np.ndarray,
    ordered_reads: pd.DataFrame,
    leaf_order: list[int],
    selected_indices: list[int],
) -> list[dict]:
    alias_by_original = {
        original_index: ordered_reads.iloc[display_position]["read_alias"]
        for display_position, original_index in enumerate(leaf_order)
    }
    label_by_original = {
        original_index: ordered_reads.iloc[display_position]["formal_label"]
        for display_position, original_index in enumerate(leaf_order)
    }
    rows: list[dict] = []
    for row_index in selected_indices:
        for column_index in selected_indices:
            row_label = str(label_by_original[row_index])
            column_label = str(label_by_original[column_index])
            rows.append(
                {
                    "row_read": str(alias_by_original[row_index]),
                    "column_read": str(alias_by_original[column_index]),
                    "row_group": f"{GROUP_META[row_label]['alias']} ({row_label})",
                    "column_group": f"{GROUP_META[column_label]['alias']} ({column_label})",
                    "distance": round(float(distance[row_index, column_index]), 6),
                }
            )
    if len(rows) != 100:
        raise AssertionError(f"Representative heatmap must have 100 cells, found {len(rows)}")
    return rows


def stage_dataset(
    connection: sqlite3.Connection,
    dataset_id: str,
    rows: list[dict],
    description: str,
    upstream_source_id: str,
) -> tuple[list[dict], dict]:
    if not rows:
        raise AssertionError(f"Dataset {dataset_id} is empty")
    columns = list(rows[0].keys())
    frame = pd.DataFrame(rows, columns=columns)
    frame.insert(0, "_row_order", np.arange(frame.shape[0], dtype=int))
    table_name = f"artifact_{dataset_id}"
    frame.to_sql(table_name, connection, if_exists="replace", index=False)
    quoted_columns = ", ".join(f'"{column}"' for column in columns)
    sql = f'SELECT {quoted_columns} FROM "{table_name}" ORDER BY "_row_order" ASC'
    reviewed = pd.read_sql_query(sql, connection)
    reviewed_rows = reviewed.to_dict(orient="records")
    return reviewed_rows, {
        "id": f"{dataset_id}_sql",
        "label": f"{description} via executed SQLite staging query",
        "path": "artifact_staging.sqlite",
        "query": {
            "description": (
                "Actual query executed against the persisted artifact staging database; "
                "returned rows are the exact artifact snapshot dataset."
            ),
            "engine": f"SQLite {sqlite3.sqlite_version}",
            "language": "SQL",
            "sql": sql,
            "tables_used": [table_name],
            "filters": ["dataset=H2009", "chrom=chr18", "focal_pos=563687", "after-peel ALT reads only"],
            "metric_definitions": [description],
        },
        "upstreamSourceId": upstream_source_id,
    }


def image_html_block(image_path: Path) -> str:
    encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
    return (
        '<figure style="margin:0;font-family:Inter,system-ui,sans-serif;color:#152734">'
        '<div style="display:flex;gap:8px;flex-wrap:wrap;margin:0 0 12px">'
        '<span style="padding:5px 9px;border-radius:999px;background:#E7F3F1;color:#12655F;'
        'font-size:12px;font-weight:700">DEMO · SINGLE LOCUS</span>'
        '<span style="padding:5px 9px;border-radius:999px;background:#FFF1DE;color:#8A5416;'
        'font-size:12px;font-weight:700">FULL 86 ALT READS</span>'
        "</div>"
        f'<img alt="H2009 chr18 563687 terminal ALT read UPGMA dendrogram aligned to a full '
        f'86 by 86 Bernoulli methylation distance heatmap" src="data:image/png;base64,{encoded}" '
        'style="display:block;width:100%;height:auto;box-sizing:border-box;'
        'border:1px solid #D7DDD9;border-radius:12px;'
        'background:#F7F8F5">'
        '<figcaption style="margin-top:12px;font-size:13px;line-height:1.65;color:#55616F">'
        "左側為 86 個 after-peel focal-ALT reads 的 average-linkage UPGMA；中間色帶為 formal M1 "
        "methyl group；右側為相同 leaf order 的完整 read×read Bernoulli distance。樹與矩陣使用同一排序，"
        "分枝可左右旋轉但不改變距離或群組。UPGMA 是表觀距離視覺化，不是 mutation ancestry。"
        "</figcaption></figure>"
    )


def source_inventory_rows(source_hashes: dict[str, str]) -> list[dict]:
    return [
        {
            "source": "all_ssnv_site_results.tsv.gz",
            "role": "87→86、K=2、81/5、ARI 與 strict-confirm status",
            "sha256": source_hashes["site_results"],
        },
        {
            "source": "stable_multigroup_read_assignments.jsonl.gz",
            "role": "精確 86 read IDs、labels 與 root split trace",
            "sha256": source_hashes["assignments"],
        },
        {
            "source": "distance/BERNOULLI/matrix.csv",
            "role": "124×124 primary distance；依 86 IDs exact subset",
            "sha256": source_hashes["distance"],
        },
        {
            "source": "H2009.topology.jsonl",
            "role": "representative VAF-best mutation topology",
            "sha256": source_hashes["topology"],
        },
        {
            "source": "H2009.candidate_factorization.jsonl",
            "role": "6/6 global-best terminal-edge census",
            "sha256": source_hashes["factorization"],
        },
    ]


def main() -> None:
    for directory in [FIGURES, DATA, RESULTS]:
        directory.mkdir(parents=True, exist_ok=True)

    site, assignment, topology, factorization, read_frame, distance = load_inputs()
    labels = read_frame["formal_label"].astype(str).to_numpy()
    group_rows, group_stats = build_group_distance_summary(distance, labels)
    z_matrix, leaf_order, upgma_stats = build_upgma(distance, labels)
    ordered_reads = assign_display_aliases(read_frame, leaf_order)

    first_split = assignment["coarse_split_trace"][0]
    if not first_split["passed"] or first_split["child_sizes"] != [5, 81]:
        raise AssertionError(f"Unexpected first split: {first_split}")
    if not math.isclose(
        group_stats["between_within_ratio"],
        float(first_split["observed_between_within"]),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise AssertionError(
            f"Recomputed ratio {group_stats['between_within_ratio']} != trace "
            f"{first_split['observed_between_within']}"
        )
    if not math.isclose(upgma_stats["root_height"], group_stats["between_mean"], abs_tol=1e-12):
        raise AssertionError("UPGMA root height differs from between-group average distance")

    expected_site_values = {
        "n_alt_raw": 87,
        "n_alt_matrix": 87,
        "n_alt_after_peel": 86,
        "n_alt_peeled": 1,
        "coarse_ng": 2,
    }
    for field, expected in expected_site_values.items():
        if int(site[field]) != expected:
            raise AssertionError(f"{field}: expected {expected}, found {site[field]}")
    if json.loads(site["cluster_sizes"]) != {"1-1": 81, "1-2": 5}:
        raise AssertionError(f"Unexpected cluster sizes: {site['cluster_sizes']}")
    if not as_bool(site["stable_null_multigroup"]) or as_bool(site["unstable"]):
        raise AssertionError("M1 stability flags do not support stable multigroup")
    if not math.isclose(float(site["modal_assignment_ari_min"]), 1.0, abs_tol=1e-12):
        raise AssertionError("Expected modal assignment minimum ARI=1.0")
    if site["strict_confirm_status"] != "NOT_RUN":
        raise AssertionError(f"Unexpected strict confirm status: {site['strict_confirm_status']}")

    terminal_edges = [
        edge
        for edge in topology["representative_best_edges"]
        if int(edge["parent_vertex"]) == 14 and int(edge["child_vertex"]) == 15
    ]
    if len(terminal_edges) != 1 or int(terminal_edges[0]["acquired_position"]) != 563687:
        raise AssertionError("Representative best topology lacks terminal 14→15 / chr18:563687 edge")
    global_terminal = [
        edge
        for edge in factorization["global_best_edge_incidence"]
        if int(edge[0]) == 14 and int(edge[1]) == 15
    ]
    best_tree_count = int(factorization["best_tree_tie_count"])
    if global_terminal != [[14, 15, "6"]] or best_tree_count != 6:
        raise AssertionError(f"Expected terminal edge in 6/6 best trees, found {global_terminal}/{best_tree_count}")
    if factorization["global_best_coarse_class_tree_counts"] != {"Direct-only": "6"}:
        raise AssertionError("Expected six Direct-only VAF-global-best trees")

    draw_upgma_heatmap(z_matrix, distance, ordered_reads, leaf_order, FIGURE_OUT)

    leaf_export = ordered_reads[
        [
            "display_order",
            "read_alias",
            "read_id",
            "read_name_sha8",
            "methyl_group",
            "formal_label",
            "latest_hp",
            "latest_ps",
            "strand",
            "start",
            "end",
            "mapq",
        ]
    ].copy()
    leaf_export.to_csv(LEAF_ORDER_OUT, sep="\t", index=False)
    pd.DataFrame(group_rows).to_csv(GROUP_DISTANCE_OUT, sep="\t", index=False)

    representatives = select_representatives(distance, ordered_reads, leaf_order)
    representative_rows = representative_distance_rows(distance, ordered_reads, leaf_order, representatives)

    summary_metrics = [
        {
            "alt_raw": int(site["n_alt_raw"]),
            "alt_after_peel": int(site["n_alt_after_peel"]),
            "peeled": int(site["n_alt_peeled"]),
            "methyl_group_count": int(site["coarse_ng"]),
            "group_a_reads": 81,
            "group_b_reads": 5,
            "between_within_ratio": round(group_stats["between_within_ratio"], 6),
            "null95_threshold": round(float(first_split["null_threshold"]), 6),
            "empirical_rank_p": round(float(first_split["empirical_p"]), 6),
            "assignment_ari_min": round(float(site["modal_assignment_ari_min"]), 6),
            "upgma_root_height": round(upgma_stats["root_height"], 6),
            "vaf_best_tree_count": best_tree_count,
            "terminal_edge_tree_count": int(global_terminal[0][2]),
        }
    ]

    source_hashes = {
        "site_results": EXPECTED_SITE_RESULTS_SHA,
        "assignments": EXPECTED_ASSIGNMENT_SHA,
        "reads": EXPECTED_PRIMARY_SHA[READS],
        "distance": EXPECTED_PRIMARY_SHA[DISTANCE],
        "methylation": EXPECTED_PRIMARY_SHA[METHYLATION],
        "topology": sha256_file(TOPOLOGY),
        "factorization": sha256_file(FACTORIZATION),
        "figure": sha256_file(FIGURE_OUT),
    }
    inventory_rows = source_inventory_rows(source_hashes)

    if STAGING_DB.exists():
        STAGING_DB.unlink()
    connection = sqlite3.connect(STAGING_DB)
    snapshot_datasets: dict[str, list[dict]] = {}
    sql_sources: list[dict] = []
    dataset_specs = [
        (
            "summary_metrics",
            summary_metrics,
            "H2009 terminal-ALT headline counts, M1 separation, stability, and VAF terminal-edge support",
            "m1_assignment",
        ),
        (
            "representative_distance_cells",
            representative_rows,
            "Deterministic 10/86 read display subset; distances retain the full 86-read source values",
            "region_distance",
        ),
        (
            "group_distance_summary",
            group_rows,
            "Full-86 read within- and between-methyl-group mean Bernoulli distances",
            "region_distance",
        ),
        (
            "source_inventory",
            inventory_rows,
            "Exact source inventory and SHA256 provenance",
            "m1_assignment",
        ),
    ]
    for dataset_id, rows, description, upstream in dataset_specs:
        reviewed_rows, source = stage_dataset(connection, dataset_id, rows, description, upstream)
        snapshot_datasets[dataset_id] = reviewed_rows
        sql_sources.append(source)
    connection.commit()
    connection.close()
    staging_sha = sha256_file(STAGING_DB)
    for source in sql_sources:
        source["sha256"] = staging_sha

    raw_sources = [
        {
            "id": "m1_site_results",
            "label": "v10 source-locked all-sSNV site results",
            "path": "all_ssnv_site_results.tsv.gz",
            "sha256": source_hashes["site_results"],
            "query": {
                "description": "Exact row extraction for H2009 chr18:563687 A>G.",
                "engine": "Python 3 / csv streaming",
                "language": "Python",
                "tables_used": ["all_ssnv_site_results.tsv.gz"],
                "filters": ["dataset=H2009", "chrom=chr18", "pos=563687"],
                "metric_definitions": [
                    "M1 groups are stable focal-ALT methylation-distance groups, not cellular clone labels."
                ],
            },
        },
        {
            "id": "m1_assignment",
            "label": "v10 stable focal-ALT read assignment",
            "path": "all_ssnv_stable_multigroup_read_assignments.jsonl.gz",
            "sha256": source_hashes["assignments"],
            "query": {
                "description": "Exact JSONL row extraction and 86-ID join.",
                "engine": "Python 3 streaming JSONL",
                "language": "Python",
                "tables_used": ["all_ssnv_stable_multigroup_read_assignments.jsonl.gz"],
                "filters": ["dataset=H2009", "chrom=chr18", "pos=563687"],
            },
        },
        {
            "id": "region_distance",
            "label": "H2009 chr18:563687 primary read-distance matrix",
            "path": "chr18_558687_568687/distance/BERNOULLI/matrix.csv",
            "sha256": source_hashes["distance"],
            "query": {
                "description": "Subset the 124x124 primary region matrix to the exact 86 after-peel focal-ALT IDs.",
                "engine": "Python 3 / pandas / scipy",
                "language": "Python",
                "tables_used": ["matrix.csv", "reads.tsv"],
                "filters": ["86 assignment read IDs", "D=max(D,D.T)", "diagonal=0", "average linkage"],
                "metric_definitions": [
                    "Bernoulli methylation distance is the M1 read-pair distance.",
                    "UPGMA branch height is a merge distance, not elapsed evolutionary time.",
                ],
            },
        },
        {
            "id": "topology_factorization",
            "label": "Exact-PS VAF-global-best topology factorization",
            "path": "H2009.candidate_factorization.jsonl",
            "sha256": source_hashes["factorization"],
            "query": {
                "description": "Exact candidate extraction for active positions 563687,567920,581932,632412.",
                "engine": "Python 3 streaming JSONL",
                "language": "Python",
                "tables_used": ["H2009.topology.jsonl", "H2009.candidate_factorization.jsonl"],
                "filters": ["sample=H2009", "chrom=chr18", "phase_set=307022", "hp_family=1"],
                "metric_definitions": [
                    "6/6 VAF-global-best trees share edge 14→15, whose acquired position is 563687.",
                    "The six trees tie in upstream mutation order; this does not make methyl UPGMA a mutation tree.",
                ],
            },
        },
    ]
    sources = raw_sources + sql_sources

    title = "H2009 chr18:563687 ALT-read 甲基距離與 UPGMA"
    figure_body = image_html_block(FIGURE_OUT)
    now = datetime.now(timezone.utc).isoformat()
    manifest = {
        "version": 1,
        "surface": "report",
        "title": title,
        "description": (
            "單一 H2009 terminal-ALT 位點的完整 86-read 甲基距離 UPGMA、同序 86×86 熱圖、"
            "10-read native detail 與 exact-PS VAF topology context。"
        ),
        "generatedAt": now,
        "sources": sources,
        "cards": [
            {
                "id": "read_scope_card",
                "description": "Formal M1 focal-ALT read scope.",
                "dataset": "summary_metrics",
                "sourceId": "summary_metrics_sql",
                "metrics": [
                    {"label": "After-peel ALT reads", "field": "alt_after_peel", "format": "number"},
                    {"label": "Raw ALT", "field": "alt_raw", "format": "number"},
                    {"label": "Peeled", "field": "peeled", "format": "number"},
                ],
            },
            {
                "id": "group_card",
                "description": "Formal stable methyl groups; K is not clone count.",
                "dataset": "summary_metrics",
                "sourceId": "summary_metrics_sql",
                "metrics": [
                    {"label": "Methyl groups K", "field": "methyl_group_count", "format": "number"},
                    {"label": "MG-A", "field": "group_a_reads", "format": "number"},
                    {"label": "MG-B", "field": "group_b_reads", "format": "number"},
                ],
            },
            {
                "id": "separation_card",
                "description": "Full-86 M1 root split versus column-shuffle null95.",
                "dataset": "summary_metrics",
                "sourceId": "summary_metrics_sql",
                "metrics": [
                    {"label": "Between / within", "field": "between_within_ratio", "format": "number"},
                    {"label": "Null95", "field": "null95_threshold", "format": "number"},
                    {"label": "Rank p", "field": "empirical_rank_p", "format": "number"},
                ],
            },
            {
                "id": "stability_card",
                "description": "Seed stability and independent VAF-topology terminal-edge census.",
                "dataset": "summary_metrics",
                "sourceId": "summary_metrics_sql",
                "metrics": [
                    {"label": "Min assignment ARI", "field": "assignment_ari_min", "format": "number"},
                    {"label": "Terminal edge trees", "field": "terminal_edge_tree_count", "format": "number"},
                    {"label": "VAF-best trees", "field": "vaf_best_tree_count", "format": "number"},
                ],
            },
        ],
        "charts": [
            {
                "id": "representative_distance_heatmap",
                "title": "10/86 representative ALT-read methylation-distance matrix",
                "subtitle": "5 deterministic MG-A medoid-nearest reads＋全部 5 MG-B reads；統計與主圖使用完整 86 reads。",
                "intent": "relationship",
                "question": "完整 86-read block structure 在可辨識的 read 子集中是否仍可直接觀察？",
                "rationale": "10×10 native heatmap保留 individual read 距離，並可在 desktop/mobile 一次看完整。",
                "type": "heatmap",
                "dataset": "representative_distance_cells",
                "sourceId": "representative_distance_cells_sql",
                "encodings": {
                    "x": {"field": "row_read", "type": "nominal", "label": "ALT read alias"},
                    "y": {"field": "distance", "type": "quantitative", "label": "Bernoulli distance"},
                    "color": {"field": "column_read", "type": "nominal", "label": "ALT read alias"},
                },
                "valueFormat": "number",
                "palette": {"kind": "sequential"},
            }
        ],
        "tables": [
            {
                "id": "group_distance_table",
                "title": "完整 86 reads 的群內／群間距離",
                "subtitle": "同群只計 unordered off-diagonal pairs；跨群計 81×5 pairs。",
                "dataset": "group_distance_summary",
                "sourceId": "group_distance_summary_sql",
                "defaultSort": {"field": "mean_distance", "direction": "desc"},
                "density": "compact",
                "columns": [
                    {"field": "row_group", "label": "Row group", "type": "text"},
                    {"field": "column_group", "label": "Column group", "type": "text"},
                    {"field": "mean_distance", "label": "Mean distance", "type": "number"},
                    {"field": "pair_count", "label": "Pair count", "type": "number"},
                ],
            },
            {
                "id": "source_inventory_table",
                "title": "來源與 SHA256",
                "subtitle": "主矩陣、M1 assignment 與 exact-PS factorization 的可追溯來源。",
                "dataset": "source_inventory",
                "sourceId": "source_inventory_sql",
                "defaultSort": {"field": "source", "direction": "asc"},
                "density": "compact",
                "columns": [
                    {"field": "source", "label": "Source", "type": "text"},
                    {"field": "role", "label": "Role", "type": "text"},
                    {"field": "sha256", "label": "SHA256", "type": "text"},
                ],
            },
        ],
        "blocks": [
            {"id": "title", "type": "markdown", "body": f"# {title}"},
            {
                "id": "scope",
                "type": "markdown",
                "sourceId": "m1_site_results",
                "body": (
                    "## 一句結論\n\n"
                    "> **DEMO / PARTIAL SCOPE — 單一位點，服務 G3。** "
                    "H2009 `chr18:563687 A>G` 的 86 個 after-peel ALT reads，在同一 HP1／PS=307022 內形成 "
                    "**81＋5** 的穩定甲基距離二群；完整 UPGMA 與同序熱圖可直接看到兩個 block。"
                    "這支持 read-level residual epigenetic substructure，但不等於兩個 cellular clones。"
                ),
            },
            {
                "id": "headline_metrics",
                "type": "metric-strip",
                "cardIds": ["read_scope_card", "group_card", "separation_card", "stability_card"],
            },
            {
                "id": "full_figure_intro",
                "type": "markdown",
                "sourceId": "region_distance",
                "body": (
                    "## 完整 86-read UPGMA＋同序距離熱圖\n\n"
                    "左樹與右矩陣由同一個 canonical average-linkage leaf order 產生；熱圖的 rows 與 columns "
                    "都依該順序排列，因此對角線兩側的 block 可與樹枝直接對齊。"
                ),
            },
            {
                "id": "full_figure",
                "type": "html",
                "sourceId": "region_distance",
                "body": figure_body,
            },
            {
                "id": "reading_guide",
                "type": "markdown",
                "sourceId": "m1_assignment",
                "body": (
                    "## 如何讀這張圖\n\n"
                    "1. 深色代表 read 間甲基 pattern 距離較大；對角線為 0。  \n"
                    "2. MG-A（81 reads）與 MG-B（5 reads）各自形成低距離方塊，跨群區塊較深。  \n"
                    "3. 完整 86 reads 的群內平均距離為 **0.174282**，群間為 **0.504450**，比值 "
                    "**2.894439**；高於 40 次 column-shuffle 的 null95 **1.381452**。  \n"
                    "4. empirical rank p=`1/41=0.02439` 是有限 null 排名描述，不重標為傳統漸近 p-value。"
                ),
            },
            {
                "id": "native_detail_intro",
                "type": "markdown",
                "sourceId": "representative_distance_cells_sql",
                "body": (
                    "## 可辨讀的 10-read native detail\n\n"
                    "Portable artifact 每個 dataset 上限 2,000 rows；完整 tidy 86²=7,396 cells 已由上方主圖呈現。"
                    "下圖固定選 MG-A medoid-nearest 5 reads，並保留全部 5 個 MG-B reads；只縮小顯示，不重算 "
                    "K、ratio 或 UPGMA 結論。"
                ),
            },
            {"id": "native_heatmap", "type": "chart", "chartId": "representative_distance_heatmap"},
            {"id": "group_distance", "type": "table", "tableId": "group_distance_table"},
            {
                "id": "topology_context",
                "type": "markdown",
                "sourceId": "topology_factorization",
                "body": (
                    "## 與突變拓撲的關係\n\n"
                    "獨立的 exact-PS／HP1 候選樹以每個位點 VAF 評分後，有 **6 個並列最可能樹**；6/6 都是 "
                    "Direct-only，且都包含 terminal edge `14→15`，此 edge 新增 `chr18:563687`。"
                    "因此可以把甲基二群標在「terminal 563687 ALT-bearing molecules」上；但上游三個突變的次序仍有 "
                    "6 種並列，甲基 UPGMA 也不能用來替它們指定父子順序。"
                ),
            },
            {
                "id": "limitations",
                "type": "markdown",
                "sourceId": "m1_site_results",
                "body": (
                    "## Claim ceiling\n\n"
                    "- 可說：同一 terminal-ALT、同一 HP1／PS 的 reads 內，有 seed-stable、null-calibrated 的甲基多群。  \n"
                    "- 不可直接說：MG-A／MG-B 就是兩個 cellular clones、甲基造成突變、突變造成甲基，或 UPGMA 是真實演化樹。  \n"
                    "- 此位點的 strict confirm status=`NOT_RUN`；要升級為 clone/subclone，仍需 linked independent marker、"
                    "matched-normal、CN/purity/multiplicity/CCF 與 single-cell/colony/spatial 正交驗證。"
                ),
            },
            {
                "id": "source_intro",
                "type": "markdown",
                "sourceId": "source_inventory_sql",
                "body": "## 可追溯來源\n\n以下 SHA256 與 exact-row filters 由建置腳本實際核對。",
            },
            {"id": "source_inventory", "type": "table", "tableId": "source_inventory_table"},
        ],
    }
    artifact = {
        "surface": "report",
        "manifest": manifest,
        "snapshot": {
            "version": 1,
            "generatedAt": now,
            "status": "ready",
            "datasets": snapshot_datasets,
        },
        "sources": sources,
    }
    ARTIFACT_OUT.write_text(json.dumps(artifact, ensure_ascii=False, indent=2), encoding="utf-8")

    validation_summary = {
        "schema_name": "intersubmod.h2009_terminal_alt_methyl_upgma_visual",
        "schema_version": "1.0.0",
        "status": "PASS",
        "generated_at": now,
        "task_type": "F_DEMO_ILLUSTRATION",
        "goal": "G3",
        "scope": {
            "sample": "H2009",
            "locus": "chr18:563687 A>G",
            "hp": "1-1",
            "phase_set": 307022,
            "alt_raw": 87,
            "alt_after_peel": 86,
            "peeled_read_id": "69",
            "methyl_groups": {"MG-A (1-1)": 81, "MG-B (1-2)": 5},
        },
        "upgma": {
            **upgma_stats,
            "linkage_method": "average",
            "linkage_rows": int(z_matrix.shape[0]),
            "same_leaf_order_for_tree_rows_columns": True,
        },
        "m1_root_split": {
            **group_stats,
            "formal_trace_ratio": float(first_split["observed_between_within"]),
            "null95_threshold": float(first_split["null_threshold"]),
            "empirical_rank_p": float(first_split["empirical_p"]),
            "valid_nulls": int(first_split["n_valid_null"]),
            "passed": bool(first_split["passed"]),
            "assignment_ari_min": float(site["modal_assignment_ari_min"]),
        },
        "mutation_topology_context": {
            "active_positions": topology["active_positions"],
            "best_tree_tie_count": best_tree_count,
            "global_best_coarse_class_tree_counts": factorization["global_best_coarse_class_tree_counts"],
            "terminal_edge": [14, 15],
            "terminal_position": 563687,
            "terminal_edge_tree_count": int(global_terminal[0][2]),
            "terminal_edge_in_all_best_trees": True,
        },
        "checks": {
            "primary_sha_match": True,
            "raw_region_matrix_shape": [124, 124],
            "formal_matrix_shape": [86, 86],
            "matrix_symmetric": True,
            "matrix_diagonal_zero": True,
            "formal_ratio_exact_match": True,
            "group_transition_count": upgma_stats["transition_count"],
            "representative_heatmap_cells": len(representative_rows),
            "strict_confirm_status": site["strict_confirm_status"],
        },
        "outputs": {
            "artifact": str(ARTIFACT_OUT.relative_to(REPO)),
            "figure": str(FIGURE_OUT.relative_to(REPO)),
            "leaf_order": str(LEAF_ORDER_OUT.relative_to(REPO)),
            "group_distance": str(GROUP_DISTANCE_OUT.relative_to(REPO)),
            "staging_db": str(STAGING_DB.relative_to(REPO)),
        },
        "sha256": {**source_hashes, "staging_db": staging_sha, "artifact": sha256_file(ARTIFACT_OUT)},
        "claim_ceiling": (
            "Supports stable read-level methylation substructure within terminal-ALT molecules; "
            "does not establish cellular clones, causality, or a unique mutation-ancestry tree."
        ),
    }
    SUMMARY_OUT.write_text(json.dumps(validation_summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(
        json.dumps(
            {
                "status": "PASS",
                "artifact": str(ARTIFACT_OUT),
                "figure": str(FIGURE_OUT),
                "summary": str(SUMMARY_OUT),
                "formal_reads": 86,
                "groups": {"1-1": 81, "1-2": 5},
                "between_within_ratio": group_stats["between_within_ratio"],
                "upgma_root_height": upgma_stats["root_height"],
                "terminal_edge_best_trees": f"{global_terminal[0][2]}/{best_tree_count}",
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
