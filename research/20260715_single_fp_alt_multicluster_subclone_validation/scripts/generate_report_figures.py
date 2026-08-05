#!/usr/bin/env python3
"""Generate aggregate and actual-case figures from the frozen report dataset."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap


TOPIC_ROOT = Path(__file__).resolve().parents[1]
COLORS = {
    "ink": "#1F2523",
    "muted": "#6B6F6C",
    "accent": "#D97757",
    "teal": "#007C83",
    "gold": "#C69214",
    "purple": "#715A8C",
    "red": "#B94A48",
    "green": "#2D7A57",
    "light": "#E8E6DF",
    "missing": "#ECEBE6",
}
GROUP_COLORS = [COLORS["teal"], COLORS["accent"], COLORS["purple"], COLORS["gold"], COLORS["green"]]


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Noto Sans CJK TC", "DejaVu Sans", "sans-serif"],
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "axes.edgecolor": "#B9B7B0",
            "axes.linewidth": 0.8,
            "xtick.color": COLORS["muted"],
            "ytick.color": COLORS["muted"],
            "text.color": COLORS["ink"],
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=190, bbox_inches="tight")
    plt.close(fig)


def pct(value: float) -> str:
    return f"{100 * value:.2f}%"


def add_source(fig: plt.Figure, text: str) -> None:
    fig.text(0.01, 0.005, text, fontsize=7, color=COLORS["muted"], ha="left", va="bottom")


def evidence_chain(report: dict[str, Any], output: Path) -> None:
    fp = report["fp_primary"]
    strict = report["strict_sensitivity"]
    labels = [
        "Latest LongPhase-S PASS truth-FP",
        "Evaluable focal ALT readsets",
        "Column-null stable multigroup",
        "Residual after measured axes",
        "HP-tag-covered high threshold",
        "Strict follow-up\n10 sites / 9 readsets / 6 components",
    ]
    values = [fp["n_sites"], fp["n_evaluable"], fp["n_stable"], fp["n_residual"], fp["n_hp_tag_covered_high_threshold"], strict["n_strict_epigenetic_followup_candidates"]]
    colors = [COLORS["muted"], COLORS["teal"], COLORS["gold"], COLORS["accent"], COLORS["purple"], COLORS["green"]]
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=colors, height=0.62)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(5, 12000)
    ax.set_title("Conservative evidence chain: heterogeneity signal contracts from 7,745 sites to 10 follow-up candidates", loc="left", fontweight="bold")
    ax.set_xlabel("Number of sites (log scale)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for bar, value in zip(bars, values):
        ax.text(value * 1.08, bar.get_y() + bar.get_height() / 2, f"{value:,}  ({value / values[0]:.2%} of all)", va="center", fontsize=8.5)
    ax.text(0.99, -0.26, "The chain is evidentiary, not a PPV/FDR calibration and not a clone count.", transform=ax.transAxes, ha="right", color=COLORS["red"], fontsize=8.5)
    add_source(fig, "Source: report_dataset_v1/report_dataset.json; latest_full_v3_frozen + strict_survival_v1")
    fig.subplots_adjust(left=0.34, right=0.94, bottom=0.22, top=0.85)
    save(fig, output)


def matched_specificity(report: dict[str, Any], output: Path) -> None:
    endpoints = [
        ("stable_null_multigroup", "Stable multigroup"),
        ("residual_unexplained_multigroup", "Residual multigroup"),
        ("phase_anchored_robust_epigenetic_candidate", "HP-tag-covered high threshold"),
    ]
    data = report["matched_tp_specificity"]["all_pairs"]
    fp = [data[key]["fp_fraction"] for key, _ in endpoints]
    tp = [data[key]["tp_fraction"] for key, _ in endpoints]
    p_values = [data[key]["exact_mcnemar_p"] for key, _ in endpoints]
    x = np.arange(len(endpoints))
    width = 0.34
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    left = ax.bar(x - width / 2, fp, width, label="Truth FP", color=COLORS["accent"])
    right = ax.bar(x + width / 2, tp, width, label="Matched truth TP", color=COLORS["teal"])
    ax.set_xticks(x, [label for _, label in endpoints])
    ax.set_ylabel("Fraction of all matched sites")
    ax.set_ylim(0, max(fp + tp) * 1.35)
    ax.set_title("The high-threshold endpoint has no truth-FP specificity", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    for bars in (left, right):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.0015, pct(bar.get_height()), ha="center", va="bottom", fontsize=8)
    for index, p_value in enumerate(p_values):
        label = f"paired p={p_value:.3g}"
        ax.text(index, max(fp[index], tp[index]) * 1.22, label, ha="center", fontsize=8, color=COLORS["muted"])
    add_source(fig, "Source: fp_matched_tp_comparison_v1; n=7,745 matched pairs")
    fig.subplots_adjust(bottom=0.21, top=0.86)
    save(fig, output)


def per_sample_risk_difference(report: dict[str, Any], output: Path) -> None:
    rows = report["per_sample"]
    samples = [row["sample"].replace("HCC1395_DORADO", "HCC1395 Dorado") for row in rows]
    stable = [100 * row["matched_stable_risk_difference"] for row in rows]
    high = [100 * row["matched_high_threshold_risk_difference"] for row in rows]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    ax.axvline(0, color=COLORS["ink"], linewidth=0.9)
    ax.scatter(stable, y - 0.13, s=52, color=COLORS["accent"], label="Stable multigroup", marker="o")
    ax.scatter(high, y + 0.13, s=55, color=COLORS["teal"], label="HP-tag-covered high threshold", marker="s")
    for index, (left, right) in enumerate(zip(stable, high)):
        ax.plot([left, right], [index - 0.13, index + 0.13], color=COLORS["light"], linewidth=1)
    ax.set_yticks(y, samples)
    ax.invert_yaxis()
    ax.set_xlabel("Paired risk difference: truth FP minus matched TP (percentage points)")
    ax.set_title("Cross-dataset directions are inconsistent; the pooled stable effect is HCC1395-sensitive", loc="left", fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(frameon=False, loc="lower right")
    ax.text(0.01, -0.17, "HCC1395 and HCC1395 Dorado are one biological sample. Removing both: stable RD = -0.153 pp, p = 0.780.", transform=ax.transAxes, color=COLORS["red"], fontsize=8.3)
    add_source(fig, "Source: fp_matched_tp_robustness_v1 + per_sample_metrics.tsv")
    fig.subplots_adjust(left=0.22, bottom=0.2, top=0.86)
    save(fig, output)


def strict_sensitivity(report: dict[str, Any], output: Path) -> None:
    strict = report["strict_sensitivity"]
    labels = [
        "Primary stable",
        "Column empirical-p pass",
        "Row-circular empirical-p pass",
        "High-threshold endpoint",
        "High-threshold + REF background",
        "High-threshold + row-circular",
        "All strict gates",
    ]
    values = [
        strict["n_primary_group_count_stable_candidates"],
        strict["n_column_empirical_p_gated_survivors"],
        strict["n_row_circular_empirical_p_gated_survivors"],
        strict["n_primary_e4_hp_covered_candidates"],
        strict["n_e4_background_surviving"],
        strict["n_e4_row_circular_assignment_robust_survivors"],
        strict["n_strict_epigenetic_followup_candidates"],
    ]
    colors = [COLORS["muted"], COLORS["gold"], COLORS["accent"], COLORS["purple"], COLORS["teal"], COLORS["accent"], COLORS["green"]]
    fig, ax = plt.subplots(figsize=(9.3, 5.0))
    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=colors, height=0.62)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 640)
    ax.set_xlabel("Candidate sites")
    ax.set_title("A row-structure-preserving null and allele background remove most nominal candidates", loc="left", fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for bar, value in zip(bars, values):
        ax.text(value + 8, bar.get_y() + bar.get_height() / 2, f"{value:,}", va="center", fontweight="bold")
    ax.text(0.99, -0.16, "The bars are overlapping sensitivity sets, not a calibrated discovery funnel.", transform=ax.transAxes, ha="right", color=COLORS["red"], fontsize=8.5)
    add_source(fig, "Source: strict_survival_v1; 199 null replicates, 10 seeds, minimum assignment ARI 0.8")
    fig.subplots_adjust(left=0.33, bottom=0.18, top=0.86)
    save(fig, output)


def topology_context(report: dict[str, Any], output: Path) -> None:
    topology = report["topology"]
    strict = report["strict_sensitivity"]["topology_context_counts"]
    source_sets = {
        "All FP sites": topology["all_sites"]["topology_context_counts"],
        "Stable": topology["stable_null_multigroup"]["topology_context_counts"],
        "High threshold": topology["phase_anchored_robust_epigenetic_candidate"]["topology_context_counts"],
        "Strict follow-up": strict,
    }
    categories = [
        ("Linear regional", ["PRIMARY_ALL_STORED_TREES_LINEAR"], COLORS["teal"]),
        ("Branching / recurrence", ["PRIMARY_INCLUDES_BRANCHING_OR_RECURRENCE"], COLORS["accent"]),
        ("Trivial one edge", ["PRIMARY_TRIVIAL_ONE_EDGE_ONLY"], COLORS["gold"]),
        ("Incomplete", ["PRIMARY_CANDIDATE_SET_INCOMPLETE"], COLORS["purple"]),
        ("No usable tree", ["NO_REGIONAL_TREE_FOR_SITE", "AUXILIARY_OR_UNPHASED_ONLY", "NO_ALT_FAMILY_UNIT_MATCH"], COLORS["light"]),
    ]
    fig, ax = plt.subplots(figsize=(9.4, 4.8))
    left = np.zeros(len(source_sets))
    x = np.arange(len(source_sets))
    for label, keys, color in categories:
        values = []
        for counts in source_sets.values():
            total = sum(counts.values())
            values.append(sum(counts.get(key, 0) for key in keys) / total)
        bars = ax.bar(x, values, bottom=left, color=color, width=0.62, label=label)
        for index, (bar, value) in enumerate(zip(bars, values)):
            if value >= 0.08:
                ax.text(bar.get_x() + bar.get_width() / 2, left[index] + value / 2, f"{value:.0%}", ha="center", va="center", fontsize=7.5)
        left += np.array(values)
    ax.set_xticks(x, list(source_sets))
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction within endpoint")
    ax.set_title("The methylation endpoint appears in both linear and branching regional contexts", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    ax.text(0.99, -0.31, "High-threshold linear-context enrichment: Fisher OR 1.217, p=0.441. Focal linear topology identifiable: 0 sites.", transform=ax.transAxes, ha="right", color=COLORS["red"], fontsize=8.3)
    add_source(fig, "Source: latest_topology_context_summary.json + strict candidate join")
    fig.subplots_adjust(bottom=0.3, top=0.86)
    save(fig, output)


def hp_tags(report: dict[str, Any], output: Path) -> None:
    scopes = [
        ("All focal ALT", report["hp_tags"]["all_fp"]),
        ("Stable", report["hp_tags"]["stable"]),
        ("High threshold", report["hp_tags"]["hp_tag_covered_high_threshold"]),
        ("Strict follow-up", report["hp_tags"]["strict_followup"]),
    ]
    categories = [
        ("First somatic haplotype", "first_somatic_haplotype", COLORS["teal"]),
        ("Ambiguous HP3", "ambiguous", COLORS["gold"]),
        ("Untagged / other", "untagged_or_other", COLORS["light"]),
        ("Explicit second haplotype", "second_somatic_haplotype", COLORS["red"]),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    left = np.zeros(len(scopes))
    x = np.arange(len(scopes))
    for label, key, color in categories:
        values = [scope["family_counts"].get(key, 0) / scope["n_alt_reads"] for _, scope in scopes]
        ax.bar(x, values, bottom=left, color=color, width=0.62, label=label)
        left += np.array(values)
    ax.set_xticks(x, [label for label, _ in scopes])
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Fraction of focal ALT reads")
    ax.set_title("No focal ALT read carries an explicit LongPhase-S second-somatic-haplotype tag", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    totals = [scope["n_alt_reads"] for _, scope in scopes]
    for index, total in enumerate(totals):
        ax.text(index, 1.005, f"n={total:,}", ha="center", va="bottom", fontsize=8)
    ax.text(0.99, -0.28, "HP tag coverage is not orthogonal clone confirmation; methylation group labels are separate from HP tags.", transform=ax.transAxes, ha="right", color=COLORS["red"], fontsize=8.3)
    add_source(fig, "Source: alt_hp_counts parsed from all 7,745 latest FP site rows")
    fig.subplots_adjust(bottom=0.28, top=0.84)
    save(fig, output)


def coverage_dependence(report: dict[str, Any], output: Path) -> None:
    rows = report["fp_primary"]["coverage_bands"]
    x = np.arange(len(rows))
    stable = [row["stable_fraction"] for row in rows]
    high = [row["hp_tag_covered_high_threshold_fraction"] for row in rows]
    width = 0.34
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    a = ax.bar(x - width / 2, stable, width, color=COLORS["accent"], label="Stable multigroup")
    b = ax.bar(x + width / 2, high, width, color=COLORS["teal"], label="High threshold")
    ax.set_xticks(x, [f"{row['band']}\n(n={row['n_evaluable']:,})" for row in rows])
    ax.set_ylabel("Fraction of evaluable sites")
    ax.set_title("Detection probability rises strongly with focal ALT read depth", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    for bars in (a, b):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.004, pct(bar.get_height()), ha="center", fontsize=8)
    ax.set_ylim(0, max(stable) * 1.32)
    add_source(fig, "Source: latest_full_v3_frozen; coverage = n_alt_after_peel")
    fig.subplots_adjust(bottom=0.2, top=0.86)
    save(fig, output)


def legacy_compatibility(report: dict[str, Any], output: Path) -> None:
    agreement = report["legacy_sensitivity"]["endpoint_agreement"]
    endpoints = [
        ("stable_null_multigroup", "Stable"),
        ("residual_unexplained_multigroup", "Residual"),
        ("phase_anchored_robust_epigenetic_candidate", "High threshold"),
    ]
    latest = [agreement[key]["latest_fraction_all_common"] for key, _ in endpoints]
    legacy = [agreement[key]["legacy_canonical_fraction_all_common"] for key, _ in endpoints]
    jaccard = [agreement[key]["jaccard"] for key, _ in endpoints]
    x = np.arange(len(endpoints))
    width = 0.34
    fig, ax = plt.subplots(figsize=(8.7, 4.8))
    a = ax.bar(x - width / 2, latest, width, color=COLORS["teal"], label="Latest frozen semantics")
    b = ax.bar(x + width / 2, legacy, width, color=COLORS["muted"], label="Legacy canonical semantics")
    ax.set_xticks(x, [label for _, label in endpoints])
    ax.set_ylabel("Fraction among 3,302 common FP sites")
    ax.set_title("Prior InterSubMod outputs reproduce the same phenomenon, but not an exchangeable rate", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    for bars in (a, b):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.004, pct(bar.get_height()), ha="center", fontsize=8)
    for index, value in enumerate(jaccard):
        ax.text(index, max(latest[index], legacy[index]) * 1.18, f"endpoint Jaccard={value:.3f}", ha="center", color=COLORS["muted"], fontsize=8)
    ax.set_ylim(0, max(legacy) * 1.38)
    ax.text(0.99, -0.18, "Methylation matrices: 100% identical. Bernoulli distance matrices: 63.0% identical.", transform=ax.transAxes, ha="right", color=COLORS["red"], fontsize=8.3)
    add_source(fig, "Source: latest_canonical_compatibility_audit.json; legacy is sensitivity context only")
    fig.subplots_adjust(bottom=0.2, top=0.86)
    save(fig, output)


def load_assignments(path: Path) -> dict[tuple[str, str, int], dict[str, Any]]:
    result = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            result[(row["sample"], row["chrom"], int(row["pos"]))] = row
    return result


def load_case_matrix(case: dict[str, Any], assignment: dict[str, Any]) -> tuple[np.ndarray, list[int], list[str], dict[str, str]]:
    matrix_path = Path(case["region_dir"]) / "methylation" / "methylation.csv"
    frame = pd.read_csv(matrix_path, na_values=["NA"])
    frame["read_id"] = frame["read_id"].astype(str)
    frame = frame.set_index("read_id")
    read_ids = [str(value) for value in assignment["all_after_peel_read_ids"]]
    labels = [str(value) for value in assignment["coarse_labels_all_after_peel"]]
    missing_reads = sorted(set(read_ids) - set(frame.index))
    if missing_reads:
        raise RuntimeError(f"Methylation matrix lacks assigned reads for {case['case_id']}: {missing_reads[:3]}")
    frame = frame.loc[read_ids]
    positions = [int(value) for value in frame.columns]
    matrix = frame.to_numpy(dtype=float)
    minimum_called = max(3, int(np.ceil(matrix.shape[0] * 0.2)))
    keep = np.sum(np.isfinite(matrix), axis=0) >= minimum_called
    matrix = matrix[:, keep]
    positions = [position for position, selected in zip(positions, keep) if selected]
    if matrix.shape[1] > 80:
        selection = np.unique(np.linspace(0, matrix.shape[1] - 1, 80).round().astype(int))
        matrix = matrix[:, selection]
        positions = [positions[index] for index in selection]

    unique_labels = sorted(set(labels))
    label_names = {label: f"Methyl group {chr(65 + index)}" for index, label in enumerate(unique_labels)}
    order = sorted(
        range(len(read_ids)),
        key=lambda index: (unique_labels.index(labels[index]), np.nanmean(matrix[index]) if np.isfinite(matrix[index]).any() else -1),
    )
    matrix = matrix[order]
    ordered_labels = [labels[index] for index in order]
    return matrix, positions, ordered_labels, label_names


def draw_case_heatmap(ax: plt.Axes, case: dict[str, Any], assignment: dict[str, Any]) -> None:
    matrix, positions, labels, label_names = load_case_matrix(case, assignment)
    cmap = mpl.colormaps["cividis"].copy()
    cmap.set_bad(COLORS["missing"])
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    unique_labels = sorted(label_names)
    numeric_labels = np.array([unique_labels.index(label) for label in labels])[:, None]
    group_ax = ax.inset_axes((-0.055, 0, 0.025, 1), transform=ax.transAxes)
    group_ax.imshow(numeric_labels, aspect="auto", interpolation="nearest", cmap=ListedColormap(GROUP_COLORS[: len(unique_labels)]), vmin=0, vmax=max(1, len(unique_labels) - 1))
    group_ax.set_axis_off()

    boundaries = []
    for index in range(1, len(labels)):
        if labels[index] != labels[index - 1]:
            boundaries.append(index - 0.5)
    for boundary in boundaries:
        ax.axhline(boundary, color="white", linewidth=1.2)
    counts = Counter(labels)
    group_text = " | ".join(f"{label_names[label]} n={counts[label]}" for label in unique_labels)
    site = f"{case['sample']}  {case['chrom']}:{case['pos']:,} {case['ref']}>{case['alt']}"
    ax.set_title(f"{site}\n{group_text}", loc="left", fontsize=9.5, fontweight="bold")
    ax.set_ylabel(f"ALT reads (n={matrix.shape[0]})")
    if positions:
        tick_indices = np.unique(np.linspace(0, len(positions) - 1, min(4, len(positions))).round().astype(int))
        ax.set_xticks(tick_indices, [f"{positions[index] / 1e6:.3f} Mb" for index in tick_indices], rotation=0)
    ax.set_xlabel("CpG genomic position")
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.spines[:].set_visible(False)
    hp_text = ", ".join(f"HP {key}: {value}" for key, value in case["alt_hp_counts"].items())
    topo = case["topology_context"].replace("PRIMARY_", "").replace("_", " ").title()
    ax.text(0, -0.22, f"AF={case['caller_af']:.3f}; normal AF={(case['normal_af'] or 0):.3f}; {hp_text}\nRegional context: {topo}", transform=ax.transAxes, fontsize=7.5, color=COLORS["muted"], va="top")


def case_heatmap_grid(report: dict[str, Any], assignments: dict[tuple[str, str, int], dict[str, Any]], case_ids: list[str], title: str, output: Path, columns: int = 2) -> None:
    cases_by_id = {case["case_id"]: case for case in report["case_catalog"]}
    rows = int(np.ceil(len(case_ids) / columns))
    figure_height = 5.8 * rows + (0.7 if rows == 1 else 0)
    fig, axes = plt.subplots(rows, columns, figsize=(7.2 * columns, figure_height), squeeze=False)
    for ax, case_id in zip(axes.flat, case_ids):
        case = cases_by_id[case_id]
        key = (case["sample"], case["chrom"], int(case["pos"]))
        if key not in assignments:
            raise RuntimeError(f"No frozen stable assignment for case {case_id}")
        draw_case_heatmap(ax, case, assignments[key])
    for ax in axes.flat[len(case_ids):]:
        ax.set_visible(False)
    fig.suptitle(title, x=0.02, y=0.985, ha="left", fontsize=14, fontweight="bold")
    colorbar_ax = fig.add_axes([0.945, 0.36, 0.009, 0.24])
    scalar = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap=mpl.colormaps["cividis"])
    colorbar = fig.colorbar(scalar, cax=colorbar_ax)
    colorbar.set_label("Methylation probability", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)
    add_source(fig, "Actual per-read InterSubMod methylation matrices; rows grouped by frozen Python null-gated assignment. Similar-looking cluster and HP labels are distinct.")
    top = 0.79 if rows == 1 else 0.89
    fig.subplots_adjust(left=0.08, right=0.925, bottom=0.14, top=top, hspace=0.58, wspace=0.24)
    save(fig, output)


def write_manifest(paths: list[Path], output: Path) -> None:
    payload = {
        "schema_name": "intersubmod.single_fp_alt_multicluster_report_figures",
        "schema_version": "1.0.0",
        "figures": [
            {"path": str(path), "size_bytes": path.stat().st_size, "exists": path.exists()} for path in paths
        ],
        "pass": all(path.exists() and path.stat().st_size > 10000 for path in paths),
    }
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dataset", type=Path, default=TOPIC_ROOT / "results" / "report_dataset_v1" / "report_dataset.json")
    parser.add_argument("--assignments", type=Path, default=TOPIC_ROOT / "results" / "focal_alt_multicluster" / "latest_full_v3_frozen" / "latest_stable_multigroup_read_assignments.jsonl")
    parser.add_argument("--output-dir", type=Path, default=TOPIC_ROOT / "figures")
    args = parser.parse_args()

    setup_style()
    report = json.loads(args.report_dataset.read_text(encoding="utf-8"))
    assignments = load_assignments(args.assignments)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = [args.output_dir / f"{index:02d}_{name}.png" for index, name in enumerate(
        [
            "evidence_chain",
            "matched_tp_specificity",
            "per_sample_risk_difference",
            "strict_sensitivity",
            "topology_context",
            "longphase_hp_tags",
            "coverage_dependence",
            "legacy_compatibility",
            "technical_replication_heatmaps",
            "topology_contrast_heatmaps",
            "confounder_heatmaps",
        ],
        start=1,
    )]

    evidence_chain(report, outputs[0])
    matched_specificity(report, outputs[1])
    per_sample_risk_difference(report, outputs[2])
    strict_sensitivity(report, outputs[3])
    topology_context(report, outputs[4])
    hp_tags(report, outputs[5])
    coverage_dependence(report, outputs[6])
    legacy_compatibility(report, outputs[7])
    case_heatmap_grid(
        report,
        assignments,
        ["replicated_chr1_hcc", "replicated_chr1_dorado", "replicated_chr9_hcc", "replicated_chr9_dorado"],
        "Technical repeat: two exact FP sites recur, but one changes from two to three methyl groups",
        outputs[8],
    )
    case_heatmap_grid(
        report,
        assignments,
        ["linear_context_not_lineage", "branching_context_same_endpoint"],
        "The same high-threshold methylation endpoint occurs in linear and branching regional contexts",
        outputs[9],
    )
    case_heatmap_grid(
        report,
        assignments,
        ["hp_axis_counterexample", "technical_axis_counterexample", "strict_branching_candidate"],
        "Case falsification: HP/technical confounds and a strict candidate with branching context",
        outputs[10],
        columns=3,
    )
    write_manifest(outputs, args.output_dir / "figure_manifest.json")


if __name__ == "__main__":
    main()
