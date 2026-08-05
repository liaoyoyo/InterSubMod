#!/usr/bin/env python3
"""Build v2 figure set for LOH Round 1 cross-sample audit."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    CURRENT_CLASSES_V2,
    LEGACY_CLASSES,
    LOH_LEGACY_CLASSES,
    UNKNOWN_CURRENT_CLASS,
    VERIFICATION_PROVENANCE_COLUMNS,
    SchemaContractError,
    extract_provenance_frame,
    read_evidence,
    select_current_view,
    select_legacy_view,
    select_loh_legacy,
)


RNG_SEED = 20260329
VIOLIN_MAX_POINTS = 4000

ROUND1_DIR_DEFAULT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260327_loh_round1_cross_sample_audit"
)
OUTPUT_DIR_DEFAULT = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/"
    "20260329_loh_round1_cross_sample_audit_v2_figures"
)

SAMPLE_ORDER = [
    "HCC1395_HKU_5kHz",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
MODE_ORDER = ["paired", "to"]

BASE_COLORS = {
    "Noise": "#6b7280",
    "Weak": "#2a9d8f",
    "Strong": "#2563eb",
    "Subclone": "#d97706",
    "Unknown": "#64748b",
    "Strong_Bidirectional": "#1d4ed8",
    "ClusterFirstOnly": "#d97706",
    "LOH-Structure": "#7c3aed",
    "MultiGroupNoLabel": "#0891b2",
    "LabelShift": "#0f766e",
    "PermanovaLocation": "#2563eb",
    "StructureNoLabel": "#4f46e5",
    "DispersionStructure": "#be123c",
    "Noise_Uniform": "#6b7280",
    "Noise_Chaotic": "#4b5563",
    "Noise_Uncorrelated": "#9ca3af",
    UNKNOWN_CURRENT_CLASS: "#64748b",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round1-dir", default=str(ROUND1_DIR_DEFAULT))
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR_DEFAULT))
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize an old Round 1 workspace without schema metadata.",
    )
    return parser.parse_args()


def ordered_pairs(df: pd.DataFrame) -> list[tuple[str, str]]:
    available = {(row["sample_label"], row["mode"]) for _, row in df[["sample_label", "mode"]].drop_duplicates().iterrows()}
    return [(sample, mode) for sample in SAMPLE_ORDER for mode in MODE_ORDER if (sample, mode) in available]


def blend(color: str, target: str, alpha: float) -> tuple[float, float, float]:
    rgb = np.array(mcolors.to_rgb(color))
    target_rgb = np.array(mcolors.to_rgb(target))
    out = (1.0 - alpha) * rgb + alpha * target_rgb
    return tuple(out.tolist())


def save_fig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def validate_round1_provenance(
    verif_df: pd.DataFrame,
    all_df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> dict[str, object]:
    """Validate derived grouping identity and the source-row provenance receipt."""

    if "VerificationProvenanceStatus" not in all_df.columns:
        if not allow_unversioned_v1:
            raise SchemaContractError(
                "LOH v2 figures: upstream schema metadata missing; explicit "
                "--allow-unversioned-v1 is required for a historical workspace"
            )
        historical = all_df[["verification_class", "loh_subtype"]].rename(
            columns={"verification_class": "VerificationClass", "loh_subtype": "LOH_Subtype"}
        )
        extract_provenance_frame(historical, allow_unversioned_v1=True)
        all_df["VerificationProvenanceStatus"] = "UNVERSIONED_DERIVED"
        verif_df["VerificationProvenanceStatus"] = "UNVERSIONED_DERIVED"
        return {
            "schema_status": "UNVERSIONED_DERIVED",
            "class_order": list(LEGACY_CLASSES),
            "authorization": "--allow-unversioned-v1",
            "provenance_fields": [],
        }

    statuses = sorted(
        all_df["VerificationProvenanceStatus"].dropna().astype(str).unique().tolist()
    )
    if len(statuses) != 1:
        raise SchemaContractError(f"LOH v2 figures: mixed source schema statuses: {statuses}")
    status = statuses[0]
    summary_statuses = sorted(
        verif_df.get("VerificationProvenanceStatus", pd.Series(dtype=str))
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    if summary_statuses != [status]:
        raise SchemaContractError(
            "LOH v2 figures: derived verification table status does not match source rows: "
            f"source={statuses}, derived={summary_statuses}"
        )

    if status == "V2":
        required_summary_metadata = [
            "VerificationSchemaVersion",
            "VerificationProvenanceSourceField",
            "LOHProvenanceSourceField",
        ]
        missing_summary_metadata = [
            field for field in required_summary_metadata if field not in verif_df.columns
        ]
        if missing_summary_metadata:
            raise SchemaContractError(
                "LOH v2 figures: derived verification table missing schema metadata: "
                + ", ".join(missing_summary_metadata)
            )
        current = select_current_view(all_df)
        select_legacy_view(all_df)
        read_evidence(all_df)
        loh = select_loh_legacy(all_df)
        missing = [field for field in VERIFICATION_PROVENANCE_COLUMNS if field not in all_df.columns]
        if missing:
            raise SchemaContractError(
                "LOH v2 figures: source rows missing provenance fields: " + ", ".join(missing)
            )
        known = current.values != UNKNOWN_CURRENT_CLASS
        if known.any():
            extract_provenance_frame(all_df.loc[known])
        if not all_df["verification_class"].astype(str).equals(current.values.astype(str)):
            raise SchemaContractError(
                "LOH v2 figures: derived verification_class is not the validated current C2 view"
            )
        if not all_df["loh_subtype"].astype(str).equals(loh.values.astype(str)):
            raise SchemaContractError(
                "LOH v2 figures: derived loh_subtype is not canonical LOH_Subtype_LegacyVC"
            )
        current_sources = sorted(
            all_df["VerificationProvenanceSourceField"].dropna().astype(str).unique().tolist()
        )
        loh_sources = sorted(
            all_df["LOHProvenanceSourceField"].dropna().astype(str).unique().tolist()
        )
        summary_current_sources = sorted(
            verif_df["VerificationProvenanceSourceField"].dropna().astype(str).unique().tolist()
        )
        summary_loh_sources = sorted(
            verif_df["LOHProvenanceSourceField"].dropna().astype(str).unique().tolist()
        )
        if current_sources != ["VerificationClass"] or summary_current_sources != current_sources:
            raise SchemaContractError(
                "LOH v2 figures: current-view provenance source must be VerificationClass"
            )
        if loh_sources != ["LOH_Subtype_LegacyVC"] or summary_loh_sources != loh_sources:
            raise SchemaContractError(
                "LOH v2 figures: LOH provenance source must be LOH_Subtype_LegacyVC"
            )
        versions = pd.to_numeric(verif_df["VerificationSchemaVersion"], errors="coerce")
        if versions.isna().any() or sorted(versions.unique().tolist()) != [2.0]:
            raise SchemaContractError(
                "LOH v2 figures: derived verification table must carry VerificationSchemaVersion=2"
            )
        invalid_classes = sorted(
            set(verif_df["verification_class"].dropna().astype(str))
            - set(CURRENT_CLASSES_V2)
            - {UNKNOWN_CURRENT_CLASS}
        )
        invalid_loh = sorted(
            set(verif_df["loh_subtype"].dropna().astype(str)) - set(LOH_LEGACY_CLASSES)
        )
        if invalid_classes or invalid_loh:
            raise SchemaContractError(
                f"LOH v2 figures: invalid derived classes current={invalid_classes}, loh={invalid_loh}"
            )
        return {
            "schema_status": status,
            "class_order": list(CURRENT_CLASSES_V2) + [UNKNOWN_CURRENT_CLASS],
            "schema_versions": [2],
            "current_selection_field": "VerificationClass",
            "loh_selection_field": "LOH_Subtype_LegacyVC",
            "provenance_fields": list(VERIFICATION_PROVENANCE_COLUMNS),
        }

    if status != "UNVERSIONED_V1":
        raise SchemaContractError(f"LOH v2 figures: unsupported schema status {status!r}")
    if not allow_unversioned_v1:
        raise SchemaContractError(
            "LOH v2 figures: UNVERSIONED_V1 input requires --allow-unversioned-v1 authorization"
        )
    historical = all_df[["verification_class", "loh_subtype"]].rename(
        columns={"verification_class": "VerificationClass", "loh_subtype": "LOH_Subtype"}
    )
    extract_provenance_frame(historical, allow_unversioned_v1=True)
    return {
        "schema_status": status,
        "class_order": list(LEGACY_CLASSES),
        "authorization": "--allow-unversioned-v1",
        "provenance_fields": ["VerificationClass", "LOH_Subtype"],
    }


def load_round1_tables(
    round1_dir: Path,
    allow_unversioned_v1: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    loh_df = pd.read_csv(round1_dir / "loh_enrichment_by_sample_mode.tsv", sep="\t")
    verif_df = pd.read_csv(round1_dir / "verificationclass_by_loh_subtype.tsv", sep="\t")
    usecols = [
        "sample_label",
        "mode",
        "truth_label",
        "effective_hp_reads",
        "hp_ratio_core",
        "num_reads",
        "verification_class",
        "loh_subtype",
        "core_loh_like",
        "allele_delta",
        "pairwise_median_dist",
        "quality_score",
        "hp_assign_rate",
        "hp0_ratio",
        "hp3_ratio",
        "caller_af",
    ]
    all_path = round1_dir / "all_region_rows.tsv.gz"
    available = pd.read_csv(all_path, sep="\t", nrows=0).columns.tolist()
    provenance_candidates = [
        *VERIFICATION_PROVENANCE_COLUMNS,
        "VerificationProvenanceStatus",
        "VerificationProvenanceSourceField",
        "VerificationProvenanceWarnings",
        "VerificationUnknownCurrentCounts",
        "LOHProvenanceSourceField",
    ]
    required_missing = [field for field in usecols if field not in available]
    if required_missing:
        raise SchemaContractError(
            "LOH v2 figures: source table missing required columns: " + ", ".join(required_missing)
        )
    selected_cols = usecols + [field for field in provenance_candidates if field in available]
    all_df = pd.read_csv(all_path, sep="\t", usecols=selected_cols, low_memory=False)
    receipt = validate_round1_provenance(
        verif_df,
        all_df,
        allow_unversioned_v1=allow_unversioned_v1,
    )
    all_df.attrs["verification_contract"] = receipt
    verif_df.attrs["verification_contract"] = receipt
    numeric_cols = [
        "effective_hp_reads",
        "hp_ratio_core",
        "num_reads",
        "allele_delta",
        "pairwise_median_dist",
        "quality_score",
        "hp_assign_rate",
        "hp0_ratio",
        "hp3_ratio",
        "caller_af",
    ]
    for col in numeric_cols:
        all_df[col] = pd.to_numeric(all_df[col], errors="coerce")
    all_df["core_loh_like"] = all_df["core_loh_like"].astype(bool)
    return loh_df, verif_df, all_df


def plot_fig01(loh_df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    bar_width = 0.38
    x = np.arange(len(SAMPLE_ORDER))

    for ax, mode in zip(axes, MODE_ORDER):
        sub = loh_df[loh_df["mode"] == mode].copy()
        sub["sample_label"] = pd.Categorical(sub["sample_label"], categories=SAMPLE_ORDER, ordered=True)
        sub = sub.sort_values("sample_label")
        tp = sub["tp_core_loh_like_frac"].to_numpy()
        fp = sub["fp_core_loh_like_frac"].to_numpy()
        tp_n = sub["tp_total"].to_numpy()
        enrich = sub["core_loh_fp_vs_tp_enrichment"].to_numpy()
        fp_n = sub["fp_total"].to_numpy()

        ax.bar(x - bar_width / 2, tp, width=bar_width, color="#60a5fa", label="TP LOH-like fraction")
        ax.bar(x + bar_width / 2, fp, width=bar_width, color="#f59e0b", label="FP LOH-like fraction")
        for xpos, frac, n in zip(x, tp, tp_n):
            ax.text(xpos - bar_width / 2, frac + 0.012, f"TP={int(n)}", ha="center", va="bottom", fontsize=7)
        for xpos, frac, ratio, n in zip(x, fp, enrich, fp_n):
            ax.text(xpos + bar_width / 2, frac + 0.012, f"{ratio:.2f}x\nFP={int(n)}", ha="center", va="bottom", fontsize=7)
        ax.set_ylim(0, max(0.8, float(np.nanmax(np.concatenate([tp, fp]))) + 0.16))
        ax.set_ylabel("fraction")
        ax.set_title(f"{mode}: TP vs FP LOH-like fraction")
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="upper left", fontsize=8)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(SAMPLE_ORDER, rotation=25, ha="right")
    fig.suptitle("Fig01 v2 LOH-like fraction overview", y=1.02)
    fig.tight_layout()
    save_fig(fig, out_path)


def plot_fig02(all_df: pd.DataFrame, out_path: Path) -> None:
    truth_order = ["TP", "FP"]
    mode_order = ["paired", "to"]
    bins = np.linspace(0, 1, 21)
    rows = len(SAMPLE_ORDER) * len(mode_order)
    fig, axes = plt.subplots(rows, 2, figsize=(12, max(2.0 * rows, 16)), squeeze=False, sharex=True)
    fill_map = {"TP": "#2563eb", "FP": "#d97706"}

    for sample_idx, sample_label in enumerate(SAMPLE_ORDER):
        sample_payload: dict[tuple[str, str], tuple[np.ndarray, int]] = {}
        sample_ymax = 0

        for mode in mode_order:
            for truth_label in truth_order:
                subset = all_df[
                    (all_df["sample_label"] == sample_label)
                    & (all_df["mode"] == mode)
                    & (all_df["truth_label"] == truth_label)
                    & (all_df["effective_hp_reads"] > 0)
                ]
                vals = subset["hp_ratio_core"].dropna().to_numpy()
                hist_counts, _ = np.histogram(vals, bins=bins)
                sample_payload[(mode, truth_label)] = (vals, len(vals))
                sample_ymax = max(sample_ymax, int(hist_counts.max()) if len(hist_counts) else 0)

        sample_ymax = max(5, int(math.ceil(sample_ymax * 1.08)))

        for mode_idx, mode in enumerate(mode_order):
            row_idx = sample_idx * len(mode_order) + mode_idx
            for col_idx, truth_label in enumerate(truth_order):
                ax = axes[row_idx][col_idx]
                vals, nvals = sample_payload[(mode, truth_label)]
                ax.hist(vals, bins=bins, color=fill_map[truth_label], alpha=0.8, edgecolor="white", linewidth=0.5)
                ax.set_ylim(0, sample_ymax)
                ax.grid(axis="y", alpha=0.2)

                if row_idx == 0:
                    ax.set_title(truth_label, fontsize=11)

                if col_idx == 0:
                    ax.set_ylabel(f"{sample_label}\n{mode}\ncount", fontsize=8)

                ax.text(
                    0.98,
                    0.92,
                    f"n={nvals}",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=8,
                    bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none", "pad": 1.5},
                )

                if row_idx == rows - 1:
                    ax.set_xlabel("hp_ratio_core")

    fig.suptitle(
        "Fig02 v2 hp_ratio_core distribution by count\nleft=TP, right=FP; paired above same-sample TO; same-sample y-scale matched",
        y=1.01,
    )
    fig.tight_layout()
    save_fig(fig, out_path)


def plot_fig03(all_df: pd.DataFrame, out_path: Path) -> None:
    pairs = ordered_pairs(all_df)
    cols = 2
    rows = math.ceil(len(pairs) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(14, max(3.4 * rows, 4)), squeeze=False)

    for ax in axes.flatten():
        ax.axis("off")

    for idx, (sample_label, mode) in enumerate(pairs):
        ax = axes[idx // cols][idx % cols]
        ax.axis("on")
        subset = all_df[(all_df["sample_label"] == sample_label) & (all_df["mode"] == mode)].copy()
        xcap = float(np.nanquantile(subset["effective_hp_reads"].fillna(0), 0.995))
        xcap = max(20.0, min(xcap, float(np.nanmax(subset["effective_hp_reads"].fillna(0)))))
        fp = subset[subset["truth_label"] == "FP"]
        tp = subset[subset["truth_label"] == "TP"]

        ax.scatter(
            tp["effective_hp_reads"],
            tp["hp_ratio_core"],
            s=5,
            alpha=0.10,
            color="#0000ff",
            label="TP",
            zorder=2,
            rasterized=True,
        )
        ax.scatter(
            fp["effective_hp_reads"],
            fp["hp_ratio_core"],
            s=5,
            alpha=0.10,
            color="#ff0000",
            label="FP",
            zorder=3,
            rasterized=True,
        )
        ax.set_xlim(0, xcap)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"{sample_label} {mode} (x<=p99.5={xcap:.0f})")
        ax.set_xlabel("effective_hp_reads")
        ax.set_ylabel("hp_ratio_core")
        ax.grid(alpha=0.15)
        ax.legend(fontsize=7, loc="upper right")

    fig.suptitle("Fig03 v2 effective HP vs hp_ratio_core", y=1.02)
    fig.tight_layout()
    save_fig(fig, out_path)


def combo_color(verification_class: str, loh_subtype: str) -> tuple[float, float, float]:
    base = BASE_COLORS.get(verification_class, BASE_COLORS["Unknown"])
    if loh_subtype == "None":
        return blend(base, "#ffffff", 0.58)
    return blend(base, "#000000", 0.12)


def plot_fig04(verif_df: pd.DataFrame, out_path: Path) -> None:
    df = verif_df.copy()
    df["sample_mode_truth"] = df["sample_label"] + "|" + df["mode"] + "|" + df["truth_label"]
    order = [f"{sample}|{mode}|{truth}" for sample in SAMPLE_ORDER for mode in MODE_ORDER for truth in ["TP", "FP"]]
    order = [x for x in order if x in set(df["sample_mode_truth"])]

    receipt = verif_df.attrs.get("verification_contract", {})
    class_order = receipt.get("class_order", list(LEGACY_CLASSES))
    combo_order = [
        (cls, subtype)
        for cls in class_order
        for subtype in LOH_LEGACY_CLASSES
    ]
    combo_order = [x for x in combo_order if ((df["verification_class"] == x[0]) & (df["loh_subtype"] == x[1])).any()]

    pivot = pd.DataFrame(index=order)
    for cls, subtype in combo_order:
        key = f"{cls}|{subtype}"
        counts = (
            df[(df["verification_class"] == cls) & (df["loh_subtype"] == subtype)]
            .groupby("sample_mode_truth")["count"]
            .sum()
            .reindex(order, fill_value=0)
        )
        pivot[key] = counts

    frac = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(max(14, 0.5 * len(order)), 6.5))
    bottom = np.zeros(len(frac.index))
    for key in frac.columns:
        cls, subtype = key.split("|", 1)
        values = frac[key].to_numpy()
        label = f"{cls} / {'non-LOH' if subtype == 'None' else subtype}"
        ax.bar(
            np.arange(len(frac.index)),
            values,
            bottom=bottom,
            color=combo_color(cls, subtype),
            label=label,
            width=0.85,
        )
        bottom += values

    ax.set_xticks(np.arange(len(frac.index)))
    ax.set_xticklabels(frac.index.tolist(), rotation=55, ha="right")
    ax.set_ylabel("fraction within sample|mode|truth")
    ax.set_title("Fig04 v2 VerificationClass base color, LOH as darker shade")
    ax.legend(fontsize=6, ncol=2, loc="upper left")
    fig.tight_layout()
    save_fig(fig, out_path)


def violin_ready(values: np.ndarray, rng: np.random.Generator, max_points: int = VIOLIN_MAX_POINTS) -> np.ndarray:
    vals = values[~np.isnan(values)]
    if len(vals) == 0:
        return np.array([0.0, 0.0])
    if len(vals) > max_points:
        vals = rng.choice(vals, size=max_points, replace=False)
    if len(vals) == 1:
        return np.array([vals[0], vals[0]])
    return vals


def plot_fig05(all_df: pd.DataFrame, out_path: Path) -> None:
    df = all_df.copy()
    rng = np.random.default_rng(RNG_SEED)
    df["group_label"] = (
        df["mode"]
        + "\n"
        + df["truth_label"]
        + "\n"
        + np.where(df["core_loh_like"], "LOH", "non-LOH")
    )
    group_order = [
        "paired\nTP\nLOH",
        "paired\nTP\nnon-LOH",
        "paired\nFP\nLOH",
        "paired\nFP\nnon-LOH",
        "to\nTP\nLOH",
        "to\nTP\nnon-LOH",
        "to\nFP\nLOH",
        "to\nFP\nnon-LOH",
    ]
    group_colors = {
        "paired\nTP\nLOH": "#2563eb",
        "paired\nTP\nnon-LOH": "#93c5fd",
        "paired\nFP\nLOH": "#1d4ed8",
        "paired\nFP\nnon-LOH": "#bfdbfe",
        "to\nTP\nLOH": "#d97706",
        "to\nTP\nnon-LOH": "#fcd34d",
        "to\nFP\nLOH": "#b45309",
        "to\nFP\nnon-LOH": "#fde68a",
    }
    features = [
        ("allele_delta", "AlleleDelta"),
        ("pairwise_median_dist", "PairwiseMedianDist"),
        ("quality_score", "Quality_Score"),
        ("hp_assign_rate", "hp_assign_rate"),
        ("hp0_ratio", "hp0_ratio"),
        ("hp3_ratio", "hp3_ratio (TO structural zero)"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(17, 11.5), squeeze=False)
    for ax, (column, title) in zip(axes.flatten(), features):
        data = [violin_ready(df[df["group_label"] == group][column].to_numpy(dtype=float), rng=rng) for group in group_order]
        parts = ax.violinplot(data, positions=np.arange(1, len(group_order) + 1), showmeans=False, showmedians=True, widths=0.85)
        for body, group in zip(parts["bodies"], group_order):
            body.set_facecolor(group_colors[group])
            body.set_edgecolor("#374151")
            body.set_alpha(0.7)
        parts["cmedians"].set_color("#111827")
        parts["cmedians"].set_linewidth(1.2)
        ax.set_title(title)
        ax.set_xticks(np.arange(1, len(group_order) + 1))
        ax.set_xticklabels(group_order, fontsize=8)
        ax.grid(axis="y", alpha=0.2)

    fig.suptitle("Fig05 v2 LOH / non-LOH feature violin plots", y=1.02)
    fig.tight_layout()
    save_fig(fig, out_path)


def assign_hp_ratio_bins_v2(df: pd.DataFrame) -> pd.Series:
    labels = ["<0.05", "0.05-0.10", "0.10-0.20", "0.20-0.40", "0.40-0.60", "0.60-0.80", "0.80-0.90", "0.90-0.95", ">0.95"]
    valid = (df["effective_hp_reads"] > 0) & df["hp_ratio_core"].notna()
    bins = pd.Series("no_eff", index=df.index, dtype="object")
    bins.loc[valid] = pd.cut(
        df.loc[valid, "hp_ratio_core"],
        bins=[-np.inf, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 0.90, 0.95, np.inf],
        labels=labels,
        right=False,
    ).astype("object")
    return bins.fillna("no_eff")


def assign_effective_hp_bins_v2(series: pd.Series) -> pd.Series:
    labels = ["0", "1-9", "10-29", "30-49", "50-74", "75-99", "100-149", ">=150"]
    bins = pd.cut(
        series.fillna(-1),
        bins=[-np.inf, 0.5, 10, 30, 50, 75, 100, 150, np.inf],
        labels=labels,
        right=False,
    )
    return bins.astype("object").fillna("0")


def heatmap_table(df: pd.DataFrame, value_col: str, bin_order: list[str], label: str) -> pd.DataFrame:
    sample_mode = pd.DataFrame(
        [(sample, mode) for sample in SAMPLE_ORDER for mode in MODE_ORDER],
        columns=["sample_label", "mode"],
    )
    totals = df.groupby(["sample_label", "mode"]).size().rename("total").reset_index()
    counts = (
        df.groupby(["sample_label", "mode", value_col])
        .size()
        .rename("count")
        .reset_index()
        .rename(columns={value_col: "bin_label"})
    )
    grid = sample_mode.assign(key=1).merge(pd.DataFrame({"bin_label": bin_order, "key": 1}), on="key").drop(columns=["key"])
    grid = grid.merge(totals, on=["sample_label", "mode"], how="left")
    grid = grid.merge(counts, on=["sample_label", "mode", "bin_label"], how="left")
    grid["sample_mode"] = grid["sample_label"] + "|" + grid["mode"]
    grid["metric"] = label
    grid["total"] = grid["total"].fillna(0).astype(int)
    grid["count"] = grid["count"].fillna(0).astype(int)
    grid["fraction"] = np.where(grid["total"] > 0, grid["count"] / grid["total"], 0.0)
    return grid[["sample_mode", "bin_label", "count", "total", "fraction", "metric"]]


def plot_single_heatmap(ax: plt.Axes, heat_df: pd.DataFrame, bin_order: list[str], title: str) -> None:
    order = [f"{sample}|{mode}" for sample in SAMPLE_ORDER for mode in MODE_ORDER if f"{sample}|{mode}" in set(heat_df["sample_mode"])]
    fraction_pivot = (
        heat_df.pivot(index="sample_mode", columns="bin_label", values="fraction")
        .reindex(index=order, columns=bin_order)
        .fillna(0.0)
    )
    count_pivot = (
        heat_df.pivot(index="sample_mode", columns="bin_label", values="count")
        .reindex(index=order, columns=bin_order)
        .fillna(0)
        .astype(int)
    )
    fraction_values = fraction_pivot.to_numpy()
    im = ax.imshow(fraction_values, aspect="auto", cmap="magma", vmin=0, vmax=max(0.35, float(fraction_values.max())))
    ax.set_title(title)
    ax.set_xticks(np.arange(len(bin_order)))
    ax.set_xticklabels(bin_order, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    for i in range(fraction_pivot.shape[0]):
        for j in range(fraction_pivot.shape[1]):
            val = fraction_pivot.iat[i, j]
            count = count_pivot.iat[i, j]
            text = f"{val * 100:.0f}%\n{count}"
            ax.text(j, i, text, ha="center", va="center", color="white" if val < 0.2 else "black", fontsize=6)
    return im


def plot_fig06(all_df: pd.DataFrame, out_path: Path) -> None:
    df = all_df.copy()
    df["hp_ratio_bin_v2"] = assign_hp_ratio_bins_v2(df)
    df["effective_hp_bin_v2"] = assign_effective_hp_bins_v2(df["effective_hp_reads"])
    hp_bins = ["no_eff", "<0.05", "0.05-0.10", "0.10-0.20", "0.20-0.40", "0.40-0.60", "0.60-0.80", "0.80-0.90", "0.90-0.95", ">0.95"]
    eff_bins = ["0", "1-9", "10-29", "30-49", "50-74", "75-99", "100-149", ">=150"]
    hp_heat = heatmap_table(df, "hp_ratio_bin_v2", hp_bins, "hp_ratio_core")
    eff_heat = heatmap_table(df, "effective_hp_bin_v2", eff_bins, "effective_hp_reads")

    fig, axes = plt.subplots(1, 2, figsize=(18, 7.5), squeeze=False)
    im1 = plot_single_heatmap(axes[0][0], hp_heat, hp_bins, "hp_ratio_core bin fraction (%)")
    im2 = plot_single_heatmap(axes[0][1], eff_heat, eff_bins, "effective_hp_reads bin fraction (%)")
    cbar1 = fig.colorbar(im1, ax=axes[0][0], fraction=0.046, pad=0.04)
    cbar2 = fig.colorbar(im2, ax=axes[0][1], fraction=0.046, pad=0.04)
    cbar1.set_label("fraction within sample|mode")
    cbar2.set_label("fraction within sample|mode")
    fig.suptitle("Fig06 v2 sample bin heatmap (cell text: % and n within sample|mode)", y=1.02)
    fig.tight_layout()
    save_fig(fig, out_path)


def main() -> None:
    args = parse_args()
    round1_dir = Path(args.round1_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    loh_df, verif_df, all_df = load_round1_tables(
        round1_dir,
        allow_unversioned_v1=args.allow_unversioned_v1,
    )
    receipt = all_df.attrs.get("verification_contract", {})
    (output_dir / "verification_schema_contract.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    plot_fig01(loh_df, output_dir / "fig01_loh_like_fraction_overview.png")
    plot_fig02(all_df, output_dir / "fig02_hp_ratio_core_distribution.png")
    plot_fig03(all_df, output_dir / "fig03_effective_hp_vs_hp_ratio_scatter.png")
    plot_fig04(verif_df, output_dir / "fig04_verificationclass_lohsubtype_structure.png")
    plot_fig05(all_df, output_dir / "fig05_loh_nonloh_feature_boxplots.png")
    plot_fig06(all_df, output_dir / "fig06_sample_bin_heatmap.png")


if __name__ == "__main__":
    main()
