#!/usr/bin/env python3
"""Figures for TP-vs-FP structure/label characterization from verified JSON assets."""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.plot_setup import setup_plot_style  # noqa: E402
from scripts.lib.verification_schema_contract import (  # noqa: E402
    LEGACY_CLASSES,
    UNKNOWN_LEGACY_CLASS,
    SchemaContractError,
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/"
    "2026/06/20260619_tp_fp_structure_association_HCC1395"
)
LEGACY_DIMENSION = "legacy_VerificationClass"
TP_COLOR, FP_COLOR = "#2c7fb8", "#d95f0e"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument(
        "--assets-dir",
        default=None,
        help="Verified asset directory (default: <output-dir>/_assets).",
    )
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize assets whose upstream legacy view used unversioned v1 input.",
    )
    return parser.parse_args()


def _read_json(path: Path):
    try:
        with path.open() as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SchemaContractError(f"figure asset {path}: cannot read valid JSON: {exc}") from exc


def _validate_taxonomy_metadata(meta, allow_unversioned_v1=False):
    if not isinstance(meta, dict):
        raise SchemaContractError("figure run_meta: top-level JSON must be an object")
    selection = meta.get("verification_selection")
    if not isinstance(selection, dict):
        raise SchemaContractError("figure run_meta: verification_selection metadata is required")
    required = {"selection_field", "schema_status", "categories", "unknown_counts", "warnings"}
    missing = sorted(required - set(selection))
    if missing:
        raise SchemaContractError(
            "figure run_meta: verification_selection is missing: " + ", ".join(missing)
        )
    if meta.get("verification_selection_values") != list(LEGACY_CLASSES):
        raise SchemaContractError(
            "figure run_meta: verification_selection_values must be the frozen legacy four-state order"
        )
    if selection["categories"] != list(LEGACY_CLASSES):
        raise SchemaContractError("figure run_meta: legacy category identity/order is not canonical")

    status = selection["schema_status"]
    field = selection["selection_field"]
    if status == "LEGACY_EXPLICIT":
        if field != "VerificationClass_Legacy":
            raise SchemaContractError(
                "figure run_meta: LEGACY_EXPLICIT must select VerificationClass_Legacy"
            )
    elif status == "UNVERSIONED_V1":
        if not allow_unversioned_v1:
            raise SchemaContractError(
                "figure assets are unversioned v1; --allow-unversioned-v1 is required"
            )
        if field != "VerificationClass" or not selection["warnings"]:
            raise SchemaContractError(
                "figure run_meta: unversioned v1 must identify VerificationClass and retain a warning"
            )
        warnings.warn(
            "UNVERSIONED: rendering explicitly authorized v1 legacy taxonomy assets",
            UserWarning,
            stacklevel=2,
        )
    else:
        raise SchemaContractError(f"figure run_meta: unsupported schema_status {status!r}")
    return selection


def load_verified_assets(assets_dir: Path, allow_unversioned_v1=False):
    """Load all three upstream artifacts and require their taxonomy contract before plotting."""
    assets_dir = Path(assets_dir)
    meta = _read_json(assets_dir / "run_meta.json")
    selection = _validate_taxonomy_metadata(meta, allow_unversioned_v1)
    crosstabs = _read_json(assets_dir / "crosstabs.json")
    confound = _read_json(assets_dir / "confound_control.json")
    if not isinstance(crosstabs, dict) or not isinstance(confound, dict):
        raise SchemaContractError("figure assets: crosstabs/confound_control must be JSON objects")

    if LEGACY_DIMENSION not in crosstabs:
        raise SchemaContractError(
            f"figure crosstabs: required explicit legacy key {LEGACY_DIMENSION!r} is missing"
        )
    try:
        matched = confound["matched_subset"]["crosstabs"]
    except (KeyError, TypeError) as exc:
        raise SchemaContractError("figure confound_control: matched_subset.crosstabs is missing") from exc
    if LEGACY_DIMENSION not in matched:
        raise SchemaContractError(
            f"figure matched crosstabs: required explicit legacy key {LEGACY_DIMENSION!r} is missing"
        )
    return crosstabs, confound, meta, selection


def cells(data, dimension):
    try:
        result = data[dimension]["cells"]
    except (KeyError, TypeError) as exc:
        raise SchemaContractError(f"figure crosstab {dimension!r}: cells are missing") from exc
    if not isinstance(result, list):
        raise SchemaContractError(f"figure crosstab {dimension!r}: cells must be a list")
    return result


def pct_pair(cell_list, order):
    if any(not isinstance(cell, dict) for cell in cell_list):
        raise SchemaContractError("figure crosstab: every cell must be an object")
    by_category = {cell.get("category"): cell for cell in cell_list}
    if len(by_category) != len(cell_list):
        raise SchemaContractError("figure crosstab: category values must be unique")
    missing = [category for category in order if category not in by_category]
    if missing:
        raise SchemaContractError("figure crosstab: missing categories: " + ", ".join(missing))
    for category in order:
        for field in ("pct_tp", "pct_fp"):
            value = by_category[category].get(field)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
                raise SchemaContractError(
                    f"figure crosstab: {category}.{field} must be numeric in [0,1]"
                )
    tp = [by_category[category]["pct_tp"] * 100 for category in order]
    fp = [by_category[category]["pct_fp"] * 100 for category in order]
    return tp, fp


def grouped(ax, order, tp, fp, title, ylab="占該類 %"):
    x = np.arange(len(order))
    width = 0.38
    ax.bar(x - width / 2, tp, width, label="TP", color=TP_COLOR)
    ax.bar(x + width / 2, fp, width, label="FP", color=FP_COLOR)
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=20, ha="right")
    ax.set_title(title)
    ax.set_ylabel(ylab)
    ax.legend(frameon=False)
    for index, (tp_value, fp_value) in enumerate(zip(tp, fp)):
        ax.text(index - width / 2, tp_value + 0.5, f"{tp_value:.1f}", ha="center", va="bottom", fontsize=7)
        ax.text(index + width / 2, fp_value + 0.5, f"{fp_value:.1f}", ha="center", va="bottom", fontsize=7)


def render_figures(crosstabs, confound, output_dir: Path):
    setup_plot_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    order = ["allele_only", "both", "hp_only", "neither"]
    tp, fp = pct_pair(cells(crosstabs, "label_axis_hp_vs_allele"), order)
    grouped(axes[0], order, tp, fp, "HP軸 vs allele軸結構（全集 TP29754/FP627）")
    matched = confound["matched_subset"]["crosstabs"]
    matched_tp, matched_fp = pct_pair(cells(matched, "label_axis_hp_vs_allele"), order)
    grouped(axes[1], order, matched_tp, matched_fp, "覆蓋配對後（TP3121/FP627, KS p≈1）")
    fig.suptitle("TP=allele軸甲基結構 / FP=HP(germline)軸結構 — 覆蓋配對後仍存活", fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_dir / "fig1_label_axis_full_vs_matched.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    legacy_order = list(LEGACY_CLASSES)
    unknown_cell = next(
        (cell for cell in cells(crosstabs, LEGACY_DIMENSION) if cell.get("category") == UNKNOWN_LEGACY_CLASS),
        None,
    )
    if unknown_cell and (unknown_cell.get("n_tp", 0) or unknown_cell.get("n_fp", 0)):
        legacy_order.append(UNKNOWN_LEGACY_CLASS)
    tp, fp = pct_pair(cells(crosstabs, LEGACY_DIMENSION), legacy_order)
    grouped(axes[0], legacy_order, tp, fp, "Verdict 分類（VerificationClass_Legacy）")
    alignment_order = ["aligned", "misaligned", "no_cluster"]
    tp, fp = pct_pair(cells(crosstabs, "alignment_state"), alignment_order)
    grouped(axes[1], alignment_order, tp, fp, "cluster×label 對齊狀態")
    fig.suptitle("FP 偏 Noise/no_cluster；TP 偏 Strong/aligned（皆富集非排他）", fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_dir / "fig2_verdict_alignment.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    delta_order = ["[0,0.05)", "[0.05,0.1)", "[0.1,0.2)", ">=0.2"]
    tp, fp = pct_pair(cells(crosstabs, "hpMergedDelta_bin"), delta_order)
    grouped(axes[0], delta_order, tp, fp, "HP軸 |Δβ| 分箱 — 大值 FP 富集(≥0.2: 7.4×)")
    tp, fp = pct_pair(cells(crosstabs, "alleleDelta_bin"), delta_order)
    grouped(axes[1], delta_order, tp, fp, "allele軸 |Δβ| 分箱 — 大值 TP 偏多")
    fig.suptitle("HP軸大Δβ = FP 簽名；allele軸大Δβ = TP 偏多", fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_dir / "fig3_delta_magnitude.png", bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    assets_dir = Path(args.assets_dir).resolve() if args.assets_dir else output_dir / "_assets"
    try:
        crosstabs, confound, _meta, _selection = load_verified_assets(
            assets_dir,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
        # Validate every category/value consumed below before creating any output.
        label_order = ["allele_only", "both", "hp_only", "neither"]
        pct_pair(cells(crosstabs, "label_axis_hp_vs_allele"), label_order)
        pct_pair(
            cells(confound["matched_subset"]["crosstabs"], "label_axis_hp_vs_allele"),
            label_order,
        )
        legacy_cells = cells(crosstabs, LEGACY_DIMENSION)
        legacy_order = list(LEGACY_CLASSES)
        unknown_cell = next(
            (cell for cell in legacy_cells if cell.get("category") == UNKNOWN_LEGACY_CLASS),
            None,
        )
        if unknown_cell and (unknown_cell.get("n_tp", 0) or unknown_cell.get("n_fp", 0)):
            legacy_order.append(UNKNOWN_LEGACY_CLASS)
        pct_pair(legacy_cells, legacy_order)
        pct_pair(cells(crosstabs, "alignment_state"), ["aligned", "misaligned", "no_cluster"])
        delta_order = ["[0,0.05)", "[0.05,0.1)", "[0.1,0.2)", ">=0.2"]
        pct_pair(cells(crosstabs, "hpMergedDelta_bin"), delta_order)
        pct_pair(cells(crosstabs, "alleleDelta_bin"), delta_order)
    except SchemaContractError as exc:
        raise SystemExit(f"[tp-fp-figures][schema-contract] {exc}") from exc

    render_figures(crosstabs, confound, output_dir)
    print("figs written to", output_dir)
    for filename in (
        "fig1_label_axis_full_vs_matched.png",
        "fig2_verdict_alignment.png",
        "fig3_delta_magnitude.png",
    ):
        path = output_dir / filename
        print(filename, path.stat().st_size if path.exists() else "MISSING")


if __name__ == "__main__":
    main()
