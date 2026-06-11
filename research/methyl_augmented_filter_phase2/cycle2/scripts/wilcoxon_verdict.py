#!/usr/bin/env python3
"""
Phase 2 Cycle 2 Step B4' — Wilcoxon paired signed-rank test for H_C1_5 verdict.

Input:
  cycle2/data/cycle2_cross_sample_delta_f1.tsv (5 samples × 3 modes)

For each mode {transfer_fixed, transfer_sweep, refit}:
  - n = 5 samples
  - Wilcoxon paired signed-rank vs 0 (exact min p with n=5 is 0.0625, achievable
    only if all 5 same direction)
  - Direction concordance (count of +/0/- and majority direction)

H_C1_5 verdict rule (per prompt):
  - POSITIVE      : >=4/5 ΔF1 > 0 AND Wilcoxon p < 0.0625 (requires 5/5 same direction)
  - DIRECTION     : >=4/5 same direction BUT Wilcoxon p >= 0.0625 (n=5 power-limited)
  - NEGATIVE      : <=2/5 same direction (majority opposite or split)
  - MIXED         : 3/5 same direction (no clear majority)

Outputs:
  cycle2/data/cycle2_b4_wilcoxon_summary.tsv  (per-mode Wilcoxon stats + verdict)
  cycle2/figures/cycle2_b4_wilcoxon.png        (per-sample ΔF1 bar with Wilcoxon p annotation)
  cycle2/intermediate/cycle2_b4_summary.json    (full JSON with per-mode verdicts)
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUT_DIR = REPO_ROOT / "research/methyl_augmented_filter_phase2/cycle2"
INPUT_TSV = OUT_DIR / "data" / "cycle2_cross_sample_delta_f1.tsv"

MODES = ["transfer_fixed", "transfer_sweep", "refit"]
SAMPLE_ORDER = ["HCC1395", "H1437", "H2009", "HCC1954", "HCC1937"]


def classify_verdict(values: np.ndarray, p_wilcoxon: float) -> tuple[str, dict]:
    """Apply H_C1_5 verdict rule."""
    n_pos = int(np.sum(values > 0))
    n_neg = int(np.sum(values < 0))
    n_zero = int(np.sum(values == 0))
    n = len(values)

    # Treat near-zero (|ΔF1| < 1e-5) as zero for direction count.
    n_pos_strict = int(np.sum(values > 1e-5))
    n_neg_strict = int(np.sum(values < -1e-5))
    n_zero_strict = n - n_pos_strict - n_neg_strict

    same_dir = max(n_pos_strict, n_neg_strict)

    if same_dir >= 4 and n_pos_strict > n_neg_strict and p_wilcoxon < 0.0625:
        verdict = "POSITIVE"
    elif same_dir >= 4 and n_pos_strict > n_neg_strict:
        verdict = "DIRECTION_POSITIVE"
    elif same_dir >= 4 and n_neg_strict > n_pos_strict:
        verdict = "NEGATIVE" if p_wilcoxon < 0.0625 else "DIRECTION_NEGATIVE"
    elif same_dir <= 2:
        verdict = "MIXED_NEGATIVE"
    else:
        verdict = "MIXED"

    return verdict, {
        "n": n,
        "n_pos_strict": n_pos_strict,
        "n_neg_strict": n_neg_strict,
        "n_zero_strict": n_zero_strict,
        "majority_direction": "+" if n_pos_strict > n_neg_strict else ("-" if n_neg_strict > n_pos_strict else "0"),
        "same_direction_count": same_dir,
    }


def wilcoxon_safe(values: np.ndarray) -> dict:
    """Run Wilcoxon paired signed-rank vs 0, handling small n / zeros."""
    # Exclude exact zeros (Wilcoxon convention "wilcox" zero-method).
    out = {
        "n_used": int(np.sum(values != 0)),
        "stat": None,
        "p_value": None,
        "method": "exact",
        "note": "",
    }
    if out["n_used"] == 0:
        out["note"] = "all zeros — Wilcoxon undefined"
        return out
    try:
        stat, p = wilcoxon(values, zero_method="wilcox", alternative="two-sided", mode="exact")
        out["stat"] = float(stat)
        out["p_value"] = float(p)
    except ValueError as e:
        try:
            stat, p = wilcoxon(values, zero_method="wilcox", alternative="two-sided")
            out["stat"] = float(stat)
            out["p_value"] = float(p)
            out["method"] = "approx"
            out["note"] = f"exact failed: {e}; fallback approximate"
        except Exception as e2:
            out["note"] = f"both modes failed: {e2}"
    return out


def make_plot(per_mode: dict, out_path: Path) -> None:
    fig, axes = plt.subplots(1, len(MODES), figsize=(16, 5), sharey=False)
    for ax, mode in zip(axes, MODES):
        data = per_mode[mode]
        samples = data["samples"]
        values = data["values"]
        colors = ["#2ca02c" if v > 1e-5 else ("#d62728" if v < -1e-5 else "#7f7f7f")
                  for v in values]
        bars = ax.bar(samples, values, color=colors)
        ax.axhline(0, color="black", linewidth=0.6)
        ax.set_title(
            f"Mode: {mode}\n"
            f"n_pos={data['direction']['n_pos_strict']}/5, "
            f"n_neg={data['direction']['n_neg_strict']}/5, "
            f"Wilcoxon p={data['wilcoxon']['p_value']:.4f}"
            if data["wilcoxon"]["p_value"] is not None else
            f"Mode: {mode}\nWilcoxon: {data['wilcoxon']['note']}",
            fontsize=10,
        )
        ax.set_ylabel("ΔF1")
        ax.tick_params(axis="x", labelrotation=30)
        ax.grid(alpha=0.3, axis="y")
        for bar, v in zip(bars, values):
            ax.annotate(
                f"{v:+.4f}",
                xy=(bar.get_x() + bar.get_width() / 2, v),
                xytext=(0, 4 if v >= 0 else -10),
                textcoords="offset points",
                ha="center", va="bottom" if v >= 0 else "top",
                fontsize=8,
            )
        # Add verdict box.
        ax.text(
            0.5, 0.95,
            f"verdict: {data['verdict']}",
            transform=ax.transAxes,
            ha="center", va="top",
            fontsize=10, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                      edgecolor="black", alpha=0.85),
        )

    fig.suptitle(
        "Step B4' — H_C1_5 cross-sample Wilcoxon verdict (n=5)\n"
        "Per-sample ΔF1 vs caller F1 baseline; Wilcoxon exact min p=0.0625 (requires 5/5 same direction)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-tsv", default=str(INPUT_TSV))
    args = parser.parse_args()

    tsv = pd.read_csv(args.input_tsv, sep="\t")
    print(f"Loaded: {args.input_tsv}, rows={len(tsv)}")

    per_mode: dict = {}
    summary_rows = []
    for mode in MODES:
        sub = tsv[tsv["mode"] == mode].copy()
        # Order samples canonically.
        sub = sub.set_index("sample").reindex(SAMPLE_ORDER).reset_index()
        values = sub["delta_F1"].to_numpy(dtype=float)
        delta_p = sub["delta_precision"].to_numpy(dtype=float)
        delta_r = sub["delta_recall"].to_numpy(dtype=float)

        wilcox_F1 = wilcoxon_safe(values)
        wilcox_P = wilcoxon_safe(delta_p)
        wilcox_R = wilcoxon_safe(delta_r)
        verdict, direction = classify_verdict(values, wilcox_F1["p_value"] or 1.0)

        per_mode[mode] = {
            "samples": SAMPLE_ORDER,
            "values": values.tolist(),
            "delta_precision": delta_p.tolist(),
            "delta_recall": delta_r.tolist(),
            "wilcoxon": wilcox_F1,
            "wilcoxon_precision": wilcox_P,
            "wilcoxon_recall": wilcox_R,
            "direction": direction,
            "verdict": verdict,
        }
        summary_rows.append({
            "mode": mode,
            "n_pos": direction["n_pos_strict"],
            "n_neg": direction["n_neg_strict"],
            "n_zero": direction["n_zero_strict"],
            "majority_direction": direction["majority_direction"],
            "same_dir_count": direction["same_direction_count"],
            "wilcoxon_p_F1": wilcox_F1["p_value"],
            "wilcoxon_stat_F1": wilcox_F1["stat"],
            "wilcoxon_method_F1": wilcox_F1["method"],
            "wilcoxon_note_F1": wilcox_F1["note"],
            "wilcoxon_p_precision": wilcox_P["p_value"],
            "wilcoxon_p_recall": wilcox_R["p_value"],
            "verdict": verdict,
            "values_per_sample": " ".join(f"{s}={v:+.5f}" for s, v in zip(SAMPLE_ORDER, values)),
            "mean_delta_F1": float(np.mean(values)),
            "median_delta_F1": float(np.median(values)),
        })

        print(
            f"\n[{mode}]"
            f"  pos={direction['n_pos_strict']}  neg={direction['n_neg_strict']}"
            f"  zero={direction['n_zero_strict']}"
            f"  Wilcoxon p (F1) = {wilcox_F1['p_value']}  -> {verdict}"
        )
        for s, v in zip(SAMPLE_ORDER, values):
            print(f"    {s:>8s}  ΔF1={v:+.5f}")

    out_tsv = OUT_DIR / "data" / "cycle2_b4_wilcoxon_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(out_tsv, sep="\t", index=False)
    print(f"\nWrote summary TSV: {out_tsv}")

    fig_path = OUT_DIR / "figures" / "cycle2_b4_wilcoxon.png"
    make_plot(per_mode, fig_path)
    print(f"Wrote figure: {fig_path}")

    json_path = OUT_DIR / "intermediate" / "cycle2_b4_summary.json"
    with open(json_path, "w") as f:
        json.dump({
            "cycle": "Phase2_Cycle2_StepB4prime",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "input": args.input_tsv,
            "per_mode": per_mode,
            "summary_table": summary_rows,
        }, f, indent=2, default=str)
    print(f"Wrote JSON: {json_path}")


if __name__ == "__main__":
    main()
