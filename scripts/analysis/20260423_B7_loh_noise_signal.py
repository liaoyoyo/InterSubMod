#!/usr/bin/env python3
"""B7 analysis: LOH_Noise + Extreme AF unexpected signal (TP 88-96%).

Investigates whether LOH_Noise (legacy-derived low-confidence LOH class
but still Potential_LOH=True per HP_Ratio) behaves like LOH_Strong under
Extreme AF; if so, LOH.bed Noise tier may be overly conservative.

Primary target: HCC1395 TO (mode label "to_pileup"); baseline TP rate ~0.711.
Then cross-sample consistency on all 7 paired_full samples.

Outputs:
- Per-sample stratified TSVs
- Heatmap fig (LOH_Subtype x AF_class x CN_tier, TP rate)
- Per-sample TP rate panels for LOH_Noise vs LOH_Strong (Extreme AF)
- Summary JSON
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import select_loh_legacy

plt.rcParams["font.family"] = "DejaVu Sans"

DATA_TSV = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/data/"
    "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
)
OUT_FIG_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/"
    "figures/20260423_B7_loh_noise"
)
OUT_DATA_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/ng_kde_rescaling/data"
)

# Coverage_Category -> CN tier lumping
CN_TIER_MAP = {
    "CNV_Loss": "CN_le1",
    "Low": "CN_le1",
    "Normal": "CN_eq2",
    "Elevated": "CN_eq3",
    "CNV_Gain": "CN_eq3",
    "High_Copy": "CN_ge4",
}
CN_TIER_ORDER = ["CN_le1", "CN_eq2", "CN_eq3", "CN_ge4"]
AF_ORDER = ["Extreme", "Intermediate", "Near-half"]
SUBTYPE_ORDER = ["LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"]
PAIRED_SAMPLES = [
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
    "H2009",
    "H1437",
    "COLO829",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical unversioned LOH_Subtype input; "
            "versioned input requires LOH_Subtype_LegacyVC and an exact alias"
        ),
    )
    return parser.parse_args()


def attach_loh_legacy_view(
    frame: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    view = select_loh_legacy(frame, allow_unversioned_v1=allow_unversioned_v1)
    prepared = frame.copy()
    prepared["_loh_subtype_legacy"] = view.values
    prepared.attrs["loh_schema_contract"] = {
        "selection_field": view.field,
        "schema_status": view.schema_status,
        "allow_unversioned_v1": allow_unversioned_v1,
        "warnings": list(view.warning_messages),
    }
    return prepared


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (center - half, center + half)


def load_main(allow_unversioned_v1: bool = False) -> pd.DataFrame:
    print(f"[B7] Loading {DATA_TSV}")
    df = attach_loh_legacy_view(
        pd.read_csv(DATA_TSV, sep="\t", low_memory=False),
        allow_unversioned_v1=allow_unversioned_v1,
    )
    df["CN_tier"] = df["Coverage_Category"].map(CN_TIER_MAP).fillna("other")
    return df


def tp_rate(sub: pd.DataFrame) -> dict:
    n = len(sub)
    k = int((sub["tp_label"] == 1).sum())
    p = k / n if n > 0 else float("nan")
    lo, hi = wilson_ci(k, n)
    return {"n": n, "k": k, "tp_rate": p, "ci_lo": lo, "ci_hi": hi}


def stratified_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, mode, subtype, af_class, cn_tier), sub in df.groupby(
        ["sample", "mode", "_loh_subtype_legacy", "AF_class", "CN_tier"], dropna=False
    ):
        r = tp_rate(sub)
        r.update(
            sample=sample,
            mode=mode,
            LOH_Subtype_LegacyVC=subtype,
            AF_class=af_class,
            CN_tier=cn_tier,
        )
        rows.append(r)
    return pd.DataFrame(rows)


def hcc1395_to_extreme_af_table(df: pd.DataFrame) -> pd.DataFrame:
    sub = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")]
    rows = []
    for (subtype, af_class, cn_tier), g in sub.groupby(
        ["_loh_subtype_legacy", "AF_class", "CN_tier"], dropna=False
    ):
        r = tp_rate(g)
        r.update(LOH_Subtype_LegacyVC=subtype, AF_class=af_class, CN_tier=cn_tier)
        rows.append(r)
    return pd.DataFrame(rows)


def per_sample_subtype_af_tp(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, mode), g in df.groupby(["sample", "mode"], dropna=False):
        for (subtype, af_class), sub in g.groupby(
            ["_loh_subtype_legacy", "AF_class"], dropna=False
        ):
            r = tp_rate(sub)
            r.update(
                sample=sample,
                mode=mode,
                LOH_Subtype_LegacyVC=subtype,
                AF_class=af_class,
            )
            rows.append(r)
    return pd.DataFrame(rows)


def heatmap_hcc1395_to(tab: pd.DataFrame) -> None:
    """Heatmap: LOH_Subtype (rows) x (CN_tier x AF_class) (cols); cell = TP rate."""
    pivot_tp = tab.pivot_table(
        index="LOH_Subtype_LegacyVC",
        columns=["CN_tier", "AF_class"],
        values="tp_rate",
        aggfunc="first",
    )
    pivot_n = tab.pivot_table(
        index="LOH_Subtype_LegacyVC",
        columns=["CN_tier", "AF_class"],
        values="n",
        aggfunc="first",
    )

    # enforce ordering
    col_order = []
    for cn in CN_TIER_ORDER:
        for af in AF_ORDER:
            if (cn, af) in pivot_tp.columns:
                col_order.append((cn, af))
    pivot_tp = pivot_tp.reindex(columns=col_order)
    pivot_n = pivot_n.reindex(columns=col_order)
    pivot_tp = pivot_tp.reindex(index=[s for s in SUBTYPE_ORDER if s in pivot_tp.index])
    pivot_n = pivot_n.reindex(index=pivot_tp.index)

    fig, ax = plt.subplots(figsize=(1.0 + 0.9 * len(col_order), 3.2))
    im = ax.imshow(pivot_tp.values, cmap="RdYlGn", vmin=0.0, vmax=1.0, aspect="auto")

    # annotate
    for i in range(pivot_tp.shape[0]):
        for j in range(pivot_tp.shape[1]):
            v = pivot_tp.values[i, j]
            n = pivot_n.values[i, j]
            if np.isnan(v):
                txt = "-"
            else:
                n_int = int(n) if not np.isnan(n) else 0
                txt = f"{v:.2f}\nn={n_int}"
            ax.text(
                j, i, txt, ha="center", va="center", fontsize=7,
                color="black" if not np.isnan(v) and 0.35 < v < 0.7 else (
                    "white" if not np.isnan(v) and v < 0.35 else "black"
                ),
            )
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels([f"{cn}\n{af}" for cn, af in col_order], fontsize=7, rotation=0)
    ax.set_yticks(range(len(pivot_tp.index)))
    ax.set_yticklabels(pivot_tp.index, fontsize=9)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    cbar.set_label("TP rate", fontsize=8)
    ax.set_title(
        "B7: HCC1395 TO — TP rate by LOH_Subtype × CN_tier × AF_class\n"
        "(dashed ref: main-track baseline TP=0.711)",
        fontsize=10,
    )
    fig.tight_layout()
    out = OUT_FIG_DIR / "fig_b7_hcc1395_to_heatmap.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[B7] heatmap -> {out}")


def panel_cross_sample_extreme(per_sample: pd.DataFrame) -> None:
    """Per-sample bar plot: TP rate for LOH_Noise vs LOH_Strong at Extreme AF,
    across all 7 paired_full samples, plus HCC1395 TO for reference."""
    tgt = per_sample[
        (per_sample["AF_class"] == "Extreme")
        & (per_sample["LOH_Subtype_LegacyVC"].isin(["LOH_Noise", "LOH_Strong"]))
    ].copy()

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)

    def draw(ax, mode_label, title):
        sub = tgt[tgt["mode"] == mode_label].copy()
        samples_here = [s for s in PAIRED_SAMPLES if s in sub["sample"].unique()]
        x = np.arange(len(samples_here))
        width = 0.38
        noise_vals, strong_vals, noise_n, strong_n = [], [], [], []
        noise_lo, noise_hi, strong_lo, strong_hi = [], [], [], []
        for s in samples_here:
            row_n = sub[(sub["sample"] == s) & (sub["LOH_Subtype_LegacyVC"] == "LOH_Noise")]
            row_s = sub[(sub["sample"] == s) & (sub["LOH_Subtype_LegacyVC"] == "LOH_Strong")]
            if len(row_n) == 1:
                noise_vals.append(float(row_n["tp_rate"].iloc[0]))
                noise_n.append(int(row_n["n"].iloc[0]))
                noise_lo.append(float(row_n["ci_lo"].iloc[0]))
                noise_hi.append(float(row_n["ci_hi"].iloc[0]))
            else:
                noise_vals.append(np.nan); noise_n.append(0)
                noise_lo.append(np.nan); noise_hi.append(np.nan)
            if len(row_s) == 1:
                strong_vals.append(float(row_s["tp_rate"].iloc[0]))
                strong_n.append(int(row_s["n"].iloc[0]))
                strong_lo.append(float(row_s["ci_lo"].iloc[0]))
                strong_hi.append(float(row_s["ci_hi"].iloc[0]))
            else:
                strong_vals.append(np.nan); strong_n.append(0)
                strong_lo.append(np.nan); strong_hi.append(np.nan)

        def _safe_err(vals, los, his):
            lo_err, hi_err = [], []
            for v, lo, hi in zip(vals, los, his):
                if np.isnan(v) or np.isnan(lo) or np.isnan(hi):
                    lo_err.append(0.0); hi_err.append(0.0)
                else:
                    lo_err.append(max(v - lo, 0.0))
                    hi_err.append(max(hi - v, 0.0))
            return [lo_err, hi_err]

        noise_err = _safe_err(noise_vals, noise_lo, noise_hi)
        strong_err = _safe_err(strong_vals, strong_lo, strong_hi)
        # replace NaN bar heights with 0 for plotting
        noise_plot = [0 if np.isnan(v) else v for v in noise_vals]
        strong_plot = [0 if np.isnan(v) else v for v in strong_vals]
        ax.bar(x - width / 2, noise_plot, width, yerr=noise_err, color="#e08b4c",
               label="LOH_Noise", capsize=2)
        ax.bar(x + width / 2, strong_plot, width, yerr=strong_err, color="#3c7a59",
               label="LOH_Strong", capsize=2)
        for xi, (nv, sv, nn, sn) in enumerate(zip(noise_vals, strong_vals, noise_n, strong_n)):
            if not np.isnan(nv):
                ax.text(xi - width / 2, nv + 0.015, f"n={nn}", ha="center", fontsize=6)
            if not np.isnan(sv):
                ax.text(xi + width / 2, sv + 0.015, f"n={sn}", ha="center", fontsize=6)
        ax.set_xticks(x)
        ax.set_xticklabels(samples_here, rotation=20, ha="right", fontsize=8)
        ax.set_ylabel("TP rate (Extreme AF)")
        ax.set_ylim(0, 1.05)
        ax.axhline(0.711, color="grey", linestyle="--", linewidth=0.8, label="HCC1395 TO baseline 0.711")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(axis="y", linestyle=":", alpha=0.5)

    draw(axes[0], "paired_full", "Paired_full: LOH_Noise vs LOH_Strong @ Extreme AF")
    draw(axes[1], "to_pileup", "TO (pileup): LOH_Noise vs LOH_Strong @ Extreme AF")
    fig.tight_layout()
    out = OUT_FIG_DIR / "fig_b7_cross_sample_extreme.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[B7] cross-sample panel -> {out}")


def build_summary(
    hcc1395_to_tab: pd.DataFrame,
    per_sample: pd.DataFrame,
) -> dict:
    summary = {
        "analysis": "B7_LOH_Noise_Extreme_AF_signal",
        "date": "2026-04-23",
        "source_tsv": str(DATA_TSV),
    }

    # 1. HCC1395 TO headline: LOH_Noise vs LOH_Strong at Extreme AF (all CN)
    tgt = hcc1395_to_tab[hcc1395_to_tab["AF_class"] == "Extreme"]
    agg = (
        tgt.groupby("LOH_Subtype_LegacyVC")[["n", "k"]]
        .sum()
        .reset_index()
    )
    agg["tp_rate"] = agg["k"] / agg["n"]
    agg["ci_lo"], agg["ci_hi"] = zip(*agg.apply(lambda r: wilson_ci(int(r["k"]), int(r["n"])), axis=1))
    summary["hcc1395_to_extreme_af_pooled_cn"] = agg.to_dict(orient="records")

    # 2. HCC1395 TO LOH_Noise stratified by CN at Extreme AF
    strat = hcc1395_to_tab[
        (hcc1395_to_tab["AF_class"] == "Extreme")
        & (hcc1395_to_tab["LOH_Subtype_LegacyVC"] == "LOH_Noise")
    ]
    summary["hcc1395_to_LOH_Noise_extreme_by_CN"] = strat.to_dict(orient="records")

    # 3. cross-sample Noise vs Strong @ Extreme AF
    cross = per_sample[
        (per_sample["AF_class"] == "Extreme")
        & (per_sample["LOH_Subtype_LegacyVC"].isin(["LOH_Noise", "LOH_Strong"]))
    ]
    summary["cross_sample_extreme_af"] = cross.to_dict(orient="records")

    # 4. Cohen's h between Noise and Strong TP rate per sample×mode
    cohen_rows = []
    for (sample, mode), g in cross.groupby(["sample", "mode"]):
        rn = g[g["LOH_Subtype_LegacyVC"] == "LOH_Noise"]
        rs = g[g["LOH_Subtype_LegacyVC"] == "LOH_Strong"]
        if len(rn) != 1 or len(rs) != 1:
            continue
        pn, ps = float(rn["tp_rate"].iloc[0]), float(rs["tp_rate"].iloc[0])
        # clip
        pn_c = min(max(pn, 1e-6), 1 - 1e-6)
        ps_c = min(max(ps, 1e-6), 1 - 1e-6)
        h = 2 * (np.arcsin(np.sqrt(pn_c)) - np.arcsin(np.sqrt(ps_c)))
        cohen_rows.append(
            {
                "sample": sample,
                "mode": mode,
                "noise_tp": pn,
                "strong_tp": ps,
                "noise_n": int(rn["n"].iloc[0]),
                "strong_n": int(rs["n"].iloc[0]),
                "cohen_h_noise_minus_strong": float(h),
            }
        )
    summary["cohen_h_per_sample"] = cohen_rows

    return summary


def main() -> int:
    args = parse_args()
    OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = load_main(allow_unversioned_v1=args.allow_unversioned_v1)
    loh_schema_contract = dict(df.attrs["loh_schema_contract"])

    # full stratified table across 7 samples
    full_tab = stratified_table(df)
    full_out = OUT_DATA_DIR / "b7_loh_noise_stratified_fullsamples.tsv"
    full_tab.to_csv(full_out, sep="\t", index=False)
    print(f"[B7] full stratified -> {full_out} ({len(full_tab)} rows)")

    # HCC1395 TO focused
    hcc_to_tab = hcc1395_to_extreme_af_table(df)
    hcc_to_out = OUT_DATA_DIR / "b7_loh_noise_hcc1395_TO.tsv"
    hcc_to_tab.to_csv(hcc_to_out, sep="\t", index=False)
    print(f"[B7] HCC1395 TO table -> {hcc_to_out}")

    # per-sample aggregated (by AF only, pooled CN)
    per_sample = per_sample_subtype_af_tp(df)
    per_sample_out = OUT_DATA_DIR / "b7_loh_noise_per_sample.tsv"
    per_sample.to_csv(per_sample_out, sep="\t", index=False)
    print(f"[B7] per-sample -> {per_sample_out}")

    # Figures
    heatmap_hcc1395_to(hcc_to_tab)
    panel_cross_sample_extreme(per_sample)

    summary = build_summary(hcc_to_tab, per_sample)
    summary["loh_schema_contract"] = loh_schema_contract
    summary_out = OUT_DATA_DIR / "b7_loh_noise_summary.json"
    with open(summary_out, "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print(f"[B7] summary -> {summary_out}")

    # Quick headline print
    print("\n=== HEADLINE (HCC1395 TO @ Extreme AF, pooled CN) ===")
    for r in summary["hcc1395_to_extreme_af_pooled_cn"]:
        print(
            f"  {r['LOH_Subtype_LegacyVC']:15s}  n={int(r['n']):6d}  "
            f"TP rate={r['tp_rate']:.3f}  CI=[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]"
        )
    print("\n=== CROSS-SAMPLE Noise vs Strong TP rate @ Extreme AF ===")
    for r in summary["cohen_h_per_sample"]:
        print(
            f"  {r['sample']:18s} {r['mode']:12s} "
            f"Noise={r['noise_tp']:.3f}(n={r['noise_n']})  "
            f"Strong={r['strong_tp']:.3f}(n={r['strong_n']})  "
            f"h={r['cohen_h_noise_minus_strong']:+.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
