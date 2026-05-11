#!/usr/bin/env python3
"""
B3 | Paired-mode control for obs18 NG=2 composition analysis (H-D3).

Hypothesis (H-D3):
  In TO mode, obs18 showed Inner same-hap vs Outer cross-het TP-rate gap
  median +0.37 across 6/6 samples, explained by LOH-constrained phasing:
  germline hets within LOH are physically forced to HP+HPS (same-hap) and
  caller keeps them as candidate somatic (AF=1).

  In paired mode, the matched-normal germline caller REMOVES germline hets.
  Therefore Outer cross-het TP rate should rise and the gap should collapse
  toward zero.

  - Gap disappears  -> H-D3 PASS, strengthens LOH-constrained phasing mechanism
  - Gap preserved   -> H-D3 FAIL, mechanism needs revision

Data:
  Per-sample significance_summary.csv (paired_full, latest run per sample)
  with HPFineN_HP1 / HPFineN_HP1S / HPFineN_HP2 / HPFineN_HP2S columns.

Output:
  - figures under docs/experiments/in_progress/2026/04/figures/20260423_B3_paired_obs18/
  - TSV summaries + Wilcoxon JSON under research/tpfp_loh_af_kde_discrimination/data/
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

# Reuse obs common helpers if available
sys.path.insert(
    0,
    "/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/scripts",
)
try:
    from _obs_common import DATA_DIR, PALETTE, apply_style, research_suptitle, takeaway_caption
except Exception:
    DATA_DIR = Path(
        "/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data"
    )
    PALETTE = {
        "dark": "#1E2A44",
        "green": "#009E73",
        "red": "#D55E00",
        "blue": "#5B8DB7",
        "grey": "#8A8A8A",
    }

    def apply_style() -> None:
        plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})

    def research_suptitle(fig, title: str, context: str = "", y: float = 0.96, fontsize: int = 16) -> None:
        fig.suptitle(title + ("\n" + context if context else ""), y=y, fontsize=fontsize)

    def takeaway_caption(fig, text: str, color: str = "blue", y: float = 0.015, fontsize: int = 11) -> None:
        fig.text(0.5, y, text, ha="center", va="bottom", fontsize=fontsize, color=color, wrap=True)


FIG_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/"
    "figures/20260423_B3_paired_obs18"
)
FIG_DIR.mkdir(parents=True, exist_ok=True)

CANON_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/output/canonical")

# Latest paired_full run per sample (verified 2026-04-23 via ls -td)
SAMPLES_PAIRED: dict[str, str] = {
    "HCC1395": "20260420_HCC1395_paired_full_full_2",
    "HCC1395_DORADO": "20260420_HCC1395_DORADO_paired_full_full",
    "HCC1937": "20260421_HCC1937_paired_full_full",
    "HCC1954": "20260421_HCC1954_paired_full_full",
    "H2009": "20260421_H2009_paired_full_full",
    "H1437": "20260421_H1437_paired_full_full",
    "COLO829": "20260421_COLO829_paired_full_full",
}

BUCKET_COLS = ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]

COMBOS_ORDER = [
    "same_HP1 (HP1 + HP1-1)",
    "same_HP2 (HP2 + HP2-1)",
    "cross_het (HP1 + HP2-1)",
    "cross_het_inv (HP1-1 + HP2)",
    "other",
]

# TO mode obs18 reference
TO_OBS18_TSV = (
    "/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/"
    "data/obs18_NG2_composition_by_sample.tsv"
)


# ---------------------------------------------------------------------------
# 1. Load + categorize
# ---------------------------------------------------------------------------
def categorize_combo(has_hp1: int, has_hp1s: int, has_hp2: int, has_hp2s: int) -> str:
    if has_hp1 and has_hp1s and not has_hp2 and not has_hp2s:
        return "same_HP1 (HP1 + HP1-1)"
    if has_hp2 and has_hp2s and not has_hp1 and not has_hp1s:
        return "same_HP2 (HP2 + HP2-1)"
    if has_hp1 and has_hp2s and not has_hp1s and not has_hp2:
        return "cross_het (HP1 + HP2-1)"
    if has_hp1s and has_hp2 and not has_hp1 and not has_hp2s:
        return "cross_het_inv (HP1-1 + HP2)"
    return "other"


def load_sample(name: str, run_id: str) -> pd.DataFrame:
    base = CANON_ROOT / name / "paired_full" / run_id
    paths = {
        "tp": base / "intersubmod_tp" / "significance_summary.csv",
        "fp": base / "intersubmod_fp" / "significance_summary.csv",
    }
    frames = []
    for label, col in [("tp", 1), ("fp", 0)]:
        p = paths[label]
        if not p.exists():
            print(f"  [!] missing {p}")
            return pd.DataFrame()
        df = pd.read_csv(
            p,
            low_memory=False,
            usecols=["RegionID", "AlleleDelta", "HPFineNGroups", "Potential_LOH", "LOH_Subtype"] + BUCKET_COLS,
        )
        df["tp_label"] = col
        df["sample"] = name
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def analyze_sample(df: pd.DataFrame) -> pd.DataFrame:
    ng2 = df[df["HPFineNGroups"] == 2].copy()
    if ng2.empty:
        return pd.DataFrame()
    for c in BUCKET_COLS:
        ng2[c] = ng2[c].fillna(0).astype(int)
    ng2["combo"] = ng2.apply(
        lambda r: categorize_combo(
            int(r["HPFineN_HP1"] > 0),
            int(r["HPFineN_HP1S"] > 0),
            int(r["HPFineN_HP2"] > 0),
            int(r["HPFineN_HP2S"] > 0),
        ),
        axis=1,
    )
    ng2["loh_side"] = ng2["Potential_LOH"].apply(
        lambda x: "Inner" if str(x).lower() in ("true", "1") else "Outer"
    )
    ng2["AF_class"] = ng2["AlleleDelta"].clip(0, 1).apply(
        lambda af: "Extreme" if (af < 0.1 or af > 0.9) else ("Near-half" if 0.4 <= af <= 0.6 else "Intermediate")
    )
    ng2_ext = ng2[ng2["AF_class"] == "Extreme"].copy()
    if ng2_ext.empty:
        return pd.DataFrame()
    out = (
        ng2_ext.groupby(["sample", "loh_side", "combo"])
        .agg(n=("tp_label", "size"), n_tp=("tp_label", "sum"))
        .reset_index()
    )
    out["tp_rate"] = out["n_tp"] / out["n"]
    out["n_fp"] = out["n"] - out["n_tp"]
    return out


# ---------------------------------------------------------------------------
# 2. Plots
# ---------------------------------------------------------------------------
def plot_heatmap(all_results: pd.DataFrame) -> Path:
    samples = list(SAMPLES_PAIRED.keys())
    fig, axes = plt.subplots(1, 2, figsize=(20, 9), sharey=True)
    research_suptitle(
        fig,
        "Paired mode NG=2 組成拆分（H-D3 對照）：Inner vs Outer × TP rate",
        context="Extreme AF 子集 · cell text = TP rate (n) · grey = n<20 · paired_full, matched-normal germline caller",
        y=0.96,
        fontsize=15,
    )

    from matplotlib.colors import Normalize

    norm = Normalize(vmin=0, vmax=1)
    cmap = plt.get_cmap("RdYlGn")
    cmap.set_bad(color=PALETTE["grey"], alpha=0.4)
    last_im = None

    for ax, side in zip(axes, ["Inner", "Outer"]):
        mat = np.full((len(samples), len(COMBOS_ORDER)), np.nan)
        ns = np.zeros((len(samples), len(COMBOS_ORDER)), dtype=int)
        for i, s in enumerate(samples):
            for j, c in enumerate(COMBOS_ORDER):
                sub = all_results[
                    (all_results["sample"] == s) & (all_results["loh_side"] == side) & (all_results["combo"] == c)
                ]
                if not sub.empty:
                    n = int(sub["n"].iloc[0])
                    ns[i, j] = n
                    if n >= 20:
                        mat[i, j] = sub["tp_rate"].iloc[0]

        masked = np.ma.masked_invalid(mat)
        im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")
        last_im = im

        ax.set_xticks(range(len(COMBOS_ORDER)))
        ax.set_xticklabels([c.replace(" (", "\n(") for c in COMBOS_ORDER], rotation=0, ha="center", fontsize=9)
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples, fontsize=11)
        ax.set_title(f"{side}（LOH {'內' if side=='Inner' else '外'}）", fontsize=13, fontweight="bold")

        for i in range(len(samples)):
            for j in range(len(COMBOS_ORDER)):
                n = ns[i, j]
                rate = mat[i, j]
                if n == 0:
                    continue
                if np.isnan(rate):
                    ax.text(j, i, f"n={n}\n<20", ha="center", va="center", fontsize=7, color=PALETTE["dark"])
                else:
                    colour = "black" if rate > 0.55 else "white"
                    ax.text(j, i, f"{rate*100:.0f}%\nn={n:,}", ha="center", va="center", fontsize=8, color=colour)

    cbar_ax = fig.add_axes([0.93, 0.12, 0.012, 0.75])
    if last_im is not None:
        cbar = fig.colorbar(last_im, cax=cbar_ax)
        cbar.set_label("TP rate", fontsize=10)

    # Takeaway
    inner_samehap = (
        all_results[
            (all_results["loh_side"] == "Inner")
            & (all_results["combo"].str.startswith("same_"))
            & (all_results["n"] >= 50)
        ]
        .groupby("sample")
        .agg(avg_tp=("tp_rate", "mean"))
    )
    outer_crosshet = (
        all_results[
            (all_results["loh_side"] == "Outer")
            & (all_results["combo"].str.startswith("cross_"))
            & (all_results["n"] >= 50)
        ]
        .groupby("sample")
        .agg(avg_tp=("tp_rate", "mean"))
    )
    inner_med = inner_samehap["avg_tp"].median() if len(inner_samehap) else np.nan
    outer_med = outer_crosshet["avg_tp"].median() if len(outer_crosshet) else np.nan
    gap = inner_med - outer_med

    takeaway = (
        f"Paired mode: Inner × same-hap TP median = {inner_med:.2f}（{len(inner_samehap)} samples）· "
        f"Outer × cross-het TP median = {outer_med:.2f}（{len(outer_crosshet)} samples）· "
        f"gap = {gap:+.2f}."
    )
    if abs(gap) < 0.10:
        takeaway += " ✓ Gap 在 paired 模式下大幅縮小（|gap|<0.10）→ H-D3 PASS：germline caller 移除 cross-het。"
        colour = "green"
    elif gap > 0.20:
        takeaway += " ⚠ Gap 在 paired 模式仍 ≥0.20 → H-D3 FAIL：機制非純 LOH-phasing，需檢視其他來源。"
        colour = "red"
    else:
        takeaway += " ◐ Gap 部分縮小（0.10-0.20）→ 機制部分解釋但仍有殘留。"
        colour = "blue"

    takeaway_caption(fig, takeaway, color=colour, y=0.015, fontsize=11)
    fig.tight_layout(rect=[0, 0.05, 0.92, 0.92])
    out = FIG_DIR / "B3_paired_NG2_composition_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_gap_comparison(paired_res: pd.DataFrame, to_res: pd.DataFrame) -> tuple[Path, pd.DataFrame, dict]:
    """Side-by-side bar chart: TO gap vs Paired gap per sample."""
    def per_sample_gap(df: pd.DataFrame) -> pd.DataFrame:
        inner_same = (
            df[(df["loh_side"] == "Inner") & (df["combo"].str.startswith("same_")) & (df["n"] >= 50)]
            .groupby("sample")
            .agg(inner_same_tp=("tp_rate", "mean"), inner_same_n=("n", "sum"))
        )
        outer_cross = (
            df[(df["loh_side"] == "Outer") & (df["combo"].str.startswith("cross_")) & (df["n"] >= 50)]
            .groupby("sample")
            .agg(outer_cross_tp=("tp_rate", "mean"), outer_cross_n=("n", "sum"))
        )
        joined = inner_same.join(outer_cross, how="outer")
        joined["gap"] = joined["inner_same_tp"] - joined["outer_cross_tp"]
        return joined.reset_index()

    paired_gap = per_sample_gap(paired_res).rename(columns={"gap": "paired_gap"})
    to_gap = per_sample_gap(to_res).rename(columns={"gap": "to_gap"})
    merged = paired_gap.merge(to_gap[["sample", "to_gap"]], on="sample", how="outer")

    # Wilcoxon on paired gaps across samples (H0: paired_gap == 0)
    paired_vals = merged["paired_gap"].dropna().values
    wil_stat: dict = {}
    if len(paired_vals) >= 5:
        try:
            stat, p = wilcoxon(paired_vals, alternative="two-sided", zero_method="wilcox")
            wil_stat["paired_gap_wilcoxon_stat"] = float(stat)
            wil_stat["paired_gap_wilcoxon_p"] = float(p)
        except Exception as e:
            wil_stat["paired_gap_wilcoxon_error"] = str(e)
    wil_stat["paired_gap_median"] = float(np.nanmedian(merged["paired_gap"]))
    wil_stat["to_gap_median"] = float(np.nanmedian(merged["to_gap"]))
    wil_stat["paired_gap_values"] = {r["sample"]: None if pd.isna(r["paired_gap"]) else float(r["paired_gap"]) for _, r in merged.iterrows()}
    wil_stat["to_gap_values"] = {r["sample"]: None if pd.isna(r["to_gap"]) else float(r["to_gap"]) for _, r in merged.iterrows()}

    # Paired vs TO difference (per sample with both values)
    both = merged.dropna(subset=["paired_gap", "to_gap"])
    if len(both) >= 5:
        try:
            stat2, p2 = wilcoxon(both["to_gap"].values, both["paired_gap"].values, alternative="two-sided", zero_method="wilcox")
            wil_stat["to_vs_paired_wilcoxon_stat"] = float(stat2)
            wil_stat["to_vs_paired_wilcoxon_p"] = float(p2)
        except Exception as e:
            wil_stat["to_vs_paired_wilcoxon_error"] = str(e)

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    samples = merged["sample"].tolist()
    x = np.arange(len(samples))
    width = 0.38
    to_vals = merged["to_gap"].fillna(0).values
    paired_vals_plot = merged["paired_gap"].fillna(0).values
    bars_to = ax.bar(x - width / 2, to_vals, width, label="TO mode (obs18 baseline)", color=PALETTE["red"], edgecolor=PALETTE["dark"])
    bars_paired = ax.bar(x + width / 2, paired_vals_plot, width, label="Paired mode (this B3)", color=PALETTE["green"], edgecolor=PALETTE["dark"])

    for i, (to_v, pa_v) in enumerate(zip(merged["to_gap"].values, merged["paired_gap"].values)):
        if not pd.isna(to_v):
            ax.text(i - width / 2, to_v + 0.01, f"{to_v:+.2f}", ha="center", va="bottom", fontsize=9)
        else:
            ax.text(i - width / 2, 0.02, "N/A", ha="center", va="bottom", fontsize=9, color=PALETTE["grey"])
        if not pd.isna(pa_v):
            ax.text(i + width / 2, pa_v + 0.01 if pa_v >= 0 else pa_v - 0.03, f"{pa_v:+.2f}", ha="center", va="bottom" if pa_v >= 0 else "top", fontsize=9)
        else:
            ax.text(i + width / 2, 0.02, "N/A", ha="center", va="bottom", fontsize=9, color=PALETTE["grey"])

    ax.axhline(0, color=PALETTE["dark"], linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=20, ha="right")
    ax.set_ylabel("Inner same-hap − Outer cross-het TP rate gap")
    ax.set_title("H-D3 對照：TO gap vs Paired gap per sample（Extreme AF, NG=2）", fontsize=13, fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    # Annotate stats
    stat_lines = [f"Paired gap median = {wil_stat['paired_gap_median']:+.3f}", f"TO gap median = {wil_stat['to_gap_median']:+.3f}"]
    if "paired_gap_wilcoxon_p" in wil_stat:
        stat_lines.append(f"Wilcoxon paired_gap vs 0: p = {wil_stat['paired_gap_wilcoxon_p']:.3g}")
    if "to_vs_paired_wilcoxon_p" in wil_stat:
        stat_lines.append(f"Wilcoxon TO vs Paired (paired): p = {wil_stat['to_vs_paired_wilcoxon_p']:.3g}")
    ax.text(
        0.02,
        0.98,
        "\n".join(stat_lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor=PALETTE["dark"], alpha=0.9),
    )

    fig.tight_layout()
    out = FIG_DIR / "B3_paired_vs_TO_gap_per_sample.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out, merged, wil_stat


def plot_proportion_stacked(all_results: pd.DataFrame) -> Path:
    samples = list(SAMPLES_PAIRED.keys())
    combo_colors = {
        "same_HP1 (HP1 + HP1-1)": "#2E7D5B",
        "same_HP2 (HP2 + HP2-1)": "#5BA889",
        "cross_het (HP1 + HP2-1)": "#C26B3A",
        "cross_het_inv (HP1-1 + HP2)": "#E19662",
        "other": "#8A8A8A",
    }
    fig, axes = plt.subplots(1, 2, figsize=(20, 9), sharey=True)
    research_suptitle(
        fig,
        "Paired mode NG=2 組成 Inner vs Outer 占比（H-D3 對照）",
        context="stacked bar：每樣本 NG=2 cells 的組成比例 · 綠系 = 同 hap · 橘系 = 跨 hap",
        y=0.96,
        fontsize=15,
    )

    for ax, side in zip(axes, ["Inner", "Outer"]):
        bottoms = np.zeros(len(samples))
        for combo in COMBOS_ORDER:
            heights = []
            for s in samples:
                sub = all_results[(all_results["sample"] == s) & (all_results["loh_side"] == side)]
                total = sub["n"].sum() if len(sub) else 0
                cnt = sub[sub["combo"] == combo]["n"].sum() if len(sub) else 0
                heights.append(cnt / total if total > 0 else 0)
            heights = np.array(heights)
            ax.bar(
                range(len(samples)),
                heights,
                bottom=bottoms,
                color=combo_colors[combo],
                edgecolor=PALETTE["dark"],
                linewidth=0.3,
                label=combo if side == "Inner" else None,
            )
            for i, (h, b) in enumerate(zip(heights, bottoms)):
                if h >= 0.05:
                    ax.text(i, b + h / 2, f"{h*100:.0f}%", ha="center", va="center", fontsize=8, color="white", fontweight="bold")
            bottoms += heights

        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=25, ha="right", fontsize=10)
        ax.set_ylabel("proportion of NG=2 cells", fontsize=10)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"{side}（LOH {'內' if side=='Inner' else '外'}）", fontsize=13, fontweight="bold")
        ax.grid(True, axis="y", alpha=0.25)

    axes[0].legend(loc="upper center", bbox_to_anchor=(1.05, -0.08), ncol=5, fontsize=9, frameon=False)
    fig.tight_layout(rect=[0, 0.05, 1, 0.92])
    out = FIG_DIR / "B3_paired_NG2_composition_proportion.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# 3. Main
# ---------------------------------------------------------------------------
def main() -> None:
    apply_style()

    # TO reference
    to_df = pd.read_csv(TO_OBS18_TSV, sep="\t")

    # Load paired samples
    all_results = []
    for name, run_id in SAMPLES_PAIRED.items():
        print(f"[B3] loading paired {name} ({run_id})")
        df = load_sample(name, run_id)
        if df.empty:
            continue
        res = analyze_sample(df)
        if not res.empty:
            all_results.append(res)

    if not all_results:
        print("[B3] no paired data loaded")
        return

    combined = pd.concat(all_results, ignore_index=True)
    out_tsv = DATA_DIR / "B3_paired_NG2_composition_by_sample.tsv"
    combined.to_csv(out_tsv, sep="\t", index=False)
    print(f"[B3] wrote {out_tsv}  ({len(combined)} rows)")

    print("\n=== Paired NG=2 組成 × Inner/Outer × TP rate（7 samples, Extreme AF） ===")
    pivot = combined.pivot_table(
        index=["sample", "loh_side"], columns="combo", values="tp_rate", aggfunc="first"
    ).round(3)
    print(pivot.to_string())

    out1 = plot_heatmap(combined)
    print(f"\n[B3] wrote {out1}")

    out2 = plot_proportion_stacked(combined)
    print(f"[B3] wrote {out2}")

    out3, merged_gap, wil = plot_gap_comparison(combined, to_df)
    print(f"[B3] wrote {out3}")

    gap_tsv = DATA_DIR / "B3_paired_vs_TO_gap_per_sample.tsv"
    merged_gap.to_csv(gap_tsv, sep="\t", index=False)
    print(f"[B3] wrote {gap_tsv}")

    wil_json = DATA_DIR / "B3_wilcoxon_gap_stats.json"
    wil_json.write_text(json.dumps(wil, indent=2, default=str))
    print(f"[B3] wrote {wil_json}")

    print("\n=== Per-sample gap comparison ===")
    print(merged_gap.round(3).to_string(index=False))
    print("\n=== Wilcoxon stats ===")
    print(json.dumps(wil, indent=2, default=str))


if __name__ == "__main__":
    main()
