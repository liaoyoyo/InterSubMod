#!/usr/bin/env python3
"""
obs18 | NG=2 composition decomposition across 6 TO samples.

Confirms the corrected mechanism: at NG=2, the Inner-vs-Outer TP rate gap
is explained by WHICH 2 buckets are populated among {HP1, HP1-1, HP2, HP2-1}.

Same-haplotype (HP1+HP1-1 or HP2+HP2-1) should dominate Inner (LOH regions)
and have high TP rate. Cross-haplotype (HP1+HP2-1 or HP1-1+HP2) should dominate
Outer and have ambiguous TP (germline het looks similar to somatic het).

Data source: per-sample significance_summary.csv (has HPFineN_HP1/S/HP2/S columns).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (
    DATA_DIR, FIG_NEW_DIR, PALETTE,
    apply_style, research_suptitle, sample_label, takeaway_caption,
)

OUT_DIR = FIG_NEW_DIR / "obs18"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Each TO sample and its significance_summary.csv paths
SAMPLES_TO = {
    "HCC1395": {
        "tp": "/tmp/ism_hp_fix_phase1/tp_off/significance_summary.csv",
        "fp": "/tmp/ism_hp_fix_phase1/fp_off/significance_summary.csv",
    },
    "HCC1395_DORADO": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "HCC1937": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "HCC1954": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "H2009": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
    "H1437": {
        "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod/intersubmod_tp/significance_summary.csv",
        "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step05_intersubmod/intersubmod_fp/significance_summary.csv",
    },
}

BUCKET_COLS = ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]


def categorize_combo(has_hp1: int, has_hp1s: int, has_hp2: int, has_hp2s: int) -> str:
    """Classify NG=2 composition into 4 scientific categories."""
    # Same-haplotype split (variant within single hap)
    if has_hp1 and has_hp1s and not has_hp2 and not has_hp2s:
        return "same_HP1 (HP1 + HP1-1)"
    if has_hp2 and has_hp2s and not has_hp1 and not has_hp1s:
        return "same_HP2 (HP2 + HP2-1)"
    # Cross-haplotype (canonical het phasing)
    if has_hp1 and has_hp2s and not has_hp1s and not has_hp2:
        return "cross_het (HP1 + HP2-1)"
    if has_hp1s and has_hp2 and not has_hp1 and not has_hp2s:
        return "cross_het_inv (HP1-1 + HP2)"
    # Other rare
    return "other"


def load_sample(name: str, paths: dict) -> pd.DataFrame:
    frames = []
    for label, col in [("tp", 1), ("fp", 0)]:
        try:
            df = pd.read_csv(paths[label], low_memory=False,
                             usecols=["RegionID", "AlleleDelta", "HPFineNGroups",
                                      "Potential_LOH", "LOH_Subtype"] + BUCKET_COLS)
            df["tp_label"] = col
            df["sample"] = name
            frames.append(df)
        except Exception as e:
            print(f"  [!] {name} {label}: {e}")
            return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def analyze_sample(df: pd.DataFrame) -> pd.DataFrame:
    """For NG=2 rows, categorize composition and group by loh_side."""
    ng2 = df[df["HPFineNGroups"] == 2].copy()
    if len(ng2) == 0:
        return pd.DataFrame()

    # Fill NaN with 0
    for c in BUCKET_COLS:
        ng2[c] = ng2[c].fillna(0).astype(int)

    ng2["combo"] = ng2.apply(
        lambda r: categorize_combo(
            int(r["HPFineN_HP1"] > 0), int(r["HPFineN_HP1S"] > 0),
            int(r["HPFineN_HP2"] > 0), int(r["HPFineN_HP2S"] > 0),
        ),
        axis=1,
    )
    ng2["loh_side"] = ng2["Potential_LOH"].apply(
        lambda x: "Inner" if str(x).lower() in ("true", "1") else "Outer"
    )
    # Filter Extreme AF only for consistency with obs17
    ng2["AF_class"] = ng2["AlleleDelta"].clip(0, 1).apply(
        lambda af: "Extreme" if (af < 0.1 or af > 0.9) else
        ("Near-half" if 0.4 <= af <= 0.6 else "Intermediate")
    )
    ng2_ext = ng2[ng2["AF_class"] == "Extreme"].copy()

    result = (ng2_ext.groupby(["sample", "loh_side", "combo"])
              .agg(n=("tp_label", "size"), n_tp=("tp_label", "sum"))
              .reset_index())
    result["tp_rate"] = result["n_tp"] / result["n"]
    result["n_fp"] = result["n"] - result["n_tp"]
    return result


def plot_composition_heatmap(all_results: pd.DataFrame) -> Path:
    """6 samples × composition × Inner/Outer heatmap."""
    combos_order = [
        "same_HP1 (HP1 + HP1-1)",
        "same_HP2 (HP2 + HP2-1)",
        "cross_het (HP1 + HP2-1)",
        "cross_het_inv (HP1-1 + HP2)",
        "other",
    ]
    samples = list(SAMPLES_TO.keys())

    fig, axes = plt.subplots(1, 2, figsize=(20, 9), sharey=True)
    research_suptitle(
        fig,
        "NG=2 組成拆分：Inner 由「同 haplotype 內分裂」支配，Outer 由「跨 haplotype」支配？",
        context="Extreme AF 子集 · cell text = TP rate (n) · grey = n<20 · 兩張圖左=Inner 右=Outer",
        y=0.96, fontsize=16,
    )

    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=1)
    cmap = plt.get_cmap("RdYlGn")
    cmap.set_bad(color=PALETTE["grey"], alpha=0.4)
    last_im = None

    for ax, side in zip(axes, ["Inner", "Outer"]):
        mat = np.full((len(samples), len(combos_order)), np.nan)
        ns = np.zeros((len(samples), len(combos_order)), dtype=int)
        for i, s in enumerate(samples):
            for j, c in enumerate(combos_order):
                sub = all_results[(all_results["sample"] == s) &
                                  (all_results["loh_side"] == side) &
                                  (all_results["combo"] == c)]
                if not sub.empty:
                    n = int(sub["n"].iloc[0])
                    ns[i, j] = n
                    if n >= 20:
                        mat[i, j] = sub["tp_rate"].iloc[0]

        masked = np.ma.masked_invalid(mat)
        im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")
        last_im = im

        ax.set_xticks(range(len(combos_order)))
        ax.set_xticklabels([c.replace(" (", "\n(") for c in combos_order],
                           rotation=0, ha="center", fontsize=9)
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples, fontsize=11)
        ax.set_title(f"{side}（LOH {'內' if side=='Inner' else '外'}）",
                     fontsize=13, fontweight="bold", color=PALETTE["dark"])

        for i in range(len(samples)):
            for j in range(len(combos_order)):
                n = ns[i, j]
                rate = mat[i, j]
                if n == 0:
                    continue
                if np.isnan(rate):
                    ax.text(j, i, f"n={n}\n<20", ha="center", va="center",
                            fontsize=7, color=PALETTE["dark"])
                else:
                    colour = "black" if rate > 0.55 else "white"
                    ax.text(j, i, f"{rate*100:.0f}%\nn={n:,}",
                            ha="center", va="center", fontsize=8, color=colour)

    cbar_ax = fig.add_axes([0.93, 0.12, 0.012, 0.75])
    if last_im is not None:
        cbar = fig.colorbar(last_im, cax=cbar_ax)
        cbar.set_label("TP rate", fontsize=10)

    # Compute takeaway
    # Aggregate Inner same-hap TP rate vs Outer cross-het TP rate
    inner_samehap = all_results[
        (all_results["loh_side"] == "Inner") &
        (all_results["combo"].str.startswith("same_")) &
        (all_results["n"] >= 50)
    ].groupby("sample").agg(avg_tp=("tp_rate", "mean"))
    outer_crosshet = all_results[
        (all_results["loh_side"] == "Outer") &
        (all_results["combo"].str.startswith("cross_")) &
        (all_results["n"] >= 50)
    ].groupby("sample").agg(avg_tp=("tp_rate", "mean"))

    inner_med = inner_samehap["avg_tp"].median() if len(inner_samehap) else np.nan
    outer_med = outer_crosshet["avg_tp"].median() if len(outer_crosshet) else np.nan

    takeaway = (
        f"Inner × same-haplotype 組合 TP rate median = {inner_med:.2f}（跨 {len(inner_samehap)} 樣本）· "
        f"Outer × cross-haplotype 組合 TP rate median = {outer_med:.2f}（跨 {len(outer_crosshet)} 樣本）· "
        f"gap = {inner_med - outer_med:+.2f}。"
    )
    if inner_med - outer_med > 0.15:
        takeaway += " ✓ 同 hap Inner >> 跨 hap Outer，確認 LOH-constrained phasing 機制。"
        colour = "green"
    else:
        takeaway += " ◐ gap 較小，可能有其他 confound。"
        colour = "blue"

    takeaway_caption(fig, takeaway, color=colour, y=0.015, fontsize=11)

    fig.tight_layout(rect=[0, 0.05, 0.92, 0.92])
    out = OUT_DIR / "obs18_NG2_composition_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_proportion_stacked(all_results: pd.DataFrame) -> Path:
    """For each sample × side, stacked bar showing which compositions dominate."""
    samples = list(SAMPLES_TO.keys())
    combos_order = ["same_HP1 (HP1 + HP1-1)", "same_HP2 (HP2 + HP2-1)",
                    "cross_het (HP1 + HP2-1)", "cross_het_inv (HP1-1 + HP2)", "other"]
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
        "NG=2 組成在 Inner vs Outer 的**占比**差異（Extreme AF）",
        context="stacked bar：每樣本 NG=2 cells 的組成比例 · 綠系 = 同 hap · 橘系 = 跨 hap",
        y=0.96, fontsize=16,
    )

    for ax, side in zip(axes, ["Inner", "Outer"]):
        bottoms = np.zeros(len(samples))
        for combo in combos_order:
            heights = []
            for s in samples:
                sub = all_results[(all_results["sample"] == s) &
                                  (all_results["loh_side"] == side)]
                total = sub["n"].sum() if len(sub) else 0
                cnt = sub[sub["combo"] == combo]["n"].sum() if len(sub) else 0
                heights.append(cnt / total if total > 0 else 0)
            heights = np.array(heights)
            ax.bar(range(len(samples)), heights, bottom=bottoms,
                   color=combo_colors[combo], edgecolor=PALETTE["dark"],
                   linewidth=0.3, label=combo if side == "Inner" else None)
            # Annotate proportions >5%
            for i, (h, b) in enumerate(zip(heights, bottoms)):
                if h >= 0.05:
                    ax.text(i, b + h / 2, f"{h*100:.0f}%", ha="center", va="center",
                            fontsize=8, color="white", fontweight="bold")
            bottoms += heights

        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=25, ha="right", fontsize=10)
        ax.set_ylabel("proportion of NG=2 cells", fontsize=10)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"{side}（LOH {'內' if side=='Inner' else '外'}）",
                     fontsize=13, fontweight="bold", color=PALETTE["dark"])
        ax.grid(True, axis="y", alpha=0.25)

    # Single legend
    axes[0].legend(loc="upper center", bbox_to_anchor=(1.05, -0.08),
                   ncol=5, fontsize=9, frameon=False)

    # Takeaway
    rows = []
    for s in samples:
        for side in ["Inner", "Outer"]:
            sub = all_results[(all_results["sample"] == s) &
                              (all_results["loh_side"] == side)]
            total = sub["n"].sum()
            same = sub[sub["combo"].str.startswith("same_")]["n"].sum()
            cross = sub[sub["combo"].str.startswith("cross_")]["n"].sum()
            rows.append({
                "sample": s, "side": side, "total": total,
                "same_pct": same / total if total else 0,
                "cross_pct": cross / total if total else 0,
            })
    comp_df = pd.DataFrame(rows)
    inner_same_med = comp_df[comp_df["side"] == "Inner"]["same_pct"].median()
    outer_cross_med = comp_df[comp_df["side"] == "Outer"]["cross_pct"].median()

    takeaway = (f"Inner NG=2 中「same-hap」占比 median = {inner_same_med*100:.0f}%（跨 6 樣本）· "
                f"Outer NG=2 中「cross-het」占比 median = {outer_cross_med*100:.0f}%。")
    if inner_same_med > 0.5 and outer_cross_med > 0.5:
        takeaway += " ✓ 組成在 Inner / Outer 幾乎完全分離 → 確認 LOH 區物理上強迫 same-hap phasing。"
        colour = "green"
    elif inner_same_med > 0.35:
        takeaway += " ◐ 主趨勢存在但不完全分離。"
        colour = "blue"
    else:
        takeaway += " ⚠ 趨勢不顯著，需檢查 phasing tool 行為。"
        colour = "red"

    takeaway_caption(fig, takeaway, color=colour, y=0.015, fontsize=11)

    fig.tight_layout(rect=[0, 0.05, 1, 0.92])
    out = OUT_DIR / "obs18_NG2_composition_proportion.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out, comp_df


def main() -> None:
    apply_style()
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []
    for name, paths in SAMPLES_TO.items():
        print(f"[obs18] loading {name}")
        df = load_sample(name, paths)
        if df.empty:
            continue
        res = analyze_sample(df)
        if not res.empty:
            all_results.append(res)

    if not all_results:
        print("[obs18] no data")
        return

    combined = pd.concat(all_results, ignore_index=True)
    out_tsv = DATA_DIR / "obs18_NG2_composition_by_sample.tsv"
    combined.to_csv(out_tsv, sep="\t", index=False)
    print(f"[obs18] wrote {out_tsv}  ({len(combined)} rows)")

    print("\n=== NG=2 組成 × Inner/Outer × TP rate（6 TO samples, Extreme AF） ===")
    pivot = combined.pivot_table(
        index=["sample", "loh_side"], columns="combo",
        values="tp_rate", aggfunc="first",
    ).round(3)
    print(pivot.to_string())

    out1 = plot_composition_heatmap(combined)
    print(f"\n[obs18] wrote {out1}")
    out2, comp_df = plot_proportion_stacked(combined)
    print(f"[obs18] wrote {out2}")

    out_tsv2 = DATA_DIR / "obs18_NG2_composition_proportion.tsv"
    comp_df.to_csv(out_tsv2, sep="\t", index=False)
    print(f"[obs18] wrote {out_tsv2}")

    print("\n=== Composition proportion summary ===")
    print(comp_df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
