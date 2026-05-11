#!/usr/bin/env python3
"""
obs17 | TO mode 4D analysis: sample × AF_class × loh_side × HPFineNGroups.

Goal: 觀察 TO 狀況下 sample × AF × LOH 內外 × 甲基分群 的特性與異常。

Outputs:
  - data/obs17_TO_4d_cube.tsv            — long-form per-cell TP rate + n + FP
  - data/obs17_TO_direction_summary.tsv  — per-sample direction vote table
  - figures/new/obs17/obs17_TO_afxng_heatmap_by_lohside.png
  - figures/new/obs17/obs17_TO_inner_vs_outer_ng_gradient.png
  - figures/new/obs17/obs17_TO_anomalies.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (
    DATA_DIR, DIM_ORDERS, FIG_NEW_DIR, PALETTE,
    SAMPLE_MODE_TO, apply_style, build_5d_key,
    compute_tp_rate_by_key, load_master, research_suptitle,
    sample_label, takeaway_caption, annotate_power_tier,
)

OUT_DIR = FIG_NEW_DIR / "obs17"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def prepare_to_df(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to TO mode, build 5D key, drop NG=0 noise."""
    to = df[df["mode"] == "to_pileup"].copy()
    to = build_5d_key(to)
    to = to[to["HPFineNGroups"].isin(DIM_ORDERS["HPFineNGroups"])]
    # Drop Unknown loh_side
    to = to[to["loh_side"].isin(["Inner", "Outer"])]
    return to


def build_4d_cube(to: pd.DataFrame) -> pd.DataFrame:
    """Per-cell 4D cube: sample × AF_class × loh_side × HPFineNGroups."""
    rows = []
    for sample in [s for s, _ in SAMPLE_MODE_TO]:
        sub = to[to["sample"] == sample]
        if len(sub) == 0:
            continue
        agg = compute_tp_rate_by_key(
            sub,
            ["AF_class", "loh_side", "HPFineNGroups"],
            min_n=20,
        )
        agg["sample"] = sample
        agg["baseline_tp"] = sub["tp_label"].mean()
        rows.append(agg)
    out = pd.concat(rows, ignore_index=True)
    # Reorder columns
    out = out[["sample", "baseline_tp", "AF_class", "loh_side", "HPFineNGroups",
               "n", "n_tp", "tp_rate", "wilson_lo", "wilson_hi"]]
    out["n_fp"] = out["n"] - out["n_tp"]
    return out


def plot_afxng_heatmap_by_lohside(cube: pd.DataFrame) -> Path:
    """6 samples × 2 loh_side = 12 heatmap panels, AF × NG colored by TP rate.

    Layout: 4 rows × 3 cols where each 2 consecutive rows belong to one sample
    (top = Inner, bottom = Outer) — but this is hard to read. Instead do
    6 rows × 2 cols (one row per sample, cols = Inner/Outer).
    """
    samples = [s for s, _ in SAMPLE_MODE_TO]
    af_order = DIM_ORDERS["AF_class"]
    ng_order = DIM_ORDERS["HPFineNGroups"]

    fig, axes = plt.subplots(6, 2, figsize=(11, 22), sharex=True, sharey=True)
    research_suptitle(
        fig,
        "TO mode：AF × HPFineNGroups (NG) 跨樣本 × LOH 內外對比",
        context="每行 = 一個 TO sample；左欄 = Inner (LOH 內部)，右欄 = Outer；"
                "cell color = TP rate；cell text = TP% / n；grey = n<20",
        y=0.985, fontsize=16,
    )

    norm = Normalize(vmin=0.0, vmax=1.0)
    cmap = plt.get_cmap("RdYlGn")
    cmap.set_bad(color=PALETTE["grey"], alpha=0.4)
    last_im = None

    for row_i, sample in enumerate(samples):
        baseline = cube[cube["sample"] == sample]["baseline_tp"].iloc[0] \
            if (cube["sample"] == sample).any() else np.nan
        for col_i, side in enumerate(["Inner", "Outer"]):
            ax = axes[row_i, col_i]
            sub = cube[(cube["sample"] == sample) & (cube["loh_side"] == side)]
            pivot_rate = sub.pivot_table(
                index="AF_class", columns="HPFineNGroups", values="tp_rate"
            ).reindex(index=af_order, columns=ng_order)
            pivot_n = sub.pivot_table(
                index="AF_class", columns="HPFineNGroups", values="n", fill_value=0
            ).reindex(index=af_order, columns=ng_order).fillna(0).astype(int)

            data = pivot_rate.values.astype(float)
            masked = np.ma.masked_invalid(data)
            im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")
            last_im = im

            ax.set_xticks(range(len(ng_order)))
            ax.set_xticklabels([f"NG={n}" for n in ng_order], fontsize=9)
            ax.set_yticks(range(len(af_order)))
            ax.set_yticklabels(af_order, fontsize=9)

            panel_title = f"{sample} · {side}" + (f"  base={baseline:.2f}" if pd.notna(baseline) else "")
            ax.set_title(panel_title, fontsize=11, fontweight="bold", color=PALETTE["dark"])

            for yi in range(len(af_order)):
                for xi in range(len(ng_order)):
                    n = int(pivot_n.iat[yi, xi])
                    rate = pivot_rate.iat[yi, xi]
                    if n == 0:
                        continue
                    if np.isnan(rate):
                        ax.text(xi, yi, f"n={n}\n<20", ha="center", va="center",
                                fontsize=7, color=PALETTE["dark"])
                    else:
                        colour = "black" if rate > 0.55 else "white"
                        ax.text(xi, yi, f"{rate*100:.0f}%\nn={n:,}",
                                ha="center", va="center", fontsize=8, color=colour)

    # shared colorbar
    cbar_ax = fig.add_axes([0.93, 0.12, 0.012, 0.78])
    if last_im is not None:
        cbar = fig.colorbar(last_im, cax=cbar_ax)
        cbar.set_label("TP rate", fontsize=10)

    # Compute takeaway
    # Direction: for each sample × side, does TP rate go up with NG within a given AF class?
    up_count, down_count, flat_count = 0, 0, 0
    for sample in samples:
        for side in ["Inner", "Outer"]:
            sub = cube[(cube["sample"] == sample) & (cube["loh_side"] == side) &
                       (cube["AF_class"] == "Extreme")]
            rates = sub.set_index("HPFineNGroups")["tp_rate"].reindex(ng_order).dropna()
            if len(rates) < 2:
                continue
            gap = rates.iloc[-1] - rates.iloc[0]
            if gap > 0.05: up_count += 1
            elif gap < -0.05: down_count += 1
            else: flat_count += 1

    total = up_count + down_count + flat_count
    takeaway = (f"NG 梯度方向（Extreme AF 內部 slice）：↑{up_count} / ↓{down_count} / "
                f"○{flat_count}（共 {total} 個 sample×side）；"
                f"Inner-Outer 分離的必要性見右側對比。")
    colour = "green" if up_count >= total * 0.6 else "blue"
    takeaway_caption(fig, takeaway, color=colour, y=0.013, fontsize=11)

    fig.tight_layout(rect=[0, 0.04, 0.92, 0.97])
    out = OUT_DIR / "obs17_TO_afxng_heatmap_by_lohside.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_inner_vs_outer_ng_gradient(cube: pd.DataFrame) -> Path:
    """For each sample, plot Inner NG gradient vs Outer NG gradient on Extreme AF."""
    samples = [s for s, _ in SAMPLE_MODE_TO]
    ng_order = DIM_ORDERS["HPFineNGroups"]

    fig, axes = plt.subplots(2, 3, figsize=(18, 11), sharey=True)
    research_suptitle(
        fig,
        "TO mode：Inner vs Outer 下 NG 梯度跨樣本比較（Extreme AF 子集）",
        context="藍 = Inner (LOH 內部)；橘 = Outer；實線 = TP rate；error bar = Wilson 95% CI；"
                "虛線 = sample baseline；只看 Extreme AF（占 TO 95%）",
        y=0.975, fontsize=17,
    )
    axes_flat = list(axes.flatten())

    for ax, sample in zip(axes_flat, samples):
        baseline = cube[cube["sample"] == sample]["baseline_tp"].iloc[0]
        sub = cube[(cube["sample"] == sample) & (cube["AF_class"] == "Extreme")]
        for side, color in [("Inner", PALETTE["blue"]), ("Outer", PALETTE["accent"])]:
            side_sub = sub[sub["loh_side"] == side].set_index("HPFineNGroups").reindex(ng_order)
            xs = np.arange(len(ng_order))
            ys = side_sub["tp_rate"].values
            lo = side_sub["wilson_lo"].values
            hi = side_sub["wilson_hi"].values
            valid = ~np.isnan(ys)
            ax.plot(xs[valid], ys[valid], "o-", color=color, label=side, linewidth=2,
                    markersize=7)
            ax.errorbar(xs[valid], ys[valid],
                        yerr=[ys[valid] - lo[valid], hi[valid] - ys[valid]],
                        fmt="none", ecolor=color, capsize=4, alpha=0.6)
            # annotate n
            for i, ng in enumerate(ng_order):
                if valid[i]:
                    n = int(side_sub["n"].iat[i])
                    ax.text(xs[i], ys[i] + 0.03,
                            f"n={n:,}", ha="center", fontsize=7, color=color)

        ax.axhline(baseline, color=PALETTE["dark"], linestyle="--", linewidth=1, alpha=0.5,
                   label=f"baseline {baseline:.2f}")
        ax.set_xticks(range(len(ng_order)))
        ax.set_xticklabels([f"NG={n}" for n in ng_order], fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("TP rate", fontsize=10)
        ax.set_title(sample_label(sample, "to_pileup"),
                     fontsize=12, fontweight="bold", color=PALETTE["dark"])
        ax.grid(True, alpha=0.25)
        ax.legend(loc="lower right", fontsize=8)
        annotate_power_tier(ax, sample, "to_pileup",
                            baseline_tp=baseline,
                            n_total=int(sub["n"].sum()))

    # Takeaway: Inner vs Outer gap direction across samples (NG=3 as reference)
    deltas = []
    for sample in samples:
        sub = cube[(cube["sample"] == sample) & (cube["AF_class"] == "Extreme") &
                   (cube["HPFineNGroups"] == 3)]
        inn = sub[sub["loh_side"] == "Inner"]["tp_rate"]
        out_ = sub[sub["loh_side"] == "Outer"]["tp_rate"]
        if len(inn) and len(out_) and pd.notna(inn.iloc[0]) and pd.notna(out_.iloc[0]):
            deltas.append((sample, inn.iloc[0] - out_.iloc[0]))

    if deltas:
        inner_higher = sum(1 for _, d in deltas if d > 0.02)
        outer_higher = sum(1 for _, d in deltas if d < -0.02)
        avg_gap = np.mean([d for _, d in deltas])
        samples_with_data = len(deltas)
        takeaway = (f"Inner–Outer gap（NG=3, Extreme AF）：Inner 高 {inner_higher}/{samples_with_data}，"
                    f"Outer 高 {outer_higher}/{samples_with_data}；"
                    f"平均 gap = {avg_gap:+.3f}。")
        if inner_higher >= samples_with_data * 0.6:
            takeaway += " ✓ Inner 方向優勢跨樣本一致。"
            colour = "green"
        elif outer_higher >= samples_with_data * 0.6:
            takeaway += " ↓ Outer 方向優勢跨樣本一致。"
            colour = "blue"
        else:
            takeaway += " ⚠ 方向不一致，需樣本別分析。"
            colour = "red"
    else:
        takeaway = "NG=3 × Extreme × Inner/Outer 資料不足以下結論"
        colour = "grey"

    takeaway_caption(fig, takeaway, color=colour, y=0.015, fontsize=11)

    fig.tight_layout(rect=[0, 0.05, 1, 0.93])
    out = OUT_DIR / "obs17_TO_inner_vs_outer_ng_gradient.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_anomalies(cube: pd.DataFrame) -> Path:
    """Identify anomaly cells: where a cell's TP rate deviates significantly from the
    cross-sample median for that (AF, loh_side, NG) combination."""
    samples = [s for s, _ in SAMPLE_MODE_TO]

    # Filter cells with n >= 50 (robust enough for comparison)
    robust = cube[cube["n"] >= 50].copy()
    # Compute cross-sample median TP rate per (AF, loh_side, NG)
    medians = (robust.groupby(["AF_class", "loh_side", "HPFineNGroups"])["tp_rate"]
                     .median().rename("median_tp").reset_index())
    robust = robust.merge(medians, on=["AF_class", "loh_side", "HPFineNGroups"])
    robust["deviation"] = robust["tp_rate"] - robust["median_tp"]
    robust["abs_deviation"] = robust["deviation"].abs()

    # Per-sample anomaly score = median of |deviation| across its cells
    per_sample = robust.groupby("sample").agg(
        n_cells=("tp_rate", "size"),
        median_abs_dev=("abs_deviation", "median"),
        max_abs_dev=("abs_deviation", "max"),
        mean_deviation=("deviation", "mean"),
    ).reindex(samples)

    # Top-15 most deviant cells
    top_anom = robust.reindex(robust["abs_deviation"].sort_values(ascending=False).index).head(15)

    fig, (ax_sum, ax_bar) = plt.subplots(1, 2, figsize=(20, 8),
                                        gridspec_kw={"width_ratios": [2, 3]})
    research_suptitle(
        fig,
        "TO mode 4D cells 的跨樣本異常偵測",
        context="左圖：各樣本在共享 (AF, loh_side, NG) 格上的 TP rate 偏離 cross-sample median 的幅度；"
                "右圖：Top-15 最偏離單元（|deviation| 最大）",
        y=0.96, fontsize=17,
    )

    # Left: bar per sample
    y_pos = np.arange(len(samples))
    mean_dev = per_sample["mean_deviation"].values
    colors = [PALETTE["green"] if v > 0.02 else (PALETTE["red"] if v < -0.02 else PALETTE["blue"])
              for v in mean_dev]
    ax_sum.barh(y_pos, mean_dev, color=colors, edgecolor=PALETTE["dark"])
    ax_sum.set_yticks(y_pos)
    ax_sum.set_yticklabels([f"{s} (n_cells={int(n)})"
                            for s, n in zip(samples, per_sample["n_cells"].values)])
    ax_sum.set_xlabel("Mean deviation from cross-sample median TP rate", fontsize=10)
    ax_sum.axvline(0, color=PALETTE["dark"], linewidth=1)
    ax_sum.set_title("各樣本偏離度總覽", fontsize=12, fontweight="bold", color=PALETTE["dark"])
    ax_sum.grid(True, axis="x", alpha=0.25)
    for i, v in enumerate(mean_dev):
        ax_sum.text(v + (0.003 if v >= 0 else -0.003), i,
                    f"{v:+.3f}", va="center",
                    ha="left" if v >= 0 else "right", fontsize=9)

    # Right: top anomalies
    ys = np.arange(len(top_anom))
    ax_bar.barh(ys, top_anom["deviation"].values,
                color=[PALETTE["green"] if v > 0 else PALETTE["red"]
                       for v in top_anom["deviation"].values],
                edgecolor=PALETTE["dark"])
    ax_bar.set_yticks(ys)
    labels = [f"{r['sample'][:10]:<10} | AF={r['AF_class'][:4]} | {r['loh_side']} | NG={r['HPFineNGroups']}"
              + f" | n={int(r['n']):,} | TP={r['tp_rate']:.2f} vs median {r['median_tp']:.2f}"
              for _, r in top_anom.iterrows()]
    ax_bar.set_yticklabels(labels, fontsize=8, family="monospace")
    ax_bar.set_xlabel("TP rate deviation from cross-sample median", fontsize=10)
    ax_bar.axvline(0, color=PALETTE["dark"], linewidth=1)
    ax_bar.invert_yaxis()
    ax_bar.grid(True, axis="x", alpha=0.25)
    ax_bar.set_title("Top-15 異常 cells（|dev| 最大）", fontsize=12, fontweight="bold",
                     color=PALETTE["dark"])

    # Takeaway: which sample is most anomalous?
    most_anomalous = per_sample["median_abs_dev"].idxmax()
    most_anomalous_val = per_sample["median_abs_dev"].max()
    least_anomalous = per_sample["median_abs_dev"].idxmin()
    takeaway = (f"最偏離樣本：{most_anomalous}（median |dev| = {most_anomalous_val:.3f}）；"
                f"最一致樣本：{least_anomalous}。Top-1 異常 cell："
                f"{top_anom.iloc[0]['sample']} {top_anom.iloc[0]['AF_class']}/"
                f"{top_anom.iloc[0]['loh_side']}/NG={top_anom.iloc[0]['HPFineNGroups']} "
                f"(deviation={top_anom.iloc[0]['deviation']:+.2f})。")
    colour = "red" if most_anomalous_val > 0.15 else ("blue" if most_anomalous_val > 0.08 else "green")
    takeaway_caption(fig, takeaway, color=colour, y=0.02, fontsize=11)

    fig.tight_layout(rect=[0, 0.06, 1, 0.93])
    out = OUT_DIR / "obs17_TO_anomalies.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out, per_sample, top_anom


def build_direction_summary(cube: pd.DataFrame) -> pd.DataFrame:
    """Per-sample direction vote: for each (AF, loh_side), what's the NG gradient direction?"""
    rows = []
    samples = [s for s, _ in SAMPLE_MODE_TO]
    ng_order = DIM_ORDERS["HPFineNGroups"]
    for sample in samples:
        for af in DIM_ORDERS["AF_class"]:
            for side in ["Inner", "Outer"]:
                sub = cube[(cube["sample"] == sample) & (cube["AF_class"] == af) &
                           (cube["loh_side"] == side)]
                rates = sub.set_index("HPFineNGroups")["tp_rate"].reindex(ng_order).dropna()
                if len(rates) < 2:
                    direction, gap = "insufficient", np.nan
                else:
                    gap = rates.iloc[-1] - rates.iloc[0]
                    if gap > 0.05:
                        direction = "up"
                    elif gap < -0.05:
                        direction = "down"
                    else:
                        direction = "flat"
                rows.append({
                    "sample": sample, "AF_class": af, "loh_side": side,
                    "NG_first": rates.iloc[0] if len(rates) else np.nan,
                    "NG_last": rates.iloc[-1] if len(rates) else np.nan,
                    "gap": gap, "direction": direction,
                    "n_total": int(sub["n"].sum()),
                })
    return pd.DataFrame(rows)


def main() -> None:
    apply_style()
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    df = load_master()
    to = prepare_to_df(df)
    print(f"[obs17] TO rows after prepare: {len(to)}")

    cube = build_4d_cube(to)
    tsv_cube = DATA_DIR / "obs17_TO_4d_cube.tsv"
    cube.to_csv(tsv_cube, sep="\t", index=False)
    print(f"[obs17] wrote {tsv_cube}  ({len(cube)} cells)")

    direction = build_direction_summary(cube)
    tsv_dir = DATA_DIR / "obs17_TO_direction_summary.tsv"
    direction.to_csv(tsv_dir, sep="\t", index=False)
    print(f"[obs17] wrote {tsv_dir}  ({len(direction)} rows)")
    print("\nDirection summary:")
    print(direction.to_string(index=False))

    out1 = plot_afxng_heatmap_by_lohside(cube)
    print(f"[obs17] wrote {out1}")
    out2 = plot_inner_vs_outer_ng_gradient(cube)
    print(f"[obs17] wrote {out2}")
    out3, per_sample, top_anom = plot_anomalies(cube)
    print(f"[obs17] wrote {out3}")
    print("\nPer-sample anomaly:")
    print(per_sample)
    print("\nTop-10 anomaly cells:")
    print(top_anom.head(10)[["sample","AF_class","loh_side","HPFineNGroups","n","tp_rate","median_tp","deviation"]].to_string(index=False))


if __name__ == "__main__":
    main()
