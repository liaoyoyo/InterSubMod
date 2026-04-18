#!/usr/bin/env python3
"""
LOH Read Threshold Analysis — 視覺化論證 ISM LOH 判定缺少最低 read 門檻的影響。

產出圖表：
  Fig01: LOH rate vs effective_hp (TO / Paired / Binomial null) — 三線比較
  Fig02: HP_Ratio 散佈圖 by ehp tier — 展示低 read 的 LOH 是隨機噪音
  Fig03: LOH excess heatmap (sample × ehp tier) — 跨樣本一致性
  Fig04: LOH 累積貢獻 — 低 read 的 LOH 佔整體多少
  Fig05: 高 read 區 TO vs Paired LOH — 論證 scaffold 品質才是主因
  Fig06: Binomial null vs observed violin — per-ehp 分佈比較
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import binom
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── Config ──
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260402_loh_read_threshold_analysis")
FIG_DIR = OUT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_SHORT = {"HCC1395": "HCC1395", "HCC1395_DORADO": "HCC1395_D", "COLO829": "COLO829",
                "H1437": "H1437", "H2009": "H2009", "HCC1937": "HCC1937", "HCC1954": "HCC1954"}

# ── Load ──
print("Loading master dataset...")
cols = ["mode", "truth_label", "HP1FamilyN", "HP2FamilyN", "HP_Ratio",
        "Potential_LOH", "caller_filter", "sample"]
df = pd.read_csv(MASTER, sep="\t", usecols=cols)
df["effective_hp"] = df["HP1FamilyN"] + df["HP2FamilyN"]

to = df[(df["mode"] == "to") & (df["caller_filter"] == "PASS")].copy()
paired = df[(df["mode"] == "paired") & (df["caller_filter"] == "PASS")].copy()
print(f"TO PASS: {len(to):,}, Paired PASS: {len(paired):,}")


def save_fig(fig, name):
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path.name}")


# ── Binomial null model ──
def binomial_null_loh_rate(n, p=0.5, threshold=0.1):
    """P(LOH) under binomial null with epsilon smoothing."""
    if n == 0:
        return 0.0
    p_loh = 0.0
    for x in range(n + 1):
        ratio = (x + 0.001) / (n + 0.002)
        if ratio < threshold or ratio > (1 - threshold):
            p_loh += binom.pmf(x, n, p)
    return p_loh


# ── Fig01: LOH rate vs effective_hp — TO / Paired / Binomial null ──
print("Fig01: LOH rate vs effective_hp...")
ehp_points = [1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300]

# Compute observed rates in rolling windows
def loh_rate_by_ehp_bins(data, centers):
    rates = []
    ns = []
    for c in centers:
        lo = max(0, int(c * 0.7))
        hi = int(c * 1.3) + 1
        if c <= 5:
            lo, hi = c, c + 1
        sub = data[(data["effective_hp"] >= lo) & (data["effective_hp"] < hi)]
        if len(sub) > 0:
            rates.append(sub["Potential_LOH"].mean() * 100)
            ns.append(len(sub))
        else:
            rates.append(np.nan)
            ns.append(0)
    return rates, ns

to_rates, to_ns = loh_rate_by_ehp_bins(to, ehp_points)
paired_rates, paired_ns = loh_rate_by_ehp_bins(paired, ehp_points)
null_rates = [binomial_null_loh_rate(n) * 100 for n in ehp_points]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={"height_ratios": [3, 1]})

ax1.plot(ehp_points, to_rates, "o-", color="#d62728", linewidth=2.5, markersize=7, label="TO observed", zorder=3)
ax1.plot(ehp_points, paired_rates, "s-", color="#1f77b4", linewidth=2.5, markersize=7, label="Paired observed", zorder=3)
ax1.plot(ehp_points, null_rates, "^--", color="#2ca02c", linewidth=2, markersize=6, label="Binomial null (p=0.5)", zorder=2)

# Shade the "no-threshold" danger zone
ax1.axvspan(0, 10, alpha=0.08, color="red", label="No min-read protection zone")
ax1.axhline(y=10, color="gray", linestyle=":", alpha=0.5)
ax1.annotate("LOH threshold = 10%", xy=(300, 11), fontsize=9, color="gray")

ax1.set_xscale("log")
ax1.set_xlim(0.8, 400)
ax1.set_ylim(-2, 105)
ax1.set_ylabel("LOH-like rate (%)", fontsize=13)
ax1.set_title("Fig01: LOH Rate vs Effective HP Reads — No Minimum Threshold Protection\n"
              "ISM is_potential_loh() has NO read count check → low-read LOH is noise",
              fontsize=13, fontweight="bold")
ax1.legend(fontsize=11, loc="upper right")
ax1.grid(True, alpha=0.3)

# Excess panel
excess = [t - n if not (np.isnan(t) or np.isnan(n)) else np.nan for t, n in zip(to_rates, null_rates)]
ax2.bar(range(len(ehp_points)), excess, color="#d62728", alpha=0.7, edgecolor="white")
ax2.set_xticks(range(len(ehp_points)))
ax2.set_xticklabels([str(x) for x in ehp_points], fontsize=9)
ax2.set_xlabel("effective_hp_reads", fontsize=13)
ax2.set_ylabel("TO excess over null (%)", fontsize=12)
ax2.axhline(0, color="black", linewidth=0.5)
ax2.set_title("Self-phasing excess: stable +36-64% across ALL read counts", fontsize=11)
ax2.grid(True, alpha=0.3, axis="y")

fig.tight_layout()
save_fig(fig, "fig01_loh_rate_vs_effective_hp")


# ── Fig02: HP_Ratio scatter by ehp tier — 4-panel ──
print("Fig02: HP_Ratio scatter by ehp tier...")
fig, axes = plt.subplots(2, 4, figsize=(20, 10))

tiers = [("ehp=1", 1, 2), ("ehp=2-4", 2, 5), ("ehp=5-9", 5, 10), ("ehp=10-19", 10, 20),
         ("ehp=20-29", 20, 30), ("ehp=30-49", 30, 50), ("ehp=50-99", 50, 100), ("ehp≥100", 100, 999999)]

for idx, (label, lo, hi) in enumerate(tiers):
    ax = axes[idx // 4, idx % 4]
    sub = to[(to["effective_hp"] >= lo) & (to["effective_hp"] < hi)].copy()
    if len(sub) > 2000:
        sub = sub.sample(2000, random_state=42)

    tp = sub[sub["truth_label"] == "TP"]
    fp = sub[sub["truth_label"] == "FP"]

    ax.scatter(tp["effective_hp"], tp["HP_Ratio"], s=8, alpha=0.4, c="#1f77b4", label=f"TP (n={len(tp)})")
    ax.scatter(fp["effective_hp"], fp["HP_Ratio"], s=8, alpha=0.4, c="#d62728", label=f"FP (n={len(fp)})")

    # LOH zones
    ax.axhline(0.1, color="orange", linestyle="--", alpha=0.6, linewidth=1)
    ax.axhline(0.9, color="orange", linestyle="--", alpha=0.6, linewidth=1)
    ax.axhspan(0, 0.1, alpha=0.06, color="red")
    ax.axhspan(0.9, 1.0, alpha=0.06, color="red")

    loh_pct = sub["Potential_LOH"].mean() * 100
    ax.set_title(f"{label}\nLOH={loh_pct:.0f}% (n={len(sub)})", fontsize=10, fontweight="bold")
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("HP_Ratio" if idx % 4 == 0 else "")
    ax.legend(fontsize=7, loc="center right")
    ax.grid(True, alpha=0.2)

fig.suptitle("Fig02: HP_Ratio Distribution by Effective HP Reads (TO mode)\n"
             "Low reads → HP_Ratio locked at 0 or 1 → trivial LOH; High reads → self-phasing LOH",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
save_fig(fig, "fig02_hp_ratio_scatter_by_ehp_tier")


# ── Fig03: LOH excess heatmap (sample × ehp tier) ──
print("Fig03: LOH excess heatmap...")
tier_bins = [(0, 5, "0-4"), (5, 10, "5-9"), (10, 20, "10-19"), (20, 30, "20-29"),
             (30, 50, "30-49"), (50, 100, "50-99"), (100, 200, "100-199"), (200, 999999, "200+")]

# TO LOH rate per sample per tier
heatmap_data = []
for sample in SAMPLE_ORDER:
    row = []
    for lo, hi, label in tier_bins:
        sub = to[(to["sample"] == sample) & (to["effective_hp"] >= lo) & (to["effective_hp"] < hi)]
        if len(sub) >= 10:
            row.append(sub["Potential_LOH"].mean() * 100)
        else:
            row.append(np.nan)
    heatmap_data.append(row)

hm_df = pd.DataFrame(heatmap_data, index=[SAMPLE_SHORT[s] for s in SAMPLE_ORDER],
                      columns=[t[2] for t in tier_bins])

# Paired average for reference
paired_ref = []
for lo, hi, label in tier_bins:
    sub = paired[(paired["effective_hp"] >= lo) & (paired["effective_hp"] < hi)]
    paired_ref.append(sub["Potential_LOH"].mean() * 100 if len(sub) >= 10 else np.nan)

fig, ax = plt.subplots(figsize=(14, 6))
im = ax.imshow(hm_df.values, cmap="YlOrRd", aspect="auto", vmin=0, vmax=80)

for i in range(len(SAMPLE_ORDER)):
    for j in range(len(tier_bins)):
        val = hm_df.values[i, j]
        if not np.isnan(val):
            color = "white" if val > 50 else "black"
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=9,
                    fontweight="bold", color=color)

ax.set_xticks(range(len(tier_bins)))
ax.set_xticklabels([t[2] for t in tier_bins], fontsize=11)
ax.set_yticks(range(len(SAMPLE_ORDER)))
ax.set_yticklabels([SAMPLE_SHORT[s] for s in SAMPLE_ORDER], fontsize=11)
ax.set_xlabel("effective_hp_reads tier", fontsize=13)
ax.set_ylabel("Sample", fontsize=13)

# Add paired reference row below
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(range(len(tier_bins)))
paired_labels = [f"P:{v:.0f}%" if not np.isnan(v) else "—" for v in paired_ref]
ax2.set_xticklabels(paired_labels, fontsize=9, color="#1f77b4")
ax2.tick_params(top=True, bottom=False)

cbar = fig.colorbar(im, ax=ax, shrink=0.8, label="LOH-like rate (%)")

ax.set_title("Fig03: TO LOH Rate by Sample × Effective HP Tier\n"
             "Top axis (blue) = Paired reference. All tiers show self-phasing excess over Paired.",
             fontsize=13, fontweight="bold", pad=30)
fig.tight_layout()
save_fig(fig, "fig03_loh_excess_heatmap_sample_ehp")


# ── Fig04: Cumulative LOH contribution ──
print("Fig04: Cumulative LOH contribution...")
to_sorted = to.sort_values("effective_hp").copy()
to_sorted["loh_int"] = to_sorted["Potential_LOH"].astype(int)

total_loh = to_sorted["loh_int"].sum()
total_regions = len(to_sorted)
to_sorted["cum_loh"] = to_sorted["loh_int"].cumsum() / total_loh * 100
to_sorted["cum_regions"] = np.arange(1, len(to_sorted) + 1) / total_regions * 100

fig, ax = plt.subplots(figsize=(12, 6))

# Mark key thresholds
for thresh, color, label in [(10, "red", "ehp<10"), (30, "orange", "ehp<30"), (50, "green", "ehp<50")]:
    mask = to_sorted["effective_hp"] < thresh
    n_below = mask.sum()
    loh_below = to_sorted.loc[mask, "loh_int"].sum()
    pct_regions = n_below / total_regions * 100
    pct_loh = loh_below / total_loh * 100
    ax.axvline(pct_regions, color=color, linestyle="--", alpha=0.6, linewidth=1.5)
    ax.annotate(f"{label}\n{pct_regions:.1f}% regions\n{pct_loh:.1f}% LOH",
                xy=(pct_regions, pct_loh), xytext=(pct_regions + 3, pct_loh + 10),
                fontsize=9, color=color, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=color, lw=1.5))

ax.plot(to_sorted["cum_regions"].values, to_sorted["cum_loh"].values,
        color="#d62728", linewidth=2)
ax.plot([0, 100], [0, 100], "k--", alpha=0.3, label="Uniform (no bias)")
ax.set_xlabel("Cumulative % of TO PASS regions (sorted by effective_hp)", fontsize=12)
ax.set_ylabel("Cumulative % of LOH-like regions", fontsize=12)
ax.set_title("Fig04: Cumulative LOH Contribution by Read Count\n"
             "Low-read regions contribute disproportionate LOH, but bulk comes from high-read self-phasing",
             fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.grid(True, alpha=0.3)
fig.tight_layout()
save_fig(fig, "fig04_cumulative_loh_contribution")


# ── Fig05: High-read LOH — TO vs Paired with truth label ──
print("Fig05: High-read LOH TO vs Paired...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for idx, (title, lo, hi) in enumerate([("ehp 50-99", 50, 100), ("ehp 100-199", 100, 200), ("ehp 200+", 200, 999999)]):
    ax = axes[idx]
    t_sub = to[(to["effective_hp"] >= lo) & (to["effective_hp"] < hi)]
    p_sub = paired[(paired["effective_hp"] >= lo) & (paired["effective_hp"] < hi)]

    categories = ["TO TP", "TO FP", "Paired TP", "Paired FP"]
    loh_rates = []
    ns = []
    for label_name, data, truth in [("TO TP", t_sub, "TP"), ("TO FP", t_sub, "FP"),
                                     ("Paired TP", p_sub, "TP"), ("Paired FP", p_sub, "FP")]:
        sub = data[data["truth_label"] == truth]
        rate = sub["Potential_LOH"].mean() * 100 if len(sub) > 0 else 0
        loh_rates.append(rate)
        ns.append(len(sub))

    colors = ["#d62728", "#ff7f0e", "#1f77b4", "#17becf"]
    bars = ax.bar(categories, loh_rates, color=colors, edgecolor="white", linewidth=1.5)

    for bar, rate, n in zip(bars, loh_rates, ns):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{rate:.1f}%\n(n={n:,})", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_ylim(0, max(loh_rates) * 1.4)
    ax.set_ylabel("LOH-like rate (%)" if idx == 0 else "")
    ax.set_title(f"{title}", fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.2, axis="y")

fig.suptitle("Fig05: LOH Rate at HIGH Read Counts — Self-Phasing Not Read-Count Noise\n"
             "Even at 200+ reads, TO LOH is 14× Paired → scaffold quality is the real problem",
             fontsize=14, fontweight="bold")
fig.tight_layout()
save_fig(fig, "fig05_high_read_loh_to_vs_paired")


# ── Fig06: Two-problem decomposition ──
print("Fig06: Two-problem decomposition...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left: Problem A — Low read noise
tier_labels = ["1", "2-4", "5-9", "10-19", "20-29", "30-49", "50-99", "100+"]
tier_ranges = [(1, 2), (2, 5), (5, 10), (10, 20), (20, 30), (30, 50), (50, 100), (100, 999999)]

to_loh_rates = []
null_loh_rates = []
for lo, hi in tier_ranges:
    sub = to[(to["effective_hp"] >= lo) & (to["effective_hp"] < hi)]
    to_loh_rates.append(sub["Potential_LOH"].mean() * 100 if len(sub) > 0 else 0)
    # Compute null rate for midpoint
    mid = (lo + hi) // 2 if hi < 999999 else 150
    null_loh_rates.append(binomial_null_loh_rate(mid) * 100)

x = np.arange(len(tier_labels))
width = 0.35
bars1 = ax1.bar(x - width / 2, null_loh_rates, width, label="Binomial null (random HP)", color="#2ca02c", alpha=0.7)
bars2 = ax1.bar(x + width / 2, to_loh_rates, width, label="TO observed", color="#d62728", alpha=0.7)

# Annotate the excess
for i in range(len(tier_labels)):
    excess = to_loh_rates[i] - null_loh_rates[i]
    if excess > 5:
        ax1.annotate(f"+{excess:.0f}%", xy=(x[i] + width / 2, to_loh_rates[i]),
                     xytext=(0, 5), textcoords="offset points", ha="center", fontsize=8,
                     color="#d62728", fontweight="bold")

ax1.set_xticks(x)
ax1.set_xticklabels(tier_labels, fontsize=10)
ax1.set_xlabel("effective_hp_reads", fontsize=12)
ax1.set_ylabel("LOH-like rate (%)", fontsize=12)
ax1.set_title("Problem A: Low-Read Noise\n"
              "ehp=1 is 100% LOH by math;\n"
              "BUT excess is stable +36-64% at ALL tiers → not just noise",
              fontsize=11, fontweight="bold")
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3, axis="y")

# Right: Problem B — Self-phasing scaffold
to_rates_b = []
paired_rates_b = []
for lo, hi in tier_ranges:
    t_sub = to[(to["effective_hp"] >= lo) & (to["effective_hp"] < hi)]
    p_sub = paired[(paired["effective_hp"] >= lo) & (paired["effective_hp"] < hi)]
    to_rates_b.append(t_sub["Potential_LOH"].mean() * 100 if len(t_sub) > 0 else 0)
    paired_rates_b.append(p_sub["Potential_LOH"].mean() * 100 if len(p_sub) > 0 else 0)

bars3 = ax2.bar(x - width / 2, paired_rates_b, width, label="Paired", color="#1f77b4", alpha=0.7)
bars4 = ax2.bar(x + width / 2, to_rates_b, width, label="TO", color="#d62728", alpha=0.7)

# Ratio annotations
for i in range(len(tier_labels)):
    if paired_rates_b[i] > 0:
        ratio = to_rates_b[i] / paired_rates_b[i]
        if ratio > 1.3:
            ax2.annotate(f"{ratio:.1f}×", xy=(x[i] + width / 2, to_rates_b[i]),
                         xytext=(0, 5), textcoords="offset points", ha="center", fontsize=8,
                         color="#d62728", fontweight="bold")

ax2.set_xticks(x)
ax2.set_xticklabels(tier_labels, fontsize=10)
ax2.set_xlabel("effective_hp_reads", fontsize=12)
ax2.set_ylabel("LOH-like rate (%)", fontsize=12)
ax2.set_title("Problem B: Self-Phasing Scaffold\n"
              "TO/Paired ratio grows at high reads (14× at 100+)\n"
              "→ Scaffold quality, not read count, is the root cause",
              fontsize=11, fontweight="bold")
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3, axis="y")

fig.suptitle("Fig06: Two Distinct Problems in TO LOH Determination",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
save_fig(fig, "fig06_two_problem_decomposition")


# ── Summary TSV ──
print("\nSaving summary TSV...")
summary_rows = []
for lo, hi, label in [(0, 1, "0"), (1, 2, "1"), (2, 5, "2-4"), (5, 10, "5-9"),
                       (10, 20, "10-19"), (20, 30, "20-29"), (30, 50, "30-49"),
                       (50, 100, "50-99"), (100, 200, "100-199"), (200, 999999, "200+")]:
    t_sub = to[(to["effective_hp"] >= lo) & (to["effective_hp"] < hi)]
    p_sub = paired[(paired["effective_hp"] >= lo) & (paired["effective_hp"] < hi)]
    mid = (lo + min(hi, 300)) // 2
    null_rate = binomial_null_loh_rate(mid) * 100

    t_tp = t_sub[t_sub["truth_label"] == "TP"]
    t_fp = t_sub[t_sub["truth_label"] == "FP"]
    p_tp = p_sub[p_sub["truth_label"] == "TP"]
    p_fp = p_sub[p_sub["truth_label"] == "FP"]

    summary_rows.append({
        "ehp_tier": label,
        "to_n": len(t_sub),
        "to_loh_rate": t_sub["Potential_LOH"].mean() if len(t_sub) > 0 else np.nan,
        "to_tp_loh_rate": t_tp["Potential_LOH"].mean() if len(t_tp) > 0 else np.nan,
        "to_fp_loh_rate": t_fp["Potential_LOH"].mean() if len(t_fp) > 0 else np.nan,
        "paired_n": len(p_sub),
        "paired_loh_rate": p_sub["Potential_LOH"].mean() if len(p_sub) > 0 else np.nan,
        "paired_tp_loh_rate": p_tp["Potential_LOH"].mean() if len(p_tp) > 0 else np.nan,
        "paired_fp_loh_rate": p_fp["Potential_LOH"].mean() if len(p_fp) > 0 else np.nan,
        "binomial_null_rate": null_rate / 100,
        "to_excess_over_null": (t_sub["Potential_LOH"].mean() - null_rate / 100) if len(t_sub) > 0 else np.nan,
    })

summary_df = pd.DataFrame(summary_rows)
summary_path = OUT / "loh_read_threshold_summary.tsv"
summary_df.to_csv(summary_path, sep="\t", index=False, float_format="%.4f")
print(f"  Saved: {summary_path.name}")

print("\n=== Done. 6 figures + 1 TSV ===")
