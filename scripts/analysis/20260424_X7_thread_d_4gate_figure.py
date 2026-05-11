#!/usr/bin/env python3
"""
X7 — Thread D 4-gate integrated figure.

Visualizes the full verification chain for "LOH-constrained phasing" mechanism:

Panel A (descriptive):  obs18 per-sample Inner same_HP1 vs Outer cross_het TP rate (6 TO samples)
Panel B (formal stats): Wilcoxon signed-rank + bootstrap 95% CI (B1 → replicated in X5)
Panel C (paired NC):    B3 paired-mode gap collapse (H-D3 confirmation)
Panel D (flag-on NC):   X3 bucket physical collapse (same_HP1 n=219→0)

Output: research/tpfp_loh_af_kde_discrimination/figures/new/x7/thread_d_4gate.png
"""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

DATA_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/figures/new/x7")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load():
    with open(DATA_DIR / "obs18_wilcoxon_B1.json") as f:
        b1 = json.load(f)
    with open(DATA_DIR / "B3_wilcoxon_gap_stats.json") as f:
        b3 = json.load(f)
    with open(DATA_DIR / "X3_flag_onoff_NC.json") as f:
        x3 = json.load(f)
    with open(DATA_DIR / "X5_crosssample_summary.json") as f:
        x5 = json.load(f)
    return b1, b3, x3, x5


def main():
    b1, b3, x3, x5 = load()

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        "Thread D — LOH-constrained Phasing: 4-gate Verification Chain\n"
        "obs18 descriptive → formal stats → 2× negative controls",
        fontsize=15, fontweight='bold', y=0.98,
    )

    # Panel A: obs18 Inner vs Outer (B1)
    ax = axes[0, 0]
    per = b1["per_sample"]
    samples = [p["sample"] for p in per]
    inner = [p["inner_same_HP1_tp_rate"] for p in per]
    outer = [p["outer_cross_het_tp_rate"] for p in per]
    x = np.arange(len(samples))
    w = 0.35
    ax.bar(x - w/2, inner, w, label="Inner same_HP1 (LOH + same-hap)",
           color="#2E86AB", edgecolor='black')
    ax.bar(x + w/2, outer, w, label="Outer cross_het (non-LOH + cross-hap)",
           color="#E63946", edgecolor='black')
    for i, (a, b) in enumerate(zip(inner, outer)):
        ax.annotate(f"Δ=+{a-b:.2f}", (i, max(a, b) + 0.03),
                    ha='center', fontsize=9, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=20, ha='right', fontsize=9)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("TP rate")
    ax.set_title("A · obs18 6/6 positive gap (descriptive)",
                 fontsize=12, fontweight='bold')
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    # Panel B: Wilcoxon + CI (B1 + X5 replication)
    ax = axes[0, 1]
    gaps = b1["gaps"]
    x5_pivot = x5.get("obs18_pivot", [])
    x5_gaps = [p["gap"] for p in x5_pivot if p["gap"] is not None]

    positions = [1, 2]
    bp = ax.boxplot([gaps, x5_gaps], positions=positions, widths=0.5,
                    patch_artist=True, showmeans=True,
                    boxprops=dict(facecolor='#8ECAE6', alpha=0.6),
                    medianprops=dict(color='#E63946', linewidth=2))
    for i, data in enumerate([gaps, x5_gaps]):
        ax.scatter([positions[i]] * len(data), data, color='black', s=40, zorder=3)
    wil = b1["wilcoxon"]
    ax.axhline(0, color='gray', linestyle='--', alpha=0.6)
    ax.set_xticks(positions)
    ax.set_xticklabels(["B1\n(pre-KDE master)\nn=6", "X5\n(post-KDE master)\nn=6"])
    ax.set_ylabel("TP gap (Inner same_HP1 − Outer cross_het)")
    ax.set_title(f"B · Wilcoxon: W={wil['statistic_W']:.0f}, p={wil['p_value']:.4f} (both n=6)",
                 fontsize=12, fontweight='bold')
    ci = b1["bootstrap_median_gap_95ci"]
    ax.axhspan(ci['ci_low'], ci['ci_high'], alpha=0.1, color='green',
               label=f"Bootstrap 95% CI [{ci['ci_low']:.2f}, {ci['ci_high']:.2f}]")
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    # Panel C: B3 paired-mode NC (gap collapse to 0)
    ax = axes[1, 0]
    b3_samples = list(b3["paired_gap_values"].keys())
    b3_paired = [b3["paired_gap_values"][s] for s in b3_samples]
    b3_to_dict = b3["to_gap_values"]
    b3_to = [b3_to_dict[s] if b3_to_dict[s] is not None else np.nan for s in b3_samples]

    x = np.arange(len(b3_samples))
    w = 0.35
    ax.bar(x - w/2, b3_to, w, label="TO mode (gap present)", color='#F77F00', edgecolor='black')
    ax.bar(x + w/2, b3_paired, w, label="Paired mode (gap ≈0)", color='#06A77D', edgecolor='black')
    ax.axhline(0, color='black', linewidth=0.8)
    ax.axhline(b3["to_gap_median"], color='#F77F00', linestyle='--', alpha=0.6,
               label=f"TO median {b3['to_gap_median']:.2f}")
    ax.axhline(b3["paired_gap_median"], color='#06A77D', linestyle='--', alpha=0.6,
               label=f"Paired median ≈0")
    ax.set_xticks(x)
    ax.set_xticklabels(b3_samples, rotation=20, ha='right', fontsize=9)
    ax.set_ylabel("TP gap")
    ax.set_title(f"C · B3 paired-mode NC: gap collapses\n(Paired Wilcoxon p={b3['paired_gap_wilcoxon_p']:.2f} n.s.)",
                 fontsize=12, fontweight='bold')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(axis='y', alpha=0.3)

    # Panel D: X3 flag-on NC bucket physical collapse
    ax = axes[1, 1]
    x3_off = x3["data"]["flag_off"]
    x3_on = x3["data"]["flag_on"]
    ng_off = x3_off["ng_distribution_top6"]
    ng_on = x3_on["ng_distribution_top6"]

    ng_keys = sorted(set(list(ng_off.keys()) + list(ng_on.keys())))
    ng_keys = [int(k) for k in ng_keys if int(k) <= 5]
    off_vals = [ng_off.get(str(k), ng_off.get(k, 0)) for k in ng_keys]
    on_vals = [ng_on.get(str(k), ng_on.get(k, 0)) for k in ng_keys]

    x = np.arange(len(ng_keys))
    w = 0.35
    b1_ = ax.bar(x - w/2, off_vals, w, label="flag=off (original)",
                 color='#4F86F7', edgecolor='black')
    b2_ = ax.bar(x + w/2, on_vals, w, label="flag=on (--germline-hp-only)",
                 color='#C94A4A', edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels([f"NG={k}" for k in ng_keys])
    ax.set_ylabel("Region count (HCC1395 TO)")
    ax.set_title(f"D · X3 flag=on NC: NG≥3 physical collapse\n(same_HP1 bucket: n=219 → n=0)",
                 fontsize=12, fontweight='bold')
    for rect, v in zip(b1_, off_vals):
        ax.annotate(f"{v:,}", (rect.get_x() + rect.get_width()/2, rect.get_height()),
                    ha='center', va='bottom', fontsize=8)
    for rect, v in zip(b2_, on_vals):
        ax.annotate(f"{v:,}", (rect.get_x() + rect.get_width()/2, rect.get_height()),
                    ha='center', va='bottom', fontsize=8)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = OUT_DIR / "thread_d_4gate.png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure → {out}")

    # Summary table
    summary = {
        "gate_A_descriptive": {"samples": 6, "all_positive": True, "median_gap": b1["descriptives"]["median_gap"]},
        "gate_B_formal_stats": {"W": wil["statistic_W"], "p": wil["p_value"],
                                "ci_95": [ci['ci_low'], ci['ci_high']], "replicated_in_x5": True},
        "gate_C_paired_NC": {"paired_median_gap": b3["paired_gap_median"],
                             "paired_wilcoxon_p": b3["paired_gap_wilcoxon_p"],
                             "interpretation": "gap collapses under matched-normal germline caller"},
        "gate_D_flag_on_NC": {"ng_ge_3_off": sum(v for k, v in ng_off.items() if int(k) >= 3),
                              "ng_ge_3_on": sum(v for k, v in ng_on.items() if int(k) >= 3),
                              "same_HP1_off_n": x3_off["full_grid"],
                              "interpretation": "same-hap bucket physically disappears"},
    }
    out_json = OUT_DIR / "thread_d_4gate_summary.json"
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"Summary → {out_json}")


if __name__ == "__main__":
    main()
