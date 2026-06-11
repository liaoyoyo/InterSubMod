#!/usr/bin/env python3
"""
41 — A6 Fig2: regime-A credible-subset gate evidence (the headline G1 figure).
Panel A: blind-ARI distribution regime-A (somatic HP-axis) vs germline-het null (baseline-allelic).
Panel B: 4-gate funnel (75 regime-A -> 62 evaluated -> 23 Tier-A pass) + gate-pass scorecard.
Panel C: LOH 3-class diagnostic stacked bar (self-phasing / candidate-subclone / CN-regression).

Reads regimeA_credible_probe.json + regimeA_hardening.json + loh_diagnostic_classifier.json.
No BAM. Single-sample HCC1395 Tier-A. Output: figures/fig2_gate_evidence.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts")
from lib.plot_setup import setup_plot_style
setup_plot_style()

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
PROBE = WS / "genome_survey_v2/regimeA_credible_probe.json"
HARD = WS / "genome_survey_v2/regimeA_hardening.json"
LOHC = WS / "genome_survey_v2/loh_diagnostic_classifier.json"
RESID = WS / "genome_survey_v2/regimeA_residual_controls.json"
FIG2 = WS / "figures/fig2_gate_evidence.png"
BLUE, RED, GREEN, ORANGE, GREY, PURPLE = "#1E3A8A", "#C2410C", "#15803D", "#D97757", "#9CA3AF", "#7C3AED"


def main():
    probe = json.load(open(PROBE)); hard = json.load(open(HARD)); lohc = json.load(open(LOHC))
    resid = json.load(open(RESID)); rsum = resid['summary']
    n_survive_both = int(rsum['n_survive_BOTH'].split("/")[0])
    ps = probe['summary']; hs = hard['summary']
    regimeA_ari = [r['ari'] for r in probe['loci']]
    het_ari = [r['ari'] for r in hard['het_null']]
    m2 = hs['M2_regimeA_vs_hetnull']

    fig, ax = plt.subplots(1, 3, figsize=(18, 6), dpi=140)

    # ---- Panel A: ARI distribution regime-A vs het-null ----
    a = ax[0]
    parts = a.violinplot([regimeA_ari, het_ari], positions=[1, 2], showmedians=True, widths=0.7)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(GREEN if i == 0 else GREY); pc.set_alpha(0.55)
    for key in ('cmedians', 'cbars', 'cmins', 'cmaxes'):
        if key in parts:
            parts[key].set_color("#333")
    # jitter points
    for i, (vals, col) in enumerate([(regimeA_ari, GREEN), (het_ari, GREY)]):
        x = np.random.normal(i + 1, 0.05, len(vals))
        a.scatter(x, vals, s=12, c=col, alpha=0.5, edgecolors="none", zorder=3)
    a.axhline(0.30, color=RED, ls="--", lw=1, alpha=0.7)
    a.text(2.45, 0.31, "ARI=0.30\ncluster 門檻", fontsize=9, color=RED, ha="right")
    a.set_xticks([1, 2])
    a.set_xticklabels([f"regime-A\n(somatic HP-axis)\nmed={m2['regimeA_median_ari']:.3f}",
                       f"germline-het null\n(baseline-allelic)\nmed={m2['hetnull_median_ari']:.3f}"], fontsize=10)
    a.set_ylabel("blind-ARI (M1, 盲分群恢復 somatic/germline split)", fontsize=11)
    a.set_title(f"A · M2 — regime-A 顯著高於 baseline\n"
                f"MW p={m2['mw_p_greater']:.1e} · Cliff δ={m2['cliffs_delta']:.2f} · "
                f"med diff CI [{m2['median_diff_ci95'][0]:.2f},{m2['median_diff_ci95'][1]:.2f}]",
                fontsize=12, fontweight="bold")

    # ---- Panel B: gate funnel + scorecard ----
    b = ax[1]
    stages = [f"regime-A\n(n={ps['regimeA_n_total']})", f"可評估\n(n={ps['n_evaluated']})",
              f"length-placebo\npass (n={ps['n_pass_tierA']})", f"全battery\n存活 (n={n_survive_both})"]
    vals = [ps['regimeA_n_total'], ps['n_evaluated'], ps['n_pass_tierA'], n_survive_both]
    cols = [GREY, BLUE, "#86C5A0", GREEN]
    b.bar(range(4), vals, color=cols, alpha=0.85, width=0.6)
    for i, v in enumerate(vals):
        b.text(i, v + 1, str(v), ha="center", fontsize=12, fontweight="bold")
    b.set_xticks(range(4)); b.set_xticklabels(stages, fontsize=9)
    b.set_ylabel("位點數", fontsize=11)
    gate_txt = (f"gate (regime-A credible 子集):\n"
                f"  ✓ M1 blind-ARI median {ps['median_ari']:.3f}\n"
                f"  ✓ M2 vs het-null: Cliff δ={m2['cliffs_delta']:.2f} (p={m2['mw_p_greater']:.0e})\n"
                f"  ✓ M3 rarefied sil {hs['M3_rarefied_silhouette']['regimeA_pass_median']:.2f} vs {hs['M3_rarefied_silhouette']['hetnull_median']:.2f}\n"
                f"  ✓ M5 coverage ρ={ps['M5_spearman_ari_vs_logcov']['rho']:.2f} (NS)\n"
                f"  ✓ M4c CpG-context: {rsum['M4c_CpGcontext_survive']}\n"
                f"  ⚠ M8-strong collider: {rsum['M8strong_random_anchor_pass']}\n"
                f"     (7 刷掉, 含 SOX2/HOTTIP/SDHAF1)")
    b.text(0.03, 0.97, gate_txt, transform=b.transAxes, fontsize=9, va="top",
           bbox=dict(boxstyle="round", facecolor="#F0FDF4", edgecolor=GREEN, alpha=0.9))
    b.set_title("B · 全 artifact gate funnel (單樣本 Tier-A; 強化collider刷掉好看基因)", fontsize=11.5, fontweight="bold")
    b.set_ylim(0, ps['regimeA_n_total'] * 1.15)

    # ---- Panel C: LOH 3-class diagnostic ----
    c = ax[2]
    ls = lohc['summary']['class_counts']
    order = ["self_phasing_artifact", "candidate_subclone", "CN_regression", "intermediate"]
    labels = ["self-phasing\nartifact", "candidate\nsubclone", "CN /\nregression", "inter"]
    cc = [GREY, GREEN, ORANGE, "#C9CDD4"]
    vals = [ls.get(k, 0) for k in order]
    tot = sum(vals)
    bars = c.bar(range(len(order)), vals, color=cc, alpha=0.85, width=0.65)
    for i, v in enumerate(vals):
        if v:
            c.text(i, v + 0.5, f"{v}\n({100*v/tot:.0f}%)", ha="center", fontsize=10, fontweight="bold")
    c.set_xticks(range(len(order))); c.set_xticklabels(labels, fontsize=10)
    c.set_ylabel("LOH-regime 位點數 (n=%d evaluated)" % tot, fontsize=11)
    c.set_title("C · G2 — LOH 表觀雙 haplotype 成因診斷\n"
                "72% = self-phasing artifact (非真 cluster)；18% candidate subclone",
                fontsize=12, fontweight="bold")
    c.set_ylim(0, max(vals) * 1.25)

    fig.suptitle("Fig2 · G1 credible-subset 有真 somatic 甲基 cluster (refine 5/31) + G2 LOH 成因診斷 · HCC1395 單樣本 Tier-A",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(FIG2, facecolor="white", bbox_inches="tight")
    print(f"Fig2 -> {FIG2} ({FIG2.stat().st_size//1024} KB)")
    print(f"  Panel A: regime-A med ARI {m2['regimeA_median_ari']:.3f} vs het-null {m2['hetnull_median_ari']:.3f}")
    print(f"  Panel B: funnel {ps['regimeA_n_total']}->{ps['n_evaluated']}->{ps['n_pass_tierA']}")
    print(f"  Panel C: LOH classes {ls}")


if __name__ == "__main__":
    main()
