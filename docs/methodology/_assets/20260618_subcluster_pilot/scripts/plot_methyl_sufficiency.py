"""[甲基輔助有效性審查圖] 3-panel: (a) 充分度漏斗 (b) power 審查 (c) by-shape rate。
讀 sm_methyl_sufficiency_audit.json (§13-A 由 JSON 注入)。English labels (避 CJK font 依賴)。
"""
import sys
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ACCENT = "#D97757"
GREY = "#9a9183"
DARK = "#141413"
RED = "#c0392b"
GREEN = "#4a7c59"


def main(audit_path, out_path):
    d = json.load(open(audit_path))
    F = d["funnel"]
    P = d["power_audit"]
    S = d["by_tree_shape"]

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))

    # (a) sufficiency funnel
    labels = ["target\n(struct,CN-clean)", "tested\n(>=2 geno-pop)", "powered\n(CpG>=10,popB>=5)", "corroborated\n(>=1 sig CpG)"]
    vals = [F["n_target_structured_CNclean"], F["n_tested(>=2 geno-pop)"],
            F["n_powered(cpg>=10 & popB>=5)"], F["n_corroborated(>=1 sig CpG)"]]
    colors = [GREY, GREY, ACCENT, GREEN]
    bars = ax[0].barh(range(len(vals))[::-1], vals, color=colors, edgecolor=DARK)
    ax[0].set_yticks(range(len(vals))[::-1])
    ax[0].set_yticklabels(labels, fontsize=9)
    for b, v in zip(bars, vals):
        ax[0].text(v + 8, b.get_y() + b.get_height() / 2, str(v), va="center", fontsize=10, fontweight="bold")
    ax[0].set_title("(a) Methylation data-sufficiency funnel", fontsize=11, fontweight="bold")
    ax[0].set_xlim(0, max(vals) * 1.18)
    ax[0].spines[["top", "right"]].set_visible(False)
    pct = 100 * F["n_corroborated(>=1 sig CpG)"] / F["n_tested(>=2 geno-pop)"]
    ax[0].text(0.98, 0.04, f"corroborated = {pct:.1f}% of tested\n(0 new partition; cis-control 0/740\nnon-evaluable => subclone-specific UNVERIFIED)",
               transform=ax[0].transAxes, ha="right", va="bottom", fontsize=8.5, color=RED,
               bbox=dict(boxstyle="round", fc="#fdf0ee", ec=RED))

    # (b) power dose-response — the key panel (corrected: signal exists, power-gated)
    D = P["dose_response"]
    bins = ["popB_n[5,8)", "popB_n[8,12)", "popB_n[12,20)", "popB_n>=20"]
    blab = ["5-7", "8-11", "12-19", ">=20"]
    rates = [D[b]["rate"] * 100 for b in bins]
    ns = [f'{D[b]["corrob"]}/{D[b]["n"]}' for b in bins]
    cmap = [GREY, GREY, ACCENT, GREEN]
    bars = ax[1].bar(range(4), rates, color=cmap, edgecolor=DARK)
    ax[1].set_xticks(range(4))
    ax[1].set_xticklabels(blab, fontsize=9)
    for b, r, n in zip(bars, rates, ns):
        ax[1].text(b.get_x() + b.get_width() / 2, r + 1, f"{r:.1f}%\n{n}", ha="center", fontsize=9, fontweight="bold")
    ax[1].set_title("(b) Power dose-response: corroboration rises steeply with\ncoverage => signal EXISTS, power-gated (not no-signal)",
                    fontsize=10.5, fontweight="bold")
    ax[1].set_xlabel("smaller-population read depth (popB_n)")
    ax[1].set_ylabel("% corroborated")
    ax[1].set_ylim(0, max(rates) * 1.25)
    ax[1].spines[["top", "right"]].set_visible(False)
    starv = D["power_starved[5,12)_share_of_powered"] * 100
    ax[1].text(0.02, 0.97, f"but only {D['popB_n>=20']['n']} regions reach popB_n>=20;\n{starv:.0f}% of 'powered' is power-starved [5,12)",
               transform=ax[1].transAxes, ha="left", va="top", fontsize=8.5, color=RED,
               bbox=dict(boxstyle="round", fc="#fdf0ee", ec=RED))

    # (c) by shape
    order = ["full_tree", "linear_nested", "sibling_only", "co_linked_lineage"]
    order = [o for o in order if o in S]
    rates = [S[o]["rate"] * 100 for o in order]
    ns = [f'{S[o]["corroborated"]}/{S[o]["tested"]}' for o in order]
    bars = ax[2].bar(range(len(order)), rates, color=ACCENT, edgecolor=DARK)
    ax[2].set_xticks(range(len(order)))
    ax[2].set_xticklabels([o.replace("_", "\n") for o in order], fontsize=8.5)
    for b, r, n in zip(bars, rates, ns):
        ax[2].text(b.get_x() + b.get_width() / 2, r + 0.2, f"{r:.1f}%\n{n}", ha="center", fontsize=8.5)
    ax[2].set_title("(c) Corroboration rate by tree shape", fontsize=11, fontweight="bold")
    ax[2].set_ylabel("% corroborated")
    ax[2].set_ylim(0, max(rates) * 1.35)
    ax[2].spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    print("saved", out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
