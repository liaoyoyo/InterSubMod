#!/usr/bin/env python3
"""
100 - Catalog observation distribution figures P1-P6 (T-OBS-1, schema section C).
"Forest view" (see-the-forest-and-the-trees): per-TAG counts + per-column distributions.

Reads landed catalog (script 98) + discovery_venn (99) + cis_scan_full. Every number on a
figure is grep-able from those files (no fabricated values). P3 uses documented 🟢 ARI controls
(somatic TP 0.135 / germline-het null 0.177 from master_draft:111; imprinting GNAS blind-ARI
0.473 from imprinting_poscontrol.out — NOT the stale schema "1.000", flagged for T-PROV).

Output: docs/paper_focus/04_figures/P{1..6}_*.png + catalog_overview.png
"""
import sys, csv, json, os
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod")
from scripts.lib.plot_setup import setup_plot_style
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
setup_plot_style()

GS = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2"
FIG = "/big7_disk/liaoyoyo2001/InterSubMod/docs/paper_focus/04_figures"
os.makedirs(FIG, exist_ok=True)
ACCENT = "#D97757"; INK = "#141413"; GREY = "#9b9b93"; BLUE = "#3b6ea5"

rows = list(csv.DictReader(open(f"{GS}/catalog/catalog_skeleton.tsv"), delimiter="\t"))
tagc = json.load(open(f"{GS}/catalog/catalog_tag_counts.json"))
venn = json.load(open(f"{GS}/catalog/discovery_venn.json"))


def fnum(x):
    try:
        v = float(x)
        return None if v != v else v
    except (ValueError, TypeError):
        return None


# ---------- P1: dbeta_max distribution ----------
def p1(ax):
    vals = [fnum(r["dbeta_max"]) for r in rows]
    vals = [v for v in vals if v is not None and v > 0]
    ax.hist(vals, bins=60, color=GREY, edgecolor="white", linewidth=0.3)
    ax.axvline(0.10, color=BLUE, ls="--", lw=1, label="medium ≥0.10")
    ax.axvline(0.20, color=ACCENT, ls="--", lw=1, label="large ≥0.20")
    ax.set_yscale("log")
    ax.set_xlabel("dbeta_max  (max |AlleleDelta|,|HPMergedDelta|)")
    ax.set_ylabel("loci (log)")
    ax.set_title("P1 · Δβ magnitude distribution", fontsize=10, loc="left")
    ax.legend(fontsize=7)


# ---------- P2: CramersV / reliability ----------
def p2(ax):
    cv = [fnum(r["cramersV_value"]) for r in rows]
    cv = [v for v in cv if v is not None and v > 0]
    ax.hist(cv, bins=50, color=ACCENT, alpha=0.8, edgecolor="white", linewidth=0.3)
    ax.axvline(0.10, color=INK, ls="--", lw=1, label="reliable gate ≥0.10")
    ax.set_yscale("log")
    ax.set_xlabel("CramersV (Cochran-gated; 0 = gated-out)")
    ax.set_ylabel("loci with CV>0 (log)")
    ax.set_title("P2 · CramersV (gated) + latent region", fontsize=10, loc="left")
    # latent annotation
    rc = {"reliable": 0, "latent": 0, "none": 0}
    for r in rows:
        rc[r["clustering_reliability"]] += 1
    ax.text(0.97, 0.95, f"reliable={rc['reliable']:,}\nlatent(PERMANOVA-recovered)={rc['latent']:,}\nnone={rc['none']:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round", fc="#FAF9F5", ec=GREY))
    ax.legend(fontsize=7, loc="lower left")


# ---------- P3: ARI controls (documented 🟢) ----------
def p3(ax):
    labels = ["somatic TP\n(median)", "germline-het\nnull", "imprinting +ctrl\nGNAS blind-ARI"]
    vals = [0.135, 0.177, 0.473]
    cols = [ACCENT, GREY, BLUE]
    bars = ax.bar(labels, vals, color=cols, edgecolor="white")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.01, f"{v:.3f}", ha="center", fontsize=8)
    ax.set_ylabel("blind-ARI")
    ax.set_ylim(0, 0.55)
    ax.set_title("P3 · ARI ruler: TP < null = germline-allelic", fontsize=10, loc="left")
    ax.text(0.5, 0.92, "somatic TP ARI (0.135) < germline-het null (0.177)\n→ clustering is germline-allelic, NOT somatic subclone (B2)",
            transform=ax.transAxes, ha="center", va="top", fontsize=6.5,
            bbox=dict(boxstyle="round", fc="#FFF6F2", ec=ACCENT))


# ---------- P4: within_d vs hp-axis_d scatter (cis vs copy) ----------
def p4(ax):
    cis = json.load(open(f"{GS}/cis_scan_full.json"))
    cp = {d["locus"]: d for d in json.load(open(f"{GS}/copy_partition_confirm.json"))["rows"]}
    xs, ys, cols, labs = [], [], [], []
    for loc, c in cp.items():
        dw = abs(c["d_within"]); dhp = abs(c["d_HP"])
        xs.append(dhp); ys.append(dw)
        is_cis = dw > 1.5 * abs(c["d_copy"])
        cols.append(ACCENT if is_cis else GREY)
        labs.append(loc)
    ax.scatter(xs, ys, c=cols, s=45, edgecolor="white", zorder=3)
    mx = max(xs + ys) * 1.1
    ax.plot([0, mx], [0, mx], color=INK, ls=":", lw=1, label="within = HP (cis line)")
    for x, y, loc in zip(xs, ys, labs):
        if loc in ("chr17:79991120", "chr13:32315128"):
            tag = "chr17/TBC1D16 (cis)" if "17" in loc else "BRCA2 (copy)"
            ax.annotate(tag, (x, y), fontsize=7, xytext=(5, 4), textcoords="offset points")
    ax.set_xlabel("HP-axis |Δβ| (d_HP)")
    ax.set_ylabel("within-subclone |Δβ| (d_within)")
    ax.set_title("P4 · within vs HP-axis: cis above line", fontsize=10, loc="left")
    ax.legend(fontsize=7)


# ---------- P5: per-TAG counts (R6 MAIN FIGURE) ----------
def p5(ax):
    order = ["G", "E", "C", "F", "B", "D", "A"]
    leg = tagc["tag_legend"]
    counts = [tagc["overall"].get(t, 0) for t in order]
    names = [f"{t}: {leg[t]}" for t in order]
    cmap = {"A": ACCENT, "B": "#c44", "C": BLUE, "D": "#d68", "E": "#6a9", "F": "#a86", "G": GREY}
    bars = ax.barh(names, counts, color=[cmap[t] for t in order], edgecolor="white")
    ax.set_xscale("log")
    for b, c in zip(bars, counts):
        ax.text(c * 1.15, b.get_y() + b.get_height() / 2, f"{c:,}", va="center", fontsize=7.5)
    ax.set_xlabel("loci (log)   ·   total = {:,}".format(tagc["n_loci"]))
    ax.set_title("P5 · per-TAG counts = Results R6 main", fontsize=10, loc="left")


# ---------- P6: coverage vs reliability ----------
def p6(ax):
    by = {"reliable": [], "latent": [], "none": []}
    for r in rows:
        c = fnum(r["coverage"])
        if c is not None:
            by[r["clustering_reliability"]].append(c)
    data = [by["reliable"], by["latent"], by["none"]]
    bp = ax.boxplot(data, labels=["reliable", "latent", "none"], showfliers=False, patch_artist=True)
    for patch, col in zip(bp["boxes"], [ACCENT, "#6a9", GREY]):
        patch.set_facecolor(col); patch.set_alpha(0.7)
    meds = [np.median(d) if d else 0 for d in data]
    for i, m in enumerate(meds, 1):
        ax.text(i, m, f"{m:.0f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_ylabel("coverage (NumReads)")
    ax.set_title("P6 · coverage vs reliability (not purely cov-driven)", fontsize=10, loc="left")


PANELS = [("P1_dbeta_distribution", p1), ("P2_cramersv_reliability", p2),
          ("P3_ari_controls", p3), ("P4_within_vs_hp_scatter", p4),
          ("P5_per_tag_counts_R6", p5), ("P6_coverage_vs_reliability", p6)]

# individual
for name, fn in PANELS:
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    fn(ax)
    fig.tight_layout()
    fig.savefig(f"{FIG}/{name}.png", dpi=140)
    plt.close(fig)

# overview 2x3
fig, axes = plt.subplots(2, 3, figsize=(15.5, 8))
for (name, fn), ax in zip(PANELS, axes.ravel()):
    fn(ax)
fig.suptitle("位點甲基分群 catalog — 觀察分佈圖 (forest view) · HCC1395+5 樣本 332,705 loci",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{FIG}/catalog_overview.png", dpi=140)
plt.close(fig)

print(f"[100] wrote 6 panels + overview to {FIG}")
for name, _ in PANELS:
    print("  ", f"{FIG}/{name}.png")
print("DONE")
