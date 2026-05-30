"""
F5_real_ism_heatmap.py — F5 figure with REAL HCC1395 ISM data (replaces synthetic)

Replaces the previous synthetic (rng.beta placeholder) F5 produced by
gen_concept_figures.py:fig5_ism_readcpg_heatmap with a heatmap of the REAL
read x CpG methylation-probability matrix from one HCC1395 ISM TP region.

Real data source (provenance):
    chr2:18,072,546 G->C somatic TP SNV (QS 0.9887), V3F_on_tp pipeline,
    region window chr2:18,067,546-18,077,546
    /big7_disk/.../filtered_snv_tp/chr2/chr2_18072546/chr2_18067546_18077546/
        methylation/methylation.csv  (43 reads x 53 CpG, prob 0..1, NA = no call)
        reads/reads.tsv              (read ids + HP tag + ALT/REF support)
        clustering/linkage_matrix.csv (scipy Z: cluster_i,cluster_j,distance,new_id,size)

Read ordering:
    Prefer the precomputed hierarchical-clustering dendrogram leaf order derived
    from clustering/linkage_matrix.csv (scipy dendrogram). Fall back to ordering
    by HP tag if the linkage matrix is missing.

Output:
    research/hku_collaboration/figures/F5_ism_readcpg_heatmap.png  (OVERWRITE)
    DPI=150, size comparable to original (13 x 7).

Repro:
    cd /big7_disk/liaoyoyo2001/InterSubMod
    python3 research/hku_collaboration/scripts/F5_real_ism_heatmap.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy.cluster.hierarchy import dendrogram  # noqa: E402

# --- CJK-safe font chain (mandatory: avoid tofu, avoid emoji in figure text) ---
# 2-font chain verified per-glyph fallback on matplotlib 3.6.2:
#   DejaVu Sans for Latin/numbers, Droid Sans Fallback for CJK.
_FONT_CHAIN = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["font.family"] = _FONT_CHAIN
matplotlib.rcParams["font.sans-serif"] = _FONT_CHAIN
matplotlib.rcParams["axes.unicode_minus"] = False

PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
REGION_DIR = (
    PROJECT_ROOT
    / "research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance"
    / "V3F_on_tp/filtered_snv_tp/chr2/chr2_18072546/chr2_18067546_18077546"
)
METH_CSV = REGION_DIR / "methylation" / "methylation.csv"
READS_TSV = REGION_DIR / "reads" / "reads.tsv"
LINKAGE_CSV = REGION_DIR / "clustering" / "linkage_matrix.csv"
LEAF_ORDER_TXT = REGION_DIR / "clustering" / "leaf_order.txt"

FIG_DIR = PROJECT_ROOT / "research" / "hku_collaboration" / "figures"
OUT_PNG = FIG_DIR / "F5_ism_readcpg_heatmap.png"

# Colors (consistent with gen_concept_figures.py palette)
COLOR_HP1 = "#5C7CFA"
COLOR_HP2 = "#FA5252"
COLOR_HP3 = "#868E96"
COLOR_HP0 = "#CED4DA"  # untagged (hp=0)
COLOR_ALT = "#FF6B6B"  # red
COLOR_REF = "#4DABF7"  # blue

SOURCE_NOTE = (
    "Source: HCC1395 ISM real data · chr2:18,072,546 region · "
    "F5_real_ism_heatmap.py · 2026-05-29"
)


def hp_to_label(hp_val) -> str:
    """Map integer HP tag in reads.tsv to a display label.

    reads.tsv encodes phasing as an integer: 0 = untagged (longphase could not
    assign), 1 = HP1 lineage, 2 = HP2 lineage. Higher subdivision tags
    (HP1-1/HP2-1/HP3) are not present in this region's reads.tsv.
    """
    try:
        v = int(hp_val)
    except (TypeError, ValueError):
        return "none"
    return {0: "none", 1: "HP1", 2: "HP2"}.get(v, f"HP{v}")


def load_methylation(path: Path):
    """Return (matrix float [n_reads x n_cpg] with NaN, read_ids list, cpg_pos list)."""
    df = pd.read_csv(path, na_values=["NA", "N", ""])
    read_ids = df["read_id"].astype(str).tolist()
    cpg_cols = [c for c in df.columns if c != "read_id"]
    mat = df[cpg_cols].to_numpy(dtype=float)
    return mat, read_ids, cpg_cols


def _rebuild_scipy_Z(ldf, n_reads):
    """The stored linkage_matrix.csv uses original observation indices that scipy's
    dendrogram rejects (a leaf id can appear as both leaf and a merged-cluster id).
    Remap each row's new_cluster_id to scipy's canonical numbering (n + step) while
    translating cluster references through that mapping, yielding a valid Z."""
    id_remap = {}  # stored new_cluster_id -> scipy cluster id
    rows = []
    for step, (_, r) in enumerate(ldf.iterrows()):
        ci, cj = int(r["cluster_i"]), int(r["cluster_j"])
        a = id_remap.get(ci, ci)
        b = id_remap.get(cj, cj)
        rows.append([a, b, float(r["distance"]), int(r["size"])])
        id_remap[int(r["new_cluster_id"])] = n_reads + step
    return np.array(rows, dtype=float)


def get_read_order(reads_df, read_ids, hp_labels):
    """Return (order_list, mode_str).

    Priority 1: precomputed dendrogram leaf_order.txt (read_names in cluster order).
    Priority 2: scipy dendrogram rebuilt from linkage_matrix.csv.
    Priority 3: group by HP tag.
    """
    n_reads = len(read_ids)
    # read_id (row index in methylation.csv) -> read_name(UUID)
    rid_to_name = {str(r): str(reads_df.set_index("read_id").loc[str(r), "read_name"])
                   for r in read_ids if str(r) in reads_df["read_id"].astype(str).values}
    name_to_pos = {name: pos for pos, rid in enumerate(read_ids)
                   for name in [rid_to_name.get(str(rid))] if name is not None}

    # Priority 1: leaf_order.txt
    if LEAF_ORDER_TXT.exists():
        names = [ln.strip() for ln in LEAF_ORDER_TXT.read_text().splitlines() if ln.strip()]
        order = [name_to_pos[n] for n in names if n in name_to_pos]
        if sorted(order) == list(range(n_reads)):
            return order, "dendrogram leaf_order.txt (precomputed clustering)"

    # Priority 2: rebuild scipy Z from linkage_matrix.csv
    if LINKAGE_CSV.exists():
        try:
            ldf = pd.read_csv(LINKAGE_CSV, sep=None, engine="python")
            Z = _rebuild_scipy_Z(ldf, n_reads)
            if Z.shape[0] == n_reads - 1:
                dd = dendrogram(Z, no_plot=True)
                leaves = [int(x) for x in dd["leaves"]]
                if sorted(leaves) == list(range(n_reads)):
                    return leaves, "dendrogram rebuilt from linkage_matrix.csv"
        except Exception as exc:  # noqa: BLE001
            print(f"[warn] linkage rebuild failed ({exc}); falling back")

    # Priority 3: order by HP tag (HP1, HP2, none), stable within group
    hp_rank = {"HP1": 0, "HP2": 1, "none": 2}
    order = sorted(range(n_reads), key=lambda i: (hp_rank.get(hp_labels[i], 9), i))
    return order, "HP-tag grouping (clustering order unavailable)"


def main():
    mat, read_ids, cpg_cols = load_methylation(METH_CSV)
    n_reads, n_cpg = mat.shape
    print(f"[load] methylation matrix: {n_reads} reads x {n_cpg} CpG (NaN = no methyl call)")

    # HP + allele labels from reads.tsv, indexed by row order of methylation.csv.
    reads_df = pd.read_csv(READS_TSV, sep="\t")
    reads_df["read_id"] = reads_df["read_id"].astype(str)
    rmap = reads_df.set_index("read_id")
    hp_labels = [hp_to_label(rmap.loc[rid, "hp"]) if rid in rmap.index else "none" for rid in read_ids]
    alt_labels = [
        str(rmap.loc[rid, "alt_support"]).upper() if rid in rmap.index else "NA" for rid in read_ids
    ]

    order, order_mode = get_read_order(reads_df, read_ids, hp_labels)
    print(f"[order] read order = {order_mode}")

    mat_o = mat[order, :]
    hp_o = [hp_labels[i] for i in order]
    alt_o = [alt_labels[i] for i in order]
    rid_o = [read_ids[i] for i in order]

    # --- figure layout: HP bar | ALT/REF bar | heatmap | colorbar ---
    fig = plt.figure(figsize=(13, 7), dpi=150)
    gs = fig.add_gridspec(
        1, 4, width_ratios=[0.28, 0.28, 6, 0.3], wspace=0.04,
        left=0.10, right=0.95, top=0.85, bottom=0.16,
    )

    hp_color_map = {"HP1": COLOR_HP1, "HP2": COLOR_HP2, "HP3": COLOR_HP3, "none": COLOR_HP0}
    alt_color_map = {"ALT": COLOR_ALT, "REF": COLOR_REF, "NA": "#E9ECEF"}

    # 1. HP sidebar
    ax_hp = fig.add_subplot(gs[0])
    for i, hp in enumerate(hp_o):
        ax_hp.add_patch(Rectangle((0, i), 1, 1, facecolor=hp_color_map.get(hp, COLOR_HP0),
                                  edgecolor="white", linewidth=0.4))
    ax_hp.set_xlim(0, 1)
    ax_hp.set_ylim(0, n_reads)
    ax_hp.invert_yaxis()
    ax_hp.set_xticks([])
    ax_hp.set_yticks([])
    ax_hp.set_title("HP", fontsize=9, pad=4)
    for sp in ax_hp.spines.values():
        sp.set_visible(False)

    # 2. ALT/REF sidebar (the dominant real signal: Cramér's V 0.87 vs HP 0.78)
    ax_al = fig.add_subplot(gs[1])
    for i, al in enumerate(alt_o):
        ax_al.add_patch(Rectangle((0, i), 1, 1, facecolor=alt_color_map.get(al, "#E9ECEF"),
                                  edgecolor="white", linewidth=0.4))
    ax_al.set_xlim(0, 1)
    ax_al.set_ylim(0, n_reads)
    ax_al.invert_yaxis()
    ax_al.set_xticks([])
    ax_al.set_yticks([])
    ax_al.set_title("Allele", fontsize=9, pad=4)
    for sp in ax_al.spines.values():
        sp.set_visible(False)

    # 3. main heatmap (NaN cells masked -> shown as light gray)
    ax_main = fig.add_subplot(gs[2])
    sns.heatmap(
        mat_o, ax=ax_main, cmap="RdYlBu_r", vmin=0, vmax=1, cbar=False,
        mask=np.isnan(mat_o),
        xticklabels=[c for c in cpg_cols], yticklabels=rid_o,
        linewidths=0.2, linecolor="#F8F9FA",
    )
    ax_main.set_facecolor("#E9ECEF")  # masked NaN cells
    ax_main.set_xlabel("CpG genomic position (5' -> 3')", fontsize=10)
    ax_main.set_ylabel("")
    # CpG positions are dense; show every 4th tick to stay readable
    xt = np.arange(n_cpg) + 0.5
    ax_main.set_xticks(xt[::4])
    ax_main.set_xticklabels([cpg_cols[i] for i in range(0, n_cpg, 4)], rotation=90, fontsize=6.5)
    ax_main.set_yticks(np.arange(n_reads) + 0.5)
    ax_main.set_yticklabels([f"read_{r}" for r in rid_o], fontsize=6, rotation=0)

    # 4. colorbar
    ax_cb = fig.add_subplot(gs[3])
    sm = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax_cb)
    cb.set_label("methylation prob P(5mC)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # sidebar legends
    hp_present = [h for h in ["HP1", "HP2", "HP3", "none"] if h in set(hp_o)]
    legend_hp = [mpatches.Patch(facecolor=hp_color_map[h], edgecolor="#ADB5BD",
                                label=f"HP={h}") for h in hp_present]
    legend_al = [mpatches.Patch(facecolor=alt_color_map[a], edgecolor="#ADB5BD", label=a)
                 for a in ["ALT", "REF"] if a in set(alt_o)]
    legend_na = [mpatches.Patch(facecolor="#E9ECEF", edgecolor="#ADB5BD", label="no methyl call (NA)")]
    fig.legend(handles=legend_hp + legend_al + legend_na, loc="lower center",
               ncol=len(legend_hp) + len(legend_al) + 1, fontsize=8.5,
               frameon=True, bbox_to_anchor=(0.5, 0.005))

    fig.suptitle(
        "F5 · HCC1395 ISM — chr2:18,072,546 region "
        f"({n_reads} reads × {n_cpg} CpG, real methylation matrix)",
        fontsize=13.5, fontweight="bold", y=0.965,
    )
    fig.text(
        0.525, 0.905,
        "Real read×CpG 5mC-probability matrix; reads ordered by hierarchical "
        "clustering dendrogram (NHD distance). HP / Allele sidebars from reads.tsv.",
        ha="center", fontsize=9, color="#495057", style="italic",
    )
    fig.text(0.99, 0.005, SOURCE_NOTE, fontsize=7, color="#ADB5BD", ha="right", va="bottom")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F5] wrote {OUT_PNG}")
    print(f"[F5] dims={n_reads}x{n_cpg}  order={order_mode}  "
          f"HP sidebar={sorted(set(hp_o))}  Allele sidebar={sorted(set(alt_o))}")


if __name__ == "__main__":
    main()
