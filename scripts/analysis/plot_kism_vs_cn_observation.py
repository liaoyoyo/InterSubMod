#!/usr/bin/env python3
"""
plot_kism_vs_cn_observation.py
從 analyze_kism_vs_cn_perread.py 的 region_table_{tp,fp}.tsv + aggregate 產觀察圖panel。
English labels (避 CJK 字型問題)。輸出 figs/kism_vs_cn_observation.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/_assets")
FIGS = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/20260621_kism_vs_cn_perread/figs")
FIGS.mkdir(parents=True, exist_ok=True)

tp = pd.read_csv(A / "region_table_tp.tsv", sep="\t")
fp = pd.read_csv(A / "region_table_fp.tsv", sep="\t")
agg_tp = json.loads((A / "aggregate_tp.json").read_text())
agg_fp = json.loads((A / "aggregate_fp.json").read_text())

fig, ax = plt.subplots(2, 3, figsize=(16, 9))
fig.suptitle("ISM methylation clustering (k_ISM) vs CN/SV — HCC1395 genome-wide (TP 29,754 / FP 627)",
             fontsize=13, fontweight="bold")

# P1: structure rate TP vs FP
a = ax[0, 0]
a.bar(["TP", "FP"], [agg_tp["pct_structure_exists"] * 100, agg_fp["pct_structure_exists"] * 100],
      color=["#e74c3c", "#95a5a6"])
a.set_ylabel("% regions with significant\nmethylation cluster structure")
a.set_title("(1) Only ~10% regions have structure\n(PERMANOVA gate; TP 2.5x > FP)")
for i, v in enumerate([agg_tp["pct_structure_exists"] * 100, agg_fp["pct_structure_exists"] * 100]):
    a.text(i, v + 0.2, f"{v:.1f}%", ha="center", fontweight="bold")

# P2: k_eff distribution by k_CN (TP structured) -> Q1 no correlation
a = ax[0, 1]
s = tp[tp["structure_exists"] == True].dropna(subset=["k_CN"])
loh = s[s["loh"] == True]["k_eff"]; het = s[s["loh"] == False]["k_eff"]
bins = np.arange(1.5, 7.5, 1)
a.hist([loh, het], bins=bins, label=[f"LOH k_CN=1 (n={len(loh)})", f"het k_CN=2 (n={len(het)})"],
       color=["#9b59b6", "#2980b9"], density=True)
a.set_xlabel("k_eff (effective cluster count)"); a.set_ylabel("density")
a.set_title(f"(2) k_ISM NOT tracking k_CN\nSpearman rho={agg_tp['Q1_spearman_rho']:.3f}")
a.legend(fontsize=8)

# P3: alignment class breakdown (TP structured)
a = ax[0, 2]
cls = agg_tp["Q3_alignment_class_pct"]
order = ["unaligned", "candidate_beyond_CN", "aligned_somatic_allele_in_LOH", "ambiguous",
         "no_usable_allele_label", "CN_explained_HP", "aligned_somatic_allele"]
order = [c for c in order if c in cls]
vals = [cls[c] for c in order]
colors = ["#7f8c8d", "#e67e22", "#f1c40f", "#bdc3c7", "#bdc3c7", "#27ae60", "#16a085"]
a.barh(range(len(order)), vals, color=colors[:len(order)])
a.set_yticks(range(len(order))); a.set_yticklabels(order, fontsize=8)
a.invert_yaxis(); a.set_xlabel("% of structured regions")
a.set_title("(3) Cluster alignment to genetic axis\nCN_explained_HP = 0.64% only")
for i, v in enumerate(vals):
    a.text(v + 0.5, i, f"{v:.1f}%", va="center", fontsize=8)

# P4: cluster x HP CramerV vs cluster x ALT CramerV scatter
a = ax[1, 0]
ss = s.dropna(subset=["cluster_alt_cramersv"])
a.scatter(ss["cluster_hp_cramersv"].fillna(-0.05), ss["cluster_alt_cramersv"],
          alpha=0.3, s=14, color="#34495e")
a.axhline(0.7, color="green", ls="--", lw=1, label="aligned thr 0.7")
a.axvline(0.7, color="green", ls="--", lw=1)
a.set_xlabel("cluster x germline-HP Cramer's V\n(-0.05 = HP axis unusable)")
a.set_ylabel("cluster x somatic ALT/REF Cramer's V")
a.set_title("(4) Most clusters align to NEITHER\nallele axis (both V < 0.7)")
a.legend(fontsize=8)

# P5: het vs LOH outcome stacked
a = ax[1, 1]
het_s = s[s["loh"] == False]; loh_s = s[s["loh"] == True]
def frac(sub, cl): return (sub["alignment_class"] == cl).mean() * 100 if len(sub) else 0
cats = ["CN_explained_HP", "aligned_somatic_allele", "aligned_somatic_allele_in_LOH",
        "candidate_beyond_CN", "unaligned", "ambiguous", "no_usable_allele_label"]
het_v = [frac(het_s, c) for c in cats]; loh_v = [frac(loh_s, c) for c in cats]
cmap = {"CN_explained_HP": "#27ae60", "aligned_somatic_allele": "#16a085",
        "aligned_somatic_allele_in_LOH": "#f1c40f", "candidate_beyond_CN": "#e67e22",
        "unaligned": "#7f8c8d", "ambiguous": "#bdc3c7", "no_usable_allele_label": "#d0d3d4"}
bottom_h = bottom_l = 0
for c, h, l in zip(cats, het_v, loh_v):
    a.bar(["het", "LOH"], [h, l], bottom=[bottom_h, bottom_l], color=cmap[c], label=c)
    bottom_h += h; bottom_l += l
a.set_ylabel("% of structured regions"); a.set_title("(5) het: 85% unaligned\nLOH: 72% candidate_beyond_CN")
a.legend(fontsize=6, loc="upper right")

# P6: SV vs HP alignment rate on SV-informative subset
a = ax[1, 2]
a.bar(["SV axis", "HP axis"], [agg_tp.get("Q2_sv_aligned_pct", 0), agg_tp.get("Q2_hp_aligned_pct_same_subset", 0)],
      color=["#c0392b", "#2980b9"])
a.set_ylabel("% aligned (CramerV>=0.7 + sig)")
a.set_title(f"(6) SV axis weak\n(SV-informative n={agg_tp.get('Q2_sv_informative_regions',0)})")
for i, v in enumerate([agg_tp.get("Q2_sv_aligned_pct", 0), agg_tp.get("Q2_hp_aligned_pct_same_subset", 0)]):
    a.text(i, v + 0.1, f"{v:.1f}%", ha="center", fontweight="bold")

plt.tight_layout()
out = FIGS / "kism_vs_cn_observation.png"
plt.savefig(out, dpi=140, bbox_inches="tight")
plt.close()
print(f"saved {out}")
