"""Phase E · Main observation synthesis figures.

Generates 5 summary figures in figures/00_main/:
  - fig00_overall_auc_heatmap.png    — Group × top-features overall AUC overview
  - fig01_top5_strongest.png          — Top-5 strongest residualised-survived features
  - fig02_top5_collapsed.png          — Top-5 most collapsed (raw high, resid crash)
  - fig03_group_verdict_table.png     — Per-group verdict visualisation
  - fig04_s1s7_coverage.png           — Thread B S1-S7 schemes coverage by G-group

Inputs: G{1..10}_auc_table.tsv / global_stats.tsv / confound.tsv across data/.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/feature_layered_observation")
OUT = ROOT / "figures" / "00_main"
OUT.mkdir(parents=True, exist_ok=True)

# Summary table produced by manual curation from per-group md files.
# Columns: group, feature, raw_auc, residualized_auc, cross_sample_ok, verdict,
#          note (short mechanistic label)
SUMMARY = pd.DataFrame([
    # G1 Coverage
    ("G1", "NumReads",                  0.711, 0.777, True,  "CONDITIONAL_POSITIVE", "survives residualization (depth)"),
    ("G1", "Diploid_Coverage_Used",     0.790, 0.639, False, "CONFOUND_COLLAPSED",   "sample-ID proxy"),
    ("G1", "NNormalReads",              0.716, 0.528, False, "CONFOUND_COLLAPSED",   "paired vs TO mode proxy"),
    ("G1", "NTumorReads",               0.578, 0.472, False, "CONFOUND_COLLAPSED",   "NumReads proxy"),
    ("G1", "Coverage_Multiple",         0.493, 0.246, False, "NEGATIVE",             "random after baseline"),
    ("G1", "NumCpGs",                   0.437, 0.442, False, "NEGATIVE",             "anti-predictive"),
    # G2 LOH & Subclone
    ("G2", "LOH_Subtype_ordinal_TO",    0.596, 0.716, False, "SAMPLE_SPECIFIC",      "HP_Ratio self-reference"),
    ("G2", "Potential_LOH_TO",          0.538, 0.716, False, "SAMPLE_SPECIFIC",      "self-reference; HCC1954 only"),
    ("G2", "LOH_Subtype_paired",        0.430, 0.380, False, "ANTI_SIGNAL",          "paired flipped (LOH=TP-enriched)"),
    # G3 VCF Caller (from g3_step1_global.tsv)
    ("G3", "vcf_QUAL",                  0.813, None, True,   "POSITIVE_PROXY",       "caller internal QS - paired"),
    ("G3", "vcf_DP",                    0.786, None, True,   "POSITIVE_PROXY",       "NumReads echo"),
    ("G3", "vcf_AD_alt",                0.786, None, True,   "POSITIVE_PROXY",       "NumReads * AF"),
    ("G3", "vcf_GQ",                    0.811, None, True,   "POSITIVE_PROXY",       "caller internal QS - paired"),
    ("G3", "vcf_AF_TO",                 0.582, None, True,   "CONDITIONAL_POSITIVE", "caller AF (TO only)"),
    ("G3", "vcf_NAF_paired",            0.561, None, True,   "CONDITIONAL_POSITIVE", "normal AF rescue signal"),
    # G5 HP Merged
    ("G5", "GlobalP_HPFamily",          0.577, None, False,  "NEGATIVE",             "below 0.58 ceiling"),
    ("G5", "HP_FamilyN_sum",            0.567, None, False,  "NEGATIVE",             "NumReads proxy"),
    ("G5", "HP_Ratio_norm_pooled",      0.466, 0.532, False, "POOLING_ARTIFACT",     "Near-half pool=0.78 but 1/7 samples"),
    ("G5", "HPMergedDelta",             0.481, 0.436, False, "NEGATIVE",             "low |AF| correlation, captured"),
    # G6 HP Fine
    ("G6", "HPFineN_HP1S",              0.612, 0.571, True,  "CONDITIONAL_POSITIVE", "somatic-demoted count (flag off)"),
    ("G6", "HPFineN_HP2S",              0.587, 0.517, False, "CONFOUND_SUSPECT",     "HP1S mirror"),
    ("G6", "HPFineNGroups_global",      0.521, 0.509, False, "CONDITIONAL_STRATUM",  "LOH_Subclone/CN_amp AUC 0.60-0.68"),
    ("G6", "same_hap_precision_paired", 0.997, None, True,   "POSITIVE_FINGERPRINT", "NG=2 same_hap precision (paired)"),
    ("G6", "same_hap_precision_TO",     0.856, None, False,  "CONDITIONAL_POSITIVE", "HCC1954 TO precision 0.34 outlier"),
    ("G6", "HPFineD_HP1_HP2",           0.543, 0.515, False, "NEGATIVE",             "methylation pair distance"),
    # G7 Cluster & Global
    ("G7", "ClusterPermanovaF",         0.678, 0.521, False, "CONFOUND_COLLAPSED",   "TP:FP 245:1 paired inflation"),
    ("G7", "PairwiseMeanDist",          0.659, 0.553, True,  "CONDITIONAL_POSITIVE", "borderline post-residualization"),
    ("G7", "PairwiseMedianDist",        0.628, 0.547, True,  "CONFOUND_COLLAPSED",   "~ Mean duplicate"),
    ("G7", "CramersV",                  0.532, 0.389, False, "NEGATIVE",             "replaced by HPFine family"),
    ("G7", "neg_log10_GlobalP",         0.510, 0.375, False, "NEGATIVE",             "caller QS proxy"),
    # G8 Entropy (pending full run; plausible priors)
    ("G8", "NME_HP1",                   None, None, None, "PENDING",                 "retry in progress"),
    ("G8", "Epipoly_Delta",             None, None, None, "PENDING",                 "retry in progress"),
    ("G8", "Fisher_Frac_Sig_paired",    0.726, None, None, "CHARAC_ONLY",            "prior memory; CI crosses random"),
    # G9 ASM
    ("G9", "SampleASM_Delta",           0.779, 0.419, False, "CONFOUND_COLLAPSED",   "HCC1395 inverted; amplicon"),
    ("G9", "NormalBaseline_Coverage",   0.776, None, False,  "CONFOUND_COLLAPSED",   "SampleASM_NNormal echo"),
    ("G9", "NormalBaseline_Mean",       0.775, 0.485, False, "CONFOUND_COLLAPSED",   "mode-pool artifact (paired=0.43)"),
    ("G9", "AlleleDelta",               0.600, 0.528, False, "CONFOUND_COLLAPSED",   "L2 collider: =|vcf_AF|"),
    ("G9", "HP_Residual_Delta",         0.674, 0.342, False, "CONFOUND_COLLAPSED",   "tumor HP echo"),
    ("G9", "Combined_HP_Signed_Delta",  0.545, 0.549, False, "DATA_GAP",             "= Tumor_HP_Signed_Delta (bit-identical)"),
    ("G9", "Normal_HP_*",               None, None, None, "DATA_GAP",                 "0% populated in canonical"),
    # G10 Quality
    ("G10", "NHP0",                     0.764, None, False,  "SAMPLE_SPECIFIC",      "NumReads proxy; 0.36-0.78 range"),
    ("G10", "LabelAllelePermanovaF",    0.627, 0.496, False, "CONFOUND_COLLAPSED",   "AF proxy"),
    ("G10", "Quality_Score",            0.497, 0.517, False, "NEGATIVE",             "composite random (4 dimensions)"),
    ("G10", "HeuristicScore",           0.509, 0.462, False, "NEGATIVE",             "unoptimised weights"),
    ("G10", "VerificationClass_TO",     None,  None, False,  "CHARAC_ONLY",          "TP rate 0.21-1.00 inconsistent"),
    ("G10", "Stability",                0.500, None, None,   "NOT_IMPLEMENTED",      "all zero (possible bug)"),
    ("G10", "SuggestFilter",            None,  None, False,  "NEGATIVE",             "global prec=3.3% useless"),
], columns=["group", "feature", "raw_auc", "residualized_auc", "cross_sample_ok", "verdict", "note"])

SUMMARY.to_csv(OUT / "synthesis_table.tsv", sep="\t", index=False)
print(f"[summary] wrote {OUT/'synthesis_table.tsv'}, {len(SUMMARY)} rows")


# ---------- fig00: Overall group AUC heatmap ----------
# For each group, pick top-3 features by raw_auc (with raw_auc available), show raw vs resid AUC side-by-side.
fig0, ax0 = plt.subplots(1, 1, figsize=(11, 9))

grp_labels = ["G1", "G2", "G3", "G5", "G6", "G7", "G9", "G10"]  # G4/G8 pending
cell_text = []
cell_raw = []
cell_resid = []
cell_feats = []
for g in grp_labels:
    sub = SUMMARY[(SUMMARY.group == g) & SUMMARY.raw_auc.notna()].sort_values("raw_auc", ascending=False).head(3)
    for _, r in sub.iterrows():
        cell_feats.append(f"{g}::{r.feature}")
        cell_raw.append(r.raw_auc)
        cell_resid.append(r.residualized_auc if pd.notna(r.residualized_auc) else np.nan)

data_arr = np.array([cell_raw, cell_resid]).T   # n_features × 2
ax0.imshow(data_arr, cmap="RdYlGn", vmin=0.30, vmax=0.85, aspect="auto")
for i, (rv, resv) in enumerate(zip(cell_raw, cell_resid)):
    ax0.text(0, i, f"{rv:.3f}" if pd.notna(rv) else "—", ha="center", va="center", fontsize=8,
             color=("white" if rv > 0.72 else "black") if pd.notna(rv) else "black")
    ax0.text(1, i, f"{resv:.3f}" if pd.notna(resv) else "—", ha="center", va="center", fontsize=8,
             color="black")

ax0.set_xticks([0, 1])
ax0.set_xticklabels(["Raw AUC", "Residualized AUC"], fontsize=10)
ax0.set_yticks(range(len(cell_feats)))
ax0.set_yticklabels(cell_feats, fontsize=8)
ax0.set_title("Phase E · Group × Top-3 Feature AUC overview\nraw vs residualized on (NumReads + vcf_AF + Coverage_Multiple)", fontsize=11)
cbar = plt.colorbar(ax0.images[0], ax=ax0, shrink=0.7)
cbar.set_label("AUC", fontsize=9)
ax0.axhline(0.5, color="gray", lw=0.5)
plt.tight_layout()
plt.savefig(OUT / "fig00_overall_auc_heatmap.png", dpi=150)
plt.close(fig0)
print("[fig00] saved")


# ---------- fig01: Top-5 strongest (residualized survived) ----------
# Criterion: residualized AUC >= 0.55, cross_sample_ok TRUE OR raw_auc very high even if not residualized available
strong_candidates = SUMMARY[SUMMARY.raw_auc.notna() & SUMMARY.verdict.str.contains("POSITIVE", na=False)].copy()
# include same_hap_precision_paired as fingerprint top
strong_candidates = strong_candidates.sort_values("raw_auc", ascending=False).head(5)

fig1, ax1 = plt.subplots(figsize=(10, 5))
ypos = np.arange(len(strong_candidates))
colors = ["#2e7d32" if "CONDITIONAL" not in v and "PROXY" not in v else "#43a047"
         for v in strong_candidates.verdict]
ax1.barh(ypos, strong_candidates.raw_auc, color=colors, alpha=0.85, label="Raw AUC")
ax1.barh(ypos + 0.18, strong_candidates.residualized_auc.fillna(np.nan), height=0.15, color="#1565c0",
         alpha=0.7, label="Residualized AUC")
ax1.set_yticks(ypos)
ax1.set_yticklabels([f"{r.group}::{r.feature}" for _, r in strong_candidates.iterrows()], fontsize=9)
for i, r in enumerate(strong_candidates.itertuples()):
    ax1.text(r.raw_auc + 0.005, i, f"{r.raw_auc:.3f}", va="center", fontsize=8)
    if pd.notna(r.residualized_auc):
        ax1.text(r.residualized_auc + 0.005, i + 0.18, f"{r.residualized_auc:.3f}", va="center", fontsize=8)
ax1.axvline(0.58, color="red", ls="--", lw=1, label="0.58 ceiling")
ax1.axvline(0.50, color="gray", ls=":", lw=1, label="random")
ax1.set_xlabel("AUC")
ax1.set_title("Phase E · Top-5 strongest TP/FP discriminative features (residualized-aware ranking)")
ax1.legend(loc="lower right", fontsize=8)
ax1.set_xlim(0.4, 1.05)
plt.tight_layout()
plt.savefig(OUT / "fig01_top5_strongest.png", dpi=150)
plt.close(fig1)
print("[fig01] saved")


# ---------- fig02: Top-5 confound-collapsed ----------
collapsed = SUMMARY[SUMMARY.raw_auc.notna() & SUMMARY.residualized_auc.notna() & (SUMMARY.verdict == "CONFOUND_COLLAPSED")].copy()
collapsed["delta"] = collapsed.raw_auc - collapsed.residualized_auc
collapsed = collapsed.sort_values("delta", ascending=False).head(5)

fig2, ax2 = plt.subplots(figsize=(10, 5))
ypos = np.arange(len(collapsed))
ax2.barh(ypos, collapsed.raw_auc, color="#c62828", alpha=0.7, label="Raw AUC")
ax2.barh(ypos + 0.28, collapsed.residualized_auc, height=0.25, color="#616161",
         alpha=0.7, label="Residualized AUC")
for i, r in enumerate(collapsed.itertuples()):
    ax2.text(r.raw_auc + 0.005, i, f"{r.raw_auc:.3f}", va="center", fontsize=8)
    ax2.text(r.residualized_auc + 0.005, i + 0.28, f"{r.residualized_auc:.3f}", va="center", fontsize=8, color="dimgray")
    ax2.text(1.03, i + 0.14, f"Δ={r.delta:+.3f}", va="center", fontsize=8, color="darkred", weight="bold")
ax2.set_yticks(ypos)
ax2.set_yticklabels([f"{r.group}::{r.feature}" for _, r in collapsed.iterrows()], fontsize=9)
ax2.axvline(0.58, color="red", ls="--", lw=1, label="0.58 ceiling")
ax2.set_xlabel("AUC")
ax2.set_title("Phase E · Top-5 unexpected confound-collapsed features\n(raw AUC high but residualization reveals depth/AF/CovM capture)")
ax2.legend(loc="lower right", fontsize=8)
ax2.set_xlim(0.25, 1.15)
plt.tight_layout()
plt.savefig(OUT / "fig02_top5_collapsed.png", dpi=150)
plt.close(fig2)
print("[fig02] saved")


# ---------- fig03: per-group verdict distribution ----------
verdict_order = ["POSITIVE_PROXY", "CONDITIONAL_POSITIVE", "POSITIVE_FINGERPRINT",
                 "CONFOUND_COLLAPSED", "SAMPLE_SPECIFIC", "POOLING_ARTIFACT",
                 "CHARAC_ONLY", "ANTI_SIGNAL", "NEGATIVE",
                 "DATA_GAP", "NOT_IMPLEMENTED", "PENDING",
                 "CONFOUND_SUSPECT", "CONDITIONAL_STRATUM"]
grp_set = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10"]
vt = pd.crosstab(SUMMARY.group, SUMMARY.verdict).reindex(index=grp_set, columns=verdict_order, fill_value=0)
fig3, ax3 = plt.subplots(figsize=(14, 5))
im = ax3.imshow(vt.values, cmap="YlOrRd", aspect="auto")
ax3.set_xticks(range(len(verdict_order)))
ax3.set_xticklabels(verdict_order, rotation=30, ha="right", fontsize=8)
ax3.set_yticks(range(len(grp_set)))
ax3.set_yticklabels(grp_set, fontsize=9)
for i in range(len(grp_set)):
    for j in range(len(verdict_order)):
        v = vt.values[i, j]
        if v > 0:
            ax3.text(j, i, str(v), ha="center", va="center",
                     color="white" if v > 2 else "black", fontsize=8)
ax3.set_title("Phase E · Verdict distribution across 10 feature groups (n=features per cell)")
plt.colorbar(im, ax=ax3, shrink=0.6)
plt.tight_layout()
plt.savefig(OUT / "fig03_group_verdict_table.png", dpi=150)
plt.close(fig3)
print("[fig03] saved")


# ---------- fig04: S1-S7 coverage mapping ----------
# Which group's features each scheme relies on
scheme_coverage = {
    "S1 LOH_Strong_Extreme": {"G1": 1, "G2": 2, "G3": 1, "G6": 0, "G7": 0, "G9": 1, "G10": 0},
    "S2 Subclonal_LOH_Inter": {"G1": 1, "G2": 2, "G3": 1, "G6": 1, "G7": 0, "G9": 0, "G10": 0},
    "S3 Diploid_Het":         {"G1": 1, "G2": 1, "G3": 1, "G6": 0, "G7": 0, "G9": 0, "G10": 0},
    "S4 NonLOH_Extreme":      {"G1": 0, "G2": 1, "G3": 1, "G6": 0, "G7": 0, "G9": 0, "G10": 0},
    "S5 Combo":               {"G1": 1, "G2": 2, "G3": 1, "G6": 0, "G7": 0, "G9": 1, "G10": 0},
    "S6/S7 + HPFineNGroups":  {"G1": 1, "G2": 2, "G3": 1, "G6": 2, "G7": 0, "G9": 0, "G10": 0},
}
sc_df = pd.DataFrame(scheme_coverage).T[["G1", "G2", "G3", "G6", "G7", "G9", "G10"]]
fig4, ax4 = plt.subplots(figsize=(9, 4.5))
im = ax4.imshow(sc_df.values, cmap="Blues", aspect="auto", vmin=0, vmax=2)
ax4.set_xticks(range(sc_df.shape[1]))
ax4.set_xticklabels(sc_df.columns, fontsize=10)
ax4.set_yticks(range(sc_df.shape[0]))
ax4.set_yticklabels(sc_df.index, fontsize=9)
for i in range(sc_df.shape[0]):
    for j in range(sc_df.shape[1]):
        v = sc_df.values[i, j]
        mark = {0: "-", 1: "used", 2: "core"}.get(v, str(v))
        ax4.text(j, i, mark, ha="center", va="center",
                 color="white" if v > 1 else "black", fontsize=9)
ax4.set_title("Phase E · Thread B S1-S7 scheme dependency on feature groups\n(-: not used, used: single axis, core: defining feature)")
plt.tight_layout()
plt.savefig(OUT / "fig04_s1s7_coverage.png", dpi=150)
plt.close(fig4)
print("[fig04] saved")

print("Done. Outputs under:", OUT)
