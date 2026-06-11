#!/usr/bin/env python3
"""
51_context_tpfp_clustering.py
=============================
C2 - USER HYPOTHESES: read-clustering x genomic context vs TP/FP (HCC1395, single sample).

Input : genome_survey_v2/cn_confound/master_o2_error.tsv  (596 BAM-pass HP-axis loci, with blind_ari)
Output: genome_survey_v2/cn_confound/c2_context_tpfp_clustering.json
        figures/cn_confound/c2_context_tpfp_clustering.png

Pre-registered user hypotheses (tested explicitly, verdict per each):
  H-ctx1: 'LOH 內有甲基分群(high blind_ari) 的可能是異常拷貝 FP'
          -> within {LOH & abnormal-CN (cn_class!=neutral) & high-clustering},
             is TP-rate LOWER than overall (i.e. FP-enriched)? Fisher exact.
  H-ctx2: '非LOH 但 coverage 高 的可能是 TP'
          -> within {nonLOH & high-coverage (n_paired_cpg top quartile) & high-clustering},
             is TP-rate HIGHER than overall? Fisher exact.

GUARDS (auc-confound-guard): any 'discriminative' context is a CONTEXT effect, NOT a
methylation-feature classifier. We report TP-rate differences with Wilson CI + Fisher exact
odds ratio + p, and explicitly DO NOT fit / claim a classifier or AUC.

POWER: FP n=43 total; 41 of 43 FP are in LOH, only 2 FP in nonLOH. Many context cells have
0 FP -> TP-rate=100% with no discriminative power. Per-cell n reported; cells with FP<5
flagged underpowered; cells with FP=0 flagged 'no-FP / uninformative'.

Each of the 596 rows is a UNIQUE (chrom,somatic_pos) -> one HP-axis record per locus in this
BAM-pass subset -> NO within-subset pseudoreplication (verified: 0 positions with >1 row).
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.proportion import proportion_confint
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# ---------------------------------------------------------------- paths
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
IN_TSV = f"{ROOT}/genome_survey_v2/cn_confound/master_o2_error.tsv"
OUT_JSON = f"{ROOT}/genome_survey_v2/cn_confound/c2_context_tpfp_clustering.json"
OUT_PNG = f"{ROOT}/figures/cn_confound/c2_context_tpfp_clustering.png"
CJK = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
try:
    cjk_fp = FontProperties(fname=CJK)
except Exception:
    cjk_fp = None

# ---------------------------------------------------------------- load
df = pd.read_csv(IN_TSV, sep="\t")
N = len(df)
n_tp = int((df["is_tp"] == 1).sum())
n_fp = int((df["is_tp"] == 0).sum())
overall_tp_rate = n_tp / N

# integrity: unique positions (no pseudoreplication within this subset)
n_unique_pos = df.groupby(["chrom", "somatic_pos"]).ngroups
pseudoreplicated = bool(n_unique_pos != N)

# ---------------------------------------------------------------- derived context vars
# clustering threshold: pre-registered as (>= median) OR (>= 0.3); use median (data-driven)
# AND report a fixed-0.3 variant for robustness.
ari_med = float(df["blind_ari"].median())
df["clust_hi_med"] = (df["blind_ari"] >= ari_med).astype(int)
df["clust_hi_03"] = (df["blind_ari"] >= 0.3).astype(int)
# primary clustering def used in cross-tab = median split (balances n); 0.3 reported alongside
df["clust_hi"] = df["clust_hi_med"]

# coverage: top quartile of n_paired_cpg
cov_q75 = float(df["n_paired_cpg"].quantile(0.75))
df["cov_hi"] = (df["n_paired_cpg"] >= cov_q75).astype(int)

# abnormal CN
df["abnormal_cn"] = (df["cn_class"] != "neutral").astype(int)


# ---------------------------------------------------------------- helpers
def wilson_ci(k, n):
    if n == 0:
        return (None, None)
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return (float(lo), float(hi))


def cell_stats(mask, label):
    """TP-rate + Wilson CI + per-cell n for an arbitrary boolean mask."""
    sub = df[mask]
    n = int(len(sub))
    tp = int((sub["is_tp"] == 1).sum())
    fp = int((sub["is_tp"] == 0).sum())
    rate = (tp / n) if n else None
    lo, hi = wilson_ci(tp, n) if n else (None, None)
    return {
        "label": label,
        "n": n,
        "tp": tp,
        "fp": fp,
        "tp_rate": rate,
        "tp_rate_ci95": [lo, hi],
        "underpowered_fp_lt5": bool(fp < 5),
        "no_fp_uninformative": bool(fp == 0),
    }


def fisher_cell_vs_rest(mask, label, direction):
    """
    2x2 Fisher: cell vs the REST of the cohort.
                       in-cell   rest
        TP (is_tp=1)     a        b
        FP (is_tp=0)     c        d
    OR>1 => cell is TP-enriched (relative to rest); OR<1 => cell is FP-enriched.
    direction: 'lower' (H-ctx1, expect FP-enriched OR<1) or 'higher' (H-ctx2, expect TP-enriched OR>1).
    """
    incell = df[mask]
    rest = df[~mask]
    a = int((incell["is_tp"] == 1).sum())
    c = int((incell["is_tp"] == 0).sum())
    b = int((rest["is_tp"] == 1).sum())
    d = int((rest["is_tp"] == 0).sum())
    table = [[a, b], [c, d]]
    # two-sided p (report), plus one-sided p in the hypothesized direction
    orr, p_two = fisher_exact(table, alternative="two-sided")
    if direction == "lower":  # FP-enriched in cell => OR<1 => 'less'
        _, p_one = fisher_exact(table, alternative="less")
    else:  # 'higher' TP-enriched => OR>1 => 'greater'
        _, p_one = fisher_exact(table, alternative="greater")
    cell_n = a + c
    cell_rate = (a / cell_n) if cell_n else None
    rest_n = b + d
    rest_rate = (b / rest_n) if rest_n else None
    lo, hi = wilson_ci(a, cell_n) if cell_n else (None, None)
    return {
        "label": label,
        "table_2x2": {"cell_TP": a, "cell_FP": c, "rest_TP": b, "rest_FP": d},
        "cell_n": cell_n,
        "cell_tp_rate": cell_rate,
        "cell_tp_rate_ci95": [lo, hi],
        "rest_n": rest_n,
        "rest_tp_rate": rest_rate,
        "odds_ratio": (float(orr) if np.isfinite(orr) else None),
        "p_two_sided": float(p_two),
        "p_one_sided_directional": float(p_one),
        "direction_tested": direction,
        "cell_fp_n": c,
        "underpowered_fp_lt5": bool(c < 5),
        "no_fp_uninformative": bool(c == 0),
    }


# ---------------------------------------------------------------- H-ctx1
# LOH & abnormal-CN & high-clustering  -> expect TP-rate LOWER (FP-enriched)
m_ctx1 = (df["loh_status"] == "LOH") & (df["abnormal_cn"] == 1) & (df["clust_hi"] == 1)
hctx1 = fisher_cell_vs_rest(m_ctx1, "LOH & abnormal-CN & high-clustering(>=median)", "lower")
# robustness with 0.3 clustering threshold
m_ctx1_03 = (df["loh_status"] == "LOH") & (df["abnormal_cn"] == 1) & (df["clust_hi_03"] == 1)
hctx1_03 = fisher_cell_vs_rest(m_ctx1_03, "LOH & abnormal-CN & clustering>=0.3", "lower")

# verdict logic
def verdict_lower(res):
    # supported = TP-rate lower (OR<1) AND directional p<0.05 AND enough FP power
    if res["cell_fp_n"] == 0:
        return "underpowered"  # no FP in cell => cannot show FP-enrichment
    if res["cell_n"] < 10:
        return "underpowered"
    enough = res["cell_fp_n"] >= 5
    eff = (res["odds_ratio"] is not None and res["odds_ratio"] < 1.0)
    sig = res["p_one_sided_directional"] < 0.05
    if eff and sig and enough:
        return "supported"
    if eff and sig and not enough:
        return "supported_but_underpowered"
    return "refuted"


def verdict_higher(res):
    if res["cell_n"] < 10:
        return "underpowered"
    # for TP-enrichment we also need the comparison to have FP somewhere; if rest has FP and
    # cell has 0 FP that *can* be TP-enrichment, but with cell FP=0 and tiny rest-FP power is poor.
    eff = (res["odds_ratio"] is not None and res["odds_ratio"] > 1.0) or (res["cell_fp_n"] == 0 and res["cell_tp_rate"] == 1.0)
    sig = res["p_one_sided_directional"] < 0.05
    # the discriminating FP pool for nonLOH is only 2 -> structurally underpowered
    if res["cell_fp_n"] == 0 and res["rest_n"] is not None:
        # whether this is informative depends on total FP available in the comparison universe
        pass
    if eff and sig:
        return "supported"
    if eff and not sig:
        return "trend_not_sig"
    return "refuted"


v_ctx1 = verdict_lower(hctx1)

# ---------------------------------------------------------------- H-ctx2
# nonLOH & high-coverage & high-clustering -> expect TP-rate HIGHER
m_ctx2 = (df["loh_status"] == "nonLOH") & (df["cov_hi"] == 1) & (df["clust_hi"] == 1)
hctx2 = fisher_cell_vs_rest(m_ctx2, "nonLOH & high-cov(top-quartile) & high-clustering(>=median)", "higher")
m_ctx2_03 = (df["loh_status"] == "nonLOH") & (df["cov_hi"] == 1) & (df["clust_hi_03"] == 1)
hctx2_03 = fisher_cell_vs_rest(m_ctx2_03, "nonLOH & high-cov & clustering>=0.3", "higher")
v_ctx2 = verdict_higher(hctx2)

# nonLOH FP availability (structural power ceiling for H-ctx2)
nonloh_fp_total = int(((df["loh_status"] == "nonLOH") & (df["is_tp"] == 0)).sum())

# ---------------------------------------------------------------- FULL see-tree-and-forest cross-tab
# TP-rate by {LOH x cn_class x clustering(hi/lo) x cov(hi/lo)} with per-cell n
crosstab_cells = []
for loh in ["LOH", "nonLOH"]:
    for cn in ["gain", "neutral", "loss"]:
        for clu in [1, 0]:
            for cov in [1, 0]:
                m = ((df["loh_status"] == loh) & (df["cn_class"] == cn)
                     & (df["clust_hi"] == clu) & (df["cov_hi"] == cov))
                if m.sum() == 0:
                    continue
                st = cell_stats(m, f"{loh}|{cn}|clust_{'hi' if clu else 'lo'}|cov_{'hi' if cov else 'lo'}")
                st.update({"loh": loh, "cn_class": cn, "clust_hi": int(clu), "cov_hi": int(cov)})
                crosstab_cells.append(st)

# also a compact marginal grid: LOH x clustering(hi/lo) (most interpretable for the figure)
marginal_grid = {}
for loh in ["LOH", "nonLOH"]:
    for clu in [1, 0]:
        m = (df["loh_status"] == loh) & (df["clust_hi"] == clu)
        st = cell_stats(m, f"{loh}|clust_{'hi' if clu else 'lo'}")
        marginal_grid[f"{loh}_clust_{'hi' if clu else 'lo'}"] = st

# ---------------------------------------------------------------- rank trees by blind_ari within is_tp x context
def top_loci(mask, k=8):
    sub = df[mask].sort_values("blind_ari", ascending=False).head(k)
    out = []
    for _, r in sub.iterrows():
        out.append({
            "chrom": r["chrom"], "pos": int(r["somatic_pos"]), "axis": r["axis"],
            "is_tp": int(r["is_tp"]), "loh": r["loh_status"], "cn_class": r["cn_class"],
            "median_cn": float(r["median_cn"]), "n_paired_cpg": int(r["n_paired_cpg"]),
            "blind_ari": float(r["blind_ari"]), "abs_delta": float(r["abs_delta"]),
        })
    return out

tree_ranks = {
    "TP_LOH_top_blind_ari": top_loci((df["is_tp"] == 1) & (df["loh_status"] == "LOH")),
    "TP_nonLOH_top_blind_ari": top_loci((df["is_tp"] == 1) & (df["loh_status"] == "nonLOH")),
    "FP_LOH_top_blind_ari": top_loci((df["is_tp"] == 0) & (df["loh_status"] == "LOH")),
    "FP_nonLOH_top_blind_ari": top_loci((df["is_tp"] == 0) & (df["loh_status"] == "nonLOH")),
}

# canonical BRCA2 cell
brca2_row = df[(df["chrom"] == "chr13") & (df["somatic_pos"] == 32315128)]
brca2 = None
if len(brca2_row):
    r = brca2_row.iloc[0]
    brca2_clu = "hi" if r["blind_ari"] >= ari_med else "lo"
    brca2_cov = "hi" if r["n_paired_cpg"] >= cov_q75 else "lo"
    brca2_cell = f"{r['loh_status']}|{r['cn_class']}|clust_{brca2_clu}|cov_{brca2_cov}"
    # peers in same context cell
    peer_mask = ((df["loh_status"] == r["loh_status"]) & (df["cn_class"] == r["cn_class"])
                 & (df["clust_hi"] == (1 if brca2_clu == "hi" else 0))
                 & (df["cov_hi"] == (1 if brca2_cov == "hi" else 0)))
    peer = cell_stats(peer_mask, brca2_cell)
    brca2 = {
        "chrom": r["chrom"], "pos": int(r["somatic_pos"]), "axis": r["axis"],
        "is_tp": int(r["is_tp"]), "loh": r["loh_status"], "cn_class": r["cn_class"],
        "median_cn": float(r["median_cn"]), "n_paired_cpg": int(r["n_paired_cpg"]),
        "blind_ari": float(r["blind_ari"]), "abs_delta": float(r["abs_delta"]),
        "context_cell": brca2_cell,
        "cell_membership": {"clust": brca2_clu, "cov": brca2_cov},
        "cell_n": peer["n"], "cell_tp": peer["tp"], "cell_fp": peer["fp"],
        "cell_tp_rate": peer["tp_rate"],
    }

# blind_ari TP vs FP within LOH / nonLOH (Mann-Whitney, descriptive)
from scipy.stats import mannwhitneyu
def mw_ari(loh):
    sub = df[df["loh_status"] == loh]
    tp = sub[sub["is_tp"] == 1]["blind_ari"].values
    fp = sub[sub["is_tp"] == 0]["blind_ari"].values
    res = {"loh": loh, "tp_n": int(len(tp)), "fp_n": int(len(fp)),
           "tp_median_ari": float(np.median(tp)) if len(tp) else None,
           "fp_median_ari": float(np.median(fp)) if len(fp) else None}
    if len(tp) >= 1 and len(fp) >= 1:
        try:
            u, p = mannwhitneyu(tp, fp, alternative="two-sided")
            res["mannwhitney_p"] = float(p)
        except Exception:
            res["mannwhitney_p"] = None
    else:
        res["mannwhitney_p"] = None
    res["underpowered_fp_lt5"] = bool(len(fp) < 5)
    return res

ari_tp_vs_fp = {"LOH": mw_ari("LOH"), "nonLOH": mw_ari("nonLOH")}

# ---------------------------------------------------------------- assemble JSON
result = {
    "analysis": "C2_context_tpfp_clustering",
    "sample": "HCC1395",
    "input": IN_TSV,
    "n_loci": N,
    "n_unique_positions": int(n_unique_pos),
    "pseudoreplicated_within_subset": pseudoreplicated,
    "n_tp": n_tp,
    "n_fp": n_fp,
    "overall_tp_rate": overall_tp_rate,
    "thresholds": {
        "blind_ari_median": ari_med,
        "blind_ari_fixed": 0.3,
        "n_paired_cpg_q75": cov_q75,
    },
    "fp_distribution": {
        "fp_in_LOH": int(((df["loh_status"] == "LOH") & (df["is_tp"] == 0)).sum()),
        "fp_in_nonLOH": nonloh_fp_total,
        "note": "41 of 43 FP are in LOH; nonLOH has only 2 FP -> any nonLOH TP-enrichment claim is structurally underpowered.",
    },
    "H_ctx1": {
        "hypothesis": "LOH 內有甲基分群(high blind_ari)+異常拷貝 -> FP-enriched (TP-rate LOWER)",
        "primary_median_split": hctx1,
        "robustness_ari_0.3": hctx1_03,
        "verdict": v_ctx1,
    },
    "H_ctx2": {
        "hypothesis": "非LOH + coverage 高 + high-clustering -> TP-rate HIGHER",
        "primary_median_split": hctx2,
        "robustness_ari_0.3": hctx2_03,
        "verdict": v_ctx2,
        "structural_power_note": f"nonLOH total FP = {nonloh_fp_total}; with cov/clustering filters the in-cell FP count collapses -> cannot statistically demonstrate TP-enrichment.",
    },
    "marginal_grid_LOH_x_clustering": marginal_grid,
    "full_crosstab_LOHxCNxClustxCov": crosstab_cells,
    "tree_ranks_by_blind_ari": tree_ranks,
    "brca2_canonical": brca2,
    "blind_ari_tp_vs_fp_within_loh": ari_tp_vs_fp,
    "guard_note": "auc-confound-guard: TP-rate differences are CONTEXT effects (LOH/CN/coverage), NOT a methylation-feature classifier. No AUC/classifier is fit or claimed. FP n=43 (43 of 596) -> severe underpower in most context cells.",
}

with open(OUT_JSON, "w") as f:
    json.dump(result, f, indent=2)
print("wrote", OUT_JSON)

# ---------------------------------------------------------------- FIGURE
fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.38, wspace=0.28)

# Panel A: TP-rate grid LOH x clustering with n labels and CI
axA = fig.add_subplot(gs[0, 0])
order = ["LOH_clust_hi", "LOH_clust_lo", "nonLOH_clust_hi", "nonLOH_clust_lo"]
labels = ["LOH\nclust-hi", "LOH\nclust-lo", "nonLOH\nclust-hi", "nonLOH\nclust-lo"]
rates, errs_lo, errs_hi, ns, fps = [], [], [], [], []
for k in order:
    c = marginal_grid[k]
    rates.append(c["tp_rate"] * 100 if c["tp_rate"] is not None else 0)
    lo, hi = c["tp_rate_ci95"]
    errs_lo.append((c["tp_rate"] - lo) * 100 if lo is not None else 0)
    errs_hi.append((hi - c["tp_rate"]) * 100 if hi is not None else 0)
    ns.append(c["n"]); fps.append(c["fp"])
colors = ["#d62728" if f >= 5 else "#ff9896" for f in fps]
bars = axA.bar(range(4), rates, yerr=[errs_lo, errs_hi], capsize=4, color=colors, edgecolor="black")
axA.axhline(overall_tp_rate * 100, ls="--", color="gray", lw=1.2, label=f"overall {overall_tp_rate*100:.1f}%")
axA.set_xticks(range(4)); axA.set_xticklabels(labels, fontsize=9)
axA.set_ylabel("TP-rate (%)"); axA.set_ylim(0, 105)
axA.set_title("A. TP-rate by LOH x clustering (Wilson 95% CI)\nbar label: n / FP=fp", fontsize=10)
for i, b in enumerate(bars):
    axA.text(b.get_x() + b.get_width() / 2, 3, f"n={ns[i]}\nFP={fps[i]}", ha="center", va="bottom", fontsize=8)
axA.legend(fontsize=8, loc="lower left")

# Panel B: blind_ari distribution TP vs FP within LOH and nonLOH (strip/box)
axB = fig.add_subplot(gs[0, 1])
positions, data, plabels, pcolors = [], [], [], []
i = 0
for loh in ["LOH", "nonLOH"]:
    for tp_flag, col, nm in [(1, "#1f77b4", "TP"), (0, "#d62728", "FP")]:
        sub = df[(df["loh_status"] == loh) & (df["is_tp"] == tp_flag)]["blind_ari"].values
        positions.append(i); data.append(sub)
        plabels.append(f"{loh}\n{nm}\n(n={len(sub)})"); pcolors.append(col)
        i += 1
bp = axB.boxplot(data, positions=positions, widths=0.6, patch_artist=True, showfliers=False)
for patch, col in zip(bp["boxes"], pcolors):
    patch.set_facecolor(col); patch.set_alpha(0.4)
for j, (pos, sub, col) in enumerate(zip(positions, data, pcolors)):
    jitter = np.random.RandomState(0).normal(0, 0.06, len(sub))
    axB.scatter(np.full(len(sub), pos) + jitter, sub, s=10, color=col, alpha=0.5, zorder=3)
axB.set_xticks(positions); axB.set_xticklabels(plabels, fontsize=8)
axB.set_ylabel("blind_ari (read clustering)")
axB.set_title("B. blind_ari: TP vs FP within LOH / nonLOH\n(nonLOH FP n=2 -> uninformative)", fontsize=10)
# mark BRCA2
if brca2:
    axB.scatter([0], [brca2["blind_ari"]], marker="*", s=260, color="gold",
                edgecolor="black", zorder=5, label="BRCA2 (TP, 0.790)")
    axB.legend(fontsize=8, loc="upper right")

# Panel C: full crosstab heatmap-ish table (LOH x cn x clust x cov), color by TP-rate, annotate n & FP
axC = fig.add_subplot(gs[1, :])
axC.axis("off")
rows = sorted(crosstab_cells, key=lambda c: (c["loh"], c["cn_class"], -c["clust_hi"], -c["cov_hi"]))
header = ["context cell", "n", "TP", "FP", "TP-rate", "95% CI", "flag"]
table_data = []
cell_colors = []
for c in rows:
    rate = c["tp_rate"]
    lo, hi = c["tp_rate_ci95"]
    ci = f"[{lo*100:.0f},{hi*100:.0f}]" if lo is not None else "-"
    flag = "FP=0 uninform" if c["no_fp_uninformative"] else ("UNDERPWR" if c["underpowered_fp_lt5"] else "ok")
    table_data.append([c["label"], c["n"], c["tp"], c["fp"], f"{rate*100:.1f}%" if rate is not None else "-", ci, flag])
    # color row by FP presence
    if c["no_fp_uninformative"]:
        rowcol = "#f0f0f0"
    elif c["fp"] >= 5:
        rowcol = "#ffd6d6"
    else:
        rowcol = "#fff0e0"
    cell_colors.append([rowcol] * len(header))
tab = axC.table(cellText=table_data, colLabels=header, loc="center",
                cellColours=cell_colors, colColours=["#cccccc"] * len(header))
tab.auto_set_font_size(False); tab.set_fontsize(7.2); tab.scale(1, 1.25)
axC.set_title("C. Full see-tree-and-forest cross-tab: TP-rate by {LOH x CN x clustering x coverage}  "
              "(gray=FP0 uninformative, orange=FP<5 underpowered, red=FP>=5)", fontsize=10, pad=14)

fig.suptitle("C2: read-clustering x genomic context vs TP/FP (HCC1395, 596 BAM-pass HP loci; FP n=43)\n"
             "CONTEXT effect not a classifier (auc-confound-guard)", fontsize=12, fontweight="bold")
fig.savefig(OUT_PNG, dpi=140, bbox_inches="tight")
print("wrote", OUT_PNG)

# ---------------------------------------------------------------- console summary
print("\n=== SUMMARY ===")
print(f"N={N} TP={n_tp} FP={n_fp} overall_TP_rate={overall_tp_rate:.3f}")
print(f"blind_ari median={ari_med:.3f}  n_paired_cpg q75={cov_q75:.1f}")
print(f"FP in LOH={result['fp_distribution']['fp_in_LOH']}  FP in nonLOH={nonloh_fp_total}")
print(f"\nH-ctx1 [{v_ctx1}] cell n={hctx1['cell_n']} TPrate={hctx1['cell_tp_rate']:.3f} "
      f"(rest={hctx1['rest_tp_rate']:.3f}) OR={hctx1['odds_ratio']} "
      f"p1={hctx1['p_one_sided_directional']:.4f} cellFP={hctx1['cell_fp_n']}")
print(f"H-ctx2 [{v_ctx2}] cell n={hctx2['cell_n']} TPrate={hctx2['cell_tp_rate']} "
      f"(rest={hctx2['rest_tp_rate']:.3f}) OR={hctx2['odds_ratio']} "
      f"p1={hctx2['p_one_sided_directional']:.4f} cellFP={hctx2['cell_fp_n']}")
if brca2:
    print(f"\nBRCA2 cell={brca2['context_cell']} blind_ari={brca2['blind_ari']:.3f} "
          f"cell_n={brca2['cell_n']} cell_TP={brca2['cell_tp']} cell_FP={brca2['cell_fp']} "
          f"cell_TPrate={brca2['cell_tp_rate']:.3f}")
print("\nDONE")
