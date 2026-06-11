#!/usr/bin/env python3
"""
32_credibility_patch.py — data-ize the red-team challenges (pure-Python, existing TSVs).

Patches (all canonical口徑: |mean_delta|>=0.1, Bonferroni alpha = 0.05/51091):
  P1  FDR vs Bonferroni full picture (+ |Δβ| sweep) + is the NEGATIVE correction-robust?
  P2  chr8 leave-one-chromosome-out jackknife of the FP/TP strong-ASM OR
  P3  coverage-stratified detection-power curve (where could ASM hide?)
  P4  rigorous germline-het null: pool + coverage-match + bootstrap CI of (TP - null)

Read-only inputs: genome_survey_v2/asm_dualaxis_{tp,fp,germline_negctrl,germline_negctrl2}.tsv
Outputs: genome_survey_v2/credibility_patch_results.json  +  figures/credibility_patch/*.png
"""
import csv, json, os, math, random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

random.seed(20260531); np.random.seed(20260531)
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
GS = os.path.join(ROOT, "genome_survey_v2")
FIGDIR = os.path.join(ROOT, "figures", "credibility_patch")
os.makedirs(FIGDIR, exist_ok=True)
try:
    import sys; sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod")
    from scripts.lib.plot_setup import setup_plot_style
    setup_plot_style(base_size=11, dpi=150)
except Exception:
    pass

def load(name):
    return list(csv.DictReader(open(os.path.join(GS, name)), delimiter="\t"))
def fl(x):
    try: return float(x)
    except: return None

tp = load("asm_dualaxis_tp.tsv")
fp = load("asm_dualaxis_fp.tsv")
g1 = load("asm_dualaxis_germline_negctrl.tsv")
g2 = load("asm_dualaxis_germline_negctrl2.tsv")
gnull = g1 + g2  # pooled germline-het null (1334 records)

N_BONF = 51091
BONF = 0.05 / N_BONF
DTHR = 0.10

def valid(rows):
    return [r for r in rows if fl(r["wilcoxon_p"]) is not None and fl(r["mean_delta"]) is not None]
def is_strong(r, alpha=BONF, dthr=DTHR):
    return fl(r["wilcoxon_p"]) < alpha and abs(fl(r["mean_delta"])) >= dthr
def strong(rows, alpha=BONF, dthr=DTHR):
    return [r for r in valid(rows) if is_strong(r, alpha, dthr)]
def bh_threshold(rows, q=0.05):
    ps = sorted(fl(r["wilcoxon_p"]) for r in valid(rows))
    m = len(ps); thr = 0.0
    for i, p in enumerate(ps, 1):
        if p <= (i / m) * q:
            thr = p
    return thr

results = {"meta": {"date": "2026-05-31", "alpha_bonf": BONF, "dthr": DTHR,
                    "n_tp": len(tp), "n_fp": len(fp), "n_gnull": len(gnull)}}

# =========================================================================
# P1 — FDR vs Bonferroni full picture
# =========================================================================
tpv = valid(tp)
bh = bh_threshold(tp)
bonf_strong = strong(tp)
fdr_strong = [r for r in tpv if fl(r["wilcoxon_p"]) <= bh and abs(fl(r["mean_delta"])) >= DTHR]
# |Δβ| sweep on Bonferroni-passers and FDR-passers
bonf_pass = [r for r in tpv if fl(r["wilcoxon_p"]) < BONF]
fdr_pass = [r for r in tpv if fl(r["wilcoxon_p"]) <= bh]
sweep = {}
for d in [0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.30]:
    sweep[d] = {"bonf": sum(1 for r in bonf_pass if abs(fl(r["mean_delta"])) >= d),
                "fdr": sum(1 for r in fdr_pass if abs(fl(r["mean_delta"])) >= d)}
# Is the anti-discriminative + null>=TP correction-robust? recompute under FDR
def strong_fdr(rows):
    bh_r = bh  # use TP's BH threshold as a common p-cut for comparability
    return [r for r in valid(rows) if fl(r["wilcoxon_p"]) <= bh_r and abs(fl(r["mean_delta"])) >= DTHR]
tp_fdr_rate = len(fdr_strong) / len(tp)
fp_fdr_rate = len(strong_fdr(fp)) / len(fp)
gn_fdr_rate = len(strong_fdr(gnull)) / len(gnull)
results["P1_fdr_vs_bonf"] = {
    "bonf_strong": len(bonf_strong), "bonf_pct": round(100*len(bonf_strong)/len(tp), 3),
    "fdr_strong": len(fdr_strong), "fdr_pct": round(100*len(fdr_strong)/len(tp), 2),
    "fold": round(len(fdr_strong)/max(1, len(bonf_strong)), 1),
    "bh_threshold": bh, "sweep": sweep,
    "correction_robust_check": {
        "TP_strong_rate_bonf": round(100*len(bonf_strong)/len(tp), 3),
        "FP_strong_rate_bonf": round(100*len(strong(fp))/len(fp), 3),
        "gnull_strong_rate_bonf": round(100*len(strong(gnull))/len(gnull), 3),
        "TP_strong_rate_fdr": round(100*tp_fdr_rate, 2),
        "FP_strong_rate_fdr": round(100*fp_fdr_rate, 2),
        "gnull_strong_rate_fdr": round(100*gn_fdr_rate, 2),
        "note": "If null>=TP holds under BOTH Bonferroni AND FDR -> NEGATIVE is correction-robust."
    }
}

# =========================================================================
# P2 — chr8 leave-one-chromosome-out jackknife of FP/TP strong-ASM OR
# =========================================================================
def or_fp_tp(tprows, fprows):
    ts, fs = strong(tprows), strong(fprows)
    a, b = len(fs), len(fprows) - len(fs)   # FP strong / not
    c, d = len(ts), len(tprows) - len(ts)   # TP strong / not
    OR = (a * d) / (b * c) if b * c else float("inf")  # FP-enrichment OR
    return OR, len(fs), len(fprows), len(ts), len(tprows)
chroms = sorted({r["chrom"] for r in tp}, key=lambda c: (len(c), c))
OR_full = or_fp_tp(tp, fp)[0]
loo = {}
for ch in chroms:
    tpx = [r for r in tp if r["chrom"] != ch]
    fpx = [r for r in fp if r["chrom"] != ch]
    loo[ch] = round(or_fp_tp(tpx, fpx)[0], 3)
loo_vals = list(loo.values())
fp_strong_by_chr = {}
for r in strong(fp):
    fp_strong_by_chr[r["chrom"]] = fp_strong_by_chr.get(r["chrom"], 0) + 1
results["P2_chr8_jackknife"] = {
    "OR_full": round(OR_full, 3),
    "OR_leave_one_chrom_out": loo,
    "OR_drop_chr8": loo.get("chr8"),
    "jackknife_min": round(min(loo_vals), 3), "jackknife_max": round(max(loo_vals), 3),
    "fp_strong_by_chrom": dict(sorted(fp_strong_by_chr.items(), key=lambda x: -x[1])),
    "fp_strong_total": len(strong(fp)),
    "chr8_share_of_fp_strong": round(fp_strong_by_chr.get("chr8", 0) / max(1, len(strong(fp))), 3),
}

# =========================================================================
# P3 — coverage-stratified detection-power curve
# =========================================================================
bins = [(5, 10), (10, 15), (15, 25), (25, 50), (50, 100), (100, 200), (200, 10**9)]
def wilcoxon_min_p(n):
    # two-sided signed-rank smallest achievable p ~ 2 / 2^n (all same-sign), capped
    return min(1.0, 2.0 / (2 ** n)) if n <= 20 else 0.0
power = []
for lo, hi in bins:
    sub = [r for r in tpv if lo <= int(r["n_paired_cpg"]) < hi]
    ss = [r for r in sub if is_strong(r)]
    sig = [r for r in sub if fl(r["wilcoxon_p"]) < BONF]
    min_strong_d = min((abs(fl(r["mean_delta"])) for r in ss), default=None)
    nmid = (lo + min(hi, lo + 5))
    power.append({
        "bin": f"{lo}-{hi if hi < 10**8 else '+'}", "n": len(sub),
        "n_bonf_sig": len(sig), "n_strong": len(ss),
        "strong_pct": round(100 * len(ss) / max(1, len(sub)), 3),
        "min_strong_abs_dbeta": round(min_strong_d, 3) if min_strong_d else None,
        "wilcoxon_floor_p_at_lo": round(wilcoxon_min_p(lo), 4),
        "can_reach_bonf": wilcoxon_min_p(lo) < BONF,
    })
n_universe = 39447; n_tested_positions = 30511
results["P3_coverage_power"] = {
    "n_universe_TP": n_universe, "n_tested_positions": n_tested_positions,
    "n_excluded": n_universe - n_tested_positions,
    "excluded_pct": round(100 * (n_universe - n_tested_positions) / n_universe, 1),
    "exclusion_reason": "no testable axis (no nearby CpG OR <5 paired CpG in any read-group)",
    "by_bin": power,
    "interpretation": "Bins where wilcoxon_floor_p_at_lo > Bonferroni alpha CANNOT reach strong-ASM regardless of true effect -> ASM at sparse loci is undetectable, not absent.",
}

# =========================================================================
# P4 — rigorous germline-het null: pool + coverage-match + bootstrap CI
# =========================================================================
def by_axis(rows, ax):
    return [r for r in rows if r["axis_type"] == ax]
def strong_rate(rows):
    v = valid(rows)
    return len(strong(rows)) / len(v) if v else 0.0
def median_abs_delta(rows):
    xs = [abs(fl(r["mean_delta"])) for r in valid(rows)]
    return float(np.median(xs)) if xs else None

# coverage match by n_paired_cpg bins
cov_bins = [(5, 15), (15, 30), (30, 60), (60, 120), (120, 10**9)]
def cov_bin_of(r):
    n = int(r["n_paired_cpg"])
    for i, (lo, hi) in enumerate(cov_bins):
        if lo <= n < hi: return i
    return len(cov_bins) - 1
matched = []
for i, (lo, hi) in enumerate(cov_bins):
    tpb = [r for r in valid(tp) if cov_bin_of(r) == i]
    gnb = [r for r in valid(gnull) if cov_bin_of(r) == i]
    matched.append({
        "cov_bin": f"{lo}-{hi if hi < 10**8 else '+'}",
        "n_tp": len(tpb), "n_null": len(gnb),
        "tp_strong_pct": round(100 * strong_rate(tpb), 3),
        "null_strong_pct": round(100 * strong_rate(gnb), 3),
        "tp_median_abs_dbeta": round(median_abs_delta(tpb), 4) if tpb else None,
        "null_median_abs_dbeta": round(median_abs_delta(gnb), 4) if gnb else None,
    })

# bootstrap CI of (TP_strong_rate - null_strong_rate), overall + coverage-stratified pooled
def boot_diff(tprows, gnrows, nboot=3000):
    tpv2, gnv2 = valid(tprows), valid(gnrows)
    tp_flags = np.array([1 if is_strong(r) else 0 for r in tpv2])
    gn_flags = np.array([1 if is_strong(r) else 0 for r in gnv2])
    diffs = []
    for _ in range(nboot):
        bt = np.random.choice(tp_flags, len(tp_flags), replace=True).mean()
        bg = np.random.choice(gn_flags, len(gn_flags), replace=True).mean()
        diffs.append(bt - bg)
    diffs = np.array(diffs)
    return float(diffs.mean()), float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))
d_all = boot_diff(tp, gnull)
d_hp = boot_diff(by_axis(tp, "HP"), by_axis(gnull, "HP"))
d_al = boot_diff(by_axis(tp, "ALLELE"), by_axis(gnull, "ALLELE"))
results["P4_rigorous_null"] = {
    "pooled_null_n": len(gnull),
    "negctrl_n": len(g1), "negctrl_strong_pct": round(100*strong_rate(g1), 3),
    "negctrl2_n": len(g2), "negctrl2_strong_pct": round(100*strong_rate(g2), 3),
    "pooled_null_strong_pct": round(100*strong_rate(gnull), 3),
    "tp_strong_pct": round(100*strong_rate(tp), 3),
    "coverage_matched": matched,
    "bootstrap_diff_TP_minus_null": {
        "all_axis": {"point": round(100*d_all[0], 3), "ci95": [round(100*d_all[1], 3), round(100*d_all[2], 3)]},
        "HP_axis": {"point": round(100*d_hp[0], 3), "ci95": [round(100*d_hp[1], 3), round(100*d_hp[2], 3)]},
        "ALLELE_axis": {"point": round(100*d_al[0], 3), "ci95": [round(100*d_al[1], 3), round(100*d_al[2], 3)]},
    },
    "verdict": "If CI of (TP-null) includes 0 or is negative -> somatic strong-ASM NOT above germline baseline (NEGATIVE upgraded to L1).",
}

json.dump(results, open(os.path.join(GS, "credibility_patch_results.json"), "w"), indent=1)
print("WROTE credibility_patch_results.json")
print(json.dumps({k: (v if k == "meta" else "...") for k, v in results.items()}, indent=1)[:200])

# =========================================================================
# FIGURES
# =========================================================================
AC = "#D97757"; GR = "#15803D"; RD = "#C2410C"; BL = "#1E3A8A"; GY = "#6B6B66"

# Fig 1 — FDR vs Bonferroni sweep
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ds = sorted(sweep.keys())
ax.plot(ds, [sweep[d]["fdr"] for d in ds], "-o", color=AC, label="FDR (BH .05) passers")
ax.plot(ds, [sweep[d]["bonf"] for d in ds], "-s", color=BL, label="Bonferroni passers")
ax.axvline(0.10, ls="--", color=GY, lw=1)
ax.annotate(f"172 (Bonf)\n3,018 (FDR)\nat |Δβ|≥0.10", xy=(0.10, 172), xytext=(0.165, 1500),
            fontsize=9, arrowprops=dict(arrowstyle="->", color=GY))
ax.set_xlabel("|Δβ| effect-size gate"); ax.set_ylabel("strong-ASM count (log)")
ax.set_yscale("log"); ax.set_title("P1 — 'rare 0.34%' is correction-dependent: Bonferroni 172 vs FDR 3,018 (17.5x)")
ax.legend(); ax.grid(True, alpha=0.25)
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "p1_fdr_vs_bonf.png"), dpi=150, bbox_inches="tight"); plt.close()

# Fig 2 — chr8 jackknife
fig, ax = plt.subplots(figsize=(7.6, 5.2))
items = sorted(loo.items(), key=lambda x: x[1])
labels = [k for k, v in items]; vals = [v for k, v in items]
colors = [RD if k == "chr8" else AC for k in labels]
ax.barh(range(len(labels)), vals, color=colors, alpha=0.85)
ax.axvline(OR_full, ls="--", color=GY, lw=1.2, label=f"full OR={OR_full:.2f}")
ax.set_yticks(range(len(labels))); ax.set_yticklabels([f"drop {k}" for k in labels], fontsize=8.5)
ax.set_xlabel("FP-enrichment OR (FP strong-ASM / TP strong-ASM)")
ax.set_title(f"P2 — leave-one-chromosome-out: dropping chr8 alone moves OR {OR_full:.2f}→{loo.get('chr8'):.2f}")
ax.legend(); ax.grid(True, axis="x", alpha=0.25)
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "p2_chr8_jackknife.png"), dpi=150, bbox_inches="tight"); plt.close()

# Fig 3 — coverage power curve
fig, ax1 = plt.subplots(figsize=(7.6, 4.8))
xb = [p["bin"] for p in power]; xi = range(len(xb))
strong_pct = [p["strong_pct"] for p in power]
ax1.bar(list(xi), strong_pct, color=AC, alpha=0.8, label="strong-ASM % of bin")
ax1.set_ylabel("strong-ASM % within coverage bin", color=AC)
ax1.set_xticks(list(xi)); ax1.set_xticklabels(xb, rotation=30, ha="right", fontsize=8.5)
ax1.set_xlabel("n_paired_cpg (coverage) bin")
for i, p in enumerate(power):
    if not p["can_reach_bonf"]:
        ax1.text(i, 0.02, "✗ can't\nreach\nBonf", ha="center", va="bottom", fontsize=7, color=RD)
ax1.set_title("P3 — detection-power: sparse-coverage loci CANNOT reach strong-ASM (ASM undetectable, not absent)")
ax1.grid(True, axis="y", alpha=0.25)
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "p3_coverage_power.png"), dpi=150, bbox_inches="tight"); plt.close()

# Fig 4 — null vs TP coverage-matched
fig, ax = plt.subplots(figsize=(7.6, 4.8))
xb = [m["cov_bin"] for m in matched]; xi = np.arange(len(xb))
w = 0.38
ax.bar(xi - w/2, [m["tp_strong_pct"] for m in matched], w, color=AC, label="somatic TP", alpha=0.85)
ax.bar(xi + w/2, [m["null_strong_pct"] for m in matched], w, color=GR, label="germline-het null", alpha=0.85)
ax.set_xticks(xi); ax.set_xticklabels(xb, fontsize=9)
ax.set_xlabel("n_paired_cpg (coverage) bin"); ax.set_ylabel("strong-ASM % (coverage-matched)")
d = results["P4_rigorous_null"]["bootstrap_diff_TP_minus_null"]["all_axis"]
ax.set_title(f"P4 — coverage-matched: somatic ≈/≤ germline baseline\nbootstrap (TP−null) all-axis = {d['point']:+.2f}% CI[{d['ci95'][0]:+.2f},{d['ci95'][1]:+.2f}]")
ax.legend(); ax.grid(True, axis="y", alpha=0.25)
plt.tight_layout(); plt.savefig(os.path.join(FIGDIR, "p4_null_matched.png"), dpi=150, bbox_inches="tight"); plt.close()

print("WROTE 4 figures to", FIGDIR)
# console summary
print("\n=== SUMMARY ===")
print(f"P1 FDR/Bonf: 172 -> {results['P1_fdr_vs_bonf']['fdr_strong']} ({results['P1_fdr_vs_bonf']['fold']}x); "
      f"correction-robust null check: TP_fdr={results['P1_fdr_vs_bonf']['correction_robust_check']['TP_strong_rate_fdr']}% "
      f"null_fdr={results['P1_fdr_vs_bonf']['correction_robust_check']['gnull_strong_rate_fdr']}%")
print(f"P2 chr8: OR full={OR_full:.2f} drop-chr8={loo.get('chr8'):.2f} jackknife range [{min(loo_vals):.2f},{max(loo_vals):.2f}]; chr8 share of FP-strong={results['P2_chr8_jackknife']['chr8_share_of_fp_strong']}")
print(f"P3 coverage: excluded {results['P3_coverage_power']['n_excluded']} ({results['P3_coverage_power']['excluded_pct']}%); bins that can't reach Bonf: "
      + ",".join(p["bin"] for p in power if not p["can_reach_bonf"]))
print(f"P4 null: pooled null strong={results['P4_rigorous_null']['pooled_null_strong_pct']}% vs TP {results['P4_rigorous_null']['tp_strong_pct']}%; "
      f"bootstrap(TP-null) all={d['point']:+.2f}% CI[{d['ci95'][0]:+.2f},{d['ci95'][1]:+.2f}] HP={results['P4_rigorous_null']['bootstrap_diff_TP_minus_null']['HP_axis']['point']:+.2f}%")
