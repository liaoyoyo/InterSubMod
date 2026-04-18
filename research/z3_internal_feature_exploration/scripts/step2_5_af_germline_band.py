#!/usr/bin/env python3
"""Step 2.5 — AF in [0.4, 0.6] germline band FP separation.

Hypothesis H-Z3d: FPs with caller_af in [0.4, 0.6] are germline het leaks,
and CN (Coverage_Multiple) x HPFineNGroups x AlleleDelta can separate
subtypes:
  - CN~1 + AF[0.4,0.6] => missed LOH (cnLOH pattern)
  - CN>2 + AF[0.4,0.6] => CNV artifact (germline AF drift)
  - CN~2 + AF[0.4,0.6] + NGroups>=2 => true somatic (subclonal het)
  - CN~2 + AF[0.4,0.6] + NGroups=1  => germline het (filter candidate)

Scope: all 7 samples x TO mode x caller_af in [0.4, 0.6].
"""
import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')

MASTER = ("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
          "20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz")
OUT_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/z3_internal_feature_exploration"
FIG_DIR = os.path.join(OUT_DIR, "figures")
DATA_DIR = os.path.join(OUT_DIR, "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

SAMPLE_ORDER = ["H2009", "H1437", "HCC1395", "HCC1395_DORADO",
                "HCC1937", "HCC1954", "COLO829"]

COLS_NEEDED = [
    "sample", "mode", "truth_label",
    "core_loh_like", "to_loh_bed_hit",
    "caller_af", "HPFineNGroups", "NumReads",
    "Coverage_Multiple", "AlleleDelta",
    "caller_ad_ref", "caller_ad_alt",
]

print("Loading master dataset...")
df = pd.read_csv(MASTER, sep='\t', usecols=COLS_NEEDED, low_memory=False)

for c in COLS_NEEDED:
    if c not in ("sample", "mode", "truth_label", "core_loh_like", "to_loh_bed_hit"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
for col in ("core_loh_like", "to_loh_bed_hit"):
    df[col] = df[col].map({"True": True, "False": False, True: True, False: False}).fillna(False)

df["is_tp"] = df["truth_label"] == "TP"
df["loh_flag"] = np.where(df["mode"] == "paired", df["core_loh_like"], df["to_loh_bed_hit"])

# AD-based AF (more precise than caller_af for germline detection)
tot_ad = df["caller_ad_ref"].fillna(0) + df["caller_ad_alt"].fillna(0)
df["ad_af"] = np.where(tot_ad > 0, df["caller_ad_alt"] / tot_ad, np.nan)

# Focus on TO mode
dto = df[df["mode"] == "to"].copy()

# Germline band: caller_af in [0.4, 0.6]
dto["in_germ_band"] = (dto["caller_af"] >= 0.4) & (dto["caller_af"] <= 0.6)

print(f"\nTO mode total: {len(dto):,}")
print(f"TO in germline band [0.4, 0.6]: {dto['in_germ_band'].sum():,}")

# ── 1. Sample-level counts ─────────────────────────────────────────
print("\n=== Per-sample AF[0.4, 0.6] distribution ===")
sample_summary = []
for sample in SAMPLE_ORDER:
    ds = dto[(dto["sample"] == sample) & dto["in_germ_band"]]
    n = len(ds)
    n_tp = ds["is_tp"].sum()
    n_fp = n - n_tp
    # Check distribution over zones is via loh+NG, but simpler: report LOH presence
    n_loh = ds["loh_flag"].sum()
    tp_rate = n_tp / n if n else np.nan
    print(f"  {sample:16s}: N={n:6d}  TP={n_tp:5d}  FP={n_fp:5d}  TP_rate={tp_rate:.3f}  LOH={n_loh:5d}")
    sample_summary.append({
        "sample": sample, "n": n, "n_tp": int(n_tp), "n_fp": int(n_fp),
        "tp_rate": tp_rate, "n_loh": int(n_loh),
    })
pd.DataFrame(sample_summary).to_csv(
    os.path.join(DATA_DIR, "z3_af_germline_band_sample_summary.tsv"), sep='\t', index=False)


# ── 2. CN x NG joint stratification ────────────────────────────────
def cn_bin(cm):
    if pd.isna(cm):
        return "NA"
    if cm < 0.75:
        return "CN~1"
    if cm <= 1.5:
        return "CN~2"
    return "CN>2"


def ng_bin(ng):
    if pd.isna(ng):
        return "NA"
    if ng <= 1:
        return "NG<=1"
    if ng <= 3:
        return "NG2-3"
    return "NG>=4"


dto["cn_bin"] = dto["Coverage_Multiple"].apply(cn_bin)
dto["ng_bin"] = dto["HPFineNGroups"].apply(ng_bin)

cn_order = ["CN~1", "CN~2", "CN>2"]
ng_order = ["NG<=1", "NG2-3", "NG>=4"]

print("\n=== CN x NG stratified TP rate (AF in [0.4, 0.6]) ===")
strat_rows = []
for sample in SAMPLE_ORDER:
    ds = dto[(dto["sample"] == sample) & dto["in_germ_band"]]
    for cnb in cn_order:
        for ngb in ng_order:
            sub = ds[(ds["cn_bin"] == cnb) & (ds["ng_bin"] == ngb)]
            n = len(sub)
            n_tp = sub["is_tp"].sum()
            n_fp = n - n_tp
            tp_rate = n_tp / n if n else np.nan
            # LOH prevalence within this cell
            loh_pct = sub["loh_flag"].mean() if n else np.nan
            strat_rows.append({
                "sample": sample, "cn_bin": cnb, "ng_bin": ngb,
                "n": n, "n_tp": int(n_tp), "n_fp": int(n_fp),
                "tp_rate": tp_rate, "loh_pct": loh_pct,
            })

strat_df = pd.DataFrame(strat_rows)
strat_df.to_csv(os.path.join(DATA_DIR, "z3_af_germline_band_cn_ng.tsv"), sep='\t', index=False)


# Fisher test for critical cell: CN~2 + NG<=1 vs CN~2 + NG>=4 within each sample
print("\n=== H-Z3d critical test: CN~2 + NG<=1 vs CN~2 + NG>=4 ===")
fisher_rows = []
for sample in SAMPLE_ORDER:
    ds = dto[(dto["sample"] == sample) & dto["in_germ_band"] & (dto["cn_bin"] == "CN~2")]
    sub_low = ds[ds["ng_bin"] == "NG<=1"]
    sub_high = ds[ds["ng_bin"] == "NG>=4"]
    tp_low, fp_low = sub_low["is_tp"].sum(), len(sub_low) - sub_low["is_tp"].sum()
    tp_high, fp_high = sub_high["is_tp"].sum(), len(sub_high) - sub_high["is_tp"].sum()
    rate_low = tp_low / (tp_low + fp_low) if (tp_low + fp_low) else np.nan
    rate_high = tp_high / (tp_high + fp_high) if (tp_high + fp_high) else np.nan
    if (tp_low + fp_low) >= 5 and (tp_high + fp_high) >= 5:
        odds, p = stats.fisher_exact(
            [[tp_high, fp_high], [tp_low, fp_low]], alternative="greater")
    else:
        odds, p = np.nan, np.nan
    delta = rate_high - rate_low if not (np.isnan(rate_high) or np.isnan(rate_low)) else np.nan
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  {sample:16s}: NG<=1 TP={rate_low:.3f}(N={tp_low+fp_low:4d})  "
          f"NG>=4 TP={rate_high:.3f}(N={tp_high+fp_high:4d})  "
          f"delta={delta if not np.isnan(delta) else float('nan'):+.3f}  OR={odds:.2f}  {sig}")
    fisher_rows.append({
        "sample": sample,
        "low_tp": int(tp_low), "low_fp": int(fp_low), "low_rate": rate_low,
        "high_tp": int(tp_high), "high_fp": int(fp_high), "high_rate": rate_high,
        "delta": delta, "odds_ratio": odds, "p_value": p,
    })
fisher_df = pd.DataFrame(fisher_rows)
fisher_df.to_csv(os.path.join(DATA_DIR, "z3_af_germline_band_fisher.tsv"), sep='\t', index=False)


# ── 3. Heatmap figure — sample x (CN x NG) TP rate ─────────────────
pivot = strat_df.pivot_table(index="sample",
                             columns=["cn_bin", "ng_bin"], values="tp_rate")
# Order columns
cols = [(c, n) for c in cn_order for n in ng_order]
pivot = pivot.reindex(columns=cols)
pivot = pivot.reindex(SAMPLE_ORDER)

# Annotate with N counts
n_pivot = strat_df.pivot_table(index="sample",
                               columns=["cn_bin", "ng_bin"], values="n")
n_pivot = n_pivot.reindex(columns=cols).reindex(SAMPLE_ORDER)

fig, ax = plt.subplots(figsize=(14, 6))
im = ax.imshow(pivot.values, cmap='RdYlGn', aspect='auto', vmin=0.0, vmax=1.0)
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([f"{c}\n{n}" for c, n in cols], fontsize=9)
ax.set_yticks(range(len(SAMPLE_ORDER)))
ax.set_yticklabels(SAMPLE_ORDER, fontsize=10)

for i in range(pivot.values.shape[0]):
    for j in range(pivot.values.shape[1]):
        v = pivot.values[i, j]
        nn = n_pivot.values[i, j]
        if not np.isnan(v):
            color = 'white' if v < 0.3 or v > 0.7 else 'black'
            ax.text(j, i, f'{v:.2f}\nn={int(nn)}', ha='center', va='center',
                    fontsize=7, color=color,
                    fontweight='bold' if v < 0.3 and nn >= 50 else 'normal')

ax.set_title("H-Z3d: TP Rate in caller_af[0.4, 0.6] band by CN x NGroups\n"
             "Red + large N = strong germline-FP candidate (drop in confidence)",
             fontweight='bold')
plt.colorbar(im, ax=ax, label='TP Rate', shrink=0.8)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "z3_af_germline_band_heatmap.png"), dpi=150, bbox_inches='tight')
plt.close()

# ── 4. AlleleDelta AUC within AF[0.4, 0.6] ────────────────────────
print("\n=== AlleleDelta AUC within AF[0.4, 0.6] per sample ===")
from sklearn.metrics import roc_auc_score
ad_rows = []
for sample in SAMPLE_ORDER:
    ds = dto[(dto["sample"] == sample) & dto["in_germ_band"]]
    vals = ds["AlleleDelta"].values
    labs = ds["is_tp"].values
    m = ~np.isnan(vals)
    if m.sum() >= 20 and labs[m].sum() >= 5 and (~labs[m]).sum() >= 5:
        try:
            auc = roc_auc_score(labs[m], vals[m])
        except ValueError:
            auc = np.nan
    else:
        auc = np.nan
    auc_abs = max(auc, 1 - auc) if not np.isnan(auc) else np.nan
    print(f"  {sample:16s}: |AUC|={auc_abs:.3f}  N={m.sum()}")
    ad_rows.append({"sample": sample, "auc": auc, "auc_abs": auc_abs, "n": int(m.sum())})
pd.DataFrame(ad_rows).to_csv(
    os.path.join(DATA_DIR, "z3_af_germline_band_alleledelta_auc.tsv"), sep='\t', index=False)


# ── 5. AD-based AF dispersion (germline signature) ─────────────────
print("\n=== ad_af distribution in AF[0.4, 0.6] — TP vs FP ===")
for sample in SAMPLE_ORDER:
    ds = dto[(dto["sample"] == sample) & dto["in_germ_band"]]
    tp_af = ds[ds["is_tp"]]["ad_af"].dropna()
    fp_af = ds[~ds["is_tp"]]["ad_af"].dropna()
    if len(tp_af) >= 10 and len(fp_af) >= 10:
        tp_mean, tp_std = tp_af.mean(), tp_af.std()
        fp_mean, fp_std = fp_af.mean(), fp_af.std()
        # Distance of FP mean from 0.5 — if ~0, suggests germline enrichment
        print(f"  {sample:16s}: TP mean_af={tp_mean:.3f}+-{tp_std:.3f}  "
              f"FP mean_af={fp_mean:.3f}+-{fp_std:.3f}  "
              f"|FP-0.5|={abs(fp_mean-0.5):.3f}")


# ── 6. Decision summary ────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY — Step 2.5 H-Z3d Decision")
print("=" * 70)

# Check critical cell: CN~2 + NG<=1 TP rate per sample
critical = strat_df[(strat_df["cn_bin"] == "CN~2") & (strat_df["ng_bin"] == "NG<=1")]
n_samples_low = (critical["tp_rate"] < 0.3).sum()
n_samples_mid = ((critical["tp_rate"] >= 0.3) & (critical["tp_rate"] < 0.5)).sum()
n_samples_high = (critical["tp_rate"] >= 0.5).sum()
valid_samples = critical[critical["n"] >= 30]
n_valid = len(valid_samples)

print(f"\nCritical cell (CN~2 + NG<=1 + AF[0.4, 0.6]):")
print(f"  Samples with n>=30: {n_valid}/7")
for _, r in critical.iterrows():
    flag = ""
    if r["n"] >= 30:
        if r["tp_rate"] < 0.3:
            flag = " *** germline-enriched"
        elif r["tp_rate"] < 0.5:
            flag = " ** partial germline"
    print(f"    {r['sample']:16s}: n={r['n']:5d}  tp_rate={r['tp_rate']:.3f}  "
          f"loh_pct={r['loh_pct']:.2f}{flag}")

# Decision
germline_samples = valid_samples[valid_samples["tp_rate"] < 0.3]["sample"].tolist()
partial_samples = valid_samples[(valid_samples["tp_rate"] >= 0.3) & (valid_samples["tp_rate"] < 0.5)]["sample"].tolist()

print("\nDecision:")
if len(germline_samples) >= 5:
    print(f"  -> POSITIVE: {len(germline_samples)}/7 samples show CN~2 + NG<=1 + AF[0.4,0.6]")
    print(f"     TP_rate < 0.3, consistent with germline het leak.")
    print(f"     samples: {germline_samples}")
elif len(germline_samples) + len(partial_samples) >= 3:
    print(f"  -> CONDITIONAL: {len(germline_samples)} strong + {len(partial_samples)} partial")
    print(f"     strong: {germline_samples}")
    print(f"     partial: {partial_samples}")
else:
    print(f"  -> NEGATIVE: only {len(germline_samples)} samples with TP rate < 0.3")
    print("     AF[0.4,0.6] FPs cannot be cleanly separated by CN x NG.")

print(f"\nOutputs in {DATA_DIR}/ and {FIG_DIR}/")
