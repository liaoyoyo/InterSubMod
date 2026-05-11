#!/usr/bin/env python3
"""Phase 3 synthesis: Per-sample S1-S7 stratified filter framework across TO + paired.

Inputs
------
research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz
    Columns of interest:
      sample, mode, tp_label, LOH_Subtype, AF_class, AF, Coverage_Multiple,
      LOH_Bed_Overlap, HPFineNGroups, ...

Schemes (Thread B definitions, week-aligned)
-------------------------------------------
  S1: LOH_Strong AND AF_class == Extreme  (AF<0.1 or >0.9)
  S2: LOH_Subclone AND AF in (0.1,0.4] U [0.6,0.9)
  S3: Diploid Het = (LOH_Subtype=='None') AND (AF in [0.4,0.6]) AND (CovM in [0.65,1.33])
      (CN tier T1/T2 — near-diploid)
  S4: (LOH_Subtype=='None') AND AF_class == Extreme          (ambiguous)
  S5: (S1 OR S2 OR S3) AND NOT S4                            (combo-clean)
  S6: S1 AND HPFineNGroups >= 3
  S7: S5 AND HPFineNGroups >= 3

Outputs
-------
docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/
  ├── s1s7_per_sample.tsv
  ├── s1s7_heatmap_tp_rate.png
  ├── s1s7_heatmap_fold.png
  ├── s3_cross_sample_wilcoxon.json
  ├── ng_gap_aggregate.json          (Thread D NG=2 effective median)
  └── summary.json                   (top-level summary for report)
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

# ---------------------------------------------------------------------------
ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
SRC = ROOT / "research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
OUT_DIR = ROOT / "docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis"
OUT_DIR.mkdir(parents=True, exist_ok=True)

B1_JSON = ROOT / "research/tpfp_loh_af_kde_discrimination/data/obs18_wilcoxon_B1.json"
B3_JSON = ROOT / "research/tpfp_loh_af_kde_discrimination/data/B3_wilcoxon_gap_stats.json"

TO_SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
PAIRED_SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]

AF_EXTREME_LO = 0.1
AF_EXTREME_HI = 0.9
AF_HET_LO = 0.4
AF_HET_HI = 0.6
COVM_DIPLOID_LO = 0.65
COVM_DIPLOID_HI = 1.33

# ---------------------------------------------------------------------------

def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score 95% CI for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    margin = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denom
    return (max(0.0, center - margin), min(1.0, center + margin))


def classify_schemes(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate boolean columns S1..S7 per row."""
    af = df["AF"].astype(float)
    af_class = df["AF_class"].astype(str)
    loh = df["LOH_Subtype"].fillna("None").astype(str)
    covm = df["Coverage_Multiple"].astype(float)
    ng = df["HPFineNGroups"].astype(float)

    extreme = af_class.eq("Extreme")  # consistent with master column
    intermediate_af = (((af > AF_EXTREME_LO) & (af <= AF_HET_LO))
                       | ((af >= AF_HET_HI) & (af < AF_EXTREME_HI)))
    near_half = af.between(AF_HET_LO, AF_HET_HI, inclusive="both")
    diploid_cn = covm.between(COVM_DIPLOID_LO, COVM_DIPLOID_HI, inclusive="both")

    s1 = loh.eq("LOH_Strong") & extreme
    s2 = loh.eq("LOH_Subclone") & intermediate_af
    s3 = loh.eq("None") & near_half & diploid_cn
    s4 = loh.eq("None") & extreme
    s5 = (s1 | s2 | s3) & ~s4
    s6 = s1 & (ng >= 3)
    s7 = s5 & (ng >= 3)

    out = df.copy()
    out["S1"] = s1
    out["S2"] = s2
    out["S3"] = s3
    out["S4"] = s4
    out["S5"] = s5
    out["S6"] = s6
    out["S7"] = s7
    return out


def per_scheme_stats(sub: pd.DataFrame, baseline_tp: float) -> dict:
    n = len(sub)
    if n == 0:
        return {"n": 0, "k_tp": 0, "tp_rate": float("nan"),
                "ci_lo": float("nan"), "ci_hi": float("nan"),
                "fold_vs_baseline": float("nan"),
                "fp_reduction_vs_baseline": float("nan")}
    k = int(sub["tp_label"].astype(int).sum())
    tp = k / n
    lo, hi = wilson_ci(k, n)
    fold = tp / baseline_tp if baseline_tp > 0 else float("nan")
    # FP reduction (precision-framed): (FP rate baseline - FP rate scheme) / FP rate baseline
    fp_base = 1.0 - baseline_tp
    fp_scheme = 1.0 - tp
    fp_red = (fp_base - fp_scheme) / fp_base if fp_base > 0 else float("nan")
    return {"n": n, "k_tp": k, "tp_rate": tp,
            "ci_lo": lo, "ci_hi": hi,
            "fold_vs_baseline": fold,
            "fp_reduction_vs_baseline": fp_red}


def build_table(df: pd.DataFrame, samples: list[str], mode: str) -> pd.DataFrame:
    rows = []
    for sample in samples:
        sub_all = df[(df["sample"] == sample) & (df["mode"] == mode)]
        if len(sub_all) == 0:
            continue
        base_tp = sub_all["tp_label"].astype(int).mean()
        base_n = len(sub_all)
        rows.append({"sample": sample, "mode": mode, "scheme": "baseline",
                     "n": base_n, "k_tp": int(sub_all["tp_label"].astype(int).sum()),
                     "tp_rate": base_tp,
                     "ci_lo": wilson_ci(int(sub_all["tp_label"].astype(int).sum()), base_n)[0],
                     "ci_hi": wilson_ci(int(sub_all["tp_label"].astype(int).sum()), base_n)[1],
                     "fold_vs_baseline": 1.0,
                     "fp_reduction_vs_baseline": 0.0})
        for s in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
            stats_d = per_scheme_stats(sub_all[sub_all[s]], base_tp)
            stats_d.update({"sample": sample, "mode": mode, "scheme": s})
            rows.append(stats_d)
    return pd.DataFrame(rows)


def heatmap(df: pd.DataFrame, out_path: Path, value_col: str, title: str,
            cmap: str, vmin: float | None, vmax: float | None,
            annotate_n: bool = True, annotate_fold: bool = False) -> None:
    pivot = df.pivot_table(index=["mode", "sample"], columns="scheme",
                           values=value_col, aggfunc="first")
    # Order rows: TO first then paired_full; samples in canonical order
    order = ([("to_pileup", s) for s in TO_SAMPLES] +
             [("paired_full", s) for s in PAIRED_SAMPLES])
    order = [o for o in order if o in pivot.index]
    pivot = pivot.loc[order]
    scheme_order = ["baseline", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]
    scheme_order = [s for s in scheme_order if s in pivot.columns]
    pivot = pivot[scheme_order]

    # pivot table for n
    n_pivot = df.pivot_table(index=["mode", "sample"], columns="scheme",
                             values="n", aggfunc="first").loc[order, scheme_order]
    fold_pivot = df.pivot_table(index=["mode", "sample"], columns="scheme",
                                values="fold_vs_baseline", aggfunc="first"
                                ).loc[order, scheme_order]

    fig, ax = plt.subplots(figsize=(9, 0.5 * len(order) + 2.5))
    im = ax.imshow(pivot.values, cmap=cmap, aspect="auto",
                   vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(scheme_order)))
    ax.set_xticklabels(scheme_order, rotation=0)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{m.split('_')[0]}·{s}" for m, s in order])

    for i in range(len(order)):
        for j in range(len(scheme_order)):
            v = pivot.iloc[i, j]
            if pd.isna(v):
                continue
            n_v = n_pivot.iloc[i, j]
            label = f"{v:.2f}"
            if annotate_n and not pd.isna(n_v):
                label += f"\n(n={int(n_v)})"
            if annotate_fold:
                f_v = fold_pivot.iloc[i, j]
                if not pd.isna(f_v):
                    label += f"\n{f_v:.2f}x"
            color = "white" if v < (vmin + vmax) / 2 else "black" if cmap != "RdBu_r" else "black"
            ax.text(j, i, label, ha="center", va="center",
                    fontsize=6.5, color=color)
    ax.set_title(title)
    ax.set_xlabel("Scheme")
    fig.colorbar(im, ax=ax, label=value_col)
    plt.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def cross_sample_wilcoxon_s3(df_stats: pd.DataFrame, samples: list[str],
                             mode: str) -> dict:
    """Wilcoxon one-sided (greater) on per-sample (S3 tp_rate - baseline tp_rate)."""
    gaps = {}
    for sample in samples:
        row_base = df_stats[(df_stats["sample"] == sample)
                            & (df_stats["mode"] == mode)
                            & (df_stats["scheme"] == "baseline")]
        row_s3 = df_stats[(df_stats["sample"] == sample)
                          & (df_stats["mode"] == mode)
                          & (df_stats["scheme"] == "S3")]
        if row_base.empty or row_s3.empty:
            continue
        bt = float(row_base.iloc[0]["tp_rate"])
        st = row_s3.iloc[0]["tp_rate"]
        if pd.isna(st):
            continue
        gaps[sample] = float(st) - bt
    vals = list(gaps.values())
    if len(vals) < 3:
        return {"gaps": gaps, "n": len(vals),
                "error": "insufficient samples for Wilcoxon"}
    try:
        res = stats.wilcoxon(vals, alternative="greater", zero_method="wilcox",
                             method="exact")
        return {"gaps": gaps,
                "n": len(vals),
                "median_gap": float(np.median(vals)),
                "mean_gap": float(np.mean(vals)),
                "min_gap": float(np.min(vals)),
                "max_gap": float(np.max(vals)),
                "wilcoxon_stat": float(res.statistic),
                "wilcoxon_p_greater": float(res.pvalue),
                "significant_0.05": bool(res.pvalue < 0.05)}
    except Exception as e:
        return {"gaps": gaps, "n": len(vals), "error": str(e)}


def ng_gap_aggregate(b1_json: dict, b2_overlay_note: str) -> dict:
    """Produce a Thread D aggregate: NG=2 gap median with/without HCC1954.

    The B2 analysis demonstrated HCC1954 Outer cross-het TP ~0.08 is caller
    FP background, not phasing failure. Effective median excludes HCC1954.
    """
    rows = b1_json["per_sample"]
    gaps_all = {r["sample"]: r["gap"] for r in rows}
    gaps_excl = {k: v for k, v in gaps_all.items() if k != "HCC1954"}
    vals_all = list(gaps_all.values())
    vals_excl = list(gaps_excl.values())

    return {
        "gaps_all_samples": gaps_all,
        "median_gap_all": float(np.median(vals_all)) if vals_all else None,
        "gaps_excl_hcc1954_caller_noise": gaps_excl,
        "median_gap_excl_hcc1954": float(np.median(vals_excl)) if vals_excl else None,
        "note": b2_overlay_note,
    }


# ---------------------------------------------------------------------------

def main() -> None:
    print(f"[load] {SRC}")
    df = pd.read_csv(SRC, sep="\t", low_memory=False)
    print(f"[load] rows={len(df)}, samples={sorted(df['sample'].unique())}, modes={sorted(df['mode'].unique())}")

    df = classify_schemes(df)

    # Per-sample tables
    to_table = build_table(df, TO_SAMPLES, "to_pileup")
    paired_table = build_table(df, PAIRED_SAMPLES, "paired_full")
    full_table = pd.concat([to_table, paired_table], ignore_index=True)
    tsv_path = OUT_DIR / "s1s7_per_sample.tsv"
    full_table.to_csv(tsv_path, sep="\t", index=False)
    print(f"[write] {tsv_path} rows={len(full_table)}")

    # Heatmaps
    heatmap(full_table, OUT_DIR / "s1s7_heatmap_tp_rate.png",
            value_col="tp_rate",
            title="S1-S7 TP rate (with n)  · 6 TO + 7 paired",
            cmap="viridis", vmin=0.0, vmax=1.0, annotate_n=True)
    heatmap(full_table, OUT_DIR / "s1s7_heatmap_fold.png",
            value_col="fold_vs_baseline",
            title="S1-S7 fold vs baseline  · 6 TO + 7 paired",
            cmap="RdBu_r", vmin=0.5, vmax=2.0,
            annotate_n=True, annotate_fold=False)

    # Cross-sample Wilcoxon on S3 (TO)
    s3_wilcoxon_to = cross_sample_wilcoxon_s3(full_table, TO_SAMPLES, "to_pileup")
    s3_wilcoxon_paired = cross_sample_wilcoxon_s3(full_table, PAIRED_SAMPLES,
                                                  "paired_full")
    wilcoxon_out = {"to_pileup": s3_wilcoxon_to, "paired_full": s3_wilcoxon_paired}
    with open(OUT_DIR / "s3_cross_sample_wilcoxon.json", "w") as fh:
        json.dump(wilcoxon_out, fh, indent=2)
    print("[write] s3_cross_sample_wilcoxon.json")

    # Thread D aggregate
    with open(B1_JSON) as fh:
        b1 = json.load(fh)
    ng_agg = ng_gap_aggregate(b1, b2_overlay_note=(
        "B2 (2026-04-23) 證實 HCC1954 TO Outer cross-het TP=0.08 為 caller-FP "
        "背景（84% all-FP regions），非 phasing failure；effective median "
        "應排除 HCC1954。"))
    with open(OUT_DIR / "ng_gap_aggregate.json", "w") as fh:
        json.dump(ng_agg, fh, indent=2)
    print("[write] ng_gap_aggregate.json")

    # Top-level summary
    summary = {
        "analysis": "Phase3_Cross_Sample_S1S7_Synthesis",
        "date": "2026-04-23",
        "source_tsv": str(SRC),
        "to_samples": TO_SAMPLES,
        "paired_samples": PAIRED_SAMPLES,
        "scheme_definitions": {
            "S1": "LOH_Strong AND AF_class=Extreme",
            "S2": "LOH_Subclone AND AF in (0.1, 0.4] U [0.6, 0.9)",
            "S3": "LOH=None AND AF in [0.4, 0.6] AND CovM in [0.65, 1.33]",
            "S4": "LOH=None AND AF_class=Extreme (ambiguous)",
            "S5": "(S1 OR S2 OR S3) AND NOT S4",
            "S6": "S1 AND HPFineNGroups>=3",
            "S7": "S5 AND HPFineNGroups>=3",
        },
        "s3_wilcoxon": wilcoxon_out,
        "thread_d_ng_gap": ng_agg,
    }
    # Add compact per-scheme summary across TO
    to_summary = {}
    for s in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
        rows = full_table[(full_table["mode"] == "to_pileup")
                          & (full_table["scheme"] == s)]
        tps = rows["tp_rate"].dropna().tolist()
        ns = rows["n"].dropna().tolist()
        to_summary[s] = {
            "tp_rate_min": float(np.min(tps)) if tps else None,
            "tp_rate_median": float(np.median(tps)) if tps else None,
            "tp_rate_max": float(np.max(tps)) if tps else None,
            "total_n": int(np.sum(ns)) if ns else 0,
            "sample_count": len(tps),
        }
    summary["to_scheme_cross_sample"] = to_summary
    # Paired
    paired_summary = {}
    for s in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
        rows = full_table[(full_table["mode"] == "paired_full")
                          & (full_table["scheme"] == s)]
        tps = rows["tp_rate"].dropna().tolist()
        ns = rows["n"].dropna().tolist()
        paired_summary[s] = {
            "tp_rate_min": float(np.min(tps)) if tps else None,
            "tp_rate_median": float(np.median(tps)) if tps else None,
            "tp_rate_max": float(np.max(tps)) if tps else None,
            "total_n": int(np.sum(ns)) if ns else 0,
            "sample_count": len(tps),
        }
    summary["paired_scheme_cross_sample"] = paired_summary

    # Decide verdict on S3 TO
    to_s3 = summary["to_scheme_cross_sample"]["S3"]
    if to_s3["sample_count"] >= 5 and to_s3["tp_rate_median"] is not None:
        if to_s3["tp_rate_median"] >= 0.85 and s3_wilcoxon_to.get("wilcoxon_p_greater", 1.0) < 0.05:
            s3_verdict = "POSITIVE"
        elif to_s3["tp_rate_median"] >= 0.75:
            s3_verdict = "CONDITIONAL_POSITIVE"
        else:
            s3_verdict = "CONDITIONAL_NEGATIVE"
    else:
        s3_verdict = "INCONCLUSIVE"
    summary["s3_verdict_TO"] = s3_verdict

    # S5 verdict parallel
    to_s5 = summary["to_scheme_cross_sample"]["S5"]
    if to_s5["sample_count"] >= 5 and to_s5["tp_rate_median"] is not None:
        if to_s5["tp_rate_median"] >= 0.85:
            s5_verdict = "POSITIVE"
        elif to_s5["tp_rate_median"] >= 0.75:
            s5_verdict = "CONDITIONAL_POSITIVE"
        else:
            s5_verdict = "CONDITIONAL_NEGATIVE"
    else:
        s5_verdict = "INCONCLUSIVE"
    summary["s5_verdict_TO"] = s5_verdict

    # NG marginal: S6 vs S1, S7 vs S5 per-sample delta
    ng_marginal = {"S6_minus_S1": [], "S7_minus_S5": []}
    for sample in TO_SAMPLES:
        base_map = {}
        for scheme in ["S1", "S5", "S6", "S7"]:
            row = full_table[(full_table["sample"] == sample)
                             & (full_table["mode"] == "to_pileup")
                             & (full_table["scheme"] == scheme)]
            if not row.empty:
                base_map[scheme] = {"tp_rate": float(row.iloc[0]["tp_rate"])
                                    if not pd.isna(row.iloc[0]["tp_rate"]) else None,
                                    "n": int(row.iloc[0]["n"])}
        if base_map.get("S1") and base_map.get("S6") and base_map["S1"]["tp_rate"] is not None:
            if base_map["S6"]["tp_rate"] is not None:
                ng_marginal["S6_minus_S1"].append({
                    "sample": sample,
                    "delta": base_map["S6"]["tp_rate"] - base_map["S1"]["tp_rate"],
                    "n_S1": base_map["S1"]["n"],
                    "n_S6": base_map["S6"]["n"]})
        if base_map.get("S5") and base_map.get("S7") and base_map["S5"]["tp_rate"] is not None:
            if base_map["S7"]["tp_rate"] is not None:
                ng_marginal["S7_minus_S5"].append({
                    "sample": sample,
                    "delta": base_map["S7"]["tp_rate"] - base_map["S5"]["tp_rate"],
                    "n_S5": base_map["S5"]["n"],
                    "n_S7": base_map["S7"]["n"]})
    for key in ["S6_minus_S1", "S7_minus_S5"]:
        deltas = [r["delta"] for r in ng_marginal[key]]
        if deltas:
            ng_marginal[f"{key}_median_delta"] = float(np.median(deltas))
            ng_marginal[f"{key}_n_samples"] = len(deltas)
    summary["ng_marginal_gain_TO"] = ng_marginal

    with open(OUT_DIR / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[write] summary.json")

    # Print compact console summary
    print("\n=== Phase 3 cross-sample summary (TO, 6 samples) ===")
    for s in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
        st = summary["to_scheme_cross_sample"][s]
        tpr = st["tp_rate_median"]
        print(f"  {s}: median TP={tpr:.3f}" if tpr is not None else f"  {s}: n/a",
              f"min={st['tp_rate_min']:.3f} max={st['tp_rate_max']:.3f}"
              if st['tp_rate_min'] is not None else "",
              f"n={st['total_n']} samples={st['sample_count']}")
    print(f"\n  S3 Wilcoxon (TO, greater): p={s3_wilcoxon_to.get('wilcoxon_p_greater', 'n/a')}")
    print(f"  S3 Verdict: {s3_verdict}")
    print(f"  S5 Verdict: {s5_verdict}")
    print(f"  Thread D NG=2 gap median (excl HCC1954): "
          f"{ng_agg['median_gap_excl_hcc1954']:.3f}")
    print(f"  NG marginal S6-S1 median delta: "
          f"{ng_marginal.get('S6_minus_S1_median_delta', 'n/a')}")
    print(f"  NG marginal S7-S5 median delta: "
          f"{ng_marginal.get('S7_minus_S5_median_delta', 'n/a')}")


if __name__ == "__main__":
    main()
