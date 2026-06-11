#!/usr/bin/env python3
"""Step 1 stages 5-7: deltas, trajectory, off/on compare, H1a/H1b/H1c judgments.

Inputs (from build_three_way_master.py):
- step1_master_three_way.tsv (wide-format, region × version × flag features)

Outputs:
- step1_delta_v5_vs_v3f.tsv / step1_delta_v6_vs_v5.tsv / step1_delta_v6_vs_v3f.tsv
- step1_trajectory.tsv (per-region 5-class A-E)
- step1_off_vs_on_compare.tsv
- step1_summary_stats.json (effect sizes for H1a/b/c)
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from statistics import mean, median

import pandas as pd
import numpy as np

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"
MASTER = STEP1 / "step1_master_three_way.tsv"

BAMS = ["V3F", "V5", "V6"]
CORE_FEATURES = ["NG", "hp_1_1", "hp_2_1", "hp_33", "n_reads"]
# In master the columns are named "1-1", "2-1", "33"; map to canonical names
COL_MAP = {"hp_1_1": "1-1", "hp_2_1": "2-1", "hp_33": "33", "NG": "NG", "n_reads": "n_reads"}


def col(bam: str, flag: str, feat: str) -> str:
    return f"{bam}_{flag}_{COL_MAP[feat]}"


def wilson_ci(p: float, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    margin = (z / denom) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, center - margin), min(1.0, center + margin))


def safe_ratio(num: float, den: float) -> float:
    if den == 0 or pd.isna(num) or pd.isna(den):
        return float("nan")
    return num / den


def main() -> int:
    print(f"[Step 1.5] loading {MASTER} ...", flush=True)
    df = pd.read_csv(MASTER, sep="\t", low_memory=False)
    print(f"  rows={len(df)}  cols={len(df.columns)}", flush=True)
    print(f"  label distribution: {df['label'].value_counts().to_dict()}", flush=True)

    # Coerce numeric columns
    numeric_prefixes = []
    for bam in BAMS:
        for flag in ["off", "on"]:
            for feat in CORE_FEATURES:
                numeric_prefixes.append(col(bam, flag, feat))
    for c in numeric_prefixes:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["caller_af"] = pd.to_numeric(df["caller_af"], errors="coerce")
    df["Coverage_Multiple"] = pd.to_numeric(df["Coverage_Multiple"], errors="coerce")
    df["Diploid_Coverage_Used"] = pd.to_numeric(df["Diploid_Coverage_Used"], errors="coerce")

    # --- LOH side derivation ----------------------------------------------
    # Prefer master.tsv loh_side; for unjoined rows mark UNKNOWN.
    df["loh_side_final"] = df["loh_side"].where(df["loh_side"].notna() & (df["loh_side"] != ""), "UNKNOWN")
    df["loh_side_final"] = df["loh_side_final"].fillna("UNKNOWN")

    # ----------------------------------------------------------------------
    # Stage 5: per-region deltas for core features (V3F → V5 → V6, flag=off)
    # ----------------------------------------------------------------------
    print("[Step 1.5] computing per-region deltas (flag=off)", flush=True)
    delta_specs = [
        ("v5_vs_v3f", "V5", "V3F"),
        ("v6_vs_v5", "V6", "V5"),
        ("v6_vs_v3f", "V6", "V3F"),
    ]
    for tag, a, b in delta_specs:
        out = STEP1 / f"step1_delta_{tag}.tsv"
        cols_keep = ["region_id", "chr", "pos", "label", "loh_side_final", "caller_af",
                     "LOH_Bed_Annotation", "Coverage_Multiple", "master_join_ok"]
        d = df[cols_keep].copy()
        for feat in CORE_FEATURES:
            d[f"{feat}_{a}"] = df[col(a, "off", feat)]
            d[f"{feat}_{b}"] = df[col(b, "off", feat)]
            d[f"delta_{feat}"] = df[col(a, "off", feat)] - df[col(b, "off", feat)]
            # ratio guard for 1-1:2-1
        d["hp_1_1_2_1_ratio_" + a] = d.apply(lambda r: safe_ratio(r[f"hp_1_1_{a}"], r[f"hp_2_1_{a}"]), axis=1)
        d["hp_1_1_2_1_ratio_" + b] = d.apply(lambda r: safe_ratio(r[f"hp_1_1_{b}"], r[f"hp_2_1_{b}"]), axis=1)
        d["delta_hp_1_1_2_1_ratio"] = d["hp_1_1_2_1_ratio_" + a] - d["hp_1_1_2_1_ratio_" + b]
        d.to_csv(out, sep="\t", index=False, na_rep="NA")
        print(f"  wrote {out.name} ({len(d)} rows)", flush=True)

        # Stats
        for feat in ["NG", "hp_1_1", "hp_2_1", "hp_33"]:
            v = d[f"delta_{feat}"].dropna()
            if len(v) == 0:
                continue
            n_pos = int((v > 0).sum())
            n_neg = int((v < 0).sum())
            n_zero = int((v == 0).sum())
            print(f"    Δ{feat} {tag}: mean={v.mean():.3f} median={v.median():.3f} "
                  f"IQR=[{v.quantile(0.25):.3f},{v.quantile(0.75):.3f}] "
                  f"frac_increase={n_pos/len(v):.3f}", flush=True)

    # ----------------------------------------------------------------------
    # Stage 6: per-region trajectory classes (V3F→V5→V6 on NG)
    # ----------------------------------------------------------------------
    print("[Step 1.5] computing per-region trajectory (NG)...", flush=True)
    traj = df[["region_id", "chr", "pos", "label", "loh_side_final", "caller_af",
               "LOH_Bed_Annotation", "Coverage_Multiple", "master_join_ok"]].copy()
    traj["NG_V3F"] = df[col("V3F", "off", "NG")]
    traj["NG_V5"] = df[col("V5", "off", "NG")]
    traj["NG_V6"] = df[col("V6", "off", "NG")]
    traj["delta_V3F_V5"] = traj["NG_V5"] - traj["NG_V3F"]
    traj["delta_V5_V6"] = traj["NG_V6"] - traj["NG_V5"]
    traj["delta_V3F_V6"] = traj["NG_V6"] - traj["NG_V3F"]

    def classify(d1: float, d2: float) -> str:
        if pd.isna(d1) or pd.isna(d2):
            return "MISSING"
        if d1 > 0 and d2 > 0:
            return "A_both_improve"
        if d1 > 0 and d2 <= 0:
            return "B_only_V5_improve"
        if d1 <= 0 and d2 > 0:
            return "C_only_V6_improve"
        if d1 == 0 and d2 == 0:
            return "D_no_change"
        if d1 < 0 or d2 < 0:
            return "E_reverse_or_decrease"
        return "D_no_change"

    # Standard 5-class: A both improve, B only V5 improve, C only V6 improve, D no change, E reverse
    def classify_strict(d1: float, d2: float) -> str:
        if pd.isna(d1) or pd.isna(d2):
            return "MISSING"
        s1 = 1 if d1 > 0 else (-1 if d1 < 0 else 0)
        s2 = 1 if d2 > 0 else (-1 if d2 < 0 else 0)
        if s1 == 1 and s2 == 1:
            return "A_both_improve"
        if s1 == 1 and s2 != 1:
            return "B_only_V5_improve"
        if s1 != 1 and s2 == 1:
            return "C_only_V6_improve"
        if s1 == 0 and s2 == 0:
            return "D_no_change"
        return "E_reverse_or_decrease"

    traj["class"] = traj.apply(lambda r: classify_strict(r["delta_V3F_V5"], r["delta_V5_V6"]), axis=1)
    traj_out = STEP1 / "step1_trajectory.tsv"
    traj.to_csv(traj_out, sep="\t", index=False, na_rep="NA")
    print(f"  wrote {traj_out.name} ({len(traj)} rows)", flush=True)

    # Class proportions
    class_counts = traj["class"].value_counts().to_dict()
    total = sum(class_counts.values())
    print("  class distribution (NG trajectory):")
    for k, v in sorted(class_counts.items()):
        print(f"    {k}: n={v} ({v/total:.1%})", flush=True)

    # By label
    class_by_label = traj.groupby(["label", "class"]).size().unstack(fill_value=0)
    print("  by label:")
    print(class_by_label.to_string(), flush=True)

    # ----------------------------------------------------------------------
    # Stage 7: off / on flag compare — for canonical phasing-weak cells
    # ----------------------------------------------------------------------
    print("[Step 1.5] computing off/on flag comparison (V5 vs V6 mask check)...", flush=True)
    cells = []
    # Use NG=2 inner same_hp1 (no LOH info before master join; approximate by loh_side_final)
    for label in ["TP", "FP"]:
        sub = df[df["label"] == label]
        n_tot = len(sub)
        cell_rows = []
        # 1. NG_off == 2 cells — global
        for bam in BAMS:
            for flag in ["off", "on"]:
                ng_col = col(bam, flag, "NG")
                n = int(sub[ng_col].notna().sum())
                ng2 = int((sub[ng_col] == 2).sum())
                ng3p = int((sub[ng_col] >= 3).sum())
                cell_rows.append({
                    "label": label, "BAM": bam, "flag": flag,
                    "n_total": n,
                    "NG_eq_2": ng2,
                    "NG_eq_2_frac": ng2 / n if n else float("nan"),
                    "NG_ge_3": ng3p,
                    "NG_ge_3_frac": ng3p / n if n else float("nan"),
                })
        cells.extend(cell_rows)

    # LOH stratified — only on master-joined rows
    df_joined = df[df["master_join_ok"] == 1]
    for label in ["TP", "FP"]:
        for side in sorted(df_joined["loh_side_final"].dropna().unique()):
            sub = df_joined[(df_joined["label"] == label) & (df_joined["loh_side_final"] == side)]
            if len(sub) < 10:
                continue
            for bam in BAMS:
                for flag in ["off", "on"]:
                    ng_col = col(bam, flag, "NG")
                    n = int(sub[ng_col].notna().sum())
                    ng2 = int((sub[ng_col] == 2).sum())
                    ng3p = int((sub[ng_col] >= 3).sum())
                    cells.append({
                        "label": label, "loh_side": side, "BAM": bam, "flag": flag,
                        "n_total": n,
                        "NG_eq_2": ng2,
                        "NG_eq_2_frac": ng2 / n if n else float("nan"),
                        "NG_ge_3": ng3p,
                        "NG_ge_3_frac": ng3p / n if n else float("nan"),
                    })

    cells_df = pd.DataFrame(cells)
    cells_out = STEP1 / "step1_off_vs_on_compare.tsv"
    cells_df.to_csv(cells_out, sep="\t", index=False, na_rep="NA")
    print(f"  wrote {cells_out.name} ({len(cells_df)} rows)", flush=True)

    # ----------------------------------------------------------------------
    # H1a/H1b/H1c judgment
    # NG=2 marker is a characterization handle; we compute TP-rate proxy at NG_off=2 cell.
    # Spec: "Inner same-hap TP gap" — proxy = (TP rate in NG=2 cells) − (TP rate in NG=3+ cells)
    # because phaseC has no per-cell HP-direction grain. Use master-joined LOH inner only.
    # ----------------------------------------------------------------------
    print("[Step 1.5] computing H1a/H1b/H1c judgment ...", flush=True)
    joined = df[df["master_join_ok"] == 1].copy()
    # Define Inner via loh_side or LOH_Bed_Annotation
    joined["loh_inner"] = joined["loh_side_final"].astype(str).str.lower().str.contains("inner")
    joined["loh_outer"] = joined["loh_side_final"].astype(str).str.lower().str.contains("outer")

    summary: dict = {}
    inner_mask = joined["loh_inner"]
    outer_mask = joined["loh_outer"]

    def tp_rate_at_ng(g: pd.DataFrame, ng_col_name: str, ng_eq: int) -> tuple[float, int, int]:
        m = g[ng_col_name] == ng_eq
        if m.sum() == 0:
            return (float("nan"), 0, 0)
        n_tp = int(((g["label"] == "TP") & m).sum())
        n_fp = int(((g["label"] == "FP") & m).sum())
        n_total = n_tp + n_fp
        return (n_tp / n_total if n_total else float("nan"), n_tp, n_fp)

    bam_metrics = {}
    for bam in BAMS:
        ng_col_name = col(bam, "off", "NG")
        inner_rate_ng2, in_tp, in_fp = tp_rate_at_ng(joined[inner_mask], ng_col_name, 2)
        outer_rate_ng2, ou_tp, ou_fp = tp_rate_at_ng(joined[outer_mask], ng_col_name, 2)
        gap = (inner_rate_ng2 - outer_rate_ng2) if (not pd.isna(inner_rate_ng2) and not pd.isna(outer_rate_ng2)) else float("nan")
        bam_metrics[bam] = {
            "inner_NG2_tp_rate": inner_rate_ng2,
            "inner_NG2_n_tp": in_tp,
            "inner_NG2_n_fp": in_fp,
            "outer_NG2_tp_rate": outer_rate_ng2,
            "outer_NG2_n_tp": ou_tp,
            "outer_NG2_n_fp": ou_fp,
            "inner_minus_outer_gap": gap,
        }
        # Wilson CI for gap is non-trivial; report point estimate + per-rate CIs
        if not pd.isna(inner_rate_ng2):
            bam_metrics[bam]["inner_NG2_tp_rate_CI95"] = wilson_ci(inner_rate_ng2, in_tp + in_fp)
        if not pd.isna(outer_rate_ng2):
            bam_metrics[bam]["outer_NG2_tp_rate_CI95"] = wilson_ci(outer_rate_ng2, ou_tp + ou_fp)

    summary["bam_metrics"] = bam_metrics

    # H1a: delta(V5 inner-outer gap) - delta(V3F inner-outer gap) ≥ 0.03
    gap_v3f = bam_metrics["V3F"]["inner_minus_outer_gap"]
    gap_v5 = bam_metrics["V5"]["inner_minus_outer_gap"]
    gap_v6 = bam_metrics["V6"]["inner_minus_outer_gap"]
    h1a = gap_v5 - gap_v3f
    h1b = gap_v6 - gap_v5
    h1c = gap_v6 - gap_v3f

    def verdict(delta: float, thr: float) -> str:
        if pd.isna(delta):
            return "UNKNOWN"
        if delta >= thr:
            return "POSITIVE"
        if delta <= -thr:
            return "NEGATIVE_REVERSE"
        return "NEGATIVE"

    summary["H1a"] = {"delta_gap_V5_minus_V3F": h1a, "threshold": 0.03, "verdict": verdict(h1a, 0.03)}
    summary["H1b"] = {"delta_gap_V6_minus_V5": h1b, "threshold": 0.03, "verdict": verdict(h1b, 0.03)}
    summary["H1c"] = {"delta_gap_V6_minus_V3F": h1c, "threshold": 0.06, "verdict": verdict(h1c, 0.06)}

    # Also: overall marker_off (NG_off ≥ 3) TP rate per BAM
    overall = {}
    for bam in BAMS:
        ng_col_name = col(bam, "off", "NG")
        all_data = df.copy()
        n_tp = int(((all_data["label"] == "TP") & (all_data[ng_col_name] >= 3)).sum())
        n_fp = int(((all_data["label"] == "FP") & (all_data[ng_col_name] >= 3)).sum())
        total = n_tp + n_fp
        rate = n_tp / total if total else float("nan")
        ci = wilson_ci(rate, total) if total else (float("nan"), float("nan"))
        overall[bam] = {
            "marker_off_n": total,
            "marker_off_n_tp": n_tp,
            "marker_off_n_fp": n_fp,
            "marker_off_tp_rate": rate,
            "marker_off_tp_rate_CI95": ci,
        }
    summary["genome_marker_off_NGge3"] = overall

    # Diploid stale check
    stale_diploid = int((df["Diploid_Coverage_Used"] == 75.0).sum())
    summary["stale_diploid75_rows"] = stale_diploid
    summary["total_rows"] = int(len(df))
    summary["master_joined_rows"] = int(df["master_join_ok"].sum())

    summary_path = STEP1 / "step1_summary_stats.json"
    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2, default=lambda o: list(o) if isinstance(o, tuple) else str(o))
    print(f"  wrote {summary_path.name}", flush=True)
    print(json.dumps(summary, indent=2, default=lambda o: list(o) if isinstance(o, tuple) else str(o)), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
