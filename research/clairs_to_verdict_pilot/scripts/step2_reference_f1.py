#!/usr/bin/env python3
"""Step 2: Reference F1 computation under hypothetical Verdict filters.

Context:
    Step 1 showed Verdict_Germline attaches only to LowQual records (never PASS),
    so the plan's nominal "remove Verdict_Germline from PASS" filter is trivially
    a no-op (ΔF1 = 0) on this subsample. To give the report a non-trivial
    reference, we enumerate several hypothetical scenarios in addition to the
    nominal one and record the resulting F1.

Scenarios (all derived from step1_verdict_vs_seqc2.tsv within HighConf BED):
    S0  baseline                — PASS set as output by ClairS-TO
    S1  PASS − Verdict_Germline — plan's nominal filter (expected no-op)
    S2  PASS ∩ Verdict_Somatic∪SubclonalSomatic  — aggressive precision-first
    S3  PASS ∪ (LowQual ∩ Verdict_Somatic)       — rescue any Verdict-somatic
        LowQual candidates (expected no-op; listed for completeness)

Outputs:
    data/step2_reference_f1.tsv
    figures/step2_reference_f1_barchart.png
"""
from __future__ import annotations

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


TRUTH_TOTAL_PATH_DEFAULT = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--step1-tsv", required=True)
    p.add_argument("--truth-total", type=int, default=39447,
                   help="Total SEQC2 sSNV truth (default 39447, from Step 1 output)")
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def f1_set(tp: int, fp: int, truth_total: int) -> dict[str, float]:
    fn = truth_total - tp
    if tp + fp == 0:
        prec = 0.0
    else:
        prec = tp / (tp + fp)
    rec = tp / truth_total if truth_total else 0.0
    if prec + rec == 0:
        f1 = 0.0
    else:
        f1 = 2 * prec * rec / (prec + rec)
    return {"tp": tp, "fp": fp, "fn": fn,
            "precision": prec, "recall": rec, "f1": f1}


def evaluate(call_df: pd.DataFrame, truth_total: int, name: str) -> dict:
    # Only in_highconf calls can be TP; outside is always non-evaluable → drop them
    # (matches standard sSNV benchmark practice)
    hc = call_df[call_df["in_highconf"] == 1]
    tp = int((hc["seqc2_label"] == "TP").sum())
    fp = int((hc["seqc2_label"] == "FP").sum())
    m = f1_set(tp, fp, truth_total)
    m["scenario"] = name
    m["n_calls"] = int(len(hc))
    return m


def main() -> None:
    args = parse_args()
    data_dir = os.path.join(args.out_dir, "data") if os.path.isdir(os.path.join(args.out_dir, "data")) else args.out_dir
    fig_dir_candidate = os.path.join(os.path.dirname(data_dir), "figures")
    fig_dir = fig_dir_candidate if os.path.isdir(fig_dir_candidate) else data_dir

    print(f"[step2] loading {args.step1_tsv}")
    df = pd.read_csv(args.step1_tsv, sep="\t", low_memory=False)
    print(f"        {len(df):,} records")

    # PASS set
    pass_df = df[df["is_pass"] == 1].copy()
    lowqual_df = df[df["filter"].str.contains("LowQual", na=False)].copy()
    print(f"        PASS n={len(pass_df):,}   LowQual n={len(lowqual_df):,}")

    results = []

    # S0 baseline
    results.append(evaluate(pass_df, args.truth_total, "S0_baseline_PASS"))

    # S1 remove Verdict_Germline from PASS (expected no-op)
    s1 = pass_df[pass_df["verdict"] != "Verdict_Germline"]
    results.append(evaluate(s1, args.truth_total, "S1_PASS_minus_Germline"))

    # S2 retain only PASS ∩ (Verdict_Somatic ∪ Verdict_SubclonalSomatic)
    s2 = pass_df[pass_df["verdict"].isin(["Verdict_Somatic", "Verdict_SubclonalSomatic"])]
    results.append(evaluate(s2, args.truth_total, "S2_PASS_only_Verdict_Somatic_or_Subclonal"))

    # S3 PASS ∪ (LowQual ∩ Verdict_Somatic) — rescue (expected empty)
    rescue = lowqual_df[lowqual_df["verdict"] == "Verdict_Somatic"]
    s3 = pd.concat([pass_df, rescue], ignore_index=True)
    results.append(evaluate(s3, args.truth_total, "S3_PASS_plus_LowQual_Verdict_Somatic"))

    out = pd.DataFrame(results)
    baseline_f1 = out.loc[out["scenario"] == "S0_baseline_PASS", "f1"].iloc[0]
    out["delta_f1_vs_baseline"] = out["f1"] - baseline_f1

    # Reorder columns for readability
    cols = ["scenario", "n_calls", "tp", "fp", "fn",
            "precision", "recall", "f1", "delta_f1_vs_baseline"]
    out = out[cols]

    out_tsv = os.path.join(data_dir, "step2_reference_f1.tsv")
    out.to_csv(out_tsv, sep="\t", index=False)
    print(f"[step2] wrote {out_tsv}")

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=(9, 5))
    xs = list(range(len(out)))
    bars = ax.bar(xs, out["f1"], color=["#1f77b4", "#2ca02c", "#d62728", "#ff7f0e"],
                  alpha=0.85, edgecolor="black")
    for x, f1, d in zip(xs, out["f1"], out["delta_f1_vs_baseline"]):
        label = f"F1={f1:.4f}"
        if x > 0:
            label += f"\nΔ={d:+.4f}"
        ax.text(x, f1 + 0.01, label, ha="center", fontsize=9)
    ax.axhline(baseline_f1, color="grey", linestyle="--", alpha=0.6, label=f"baseline F1={baseline_f1:.4f}")
    ax.set_xticks(xs)
    ax.set_xticklabels(out["scenario"], rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("F1")
    ax.set_ylim(0, max(out["f1"]) * 1.25)
    ax.set_title("Hypothetical reference F1 under Verdict-based filter scenarios\n"
                 "HCC1395 t20_n30 · SEQC2 HighConf · PASS baseline")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "step2_reference_f1_barchart.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[step2] wrote {fig_path}")

    # ---- Stdout summary ----
    print("\n" + "=" * 70)
    print("Step 2 summary — hypothetical Verdict filter reference F1")
    print("=" * 70)
    for _, r in out.iterrows():
        print(f"  {r['scenario']:<42}  n={r['n_calls']:>6}  TP={r['tp']:>5}  FP={r['fp']:>5}  "
              f"P={r['precision']:.3f}  R={r['recall']:.3f}  F1={r['f1']:.4f}  Δ={r['delta_f1_vs_baseline']:+.4f}")
    print("=" * 70)

    # Gate for the plan's judgment path
    s1_delta = out.loc[out["scenario"] == "S1_PASS_minus_Germline", "delta_f1_vs_baseline"].iloc[0]
    print(f"\n  Plan's nominal filter (S1) ΔF1 = {s1_delta:+.4f}")
    if s1_delta >= 0.01:
        print("  Gate: ≥+0.01 → strong POSITIVE for direct F1 gain")
    elif s1_delta >= 0.005:
        print("  Gate: +0.005..+0.01 → CONDITIONAL")
    else:
        print("  Gate: <+0.005 → NEGATIVE for direct F1 gain (Verdict_Germline already excluded from PASS)")


if __name__ == "__main__":
    main()
