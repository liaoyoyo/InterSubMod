#!/usr/bin/env python3
"""Cycle 3 Step 1 — Apply caller-F1-headroom-aware gated filter on cycle 2 5 samples.

Gate rule (from cycle3_caller_f1_gate.json):
    IF caller_F1 < 0.80 AND FP_density > 0.10:
        action = apply per-sample re-fit ΔF1 (from cycle 2 B3' refit column)
    ELSE:
        action = skip filter (ΔF1 = 0, caller F1 preserved)

Reads cycle 2 cross_sample_delta_f1.tsv to extract per-sample re-fit ΔF1 and
caller F1 baseline; FP density computed from tp_total_used + fp_total_used.

Output:
    cycle3/data/cycle3_step1_gated_delta_f1.tsv
    cycle3/figures/cycle3_step1_gated_delta_f1.png
    cycle3/intermediate/cycle3_step1_summary.json

Pre-reg verdicts:
    H_C3_1: qualifying mean ΔF1 ≥ +0.01 (HCC1395 + HCC1937 expected)
    H_C3_3: non-qualifying |drift| < 0.001 (H1437/H2009/HCC1954 expected drift=0)
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

CYCLE3_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/methyl_augmented_filter_phase2/cycle3")
CYCLE2_DELTA_TSV = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/methyl_augmented_filter_phase2/cycle2/data/cycle2_cross_sample_delta_f1.tsv"
)
GATE_CONFIG = CYCLE3_DIR / "cycle3_caller_f1_gate.json"

OUT_TSV = CYCLE3_DIR / "data" / "cycle3_step1_gated_delta_f1.tsv"
OUT_FIG = CYCLE3_DIR / "figures" / "cycle3_step1_gated_delta_f1.png"
OUT_SUMMARY = CYCLE3_DIR / "intermediate" / "cycle3_step1_summary.json"

SAMPLE_ORDER = ["HCC1395", "H1437", "H2009", "HCC1954", "HCC1937"]


def main():
    gate = json.loads(GATE_CONFIG.read_text())
    f1_max = gate["gate_conditions"]["caller_F1_max"]
    fp_density_min = gate["gate_conditions"]["FP_density_min"]
    pass_thr = gate["decision_thresholds"]["H_C3_1_qualifying_mean_delta_F1_pass"]
    drift_thr = gate["decision_thresholds"]["H_C3_3_non_qualifying_drift_max"]

    df = pd.read_csv(CYCLE2_DELTA_TSV, sep="\t")
    refit = df[df["mode"] == "refit"].set_index("sample")

    rows = []
    for sample in SAMPLE_ORDER:
        r = refit.loc[sample]
        tp = float(r["tp_total_used"])
        fp = float(r["fp_total_used"])
        fp_density = fp / (tp + fp)
        caller_f1 = float(r["caller_F1"])

        gate_pass = (caller_f1 < f1_max) and (fp_density > fp_density_min)

        if gate_pass:
            gated_delta_f1 = float(r["delta_F1"])
            gated_tau = float(r["tau"])
            action = "apply_refit"
            f1_post = float(r["F1_post"])
        else:
            gated_delta_f1 = 0.0
            gated_tau = None
            action = "skip"
            f1_post = caller_f1

        rows.append({
            "sample": sample,
            "caller_F1": caller_f1,
            "fp_density": fp_density,
            "tp_total": int(tp),
            "fp_total": int(fp),
            "gate_caller_F1_lt_080": caller_f1 < f1_max,
            "gate_fp_density_gt_010": fp_density > fp_density_min,
            "gate_pass": gate_pass,
            "action": action,
            "gated_tau": gated_tau,
            "gated_F1_post": f1_post,
            "gated_delta_F1": gated_delta_f1,
            "refit_delta_F1_reference": float(r["delta_F1"]),
        })

    out = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUT_TSV}")

    qualifying = out[out["gate_pass"]]
    non_qualifying = out[~out["gate_pass"]]

    qual_mean = float(qualifying["gated_delta_F1"].mean()) if len(qualifying) else float("nan")
    qual_n = int(len(qualifying))
    qual_samples = qualifying["sample"].tolist()

    non_qual_max_abs_drift = float(non_qualifying["gated_delta_F1"].abs().max()) if len(non_qualifying) else 0.0
    non_qual_n = int(len(non_qualifying))
    non_qual_samples = non_qualifying["sample"].tolist()

    h_c3_1_verdict = "PASS" if qual_mean >= pass_thr else ("FALSIFY" if qual_mean < gate["decision_thresholds"]["H_C3_1_qualifying_mean_delta_F1_falsify"] else "MARGINAL")
    h_c3_3_verdict = "PASS" if non_qual_max_abs_drift < drift_thr else "FAIL"

    summary = {
        "cycle_id": gate["cycle_id"],
        "gate_rule": {
            "caller_F1_max": f1_max,
            "FP_density_min": fp_density_min,
        },
        "qualifying": {
            "n": qual_n,
            "samples": qual_samples,
            "mean_gated_delta_F1": qual_mean,
            "threshold_pass": pass_thr,
        },
        "non_qualifying": {
            "n": non_qual_n,
            "samples": non_qual_samples,
            "max_abs_drift": non_qual_max_abs_drift,
            "threshold_drift_max": drift_thr,
        },
        "verdicts": {
            "H_C3_1_qualifying_mean_delta_F1": {
                "value": qual_mean,
                "threshold": pass_thr,
                "verdict": h_c3_1_verdict,
            },
            "H_C3_3_non_qualifying_drift": {
                "value": non_qual_max_abs_drift,
                "threshold": drift_thr,
                "verdict": h_c3_3_verdict,
            },
            "H_C3_2_panel_expansion": {
                "verdict": "PENDING_STEP_2",
                "note": "Panel expansion to n≥4 low-F1 samples blocked until Step 2 survey complete",
            },
        },
        "naive_5sample_mean_for_comparison": float(out["gated_delta_F1"].mean()),
    }

    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"Wrote {OUT_SUMMARY}")

    # Figure: bar chart 5 sample gated ΔF1 + threshold lines
    fig, ax = plt.subplots(figsize=(10, 5.5))
    colors = ["#2E7D32" if g else "#9E9E9E" for g in out["gate_pass"]]
    bars = ax.bar(out["sample"], out["gated_delta_F1"], color=colors, edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(pass_thr, color="#1976D2", linestyle="--", linewidth=1, label=f"H_C3_1 pass = +{pass_thr}")
    ax.axhline(drift_thr, color="#D32F2F", linestyle=":", linewidth=1, label=f"H_C3_3 drift cap = ±{drift_thr}")
    ax.axhline(-drift_thr, color="#D32F2F", linestyle=":", linewidth=1)

    for i, (sample, df1, gate_pass) in enumerate(zip(out["sample"], out["gated_delta_F1"], out["gate_pass"])):
        label = f"{df1:+.5f}"
        if gate_pass:
            label += "\n(apply re-fit)"
        else:
            label += "\n(skip)"
        ax.text(i, df1 + (0.0015 if df1 >= 0 else -0.0015), label,
                ha="center", va="bottom" if df1 >= 0 else "top", fontsize=8.5)

    ax.set_ylabel("Gated ΔF1 (F1_post − caller_F1)", fontsize=11)
    ax.set_title("Cycle 3 Step 1 — Caller-F1-headroom-aware gated filter\n"
                 "Green = qualify (caller F1<0.80 AND FP density>0.10) | Grey = skip",
                 fontsize=11)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=120)
    plt.close(fig)
    print(f"Wrote {OUT_FIG}")

    # Console verdict
    print()
    print(f"=== Cycle 3 Step 1 Verdicts ===")
    print(f"  Qualifying samples (n={qual_n}): {qual_samples}")
    print(f"  Qualifying mean ΔF1 = {qual_mean:+.5f} (threshold +{pass_thr})  → H_C3_1: {h_c3_1_verdict}")
    print(f"  Non-qualifying samples (n={non_qual_n}): {non_qual_samples}")
    print(f"  Non-qualifying max |drift| = {non_qual_max_abs_drift:.6f} (cap {drift_thr})  → H_C3_3: {h_c3_3_verdict}")
    print(f"  H_C3_2 panel expansion: PENDING_STEP_2")

    return 0


if __name__ == "__main__":
    sys.exit(main())
