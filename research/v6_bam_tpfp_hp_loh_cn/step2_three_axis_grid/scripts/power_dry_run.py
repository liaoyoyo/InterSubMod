"""Step 1.5 Power dry-run gate — see 00_PLAN.md §Step 2 / 02_methodology_notes §2.

Builds the 3-axis (LOH × HP × CN) preview, counts cells by power tier, audits
collinearity, and writes power_dry_run.md + power_dry_run.tsv. Gate threshold:
powered cells (n≥50) ≥ 15 → proceed to full 50-cell grid; otherwise downscale.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _common import (
    COV_LEVELS,
    HP_LEVELS,
    HP_PREFIX_PRIMARY,
    HP_PREFIX_SENSITIVITY,
    LOH_LEVELS,
    MARGINAL_THRESHOLD,
    POWERED_THRESHOLD,
    STEP2_DIR,
    annotate_axes,
    cell_id,
    joined_subset,
    load_master,
    log_msg,
    wilson_ci,
)


def cramers_v(a: pd.Series, b: pd.Series) -> float:
    cont = pd.crosstab(a, b)
    chi2 = stats.chi2_contingency(cont)[0]
    n = cont.values.sum()
    r, k = cont.shape
    if min(r - 1, k - 1) == 0 or n == 0:
        return float("nan")
    return float(np.sqrt(chi2 / (n * min(r - 1, k - 1))))


def vif_for(df: pd.DataFrame, target: str, features: list[str]) -> float:
    """Simple R^2-based VIF: regress target on features (numeric, OHE for categorical)."""
    from sklearn.linear_model import LinearRegression

    X = pd.get_dummies(df[features], drop_first=True).astype(float).values
    y = df[target].astype(float).values
    if X.shape[1] == 0:
        return float("nan")
    reg = LinearRegression().fit(X, y)
    yhat = reg.predict(X)
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return 1.0 / max(1e-6, 1.0 - r2)


def main():
    log_msg("=== Step 1.5 power dry-run START ===")

    df_all = load_master()
    log_msg(f"Loaded master TSV: {len(df_all):,} rows × {len(df_all.columns)} cols")

    df = joined_subset(df_all)
    log_msg(
        f"Subset master_join_ok==1: {len(df):,} rows "
        f"(TP={int((df['label']=='TP').sum()):,}, FP={int((df['label']=='FP').sum()):,})"
    )

    # ---- Diploid_Coverage_Used audit ----
    dcu_median = float(df["Diploid_Coverage_Used"].median())
    dcu_audit = {
        "median": dcu_median,
        "mean": float(df["Diploid_Coverage_Used"].mean()),
        "expected_HCC1395_TO_KDE": 61.0,
        "actual_paired_full_value": 115.0,
        "matches_HCC1395_TO_KDE_61": abs(dcu_median - 61.0) < 1.0,
        "flag_stale_75": abs(dcu_median - 75.0) < 1.0,
    }
    log_msg(f"Diploid_Coverage_Used median = {dcu_median} (HCC1395 paired_full baseline ≈ 115)")

    # ---- Annotate axes (primary V6_off) ----
    df = annotate_axes(df, hp_prefix=HP_PREFIX_PRIMARY)
    df_on = annotate_axes(joined_subset(df_all), hp_prefix=HP_PREFIX_SENSITIVITY)

    # ---- Build cell counts (full 75 raw combos, conceptual 50 after pruning UNKNOWN) ----
    rows = []
    for loh in LOH_LEVELS:
        for hp in HP_LEVELS:
            for cov in COV_LEVELS:
                cell = df[
                    (df["loh_side_norm"] == loh)
                    & (df["hp_bucket"] == hp)
                    & (df["cov_bin"] == cov)
                ]
                cell_on = df_on[
                    (df_on["loh_side_norm"] == loh)
                    & (df_on["hp_bucket"] == hp)
                    & (df_on["cov_bin"] == cov)
                ]
                n = len(cell)
                n_tp = int((cell["label"] == "TP").sum())
                n_fp = int((cell["label"] == "FP").sum())
                lo, hi = wilson_ci(n_tp, n) if n > 0 else (float("nan"), float("nan"))
                rows.append(
                    {
                        "cell_id": cell_id(loh, hp, cov),
                        "loh_side": loh,
                        "hp_bucket": hp,
                        "cov_bin": cov,
                        "n": n,
                        "n_TP": n_tp,
                        "n_FP": n_fp,
                        "TP_rate": (n_tp / n) if n else float("nan"),
                        "wilson_lo": lo,
                        "wilson_hi": hi,
                        "wilson_width": (hi - lo) if n else float("nan"),
                        "powered": n >= POWERED_THRESHOLD,
                        "marginal": MARGINAL_THRESHOLD <= n < POWERED_THRESHOLD,
                        "underpowered": 0 < n < MARGINAL_THRESHOLD,
                        "empty": n == 0,
                        "n_on_sensitivity": len(cell_on),
                    }
                )

    grid = pd.DataFrame(rows)

    n_non_empty = int((grid["n"] >= 1).sum())
    n_underpow = int(grid["underpowered"].sum())
    n_marg = int(grid["marginal"].sum())
    n_powered = int(grid["powered"].sum())
    n_unknown = int((grid["loh_side"] == "UNKNOWN").sum())

    # Cells excluding UNKNOWN loh row (conceptual 50)
    grid_main = grid[grid["loh_side"] != "UNKNOWN"].copy()
    n_main_powered = int(grid_main["powered"].sum())
    n_main_marg = int(grid_main["marginal"].sum())
    n_main_non_empty = int((grid_main["n"] >= 1).sum())

    log_msg(
        f"Cells (75 raw / 50 main): non_empty={n_non_empty}, "
        f"powered={n_powered} (main={n_main_powered}), marginal={n_marg}, underpowered={n_underpow}"
    )

    # ---- Collinearity audit (5 axes incl. NG/AF/NumReads covariates) ----
    df["AF_bin"] = pd.cut(
        df["caller_af"], bins=[-0.01, 0.1, 0.25, 0.4, 0.6, 1.01], labels=["af1", "af2", "af3", "af4", "af5"]
    ).astype(object)
    df["NR_bin"] = pd.qcut(df["n_reads_v6off"].clip(lower=1), q=5, duplicates="drop").astype(str)
    df["NG_bin"] = df["NG_v6off"].astype(str)

    pairs = [
        ("loh_side_norm", "hp_bucket"),
        ("loh_side_norm", "cov_bin"),
        ("hp_bucket", "cov_bin"),
        ("loh_side_norm", "NG_bin"),
        ("hp_bucket", "NG_bin"),
        ("cov_bin", "NG_bin"),
        ("loh_side_norm", "AF_bin"),
        ("hp_bucket", "AF_bin"),
        ("cov_bin", "AF_bin"),
        ("NG_bin", "AF_bin"),
        ("loh_side_norm", "NR_bin"),
        ("hp_bucket", "NR_bin"),
        ("cov_bin", "NR_bin"),
        ("NG_bin", "NR_bin"),
        ("AF_bin", "NR_bin"),
    ]
    collinear = []
    for a, b in pairs:
        v = cramers_v(df[a], df[b])
        collinear.append({"axis_a": a, "axis_b": b, "cramers_v": v})
    coll_df = pd.DataFrame(collinear).sort_values("cramers_v", ascending=False)

    # VIF for numeric covariates against the 3 grid axes
    vif_records = []
    feats_grid = ["loh_side_norm", "hp_bucket", "cov_bin"]
    for tgt in ["NG_v6off", "caller_af", "n_reads_v6off"]:
        vif_records.append({"target": tgt, "vif_vs_grid": vif_for(df, tgt, feats_grid)})
    vif_df = pd.DataFrame(vif_records)

    # ---- Decision ----
    gate_pass = n_main_powered >= 15
    if gate_pass:
        decision = "PROCEED full 3-axis (50 main cells, plus UNKNOWN row as auxiliary)"
        downscale = None
    else:
        # Fallback: 2-axis (LOH × HP) — 10 cells (excluding UNKNOWN)
        decision = "DOWNSCALE to 2-axis (LOH × HP) — 10 cells"
        downscale = "loh_x_hp_10cell"

    # ---- Write outputs ----
    STEP2_DIR.mkdir(parents=True, exist_ok=True)
    out_tsv = STEP2_DIR / "power_dry_run.tsv"
    grid.to_csv(out_tsv, sep="\t", index=False, float_format="%.6f")
    log_msg(f"Wrote {out_tsv} ({out_tsv.stat().st_size:,} bytes)")

    coll_df.to_csv(STEP2_DIR / "intermediate" / "collinearity_audit.tsv", sep="\t", index=False, float_format="%.4f")
    vif_df.to_csv(STEP2_DIR / "intermediate" / "vif_audit.tsv", sep="\t", index=False, float_format="%.4f")

    summary = {
        "n_rows_total": int(len(df_all)),
        "n_rows_joined": int(len(df)),
        "n_tp_joined": int((df["label"] == "TP").sum()),
        "n_fp_joined": int((df["label"] == "FP").sum()),
        "diploid_coverage_used_audit": dcu_audit,
        "cells_75_raw": {
            "non_empty": n_non_empty,
            "empty": int(grid["empty"].sum()),
            "powered": n_powered,
            "marginal": n_marg,
            "underpowered": n_underpow,
        },
        "cells_50_main_no_unknown": {
            "non_empty": n_main_non_empty,
            "powered": n_main_powered,
            "marginal": n_main_marg,
        },
        "unknown_loh_cells": n_unknown,
        "gate_threshold_main_powered": 15,
        "gate_pass": bool(gate_pass),
        "decision": decision,
        "downscale_option": downscale,
        "fp_constraint_note": (
            "Master-join restricts FP universe to 506 / 4842 rows. FP analyses limited to this subset."
        ),
    }
    with open(STEP2_DIR / "intermediate" / "power_dry_run_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ---- Markdown report ----
    md_lines = []
    md_lines.append("# Step 1.5 — Power dry-run gate")
    md_lines.append("")
    md_lines.append(f"- Master rows (joined / total): **{len(df):,} / {len(df_all):,}**")
    md_lines.append(
        f"- TP / FP in joined subset: **{int((df['label']=='TP').sum()):,} / {int((df['label']=='FP').sum()):,}**"
    )
    md_lines.append("- Primary BAM/flag: **V6, --germline-hp-only=off**")
    md_lines.append("- HP bucket prefix: `V6_off` (HP1, HP1-1, HP2, HP2-1)")
    md_lines.append("")
    md_lines.append("## Diploid_Coverage_Used audit")
    md_lines.append("")
    md_lines.append(f"- Median = **{dcu_median:.1f}** (mean {dcu_audit['mean']:.2f})")
    md_lines.append("- HCC1395 paired_full pipeline baseline ≈ 115 (TO-mode KDE baseline would be 61)")
    md_lines.append(
        f"- Stale 75.0 flag triggered: **{dcu_audit['flag_stale_75']}** "
        f"(matches paired_full 115: **{abs(dcu_median - 115.0) < 1.0}**)"
    )
    md_lines.append("")
    md_lines.append("Note: Step 1 master is built on `paired_full` mode KDE-corrected master.tsv.gz "
                    "(see 01_data_provenance §7). The 61 value mentioned in 02_methodology_notes §1.1 "
                    "refers to HCC1395 TO mode; paired_full uses different diploid coverage (~115). "
                    "Coverage_Multiple is still scale-invariant ratio (Step 1 verified no 75.0 stale rows).")
    md_lines.append("")
    md_lines.append("## Cell power summary")
    md_lines.append("")
    md_lines.append(f"- Raw 75-cell grid (3 LOH × 5 HP × 5 CN incl. UNKNOWN LOH):")
    md_lines.append(f"  - Non-empty: **{n_non_empty}**  /  Empty: **{int(grid['empty'].sum())}**")
    md_lines.append(
        f"  - Powered (n≥50): **{n_powered}**  /  Marginal (30-49): **{n_marg}**  /  "
        f"Underpowered (1-29): **{n_underpow}**"
    )
    md_lines.append(f"- Main 50-cell grid (LOH ∈ {{Inner, Outer}}):")
    md_lines.append(
        f"  - Non-empty: **{n_main_non_empty}**  Powered: **{n_main_powered}**  Marginal: **{n_main_marg}**"
    )
    md_lines.append(f"- UNKNOWN LOH cells (auxiliary): **{n_unknown}**")
    md_lines.append("")
    md_lines.append("## Gate decision")
    md_lines.append("")
    md_lines.append(f"- Threshold: main_powered ≥ 15")
    md_lines.append(f"- Observed main_powered = **{n_main_powered}**")
    md_lines.append(f"- **Decision: {decision}**")
    md_lines.append("")
    md_lines.append("## FP universe constraint (carried from Step 1)")
    md_lines.append("")
    md_lines.append("- FP rows in master TSV: 4,842 total → only **506 (10.5%)** join master.tsv.gz")
    md_lines.append("- FP rate / log-odds metrics are restricted to this 506-FP subset")
    md_lines.append("- TP universe is robust: 29,530 / 30,490 joined (96.8%)")
    md_lines.append("")
    md_lines.append("## Collinearity audit (Cramér's V, top 10)")
    md_lines.append("")
    md_lines.append("| axis_a | axis_b | Cramér's V |")
    md_lines.append("|---|---|---|")
    for _, r in coll_df.head(10).iterrows():
        md_lines.append(f"| {r['axis_a']} | {r['axis_b']} | {r['cramers_v']:.3f} |")
    md_lines.append("")
    md_lines.append("Full audit: `intermediate/collinearity_audit.tsv`")
    md_lines.append("")
    md_lines.append("VIF (R²-based, target ≤ 5 considered safe):")
    md_lines.append("")
    md_lines.append("| target | VIF vs (loh, hp, cov) |")
    md_lines.append("|---|---|")
    for _, r in vif_df.iterrows():
        md_lines.append(f"| {r['target']} | {r['vif_vs_grid']:.2f} |")
    md_lines.append("")
    md_lines.append("## Top-15 most populated cells (preview)")
    md_lines.append("")
    md_lines.append("| cell_id | n | n_TP | n_FP | TP_rate | Wilson 95% CI |")
    md_lines.append("|---|---|---|---|---|---|")
    top = grid.sort_values("n", ascending=False).head(15)
    for _, r in top.iterrows():
        if r["n"] == 0:
            continue
        md_lines.append(
            f"| {r['cell_id']} | {int(r['n']):,} | {int(r['n_TP']):,} | {int(r['n_FP'])} | "
            f"{r['TP_rate']:.4f} | [{r['wilson_lo']:.4f}, {r['wilson_hi']:.4f}] |"
        )
    md_lines.append("")
    md_lines.append("---")
    md_lines.append(
        "Generated by `step2_three_axis_grid/scripts/power_dry_run.py`. "
        "Build log: `intermediate/build_log.txt`. JSON summary: `intermediate/power_dry_run_summary.json`."
    )

    (STEP2_DIR / "power_dry_run.md").write_text("\n".join(md_lines) + "\n")
    log_msg(f"Wrote power_dry_run.md")
    log_msg(f"=== Gate decision: {decision} ===")
    log_msg("=== Step 1.5 power dry-run END ===")

    return summary


if __name__ == "__main__":
    main()
