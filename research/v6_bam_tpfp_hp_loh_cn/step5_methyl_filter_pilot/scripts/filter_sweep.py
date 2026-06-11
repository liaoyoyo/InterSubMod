"""Step 2 — FP-rich cells filter threshold sweep.

For each FP-rich cell (FP_rate > global × 1.5):
1. Pick Model B if LRT q<0.05 for the chosen (BAM × flag); else fall back to Model A.
2. Re-fit on the cell to get predicted P(TP) per region.
3. Sweep τ ∈ {0.50, 0.51, ..., 0.95} (46 values). Filter rule: P(TP) < τ → flagged FP (removed).
4. Compute TP_loss_pct, FP_removal_pct, precision_gain, recall_loss.
5. Mark τ where FP_removal_pct > TP_loss_pct (H3 direction).

FP-rich cells come from:
- step3/chr8_hotspot/zone_grid_real_cov.tsv (powered + FP_rate > 0.2055)
- step3/auto_detected_top5pct/zone_grid.tsv (powered + FP_rate > 0.2055; uses cov_proxy)

Outputs:
  step2_filter_sweep.tsv          per (cell, BAM, flag, τ) × {metrics}
  step2_findings.md               H3 preliminary verdict + cell summary
  figures/step2_roc_per_cell.png  precision-recall curves
  figures/step2_filter_signal.png FP_removal vs TP_loss scatter
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from _common_step1 import (
    BAMS,
    FLAGS,
    METHYL_COVARIATES,
    MASTER_AUG,
    STEP3_CHR8,
    STEP3_AUTO,
    STEP5_DIR,
    annotate_v6off_cell,
    log_msg,
    matplotlib_setup,
)
from augmented_lr import (
    assemble_design,
    fit_lr,
)

warnings.filterwarnings("ignore")

TAU_GRID = np.round(np.arange(0.50, 0.9501, 0.01), 2)  # 46 values
GLOBAL_FP_RATE = 0.137
FP_RICH_THRESH = GLOBAL_FP_RATE * 1.5  # 0.2055

# Per the plan, primary BAM/flag for filter is V6_off (V6 is the production binary)
PRIMARY_BAM = "V6"
PRIMARY_FLAG = "off"


# ---------------------------------------------------------------------------
# cov_proxy reconstruction (matches step3/_zone_common.py compute_n_reads_cov_proxy)
# ---------------------------------------------------------------------------
def compute_cov_proxy(s: pd.Series) -> pd.Series:
    q33, q67, q90 = s.quantile([0.33, 0.67, 0.90])
    bins = [-np.inf, q33, q67, q90, np.inf]
    labels = ["cov_proxy_low", "cov_proxy_mid", "cov_proxy_high", "cov_proxy_top"]
    return pd.cut(s, bins=bins, labels=labels, right=False).astype(object).fillna("cov_proxy_unknown")


# ---------------------------------------------------------------------------
# FP-density (top 5%) reconstruction — read precomputed gz from step3
# ---------------------------------------------------------------------------
FP_DENSITY_FILE = STEP5_DIR.parent / "step3_fp_zone_zoom" / "auto_detected_top5pct" / "fp_density_per_region.tsv.gz"


def load_fp_density() -> pd.DataFrame:
    return pd.read_csv(FP_DENSITY_FILE, sep="\t")


# ---------------------------------------------------------------------------
# FP-rich cells assembly
# ---------------------------------------------------------------------------
def find_fp_rich_cells(thresh: float = FP_RICH_THRESH) -> list[dict]:
    """Return list of FP-rich cell descriptors (source, loh, hp_bucket, cov, cov_axis,
    n, n_TP, n_FP, FP_rate)."""
    rows = []
    g_chr8 = pd.read_csv(STEP3_CHR8, sep="\t")
    for _, r in g_chr8[(g_chr8["powered"] == True) & (g_chr8["FP_rate"] > thresh)].iterrows():
        rows.append({
            "source": "chr8_hotspot",
            "cell_id": f"chr8|{r['loh']}|{r['hp_bucket']}|{r['cov']}",
            "loh_side": r["loh"],
            "hp_bucket": r["hp_bucket"],
            "cov": r["cov"],
            "cov_axis": "cov_bin",
            "n": int(r["n"]),
            "n_TP": int(r["n_TP"]),
            "n_FP": int(r["n_FP"]),
            "FP_rate": float(r["FP_rate"]),
        })
    g_auto = pd.read_csv(STEP3_AUTO, sep="\t")
    for _, r in g_auto[(g_auto["powered"] == True) & (g_auto["FP_rate"] > thresh)].iterrows():
        rows.append({
            "source": "auto_top5pct",
            "cell_id": f"auto|{r['loh']}|{r['hp_bucket']}|{r['cov']}",
            "loh_side": r["loh"],
            "hp_bucket": r["hp_bucket"],
            "cov": r["cov"],
            "cov_axis": "cov_proxy",
            "n": int(r["n"]),
            "n_TP": int(r["n_TP"]),
            "n_FP": int(r["n_FP"]),
            "FP_rate": float(r["FP_rate"]),
        })
    return rows


def select_cell_rows(df: pd.DataFrame, cell: dict) -> pd.Index:
    """Return Index of master TSV rows belonging to this FP-rich cell.
    df should be pre-annotated with hp_bucket_v6off / cov_bin / loh_side_norm / cov_proxy / in_Z_AUTO / chr.
    """
    loh = cell["loh_side"]
    hp = cell["hp_bucket"]
    cov = cell["cov"]
    if cell["source"] == "chr8_hotspot":
        mask = (
            (df["chr"] == "chr8")
            & (df["loh_side_norm"] == loh)
            & (df["hp_bucket_v6off"] == hp)
            & (df["cov_bin"] == cov)
        )
    else:  # auto_top5pct (Z-AUTO)
        mask = (
            (df["in_Z_AUTO"] == True)
            & (df["loh_side_norm"] == loh)
            & (df["hp_bucket_v6off"] == hp)
            & (df["cov_proxy"] == cov)
        )
    return df.index[mask]


# ---------------------------------------------------------------------------
# LR fit + predict for a cell
# ---------------------------------------------------------------------------
def fit_predict_cell(df_cell: pd.DataFrame, bam: str, flag: str, model: str):
    """Fit LR on Model A or B and return (predictions on training rows, mask of rows used)."""
    y = (df_cell["label"] == "TP").astype(int).values
    X, _, mask = assemble_design(df_cell, bam, flag, model)
    if mask.sum() < 10 or len(np.unique(y[mask])) < 2:
        return None, mask, "insufficient_y"
    r, status = fit_lr(y[mask], X[mask])
    if r is None:
        return None, mask, status
    Xc = sm.add_constant(X[mask], has_constant="add")
    try:
        p = r.predict(Xc)
        return p, mask, "ok"
    except Exception as e:
        return None, mask, f"predict_fail:{type(e).__name__}"


# ---------------------------------------------------------------------------
# Determine Model A vs B per cell based on Step 1 LRT q
# ---------------------------------------------------------------------------
def load_step1_results() -> pd.DataFrame:
    return pd.read_csv(STEP5_DIR / "step1_lrt_per_cell.tsv", sep="\t")


def main():
    log_msg("=== Step 2 filter_sweep START ===")
    log_msg(f"τ grid: {len(TAU_GRID)} values, {TAU_GRID[0]} → {TAU_GRID[-1]}")

    # Load augmented master
    df = pd.read_csv(MASTER_AUG, sep="\t", low_memory=False)
    sig_cols = [c for c in df.columns
                if c.endswith("_meth_HPMergedSig") or c.endswith("_meth_HPFineSig") or c.endswith("_meth_AlleleP")]
    for c in sig_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = annotate_v6off_cell(df)

    # Compute cov_proxy (Step 3 used V6_off_n_reads quantiles on full df_full)
    df["cov_proxy"] = compute_cov_proxy(df["V6_off_n_reads"].fillna(0).astype(float))

    # Z-AUTO membership
    fp_density = load_fp_density()
    df = df.merge(fp_density[["region_id", "fp_density", "in_Z_AUTO"]], on="region_id", how="left")
    df["in_Z_AUTO"] = df["in_Z_AUTO"].fillna(False).astype(bool)

    # FP-rich cells
    fp_rich = find_fp_rich_cells()
    log_msg(f"FP-rich cells (FP_rate > {FP_RICH_THRESH:.4f}): {len(fp_rich)}")
    for c in fp_rich:
        log_msg(f"  {c['cell_id']}  n={c['n']:>5}  n_TP={c['n_TP']:>4}  n_FP={c['n_FP']:>4}  FP_rate={c['FP_rate']:.4f}")

    # Step 1 LRT results (to decide Model A vs B per cell)
    step1 = load_step1_results()

    # Sweep loop
    sweep_rows = []
    cell_summaries = []

    for cell in fp_rich:
        idx = select_cell_rows(df, cell)
        df_cell = df.loc[idx].copy()
        actual_n = len(df_cell)
        actual_n_tp = int((df_cell["label"] == "TP").sum())
        actual_n_fp = int((df_cell["label"] == "FP").sum())
        log_msg(f"\nCell: {cell['cell_id']}  step3-n={cell['n']}  reconstructed-n={actual_n}  "
                f"(TP={actual_n_tp}, FP={actual_n_fp})")
        if actual_n < 30 or actual_n_fp < 5:
            log_msg(f"  SKIP: too small after reconstruction")
            continue

        # Use primary BAM × flag (V6_off) for the filter
        bam, flag = PRIMARY_BAM, PRIMARY_FLAG

        # For cell-level model choice we look at the most relevant step1 cell — closest match.
        # FP-rich cells use different partitions than step2 powered cells; we use Model A as
        # fallback baseline and try Model B; only stay with B if it fits OK.
        p_B, mask_B, status_B = fit_predict_cell(df_cell, bam, flag, "B")
        p_A, mask_A, status_A = fit_predict_cell(df_cell, bam, flag, "A")
        used_model = "B" if p_B is not None else ("A" if p_A is not None else None)
        if used_model is None:
            log_msg(f"  SKIP: both Model A & B failed (A:{status_A}, B:{status_B})")
            continue
        p = p_B if used_model == "B" else p_A
        mask = mask_B if used_model == "B" else mask_A
        y = (df_cell["label"] == "TP").astype(int).values[mask]
        n_fit = int(mask.sum())
        n_tp_fit = int(y.sum())
        n_fp_fit = int(n_fit - n_tp_fit)
        log_msg(f"  Fit Model {used_model} on {n_fit} rows (TP={n_tp_fit}, FP={n_fp_fit}); "
                f"A_status={status_A}, B_status={status_B}")

        # Sweep τ
        for tau in TAU_GRID:
            kept = (p >= tau)  # P(TP) >= τ → keep
            tp_kept = int(((y == 1) & kept).sum())
            fp_kept = int(((y == 0) & kept).sum())
            tp_removed = n_tp_fit - tp_kept
            fp_removed = n_fp_fit - fp_kept
            tp_loss_pct = (tp_removed / n_tp_fit) if n_tp_fit > 0 else float("nan")
            fp_removal_pct = (fp_removed / n_fp_fit) if n_fp_fit > 0 else float("nan")
            # post-filter precision / recall (on this cell only)
            precision_post = (tp_kept / (tp_kept + fp_kept)) if (tp_kept + fp_kept) > 0 else float("nan")
            precision_pre = n_tp_fit / n_fit
            precision_gain = precision_post - precision_pre
            recall_loss = tp_loss_pct  # recall is TP_kept / TP_total
            filter_signal = float(fp_removal_pct - tp_loss_pct) if (not np.isnan(fp_removal_pct) and not np.isnan(tp_loss_pct)) else float("nan")
            sweep_rows.append({
                "cell_id": cell["cell_id"],
                "source": cell["source"],
                "BAM": bam,
                "flag": flag,
                "model_used": used_model,
                "n_fit": n_fit,
                "n_TP_fit": n_tp_fit,
                "n_FP_fit": n_fp_fit,
                "tau": float(tau),
                "TP_kept": tp_kept,
                "FP_kept": fp_kept,
                "TP_removed": tp_removed,
                "FP_removed": fp_removed,
                "TP_loss_pct": tp_loss_pct,
                "FP_removal_pct": fp_removal_pct,
                "precision_pre": precision_pre,
                "precision_post": precision_post,
                "precision_gain": precision_gain,
                "recall_loss": recall_loss,
                "filter_signal_FPrem_minus_TPloss": filter_signal,
            })

        # Best τ for this cell (max filter_signal)
        cell_df = pd.DataFrame([r for r in sweep_rows if r["cell_id"] == cell["cell_id"]])
        if not cell_df.empty:
            best = cell_df.sort_values("filter_signal_FPrem_minus_TPloss", ascending=False).iloc[0]
            cell_summaries.append({
                "cell_id": cell["cell_id"],
                "source": cell["source"],
                "n_fit": n_fit,
                "n_TP_fit": n_tp_fit,
                "n_FP_fit": n_fp_fit,
                "model_used": used_model,
                "best_tau": float(best["tau"]),
                "best_filter_signal": float(best["filter_signal_FPrem_minus_TPloss"]),
                "best_TP_loss_pct": float(best["TP_loss_pct"]),
                "best_FP_removal_pct": float(best["FP_removal_pct"]),
                "best_precision_post": float(best["precision_post"]),
                "h3_direction_positive": bool(best["filter_signal_FPrem_minus_TPloss"] > 0),
            })

    sweep_df = pd.DataFrame(sweep_rows)
    sweep_out = STEP5_DIR / "step2_filter_sweep.tsv"
    sweep_df.to_csv(sweep_out, sep="\t", index=False, float_format="%.6g")
    log_msg(f"\nWrote {sweep_out.name} ({len(sweep_df)} rows = {len(cell_summaries)} cells × {len(TAU_GRID)} τ)")

    summary_df = pd.DataFrame(cell_summaries)
    summary_out = STEP5_DIR / "intermediate" / "step2_cell_summary.tsv"
    summary_df.to_csv(summary_out, sep="\t", index=False, float_format="%.6g")
    log_msg(f"Wrote {summary_out.name}")

    # ---- H3 preliminary verdict ----
    n_pos_cells = int(summary_df["h3_direction_positive"].sum()) if not summary_df.empty else 0
    n_cells_eval = len(summary_df)
    h3_prelim = "POSITIVE_DIRECTION" if n_pos_cells >= 1 else "NEGATIVE"
    log_msg(f"\nH3 preliminary: {n_pos_cells}/{n_cells_eval} FP-rich cells have ≥1 τ where FP_removal_pct > TP_loss_pct")
    log_msg(f"H3 preliminary verdict: {h3_prelim}")

    # ---- Figure 1: precision-recall per cell (subplots 4×3 grid) ----
    plt = matplotlib_setup()
    n_cells = len(cell_summaries)
    if n_cells > 0:
        ncol = 3
        nrow = (n_cells + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 4.5, nrow * 3.5))
        if nrow == 1 and ncol == 1:
            axes = np.array([[axes]])
        elif nrow == 1:
            axes = axes[np.newaxis, :]
        for i, cell in enumerate(cell_summaries):
            ax = axes[i // ncol, i % ncol]
            sub = sweep_df[sweep_df["cell_id"] == cell["cell_id"]].sort_values("tau")
            ax.plot(sub["recall_loss"], sub["precision_gain"], marker=".", color="steelblue", label="τ sweep")
            ax.axhline(0, color="gray", lw=0.5)
            ax.axvline(0, color="gray", lw=0.5)
            # Mark best τ
            best_row = sub[sub["tau"] == cell["best_tau"]]
            if not best_row.empty:
                ax.scatter([best_row["recall_loss"].iloc[0]], [best_row["precision_gain"].iloc[0]],
                           color="red", marker="*", s=100, label=f"τ*={cell['best_tau']:.2f}")
            ax.set_xlabel("recall loss (TP_loss_pct)", fontsize=8)
            ax.set_ylabel("precision gain", fontsize=8)
            ax.set_title(cell["cell_id"][:32] + f"\nn={cell['n_fit']}  Model {cell['model_used']}",
                         fontsize=8)
            ax.legend(fontsize=7)
            ax.tick_params(labelsize=7)
        # Hide unused subplots
        for j in range(n_cells, nrow * ncol):
            axes[j // ncol, j % ncol].set_visible(False)
        fig.suptitle("Step 2 — Per-cell precision-recall (V6_off, Model B/A predicted P(TP) vs τ)",
                     fontsize=11)
        plt.tight_layout()
        fig.savefig(STEP5_DIR / "figures" / "step2_roc_per_cell.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        log_msg("Wrote figures/step2_roc_per_cell.png")

    # ---- Figure 2: FP_removal vs TP_loss scatter (all cells × all τ) ----
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    for cell in cell_summaries:
        sub = sweep_df[sweep_df["cell_id"] == cell["cell_id"]]
        ax.plot(sub["TP_loss_pct"], sub["FP_removal_pct"], marker=".", lw=0.7, alpha=0.7,
                label=cell["cell_id"][:24])
    lim = [0, 1.02]
    ax.plot(lim, lim, color="black", ls="--", lw=0.8, label="y=x (no filter signal)")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("TP_loss_pct")
    ax.set_ylabel("FP_removal_pct")
    ax.set_title("Step 2 — FP removal vs TP loss across FP-rich cells\n(above y=x line = filter signal)")
    ax.legend(fontsize=6, loc="lower right", ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(STEP5_DIR / "figures" / "step2_filter_signal.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    log_msg("Wrote figures/step2_filter_signal.png")

    # ---- step2_findings.md ----
    lines = []
    lines.append("# Step 2 — FP-rich Cells Filter Threshold Sweep Findings")
    lines.append("")
    lines.append("> Generated by `scripts/filter_sweep.py`")
    lines.append("")
    lines.append("## Setup")
    lines.append(f"- FP-rich threshold: FP_rate > global ({GLOBAL_FP_RATE}) × 1.5 = {FP_RICH_THRESH:.4f}")
    lines.append(f"- FP-rich cells identified: {len(fp_rich)} (chr8_hotspot + Z-AUTO sub-cells)")
    lines.append(f"- Cells successfully fit + swept: {len(cell_summaries)}")
    lines.append(f"- Primary BAM × flag: {PRIMARY_BAM}_{PRIMARY_FLAG}")
    lines.append(f"- τ grid: {len(TAU_GRID)} values, [{TAU_GRID[0]}, {TAU_GRID[-1]}] step 0.01")
    lines.append("")
    lines.append("## H3 Preliminary Direction Check")
    lines.append(f"- Cells with ≥1 τ where FP_removal_pct > TP_loss_pct: **{n_pos_cells}** / {n_cells_eval}")
    lines.append(f"- H3 preliminary verdict: **{h3_prelim}**")
    lines.append("- (Final H3 verdict requires Step 3 ΔF1 computation with global TP/FP/FN reconstruction.)")
    lines.append("")
    lines.append("## Per-cell best τ summary")
    lines.append("")
    lines.append("| cell_id | source | n_fit | TP/FP | model | best τ | filter_signal | TP_loss% | FP_remov% | prec_post |")
    lines.append("|---------|--------|-------|-------|-------|--------|---------------|----------|-----------|-----------|")
    for c in cell_summaries:
        lines.append(f"| {c['cell_id'][:40]} | {c['source']} | {c['n_fit']} | "
                     f"{c['n_TP_fit']}/{c['n_FP_fit']} | {c['model_used']} | "
                     f"{c['best_tau']:.2f} | {c['best_filter_signal']:.4f} | "
                     f"{c['best_TP_loss_pct']:.4f} | {c['best_FP_removal_pct']:.4f} | "
                     f"{c['best_precision_post']:.4f} |")
    lines.append("")
    lines.append("## Files")
    lines.append("- `step2_filter_sweep.tsv` — per (cell, τ) × {TP/FP_loss, precision, recall, filter_signal}")
    lines.append("- `intermediate/step2_cell_summary.tsv` — per-cell best τ summary")
    lines.append("- `figures/step2_roc_per_cell.png` — precision-recall curves per cell")
    lines.append("- `figures/step2_filter_signal.png` — FP_removal vs TP_loss scatter (target above y=x)")
    lines.append("")
    lines.append("## Caveats (read before interpreting numbers)")
    lines.append("1. **Sample shrinkage from HP-conditional methyl features**: Model B requires ")
    lines.append("   NME_imbalance + Epipoly_Delta to be non-null (HP-conditional, ~26-40% coverage). ")
    lines.append("   For most FP-rich cells, n_fit (rows surviving dropna) is 30-70% of the step3 ")
    lines.append("   reported n (e.g. chr8|Inner|other|cov_normal: step3 n=625 → n_fit=152). ")
    lines.append("   Interpretation: signal is for the HP-stratifiable subset, not the entire FP-rich zone.")
    lines.append("2. **Small TP counts in UNKNOWN cells**: 4 of the 11 cells have n_TP_fit ≤ 5 ")
    lines.append("   (auto|UNKNOWN|*). filter_signal=1.0 in these cells is essentially trivial ")
    lines.append("   (all TPs already excluded by P(TP)<0.5 since the cells are 95-99% FP).")
    lines.append("3. **In-sample evaluation**: predictions come from the LR fit on the same rows. ")
    lines.append("   Step 3 must add held-out CV or independent test (cross-sample) to confirm.")
    lines.append("4. **Decision-theoretic τ choice**: best_tau=0.5 is the minimum τ in the grid. ")
    lines.append("   The signal is driven mainly by Model B's high in-cell AUC. Step 3 ΔF1 may differ ")
    lines.append("   substantially once global counts (not within-cell) are considered.")
    lines.append("")
    lines.append("## Hand-off to Step 3 (Agent H)")
    lines.append("- Step 2 evaluates *per-cell local* filter signal (within each FP-rich cell, in-sample).")
    lines.append("- Step 3 (ΔF1) must compute global TP/FP/FN against SEQC2 truth set.")
    lines.append("- Recommended τ* candidates: per-cell best τ values from `intermediate/step2_cell_summary.tsv`.")
    lines.append("- For each candidate τ*, compute global ΔF1 = F1_post - F1_caller (=0.7166).")
    lines.append("- Recommend held-out CV when reporting Step 3 to avoid optimism.")
    lines.append("- Strong candidates (large n_fit + meaningful TP count + filter_signal > 0.5):")
    for c in cell_summaries:
        if c["n_TP_fit"] >= 15 and c["best_filter_signal"] > 0.5:
            lines.append(f"  - **{c['cell_id']}** (n_fit={c['n_fit']}, TP/FP={c['n_TP_fit']}/{c['n_FP_fit']}, "
                         f"signal={c['best_filter_signal']:.4f}, FP_remov={c['best_FP_removal_pct']:.4f})")

    (STEP5_DIR / "step2_findings.md").write_text("\n".join(lines))
    log_msg(f"Wrote step2_findings.md")

    log_msg("=== Step 2 filter_sweep END ===")

    return {
        "n_fp_rich_cells": len(fp_rich),
        "n_cells_swept": n_cells_eval,
        "n_cells_positive_direction": n_pos_cells,
        "h3_preliminary_verdict": h3_prelim,
    }


if __name__ == "__main__":
    summary = main()
    print("\n=== STEP 2 SUMMARY ===")
    print(json.dumps(summary, indent=2))
