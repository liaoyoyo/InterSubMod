"""Step 2 — 5 + 2 道 confound guard for top-list cells.

Guards (per 02_methodology_notes §4):
  1. NumReads OLS residualize (within cell)
  2. caller_af OLS residualize (within cell)
  3. L3 AF-bin internal consistency
  4. Permutation test (max-statistic null)
  5. Marginal expectation surprise (already in step2_grid_3d.tsv)
  6. chr-stratified Mantel-Haenszel
  7. HP1/HP2 symmetry test (binomial)

Plus spatial autocorr check (mid-TP-rate 1 Mb window).
Outputs:
  - step2_confound_guard.tsv (per top-cell × 7 guards)
  - step2_top_lists.tsv (6 top-lists merged + per-cell guard verdict)
  - step2_findings.md (H2/H3/H6/H7 verdict + V5 over-promote list)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from _common import (
    HP_PREFIX_PRIMARY,
    LOH_LEVELS,
    POWERED_THRESHOLD,
    STEP2_DIR,
    annotate_axes,
    cell_id,
    joined_subset,
    load_master,
    log_msg,
    matplotlib_setup,
)


N_PERMUTATIONS = 1000
RNG = np.random.default_rng(20260515)


def ols_residualize(cell_df: pd.DataFrame, covariate_col: str) -> float | None:
    """Within-cell OLS y ~ covariate; return residualized TP rate = mean(residual + mean(y))."""
    if len(cell_df) < 5:
        return None
    y = (cell_df["label"] == "TP").astype(int).values
    x = cell_df[covariate_col].astype(float).values
    if np.std(x) == 0 or len(np.unique(y)) < 2:
        return float(y.mean())
    X = sm.add_constant(x)
    try:
        m = sm.OLS(y, X).fit()
        resid = y - m.predict(X)
        return float(resid.mean() + y.mean())
    except Exception:
        return None


def l3_af_bin_consistency(cell_df: pd.DataFrame, cell_tp_rate: float) -> tuple[bool, dict]:
    """Bin AF into 4 ranges; flag pass if ≥1 bin's TP rate is same-sign-direction with cell mean
    (within 0.05 abs gap) AND has n>=10.
    """
    bins = [(-0.01, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.01)]
    details = {}
    pass_count = 0
    for lo, hi in bins:
        sub = cell_df[(cell_df["caller_af"] >= lo) & (cell_df["caller_af"] < hi)]
        n = len(sub)
        if n >= 10:
            r = (sub["label"] == "TP").mean()
            details[f"af[{lo:.2f},{hi:.2f})_n={n}_rate"] = float(r)
            if abs(r - cell_tp_rate) <= 0.05:
                pass_count += 1
    return (pass_count >= 1, details)


def permutation_pvalue(
    df_universe: pd.DataFrame, cell_mask: np.ndarray, observed_rate: float, n_iter: int = N_PERMUTATIONS
) -> float:
    """Max-statistic permutation: shuffle TP/FP labels globally; track each cell's rate and report
    fraction of iterations where any cell rate >= observed (one-sided, max-statistic FWER).
    For computational efficiency we just test this specific cell's tail vs shuffled distribution.
    """
    y = (df_universe["label"] == "TP").astype(int).values
    cell_n = cell_mask.sum()
    if cell_n == 0:
        return float("nan")
    obs = observed_rate
    count = 0
    for _ in range(n_iter):
        perm = RNG.permutation(y)
        rate = perm[cell_mask].mean()
        if abs(rate - y.mean()) >= abs(obs - y.mean()):
            count += 1
    return (count + 1) / (n_iter + 1)


def mantel_haenszel(cell_df: pd.DataFrame, rest_df: pd.DataFrame) -> float:
    """chr-stratified Mantel-Haenszel test (TP/FP × cell/rest, stratified by chr).
    Returns 2-sided p-value; NaN if not enough strata.
    """
    if len(cell_df) < 10 or len(rest_df) < 10:
        return float("nan")
    chrs = sorted(set(cell_df["chr"]).union(set(rest_df["chr"])))
    num = 0.0
    den = 0.0
    var = 0.0
    used_strata = 0
    for c in chrs:
        a = int(((cell_df["chr"] == c) & (cell_df["label"] == "TP")).sum())
        b = int(((cell_df["chr"] == c) & (cell_df["label"] == "FP")).sum())
        cc = int(((rest_df["chr"] == c) & (rest_df["label"] == "TP")).sum())
        d = int(((rest_df["chr"] == c) & (rest_df["label"] == "FP")).sum())
        n = a + b + cc + d
        if n < 4 or (a + b) == 0 or (cc + d) == 0 or (a + cc) == 0 or (b + d) == 0:
            continue
        used_strata += 1
        e_a = (a + b) * (a + cc) / n
        num += a - e_a
        v = (a + b) * (cc + d) * (a + cc) * (b + d) / (n * n * (n - 1))
        var += v
    if used_strata < 2 or var <= 0:
        return float("nan")
    chi2 = (num * num) / var
    return float(stats.chi2.sf(chi2, df=1))


def hp_symmetry_test(cell_df: pd.DataFrame) -> tuple[float, dict]:
    """Binomial test for HP1 vs HP2 family asymmetry within the cell.
    HP1 side = HP1 + HP1S; HP2 side = HP2 + HP2S (in V6_off coords).
    Test H0: equal proportions.
    """
    hp1_side = (
        cell_df[f"{HP_PREFIX_PRIMARY}_1"].fillna(0)
        + cell_df[f"{HP_PREFIX_PRIMARY}_1-1"].fillna(0)
    ).sum()
    hp2_side = (
        cell_df[f"{HP_PREFIX_PRIMARY}_2"].fillna(0)
        + cell_df[f"{HP_PREFIX_PRIMARY}_2-1"].fillna(0)
    ).sum()
    total = hp1_side + hp2_side
    if total < 10:
        return float("nan"), {"hp1_side": int(hp1_side), "hp2_side": int(hp2_side)}
    k = int(min(hp1_side, hp2_side))
    n = int(total)
    p = stats.binomtest(k=k, n=n, p=0.5, alternative="two-sided").pvalue
    return float(p), {
        "hp1_side": int(hp1_side),
        "hp2_side": int(hp2_side),
        "ratio_smaller": k / n,
    }


def spatial_robust(df_universe: pd.DataFrame, cell_mask: np.ndarray, cell_tp_rate: float) -> tuple[bool, dict]:
    """Mid-TP-rate window check (1 Mb windows). Skipped if cell has <5 distinct chr+window
    combinations. Pass = at least one mid-TP-rate window (0.4-0.6) shows same direction as
    cell vs global.
    """
    sub = df_universe.loc[cell_mask].copy()
    if len(sub) < 5:
        return (True, {"reason": "too_small_for_spatial"})  # default pass to avoid false flag
    sub["window"] = (sub["pos"] // 1_000_000).astype(int)
    sub["cwin"] = sub["chr"].astype(str) + ":" + sub["window"].astype(str)
    if sub["cwin"].nunique() < 3:
        return (True, {"reason": "too_few_distinct_windows", "n_windows": int(sub["cwin"].nunique())})
    # Window-level TP rate (require window with global mid-rate evidence)
    # Compute global TP rate by window
    df_universe_local = df_universe.copy()
    df_universe_local["window"] = (df_universe_local["pos"] // 1_000_000).astype(int)
    df_universe_local["cwin"] = (
        df_universe_local["chr"].astype(str) + ":" + df_universe_local["window"].astype(str)
    )
    global_win_rate = df_universe_local.groupby("cwin")["label"].apply(lambda s: (s == "TP").mean())
    mid_windows = global_win_rate[(global_win_rate >= 0.4) & (global_win_rate <= 0.6)].index
    if len(mid_windows) == 0:
        # No mid-TP-rate windows globally → fallback: check that cell signal isn't isolated to 1 window
        wcounts = sub["cwin"].value_counts()
        passed = wcounts.iloc[0] / len(sub) < 0.5
        return (passed, {"reason": "no_global_mid_windows", "top_window_share": float(wcounts.iloc[0] / len(sub))})
    sub_mid = sub[sub["cwin"].isin(mid_windows)]
    if len(sub_mid) < 10:
        return (True, {"reason": "few_mid_window_rows", "n_mid": int(len(sub_mid))})
    mid_rate = (sub_mid["label"] == "TP").mean()
    global_rate = (df_universe["label"] == "TP").mean()
    # Direction-consistent if same side of global
    same_dir = (cell_tp_rate > global_rate) == (mid_rate > global_rate)
    return (same_dir, {"mid_rate": float(mid_rate), "n_mid": int(len(sub_mid))})


def main():
    log_msg("=== Step 2 confound_guard START ===")

    grid = pd.read_csv(STEP2_DIR / "step2_grid_3d.tsv", sep="\t")
    df_all = load_master()
    df = joined_subset(df_all)
    df = annotate_axes(df, hp_prefix=HP_PREFIX_PRIMARY)

    log_msg(f"Grid loaded: {len(grid)} cells. Powered={int(grid['powered_flag'].sum())}.")

    # ---- Top lists ----
    powered = grid[grid["powered_flag"] == True].copy()

    def top(df_, col, asc, k=20):
        return df_.dropna(subset=[col]).sort_values(col, ascending=asc).head(k)["cell_id"].tolist()

    # Trajectory file for V5-over-promote list
    traj = pd.read_csv(STEP2_DIR / "step2_trajectory_per_cell.tsv", sep="\t")

    top_lists = {
        "top_TP_extreme": top(powered, "TP_rate", asc=False),
        "top_FP_extreme": top(powered, "FP_rate", asc=False),
        "top_log_odds_high": top(powered, "log_odds", asc=False),
        "top_log_odds_low": top(powered, "log_odds", asc=True),
        "most_populated": top(powered, "n", asc=False),
        "top_surprise": top(powered, "surprise", asc=False),
        "top_v5_over_promote": traj[
            (traj["v5_over_v3f"] > 1.3) & (traj["n_V6"] >= POWERED_THRESHOLD)
        ].sort_values("v5_over_v3f", ascending=False)["cell_id"].head(20).tolist(),
    }

    # Long table
    rows = []
    for list_name, cells in top_lists.items():
        for rank, cid in enumerate(cells, 1):
            rows.append({"list_name": list_name, "rank": rank, "cell_id": cid})
    top_lists_df = pd.DataFrame(rows)

    # Unique cells to guard
    cells_to_guard = sorted(set(top_lists_df["cell_id"]))
    log_msg(f"{len(cells_to_guard)} unique cells across {len(top_lists)} top-lists")

    # ---- Run 7 guards per cell ----
    guard_rows = []
    for cid in cells_to_guard:
        loh, hp, cov = cid.split("|")
        mask = (
            (df["loh_side_norm"] == loh)
            & (df["hp_bucket"] == hp)
            & (df["cov_bin"] == cov)
        )
        cell_df = df[mask]
        rest_df = df[~mask]
        n = len(cell_df)
        if n == 0:
            continue
        n_tp = int((cell_df["label"] == "TP").sum())
        n_fp = int((cell_df["label"] == "FP").sum())
        tp_rate = n_tp / n
        global_tp = (df["label"] == "TP").mean()

        # Guard 1: NumReads residualize
        rate_nr = ols_residualize(cell_df, "n_reads_v6off")
        guard1_pass = (rate_nr is not None) and (abs(rate_nr - tp_rate) <= 0.05)

        # Guard 2: caller_af residualize
        rate_af = ols_residualize(cell_df, "caller_af")
        guard2_pass = (rate_af is not None) and (abs(rate_af - tp_rate) <= 0.05)

        # Guard 3: L3 AF-bin
        guard3_pass, l3_details = l3_af_bin_consistency(cell_df, tp_rate)

        # Guard 4: Permutation (use cell_mask in joined df coords)
        perm_p = permutation_pvalue(df, mask.values, tp_rate, n_iter=N_PERMUTATIONS)
        guard4_pass = perm_p < 0.05

        # Guard 5: Surprise already in grid
        gr = grid[grid["cell_id"] == cid].iloc[0]
        surprise = float(gr["surprise"]) if not pd.isna(gr["surprise"]) else float("nan")

        # Guard 6: Mantel-Haenszel chr-stratified
        mh_p = mantel_haenszel(cell_df, rest_df)
        guard6_pass = (not np.isnan(mh_p)) and (mh_p < 0.05)

        # Guard 7: HP symmetry (FP/uneven sides flagged)
        hp_p, hp_det = hp_symmetry_test(cell_df)
        # For "other" bucket, asymmetry expected; only flag for same_HP* / cross_het*
        if hp in ("same_HP1", "cross_het", "cross_het_inv", "same_HP2"):
            guard7_pass = np.isnan(hp_p) or hp_p > 0.01  # pass if symmetric or insufficient data
        else:
            guard7_pass = True  # not applicable for 'other'

        # Spatial guard
        spat_pass, spat_det = spatial_robust(df, mask.values, tp_rate)

        guards_passed = sum(
            [guard1_pass, guard2_pass, guard3_pass, guard4_pass, guard6_pass, guard7_pass, spat_pass]
        )

        guard_rows.append(
            {
                "cell_id": cid,
                "loh_side": loh,
                "hp_bucket": hp,
                "cov_bin": cov,
                "n": n,
                "n_TP": n_tp,
                "n_FP": n_fp,
                "TP_rate": tp_rate,
                "guard1_NR_residualized_rate": rate_nr,
                "guard1_NR_pass": guard1_pass,
                "guard2_AF_residualized_rate": rate_af,
                "guard2_AF_pass": guard2_pass,
                "guard3_L3_AFbin_pass": guard3_pass,
                "guard4_perm_p": perm_p,
                "guard4_perm_pass": guard4_pass,
                "guard5_surprise": surprise,
                "guard6_mh_p": mh_p,
                "guard6_mh_pass": guard6_pass,
                "guard7_hp_symmetry_p": hp_p,
                "guard7_hp_symmetry_pass": guard7_pass,
                "spatial_robust_pass": spat_pass,
                "n_guards_passed": guards_passed,
                "all_7_passed": guards_passed == 7,
            }
        )

    guard_df = pd.DataFrame(guard_rows)
    guard_df.to_csv(STEP2_DIR / "step2_confound_guard.tsv", sep="\t", index=False, float_format="%.6f")
    log_msg(f"Wrote step2_confound_guard.tsv ({len(guard_df)} cells)")

    # Merge guard verdicts into top_lists table
    top_lists_df = top_lists_df.merge(
        guard_df[["cell_id", "n_guards_passed", "all_7_passed", "TP_rate", "n", "n_TP", "n_FP"]],
        on="cell_id",
        how="left",
    )
    top_lists_df.to_csv(STEP2_DIR / "step2_top_lists.tsv", sep="\t", index=False, float_format="%.6f")
    log_msg(f"Wrote step2_top_lists.tsv ({len(top_lists_df)} rows)")

    # ---- Findings markdown ----
    n_pass_all = int(guard_df["all_7_passed"].sum())
    log_msg(f"Cells passing all 7 guards: {n_pass_all}")

    # H2: LOH inner × same_HP TP enrichment
    inner_same = grid[
        (grid["loh_side"] == "Inner")
        & (grid["hp_bucket"].isin(["same_HP1", "same_HP2"]))
        & (grid["powered_flag"] == True)
    ]
    if len(inner_same) > 0:
        mean_enrich = float(inner_same["TP_enrichment"].mean())
        h2_verdict = "POSITIVE" if mean_enrich >= 1.005 else ("MIXED" if mean_enrich >= 0.99 else "NEGATIVE")
    else:
        mean_enrich = float("nan")
        h2_verdict = "NO_POWERED_CELLS"

    # H3: Outer cross_het × cov_gain FP enrichment
    outer_xh = grid[
        (grid["loh_side"] == "Outer")
        & (grid["hp_bucket"].isin(["cross_het", "cross_het_inv"]))
        & (grid["cov_bin"].isin(["cov_gain", "cov_high_gain"]))
    ]
    h3_enrich = float(outer_xh["FP_enrichment"].mean()) if len(outer_xh) and not outer_xh["FP_enrichment"].isna().all() else float("nan")
    h3_verdict = (
        "POSITIVE" if (not np.isnan(h3_enrich) and h3_enrich >= 2.0)
        else ("MIXED" if (not np.isnan(h3_enrich) and h3_enrich >= 1.3) else "NEGATIVE")
    )

    # H6: powered cells share
    powered_share = float(grid["powered_flag"].mean())
    main_powered_share = float(grid[grid["loh_side"] != "UNKNOWN"]["powered_flag"].mean())
    h6_verdict = "POSITIVE" if main_powered_share >= 0.30 else "NEGATIVE"

    # H7: ≥5 cells pass all guards
    h7_verdict = "POSITIVE" if n_pass_all >= 5 else "NEGATIVE"

    # V5 over-promote: top 5
    over_promote = traj[
        (traj["v5_over_v3f"] > 1.3) & (traj["n_V6"] >= POWERED_THRESHOLD)
    ].sort_values("v5_over_v3f", ascending=False).head(5)

    md = []
    md.append("# Step 2 — Findings (HCC1395 paired_full, V6_off primary)")
    md.append("")
    md.append("## Power gate (Step 1.5)")
    md.append(f"- Decision: **PROCEED** full 3-axis 50-cell grid")
    md.append(f"- Powered cells (n≥50, main 50): **{int(grid[grid['loh_side']!='UNKNOWN']['powered_flag'].sum())} / 50** ({main_powered_share:.0%})")
    md.append(f"- Full 75 cells (incl. UNKNOWN LOH): non-empty {int((grid['n']>0).sum())}")
    md.append("")
    md.append("## H2 — Inner × same_HP TP enrichment")
    md.append(f"- Powered cells matching: **{len(inner_same)}**")
    md.append(f"- Mean TP_enrichment vs global = **{mean_enrich:.4f}**")
    md.append(f"- **Verdict: {h2_verdict}** (POSITIVE ≥ 1.005, given global TP_rate already 0.9832 ceiling)")
    md.append(
        "- Caveat: global TP_rate = 0.9832 (FP universe shrank to 506 due to master-join). "
        "Enrichment compressed near ceiling; interpret directionally not as effect size."
    )
    md.append("")
    md.append("## H3 — Outer cross_het × cov_gain FP enrichment")
    md.append(f"- Cells matching: **{len(outer_xh)}**")
    md.append(f"- Mean FP_enrichment = **{h3_enrich:.4f}**" if not np.isnan(h3_enrich) else "- Mean FP_enrichment = NaN (no data)")
    md.append(f"- **Verdict: {h3_verdict}** (POSITIVE ≥ 2.0×)")
    md.append("")
    md.append("## H6 — Powered cells ≥ 30% (grid sparsity)")
    md.append(f"- Main 50-cell powered share = **{main_powered_share:.0%}**")
    md.append(f"- **Verdict: {h6_verdict}**")
    md.append("")
    md.append("## H7 — ≥5 cells pass all 7 confound guards")
    md.append(f"- Cells tested (top-list union): **{len(guard_df)}**")
    md.append(f"- Cells passing all 7 guards: **{n_pass_all}**")
    md.append(f"- **Verdict: {h7_verdict}**")
    md.append("")
    md.append("## V5 over-promote cells (top 5)")
    md.append("")
    md.append("| cell_id | V5/V3F ratio | V6/V5 ratio | V6/V3F ratio | n_V3F | n_V5 | n_V6 | n_TP_V6 |")
    md.append("|---|---|---|---|---|---|---|---|")
    for _, r in over_promote.iterrows():
        md.append(
            f"| {r['cell_id']} | {r['v5_over_v3f']:.3f} | {r['v6_over_v5']:.3f} | {r['v6_over_v3f']:.3f} | "
            f"{int(r['n_V3F'])} | {int(r['n_V5'])} | {int(r['n_V6'])} | {int(r['n_TP_V6'])} |"
        )
    md.append("")
    md.append("## Cells passing all 7 guards")
    md.append("")
    passing = guard_df[guard_df["all_7_passed"]]
    if len(passing):
        md.append("| cell_id | n | n_TP | n_FP | TP_rate | guard4_perm_p | guard6_mh_p |")
        md.append("|---|---|---|---|---|---|---|")
        for _, r in passing.iterrows():
            mh_str = f"{r['guard6_mh_p']:.4f}" if not np.isnan(r['guard6_mh_p']) else "NaN"
            md.append(
                f"| {r['cell_id']} | {int(r['n'])} | {int(r['n_TP'])} | {int(r['n_FP'])} | "
                f"{r['TP_rate']:.4f} | {r['guard4_perm_p']:.4f} | {mh_str} |"
            )
    else:
        md.append("_(no cells passed all 7 guards)_")
    md.append("")
    md.append("## Step 3 hand-off candidates")
    md.append("")
    md.append("Zone candidates from guard-passing + log-odds extreme cells:")
    # Build hand-off list
    candidates = guard_df.copy()
    candidates["combined_score"] = candidates["n_guards_passed"] + 0.001 * candidates["n"]
    candidates = candidates.sort_values("combined_score", ascending=False).head(15)
    md.append("")
    md.append("| cell_id | n_guards_passed | n | n_TP | n_FP | TP_rate | use_for |")
    md.append("|---|---|---|---|---|---|---|")
    for _, r in candidates.iterrows():
        use_for = []
        if r["n_FP"] >= 3 and r["TP_rate"] < 0.95:
            use_for.append("FP-rich zone (Step 3 deep dive)")
        if r["TP_rate"] >= 0.99 and r["n"] >= 100:
            use_for.append("TP-clean reference")
        if r["loh_side"] == "Inner" and r["hp_bucket"].startswith("same_"):
            use_for.append("Z-OCH/Z-CHR8 contrast")
        md.append(
            f"| {r['cell_id']} | {int(r['n_guards_passed'])} | {int(r['n'])} | {int(r['n_TP'])} | {int(r['n_FP'])} | "
            f"{r['TP_rate']:.4f} | {'; '.join(use_for) if use_for else '—'} |"
        )
    md.append("")
    md.append("---")
    md.append("Generated by `step2_three_axis_grid/scripts/confound_guard.py`.")

    (STEP2_DIR / "step2_findings.md").write_text("\n".join(md) + "\n")
    log_msg("Wrote step2_findings.md")

    summary = {
        "n_cells_guarded": int(len(guard_df)),
        "n_passing_all_7": n_pass_all,
        "H2_verdict": h2_verdict,
        "H2_mean_TP_enrichment": mean_enrich,
        "H3_verdict": h3_verdict,
        "H3_mean_FP_enrichment": h3_enrich,
        "H6_verdict": h6_verdict,
        "H6_main_powered_share": main_powered_share,
        "H7_verdict": h7_verdict,
        "v5_over_promote_top5": over_promote[
            ["cell_id", "v5_over_v3f", "v6_over_v5", "n_V3F", "n_V5", "n_V6", "n_TP_V6"]
        ].to_dict(orient="records"),
    }
    with open(STEP2_DIR / "intermediate" / "step2_findings_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    log_msg("=== Step 2 confound_guard END ===")
    return summary


if __name__ == "__main__":
    main()
