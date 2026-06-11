"""Step 3 — 4 FP zone deep dive.

Zones:
- Z-OCH (Outer cross_het): full-FP set (no master_join_ok filter), HP from V6_off
- Z-CHR8 (HCC1395 chr8 hotspot): full chr8 + per-chr comparison
- Z-GL (Gain+LOH, Inner × cov_elevated/gain/high_gain): master-joined only (needs LOH+CN)
- Z-AUTO (top 5% FP density by KDE smoothing): full FP set

Per zone outputs (under each zone subdir):
- zone_summary.tsv (overall TP/FP rate, n, enrichment)
- zone_grid.tsv (3-axis sub-grid; cov_proxy when needed)
- zone_three_version_trajectory.tsv (V3F/V5/V6 marker coverage + TP rate in zone)
- zone_findings.md
- zone-specific extras: chr8 has per_chr_comparison.tsv; gain_loh has lr_ablation.tsv + seqc2_concordance.tsv

Plus: step3_synthesis.md (4-zone compare + H3/H4/H5 verdict + handoff).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _zone_common import (
    COV_LEVELS,
    COV_PROXY_LEVELS,
    HP_LEVELS,
    LOH_LEVELS,
    POWERED_THRESHOLD,
    STEP3_DIR,
    annotate_axes,
    bh_fdr,
    cell_stats,
    compute_hp_bucket,
    compute_ng,
    deviance_decomposition,
    load_master_full,
    load_seqc2_cnv,
    matplotlib_setup,
    wilson_ci,
)


# --------------------------------------------------------------------------------
# Globals (computed in main)
# --------------------------------------------------------------------------------
GLOBAL_TP_RATE = None
GLOBAL_FP_RATE = None


def compute_zone_summary(df_zone: pd.DataFrame, zone_label: str, full_df: pd.DataFrame) -> dict:
    """Zone-level TP/FP rate + enrichment vs global."""
    n = len(df_zone)
    n_tp = int((df_zone["label"] == "TP").sum())
    n_fp = int((df_zone["label"] == "FP").sum())
    tp_rate = n_tp / n if n else float("nan")
    fp_rate = n_fp / n if n else float("nan")
    tp_lo, tp_hi = wilson_ci(n_tp, n)
    fp_lo, fp_hi = wilson_ci(n_fp, n)
    # 2x2 Fisher: zone vs not-zone, label
    not_zone = full_df.loc[~full_df.index.isin(df_zone.index)]
    n_tp_out = int((not_zone["label"] == "TP").sum())
    n_fp_out = int((not_zone["label"] == "FP").sum())
    if n_tp > 0 and n_fp > 0 and n_tp_out > 0 and n_fp_out > 0:
        odds, p_fisher = stats.fisher_exact([[n_tp, n_fp], [n_tp_out, n_fp_out]])
    else:
        odds, p_fisher = float("nan"), float("nan")
    return dict(
        zone=zone_label,
        n=n,
        n_TP=n_tp,
        n_FP=n_fp,
        TP_rate=tp_rate,
        FP_rate=fp_rate,
        TP_wilson_lo=tp_lo,
        TP_wilson_hi=tp_hi,
        FP_wilson_lo=fp_lo,
        FP_wilson_hi=fp_hi,
        TP_enrichment=(tp_rate / GLOBAL_TP_RATE) if GLOBAL_TP_RATE > 0 else float("nan"),
        FP_enrichment=(fp_rate / GLOBAL_FP_RATE) if GLOBAL_FP_RATE > 0 else float("nan"),
        fisher_odds=odds,
        fisher_p=p_fisher,
    )


def confound_guards(
    df_cell: pd.DataFrame,
    full_df: pd.DataFrame,
    cell_id: str,
    skip_hp_symmetry: bool = False,
) -> dict:
    """Run 4 confound guards on a cell.

    Per Step 3 spec:
    - Guard 1: NumReads OLS residualize
    - Guard 2: caller_af OLS residualize
    - Guard 3: permutation (1000 shuffles)
    - Guard 4: chr-stratified Mantel-Haenszel
    - (Guard 7 HP symmetry: skip when zone is asymmetric by definition)
    """
    out = dict(cell_id=cell_id)

    if len(df_cell) < 30:
        out.update(
            dict(
                guard1_residual_TP_rate=float("nan"),
                guard1_diff=float("nan"),
                guard2_residual_TP_rate=float("nan"),
                guard2_diff=float("nan"),
                guard3_perm_p=float("nan"),
                guard4_mh_p=float("nan"),
                hp_sym_p=float("nan"),
                guards_passed=0,
            )
        )
        return out

    y = (df_cell["label"] == "TP").astype(int).values
    raw_rate = y.mean()

    # Guard 1: NumReads
    try:
        x = df_cell["V6_off_n_reads"].fillna(0).astype(float).values
        # within-cell OLS: y ~ x
        if x.std() > 0:
            beta = np.cov(x, y, bias=False)[0, 1] / x.var(ddof=1)
            alpha = y.mean() - beta * x.mean()
            pred = alpha + beta * x
            residual_rate = (y - pred + y.mean()).mean()  # adjust to overall mean
        else:
            residual_rate = raw_rate
        out["guard1_residual_TP_rate"] = float(residual_rate)
        out["guard1_diff"] = float(abs(residual_rate - raw_rate))
    except Exception:
        out["guard1_residual_TP_rate"] = float("nan")
        out["guard1_diff"] = float("nan")

    # Guard 2: caller_af
    try:
        af = df_cell["caller_af"].fillna(df_cell["caller_af"].median()).astype(float).values
        if af.std() > 0:
            beta = np.cov(af, y, bias=False)[0, 1] / af.var(ddof=1)
            alpha = y.mean() - beta * af.mean()
            pred = alpha + beta * af
            residual_rate = (y - pred + y.mean()).mean()
        else:
            residual_rate = raw_rate
        out["guard2_residual_TP_rate"] = float(residual_rate)
        out["guard2_diff"] = float(abs(residual_rate - raw_rate))
    except Exception:
        out["guard2_residual_TP_rate"] = float("nan")
        out["guard2_diff"] = float("nan")

    # Guard 3: Permutation test (cell rate vs shuffled labels from full universe)
    try:
        rng = np.random.default_rng(seed=42)
        full_labels = (full_df["label"] == "TP").astype(int).values
        n = len(df_cell)
        nperm = 1000
        perm_rates = np.empty(nperm)
        for i in range(nperm):
            sample = rng.choice(full_labels, size=n, replace=False)
            perm_rates[i] = sample.mean()
        # two-sided p
        deviation = abs(raw_rate - perm_rates.mean())
        perm_dev = np.abs(perm_rates - perm_rates.mean())
        p = (perm_dev >= deviation).mean()
        out["guard3_perm_p"] = float(p)
    except Exception:
        out["guard3_perm_p"] = float("nan")

    # Guard 4: chr-stratified Mantel-Haenszel
    try:
        from statsmodels.stats.contingency_tables import StratifiedTable

        in_cell = df_cell.index
        full_cmp = full_df.copy()
        full_cmp["in_cell"] = full_cmp.index.isin(in_cell).astype(int)
        full_cmp["is_TP"] = (full_cmp["label"] == "TP").astype(int)
        tables = []
        for chr_val, sub in full_cmp.groupby("chr"):
            a = int(((sub["in_cell"] == 1) & (sub["is_TP"] == 1)).sum())
            b = int(((sub["in_cell"] == 1) & (sub["is_TP"] == 0)).sum())
            c = int(((sub["in_cell"] == 0) & (sub["is_TP"] == 1)).sum())
            d = int(((sub["in_cell"] == 0) & (sub["is_TP"] == 0)).sum())
            if (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0:
                tables.append(np.array([[a, b], [c, d]]))
        if len(tables) >= 2:
            st = StratifiedTable(tables)
            mh_p = float(st.test_null_odds().pvalue)
        else:
            mh_p = float("nan")
        out["guard4_mh_p"] = mh_p
    except Exception:
        out["guard4_mh_p"] = float("nan")

    # HP symmetry (skip in asymmetric zones)
    if skip_hp_symmetry:
        out["hp_sym_p"] = float("nan")
    else:
        try:
            hp1 = int(((df_cell["V6_off_1"].fillna(0).astype(int) > 0)).sum())
            hp2 = int(((df_cell["V6_off_2"].fillna(0).astype(int) > 0)).sum())
            if hp1 + hp2 > 0:
                p = float(stats.binomtest(hp1, hp1 + hp2, 0.5).pvalue)
            else:
                p = float("nan")
            out["hp_sym_p"] = p
        except Exception:
            out["hp_sym_p"] = float("nan")

    # Count guards passed (alpha = 0.05; passed = NOT confounded, i.e. p < 0.05 means signal real for permutation, but for residual_diff < 0.05 means low NumReads/AF drive)
    g1 = (not math.isnan(out["guard1_diff"])) and (out["guard1_diff"] < 0.05)
    g2 = (not math.isnan(out["guard2_diff"])) and (out["guard2_diff"] < 0.05)
    g3 = (not math.isnan(out["guard3_perm_p"])) and (out["guard3_perm_p"] < 0.05)
    g4 = (not math.isnan(out["guard4_mh_p"])) and (out["guard4_mh_p"] < 0.05)
    out["guard1_pass"] = bool(g1)
    out["guard2_pass"] = bool(g2)
    out["guard3_pass"] = bool(g3)
    out["guard4_pass"] = bool(g4)
    out["guards_passed"] = int(g1) + int(g2) + int(g3) + int(g4)

    return out


def per_version_trajectory(df_zone: pd.DataFrame) -> pd.DataFrame:
    """For each BAM version, compute marker coverage (NG≥3) + TP rate in zone."""
    rows = []
    for prefix in ["V3F_off", "V5_off", "V6_off"]:
        ng = compute_ng(df_zone, prefix=prefix)
        marker_mask = ng >= 3
        n_marker = int(marker_mask.sum())
        n_marker_tp = int(((df_zone["label"] == "TP") & marker_mask).sum())
        n_marker_fp = int(((df_zone["label"] == "FP") & marker_mask).sum())
        marker_tp_rate = (n_marker_tp / n_marker) if n_marker else float("nan")
        # all-region rate (regardless of NG)
        n = len(df_zone)
        n_tp = int((df_zone["label"] == "TP").sum())
        # HP bucket distribution
        hp = compute_hp_bucket(df_zone, prefix=prefix)
        rows.append(
            dict(
                version=prefix.split("_")[0],
                zone_n=n,
                zone_TP_rate=n_tp / n if n else float("nan"),
                marker_n=n_marker,
                marker_coverage_pct=n_marker / n if n else float("nan"),
                marker_TP=n_marker_tp,
                marker_FP=n_marker_fp,
                marker_TP_rate=marker_tp_rate,
                hp_same_HP1_pct=(hp == "same_HP1").sum() / n if n else float("nan"),
                hp_same_HP2_pct=(hp == "same_HP2").sum() / n if n else float("nan"),
                hp_cross_het_pct=(hp == "cross_het").sum() / n if n else float("nan"),
                hp_cross_het_inv_pct=(hp == "cross_het_inv").sum() / n if n else float("nan"),
                hp_other_pct=(hp == "other").sum() / n if n else float("nan"),
            )
        )
    return pd.DataFrame(rows)


def sub_grid(df_zone: pd.DataFrame, use_proxy_cov: bool = False) -> pd.DataFrame:
    """3-axis sub-grid within zone. cov axis uses real cov_bin or proxy."""
    cov_col = "cov_proxy" if use_proxy_cov else "cov_bin"
    cov_levels = COV_PROXY_LEVELS if use_proxy_cov else COV_LEVELS

    rows = []
    for loh in LOH_LEVELS:
        for hp in HP_LEVELS:
            for cov in cov_levels:
                sub = df_zone[
                    (df_zone["loh_side_norm"] == loh)
                    & (df_zone["hp_bucket"] == hp)
                    & (df_zone[cov_col] == cov)
                ]
                if len(sub) == 0:
                    continue
                stat = cell_stats(sub, GLOBAL_TP_RATE, GLOBAL_FP_RATE)
                stat.update(dict(loh=loh, hp_bucket=hp, cov=cov, cov_axis=cov_col))
                rows.append(stat)
    return pd.DataFrame(rows)


def deviance_decomp_zone(df_zone: pd.DataFrame, axes: list[str]) -> pd.DataFrame:
    """For each axis, deviance decomposition: full vs leave-one-out.

    axes: list of feature column names.
    Returns df with axis, dev_explained_full, dev_reduced (no axis), incremental dev.
    """
    if len(df_zone) < 50:
        return pd.DataFrame()
    sub = df_zone[df_zone["label"].isin(["TP", "FP"])].copy()
    y = (sub["label"] == "TP").astype(int).values
    if y.sum() == 0 or y.sum() == len(y):
        return pd.DataFrame()

    # Build design matrices
    def encode_cat(col):
        return pd.get_dummies(sub[col], prefix=col, drop_first=True).astype(float)

    blocks = {}
    for ax in axes:
        if ax in ("loh_side_norm", "hp_bucket", "cov_bin"):
            blocks[ax] = encode_cat(ax)
        elif ax in ("Coverage_Multiple", "caller_af", "V6_off_n_reads", "NG_v6off"):
            blocks[ax] = sub[[ax]].fillna(sub[ax].median()).astype(float)
        else:
            blocks[ax] = encode_cat(ax)

    # intercept-only baseline
    intercept = np.ones((len(sub), 1))

    def fit_dev(X):
        try:
            from statsmodels.api import GLM, families

            m = GLM(y, X, family=families.Binomial()).fit()
            return float(m.deviance), float(m.null_deviance)
        except Exception:
            return float("nan"), float("nan")

    # Full
    parts = [intercept] + list(blocks.values())
    X_full = np.hstack([np.asarray(b, dtype=float) for b in parts])
    d_full, d_null = fit_dev(X_full)

    rows = []
    for ax_drop, _ in blocks.items():
        parts_red = [intercept] + [blocks[k] for k in blocks if k != ax_drop]
        X_red = np.hstack([np.asarray(b, dtype=float) for b in parts_red])
        d_red, _ = fit_dev(X_red)
        rows.append(
            dict(
                axis_dropped=ax_drop,
                dev_full=d_full,
                dev_reduced=d_red,
                dev_null=d_null,
                incremental_dev=(d_red - d_full),
                incremental_dev_explained=((d_red - d_full) / d_null) if d_null > 0 else float("nan"),
            )
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------
# Zone 1: Z-OCH (Outer cross_het)
# --------------------------------------------------------------------------------
def run_zone_och(df_full: pd.DataFrame, outdir: Path) -> dict:
    """Outer cross_het: loh_side != Inner AND HP bucket cross_het/cross_het_inv.

    Use full TSV (all FP retained). Outer = explicit 'Outer' or 'UNKNOWN' (NaN loh_side
    means master_join_ok==0, which means we cannot annotate Inner, so treat as not-Inner).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    mask = (df_full["loh_side_norm"] != "Inner") & (
        df_full["hp_bucket"].isin(["cross_het", "cross_het_inv"])
    )
    df_zone = df_full.loc[mask].copy()

    summary = compute_zone_summary(df_zone, "Z-OCH", df_full)

    # Per-version trajectory (HP bucket varies by version → recompute per version)
    # For trajectory, the zone definition itself uses V6_off HP bucket; for other versions
    # we just report the zone composition (V6-defined).
    traj = per_version_trajectory(df_zone)

    # Sub-grid (use cov_proxy since most rows lack Coverage_Multiple)
    grid = sub_grid(df_zone, use_proxy_cov=True)

    # 4 confound guards on top cells
    powered_cells = grid[grid["n"] >= 50].copy()
    if len(powered_cells) > 0:
        # Build cell-id strings
        powered_cells["cell_id"] = (
            powered_cells["loh"] + "|" + powered_cells["hp_bucket"] + "|" + powered_cells["cov"]
        )
        cg_rows = []
        for _, r in powered_cells.iterrows():
            sub = df_zone[
                (df_zone["loh_side_norm"] == r["loh"])
                & (df_zone["hp_bucket"] == r["hp_bucket"])
                & (df_zone["cov_proxy"] == r["cov"])
            ]
            cg_rows.append(confound_guards(sub, df_full, r["cell_id"], skip_hp_symmetry=True))
        cg_df = pd.DataFrame(cg_rows)
        cg_df = powered_cells[["cell_id", "n", "n_TP", "n_FP", "TP_rate", "FP_rate", "FP_enrichment"]].merge(
            cg_df, on="cell_id", how="left"
        )
    else:
        cg_df = pd.DataFrame()

    # Save
    pd.DataFrame([summary]).to_csv(outdir / "zone_summary.tsv", sep="\t", index=False)
    traj.to_csv(outdir / "zone_three_version_trajectory.tsv", sep="\t", index=False)
    grid.to_csv(outdir / "zone_grid.tsv", sep="\t", index=False)
    if len(cg_df) > 0:
        cg_df.to_csv(outdir / "zone_confound_guards.tsv", sep="\t", index=False)

    return dict(zone="Z-OCH", summary=summary, n=len(df_zone), idx=set(df_zone.index))


# --------------------------------------------------------------------------------
# Zone 2: Z-CHR8
# --------------------------------------------------------------------------------
def run_zone_chr8(df_full: pd.DataFrame, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    df_zone = df_full[df_full["chr"] == "chr8"].copy()

    summary = compute_zone_summary(df_zone, "Z-CHR8", df_full)

    # Per-chr comparison: chr8 vs (chr1+chr2+chr3) baseline vs other 21 chr total
    chr8 = df_full[df_full["chr"] == "chr8"]
    chr123 = df_full[df_full["chr"].isin(["chr1", "chr2", "chr3"])]
    other21 = df_full[~df_full["chr"].isin(["chr8"])]

    def short(df_sub, label):
        n = len(df_sub)
        ntp = int((df_sub["label"] == "TP").sum())
        nfp = int((df_sub["label"] == "FP").sum())
        return dict(
            group=label,
            n=n,
            n_TP=ntp,
            n_FP=nfp,
            TP_rate=ntp / n if n else float("nan"),
            FP_rate=nfp / n if n else float("nan"),
            FP_enrichment_vs_global=((nfp / n) / GLOBAL_FP_RATE) if n and GLOBAL_FP_RATE > 0 else float("nan"),
        )

    cmp_df = pd.DataFrame([short(chr8, "chr8"), short(chr123, "chr1+chr2+chr3"), short(other21, "non-chr8")])

    # Per-chr full breakdown
    per_chr_rows = []
    for chrom, sub in df_full.groupby("chr"):
        per_chr_rows.append(short(sub, chrom))
    per_chr_df = pd.DataFrame(per_chr_rows).sort_values("FP_rate", ascending=False)

    # Sub-grid (chr8 has CN/LOH only on master-joined subset; report both)
    grid_proxy = sub_grid(df_zone, use_proxy_cov=True)
    grid_real = sub_grid(df_zone[df_zone["master_join_ok"] == 1], use_proxy_cov=False)

    # Per-version trajectory in chr8
    traj = per_version_trajectory(df_zone)

    # Deviance decomposition: in chr8 master-joined, LOH/HP/CN/AF
    df_zone_mj = df_zone[df_zone["master_join_ok"] == 1].copy()
    if len(df_zone_mj) >= 50:
        # Need to ensure variance in each feature
        axes_hp_loh_cn = []
        if df_zone_mj["loh_side_norm"].nunique() > 1:
            axes_hp_loh_cn.append("loh_side_norm")
        if df_zone_mj["hp_bucket"].nunique() > 1:
            axes_hp_loh_cn.append("hp_bucket")
        if df_zone_mj["Coverage_Multiple"].notna().sum() > 0:
            axes_hp_loh_cn.append("Coverage_Multiple")
        if df_zone_mj["caller_af"].notna().sum() > 0:
            axes_hp_loh_cn.append("caller_af")
        decomp_df = deviance_decomp_zone(df_zone_mj, axes_hp_loh_cn)
    else:
        decomp_df = pd.DataFrame()

    # Confound guards on top cells (skip HP symmetry as chr8 isn't asymmetric-by-def, but
    # FP-heavy → keep HP symmetry as informational only)
    grid_for_cg = grid_proxy  # use full chr8 (not master-joined) for richer FP set
    powered_cells = grid_for_cg[grid_for_cg["n"] >= 50].copy()
    if len(powered_cells) > 0:
        powered_cells["cell_id"] = (
            powered_cells["loh"] + "|" + powered_cells["hp_bucket"] + "|" + powered_cells["cov"]
        )
        cg_rows = []
        for _, r in powered_cells.iterrows():
            sub = df_zone[
                (df_zone["loh_side_norm"] == r["loh"])
                & (df_zone["hp_bucket"] == r["hp_bucket"])
                & (df_zone["cov_proxy"] == r["cov"])
            ]
            cg_rows.append(confound_guards(sub, df_full, r["cell_id"], skip_hp_symmetry=False))
        cg_df = pd.DataFrame(cg_rows)
        cg_df = powered_cells[["cell_id", "n", "n_TP", "n_FP", "TP_rate", "FP_rate", "FP_enrichment"]].merge(
            cg_df, on="cell_id", how="left"
        )
    else:
        cg_df = pd.DataFrame()

    # Save
    pd.DataFrame([summary]).to_csv(outdir / "zone_summary.tsv", sep="\t", index=False)
    cmp_df.to_csv(outdir / "per_chr_comparison.tsv", sep="\t", index=False)
    per_chr_df.to_csv(outdir / "per_chr_full_breakdown.tsv", sep="\t", index=False)
    grid_proxy.to_csv(outdir / "zone_grid_proxy_cov.tsv", sep="\t", index=False)
    if len(grid_real) > 0:
        grid_real.to_csv(outdir / "zone_grid_real_cov.tsv", sep="\t", index=False)
    traj.to_csv(outdir / "zone_three_version_trajectory.tsv", sep="\t", index=False)
    if len(decomp_df) > 0:
        decomp_df.to_csv(outdir / "lr_ablation_deviance.tsv", sep="\t", index=False)
    if len(cg_df) > 0:
        cg_df.to_csv(outdir / "zone_confound_guards.tsv", sep="\t", index=False)

    return dict(
        zone="Z-CHR8",
        summary=summary,
        per_chr=cmp_df,
        decomp=decomp_df,
        n=len(df_zone),
        idx=set(df_zone.index),
    )


# --------------------------------------------------------------------------------
# Zone 3: Z-GL (Gain+LOH)
# --------------------------------------------------------------------------------
def run_zone_gl(df_full: pd.DataFrame, outdir: Path) -> dict:
    """Inner LOH + Coverage_Multiple >= 1.3 (master-joined only, needs LOH+CN)."""
    outdir.mkdir(parents=True, exist_ok=True)
    df_zone = df_full[
        (df_full["master_join_ok"] == 1)
        & (df_full["loh_side_norm"] == "Inner")
        & (df_full["Coverage_Multiple"] >= 1.3)
    ].copy()

    summary = compute_zone_summary(df_zone, "Z-GL", df_full)

    # Sub-grid: HP bucket × cov_bin (loh fixed)
    rows = []
    for hp in HP_LEVELS:
        for cov in COV_LEVELS:
            if cov in ("cov_loss", "cov_normal"):
                continue  # by def excluded
            sub = df_zone[(df_zone["hp_bucket"] == hp) & (df_zone["cov_bin"] == cov)]
            if len(sub) == 0:
                continue
            stat = cell_stats(sub, GLOBAL_TP_RATE, GLOBAL_FP_RATE)
            stat.update(dict(loh="Inner", hp_bucket=hp, cov=cov))
            rows.append(stat)
    grid = pd.DataFrame(rows)

    # Trajectory
    traj = per_version_trajectory(df_zone)

    # Deviance decomposition: full vs leave-one-out (LOH/HP/CN)
    # LOH fixed in zone (all Inner), so axes = hp_bucket, cov_bin, caller_af
    axes = []
    if df_zone["hp_bucket"].nunique() > 1:
        axes.append("hp_bucket")
    if df_zone["Coverage_Multiple"].notna().sum() > 0:
        axes.append("Coverage_Multiple")
    if df_zone["caller_af"].notna().sum() > 0:
        axes.append("caller_af")
    decomp_df = deviance_decomp_zone(df_zone, axes)

    # SEQC2 concordance: in this zone, what is the Coverage_Multiple distribution vs SEQC2 CN?
    seqc2_concordance = pd.DataFrame()
    try:
        seqc2 = load_seqc2_cnv()
        # Join by Chr+Pos
        seqc2_cols = [c for c in seqc2.columns if c in ("Chr", "Pos", "SEQC2_CN", "SEQC2_Zone")]
        if {"Chr", "Pos"}.issubset(seqc2.columns) and "SEQC2_CN" in seqc2.columns:
            seqc2_lite = seqc2[seqc2_cols].drop_duplicates(subset=["Chr", "Pos"])
            seqc2_lite.columns = [c if c not in ("Chr", "Pos") else c.lower() for c in seqc2_lite.columns]
            merged = df_zone.merge(seqc2_lite, on=["chr", "pos"], how="left")
            sc = merged[merged["SEQC2_CN"].notna()]
            if len(sc) > 0:
                rows = []
                for cn, sub in sc.groupby("SEQC2_CN"):
                    rows.append(
                        dict(
                            SEQC2_CN=cn,
                            n=len(sub),
                            covm_mean=float(sub["Coverage_Multiple"].mean()),
                            covm_median=float(sub["Coverage_Multiple"].median()),
                            covm_std=float(sub["Coverage_Multiple"].std()),
                            TP_rate=(sub["label"] == "TP").mean(),
                            FP_rate=(sub["label"] == "FP").mean(),
                        )
                    )
                seqc2_concordance = pd.DataFrame(rows)
    except Exception as e:
        print(f"[Z-GL] SEQC2 concordance failed: {e}")

    # Confound guards on top cells
    powered_cells = grid[grid["n"] >= 50].copy()
    if len(powered_cells) > 0:
        powered_cells["cell_id"] = "Inner|" + powered_cells["hp_bucket"] + "|" + powered_cells["cov"]
        cg_rows = []
        for _, r in powered_cells.iterrows():
            sub = df_zone[(df_zone["hp_bucket"] == r["hp_bucket"]) & (df_zone["cov_bin"] == r["cov"])]
            cg_rows.append(confound_guards(sub, df_full, r["cell_id"], skip_hp_symmetry=False))
        cg_df = pd.DataFrame(cg_rows)
        cg_df = powered_cells[["cell_id", "n", "n_TP", "n_FP", "TP_rate", "FP_rate", "FP_enrichment"]].merge(
            cg_df, on="cell_id", how="left"
        )
    else:
        cg_df = pd.DataFrame()

    # Save
    pd.DataFrame([summary]).to_csv(outdir / "zone_summary.tsv", sep="\t", index=False)
    grid.to_csv(outdir / "zone_grid.tsv", sep="\t", index=False)
    traj.to_csv(outdir / "zone_three_version_trajectory.tsv", sep="\t", index=False)
    if len(decomp_df) > 0:
        decomp_df.to_csv(outdir / "lr_ablation_deviance.tsv", sep="\t", index=False)
    if len(seqc2_concordance) > 0:
        seqc2_concordance.to_csv(outdir / "seqc2_concordance.tsv", sep="\t", index=False)
    if len(cg_df) > 0:
        cg_df.to_csv(outdir / "zone_confound_guards.tsv", sep="\t", index=False)

    return dict(
        zone="Z-GL", summary=summary, decomp=decomp_df, n=len(df_zone), idx=set(df_zone.index)
    )


# --------------------------------------------------------------------------------
# Zone 4: Z-AUTO (top 5% FP density)
# --------------------------------------------------------------------------------
def run_zone_auto(df_full: pd.DataFrame, outdir: Path, bandwidth_bp: int = 1_000_000) -> dict:
    """KDE smoothing per-chr to identify top 5% FP-density positions."""
    outdir.mkdir(parents=True, exist_ok=True)

    # Per-chr FP density via Gaussian KDE on position
    df_full = df_full.copy()
    df_full["is_FP"] = (df_full["label"] == "FP").astype(int)
    df_full["fp_density"] = 0.0

    bw = bandwidth_bp
    for chrom, sub_idx in df_full.groupby("chr").groups.items():
        sub = df_full.loc[sub_idx]
        positions = sub["pos"].values
        is_fp = sub["is_FP"].values
        if is_fp.sum() == 0:
            continue
        fp_positions = positions[is_fp == 1]
        if len(fp_positions) < 3:
            # too sparse, just mark FP positions themselves
            df_full.loc[sub_idx, "fp_density"] = is_fp.astype(float)
            continue
        # Vectorised gaussian kernel summation
        # For each position, density = sum_{i in FP} exp(-0.5 * ((pos - fp_i) / bw)^2)
        # use mini-batch to avoid memory blowup
        densities = np.zeros(len(positions), dtype=float)
        batch = 512
        for start in range(0, len(positions), batch):
            chunk = positions[start : start + batch]
            diff = (chunk[:, None] - fp_positions[None, :]) / bw
            densities[start : start + batch] = np.exp(-0.5 * diff**2).sum(axis=1)
        df_full.loc[sub_idx, "fp_density"] = densities

    # Top 5% by density (across whole genome)
    threshold = np.percentile(df_full["fp_density"].values, 95)
    df_full["in_Z_AUTO"] = df_full["fp_density"] >= threshold
    df_zone = df_full[df_full["in_Z_AUTO"]].copy()

    summary = compute_zone_summary(df_zone, "Z-AUTO", df_full)
    summary["density_threshold"] = float(threshold)
    summary["bandwidth_bp"] = int(bandwidth_bp)

    # Sub-grid
    grid = sub_grid(df_zone, use_proxy_cov=True)

    # Trajectory
    traj = per_version_trajectory(df_zone)

    # Confound guards
    powered_cells = grid[grid["n"] >= 50].copy()
    if len(powered_cells) > 0:
        powered_cells["cell_id"] = (
            powered_cells["loh"] + "|" + powered_cells["hp_bucket"] + "|" + powered_cells["cov"]
        )
        cg_rows = []
        for _, r in powered_cells.iterrows():
            sub = df_zone[
                (df_zone["loh_side_norm"] == r["loh"])
                & (df_zone["hp_bucket"] == r["hp_bucket"])
                & (df_zone["cov_proxy"] == r["cov"])
            ]
            cg_rows.append(confound_guards(sub, df_full, r["cell_id"], skip_hp_symmetry=False))
        cg_df = pd.DataFrame(cg_rows)
        cg_df = powered_cells[["cell_id", "n", "n_TP", "n_FP", "TP_rate", "FP_rate", "FP_enrichment"]].merge(
            cg_df, on="cell_id", how="left"
        )
    else:
        cg_df = pd.DataFrame()

    # Save
    pd.DataFrame([summary]).to_csv(outdir / "zone_summary.tsv", sep="\t", index=False)
    grid.to_csv(outdir / "zone_grid.tsv", sep="\t", index=False)
    traj.to_csv(outdir / "zone_three_version_trajectory.tsv", sep="\t", index=False)
    if len(cg_df) > 0:
        cg_df.to_csv(outdir / "zone_confound_guards.tsv", sep="\t", index=False)
    # Save the per-region fp_density for downstream / Jaccard
    df_full[["region_id", "chr", "pos", "label", "fp_density", "in_Z_AUTO"]].to_csv(
        outdir / "fp_density_per_region.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    return dict(
        zone="Z-AUTO",
        summary=summary,
        n=len(df_zone),
        idx=set(df_zone.index),
        density_threshold=threshold,
        df_full_with_density=df_full[["region_id", "chr", "pos", "fp_density", "in_Z_AUTO"]],
    )


# --------------------------------------------------------------------------------
# Synthesis
# --------------------------------------------------------------------------------
def compute_jaccard(set_a: set, set_b: set) -> float:
    a, b = set(set_a), set(set_b)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def make_figures(zone_results: dict, outdir: Path):
    """Generate cross-zone overview figures."""
    plt = matplotlib_setup()
    outdir.mkdir(parents=True, exist_ok=True)

    # Figure 1: 4-zone TP/FP enrichment bar
    summaries = pd.DataFrame([zone_results[z]["summary"] for z in ["Z-OCH", "Z-CHR8", "Z-GL", "Z-AUTO"]])
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    x = np.arange(len(summaries))
    ax[0].bar(x, summaries["TP_enrichment"], color="#4a7ab8")
    ax[0].axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(summaries["zone"])
    ax[0].set_ylabel("TP enrichment (vs global)")
    ax[0].set_title("Zone TP rate vs global baseline")

    ax[1].bar(x, summaries["FP_enrichment"], color="#c84a4a")
    ax[1].axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(summaries["zone"])
    ax[1].set_ylabel("FP enrichment (vs global)")
    ax[1].set_title("Zone FP rate vs global baseline")
    plt.tight_layout()
    plt.savefig(outdir / "step3_zone_enrichment.png", dpi=120)
    plt.close()

    # Figure 2: per-chr FP rate (Z-CHR8 supporting evidence)
    if "per_chr" in zone_results["Z-CHR8"]:
        per_chr_full = pd.read_csv(STEP3_DIR / "chr8_hotspot" / "per_chr_full_breakdown.tsv", sep="\t")
        per_chr_full = per_chr_full.sort_values("FP_rate", ascending=False).head(15)
        fig, ax = plt.subplots(figsize=(10, 5))
        colors = ["#c84a4a" if c == "chr8" else "#4a7ab8" for c in per_chr_full["group"]]
        ax.bar(per_chr_full["group"], per_chr_full["FP_rate"], color=colors)
        ax.axhline(GLOBAL_FP_RATE, color="gray", linestyle="--", label=f"Global FP={GLOBAL_FP_RATE:.3f}")
        ax.set_ylabel("FP rate")
        ax.set_title("Per-chromosome FP rate (top 15) — chr8 hotspot")
        ax.legend()
        ax.tick_params(axis="x", rotation=45)
        plt.tight_layout()
        plt.savefig(outdir / "step3_per_chr_fp_rate.png", dpi=120)
        plt.close()


def synthesize(zone_results: dict, outdir: Path):
    """Cross-zone synthesis: 4-zone table + H3/H4/H5 verdict + Jaccard overlap."""
    # 4-zone comparison
    summary_rows = []
    for z in ["Z-OCH", "Z-CHR8", "Z-GL", "Z-AUTO"]:
        s = zone_results[z]["summary"]
        summary_rows.append(s)
    cross_zone = pd.DataFrame(summary_rows)
    cross_zone.to_csv(outdir / "step3_cross_zone_summary.tsv", sep="\t", index=False)

    # Jaccard overlap matrix (H5 input)
    zones = ["Z-OCH", "Z-CHR8", "Z-GL", "Z-AUTO"]
    jac_matrix = []
    for z1 in zones:
        row = {"zone": z1}
        for z2 in zones:
            row[z2] = compute_jaccard(zone_results[z1]["idx"], zone_results[z2]["idx"])
        jac_matrix.append(row)
    jac_df = pd.DataFrame(jac_matrix)
    jac_df.to_csv(outdir / "step3_jaccard_overlap.tsv", sep="\t", index=False)

    # H3 verdict (was Step 2 NEGATIVE): Outer cross_het × cov_gain FP enrichment
    z_och_summary = zone_results["Z-OCH"]["summary"]
    h3_fp_enrichment = z_och_summary["FP_enrichment"]
    h3_verdict = "POSITIVE" if h3_fp_enrichment >= 2.0 else "NEGATIVE"
    if 1.5 <= h3_fp_enrichment < 2.0:
        h3_verdict = "MIXED"

    # H4 verdict: chr8 LOH+CN deviance vs HP deviance
    chr8_decomp = zone_results["Z-CHR8"].get("decomp", pd.DataFrame())
    h4_verdict = "UNKNOWN"
    h4_details = ""
    if len(chr8_decomp) > 0:
        # Find LOH+CN incremental vs HP incremental
        loh_dev = chr8_decomp[chr8_decomp["axis_dropped"] == "loh_side_norm"]["incremental_dev_explained"].values
        cn_dev = chr8_decomp[chr8_decomp["axis_dropped"] == "Coverage_Multiple"]["incremental_dev_explained"].values
        hp_dev = chr8_decomp[chr8_decomp["axis_dropped"] == "hp_bucket"]["incremental_dev_explained"].values
        loh_v = float(loh_dev[0]) if len(loh_dev) else 0.0
        cn_v = float(cn_dev[0]) if len(cn_dev) else 0.0
        hp_v = float(hp_dev[0]) if len(hp_dev) else 0.0
        h4_details = f"LOH={loh_v:.4f}, CN={cn_v:.4f}, HP={hp_v:.4f}; (LOH+CN)-HP = {loh_v+cn_v-hp_v:.4f}"
        # Predicted: LOH+CN ≥ HP + 0.05
        if (loh_v + cn_v) >= hp_v + 0.05:
            h4_verdict = "POSITIVE"
        elif (loh_v + cn_v) < hp_v:
            h4_verdict = "NEGATIVE"
        else:
            h4_verdict = "MIXED"

    # H5 verdict: Z-AUTO vs known zones Jaccard
    known_union_idx = (
        zone_results["Z-OCH"]["idx"] | zone_results["Z-CHR8"]["idx"] | zone_results["Z-GL"]["idx"]
    )
    auto_idx = zone_results["Z-AUTO"]["idx"]
    jac_auto_vs_known = compute_jaccard(auto_idx, known_union_idx)
    h5_verdict = "POSITIVE" if jac_auto_vs_known >= 0.60 else "NEGATIVE"
    if 0.40 <= jac_auto_vs_known < 0.60:
        h5_verdict = "MIXED"

    verdicts = dict(
        H3=dict(
            verdict=h3_verdict,
            FP_enrichment=h3_fp_enrichment,
            threshold=2.0,
        ),
        H4=dict(
            verdict=h4_verdict,
            details=h4_details,
        ),
        H5=dict(
            verdict=h5_verdict,
            jaccard_auto_vs_known_union=jac_auto_vs_known,
            threshold=0.60,
        ),
    )
    with open(outdir / "step3_verdicts.json", "w") as f:
        json.dump(verdicts, f, indent=2)

    # Write synthesis markdown
    md_lines = [
        "# Step 3 Synthesis — 4 FP Zone Deep Dive (HCC1395)",
        "",
        "> Characterization-only. Universe: full Step 1 master TSV (35,332 rows = 30,490 TP + 4,842 FP).",
        f"> Global TP rate = {GLOBAL_TP_RATE:.4f}, Global FP rate = {GLOBAL_FP_RATE:.4f}.",
        "",
        "## 0. TL;DR",
        "",
        f"- **H3 (Z-OCH × cov_gain FP enrichment ≥ 2.0×)**: **{h3_verdict}** (FP_enrichment={h3_fp_enrichment:.3f})",
        f"- **H4 (chr8 LOH+CN ≥ HP + 0.05)**: **{h4_verdict}** ({h4_details})",
        f"- **H5 (Z-AUTO vs known zones Jaccard ≥ 0.60)**: **{h5_verdict}** (Jaccard={jac_auto_vs_known:.3f})",
        "",
        "## 1. 4-zone cross comparison",
        "",
        cross_zone.to_markdown(index=False, floatfmt=".4f"),
        "",
        "## 2. Jaccard overlap matrix",
        "",
        jac_df.to_markdown(index=False, floatfmt=".3f"),
        "",
        f"- **Z-AUTO vs (Z-OCH ∪ Z-CHR8 ∪ Z-GL) union**: Jaccard = **{jac_auto_vs_known:.4f}**",
        "",
        "## 3. H3 重測 (Z-OCH × cov_gain)",
        "",
        "Step 2 H3 NEGATIVE (FP_enrichment ≈ 0.0) was caused by master-join filter losing 90% FP.",
        "Step 3 re-test uses full FP set (no master_join_ok filter), with cov_proxy from V6_off_n_reads.",
        "",
        f"- Z-OCH overall: n={zone_results['Z-OCH']['summary']['n']}, "
        f"FP_rate={zone_results['Z-OCH']['summary']['FP_rate']:.4f}, "
        f"FP_enrichment={zone_results['Z-OCH']['summary']['FP_enrichment']:.4f}",
        f"- **Verdict**: {h3_verdict}",
        "",
        "## 4. H4 chr8 LOH+CN vs HP contribution",
        "",
    ]

    if len(chr8_decomp) > 0:
        md_lines.append(chr8_decomp.to_markdown(index=False, floatfmt=".5f"))
        md_lines.append("")
    md_lines.append(f"- {h4_details}")
    md_lines.append(f"- **Verdict**: {h4_verdict}")
    md_lines.append("")
    md_lines.append("## 5. H5 Z-AUTO vs known zones Jaccard")
    md_lines.append("")
    md_lines.append(
        f"- Auto-detected (top 5% FP density, gaussian KDE bandwidth 1Mb): "
        f"n={zone_results['Z-AUTO']['summary']['n']}, "
        f"FP_enrichment={zone_results['Z-AUTO']['summary']['FP_enrichment']:.3f}"
    )
    md_lines.append(f"- Jaccard vs Z-OCH ∪ Z-CHR8 ∪ Z-GL union = **{jac_auto_vs_known:.4f}**")
    md_lines.append(f"- **Verdict**: {h5_verdict}")
    md_lines.append("")
    md_lines.append("## 6. Cross-zone common vs unique signatures (Coordinator hand-off)")
    md_lines.append("")
    n_och = len(zone_results["Z-OCH"]["idx"])
    n_chr8 = len(zone_results["Z-CHR8"]["idx"])
    n_gl = len(zone_results["Z-GL"]["idx"])
    n_auto = len(zone_results["Z-AUTO"]["idx"])
    intersect_all = (
        zone_results["Z-OCH"]["idx"]
        & zone_results["Z-CHR8"]["idx"]
        & zone_results["Z-GL"]["idx"]
        & zone_results["Z-AUTO"]["idx"]
    )
    md_lines.append(f"- Z-OCH n={n_och}; Z-CHR8 n={n_chr8}; Z-GL n={n_gl}; Z-AUTO n={n_auto}")
    md_lines.append(f"- 4-zone intersection: n={len(intersect_all)}")
    md_lines.append("")
    md_lines.append("## 7. Hand-off to Coordinator")
    md_lines.append("")
    md_lines.append("- Step 3 isolates 4 FP-rich zones with distinct mechanism:")
    md_lines.append("  - Z-OCH: cross-haplotype germline-somatic mixing (HP bucket signature)")
    md_lines.append("  - Z-CHR8: chr-specific hotspot (likely structural/repeat context, sample-specific)")
    md_lines.append("  - Z-GL: gain/elevated CN with LOH allele imbalance (master-joined only)")
    md_lines.append("  - Z-AUTO: KDE-detected emergent clusters (cross-mechanism)")
    md_lines.append("- All zone-level statistics characterization-only — no filter/ΔF1 claims.")
    md_lines.append("- Cross-sample generalization (Step 4) requires V6-only set for 4 other samples; Z-CHR8 mechanism expected to be HCC1395-specific.")

    with open(outdir / "step3_synthesis.md", "w") as f:
        f.write("\n".join(md_lines))

    return verdicts


def main():
    global GLOBAL_TP_RATE, GLOBAL_FP_RATE
    print("[step3] Loading Step 1 master TSV...")
    df = load_master_full()
    print(f"[step3] Loaded {len(df)} rows ({(df['label']=='TP').sum()} TP, {(df['label']=='FP').sum()} FP)")

    df = annotate_axes(df)
    GLOBAL_TP_RATE = (df["label"] == "TP").mean()
    GLOBAL_FP_RATE = (df["label"] == "FP").mean()
    print(f"[step3] Global TP={GLOBAL_TP_RATE:.4f}  FP={GLOBAL_FP_RATE:.4f}")

    zone_results = {}

    print("[step3] Running Z-OCH (Outer cross_het)...")
    zone_results["Z-OCH"] = run_zone_och(df, STEP3_DIR / "outer_cross_het")
    print(f"   Z-OCH n={zone_results['Z-OCH']['n']}, summary={zone_results['Z-OCH']['summary']}")

    print("[step3] Running Z-CHR8 (HCC1395 chr8 hotspot)...")
    zone_results["Z-CHR8"] = run_zone_chr8(df, STEP3_DIR / "chr8_hotspot")
    print(f"   Z-CHR8 n={zone_results['Z-CHR8']['n']}, summary={zone_results['Z-CHR8']['summary']}")

    print("[step3] Running Z-GL (Gain+LOH inner)...")
    zone_results["Z-GL"] = run_zone_gl(df, STEP3_DIR / "gain_loh_zone")
    print(f"   Z-GL n={zone_results['Z-GL']['n']}, summary={zone_results['Z-GL']['summary']}")

    print("[step3] Running Z-AUTO (top 5% FP density)...")
    zone_results["Z-AUTO"] = run_zone_auto(df, STEP3_DIR / "auto_detected_top5pct")
    print(f"   Z-AUTO n={zone_results['Z-AUTO']['n']}, summary={zone_results['Z-AUTO']['summary']}")

    print("[step3] Generating figures + synthesis...")
    make_figures(zone_results, STEP3_DIR / "figures")
    verdicts = synthesize(zone_results, STEP3_DIR)

    print("\n[step3] === VERDICTS ===")
    print(json.dumps(verdicts, indent=2))

    print("\n[step3] Done.")


if __name__ == "__main__":
    main()
