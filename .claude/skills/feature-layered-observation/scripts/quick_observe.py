#!/usr/bin/env python3
"""
Quick single-feature observer — smoke test and Step 1 / Step 5 shortcut.

Usage:
    python quick_observe.py --feature <feature_id> [--mode paired|to|all]
                            [--master /path/to/master.tsv.gz]
                            [--stage full|step1|step5|step1+5]

Default master: research/feature_layered_observation/data/merged_with_vcf.tsv.gz
Default stage : step1+5 (fast smoke: global AUC + confound guard only)

Output:
    stdout TSV-style summary (feature, raw_auc [95% CI], resid_auc, af_bin_aucs,
                              perm_p, verdict_hint)
    --out DIR optional; writes summary tsv
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

DEFAULT_MASTER = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "feature_layered_observation/data/merged_with_vcf.tsv.gz"
)
CONFOUNDS = ["NumReads", "vcf_AF", "Coverage_Multiple"]


def hanley_mcneil_ci(auc: float, n_pos: int, n_neg: int, alpha: float = 0.05) -> tuple[float, float]:
    if n_pos < 5 or n_neg < 5:
        return (np.nan, np.nan)
    q1 = auc / (2 - auc)
    q2 = 2 * auc * auc / (1 + auc)
    se2 = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc * auc) +
           (n_neg - 1) * (q2 - auc * auc)) / (n_pos * n_neg)
    se = np.sqrt(max(se2, 0.0))
    z = stats.norm.ppf(1 - alpha / 2)
    return (max(auc - z * se, 0.0), min(auc + z * se, 1.0))


def safe_auc(label: np.ndarray, score: np.ndarray) -> tuple[float, int, int]:
    m = ~np.isnan(score) & ~np.isnan(label)
    if m.sum() < 10:
        return (np.nan, 0, 0)
    y = label[m].astype(int)
    s = score[m]
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos < 5 or n_neg < 5:
        return (np.nan, n_pos, n_neg)
    return (float(roc_auc_score(y, s)), n_pos, n_neg)


def step1_global(df: pd.DataFrame, feature: str) -> dict:
    y = df["tp_label"].values
    s = df[feature].astype(float).values
    auc, n_pos, n_neg = safe_auc(y, s)
    lo, hi = hanley_mcneil_ci(auc, n_pos, n_neg)
    m = ~np.isnan(s)
    tp_mask = (y == 1) & m
    fp_mask = (y == 0) & m
    if tp_mask.sum() >= 5 and fp_mask.sum() >= 5:
        mean_tp = float(s[tp_mask].mean())
        mean_fp = float(s[fp_mask].mean())
        pooled_sd = float(np.sqrt(((s[tp_mask].var(ddof=1) * (tp_mask.sum() - 1)) +
                                    (s[fp_mask].var(ddof=1) * (fp_mask.sum() - 1))) /
                                   (tp_mask.sum() + fp_mask.sum() - 2)))
        d = (mean_tp - mean_fp) / pooled_sd if pooled_sd > 0 else np.nan
        try:
            u, p = stats.mannwhitneyu(s[tp_mask], s[fp_mask], alternative="two-sided")
        except ValueError:
            p = np.nan
    else:
        mean_tp = mean_fp = d = p = np.nan
    return {
        "feature": feature,
        "auc": auc,
        "auc_lo": lo,
        "auc_hi": hi,
        "cohen_d": d,
        "mwu_p": float(p) if p is not None and not np.isnan(p) else np.nan,
        "mean_tp": mean_tp,
        "mean_fp": mean_fp,
        "n_tp": n_pos,
        "n_fp": n_neg,
    }


def within_group_ols_residual(df: pd.DataFrame, feature: str,
                               confounds: list[str]) -> np.ndarray:
    """Within-group OLS: fit separately on TP and FP, return residuals."""
    res = np.full(len(df), np.nan)
    y_all = df[feature].astype(float).values
    for label in [0, 1]:
        mask = (df["tp_label"].values == label)
        if mask.sum() < 10:
            continue
        X = df.loc[mask, confounds].astype(float).values
        y = y_all[mask]
        valid = ~np.isnan(y) & ~np.isnan(X).any(axis=1)
        if valid.sum() < 10:
            continue
        lr = LinearRegression().fit(X[valid], y[valid])
        pred = lr.predict(X[valid])
        sub = np.full(mask.sum(), np.nan)
        sub[valid] = y[valid] - pred
        res[mask] = sub
    return res


def step5_confound(df: pd.DataFrame, feature: str, confounds: list[str],
                   n_perm: int = 1000) -> dict:
    y = df["tp_label"].values
    s = df[feature].astype(float).values
    raw_auc, _, _ = safe_auc(y, s)

    # Gate 1: within-group OLS
    resid = within_group_ols_residual(df, feature, confounds)
    resid_auc, _, _ = safe_auc(y, resid)

    # Gate 2: AF-bin
    af_bin_aucs: dict[str, float] = {}
    if "AF_class" in df.columns:
        for name in ["Extreme", "Near-half", "Intermediate"]:
            mask = (df["AF_class"] == name).values
            if mask.sum() < 50:
                continue
            a, _, _ = safe_auc(y[mask], s[mask])
            af_bin_aucs[name] = a

    # Gate 3: permutation (on raw feature, for stability)
    valid = ~np.isnan(s)
    if valid.sum() > 100 and not np.isnan(raw_auc):
        rng = np.random.default_rng(42)
        s_v = s[valid]
        y_v = y[valid].astype(int)
        null = np.empty(n_perm)
        for i in range(n_perm):
            null[i] = roc_auc_score(rng.permutation(y_v), s_v)
        perm_p = float((null >= raw_auc).mean())
    else:
        perm_p = np.nan

    af_range = (max(af_bin_aucs.values()) - min(af_bin_aucs.values())
                if af_bin_aucs else np.nan)

    return {
        "feature": feature,
        "raw_auc": raw_auc,
        "resid_auc": resid_auc,
        "af_bin_aucs": af_bin_aucs,
        "af_bin_range": af_range,
        "perm_p": perm_p,
    }


def hint_verdict(global_res: dict, confound_res: dict) -> str:
    raw = global_res["auc"]
    if np.isnan(raw):
        return "NO_DATA"
    if raw < 0.50:
        return "ANTI_SIGNAL"
    if raw < 0.53:
        return "NEGATIVE (near-random)"
    if raw < 0.58:
        return "NEGATIVE (below Beyond-AUC ceiling)"
    resid = confound_res["resid_auc"]
    if np.isnan(resid) or resid < 0.55:
        return "CONFOUND_COLLAPSED (resid<0.55)"
    af_range = confound_res["af_bin_range"]
    if not np.isnan(af_range) and af_range >= 0.10:
        return "CONFOUND_COLLAPSED (AF-driven)"
    perm_p = confound_res["perm_p"]
    if not np.isnan(perm_p) and perm_p >= 0.05:
        return "NEGATIVE (perm p>=0.05)"
    # passes step1/5 only; full verdict needs step3+6
    if raw >= 0.65:
        return "POSITIVE_CANDIDATE (needs step3+6)"
    return "CONDITIONAL_POSITIVE_CANDIDATE (needs step3+6)"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--feature", required=True)
    ap.add_argument("--mode", choices=["paired", "to", "all"], default="all")
    ap.add_argument("--master", default=str(DEFAULT_MASTER))
    ap.add_argument("--stage", choices=["full", "step1", "step5", "step1+5"],
                    default="step1+5")
    ap.add_argument("--sample-subset", default=None,
                    help="comma-separated samples; default all 7")
    ap.add_argument("--out", default=None, help="optional output TSV path")
    ap.add_argument("--n-perm", type=int, default=1000)
    args = ap.parse_args()

    master_path = Path(args.master)
    if not master_path.exists():
        print(f"[err] master not found: {master_path}", file=sys.stderr)
        return 2

    usecols_probe = pd.read_csv(master_path, sep="\t", nrows=1, low_memory=False)
    available = set(usecols_probe.columns)
    if args.feature not in available:
        print(f"[err] feature '{args.feature}' not in master columns", file=sys.stderr)
        print(f"      master has {len(available)} cols; sample: "
              f"{sorted(list(available))[:20]}", file=sys.stderr)
        return 3

    usecols = ["sample", "mode", "tp_label", "AF_class", args.feature]
    for c in CONFOUNDS:
        if c in available and c != args.feature:
            usecols.append(c)

    print(f"[load] {master_path.name} cols={usecols}", file=sys.stderr)
    df = pd.read_csv(master_path, sep="\t", usecols=usecols, low_memory=False)
    print(f"[load] rows={len(df):,}", file=sys.stderr)

    if args.mode != "all":
        mode_map = {"paired": ["paired_full", "paired_pileup", "paired"],
                    "to": ["to_pileup", "to"]}
        df = df[df["mode"].isin(mode_map[args.mode])].copy()
        print(f"[filter] mode={args.mode} rows={len(df):,}", file=sys.stderr)

    if args.sample_subset:
        samples = [s.strip() for s in args.sample_subset.split(",")]
        df = df[df["sample"].isin(samples)].copy()
        print(f"[filter] samples={samples} rows={len(df):,}", file=sys.stderr)

    df["tp_label"] = df["tp_label"].astype(int)
    df = df.dropna(subset=[args.feature])
    print(f"[after dropna {args.feature}] rows={len(df):,}", file=sys.stderr)

    g = step1_global(df, args.feature)
    c = (step5_confound(df, args.feature, [x for x in CONFOUNDS if x in df.columns],
                        n_perm=args.n_perm)
         if args.stage in ("step5", "step1+5", "full") else None)

    print()
    print(f"=== quick_observe · feature={args.feature} · mode={args.mode} ===")
    print(f"Rows analyzed : {len(df):,}   TP={g['n_tp']:,}   FP={g['n_fp']:,}")
    print(f"Mean TP / FP  : {g['mean_tp']:.4f} / {g['mean_fp']:.4f}   "
          f"Cohen's d = {g['cohen_d']:+.3f}")
    print(f"MW-U p        : {g['mwu_p']:.3e}")
    print(f"RAW AUC       : {g['auc']:.4f}  [95% CI {g['auc_lo']:.4f}, {g['auc_hi']:.4f}]")
    if c is not None:
        print(f"Resid AUC     : {c['resid_auc']:.4f}   "
              f"(confounds: {[x for x in CONFOUNDS if x in df.columns]})")
        print(f"AF-bin AUC    : {c['af_bin_aucs']}")
        print(f"AF-bin range  : {c['af_bin_range']:.4f}"
              if not np.isnan(c['af_bin_range']) else "AF-bin range  : n/a")
        print(f"Permutation p : {c['perm_p']:.4f}" if not np.isnan(c['perm_p'])
              else "Permutation p : n/a")
        print(f"Verdict hint  : {hint_verdict(g, c)}")
    print()

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        rec = {**g}
        if c is not None:
            rec.update({
                "raw_auc": c["raw_auc"],
                "resid_auc": c["resid_auc"],
                "af_bin_aucs": str(c["af_bin_aucs"]),
                "af_bin_range": c["af_bin_range"],
                "perm_p": c["perm_p"],
                "verdict_hint": hint_verdict(g, c),
            })
        pd.DataFrame([rec]).to_csv(out, sep="\t", index=False)
        print(f"[out] wrote {out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
