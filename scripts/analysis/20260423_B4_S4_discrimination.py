#!/usr/bin/env python3
"""B4 | HCC1395 TO S4 (LOH=None + Extreme AF) secondary-discrimination pilot.

Hypothesis: The S4 bucket (LOH_Subtype=None + AF<0.1 or AF>0.9) has no
first-level discriminative power (TP rate = baseline ~71%). We test whether
Thread D's "Inner same-hap fingerprint" (NG=2 with single-haplotype bucket
occupancy {HP1+HP1S} or {HP2+HP2S}) plus auxiliary features can pick a
high-precision subset inside S4.

Pass criteria:
    LR AUC >= 0.65 -> POSITIVE (secondary discrimination possible)
    LR AUC <  0.60 -> NEGATIVE

Outputs:
    docs/experiments/in_progress/2026/04/figures/20260423_B4_S4/
        feature_auc_bar.png
        lr_roc_curve.png
        samehap_fp_tp_stack.png
    research/tpfp_loh_af_kde_discrimination/data/B4_S4_auc_summary.tsv
    research/tpfp_loh_af_kde_discrimination/data/B4_S4_lr_cv.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# ---------- Paths ----------
PHASE1_DIR = Path("/tmp/ism_hp_fix_phase1")
TP_CSV = PHASE1_DIR / "tp_off" / "significance_summary.csv"
FP_CSV = PHASE1_DIR / "fp_off" / "significance_summary.csv"
MASTER_TSV = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "tpfp_loh_af_kde_discrimination/data/master_extended.tsv.gz"
)
FIG_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/04/"
    "figures/20260423_B4_S4"
)
DATA_DIR = MASTER_TSV.parent

FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---------- Load & filter ----------
def load_s4_data() -> pd.DataFrame:
    """Load HCC1395 TO, merge AF from master_extended, restrict to S4."""
    tp = pd.read_csv(TP_CSV)
    tp["tp_label"] = 1
    fp = pd.read_csv(FP_CSV)
    fp["tp_label"] = 0
    sig = pd.concat([tp, fp], axis=0, ignore_index=True)

    master = pd.read_csv(MASTER_TSV, sep="\t", low_memory=False)
    master_sub = master[(master["sample"] == "HCC1395") & (master["mode"] == "to_pileup")].copy()

    # Merge AF and LOH_Subtype on (Chr, Pos). sig has LOH_Subtype already,
    # but use master's AF field (sig has AlleleDelta not caller AF).
    key_cols = ["Chr", "Pos"]
    merged = sig.merge(
        master_sub[key_cols + ["AF", "AF_class"]],
        on=key_cols,
        how="inner",
    )

    # S4 filter: LOH_Subtype=None (string "None") + AF extreme
    s4 = merged[
        (merged["LOH_Subtype"] == "None") & ((merged["AF"] < 0.1) | (merged["AF"] > 0.9))
    ].copy()
    return merged, s4


# ---------- Feature engineering ----------
def build_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add derived features incl. same-hap fingerprint."""
    df = df.copy()

    # Bucket count booleans (0/1 for having reads in that bucket)
    for b in ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]:
        df[f"has_{b}"] = (df[b].fillna(0) > 0).astype(int)

    # Same-hap fingerprint: NG>=2 AND (only HP1-family has reads) OR (only HP2-family has reads)
    # i.e. at least 2 buckets populated but all from the same haplotype family
    same_hp1 = (df["has_HPFineN_HP1"] == 1) & (df["has_HPFineN_HP1S"] == 1) & \
               (df["has_HPFineN_HP2"] == 0) & (df["has_HPFineN_HP2S"] == 0)
    same_hp2 = (df["has_HPFineN_HP2"] == 1) & (df["has_HPFineN_HP2S"] == 1) & \
               (df["has_HPFineN_HP1"] == 0) & (df["has_HPFineN_HP1S"] == 0)
    df["SameHapFingerprint"] = (same_hp1 | same_hp2).astype(int)

    # Ratio of HP1-family reads vs HP2-family (informative even at NG>=2)
    hp1_total = df["HPFineN_HP1"].fillna(0) + df["HPFineN_HP1S"].fillna(0)
    hp2_total = df["HPFineN_HP2"].fillna(0) + df["HPFineN_HP2S"].fillna(0)
    total = hp1_total + hp2_total
    with np.errstate(divide="ignore", invalid="ignore"):
        df["HP1_family_frac"] = np.where(total > 0, hp1_total / total, 0.5)
    # Distance from balanced (0.5) -- extreme fraction = single-hap dominance
    df["HP_family_imbalance"] = np.abs(df["HP1_family_frac"] - 0.5)

    return df


# ---------- Analysis ----------
# NOTE: Coverage_Multiple is EXCLUDED — tp_off run uses Diploid_Coverage_Used=115
# and fp_off run uses 113, causing deterministic run-level label leakage when
# combined with NumReads (LR AUC -> 1.000). Diploid_Coverage_Used, NumReads
# (when paired with any run-level normaliser), CovM_used are all contaminated.
# Keep only per-region biological features that do not depend on the run-level
# baseline constant.
CANDIDATE_FEATURES = [
    "HPFineNGroups",
    "SameHapFingerprint",
    "HP_family_imbalance",
    "HP1_family_frac",
    "HP_Ratio",
    "Quality_Score",
    "AlleleDelta",
    "HPFine_NGroups_CF",
    "AF",
    "NTumorReads",
    "NNormalReads",
    "NumCpGs",
    "Fisher_Frac_Sig",
    "Entropy_Imbalance",
    "HPMergedDelta",
    "HPFineF",
    "CramersV_HPFine",
]


def single_feature_auc(df: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    rows = []
    for f in CANDIDATE_FEATURES:
        if f not in df.columns:
            continue
        x = pd.to_numeric(df[f], errors="coerce")
        mask = x.notna()
        if mask.sum() < 50 or len(np.unique(y[mask])) < 2:
            continue
        x_valid = x[mask].values
        y_valid = y[mask]
        try:
            auc_pos = roc_auc_score(y_valid, x_valid)
        except ValueError:
            continue
        # Report MAX(AUC, 1-AUC) with direction
        auc = max(auc_pos, 1 - auc_pos)
        direction = "+" if auc_pos >= 0.5 else "-"
        rows.append({"feature": f, "auc": auc, "direction": direction, "n": int(mask.sum())})
    return pd.DataFrame(rows).sort_values("auc", ascending=False).reset_index(drop=True)


def lr_cv_auc(df: pd.DataFrame, y: np.ndarray, features: list[str]) -> dict:
    feats = [f for f in features if f in df.columns]
    X = df[feats].apply(pd.to_numeric, errors="coerce")
    valid = X.notna().all(axis=1)
    X = X.loc[valid].values
    y = y[valid.values]

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fold_aucs = []
    all_y, all_p = [], []
    for tr, te in skf.split(X, y):
        scaler = StandardScaler().fit(X[tr])
        Xtr, Xte = scaler.transform(X[tr]), scaler.transform(X[te])
        clf = LogisticRegression(max_iter=1000, class_weight="balanced", random_state=42)
        clf.fit(Xtr, y[tr])
        p = clf.predict_proba(Xte)[:, 1]
        fold_aucs.append(roc_auc_score(y[te], p))
        all_y.append(y[te])
        all_p.append(p)
    all_y = np.concatenate(all_y)
    all_p = np.concatenate(all_p)
    pooled_auc = roc_auc_score(all_y, all_p)
    return {
        "features_used": feats,
        "n_total": int(len(y)),
        "n_valid": int(valid.sum()),
        "fold_aucs": [float(a) for a in fold_aucs],
        "mean_auc": float(np.mean(fold_aucs)),
        "std_auc": float(np.std(fold_aucs)),
        "pooled_auc": float(pooled_auc),
        "y_pooled": all_y.tolist(),
        "p_pooled": all_p.tolist(),
    }


# ---------- Plots ----------
def plot_feature_auc_bar(auc_df: pd.DataFrame, out_path: Path) -> None:
    top = auc_df.head(10).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["tab:green" if d == "+" else "tab:orange" for d in top["direction"]]
    ax.barh(top["feature"], top["auc"], color=colors)
    ax.axvline(0.5, color="grey", ls="--", lw=1, label="Random")
    ax.axvline(0.6, color="red", ls=":", lw=1, label="0.60 threshold")
    ax.axvline(0.65, color="darkred", ls=":", lw=1, label="0.65 threshold")
    ax.set_xlabel("AUC (max of positive/negative direction)")
    ax.set_title("B4 | Single-feature AUC inside HCC1395 TO S4 bucket")
    ax.set_xlim(0.45, max(0.75, top["auc"].max() + 0.02))
    for i, (f, a) in enumerate(zip(top["feature"], top["auc"])):
        ax.text(a + 0.005, i, f"{a:.3f}", va="center", fontsize=8)
    ax.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_lr_roc(lr_result: dict, out_path: Path) -> None:
    y = np.asarray(lr_result["y_pooled"])
    p = np.asarray(lr_result["p_pooled"])
    fpr, tpr, _ = roc_curve(y, p)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr, tpr, lw=2, label=f"LR (pooled AUC={lr_result['pooled_auc']:.3f})")
    ax.plot([0, 1], [0, 1], "--", color="grey", lw=1, label="Random")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(
        f"B4 | Logistic Regression 5-fold CV ROC\n"
        f"mean±std AUC = {lr_result['mean_auc']:.3f}±{lr_result['std_auc']:.3f}"
    )
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_samehap_stack(df: pd.DataFrame, out_path: Path) -> None:
    # Stack: count by SameHapFingerprint x TP/FP
    ct = df.groupby(["SameHapFingerprint", "tp_label"]).size().unstack(fill_value=0)
    ct = ct.reindex([0, 1])
    labels = ["No (mixed-hap / single-bucket)", "Yes (same-hap fingerprint)"]
    tp_vals = ct.get(1, pd.Series([0, 0]))
    fp_vals = ct.get(0, pd.Series([0, 0]))
    x = np.arange(len(labels))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    w = 0.6
    ax1.bar(x, tp_vals.values, w, label="TP", color="tab:green")
    ax1.bar(x, fp_vals.values, w, bottom=tp_vals.values, label="FP", color="tab:orange")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=10)
    ax1.set_ylabel("Region count")
    ax1.set_title("Count stack")
    for i, (t, f) in enumerate(zip(tp_vals.values, fp_vals.values)):
        total = t + f
        ax1.text(i, total, f"n={total}\nTP={t}", ha="center", va="bottom", fontsize=8)
    ax1.legend()

    total_vals = tp_vals.values + fp_vals.values
    tp_rate = np.divide(tp_vals.values, total_vals, out=np.zeros_like(tp_vals.values, dtype=float),
                        where=total_vals > 0)
    ax2.bar(x, tp_rate, w, color="tab:blue")
    baseline = df["tp_label"].mean()
    ax2.axhline(baseline, color="red", ls="--", label=f"S4 baseline TP rate = {baseline:.3f}")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=10)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("TP rate")
    ax2.set_title("TP rate by same-hap fingerprint")
    for i, r in enumerate(tp_rate):
        ax2.text(i, r + 0.01, f"{r:.3f}", ha="center", fontsize=9)
    ax2.legend()

    fig.suptitle("B4 | Same-hap fingerprint TP/FP composition (HCC1395 TO, S4)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------- Main ----------
def main() -> int:
    print("[B4] Loading HCC1395 TO significance_summary + master_extended...")
    full, s4 = load_s4_data()
    print(f"  Full HCC1395 TO: n={len(full)} TP rate={full['tp_label'].mean():.4f}")
    print(
        f"  S4 (LOH=None & AF extreme): n={len(s4)}, "
        f"TP={s4['tp_label'].sum()}, FP={(s4['tp_label']==0).sum()}, "
        f"TP rate={s4['tp_label'].mean():.4f}"
    )

    s4 = build_features(s4)
    y = s4["tp_label"].values

    # Single-feature AUC
    auc_df = single_feature_auc(s4, y)
    print("\n[B4] Single-feature AUC (top 10):")
    print(auc_df.head(10).to_string(index=False))
    auc_df.to_csv(DATA_DIR / "B4_S4_auc_summary.tsv", sep="\t", index=False)

    # LR CV
    print("\n[B4] Logistic Regression 5-fold CV...")
    lr_result = lr_cv_auc(s4, y, CANDIDATE_FEATURES)
    print(
        f"  Features used ({len(lr_result['features_used'])}): {lr_result['features_used']}"
    )
    print(f"  n_valid: {lr_result['n_valid']}/{lr_result['n_total']}")
    print(f"  Fold AUCs: {[f'{a:.3f}' for a in lr_result['fold_aucs']]}")
    print(f"  Mean AUC: {lr_result['mean_auc']:.3f} ± {lr_result['std_auc']:.3f}")
    print(f"  Pooled AUC: {lr_result['pooled_auc']:.3f}")

    # Decision
    mean_auc = lr_result["mean_auc"]
    if mean_auc >= 0.65:
        verdict = "POSITIVE"
    elif mean_auc < 0.60:
        verdict = "NEGATIVE"
    else:
        verdict = "INCONCLUSIVE"
    print(f"\n[B4] VERDICT: {verdict} (mean LR AUC = {mean_auc:.3f})")

    # Save LR result (trim y/p for JSON)
    lr_json = {k: v for k, v in lr_result.items() if k not in ("y_pooled", "p_pooled")}
    lr_json["verdict"] = verdict
    lr_json["s4_n"] = int(len(s4))
    lr_json["s4_tp"] = int(s4["tp_label"].sum())
    lr_json["s4_fp"] = int((s4["tp_label"] == 0).sum())
    lr_json["s4_tp_rate"] = float(s4["tp_label"].mean())
    with (DATA_DIR / "B4_S4_lr_cv.json").open("w") as fh:
        json.dump(lr_json, fh, indent=2)

    # Plots
    print("\n[B4] Generating figures...")
    plot_feature_auc_bar(auc_df, FIG_DIR / "feature_auc_bar.png")
    plot_lr_roc(lr_result, FIG_DIR / "lr_roc_curve.png")
    plot_samehap_stack(s4, FIG_DIR / "samehap_fp_tp_stack.png")
    print(f"  Figures -> {FIG_DIR}")

    # Same-hap fingerprint breakdown
    same_hap_stats = s4.groupby("SameHapFingerprint").agg(
        n=("tp_label", "size"),
        tp=("tp_label", "sum"),
        tp_rate=("tp_label", "mean"),
    )
    print("\n[B4] Same-hap fingerprint breakdown:")
    print(same_hap_stats.to_string())

    return 0


if __name__ == "__main__":
    sys.exit(main())
