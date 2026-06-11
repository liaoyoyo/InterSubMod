#!/usr/bin/env python3
"""
44_q3_reverse_predict.py — Q3: reverse-predict ploidy/CN from methylation features.

Pre-reg H-B: methylation features CANNOT predict CN beyond coverage.
  Falsifier:  (full-model AUC) - (coverage-only AUC) > 0.05
  Threshold:  dAUC <= 0.02  ->  methylation does NOT encode CN.

Decisive task: BINARY  neutral (CN=2)  vs  gain (CN>2).  loss dropped (too few).
Sites:  primarily HP-axis significant (axis_type==HP & wilcoxon_p<0.05);
        also report all-HP (no sig filter) as robustness.

Features = [mean_delta, abs_delta, mean_germ_beta, mean_som_beta, n_paired_cpg]
           (+ blind_ari, silhouette where joinable from master_o2_error).
Coverage-only baseline = n_paired_cpg alone.

Models: LogisticRegression + RandomForest. 5-fold CV AUC for both.
dAUC = full_AUC - coverage_only_AUC.

auc-confound-guard:
  - within-coverage-stratum AUC (4 quartile bins of n_paired_cpg)
  - permutation test (1000x label shuffle) for methylation-ONLY model
    (drop coverage feature -> [mean_delta, abs_delta, mean_germ_beta, mean_som_beta])

Outputs (ALL detail to disk):
  json: genome_survey_v2/cn_confound/q3_reverse_predict.json
  fig:  figures/cn_confound/q3_reverse_predict.png

Conventions match:
  - HP-sig filter  -> 40_cn_annotate.py:115 (axis_type=="HP" & wilcoxon_p<0.05)
  - o2 join keys   -> (chrom, somatic_pos, axis)
Single sample (HCC1395), A pilot.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CN_DIR = f"{ROOT}/genome_survey_v2/cn_confound"
FIG_DIR = f"{ROOT}/figures/cn_confound"
MASTER_O1 = f"{CN_DIR}/master_o1_cn.tsv"
MASTER_O2 = f"{CN_DIR}/master_o2_error.tsv"
OUT_JSON = f"{CN_DIR}/q3_reverse_predict.json"
OUT_FIG = f"{FIG_DIR}/q3_reverse_predict.png"

os.makedirs(CN_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

RNG = 42
N_PERM = 1000
N_FOLDS = 5

BASE_FEATS = ["mean_delta", "abs_delta", "mean_germ_beta", "mean_som_beta", "n_paired_cpg"]
METH_ONLY_FEATS = ["mean_delta", "abs_delta", "mean_germ_beta", "mean_som_beta"]  # drop coverage
COVERAGE_FEAT = ["n_paired_cpg"]
O2_FEATS = ["blind_ari", "silhouette"]


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def make_models(fast: bool = False):
    """Return dict name -> sklearn estimator (fresh, unfitted).

    fast=True uses a lighter RF (fewer trees) for the 1000x permutation loop —
    permutation only needs a stable AUC estimate, not the production tree count.
    """
    lr = Pipeline(
        [
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=RNG)),
        ]
    )
    rf = RandomForestClassifier(
        n_estimators=120 if fast else 400,
        max_depth=12 if fast else None,
        min_samples_leaf=10 if fast else 5,
        class_weight="balanced",
        n_jobs=-1,
        random_state=RNG,
    )
    return {"logistic": lr, "rf": rf}


def cv_auc(X: pd.DataFrame, y: np.ndarray, estimator) -> float:
    """5-fold stratified CV AUC via out-of-fold predicted probabilities."""
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG)
    # clone-safe: cross_val_predict refits per fold
    proba = cross_val_predict(estimator, X, y, cv=skf, method="predict_proba", n_jobs=None)[:, 1]
    return float(roc_auc_score(y, proba))


def perm_model(kind: str):
    """Lightweight, single-threaded estimator for the parallel permutation loop.

    n_jobs=1 so that joblib parallelizes ACROSS permutations (not nested), which is
    far more efficient for the 1000x label-shuffle loop. RF is leaner (60 trees,
    depth-capped) — permutation only needs a stable null AUC estimate.
    """
    if kind == "logistic":
        return Pipeline(
            [
                ("scale", StandardScaler()),
                ("clf", LogisticRegression(max_iter=1000, class_weight="balanced", random_state=RNG)),
            ]
        )
    return RandomForestClassifier(
        n_estimators=60,
        max_depth=10,
        min_samples_leaf=15,
        class_weight="balanced",
        n_jobs=1,
        random_state=RNG,
    )


def _fast_cv_auc(Xarr: np.ndarray, y: np.ndarray, estimator, splits) -> float:
    """Manual 5-fold OOF AUC (no cross_val_predict overhead). estimator must be n_jobs=1."""
    oof = np.empty(len(y), dtype=float)
    for tr, te in splits:
        est = clone(estimator)
        est.fit(Xarr[tr], y[tr])
        oof[te] = est.predict_proba(Xarr[te])[:, 1]
    return float(roc_auc_score(y, oof))


def oof_proba(X: pd.DataFrame, y: np.ndarray, estimator) -> np.ndarray:
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG)
    return cross_val_predict(estimator, X, y, cv=skf, method="predict_proba", n_jobs=None)[:, 1]


def permutation_test(X: pd.DataFrame, y: np.ndarray, kind: str, n_perm: int = N_PERM, n_jobs: int = 24):
    """
    Permutation test for methylation-only model: shuffle labels n_perm times,
    recompute 5-fold CV AUC each time. Parallelized ACROSS permutations via joblib
    (each worker uses a single-threaded perm_model). Returns (obs, perm_aucs, p).
    p = (#{perm_auc >= observed} + 1) / (n_perm + 1).
    kind in {"logistic","rf"}.
    """
    Xarr = np.asarray(X, dtype=float)
    y = np.asarray(y)
    # fixed fold structure (same splits for obs + every perm, so AUCs are comparable)
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG)
    splits = list(skf.split(Xarr, y))

    est = perm_model(kind)
    obs = _fast_cv_auc(Xarr, y, est, splits)

    # pre-generate permuted label vectors with a single seeded RNG (reproducible)
    rng = np.random.default_rng(RNG)
    perm_labels = [rng.permutation(y) for _ in range(n_perm)]

    def _one(yp):
        # note: shuffling labels breaks stratification balance only trivially; we
        # reuse the original splits' index structure on permuted labels (standard).
        return _fast_cv_auc(Xarr, yp, perm_model(kind), splits)

    perm = Parallel(n_jobs=n_jobs, backend="loky")(delayed(_one)(yp) for yp in perm_labels)
    perm = np.asarray(perm, dtype=float)
    p = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    return obs, perm, float(p)


def stratum_auc(df: pd.DataFrame, y: np.ndarray, feats, estimator, n_bins: int = 4):
    """
    within-coverage-stratum AUC: bin by n_paired_cpg quartiles, compute CV AUC
    of the methylation-only feature set INSIDE each coverage bin (controls coverage).
    Uses methylation-only feats (excludes n_paired_cpg) so coverage is held ~constant.
    """
    cov = df["n_paired_cpg"].values
    # quartile edges (unique to avoid degenerate bins)
    qs = np.quantile(cov, np.linspace(0, 1, n_bins + 1))
    qs = np.unique(qs)
    bins = np.clip(np.digitize(cov, qs[1:-1], right=True), 0, len(qs) - 2)
    out = []
    meth_feats = [f for f in feats if f != "n_paired_cpg"]
    for b in range(len(qs) - 1):
        mask = bins == b
        yb = y[mask]
        Xb = df.loc[mask, meth_feats]
        n = int(mask.sum())
        n_pos = int(yb.sum())
        rec = {
            "bin": b,
            "cov_lo": float(qs[b]),
            "cov_hi": float(qs[b + 1]),
            "n": n,
            "n_pos_gain": n_pos,
            "n_neg_neutral": int(n - n_pos),
        }
        # need both classes & enough samples for stratified 5-fold
        if n_pos >= N_FOLDS and (n - n_pos) >= N_FOLDS and len(np.unique(yb)) == 2:
            try:
                rec["auc"] = cv_auc(Xb, yb, estimator)
            except Exception as e:  # noqa: BLE001
                rec["auc"] = None
                rec["err"] = str(e)
        else:
            rec["auc"] = None
            rec["err"] = "insufficient_per_class_for_5fold"
        out.append(rec)
    return out


def run_block(df: pd.DataFrame, label: str, has_o2: bool, run_perm: bool = True):
    """
    Run the full Q3 analysis on a dataframe already filtered to neutral/gain binary.
    Returns a dict of results.  has_o2=True adds blind_ari/silhouette to full feats.
    run_perm=False skips the 1000x permutation test (used for the large all-HP
    robustness block where dAUC + meth-only AUC already settle the verdict and the
    n=25K permutation loop would be prohibitively slow).
    """
    y = (df["cn_class"].values == "gain").astype(int)  # gain=1 (positive), neutral=0
    n = len(df)
    n_gain = int(y.sum())
    n_neutral = int(n - n_gain)

    full_feats = list(BASE_FEATS)
    if has_o2:
        full_feats = full_feats + O2_FEATS

    X_full = df[full_feats].copy()
    X_cov = df[COVERAGE_FEAT].copy()
    X_meth = df[METH_ONLY_FEATS].copy()

    models = make_models()
    res = {
        "label": label,
        "n": n,
        "n_gain": n_gain,
        "n_neutral": n_neutral,
        "prevalence_gain": float(n_gain / n),
        "full_features": full_feats,
        "coverage_feature": COVERAGE_FEAT,
        "meth_only_features": METH_ONLY_FEATS,
        "models": {},
        "dAUC": {},
        "confound_guard": {},
    }

    # --- full vs coverage-only AUC for both models ---
    for mname in ("logistic", "rf"):
        full_auc = cv_auc(X_full, y, make_models()[mname])
        cov_auc = cv_auc(X_cov, y, make_models()[mname])
        meth_auc = cv_auc(X_meth, y, make_models()[mname])
        dauc = full_auc - cov_auc
        res["models"][mname] = {
            "full_auc": full_auc,
            "coverage_only_auc": cov_auc,
            "meth_only_auc": meth_auc,
            "dAUC_full_minus_cov": dauc,
        }
        res["dAUC"][mname] = dauc

    # --- feature importances (RF on full feats, fit on all data) ---
    rf_full = make_models()["rf"]
    rf_full.fit(X_full, y)
    res["feature_importance_rf"] = dict(
        sorted(
            zip(full_feats, [float(v) for v in rf_full.feature_importances_]),
            key=lambda kv: -kv[1],
        )
    )
    # logistic coefficients (standardized) as secondary importance
    lr_full = make_models()["logistic"]
    lr_full.fit(X_full, y)
    coefs = lr_full.named_steps["clf"].coef_[0]
    res["logistic_coef_std"] = dict(
        sorted(zip(full_feats, [float(c) for c in coefs]), key=lambda kv: -abs(kv[1]))
    )

    # --- confound guard 1: within-coverage-stratum AUC (4 bins) ---
    # use RF (more flexible) on meth-only feats inside each coverage bin
    res["confound_guard"]["within_coverage_stratum_auc_rf"] = stratum_auc(
        df, y, BASE_FEATS, make_models()["rf"], n_bins=4
    )
    res["confound_guard"]["within_coverage_stratum_auc_logistic"] = stratum_auc(
        df, y, BASE_FEATS, make_models()["logistic"], n_bins=4
    )

    # --- confound guard 2: permutation test for METHYLATION-ONLY model ---
    # (drop coverage) -> if this AUC ~0.5 and perm p high, methylation has no CN signal.
    # Uses fast RF (120 trees) for the 1000x loop; logistic is already cheap.
    if run_perm:
        obs_lr, perm_lr, p_lr = permutation_test(X_meth, y, "logistic", N_PERM)
        obs_rf, perm_rf, p_rf = permutation_test(X_meth, y, "rf", N_PERM)
        res["confound_guard"]["permutation_meth_only"] = {
            "logistic": {
                "observed_auc": obs_lr,
                "perm_auc_mean": float(np.mean(perm_lr)),
                "perm_auc_std": float(np.std(perm_lr)),
                "perm_auc_p95": float(np.quantile(perm_lr, 0.95)),
                "p_value": p_lr,
                "n_perm": N_PERM,
            },
            "rf": {
                "observed_auc": obs_rf,
                "perm_auc_mean": float(np.mean(perm_rf)),
                "perm_auc_std": float(np.std(perm_rf)),
                "perm_auc_p95": float(np.quantile(perm_rf, 0.95)),
                "p_value": p_rf,
                "n_perm": N_PERM,
                "note": "lean RF (60 trees, depth 10) single-threaded, joblib-parallel across perms",
            },
        }
        res["_perm_rf"] = perm_rf.tolist()
        res["_perm_obs_rf"] = obs_rf
    else:
        # robustness block: report meth-only AUC without the heavy 1000x loop
        mo_lr = cv_auc(X_meth, y, make_models()["logistic"])
        mo_rf = cv_auc(X_meth, y, make_models()["rf"])
        res["confound_guard"]["permutation_meth_only"] = {
            "logistic": {"observed_auc": mo_lr, "p_value": None, "note": "perm skipped (large-n robustness block)"},
            "rf": {"observed_auc": mo_rf, "p_value": None, "note": "perm skipped (large-n robustness block)"},
        }

    # --- store OOF proba for ROC plotting (primary = HP-sig, logistic + rf) ---
    res["_roc"] = {
        "y": y.tolist(),
        "full_logistic": oof_proba(X_full, y, make_models()["logistic"]).tolist(),
        "cov_logistic": oof_proba(X_cov, y, make_models()["logistic"]).tolist(),
        "full_rf": oof_proba(X_full, y, make_models()["rf"]).tolist(),
        "cov_rf": oof_proba(X_cov, y, make_models()["rf"]).tolist(),
    }
    return res


# ----------------------------------------------------------------------------
# Verdict logic (per H-B)
# ----------------------------------------------------------------------------
def verdict_for(res: dict) -> dict:
    """
    H-B threshold: dAUC <= 0.02 -> methylation does NOT encode CN.
    Falsifier:     dAUC > 0.05  -> methylation DOES encode CN beyond coverage.
    Use the MAX dAUC across logistic & rf as the decisive (most generous to falsifier).
    Also confirm methylation-only permutation AUC ~ 0.5.
    """
    dauc_lr = res["dAUC"]["logistic"]
    dauc_rf = res["dAUC"]["rf"]
    dauc_max = max(dauc_lr, dauc_rf)
    perm = res["confound_guard"]["permutation_meth_only"]
    meth_only_auc_max = max(perm["logistic"]["observed_auc"], perm["rf"]["observed_auc"])
    p_vals = [perm["logistic"]["p_value"], perm["rf"]["p_value"]]
    p_vals = [p for p in p_vals if p is not None]
    perm_p_min = min(p_vals) if p_vals else None

    if dauc_max > 0.05:
        v = "REFUTE_HB_methylation_encodes_CN"  # falsifier triggered
    elif dauc_max <= 0.02:
        v = "CONFIRM_HB_methylation_does_NOT_encode_CN"
    else:
        v = "INCONCLUSIVE_gray_zone_0.02_to_0.05"

    return {
        "dAUC_logistic": dauc_lr,
        "dAUC_rf": dauc_rf,
        "dAUC_max": dauc_max,
        "meth_only_auc_max": meth_only_auc_max,
        "meth_only_perm_p_min": perm_p_min,
        "falsifier_threshold_gt": 0.05,
        "confirm_threshold_le": 0.02,
        "verdict": v,
    }


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def make_figure(primary: dict, all_hp: dict):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod")
    try:
        from scripts.lib.plot_setup import setup_plot_style

        setup_plot_style()
    except Exception:  # noqa: BLE001
        pass

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # (A) ROC: full vs coverage-only (HP-sig primary), logistic + rf -------------
    ax = axes[0, 0]
    roc = primary["_roc"]
    y = np.array(roc["y"])
    for key, lab, col, ls in [
        ("full_rf", "Full (RF)", "C0", "-"),
        ("cov_rf", "Coverage-only (RF)", "C0", "--"),
        ("full_logistic", "Full (LogReg)", "C1", "-"),
        ("cov_logistic", "Coverage-only (LogReg)", "C1", "--"),
    ]:
        fpr, tpr, _ = roc_curve(y, roc[key])
        auc = roc_auc_score(y, roc[key])
        ax.plot(fpr, tpr, col, linestyle=ls, lw=1.8, label=f"{lab}: AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], "k:", lw=1, alpha=0.5)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("(A) ROC: gain vs neutral  [HP-sig]")
    ax.legend(fontsize=8, loc="lower right")

    # (B) within-coverage-stratum AUC (RF, HP-sig) -------------------------------
    ax = axes[0, 1]
    strata = primary["confound_guard"]["within_coverage_stratum_auc_rf"]
    xs = [f"Q{s['bin']+1}\n[{s['cov_lo']:.0f},{s['cov_hi']:.0f}]\nn={s['n']}" for s in strata]
    aucs = [s["auc"] if s["auc"] is not None else np.nan for s in strata]
    bars = ax.bar(range(len(xs)), aucs, color="C2", alpha=0.8)
    for i, (b, a) in enumerate(zip(bars, aucs)):
        if not np.isnan(a):
            ax.text(b.get_x() + b.get_width() / 2, a + 0.01, f"{a:.2f}", ha="center", fontsize=8)
    ax.axhline(0.5, color="r", ls="--", lw=1, label="chance (0.5)")
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels(xs, fontsize=7)
    ax.set_ylabel("AUC (meth-only, within coverage bin)")
    ax.set_ylim(0, 1)
    ax.set_title("(B) Within-coverage-stratum AUC  [confound guard]")
    ax.legend(fontsize=8)

    # (C) feature importance (RF, HP-sig full model) -----------------------------
    ax = axes[1, 0]
    fi = primary["feature_importance_rf"]
    names = list(fi.keys())
    vals = list(fi.values())
    ypos = np.arange(len(names))
    ax.barh(ypos, vals, color="C3", alpha=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("RF feature importance")
    ax.set_title("(C) Feature importance  [HP-sig full model]")
    # annotate coverage feature
    for i, nm in enumerate(names):
        if nm == "n_paired_cpg":
            ax.text(vals[i] + 0.005, i, "<-coverage", va="center", fontsize=8, color="darkred")

    # (D) dAUC bar: full - coverage-only, both models, both site sets ------------
    ax = axes[1, 1]
    groups = ["HP-sig\nLogReg", "HP-sig\nRF", "all-HP\nLogReg", "all-HP\nRF"]
    dvals = [
        primary["dAUC"]["logistic"],
        primary["dAUC"]["rf"],
        all_hp["dAUC"]["logistic"],
        all_hp["dAUC"]["rf"],
    ]
    cols = ["C0" if d <= 0.02 else ("orange" if d <= 0.05 else "red") for d in dvals]
    bars = ax.bar(range(len(groups)), dvals, color=cols, alpha=0.85)
    for b, d in zip(bars, dvals):
        ax.text(
            b.get_x() + b.get_width() / 2,
            d + (0.001 if d >= 0 else -0.003),
            f"{d:+.3f}",
            ha="center",
            va="bottom" if d >= 0 else "top",
            fontsize=8,
        )
    ax.axhline(0.02, color="green", ls="--", lw=1, label="confirm thr (0.02)")
    ax.axhline(0.05, color="red", ls="--", lw=1, label="falsifier thr (0.05)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, fontsize=8)
    ax.set_ylabel("dAUC = full - coverage-only")
    ax.set_title("(D) dAUC vs H-B thresholds")
    ax.legend(fontsize=8)

    fig.suptitle(
        "Q3 reverse-predict CN (gain vs neutral) from methylation | HCC1395 single-sample | H-B",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_FIG, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    df = pd.read_csv(MASTER_O1, sep="\t")

    # HP-axis subset
    hp = df[df["axis_type"] == "HP"].copy()
    hp_sig = hp[hp["wilcoxon_p"] < 0.05].copy()

    # join o2 (blind_ari, silhouette) on (chrom, somatic_pos, axis)
    o2 = pd.read_csv(MASTER_O2, sep="\t")[["chrom", "somatic_pos", "axis"] + O2_FEATS].copy()

    def prep_binary(sub, name):
        b = sub[sub["cn_class"].isin(["neutral", "gain"])].copy()
        # ensure no NaN in base feats (sanity)
        b = b.dropna(subset=BASE_FEATS)
        return b

    hp_sig_bin = prep_binary(hp_sig, "hp_sig")
    all_hp_bin = prep_binary(hp, "all_hp")

    # o2-joined version of HP-sig (subset where blind_ari/silhouette available)
    hp_sig_o2 = hp_sig_bin.merge(o2, on=["chrom", "somatic_pos", "axis"], how="inner")
    hp_sig_o2 = hp_sig_o2.dropna(subset=BASE_FEATS + O2_FEATS)

    print(f"[data] HP-sig binary n={len(hp_sig_bin)} "
          f"(gain={sum(hp_sig_bin.cn_class=='gain')}, neutral={sum(hp_sig_bin.cn_class=='neutral')})")
    print(f"[data] all-HP binary n={len(all_hp_bin)} "
          f"(gain={sum(all_hp_bin.cn_class=='gain')}, neutral={sum(all_hp_bin.cn_class=='neutral')})")
    print(f"[data] HP-sig+o2 (blind_ari/silhouette joined) n={len(hp_sig_o2)} "
          f"(gain={sum(hp_sig_o2.cn_class=='gain')}, neutral={sum(hp_sig_o2.cn_class=='neutral')})")

    # ---- run blocks ----
    print("[run] primary = HP-sig (5 meth features) + 1000x perm ...", flush=True)
    primary = run_block(hp_sig_bin, "HP_sig_5feat", has_o2=False, run_perm=True)
    print("[done] primary block", flush=True)

    print("[run] all-HP (robustness, perm skipped for n=25K) ...", flush=True)
    all_hp_res = run_block(all_hp_bin, "all_HP_5feat", has_o2=False, run_perm=False)
    print("[done] all-HP block", flush=True)

    o2_res = None
    # only run o2-joined block if both classes present with enough samples for 5-fold
    if len(hp_sig_o2) >= 50 and sum(hp_sig_o2.cn_class == "neutral") >= N_FOLDS \
            and sum(hp_sig_o2.cn_class == "gain") >= N_FOLDS:
        print(f"[run] HP-sig + o2 (7 features: 5 meth + blind_ari + silhouette) + 1000x perm ...", flush=True)
        o2_res = run_block(hp_sig_o2, "HP_sig_7feat_with_o2", has_o2=True, run_perm=True)
        print("[done] o2 block", flush=True)
    else:
        print(f"[skip] o2-joined block: insufficient (n={len(hp_sig_o2)})")

    # ---- verdicts ----
    v_primary = verdict_for(primary)
    v_all = verdict_for(all_hp_res)
    v_o2 = verdict_for(o2_res) if o2_res else None

    # ---- assemble output (strip heavy _roc/_perm arrays from json detail but keep summary) ----
    def strip_heavy(r):
        r2 = {k: v for k, v in r.items() if not k.startswith("_")}
        return r2

    out = {
        "task": "Q3_reverse_predict_CN_from_methylation",
        "sample": "HCC1395",
        "scope": "single_sample_A_pilot",
        "hypothesis": "H-B: methylation features CANNOT predict CN beyond coverage",
        "pre_reg": {
            "falsifier": "(full-model AUC) - (coverage-only AUC) > 0.05",
            "confirm_threshold": "dAUC <= 0.02 -> methylation does NOT encode CN",
            "target": "binary neutral(CN=2) vs gain(CN>2); loss dropped (too few)",
            "positive_class": "gain (CN>2) = 1",
            "features_full": BASE_FEATS,
            "features_coverage_only": COVERAGE_FEAT,
            "features_meth_only_for_perm": METH_ONLY_FEATS,
            "sig_site_def": "axis_type=='HP' & wilcoxon_p<0.05 (matches 40_cn_annotate.py:115)",
            "cv": f"{N_FOLDS}-fold stratified, OOF proba AUC",
            "n_perm": N_PERM,
        },
        "n_loss_dropped": {
            "hp_sig": int((hp_sig["cn_class"] == "loss").sum()),
            "all_hp": int((hp["cn_class"] == "loss").sum()),
        },
        "results": {
            "primary_HP_sig": strip_heavy(primary),
            "robustness_all_HP": strip_heavy(all_hp_res),
            "robustness_HP_sig_with_o2": strip_heavy(o2_res) if o2_res else None,
        },
        "verdicts": {
            "primary_HP_sig": v_primary,
            "robustness_all_HP": v_all,
            "robustness_HP_sig_with_o2": v_o2,
        },
        "decision_rule_result": {
            "decisive_block": "primary_HP_sig",
            "dAUC_max": v_primary["dAUC_max"],
            "confirm_threshold_le": 0.02,
            "falsifier_threshold_gt": 0.05,
            "met": v_primary["dAUC_max"] <= 0.02,
            "verdict": v_primary["verdict"],
            "meth_only_perm_check": (
                f"meth-only AUC max={v_primary['meth_only_auc_max']:.3f}, "
                f"perm p_min={v_primary['meth_only_perm_p_min']:.4f}"
            ),
        },
        "outputs": {"json": OUT_JSON, "figure": OUT_FIG},
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[write] {OUT_JSON}")

    make_figure(primary, all_hp_res)
    print(f"[write] {OUT_FIG}")

    # ---- compact stdout summary ----
    print("\n" + "=" * 60)
    print("COMPACT_SUMMARY")
    print("=" * 60)
    print(json.dumps(
        {
            "primary_dAUC_logistic": round(v_primary["dAUC_logistic"], 4),
            "primary_dAUC_rf": round(v_primary["dAUC_rf"], 4),
            "primary_dAUC_max": round(v_primary["dAUC_max"], 4),
            "primary_full_auc_rf": round(primary["models"]["rf"]["full_auc"], 4),
            "primary_cov_auc_rf": round(primary["models"]["rf"]["coverage_only_auc"], 4),
            "primary_full_auc_lr": round(primary["models"]["logistic"]["full_auc"], 4),
            "primary_cov_auc_lr": round(primary["models"]["logistic"]["coverage_only_auc"], 4),
            "meth_only_auc_rf": round(primary["models"]["rf"]["meth_only_auc"], 4),
            "meth_only_auc_lr": round(primary["models"]["logistic"]["meth_only_auc"], 4),
            "meth_only_perm_p_min": round(v_primary["meth_only_perm_p_min"], 4),
            "verdict": v_primary["verdict"],
            "all_HP_dAUC_max": round(v_all["dAUC_max"], 4),
            "all_HP_verdict": v_all["verdict"],
            "o2_dAUC_max": round(v_o2["dAUC_max"], 4) if v_o2 else None,
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
