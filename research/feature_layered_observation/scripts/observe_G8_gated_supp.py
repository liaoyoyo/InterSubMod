#!/usr/bin/env python3
"""G8 supplementary gated analysis (task retry minimal version).

Adds PerCpgASM_Valid==True gate (required for Fisher_* validity) and produces
the four core figures requested by the task:
  fig01_G8_global_auc_bar.png
  fig02_G8_per_sample_heatmap.png
  fig05_G8_confound_residualized.png
  fig_fisher_frac_sig_paired.png

Reads pre-enriched tsv produced by observe_G8_entropy.py
(data/G8/g8_enriched.tsv.gz) so we avoid re-loading all raw CSVs.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
DATA_DIR = ROOT / "research/feature_layered_observation/data/G8"
FIG_DIR = ROOT / "research/feature_layered_observation/figures/G8_entropy"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ENRICHED = DATA_DIR / "g8_enriched.tsv.gz"

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]

G8_FEATURES = [
    "NME_HP1", "NME_HP2", "Entropy_Imbalance",
    "Epipoly_HP1", "Epipoly_HP2", "Epipoly_Delta",
    "PerCpgASM_Valid",
    "Fisher_N_Sig", "Fisher_Frac_Sig", "Fisher_N_Tested", "Fisher_MaxNegLogFDR",
]
CONTINUOUS = [f for f in G8_FEATURES if f != "PerCpgASM_Valid"]

PALETTE = {"tp": "#5B8DB7", "fp": "#D55E00", "accent": "#A85540", "grey": "#8A8A8A"}


# ----- stat helpers -----

def _auc_hanley(y: np.ndarray, s: np.ndarray) -> Tuple[float, float, float]:
    y = np.asarray(y); s = np.asarray(s, dtype=float)
    mask = ~np.isnan(s) & ~pd.isna(y)
    y = y[mask]; s = s[mask]
    if len(y) < 5 or y.sum() == 0 or (1 - y).sum() == 0:
        return (np.nan, np.nan, np.nan)
    pos = s[y == 1]; neg = s[y == 0]
    x = np.concatenate([pos, neg])
    r = stats.rankdata(x)
    r_pos = r[: len(pos)].sum()
    n1 = len(pos); n0 = len(neg)
    auc = (r_pos - n1 * (n1 + 1) / 2) / (n1 * n0)
    q1 = auc / (2 - auc); q2 = (2 * auc * auc) / (1 + auc)
    var = (auc * (1 - auc) + (n1 - 1) * (q1 - auc ** 2) + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0)
    se = np.sqrt(max(var, 0))
    return (auc, auc - 1.96 * se, auc + 1.96 * se)


def _within_group_resid(d: pd.DataFrame, feat: str, covars: List[str], keys: List[str]) -> np.ndarray:
    out = pd.Series(np.nan, index=d.index, dtype=float)
    y_full = pd.to_numeric(d[feat], errors="coerce")
    for _, g in d.groupby(keys, observed=True):
        idx = g.index
        X = g[covars].apply(pd.to_numeric, errors="coerce")
        y = y_full.loc[idx]
        m = X.notna().all(axis=1) & y.notna()
        if m.sum() < 30 or X.loc[m].std().min() == 0:
            continue
        Xm = np.column_stack([np.ones(m.sum()), X.loc[m].values])
        ym = y.loc[m].values
        try:
            beta, *_ = np.linalg.lstsq(Xm, ym, rcond=None)
            out.loc[m[m].index] = ym - Xm @ beta
        except np.linalg.LinAlgError:
            continue
    return out.values


# ----- main -----

def main():
    print("[gated] loading enriched ...")
    df = pd.read_csv(ENRICHED, sep="\t", low_memory=False)
    print(f"[gated] shape: {df.shape}")

    # Gate: PerCpgASM_Valid == True
    if "PerCpgASM_Valid" in df.columns:
        mask = df["PerCpgASM_Valid"].map({True: True, False: False,
                                          "True": True, "False": False})
        g = df[mask == True].copy()
    else:
        g = df.copy()
    print(f"[gated] after PerCpgASM_Valid=True: {len(g)} / {len(df)} rows "
          f"({len(g) / max(len(df), 1):.1%})")

    # -------- Step 1: gated global AUC per feature --------
    rows = []
    for f in G8_FEATURES:
        if f not in g.columns:
            continue
        if f == "PerCpgASM_Valid":
            continue  # the gate itself
        s = pd.to_numeric(g[f], errors="coerce")
        y = g["tp_label"].astype(int).values
        a, lo, hi = _auc_hanley(y, s.values)
        pos = s[y == 1].dropna(); neg = s[y == 0].dropna()
        try:
            _, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
        except ValueError:
            p = np.nan
        rows.append(dict(feature=f, n_tp=int(len(pos)), n_fp=int(len(neg)),
                         mean_tp=float(pos.mean()), mean_fp=float(neg.mean()),
                         auc=a, auc_lo=lo, auc_hi=hi, mwu_p=p))
    auc_tbl = pd.DataFrame(rows).sort_values("auc", ascending=False)
    auc_tbl.to_csv(DATA_DIR / "G8_auc_table_gated.tsv", sep="\t", index=False)
    print(auc_tbl.to_string(index=False))

    # fig01: bar chart
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=120)
    d = auc_tbl.dropna(subset=["auc"]).sort_values("auc")
    y_pos = np.arange(len(d))
    xerr = np.vstack([(d["auc"] - d["auc_lo"]).clip(lower=0).values,
                      (d["auc_hi"] - d["auc"]).clip(lower=0).values])
    cols = [PALETTE["accent"] if v >= 0.58 else PALETTE["grey"] for v in d["auc"]]
    ax.barh(y_pos, d["auc"], xerr=xerr, color=cols, ecolor="#333", capsize=3)
    ax.set_yticks(y_pos); ax.set_yticklabels(d["feature"])
    ax.axvline(0.5, ls="--", c="#888", lw=1)
    ax.axvline(0.58, ls="--", c="#c00", lw=1, label="Beyond-AUC 0.58")
    ax.set_xlabel("AUC (TP vs FP, PerCpgASM_Valid=True gate, Hanley-McNeil 95% CI)")
    ax.set_xlim(0.35, 0.85)
    ax.set_title("G8 · Global AUC (gated on PerCpgASM_Valid)")
    ax.legend(loc="lower right")
    for i, (a_, f) in enumerate(zip(d["auc"].values, d["feature"].values)):
        ax.text(a_ + 0.005, i, f"{a_:.3f}", va="center", fontsize=7)
    fig.tight_layout(); fig.savefig(FIG_DIR / "fig01_G8_global_auc_bar.png"); plt.close(fig)

    # -------- Step 4: per-sample × mode AUC heatmap --------
    ps_rows = []
    for f in CONTINUOUS:
        if f not in g.columns:
            continue
        for (samp, mode), sub in g.groupby(["sample", "mode"]):
            s = pd.to_numeric(sub[f], errors="coerce")
            if s.notna().sum() < 30:
                continue
            y = sub["tp_label"].astype(int).values
            n_tp = int((y == 1).sum()); n_fp = int((y == 0).sum())
            a, lo, hi = _auc_hanley(y, s.values)
            ps_rows.append(dict(feature=f, sample=samp, mode=mode,
                                n=len(sub), n_tp=n_tp, n_fp=n_fp,
                                auc=a, auc_lo=lo, auc_hi=hi))
    ps = pd.DataFrame(ps_rows)
    ps.to_csv(DATA_DIR / "step3_per_sample_auc_gated.tsv", sep="\t", index=False)

    # Heatmap: rows = sample×mode, columns = features
    ps["cell_id"] = ps["sample"].astype(str) + " · " + ps["mode"]
    cell_order = [f"{s} · {m}" for m in ["paired_full", "to_pileup"] for s in SAMPLE_ORDER]
    feats_plot = [f for f in ["Fisher_Frac_Sig", "Fisher_N_Sig", "Fisher_MaxNegLogFDR",
                              "Fisher_N_Tested", "Epipoly_HP1", "Epipoly_HP2",
                              "Epipoly_Delta", "NME_HP1", "NME_HP2",
                              "Entropy_Imbalance"] if f in ps["feature"].unique()]
    pivot = ps.pivot_table(index="cell_id", columns="feature", values="auc", aggfunc="first")
    pivot = pivot.reindex([c for c in cell_order if c in pivot.index])[feats_plot]
    fig, ax = plt.subplots(figsize=(1.1 * len(feats_plot) + 3, 0.4 * len(pivot) + 2), dpi=120)
    arr = pivot.values
    im = ax.imshow(arr, cmap="RdBu_r", vmin=0.35, vmax=0.75, aspect="auto")
    ax.set_xticks(range(len(feats_plot))); ax.set_xticklabels(feats_plot, rotation=40, ha="right")
    ax.set_yticks(range(len(pivot.index))); ax.set_yticklabels(pivot.index, fontsize=7)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if (v < 0.45 or v > 0.65) else "#111")
    plt.colorbar(im, ax=ax, shrink=0.7, label="AUC")
    ax.set_title("G8 · Per-(sample × mode) AUC heatmap (gated on PerCpgASM_Valid=True)")
    fig.tight_layout(); fig.savefig(FIG_DIR / "fig02_G8_per_sample_heatmap.png"); plt.close(fig)

    # -------- Step 5: residualization on vcf_AF + NumReads + Coverage_Multiple --------
    covars = [c for c in ["vcf_AF", "NumReads", "Coverage_Multiple"] if c in g.columns]
    for c in covars:
        g[c] = pd.to_numeric(g[c], errors="coerce")
    g["LOH_Subtype"] = g.get("LOH_Subtype", "None").fillna("None")

    conf_rows = []
    only_signif = auc_tbl[auc_tbl["auc"] >= 0.58]["feature"].tolist()
    # Per task: residualize ONLY for continuous features with AUC >= 0.58, but also
    # include Entropy_Imbalance (task asks specifically) and Fisher_Frac_Sig regardless.
    resid_targets = sorted(set(only_signif) | {"Entropy_Imbalance", "Fisher_Frac_Sig", "Epipoly_Delta",
                                               "Epipoly_HP1", "Epipoly_HP2"})
    resid_targets = [f for f in resid_targets if f in g.columns and f != "PerCpgASM_Valid"]

    for f in resid_targets:
        s = pd.to_numeric(g[f], errors="coerce")
        m = s.notna()
        if m.sum() < 200:
            continue
        a_raw, r_lo, r_hi = _auc_hanley(g.loc[m, "tp_label"].values, s.loc[m].values)
        resid = _within_group_resid(g, f, covars, ["sample", "mode", "LOH_Subtype"])
        ok = ~np.isnan(resid)
        a_res, res_lo, res_hi = _auc_hanley(g.loc[ok, "tp_label"].values, resid[ok])
        conf_rows.append(dict(feature=f, n_raw=int(m.sum()), n_res=int(ok.sum()),
                              auc_raw=a_raw, auc_raw_lo=r_lo, auc_raw_hi=r_hi,
                              auc_resid=a_res, auc_resid_lo=res_lo, auc_resid_hi=res_hi))
    conf = pd.DataFrame(conf_rows).sort_values("auc_raw", ascending=False)
    conf.to_csv(DATA_DIR / "step5_confound_gated.tsv", sep="\t", index=False)
    print("\n[gated] Step 5 confound:")
    print(conf.to_string(index=False))

    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=120)
    d = conf.dropna(subset=["auc_raw"]).sort_values("auc_raw")
    y_pos = np.arange(len(d))
    bh = 0.38
    ax.barh(y_pos - bh / 2, d["auc_raw"], bh, color=PALETTE["accent"], label="raw AUC (gated)")
    ax.barh(y_pos + bh / 2, d["auc_resid"], bh, color=PALETTE["grey"],
            label="within-(sample,mode,LOH) resid on vcf_AF+NumReads+CovM")
    ax.axvline(0.5, ls="--", c="#888"); ax.axvline(0.58, ls="--", c="#c00")
    ax.set_yticks(y_pos); ax.set_yticklabels(d["feature"], fontsize=9)
    ax.set_xlim(0.3, 0.85); ax.set_xlabel("AUC")
    ax.legend(loc="lower right", fontsize=8)
    for i, (r, rr) in enumerate(zip(d["auc_raw"].values, d["auc_resid"].values)):
        ax.text(r + 0.003, i - bh / 2, f"{r:.3f}", va="center", fontsize=7)
        if not np.isnan(rr):
            ax.text(rr + 0.003, i + bh / 2, f"{rr:.3f}", va="center", fontsize=7)
    ax.set_title("G8 · Confound guard: raw vs within-group residualized (gated)")
    fig.tight_layout(); fig.savefig(FIG_DIR / "fig05_G8_confound_residualized.png"); plt.close(fig)

    # -------- Fisher_Frac_Sig paired 7-sample focus fig --------
    sub = ps[(ps["feature"] == "Fisher_Frac_Sig") & (ps["mode"] == "paired_full")]
    sub = sub.set_index("sample").reindex(SAMPLE_ORDER)
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=120)
    x = np.arange(len(SAMPLE_ORDER))
    aucs = sub["auc"].values
    err = np.vstack([(sub["auc"] - sub["auc_lo"]).clip(lower=0).fillna(0).values,
                     (sub["auc_hi"] - sub["auc"]).clip(lower=0).fillna(0).values])
    cols = [PALETTE["accent"] if (v is not None and not np.isnan(v) and v >= 0.58)
            else PALETTE["grey"] for v in aucs]
    bars = ax.bar(x, aucs, 0.6, color=cols, yerr=err, ecolor="#333", capsize=3)
    for xi, v, ntp, nfp in zip(x, aucs, sub["n_tp"].values, sub["n_fp"].values):
        if not np.isnan(v):
            ax.text(xi, v + 0.015, f"{v:.3f}\nTP={int(ntp)} FP={int(nfp)}",
                    ha="center", fontsize=7)
    ax.axhline(0.5, ls="--", c="#888", lw=0.8)
    ax.axhline(0.58, ls="--", c="#c00", lw=0.8)
    ax.axhline(0.65, ls="--", c="#080", lw=0.8, label="Question threshold 0.65")
    ax.set_ylim(0.3, 0.9)
    ax.set_xticks(x); ax.set_xticklabels(SAMPLE_ORDER, rotation=30, ha="right")
    ax.set_ylabel("Fisher_Frac_Sig AUC (gated)")
    n_above = int((pd.Series(aucs) >= 0.65).sum())
    total = int(pd.Series(aucs).notna().sum())
    ax.set_title(f"G8 · Fisher_Frac_Sig paired_full 7-sample AUC\n"
                 f"{n_above}/{total} samples reach ≥0.65 threshold")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(FIG_DIR / "fig_fisher_frac_sig_paired.png"); plt.close(fig)

    print("[gated] DONE")


if __name__ == "__main__":
    main()
