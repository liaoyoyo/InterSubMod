#!/usr/bin/env python3
"""Per-sample AUC of HP_Ratio_norm within Near-half / Intermediate AF."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/feature_layered_observation")
DATA = ROOT / "data/G5/master_g5.tsv.gz"
OUT  = ROOT / "data/G5/G5_HP_Ratio_norm_persample_AF.tsv"

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]
AF_ORDER = ["Extreme", "Near-half", "Intermediate"]


def auc(y, s):
    y = np.asarray(y); s = np.asarray(s, dtype=float)
    m = ~np.isnan(s)
    y = y[m]; s = s[m]
    if len(y) < 5 or (y == 1).sum() < 5 or (y == 0).sum() < 5:
        return np.nan, 0, 0
    r = stats.rankdata(s)
    n_pos = int((y == 1).sum()); n_neg = int((y == 0).sum())
    a = (r[y == 1].sum() - n_pos*(n_pos+1)/2) / (n_pos*n_neg)
    return a, n_pos, n_neg


def main():
    df = pd.read_csv(DATA, sep="\t", usecols=["sample", "mode", "AF_class",
                                              "tp_label", "HP_Ratio",
                                              "LOH_Subtype"],
                     low_memory=False)
    df = df.dropna(subset=["HP_Ratio", "tp_label"]).copy()
    df["HP_Ratio_norm"] = (df["HP_Ratio"] - 0.5).abs()
    rows = []
    for s in SAMPLE_ORDER:
        for mode in ["paired_full", "to_pileup"]:
            for af in AF_ORDER:
                m = ((df["sample"] == s) & (df["mode"] == mode) &
                     (df["AF_class"] == af))
                if m.sum() < 100:
                    continue
                sub = df.loc[m]
                a_raw, npp, nnn = auc(sub["tp_label"],
                                      sub["HP_Ratio"].to_numpy())
                a_norm, _, _ = auc(sub["tp_label"],
                                   sub["HP_Ratio_norm"].to_numpy())
                rows.append({"sample": s, "mode": mode, "AF_class": af,
                             "n": int(m.sum()),
                             "tp_rate": float(sub["tp_label"].mean()),
                             "auc_HP_Ratio": a_raw,
                             "auc_HP_Ratio_norm": a_norm,
                             "n_pos": npp, "n_neg": nnn})
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    # Also pivot direction consistency
    print("--- HP_Ratio_norm AUC by sample × mode × AF (n>=100) ---")
    piv = out.pivot_table(index=["sample", "mode"], columns="AF_class",
                          values="auc_HP_Ratio_norm", aggfunc="first")
    print(piv.to_string(float_format=lambda x: f"{x:.3f}"))
    # Count (AUC>0.65 in Near-half) / total samples
    near = out[out["AF_class"] == "Near-half"]
    inter = out[out["AF_class"] == "Intermediate"]
    def summary(frame, name):
        if len(frame) == 0:
            print(f"{name}: empty"); return
        n_total = len(frame)
        n_gt58 = (frame["auc_HP_Ratio_norm"] > 0.58).sum()
        n_gt65 = (frame["auc_HP_Ratio_norm"] > 0.65).sum()
        print(f"{name}: {n_gt58}/{n_total} AUC_norm>0.58, {n_gt65}/{n_total} AUC_norm>0.65")
        print(f"  median AUC_norm: {frame['auc_HP_Ratio_norm'].median():.3f}")
    summary(near, "Near-half")
    summary(inter, "Intermediate")
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
