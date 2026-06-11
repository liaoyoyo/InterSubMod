#!/usr/bin/env python3
"""obs18c — add COLO829 to make the n=7 cross-sample NG=2 composition table.

COLO829's KDE-corrected TO ISM was really run 2026-04-23 (run.log '[TO pilot] Completed')
but was never added to obs18 SAMPLES_TO, so the Grade-B Wilcoxon was n=6 (5 cell lines +
HCC1395_DORADO replicate). This script computes COLO829's NG=2 Inner/Outer composition with
the SAME obs18 logic and appends it to the cached n=6 composition TSV, producing the n=7 TSV
that obs18b_wilcoxon_ng2_gap.py then tests. No C++ rerun — pure aggregation on existing CSVs.

Output: data/obs18_NG2_composition_by_sample_n7.tsv  (n=6 cached rows + COLO829 rows)
"""
import pandas as pd
from pathlib import Path

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination")
N6_TSV = ROOT / "data/obs18_NG2_composition_by_sample.tsv"
N7_TSV = ROOT / "data/obs18_NG2_composition_by_sample_n7.tsv"

COLO829 = {
    "tp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod/intersubmod_tp/significance_summary.csv",
    "fp": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod/intersubmod_fp/significance_summary.csv",
}
BUCKET_COLS = ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]


def categorize_combo(h1, h1s, h2, h2s):
    if h1 and h1s and not h2 and not h2s:
        return "same_HP1 (HP1 + HP1-1)"
    if h2 and h2s and not h1 and not h1s:
        return "same_HP2 (HP2 + HP2-1)"
    if h1 and h2s and not h1s and not h2:
        return "cross_het (HP1 + HP2-1)"
    if h1s and h2 and not h1 and not h2s:
        return "cross_het_inv (HP1-1 + HP2)"
    return "other"


def analyze_sample(name, paths):
    frames = []
    for label, col in [("tp", 1), ("fp", 0)]:
        df = pd.read_csv(paths[label], low_memory=False,
                         usecols=["RegionID", "AlleleDelta", "HPFineNGroups",
                                  "Potential_LOH"] + BUCKET_COLS)
        df["tp_label"] = col
        df["sample"] = name
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    ng2 = df[df["HPFineNGroups"] == 2].copy()
    for c in BUCKET_COLS:
        ng2[c] = ng2[c].fillna(0).astype(int)
    ng2["combo"] = ng2.apply(lambda r: categorize_combo(
        r["HPFineN_HP1"] > 0, r["HPFineN_HP1S"] > 0,
        r["HPFineN_HP2"] > 0, r["HPFineN_HP2S"] > 0), axis=1)
    ng2["loh_side"] = ng2["Potential_LOH"].apply(
        lambda x: "Inner" if str(x).lower() in ("true", "1") else "Outer")
    ng2["AF_class"] = ng2["AlleleDelta"].clip(0, 1).apply(
        lambda af: "Extreme" if (af < 0.1 or af > 0.9) else "x")
    ext = ng2[ng2["AF_class"] == "Extreme"].copy()
    res = (ext.groupby(["sample", "loh_side", "combo"])
           .agg(n=("tp_label", "size"), n_tp=("tp_label", "sum")).reset_index())
    res["tp_rate"] = res["n_tp"] / res["n"]
    res["n_fp"] = res["n"] - res["n_tp"]
    return res


def main():
    n6 = pd.read_csv(N6_TSV, sep="\t")
    colo = analyze_sample("COLO829", COLO829)
    n7 = pd.concat([n6, colo[n6.columns]], ignore_index=True)
    n7.to_csv(N7_TSV, sep="\t", index=False)

    # report COLO829 gap
    def get(side, combo):
        r = colo[(colo.loh_side == side) & (colo.combo == combo)]
        return (float(r.tp_rate.iloc[0]), int(r.n.iloc[0])) if not r.empty else (None, 0)
    inner_rate, inner_n = get("Inner", "same_HP1 (HP1 + HP1-1)")
    outer_rate, outer_n = get("Outer", "cross_het (HP1 + HP2-1)")
    print(f"COLO829 Inner same_HP1 TP rate = {inner_rate:.4f} (n={inner_n})")
    print(f"COLO829 Outer cross_het TP rate = {outer_rate:.4f} (n={outer_n})")
    if inner_rate is not None and outer_rate is not None:
        print(f"COLO829 gap = {inner_rate - outer_rate:.4f}")
    print(f"n7 TSV written: {N7_TSV} ({len(n7)} rows, samples={sorted(n7['sample'].unique())})")


if __name__ == "__main__":
    main()
