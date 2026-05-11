#!/usr/bin/env python3
"""
Step 4b: TP/FP discrimination on TO mode (HCC1395 only, 40k rows).

Rationale: paired_full has TP 95.9% which leaves minimal FP for discrimination
analysis. TO mode has TP 71.1% -- that's the real test bed. This step mirrors
step4 schemes on HCC1395 TO mode.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import DATA_DIR, CN_STRATEGIES, assign_cn_tier
from step4_tpfp_discrimination import (
    wilson_ci, tp_enrichment, compute_cube, add_enrichment,
)


MASTER_TSV = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
CN_F_LABEL = {
    "T0": "CN<1 (del/LOH)",
    "T1": "CN~2 (diploid)",
    "T2": "CN~3 (gain)",
    "T3": "CN~4",
    "T4": "CN>=5 (amp)",
}


def main() -> None:
    df = pd.read_csv(MASTER_TSV, sep="\t", compression="gzip", low_memory=False)
    to = df[df["mode"] == "to_pileup"].copy()
    print(f"[step4b] TO rows (HCC1395): {len(to)}")
    print(f"[step4b] TP rate: {to['tp_label'].mean():.4f} (TP={to['tp_label'].sum()}, FP={(1-to['tp_label']).sum()})")

    to["cn_tier_F"] = assign_cn_tier(to["CovM_used"], CN_STRATEGIES["F_SEQC2_grounded"])
    to["cn_tier_F_label"] = to["cn_tier_F"].astype(str).map(CN_F_LABEL)
    to["LOH_Subtype"] = to["LOH_Subtype"].fillna("None")

    totals = pd.DataFrame([dict(
        sample="HCC1395_TO",
        N=len(to), N_tp=int(to["tp_label"].sum()),
        N_fp=int((to["tp_label"]==0).sum()),
        overall_tp_rate=to["tp_label"].mean(),
    )])
    print("\n[step4b] totals:")
    print(totals.to_string(index=False))

    to["sample"] = "HCC1395_TO"
    cube_coarse = compute_cube(to, ["sample","LOH_Subtype","cn_tier_F","AF_class"])
    cube_coarse = add_enrichment(cube_coarse, totals[["sample","N_tp","N_fp","overall_tp_rate"]])
    cube_coarse["cn_tier_F_label"] = cube_coarse["cn_tier_F"].astype(str).map(CN_F_LABEL)
    cube_coarse.to_csv(DATA_DIR / "tpfp_cube_coarse_TO.tsv", sep="\t", index=False)
    print(f"[step4b] wrote tpfp_cube_coarse_TO.tsv: {cube_coarse.shape}")

    cube_fine = compute_cube(to, ["sample","LOH_Subtype","cn_tier_F","AF_bin10"])
    cube_fine = add_enrichment(cube_fine, totals[["sample","N_tp","N_fp","overall_tp_rate"]])
    cube_fine.to_csv(DATA_DIR / "tpfp_cube_fine_TO.tsv", sep="\t", index=False)
    print(f"[step4b] wrote tpfp_cube_fine_TO.tsv: {cube_fine.shape}")

    # Top TP cells
    MIN_N = 200
    high_tp = cube_coarse[(cube_coarse["n"] >= MIN_N)].copy()
    high_tp["lift"] = high_tp["tp_rate"] - high_tp["overall_tp_rate"]
    high_tp = high_tp.sort_values("tp_rate", ascending=False)
    print(f"\n[step4b] Top 15 TP-enriched cells (n>={MIN_N}):")
    print(high_tp[["LOH_Subtype","cn_tier_F_label","AF_class","n","n_tp","n_fp","tp_rate","lift","tp_enrichment"]].head(15).to_string(index=False))

    high_fp = cube_coarse[(cube_coarse["n"] >= MIN_N)].copy()
    high_fp["fp_rate"] = high_fp["n_fp"] / high_fp["n"]
    high_fp = high_fp.sort_values("fp_rate", ascending=False)
    print(f"\n[step4b] Top 15 FP-enriched cells (n>={MIN_N}):")
    print(high_fp[["LOH_Subtype","cn_tier_F_label","AF_class","n","n_tp","n_fp","fp_rate","tp_enrichment"]].head(15).to_string(index=False))

    high_tp.to_csv(DATA_DIR / "tpfp_top_high_tp_cells_TO.tsv", sep="\t", index=False)
    high_fp.to_csv(DATA_DIR / "tpfp_top_high_fp_cells_TO.tsv", sep="\t", index=False)

    # Filter schemes on TO mode
    def metrics(mask, name, rationale):
        total_tp = int(to["tp_label"].sum())
        total_fp = int((to["tp_label"]==0).sum())
        sub = to[mask]
        n = len(sub); n_tp = int(sub["tp_label"].sum()); n_fp = n - n_tp
        return dict(
            scheme=name, rationale=rationale,
            n_included=n, frac=n/len(to),
            tp_rate=n_tp/n if n>0 else 0.0,
            tp_recall=n_tp/total_tp if total_tp>0 else 0.0,
            fp_recall=n_fp/total_fp if total_fp>0 else 0.0,
            fp_reduction=1-(n_fp/total_fp) if total_fp>0 else 0.0,
            n_tp=n_tp, n_fp=n_fp,
        )

    loh = to["LOH_Subtype"]; af = to["AF_class"]; cn = to["cn_tier_F"].astype(str); ng = to["HPFineNGroups"]
    m1 = (loh=="LOH_Strong") & (af=="Extreme")
    m2 = (loh.isin(["LOH_Subclone","LOH_Weak"])) & (af=="Intermediate")
    m3 = (loh=="None") & (af=="Near-half") & (cn.isin(["T1","T2"]))
    m4 = (loh=="None") & (af=="Extreme")

    schemes = [
        metrics(m1, "S1_LOH_Strong_Extreme", "Pure-haplotype LOH; expect highest TP."),
        metrics(m2, "S2_Subclonal_LOH_Intermediate", "Subclonal/weak LOH intermediate AF."),
        metrics(m3, "S3_Diploid_Het", "No LOH + near-half AF + CN~2-3; canonical het somatic."),
        metrics(m4, "S4_NonLOH_Extreme", "No LOH + extreme AF; germline leak risk."),
        metrics((m1|m2|m3) & ~m4, "S5_Combo", "Union S1+S2+S3, excluding S4."),
        metrics(m1 & (ng>=3), "S6_S1_NG>=3", "S1 + NG>=3 subclone marker."),
        metrics((m1|m2|m3) & (ng>=3) & ~m4, "S7_Combo_NG>=3", "S5 + NG>=3."),
        metrics((m1|m3) & (ng==4) & (to["NumReads"]>=80), "S8_Pure_NG4_NR80", "Strictest: S1orS3 + NG=4 + NR>=80."),
    ]
    schemes_df = pd.DataFrame(schemes)
    schemes_df.to_csv(DATA_DIR / "tpfp_stratified_filter_schemes_TO.tsv", sep="\t", index=False)
    print("\n[step4b] TO Filter scheme metrics:")
    print(schemes_df[["scheme","n_included","tp_rate","tp_recall","fp_recall","fp_reduction","n_tp","n_fp"]].to_string(index=False))


if __name__ == "__main__":
    main()
