#!/usr/bin/env python3
"""
HD-4: Attribution of the AF -> NGroups association.

Question: is the LOH intermediate-AF -> higher HPFineNGroups association
  (a) a methylation-biology signal, or
  (b) a phasing/haplotype-admixture artifact (NGroups = count of populated HP
      sub-families {HP1, HP1-1, HP2, HP2-1}, mechanically driven by allele balance)?

Strategy (all on landed master data, paired mode, LOH-bed TP, CN1 to match the
original H1p analysis):
  0. Reproduce baseline NGroups~AF (point-biserial via |AF-0.5|).
  1. Partial correlation NGroups~AF controlling for METHYLATION features.
  2. Partial correlation NGroups~AF controlling for PHASING-OCCUPANCY proxies
     (HP-family balance, assignment/unphased/conflict ratios).
  3. Decompose NGroups into its HP sub-family components to show what AF moves.

We use the SAME filters/thresholds as step2 (utils.classify_af, cn_tier, LOH.bed).
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import pandas as pd
from scipy import stats as st
from utils import (MASTER, classify_af, cn_tier, load_loh_beds, annotate_loh_bed)

DATADIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_subclone_af_paired/data")

EXTRA = ["NHP3", "NHP0", "hp0_ratio", "hp3_ratio", "hp_assign_rate"]
BASE = ["Chr", "Pos", "sample", "mode", "truth_label", "caller_af",
        "HPFineNGroups", "HPFineF", "HPFineP", "Coverage_Multiple", "NumReads",
        "HP_Ratio", "AlleleDelta", "PairwiseMeanDist", "CramersV",
        "ClusterPermanovaF", "HP1FamilyN", "HP2FamilyN"]


def partial_spearman(df, x, y, controls):
    """Spearman partial correlation of x,y controlling for `controls`.
    Method: rank-transform all, OLS-residualize x and y on controls, correlate residuals."""
    sub = df[[x, y] + controls].dropna()
    n = len(sub)
    if n < 30:
        return dict(n=n, r=np.nan, p=np.nan, note="insufficient n")
    R = sub.rank()
    C = np.column_stack([np.ones(n)] + [R[c].values for c in controls])
    def resid(v):
        beta, *_ = np.linalg.lstsq(C, v, rcond=None)
        return v - C @ beta
    rx, ry = resid(R[x].values), resid(R[y].values)
    if rx.std() == 0 or ry.std() == 0:
        return dict(n=n, r=np.nan, p=np.nan, note="zero-variance residual")
    r, p = st.pearsonr(rx, ry)  # pearson on residual ranks = partial spearman
    return dict(n=int(n), r=float(r), p=float(p))


def main():
    print("Loading master...")
    df = pd.read_csv(MASTER, sep="\t", usecols=BASE + EXTRA, low_memory=False)
    dfp = df[df["mode"] == "paired"].copy()
    print(f"  paired rows: {len(dfp):,}")

    beds = load_loh_beds()
    dfp["loh_bed_hit"] = annotate_loh_bed(dfp, beds)
    dfp["af_class"] = dfp["caller_af"].apply(classify_af)
    dfp["cn_tier"] = dfp["Coverage_Multiple"].apply(cn_tier)

    # Original H1p analysis cohort: LOH-bed TP, CN1
    coh = dfp[(dfp["loh_bed_hit"]) & (dfp["truth_label"] == "TP") &
              (dfp["cn_tier"] == "CN1")].copy()
    coh = coh.dropna(subset=["HPFineNGroups", "caller_af"])
    print(f"  CN1 LOH TP cohort: {len(coh):,}")

    # AF-centrality: distance from 0.5 (extreme=large, intermediate=small).
    # Use af_centrality = 0.5 - |AF-0.5|  so HIGHER = more intermediate (matches ΔNG>0).
    coh["af_dist_half"] = (coh["caller_af"] - 0.5).abs()
    coh["af_centrality"] = 0.5 - coh["af_dist_half"]   # 0=extreme, 0.5=exactly half

    # Phasing-admixture proxies (NO methylation):
    hp1 = coh["HP1FamilyN"].fillna(0)
    hp2 = coh["HP2FamilyN"].fillna(0)
    tot = hp1 + hp2
    coh["hp_balance"] = np.where(tot > 0, np.minimum(hp1, hp2) / np.maximum(hp1, hp2).replace(0, np.nan), 0)
    coh["hp_balance"] = coh["hp_balance"].fillna(0)
    coh["hp_assigned_total"] = tot
    coh["hp_minority_frac"] = np.where(tot > 0, np.minimum(hp1, hp2) / tot, 0)

    out = {"cohort": "paired, LOH-bed TP, CN1 (matches step2 H1p)",
           "n_cohort": int(len(coh)),
           "ngroups_definition": ("HPFineNGroups = LabelTest::hp_to_fine_labels -> "
               "count of distinct populated sub-families among {HP1,HP1-1,HP2,HP2-1} "
               "(max 4), thresholded by min_reads_per_group. Computed from BAM HP tag "
               "STRINGS only (ReadParser integer 1/2/11/21/33). Methylation distance "
               "matrix is used only for the PERMANOVA F-stat (HPFineF), NOT for the count. "
               "Source: src/core/LabelTest.cpp:265-305, :607-633."),
           "results": {}}

    # ---- 0. Baseline: NGroups ~ af_centrality (Spearman, no control) ----
    r0, p0 = st.spearmanr(coh["af_centrality"], coh["HPFineNGroups"])
    out["results"]["baseline_ngroups_vs_af_centrality"] = dict(
        spearman_r=float(r0), p=float(p0), n=int(len(coh)),
        note="positive r = more central(intermediate) AF -> higher NGroups (reproduces H1p)")
    print(f"\n0. baseline NGroups~af_centrality: r={r0:.4f} p={p0:.2e}")

    # ---- 1. Control for METHYLATION features ----
    meth_sets = {
        "HPFineF_only": ["HPFineF"],
        "AlleleDelta_only": ["AlleleDelta"],
        "PairwiseMeanDist_only": ["PairwiseMeanDist"],
        "ClusterPermanovaF_only": ["ClusterPermanovaF"],
        "all_methylation": ["HPFineF", "AlleleDelta", "PairwiseMeanDist",
                            "CramersV", "ClusterPermanovaF"],
    }
    out["results"]["control_for_methylation"] = {}
    for name, ctrl in meth_sets.items():
        res = partial_spearman(coh, "af_centrality", "HPFineNGroups", ctrl)
        out["results"]["control_for_methylation"][name] = dict(controls=ctrl, **res)
        print(f"1. NGroups~AF | {name}: r={res.get('r'):.4f} (n={res.get('n')})" if not np.isnan(res.get('r', np.nan)) else f"1. {name}: {res}")

    # ---- 2. Control for PHASING-OCCUPANCY proxies ----
    phas_sets = {
        "hp_balance_only": ["hp_balance"],
        "hp_minority_frac_only": ["hp_minority_frac"],
        "unphased_conflict": ["hp0_ratio", "hp3_ratio", "hp_assign_rate"],
        "all_phasing_occupancy": ["hp_balance", "hp_minority_frac",
                                  "hp_assigned_total", "hp0_ratio",
                                  "hp3_ratio", "hp_assign_rate"],
    }
    out["results"]["control_for_phasing_occupancy"] = {}
    for name, ctrl in phas_sets.items():
        res = partial_spearman(coh, "af_centrality", "HPFineNGroups", ctrl)
        out["results"]["control_for_phasing_occupancy"][name] = dict(controls=ctrl, **res)
        print(f"2. NGroups~AF | {name}: r={res.get('r'):.4f} (n={res.get('n')})" if not np.isnan(res.get('r', np.nan)) else f"2. {name}: {res}")

    # ---- 2b. Reverse framing: AF -> hp_balance (does AF drive admixture?) ----
    rb, pb = st.spearmanr(coh["af_centrality"], coh["hp_balance"])
    rbn, pbn = st.spearmanr(coh["hp_balance"], coh["HPFineNGroups"])
    out["results"]["mechanism_chain"] = dict(
        af_centrality_to_hp_balance=dict(spearman_r=float(rb), p=float(pb)),
        hp_balance_to_ngroups=dict(spearman_r=float(rbn), p=float(pbn)),
        note="if AF drives hp_balance AND hp_balance drives NGroups, the chain is phasing-mediated")
    print(f"\n2b. af_centrality->hp_balance: r={rb:.4f}; hp_balance->NGroups: r={rbn:.4f}")

    # ---- 3. What does NGroups actually count? pct with >=2 groups by AF class ----
    by_class = coh.groupby("af_class").agg(
        n=("HPFineNGroups", "size"),
        mean_ng=("HPFineNGroups", "mean"),
        pct_ge2=("HPFineNGroups", lambda s: 100*(s >= 2).mean()),
        pct_ge3=("HPFineNGroups", lambda s: 100*(s >= 3).mean()),
        mean_hp_balance=("hp_balance", "mean"),
        mean_hp_minority_frac=("hp_minority_frac", "mean"),
    ).reset_index()
    out["results"]["ngroups_decomposition_by_af_class"] = by_class.round(4).to_dict("records")
    print("\n3. By AF class:\n", by_class.round(3).to_string(index=False))

    DATADIR.mkdir(exist_ok=True)
    with open(DATADIR / "hd4_attribution.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {DATADIR/'hd4_attribution.json'}")


if __name__ == "__main__":
    main()
