#!/usr/bin/env python3
"""
Step 4: TP/FP discrimination re-analysis (biology-module framing).

Pivot from v1 report which focused on NG effect size (off-topic).
This step answers:
  - Under new-KDE CN tiering, do LOH x AF combinations separate TP from FP?
  - Which cells (biology modules) yield high TP_rate with low FP enrichment?
  - What stratified filter schemes can be proposed per biology module?

Biology modules (from LOH_Subtype x AF_class):
  None + Extreme         -> germline-like / extreme somatic
  None + Near-half       -> diploid heterozygous somatic
  None + Intermediate    -> subclonal somatic in balanced CN
  LOH_Weak + Extreme     -> weak LOH with haplotype bias
  LOH_Strong + Extreme   -> deletion/cnLOH pure haplotype
  LOH_Strong + Intermediate -> residual admixture in LOH
  LOH_Subclone + *       -> subclonal LOH (mixed populations)
  LOH_Noise              -> ambiguous LOH annotation
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import DATA_DIR, OBS_DIR, SAMPLE_ORDER, CN_STRATEGIES, assign_cn_tier


MASTER_TSV = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"


# ----------------------------------------------------------------------------
def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    if n == 0:
        return (0.0, 0.0, 1.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z / denom) * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return p, max(0.0, center - half), min(1.0, center + half)


def tp_enrichment(n_tp: int, n_fp: int, N_tp: int, N_fp: int) -> float:
    """(n_tp/N_tp) / (n_fp/N_fp); returns +inf if no FP, 0 if no TP."""
    if n_fp == 0:
        return float("inf") if n_tp > 0 else 1.0
    if N_tp == 0 or N_fp == 0:
        return 1.0
    pt = n_tp / N_tp
    pf = n_fp / N_fp
    return pt / pf if pf > 0 else float("inf")


# ----------------------------------------------------------------------------
def compute_cube(df: pd.DataFrame, dims: list[str]) -> pd.DataFrame:
    """Compute TP/FP counts + rates along given dimensions."""
    g = df.groupby(dims, observed=True, dropna=False).agg(
        n=("tp_label", "size"),
        n_tp=("tp_label", "sum"),
    ).reset_index()
    g["n_fp"] = g["n"] - g["n_tp"]
    g["tp_rate"] = g["n_tp"] / g["n"]
    return g


def add_enrichment(cube: pd.DataFrame, per_sample_totals: pd.DataFrame) -> pd.DataFrame:
    """Join per-sample TP/FP totals and compute enrichment + Wilson CI."""
    out = cube.merge(per_sample_totals, on="sample", how="left")
    out["tp_enrichment"] = out.apply(
        lambda r: tp_enrichment(int(r["n_tp"]), int(r["n_fp"]), int(r["N_tp"]), int(r["N_fp"])),
        axis=1,
    )
    ci = out.apply(lambda r: wilson_ci(int(r["n_tp"]), int(r["n"])), axis=1)
    out["tp_rate_lo"] = [c[1] for c in ci]
    out["tp_rate_hi"] = [c[2] for c in ci]
    out["tp_rate_half_width"] = (out["tp_rate_hi"] - out["tp_rate_lo"]) / 2.0
    return out


# ----------------------------------------------------------------------------
def main() -> None:
    print(f"[step4] loading master: {MASTER_TSV}")
    df = pd.read_csv(MASTER_TSV, sep="\t", compression="gzip", low_memory=False)
    print(f"[step4] shape: {df.shape}")

    # Restrict to paired_full primary frame (use to_pileup separately for HCC1395)
    primary = df[df["mode"] == "paired_full"].copy()
    print(f"[step4] paired_full rows: {len(primary)}")

    # Assign CN tier F (SEQC2-grounded empirical midpoints 0.65/0.99/1.33/1.82)
    primary["cn_tier_F"] = assign_cn_tier(primary["CovM_used"], CN_STRATEGIES["F_SEQC2_grounded"])
    # Human label for CN tier F
    CN_F_LABEL = {
        "T0": "CN<1 (del/LOH)",
        "T1": "CN~2 (diploid)",
        "T2": "CN~3 (gain)",
        "T3": "CN~4",
        "T4": "CN>=5 (amp)",
    }
    primary["cn_tier_F_label"] = primary["cn_tier_F"].astype(str).map(CN_F_LABEL)

    # Ensure categorical ordering
    primary["LOH_Subtype"] = primary["LOH_Subtype"].fillna("None")

    # Per-sample totals
    totals = primary.groupby("sample").agg(
        N=("tp_label", "size"),
        N_tp=("tp_label", "sum"),
    ).reset_index()
    totals["N_fp"] = totals["N"] - totals["N_tp"]
    totals["overall_tp_rate"] = totals["N_tp"] / totals["N"]
    print("\n[step4] per-sample totals:")
    print(totals.to_string(index=False))

    # ----------------------------------------------------------------------
    # A1: 3D cube per-sample: LOH_Subtype x CN_F x AF_class x tp_label
    # (use AF_class for coarse, AF_bin10 for fine)
    # ----------------------------------------------------------------------
    cube_coarse = compute_cube(primary, ["sample", "LOH_Subtype", "cn_tier_F", "AF_class"])
    cube_coarse = add_enrichment(cube_coarse, totals[["sample", "N_tp", "N_fp", "overall_tp_rate"]])
    # Add human label
    cube_coarse["cn_tier_F_label"] = cube_coarse["cn_tier_F"].astype(str).map(CN_F_LABEL)
    cube_coarse.to_csv(DATA_DIR / "tpfp_cube_coarse.tsv", sep="\t", index=False)
    print(f"[step4] wrote tpfp_cube_coarse.tsv: {cube_coarse.shape}")

    cube_fine = compute_cube(primary, ["sample", "LOH_Subtype", "cn_tier_F", "AF_bin10"])
    cube_fine = add_enrichment(cube_fine, totals[["sample", "N_tp", "N_fp", "overall_tp_rate"]])
    cube_fine.to_csv(DATA_DIR / "tpfp_cube_fine.tsv", sep="\t", index=False)
    print(f"[step4] wrote tpfp_cube_fine.tsv: {cube_fine.shape}")

    # ----------------------------------------------------------------------
    # A2: Top cells with high TP_rate AND high FP enrichment (opposites)
    # ----------------------------------------------------------------------
    MIN_N = 100
    # High-TP cells (high TP rate vs baseline + enough n)
    high_tp = cube_coarse[(cube_coarse["n"] >= MIN_N)].copy()
    high_tp["lift"] = high_tp["tp_rate"] - high_tp["overall_tp_rate"]
    high_tp = high_tp.sort_values("lift", ascending=False)
    high_tp.to_csv(DATA_DIR / "tpfp_top_high_tp_cells.tsv", sep="\t", index=False)
    print(f"[step4] top 10 TP-enriched cells (n>={MIN_N}):")
    print(high_tp[["sample","LOH_Subtype","cn_tier_F_label","AF_class","n","n_tp","n_fp","tp_rate","tp_enrichment","lift"]].head(10).to_string(index=False))

    # High-FP cells
    high_fp = cube_coarse[(cube_coarse["n"] >= MIN_N) & (cube_coarse["n_fp"] >= 10)].copy()
    high_fp["fp_rate"] = high_fp["n_fp"] / high_fp["n"]
    high_fp["fp_enrichment"] = high_fp.apply(
        lambda r: (r["n_fp"]/r["N_fp"]) / (r["n_tp"]/r["N_tp"]) if (r["n_tp"]>0 and r["N_fp"]>0) else float("inf"),
        axis=1,
    )
    high_fp = high_fp.sort_values("fp_enrichment", ascending=False)
    high_fp.to_csv(DATA_DIR / "tpfp_top_high_fp_cells.tsv", sep="\t", index=False)
    print(f"\n[step4] top 10 FP-enriched cells (n>={MIN_N}, n_fp>=10):")
    print(high_fp[["sample","LOH_Subtype","cn_tier_F_label","AF_class","n","n_tp","n_fp","fp_rate","fp_enrichment"]].head(10).to_string(index=False))

    # ----------------------------------------------------------------------
    # A3: NG effect within biology modules (does NG>=3 lift TP rate?)
    # Strata: LOH_Subtype x AF_class x sample, NG levels: all / NG>=3 / NG=4
    # ----------------------------------------------------------------------
    ng_rows = []
    for (samp, loh_sub, af_cls), g in primary.groupby(["sample", "LOH_Subtype", "AF_class"], observed=True):
        if len(g) < 50:
            continue
        for ng_filter, subset in [
            ("all", g),
            ("NG>=3", g[g["HPFineNGroups"] >= 3]),
            ("NG=4", g[g["HPFineNGroups"] == 4]),
            ("NG=4_NR>=80", g[(g["HPFineNGroups"] == 4) & (g["NumReads"] >= 80)]),
        ]:
            n = len(subset); n_tp = int(subset["tp_label"].sum())
            p, lo, hi = wilson_ci(n_tp, n)
            ng_rows.append(dict(
                sample=samp, LOH_Subtype=loh_sub, AF_class=af_cls,
                ng_filter=ng_filter, n=n, n_tp=n_tp, n_fp=n-n_tp,
                tp_rate=p, tp_rate_lo=lo, tp_rate_hi=hi,
            ))
    ng_module = pd.DataFrame(ng_rows)
    ng_module.to_csv(DATA_DIR / "tpfp_ng_within_biology_modules.tsv", sep="\t", index=False)
    print(f"\n[step4] wrote tpfp_ng_within_biology_modules.tsv: {ng_module.shape}")

    # ----------------------------------------------------------------------
    # A4: Stratified filter scheme proposals (4 biology-informed modules)
    # For each scheme, compute across all samples:
    #   - n included, TP rate, TP recall (vs all TP), FP reduction
    # ----------------------------------------------------------------------
    def scheme_metrics(primary_df: pd.DataFrame, mask: pd.Series, name: str, rationale: str) -> dict:
        total_tp = int(primary_df["tp_label"].sum())
        total_fp = int((primary_df["tp_label"] == 0).sum())
        sub = primary_df[mask]
        n = len(sub)
        n_tp = int(sub["tp_label"].sum())
        n_fp = n - n_tp
        return dict(
            scheme=name,
            rationale=rationale,
            n_included=n,
            frac_of_all_regions=n / len(primary_df),
            tp_rate=n_tp / n if n > 0 else 0.0,
            tp_recall=n_tp / total_tp if total_tp > 0 else 0.0,
            fp_recall=n_fp / total_fp if total_fp > 0 else 0.0,
            fp_reduction=1 - (n_fp / total_fp) if total_fp > 0 else 0.0,
            n_tp=n_tp, n_fp=n_fp,
        )

    schemes = []

    # S1: Pure-haplotype LOH module (strong signal, high TP expected)
    m1 = (primary["LOH_Subtype"] == "LOH_Strong") & (primary["AF_class"] == "Extreme")
    schemes.append(scheme_metrics(primary, m1, "S1_LOH_Strong_Extreme_AF",
        "Deletion/cnLOH + pure-haplotype AF; expect highest TP purity (single surviving allele)."))

    # S2: Subclonal LOH module (intermediate AF; admixture signal)
    m2 = (primary["LOH_Subtype"].isin(["LOH_Subclone", "LOH_Weak"])) & (primary["AF_class"] == "Intermediate")
    schemes.append(scheme_metrics(primary, m2, "S2_Subclonal_LOH_Intermediate_AF",
        "Subclonal/weak LOH with intermediate AF; admixed tumor populations, NG can further resolve."))

    # S3: Diploid heterozygous somatic
    m3 = (primary["LOH_Subtype"] == "None") & (primary["AF_class"] == "Near-half") & (primary["cn_tier_F"].astype(str).isin(["T1","T2"]))
    schemes.append(scheme_metrics(primary, m3, "S3_Diploid_Het_Somatic",
        "No LOH + near-half AF + CN~2-3; canonical heterozygous somatic signal."))

    # S4: High-risk germline-leak (extreme AF, non-LOH) -> expected high FP
    m4 = (primary["LOH_Subtype"] == "None") & (primary["AF_class"] == "Extreme")
    schemes.append(scheme_metrics(primary, m4, "S4_NonLOH_Extreme_AF_HIGH_FP_RISK",
        "No LOH + extreme AF; suspect germline leak or mapping artifact; expected high FP."))

    # S5: Filter combo = S1 OR S2 OR S3 (exclude S4 high-risk)
    m5 = (m1 | m2 | m3) & ~m4
    schemes.append(scheme_metrics(primary, m5, "S5_Combined_S1_S2_S3_EXCLUDE_S4",
        "Union of high-confidence modules; excludes S4 high-FP-risk."))

    # S6: S1 + NG>=3 subclone confirmation
    m6 = m1 & (primary["HPFineNGroups"] >= 3)
    schemes.append(scheme_metrics(primary, m6, "S6_S1_plus_NG>=3",
        "S1 + NG>=3 subclone marker; requires both LOH biology + read haplotype structure."))

    schemes_df = pd.DataFrame(schemes)
    schemes_df.to_csv(DATA_DIR / "tpfp_stratified_filter_schemes.tsv", sep="\t", index=False)
    print("\n[step4] Filter scheme metrics:")
    print(schemes_df[["scheme","n_included","tp_rate","tp_recall","fp_recall","fp_reduction","n_tp","n_fp"]].to_string(index=False))

    # ----------------------------------------------------------------------
    # A5: Per-sample module-wise metrics (for small-multiple heatmaps)
    # ----------------------------------------------------------------------
    def per_sample_scheme(samp_df: pd.DataFrame, mask: pd.Series) -> dict:
        total_tp = int(samp_df["tp_label"].sum())
        total_fp = int((samp_df["tp_label"] == 0).sum())
        sub = samp_df[mask]
        n = len(sub); n_tp = int(sub["tp_label"].sum())
        return dict(
            n=n, tp_rate=(n_tp / n if n > 0 else 0.0),
            tp_recall=(n_tp / total_tp if total_tp > 0 else 0.0),
            fp_recall=((n - n_tp) / total_fp if total_fp > 0 else 0.0),
        )

    rows_ps = []
    for samp in SAMPLE_ORDER:
        sdf = primary[primary["sample"] == samp]
        if sdf.empty:
            continue
        loh = sdf["LOH_Subtype"]; af = sdf["AF_class"]; cn = sdf["cn_tier_F"].astype(str)
        ng = sdf["HPFineNGroups"]
        for scheme_name, mask in [
            ("S1_LOH_Strong_Extreme", (loh=="LOH_Strong") & (af=="Extreme")),
            ("S2_Subclonal_LOH_Inter", (loh.isin(["LOH_Subclone","LOH_Weak"])) & (af=="Intermediate")),
            ("S3_Diploid_Het", (loh=="None") & (af=="Near-half") & (cn.isin(["T1","T2"]))),
            ("S4_NonLOH_Extreme", (loh=="None") & (af=="Extreme")),
            ("S5_Combo_S1_S2_S3", (((loh=="LOH_Strong") & (af=="Extreme")) |
                                    ((loh.isin(["LOH_Subclone","LOH_Weak"])) & (af=="Intermediate")) |
                                    ((loh=="None") & (af=="Near-half") & (cn.isin(["T1","T2"])))) &
                                   ~((loh=="None") & (af=="Extreme"))),
            ("S6_S1_NG>=3", (loh=="LOH_Strong") & (af=="Extreme") & (ng>=3)),
        ]:
            metrics = per_sample_scheme(sdf, mask)
            rows_ps.append(dict(sample=samp, scheme=scheme_name, **metrics))
    per_sample_schemes = pd.DataFrame(rows_ps)
    per_sample_schemes.to_csv(DATA_DIR / "tpfp_per_sample_schemes.tsv", sep="\t", index=False)
    print(f"\n[step4] wrote tpfp_per_sample_schemes.tsv: {per_sample_schemes.shape}")

    print("\n[step4] DONE")


if __name__ == "__main__":
    main()
