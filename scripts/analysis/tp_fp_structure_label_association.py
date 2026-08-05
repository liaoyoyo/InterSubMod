#!/usr/bin/env python3
"""
TP vs FP × structure-label association characterization (HCC1395 paired_full, single-sample).

Goal (characterization, NOT a filter — that direction is NEGATIVE/DEAD):
  Describe how structure-label association patterns (verdict class, cluster×label
  alignment, HP-axis vs allele-axis PERMANOVA, Δβ magnitude, CramerV) distribute
  across TP vs FP somatic-SNV regions. Answer: "which structure is FP-common?",
  "does any pattern appear (near-)only in TP?".

Hard guardrails baked in:
  * The historical verdict panel explicitly reads the legacy four-state taxonomy;
    it answers a different question than the raw PERMANOVA axes -> kept as a SEPARATE
    dimension, never merged with raw-axis dimensions or current v2 classes.
  * Coverage is the #1 confound (FP has fewer CpGs). Every key categorical dimension
    is re-tested with (a) CMH stratified by NumCpGs x NumReads tertiles, and
    (b) a coverage+LOH matched TP subset. Both reported alongside the full set.
  * "Only in TP" claims are bounded by the rule-of-three (FP n=627) + Wilson CIs.
  * PERMDISP column is degenerate in this run (all dispersion-warn=false) -> any
    PERMANOVA-sig cell is location-OR-dispersion; flagged in output, not separated.

Inputs : intersubmod_{tp,fp}/significance_summary.csv (114 cols, 1 row = 1 somatic SNV)
Outputs: _assets/crosstabs.json, _assets/confound_control.json, _assets/run_meta.json
Run    : python3 scripts/analysis/tp_fp_structure_label_association.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency, ks_2samp
from scipy.spatial import cKDTree
from statsmodels.stats.proportion import proportion_confint
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import StratifiedTable

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    LEGACY_CLASSES,
    UNKNOWN_LEGACY_CLASS,
    SchemaContractError,
    select_legacy_view,
)

# ──────────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────────
RUN_DIR = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/"
           "paired_full/20260420_HCC1395_paired_full_full_2")
TP_CSV = os.path.join(RUN_DIR, "intersubmod_tp", "significance_summary.csv")
FP_CSV = os.path.join(RUN_DIR, "intersubmod_fp", "significance_summary.csv")

OUT_DIR = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/"
           "2026/06/20260619_tp_fp_structure_association_HCC1395")
ASSETS = os.path.join(OUT_DIR, "_assets")

BOOL_COLS = [
    "PassedGating", "ClusterPermanovaValid", "ClusterDispersionWarn",
    "HPMergedSig", "HPFineSig", "AlleleSig", "LabelHPPermanovaValid",
    "LabelHPDispersionWarn", "LabelAllelePermanovaValid",
    "LabelAlleleDispersionWarn", "Potential_LOH", "PerCpgASM_Valid",
    "Significant", "SuggestFilter", "Tumor_HP_Valid", "Normal_HP_Valid",
    "SampleASM_Sig", "HP_Residual_Sig",
    "WithinHP_CleanMultigroup",
]
SIG = 0.05
VERIFICATION_SELECTED = "_VerificationClass_Legacy_Selected"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Historical TP/FP structure-label analysis using an explicit legacy verification view."
    )
    parser.add_argument("--tp-csv", default=TP_CSV)
    parser.add_argument("--fp-csv", default=FP_CSV)
    parser.add_argument("--output-dir", default=OUT_DIR)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize historical unversioned VerificationClass as the legacy four-state input.",
    )
    return parser.parse_args()


# ──────────────────────────────────────────────────────────────────────────
# Loading / typing
# ──────────────────────────────────────────────────────────────────────────
def _to_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False}
    )


def load(tp_csv, fp_csv, allow_unversioned_v1=False):
    tp = pd.read_csv(tp_csv, low_memory=False)
    fp = pd.read_csv(fp_csv, low_memory=False)
    if list(tp.columns) != list(fp.columns):
        raise SchemaContractError("TP/FP column schemas differ")
    tp["source_scope"] = "tp"
    fp["source_scope"] = "fp"
    df = pd.concat([tp, fp], ignore_index=True)
    view = select_legacy_view(
        df,
        allow_unversioned_v1=allow_unversioned_v1,
        unversioned_unknown_policy="exclude",
    )
    df[VERIFICATION_SELECTED] = view.values.fillna(UNKNOWN_LEGACY_CLASS)
    for c in BOOL_COLS:
        if c not in df.columns:
            raise SchemaContractError(f"missing required boolean field: {c}")
        df[c] = _to_bool(df[c])
        if df[c].isna().any():
            raise SchemaContractError(f"invalid or missing boolean values in {c}")
    return df, view


# ──────────────────────────────────────────────────────────────────────────
# Derived structural dimensions (each independent; do not merge verdict & raw)
# ──────────────────────────────────────────────────────────────────────────
def derive(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()

    # Alignment state (PAIRED cluster x label axis; def from 20260617 inventory)
    def align_state(r):
        if not r["PassedGating"]:
            return "no_cluster"
        return "aligned" if r["GlobalP"] <= SIG else "misaligned"
    d["alignment_state"] = d.apply(align_state, axis=1)

    # Label-PERMANOVA HP-axis vs allele-axis (conditioned on validity)
    hp_sig = d["LabelHPPermanovaValid"] & (d["LabelHPPermanovaP"] <= SIG)
    al_sig = d["LabelAllelePermanovaValid"] & (d["LabelAllelePermanovaP"] <= SIG)
    d["label_axis"] = np.select(
        [hp_sig & al_sig, hp_sig & ~al_sig, ~hp_sig & al_sig],
        ["both", "hp_only", "allele_only"], default="neither")

    # Raw cluster PERMANOVA significance (only meaningful where valid)
    d["cluster_permanova"] = np.select(
        [~d["ClusterPermanovaValid"],
         d["ClusterPermanovaValid"] & (d["ClusterPermanovaP"] <= SIG)],
        ["invalid", "sig"], default="nonsig")

    # Exact canonical per-row within-HP state. The deprecated report stub is never parsed.
    d["within_hp"] = np.where(
        d["WithinHP_CleanMultigroup"], "multi", "single")

    # Magnitude bins (|Δβ|), CramerV bins
    def absbin(col, edges):
        v = d[col].abs()
        labels = [f"[{edges[i]},{edges[i+1]})" for i in range(len(edges) - 2)]
        labels.append(f">={edges[-2]}")
        return pd.cut(v, bins=edges, labels=labels, right=False,
                      include_lowest=True).astype(str)

    db_edges = [0, 0.05, 0.10, 0.20, 1.01]
    d["alleleDelta_bin"] = absbin("AlleleDelta", db_edges)
    d["hpMergedDelta_bin"] = absbin("HPMergedDelta", db_edges)
    d["hpResidualDelta_bin"] = absbin("HP_Residual_Delta", db_edges)  # somatic-controlled
    d["sampleASM_bin"] = absbin("SampleASM_Delta", db_edges)          # tumor-vs-normal read-set
    d["cramersV_bin"] = pd.cut(
        d["CramersV"], bins=[0, 0.1, 0.2, 0.3, 1.01],
        labels=["[0,.1)", "[.1,.2)", "[.2,.3)", ">=.3"],
        right=False, include_lowest=True).astype(str)
    return d


# ──────────────────────────────────────────────────────────────────────────
# Per-dimension cross-tab + cell stats
# ──────────────────────────────────────────────────────────────────────────
def wilson(k, n):
    if n == 0:
        return [None, None]
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return [round(lo, 6), round(hi, 6)]


def crosstab_dimension(d, col, order=None):
    """One-vs-rest Fisher per category + props + Wilson CI + rule-of-three."""
    n_tp = int((d.source_scope == "tp").sum())
    n_fp = int((d.source_scope == "fp").sum())
    cats = order or sorted(d[col].dropna().unique().tolist())
    rows, pvals = [], []
    for c in cats:
        in_c = d[col] == c
        a = int((in_c & (d.source_scope == "fp")).sum())   # FP in cat
        b = int((in_c & (d.source_scope == "tp")).sum())   # TP in cat
        cfp, dtp = n_fp - a, n_tp - b
        # one-vs-rest 2x2: rows = scope(fp/tp), cols = in_cat / not
        odds, p = fisher_exact([[a, cfp], [b, dtp]])
        pct_fp = a / n_fp if n_fp else None
        pct_tp = b / n_tp if n_tp else None
        enr = (pct_fp / pct_tp) if (pct_tp and pct_tp > 0) else None
        rule3 = (3.0 / n_fp) if a == 0 else None  # 95% upper bound on FP rate
        rows.append({
            "category": str(c),
            "n_tp": b, "n_fp": a,
            "pct_tp": round(pct_tp, 6) if pct_tp is not None else None,
            "pct_fp": round(pct_fp, 6) if pct_fp is not None else None,
            "enrichment_fp_over_tp": round(enr, 4) if enr is not None else None,
            "wilson95_pct_tp": wilson(b, n_tp),
            "wilson95_pct_fp": wilson(a, n_fp),
            "fisher_or_fp_in_cat": (round(odds, 4) if np.isfinite(odds) else None),
            "fisher_p": p,
            "rule_of_three_fp_upper": round(rule3, 6) if rule3 is not None else None,
        })
        pvals.append(p)
    # omnibus chi2 on full scope x category table
    ct = pd.crosstab(d.source_scope, d[col])
    try:
        chi2, chi_p, dofree, _ = chi2_contingency(ct.values)
        omni = {"chi2": round(float(chi2), 4), "p": float(chi_p), "dof": int(dofree)}
    except Exception as e:
        omni = {"error": str(e)}
    return {"n_tp": n_tp, "n_fp": n_fp, "cells": rows,
            "omnibus_chi2": omni, "_fisher_pvals": pvals}


# ──────────────────────────────────────────────────────────────────────────
# Coverage confound control: (a) CMH stratified, (b) matched subset
# ──────────────────────────────────────────────────────────────────────────
def tertile_label(s):
    q = s.quantile([1/3, 2/3]).values
    return pd.cut(s, bins=[-np.inf, q[0], q[1], np.inf],
                  labels=["lo", "mid", "hi"]).astype(str)


def cmh_for_category(d, col, cat):
    """CMH common OR + Breslow-Day for (FP-in-cat vs rest) across CpG x reads strata."""
    dd = d.copy()
    dd["cpg_t"] = tertile_label(dd["NumCpGs"])
    dd["rd_t"] = tertile_label(dd["NumReads"])
    dd["in_cat"] = (dd[col] == cat)
    tables = []
    for (_, _), g in dd.groupby(["cpg_t", "rd_t"]):
        a = int(((g.source_scope == "fp") & g.in_cat).sum())
        b = int(((g.source_scope == "fp") & ~g.in_cat).sum())
        c = int(((g.source_scope == "tp") & g.in_cat).sum())
        e = int(((g.source_scope == "tp") & ~g.in_cat).sum())
        if (a + b) > 0 and (c + e) > 0 and (a + c) > 0 and (b + e) > 0:
            tables.append([[a, b], [c, e]])
    if len(tables) < 2:
        return {"strata_used": len(tables), "note": "insufficient strata"}
    arr = np.array(tables).transpose(1, 2, 0)  # 2x2xK
    st = StratifiedTable(arr)
    cmh = st.test_null_odds()
    bd = st.test_equal_odds()
    return {
        "strata_used": len(tables),
        "cmh_common_oddsratio_fp": round(float(st.oddsratio_pooled), 4),
        "cmh_statistic": round(float(cmh.statistic), 4),
        "cmh_p": float(cmh.pvalue),
        "breslow_day_stat": round(float(bd.statistic), 4),
        "breslow_day_p": float(bd.pvalue),
    }


def matched_subset(d, k=5, seed=42):
    """Greedy NN match each FP to k TP on (NumReads, NumCpGs) within LOH stratum."""
    rng = np.random.default_rng(seed)
    keep_idx = []
    feats = ["NumReads", "NumCpGs"]
    mu = d[feats].mean().values
    sd = d[feats].std().values
    sd[sd == 0] = 1.0
    for loh_val in [True, False]:
        sub = d[d.Potential_LOH == loh_val]
        tp = sub[sub.source_scope == "tp"]
        fp = sub[sub.source_scope == "fp"]
        if len(fp) == 0 or len(tp) == 0:
            continue
        tp_z = ((tp[feats].values - mu) / sd)
        fp_z = ((fp[feats].values - mu) / sd)
        tree = cKDTree(tp_z)
        used = set()
        for i in range(len(fp_z)):
            # query enough neighbors to fill k unused
            dist, nb = tree.query(fp_z[i], k=min(k * 6, len(tp_z)))
            nb = np.atleast_1d(nb)
            picked = 0
            for j in nb:
                if j not in used:
                    used.add(j)
                    picked += 1
                    if picked >= k:
                        break
        keep_idx.extend(tp.index[list(used)].tolist())
        keep_idx.extend(fp.index.tolist())
    matched = d.loc[sorted(set(keep_idx))].copy()
    # balance check: KS on NumReads/NumCpGs between FP and matched-TP
    m_tp = matched[matched.source_scope == "tp"]
    m_fp = matched[matched.source_scope == "fp"]
    bal = {}
    for f in feats:
        ks = ks_2samp(m_fp[f].values, m_tp[f].values)
        bal[f] = {"ks_stat": round(float(ks.statistic), 4), "ks_p": float(ks.pvalue),
                  "median_fp": float(m_fp[f].median()), "median_tp": float(m_tp[f].median())}
    return matched, {"n_tp_matched": int(len(m_tp)), "n_fp": int(len(m_fp)),
                     "k": k, "balance_ks": bal}


# ──────────────────────────────────────────────────────────────────────────
def main():
    global ASSETS
    args = parse_args()
    try:
        df, verification_view = load(
            args.tp_csv,
            args.fp_csv,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
    except SchemaContractError as error:
        raise SystemExit(f"FATAL: verification schema contract failed: {error}") from error
    ASSETS = os.path.join(args.output_dir, "_assets")
    os.makedirs(ASSETS, exist_ok=True)
    d = derive(df)

    n_tp = int((d.source_scope == "tp").sum())
    n_fp = int((d.source_scope == "fp").sum())

    # degeneracy / coverage sanity (recorded as ground-truth context)
    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tp_csv": str(Path(args.tp_csv).resolve()),
        "fp_csv": str(Path(args.fp_csv).resolve()),
        "verification_selection": verification_view.metadata(),
        "verification_selection_values": list(LEGACY_CLASSES),
        "n_tp": n_tp, "n_fp": n_fp,
        "permdisp_degenerate": {
            "ClusterDispersionWarn_true_count": int(d.ClusterDispersionWarn.sum()),
            "ClusterDispersionP_lt05_count": int((d.ClusterDispersionP < SIG).sum()),
            "note": "PERMDISP carries no signal in this run; PERMANOVA-sig = location-OR-dispersion",
        },
        "cluster_permanova_valid": {
            "tp_true": int(((d.source_scope == "tp") & d.ClusterPermanovaValid).sum()),
            "fp_true": int(((d.source_scope == "fp") & d.ClusterPermanovaValid).sum()),
        },
        "coverage_medians": {
            "NumReads_tp": float(d[d.source_scope == "tp"].NumReads.median()),
            "NumReads_fp": float(d[d.source_scope == "fp"].NumReads.median()),
            "NumCpGs_tp": float(d[d.source_scope == "tp"].NumCpGs.median()),
            "NumCpGs_fp": float(d[d.source_scope == "fp"].NumCpGs.median()),
        },
        "potential_loh_rate": {
            "tp": float(d[d.source_scope == "tp"].Potential_LOH.mean()),
            "fp": float(d[d.source_scope == "fp"].Potential_LOH.mean()),
        },
    }

    # ── Dimensions ──
    dims = {
        "legacy_VerificationClass": (
            VERIFICATION_SELECTED,
            list(LEGACY_CLASSES) + [UNKNOWN_LEGACY_CLASS],
        ),
        "alignment_state": ("alignment_state",
                            ["aligned", "misaligned", "no_cluster"]),
        "label_axis_hp_vs_allele": ("label_axis",
                                    ["both", "hp_only", "allele_only", "neither"]),
        "cluster_permanova": ("cluster_permanova", ["sig", "nonsig", "invalid"]),
        "cramersV_bin": ("cramersV_bin", None),
        "alleleDelta_bin": ("alleleDelta_bin", None),
        "hpMergedDelta_bin": ("hpMergedDelta_bin", None),
        "hpResidualDelta_bin_somatic_controlled": ("hpResidualDelta_bin", None),
        "sampleASM_bin_tumor_vs_normal": ("sampleASM_bin", None),
        "within_hp_multigroup_proxy": ("within_hp", ["multi", "single"]),
    }
    crosstabs, all_p, p_index = {}, [], []
    for name, (col, order) in dims.items():
        res = crosstab_dimension(d, col, order)
        for i, _ in enumerate(res.pop("_fisher_pvals")):
            p_index.append((name, i))
        crosstabs[name] = res

    # collect p-values for family-wide BH-FDR (one-vs-rest Fisher across ALL cells)
    for name in crosstabs:
        for cell in crosstabs[name]["cells"]:
            all_p.append(cell["fisher_p"])
    rej, q, _, _ = multipletests(all_p, alpha=0.05, method="fdr_bh")
    qi = 0
    for name in crosstabs:
        for cell in crosstabs[name]["cells"]:
            cell["fisher_q_bh"] = float(q[qi])
            cell["fdr_sig_q05"] = bool(rej[qi])
            qi += 1
    meta["fdr_family_size"] = len(all_p)

    # ── Confound control on key categorical dims ──
    cmh = {}
    for dim_col, cats in [
        (VERIFICATION_SELECTED, list(LEGACY_CLASSES) + [UNKNOWN_LEGACY_CLASS]),
        ("alignment_state", ["aligned", "misaligned", "no_cluster"]),
        ("label_axis", ["both", "hp_only", "allele_only"]),
    ]:
        cmh[dim_col] = {c: cmh_for_category(d, dim_col, c) for c in cats}

    matched, match_meta = matched_subset(d, k=5)
    matched_crosstabs = {}
    for name, col, order in [
        ("legacy_VerificationClass", VERIFICATION_SELECTED,
         list(LEGACY_CLASSES) + [UNKNOWN_LEGACY_CLASS]),
        ("alignment_state", "alignment_state",
         ["aligned", "misaligned", "no_cluster"]),
        ("label_axis_hp_vs_allele", "label_axis",
         ["both", "hp_only", "allele_only", "neither"]),
    ]:
        r = crosstab_dimension(matched, col, order)
        r.pop("_fisher_pvals", None)
        matched_crosstabs[name] = r

    # ── Write ──
    with open(os.path.join(ASSETS, "run_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(ASSETS, "crosstabs.json"), "w") as f:
        json.dump(crosstabs, f, indent=2)
    with open(os.path.join(ASSETS, "confound_control.json"), "w") as f:
        json.dump({"cmh": cmh, "matched_subset": {"meta": match_meta,
                   "crosstabs": matched_crosstabs}}, f, indent=2)

    print("OK  n_tp=%d n_fp=%d  dims=%d  fdr_family=%d"
          % (n_tp, n_fp, len(crosstabs), len(all_p)))
    print("matched: TP=%d FP=%d  CpG KS p=%.3g  reads KS p=%.3g"
          % (match_meta["n_tp_matched"], match_meta["n_fp"],
             match_meta["balance_ks"]["NumCpGs"]["ks_p"],
             match_meta["balance_ks"]["NumReads"]["ks_p"]))


if __name__ == "__main__":
    main()
