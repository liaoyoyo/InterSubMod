#!/usr/bin/env python3
"""Consolidate every computed result into one report-ready data file.

Every number quoted in the write-up must come from this file, so the report can
be assembled by injection rather than by hand (project rule 13-A).  Nothing here
recomputes science; it reads the per-step outputs, adds the per-chromosome
distribution view, and picks concrete worked examples for the "see the forest
and the trees" requirement.
"""

from __future__ import annotations

import json
import os
from collections import defaultdict

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OUT = os.path.join(DATA, "consolidated_findings.json")

ANNOT = os.path.join(DATA, "hcc1395_unit_cn_annotation.jsonl")


def load(name):
    with open(os.path.join(DATA, name)) as fh:
        return json.load(fh)


def main():
    cn_sum = load("cn_annotation_summary.json")
    af_exp = load("af_vs_cn_expectation.json")
    adj = load("savana_vs_seqc2_adjudication.json")
    refit = load("savana_refit_grid.json")
    purity = load("purity_selfconsistency_audit.json")
    res = load("resolution_vs_cn_stratified.json")
    mech1 = load("mechanism_and_intervals.json")
    mech2 = load("mechanism_v2_distinctness.json")

    # ---- per-chromosome distribution ----
    chrom_rows = []
    for chrom, classes in cn_sum["unit_class_by_chrom"].items():
        tot = sum(classes.values())
        neu = classes.get("neutral", 0)
        chrom_rows.append(
            {
                "chrom": chrom,
                "units": tot,
                "neutral": neu,
                "neutral_percent": round(100 * neu / tot, 2) if tot else None,
                "cn_altered_percent": round(100 * (tot - neu) / tot, 2) if tot else None,
                "gain": classes.get("gain", 0),
                "gain_loh": classes.get("gain_loh", 0),
                "cnloh": classes.get("cnloh", 0),
                "loss": classes.get("loss", 0) + classes.get("loss_loh", 0),
            }
        )

    def chrom_key(r):
        c = r["chrom"].replace("chr", "")
        return int(c) if c.isdigit() else 99

    chrom_rows.sort(key=chrom_key)

    with_units = [r for r in chrom_rows if r["units"] >= 50]
    cleanest = sorted(with_units, key=lambda r: -r["neutral_percent"])[:5]
    dirtiest = sorted(with_units, key=lambda r: r["neutral_percent"])[:5]
    biggest = sorted(chrom_rows, key=lambda r: -r["units"])[:5]

    # ---- concrete worked examples ----
    examples = {"neutral_unique": None, "gain_loh_af_broke_tie": None, "identical_af_tied": None}
    neutral_pool, gainloh_pool, ident_pool = [], [], []
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") != "ANNOTATED":
                continue
            if r.get("unit_status") != "ranked":
                continue
            try:
                ttc = int(r.get("total_tree_count") or 0)
            except (TypeError, ValueError):
                ttc = 0
            slim = {
                "region_id": r["region_id"],
                "chrom": r["chrom"],
                "span_start": r.get("span_start"),
                "span_end": r.get("span_end"),
                "hp_family": r.get("hp_family"),
                "active_bit_count": r.get("active_bit_count"),
                "seqc2_class": r["seqc2_class"],
                "seqc2_total_cn_weighted": r.get("seqc2_total_cn_weighted"),
                "seqc2_loh_fraction": r.get("seqc2_loh_fraction"),
                "resolution_class": r.get("resolution_class"),
                "total_tree_count": ttc,
                "best_tree_unique": r.get("best_tree_unique"),
            }
            if r["seqc2_class"] == "neutral" and ttc == 1:
                neutral_pool.append(slim)
            if (
                r["seqc2_class"] == "gain_loh"
                and ttc > 1
                and r.get("best_tree_unique")
            ):
                gainloh_pool.append(slim)

    if neutral_pool:
        examples["neutral_unique"] = neutral_pool[0]
    if gainloh_pool:
        gainloh_pool.sort(key=lambda x: -(x["seqc2_total_cn_weighted"] or 0))
        examples["gain_loh_af_broke_tie"] = gainloh_pool[0]

    # ---- headline synthesis ----
    contested_neu = res["marginal_test"]["neutral"]
    contested_alt = res["marginal_test"]["cn_altered"]

    synthesis = {
        "question_1_how_much_sits_on_cn_altered_ground": {
            "unit_level": cn_sum["unit_cn_altered"],
            "site_level": cn_sum["site_cn_altered"],
            "unit_class_percent": cn_sum["unit_class_percent"],
            "note": "denominator is canonical 2026-07-24 ranked+unranked located units for HCC1395; the 2026-07-10 figure of 94.2% used the superseded 23,810-sSNV backbone",
        },
        "question_2_is_the_topology_layer_affected": {
            "structure_unique_percent_neutral": res["headline_rates"]["cn_neutral_only"][
                "structure_unique_percent"
            ],
            "structure_unique_percent_altered": res["headline_rates"]["cn_altered_only"][
                "structure_unique_percent"
            ],
            "reading": "structural uniqueness does not rise on CN-altered ground (it is slightly lower), so the co-occurrence/topology layer shows no CN inflation",
        },
        "question_3_is_the_read_af_layer_affected": {
            "af_broke_tie_percent_neutral": contested_neu,
            "af_broke_tie_percent_altered": contested_alt,
            "fisher": {
                "odds_ratio": res["marginal_test"]["fisher_odds_ratio_altered_vs_neutral"],
                "p_value": res["marginal_test"]["fisher_p_value"],
            },
            "cmh_stratified_tree_count": res["stratified_by_candidate_tree_count"]["cmh"],
            "cmh_stratified_k": res["stratified_by_k"]["cmh"],
            "cmh_stratified_both": res["stratified_by_both"]["cmh"],
            "reading": "the read-AF tie-breaking layer is strongly CN-inflated and the effect survives stratification on candidate-tree count and k",
        },
        "question_4_mechanism": {
            "spread_hypothesis": {
                "claim": "CN widens within-unit AF spread",
                "result": mech1["mechanism"]["test_altered_vs_neutral"],
                "verdict": "initially REFUTED against a site set that wrongly included inactive all-reference rows; SUPPORTED once restricted to active sites",
            },
            "supported_hypothesis": {
                "claim": "CN determines whether the active read-AF values are arithmetically distinguishable at all; identical values force a tie",
                "prediction_1": mech2["prediction_1_identical_af_more_common_on_neutral"],
                "prediction_2": mech2["prediction_2_identical_af_units_stay_tied"],
                "mediation": mech2["mediation_check"]["cmh_controlling_distinctness_and_tree_count"],
                "non_degenerate_contrast": mech2["non_degenerate_contrast"],
                "arithmetically_unresolvable_share": mech2["arithmetically_unresolvable_share"],
                "verdict": "SUPPORTED and fully mediating: once arithmetic resolvability is controlled the CN effect is not detectable",
            },
            "correction_note": {
                "what_was_wrong": "the first run computed distinctness and spread over every af_coverage row, including inactive all-reference sites (fraction 0/1) that the score never sees",
                "found_by": "adversarial review, 2026-08-10",
                "effect": "reversed the mediation verdict (partial -> full) and reversed the spread test (refuted -> supported); headline CN association, excess estimate and SAVANA results were unaffected",
            },
            "cn_explainability_of_af": af_exp["cn_explainability"],
        },
        "question_5_how_much_of_the_conclusion_is_cn_attributable": mech1[
            "cn_attributable_uniqueness_excess"
        ],
        "question_6_savana_as_pipeline_validation": {
            "published_fit": refit["published_fit"],
            "grid_optimum": refit["best_by_state_agreement"],
            "seqc2_implied_ploidy": adj["seqc2_implied_ploidy_autosomal_bp_weighted"],
            "state_agreement_published": adj["state_agreement"],
            "regression_nonneutral": adj["regression_savana_on_seqc2_nonneutral_only"],
            "per_chrom_bias": adj["per_chrom_mean_bias_savana_minus_seqc2"],
            "verdict": "signal and segmentation are sound; the purity/ploidy fit is wrong; error is a calibration failure, which is correctable in principle but only where a truth set or a valid prior exists",
        },
        "question_7_other_samples": {
            "calibration": purity["calibration_case"],
            "per_sample": {
                s: {
                    "fitted_purity": v.get("fitted_purity"),
                    "fitted_ploidy": v.get("fitted_ploidy"),
                    "violation_rate_percent": v.get("violation_rate_percent"),
                    "baf_median": v.get("baf_median"),
                }
                for s, v in purity["samples"].items()
                if isinstance(v, dict) and "fitted_purity" in v
            },
        },
    }

    out = {
        "generated_by": os.path.basename(__file__),
        "sample": "HCC1395",
        "scope": "chr1-22",
        "canonical_authority": (
            "docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json"
        ),
        "synthesis": synthesis,
        "per_chromosome_distribution": chrom_rows,
        "chromosome_extremes": {
            "cleanest_by_neutral_percent": cleanest,
            "most_cn_altered": dirtiest,
            "largest_by_unit_count": biggest,
        },
        "worked_examples": examples,
        "source_files": {
            "cn_annotation_summary": "cn_annotation_summary.json",
            "af_vs_cn_expectation": "af_vs_cn_expectation.json",
            "savana_vs_seqc2_adjudication": "savana_vs_seqc2_adjudication.json",
            "savana_refit_grid": "savana_refit_grid.json",
            "purity_selfconsistency_audit": "purity_selfconsistency_audit.json",
            "resolution_vs_cn_stratified": "resolution_vs_cn_stratified.json",
            "mechanism_spread_refuted": "mechanism_and_intervals.json",
            "mechanism_distinctness_supported": "mechanism_v2_distinctness.json",
        },
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print("per-chromosome (top by unit count):")
    for r in biggest:
        print(
            f"  {r['chrom']:6s} units={r['units']:5d} neutral={r['neutral']:4d} "
            f"({r['neutral_percent']:5.2f}%) altered={r['cn_altered_percent']:5.2f}%"
        )
    print("\ncleanest chromosomes (highest neutral %):")
    for r in cleanest:
        print(
            f"  {r['chrom']:6s} units={r['units']:5d} neutral%={r['neutral_percent']:5.2f}"
        )
    print("\nmost CN-altered:")
    for r in dirtiest:
        print(
            f"  {r['chrom']:6s} units={r['units']:5d} neutral%={r['neutral_percent']:5.2f}"
        )
    print(f"\nworked examples: {json.dumps(examples, indent=1, ensure_ascii=False)}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
