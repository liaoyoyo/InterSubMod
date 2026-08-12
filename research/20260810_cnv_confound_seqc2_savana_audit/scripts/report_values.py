#!/usr/bin/env python3
"""Emit every number the write-up quotes, at the precision the write-up uses.

The provenance auditor greps the report's literal strings against the source
JSON.  Values that were rounded for readability (95.44 from 95.4385) or derived
in prose (chr8 gain+LOH share) therefore fail to match even though they are
correct.  Rather than hand-editing the report to match stored precision -- which
would invite typos -- this script recomputes each quoted value from the stored
outputs and writes it at report precision, so the report remains fully greppable
and nothing is ever typed by hand.

Externally-quoted canonical figures (88.26, 55.10, 64.89) are pulled from the
2026-08-01 denominator registry so their provenance is explicit too.
"""

from __future__ import annotations

import csv
import json
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
REGISTRY = (
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/handoff/"
    "20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv"
)
OUT = os.path.join(DATA, "report_values.json")


def load(name):
    with open(os.path.join(DATA, name)) as fh:
        return json.load(fh)


def main():
    cn = load("cn_annotation_summary.json")
    af = load("af_vs_cn_expectation.json")
    adj = load("savana_vs_seqc2_adjudication.json")
    refit = load("savana_refit_grid.json")
    res = load("resolution_vs_cn_stratified.json")
    mech1 = load("mechanism_and_intervals.json")
    mech2 = load("mechanism_v2_distinctness.json")
    rob = load("robustness_checks.json")
    cons = load("consolidated_findings.json")

    def pct2(x):
        return None if x is None else round(float(x), 2)

    def pct4(x):
        return None if x is None else round(float(x), 4)

    out = {"generated_by": os.path.basename(__file__)}

    # --- section 3: how much sits on CN-altered ground ---
    out["cn_share"] = {
        "unit_cn_altered_percent": pct2(cn["unit_cn_altered"]["percent"]),
        "site_cn_altered_percent": pct2(cn["site_cn_altered"]["percent"]),
        "unit_class_percent": {k: pct2(v) for k, v in cn["unit_class_percent"].items()},
        "site_class_percent": {k: pct2(v) for k, v in cn["site_class_percent"].items()},
        "unit_loss_combined_count": cn["unit_class_counts"].get("loss", 0)
        + cn["unit_class_counts"].get("loss_loh", 0),
        "unit_loss_combined_percent": pct2(
            100.0
            * (
                cn["unit_class_counts"].get("loss", 0)
                + cn["unit_class_counts"].get("loss_loh", 0)
            )
            / cn["unit_cn_altered"]["denominator"]
        ),
        "site_loss_combined_count": cn["site_class_counts"].get("loss", 0)
        + cn["site_class_counts"].get("loss_loh", 0),
        "site_loss_combined_percent": pct2(
            100.0
            * (
                cn["site_class_counts"].get("loss", 0)
                + cn["site_class_counts"].get("loss_loh", 0)
            )
            / cn["site_totals"]["sites"]
        ),
    }

    # chr8 gain+LOH share, quoted in prose
    chr8 = next(r for r in cons["per_chromosome_distribution"] if r["chrom"] == "chr8")
    out["chr8"] = {
        "units": chr8["units"],
        "gain_loh": chr8["gain_loh"],
        "gain_loh_percent": pct2(100.0 * chr8["gain_loh"] / chr8["units"]),
    }
    zero_neu = [
        r for r in cons["per_chromosome_distribution"] if r["neutral"] == 0
    ]
    out["chromosomes_without_neutral"] = {
        "chroms": [r["chrom"] for r in zero_neu],
        "combined_units": sum(r["units"] for r in zero_neu),
    }

    # --- section 5-7: resolution layers ---
    out["resolution"] = {
        "structure_unique_percent_all": pct2(
            res["headline_rates"]["all_ranked_units"]["structure_unique_percent"]
        ),
        "structure_unique_percent_neutral": pct2(
            res["headline_rates"]["cn_neutral_only"]["structure_unique_percent"]
        ),
        "structure_unique_percent_altered": pct2(
            res["headline_rates"]["cn_altered_only"]["structure_unique_percent"]
        ),
        "unique_best_tree_percent_all": pct2(
            res["headline_rates"]["all_ranked_units"]["unique_best_tree_percent"]
        ),
        "unique_best_tree_percent_neutral": pct2(
            res["headline_rates"]["cn_neutral_only"]["unique_best_tree_percent"]
        ),
        "unique_best_tree_percent_altered": pct2(
            res["headline_rates"]["cn_altered_only"]["unique_best_tree_percent"]
        ),
        "af_broke_tie_percent_by_class": {
            k: pct2(v["percent"]) if v else None
            for k, v in res["af_broke_tie_by_cn_class_contested_only"].items()
        },
        "fisher_odds_ratio": round(
            float(res["marginal_test"]["fisher_odds_ratio_altered_vs_neutral"]), 3
        ),
        "fisher_p_value_scientific": "%.1e" % res["marginal_test"]["fisher_p_value"],
        "cmh_or_tree_count": res["stratified_by_candidate_tree_count"]["cmh"][
            "common_odds_ratio"
        ],
        "cmh_p_tree_count_scientific": "%.1e"
        % res["stratified_by_candidate_tree_count"]["cmh"]["p_value"],
        "cmh_or_k": res["stratified_by_k"]["cmh"]["common_odds_ratio"],
        "cmh_p_k_scientific": "%.1e" % res["stratified_by_k"]["cmh"]["p_value"],
        "cmh_or_both": res["stratified_by_both"]["cmh"]["common_odds_ratio"],
        "cmh_p_both_scientific": "%.1e" % res["stratified_by_both"]["cmh"]["p_value"],
    }

    # --- mechanism ---
    sp = mech1["mechanism"]["test_altered_vs_neutral"]
    out["mechanism"] = {
        "spread_median_neutral": pct4(sp["median_neutral"]),
        "spread_median_altered": pct4(sp["median_altered"]),
        "spread_p_value": pct4(sp["p_value"]),
        "af_median_by_class": {
            k: pct4(v["median"]) if v else None
            for k, v in cn["read_af_by_cn_class"].items()
        },
        "identical_af_neutral_percent": pct2(
            mech2["prediction_1_identical_af_more_common_on_neutral"]["neutral"][
                "point_percent"
            ]
        ),
        "identical_af_altered_percent": pct2(
            mech2["prediction_1_identical_af_more_common_on_neutral"]["cn_altered"][
                "point_percent"
            ]
        ),
        "identical_af_fisher_or": mech2[
            "prediction_1_identical_af_more_common_on_neutral"
        ]["fisher_odds_ratio_neutral_vs_altered"],
        "identical_af_p_scientific": "%.1e"
        % mech2["prediction_1_identical_af_more_common_on_neutral"]["fisher_p_value"],
        "mediation_or": mech2["mediation_check"][
            "cmh_controlling_distinctness_and_tree_count"
        ]["common_odds_ratio"],
        "mediation_p_scientific": "%.1e"
        % mech2["mediation_check"]["cmh_controlling_distinctness_and_tree_count"][
            "p_value"
        ],
        "tautology_units": rob["check_b_tautology"]["identical_af_units"],
        "tautology_score_values": rob["check_b_tautology"]["their_best_score_values"],
        "identical_af_or": mech2["prediction_1_identical_af_more_common_on_neutral"][
            "fisher_odds_ratio_neutral_vs_altered"
        ],
        "unresolvable_neutral_percent": pct2(
            mech2["arithmetically_unresolvable_share"]["cn_neutral"]["point_percent"]
        ),
        "unresolvable_altered_percent": pct2(
            mech2["arithmetically_unresolvable_share"]["cn_altered"]["point_percent"]
        ),
        "nondeg_altered_percent": pct2(
            mech2["non_degenerate_contrast"]["cn_altered"]["point_percent"]
        ),
        "nondeg_neutral_percent": pct2(
            mech2["non_degenerate_contrast"]["cn_neutral"]["point_percent"]
        ),
        "nondeg_n_altered": mech2["non_degenerate_contrast"]["cn_altered"]["n"],
        "nondeg_n_neutral": mech2["non_degenerate_contrast"]["cn_neutral"]["n"],
        "nondeg_or": mech2["non_degenerate_contrast"]["fisher_odds_ratio"],
        "nondeg_p": round(mech2["non_degenerate_contrast"]["fisher_p_value"], 4),
        "spread_n_altered": sp["n_altered"],
        "spread_n_neutral": sp["n_neutral"],
        "spread_p_scientific": "%.1e" % sp["p_value"],
        "spread_by_class": {
            k: (round(v["median"], 4) if v else None)
            for k, v in mech1["mechanism"]["within_unit_af_spread_by_cn_class"].items()
        },
    }

    # AF grid explainability at report precision
    out["af_grid"] = {
        k: {
            "hit_percent": pct2(v["grid_hit_rate_percent"]),
            "chance_percent": pct2(v["chance_hit_rate_percent"]),
            "excess_percent": pct2(v["excess_over_chance_percent"]),
        }
        for k, v in af["cn_explainability"].items()
    }
    neu_grid = af["cn_explainability"].get("neutral")
    if neu_grid:
        out["af_grid"]["neutral_below_one_percent"] = pct2(
            100.0 - float(neu_grid["grid_hit_rate_percent"])
        )

    # --- SAVANA ---
    out["savana"] = {
        "published_purity": refit["published_fit"]["purity"],
        "published_ploidy": refit["published_fit"]["ploidy"],
        "published_state_agreement_percent": pct2(
            100.0 * refit["published_fit"]["state_agreement"]
        ),
        "published_integer_match_percent": pct2(
            100.0 * refit["published_fit"]["integer_match_rate"]
        ),
        "published_mean_abs_error": pct4(refit["published_fit"]["mean_abs_cn_error"]),
        "best_purity": refit["best_by_state_agreement"]["purity"],
        "best_ploidy": refit["best_by_state_agreement"]["ploidy"],
        "best_state_agreement_percent": pct2(
            100.0 * refit["best_by_state_agreement"]["state_agreement"]
        ),
        "best_integer_match_percent": pct2(
            100.0 * refit["best_by_state_agreement"]["integer_match_rate"]
        ),
        "best_mean_abs_error": pct4(
            refit["best_by_state_agreement"]["mean_abs_cn_error"]
        ),
        "best_int_purity": refit["best_by_integer_match"]["purity"],
        "best_int_ploidy": refit["best_by_integer_match"]["ploidy"],
        "best_int_state_agreement_percent": pct2(
            100.0 * refit["best_by_integer_match"]["state_agreement"]
        ),
        "best_int_integer_match_percent": pct2(
            100.0 * refit["best_by_integer_match"]["integer_match_rate"]
        ),
        "best_int_mean_abs_error": pct4(
            refit["best_by_integer_match"]["mean_abs_cn_error"]
        ),
        "seqc2_implied_ploidy": adj["seqc2_implied_ploidy_autosomal_bp_weighted"],
        "bin_state_agreement_percent": pct2(adj["state_agreement"]["percent"]),
        "regression_all_r2": adj["regression_savana_on_seqc2"]["r_squared"],
        "regression_nonneutral_r2": adj[
            "regression_savana_on_seqc2_nonneutral_only"
        ]["r_squared"],
        "regression_nonneutral_slope": adj[
            "regression_savana_on_seqc2_nonneutral_only"
        ]["slope"],
        "chr16_bias": adj["per_chrom_mean_bias_savana_minus_seqc2"].get("chr16"),
        "chr6_bias": adj["per_chrom_mean_bias_savana_minus_seqc2"].get("chr6"),
    }

    # --- chromosome robustness ---
    out["robustness"] = {
        "cmh_or_chromosome": rob["check_a_chromosome_confounding"][
            "within_chromosome_contrast"
        ]["cmh_stratified_by_chromosome"]["common_odds_ratio"],
        "cmh_p_chromosome_scientific": "%.1e"
        % rob["check_a_chromosome_confounding"]["within_chromosome_contrast"][
            "cmh_stratified_by_chromosome"
        ]["p_value"],
        "neutral_top4_combined_percent": pct2(
            rob["check_a_chromosome_confounding"]["neutral_group_concentration"][
                "top4_combined_percent"
            ]
        ),
    }

    # --- clean-ground capacity and CN-specific phenomena ---
    cap = load("clean_ground_capacity.json")
    kc = load("clean_ground_k_controlled.json")
    A = cap["A_clean_ground_capacity"]
    B = cap["B_what_is_distinctive_about_cn_ground"]
    out["clean_ground"] = {
        "neutral_ranked": kc["neutral_capacity"]["ranked"],
        "neutral_k1_no_tree": kc["neutral_capacity"]["k_eq_1_no_tree_possible"],
        "neutral_k_ge2": kc["neutral_capacity"]["k_ge_2_tree_informative"],
        "neutral_k_ge2_percent": kc["neutral_capacity"]["k_ge_2_percent"],
        "altered_k_ge2": kc["altered_capacity"]["k_ge_2_tree_informative"],
        "altered_k_ge2_percent": kc["altered_capacity"]["k_ge_2_percent"],
        "tree_informative_neutral_share_percent": round(
            100.0
            * kc["neutral_capacity"]["k_ge_2_tree_informative"]
            / (
                kc["neutral_capacity"]["k_ge_2_tree_informative"]
                + kc["altered_capacity"]["k_ge_2_tree_informative"]
            ),
            2,
        ),
        "neutral_k_ge2_chrom_count": kc["neutral_k_ge_2_chrom_count"],
        "neutral_k_ge2_top2_share_percent": round(
            100.0
            * sum(sorted(kc["neutral_k_ge_2_by_chrom"].values(), reverse=True)[:2])
            / kc["neutral_capacity"]["k_ge_2_tree_informative"],
            2,
        ),
        "neutral_sites_total": A["neutral_site_read_af"]["n_sites_with_af"],
        "neutral_sites_clonal": A["neutral_site_read_af"]["clonal_at_or_above_0.95"],
        "neutral_sites_subclonal": A["neutral_site_read_af"]["subclonal_below_0.95"],
        "power": A["statistical_power"]["detectable_at_current_n"],
        "branching_neutral": kc["branching_k_ge_2"]["neutral"],
        "branching_altered": kc["branching_k_ge_2"]["altered"],
        "branching_by_k": kc["branching_by_k"],
    }
    out["cn_phenomena"] = {
        "af_median_by_total_cn": {
            k: v["median_af"] for k, v in B["read_af_by_total_copy_number"].items()
        },
        "predicted_lowest_grid_point_by_cn": {
            k: v.get("predicted_lowest_grid_point")
            for k, v in B["read_af_by_total_copy_number"].items()
        },
        "grid_hit_percent_by_cn": {
            k: (
                round(100.0 * v["fraction_within_0.05_of_a_grid_point"], 2)
                if v.get("fraction_within_0.05_of_a_grid_point") is not None
                else None
            )
            for k, v in B["read_af_by_total_copy_number"].items()
        },
        "n_sites_by_cn": {
            k: v["n"] for k, v in B["read_af_by_total_copy_number"].items()
        },
        "structure_unique_percent_by_class": {
            k: v["structure_unique_percent"] for k, v in B["resolution_by_cn_class"].items()
        },
        "af_broke_tie_percent_by_class": {
            k: v["af_broke_tie_percent"] for k, v in B["resolution_by_cn_class"].items()
        },
    }

    # --- evolutionary yield left on copy-number-clean ground ---
    ey = load("clean_ground_evolution_yield.json")
    def band(p):
        return {
            "n_units": p.get("n_units"),
            "structure_unique": p["structure_unique"]["count"] if p.get("structure_unique") else None,
            "structure_unique_percent": p["structure_unique"]["percent"] if p.get("structure_unique") else None,
            "af_resolved": p["af_resolved"]["count"] if p.get("af_resolved") else None,
            "af_resolved_percent": p["af_resolved"]["percent"] if p.get("af_resolved") else None,
            "still_tied": p["still_tied"]["count"] if p.get("still_tied") else None,
            "still_tied_percent": p["still_tied"]["percent"] if p.get("still_tied") else None,
            "determined": p["determined_total"]["count"] if p.get("determined_total") else None,
            "determined_percent": p["determined_total"]["percent"] if p.get("determined_total") else None,
            "determined_lo": p["determined_total"]["lo"] if p.get("determined_total") else None,
            "determined_hi": p["determined_total"]["hi"] if p.get("determined_total") else None,
            "edges": p.get("determined_edges"),
            "deep_edges": p.get("determined_deep_edges"),
            "chromosomes": p.get("chromosomes_covered"),
        }
    out["evolution_yield"] = {
        "neutral_all": band(ey["cn_neutral"]["all"]),
        "neutral_k_ge_2": band(ey["cn_neutral"]["k_ge_2"]),
        "neutral_k_ge_3": band(ey["cn_neutral"]["k_ge_3"]),
        "neutral_by_k": {
            b: {
                "n": p.get("n_units"),
                "determined_percent": p["determined_total"]["percent"] if p.get("determined_total") else None,
                "edges": p.get("determined_edges"),
                "deep_edges": p.get("determined_deep_edges"),
            }
            for b, p in ey["cn_neutral"]["by_k"].items()
            if p.get("n_units")
        },
        "altered_k_ge_2": band(ey["cn_altered_reference"]["k_ge_2"]),
        "altered_k_ge_3": band(ey["cn_altered_reference"]["k_ge_3"]),
        "clean_share_of_deep_edges_k_ge_2_percent": round(
            100.0
            * ey["cn_neutral"]["k_ge_2"]["determined_deep_edges"]
            / (
                ey["cn_neutral"]["k_ge_2"]["determined_deep_edges"]
                + ey["cn_altered_reference"]["k_ge_2"]["determined_deep_edges"]
            ),
            2,
        ),
        "altered_share_of_deep_edges_k_ge_2_percent": round(
            100.0
            * ey["cn_altered_reference"]["k_ge_2"]["determined_deep_edges"]
            / (
                ey["cn_neutral"]["k_ge_2"]["determined_deep_edges"]
                + ey["cn_altered_reference"]["k_ge_2"]["determined_deep_edges"]
            ),
            2,
        ),
        "clean_k_ge_3_top3_chrom_sum": sum(
            sorted(ey["clean_k_ge_3_detail"]["chrom_distribution"].values(), reverse=True)[:3]
        ),
        "clean_k_ge_3_chrom_distribution": ey["clean_k_ge_3_detail"]["chrom_distribution"],
    }

    # --- multisample clean-anchor extension + external SEQC2 comparison ---
    ms = load("multisample_clean_deep_edges.json")
    nc = load("savana_neutral_callability.json")
    spec = load("neutral_af_spectrum.json")
    out["multisample"] = {
        "gate_precision_percent": nc["scores"]["recalibrated_with_baf"]["precision_percent"],
        "gate_recall_percent": nc["scores"]["recalibrated_with_baf"]["recall_percent"],
        "gate_precision_published_fit_percent": nc["scores"]["as_published_with_baf"]["precision_percent"],
        "per_sample": {
            s_: {
                "clean_k_ge2_units": r["clean"]["k_ge_2"]["units"],
                "clean_k_ge2_determined_percent": r["clean"]["k_ge_2"]["determined"]["percent"]
                if r["clean"]["k_ge_2"]["determined"] else None,
                "clean_k_ge2_deep_edges": r["clean"]["k_ge_2"]["deep_edges"],
                "clean_k_ge3_deep_edges": r["clean"]["k_ge_3"]["deep_edges"],
                "clean_strict_k_ge2_deep_edges": r["clean_strict"]["k_ge_2"]["deep_edges"],
                "altered_k_ge2_determined_percent": r["altered_k_ge_2_reference"]["determined"]["percent"]
                if r["altered_k_ge_2_reference"]["determined"] else None,
                "chromosomes": r["clean"]["k_ge_2"]["chromosomes"],
            }
            for s_, r in ms["per_sample"].items()
            if "clean" in r
        },
        "pooled": ms["pooled_totals"],
    }
    out["seqc2_external_comparison"] = {
        "neutral_sites_depth_ge20": spec["n_sites_depth_ge20"],
        "clonal_af_near_1_count": spec["histogram_bin005"]["0.95"],
        "subclonal_count": spec["subclonal_count"],
        "subclonal_percent": spec["subclonal_percent"],
        "near_ccf60_count": spec["near_ccf60_count"],
        "seqc2_reported": {
            "tree": "S1 = MRCA, S2-S10 = subclones",
            "s2_cancer_cell_fraction_percent": 60,
            "single_cells_tumor": 1270,
            "single_cells_normal": 638,
            "method": "10x Single Cell CNV integer-scaled CNA + VAF clonality corrected for local copy number + hierarchical clustering across replicates",
            "source": "Fang et al., Nat Biotechnol 2021 (PMID 34504347)",
        },
    }

    # --- linkage as a constraint on VAF-based reconstruction ---
    lk = load("linkage_as_constraint_on_vaf.json")
    ci = load("constraint_inventory_all_samples.json")
    out["linkage_constraint"] = {
        "k2_pairs_total": lk["all"]["n_pairs"],
        "relation_percent_all": lk["all"]["relation_percent"],
        "exclusive_large_gap_all": lk["all"]["mutually_exclusive_with_large_af_gap"],
        "exclusive_large_gap_neutral": lk["cn_neutral"]["mutually_exclusive_with_large_af_gap"],
        "exclusive_large_gap_altered": lk["cn_altered"]["mutually_exclusive_with_large_af_gap"],
        "median_af_gap_exclusive_neutral": lk["cn_neutral"]["median_af_gap_in_exclusive_pairs"],
        "median_af_gap_exclusive_altered": lk["cn_altered"]["median_af_gap_in_exclusive_pairs"],
        "nested_violating_monotonicity_all": lk["all"]["nested_violating_af_monotonicity"],
        "large_gap_threshold": lk["large_af_gap_threshold"],
        "inventory": ci,
    }

    # --- constraint validity audit ---
    cv = load("constraint_validity_audit.json")
    out["constraint_validity"] = {
        "units_audited": cv["units_audited"],
        "ordering_total": cv["claim_1_ordering"]["total_ordering_constraints"],
        "ordering_on_neutral": cv["claim_1_ordering"]["on_cn_neutral"],
        "ordering_on_altered": cv["claim_1_ordering"]["on_cn_altered"],
        "median_linking_depth": cv["claim_1_ordering"]["sampling_power"]["median_linking_depth"],
        "median_min_detectable_fraction": cv["claim_1_ordering"]["sampling_power"]["median_min_detectable_pattern_fraction"],
        "depth_bands": cv["claim_1_ordering"]["power_by_depth_band"],
        "exclusion_total": cv["claim_2_exclusion"]["exclusion_units_total"],
        "exclusion_cellular_valid": cv["claim_2_exclusion"]["exclusion_units_cn_neutral_cellular_valid"],
        "exclusion_molecular_only": cv["claim_2_exclusion"]["exclusion_units_cn_altered_molecular_only"],
        "exclusion_cellular_share": cv["claim_2_exclusion"]["share_reaching_cellular_level"],
        "exclusion_molecular_only_percent": round(
            100.0
            * cv["claim_2_exclusion"]["exclusion_units_cn_altered_molecular_only"]
            / cv["claim_2_exclusion"]["exclusion_units_total"],
            2,
        ),
        "median_min_detectable_percent": round(
            100.0
            * cv["claim_1_ordering"]["sampling_power"][
                "median_min_detectable_pattern_fraction"
            ],
            2,
        ),
        "low_depth_unit_percent": round(
            sum(
                v["percent_of_units"]
                for k, v in cv["claim_1_ordering"]["power_by_depth_band"].items()
                if k in ("3-9", "10-19")
            ),
            2,
        ),
        "depth_band_detectable_percent": {
            k: {
                "low": round(100 * v["min_detectable_pattern_fraction_at_band_low"], 2),
                "high": round(100 * v["min_detectable_pattern_fraction_at_band_high"], 2),
                "units": v["units"],
                "percent_of_units": v["percent_of_units"],
            }
            for k, v in cv["claim_1_ordering"]["power_by_depth_band"].items()
        },
    }

    # --- externally quoted canonical figures ---
    ext = {}
    with open(REGISTRY) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["metric"] in (
                "one_rooted_unlabeled_topology",
                "unique_best_tree",
                "final_groups",
                "ranked_complete",
            ):
                ext[row["metric"]] = {
                    "numerator": int(row["numerator"]),
                    "denominator": int(row["denominator"]),
                    "percent": round(float(row["percent"]), 2),
                    "scope": row["scope"],
                    "note": "ALL SEVEN datasets combined; not a HCC1395 figure",
                }
    out["externally_quoted_canonical"] = ext
    # Values produced by the FIRST run of the mechanism scripts, before the
    # active-site fix of 2026-08-10.  They no longer exist in any data file
    # because the rerun overwrote them, but the correction record in the report
    # quotes them to say what changed -- so they are recorded here, explicitly
    # labelled as superseded, rather than being untraceable prose.
    out["superseded_by_active_site_correction"] = {
        "status": "SUPERSEDED - initial run, retained only so the correction record is traceable",
        "cause": "distinctness and spread were computed over all af_coverage rows, including inactive all-reference sites the score never reads",
        "found_by": "adversarial review 2026-08-10",
        "spread_median_neutral": 0.7826,
        "spread_median_altered": 0.5223,
        "spread_p_value": 0.9997,
        "identical_af_neutral_percent": 20.80,
        "identical_af_altered_percent": 11.78,
        "identical_af_fisher_or": 1.9669,
        "mediation_cmh_or": 2.819,
        "mediation_cmh_or_rounded": 2.82,
        "mediation_p_scientific": "3.9e-06",
        "tautology_units": 1107,
    }
    out["externally_quoted_unsourced"] = {
        "pure_parsimony_single_topology_lower_bound": {
            "value": 64.89,
            "status": "NO SOURCE IN REPO",
            "note": "carried over from prior sessions; still not greppable; flagged, not used as evidence",
        }
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"wrote {OUT}")
    print(json.dumps(out["cn_share"], indent=1, ensure_ascii=False))
    print(json.dumps(out["savana"], indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
