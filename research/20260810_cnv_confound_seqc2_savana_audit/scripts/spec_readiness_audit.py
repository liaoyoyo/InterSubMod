#!/usr/bin/env python3
"""What does the sr2b/sr2c joint-inference spec need, and what do we already have?

The lab tutorial (sr2b: CCF + grouping, sr2c: local molecular phylogeny) defines
a joint model

    theta_mk = rho * c_k * mu_m / (rho * CN_m + 2(1 - rho))

with mutation-to-group assignment z, tree topology T, copy genotypes G and cell
fractions w optimised together; S1 (single-variant windows) fitted by
beta-binomial, S2 (multi-variant windows) by a multinomial over observed read
states, and a one-way prior from S1 into S2 with escape mass lambda.

sr2c states plainly that it is "a specification, not yet validated".  This
script audits what the specification would need against what the frozen
canonical outputs already contain, and computes in place the quantities that
can be computed today:

  * the molecular state distribution nu_{u,h} per unit, from populations_by_hp
  * the marginalisation burden: how much of the evidence is partial (carries X)
    rather than fully covering, since sr2c requires unobserved positions to be
    marginalised rather than assumed reference
  * the observable support behind each unit, which bounds the error floor that
    sr2c says is the only one measurable in place

It does NOT fit the model.  The purpose is to establish which inputs exist,
which are derivable, and which are genuinely missing.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
ANNOT = os.path.join(DATA, "hcc1395_unit_cn_annotation.jsonl")
MLHP = (
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/"
    "samples/HCC1395/HCC1395.exact_ps_mlhp.json"
)
SAVANA = (
    "/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna/"
    "HCC1395_segmented_absolute_copy_number.tsv"
)
OUT = os.path.join(DATA, "spec_readiness_audit.json")


def main():
    cn_by_region = {}
    with open(ANNOT) as fh:
        for ln in fh:
            r = json.loads(ln)
            if r.get("cn_status") == "ANNOTATED":
                cn_by_region[r["region_id"]] = r["seqc2_class"]

    d = json.load(open(MLHP))
    groups = d.get("groups") or []

    n_groups = len(groups)
    cn_available = sum(1 for g in groups if g.get("cn_data_available"))
    vaf_eligible = sum(1 for g in groups if g.get("vaf_eligible"))

    full_mol = 0
    partial_mol = 0
    units_with_full = 0
    units_multi_state = 0
    state_counts = []
    full_share = []
    nu_top = []
    support_per_unit = []
    x_positions_share = []

    by_cn = defaultdict(lambda: {"units": 0, "full": 0, "partial": 0})

    for g in groups:
        pops = g.get("populations_by_hp") or {}
        subs = g.get("subread_groups_by_hp") or {}
        k = int(g.get("n_sSNV") or 0)
        if k < 2:
            continue
        f = sum(sum(v.values()) for v in pops.values())
        p = sum(sum(v.values()) for v in subs.values())
        full_mol += f
        partial_mol += p
        if f:
            units_with_full += 1
        states = sum(len(v) for v in pops.values())
        state_counts.append(states)
        if states > 1:
            units_multi_state += 1
        tot = f + p
        if tot:
            full_share.append(f / tot)
        support_per_unit.append(tot)

        # nu: molecular proportions over fully-observed states
        if f:
            allv = []
            for v in pops.values():
                allv.extend(v.values())
            allv.sort(reverse=True)
            nu_top.append(allv[0] / f)

        # how many position slots are X (unobserved) among partial patterns
        xs = tot_slots = 0
        for v in subs.values():
            for pat, c in v.items():
                xs += pat.count("X") * c
                tot_slots += len(pat) * c
        if tot_slots:
            x_positions_share.append(xs / tot_slots)

        cls = cn_by_region.get(g.get("region_id"))
        if cls:
            b = by_cn[cls]
            b["units"] += 1
            b["full"] += f
            b["partial"] += p

    def stats(v, name):
        if not v:
            return None
        s = sorted(v)
        n = len(s)
        return {
            "n": n,
            "median": round(s[n // 2], 4),
            "q1": round(s[n // 4], 4),
            "q3": round(s[(3 * n) // 4], 4),
            "mean": round(sum(s) / n, 4),
        }

    # ---- allele-specific CN availability from SAVANA ----
    seg_tot = seg_minor = seg_het = 0
    with open(SAVANA) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            seg_tot += 1
            try:
                mi = float(p[col["minorAlleleCopyNumber"]])
            except (ValueError, IndexError, KeyError):
                continue
            seg_minor += 1
            if mi >= 0.5:
                seg_het += 1

    inputs = {
        "rho_purity": {
            "spec_symbol": "rho",
            "status": "AVAILABLE_WITH_CAVEAT",
            "what_we_have": "HCC1395 recalibrated to 1.0 and externally corroborated (SEQC2-implied ploidy 2.9104, refit grid 2.95, KB 2.85); H1437/H2009/HCC1954 use SAVANA published purity, passed BAF self-consistency",
            "gap": "no external purity truth for the three non-HCC1395 samples",
        },
        "CN_m_total": {
            "spec_symbol": "CN_m",
            "status": "AVAILABLE",
            "what_we_have": "SEQC2 integer CN for HCC1395 (gain 1502 + loss 63 segments); SAVANA copyNumber elsewhere",
            "gap": "none for total CN",
        },
        "allele_specific_CN": {
            "spec_symbol": "G (copy genotypes)",
            "status": "AVAILABLE_UNVALIDATED",
            "what_we_have": "SAVANA minorAlleleCopyNumber present on %.1f%% of segments (%d/%d); %d segments have minor CN >= 0.5"
            % (100 * seg_minor / seg_tot, seg_minor, seg_tot, seg_het),
            "gap": "never validated against truth; SEQC2 publishes no allele-specific CN to check it against",
        },
        "hp_to_allele_mapping": {
            "spec_symbol": "needed to turn CN_m into the copies carried by the observed haplotype",
            "status": "MISSING_BUT_DERIVABLE",
            "what_we_have": "germline phasing already defines HP families and the phased VCF supplies 1.69M autosomal phased het SITE POSITIONS with each allele's HP assignment; SAVANA gives major/minor CN per segment",
            "gap": "the link between the two is not computed: which HP family carries the major allele. A pilot showed the phased germline VCF ALONE is insufficient - its AD field comes from the matched normal (HCC1395BL), so every het site is ~50:50 by construction and LOH separation AUC was 0.4341. Per-haplotype DEPTHS must be counted on the HP-tagged TUMOUR bam; the existing HP tags mean no re-phasing is needed",
            "why_it_matters": "mu_m / CN cannot be turned into a per-haplotype multiplicity without it; this is the single blocking input",
        },
        "mu_m_multiplicity": {
            "spec_symbol": "mu_m",
            "status": "NOT_YET_ESTIMATED",
            "what_we_have": "strong empirical evidence it is operating: read-AF falls monotonically with CN (1.0 -> 0.19 from CN 1 to 7) and concentrates on the m/c grid (85.37% within 0.05 at CN=8)",
            "gap": "requires hp_to_allele_mapping; currently only the aggregate signature is measured, not per-mutation mu",
        },
        "read_states_S2": {
            "spec_symbol": "D_v, X_u",
            "status": "AVAILABLE",
            "what_we_have": "populations_by_hp gives fully-covering molecule patterns with counts; subread_groups_by_hp gives partial patterns carrying X",
            "gap": "none - this is exactly the S2 evidence the spec consumes",
        },
        "alt_ref_counts_S1": {
            "spec_symbol": "a_m, d_m",
            "status": "AVAILABLE",
            "what_we_have": "col_coverage_by_hp per site [ref, alt] within each HP family",
            "gap": "these are HP-family-specific, not caller VAF; the spec's theta is written for sample-level VAF and must be re-derived for the HP-conditional denominator",
        },
        "error_floor_delta": {
            "spec_symbol": "delta_u = f(epsilon, eta)",
            "status": "MEASURABLE_NOT_MEASURED",
            "what_we_have": "the spec states only the error floor is measurable in place; the partial/full pattern inventory needed to measure it exists",
            "gap": "eta (whole-molecule family mislabeling) not yet estimated from data; spec notes it dominates point-wise error at k=3 (3e-3 vs 1e-7)",
        },
        "lambda_escape_mass": {
            "spec_symbol": "lambda",
            "status": "JUDGEMENT_PARAMETER",
            "what_we_have": "spec calls it the only parameter chosen by judgement and requires sensitivity analysis",
            "gap": "our CN-neutral AF spectrum is continuous with no discrete CCF atoms (267 subclonal sites, median 0.4583), which argues against a small lambda; this is evidence for the sensitivity analysis, not a value",
        },
        "tree_topology_T_and_ancestry_E": {
            "spec_symbol": "T, E_u",
            "status": "ALREADY_COMPUTED",
            "what_we_have": "7,554 ordering constraints for HCC1395 (78,655 across 7 samples); pairwise four-gamete compatibility already tested",
            "gap": "none; this is what we contribute to the spec rather than need from it",
        },
    }

    out = {
        "generated_by": os.path.basename(__file__),
        "spec_source": [
            "https://ccu-bioinformatics-lab.github.io/lab-tutorial/sr2b.html",
            "https://ccu-bioinformatics-lab.github.io/lab-tutorial/sr2c.html",
        ],
        "spec_status_per_its_own_text": "sr2c states: 'This is a specification, not yet validated'",
        "canonical_mlhp_inventory": {
            "groups_total": n_groups,
            "groups_with_cn_data": cn_available,
            "groups_vaf_eligible": vaf_eligible,
            "note": "cn_data_available is False for every group in the frozen run - copy number was never integrated into the unit-level evidence",
        },
        "s2_evidence_scale": {
            "units_k_ge_2": len(state_counts),
            "fully_covering_molecules": full_mol,
            "partial_molecules_with_X": partial_mol,
            "partial_share_percent": round(
                100 * partial_mol / (full_mol + partial_mol), 2
            )
            if (full_mol + partial_mol)
            else None,
            "units_with_any_full_molecule": units_with_full,
            "units_with_more_than_one_observed_state": units_multi_state,
            "observed_states_per_unit": stats(state_counts, "states"),
            "full_evidence_share_per_unit": stats(full_share, "full_share"),
            "dominant_state_proportion_nu": stats(nu_top, "nu_top"),
            "molecules_per_unit": stats(support_per_unit, "support"),
            "x_slot_share_in_partial_patterns": stats(x_positions_share, "x_share"),
        },
        "s2_evidence_by_cn_class": {
            k: {
                "units": v["units"],
                "full_molecules": v["full"],
                "partial_molecules": v["partial"],
                "partial_share_percent": round(
                    100 * v["partial"] / (v["full"] + v["partial"]), 2
                )
                if (v["full"] + v["partial"])
                else None,
            }
            for k, v in sorted(by_cn.items())
        },
        "required_inputs": inputs,
    }

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(f"mlhp groups: {n_groups}; cn_data_available: {cn_available}; vaf_eligible: {vaf_eligible}")
    print(f"S2 evidence: full={full_mol:,} partial(X)={partial_mol:,} "
          f"partial share={out['s2_evidence_scale']['partial_share_percent']}%")
    print(f"units k>=2 with S2: {len(state_counts)}; multi-state: {units_multi_state}")
    print(f"states/unit: {out['s2_evidence_scale']['observed_states_per_unit']}")
    print(f"full-evidence share/unit: {out['s2_evidence_scale']['full_evidence_share_per_unit']}")
    print(f"X slot share in partial: {out['s2_evidence_scale']['x_slot_share_in_partial_patterns']}")
    print(f"\nallele-specific CN: minor present on {seg_minor}/{seg_tot} segments ({100*seg_minor/seg_tot:.1f}%)")
    print("\nrequired-input status:")
    for k, v in inputs.items():
        print(f"  {k:28s} {v['status']}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
