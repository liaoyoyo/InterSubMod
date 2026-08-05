#!/usr/bin/env python3
"""Verify all pre-registered layered reconstruction v2 invariants."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pysam


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def add_check(checks, name, passed, observed, expected):
    checks.append({
        "name": name,
        "pass": bool(passed),
        "observed": observed,
        "expected": expected,
    })


def verify_sample(sample_meta, sample_dir):
    sample = sample_meta["sample"]
    layered_path = sample_dir / f"layered_reconstruction_{sample}.json"
    region_path = sample_dir / f"layered_region_view_{sample}.json"
    site_summary_path = sample_dir / f"ssnv_site_ledger_{sample}.summary.json"
    site_ledger_path = sample_dir / f"ssnv_site_ledger_{sample}.tsv.gz"
    site_index_path = Path(str(site_ledger_path) + ".tbi")
    checks = []
    add_check(checks, "layered_exists", layered_path.is_file(), str(layered_path), "existing file")
    add_check(checks, "region_view_exists", region_path.is_file(), str(region_path), "existing file")
    add_check(checks, "site_ledger_exists", site_ledger_path.is_file(), str(site_ledger_path), "existing file")
    add_check(checks, "site_ledger_index_exists", site_index_path.is_file(), str(site_index_path), "existing file")
    add_check(checks, "site_ledger_summary_exists", site_summary_path.is_file(), str(site_summary_path), "existing file")
    if not layered_path.is_file() or not region_path.is_file():
        return {"sample": sample, "pass": False, "checks": checks}

    layered = json.loads(layered_path.read_text(encoding="utf-8"))
    view = json.loads(region_path.read_text(encoding="utf-8"))
    detail = layered.get("detail", [])
    census = view.get("census", {})
    l1 = layered.get("L1_ssnv_algorithm", {})
    funnel = census.get("funnel", {})
    tag_census = layered.get("read_tag_census", {})
    site_summary = json.loads(site_summary_path.read_text(encoding="utf-8")) if site_summary_path.is_file() else {}
    ledger_queryable = False
    if site_ledger_path.is_file() and site_index_path.is_file():
        try:
            with pysam.TabixFile(str(site_ledger_path)) as indexed:
                probe_contig = next(iter(indexed.contigs), None)
                ledger_queryable = bool(probe_contig and next(indexed.fetch(probe_contig), None))
        except (OSError, ValueError):
            ledger_queryable = False
    add_check(checks, "site_ledger_tabix_queryable", ledger_queryable, ledger_queryable, True)

    add_check(checks, "schema_v2", layered.get("schema_version") == "2.0" and view.get("schema_version") == "2.0",
              [layered.get("schema_version"), view.get("schema_version")], ["2.0", "2.0"])
    ref_in_primary = sum(1 for u in detail if u.get("reference_only") and u.get("is_primary_lineage"))
    h3_in_primary = sum(1 for u in detail if u.get("family") == "3" and u.get("is_primary_lineage"))
    h4_in_primary = sum(1 for u in detail if u.get("family") == "4" and u.get("is_primary_lineage"))
    add_check(checks, "H1_reference_not_primary", ref_in_primary == 0, ref_in_primary, 0)
    add_check(checks, "H2_H3_not_primary", h3_in_primary == 0, h3_in_primary, 0)
    add_check(checks, "HP4_not_primary", h4_in_primary == 0, h4_in_primary, 0)

    role_counts = {
        "primary": sum(1 for u in detail if u.get("is_primary_lineage")),
        "reference": sum(1 for u in detail if u.get("reference_only")),
        "H3_auxiliary": sum(1 for u in detail if u.get("is_h3_auxiliary") and u.get("mutation_bearing")),
        "H4_auxiliary": sum(1 for u in detail if u.get("is_h4_auxiliary") and u.get("mutation_bearing")),
        "unphased": sum(1 for u in detail if u.get("family") == "none"),
    }
    add_check(checks, "role_count_primary_matches_L1", role_counts["primary"] == l1.get("n_primary_lineage_units"),
              role_counts["primary"], l1.get("n_primary_lineage_units"))
    add_check(checks, "role_count_reference_matches_L1", role_counts["reference"] == l1.get("n_reference_only_controls"),
              role_counts["reference"], l1.get("n_reference_only_controls"))

    eligible_skipped = sum(1 for u in detail if not u.get("capped") and u.get("verification_skipped"))
    verification_fail = sum(1 for u in detail if u.get("verification_status") == "fail")
    add_check(checks, "H4_zero_eligible_skipped", eligible_skipped == 0, eligible_skipped, 0)
    add_check(checks, "H4_zero_verification_fail", verification_fail == 0, verification_fail, 0)
    add_check(checks, "H4_all_eligible_pass", l1.get("all_eligible_V1V7_pass") is True,
              l1.get("all_eligible_V1V7_pass"), True)

    noncapped = [u for u in detail if not u.get("capped")]
    incomplete = sum(1 for u in noncapped if not u.get("analysis_candidate_set_complete"))
    count_mismatch = sum(1 for u in noncapped if u.get("analysis_trees_generated") != u.get("n_trees"))
    missing_exact_shape = sum(1 for u in noncapped if u.get("n_distinct_shapes_exact") is None)
    add_check(checks, "H6_complete_candidate_analysis", incomplete == 0, incomplete, 0)
    add_check(checks, "H6_candidate_count_matches", count_mismatch == 0, count_mismatch, 0)
    add_check(checks, "H6_exact_shape_present", missing_exact_shape == 0, missing_exact_shape, 0)
    negative_tree_hidden = sum(
        1 for u in detail for tree in u.get("trees", []) if tree.get("n_hidden", 0) < 0)
    add_check(checks, "no_negative_tree_hidden", negative_tree_hidden == 0, negative_tree_hidden, 0)

    add_check(checks, "H5_upstream_funnel_conservation", funnel.get("check_scope_conservation") is True,
              funnel.get("check_scope_conservation"), True)
    add_check(checks, "H5_autosomal_matches_upstream", funnel.get("check_autosomal_matches_upstream") is True,
              funnel.get("check_autosomal_matches_upstream"), True)
    add_check(checks, "H5_six_branch_conservation", funnel.get("check_six_branch_conservation") is True,
              funnel.get("check_six_branch_conservation"), True)
    add_check(checks, "read_tag_sidecar_exact_join", tag_census.get("check_exact_sidecar_join") is True,
              tag_census, "missing=0, conflicts=0, exact_matches=alignment exposures")
    phase_regions = tag_census.get("phase_set_region_counts", {})
    add_check(
        checks,
        "phase_set_region_census_conserves",
        sum(phase_regions.values()) == tag_census.get("n_regions_with_phase_set_census"),
        {"counts": phase_regions, "total": tag_census.get("n_regions_with_phase_set_census")},
        "none+one+multiple equals all retained regions",
    )
    add_check(
        checks,
        "phase_set_is_QC_not_topology",
        layered.get("analysis_contract", {}).get("PS")
        == "preserved for phase-block QC; not used as a topology edge or lineage label",
        layered.get("analysis_contract", {}).get("PS"),
        "preserved for phase-block QC; not used as a topology edge or lineage label",
    )
    add_check(checks, "site_ledger_all_records_conserved", site_summary.get("pass") is True,
              site_summary.get("checks"), "all site-ledger checks true")

    missing_cn = sample_meta.get("cn_source") == "unavailable"
    if missing_cn:
        measured_cn = sum(1 for u in detail if u.get("cn") != "unavailable")
        candidate_cn = sum(1 for u in detail if (u.get("L2_m_channel") or {}).get("verdict") == "candidate_keep")
        add_check(checks, "H3_missing_CN_all_unavailable", measured_cn == 0, measured_cn, 0)
        add_check(checks, "H3_missing_CN_zero_candidate", candidate_cn == 0, candidate_cn, 0)

    region_mismatch = 0
    for region in view.get("regions", []):
        expected = len({
            u["family"] for u in region.get("lineages", []) if u.get("is_primary_lineage")
        })
        if expected != region.get("hp_multiplicity"):
            region_mismatch += 1
    add_check(checks, "region_primary_multiplicity_recomputes", region_mismatch == 0, region_mismatch, 0)

    spans = [r.get("end", 0) - r.get("start", 0) for r in view.get("regions", [])]
    metrics = {
        "roles": role_counts,
        "n_regions": census.get("n_regions"),
        "hp_multiplicity": census.get("hp_multiplicity"),
        "region_determinacy": census.get("region_determinacy"),
        "verification_status": l1.get("verification_status"),
        "n_capped": l1.get("n_verification_not_applicable_capped"),
        "funnel": funnel,
        "read_tag_census": tag_census,
        "site_ledger": site_summary,
        "site_ledger_sha256": sha256(site_ledger_path) if site_ledger_path.is_file() else None,
        "site_ledger_index_sha256": sha256(site_index_path) if site_index_path.is_file() else None,
        "region_span_gt_50kb": sum(1 for x in spans if x > 50000),
        "region_span_max": max(spans, default=0),
        "layered_sha256": sha256(layered_path),
        "region_view_sha256": sha256(region_path),
    }
    return {"sample": sample, "biological_id": sample_meta.get("biological_id"),
            "pass": all(c["pass"] for c in checks), "checks": checks, "metrics": metrics}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.run_root / "verification_summary.json"
    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    results = [verify_sample(meta, args.run_root / "samples" / meta["sample"])
               for meta in manifest["samples"]]
    summary = {
        "schema_version": "2.0",
        "run_root": str(args.run_root.resolve()),
        "input_manifest": str(args.input_manifest.resolve()),
        "dataset_count": len(results),
        "biological_sample_count": len({r.get("biological_id") for r in results}),
        "all_pass": all(r["pass"] for r in results),
        "n_pass": sum(r["pass"] for r in results),
        "n_fail": sum(not r["pass"] for r in results),
        "samples": results,
    }
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tsv = output.with_suffix(".tsv")
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "pass", "primary", "reference", "H3_auxiliary", "H4_auxiliary", "regions", "capped", "failed_checks"])
        for result in results:
            metrics = result.get("metrics", {})
            roles = metrics.get("roles", {})
            failed = ",".join(c["name"] for c in result["checks"] if not c["pass"])
            writer.writerow([result["sample"], result["pass"], roles.get("primary"), roles.get("reference"),
                             roles.get("H3_auxiliary"), roles.get("H4_auxiliary"), metrics.get("n_regions"),
                             metrics.get("n_capped"), failed])
    print(f"VERIFY layered-v2: {summary['n_pass']}/{summary['dataset_count']} pass -> {output}")
    for result in results:
        failed = [c["name"] for c in result["checks"] if not c["pass"]]
        print(f"  {result['sample']}: {'PASS' if result['pass'] else 'FAIL'}" + (f" ({','.join(failed)})" if failed else ""))
    raise SystemExit(0 if summary["all_pass"] else 1)


if __name__ == "__main__":
    main()
