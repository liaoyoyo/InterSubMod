#!/usr/bin/env python3
"""Compare canonical LongPhase-S PASS and ClairS-input tree backbones across all datasets."""

import argparse
import csv
import json
from pathlib import Path


def load_view(run_root, sample):
    path = run_root / "samples" / sample / f"layered_region_view_{sample}.json"
    return json.loads(path.read_text(encoding="utf-8"))


def retained_positions(run_root, sample):
    values = set()
    for path in (run_root / "samples" / sample).glob("mlhp_part_*.json"):
        doc = json.loads(path.read_text(encoding="utf-8"))
        for group in doc.get("groups", []):
            values.update((group["chrom"], int(pos)) for pos in group.get("positions", []))
    return values


def primary_units(run_root, sample):
    path = run_root / "samples" / sample / f"layered_reconstruction_{sample}.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    values = {}
    for unit in document.get("detail", []):
        if not unit.get("is_primary_lineage") or unit.get("capped"):
            continue
        key = (unit["chrom"], int(unit["start"]), int(unit["end"]), str(unit["family"]))
        if key in values:
            raise RuntimeError(f"duplicate primary unit key for {sample}: {key}")
        digest = unit.get("analysis_tree_digest_sha256")
        if not digest:
            raise RuntimeError(f"missing complete candidate-tree digest for {sample}: {key}")
        values[key] = digest
    return values


def metrics(view):
    census = view["census"]
    l1 = census["L1"]
    roles = {
        "primary": l1["n_primary_lineage_units"],
        "reference": l1["n_reference_only_controls"],
        "H3_auxiliary": l1["n_unresolved_H3_auxiliary"],
        "H4_auxiliary": l1["n_shared_H4_auxiliary"],
    }
    primary = roles["primary"] or 1
    determined = l1["determinacy_primary_lineage"].get("determined", 0)
    n_regions = census["n_regions"] or 1
    multi = census["hp_multiplicity"].get("2", 0)
    all_det = census["region_determinacy"].get("all_determined", 0)
    with_primary = census["n_regions_with_primary_lineage"] or 1
    return {
        "universe": census["U1_sSNV_somatic_total"],
        "autosomal": census["U1_sSNV_scope_chr1_22"],
        "retained_sSNV": census["funnel"]["L6_retained_sSNV"],
        "n_regions": census["n_regions"],
        **roles,
        "determined_primary": determined,
        "determined_primary_pct": 100 * determined / primary,
        "multiHP_regions": multi,
        "multiHP_pct_all_regions": 100 * multi / n_regions,
        "all_determined_regions": all_det,
        "all_determined_pct_primary_regions": 100 * all_det / with_primary,
    }


def jaccard(a, b):
    return len(a & b) / len(a | b) if a or b else 1.0


def verdict(deltas, retained_jaccard=1.0, unit_jaccard=1.0, topology_concordance=1.0):
    maximum = max(abs(value) for value in deltas)
    if (maximum < 5 and retained_jaccard >= 0.95 and unit_jaccard >= 0.90
            and topology_concordance >= 0.95):
        return "robust_all_gates"
    if (maximum < 10 and retained_jaccard >= 0.80 and unit_jaccard >= 0.70
            and topology_concordance >= 0.80):
        return "moderately_sensitive"
    return "backbone_sensitive"


def validate_run_contract(run_root, expected_samples, expected_contract):
    verification_path = run_root / "verification_summary.json"
    success_path = run_root / "_SUCCESS"
    if not verification_path.is_file() or not success_path.is_file():
        raise SystemExit(f"run is incomplete: {run_root}")
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    verified_samples = [item.get("sample") for item in verification.get("samples", [])]
    if (verification.get("all_pass") is not True or verification.get("n_pass") != 7
            or len(verified_samples) != 7 or set(verified_samples) != set(expected_samples)):
        raise SystemExit(f"run verification/sample set failed: {run_root}")
    frozen_manifest_path = run_root / "input_manifest.json"
    frozen_lock_path = run_root / "frozen_input_lock.json"
    if frozen_manifest_path.is_file():
        frozen = json.loads(frozen_manifest_path.read_text(encoding="utf-8"))
        frozen_samples = [item.get("sample") for item in frozen.get("samples", [])]
        if frozen.get("tree_input_contract") != expected_contract or frozen_samples != expected_samples:
            raise SystemExit(f"v2 run tree contract/sample order failed: {run_root}")
    elif frozen_lock_path.is_file():
        frozen = json.loads(frozen_lock_path.read_text(encoding="utf-8"))
        frozen_samples = [item.get("sample") for item in frozen.get("samples", [])]
        expected_profiles = {
            "longphase_recalibrated_PASS": (
                "comprehensive_validation",
                "longphase_s_recalibrated_FILTER_PASS",
                "longphase_recalibrated_pass_vcf",
            ),
            "clairs_PASS_input": (
                "backbone_sensitivity",
                "clairs_FILTER_PASS_sensitivity",
                "caller_pass_baseline_vcf",
            ),
        }
        if expected_contract not in expected_profiles:
            raise SystemExit(f"unsupported expected v3 tree contract: {expected_contract}")
        expected_task, expected_tree, selected_role = expected_profiles[expected_contract]
        exact_roles = all(
            set(item.get("somatic", {}))
            == {
                "caller_raw_vcf", "longphase_input_vcf", "caller_pass_baseline_vcf",
                "longphase_recalibrated_all_vcf", "longphase_recalibrated_pass_vcf", "tree_vcf",
            }
            and item["somatic"]["tree_vcf"] == item["somatic"][selected_role]
            and item["somatic"]["longphase_recalibrated_pass_vcf"]
            != item["somatic"]["longphase_recalibrated_all_vcf"]
            for item in frozen.get("samples", [])
        )
        if (len(frozen_samples) != 7 or set(frozen_samples) != set(expected_samples)
                or frozen.get("all_pass") is not True
                or frozen.get("analysis_contract", {}).get("production_filter_policy", {}).get("name")
                != "production_default_filter"
                or frozen.get("analysis_contract", {}).get("task_type") != expected_task
                or frozen.get("analysis_contract", {}).get("tree_input_contract") != expected_tree
                or not exact_roles):
            raise SystemExit(f"v3 frozen lock/tree role failed: {run_root}")
    else:
        raise SystemExit(f"recognized frozen manifest/lock missing: {run_root}")
    for sample in expected_samples:
        sample_root = run_root / "samples" / sample
        layered_path = sample_root / f"layered_reconstruction_{sample}.json"
        site_path = sample_root / f"ssnv_site_ledger_{sample}.summary.json"
        if not layered_path.is_file() or not site_path.is_file():
            raise SystemExit(f"{sample} backbone scientific outputs are incomplete: {run_root}")
        layered = json.loads(layered_path.read_text(encoding="utf-8"))
        tags = layered.get("read_tag_census", {})
        phase_counts = tags.get("phase_set_region_counts", {})
        site = json.loads(site_path.read_text(encoding="utf-8"))
        expected_site_tree = (
            "longphase_recalibrated_PASS"
            if expected_contract == "longphase_recalibrated_PASS"
            else "clairs_PASS_input"
        )
        if (tags.get("check_exact_sidecar_join") is not True
                or sum(phase_counts.values()) != tags.get("n_regions_with_phase_set_census")
                or layered.get("analysis_contract", {}).get("PS")
                != "preserved for phase-block QC; not used as a topology edge or lineage label"
                or site.get("pass") is not True
                or site.get("longphase_input_contract") != "clairs_raw_all"
                or site.get("tree_contract") != expected_site_tree):
            raise SystemExit(f"{sample} backbone science/tag/site-ledger contract failed: {run_root}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-run", required=True, type=Path)
    parser.add_argument("--clairs-run", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    canonical_v2 = (
        manifest.get("schema_version") == "2.1"
        and manifest.get("tree_input_contract") == "longphase_recalibrated_PASS"
    )
    canonical_v3 = (
        manifest.get("schema_version") == "3.0.0"
        and manifest.get("analysis_contract", {}).get("task_type") == "comprehensive_validation"
        and manifest.get("analysis_contract", {}).get("tree_input_contract")
        == "longphase_s_recalibrated_FILTER_PASS"
    )
    if (manifest.get("dataset_count") != 7 or not (canonical_v2 or canonical_v3)):
        raise SystemExit("input manifest is not the full canonical LongPhase-S PASS manifest")
    expected_samples = [item["sample"] for item in manifest.get("samples", [])]
    if len(expected_samples) != 7 or len(set(expected_samples)) != 7:
        raise SystemExit("input manifest must contain seven unique datasets")
    validate_run_contract(args.base_run, expected_samples, "longphase_recalibrated_PASS")
    validate_run_contract(args.clairs_run, expected_samples, "clairs_PASS_input")
    rows = []
    for meta in manifest["samples"]:
        sample = meta["sample"]
        label = "clairs_PASS_input_vs_longphase_PASS"
        alt_run = args.clairs_run
        base_view = load_view(args.base_run, sample)
        alt_view = load_view(alt_run, sample)
        base_metrics, alt_metrics = metrics(base_view), metrics(alt_view)
        deltas = {
            "determined_primary_pp": alt_metrics["determined_primary_pct"] - base_metrics["determined_primary_pct"],
            "multiHP_pp": alt_metrics["multiHP_pct_all_regions"] - base_metrics["multiHP_pct_all_regions"],
            "all_determined_region_pp": alt_metrics["all_determined_pct_primary_regions"] - base_metrics["all_determined_pct_primary_regions"],
        }
        base_positions = retained_positions(args.base_run, sample)
        alt_positions = retained_positions(alt_run, sample)
        retained_jaccard = jaccard(base_positions, alt_positions)
        base_units = primary_units(args.base_run, sample)
        alt_units = primary_units(alt_run, sample)
        shared_units = set(base_units) & set(alt_units)
        unit_jaccard = jaccard(set(base_units), set(alt_units))
        topology_concordance = (
            sum(base_units[key] == alt_units[key] for key in shared_units) / len(shared_units)
            if shared_units else (1.0 if not base_units and not alt_units else 0.0)
        )
        rows.append({"label": label, "sample": sample, "base": base_metrics, "alternative": alt_metrics,
                     "delta": deltas, "retained_position_jaccard": retained_jaccard,
                     "primary_unit_key_jaccard": unit_jaccard,
                     "n_base_primary_units": len(base_units), "n_alternative_primary_units": len(alt_units),
                     "n_shared_primary_units": len(shared_units),
                     "shared_unit_topology_digest_concordance": topology_concordance,
                     "verdict": verdict(deltas.values(), retained_jaccard, unit_jaccard,
                                        topology_concordance)})
    maximum = max(max(abs(value) for value in row["delta"].values()) for row in rows)
    output = {"schema_version": "2.1", "scope": "chr1-22 x 7 datasets / 6 biological samples",
              "canonical": "LongPhase-S _sc.vcf FILTER=PASS",
              "alternative": "ClairS paired FILTER=PASS selected tree using the same raw-all LongPhase-S HP/PS tags",
              "decision_thresholds": {"robust": "max abs delta <5pp",
              "moderate": "5-10pp", "sensitive": ">=10pp",
              "retained_position_jaccard": "robust >=0.95; moderate >=0.80",
              "primary_unit_key_jaccard": "robust >=0.90; moderate >=0.70",
              "shared_topology_digest_concordance": "robust >=0.95; moderate >=0.80"},
              "aggregate": {"max_abs_delta_pp": maximum,
                            "min_retained_position_jaccard": min(row["retained_position_jaccard"] for row in rows),
                            "min_primary_unit_key_jaccard": min(row["primary_unit_key_jaccard"] for row in rows),
                            "min_shared_topology_digest_concordance": min(
                                row["shared_unit_topology_digest_concordance"] for row in rows),
                            "verdict": verdict(
                                [value for row in rows for value in row["delta"].values()],
                                min(row["retained_position_jaccard"] for row in rows),
                                min(row["primary_unit_key_jaccard"] for row in rows),
                                min(row["shared_unit_topology_digest_concordance"] for row in rows),
                            )},
              "comparisons": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with args.output.with_suffix(".tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["label", "sample", "retained_jaccard", "primary_unit_jaccard",
                         "shared_topology_digest_concordance", "determined_delta_pp", "multiHP_delta_pp",
                         "region_determined_delta_pp", "verdict"])
        for row in rows:
            writer.writerow([row["label"], row["sample"], row["retained_position_jaccard"],
                             row["primary_unit_key_jaccard"],
                             row["shared_unit_topology_digest_concordance"],
                             row["delta"]["determined_primary_pp"], row["delta"]["multiHP_pp"],
                             row["delta"]["all_determined_region_pp"], row["verdict"]])
    print(f"BACKBONE SENSITIVITY -> {args.output}")
    print(f"  aggregate: {output['aggregate']}")
    for row in rows:
        print(f"  {row['sample']}: siteJ={row['retained_position_jaccard']:.3f} "
              f"unitJ={row['primary_unit_key_jaccard']:.3f} "
              f"topology={row['shared_unit_topology_digest_concordance']:.3f} "
              f"{row['delta']} {row['verdict']}")


if __name__ == "__main__":
    main()
