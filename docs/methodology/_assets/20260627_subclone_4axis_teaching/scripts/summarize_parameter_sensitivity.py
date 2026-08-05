#!/usr/bin/env python3
"""Summarize all-dataset MAPQ/BaseQ/MINREAD/MAX_SNV sensitivity runs."""

import argparse
import csv
import json
import statistics
from pathlib import Path


VARIANTS = ["mapq30", "baseq10", "minread4", "maxsnv6"]


def load_view(run_root, sample):
    path = run_root / "samples" / sample / f"layered_region_view_{sample}.json"
    return json.loads(path.read_text(encoding="utf-8"))


def positions(run_root, sample):
    result = set()
    for path in (run_root / "samples" / sample).glob("mlhp_part_*.json"):
        for group in json.loads(path.read_text(encoding="utf-8")).get("groups", []):
            result.update((group["chrom"], int(pos)) for pos in group.get("positions", []))
    return result


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
            raise RuntimeError(f"missing candidate-tree digest for {sample}: {key}")
        values[key] = digest
    return values


def metrics(view):
    census = view["census"]
    l1 = census["L1"]
    primary = l1["n_primary_lineage_units"] or 1
    determined = l1["determinacy_primary_lineage"].get("determined", 0)
    regions = census["n_regions"] or 1
    with_primary = census["n_regions_with_primary_lineage"] or 1
    multi = census["hp_multiplicity"].get("2", 0)
    all_det = census["region_determinacy"].get("all_determined", 0)
    return {"retained_sSNV": census["funnel"]["L6_retained_sSNV"], "regions": census["n_regions"],
            "primary": l1["n_primary_lineage_units"], "determined_pct": 100 * determined / primary,
            "multiHP_pct": 100 * multi / regions, "region_determined_pct": 100 * all_det / with_primary}


def jaccard(a, b):
    return len(a & b) / len(a | b) if a or b else 1.0


def verdict(maximum, site_jaccard, unit_jaccard, topology_concordance):
    if (maximum < 5 and site_jaccard >= 0.95 and unit_jaccard >= 0.90
            and topology_concordance >= 0.95):
        return "robust_all_gates"
    if (maximum < 10 and site_jaccard >= 0.80 and unit_jaccard >= 0.70
            and topology_concordance >= 0.80):
        return "moderately_sensitive"
    return "sensitive"


def validate_run(run_root, expected_samples, expected_contract="longphase_recalibrated_PASS"):
    required = [run_root / "_SUCCESS", run_root / "verification_summary.json"]
    if any(not path.is_file() for path in required):
        raise SystemExit(f"incomplete sensitivity run: {run_root}")
    verification = json.loads(required[1].read_text(encoding="utf-8"))
    verified_samples = [item.get("sample") for item in verification.get("samples", [])]
    if (verification.get("all_pass") is not True or verification.get("n_pass") != 7
            or len(verified_samples) != 7 or set(verified_samples) != set(expected_samples)):
        raise SystemExit(f"sensitivity run contract/sample set failed: {run_root}")
    manifest_path = run_root / "input_manifest.json"
    lock_path = run_root / "frozen_input_lock.json"
    if manifest_path.is_file():
        frozen = json.loads(manifest_path.read_text(encoding="utf-8"))
        frozen_samples = [item.get("sample") for item in frozen.get("samples", [])]
        if frozen_samples != expected_samples or frozen.get("tree_input_contract") != expected_contract:
            raise SystemExit(f"v2 sensitivity frozen manifest failed: {run_root}")
    elif lock_path.is_file() and expected_contract == "longphase_recalibrated_PASS":
        frozen = json.loads(lock_path.read_text(encoding="utf-8"))
        frozen_samples = [item.get("sample") for item in frozen.get("samples", [])]
        canonical_roles = all(
            set(item.get("somatic", {}))
            == {"caller_raw_vcf", "longphase_input_vcf", "longphase_recalibrated_all_vcf", "tree_vcf"}
            and item["somatic"]["tree_vcf"].get("path")
            != item["somatic"]["longphase_recalibrated_all_vcf"].get("path")
            for item in frozen.get("samples", [])
        )
        if (len(frozen_samples) != 7 or set(frozen_samples) != set(expected_samples)
                or frozen.get("all_pass") is not True or not canonical_roles
                or frozen.get("analysis_contract", {}).get("production_filter_policy", {}).get("name")
                != "production_default_filter"):
            raise SystemExit(f"v3 canonical frozen lock failed: {run_root}")
    else:
        raise SystemExit(f"recognized frozen manifest/lock missing: {run_root}")
    for sample in expected_samples:
        sample_root = run_root / "samples" / sample
        layered_path = sample_root / f"layered_reconstruction_{sample}.json"
        site_path = sample_root / f"ssnv_site_ledger_{sample}.summary.json"
        if not layered_path.is_file() or not site_path.is_file():
            raise SystemExit(f"{sample} sensitivity scientific outputs are incomplete: {run_root}")
        layered = json.loads(layered_path.read_text(encoding="utf-8"))
        tags = layered.get("read_tag_census", {})
        phase_counts = tags.get("phase_set_region_counts", {})
        site = json.loads(site_path.read_text(encoding="utf-8"))
        if (tags.get("check_exact_sidecar_join") is not True
                or sum(phase_counts.values()) != tags.get("n_regions_with_phase_set_census")
                or layered.get("analysis_contract", {}).get("PS")
                != "preserved for phase-block QC; not used as a topology edge or lineage label"
                or site.get("pass") is not True):
            raise SystemExit(f"{sample} sensitivity science/tag/site-ledger contract failed: {run_root}")


def validate_variant_params(run_root, variant):
    path = run_root / "params.json"
    if path.is_file():
        document = json.loads(path.read_text(encoding="utf-8"))
        params = document.get("params", document)
    else:
        receipt_path = run_root / "launch_receipt.json"
        if not receipt_path.is_file():
            raise SystemExit(f"variant params/launch receipt missing: {run_root}")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        params = receipt.get("extra", {}).get("analysis_params", {})
    expected = {
        "mapq30": ("MAPQ_MIN", 30),
        "baseq10": ("BASEQ_MIN", 10),
        "minread4": ("MINREAD", 4),
        "maxsnv6": ("MAX_SNV", 6),
    }[variant]
    if params.get(expected[0]) != expected[1]:
        raise SystemExit(f"{variant} run does not declare {expected[0]}={expected[1]}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-run", required=True, type=Path)
    parser.add_argument("--parameter-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    expected_samples = [item["sample"] for item in manifest.get("samples", [])]
    if (manifest.get("schema_version") != "2.1"
            or manifest.get("tree_input_contract") != "longphase_recalibrated_PASS"
            or len(expected_samples) != 7 or len(set(expected_samples)) != 7):
        raise SystemExit("parameter sensitivity requires the full canonical clean manifest")
    validate_run(args.base_run, expected_samples)
    rows = []
    for variant in VARIANTS:
        variant_root = args.parameter_root / variant / "run"
        validate_run(variant_root, expected_samples)
        validate_variant_params(variant_root, variant)
        for meta in manifest["samples"]:
            sample = meta["sample"]
            base = metrics(load_view(args.base_run, sample))
            alt = metrics(load_view(variant_root, sample))
            delta = {"determined_pp": alt["determined_pct"] - base["determined_pct"],
                     "multiHP_pp": alt["multiHP_pct"] - base["multiHP_pct"],
                     "region_determined_pp": alt["region_determined_pct"] - base["region_determined_pct"]}
            site_jaccard = jaccard(positions(args.base_run, sample), positions(variant_root, sample))
            base_units = primary_units(args.base_run, sample)
            alt_units = primary_units(variant_root, sample)
            shared = set(base_units) & set(alt_units)
            unit_jaccard = jaccard(set(base_units), set(alt_units))
            topology = (sum(base_units[key] == alt_units[key] for key in shared) / len(shared)
                        if shared else (1.0 if not base_units and not alt_units else 0.0))
            rows.append({"variant": variant, "sample": sample, "base": base, "alternative": alt, "delta": delta,
                         "retained_position_jaccard": site_jaccard,
                         "primary_unit_key_jaccard": unit_jaccard,
                         "shared_unit_topology_digest_concordance": topology,
                         "n_shared_primary_units": len(shared)})
    aggregate = {}
    for variant in VARIANTS:
        subset = [row for row in rows if row["variant"] == variant]
        abs_deltas = [max(abs(value) for value in row["delta"].values()) for row in subset]
        maximum = max(abs_deltas)
        site_jaccard = min(row["retained_position_jaccard"] for row in subset)
        unit_jaccard = min(row["primary_unit_key_jaccard"] for row in subset)
        topology = min(row["shared_unit_topology_digest_concordance"] for row in subset)
        aggregate[variant] = {"max_abs_delta_pp": maximum,
                              "median_abs_delta_pp": statistics.median(abs_deltas),
                              "min_retained_position_jaccard": site_jaccard,
                              "min_primary_unit_key_jaccard": unit_jaccard,
                              "min_shared_topology_digest_concordance": topology,
                              "verdict": verdict(maximum, site_jaccard, unit_jaccard, topology)}
    output = {"schema_version": "2.0", "scope": "7 datasets / chr1-22",
              "thresholds": {"robust": "<5pp + site J>=0.95 + unit J>=0.90 + topology>=0.95",
                             "moderate": "<10pp + site J>=0.80 + unit J>=0.70 + topology>=0.80",
                             "sensitive": "otherwise"},
              "aggregate": aggregate, "rows": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with args.output.with_suffix(".tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant", "sample", "retained_jaccard", "primary_unit_jaccard",
                         "shared_topology_digest_concordance", "determined_delta_pp", "multiHP_delta_pp",
                         "region_determined_delta_pp"])
        for row in rows:
            writer.writerow([row["variant"], row["sample"], row["retained_position_jaccard"],
                             row["primary_unit_key_jaccard"],
                             row["shared_unit_topology_digest_concordance"],
                             row["delta"]["determined_pp"], row["delta"]["multiHP_pp"], row["delta"]["region_determined_pp"]])
    print(f"PARAMETER SENSITIVITY -> {args.output}")
    for variant, summary in aggregate.items():
        print(f"  {variant}: {summary}")


if __name__ == "__main__":
    main()
