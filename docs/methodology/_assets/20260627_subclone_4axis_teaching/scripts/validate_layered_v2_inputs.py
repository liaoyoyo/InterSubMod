#!/usr/bin/env python3
"""Preflight and lock the seven-dataset layered v2 input manifest."""

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import pysam


AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
CN_STATES = {"gain", "loss", "loh", "neutral"}
EXPECTED_BINDING = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def digest_json(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def inspect_vcf(path):
    header = []
    total = autosomal = out_scope = naf0 = 0
    filters = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                header.append(line.rstrip())
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10 or len(fields[3]) != 1 or len(fields[4]) != 1:
                continue
            total += 1
            filters.add(fields[6])
            if fields[0] in AUTOSOMES:
                autosomal += 1
            else:
                out_scope += 1
            fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
            try:
                naf0 += float(fmt.get("NAF", "nan")) == 0.0
            except ValueError:
                pass
    source = [x for x in header if x.startswith("##source=")]
    version = [x for x in header if x.startswith("##clairs_version=")]
    cmdline = [x for x in header if "cmdline" in x.lower()]
    return {
        "sha256": sha256(path),
        "snv_total": total,
        "autosomal_chr1_22": autosomal,
        "out_of_scope": out_scope,
        "filters": sorted(filters),
        "NAF_zero": naf0,
        "NAF_zero_fraction": round(naf0 / total, 6) if total else None,
        "source_header": source,
        "version_header": version,
        "cmdline_header": cmdline[:2],
    }


def vcf_key_counters(path):
    all_keys = Counter()
    pass_keys = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            record_key = (record.contig, int(record.pos), record.ref, tuple(record.alts or ()))
            all_keys[record_key] += 1
            if set(record.filter.keys()) == {"PASS"}:
                pass_keys[record_key] += 1
    return all_keys, pass_keys


def inspect_vcf_flow(raw_vcf, longphase_input, recalibrated_all, tree_input, tree_contract):
    raw_all, raw_pass = vcf_key_counters(raw_vcf)
    lps_all, _ = vcf_key_counters(longphase_input)
    recal_all, recal_pass = vcf_key_counters(recalibrated_all)
    tree_all, tree_pass = vcf_key_counters(tree_input)
    checks = {
        "raw_PASS_equals_longphase_input": raw_pass == lps_all,
        "longphase_input_equals_recalibrated_all": lps_all == recal_all,
        "tree_input_all_FILTER_PASS": tree_all == tree_pass,
    }
    if tree_contract == "longphase_recalibrated_PASS":
        checks["recalibrated_PASS_equals_tree_input"] = recal_pass == tree_all
    elif tree_contract == "clairs_PASS_input":
        checks["longphase_input_equals_tree_input"] = lps_all == tree_all
    else:
        checks["recognized_tree_contract"] = False
    return {
        "counts": {
            "raw_all": sum(raw_all.values()),
            "raw_PASS": sum(raw_pass.values()),
            "longphase_input": sum(lps_all.values()),
            "recalibrated_all": sum(recal_all.values()),
            "recalibrated_PASS": sum(recal_pass.values()),
            "tree_input": sum(tree_all.values()),
        },
        "checks": checks,
        "pass": all(checks.values()),
    }


def inspect_bam(path):
    with pysam.AlignmentFile(str(path), "rb") as bam:
        header = bam.header.to_dict()
        references = list(zip(bam.references, bam.lengths))
    index = Path(str(path) + ".bai")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
        "header_sha256": digest_json(header),
        "reference_dictionary_sha256": digest_json(references),
        "index_path": str(index),
        "index_size_bytes": index.stat().st_size,
        "index_sha256": sha256(index),
    }


def inspect_cn(path):
    states = set()
    n = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 4:
            continue
        int(fields[1]); int(fields[2])
        states.add(fields[3].lower())
        n += 1
    return {"sha256": sha256(path), "n_intervals": n, "states": sorted(states),
            "states_recognized": states <= CN_STATES}


def inspect_production_inventory(path, tumor_bam):
    with path.open(encoding="utf-8") as handle:
        rows = {row["role"]: row for row in csv.DictReader(handle, delimiter="\t")}
    recorded = rows.get("tumor_bam")
    errors = []
    if not recorded:
        return {"path": str(path), "match": False, "errors": ["tumor_bam row missing"]}
    if Path(recorded["path"]).resolve() != tumor_bam.resolve():
        errors.append("tumor BAM path differs from production inventory")
    logical_stat = tumor_bam.lstat()
    if int(recorded["size_bytes"]) != logical_stat.st_size:
        errors.append("tumor BAM size differs from production inventory")
    try:
        date_part, zone = recorded["mtime"].rsplit(" ", 1)
        whole, fraction = date_part.rsplit(".", 1)
        seconds = int(dt.datetime.strptime(f"{whole} {zone}", "%Y-%m-%d %H:%M:%S %z").timestamp())
        recorded_mtime_ns = seconds * 1_000_000_000 + int(fraction[:9].ljust(9, "0"))
        if recorded_mtime_ns != logical_stat.st_mtime_ns:
            errors.append("tumor BAM mtime differs from production inventory")
    except (ValueError, TypeError):
        errors.append("tumor BAM production mtime is not parseable")
    return {"path": str(path), "sha256": sha256(path), "match": not errors,
            "recorded_tumor_bam": recorded, "errors": errors}


def inspect_bam_fingerprint(expected, tumor_bam):
    logical = tumor_bam.lstat()
    resolved = tumor_bam.resolve()
    target = resolved.stat()
    observed = {
        "logical_path": str(tumor_bam),
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "resolved_path": str(resolved),
        "resolved_size_bytes": target.st_size,
        "resolved_mtime_ns": target.st_mtime_ns,
    }
    required = set(observed)
    errors = [f"{key} missing" for key in sorted(required - set(expected))]
    errors.extend(f"{key} changed" for key, value in expected.items() if observed.get(key) != value)
    return {"match": not errors, "expected": expected, "observed": observed, "errors": errors}


def validate_production_closeout(contract, samples):
    fields = {
        "closeout": ("production_closeout", "production_closeout_sha256"),
        "success": ("production_success", "production_success_sha256"),
        "artifacts": ("production_artifacts_manifest", "production_artifacts_manifest_sha256"),
    }
    paths = {}
    for role, (path_field, hash_field) in fields.items():
        path = Path(contract.get(path_field, ""))
        if not path.is_file() or contract.get(hash_field) != sha256(path):
            raise SystemExit(f"production {role} path/hash contract failed: {path}")
        paths[role] = path
    closeout = json.loads(paths["closeout"].read_text(encoding="utf-8"))
    success = json.loads(paths["success"].read_text(encoding="utf-8"))
    if (closeout.get("status") != "PASS" or closeout.get("dataset_count") != 7
            or closeout.get("n_pass") != 7 or closeout.get("truth_flags") is not False
            or closeout.get("tree_backbone") != "LongPhase-S _sc.vcf FILTER=PASS"):
        raise SystemExit("production closeout does not prove the canonical seven-dataset contract")
    if (success.get("status") != "SUCCESS"
            or Path(success.get("closeout_receipt", "")).resolve() != paths["closeout"].resolve()
            or success.get("closeout_receipt_sha256") != sha256(paths["closeout"])
            or Path(success.get("artifacts_manifest", "")).resolve() != paths["artifacts"].resolve()
            or success.get("artifacts_manifest_sha256") != sha256(paths["artifacts"])):
        raise SystemExit("production _SUCCESS receipt does not bind the declared closeout")
    if any(Path(item.get("longphase_production_closeout", "")).resolve() != paths["closeout"].resolve()
           for item in samples):
        raise SystemExit("one or more samples are not bound to the production closeout")
    return {role: str(path) for role, path in paths.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    clean_tag_contract = manifest.get("schema_version") == "2.1"
    if not clean_tag_contract:
        raise SystemExit("clean layered validation requires manifest schema_version 2.1")
    contract = manifest.get("tag_contract", {})
    tree_contract = manifest.get("tree_input_contract")
    canonical_contract = (
        manifest.get("task_type") == "B_COMPREHENSIVE_VALIDATION"
        and
        tree_contract == "longphase_recalibrated_PASS"
        and contract.get("tree_backbone") == "LongPhase-S _sc.vcf FILTER=PASS"
    )
    sensitivity_contract = (
        manifest.get("task_type") == "B_BACKBONE_SENSITIVITY"
        and tree_contract == "clairs_PASS_input"
        and contract.get("tree_backbone") == "ClairS PASS sensitivity"
    )
    if (contract.get("truth_flags") is not False or contract.get("PS_preserved") is not True
            or contract.get("bam_identity_locked") is not True
            or contract.get("longphase_filtering_policy") != "production_default_filter"
            or not (canonical_contract or sensitivity_contract)):
        raise SystemExit(
            "manifest must declare no-truth exact tags and a recognized canonical/sensitivity tree contract"
        )
    samples = manifest.get("samples", [])
    observed_binding = {item.get("sample"): item.get("biological_id") for item in samples}
    if (manifest.get("dataset_count") != 7 or manifest.get("biological_sample_count") != 6
            or len(samples) != 7 or len(observed_binding) != 7 or observed_binding != EXPECTED_BINDING):
        raise SystemExit(f"manifest does not match the canonical 7-dataset/6-biological binding: {observed_binding}")
    closeout_info = validate_production_closeout(contract, samples)
    results = []
    for meta in samples:
        errors = []
        bam = Path(meta["tumor_bam"])
        vcf = Path(meta["somatic_vcf"])
        for label, path in (("tumor_bam", bam), ("tumor_bam_index", Path(str(bam) + ".bai")), ("somatic_vcf", vcf)):
            if not path.is_file():
                errors.append(f"{label} missing: {path}")
        vcf_indexes = [Path(str(vcf) + suffix) for suffix in (".csi", ".tbi")]
        if not any(path.is_file() for path in vcf_indexes):
            errors.append(f"somatic_vcf index missing: {vcf}")
        cn_path = Path(meta["cn_bed"]) if meta.get("cn_bed") else None
        if cn_path and not cn_path.is_file():
            errors.append(f"cn_bed missing: {cn_path}")
        sidecar = Path(meta["read_tag_sidecar"]) if meta.get("read_tag_sidecar") else None
        sidecar_index = Path(meta["read_tag_sidecar_index"]) if meta.get("read_tag_sidecar_index") else None
        tag_validation = Path(meta["read_tag_validation"]) if meta.get("read_tag_validation") else None
        input_inventory = Path(meta["longphase_input_inventory"]) if meta.get("longphase_input_inventory") else None
        raw_vcf = Path(meta["caller_raw_vcf"]) if meta.get("caller_raw_vcf") else None
        longphase_input = Path(meta["longphase_input_vcf"]) if meta.get("longphase_input_vcf") else None
        recal_vcf = Path(meta["longphase_recalibrated_all_vcf"]) if meta.get("longphase_recalibrated_all_vcf") else None
        for label, path in (("read_tag_sidecar", sidecar), ("read_tag_sidecar_index", sidecar_index),
                            ("read_tag_validation", tag_validation), ("caller_raw_vcf", raw_vcf),
                            ("longphase_input_vcf", longphase_input),
                            ("longphase_recalibrated_all_vcf", recal_vcf),
                            ("longphase_input_inventory", input_inventory)):
            if path is None or not path.is_file():
                errors.append(f"{label} missing: {path}")
        if errors:
            results.append({"sample": meta["sample"], "pass": False, "errors": errors})
            continue
        vcf_info = inspect_vcf(vcf)
        cn_info = inspect_cn(cn_path) if cn_path else {"status": "unavailable"}
        tag_info = None
        inventory_info = None
        flow_info = inspect_vcf_flow(raw_vcf, longphase_input, recal_vcf, vcf, tree_contract)
        if not flow_info["pass"]:
            errors.append(f"ClairS/LongPhase-S/tree VCF flow failed: {flow_info['checks']}")
        if sidecar and sidecar_index and tag_validation:
            validation = json.loads(tag_validation.read_text(encoding="utf-8"))
            try:
                with pysam.TabixFile(str(sidecar)) as tags:
                    probe_contig = next(iter(tags.contigs), None)
                    first_tag_record = next(tags.fetch(probe_contig), None) if probe_contig else None
            except (OSError, ValueError) as error:
                probe_contig = None
                first_tag_record = None
                errors.append(f"read tag sidecar tabix error: {error}")
            if first_tag_record is None:
                errors.append("read tag sidecar contains no fetchable records")
            if validation.get("pass") is not True:
                errors.append("production read tag validation is not PASS")
            if Path(validation.get("sidecar", "")).resolve() != Path(str(sidecar).removesuffix(".gz")).resolve():
                errors.append("production read tag validation is bound to a different sidecar")
            if validation.get("checks", {}).get("truth_flags_absent") is not True:
                errors.append("production read tag source contains truth flags")
            if validation.get("region") != "all":
                errors.append(f"production read tag validation region is {validation.get('region')}")
            if validation.get("duplicate_exact_alignment_rows") != 0:
                errors.append("exact alignment identity is not one-to-one")
            required_checks = (
                "parser_count_matches_input",
                "sidecar_row_count_matches_capture",
                "tagged_count_matches_execution",
                "sidecar_no_unknown_HP",
                "sidecar_no_exact_identity_conflicts",
                "recalibrated_preserves_all_input_keys",
            )
            failed_checks = [name for name in required_checks if validation.get("checks", {}).get(name) is not True]
            if failed_checks:
                errors.append(f"production read tag validation checks failed: {failed_checks}")
            tag_info = {
                "path": str(sidecar), "size_bytes": sidecar.stat().st_size,
                "index_path": str(sidecar_index), "index_sha256": sha256(sidecar_index),
                "validation_path": str(tag_validation), "validation_sha256": sha256(tag_validation),
                "validation_pass": validation.get("pass"), "probe_contig": probe_contig,
                "probe_fetch_nonempty": first_tag_record is not None,
                "HP_counts": validation.get("HP_counts", {}),
            }
        if input_inventory:
            inventory_info = inspect_production_inventory(input_inventory, bam)
            errors.extend(inventory_info["errors"])
        fingerprint_info = inspect_bam_fingerprint(meta.get("tumor_bam_fingerprint", {}), bam)
        errors.extend(fingerprint_info["errors"])
        if vcf_info["filters"] != ["PASS"]:
            errors.append(f"VCF FILTER set is {vcf_info['filters']}")
        if not any("ClairS" in x for x in vcf_info["source_header"]):
            errors.append("VCF source header is not ClairS")
        if cn_path and not cn_info["states_recognized"]:
            errors.append(f"unrecognized CN states: {cn_info['states']}")
        results.append({
            "sample": meta["sample"],
            "biological_id": meta["biological_id"],
            "pass": not errors,
            "errors": errors,
            "bam": inspect_bam(bam),
            "vcf": vcf_info,
            "vcf_flow": flow_info,
            "vcf_index": next({"path": str(p), "sha256": sha256(p), "size_bytes": p.stat().st_size}
                              for p in vcf_indexes if p.is_file()),
            "cn_source": meta["cn_source"],
            "cn": cn_info,
            "read_tags": tag_info,
            "longphase_input_inventory": inventory_info,
            "tumor_bam_fingerprint": fingerprint_info,
        })
    output = {
        "schema_version": "2.1" if any(x.get("read_tags") for x in results) else "2.0",
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": sha256(args.manifest),
        "dataset_count": len(results),
        "biological_sample_count": len({x.get("biological_id") for x in results}),
        "production_closeout": closeout_info,
        "all_pass": all(x["pass"] for x in results),
        "samples": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"INPUT PREFLIGHT: {sum(x['pass'] for x in results)}/{len(results)} pass -> {args.output}")
    for result in results:
        print(f"  {result['sample']}: {'PASS' if result['pass'] else 'FAIL'}"
              + (f" ({'; '.join(result['errors'])})" if result["errors"] else ""))
    raise SystemExit(0 if output["all_pass"] else 1)


if __name__ == "__main__":
    main()
