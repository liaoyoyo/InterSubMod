#!/usr/bin/env python3
"""Derive an exact MINREAD=4 run from MINREAD=3 mlhp outputs without rereading BAM."""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE_FILES = (
    "derive_minread4_run.py",
    "layered_tree_reconstruction.py",
    "build_region_view.py",
    "build_ssnv_site_ledger.py",
    "verify_layered_v2.py",
)


def filtered(values, minimum=4):
    return {key: value for key, value in (values or {}).items() if value >= minimum}


def transform_group(group):
    group = dict(group)
    group["populations"] = filtered(group.get("populations"))
    group["subread_groups"] = filtered(group.get("subread_groups"))
    group["populations_by_hp"] = {
        family: kept for family, values in (group.get("populations_by_hp") or {}).items()
        if (kept := filtered(values))
    }
    group["subread_groups_by_hp"] = {
        family: kept for family, values in (group.get("subread_groups_by_hp") or {}).items()
        if (kept := filtered(values))
    }
    active = set(group["populations_by_hp"]) | set(group["subread_groups_by_hp"])
    group["col_coverage_by_hp"] = {
        family: values for family, values in (group.get("col_coverage_by_hp") or {}).items() if family in active
    }
    group["n_populations"] = len(group["populations"])
    group["n_populations_with_ALT"] = sum("A" in genotype for genotype in group["populations"])
    keep = group.get("n_full_cov_reads", 0) >= 4 or bool(group["subread_groups"])
    return group, keep


def aggregate(groups, allowed_cn=None):
    out = Counter()
    for group in groups:
        if allowed_cn is not None and group.get("cn") not in allowed_cn:
            continue
        n = group["n_populations"]
        out[f"npop_{min(n, 5)}"] += 1
        out[f"nsnv_{min(group['n_sSNV'], 6)}"] += 1
        if n >= 2:
            out["groups_ge2_pop"] += 1
        if n >= 3:
            out["groups_ge3_pop"] += 1
    return dict(out)


def transform_doc(doc):
    groups = []
    dropped_groups = dropped_snvs = 0
    for raw in doc.get("groups", []):
        group, keep = transform_group(raw)
        if keep:
            groups.append(group)
        else:
            dropped_groups += 1
            dropped_snvs += raw["n_sSNV"]
    funnel = dict(doc["input_funnel"])
    funnel["n_groups_read_unsupported"] += dropped_groups
    funnel["n_sSNV_read_unsupported"] += dropped_snvs
    funnel["n_groups_retained"] -= dropped_groups
    funnel["n_sSNV_retained"] -= dropped_snvs
    funnel["n_sSNV_accounted"] = (funnel["n_positional_singleton"] + funnel["n_sSNV_cap_excluded"]
                                     + funnel["n_sSNV_read_unsupported"] + funnel["n_sSNV_retained"])
    funnel["check_scope_conservation"] = funnel["n_sSNV_scope_input"] == funnel["n_sSNV_accounted"]
    params = dict(doc["params"])
    params["MINREAD"] = 4
    params["derived_from_MINREAD"] = 3
    params["derivation"] = "exact upward threshold filter; no BAM reread"
    tag_census = dict(doc.get("read_tag_census", {}))
    phase_set_region_counts = Counter()
    for group in groups:
        n_phase_sets = int(group.get("n_unique_phase_sets", 0))
        phase_set_region_counts["none" if n_phase_sets == 0 else ("one" if n_phase_sets == 1 else "multiple")] += 1
    tag_census["phase_set_region_counts"] = dict(phase_set_region_counts)
    tag_census["n_regions_with_phase_set_census"] = len(groups)
    doc = dict(doc)
    doc.update({"params": params, "input_funnel": funnel, "groups": groups,
                "read_tag_census": tag_census,
                "n_groups_analyzed": len(groups), "aggregate_all": aggregate(groups),
                "aggregate_clean_loh_neutral": aggregate(groups, {"neutral", "loh"}),
                "aggregate_strict_neutral": aggregate(groups, {"neutral"})})
    return doc


def run(command, env, log):
    proc = subprocess.run(command, env=env, text=True, capture_output=True, check=False)
    log.append({"command": command, "exit_code": proc.returncode,
                "stdout": proc.stdout[-8000:], "stderr": proc.stderr[-8000:]})
    if proc.returncode:
        raise RuntimeError(f"command failed: {' '.join(command)}\n{proc.stderr}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-run", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    verification_path = args.base_run / "verification_summary.json"
    if (manifest.get("schema_version") != "2.1"
            or manifest.get("tag_contract", {}).get("truth_flags") is not False
            or manifest.get("tag_contract", {}).get("tree_backbone") != "LongPhase-S _sc.vcf FILTER=PASS"
            or manifest.get("tree_input_contract") != "longphase_recalibrated_PASS"):
        raise SystemExit("MINREAD=4 derivation requires the clean no-truth schema 2.1 manifest")
    if not verification_path.is_file() or not (args.base_run / "_SUCCESS").is_file():
        raise SystemExit(f"base verification/_SUCCESS missing: {args.base_run}")
    base_verification = json.loads(verification_path.read_text(encoding="utf-8"))
    if base_verification.get("all_pass") is not True or base_verification.get("n_pass") != len(manifest["samples"]):
        raise SystemExit("MINREAD=4 derivation requires an all-pass clean base run")
    for result in base_verification.get("samples", []):
        sample = result.get("sample")
        exact_join = result.get("metrics", {}).get("read_tag_census", {}).get("check_exact_sidecar_join")
        if exact_join is None:
            layered_path = args.base_run / "samples" / sample / f"layered_reconstruction_{sample}.json"
            exact_join = json.loads(layered_path.read_text(encoding="utf-8")).get(
                "read_tag_census", {}).get("check_exact_sidecar_join")
        if exact_join is not True:
            raise SystemExit(f"{sample} base exact sidecar join is not PASS")
    for meta in manifest["samples"]:
        sample = meta["sample"]
        parts = sorted((args.base_run / "samples" / sample).glob("mlhp_part_*.json"))
        if len(parts) != 5:
            raise SystemExit(f"{sample} base run has {len(parts)} mlhp parts, expected 5")
        if any(json.loads(path.read_text(encoding="utf-8"))["params"].get("MINREAD") != 3 for path in parts):
            raise SystemExit(f"{sample} base run is not MINREAD=3")
    if args.output_root.exists():
        raise SystemExit(f"immutable output root exists: {args.output_root}")
    args.output_root.mkdir(parents=True)
    (args.output_root / "samples").mkdir()
    source_root = args.output_root / "source"
    source_root.mkdir()
    for name in SOURCE_FILES:
        shutil.copy2(HERE / name, source_root / name)
    shutil.copy2(args.input_manifest, args.output_root / "input_manifest.json")
    with (args.output_root / "source_bundle.sha256").open("w", encoding="utf-8") as handle:
        for path in sorted(source_root.iterdir()):
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n")
    (args.output_root / "params.json").write_text(json.dumps({
        "variant": "minread4", "MINREAD": 4, "derived_from_MINREAD": 3,
        "derivation": "exact upward threshold filter; no BAM reread",
        "base_run": str(args.base_run.resolve()),
    }, indent=2) + "\n", encoding="utf-8")
    python = sys.executable
    command_log = []
    try:
        for meta in manifest["samples"]:
            sample = meta["sample"]
            source_dir = args.base_run / "samples" / sample
            target_dir = args.output_root / "samples" / sample
            target_dir.mkdir()
            for path in sorted(source_dir.glob("mlhp_part_*.json")):
                transformed = transform_doc(json.loads(path.read_text(encoding="utf-8")))
                (target_dir / path.name).write_text(json.dumps(transformed, ensure_ascii=False) + "\n", encoding="utf-8")
            env = os.environ.copy()
            env.update({
                "SM_ML_GLOB": str(target_dir / "mlhp_part_*.json"),
                "SM_OUT": str(target_dir / f"layered_reconstruction_{sample}.json"),
                "SM_VERIFY_EVERY": "1",
                "SM_ANALYSIS_TREE_CAP": "0",
                "SM_DISPLAY_TREE_CAP": "32",
                "SM_CN_INT_GAIN": meta.get("cn_int_gain") or "",
                "SM_CN_INT_LOSS": meta.get("cn_int_loss") or "",
            })
            run([python, str(source_root / "layered_tree_reconstruction.py")], env, command_log)
            env.update({
                "SM_LAYERED": str(target_dir / f"layered_reconstruction_{sample}.json"),
                "SM_OUT": str(target_dir / f"layered_region_view_{sample}.json"),
                "SM_SAMPLE": sample,
                "SM_SOMATIC_VCF": meta["somatic_vcf"],
                "SM_INTEGRATION": meta.get("integration_json") or "",
                "SM_BACKBONE_SOURCE": meta.get("somatic_vcf_role") or "LongPhase-S recalibrated FILTER=PASS tree input",
            })
            run([python, str(source_root / "build_region_view.py")], env, command_log)
            run([python, str(source_root / "build_ssnv_site_ledger.py"), "--sample", sample,
                 "--caller-raw-vcf", meta["caller_raw_vcf"],
                 "--longphase-input-vcf", meta["longphase_input_vcf"],
                 "--tree-input-vcf", meta["somatic_vcf"],
                 "--recalibrated-vcf", meta["longphase_recalibrated_all_vcf"],
                 "--tree-contract", "longphase_recalibrated_PASS",
                 "--mlhp-glob", str(target_dir / "mlhp_part_*.json"),
                 "--output-tsv-gz", str(target_dir / f"ssnv_site_ledger_{sample}.tsv.gz"),
                 "--output-summary", str(target_dir / f"ssnv_site_ledger_{sample}.summary.json")],
                os.environ.copy(), command_log)
            print(f"{sample}: MINREAD=4 derived and reconstructed")
        run([python, str(source_root / "verify_layered_v2.py"), "--run-root", str(args.output_root),
             "--input-manifest", str(args.input_manifest), "--output", str(args.output_root / "verification_summary.json")],
            os.environ.copy(), command_log)
        verification = json.loads((args.output_root / "verification_summary.json").read_text(encoding="utf-8"))
        if verification.get("all_pass") is not True or verification.get("n_pass") != 7:
            raise RuntimeError("derived MINREAD=4 verification is not 7/7 all-pass")
        (args.output_root / "derive_commands.json").write_text(
            json.dumps(command_log, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        success_pending = args.output_root / f"_SUCCESS.pending.{os.getpid()}"
        success_pending.write_text(json.dumps({
            "status": "SUCCESS", "variant": "minread4", "dataset_count": 7,
            "verification_sha256": hashlib.sha256(
                (args.output_root / "verification_summary.json").read_bytes()).hexdigest(),
            "source_bundle_sha256": hashlib.sha256(
                (args.output_root / "source_bundle.sha256").read_bytes()).hexdigest(),
        }, indent=2) + "\n", encoding="utf-8")
        os.replace(success_pending, args.output_root / "_SUCCESS")
    except Exception as error:
        failed_pending = args.output_root / f"_FAILED.pending.{os.getpid()}"
        failed_pending.write_text(json.dumps({"status": "FAILED", "error": str(error)}, indent=2) + "\n",
                                  encoding="utf-8")
        os.replace(failed_pending, args.output_root / "_FAILED")
        raise


if __name__ == "__main__":
    main()
