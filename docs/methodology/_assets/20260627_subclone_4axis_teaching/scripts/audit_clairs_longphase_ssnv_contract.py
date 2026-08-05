#!/usr/bin/env python3
"""Audit ClairS -> LongPhase-S -> sSNV site and read-tag contracts."""

import argparse
import bisect
import csv
import json
import re
import shlex
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam


EXPECTED_HP_TAGS = [".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"]
AUTOSOMES = {f"chr{chrom}" for chrom in range(1, 23)}


def record_key(record):
    return (record.contig, int(record.pos), record.ref, tuple(record.alts or ()))


def is_pass(record):
    return "PASS" in set(record.filter.keys())


def is_biallelic_snv(record):
    return len(record.ref) == 1 and len(record.alts or ()) == 1 and len(record.alts[0]) == 1


def scan_vcf(path):
    counts = Counter({
        name: 0 for name in (
            "records", "autosomal_records", "non_autosomal_records", "biallelic_snvs",
            "pass_records", "pass_autosomal_records", "pass_non_autosomal_records",
            "pass_biallelic_snvs",
        )
    })
    keys = Counter()
    pass_keys = Counter()
    filters = Counter()
    headers = []
    with pysam.VariantFile(str(path)) as vcf:
        headers = str(vcf.header).splitlines()
        for record in vcf:
            key = record_key(record)
            keys[key] += 1
            counts["records"] += 1
            scope = "autosomal" if record.contig in AUTOSOMES else "non_autosomal"
            counts[f"{scope}_records"] += 1
            if is_biallelic_snv(record):
                counts["biallelic_snvs"] += 1
            labels = list(record.filter.keys()) or ["."]
            for label in labels:
                filters[label] += 1
            if is_pass(record):
                pass_keys[key] += 1
                counts["pass_records"] += 1
                counts[f"pass_{scope}_records"] += 1
                if is_biallelic_snv(record):
                    counts["pass_biallelic_snvs"] += 1
    return {
        "path": str(path),
        "counts": dict(counts),
        "filters": dict(filters),
        "keys": keys,
        "pass_keys": pass_keys,
        "source": next((line.split("=", 1)[1] for line in headers if line.startswith("##source=")), None),
        "version": next((line.split("=", 1)[1] for line in headers if line.startswith("##clairs_version=")), None),
    }


def counter_diff(left, right):
    missing = left - right
    extra = right - left
    return {
        "equal": not missing and not extra,
        "missing_count": sum(missing.values()),
        "extra_count": sum(extra.values()),
        "missing_examples": [list(key) for key in list(missing)[:5]],
        "extra_examples": [list(key) for key in list(extra)[:5]],
    }


def longphase_program(bam_path):
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        programs = bam.header.to_dict().get("PG", [])
    candidates = [item for item in programs if "somatic_haplotag" in item.get("CL", "")]
    if not candidates:
        return {"program": None, "command": None, "args": []}
    item = candidates[-1]
    command = item.get("CL", "")
    try:
        args = shlex.split(command)
    except ValueError:
        args = command.split()
    return {"program": item, "command": command, "args": args}


def option(args, name):
    if name not in args:
        return None
    index = args.index(name)
    return args[index + 1] if index + 1 < len(args) else True


def load_bed(path):
    intervals = defaultdict(list)
    if not path or not Path(path).is_file():
        return intervals
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            intervals[fields[0]].append((start, end))
    merged = {}
    for chrom, values in intervals.items():
        output = []
        for start, end in sorted(values):
            if output and start <= output[-1][1]:
                output[-1] = (output[-1][0], max(output[-1][1], end))
            else:
                output.append((start, end))
        merged[chrom] = {
            "starts": [value[0] for value in output],
            "ends": [value[1] for value in output],
        }
    return merged


def in_bed(intervals, chrom, pos1):
    values = intervals.get(chrom)
    if not values:
        return False
    pos0 = pos1 - 1
    index = bisect.bisect_right(values["starts"], pos0) - 1
    return index >= 0 and pos0 < values["ends"][index]


def bed_overlap(path, bed_path):
    intervals = load_bed(bed_path)
    if not intervals:
        return None
    total = overlap = 0
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            total += 1
            overlap += int(in_bed(intervals, record.contig, record.pos))
    return {
        "total": total,
        "inside_truth_bed": overlap,
        "outside_truth_bed": total - overlap,
        "inside_fraction": overlap / total if total else None,
    }


def parse_log(path):
    text = path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""
    result = {"path": str(path)}
    patterns = {
        "tumor_snp_count": r"Tumor SNP count:\s*(\d+)",
        "somatic_variant_flag_count": r"somatic variant count\(Flag\):\s*(\d+)",
        "total_alignment": r"total alignment\s*:\s*(\d+)",
        "total_tagged_alignments": r"total tagged alignments\s*:\s*(\d+)",
        "total_untagged": r"total untagged\s*:\s*(\d+)",
    }
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        result[name] = int(match.group(1)) if match else None
    result["truth_bed_removal_executed"] = "removing tumor & truth somatic variants outside bed regions" in text
    return result


def consumer_contract(script_dir):
    source = (script_dir / "sm_multilocus_combinations.py").read_text(encoding="utf-8")
    sys.path.insert(0, str(script_dir))
    from sm_multilocus_combinations import germ_family
    mapping = {tag: germ_family(tag) for tag in EXPECTED_HP_TAGS}
    return {
        "expected_longphase_hp_vocabulary": EXPECTED_HP_TAGS,
        "current_family_mapping": mapping,
        "H4_is_explicit": mapping["4"] == "4",
        "PS_is_persisted_in_mlhp_output": "n_unique_phase_sets" in source and "raw_HP_with_PS_counts" in source,
        "fine_hp_counts_are_persisted": '"raw_hp' in source.lower(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    rows = []
    for meta in manifest["samples"]:
        sample = meta["sample"]
        tagged_bam = Path(meta["tumor_bam"])
        longphase_dir = tagged_bam.parent
        run_root = longphase_dir.parent
        context_path = run_root / "run_context.json"
        context = json.loads(context_path.read_text(encoding="utf-8"))
        raw_path = Path(context["somatic_vcf"])
        input_path = Path(meta["somatic_vcf"])
        sc_path = longphase_dir / f"{sample}_tagged_sc.vcf"
        raw = scan_vcf(raw_path)
        caller_input = scan_vcf(input_path)
        recalibrated = scan_vcf(sc_path)
        pg = longphase_program(tagged_bam)
        truth_bed = option(pg["args"], "--truth-bed")
        truth_vcf = option(pg["args"], "--truth-vcf")
        log = parse_log(longphase_dir / "longphase_s.log")
        raw_to_input = counter_diff(raw["pass_keys"], caller_input["keys"])
        input_to_output = counter_diff(caller_input["keys"], recalibrated["keys"])
        row = {
            "sample": sample,
            "biological_id": meta["biological_id"],
            "raw_clairs": {key: value for key, value in raw.items() if key not in {"keys", "pass_keys"}},
            "longphase_input": {key: value for key, value in caller_input.items() if key not in {"keys", "pass_keys"}},
            "longphase_recalibrated_output": {key: value for key, value in recalibrated.items() if key not in {"keys", "pass_keys"}},
            "raw_PASS_to_longphase_input": raw_to_input,
            "longphase_input_to_recalibrated_all": input_to_output,
            "longphase_parser_count_matches_input": log["tumor_snp_count"] == caller_input["counts"]["biallelic_snvs"],
            "longphase_log": log,
            "tagging_command": pg["command"],
            "tag_supplementary": "--tagSupplementary" in pg["args"],
            "mapq": option(pg["args"], "-q"),
            "truth_vcf": truth_vcf,
            "truth_bed": truth_bed,
            "truth_bed_scope": bed_overlap(input_path, truth_bed),
            "production_genomewide_tagging_contract": not bool(truth_bed),
        }
        row["pass"] = all([
            raw_to_input["equal"],
            input_to_output["equal"],
            row["longphase_parser_count_matches_input"],
            row["production_genomewide_tagging_contract"],
            row["tag_supplementary"],
            str(row["mapq"]) == "20",
        ])
        rows.append(row)
        print(f"{sample}: rawPASS->input={raw_to_input['equal']} input->sc={input_to_output['equal']} "
              f"parser={row['longphase_parser_count_matches_input']} production_scope={row['production_genomewide_tagging_contract']}")
    consumer = consumer_contract(Path(__file__).resolve().parent)
    output = {
        "schema_version": "1.0",
        "task_type": "B_comprehensive_validation",
        "claim": "All ClairS PASS records must reach LongPhase-S and all sSNV-consumed reads must use an explicit HP/PS contract",
        "dataset_count": len(rows),
        "n_pass": sum(row["pass"] for row in rows),
        "all_pass": all(row["pass"] for row in rows) and consumer["H4_is_explicit"]
                    and consumer["PS_is_persisted_in_mlhp_output"] and consumer["fine_hp_counts_are_persisted"],
        "consumer_contract": consumer,
        "samples": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with args.output.with_suffix(".tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "raw_records", "raw_PASS", "raw_non_autosomal", "raw_PASS_non_autosomal",
                         "lps_input", "lps_input_non_autosomal", "lps_sc_all", "lps_sc_PASS",
                         "rawPASS_input_equal", "input_sc_equal", "truth_bed", "inside_bed", "outside_bed",
                         "production_scope", "pass"])
        for row in rows:
            scope = row["truth_bed_scope"] or {}
            writer.writerow([
                row["sample"], row["raw_clairs"]["counts"]["records"],
                row["raw_clairs"]["counts"]["pass_records"],
                row["raw_clairs"]["counts"]["non_autosomal_records"],
                row["raw_clairs"]["counts"]["pass_non_autosomal_records"],
                row["longphase_input"]["counts"]["records"],
                row["longphase_input"]["counts"]["non_autosomal_records"],
                row["longphase_recalibrated_output"]["counts"]["records"],
                row["longphase_recalibrated_output"]["counts"]["pass_records"],
                row["raw_PASS_to_longphase_input"]["equal"], row["longphase_input_to_recalibrated_all"]["equal"],
                row["truth_bed"] or "", scope.get("inside_truth_bed", ""), scope.get("outside_truth_bed", ""),
                row["production_genomewide_tagging_contract"], row["pass"],
            ])
    print(f"AUDIT -> {args.output} ({output['n_pass']}/{output['dataset_count']} sample contracts pass; all_pass={output['all_pass']})")
    raise SystemExit(0 if output["all_pass"] else 1)


if __name__ == "__main__":
    main()
