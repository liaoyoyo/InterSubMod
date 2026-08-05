#!/usr/bin/env python3
"""Build a lossless per-site ClairS/LongPhase-S/sSNV disposition ledger."""

import argparse
import csv
import glob
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam


AUTOSOMES = {f"chr{value}" for value in range(1, 23)}


def scalar(value):
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    if isinstance(value, tuple):
        return [scalar(item) for item in value]
    if isinstance(value, dict):
        return {str(key): scalar(item) for key, item in value.items()}
    return value


def key(record):
    return (record.contig, int(record.pos), record.ref, tuple(record.alts or ()))


def is_snv(record):
    return len(record.ref) == 1 and len(record.alts or ()) == 1 and len(record.alts[0]) == 1


def make_groups(records, distance):
    groups = []
    current = []
    for record in records:
        if current and record.pos - current[-1].pos > distance:
            groups.append(current)
            current = []
        current.append(record)
    if current:
        groups.append(current)
    return groups


def recalibrated_records(path):
    result = {}
    keys = Counter()
    pass_keys = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            record_key = key(record)
            sample_values = {
                sample: {format_key: scalar(call.get(format_key)) for format_key in record.format.keys()}
                for sample, call in record.samples.items()
            }
            result[record_key] = {
                "id": record.id or ".",
                "qual": None if record.qual is None else record.qual,
                "filter": ";".join(record.filter.keys()) or ".",
                "info": {name: scalar(value) for name, value in record.info.items()},
                "format_keys": list(record.format.keys()),
                "sample_values": sample_values,
            }
            keys[record_key] += 1
            if set(record.filter.keys()) == {"PASS"}:
                pass_keys[record_key] += 1
    return result, keys, pass_keys


def load_mlhp(paths):
    groups = {}
    params = []
    for path in paths:
        doc = json.loads(path.read_text(encoding="utf-8"))
        params.append(doc.get("params", {}))
        for group in doc.get("groups", []):
            positions = tuple(int(value) for value in group.get("positions", []))
            groups[(group["chrom"], positions)] = group
    return groups, params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--longphase-input-vcf", required=True, type=Path,
                        help="Normalized ClairS VCF supplied to LongPhase-S")
    parser.add_argument("--longphase-input-contract",
                        choices=("clairs_PASS_only", "clairs_raw_all"),
                        default="clairs_PASS_only")
    parser.add_argument("--tree-input-vcf", required=True, type=Path,
                        help="Explicit VCF supplied to sSNV reconstruction; role is fixed by --tree-contract")
    parser.add_argument("--tree-contract",
                        choices=("longphase_recalibrated_PASS", "clairs_PASS_input", "clairs_PASS_backbone"),
                        default="longphase_recalibrated_PASS")
    parser.add_argument("--caller-raw-vcf", type=Path,
                        help="Raw ClairS output; every record is retained in the ledger")
    parser.add_argument("--recalibrated-vcf", required=True, type=Path)
    parser.add_argument("--mlhp-glob", required=True)
    parser.add_argument("--output-tsv-gz", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    args = parser.parse_args()
    mlhp_paths = [Path(path) for path in sorted(glob.glob(args.mlhp_glob))]
    if not mlhp_paths:
        raise SystemExit(f"no mlhp inputs: {args.mlhp_glob}")
    mlhp_groups, params = load_mlhp(mlhp_paths)
    tier_values = {int(item.get("TIER_R", 50000)) for item in params}
    cap_values = {int(item.get("MAX_SNV", 8)) for item in params}
    if len(tier_values) != 1 or len(cap_values) != 1:
        raise SystemExit(f"inconsistent mlhp params: TIER_R={tier_values} MAX_SNV={cap_values}")
    tier_r, max_snv = tier_values.pop(), cap_values.pop()
    recal, recal_keys, recal_pass_keys = recalibrated_records(args.recalibrated_vcf)
    with pysam.VariantFile(str(args.longphase_input_vcf)) as vcf:
        longphase_records = [record.copy() for record in vcf]
    with pysam.VariantFile(str(args.tree_input_vcf)) as vcf:
        tree_records = [record.copy() for record in vcf]
    records = longphase_records
    if args.caller_raw_vcf:
        with pysam.VariantFile(str(args.caller_raw_vcf)) as vcf:
            records = [record.copy() for record in vcf]
    longphase_keys = Counter(key(record) for record in longphase_records)
    tree_keys = Counter(key(record) for record in tree_records)
    raw_keys = Counter(key(record) for record in records)
    raw_pass_keys = Counter(key(record) for record in records if set(record.filter.keys()) == {"PASS"})
    duplicate_counts = {
        "raw_clairs": sum(value - 1 for value in Counter(key(record) for record in records).values()),
        "longphase_input": sum(value - 1 for value in longphase_keys.values()),
        "longphase_recalibrated": sum(value - 1 for value in recal_keys.values()),
        "tree_input": sum(value - 1 for value in tree_keys.values()),
    }
    status = {}
    seen_mlhp_keys = set()
    branch_counts = Counter()
    by_chrom = defaultdict(list)
    for record in tree_records:
        if not is_snv(record):
            status[key(record)] = {"branch": "unsupported_non_biallelic", "selected": False}
            branch_counts["unsupported_non_biallelic"] += 1
        elif record.contig not in AUTOSOMES:
            status[key(record)] = {"branch": "out_of_scope_non_autosomal", "selected": False}
            branch_counts["out_of_scope_non_autosomal"] += 1
        else:
            by_chrom[record.contig].append(record)
    for chrom, chrom_records in by_chrom.items():
        for component_index, component in enumerate(make_groups(chrom_records, tier_r), start=1):
            component_positions = [record.pos for record in component]
            component_id = f"{chrom}:{component_positions[0]}-{component_positions[-1]}"
            if len(component) < 2:
                record = component[0]
                status[key(record)] = {"branch": "positional_singleton", "selected": False,
                                       "component_id": component_id, "component_size": 1}
                branch_counts["positional_singleton"] += 1
                continue
            selected = component
            if len(component) > max_snv:
                best = min(range(len(component) - max_snv + 1),
                           key=lambda index: component[index + max_snv - 1].pos - component[index].pos)
                selected = component[best:best + max_snv]
            selected_keys = {key(record) for record in selected}
            selected_positions = tuple(record.pos for record in selected)
            selected_id = f"{chrom}:{selected_positions[0]}-{selected_positions[-1]}"
            raw = mlhp_groups.get((chrom, selected_positions))
            for record in component:
                record_key = key(record)
                common = {"component_id": component_id, "component_size": len(component),
                          "selected_group_id": selected_id, "selected_group_size": len(selected)}
                if record_key not in selected_keys:
                    status[record_key] = {**common, "branch": "max_snv_excluded", "selected": False}
                    branch_counts["max_snv_excluded"] += 1
                elif raw is None:
                    status[record_key] = {**common, "branch": "read_unsupported", "selected": True}
                    branch_counts["read_unsupported"] += 1
                else:
                    seen_mlhp_keys.add((chrom, selected_positions))
                    coverage = (raw.get("col_coverage") or {}).get(str(record.pos), [0, 0])
                    family_coverage = {
                        family: (values or {}).get(str(record.pos), [0, 0])
                        for family, values in (raw.get("col_coverage_by_hp") or {}).items()
                    }
                    status[record_key] = {
                        **common, "branch": "retained", "selected": True,
                        "pooled_ref_reads": coverage[0], "pooled_alt_reads": coverage[1],
                        "family_coverage": family_coverage,
                        "raw_HP_counts": raw.get("raw_HP_counts", {}),
                        "raw_HP_with_PS_counts": raw.get("raw_HP_with_PS_counts", {}),
                        "n_unique_phase_sets": raw.get("n_unique_phase_sets", 0),
                        "phase_set_counts": raw.get("phase_set_counts", {}),
                        "phase_set_HP_counts": raw.get("phase_set_HP_counts", {}),
                    }
                    branch_counts["retained"] += 1
    args.output_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample", "chrom", "pos", "id", "ref", "alt", "qual", "caller_filter",
        "info_json", "format_keys", "sample_values_json",
        "longphase_id", "longphase_qual", "longphase_recalibrated_filter", "longphase_info_json",
        "longphase_filter_transition",
        "longphase_format_keys", "longphase_sample_values_json",
        "ssnv_branch", "selected_for_read_census", "component_id", "component_size",
        "selected_group_id", "selected_group_size", "pooled_ref_reads", "pooled_alt_reads",
        "family_coverage_json", "raw_HP_counts_json", "raw_HP_with_PS_counts_json", "n_unique_phase_sets",
        "phase_set_counts_json", "phase_set_HP_counts_json",
    ]
    recal_missing = 0
    with pysam.BGZFile(str(args.output_tsv_gz), "w") as compressed, \
            io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in records:
            record_key = key(record)
            is_longphase_input = record_key in longphase_keys
            is_tree_input = record_key in tree_keys
            if not is_longphase_input:
                branch = ("excluded_before_longphase_nonPASS"
                          if args.longphase_input_contract == "clairs_PASS_only"
                          else "missing_from_longphase_input")
                disposition = {"branch": branch, "selected": False}
                branch_counts[branch] += 1
            elif not is_tree_input:
                disposition = {"branch": "excluded_by_longphase_filter", "selected": False}
                branch_counts["excluded_by_longphase_filter"] += 1
            else:
                disposition = status[record_key]
            recal_record = recal.get(record_key) if is_longphase_input else None
            recal_missing += int(is_longphase_input and recal_record is None)
            sample_values = {
                sample: {format_key: scalar(call.get(format_key)) for format_key in record.format.keys()}
                for sample, call in record.samples.items()
            }
            writer.writerow({
                "sample": args.sample, "chrom": record.contig, "pos": record.pos,
                "id": record.id or ".", "ref": record.ref, "alt": ",".join(record.alts or ()),
                "qual": "." if record.qual is None else record.qual,
                "caller_filter": ";".join(record.filter.keys()) or ".",
                "info_json": json.dumps({name: scalar(value) for name, value in record.info.items()},
                                        ensure_ascii=False, separators=(",", ":")),
                "format_keys": ":".join(record.format.keys()),
                "sample_values_json": json.dumps(sample_values, ensure_ascii=False, separators=(",", ":")),
                "longphase_id": recal_record["id"] if recal_record else "",
                "longphase_qual": ("" if not recal_record or recal_record["qual"] is None
                                    else recal_record["qual"]),
                "longphase_recalibrated_filter": (
                    recal_record["filter"] if recal_record else ("MISSING" if is_longphase_input else "NOT_INPUT")
                ),
                "longphase_filter_transition": (
                    (";".join(record.filter.keys()) or ".") + "->" +
                    (recal_record["filter"] if recal_record else
                     ("MISSING" if is_longphase_input else "NOT_INPUT"))
                ),
                "longphase_info_json": json.dumps(
                    recal_record["info"] if recal_record else {}, ensure_ascii=False, separators=(",", ":")
                ),
                "longphase_format_keys": ":".join(recal_record["format_keys"]) if recal_record else "",
                "longphase_sample_values_json": json.dumps(
                    recal_record["sample_values"] if recal_record else {},
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
                "ssnv_branch": disposition["branch"],
                "selected_for_read_census": str(disposition.get("selected", False)).lower(),
                "component_id": disposition.get("component_id", ""),
                "component_size": disposition.get("component_size", ""),
                "selected_group_id": disposition.get("selected_group_id", ""),
                "selected_group_size": disposition.get("selected_group_size", ""),
                "pooled_ref_reads": disposition.get("pooled_ref_reads", ""),
                "pooled_alt_reads": disposition.get("pooled_alt_reads", ""),
                "family_coverage_json": json.dumps(disposition.get("family_coverage", {}), separators=(",", ":")),
                "raw_HP_counts_json": json.dumps(disposition.get("raw_HP_counts", {}), separators=(",", ":")),
                "raw_HP_with_PS_counts_json": json.dumps(disposition.get("raw_HP_with_PS_counts", {}), separators=(",", ":")),
                "n_unique_phase_sets": disposition.get("n_unique_phase_sets", ""),
                "phase_set_counts_json": json.dumps(disposition.get("phase_set_counts", {}), separators=(",", ":")),
                "phase_set_HP_counts_json": json.dumps(disposition.get("phase_set_HP_counts", {}), separators=(",", ":")),
            })
    pysam.tabix_index(
        str(args.output_tsv_gz),
        seq_col=1,
        start_col=2,
        end_col=2,
        line_skip=1,
        zerobased=False,
        force=True,
    )
    autosomal_snv = sum(1 for record in tree_records if is_snv(record) and record.contig in AUTOSOMES)
    autosomal_accounted = sum(branch_counts[name] for name in
                              ("positional_singleton", "max_snv_excluded", "read_unsupported", "retained"))
    checks = {
        "all_raw_records_written": sum(branch_counts.values()) == len(records),
        "longphase_input_equals_recalibrated_all": recal_missing == 0 and longphase_keys == recal_keys,
        "autosomal_snv_conservation": autosomal_accounted == autosomal_snv,
        "all_mlhp_retained_groups_joined": seen_mlhp_keys == set(mlhp_groups),
        "record_keys_are_unique_at_all_four_layers": all(value == 0 for value in duplicate_counts.values()),
    }
    if args.longphase_input_contract == "clairs_raw_all":
        checks["raw_all_equals_longphase_input"] = raw_keys == longphase_keys
    else:
        checks["raw_PASS_equals_longphase_input"] = raw_pass_keys == longphase_keys
    if args.tree_contract == "longphase_recalibrated_PASS":
        checks["recalibrated_PASS_equals_tree_input"] = recal_pass_keys == tree_keys
    elif args.tree_contract == "clairs_PASS_input":
        if args.longphase_input_contract == "clairs_raw_all":
            checks["caller_raw_PASS_equals_tree_input"] = raw_pass_keys == tree_keys
        else:
            checks["longphase_input_equals_tree_input"] = longphase_keys == tree_keys
    else:
        checks["caller_raw_PASS_equals_tree_input"] = raw_pass_keys == tree_keys
    transition_counts = Counter()
    for record in records:
        recal_record = recal.get(key(record))
        caller_filter = ";".join(record.filter.keys()) or "."
        output_filter = recal_record["filter"] if recal_record else "NOT_INPUT"
        transition_counts[f"{caller_filter}->{output_filter}"] += 1
    summary = {
        "schema_version": "2.0", "sample": args.sample,
        "longphase_input_contract": args.longphase_input_contract,
        "tree_contract": args.tree_contract,
        "caller_raw_vcf": str(args.caller_raw_vcf) if args.caller_raw_vcf else str(args.longphase_input_vcf),
        "longphase_input_vcf": str(args.longphase_input_vcf),
        "tree_input_vcf": str(args.tree_input_vcf),
        "recalibrated_vcf": str(args.recalibrated_vcf),
        "mlhp_inputs": [str(path) for path in mlhp_paths], "output_tsv_gz": str(args.output_tsv_gz),
        "output_index": str(args.output_tsv_gz) + ".tbi",
        "TIER_R": tier_r, "MAX_SNV": max_snv, "raw_clairs_records": len(records),
        "longphase_input_records": len(longphase_records),
        "longphase_recalibrated_records": sum(recal_keys.values()),
        "duplicate_record_key_excess": duplicate_counts,
        "tree_input_records": len(tree_records),
        "autosomal_biallelic_snvs": autosomal_snv, "branch_counts": dict(branch_counts),
        "filter_transition_counts": dict(sorted(transition_counts.items())),
        "recalibrated_key_missing": recal_missing, "checks": checks, "pass": all(checks.values()),
    }
    args.output_summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"SITE LEDGER [{args.sample}]: rows={len(records)} branches={dict(branch_counts)} -> {args.output_tsv_gz}")
    raise SystemExit(0 if summary["pass"] else 1)


if __name__ == "__main__":
    main()
