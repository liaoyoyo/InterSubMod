#!/usr/bin/env python3
"""Build fail-closed PyClone-VI inputs for the cross-source validation.

The script deliberately uses exact CHROM:POS:REF:ALT joins, caller DP/AF counts,
and measured SAVANA allele-specific CN. It never fills missing tumour CN with 2.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, NamedTuple, Optional, Sequence, Set, Tuple


AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
Key = Tuple[str, int, str, str]


class VcfDatum(NamedTuple):
    dp: int
    af: float
    ref_counts: int
    alt_counts: int


class Segment(NamedTuple):
    chrom: str
    start: int
    end: int
    segment_id: str
    total_raw: float
    minor_raw: float
    total_i: int
    minor_i: int
    major_i: int
    near_integer: bool


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_gzip_text(path: Path, rows: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                for row in rows:
                    text.write(row)


def run_bcftools_query(vcf: Path) -> Iterator[Tuple[Key, str, str]]:
    """Yield exact biallelic autosomal SNVs with raw DP and AF strings."""
    command = [
        "bcftools",
        "query",
        "-f",
        "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DP\\t%AF]\\n",
        str(vcf),
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        chrom, pos_text, ref, alt, dp_text, af_text = fields[:6]
        if chrom not in AUTOSOMES or len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        yield (chrom, int(pos_text), ref.upper(), alt.upper()), dp_text, af_text
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"bcftools query failed ({return_code}) for {vcf}: {stderr.strip()}")


def run_bcftools_keys(vcf: Path) -> Iterator[Key]:
    """Yield exact biallelic autosomal SNV keys, including sample-less truth VCFs."""
    command = ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(vcf)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        chrom, pos_text, ref, alt = fields[:4]
        if chrom not in AUTOSOMES or len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        yield chrom, int(pos_text), ref.upper(), alt.upper()
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"bcftools query failed ({return_code}) for {vcf}: {stderr.strip()}")


def parse_single_number(text: str) -> Optional[float]:
    if text in {"", "."}:
        return None
    value_text = text.split(",", 1)[0]
    try:
        value = float(value_text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def read_vcf(vcf: Path, require_counts: bool) -> Tuple[Dict[Key, VcfDatum], Set[Key], Mapping[str, int]]:
    data: Dict[Key, VcfDatum] = {}
    keys: Set[Key] = set()
    qa: MutableMapping[str, int] = {
        "records_autosomal_biallelic_snv": 0,
        "records_duplicate_exact_key": 0,
        "records_missing_or_invalid_dp_af": 0,
        "records_count_conservation_fail": 0,
    }
    for key, dp_text, af_text in run_bcftools_query(vcf):
        qa["records_autosomal_biallelic_snv"] += 1
        if key in keys:
            qa["records_duplicate_exact_key"] += 1
            continue
        keys.add(key)
        dp_value = parse_single_number(dp_text)
        af_value = parse_single_number(af_text)
        if dp_value is None or af_value is None or dp_value <= 0 or not 0 <= af_value <= 1:
            qa["records_missing_or_invalid_dp_af"] += 1
            continue
        dp = int(math.floor(dp_value + 0.5))
        alt_counts = int(math.floor(dp * af_value + 0.5))
        alt_counts = min(max(alt_counts, 0), dp)
        ref_counts = dp - alt_counts
        if ref_counts + alt_counts != dp:
            qa["records_count_conservation_fail"] += 1
            continue
        data[key] = VcfDatum(dp, af_value, ref_counts, alt_counts)
    if require_counts and not data:
        raise ValueError(f"No usable DP/AF records in {vcf}")
    return data, keys, qa


def read_truth_keys(vcf: Path) -> Tuple[Set[Key], Mapping[str, int]]:
    keys: Set[Key] = set()
    qa: MutableMapping[str, int] = {
        "records_autosomal_biallelic_snv": 0,
        "records_duplicate_exact_key": 0,
    }
    for key in run_bcftools_keys(vcf):
        qa["records_autosomal_biallelic_snv"] += 1
        if key in keys:
            qa["records_duplicate_exact_key"] += 1
            continue
        keys.add(key)
    return keys, qa


def round_half_up(value: float) -> int:
    return int(math.floor(value + 0.5))


class SegmentIndex:
    def __init__(self, path: Path):
        self.path = path
        self.by_chrom: Dict[str, List[Segment]] = {}
        self.starts: Dict[str, List[int]] = {}
        self.qa: MutableMapping[str, int] = {
            "segments_total_rows": 0,
            "segments_autosomal": 0,
            "segments_missing_total_cn": 0,
            "segments_missing_minor_cn": 0,
            "segments_invalid_coordinates": 0,
            "segments_invalid_discrete_cn": 0,
            "segments_overlap": 0,
            "segments_near_integer": 0,
        }
        self._load()

    def _load(self) -> None:
        with self.path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"chromosome", "start", "end", "segment_id", "copyNumber", "minorAlleleCopyNumber"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(f"CN file {self.path} misses columns: {sorted(missing)}")
            for row in reader:
                self.qa["segments_total_rows"] += 1
                chrom = row["chromosome"]
                if chrom not in AUTOSOMES:
                    continue
                self.qa["segments_autosomal"] += 1
                try:
                    start, end = int(row["start"]), int(row["end"])
                except ValueError:
                    self.qa["segments_invalid_coordinates"] += 1
                    continue
                if start < 1 or end < start:
                    self.qa["segments_invalid_coordinates"] += 1
                    continue
                try:
                    total_raw = float(row["copyNumber"])
                except ValueError:
                    self.qa["segments_missing_total_cn"] += 1
                    continue
                try:
                    minor_raw = float(row["minorAlleleCopyNumber"])
                except ValueError:
                    self.qa["segments_missing_minor_cn"] += 1
                    continue
                if not math.isfinite(total_raw) or not math.isfinite(minor_raw):
                    self.qa["segments_missing_minor_cn"] += 1
                    continue
                total_i = round_half_up(total_raw)
                minor_nearest = round_half_up(minor_raw)
                minor_i = min(minor_nearest, total_i // 2)
                major_i = total_i - minor_i
                if total_i < 1 or minor_i < 0 or major_i < 1 or major_i < minor_i:
                    self.qa["segments_invalid_discrete_cn"] += 1
                    continue
                near_integer = abs(total_raw - total_i) <= 0.25 and abs(minor_raw - minor_nearest) <= 0.25
                if near_integer:
                    self.qa["segments_near_integer"] += 1
                segment = Segment(
                    chrom,
                    start,
                    end,
                    row["segment_id"],
                    total_raw,
                    minor_raw,
                    total_i,
                    minor_i,
                    major_i,
                    near_integer,
                )
                self.by_chrom.setdefault(chrom, []).append(segment)
        for chrom, segments in self.by_chrom.items():
            segments.sort(key=lambda value: (value.start, value.end))
            previous_end = 0
            for segment in segments:
                if segment.start <= previous_end:
                    self.qa["segments_overlap"] += 1
                previous_end = max(previous_end, segment.end)
            self.starts[chrom] = [segment.start for segment in segments]
        if self.qa["segments_overlap"]:
            raise ValueError(f"Overlapping SAVANA segments are forbidden: {self.path}")

    def find(self, key: Key) -> Optional[Segment]:
        chrom, pos, _, _ = key
        starts = self.starts.get(chrom)
        if not starts:
            return None
        index = bisect.bisect_right(starts, pos) - 1
        if index < 0:
            return None
        segment = self.by_chrom[chrom][index]
        return segment if segment.start <= pos <= segment.end else None


def mutation_id(key: Key) -> str:
    chrom, pos, ref, alt = key
    return f"{chrom}:{pos}:{ref}>{alt}"


def write_input_bundle(
    output_dir: Path,
    bundle_id: str,
    sample_data: Mapping[str, Mapping[Key, VcfDatum]],
    universe: Set[Key],
    cn_index: SegmentIndex,
    purity_by_sample: Mapping[str, float],
    selection_contract: str,
    near_integer_only: bool,
    source_paths: Mapping[str, str],
    warnings: Sequence[str],
) -> Mapping[str, object]:
    samples = list(sample_data)
    counters: MutableMapping[str, int] = {
        "universe_exact_keys": len(universe),
        "excluded_missing_counts_any_sample": 0,
        "excluded_no_measured_cn_segment": 0,
        "excluded_not_near_integer_cn": 0,
        "included_mutations": 0,
        "included_rows": 0,
        "count_conservation_fail": 0,
        "duplicate_mutation_sample_rows": 0,
    }
    selected: List[Tuple[Key, Segment]] = []
    for key in sorted(universe, key=lambda item: (int(item[0][3:]), item[1], item[2], item[3])):
        if any(key not in sample_data[sample] for sample in samples):
            counters["excluded_missing_counts_any_sample"] += 1
            continue
        segment = cn_index.find(key)
        if segment is None:
            counters["excluded_no_measured_cn_segment"] += 1
            continue
        if near_integer_only and not segment.near_integer:
            counters["excluded_not_near_integer_cn"] += 1
            continue
        selected.append((key, segment))

    input_path = output_dir / f"{bundle_id}.pyclone_input.tsv.gz"
    metadata_path = output_dir / f"{bundle_id}.site_metadata.tsv.gz"

    def input_rows() -> Iterator[str]:
        yield "mutation_id\tsample_id\tref_counts\talt_counts\tnormal_cn\tmajor_cn\tminor_cn\ttumour_content\n"
        seen: Set[Tuple[str, str]] = set()
        for key, segment in selected:
            mid = mutation_id(key)
            for sample in samples:
                pair = (mid, sample)
                if pair in seen:
                    counters["duplicate_mutation_sample_rows"] += 1
                    continue
                seen.add(pair)
                datum = sample_data[sample][key]
                if datum.ref_counts + datum.alt_counts != datum.dp:
                    counters["count_conservation_fail"] += 1
                    continue
                yield (
                    f"{mid}\t{sample}\t{datum.ref_counts}\t{datum.alt_counts}\t2\t"
                    f"{segment.major_i}\t{segment.minor_i}\t{purity_by_sample[sample]:.6f}\n"
                )
                counters["included_rows"] += 1

    def metadata_rows() -> Iterator[str]:
        count_columns = "".join(
            f"\t{sample}_DP\t{sample}_AF\t{sample}_ref_counts\t{sample}_alt_counts" for sample in samples
        )
        yield (
            "mutation_id\tchrom\tpos\tref\talt\tsegment_id\tsegment_start\tsegment_end\t"
            "total_cn_raw\tminor_cn_raw\ttotal_cn_discrete\tmajor_cn\tminor_cn\tnear_integer"
            f"{count_columns}\n"
        )
        for key, segment in selected:
            chrom, pos, ref, alt = key
            values = [
                mutation_id(key), chrom, str(pos), ref, alt, segment.segment_id,
                str(segment.start), str(segment.end), f"{segment.total_raw:.8g}",
                f"{segment.minor_raw:.8g}", str(segment.total_i), str(segment.major_i),
                str(segment.minor_i), "1" if segment.near_integer else "0",
            ]
            for sample in samples:
                datum = sample_data[sample][key]
                values.extend([
                    str(datum.dp), f"{datum.af:.10g}", str(datum.ref_counts), str(datum.alt_counts)
                ])
            yield "\t".join(values) + "\n"

    stable_gzip_text(input_path, input_rows())
    stable_gzip_text(metadata_path, metadata_rows())
    counters["included_mutations"] = len(selected)
    expected_rows = counters["included_mutations"] * len(samples)
    status = "PASS"
    failures: List[str] = []
    if counters["included_rows"] != expected_rows:
        status = "FAIL"
        failures.append(f"included_rows={counters['included_rows']} expected={expected_rows}")
    if counters["count_conservation_fail"] or counters["duplicate_mutation_sample_rows"]:
        status = "FAIL"
        failures.append("count conservation or mutation/sample uniqueness failed")
    if not selected:
        status = "FAIL"
        failures.append("no mutation survived fail-closed gates")
    qa: Dict[str, object] = {
        "schema_name": "intersubmod.pyclone_input_qa",
        "schema_version": "1.0.0",
        "bundle_id": bundle_id,
        "status": status,
        "failures": failures,
        "samples": samples,
        "selection_contract": selection_contract,
        "near_integer_only": near_integer_only,
        "copy_number_policy": (
            "SAVANA continuous allele-specific CN discretized with half-up rule; uncovered/missing/invalid CN excluded"
        ),
        "count_policy": "alt=floor(DP*AF+0.5); ref=DP-alt; exact conservation required",
        "counters": dict(counters),
        "source_paths": dict(source_paths),
        "warnings": list(warnings),
        "artifacts": {
            "pyclone_input": str(input_path.resolve()),
            "pyclone_input_sha256": sha256_file(input_path),
            "site_metadata": str(metadata_path.resolve()),
            "site_metadata_sha256": sha256_file(metadata_path),
        },
    }
    qa_path = output_dir / f"{bundle_id}.qa.json"
    qa_path.write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n")
    qa["artifacts"]["qa_json"] = str(qa_path.resolve())
    return qa


def file_receipt(path: Path) -> Mapping[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def parse_args() -> argparse.Namespace:
    topic = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=topic / "config" / "pyclone_validation_config.json")
    parser.add_argument("--output-dir", type=Path, default=topic / "data" / "pyclone_inputs")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = json.loads(args.config.read_text())
    manifest_path = Path(config["source_manifest"])
    manifest = json.loads(manifest_path.read_text())
    sample_manifest = {entry["sample"]: entry for entry in manifest["samples"]}
    args.output_dir.mkdir(parents=True, exist_ok=True)

    source_receipts: Dict[str, Mapping[str, object]] = {
        "config": file_receipt(args.config),
        "source_manifest": file_receipt(manifest_path),
    }
    vcf_cache: Dict[Tuple[str, bool], Tuple[Dict[Key, VcfDatum], Set[Key], Mapping[str, int]]] = {}
    truth_cache: Dict[str, Tuple[Set[Key], Mapping[str, int]]] = {}
    cn_cache: Dict[str, SegmentIndex] = {}

    def load_vcf(path_text: str, require_counts: bool = True):
        key = (path_text, require_counts)
        if key not in vcf_cache:
            path = Path(path_text)
            vcf_cache[key] = read_vcf(path, require_counts=require_counts)
            source_receipts.setdefault(path_text, file_receipt(path))
        return vcf_cache[key]

    def load_truth(path_text: str):
        if path_text not in truth_cache:
            path = Path(path_text)
            truth_cache[path_text] = read_truth_keys(path)
            source_receipts.setdefault(path_text, file_receipt(path))
        return truth_cache[path_text]

    def load_cn(path_text: str):
        if path_text not in cn_cache:
            path = Path(path_text)
            cn_cache[path_text] = SegmentIndex(path)
            source_receipts.setdefault(path_text, file_receipt(path))
        return cn_cache[path_text]

    bundles: List[Mapping[str, object]] = []
    source_qa: Dict[str, object] = {}

    pair = config["hcc1395_pair"]
    first, second = pair["samples"]
    first_tree_path = sample_manifest[first]["somatic"]["tree_vcf"]["path"]
    second_tree_path = sample_manifest[second]["somatic"]["tree_vcf"]["path"]
    first_all_path = sample_manifest[first]["somatic"]["longphase_recalibrated_all_vcf"]["path"]
    second_all_path = sample_manifest[second]["somatic"]["longphase_recalibrated_all_vcf"]["path"]
    first_tree_data, first_tree_keys, first_tree_qa = load_vcf(first_tree_path)
    second_tree_data, second_tree_keys, second_tree_qa = load_vcf(second_tree_path)
    first_all_data, _, first_all_qa = load_vcf(first_all_path)
    second_all_data, _, second_all_qa = load_vcf(second_all_path)
    truth_keys, pair_truth_qa = load_truth(pair["truth_vcf"])
    pair_cn = load_cn(pair["copy_number_tsv"])
    primary_universe = first_tree_keys & second_tree_keys & truth_keys
    union_universe = (first_tree_keys | second_tree_keys) & truth_keys
    pair_purity = {first: float(pair["purity"]), second: float(pair["purity"])}
    pair_warnings = [pair["dorado_copy_number_policy"]]
    source_qa["hcc1395_pair"] = {
        "first_tree": first_tree_qa,
        "second_tree": second_tree_qa,
        "first_all": first_all_qa,
        "second_all": second_all_qa,
        "truth": pair_truth_qa,
        "primary_shared_tree_truth_exact_keys": len(primary_universe),
        "sensitivity_pass_union_truth_exact_keys_before_both_count_gate": len(union_universe),
        "cn": dict(pair_cn.qa),
    }
    pair_contract = "exact (tree PASS HCC1395 intersect tree PASS HCC1395_DORADO intersect SEQC2 high-confidence sSNV)"
    union_contract = "exact ((tree PASS HCC1395 union tree PASS HCC1395_DORADO) intersect SEQC2 high-confidence sSNV), counts required in both recalibrated-all VCFs; not a BAM recount"
    pair_sources = {
        "HCC1395_tree_vcf": first_tree_path,
        "HCC1395_DORADO_tree_vcf": second_tree_path,
        "HCC1395_recalibrated_all_vcf": first_all_path,
        "HCC1395_DORADO_recalibrated_all_vcf": second_all_path,
        "truth_vcf": pair["truth_vcf"],
        "copy_number_tsv": pair["copy_number_tsv"],
    }
    for near in (False, True):
        suffix = "near_integer" if near else "main"
        bundles.append(write_input_bundle(
            args.output_dir,
            f"hcc1395_pair_primary_joint_{suffix}",
            {first: first_tree_data, second: second_tree_data},
            primary_universe,
            pair_cn,
            pair_purity,
            pair_contract,
            near,
            pair_sources,
            pair_warnings,
        ))
        for sample, data in ((first, first_tree_data), (second, second_tree_data)):
            bundles.append(write_input_bundle(
                args.output_dir,
                f"hcc1395_pair_primary_separate_{sample}_{suffix}",
                {sample: data},
                primary_universe,
                pair_cn,
                {sample: pair_purity[sample]},
                pair_contract + "; separate fit on the identical primary mutation universe",
                near,
                pair_sources,
                pair_warnings,
            ))
        bundles.append(write_input_bundle(
            args.output_dir,
            f"hcc1395_pair_pass_union_joint_{suffix}",
            {first: first_all_data, second: second_all_data},
            union_universe,
            pair_cn,
            pair_purity,
            union_contract,
            near,
            pair_sources,
            pair_warnings + ["Sensitivity counts come from VCF DP/AF; this is explicitly not a complete BAM recount."],
        ))

    for sample, sample_config in config["individual_samples"].items():
        tree_path = sample_manifest[sample]["somatic"]["tree_vcf"]["path"]
        tree_data, tree_keys, tree_qa = load_vcf(tree_path)
        sample_truth_keys, sample_truth_qa = load_truth(sample_config["truth_vcf"])
        sample_cn = load_cn(sample_config["copy_number_tsv"])
        universe = tree_keys & sample_truth_keys
        source_qa[sample] = {
            "tree": tree_qa,
            "truth": sample_truth_qa,
            "tree_truth_exact_keys": len(universe),
            "cn": dict(sample_cn.qa),
        }
        sources = {
            "tree_vcf": tree_path,
            "truth_vcf": sample_config["truth_vcf"],
            "copy_number_tsv": sample_config["copy_number_tsv"],
        }
        contract = f"exact ({sample} tree PASS intersect {sample_config['truth_label']})"
        for near in (False, True):
            suffix = "near_integer" if near else "main"
            bundles.append(write_input_bundle(
                args.output_dir,
                f"{sample}_individual_{suffix}",
                {sample: tree_data},
                universe,
                sample_cn,
                {sample: float(sample_config["purity"])},
                contract,
                near,
                sources,
                [sample_config["truth_label"]],
            ))

    blocked_path = args.output_dir / "blocked_samples.tsv"
    blocked_rows = ["sample\tstatus\treason\n"] + [
        f"{sample}\tBLOCKED\t{reason}\n" for sample, reason in config["blocked_samples"].items()
    ]
    stable_gzip_text(blocked_path.with_suffix(".tsv.gz"), blocked_rows)

    summary_rows = [
        "bundle_id\tstatus\tsamples\tnear_integer_only\tuniverse_exact_keys\tincluded_mutations\t"
        "included_rows\texcluded_missing_counts_any_sample\texcluded_no_measured_cn_segment\t"
        "excluded_not_near_integer_cn\tinput_path\n"
    ]
    for qa in bundles:
        counters = qa["counters"]
        summary_rows.append(
            "\t".join([
                str(qa["bundle_id"]), str(qa["status"]), ",".join(qa["samples"]),
                str(int(bool(qa["near_integer_only"]))), str(counters["universe_exact_keys"]),
                str(counters["included_mutations"]), str(counters["included_rows"]),
                str(counters["excluded_missing_counts_any_sample"]),
                str(counters["excluded_no_measured_cn_segment"]),
                str(counters["excluded_not_near_integer_cn"]),
                str(qa["artifacts"]["pyclone_input"]),
            ]) + "\n"
        )
    summary_path = args.output_dir / "input_qa_summary.tsv"
    summary_path.write_text("".join(summary_rows))
    provenance = {
        "schema_name": "intersubmod.pyclone_input_build_provenance",
        "schema_version": "1.0.0",
        "builder": str(Path(__file__).resolve()),
        "builder_sha256": sha256_file(Path(__file__).resolve()),
        "command": [str(Path(__file__).resolve()), "--config", str(args.config.resolve()), "--output-dir", str(args.output_dir.resolve())],
        "python_version": sys.version,
        "bcftools_version": subprocess.run(["bcftools", "--version"], capture_output=True, text=True, check=True).stdout.splitlines()[0],
        "config": config,
        "source_receipts": source_receipts,
        "source_qa": source_qa,
        "bundles": bundles,
        "blocked_samples": config["blocked_samples"],
        "summary_tsv": str(summary_path.resolve()),
    }
    provenance_path = args.output_dir / "provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    failed = [qa["bundle_id"] for qa in bundles if qa["status"] != "PASS"]
    print(f"Built {len(bundles)} PyClone-VI bundles in {args.output_dir}")
    print(f"Primary shared HCC1395 exact universe: {len(primary_universe)}")
    print(f"PASS-union HCC1395 sensitivity universe before count/CN gates: {len(union_universe)}")
    print(f"QA summary: {summary_path}")
    print(f"Provenance: {provenance_path}")
    if failed:
        print(f"FAIL bundles: {', '.join(failed)}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
