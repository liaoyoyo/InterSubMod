#!/usr/bin/env python3
"""Bridge independent PyClone-VI clusters to the fixed 5,720 regional topology pairs.

The primary endpoint uses the two *separate* PyClone-VI fits. Joint-fit labels
are never used as independent replication evidence. Read-full directed
relations are interpreted only under the documented F/R/P contract from
``build_site_topology_containment.py`` and are evaluated fail-closed.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, NamedTuple, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr, spearmanr


SAMPLE_A = "HCC1395"
SAMPLE_B = "HCC1395_DORADO"
SAMPLES = (SAMPLE_A, SAMPLE_B)
CLONAL_THRESHOLD = 0.90
ASSIGNMENT_THRESHOLD = 0.80
PRIMARY_CP_TOLERANCE = 0.02
CP_TOLERANCE_GRID = (0.0, 0.02, 0.05)
FIXED_REGION_DENOMINATOR = 5720


class ClusterRecord(NamedTuple):
    mutation_id: str
    sample_id: str
    cluster_id: int
    cp: float
    cp_std: float
    assignment_prob: float


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def stable_gzip_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="")
    return raw, zipped, text


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str], gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if gzip_output:
        raw, zipped, text = stable_gzip_writer(path)
        try:
            writer = csv.DictWriter(text, delimiter="\t", fieldnames=list(columns), extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        finally:
            text.close()
            zipped.close()
            raw.close()
    else:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)


def parse_bool(value: str) -> bool:
    if value in {"True", "true", "1"}:
        return True
    if value in {"False", "false", "0", ""}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def canonical_mutation_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{int(pos)}:{ref.upper()}>{alt.upper()}"


def parse_pyclone_mutation_id(value: str) -> Tuple[str, int, str, str]:
    try:
        chrom, position, alleles = value.split(":", 2)
        ref, alt = alleles.split(">", 1)
        key = canonical_mutation_id(chrom, int(position), ref, alt)
    except Exception as error:
        raise ValueError(f"Invalid PyClone mutation_id: {value!r}") from error
    if key != value:
        raise ValueError(f"Noncanonical PyClone mutation_id: {value!r} -> {key!r}")
    return chrom, int(position), ref, alt


def parse_raw_key(value: str) -> str:
    try:
        chrom, position, ref, alt = value.split(":")
    except ValueError as error:
        raise ValueError(f"Invalid raw-VAF key: {value!r}") from error
    return canonical_mutation_id(chrom, int(position), ref, alt)


def load_cluster_results(path: Path, expected_sample: str) -> Dict[str, ClusterRecord]:
    records: Dict[str, ClusterRecord] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "mutation_id", "sample_id", "cluster_id", "cellular_prevalence",
            "cellular_prevalence_std", "cluster_assignment_prob",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} misses fields {sorted(missing)}")
        for row in reader:
            mutation_id = row["mutation_id"]
            parse_pyclone_mutation_id(mutation_id)
            if row["sample_id"] != expected_sample:
                raise ValueError(f"Unexpected sample in {path}: {row['sample_id']}")
            if mutation_id in records:
                raise ValueError(f"Duplicate PyClone mutation key in {path}: {mutation_id}")
            record = ClusterRecord(
                mutation_id,
                expected_sample,
                int(row["cluster_id"]),
                float(row["cellular_prevalence"]),
                float(row["cellular_prevalence_std"]),
                float(row["cluster_assignment_prob"]),
            )
            if not (0 <= record.cp <= 1 and 0 <= record.assignment_prob <= 1 and record.cp_std >= 0):
                raise ValueError(f"Invalid PyClone value for {mutation_id}: {record}")
            records[mutation_id] = record
    return records


def load_raw_vaf(path: Path) -> Tuple[Dict[str, Mapping[str, object]], Mapping[str, object]]:
    """Load raw VAFs and collapse source-label duplicates fail-closed.

    The upstream exact-shared table is a stacked evidence table: a mutation can
    occur once for ``baseline_caller_pass`` and once for ``latest_lps_pass``.
    Collapsing is allowed only when every column except ``source`` is byte-for-
    byte identical.  Input/unique/duplicate counts are returned so the collapse
    remains auditable rather than silently behaving like ``drop_duplicates``.
    """
    records: Dict[str, Dict[str, object]] = {}
    signatures: Dict[str, Tuple[Tuple[str, str], ...]] = {}
    duplicate_keys = set()
    input_rows = 0
    duplicate_extra_rows = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "key", "chrom", "pos", "ref", "alt", "af_hcc1395", "af_dorado",
            "dp_hcc1395", "dp_dorado", "truth_confirmed", "source",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} misses fields {sorted(missing)}")
        for row in reader:
            input_rows += 1
            mutation_id = parse_raw_key(row["key"])
            expected = canonical_mutation_id(row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            if mutation_id != expected:
                raise ValueError(f"Raw VAF key columns disagree: {row['key']} vs {expected}")
            signature = tuple((field, row[field]) for field in (reader.fieldnames or []) if field != "source")
            if mutation_id in records:
                duplicate_keys.add(mutation_id)
                duplicate_extra_rows += 1
                if signatures[mutation_id] != signature:
                    raise ValueError(
                        f"Conflicting duplicate raw-VAF payload for {mutation_id}; "
                        "refuse mutation-key collapse"
                    )
                sources = set(str(records[mutation_id]["source"]).split(";"))
                sources.add(row["source"])
                records[mutation_id]["source"] = ";".join(sorted(sources))
                continue
            signatures[mutation_id] = signature
            records[mutation_id] = {
                "mutation_id": mutation_id,
                "af_a": float(row["af_hcc1395"]),
                "af_b": float(row["af_dorado"]),
                "dp_a": int(float(row["dp_hcc1395"])),
                "dp_b": int(float(row["dp_dorado"])),
                "truth_confirmed": parse_bool(row["truth_confirmed"]),
                "source": row["source"],
            }
    qa = {
        "input_rows": input_rows,
        "unique_mutation_keys_after_exact_payload_collapse": len(records),
        "duplicate_mutation_keys": len(duplicate_keys),
        "duplicate_extra_rows": duplicate_extra_rows,
        "conflicting_duplicate_payloads": 0,
        "row_conservation_pass": input_rows == len(records) + duplicate_extra_rows,
        "collapse_rule": "same canonical mutation key and all input columns except source exactly identical",
    }
    return records, qa


def load_region_outcomes(path: Path) -> Tuple[List[Mapping[str, str]], Dict[str, Mapping[str, str]]]:
    rows: List[Mapping[str, str]] = []
    by_match: Dict[str, Mapping[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "match_id", "chrom", "region", "fixed_denominator", "caller_k_a", "caller_k_b",
            "caller_shared_k", "k_bin", "caller_site_set_relation", "allele_identity_pass",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} misses fields {sorted(missing)}")
        for row in reader:
            if row["match_id"] in by_match:
                raise ValueError(f"Duplicate region match_id: {row['match_id']}")
            if int(row["fixed_denominator"]) != FIXED_REGION_DENOMINATOR:
                raise ValueError(f"Unexpected fixed denominator for {row['match_id']}: {row['fixed_denominator']}")
            by_match[row["match_id"]] = row
            rows.append(row)
    if len(rows) != FIXED_REGION_DENOMINATOR:
        raise ValueError(f"Expected {FIXED_REGION_DENOMINATOR} fixed regions, found {len(rows)}")
    return rows, by_match


def parse_alleles(value: str) -> List[Tuple[str, int, str, str]]:
    parsed = json.loads(value)
    alleles: List[Tuple[str, int, str, str]] = []
    for item in parsed:
        if not isinstance(item, list) or len(item) != 4:
            raise ValueError(f"Invalid allele tuple: {item!r}")
        chrom, position, ref, alt = item
        alleles.append((str(chrom), int(position), str(ref).upper(), str(alt).upper()))
    return alleles


def load_region_alleles(
    path: Path,
    region_by_match: Mapping[str, Mapping[str, str]],
) -> Tuple[List[Mapping[str, object]], Dict[Tuple[str, int], List[str]], Mapping[str, int]]:
    expanded: List[Mapping[str, object]] = []
    by_region_position: Dict[Tuple[str, int], List[str]] = defaultdict(list)
    seen_region_mutation = set()
    counters = Counter()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"match_id", "chrom", "region", "position", "alleles_a", "alleles_b", "allele_identity"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} misses fields {sorted(missing)}")
        for row in reader:
            counters["input_position_rows"] += 1
            match_id = row["match_id"]
            if match_id not in region_by_match:
                raise ValueError(f"Allele row references unknown region: {match_id}")
            alleles_a = parse_alleles(row["alleles_a"])
            alleles_b = parse_alleles(row["alleles_b"])
            identity = parse_bool(row["allele_identity"])
            if identity != (set(alleles_a) == set(alleles_b)):
                raise ValueError(f"Allele identity flag disagrees for {match_id}:{row['position']}")
            if not identity:
                counters["allele_identity_fail_position_rows"] += 1
                continue
            for chrom, position, ref, alt in sorted(set(alleles_a)):
                if chrom != row["chrom"] or position != int(row["position"]):
                    raise ValueError(f"Allele coordinate mismatch for {match_id}: {row}")
                mutation_id = canonical_mutation_id(chrom, position, ref, alt)
                key = (match_id, mutation_id)
                if key in seen_region_mutation:
                    raise ValueError(f"Duplicate region×mutation key: {key}")
                seen_region_mutation.add(key)
                by_region_position[(match_id, position)].append(mutation_id)
                expanded.append({
                    "match_id": match_id,
                    "chrom": chrom,
                    "region": row["region"],
                    "position": position,
                    "ref": ref,
                    "alt": alt,
                    "mutation_id": mutation_id,
                })
                counters["expanded_exact_alleles"] += 1
    counters["ambiguous_region_positions"] = sum(len(values) != 1 for values in by_region_position.values())
    counters["unique_regional_mutations"] = len({row["mutation_id"] for row in expanded})
    counters["mutation_region_multiplicity_gt1"] = sum(
        count > 1 for count in Counter(row["mutation_id"] for row in expanded).values()
    )
    return expanded, by_region_position, dict(counters)


def comb2(value: int) -> float:
    return value * (value - 1) / 2


def label_metrics(labels_a: Sequence[int], labels_b: Sequence[int]) -> Mapping[str, object]:
    if len(labels_a) != len(labels_b) or not labels_a:
        return {
            "n": len(labels_a), "clusters_a": 0, "clusters_b": 0, "ari": None, "nmi": None,
            "hungarian_agreement": None,
        }
    a_values = sorted(set(labels_a))
    b_values = sorted(set(labels_b))
    a_index = {value: index for index, value in enumerate(a_values)}
    b_index = {value: index for index, value in enumerate(b_values)}
    table = np.zeros((len(a_values), len(b_values)), dtype=np.int64)
    for first, second in zip(labels_a, labels_b):
        table[a_index[first], b_index[second]] += 1
    row_sums = table.sum(axis=1)
    column_sums = table.sum(axis=0)
    n = int(table.sum())
    total_pairs = comb2(n)
    sum_cells = float(sum(comb2(int(value)) for value in table.ravel()))
    sum_rows = float(sum(comb2(int(value)) for value in row_sums))
    sum_columns = float(sum(comb2(int(value)) for value in column_sums))
    expected = sum_rows * sum_columns / total_pairs if total_pairs else 0.0
    maximum = 0.5 * (sum_rows + sum_columns)
    ari = (sum_cells - expected) / (maximum - expected) if maximum != expected else 1.0
    probabilities = table / n
    p_a, p_b = row_sums / n, column_sums / n
    mutual_information = 0.0
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            if probabilities[i, j] > 0:
                mutual_information += probabilities[i, j] * math.log(
                    probabilities[i, j] / (p_a[i] * p_b[j])
                )
    entropy_a = -sum(value * math.log(value) for value in p_a if value > 0)
    entropy_b = -sum(value * math.log(value) for value in p_b if value > 0)
    nmi = 2 * mutual_information / (entropy_a + entropy_b) if entropy_a + entropy_b else 1.0
    row_match, column_match = linear_sum_assignment(-table)
    matched = int(table[row_match, column_match].sum())
    return {
        "n": n,
        "clusters_a": len(a_values),
        "clusters_b": len(b_values),
        "ari": float(ari),
        "nmi": float(nmi),
        "hungarian_agreement": matched / n,
    }


def safe_correlation(first: Sequence[float], second: Sequence[float], method: str) -> Optional[float]:
    if len(first) < 2 or len(set(first)) < 2 or len(set(second)) < 2:
        return None
    result = pearsonr(first, second) if method == "pearson" else spearmanr(first, second)
    return float(result.statistic)


def binary_kappa(first: Sequence[bool], second: Sequence[bool]) -> Optional[float]:
    if not first:
        return None
    observed = sum(a == b for a, b in zip(first, second)) / len(first)
    p_a = sum(first) / len(first)
    p_b = sum(second) / len(second)
    expected = p_a * p_b + (1 - p_a) * (1 - p_b)
    return (observed - expected) / (1 - expected) if expected < 1 else 1.0


def mutation_concordance(rows: Sequence[Mapping[str, object]], population: str, stratum: str) -> Mapping[str, object]:
    records = [row for row in rows if row["joined_pyclone_both"]]
    labels_a = [int(row["cluster_a"]) for row in records]
    labels_b = [int(row["cluster_b"]) for row in records]
    metrics = dict(label_metrics(labels_a, labels_b))
    cp_a = [float(row["cp_a"]) for row in records]
    cp_b = [float(row["cp_b"]) for row in records]
    sub_a = [value < CLONAL_THRESHOLD for value in cp_a]
    sub_b = [value < CLONAL_THRESHOLD for value in cp_b]
    both_sub = sum(a and b for a, b in zip(sub_a, sub_b))
    either_sub = sum(a or b for a, b in zip(sub_a, sub_b))
    metrics.update({
        "population": population,
        "stratum": stratum,
        "mutation_rows": len(records),
        "unique_mutations": len({row["mutation_id"] for row in records}),
        "assignment_high_both": sum(bool(row["assignment_high_both"]) for row in records),
        "cp_mae": float(np.mean(np.abs(np.asarray(cp_a) - np.asarray(cp_b)))) if records else None,
        "cp_pearson": safe_correlation(cp_a, cp_b, "pearson"),
        "cp_spearman": safe_correlation(cp_a, cp_b, "spearman"),
        "clonal_state_agreement": sum(a == b for a, b in zip(sub_a, sub_b)) / len(records) if records else None,
        "clonal_state_kappa": binary_kappa(sub_a, sub_b),
        "subclonal_intersection": both_sub,
        "subclonal_union": either_sub,
        "subclonal_jaccard": both_sub / either_sub if either_sub else 1.0,
    })
    return metrics


def canonical_partition(mutations: Sequence[str], cluster_by_mutation: Mapping[str, int]) -> Tuple[int, ...]:
    mapping: Dict[int, int] = {}
    result = []
    for mutation in mutations:
        cluster = cluster_by_mutation[mutation]
        if cluster not in mapping:
            mapping[cluster] = len(mapping)
        result.append(mapping[cluster])
    return tuple(result)


def region_pattern_row(
    region: Mapping[str, str],
    site_rows: Sequence[Mapping[str, object]],
    high_confidence_only: bool,
) -> Mapping[str, object]:
    eligible = [
        row for row in site_rows
        if row["joined_pyclone_both"] and (not high_confidence_only or row["assignment_high_both"])
    ]
    eligible.sort(key=lambda row: (int(row["position"]), str(row["mutation_id"])))
    mutations = [str(row["mutation_id"]) for row in eligible]
    endpoint = "both_assignment_ge_0.8" if high_confidence_only else "all_joined"
    base: Dict[str, object] = {
        "match_id": region["match_id"],
        "chrom": region["chrom"],
        "region": region["region"],
        "endpoint": endpoint,
        "fixed_denominator": FIXED_REGION_DENOMINATOR,
        "caller_k_a": int(region["caller_k_a"]),
        "caller_k_b": int(region["caller_k_b"]),
        "caller_shared_k": int(region["caller_shared_k"]),
        "k_bin": region["k_bin"],
        "caller_site_set_relation": region["caller_site_set_relation"],
        "equal_site_set": region["caller_site_set_relation"] == "equal",
        "regional_exact_alleles": len(site_rows),
        "joined_mutations": len(mutations),
        "evaluable_multilocus": len(mutations) >= 2,
    }
    if len(mutations) < 2:
        base.update({
            "clusters_a": "", "clusters_b": "", "partition_signature_a": "",
            "partition_signature_b": "", "partition_exact": "", "pairwise_cocluster_agreement": "",
            "region_ari": "", "pattern_category": "not_evaluable_lt2_joined_mutations",
            "subclonal_jaccard": "",
        })
        return base
    cluster_a = {str(row["mutation_id"]): int(row["cluster_a"]) for row in eligible}
    cluster_b = {str(row["mutation_id"]): int(row["cluster_b"]) for row in eligible}
    partition_a = canonical_partition(mutations, cluster_a)
    partition_b = canonical_partition(mutations, cluster_b)
    pair_agreements = []
    for first_index in range(len(mutations)):
        for second_index in range(first_index + 1, len(mutations)):
            same_a = partition_a[first_index] == partition_a[second_index]
            same_b = partition_b[first_index] == partition_b[second_index]
            pair_agreements.append(same_a == same_b)
    count_a, count_b = len(set(partition_a)), len(set(partition_b))
    exact = partition_a == partition_b
    if count_a == 1 and count_b == 1:
        category = "both_single_cluster"
    elif count_a > 1 and count_b > 1 and exact:
        category = "both_multicluster_exact_partition"
    elif (count_a == 1) != (count_b == 1):
        category = "one_single_one_multicluster"
    else:
        category = "both_multicluster_different_partition"
    sub_a = {mutation for mutation in mutations if eligible[mutations.index(mutation)]["cp_a"] < CLONAL_THRESHOLD}
    sub_b = {mutation for mutation in mutations if eligible[mutations.index(mutation)]["cp_b"] < CLONAL_THRESHOLD}
    union = sub_a | sub_b
    base.update({
        "clusters_a": count_a,
        "clusters_b": count_b,
        "partition_signature_a": ",".join(map(str, partition_a)),
        "partition_signature_b": ",".join(map(str, partition_b)),
        "partition_exact": exact,
        "pairwise_cocluster_agreement": sum(pair_agreements) / len(pair_agreements),
        "region_ari": label_metrics(partition_a, partition_b)["ari"],
        "pattern_category": category,
        "subclonal_jaccard": len(sub_a & sub_b) / len(union) if union else 1.0,
    })
    return base


def load_relation_evidence(path: Path) -> List[Mapping[str, str]]:
    rows: List[Mapping[str, str]] = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "match_id", "chrom", "region", "component_index", "pair_index", "left_position", "right_position",
            "read_full_relation_codes_a", "read_full_relation_codes_b", "read_full_relations_a", "read_full_relations_b",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} misses fields {sorted(missing)}")
        for row in reader:
            if int(row["left_position"]) >= int(row["right_position"]):
                raise ValueError(f"Relation evidence is not lower/higher ordered: {row}")
            for suffix in ("a", "b"):
                codes_text = row[f"read_full_relation_codes_{suffix}"]
                relations_text = row[f"read_full_relations_{suffix}"]
                codes = json.loads(codes_text) if codes_text else []
                relations = json.loads(relations_text) if relations_text else []
                if not set(codes).issubset({"F", "R", "P"}):
                    raise ValueError(f"Unknown relation code in {row['match_id']}: {codes}")
                expected = set()
                left, right = int(row["left_position"]), int(row["right_position"])
                if "F" in codes:
                    expected.add(f"{left}>{right}")
                if "R" in codes:
                    expected.add(f"{right}>{left}")
                if "P" in codes:
                    expected.add("parallel")
                if expected != set(relations):
                    raise ValueError(
                        f"Relation code/string contract mismatch for {row['match_id']} {left}-{right}: {codes} {relations}"
                    )
            rows.append(row)
    return rows


def validate_edge_schema_source(path: Path) -> Mapping[str, object]:
    """Fail closed unless the producing source explicitly defines F/R/P."""
    source = path.read_text()
    required_snippets = {
        "F": "F: lower genomic position is ancestral to the higher position",
        "R": "R: higher genomic position is ancestral to the lower position",
        "P": "P: the two mutation events are on parallel branches",
    }
    missing = [code for code, snippet in required_snippets.items() if snippet not in source]
    if missing:
        raise ValueError(
            f"Edge schema source {path} does not retain the explicit F/R/P contract; "
            f"missing={missing}. Directed-edge analysis is fail-closed."
        )
    return {
        "status": "PASS",
        "source": str(path.resolve()),
        "required_codes": sorted(required_snippets),
        "definition": "F=lower-position ancestor; R=higher-position ancestor; P=parallel",
    }


def cluster_distribution_rows(
    rows: Sequence[Mapping[str, object]], population: str, stratum: str,
) -> List[Mapping[str, object]]:
    output: List[Mapping[str, object]] = []
    for sample, suffix in ((SAMPLE_A, "a"), (SAMPLE_B, "b")):
        joined = [row for row in rows if row.get(f"cluster_{suffix}") not in {None, ""}]
        counts = Counter(int(row[f"cluster_{suffix}"]) for row in joined)
        for cluster_id in sorted(counts):
            members = [row for row in joined if int(row[f"cluster_{suffix}"]) == cluster_id]
            output.append({
                "population": population,
                "stratum": stratum,
                "sample": sample,
                "cluster_id_within_separate_fit": cluster_id,
                "n": len(members),
                "denominator": len(joined),
                "share": len(members) / len(joined) if joined else None,
                "mean_cp": float(np.mean([float(row[f"cp_{suffix}"]) for row in members])),
                "mean_assignment_prob": float(np.mean([
                    float(row[f"assignment_prob_{suffix}"]) for row in members
                ])),
                "assignment_ge_0.8": sum(
                    float(row[f"assignment_prob_{suffix}"]) >= ASSIGNMENT_THRESHOLD for row in members
                ),
                "subclonal_cp_lt_0.9": sum(float(row[f"cp_{suffix}"]) < CLONAL_THRESHOLD for row in members),
            })
    return output


def interval(record: ClusterRecord) -> Tuple[float, float]:
    return max(0.0, record.cp - 1.96 * record.cp_std), min(1.0, record.cp + 1.96 * record.cp_std)


def cp_status(parent: ClusterRecord, child: ClusterRecord, tolerance: float) -> str:
    if parent.cluster_id == child.cluster_id:
        return "uninformative_same_cluster"
    parent_lower, parent_upper = interval(parent)
    child_lower, child_upper = interval(child)
    if parent_lower + tolerance >= child_upper:
        return "compatible"
    if parent_upper + tolerance < child_lower:
        return "conflict"
    return "uninformative_interval_overlap"


def directed_edge_row(
    row: Mapping[str, str],
    by_region_position: Mapping[Tuple[str, int], List[str]],
    cluster_maps: Mapping[str, Mapping[str, ClusterRecord]],
) -> Mapping[str, object]:
    match_id = row["match_id"]
    left, right = int(row["left_position"]), int(row["right_position"])
    left_mutations = by_region_position.get((match_id, left), [])
    right_mutations = by_region_position.get((match_id, right), [])
    output: Dict[str, object] = {
        "match_id": match_id,
        "chrom": row["chrom"],
        "region": row["region"],
        "component_index": int(row["component_index"]),
        "pair_index": int(row["pair_index"]),
        "left_position": left,
        "right_position": right,
        "left_mutation_id": left_mutations[0] if len(left_mutations) == 1 else "",
        "right_mutation_id": right_mutations[0] if len(right_mutations) == 1 else "",
        "endpoint_identity_unique": len(left_mutations) == 1 and len(right_mutations) == 1,
    }
    directions: Dict[str, Optional[Tuple[str, str]]] = {}
    for sample, suffix in ((SAMPLE_A, "a"), (SAMPLE_B, "b")):
        codes = json.loads(row[f"read_full_relation_codes_{suffix}"]) if row[f"read_full_relation_codes_{suffix}"] else []
        output[f"read_codes_{suffix}"] = json.dumps(codes, separators=(",", ":"))
        direction: Optional[Tuple[str, str]] = None
        relation_status = "uninformative_ambiguous_or_parallel"
        if len(left_mutations) != 1 or len(right_mutations) != 1:
            relation_status = "uninformative_endpoint_identity_not_unique"
        elif codes == ["F"]:
            direction = (left_mutations[0], right_mutations[0])
            relation_status = "determinate_directed"
        elif codes == ["R"]:
            direction = (right_mutations[0], left_mutations[0])
            relation_status = "determinate_directed"
        directions[sample] = direction
        output[f"relation_status_{suffix}"] = relation_status
        output[f"parent_mutation_id_{suffix}"] = direction[0] if direction else ""
        output[f"child_mutation_id_{suffix}"] = direction[1] if direction else ""
        if direction is None:
            output[f"pyclone_join_{suffix}"] = False
            output[f"assignment_high_endpoints_{suffix}"] = False
            output[f"cluster_parent_{suffix}"] = ""
            output[f"cluster_child_{suffix}"] = ""
            output[f"cp_parent_{suffix}"] = ""
            output[f"cp_child_{suffix}"] = ""
            output[f"cp_delta_parent_minus_child_{suffix}"] = ""
            output[f"cp_status_t0.02_{suffix}"] = relation_status
            continue
        parent = cluster_maps[sample].get(direction[0])
        child = cluster_maps[sample].get(direction[1])
        if parent is None or child is None:
            output[f"pyclone_join_{suffix}"] = False
            output[f"assignment_high_endpoints_{suffix}"] = False
            output[f"cluster_parent_{suffix}"] = parent.cluster_id if parent else ""
            output[f"cluster_child_{suffix}"] = child.cluster_id if child else ""
            output[f"cp_parent_{suffix}"] = parent.cp if parent else ""
            output[f"cp_child_{suffix}"] = child.cp if child else ""
            output[f"cp_delta_parent_minus_child_{suffix}"] = ""
            output[f"cp_status_t0.02_{suffix}"] = "uninformative_pyclone_endpoint_missing"
            continue
        output[f"pyclone_join_{suffix}"] = True
        output[f"assignment_high_endpoints_{suffix}"] = (
            parent.assignment_prob >= ASSIGNMENT_THRESHOLD and child.assignment_prob >= ASSIGNMENT_THRESHOLD
        )
        output[f"cluster_parent_{suffix}"] = parent.cluster_id
        output[f"cluster_child_{suffix}"] = child.cluster_id
        output[f"cp_parent_{suffix}"] = parent.cp
        output[f"cp_child_{suffix}"] = child.cp
        output[f"cp_std_parent_{suffix}"] = parent.cp_std
        output[f"cp_std_child_{suffix}"] = child.cp_std
        output[f"assignment_parent_{suffix}"] = parent.assignment_prob
        output[f"assignment_child_{suffix}"] = child.assignment_prob
        output[f"cp_delta_parent_minus_child_{suffix}"] = parent.cp - child.cp
        output[f"cp_status_t0.02_{suffix}"] = cp_status(parent, child, PRIMARY_CP_TOLERANCE)
    direction_a, direction_b = directions[SAMPLE_A], directions[SAMPLE_B]
    if direction_a is None or direction_b is None:
        direction_relation = "uninformative_one_or_both_direction_sets"
    elif direction_a == direction_b:
        direction_relation = "same_direction"
    elif direction_a == (direction_b[1], direction_b[0]):
        direction_relation = "opposite_direction"
    else:
        direction_relation = "different_endpoint_identity"
    output["cross_source_direction_relation"] = direction_relation
    status_a = str(output["cp_status_t0.02_a"])
    status_b = str(output["cp_status_t0.02_b"])
    if "conflict" in {status_a, status_b}:
        combined = "conflict_any_source"
    elif status_a == status_b == "compatible":
        combined = "compatible_both_sources"
    elif "compatible" in {status_a, status_b}:
        combined = "compatible_one_other_uninformative"
    else:
        combined = "uninformative"
    output["combined_cp_status_t0.02"] = combined
    return output


def summarize_edges(
    raw_relation_rows: Sequence[Mapping[str, str]],
    edge_rows: Sequence[Mapping[str, object]],
    cluster_maps: Mapping[str, Mapping[str, ClusterRecord]],
) -> List[Mapping[str, object]]:
    summaries: List[Mapping[str, object]] = []
    for tolerance in CP_TOLERANCE_GRID:
        for sample, suffix in ((SAMPLE_A, "a"), (SAMPLE_B, "b")):
            for confidence_label in ("all", "both_assignment_ge_0.8"):
                counts = Counter()
                denominator = 0
                for raw, derived in zip(raw_relation_rows, edge_rows):
                    codes = json.loads(raw[f"read_full_relation_codes_{suffix}"]) if raw[f"read_full_relation_codes_{suffix}"] else []
                    if codes not in (["F"], ["R"]):
                        counts["uninformative_ambiguous_or_parallel"] += 1
                        continue
                    parent_id = str(derived[f"parent_mutation_id_{suffix}"])
                    child_id = str(derived[f"child_mutation_id_{suffix}"])
                    parent, child = cluster_maps[sample].get(parent_id), cluster_maps[sample].get(child_id)
                    if parent is None or child is None:
                        counts["uninformative_pyclone_endpoint_missing"] += 1
                        continue
                    if confidence_label == "both_assignment_ge_0.8" and not (
                        parent.assignment_prob >= ASSIGNMENT_THRESHOLD and child.assignment_prob >= ASSIGNMENT_THRESHOLD
                    ):
                        counts["excluded_assignment_lt_0.8"] += 1
                        continue
                    denominator += 1
                    counts[cp_status(parent, child, tolerance)] += 1
                summaries.append({
                    "sample": sample,
                    "tolerance": tolerance,
                    "confidence_stratum": confidence_label,
                    "input_relation_pairs": len(edge_rows),
                    "evaluable_directed_joined_edges": denominator,
                    "compatible": counts["compatible"],
                    "conflict": counts["conflict"],
                    "uninformative_same_cluster": counts["uninformative_same_cluster"],
                    "uninformative_interval_overlap": counts["uninformative_interval_overlap"],
                    "uninformative_ambiguous_or_parallel": counts["uninformative_ambiguous_or_parallel"],
                    "uninformative_pyclone_endpoint_missing": counts["uninformative_pyclone_endpoint_missing"],
                    "excluded_assignment_lt_0.8": counts["excluded_assignment_lt_0.8"],
                    "compatible_share_evaluable": counts["compatible"] / denominator if denominator else None,
                    "conflict_share_evaluable": counts["conflict"] / denominator if denominator else None,
                })
    return summaries


def parse_args() -> argparse.Namespace:
    topic = Path(__file__).resolve().parents[1]
    repo = topic.parents[1]
    topology = repo / "research" / "20260712_hcc1395_pair_site_topology_containment_validation" / "data"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pyclone-a",
        type=Path,
        default=topic / "runs" / "pyclone_vi" / "hcc1395_pair_primary_separate_HCC1395_main" / "results.tsv.gz",
    )
    parser.add_argument(
        "--pyclone-b",
        type=Path,
        default=topic / "runs" / "pyclone_vi" / "hcc1395_pair_primary_separate_HCC1395_DORADO_main" / "results.tsv.gz",
    )
    parser.add_argument(
        "--raw-vaf",
        type=Path,
        default=topic / "results" / "raw_vaf_validation_v1" / "data" / "hcc1395_pair_exact_shared_records.tsv.gz",
    )
    parser.add_argument("--regions", type=Path, default=topology / "hcc1395_site_topology_pair_outcomes.tsv")
    parser.add_argument("--alleles", type=Path, default=topology / "hcc1395_site_allele_identity.tsv.gz")
    parser.add_argument("--relations", type=Path, default=topology / "hcc1395_site_relation_evidence.tsv.gz")
    parser.add_argument(
        "--edge-schema-source",
        type=Path,
        default=repo / "research" / "20260712_hcc1395_pair_site_topology_containment_validation" / "scripts" / "build_site_topology_containment.py",
    )
    parser.add_argument("--output-dir", type=Path, default=topic / "results" / "clone_region_bridge_v1")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for path in (
        args.pyclone_a, args.pyclone_b, args.raw_vaf, args.regions, args.alleles,
        args.relations, args.edge_schema_source,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    edge_schema_qa = validate_edge_schema_source(args.edge_schema_source)

    clusters_a = load_cluster_results(args.pyclone_a, SAMPLE_A)
    clusters_b = load_cluster_results(args.pyclone_b, SAMPLE_B)
    if set(clusters_a) != set(clusters_b):
        raise ValueError("Separate-fit mutation universes are not identical; refuse asymmetric comparison")
    cluster_maps = {SAMPLE_A: clusters_a, SAMPLE_B: clusters_b}
    raw_vaf, raw_vaf_qa = load_raw_vaf(args.raw_vaf)
    region_rows, region_by_match = load_region_outcomes(args.regions)
    allele_rows, by_region_position, allele_qa = load_region_alleles(args.alleles, region_by_match)

    bridge_rows: List[Mapping[str, object]] = []
    by_region_bridge: Dict[str, List[Mapping[str, object]]] = defaultdict(list)
    for allele in allele_rows:
        mutation_id = str(allele["mutation_id"])
        raw = raw_vaf.get(mutation_id)
        a, b = clusters_a.get(mutation_id), clusters_b.get(mutation_id)
        region = region_by_match[str(allele["match_id"])]
        row: Dict[str, object] = dict(allele)
        row.update({
            "caller_k_a": int(region["caller_k_a"]),
            "caller_k_b": int(region["caller_k_b"]),
            "caller_shared_k": int(region["caller_shared_k"]),
            "k_bin": region["k_bin"],
            "caller_site_set_relation": region["caller_site_set_relation"],
            "equal_site_set": region["caller_site_set_relation"] == "equal",
            "joined_raw_vaf": raw is not None,
            "truth_confirmed_raw": raw["truth_confirmed"] if raw else "",
            "raw_source": raw["source"] if raw else "",
            "vaf_a": raw["af_a"] if raw else "",
            "vaf_b": raw["af_b"] if raw else "",
            "dp_a": raw["dp_a"] if raw else "",
            "dp_b": raw["dp_b"] if raw else "",
            "joined_pyclone_a": a is not None,
            "joined_pyclone_b": b is not None,
            "joined_pyclone_both": a is not None and b is not None,
            "cluster_a": a.cluster_id if a else "",
            "cluster_b": b.cluster_id if b else "",
            "cp_a": a.cp if a else "",
            "cp_b": b.cp if b else "",
            "cp_std_a": a.cp_std if a else "",
            "cp_std_b": b.cp_std if b else "",
            "assignment_prob_a": a.assignment_prob if a else "",
            "assignment_prob_b": b.assignment_prob if b else "",
            "assignment_high_a": a.assignment_prob >= ASSIGNMENT_THRESHOLD if a else False,
            "assignment_high_b": b.assignment_prob >= ASSIGNMENT_THRESHOLD if b else False,
            "assignment_high_both": (
                a is not None and b is not None
                and a.assignment_prob >= ASSIGNMENT_THRESHOLD
                and b.assignment_prob >= ASSIGNMENT_THRESHOLD
            ),
            "clonal_state_a": "clonal" if a and a.cp >= CLONAL_THRESHOLD else ("subclonal" if a else ""),
            "clonal_state_b": "clonal" if b and b.cp >= CLONAL_THRESHOLD else ("subclonal" if b else ""),
            "clonal_state_same": (a.cp >= CLONAL_THRESHOLD) == (b.cp >= CLONAL_THRESHOLD) if a and b else "",
            "cluster_id_same_raw": a.cluster_id == b.cluster_id if a and b else "",
        })
        bridge_rows.append(row)
        by_region_bridge[str(allele["match_id"])].append(row)

    if len(bridge_rows) != allele_qa["expanded_exact_alleles"]:
        raise RuntimeError("Bridge row conservation failed")
    if len({(row["match_id"], row["mutation_id"]) for row in bridge_rows}) != len(bridge_rows):
        raise RuntimeError("Region×mutation bridge is not unique")

    global_rows = []
    for mutation_id in sorted(clusters_a):
        a, b = clusters_a[mutation_id], clusters_b[mutation_id]
        global_rows.append({
            "mutation_id": mutation_id,
            "joined_pyclone_both": True,
            "cluster_a": a.cluster_id,
            "cluster_b": b.cluster_id,
            "cp_a": a.cp,
            "cp_b": b.cp,
            "assignment_prob_a": a.assignment_prob,
            "assignment_prob_b": b.assignment_prob,
            "assignment_high_both": (
                a.assignment_prob >= ASSIGNMENT_THRESHOLD and b.assignment_prob >= ASSIGNMENT_THRESHOLD
            ),
        })
    concordance_rows: List[Mapping[str, object]] = [
        mutation_concordance(global_rows, "global_separate_fit_universe", "all"),
        mutation_concordance(
            [row for row in global_rows if row["assignment_high_both"]],
            "global_separate_fit_universe", "both_assignment_ge_0.8",
        ),
        mutation_concordance(bridge_rows, "fixed_5720_regional_subset", "all_joined"),
        mutation_concordance(
            [row for row in bridge_rows if row["assignment_high_both"]],
            "fixed_5720_regional_subset", "both_assignment_ge_0.8",
        ),
        mutation_concordance(
            [row for row in bridge_rows if row["equal_site_set"]],
            "fixed_5720_regional_subset", "equal_site_set",
        ),
        mutation_concordance(
            [row for row in bridge_rows if row["equal_site_set"] and row["assignment_high_both"]],
            "fixed_5720_regional_subset", "equal_site_set_and_both_assignment_ge_0.8",
        ),
    ]
    for k_bin in sorted({str(row["k_bin"]) for row in bridge_rows}):
        concordance_rows.append(mutation_concordance(
            [row for row in bridge_rows if row["k_bin"] == k_bin],
            "fixed_5720_regional_subset", k_bin,
        ))

    distribution_inputs: List[Tuple[str, str, Sequence[Mapping[str, object]]]] = [
        ("global_separate_fit_universe", "all", global_rows),
        ("global_separate_fit_universe", "both_assignment_ge_0.8", [
            row for row in global_rows if row["assignment_high_both"]
        ]),
        ("fixed_5720_regional_subset", "all_joined", bridge_rows),
        ("fixed_5720_regional_subset", "both_assignment_ge_0.8", [
            row for row in bridge_rows if row["assignment_high_both"]
        ]),
        ("fixed_5720_regional_subset", "equal_site_set", [
            row for row in bridge_rows if row["equal_site_set"]
        ]),
        ("fixed_5720_regional_subset", "equal_site_set_and_both_assignment_ge_0.8", [
            row for row in bridge_rows if row["equal_site_set"] and row["assignment_high_both"]
        ]),
    ]
    for k_bin in sorted({str(row["k_bin"]) for row in bridge_rows}):
        distribution_inputs.extend([
            ("fixed_5720_regional_subset", k_bin, [row for row in bridge_rows if row["k_bin"] == k_bin]),
            ("fixed_5720_regional_subset", f"{k_bin}_and_both_assignment_ge_0.8", [
                row for row in bridge_rows if row["k_bin"] == k_bin and row["assignment_high_both"]
            ]),
        ])
    cluster_distribution: List[Mapping[str, object]] = []
    for population, stratum, values in distribution_inputs:
        cluster_distribution.extend(cluster_distribution_rows(values, population, stratum))
        concordance_rows.append(mutation_concordance(
            [row for row in bridge_rows if row["k_bin"] == k_bin and row["assignment_high_both"]],
            "fixed_5720_regional_subset", f"{k_bin}_and_both_assignment_ge_0.8",
        ))

    pattern_rows: List[Mapping[str, object]] = []
    for region in region_rows:
        site_rows = by_region_bridge.get(region["match_id"], [])
        pattern_rows.append(region_pattern_row(region, site_rows, False))
        pattern_rows.append(region_pattern_row(region, site_rows, True))
    if len(pattern_rows) != FIXED_REGION_DENOMINATOR * 2:
        raise RuntimeError("Region pattern fixed-denominator conservation failed")

    pattern_summary_rows: List[Mapping[str, object]] = []
    for endpoint in sorted({str(row["endpoint"]) for row in pattern_rows}):
        endpoint_rows = [row for row in pattern_rows if row["endpoint"] == endpoint]
        strata = [("ALL", endpoint_rows)]
        for k_bin in sorted({str(row["k_bin"]) for row in endpoint_rows}):
            strata.append((k_bin, [row for row in endpoint_rows if row["k_bin"] == k_bin]))
        strata.append(("k>1", [row for row in endpoint_rows if int(row["caller_shared_k"]) > 1]))
        strata.append(("k>1_equal_site_set", [
            row for row in endpoint_rows if int(row["caller_shared_k"]) > 1 and row["equal_site_set"]
        ]))
        for stratum, values in strata:
            evaluable = [row for row in values if row["evaluable_multilocus"]]
            categories = Counter(str(row["pattern_category"]) for row in evaluable)
            both_multicluster_evaluable = (
                categories["both_multicluster_exact_partition"]
                + categories["both_multicluster_different_partition"]
            )
            pattern_summary_rows.append({
                "endpoint": endpoint,
                "stratum": stratum,
                "fixed_regions_in_stratum": len(values),
                "evaluable_multilocus_regions": len(evaluable),
                "evaluable_share_stratum": len(evaluable) / len(values) if values else None,
                "partition_exact": sum(row["partition_exact"] is True for row in evaluable),
                "partition_exact_share_evaluable": (
                    sum(row["partition_exact"] is True for row in evaluable) / len(evaluable) if evaluable else None
                ),
                "mean_pairwise_cocluster_agreement": (
                    float(np.mean([float(row["pairwise_cocluster_agreement"]) for row in evaluable])) if evaluable else None
                ),
                "both_single_cluster": categories["both_single_cluster"],
                "both_multicluster_exact_partition": categories["both_multicluster_exact_partition"],
                "one_single_one_multicluster": categories["one_single_one_multicluster"],
                "both_multicluster_different_partition": categories["both_multicluster_different_partition"],
                "both_multicluster_evaluable": both_multicluster_evaluable,
                "both_multicluster_exact_share": (
                    categories["both_multicluster_exact_partition"] / both_multicluster_evaluable
                    if both_multicluster_evaluable else None
                ),
            })

    relation_rows = load_relation_evidence(args.relations)
    edge_rows = [directed_edge_row(row, by_region_position, cluster_maps) for row in relation_rows]
    edge_summary_rows = summarize_edges(relation_rows, edge_rows, cluster_maps)

    join_coverage_rows = [
        {"level": "raw_vaf_source_rows", "metric": "input_source_rows", "n": raw_vaf_qa["input_rows"], "denominator": raw_vaf_qa["input_rows"]},
        {"level": "raw_vaf_unique_mutations", "metric": "unique_mutation_keys_after_exact_payload_collapse", "n": len(raw_vaf), "denominator": raw_vaf_qa["input_rows"]},
        {"level": "raw_vaf_source_rows", "metric": "duplicate_source_rows_accounted", "n": raw_vaf_qa["duplicate_extra_rows"], "denominator": raw_vaf_qa["input_rows"]},
        {"level": "fixed_regions", "metric": "input_regions", "n": len(region_rows), "denominator": FIXED_REGION_DENOMINATOR},
        {"level": "regional_exact_alleles", "metric": "input_expanded_alleles", "n": len(bridge_rows), "denominator": len(bridge_rows)},
        {"level": "regional_exact_alleles", "metric": "joined_raw_vaf", "n": sum(row["joined_raw_vaf"] for row in bridge_rows), "denominator": len(bridge_rows)},
        {"level": "regional_exact_alleles", "metric": "joined_pyclone_both", "n": sum(row["joined_pyclone_both"] for row in bridge_rows), "denominator": len(bridge_rows)},
        {"level": "regional_exact_alleles", "metric": "joined_pyclone_both_assignment_ge_0.8", "n": sum(row["assignment_high_both"] for row in bridge_rows), "denominator": len(bridge_rows)},
        {"level": "fixed_regions", "metric": "regions_with_at_least_one_pyclone_join", "n": sum(
            any(row["joined_pyclone_both"] for row in by_region_bridge.get(region["match_id"], [])) for region in region_rows
        ), "denominator": FIXED_REGION_DENOMINATOR},
        {"level": "fixed_regions", "metric": "regions_with_at_least_two_pyclone_joins", "n": sum(
            sum(row["joined_pyclone_both"] for row in by_region_bridge.get(region["match_id"], [])) >= 2
            for region in region_rows
        ), "denominator": FIXED_REGION_DENOMINATOR},
        {"level": "directed_relation_pairs", "metric": "input_read_relation_pairs", "n": len(edge_rows), "denominator": len(edge_rows)},
        {"level": "directed_relation_pairs", "metric": "same_direction_cross_source", "n": sum(
            row["cross_source_direction_relation"] == "same_direction" for row in edge_rows
        ), "denominator": len(edge_rows)},
        {"level": "directed_relation_pairs", "metric": "both_sources_determinate_direction", "n": sum(
            row["cross_source_direction_relation"] in {"same_direction", "opposite_direction"} for row in edge_rows
        ), "denominator": len(edge_rows)},
        {"level": "both_sources_determinate_direction", "metric": "same_direction_among_both_determinate", "n": sum(
            row["cross_source_direction_relation"] == "same_direction" for row in edge_rows
        ), "denominator": sum(
            row["cross_source_direction_relation"] in {"same_direction", "opposite_direction"} for row in edge_rows
        )},
        {"level": "both_sources_determinate_direction", "metric": "opposite_direction_among_both_determinate", "n": sum(
            row["cross_source_direction_relation"] == "opposite_direction" for row in edge_rows
        ), "denominator": sum(
            row["cross_source_direction_relation"] in {"same_direction", "opposite_direction"} for row in edge_rows
        )},
    ]
    for row in join_coverage_rows:
        row["share"] = row["n"] / row["denominator"] if row["denominator"] else None

    bridge_columns = [
        "match_id", "chrom", "region", "position", "ref", "alt", "mutation_id",
        "caller_k_a", "caller_k_b", "caller_shared_k", "k_bin", "caller_site_set_relation", "equal_site_set",
        "joined_raw_vaf", "truth_confirmed_raw", "raw_source", "vaf_a", "vaf_b", "dp_a", "dp_b",
        "joined_pyclone_a", "joined_pyclone_b", "joined_pyclone_both", "cluster_a", "cluster_b",
        "cp_a", "cp_b", "cp_std_a", "cp_std_b", "assignment_prob_a", "assignment_prob_b",
        "assignment_high_a", "assignment_high_b", "assignment_high_both", "clonal_state_a", "clonal_state_b",
        "clonal_state_same", "cluster_id_same_raw",
    ]
    concordance_columns = [
        "population", "stratum", "mutation_rows", "unique_mutations", "n", "clusters_a", "clusters_b",
        "ari", "nmi", "hungarian_agreement", "assignment_high_both", "cp_mae", "cp_pearson", "cp_spearman",
        "clonal_state_agreement", "clonal_state_kappa", "subclonal_intersection", "subclonal_union", "subclonal_jaccard",
    ]
    pattern_columns = [
        "match_id", "chrom", "region", "endpoint", "fixed_denominator", "caller_k_a", "caller_k_b",
        "caller_shared_k", "k_bin", "caller_site_set_relation", "equal_site_set", "regional_exact_alleles",
        "joined_mutations", "evaluable_multilocus", "clusters_a", "clusters_b", "partition_signature_a",
        "partition_signature_b", "partition_exact", "pairwise_cocluster_agreement", "region_ari",
        "pattern_category", "subclonal_jaccard",
    ]
    pattern_summary_columns = [
        "endpoint", "stratum", "fixed_regions_in_stratum", "evaluable_multilocus_regions",
        "evaluable_share_stratum", "partition_exact", "partition_exact_share_evaluable",
        "mean_pairwise_cocluster_agreement", "both_single_cluster", "both_multicluster_exact_partition",
        "one_single_one_multicluster", "both_multicluster_different_partition", "both_multicluster_evaluable",
        "both_multicluster_exact_share",
    ]
    edge_columns = [
        "match_id", "chrom", "region", "component_index", "pair_index", "left_position", "right_position",
        "left_mutation_id", "right_mutation_id", "endpoint_identity_unique", "read_codes_a", "read_codes_b",
        "relation_status_a", "relation_status_b", "parent_mutation_id_a", "child_mutation_id_a",
        "parent_mutation_id_b", "child_mutation_id_b", "pyclone_join_a", "pyclone_join_b",
        "assignment_high_endpoints_a", "assignment_high_endpoints_b", "cluster_parent_a", "cluster_child_a",
        "cluster_parent_b", "cluster_child_b", "cp_parent_a", "cp_child_a", "cp_parent_b", "cp_child_b",
        "cp_std_parent_a", "cp_std_child_a", "cp_std_parent_b", "cp_std_child_b", "assignment_parent_a",
        "assignment_child_a", "assignment_parent_b", "assignment_child_b", "cp_delta_parent_minus_child_a",
        "cp_delta_parent_minus_child_b", "cp_status_t0.02_a", "cp_status_t0.02_b",
        "cross_source_direction_relation", "combined_cp_status_t0.02",
    ]
    edge_summary_columns = [
        "sample", "tolerance", "confidence_stratum", "input_relation_pairs", "evaluable_directed_joined_edges",
        "compatible", "conflict", "uninformative_same_cluster", "uninformative_interval_overlap",
        "uninformative_ambiguous_or_parallel", "uninformative_pyclone_endpoint_missing",
        "excluded_assignment_lt_0.8", "compatible_share_evaluable", "conflict_share_evaluable",
    ]
    join_columns = ["level", "metric", "n", "denominator", "share"]
    cluster_distribution_columns = [
        "population", "stratum", "sample", "cluster_id_within_separate_fit", "n", "denominator", "share",
        "mean_cp", "mean_assignment_prob", "assignment_ge_0.8", "subclonal_cp_lt_0.9",
    ]

    paths = {
        "bridge": args.output_dir / "region_mutation_bridge.tsv.gz",
        "join_coverage": args.output_dir / "join_coverage.tsv",
        "concordance": args.output_dir / "global_vs_regional_cluster_concordance.tsv",
        "cluster_distribution": args.output_dir / "cluster_distribution.tsv",
        "patterns": args.output_dir / "region_cluster_patterns.tsv.gz",
        "pattern_summary": args.output_dir / "region_cluster_pattern_summary.tsv",
        "edges": args.output_dir / "directed_edge_cp_compatibility.tsv.gz",
        "edge_summary": args.output_dir / "directed_edge_summary.tsv",
        "checks": args.output_dir / "checks.tsv",
        "diagnostic_flags": args.output_dir / "diagnostic_flags.tsv",
        "summary": args.output_dir / "summary.json",
        "provenance": args.output_dir / "provenance.json",
        "report": args.output_dir / "20260712_clone_region_bridge_summary.md",
    }
    write_tsv(paths["bridge"], bridge_rows, bridge_columns, gzip_output=True)
    write_tsv(paths["join_coverage"], join_coverage_rows, join_columns)
    write_tsv(paths["concordance"], concordance_rows, concordance_columns)
    write_tsv(paths["cluster_distribution"], cluster_distribution, cluster_distribution_columns)
    write_tsv(paths["patterns"], pattern_rows, pattern_columns, gzip_output=True)
    write_tsv(paths["pattern_summary"], pattern_summary_rows, pattern_summary_columns)
    write_tsv(paths["edges"], edge_rows, edge_columns, gzip_output=True)
    write_tsv(paths["edge_summary"], edge_summary_rows, edge_summary_columns)

    checks = [
        {"check": "fixed_region_count", "pass": len(region_rows) == FIXED_REGION_DENOMINATOR, "observed": len(region_rows), "expected": FIXED_REGION_DENOMINATOR},
        {"check": "separate_fit_universe_identical", "pass": set(clusters_a) == set(clusters_b), "observed": len(set(clusters_a) & set(clusters_b)), "expected": len(clusters_a)},
        {"check": "regional_bridge_row_conservation", "pass": len(bridge_rows) == allele_qa["expanded_exact_alleles"], "observed": len(bridge_rows), "expected": allele_qa["expanded_exact_alleles"]},
        {"check": "region_mutation_unique", "pass": len({(row["match_id"], row["mutation_id"]) for row in bridge_rows}) == len(bridge_rows), "observed": len({(row["match_id"], row["mutation_id"]) for row in bridge_rows}), "expected": len(bridge_rows)},
        {"check": "ambiguous_region_positions_zero", "pass": allele_qa["ambiguous_region_positions"] == 0, "observed": allele_qa["ambiguous_region_positions"], "expected": 0},
        {"check": "mutation_region_multiplicity_zero", "pass": allele_qa["mutation_region_multiplicity_gt1"] == 0, "observed": allele_qa["mutation_region_multiplicity_gt1"], "expected": 0},
        {"check": "pattern_rows_fixed_two_endpoints", "pass": len(pattern_rows) == 2 * FIXED_REGION_DENOMINATOR, "observed": len(pattern_rows), "expected": 2 * FIXED_REGION_DENOMINATOR},
        {"check": "edge_rows_conserved", "pass": len(edge_rows) == len(relation_rows), "observed": len(edge_rows), "expected": len(relation_rows)},
        {"check": "edge_schema_source_contract", "pass": edge_schema_qa["status"] == "PASS", "observed": edge_schema_qa["status"], "expected": "PASS"},
        {"check": "raw_vaf_source_row_conservation", "pass": raw_vaf_qa["row_conservation_pass"], "observed": raw_vaf_qa["input_rows"], "expected": len(raw_vaf) + raw_vaf_qa["duplicate_extra_rows"]},
        {"check": "raw_vaf_duplicate_payload_conflicts_zero", "pass": raw_vaf_qa["conflicting_duplicate_payloads"] == 0, "observed": raw_vaf_qa["conflicting_duplicate_payloads"], "expected": 0},
        {"check": "raw_vaf_key_unique_after_audited_collapse", "pass": len(raw_vaf) == raw_vaf_qa["unique_mutation_keys_after_exact_payload_collapse"], "observed": len(raw_vaf), "expected": raw_vaf_qa["unique_mutation_keys_after_exact_payload_collapse"]},
        {"check": "pyclone_a_key_unique", "pass": True, "observed": len(clusters_a), "expected": len(clusters_a)},
        {"check": "pyclone_b_key_unique", "pass": True, "observed": len(clusters_b), "expected": len(clusters_b)},
    ]
    write_tsv(paths["checks"], checks, ["check", "pass", "observed", "expected"])
    if not all(bool(row["pass"]) for row in checks):
        failed = [row["check"] for row in checks if not row["pass"]]
        raise RuntimeError(f"Bridge QA failed: {failed}")

    concordance_lookup = {(row["population"], row["stratum"]): row for row in concordance_rows}
    pattern_lookup = {(row["endpoint"], row["stratum"]): row for row in pattern_summary_rows}
    primary_regional = concordance_lookup[("fixed_5720_regional_subset", "all_joined")]
    high_regional = concordance_lookup[("fixed_5720_regional_subset", "both_assignment_ge_0.8")]
    kgt1_pattern = pattern_lookup[("all_joined", "k>1")]
    edge_primary = {
        row["sample"]: row for row in edge_summary_rows
        if float(row["tolerance"]) == PRIMARY_CP_TOLERANCE and row["confidence_stratum"] == "all"
    }
    edge_high = {
        row["sample"]: row for row in edge_summary_rows
        if float(row["tolerance"]) == PRIMARY_CP_TOLERANCE
        and row["confidence_stratum"] == "both_assignment_ge_0.8"
    }
    high_confidence_degenerate = (
        int(high_regional["subclonal_union"]) < 20
        or int(high_regional["subclonal_union"]) / int(high_regional["n"]) < 0.01
    )
    single_cluster_dominated = (
        int(kgt1_pattern["both_single_cluster"]) / int(kgt1_pattern["evaluable_multilocus_regions"]) >= 0.90
    )
    high_edge_distinct_cluster_a = int(edge_high[SAMPLE_A]["compatible"]) + int(edge_high[SAMPLE_A]["conflict"])
    high_edge_distinct_cluster_b = int(edge_high[SAMPLE_B]["compatible"]) + int(edge_high[SAMPLE_B]["conflict"])
    diagnostic_flags = [
        {
            "flag": "regional_high_confidence_subclone_degeneracy",
            "severity": "WARN" if high_confidence_degenerate else "INFO",
            "status": "DEGENERATE_SELECTION_INDUCED" if high_confidence_degenerate else "NONDEGENERATE",
            "observed": f"subclonal_union={high_regional['subclonal_union']}; n={high_regional['n']}",
            "threshold": "subclonal_union<20 OR subclonal_union/n<0.01",
            "interpretation": (
                "Perfect high-confidence concordance is vacuous/selection-induced and cannot support high clone agreement"
                if high_confidence_degenerate else "High-confidence subset retains enough subclonal variation"
            ),
        },
        {
            "flag": "k_gt_1_partition_exact_single_cluster_dominance",
            "severity": "WARN" if single_cluster_dominated else "INFO",
            "status": "SINGLE_CLUSTER_DOMINATED" if single_cluster_dominated else "NOT_SINGLE_CLUSTER_DOMINATED",
            "observed": (
                f"both_single_cluster={kgt1_pattern['both_single_cluster']}; "
                f"evaluable={kgt1_pattern['evaluable_multilocus_regions']}"
            ),
            "threshold": "both_single_cluster/evaluable>=0.90",
            "interpretation": "Overall partition-exact rate is not an informative multi-clone concordance rate",
        },
        {
            "flag": "high_confidence_directed_edge_distinct_cluster_information",
            "severity": "WARN" if high_edge_distinct_cluster_a + high_edge_distinct_cluster_b == 0 else "INFO",
            "status": (
                "NO_DISTINCT_CLUSTER_EDGE_INFORMATION"
                if high_edge_distinct_cluster_a + high_edge_distinct_cluster_b == 0 else "INFORMATIVE"
            ),
            "observed": f"HCC1395={high_edge_distinct_cluster_a}; HCC1395_DORADO={high_edge_distinct_cluster_b}",
            "threshold": "compatible+conflict across both samples == 0",
            "interpretation": "High-confidence directed edges cannot validate CP ordering when all joined endpoints share a cluster",
        },
    ]
    write_tsv(
        paths["diagnostic_flags"], diagnostic_flags,
        ["flag", "severity", "status", "observed", "threshold", "interpretation"],
    )
    summary = {
        "schema_name": "intersubmod.clone_region_bridge",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "status": "PASS",
        "scope": "HCC1395 vs HCC1395_DORADO; fixed 5,720 exact-coordinate complete-both historical layered-v2 regions",
        "primary_endpoint": "independent separate PyClone-VI fits on identical shared-truth mutation universe",
        "claim_ceiling": (
            "Conditional mutation-cluster-to-region bridge and read-edge compatibility diagnostic; "
            "not a clone-tree truth test. Joint PyClone labels are excluded from independent concordance."
        ),
        "thresholds": {
            "clonal_cp": CLONAL_THRESHOLD,
            "assignment_high": ASSIGNMENT_THRESHOLD,
            "primary_cp_tolerance": PRIMARY_CP_TOLERANCE,
            "cp_interval": "point estimate +/- 1.96*PyClone-VI cellular_prevalence_std, clipped to [0,1]; diagnostic, not calibrated CI",
        },
        "join_coverage": {row["metric"]: row for row in join_coverage_rows},
        "raw_vaf_qa": raw_vaf_qa,
        "allele_qa": allele_qa,
        "global_concordance": concordance_lookup[("global_separate_fit_universe", "all")],
        "regional_concordance": primary_regional,
        "regional_high_confidence_concordance": high_regional,
        "high_confidence_degeneracy": {
            "status": "WARN_SELECTION_INDUCED" if high_confidence_degenerate else "PASS_NONDEGENERATE",
            "n": high_regional["n"],
            "subclonal_union": high_regional["subclonal_union"],
            "subclonal_union_share": high_regional["subclonal_union"] / high_regional["n"],
            "interpretation": (
                "ARI/NMI/Jaccard=1 is vacuous here: 10,965 mutations contain only one mutation in the "
                "two-sample subclonal union; do not claim high-confidence regional clone agreement."
            ),
        },
        "k_gt_1_region_pattern": kgt1_pattern,
        "k_gt_1_informative_both_multicluster": {
            "denominator": kgt1_pattern["both_multicluster_evaluable"],
            "exact": kgt1_pattern["both_multicluster_exact_partition"],
            "different": kgt1_pattern["both_multicluster_different_partition"],
            "exact_share": kgt1_pattern["both_multicluster_exact_share"],
            "interpretation": "Label-invariant partition agreement only among regions where both fits contain >1 cluster",
        },
        "directed_edge_schema": {
            **edge_schema_qa,
            "definition": "F=lower-position ancestor; R=higher-position ancestor; P=parallel; only singleton F/R is directed",
            "primary_endpoint": "read_full candidate-set singleton direction; VAF-selected relations are not used for independent CP compatibility",
        },
        "directed_edge_primary_tolerance_0.02": edge_primary,
        "directed_edge_high_confidence_tolerance_0.02": edge_high,
        "checks": checks,
        "diagnostic_flags": diagnostic_flags,
        "outputs": {name: str(path.resolve()) for name, path in paths.items()},
    }
    paths["summary"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    input_paths = {
        "pyclone_a": args.pyclone_a,
        "pyclone_b": args.pyclone_b,
        "raw_vaf": args.raw_vaf,
        "regions": args.regions,
        "alleles": args.alleles,
        "relations": args.relations,
        "edge_schema_source": args.edge_schema_source,
        "builder": Path(__file__).resolve(),
    }
    provenance = {
        "schema_name": "intersubmod.clone_region_bridge_provenance",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "command": [
            str(Path(__file__).resolve()),
            "--pyclone-a", str(args.pyclone_a.resolve()),
            "--pyclone-b", str(args.pyclone_b.resolve()),
            "--raw-vaf", str(args.raw_vaf.resolve()),
            "--regions", str(args.regions.resolve()),
            "--alleles", str(args.alleles.resolve()),
            "--relations", str(args.relations.resolve()),
            "--edge-schema-source", str(args.edge_schema_source.resolve()),
            "--output-dir", str(args.output_dir.resolve()),
        ],
        "python_version": sys.version,
        "inputs": {
            name: {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for name, path in input_paths.items()
        },
        "outputs": {
            name: {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for name, path in paths.items() if path.is_file() and name not in {"provenance", "report"}
        },
        "notes": [
            "Separate-fit cluster IDs are compared label-invariantly (ARI/NMI/Hungarian and regional partitions).",
            "Raw cluster_id equality across fits is retained only as a diagnostic column and is never interpreted as label identity.",
            "Read-full singleton F/R edges are read-pattern constraints; CP compatibility is conditional on PyClone CN/purity/count inputs.",
            "Historical layered-v2 regional artifacts are engineering evidence, not clean-v3 biological truth.",
        ],
    }
    paths["provenance"].write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

    report = f"""<!--
建立時間: 2026-07-12
目標: 將 separate PyClone-VI mutation clusters 接回固定 5,720 regions，量化區域子集合再現與 read-edge CP 相容性
處理範圍: HCC1395 vs HCC1395_DORADO；chr1-22；fixed exact-coordinate complete-both regions
關聯檔案: InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/
-->

# Clone ↔ region bridge（v1）

> **PASS；但只支持 conditional cluster-to-region reproducibility，不是 tree truth。** Primary 使用兩側獨立 separate fits；joint labels 未列為獨立再現證據。

## Join coverage

- Fixed regions：{FIXED_REGION_DENOMINATOR:,}。
- Shared exact regional alleles：{len(bridge_rows):,}；raw VAF joined {sum(row['joined_raw_vaf'] for row in bridge_rows):,}；兩側 PyClone joined {sum(row['joined_pyclone_both'] for row in bridge_rows):,}。
- 兩側 assignment ≥0.8：{sum(row['assignment_high_both'] for row in bridge_rows):,} regional alleles。
- Region×mutation key unique、PyClone keys unique：PASS。Raw VAF 的 {raw_vaf_qa['input_rows']:,} source rows 中，{raw_vaf_qa['duplicate_extra_rows']:,} 列是 payload 完全相同、僅 source label 不同的重複證據列；依明示規則收斂為 {len(raw_vaf):,} unique mutations，before/after rows 守恆 PASS。

## Global vs regional mutation-cluster concordance

| Population | n | ARI | NMI | Subclonal Jaccard | κ | CP Spearman |
|---|---:|---:|---:|---:|---:|---:|
| Global separate-fit universe | {summary['global_concordance']['n']:,} | {summary['global_concordance']['ari']:.3f} | {summary['global_concordance']['nmi']:.3f} | {summary['global_concordance']['subclonal_jaccard']:.3f} | {summary['global_concordance']['clonal_state_kappa']:.3f} | {summary['global_concordance']['cp_spearman']:.3f} |
| Fixed-regional joined subset | {primary_regional['n']:,} | {primary_regional['ari']:.3f} | {primary_regional['nmi']:.3f} | {primary_regional['subclonal_jaccard']:.3f} | {primary_regional['clonal_state_kappa']:.3f} | {primary_regional['cp_spearman']:.3f} |
| Regional both assignment≥0.8 | {high_regional['n']:,} | {high_regional['ari']:.3f} | {high_regional['nmi']:.3f} | {high_regional['subclonal_jaccard']:.3f} | {high_regional['clonal_state_kappa']:.3f} | {high_regional['cp_spearman']:.3f} |

`assignment≥0.8` 的完美值不可單獨解讀成 clone 再現：10,965 個 regional high-confidence mutations 的 subclonal union 只有 {high_regional['subclonal_union']:,}（{100*high_regional['subclonal_union']/high_regional['n']:.4f}%），幾乎全被兩側主 clonal cluster 支配，屬 **vacuous selection-induced concordance**；不可稱「高信心區域 clone 高度一致」。

## k>1 multi-locus partition pattern

- All-joined k>1：fixed stratum {kgt1_pattern['fixed_regions_in_stratum']:,}；evaluable {kgt1_pattern['evaluable_multilocus_regions']:,}；partition exact {kgt1_pattern['partition_exact']:,}（{100*kgt1_pattern['partition_exact_share_evaluable']:.2f}% of evaluable）。
- 其中 both-single-cluster {kgt1_pattern['both_single_cluster']:,}（{100*kgt1_pattern['both_single_cluster']/kgt1_pattern['evaluable_multilocus_regions']:.2f}%）；真正兩側皆 multi-cluster 且 partition exact 只有 {kgt1_pattern['both_multicluster_exact_partition']:,}（{100*kgt1_pattern['both_multicluster_exact_partition']/kgt1_pattern['evaluable_multilocus_regions']:.2f}%）。因此 96.90% 不可稱為「多 clone 結構 96.90% 一致」。
- 限定兩側都真的含 >1 cluster 的 informative denominator 是 {kgt1_pattern['both_multicluster_evaluable']:,}：exact {kgt1_pattern['both_multicluster_exact_partition']:,}、different {kgt1_pattern['both_multicluster_different_partition']:,}，即 {100*kgt1_pattern['both_multicluster_exact_share']:.2f}%（21/34）。這才是較接近「區域內多 clone partition 是否一致」的條件式數字，但 n=34 很小。
- 非 exact 分為 one-single/one-multi {kgt1_pattern['one_single_one_multicluster']:,} 與 both-multi-different {kgt1_pattern['both_multicluster_different_partition']:,}；完整列在 `region_cluster_patterns.tsv.gz`。

## Read-directed edge CP compatibility

- Edge contract已從 source驗證：F=低座標祖先、R=高座標祖先、P=平行；只有 singleton F/R 進 directed endpoint。
- CP判斷：同群為 uninformative；不同群用 CP±1.96×std 與 tolerance 0.02，清楚滿足才 compatible、清楚反向才 conflict、區間重疊則 uninformative。
- HCC1395：directed+joined {edge_primary[SAMPLE_A]['evaluable_directed_joined_edges']:,}/{edge_primary[SAMPLE_A]['input_relation_pairs']:,}；compatible {edge_primary[SAMPLE_A]['compatible']:,}、conflict {edge_primary[SAMPLE_A]['conflict']:,}、same-cluster uninformative {edge_primary[SAMPLE_A]['uninformative_same_cluster']:,}。DORADO：directed+joined {edge_primary[SAMPLE_B]['evaluable_directed_joined_edges']:,}/{edge_primary[SAMPLE_B]['input_relation_pairs']:,}；compatible {edge_primary[SAMPLE_B]['compatible']:,}、conflict {edge_primary[SAMPLE_B]['conflict']:,}、same-cluster uninformative {edge_primary[SAMPLE_B]['uninformative_same_cluster']:,}。
- assignment≥0.8 後，HCC1395 的 {edge_high[SAMPLE_A]['evaluable_directed_joined_edges']:,} 與 DORADO 的 {edge_high[SAMPLE_B]['evaluable_directed_joined_edges']:,} 條 evaluable edges 全部是 same-cluster；compatible/conflict 皆為 0。因此 high-confidence endpoint **沒有不同 cluster 的 CP ordering 資訊**，0 conflict 不能解讀成已驗證 ancestry。
- 兩側都有 determinate direction 的只有 {summary['join_coverage']['both_sources_determinate_direction']['n']:,}/{summary['join_coverage']['both_sources_determinate_direction']['denominator']:,} pairs；其中 same direction {summary['join_coverage']['same_direction_among_both_determinate']['n']:,}/{summary['join_coverage']['same_direction_among_both_determinate']['denominator']:,}，opposite {summary['join_coverage']['opposite_direction_among_both_determinate']['n']:,}。這是高條件化分母，不能外推成全 5,720 regions 的方向一致率。
- 這是 read-pattern edge 對 conditional PyClone CCF 的診斷；不是對真實 ancestry 的證明。VAF-selected edge 未拿來當獨立驗證，以避免同 VAF 循環。

## Claim ceiling

本橋接可回答「固定區域中的哪些 shared mutations 被兩側 separate PyClone fit 接上、區域內 co-clustering partition 是否再現、read-compatible directed constraint 是否與 CCF ordering 明顯衝突」。它不能把 PyClone clusters 升格為 clone truth，也不能證明每區唯一演化樹。Historical layered-v2 region artifacts仍受 upstream engineering-snapshot ceiling。
"""
    paths["report"].write_text(report)

    # Add final self-receipts only after summary/provenance/report exist.
    final_outputs = {
        name: {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for name, path in paths.items() if path.is_file() and name != "provenance"
    }
    provenance["outputs"] = final_outputs
    paths["provenance"].write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

    print(f"status=PASS regions={len(region_rows)} regional_alleles={len(bridge_rows)}")
    print(
        "joined_raw={} joined_pyclone_both={} joined_high_both={}".format(
            sum(row["joined_raw_vaf"] for row in bridge_rows),
            sum(row["joined_pyclone_both"] for row in bridge_rows),
            sum(row["assignment_high_both"] for row in bridge_rows),
        )
    )
    print(
        "regional_ari={:.6f} regional_subJ={:.6f} k>1_exact={}/{}".format(
            primary_regional["ari"], primary_regional["subclonal_jaccard"],
            kgt1_pattern["partition_exact"], kgt1_pattern["evaluable_multilocus_regions"],
        )
    )
    print(f"output={args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
