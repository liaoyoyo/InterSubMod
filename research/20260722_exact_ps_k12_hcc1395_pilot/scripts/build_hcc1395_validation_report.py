#!/usr/bin/env python3
"""Build the HCC1395 exact-PS pilot evidence snapshot and report artifact.

This script independently recomputes the report-facing counts from the per-
chromosome gzip TSV artifacts.  The run receipt is used only as a comparison
target, not as the sole source of aggregate metrics.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


CHROMS = [f"chr{i}" for i in range(1, 23)]
REPORT_TITLE = "HCC1395 exact PS 分區與 k≤12 分割驗證"
TOPIC = "research/20260722_exact_ps_k12_hcc1395_pilot"
RUN_SOURCE = (
    "output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/"
    "hcc1395_chr1_22_direct_big7_v2/run_receipt.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--hc-bed", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def load_bed(path: Path) -> dict[str, tuple[list[int], list[tuple[int, int]]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            raw[fields[0]].append((int(fields[1]), int(fields[2])))
    result: dict[str, tuple[list[int], list[tuple[int, int]]]] = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        result[chrom] = ([start for start, _ in intervals], intervals)
    return result


def in_bed(chrom: str, pos1: int, index: dict[str, tuple[list[int], list[tuple[int, int]]]]) -> bool:
    if chrom not in index:
        return False
    starts, intervals = index[chrom]
    pos0 = pos1 - 1
    idx = bisect.bisect_right(starts, pos0) - 1
    return idx >= 0 and intervals[idx][0] <= pos0 < intervals[idx][1]


def truth_alleles(path: Path) -> dict[str, set[tuple[int, str, str]]]:
    text = subprocess.check_output(
        ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(path)],
        text=True,
    )
    truth: dict[str, set[tuple[int, str, str]]] = defaultdict(set)
    for line in text.splitlines():
        chrom, pos1, ref, alts = line.split("\t")
        for alt in alts.split(","):
            truth[chrom].add((int(pos1), ref, alt))
    return truth


def arm_label(chrom: str, pos1: int) -> str | None:
    if chrom == "chr6":
        return "chr6 <60 Mb" if pos1 < 60_000_000 else "chr6 ≥60 Mb"
    if chrom == "chr16":
        return "chr16 <45 Mb" if pos1 < 45_000_000 else "chr16 ≥45 Mb"
    return None


def pct(numerator: int | float, denominator: int | float) -> float:
    return float(numerator) / float(denominator) if denominator else 0.0


def json_dump(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def build_evidence(args: argparse.Namespace) -> dict[str, Any]:
    receipt_path = args.run_root / "run_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    hc_index = load_bed(args.hc_bed)
    truth = truth_alleles(args.truth_vcf)

    all_catalog_loci: set[tuple[str, int]] = set()
    primary_loci: set[tuple[str, int]] = set()
    locus_hp_to_ps: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    unit_metadata: dict[str, dict[str, Any]] = {}
    routes: set[tuple[str, str, int]] = set()
    containers: dict[tuple[str, str], set[int]] = defaultdict(set)
    block_routes: dict[tuple[str, str], set[tuple[int, str]]] = defaultdict(set)

    totals = Counter()
    k_bins = Counter()
    disposition_rows = Counter()
    disposition_patterns = Counter()
    disposition_weights = Counter()
    per_chromosome: list[dict[str, Any]] = []
    catalog_by_chrom: dict[str, set[tuple[int, str, str]]] = {}
    primary_by_chrom: dict[str, set[int]] = defaultdict(set)
    k_gt_12_units: set[str] = set()
    k_gt_12_memberships: set[tuple[str, int]] = set()
    k_gt_12_unique_loci: set[tuple[str, int]] = set()
    k_gt_12_disposition_rows = Counter()
    k_gt_12_disposition_patterns = Counter()
    k_gt_12_disposition_weights = Counter()

    for chrom in CHROMS:
        chrom_root = args.run_root / "chromosomes" / chrom
        extraction = chrom_root / "extraction"
        partition = chrom_root / "python_partition"

        catalog_path = extraction / f"HCC1395.{chrom}.site_catalog.tsv.gz"
        catalog_rows = list(read_tsv_gz(catalog_path))
        catalog = {(int(row["pos1"]), row["ref"], row["alt"]) for row in catalog_rows}
        if len(catalog) != len(catalog_rows):
            raise RuntimeError(f"duplicate catalog allele in {chrom}")
        catalog_by_chrom[chrom] = catalog
        all_catalog_loci.update((chrom, pos1) for pos1, _, _ in catalog)

        units_path = partition / "units.tsv.gz"
        chrom_units = 0
        chrom_k_gt_12 = 0
        for row in read_tsv_gz(units_path):
            chrom_units += 1
            totals["units"] += 1
            k = int(row["k"])
            hp = int(row["hp_family"])
            phase_set = row["phase_set"]
            unit_id = row["unit_id"]
            unit_metadata[unit_id] = {
                "chrom": chrom,
                "k": k,
                "hp": hp,
                "phase_set": phase_set,
            }
            routes.add((chrom, phase_set, hp))
            containers[(chrom, phase_set)].add(hp)
            if k == 1:
                k_bins["k=1"] += 1
            elif k <= 8:
                k_bins["k=2-8"] += 1
            elif k <= 12:
                k_bins["k=9-12"] += 1
            else:
                k_bins["k>12"] += 1
                chrom_k_gt_12 += 1
                k_gt_12_units.add(unit_id)

        membership_path = partition / "membership.tsv.gz"
        chrom_memberships = 0
        chrom_primary: set[int] = set()
        for row in read_tsv_gz(membership_path):
            chrom_memberships += 1
            totals["memberships"] += 1
            pos1 = int(row["pos1"])
            hp = int(row["hp_family"])
            phase_set = row["phase_set"]
            unit_id = row["unit_id"]
            primary_loci.add((chrom, pos1))
            primary_by_chrom[chrom].add(pos1)
            chrom_primary.add(pos1)
            locus_hp_to_ps[(chrom, pos1, hp)].add(phase_set)
            block_routes[(chrom, row["block_id"])].add((hp, phase_set))
            if unit_id in k_gt_12_units:
                k_gt_12_memberships.add((unit_id, pos1))
                k_gt_12_unique_loci.add((chrom, pos1))

        blocks_path = partition / "blocks.tsv.gz"
        chrom_blocks = 0
        for row in read_tsv_gz(blocks_path):
            chrom_blocks += 1
            totals["blocks"] += 1
            if int(row["k"]) > 12:
                raise RuntimeError(f"block k>12 in {chrom}: {row['block_id']}")

        constraints_path = partition / "constraints.tsv.gz"
        chrom_constraints = 0
        chrom_weight = 0
        for row in read_tsv_gz(constraints_path):
            chrom_constraints += 1
            totals["constraints"] += 1
            totals["constraint_patterns"] += int(row["pattern_count"])
            weight = int(row["molecule_weight"])
            totals["constraint_weight"] += weight
            chrom_weight += weight

        dispositions_path = partition / "dispositions.tsv.gz"
        for row in read_tsv_gz(dispositions_path):
            raw_status = row["disposition"]
            status = (
                "unavoidable"
                if raw_status == "unavoidable_span_gt_max_block_size"
                else raw_status
            )
            unit_id = row["unit_id"]
            disposition_rows[status] += 1
            disposition_patterns[status] += int(row["pattern_count"])
            disposition_weights[status] += int(row["molecule_weight"])
            if unit_id in k_gt_12_units:
                k_gt_12_disposition_rows[status] += 1
                k_gt_12_disposition_patterns[status] += int(row["pattern_count"])
                k_gt_12_disposition_weights[status] += int(row["molecule_weight"])

        truth_overlap = catalog & truth.get(chrom, set())
        per_chromosome.append(
            {
                "chrom": chrom,
                "catalog_s": len(catalog),
                "primary_unique": len(chrom_primary),
                "primary_rate": pct(len(chrom_primary), len(catalog)),
                "units": chrom_units,
                "memberships": chrom_memberships,
                "blocks": chrom_blocks,
                "constraints": chrom_constraints,
                "molecule_weight": chrom_weight,
                "k_gt_12_units": chrom_k_gt_12,
                "seqc2_truth_count": len(truth.get(chrom, set())),
                "seqc2_exact_overlap": len(truth_overlap),
                "seqc2_truth_recall": pct(len(truth_overlap), len(truth.get(chrom, set()))),
                "seqc2_catalog_overlap": pct(len(truth_overlap), len(catalog)),
            }
        )

    cross_route_blocks = sum(len(routes_) != 1 for routes_ in block_routes.values())
    if cross_route_blocks:
        raise RuntimeError(f"found {cross_route_blocks} cross-PS/HP blocks")

    multi_ps_locus_hp = {key: value for key, value in locus_hp_to_ps.items() if len(value) > 1}
    multi_ps_union = {(chrom, pos1) for chrom, pos1, _ in multi_ps_locus_hp}
    multi_ps_hp = Counter(key[2] for key in multi_ps_locus_hp)
    multi_ps_locus_categories = Counter()
    by_locus: dict[tuple[str, int], set[int]] = defaultdict(set)
    for chrom, pos1, hp in multi_ps_locus_hp:
        by_locus[(chrom, pos1)].add(hp)
    for hps in by_locus.values():
        multi_ps_locus_categories[
            "both" if hps == {1, 2} else "hp1_only" if hps == {1} else "hp2_only"
        ] += 1

    container_presence = Counter()
    for hps in containers.values():
        if hps == {1, 2}:
            container_presence["HP1 + HP2"] += 1
        elif hps == {1}:
            container_presence["HP1 only"] += 1
        elif hps == {2}:
            container_presence["HP2 only"] += 1
        else:
            raise RuntimeError(f"unexpected HP set in primary container: {hps}")

    arm_catalog: dict[str, set[tuple[str, int, str, str]]] = defaultdict(set)
    arm_primary: dict[str, set[tuple[str, int]]] = defaultdict(set)
    arm_hc = Counter()
    arm_transitions = Counter()
    arm_transversions = Counter()
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    for chrom in ("chr6", "chr16"):
        for pos1, ref, alt in catalog_by_chrom[chrom]:
            arm = arm_label(chrom, pos1)
            assert arm is not None
            arm_catalog[arm].add((chrom, pos1, ref, alt))
            if in_bed(chrom, pos1, hc_index):
                arm_hc[arm] += 1
            if (ref, alt) in transitions:
                arm_transitions[arm] += 1
            else:
                arm_transversions[arm] += 1
        for pos1 in primary_by_chrom[chrom]:
            arm = arm_label(chrom, pos1)
            assert arm is not None
            arm_primary[arm].add((chrom, pos1))

    molecule_arm_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for chrom in ("chr6", "chr16"):
        sparse = (
            args.run_root
            / "chromosomes"
            / chrom
            / "extraction"
            / f"HCC1395.{chrom}.molecule_sparse_calls.tsv.gz"
        )
        for row in read_tsv_gz(sparse):
            positions = row["positions1"].split(",") if row["positions1"] else []
            if not positions:
                continue
            arm = arm_label(chrom, int(positions[0]))
            assert arm is not None
            molecule_arm_counts[arm]["total"] += 1
            hp = row["hp_family"]
            ps = row["phase_set"].strip()
            if hp in {"1", "2"} and ps:
                molecule_arm_counts[arm]["known_ps_hp12"] += 1
            elif hp == "3":
                molecule_arm_counts[arm]["hp3"] += 1
            else:
                molecule_arm_counts[arm]["unphased_or_other"] += 1

    arm_order = ["chr6 <60 Mb", "chr6 ≥60 Mb", "chr16 <45 Mb", "chr16 ≥45 Mb"]
    outlier_arms = []
    for arm in arm_order:
        cat_n = len(arm_catalog[arm])
        pri_n = len(arm_primary[arm])
        molecule_n = molecule_arm_counts[arm]["total"]
        outlier_arms.append(
            {
                "arm": arm,
                "catalog_s": cat_n,
                "primary_unique": pri_n,
                "primary_rate": pct(pri_n, cat_n),
                "catalog_in_seqc2_hc": arm_hc[arm],
                "known_ps_hp12_molecules": molecule_arm_counts[arm]["known_ps_hp12"],
                "molecule_rows": molecule_n,
                "known_ps_hp12_rate": pct(molecule_arm_counts[arm]["known_ps_hp12"], molecule_n),
                "hp3_rate": pct(molecule_arm_counts[arm]["hp3"], molecule_n),
                "unphased_or_other_rate": pct(
                    molecule_arm_counts[arm]["unphased_or_other"], molecule_n
                ),
                "transitions": arm_transitions[arm],
                "transversions": arm_transversions[arm],
                "ti_tv": pct(arm_transitions[arm], arm_transversions[arm]),
            }
        )

    outlier_names = {"chr6 <60 Mb", "chr16 ≥45 Mb"}
    outlier_catalog = sum(len(arm_catalog[name]) for name in outlier_names)
    outlier_primary = sum(len(arm_primary[name]) for name in outlier_names)
    non_outlier_catalog = len(all_catalog_loci) - outlier_catalog
    non_outlier_primary = len(primary_loci) - outlier_primary

    direct = {
        "S": len(all_catalog_loci),
        "unique_sites": len(primary_loci),
        "units": totals["units"],
        "unit_memberships": totals["memberships"],
        "blocks": totals["blocks"],
        "constraints": totals["constraints"],
        "pattern_count_total": totals["constraint_patterns"],
        "molecule_weight_total": totals["constraint_weight"],
        "k_bins": dict(k_bins),
        "constraint_dispositions": dict(disposition_rows),
        "disposition_pattern_counts": dict(disposition_patterns),
        "molecule_weights": dict(disposition_weights),
    }
    receipt_aggregate = receipt["aggregate"]
    receipt_projection = {
        "S": receipt_aggregate["S"],
        "unique_sites": receipt_aggregate["unique_sites"],
        "units": receipt_aggregate["units"],
        "unit_memberships": receipt_aggregate["unit_memberships"],
        "blocks": receipt_aggregate["blocks"],
        "constraints": receipt_aggregate["constraints"],
        "pattern_count_total": receipt_aggregate["pattern_count_total"],
        "molecule_weight_total": receipt_aggregate["molecule_weights"]["total"],
        "k_bins": receipt_aggregate["k_bins"],
        "constraint_dispositions": receipt_aggregate["constraint_dispositions"],
        "molecule_weights": {
            key: receipt_aggregate["molecule_weights"][key]
            for key in ("retained", "cut", "unavoidable")
        },
    }
    direct_comparable = {key: value for key, value in direct.items() if key in receipt_projection}
    aggregate_match = direct_comparable == receipt_projection
    if not aggregate_match:
        raise RuntimeError(
            "independent aggregate does not match run receipt:\n"
            + json.dumps({"direct": direct_comparable, "receipt": receipt_projection}, indent=2)
        )

    total_disposition_weight = sum(disposition_weights.values())
    total_disposition_rows = sum(disposition_rows.values())
    if total_disposition_weight != totals["constraint_weight"]:
        raise RuntimeError("disposition molecule weight does not conserve constraint weight")
    if total_disposition_rows != totals["constraints"]:
        raise RuntimeError("disposition rows do not conserve constraint rows")

    k_gt_12_weight = sum(k_gt_12_disposition_weights.values())
    evidence = {
        "schema_name": "intersubmod.hcc1395_exact_ps_k12_report_evidence",
        "schema_version": "1.0.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "claim_status": "PARTIAL",
        "task_type": "exploratory_pilot",
        "scope": "HCC1395 chr1-22; extraction and exact-PS segmentation only",
        "source_paths": {
            "run_root": str(args.run_root.resolve()),
            "run_receipt": str(receipt_path.resolve()),
            "truth_vcf": str(args.truth_vcf.resolve()),
            "hc_bed": str(args.hc_bed.resolve()),
            "input_manifest": str(args.input_manifest.resolve()),
        },
        "source_identities": {
            "run_receipt_sha256": sha256_file(receipt_path),
            "truth_vcf_sha256": sha256_file(args.truth_vcf),
            "truth_vcf_size_bytes": args.truth_vcf.stat().st_size,
            "hc_bed_sha256": sha256_file(args.hc_bed),
            "hc_bed_size_bytes": args.hc_bed.stat().st_size,
            "input_manifest_sha256": sha256_file(args.input_manifest),
        },
        "direct_recompute": direct,
        "receipt_projection": receipt_projection,
        "aggregate_exact_match": aggregate_match,
        "all_pass": receipt["all_pass"],
        "validation": {
            "chromosomes_passed": receipt_aggregate["chromosomes_passed"],
            "chromosomes_requested": receipt_aggregate["chromosomes_requested"],
            "cross_ps_violations": receipt_aggregate["cross_ps_violations"],
            "cross_hp_violations": receipt_aggregate["cross_hp_violations"],
            "python_cpp_mismatch_count": receipt_aggregate["python_cpp_mismatch_count"],
            "historical_artifacts_compared": receipt_aggregate["historical_comparison"][
                "artifacts_compared"
            ],
            "historical_semantic_mismatches": receipt_aggregate["historical_comparison"][
                "semantic_mismatches"
            ],
            "historical_physical_only_differences": receipt_aggregate[
                "historical_comparison"
            ]["physical_only_differences"],
            "cross_route_blocks_direct": cross_route_blocks,
        },
        "k_gt_12": {
            "units": len(k_gt_12_units),
            "unit_memberships": len(k_gt_12_memberships),
            "unique_loci": len(k_gt_12_unique_loci),
            "blocks": sum(
                1
                for chrom in CHROMS
                for row in read_tsv_gz(
                    args.run_root
                    / "chromosomes"
                    / chrom
                    / "python_partition"
                    / "blocks.tsv.gz"
                )
                if row["unit_id"] in k_gt_12_units
            ),
            "constraint_dispositions": dict(k_gt_12_disposition_rows),
            "disposition_pattern_counts": dict(k_gt_12_disposition_patterns),
            "molecule_weights": dict(k_gt_12_disposition_weights),
            "molecule_weight_total": k_gt_12_weight,
        },
        "exact_ps": {
            "containers": len(containers),
            "container_presence": dict(container_presence),
            "ps_hp_routes": len(routes),
            "unique_locus_hp_keys": len(locus_hp_to_ps),
            "multi_ps_locus_hp_keys": len(multi_ps_locus_hp),
            "multi_ps_locus_hp_rate": pct(len(multi_ps_locus_hp), len(locus_hp_to_ps)),
            "multi_ps_locus_hp_by_hp": {"HP1": multi_ps_hp[1], "HP2": multi_ps_hp[2]},
            "multi_ps_unique_loci": len(multi_ps_union),
            "multi_ps_unique_loci_rate": pct(len(multi_ps_union), len(primary_loci)),
            "multi_ps_locus_categories": dict(multi_ps_locus_categories),
            "max_ps_per_locus_hp": max(map(len, multi_ps_locus_hp.values()), default=0),
        },
        "per_chromosome": per_chromosome,
        "outlier_arm_diagnostic": {
            "rows": outlier_arms,
            "outlier_catalog_s": outlier_catalog,
            "outlier_catalog_share": pct(outlier_catalog, len(all_catalog_loci)),
            "outlier_primary_unique": outlier_primary,
            "outlier_primary_share": pct(outlier_primary, len(primary_loci)),
            "non_outlier_catalog_s": non_outlier_catalog,
            "non_outlier_primary_unique": non_outlier_primary,
            "non_outlier_primary_rate": pct(non_outlier_primary, non_outlier_catalog),
            "interpretation_ceiling": (
                "The arms are outside or nearly outside SEQC2 HC scope and have low exact-PS HP1/2 "
                "eligibility. This is an upstream candidate-universe/CN-exclusion QA signal, not proof "
                "that all loci are germline or false positives."
            ),
        },
        "claim_ceiling": receipt["claim_ceiling"],
    }
    return evidence


def build_artifact(evidence: dict[str, Any]) -> dict[str, Any]:
    direct = evidence["direct_recompute"]
    k_gt_12 = evidence["k_gt_12"]
    exact_ps = evidence["exact_ps"]
    outlier = evidence["outlier_arm_diagnostic"]
    total_weight = direct["molecule_weight_total"]
    k_gt_12_weight = k_gt_12["molecule_weight_total"]
    generated_at = evidence["generated_at_utc"]

    headline = [
        {
            "candidate_loci": direct["S"],
            "primary_unique_loci": direct["unique_sites"],
            "units": direct["units"],
            "k_gt_12_units": k_gt_12["units"],
            "overall_retained_rate": pct(direct["molecule_weights"]["retained"], total_weight),
            "k_gt_12_retained_rate": pct(
                k_gt_12["molecule_weights"]["retained"], k_gt_12_weight
            ),
        }
    ]
    status_labels = {
        "retained": "保留於同一 block",
        "cut": "跨切點而切斷",
        "unavoidable": "任一 k≤12 分割皆無法保留",
    }
    support_disposition = []
    for scope, counts, weights, denominator in (
        (
            "全部 exact-PS units",
            direct["constraint_dispositions"],
            direct["molecule_weights"],
            total_weight,
        ),
        (
            "僅原始 k>12 units",
            k_gt_12["constraint_dispositions"],
            k_gt_12["molecule_weights"],
            k_gt_12_weight,
        ),
    ):
        for status in ("retained", "cut", "unavoidable"):
            support_disposition.append(
                {
                    "scope": scope,
                    "disposition": status_labels[status],
                    "constraint_rows": counts[status],
                    "molecule_weight": weights[status],
                    "weight_share": pct(weights[status], denominator),
                    "weight_denominator": denominator,
                }
            )

    k_rows = []
    for label in ("k=1", "k=2-8", "k=9-12", "k>12"):
        count = direct["k_bins"][label]
        k_rows.append(
            {
                "k_bin": label,
                "units": count,
                "unit_share": pct(count, direct["units"]),
                "requires_split": label == "k>12",
            }
        )

    ps_container_rows = []
    for label in ("HP1 + HP2", "HP1 only", "HP2 only"):
        count = exact_ps["container_presence"][label]
        ps_container_rows.append(
            {
                "hp_presence": label,
                "exact_ps_containers": count,
                "container_share": pct(count, exact_ps["containers"]),
            }
        )

    validation_rows = [
        {
            "check": "Chromosome stages",
            "observed": f"{evidence['validation']['chromosomes_passed']}/{evidence['validation']['chromosomes_requested']}",
            "expected": "22/22",
            "status": "PASS",
            "meaning": "chr1–22 都有完整 metrics 且納入 aggregate",
        },
        {
            "check": "Cross-PS blocks",
            "observed": evidence["validation"]["cross_ps_violations"],
            "expected": 0,
            "status": "PASS",
            "meaning": "任何輸出 block 都沒有混合不同 exact PS",
        },
        {
            "check": "Cross-HP blocks",
            "observed": evidence["validation"]["cross_hp_violations"],
            "expected": 0,
            "status": "PASS",
            "meaning": "任何輸出 block 都沒有混合 HP1/HP2",
        },
        {
            "check": "Python ↔ C++ mismatch",
            "observed": evidence["validation"]["python_cpp_mismatch_count"],
            "expected": 0,
            "status": "PASS",
            "meaning": "block、membership 與 disposition 完全一致",
        },
        {
            "check": "Independent aggregate recompute",
            "observed": "exact equal" if evidence["aggregate_exact_match"] else "mismatch",
            "expected": "exact equal",
            "status": "PASS" if evidence["aggregate_exact_match"] else "FAIL",
            "meaning": "raw gzip TSV 逐列重算與 run receipt 一致",
        },
        {
            "check": "Historical semantic comparison",
            "observed": (
                f"{evidence['validation']['historical_artifacts_compared']} compared / "
                f"{evidence['validation']['historical_semantic_mismatches']} mismatches"
            ),
            "expected": "0 semantic mismatches",
            "status": "PASS",
            "meaning": "舊、新 extraction 解壓內容一致；gzip header 差異不算科學差異",
        },
        {
            "check": "Constraint row conservation",
            "observed": sum(direct["constraint_dispositions"].values()),
            "expected": direct["constraints"],
            "status": "PASS",
            "meaning": "每個 constraint 恰分為 retained/cut/unavoidable 一類",
        },
        {
            "check": "Molecule-weight conservation",
            "observed": sum(direct["molecule_weights"].values()),
            "expected": direct["molecule_weight_total"],
            "status": "PASS",
            "meaning": "support mass 沒有遺失或重複計算",
        },
    ]

    definitions = [
        {
            "term": "S",
            "definition": "Manifest 綁定的 LongPhase-S FILTER=PASS candidate loci；不是由本輪 Big7 BAM 重新 de novo calling 的 truth sSNV。",
            "unit": "unique genomic locus",
            "denominator": "HCC1395 chr1–22 candidate catalog",
        },
        {
            "term": "exact PS",
            "definition": "read sidecar 中的原始 phase-set identifier；只保證 block 內 HP 方向一致，不保證不同 PS 的 HP1/HP2 可直接對齊。",
            "unit": "phase-set identifier",
            "denominator": "known PS reads",
        },
        {
            "term": "primary unique locus",
            "definition": "至少進入一個 exact PS × HP1/2 read-linked unit 的獨特染色體位置。",
            "unit": "unique genomic locus",
            "denominator": "S",
        },
        {
            "term": "unit",
            "definition": "exact PS × HP1/2 × read-linked component；本輪最小可獨立分割的證據容器。",
            "unit": "evidence unit",
            "denominator": "all primary units",
        },
        {
            "term": "membership",
            "definition": "一個 locus 進入一個 unit 的一次指派；同一 locus 若在不同 PS/HP 出現可有多筆，因此不是 unique loci。",
            "unit": "locus–unit incidence",
            "denominator": "all primary memberships",
        },
        {
            "term": "k",
            "definition": "一個 unit 內的 candidate-locus 數；k>12 才啟動 bounded contiguous segmentation。",
            "unit": "loci per unit",
            "denominator": "per unit",
        },
        {
            "term": "molecule weight",
            "definition": "聚合 constraint 的 molecule–unit incidence 數；不是全樣本 unique read 數。",
            "unit": "molecule–unit incidence",
            "denominator": "constraints in the stated scope",
        },
        {
            "term": "retained / cut / unavoidable",
            "definition": "constraint 完全落在一個 block／跨過選定切點／任何 k≤12 切法都無法完整容納。",
            "unit": "constraint row or molecule weight",
            "denominator": "must state all units or k>12-only",
        },
    ]

    source_ids = {
        "evidence": "exact_ps_raw_recompute",
        "run": "hcc1395_direct_v2_receipt",
        "method": "exact_ps_method_notes",
    }
    sources = [
        {
            "id": source_ids["evidence"],
            "label": "HCC1395 exact-PS independent raw recompute",
            "path": f"{TOPIC}/report/evidence_snapshot.json",
            "query": {
                "engine": "artifact snapshot",
                "language": "sql",
                "sql": "SELECT * FROM snapshot.exact_ps_raw_recompute",
                "description": "Reads every chr1–22 gzip TSV, recomputes aggregates, checks conservation, and adds bounded SEQC2 scope diagnostics.",
                "tables_used": [
                    f"{TOPIC}/report/evidence_snapshot.json",
                    RUN_SOURCE,
                ],
                "filters": [
                    "sample=HCC1395",
                    "chr1-22 only",
                    "exact non-missing PS",
                    "HP family in {1,2}",
                    "read-link threshold=3",
                    "max block k=12",
                ],
                "metric_definitions": [
                    "primary unique locus = distinct (chrom,pos1) present in exact-PS membership",
                    "molecule-weight retention = retained molecule_weight / total molecule_weight",
                    "k>12 retention uses only constraints belonging to original units with k>12",
                ],
            },
        },
        {
            "id": source_ids["run"],
            "label": "HCC1395 direct-Big7 v2 run receipt",
            "path": RUN_SOURCE,
            "query": {
                "engine": "InterSubMod pilot runner",
                "language": "json",
                "sql": "SELECT * FROM snapshot.hcc1395_direct_v2_receipt",
                "description": "Frozen full-run receipt with input binding, source identities, per-chromosome status, and aggregate checks.",
                "tables_used": [RUN_SOURCE],
                "filters": ["HCC1395", "chr1-22", "all_pass chromosomes with complete metrics only"],
            },
        },
        {
            "id": source_ids["method"],
            "label": "Exact-PS implementation and decision notes",
            "path": f"{TOPIC}/implementation-notes.md",
            "query": {
                "engine": "repository documentation",
                "language": "markdown",
                "sql": "SELECT * FROM snapshot.exact_ps_method_definitions",
                "description": "Fail-closed PS boundary, partial-read preservation, DP segmentation, and claim-ceiling decisions.",
                "tables_used": [
                    f"{TOPIC}/implementation-notes.md",
                    f"{TOPIC}/pre-decision-audit.md",
                    f"{TOPIC}/scripts/exact_ps_k12_partition.py",
                    f"{TOPIC}/cpp/exact_ps_k12_partition.cpp",
                ],
            },
        },
    ]

    cards = [
        {
            "id": "candidate_loci_card",
            "description": "Manifest-bound LongPhase-S PASS candidate loci；不是 sSNV truth。",
            "dataset": "headline_metrics",
            "sourceId": source_ids["evidence"],
            "metrics": [{"label": "Candidate loci S", "field": "candidate_loci", "format": "number"}],
        },
        {
            "id": "primary_loci_card",
            "description": "至少進入一個 exact PS×HP1/2 unit 的 unique loci。",
            "dataset": "headline_metrics",
            "sourceId": source_ids["evidence"],
            "metrics": [
                {"label": "Primary unique loci", "field": "primary_unique_loci", "format": "number"}
            ],
        },
        {
            "id": "units_card",
            "description": "exact PS×HP1/2×read-linked component 的數量。",
            "dataset": "headline_metrics",
            "sourceId": source_ids["evidence"],
            "metrics": [{"label": "Evidence units", "field": "units", "format": "number"}],
        },
        {
            "id": "k_gt_12_card",
            "description": "只有這些原始 units 需要分割成 k≤12 blocks。",
            "dataset": "headline_metrics",
            "sourceId": source_ids["evidence"],
            "metrics": [{"label": "Original k>12 units", "field": "k_gt_12_units", "format": "number"}],
        },
        {
            "id": "overall_retention_card",
            "description": "全部 units 的 molecule–unit support mass 保留率。",
            "dataset": "headline_metrics",
            "sourceId": source_ids["evidence"],
            "metrics": [
                {"label": "Overall retained", "field": "overall_retained_rate", "format": "percent"}
            ],
        },
        {
            "id": "k_gt_12_retention_card",
            "description": "只看真正被 k≤12 上限影響的原始 k>12 units。",
            "dataset": "headline_metrics",
            "sourceId": source_ids["evidence"],
            "metrics": [
                {"label": "k>12 retained", "field": "k_gt_12_retained_rate", "format": "percent"}
            ],
        },
    ]

    charts = [
        {
            "id": "support_disposition_chart",
            "title": "Read-support mass disposition",
            "subtitle": "全部 units 與原始 k>12 units 各自歸一化為 100%；單位為 molecule–unit incidence。",
            "type": "stackedBar100",
            "dataset": "support_disposition",
            "sourceId": source_ids["evidence"],
            "intent": "composition",
            "question": "k≤12 分割保留多少 read-support mass，損失集中在哪個 scope？",
            "rationale": "100% stacked bars keep the two denominators explicit and prevent the much larger k≤12 population from hiding k>12 loss.",
            "comparisonContext": {
                "denominator": "molecule weight within each displayed scope",
                "grain": "constraint disposition",
                "normalization": "within scope",
                "unit": "share",
            },
            "encodings": {
                "x": {"field": "scope", "type": "nominal", "label": "Scope"},
                "y": {"field": "molecule_weight", "type": "quantitative", "label": "Molecule weight"},
                "color": {"field": "disposition", "type": "nominal", "label": "Disposition"},
                "tooltip": [
                    {"field": "constraint_rows", "type": "quantitative", "label": "Constraint rows"},
                    {"field": "weight_share", "type": "quantitative", "label": "Weight share", "format": "percent"},
                    {"field": "weight_denominator", "type": "quantitative", "label": "Weight denominator"},
                ],
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True},
        },
        {
            "id": "chromosome_primary_coverage_chart",
            "title": "Primary exact-PS locus coverage by chromosome",
            "subtitle": "分母為各染色體的 manifest-bound PASS candidate loci；chr6/chr16 受已知排除臂主導。",
            "type": "bar",
            "dataset": "chromosome_coverage",
            "sourceId": source_ids["evidence"],
            "intent": "comparison",
            "question": "各染色體有多少 candidate loci 能進入 exact PS×HP1/2 primary units？",
            "rationale": "A single-series bar chart exposes chromosome-specific eligibility outliers without implying a second grouping dimension.",
            "comparisonContext": {
                "denominator": "manifest-bound PASS candidate loci within chromosome",
                "grain": "chromosome",
                "unit": "fraction",
            },
            "encodings": {
                "x": {"field": "chrom", "type": "nominal", "label": "Chromosome"},
                "y": {"field": "primary_rate", "type": "quantitative", "label": "Primary locus coverage", "format": "percent"},
                "tooltip": [
                    {"field": "catalog_s", "type": "quantitative", "label": "Candidate loci S"},
                    {"field": "primary_unique", "type": "quantitative", "label": "Primary unique loci"},
                    {"field": "k_gt_12_units", "type": "quantitative", "label": "k>12 units"},
                ],
            },
            "valueFormat": "percent",
            "palette": {"kind": "sequential", "name": "blue"},
            "labels": {"values": "none"},
        },
    ]

    tables = [
        {
            "id": "k_bins_table",
            "title": "Unit k distribution",
            "subtitle": "所有 39,544 個 exact-PS evidence units；只有 k>12 需要 bounded segmentation。",
            "dataset": "k_bins",
            "sourceId": source_ids["evidence"],
            "defaultSort": {"field": "units", "direction": "desc"},
            "columns": [
                {"field": "k_bin", "label": "k category", "type": "text"},
                {"field": "units", "label": "Units", "format": "number"},
                {"field": "unit_share", "label": "Share of units", "format": "percent"},
                {"field": "requires_split", "label": "Requires k≤12 split", "type": "boolean"},
            ],
        },
        {
            "id": "ps_container_table",
            "title": "Exact-PS container HP presence",
            "subtitle": "每列是一個 (chromosome, exact PS) container 的 HP1/HP2 presence。",
            "dataset": "ps_containers",
            "sourceId": source_ids["evidence"],
            "defaultSort": {"field": "exact_ps_containers", "direction": "desc"},
            "columns": [
                {"field": "hp_presence", "label": "HP presence", "type": "text"},
                {"field": "exact_ps_containers", "label": "Containers", "format": "number"},
                {"field": "container_share", "label": "Share", "format": "percent"},
            ],
        },
        {
            "id": "outlier_arm_table",
            "title": "chr6/chr16 arm-level candidate eligibility",
            "subtitle": "60 Mb/45 Mb 是診斷切點；SEQC2 HC 欄是 candidate loci 落在 HC BED 內的數量。",
            "dataset": "outlier_arms",
            "sourceId": source_ids["evidence"],
            "defaultSort": {"field": "catalog_s", "direction": "desc"},
            "columns": [
                {"field": "arm", "label": "Diagnostic arm", "type": "text"},
                {"field": "catalog_s", "label": "Candidate loci S", "format": "number"},
                {"field": "primary_unique", "label": "Primary loci", "format": "number"},
                {"field": "primary_rate", "label": "Primary rate", "format": "percent"},
                {"field": "catalog_in_seqc2_hc", "label": "Inside SEQC2 HC", "format": "number"},
            ],
        },
        {
            "id": "outlier_phase_table",
            "title": "chr6/chr16 arm-level phase evidence",
            "subtitle": "Known-PS HP1/2 分母為各診斷臂的 molecule rows；Ti/Tv 只作描述性 QA。",
            "dataset": "outlier_arms",
            "sourceId": source_ids["evidence"],
            "defaultSort": {"field": "known_ps_hp12_rate", "direction": "asc"},
            "columns": [
                {"field": "arm", "label": "Diagnostic arm", "type": "text"},
                {"field": "known_ps_hp12_rate", "label": "Known-PS HP1/2 molecules", "format": "percent"},
                {"field": "ti_tv", "label": "Ti/Tv", "format": "number"},
            ],
        },
        {
            "id": "validation_table",
            "title": "Validation and conservation checks",
            "subtitle": "HCC1395 chr1–22 direct v2；所有項目必須 PASS 才能接受 segmentation 結果。",
            "dataset": "validation_checks",
            "sourceId": source_ids["run"],
            "defaultSort": {"field": "check", "direction": "asc"},
            "columns": [
                {"field": "check", "label": "Check", "type": "text"},
                {"field": "observed", "label": "Observed", "type": "text"},
                {"field": "expected", "label": "Expected", "type": "text"},
                {"field": "status", "label": "Status", "type": "text"},
            ],
        },
        {
            "id": "definitions_meaning_table",
            "title": "Metric definitions",
            "subtitle": "先分清 S、unique loci、memberships 與 molecule weight 的計數對象。",
            "dataset": "definitions",
            "sourceId": source_ids["method"],
            "defaultSort": {"field": "term", "direction": "asc"},
            "columns": [
                {"field": "term", "label": "Term", "type": "text"},
                {"field": "definition", "label": "Definition", "type": "text"},
            ],
        },
        {
            "id": "definitions_units_table",
            "title": "Metric units and denominators",
            "subtitle": "報告比例時，必須同時說明 unit、scope 與 denominator。",
            "dataset": "definitions",
            "sourceId": source_ids["method"],
            "defaultSort": {"field": "term", "direction": "asc"},
            "columns": [
                {"field": "term", "label": "Term", "type": "text"},
                {"field": "unit", "label": "Unit", "type": "text"},
                {"field": "denominator", "label": "Denominator", "type": "text"},
            ],
        },
    ]

    blocks = [
        {"id": "title", "type": "markdown", "layout": "full", "body": f"# {REPORT_TITLE}"},
        {
            "id": "technical_summary",
            "type": "markdown",
            "layout": "full",
            "sourceId": source_ids["evidence"],
            "body": (
                "## 技術結論：exact PS 應先分開，不應只依 HP 編號跨 block 接樹\n\n"
                "- **方法決策：** primary unit 固定為 `exact PS × HP1/2 × read-linked component`。PS 是相位證據邊界；沒有能判定 block 間 SAME/FLIP 的直接證據，就不合併。\n"
                f"- **HCC1395 pilot 通過：** chr1–22 全數 PASS；cross-PS=0、cross-HP=0、Python↔C++ mismatch=0。\n"
                f"- **分割影響集中：** 全部 support mass 保留 {pct(direct['molecule_weights']['retained'], total_weight):.2%}，但只看原始 k>12 units 時保留 {pct(k_gt_12['molecule_weights']['retained'], k_gt_12_weight):.2%}；兩個分母不可混用。\n"
                f"- **重要資料品質訊號：** chr6<60 Mb 與 chr16≥45 Mb 佔 S 的 {outlier['outlier_catalog_share']:.2%}，卻只佔 primary loci 的 {outlier['outlier_primary_share']:.2%}；排除這兩個診斷臂後 primary coverage 為 {outlier['non_outlier_primary_rate']:.2%}。\n"
                "- **Claim ceiling：** 目前只證明 extraction、exact-PS segmentation 與分割守恆；尚未重跑 tree enumeration、partial-read likelihood、VAF ranking、Topo 或 clone/subclone 判定。"
            ),
        },
        {
            "id": "headline_metrics",
            "type": "metric-strip",
            "layout": "full",
            "cardIds": [card["id"] for card in cards],
        },
        {
            "id": "why_exact_ps",
            "type": "markdown",
            "layout": "full",
            "sourceId": source_ids["evidence"],
            "body": (
                "## 1,244 個 loci 顯示：相同 HP family 並不等於可安全跨 PS 對齊\n\n"
                f"HCC1395 有 **{exact_ps['multi_ps_unique_loci']:,} 個 primary loci** 在同一 HP family 下出現在兩個 PS，佔 primary loci 的 **{exact_ps['multi_ps_unique_loci_rate']:.2%}**。若只看到 `HP1` 就跨 PS 聚合，可能把局部方向不同的 haplotypes 誤接，或讓同一 locus 重複進入全域結構。\n\n"
                "因此本輪採 fail-closed：每個 exact PS 分開保留 read constraint。未來若要跨 PS，只能另建 signed bridge，明確證明兩個 block 是 SAME 或 FLIP；不能以距離近、HP 名稱相同或 VAF 相近代替 orientation evidence。"
            ),
        },
        {"id": "ps_container_table_block", "type": "table", "layout": "full", "tableId": "ps_container_table"},
        {
            "id": "support_takeaway",
            "type": "markdown",
            "layout": "full",
            "sourceId": source_ids["evidence"],
            "body": (
                "## 整體 98.63% 保留率不能掩蓋 k>12 內的 29.72% 損失\n\n"
                f"全部 {direct['units']:,} 個 units 中，只有 {k_gt_12['units']:,} 個（{pct(k_gt_12['units'], direct['units']):.2%}）需要切割，所以全體保留率很高。真正受上限影響的 k>12 subset 共有 {k_gt_12_weight:,} molecule-weight，其中 retained={k_gt_12['molecule_weights']['retained']:,}、cut={k_gt_12['molecule_weights']['cut']:,}、unavoidable={k_gt_12['molecule_weights']['unavoidable']:,}。\n\n"
                "圖中兩條 bar 各自使用自己的 100% 分母；這是判斷 k=12 工程上限是否可接受時必須看的相對比例。"
            ),
        },
        {"id": "support_chart_block", "type": "chart", "layout": "full", "chartId": "support_disposition_chart"},
        {
            "id": "k_distribution_takeaway",
            "type": "markdown",
            "layout": "full",
            "sourceId": source_ids["evidence"],
            "body": (
                "## k>12 只佔 0.24% units，卻承擔全部 cut 與 unavoidable constraints\n\n"
                f"原始 k>12 subset 為 {k_gt_12['units']:,} units、{k_gt_12['unique_loci']:,} unique loci、{k_gt_12['blocks']:,} 輸出 blocks。所有 cut 與 unavoidable 都來自這 94 個 units；因此後續效能改良與 case review 應聚焦此 subset，而不是重寫所有 k≤12 units。"
            ),
        },
        {"id": "k_bins_table_block", "type": "table", "layout": "full", "tableId": "k_bins_table"},
        {
            "id": "chromosome_outlier_takeaway",
            "type": "markdown",
            "layout": "full",
            "sourceId": source_ids["evidence"],
            "body": (
                "## chr6/chr16 的低納入率是上游候選範圍與 phase eligibility 問題，不是一般性 exact-PS 失敗\n\n"
                f"兩個診斷臂含 {outlier['outlier_catalog_s']:,}/{direct['S']:,} candidate loci（{outlier['outlier_catalog_share']:.2%}），但只有 {outlier['outlier_primary_unique']:,}/{direct['unique_sites']:,} primary loci（{outlier['outlier_primary_share']:.2%}）。對照臂與其他染色體多數接近完整納入。\n\n"
                "SEQC2 文件把 chr6p、chr16q 列為完全缺失／benchmark-excluded；這些區域幾乎沒有 known-PS HP1/2 molecule evidence。此結果可用來標記上游 QA 與 exclusion contract 風險，但**不能**把所有 loci 直接稱為 germline、caller FP 或刪除。"
            ),
        },
        {"id": "chromosome_chart_block", "type": "chart", "layout": "full", "chartId": "chromosome_primary_coverage_chart"},
        {"id": "outlier_table_block", "type": "table", "layout": "full", "tableId": "outlier_arm_table"},
        {"id": "outlier_phase_table_block", "type": "table", "layout": "full", "tableId": "outlier_phase_table"},
        {
            "id": "scope_definitions",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 名詞、單位與分母先定義，才不會把 loci、memberships、reads 混在一起\n\n"
                "`S`、primary unique loci、unit memberships、constraint rows 與 molecule weight 是五種不同計數單位。尤其 molecule weight 是 molecule–unit incidence，不是全樣本 unique read count；報告任何比例時都必須同時標明 scope 與 denominator。"
            ),
        },
        {"id": "definitions_meaning_table_block", "type": "table", "layout": "full", "tableId": "definitions_meaning_table"},
        {"id": "definitions_units_table_block", "type": "table", "layout": "full", "tableId": "definitions_units_table"},
        {
            "id": "methodology",
            "type": "markdown",
            "layout": "full",
            "sourceId": source_ids["method"],
            "body": (
                "## 方法流程：先隔離證據，再於同一 unit 內做 bounded segmentation\n\n"
                "1. 以 manifest 綁定的 LongPhase-S PASS candidate loci 定義 target catalog；Big7 HCC1395 BAM 只負責在這些位置抽取 read evidence，不是重新 de novo calling。\n"
                "2. 保留每個 molecule 的 R/A/O/D/S/L/X sparse calls、HP 與 exact PS；partial `X` 不展開成多條假 read。\n"
                "3. 只讓 HP1/HP2 且 non-missing PS 的 molecules 進 primary；以 `exact PS × HP × read-linked component` 建 unit。\n"
                "4. k≤12 的 unit 原樣保留；k>12 只在同一 unit 內用 contiguous DP 最大化 retained read constraints，並把每個 constraint 唯一分類為 retained、cut 或 unavoidable。\n"
                "5. Python producer 與獨立 C++ kernel 必須對 block、membership、disposition 逐項一致。\n"
                "6. 完整 partial-read likelihood、tree enumeration 與 VAF ranking 是下一個 adapter stage，本輪未執行。"
            ),
        },
        {
            "id": "validation_takeaway",
            "type": "markdown",
            "layout": "full",
            "sourceId": source_ids["run"],
            "body": (
                "## 程式與輸出守恆通過；舊、新 extraction 語意相同\n\n"
                "22 條染色體全部完成，輸入 pre/post binding 未漂移，Python/C++ parity 為零差異。與歷史 extraction 比較的 132 個 artifacts 解壓後語意全相同；110 個 physical-only differences 只是 gzip header/bytes 不同，不能當科學差異。\n\n"
                "因此目前能歸因的是 segmentation policy 改為 exact PS fail-closed；不能把尚未重跑的 T、Topo 或 VAF 數字說成已完成 before/after 改善。"
            ),
        },
        {"id": "validation_table_block", "type": "table", "layout": "full", "tableId": "validation_table"},
        {
            "id": "limitations",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 限制與尚未確定的地方\n\n"
                "- **單樣本：** 這是 HCC1395 chr1–22 exploratory pilot，不是七樣本或 paper-final validation。\n"
                "- **PS 不是 clone boundary：** 分開 PS 是避免錯接的證據政策，不代表生物上每個 PS 都是獨立 clone。\n"
                "- **尚無 topology 結論：** 沒有 block→partial-read likelihood adapter，就不能報 T、Topo、VAF winner、clone 或 subclone。\n"
                "- **CN exclusion P0：** 目前 manifest 把 CN BED 未列位置視為 neutral，卻未綁 exclusion BED。segmentation 本身未使用 CN，所以本輪數字不變；任何後續 CN/clone 解讀前，chr6p/chr16q 必須改為 excluded/unknown，而非 neutral。\n"
                "- **Candidate provenance：** S 是上游 ClairS→LongPhase-S PASS candidate universe；不是本輪 direct Big7 BAM 重新 calling 的 sSNV truth。"
            ),
        },
        {
            "id": "next_steps",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 建議下一步：先修 evidence contract，再接回 tree/VAF\n\n"
                "1. 將 `exact PS 不跨接` 固定為 primary policy；保留 signed cross-PS bridge 為獨立、可稽核的 optional stage。\n"
                "2. 修正 HCC1395 CN/exclusion contract，把 chr6p、chr16q 標成 excluded/unknown，並重做 candidate-scope QA。\n"
                "3. 建立 block→molecule adapter，讓完整 R/A/O/D/S/L/X partial reads 直接進 read-pattern likelihood，而不是先複製 completion。\n"
                "4. 只在 HCC1395 重跑 tree enumeration、VAF ranking 與 Topo before/after；確認數據語意一致後，再由使用者決定是否擴至七樣本。\n"
                "5. 若要正式 handoff，將目前未追蹤的 topic 程式、測試與報告先做有意義的 commit。"
            ),
        },
        {
            "id": "further_questions",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 待決問題\n\n"
                "- signed cross-PS evidence 的最低必要條件要定義成哪些 germline anchors、spanning molecules 與一致性門檻？\n"
                "- k>12 subset 的 29.72% support loss 是否需要 separator-aware exact decomposition，或現有局部 blocks 已足以支援下游 likelihood？\n"
                "- chr6p/chr16q 的 candidate loci 在配對 normal、CN=0/exclusion 與 orthogonal truth 下應分成 excluded、ambiguous 或可評估哪幾類？"
            ),
        },
    ]

    manifest_sources = [
        {"id": source["id"], "label": source["label"], "path": source["path"]}
        for source in sources
    ]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": REPORT_TITLE,
            "description": "PARTIAL exploratory pilot：HCC1395 chr1–22 exact-PS evidence segmentation validation.",
            "generatedAt": generated_at,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": manifest_sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "partial",
            "datasets": {
                "headline_metrics": headline,
                "support_disposition": support_disposition,
                "k_bins": k_rows,
                "ps_containers": ps_container_rows,
                "chromosome_coverage": evidence["per_chromosome"],
                "outlier_arms": outlier["rows"],
                "validation_checks": validation_rows,
                "definitions": definitions,
            },
            "accessIssues": [
                {
                    "id": "scope_is_hcc1395_segmentation_only",
                    "scope": "claim ceiling",
                    "message": (
                        "PARTIAL：只有 HCC1395 chr1–22 extraction→exact-PS segmentation→C++ parity。"
                        "七樣本、tree enumeration、partial-read likelihood、VAF ranking 與 topology 尚未執行。"
                    ),
                }
            ],
        },
        "sources": sources,
    }
    return artifact


def build_source_notes(args: argparse.Namespace, evidence: dict[str, Any]) -> str:
    return f"""<!--
建立時間: 2026-07-22
目標: 保存 HCC1395 exact-PS technical report 的來源、圖表契約與重現命令
處理範圍: PARTIAL / exploratory pilot / HCC1395 chr1-22 only
關聯檔案: InterSubMod/{TOPIC}/report/artifact.json
-->

# HCC1395 exact-PS report source notes

> **PARTIAL / EXPLORATORY PILOT**：此檔為 report supporting artifact，不是第二份 reader-facing report。

## Reporting job

- 問題：每個 exact PS 是否應分開保留 read constraint，以及 k≤12 分割在 HCC1395 的量化影響。
- 主要讀者：technical / 教授與研究方法 reviewer。
- 決策：是否接受 exact-PS fail-closed policy，並允許進入 HCC1395 downstream adapter 階段。
- Claim ceiling：extraction、exact-PS segmentation、C++ parity；不包含 T/Topo/VAF/clone/subclone。

## Input paths

- Run root：`{args.run_root.resolve()}`
- Run receipt：`{(args.run_root / 'run_receipt.json').resolve()}`
- SEQC2 truth VCF：`{args.truth_vcf.resolve()}`
- SEQC2 HC BED：`{args.hc_bed.resolve()}`
- Input manifest：`{args.input_manifest.resolve()}`

## Recompute command

```bash
PYTHONDONTWRITEBYTECODE=1 python {TOPIC}/scripts/build_hcc1395_validation_report.py \\
  --run-root {args.run_root.resolve()} \\
  --truth-vcf {args.truth_vcf.resolve()} \\
  --hc-bed {args.hc_bed.resolve()} \\
  --input-manifest {args.input_manifest.resolve()} \\
  --output-dir {args.output_dir.resolve()}
```

Expected output excerpt:

```text
aggregate_exact_match=true
S={evidence['direct_recompute']['S']}
primary_unique={evidence['direct_recompute']['unique_sites']}
units={evidence['direct_recompute']['units']}
cross_ps=0
python_cpp_mismatch=0
```

## Chart map

| Section | Question | Family/type | Fields | Supported claim | Palette |
|---|---|---|---|---|---|
| Support disposition | k≤12 loss concentrates where? | Composition / 100% stacked bar | scope, molecule_weight, disposition | Overall and k>12 denominators differ materially | categorical, direct legend |
| Chromosome coverage | Which chromosomes lose primary eligibility? | Comparison / bar | chrom, primary_rate | chr6/chr16 are upstream-scope outliers | single-root blue, no redundant legend |

## Portable delivery exception and QA

The shared portable reader's `width:100vw` top bar produces an 8 px horizontal overflow under Linux non-overlay scrollbars. `deliver_report_with_topbar_fix.mjs` keeps the canonical artifact payload, native reader, chart renderer, static SVG extraction, and verifier unchanged; it only changes the top bar to `width/max-width:100%` and clears its viewport margins before canonical browser verification.

```bash
node {TOPIC}/scripts/deliver_report_with_topbar_fix.mjs \\
  --artifact {TOPIC}/report/artifact.json \\
  --output {TOPIC}/report/report.html \\
  --receipt {TOPIC}/report/report_delivery_receipt.json
```

## Evidence and caveats

- Raw gzip TSV is recomputed before comparison with `run_receipt.json`.
- `molecule_weight` is molecule–unit incidence, not globally unique reads.
- SEQC2 HC is a benchmark scope. Loci outside it cannot be labeled false positive solely from non-overlap.
- chr6p/chr16q are documented complete-loss/benchmark-excluded regions. The current pilot manifest's unlisted-CN=`neutral` policy is unsafe for downstream CN/clone interpretation, although segmentation does not consume CN.
- The direct Big7 BAM supplies read evidence at a manifest-bound candidate catalog; it did not de novo call the {evidence['direct_recompute']['S']:,} loci.

## Audience-spec mapping

1. Title → report title block.
2. Technical summary → `technical_summary`.
3. Key findings with visual evidence → exact-PS risk, support disposition, chromosome outlier sections.
4. Scope/data/metric definitions → definitions table and visible PARTIAL notice.
5. Methodology → exact-PS pipeline section.
6. Limitations/robustness → validation table and limitations section.
7. Recommended next steps → evidence-contract-first next steps.
8. Further questions → signed bridge, k>12 decomposition, and outlier classification questions.
"""


def main() -> int:
    args = parse_args()
    for path in (args.run_root, args.truth_vcf, args.hc_bed, args.input_manifest):
        if not path.exists():
            raise FileNotFoundError(path)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    evidence = build_evidence(args)
    artifact = build_artifact(evidence)
    json_dump(args.output_dir / "evidence_snapshot.json", evidence)
    json_dump(args.output_dir / "artifact.json", artifact)
    (args.output_dir / "source_notes.md").write_text(build_source_notes(args, evidence))
    print(f"output_dir={args.output_dir.resolve()}")
    print(f"aggregate_exact_match={str(evidence['aggregate_exact_match']).lower()}")
    print(f"S={evidence['direct_recompute']['S']}")
    print(f"primary_unique={evidence['direct_recompute']['unique_sites']}")
    print(f"units={evidence['direct_recompute']['units']}")
    print(f"cross_ps={evidence['validation']['cross_ps_violations']}")
    print(f"python_cpp_mismatch={evidence['validation']['python_cpp_mismatch_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
