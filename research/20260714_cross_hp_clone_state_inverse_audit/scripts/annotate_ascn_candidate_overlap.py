#!/usr/bin/env python3
"""Annotate selected cross-HP candidates with a conservative SAVANA 1+1 magnitude screen."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SAMPLE_ALIAS = {"HCC1395_DORADO": "HCC1395"}


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def solution_stability(ranked_path: Path) -> dict[str, Any]:
    if not ranked_path.exists():
        return {"available": False, "pass": False, "reason": "ranked_solution_unavailable"}
    rows = sorted(read_tsv(ranked_path), key=lambda row: int(row["rank"]))
    best = rows[0]
    margin = None if len(rows) == 1 else float(rows[1]["distance"]) - float(best["distance"])
    passed = len(rows) == 1 or (margin is not None and margin >= 0.02)
    return {
        "available": True,
        "pass": passed,
        "n_solutions": len(rows),
        "best_purity": float(best["purity"]),
        "best_ploidy": float(best["ploidy"]),
        "best_distance": float(best["distance"]),
        "best_second_distance_margin": margin,
        "rule": "one solution or second_distance-best_distance >= 0.02",
    }


def segment_screen(record: dict[str, Any], rows: list[dict[str, str]]) -> dict[str, Any]:
    chrom = str(record["chrom"])
    query_start = int(record["start"]) - 10_000
    query_end = int(record["end"]) + 10_000
    containing = [
        row
        for row in rows
        if row["chromosome"] == chrom
        and int(float(row["start"])) <= query_start
        and int(float(row["end"])) >= query_end
    ]
    if len(containing) != 1:
        return {
            "pass": False,
            "reason": "region_plus_10kb_not_in_exactly_one_segment",
            "containing_segment_count": len(containing),
        }
    row = containing[0]
    total = float(row["copyNumber"])
    minor = float(row["minorAlleleCopyNumber"])
    major = total - minor
    length = int(float(row["end"])) - int(float(row["start"])) + 1
    bin_count = int(float(row["bin_count"]))
    het = int(float(row["no_hetSNPs"]))
    checks = {
        "total_cn": abs(total - 2.0) <= 0.25,
        "minor_cn": abs(minor - 1.0) <= 0.25,
        "major_cn": abs(major - 1.0) <= 0.25,
        "segment_length": length >= 100_000,
        "bin_count": bin_count >= 10,
        "het_snps": het >= 20,
    }
    return {
        "pass": all(checks.values()),
        "reason": "ok" if all(checks.values()) else "one_or_more_segment_checks_failed",
        "segment_id": row["segment_id"],
        "segment_start": int(float(row["start"])),
        "segment_end": int(float(row["end"])),
        "total_cn": total,
        "minor_cn": minor,
        "major_cn": major,
        "segment_length": length,
        "bin_count": bin_count,
        "no_hetSNPs": het,
        "checks": checks,
        "orientation_caveat": "major/minor magnitudes are not oriented to HP1 versus HP2",
    }


def category_counts(records: list[dict[str, Any]], key: str) -> dict[str, int]:
    predicates = {
        "topology_invariant": lambda row: bool(row["direct_sister_shape_invariant"]),
        "tree_unique": lambda row: bool(row["direct_sister_tree_unique"]),
        "topology_invariant_analysis_complete": lambda row: bool(
            row["direct_sister_shape_invariant_analysis_complete"]
        ),
        "tree_unique_analysis_complete": lambda row: bool(
            row["direct_sister_tree_unique_analysis_complete"]
        ),
        "observed_hp1_direct_only_hp2_sister_only": lambda row: bool(
            row.get("observed_hp1_direct_only_hp2_sister_only")
        ),
        "topology_invariant_any_collision": lambda row: bool(
            row["direct_sister_shape_invariant"] and row["collision_positions"]
        ),
        "tree_unique_any_collision": lambda row: bool(
            row["direct_sister_tree_unique"] and row["collision_positions"]
        ),
    }
    return {
        name: sum(1 for row in records if predicate(row) and row["ascn_screen"][key])
        for name, predicate in predicates.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--cn-root", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    args = parser.parse_args()

    audit = load_json(args.audit)
    records = audit["candidate_records"]
    cache: dict[
        str, tuple[list[dict[str, str]], dict[str, Any], str, str, str | None, str | None]
    ] = {}
    annotated: list[dict[str, Any]] = []
    for record in records:
        dataset = str(record["dataset"])
        source_sample = SAMPLE_ALIAS.get(dataset, dataset)
        if source_sample not in cache:
            base = args.cn_root / source_sample / "savana_wgs" / "cna_normalhet"
            segment_path = base / f"{source_sample}_segmented_absolute_copy_number.tsv"
            ranked_path = base / f"{source_sample}_ranked_solutions.tsv"
            rows = read_tsv(segment_path) if segment_path.exists() else []
            cache[source_sample] = (
                rows,
                solution_stability(ranked_path),
                str(segment_path.resolve()) if segment_path.exists() else str(segment_path),
                str(ranked_path.resolve()) if ranked_path.exists() else str(ranked_path),
                sha256_file(segment_path),
                sha256_file(ranked_path),
            )
        rows, stability, segment_path, ranked_path, _, _ = cache[source_sample]
        segment = segment_screen(record, rows) if rows else {
            "pass": False,
            "reason": "segment_source_unavailable",
        }
        canonical_neutral_required = dataset in {"HCC1395", "HCC1395_DORADO"}
        canonical_neutral_pass = (
            record["cn_label"] == "neutral" if canonical_neutral_required else True
        )
        magnitude_pass = bool(segment["pass"] and canonical_neutral_pass)
        stable_pass = bool(magnitude_pass and stability["pass"])
        annotated.append(
            {
                **record,
                "ascn_screen": {
                    "source_sample": source_sample,
                    "segment_path": segment_path,
                    "ranked_solution_path": ranked_path,
                    "segment": segment,
                    "solution_stability": stability,
                    "canonical_neutral_required": canonical_neutral_required,
                    "canonical_neutral_pass": canonical_neutral_pass,
                    "magnitude_pass": magnitude_pass,
                    "stable_pass": stable_pass,
                    "screening_only": True,
                    "not_proven": [
                        "HP1/HP2 orientation",
                        "absence of subclonal CNA",
                        "cell-level homolog pairing",
                    ],
                },
            }
        )

    # Recompute the declared gate over the complete 22,779-region dual-HP mother set,
    # not only over the compact candidate table.
    overall_by_dataset: dict[str, dict[str, int]] = {}
    run_dir = Path(audit["input_run"]["path"])
    overall_total = 0
    overall_magnitude = 0
    overall_stable = 0
    for sample_dir in sorted((run_dir / "samples").iterdir()):
        if not sample_dir.is_dir():
            continue
        dataset = sample_dir.name
        source_sample = SAMPLE_ALIAS.get(dataset, dataset)
        if source_sample not in cache:
            base = args.cn_root / source_sample / "savana_wgs" / "cna_normalhet"
            segment_path_obj = base / f"{source_sample}_segmented_absolute_copy_number.tsv"
            ranked_path_obj = base / f"{source_sample}_ranked_solutions.tsv"
            rows = read_tsv(segment_path_obj) if segment_path_obj.exists() else []
            cache[source_sample] = (
                rows,
                solution_stability(ranked_path_obj),
                str(segment_path_obj.resolve()) if segment_path_obj.exists() else str(segment_path_obj),
                str(ranked_path_obj.resolve()) if ranked_path_obj.exists() else str(ranked_path_obj),
                sha256_file(segment_path_obj),
                sha256_file(ranked_path_obj),
            )
        rows, stability, _, _, _, _ = cache[source_sample]
        region_view = load_json(sample_dir / f"layered_region_view_{dataset}.json")
        dataset_counts = {"dual_hp_total": 0, "magnitude_pass": 0, "stable_pass": 0}
        for region in region_view.get("regions", []):
            primary_families = {
                str(lineage.get("family"))
                for lineage in region.get("lineages", [])
                if lineage.get("is_primary_lineage") and lineage.get("mutation_bearing")
            }
            if not {"1", "2"}.issubset(primary_families):
                continue
            dataset_counts["dual_hp_total"] += 1
            overall_total += 1
            screen_record = {
                "chrom": region["chrom"],
                "start": int(region["start"]),
                "end": int(region["end"]),
            }
            segment = segment_screen(screen_record, rows) if rows else {
                "pass": False,
                "reason": "segment_source_unavailable",
            }
            canonical_neutral_required = dataset in {"HCC1395", "HCC1395_DORADO"}
            canonical_neutral_pass = (
                region.get("cn") == "neutral" if canonical_neutral_required else True
            )
            magnitude_pass = bool(segment["pass"] and canonical_neutral_pass)
            stable_pass = bool(magnitude_pass and stability["pass"])
            if magnitude_pass:
                dataset_counts["magnitude_pass"] += 1
                overall_magnitude += 1
            if stable_pass:
                dataset_counts["stable_pass"] += 1
                overall_stable += 1
        overall_by_dataset[dataset] = dataset_counts

    source_manifest = {
        sample: {
            "segment_path": values[2],
            "ranked_solution_path": values[3],
            "segment_sha256": values[4],
            "ranked_solution_sha256": values[5],
        }
        for sample, values in sorted(cache.items())
    }
    payload = {
        "schema_name": "intersubmod.cross_hp_candidate_ascn_overlap",
        "schema_version": "1.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_audit": str(args.audit.resolve()),
        "cn_root": str(args.cn_root.resolve()),
        "gate": {
            "region_margin_bp": 10_000,
            "total_cn_tolerance": 0.25,
            "minor_cn_tolerance": 0.25,
            "major_cn_tolerance": 0.25,
            "minimum_segment_bp": 100_000,
            "minimum_bins": 10,
            "minimum_het_snps": 20,
            "solution_stability": "one solution or best-second distance margin >=0.02",
            "HCC1395_extra_gate": "canonical SEQC2-neutral label",
        },
        "category_magnitude_pass_counts": category_counts(annotated, "magnitude_pass"),
        "category_stable_pass_counts": category_counts(annotated, "stable_pass"),
        "dual_hp_mother_set_screen": {
            "total": overall_total,
            "magnitude_pass": overall_magnitude,
            "stable_pass": overall_stable,
            "by_dataset": overall_by_dataset,
            "unit": "dataset-region; HCC1395 and HCC1395_DORADO are not independent biological samples",
        },
        "source_manifest": source_manifest,
        "H2009_chr1_120007237_120040749": next(
            (
                row["ascn_screen"]
                for row in annotated
                if row["dataset"] == "H2009" and row["region"] == "chr1:120007237-120040749"
            ),
            None,
        ),
        "candidate_records": annotated,
        "claim_ceiling": (
            "SAVANA major/minor copy-number magnitude screen only; it does not orient CN to HP1/HP2 "
            "or establish biological clone pairing; ranked-output row count is not proof of global "
            "purity/ploidy solution-space uniqueness."
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    fields = [
        "dataset", "region", "n_sSNV", "direct_sister_shape_invariant",
        "direct_sister_tree_unique", "observed_hp1_direct_only_hp2_sister_only",
        "collision_positions", "cn_label", "magnitude_pass", "stable_pass",
        "segment_id", "total_cn", "minor_cn", "major_cn", "solution_margin",
    ]
    with args.out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in annotated:
            screen = row["ascn_screen"]
            segment = screen["segment"]
            stability = screen["solution_stability"]
            writer.writerow(
                {
                    "dataset": row["dataset"],
                    "region": row["region"],
                    "n_sSNV": row["n_sSNV"],
                    "direct_sister_shape_invariant": row["direct_sister_shape_invariant"],
                    "direct_sister_tree_unique": row["direct_sister_tree_unique"],
                    "observed_hp1_direct_only_hp2_sister_only": row.get(
                        "observed_hp1_direct_only_hp2_sister_only"
                    ),
                    "collision_positions": json.dumps(row["collision_positions"]),
                    "cn_label": row["cn_label"],
                    "magnitude_pass": screen["magnitude_pass"],
                    "stable_pass": screen["stable_pass"],
                    "segment_id": segment.get("segment_id"),
                    "total_cn": segment.get("total_cn"),
                    "minor_cn": segment.get("minor_cn"),
                    "major_cn": segment.get("major_cn"),
                    "solution_margin": stability.get("best_second_distance_margin"),
                }
            )
    print(json.dumps(
        {
            "magnitude_pass": payload["category_magnitude_pass_counts"],
            "stable_pass": payload["category_stable_pass_counts"],
            "dual_hp_mother_set_screen": payload["dual_hp_mother_set_screen"],
            "H2009_candidate": payload["H2009_chr1_120007237_120040749"],
            "outputs": [str(args.out_json.resolve()), str(args.out_tsv.resolve())],
        },
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()
