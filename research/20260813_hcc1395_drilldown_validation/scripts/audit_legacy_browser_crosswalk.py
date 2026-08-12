#!/usr/bin/env python3
"""Audit the 2026-07-26 methyl browser against a current drilldown bundle.

This script deliberately treats the two axis systems as non-equivalent.  It
reports coordinate overlap and a descriptive cross-tab, but it never labels
one classification as ground truth for the other.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


LEGACY_A = "A_ALLELE_clean"
LEGACY_CLASSES = (
    LEGACY_A,
    "B_HP_ASM",
    "C_other",
    "D_split_unaligned",
    "E_unstable",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def boot_data(index_html: Path) -> dict:
    text = index_html.read_text(encoding="utf-8")
    match = re.search(
        r'<script type="application/json" id="bootData">(.*?)</script>',
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise ValueError(f"bootData not found: {index_html}")
    return json.loads(match.group(1))


def current_loci(l1: dict) -> dict[tuple[str, int], int]:
    loci: dict[tuple[str, int], int] = {}
    previous_chrom = -1
    position = 0
    for idx, chrom_idx in enumerate(l1["chrom"]):
        if chrom_idx != previous_chrom:
            position = int(l1["dpos"][idx])
            previous_chrom = chrom_idx
        else:
            position += int(l1["dpos"][idx])
        loci[(l1["chroms"][chrom_idx], position)] = int(l1["axis"][idx])
    if len(loci) != int(l1["n"]):
        raise ValueError("decoded current loci are not unique")
    return loci


def legacy_loci(tsv: Path) -> dict[tuple[str, int], dict[str, str]]:
    with tsv.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    loci = {(row["chrom"], int(row["snv_pos"])): row for row in rows}
    if len(loci) != len(rows):
        raise ValueError("legacy locus coordinates are not unique")
    return loci


def pct(num: int, den: int) -> float | None:
    return round(num / den * 100, 6) if den else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-dir", required=True, type=Path)
    parser.add_argument("--current-bundle", required=True, type=Path)
    parser.add_argument("--scratchpad-root", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    legacy_tsv = args.legacy_dir / "methyl_class_x_linkage_annotation.tsv"
    legacy_data_path = args.legacy_dir / "20260726_ssnv_branch_x_methyl_browser.data.json"
    legacy_build = args.legacy_dir / "build_browser.py"
    legacy_html = args.legacy_dir / "20260726_ssnv_branch_x_methyl_browser.standalone.html"
    current_index = args.current_bundle / "index.html"

    legacy = legacy_loci(legacy_tsv)
    boot = boot_data(current_index)
    current = current_loci(boot["l1"])
    legacy_data = json.loads(legacy_data_path.read_text(encoding="utf-8"))

    legacy_keys = set(legacy)
    current_keys = set(current)
    shared = legacy_keys & current_keys

    cross_rows = []
    for legacy_class in LEGACY_CLASSES:
        keys = [key for key in shared if legacy[key]["methyl_class"] == legacy_class]
        cross_rows.append(
            {
                "legacy_class": legacy_class,
                "shared_loci": len(keys),
                "current_no_ism_axis_code_0": sum(current[key] == 0 for key in keys),
                "current_tested_none_axis_code_8": sum(current[key] == 8 for key in keys),
                "current_hp_raw_significant": sum(bool(current[key] & 1) for key in keys),
                "current_alt_raw_significant": sum(bool(current[key] & 2) for key in keys),
                "current_lineage_raw_significant": sum(bool(current[key] & 4) for key in keys),
            }
        )

    legacy_a_keys = {key for key in shared if legacy[key]["methyl_class"] == LEGACY_A}
    current_alt_keys = {key for key in shared if current[key] & 2}
    both = legacy_a_keys & current_alt_keys
    union = legacy_a_keys | current_alt_keys
    legacy_a_no_l4 = sum(current[key] == 0 for key in legacy_a_keys)
    legacy_a_l4_available = len(legacy_a_keys) - legacy_a_no_l4

    regions_total = int(legacy_data["regions_total"])
    regions_with_a = int(legacy_data["regions_with_A"])
    candidate_regions = int(legacy_data["n_branch_methyl_regions"])
    displayed_cases = int(legacy_data["n_cases"])
    funnel_rows = [
        {"stage": "legacy_regions_with_methyl_annotation", "n": regions_total,
         "pct_of_first_stage": 100.0, "pct_of_previous_stage": 100.0},
        {"stage": "legacy_regions_with_at_least_one_A", "n": regions_with_a,
         "pct_of_first_stage": pct(regions_with_a, regions_total),
         "pct_of_previous_stage": pct(regions_with_a, regions_total)},
        {"stage": "legacy_branch_plus_methyl_candidates", "n": candidate_regions,
         "pct_of_first_stage": pct(candidate_regions, regions_total),
         "pct_of_previous_stage": pct(candidate_regions, regions_with_a)},
        {"stage": "legacy_displayed_heatmap_cases", "n": displayed_cases,
         "pct_of_first_stage": pct(displayed_cases, regions_total),
         "pct_of_previous_stage": pct(displayed_cases, candidate_regions)},
    ]

    build_text = legacy_build.read_text(encoding="utf-8")
    absolute_tmp_inputs = bool(re.search(r'Path\("/tmp/', build_text))
    input_hashes = {}
    if args.scratchpad_root:
        for name in ("methyl_backbone.json", "branch_methyl_cases.json", "methyl_regions.json"):
            path = args.scratchpad_root / name
            input_hashes[name] = {
                "path": str(path),
                "exists": path.is_file(),
                "bytes": path.stat().st_size if path.is_file() else None,
                "sha256": sha256(path) if path.is_file() else None,
            }

    result = {
        "status": "NON_EQUIVALENT_SCHEMAS",
        "claim_ceiling": "descriptive coordinate crosswalk only",
        "sample": {"legacy_inferred": "HCC1395", "current": boot.get("sample")},
        "inputs": {
            "legacy_tsv": {"path": str(legacy_tsv), "sha256": sha256(legacy_tsv)},
            "legacy_data": {"path": str(legacy_data_path), "sha256": sha256(legacy_data_path)},
            "legacy_build": {"path": str(legacy_build), "sha256": sha256(legacy_build)},
            "legacy_html": {"path": str(legacy_html), "bytes": legacy_html.stat().st_size,
                            "sha256": sha256(legacy_html)},
            "current_index": {"path": str(current_index), "sha256": sha256(current_index)},
            "scratchpad": input_hashes,
        },
        "reproducibility": {
            "build_uses_absolute_tmp_inputs": absolute_tmp_inputs,
            "stable_manifest_embedded": False,
            "assessment": "legacy HTML is viewable, but its builder does not use a stable repository input manifest",
        },
        "coordinate_crosswalk": {
            "legacy_unique_loci": len(legacy_keys),
            "current_unique_loci": len(current_keys),
            "shared_loci": len(shared),
            "legacy_only_loci": len(legacy_keys - current_keys),
            "current_only_loci": len(current_keys - legacy_keys),
            "shared_pct_of_legacy": pct(len(shared), len(legacy_keys)),
            "shared_pct_of_current": pct(len(shared), len(current_keys)),
            "coordinate_jaccard": round(
                len(shared) / len(legacy_keys | current_keys), 6
            ) if legacy_keys or current_keys else None,
        },
        "current_axis_encoding_provenance": {
            "source": "immutable pre-hardening HCC1395_v1 bundle",
            "observed_code_range": [min(current.values()), max(current.values())],
            "axis_untested_code_16_available": False,
            "note": (
                "Codes in this crosswalk come from the immutable legacy v1 index (0-8). "
                "The source-level AXIS_UNTESTED=16 third state was added later and is not "
                "retroactively present in this bundle."
            ),
        },
        "legacy_class_counts_all": dict(Counter(row["methyl_class"] for row in legacy.values())),
        "current_axis_code_counts_all": {
            str(key): value for key, value in sorted(Counter(current.values()).items())
        },
        "legacy_A_vs_current_ALT_descriptive_only": {
            "legacy_A_shared": len(legacy_a_keys),
            "legacy_A_no_L4_axis_code_0": legacy_a_no_l4,
            "legacy_A_L4_available": legacy_a_l4_available,
            "current_ALT_shared": len(current_alt_keys),
            "intersection": len(both),
            "jaccard": round(len(both) / len(union), 6) if union else None,
            "denominator_note": (
                f"Main-report ALT/HP percentages use {legacy_a_l4_available} L4-available "
                f"legacy-A loci = {len(legacy_a_keys)} shared legacy-A loci - "
                f"{legacy_a_no_l4} current axis-code-0 loci. Jaccard uses the full "
                "shared-coordinate sets and therefore has a different denominator."
            ),
            "warning": (
                "Different locus universes, gates and statistics; this is not accuracy, "
                "validation, or a class-to-axis mapping."
            ),
        },
        "legacy_funnel": funnel_rows,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "legacy_browser_crosswalk.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    cross_path = args.output_dir / "legacy_current_axis_crosswalk.csv"
    with cross_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cross_rows[0]))
        writer.writeheader()
        writer.writerows(cross_rows)

    funnel_path = args.output_dir / "legacy_browser_funnel.csv"
    with funnel_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(funnel_rows[0]))
        writer.writeheader()
        writer.writerows(funnel_rows)

    summary = {
        "status": result["status"],
        **result["coordinate_crosswalk"],
        "legacy_A_current_ALT_intersection": len(both),
        "legacy_A_current_ALT_jaccard": result["legacy_A_vs_current_ALT_descriptive_only"]["jaccard"],
        "coordinate_jaccard": result["coordinate_crosswalk"]["coordinate_jaccard"],
        "legacy_displayed_cases": displayed_cases,
        "outputs": [str(json_path), str(cross_path), str(funnel_path)],
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
