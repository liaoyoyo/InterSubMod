#!/usr/bin/env python3
"""Build the current InterSubMod layered-reconstruction observation report.

The report is deliberately a snapshot of two distinct evidence tracks:

1. the hash-verified 2026-07-10 historical engineering snapshot, and
2. the clean LongPhase-S production-sidecar gate that may still be running.

It never upgrades missing read-AF or methylation evidence to zero and never
describes regional mutation-state candidates as confirmed cell clones.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import html
import importlib.util
import json
import re
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from math import prod
from pathlib import Path


SAMPLE_ORDER = [
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
]

TREE_CLASS_LABELS = {
    "determined": "模型內唯一",
    "ambiguous_structure(多完成/多結構)": "多候選／歧義",
    "capped(太密;枚舉未完)": "solver capped",
    "recurrence_required": "需要 recurrence",
}

SOURCE_LABELS = {
    "verification": "Source: InterSubMod hash-verified historical engineering snapshot<br>File: verification_summary.json",
    "region": "Source: InterSubMod layered-v2 region view<br>File: layered_region_view_&lt;sample&gt;.json",
    "layered": "Source: InterSubMod layered-v2 reconstruction units<br>File: layered_reconstruction_&lt;sample&gt;.json",
    "mlhp": "Source: InterSubMod layered-v2 read observations<br>Files: mlhp_part_1.json … mlhp_part_5.json",
    "vcf": "Source: historical snapshot input: ClairS paired PASS call set<br>File: somatic_pass.vcf.gz; superseded for future canonical tree input",
    "production": "Source: normalized raw-all LongPhase-S production-sidecar evidence<br>Files: run_status.tsv, verification_summary.json, _SUCCESS and/or _FAILED",
    "read_af": "Source: report input audit<br>Finding: no canonical read-AF ordering artifact under the run root",
    "method": "Source: layered-v2 historical run contract and solver configuration<br>Files: input_manifest.json and run-scoped method scripts",
    "seqc2": "Source: SEQC2 HCC1395 benchmark artifacts<br>Files: High-Confidence_Regions_v1.2.bed and CN benchmark BED",
    "scope": "Source: historical snapshot validation boundary<br>File: VALIDATION_SCOPE.md",
    "code_hash": "Source: historical run code provenance check<br>File: code.sha256",
    "derived": "Source: report builder deterministic recomputation<br>Inputs: hash-verified historical region, layered, MLHP and VCF artifacts",
}

REPORT_STEM = "20260711_分層重建全景數據觀察_01"
DATA_SNAPSHOT_NAME = f"{REPORT_STEM}.data.json"
SOURCE_NOTES_NAME = f"{REPORT_STEM}.source_notes.json"
VALIDATION_NAME = f"{REPORT_STEM}.validation.json"
BROWSER_QA_NAME = "audit/after/metrics.json"


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def dump_json(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def pct(numerator: float, denominator: float):
    return numerator / denominator if denominator else None


def median(values):
    return statistics.median(values) if values else None


def percentile(values, q):
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def format_value(value, kind="integer"):
    if value is None:
        return "未評估"
    if kind == "percent":
        return f"{value * 100:.2f}%"
    if kind == "percent1":
        return f"{value * 100:.1f}%"
    if kind == "decimal2":
        return f"{value:.2f}"
    if kind == "decimal1":
        return f"{value:.1f}"
    return f"{int(value):,}"


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha_manifest(manifest: Path):
    checked = []
    if not manifest.is_file():
        return {"manifest": str(manifest), "pass": False, "checked": checked, "error": "missing"}
    for raw in manifest.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        relative = relative.lstrip("* ")
        target = manifest.parent / relative
        actual = sha256(target) if target.is_file() else None
        checked.append({
            "file": relative,
            "expected": expected,
            "actual": actual,
            "pass": actual == expected,
        })
    return {
        "manifest": str(manifest),
        "pass": bool(checked) and all(item["pass"] for item in checked),
        "checked": checked,
    }


def count_vcf_snv_contigs(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    counts = Counter()
    nonpass_snv = 0
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) > 4 and len(fields[3]) == 1 and len(fields[4]) == 1:
                counts[fields[0]] += 1
                if len(fields) < 7 or fields[6] != "PASS":
                    nonpass_snv += 1
    return {"contigs": dict(counts), "nonpass_snv": nonpass_snv}


def load_mlhp_groups(sample_dir: Path):
    groups = []
    params = None
    parts = sorted(sample_dir.glob("mlhp_part_*.json"))
    for part in parts:
        payload = load_json(part)
        groups.extend(payload.get("groups", []))
        if params is None:
            params = payload.get("params", {})
    return groups, parts, params or {}


def region_key_from_group(group):
    return f"{group['chrom']}:{group['start']}-{group['end']}"


def c_bin(value):
    """Return the professor-facing joint candidate-combination bin."""
    if value is None:
        return "incomplete"
    return f"C={value}" if 1 <= value <= 6 else "C>6"


def region_candidate_topology(regions, detail_by_key):
    """Derive region-level HP/H3, exact-combination and topology summaries.

    Primary HP1/HP2 reconstruction units are solved independently.  For an
    analysis-complete region, their joint exact-tree combination count is the
    Cartesian product of per-unit ``n_trees``; the analogous topology count is
    the product of exact unlabeled-shape counts.  Any incomplete/capped primary
    unit makes both region-level counts unavailable rather than zero.
    """
    hp_h3 = Counter()
    c_counts = {scope: Counter() for scope in ("single", "double", "all")}
    topology = Counter()
    topology_by_hp = {scope: Counter() for scope in ("single", "double")}
    hidden = Counter()
    hidden_by_topology = defaultdict(Counter)
    hidden_by_topology_by_hp = {
        scope: defaultdict(Counter) for scope in ("single", "double")
    }
    complete_joint_c = []
    complete_joint_topology = []
    invalid = []
    capped_label_discordance = []

    for region in regions:
        primary = [line for line in region.get("lineages", []) if line.get("is_primary_lineage")]
        primary_count = len(primary)
        # Canonical H3? denominator counts emitted mutation-bearing H3 auxiliary
        # units only.  A reference-only family=3 record must remain H3- here.
        has_h3 = int(region.get("n_H3_auxiliary", 0)) > 0
        primary_families = {str(line.get("family")) for line in primary}
        any_primary_capped = any(line.get("capped") for line in primary)
        if any_primary_capped and region.get("region_determinacy") != "has_capped":
            capped_label_discordance.append({
                "region": region.get("region"),
                "legacy_region_determinacy": region.get("region_determinacy"),
            })
        if primary_count == 0:
            scope = "none"
        elif primary_count == 1 and primary_families <= {"1", "2"}:
            scope = "single"
        elif primary_count == 2 and primary_families == {"1", "2"}:
            scope = "double"
        else:
            scope = "other"
            invalid.append({"region": region.get("region"), "reason": f"primary_count={primary_count}"})
        hp_h3[f"{scope}_{'with_h3' if has_h3 else 'without_h3'}"] += 1
        if scope not in {"single", "double"}:
            continue

        complete = True
        for line in primary:
            detail = detail_by_key.get((region.get("region"), str(line.get("family"))))
            line_complete = bool(
                detail
                and line.get("analysis_candidate_set_complete") is True
                and not line.get("capped")
                and line.get("verification_status") == "full_pass"
                and line.get("verification_complete") is True
                and detail.get("analysis_trees_generated") == line.get("n_trees")
                and isinstance(line.get("n_trees"), int)
                and line.get("n_trees", 0) >= 1
                and isinstance(line.get("n_distinct_shapes_exact"), int)
                and line.get("n_distinct_shapes_exact", 0) >= 1
            )
            complete = complete and line_complete
        if not complete:
            c_counts[scope]["incomplete"] += 1
            c_counts["all"]["incomplete"] += 1
            topology["incomplete"] += 1
            topology_by_hp[scope]["incomplete"] += 1
            hidden["incomplete"] += 1
            continue

        joint_c = prod(int(line["n_trees"]) for line in primary)
        joint_topology = prod(int(line["n_distinct_shapes_exact"]) for line in primary)
        hidden_nodes = sum(int(line.get("n_hidden", 0)) for line in primary)
        if not (joint_c >= joint_topology >= 1):
            invalid.append({
                "region": region.get("region"),
                "reason": f"joint_C={joint_c}, joint_Topo={joint_topology}",
            })
            continue

        bucket = c_bin(joint_c)
        c_counts[scope][bucket] += 1
        c_counts["all"][bucket] += 1
        complete_joint_c.append(joint_c)
        complete_joint_topology.append(joint_topology)

        if joint_c == 1 and joint_topology == 1:
            topology_class = "exact_and_topology_unique"
        elif joint_c > 1 and joint_topology == 1:
            topology_class = "topology_unique_exact_multiple"
        elif joint_c > 1 and joint_topology > 1:
            topology_class = "topology_multiple_exact_multiple"
        else:
            topology_class = "impossible_exact_unique_topology_multiple"
            invalid.append({
                "region": region.get("region"),
                "reason": f"impossible class joint_C={joint_c}, joint_Topo={joint_topology}",
            })
        topology[topology_class] += 1
        topology_by_hp[scope][topology_class] += 1
        hidden_class = "hidden_zero" if hidden_nodes == 0 else "hidden_positive"
        hidden[hidden_class] += 1
        hidden_by_topology[topology_class][hidden_class] += 1
        hidden_by_topology_by_hp[scope][topology_class][hidden_class] += 1

    c_order = ["C=1", "C=2", "C=3", "C=4", "C=5", "C=6", "C>6", "incomplete"]
    topology_order = [
        "exact_and_topology_unique",
        "topology_unique_exact_multiple",
        "topology_multiple_exact_multiple",
        "incomplete",
        "impossible_exact_unique_topology_multiple",
    ]
    hp_h3_order = [
        "single_without_h3", "single_with_h3", "double_without_h3", "double_with_h3",
        "none_without_h3", "none_with_h3", "other_without_h3", "other_with_h3",
    ]
    return {
        "definitions": {
            "C_region": "product of n_trees across analysis-complete primary HP1/HP2 units",
            "Topo_region": "product of n_distinct_shapes_exact across those units",
            "hidden_region": "sum of per-primary-unit minimal parsimony-inferred mutation-state nodes",
        },
        "hp_h3": {key: int(hp_h3.get(key, 0)) for key in hp_h3_order},
        "C_bins": {
            scope: {key: int(c_counts[scope].get(key, 0)) for key in c_order}
            for scope in ("single", "double", "all")
        },
        "topology_classes": {key: int(topology.get(key, 0)) for key in topology_order},
        "topology_by_hp": {
            scope: {key: int(topology_by_hp[scope].get(key, 0)) for key in topology_order}
            for scope in ("single", "double")
        },
        "hidden_classes": {
            "hidden_zero": int(hidden.get("hidden_zero", 0)),
            "hidden_positive": int(hidden.get("hidden_positive", 0)),
            "incomplete": int(hidden.get("incomplete", 0)),
        },
        "hidden_by_topology": {
            key: {
                "hidden_zero": int(hidden_by_topology[key].get("hidden_zero", 0)),
                "hidden_positive": int(hidden_by_topology[key].get("hidden_positive", 0)),
            }
            for key in topology_order[:3]
        },
        "hidden_by_topology_by_hp": {
            scope: {
                key: {
                    "hidden_zero": int(hidden_by_topology_by_hp[scope][key].get("hidden_zero", 0)),
                    "hidden_positive": int(hidden_by_topology_by_hp[scope][key].get("hidden_positive", 0)),
                }
                for key in topology_order[:3]
            }
            for scope in ("single", "double")
        },
        "complete_regions": len(complete_joint_c),
        "incomplete_regions": int(topology.get("incomplete", 0)),
        "joint_exact_combinations_total": sum(complete_joint_c),
        "joint_topology_combinations_total": sum(complete_joint_topology),
        "max_C_region": max(complete_joint_c, default=0),
        "max_Topo_region": max(complete_joint_topology, default=0),
        "ranking_status": "not_evaluated",
        "ranked_regions": 0,
        "legacy_has_capped_label_count": sum(1 for region in regions if region.get("region_determinacy") == "has_capped"),
        "capped_label_discordance": capped_label_discordance,
        "invalid_states": invalid,
    }


def sample_metrics(sample_meta, verification_result, run_root: Path):
    sample = sample_meta["sample"]
    sample_dir = run_root / "samples" / sample
    region_path = sample_dir / f"layered_region_view_{sample}.json"
    layered_path = sample_dir / f"layered_reconstruction_{sample}.json"
    region_view = load_json(region_path)
    layered = load_json(layered_path)
    detail_by_key = {
        (unit.get("region"), str(unit.get("family"))): unit
        for unit in layered["detail"]
    }
    groups, group_parts, mlhp_params = load_mlhp_groups(sample_dir)
    regions = region_view["regions"]
    primary_units = [unit for unit in layered["detail"] if unit.get("is_primary_lineage")]
    noncapped_primary = [
        unit for unit in primary_units
        if not unit.get("capped") and unit.get("analysis_candidate_set_complete") is True
    ]

    funnel = verification_result["metrics"]["funnel"]
    roles = verification_result["metrics"]["roles"]
    region_determinacy = verification_result["metrics"]["region_determinacy"]
    hp_multiplicity = {str(key): value for key, value in verification_result["metrics"]["hp_multiplicity"].items()}

    k_counts = Counter(int(group["n_sSNV"]) for group in groups)
    k_sites = Counter()
    for key, count in k_counts.items():
        k_sites[key] = key * count
    natural_k8 = sum(1 for group in groups if group["n_sSNV"] == 8 and group.get("pre_cap_n_sSNV") == 8)
    compressed_k8 = sum(1 for group in groups if group["n_sSNV"] == 8 and group.get("pre_cap_n_sSNV", 8) > 8)

    raw_full_observations = sum(int(group.get("n_full_cov_reads", 0)) for group in groups)
    full_pattern_observations = sum(
        sum(int(value) for value in (group.get("populations") or {}).values()) for group in groups
    )
    partial_pattern_observations = sum(
        int(value)
        for group in groups
        for pattern, value in (group.get("subread_groups") or {}).items()
        if "X" in pattern
    )
    broad_hp_observations = sum(
        int(value) for group in groups for value in (group.get("reads_by_hp") or {}).values()
    )
    regions_with_full = sum(1 for group in groups if group.get("n_full_cov_reads", 0) > 0)
    regions_with_alt_population = sum(1 for group in groups if group.get("n_populations_with_ALT", 0) > 0)

    primary_full_partial = sum(
        1 for unit in primary_units if unit.get("n_full_pops", 0) > 0 and unit.get("n_partial", 0) > 0
    )
    primary_partial_only = sum(
        1 for unit in primary_units if unit.get("n_full_pops", 0) == 0 and unit.get("n_partial", 0) > 0
    )
    primary_full_only = sum(
        1 for unit in primary_units if unit.get("n_full_pops", 0) > 0 and unit.get("n_partial", 0) == 0
    )

    l1_counts = Counter(unit["L1_class"] for unit in primary_units)
    exact_tree_unique = sum(1 for unit in noncapped_primary if unit.get("L1_class") == "determined")
    shape_unique = sum(1 for unit in noncapped_primary if unit.get("n_distinct_shapes_exact") == 1)
    total_candidate_trees = sum(int(unit.get("n_trees", 0)) for unit in noncapped_primary)
    hidden_total_all_primary = sum(int(unit.get("n_hidden", 0)) for unit in primary_units)
    hidden_total_complete = sum(int(unit.get("n_hidden", 0)) for unit in noncapped_primary)
    hidden_unit_bins = Counter(
        str(value) if value < 5 else "5+"
        for value in (int(unit.get("n_hidden", 0)) for unit in noncapped_primary)
    )
    candidate_counts = [int(unit.get("n_trees", 0)) for unit in noncapped_primary]

    regions_with_primary = [region for region in regions if region.get("n_primary_lineages", 0) > 0]
    candidate_topology = region_candidate_topology(regions, detail_by_key)
    region_all_shape_unique = 0
    for region in regions_with_primary:
        lineages = [line for line in region.get("lineages", []) if line.get("is_primary_lineage")]
        if lineages and all(
            line.get("analysis_candidate_set_complete") is True and line.get("n_distinct_shapes_exact") == 1
            for line in lineages
        ):
            region_all_shape_unique += 1

    mlhp_region_keys = {region_key_from_group(group) for group in groups}
    view_region_keys = {region["region"] for region in regions}
    missing_regions = sorted(mlhp_region_keys - view_region_keys)
    missing_region_detail = [
        {
            "region": region_key_from_group(group),
            "k": int(group["n_sSNV"]),
            "n_full_cov_reads": int(group.get("n_full_cov_reads", 0)),
            "n_populations": int(group.get("n_populations", 0)),
        }
        for group in groups
        if region_key_from_group(group) in set(missing_regions)
    ]
    extra_regions = sorted(view_region_keys - mlhp_region_keys)
    mlhp_site_sum = sum(int(group["n_sSNV"]) for group in groups)
    view_site_sum = sum(int(region["n_sSNV"]) for region in regions)

    vcf_path = Path(sample_meta["somatic_vcf"])
    vcf_census = count_vcf_snv_contigs(vcf_path)
    contig_counts = vcf_census["contigs"]
    chr_x = int(contig_counts.get("chrX", 0))
    chr_y = int(contig_counts.get("chrY", 0))
    nonautosomal_other = int(sum(
        value for chrom, value in contig_counts.items()
        if chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    ))
    well_explained_case = None
    if sample == "HCC1395":
        selected = next(
            (region for region in regions if region.get("region") == "chr2:221958321-222018361"),
            None,
        )
        if selected:
            primary = [line for line in selected.get("lineages", []) if line.get("is_primary_lineage")]
            well_explained_case = {
                "sample": sample,
                "region": selected["region"],
                "k": int(selected["n_sSNV"]),
                "hp_multiplicity": int(selected["hp_multiplicity"]),
                "n_full_cov_reads": int(selected.get("n_full_cov_reads", 0)),
                "families": [
                    {
                        "family": line["family"],
                        "n_trees": int(line["n_trees"]),
                        "n_hidden": int(line["n_hidden"]),
                        "n_reads": int(line["n_reads"]),
                        "n_full_pops": int(line["n_full_pops"]),
                        "n_partial": int(line["n_partial"]),
                    }
                    for line in primary
                ],
            }

    return {
        "sample": sample,
        "biological_id": verification_result["biological_id"],
        "pass": verification_result["pass"],
        "paths": {
            "region_view": str(region_path),
            "layered": str(layered_path),
            "mlhp_parts": [str(path) for path in group_parts],
            "vcf": str(vcf_path),
        },
        "site_funnel": {
            "universe": funnel["L1_all_pass_universe"],
            "autosomal": funnel["autosomal_chr1_22"],
            "out_of_scope": funnel["L2_out_of_scope_non_autosomal"],
            "singleton": funnel["L3_positional_singleton"],
            "multilocus_pre_cap": funnel["L4_multilocus_pre_cap_sSNV"],
            "cap_excluded": funnel["L5_cap_excluded_sSNV"],
            "read_unsupported": funnel["L5_read_unsupported_sSNV"],
            "retained": funnel["L6_retained_sSNV"],
            "scope_conservation": funnel["check_six_branch_conservation"],
        },
        "sex_chromosome_census": {
            "chrX": chr_x,
            "chrY": chr_y,
            "other": nonautosomal_other,
            "nonautosomal_fraction": pct(chr_x + chr_y + nonautosomal_other, funnel["L1_all_pass_universe"]),
            "vcf_total_matches": sum(contig_counts.values()) == funnel["L1_all_pass_universe"],
            "vcf_nonpass_snv": vcf_census["nonpass_snv"],
        },
        "regions": {
            "W_all_pre": funnel["L3_positional_singleton"] + funnel["n_multilocus_pre_cap_groups"],
            "W_k1": funnel["L3_positional_singleton"],
            "W_k_gt1": funnel["n_multilocus_pre_cap_groups"],
            "W_pre": funnel["n_multilocus_pre_cap_groups"],
            "W_ret": funnel["n_groups_retained"],
            "W_tree": verification_result["metrics"]["n_regions"],
            "W_primary": region_view["census"]["n_regions_with_primary_lineage"],
            "multi_HP": int(hp_multiplicity.get("2", 0)),
            "single_primary_HP": int(hp_multiplicity.get("1", 0)),
            "no_primary": int(hp_multiplicity.get("0", 0)),
            "missing_from_region_view": missing_regions,
            "missing_from_region_view_detail": missing_region_detail,
            "extra_in_region_view": extra_regions,
            "mlhp_site_sum": mlhp_site_sum,
            "region_view_site_sum": view_site_sum,
            "site_sum_delta": mlhp_site_sum - view_site_sum,
        },
        "k_distribution": {
            "counts": {str(k): int(k_counts.get(k, 0)) for k in range(2, 9)},
            "sites": {str(k): int(k_sites.get(k, 0)) for k in range(2, 9)},
            "natural_k8": natural_k8,
            "compressed_k8": compressed_k8,
            "mean_k": pct(mlhp_site_sum, len(groups)),
            "median_k": median([group["n_sSNV"] for group in groups]),
            "raw_pre_cap_max": max((int(group.get("pre_cap_n_sSNV", group["n_sSNV"])) for group in groups), default=0),
        },
        "read_evidence": {
            "regions_with_full": regions_with_full,
            "regions_without_full": len(groups) - regions_with_full,
            "regions_with_alt_population": regions_with_alt_population,
            "regions_without_alt_population": len(groups) - regions_with_alt_population,
            "raw_full_observations": raw_full_observations,
            "full_pattern_observations": full_pattern_observations,
            "partial_pattern_observations": partial_pattern_observations,
            "retained_pattern_observations": full_pattern_observations + partial_pattern_observations,
            "broad_hp_observations": broad_hp_observations,
            "primary_full_partial": primary_full_partial,
            "primary_partial_only": primary_partial_only,
            "primary_full_only": primary_full_only,
        },
        "units": {
            "primary": roles["primary"],
            "reference": roles["reference"],
            "H3_auxiliary": roles["H3_auxiliary"],
            "unphased_auxiliary": roles.get("unphased", 0),
            "determined": int(l1_counts.get("determined", 0)),
            "ambiguous": int(l1_counts.get("ambiguous_structure(多完成/多結構)", 0)),
            "solver_capped": int(l1_counts.get("capped(太密;枚舉未完)", 0)),
            "recurrence": int(l1_counts.get("recurrence_required", 0)),
            "noncapped_complete": len(noncapped_primary),
            "exact_tree_unique": exact_tree_unique,
            "shape_unique": shape_unique,
            "total_candidate_trees": total_candidate_trees,
            "candidate_tree_median": median(candidate_counts),
            "candidate_tree_p90": percentile(candidate_counts, 0.9),
            "hidden_total_all_primary": hidden_total_all_primary,
            "hidden_total_complete": hidden_total_complete,
            "hidden_unit_bins_complete": {
                key: int(hidden_unit_bins.get(key, 0)) for key in ("0", "1", "2", "3", "4", "5+")
            },
            "region_all_exact": int(region_determinacy.get("all_determined", 0)),
            "region_all_shape": region_all_shape_unique,
        },
        "region_candidate_topology": candidate_topology,
        "L2_cn": layered["L2_cn"],
        "cn_contract": {
            key: mlhp_params.get(key)
            for key in ("cn_data_available", "cn_source", "cn_interpretation")
        },
        "L3_methyl": layered["L3_methyl"],
        "verification_status": verification_result["metrics"]["verification_status"],
        "verification_failure": layered["L1_ssnv_algorithm"].get("n_verification_fail", 0),
        "method_params": {
            key: mlhp_params.get(key)
            for key in ("MINREAD", "MAX_SNV", "TIER_R", "MAPQ_MIN", "BASEQ_MIN")
        },
        "well_explained_case": well_explained_case,
    }


def read_production_status(production_root: Path):
    manifest_path = production_root / "input_manifest.json"
    status_path = production_root / "run_status.tsv"
    expected = 0
    manifest_sample_names = []
    if manifest_path.is_file():
        manifest = load_json(manifest_path)
        expected = int(manifest.get("dataset_count", 0))
        for item in manifest.get("samples", []):
            manifest_sample_names.append(item.get("sample") if isinstance(item, dict) else str(item))
        manifest_sample_names = [name for name in manifest_sample_names if name]
    rows = []
    if status_path.is_file():
        with status_path.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    allowed_samples = set(manifest_sample_names)
    sample_rows = [row for row in rows if row.get("sample") in allowed_samples]
    ignored_rows = [row for row in rows if row.get("sample") not in allowed_samples]
    latest = {}
    for row in sample_rows:
        latest[row["sample"]] = row
    counts = Counter(row["status"] for row in latest.values())
    aggregate_path = production_root / "verification_summary.json"
    aggregate = load_json(aggregate_path) if aggregate_path.is_file() else None
    failure_path = production_root / "_FAILED"
    failure = None
    if failure_path.is_file():
        try:
            failure = load_json(failure_path)
        except json.JSONDecodeError:
            failure = {"status": "FAILED", "reason": failure_path.read_text(encoding="utf-8").strip()}
    aggregate_all_pass = aggregate.get("all_pass") if isinstance(aggregate, dict) else None
    aggregate_complete = bool(
        isinstance(aggregate, dict)
        and aggregate_all_pass is True
        and int(aggregate.get("dataset_count", -1)) == expected
        and int(aggregate.get("n_pass", -1)) == expected
    )
    sample_status_complete = bool(
        expected > 0
        and set(latest) == allowed_samples
        and len(latest) == expected
        and all(row.get("status") == "PASS" for row in latest.values())
    )
    if failure is not None or counts.get("FAIL", 0) or aggregate_all_pass is False:
        state = "failed"
    elif aggregate_complete and sample_status_complete:
        state = "complete"
    elif rows or aggregate is not None:
        state = "in_progress"
    else:
        state = "not_started"
    return {
        "root": str(production_root),
        "expected_datasets": expected,
        "status_rows": rows,
        "sample_status_rows": sample_rows,
        "ignored_non_sample_status_rows": ignored_rows,
        "manifest_sample_names": manifest_sample_names,
        "latest_by_sample": latest,
        "n_started_or_latest": len(latest),
        "n_pass": int(counts.get("PASS", 0)),
        "n_fail": int(counts.get("FAIL", 0)),
        "status_counts": dict(counts),
        "aggregate_verification": aggregate,
        "aggregate_exists": aggregate is not None,
        "aggregate_all_pass": aggregate_all_pass,
        "failure_marker": str(failure_path) if failure_path.is_file() else None,
        "failure": failure,
        "aggregate_complete": aggregate_complete,
        "sample_status_complete": sample_status_complete,
        "state": state,
    }


class HtmlBuilder:
    def __init__(self):
        self.tooltip_index = 0
        self.table_index = 0

    def tooltip(self, shown, source_key="derived"):
        """Keep values visually clean while retaining machine-readable provenance."""
        self.tooltip_index += 1
        return (
            f'<span class="sourced-value" data-source-key="{html.escape(source_key)}">'
            f'{html.escape(str(shown))}</span>'
        )

    def number(self, value, kind="integer", source_key="derived"):
        return self.tooltip(format_value(value, kind), source_key)

    def chart_source(self, tooltip_id, source_keys, chart_title):
        return (
            f'<a class="source-link chart-source" href="{html.escape(SOURCE_NOTES_NAME)}" '
            f'aria-label="開啟「{html.escape(chart_title)}」的資料來源紀錄">來源</a>'
        )

    def exact_table(self, headers, rows, source_keys, numeric_columns=None, first_col_left=True):
        self.table_index += 1
        numeric_columns = set(range(1, len(headers))) if numeric_columns is None else set(numeric_columns)
        caption = "資料表 " + str(self.table_index) + "：" + "、".join(str(header) for header in headers[:3])
        thead = "".join(f'<th scope="col">{html.escape(str(header))}</th>' for header in headers)
        body_rows = []
        for row_index, row in enumerate(rows):
            cells = []
            for index, value in enumerate(row):
                shown = html.escape(str(value))
                if index in numeric_columns:
                    if isinstance(source_keys, list) and source_keys and isinstance(source_keys[0], list):
                        source_key = source_keys[row_index][index]
                    elif isinstance(source_keys, list):
                        source_key = source_keys[index]
                    else:
                        source_key = source_keys
                    shown = self.tooltip(value, source_key)
                if first_col_left and index == 0:
                    cells.append(f'<th scope="row" class="left">{shown}</th>')
                else:
                    cells.append(f"<td>{shown}</td>")
            body_rows.append("<tr>" + "".join(cells) + "</tr>")
        return (
            '<div class="table-scroll" tabindex="0" role="region" aria-label="' + html.escape(caption) + '"><table><caption>' + html.escape(caption) +
            '</caption><thead><tr>' + thead + "</tr></thead><tbody>" +
            "".join(body_rows) + "</tbody></table></div>"
        )


def chart_payload(samples):
    charts = []

    funnel_data = []
    for sample in samples:
        funnel = sample["site_funnel"]
        for stage, field in [
            ("chrX/Y/其他", "out_of_scope"),
            ("singleton", "singleton"),
            ("densest-8 排除", "cap_excluded"),
            ("read 不足", "read_unsupported"),
            ("retained", "retained"),
        ]:
            funnel_data.append({"dataset": sample["sample"], "stage": stage, "value": funnel[field]})
    charts.append({
        "id": "site-funnel-composition",
        "height": 380,
        "type": "bar",
        "dataset": {
            "id": "site-funnel-composition",
            "title": "PASS sSNV universe composition",
            "data": funnel_data,
            "chart_spec": {
                "id": "site-funnel-composition",
                "dataset": "site-funnel-composition",
                "title": "PASS sSNV universe composition",
                "type": "bar",
                "palette": {"kind": "categorical", "roots": ["blue", "gold", "orange", "olive", "pink"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "value", "type": "quantitative", "label": "sSNV 數"},
                    "color": {"field": "stage", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "sSNV 數",
                "valueFormat": "number",
                "settings": {"groupMode": "stacked"},
            },
        },
    })

    region_stage_data = []
    for sample in samples:
        region = sample["regions"]
        for label, field in [
            ("W_pre", "W_pre"),
            ("W_ret", "W_ret"),
            ("W_tree", "W_tree"),
            ("W_primary", "W_primary"),
        ]:
            region_stage_data.append({"dataset": sample["sample"], "stage": label, "value": region[field]})
    charts.append({
        "id": "region-stage-counts",
        "height": 380,
        "type": "bar",
        "dataset": {
            "id": "region-stage-counts",
            "title": "Region-stage counts",
            "data": region_stage_data,
            "chart_spec": {
                "id": "region-stage-counts",
                "dataset": "region-stage-counts",
                "title": "Region-stage counts",
                "type": "bar",
                "palette": {"kind": "categorical", "roots": ["blue", "gold", "orange", "olive"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "value", "type": "quantitative", "label": "region/group 數"},
                    "color": {"field": "stage", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "region/group 數",
                "valueFormat": "number",
                "settings": {"groupMode": "grouped"},
            },
        },
    })

    k_data = []
    for sample in samples:
        total = sample["regions"]["W_ret"]
        for k in range(2, 9):
            k_data.append({
                "dataset": sample["sample"],
                "k": f"k={k}",
                "share": sample["k_distribution"]["counts"][str(k)] / total if total else 0,
            })
    charts.append({
        "id": "k-distribution",
        "height": 380,
        "type": "bar",
        "dataset": {
            "id": "k-distribution",
            "title": "Solver input k distribution",
            "data": k_data,
            "chart_spec": {
                "id": "k-distribution",
                "dataset": "k-distribution",
                "title": "Solver input k distribution",
                "type": "bar",
                "palette": {"kind": "sequential", "roots": ["blue"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "區域比例"},
                    "color": {"field": "k", "type": "ordinal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "區域比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "stacked"},
            },
        },
    })

    read_data = []
    for sample in samples:
        read = sample["read_evidence"]
        denominator = read["retained_pattern_observations"]
        for label, field in [("完整 genotype pattern", "full_pattern_observations"), ("partial-X pattern", "partial_pattern_observations")]:
            read_data.append({
                "dataset": sample["sample"],
                "pattern": label,
                "share": read[field] / denominator if denominator else 0,
            })
    charts.append({
        "id": "read-pattern-composition",
        "height": 360,
        "type": "bar",
        "dataset": {
            "id": "read-pattern-composition",
            "title": "Full versus partial read-pattern composition",
            "data": read_data,
            "chart_spec": {
                "id": "read-pattern-composition",
                "dataset": "read-pattern-composition",
                "title": "Full versus partial read-pattern composition",
                "type": "bar",
                "palette": {"kind": "semantic", "roots": ["blue", "gold"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "pattern observation 比例"},
                    "color": {"field": "pattern", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "pattern observation 比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "stacked"},
            },
        },
    })

    hp_data = []
    for sample in samples:
        region = sample["regions"]
        denominator = region["W_tree"]
        for label, field in [("單一 primary HP", "single_primary_HP"), ("HP1+HP2", "multi_HP"), ("無 primary", "no_primary")]:
            hp_data.append({
                "dataset": sample["sample"],
                "class": label,
                "share": region[field] / denominator if denominator else 0,
            })
    charts.append({
        "id": "hp-multiplicity",
        "height": 360,
        "type": "bar",
        "dataset": {
            "id": "hp-multiplicity",
            "title": "Primary HP multiplicity",
            "data": hp_data,
            "chart_spec": {
                "id": "hp-multiplicity",
                "dataset": "hp-multiplicity",
                "title": "Primary HP multiplicity",
                "type": "bar",
                "palette": {"kind": "categorical", "roots": ["blue", "gold", "olive"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "W_tree 比例"},
                    "color": {"field": "class", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "W_tree 比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "stacked"},
            },
        },
    })

    outcome_data = []
    for sample in samples:
        units = sample["units"]
        denominator = units["primary"]
        for label, field in [
            ("模型內唯一", "determined"),
            ("多候選／歧義", "ambiguous"),
            ("solver capped", "solver_capped"),
            ("需要 recurrence", "recurrence"),
        ]:
            outcome_data.append({
                "dataset": sample["sample"],
                "class": label,
                "share": units[field] / denominator if denominator else 0,
            })
    charts.append({
        "id": "tree-outcomes",
        "height": 380,
        "type": "bar",
        "dataset": {
            "id": "tree-outcomes",
            "title": "Primary reconstruction-unit outcomes",
            "data": outcome_data,
            "chart_spec": {
                "id": "tree-outcomes",
                "dataset": "tree-outcomes",
                "title": "Primary reconstruction-unit outcomes",
                "type": "bar",
                "palette": {"kind": "semantic", "roots": ["blue", "gold", "orange", "pink"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "primary unit 比例"},
                    "color": {"field": "class", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "primary unit 比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "stacked"},
            },
        },
    })

    determinacy_data = []
    for sample in samples:
        units = sample["units"]
        levels = [
            ("Exact tree unique", pct(units["exact_tree_unique"], units["noncapped_complete"])),
            ("Shape unique", pct(units["shape_unique"], units["noncapped_complete"])),
        ]
        for label, value in levels:
            determinacy_data.append({"dataset": sample["sample"], "level": label, "share": value or 0})
    charts.append({
        "id": "determinacy-levels",
        "height": 380,
        "type": "bar",
        "dataset": {
            "id": "determinacy-levels",
            "title": "Exact-tree versus shape uniqueness",
            "data": determinacy_data,
            "chart_spec": {
                "id": "determinacy-levels",
                "dataset": "determinacy-levels",
                "title": "Exact-tree versus shape uniqueness",
                "type": "bar",
                "palette": {"kind": "semantic", "roots": ["blue", "gold"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "比例"},
                    "color": {"field": "level", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "grouped"},
            },
        },
    })

    cap_data = [{
        "dataset": sample["sample"],
        "share": pct(sample["units"]["solver_capped"], sample["units"]["primary"]) or 0,
    } for sample in samples]
    charts.append({
        "id": "solver-cap-rate",
        "height": 380,
        "type": "bar",
        "dataset": {
            "id": "solver-cap-rate",
            "title": "Solver-capped rate",
            "data": cap_data,
            "chart_spec": {
                "id": "solver-cap-rate",
                "dataset": "solver-cap-rate",
                "title": "Solver-capped rate",
                "type": "bar",
                "palette": {"kind": "sequential", "roots": ["orange"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "primary unit 比例"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "primary unit 比例",
                "valueFormat": "percent",
                "settings": {"orientation": "horizontal", "groupMode": "grouped"},
            },
        },
    })

    chromosome_data = [{
        "dataset": sample["sample"],
        "share": sample["sex_chromosome_census"]["nonautosomal_fraction"] or 0,
    } for sample in samples]
    charts.append({
        "id": "sex-chromosome-share",
        "height": 380,
        "type": "bar",
        "dataset": {
            "id": "sex-chromosome-share",
            "title": "Non-autosomal PASS sSNV share",
            "data": chromosome_data,
            "chart_spec": {
                "id": "sex-chromosome-share",
                "dataset": "sex-chromosome-share",
                "title": "Non-autosomal PASS sSNV share",
                "type": "bar",
                "palette": {"kind": "sequential", "roots": ["gold"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "全部 PASS sSNV 比例"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "全部 PASS sSNV 比例",
                "valueFormat": "percent",
                "settings": {"orientation": "horizontal", "groupMode": "grouped"},
            },
        },
    })

    hp_h3_data = []
    hp_h3_labels = [
        ("single_without_h3", "單一 primary · 無 H3?"),
        ("single_with_h3", "單一 primary · 有 H3?"),
        ("double_without_h3", "HP1+HP2 · 無 H3?"),
        ("double_with_h3", "HP1+HP2 · 有 H3?"),
        ("none_without_h3", "無 primary · 無 H3?"),
        ("none_with_h3", "無 primary · 有 H3?"),
    ]
    for sample in samples:
        counts = sample["region_candidate_topology"]["hp_h3"]
        denominator = sample["regions"]["W_tree"]
        for key, label in hp_h3_labels:
            hp_h3_data.append({
                "dataset": sample["sample"],
                "class": label,
                "share": counts[key] / denominator if denominator else 0,
            })
    charts.append({
        "id": "hp-h3-strata",
        "height": 390,
        "type": "bar",
        "dataset": {
            "id": "hp-h3-strata",
            "title": "Primary HP multiplicity × H3? auxiliary",
            "data": hp_h3_data,
            "chart_spec": {
                "id": "hp-h3-strata",
                "dataset": "hp-h3-strata",
                "title": "Primary HP multiplicity × H3? auxiliary",
                "type": "bar",
                "palette": {"kind": "categorical", "roots": ["blue", "gold", "orange", "olive", "pink"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "W_tree 比例"},
                    "color": {"field": "class", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "W_tree 比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "stacked"},
            },
        },
    })

    c_labels = ["C=1", "C=2", "C=3", "C=4", "C=5", "C=6", "C>6", "incomplete"]
    for scope, chart_id, title in [
        ("single", "candidate-combinations-single", "單一 primary HP 的 joint exact-tree combinations"),
        ("double", "candidate-combinations-double", "HP1+HP2 的 joint exact-tree combinations"),
    ]:
        c_data = []
        for sample in samples:
            counts = sample["region_candidate_topology"]["C_bins"][scope]
            denominator = sum(counts.values())
            for label in c_labels:
                c_data.append({
                    "dataset": sample["sample"],
                    "C_bin": label,
                    "share": counts[label] / denominator if denominator else 0,
                })
        charts.append({
            "id": chart_id,
            "height": 390,
            "type": "bar",
            "dataset": {
                "id": chart_id,
                "title": title,
                "data": c_data,
                "chart_spec": {
                    "id": chart_id,
                    "dataset": chart_id,
                    "title": title,
                    "type": "bar",
                    "palette": {"kind": "sequential", "roots": ["blue"]},
                    "encodings": {
                        "x": {"field": "dataset", "type": "nominal"},
                        "y": {"field": "share", "type": "quantitative", "label": f"{scope} primary region 比例"},
                        "color": {"field": "C_bin", "type": "ordinal"},
                    },
                    "xAxisTitle": "",
                    "yAxisTitle": f"{scope} primary region 比例",
                    "valueFormat": "percent",
                    "settings": {"groupMode": "stacked"},
                },
            },
        })

    topology_labels = [
        ("exact_and_topology_unique", "C=1 · Topo=1"),
        ("topology_unique_exact_multiple", "C>1 · Topo=1"),
        ("topology_multiple_exact_multiple", "C>1 · Topo>1"),
        ("incomplete", "incomplete / capped"),
    ]
    topology_data = []
    hidden_data = []
    for sample in samples:
        metrics = sample["region_candidate_topology"]
        denominator = sample["regions"]["W_primary"]
        for key, label in topology_labels:
            topology_data.append({
                "dataset": sample["sample"],
                "class": label,
                "share": metrics["topology_classes"][key] / denominator if denominator else 0,
            })
        for key, label in [
            ("hidden_zero", "0 hidden nodes"),
            ("hidden_positive", ">0 hidden nodes"),
            ("incomplete", "incomplete / unavailable"),
        ]:
            hidden_data.append({
                "dataset": sample["sample"],
                "class": label,
                "share": metrics["hidden_classes"][key] / denominator if denominator else 0,
            })
    charts.append({
        "id": "region-topology-classes",
        "height": 390,
        "type": "bar",
        "dataset": {
            "id": "region-topology-classes",
            "title": "Region-level exact/topology uniqueness classes",
            "data": topology_data,
            "chart_spec": {
                "id": "region-topology-classes",
                "dataset": "region-topology-classes",
                "title": "Region-level exact/topology uniqueness classes",
                "type": "bar",
                "palette": {"kind": "semantic", "roots": ["blue", "gold", "orange", "pink"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "W_primary 比例"},
                    "color": {"field": "class", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "W_primary 比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "stacked"},
            },
        },
    })
    charts.append({
        "id": "hidden-node-classes",
        "height": 370,
        "type": "bar",
        "dataset": {
            "id": "hidden-node-classes",
            "title": "Region-level hidden mutation-state node status",
            "data": hidden_data,
            "chart_spec": {
                "id": "hidden-node-classes",
                "dataset": "hidden-node-classes",
                "title": "Region-level hidden mutation-state node status",
                "type": "bar",
                "palette": {"kind": "semantic", "roots": ["blue", "gold", "orange"]},
                "encodings": {
                    "x": {"field": "dataset", "type": "nominal"},
                    "y": {"field": "share", "type": "quantitative", "label": "W_primary 比例"},
                    "color": {"field": "class", "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "W_primary 比例",
                "valueFormat": "percent",
                "settings": {"groupMode": "stacked"},
            },
        },
    })
    absolute_composition_ids = {"site-funnel-composition"}
    normalized_composition_ids = {
        "k-distribution", "read-pattern-composition", "hp-multiplicity", "tree-outcomes",
        "hp-h3-strata", "candidate-combinations-single", "candidate-combinations-double",
        "region-topology-classes", "hidden-node-classes",
    }
    horizontal_comparison_ids = {"solver-cap-rate", "sex-chromosome-share"}
    for chart in charts:
        spec = chart["dataset"]["chart_spec"]
        if chart["id"] in absolute_composition_ids:
            chart["type"] = spec["type"] = "stackedBar"
            spec.setdefault("settings", {})["groupMode"] = "stacked"
        elif chart["id"] in normalized_composition_ids:
            chart["type"] = spec["type"] = "stackedBar100"
            spec.setdefault("settings", {})["groupMode"] = "stacked100"
        elif chart["id"] in horizontal_comparison_ids:
            chart["type"] = spec["type"] = "horizontalBar"
    chart_questions = {
        "site-funnel-composition": ("composition", "Where do historical PASS sSNVs leave the analysis scope?", "Mutually exclusive site-level branches require a stacked composition view."),
        "region-stage-counts": ("comparison", "How many region/groups remain at each W stage?", "Grouped bars preserve the fact that these are nested counts, not additive branches."),
        "k-distribution": ("composition", "What k characterizes retained MLHP groups?", "A 100% composition view compares group-size profiles across datasets."),
        "read-pattern-composition": ("composition", "How much observation multiplicity is full versus partial-X?", "A common observation denominator supports a stacked composition view."),
        "hp-multiplicity": ("composition", "How many primary reconstruction units occur per emitted region?", "Mutually exclusive W_tree classes sum to one."),
        "tree-outcomes": ("composition", "What is the L1 model outcome per primary reconstruction unit?", "Mutually exclusive L1 classes share one denominator."),
        "determinacy-levels": ("comparison", "How does exact-tree uniqueness differ from shape uniqueness?", "Both measures share the non-capped analysis-complete primary-unit denominator."),
        "solver-cap-rate": ("comparison", "Which dataset has the largest incomplete-search burden?", "A single-rate horizontal comparison exposes the outlier."),
        "sex-chromosome-share": ("comparison", "How much of the historical caller universe is non-autosomal?", "A single-rate horizontal comparison shows why chrX/Y need a separate contract."),
        "hp-h3-strata": ("composition", "How do primary HP multiplicity and H3? auxiliary evidence intersect?", "Six mutually exclusive W_tree strata expose H3? without adding it to the primary denominator."),
        "candidate-combinations-single": ("composition", "How many exact candidates remain when one primary HP is mutation-bearing?", "C bins plus incomplete cases partition the one-primary region denominator."),
        "candidate-combinations-double": ("composition", "How many joint exact combinations remain when HP1 and HP2 are both mutation-bearing?", "The Cartesian-product C bins plus incomplete cases partition the two-primary region denominator."),
        "region-topology-classes": ("composition", "Which exact/tree-shape uniqueness state remains per primary region?", "Three feasible classes plus incomplete search partition W_primary."),
        "hidden-node-classes": ("composition", "How often do complete regional candidates require unobserved intermediate mutation states?", "Zero, positive and unavailable hidden-node states partition W_primary without calling them clones."),
    }
    for chart in charts:
        intent, question, rationale = chart_questions[chart["id"]]
        chart["dataset"]["chart_spec"].update({
            "intent": intent,
            "question": question,
            "rationale": rationale,
        })
    return {"charts": charts}


def validate_official_chart_contract(payload, contract_path: Path):
    """Run the Data Analytics plugin's canonical chart contract fail-closed."""
    result = {"path": str(contract_path), "errors": [], "warnings": []}
    if not contract_path.is_file():
        result["errors"].append(f"official chart contract missing: {contract_path}")
        return result
    module_name = "intersubmod_data_analytics_chart_contract"
    module_spec = importlib.util.spec_from_file_location(module_name, contract_path)
    if module_spec is None or module_spec.loader is None:
        result["errors"].append(f"official chart contract cannot be loaded: {contract_path}")
        return result
    module = importlib.util.module_from_spec(module_spec)
    sys.path.insert(0, str(contract_path.parent))
    try:
        module_spec.loader.exec_module(module)
    except Exception as exc:  # fail closed on contract import/runtime drift
        result["errors"].append(f"official chart contract import failed: {exc}")
        return result
    finally:
        sys.path.pop(0)

    sibling_specs = [chart["dataset"]["chart_spec"] for chart in payload["charts"]]
    for index, chart in enumerate(payload["charts"]):
        path = f"charts[{index}]"
        spec = chart["dataset"]["chart_spec"]
        rows = chart["dataset"]["data"]
        try:
            chart_type = module.validate_chart_type(spec.get("type"), path)
            if chart.get("type") != chart_type:
                raise module.ContractError(f"{path} top-level type differs from chart_spec.type")
            module.validate_chart_intent_metadata(spec, path, require=True)
            module.validate_chart_presentation_options(spec, path)
            module.validate_chart_field_compatibility(
                chart_type,
                module.chart_encoding_field(spec, "x"),
                [field for field in module.chart_y_fields(spec) if field],
                path,
                spec,
            )
            module.validate_chart_data_compatibility(spec, rows, path)
            for warning in module.chart_quality_warnings(spec, rows, sibling_charts=sibling_specs):
                result["warnings"].append(f"{chart['id']}: {warning}")
        except Exception as exc:
            result["errors"].append(f"{chart['id']}: {exc}")
    return result


def chart_block(builder, chart_id, title, subtitle, fallback_table, exact_table, source_keys, note):
    # The same semantic table appears once as a no-JS fallback and once as an
    # expandable exact-data table.  Give the second copy unique tooltip IDs.
    exact_table = exact_table.replace("source-tooltip-", f"source-tooltip-{chart_id}-exact-")
    return f'''
    <figure class="card source-figure report-chart" aria-labelledby="chart-title-{html.escape(chart_id)}">
      <div class="card-head"><h3 id="chart-title-{html.escape(chart_id)}">{html.escape(title)}</h3><p>{html.escape(subtitle)}</p><span class="historical-chip">Historical snapshot · not current canonical result</span></div>
      <p class="chart-scroll-hint">左右滑動查看全部 7 datasets</p>
      <div class="chart-wrap" role="region" tabindex="0" aria-label="{html.escape(title)}圖表；可左右捲動">
        <div data-recharts-chart="{html.escape(chart_id)}">
          <div class="chart-fallback" data-recharts-fallback>{fallback_table}</div>
          <div data-recharts-live aria-hidden="true"></div>
        </div>
      </div>
      <figcaption class="chart-note">{html.escape(note)}</figcaption>
      {builder.chart_source(f"chart-source-{chart_id}", source_keys, title)}
    </figure>
    <details class="exact-data"><summary>展開「{html.escape(title)}」的精確數據、分子與分母</summary>{exact_table}</details>
    '''


def build_report_html(snapshot, shell_template: Path):
    builder = HtmlBuilder()
    samples = snapshot["samples"]
    by_name = {sample["sample"]: sample for sample in samples}
    hcc = by_name["HCC1395"]
    h2009 = by_name["H2009"]
    verification = snapshot["baseline_verification"]
    production = snapshot["clean_production"]
    auxiliary = snapshot["auxiliary_evidence"]
    method = snapshot["method_contract"]
    case = hcc["well_explained_case"]
    generated = snapshot["generated_at"]
    production_failure = production.get("failure") or {}
    production_failure_label = production_failure.get("reason_code") or production_failure.get("status") or "none"
    production_latest = production.get("latest_by_sample", {})
    production_passed = [
        sample for sample in SAMPLE_ORDER
        if production_latest.get(sample, {}).get("status") == "PASS"
    ]
    production_active = [
        sample for sample in SAMPLE_ORDER
        if production_latest.get(sample, {}).get("status") in {"START", "RUNNING"}
    ]
    production_failed = [
        sample for sample in SAMPLE_ORDER
        if production_latest.get(sample, {}).get("status") == "FAIL"
    ]
    production_root_name = Path(production["root"]).name

    if production["state"] == "complete":
        production_readout = "normalized raw-all production 已完成 7/7，aggregate verification 通過"
        production_boundary = "下游 canonical tree 結果仍須由明示的 completed layered run 另行發布。"
    elif production["state"] == "in_progress":
        active_label = "、".join(production_active) if production_active else "等待下一個 dataset"
        production_readout = (
            f"normalized raw-all production 目前 {production['n_pass']}/{production['expected_datasets']} PASS；"
            f"進行中：{active_label}"
        )
        production_boundary = "aggregate _SUCCESS 尚不存在；本頁歷史圖表不得替代尚未完成的 canonical 結果。"
    elif production["state"] == "failed":
        production_readout = f"production 已 fail-closed：{production_failure_label}"
        production_boundary = "失敗 run 不產生 canonical 結果，也不得由 historical snapshot 補位。"
    else:
        production_readout = "normalized raw-all production 尚未開始"
        production_boundary = "目前只有 historical engineering baseline，可用於方法觀察而非 current result。"

    production_passed_label = "、".join(production_passed) if production_passed else "尚無"
    production_failed_label = "、".join(production_failed) if production_failed else "無"
    generated_display = generated.replace("T", " ")[:16] + " · Asia/Taipei"

    funnel_rows = []
    for sample in samples:
        f = sample["site_funnel"]
        funnel_rows.append([
            sample["sample"], f"{f['universe']:,}", f"{f['out_of_scope']:,}", f"{f['singleton']:,}",
            f"{f['cap_excluded']:,}", f"{f['read_unsupported']:,}", f"{f['retained']:,}",
        ])
    funnel_table = builder.exact_table(
        ["Dataset", "PASS U", "non-autosomal", "singleton", "cap-excluded", "read-unsupported", "retained"],
        funnel_rows, "verification")

    region_rows = []
    for sample in samples:
        region = sample["regions"]
        region_rows.append([
            sample["sample"], f"{region['W_pre']:,}", f"{region['W_ret']:,}", f"{region['W_tree']:,}",
            f"{region['W_primary']:,}", f"{region['W_ret'] - region['W_tree']:,}",
        ])
    region_table = builder.exact_table(
        ["Dataset", "W_pre", "W_ret", "W_tree", "W_primary", "W_ret−W_tree"],
        region_rows, "derived")

    k_rows = []
    for sample in samples:
        total = sample["regions"]["W_ret"]
        k_rows.append([
            sample["sample"],
            *[f"{sample['k_distribution']['counts'][str(k)]:,} ({format_value(sample['k_distribution']['counts'][str(k)] / total, 'percent1')})" for k in range(2, 9)],
            f"{sample['k_distribution']['natural_k8']:,}",
            f"{sample['k_distribution']['compressed_k8']:,}",
            format_value(sample['k_distribution']['mean_k'], "decimal2"),
            format_value(sample['k_distribution']['median_k'], "decimal1"),
        ])
    k_table = builder.exact_table(
        ["Dataset", *[f"k={k}" for k in range(2, 9)], "natural k8", "cap-compressed k8", "mean k", "median k"],
        k_rows, "mlhp")

    read_rows = []
    for sample in samples:
        read = sample["read_evidence"]
        denom = read["retained_pattern_observations"]
        read_rows.append([
            sample["sample"], f"{read['full_pattern_observations']:,}",
            format_value(pct(read["full_pattern_observations"], denom), "percent"),
            f"{read['partial_pattern_observations']:,}",
            format_value(pct(read["partial_pattern_observations"], denom), "percent"),
            f"{denom:,}",
            f"{read['regions_with_full']:,}/{sample['regions']['W_ret']:,}",
            format_value(pct(read['regions_with_full'], sample['regions']['W_ret']), "percent"),
            f"{read['regions_with_alt_population']:,}/{sample['regions']['W_ret']:,}",
            format_value(pct(read['regions_with_alt_population'], sample['regions']['W_ret']), "percent"),
        ])
    read_table = builder.exact_table(
        ["Dataset", "full patterns", "full %", "partial-X patterns", "partial-X %", "pattern denominator",
         "regions with full/W_ret", "full-region %", "regions with ALT full pop/W_ret", "ALT-region %"],
        read_rows, "mlhp")

    hp_rows = []
    for sample in samples:
        region = sample["regions"]
        hp_rows.append([
            sample["sample"], f"{region['single_primary_HP']:,}", f"{region['multi_HP']:,}",
            format_value(pct(region["multi_HP"], region["W_tree"]), "percent"),
            f"{region['no_primary']:,}", f"{region['W_tree']:,}",
        ])
    hp_table = builder.exact_table(
        ["Dataset", "one primary reconstruction unit", "HP1 & HP2 mutation-bearing units", "two-unit %", "no primary unit", "W_tree"], hp_rows, "region")

    outcome_rows = []
    for sample in samples:
        units = sample["units"]
        outcome_rows.append([
            sample["sample"], f"{units['determined']:,} ({format_value(pct(units['determined'], units['primary']), 'percent')})",
            f"{units['ambiguous']:,} ({format_value(pct(units['ambiguous'], units['primary']), 'percent')})",
            f"{units['solver_capped']:,} ({format_value(pct(units['solver_capped'], units['primary']), 'percent')})",
            f"{units['recurrence']:,} ({format_value(pct(units['recurrence'], units['primary']), 'percent')})",
            f"{units['primary']:,}",
        ])
    outcome_table = builder.exact_table(
        ["Dataset", "模型內候選唯一", "多候選／歧義", "solver capped", "model recurrence-required", "primary reconstruction units"],
        outcome_rows, "layered")

    determinacy_rows = []
    for sample in samples:
        units = sample["units"]
        determinacy_rows.append([
            sample["sample"],
            f"{units['exact_tree_unique']:,}/{units['noncapped_complete']:,} ({format_value(pct(units['exact_tree_unique'], units['noncapped_complete']), 'percent')})",
            f"{units['shape_unique']:,}/{units['noncapped_complete']:,} ({format_value(pct(units['shape_unique'], units['noncapped_complete']), 'percent')})",
            f"{units['region_all_exact']:,}/{sample['regions']['W_primary']:,} ({format_value(pct(units['region_all_exact'], sample['regions']['W_primary']), 'percent')})",
            f"{units['region_all_shape']:,}/{sample['regions']['W_primary']:,} ({format_value(pct(units['region_all_shape'], sample['regions']['W_primary']), 'percent')})",
        ])
    determinacy_table = builder.exact_table(
        ["Dataset", "exact-tree unique", "shape unique", "region all-exact", "region all-shape"], determinacy_rows, "derived")

    complexity_rows = []
    for sample in samples:
        units = sample["units"]
        complexity_rows.append([
            sample["sample"], f"{units['noncapped_complete']:,}", f"{units['total_candidate_trees']:,}",
            format_value(pct(units['total_candidate_trees'], units['noncapped_complete']), "decimal2"),
            format_value(units['candidate_tree_median'], "decimal1"),
            format_value(units['candidate_tree_p90'], "decimal1"),
            f"{units['hidden_total_complete']:,}", format_value(pct(units['hidden_total_complete'], units['noncapped_complete']), "decimal2"),
            "; ".join(f"{key}:{value:,}" for key, value in units["hidden_unit_bins_complete"].items()),
        ])
    complexity_table = builder.exact_table(
        ["Dataset", "analysis-complete non-capped units", "candidate trees", "mean trees/unit", "median", "P90",
         "complete-set hidden nodes", "mean nodes/complete unit", "unit bins 0/1/2/3/4/5+"], complexity_rows, "layered")

    l2_rows = []
    for sample in samples:
        l2 = sample["L2_cn"]
        split = "; ".join(f"{key}: {value:,}" for key, value in l2.get("primary_cn_split", {}).items()) or "none"
        l2_rows.append([
            sample["sample"], f"{l2['n_primary_recurrence_sent_to_cn']:,}", f"{sample['units']['primary']:,}",
            split, str(sample["cn_contract"]["cn_source"] or "unavailable"),
        ])
    l2_table = builder.exact_table(
        ["Dataset", "model recurrence-required sent to L2", "primary units", "post-hoc CN split", "CN source"],
        l2_rows, "layered", numeric_columns={1, 2})

    cap_rows = [[
        sample["sample"], f"{sample['units']['solver_capped']:,}", f"{sample['units']['primary']:,}",
        format_value(pct(sample["units"]["solver_capped"], sample["units"]["primary"]), "percent"),
    ] for sample in samples]
    cap_table = builder.exact_table(["Dataset", "solver-capped", "primary units", "rate"], cap_rows, "layered")

    chromosome_rows = [[
        sample["sample"], f"{sample['sex_chromosome_census']['chrX']:,}", f"{sample['sex_chromosome_census']['chrY']:,}",
        f"{sample['site_funnel']['universe']:,}", format_value(sample["sex_chromosome_census"]["nonautosomal_fraction"], "percent"),
    ] for sample in samples]
    chromosome_table = builder.exact_table(["Dataset", "chrX", "chrY", "PASS U", "non-autosomal %"], chromosome_rows, "vcf")

    def dual_funnel_rows(label, funnel, region):
        site_u = funnel["universe"]
        autosomal = funnel["autosomal"]
        multilocus_sites = funnel["multilocus_pre_cap"]
        wall = region["W_all_pre"]
        return [
            [label, "Site", "S / caller-PASS sSNV universe", "sSNV", f"{site_u:,}", f"{site_u:,} · 100.00%", f"{site_u:,} · 100.00%"],
            [label, "Site", "chr1–22 in scope", "sSNV", f"{autosomal:,}", f"{site_u:,} · {format_value(pct(autosomal, site_u), 'percent')}", f"{site_u:,} · {format_value(pct(autosomal, site_u), 'percent')}"],
            [label, "Site", "chrX/Y/other out of scope", "sSNV", f"{funnel['out_of_scope']:,}", f"{site_u:,} · {format_value(pct(funnel['out_of_scope'], site_u), 'percent')}", f"{site_u:,} · {format_value(pct(funnel['out_of_scope'], site_u), 'percent')}"],
            [label, "Site", "k=1 positional singleton", "sSNV = component", f"{funnel['singleton']:,}", f"{site_u:,} · {format_value(pct(funnel['singleton'], site_u), 'percent')}", f"{autosomal:,} · {format_value(pct(funnel['singleton'], autosomal), 'percent')}"],
            [label, "Site", "k>1 multi-locus pre-cap", "sSNV", f"{multilocus_sites:,}", f"{site_u:,} · {format_value(pct(multilocus_sites, site_u), 'percent')}", f"{autosomal:,} · {format_value(pct(multilocus_sites, autosomal), 'percent')}"],
            [label, "Site", "densest-8 cap-excluded", "sSNV", f"{funnel['cap_excluded']:,}", f"{site_u:,} · {format_value(pct(funnel['cap_excluded'], site_u), 'percent')}", f"{multilocus_sites:,} · {format_value(pct(funnel['cap_excluded'], multilocus_sites), 'percent')}"],
            [label, "Site", "read-unsupported", "sSNV", f"{funnel['read_unsupported']:,}", f"{site_u:,} · {format_value(pct(funnel['read_unsupported'], site_u), 'percent')}", f"{multilocus_sites:,} · {format_value(pct(funnel['read_unsupported'], multilocus_sites), 'percent')}"],
            [label, "Site", "retained", "sSNV", f"{funnel['retained']:,}", f"{site_u:,} · {format_value(pct(funnel['retained'], site_u), 'percent')}", f"{multilocus_sites:,} · {format_value(pct(funnel['retained'], multilocus_sites), 'percent')}"],
            [label, "Unit bridge", "W_all,pre = k1 components + k>1 groups", "region/group", f"{wall:,}", f"{wall:,} · 100.00%", f"{wall:,} · 100.00%"],
            [label, "Region", "k=1 components (not sent to tree solver)", "region", f"{region['W_k1']:,}", f"{wall:,} · {format_value(pct(region['W_k1'], wall), 'percent')}", f"{wall:,} · {format_value(pct(region['W_k1'], wall), 'percent')}"],
            [label, "Region", "k>1 W_pre", "group", f"{region['W_k_gt1']:,}", f"{wall:,} · {format_value(pct(region['W_k_gt1'], wall), 'percent')}", f"{wall:,} · {format_value(pct(region['W_k_gt1'], wall), 'percent')}"],
            [label, "Region", "W_ret", "group", f"{region['W_ret']:,}", f"{wall:,} · {format_value(pct(region['W_ret'], wall), 'percent')}", f"{region['W_pre']:,} · {format_value(pct(region['W_ret'], region['W_pre']), 'percent')}"],
            [label, "Region", "W_tree", "emitted region", f"{region['W_tree']:,}", f"{wall:,} · {format_value(pct(region['W_tree'], wall), 'percent')}", f"{region['W_ret']:,} · {format_value(pct(region['W_tree'], region['W_ret']), 'percent')}"],
            [label, "Region", "W_primary", "primary region", f"{region['W_primary']:,}", f"{wall:,} · {format_value(pct(region['W_primary'], wall), 'percent')}", f"{region['W_tree']:,} · {format_value(pct(region['W_primary'], region['W_tree']), 'percent')}"],
        ]

    pooled_funnel = {
        key: sum(sample["site_funnel"][key] for sample in samples)
        for key in ("universe", "autosomal", "out_of_scope", "singleton", "multilocus_pre_cap", "cap_excluded", "read_unsupported", "retained")
    }
    pooled_regions = {
        key: sum(sample["regions"][key] for sample in samples)
        for key in ("W_all_pre", "W_k1", "W_k_gt1", "W_pre", "W_ret", "W_tree", "W_primary")
    }
    pooled_funnel_table = builder.exact_table(
        ["Scope", "Layer", "Stage", "Unit", "Count", "Total denominator · share", "Parent denominator · relative share"],
        dual_funnel_rows("7-dataset row aggregate", pooled_funnel, pooled_regions),
        "derived",
        numeric_columns={4, 5, 6},
    )

    hp_h3_rows = []
    hp_h3_key_labels = [
        ("single_without_h3", "HP1 xor HP2", "not H3?"),
        ("single_with_h3", "HP1 xor HP2", "with H3?"),
        ("double_without_h3", "HP1 and HP2", "not H3?"),
        ("double_with_h3", "HP1 and HP2", "with H3?"),
        ("none_without_h3", "no primary", "not H3?"),
        ("none_with_h3", "no primary", "with H3?"),
    ]
    for sample in samples:
        metric = sample["region_candidate_topology"]
        counts = metric["hp_h3"]
        parent_totals = {
            "HP1 xor HP2": counts["single_without_h3"] + counts["single_with_h3"],
            "HP1 and HP2": counts["double_without_h3"] + counts["double_with_h3"],
            "no primary": counts["none_without_h3"] + counts["none_with_h3"],
        }
        for key, hp_class, h3_class in hp_h3_key_labels:
            count = counts[key]
            hp_h3_rows.append([
                sample["sample"], hp_class, h3_class, f"{count:,}",
                format_value(pct(count, sample["regions"]["W_tree"]), "percent"),
                format_value(pct(count, parent_totals[hp_class]), "percent"),
                f"W_tree={sample['regions']['W_tree']:,}", f"parent={parent_totals[hp_class]:,}",
            ])
    pooled_hp_h3 = Counter()
    for sample in samples:
        pooled_hp_h3.update(sample["region_candidate_topology"]["hp_h3"])
    pooled_wtree = sum(sample["regions"]["W_tree"] for sample in samples)
    pooled_parent_totals = {
        "HP1 xor HP2": pooled_hp_h3["single_without_h3"] + pooled_hp_h3["single_with_h3"],
        "HP1 and HP2": pooled_hp_h3["double_without_h3"] + pooled_hp_h3["double_with_h3"],
        "no primary": pooled_hp_h3["none_without_h3"] + pooled_hp_h3["none_with_h3"],
    }
    for key, hp_class, h3_class in hp_h3_key_labels:
        count = pooled_hp_h3[key]
        hp_h3_rows.append([
            "SUM 7 dataset rows", hp_class, h3_class, f"{count:,}",
            format_value(pct(count, pooled_wtree), "percent"),
            format_value(pct(count, pooled_parent_totals[hp_class]), "percent"),
            f"W_tree={pooled_wtree:,}", f"parent={pooled_parent_totals[hp_class]:,}",
        ])
    hp_h3_table = builder.exact_table(
        ["Dataset", "Primary class", "H3?", "Count", "Total share / W_tree", "Relative share / HP class", "Total denominator", "Parent denominator"],
        hp_h3_rows, "region", numeric_columns={3, 4, 5, 6, 7})

    c_tables = {}
    for scope, scope_label in [("single", "one primary HP"), ("double", "HP1 and HP2")]:
        rows = []
        for sample in samples:
            metrics = sample["region_candidate_topology"]
            counts = metrics["C_bins"][scope]
            parent = sum(counts.values())
            complete_parent = parent - counts["incomplete"]
            for bucket in ("C=1", "C=2", "C=3", "C=4", "C=5", "C=6", "C>6", "incomplete"):
                count = counts[bucket]
                complete_share = "不適用" if bucket == "incomplete" else format_value(pct(count, complete_parent), "percent")
                rows.append([
                    sample["sample"], scope_label, bucket, f"{count:,}",
                    format_value(pct(count, sample["regions"]["W_tree"]), "percent"),
                    format_value(pct(count, parent), "percent"), complete_share,
                    f"W_tree={sample['regions']['W_tree']:,}", f"HP parent={parent:,}", f"complete={complete_parent:,}",
                ])
        pooled_counts = Counter()
        for sample in samples:
            pooled_counts.update(sample["region_candidate_topology"]["C_bins"][scope])
        pooled_parent = sum(pooled_counts.values())
        pooled_complete = pooled_parent - pooled_counts["incomplete"]
        for bucket in ("C=1", "C=2", "C=3", "C=4", "C=5", "C=6", "C>6", "incomplete"):
            count = pooled_counts[bucket]
            rows.append([
                "SUM 7 dataset rows", scope_label, bucket, f"{count:,}",
                format_value(pct(count, pooled_wtree), "percent"),
                format_value(pct(count, pooled_parent), "percent"),
                "不適用" if bucket == "incomplete" else format_value(pct(count, pooled_complete), "percent"),
                f"W_tree={pooled_wtree:,}", f"HP parent={pooled_parent:,}", f"complete={pooled_complete:,}",
            ])
        c_tables[scope] = builder.exact_table(
            ["Dataset", "HP scope", "C bin", "Count", "Total share / W_tree", "Relative share / HP scope", "Complete-only share", "Total denominator", "Parent denominator", "Complete denominator"],
            rows, "derived", numeric_columns={3, 4, 5, 6, 7, 8, 9})

    topology_rows = []
    topology_names = [
        ("exact_and_topology_unique", "C=1, Topo=1", "exact tree unique + topology unique"),
        ("topology_unique_exact_multiple", "C>1, Topo=1", "topology unique; exact tree non-unique"),
        ("topology_multiple_exact_multiple", "C>1, Topo>1", "topology and exact tree non-unique"),
        ("incomplete", "incomplete", "solver-capped/incomplete; C and Topo unavailable"),
    ]
    hidden_rows = []
    topology_hidden_rows = []
    topology_sources = list(samples)
    for sample in topology_sources:
        metrics = sample["region_candidate_topology"]
        wtree = sample["regions"]["W_tree"]
        wprimary = sample["regions"]["W_primary"]
        complete = metrics["complete_regions"]
        for scope, scope_label in [("all", "all primary regions"), ("single", "one primary HP"), ("double", "HP1 and HP2")]:
            counts = metrics["topology_classes"] if scope == "all" else metrics["topology_by_hp"][scope]
            parent = wprimary if scope == "all" else sum(metrics["C_bins"][scope].values())
            complete_parent = parent - counts["incomplete"]
            for key, formula, meaning in topology_names:
                count = counts[key]
                topology_rows.append([
                    sample["sample"], scope_label, formula, meaning, f"{count:,}",
                    format_value(pct(count, wtree), "percent"),
                    format_value(pct(count, parent), "percent"),
                    "不適用" if key == "incomplete" else format_value(pct(count, complete_parent), "percent"),
                    f"W_tree={wtree:,}", f"parent={parent:,}", f"complete={complete_parent:,}",
                ])
        for key, label in [("hidden_zero", "n_hidden,region=0"), ("hidden_positive", "n_hidden,region>0"), ("incomplete", "incomplete / unavailable")]:
            count = metrics["hidden_classes"][key]
            hidden_rows.append([
                sample["sample"], label, f"{count:,}",
                format_value(pct(count, wtree), "percent"),
                format_value(pct(count, wprimary), "percent"),
                "不適用" if key == "incomplete" else format_value(pct(count, complete), "percent"),
                f"W_tree={wtree:,}", f"W_primary={wprimary:,}", f"complete={complete:,}",
            ])
        for scope, scope_label in [("single", "one primary HP"), ("double", "HP1 and HP2")]:
            for topo_key, topo_formula, _ in topology_names[:3]:
                for hidden_key, hidden_label in [("hidden_zero", "H=0"), ("hidden_positive", "H>0")]:
                    count = metrics["hidden_by_topology_by_hp"][scope][topo_key][hidden_key]
                    topology_hidden_rows.append([
                        sample["sample"], scope_label, topo_formula, hidden_label, f"{count:,}",
                        format_value(pct(count, wtree), "percent"),
                        format_value(pct(count, metrics["complete_regions"]), "percent"),
                    ])
    pooled_wprimary = sum(sample["regions"]["W_primary"] for sample in samples)
    pooled_complete_regions = sum(sample["region_candidate_topology"]["complete_regions"] for sample in samples)
    for scope, scope_label in [("all", "all primary regions"), ("single", "one primary HP"), ("double", "HP1 and HP2")]:
        counts = Counter()
        for sample in samples:
            metric = sample["region_candidate_topology"]
            counts.update(metric["topology_classes"] if scope == "all" else metric["topology_by_hp"][scope])
        parent = pooled_wprimary if scope == "all" else sum(
            sum(sample["region_candidate_topology"]["C_bins"][scope].values()) for sample in samples
        )
        complete_parent = parent - counts["incomplete"]
        for key, formula, meaning in topology_names:
            count = counts[key]
            topology_rows.append([
                "SUM 7 dataset rows", scope_label, formula, meaning, f"{count:,}",
                format_value(pct(count, pooled_wtree), "percent"),
                format_value(pct(count, parent), "percent"),
                "不適用" if key == "incomplete" else format_value(pct(count, complete_parent), "percent"),
                f"W_tree={pooled_wtree:,}", f"parent={parent:,}", f"complete={complete_parent:,}",
            ])
    pooled_hidden = Counter()
    for sample in samples:
        pooled_hidden.update(sample["region_candidate_topology"]["hidden_classes"])
    for key, label in [("hidden_zero", "n_hidden,region=0"), ("hidden_positive", "n_hidden,region>0"), ("incomplete", "incomplete / unavailable")]:
        count = pooled_hidden[key]
        hidden_rows.append([
            "SUM 7 dataset rows", label, f"{count:,}",
            format_value(pct(count, pooled_wtree), "percent"),
            format_value(pct(count, pooled_wprimary), "percent"),
            "不適用" if key == "incomplete" else format_value(pct(count, pooled_complete_regions), "percent"),
            f"W_tree={pooled_wtree:,}", f"W_primary={pooled_wprimary:,}", f"complete={pooled_complete_regions:,}",
        ])
    for scope, scope_label in [("single", "one primary HP"), ("double", "HP1 and HP2")]:
        for topo_key, topo_formula, _ in topology_names[:3]:
            for hidden_key, hidden_label in [("hidden_zero", "H=0"), ("hidden_positive", "H>0")]:
                count = sum(
                    sample["region_candidate_topology"]["hidden_by_topology_by_hp"][scope][topo_key][hidden_key]
                    for sample in samples
                )
                topology_hidden_rows.append([
                    "SUM 7 dataset rows", scope_label, topo_formula, hidden_label, f"{count:,}",
                    format_value(pct(count, pooled_wtree), "percent"),
                    format_value(pct(count, pooled_complete_regions), "percent"),
                ])
    topology_table = builder.exact_table(
        ["Dataset", "HP scope", "Formal state", "Meaning", "Count", "Total share / W_tree", "Relative share / HP scope", "Complete-only share", "Total denominator", "Parent denominator", "Complete denominator"],
        topology_rows, "derived", numeric_columns={4, 5, 6, 7, 8, 9, 10})
    hidden_table = builder.exact_table(
        ["Dataset", "Hidden state", "Count", "Total share / W_tree", "Relative share / W_primary", "Complete-only share", "Total denominator", "Parent denominator", "Complete denominator"],
        hidden_rows, "derived", numeric_columns={2, 3, 4, 5, 6, 7, 8})
    topology_hidden_table = builder.exact_table(
        ["Dataset", "HP scope", "Topology class", "Hidden state", "Count", "Total share / W_tree", "Share / all complete regions"],
        topology_hidden_rows, "derived", numeric_columns={4, 5, 6})

    charts = "".join([
        chart_block(builder, "site-funnel-composition", "7 datasets 的 sSNV 六分支漏斗",
                    "Absolute sSNV counts；non-autosomal、singleton、densest-8 排除、read-unsupported、retained 為互斥分支。",
                    funnel_table, funnel_table, ["verification"],
                    "這是位點單位；不能與後面的 region W 直接接成同一條 Sankey。"),
        chart_block(builder, "region-stage-counts", "7 datasets 的 W_pre → W_ret → W_tree → W_primary",
                    "Absolute region/group counts；各階段為巢狀分母，不是相加後等於總數的互斥分支。",
                    region_table, region_table, ["verification", "region"],
                    "HCC1395 的 W_ret−W_tree=1；其餘六個 datasets 為 0。"),
        chart_block(builder, "k-distribution", "Retained MLHP group k 的組成",
                    "各 dataset 以 W_ret 為分母；k=2–8，並在精確表區分 natural 與 cap-compressed k8。",
                    k_table, k_table, ["mlhp"],
                    "k 是輸入 sSNV 數，不是樹節點數、clone 數或拓撲大小。"),
        chart_block(builder, "read-pattern-composition", "完整與 partial-X read-pattern 的比例",
                    "分母為 retained full-genotype patterns + 含 X 的 partial patterns；單位是 read×region observation。",
                    read_table, read_table, ["mlhp"],
                    "同一 read 可在不同 region 重複出現；這不是 unique molecules，更不是 cell count。"),
        chart_block(builder, "hp-multiplicity", "Region 的 primary reconstruction-unit multiplicity",
                    "分母 W_tree；HP1、HP2 都形成 mutation-bearing reconstruction unit 時列入 two-unit。",
                    hp_table, hp_table, ["region"],
                    "multi-HP 不是 multi-clone；HP 是等位 family，不是細胞族群標籤。"),
        chart_block(builder, "hp-h3-strata", "Primary HP × H3? 的完整六分支",
                    "分母 W_tree；先分 HP1 xor HP2、HP1 and HP2、no primary，再分 with/without H3?。",
                    hp_h3_table, hp_h3_table, ["region", "derived"],
                    "H3? 是 unresolved auxiliary family，不可加入 primary HP 分母，也不能當第三個 clone。"),
        chart_block(builder, "candidate-combinations-single", "單一 primary HP：C=1…6、>6 與 incomplete",
                    "C 等於該 primary unit 的 n_trees；圖中比例分母為所有 one-primary regions，精確表另列 W_tree 總比例與 complete-only 比例。",
                    c_tables["single"], c_tables["single"], ["region", "derived"],
                    "C 是 exact candidate-tree combination count，不是群落數、clone 數或可信度分數。"),
        chart_block(builder, "candidate-combinations-double", "HP1+HP2：joint C=1…6、>6 與 incomplete",
                    "C_region=n_trees(HP1)×n_trees(HP2)；圖中分母為所有 two-primary regions。",
                    c_tables["double"], c_tables["double"], ["region", "derived"],
                    "乘積代表兩個獨立 primary reconstruction units 的 Cartesian combinations；任一 unit capped 時 C 不可用。"),
        chart_block(builder, "region-topology-classes", "Region 最終可能出現的三種可行拓撲判定",
                    "分母 W_primary；三個可行狀態加上 incomplete，精確表同時列 W_tree、W_primary 與 complete-only 比例。",
                    topology_table, topology_table, ["region", "derived"],
                    "C=1 必然 Topo=1；『exact tree 唯一但 topology 不唯一』在數學上不可能。"),
        chart_block(builder, "hidden-node-classes", "Hidden mutation-state nodes：0、>0 與 unavailable",
                    "n_hidden,region 為完整 primary HP units 的 minimal-hidden node 數加總；分母 W_primary。",
                    hidden_table, hidden_table + '<h4 class="subtable-title">依 HP scope × topology class 的 hidden 交叉</h4>' + topology_hidden_table, ["layered", "derived"],
                    "hidden=0 不是單次突變；hidden>0 也不是 hidden clone，只是未被 read pattern 直接觀察的中間 mutation state。"),
        chart_block(builder, "tree-outcomes", "Primary reconstruction unit 的 L1 結果",
                    "分母為 primary mutation-bearing HP1/HP2 units。",
                    outcome_table, outcome_table, ["layered"],
                    "模型內唯一只代表目前限制下候選樹唯一，不等於真實演化史被確認。"),
        chart_block(builder, "determinacy-levels", "Exact tree 與 shape 的模型內唯一性",
                    "圖中兩者共同以 analysis-complete non-capped primary units 為分母；region-level 指標只列在精確表。",
                    determinacy_table, determinacy_table, ["layered", "region"],
                    "Shape unique 只代表候選樹共享無標記形狀；不代表突變順序或細胞演化史唯一。"),
        chart_block(builder, "solver-cap-rate", "Solver-capped 比例",
                    "各 dataset 以 primary units 為分母；橫向比較搜索未完成的負擔。",
                    cap_table, cap_table, ["layered"],
                    "H2009 是明顯 outlier；原因仍需 k、coverage 與 partial evidence diagnostics，不能先下生物學結論。"),
        chart_block(builder, "sex-chromosome-share", "非 autosomal PASS sSNV 比例",
                    "分母為每個 caller PASS sSNV universe；目前主要由 chrX 構成。",
                    chromosome_table, chromosome_table, ["vcf"],
                    "樣本間差異極大，而 manifest 尚無 sex/karyotype/PAR contract，因此必須與 chr1–22 分層。"),
    ])

    hcc_f = hcc["site_funnel"]
    hcc_r = hcc["regions"]
    hcc_u = hcc["units"]
    hcc_read = hcc["read_evidence"]
    missing_region = hcc_r["missing_from_region_view"][0] if hcc_r["missing_from_region_view"] else "無"

    status_table = builder.exact_table(
        ["Evidence track", "Expected", "PASS", "Current artifact state", "Release role"],
        [
            ["Hash-verified historical snapshot", str(verification["dataset_count"]), str(verification["n_pass"]),
             "7/7 outputs + hashes", "可作內部數據觀察；不可作 clean publication evidence"],
            ["Normalized raw-all LongPhase-S production", str(production["expected_datasets"]), str(production["n_pass"]),
             f"state={production['state']}; PASS={production_passed_label}; active={'、'.join(production_active) or 'none'}; FAIL={production_failed_label}; aggregate={production['aggregate_exists']}/{production['aggregate_all_pass']}",
             {"complete": "release candidate", "failed": "FAIL／需修正", "in_progress": "進行中／PROBE", "not_started": "尚未開始"}[production["state"]]],
            ["candidate-tree read-AF ordering", str(auxiliary["read_af"]["expected_datasets"]),
             str(auxiliary["read_af"]["dataset_artifacts"]), "0/7 dataset artifacts；not_generated",
             "僅能作候選相對權重；不能當 CCF／矯正／確認"],
            ["L3 methylation", str(auxiliary["methylation"]["expected_datasets"]),
             str(auxiliary["methylation"]["evaluated_datasets"]), "7/7 not_evaluated；0 units evaluated",
             "介面保留 bounded auxiliary；本 run 沒有甲基結果"],
        ],
        [
            ["verification", "verification", "verification", "verification", "verification"],
            ["production", "production", "production", "production", "production"],
            ["read_af", "read_af", "read_af", "read_af", "read_af"],
            ["layered", "layered", "layered", "layered", "layered"],
        ],
        numeric_columns={1, 2},
    )

    unit_table = builder.exact_table(
        ["Layer", "Canonical unit", "Primary denominator", "Correct reading", "Do not read as"],
        [
            ["U1 site", "biallelic caller-PASS sSNV", "PASS U or chr1–22 U", "operational site universe", "biological truth"],
            ["W_pre", "multi-locus connected component", "all pre-cap groups", "region grouping before read support", "clone"],
            ["W_ret", "retained MLHP group", "W_pre", "group that passed cap/read gates", "emitted tree unit in every case"],
            ["W_tree", "emitted region-view record", "W_ret", "region with a detail record", "same denominator as W_ret when a gap exists"],
            ["W_primary", "region with ≥1 primary unit", "W_tree", "mutation-bearing HP1/2 reconstruction present", "confirmed subclone region"],
            ["k", "selected sSNVs per retained group", "W_ret", "solver/problem input size", "tree nodes or clone count"],
            ["Read", "read×region R/A/X observation", "full+partial observations or W_ret", "same-read pattern evidence", "unique molecule or cell count"],
            ["Primary unit", "region × mutation-bearing HP1/2", "all primary reconstruction units", "one reconstruction problem", "cell lineage or clone"],
            ["Candidate tree", "mutation-state graph candidate", "analysis-complete non-capped unit", "model-consistent solution", "posterior or true history"],
            ["L2 CN", "post-hoc annotation of recurrence-required unit", "recurrence-required units", "confound interpretation", "L1 solver input or causal proof"],
            ["read-AF", "candidate-set relative ordering artifact", "digest-matched complete candidate set", "future relative weighting", "caller VAF, family fraction, CCF or correction"],
            ["L3 methylation", "bounded auxiliary interface", "evaluated units only", "future residual/negative flag", "current evidence, ranking or confirmation"],
        ],
        "derived",
        numeric_columns=set(),
    )

    hcc_site_table = builder.exact_table(
        ["Branch", "Count", "Denominator", "Share"],
        [
            ["PASS sSNV universe U", f"{hcc_f['universe']:,}", f"{hcc_f['universe']:,}", "100.00%"],
            ["chr1–22", f"{hcc_f['autosomal']:,}", f"{hcc_f['universe']:,}", format_value(pct(hcc_f['autosomal'], hcc_f['universe']), 'percent')],
            ["chrX/Y out-of-scope", f"{hcc_f['out_of_scope']:,}", f"{hcc_f['universe']:,}", format_value(pct(hcc_f['out_of_scope'], hcc_f['universe']), 'percent')],
            ["positional singleton", f"{hcc_f['singleton']:,}", f"{hcc_f['autosomal']:,}", format_value(pct(hcc_f['singleton'], hcc_f['autosomal']), 'percent')],
            ["densest-8 cap-excluded", f"{hcc_f['cap_excluded']:,}", f"{hcc_f['autosomal']:,}", format_value(pct(hcc_f['cap_excluded'], hcc_f['autosomal']), 'percent')],
            ["read-unsupported", f"{hcc_f['read_unsupported']:,}", f"{hcc_f['autosomal']:,}", format_value(pct(hcc_f['read_unsupported'], hcc_f['autosomal']), 'percent')],
            ["retained", f"{hcc_f['retained']:,}", f"{hcc_f['autosomal']:,}", format_value(pct(hcc_f['retained'], hcc_f['autosomal']), 'percent')],
        ], "verification")

    hcc_region_table = builder.exact_table(
        ["Region stage", "Count", "Conservation"],
        [
            ["W_pre", f"{hcc_r['W_pre']:,}", f"{hcc_r['W_pre']:,} = {hcc_r['W_pre']-hcc_r['W_ret']:,} unsupported + {hcc_r['W_ret']:,} retained"],
            ["W_ret", f"{hcc_r['W_ret']:,}", f"{hcc_r['W_ret']:,} = {hcc_r['W_tree']:,} emitted + {len(hcc_r['missing_from_region_view']):,} missing"],
            ["W_tree", f"{hcc_r['W_tree']:,}", f"{hcc_r['W_tree']:,} = {hcc_r['W_primary']:,} primary + {hcc_r['no_primary']:,} no-primary"],
            ["W_primary", f"{hcc_r['W_primary']:,}", "mutation-bearing HP1/HP2 region"],
        ], "derived")

    hcc_layer_rows = [
        ["Site retained", f"{hcc_f['retained']:,}", f"{hcc_f['autosomal']:,}",
         format_value(pct(hcc_f['retained'], hcc_f['autosomal']), "percent"), "caller-PASS sSNV / chr1–22"],
        ["Full R/A patterns", f"{hcc_read['full_pattern_observations']:,}", f"{hcc_read['retained_pattern_observations']:,}",
         format_value(pct(hcc_read['full_pattern_observations'], hcc_read['retained_pattern_observations']), "percent"), "read×region observations"],
        ["Partial-X patterns", f"{hcc_read['partial_pattern_observations']:,}", f"{hcc_read['retained_pattern_observations']:,}",
         format_value(pct(hcc_read['partial_pattern_observations'], hcc_read['retained_pattern_observations']), "percent"), "read×region observations"],
        ["Regions with full read", f"{hcc_read['regions_with_full']:,}", f"{hcc_r['W_ret']:,}",
         format_value(pct(hcc_read['regions_with_full'], hcc_r['W_ret']), "percent"), "region denominator"],
        ["Regions with ALT full population", f"{hcc_read['regions_with_alt_population']:,}", f"{hcc_r['W_ret']:,}",
         format_value(pct(hcc_read['regions_with_alt_population'], hcc_r['W_ret']), "percent"), "region denominator"],
        ["HP region partition", f"{hcc_r['single_primary_HP']:,} one; {hcc_r['multi_HP']:,} two; {hcc_r['no_primary']:,} none",
         f"{hcc_r['W_tree']:,}", "100.00%", "reconstruction-unit multiplicity; not clones"],
        ["L1 outcome partition", f"{hcc_u['determined']:,} unique; {hcc_u['ambiguous']:,} ambiguous; {hcc_u['solver_capped']:,} capped; {hcc_u['recurrence']:,} recurrence-required",
         f"{hcc_u['primary']:,}", "100.00%", "primary reconstruction units"],
        ["Exact-tree unique", f"{hcc_u['exact_tree_unique']:,}", f"{hcc_u['noncapped_complete']:,}",
         format_value(pct(hcc_u['exact_tree_unique'], hcc_u['noncapped_complete']), "percent"), "analysis-complete non-capped units"],
        ["Shape unique", f"{hcc_u['shape_unique']:,}", f"{hcc_u['noncapped_complete']:,}",
         format_value(pct(hcc_u['shape_unique'], hcc_u['noncapped_complete']), "percent"), "shape only; order may differ"],
        ["Region all-exact", f"{hcc_u['region_all_exact']:,}", f"{hcc_r['W_primary']:,}",
         format_value(pct(hcc_u['region_all_exact'], hcc_r['W_primary']), "percent"), "all primary units exact within region"],
        ["Candidate-tree burden", f"{hcc_u['total_candidate_trees']:,}", f"{hcc_u['noncapped_complete']:,}",
         format_value(pct(hcc_u['total_candidate_trees'], hcc_u['noncapped_complete']), "decimal2"), "mean candidates/unit; not clone count"],
        ["L2 model recurrence-required", f"{hcc['L2_cn']['n_primary_recurrence_sent_to_cn']:,}", f"{hcc_u['primary']:,}",
         f"20 CN-amp-confounded; 7 LOH unresolved", "post-hoc annotation; not solver input"],
        ["chrX caller branch", f"{hcc['sex_chromosome_census']['chrX']:,}", f"{hcc_f['universe']:,}",
         format_value(pct(hcc['sex_chromosome_census']['chrX'], hcc_f['universe']), "percent"), "historical caller universe; not mutation burden"],
    ]
    hcc_layer_sources = []
    for row_index in range(len(hcc_layer_rows)):
        key = "mlhp" if 1 <= row_index <= 4 else "region" if row_index == 5 else "layered"
        if row_index == 0:
            key = "verification"
        elif row_index == 12:
            key = "vcf"
        hcc_layer_sources.append([key] * 5)
    hcc_layer_table = builder.exact_table(
        ["Layer", "Numerator/detail", "Denominator", "Rate/mean", "Correct interpretation"],
        hcc_layer_rows, hcc_layer_sources, numeric_columns={1, 2, 3})

    pooled_topology = Counter()
    pooled_hidden_classes = Counter()
    pooled_c_bins = Counter()
    for sample in samples:
        pooled_topology.update(sample["region_candidate_topology"]["topology_classes"])
        pooled_hidden_classes.update(sample["region_candidate_topology"]["hidden_classes"])
        pooled_c_bins.update(sample["region_candidate_topology"]["C_bins"]["all"])

    logic_table = builder.exact_table(
        ["Formal state", "Exact tree", "Topology", "Hidden nodes", "Most-likely wording", "Feasibility"],
        [
            ["C=1, Topo=1", "唯一", "唯一", "可為 0 或 >0", "sole model-compatible candidate；不是 calibrated most-likely", "可行"],
            ["C>1, Topo=1", "不唯一", "唯一", "可為 0 或 >0", "需要額外 ranking 才能挑 exact tree；本輪未評估", "可行"],
            ["C>1, Topo>1", "不唯一", "不唯一", "可為 0 或 >0", "需要額外 ranking 才能挑 topology；本輪未評估", "可行"],
            ["C=1, Topo>1", "唯一", "不唯一", "—", "—", "數學上不可能"],
            ["incomplete / capped", "不可判", "不可判", "不可比較", "ranking 不適用", "候選全集未完成"],
        ],
        "derived", numeric_columns=set())

    definition_details = "".join([
        '<details id="definition-s" class="concept-detail"><summary>S／U：sSNV 位點宇宙</summary><p>本報告的 S/U 是 historical caller-PASS biallelic sSNV。它是操作性輸入，不是 biological truth；raw ClairS-all 是另一個上游契約。</p></details>',
        '<details id="definition-w" class="concept-detail"><summary>W：region／group 的各階段</summary><p>W_all,pre 先包含 k=1 singleton components 與 k&gt;1 groups；W_pre/W_ret/W_tree/W_primary 只沿 k&gt;1 多位點重建路徑。S→W 是單位轉換，不能當 site retention rate。</p></details>',
        '<details id="definition-k" class="concept-detail"><summary>k：一個 group 選入的 sSNV 數</summary><p>k 是重建問題大小；k=1 不進多位點樹 solver，k=2–8 進 retained MLHP。k 不是樹節點數、clone 數或 subclone 數。</p></details>',
        '<details id="definition-hp" class="concept-detail"><summary>HP1／HP2：primary reconstruction units</summary><p>一個 primary unit = region × mutation-bearing HP1/2。HP 是等位 family context；HP1、HP2 都有 unit 只代表兩個 allele-specific reconstruction problems，不代表兩個 cell clones。</p></details>',
        '<details id="definition-h3" class="concept-detail"><summary>H3?：unresolved auxiliary HP context</summary><p>H3? 可與 0、1、2 個 primary units 共存；它不進 primary、C 或 determinacy 分母，也不是第三個 haplotype、clone 或 subclone。</p></details>',
        '<details id="definition-c" class="concept-detail"><summary>C：region 的 joint exact-tree combinations</summary><p>單一 primary 時 C=n_trees；HP1+HP2 時 C=n_trees(HP1)×n_trees(HP2)。任一 unit incomplete/capped 時 C 不可用，不能記為 0。</p></details>',
        '<details id="definition-topo" class="concept-detail"><summary>Topo：joint unlabeled-shape combinations</summary><p>Topo_region 是各 primary unit 的 n_distinct_shapes_exact 乘積。Topo=1 只代表候選共享同一無標記形狀，不保證 mutation labels、順序或真實演化史唯一。</p></details>',
        '<details id="definition-hidden" class="concept-detail"><summary>Hidden：parsimony-inferred mutation-state nodes</summary><p>region 值是完整 primary units 的 minimal-hidden node 數加總。hidden=0 不等於單次突變；hidden&gt;0 也不是 hidden clones。</p></details>',
        '<details id="definition-rate" class="concept-detail"><summary>總比例與父層相對比例</summary><p>總比例回答「占全體多少」；父層相對比例回答「通過上一層後占多少」。每列都必須同時標 numerator、denominator 與 unit。</p></details>',
        '<details id="definition-ranking" class="concept-detail"><summary>Most-likely／read-AF：目前未評估</summary><p>candidate-tree read-AF artifacts=0/7，ranked regions=0；沒有 calibrated likelihood、posterior 或 winner tree。caller VAF、family read fraction 與 candidate-tree ordering 是三個不同概念。</p></details>',
        '<details id="definition-methyl" class="concept-detail"><summary>L3 methylation：介面保留、結果未評估</summary><p>7/7 datasets 為 not_evaluated，0 units evaluated；本輪甲基沒有參與候選樹生成、排序、修正或確認。</p></details>',
    ])

    pipeline_svg = '''
    <figure class="svg-explainer" aria-labelledby="pipeline-svg-title">
      <figcaption id="pipeline-svg-title"><strong>S → W → HP → C → Topo：</strong>每個箭頭的單位、分母與證據角色</figcaption>
      <svg viewBox="0 0 1160 390" role="img" aria-labelledby="pipeline-svg-title pipeline-svg-desc">
        <desc id="pipeline-svg-desc">Observed sSNV and read evidence becomes region groups, allele-specific primary reconstruction units, exact candidate combinations and topology classes. Read-AF and methylation remain unavailable auxiliary evidence.</desc>
        <defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" class="svg-arrow-head"/></marker></defs>
        <path d="M184 126H254" class="svg-arrow"/><path d="M424 126H494" class="svg-arrow"/><path d="M664 126H734" class="svg-arrow"/><path d="M904 126H974" class="svg-arrow"/>
        <g class="svg-node observed" transform="translate(24 70)"><rect width="160" height="112" rx="18"/><text x="18" y="30" class="svg-step">01 · OBSERVED</text><text x="18" y="58" class="svg-title">S / U</text><text x="18" y="82">caller-PASS sSNV</text><text x="18" y="101">單位：site</text></g>
        <g class="svg-node derived" transform="translate(254 70)"><rect width="170" height="112" rx="18"/><text x="18" y="30" class="svg-step">02 · GROUPING</text><text x="18" y="58" class="svg-title">W + k</text><text x="18" y="82">k1 vs k&gt;1 groups</text><text x="18" y="101">單位：region/group</text></g>
        <g class="svg-node derived" transform="translate(494 70)"><rect width="170" height="112" rx="18"/><text x="18" y="30" class="svg-step">03 · ALLELIC</text><text x="18" y="58" class="svg-title">HP1 / HP2</text><text x="18" y="82">0 / 1 / 2 primary units</text><text x="18" y="101">H3? auxiliary</text></g>
        <g class="svg-node model" transform="translate(734 70)"><rect width="170" height="112" rx="18"/><text x="18" y="30" class="svg-step">04 · ENUMERATE</text><text x="18" y="58" class="svg-title">C</text><text x="18" y="82">joint exact candidates</text><text x="18" y="101">1…N / incomplete</text></g>
        <g class="svg-node model" transform="translate(974 70)"><rect width="160" height="112" rx="18"/><text x="18" y="30" class="svg-step">05 · CLASSIFY</text><text x="18" y="58" class="svg-title">Topo</text><text x="18" y="82">unique / multiple</text><text x="18" y="101">+ hidden state</text></g>
        <line x1="334" y1="205" x2="334" y2="240" class="unit-divider"/><text x="348" y="228" class="svg-note">S→W：many sites → one component；比例在此重新起算</text>
        <g class="aux-rail"><rect x="494" y="272" width="640" height="86" rx="18"/><text x="516" y="298" class="svg-step">AUXILIARY RAIL · 不進入主鏈的候選生成</text><text x="516" y="325" class="svg-title-sm">CN：post-hoc confound only</text><text x="770" y="325" class="svg-title-sm">read-AF：0/7 · NOT EVALUATED</text><text x="770" y="347" class="svg-title-sm">Methylation：0 units · NOT EVALUATED</text></g>
        <g class="svg-legend" transform="translate(24 300)"><circle cx="8" cy="0" r="7" class="legend-observed"/><text x="22" y="5">Observed evidence</text><circle cx="8" cy="28" r="7" class="legend-derived"/><text x="22" y="33">Derived grouping</text><circle cx="8" cy="56" r="7" class="legend-model"/><text x="22" y="61">Model-derived candidate</text></g>
      </svg>
    </figure>'''

    topology_svg = '''
    <figure class="svg-explainer topology-atlas" aria-labelledby="topology-svg-title">
      <figcaption id="topology-svg-title"><strong>Topology atlas：</strong>三個可行狀態與一個不可能狀態</figcaption>
      <svg viewBox="0 0 1160 330" role="img" aria-labelledby="topology-svg-title topology-svg-desc">
        <desc id="topology-svg-desc">Four panels show exact and topology uniqueness: unique exact tree, multiple exact trees with one shape, multiple shapes, and the impossible state of one exact tree with multiple topologies.</desc>
        <g class="topo-panel ok" transform="translate(18 28)"><rect width="265" height="270" rx="20"/><text x="20" y="30" class="svg-step">A · C=1, Topo=1</text><path d="M66 83L132 135L198 83" class="tree-edge"/><circle cx="66" cy="83" r="12"/><circle cx="132" cy="135" r="12"/><circle cx="198" cy="83" r="12"/><text x="20" y="190" class="svg-title-sm">Exact + topology unique</text><text x="20" y="217">sole model-compatible candidate</text><text x="20" y="242">hidden 可為 0 或 &gt;0</text></g>
        <g class="topo-panel ok" transform="translate(303 28)"><rect width="265" height="270" rx="20"/><text x="20" y="30" class="svg-step">B · C&gt;1, Topo=1</text><path d="M50 83L116 135L182 83" class="tree-edge"/><path d="M88 75L154 127L220 75" class="tree-edge ghost"/><circle cx="50" cy="83" r="10"/><circle cx="116" cy="135" r="10"/><circle cx="182" cy="83" r="10"/><text x="20" y="190" class="svg-title-sm">Shape-only unique</text><text x="20" y="217">labels/order differ；shape same</text><text x="20" y="242">winner 未評估</text></g>
        <g class="topo-panel warn" transform="translate(588 28)"><rect width="265" height="270" rx="20"/><text x="20" y="30" class="svg-step">C · C&gt;1, Topo&gt;1</text><path d="M54 83L116 135L178 83" class="tree-edge"/><path d="M204 76L204 140M204 105L236 140" class="tree-edge ghost"/><circle cx="54" cy="83" r="10"/><circle cx="116" cy="135" r="10"/><circle cx="178" cy="83" r="10"/><text x="20" y="190" class="svg-title-sm">Topology non-unique</text><text x="20" y="217">多個 exact + 多個 shapes</text><text x="20" y="242">winner 未評估</text></g>
        <g class="topo-panel impossible" transform="translate(873 28)"><rect width="265" height="270" rx="20"/><text x="20" y="30" class="svg-step">D · C=1, Topo&gt;1</text><path d="M58 76L205 150M205 76L58 150" class="cross"/><text x="20" y="190" class="svg-title-sm">不可能狀態</text><text x="20" y="217">一棵 exact tree 只能對應</text><text x="20" y="239">一個 unlabeled shape</text></g>
      </svg>
    </figure>'''

    dataset_inspections = []
    for sample in samples:
        metrics = sample["region_candidate_topology"]
        region = sample["regions"]
        funnel_detail = builder.exact_table(
            ["Dataset", "Layer", "Stage", "Unit", "Count", "Total denominator · share", "Parent denominator · relative share"],
            dual_funnel_rows(sample["sample"], sample["site_funnel"], region),
            "derived", numeric_columns={4, 5, 6})
        hp_counts = metrics["hp_h3"]
        classification_rows = [
            ["HP×H3", "single H3− / H3+", f"{hp_counts['single_without_h3']:,} / {hp_counts['single_with_h3']:,}", f"single parent={region['single_primary_HP']:,}", "H3? 不進 primary 分母"],
            ["HP×H3", "double H3− / H3+", f"{hp_counts['double_without_h3']:,} / {hp_counts['double_with_h3']:,}", f"double parent={region['multi_HP']:,}", "HP1+HP2 不是兩個 clones"],
            ["C", "one primary", "; ".join(f"{key}:{value:,}" for key, value in metrics["C_bins"]["single"].items()), f"parent={region['single_primary_HP']:,}", "exact candidate combinations"],
            ["C", "HP1 and HP2", "; ".join(f"{key}:{value:,}" for key, value in metrics["C_bins"]["double"].items()), f"parent={region['multi_HP']:,}", "joint Cartesian combinations"],
            ["Topo", "C1/T1 · C&gt;1/T1 · C&gt;1/T&gt;1 · incomplete", f"{metrics['topology_classes']['exact_and_topology_unique']:,} · {metrics['topology_classes']['topology_unique_exact_multiple']:,} · {metrics['topology_classes']['topology_multiple_exact_multiple']:,} · {metrics['topology_classes']['incomplete']:,}", f"W_primary={region['W_primary']:,}", "model-internal uniqueness"],
            ["Hidden", "H=0 · H&gt;0 · incomplete", f"{metrics['hidden_classes']['hidden_zero']:,} · {metrics['hidden_classes']['hidden_positive']:,} · {metrics['hidden_classes']['incomplete']:,}", f"W_primary={region['W_primary']:,}", "mutation-state nodes；不是 clones"],
            ["Most-likely", "candidate-tree ranking", "NOT EVALUATED", "0 ranked regions", "沒有 winner/posterior"],
            ["read-AF", "candidate-tree ordering artifact", "0/7 datasets", "0 evaluated units", "caller VAF 與此不同"],
            ["L3", "methylation", "not_evaluated", "0 evaluated units", "未參與排序或修正"],
            ["chrX", "caller census only", f"{sample['sex_chromosome_census']['chrX']:,}", f"PASS U={sample['site_funnel']['universe']:,}", "不進 chr1–22 topology"],
        ]
        classification_table = builder.exact_table(
            ["Layer", "Class", "Count/detail", "Denominator", "Correct interpretation"],
            classification_rows, "derived", numeric_columns={2, 3})
        duplicate_note = "；與 HCC1395 為同一 biological sample 的技術處理版本" if sample["sample"] == "HCC1395_DORADO" else ""
        dataset_inspections.append(f'''
        <details class="dataset-detail" id="dataset-{html.escape(sample['sample'].lower().replace('_', '-'))}">
          <summary>{html.escape(sample['sample'])}：展開完整 S→W→HP→C→Topo 檢視</summary>
          <div class="dataset-meta"><span>Biological ID: {html.escape(str(sample['biological_id']))}{duplicate_note}</span><span>Historical snapshot · 6/7 truth-BED-conditioned context</span></div>
          <h3>完整雙單位 funnel</h3>{funnel_detail}
          <h3>HP、C、Topo、hidden 與未評估證據</h3>{classification_table}
          <p class="legal-reading"><strong>合法結論：</strong>可描述此 dataset 的 regional mutation-state candidate 組成與模型內唯一性。<br><strong>禁止結論：</strong>不可把 HP、C、Topo 或 hidden counts 直接稱為 clone/subclone 數、比例或真實祖先關係。</p>
        </details>''')
    dataset_inspections = "".join(dataset_inspections)

    production_state = {
        "complete": "完成（aggregate all_pass）",
        "failed": "失敗／已中止",
        "in_progress": "進行中",
        "not_started": "尚未開始",
    }[production["state"]]
    data_drawer = f'''
    <details class="data-drawer">
      <summary>資料與驗證</summary>
      <div class="data-link-grid">
        <a href="{html.escape(DATA_SNAPSHOT_NAME)}" aria-label="開啟可重算數據快照 JSON"><strong>可重算數據</strong><span>本頁所有展示值</span></a>
        <a href="{html.escape(SOURCE_NOTES_NAME)}" aria-label="開啟圖表與來源對照 JSON"><strong>圖表與來源</strong><span>claim-to-source mapping</span></a>
        <a href="{html.escape(VALIDATION_NAME)}" aria-label="開啟 artifact validation JSON"><strong>Artifact gate</strong><span>結構與數據驗證</span></a>
        <a href="{html.escape(BROWSER_QA_NAME)}" aria-label="開啟 browser QA JSON"><strong>Browser QA</strong><span>桌機／手機驗證</span></a>
      </div>
    </details>'''
    main_content = f'''
    <a class="skip-link" href="#report-main">跳到主要內容</a>
    <div class="shell">
      <header class="topbar">
        <div class="brand"><span class="mark" aria-hidden="true"></span>InterSubMod Research</div>
        <div class="meta">資料截止 · {html.escape(generated_display)}</div>
      </header>
      <main id="report-main" tabindex="-1" data-report-audience="technical">
        <article class="reading report-hero">
          <div class="validation-ribbon">PARTIAL SCIENTIFIC VALIDATION · PUBLICATION NO-GO</div>
          <div class="blossom" aria-hidden="true"><i></i><i></i><i></i><i></i><i></i></div>
          <div class="kicker">分層重建 · 目前進度 + 歷史基線 · 技術報告</div>
          <header data-contract-section="title"><h1>InterSubMod 分層重建：目前進度、歷史基線與證據邊界</h1></header>
          <p class="deck">先回答目前資料能不能用，再把 historical engineering snapshot 降為方法與分母的參考附錄；任何舊數字都不能替代尚未完成的 normalized raw-all canonical run。</p>
          <section class="summary" data-contract-section="technical-summary">
            <div class="summary-label">Technical Summary</div>
            <div class="summary-body">
              <p><strong>目前 canonical production：</strong>{html.escape(production_readout)}。{html.escape(production_boundary)}</p>
              <p><strong>歷史數據只能作工程基線：</strong>{builder.number(verification['n_pass'], source_key='verification')}/{builder.number(verification['dataset_count'], source_key='verification')} dataset outputs 的 artifact hashes 通過，但 {builder.number(6, source_key='scope')}/{builder.number(7, source_key='scope')} 舊 HP/PS inputs 受 truth-BED conditioning，且目前 code snapshot 只有 {builder.number(verification['code_hash_manifest_pass'], source_key='code_hash')}/{builder.number(verification['code_hash_manifest_total'], source_key='code_hash')} 符合記錄。</p>
              <p><strong>未閉合證據：</strong>candidate-tree read-AF 為 {builder.number(auxiliary['read_af']['dataset_artifacts'], source_key='read_af')}/{builder.number(auxiliary['read_af']['expected_datasets'], source_key='read_af')} artifacts；L3 為 {builder.number(auxiliary['methylation']['evaluated_datasets'], source_key='layered')}/{builder.number(auxiliary['methylation']['expected_datasets'], source_key='layered')} datasets evaluated；HCC1395 尚有 {builder.number(len(hcc_r['missing_from_region_view']), source_key='derived')} region／{builder.number(hcc_r['site_sum_delta'], source_key='derived')} sSNV 的 region-view 缺口。</p>
              <p><strong>主張上限：</strong>本 snapshot 產生 regional mutation-state candidate sets；部分結構與樣本內異質性相容，但不能建立細胞 clone/subclone 的存在、比例、祖先關係或完整演化史。</p>
            </div>
          </section>
          <section class="metrics">
            <div class="metric current-metric"><div class="metric-label">Raw-all production PASS</div><div class="metric-value">{builder.number(production['n_pass'], source_key='production')}/{builder.number(production['expected_datasets'], source_key='production')}</div><div class="metric-note">{html.escape('進行中：' + '、'.join(production_active) if production_active else production_state)}</div></div>
            <div class="metric"><div class="metric-label">Historical artifact outputs</div><div class="metric-value">{builder.number(verification['n_pass'], source_key='verification')}/{builder.number(verification['dataset_count'], source_key='verification')}</div><div class="metric-note">僅 historical engineering baseline</div></div>
            <div class="metric"><div class="metric-label">Biological samples</div><div class="metric-value">{builder.number(verification['biological_sample_count'], source_key='verification')}</div><div class="metric-note">7 dataset rows；HCC1395/DORADO 同一生物樣本</div></div>
            <div class="metric"><div class="metric-label">Scientific release</div><div class="metric-value status-no-go">NO-GO</div><div class="metric-note">可內部分享；不可寫成 biological confirmation</div></div>
          </section>
          {data_drawer}
        </article>

        <section class="reading current-gate" id="current-status" aria-labelledby="current-status-title">
          <div>
            <span class="gate-label">CURRENT · NORMALIZED RAW-ALL</span>
            <h2 id="current-status-title">{html.escape(production_readout)}</h2>
            <p>{html.escape(production_boundary)} 截止 {html.escape(generated_display)}；本頁是靜態 snapshot。</p>
          </div>
          <div class="gate-progress">
            <progress value="{production['n_pass']}" max="{production['expected_datasets']}" aria-label="normalized raw-all production PASS 進度"></progress>
            <span>PASS：{html.escape(production_passed_label)}</span>
            <span>Active：{html.escape('、'.join(production_active) or '無')}</span>
          </div>
        </section>

        <nav class="reading mini-toc" aria-label="報告章節">
          <strong>快速閱讀</strong>
          <a href="#current-status">目前狀態</a><a href="#professor-flow">方法主鏈</a><a href="#concept-glossary">名詞定義</a>
          <a href="#topology-atlas">拓撲判斷</a><a href="#historical-atlas">歷史附錄</a><a href="#scope-limitations">限制</a>
        </nav>

        <article class="reading" id="professor-flow">
          <section class="narrative">
            <span class="evidence-badge observed">Observed evidence</span><span class="evidence-badge model">Model-derived candidate</span><span class="evidence-badge unavailable">Unavailable / not evaluated</span>
            <h2>先看懂主鏈：資料不是一路變成 clone，而是逐層改變分析單位</h2>
            <p>單位依序是 site → region/group → primary reconstruction unit → joint candidate combination → topology class。CN、candidate-tree read-AF 與甲基化位在 auxiliary rail；它們不生成候選樹，其中 read-AF 與甲基在本輪尚未評估。</p>
          </section>
          {pipeline_svg}
          <details class="technical-appendix compact-appendix" id="complete-funnel">
            <summary><span>Historical aggregate 指標與雙分母漏斗</span><small>7 dataset rows／6 biological samples；只作描述性 inventory</small></summary>
            <div class="appendix-body">
              <section class="metrics topology-metrics">
                <div class="metric"><div class="metric-label">W_primary · 7 dataset rows</div><div class="metric-value">{builder.number(pooled_wprimary, source_key='derived')}</div><div class="metric-note">dataset-row aggregate；不是 7 個獨立生物樣本</div></div>
                <div class="metric"><div class="metric-label">Complete candidate regions</div><div class="metric-value">{builder.number(pooled_complete_regions, source_key='derived')}</div><div class="metric-note">{builder.number(pct(pooled_complete_regions, pooled_wprimary), 'percent', 'derived')} of W_primary</div></div>
                <div class="metric"><div class="metric-label">C=1 · Topo=1</div><div class="metric-value">{builder.number(pooled_topology['exact_and_topology_unique'], source_key='derived')}</div><div class="metric-note">sole model-compatible exact combination</div></div>
                <div class="metric"><div class="metric-label">C&gt;1 · Topo&gt;1</div><div class="metric-value">{builder.number(pooled_topology['topology_multiple_exact_multiple'], source_key='derived')}</div><div class="metric-note">多 exact candidates 且多 shapes</div></div>
              </section>
              <section class="narrative">
                <h2>Site 與 region 必須使用兩套分母</h2>
                <p>總比例從該單位全體出發；父層相對比例只回答上一層之後留下多少。S→W 處會重設分母，不能把 autosomal sites 與 components 串成同一守恆流。</p>
              </section>
              <section class="card table-card"><div class="card-head"><h3>Historical 7-row aggregate：S/U → W → primary</h3><p>HCC1395/DORADO 重複同一 biological sample，故只作描述性加總。</p></div>{pooled_funnel_table}</section>
            </div>
          </details>

          <section id="concept-glossary" class="narrative">
            <h2>名詞與概念：需要細節時再展開</h2>
            <p>每個定義都同時說明「是什麼」與「不能解讀成什麼」，避免 HP、candidate tree、hidden node 被誤寫成 clone/subclone。</p>
          </section>
          <div class="concept-grid">{definition_details}</div>

          <section id="topology-atlas" class="narrative">
            <h2>最後會出現三種可行拓撲狀態；第四種在數學上不可能</h2>
            <p>在完整候選集內必有 C≥Topo≥1。`C=1` 時 exact candidate 只有一棵，因此它只能對應一種 shape；hidden=0 或 hidden&gt;0 都不會改變這個邏輯。</p>
          </section>
          {topology_svg}
          <section class="card table-card"><div class="card-head"><h3>C、Topo、hidden 與 most-likely 的判斷矩陣</h3><p>「sole candidate」是集合大小；「most likely」是排序結果。兩者不能互換。</p></div>{logic_table}</section>
          <section class="na-panel" aria-label="Most-likely topology 未評估">
            <span class="evidence-badge unavailable">Unavailable / not evaluated</span>
            <div><strong>Most-likely topology：NOT EVALUATED</strong><p>candidate-tree read-AF artifacts = {builder.number(auxiliary['read_af']['dataset_artifacts'], source_key='read_af')}/{builder.number(auxiliary['read_af']['expected_datasets'], source_key='read_af')} datasets；ranked regions = {builder.number(0, source_key='read_af')}。沒有 calibrated ranking、posterior 或 winner tree。</p></div>
          </section>
        </article>

        <article class="reading">
          <section class="narrative" data-contract-section="key-findings">
            <h2>Current 與 historical 是兩條不同 evidence track</h2>
            <p>{html.escape(production_readout)}。Historical 圖表仍保留作方法與分母稽核，但 6/7 HP/PS inputs 受 truth-BED conditioning，且目前只有 {builder.number(verification['code_hash_manifest_pass'], source_key='code_hash')}/{builder.number(verification['code_hash_manifest_total'], source_key='code_hash')} <code>code.sha256</code> entries 符合記錄；兩者不得合併成同一組 current rates。</p>
          </section>
          <section class="card table-card"><div class="card-head"><h3>Evidence-track status</h3><p>截止時間、完成狀態與合法用途。</p></div>{status_table}</section>
          <section class="narrative latest-flow">
            <h2>Canonical 資料流：normalized raw-all 進 LongPhase-S，recalibrated PASS 才進 tree</h2>
            <div class="flow-grid">
              <div><span>1</span><strong>raw ClairS all</strong><p>先正規化，同時保留 lossless site ledger 與原始 disposition。</p></div>
              <div><span>2</span><strong>normalized raw-all → LongPhase-S</strong><p>PASS 與 non-PASS 都進 native recalibration，rescue/remove 分支才可被驗證。</p></div>
              <div><span>3</span><strong>LongPhase-S <code>_sc.vcf</code> all</strong><p>record keys 守恆、保存 recalibrated FILTER，並產出 clean HP/PS sidecar。</p></div>
              <div><span>4</span><strong><code>_sc.vcf</code> PASS</strong><p>最新 primary sSNV tree backbone；需 7-dataset record-key gates。</p></div>
              <div><span>5</span><strong>Original BAM + exact sidecar join</strong><p>輸出 R/A/X 與 HP family evidence；equivalence 未過前 fail closed。</p></div>
              <div><span>6</span><strong>Region → candidate trees</strong><p>W、k、read、HP、L1/L2/L3 逐層保留單位與分母。</p></div>
            </div>
            <div class="caveat"><strong>圖表資料邊界：</strong>下列 7-dataset 圖仍是 historical ClairS-PASS snapshot，只用來說明方法與量化分布。Active production root <code>{html.escape(production_root_name)}</code> 目前為 {builder.number(production['n_pass'], source_key='production')}/{builder.number(production['expected_datasets'], source_key='production')} PASS；aggregate `_SUCCESS` 與 fresh downstream canonical metrics 尚未完成，不能用舊數字代替。</div>
          </section>
        </article>

        <details class="wide technical-appendix chart-atlas" id="historical-atlas">
          <summary><span>Historical 圖表附錄</span><small>14 張圖與精確表；非 current canonical result</small></summary>
          <div class="appendix-body chart-stack">{charts}</div>
        </details>

        <article class="reading" id="dataset-views">
          <section class="narrative"><h2>7 個 dataset 的單獨檢視</h2><p>每個摺疊區塊都可獨立重讀完整 S→W→HP→C→Topo，並保留自己的分母、H3? 交叉、未評估證據與合法解釋上限。</p></section>
          <div class="dataset-list">{dataset_inspections}</div>
        </article>

        <article class="reading">
          <section class="narrative">
            <h2>候選樹負擔、hidden mutation-state nodes 與 L2 CN 必須另列</h2>
            <p>候選樹與 hidden nodes 是模型狀態，不是 clone 數；L2 CN 是對 model recurrence-required units 的事後分層，不進入 L1 solver。</p>
          </section>
          <details class="technical-appendix">
            <summary><span>展開 candidate-tree／hidden／L2 精確表</span><small>Historical audit detail</small></summary>
            <div class="appendix-body">
              <section class="card table-card"><div class="card-head"><h3>Candidate-tree and hidden-node burden</h3><p>只用 analysis-complete non-capped primary units；capped fallback 不混入可比較分母。</p></div>{complexity_table}</section>
              <section class="card table-card"><div class="card-head"><h3>L2 post-hoc CN classification</h3><p>各 dataset 的 CN source 不同；不可 pooled 成同一個 biological proportion。</p></div>{l2_table}</section>
            </div>
          </details>

          <section class="narrative">
            <h2>HCC1395 的雙漏斗揭露真正缺口：位點守恆，但 retained region 並未完全 emit</h2>
            <p>位點漏斗與 region 漏斗使用不同單位，必須分開閱讀。最大的位點損失來自 densest-8 工程 cap；region 層只有一個 retained group 沒有輸出 family/detail unit，但這一個缺口使 retained sSNV 與 region-view Σk 相差四。</p>
          </section>
          <details class="technical-appendix">
            <summary><span>展開 HCC1395 雙漏斗與全層精確表</span><small>1 region／4 sSNV gap 的可重算證據</small></summary>
            <div class="appendix-body">
              <section class="card table-card"><div class="card-head"><h3>HCC1395 site funnel</h3><p>分母在每列明列；不同分母不可直接比較百分比。</p></div>{hcc_site_table}</section>
              <section class="card table-card"><div class="card-head"><h3>HCC1395 region W funnel</h3><p>W_pre → W_ret → W_tree → W_primary。</p></div>{hcc_region_table}</section>
              <section class="card table-card"><div class="card-head"><h3>HCC1395 all-layer quantitative audit</h3><p>read、HP、L1、候選樹、L2 與 chrX 各自保留分母與單位。</p></div>{hcc_layer_table}</section>
            </div>
          </details>
          <div class="caveat"><strong>QA warning.</strong> 缺少的 region 是 <code>{html.escape(missing_region)}</code>；MLHP retained sites={builder.number(hcc_r['mlhp_site_sum'], source_key='mlhp')}，region-view Σk={builder.number(hcc_r['region_view_site_sum'], source_key='region')}。現有 final verifier 沒有攔截這個差異。</div>

          <section class="narrative">
            <h2>四層案例讓整體、重現性、極端值與可解釋區域同時可見</h2>
            <p>Aggregate 用來看全局；HCC1395/DORADO 是同一生物樣本的技術比較；H2009 是 solver-cap outlier；chr2 個案則示範模型內唯一但仍受 partial-read 限制的區域。</p>
          </section>
          <section class="case-grid">
            <div class="case-card"><span class="case-label">Aggregate</span><h3>7 datasets / 6 biological samples</h3><p>{builder.number(sum(s['site_funnel']['retained'] for s in samples), source_key='derived')} retained sSNV，{builder.number(sum(s['regions']['W_tree'] for s in samples), source_key='derived')} emitted regions，{builder.number(sum(s['units']['primary'] for s in samples), source_key='derived')} primary units。這是 dataset-level 加總，不是獨立生物樣本加權。</p></div>
            <div class="case-card"><span class="case-label">Technical replicate</span><h3>HCC1395 vs DORADO</h3><p>multi-HP marginal rate={builder.number(pct(hcc_r['multi_HP'], hcc_r['W_tree']), 'percent', 'region')} vs {builder.number(pct(by_name['HCC1395_DORADO']['regions']['multi_HP'], by_name['HCC1395_DORADO']['regions']['W_tree']), 'percent', 'region')}；model-unique rate={builder.number(pct(hcc_u['determined'], hcc_u['primary']), 'percent', 'layered')} vs {builder.number(pct(by_name['HCC1395_DORADO']['units']['determined'], by_name['HCC1395_DORADO']['units']['primary']), 'percent', 'layered')}。這只是 marginal-rate similarity，尚未做 locus-matched reproducibility，也不是跨生物樣本驗證。</p></div>
            <div class="case-card"><span class="case-label">Extreme outlier</span><h3>H2009 solver capped</h3><p>{builder.number(h2009['units']['solver_capped'], source_key='layered')}/{builder.number(h2009['units']['primary'], source_key='layered')}={builder.number(pct(h2009['units']['solver_capped'], h2009['units']['primary']), 'percent', 'layered')} primary units。以 unit flags 重算有 {builder.number(h2009['region_candidate_topology']['incomplete_regions'], source_key='derived')} incomplete regions；舊 <code>region_determinacy=has_capped</code> 只有 {builder.number(h2009['region_candidate_topology']['legacy_has_capped_label_count'], source_key='region')}，因 {builder.number(len(h2009['region_candidate_topology']['capped_label_discordance']), source_key='derived')} 個 recurrence+capped 混合區被優先標成 recurrence。先做 coverage/k/partial diagnostics，不能直接歸因為腫瘤較複雜。</p></div>
            <div class="case-card"><span class="case-label">Well-explained · HCC1395</span><h3>{html.escape(case['region'])}</h3><p>k={builder.number(case['k'], source_key='region')}、HP multiplicity={builder.number(case['hp_multiplicity'], source_key='region')}；兩個 mutation-bearing family 都各有 {builder.number(case['families'][0]['n_trees'], source_key='region')} 棵候選樹，且各含 {builder.number(case['families'][0]['n_hidden'], source_key='region')} 個 parsimony-inferred mutation-state nodes（不是 hidden clones）。兩 family 的 full populations 都為 {builder.number(sum(item['n_full_pops'] for item in case['families']), source_key='region')}，partial pattern types 為 {builder.number(case['families'][0]['n_partial'], source_key='region')}／{builder.number(case['families'][1]['n_partial'], source_key='region')}，仍主要由 partial evidence 約束。</p></div>
          </section>
        </article>

        <article class="reading">
          <section data-contract-section="scope-data-and-metric-definitions" class="narrative">
            <h2>Scope 與單位先於解釋：chrX 必須另建 ploidy／PAR／CN 契約</h2>
            <p>所有圖表都必須先讀「單位」再讀百分比；跨層的數字不能直接相除或串成單一 biological funnel。</p>
            <details class="technical-appendix compact-appendix"><summary><span>展開跨層 unit dictionary</span><small>U、W、k、read、HP、tree 與 auxiliary evidence</small></summary><div class="appendix-body"><div class="card table-card"><div class="card-head"><h3>跨層 unit dictionary</h3><p>每一層同時列正確讀法與禁止讀法。</p></div>{unit_table}</div></div></details>
            <p>目前把 chr1–22 與 chrX/Y 分開是合理的 scope control，但不是 solver 無法處理 chrX。VCF 規格要求 male non-PAR X 使用 haploid genotype；GRCh38 的 PAR1/PAR2 又是 X/Y 共享區域。癌症還可能有 X loss、gain、LOH，而 XX 的 X-inactivation 會產生既有甲基化差異。</p>
            <div class="definition-grid">
              <div><strong>Autosome</strong><span>正常基準通常為兩份；癌症 CN 仍可改變。</span></div>
              <div><strong>XX non-PAR X</strong><span>DNA 通常兩份；Xa/Xi 甲基基準不同。</span></div>
              <div><strong>XY non-PAR X</strong><span>通常一份；不能沿用兩個 homolog 的 HP 分母。</span></div>
              <div><strong>PAR1/PAR2</strong><span>XY 中由 X+Y 提供兩份同源序列；mapping 與 ploidy 另算。</span></div>
            </div>
            <p>HCC1395 的 SEQC2 high-confidence BED 與 CN BED 都存在，但其中 chrX interval 皆為 {builder.number(0, source_key='seqc2')}；這代表 benchmark 排除，不代表 chrX 生物缺失或 neutral。HCC1395/HCC1395BL 有複雜 karyotype，不能只由 donor sex 推定 copy complement；current manifest 也沒有任何 sample 的 sex、karyotype、PAR policy 或 X-specific CN contract。若不修改 contract 就直接納入 chrX，現行 CN lookup 在未命中時可能誤回 neutral。故正確方向是新增 <code>sex_chromosome_exploratory</code> 分層，而不是把 chrX 永久丟棄。</p>
            <p class="external-links">科學依據：<a href="https://samtools.github.io/hts-specs/VCFv4.5.pdf">GA4GH VCF 4.5</a> · <a href="https://www.ncbi.nlm.nih.gov/grc/human">GRC GRCh38 PAR</a> · <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360035890511-Assigning-per-sample-genotypes-HaplotypeCaller">GATK ploidy</a> · <a href="https://sites.google.com/view/seqc2/home/identify-high-confidence-regions">SEQC2 high-confidence scope</a> · <a href="https://www.nature.com/articles/s41597-021-01077-5">HCC1395/HCC1395BL karyotype</a></p>
          </section>

          <section data-contract-section="methodology" class="narrative">
            <h2>方法是「same-read observation → family stratification → candidate-tree enumeration」</h2>
            <ol class="method-steps">
              <li><strong>Historical chart universe：</strong>paired ClairS FILTER=PASS biallelic sSNV；這是舊工程 snapshot，不是 biological truth，也不是 04:50 決策後的 canonical tree backbone。</li>
              <li><strong>Current canonical universe：</strong>raw ClairS all → normalized raw-all LongPhase-S input → <code>_sc.vcf</code> all → <code>_sc.vcf</code> FILTER=PASS tree input；record keys、FILTER transitions 與 disposition 必須守恆。</li>
              <li><strong>Region：</strong>chr1–22 上相鄰 gap≤{builder.number(method['TIER_R'] / 1000, 'decimal1', 'method')} kb 的 connected component；總 span 可大於該 gap threshold。</li>
              <li><strong>Engineering cap：</strong>每個 component 最多取最密集連續 {builder.number(method['MAX_SNV'], source_key='method')} 個 sSNV；其餘仍記入漏斗。</li>
              <li><strong>Read observation：</strong>ONT molecule 轉為 R/A/X；MAPQ≥{builder.number(method['MAPQ_MIN'], source_key='method')}、BASEQ≥{builder.number(method['BASEQ_MIN'], source_key='method')}、MINREAD={builder.number(method['MINREAD'], source_key='method')}。</li>
              <li><strong>HP family：</strong>HP1/HP2 分開建 primary mutation-state tree；H3?、reference-only、unphased 另列。</li>
              <li><strong>Solver：</strong>枚舉 minimal-hidden、root-connected、unit-flip candidate trees；capped 代表搜索未完成。</li>
              <li><strong>Auxiliary evidence：</strong>CN 只解釋 recurrence-like model outcomes；candidate-tree read-AF 即使產出也只能給既有候選相對權重，不能當 caller VAF、family read fraction、CCF 或 topology confirmation；甲基介面只保留 bounded residual role，本 run 沒有甲基結果。</li>
            </ol>
          </section>

          <section id="scope-limitations" data-contract-section="limitations-uncertainty-and-robustness-checks" class="narrative">
            <h2>目前最重要的限制不是圖不夠，而是證據層尚未閉合</h2>
            <div class="status-matrix">
              <div><span class="pill">已驗</span><strong>Baseline artifact integrity</strong><p>{builder.number(verification['n_pass'], source_key='verification')}/{builder.number(verification['dataset_count'], source_key='verification')} outputs；root 與 {builder.number(verification['sample_hash_manifests_pass'], source_key='verification')} sample hash manifests 通過。</p></div>
              <div><span class="pill warn">缺口</span><strong>HCC1395 denominator</strong><p>{builder.number(len(hcc_r['missing_from_region_view']), source_key='derived')} region／{builder.number(hcc_r['site_sum_delta'], source_key='derived')} sSNV 未在 region-view emit。</p></div>
              <div><span class="pill warn">未產出</span><strong>candidate-tree read-AF ordering</strong><p>{builder.number(auxiliary['read_af']['dataset_artifacts'], source_key='read_af')}/{builder.number(auxiliary['read_af']['expected_datasets'], source_key='read_af')} dataset artifacts；與 caller VAF、family read fraction、CCF 分開，不能聲稱已矯正或確認樹。</p></div>
              <div><span class="pill warn">未評估</span><strong>Methylation L3</strong><p>{builder.number(auxiliary['methylation']['evaluated_datasets'], source_key='layered')}/{builder.number(auxiliary['methylation']['expected_datasets'], source_key='layered')} datasets、{builder.number(auxiliary['methylation']['evaluated_units'], source_key='layered')} units evaluated；7/7 都是 not_evaluated，不能當成陰性、0 或 confirmation。</p></div>
              <div><span class="pill warn">未鎖定</span><strong>Historical code snapshot</strong><p>{builder.number(verification['code_hash_manifest_pass'], source_key='code_hash')}/{builder.number(verification['code_hash_manifest_total'], source_key='code_hash')} entries 目前仍符合 <code>code.sha256</code>；輸出 hashes 通過不等於 producer code 可完整重現。</p></div>
              <div><span class="pill warn">外部資料</span><strong>Biological truth</strong><p>缺 single-cell／multi-region orthogonal truth；bulk molecule 不等於 cell。</p></div>
              <div><span class="pill warn">HP census</span><strong>HP4 contract 未獨立閉合</strong><p>本 snapshot 明列 HP1/HP2、H3? 與 none/unphased；現行 driver 沒有獨立 HP4 denominator，未知 tag 不可默認為生物學 unphased。</p></div>
              <div><span class="pill warn">{production_state}</span><strong>Normalized raw-all production</strong><p>{builder.number(production['n_pass'], source_key='production')}/{builder.number(production['expected_datasets'], source_key='production')} PASS；active={html.escape('、'.join(production_active) or '無')}，aggregate gate 未通過前不能 promotion。</p></div>
            </div>
            <details class="provenance"><summary>輸入、命令與可重現性資料</summary>
              <p><strong>Engineering baseline：</strong><code>{html.escape(snapshot['paths']['baseline_run_root'])}</code></p>
              <p><strong>Normalized raw-all production：</strong><code>{html.escape(snapshot['paths']['production_run_root'])}</code></p>
              <p><strong>Input manifest：</strong><a href="{html.escape(Path(snapshot['paths']['baseline_manifest']).as_uri())}" aria-label="開啟 historical input manifest JSON">開啟歷史輸入清單</a></p>
              <p><strong>核心驗證：</strong><code>sha256sum -c verification.sha256</code>；每個 sample 執行 <code>sha256sum -c output.sha256</code>。</p>
              <p><strong>Report builder：</strong><code>{html.escape(snapshot['paths']['builder'])}</code></p>
            </details>
          </section>

          <section data-contract-section="recommended-next-steps" class="narrative">
            <h2>下一步按 gate 推進，不再把缺失資料用內部指標補掉</h2>
            <ol>
              <li>讓 normalized raw-all LongPhase-S production 完成 7/7 validation 與 aggregate `_SUCCESS`，鎖定 FILTER transitions、input/output hashes 與真正 consumed manifest。</li>
              <li>完成 direct clean tagged-BAM vs original BAM+sidecar 的 bounded exact-equivalence 與 adversarial fail-closed tests。</li>
              <li>修正 HCC1395 retained→region-view 守恆 gate，讓 1 region／4 sSNV 缺口能被 final verifier 阻擋或明確分類。</li>
              <li>另建 chrX/PAR exploratory contract；在 manifest 加 sex-chromosome complement、tumor/normal local CN、PAR 與 benchmark BED。</li>
              <li>若要升級 clone claim，取得 single-cell 或 multi-region truth；若要稱 CCF，需 purity + allele-specific integer CN + mutation multiplicity。</li>
            </ol>
          </section>

          <section data-contract-section="further-questions" class="narrative">
            <h2>仍會改變結論的問題</h2>
            <ul>
              <li>H2009 的高 solver-cap rate 是否與較高 k、partial coverage 或 pattern complexity 相關？完成 solver diagnostics 後，再觀察 CN-stratum association；CN 不進 L1 solver。</li>
              <li>HCC1395 chrX 大量 PASS records 在 paired caller、normal chrX copy state 與 benchmark scope之間如何分解？</li>
              <li>哪些模型內唯一 region 能在 single cells 中觀察到相同 mutation-state coexistence？</li>
            </ul>
          </section>

          <section class="caveat"><strong>Interpretation ceiling.</strong> 本報告是完整的內部數據觀察與方法稽核快照，不是 biological validation release。任何圖中的「determined」「multi-HP」「shape unique」都必須保留其模型與分母限定。</section>
        </article>
        <footer class="report-footer">Run <code>{html.escape(Path(snapshot['paths']['baseline_run_root']).name)}</code> · primary scope chr1–22 · report state PARTIAL/PROBE · output hash recorded in companion validation JSON · claim ceiling: regional mutation-state candidate sets, not confirmed clones/subclones.</footer>
      </main>
    </div>
    <!-- DATA_ANALYTICS_HTML_REPORT_RUNTIME -->
    '''

    template = shell_template.read_text(encoding="utf-8")
    extra_css = r'''
    :root { --muted:#616161; --warning:#9b3900; --observed:#0a67b5; --derived:#1f7a58; --model:#9a5a00; --unavailable:#5f6368; }
    @media(prefers-color-scheme:dark){:root{--muted:#c2c2c2;--warning:#ffb48c;--observed:#73baff;--derived:#80d7ae;--model:#ffc46b;--unavailable:#d0d0d0}}
    body { font-size:15px; }
    p, li { line-height:1.62; }
    h3 { font-size:16px; line-height:1.45; }
    h4 { margin:26px 0 10px; font-size:15px; }
    .skip-link { position:fixed; z-index:2000; top:8px; left:8px; width:1px; height:1px; padding:0; overflow:hidden; clip-path:inset(50%); border-radius:10px; background:var(--text); color:var(--surface); }
    .skip-link:focus { width:auto; height:auto; padding:10px 14px; overflow:visible; clip-path:none; outline:3px solid var(--purple); }
    .report-hero { position: relative; }
    .validation-ribbon { display:inline-flex; margin:0 0 18px; padding:7px 10px; border:1px solid color-mix(in srgb,var(--warning) 55%,var(--border)); border-radius:999px; color:var(--warning); background:color-mix(in srgb,var(--warning) 8%,var(--surface)); font-size:11px; font-weight:800; letter-spacing:.045em; }
    .blossom { position: absolute; top: 2px; right: 2px; width: 74px; height: 74px; opacity: .5; }
    .blossom i { position: absolute; left: 31px; top: 8px; width: 12px; height: 30px; border: 1px solid var(--border-strong); border-radius: 999px 999px 40% 40%; transform-origin: 6px 29px; background: color-mix(in srgb, var(--blue) 10%, transparent); }
    .blossom i:nth-child(2){transform:rotate(72deg)} .blossom i:nth-child(3){transform:rotate(144deg)} .blossom i:nth-child(4){transform:rotate(216deg)} .blossom i:nth-child(5){transform:rotate(288deg)}
    .status-no-go { color: var(--warning); }
    .sourced-value { font-variant-numeric:tabular-nums; text-decoration:none; }
    .data-drawer,.technical-appendix { margin:18px 0 34px; overflow:hidden; border:1px solid var(--border); border-radius:18px; background:var(--surface); }
    .data-drawer > summary,.technical-appendix > summary { display:flex; cursor:pointer; list-style:none; align-items:center; justify-content:space-between; gap:16px; padding:15px 18px; color:var(--text); font-weight:750; }
    .data-drawer > summary::-webkit-details-marker,.technical-appendix > summary::-webkit-details-marker { display:none; }
    .data-drawer > summary::after,.technical-appendix > summary::after { content:"展開"; color:var(--muted); font-size:11px; font-weight:650; }
    .data-drawer[open] > summary::after,.technical-appendix[open] > summary::after { content:"收合"; }
    .data-drawer[open] > summary,.technical-appendix[open] > summary { border-bottom:1px solid var(--border); background:var(--surface-tertiary); }
    .technical-appendix summary span { display:block; }
    .technical-appendix summary small { margin-left:auto; color:var(--muted); font-size:11px; font-weight:500; text-align:right; }
    .appendix-body { padding:18px; }
    .appendix-body > :first-child { margin-top:0; }
    .data-link-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:10px; padding:14px; }
    .data-link-grid a { display:grid; gap:3px; min-height:68px; padding:12px 14px; border:1px solid var(--border); border-radius:12px; color:var(--text); text-decoration:none; }
    .data-link-grid a:hover,.data-link-grid a:focus-visible { border-color:var(--blue); outline:2px solid var(--blue); outline-offset:1px; }
    .data-link-grid span { color:var(--muted); font-size:11px; }
    .current-gate { display:grid; grid-template-columns:minmax(0,1fr) minmax(240px,320px); gap:28px; align-items:center; margin-top:22px; padding:22px 24px; border:1px solid color-mix(in srgb,var(--blue) 35%,var(--border)); border-radius:20px; background:color-mix(in srgb,var(--blue) 5%,var(--surface)); }
    .current-gate h2 { margin:7px 0; font-size:21px; line-height:1.35; }
    .current-gate p { margin:0; color:var(--secondary); font-size:12px; }
    .gate-label { color:var(--blue); font-size:10px; font-weight:800; letter-spacing:.09em; }
    .gate-progress { display:grid; gap:8px; color:var(--secondary); font-size:11px; }
    .gate-progress progress { width:100%; height:12px; accent-color:var(--blue); }
    .report-chart { margin: 42px 0 12px; }
    .report-chart .chart-wrap { min-height: 360px; }
    .exact-data { margin: 0 0 42px; border: 1px solid var(--border); border-radius: 14px; background: var(--surface); }
    .exact-data summary, .provenance summary { cursor: pointer; padding: 12px 16px; color: var(--secondary); font-weight: 600; }
    .exact-data[open] summary, .provenance[open] summary { border-bottom: 1px solid var(--border); }
    th.left { color:var(--text); font-weight:600; text-align:left; }
    caption { position:absolute; width:1px; height:1px; padding:0; margin:-1px; overflow:hidden; clip:rect(0,0,0,0); white-space:nowrap; border:0; }
    .table-scroll { position:relative; scrollbar-color:var(--border-strong) transparent; }
    .table-scroll:focus-visible { outline:3px solid var(--purple); outline-offset:2px; }
    tbody tr:last-child td, tbody tr:last-child th { border-bottom:0; }
    .mini-toc { position:sticky; z-index:8; top:0; display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin-top:18px; padding:12px 14px; border:1px solid var(--border); border-radius:16px; background:color-mix(in srgb,var(--surface) 94%,transparent); backdrop-filter:blur(12px); }
    .mini-toc strong { margin-right:8px; }
    .mini-toc a { padding:5px 8px; border-radius:8px; color:var(--secondary); font-size:12px; text-decoration:none; }
    .mini-toc a:hover,.mini-toc a:focus-visible { background:var(--surface-tertiary); color:var(--text); outline:2px solid var(--purple); }
    .evidence-badge,.historical-chip { display:inline-flex; align-items:center; margin:0 6px 8px 0; padding:4px 9px; border:1px solid currentColor; border-radius:999px; font-size:11px; font-weight:750; letter-spacing:.025em; }
    .evidence-badge.observed { color:var(--observed); }.evidence-badge.model { color:var(--model); }.evidence-badge.unavailable { color:var(--unavailable); border-style:dashed; }
    .historical-chip { margin:10px 0 0; color:var(--warning); border-style:dashed; }
    .concept-grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:10px; margin:16px 0 46px; }
    .concept-detail,.dataset-detail { border:1px solid var(--border); border-radius:16px; background:var(--surface); }
    .concept-detail summary,.dataset-detail summary { cursor:pointer; padding:15px 17px; color:var(--text); font-weight:700; }
    .concept-detail p { margin:0; padding:0 17px 17px; color:var(--secondary); }
    .dataset-list { display:grid; gap:14px; margin-bottom:50px; }
    .dataset-detail { scroll-margin-top:80px; overflow:hidden; }
    .dataset-detail[open] summary { border-bottom:1px solid var(--border); background:var(--surface-tertiary); }
    .dataset-detail > :not(summary) { margin-left:18px; margin-right:18px; }
    .dataset-meta { display:flex; flex-wrap:wrap; gap:8px 18px; padding:16px 0 4px; color:var(--secondary); font-size:12px; }
    .legal-reading { margin-bottom:20px!important; padding:14px 16px; border-left:4px solid var(--model); background:var(--surface-tertiary); color:var(--secondary); }
    .na-panel { display:grid; grid-template-columns:auto 1fr; gap:16px; align-items:start; margin:28px 0 52px; padding:20px; border:1px dashed var(--unavailable); border-radius:18px; background:var(--surface-tertiary); }
    .na-panel p { margin:6px 0 0; color:var(--secondary); }
    .svg-explainer { margin:18px 0 38px; padding:18px; overflow-x:auto; border:1px solid var(--border); border-radius:22px; background:linear-gradient(145deg,var(--surface),var(--surface-tertiary)); }
    .svg-explainer figcaption { margin-bottom:12px; color:var(--secondary); }
    .svg-explainer svg { display:block; width:100%; min-width:760px; height:auto; font-family:var(--ds-font); }
    .svg-node rect,.topo-panel rect { fill:var(--surface); stroke-width:1.5; }
    .svg-node.observed rect { stroke:var(--observed); }.svg-node.derived rect { stroke:var(--derived); }.svg-node.model rect { stroke:var(--model); }
    .svg-node text,.topo-panel text,.svg-legend text,.aux-rail text { fill:var(--secondary); font-size:12px; }
    .svg-step { fill:var(--muted)!important; font-size:10px!important; font-weight:800; letter-spacing:.08em; }
    .svg-title { fill:var(--text)!important; font-size:22px!important; font-weight:750; }.svg-title-sm { fill:var(--text)!important; font-size:13px!important; font-weight:700; }
    .svg-arrow { stroke:var(--border-strong); stroke-width:2; marker-end:url(#arrow); }.svg-arrow-head { fill:var(--secondary); }.unit-divider { stroke:var(--warning); stroke-width:2; stroke-dasharray:5 5; }.svg-note { fill:var(--warning); font-size:12px; }
    .aux-rail rect { fill:color-mix(in srgb,var(--unavailable) 7%,var(--surface)); stroke:var(--unavailable); stroke-dasharray:6 5; }
    .legend-observed { fill:var(--observed); }.legend-derived { fill:var(--derived); }.legend-model { fill:var(--model); }
    .topo-panel.ok rect { stroke:var(--derived); }.topo-panel.warn rect { stroke:var(--model); }.topo-panel.impossible rect { stroke:var(--warning); stroke-dasharray:6 4; }
    .tree-edge { fill:none; stroke:var(--observed); stroke-width:3; }.tree-edge.ghost { stroke:var(--model); stroke-dasharray:6 4; }.topo-panel circle { fill:var(--surface); stroke:var(--observed); stroke-width:2; }.cross { stroke:var(--warning); stroke-width:5; }
    .topology-metrics { margin-top:10px; }.subtable-title { padding:16px 18px 4px; }
    .chart-scroll-hint { display:none; margin:8px 18px 0; color:var(--muted); font-size:12px; }
    .source-figure { position:relative; }
    .source-figure .card-head { padding-right:86px; }
    .chart-source { position:absolute; top:14px; right:16px; padding:7px 10px; border:1px solid var(--border); border-radius:999px; background:var(--surface); color:var(--secondary); font-size:11px; font-weight:700; text-decoration:none; }
    .chart-source:hover,.chart-source:focus-visible { border-color:var(--blue); color:var(--blue); outline:2px solid var(--blue); }
    [data-recharts-live] svg:focus-visible,.chart-wrap:focus-visible { outline:3px solid var(--purple); outline-offset:2px; }
    .case-grid, .status-matrix, .definition-grid, .flow-grid { display: grid; grid-template-columns: repeat(2,minmax(0,1fr)); gap: 12px; margin: 20px 0 46px; }
    .case-card, .status-matrix > div, .definition-grid > div, .flow-grid > div { padding: 18px; border: 1px solid var(--border); border-radius: 18px; background: var(--surface); }
    .case-card p, .status-matrix p, .definition-grid span, .flow-grid p { display: block; margin: 8px 0 0; color: var(--secondary); font-size: 12px; line-height: 19px; }
    .flow-grid > div > span { display:inline-grid; place-items:center; width:24px; height:24px; margin-right:8px; border-radius:50%; background:color-mix(in srgb,var(--blue) 12%,var(--surface)); color:var(--blue); font-size:11px; font-weight:800; }
    .case-label { color: var(--blue); font-size: 11px; font-weight: 700; letter-spacing: .05em; text-transform: uppercase; }
    .method-steps { padding-left: 22px; }
    .method-steps li { margin: 10px 0; color: var(--secondary); }
    .method-steps strong { color: var(--text); }
    .external-links { font-size: 12px; }
    a { color: var(--blue); }
    .provenance { margin: 20px 0; border: 1px solid var(--border); border-radius: 14px; }
    .provenance p { margin: 12px 16px; }
    .report-footer { max-width:var(--reading); margin:34px auto 0; padding:20px 0 8px; border-top:1px solid var(--border); color:var(--muted); font-size:11px; line-height:18px; }
    code { font-family: ui-monospace, SFMono-Regular, Consolas, monospace; overflow-wrap: anywhere; }
    @media(max-width:800px){
      .case-grid,.status-matrix,.definition-grid,.flow-grid,.concept-grid,.data-link-grid{grid-template-columns:1fr}.blossom{display:none}
      .topbar{align-items:flex-start;flex-wrap:wrap}.topbar .meta{width:100%;max-width:none;overflow:visible;text-overflow:clip;white-space:normal;padding-left:27px;line-height:1.35}
      .current-gate{grid-template-columns:1fr;margin-left:18px;margin-right:18px}.technical-appendix summary{align-items:flex-start;flex-wrap:wrap}.technical-appendix summary small{width:100%;margin:0;text-align:left}
      .mini-toc{position:static}.svg-explainer{margin-left:-18px;margin-right:-18px;border-radius:0;border-left:0;border-right:0}
      .report-chart .chart-wrap{overflow-x:auto;padding:12px 14px}.report-chart [data-recharts-chart]{min-width:760px}
      .chart-scroll-hint{display:block}.source-figure .card-head{padding-right:20px}
      .chart-source{position:static;display:inline-flex;margin:10px 18px 14px;padding:8px 12px}
      .na-panel{grid-template-columns:1fr}
    }
    @media(prefers-reduced-motion:reduce){*,*::before,*::after{scroll-behavior:auto!important;animation-duration:.01ms!important;animation-iteration-count:1!important;transition-duration:.01ms!important}}
    @media print{.mini-toc,.chart-source,.chart-scroll-hint,.data-drawer{display:none!important}.shell{box-shadow:none;border:0}.concept-detail>*,.dataset-detail>*{display:block!important}.exact-data{break-inside:avoid}.svg-explainer{break-inside:avoid}}
    '''
    template = template.replace('lang="en"', 'lang="zh-Hant"')
    template = template.replace("{{TITLE}}", "InterSubMod 分層重建：目前進度、歷史基線與證據邊界")
    template = template.replace("</style>", extra_css + "\n</style>", 1)
    body_start = template.index("<body>") + len("<body>")
    body_end = template.rindex("</body>")
    return template[:body_start] + main_content + template[body_end:]


def build_chart_map(output_html: Path, snapshot):
    run_root = Path(snapshot["paths"]["baseline_run_root"])
    production_root = Path(snapshot["paths"]["production_run_root"])
    return {
        "delivery_mode": "html",
        "audience": "technical",
        "palette_policy": "stacked composition uses relaxed multi-category; binary and single-measure charts use one/two roots plus neutrals",
        "charts": [
            {"id": "site-funnel-composition", "question": "Where do PASS sSNVs leave the solver?", "family": "composition", "type": "stacked bar", "denominator": "PASS sSNV universe by dataset"},
            {"id": "region-stage-counts", "question": "What remains at each W stage?", "family": "comparison", "type": "grouped bar", "denominator": "explicit nested region/group stage"},
            {"id": "k-distribution", "question": "What k enters the solver?", "family": "composition", "type": "100% stacked bar", "denominator": "W_ret by dataset"},
            {"id": "read-pattern-composition", "question": "How much evidence is full versus partial?", "family": "composition", "type": "100% stacked bar", "denominator": "retained pattern observations"},
            {"id": "hp-multiplicity", "question": "How many primary HP families occur per region?", "family": "composition", "type": "100% stacked bar", "denominator": "W_tree"},
            {"id": "hp-h3-strata", "question": "How do HP multiplicity and H3? auxiliary evidence intersect?", "family": "composition", "type": "100% stacked bar", "denominator": "W_tree"},
            {"id": "candidate-combinations-single", "question": "How many exact candidates remain in one-primary regions?", "family": "composition", "type": "100% stacked bar", "denominator": "one-primary regions"},
            {"id": "candidate-combinations-double", "question": "How many joint exact combinations remain in two-primary regions?", "family": "composition", "type": "100% stacked bar", "denominator": "two-primary regions"},
            {"id": "region-topology-classes", "question": "Which feasible exact/topology state remains per primary region?", "family": "composition", "type": "100% stacked bar", "denominator": "W_primary"},
            {"id": "hidden-node-classes", "question": "Do complete region candidates require hidden mutation-state nodes?", "family": "composition", "type": "100% stacked bar", "denominator": "W_primary"},
            {"id": "tree-outcomes", "question": "What is the L1 inference outcome?", "family": "composition", "type": "100% stacked bar", "denominator": "primary reconstruction units"},
            {"id": "determinacy-levels", "question": "How do exact-tree and shape uniqueness differ?", "family": "comparison", "type": "grouped bar", "denominator": "analysis-complete non-capped primary reconstruction units"},
            {"id": "solver-cap-rate", "question": "Which dataset has the highest incomplete-search burden?", "family": "ranking", "type": "horizontal bar", "denominator": "primary reconstruction units"},
            {"id": "sex-chromosome-share", "question": "How large is the non-autosomal caller branch?", "family": "ranking", "type": "horizontal bar", "denominator": "PASS sSNV universe"},
        ],
        "repeated_type_rationale": "Stacked bars intentionally keep the same denominator-composition grammar across sequential biological/engineering layers; each chart changes the unit and states it explicitly.",
        "omissions": {
            "mixed_unit_sankey": "site, region and reconstruction-unit denominators must not be connected as one flow",
            "raw_component_density": "raw component boundaries/spans are not stored; pre_cap_n divided by selected densest-8 span would be misleading",
            "read_af": "canonical output not generated",
            "methylation": "not evaluated",
        },
        "source_contracts": {
            "historical_scope": {
                "path": str(run_root / "VALIDATION_SCOPE.md"),
                "role": "claim boundary; 6/7 truth-BED-conditioned HP/PS inputs",
            },
            "verification": {
                "paths": [str(run_root / "verification_summary.json"), str(run_root / "verification.sha256")],
                "role": "historical artifact integrity and funnel census",
            },
            "sample_outputs": {
                "path_template": str(run_root / "samples/<sample>/{mlhp_part_*.json,layered_reconstruction_<sample>.json,layered_region_view_<sample>.json,output.sha256}"),
                "role": "read, region, HP, candidate-tree, L2 and L3 metrics",
            },
            "production": {
                "paths": [str(production_root / "input_manifest.json"), str(production_root / "run_status.tsv"), str(production_root / "verification_summary.json")],
                "role": "latest truth-unrestricted production probe status; missing aggregate remains pending",
            },
            "report_builder": {
                "path": snapshot["paths"]["builder"],
                "role": "deterministic recomputation and HTML generation",
            },
        },
        "claim_boundaries": [
            "Historical charts are not the current canonical LongPhase-S recalibrated PASS tree run.",
            "read-AF is not generated; methylation is not evaluated.",
            "Regional mutation-state candidates are not confirmed cell clones or subclones.",
            "chrX/Y are census-only until a sex/PAR/ploidy/CN contract exists.",
        ],
        "final_artifact": str(output_html),
    }


def build_snapshot(args):
    verification_path = args.run_root / "verification_summary.json"
    manifest_path = args.run_root / "input_manifest.json"
    verification = load_json(verification_path)
    manifest = load_json(manifest_path)
    method_contract = load_json(args.run_root / "params.json")["params"]
    verification_by_sample = {item["sample"]: item for item in verification["samples"]}
    manifest_by_sample = {item["sample"]: item for item in manifest["samples"]}
    samples = []
    for sample in SAMPLE_ORDER:
        samples.append(sample_metrics(manifest_by_sample[sample], verification_by_sample[sample], args.run_root))

    root_hash = verify_sha_manifest(args.run_root / "verification.sha256")
    sample_hash_results = [verify_sha_manifest(args.run_root / "samples" / sample / "output.sha256") for sample in SAMPLE_ORDER]
    code_hash = verify_sha_manifest(args.run_root / "code.sha256")
    production = read_production_status(args.production_root)
    chart_contract_path = args.embed_helper.parent / "chart_contract_v0_2_6.py"

    hcc_hc_bed_chr_x = None
    if args.hcc_hc_bed.is_file():
        with args.hcc_hc_bed.open(encoding="utf-8") as handle:
            hcc_hc_bed_chr_x = sum(1 for line in handle if line.split("\t", 1)[0] in {"chrX", "X"})
    hcc_cn_bed_chr_x = None
    if args.hcc_cn_bed.is_file():
        with args.hcc_cn_bed.open(encoding="utf-8") as handle:
            hcc_cn_bed_chr_x = sum(1 for line in handle if line.split("\t", 1)[0] in {"chrX", "X"})

    return {
        "schema_version": "1.0",
        "generated_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "report_scope": {
            "task_type": "B comprehensive observation report",
            "primary_genome": "chr1-22",
            "sex_chromosome_scope": "out-of-scope census + rationale; no topology reconstruction",
            "dataset_count": verification["dataset_count"],
            "biological_sample_count": verification["biological_sample_count"],
            "claim_ceiling": "regional mutation-state candidate trees; not confirmed cell clones",
        },
        "paths": {
            "baseline_run_root": str(args.run_root),
            "production_run_root": str(args.production_root),
            "baseline_manifest": str(manifest_path),
            "builder": str(Path(__file__).resolve()),
            "validation_scope": str(args.run_root / "VALIDATION_SCOPE.md"),
            "params": str(args.run_root / "params.json"),
        },
        "report_generation_provenance": {
            "builder": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
            "shell_template": {"path": str(args.shell_template), "sha256": sha256(args.shell_template)},
            "embed_helper": {"path": str(args.embed_helper), "sha256": sha256(args.embed_helper)},
            "chart_contract": {
                "path": str(chart_contract_path),
                "sha256": sha256(chart_contract_path) if chart_contract_path.is_file() else None,
            },
            "baseline_manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
            "baseline_params": {"path": str(args.run_root / "params.json"), "sha256": sha256(args.run_root / "params.json")},
        },
        "baseline_verification": {
            "schema_version": verification.get("schema_version"),
            "dataset_count": verification["dataset_count"],
            "biological_sample_count": verification["biological_sample_count"],
            "all_pass": verification["all_pass"],
            "n_pass": verification["n_pass"],
            "n_fail": verification["n_fail"],
            "root_hash_manifest_pass": root_hash["pass"],
            "sample_hash_manifests_pass": sum(1 for item in sample_hash_results if item["pass"]),
            "sample_hash_manifests_total": len(sample_hash_results),
            "code_hash_manifest_pass": sum(1 for item in code_hash["checked"] if item["pass"]),
            "code_hash_manifest_total": len(code_hash["checked"]),
            "root_hash_detail": root_hash,
            "sample_hash_detail": sample_hash_results,
            "code_hash_detail": code_hash,
            "manifest_sample_names": [item["sample"] for item in manifest["samples"]],
            "verification_sample_names": [item["sample"] for item in verification["samples"]],
            "scientific_role": "hash-verified historical engineering snapshot only; 6/7 truth-BED-conditioned HP/PS inputs",
        },
        "method_contract": method_contract,
        "auxiliary_evidence": {
            "read_af": {
                "status": "not_generated",
                "expected_datasets": verification["dataset_count"],
                "dataset_artifacts": 1 if (args.run_root / "read_af_tree_ordering.json").is_file() else 0,
                "evaluated_units": 0,
            },
            "methylation": {
                "status": "not_evaluated",
                "expected_datasets": verification["dataset_count"],
                "evaluated_datasets": sum(1 for sample in samples if sample["L3_methyl"].get("status") != "not_evaluated"),
                "evaluated_units": 0,
            },
        },
        "clean_production": production,
        "seqc2_chrX_contract": {
            "high_confidence_bed_exists": args.hcc_hc_bed.is_file(),
            "high_confidence_bed_chrX_intervals": hcc_hc_bed_chr_x,
            "cn_bed_exists": args.hcc_cn_bed.is_file(),
            "cn_bed_chrX_intervals": hcc_cn_bed_chr_x,
            "manifest_samples_with_sex_or_karyotype": sum(
                1 for item in manifest["samples"]
                if any(key in item for key in ("sex", "karyotype", "sex_chromosome_complement"))
            ),
        },
        "samples": samples,
    }


def validate_snapshot(snapshot):
    errors = []
    baseline = snapshot["baseline_verification"]
    if baseline["schema_version"] != "2.0":
        errors.append("historical observation profile requires verification schema 2.0")
    if baseline["all_pass"] is not True or baseline["n_pass"] != baseline["dataset_count"]:
        errors.append("engineering baseline is not 7/7 all-pass")
    if baseline["root_hash_manifest_pass"] is not True:
        errors.append("root verification.sha256 failed")
    if baseline["sample_hash_manifests_pass"] != baseline["sample_hash_manifests_total"]:
        errors.append("one or more sample output.sha256 manifests failed")
    if len(snapshot["samples"]) != 7:
        errors.append("sample count is not seven")
    if set(baseline["manifest_sample_names"]) != set(SAMPLE_ORDER) or set(baseline["verification_sample_names"]) != set(SAMPLE_ORDER):
        errors.append("manifest/verification sample set differs from the declared seven datasets")
    if len(set(baseline["manifest_sample_names"])) != 7 or len(set(baseline["verification_sample_names"])) != 7:
        errors.append("manifest/verification sample set contains duplicates")
    production = snapshot["clean_production"]
    if production["expected_datasets"] != 7 or len(production["manifest_sample_names"]) != 7:
        errors.append("clean production manifest does not declare exactly seven datasets")
    if production["state"] not in {"not_started", "in_progress", "failed", "complete"}:
        errors.append("clean production state is not a valid tri-state-plus-failure value")
    if production["state"] == "complete" and not (production["aggregate_complete"] and production["sample_status_complete"]):
        errors.append("production completion lacks exact sample-status and aggregate all-pass gates")
    if not snapshot["seqc2_chrX_contract"]["high_confidence_bed_exists"] or not snapshot["seqc2_chrX_contract"]["cn_bed_exists"]:
        errors.append("one or more required SEQC2 chrX-scope source BEDs are missing")
    if snapshot["auxiliary_evidence"]["read_af"]["dataset_artifacts"] != 0:
        errors.append("unexpected read-AF artifact found; report wording/profile must be upgraded")
    if snapshot["auxiliary_evidence"]["methylation"]["evaluated_datasets"] != 0:
        errors.append("unexpected evaluated methylation dataset found; report wording/profile must be upgraded")
    for sample in snapshot["samples"]:
        funnel = sample["site_funnel"]
        if funnel["universe"] != funnel["out_of_scope"] + funnel["singleton"] + funnel["cap_excluded"] + funnel["read_unsupported"] + funnel["retained"]:
            errors.append(f"{sample['sample']}: site funnel does not conserve")
        if sum(sample["k_distribution"]["counts"].values()) != sample["regions"]["W_ret"]:
            errors.append(f"{sample['sample']}: k distribution does not equal W_ret")
        units = sample["units"]
        if units["primary"] != units["determined"] + units["ambiguous"] + units["solver_capped"] + units["recurrence"]:
            errors.append(f"{sample['sample']}: primary L1 outcomes do not conserve")
        if sample["sex_chromosome_census"]["vcf_total_matches"] is not True:
            errors.append(f"{sample['sample']}: VCF per-contig total does not match verification universe")
        if sample["sex_chromosome_census"]["vcf_nonpass_snv"] != 0:
            errors.append(f"{sample['sample']}: historical somatic VCF contains non-PASS biallelic SNVs")
        if sample["L3_methyl"].get("status") != "not_evaluated":
            errors.append(f"{sample['sample']}: unexpected L3 status")
        region = sample["regions"]
        if region["W_all_pre"] != region["W_k1"] + region["W_k_gt1"] or region["W_k_gt1"] != region["W_pre"]:
            errors.append(f"{sample['sample']}: k=1/k>1 region bridge does not conserve")
        if region["W_tree"] != region["single_primary_HP"] + region["multi_HP"] + region["no_primary"]:
            errors.append(f"{sample['sample']}: HP multiplicity does not conserve W_tree")
        if region["mlhp_site_sum"] != sample["site_funnel"]["retained"]:
            errors.append(f"{sample['sample']}: MLHP site sum does not equal retained sSNV")
        if region["W_ret"] - region["W_tree"] != len(region["missing_from_region_view"]):
            errors.append(f"{sample['sample']}: W gap does not match explicitly missing region count")
        if region["site_sum_delta"] != sum(item["k"] for item in region["missing_from_region_view_detail"]):
            errors.append(f"{sample['sample']}: site gap does not match missing-region k sum")
        if sample["L2_cn"]["n_primary_recurrence_sent_to_cn"] != sample["units"]["recurrence"]:
            errors.append(f"{sample['sample']}: L2 primary count does not match model recurrence-required units")
        if sample["method_params"] != {key: snapshot["method_contract"][key] for key in sample["method_params"]}:
            errors.append(f"{sample['sample']}: method parameters differ from run-level params.json")
        topology = sample["region_candidate_topology"]
        if sum(topology["hp_h3"].values()) != region["W_tree"]:
            errors.append(f"{sample['sample']}: HP×H3 strata do not conserve W_tree")
        if sum(topology["C_bins"]["single"].values()) != region["single_primary_HP"]:
            errors.append(f"{sample['sample']}: single-primary C bins do not conserve")
        if sum(topology["C_bins"]["double"].values()) != region["multi_HP"]:
            errors.append(f"{sample['sample']}: double-primary C bins do not conserve")
        if sum(topology["C_bins"]["all"].values()) != region["W_primary"]:
            errors.append(f"{sample['sample']}: all-primary C bins do not conserve W_primary")
        if topology["complete_regions"] + topology["incomplete_regions"] != region["W_primary"]:
            errors.append(f"{sample['sample']}: complete/incomplete regions do not conserve W_primary")
        if sum(topology["topology_classes"].values()) != region["W_primary"]:
            errors.append(f"{sample['sample']}: topology classes do not conserve W_primary")
        if topology["topology_classes"]["impossible_exact_unique_topology_multiple"] != 0:
            errors.append(f"{sample['sample']}: mathematically impossible C=1/Topo>1 state found")
        if sum(topology["hidden_classes"].values()) != region["W_primary"]:
            errors.append(f"{sample['sample']}: hidden classes do not conserve W_primary")
        if topology["invalid_states"]:
            errors.append(f"{sample['sample']}: invalid region candidate/topology states found")
    hcc = next(sample for sample in snapshot["samples"] if sample["sample"] == "HCC1395")
    if not hcc["well_explained_case"] or len(hcc["well_explained_case"]["families"]) != 2:
        errors.append("HCC1395 well-explained case could not be read from the region-view artifact")
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--production-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--shell-template", required=True, type=Path)
    parser.add_argument("--embed-helper", required=True, type=Path)
    parser.add_argument("--hcc-hc-bed", required=True, type=Path)
    parser.add_argument("--hcc-cn-bed", required=True, type=Path)
    args = parser.parse_args()

    if "validated" in args.output_dir.parts or "in_progress" not in args.output_dir.parts:
        raise SystemExit("baseline observation reports must publish under docs/reports/in_progress, never validated")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_html = args.output_dir / f"{REPORT_STEM}.standalone.html"
    output_data = args.output_dir / DATA_SNAPSHOT_NAME
    output_notes = args.output_dir / SOURCE_NOTES_NAME
    output_validation = args.output_dir / VALIDATION_NAME

    snapshot = build_snapshot(args)
    errors = validate_snapshot(snapshot)
    payload = chart_payload(snapshot["samples"])
    chart_contract = validate_official_chart_contract(
        payload,
        Path(snapshot["report_generation_provenance"]["chart_contract"]["path"]),
    )
    errors.extend(chart_contract["errors"])
    allowed_palette_kinds = {"categorical", "sequential", "diverging", "semantic", "identity"}
    chart_ids = [chart["id"] for chart in payload["charts"]]
    if len(chart_ids) != len(set(chart_ids)):
        errors.append("chart IDs are not unique")
    for chart in payload["charts"]:
        spec = chart["dataset"]["chart_spec"]
        if spec.get("palette", {}).get("kind") not in allowed_palette_kinds:
            errors.append(f"{chart['id']}: unsupported palette kind")
        for key in ("intent", "question", "rationale"):
            if not spec.get(key):
                errors.append(f"{chart['id']}: missing chart-spec {key}")
    report_shell = build_report_html(snapshot, args.shell_template)
    chart_map = build_chart_map(output_html, snapshot)
    dump_json(output_data, snapshot)
    dump_json(output_notes, chart_map)

    candidate_html = None
    with tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        shell_path = directory / "report-shell.html"
        payload_path = directory / "report-payload.json"
        embedded_path = directory / "embedded-report.html"
        shell_path.write_text(report_shell, encoding="utf-8")
        dump_json(payload_path, payload)
        process = subprocess.run(
            [sys.executable, str(args.embed_helper), "--input", str(shell_path), "--payload", str(payload_path), "--output", str(embedded_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        if process.returncode:
            errors.append("Recharts runtime embedding failed: " + process.stderr.strip())
        else:
            embedded_html = embedded_path.read_text(encoding="utf-8")
            candidate_html = "\n".join(line.rstrip() for line in embedded_html.splitlines()) + "\n"

    if candidate_html is not None:
        ids = re.findall(r'\bid="([^"]+)"', candidate_html)
        duplicate_ids = sorted(key for key, value in Counter(ids).items() if value > 1)
        if duplicate_ids:
            errors.append(f"duplicate HTML IDs: {len(duplicate_ids)}")
        if candidate_html.count("data-recharts-chart=") != len(payload["charts"]):
            errors.append("Recharts host count differs from payload chart count")
        if len(re.findall(r'\sdata-recharts-fallback(?:\s|>)', candidate_html)) != len(payload["charts"]):
            errors.append("static fallback count differs from payload chart count")
        if candidate_html.count('class="card source-figure report-chart"') != len(payload["charts"]):
            errors.append("source-figure count differs from payload chart count")
        if candidate_html.count('class="source-link chart-source"') != len(payload["charts"]):
            errors.append("chart-level hidden source-link count differs from payload chart count")
        for companion_name in (DATA_SNAPSHOT_NAME, SOURCE_NOTES_NAME, VALIDATION_NAME, BROWSER_QA_NAME):
            if f'href="{companion_name}"' not in candidate_html:
                errors.append(f"missing collapsed companion-data link: {companion_name}")
        if "DATA_ANALYTICS_HTML_REPORT_RUNTIME" in candidate_html:
            errors.append("HTML runtime marker remains after embedding")
        if re.search(r'<script[^>]+src=["\']https?://', candidate_html, re.I) or re.search(r'<link[^>]+href=["\']https?://', candidate_html, re.I):
            errors.append("remote JavaScript or stylesheet dependency found")
        figure_tags = re.findall(r"<figure\b[^>]*>", candidate_html, re.I)
        if any("aria-labelledby=" not in tag for tag in figure_tags):
            errors.append("one or more figures lack an accessible name")
        table_count = len(re.findall(r"<table\b", candidate_html, re.I))
        caption_count = len(re.findall(r"<caption\b", candidate_html, re.I))
        if table_count != caption_count:
            errors.append("one or more tables lack a caption")
        th_tags = re.findall(r"<th\b[^>]*>", candidate_html, re.I)
        if any("scope=" not in tag for tag in th_tags):
            errors.append("one or more table headers lack scope")
        summary_labels = [
            re.sub(r"<[^>]+>", "", value).strip()
            for value in re.findall(r"<summary\b[^>]*>(.*?)</summary>", candidate_html, re.I | re.S)
        ]
        duplicate_summaries = sorted(key for key, value in Counter(summary_labels).items() if value > 1)
        if duplicate_summaries:
            errors.append(f"duplicate details summary labels: {len(duplicate_summaries)}")
        if 'class="skip-link" href="#report-main"' not in candidate_html or 'id="report-main"' not in candidate_html:
            errors.append("accessible skip link/main target is missing")
        if candidate_html.count('class="source-tooltip inline-source"'):
            errors.append("inline numeric source tooltips must be removed from the reader-facing number flow")
        if candidate_html.count('class="sourced-value"') < 1:
            errors.append("machine-readable sourced-value annotations are missing")
        for section in (
            "title", "technical-summary", "key-findings", "scope-data-and-metric-definitions", "methodology",
            "limitations-uncertainty-and-robustness-checks", "recommended-next-steps", "further-questions",
        ):
            if f'data-contract-section="{section}"' not in candidate_html:
                errors.append(f"missing required technical report section: {section}")
        for required_text in ("PARTIAL SCIENTIFIC VALIDATION", "PUBLICATION NO-GO", "read-AF", "not_evaluated", "S → W → HP → C → Topo", "C=1, Topo=1", "NOT EVALUATED", "7,928", "7,927", "25,639", "25,635"):
            if required_text not in candidate_html:
                errors.append(f"missing required report evidence text: {required_text}")

    if not errors and candidate_html is not None:
        output_html.write_text(candidate_html, encoding="utf-8")

    validation = {
        "generated_at": snapshot["generated_at"],
        "status": "pass" if not errors else "fail",
        "artifact_validation_status": "PASS" if not errors else "FAIL",
        "scientific_release_status": "NO_GO",
        "scientific_release_reasons": [
            "historical 6/7 truth-BED-conditioned HP/PS inputs",
            f"historical code.sha256 currently matches {snapshot['baseline_verification']['code_hash_manifest_pass']}/{snapshot['baseline_verification']['code_hash_manifest_total']} entries",
            f"normalized raw-all production state={snapshot['clean_production']['state']}; pass={snapshot['clean_production']['n_pass']}/{snapshot['clean_production']['expected_datasets']}; aggregate_complete={snapshot['clean_production']['aggregate_complete']}",
            "candidate-tree read-AF not generated and L3 methylation not evaluated",
            f"HCC1395 historical region-view gap: {len(next(sample for sample in snapshot['samples'] if sample['sample'] == 'HCC1395')['regions']['missing_from_region_view'])} region(s) / {next(sample for sample in snapshot['samples'] if sample['sample'] == 'HCC1395')['regions']['site_sum_delta']} sSNV(s)",
            "no single-cell or multi-region orthogonal truth",
        ],
        "errors": errors,
        "output_html": str(output_html),
        "output_html_sha256": sha256(output_html) if not errors and output_html.is_file() else None,
        "output_data": str(output_data),
        "output_data_sha256": sha256(output_data),
        "chart_count": len(payload["charts"]),
        "official_chart_contract_status": "PASS" if not chart_contract["errors"] else "FAIL",
        "official_chart_contract_path": chart_contract["path"],
        "official_chart_contract_errors": chart_contract["errors"],
        "official_chart_contract_warnings": chart_contract["warnings"],
        "source_tooltip_count": sum(
            1 for classes in re.findall(r'class="([^"]*)"', candidate_html)
            if "source-tooltip" in classes.split()
        ) if candidate_html else 0,
        "inline_source_tooltip_count": candidate_html.count('class="source-tooltip inline-source"') if candidate_html else 0,
        "focusable_chart_source_control_count": candidate_html.count('class="source-link chart-source"') if candidate_html else 0,
        "hidden_companion_data_link_count": sum(
            candidate_html.count(f'href="{name}"')
            for name in (DATA_SNAPSHOT_NAME, SOURCE_NOTES_NAME, VALIDATION_NAME, BROWSER_QA_NAME)
        ) if candidate_html else 0,
        "sourced_value_annotation_count": candidate_html.count('class="sourced-value"') if candidate_html else 0,
        "recharts_host_count": candidate_html.count("data-recharts-chart=") if candidate_html else 0,
        "fallback_table_count": len(re.findall(r'\sdata-recharts-fallback(?:\s|>)', candidate_html)) if candidate_html else 0,
        "duplicate_id_count": len(duplicate_ids) if candidate_html is not None else None,
        "unnamed_figure_count": sum(
            1 for tag in re.findall(r"<figure\b[^>]*>", candidate_html, re.I) if "aria-labelledby=" not in tag
        ) if candidate_html else None,
        "tables_without_caption": (
            len(re.findall(r"<table\b", candidate_html, re.I)) - len(re.findall(r"<caption\b", candidate_html, re.I))
        ) if candidate_html else None,
        "table_headers_without_scope": sum(
            1 for tag in re.findall(r"<th\b[^>]*>", candidate_html, re.I) if "scope=" not in tag
        ) if candidate_html else None,
        "skip_link_present": 'class="skip-link" href="#report-main"' in candidate_html if candidate_html else False,
        "runtime_marker_remaining": "DATA_ANALYTICS_HTML_REPORT_RUNTIME" in candidate_html if candidate_html else True,
    }
    dump_json(output_validation, validation)
    if errors:
        raise SystemExit("report validation failed: " + "; ".join(errors))
    print(f"REPORT HTML -> {output_html}")
    print(f"DATA SNAPSHOT -> {output_data}")
    print(f"SOURCE NOTES -> {output_notes}")
    print(f"VALIDATION -> {output_validation}")


if __name__ == "__main__":
    main()
