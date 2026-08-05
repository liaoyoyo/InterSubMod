#!/usr/bin/env python3
"""Rebuild fixed-denominator HCC1395 regional/site/topology precision strata.

The 5,720 exact-coordinate regions complete in both HCC1395 technical datasets
are the immutable spine.  Every joined endpoint is recomputed from row-level
artifacts; the integrated context is used only as an independent conservation
gate.  Missing columns or non-unique keys fail closed.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import sys
import tempfile
from typing import Any, Callable, Iterable, Mapping, Sequence


SCRIPT = Path(__file__).resolve()
TOPIC = SCRIPT.parents[1]
REPO = SCRIPT.parents[3]
FIXED_N = 5720
SAMPLE_A = "HCC1395"
SAMPLE_B = "HCC1395_DORADO"
K_STRATA = ("ALL", "k=1", "k=2", "k=3", "k>=4")
COARSE_STATES = (
    "no_within_hp_relation",
    "sister_only",
    "direct_only",
    "sister_and_direct",
    "topology_multiple_unresolved",
)


@dataclass(frozen=True)
class MetricSpec:
    group: str
    metric_id: str
    label: str
    gate: Callable[[Mapping[str, Any]], bool]
    success: Callable[[Mapping[str, Any]], bool]
    source_fields: str
    claim_ceiling: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coarse-dir",
        type=Path,
        default=REPO / "research" / "20260712_hcc1395_pair_coarse_topology_gene_drug_validation" / "data",
    )
    parser.add_argument(
        "--site-dir",
        type=Path,
        default=REPO / "research" / "20260712_hcc1395_pair_site_topology_containment_validation" / "data",
    )
    parser.add_argument(
        "--integrated-context",
        type=Path,
        default=TOPIC / "results" / "integrated_topology_context_v1.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=TOPIC / "results" / "regional_precision_validation_v1",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_tsv(path: Path, required: Iterable[str]) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = set(required).difference(fields)
        if missing:
            raise ValueError(f"{path} lacks required fields: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} is empty")
    return rows, fields


def boolean(value: str, *, allow_empty: bool = False) -> bool | None:
    if value == "True":
        return True
    if value == "False":
        return False
    if allow_empty and value == "":
        return None
    raise ValueError(f"Invalid boolean token: {value!r}")


def integer(value: str) -> int:
    converted = int(value)
    if str(converted) != value and not (value.startswith("+") and str(converted) == value[1:]):
        raise ValueError(f"Non-canonical integer token: {value!r}")
    return converted


def k_stratum(k: int) -> str:
    if k == 1:
        return "k=1"
    if k == 2:
        return "k=2"
    if k == 3:
        return "k=3"
    if k >= 4:
        return "k>=4"
    raise ValueError(f"caller_shared_k must be positive, got {k}")


def unlabeled_signature(value: str) -> tuple[str, ...]:
    if not value:
        return tuple()
    components = []
    for token in value.split("|"):
        if "=" not in token:
            raise ValueError(f"Malformed HP signature: {value}")
        components.append(token.split("=", 1)[1])
    return tuple(sorted(components))


def unique_index(rows: Sequence[Mapping[str, str]], key: str, label: str) -> dict[str, Mapping[str, str]]:
    result: dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = row[key]
        if value in result:
            raise ValueError(f"Duplicate {label} key {key}={value}")
        result[value] = row
    return result


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def write_json(path: Path, value: Any) -> None:
    atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n")


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, delete=False) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def write_gzip_tsv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=path.parent, delete=False) as raw_handle:
        temporary = Path(raw_handle.name)
    try:
        with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="raise", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def main() -> int:
    args = parse_args()
    paths = {
        "pair_matches": args.coarse_dir / "hcc1395_pair_matches.tsv",
        "vaf_pair_regions": args.coarse_dir / "hcc1395_vaf_pair_regions.tsv",
        "site_pair_outcomes": args.site_dir / "hcc1395_site_topology_pair_outcomes.tsv",
        "site_regions": args.site_dir / "hcc1395_site_topology_regions.tsv",
        "allele_identity": args.site_dir / "hcc1395_site_allele_identity.tsv.gz",
        "integrated_context": args.integrated_context,
    }

    outcome_rows, outcome_fields = load_tsv(
        paths["site_pair_outcomes"],
        {
            "match_id", "chrom", "region", "fixed_denominator", "caller_k_a", "caller_k_b",
            "caller_shared_k", "caller_site_set_relation", "shared_allele_identity_count",
            "shared_allele_mismatch_count", "allele_identity_pass", "read_full_evaluable",
            "read_full_candidate_relation", "read_full_outcome", "vaf_official_evaluable",
            "vaf_official_candidate_relation", "vaf_official_outcome",
        },
    )
    site_rows, site_fields = load_tsv(
        paths["site_regions"],
        {
            "match_id", "region", "swap_tolerant_read_full_status", "swap_tolerant_read_full_category",
            "swap_tolerant_read_full_compatible", "swap_tolerant_read_full_unique_exact",
            "swap_tolerant_vaf_official_status", "swap_tolerant_vaf_official_category",
            "swap_tolerant_vaf_official_compatible", "swap_tolerant_vaf_official_unique_exact",
        },
    )
    pair_rows_all, pair_fields = load_tsv(
        paths["pair_matches"],
        {
            "scenario", "match_id", "chrom", "region_a", "region_b", "complete_a", "complete_b",
            "category_a", "category_b", "category_agree", "tree_digest_signature_a", "tree_digest_signature_b",
        },
    )
    pair_rows = [
        row for row in pair_rows_all
        if row["scenario"] == "exact_coordinate" and boolean(row["complete_a"]) and boolean(row["complete_b"])
    ]
    vaf_rows, vaf_fields = load_tsv(
        paths["vaf_pair_regions"],
        {
            "match_id", "region", "uses_vaf_a", "uses_vaf_b", "unique_exact_a", "unique_exact_b",
            "candidate_ordered_agree", "candidate_swap_tolerant_agree", "both_shape_selectable",
            "shape_ordered_agree", "shape_swap_tolerant_agree", "both_actual_vaf_evaluable",
            "both_actual_vaf_single_shape",
        },
    )
    allele_rows, allele_fields = load_tsv(
        paths["allele_identity"],
        {"match_id", "chrom", "region", "position", "alleles_a", "alleles_b", "allele_identity"},
    )
    integrated = load_json(paths["integrated_context"])

    outcome_by_id = unique_index(outcome_rows, "match_id", "site outcome")
    site_by_id = unique_index(site_rows, "match_id", "site topology")
    pair_by_id = unique_index(pair_rows, "match_id", "coarse pair")
    vaf_by_id = unique_index(vaf_rows, "match_id", "VAF pair")
    allele_keys = set()
    allele_by_match: dict[str, dict[str, int]] = {}
    for row in allele_rows:
        key = (row["match_id"], row["position"])
        if key in allele_keys:
            raise ValueError(f"Duplicate allele identity key: {key}")
        allele_keys.add(key)
        bucket = allele_by_match.setdefault(row["match_id"], {"identity": 0, "mismatch": 0})
        bucket["identity" if boolean(row["allele_identity"]) else "mismatch"] += 1

    keyset = set(outcome_by_id)
    if len(outcome_rows) != FIXED_N or len(keyset) != FIXED_N:
        raise ValueError(f"Fixed spine must be {FIXED_N} unique rows, observed {len(outcome_rows)}/{len(keyset)}")
    for label, mapping in (("site topology", site_by_id), ("coarse pair", pair_by_id), ("VAF pair", vaf_by_id)):
        if set(mapping) != keyset:
            missing = sorted(keyset.difference(mapping))[:5]
            extra = sorted(set(mapping).difference(keyset))[:5]
            raise ValueError(f"{label} keyset mismatch; missing={missing}, extra={extra}")
    if set(allele_by_match) != keyset:
        raise ValueError("Allele identity evidence does not cover every fixed region")

    regional_rows: list[dict[str, Any]] = []
    for match_id in sorted(keyset):
        outcome = outcome_by_id[match_id]
        site = site_by_id[match_id]
        pair = pair_by_id[match_id]
        vaf = vaf_by_id[match_id]
        if len({outcome["region"], site["region"], pair["region_a"], pair["region_b"], vaf["region"]}) != 1:
            raise ValueError(f"Region mismatch for {match_id}")
        if integer(outcome["fixed_denominator"]) != FIXED_N:
            raise ValueError(f"Bad fixed denominator for {match_id}")
        k_a = integer(outcome["caller_k_a"])
        k_b = integer(outcome["caller_k_b"])
        k_shared = integer(outcome["caller_shared_k"])
        identity_count = integer(outcome["shared_allele_identity_count"])
        mismatch_count = integer(outcome["shared_allele_mismatch_count"])
        evidence_counts = allele_by_match[match_id]
        if (identity_count, mismatch_count) != (evidence_counts["identity"], evidence_counts["mismatch"]):
            raise ValueError(f"Allele evidence count mismatch for {match_id}")
        if identity_count + mismatch_count != k_shared:
            raise ValueError(f"Shared allele conservation failed for {match_id}")

        read_evaluable = boolean(outcome["read_full_evaluable"])
        vaf_evaluable = boolean(outcome["vaf_official_evaluable"])
        if read_evaluable != (site["swap_tolerant_read_full_status"] == "evaluable"):
            raise ValueError(f"Read evaluability mismatch for {match_id}")
        if vaf_evaluable != (site["swap_tolerant_vaf_official_status"] == "evaluable"):
            raise ValueError(f"VAF evaluability mismatch for {match_id}")
        if outcome["read_full_candidate_relation"] != site["swap_tolerant_read_full_category"]:
            raise ValueError(f"Read relation mismatch for {match_id}")
        if outcome["vaf_official_candidate_relation"] != site["swap_tolerant_vaf_official_category"]:
            raise ValueError(f"VAF relation mismatch for {match_id}")

        read_compatible = boolean(site["swap_tolerant_read_full_compatible"], allow_empty=True)
        read_unique = boolean(site["swap_tolerant_read_full_unique_exact"], allow_empty=True)
        vaf_compatible = boolean(site["swap_tolerant_vaf_official_compatible"], allow_empty=True)
        vaf_unique = boolean(site["swap_tolerant_vaf_official_unique_exact"], allow_empty=True)
        if read_evaluable != (read_compatible is not None and read_unique is not None):
            raise ValueError(f"Read compatibility nullability mismatch for {match_id}")
        if vaf_evaluable != (vaf_compatible is not None and vaf_unique is not None):
            raise ValueError(f"VAF compatibility nullability mismatch for {match_id}")

        tree_a = pair["tree_digest_signature_a"]
        tree_b = pair["tree_digest_signature_b"]
        if not tree_a or not tree_b:
            raise ValueError(f"Missing unranked tree digest for {match_id}")
        unique_vaf_exact_gate = all(
            boolean(vaf[field]) for field in ("uses_vaf_a", "uses_vaf_b", "unique_exact_a", "unique_exact_b")
        )
        regional_rows.append(
            {
                "match_id": match_id,
                "chrom": outcome["chrom"],
                "region": outcome["region"],
                "fixed_denominator": FIXED_N,
                "caller_k_a": k_a,
                "caller_k_b": k_b,
                "caller_shared_k": k_shared,
                "k_stratum": k_stratum(k_shared),
                "caller_site_set_relation": outcome["caller_site_set_relation"],
                "shared_allele_identity_count": identity_count,
                "shared_allele_mismatch_count": mismatch_count,
                "all_shared_alleles_identical": boolean(outcome["allele_identity_pass"]),
                "coarse_category_a": pair["category_a"],
                "coarse_category_b": pair["category_b"],
                "coarse_agree": boolean(pair["category_agree"]),
                "unranked_tree_digest_ordered_agree": tree_a == tree_b,
                "unranked_tree_digest_hp_swap_agree": unlabeled_signature(tree_a) == unlabeled_signature(tree_b),
                "read_evaluable": read_evaluable,
                "read_candidate_relation": outcome["read_full_candidate_relation"],
                "read_compatible": read_compatible,
                "read_unique_exact": read_unique,
                "read_outcome": outcome["read_full_outcome"],
                "vaf_official_evaluable": vaf_evaluable,
                "vaf_official_candidate_relation": outcome["vaf_official_candidate_relation"],
                "vaf_official_compatible": vaf_compatible,
                "vaf_official_unique_exact": vaf_unique,
                "vaf_official_outcome": outcome["vaf_official_outcome"],
                "structure_first_vaf_shape_evaluable": boolean(vaf["both_shape_selectable"]),
                "structure_first_vaf_shape_ordered_agree": boolean(vaf["shape_ordered_agree"], allow_empty=True),
                "structure_first_vaf_shape_hp_swap_agree": boolean(vaf["shape_swap_tolerant_agree"], allow_empty=True),
                "both_actual_vaf_evaluable": boolean(vaf["both_actual_vaf_evaluable"]),
                "both_actual_vaf_single_shape": boolean(vaf["both_actual_vaf_single_shape"]),
                "vaf_unique_mutation_labeled_exact_forest_evaluable": unique_vaf_exact_gate,
                "vaf_unique_exact_forest_ordered_agree": (
                    boolean(vaf["candidate_ordered_agree"], allow_empty=True) if unique_vaf_exact_gate else None
                ),
                "vaf_unique_exact_forest_hp_swap_agree": (
                    boolean(vaf["candidate_swap_tolerant_agree"], allow_empty=True) if unique_vaf_exact_gate else None
                ),
            }
        )

    by_stratum = {
        "ALL": regional_rows,
        **{name: [row for row in regional_rows if row["k_stratum"] == name] for name in K_STRATA if name != "ALL"},
    }
    if sum(len(by_stratum[name]) for name in K_STRATA if name != "ALL") != FIXED_N:
        raise ValueError("k-stratum subtotal does not conserve to 5,720")

    metric_specs = [
        MetricSpec("inventory", "region_count", "Regions in stratum", lambda row: True, lambda row: True,
                   "caller_shared_k", "Region inventory only"),
        MetricSpec("site_set", "site_set_equal", "Exact caller site set equal", lambda row: True,
                   lambda row: row["caller_site_set_relation"] == "equal", "caller_site_set_relation", "Site membership only"),
        MetricSpec("site_set", "site_set_a_subset_b", "HCC1395 proper subset of DORADO", lambda row: True,
                   lambda row: row["caller_site_set_relation"] == "a_proper_subset_b", "caller_site_set_relation", "Site membership only"),
        MetricSpec("site_set", "site_set_b_subset_a", "DORADO proper subset of HCC1395", lambda row: True,
                   lambda row: row["caller_site_set_relation"] == "b_proper_subset_a", "caller_site_set_relation", "Site membership only"),
        MetricSpec("site_set", "site_set_partial_overlap", "Partial-overlap caller site set", lambda row: True,
                   lambda row: row["caller_site_set_relation"] == "partial_overlap", "caller_site_set_relation", "Site membership only"),
        MetricSpec("site_set", "site_set_disjoint", "Disjoint caller site set", lambda row: True,
                   lambda row: row["caller_site_set_relation"] == "disjoint", "caller_site_set_relation", "Site membership only"),
        MetricSpec("allele", "all_shared_alleles_identical_region", "Regions with identity at every shared allele", lambda row: True,
                   lambda row: row["all_shared_alleles_identical"], "allele_identity_pass + allele evidence", "Allele identity; not topology"),
        MetricSpec("coarse", "coarse_five_state_agreement", "Coarse five-state agreement", lambda row: True,
                   lambda row: row["coarse_agree"], "category_a/category_b", "Agreement, not accuracy"),
        MetricSpec("unranked_exact_tree", "unranked_tree_digest_ordered_agreement", "Unranked mutation-labeled tree-set digest ordered agreement",
                   lambda row: True, lambda row: row["unranked_tree_digest_ordered_agree"], "tree_digest_signature_a/b", "Feasible candidate-set identity; not selected truth"),
        MetricSpec("unranked_exact_tree", "unranked_tree_digest_hp_swap_agreement", "Unranked mutation-labeled tree-set digest HP-swap agreement",
                   lambda row: True, lambda row: row["unranked_tree_digest_hp_swap_agree"], "tree_digest_signature_a/b", "Feasible candidate-set identity; not selected truth"),
        MetricSpec("read_constraint", "read_evaluable", "Read-only constraint evaluable", lambda row: True,
                   lambda row: row["read_evaluable"], "swap_tolerant_read_full_status", "Evaluability; same-read evidence"),
        MetricSpec("read_constraint", "read_compatible", "Read-only candidate spaces compatible", lambda row: row["read_evaluable"],
                   lambda row: row["read_compatible"], "swap_tolerant_read_full_compatible", "Candidate compatibility; ambiguity retained"),
        MetricSpec("read_constraint", "read_unique_exact", "Read-only unique exact projection", lambda row: row["read_evaluable"],
                   lambda row: row["read_unique_exact"], "swap_tolerant_read_full_unique_exact", "Conditional small subset"),
        MetricSpec("read_constraint", "read_strict_full_exact", "Read strict full exact", lambda row: row["read_evaluable"],
                   lambda row: row["read_outcome"] == "strict_full_exact", "read_full_outcome", "Full candidate relation equality"),
        MetricSpec("read_constraint", "read_a_induced_substructure", "Read HCC1395 induced substructure", lambda row: row["read_evaluable"],
                   lambda row: row["read_outcome"] == "A_induced_substructure", "read_full_outcome", "True induced substructure only"),
        MetricSpec("read_constraint", "read_b_induced_substructure", "Read DORADO induced substructure", lambda row: row["read_evaluable"],
                   lambda row: row["read_outcome"] == "B_induced_substructure", "read_full_outcome", "True induced substructure only"),
        MetricSpec("read_constraint", "read_strict_plus_induced", "Read strict exact + true induced substructure", lambda row: row["read_evaluable"],
                   lambda row: row["read_outcome"] in {"strict_full_exact", "A_induced_substructure", "B_induced_substructure"},
                   "read_full_outcome", "Structural agreement; excludes resolution-only/shared-core/overlap/conflict"),
        MetricSpec("vaf_constraint", "vaf_official_evaluable", "Official VAF constraint evaluable", lambda row: True,
                   lambda row: row["vaf_official_evaluable"], "swap_tolerant_vaf_official_status", "Same-read VAF heuristic"),
        MetricSpec("vaf_constraint", "vaf_official_compatible", "Official VAF candidate spaces compatible", lambda row: row["vaf_official_evaluable"],
                   lambda row: row["vaf_official_compatible"], "swap_tolerant_vaf_official_compatible", "Same-read VAF heuristic"),
        MetricSpec("vaf_constraint", "vaf_official_conflict", "Official VAF candidate-space conflict", lambda row: row["vaf_official_evaluable"],
                   lambda row: row["vaf_official_outcome"] == "conflict", "vaf_official_outcome", "Disjoint VAF-selected candidate spaces"),
        MetricSpec("vaf_constraint", "vaf_official_strict_full_exact", "Official VAF strict full exact", lambda row: row["vaf_official_evaluable"],
                   lambda row: row["vaf_official_outcome"] == "strict_full_exact", "vaf_official_outcome", "Same-read VAF heuristic"),
        MetricSpec("vaf_constraint", "vaf_official_a_induced_substructure", "Official VAF HCC1395 induced substructure", lambda row: row["vaf_official_evaluable"],
                   lambda row: row["vaf_official_outcome"] == "A_induced_substructure", "vaf_official_outcome", "True induced substructure only"),
        MetricSpec("vaf_constraint", "vaf_official_b_induced_substructure", "Official VAF DORADO induced substructure", lambda row: row["vaf_official_evaluable"],
                   lambda row: row["vaf_official_outcome"] == "B_induced_substructure", "vaf_official_outcome", "True induced substructure only"),
        MetricSpec("vaf_constraint", "vaf_official_strict_plus_induced", "Official VAF strict exact + true induced substructure", lambda row: row["vaf_official_evaluable"],
                   lambda row: row["vaf_official_outcome"] in {"strict_full_exact", "A_induced_substructure", "B_induced_substructure"},
                   "vaf_official_outcome", "Same-read VAF heuristic; excludes resolution-only/shared-core/overlap/conflict"),
        MetricSpec("structure_first_vaf_shape", "shape_evaluable", "Structure-first + VAF single-shape evaluable", lambda row: True,
                   lambda row: row["structure_first_vaf_shape_evaluable"], "both_shape_selectable", "Unlabeled shape; mutation direction removed"),
        MetricSpec("structure_first_vaf_shape", "shape_ordered_agreement", "Structure-first + VAF ordered-shape agreement",
                   lambda row: row["structure_first_vaf_shape_evaluable"], lambda row: row["structure_first_vaf_shape_ordered_agree"],
                   "shape_ordered_agree", "Unlabeled shape; mutation direction removed"),
        MetricSpec("structure_first_vaf_shape", "shape_hp_swap_agreement", "Structure-first + VAF HP-swap shape agreement",
                   lambda row: row["structure_first_vaf_shape_evaluable"], lambda row: row["structure_first_vaf_shape_hp_swap_agree"],
                   "shape_swap_tolerant_agree", "Unlabeled shape; mutation direction removed"),
        MetricSpec("vaf_unique_exact_forest", "unique_exact_forest_evaluable", "Both sides use VAF and each has one mutation-labeled exact forest",
                   lambda row: True, lambda row: row["vaf_unique_mutation_labeled_exact_forest_evaluable"],
                   "uses_vaf_a/b + unique_exact_a/b", "Selected same-read VAF subset"),
        MetricSpec("vaf_unique_exact_forest", "unique_exact_forest_ordered_agreement", "VAF-unique mutation-labeled exact forest ordered agreement",
                   lambda row: row["vaf_unique_mutation_labeled_exact_forest_evaluable"], lambda row: row["vaf_unique_exact_forest_ordered_agree"],
                   "candidate_ordered_agree", "Selected same-read VAF exact forest; not truth"),
        MetricSpec("vaf_unique_exact_forest", "unique_exact_forest_hp_swap_agreement", "VAF-unique mutation-labeled exact forest HP-swap agreement",
                   lambda row: row["vaf_unique_mutation_labeled_exact_forest_evaluable"], lambda row: row["vaf_unique_exact_forest_hp_swap_agree"],
                   "candidate_swap_tolerant_agree", "Selected same-read VAF exact forest; not truth"),
    ]

    metric_rows: list[dict[str, Any]] = []
    for stratum in K_STRATA:
        subset = by_stratum[stratum]
        for spec in metric_specs:
            evaluable = [row for row in subset if spec.gate(row)]
            numerator = sum(bool(spec.success(row)) for row in evaluable)
            availability = "AVAILABLE" if subset else "AVAILABLE_ZERO_STRATUM"
            metric_rows.append(
                {
                    "stratum": stratum,
                    "k_definition": "caller_shared_k across exact-coordinate region pair",
                    "metric_group": spec.group,
                    "metric_id": spec.metric_id,
                    "metric_label": spec.label,
                    "unit": "regions",
                    "global_fixed_denominator": FIXED_N,
                    "stratum_fixed_denominator": len(subset),
                    "evaluable_denominator": len(evaluable),
                    "numerator": numerator,
                    "rate_global_fixed": rate(numerator, FIXED_N),
                    "rate_stratum_fixed": rate(numerator, len(subset)),
                    "rate_evaluable": rate(numerator, len(evaluable)),
                    "availability": availability,
                    "source_fields": spec.source_fields,
                    "claim_ceiling": spec.claim_ceiling,
                }
            )

    k_rows = []
    for stratum in K_STRATA[1:]:
        subset = by_stratum[stratum]
        values = [row["caller_shared_k"] for row in subset]
        k_rows.append(
            {
                "k_stratum": stratum,
                "k_basis": "caller_shared_k",
                "region_count": len(subset),
                "rate_fixed_5720": len(subset) / FIXED_N,
                "observed_min_k": min(values) if values else None,
                "observed_max_k": max(values) if values else None,
                "availability": "AVAILABLE" if subset else "AVAILABLE_ZERO_STRATUM",
            }
        )

    relations = ("equal", "a_proper_subset_b", "b_proper_subset_a", "partial_overlap", "disjoint")
    site_set_rows = []
    allele_summary_rows = []
    for stratum in K_STRATA:
        subset = by_stratum[stratum]
        for relation in relations:
            count = sum(row["caller_site_set_relation"] == relation for row in subset)
            site_set_rows.append(
                {
                    "stratum": stratum,
                    "relation": relation,
                    "regions": count,
                    "stratum_denominator": len(subset),
                    "rate_within_stratum": rate(count, len(subset)),
                    "rate_fixed_5720": count / FIXED_N,
                }
            )
        identity = sum(row["shared_allele_identity_count"] for row in subset)
        mismatch = sum(row["shared_allele_mismatch_count"] for row in subset)
        allele_summary_rows.append(
            {
                "stratum": stratum,
                "regions": len(subset),
                "shared_alleles": identity + mismatch,
                "identity_alleles": identity,
                "mismatch_alleles": mismatch,
                "identity_rate": rate(identity, identity + mismatch),
                "all_identity_regions": sum(row["all_shared_alleles_identical"] for row in subset),
                "all_identity_region_rate": rate(sum(row["all_shared_alleles_identical"] for row in subset), len(subset)),
            }
        )

    confusion_rows = []
    for stratum in K_STRATA:
        subset = by_stratum[stratum]
        for category_a in COARSE_STATES:
            for category_b in COARSE_STATES:
                count = sum(
                    row["coarse_category_a"] == category_a and row["coarse_category_b"] == category_b
                    for row in subset
                )
                confusion_rows.append(
                    {
                        "stratum": stratum,
                        "category_a": category_a,
                        "category_b": category_b,
                        "regions": count,
                        "stratum_denominator": len(subset),
                        "rate_within_stratum": rate(count, len(subset)),
                        "diagonal": category_a == category_b,
                    }
                )

    availability_rows = [
        {"endpoint": "k=1/k=2/k=3/k>=4 region counts", "status": "AVAILABLE_ZERO_K1", "source": "caller_shared_k", "reason": "All 5,720 rows carry integer shared-k; k=1 is a supported zero, not missing data."},
        {"endpoint": "exact site-set equal/subset/partial", "status": "AVAILABLE", "source": "caller_site_set_relation", "reason": "Mutually exclusive region-level relation."},
        {"endpoint": "shared allele identity", "status": "AVAILABLE", "source": "allele identity TSV + pair outcome counts", "reason": "Unique match_id+position evidence reconciles to caller_shared_k."},
        {"endpoint": "coarse five-state agreement", "status": "AVAILABLE", "source": "exact complete-both pair rows", "reason": "Both category labels present for every fixed region."},
        {"endpoint": "structure-first + VAF shape", "status": "AVAILABLE", "source": "both_shape_selectable + ordered/swap agreement", "reason": "Conditional unlabeled-shape endpoint; direction removed."},
        {"endpoint": "mutation-labeled exact and HP-swap", "status": "AVAILABLE", "source": "unranked tree digests + VAF unique exact forests", "reason": "Candidate-set and selected-forest endpoints are reported separately."},
        {"endpoint": "induced substructure", "status": "AVAILABLE", "source": "read_full_outcome / vaf_official_outcome", "reason": "Only A/B_induced_substructure is counted; resolution-only/shared-core excluded."},
        {"endpoint": "read-only constraint compatibility", "status": "AVAILABLE", "source": "swap_tolerant_read_full_compatible", "reason": "Evaluability and compatible denominators are explicit."},
        {"endpoint": "unique biological true-tree accuracy", "status": "UNAVAILABLE_NO_ORTHOGONAL_TRUTH", "source": "none", "reason": "No single-cell, lineage, multi-region, or synthetic exact-tree truth in these artifacts."},
        {"endpoint": "independent VAF validation", "status": "UNAVAILABLE_SAME_READ_EVIDENCE", "source": "VAF ranking", "reason": "VAF-selected endpoints reuse the same reads and are heuristic, not independent truth."},
    ]

    checks: list[dict[str, Any]] = []

    def check(name: str, expected: Any, observed: Any, evidence: str) -> None:
        status = "PASS" if expected == observed else "FAIL"
        checks.append({"check": name, "expected": expected, "observed": observed, "status": status, "evidence": evidence})

    check("fixed spine rows", FIXED_N, len(regional_rows), "site topology pair outcomes")
    check("fixed spine unique match_id", FIXED_N, len({row["match_id"] for row in regional_rows}), "joined output")
    check("k subtotal", FIXED_N, sum(row["region_count"] for row in k_rows), "k_region_counts.tsv")
    check("site-set subtotal", FIXED_N, sum(row["regions"] for row in site_set_rows if row["stratum"] == "ALL"), "all five relations")
    check("coarse confusion subtotal", FIXED_N, sum(row["regions"] for row in confusion_rows if row["stratum"] == "ALL"), "5x5 confusion")
    check("shared allele conservation", sum(row["caller_shared_k"] for row in regional_rows),
          sum(row["shared_allele_identity_count"] + row["shared_allele_mismatch_count"] for row in regional_rows), "row-level allele counts")

    pair_context = integrated["hcc1395_pair"]
    check("integrated fixed denominator", pair_context["fixed_denominator"], len(regional_rows), "integrated context")
    check("integrated site equal", pair_context["site_inventory"]["caller_equal_site_regions"],
          sum(row["caller_site_set_relation"] == "equal" for row in regional_rows), "row-level site relation")
    check("integrated site A subset B", pair_context["site_inventory"]["caller_a_proper_subset_b_regions"],
          sum(row["caller_site_set_relation"] == "a_proper_subset_b" for row in regional_rows), "row-level site relation")
    check("integrated site B subset A", pair_context["site_inventory"]["caller_b_proper_subset_a_regions"],
          sum(row["caller_site_set_relation"] == "b_proper_subset_a" for row in regional_rows), "row-level site relation")
    check("integrated site partial overlap", pair_context["site_inventory"]["caller_partial_overlap_regions"],
          sum(row["caller_site_set_relation"] == "partial_overlap" for row in regional_rows), "row-level site relation")
    check("integrated shared allele identity", pair_context["site_inventory"]["shared_allele_identity"],
          sum(row["shared_allele_identity_count"] for row in regional_rows), "allele evidence")
    check("integrated shared allele mismatch", pair_context["site_inventory"]["shared_allele_mismatch"],
          sum(row["shared_allele_mismatch_count"] for row in regional_rows), "allele evidence")
    check("integrated coarse agreement", pair_context["coarse_topology"]["agree_n"],
          sum(row["coarse_agree"] for row in regional_rows), "coarse categories")
    check("integrated unranked exact tree ordered", pair_context["unranked_exact_candidate_tree_space"]["ordered_agree_n"],
          sum(row["unranked_tree_digest_ordered_agree"] for row in regional_rows), "tree digest signatures")
    check("integrated unranked exact tree HP-swap", pair_context["unranked_exact_candidate_tree_space"]["hp_swap_tolerant_agree_n"],
          sum(row["unranked_tree_digest_hp_swap_agree"] for row in regional_rows), "tree digest signatures")
    read_context = pair_context["read_full_candidate_constraints"]
    check("integrated read evaluable", read_context["evaluable_n"], sum(row["read_evaluable"] for row in regional_rows), "read status")
    check("integrated read strict exact", read_context["strict_full_exact_n"],
          sum(row["read_outcome"] == "strict_full_exact" for row in regional_rows), "read outcomes")
    check("integrated read induced", read_context["true_induced_substructure_n"],
          sum(row["read_outcome"] in {"A_induced_substructure", "B_induced_substructure"} for row in regional_rows), "read outcomes")
    check("integrated read strict+induced", read_context["strict_plus_true_induced_n"],
          sum(row["read_outcome"] in {"strict_full_exact", "A_induced_substructure", "B_induced_substructure"} for row in regional_rows), "read outcomes")
    check("read compatible", read_context["evaluable_n"] - read_context["conflict_n"],
          sum(row["read_compatible"] is True for row in regional_rows), "read compatibility")
    vaf_context = pair_context["vaf_official_candidate_constraints"]
    check("integrated VAF evaluable", vaf_context["evaluable_n"], sum(row["vaf_official_evaluable"] for row in regional_rows), "VAF status")
    check("integrated VAF strict exact", vaf_context["strict_full_exact_n"],
          sum(row["vaf_official_outcome"] == "strict_full_exact" for row in regional_rows), "VAF outcomes")
    check("integrated VAF induced", vaf_context["true_induced_substructure_n"],
          sum(row["vaf_official_outcome"] in {"A_induced_substructure", "B_induced_substructure"} for row in regional_rows), "VAF outcomes")
    check("integrated VAF strict+induced", vaf_context["strict_plus_true_induced_n"],
          sum(row["vaf_official_outcome"] in {"strict_full_exact", "A_induced_substructure", "B_induced_substructure"} for row in regional_rows), "VAF outcomes")
    check("integrated VAF conflicts", vaf_context["conflict_n"],
          sum(row["vaf_official_outcome"] == "conflict" for row in regional_rows), "VAF outcomes")
    selected_context = pair_context["vaf_selected_tree_and_shape"]
    unique_context = selected_context["unique_mutation_labeled_exact_forest"]
    check("integrated VAF unique exact forest evaluable", unique_context["n"],
          sum(row["vaf_unique_mutation_labeled_exact_forest_evaluable"] for row in regional_rows), "VAF exact forest gate")
    check("integrated VAF unique exact forest ordered", unique_context["ordered_agree_n"],
          sum(row["vaf_unique_exact_forest_ordered_agree"] is True for row in regional_rows), "VAF exact forest")
    check("integrated VAF unique exact forest HP-swap", unique_context["hp_swap_tolerant_agree_n"],
          sum(row["vaf_unique_exact_forest_hp_swap_agree"] is True for row in regional_rows), "VAF exact forest")
    shape_context = selected_context["structure_first_plus_vaf_single_shape"]
    check("integrated structure-first+VAF shape evaluable", shape_context["n"],
          sum(row["structure_first_vaf_shape_evaluable"] for row in regional_rows), "shape gate")
    check("integrated structure-first+VAF shape ordered", shape_context["ordered_agree_n"],
          sum(row["structure_first_vaf_shape_ordered_agree"] is True for row in regional_rows), "shape endpoint")
    check("integrated structure-first+VAF shape HP-swap", shape_context["hp_swap_tolerant_agree_n"],
          sum(row["structure_first_vaf_shape_hp_swap_agree"] is True for row in regional_rows), "shape endpoint")

    failed = [row for row in checks if row["status"] != "PASS"]
    if failed:
        raise RuntimeError("Regional precision validation failed: " + "; ".join(row["check"] for row in failed))

    input_inventory = []
    field_map = {
        "pair_matches": pair_fields,
        "vaf_pair_regions": vaf_fields,
        "site_pair_outcomes": outcome_fields,
        "site_regions": site_fields,
        "allele_identity": allele_fields,
        "integrated_context": sorted(integrated.keys()),
    }
    row_map = {
        "pair_matches": len(pair_rows_all),
        "vaf_pair_regions": len(vaf_rows),
        "site_pair_outcomes": len(outcome_rows),
        "site_regions": len(site_rows),
        "allele_identity": len(allele_rows),
        "integrated_context": 1,
    }
    key_map = {
        "pair_matches": "scenario+match_id; selected exact complete-both key=match_id",
        "vaf_pair_regions": "match_id",
        "site_pair_outcomes": "match_id",
        "site_regions": "match_id",
        "allele_identity": "match_id+position",
        "integrated_context": "singleton snapshot",
    }
    for source_id, path in paths.items():
        input_inventory.append(
            {
                "source_id": source_id,
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "rows": row_map[source_id],
                "columns": len(field_map[source_id]),
                "primary_key": key_map[source_id],
                "selected_fixed_rows": len(pair_rows) if source_id == "pair_matches" else (FIXED_N if source_id not in {"allele_identity", "integrated_context"} else ""),
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    regional_path = args.output_dir / "regional_precision_outcomes.tsv.gz"
    metrics_path = args.output_dir / "regional_precision_metrics.tsv"
    k_path = args.output_dir / "k_region_counts.tsv"
    site_path = args.output_dir / "site_set_relations_by_k.tsv"
    allele_path = args.output_dir / "allele_identity_by_k.tsv"
    confusion_path = args.output_dir / "coarse_confusion_by_k.tsv"
    availability_path = args.output_dir / "availability.tsv"
    inventory_path = args.output_dir / "input_schema_inventory.tsv"
    checks_path = args.output_dir / "validation_checks.tsv"
    summary_path = args.output_dir / "regional_precision_summary.json"
    record_path = args.output_dir / "execution_record.md"

    regional_fields = list(regional_rows[0].keys())
    metric_fields = list(metric_rows[0].keys())
    write_gzip_tsv(regional_path, regional_rows, regional_fields)
    write_tsv(metrics_path, metric_rows, metric_fields)
    write_tsv(k_path, k_rows, list(k_rows[0].keys()))
    write_tsv(site_path, site_set_rows, list(site_set_rows[0].keys()))
    write_tsv(allele_path, allele_summary_rows, list(allele_summary_rows[0].keys()))
    write_tsv(confusion_path, confusion_rows, list(confusion_rows[0].keys()))
    write_tsv(availability_path, availability_rows, list(availability_rows[0].keys()))
    write_tsv(inventory_path, input_inventory, list(input_inventory[0].keys()))
    write_tsv(checks_path, checks, list(checks[0].keys()))

    global_metric_map = {row["metric_id"]: row for row in metric_rows if row["stratum"] == "ALL"}
    summary = {
        "schema_name": "intersubmod.hcc1395_regional_precision_validation",
        "schema_version": "1.0.0",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "task_type": "B_comprehensive_validation",
        "goals": ["G4", "G5"],
        "status": "PASS_WITH_CLAIM_CEILING",
        "scope": {
            "genome": "GRCh38 chr1-22",
            "samples": [SAMPLE_A, SAMPLE_B],
            "fixed_population": FIXED_N,
            "fixed_population_definition": "exact-coordinate regions complete in both HCC1395 technical datasets",
            "k_definition": "caller_shared_k across the exact-coordinate pair",
            "k_counts": {row["k_stratum"]: row["region_count"] for row in k_rows},
        },
        "denominator_contract": {
            "global_fixed_denominator": FIXED_N,
            "stratum_fixed_denominator": "all regions in the selected k stratum, including not-evaluable endpoints",
            "evaluable_denominator": "regions satisfying the endpoint-specific gate",
            "allele_denominator": "shared CHROM:POS allele positions, separately reported from region denominators",
        },
        "headline": {
            "site_set_equal": global_metric_map["site_set_equal"],
            "shared_allele_identity": allele_summary_rows[0],
            "coarse_five_state_agreement": global_metric_map["coarse_five_state_agreement"],
            "read_compatible": global_metric_map["read_compatible"],
            "read_strict_plus_induced": global_metric_map["read_strict_plus_induced"],
            "vaf_official_strict_plus_induced": global_metric_map["vaf_official_strict_plus_induced"],
            "vaf_official_conflict": global_metric_map["vaf_official_conflict"],
            "shape_hp_swap_agreement": global_metric_map["shape_hp_swap_agreement"],
            "unique_exact_forest_hp_swap_agreement": global_metric_map["unique_exact_forest_hp_swap_agreement"],
        },
        "availability": availability_rows,
        "validation": {"checks": checks, "pass": len(checks), "fail": 0},
        "inputs": input_inventory,
        "outputs": {},
        "claim_ceiling": {
            "allowed": "fixed-population engineering coherence and cross-technical-source structural reproducibility at explicitly reported resolution",
            "not_allowed": "unique biological true-tree accuracy or independent validation by VAF/read evidence reused from the same observations",
        },
    }
    output_paths = [
        regional_path, metrics_path, k_path, site_path, allele_path, confusion_path,
        availability_path, inventory_path, checks_path,
    ]
    for output in output_paths:
        summary["outputs"][output.name] = {
            "path": str(output.resolve()), "bytes": output.stat().st_size, "sha256": sha256(output)
        }

    command = shlex.join([sys.executable, str(SCRIPT.resolve()), *sys.argv[1:]])
    record = f"""<!--
建立時間: {summary['generated_at']}
目標: 重建 HCC1395 兩技術來源固定 5,720-region 的位點、區域與拓撲精確分層
處理範圍: GRCh38 chr1-22; exact-coordinate complete-both regions; historical layered-v2 context
關聯檔案: InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/scripts/build_regional_precision_validation.py
-->

# Regional precision validation execution record

## 結論

**PASS WITH CLAIM CEILING** — {len(checks)}/{len(checks)} conservation checks PASS。本輪支持固定 population 的技術再現性，不支持唯一真實演化樹 accuracy。

## 輸入路徑

""" + "\n".join(f"- `{row['path']}` (sha256 `{row['sha256']}`)" for row in input_inventory) + f"""

## 執行命令

```bash
{command}
```

## 輸出路徑

""" + "\n".join(f"- `{path.resolve()}`" for path in [*output_paths, summary_path, record_path]) + f"""

## 實際輸出片段

```text
k=1: {len(by_stratum['k=1']):,}
k=2: {len(by_stratum['k=2']):,}
k=3: {len(by_stratum['k=3']):,}
k>=4: {len(by_stratum['k>=4']):,}
site-set equal: {global_metric_map['site_set_equal']['numerator']:,}/{FIXED_N:,}
shared allele identity: {allele_summary_rows[0]['identity_alleles']:,}/{allele_summary_rows[0]['shared_alleles']:,}
coarse 5-state: {global_metric_map['coarse_five_state_agreement']['numerator']:,}/{FIXED_N:,}
read compatible: {global_metric_map['read_compatible']['numerator']:,}/{global_metric_map['read_compatible']['evaluable_denominator']:,}
read strict+induced: {global_metric_map['read_strict_plus_induced']['numerator']:,}/{FIXED_N:,} fixed; /{global_metric_map['read_strict_plus_induced']['evaluable_denominator']:,} evaluable
VAF strict+induced: {global_metric_map['vaf_official_strict_plus_induced']['numerator']:,}/{FIXED_N:,} fixed; /{global_metric_map['vaf_official_strict_plus_induced']['evaluable_denominator']:,} evaluable
structure-first+VAF shape swap: {global_metric_map['shape_hp_swap_agreement']['numerator']:,}/{FIXED_N:,} fixed; /{global_metric_map['shape_hp_swap_agreement']['evaluable_denominator']:,} evaluable
VAF-unique exact forest swap: {global_metric_map['unique_exact_forest_hp_swap_agreement']['numerator']:,}/{FIXED_N:,} fixed; /{global_metric_map['unique_exact_forest_hp_swap_agreement']['evaluable_denominator']:,} evaluable
validation: {len(checks)}/{len(checks)} PASS
```

## 分母與可用性說明

- `global_fixed_denominator=5,720` 是固定 population；不因 endpoint 可否評估而改變。
- `stratum_fixed_denominator` 是 k 分層內全部 regions；`evaluable_denominator` 才是 endpoint-specific gate 後分母。
- k 以 `caller_shared_k` 定義。k=1 在這 5,720 regions 中是資料支持的 0，不是缺值。
- VAF/read endpoints 復用同批 reads；shape 移除 mutation labels/direction；二者都不是 orthogonal truth。
"""
    atomic_text(record_path, record)

    summary["outputs"][summary_path.name] = {
        "path": str(summary_path.resolve()),
        "receipt_note": "Self size/hash intentionally omitted to avoid a recursive self-receipt cycle",
    }
    summary["outputs"][record_path.name] = {
        "path": str(record_path.resolve()), "bytes": record_path.stat().st_size, "sha256": sha256(record_path)
    }
    write_json(summary_path, summary)

    print(json.dumps({
        "status": "PASS",
        "fixed_regions": FIXED_N,
        "k_counts": summary["scope"]["k_counts"],
        "validation_pass": len(checks),
        "validation_fail": 0,
        "output_dir": str(args.output_dir.resolve()),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
