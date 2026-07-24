#!/usr/bin/env python3
"""Build the exact-PS authority layered workstation.

The builder fails closed unless the 2026-07-24 all-seven exact-PS topology run,
the exact topology-signature census, and the k/HP funnel all pass their
hash/schema/arithmetic contracts.  It writes a cohort index plus seven
standalone sample pages to ``docs/methodology/_assets/layered_workstation``.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping


HERE = Path(__file__).resolve().parent
TOPIC_ROOT = HERE.parent
REPO_ROOT = HERE.parents[2]
OUTDIR = REPO_ROOT / "docs" / "methodology" / "_assets" / "layered_workstation"
TOPOLOGY_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1"
)
CENSUS_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260724_exact_ps_cpp_topology_signature_census/all7_v1"
)
FUNNEL_PATH = (
    TOPIC_ROOT / "data" / "20260724_exactPS_k_hp_funnel_census_01.json"
)
EXPECTED_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
UI_CONTRACT = "layered-workstation-exact-ps-v1"
AUTHORITY_NAME = "20260724-exact-ps-hp-strict-read-linkage"
CHROM_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
}
TARGET_REGION = ("chr10", 87818272, 87928739)


class AuthorityError(RuntimeError):
    """Raised when the frozen data authority cannot be proven."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuthorityError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuthorityError(f"cannot parse {path}: {exc}") from exc
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def iter_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    require(path.is_file(), f"missing JSONL: {path}")
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise AuthorityError(
                    f"invalid JSONL {path}:{line_number}: {exc}"
                ) from exc
            require(isinstance(row, dict), f"non-object row {path}:{line_number}")
            yield row


def all_true(value: Any, label: str) -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            all_true(child, f"{label}.{key}")
    elif isinstance(value, bool):
        require(value, f"false authority check: {label}")


def verify_identity(spec: Mapping[str, Any], label: str) -> Path:
    path = Path(str(spec.get("path", "")))
    require(path.is_file(), f"missing bound artifact {label}: {path}")
    expected_size = spec.get("size_bytes")
    require(
        isinstance(expected_size, int) and path.stat().st_size == expected_size,
        f"size drift for {label}: {path}",
    )
    expected_sha = spec.get("sha256")
    require(
        isinstance(expected_sha, str) and sha256_file(path) == expected_sha,
        f"SHA-256 drift for {label}: {path}",
    )
    return path


def verify_nested_identities(value: Any, label: str) -> int:
    count = 0
    if isinstance(value, Mapping):
        if {"path", "sha256", "size_bytes"} <= set(value):
            verify_identity(value, label)
            count += 1
        else:
            for key, child in value.items():
                count += verify_nested_identities(child, f"{label}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            count += verify_nested_identities(child, f"{label}[{index}]")
    return count


def pct(numerator: int | float, denominator: int | float) -> float:
    return 0.0 if not denominator else 100.0 * float(numerator) / float(denominator)


def json_for_html(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def atomic_write(path: Path, document: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = path.parent / f".{path.name}.stage.{os.getpid()}"
    with staging.open("w", encoding="utf-8") as handle:
        handle.write(document)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(staging, path)


def load_authority() -> dict[str, Any]:
    cohort_receipt_path = TOPOLOGY_ROOT / "cohort_receipt.json"
    cohort_receipt = load_json(cohort_receipt_path)
    require(
        cohort_receipt.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_cohort.cohort_receipt",
        "unexpected cohort receipt schema",
    )
    require(
        cohort_receipt.get("all_pass") is True
        and cohort_receipt.get("technical_all_pass") is True,
        "cohort technical gate did not pass",
    )
    scope = cohort_receipt.get("scope") or {}
    require(scope.get("autosomes_complete") is True, "cohort autosomes incomplete")
    require(scope.get("canonical_all7") is True, "cohort is not canonical all-seven")
    require(
        tuple(scope.get("samples") or ()) == EXPECTED_SAMPLES,
        "cohort sample order/scope mismatch",
    )
    require(
        tuple(scope.get("chromosomes") or ()) == tuple(CHROM_LENGTHS),
        "cohort chromosome scope mismatch",
    )

    all7_summary_path = TOPOLOGY_ROOT / "summary" / "all7_summary.json"
    all7_summary = load_json(all7_summary_path)
    require(
        all7_summary.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_af.cohort_summary",
        "unexpected all-seven topology summary schema",
    )
    require(all7_summary.get("all_pass") is True, "all-seven summary not PASS")
    all_true(all7_summary.get("checks") or {}, "all7_summary.checks")
    verify_identity(
        all7_summary.get("cohort_receipt_identity") or {},
        "all7_summary.cohort_receipt",
    )

    pipeline_identity_count = 0
    derived_funnel_samples: dict[str, dict[str, Any]] = {}
    for sample in EXPECTED_SAMPLES:
        receipt_specs = (cohort_receipt.get("sample_receipts") or {}).get(sample)
        require(isinstance(receipt_specs, Mapping), f"missing {sample} receipts")
        verify_identity(receipt_specs.get("partition") or {}, f"{sample}.partition")
        pipeline_path = verify_identity(
            receipt_specs.get("pipeline") or {}, f"{sample}.pipeline"
        )
        pipeline = load_json(pipeline_path)
        require(
            pipeline.get("all_pass") is True
            and pipeline.get("technical_all_pass") is True,
            f"{sample} pipeline did not pass",
        )
        require(pipeline.get("sample") == sample, f"{sample} pipeline sample mismatch")
        pipeline_identity_count += verify_nested_identities(
            pipeline.get("bound_artifacts") or [], f"{sample}.bound_artifacts"
        )

        mlhp_receipt_path = (
            TOPOLOGY_ROOT
            / "samples"
            / sample
            / f"{sample}.exact_ps_mlhp.json.receipt.json"
        )
        mlhp_receipt = load_json(mlhp_receipt_path)
        require(mlhp_receipt.get("all_pass") is True, f"{sample} MLHP receipt not PASS")
        funnel = mlhp_receipt.get("funnel") or {}
        require(
            funnel.get("check_cross_ps_zero") is True
            and funnel.get("check_cross_hp_zero") is True,
            f"{sample} MLHP contains cross-PS/HP leakage",
        )
        require(
            int((mlhp_receipt.get("counts") or {}).get("python_cpp_mismatch_count", -1))
            == 0,
            f"{sample} Python/C++ partition mismatch",
        )
        verify_identity(mlhp_receipt.get("output") or {}, f"{sample}.mlhp_output")
        derived_funnel_samples[sample] = {
            "funnel": {
                "source_W": int(funnel.get("tree_eligible_read_linked_regions", -1)),
                "bounded_blocks": int(funnel.get("bounded_blocks", -1)),
                "final_groups": int(funnel.get("tree_input_groups", -1)),
            },
            "source_overview": {
                "all_components": int(funnel.get("strict_components_all", -1))
            },
        }
        require(
            all(
                value >= 0
                for section in derived_funnel_samples[sample].values()
                for value in section.values()
            ),
            f"{sample} MLHP funnel fields are incomplete",
        )

    census_receipt_path = CENSUS_ROOT / "receipt.v2.json"
    census_receipt = load_json(census_receipt_path)
    require(
        census_receipt.get("schema_name")
        == "intersubmod.exact_ps_cpp_topology_signature_census.receipt",
        "unexpected census receipt schema",
    )
    all_true(census_receipt.get("checks") or {}, "census_receipt.checks")
    census_scope = census_receipt.get("scope") or {}
    require(
        tuple(census_scope.get("datasets") or ()) == EXPECTED_SAMPLES,
        "census sample scope mismatch",
    )
    require(
        int(census_scope.get("ranked_units", -1)) == 71955
        and int(census_scope.get("global_best_trees_enumerated", -1)) == 680527,
        "census scope totals mismatch",
    )
    census_identity_count = verify_nested_identities(
        {
            "sample_outputs": census_receipt.get("sample_outputs") or {},
            "summary": census_receipt.get("summary") or {},
            "implementation": census_receipt.get("implementation") or {},
        },
        "census_receipt",
    )
    require(census_identity_count >= 30, "census receipt identity coverage too small")

    census_summary = load_json(CENSUS_ROOT / "summary.json")
    cohort_census = census_summary.get("cohort") or {}
    all_true(cohort_census.get("checks") or {}, "census_summary.cohort.checks")
    require(
        int(cohort_census.get("ranked_units", -1)) == 71955,
        "census summary ranked denominator mismatch",
    )
    resolution = cohort_census.get("resolution") or {}
    require(
        sum(int((resolution.get(key) or {}).get("n", -1)) for key in (
            "UNIQUE_TREE",
            "TIED_SAME_TOPOLOGY",
            "TIED_CROSS_TOPOLOGY",
        ))
        == 71955,
        "census resolution conservation failed",
    )

    # The research data directory is intentionally git-ignored.  Derive the
    # renderer funnel from hash-bound MLHP receipts, and use the local census
    # JSON only as an additional cross-check when it is present.
    funnel = {"samples": derived_funnel_samples}
    funnel_sha256 = None
    if FUNNEL_PATH.is_file():
        local_funnel = load_json(FUNNEL_PATH)
        all_true(
            ((local_funnel.get("checks") or {}).get("cohort") or {}),
            "funnel.cohort",
        )
        per_sample_checks = (
            (local_funnel.get("checks") or {}).get("per_sample") or []
        )
        require(len(per_sample_checks) == 7, "funnel sample check count mismatch")
        for row in per_sample_checks:
            all_true(row, f"funnel.{row.get('dataset', 'unknown')}")
        local_samples = local_funnel.get("samples") or {}
        require(
            set(local_samples) == set(EXPECTED_SAMPLES),
            "local funnel sample set mismatch",
        )
        for sample in EXPECTED_SAMPLES:
            local = local_samples[sample]
            for key in ("source_W", "bounded_blocks", "final_groups"):
                require(
                    int((local.get("funnel") or {}).get(key, -1))
                    == int(derived_funnel_samples[sample]["funnel"][key]),
                    f"{sample} local/MLHP funnel mismatch: {key}",
                )
            require(
                int((local.get("source_overview") or {}).get("all_components", -1))
                == int(
                    derived_funnel_samples[sample]["source_overview"][
                        "all_components"
                    ]
                ),
                f"{sample} local/MLHP component total mismatch",
            )
        funnel_sha256 = sha256_file(FUNNEL_PATH)

    totals = all7_summary.get("totals") or {}
    require(int(totals.get("groups_total", -1)) == 98955, "group total mismatch")
    require(int(totals.get("ranked_units", -1)) == 71955, "ranked total mismatch")
    require(
        int(totals.get("mutation_family_abstain_units", -1)) == 10717,
        "ABSTAIN total mismatch",
    )

    return {
        "cohort_receipt": cohort_receipt,
        "cohort_receipt_path": cohort_receipt_path,
        "cohort_receipt_sha256": sha256_file(cohort_receipt_path),
        "all7_summary": all7_summary,
        "all7_summary_path": all7_summary_path,
        "all7_summary_sha256": sha256_file(all7_summary_path),
        "census_receipt": census_receipt,
        "census_receipt_path": census_receipt_path,
        "census_receipt_sha256": sha256_file(census_receipt_path),
        "census_summary": census_summary,
        "census_summary_path": CENSUS_ROOT / "summary.json",
        "census_summary_sha256": sha256_file(CENSUS_ROOT / "summary.json"),
        "funnel": funnel,
        "funnel_path": FUNNEL_PATH if FUNNEL_PATH.is_file() else None,
        "funnel_sha256": funnel_sha256,
        "pipeline_identity_count": pipeline_identity_count,
        "census_identity_count": census_identity_count,
    }


def sample_record_map(rows: Iterable[Mapping[str, Any]]) -> dict[str, Mapping[str, Any]]:
    result = {}
    for row in rows:
        sample = str(row.get("sample", ""))
        require(sample and sample not in result, f"duplicate/missing sample summary: {sample}")
        result[sample] = row
    require(set(result) == set(EXPECTED_SAMPLES), "sample summary set mismatch")
    return result


def resolution_key(topology: Mapping[str, Any], census: Mapping[str, Any] | None) -> str:
    if census:
        return str(census["resolution_class"])
    status = str(topology.get("read_af_status", "unknown"))
    return {
        "not_applicable_no_active_alt": "NO_ACTIVE_ALT",
        "zero_denominator": "ZERO_DENOMINATOR",
        "not_evaluable_family_incomplete": "RESOURCE_ABSTAIN",
        "recurrence_not_evaluable": "RECURRENCE_REQUIRED",
    }.get(status, "NOT_EVALUABLE")


def coarse_key(census: Mapping[str, Any] | None) -> str:
    if not census:
        return "Not evaluable"
    if int(census.get("coarse_class_count", 0)) != 1:
        return "Cross-class unresolved"
    counts = census.get("coarse_class_tree_counts") or {}
    present = [str(key) for key, count in counts.items() if int(count) > 0]
    return present[0] if len(present) == 1 else "Cross-class unresolved"


def compact_region(
    group: Mapping[str, Any],
    topology: Mapping[str, Any],
    census: Mapping[str, Any] | None,
) -> dict[str, Any]:
    positions = [int(value) for value in group.get("positions") or []]
    require(positions, f"empty region positions: {topology.get('region_id')}")
    hp = str(group.get("hp_family"))
    require(hp in {"1", "2"}, f"non-primary HP group: {topology.get('region_id')}")
    col_coverage = ((group.get("col_coverage_by_hp") or {}).get(hp) or {})
    coverage_fallback = [
        {
            "position": int(position),
            "ref_reads": int((col_coverage.get(str(position)) or [0, 0])[0]),
            "alt_reads": int((col_coverage.get(str(position)) or [0, 0])[1]),
        }
        for position in positions
    ]
    signatures = []
    if census:
        for signature in census.get("topology_signatures") or []:
            signatures.append(
                {
                    "shape": signature.get("shape_signature"),
                    "coarse": signature.get("coarse_class"),
                    "trees": str(signature.get("tree_count")),
                    "sha256": signature.get("shape_sha256"),
                }
            )
    populations = ((group.get("populations_by_hp") or {}).get(hp) or {})
    subreads = ((group.get("subread_groups_by_hp") or {}).get(hp) or {})
    return {
        "i": int(topology["group_index"]),
        "id": str(topology["region_id"]),
        "chrom": str(group["chrom"]),
        "start": int(group["start"]),
        "end": int(group["end"]),
        "mid": (int(group["start"]) + int(group["end"])) // 2,
        "span": int(group.get("span", int(group["end"]) - int(group["start"]))),
        "ps": str(group.get("phase_set")),
        "hp": hp,
        "unit": str(group.get("unit_id")),
        "block": str(group.get("block_id")),
        "component": str(group.get("component_id")),
        "n": int(group.get("n_sSNV", len(positions))),
        "activeK": int(topology.get("active_bit_count", 0)),
        "positions": positions,
        "resolution": resolution_key(topology, census),
        "coarse": coarse_key(census),
        "readStatus": str(topology.get("read_af_status", "")),
        "unitStatus": str(topology.get("unit_status", "")),
        "familyStatus": str(topology.get("family_status", "")),
        "objectiveStatus": str(topology.get("objective_status", "")),
        "reason": str(topology.get("reason", "")),
        "message": str(topology.get("message", "")),
        "score": str(topology.get("best_score_fraction", "")),
        "tieCount": str(topology.get("best_tree_tie_count", "0")),
        "vertexSetCount": int(topology.get("best_vertex_set_count", 0)),
        "treeUnique": bool(topology.get("best_tree_unique", False)),
        "recurrenceRequired": bool(topology.get("recurrence_required", False)),
        "searchNodes": int(topology.get("search_nodes", 0)),
        "coverage": topology.get("af_coverage") or coverage_fallback,
        "vertices": topology.get("representative_best_vertices") or [],
        "edges": topology.get("representative_best_edges") or [],
        "signatures": signatures,
        "exactShapeCount": int(census.get("topology_signature_count", 0))
        if census
        else 0,
        "coarseCount": int(census.get("coarse_class_count", 0)) if census else 0,
        "enumeratedTrees": str(census.get("enumerated_best_tree_count", "0"))
        if census
        else "0",
        "patterns": populations,
        "subreads": subreads,
        "fullReads": int(group.get("n_full_cov_reads", 0)),
        "hpReads": int((group.get("reads_by_hp") or {}).get(hp, 0)),
        "projectedRows": int(group.get("projected_molecule_rows", 0)),
        "treeMolecules": int(
            group.get("tree_supported_molecule_block_incidences", 0)
        ),
        "patternCount": int(group.get("tree_supported_pattern_count", 0)),
        "coverageMeaning": str(group.get("coverage_interpretation", "")),
    }


def load_sample(
    sample: str,
    authority: Mapping[str, Any],
    topology_summary: Mapping[str, Any],
    census_summary: Mapping[str, Any],
    funnel: Mapping[str, Any],
) -> dict[str, Any]:
    sample_root = TOPOLOGY_ROOT / "samples" / sample
    mlhp_path = sample_root / f"{sample}.exact_ps_mlhp.json"
    topology_path = sample_root / f"{sample}.topology.jsonl"
    census_path = CENSUS_ROOT / f"{sample}.census.jsonl"
    topology_receipt = load_json(sample_root / f"{sample}.topology.receipt.json")
    require(topology_receipt.get("all_pass") is True, f"{sample} topology not PASS")
    verify_identity(topology_receipt.get("input") or {}, f"{sample}.topology.input")
    verify_identity(topology_receipt.get("output") or {}, f"{sample}.topology.output")

    mlhp = load_json(mlhp_path)
    require(
        mlhp.get("schema_name") == "intersubmod.exact_ps_layered_tree_input",
        f"{sample} unexpected MLHP schema",
    )
    require(mlhp.get("sample") == sample, f"{sample} MLHP sample mismatch")
    groups = mlhp.get("groups") or []
    topology_rows = list(iter_jsonl(topology_path))
    census_rows = list(iter_jsonl(census_path))

    expected_groups = int(topology_summary.get("groups_total", -1))
    expected_ranked = int(census_summary.get("ranked_units", -1))
    require(
        len(groups) == len(topology_rows) == expected_groups,
        f"{sample} MLHP/topology/group total mismatch",
    )
    require(len(census_rows) == expected_ranked, f"{sample} census row total mismatch")

    census_by_index: dict[int, Mapping[str, Any]] = {}
    for row in census_rows:
        index = int(row.get("group_index", -1))
        require(index not in census_by_index, f"{sample} duplicate census group {index}")
        require(
            row.get("canonical_reproduction_pass") is True,
            f"{sample} census row failed reproduction: {index}",
        )
        census_by_index[index] = row

    rows = []
    ranked_indices = set()
    for index, (group, topology) in enumerate(zip(groups, topology_rows)):
        require(
            topology.get("schema_name")
            == "intersubmod.exact_ps_cpp_topology_af.unit"
            and topology.get("schema_version") == "1.0.0",
            f"{sample} unexpected topology schema at {index}",
        )
        require(
            int(topology.get("group_index", -1)) == index,
            f"{sample} topology order mismatch at {index}",
        )
        require(
            group.get("region_id") == topology.get("region_id"),
            f"{sample} MLHP/topology region mismatch at {index}",
        )
        require(
            group.get("unit_id") == topology.get("unit_id")
            and group.get("block_id") == topology.get("block_id")
            and group.get("chrom") == topology.get("chrom")
            and str(group.get("phase_set")) == str(topology.get("phase_set"))
            and str(group.get("hp_family")) == str(topology.get("hp_family")),
            f"{sample} MLHP/topology identity mismatch at {index}",
        )
        require(
            len(group.get("positions") or [])
            == int(group.get("n_sSNV", -1))
            == int(topology.get("original_bit_count", -2)),
            f"{sample} position/original-k mismatch at {index}",
        )
        census = census_by_index.get(index)
        if topology.get("read_af_status") == "ranked_complete":
            require(census is not None, f"{sample} ranked row missing census: {index}")
            ranked_indices.add(index)
        else:
            require(census is None, f"{sample} unranked row unexpectedly in census: {index}")
        if census:
            require(
                census.get("schema_name")
                == "intersubmod.exact_ps_cpp_topology_signature_census.unit"
                and census.get("schema_version") == "1.0.0",
                f"{sample} unexpected census schema at {index}",
            )
            require(
                census.get("region_id") == topology.get("region_id")
                and census.get("unit_id") == topology.get("unit_id")
                and census.get("block_id") == topology.get("block_id"),
                f"{sample} topology/census key mismatch at {index}",
            )
            require(
                int(census.get("active_bit_count", -1))
                == int(topology.get("active_bit_count", -2))
                and str(census.get("best_score_fraction"))
                == str(topology.get("best_score_fraction"))
                and str(census.get("best_tree_tie_count"))
                == str(topology.get("best_tree_tie_count")),
                f"{sample} topology/census AF result mismatch at {index}",
            )
            tie_count = int(census.get("best_tree_tie_count", -1))
            require(
                sum(
                    int(signature.get("tree_count", 0))
                    for signature in census.get("topology_signatures") or []
                )
                == tie_count
                == sum(
                    int(count)
                    for count in (census.get("coarse_class_tree_counts") or {}).values()
                ),
                f"{sample} signature/coarse tree conservation failed at {index}",
            )
        rows.append(compact_region(group, topology, census))

    require(
        ranked_indices == set(census_by_index),
        f"{sample} census is not exactly the ranked subset",
    )
    derived_resolution = Counter(row["resolution"] for row in rows)
    expected_resolution = census_summary.get("resolution") or {}
    for key in ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY"):
        require(
            derived_resolution[key] == int((expected_resolution.get(key) or {}).get("n", -1)),
            f"{sample} {key} derived count mismatch",
        )
    require(
        derived_resolution["RESOURCE_ABSTAIN"]
        == int(topology_summary.get("mutation_family_abstain_units", -1)),
        f"{sample} ABSTAIN count mismatch",
    )
    require(
        derived_resolution["RECURRENCE_REQUIRED"]
        == int(topology_summary.get("recurrence_required_units", -1)),
        f"{sample} recurrence-required count mismatch",
    )
    require(
        sum(1 for row in rows if row["chrom"] in CHROM_LENGTHS) == len(rows),
        f"{sample} contains non-autosomal final groups",
    )
    if sample in {"HCC1395", "HCC1395_DORADO"}:
        chrom, target_start, target_end = TARGET_REGION
        overlaps = [
            row
            for row in rows
            if row["chrom"] == chrom
            and row["start"] <= target_end
            and row["end"] >= target_start
        ]
        if sample == "HCC1395":
            require(len(overlaps) == 1, "HCC1395 target region must have one exact group")
            target = overlaps[0]
            require(
                target["positions"] == [87818272, 87840023]
                and target["resolution"] == "UNIQUE_TREE"
                and target["coarse"] == "Direct-only",
                "HCC1395 target exact-PS regression failed",
            )
        else:
            require(
                not overlaps,
                "DORADO target region must have no final exact-PS topology group",
            )

    compact_summary = {
        "groups": len(rows),
        "mutation": int(topology_summary["mutation_bearing_units"]),
        "ranked": len(census_rows),
        "abstain": int(topology_summary["mutation_family_abstain_units"]),
        "zero": int(topology_summary["zero_denominator_units"]),
        "noAlt": int(topology_summary["no_active_alt_units"]),
        "recurrence": int(topology_summary["recurrence_required_units"]),
        "uniqueTree": int((expected_resolution["UNIQUE_TREE"])["n"]),
        "sameTopology": int((expected_resolution["TIED_SAME_TOPOLOGY"])["n"]),
        "crossTopology": int((expected_resolution["TIED_CROSS_TOPOLOGY"])["n"]),
        "oneExact": int((census_summary.get("one_exact_topology") or {})["n"]),
        "oneCoarse": int((census_summary.get("one_coarse_class") or {})["n"]),
        "bestTrees": int(census_summary["best_tree_total"]),
    }
    return {
        "sample": sample,
        "rows": rows,
        "summary": compact_summary,
        "topology_summary": topology_summary,
        "census_summary": census_summary,
        "funnel": funnel,
        "paths": {
            "mlhp": str(mlhp_path),
            "topology": str(topology_path),
            "census": str(census_path),
            "pipelineReceipt": str(sample_root / "pipeline_receipt.json"),
            "censusReceipt": str(authority["census_receipt_path"]),
        },
        "hashes": {
            "mlhp": sha256_file(mlhp_path),
            "topology": sha256_file(topology_path),
            "census": sha256_file(census_path),
        },
    }


COMMON_CSS = r"""
:root{
  --paper:#f3f0e6;--panel:#fffdf7;--ink:#17231e;--muted:#667069;
  --line:#c9c8bd;--accent:#176b58;--accent2:#285f8f;--amber:#b66e20;
  --danger:#a94336;--purple:#6b5592;--soft:#e7e5da;--shadow:0 10px 32px rgba(30,42,35,.08);
}
*{box-sizing:border-box}
html{scroll-behavior:smooth;background:var(--paper)}
body{margin:0;color:var(--ink);background:var(--paper);font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI","Noto Sans TC",sans-serif;line-height:1.55}
a{color:var(--accent);text-underline-offset:3px}
button,input,select{font:inherit}
button{color:inherit}
.authority{position:sticky;top:0;z-index:30;display:flex;gap:.8rem;align-items:center;justify-content:center;padding:.55rem 1rem;background:#123c34;color:#fff;font-size:.82rem;letter-spacing:.02em}
.authority b{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.shell{width:min(1480px,calc(100% - 32px));margin:auto}
.hero{padding:3.5rem 0 2rem;border-bottom:1px solid var(--line)}
.eyebrow{font-size:.78rem;font-weight:800;letter-spacing:.12em;text-transform:uppercase;color:var(--accent)}
h1{font-family:Georgia,"Noto Serif TC",serif;font-size:clamp(2.2rem,5vw,5.4rem);line-height:.98;letter-spacing:-.045em;margin:.35rem 0 1rem;max-width:1100px}
h2{font-family:Georgia,"Noto Serif TC",serif;font-size:clamp(1.55rem,3vw,2.5rem);line-height:1.1;margin:0 0 .7rem}
h3{font-size:1rem;margin:0 0 .45rem}
p{margin:.35rem 0 .9rem}
.lead{font-size:clamp(1rem,2vw,1.3rem);max-width:900px;color:#3e4a44}
.boundary{display:grid;grid-template-columns:auto 1fr;gap:.7rem 1rem;margin-top:1.2rem;padding:1rem 1.1rem;border-left:5px solid var(--amber);background:#fff7e7;max-width:1100px}
.boundary strong{color:#7a4614}
.nav{display:flex;gap:.45rem;flex-wrap:wrap;padding:1rem 0}
.nav a,.back{display:inline-flex;align-items:center;min-height:40px;padding:.45rem .75rem;border:1px solid var(--line);background:var(--panel);text-decoration:none;color:var(--ink)}
main{padding-bottom:5rem}
.section{padding:2.4rem 0;border-bottom:1px solid var(--line)}
.section-head{display:flex;justify-content:space-between;gap:1rem;align-items:end;margin-bottom:1.15rem}
.section-head p{max-width:760px;color:var(--muted);margin:0}
.metrics{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:1px;background:var(--line);border:1px solid var(--line)}
.metric{background:var(--panel);padding:1.15rem;min-height:128px}
.metric .label{display:block;color:var(--muted);font-size:.78rem;min-height:2.4em}
.metric .value{display:block;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:clamp(1.45rem,2.5vw,2.45rem);font-weight:800;line-height:1.05;margin:.35rem 0}
.metric .note{display:block;font-size:.76rem;color:var(--muted)}
.grid-2{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem}
.grid-3{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:1rem}
.panel{background:var(--panel);border:1px solid var(--line);padding:1.1rem;box-shadow:var(--shadow)}
.panel.flat{box-shadow:none}
.denom{font-size:.76rem;color:var(--muted)}
.stack{display:flex;height:22px;overflow:hidden;border:1px solid rgba(0,0,0,.15);margin:.8rem 0 .6rem;background:#ddd}
.stack span{min-width:1px}
.dist{display:grid;gap:.45rem}
.dist-row{display:grid;grid-template-columns:minmax(120px,1fr) auto 70px;gap:.5rem;align-items:center;font-size:.82rem}
.swatch{width:.75rem;height:.75rem;display:inline-block;margin-right:.4rem;vertical-align:-.05rem}
.count{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;text-align:right}
.pct{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--muted);text-align:right}
.funnel{display:grid;grid-template-columns:repeat(4,1fr);gap:1.5rem;align-items:stretch}
.funnel-step{position:relative;padding:1.1rem;background:var(--panel);border:1px solid var(--line)}
.funnel-step:not(:last-child)::after{content:"→";position:absolute;right:-1.25rem;top:40%;font-size:1.5rem;color:var(--muted)}
.funnel-step b{display:block;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:1.45rem}
.table-wrap{overflow:auto;border:1px solid var(--line);background:var(--panel)}
table{border-collapse:collapse;width:100%;min-width:720px}
th,td{padding:.7rem .75rem;border-bottom:1px solid var(--line);text-align:left;vertical-align:top;font-size:.82rem}
th{position:sticky;top:0;background:#e9e7dd;font-size:.72rem;letter-spacing:.04em;text-transform:uppercase;z-index:2}
tr:last-child td{border-bottom:0}
td.num{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;text-align:right}
.sample-link{font-weight:800;font-size:1rem}
.barline{height:7px;background:var(--soft);min-width:100px;margin-top:.25rem}
.barline span{display:block;height:100%;background:var(--accent)}
.callout{padding:1.2rem;border:1px solid var(--line);background:#edf5f1}
.tag{display:inline-block;padding:.15rem .45rem;border:1px solid currentColor;font-size:.7rem;font-weight:800;letter-spacing:.04em;text-transform:uppercase;margin:.1rem .15rem .1rem 0}
.controls{display:grid;grid-template-columns:1fr auto;gap:1rem;align-items:end;margin:1rem 0}
.mode-buttons,.legend{display:flex;flex-wrap:wrap;gap:.35rem}
.mode-btn,.legend-btn,.action{min-height:40px;padding:.42rem .7rem;border:1px solid var(--line);background:var(--panel);cursor:pointer}
.mode-btn[aria-pressed=true]{background:var(--ink);color:#fff;border-color:var(--ink)}
.legend-btn{display:inline-flex;align-items:center;gap:.35rem}
.legend-btn[aria-pressed=true]{outline:3px solid color-mix(in srgb,var(--c),transparent 55%);border-color:var(--c);font-weight:800}
.legend-dot{width:.72rem;height:.72rem;background:var(--c);border-radius:50%}
.search{display:flex;gap:.4rem;min-width:min(100%,420px)}
.search input{width:100%;min-height:42px;padding:.5rem .65rem;border:1px solid var(--line);background:#fff}
.genome-wrap{border:1px solid var(--line);background:var(--panel);overflow:hidden}
#genomeCanvas{display:block;width:100%;height:660px;touch-action:manipulation;cursor:crosshair}
.genome-status{display:flex;justify-content:space-between;gap:1rem;padding:.55rem .75rem;border-top:1px solid var(--line);font-size:.78rem;color:var(--muted)}
.browser{display:grid;grid-template-columns:minmax(300px,.8fr) minmax(0,1.7fr);gap:1rem;margin-top:1rem;align-items:start}
.region-list{max-height:780px;overflow:auto;border:1px solid var(--line);background:var(--panel)}
.region-button{display:grid;width:100%;grid-template-columns:1fr auto;gap:.25rem .7rem;padding:.7rem;border:0;border-bottom:1px solid var(--line);background:transparent;text-align:left;cursor:pointer}
.region-button:hover,.region-button[aria-current=true]{background:#e8f1ed}
.region-button strong{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.78rem;overflow-wrap:anywhere}
.region-button small{color:var(--muted)}
.detail{min-height:460px;scroll-margin-top:58px}
.detail-head{display:flex;justify-content:space-between;gap:1rem;align-items:start}
.detail-id{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.92rem;overflow-wrap:anywhere}
.evidence-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1px;background:var(--line);border:1px solid var(--line);margin:.8rem 0}
.evidence-cell{background:var(--panel);padding:.7rem}
.evidence-cell span{display:block;color:var(--muted);font-size:.7rem}
.evidence-cell b{display:block;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere}
.detail-grid{display:grid;grid-template-columns:1.15fr .85fr;gap:1rem}
.tree{min-height:290px;border:1px solid var(--line);background:#fbfaf4;overflow:auto}
.tree svg{display:block;min-width:520px;width:100%;height:290px}
.mini-table{min-width:0}
.mini-table th,.mini-table td{font-size:.75rem;padding:.45rem}
.empty{padding:3rem 1rem;text-align:center;color:var(--muted)}
details{border:1px solid var(--line);background:var(--panel);margin:.7rem 0}
summary{cursor:pointer;padding:.85rem 1rem;font-weight:800}
.details-body{padding:0 1rem 1rem}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere}
.path{font-size:.75rem;display:block;margin:.4rem 0}
.footer{padding:2rem 0;color:var(--muted);font-size:.78rem}
.danger{color:var(--danger)}
.accent{color:var(--accent)}
.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
@media(max-width:980px){
  .metrics{grid-template-columns:repeat(2,minmax(0,1fr))}
  .grid-3,.funnel{grid-template-columns:repeat(2,minmax(0,1fr))}
  .funnel-step::after{display:none}
  .browser,.detail-grid{grid-template-columns:1fr}
  .region-list{max-height:360px}
  .evidence-grid{grid-template-columns:repeat(2,minmax(0,1fr))}
}
@media(max-width:640px){
  .shell{width:min(100% - 20px,1480px)}
  .authority{position:static;justify-content:flex-start;flex-wrap:wrap;white-space:normal;font-size:.7rem}
  .hero{padding:2rem 0 1.4rem}
  h1{font-size:2.45rem}
  .boundary{grid-template-columns:1fr}
  .metrics,.grid-2,.grid-3,.funnel{grid-template-columns:1fr}
  .metric{min-height:0}
  .section{padding:1.7rem 0}
  .section-head,.detail-head{display:block}
  .controls{grid-template-columns:1fr}
  .search{min-width:0}
  #genomeCanvas{height:610px}
  .genome-status{display:block}
  .evidence-grid{grid-template-columns:1fr}
}
"""


RESOLUTION_COLORS = {
    "UNIQUE_TREE": "#176b58",
    "TIED_SAME_TOPOLOGY": "#285f8f",
    "TIED_CROSS_TOPOLOGY": "#a94336",
    "NO_ACTIVE_ALT": "#b9b8ad",
    "ZERO_DENOMINATOR": "#736f68",
    "RESOURCE_ABSTAIN": "#b66e20",
    "RECURRENCE_REQUIRED": "#7d4f2f",
    "NOT_EVALUABLE": "#3f4441",
}
RESOLUTION_LABELS = {
    "UNIQUE_TREE": "唯一最佳樹",
    "TIED_SAME_TOPOLOGY": "多樹・同 topology",
    "TIED_CROSS_TOPOLOGY": "多樹・跨 topology",
    "NO_ACTIVE_ALT": "無 active ALT",
    "ZERO_DENOMINATOR": "read-AF 分母為 0",
    "RESOURCE_ABSTAIN": "資源 guard・ABSTAIN",
    "RECURRENCE_REQUIRED": "recurrence-required・不排序",
    "NOT_EVALUABLE": "不可評估",
}
COARSE_COLORS = {
    "Single-only": "#66717f",
    "Direct-only": "#176b58",
    "Sister-only": "#6b5592",
    "Sister+direct": "#b66e20",
    "Cross-class unresolved": "#a94336",
    "Not evaluable": "#b9b8ad",
}
COARSE_LABELS = {
    "Single-only": "單支",
    "Direct-only": "直系",
    "Sister-only": "旁系",
    "Sister+direct": "直系＋旁系",
    "Cross-class unresolved": "跨形態未定",
    "Not evaluable": "不可評估",
}


def distribution_html(
    counts: Mapping[str, int],
    denominator: int,
    colors: Mapping[str, str],
    labels: Mapping[str, str],
    order: Iterable[str],
) -> str:
    segments = []
    rows = []
    for key in order:
        count = int(counts.get(key, 0))
        if count:
            segments.append(
                f'<span title="{esc(labels.get(key, key))}: {count:,}" '
                f'style="width:{pct(count, denominator):.8f}%;background:{colors[key]}"></span>'
            )
        rows.append(
            '<div class="dist-row">'
            f'<span><i class="swatch" style="background:{colors[key]}"></i>{esc(labels.get(key, key))}</span>'
            f'<b class="count">{count:,}</b>'
            f'<span class="pct">{pct(count, denominator):.1f}%</span>'
            "</div>"
        )
    return (
        f'<div class="stack" aria-label="distribution">{"".join(segments)}</div>'
        f'<div class="dist">{"".join(rows)}</div>'
    )


def sample_page(authority: Mapping[str, Any], data: Mapping[str, Any]) -> str:
    sample = str(data["sample"])
    summary = data["summary"]
    rows = data["rows"]
    resolution_counts = Counter(row["resolution"] for row in rows)
    coarse_counts = Counter(row["coarse"] for row in rows)
    hp_counts = Counter(row["hp"] for row in rows)
    funnel = data["funnel"]["funnel"]
    exact_pct = pct(summary["oneExact"], summary["ranked"])
    ranked_pct = pct(summary["ranked"], summary["mutation"])

    page_payload = {
        "sample": sample,
        "summary": summary,
        "rows": rows,
        "chromLengths": CHROM_LENGTHS,
        "definitions": {
            "resolution": {
                key: {"label": RESOLUTION_LABELS[key], "color": color}
                for key, color in RESOLUTION_COLORS.items()
            },
            "coarse": {
                key: {"label": COARSE_LABELS[key], "color": color}
                for key, color in COARSE_COLORS.items()
            },
            "hp": {
                "1": {"label": "HP1", "color": "#2d6f91"},
                "2": {"label": "HP2", "color": "#a35d73"},
            },
            "status": {
                "ranked_complete": {"label": "ranked complete", "color": "#176b58"},
                "not_evaluable_family_incomplete": {
                    "label": "resource ABSTAIN",
                    "color": "#b66e20",
                },
                "zero_denominator": {"label": "zero denominator", "color": "#736f68"},
                "not_applicable_no_active_alt": {
                    "label": "no active ALT",
                    "color": "#b9b8ad",
                },
                "recurrence_not_evaluable": {
                    "label": "recurrence required",
                    "color": "#7d4f2f",
                },
            },
            "k": {
                str(key): {
                    "label": f"active k={key}" if key < 6 else "active k≥6",
                    "color": color,
                }
                for key, color in {
                    0: "#b9b8ad",
                    1: "#6483a3",
                    2: "#176b58",
                    3: "#4d8c69",
                    4: "#8a8b3f",
                    5: "#b66e20",
                    6: "#a94336",
                }.items()
            },
        },
    }
    is_hcc_pair = sample in {"HCC1395", "HCC1395_DORADO"}
    pair_note = (
        '<span class="tag">同一生物樣本・技術資料集</span>'
        if is_hcc_pair
        else '<span class="tag">獨立生物樣本</span>'
    )
    resolution_panel = distribution_html(
        resolution_counts,
        summary["ranked"],
        RESOLUTION_COLORS,
        RESOLUTION_LABELS,
        ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY"),
    )
    coarse_panel = distribution_html(
        coarse_counts,
        summary["ranked"],
        COARSE_COLORS,
        COARSE_LABELS,
        (
            "Single-only",
            "Direct-only",
            "Sister-only",
            "Sister+direct",
            "Cross-class unresolved",
        ),
    )
    status_panel = distribution_html(
        resolution_counts,
        summary["groups"],
        RESOLUTION_COLORS,
        RESOLUTION_LABELS,
        (
            "UNIQUE_TREE",
            "TIED_SAME_TOPOLOGY",
            "TIED_CROSS_TOPOLOGY",
            "ZERO_DENOMINATOR",
            "RESOURCE_ABSTAIN",
            "RECURRENCE_REQUIRED",
            "NO_ACTIVE_ALT",
        ),
    )
    provenance_rows = "".join(
        f'<span class="path"><b>{esc(label)}</b> <code>{esc(path)}</code></span>'
        for label, path in data["paths"].items()
    )
    hash_rows = "".join(
        f'<span class="path"><b>{esc(label)} SHA-256</b> <code>{esc(value)}</code></span>'
        for label, value in data["hashes"].items()
    )
    boundary_case = ""
    if sample == "HCC1395":
        boundary_case = """
<section class="section" id="boundary-case"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Boundary regression exemplar</div><h2>chr10:87818272–87928739 已不再是四點同一區</h2></div><p>這個例子直接展示「同 PS」仍不足以建立 topology unit；還必須有嚴格 read bridge。</p></div>
  <div class="grid-2">
    <article class="panel"><h3>新版 tree input</h3><p><b>{87818272, 87840023}</b> · exact PS=87549798 · HP1 · bridge support=8。結果為 <b>UNIQUE_TREE / Direct-only</b>，代表路徑 ROOT → H_AR → H_AA。</p><p><a href="#genome" onclick="setTimeout(()=>{const q=document.getElementById('regionSearch');q.value='chr10:87818272-87928739';document.getElementById('searchForm').requestSubmit()},50)">在全基因圖定位 →</a></p></article>
    <article class="panel"><h3>未進入同一 tree 的兩點</h3><p>S3=87888228、S4=87928739 與相鄰點 bridge support 都是 0，分別停在 source singleton abstain；頁面不再重現舊四點 RRAA→ARAA 路徑。</p><p class="denom">這證明 exact PS 是必要邊界，但同 PS 本身不是 topology edge。</p></article>
  </div>
</div></section>"""
    elif sample == "HCC1395_DORADO":
        boundary_case = """
<section class="section" id="boundary-case"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Technical pair boundary</div><h2>chr10:87818272–87928739 沒有可比較的 DORADO tree</h2></div><p>四個共同座標在 DORADO 都是 strict source singleton，沒有 final MLHP/topology/census row。</p></div>
  <div class="callout"><b>正確判讀：</b>可比較共同座標與邊際 read-AF，但 topology 應標示「不可比較」，不能寫成 topology=0，也不能用 HCC1395 的 HP1/PS 對位 DORADO 的 HP 標號。</div>
</div></section>"""
    document = f"""<!doctype html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<meta name="intersubmod-authority" content="{AUTHORITY_NAME}">
<meta name="intersubmod-ui-contract" content="{UI_CONTRACT}">
<meta name="intersubmod-canonical-sample" content="{esc(sample)}">
<meta name="intersubmod-cohort-receipt-sha256" content="{authority['cohort_receipt_sha256']}">
<meta name="intersubmod-census-receipt-sha256" content="{authority['census_receipt_sha256']}">
<title>Exact-PS 分層拓撲工作站 · {esc(sample)}</title>
<style>{COMMON_CSS}</style>
</head>
<body>
<div class="authority"><b>PRIMARY · 2026-07-24 exact PS × HP</b><span>strict endpoint read-linkage ≥3 · chr1–22 · receipt-bound</span></div>
<header class="hero">
  <div class="shell">
    <a class="back" href="index.html">← cohort 總覽</a>
    <div class="eyebrow">Regional mathematical topology · {esc(sample)}</div>
    <h1>精確 PS 分層<br>拓撲工作站</h1>
    <p class="lead">每個點都是同一 exact PS、同一 primary HP、且以嚴格 read-linkage 建立的 bounded block。read-AF 用於排序 recurrence-allowed 最小圖模型；不是 caller VAF、CCF 或 clone posterior。</p>
    <div>{pair_note}<span class="tag">7 datasets / 6 biological samples</span><span class="tag">GRCh38 chr1–22</span></div>
    <div class="boundary">
      <strong>判讀上限</strong>
      <span>「唯一最佳樹」只表示目前數學模型與 read-AF score 下只有一棵 global-best tree；不證明唯一細胞 clone、真實生物祖先、CN/LOH 正確或 posterior 已校準。代表樹只作一棵 exemplar；exact signature census 才是 tie 判定權威。</span>
    </div>
    <nav class="nav" aria-label="頁面導覽">
      <a href="#overview">樣本全貌</a><a href="#distributions">拓撲比例</a><a href="#genome">全基因組</a><a href="#method">方法與證據</a>
    </nav>
  </div>
</header>
<main>
<section class="section" id="overview"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Evidence ladder</div><h2>先看分母，再看拓撲</h2></div><p>同一頁同時保留 final groups、mutation-bearing、ranked complete 與 resource ABSTAIN；數字不跨分母混算。</p></div>
  <div class="metrics" data-testid="sample-metrics">
    <div class="metric"><span class="label">final topology groups</span><span class="value">{summary['groups']:,}</span><span class="note">exact PS × HP bounded blocks</span></div>
    <div class="metric"><span class="label">mutation-bearing</span><span class="value">{summary['mutation']:,}</span><span class="note">{pct(summary['mutation'],summary['groups']):.1f}% of final groups</span></div>
    <div class="metric"><span class="label">read-AF ranked complete</span><span class="value">{summary['ranked']:,}</span><span class="note">{ranked_pct:.1f}% of mutation-bearing</span></div>
    <div class="metric"><span class="label">單一 exact topology</span><span class="value">{exact_pct:.1f}%</span><span class="note">{summary['oneExact']:,} / {summary['ranked']:,}</span></div>
    <div class="metric"><span class="label">resource-limit ABSTAIN</span><span class="value danger">{summary['abstain']:,}</span><span class="note">不可當作 unresolved=negative</span></div>
  </div>
  <div class="funnel" style="margin-top:1rem">
    <div class="funnel-step"><span>source read-linked W</span><b>{int(funnel['source_W']):,}</b><small>strict endpoint components k≥2</small></div>
    <div class="funnel-step"><span>bounded blocks</span><b>{int(funnel['bounded_blocks']):,}</b><small>長區域依證據切分</small></div>
    <div class="funnel-step"><span>final groups</span><b>{int(funnel['final_groups']):,}</b><small>排除 k1 / pattern unsupported</small></div>
    <div class="funnel-step"><span>ranked complete</span><b>{summary['ranked']:,}</b><small>可做 exact signature census</small></div>
  </div>
</div></section>
<section class="section" id="distributions"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Sample-wide census</div><h2>三個問題、三個分母</h2></div><p>狀態圖以全部 final groups 為分母；determinacy 與形態只以 ranked complete 為分母。</p></div>
  <div class="grid-3">
    <article class="panel"><h3>能否完成 read-AF 判定？</h3><div class="denom">分母 = {summary['groups']:,} final groups</div>{status_panel}</article>
    <article class="panel"><h3>Determinacy｜目前辨識到哪一層</h3><div class="denom">分母 = {summary['ranked']:,} ranked complete</div>{resolution_panel}</article>
    <article class="panel"><h3>🌳 分子累積形狀</h3><div class="denom">分母 = {summary['ranked']:,} ranked complete</div>{coarse_panel}</article>
  </div>
</div></section>
{boundary_case}
<section class="section" id="genome"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">GRCh38 genome atlas</div><h2>chr1–22 全基因組分布</h2></div><p>圖例支援聯集多選：點一次選取，再點一次取消；沒有任何類別被選取時顯示全部。點圖上位置或下方清單檢視 exact PS region。</p></div>
  <div class="controls">
    <div>
      <span class="denom">著色維度</span>
      <div class="mode-buttons" id="modeButtons" data-testid="mode-buttons"></div>
    </div>
    <form class="search" id="searchForm">
      <label class="sr-only" for="regionSearch">搜尋座標</label>
      <input id="regionSearch" placeholder="chr10:87818272-87928739" autocomplete="off">
      <button class="action" type="submit">搜尋</button>
      <button class="action" type="button" id="clearSearch">清除</button>
    </form>
  </div>
  <div class="legend" id="legend" data-testid="legend"></div>
  <div class="genome-wrap">
    <canvas id="genomeCanvas" data-testid="genome-canvas" aria-label="GRCh38 chr1 到 chr22 拓撲區域分布"></canvas>
    <div class="genome-status"><span id="filterStatus"></span><span>座標按 GRCh38 chromosome bp 比例；點大小僅為可見性，不代表區域長度。</span></div>
  </div>
  <div class="browser">
    <aside>
      <div class="region-list" id="regionList" data-testid="region-list"></div>
    </aside>
    <article class="panel detail" id="regionDetail" data-testid="region-detail"><div class="empty">從全基因圖、座標搜尋或清單選擇一個 exact-PS region。</div></article>
  </div>
</div></section>
<section class="section" id="method"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Interpretation contract</div><h2>名詞與方法邊界</h2></div></div>
  <div class="grid-2">
    <article class="panel flat"><h3>🎯 Determinacy</h3><p><b>唯一最佳樹</b>：AF-global-best 完整樹只有一棵。<b>多樹・同 topology</b>：多個最佳 parent assignments 仍屬同一 rooted-unlabeled shape。<b>多樹・跨 topology</b>：最佳集合跨 shape，不能用 deterministic representative 冒充唯一形狀。</p></article>
    <article class="panel flat"><h3>🌳 拓撲型態</h3><p><b>單支</b>只有單一 mutation-state；<b>直系</b>有 root-to-node depth≥2；<b>旁系</b>有 outdegree≥2；<b>直系＋旁系</b>同時具有兩種圖形。這些是分子狀態圖幾何，不是 clone 身分。</p></article>
    <article class="panel flat"><h3>📍 基因體位置</h3><p>region 必須落在同一 exact PS 與同一 HP1/HP2，並由 threshold=3 嚴格 endpoint read-linkage component 支持；不再用相鄰 gap≤50 kb 傳遞合併。</p></article>
    <article class="panel flat"><h3>代表樹與完整候選</h3><p>detail network 顯示 frozen topology output 的一棵 deterministic global-best representative。完整最佳集合的權威是 signature 與 tree-count census；目前輸出未保存每棵 tied tree 的所有 edge list，因此頁面不捏造「完整 edge 聯集」。</p></article>
  </div>
  <details>
    <summary>機器證據與原始 JSON（預設收合）</summary>
    <div class="details-body">{provenance_rows}{hash_rows}
      <span class="path"><b>cohort receipt SHA-256</b> <code>{authority['cohort_receipt_sha256']}</code></span>
      <span class="path"><b>census receipt SHA-256</b> <code>{authority['census_receipt_sha256']}</code></span>
    </div>
  </details>
  <details>
    <summary>歷史 baseline（legacy 50 kb；不進入本頁統計）</summary>
    <div class="details-body"><p>2026-07-13 canonical-v5 頁面曾用 <code>adjacent_gap&lt;=50000</code> 傳遞合併與舊 92.18% topology 口徑。它保留為方法史，不是本頁 primary authority，也不能與新版 88.2579% exact rooted-unlabeled topology 比例直接混用。</p></div>
  </details>
</div></section>
</main>
<footer class="footer"><div class="shell">InterSubMod · {AUTHORITY_NAME} · UI contract {UI_CONTRACT}</div></footer>
<script type="application/json" id="pageData">{json_for_html(page_payload)}</script>
<script>{SAMPLE_JS}</script>
</body></html>"""
    return document


SAMPLE_JS = r"""
(() => {
  "use strict";
  const data = JSON.parse(document.getElementById("pageData").textContent);
  const rows = data.rows;
  const chroms = Object.keys(data.chromLengths);
  const state = {mode:"resolution", selected:new Set(), query:null, current:null, filtered:rows};
  const modes = [
    ["resolution","Determinacy"],["coarse","拓撲形態"],["hp","HP family"],
    ["status","判定狀態"],["k","active k"]
  ];
  const modeButtons = document.getElementById("modeButtons");
  const legend = document.getElementById("legend");
  const canvas = document.getElementById("genomeCanvas");
  const ctx = canvas.getContext("2d");
  const list = document.getElementById("regionList");
  const detail = document.getElementById("regionDetail");
  const status = document.getElementById("filterStatus");
  let hitPoints = [];

  function modeValue(row) {
    if (state.mode === "status") return row.readStatus;
    if (state.mode === "k") return String(Math.min(6,row.activeK));
    return String(row[state.mode]);
  }
  function defs(){ return data.definitions[state.mode]; }
  function parseQuery(raw) {
    const match = raw.trim().replace(/,/g,"").match(/^(chr(?:[1-9]|1\d|2[0-2]))(?::(\d+)(?:-(\d+))?)?$/i);
    if(!match) return null;
    const chrom = "chr" + match[1].slice(3);
    const start = match[2] ? Number(match[2]) : 1;
    const end = match[3] ? Number(match[3]) : start;
    return {chrom,start:Math.min(start,end),end:Math.max(start,end)};
  }
  function matchesQuery(row){
    if(!state.query) return true;
    return row.chrom===state.query.chrom && row.start<=state.query.end && row.end>=state.query.start;
  }
  function refresh() {
    state.filtered = rows.filter(row => matchesQuery(row) && (!state.selected.size || state.selected.has(modeValue(row))));
    renderCanvas(); renderList();
    status.textContent = `顯示 ${state.filtered.length.toLocaleString()} / ${rows.length.toLocaleString()} regions` +
      (state.selected.size ? ` · ${state.selected.size} 類聯集` : " · 全部類別") +
      (state.query ? ` · ${state.query.chrom}:${state.query.start.toLocaleString()}–${state.query.end.toLocaleString()}` : "");
  }
  function renderModes(){
    modeButtons.innerHTML="";
    for(const [key,label] of modes){
      const button=document.createElement("button");
      button.type="button";button.className="mode-btn";button.textContent=label;
      button.setAttribute("aria-pressed",String(state.mode===key));
      button.addEventListener("click",()=>{if(state.mode===key)return;state.mode=key;state.selected.clear();renderModes();renderLegend();refresh();});
      modeButtons.append(button);
    }
  }
  function renderLegend(){
    legend.innerHTML="";
    const counts=new Map();
    rows.filter(matchesQuery).forEach(row=>counts.set(modeValue(row),(counts.get(modeValue(row))||0)+1));
    for(const [key,spec] of Object.entries(defs())){
      if(!counts.get(key)) continue;
      const button=document.createElement("button");button.type="button";button.className="legend-btn";
      button.style.setProperty("--c",spec.color);button.setAttribute("aria-pressed",String(state.selected.has(key)));
      button.innerHTML=`<i class="legend-dot"></i><span>${spec.label}</span><b class="count">${counts.get(key).toLocaleString()}</b>`;
      button.addEventListener("click",()=>{state.selected.has(key)?state.selected.delete(key):state.selected.add(key);renderLegend();refresh();});
      legend.append(button);
    }
  }
  function resizeCanvas(){
    const rect=canvas.getBoundingClientRect();const dpr=Math.min(window.devicePixelRatio||1,2);
    canvas.width=Math.max(320,Math.floor(rect.width*dpr));canvas.height=Math.floor(rect.height*dpr);
    ctx.setTransform(dpr,0,0,dpr,0,0);
    return {w:rect.width,h:rect.height};
  }
  function renderCanvas(){
    const {w,h}=resizeCanvas();ctx.clearRect(0,0,w,h);hitPoints=[];
    const left=w<560?48:72,right=16,top=18,rowH=(h-top-12)/chroms.length;
    const rowSet=new Set(state.filtered);
    ctx.font=`${w<560?10:11}px ui-monospace,monospace`;ctx.textBaseline="middle";
    chroms.forEach((chrom,index)=>{
      const y=top+index*rowH+rowH/2;
      ctx.fillStyle="#425048";ctx.textAlign="right";ctx.fillText(chrom,left-8,y);
      ctx.strokeStyle="#cfcec4";ctx.lineWidth=Math.max(5,rowH*.34);ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(w-right,y);ctx.stroke();
      ctx.strokeStyle="#8c918c";ctx.lineWidth=.6;ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(w-right,y);ctx.stroke();
    });
    ctx.globalAlpha=.72;
    for(const row of state.filtered){
      const index=chroms.indexOf(row.chrom);if(index<0)continue;
      const y=top+index*rowH+rowH/2;
      const x=left+(w-left-right)*(row.mid/data.chromLengths[row.chrom]);
      const spec=defs()[modeValue(row)]||{color:"#444"};
      ctx.fillStyle=spec.color;ctx.beginPath();ctx.arc(x,y,w<560?1.7:2.2,0,Math.PI*2);ctx.fill();
      hitPoints.push({x,y,row});
    }
    ctx.globalAlpha=1;
    if(state.current && rowSet.has(state.current)){
      const index=chroms.indexOf(state.current.chrom);const y=top+index*rowH+rowH/2;
      const x=left+(w-left-right)*(state.current.mid/data.chromLengths[state.current.chrom]);
      ctx.strokeStyle="#111";ctx.lineWidth=2;ctx.beginPath();ctx.arc(x,y,6,0,Math.PI*2);ctx.stroke();
    }
  }
  function renderList(){
    list.innerHTML="";
    const shown=state.filtered.slice(0,160);
    if(!shown.length){list.innerHTML='<div class="empty">沒有符合目前聯集與座標條件的 region。</div>';return;}
    for(const row of shown){
      const button=document.createElement("button");button.type="button";button.className="region-button";
      button.setAttribute("aria-current",String(state.current===row));
      const spec=defs()[modeValue(row)]||{label:modeValue(row),color:"#555"};
      button.innerHTML=`<strong>${row.chrom}:${row.start.toLocaleString()}–${row.end.toLocaleString()}</strong><span class="tag" style="color:${spec.color}">${spec.label}</span><small>PS ${row.ps} · HP${row.hp} · n=${row.n} · active k=${row.activeK}</small><small>${row.span.toLocaleString()} bp</small>`;
      button.addEventListener("click",()=>selectRow(row,true));list.append(button);
    }
    if(state.filtered.length>shown.length){
      const note=document.createElement("div");note.className="denom";note.style.padding=".7rem";
      note.textContent=`清單只顯示前 ${shown.length} 筆；全基因圖仍包含全部 ${state.filtered.length.toLocaleString()} 筆。`;
      list.append(note);
    }
  }
  function escapeHtml(value){
    return String(value??"").replace(/[&<>"']/g,ch=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[ch]));
  }
  function fractionPct(frac){
    if(!frac||!String(frac).includes("/"))return "N/A";
    const [a,b]=String(frac).split("/").map(Number);return b?`${(100*a/b).toFixed(1)}%`:"N/A";
  }
  function treeSvg(row){
    if(!row.vertices.length||!row.edges.length)return '<div class="empty">此 group 沒有可展示的 ranked representative tree。</div>';
    const vertices=new Map(row.vertices.map(v=>[String(v.vertex),v]));
    const parent=new Map(row.edges.map(e=>[String(e.child_vertex),String(e.parent_vertex)]));
    const depth=id=>{let d=0,cur=id,seen=new Set();while(parent.has(cur)&&!seen.has(cur)){seen.add(cur);cur=parent.get(cur);d++;}return d;};
    const levels=new Map();
    for(const [id,v] of vertices){const d=depth(id);if(!levels.has(d))levels.set(d,[]);levels.get(d).push([id,v]);}
    const maxDepth=Math.max(...levels.keys());const pos=new Map();
    for(const [d,items] of levels){items.sort((a,b)=>Number(a[0])-Number(b[0]));items.forEach(([id],i)=>pos.set(id,{x:70+d*(420/Math.max(1,maxDepth)),y:items.length===1?145:48+i*(190/Math.max(1,items.length-1))}));}
    const edges=row.edges.map(e=>{const a=pos.get(String(e.parent_vertex)),b=pos.get(String(e.child_vertex));return `<g><line x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" stroke="#52625a" stroke-width="2" marker-end="url(#arrow)"/><text x="${(a.x+b.x)/2}" y="${(a.y+b.y)/2-6}" text-anchor="middle" font-size="10" fill="#667069">${escapeHtml(e.acquired_position??"")}</text></g>`}).join("");
    const nodes=[...vertices].map(([id,v])=>{const p=pos.get(id);return `<g><circle cx="${p.x}" cy="${p.y}" r="19" fill="${id==="0"?"#17231e":"#fffdf7"}" stroke="#176b58" stroke-width="2"/><text x="${p.x}" y="${p.y+4}" text-anchor="middle" font-size="10" fill="${id==="0"?"#fff":"#17231e"}">${escapeHtml(v.label)}</text></g>`}).join("");
    return `<svg viewBox="0 0 520 290" role="img" aria-label="一棵 deterministic global-best representative tree"><defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#52625a"/></marker></defs>${edges}${nodes}</svg>`;
  }
  function patternRows(obj,label){
    const entries=Object.entries(obj||{}).sort((a,b)=>b[1]-a[1]);
    if(!entries.length)return `<tr><td>${label}</td><td colspan="2">無</td></tr>`;
    return entries.slice(0,24).map(([p,n])=>`<tr><td>${label}</td><td class="mono">${escapeHtml(p)}</td><td class="num">${Number(n).toLocaleString()}</td></tr>`).join("");
  }
  function renderDetail(row){
    const rSpec=data.definitions.resolution[row.resolution]||{label:row.resolution,color:"#444"};
    const cSpec=data.definitions.coarse[row.coarse]||{label:row.coarse,color:"#444"};
    const coverage=(row.coverage||[]).map(c=>`<tr><td class="mono">${Number(c.position).toLocaleString()}</td><td class="num">${Number(c.ref_reads||0).toLocaleString()}</td><td class="num">${Number(c.alt_reads||0).toLocaleString()}</td><td class="num">${c.fraction?escapeHtml(c.fraction):"—"}</td><td class="num">${c.fraction?fractionPct(c.fraction):((Number(c.ref_reads||0)+Number(c.alt_reads||0))?`${(100*Number(c.alt_reads||0)/(Number(c.ref_reads||0)+Number(c.alt_reads||0))).toFixed(1)}%`:"N/A")}</td></tr>`).join("");
    const signatures=(row.signatures||[]).map(s=>`<tr><td class="mono">${escapeHtml(s.shape)}</td><td>${escapeHtml(s.coarse)}</td><td class="num">${escapeHtml(s.trees)}</td></tr>`).join("")||'<tr><td colspan="3">只對 ranked complete groups 提供 exact signature census。</td></tr>';
    detail.innerHTML=`
      <div class="detail-head"><div><div class="eyebrow">Exact-PS regional evidence</div><h3 class="detail-id">${escapeHtml(row.id)}</h3></div><div><span class="tag" style="color:${rSpec.color}">${rSpec.label}</span><span class="tag" style="color:${cSpec.color}">${cSpec.label}</span></div></div>
      <div class="evidence-grid">
        <div class="evidence-cell"><span>GRCh38</span><b>${row.chrom}:${row.start.toLocaleString()}–${row.end.toLocaleString()}</b></div>
        <div class="evidence-cell"><span>exact phase set</span><b>PS ${escapeHtml(row.ps)}</b></div>
        <div class="evidence-cell"><span>primary family</span><b>HP${row.hp}</b></div>
        <div class="evidence-cell"><span>sites</span><b>n=${row.n} · active k=${row.activeK}</b></div>
        <div class="evidence-cell"><span>best tree ties</span><b>${escapeHtml(row.tieCount)}</b></div>
        <div class="evidence-cell"><span>exact shape count</span><b>${row.exactShapeCount||"—"}</b></div>
        <div class="evidence-cell"><span>AF ordering score</span><b>${escapeHtml(row.score||"N/A")}</b></div>
        <div class="evidence-cell"><span>status / reason</span><b>${escapeHtml(row.readStatus)} · ${escapeHtml(row.reason)}</b></div>
      </div>
      <div class="callout"><b>邊界證據：</b>同一 PS=${escapeHtml(row.ps)}、同一 HP${row.hp}、threshold=3 strict endpoint component；unit <code>${escapeHtml(row.unit)}</code> / block <code>${escapeHtml(row.block)}</code>。此 group 不由 50 kb gap 規則建立。</div>
      <div class="detail-grid" style="margin-top:1rem">
        <div><h3>代表最佳樹</h3><div class="tree">${treeSvg(row)}</div><p class="denom">僅一棵 deterministic representative；若 tieCount&gt;1，完整同形／跨形判定請看右側 signature census。</p></div>
        <div><h3>Exact topology signatures</h3><div class="table-wrap"><table class="mini-table"><thead><tr><th>rooted shape</th><th>coarse</th><th>best trees</th></tr></thead><tbody>${signatures}</tbody></table></div></div>
      </div>
      <details open><summary>每個 sSNV 的 HP-specific read-ALT 分層</summary><div class="details-body"><div class="table-wrap"><table class="mini-table"><thead><tr><th>position</th><th>REF</th><th>ALT</th><th>exact fraction</th><th>ALT/(R+A)</th></tr></thead><tbody>${coverage}</tbody></table></div><p class="denom">${escapeHtml(row.coverageMeaning)}</p></div></details>
      <details><summary>分子 pattern 與 partial-read 證據</summary><div class="details-body"><p>full coverage=${row.fullReads.toLocaleString()} · HP reads=${row.hpReads.toLocaleString()} · projected rows=${row.projectedRows.toLocaleString()} · tree-supported molecule incidences=${row.treeMolecules.toLocaleString()}</p><div class="table-wrap"><table class="mini-table"><thead><tr><th>來源</th><th>pattern</th><th>reads</th></tr></thead><tbody>${patternRows(row.patterns,"full")}${patternRows(row.subreads,"partial")}</tbody></table></div></div></details>
      <details><summary>識別鍵與 solver 狀態</summary><div class="details-body"><span class="path"><b>component</b> <code>${escapeHtml(row.component)}</code></span><span class="path"><b>unit</b> <code>${escapeHtml(row.unit)}</code></span><span class="path"><b>block</b> <code>${escapeHtml(row.block)}</code></span><span class="path"><b>family/objective</b> <code>${escapeHtml(row.familyStatus)} / ${escapeHtml(row.objectiveStatus)}</code></span><span class="path"><b>search nodes</b> <code>${row.searchNodes.toLocaleString()}</code></span><p>${escapeHtml(row.message)}</p></div></details>`;
  }
  function selectRow(row,scroll){
    state.current=row;renderDetail(row);renderCanvas();renderList();
    if(scroll&&window.innerWidth<980) detail.scrollIntoView({behavior:"smooth",block:"start"});
  }
  canvas.addEventListener("click",event=>{
    const rect=canvas.getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top;
    let best=null,dist=Infinity;for(const hit of hitPoints){const d=Math.hypot(hit.x-x,hit.y-y);if(d<dist){best=hit;dist=d;}}
    if(best&&dist<12)selectRow(best.row,false);
  });
  document.getElementById("searchForm").addEventListener("submit",event=>{
    event.preventDefault();const raw=document.getElementById("regionSearch").value;
    const parsed=parseQuery(raw);if(!parsed){document.getElementById("regionSearch").setCustomValidity("格式例：chr10:87818272-87928739");document.getElementById("regionSearch").reportValidity();return;}
    document.getElementById("regionSearch").setCustomValidity("");state.query=parsed;renderLegend();refresh();if(state.filtered.length)selectRow(state.filtered[0],false);
  });
  document.getElementById("clearSearch").addEventListener("click",()=>{state.query=null;document.getElementById("regionSearch").value="";renderLegend();refresh();});
  let resizeTimer;window.addEventListener("resize",()=>{clearTimeout(resizeTimer);resizeTimer=setTimeout(renderCanvas,120);});
  renderModes();renderLegend();refresh();
})();
"""


def js_divergence_similarity(a: list[float], b: list[float]) -> float:
    total_a, total_b = sum(a), sum(b)
    if total_a <= 0 or total_b <= 0:
        return 0.0
    p = [value / total_a for value in a]
    q = [value / total_b for value in b]
    m = [(x + y) / 2 for x, y in zip(p, q)]

    def kl(x: list[float], y: list[float]) -> float:
        return sum(v * math.log(v / w, 2) for v, w in zip(x, y) if v > 0 and w > 0)

    distance = math.sqrt(max(0.0, (kl(p, m) + kl(q, m)) / 2))
    return max(0.0, 1.0 - distance)


def index_page(authority: Mapping[str, Any], samples: list[Mapping[str, Any]]) -> str:
    totals = authority["all7_summary"]["totals"]
    census = authority["census_summary"]["cohort"]
    funnel_samples = authority["funnel"]["samples"]
    source_components = sum(
        int(row["source_overview"]["all_components"]) for row in funnel_samples.values()
    )
    source_w = sum(int(row["funnel"]["source_W"]) for row in funnel_samples.values())
    bounded = sum(int(row["funnel"]["bounded_blocks"]) for row in funnel_samples.values())
    resolution = {
        key: int((census["resolution"][key])["n"])
        for key in ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY")
    }
    coarse = {
        key: int((census["resolved_coarse_class"][key])["n"])
        for key in ("Single-only", "Direct-only", "Sister-only", "Sister+direct")
    }
    coarse["Cross-class unresolved"] = int(census["cross_coarse_class"]["n"])
    topology_map = sample_record_map(authority["all7_summary"]["samples"])
    census_map = sample_record_map(authority["census_summary"]["samples"])
    sample_rows = []
    for item in samples:
        sample = item["sample"]
        t = topology_map[sample]
        c = census_map[sample]
        ranked = int(c["ranked_units"])
        unique = int(c["resolution"]["UNIQUE_TREE"]["n"])
        exact = int(c["one_exact_topology"]["n"])
        abstain = int(t["mutation_family_abstain_units"])
        sample_rows.append(
            f"""<tr>
<td><a class="sample-link" href="{esc(sample)}.html">{esc(sample)}</a>{'<br><span class="tag">same biological sample</span>' if sample in {'HCC1395','HCC1395_DORADO'} else ''}</td>
<td class="num">{int(t['groups_total']):,}</td><td class="num">{ranked:,}</td>
<td><b>{pct(unique,ranked):.1f}%</b><div class="barline"><span style="width:{pct(unique,ranked):.4f}%"></span></div></td>
<td class="num">{pct(exact,ranked):.1f}%</td><td class="num danger">{abstain:,}</td>
<td><a href="{esc(sample)}.html#genome">開啟 GRCh38 →</a></td></tr>"""
        )

    h1, h2 = census_map["HCC1395"], census_map["HCC1395_DORADO"]
    res_keys = ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY")
    hcc_similarity = js_divergence_similarity(
        [int(h1["resolution"][key]["n"]) for key in res_keys],
        [int(h2["resolution"][key]["n"]) for key in res_keys],
    )
    resolution_panel = distribution_html(
        resolution,
        int(census["ranked_units"]),
        RESOLUTION_COLORS,
        RESOLUTION_LABELS,
        res_keys,
    )
    coarse_panel = distribution_html(
        coarse,
        int(census["ranked_units"]),
        COARSE_COLORS,
        COARSE_LABELS,
        (
            "Single-only",
            "Direct-only",
            "Sister-only",
            "Sister+direct",
            "Cross-class unresolved",
        ),
    )
    evidence_paths = {
        "topology cohort receipt": authority["cohort_receipt_path"],
        "topology all7 summary": authority["all7_summary_path"],
        "signature census receipt": authority["census_receipt_path"],
        "signature census summary": authority["census_summary_path"],
        "k / HP funnel": authority["funnel_path"]
        or "derived from hash-bound per-sample MLHP receipts",
    }
    path_html = "".join(
        f'<span class="path"><b>{esc(label)}</b> <code>{esc(path)}</code></span>'
        for label, path in evidence_paths.items()
    )
    return f"""<!doctype html>
<html lang="zh-Hant"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<meta name="intersubmod-authority" content="{AUTHORITY_NAME}">
<meta name="intersubmod-ui-contract" content="{UI_CONTRACT}">
<meta name="intersubmod-cohort-receipt-sha256" content="{authority['cohort_receipt_sha256']}">
<meta name="intersubmod-census-receipt-sha256" content="{authority['census_receipt_sha256']}">
<title>Exact-PS layered workstation · cohort command center</title>
<style>{COMMON_CSS}</style></head>
<body>
<div class="authority"><b>PRIMARY · 2026-07-24 exact PS × HP</b><span>7 datasets / 6 biological samples · chr1–22 · strict endpoint read-linkage ≥3</span></div>
<header class="hero"><div class="shell">
  <div class="eyebrow">InterSubMod evidence atlas</div>
  <h1>從 read-linkage<br>到 exact topology</h1>
  <p class="lead">新的 layered workstation 以 exact phase-set、primary HP 與嚴格 read-linkage 重建區域；舊 50 kb proximity grouping 不再是預設資料。</p>
  <div class="boundary"><strong>最重要限制</strong><span>7/7 pipeline 是 technical PASS，但全 cohort 仍有 {int(totals['mutation_family_abstain_units']):,} 個 mutation-bearing groups 因資源 guard ABSTAIN，因此不是 topology-complete。所有比例均顯示明確分母。</span></div>
  <nav class="nav"><a href="#evidence">資料漏斗</a><a href="#topology">拓撲全貌</a><a href="#samples">7 組資料</a><a href="#pair">HCC1395 技術驗證</a><a href="#provenance">證據</a></nav>
</div></header>
<main>
<section class="section" id="evidence"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Authority funnel</div><h2>先分清楚「區域」與「可判定區域」</h2></div><p>source components 包含 k=1 singleton；source W 才是 k≥2 strict read-linked components。bounded block 與 final group 由證據切分，不使用 50 kb 傳遞合併。</p></div>
  <div class="metrics">
    <div class="metric"><span class="label">candidate sSNV universe</span><span class="value">469,849</span><span class="note">7 technical datasets</span></div>
    <div class="metric"><span class="label">source components</span><span class="value">{source_components:,}</span><span class="note">含 170,131 k=1 singleton</span></div>
    <div class="metric"><span class="label">final topology groups</span><span class="value">{int(totals['groups_total']):,}</span><span class="note">439,685 sSNV memberships</span></div>
    <div class="metric"><span class="label">ranked complete</span><span class="value">{int(totals['ranked_units']):,}</span><span class="note">{pct(totals['ranked_units'],totals['mutation_bearing_units']):.1f}% of mutation-bearing</span></div>
    <div class="metric"><span class="label">global-best trees enumerated</span><span class="value">{int(census['best_tree_total']):,}</span><span class="note">canonical score/tie reproduced</span></div>
  </div>
  <div class="funnel" style="margin-top:1rem">
    <div class="funnel-step"><span>source components</span><b>{source_components:,}</b><small>PS × HP graph components</small></div>
    <div class="funnel-step"><span>strict read-linked W</span><b>{source_w:,}</b><small>k≥2；singleton 不進 tree</small></div>
    <div class="funnel-step"><span>bounded blocks</span><b>{bounded:,}</b><small>證據邊界切分</small></div>
    <div class="funnel-step"><span>final groups</span><b>{int(totals['groups_total']):,}</b><small>新全基因工作站 universe</small></div>
  </div>
</div></section>
<section class="section" id="topology"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Cohort topology census</div><h2>88.26% 有單一 exact shape；不等於唯一 clone tree</h2></div><p>71,955 ranked-complete groups 為分母。舊 legacy 92.18% 口徑不能沿用。</p></div>
  <div class="grid-2">
    <article class="panel"><h3>🎯 Determinacy｜最佳樹同形到哪一層</h3><div class="denom">n={int(census['ranked_units']):,}</div>{resolution_panel}</article>
    <article class="panel"><h3>🌳 拓撲型態｜分子累積形狀</h3><div class="denom">n={int(census['ranked_units']):,}</div>{coarse_panel}</article>
  </div>
  <div class="callout" style="margin-top:1rem"><b>精確結論：</b>單一 rooted-unlabeled topology = {int(census['one_exact_topology']['n']):,}/{int(census['ranked_units']):,} = {float(census['one_exact_topology']['pct_ranked']):.4f}%；單一四類 coarse geometry = {int(census['one_coarse_class']['n']):,}/{int(census['ranked_units']):,} = {float(census['one_coarse_class']['pct_ranked']):.4f}%。兩者不能替代 cellular clone / ancestry 驗證。</div>
</div></section>
<section class="section" id="samples"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Seven technical datasets</div><h2>每組資料直接進入全基因視圖</h2></div><p>HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 technical datasets；因此 7 組資料只代表 6 個生物樣本。</p></div>
  <div class="table-wrap"><table><thead><tr><th>dataset</th><th>final groups</th><th>ranked</th><th>唯一最佳樹</th><th>單一 exact shape</th><th>ABSTAIN</th><th>工作站</th></tr></thead><tbody>{''.join(sample_rows)}</tbody></table></div>
</div></section>
<section class="section" id="pair"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Technical concordance</div><h2>HCC1395 × HCC1395_DORADO</h2></div><p>這是同一生物樣本在不同 basecalling / technical dataset 下的組成比較，不能算兩個獨立 biological replicates。</p></div>
  <div class="grid-2">
    <article class="panel"><h3>Resolution composition JS similarity</h3><div class="metric" style="padding:0;min-height:0"><span class="value">{100*hcc_similarity:.1f}%</span><span class="note">1 − Jensen–Shannon distance；只比較三類 determinacy composition</span></div></article>
    <article class="panel"><h3>判讀</h3><p>兩組資料的 determinacy 組成高度接近，可作技術穩健性訊號；但 PS 邊界、coverage 與 region universe 不相同，不能由這個單一比例宣稱逐區域 topology 完全一致。</p><p><a href="HCC1395.html">HCC1395 →</a>　<a href="HCC1395_DORADO.html">DORADO →</a></p></article>
  </div>
</div></section>
<section class="section" id="provenance"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Fail-closed provenance</div><h2>JSON 路徑不干擾主要數字</h2></div><p>原始 JSON、receipt 與 SHA-256 收在下方；builder 另逐列核對 MLHP/topology/census keys、row conservation 與 source identities。</p></div>
  <details><summary>機器證據與原始 JSON（預設收合）</summary><div class="details-body">{path_html}<span class="path"><b>cohort receipt SHA-256</b> <code>{authority['cohort_receipt_sha256']}</code></span><span class="path"><b>census receipt SHA-256</b> <code>{authority['census_receipt_sha256']}</code></span><span class="path"><b>verified pipeline bound identities</b> <code>{authority['pipeline_identity_count']}</code></span><span class="path"><b>verified census identities</b> <code>{authority['census_identity_count']}</code></span></div></details>
  <details><summary>歷史 50 kb baseline（不進入新統計）</summary><div class="details-body"><p>舊 canonical-v5 工作站用 proximity/transitive grouping 與不同 candidate-tree schema。它僅保留作方法演進說明；所有本頁 KPI、分布與樣本頁皆來自 2026-07-24 exact-PS authority。</p></div></details>
</div></section>
</main>
<footer class="footer"><div class="shell">InterSubMod · {AUTHORITY_NAME} · UI contract {UI_CONTRACT}</div></footer>
</body></html>"""


def build(index_only: bool = False, verify_only: bool = False) -> dict[str, Any]:
    authority = load_authority()
    topology_map = sample_record_map(authority["all7_summary"]["samples"])
    census_map = sample_record_map(authority["census_summary"]["samples"])
    funnel_samples = authority["funnel"]["samples"]
    require(set(funnel_samples) == set(EXPECTED_SAMPLES), "funnel sample set mismatch")

    samples = []
    for sample in EXPECTED_SAMPLES:
        item = load_sample(
            sample,
            authority,
            topology_map[sample],
            census_map[sample],
            funnel_samples[sample],
        )
        samples.append(item)
        if not verify_only and not index_only:
            atomic_write(OUTDIR / f"{sample}.html", sample_page(authority, item))

    require(
        sum(item["summary"]["groups"] for item in samples) == 98955,
        "all-page group count does not conserve 98,955",
    )
    require(
        sum(item["summary"]["ranked"] for item in samples) == 71955,
        "all-page ranked count does not conserve 71,955",
    )
    if not verify_only:
        if index_only:
            for sample in EXPECTED_SAMPLES:
                page = OUTDIR / f"{sample}.html"
                require(page.is_file(), f"index-only requires existing {page}")
                prefix = page.read_text(encoding="utf-8", errors="replace")[:12000]
                require(
                    f'<meta name="intersubmod-authority" content="{AUTHORITY_NAME}">'
                    in prefix,
                    f"stale sample page authority: {page}",
                )
                require(
                    f'<meta name="intersubmod-census-receipt-sha256" content="{authority["census_receipt_sha256"]}">'
                    in prefix,
                    f"stale census binding: {page}",
                )
        atomic_write(OUTDIR / "index.html", index_page(authority, samples))
    return {"authority": authority, "samples": samples}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index-only", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    try:
        result = build(index_only=args.index_only, verify_only=args.verify_only)
    except AuthorityError as exc:
        print(f"FAIL CLOSED: {exc}", file=sys.stderr)
        return 2
    groups = sum(item["summary"]["groups"] for item in result["samples"])
    ranked = sum(item["summary"]["ranked"] for item in result["samples"])
    print(f"OK exact-PS authority verified: samples=7 groups={groups} ranked={ranked}")
    print(f"INPUT_TOPOLOGY_ROOT={TOPOLOGY_ROOT}")
    print(f"INPUT_CENSUS_ROOT={CENSUS_ROOT}")
    if not args.verify_only:
        print(f"OUTPUT={OUTDIR / 'index.html'}")
        print(f"OUTPUT_SAMPLE_PAGES=7")
    print(
        "CENSUS_RECEIPT_SHA256="
        f"{result['authority']['census_receipt_sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
