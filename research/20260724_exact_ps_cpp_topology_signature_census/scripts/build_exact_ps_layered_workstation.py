#!/usr/bin/env python3
"""Build the exact-PS authority layered workstation.

The builder fails closed unless the 2026-07-24 all-seven exact-PS topology run,
the exact topology-signature census, and the k/HP funnel all pass their
hash/schema/arithmetic contracts.  It writes a cohort index plus seven
standalone sample pages to ``docs/methodology/_assets/layered_workstation``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import os
import re
import sys
from bisect import bisect_right
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
CANDIDATE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260725_exact_ps_candidate_factorization/all7_v2"
)
SIMILARITY_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260725_exact_ps_cohort_similarity/all7_v1/cohort_similarity.v1.json"
)
GENE_DRUG_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260725_exact_ps_gene_drug_annotation/all7_v2"
)
GENE_DRUG_PATH = GENE_DRUG_ROOT / "annotation.v1.json"
GENE_DRUG_RECEIPT_PATH = GENE_DRUG_ROOT / "receipt.json"
GENE_DRUG_SOURCE_PATHS = {
    "gencode": Path(
        "/big7_disk/liaoyoyo2001/gene_annotation/"
        "gencode.v46.basic.annotation.gtf.gz"
    ),
    "cgc": Path(
        "/big7_disk/liaoyoyo2001/gene_annotation/"
        "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz"
    ),
    "cosmic_genes": Path(
        "/big7_disk/liaoyoyo2001/gene_annotation/cosmic_v104/"
        "Cosmic_Genes_v104_GRCh38.tsv.gz"
    ),
    "dgidb": Path(
        "/big7_disk/liaoyoyo2001/gene_annotation/dgidb_interactions.tsv"
    ),
}
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
# Presentation aliases/order only; authority validation, hashes, sidecars, and
# filenames must continue to use EXPECTED_SAMPLES' internal keys.
DISPLAY_ORDER = (
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
    "H1437",
    "H2009",
    "COLO829",
)
DISPLAY_LABELS = {
    "HCC1395": "HCC1395_HKU",
    "HCC1395_DORADO": "HCC1395_NYGC",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
    "H1437": "H1437",
    "H2009": "H2009",
    "COLO829": "COLO829",
}
UI_CONTRACT = "layered-workstation-exact-ps-v3"
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
EXPECTED_GENE_DRUG_REGION_COUNTS = {
    "HCC1395": {"cgc": 480, "cgc_drug": 174},
    "HCC1395_DORADO": {"cgc": 240, "cgc_drug": 91},
    "HCC1937": {"cgc": 189, "cgc_drug": 66},
    "HCC1954": {"cgc": 273, "cgc_drug": 90},
    "H1437": {"cgc": 537, "cgc_drug": 167},
    "H2009": {"cgc": 1379, "cgc_drug": 503},
    "COLO829": {"cgc": 456, "cgc_drug": 161},
}


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


def load_gene_drug_authority() -> dict[str, Any]:
    receipt = load_json(GENE_DRUG_RECEIPT_PATH)
    require(
        receipt.get("schema_name")
        == "intersubmod.exact_ps_gene_drug_annotation.receipt"
        and receipt.get("schema_version") == "1.1.0"
        and receipt.get("all_pass") is True,
        "unexpected/non-PASS gene-drug annotation receipt",
    )
    all_true(receipt.get("checks") or {}, "gene_drug_receipt.checks")
    output = receipt.get("output") or {}
    require(
        output.get("path") == GENE_DRUG_PATH.name
        and int(output.get("size_bytes", -1)) == GENE_DRUG_PATH.stat().st_size
        and output.get("sha256") == sha256_file(GENE_DRUG_PATH),
        "gene-drug annotation output identity mismatch",
    )
    raw_sources = receipt.get("raw_sources") or {}
    require(
        set(raw_sources) == set(GENE_DRUG_SOURCE_PATHS),
        "gene-drug raw-source set mismatch",
    )
    for key, path in GENE_DRUG_SOURCE_PATHS.items():
        identity = raw_sources[key]
        require(
            path.is_file()
            and identity.get("file_name") == path.name
            and int(identity.get("size_bytes", -1)) == path.stat().st_size
            and identity.get("sha256") == sha256_file(path),
            f"gene-drug raw-source identity mismatch: {key}",
        )

    annotation = load_json(GENE_DRUG_PATH)
    require(
        annotation.get("schema_name")
        == "intersubmod.exact_ps_gene_drug_annotation"
        and annotation.get("schema_version") == "1.1.0",
        "unexpected gene-drug annotation schema",
    )
    scope = annotation.get("scope") or {}
    require(
        tuple(scope.get("chromosomes") or ()) == tuple(CHROM_LENGTHS)
        and scope.get("genome_build") == "GRCh38"
        and scope.get("coordinate_system")
        == "GENCODE_v46_gene_body_1_based_inclusive"
        and scope.get("partial") is False,
        "gene-drug annotation scope mismatch",
    )
    require(
        scope.get("intersection_contract")
        == (
            "position_in_GENCODE_gene_body AND same_gene_in_CGC; "
            "drug subset additionally requires DGIdb approved=TRUE AND "
            "anti_neoplastic=TRUE"
        ),
        "gene-drug intersection contract mismatch",
    )
    counts = annotation.get("counts") or {}
    require(
        int(counts.get("cgc_source_rows", -1)) == 768
        and int(counts.get("cgc_autosomal_rows", -1)) == 727
        and int(counts.get("cgc_genes", -1)) == 721
        and int(counts.get("cgc_with_approved_antineoplastic", -1)) == 279
        and int(counts.get("dgidb_source_rows", -1)) == 98239
        and int(counts.get("dgidb_approved_and_antineoplastic_rows", -1))
        == 9966,
        "gene-drug source/count regression failed",
    )
    require(
        int(counts.get("dgidb_selected_rows_indexed_by_symbol_only", -1))
        == 224
        and int(
            counts.get(
                "dgidb_symbol_only_rows_accepted_after_gencode_gate",
                -1,
            )
        )
        == 0
        and int(
            counts.get("dgidb_symbol_only_keys_overlapping_cgc", -1)
        )
        == 0
        and int(counts.get("dgidb_symbol_fallback_gene_matches", -1))
        == 0,
        "unsafe/free-text DGIdb symbol fallback entered production joins",
    )
    claim_ceiling = annotation.get("claim_ceiling") or {}
    require(
        claim_ceiling.get("level") == "gene_level_annotation_only"
        and "clinical_actionability"
        in (claim_ceiling.get("prohibited_claims") or []),
        "gene-drug claim ceiling is incomplete",
    )

    genes = annotation.get("genes") or []
    require(len(genes) == 721, "gene-drug gene row count mismatch")
    seen_gene_ids: set[str] = set()
    by_chrom: dict[str, list[Mapping[str, Any]]] = {
        chrom: [] for chrom in CHROM_LENGTHS
    }
    previous_key: tuple[int, int, int, str] | None = None
    for gene in genes:
        chrom = str(gene.get("chrom", ""))
        gene_id = str(gene.get("gene_id", ""))
        start = int(gene.get("start", -1))
        end = int(gene.get("end", -1))
        require(
            chrom in CHROM_LENGTHS
            and gene_id
            and gene_id not in seen_gene_ids
            and 1 <= start <= end
            and str(gene.get("hgnc_id", "")).startswith("HGNC:")
            and gene.get("gene_match_basis") == "cosg_to_hgnc",
            f"invalid gene-drug gene record: {gene_id or gene}",
        )
        key = (int(chrom[3:]), start, end, gene_id)
        require(
            previous_key is None or key > previous_key,
            f"gene-drug genes are not strictly sorted at {gene_id}",
        )
        previous_key = key
        seen_gene_ids.add(gene_id)
        for drug in gene.get("approved_antineoplastic_drugs") or []:
            require(
                str(drug.get("claim_name", "")).strip()
                and all(
                    str(source.get("name", "")).strip()
                    and str(source.get("version", "")).strip()
                    for source in drug.get("sources") or []
                ),
                f"invalid DGIdb drug claim at {gene_id}",
            )
        by_chrom[chrom].append(gene)
    interval_index = {
        chrom: {
            "genes": rows,
            "starts": [int(row["start"]) for row in rows],
        }
        for chrom, rows in by_chrom.items()
    }
    return {
        "annotation": annotation,
        "receipt": receipt,
        "interval_index": interval_index,
        "sha256": str(output["sha256"]),
        "receipt_sha256": sha256_file(GENE_DRUG_RECEIPT_PATH),
    }


def gene_drug_hits(
    chrom: str,
    positions: Iterable[int],
    interval_index: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    chromosome = interval_index.get(chrom) or {"genes": [], "starts": []}
    genes = chromosome["genes"]
    starts = chromosome["starts"]
    hits: dict[str, dict[str, Any]] = {}
    for position in sorted({int(value) for value in positions}):
        for gene in genes[: bisect_right(starts, position)]:
            if int(gene["end"]) < position:
                continue
            gene_id = str(gene["gene_id"])
            record = hits.setdefault(
                gene_id,
                {
                    "geneId": gene_id,
                    "symbol": str(gene["symbol"]),
                    "hgncId": str(gene["hgnc_id"]),
                    "start": int(gene["start"]),
                    "end": int(gene["end"]),
                    "tier": int(gene["cgc_tier"]),
                    "role": str(gene.get("cgc_role") or ""),
                    "tumour": str(gene.get("cgc_tumour") or ""),
                    "loci": [],
                    "drugs": [
                        {
                            "claimName": str(drug["claim_name"]),
                            "sources": [
                                {
                                    "name": str(source["name"]),
                                    "version": str(source["version"]),
                                }
                                for source in drug.get("sources") or []
                            ],
                        }
                        for drug in gene.get("approved_antineoplastic_drugs")
                        or []
                    ],
                },
            )
            record["loci"].append(position)
    return sorted(
        hits.values(),
        key=lambda row: (row["start"], row["end"], row["geneId"]),
    )


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

    candidate_receipt_path = CANDIDATE_ROOT / "receipt.json"
    candidate_receipt = load_json(candidate_receipt_path)
    require(
        candidate_receipt.get("schema_name")
        == "intersubmod.exact_ps_topology_candidate_factorization.receipt"
        and candidate_receipt.get("schema_version") == "1.1.0",
        "unexpected candidate-factorization receipt schema",
    )
    require(
        candidate_receipt.get("all_pass") is True
        and candidate_receipt.get("technical_all_pass") is True
        and candidate_receipt.get("validation_evidence_eligible") is False
        and candidate_receipt.get("all_mutation_bearing_families_complete")
        is False,
        "candidate-factorization cohort did not pass",
    )
    all_true(candidate_receipt.get("checks") or {}, "candidate_receipt.checks")
    candidate_scope = candidate_receipt.get("scope") or {}
    require(
        tuple(candidate_scope.get("datasets") or ()) == EXPECTED_SAMPLES
        and tuple(candidate_scope.get("chromosomes") or ())
        == tuple(CHROM_LENGTHS),
        "candidate-factorization all-seven/autosome scope mismatch",
    )
    candidate_totals = candidate_receipt.get("totals") or {}
    require(
        int(candidate_totals.get("ranked_units", -1)) == 71955
        and int(candidate_totals.get("minimum_trees", -1)) == 972592
        and int(candidate_totals.get("global_best_trees", -1)) == 680527,
        "candidate-factorization cohort totals mismatch",
    )
    candidate_identity_count = verify_nested_identities(
        {
            "upstream_authority": candidate_receipt.get("upstream_authority") or {},
            "implementation": candidate_receipt.get("implementation") or {},
            "samples": candidate_receipt.get("samples") or {},
        },
        "candidate_receipt",
    )
    require(
        candidate_identity_count >= 38,
        "candidate-factorization identity coverage too small",
    )

    similarity = load_json(SIMILARITY_PATH)
    require(
        similarity.get("schema_name") == "intersubmod.exact_ps_cohort_similarity"
        and similarity.get("schema_version") == "1.0.0",
        "unexpected cohort-similarity schema",
    )
    all_true(similarity.get("checks") or {}, "similarity.checks")
    similarity_scope = similarity.get("scope") or {}
    require(
        tuple(similarity_scope.get("datasets") or ()) == EXPECTED_SAMPLES
        and tuple(similarity_scope.get("chromosomes") or ())
        == tuple(CHROM_LENGTHS)
        and similarity_scope.get("partial") is False,
        "cohort-similarity all-seven/autosome scope mismatch",
    )
    require(
        len(similarity.get("pairwise_profiles") or []) == 21,
        "cohort-similarity pair count mismatch",
    )
    similarity_identity_count = verify_nested_identities(
        similarity.get("source_identities") or {},
        "similarity.source_identities",
    )
    require(
        similarity_identity_count >= 46,
        "cohort-similarity identity coverage too small",
    )
    hcc_profile = next(
        (
            row
            for row in similarity.get("pairwise_profiles") or []
            if row.get("pair_id") == "HCC1395__HCC1395_DORADO"
        ),
        None,
    )
    require(
        hcc_profile is not None
        and int(hcc_profile.get("similarity_rank_of_21", -1)) == 1
        and math.isclose(
            float(hcc_profile.get("profile_similarity", -1)),
            0.926381219637838,
            abs_tol=1e-15,
        ),
        "HCC technical-pair profile regression failed",
    )
    gene_drug_authority = load_gene_drug_authority()

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
        "candidate_receipt": candidate_receipt,
        "candidate_receipt_path": candidate_receipt_path,
        "candidate_receipt_sha256": sha256_file(candidate_receipt_path),
        "candidate_identity_count": candidate_identity_count,
        "similarity": similarity,
        "similarity_path": SIMILARITY_PATH,
        "similarity_sha256": sha256_file(SIMILARITY_PATH),
        "similarity_identity_count": similarity_identity_count,
        "gene_drug": gene_drug_authority["annotation"],
        "gene_drug_receipt": gene_drug_authority["receipt"],
        "gene_drug_interval_index": gene_drug_authority["interval_index"],
        "gene_drug_path": GENE_DRUG_PATH,
        "gene_drug_receipt_path": GENE_DRUG_RECEIPT_PATH,
        "gene_drug_sha256": gene_drug_authority["sha256"],
        "gene_drug_receipt_sha256": gene_drug_authority[
            "receipt_sha256"
        ],
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


def compact_candidate_factorization(
    factorization: Mapping[str, Any],
    topology: Mapping[str, Any],
    census: Mapping[str, Any],
) -> dict[str, Any]:
    require(
        factorization.get("schema_name")
        == "intersubmod.exact_ps_topology_candidate_factorization.unit"
        and factorization.get("schema_version") == "1.0.0",
        f"unexpected candidate schema at {topology.get('region_id')}",
    )
    for key in ("group_index", "region_id", "unit_id", "block_id"):
        require(
            str(factorization.get(key)) == str(topology.get(key)),
            f"candidate/topology {key} mismatch at {topology.get('region_id')}",
        )
    for key in (
        "active_bit_count",
        "minimum_vertex_set_count",
        "best_vertex_set_count",
    ):
        require(
            int(factorization.get(key, -1)) == int(topology.get(key, -2)),
            f"candidate/topology {key} mismatch at {topology.get('region_id')}",
        )
    for key in (
        "minimum_family_sha256",
        "total_tree_count",
        "best_score_fraction",
        "best_tree_tie_count",
    ):
        require(
            str(factorization.get(key)) == str(topology.get(key)),
            f"candidate/topology {key} mismatch at {topology.get('region_id')}",
        )
    require(
        [int(value) for value in factorization.get("active_original_indices") or []]
        == [int(value) for value in topology.get("active_original_indices") or []]
        and [int(value) for value in factorization.get("active_positions") or []]
        == [int(value) for value in topology.get("active_positions") or []],
        f"candidate/topology active loci mismatch at {topology.get('region_id')}",
    )
    require(
        factorization.get("resolution_class") == census.get("resolution_class"),
        f"candidate/census resolution mismatch at {topology.get('region_id')}",
    )
    require(
        factorization.get("canonical_reproduction_pass") is True
        and factorization.get("census_reproduction_pass") is True,
        f"candidate reproduction failed at {topology.get('region_id')}",
    )

    compact_sets: list[list[Any]] = []
    global_best_indices = []
    observed_vertices = [
        int(value) for value in factorization.get("observed_vertices") or []
    ]
    require(
        observed_vertices == sorted(set(observed_vertices)),
        f"candidate observed vertices invalid at {topology.get('region_id')}",
    )
    reproduced_all_edges: Counter[tuple[int, int]] = Counter()
    reproduced_best_edges: Counter[tuple[int, int]] = Counter()
    reproduced_all_edge_total = 0
    reproduced_best_edge_total = 0
    for index, item in enumerate(factorization.get("minimum_vertex_sets") or []):
        vertices = [int(value) for value in item.get("vertices") or []]
        require(
            vertices == sorted(set(vertices))
            and vertices
            and vertices[0] == 0
            and set(observed_vertices) <= set(vertices),
            f"candidate vertices invalid at {topology.get('region_id')}:{index}",
        )
        parent_rows = []
        child_values = set()
        tree_count = 1
        best_tree_count = 1
        for raw in item.get("parent_factorization") or []:
            require(
                isinstance(raw, list) and len(raw) == 3,
                f"candidate parent row invalid at {topology.get('region_id')}:{index}",
            )
            child = int(raw[0])
            legal = [int(value) for value in raw[1]]
            best = [int(value) for value in raw[2]]
            require(
                child in vertices
                and legal == sorted(set(legal))
                and best == sorted(set(best))
                and legal
                and best
                and set(best) <= set(legal)
                and all(
                    parent in vertices
                    and parent < child
                    and (parent & child) == parent
                    and bin(parent ^ child).count("1") == 1
                    for parent in legal
                ),
                f"candidate parent choices invalid at {topology.get('region_id')}:{index}",
            )
            require(
                child != 0 and child not in child_values,
                f"candidate child duplicate/root at {topology.get('region_id')}:{index}",
            )
            child_values.add(child)
            tree_count *= len(legal)
            best_tree_count *= len(best)
            parent_rows.append([child, legal, best])
        require(
            tree_count == int(item.get("total_tree_count", -1))
            and best_tree_count == int(item.get("best_tree_count", -1)),
            f"candidate parent product mismatch at {topology.get('region_id')}:{index}",
        )
        require(
            child_values == set(vertices) - {0},
            f"candidate child coverage mismatch at {topology.get('region_id')}:{index}",
        )
        is_global = item.get("is_global_best") is True
        require(
            is_global
            == (
                str(item.get("best_score_fraction"))
                == str(factorization.get("best_score_fraction"))
            ),
            f"candidate global-best score mismatch at {topology.get('region_id')}:{index}",
        )
        for child, legal, best in parent_rows:
            for parent in legal:
                reproduced_all_edges[(parent, child)] += tree_count // len(legal)
            if is_global:
                for parent in best:
                    reproduced_best_edges[(parent, child)] += (
                        best_tree_count // len(best)
                    )
        reproduced_all_edge_total += tree_count * (len(vertices) - 1)
        if is_global:
            reproduced_best_edge_total += best_tree_count * (len(vertices) - 1)
        if is_global:
            global_best_indices.append(index)
        compact_sets.append(
            [
                vertices,
                str(item.get("best_score_fraction")),
                str(item.get("total_tree_count")),
                str(item.get("best_tree_count")),
                1 if is_global else 0,
                parent_rows,
            ]
        )
    require(
        len(compact_sets) == int(factorization["minimum_vertex_set_count"])
        and len(global_best_indices) == int(factorization["best_vertex_set_count"]),
        f"candidate set conservation failed at {topology.get('region_id')}",
    )

    def compact_edges(key: str) -> list[list[Any]]:
        rows = []
        seen = set()
        for raw in factorization.get(key) or []:
            require(
                isinstance(raw, list) and len(raw) == 3,
                f"{key} row invalid at {topology.get('region_id')}",
            )
            parent, child, count = int(raw[0]), int(raw[1]), str(raw[2])
            require(
                count.isdigit()
                and int(count) > 0
                and (parent, child) not in seen,
                f"{key} count/duplicate invalid at {topology.get('region_id')}",
            )
            seen.add((parent, child))
            rows.append([parent, child, count])
        return rows

    all_edges = compact_edges("all_minimum_tree_edge_incidence")
    best_edges = compact_edges("global_best_edge_incidence")
    require(
        Counter(
            {
                (int(parent), int(child)): int(count)
                for parent, child, count in all_edges
            }
        )
        == reproduced_all_edges
        and Counter(
            {
                (int(parent), int(child)): int(count)
                for parent, child, count in best_edges
            }
        )
        == reproduced_best_edges
        and sum(int(row[2]) for row in all_edges)
        == reproduced_all_edge_total
        and sum(int(row[2]) for row in best_edges)
        == reproduced_best_edge_total,
        f"candidate edge incidence mismatch at {topology.get('region_id')}",
    )
    all_edge_keys = {(row[0], row[1]) for row in all_edges}
    require(
        {(row[0], row[1]) for row in best_edges} <= all_edge_keys,
        f"global-best edge is outside candidate union at {topology.get('region_id')}",
    )
    representative_vertices = sorted(
        int(item["vertex"])
        for item in topology.get("representative_best_vertices") or []
    )
    representative_edges = {
        (int(item["parent_vertex"]), int(item["child_vertex"]))
        for item in topology.get("representative_best_edges") or []
    }
    selected_set_indices = [
        index
        for index in global_best_indices
        if compact_sets[index][0] == representative_vertices
    ]
    require(
        len(selected_set_indices) == 1
        and representative_edges <= {(row[0], row[1]) for row in best_edges},
        f"representative is absent from global-best factorization at {topology.get('region_id')}",
    )
    signatures = [
        [
            str(item["shape_signature"]),
            str(item["coarse_class"]),
            str(item["tree_count"]),
            [int(value) for value in item["exemplar_vertices"]],
            [
                [int(edge[0]), int(edge[1])]
                for edge in item["exemplar_edges"]
            ],
        ]
        for item in factorization.get("global_best_signatures") or []
    ]
    require(
        sum(int(item[2]) for item in signatures)
        == int(factorization["best_tree_tie_count"]),
        f"candidate signature conservation failed at {topology.get('region_id')}",
    )
    return {
        "minimumSetCount": int(factorization["minimum_vertex_set_count"]),
        "totalTreeCount": str(factorization["total_tree_count"]),
        "bestSetCount": int(factorization["best_vertex_set_count"]),
        "bestTreeCount": str(factorization["best_tree_tie_count"]),
        "bestScore": str(factorization["best_score_fraction"]),
        "selectedSetIndex": selected_set_indices[0],
        "observedVertices": observed_vertices,
        "sets": compact_sets,
        "allEdges": all_edges,
        "bestEdges": best_edges,
        "signatures": signatures,
    }


def compact_region(
    group: Mapping[str, Any],
    topology: Mapping[str, Any],
    census: Mapping[str, Any] | None,
    candidate: Mapping[str, Any] | None,
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
        "activeIndices": [
            int(value) for value in topology.get("active_original_indices") or []
        ],
        "activePositions": [
            int(value) for value in topology.get("active_positions") or []
        ],
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
        "minimumSetCount": int(topology.get("minimum_vertex_set_count") or 0),
        "totalTreeCount": str(topology.get("total_tree_count") or "0"),
        "vertexSetCount": int(topology.get("best_vertex_set_count", 0)),
        "familySha256": str(topology.get("minimum_family_sha256") or ""),
        "recurrenceCompatibleSets": int(
            topology.get("recurrence_compatible_vertex_set_count") or 0
        ),
        "recurrenceIncompatibleSets": int(
            topology.get("recurrence_incompatible_vertex_set_count") or 0
        ),
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
        "candidate": compact_candidate_factorization(
            candidate, topology, census
        )
        if candidate is not None and census is not None
        else None,
    }


def load_target_boundary_evidence(mlhp: Mapping[str, Any]) -> dict[str, Any]:
    chromosome_inputs = (
        ((mlhp.get("inputs") or {}).get("chromosomes") or {}).get("chr10") or {}
    )
    receipt_path = verify_identity(
        chromosome_inputs.get("strict_endpoint_receipt") or {},
        "HCC1395.chr10.strict_endpoint_receipt",
    )
    receipt = load_json(receipt_path)
    require(
        receipt.get("schema_name") == "intersubmod.strict_exact_ps_hp_regions"
        and receipt.get("schema_version") == "1.1.0"
        and receipt.get("all_pass") is True,
        "HCC1395 chr10 strict-endpoint authority is not PASS",
    )
    all_true(receipt.get("checks") or {}, "HCC1395.chr10.strict_endpoint.checks")
    outputs = receipt.get("outputs") or {}
    edge_path = verify_identity(
        outputs.get("edges") or {}, "HCC1395.chr10.strict_endpoint.edges"
    )
    membership_path = verify_identity(
        outputs.get("membership") or {},
        "HCC1395.chr10.strict_endpoint.membership",
    )
    target_edge = None
    with gzip.open(edge_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row.get("phase_set") == "87549798"
                and row.get("hp_family") == "1"
                and int(row.get("pos_i1", -1)) == 87818272
                and int(row.get("pos_j1", -1)) == 87840023
            ):
                require(target_edge is None, "duplicate target strict endpoint edge")
                target_edge = row
    require(target_edge is not None, "target strict endpoint edge is missing")
    support = {
        key: int(target_edge[f"support_{key}"])
        for key in ("total", "RR", "RA", "AR", "AA")
    }
    require(
        support == {"total": 8, "RR": 6, "RA": 0, "AR": 0, "AA": 2}
        and target_edge.get("passes_primary_threshold") == "true",
        f"target strict endpoint support drifted: {support}",
    )
    target_membership = {}
    with gzip.open(membership_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row.get("phase_set") == "87549798"
                and row.get("linkage_basis") == "PS_HP1"
                and int(row.get("threshold", -1)) == 3
                and int(row.get("pos1", -1))
                in {87818272, 87840023, 87888228, 87928739}
            ):
                position = int(row["pos1"])
                require(
                    position not in target_membership,
                    f"duplicate target membership at {position}",
                )
                target_membership[position] = row
    require(
        set(target_membership)
        == {87818272, 87840023, 87888228, 87928739},
        "target strict membership set drifted",
    )
    require(
        all(
            target_membership[position].get("tree_eligible") == "true"
            and target_membership[position].get("component_class")
            == "READ_LINKED_MULTISITE_REGION"
            for position in (87818272, 87840023)
        )
        and all(
            target_membership[position].get("tree_eligible") == "false"
            and target_membership[position].get("inference_role")
            == "ABSTAIN_SINGLETON_UNLINKED"
            for position in (87888228, 87928739)
        ),
        "target strict component roles drifted",
    )
    return {
        "queryWindow": [87818272, 87928739],
        "includedPositions": [87818272, 87840023],
        "excludedSingletons": [87888228, 87928739],
        "threshold": 3,
        "endpointSupport": support,
        "strictReceipt": str(receipt_path),
        "strictReceiptSha256": sha256_file(receipt_path),
        "edgeSource": str(edge_path),
        "membershipSource": str(membership_path),
        "interpretation": (
            "8 molecules fixed at both endpoints (RR=6, AA=2) establish "
            "the physical read bridge. AA=2 is below MLHP MINREAD=3 as a "
            "genotype population, so it establishes connectivity but is "
            "not a solver-visible full pattern."
        ),
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
    candidate_path = CANDIDATE_ROOT / f"{sample}.candidate_factorization.jsonl"
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
    target_boundary = (
        load_target_boundary_evidence(mlhp) if sample == "HCC1395" else None
    )
    groups = mlhp.get("groups") or []
    topology_rows = list(iter_jsonl(topology_path))
    census_rows = list(iter_jsonl(census_path))
    candidate_iterator = iter(iter_jsonl(candidate_path))

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
    candidate_count = 0
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
            try:
                candidate = next(candidate_iterator)
            except StopIteration as exc:
                raise AuthorityError(
                    f"{sample} candidate sidecar ended before ranked row {index}"
                ) from exc
            require(
                int(candidate.get("group_index", -1)) == index,
                f"{sample} candidate order mismatch at {index}",
            )
            candidate_count += 1
            ranked_indices.add(index)
        else:
            require(census is None, f"{sample} unranked row unexpectedly in census: {index}")
            candidate = None
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
        compact = compact_region(group, topology, census, candidate)
        annotations = gene_drug_hits(
            compact["chrom"],
            compact["positions"],
            authority["gene_drug_interval_index"],
        )
        compact["geneDrug"] = annotations
        compact["hasCgc"] = bool(annotations)
        compact["hasCgcDrug"] = any(
            bool(hit["drugs"]) for hit in annotations
        )
        active_positions = {int(value) for value in compact["activePositions"]}
        compact["hasCgcActive"] = any(
            any(int(position) in active_positions for position in hit["loci"])
            for hit in annotations
        )
        compact["hasCgcDrugActive"] = any(
            bool(hit["drugs"])
            and any(
                int(position) in active_positions for position in hit["loci"]
            )
            for hit in annotations
        )
        compact["annotationClass"] = (
            "CGC_DRUG"
            if compact["hasCgcDrug"]
            else "CGC_ONLY"
            if compact["hasCgc"]
            else "NO_CGC"
        )
        rows.append(compact)

    try:
        extra_candidate = next(candidate_iterator)
    except StopIteration:
        extra_candidate = None
    require(
        extra_candidate is None and candidate_count == expected_ranked,
        f"{sample} candidate sidecar row count mismatch",
    )

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
    cgc_regions = sum(bool(row["hasCgc"]) for row in rows)
    cgc_drug_regions = sum(bool(row["hasCgcDrug"]) for row in rows)
    cgc_active_regions = sum(bool(row["hasCgcActive"]) for row in rows)
    cgc_drug_active_regions = sum(
        bool(row["hasCgcDrugActive"]) for row in rows
    )
    expected_annotation = EXPECTED_GENE_DRUG_REGION_COUNTS[sample]
    annotation_partition = Counter(row["annotationClass"] for row in rows)
    require(
        cgc_regions == expected_annotation["cgc"]
        and cgc_drug_regions == expected_annotation["cgc_drug"]
        and annotation_partition["CGC_DRUG"] == cgc_drug_regions
        and annotation_partition["CGC_ONLY"]
        == cgc_regions - cgc_drug_regions
        and annotation_partition["NO_CGC"] == len(rows) - cgc_regions
        and sum(annotation_partition.values()) == len(rows),
        f"{sample} exact-locus gene/drug count regression failed: "
        f"CGC={cgc_regions}, CGC∩drug={cgc_drug_regions}",
    )
    cgc_gene_ids = {
        hit["geneId"] for row in rows for hit in row["geneDrug"]
    }
    cgc_drug_gene_ids = {
        hit["geneId"]
        for row in rows
        for hit in row["geneDrug"]
        if hit["drugs"]
    }
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
            require(
                target["annotationClass"] == "NO_CGC"
                and target["geneDrug"] == [],
                "HCC1395 historical PTEN envelope leaked into exact-locus annotation",
            )
            target["boundaryContext"] = target_boundary
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
        "cgcRegions": cgc_regions,
        "cgcDrugRegions": cgc_drug_regions,
        "cgcActiveRegions": cgc_active_regions,
        "cgcDrugActiveRegions": cgc_drug_active_regions,
        "cgcGenes": len(cgc_gene_ids),
        "cgcDrugGenes": len(cgc_drug_gene_ids),
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
            "candidateFactorization": str(candidate_path),
            "pipelineReceipt": str(sample_root / "pipeline_receipt.json"),
            "censusReceipt": str(authority["census_receipt_path"]),
            "candidateReceipt": str(authority["candidate_receipt_path"]),
            "geneDrugAnnotation": str(authority["gene_drug_path"]),
            "geneDrugReceipt": str(authority["gene_drug_receipt_path"]),
            **(
                {
                    "targetStrictEndpointReceipt": str(
                        target_boundary["strictReceipt"]
                    ),
                    "targetStrictEndpointEdges": str(
                        target_boundary["edgeSource"]
                    ),
                    "targetStrictMembership": str(
                        target_boundary["membershipSource"]
                    ),
                }
                if target_boundary
                else {}
            ),
        },
        "hashes": {
            "mlhp": sha256_file(mlhp_path),
            "topology": sha256_file(topology_path),
            "census": sha256_file(census_path),
            "candidateFactorization": sha256_file(candidate_path),
            "geneDrugAnnotation": authority["gene_drug_sha256"],
            "geneDrugReceipt": authority["gene_drug_receipt_sha256"],
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
.panel,.browser>*,.detail-grid>*,.candidate-grid>*{min-width:0}
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
.annotation-overview{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1px;margin-top:1rem;border:1px solid var(--line);background:var(--line)}
.annotation-overview>div{min-width:0;padding:.8rem;background:#f7f3ef}
.annotation-overview span{display:block;color:var(--muted);font-size:.68rem}
.annotation-overview b{display:block;font:800 1rem/1.35 ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere}
.annotation-overview .annotation-contract{grid-column:span 2;background:#fff1ee}
.annotation-overview .annotation-contract b{font-family:inherit;font-size:.78rem;color:#7f2f25}
.cohort-dashboard{display:grid;gap:1rem}
.cohort-panel{min-width:0;background:var(--panel);border:1px solid var(--line);box-shadow:var(--shadow);padding:1.1rem}
.cohort-panel-head{display:flex;justify-content:space-between;gap:1rem;align-items:start;margin-bottom:.9rem}
.cohort-panel-head p{max-width:720px;color:var(--muted);font-size:.8rem;margin:0}
.compare-tabs{display:flex;flex-wrap:wrap;gap:.35rem}
.compare-btn{min-height:40px;padding:.42rem .7rem;border:1px solid var(--line);background:var(--panel);cursor:pointer}
.compare-btn[aria-pressed=true]{background:var(--ink);border-color:var(--ink);color:#fff}
.profile-chart{display:grid;gap:.65rem}
.profile-row{display:grid;grid-template-columns:minmax(150px,.55fr) minmax(280px,2.2fr) minmax(150px,.65fr);gap:.75rem;align-items:center}
.profile-sample{font-weight:800;overflow-wrap:anywhere}
.profile-sample small{display:block;color:var(--muted);font-size:.65rem;font-weight:700;letter-spacing:.04em;text-transform:uppercase}
.profile-stack{display:flex;height:30px;overflow:hidden;border:1px solid rgba(0,0,0,.2);background:var(--soft)}
.profile-stack span{display:block;min-width:1px;border-right:1px solid rgba(255,255,255,.32)}
.profile-dominant{font-size:.72rem;color:var(--muted);text-align:right}
.profile-dominant b{display:block;color:var(--ink);font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.82rem}
.profile-legend{display:flex;flex-wrap:wrap;gap:.3rem .75rem;margin-top:.9rem;padding-top:.75rem;border-top:1px solid var(--line);font-size:.72rem;color:var(--muted)}
.profile-legend span{display:inline-flex;align-items:center;gap:.35rem}
.profile-legend i{width:.72rem;height:.72rem;display:inline-block;border:1px solid rgba(0,0,0,.12)}
.heatmap-layout{display:grid;grid-template-columns:minmax(0,1.45fr) minmax(280px,.55fr);gap:1rem;align-items:start}
.heatmap-layout>*,.hcc-scoreboard>*{min-width:0}
.heatmap-scroll{overflow:auto;border:1px solid var(--line);background:var(--panel)}
.heatmap-table{min-width:840px;table-layout:fixed}
.heatmap-table th{position:static;text-transform:none;letter-spacing:0;font-size:.68rem;line-height:1.15;white-space:nowrap}
.heatmap-table th:first-child{width:150px}
.heatmap-table td{padding:.55rem .35rem;text-align:center;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.76rem;font-weight:800}
.heatmap-table td.technical-pair{outline:3px solid var(--amber);outline-offset:-3px}
.heatmap-caption{margin-top:.55rem;color:var(--muted);font-size:.73rem}
.pair-ranking{margin:0;padding:0;list-style:none;border:1px solid var(--line)}
.pair-rank{display:grid;grid-template-columns:2rem minmax(0,1fr) auto;gap:.65rem;align-items:center;padding:.68rem;border-bottom:1px solid var(--line)}
.pair-rank:last-child{border-bottom:0}
.pair-rank .rank{font:800 1rem/1 ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--muted)}
.pair-rank .pair{font-size:.76rem;font-weight:800;overflow-wrap:anywhere}
.pair-rank .pair small{display:block;margin-top:.15rem;color:var(--muted);font-size:.64rem;font-weight:600}
.pair-rank .score{font:800 .88rem/1 ui-monospace,SFMono-Regular,Menlo,monospace}
.pair-rank.technical{background:#fff7e7}
.hcc-scoreboard{display:grid;grid-template-columns:minmax(240px,.45fr) minmax(0,1.55fr);gap:1rem;align-items:stretch}
.hcc-score{display:flex;flex-direction:column;justify-content:space-between;padding:1.25rem;background:#123c34;color:#fff}
.hcc-score .label{color:#d5e8e1;font-size:.74rem;letter-spacing:.08em;text-transform:uppercase}
.hcc-score .value{display:block;margin:.5rem 0;font:800 clamp(2rem,4vw,3.4rem)/.95 ui-monospace,SFMono-Regular,Menlo,monospace;letter-spacing:-.04em}
.hcc-score .interval{font:700 .82rem/1.5 ui-monospace,SFMono-Regular,Menlo,monospace}
.hcc-score .note{margin-top:1rem;color:#d5e8e1;font-size:.72rem}
.evidence-ladder{display:grid;gap:.55rem;padding:.25rem 0}
.ladder-step{display:grid;grid-template-columns:minmax(180px,.85fr) minmax(140px,1.15fr) auto;gap:.7rem;align-items:center}
.ladder-label{font-size:.76rem;font-weight:800}
.ladder-label small{display:block;color:var(--muted);font-size:.66rem;font-weight:500}
.ladder-track{height:10px;background:var(--soft);border:1px solid rgba(0,0,0,.1)}
.ladder-track span{display:block;height:100%;background:var(--accent)}
.ladder-value{min-width:104px;text-align:right;font:800 .76rem/1 ui-monospace,SFMono-Regular,Menlo,monospace}
.claim-boundary{margin-top:1rem;padding:1rem 1.1rem;border-left:5px solid var(--danger);background:#fff1ee}
.claim-boundary strong{display:block;margin-bottom:.25rem;color:#7f2f25}
.controls{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:1rem;align-items:end;margin:1rem 0}
.mode-buttons,.legend{display:flex;flex-wrap:wrap;gap:.35rem}
.mode-btn,.legend-btn,.action{min-height:40px;padding:.42rem .7rem;border:1px solid var(--line);background:var(--panel);cursor:pointer}
.mode-btn[aria-pressed=true]{background:var(--ink);color:#fff;border-color:var(--ink)}
.annotation-filters{display:flex;flex-wrap:wrap;gap:.35rem;margin-top:.75rem}
.annotation-filter{min-height:40px;padding:.42rem .7rem;border:1px solid var(--line);background:var(--panel);cursor:pointer}
.annotation-filter[aria-pressed=true]{background:#7f2f25;border-color:#7f2f25;color:#fff}
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
.detail-grid{display:grid;grid-template-columns:minmax(0,1fr);gap:1rem}
.detail-section{margin-top:1rem;padding:1rem;border:1px solid var(--line);background:var(--panel)}
.detail-section>h3{font-size:1.08rem}
.section-number{display:inline-grid;place-items:center;width:1.65rem;height:1.65rem;margin-right:.45rem;border-radius:50%;background:var(--ink);color:#fff;font:800 .78rem/1 ui-monospace,SFMono-Regular,Menlo,monospace}
.locus-strip{display:grid;grid-auto-flow:column;grid-auto-columns:minmax(148px,1fr);gap:1px;overflow-x:auto;margin-top:.7rem;border:1px solid var(--line);background:var(--line)}
.locus-card{position:relative;min-height:142px;padding:.72rem;background:#f6f4ec}
.locus-card.active{background:#e8f2ed}
.locus-card.cancer-locus{box-shadow:inset 0 0 0 2px var(--purple)}
.locus-card.drug-locus{box-shadow:inset 0 0 0 3px var(--danger)}
.locus-card .site{display:flex;justify-content:space-between;gap:.5rem;font-size:.72rem;color:var(--muted)}
.locus-card .position{display:block;margin:.2rem 0 .45rem;font:800 .86rem/1.25 ui-monospace,SFMono-Regular,Menlo,monospace}
.locus-card .af{font:800 1.2rem/1 ui-monospace,SFMono-Regular,Menlo,monospace}
.locus-card .reads{display:block;margin-top:.35rem;font-size:.7rem;color:var(--muted)}
.af-track{height:6px;margin-top:.55rem;background:#d7d6cd}
.af-track>span{display:block;height:100%;background:var(--accent)}
.locus-genes{display:flex;flex-wrap:wrap;gap:.2rem;margin-top:.5rem}
.locus-gene{display:inline-block;padding:.12rem .35rem;border:1px solid var(--purple);color:#4f3b74;background:#f3effa;font-size:.62rem;font-weight:800}
.locus-gene.drug{border-color:var(--danger);color:#7f2f25;background:#fff1ee}
.gene-drug-panel{border-left:5px solid var(--danger);background:#fffaf8}
.gene-drug-panel>summary{display:flex;align-items:center;gap:.45rem;flex-wrap:wrap}
.gene-drug-table td:first-child{font-weight:800}
.gene-drug-table .drug-list{max-width:620px;line-height:1.7;overflow-wrap:anywhere}
.gene-drug-table .source-list{display:block;margin-top:.25rem;color:var(--muted);font-size:.68rem}
.claim-ceiling{padding:.75rem .85rem;border-left:4px solid var(--amber);background:#fff7e7;color:#62451f;font-size:.75rem}
.matrix-wrap{overflow:auto;border:1px solid var(--line);background:var(--panel)}
.state-matrix{width:max-content;min-width:100%;border-collapse:separate;border-spacing:0}
.state-matrix th,.state-matrix td{min-width:44px;padding:.4rem;text-align:center}
.state-matrix th:first-child,.state-matrix td:first-child{position:sticky;left:0;z-index:3;min-width:88px;text-align:left;background:#f1efe7}
.state-matrix th:last-child,.state-matrix td:last-child{min-width:72px}
.state-cell{font:800 .78rem/1 ui-monospace,SFMono-Regular,Menlo,monospace}
.state-cell.A{background:#cbe7da;color:#0f5f4c}
.state-cell.R{background:#eceae1;color:#4d5751}
.state-cell.X{background:#fff2d8;color:#8b581b}
.candidate-summary{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:1px;margin:.75rem 0;background:var(--line);border:1px solid var(--line)}
.candidate-summary>div{padding:.7rem;background:#f7f5ed}
.candidate-summary span{display:block;color:var(--muted);font-size:.68rem}
.candidate-summary b{display:block;font:800 1rem/1.25 ui-monospace,SFMono-Regular,Menlo,monospace}
.candidate-toolbar{display:flex;align-items:center;gap:.45rem;flex-wrap:wrap;margin:.6rem 0}
.candidate-toolbar input{width:7rem;min-height:40px;padding:.4rem;border:1px solid var(--line);background:#fff}
.candidate-network{min-height:320px;border:1px solid var(--line);background:#fbfaf4;overflow:auto}
.candidate-network svg{display:block;min-width:620px;width:100%;height:auto;min-height:320px}
.candidate-grid{display:grid;grid-template-columns:minmax(0,1fr);gap:1rem}
.edge-key{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.edge-selected td{background:#e8f2ed}
.edge-global-best:not(.edge-selected) td{background:#edf2f7}
.edge-forced{color:var(--accent)}
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
  .browser,.detail-grid,.candidate-grid{grid-template-columns:1fr}
  .heatmap-layout,.hcc-scoreboard{grid-template-columns:1fr}
  .region-list{max-height:360px}
  .evidence-grid{grid-template-columns:repeat(2,minmax(0,1fr))}
  .candidate-summary{grid-template-columns:repeat(3,minmax(0,1fr))}
  .annotation-overview{grid-template-columns:repeat(2,minmax(0,1fr))}
  .annotation-overview .annotation-contract{grid-column:span 2}
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
  .candidate-summary{grid-template-columns:repeat(2,minmax(0,1fr))}
  .annotation-overview{grid-template-columns:1fr}
  .annotation-overview .annotation-contract{grid-column:auto}
  .cohort-panel-head{display:block}
  .cohort-panel-head .compare-tabs{margin-top:.75rem}
  .profile-row{grid-template-columns:minmax(115px,.65fr) minmax(190px,1.35fr)}
  .profile-dominant{grid-column:2;text-align:left}
  .ladder-step{grid-template-columns:1fr auto}
  .ladder-track{grid-column:1/-1;grid-row:2}
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
ANNOTATION_COLORS = {
    "CGC_DRUG": "#a94336",
    "CGC_ONLY": "#6b5592",
    "NO_CGC": "#b9b8ad",
}
ANNOTATION_LABELS = {
    "CGC_DRUG": "癌症基因 ∩ approved 抗腫瘤藥",
    "CGC_ONLY": "癌症基因・無 approved 抗腫瘤藥交集",
    "NO_CGC": "已評估・無癌症基因 locus 命中",
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
    cgc_regions = int(summary.get("cgcRegions", 0))
    cgc_drug_regions = int(summary.get("cgcDrugRegions", 0))
    cgc_active_regions = int(summary.get("cgcActiveRegions", 0))
    cgc_drug_active_regions = int(summary.get("cgcDrugActiveRegions", 0))
    cgc_genes = int(summary.get("cgcGenes", 0))
    cgc_drug_genes = int(summary.get("cgcDrugGenes", 0))

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
            "annotation": {
                key: {"label": ANNOTATION_LABELS[key], "color": color}
                for key, color in ANNOTATION_COLORS.items()
            },
        },
    }
    display_label = DISPLAY_LABELS[sample]
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
    <article class="panel"><h3>新版 tree input</h3><p><b>{87818272, 87840023}</b> · exact PS=87549798 · HP1 · fixed-endpoint bridge support=8（RR=6、AA=2）。結果為 <b>UNIQUE_TREE / Direct-only</b>，selected path 是 ROOT → H_AR → H_AA。</p><p class="denom">AA=2 低於 MLHP MINREAD=3，因此只支持區域連通，未成為 solver-visible full population；候選 family 仍由 partial ALT constraints 求解。</p><p><a href="#genome" onclick="setTimeout(()=>{const q=document.getElementById('regionSearch');q.value='chr10:87818272-87928739';document.getElementById('searchForm').requestSubmit()},50)">在全基因圖定位 →</a></p></article>
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
<meta name="intersubmod-gene-drug-sha256" content="{authority['gene_drug_sha256']}">
<meta name="intersubmod-cgc-locus-regions" content="{cgc_regions}">
<meta name="intersubmod-cgc-drug-locus-regions" content="{cgc_drug_regions}">
<title>Exact-PS 分層拓撲工作站 · {esc(display_label)}</title>
<style>{COMMON_CSS}</style>
</head>
<body>
<div class="authority"><b>PRIMARY · 2026-07-24 exact PS × HP</b><span>strict endpoint read-linkage ≥3 · chr1–22 · receipt-bound</span></div>
<header class="hero">
  <div class="shell">
    <a class="back" href="index.html">← cohort 總覽</a>
    <div class="eyebrow">Regional mathematical topology · {esc(display_label)}</div>
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
  <div class="annotation-overview" data-testid="annotation-overview">
    <div><span>CGC locus-hit regions</span><b>{cgc_regions:,}</b></div>
    <div><span>CGC ∩ approved antineoplastic</span><b>{cgc_drug_regions:,}</b></div>
    <div><span>命中 CGC genes</span><b>{cgc_genes:,}</b></div>
    <div><span>同基因亦有 approved 抗腫瘤藥</span><b>{cgc_drug_genes:,}</b></div>
    <div class="annotation-contract"><span>交集契約</span><b>actual sSNV locus ∈ GENCODE v46 gene body；同一 HGNC gene 同時進入 COSMIC CGC v104 與 DGIdb approved=true ∧ anti_neoplastic=true。這是資料庫脈絡，不是拓撲排序或臨床可用性。</b></div>
  </div>
  <p class="denom">分母 = {summary['groups']:,} final groups；主計數納入 group 內全部 exact sSNV loci（active 與 non-active）。只看 active topology loci 時，CGC={cgc_active_regions:,}、同基因 CGC∩drug={cgc_drug_active_regions:,}。</p>
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
      <span class="denom" id="modeButtonsLabel">著色維度</span>
      <div class="mode-buttons" id="modeButtons" data-testid="mode-buttons" role="group" aria-labelledby="modeButtonsLabel"></div>
      <span class="denom" id="annotationFiltersLabel" style="display:block;margin-top:.7rem">癌症基因／藥物快速篩選（與目前拓撲類別做 AND）</span>
      <div class="annotation-filters" id="annotationFilters" data-testid="annotation-filters" role="group" aria-labelledby="annotationFiltersLabel"></div>
    </div>
    <form class="search" id="searchForm">
      <label class="sr-only" for="regionSearch">搜尋座標</label>
      <input id="regionSearch" placeholder="chr10:87818272-87928739" autocomplete="off">
      <button class="action" type="submit">搜尋</button>
      <button class="action" type="button" id="clearSearch">清除</button>
    </form>
  </div>
  <div class="legend" id="legend" data-testid="legend" role="group" aria-label="目前著色維度的類別圖例；可多選聯集"></div>
  <div class="genome-wrap">
    <canvas id="genomeCanvas" data-testid="genome-canvas" role="img" aria-label="GRCh38 chr1 到 chr22 拓撲區域分布；下方 region 清單提供鍵盤可操作的等價內容"></canvas>
    <div class="genome-status"><span id="filterStatus" aria-live="polite"></span><span>座標按 GRCh38 chromosome bp 比例；點大小僅為可見性，不代表區域長度。</span></div>
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
    <article class="panel flat"><h3>代表樹與完整候選</h3><p>完整 edge union 來自 v2 factorization sidecar：逐 minimum vertex set 保存 legal parent choices，並重算所有 minimum tree 與 AF-global-best tree 的 edge incidence。藍色是 global-best edge union；綠色只是一棵 deterministic selected exemplar。edge incidence 是候選樹計數，不是 read support 或 probability。</p></article>
    <article class="panel flat"><h3>癌症基因 ∩ 藥物交集</h3><p>只有實際 sSNV locus 落入 GENCODE v46 gene body 才算 gene hit；再要求同一 HGNC gene 是 COSMIC CGC v104 且 DGIdb 關聯列同時標示 approved 與 anti-neoplastic。region span 只掃過基因不算命中，藥物關聯也不表示此突變敏感、有效或可治療。</p></article>
  </div>
  <details>
    <summary>機器證據與原始 JSON（預設收合）</summary>
    <div class="details-body">{provenance_rows}{hash_rows}
      <span class="path"><b>cohort receipt SHA-256</b> <code>{authority['cohort_receipt_sha256']}</code></span>
      <span class="path"><b>census receipt SHA-256</b> <code>{authority['census_receipt_sha256']}</code></span>
      <span class="path"><b>gene/drug annotation SHA-256</b> <code>{authority['gene_drug_sha256']}</code></span>
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
  const state = {mode:"resolution", selected:new Set(), annotationFilter:"all", query:null, current:null, filtered:rows, candidateIndex:0, shapeIndex:0};
  const modes = [
    ["resolution","Determinacy"],["coarse","拓撲形態"],["hp","HP family"],
    ["status","判定狀態"],["k","active k"],["annotation","癌症基因／藥物"]
  ];
  const annotationFilters = [
    ["all","全部 regions"],
    ["cgc","只看癌症基因 locus"],
    ["cgc_drug","只看癌症基因 ∩ approved 抗腫瘤藥"]
  ];
  const modeButtons = document.getElementById("modeButtons");
  const annotationFilterButtons = document.getElementById("annotationFilters");
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
    if (state.mode === "annotation") return row.annotationClass;
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
  function matchesAnnotation(row){
    if(state.annotationFilter==="cgc")return Boolean(row.hasCgc);
    if(state.annotationFilter==="cgc_drug")return Boolean(row.hasCgcDrug);
    return true;
  }
  function annotationFilterLabel(){
    return annotationFilters.find(item=>item[0]===state.annotationFilter)?.[1]||"全部 regions";
  }
  function refresh() {
    state.filtered = rows.filter(row => matchesQuery(row) && matchesAnnotation(row) && (!state.selected.size || state.selected.has(modeValue(row))));
    renderCanvas(); renderList();
    status.textContent = `顯示 ${state.filtered.length.toLocaleString()} / ${rows.length.toLocaleString()} regions` +
      (state.selected.size ? ` · ${state.selected.size} 類聯集` : " · 全部類別") +
      ` · ${annotationFilterLabel()}` +
      (state.query ? ` · ${state.query.chrom}:${state.query.start.toLocaleString()}–${state.query.end.toLocaleString()}` : "");
  }
  function renderModes(){
    modeButtons.innerHTML="";
    for(const [key,label] of modes){
      const button=document.createElement("button");
      button.type="button";button.className="mode-btn";button.textContent=label;
      button.dataset.mode=key;
      button.setAttribute("aria-pressed",String(state.mode===key));
      button.addEventListener("click",()=>{if(state.mode===key)return;state.mode=key;state.selected.clear();renderModes();renderLegend();refresh();});
      modeButtons.append(button);
    }
  }
  function renderAnnotationFilters(){
    annotationFilterButtons.innerHTML="";
    for(const [key,label] of annotationFilters){
      const button=document.createElement("button");
      button.type="button";button.className="annotation-filter";button.textContent=label;
      button.dataset.annotationFilter=key;
      button.setAttribute("aria-pressed",String(state.annotationFilter===key));
      button.addEventListener("click",()=>{
        if(state.annotationFilter===key)return;
        state.annotationFilter=key;renderAnnotationFilters();renderLegend();refresh();
      });
      annotationFilterButtons.append(button);
    }
  }
  function renderLegend(){
    legend.innerHTML="";
    const counts=new Map();
    rows.filter(row=>matchesQuery(row)&&matchesAnnotation(row)).forEach(row=>counts.set(modeValue(row),(counts.get(modeValue(row))||0)+1));
    for(const key of [...state.selected]){
      if(!counts.get(key))state.selected.delete(key);
    }
    for(const [key,spec] of Object.entries(defs())){
      if(!counts.get(key)) continue;
      const button=document.createElement("button");button.type="button";button.className="legend-btn";
      button.dataset.legendKey=key;
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
      const annotation=row.hasCgcDrug?'<span class="tag" style="color:#a94336">CGC ∩ drug</span>':row.hasCgc?'<span class="tag" style="color:#6b5592">CGC locus</span>':"";
      const genes=(row.geneDrug||[]).map(hit=>hit.symbol).filter(Boolean);
      button.innerHTML=`<strong>${row.chrom}:${row.start.toLocaleString()}–${row.end.toLocaleString()}</strong><span><span class="tag" style="color:${spec.color}">${spec.label}</span>${annotation}</span><small>PS ${row.ps} · HP${row.hp} · n=${row.n} · active k=${row.activeK}</small><small>${genes.length?escapeHtml([...new Set(genes)].join(" · ")):row.span.toLocaleString()+" bp"}</small>`;
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
  function coveragePct(item){
    if(item?.fraction)return fractionPct(item.fraction);
    const ref=Number(item?.ref_reads||0),alt=Number(item?.alt_reads||0);
    return ref+alt?`${(100*alt/(ref+alt)).toFixed(1)}%`:"N/A";
  }
  function locusStrip(row){
    const active=new Set((row.activePositions||[]).map(Number));
    const coverage=new Map((row.coverage||[]).map(item=>[Number(item.position),item]));
    const annotations=new Map();
    for(const hit of row.geneDrug||[]){
      for(const position of hit.loci||[]){
        const key=Number(position);if(!annotations.has(key))annotations.set(key,[]);
        annotations.get(key).push(hit);
      }
    }
    return `<div class="locus-strip" role="list" aria-label="依 GRCh38 座標排序的 locus 與 HP-specific read ALT">
      ${row.positions.map((position,index)=>{
        const item=coverage.get(Number(position))||{};
        const value=coveragePct(item);
        const numeric=value==="N/A"?0:Number(value.replace("%",""));
        const geneHits=annotations.get(Number(position))||[];
        const hasDrug=geneHits.some(hit=>(hit.drugs||[]).length);
        const geneTags=geneHits.map(hit=>`<span class="locus-gene ${hasDrug&&(hit.drugs||[]).length?"drug":""}">${escapeHtml(hit.symbol)}</span>`).join("");
        return `<article class="locus-card ${active.has(Number(position))?"active":""} ${geneHits.length?"cancer-locus":""} ${hasDrug?"drug-locus":""}" role="listitem">
          <div class="site"><b>S${index+1}</b><span>${active.has(Number(position))?"active ALT":"非 active"}</span></div>
          <span class="position">${Number(position).toLocaleString()}</span>
          <span class="af">${escapeHtml(value)}</span>
          <span class="reads">ALT ${Number(item.alt_reads||0).toLocaleString()} / REF ${Number(item.ref_reads||0).toLocaleString()}</span>
          <div class="af-track" aria-hidden="true"><span style="width:${Math.max(0,Math.min(100,numeric))}%"></span></div>
          ${geneTags?`<div class="locus-genes" aria-label="此 locus 命中的癌症基因">${geneTags}</div>`:""}
        </article>`;
      }).join("")}
    </div>`;
  }
  function drugClaimName(drug){
    return String(drug?.claimName||"");
  }
  function drugSources(drug){
    return (drug?.sources||[]).map(source=>[source?.name,source?.version].filter(Boolean).join(" ")).filter(Boolean);
  }
  function drugCell(hit){
    const drugs=hit.drugs||[];
    if(!drugs.length)return '<span class="denom">此 CGC gene 在本地 DGIdb snapshot 無 approved=true ∧ anti_neoplastic=true 關聯。</span>';
    const names=[...new Set(drugs.map(drugClaimName).filter(Boolean))].sort();
    const sources=[...new Set(drugs.flatMap(drugSources))].sort();
    const first=names.slice(0,12).map(escapeHtml).join(" · ");
    const all=names.map(escapeHtml).join(" · ");
    return `<div class="drug-list">${first}${names.length>12?`<details><summary>查看全部 ${names.length} 個 claim names</summary><div class="details-body">${all}</div></details>`:""}${sources.length?`<span class="source-list">DGIdb sources/version：${sources.map(escapeHtml).join(" · ")}</span>`:""}</div>`;
  }
  function geneDrugPanel(row){
    const hits=row.geneDrug||[];
    if(!hits.length)return `<details class="gene-drug-panel" data-testid="gene-drug-panel"><summary><span class="tag">NO_HIT_EVALUATED</span>癌症基因／藥物交集：已評估，無 actual sSNV locus 命中 CGC gene body</summary><div class="details-body"><p class="denom">NO_HIT_EVALUATED：region span 掃過基因不算命中；必須有本 group 的實際 sSNV 座標落在 GENCODE v46 gene body。</p></div></details>`;
    const drugHits=hits.filter(hit=>(hit.drugs||[]).length);
    const rows=hits.map(hit=>`<tr class="${(hit.drugs||[]).length?"has-drug":""}">
      <td>${escapeHtml(hit.symbol)}<span class="source-list">${escapeHtml(hit.geneId||"")} · ${escapeHtml(hit.hgncId||"")}</span></td>
      <td class="mono">${(hit.loci||[]).map(position=>Number(position).toLocaleString()).join(" · ")}</td>
      <td>CGC T${escapeHtml(hit.tier||"—")} · ${escapeHtml(hit.role||"未註明")}<span class="source-list">${escapeHtml(hit.tumour||"")}</span></td>
      <td>${drugCell(hit)}</td>
    </tr>`).join("");
    return `<details class="gene-drug-panel" data-testid="gene-drug-panel" open>
      <summary>癌症基因 locus 命中 ${hits.length} gene${drugHits.length?` · 同基因有 approved 抗腫瘤藥 ${drugHits.length}`:""}</summary>
      <div class="details-body">
        <p class="claim-ceiling"><b>判讀上限：</b>actual sSNV locus ∈ GENCODE gene body，再以同一 HGNC gene 連到 COSMIC CGC 與 DGIdb。DGIdb 的 approved ∧ anti-neoplastic 關聯不是 mutation-specific sensitivity、不是本癌別適應症，也不是臨床 actionability；不參與 topology candidate 排序。</p>
        <div class="table-wrap" role="region" tabindex="0" aria-label="癌症基因與 approved 抗腫瘤藥交集表"><table class="mini-table gene-drug-table"><thead><tr><th>CGC gene</th><th>命中 sSNV loci</th><th>CGC context</th><th>DGIdb approved ∩ anti-neoplastic</th></tr></thead><tbody>${rows}</tbody></table></div>
      </div>
    </details>`;
  }
  function patternEntries(row){
    const entries=[];
    for(const [pattern,count] of Object.entries(row.patterns||{}))entries.push({source:"full",pattern,count:Number(count)});
    for(const [pattern,count] of Object.entries(row.subreads||{}))entries.push({source:"partial",pattern,count:Number(count)});
    return entries.sort((a,b)=>b.count-a.count||a.source.localeCompare(b.source)||a.pattern.localeCompare(b.pattern));
  }
  function readStateMatrix(row){
    const entries=patternEntries(row);
    if(!entries.length)return '<div class="empty">此 group 沒有通過 threshold 的 solver-visible pattern。</div>';
    const header=row.positions.map((position,index)=>`<th title="${Number(position).toLocaleString()}">S${index+1}</th>`).join("");
    const body=entries.map(item=>`<tr><td>${item.source}</td>${[...item.pattern].map(value=>`<td class="state-cell ${escapeHtml(value)}">${escapeHtml(value)}</td>`).join("")}<td class="num">${item.count.toLocaleString()}</td></tr>`).join("");
    return `<div class="matrix-wrap"><table class="state-matrix"><thead><tr><th>來源</th>${header}<th>reads</th></tr></thead><tbody>${body}</tbody></table></div>`;
  }
  function popcount(value){
    let count=0,mask=Number(value);while(mask){count+=mask&1;mask>>>=1;}return count;
  }
  function patternForMask(mask,row){
    const values=Array(row.n).fill("R");
    (row.activeIndices||[]).forEach((originalIndex,bit)=>{if(Number(mask)&(1<<bit))values[Number(originalIndex)]="A";});
    return values.join("");
  }
  function observedMasks(row){
    if(row.candidate&&Array.isArray(row.candidate.observedVertices)){
      return new Set(row.candidate.observedVertices.map(Number));
    }
    const result=new Set([0]);
    for(const pattern of Object.keys(row.patterns||{})){
      let mask=0;(row.activeIndices||[]).forEach((originalIndex,bit)=>{if(pattern[Number(originalIndex)]==="A")mask|=1<<bit;});
      result.add(mask);
    }
    return result;
  }
  function vertexLabel(mask,row,compact=false){
    if(Number(mask)===0)return "ROOT";
    const observed=observedMasks(row).has(Number(mask));
    if(compact){
      const sites=[];(row.activeIndices||[]).forEach((originalIndex,bit)=>{if(Number(mask)&(1<<bit))sites.push(`S${Number(originalIndex)+1}`);});
      const state=sites.length<=3?sites.join("+"):`${sites.slice(0,2).join("+")}+…(${sites.length})`;
      return `${observed?"":"H:"}${state}`;
    }
    const pattern=patternForMask(Number(mask),row);
    return observed?pattern:`H_${pattern}`;
  }
  function acquiredPosition(parent,child,row){
    const delta=Number(parent)^Number(child);let bit=0,value=delta;while(value>1){value>>>=1;bit++;}
    const originalIndex=Number((row.activeIndices||[])[bit]);return Number(row.positions[originalIndex]);
  }
  function factoredTreeRow(row,set){
    const vertices=(set?.[0]||[]).map(mask=>({
      vertex:Number(mask),
      label:vertexLabel(mask,row),
      display_label:vertexLabel(mask,row,true)
    }));
    const edges=(set?.[5]||[]).map(item=>{
      const child=Number(item[0]),parent=Number((item[2]||item[1]||[])[0]);
      return {parent_vertex:parent,child_vertex:child,parent_label:vertexLabel(parent,row),child_label:vertexLabel(child,row),acquired_position:acquiredPosition(parent,child,row)};
    });
    return {vertices,edges};
  }
  function signatureTreeRow(row,signature){
    return {
      vertices:(signature?.[3]||[]).map(mask=>({
        vertex:Number(mask),
        label:vertexLabel(mask,row),
        display_label:vertexLabel(mask,row,true)
      })),
      edges:(signature?.[4]||[]).map(edge=>({parent_vertex:Number(edge[0]),child_vertex:Number(edge[1]),parent_label:vertexLabel(edge[0],row),child_label:vertexLabel(edge[1],row),acquired_position:acquiredPosition(edge[0],edge[1],row)}))
    };
  }
  function candidateUnionSvg(row,candidate){
    const edges=candidate.allEdges||[];
    const globalBest=new Map((candidate.bestEdges||[]).map(item=>[`${item[0]}>${item[1]}`,String(item[2])]));
    const selectedExemplar=new Set((row.edges||[]).map(item=>`${Number(item.parent_vertex)}>${Number(item.child_vertex)}`));
    const vertices=new Set([0]);edges.forEach(item=>{vertices.add(Number(item[0]));vertices.add(Number(item[1]));});
    if(vertices.size>96||edges.length>220)return `<div class="empty">完整聯集含 ${vertices.size.toLocaleString()} nodes / ${edges.length.toLocaleString()} edges，超過可讀網路上限；下方 edge census 仍完整列出所有邊，右側保留 global-best shape exemplar。</div>`;
    const levels=new Map();[...vertices].sort((a,b)=>a-b).forEach(mask=>{const level=popcount(mask);if(!levels.has(level))levels.set(level,[]);levels.get(level).push(mask);});
    const maxLevel=Math.max(1,...levels.keys()),maxLane=Math.max(1,...[...levels.values()].map(items=>items.length));
    const width=Math.max(720,180*(maxLevel+1)),height=Math.max(340,64*maxLane+70),positions=new Map();
    for(const [level,items] of levels){items.forEach((mask,index)=>positions.set(mask,{x:62+level*((width-124)/maxLevel),y:height/2+(index-(items.length-1)/2)*58}));}
    const edgeMarkup=edges.map(item=>{
      const parent=Number(item[0]),child=Number(item[1]),count=String(item[2]),key=`${parent}>${child}`,inGlobalBest=globalBest.has(key),selected=selectedExemplar.has(key),forced=count===String(candidate.totalTreeCount),a=positions.get(parent),b=positions.get(child);
      const dash=forced?"":' stroke-dasharray="5 4"';
      return `<g class="candidate-edge ${forced?"edge-forced":"edge-variable"} ${inGlobalBest?"edge-global-best":""} ${selected?"edge-selected":""}" data-parent="${parent}" data-child="${child}" data-candidate-count="${escapeHtml(count)}" data-candidate-total="${escapeHtml(candidate.totalTreeCount)}" data-best-count="${escapeHtml(globalBest.get(key)||"0")}" data-best-total="${escapeHtml(candidate.bestTreeCount)}" data-selected-exemplar="${selected?"true":"false"}"><line x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" stroke="#9b9d97" stroke-width="1.6"${dash} marker-end="url(#candidateArrow)"/>${inGlobalBest?`<line x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" stroke="#285f8f" stroke-width="2.7" opacity=".58" marker-end="url(#candidateArrowBest)"/>`:""}${selected?`<line x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" stroke="#176b58" stroke-width="4.2" opacity=".9" marker-end="url(#candidateArrowSelected)"/>`:""}<title>${escapeHtml(vertexLabel(parent,row))} → ${escapeHtml(vertexLabel(child,row))}; minimum trees ${escapeHtml(count)}/${escapeHtml(candidate.totalTreeCount)}; global-best union ${escapeHtml(globalBest.get(key)||"0")}/${escapeHtml(candidate.bestTreeCount)}; selected deterministic exemplar ${selected?"yes":"no"}</title></g>`;
    }).join("");
    const observed=observedMasks(row);
    const nodeMarkup=[...vertices].sort((a,b)=>a-b).map(mask=>{const pos=positions.get(mask),isObserved=mask===0||observed.has(mask),label=vertexLabel(mask,row,true);return `<g class="candidate-node ${isObserved?"observed":"hidden"}" data-vertex="${mask}"><circle cx="${pos.x}" cy="${pos.y}" r="22" fill="${isObserved?"#fffdf7":"#f5ead5"}" stroke="${isObserved?"#176b58":"#b66e20"}" stroke-width="2" ${isObserved?"":'stroke-dasharray="4 3"'}/><text x="${pos.x}" y="${pos.y+4}" text-anchor="middle" font-size="10" fill="#17231e">${escapeHtml(label)}</text><title>${escapeHtml(vertexLabel(mask,row))} · ${isObserved?"full-read observed":"solver-hidden/intermediate"}</title></g>`;}).join("");
    return `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="完整 minimum candidate edge union、global-best edge union 與 selected deterministic exemplar"><defs><marker id="candidateArrow" markerWidth="7" markerHeight="7" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L7,3 z" fill="#9b9d97"/></marker><marker id="candidateArrowBest" markerWidth="7" markerHeight="7" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L7,3 z" fill="#285f8f"/></marker><marker id="candidateArrowSelected" markerWidth="7" markerHeight="7" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L7,3 z" fill="#176b58"/></marker></defs>${edgeMarkup}${nodeMarkup}</svg>`;
  }
  function candidateEdgeTable(row,candidate){
    const globalBest=new Map((candidate.bestEdges||[]).map(item=>[`${item[0]}>${item[1]}`,String(item[2])]));
    const selectedExemplar=new Set((row.edges||[]).map(item=>`${Number(item.parent_vertex)}>${Number(item.child_vertex)}`));
    return (candidate.allEdges||[]).map(item=>{
      const parent=Number(item[0]),child=Number(item[1]),count=String(item[2]),key=`${parent}>${child}`,top=globalBest.get(key)||"0",selected=selectedExemplar.has(key),forced=count===String(candidate.totalTreeCount);
      return `<tr class="${top!=="0"?"edge-global-best":""} ${selected?"edge-selected":""}" data-edge-key="${key}" data-candidate-count="${escapeHtml(count)}" data-candidate-total="${escapeHtml(candidate.totalTreeCount)}" data-top-count="${escapeHtml(top)}" data-top-total="${escapeHtml(candidate.bestTreeCount)}" data-selected-exemplar="${selected?"true":"false"}"><td class="edge-key">${escapeHtml(vertexLabel(parent,row))} → ${escapeHtml(vertexLabel(child,row))}</td><td>${forced?"所有 minimum trees":"部分 minimum trees"}${top!=="0"?" · global-best union":""}${selected?" · selected exemplar":""}</td><td class="num">${escapeHtml(count)} / ${escapeHtml(candidate.totalTreeCount)}</td><td class="num">${escapeHtml(top)} / ${escapeHtml(candidate.bestTreeCount)}</td></tr>`;
    }).join("");
  }
  function candidateExplorer(row){
    if(!row.candidate){
      return `<section class="detail-section"><h3><span class="section-number">3</span>Minimum candidate family</h3><div class="empty">${row.readStatus==="ranked_complete"?"候選 factorization sidecar 尚未連接；不得從 representative tree 捏造完整候選聯集。":"此 group 未完成 AF ranking，因此沒有可宣稱的完整 ranked candidate explorer。"}</div></section>`;
    }
    const candidate=row.candidate,setCount=Number(candidate.minimumSetCount),index=Math.max(0,Math.min(setCount-1,Number(state.candidateIndex)||0));state.candidateIndex=index;
    const set=candidate.sets[index],isGlobal=Number(set[4])===1,shapeCount=(candidate.signatures||[]).length,shapeIndex=Math.max(0,Math.min(shapeCount-1,Number(state.shapeIndex)||0));state.shapeIndex=shapeIndex;
    const shape=candidate.signatures[shapeIndex],candidateTree=factoredTreeRow(row,set),shapeTree=signatureTreeRow(row,shape);
    const shapeOptions=(candidate.signatures||[]).map((item,itemIndex)=>`<option value="${itemIndex}" ${itemIndex===shapeIndex?"selected":""}>shape ${itemIndex+1} · ${escapeHtml(item[1])} · ${escapeHtml(item[2])} trees</option>`).join("");
    return `<details id="candidateExplorer" class="candidate-space" ${setCount>1||row.resolution==="TIED_CROSS_TOPOLOGY"?"open":""}>
      <summary><span class="section-number">3</span>完整 minimum candidate space 與 AF 排序</summary>
      <div class="details-body">
        <div class="callout"><b>如何讀：</b>灰線是所有 minimum structural trees 的完整 edge union；虛線只出現在部分候選；藍線是所有 AF-global-best trees 的 edge union；綠線只疊目前頁面採用的一棵 deterministic selected exemplar。任何聯集都不是一棵樹，edge 次數是 candidate membership，不是 probability 或 read support。</div>
        <div class="candidate-summary">
          <div><span>minimum vertex sets</span><b>${setCount.toLocaleString()}</b></div>
          <div><span>all minimum trees</span><b>${escapeHtml(candidate.totalTreeCount)}</b></div>
          <div><span>global-best sets</span><b>${Number(candidate.bestSetCount).toLocaleString()}</b></div>
          <div><span>global-best trees</span><b>${escapeHtml(candidate.bestTreeCount)}</b></div>
          <div><span>exact best shapes</span><b>${shapeCount.toLocaleString()}</b></div>
        </div>
        <div class="candidate-grid">
          <div><h3>完整 edge union＋global-best union＋selected tree</h3><div class="candidate-network" data-testid="candidate-union">${candidateUnionSvg(row,candidate)}</div><p class="denom">所有 ${escapeHtml(candidate.totalTreeCount)} minimum trees 均已計入；藍色 global-best union 可含多棵同分樹，綠色 selected exemplar 才是一棵樹，且其 edges 必須是兩層 union 的子集。</p></div>
          <div><h3>Global-best shape exemplar</h3><label class="denom" for="shapeSelect">shape census（exemplar 不是唯一 tree）</label><select id="shapeSelect" class="action" style="width:100%;margin:.35rem 0">${shapeOptions}</select><div class="tree">${treeSvg(shapeTree)}</div><p class="denom">${escapeHtml(shape?.[0]||"—")} · ${escapeHtml(shape?.[1]||"—")} · ${escapeHtml(shape?.[2]||"0")} global-best trees</p></div>
        </div>
        <details><summary>完整候選 edge census（所有合理邊）</summary><div class="details-body"><div class="table-wrap"><table class="mini-table candidate-edge-table"><thead><tr><th>edge</th><th>角色</th><th>minimum tree incidence</th><th>global-best incidence</th></tr></thead><tbody>${candidateEdgeTable(row,candidate)}</tbody></table></div></div></details>
        <h3 style="margin-top:1rem">逐 minimum vertex set 檢視</h3>
        <div class="candidate-toolbar"><button type="button" class="action" data-candidate-step="-1" ${index===0?"disabled":""}>← 前一組</button><label>候選 <input id="candidateIndex" type="number" min="1" max="${setCount}" value="${index+1}"></label><button type="button" class="action" id="candidateGo">前往</button><button type="button" class="action" data-candidate-step="1" ${index===setCount-1?"disabled":""}>下一組 →</button><span class="tag" style="color:${isGlobal?"#176b58":"#667069"}">${isGlobal?"AF global-best":"非 global-best"}</span></div>
        <div class="candidate-grid"><div><div class="tree" data-testid="candidate-tree">${treeSvg(candidateTree)}</div></div><div><div class="evidence-cell"><span>vertex set</span><b>${(set[0]||[]).map(mask=>escapeHtml(vertexLabel(mask,row))).join(" · ")}</b></div><div class="evidence-cell"><span>AF-best score within set</span><b>${escapeHtml(set[1])}</b></div><div class="evidence-cell"><span>legal / within-set best trees</span><b>${escapeHtml(set[2])} / ${escapeHtml(set[3])}</b></div><p class="denom">畫面選每個 child 的第一個 exact best parent；若 within-set best trees&gt;1，這仍只是一棵 exemplar。</p></div></div>
      </div>
    </details>`;
  }
  function bindCandidateControls(row){
    if(!row.candidate)return;
    const refresh=()=>{const current=document.getElementById("candidateExplorer");if(!current)return;const open=current.open;current.outerHTML=candidateExplorer(row);const replacement=document.getElementById("candidateExplorer");if(replacement)replacement.open=open;bindCandidateControls(row);};
    document.querySelectorAll("#candidateExplorer [data-candidate-step]").forEach(button=>button.addEventListener("click",()=>{state.candidateIndex=Math.max(0,Math.min(row.candidate.sets.length-1,state.candidateIndex+Number(button.dataset.candidateStep)));refresh();}));
    document.getElementById("candidateGo")?.addEventListener("click",()=>{const input=document.getElementById("candidateIndex");state.candidateIndex=Math.max(0,Math.min(row.candidate.sets.length-1,Number(input?.value||1)-1));refresh();});
    document.getElementById("candidateIndex")?.addEventListener("keydown",event=>{if(event.key==="Enter"){event.preventDefault();document.getElementById("candidateGo")?.click();}});
    document.getElementById("shapeSelect")?.addEventListener("change",event=>{state.shapeIndex=Number(event.target.value||0);refresh();});
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
    const nodes=[...vertices].map(([id,v])=>{
      const p=pos.get(id),hidden=id!=="0"&&String(v.label||"").startsWith("H_");
      const visibleLabel=v.display_label||(Array.isArray(row.activeIndices)?vertexLabel(Number(id),row,true):v.label);
      return `<g class="${hidden?"hidden":"observed"}"><circle cx="${p.x}" cy="${p.y}" r="19" fill="${id==="0"?"#17231e":hidden?"#f5ead5":"#fffdf7"}" stroke="${hidden?"#b66e20":"#176b58"}" stroke-width="2"${hidden?' stroke-dasharray="4 3"':""}/><text x="${p.x}" y="${p.y+4}" text-anchor="middle" font-size="10" fill="${id==="0"?"#fff":"#17231e"}">${escapeHtml(visibleLabel)}</text><title>${escapeHtml(v.label)} · ${hidden?"solver-required hidden ancestor":"full-read observed/root"}</title></g>`;
    }).join("");
    return `<svg viewBox="0 0 520 290" role="img" aria-label="AF-global-best representative tree；虛線節點為 solver-required hidden ancestor"><defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#52625a"/></marker></defs>${edges}${nodes}</svg>`;
  }
  function patternRows(obj,label){
    const entries=Object.entries(obj||{}).sort((a,b)=>b[1]-a[1]);
    if(!entries.length)return `<tr><td>${label}</td><td colspan="2">無</td></tr>`;
    return entries.slice(0,24).map(([p,n])=>`<tr><td>${label}</td><td class="mono">${escapeHtml(p)}</td><td class="num">${Number(n).toLocaleString()}</td></tr>`).join("");
  }
  function renderDetail(row){
    const rSpec=data.definitions.resolution[row.resolution]||{label:row.resolution,color:"#444"};
    const cSpec=data.definitions.coarse[row.coarse]||{label:row.coarse,color:"#444"};
    const aSpec=data.definitions.annotation[row.annotationClass]||{label:row.annotationClass,color:"#444"};
    const coverage=(row.coverage||[]).map(c=>`<tr><td class="mono">${Number(c.position).toLocaleString()}</td><td class="num">${Number(c.ref_reads||0).toLocaleString()}</td><td class="num">${Number(c.alt_reads||0).toLocaleString()}</td><td class="num">${c.fraction?escapeHtml(c.fraction):"—"}</td><td class="num">${c.fraction?fractionPct(c.fraction):((Number(c.ref_reads||0)+Number(c.alt_reads||0))?`${(100*Number(c.alt_reads||0)/(Number(c.ref_reads||0)+Number(c.alt_reads||0))).toFixed(1)}%`:"N/A")}</td></tr>`).join("");
    const signatures=(row.signatures||[]).map(s=>`<tr><td class="mono">${escapeHtml(s.shape)}</td><td>${escapeHtml(s.coarse)}</td><td class="num">${escapeHtml(s.trees)}</td></tr>`).join("")||'<tr><td colspan="3">只對 ranked complete groups 提供 exact signature census。</td></tr>';
    const boundary=row.boundaryContext;
    const boundaryPanel=boundary?`<div class="callout"><b>查詢窗拆分證據：</b>
      S1=${Number(boundary.includedPositions[0]).toLocaleString()} 與 S2=${Number(boundary.includedPositions[1]).toLocaleString()} 有 ${Number(boundary.endpointSupport.total).toLocaleString()} 條 fixed-endpoint molecules（RR=${Number(boundary.endpointSupport.RR)}, AA=${Number(boundary.endpointSupport.AA)}），因此在 threshold=${Number(boundary.threshold)} 形成同一 strict component。
      S3=${Number(boundary.excludedSingletons[0]).toLocaleString()}、S4=${Number(boundary.excludedSingletons[1]).toLocaleString()} 沒有達門檻的 incident edge，保留為 singleton ABSTAIN，不進本樹。
      <span class="denom">AA=2 低於 MLHP MINREAD=3，故可支持物理連通，卻不作 solver-visible full genotype population；同 PS 本身仍不足以合併。</span></div>`
      :`<div class="callout"><b>邊界證據：</b>同一 PS=${escapeHtml(row.ps)}、同一 HP${row.hp}、threshold=3 strict endpoint component；unit <code>${escapeHtml(row.unit)}</code> / block <code>${escapeHtml(row.block)}</code>。此 group 不由 50 kb gap 規則建立。</div>`;
    detail.innerHTML=`
      <div class="detail-head"><div><div class="eyebrow">Exact-PS regional evidence</div><h3 class="detail-id">${escapeHtml(row.id)}</h3></div><div><span class="tag" style="color:${rSpec.color}">${rSpec.label}</span><span class="tag" style="color:${cSpec.color}">${cSpec.label}</span><span class="tag" style="color:${aSpec.color}">${aSpec.label}</span></div></div>
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
      ${boundaryPanel}
      ${geneDrugPanel(row)}
      <section class="detail-section"><h3><span class="section-number">1</span>Locus 排列與 read-ALT</h3><p class="denom">依 GRCh38 座標排列；綠底是進入 compact topology 的 active loci。百分比為同一 PS × HP family 的 ALT/(REF+ALT)，不是 caller VAF。</p>${locusStrip(row)}</section>
      <details open><summary><span class="section-number">2</span>Solver-visible read state matrix</summary><div class="details-body"><p>full coverage=${row.fullReads.toLocaleString()} · HP reads=${row.hpReads.toLocaleString()} · projected rows=${row.projectedRows.toLocaleString()} · tree-supported molecule incidences=${row.treeMolecules.toLocaleString()}</p>${readStateMatrix(row)}<p class="denom">只列 count≥3 的 solver-visible full/partial patterns；低於 threshold 的 molecule pattern 不在 MLHP group payload，不能從本表推成「不存在」。</p></div></details>
      ${candidateExplorer(row)}
      <section class="detail-section"><h3><span class="section-number">4</span>AF-global-best 結果與 determinacy</h3>
        <div class="detail-grid">
          <div><h3>${row.treeUnique?"唯一 selected tree":"一棵 selected exemplar"}</h3><div class="tree">${treeSvg(row)}</div><p class="denom">${row.treeUnique?"此樹是目前模型下唯一 global-best tree；仍不是 biological truth。":"這只是 deterministic exemplar；所有同分結果需由候選集合與 signature census 判讀。"}</p></div>
          <div><h3>Exact topology signatures</h3><div class="table-wrap"><table class="mini-table"><thead><tr><th>rooted shape</th><th>coarse</th><th>best trees</th></tr></thead><tbody>${signatures}</tbody></table></div></div>
        </div>
      </section>
      <details open><summary>排序依據｜每個 sSNV 的 HP-specific read-ALT</summary><div class="details-body"><div class="table-wrap"><table class="mini-table"><thead><tr><th>position</th><th>REF</th><th>ALT</th><th>exact fraction</th><th>ALT/(R+A)</th></tr></thead><tbody>${coverage}</tbody></table></div><p class="denom">${escapeHtml(row.coverageMeaning)}</p></div></details>
      <details><summary>Raw pattern counts（文字表）</summary><div class="details-body"><div class="table-wrap"><table class="mini-table"><thead><tr><th>來源</th><th>pattern</th><th>reads</th></tr></thead><tbody>${patternRows(row.patterns,"full")}${patternRows(row.subreads,"partial")}</tbody></table></div></div></details>
      <details><summary>識別鍵與 solver 狀態</summary><div class="details-body"><span class="path"><b>component</b> <code>${escapeHtml(row.component)}</code></span><span class="path"><b>unit</b> <code>${escapeHtml(row.unit)}</code></span><span class="path"><b>block</b> <code>${escapeHtml(row.block)}</code></span><span class="path"><b>family/objective</b> <code>${escapeHtml(row.familyStatus)} / ${escapeHtml(row.objectiveStatus)}</code></span><span class="path"><b>search nodes</b> <code>${row.searchNodes.toLocaleString()}</code></span><p>${escapeHtml(row.message)}</p></div></details>`;
    bindCandidateControls(row);
  }
  function selectRow(row,scroll){
    state.current=row;
    state.candidateIndex=row.candidate?Number(row.candidate.selectedSetIndex||0):0;
    state.shapeIndex=0;
    renderDetail(row);renderCanvas();renderList();
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
  renderModes();renderAnnotationFilters();renderLegend();refresh();
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
    similarity = authority["similarity"]
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
    samples_by_name = {item["sample"]: item for item in samples}
    require(
        set(samples_by_name) == set(DISPLAY_ORDER),
        "cohort display-order sample set mismatch",
    )
    sample_rows = []
    for sample in DISPLAY_ORDER:
        item = samples_by_name[sample]
        sample = item["sample"]
        t = topology_map[sample]
        c = census_map[sample]
        ranked = int(c["ranked_units"])
        unique = int(c["resolution"]["UNIQUE_TREE"]["n"])
        exact = int(c["one_exact_topology"]["n"])
        abstain = int(t["mutation_family_abstain_units"])
        cgc_regions = item["summary"].get("cgcRegions")
        cgc_drug_regions = item["summary"].get("cgcDrugRegions")
        cgc_text = f"{int(cgc_regions):,}" if cgc_regions is not None else "—"
        cgc_drug_text = (
            f"{int(cgc_drug_regions):,}" if cgc_drug_regions is not None else "—"
        )
        sample_rows.append(
            f"""<tr>
<td><a class="sample-link" href="{esc(sample)}.html" title="authority key: {esc(sample)}">{esc(DISPLAY_LABELS[sample])}</a>{'<br><span class="tag">same biological sample</span>' if sample in {'HCC1395','HCC1395_DORADO'} else ''}</td>
<td class="num">{int(t['groups_total']):,}</td><td class="num">{ranked:,}</td>
<td><b>{pct(unique,ranked):.1f}%</b><div class="barline"><span style="width:{pct(unique,ranked):.4f}%"></span></div></td>
<td class="num">{pct(exact,ranked):.1f}%</td><td class="num danger">{abstain:,}</td><td class="num">{cgc_text}</td><td class="num">{cgc_drug_text}</td>
<td><a href="{esc(sample)}.html#genome">開啟 GRCh38 →</a></td></tr>"""
        )

    res_keys = ("UNIQUE_TREE", "TIED_SAME_TOPOLOGY", "TIED_CROSS_TOPOLOGY")
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
    hcc_pair = next(
        row
        for row in similarity["pairwise_profiles"]
        if row["pair_id"] == "HCC1395__HCC1395_DORADO"
    )
    bootstrap = similarity["hcc1395_dorado_chromosome_block_bootstrap"]
    strict = similarity["hcc1395_dorado_strict_comparison"]
    ordered_overlap = strict["ordered_locus_key_overlap"]
    strict_compare = strict["unambiguous_comparison"]
    ladder_specs = (
        (
            "可嚴格 1:1 配對",
            int(strict_compare["pair_count"]),
            int(ordered_overlap["intersection"]),
            f"{int(ordered_overlap['ambiguous_common_keys_not_paired']):,} 個重複 key 不強配 HP",
        ),
        (
            "兩邊皆 ranked",
            int(strict_compare["both_ranked"]["numerator"]),
            int(strict_compare["both_ranked"]["denominator"]),
            "相同 ordered-locus key 的 technical pair",
        ),
        (
            "active loci 完全相同",
            int(strict_compare["same_active_loci_given_both_ranked"]["numerator"]),
            int(strict_compare["same_active_loci_given_both_ranked"]["denominator"]),
            "先要求兩邊皆 ranked",
        ),
        (
            "determinacy class 相同",
            int(strict_compare["same_resolution_class"]["numerator"]),
            int(strict_compare["same_resolution_class"]["denominator"]),
            "在 active loci 可比較集合",
        ),
        (
            "exact shape 集合相同",
            int(strict_compare["same_shape_signature_set"]["numerator"]),
            int(strict_compare["same_shape_signature_set"]["denominator"]),
            "rooted-unlabeled signature set",
        ),
        (
            "兩邊皆為唯一最佳樹",
            int(strict_compare["both_unique_tree"]["numerator"]),
            int(strict_compare["both_unique_tree"]["denominator"]),
            "selected-tree 比較的必要 gate",
        ),
        (
            "selected labeled tree 相同",
            int(strict_compare["same_selected_labeled_tree"]["numerator"]),
            int(strict_compare["same_selected_labeled_tree"]["denominator"]),
            "最嚴格逐區比較；不是 clone identity",
        ),
    )
    ladder_html = "".join(
        '<div class="ladder-step">'
        f'<div class="ladder-label">{esc(label)}<small>{esc(note)}</small></div>'
        f'<div class="ladder-track" aria-hidden="true"><span style="width:{pct(numerator, denominator):.6f}%"></span></div>'
        f'<div class="ladder-value">{numerator:,}/{denominator:,}<br>{pct(numerator, denominator):.1f}%</div>'
        "</div>"
        for label, numerator, denominator, note in ladder_specs
    )
    cohort_payload = {
        "sampleProfiles": similarity["sample_profiles"],
        "pairs": similarity["pairwise_profiles"],
        "sampleOrder": DISPLAY_ORDER,
        "displayLabels": DISPLAY_LABELS,
        "metricContract": similarity["metric_contract"],
    }
    cohort_payload_json = json.dumps(
        cohort_payload, ensure_ascii=False, separators=(",", ":")
    ).replace("</", "<\\/")
    evidence_paths = {
        "topology cohort receipt": authority["cohort_receipt_path"],
        "topology all7 summary": authority["all7_summary_path"],
        "signature census receipt": authority["census_receipt_path"],
        "signature census summary": authority["census_summary_path"],
        "cohort similarity sidecar": authority["similarity_path"],
        "gene/drug annotation sidecar": authority["gene_drug_path"],
        "gene/drug annotation receipt": authority["gene_drug_receipt_path"],
        "k / HP funnel": authority["funnel_path"]
        or "derived from hash-bound per-sample MLHP receipts",
    }
    path_html = "".join(
        f'<span class="path"><b>{esc(label)}</b> <code>{esc(path)}</code></span>'
        for label, path in evidence_paths.items()
    )
    cohort_js = r"""
(() => {
  const data = JSON.parse(document.getElementById("cohortData").textContent);
  const profileModes = {
    resolution: {
      dimension: "resolution",
      label: "Determinacy",
      denominator: "ranked complete",
      labels: {
        UNIQUE_TREE: "唯一最佳樹",
        TIED_SAME_TOPOLOGY: "多樹・同 topology",
        TIED_CROSS_TOPOLOGY: "多樹・跨 topology"
      },
      colors: {
        UNIQUE_TREE: "#176b58",
        TIED_SAME_TOPOLOGY: "#285f8f",
        TIED_CROSS_TOPOLOGY: "#a94336"
      }
    },
    morphology: {
      dimension: "morphology",
      label: "Topology / coarse",
      denominator: "ranked complete",
      labels: {
        "Single-only": "單支",
        "Sister-only": "旁系",
        "Direct-only": "直系",
        "Sister+direct": "直系＋旁系",
        "Unresolved-cross-coarse": "跨形態未定"
      },
      colors: {
        "Single-only": "#66717f",
        "Sister-only": "#6b5592",
        "Direct-only": "#176b58",
        "Sister+direct": "#b66e20",
        "Unresolved-cross-coarse": "#a94336"
      }
    },
    active_k: {
      dimension: "active_k",
      label: "Complexity / active k",
      denominator: "all final groups",
      labels: Object.fromEntries(Array.from({length: 13}, (_, i) => [String(i), `active k=${i}`])),
      colors: {
        "0":"#b9b8ad","1":"#6483a3","2":"#176b58","3":"#4d8c69",
        "4":"#8a8b3f","5":"#b66e20","6":"#a94336","7":"#8f4057",
        "8":"#6b5592","9":"#4f477d","10":"#354768","11":"#293b50","12":"#17231e"
      }
    }
  };
  const similarityModes = {
    composite: {label: "Composite", component: null},
    resolution: {label: "Determinacy", component: "resolution"},
    morphology: {label: "Topology / coarse", component: "morphology"},
    active_k: {label: "Complexity / active k", component: "active_k"},
    status: {label: "State", component: "status"},
    hp_symmetric: {label: "HP balance", component: "hp_symmetric"}
  };
  let profileMode = "resolution";
  let similarityMode = "composite";
  const escapeText = value => String(value).replace(/[&<>"']/g, char => ({
    "&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#039;"
  })[char]);
  const percent = value => `${(100 * Number(value)).toFixed(1)}%`;
  const sampleLabel = sample => data.displayLabels[sample] || sample;
  const pairValue = (row, mode) => mode === "composite"
    ? Number(row.profile_similarity)
    : Number(row.components[similarityModes[mode].component].similarity);
  const pairMap = new Map();
  data.pairs.forEach(row => {
    pairMap.set(`${row.sample_a}|${row.sample_b}`, row);
    pairMap.set(`${row.sample_b}|${row.sample_a}`, row);
  });

  function renderProfiles() {
    const mode = profileModes[profileMode];
    const chart = document.getElementById("profileChart");
    const first = data.sampleProfiles[data.sampleOrder[0]].dimensions[mode.dimension];
    const bins = first.bins;
    chart.innerHTML = data.sampleOrder.map(sample => {
      const dimension = data.sampleProfiles[sample].dimensions[mode.dimension];
      const dominant = bins.reduce(
        (best, bin) => Number(dimension.proportions[bin] || 0) > Number(dimension.proportions[best] || 0) ? bin : best,
        bins[0]
      );
      const aria = bins.map(bin => `${mode.labels[bin] || bin} ${percent(dimension.proportions[bin] || 0)}`).join("；");
      const segments = bins.map(bin => {
        const value = Number(dimension.proportions[bin] || 0);
        return `<span style="width:${100 * value}%;background:${mode.colors[bin]}" title="${escapeText(mode.labels[bin] || bin)}：${percent(value)}（${Number(dimension.counts[bin] || 0).toLocaleString()}）"></span>`;
      }).join("");
      const technical = sample === "HCC1395" || sample === "HCC1395_DORADO"
        ? "<small>same biological sample</small>" : "<small>independent sample</small>";
      return `<div class="profile-row">
        <a class="profile-sample" href="${escapeText(sample)}.html" title="authority key: ${escapeText(sample)}">${escapeText(sampleLabel(sample))}${technical}</a>
        <div class="profile-stack" role="img" aria-label="${escapeText(sampleLabel(sample))} ${escapeText(mode.label)}：${escapeText(aria)}">${segments}</div>
        <div class="profile-dominant"><b>${escapeText(mode.labels[dominant] || dominant)} ${percent(dimension.proportions[dominant] || 0)}</b>n=${Number(dimension.denominator).toLocaleString()} ${escapeText(mode.denominator)}</div>
      </div>`;
    }).join("");
    document.getElementById("profileLegend").innerHTML = bins.map(bin =>
      `<span><i style="background:${mode.colors[bin]}"></i>${escapeText(mode.labels[bin] || bin)}</span>`
    ).join("");
    document.getElementById("profileDenominator").textContent =
      `${mode.label}；每列獨立正規化，分母為 ${mode.denominator}。`;
    document.querySelectorAll("[data-profile-mode]").forEach(button => {
      button.setAttribute("aria-pressed", String(button.dataset.profileMode === profileMode));
    });
  }

  function heatStyle(value) {
    const strength = Math.max(0, Math.min(1, (value - .55) / .45));
    const lightness = 95 - 52 * strength;
    const saturation = 18 + 38 * strength;
    return `background:hsl(159 ${saturation}% ${lightness}%);color:${strength > .58 ? "#fff" : "#17231e"}`;
  }

  function renderSimilarity() {
    const order = data.sampleOrder;
    const label = similarityModes[similarityMode].label;
    const header = order.map(sample =>
      `<th scope="col" title="authority key: ${escapeText(sample)}">${escapeText(sampleLabel(sample))}</th>`
    ).join("");
    const body = order.map(sampleA => {
      const cells = order.map(sampleB => {
        const row = pairMap.get(`${sampleA}|${sampleB}`);
        const value = sampleA === sampleB ? 1 : pairValue(row, similarityMode);
        const technical = (
          (sampleA === "HCC1395" && sampleB === "HCC1395_DORADO") ||
          (sampleB === "HCC1395" && sampleA === "HCC1395_DORADO")
        );
        return `<td class="${technical ? "technical-pair" : ""}" style="${heatStyle(value)}" title="${escapeText(sampleLabel(sampleA))} × ${escapeText(sampleLabel(sampleB))} · ${escapeText(label)} ${value.toFixed(6)}" aria-label="${escapeText(sampleLabel(sampleA))} 與 ${escapeText(sampleLabel(sampleB))} ${escapeText(label)} 相似度 ${value.toFixed(6)}">${value.toFixed(3)}</td>`;
      }).join("");
      return `<tr><th scope="row" title="authority key: ${escapeText(sampleA)}">${escapeText(sampleLabel(sampleA))}</th>${cells}</tr>`;
    }).join("");
    document.getElementById("similarityHeatmap").innerHTML =
      `<table class="heatmap-table"><caption class="sr-only">7 組 technical datasets 的 ${escapeText(label)} 相似度矩陣</caption><thead><tr><th scope="col">${escapeText(label)}</th>${header}</tr></thead><tbody>${body}</tbody></table>`;
    const ranked = [...data.pairs].sort((left, right) =>
      pairValue(right, similarityMode) - pairValue(left, similarityMode)
    );
    document.getElementById("topPairs").innerHTML = ranked.slice(0, 5).map((row, index) => {
      const technical = row.relationship === "same_biological_sample_technical_pair";
      return `<li class="pair-rank ${technical ? "technical" : ""}">
        <span class="rank">${index + 1}</span>
        <span class="pair">${escapeText(sampleLabel(row.sample_a))} × ${escapeText(sampleLabel(row.sample_b))}<small>${technical ? "same biological sample · technical pair" : "cross-biological distribution only"}</small></span>
        <span class="score">${pairValue(row, similarityMode).toFixed(3)}</span>
      </li>`;
    }).join("");
    document.getElementById("heatmapCaption").textContent =
      `${label}：1 − Jensen–Shannon distance；對角線為同資料集自比。橘框標示 HCC1395_HKU × HCC1395_NYGC technical pair。`;
    document.getElementById("pairRankingTitle").textContent = `${label} · top 5 pairs`;
    document.querySelectorAll("[data-similarity-mode]").forEach(button => {
      button.setAttribute("aria-pressed", String(button.dataset.similarityMode === similarityMode));
    });
  }

  document.querySelectorAll("[data-profile-mode]").forEach(button => {
    button.addEventListener("click", () => {
      profileMode = button.dataset.profileMode;
      renderProfiles();
    });
  });
  document.querySelectorAll("[data-similarity-mode]").forEach(button => {
    button.addEventListener("click", () => {
      similarityMode = button.dataset.similarityMode;
      renderSimilarity();
    });
  });
  renderProfiles();
  renderSimilarity();
})();
"""
    return f"""<!doctype html>
<html lang="zh-Hant"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<meta name="intersubmod-authority" content="{AUTHORITY_NAME}">
<meta name="intersubmod-ui-contract" content="{UI_CONTRACT}">
<meta name="intersubmod-cohort-receipt-sha256" content="{authority['cohort_receipt_sha256']}">
<meta name="intersubmod-census-receipt-sha256" content="{authority['census_receipt_sha256']}">
<meta name="intersubmod-similarity-sha256" content="{authority['similarity_sha256']}">
<meta name="intersubmod-gene-drug-sha256" content="{authority['gene_drug_sha256']}">
<title>Exact-PS layered workstation · cohort command center</title>
<style>{COMMON_CSS}</style></head>
<body>
<div class="authority"><b>PRIMARY · 2026-07-24 exact PS × HP</b><span>7 datasets / 6 biological samples · chr1–22 · strict endpoint read-linkage ≥3</span></div>
<header class="hero"><div class="shell">
  <div class="eyebrow">InterSubMod evidence atlas</div>
  <h1>從 read-linkage<br>到 exact topology</h1>
  <p class="lead">新的 layered workstation 以 exact phase-set、primary HP 與嚴格 read-linkage 重建區域；舊 50 kb proximity grouping 不再是預設資料。</p>
  <div class="boundary"><strong>最重要限制</strong><span>7/7 pipeline 是 technical PASS，但全 cohort 仍有 {int(totals['mutation_family_abstain_units']):,} 個 mutation-bearing groups 因資源 guard ABSTAIN，因此不是 topology-complete。所有比例均顯示明確分母。</span></div>
  <nav class="nav"><a href="#evidence">資料漏斗</a><a href="#topology">拓撲全貌</a><a href="#compare">跨樣本比較</a><a href="#samples">7 組資料</a><a href="#pair">HCC1395 技術驗證</a><a href="#provenance">證據</a></nav>
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
    <article class="panel"><h3>Determinacy｜最佳樹同形到哪一層</h3><div class="denom">n={int(census['ranked_units']):,}</div>{resolution_panel}</article>
    <article class="panel"><h3>拓撲型態｜分子累積形狀</h3><div class="denom">n={int(census['ranked_units']):,}</div>{coarse_panel}</article>
  </div>
  <div class="callout" style="margin-top:1rem"><b>精確結論：</b>單一 rooted-unlabeled topology = {int(census['one_exact_topology']['n']):,}/{int(census['ranked_units']):,} = {float(census['one_exact_topology']['pct_ranked']):.4f}%；單一四類 coarse geometry = {int(census['one_coarse_class']['n']):,}/{int(census['ranked_units']):,} = {float(census['one_coarse_class']['pct_ranked']):.4f}%。兩者不能替代 cellular clone / ancestry 驗證。</div>
</div></section>
<section class="section" id="compare"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Cross-dataset state atlas</div><h2>七組資料的狀態、拓撲與複雜度並排比較</h2></div><p>先看各資料集的組成輪廓，再以同一維度檢查兩兩相似度。每一列獨立正規化；ranked-only 與 all-groups 分母不混用。</p></div>
  <div class="cohort-dashboard">
    <article class="cohort-panel" aria-labelledby="profile-title">
      <div class="cohort-panel-head">
        <div><h3 id="profile-title">Sample profile composition</h3><p id="profileDenominator">Determinacy；每列獨立正規化，分母為 ranked complete。</p></div>
        <div class="compare-tabs" role="group" aria-label="切換樣本組成維度">
          <button class="compare-btn" type="button" data-profile-mode="resolution" aria-pressed="true">Determinacy</button>
          <button class="compare-btn" type="button" data-profile-mode="morphology" aria-pressed="false">Topology / coarse</button>
          <button class="compare-btn" type="button" data-profile-mode="active_k" aria-pressed="false">Complexity / active k</button>
        </div>
      </div>
      <div class="profile-chart" id="profileChart" aria-live="polite"></div>
      <div class="profile-legend" id="profileLegend" aria-label="組成圖例"></div>
    </article>
    <article class="cohort-panel" aria-labelledby="similarity-title">
      <div class="cohort-panel-head">
        <div><h3 id="similarity-title">Pairwise profile similarity</h3><p>Composite 等權整合 state、determinacy、topology、active k 與對稱化 HP balance；切換維度可確認相似性由何處構成。</p></div>
        <div class="compare-tabs" role="group" aria-label="切換相似度維度">
          <button class="compare-btn" type="button" data-similarity-mode="composite" aria-pressed="true">Composite</button>
          <button class="compare-btn" type="button" data-similarity-mode="resolution" aria-pressed="false">Determinacy</button>
          <button class="compare-btn" type="button" data-similarity-mode="morphology" aria-pressed="false">Topology / coarse</button>
          <button class="compare-btn" type="button" data-similarity-mode="active_k" aria-pressed="false">Active k</button>
          <button class="compare-btn" type="button" data-similarity-mode="status" aria-pressed="false">State</button>
          <button class="compare-btn" type="button" data-similarity-mode="hp_symmetric" aria-pressed="false">HP balance</button>
        </div>
      </div>
      <div class="heatmap-layout">
        <div>
          <div class="heatmap-scroll" id="similarityHeatmap" aria-live="polite"></div>
          <p class="heatmap-caption" id="heatmapCaption"></p>
        </div>
        <aside>
          <h3 id="pairRankingTitle">Composite · top 5 pairs</h3>
          <ol class="pair-ranking" id="topPairs"></ol>
          <p class="heatmap-caption">排名只表示分布輪廓接近；cross-biological pair 不可據此視為同一 clone 或共同 ancestry。</p>
        </aside>
      </div>
    </article>
  </div>
</div></section>
<section class="section" id="samples"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Seven technical datasets</div><h2>每組資料直接進入全基因視圖</h2></div><p>HCC1395_HKU 與 HCC1395_NYGC 是同一 biological sample 的兩個 technical datasets；因此 7 組資料只代表 6 個生物樣本。括號外顯示名稱不改變底層 authority key 與檔名。</p></div>
  <div class="table-wrap"><table><thead><tr><th>dataset</th><th>final groups</th><th>ranked</th><th>唯一最佳樹</th><th>單一 exact shape</th><th>ABSTAIN</th><th>CGC locus</th><th>CGC ∩ drug</th><th>工作站</th></tr></thead><tbody>{''.join(sample_rows)}</tbody></table></div>
  <p class="denom">CGC locus 與 CGC ∩ drug 都以全部 final groups 為分母：至少一個實際 sSNV position 落在 GENCODE v46 gene body；後者另要求同一 HGNC gene 為 COSMIC CGC v104 且有 DGIdb approved=true ∧ anti_neoplastic=true 關聯。不是臨床 actionability。</p>
</div></section>
<section class="section" id="pair"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Technical concordance</div><h2>HCC1395_HKU × HCC1395_NYGC</h2></div><p>這是同一生物樣本在不同 technical datasets 下的組成比較，不能算兩個獨立 biological replicates；底層 authority keys 分別為 HCC1395 與 HCC1395_DORADO。</p></div>
  <div class="hcc-scoreboard">
    <article class="hcc-score">
      <div><span class="label">five-dimension profile similarity</span><span class="value">{float(hcc_pair['profile_similarity']):.6f}</span><span class="interval">chr-block bootstrap 95% CI<br>{float(bootstrap['similarity']['p2_5']):.6f}–{float(bootstrap['similarity']['p97_5']):.6f}</span></div>
      <div class="note">point rank #{int(bootstrap['rank_of_21']['point'])}/21 · rank #1 in {100*float(bootstrap['rank_of_21']['rank_1_rate']):.1f}% · top 3 in {100*float(bootstrap['rank_of_21']['top_3_rate']):.1f}% of {int(bootstrap['replicates']):,} chromosome-block replicates</div>
    </article>
    <article class="cohort-panel">
      <div class="cohort-panel-head"><div><h3>從共同座標到相同 selected tree</h3><p>每一階使用下方標明的條件分母；重複 ordered-locus key 不以 HP1/HP2 標號硬配。</p></div></div>
      <div class="evidence-ladder">{ladder_html}</div>
    </article>
  </div>
  <div class="claim-boundary"><strong>Profile 相似不等於同一棵樹，更不等於同一 clone。</strong>0.926381 是五個 cohort-level composition 維度的分布相似度；嚴格逐區比較只在兩邊皆為唯一最佳樹時成立，selected labeled tree 相同為 {int(strict_compare['same_selected_labeled_tree']['numerator']):,}/{int(strict_compare['same_selected_labeled_tree']['denominator']):,}（{100*float(strict_compare['same_selected_labeled_tree']['rate']):.1f}%）。這支持 technical concordance，但不證明 cellular clone identity 或 biological ancestry。<p><a href="HCC1395.html">開啟 HCC1395_HKU 工作站 →</a>　<a href="HCC1395_DORADO.html">開啟 HCC1395_NYGC 工作站 →</a></p></div>
</div></section>
<section class="section" id="provenance"><div class="shell">
  <div class="section-head"><div><div class="eyebrow">Fail-closed provenance</div><h2>JSON 路徑不干擾主要數字</h2></div><p>原始 JSON、receipt 與 SHA-256 收在下方；builder 另逐列核對 MLHP/topology/census keys、row conservation 與 source identities。</p></div>
  <details><summary>機器證據與原始 JSON（預設收合）</summary><div class="details-body">{path_html}<span class="path"><b>cohort receipt SHA-256</b> <code>{authority['cohort_receipt_sha256']}</code></span><span class="path"><b>census receipt SHA-256</b> <code>{authority['census_receipt_sha256']}</code></span><span class="path"><b>similarity sidecar SHA-256</b> <code>{authority['similarity_sha256']}</code></span><span class="path"><b>gene/drug sidecar SHA-256</b> <code>{authority['gene_drug_sha256']}</code></span><span class="path"><b>verified pipeline bound identities</b> <code>{authority['pipeline_identity_count']}</code></span><span class="path"><b>verified census identities</b> <code>{authority['census_identity_count']}</code></span><span class="path"><b>verified similarity identities</b> <code>{authority['similarity_identity_count']}</code></span></div></details>
  <details><summary>歷史 50 kb baseline（不進入新統計）</summary><div class="details-body"><p>舊 canonical-v5 工作站用 proximity/transitive grouping 與不同 candidate-tree schema。它僅保留作方法演進說明；所有本頁 KPI、分布與樣本頁皆來自 2026-07-24 exact-PS authority。</p></div></details>
</div></section>
</main>
<footer class="footer"><div class="shell">InterSubMod · {AUTHORITY_NAME} · UI contract {UI_CONTRACT}</div></footer>
<script type="application/json" id="cohortData">{cohort_payload_json}</script>
<script>{cohort_js}</script>
</body></html>"""


def annotation_summary_from_sample_prefix(page: Path) -> dict[str, int]:
    prefix = page.read_text(encoding="utf-8", errors="replace")[:16000]
    values: dict[str, int] = {}
    for key, meta_name in (
        ("cgcRegions", "intersubmod-cgc-locus-regions"),
        ("cgcDrugRegions", "intersubmod-cgc-drug-locus-regions"),
    ):
        match = re.search(
            rf'<meta name="{re.escape(meta_name)}" content="(\d+)">',
            prefix,
        )
        require(match is not None, f"stale/missing annotation meta in {page}: {meta_name}")
        values[key] = int(match.group(1))
    return values


def build(
    index_only: bool = False,
    verify_only: bool = False,
    retain_rows: bool = False,
) -> dict[str, Any]:
    authority = load_authority()
    topology_map = sample_record_map(authority["all7_summary"]["samples"])
    census_map = sample_record_map(authority["census_summary"]["samples"])
    funnel_samples = authority["funnel"]["samples"]
    require(set(funnel_samples) == set(EXPECTED_SAMPLES), "funnel sample set mismatch")

    if index_only:
        samples = []
        for sample in EXPECTED_SAMPLES:
            page = OUTDIR / f"{sample}.html"
            require(page.is_file(), f"index-only requires existing {page}")
            samples.append(
                {
                    "sample": sample,
                    "summary": {
                        "groups": int(topology_map[sample]["groups_total"]),
                        "ranked": int(census_map[sample]["ranked_units"]),
                        **annotation_summary_from_sample_prefix(page),
                    },
                }
            )
    else:
        samples = []
        for sample in EXPECTED_SAMPLES:
            item = load_sample(
                sample,
                authority,
                topology_map[sample],
                census_map[sample],
                funnel_samples[sample],
            )
            if not verify_only:
                atomic_write(OUTDIR / f"{sample}.html", sample_page(authority, item))
            samples.append(
                item
                if retain_rows
                else {"sample": item["sample"], "summary": item["summary"]}
            )

    require(
        sum(item["summary"]["groups"] for item in samples) == 98955,
        "all-page group count does not conserve 98,955",
    )
    require(
        sum(item["summary"]["ranked"] for item in samples) == 71955,
        "all-page ranked count does not conserve 71,955",
    )
    require(
        sum(int(item["summary"]["cgcRegions"]) for item in samples) == 3554
        and sum(int(item["summary"]["cgcDrugRegions"]) for item in samples)
        == 1252,
        "all-page exact-locus CGC / same-gene CGC∩drug counts do not conserve "
        "3,554 / 1,252",
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
    print(f"INPUT_GENE_DRUG={GENE_DRUG_PATH}")
    if not args.verify_only:
        print(f"OUTPUT={OUTDIR / 'index.html'}")
        if args.index_only:
            print("OUTPUT_SAMPLE_PAGES=0")
            print("VERIFIED_EXISTING_SAMPLE_PAGES=7")
        else:
            print("OUTPUT_SAMPLE_PAGES=7")
    print(
        "CENSUS_RECEIPT_SHA256="
        f"{result['authority']['census_receipt_sha256']}"
    )
    print(
        "GENE_DRUG_SHA256="
        f"{result['authority']['gene_drug_sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
