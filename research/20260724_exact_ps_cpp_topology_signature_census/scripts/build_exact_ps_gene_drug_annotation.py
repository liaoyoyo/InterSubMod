#!/usr/bin/env python3
"""Build a deterministic exact-PS cancer-gene/drug annotation sidecar.

The positional authority is the GENCODE v46 GRCh38 gene body. COSMIC Cancer
Gene Census coordinates are retained only as a cross-check. Drug rows are
limited to the strict DGIdb approved=TRUE AND anti_neoplastic=TRUE subset.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence


SCHEMA_NAME = "intersubmod.exact_ps_gene_drug_annotation"
SCHEMA_VERSION = "1.1.0"
RECEIPT_SCHEMA_NAME = f"{SCHEMA_NAME}.receipt"
RECEIPT_SCHEMA_VERSION = "1.1.0"
EXPECTED_CHROMOSOMES = tuple(f"chr{value}" for value in range(1, 23))

EXPECTED_FILE_NAMES = {
    "gencode": "gencode.v46.basic.annotation.gtf.gz",
    "cgc": "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz",
    "cosmic_genes": "Cosmic_Genes_v104_GRCh38.tsv.gz",
    "dgidb": "dgidb_interactions.tsv",
}

CGC_HEADER = (
    "GENE_SYMBOL",
    "NAME",
    "COSMIC_GENE_ID",
    "CHROMOSOME",
    "GENOME_START",
    "GENOME_STOP",
    "CHR_BAND",
    "SOMATIC",
    "GERMLINE",
    "TUMOUR_TYPES_SOMATIC",
    "TUMOUR_TYPES_GERMLINE",
    "CANCER_SYNDROME",
    "TISSUE_TYPE",
    "MOLECULAR_GENETICS",
    "ROLE_IN_CANCER",
    "MUTATION_TYPES",
    "TRANSLOCATION_PARTNER",
    "OTHER_GERMLINE_MUT",
    "OTHER_SYNDROME",
    "TIER",
    "SYNONYMS",
)

COSMIC_GENES_HEADER = (
    "COSMIC_GENE_ID",
    "LEGACY_GENE_ID",
    "GENE_SYMBOL",
    "GENE_ACCESSION",
    "ENTREZ_ID",
    "HGNC_ID",
    "IN_CANCER_CENSUS",
    "IS_EXPERT_CURATED",
)

DGIDB_HEADER = (
    "gene_claim_name",
    "gene_concept_id",
    "gene_name",
    "interaction_source_db_name",
    "interaction_source_db_version",
    "interaction_type",
    "interaction_score",
    "drug_claim_name",
    "drug_concept_id",
    "drug_name",
    "approved",
    "immunotherapy",
    "anti_neoplastic",
)

CLAIM_CEILING = {
    "level": "gene_level_annotation_only",
    "statement": (
        "This sidecar shows a gene-level intersection between COSMIC Cancer "
        "Gene Census membership and DGIdb rows marked approved=TRUE and "
        "anti_neoplastic=TRUE. It does not establish variant-specific or "
        "tumour-specific clinical actionability."
    ),
    "prohibited_claims": [
        "clinical_actionability",
        "treatment_recommendation",
        "variant_specific_sensitivity_or_resistance",
        "tumour_specific_efficacy",
        "dose_or_eligibility",
    ],
}


class AnnotationBuildError(RuntimeError):
    """Raised when an annotation source or output violates the contract."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AnnotationBuildError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_identity(path: Path, source_key: str) -> dict[str, Any]:
    require(source_key in EXPECTED_FILE_NAMES, f"unknown source key: {source_key}")
    resolved = path.resolve(strict=True)
    require(resolved.is_file(), f"source is not a regular file: {resolved}")
    expected_name = EXPECTED_FILE_NAMES[source_key]
    require(
        path.name == expected_name,
        f"{source_key} source basename must be {expected_name}, got {path.name}",
    )
    stat = resolved.stat()
    return {
        "file_name": path.name,
        "size_bytes": stat.st_size,
        "sha256": sha256_path(resolved),
    }


def encode_json(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def normalize_hgnc(value: str) -> str | None:
    cleaned = value.strip()
    if cleaned in {"", "NULL"}:
        return None
    match = re.fullmatch(r"(?:HGNC:|hgnc:)?([0-9]+)", cleaned)
    require(match is not None, f"invalid HGNC identifier: {value!r}")
    return f"HGNC:{match.group(1)}"


def normalize_symbol(value: str) -> str:
    return value.strip().upper()


def looks_like_hgnc_gene_symbol(value: str) -> bool:
    """Reject pathway/family/free-text claims before symbol fallback."""

    cleaned = value.strip()
    return (
        1 <= len(cleaned) <= 40
        and re.fullmatch(r"[A-Za-z][A-Za-z0-9.-]*", cleaned) is not None
    )


def is_autosomal_chromosome(chrom: str) -> bool:
    return chrom.startswith("chr") and chrom[3:].isdigit() and (
        1 <= int(chrom[3:]) <= 22
    )


def parse_gtf_attributes(raw: str, line_number: int) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for piece in raw.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        fields = piece.split(None, 1)
        require(
            len(fields) == 2,
            f"GENCODE line {line_number}: malformed attribute {piece!r}",
        )
        key, value = fields
        value = value.strip()
        if value.startswith('"'):
            require(
                len(value) >= 2 and value.endswith('"'),
                f"GENCODE line {line_number}: unterminated quoted attribute",
            )
            value = value[1:-1]
        attributes.setdefault(key, value)
    return attributes


def load_gencode(path: Path) -> tuple[list[dict[str, Any]], dict[str, int]]:
    expected_description = (
        "##description: evidence-based annotation of the human genome "
        "(GRCh38), version 46 (Ensembl 112)"
    )
    headers: list[str] = []
    records: list[dict[str, Any]] = []
    seen_gene_ids: set[str] = set()
    total_gene_rows = 0
    excluded_non_autosomal = 0

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            require(
                line.endswith("\n"),
                f"GENCODE line {line_number}: missing terminal newline",
            )
            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
                continue
            fields = line.rstrip("\n").split("\t")
            require(
                len(fields) == 9,
                f"GENCODE line {line_number}: expected 9 GTF fields",
            )
            chrom, _, feature, start_raw, end_raw, _, strand, _, raw_attributes = (
                fields
            )
            try:
                start = int(start_raw)
                end = int(end_raw)
            except ValueError as exc:
                raise AnnotationBuildError(
                    f"GENCODE line {line_number}: non-integer coordinate"
                ) from exc
            require(
                1 <= start <= end,
                f"GENCODE line {line_number}: invalid coordinate interval",
            )
            require(
                strand in {"+", "-", "."},
                f"GENCODE line {line_number}: invalid strand {strand!r}",
            )
            if feature != "gene":
                continue
            total_gene_rows += 1
            if not is_autosomal_chromosome(chrom):
                excluded_non_autosomal += 1
                continue
            attributes = parse_gtf_attributes(raw_attributes, line_number)
            for required_key in ("gene_id", "gene_name", "gene_type"):
                require(
                    attributes.get(required_key, "").strip() != "",
                    f"GENCODE line {line_number}: missing {required_key}",
                )
            gene_id = attributes["gene_id"].strip()
            require(
                gene_id not in seen_gene_ids,
                f"GENCODE duplicate autosomal gene_id: {gene_id}",
            )
            seen_gene_ids.add(gene_id)
            records.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": attributes["gene_name"].strip(),
                    "gene_type": attributes["gene_type"].strip(),
                    "hgnc_id": normalize_hgnc(attributes.get("hgnc_id", "")),
                }
            )

    require(headers, "GENCODE header is missing")
    require(
        headers[0] == expected_description,
        "unsupported GENCODE description/version/build",
    )
    require("##provider: GENCODE" in headers, "GENCODE provider header mismatch")
    require("##format: gtf" in headers, "GENCODE format header mismatch")
    require(records, "GENCODE has no autosomal gene-body records")
    records.sort(
        key=lambda row: (
            int(row["chrom"][3:]),
            row["start"],
            row["end"],
            row["gene_id"],
        )
    )
    return records, {
        "gencode_gene_rows": total_gene_rows,
        "gencode_autosomal_gene_bodies": len(records),
        "gencode_excluded_non_autosomal_gene_bodies": excluded_non_autosomal,
    }


def iter_tsv(
    path: Path,
    expected_header: Sequence[str],
    source_name: str,
    *,
    compressed: bool,
) -> Iterable[tuple[int, dict[str, str]]]:
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(
            tuple(reader.fieldnames or ()) == tuple(expected_header),
            f"{source_name} schema mismatch: got {reader.fieldnames!r}",
        )
        for line_number, row in enumerate(reader, 2):
            require(
                None not in row,
                f"{source_name} line {line_number}: unexpected extra fields",
            )
            require(
                all(value is not None for value in row.values()),
                f"{source_name} line {line_number}: missing fields",
            )
            yield line_number, row


def load_cosmic_genes(path: Path) -> dict[str, dict[str, str | None]]:
    bridge: dict[str, dict[str, str | None]] = {}
    for line_number, row in iter_tsv(
        path,
        COSMIC_GENES_HEADER,
        "COSMIC Genes",
        compressed=True,
    ):
        cosmic_gene_id = row["COSMIC_GENE_ID"].strip()
        symbol = row["GENE_SYMBOL"].strip()
        require(
            re.fullmatch(r"COSG[0-9]+", cosmic_gene_id) is not None,
            f"COSMIC Genes line {line_number}: invalid COSMIC_GENE_ID",
        )
        require(
            symbol != "",
            f"COSMIC Genes line {line_number}: empty GENE_SYMBOL",
        )
        require(
            cosmic_gene_id not in bridge,
            f"COSMIC Genes duplicate COSMIC_GENE_ID: {cosmic_gene_id}",
        )
        bridge[cosmic_gene_id] = {
            "symbol": symbol,
            "hgnc_id": normalize_hgnc(row["HGNC_ID"]),
        }
    require(bridge, "COSMIC Genes bridge is empty")
    return bridge


def load_cgc(
    path: Path,
    bridge: dict[str, dict[str, str | None]],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    rows: list[dict[str, Any]] = []
    total_rows = 0
    excluded_non_autosomal = 0
    missing_coordinates = 0
    seen_cosmic_ids: set[str] = set()

    for line_number, row in iter_tsv(
        path,
        CGC_HEADER,
        "COSMIC Cancer Gene Census",
        compressed=True,
    ):
        total_rows += 1
        cosmic_gene_id = row["COSMIC_GENE_ID"].strip()
        symbol = row["GENE_SYMBOL"].strip()
        chromosome_raw = row["CHROMOSOME"].strip()
        require(
            re.fullmatch(r"COSG[0-9]+", cosmic_gene_id) is not None,
            f"CGC line {line_number}: invalid COSMIC_GENE_ID",
        )
        require(symbol != "", f"CGC line {line_number}: empty GENE_SYMBOL")
        require(
            cosmic_gene_id not in seen_cosmic_ids,
            f"CGC duplicate COSMIC_GENE_ID: {cosmic_gene_id}",
        )
        seen_cosmic_ids.add(cosmic_gene_id)
        require(
            cosmic_gene_id in bridge,
            f"CGC line {line_number}: COSG missing from COSMIC Genes bridge",
        )
        bridge_row = bridge[cosmic_gene_id]
        require(
            normalize_symbol(str(bridge_row["symbol"])) == normalize_symbol(symbol),
            f"CGC line {line_number}: COSG symbol mismatch",
        )
        try:
            tier = int(row["TIER"])
        except ValueError as exc:
            raise AnnotationBuildError(
                f"CGC line {line_number}: invalid TIER"
            ) from exc
        require(tier in {1, 2}, f"CGC line {line_number}: unsupported TIER")

        if not chromosome_raw.isdigit() or not (
            1 <= int(chromosome_raw) <= 22
        ):
            excluded_non_autosomal += 1
            continue
        start_raw = row["GENOME_START"].strip()
        end_raw = row["GENOME_STOP"].strip()
        require(
            (start_raw == "") == (end_raw == ""),
            f"CGC line {line_number}: incomplete coordinate pair",
        )
        cgc_start: int | None = None
        cgc_end: int | None = None
        if start_raw == "":
            missing_coordinates += 1
        else:
            try:
                cgc_start = int(start_raw)
                cgc_end = int(end_raw)
            except ValueError as exc:
                raise AnnotationBuildError(
                    f"CGC line {line_number}: non-integer coordinate"
                ) from exc
            require(
                1 <= cgc_start <= cgc_end,
                f"CGC line {line_number}: invalid coordinate interval",
            )
        rows.append(
            {
                "symbol": symbol,
                "cosmic_gene_id": cosmic_gene_id,
                "hgnc_id": bridge_row["hgnc_id"],
                "chrom": f"chr{chromosome_raw}",
                "cgc_start": cgc_start,
                "cgc_end": cgc_end,
                "cgc_tier": tier,
                "cgc_role": row["ROLE_IN_CANCER"].strip(),
                "cgc_tumour": row["TUMOUR_TYPES_SOMATIC"].strip(),
            }
        )
    require(rows, "CGC has no autosomal rows")
    return rows, {
        "cgc_source_rows": total_rows,
        "cgc_autosomal_rows": len(rows),
        "cgc_excluded_non_autosomal": excluded_non_autosomal,
        "cgc_missing_coordinates": missing_coordinates,
    }


def load_dgidb(
    path: Path,
    allowed_gencode_symbols: set[str],
) -> tuple[
    dict[str, list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
    dict[str, int],
]:
    by_hgnc: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    total_rows = 0
    selected_rows = 0
    selected_without_indexable_gene = 0
    selected_rows_without_hgnc_id = 0
    symbol_only_rows_accepted = 0
    allowed_boolean_values = {"TRUE", "FALSE", "NULL"}

    for line_number, row in iter_tsv(
        path,
        DGIDB_HEADER,
        "DGIdb interactions",
        compressed=False,
    ):
        total_rows += 1
        for key in ("approved", "immunotherapy", "anti_neoplastic"):
            require(
                row[key] in allowed_boolean_values,
                f"DGIdb line {line_number}: invalid {key} boolean",
            )
        if not (
            row["approved"] == "TRUE"
            and row["anti_neoplastic"] == "TRUE"
        ):
            continue
        selected_rows += 1
        for key in (
            "drug_claim_name",
            "drug_concept_id",
            "drug_name",
            "interaction_source_db_name",
            "interaction_source_db_version",
        ):
            require(
                row[key].strip() not in {"", "NULL"},
                f"DGIdb line {line_number}: selected row lacks {key}",
            )
        evidence = {
            "claim_name": row["drug_claim_name"].strip(),
            "canonical_name": row["drug_name"].strip(),
            "concept_id": row["drug_concept_id"].strip(),
            "source_name": row["interaction_source_db_name"].strip(),
            "source_version": row["interaction_source_db_version"].strip(),
        }
        gene_concept = row["gene_concept_id"].strip()
        hgnc_match = re.fullmatch(r"hgnc:([0-9]+)", gene_concept, re.IGNORECASE)
        if hgnc_match:
            by_hgnc[f"HGNC:{hgnc_match.group(1)}"].append(evidence)
            continue
        selected_rows_without_hgnc_id += 1
        gene_name = row["gene_name"].strip()
        gene_claim_name = row["gene_claim_name"].strip()
        symbol = (
            gene_name
            if gene_name not in {"", "NULL"}
            else gene_claim_name
        )
        if (
            symbol in {"", "NULL"}
            or not looks_like_hgnc_gene_symbol(symbol)
            or normalize_symbol(symbol) not in allowed_gencode_symbols
        ):
            selected_without_indexable_gene += 1
            continue
        symbol_only_rows_accepted += 1
        by_symbol[normalize_symbol(symbol)].append(evidence)

    require(total_rows > 0, "DGIdb interactions source is empty")
    return by_hgnc, by_symbol, {
        "dgidb_source_rows": total_rows,
        "dgidb_approved_and_antineoplastic_rows": selected_rows,
        # This is the raw non-HGNC subset evaluated for symbol-only indexing.
        "dgidb_selected_rows_indexed_by_symbol_only": (
            selected_rows_without_hgnc_id
        ),
        "dgidb_selected_rows_without_hgnc_id": (
            selected_rows_without_hgnc_id
        ),
        "dgidb_symbol_only_rows_accepted_after_gencode_gate": (
            symbol_only_rows_accepted
        ),
        "dgidb_symbol_only_keys_accepted_after_gencode_gate": len(by_symbol),
        "dgidb_selected_rows_without_indexable_gene": (
            selected_without_indexable_gene
        ),
    }


def group_drugs(evidence_rows: Iterable[dict[str, str]]) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, set[Any]]] = {}
    for row in evidence_rows:
        claim_name = row["claim_name"]
        bucket = grouped.setdefault(
            claim_name,
            {
                "canonical_names": set(),
                "concept_ids": set(),
                "sources": set(),
            },
        )
        bucket["canonical_names"].add(row["canonical_name"])
        bucket["concept_ids"].add(row["concept_id"])
        bucket["sources"].add((row["source_name"], row["source_version"]))

    result: list[dict[str, Any]] = []
    for claim_name in sorted(grouped, key=lambda value: (value.casefold(), value)):
        bucket = grouped[claim_name]
        result.append(
            {
                "claim_name": claim_name,
                "canonical_names": sorted(
                    bucket["canonical_names"],
                    key=lambda value: (value.casefold(), value),
                ),
                "concept_ids": sorted(bucket["concept_ids"]),
                "sources": [
                    {"name": name, "version": version}
                    for name, version in sorted(
                        bucket["sources"],
                        key=lambda value: (
                            value[0].casefold(),
                            value[0],
                            value[1],
                        ),
                    )
                ],
            }
        )
    return result


def coordinate_crosscheck(
    cgc: dict[str, Any],
    gencode: dict[str, Any],
) -> dict[str, Any]:
    if cgc["cgc_start"] is None:
        return {"status": "cgc_coordinate_missing"}
    same_chromosome = cgc["chrom"] == gencode["chrom"]
    exact = (
        same_chromosome
        and cgc["cgc_start"] == gencode["start"]
        and cgc["cgc_end"] == gencode["end"]
    )
    overlaps = (
        same_chromosome
        and max(cgc["cgc_start"], gencode["start"])
        <= min(cgc["cgc_end"], gencode["end"])
    )
    status = (
        "exact"
        if exact
        else "overlap_nonidentical"
        if overlaps
        else "chromosome_mismatch"
        if not same_chromosome
        else "disjoint"
    )
    return {
        "status": status,
        "chrom": cgc["chrom"],
        "start": cgc["cgc_start"],
        "end": cgc["cgc_end"],
    }


def derive_receipt_checks(
    annotation: dict[str, Any],
    source_identities: dict[str, dict[str, Any]],
    output_identity: dict[str, Any],
) -> dict[str, bool]:
    """Derive every publication check from the built payload."""

    scope = annotation["scope"]
    counts = annotation["counts"]
    genes = annotation["genes"]
    expected_source_keys = set(EXPECTED_FILE_NAMES)
    source_identity_valid = (
        set(source_identities) == expected_source_keys
        and all(
            identity.get("file_name") == EXPECTED_FILE_NAMES[key]
            and isinstance(identity.get("size_bytes"), int)
            and identity["size_bytes"] > 0
            and re.fullmatch(r"[0-9a-f]{64}", str(identity.get("sha256")))
            is not None
            for key, identity in source_identities.items()
        )
    )
    coordinate_partition = sum(
        counts[key]
        for key in (
            "cgc_gencode_coordinate_exact",
            "cgc_gencode_coordinate_overlap_nonidentical",
            "cgc_gencode_coordinate_missing",
        )
    )
    annotation_bytes = encode_json(annotation)
    output_identity_valid = (
        output_identity.get("path") == "annotation.v1.json"
        and output_identity.get("size_bytes") == len(annotation_bytes)
        and output_identity.get("sha256")
        == hashlib.sha256(annotation_bytes).hexdigest()
    )
    return {
        "schema_supported": (
            annotation.get("schema_name") == SCHEMA_NAME
            and annotation.get("schema_version") == SCHEMA_VERSION
        ),
        "scope_is_chr1_through_chr22": (
            scope.get("genome_build") == "GRCh38"
            and scope.get("chromosomes") == list(EXPECTED_CHROMOSOMES)
            and scope.get("partial") is False
        ),
        "gencode_v46_gene_body_is_position_authority": (
            scope.get("coordinate_system")
            == "GENCODE_v46_gene_body_1_based_inclusive"
            and len(genes) == counts["cgc_genes"]
            and all(
                row["chrom"] in EXPECTED_CHROMOSOMES
                and 1 <= row["start"] <= row["end"]
                and row["gene_id"]
                and row["gene_name"]
                for row in genes
            )
        ),
        "cosg_to_hgnc_primary_join": (
            counts["cgc_genes"]
            == counts["cgc_hgnc_primary_gene_matches"]
            + counts["cgc_symbol_fallback_gene_matches"]
            and counts["cgc_autosomal_rows"]
            == counts["cgc_genes"]
            + counts["cgc_without_unambiguous_gencode_gene_body"]
        ),
        "symbol_fallback_only_when_hgnc_missing": (
            counts["dgidb_selected_rows_indexed_by_symbol_only"]
            == counts["dgidb_symbol_only_rows_accepted_after_gencode_gate"]
            + counts["dgidb_selected_rows_without_indexable_gene"]
            and counts["dgidb_symbol_only_keys_overlapping_cgc"]
            <= counts["dgidb_symbol_only_keys_accepted_after_gencode_gate"]
            and all(
                row["hgnc_id"] is None
                for row in genes
                if row["drug_match_basis"]
                == "symbol_fallback_missing_hgnc"
            )
        ),
        "cgc_coordinate_crosscheck_passed": (
            coordinate_partition == counts["cgc_genes"]
            and all(
                row["cgc_coordinate_crosscheck"]["status"]
                in {"exact", "overlap_nonidentical", "cgc_coordinate_missing"}
                for row in genes
            )
        ),
        "dgidb_approved_and_antineoplastic_and_filter": (
            0
            <= counts["dgidb_approved_and_antineoplastic_rows"]
            <= counts["dgidb_source_rows"]
            and counts["cgc_with_approved_antineoplastic"]
            <= counts["cgc_genes"]
            and all(
                drug.get("claim_name")
                and drug.get("sources")
                and all(
                    source.get("name") and source.get("version")
                    for source in drug["sources"]
                )
                for row in genes
                for drug in row["approved_antineoplastic_drugs"]
            )
        ),
        "raw_source_identities_recorded": source_identity_valid,
        "clinical_claim_ceiling_present": (
            annotation.get("claim_ceiling", {}).get("level")
            == "gene_level_annotation_only"
            and "clinical_actionability"
            in annotation.get("claim_ceiling", {}).get(
                "prohibited_claims",
                [],
            )
        ),
        "deterministic_serialization_contract": output_identity_valid,
    }


def build(
    *,
    gencode_path: Path,
    cgc_path: Path,
    cosmic_genes_path: Path,
    dgidb_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    source_paths = {
        "gencode": gencode_path,
        "cgc": cgc_path,
        "cosmic_genes": cosmic_genes_path,
        "dgidb": dgidb_path,
    }
    source_identities = {
        key: source_identity(path, key)
        for key, path in source_paths.items()
    }

    gencode_rows, gencode_counts = load_gencode(gencode_path)
    bridge = load_cosmic_genes(cosmic_genes_path)
    cgc_rows, cgc_counts = load_cgc(cgc_path, bridge)
    allowed_gencode_symbols = {
        normalize_symbol(row["gene_name"]) for row in gencode_rows
    }
    dgidb_hgnc, dgidb_symbol, dgidb_counts = load_dgidb(
        dgidb_path,
        allowed_gencode_symbols,
    )

    gencode_by_hgnc: dict[str, list[dict[str, Any]]] = defaultdict(list)
    gencode_by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in gencode_rows:
        if row["hgnc_id"] is not None:
            gencode_by_hgnc[row["hgnc_id"]].append(row)
        gencode_by_symbol[normalize_symbol(row["gene_name"])].append(row)

    genes: list[dict[str, Any]] = []
    seen_gencode_ids: set[str] = set()
    match_counts = {
        "cgc_hgnc_primary_gene_matches": 0,
        "cgc_symbol_fallback_gene_matches": 0,
        "cgc_without_unambiguous_gencode_gene_body": 0,
        "gencode_cgc_symbol_mismatches": 0,
        "cgc_gencode_coordinate_exact": 0,
        "cgc_gencode_coordinate_overlap_nonidentical": 0,
        "cgc_gencode_coordinate_missing": 0,
        "dgidb_hgnc_primary_gene_matches": 0,
        "dgidb_symbol_fallback_gene_matches": 0,
    }
    matched_evidence_rows = 0

    for cgc in cgc_rows:
        if cgc["hgnc_id"] is not None:
            candidates = gencode_by_hgnc.get(cgc["hgnc_id"], [])
            match_basis = "cosg_to_hgnc"
        else:
            candidates = gencode_by_symbol.get(
                normalize_symbol(cgc["symbol"]),
                [],
            )
            match_basis = "symbol_fallback_missing_hgnc"
        if len(candidates) != 1:
            match_counts[
                "cgc_without_unambiguous_gencode_gene_body"
            ] += 1
            continue
        gencode = candidates[0]
        if match_basis == "cosg_to_hgnc":
            match_counts["cgc_hgnc_primary_gene_matches"] += 1
        else:
            match_counts["cgc_symbol_fallback_gene_matches"] += 1
        if normalize_symbol(cgc["symbol"]) != normalize_symbol(
            gencode["gene_name"]
        ):
            match_counts["gencode_cgc_symbol_mismatches"] += 1

        crosscheck = coordinate_crosscheck(cgc, gencode)
        require(
            crosscheck["status"]
            in {"exact", "overlap_nonidentical", "cgc_coordinate_missing"},
            (
                "CGC/GENCODE coordinate cross-check failed for "
                f"{cgc['symbol']}: {crosscheck['status']}"
            ),
        )
        coordinate_count_key = {
            "exact": "cgc_gencode_coordinate_exact",
            "overlap_nonidentical": (
                "cgc_gencode_coordinate_overlap_nonidentical"
            ),
            "cgc_coordinate_missing": "cgc_gencode_coordinate_missing",
        }[crosscheck["status"]]
        match_counts[coordinate_count_key] += 1
        require(
            gencode["gene_id"] not in seen_gencode_ids,
            f"multiple CGC rows map to GENCODE gene {gencode['gene_id']}",
        )
        seen_gencode_ids.add(gencode["gene_id"])

        effective_hgnc = cgc["hgnc_id"] or gencode["hgnc_id"]
        if effective_hgnc is not None:
            evidence = dgidb_hgnc.get(effective_hgnc, [])
            drug_match_basis = "hgnc_id"
            if evidence:
                match_counts["dgidb_hgnc_primary_gene_matches"] += 1
        else:
            evidence = dgidb_symbol.get(
                normalize_symbol(gencode["gene_name"]),
                [],
            )
            drug_match_basis = "symbol_fallback_missing_hgnc"
            if evidence:
                match_counts["dgidb_symbol_fallback_gene_matches"] += 1
        matched_evidence_rows += len(evidence)
        drugs = group_drugs(evidence)

        genes.append(
            {
                "chrom": gencode["chrom"],
                "start": gencode["start"],
                "end": gencode["end"],
                "gene_id": gencode["gene_id"],
                "gene_name": gencode["gene_name"],
                "gene_type": gencode["gene_type"],
                "strand": gencode["strand"],
                "hgnc_id": effective_hgnc,
                "symbol": cgc["symbol"],
                "cgc_cosmic_gene_id": cgc["cosmic_gene_id"],
                "cgc_tier": cgc["cgc_tier"],
                "cgc_role": cgc["cgc_role"],
                "cgc_tumour": cgc["cgc_tumour"],
                "gene_match_basis": match_basis,
                "cgc_coordinate_crosscheck": crosscheck,
                "drug_match_basis": drug_match_basis,
                "approved_antineoplastic_drugs": drugs,
            }
        )

    genes.sort(
        key=lambda row: (
            int(row["chrom"][3:]),
            row["start"],
            row["end"],
            row["gene_id"],
        )
    )
    with_drugs = sum(
        bool(row["approved_antineoplastic_drugs"]) for row in genes
    )
    cgc_gencode_symbols = {
        normalize_symbol(row["gene_name"]) for row in genes
    }
    counts = {
        **gencode_counts,
        **cgc_counts,
        **dgidb_counts,
        **match_counts,
        "cgc_genes": len(genes),
        "cgc_with_approved_antineoplastic": with_drugs,
        "cgc_without_approved_antineoplastic": len(genes) - with_drugs,
        "matched_dgidb_evidence_rows_before_collapse": matched_evidence_rows,
        "gene_drug_claim_pairs_after_collapse": sum(
            len(row["approved_antineoplastic_drugs"]) for row in genes
        ),
        "dgidb_symbol_only_keys_overlapping_cgc": len(
            set(dgidb_symbol) & cgc_gencode_symbols
        ),
    }
    sources = {
        "gencode": {
            **source_identities["gencode"],
            "name": "GENCODE basic gene annotation",
            "version": "v46 / Ensembl 112",
            "role": "GRCh38 1-based inclusive gene-body positional authority",
        },
        "cgc": {
            **source_identities["cgc"],
            "name": "COSMIC Cancer Gene Census",
            "version": "v104",
            "role": "cancer-gene membership and descriptive metadata",
        },
        "cosmic_genes": {
            **source_identities["cosmic_genes"],
            "name": "COSMIC Genes",
            "version": "v104",
            "role": "COSG-to-HGNC identifier bridge",
        },
        "dgidb": {
            **source_identities["dgidb"],
            "name": "DGIdb interactions export",
            "version": "row-level source versions preserved per drug claim",
            "role": "approved AND anti-neoplastic gene-drug association subset",
        },
    }
    annotation = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "scope": {
            "genome_build": "GRCh38",
            "chromosomes": list(EXPECTED_CHROMOSOMES),
            "coordinate_system": "GENCODE_v46_gene_body_1_based_inclusive",
            "partial": False,
            "intersection_contract": (
                "position_in_GENCODE_gene_body AND same_gene_in_CGC; "
                "drug subset additionally requires DGIdb approved=TRUE AND "
                "anti_neoplastic=TRUE"
            ),
            "symbol_fallback_contract": (
                "only_when_HGNC_is_missing AND token_looks_like_HGNC_symbol "
                "AND token_exists_as_an_exact_GENCODE_v46_gene_name"
            ),
        },
        "sources": sources,
        "claim_ceiling": CLAIM_CEILING,
        "counts": counts,
        "genes": genes,
    }
    annotation_bytes = encode_json(annotation)
    output_identity = {
        "path": "annotation.v1.json",
        "size_bytes": len(annotation_bytes),
        "sha256": hashlib.sha256(annotation_bytes).hexdigest(),
    }
    checks = derive_receipt_checks(
        annotation,
        source_identities,
        output_identity,
    )
    require(
        all(checks.values()),
        "derived publication checks failed: "
        + ",".join(key for key, value in checks.items() if not value),
    )
    receipt = {
        "schema_name": RECEIPT_SCHEMA_NAME,
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "all_pass": all(checks.values()),
        "checks": checks,
        "scope": annotation["scope"],
        "raw_sources": source_identities,
        "counts": counts,
        "claim_ceiling": CLAIM_CEILING,
        "output": output_identity,
        "outputs": {
            "annotation": output_identity,
        },
        "determinism": {
            "json_keys_sorted": True,
            "array_order": (
                "genes by chromosome/start/end/gene_id; drugs by exact "
                "drug_claim_name; evidence sources by name/version"
            ),
            "timestamp_omitted": True,
        },
    }
    return annotation, receipt


def publish(
    output_dir: Path,
    annotation: dict[str, Any],
    receipt: dict[str, Any],
) -> None:
    require(not output_dir.exists(), f"output directory already exists: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    outputs = {
        output_dir / "annotation.v1.json": encode_json(annotation),
        output_dir / "receipt.json": encode_json(receipt),
    }
    for path, payload in outputs.items():
        with path.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gencode", type=Path, required=True)
    parser.add_argument("--cgc", type=Path, required=True)
    parser.add_argument("--cosmic-genes", type=Path, required=True)
    parser.add_argument("--dgidb", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    annotation, receipt = build(
        gencode_path=args.gencode,
        cgc_path=args.cgc,
        cosmic_genes_path=args.cosmic_genes,
        dgidb_path=args.dgidb,
    )
    publish(args.output_dir, annotation, receipt)
    counts = annotation["counts"]
    print(
        "OK "
        f"output={args.output_dir.resolve()} "
        f"cgc_genes={counts['cgc_genes']} "
        "cgc_with_approved_antineoplastic="
        f"{counts['cgc_with_approved_antineoplastic']} "
        f"annotation_sha256={receipt['output']['sha256']}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AnnotationBuildError, OSError, UnicodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
