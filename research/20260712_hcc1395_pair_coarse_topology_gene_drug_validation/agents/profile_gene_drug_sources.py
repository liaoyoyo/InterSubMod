#!/usr/bin/env python3
"""Profile local gene/drug sources and join the latest HCC1395 pair regions.

This script is intentionally read-only outside its output directory.  It uses the
historical layered-v2 exact topology catalog selected by the parent analysis and
does not treat any annotation database as topology ground truth.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUT = ROOT / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents"
TOPOLOGY = ROOT / "research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json"
REGION_TSV = ROOT / "research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_regions.tsv"
COARSE_REGIONS = ROOT / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv"
PAIR_MATCHES = ROOT / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv"
ANN = Path("/big7_disk/liaoyoyo2001/gene_annotation")
GTF = ANN / "gencode.v46.basic.annotation.gtf.gz"
CGC = ANN / "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz"
COSMIC_GENES = ANN / "cosmic_v104/Cosmic_Genes_v104_GRCh38.tsv.gz"
RESISTANCE = ANN / "cosmic_v104/Cosmic_ResistanceMutations_v104_GRCh38.tsv.gz"
CLP_MUTANT = ANN / "cosmic_v104/CellLinesProject_GenomeScreensMutant_v104_GRCh38.tsv.gz"
DGIDB = ANN / "dgidb_interactions.tsv"
WAKHAN = Path(
    "/big7_disk/liaoyoyo2001/external_validation/axis2_longread_phasing_haplotag/"
    "wakhan-cancer-cna-loh-haplotype-2025/repo/src/annotations/COSMIC_cancer_genes.tsv"
)
OLD_REGION_ANN = {
    "HCC1395": ROOT / "docs/methodology/_assets/20260627_subclone_4axis_teaching/data/region_gene_annotation.json",
    "HCC1395_DORADO": Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/"
        "HCC1395_DORADO/region_gene_annotation.json"
    ),
}
SAMPLES = ("HCC1395", "HCC1395_DORADO")
BIN_SIZE = 1_000_000
ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def clean(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.upper() == "NULL" else value


def boolish(value: str | None) -> bool:
    return clean(value).upper() == "TRUE"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).astimezone().isoformat()


def pct(num: int, den: int) -> float:
    return round(100.0 * num / den, 4) if den else 0.0


def quantile(values: list[int], q: float) -> float:
    if not values:
        return 0.0
    vals = sorted(values)
    pos = (len(vals) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(vals) - 1)
    frac = pos - lo
    return vals[lo] * (1 - frac) + vals[hi] * frac


def load_cosmic_genes() -> tuple[dict[str, dict[str, str]], dict]:
    by_cosg: dict[str, dict[str, str]] = {}
    rows = 0
    hgnc = set()
    symbols = set()
    duplicate_cosg = 0
    with gzip.open(COSMIC_GENES, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            cosg = clean(row.get("COSMIC_GENE_ID"))
            if cosg in by_cosg:
                duplicate_cosg += 1
            hid = clean(row.get("HGNC_ID"))
            hkey = f"hgnc:{hid}" if hid else ""
            symbol = clean(row.get("GENE_SYMBOL")).upper()
            by_cosg[cosg] = {"hgnc": hkey, "symbol": symbol, "ensembl": clean(row.get("GENE_ACCESSION"))}
            if hkey:
                hgnc.add(hkey)
            if symbol:
                symbols.add(symbol)
    return by_cosg, {
        "rows": rows,
        "unique_cosmic_gene_id": len(by_cosg),
        "duplicate_cosmic_gene_id_rows": duplicate_cosg,
        "unique_hgnc": len(hgnc),
        "unique_symbols": len(symbols),
    }


def load_cgc(cosmic_genes: dict[str, dict[str, str]]) -> tuple[dict[str, dict], dict[str, dict], dict]:
    by_hgnc: dict[str, dict] = {}
    by_symbol: dict[str, dict] = {}
    rows = 0
    duplicate_symbols = 0
    duplicate_hgnc = 0
    missing_hgnc = 0
    roles = Counter()
    tiers = Counter()
    aliases = set()
    with gzip.open(CGC, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            symbol = clean(row.get("GENE_SYMBOL")).upper()
            cosg = clean(row.get("COSMIC_GENE_ID"))
            hkey = cosmic_genes.get(cosg, {}).get("hgnc", "")
            record = {
                "symbol": symbol,
                "cosmic_gene_id": cosg,
                "hgnc": hkey,
                "role": clean(row.get("ROLE_IN_CANCER")),
                "tier": clean(row.get("TIER")),
                "tumour": clean(row.get("TUMOUR_TYPES_SOMATIC")),
                "chrom": clean(row.get("CHROMOSOME")),
                "start": clean(row.get("GENOME_START")),
                "end": clean(row.get("GENOME_STOP")),
            }
            if symbol in by_symbol:
                duplicate_symbols += 1
            by_symbol[symbol] = record
            if hkey:
                if hkey in by_hgnc:
                    duplicate_hgnc += 1
                by_hgnc[hkey] = record
            else:
                missing_hgnc += 1
            tiers[record["tier"] or "(missing)"] += 1
            for role in [x.strip() for x in record["role"].split(",") if x.strip()]:
                roles[role] += 1
            for alias in [x.strip().upper() for x in clean(row.get("SYNONYMS")).split(",") if x.strip()]:
                aliases.add(alias)
    profile = {
        "rows": rows,
        "unique_symbols": len(by_symbol),
        "duplicate_symbol_rows": duplicate_symbols,
        "unique_hgnc": len(by_hgnc),
        "duplicate_hgnc_rows": duplicate_hgnc,
        "missing_hgnc": missing_hgnc,
        "distinct_alias_tokens": len(aliases),
        "tiers": dict(sorted(tiers.items())),
        "role_tokens": dict(sorted(roles.items())),
    }
    return by_hgnc, by_symbol, profile


def load_dgidb() -> tuple[dict[str, dict], dict[str, dict], dict]:
    by_hgnc: dict[str, dict] = defaultdict(
        lambda: {
            "rows": 0,
            "drugs": {},
            "approved_drugs": {},
            "approved_antineoplastic_drugs": {},
            "sources": set(),
            "symbols": set(),
        }
    )
    by_symbol: dict[str, dict] = defaultdict(
        lambda: {
            "rows": 0,
            "drugs": {},
            "approved_drugs": {},
            "approved_antineoplastic_drugs": {},
            "sources": set(),
            "symbols": set(),
        }
    )
    rows = 0
    exact_seen = set()
    composite_seen = set()
    duplicate_composite_rows = 0
    genes = set()
    hgnc = set()
    drug_keys = set()
    sources = Counter()
    source_versions: dict[str, set[str]] = defaultdict(set)
    missing = Counter()
    approved_rows = anti_rows = approved_anti_rows = 0
    approved_genes = set()
    approved_drugs = set()
    approved_pairs = set()
    anti_genes = set()
    anti_drugs = set()
    anti_pairs = set()
    approved_anti_genes = set()
    approved_anti_drugs = set()
    approved_anti_pairs = set()
    with DGIDB.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        for row in reader:
            rows += 1
            full = tuple(row.get(field, "") for field in fields)
            exact_seen.add(full)
            symbol = clean(row.get("gene_name")).upper()
            hkey = clean(row.get("gene_concept_id")).lower()
            source = clean(row.get("interaction_source_db_name"))
            source_version = clean(row.get("interaction_source_db_version"))
            drug_id = clean(row.get("drug_concept_id")).lower()
            canonical_name = clean(row.get("drug_name"))
            claim_name = clean(row.get("drug_claim_name"))
            display_name = claim_name or canonical_name
            dkey = drug_id or f"name:{(canonical_name or claim_name).upper()}"
            approved = boolish(row.get("approved"))
            anti = boolish(row.get("anti_neoplastic"))
            if not symbol:
                missing["gene_name"] += 1
            if not hkey:
                missing["gene_concept_id"] += 1
            if not dkey or dkey == "name:":
                missing["drug_key"] += 1
            if not canonical_name:
                missing["drug_name"] += 1
            if not claim_name:
                missing["drug_claim_name"] += 1
            if symbol:
                genes.add(symbol)
            if hkey:
                hgnc.add(hkey)
            if dkey and dkey != "name:":
                drug_keys.add(dkey)
            if source:
                sources[source] += 1
                if source_version:
                    source_versions[source].add(source_version)
            if approved:
                approved_rows += 1
            if anti:
                anti_rows += 1
            if approved and anti:
                approved_anti_rows += 1
            joinable_gene = hkey or (f"symbol:{symbol}" if symbol else "")
            if approved and joinable_gene and dkey and dkey != "name:":
                approved_genes.add(joinable_gene)
                approved_drugs.add(dkey)
                approved_pairs.add((joinable_gene, dkey))
            if anti and joinable_gene and dkey and dkey != "name:":
                anti_genes.add(joinable_gene)
                anti_drugs.add(dkey)
                anti_pairs.add((joinable_gene, dkey))
            if approved and anti and joinable_gene and dkey and dkey != "name:":
                approved_anti_genes.add(joinable_gene)
                approved_anti_drugs.add(dkey)
                approved_anti_pairs.add((joinable_gene, dkey))
            composite = (hkey or f"symbol:{symbol}", dkey, source, clean(row.get("interaction_type")))
            if composite in composite_seen:
                duplicate_composite_rows += 1
            composite_seen.add(composite)
            payload = {
                "id": dkey,
                "claim_name": display_name,
                "canonical_name": canonical_name,
                "source": source,
            }
            targets = []
            if hkey:
                targets.append(by_hgnc[hkey])
            if symbol:
                targets.append(by_symbol[symbol])
            for target in targets:
                target["rows"] += 1
                target["symbols"].add(symbol)
                if source:
                    target["sources"].add(source)
                if dkey and dkey != "name:":
                    target["drugs"][dkey] = payload
                    if approved:
                        target["approved_drugs"][dkey] = payload
                    if approved and anti:
                        target["approved_antineoplastic_drugs"][dkey] = payload
    profile = {
        "rows": rows,
        "exact_duplicate_rows": rows - len(exact_seen),
        "duplicate_composite_rows": duplicate_composite_rows,
        "unique_gene_symbols": len(genes),
        "unique_gene_concept_ids": len(hgnc),
        "unique_drug_keys": len(drug_keys),
        "approved_rows": approved_rows,
        "approved_unique_genes": len(approved_genes),
        "approved_unique_drugs": len(approved_drugs),
        "approved_unique_gene_drug_pairs": len(approved_pairs),
        "anti_neoplastic_rows": anti_rows,
        "anti_neoplastic_unique_genes": len(anti_genes),
        "anti_neoplastic_unique_drugs": len(anti_drugs),
        "anti_neoplastic_unique_gene_drug_pairs": len(anti_pairs),
        "approved_antineoplastic_rows": approved_anti_rows,
        "approved_antineoplastic_unique_genes": len(approved_anti_genes),
        "approved_antineoplastic_unique_drugs": len(approved_anti_drugs),
        "approved_antineoplastic_unique_gene_drug_pairs": len(approved_anti_pairs),
        "missing": dict(missing),
        "source_db_count": len(sources),
        "source_rows": dict(sources.most_common()),
        "source_versions": {key: sorted(value) for key, value in sorted(source_versions.items())},
    }
    return by_hgnc, by_symbol, profile


def load_gtf() -> tuple[list[dict], dict]:
    genes: list[dict] = []
    total_features = 0
    headers = []
    gene_ids = set()
    symbols = Counter()
    hgnc = set()
    gene_types = Counter()
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                if len(headers) < 8:
                    headers.append(line.rstrip())
                continue
            if line.startswith("#"):
                continue
            total_features += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = dict(ATTR_RE.findall(parts[8]))
            gene_id = clean(attrs.get("gene_id")).split(".")[0]
            symbol = clean(attrs.get("gene_name")).upper()
            hkey = clean(attrs.get("hgnc_id")).lower()
            gene_type = clean(attrs.get("gene_type"))
            rec = {
                "chrom": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "gene_id": gene_id,
                "symbol": symbol,
                "hgnc": hkey,
                "gene_type": gene_type,
            }
            genes.append(rec)
            gene_ids.add(gene_id)
            if symbol:
                symbols[symbol] += 1
            if hkey:
                hgnc.add(hkey)
            gene_types[gene_type or "(missing)"] += 1
    profile = {
        "feature_rows": total_features,
        "gene_rows": len(genes),
        "unique_gene_ids": len(gene_ids),
        "unique_symbols": len(symbols),
        "duplicate_symbol_gene_rows": sum(count - 1 for count in symbols.values() if count > 1),
        "symbols_with_multiple_gene_ids": sum(1 for count in symbols.values() if count > 1),
        "unique_hgnc": len(hgnc),
        "missing_hgnc_gene_rows": sum(1 for gene in genes if not gene["hgnc"]),
        "gene_types": dict(gene_types.most_common()),
        "header": headers,
    }
    return genes, profile


def load_resistance() -> dict:
    rows = 0
    full_seen = set()
    sample_mut_drug = set()
    genomic_drug = set()
    genes = set()
    drugs = set()
    samples = set()
    hcc1395_rows = 0
    with gzip.open(RESISTANCE, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        for row in reader:
            rows += 1
            full_seen.add(tuple(row.get(field, "") for field in fields))
            sample = clean(row.get("SAMPLE_NAME"))
            gene = clean(row.get("GENE_SYMBOL")).upper()
            drug = clean(row.get("DRUG_NAME")).upper()
            mut = clean(row.get("GENOMIC_MUTATION_ID"))
            samples.add(sample)
            genes.add(gene)
            drugs.add(drug)
            sample_mut_drug.add((sample, mut, drug))
            genomic_drug.add(
                (
                    clean(row.get("CHROMOSOME")),
                    clean(row.get("GENOME_START")),
                    clean(row.get("GENOME_STOP")),
                    clean(row.get("GENOMIC_WT_ALLELE")),
                    clean(row.get("GENOMIC_MUT_ALLELE")),
                    drug,
                )
            )
            if sample.upper().replace("-", "") == "HCC1395":
                hcc1395_rows += 1
    return {
        "rows": rows,
        "exact_duplicate_rows": rows - len(full_seen),
        "unique_sample_mutation_drug": len(sample_mut_drug),
        "unique_genomic_allele_drug": len(genomic_drug),
        "transcript_or_sample_expansion_ratio_vs_genomic_drug": round(rows / max(1, len(genomic_drug)), 4),
        "unique_genes": len(genes),
        "unique_drugs": len(drugs),
        "unique_samples": len(samples),
        "hcc1395_named_rows": hcc1395_rows,
    }


def topology_regions() -> tuple[dict[str, list[dict]], dict[str, str]]:
    payload = json.loads(TOPOLOGY.read_text())
    source_by_sample = {
        sample["sample"]: sample["source"]
        for sample in payload["samples"]
        if sample["sample"] in SAMPLES
    }
    seen: dict[str, dict[str, dict]] = {sample: {} for sample in SAMPLES}
    with REGION_TSV.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = row["sample"]
            if name not in seen:
                continue
            region = row["region"]
            if region in seen[name]:
                raise RuntimeError(f"duplicate region key in region TSV: {name} {region}")
            chrom, coords = region.split(":", 1)
            start, end = [int(value) for value in coords.split("-", 1)]
            seen[name][region] = {
                "sample": name,
                "region": region,
                "chrom": chrom,
                "start": start,
                "end": end,
                "span_bp": end - start + 1,
                "n_primary_hp_units": int(row["primary_HP_units"]),
                "primary_hp_families": row["primary_HP_families"],
                "structural_class": row["structural_class"],
            }
    regions_by_sample = {
        sample: sorted(rows.values(), key=lambda row: (row["chrom"], row["start"], row["end"]))
        for sample, rows in seen.items()
    }
    return regions_by_sample, source_by_sample


def make_interval_index(genes: list[dict]) -> tuple[dict, dict]:
    body: dict[tuple[str, int], set[int]] = defaultdict(set)
    promoter: dict[tuple[str, int], set[int]] = defaultdict(set)
    for idx, gene in enumerate(genes):
        for bin_id in range((gene["start"] - 1) // BIN_SIZE, (gene["end"] - 1) // BIN_SIZE + 1):
            body[(gene["chrom"], bin_id)].add(idx)
        tss = gene["start"] if gene["strand"] == "+" else gene["end"]
        pstart = max(1, tss - 2000)
        pend = tss + 2000
        for bin_id in range((pstart - 1) // BIN_SIZE, (pend - 1) // BIN_SIZE + 1):
            promoter[(gene["chrom"], bin_id)].add(idx)
    return body, promoter


def candidate_indices(index: dict, chrom: str, start: int, end: int) -> set[int]:
    result: set[int] = set()
    for bin_id in range((start - 1) // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
        result.update(index.get((chrom, bin_id), ()))
    return result


def get_annotation(gene: dict, cgc_h: dict, cgc_s: dict, dg_h: dict, dg_s: dict) -> tuple[dict | None, dict | None, str, str]:
    cancer = cgc_h.get(gene["hgnc"]) if gene["hgnc"] else None
    cancer_method = "hgnc" if cancer else ""
    if not cancer and gene["symbol"] in cgc_s:
        cancer = cgc_s[gene["symbol"]]
        cancer_method = "symbol_fallback"
    drug = dg_h.get(gene["hgnc"]) if gene["hgnc"] else None
    drug_method = "hgnc" if drug else ""
    if not drug and gene["symbol"] in dg_s:
        drug = dg_s[gene["symbol"]]
        drug_method = "symbol_fallback"
    return cancer, drug, cancer_method, drug_method


def join_regions(
    regions_by_sample: dict[str, list[dict]],
    genes: list[dict],
    cgc_h: dict,
    cgc_s: dict,
    dg_h: dict,
    dg_s: dict,
) -> tuple[list[dict], dict]:
    body_index, promoter_index = make_interval_index(genes)
    long_rows: list[dict] = []
    summary: dict[str, dict] = {}
    for sample, regions in regions_by_sample.items():
        region_gene_counts = []
        region_drug_counts = []
        body_regions = promoter_regions = cancer_regions = drug_regions = approved_anti_regions = 0
        body_pairs = promoter_pairs = cancer_pairs = drug_pairs = approved_anti_pairs = 0
        exploded_drug_rows = 0
        gene_lengths = []
        symbol_fallback_cancer = symbol_fallback_drug = 0
        for region in regions:
            body_ids = {
                idx
                for idx in candidate_indices(body_index, region["chrom"], region["start"], region["end"])
                if genes[idx]["start"] <= region["end"] and genes[idx]["end"] >= region["start"]
            }
            promoter_ids = set()
            for idx in candidate_indices(promoter_index, region["chrom"], region["start"], region["end"]):
                gene = genes[idx]
                tss = gene["start"] if gene["strand"] == "+" else gene["end"]
                if max(1, tss - 2000) <= region["end"] and (tss + 2000) >= region["start"]:
                    promoter_ids.add(idx)
            if body_ids:
                body_regions += 1
            if promoter_ids:
                promoter_regions += 1
            body_pairs += len(body_ids)
            promoter_pairs += len(promoter_ids)
            region_gene_counts.append(len(body_ids | promoter_ids))
            this_cancer = this_drug = this_approved_anti = False
            this_drug_count = 0
            for idx in sorted(body_ids | promoter_ids, key=lambda value: (genes[value]["start"], genes[value]["symbol"])):
                gene = genes[idx]
                cancer, drug, cancer_method, drug_method = get_annotation(gene, cgc_h, cgc_s, dg_h, dg_s)
                if cancer:
                    this_cancer = True
                    cancer_pairs += 1
                    if cancer_method == "symbol_fallback":
                        symbol_fallback_cancer += 1
                if drug:
                    this_drug = True
                    drug_pairs += 1
                    if drug_method == "symbol_fallback":
                        symbol_fallback_drug += 1
                approved_anti = drug["approved_antineoplastic_drugs"] if drug else {}
                approved = drug["approved_drugs"] if drug else {}
                if approved_anti:
                    this_approved_anti = True
                    approved_anti_pairs += 1
                n_drugs = len(drug["drugs"]) if drug else 0
                exploded_drug_rows += n_drugs
                this_drug_count += n_drugs
                gene_lengths.append(gene["end"] - gene["start"] + 1)
                approved_anti_names = sorted(
                    {payload["claim_name"] for payload in approved_anti.values() if payload["claim_name"]}
                )
                long_rows.append(
                    {
                        **region,
                        "gene_id": gene["gene_id"],
                        "gene_symbol": gene["symbol"],
                        "hgnc_id": gene["hgnc"],
                        "gene_type": gene["gene_type"],
                        "gene_start": gene["start"],
                        "gene_end": gene["end"],
                        "gene_span_bp": gene["end"] - gene["start"] + 1,
                        "body_overlap": int(idx in body_ids),
                        "promoter_overlap_tss_2kb": int(idx in promoter_ids),
                        "cgc": int(bool(cancer)),
                        "cgc_join_method": cancer_method,
                        "cgc_role": cancer["role"] if cancer else "",
                        "cgc_tier": cancer["tier"] if cancer else "",
                        "cgc_tumour": cancer["tumour"] if cancer else "",
                        "dgidb": int(bool(drug)),
                        "dgidb_join_method": drug_method,
                        "dgidb_interaction_rows": drug["rows"] if drug else 0,
                        "dgidb_unique_drugs": n_drugs,
                        "dgidb_unique_sources": len(drug["sources"]) if drug else 0,
                        "approved_drug_count": len(approved),
                        "approved_antineoplastic_drug_count": len(approved_anti),
                        "approved_antineoplastic_drug_claim_names": "|".join(approved_anti_names),
                    }
                )
            if this_cancer:
                cancer_regions += 1
            if this_drug:
                drug_regions += 1
            if this_approved_anti:
                approved_anti_regions += 1
            region_drug_counts.append(this_drug_count)
        n_regions = len(regions)
        complete_regions = sum(region["structural_class"] != "incomplete" for region in regions)
        summary[sample] = {
            "regions_before_join": n_regions,
            "complete_regions_before_join": complete_regions,
            "incomplete_regions_before_join": n_regions - complete_regions,
            "region_gene_rows_after_join": sum(region_gene_counts),
            "distinct_regions_after_any_gene_join": sum(1 for value in region_gene_counts if value > 0),
            "regions_without_gene_or_promoter": sum(1 for value in region_gene_counts if value == 0),
            "body_overlap_region_gene_rows": body_pairs,
            "distinct_regions_with_body_gene": body_regions,
            "promoter_overlap_region_gene_rows": promoter_pairs,
            "distinct_regions_with_promoter": promoter_regions,
            "cgc_region_gene_rows": cancer_pairs,
            "distinct_regions_with_cgc_gene": cancer_regions,
            "cgc_region_coverage_pct": pct(cancer_regions, n_regions),
            "dgidb_region_gene_rows": drug_pairs,
            "distinct_regions_with_dgidb_gene": drug_regions,
            "dgidb_region_coverage_pct": pct(drug_regions, n_regions),
            "approved_antineoplastic_region_gene_rows": approved_anti_pairs,
            "distinct_regions_with_approved_antineoplastic_gene": approved_anti_regions,
            "approved_antineoplastic_region_coverage_pct": pct(approved_anti_regions, n_regions),
            "region_gene_join_blowup": round(sum(region_gene_counts) / max(1, n_regions), 4),
            "region_gene_drug_exploded_rows": exploded_drug_rows,
            "region_gene_drug_blowup_vs_regions": round(exploded_drug_rows / max(1, n_regions), 4),
            "region_gene_drug_blowup_vs_region_gene": round(exploded_drug_rows / max(1, sum(region_gene_counts)), 4),
            "max_genes_per_region": max(region_gene_counts, default=0),
            "p95_genes_per_region": round(quantile(region_gene_counts, 0.95), 4),
            "max_drugs_per_region_after_explode": max(region_drug_counts, default=0),
            "p95_drugs_per_region_after_explode": round(quantile(region_drug_counts, 0.95), 4),
            "median_joined_gene_length_bp": round(statistics.median(gene_lengths), 2) if gene_lengths else 0,
            "cgc_symbol_fallback_rows": symbol_fallback_cancer,
            "dgidb_symbol_fallback_rows": symbol_fallback_drug,
        }
    return long_rows, summary


def build_region_flags(regions_by_sample: dict[str, list[dict]], long_rows: list[dict]) -> dict[tuple[str, str], dict]:
    """Collapse the one-to-many join to one record per region.

    Gene-body overlap is the primary annotation.  TSS±2 kb is retained as an
    explicitly secondary promoter view, per the parent analysis contract.
    """
    flags: dict[tuple[str, str], dict] = {}
    set_fields = (
        "body_gene_symbols",
        "promoter_gene_symbols",
        "cgc_body_symbols",
        "cgc_promoter_symbols",
        "dgidb_interaction_body_symbols",
        "dgidb_interaction_promoter_symbols",
        "dgidb_approved_body_symbols",
        "dgidb_approved_promoter_symbols",
        "dgidb_approved_antineoplastic_body_symbols",
        "dgidb_approved_antineoplastic_promoter_symbols",
        "approved_antineoplastic_body_drug_claim_names",
    )
    for sample, regions in regions_by_sample.items():
        for region in regions:
            record = {field: set() for field in set_fields}
            record.update({
                "sample": sample,
                "region": region["region"],
                "chrom": region["chrom"],
                "start": region["start"],
                "end": region["end"],
            })
            flags[(sample, region["region"])] = record
    for row in long_rows:
        record = flags[(row["sample"], row["region"])]
        symbol = row["gene_symbol"]
        body = bool(row["body_overlap"])
        promoter = bool(row["promoter_overlap_tss_2kb"])
        if body:
            record["body_gene_symbols"].add(symbol)
        if promoter:
            record["promoter_gene_symbols"].add(symbol)
        if row["cgc"] and body:
            record["cgc_body_symbols"].add(symbol)
        if row["cgc"] and promoter:
            record["cgc_promoter_symbols"].add(symbol)
        if row["dgidb"] and body:
            record["dgidb_interaction_body_symbols"].add(symbol)
        if row["dgidb"] and promoter:
            record["dgidb_interaction_promoter_symbols"].add(symbol)
        if int(row["approved_drug_count"]) > 0 and body:
            record["dgidb_approved_body_symbols"].add(symbol)
        if int(row["approved_drug_count"]) > 0 and promoter:
            record["dgidb_approved_promoter_symbols"].add(symbol)
        if int(row["approved_antineoplastic_drug_count"]) > 0 and body:
            record["dgidb_approved_antineoplastic_body_symbols"].add(symbol)
            record["approved_antineoplastic_body_drug_claim_names"].update(
                value for value in row["approved_antineoplastic_drug_claim_names"].split("|") if value
            )
        if int(row["approved_antineoplastic_drug_count"]) > 0 and promoter:
            record["dgidb_approved_antineoplastic_promoter_symbols"].add(symbol)
    return flags


def enrich_body_primary_join_summary(summary: dict, flags: dict[tuple[str, str], dict]) -> None:
    feature_to_field = {
        "body_gene": "body_gene_symbols",
        "cgc_body": "cgc_body_symbols",
        "dgidb_interaction_body": "dgidb_interaction_body_symbols",
        "dgidb_approved_body": "dgidb_approved_body_symbols",
        "dgidb_approved_antineoplastic_body": "dgidb_approved_antineoplastic_body_symbols",
        "promoter_secondary": "promoter_gene_symbols",
    }
    for sample in SAMPLES:
        records = [record for (name, _), record in flags.items() if name == sample]
        for label, field in feature_to_field.items():
            n = sum(bool(record[field]) for record in records)
            summary[sample][f"distinct_regions_with_{label}"] = n
            summary[sample][f"{label}_region_coverage_pct"] = pct(n, len(records))


def cohens_kappa(rows: list[dict]) -> float | None:
    if not rows:
        return None
    cats = sorted({row["category_a"] for row in rows} | {row["category_b"] for row in rows})
    n = len(rows)
    po = sum(row["category_agree"] for row in rows) / n
    ca = Counter(row["category_a"] for row in rows)
    cb = Counter(row["category_b"] for row in rows)
    pe = sum((ca[cat] / n) * (cb[cat] / n) for cat in cats)
    return (po - pe) / (1 - pe) if pe < 1 else None


def exact_complete_pair_annotation(
    flags: dict[tuple[str, str], dict]
) -> tuple[list[dict], list[dict], dict]:
    """Attach body-primary annotations to the canonical exact complete-both pair."""
    pair_rows: list[dict] = []
    with PAIR_MATCHES.open() as handle:
        for raw in csv.DictReader(handle, delimiter="\t"):
            if raw["scenario"] != "exact_coordinate":
                continue
            if raw["complete_a"] != "True" or raw["complete_b"] != "True":
                continue
            if raw["region_a"] != raw["region_b"]:
                raise RuntimeError("exact-coordinate match has different region keys")
            annotation = flags[("HCC1395", raw["region_a"])]
            row = dict(raw)
            row["category_agree"] = raw["category_agree"] == "True"
            for field, value in annotation.items():
                if isinstance(value, set):
                    row[field] = "|".join(sorted(value))
            for field in (
                "body_gene_symbols",
                "promoter_gene_symbols",
                "cgc_body_symbols",
                "cgc_promoter_symbols",
                "dgidb_interaction_body_symbols",
                "dgidb_interaction_promoter_symbols",
                "dgidb_approved_body_symbols",
                "dgidb_approved_promoter_symbols",
                "dgidb_approved_antineoplastic_body_symbols",
                "dgidb_approved_antineoplastic_promoter_symbols",
                "approved_antineoplastic_body_drug_claim_names",
            ):
                row[field.replace("_symbols", "_any").replace("_drug_claim_names", "_any")] = bool(row.get(field))
            pair_rows.append(row)
    if len(pair_rows) != 5720:
        raise RuntimeError(f"exact complete-both pair expected 5720 rows, observed {len(pair_rows)}")

    features = [
        ("primary_gene_body", "body_gene_any"),
        ("primary_COSMIC_CGC_body", "cgc_body_any"),
        ("primary_DGIdb_interaction_body", "dgidb_interaction_body_any"),
        ("primary_DGIdb_approved_body", "dgidb_approved_body_any"),
        ("primary_DGIdb_approved_antineoplastic_body", "dgidb_approved_antineoplastic_body_any"),
        ("secondary_promoter_TSS_2kb", "promoter_gene_any"),
        ("secondary_promoter_CGC", "cgc_promoter_any"),
        ("secondary_promoter_DGIdb_interaction", "dgidb_interaction_promoter_any"),
        ("secondary_promoter_DGIdb_approved", "dgidb_approved_promoter_any"),
        (
            "secondary_promoter_DGIdb_approved_antineoplastic",
            "dgidb_approved_antineoplastic_promoter_any",
        ),
    ]

    def summarize(label: str, stratum: str, rows: list[dict]) -> dict:
        agree = sum(row["category_agree"] for row in rows)
        return {
            "feature": label,
            "stratum": stratum,
            "n_exact_complete_pairs": len(rows),
            "category_agree_n": agree,
            "category_agreement_pct": pct(agree, len(rows)),
            "cohens_kappa": round(cohens_kappa(rows), 6) if cohens_kappa(rows) is not None else "",
            "category_counts_a": json.dumps(Counter(row["category_a"] for row in rows), ensure_ascii=False, sort_keys=True),
            "category_counts_b": json.dumps(Counter(row["category_b"] for row in rows), ensure_ascii=False, sort_keys=True),
            "difference_present_minus_absent_pp": "",
            "interpretation": "descriptive only; same coordinate implies the annotation flag itself is shared",
        }

    baseline = summarize("ALL", "all", pair_rows)
    stratified = [baseline]
    for label, field in features:
        present = [row for row in pair_rows if row[field]]
        absent = [row for row in pair_rows if not row[field]]
        present_row = summarize(label, "present", present)
        absent_row = summarize(label, "absent", absent)
        delta = round(present_row["category_agreement_pct"] - absent_row["category_agreement_pct"], 4)
        present_row["difference_present_minus_absent_pp"] = delta
        absent_row["difference_present_minus_absent_pp"] = delta
        stratified.extend((present_row, absent_row))
    checks = {
        "exact_complete_pair_rows": len(pair_rows),
        "baseline_agree_n": sum(row["category_agree"] for row in pair_rows),
        "baseline_agreement_pct": pct(sum(row["category_agree"] for row in pair_rows), len(pair_rows)),
        "baseline_kappa": round(cohens_kappa(pair_rows), 6) if cohens_kappa(pair_rows) is not None else None,
        "all_strata_conserve_5720": all(
            sum(row["n_exact_complete_pairs"] for row in stratified if row["feature"] == label) == 5720
            for label, _ in features
        ),
    }
    return pair_rows, stratified, checks


def clp_hcc1395_sensitivity(pair_rows: list[dict]) -> tuple[list[dict], list[dict], dict]:
    """Extract HCC1395 CLP v104 alleles and overlap the exact complete pair.

    The CLP catalogue is a sample-specific sensitivity annotation only.  It is
    not a topology label and may represent a different assay, passage, or
    library from either ONT processing dataset.
    """
    region_bins: dict[tuple[str, int], set[int]] = defaultdict(set)
    for idx, row in enumerate(pair_rows):
        start = int(row["region_a"].split(":", 1)[1].split("-", 1)[0])
        end = int(row["region_a"].rsplit("-", 1)[1])
        for bin_id in range((start - 1) // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
            region_bins[(row["chrom"], bin_id)].add(idx)
        row["clp_all_variant_count"] = 0
        row["clp_confirmed_somatic_variant_count"] = 0
        row["clp_reported_elsewhere_somatic_variant_count"] = 0
        row["clp_unknown_origin_variant_count"] = 0

    variants: dict[tuple, dict] = {}
    source_rows = 0
    hcc_rows = 0
    with gzip.open(CLP_MUTANT, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source_rows += 1
            if clean(row.get("SAMPLE_NAME")).upper() != "HCC1395":
                continue
            chrom_raw = clean(row.get("CHROMOSOME"))
            if not chrom_raw.isdigit() or not 1 <= int(chrom_raw) <= 22:
                continue
            hcc_rows += 1
            chrom = f"chr{int(chrom_raw)}"
            start = int(clean(row.get("GENOME_START")))
            end = int(clean(row.get("GENOME_STOP")))
            wt = clean(row.get("GENOMIC_WT_ALLELE"))
            alt = clean(row.get("GENOMIC_MUT_ALLELE"))
            key = (chrom, start, end, wt, alt)
            record = variants.setdefault(
                key,
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "ref": wt,
                    "alt": alt,
                    "sample_name": "HCC1395",
                    "cosmic_sample_ids": set(),
                    "cosmic_genomic_mutation_ids": set(),
                    "gene_symbols": set(),
                    "somatic_statuses": set(),
                    "transcript_rows": 0,
                    "overlapping_exact_complete_regions": set(),
                },
            )
            record["cosmic_sample_ids"].add(clean(row.get("COSMIC_SAMPLE_ID")))
            record["cosmic_genomic_mutation_ids"].add(clean(row.get("GENOMIC_MUTATION_ID")))
            record["gene_symbols"].add(clean(row.get("GENE_SYMBOL")).upper())
            record["somatic_statuses"].add(clean(row.get("MUTATION_SOMATIC_STATUS")) or "(missing)")
            record["transcript_rows"] += 1

    region_variant_keys: dict[int, set[tuple]] = defaultdict(set)
    for key, variant in variants.items():
        candidates: set[int] = set()
        for bin_id in range((variant["start"] - 1) // BIN_SIZE, (variant["end"] - 1) // BIN_SIZE + 1):
            candidates.update(region_bins.get((variant["chrom"], bin_id), ()))
        for idx in candidates:
            region = pair_rows[idx]
            start = int(region["region_a"].split(":", 1)[1].split("-", 1)[0])
            end = int(region["region_a"].rsplit("-", 1)[1])
            if variant["start"] <= end and variant["end"] >= start:
                region_variant_keys[idx].add(key)
                variant["overlapping_exact_complete_regions"].add(region["region_a"])

    confirmed_label = "Confirmed somatic variant"
    elsewhere_label = "Reported in another cancer sample as somatic"
    unknown_label = "Variant of unknown origin"
    for idx, keys in region_variant_keys.items():
        pair_rows[idx]["clp_all_variant_count"] = len(keys)
        pair_rows[idx]["clp_confirmed_somatic_variant_count"] = sum(
            confirmed_label in variants[key]["somatic_statuses"] for key in keys
        )
        pair_rows[idx]["clp_reported_elsewhere_somatic_variant_count"] = sum(
            elsewhere_label in variants[key]["somatic_statuses"] for key in keys
        )
        pair_rows[idx]["clp_unknown_origin_variant_count"] = sum(
            unknown_label in variants[key]["somatic_statuses"] for key in keys
        )

    extracted = []
    for key in sorted(variants):
        variant = variants[key]
        extracted.append(
            {
                "chrom": variant["chrom"],
                "start": variant["start"],
                "end": variant["end"],
                "ref": variant["ref"],
                "alt": variant["alt"],
                "sample_name": variant["sample_name"],
                "cosmic_sample_ids": "|".join(sorted(value for value in variant["cosmic_sample_ids"] if value)),
                "cosmic_genomic_mutation_ids": "|".join(
                    sorted(value for value in variant["cosmic_genomic_mutation_ids"] if value)
                ),
                "gene_symbols": "|".join(sorted(value for value in variant["gene_symbols"] if value)),
                "somatic_statuses": "|".join(sorted(variant["somatic_statuses"])),
                "confirmed_somatic": confirmed_label in variant["somatic_statuses"],
                "reported_elsewhere_somatic": elsewhere_label in variant["somatic_statuses"],
                "unknown_origin": unknown_label in variant["somatic_statuses"],
                "transcript_rows": variant["transcript_rows"],
                "overlapping_exact_complete_regions": "|".join(
                    sorted(variant["overlapping_exact_complete_regions"])
                ),
            }
        )

    def summarize(label: str, rows: list[dict]) -> dict:
        agree = sum(row["category_agree"] for row in rows)
        return {
            "feature": label,
            "n_exact_complete_pairs": len(rows),
            "category_agree_n": agree,
            "category_agreement_pct": pct(agree, len(rows)),
            "cohens_kappa": round(cohens_kappa(rows), 6) if cohens_kappa(rows) is not None else "",
            "category_counts_a": json.dumps(Counter(row["category_a"] for row in rows), sort_keys=True),
            "category_counts_b": json.dumps(Counter(row["category_b"] for row in rows), sort_keys=True),
            "difference_present_minus_absent_pp": "",
            "claim_ceiling": "sample-specific catalogue sensitivity; not topology/clone truth",
        }

    strata = []
    for label, field in (
        ("CLP_all_genomic_allele_present", "clp_all_variant_count"),
        ("CLP_confirmed_somatic_allele_present", "clp_confirmed_somatic_variant_count"),
        ("CLP_unknown_origin_allele_present", "clp_unknown_origin_variant_count"),
    ):
        present = [row for row in pair_rows if int(row[field]) > 0]
        absent = [row for row in pair_rows if int(row[field]) == 0]
        p_row = summarize(label, present)
        a_row = summarize(label.replace("_present", "_absent"), absent)
        delta = round(p_row["category_agreement_pct"] - a_row["category_agreement_pct"], 4)
        p_row["difference_present_minus_absent_pp"] = delta
        a_row["difference_present_minus_absent_pp"] = delta
        strata.extend((p_row, a_row))

    status_counts = Counter()
    for variant in variants.values():
        for status in variant["somatic_statuses"]:
            status_counts[status] += 1
    profile = {
        "source_rows_all_samples": source_rows,
        "hcc1395_chr1_22_transcript_rows": hcc_rows,
        "hcc1395_chr1_22_unique_genomic_alleles": len(variants),
        "cosmic_sample_ids": sorted(
            {value for variant in variants.values() for value in variant["cosmic_sample_ids"] if value}
        ),
        "unique_alleles_by_somatic_status_membership": dict(status_counts),
        "confirmed_somatic_unique_alleles": sum(
            confirmed_label in variant["somatic_statuses"] for variant in variants.values()
        ),
        "unknown_origin_unique_alleles": sum(
            unknown_label in variant["somatic_statuses"] for variant in variants.values()
        ),
        "reported_elsewhere_somatic_unique_alleles": sum(
            elsewhere_label in variant["somatic_statuses"] for variant in variants.values()
        ),
        "all_unique_alleles_overlapping_exact_complete_pair": sum(
            bool(variant["overlapping_exact_complete_regions"]) for variant in variants.values()
        ),
        "confirmed_somatic_unique_alleles_overlapping_exact_complete_pair": sum(
            confirmed_label in variant["somatic_statuses"]
            and bool(variant["overlapping_exact_complete_regions"])
            for variant in variants.values()
        ),
        "exact_complete_regions_with_any_clp_allele": sum(
            int(row["clp_all_variant_count"]) > 0 for row in pair_rows
        ),
        "exact_complete_regions_with_confirmed_somatic_clp_allele": sum(
            int(row["clp_confirmed_somatic_variant_count"]) > 0 for row in pair_rows
        ),
        "exact_complete_regions_with_unknown_origin_clp_allele": sum(
            int(row["clp_unknown_origin_variant_count"]) > 0 for row in pair_rows
        ),
        "strata_conserve_5720": all(
            strata[idx]["n_exact_complete_pairs"] + strata[idx + 1]["n_exact_complete_pairs"] == 5720
            for idx in range(0, len(strata), 2)
        ),
    }
    return extracted, strata, profile


def exact_pair_gene_summary(
    pair_rows: list[dict], cgc_s: dict[str, dict], dg_s: dict[str, dict]
) -> list[dict]:
    notable = {"BRCA2", "TBC1D16", "ERBB2", "MYC"}
    selected = set(notable)
    for row in pair_rows:
        selected.update(value for value in row["cgc_body_symbols"].split("|") if value)
        selected.update(
            value for value in row["dgidb_approved_antineoplastic_body_symbols"].split("|") if value
        )
    output = []
    for symbol in sorted(selected):
        rows = [
            row for row in pair_rows if symbol in set(value for value in row["body_gene_symbols"].split("|") if value)
        ]
        cancer = cgc_s.get(symbol)
        drug = dg_s.get(symbol)
        approved = drug["approved_drugs"] if drug else {}
        approved_anti = drug["approved_antineoplastic_drugs"] if drug else {}
        approved_names = sorted(
            {value["claim_name"] for value in approved.values() if value.get("claim_name")}
        )
        approved_anti_names = sorted(
            {value["claim_name"] for value in approved_anti.values() if value.get("claim_name")}
        )
        output.append(
            {
                "gene_symbol": symbol,
                "user_notable_gene": symbol in notable,
                "n_exact_complete_regions_with_body_overlap": len(rows),
                "category_agree_n": sum(row["category_agree"] for row in rows),
                "category_agreement_pct": pct(sum(row["category_agree"] for row in rows), len(rows)),
                "cohens_kappa": round(cohens_kappa(rows), 6) if cohens_kappa(rows) is not None else "",
                "cgc_v104": bool(cancer),
                "cgc_role": cancer["role"] if cancer else "",
                "cgc_tier": cancer["tier"] if cancer else "",
                "dgidb_interaction_rows": drug["rows"] if drug else 0,
                "dgidb_unique_drugs": len(drug["drugs"]) if drug else 0,
                "dgidb_approved_unique_drugs": len(approved),
                "dgidb_approved_antineoplastic_unique_drugs": len(approved_anti),
                "approved_drug_claim_names": "|".join(approved_names),
                "approved_antineoplastic_drug_claim_names": "|".join(approved_anti_names),
                "region_keys": "|".join(row["region_a"] for row in rows),
                "claim_ceiling": "descriptive annotation-stratified technical agreement; not clinical actionability",
            }
        )
    return output


def notable_gene_locus_sensitivity(long_rows: list[dict]) -> list[dict]:
    notable = ("BRCA2", "TBC1D16", "ERBB2", "MYC")
    gene_regions: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in long_rows:
        if row["gene_symbol"] in notable and int(row["body_overlap"]):
            gene_regions[(row["sample"], row["gene_symbol"])].add(row["region"])

    coarse = {}
    with COARSE_REGIONS.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] in SAMPLES:
                coarse[(row["sample"], row["region"])] = row

    matches = []
    with PAIR_MATCHES.open() as handle:
        matches = list(csv.DictReader(handle, delimiter="\t"))
    scenario_rank = {
        "exact_coordinate": 0,
        "reciprocal_overlap_0.80": 1,
        "reciprocal_overlap_0.50": 2,
        "endpoint_anchor_5kb": 3,
    }
    output = []
    for gene in notable:
        a_regions = sorted(gene_regions.get(("HCC1395", gene), set()))
        b_regions = sorted(gene_regions.get(("HCC1395_DORADO", gene), set()))
        linked = [
            row
            for row in matches
            if row["region_a"] in a_regions and row["region_b"] in b_regions
        ]
        linked.sort(key=lambda row: (scenario_rank.get(row["scenario"], 99), -float(row["reciprocal_min"])))
        best = linked[0] if linked else {}
        a_meta = [coarse[("HCC1395", region)] for region in a_regions]
        b_meta = [coarse[("HCC1395_DORADO", region)] for region in b_regions]
        output.append(
            {
                "gene_symbol": gene,
                "hcc1395_body_regions": "|".join(a_regions),
                "hcc1395_structural_classes": "|".join(row["structural_class"] for row in a_meta),
                "hcc1395_coarse_categories": "|".join(row["coarse_category"] for row in a_meta),
                "dorado_body_regions": "|".join(b_regions),
                "dorado_structural_classes": "|".join(row["structural_class"] for row in b_meta),
                "dorado_coarse_categories": "|".join(row["coarse_category"] for row in b_meta),
                "best_pair_scenario": best.get("scenario", "unmatched"),
                "best_pair_reciprocal_overlap": best.get("reciprocal_min", ""),
                "best_pair_complete_a": best.get("complete_a", ""),
                "best_pair_complete_b": best.get("complete_b", ""),
                "best_pair_category_a": best.get("category_a", ""),
                "best_pair_category_b": best.get("category_b", ""),
                "best_pair_category_agree": best.get("category_agree", ""),
                "interpretation": (
                    "coarse topology locus sensitivity; exact tree/shape agreement and clinical actionability are not implied"
                ),
            }
        )
    return output


def profile_old_annotations(regions_by_sample: dict[str, list[dict]]) -> dict:
    result = {}
    for sample, path in OLD_REGION_ANN.items():
        payload = json.loads(path.read_text())
        old_keys = set(payload.get("regions", {}))
        new_keys = {row["region"] for row in regions_by_sample[sample]}
        result[sample] = {
            "path": str(path),
            "reported_n_regions": payload.get("summary", {}).get("n_regions"),
            "stored_region_keys": len(old_keys),
            "latest_region_keys": len(new_keys),
            "exact_key_overlap": len(old_keys & new_keys),
            "latest_coverage_pct": pct(len(old_keys & new_keys), len(new_keys)),
            "old_keys_not_in_latest": len(old_keys - new_keys),
            "latest_keys_missing_old_annotation": len(new_keys - old_keys),
            "sha256": sha256(path),
        }
    return result


def profile_wakhan(genes: list[dict], cgc_symbols: dict[str, dict]) -> dict:
    by_symbol: dict[str, list[dict]] = defaultdict(list)
    for gene in genes:
        by_symbol[gene["symbol"]].append(gene)
    rows = []
    with WAKHAN.open() as handle:
        for line in handle:
            chrom, start, end, symbol = line.rstrip("\n").split("\t")
            rows.append((chrom, int(start), int(end), symbol.upper()))
    symbols = Counter(row[3] for row in rows)
    exact_gene_coords = 0
    gencode_overlap_same_symbol = 0
    ratios = []
    for chrom, start, end, symbol in rows:
        for gene in by_symbol.get(symbol, []):
            if gene["chrom"] != chrom:
                continue
            if gene["start"] == start and gene["end"] == end:
                exact_gene_coords += 1
            if gene["start"] <= end and gene["end"] >= start:
                gencode_overlap_same_symbol += 1
                ratios.append((end - start + 1) / max(1, gene["end"] - gene["start"] + 1))
                break
    return {
        "rows": len(rows),
        "unique_symbols": len(symbols),
        "duplicate_symbol_rows": sum(value - 1 for value in symbols.values() if value > 1),
        "cgc_v104_member_symbols": sum(1 for symbol in symbols if symbol in cgc_symbols),
        "exact_gencode_v46_coordinate_rows": exact_gene_coords,
        "same_symbol_body_contained_or_overlapped_rows": gencode_overlap_same_symbol,
        "median_interval_bp": statistics.median(end - start + 1 for _, start, end, _ in rows),
        "max_interval_bp": max(end - start + 1 for _, start, end, _ in rows),
        "median_interval_to_gene_length_ratio": round(statistics.median(ratios), 4) if ratios else 0,
    }


def source_inventory(profiles: dict, source_by_sample: dict[str, str]) -> list[dict]:
    rows = [
        {
            "source_id": "gencode_v46_basic",
            "absolute_path": str(GTF),
            "version_or_snapshot": "GENCODE v46 / Ensembl 112; annotation date 2024-03-26",
            "genome_build": "GRCh38",
            "grain": "one GTF gene feature per Ensembl gene ID (profile excludes non-gene features)",
            "row_count": profiles["gencode"]["gene_rows"],
            "unique_key_count": profiles["gencode"]["unique_gene_ids"],
            "duplicate_key_rows": profiles["gencode"]["gene_rows"] - profiles["gencode"]["unique_gene_ids"],
            "join_key": "HGNC ID preferred; Ensembl gene ID; symbol only fallback",
            "license_or_terms": "license file not bundled in local annotation directory; redistribution status not re-verified",
            "coverage_caveat": "basic annotation; symbol is non-unique; coordinates are 1-based inclusive",
            "decision": "USE for region→gene coordinates",
        },
        {
            "source_id": "cosmic_cgc_v104",
            "absolute_path": str(CGC),
            "version_or_snapshot": "COSMIC Cancer Gene Census v104; local file mtime 2026-02-20",
            "genome_build": "GRCh38",
            "grain": "one curated cancer-gene record",
            "row_count": profiles["cgc"]["rows"],
            "unique_key_count": profiles["cgc"]["unique_symbols"],
            "duplicate_key_rows": profiles["cgc"]["duplicate_symbol_rows"],
            "join_key": "COSMIC_GENE_ID→Cosmic_Genes.HGNC_ID preferred; symbol fallback",
            "license_or_terms": "COSMIC README contains schema only; redistribution terms not bundled/re-verified",
            "coverage_caveat": "curated census, not all cancer-relevant genes; TBC1D16 is not a CGC v104 member",
            "decision": "USE as cancer-gene annotation, not topology truth",
        },
        {
            "source_id": "cosmic_genes_v104",
            "absolute_path": str(COSMIC_GENES),
            "version_or_snapshot": "COSMIC Genes v104; local file mtime 2026-02-20",
            "genome_build": "GRCh38",
            "grain": "one COSMIC gene mapping record",
            "row_count": profiles["cosmic_genes"]["rows"],
            "unique_key_count": profiles["cosmic_genes"]["unique_cosmic_gene_id"],
            "duplicate_key_rows": profiles["cosmic_genes"]["duplicate_cosmic_gene_id_rows"],
            "join_key": "COSMIC_GENE_ID; HGNC_ID/Ensembl are bridge keys",
            "license_or_terms": "COSMIC README contains schema only; redistribution terms not bundled/re-verified",
            "coverage_caveat": "mapping bridge, not cancer-role evidence",
            "decision": "USE to avoid symbol-only CGC join",
        },
        {
            "source_id": "dgidb_interactions_local_snapshot",
            "absolute_path": str(DGIDB),
            "version_or_snapshot": "overall DGIdb release absent; local download mtime 2026-06-29; embedded source versions vary",
            "genome_build": "N/A (gene/drug concepts)",
            "grain": "gene–drug–source interaction row",
            "row_count": profiles["dgidb"]["rows"],
            "unique_key_count": profiles["dgidb"]["unique_gene_concept_ids"],
            "duplicate_key_rows": profiles["dgidb"]["duplicate_composite_rows"],
            "join_key": "gene_concept_id (HGNC) preferred; gene_name fallback",
            "license_or_terms": "no README/license/provenance receipt beside local TSV; redistribution/currentness unverified",
            "coverage_caveat": "many-to-many; source versions heterogeneous; interaction≠benefit; canonical drug labels contain observed normalization anomalies",
            "decision": "USE only as exploratory interaction annotation; strongest display subset=approved AND anti_neoplastic",
        },
        {
            "source_id": "cosmic_resistance_v104",
            "absolute_path": str(RESISTANCE),
            "version_or_snapshot": "COSMIC Resistance Mutations v104; local file mtime 2026-02-20",
            "genome_build": "GRCh38",
            "grain": "sample×mutation×transcript×drug resistance observation",
            "row_count": profiles["resistance"]["rows"],
            "unique_key_count": profiles["resistance"]["unique_genomic_allele_drug"],
            "duplicate_key_rows": profiles["resistance"]["rows"] - profiles["resistance"]["unique_genomic_allele_drug"],
            "join_key": "exact CHROM/START/STOP/WT/ALT + drug; never gene-only",
            "license_or_terms": "COSMIC README contains schema only; redistribution terms not bundled/re-verified",
            "coverage_caveat": "transcript/sample expansion; current coarse topology artifact lacks exact variant allele keys",
            "decision": "BLOCK gene-only join; exact variant join unavailable in this artifact",
        },
        {
            "source_id": "cosmic_clp_genome_screens_mutant_v104",
            "absolute_path": str(CLP_MUTANT),
            "version_or_snapshot": "COSMIC Cell Lines Project Genome Screens Mutant v104; local file mtime 2026-02-20",
            "genome_build": "GRCh38",
            "grain": "all-sample coding mutation transcript row; HCC1395 sensitivity extracted separately",
            "row_count": profiles["clp_hcc1395"]["source_rows_all_samples"],
            "unique_key_count": "not_profiled_globally",
            "duplicate_key_rows": "not_profiled_globally",
            "join_key": "SAMPLE_NAME/COSMIC_SAMPLE_ID + exact genomic allele; transcript rows must be collapsed",
            "license_or_terms": "COSMIC README contains schema only; redistribution terms not bundled/re-verified",
            "coverage_caveat": "different assay/passage possible; coding-screen ascertainment; somatic status includes unknown origin",
            "decision": "USE only for HCC1395 sample-specific sensitivity, never as topology truth",
        },
        {
            "source_id": "cosmic_clp_hcc1395_chr1_22_alleles",
            "absolute_path": str(OUT / "hcc1395_cosmic_clp_v104_chr1_22_alleles.tsv.gz"),
            "version_or_snapshot": "derived from COSMIC CLP v104 in this run",
            "genome_build": "GRCh38 chr1-22",
            "grain": "one HCC1395 exact genomic allele after transcript collapse",
            "row_count": profiles["clp_hcc1395"]["hcc1395_chr1_22_unique_genomic_alleles"],
            "unique_key_count": profiles["clp_hcc1395"]["hcc1395_chr1_22_unique_genomic_alleles"],
            "duplicate_key_rows": 0,
            "join_key": "chrom,start,end,ref,alt",
            "license_or_terms": "derived from COSMIC; redistribution terms not bundled/re-verified",
            "coverage_caveat": (
                f"collapsed from {profiles['clp_hcc1395']['hcc1395_chr1_22_transcript_rows']} HCC1395 transcript rows; "
                "Confirmed somatic and unknown-origin reported separately"
            ),
            "decision": "USE as sample-specific positional sensitivity only",
        },
        {
            "source_id": "wakhan_bundled_cosmic_named_subset",
            "absolute_path": str(WAKHAN),
            "version_or_snapshot": "WAKHAN repo commit 013b52e (0.4.3); file itself has no CGC version/header",
            "genome_build": "appears chr-prefixed GRCh38 but not self-declared in file",
            "grain": "four-column interval×symbol",
            "row_count": profiles["wakhan"]["rows"],
            "unique_key_count": profiles["wakhan"]["unique_symbols"],
            "duplicate_key_rows": profiles["wakhan"]["duplicate_symbol_rows"],
            "join_key": "none safe for this analysis",
            "license_or_terms": "WAKHAN code is MIT; that does not establish license/provenance of bundled COSMIC-named data",
            "coverage_caveat": "99 rows/98 symbols; intervals frequently multi-megabase and do not equal gene bodies; only a subset overlaps CGC v104",
            "decision": "DO NOT USE as canonical cancer-gene truth",
        },
        {
            "source_id": "coarse_topology_all_regions",
            "absolute_path": str(COARSE_REGIONS),
            "version_or_snapshot": "2026-07-12 canonical five-class derivative of historical layered-v2",
            "genome_build": "GRCh38/chr1-22 analysis scope",
            "grain": "one dataset×region row",
            "row_count": profiles["topology_pair"]["coarse_region_rows"],
            "unique_key_count": profiles["topology_pair"]["coarse_region_rows"],
            "duplicate_key_rows": 0,
            "join_key": "sample + region",
            "license_or_terms": "project-generated",
            "coverage_caveat": "historical layered-v2 engineering snapshot; not clean-v3 scientific release",
            "decision": "USE as mother table for annotation join",
        },
        {
            "source_id": "hcc1395_pair_matches",
            "absolute_path": str(PAIR_MATCHES),
            "version_or_snapshot": "2026-07-12 exact/overlap pair derivative",
            "genome_build": "GRCh38/chr1-22 analysis scope",
            "grain": "one scenario×matched-region-pair row",
            "row_count": profiles["topology_pair"]["all_match_rows"],
            "unique_key_count": profiles["topology_pair"]["all_match_rows"],
            "duplicate_key_rows": 0,
            "join_key": "scenario + match_id",
            "license_or_terms": "project-generated",
            "coverage_caveat": "annotation-stratified result is restricted to exact_coordinate + complete_a + complete_b (n=5,720)",
            "decision": "USE as pair mother table",
        },
    ]
    for sample in SAMPLES:
        rows.append(
            {
                "source_id": f"layered_v2_regions_{sample}",
                "absolute_path": source_by_sample[sample],
                "version_or_snapshot": "20260710_232501 layered_reconstruction_v2 historical engineering snapshot",
                "genome_build": "GRCh38/chr1-22 analysis scope",
                "grain": "region×HP family; deduplicated to region before annotation join",
                "row_count": profiles["join"][sample]["regions_before_join"],
                "unique_key_count": profiles["join"][sample]["regions_before_join"],
                "duplicate_key_rows": 0,
                "join_key": "region string chrom:start-end",
                "license_or_terms": "project-generated",
                "coverage_caveat": "historical layered-v2, not clean-v3 scientific release; same biological HCC1395 sample",
                "decision": "USE for this engineering comparison with PARTIAL/NO-GO label",
            }
        )
    for row in rows:
        path = Path(row["absolute_path"])
        if path.exists():
            row["file_bytes"] = path.stat().st_size
            row["file_mtime"] = mtime(path)
            row["sha256"] = sha256(path)
        else:
            row["file_bytes"] = ""
            row["file_mtime"] = ""
            row["sha256"] = ""
    return rows


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None, gz: bool = False) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty TSV: {path}")
    fields = fields or list(rows[0])
    opener = gzip.open if gz else open
    mode = "wt" if gz else "w"
    with opener(path, mode, newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cosmic_genes, cosmic_genes_profile = load_cosmic_genes()
    cgc_h, cgc_s, cgc_profile = load_cgc(cosmic_genes)
    dg_h, dg_s, dgidb_profile = load_dgidb()
    genes, gtf_profile = load_gtf()
    resistance_profile = load_resistance()
    regions_by_sample, source_by_sample = topology_regions()
    long_rows, join_summary = join_regions(regions_by_sample, genes, cgc_h, cgc_s, dg_h, dg_s)
    region_flags = build_region_flags(regions_by_sample, long_rows)
    enrich_body_primary_join_summary(join_summary, region_flags)
    exact_pair_rows, annotation_agreement, pair_checks = exact_complete_pair_annotation(region_flags)
    clp_variants, clp_sensitivity, clp_profile = clp_hcc1395_sensitivity(exact_pair_rows)
    pair_gene_summary = exact_pair_gene_summary(exact_pair_rows, cgc_s, dg_s)
    notable_loci = notable_gene_locus_sensitivity(long_rows)
    write_tsv(OUT / "hcc1395_cosmic_clp_v104_chr1_22_alleles.tsv.gz", clp_variants, gz=True)
    coarse_region_rows = sum(1 for _ in COARSE_REGIONS.open()) - 1
    all_match_rows = 0
    exact_match_rows = 0
    with PAIR_MATCHES.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            all_match_rows += 1
            exact_match_rows += int(row["scenario"] == "exact_coordinate")
    if coarse_region_rows != 47377:
        raise RuntimeError(f"coarse mother table expected 47377 rows, observed {coarse_region_rows}")
    if exact_match_rows != 6252:
        raise RuntimeError(f"exact-coordinate pair expected 6252 rows, observed {exact_match_rows}")
    old_profile = profile_old_annotations(regions_by_sample)
    wakhan_profile = profile_wakhan(genes, cgc_s)
    profiles = {
        "generated_at": datetime.now().astimezone().isoformat(),
        "grain": "source profile plus deduplicated region→GENCODE gene/promoter→CGC/DGIdb join",
        "gencode": gtf_profile,
        "cosmic_genes": cosmic_genes_profile,
        "cgc": cgc_profile,
        "dgidb": dgidb_profile,
        "resistance": resistance_profile,
        "clp_hcc1395": clp_profile,
        "wakhan": wakhan_profile,
        "old_region_annotations": old_profile,
        "join": join_summary,
        "topology_pair": {
            "coarse_region_rows": coarse_region_rows,
            "all_match_rows": all_match_rows,
            "exact_coordinate_rows": exact_match_rows,
            **pair_checks,
        },
    }
    inventory = source_inventory(profiles, source_by_sample)
    write_tsv(OUT / "source_inventory.tsv", inventory)
    write_tsv(OUT / "hcc1395_pair_region_gene_join.tsv.gz", long_rows, gz=True)
    important_rows = [
        row
        for row in long_rows
        if row["cgc"] or row["dgidb"] or row["gene_symbol"] in {"BRCA2", "TBC1D16", "ERBB2", "MYC"}
    ]
    write_tsv(OUT / "hcc1395_pair_cancer_drug_region_hits.tsv", important_rows)
    write_tsv(OUT / "hcc1395_exact_complete_pair_gene_drug_flags.tsv", exact_pair_rows)
    write_tsv(OUT / "hcc1395_exact_complete_annotation_agreement.tsv", annotation_agreement)
    write_tsv(OUT / "hcc1395_exact_complete_clp_sensitivity.tsv", clp_sensitivity)
    write_tsv(OUT / "hcc1395_exact_complete_gene_summary.tsv", pair_gene_summary)
    write_tsv(OUT / "hcc1395_notable_gene_locus_sensitivity.tsv", notable_loci)
    qa_rows = []
    for sample, row in join_summary.items():
        for metric, value in row.items():
            qa_rows.append({"sample": sample, "metric": metric, "value": value})
    write_tsv(OUT / "hcc1395_pair_region_gene_drug_join_qa.tsv", qa_rows)
    (OUT / "gene_drug_source_profile.json").write_text(json.dumps(profiles, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"outputs": [
        str(OUT / "source_inventory.tsv"),
        str(OUT / "hcc1395_pair_region_gene_join.tsv.gz"),
        str(OUT / "hcc1395_pair_cancer_drug_region_hits.tsv"),
        str(OUT / "hcc1395_exact_complete_pair_gene_drug_flags.tsv"),
        str(OUT / "hcc1395_exact_complete_annotation_agreement.tsv"),
        str(OUT / "hcc1395_cosmic_clp_v104_chr1_22_alleles.tsv.gz"),
        str(OUT / "hcc1395_exact_complete_clp_sensitivity.tsv"),
        str(OUT / "hcc1395_exact_complete_gene_summary.tsv"),
        str(OUT / "hcc1395_notable_gene_locus_sensitivity.tsv"),
        str(OUT / "hcc1395_pair_region_gene_drug_join_qa.tsv"),
        str(OUT / "gene_drug_source_profile.json"),
    ], "join": join_summary, "pair_checks": pair_checks, "clp_hcc1395": clp_profile,
       "old_annotations": old_profile}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
