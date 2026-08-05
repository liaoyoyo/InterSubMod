#!/usr/bin/env python3
"""Fail-closed validation for the gene/drug inventory subtask outputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUT = ROOT / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    profile = json.loads((OUT / "gene_drug_source_profile.json").read_text())
    assert profile["topology_pair"]["coarse_region_rows"] == 47377
    assert profile["topology_pair"]["exact_coordinate_rows"] == 6252
    assert profile["topology_pair"]["exact_complete_pair_rows"] == 5720
    assert profile["topology_pair"]["baseline_agree_n"] == 3969
    assert profile["topology_pair"]["all_strata_conserve_5720"] is True
    assert profile["clp_hcc1395"]["hcc1395_chr1_22_unique_genomic_alleles"] == 1141
    assert profile["clp_hcc1395"]["confirmed_somatic_unique_alleles"] == 33
    assert profile["clp_hcc1395"]["unknown_origin_unique_alleles"] == 1108
    assert profile["clp_hcc1395"]["all_unique_alleles_overlapping_exact_complete_pair"] == 381

    inventory = read_tsv(OUT / "source_inventory.tsv")
    assert len(inventory) == len({row["source_id"] for row in inventory})
    for row in inventory:
        path = Path(row["absolute_path"])
        assert path.is_file(), path
        assert int(row["file_bytes"]) == path.stat().st_size, path
        assert row["sha256"] == sha256(path), path

    flags = read_tsv(OUT / "hcc1395_exact_complete_pair_gene_drug_flags.tsv")
    assert len(flags) == 5720
    assert len({row["match_id"] for row in flags}) == 5720
    assert sum(row["category_agree"] == "True" for row in flags) == 3969

    agreement = read_tsv(OUT / "hcc1395_exact_complete_annotation_agreement.tsv")
    by_feature = defaultdict(list)
    for row in agreement:
        by_feature[row["feature"]].append(row)
    assert len(by_feature["ALL"]) == 1
    assert int(by_feature["ALL"][0]["n_exact_complete_pairs"]) == 5720
    for feature, rows in by_feature.items():
        if feature == "ALL":
            continue
        assert len(rows) == 2, feature
        assert sum(int(row["n_exact_complete_pairs"]) for row in rows) == 5720, feature

    with gzip.open(OUT / "hcc1395_cosmic_clp_v104_chr1_22_alleles.tsv.gz", "rt") as handle:
        clp = list(csv.DictReader(handle, delimiter="\t"))
    assert len(clp) == 1141
    assert sum(row["confirmed_somatic"] == "True" for row in clp) == 33
    assert sum(row["unknown_origin"] == "True" for row in clp) == 1108
    assert sum(bool(row["overlapping_exact_complete_regions"]) for row in clp) == 381
    assert {row["cosmic_sample_ids"] for row in clp} == {"COSS749712"}

    clp_sensitivity = read_tsv(OUT / "hcc1395_exact_complete_clp_sensitivity.tsv")
    assert len(clp_sensitivity) == 6
    for idx in range(0, len(clp_sensitivity), 2):
        assert sum(int(row["n_exact_complete_pairs"]) for row in clp_sensitivity[idx : idx + 2]) == 5720

    notable = read_tsv(OUT / "hcc1395_notable_gene_locus_sensitivity.tsv")
    assert {row["gene_symbol"] for row in notable} == {"BRCA2", "TBC1D16", "ERBB2", "MYC"}
    brca2 = next(row for row in notable if row["gene_symbol"] == "BRCA2")
    assert brca2["best_pair_scenario"] == "reciprocal_overlap_0.80"
    assert brca2["best_pair_category_agree"] == "True"

    report = (OUT / "gene_drug_inventory.md").read_text()
    for token in (
        "Interaction",
        "Approved ∩ anti-neoplastic",
        "6,055,938",
        "COSS749712",
        "69.3881%",
        "Final verdict：資料工程 PASS；生物／臨床驗證 NO-GO",
    ):
        assert token in report, token

    print(
        json.dumps(
            {
                "status": "PASS",
                "inventory_sources": len(inventory),
                "exact_complete_pairs": len(flags),
                "annotation_features": len(by_feature) - 1,
                "clp_unique_alleles": len(clp),
                "notable_genes": len(notable),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
