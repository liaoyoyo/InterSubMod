#!/usr/bin/env python3
"""Regression tests for the exact-PS cancer-gene/drug annotation sidecar."""

from __future__ import annotations

import csv
import copy
import gzip
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
PRODUCER_PATH = (
    REPO_ROOT
    / "research"
    / "20260724_exact_ps_cpp_topology_signature_census"
    / "scripts"
    / "build_exact_ps_gene_drug_annotation.py"
)
SPEC = importlib.util.spec_from_file_location(
    "exact_ps_gene_drug_annotation",
    PRODUCER_PATH,
)
assert SPEC and SPEC.loader
PRODUCER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRODUCER)


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ExactPsGeneDrugAnnotationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.gencode = (
            self.root / "gencode.v46.basic.annotation.gtf.gz"
        )
        self.cgc = (
            self.root / "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz"
        )
        self.cosmic_genes = (
            self.root / "Cosmic_Genes_v104_GRCh38.tsv.gz"
        )
        self.dgidb = self.root / "dgidb_interactions.tsv"
        self._write_gencode()
        self._write_cosmic_genes()
        self._write_cgc()
        self._write_dgidb()

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _write_gencode(self) -> None:
        lines = [
            (
                "##description: evidence-based annotation of the human genome "
                "(GRCh38), version 46 (Ensembl 112)\n"
            ),
            "##provider: GENCODE\n",
            "##format: gtf\n",
            (
                'chr1\tHAVANA\tgene\t90\t210\t.\t+\t.\tgene_id '
                '"ENSG_TEST_TP53.1"; gene_type "protein_coding"; gene_name '
                '"TP53"; hgnc_id "HGNC:11998";\n'
            ),
            (
                'chr2\tHAVANA\tgene\t300\t450\t.\t-\t.\tgene_id '
                '"ENSG_TEST_FALLBACK.1"; gene_type "protein_coding"; '
                'gene_name "FALLBACK";\n'
            ),
            (
                'chr3\tHAVANA\tgene\t500\t600\t.\t+\t.\tgene_id '
                '"ENSG_TEST_NOT_CGC.1"; gene_type "protein_coding"; '
                'gene_name "NOTCGC"; hgnc_id "HGNC:333";\n'
            ),
            (
                'chr4\tHAVANA\tgene\t700\t800\t.\t+\t.\tgene_id '
                '"ENSG_TEST_NODRUG.1"; gene_type "protein_coding"; '
                'gene_name "NODRUG"; hgnc_id "HGNC:444";\n'
            ),
            (
                'chrX\tHAVANA\tgene\t900\t1000\t.\t+\t.\tgene_id '
                '"ENSG_TEST_X.1"; gene_type "protein_coding"; gene_name '
                '"XGENE"; hgnc_id "HGNC:999";\n'
            ),
        ]
        with gzip.open(self.gencode, "wt", encoding="utf-8", newline="") as handle:
            handle.writelines(lines)

    def _write_cosmic_genes(self) -> None:
        rows = [
            ["COSG1", "", "TP53", "", "", "11998", "y", "y"],
            ["COSG2", "", "FALLBACK", "", "", "", "y", "y"],
            ["COSG4", "", "NODRUG", "", "", "444", "y", "y"],
            ["COSG9", "", "XGENE", "", "", "999", "y", "y"],
        ]
        with gzip.open(
            self.cosmic_genes,
            "wt",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(PRODUCER.COSMIC_GENES_HEADER)
            writer.writerows(rows)

    def _cgc_row(
        self,
        symbol: str,
        cosmic_gene_id: str,
        chrom: str,
        start: str,
        end: str,
    ) -> list[str]:
        values = {field: "" for field in PRODUCER.CGC_HEADER}
        values.update(
            {
                "GENE_SYMBOL": symbol,
                "NAME": f"{symbol} name",
                "COSMIC_GENE_ID": cosmic_gene_id,
                "CHROMOSOME": chrom,
                "GENOME_START": start,
                "GENOME_STOP": end,
                "SOMATIC": "y",
                "GERMLINE": "n",
                "TUMOUR_TYPES_SOMATIC": "test tumour",
                "ROLE_IN_CANCER": "oncogene",
                "TIER": "1",
            }
        )
        return [values[field] for field in PRODUCER.CGC_HEADER]

    def _write_cgc(self) -> None:
        rows = [
            self._cgc_row("TP53", "COSG1", "1", "100", "200"),
            self._cgc_row("FALLBACK", "COSG2", "2", "310", "440"),
            self._cgc_row("NODRUG", "COSG4", "4", "700", "800"),
            self._cgc_row("XGENE", "COSG9", "X", "900", "1000"),
        ]
        with gzip.open(
            self.cgc,
            "wt",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(PRODUCER.CGC_HEADER)
            writer.writerows(rows)

    def _dgidb_row(
        self,
        *,
        gene_claim_name: str,
        gene_concept_id: str,
        gene_name: str,
        source: str,
        version: str,
        drug_claim_name: str,
        drug_concept_id: str,
        drug_name: str,
        approved: str,
        anti_neoplastic: str,
    ) -> list[str]:
        values = {
            "gene_claim_name": gene_claim_name,
            "gene_concept_id": gene_concept_id,
            "gene_name": gene_name,
            "interaction_source_db_name": source,
            "interaction_source_db_version": version,
            "interaction_type": "NULL",
            "interaction_score": "1.0",
            "drug_claim_name": drug_claim_name,
            "drug_concept_id": drug_concept_id,
            "drug_name": drug_name,
            "approved": approved,
            "immunotherapy": "FALSE",
            "anti_neoplastic": anti_neoplastic,
        }
        return [values[field] for field in PRODUCER.DGIDB_HEADER]

    def _write_dgidb(self) -> None:
        selected = self._dgidb_row(
            gene_claim_name="TP53",
            gene_concept_id="hgnc:11998",
            gene_name="TP53",
            source="SourceA",
            version="1",
            drug_claim_name="CLAIM-A",
            drug_concept_id="drug:1",
            drug_name="CANONICAL-X",
            approved="TRUE",
            anti_neoplastic="TRUE",
        )
        rows = [
            selected,
            selected,
            self._dgidb_row(
                gene_claim_name="TP53",
                gene_concept_id="hgnc:11998",
                gene_name="TP53",
                source="SourceB",
                version="2",
                drug_claim_name="CLAIM-A",
                drug_concept_id="drug:1",
                drug_name="CANONICAL-X",
                approved="TRUE",
                anti_neoplastic="TRUE",
            ),
            self._dgidb_row(
                gene_claim_name="TP53",
                gene_concept_id="hgnc:11998",
                gene_name="TP53",
                source="SourceA",
                version="1",
                drug_claim_name="CLAIM-B",
                drug_concept_id="drug:1",
                drug_name="CANONICAL-X",
                approved="TRUE",
                anti_neoplastic="TRUE",
            ),
            self._dgidb_row(
                gene_claim_name="TP53",
                gene_concept_id="hgnc:11998",
                gene_name="TP53",
                source="RejectedA",
                version="1",
                drug_claim_name="APPROVED-ONLY",
                drug_concept_id="drug:2",
                drug_name="APPROVED-ONLY",
                approved="TRUE",
                anti_neoplastic="FALSE",
            ),
            self._dgidb_row(
                gene_claim_name="TP53",
                gene_concept_id="hgnc:11998",
                gene_name="TP53",
                source="RejectedB",
                version="1",
                drug_claim_name="ANTINEOPLASTIC-ONLY",
                drug_concept_id="drug:3",
                drug_name="ANTINEOPLASTIC-ONLY",
                approved="FALSE",
                anti_neoplastic="TRUE",
            ),
            self._dgidb_row(
                gene_claim_name="FALLBACK",
                gene_concept_id="NULL",
                gene_name="FALLBACK",
                source="SourceC",
                version="3",
                drug_claim_name="FALLBACK-DRUG",
                drug_concept_id="drug:4",
                drug_name="FALLBACK-DRUG",
                approved="TRUE",
                anti_neoplastic="TRUE",
            ),
            self._dgidb_row(
                gene_claim_name="NOTCGC",
                gene_concept_id="hgnc:333",
                gene_name="NOTCGC",
                source="SourceD",
                version="4",
                drug_claim_name="NOT-CGC-DRUG",
                drug_concept_id="drug:5",
                drug_name="NOT-CGC-DRUG",
                approved="TRUE",
                anti_neoplastic="TRUE",
            ),
            self._dgidb_row(
                gene_claim_name="RAS",
                gene_concept_id="NULL",
                gene_name="NULL",
                source="RejectedFamily",
                version="1",
                drug_claim_name="RAS-FAMILY-DRUG",
                drug_concept_id="drug:6",
                drug_name="RAS-FAMILY-DRUG",
                approved="TRUE",
                anti_neoplastic="TRUE",
            ),
            self._dgidb_row(
                gene_claim_name="PI3K/AKT/MTOR PATHWAY",
                gene_concept_id="NULL",
                gene_name="NULL",
                source="RejectedPathway",
                version="1",
                drug_claim_name="PATHWAY-DRUG",
                drug_concept_id="drug:7",
                drug_name="PATHWAY-DRUG",
                approved="TRUE",
                anti_neoplastic="TRUE",
            ),
            self._dgidb_row(
                gene_claim_name="FREE TEXT TARGET",
                gene_concept_id="NULL",
                gene_name="NULL",
                source="RejectedFreeText",
                version="1",
                drug_claim_name="FREE-TEXT-DRUG",
                drug_concept_id="drug:8",
                drug_name="FREE-TEXT-DRUG",
                approved="TRUE",
                anti_neoplastic="TRUE",
            ),
        ]
        with self.dgidb.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(PRODUCER.DGIDB_HEADER)
            writer.writerows(rows)

    def _build(self) -> tuple[dict, dict]:
        return PRODUCER.build(
            gencode_path=self.gencode,
            cgc_path=self.cgc,
            cosmic_genes_path=self.cosmic_genes,
            dgidb_path=self.dgidb,
        )

    def test_schema_scope_sources_and_claim_ceiling(self) -> None:
        annotation, receipt = self._build()
        self.assertEqual(
            annotation["schema_name"],
            "intersubmod.exact_ps_gene_drug_annotation",
        )
        self.assertEqual(annotation["schema_version"], "1.1.0")
        self.assertEqual(
            receipt["schema_name"],
            "intersubmod.exact_ps_gene_drug_annotation.receipt",
        )
        self.assertEqual(receipt["schema_version"], "1.1.0")
        self.assertTrue(receipt["all_pass"])
        self.assertTrue(all(receipt["checks"].values()))
        self.assertEqual(receipt["output"]["path"], "annotation.v1.json")
        self.assertEqual(annotation["scope"]["chromosomes"][0], "chr1")
        self.assertEqual(annotation["scope"]["chromosomes"][-1], "chr22")
        self.assertFalse(annotation["scope"]["partial"])
        self.assertEqual(
            annotation["scope"]["coordinate_system"],
            "GENCODE_v46_gene_body_1_based_inclusive",
        )
        self.assertEqual(set(annotation["sources"]), {
            "gencode",
            "cgc",
            "cosmic_genes",
            "dgidb",
        })
        for source in receipt["raw_sources"].values():
            self.assertEqual(len(source["sha256"]), 64)
            self.assertGreater(source["size_bytes"], 0)
        self.assertEqual(
            annotation["claim_ceiling"]["level"],
            "gene_level_annotation_only",
        )
        self.assertIn(
            "clinical_actionability",
            annotation["claim_ceiling"]["prohibited_claims"],
        )

    def test_gencode_body_is_position_authority_and_same_gene_and(self) -> None:
        annotation, _ = self._build()
        genes = {row["symbol"]: row for row in annotation["genes"]}
        self.assertEqual(set(genes), {"TP53", "FALLBACK", "NODRUG"})
        tp53 = genes["TP53"]
        self.assertEqual((tp53["chrom"], tp53["start"], tp53["end"]), (
            "chr1",
            90,
            210,
        ))
        self.assertEqual(
            tp53["cgc_coordinate_crosscheck"]["status"],
            "overlap_nonidentical",
        )

        def matching_symbols(chrom: str, position: int) -> set[str]:
            return {
                row["symbol"]
                for row in annotation["genes"]
                if row["chrom"] == chrom
                and row["start"] <= position <= row["end"]
            }

        self.assertEqual(matching_symbols("chr1", 95), {"TP53"})
        self.assertEqual(matching_symbols("chr1", 210), {"TP53"})
        self.assertEqual(matching_symbols("chr1", 211), set())
        self.assertNotIn("NOTCGC", genes)
        self.assertEqual(
            genes["NODRUG"]["approved_antineoplastic_drugs"],
            [],
        )

    def test_strict_approved_and_antineoplastic_filter_and_duplicate_collapse(
        self,
    ) -> None:
        annotation, _ = self._build()
        tp53 = next(row for row in annotation["genes"] if row["symbol"] == "TP53")
        drugs = tp53["approved_antineoplastic_drugs"]
        self.assertEqual([row["claim_name"] for row in drugs], [
            "CLAIM-A",
            "CLAIM-B",
        ])
        self.assertEqual(
            [row["canonical_names"] for row in drugs],
            [["CANONICAL-X"], ["CANONICAL-X"]],
        )
        self.assertEqual(
            drugs[0]["sources"],
            [
                {"name": "SourceA", "version": "1"},
                {"name": "SourceB", "version": "2"},
            ],
        )
        self.assertEqual(annotation["counts"]["cgc_genes"], 3)
        self.assertEqual(
            annotation["counts"]["cgc_with_approved_antineoplastic"],
            2,
        )
        self.assertEqual(
            annotation["counts"]["dgidb_approved_and_antineoplastic_rows"],
            9,
        )
        self.assertEqual(
            annotation["counts"]["gene_drug_claim_pairs_after_collapse"],
            3,
        )

    def test_hgnc_primary_and_symbol_fallback_only_when_hgnc_missing(self) -> None:
        annotation, _ = self._build()
        genes = {row["symbol"]: row for row in annotation["genes"]}
        self.assertEqual(genes["TP53"]["gene_match_basis"], "cosg_to_hgnc")
        self.assertEqual(genes["TP53"]["drug_match_basis"], "hgnc_id")
        self.assertEqual(
            genes["FALLBACK"]["gene_match_basis"],
            "symbol_fallback_missing_hgnc",
        )
        self.assertEqual(
            genes["FALLBACK"]["drug_match_basis"],
            "symbol_fallback_missing_hgnc",
        )
        self.assertEqual(
            genes["FALLBACK"]["approved_antineoplastic_drugs"][0][
                "claim_name"
            ],
            "FALLBACK-DRUG",
        )

    def test_family_pathway_and_free_text_symbol_fallback_are_rejected(
        self,
    ) -> None:
        annotation, _ = self._build()
        counts = annotation["counts"]
        self.assertEqual(
            counts["dgidb_selected_rows_indexed_by_symbol_only"],
            4,
        )
        self.assertEqual(
            counts["dgidb_symbol_only_rows_accepted_after_gencode_gate"],
            1,
        )
        self.assertEqual(
            counts["dgidb_selected_rows_without_indexable_gene"],
            3,
        )
        self.assertEqual(
            counts["dgidb_symbol_only_keys_overlapping_cgc"],
            1,
        )
        all_claims = {
            drug["claim_name"]
            for gene in annotation["genes"]
            for drug in gene["approved_antineoplastic_drugs"]
        }
        self.assertTrue(
            {
                "RAS-FAMILY-DRUG",
                "PATHWAY-DRUG",
                "FREE-TEXT-DRUG",
            }.isdisjoint(all_claims)
        )

    def test_receipt_checks_are_derived_from_payload(self) -> None:
        annotation, receipt = self._build()
        derived = PRODUCER.derive_receipt_checks(
            annotation,
            receipt["raw_sources"],
            receipt["output"],
        )
        self.assertEqual(derived, receipt["checks"])
        self.assertTrue(all(derived.values()))

        corrupted = copy.deepcopy(annotation)
        corrupted["counts"][
            "dgidb_selected_rows_indexed_by_symbol_only"
        ] += 1
        corrupted_checks = PRODUCER.derive_receipt_checks(
            corrupted,
            receipt["raw_sources"],
            receipt["output"],
        )
        self.assertFalse(
            corrupted_checks["symbol_fallback_only_when_hgnc_missing"]
        )
        self.assertFalse(
            corrupted_checks["deterministic_serialization_contract"]
        )

    def test_chr_x_is_excluded(self) -> None:
        annotation, _ = self._build()
        self.assertNotIn(
            "chrX",
            {row["chrom"] for row in annotation["genes"]},
        )
        self.assertEqual(
            annotation["counts"]["cgc_excluded_non_autosomal"],
            1,
        )
        self.assertEqual(
            annotation["counts"][
                "gencode_excluded_non_autosomal_gene_bodies"
            ],
            1,
        )

    def test_outputs_and_receipt_are_byte_deterministic(self) -> None:
        annotation_a, receipt_a = self._build()
        annotation_b, receipt_b = self._build()
        out_a = self.root / "out-a"
        out_b = self.root / "out-b"
        PRODUCER.publish(out_a, annotation_a, receipt_a)
        PRODUCER.publish(out_b, annotation_b, receipt_b)
        for file_name in ("annotation.v1.json", "receipt.json"):
            self.assertEqual(
                (out_a / file_name).read_bytes(),
                (out_b / file_name).read_bytes(),
            )
            self.assertEqual(
                sha256_path(out_a / file_name),
                sha256_path(out_b / file_name),
            )
        parsed_receipt = json.loads((out_a / "receipt.json").read_text())
        annotation_path = out_a / "annotation.v1.json"
        self.assertEqual(
            parsed_receipt["output"]["sha256"],
            sha256_path(annotation_path),
        )
        self.assertEqual(
            parsed_receipt["output"]["size_bytes"],
            annotation_path.stat().st_size,
        )

    def test_bad_source_schema_fails_closed(self) -> None:
        with gzip.open(
            self.cosmic_genes,
            "wt",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(PRODUCER.COSMIC_GENES_HEADER[:-1])
            writer.writerow(["COSG1"] * (len(PRODUCER.COSMIC_GENES_HEADER) - 1))
        with self.assertRaisesRegex(
            PRODUCER.AnnotationBuildError,
            "schema mismatch",
        ):
            self._build()

    def test_bad_boolean_fails_closed(self) -> None:
        with self.dgidb.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        approved_index = list(PRODUCER.DGIDB_HEADER).index("approved")
        rows[1][approved_index] = "YES"
        with self.dgidb.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(rows)
        with self.assertRaisesRegex(
            PRODUCER.AnnotationBuildError,
            "invalid approved boolean",
        ):
            self._build()


class ExactPsGeneDrugAnnotationFullSourceContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        root = Path("/big7_disk/liaoyoyo2001/gene_annotation")
        cls.annotation, cls.receipt = PRODUCER.build(
            gencode_path=root / "gencode.v46.basic.annotation.gtf.gz",
            cgc_path=root / "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz",
            cosmic_genes_path=(
                root
                / "cosmic_v104"
                / "Cosmic_Genes_v104_GRCh38.tsv.gz"
            ),
            dgidb_path=root / "dgidb_interactions.tsv",
        )

    def test_full_source_symbol_only_contract_and_stable_gene_counts(
        self,
    ) -> None:
        counts = self.annotation["counts"]
        self.assertEqual(
            counts["dgidb_selected_rows_indexed_by_symbol_only"],
            224,
        )
        self.assertEqual(
            counts["dgidb_symbol_only_keys_overlapping_cgc"],
            0,
        )
        self.assertEqual(counts["cgc_genes"], 721)
        self.assertEqual(
            counts["cgc_with_approved_antineoplastic"],
            279,
        )
        self.assertTrue(self.receipt["all_pass"])
        self.assertTrue(all(self.receipt["checks"].values()))


if __name__ == "__main__":
    unittest.main()
