#!/usr/bin/env python3
"""Regression tests for the exact-PS layered workstation authority and panels."""

from __future__ import annotations

import importlib.util
import json
import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
BUILDER_PATH = (
    REPO_ROOT
    / "research"
    / "20260724_exact_ps_cpp_topology_signature_census"
    / "scripts"
    / "build_exact_ps_layered_workstation.py"
)
SPEC = importlib.util.spec_from_file_location("exact_ps_workstation_builder", BUILDER_PATH)
assert SPEC and SPEC.loader
BUILDER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILDER)


class ExactPsLayeredWorkstationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authority = BUILDER.load_authority()
        topology = BUILDER.sample_record_map(
            cls.authority["all7_summary"]["samples"]
        )
        census = BUILDER.sample_record_map(
            cls.authority["census_summary"]["samples"]
        )
        cls.topology = topology
        cls.census = census
        cls.hcc = BUILDER.load_sample(
            "HCC1395",
            cls.authority,
            topology["HCC1395"],
            census["HCC1395"],
            cls.authority["funnel"]["samples"]["HCC1395"],
        )
        cls.index_samples = [
            {
                "sample": sample,
                "summary": {
                    "cgcRegions": BUILDER.EXPECTED_GENE_DRUG_REGION_COUNTS[
                        sample
                    ]["cgc"],
                    "cgcDrugRegions": BUILDER.EXPECTED_GENE_DRUG_REGION_COUNTS[
                        sample
                    ]["cgc_drug"],
                },
            }
            for sample in BUILDER.EXPECTED_SAMPLES
        ]

    def test_all_seven_denominators_conserve(self) -> None:
        self.assertEqual(
            sum(int(self.topology[sample]["groups_total"]) for sample in BUILDER.EXPECTED_SAMPLES),
            98955,
        )
        self.assertEqual(
            sum(int(self.census[sample]["ranked_units"]) for sample in BUILDER.EXPECTED_SAMPLES),
            71955,
        )
        self.assertEqual(
            sum(
                int(self.topology[sample]["mutation_family_abstain_units"])
                for sample in BUILDER.EXPECTED_SAMPLES
            ),
            10717,
        )
        self.assertEqual(
            sum(
                int(self.topology[sample]["recurrence_required_units"])
                for sample in BUILDER.EXPECTED_SAMPLES
            ),
            45,
        )

    def test_status_partition_and_ranked_resolution(self) -> None:
        totals = self.authority["all7_summary"]["totals"]
        cohort = self.authority["census_summary"]["cohort"]
        self.assertEqual(
            {
                key: int(cohort["resolution"][key]["n"])
                for key in (
                    "UNIQUE_TREE",
                    "TIED_SAME_TOPOLOGY",
                    "TIED_CROSS_TOPOLOGY",
                )
            },
            {
                "UNIQUE_TREE": 39648,
                "TIED_SAME_TOPOLOGY": 23858,
                "TIED_CROSS_TOPOLOGY": 8449,
            },
        )
        self.assertEqual(int(totals["no_active_alt_units"]), 13014)
        self.assertEqual(int(totals["zero_denominator_units"]), 3224)
        self.assertEqual(int(totals["mutation_family_abstain_units"]), 10717)
        self.assertEqual(int(totals["recurrence_required_units"]), 45)

    def test_target_chr10_old_four_site_region_is_split(self) -> None:
        hcc = [
            row
            for row in self.hcc["rows"]
            if row["chrom"] == "chr10"
            and row["start"] <= 87928739
            and row["end"] >= 87818272
        ]
        self.assertEqual(len(hcc), 1)
        self.assertEqual(hcc[0]["positions"], [87818272, 87840023])
        self.assertEqual(hcc[0]["resolution"], "UNIQUE_TREE")
        self.assertEqual(hcc[0]["coarse"], "Direct-only")
        self.assertEqual(
            [(edge["parent_label"], edge["child_label"]) for edge in hcc[0]["edges"]],
            [("ROOT", "H_AR"), ("H_AR", "H_AA")],
        )
        self.assertNotIn(87888228, hcc[0]["positions"])
        self.assertNotIn(87928739, hcc[0]["positions"])
        self.assertEqual(
            hcc[0]["boundaryContext"]["endpointSupport"],
            {"total": 8, "RR": 6, "RA": 0, "AR": 0, "AA": 2},
        )
        self.assertEqual(hcc[0]["annotationClass"], "NO_CGC")
        self.assertFalse(hcc[0]["hasCgc"])
        self.assertFalse(hcc[0]["hasCgcDrug"])
        self.assertEqual(hcc[0]["geneDrug"], [])
        candidate = hcc[0]["candidate"]
        self.assertEqual(candidate["minimumSetCount"], 3)
        self.assertEqual(candidate["observedVertices"], [0])
        self.assertEqual(
            candidate["allEdges"],
            [[0, 1, "2"], [0, 2, "2"], [1, 3, "1"], [2, 3, "1"]],
        )
        self.assertEqual(
            candidate["bestEdges"],
            [[0, 1, "1"], [1, 3, "1"]],
        )

    def test_exact_locus_gene_drug_intersection_and_counts(self) -> None:
        self.assertEqual(
            (
                self.hcc["summary"]["cgcRegions"],
                self.hcc["summary"]["cgcDrugRegions"],
            ),
            (480, 174),
        )
        self.assertEqual(
            (
                sum(
                    row["cgc"]
                    for row in BUILDER.EXPECTED_GENE_DRUG_REGION_COUNTS.values()
                ),
                sum(
                    row["cgc_drug"]
                    for row in BUILDER.EXPECTED_GENE_DRUG_REGION_COUNTS.values()
                ),
            ),
            (3554, 1252),
        )
        hit = next(row for row in self.hcc["rows"] if row["hasCgcDrug"])
        self.assertEqual(hit["annotationClass"], "CGC_DRUG")
        for gene in hit["geneDrug"]:
            for position in gene["loci"]:
                self.assertIn(position, hit["positions"])
                self.assertLessEqual(gene["start"], position)
                self.assertLessEqual(position, gene["end"])
            for drug in gene["drugs"]:
                self.assertTrue(drug["claimName"])
                self.assertTrue(drug["sources"])
        self.assertEqual(
            self.authority["gene_drug_sha256"],
            "148f74e9f958cccffe14a2297e9453e28b2346a424adcb30e350f79e91ae4f50",
        )
        self.assertEqual(
            self.authority["gene_drug_sha256"],
            self.authority["gene_drug_receipt"]["output"]["sha256"],
        )

    def test_sample_html_embeds_one_authority_payload_for_all_panels(self) -> None:
        document = BUILDER.sample_page(
            self.authority, self.hcc
        )
        self.assertIn(
            f'<meta name="intersubmod-authority" content="{BUILDER.AUTHORITY_NAME}">',
            document,
        )
        index = BUILDER.index_page(self.authority, self.index_samples)
        self.assertIn("98,955", index)
        self.assertIn("七組資料的狀態、拓撲與複雜度並排比較", index)
        self.assertIn("CGC ∩ drug", index)
        self.assertIn("cohortData", index)
        cohort_match = re.search(
            r'<script type="application/json" id="cohortData">(.*?)</script>',
            index,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(cohort_match)
        cohort_payload = json.loads(cohort_match.group(1))
        self.assertEqual(
            cohort_payload["sampleOrder"],
            list(BUILDER.DISPLAY_ORDER),
        )
        self.assertEqual(
            cohort_payload["displayLabels"],
            BUILDER.DISPLAY_LABELS,
        )
        sample_links = re.findall(
            r'<a class="sample-link"[^>]*>([^<]+)</a>',
            index,
        )
        self.assertEqual(
            sample_links,
            [BUILDER.DISPLAY_LABELS[sample] for sample in BUILDER.DISPLAY_ORDER],
        )
        match = re.search(
            r'<script type="application/json" id="pageData">(.*?)</script>',
            document,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(match)
        payload = json.loads(match.group(1))
        self.assertEqual(len(payload["rows"]), 11590)
        self.assertEqual(payload["summary"]["ranked"], 9130)
        self.assertEqual(
            sum(row["resolution"] == "UNIQUE_TREE" for row in payload["rows"]),
            7047,
        )
        self.assertIn("不再是四點同一區", document)
        self.assertIn("完整 minimum candidate space", document)
        self.assertIn("癌症基因 ∩ approved 抗腫瘤藥", document)
        self.assertIn('data-testid="annotation-filters"', document)
        self.assertIn('data-testid="gene-drug-panel"', document)
        self.assertIn("S3=", document)
        self.assertIn("機器證據與原始 JSON（預設收合）", document)

    def test_generated_pages_are_bound_when_present(self) -> None:
        index = BUILDER.OUTDIR / "index.html"
        if not index.is_file():
            self.skipTest("generated pages have not been built yet")
        index_prefix = index.read_text(encoding="utf-8")[:12000]
        if BUILDER.AUTHORITY_NAME not in index_prefix:
            self.skipTest("generated pages still use a different authority; rebuild first")
        for sample in BUILDER.EXPECTED_SAMPLES:
            page = BUILDER.OUTDIR / f"{sample}.html"
            self.assertTrue(page.is_file(), page)
            prefix = page.read_text(encoding="utf-8")[:12000]
            self.assertIn(BUILDER.AUTHORITY_NAME, prefix)
            self.assertIn(BUILDER.UI_CONTRACT, prefix)
            self.assertIn("intersubmod-gene-drug-sha256", prefix)


if __name__ == "__main__":
    unittest.main()
