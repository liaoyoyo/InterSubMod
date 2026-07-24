#!/usr/bin/env python3
"""Regression tests for the exact-PS layered workstation authority and panels."""

from __future__ import annotations

import importlib.util
import json
import re
import unittest
from collections import Counter
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
        cls.result = BUILDER.build(verify_only=True)
        cls.by_sample = {
            item["sample"]: item for item in cls.result["samples"]
        }

    def test_all_seven_denominators_conserve(self) -> None:
        self.assertEqual(set(self.by_sample), set(BUILDER.EXPECTED_SAMPLES))
        self.assertEqual(
            sum(item["summary"]["groups"] for item in self.by_sample.values()),
            98955,
        )
        self.assertEqual(
            sum(item["summary"]["ranked"] for item in self.by_sample.values()),
            71955,
        )
        self.assertEqual(
            sum(item["summary"]["abstain"] for item in self.by_sample.values()),
            10717,
        )
        self.assertEqual(
            sum(item["summary"]["recurrence"] for item in self.by_sample.values()),
            45,
        )

    def test_status_partition_and_ranked_resolution(self) -> None:
        status = Counter()
        for item in self.by_sample.values():
            status.update(row["resolution"] for row in item["rows"])
        self.assertEqual(
            status,
            Counter(
                {
                    "UNIQUE_TREE": 39648,
                    "TIED_SAME_TOPOLOGY": 23858,
                    "TIED_CROSS_TOPOLOGY": 8449,
                    "NO_ACTIVE_ALT": 13014,
                    "ZERO_DENOMINATOR": 3224,
                    "RESOURCE_ABSTAIN": 10717,
                    "RECURRENCE_REQUIRED": 45,
                }
            ),
        )

    def test_target_chr10_old_four_site_region_is_split(self) -> None:
        def overlaps(item: dict) -> list[dict]:
            return [
                row
                for row in item["rows"]
                if row["chrom"] == "chr10"
                and row["start"] <= 87928739
                and row["end"] >= 87818272
            ]

        hcc = overlaps(self.by_sample["HCC1395"])
        dorado = overlaps(self.by_sample["HCC1395_DORADO"])
        self.assertEqual(len(hcc), 1)
        self.assertEqual(hcc[0]["positions"], [87818272, 87840023])
        self.assertEqual(hcc[0]["resolution"], "UNIQUE_TREE")
        self.assertEqual(hcc[0]["coarse"], "Direct-only")
        self.assertEqual(
            [(edge["parent_label"], edge["child_label"]) for edge in hcc[0]["edges"]],
            [("ROOT", "H_AR"), ("H_AR", "H_AA")],
        )
        self.assertEqual(dorado, [])
        self.assertNotIn(87888228, hcc[0]["positions"])
        self.assertNotIn(87928739, hcc[0]["positions"])

    def test_sample_html_embeds_one_authority_payload_for_all_panels(self) -> None:
        document = BUILDER.sample_page(
            self.result["authority"], self.by_sample["HCC1395"]
        )
        self.assertIn(
            f'<meta name="intersubmod-authority" content="{BUILDER.AUTHORITY_NAME}">',
            document,
        )
        self.assertIn("98,955", BUILDER.index_page(
            self.result["authority"], list(self.by_sample.values())
        ))
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
            Counter(row["resolution"] for row in payload["rows"])["UNIQUE_TREE"],
            7047,
        )
        self.assertIn("不再是四點同一區", document)
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


if __name__ == "__main__":
    unittest.main()
