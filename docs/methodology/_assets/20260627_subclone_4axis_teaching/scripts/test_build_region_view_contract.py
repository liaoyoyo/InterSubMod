#!/usr/bin/env python3
"""Regression contract for zero-population multilocus regions."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent


class BuildRegionViewContractTest(unittest.TestCase):
    def test_zero_population_group_remains_in_region_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            layered = root / "layered.json"
            part = root / "mlhp_part_4.json"
            vcf = root / "tree.vcf"
            output = root / "region.json"

            layered.write_text(
                json.dumps({
                    "detail": [],
                    "L0_hp_family": {"regions_total": 1},
                    "L1_ssnv_algorithm": {
                        "n_units_total_including_unphased": 0,
                        "n_primary_lineage_units": 0,
                        "n_reference_only_controls": 0,
                        "n_verification_fail": 0,
                    },
                    "L2_cn": {},
                    "L3_methyl": {},
                    "analysis_contract": {},
                    "read_tag_census": {},
                    "input_funnel": {
                        "n_sSNV_scope_input": 8,
                        "n_positional_singleton": 0,
                        "n_multilocus_pre_cap_groups": 1,
                        "n_multilocus_pre_cap_sSNV": 8,
                        "n_groups_capped_by_MAX_SNV": 0,
                        "n_sSNV_cap_excluded": 0,
                        "n_groups_read_unsupported": 0,
                        "n_sSNV_read_unsupported": 0,
                        "n_groups_retained": 1,
                        "n_sSNV_retained": 8,
                        "n_sSNV_accounted": 8,
                        "check_scope_conservation": True,
                        "grouping_rule": "adjacent_gap<=50000; total_span_not_bounded",
                    },
                }),
                encoding="utf-8",
            )
            part.write_text(
                json.dumps({
                    "groups": [{
                        "chrom": "chr9",
                        "start": 61184408,
                        "end": 61184992,
                        "n_sSNV": 8,
                        "cn": "unavailable",
                        "positions": [
                            61184408, 61184547, 61184556, 61184682,
                            61184811, 61184918, 61184970, 61184992,
                        ],
                        "n_full_cov_reads": 5,
                        "populations_by_hp": {},
                        "subread_groups_by_hp": {},
                        "col_coverage_by_hp": {},
                    }],
                }),
                encoding="utf-8",
            )
            records = "".join(
                f"chr9\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n"
                for position in range(61184408, 61184416)
            )
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                + records,
                encoding="utf-8",
            )
            environment = os.environ.copy()
            environment.update({
                "SM_LAYERED": str(layered),
                "SM_OUT": str(output),
                "SM_SAMPLE": "HCC1937_ZERO_POPULATION_FIXTURE",
                "SM_ML_GLOB": str(root / "mlhp_part_*.json"),
                "SM_SOMATIC_VCF": str(vcf),
            })

            process = subprocess.run(
                [sys.executable, str(HERE / "build_region_view.py")],
                env=environment,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(process.returncode, 0, process.stderr)
            document = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(document["census"]["n_regions"], 1)
            self.assertEqual(document["census"]["n_regions_without_primary_lineage"], 1)
            self.assertTrue(document["census"]["funnel"]["check_six_branch_conservation"])
            self.assertEqual(len(document["regions"]), 1)
            region = document["regions"][0]
            self.assertEqual(
                (region["chrom"], region["start"], region["end"]),
                ("chr9", 61184408, 61184992),
            )
            self.assertEqual(region["hp_multiplicity"], 0)
            self.assertEqual(region["region_determinacy"], "no_primary_lineage")
            self.assertEqual(region["lineages"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
