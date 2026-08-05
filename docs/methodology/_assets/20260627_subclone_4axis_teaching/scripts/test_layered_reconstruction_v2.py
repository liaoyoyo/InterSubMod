#!/usr/bin/env python3
"""Regression tests for layered reconstruction v2 semantics."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


HERE = Path(__file__).resolve().parent
os.environ.setdefault("SM_CN_INT_GAIN", "")
os.environ.setdefault("SM_CN_INT_LOSS", "")
os.environ.setdefault("SM_VERIFY_EVERY", "1")
os.environ.setdefault("SM_ANALYSIS_TREE_CAP", "0")
os.environ.setdefault("SM_DISPLAY_TREE_CAP", "1")
sys.path.insert(0, str(HERE))

import layered_tree_reconstruction as L  # noqa: E402
import derive_minread4_run as D  # noqa: E402
import read_af_tree_ordering_multisample as RA  # noqa: E402
import sm_multilocus_combinations as ML  # noqa: E402
import tree_enumeration_solver as S  # noqa: E402


def synthetic_group():
    positions = [100, 200, 300]
    return {
        "chrom": "chr1",
        "start": 100,
        "end": 300,
        "span": 200,
        "n_sSNV": 3,
        "positions": positions,
        "n_full_cov_reads": 20,
        "cn": "unavailable",
        "populations": {"RRR": 5, "ARR": 5, "RAR": 5, "RRA": 5},
        "subread_groups": {"RRR": 5, "ARR": 5, "RAR": 5, "RRA": 5},
        "populations_by_hp": {
            "1": {"RRR": 5},
            "2": {"ARR": 5},
            "3": {"RAR": 5},
            "4": {"AAR": 5},
            "none": {"RRA": 5},
        },
        "subread_groups_by_hp": {},
        "col_coverage_by_hp": {
            "1": {"100": [5, 0], "200": [5, 0], "300": [5, 0]},
            "2": {"100": [0, 5], "200": [5, 0], "300": [5, 0]},
            "3": {"100": [5, 0], "200": [0, 5], "300": [5, 0]},
            "4": {"100": [0, 5], "200": [0, 5], "300": [5, 0]},
            "none": {"100": [5, 0], "200": [5, 0], "300": [0, 5]},
        },
        "reads_by_hp": {"1": 5, "2": 5, "3": 5, "4": 5, "none": 5},
    }


def synthetic_input_doc():
    return {
        "schema_version": "2.0",
        "params": {
            "MINREAD": 3,
            "MAX_SNV": 8,
            "TIER_R": 50000,
            "MAPQ_MIN": 20,
            "BASEQ_MIN": 0,
            "cn_data_available": False,
            "cn_source": "unavailable",
        },
        "input_funnel": {
            "n_sSNV_scope_input": 6,
            "n_positional_singleton": 2,
            "n_multilocus_pre_cap_groups": 2,
            "n_multilocus_pre_cap_sSNV": 4,
            "n_groups_capped_by_MAX_SNV": 1,
            "n_sSNV_cap_excluded": 1,
            "n_groups_read_unsupported": 1,
            "n_sSNV_read_unsupported": 0,
            "n_groups_retained": 1,
            "n_sSNV_retained": 3,
            "n_sSNV_accounted": 6,
            "check_scope_conservation": True,
            "grouping_rule": "adjacent_gap<=50000; total_span_not_bounded",
        },
        "groups": [synthetic_group()],
    }


class LayeredV2Test(unittest.TestCase):
    def setUp(self):
        L._unit_idx[0] = 0

    def test_root_only_tree_has_zero_hidden(self):
        result = S.enumerate_min_trees({"RRR": 5}, [], 3, tree_cap=0)
        self.assertEqual(result["n_hidden"], 0)
        self.assertTrue(result["trees_complete"])
        self.assertEqual(result["trees"][0]["n_hidden"], 0)

    def test_missing_cn_is_unavailable(self):
        self.assertEqual(ML.cn_state(None, "chr1", 100), "unavailable")
        label, verdict = L.m_channel_split("chr1", 100, "unavailable", 0.5)
        self.assertEqual(label, "m通道不可用")
        self.assertEqual(verdict["verdict"], "unavailable")

    def test_family_roles_and_full_verification(self):
        units = {u["family"]: u for u in L.process_group(synthetic_group())}
        self.assertTrue(units["1"]["reference_only"])
        self.assertFalse(units["1"]["is_primary_lineage"])
        self.assertEqual(units["1"]["unit_role"], "reference_only_control")
        self.assertTrue(units["2"]["is_primary_lineage"])
        self.assertTrue(units["3"]["is_h3_auxiliary"])
        self.assertFalse(units["3"]["is_primary_lineage"])
        self.assertEqual(units["3"]["unit_role"], "unresolved_H3_auxiliary")
        self.assertTrue(units["4"]["is_h4_auxiliary"])
        self.assertFalse(units["4"]["is_primary_lineage"])
        self.assertEqual(units["4"]["unit_role"], "shared_H4_auxiliary")
        self.assertEqual({u["verification_status"] for u in units.values()}, {"full_pass"})

    def test_complete_longphase_hp_vocabulary_mapping(self):
        expected = {
            None: "none", ".": "none", "1": "1", "1-1": "1", "1-2": "1",
            "2": "2", "2-1": "2", "2-2": "2", "3": "3", "4": "4",
        }
        self.assertEqual({tag: ML.germ_family(tag) for tag in expected}, expected)

    def test_analysis_uses_all_trees_while_display_can_be_capped(self):
        group = synthetic_group()
        group["positions"] = [100, 200]
        group["start"], group["end"], group["n_sSNV"] = 100, 200, 2
        group["populations_by_hp"] = {"2": {"AA": 5}}
        group["subread_groups_by_hp"] = {}
        group["col_coverage_by_hp"] = {"2": {"100": [0, 5], "200": [0, 5]}}
        group["reads_by_hp"] = {"2": 5}
        unit = L.process_group(group)[0]
        self.assertEqual(unit["n_trees"], 2)
        self.assertEqual(unit["analysis_trees_generated"], 2)
        self.assertTrue(unit["analysis_candidate_set_complete"])
        self.assertEqual(unit["n_trees_stored"], 1)
        self.assertFalse(unit["display_trees_complete"])
        self.assertEqual(unit["n_distinct_shapes_exact"], 1)

    def test_read_af_ordering_prefers_higher_ancestor_fraction(self):
        read_af = {0: 0.6, 1: 0.2}
        forward = [("ROOT", "AR"), ("AR", "AA")]
        reverse = [("ROOT", "RA"), ("RA", "AA")]
        self.assertGreater(RA.ordering_score(forward, read_af)[0], RA.ordering_score(reverse, read_af)[0])
        post = RA.posterior([RA.ordering_score(forward, read_af)[0],
                             RA.ordering_score(reverse, read_af)[0]], 0.05)
        self.assertEqual(post.index(max(post)), 0)

    def test_minread4_derivation_preserves_funnel_conservation(self):
        doc = synthetic_input_doc()
        derived = D.transform_doc(doc)
        self.assertEqual(derived["params"]["MINREAD"], 4)
        self.assertTrue(derived["input_funnel"]["check_scope_conservation"])
        self.assertEqual(len(derived["groups"]), 1)
        self.assertEqual(derived["read_tag_census"]["n_regions_with_phase_set_census"], 1)
        self.assertEqual(derived["read_tag_census"]["phase_set_region_counts"], {"none": 1})

    def test_driver_and_region_view_use_primary_denominators(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            ml = td / "ml.json"
            layered = td / "layered.json"
            region = td / "region.json"
            vcf = td / "somatic_pass.vcf"
            ml.write_text(json.dumps(synthetic_input_doc()), encoding="utf-8")
            records = [
                f"chr1\t{i}\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n" for i in range(1, 7)
            ] + [
                f"chrX\t{i}\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n" for i in range(1, 3)
            ]
            vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + "".join(records), encoding="utf-8")
            with mock.patch.dict(os.environ, {"SM_ML": str(ml), "SM_OUT": str(layered)}, clear=False):
                L.main()
            out = json.loads(layered.read_text(encoding="utf-8"))
            l1 = out["L1_ssnv_algorithm"]
            self.assertEqual(l1["n_primary_lineage_units"], 1)
            self.assertEqual(l1["n_reference_only_controls"], 1)
            self.assertEqual(l1["n_unresolved_H3_auxiliary"], 1)
            self.assertEqual(l1["n_shared_H4_auxiliary"], 1)
            self.assertTrue(l1["all_eligible_V1V7_pass"])
            env = os.environ.copy()
            env.update({
                "SM_LAYERED": str(layered),
                "SM_OUT": str(region),
                "SM_SAMPLE": "SYNTHETIC",
                "SM_SOMATIC_VCF": str(vcf),
            })
            proc = subprocess.run(
                [sys.executable, str(HERE / "build_region_view.py")],
                env=env,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            view = json.loads(region.read_text(encoding="utf-8"))
            self.assertEqual(view["regions"][0]["hp_multiplicity"], 1)
            self.assertEqual(view["regions"][0]["n_reference_only_controls"], 1)
            self.assertEqual(view["regions"][0]["n_H3_auxiliary"], 1)
            self.assertEqual(view["regions"][0]["n_H4_auxiliary"], 1)
            self.assertTrue(view["census"]["funnel"]["check_six_branch_conservation"])
            self.assertEqual(view["census"]["U1_sSNV_out_of_scope"], 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
