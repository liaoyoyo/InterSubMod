#!/usr/bin/env python3
"""drilldown 跨樣本 fail-closed 與 provenance regression tests。"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
for rel in ("scripts/drilldown", "scripts/drilldown/sources",
            "scripts/drilldown/panels", "scripts/drilldown/emit", "scripts"):
    sys.path.insert(0, str(ROOT / rel))

CLI = ROOT / "scripts" / "build_drilldown_dashboard.py"

import capability  # noqa: E402
import ism  # noqa: E402
import lca_ab  # noqa: E402
import lineage_paths  # noqa: E402
import mlhp  # noqa: E402
import selfcheck  # noqa: E402
import strict_edges  # noqa: E402
import topology  # noqa: E402
from bake import panel_byte_totals  # noqa: E402


def _topology_record(sample: str = "HCC1395") -> dict:
    return {
        "schema_name": topology.SCHEMA,
        "schema_version": "1.0.0",
        "sample": sample,
        "region_id": "chr1:100-100:HP1:PS1:U1:B1",
        "chrom": "chr1",
        "active_positions": [100],
        "hp_family": "HP1",
        "phase_set": "PS1",
        "unit_id": "U1",
        "block_id": "B1",
        "active_bit_count": 1,
        "unit_status": "no_active_alt",
    }


class TestLinkageThreshold(unittest.TestCase):
    def test_default_remains_observational(self):
        cap = capability.Capability("x", capability.EXTENSION, "x")
        self.assertTrue(capability.probe_linkage(cap, 1, 100, "legacy"))
        self.assertNotEqual(cap.state, capability.MALFORMED)

    def test_threshold_can_fail_closed(self):
        cap = capability.Capability("x", capability.EXTENSION, "x")
        ok = capability.probe_linkage(
            cap, 1, 100, "sample join", min_rate=0.50,
            fail_state=capability.MALFORMED)
        self.assertFalse(ok)
        self.assertEqual(cap.state, capability.MALFORMED)
        self.assertFalse(cap.usable)

    def test_summary_does_not_claim_scientific_completeness(self):
        reg = capability.Registry()
        cap = capability.make(reg, "x", "x")
        cap.state = capability.PRESENT
        self.assertNotIn("完整", reg.summary_line())
        self.assertIn("可用", reg.summary_line())


class TestSampleIdentityFailClosed(unittest.TestCase):
    def test_token_match_is_not_substring_match(self):
        self.assertTrue(capability.sample_in_path("/run/HCC1395_v2/HCC1395.bam", "HCC1395"))
        self.assertFalse(capability.sample_in_path("/normal/HCC1395BL.bam", "HCC1395"))

    def test_topology_wrong_sample_is_unusable_core(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "HCC1395.topology.jsonl"
            path.write_text(json.dumps(_topology_record("HCC1395")) + "\n", encoding="utf-8")
            reg = capability.Registry()
            cap = topology.load(reg, path, expected_sample="COLO829")
            self.assertEqual(cap.state, capability.MALFORMED)
            with self.assertRaises(capability.Refuse):
                reg.require_core()

    def test_mlhp_wrong_sample_is_unusable(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "HCC1395.exact_ps_mlhp.json"
            path.write_text(json.dumps({
                "sample": "HCC1395",
                "groups": [{
                    "region_id": "r1", "positions": [100, 200],
                    "populations_by_hp": {}, "subread_groups_by_hp": {},
                }],
            }), encoding="utf-8")
            reg = capability.Registry()
            cap = mlhp.load(reg, path, {"r1"}, expected_sample="COLO829")
            self.assertEqual(cap.state, capability.MALFORMED)
            self.assertFalse(cap.usable)

    def test_ism_wrong_tumor_input_is_unusable(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            chrom = root / "chr1"
            chrom.mkdir()
            (chrom / "significance_summary.csv").write_text(
                "Chr,Pos,HPMergedP,HPMergedDelta\nchr1,100,0.5,0.1\n", encoding="utf-8")
            (chrom / "run_params.json").write_text(json.dumps({
                "schema_name": "intersubmod.run_params",
                "io": {"tumor_bam_path": "/runs/HCC1395_v2/HCC1395.bam"},
                "region": {"window_size_bp": 1000},
            }), encoding="utf-8")
            l1 = {"n": 1, "chroms": ["chr1"], "chrom": [0], "dpos": [100]}
            reg = capability.Registry()
            cap = ism.load(reg, root, l1, ["chr1"], expected_sample="COLO829")
            self.assertEqual(cap.state, capability.MALFORMED)
            self.assertFalse(cap.usable)

    def test_lineage_directory_cannot_lend_other_sample_file(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "HCC1395.unit_lineage_paths.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as fh:
                fh.write("region_id\tvertex\tvertex_label\tis_hidden\tdepth\tlineage_path\n")
            reg = capability.Registry()
            cap = lineage_paths.load(reg, td, {"r1"}, expected_sample="COLO829")
            self.assertEqual(cap.state, capability.MALFORMED)
            self.assertFalse(cap.usable)

    def test_lca_directory_cannot_lend_other_sample_receipts(self):
        with tempfile.TemporaryDirectory() as td:
            pre, post = Path(td) / "pre", Path(td) / "post"
            pre.mkdir()
            post.mkdir()
            doc = {"topology": "/x/HCC1395/HCC1395.topology.jsonl", "stats": {"lv_written": 1}}
            for root in (pre, post):
                (root / "HCC1395.chr1.tag_bam.receipt.json").write_text(
                    json.dumps(doc), encoding="utf-8")
            reg = capability.Registry()
            cap = lca_ab.load(reg, pre, post, expected_sample="COLO829")
            self.assertEqual(cap.state, capability.MALFORMED)
            self.assertFalse(cap.usable)

    def test_strict_edge_root_cannot_lend_other_sample_files(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / "HCC1395" / "chromosomes"
            root.mkdir(parents=True)
            reg = capability.Registry()
            cap = strict_edges.load(reg, root, [], expected_sample="COLO829")
            self.assertEqual(cap.state, capability.MALFORMED)
            self.assertFalse(cap.usable)


class TestLcaPerChromContract(unittest.TestCase):
    @staticmethod
    def _receipt(sample: str, chrom: str, lv: int, reads: int) -> dict:
        return {
            "topology": f"/x/{sample}/{sample}.topology.jsonl",
            "in_bam": f"/x/{sample}/{sample}.bam",
            "assignments": f"/x/{sample}/{sample}.assign.tsv.gz",
            "out_bam": f"/x/{sample}/{sample}.{chrom}.bam",
            "region": chrom,
            "stats": {"lv_written": lv, "reads_total": reads},
        }

    def test_aggregate_cancellation_does_not_hide_per_chrom_mismatch(self):
        with tempfile.TemporaryDirectory() as td:
            pre, post = Path(td) / "pre", Path(td) / "post"
            pre.mkdir()
            post.mkdir()
            # reads_total aggregate 都是 30，但 chr1/chr2 各自不同。
            for chrom, a, b in (("chr1", 10, 11), ("chr2", 20, 19)):
                (pre / f"HCC1395.{chrom}.tag_bam.receipt.json").write_text(
                    json.dumps(self._receipt("HCC1395", chrom, 1, a)), encoding="utf-8")
                (post / f"HCC1395.{chrom}.tag_bam.receipt.json").write_text(
                    json.dumps(self._receipt("HCC1395", chrom, 2, b)), encoding="utf-8")
            reg = capability.Registry()
            cap = lca_ab.load(reg, pre, post, expected_sample="HCC1395")
            self.assertEqual(cap.state, capability.PARTIAL)
            self.assertEqual(set(cap.payload["per_chrom_unexpected"]), {"chr1", "chr2"})
            self.assertEqual(cap.payload["pre_files"], 2)
            self.assertEqual(cap.payload["post_files"], 2)
            self.assertEqual(cap.payload["shared_files"], 2)

    def test_changed_input_identity_blocks_causal_lca_claim(self):
        with tempfile.TemporaryDirectory() as td:
            pre, post = Path(td) / "pre", Path(td) / "post"
            pre.mkdir()
            post.mkdir()
            a = self._receipt("HCC1395", "chr1", 1, 10)
            b = self._receipt("HCC1395", "chr1", 2, 10)
            b["in_bam"] = "/other/HCC1395/HCC1395_tagged.bam"
            a["threads"] = 12
            b["threads"] = 16
            (pre / "HCC1395.chr1.tag_bam.receipt.json").write_text(
                json.dumps(a), encoding="utf-8")
            (post / "HCC1395.chr1.tag_bam.receipt.json").write_text(
                json.dumps(b), encoding="utf-8")
            reg = capability.Registry()
            cap = lca_ab.load(reg, pre, post, expected_sample="HCC1395")
            self.assertEqual(cap.state, capability.PARTIAL)
            self.assertEqual(cap.payload["per_chrom_input_mismatch"],
                             {"chr1": ["in_bam", "threads"]})
            topo = capability.make(reg, "topology", "topology", tier=capability.CORE)
            topo.state = capability.PRESENT
            topo.payload = {"regions": [],
                            "l1": {"n": 0, "chrom": [], "csrVal": []}}
            c11 = next(c for c in selfcheck.run(reg)["checks"] if c["id"] == "C11")
            self.assertEqual(c11["status"], selfcheck.FAIL)


class TestPanelReceiptBytes(unittest.TestCase):
    def test_base_and_tumor_are_both_counted(self):
        self.assertEqual(panel_byte_totals({"bytes": 100, "tumorOnly": {"bytes": 60}}),
                         (100, 60, 160))
        self.assertEqual(panel_byte_totals({"bytes": 100}), (100, 0, 100))


class TestCliAndReceipt(unittest.TestCase):
    def test_unimplemented_fig_modes_are_rejected(self):
        result = subprocess.run(
            [sys.executable, str(CLI), "--sample", "HCC1395", "--figs-mode", "link",
             "--probe-only", "--topology", "/missing"],
            capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 2)
        self.assertIn("invalid choice", result.stderr)

    def test_minimal_build_emits_versioned_reproducible_receipt(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            topo = root / "HCC1395.topology.jsonl"
            topo.write_text(json.dumps(_topology_record()) + "\n", encoding="utf-8")
            out = root / "out"
            missing = str(root / "missing")
            result = subprocess.run([
                sys.executable, str(CLI), "--sample", "HCC1395",
                "--topology", str(topo), "--topology-receipt", missing,
                "--mlhp", missing, "--ism-root", missing, "--chrom-root", missing,
                "--lca-pre", missing, "--lca-post", missing, "--paths-dir", missing,
                "--out", str(out), "--figs-mode", "copy",
            ], capture_output=True, text=True, timeout=120, cwd=ROOT)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            receipt = json.loads((out / "receipt.json").read_text(encoding="utf-8"))
            self.assertEqual(receipt["schema_name"], "intersubmod.drilldown_dashboard_receipt")
            self.assertEqual(receipt["schema_version"], "2.0.0")
            self.assertEqual(receipt["claim_ceiling"], "OBSERVATION_ONLY_NO_TRUTH_SET")
            self.assertIn(receipt["validation_status"], {
                "INPUT_IDENTITY_OR_SCHEMA_FAILED", "SELF_CHECK_FAILED",
                "SOURCE_STATE_BLOCKED", "PROVENANCE_INCOMPLETE",
                "SELF_CHECK_INCOMPLETE", "SELF_CHECK_PASS_INTERNAL"})
            self.assertFalse(receipt["validation"]["scientifically_validated"])
            self.assertFalse(receipt["validation"]["gates"]["truth_set_validation"])
            self.assertFalse(receipt["validation"]["publishable_internal_qa"])
            self.assertTrue(receipt["generator"]["git_commit"])
            self.assertEqual(len(receipt["generator"]["script_sha256"]), 64)
            self.assertIn("--sample HCC1395", receipt["command"])
            self.assertGreater(receipt["output_inventory"]["total_count"], 0)
            self.assertIn("data_js", receipt["output_inventory"]["by_type"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
