#!/usr/bin/env python3
"""Tests for the source-bound historical-v5 cross-sample panel exporter."""

from __future__ import annotations

import importlib.util
import pathlib
import sys
import tempfile
import unittest


REPO = pathlib.Path("/big7_disk/liaoyoyo2001/InterSubMod")
SCRIPT = (
    REPO
    / "research/20260718_solver_methyl_edge_probe/scripts"
    / "export_cross_sample_solver_stress_panel.py"
)


def load(path: pathlib.Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


exporter = load(SCRIPT, "cross_sample_panel_exporter_test")


class StructuralPayloadTests(unittest.TestCase):
    def test_filters_below_minread_and_all_unknown(self):
        payload = exporter.structural_payload(
            {
                "AXR": 3,
                "AAR": 5,
                "XXX": 99,
                "RAA": 2,
            },
            k=3,
        )
        self.assertEqual(payload["full_patterns"], ["AAR"])
        self.assertEqual(payload["partial_patterns"], ["AXR"])
        self.assertEqual(payload["structural_alt_universe_mask"], 3)

    def test_rejects_invalid_pattern_or_count(self):
        with self.assertRaises(exporter.CrossSamplePanelError):
            exporter.structural_payload({"AB": 3}, k=2)
        with self.assertRaises(exporter.CrossSamplePanelError):
            exporter.structural_payload({"AA": True}, k=2)
        with self.assertRaises(exporter.CrossSamplePanelError):
            exporter.structural_payload({"AA": -1}, k=2)

    def test_structural_hash_contract_is_order_invariant(self):
        left = exporter.structural_payload(
            {"AXR": 3, "AAR": 5},
            k=3,
        )
        right = exporter.structural_payload(
            {"AAR": 5, "AXR": 3},
            k=3,
        )
        self.assertEqual(left, right)
        self.assertEqual(
            exporter.semantic_sha256(left),
            exporter.semantic_sha256(right),
        )


class SelectionTests(unittest.TestCase):
    @staticmethod
    def row(sample, digest, q, cap="LEVEL_BUDGET_E3"):
        return {
            "sample": sample,
            "k": 5,
            "effective_m": 4,
            "reduced_q": q,
            "historical_cap_class": cap,
            "structural_key_sha256": digest,
        }

    def test_selects_all_q_tail_and_lexicographic_stratum_minimum(self):
        rows = [
            self.row("H1437", "f" * 64, 9),
            self.row("H1437", "b" * 64, 4),
            self.row("H1437", "a" * 64, 4),
            self.row("H1437", "c" * 64, 4, "LEVEL_BUDGET_E4"),
        ]
        selected, summary = exporter.select_cases(rows)
        self.assertEqual(
            {row["structural_key_sha256"] for row in selected},
            {"f" * 64, "a" * 64, "c" * 64},
        )
        self.assertEqual(summary["q_gt_8_memberships"], 1)
        self.assertEqual(summary["strata"], 2)

    def test_rejects_selected_cross_sample_hash_overlap(self):
        shared = "a" * 64
        with self.assertRaises(exporter.CrossSamplePanelError):
            exporter.select_cases(
                [
                    self.row("H1437", shared, 9),
                    self.row("H2009", shared, 9),
                ]
            )

    def test_cap_reason_mapping_is_fail_closed(self):
        self.assertEqual(
            exporter.historical_cap_class(
                "no feasible N within extra_cap=4(pool=31)"
            ),
            "EXTRA_CAP_4",
        )
        self.assertEqual(
            exporter.historical_cap_class(
                "level e=3: C(127,3)=333375 > budget 150000(太密)"
            ),
            "LEVEL_BUDGET_E3",
        )
        self.assertEqual(
            exporter.historical_cap_class(
                "level e=4: C(63,4)=595665 > budget 150000(太密)"
            ),
            "LEVEL_BUDGET_E4",
        )
        with self.assertRaises(exporter.CrossSamplePanelError):
            exporter.historical_cap_class("unexpected")


class ManifestTests(unittest.TestCase):
    def test_manifest_hash_excludes_only_its_own_value(self):
        document = {
            "x": 1,
            "integrity": {"manifest_content_sha256": "old"},
        }
        digest = exporter.manifest_content_sha256(document)
        document["integrity"]["manifest_content_sha256"] = digest
        self.assertEqual(exporter.manifest_content_sha256(document), digest)
        document["x"] = 2
        self.assertNotEqual(exporter.manifest_content_sha256(document), digest)

    def test_write_once_refuses_existing_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = pathlib.Path(temporary) / "authority"
            document = {
                "authority": {
                    "claim_scope": "HISTORICAL_V5_SOLVER_CORE_ONLY"
                }
            }
            exporter.write_once_manifest(document, output)
            self.assertTrue((output / "manifest.json").is_file())
            self.assertTrue((output / "manifest.json.sha256").is_file())
            with self.assertRaises(exporter.CrossSamplePanelError):
                exporter.write_once_manifest(document, output)


if __name__ == "__main__":
    unittest.main()
