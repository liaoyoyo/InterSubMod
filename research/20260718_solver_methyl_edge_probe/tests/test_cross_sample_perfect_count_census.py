#!/usr/bin/env python3
"""Tests for the bounded cross-sample perfect-family count census."""

from __future__ import annotations

import copy
import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest


BASE = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = BASE / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def load(path: pathlib.Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


runner = load(
    SCRIPTS / "run_cross_sample_perfect_count_census.py",
    "cross_sample_perfect_count_census_test",
)
counter = load(
    SCRIPTS / "perfect_family_count.py",
    "perfect_family_count_test_for_census",
)


def make_case(
    digest_seed: str,
    *,
    partial_patterns,
    sample="H1437",
):
    structural = {
        "k": 2,
        "structural_alt_universe_mask": 3,
        "full_patterns": [],
        "partial_patterns": list(partial_patterns),
    }
    digest = runner.structural_case_sha256(structural)
    return {
        "detail_order": 1,
        "effective_m": 2,
        "historical_cap_class": "LEVEL_BUDGET_E3",
        "k": 2,
        "reduced_q": 2,
        "sample": sample,
        "selection_reason": digest_seed,
        "structural_input": structural,
        "structural_key_sha256": digest,
        "unit_identity": {
            "chrom": "chr1",
            "start": 1,
            "end": 2,
            "family": "1",
        },
    }


class CaseContractTests(unittest.TestCase):
    def test_exact_and_abstain_rows_are_total_and_never_ranked(self):
        cases = [
            make_case("exact", partial_patterns=["AX", "XA"]),
            make_case(
                "abstain",
                partial_patterns=["AA", "AR", "RA"],
                sample="H2009",
            ),
        ]
        # The two fixtures have different structural hashes despite sharing k/m.
        rows = runner.run_cases(cases, counter=counter, max_m=20)
        census = runner.summarize_rows(
            rows,
            counter=counter,
            expected_cases=2,
        )
        self.assertEqual(census["total_cases"], 2)
        self.assertEqual(census["exact_cases"], 1)
        self.assertEqual(census["abstain_cases"], 1)
        self.assertTrue(census["all_rows_total"])
        self.assertTrue(census["all_elapsed_finite"])
        self.assertFalse(census["ranking_allowed"])
        self.assertFalse(census["formal_speed_claim"])
        self.assertTrue(all(row["ranking_allowed"] is False for row in rows))

    def test_structural_hash_or_dimension_drift_fails_closed(self):
        case = make_case("exact", partial_patterns=["AX", "XA"])
        bad_hash = copy.deepcopy(case)
        bad_hash["structural_key_sha256"] = "0" * 64
        with self.assertRaises(runner.PerfectCountCensusError):
            runner.validate_case(bad_hash)
        bad_m = copy.deepcopy(case)
        bad_m["effective_m"] = 1
        with self.assertRaises(runner.PerfectCountCensusError):
            runner.validate_case(bad_m)

    def test_extra_structural_field_is_incompatible(self):
        case = make_case("exact", partial_patterns=["AX", "XA"])
        case["structural_input"]["unexpected"] = True
        with self.assertRaises(runner.PerfectCountCensusError):
            runner.validate_case(case)


class AuthorityAndReceiptTests(unittest.TestCase):
    def test_input_manifest_byte_sidecar_and_semantic_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            path = root / "manifest.json"
            case = make_case("exact", partial_patterns=["AX", "XA"])
            document = {
                "schema_name": runner.EXPECTED_MANIFEST_SCHEMA,
                "authority": {
                    "status": "AUTHORITATIVE",
                    "claim_scope": runner.CLAIM_SCOPE,
                },
                "checks": {"fixture": True},
                "selection": {},
                "bindings": {},
                "cases": [case],
                "integrity": {},
            }
            keys = [case["structural_key_sha256"]]
            payloads = [case["structural_input"]]
            document["selection"] = {
                "selected_structural_keys": keys,
                "key_list_sha256": runner.semantic_sha256(keys),
                "payload_list_sha256": runner.semantic_sha256(payloads),
            }
            document["integrity"]["manifest_content_sha256"] = (
                runner.document_content_sha256(
                    document,
                    field_name="manifest_content_sha256",
                )
            )
            encoded = (
                json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True)
                + "\n"
            ).encode()
            path.write_bytes(encoded)
            file_sha = __import__("hashlib").sha256(encoded).hexdigest()
            path.with_name("manifest.json.sha256").write_text(
                f"{file_sha}  manifest.json\n"
            )
            observed, binding = runner.verify_authority_manifest(
                path,
                expected_cases=1,
                verify_bound_sources=False,
            )
            self.assertEqual(observed["cases"], [case])
            self.assertEqual(binding["file_sha256"], file_sha)
            self.assertEqual(
                binding["semantic_content_sha256"],
                document["integrity"]["manifest_content_sha256"],
            )

            document["checks"]["fixture"] = False
            path.write_text(json.dumps(document))
            with self.assertRaises(runner.PerfectCountCensusError):
                runner.verify_authority_manifest(
                    path,
                    expected_cases=1,
                    verify_bound_sources=False,
                )

    def test_receipt_content_hash_and_write_once(self):
        receipt = {
            "authority": {"claim_scope": runner.CLAIM_SCOPE},
            "integrity": {},
        }
        digest = runner.document_content_sha256(
            receipt,
            field_name="receipt_content_sha256",
        )
        receipt["integrity"]["receipt_content_sha256"] = digest
        self.assertEqual(
            runner.document_content_sha256(
                receipt,
                field_name="receipt_content_sha256",
            ),
            digest,
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = pathlib.Path(temporary) / "receipt"
            runner.write_once_receipt(receipt, output)
            self.assertTrue((output / "receipt.json").is_file())
            self.assertTrue((output / "receipt.json.sha256").is_file())
            with self.assertRaises(runner.PerfectCountCensusError):
                runner.write_once_receipt(receipt, output)


if __name__ == "__main__":
    unittest.main()
