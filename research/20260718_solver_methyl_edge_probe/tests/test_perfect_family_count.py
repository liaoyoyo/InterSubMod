#!/usr/bin/env python3

from __future__ import annotations

import copy
import hashlib
import json
import pathlib
import sys
import tempfile
import unittest
from unittest import mock


BASE = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = BASE / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from perfect_family_count import (  # noqa: E402
    ABSTAIN_STATUS,
    EXACT_MANDATORY_STATUS,
    EXACT_STATUS,
    EXACT_STATUSES,
    count_manifest_case,
    count_perfect_families,
)
import run_perfect_family_count_validation as validator  # noqa: E402


class PerfectFamilyCountUnitTest(unittest.TestCase):
    def test_two_unconstrained_mutations_have_three_rooted_forests(self):
        result = count_perfect_families(
            full_patterns=[],
            partial_patterns=["AX", "XA"],
            k=2,
            structural_alt_universe_mask=3,
        )
        self.assertEqual(result.status, EXACT_STATUS)
        self.assertEqual(result.objective, 2)
        self.assertEqual(result.perfect_family_count, 3)

    def test_two_co_alt_mutations_have_two_chain_directions(self):
        result = count_perfect_families(
            full_patterns=[],
            partial_patterns=["AA"],
            k=2,
            structural_alt_universe_mask=3,
        )
        self.assertEqual(result.status, EXACT_STATUS)
        self.assertEqual(result.perfect_family_count, 2)

    def test_opposite_ref_constraints_force_sister_events(self):
        result = count_perfect_families(
            full_patterns=[],
            partial_patterns=["AR", "RA"],
            k=2,
            structural_alt_universe_mask=3,
        )
        self.assertEqual(result.status, EXACT_STATUS)
        self.assertEqual(result.perfect_family_count, 1)

    def test_recurrence_required_case_abstains(self):
        result = count_perfect_families(
            full_patterns=[],
            partial_patterns=["AA", "AR", "RA"],
            k=2,
            structural_alt_universe_mask=3,
        )
        self.assertEqual(result.status, ABSTAIN_STATUS)
        self.assertFalse(result.exact)
        self.assertIsNone(result.objective)
        self.assertEqual(result.perfect_family_count, 0)
        self.assertFalse(result.ranking_allowed)

    def test_full_read_mandatory_state_uses_m_minus_d_lower_bound(self):
        result = count_perfect_families(
            full_patterns=["AA"],
            partial_patterns=[],
            k=2,
            structural_alt_universe_mask=3,
        )
        self.assertEqual(result.status, EXACT_MANDATORY_STATUS)
        self.assertEqual(result.mandatory_nonroot_state_count, 1)
        self.assertEqual(result.objective, 1)
        self.assertEqual(result.perfect_family_count, 2)

    def test_two_mandatory_states_use_m_minus_d_lower_bound(self):
        result = count_perfect_families(
            full_patterns=["AR", "AA"],
            partial_patterns=[],
            k=2,
            structural_alt_universe_mask=3,
        )
        self.assertEqual(result.status, EXACT_MANDATORY_STATUS)
        self.assertEqual(result.mandatory_nonroot_state_count, 2)
        self.assertEqual(result.objective, 0)
        self.assertEqual(result.perfect_family_count, 1)

    def test_gapped_active_coordinates_are_compacted_without_changing_count(self):
        result = count_perfect_families(
            full_patterns=[],
            partial_patterns=["AXXA"],
            k=4,
            structural_alt_universe_mask=0b1001,
        )
        self.assertEqual(result.status, EXACT_STATUS)
        self.assertEqual(result.effective_m, 2)
        self.assertEqual(result.objective, 2)
        self.assertEqual(result.perfect_family_count, 2)

    def test_structural_alt_mask_mismatch_is_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "differs"):
            count_perfect_families(
                full_patterns=[],
                partial_patterns=["AX"],
                k=2,
                structural_alt_universe_mask=3,
            )

    def test_zero_active_mutations_are_an_exact_root_only_family(self):
        result = count_perfect_families(
            full_patterns=["RR"],
            partial_patterns=["RX"],
            k=2,
            structural_alt_universe_mask=0,
        )
        self.assertEqual(result.status, EXACT_STATUS)
        self.assertEqual(result.effective_m, 0)
        self.assertEqual(result.objective, 0)
        self.assertEqual(result.perfect_family_count, 1)

    def test_noninteger_public_parameters_fail_closed(self):
        for field, value in (
            ("k", 2.5),
            ("max_m", 2.5),
            ("structural_alt_universe_mask", 1.0),
        ):
            kwargs = {
                "full_patterns": [],
                "partial_patterns": ["AX"],
                "k": 2,
                "structural_alt_universe_mask": 1,
                "max_m": 20,
            }
            kwargs[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                count_perfect_families(**kwargs)


class FrozenStressPanelOracleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        manifest_path = (
            BASE
            / "results"
            / "solver_stress_panel_v1"
            / "authoritative_r3"
            / "manifest.json"
        )
        cls.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        cls.cases = {
            case["case_id"]: case for case in cls.manifest["cases"]
        }
        cls.results = {
            case_id: count_manifest_case(case)
            for case_id, case in cls.cases.items()
        }

    def test_all_46_frozen_keys_are_total(self):
        self.assertEqual(len(self.results), 46)
        self.assertTrue(
            all(
                result.status in EXACT_STATUSES | {ABSTAIN_STATUS}
                for result in self.results.values()
            )
        )

    def test_three_hard_perfect_family_counts(self):
        expected = {
            "k10_m10_a0944cdd1ac8fa9d": 104_640,
            "k11_m11_09b1da787e58efed": 122_281_152,
            "k14_m10_cecad6897b192d47": 27_360,
        }
        for case_id, count in expected.items():
            with self.subTest(case_id=case_id):
                result = self.results[case_id]
                self.assertEqual(result.status, EXACT_STATUS)
                self.assertEqual(result.perfect_family_count, count)
                self.assertEqual(result.objective, result.effective_m)

    def test_e5b_is_explicit_recurrence_abstention(self):
        result = self.results["k10_m9_e5b33e46c7c23c0f"]
        self.assertEqual(result.status, ABSTAIN_STATUS)
        self.assertEqual(result.perfect_family_count, 0)
        self.assertIsNone(result.objective)

    def test_eight_complete_control_oracles(self):
        units_by_case = {}
        for unit in self.manifest["units"]:
            units_by_case.setdefault(unit["case_id"], []).append(unit)
        controls = {
            case_id: units
            for case_id, units in units_by_case.items()
            if any(
                unit["source_classification"] == "complete"
                for unit in units
            )
        }
        self.assertEqual(len(controls), 8)
        exact_matches = 0
        abstentions = 0
        for case_id, units in controls.items():
            expected = {
                unit["source_candidate_vertex_sets"]
                for unit in units
                if unit["source_classification"] == "complete"
            }
            self.assertEqual(len(expected), 1)
            result = self.results[case_id]
            if result.status in EXACT_STATUSES:
                self.assertEqual(
                    result.perfect_family_count, next(iter(expected))
                )
                exact_matches += 1
            else:
                self.assertEqual(case_id, "k6_m6_bd8aa7beafa3f719")
                abstentions += 1
        self.assertEqual(exact_matches, 7)
        self.assertEqual(abstentions, 1)


class PerfectFamilyValidationBindingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel_root = (
            BASE / "results" / "solver_stress_panel_v1"
        )
        cls.pointer_path = cls.panel_root / "AUTHORITATIVE_MANIFEST.json"
        cls.r3_manifest_path = (
            cls.panel_root / "authoritative_r3" / "manifest.json"
        )
        cls.r2_manifest_path = (
            cls.panel_root / "authoritative_r2" / "manifest.json"
        )
        cls.r3_manifest = json.loads(
            cls.r3_manifest_path.read_text(encoding="utf-8")
        )
        cls.cases = {
            case["case_id"]: case for case in cls.r3_manifest["cases"]
        }
        cls.case = cls.cases["k10_m10_a0944cdd1ac8fa9d"]
        cls.real_evidence_dir = (
            cls.panel_root
            / "runs"
            / "optimized_primary_unbounded_q9"
            / "raw"
        )

    @staticmethod
    def _write_with_sidecar(
        path: pathlib.Path,
        document: dict,
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        encoded = (
            json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True)
            + "\n"
        ).encode("utf-8")
        path.write_bytes(encoded)
        sidecar = path.with_name(path.name + ".sha256")
        sidecar.write_text(
            f"{hashlib.sha256(encoded).hexdigest()}  {path.name}\n",
            encoding="ascii",
        )

    @classmethod
    def _manifest_document(
        cls,
        case: dict,
        *,
        authority_status: str | None,
    ) -> dict:
        document = {
            "schema_name": validator.MANIFEST_SCHEMA,
            "schema_version": "1.0.0",
            "checks": {"fixture": True},
            "cases": [copy.deepcopy(case)],
            "units": [{"case_id": case["case_id"]}],
            "integrity": {},
        }
        if authority_status is not None:
            document["authority"] = {"status": authority_status}
        document["integrity"]["manifest_content_sha256"] = (
            validator._manifest_content_sha256(document)
        )
        return document

    @classmethod
    def _evidence_document(
        cls,
        case: dict,
        manifest_content_sha256: str,
    ) -> dict:
        structural = case["structural_input"]
        return {
            "case_id": case["case_id"],
            "structural_key_sha256": case["structural_key_sha256"],
            "k": structural["k"],
            "effective_m": bin(
                structural["structural_alt_universe_mask"]
            ).count("1"),
            "backend": validator.EVIDENCE_BACKEND,
            "schema": validator.EVIDENCE_SCHEMA,
            "objective_status": validator.EVIDENCE_OBJECTIVE_STATUS,
            "objective_certified": True,
            "objective": bin(
                structural["structural_alt_universe_mask"]
            ).count("1"),
            "complete": False,
            "status": "CANDIDATE_SET_INCOMPLETE_DEADLINE",
            "optimal_family_count": 1,
            "manifest_content_sha256": manifest_content_sha256,
        }

    def _write_evidence_fixture(
        self,
        root: pathlib.Path,
        *,
        source_case: dict | None = None,
        evidence_mutation=None,
        unknown_source_digest: bool = False,
    ) -> tuple[pathlib.Path, pathlib.Path, dict, dict]:
        current_case = copy.deepcopy(self.case)
        current_manifest = self._manifest_document(
            current_case,
            authority_status=validator.REQUIRED_AUTHORITY_STATUS,
        )
        current_manifest_path = root / "authoritative_r3" / "manifest.json"
        self._write_with_sidecar(current_manifest_path, current_manifest)

        historical_manifest = self._manifest_document(
            copy.deepcopy(source_case or current_case),
            authority_status=None,
        )
        historical_manifest_path = root / "manifest.json"
        self._write_with_sidecar(historical_manifest_path, historical_manifest)
        source_digest = historical_manifest["integrity"][
            "manifest_content_sha256"
        ]
        if unknown_source_digest:
            source_digest = "f" * 64
        evidence = self._evidence_document(current_case, source_digest)
        if evidence_mutation is not None:
            evidence_mutation(evidence)
        evidence_dir = root / "runs" / "probe" / "raw"
        evidence_dir.mkdir(parents=True)
        evidence_path = (
            evidence_dir
            / f"{current_case['case_id']}.r1.optimized.json"
        )
        evidence_path.write_text(
            json.dumps(evidence, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return evidence_dir, current_manifest_path, current_manifest, current_case

    def test_authority_pointer_resolves_and_binds_r3(self):
        path, manifest, authority = validator._resolve_authority_pointer(
            self.pointer_path
        )
        self.assertEqual(path, self.r3_manifest_path.resolve())
        self.assertEqual(
            manifest["authority"]["status"],
            validator.REQUIRED_AUTHORITY_STATUS,
        )
        self.assertEqual(
            authority["status"],
            validator.REQUIRED_AUTHORITY_STATUS,
        )

    def test_direct_r3_is_accepted_but_direct_r2_is_rejected(self):
        observed = validator._load_frozen_manifest(self.r3_manifest_path)
        self.assertEqual(
            observed["authority"]["status"],
            validator.REQUIRED_AUTHORITY_STATUS,
        )
        with self.assertRaisesRegex(
            validator.ValidationError,
            "authority status mismatch",
        ):
            validator._load_frozen_manifest(self.r2_manifest_path)

    def test_pointer_sidecar_tamper_and_path_escape_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            pointer = {
                "schema": validator.POINTER_SCHEMA,
                "status": validator.REQUIRED_AUTHORITY_STATUS,
                "authoritative_manifest": "../outside/manifest.json",
                "authoritative_manifest_file_sha256": "0" * 64,
                "authoritative_manifest_content_sha256": "0" * 64,
            }
            pointer_path = root / "panel" / "AUTHORITATIVE_MANIFEST.json"
            self._write_with_sidecar(pointer_path, pointer)
            with self.assertRaisesRegex(
                validator.ValidationError,
                "escapes pointer directory",
            ):
                validator._resolve_authority_pointer(pointer_path)

            pointer_path.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(
                validator.ValidationError,
                "byte digest mismatch",
            ):
                validator._resolve_authority_pointer(pointer_path)

    def test_historical_evidence_requires_verified_structural_equivalence(self):
        record = validator._evidence_record(
            self.real_evidence_dir,
            self.case["case_id"],
            expected_case=self.case,
            authoritative_manifest=self.r3_manifest,
            authoritative_manifest_path=self.r3_manifest_path,
        )
        self.assertIsNotNone(record)
        self.assertEqual(
            record["source_manifest"]["binding"],
            "VERIFIED_STRUCTURAL_EQUIVALENCE",
        )
        self.assertEqual(
            record["structural_key_sha256"],
            self.case["structural_key_sha256"],
        )

    def test_stale_or_unresolvable_evidence_manifest_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            evidence_dir, current_path, current, current_case = (
                self._write_evidence_fixture(
                    root,
                    unknown_source_digest=True,
                )
            )
            with self.assertRaisesRegex(
                validator.ValidationError,
                "source manifest is missing",
            ):
                validator._evidence_record(
                    evidence_dir,
                    current_case["case_id"],
                    expected_case=current_case,
                    authoritative_manifest=current,
                    authoritative_manifest_path=current_path,
                )

        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            stale_case = copy.deepcopy(self.cases[
                "k10_m9_e5b33e46c7c23c0f"
            ])
            stale_case["case_id"] = self.case["case_id"]
            evidence_dir, current_path, current, current_case = (
                self._write_evidence_fixture(
                    root,
                    source_case=stale_case,
                )
            )
            with self.assertRaisesRegex(
                validator.ValidationError,
                "not structurally equivalent",
            ):
                validator._evidence_record(
                    evidence_dir,
                    current_case["case_id"],
                    expected_case=current_case,
                    authoritative_manifest=current,
                    authoritative_manifest_path=current_path,
                )

    def test_evidence_identity_schema_objective_and_types_fail_closed(self):
        mutations = {
            "case_id": lambda row: row.__setitem__("case_id", "stale"),
            "structural key": lambda row: row.__setitem__(
                "structural_key_sha256", "0" * 64
            ),
            "k": lambda row: row.__setitem__("k", True),
            "effective_m": lambda row: row.__setitem__("effective_m", True),
            "backend": lambda row: row.__setitem__("backend", "current"),
            "schema": lambda row: row.__setitem__("schema", "stale"),
            "objective_status": lambda row: row.__setitem__(
                "objective_status", "HEURISTIC"
            ),
            "objective_certified": lambda row: row.__setitem__(
                "objective_certified", 1
            ),
            "objective type": lambda row: row.__setitem__("objective", True),
            "complete type": lambda row: row.__setitem__("complete", 0),
        }
        for label, mutation in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                root = pathlib.Path(tmp)
                evidence_dir, current_path, current, current_case = (
                    self._write_evidence_fixture(
                        root,
                        evidence_mutation=mutation,
                    )
                )
                with self.assertRaises(validator.ValidationError):
                    validator._evidence_record(
                        evidence_dir,
                        current_case["case_id"],
                        expected_case=current_case,
                        authoritative_manifest=current,
                        authoritative_manifest_path=current_path,
                    )

    def test_hard_receipt_rejects_certified_but_wrong_objective(self):
        original = validator._evidence_record

        def wrong_objective(*args, **kwargs):
            record = original(*args, **kwargs)
            if record is not None and args[1] == self.case["case_id"]:
                record = copy.deepcopy(record)
                record["objective"] += 1
            return record

        with mock.patch.object(
            validator,
            "_evidence_record",
            side_effect=wrong_objective,
        ):
            with self.assertRaisesRegex(
                validator.ValidationError,
                "objective evidence differs",
            ):
                validator.build_receipt(
                    authority_pointer_path=self.pointer_path,
                    evidence_dir=self.real_evidence_dir,
                    source_path=SCRIPTS / "perfect_family_count.py",
                    runner_path=SCRIPTS
                    / "run_perfect_family_count_validation.py",
                )

    def test_gapped_units_and_manifest_value_types_fail_closed(self):
        document = self._manifest_document(
            copy.deepcopy(self.case),
            authority_status=validator.REQUIRED_AUTHORITY_STATUS,
        )
        document["units"] = []
        with self.assertRaisesRegex(
            validator.ValidationError,
            "case\\(s\\) with no unit",
        ):
            validator._index_manifest_cases(document)

        bad_case = copy.deepcopy(self.case)
        bad_case["structural_input"]["k"] = True
        bad_case["structural_key_sha256"] = validator._semantic_sha256(
            bad_case["structural_input"]
        )
        with self.assertRaisesRegex(
            validator.ValidationError,
            "positive int",
        ):
            validator._validate_case_record(bad_case, label="bad")

        bad_case = copy.deepcopy(self.case)
        bad_case["structural_input"]["partial_patterns"] = "AX"
        bad_case["structural_key_sha256"] = validator._semantic_sha256(
            bad_case["structural_input"]
        )
        with self.assertRaisesRegex(
            validator.ValidationError,
            "must be a list",
        ):
            validator._validate_case_record(bad_case, label="bad")


if __name__ == "__main__":
    unittest.main()
