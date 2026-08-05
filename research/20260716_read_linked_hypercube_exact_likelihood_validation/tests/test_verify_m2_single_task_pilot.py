#!/usr/bin/env python3
"""Synthetic, no-BAM tests for the single-task M2 pilot gate."""

from __future__ import annotations

import csv
import gzip
import json
import os
import pathlib
import re
import stat
import subprocess
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build_m2_patterns_and_rank as ranker  # noqa: E402
import run_full_m2_extraction as resource_producer  # noqa: E402
import verify_full_m2_receipts as independent  # noqa: E402
import verify_m2_single_task_pilot as gate  # noqa: E402


DATASET = "HCC1954"
CHROM = "chr22"


def write_gzip_tsv(path: pathlib.Path, fields, rows) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def identity(path: pathlib.Path) -> dict:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": independent.sha256_path(path),
    }


def rewrite_receipt(path: pathlib.Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    path.with_name(f"{path.name}.sha256").write_text(
        f"{independent.sha256_path(path)}  {path.name}\n", encoding="ascii"
    )


def rewrite_immutable_receipt(path: pathlib.Path, payload: dict) -> None:
    sidecar = path.with_name(f"{path.name}.sha256")
    path.chmod(0o644)
    sidecar.chmod(0o644)
    rewrite_receipt(path, payload)
    path.chmod(0o444)
    sidecar.chmod(0o444)


def write_time(path: pathlib.Path, elapsed: str = "0:01.50", rss_kib: int = 100_000) -> None:
    path.write_text(
        "\n".join((
            "Command being timed: synthetic-no-bam-fixture",
            f"Elapsed (wall clock) time (h:mm:ss or m:ss): {elapsed}",
            f"Maximum resident set size (kbytes): {rss_kib}",
            "Exit status: 0",
            "",
        )),
        encoding="utf-8",
    )


def write_synthetic_frozen_release(
    path: pathlib.Path, canonical_manifest: pathlib.Path
) -> dict:
    relative = canonical_manifest.resolve().relative_to(path.parent.resolve())
    payload = {
        "schema_name": independent.RELEASE_SCHEMA_NAME,
        "schema_version": independent.RELEASE_SCHEMA_VERSION,
        "authority_mode": "CANONICAL_V5_FROZEN",
        "validation_evidence_eligible": True,
        "all_pass": True,
        "canonical_manifest": {
            "immutable_copy": {
                "path": relative.as_posix(),
                "sha256": independent.sha256_path(canonical_manifest),
            },
        },
        "receipt_integrity": resource_producer._receipt_integrity(path),
    }
    resource_producer.write_immutable_json_exclusive(path, payload)
    sidecar = path.with_name(f"{path.name}.sha256")
    return {
        "path": str(path.resolve()),
        "sha256": independent.sha256_path(path),
        "semantic_sha256": independent.semantic_json_sha256(payload),
        "sidecar": {
            "path": str(sidecar.resolve()),
            "sha256": independent.sha256_path(sidecar),
        },
    }


def make_pilot(
    root: pathlib.Path,
    *,
    ranking_dirname: str = gate.RANKING_DIRNAME,
    bootstrap_replicates: int = 0,
) -> pathlib.Path:
    pilot = root / "pilot"
    extraction = pilot / gate.EXTRACTION_DIRNAME
    ranking = pilot / ranking_dirname
    extraction.mkdir(parents=True)

    calls = extraction / f"{DATASET}.{CHROM}.molecule_sparse_calls.tsv.gz"
    sites = extraction / f"{DATASET}.{CHROM}.site_catalog.tsv.gz"
    components = extraction / f"{DATASET}.{CHROM}.components.tsv.gz"
    membership = extraction / f"{DATASET}.{CHROM}.site_component_membership.tsv.gz"
    write_gzip_tsv(
        sites,
        ("site_index", "chrom", "pos1", "ref", "alt"),
        [
            {"site_index": index, "chrom": CHROM, "pos1": 10 + index, "ref": "C", "alt": "T"}
            for index in range(2)
        ],
    )
    write_gzip_tsv(
        components,
        (
            "dataset", "chrom", "component_basis", "phase_set", "phase_set_status",
            "inference_role", "threshold", "component_id", "start1", "end1", "k",
        ),
        [{
            "dataset": DATASET, "chrom": CHROM, "component_basis": "PS_HP1",
            "phase_set": "100", "phase_set_status": "KNOWN_PS_PRIMARY",
            "inference_role": "PRIMARY_PS_AWARE", "threshold": 3,
            "component_id": "SYNTHETIC_C", "start1": 10, "end1": 11, "k": 2,
        }],
    )
    write_gzip_tsv(
        membership,
        (
            "dataset", "chrom", "component_basis", "phase_set", "phase_set_status",
            "inference_role", "threshold", "site_index", "pos1", "component_id",
        ),
        [
            {
                "dataset": DATASET, "chrom": CHROM, "component_basis": "PS_HP1",
                "phase_set": "100", "phase_set_status": "KNOWN_PS_PRIMARY",
                "inference_role": "PRIMARY_PS_AWARE", "threshold": 3,
                "site_index": index, "pos1": 10 + index, "component_id": "SYNTHETIC_C",
            }
            for index in range(2)
        ],
    )
    # AA alone has two minimum compatible vertex sets, so B20 exercises the
    # real conditional-bootstrap branch instead of only the V=1 not-applicable path.
    call_rows = [
        {
            "molecule_id": f"AA{repeat}", "hp_family": "1", "phase_set": "100",
            "site_indices": "0,1", "call_codes": "AA", "base_qualities": "30,30",
        }
        for repeat in range(3)
    ]
    write_gzip_tsv(
        calls,
        ("molecule_id", "hp_family", "phase_set", "site_indices", "call_codes", "base_qualities"),
        call_rows,
    )

    manifest = pilot / "synthetic_manifest.json"
    extractor = pilot / "synthetic_extractor.py"
    manifest.write_text('{"fixture":"no BAM is opened"}\n', encoding="utf-8")
    extractor.write_text("# synthetic producer identity\n", encoding="utf-8")
    extraction_receipt = {
        "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
        "schema_version": "1.2.0",
        "scope": {"dataset": DATASET, "chrom": CHROM, "n_sSNV": 2},
        "provenance": {
            "manifest": {"path": str(manifest.resolve()), "sha256": independent.sha256_path(manifest)},
            "extractor": {"path": str(extractor.resolve()), "sha256": independent.sha256_path(extractor)},
        },
        "parameters": dict(gate.EXPECTED_EXTRACTION_PARAMETERS),
        "counts": {
            "raw_overlapping_alignments": 3,
            "alignment_class_primary": 3,
            "alignment_class_secondary": 0,
            "alignment_class_supplementary": 0,
            "alignment_class_unmapped": 0,
            "excluded_by_flag": 0,
            "mapq_rejected_after_flag": 0,
            "canonical_eligible_alignments": 3,
            "molecule_sparse_rows_written": 3,
            "unique_molecule_ids": 3,
            "sidecar_exact_matches": 3,
            "sidecar_missing": 0,
            "site_call_rows_sparse": 6,
            "fixed_ra_calls": 6,
            "alt_calls": 6,
        },
        "checks": {key: True for key in gate.REQUIRED_EXTRACTION_CHECKS},
        "all_pass": True,
        "outputs": {path.name: identity(path) for path in (calls, sites, components, membership)},
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": "receipt.json.sha256",
            "covers": "receipt.json",
        },
    }
    rewrite_receipt(extraction / "receipt.json", extraction_receipt)
    result = ranker.run(
        extraction,
        ranking,
        thresholds={1, 2, 3, 5},
        component_bases={"PS_HP1", "PS_HP2"},
        families=("1", "2"),
        structural_exact_pattern_minreads=(1, 2, 3, 5),
        primary_structural_exact_pattern_minread=3,
        exact_k_max=12,
        max_vertex_sets=256,
        solver_time_limit_seconds=30.0,
        fixed_error_grid=(0.005, 0.01, 0.02, 0.05),
        minimum_bq_error_rate=0.000001,
        maximum_bq_error_rate=0.25,
        conditional_candidate_ranking_bootstrap_replicates=bootstrap_replicates,
        conditional_candidate_ranking_bootstrap_seed=20260716,
        tie_tolerance=0.000001,
    )
    if not result["all_pass"]:
        raise AssertionError("synthetic production ranker fixture did not pass")
    release_path = pilot / "synthetic_release_manifest.json"
    release_identity = write_synthetic_frozen_release(release_path, manifest)
    gate_dir = pilot / "resource_gates"
    for stage, label, output_root in (
        ("extraction", "extraction", extraction),
        ("ranking", ranking_dirname, ranking),
    ):
        observed = output_root.stat()
        target = {
            "task_type": "B_COMPREHENSIVE_VALIDATION_RELEASE_PILOT",
            "dataset": DATASET,
            "chrom": CHROM,
            "gate_label": label,
            "output_root": {
                "path": str(output_root.resolve()), "st_dev": int(observed.st_dev),
                "st_ino": int(observed.st_ino),
            },
            "manifest": {
                "path": str(manifest.resolve()),
                "sha256": independent.sha256_path(manifest),
            },
            "release_manifest": release_identity,
        }
        resource_producer.create_resource_gate_receipt(
            output_root, stage=stage, gate_scope="pilot", target=target,
            producer_source={
                "path": str(resource_producer.RUNNER.resolve()),
                "sha256": independent.sha256_path(resource_producer.RUNNER),
            },
            conflicts={"process_count": 0, "root_count": 0, "representatives": []},
            receipt_path=gate_dir / f"{label}.json",
        )
    write_time(pilot / "extraction.time.txt")
    write_time(pilot / f"{ranking_dirname}.time.txt")
    return pilot


class PositiveAndTamperTests(unittest.TestCase):
    def test_valid_pilot_is_go_and_writes_authenticated_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            receipt = gate.verify_pilot(pilot, DATASET, CHROM)
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(receipt["release_gate"]["verdict"], "GO")
            self.assertFalse(receipt["verification_independence"]["reads_bam"])
            self.assertGreater(
                receipt["profile_likelihood_independent_recomputation"]["n_candidates"], 0
            )
            self.assertEqual(
                receipt["unit_audit"]["bootstrap_status_counts"],
                {"NOT_RUN_REQUESTED_0": 1},
            )
            self.assertEqual(
                set(receipt["resource_gate_receipts"]),
                {"extraction", gate.RANKING_DIRNAME},
            )
            self.assertEqual(
                receipt["resource_gate_receipts"]["extraction"]["required_reserve_bytes"],
                300 * 1024**3,
            )
            self.assertTrue(
                receipt["resource_gate_cross_binding"][
                    "both_gate_manifests_equal_extraction_child_provenance"
                ]
            )
            self.assertEqual(
                receipt["resource_gate_receipts"]["extraction"][
                    "release_manifest_identity"
                ],
                receipt["resource_gate_receipts"][gate.RANKING_DIRNAME][
                    "release_manifest_identity"
                ],
            )
            output = pilot / "pilot_gate_verification_receipt.json"
            sha = gate.write_receipt(output, receipt)
            loaded, authenticated_sha = independent.authenticate_receipt(
                output,
                "intersubmod.m2_single_task_pilot_independent_verification",
                "1.0.0",
            )
            self.assertEqual(sha, authenticated_sha)
            self.assertTrue(loaded["all_pass"])
            self.assertEqual(stat.S_IMODE(os.lstat(output).st_mode), 0o444)
            self.assertEqual(os.lstat(output).st_nlink, 1)
            sidecar = output.with_name(f"{output.name}.sha256")
            self.assertEqual(stat.S_IMODE(os.lstat(sidecar).st_mode), 0o444)
            self.assertEqual(os.lstat(sidecar).st_nlink, 1)

    def test_receipt_writer_rejects_existing_and_partial_preseed_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            complete = root / "complete.json"
            gate.write_receipt(complete, {"all_pass": True})
            with self.assertRaisesRegex(FileExistsError, "preseeded"):
                gate.write_receipt(complete, {"all_pass": False})

            receipt_only = root / "receipt_only.json"
            receipt_only.write_text("prior receipt bytes\n", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "preseeded"):
                gate.write_receipt(receipt_only, {"all_pass": True})
            self.assertEqual(
                receipt_only.read_text(encoding="utf-8"), "prior receipt bytes\n"
            )
            self.assertFalse(
                receipt_only.with_name(f"{receipt_only.name}.sha256").exists()
            )

            sidecar_only = root / "sidecar_only.json"
            preseeded_sidecar = sidecar_only.with_name(f"{sidecar_only.name}.sha256")
            preseeded_sidecar.write_text("prior sidecar bytes\n", encoding="ascii")
            with self.assertRaisesRegex(FileExistsError, "preseeded"):
                gate.write_receipt(sidecar_only, {"all_pass": True})
            self.assertFalse(sidecar_only.exists())
            self.assertEqual(
                preseeded_sidecar.read_text(encoding="ascii"), "prior sidecar bytes\n"
            )
    def test_pilot_resource_gate_tamper_and_path_swap_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            extraction_gate = pilot / "resource_gates/extraction.json"
            ranking_gate = pilot / f"resource_gates/{gate.RANKING_DIRNAME}.json"
            with self.assertRaisesRegex(gate.PilotVerificationError, "path swap"):
                gate.verify_pilot(
                    pilot,
                    DATASET,
                    CHROM,
                    ranking_resource_gate_receipt=extraction_gate,
                )

            original = json.loads(ranking_gate.read_text(encoding="utf-8"))
            tampered = json.loads(json.dumps(original))
            tampered["process_snapshot"]["semantic_sha256"] = "0" * 64
            rewrite_immutable_receipt(ranking_gate, tampered)
            with self.assertRaisesRegex(gate.PilotVerificationError, "gate_id"):
                gate.verify_pilot(pilot, DATASET, CHROM)
            rewrite_immutable_receipt(ranking_gate, original)

            producer_tampered = json.loads(json.dumps(original))
            producer_tampered["producer_source"]["sha256"] = "0" * 64
            producer_tampered["gate_id"] = independent.semantic_json_sha256({
                key: value for key, value in producer_tampered.items()
                if key not in {"gate_id", "receipt_integrity"}
            })
            rewrite_immutable_receipt(ranking_gate, producer_tampered)
            with self.assertRaisesRegex(gate.PilotVerificationError, "source SHA"):
                gate.verify_pilot(pilot, DATASET, CHROM)

    def test_self_consistent_alternate_manifest_or_release_cannot_cross_bind(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            ranking_gate_path = pilot / f"resource_gates/{gate.RANKING_DIRNAME}.json"
            original_gate = json.loads(ranking_gate_path.read_text(encoding="utf-8"))

            alternate_manifest = pilot / "alternate_manifest.json"
            alternate_manifest.write_text(
                '{"fixture":"different but internally valid"}\n', encoding="utf-8"
            )
            alternate_manifest.chmod(0o444)
            alternate_release = write_synthetic_frozen_release(
                pilot / "alternate_release_manifest.json", alternate_manifest
            )
            alternate_gate = json.loads(json.dumps(original_gate))
            alternate_gate["target"]["manifest"] = {
                "path": str(alternate_manifest.resolve()),
                "sha256": independent.sha256_path(alternate_manifest),
            }
            alternate_gate["target"]["release_manifest"] = alternate_release
            alternate_gate["gate_id"] = independent.semantic_json_sha256({
                key: value for key, value in alternate_gate.items()
                if key not in {"gate_id", "receipt_integrity"}
            })
            rewrite_immutable_receipt(ranking_gate_path, alternate_gate)
            self_consistent = gate.verify_pilot_resource_gate(
                ranking_gate_path,
                pilot_root=pilot,
                output_root=pilot / gate.RANKING_DIRNAME,
                dataset=DATASET,
                chrom=CHROM,
                stage="ranking",
                gate_label=gate.RANKING_DIRNAME,
            )
            self.assertEqual(
                self_consistent["canonical_manifest_identity"]["sha256"],
                independent.sha256_path(alternate_manifest),
            )
            with self.assertRaisesRegex(
                gate.PilotVerificationError, "do not share the exact canonical manifest"
            ):
                gate.verify_pilot(pilot, DATASET, CHROM)

            rewrite_immutable_receipt(ranking_gate_path, original_gate)
            original_manifest = pilot / "synthetic_manifest.json"
            second_release = write_synthetic_frozen_release(
                pilot / "second_release_same_manifest.json", original_manifest
            )
            release_only_gate = json.loads(json.dumps(original_gate))
            release_only_gate["target"]["release_manifest"] = second_release
            release_only_gate["gate_id"] = independent.semantic_json_sha256({
                key: value for key, value in release_only_gate.items()
                if key not in {"gate_id", "receipt_integrity"}
            })
            rewrite_immutable_receipt(ranking_gate_path, release_only_gate)
            gate.verify_pilot_resource_gate(
                ranking_gate_path,
                pilot_root=pilot,
                output_root=pilot / gate.RANKING_DIRNAME,
                dataset=DATASET,
                chrom=CHROM,
                stage="ranking",
                gate_label=gate.RANKING_DIRNAME,
            )
            with self.assertRaisesRegex(
                gate.PilotVerificationError, "exact frozen release binding"
            ):
                gate.verify_pilot(pilot, DATASET, CHROM)

    def test_extraction_output_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            calls = next((pilot / gate.EXTRACTION_DIRNAME).glob("*.molecule_sparse_calls.tsv.gz"))
            with calls.open("ab") as handle:
                handle.write(b"tamper")
            with self.assertRaisesRegex(gate.PilotVerificationError, "size mismatch|SHA-256 mismatch"):
                gate.verify_pilot(pilot, DATASET, CHROM)

    def test_receipt_sidecar_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            rank_receipt = pilot / gate.RANKING_DIRNAME / "receipt.json"
            with rank_receipt.open("a", encoding="utf-8") as handle:
                handle.write(" \n")
            with self.assertRaisesRegex(gate.PilotVerificationError, "sidecar mismatch"):
                gate.verify_pilot(pilot, DATASET, CHROM)

    def test_parameter_tamper_with_fresh_sidecar_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            path = pilot / gate.RANKING_DIRNAME / "receipt.json"
            receipt = json.loads(path.read_text(encoding="utf-8"))
            receipt["parameters"]["exact_k_max"] = 13
            rewrite_receipt(path, receipt)
            with self.assertRaisesRegex(gate.PilotVerificationError, "exact_k_max"):
                gate.verify_pilot(pilot, DATASET, CHROM)

    def test_profile_likelihood_tamper_is_rejected_after_rebinding_output_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            ranking_dir = pilot / gate.RANKING_DIRNAME
            candidate = ranking_dir / independent.CANDIDATE_SOURCE_NAME
            with gzip.open(candidate, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields, rows = tuple(reader.fieldnames or ()), list(reader)
            primary_row = next(
                row for row in rows if int(row["structural_exact_pattern_minread"]) == 3
            )
            primary_row["primary_log_likelihood"] = str(
                float(primary_row["primary_log_likelihood"]) + 0.25
            )
            write_gzip_tsv(candidate, fields, rows)
            receipt_path = ranking_dir / "receipt.json"
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["outputs"][candidate.name].update(identity(candidate))
            rewrite_receipt(receipt_path, receipt)
            with self.assertRaisesRegex(gate.PilotVerificationError, "profile likelihood"):
                gate.verify_pilot(pilot, DATASET, CHROM)

    def test_five_hour_wall_is_probe_and_nine_hour_wall_is_no_go(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(pathlib.Path(tmp))
            write_time(pilot / "ranking_bootstrap0.time.txt", elapsed="5:00:00")
            probe = gate.verify_pilot(pilot, DATASET, CHROM)
            self.assertEqual(probe["release_gate"]["verdict"], "PROBE")
            self.assertFalse(probe["all_pass"])
            write_time(pilot / "ranking_bootstrap0.time.txt", elapsed="9:00:00")
            no_go = gate.verify_pilot(pilot, DATASET, CHROM)
            self.assertEqual(no_go["release_gate"]["verdict"], "NO_GO")
            self.assertFalse(no_go["all_pass"])

    def test_bootstrap20_directory_and_status_contract_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(
                pathlib.Path(tmp),
                ranking_dirname="ranking_bootstrap20",
                bootstrap_replicates=20,
            )
            receipt = gate.verify_pilot(
                pilot,
                DATASET,
                CHROM,
                ranking_dirname="ranking_bootstrap20",
                expected_bootstrap_replicates=20,
            )
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(receipt["release_gate"]["verdict"], "GO")
            self.assertEqual(receipt["unit_audit"]["expected_bootstrap_replicates"], 20)
            self.assertEqual(receipt["unit_audit"]["bootstrap_status_counts"], {"COMPLETE": 1})
            self.assertEqual(receipt["unit_audit"]["bootstrap_evaluated_units"], 1)
            self.assertIn("ranking_bootstrap20", receipt["process_resources"])

    def test_bootstrap20_unit_replicate_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            pilot = make_pilot(
                pathlib.Path(tmp),
                ranking_dirname="ranking_bootstrap20",
                bootstrap_replicates=20,
            )
            ranking_dir = pilot / "ranking_bootstrap20"
            units = ranking_dir / gate.UNIT_SOURCE_NAME
            with gzip.open(units, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields, rows = tuple(reader.fieldnames or ()), list(reader)
            primary = next(row for row in rows if row["structural_minread_role"] == "PRIMARY")
            primary["conditional_candidate_ranking_bootstrap_replicates"] = "0"
            write_gzip_tsv(units, fields, rows)
            receipt_path = ranking_dir / "receipt.json"
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["outputs"][units.name].update(identity(units))
            rewrite_receipt(receipt_path, receipt)
            with self.assertRaisesRegex(
                gate.PilotVerificationError,
                "requested bootstrap unit has wrong status/replicate/frequency fields",
            ):
                gate.verify_pilot(
                    pilot,
                    DATASET,
                    CHROM,
                    ranking_dirname="ranking_bootstrap20",
                    expected_bootstrap_replicates=20,
                )


class RunbookOperationalClosureTests(unittest.TestCase):
    RUNBOOK = ROOT / "20260719_M2_exact_preserving資源閘門修補正式執行Runbook_04.md"

    @classmethod
    def function_body(cls, name: str) -> str:
        text = cls.RUNBOOK.read_text(encoding="utf-8")
        match = re.search(
            rf"(?ms)^{re.escape(name)}\(\) \{{\n(?P<body>.*?)^\}}$", text
        )
        if match is None:
            raise AssertionError(f"missing runbook function: {name}")
        return match.group("body")

    def test_direct_receipt_sealing_precedes_authentication_and_resume_reauthenticates(self):
        helper = self.function_body("seal_direct_receipt")
        required_fragments = (
            'test -f "$target" && test ! -L "$target"',
            'test -f "$sidecar" && test ! -L "$sidecar"',
            'test "$(stat -c %h "$target")" -eq 1',
            'test "$(stat -c %h "$sidecar")" -eq 1',
            "grep -Eq '^[0-9a-f]{64}  [^[:space:]]+$'",
            'test "$covered" = "$(basename "$target")"',
            'sha256sum --strict --status -c',
            'chmod 0444 -- "$target" "$sidecar"',
            'authenticate_sidecar "$target"',
        )
        for fragment in required_fragments:
            self.assertIn(fragment, helper)
        self.assertLess(
            helper.index('sha256sum --strict --status -c'),
            helper.index('chmod 0444 -- "$target" "$sidecar"'),
        )
        self.assertLess(
            helper.index('chmod 0444 -- "$target" "$sidecar"'),
            helper.index('authenticate_sidecar "$target"'),
        )

        extraction = self.function_body("run_pilot_extraction")
        extraction_seal = extraction.index('seal_direct_receipt "$receipt"')
        self.assertLess(extraction.index('expect_rc "$rc" 0'), extraction_seal)
        self.assertGreater(
            extraction.find('authenticate_sidecar "$receipt"', extraction_seal),
            extraction_seal,
        )

        ranking = self.function_body("run_pilot_ranking")
        ranking_seal = ranking.index('seal_direct_receipt "$receipt"')
        self.assertLess(ranking.index('expect_rc "$rc" 0'), ranking_seal)
        self.assertGreater(
            ranking.find('authenticate_sidecar "$receipt"', ranking_seal),
            ranking_seal,
        )

        runbook = self.RUNBOOK.read_text(encoding="utf-8")
        self.assertEqual(runbook.count('seal_direct_receipt "$PRE"'), 1)
        self.assertEqual(runbook.count('seal_direct_receipt "$POST"'), 1)
        resume = self.function_body("validate_pilot_resume_state")
        for path in (
            "$pilot/extraction/receipt.json",
            "$pilot/ranking_bootstrap0/receipt.json",
            "$pilot/ranking_bootstrap20/receipt.json",
            "$pilot/pilot_gate_verification_receipt.json",
            "$pilot/pilot_gate_verification_receipt.ranking_bootstrap20.json",
        ):
            self.assertIn(f'authenticate_sidecar "{path}"', resume)

    def test_direct_receipt_shell_sealer_accepts_valid_and_rejects_hardlink(self):
        authenticate = self.function_body("authenticate_sidecar")
        seal = self.function_body("seal_direct_receipt")
        script = (
            "set -euo pipefail\n"
            f"authenticate_sidecar() {{\n{authenticate}}}\n"
            f"seal_direct_receipt() {{\n{seal}}}\n"
            'seal_direct_receipt "$1"\n'
        )

        def invoke(path: pathlib.Path) -> subprocess.CompletedProcess[str]:
            return subprocess.run(
                ["bash", "-c", script, "runbook-sealer-test", str(path)],
                check=False,
                text=True,
                capture_output=True,
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            valid = root / "valid.json"
            valid.write_text('{"all_pass":true}\n', encoding="utf-8")
            valid_sidecar = valid.with_name(f"{valid.name}.sha256")
            valid_sidecar.write_text(
                f"{independent.sha256_path(valid)}  {valid.name}\n", encoding="ascii"
            )
            accepted = invoke(valid)
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            self.assertEqual(stat.S_IMODE(os.lstat(valid).st_mode), 0o444)
            self.assertEqual(stat.S_IMODE(os.lstat(valid_sidecar).st_mode), 0o444)
            self.assertEqual(os.lstat(valid).st_nlink, 1)
            self.assertEqual(os.lstat(valid_sidecar).st_nlink, 1)

            linked = root / "linked.json"
            linked.write_text('{"all_pass":true}\n', encoding="utf-8")
            linked_sidecar = linked.with_name(f"{linked.name}.sha256")
            linked_sidecar.write_text(
                f"{independent.sha256_path(linked)}  {linked.name}\n", encoding="ascii"
            )
            os.link(linked, root / "linked.alias.json")
            rejected = invoke(linked)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertNotEqual(stat.S_IMODE(os.lstat(linked).st_mode), 0o444)

    def test_runbook_exact11_allowlist_matches_current_source_bytes(self):
        runbook = self.RUNBOOK.read_text(encoding="utf-8").splitlines()
        marker = "sha256sum -c <<'SHA256_ALLOWLIST'"
        start = next(index for index, line in enumerate(runbook) if marker in line) + 1
        end = runbook.index("SHA256_ALLOWLIST", start)
        rows = [line.split(maxsplit=1) for line in runbook[start:end] if line.strip()]
        self.assertEqual(len(rows), 11)
        self.assertEqual(len({path for _, path in rows}), 11)
        repo_root = ROOT.parents[1]
        for expected_sha256, relative_path in rows:
            self.assertRegex(expected_sha256, r"^[0-9a-f]{64}$")
            source = repo_root / relative_path
            self.assertTrue(source.is_file(), relative_path)
            self.assertEqual(independent.sha256_path(source), expected_sha256)


class GateBoundaryTests(unittest.TestCase):
    @staticmethod
    def evaluate(incomplete: float, coverage: float, p99: int = 10, maximum: int = 20):
        time = {
            "wall_seconds": 1.0,
            "wall_hours": 1 / 3600,
            "maximum_resident_set_size_kib": 100,
            "maximum_resident_set_size_gib": 100 / 1024 / 1024,
        }
        unit = {
            "exact_enumeration_incomplete_rate": incomplete,
            "exact_limit_coverage": coverage,
            "n_incomplete_not_abstain": 0,
        }
        candidate = {
            "n_candidates": 10,
            "n_globally_certified_candidates": 10,
            "globally_certified_fraction": 1.0,
            "n_uncertified_winners": 0,
            "n_candidates_in_primary_incomplete_units": 0,
            "refinement_iterations": {"p99_nearest_rank": p99, "max": maximum},
        }
        return gate.evaluate_gates(time, time, unit, candidate)

    def test_incomplete_and_coverage_boundaries(self):
        self.assertEqual(self.evaluate(0.01, 0.90)["verdict"], "GO")
        self.assertEqual(self.evaluate(0.02, 0.85)["verdict"], "PROBE")
        self.assertEqual(self.evaluate(0.051, 0.95)["verdict"], "NO_GO")
        self.assertEqual(self.evaluate(0.0, 0.799)["verdict"], "NO_GO")

    def test_refinement_excess_is_probe_not_silent_go(self):
        result = self.evaluate(0.0, 1.0, p99=101, maximum=1001)
        self.assertEqual(result["verdict"], "PROBE")
        self.assertEqual(result["metrics"]["refinement_iterations"]["status"], "PROBE")


if __name__ == "__main__":
    unittest.main()
