#!/usr/bin/env python3
"""Synthetic contract tests for split pattern-methyl run merging."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping, Sequence

from statsmodels.stats.multitest import multipletests


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "research"
    / "20260727_multisite_pattern_methyl_topology_validation"
    / "scripts"
    / "merge_pattern_methyl_runs.py"
)
SPEC = importlib.util.spec_from_file_location("merge_pattern_methyl_runs_tested", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
merge = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = merge
SPEC.loader.exec_module(merge)
analysis = merge.ANALYSIS


def write_gzip_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    with merge.deterministic_gzip_text_writer(path, newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=analysis.RESULT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: analysis.output_value(row.get(field, ""))
                    for field in analysis.RESULT_FIELDS
                }
            )


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_details(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    with merge.deterministic_gzip_text_writer(path) as handle:
        for row in rows:
            handle.write(
                json.dumps(
                    row,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                + "\n"
            )


def evidence_row(
    region_id: str,
    *,
    confirmatory: bool,
    p_value: float,
    permutations: int,
) -> dict[str, Any]:
    row = {field: "" for field in analysis.RESULT_FIELDS}
    row.update(
        {
            "schema_version": analysis.SCHEMA_VERSION,
            "dataset": "SYNTH",
            "chrom": "chr1",
            "region_id": region_id,
            "unit_id": f"unit-{region_id}",
            "phase_set": "100",
            "hp_family": "1",
            "hp_raw": "1-1",
            "active_positions": "100,200",
            "n_active_bits": "2",
            "pair_full4": "true" if confirmatory else "false",
            "k_ge_3": "false",
            "input_n_complete": "40",
            "input_state_counts_json": '{"AA":20,"RR":20}',
            "analysis_n": "40",
            "analysis_state_counts_json": '{"AA":20,"RR":20}',
            "n_common_cpg": "20",
            "n_distal_cpg": "12",
            "qname_join_fraction_min": "1",
            "tile_overlap_conflicts": "0",
            "exchangeable_strata": "2",
            "exchangeable_n": "40",
            "permanova_pseudo_f": "12",
            "permanova_r2": "0.2",
            "permanova_p": str(p_value),
            "permanova_permutations_requested": str(permutations),
            "permanova_permutations_realized": str(permutations),
            "permdisp_f": "0.5",
            "permdisp_p": "0.5",
            "best_pair": "AA|RR",
            "best_pair_hamming": "2",
            "best_pair_between_mean": "0.6",
            "best_pair_pooled_within_mean": "0.2",
            "best_pair_distance_contrast": "0.4",
            "best_pair_standardized_effect": "1.0",
            "best_pair_topology_relation": "PAIR_BAND_ONLY_HAMMING_GT1",
            "max_geometry_smd": "0.1",
            "geometry_feature": "read_length:AA|RR",
            "all_states_n8": "true",
            "all_states_n10": "true",
            "equal_n_r2": "0.18",
            "equal_n_retention": "0.9",
            "rarefaction_median_r2": "0.17",
            "rarefaction_retention": "0.85",
            "distal_r2": "0.16",
            "distal_retention": "0.8",
            "multiplicity_family": (
                merge.CONFIRMATORY_FAMILY
                if confirmatory
                else merge.SECONDARY_FAMILY
            ),
            "q_by": "",
            "p_holm": "",
            "assessment": "PENDING_MULTIPLICITY",
            "evaluation_status": "EVALUABLE",
            "invalid_reason": "",
        }
    )
    return row


class MergePatternMethylRunsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="pattern-methyl-merge-")
        self.root = Path(self.temporary.name)
        self.confirm_rows = [
            evidence_row(
                "chr1:100-200:ps100:hpf1",
                confirmatory=True,
                p_value=0.001,
                permutations=999,
            ),
            evidence_row(
                "chr1:300-400:ps100:hpf1",
                confirmatory=True,
                p_value=0.01,
                permutations=999,
            ),
        ]
        self.secondary_rows = [
            evidence_row(
                "chr1:500-600:ps100:hpf1",
                confirmatory=False,
                p_value=0.0001,
                permutations=199,
            )
        ]
        self.confirm = self._write_bundle(
            "confirmatory", self.confirm_rows, permutations=999
        )
        self.secondary = self._write_bundle(
            "secondary", self.secondary_rows, permutations=199
        )
        loaded, _ = merge.load_evidence(
            self.confirm["evidence"],
            label="fixture_confirmatory",
            expected_family=merge.CONFIRMATORY_FAMILY,
        )
        self.refined_row = dict(loaded[0])
        self.refined_row["permanova_p"] = "0.00002"
        self.refined_row["permdisp_p"] = "0.55"
        self.refined_row["permanova_permutations_requested"] = "49999"
        self.refined_row["permanova_permutations_realized"] = "49999"
        self.refined = self._write_bundle(
            "refined", [self.refined_row], permutations=49999
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write_bundle(
        self,
        name: str,
        source_rows: Sequence[Mapping[str, Any]],
        *,
        permutations: int,
    ) -> dict[str, Path]:
        directory = self.root / name
        directory.mkdir()
        rows = [dict(row) for row in source_rows]
        analysis.adjust_p_values(rows)
        analysis.classify_rows(rows)
        details = [
            {
                "schema_name": f"{analysis.SCHEMA_NAME}.detail",
                "schema_version": analysis.SCHEMA_VERSION,
                "dataset": row["dataset"],
                "chrom": row["chrom"],
                "region_id": row["region_id"],
                "hp_raw": row["hp_raw"],
                "assessment": row["assessment"],
                "synthetic_payload": row["unit_id"],
            }
            for row in rows
        ]
        evidence_path = directory / "evidence.tsv.gz"
        details_path = directory / "details.jsonl.gz"
        summary_path = directory / "summary.json"
        write_gzip_tsv(evidence_path, rows)
        write_details(details_path, details)
        refinement_sources: dict[str, Any] = {
            "unit_key_file": None,
            "unit_key_receipt": None,
            "screen_evidence": None,
        }
        if name.startswith("refined"):
            unit_key_path = directory / "selected.tsv"
            with unit_key_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=analysis.KEY_FIELDS
                    if hasattr(analysis, "KEY_FIELDS")
                    else ("dataset", "chrom", "region_id", "hp_raw"),
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="ignore",
                )
                writer.writeheader()
                writer.writerows(rows)
            screen_identity = merge.file_identity(self.confirm["evidence"])
            receipt_path = directory / "selected.receipt.json"
            receipt = {
                "schema_name": analysis.ADAPTIVE_SELECTION_SCHEMA_NAME,
                "schema_version": analysis.ADAPTIVE_SELECTION_SCHEMA_VERSION,
                "all_pass": True,
                "contract": {
                    "family": merge.CONFIRMATORY_FAMILY,
                    "screen_permutations": 999,
                    "screen_floor": 0.001,
                    "refinement_permutations": 49999,
                },
                "counts": {"selected_for_refinement": len(rows)},
                "inputs": {"screen_evidence": screen_identity},
                "outputs": {
                    "unit_keys": merge.file_identity(unit_key_path),
                },
            }
            receipt_path.write_text(
                json.dumps(receipt, sort_keys=True, indent=2) + "\n",
                encoding="utf-8",
            )
            refinement_sources = {
                "unit_key_file": merge.file_identity(unit_key_path),
                "unit_key_receipt": merge.file_identity(receipt_path),
                "screen_evidence": screen_identity,
            }
        summary = {
            "schema_name": f"{analysis.SCHEMA_NAME}.summary",
            "schema_version": analysis.SCHEMA_VERSION,
            "result_status": (
                merge.PROVISIONAL_RESULT_STATUS
                if name.startswith("refined")
                else merge.AUTHORITATIVE_RESULT_STATUS
            ),
            "config": {"permutations": permutations, "seed": 20260727},
            "counts": merge.expected_counts(rows, len(details)),
            "sources": {
                "pattern_counts": {"path": "/frozen/counts", "sha256": "a" * 64},
                "artifact_catalog": {"path": "/frozen/catalog", "sha256": "b" * 64},
                "candidate_shards": [
                    {"path": "/frozen/shard", "sha256": "c" * 64, "size_bytes": 1}
                ],
                "candidate_manifest": {
                    "path": "/frozen/manifest",
                    "sha256": "d" * 64,
                },
                **refinement_sources,
                "topology_root": "/frozen/topology",
                "topology_jsonl": [
                    {
                        "dataset": "SYNTH",
                        "path": "/frozen/topology/SYNTH.jsonl",
                        "size_bytes": 1,
                        "sha256": "f" * 64,
                    }
                ],
            },
            "outputs": {
                "evidence": merge.file_identity(evidence_path),
                "details": merge.file_identity(details_path),
            },
        }
        summary_path.write_text(
            json.dumps(summary, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        return {
            "evidence": evidence_path,
            "details": details_path,
            "summary": summary_path,
        }

    def _args(
        self,
        output: Path,
        *,
        confirm: Mapping[str, Path] | None = None,
        secondary: Mapping[str, Path] | None = None,
        refined: Mapping[str, Path] | None | object = ...,
    ) -> SimpleNamespace:
        confirm = confirm or self.confirm
        secondary = secondary or self.secondary
        refined_bundle = self.refined if refined is ... else refined
        if refined_bundle is None:
            refined_evidence = None
            refined_details = None
            refined_summary = None
        else:
            refined_evidence = refined_bundle["evidence"]  # type: ignore[index]
            refined_details = refined_bundle["details"]  # type: ignore[index]
            refined_summary = refined_bundle["summary"]  # type: ignore[index]
        return SimpleNamespace(
            confirmatory_evidence=confirm["evidence"],
            confirmatory_details=confirm["details"],
            confirmatory_summary=confirm["summary"],
            secondary_evidence=secondary["evidence"],
            secondary_details=secondary["details"],
            secondary_summary=secondary["summary"],
            refined_confirmatory_evidence=refined_evidence,
            refined_confirmatory_details=refined_details,
            refined_confirmatory_summary=refined_summary,
            output_dir=output,
        )

    def _rewrite_refined_receipt(
        self,
        mutate: Any,
        *,
        bundle: Mapping[str, Path] | None = None,
    ) -> None:
        target = bundle or self.refined
        summary = json.loads(target["summary"].read_text(encoding="utf-8"))
        receipt_path = Path(summary["sources"]["unit_key_receipt"]["path"])
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        mutate(receipt)
        receipt_path.write_text(
            json.dumps(receipt, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        summary["sources"]["unit_key_receipt"] = merge.file_identity(receipt_path)
        target["summary"].write_text(
            json.dumps(summary, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )

    def test_merge_refines_and_recomputes_each_family_deterministically(self) -> None:
        first = self.root / "combined-first"
        second = self.root / "combined-second"
        result = merge.execute(self._args(first))
        merge.execute(self._args(second))

        rows = read_gzip_tsv(first / merge.OUTPUT_EVIDENCE)
        by_region = {row["region_id"]: row for row in rows}
        refined = by_region["chr1:100-200:ps100:hpf1"]
        other_confirm = by_region["chr1:300-400:ps100:hpf1"]
        secondary = by_region["chr1:500-600:ps100:hpf1"]
        expected_by = multipletests(
            [0.00002, 0.01], alpha=0.05, method="fdr_by"
        )[1]
        expected_holm = multipletests(
            [0.00002, 0.01], alpha=0.05, method="holm"
        )[1]

        self.assertEqual(refined["permanova_permutations_realized"], "49999")
        self.assertAlmostEqual(float(refined["q_by"]), float(expected_by[0]))
        self.assertAlmostEqual(float(other_confirm["q_by"]), float(expected_by[1]))
        self.assertAlmostEqual(float(refined["p_holm"]), float(expected_holm[0]))
        self.assertAlmostEqual(float(secondary["q_by"]), 0.0001)
        self.assertEqual(refined["assessment"], "ROBUST_ASSOCIATION")
        self.assertEqual(
            secondary["assessment"], "EVALUABLE_NO_ROBUST_ASSOCIATION"
        )

        receipt = json.loads((first / merge.OUTPUT_RECEIPT).read_text())
        self.assertTrue(receipt["pass"])
        self.assertEqual(receipt["counts"]["units_total"], 3)
        self.assertEqual(receipt["counts"]["refined_confirmatory_units"], 1)
        for artifact in ("evidence", "details", "summary"):
            self.assertEqual(
                receipt["inputs"]["refined_confirmatory"][artifact]["sha256"],
                merge.sha256_file(self.refined[artifact]),
            )
        self.assertEqual(
            receipt["outputs"]["evidence"]["sha256"],
            merge.sha256_file(first / merge.OUTPUT_EVIDENCE),
        )
        expected_sidecar = (
            f"{merge.sha256_file(first / merge.OUTPUT_RECEIPT)}  "
            f"{merge.OUTPUT_RECEIPT}\n"
        )
        self.assertEqual(
            (first / merge.OUTPUT_RECEIPT_SHA256).read_text(),
            expected_sidecar,
        )
        self.assertEqual(result["counts"], receipt["counts"])

        for filename in (
            merge.OUTPUT_EVIDENCE,
            merge.OUTPUT_DETAILS,
            merge.OUTPUT_SUMMARY,
            merge.OUTPUT_RECEIPT,
            merge.OUTPUT_RECEIPT_SHA256,
        ):
            self.assertEqual(
                (first / filename).read_bytes(),
                (second / filename).read_bytes(),
                filename,
            )

    def test_refinement_fails_on_identity_drift(self) -> None:
        drift = dict(self.refined_row)
        drift["analysis_n"] = "41"
        refined = self._write_bundle(
            "refined-drift", [drift], permutations=49999
        )
        output = self.root / "identity-drift-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "identity/statistic drift"
        ):
            merge.execute(self._args(output, refined=refined))
        self.assertFalse(output.exists())

    def test_merge_without_optional_refinement_retains_base_budget(self) -> None:
        noneligible_rows = [dict(row) for row in self.confirm_rows]
        noneligible_rows[0]["permanova_p"] = "0.002"
        noneligible_confirm = self._write_bundle(
            "confirmatory-no-refinement",
            noneligible_rows,
            permutations=999,
        )
        output = self.root / "combined-without-refinement"
        result = merge.execute(
            self._args(output, confirm=noneligible_confirm, refined=None)
        )
        rows = read_gzip_tsv(output / merge.OUTPUT_EVIDENCE)
        by_region = {row["region_id"]: row for row in rows}
        self.assertEqual(
            by_region["chr1:100-200:ps100:hpf1"][
                "permanova_permutations_realized"
            ],
            "999",
        )
        self.assertEqual(result["counts"]["refined_confirmatory_units"], 0)
        receipt = json.loads((output / merge.OUTPUT_RECEIPT).read_text())
        self.assertIsNone(receipt["inputs"]["refined_confirmatory"])

    def test_zero_primary_r2_accepts_undefined_retention_as_non_robust(self) -> None:
        zero = dict(self.secondary_rows[0])
        zero["permanova_r2"] = "0"
        zero["permanova_p"] = "1"
        zero["equal_n_r2"] = "0"
        zero["equal_n_retention"] = ""
        zero["rarefaction_median_r2"] = "0"
        zero["rarefaction_retention"] = ""
        secondary = self._write_bundle(
            "secondary-zero-r2", [zero], permutations=199
        )
        output = self.root / "combined-zero-r2"
        merge.execute(
            self._args(output, secondary=secondary)
        )
        rows = read_gzip_tsv(output / merge.OUTPUT_EVIDENCE)
        row = next(item for item in rows if item["region_id"] == zero["region_id"])
        self.assertEqual(row["permanova_r2"], "0")
        self.assertEqual(row["equal_n_retention"], "")
        self.assertEqual(
            row["assessment"], "EVALUABLE_NO_ROBUST_ASSOCIATION"
        )

    def test_refinement_fails_when_realized_permutations_do_not_increase(self) -> None:
        stale = dict(self.refined_row)
        stale["permanova_permutations_requested"] = "999"
        stale["permanova_permutations_realized"] = "999"
        refined = self._write_bundle(
            "refined-stale", [stale], permutations=999
        )
        output = self.root / "stale-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "frozen permutation budget"
        ):
            merge.execute(self._args(output, refined=refined))
        self.assertFalse(output.exists())

    def test_refinement_requires_complete_three_artifact_bundle(self) -> None:
        output = self.root / "partial-refinement-output"
        args = self._args(output)
        args.refined_confirmatory_details = None
        with self.assertRaisesRegex(
            merge.MergeContractError, "must be provided together"
        ):
            merge.execute(args)
        self.assertFalse(output.exists())
        self.assertFalse((self.root / "_failed_staging_archive").exists())

    def test_refinement_summary_binds_details_sha256(self) -> None:
        with gzip.open(
            self.refined["details"], "rt", encoding="utf-8"
        ) as handle:
            details = [json.loads(line) for line in handle]
        details[0]["synthetic_payload"] = "tampered"
        write_details(self.refined["details"], details)

        output = self.root / "refined-details-hash-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "summary output sha256 mismatch"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_refinement_fails_on_detail_identity_drift(self) -> None:
        with gzip.open(
            self.refined["details"], "rt", encoding="utf-8"
        ) as handle:
            details = [json.loads(line) for line in handle]
        details[0]["synthetic_payload"] = "identity-drift"
        write_details(self.refined["details"], details)
        summary = json.loads(self.refined["summary"].read_text())
        summary["outputs"]["details"] = merge.file_identity(
            self.refined["details"]
        )
        self.refined["summary"].write_text(
            json.dumps(summary, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )

        output = self.root / "refined-detail-identity-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "detail identity/statistic drift"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_refinement_config_must_match_base_except_permutations(self) -> None:
        summary = json.loads(self.refined["summary"].read_text())
        summary["config"]["seed"] = 7
        self.refined["summary"].write_text(
            json.dumps(summary, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )

        output = self.root / "refined-config-drift-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "config drift outside permutations"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_runtime_failure_archives_staging_without_deletion(self) -> None:
        output = self.root / "runtime-failure-output"
        with mock.patch.object(
            merge,
            "write_details",
            side_effect=RuntimeError("synthetic write failure"),
        ):
            with self.assertRaisesRegex(
                merge.MergeContractError, "staging archived"
            ):
                merge.execute(self._args(output))

        self.assertFalse(output.exists())
        archive_root = self.root / "_failed_staging_archive"
        archived = list(archive_root.iterdir())
        self.assertEqual(len(archived), 1)
        self.assertTrue((archived[0] / merge.OUTPUT_EVIDENCE).is_file())
        self.assertFalse((archived[0] / merge.OUTPUT_DETAILS).exists())

    def test_base_family_exact_key_overlap_fails_closed(self) -> None:
        overlapping_secondary_row = evidence_row(
            self.confirm_rows[0]["region_id"],
            confirmatory=False,
            p_value=0.01,
            permutations=199,
        )
        overlapping_secondary = self._write_bundle(
            "secondary-overlap",
            [overlapping_secondary_row],
            permutations=199,
        )
        output = self.root / "overlap-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "exact keys overlap"
        ):
            merge.execute(
                self._args(
                    output,
                    secondary=overlapping_secondary,
                )
            )
        self.assertFalse(output.exists())

    def test_summary_count_drift_is_rejected_before_publish(self) -> None:
        invalid_summary = self.root / "invalid-summary.json"
        payload = json.loads(self.secondary["summary"].read_text())
        payload["counts"]["units_total"] = 999
        invalid_summary.write_text(
            json.dumps(payload, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        invalid_secondary = dict(self.secondary)
        invalid_secondary["summary"] = invalid_summary
        output = self.root / "summary-drift-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "summary counts mismatch"
        ):
            merge.execute(
                self._args(
                    output,
                    secondary=invalid_secondary,
                    refined=None,
                )
            )
        self.assertFalse(output.exists())

    def test_base_source_drift_is_rejected(self) -> None:
        payload = json.loads(self.secondary["summary"].read_text())
        payload["sources"]["artifact_catalog"]["sha256"] = "0" * 64
        self.secondary["summary"].write_text(
            json.dumps(payload, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        output = self.root / "source-drift-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "source contract drift"
        ):
            merge.execute(self._args(output, refined=None))
        self.assertFalse(output.exists())

    def test_base_config_drift_is_rejected(self) -> None:
        payload = json.loads(self.secondary["summary"].read_text())
        payload["config"]["seed"] = 7
        self.secondary["summary"].write_text(
            json.dumps(payload, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        output = self.root / "base-config-drift-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "config drift outside permutations"
        ):
            merge.execute(self._args(output, refined=None))
        self.assertFalse(output.exists())

    def test_refinement_requires_provisional_status(self) -> None:
        payload = json.loads(self.refined["summary"].read_text())
        payload["result_status"] = merge.AUTHORITATIVE_RESULT_STATUS
        self.refined["summary"].write_text(
            json.dumps(payload, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        output = self.root / "refined-status-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "not marked provisional"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_frozen_family_budgets_are_enforced(self) -> None:
        wrong_rows = [dict(row) for row in self.confirm_rows]
        for row in wrong_rows:
            row["permanova_permutations_requested"] = "998"
            row["permanova_permutations_realized"] = "998"
        wrong = self._write_bundle(
            "confirmatory-wrong-budget",
            wrong_rows,
            permutations=998,
        )
        output = self.root / "wrong-budget-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "frozen permutation budget 999"
        ):
            merge.execute(self._args(output, confirm=wrong))
        self.assertFalse(output.exists())

    def test_base_sources_must_not_claim_refinement_bindings(self) -> None:
        summary = json.loads(self.confirm["summary"].read_text(encoding="utf-8"))
        summary["sources"]["unit_key_file"] = {
            "path": "/forged/selection.tsv",
            "size_bytes": 1,
            "sha256": "0" * 64,
        }
        self.confirm["summary"].write_text(
            json.dumps(summary, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        output = self.root / "base-refinement-source-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "base sources must set refinement bindings"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_refined_receipt_screen_binding_must_equal_confirmatory_evidence(
        self,
    ) -> None:
        self._rewrite_refined_receipt(
            lambda receipt: receipt["inputs"]["screen_evidence"].update(
                {"sha256": "0" * 64}
            )
        )
        output = self.root / "screen-binding-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "screen_evidence sha256 mismatch"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_refined_receipt_selected_count_must_equal_refined_keys(self) -> None:
        self._rewrite_refined_receipt(
            lambda receipt: receipt["counts"].update(
                {"selected_for_refinement": 2}
            )
        )
        output = self.root / "selected-count-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "selected count does not equal refined keys"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_ineligible_refined_key_is_rejected_by_recomputed_screen_gate(
        self,
    ) -> None:
        rows = read_gzip_tsv(self.confirm["evidence"])
        rows[0]["equal_n_retention"] = "0.1"
        write_gzip_tsv(self.confirm["evidence"], rows)
        summary = json.loads(self.confirm["summary"].read_text(encoding="utf-8"))
        summary["outputs"]["evidence"] = merge.file_identity(
            self.confirm["evidence"]
        )
        self.confirm["summary"].write_text(
            json.dumps(summary, sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        output = self.root / "ineligible-refinement-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "do not equal frozen adaptive-gate selection"
        ):
            merge.execute(self._args(output))
        self.assertFalse(output.exists())

    def test_eligible_screen_cannot_silently_skip_refinement(self) -> None:
        output = self.root / "missing-required-refinement-output"
        with self.assertRaisesRegex(
            merge.MergeContractError, "no refined confirmatory bundle"
        ):
            merge.execute(self._args(output, refined=None))
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
