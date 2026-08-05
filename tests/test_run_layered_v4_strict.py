"""Contract tests for the layered-v4 strict production launcher."""

from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts/run_layered_v4_strict.py"
SPEC = importlib.util.spec_from_file_location("run_layered_v4_strict", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
V4 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V4)


FIELDS = (
    "dataset",
    "chrom",
    "linkage_basis",
    "phase_set",
    "phase_set_status",
    "inference_role",
    "component_class",
    "tree_eligible",
    "threshold",
    "site_index",
    "pos1",
    "component_id",
    "linkage_rule",
)


def write_gzip_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row(
    component_id: str,
    site_index: int,
    *,
    role: str,
    component_class: str,
    tree_eligible: str,
) -> dict[str, str]:
    return {
        "dataset": "HCC1395",
        "chrom": "chr1",
        "linkage_basis": "PS_HP1",
        "phase_set": "PS10",
        "phase_set_status": (
            "KNOWN_PS_PRIMARY" if tree_eligible == "true" else "KNOWN_PS_SINGLETON_ABSTAIN"
        ),
        "inference_role": role,
        "component_class": component_class,
        "tree_eligible": tree_eligible,
        "threshold": "3",
        "site_index": str(site_index),
        "pos1": str(100 + site_index),
        "component_id": component_id,
        "linkage_rule": "strict_fixed_ra_endpoint_pair",
    }


class ParserAndScopeTests(unittest.TestCase):
    def parse(self, *extra: str):
        return V4.parse_args(["--manifest", "manifest.json", "--output-root", "out", *extra])

    def test_defaults_are_full_production_scope_and_validated_partition_boundary(self) -> None:
        args = self.parse()
        self.assertEqual(tuple(args.datasets), V4.CANONICAL_DATASETS)
        self.assertEqual(tuple(args.chromosomes), V4.AUTOSOMES)
        self.assertEqual(args.stage_through, "partition")
        self.assertTrue(V4.validate_scope(args))

    def test_subset_requires_explicit_partial_flag(self) -> None:
        args = self.parse("--datasets", "HCC1395", "--chromosomes", "chr1")
        with self.assertRaisesRegex(V4.RunnerError, "allow-partial-scope"):
            V4.validate_scope(args)
        args = self.parse(
            "--datasets", "HCC1395", "--chromosomes", "chr1", "--allow-partial-scope"
        )
        self.assertFalse(V4.validate_scope(args))

    def test_cache_pattern_must_bind_chromosome(self) -> None:
        args = self.parse("--extraction-cache-pattern", "/cache/{sample}")
        with self.assertRaisesRegex(V4.RunnerError, "chrom"):
            V4.validate_scope(args)

    def test_nonempty_or_existing_output_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "already_exists"
            output.mkdir()
            (output / "user.txt").write_text("preserve", encoding="utf-8")
            with self.assertRaisesRegex(V4.RunnerError, "must not already exist"):
                V4.prepare_output_root(output)
            self.assertEqual((output / "user.txt").read_text(encoding="utf-8"), "preserve")


class StrictMembershipGateTests(unittest.TestCase):
    def test_valid_singleton_is_funnel_only_and_multisite_is_primary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "membership.tsv.gz"
            rows = [
                row(
                    "singleton",
                    0,
                    role="ABSTAIN_SINGLETON_UNLINKED",
                    component_class="UNLINKED_SINGLETON_COMPONENT",
                    tree_eligible="false",
                ),
                row(
                    "linked",
                    1,
                    role="PRIMARY_PS_AWARE",
                    component_class="READ_LINKED_MULTISITE_REGION",
                    tree_eligible="true",
                ),
                row(
                    "linked",
                    2,
                    role="PRIMARY_PS_AWARE",
                    component_class="READ_LINKED_MULTISITE_REGION",
                    tree_eligible="true",
                ),
            ]
            write_gzip_tsv(path, FIELDS, rows)
            counts = V4.validate_strict_membership(
                path, sample="HCC1395", chrom="chr1"
            )
            self.assertEqual(counts["singleton_abstain_regions"], 1)
            self.assertEqual(counts["primary_multisite_regions"], 1)
            self.assertEqual(counts["primary_site_memberships"], 2)

    def test_singleton_marked_primary_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "membership.tsv.gz"
            write_gzip_tsv(
                path,
                FIELDS,
                [
                    row(
                        "bad",
                        0,
                        role="PRIMARY_PS_AWARE",
                        component_class="READ_LINKED_MULTISITE_REGION",
                        tree_eligible="true",
                    )
                ],
            )
            with self.assertRaisesRegex(V4.RunnerError, "singleton incorrectly enters"):
                V4.validate_strict_membership(path, sample="HCC1395", chrom="chr1")

    def test_legacy_membership_without_strict_columns_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy.tsv.gz"
            legacy_fields = tuple(field for field in FIELDS if field not in {"linkage_rule", "tree_eligible"})
            legacy = row(
                "legacy",
                0,
                role="PRIMARY_PS_AWARE",
                component_class="READ_LINKED_MULTISITE_REGION",
                tree_eligible="true",
            )
            write_gzip_tsv(path, legacy_fields, [{key: legacy[key] for key in legacy_fields}])
            with self.assertRaisesRegex(V4.RunnerError, "legacy/non-strict membership"):
                V4.validate_strict_membership(path, sample="HCC1395", chrom="chr1")

    def test_duplicate_site_ownership_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "membership.tsv.gz"
            first = row(
                "c1",
                0,
                role="ABSTAIN_SINGLETON_UNLINKED",
                component_class="UNLINKED_SINGLETON_COMPONENT",
                tree_eligible="false",
            )
            second = dict(first, component_id="c2")
            write_gzip_tsv(path, FIELDS, [first, second])
            with self.assertRaisesRegex(V4.RunnerError, "duplicate site ownership"):
                V4.validate_strict_membership(path, sample="HCC1395", chrom="chr1")


if __name__ == "__main__":
    unittest.main()
