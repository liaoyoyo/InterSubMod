#!/usr/bin/env python3
import csv
import gzip
import pathlib
import sys
import tempfile
import unittest
from collections import Counter

import pysam

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from extract_lossless_read_linkage import (  # noqa: E402
    CoordinateSidecarJoiner,
    Variant,
    apply_bridge,
    call_alignment,
    cigar_digest,
    extract_one,
    finish_difference,
    full_alignment_key,
    hp_family,
    prepare_output_directory,
    sparse_any_phase_support_at_active_boundaries,
    sparse_bridge_events,
    sparse_support_at_active_boundaries,
)


def fixture_alignment(header, name, start, cigar, sequence, *, mapq=60, flag=0, qualities=None):
    row = pysam.AlignedSegment(header)
    row.query_name = name
    row.query_sequence = sequence
    row.flag = flag
    row.reference_id = 0
    row.reference_start = start
    row.mapping_quality = mapq
    row.cigarstring = cigar
    row.query_qualities = pysam.qualitystring_to_array(qualities or ("I" * len(sequence)))
    row.set_tag("RG", "rg1")
    return row


def build_extraction_fixture(root: pathlib.Path, *, duplicate_primary=False, cross_phase=False):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}], "RG": [{"ID": "rg1"}]}
    )
    rows = [
        fixture_alignment(header, "hp1", 9, "3M", "ACT", qualities="I!I"),
        fixture_alignment(header, "reason", 9, "1M1D1N2M", "CCC"),
        fixture_alignment(header, "lowmap", 9, "5M", "AAAAA", mapq=5),
        fixture_alignment(header, "duplicate", 9, "5M", "AAAAA", flag=0x400),
        fixture_alignment(header, "qcfail", 9, "5M", "AAAAA", flag=0x200),
        fixture_alignment(header, "secondary", 9, "5M", "AAAAA", flag=0x100),
        fixture_alignment(header, "supplementary", 9, "5M", "AAAAA", flag=0x800),
        fixture_alignment(header, "hp2", 11, "3M", "GAG"),
    ]
    if duplicate_primary:
        rows.append(fixture_alignment(header, "hp1", 11, "3M", "GAG"))
    if cross_phase:
        rows.append(fixture_alignment(header, "hp1b", 11, "3M", "GAG"))
    rows.sort(key=lambda row: (row.reference_start, row.query_name, row.flag))

    bam = root / "tiny.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
        for row in rows:
            handle.write(row)
    pysam.index(str(bam))

    plain_vcf = root / "tree.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##FILTER=<ID=PASS,Description=All filters passed>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tT\t.\tPASS\t.\n"
        "chr1\t11\t.\tG\tC\t.\tPASS\t.\n"
        "chr1\t12\t.\tG\tT\t.\tPASS\t.\n"
        "chr1\t13\t.\tA\tG\t.\tPASS\t.\n"
        "chr1\t14\t.\tT\tG\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    vcf = root / "tree.vcf.gz"
    pysam.tabix_compress(str(plain_vcf), str(vcf), force=True)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)

    plain_sidecar = root / "tags.tsv"
    with plain_sidecar.open("w", encoding="utf-8") as handle:
        handle.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
        for row in rows:
            hp = "1-1" if row.query_name.startswith("hp1") else "2-1" if row.query_name == "hp2" else "."
            ps = "200" if row.query_name == "hp1b" else "100" if hp != "." else "."
            handle.write(
                f"chr1\t{row.reference_start}\t{row.reference_end}\t{row.query_name}\t{row.flag}\t"
                f"{row.mapping_quality}\t{cigar_digest(row.cigarstring)}\t{hp}\t{ps}\n"
            )
    sidecar = root / "tags.tsv.gz"
    pysam.tabix_compress(str(plain_sidecar), str(sidecar), force=True)
    pysam.tabix_index(str(sidecar), preset="bed", force=True)
    sidecar_index = pathlib.Path(f"{sidecar}.tbi")

    entry = {
        "sample": "COLO829",
        "alignment_payload": {"path": str(bam)},
        "somatic": {"tree_vcf": {"path": str(vcf)}},
        "read_tags": {"sidecar": {"path": str(sidecar)}, "index": {"path": str(sidecar_index)}},
    }
    return entry


def alignment(seq="ACGTACGT", cigar="8M", start=9, flag=0, qualities=None):
    row = pysam.AlignedSegment()
    row.query_name = "read1"
    row.query_sequence = seq
    row.flag = flag
    row.reference_id = 0
    row.reference_start = start
    row.mapping_quality = 60
    row.cigarstring = cigar
    row.query_qualities = pysam.qualitystring_to_array(qualities or ("I" * len(seq)))
    # Unit-test segments do not carry a header, so reference_name is patched by
    # assigning a minimal header-backed clone in the sidecar test only.
    return row


class CigarCallingTests(unittest.TestCase):
    def test_match_other_and_low_quality(self):
        aln = alignment(seq="ACGT", cigar="4M", start=9, qualities="I!II")
        variants = [Variant("chr1", 10, "A", "T"), Variant("chr1", 11, "G", "C"), Variant("chr1", 12, "G", "A")]
        indices, calls, qualities = call_alignment(aln, variants, [9, 10, 11], 20)
        self.assertEqual(indices, (0, 1, 2))
        self.assertEqual(calls, ("R", "L", "R"))
        self.assertEqual(qualities[1], 0)

    def test_deletion_and_refskip_are_distinct(self):
        aln = alignment(seq="ACGT", cigar="2M1D1N2M", start=9)
        variants = [
            Variant("chr1", 10, "A", "T"),
            Variant("chr1", 12, "A", "T"),
            Variant("chr1", 13, "A", "T"),
            Variant("chr1", 14, "G", "A"),
        ]
        _, calls, _ = call_alignment(aln, variants, [9, 11, 12, 13], 20)
        self.assertEqual(calls, ("R", "D", "S", "R"))

    def test_insertion_does_not_shift_reference_calls(self):
        aln = alignment(seq="ACTGT", cigar="2M1I2M", start=9)
        variants = [Variant("chr1", 10, "A", "T"), Variant("chr1", 12, "G", "A")]
        _, calls, _ = call_alignment(aln, variants, [9, 11], 20)
        self.assertEqual(calls, ("R", "R"))

    def test_other_and_uncovered_codes_are_explicit(self):
        other = alignment(seq="C", cigar="1M", start=9)
        _, calls, qualities = call_alignment(other, [Variant("chr1", 10, "A", "T")], [9], 20)
        self.assertEqual(calls, ("O",))
        self.assertEqual(qualities, (40,))

        uncovered = alignment(seq="A", cigar="1M", start=9)
        uncovered.query_sequence = None
        _, calls, qualities = call_alignment(uncovered, [Variant("chr1", 10, "A", "T")], [9], 20)
        self.assertEqual(calls, ("X",))
        self.assertEqual(qualities, (None,))


class JoinAndBridgeTests(unittest.TestCase):
    def test_coordinate_sidecar_exact_join(self):
        header = pysam.AlignmentHeader.from_references(["chr1"], [1000])
        aln = pysam.AlignedSegment(header)
        aln.query_name = "read1"
        aln.query_sequence = "ACGT"
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = 9
        aln.mapping_quality = 60
        aln.cigarstring = "4M"
        aln.query_qualities = pysam.qualitystring_to_array("IIII")
        digest = cigar_digest("4M")
        line = f"chr1\t9\t13\tread1\t0\t60\t{digest}\t1-1\t100"
        joiner = CoordinateSidecarJoiner(iter([line]), "chr1")
        tags = joiner.lookup(aln)
        self.assertEqual(tags.hp, "1-1")
        self.assertEqual(tags.ps, "100")
        self.assertEqual(joiner.matched, 1)
        self.assertEqual(full_alignment_key(aln)[-1], digest)

    def test_coordinate_sidecar_conflicting_duplicate_fails_closed(self):
        header = pysam.AlignmentHeader.from_references(["chr1"], [1000])
        aln = pysam.AlignedSegment(header)
        aln.query_name = "read1"
        aln.query_sequence = "ACGT"
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = 9
        aln.mapping_quality = 60
        aln.cigarstring = "4M"
        aln.query_qualities = pysam.qualitystring_to_array("IIII")
        digest = cigar_digest("4M")
        lines = iter(
            [
                f"chr1\t9\t13\tread1\t0\t60\t{digest}\t1\t100",
                f"chr1\t9\t13\tread1\t0\t60\t{digest}\t2\t100",
            ]
        )
        with self.assertRaisesRegex(RuntimeError, "conflicting sidecar tags"):
            CoordinateSidecarJoiner(lines, "chr1").lookup(aln)

    def test_bridge_difference_counts_one_molecule_once_per_cut(self):
        difference = [0, 0, 0, 0]
        apply_bridge(difference, [0, 1, 2, 3])
        self.assertEqual(finish_difference(difference, 4), (1, 1, 1))

    def test_sparse_ps_support_matches_dense_reference(self):
        n_sites = 12
        bridges = ((0, 4), (0, 7), (2, 7), (7, 11))
        dense = [0] * n_sites
        events = Counter()
        for bridge in bridges:
            apply_bridge(dense, bridge)
            sparse_bridge_events(events, bridge)
        active = (0, 2, 4, 7, 11)
        dense_support = finish_difference(dense, n_sites)
        expected = tuple(dense_support[index] for index in active[:-1])
        self.assertEqual(
            sparse_support_at_active_boundaries(active, events, n_sites),
            expected,
        )

    def test_sparse_ps_support_does_not_scale_with_chromosome_catalog_length(self):
        # A billion-site catalog would be infeasible for the former range(n_sites)
        # scan, but sparse event/query work remains constant for this fixture.
        n_sites = 1_000_000_000
        active = (3, n_sites - 1)
        events = Counter({3: 2, n_sites - 1: -2})
        self.assertEqual(
            sparse_support_at_active_boundaries(active, events, n_sites),
            (2,),
        )

    def test_sparse_any_known_ps_threshold_union_matches_dense_reference(self):
        n_sites = 12
        active = (0, 2, 4, 7, 11)
        phase_bridges = (
            ((0, 4), (0, 7)),
            ((2, 7), (7, 11)),
        )
        phase_events = []
        dense_supports = []
        for bridges in phase_bridges:
            events = Counter()
            dense = [0] * n_sites
            for bridge in bridges:
                sparse_bridge_events(events, bridge)
                apply_bridge(dense, bridge)
            phase_events.append(events)
            dense_supports.append(finish_difference(dense, n_sites))
        observed = sparse_any_phase_support_at_active_boundaries(
            active, phase_events, (1, 2, 3), n_sites
        )
        for threshold in (1, 2, 3):
            expected = tuple(
                any(values[cut] >= threshold for values in dense_supports)
                for cut in active[:-1]
            )
            self.assertEqual(observed[threshold], expected)

    def test_hp_family_contract(self):
        self.assertEqual(hp_family("1-2"), "1")
        self.assertEqual(hp_family("2-1"), "2")
        self.assertEqual(hp_family("3"), "3")
        self.assertEqual(hp_family("."), "none")
        for invalid in ("0", "1foo", "2-3", "5"):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "unexpected LongPhase-S HP"):
                hp_family(invalid)


class SyntheticExtractionTests(unittest.TestCase):
    def test_preflight_created_output_directory_contract_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            missing = root / "missing"
            with self.assertRaisesRegex(FileExistsError, "unavailable"):
                prepare_output_directory(missing, require_existing_empty=True)

            empty = root / "empty"
            empty.mkdir()
            prepare_output_directory(empty, require_existing_empty=True)
            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                prepare_output_directory(empty)

            (empty / "prior-output").write_text("x", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "not empty"):
                prepare_output_directory(empty, require_existing_empty=True)

            target = root / "target"
            target.mkdir()
            symlink = root / "symlink"
            symlink.symlink_to(target, target_is_directory=True)
            with self.assertRaisesRegex(FileExistsError, "not a real directory"):
                prepare_output_directory(symlink, require_existing_empty=True)

    def test_extraction_preserves_preflight_output_directory_inode(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            entry = build_extraction_fixture(root)
            output = root / "output"
            output.mkdir()
            inode_before = output.stat().st_ino
            receipt = extract_one(
                entry,
                "chr1",
                output,
                mapq_min=20,
                baseq_min=20,
                thresholds=(1, 2),
                samtools_threads=0,
                require_existing_empty_output_dir=True,
            )
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(output.stat().st_ino, inode_before)
            self.assertTrue((output / "receipt.json").is_file())

    def test_funnel_reason_mass_receipt_and_per_hp_components(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            entry = build_extraction_fixture(root)
            output = root / "output"
            receipt = extract_one(
                entry,
                "chr1",
                output,
                mapq_min=20,
                baseq_min=20,
                thresholds=(1, 2),
                samtools_threads=0,
                provenance={
                    "manifest": {"path": "/synthetic/manifest.json", "sha256": "m" * 64},
                    "extractor": {"path": "/synthetic/extractor.py", "sha256": "e" * 64},
                },
            )
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(receipt["schema_version"], "1.2.0")
            self.assertTrue(all(receipt["checks"].values()))
            counts = receipt["counts"]
            self.assertEqual(counts["raw_overlapping_alignments"], 8)
            self.assertEqual(counts["excluded_by_flag"], 4)
            self.assertEqual(counts["mapq_rejected_after_flag"], 1)
            self.assertEqual(counts["canonical_eligible_alignments"], 3)
            self.assertEqual(counts["molecule_sparse_rows_written"], 3)
            self.assertEqual(counts["sidecar_exact_matches"], 3)

            sites_path = output / "COLO829.chr1.site_catalog.tsv.gz"
            with gzip.open(sites_path, "rt", encoding="utf-8", newline="") as handle:
                site_rows = list(csv.DictReader(handle, delimiter="\t"))
            reason_totals = {
                code: sum(int(row[code]) for row in site_rows) for code in ("R", "A", "O", "D", "S", "L", "X")
            }
            for observed in ("R", "A", "O", "D", "S", "L"):
                self.assertGreater(reason_totals[observed], 0, observed)
            self.assertEqual(sum(reason_totals.values()), counts["site_call_rows_sparse"])

            summary = receipt["component_summary_by_linkage_basis"]
            self.assertEqual(summary["pooled"]["1"]["max_k"], 5)
            self.assertEqual(summary["HP1"]["1"]["max_k"], 3)
            self.assertEqual(summary["HP2"]["1"]["max_k"], 3)
            self.assertNotEqual(
                receipt["component_digests_by_linkage_basis"]["pooled"]["1"],
                receipt["component_digests_by_linkage_basis"]["HP1"]["1"],
            )
            components_path = output / "COLO829.chr1.components.tsv.gz"
            with gzip.open(components_path, "rt", encoding="utf-8", newline="") as handle:
                component_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {row["linkage_basis"] for row in component_rows},
                {"pooled", "HP1", "HP2", "HP3", "HP4", "unphased", "PS_HP1", "PS_HP2"},
            )
            primary = [row for row in component_rows if row["inference_role"] == "PRIMARY_PS_AWARE"]
            self.assertTrue(primary)
            self.assertTrue(all(row["phase_set"] == "100" for row in primary))
            ps_audit = receipt["legacy_cross_phase_set_aggregation_audit"]["1"]["1"]
            self.assertIn("legacy_active_site_components", ps_audit)
            self.assertIn("known_ps_components_sum", ps_audit)
            self.assertTrue((output / "receipt.json.sha256").is_file())

    def test_same_dataset_read_group_qname_with_two_primary_alignments_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            entry = build_extraction_fixture(root, duplicate_primary=True)
            with self.assertRaisesRegex(RuntimeError, "multiple canonical primary alignments for molecule"):
                extract_one(
                    entry,
                    "chr1",
                    root / "output",
                    mapq_min=20,
                    baseq_min=20,
                    thresholds=(1,),
                    samtools_threads=0,
                )

    def test_same_site_can_belong_to_two_phase_sets_without_primary_support_mixing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            receipt = extract_one(
                build_extraction_fixture(root, cross_phase=True),
                "chr1",
                root / "output",
                mapq_min=20,
                baseq_min=20,
                thresholds=(1,),
                samtools_threads=0,
            )
            self.assertTrue(receipt["all_pass"])
            membership_path = root / "output" / "COLO829.chr1.site_component_membership.tsv.gz"
            with gzip.open(membership_path, "rt", encoding="utf-8", newline="") as handle:
                rows = [
                    row for row in csv.DictReader(handle, delimiter="\t")
                    if row["linkage_basis"] == "PS_HP1" and row["site_index"] == "2"
                ]
            self.assertEqual({row["phase_set"] for row in rows}, {"100", "200"})
            self.assertEqual(len(rows), 2)
            self.assertEqual(
                receipt["phase_set_contract_counts"]["known_phase_sets_by_hp_family"]["1"], 2
            )


if __name__ == "__main__":
    unittest.main()
