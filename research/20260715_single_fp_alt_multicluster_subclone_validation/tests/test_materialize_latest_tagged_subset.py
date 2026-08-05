#!/usr/bin/env python3
"""Synthetic contract tests for bounded latest-sidecar BAM materialization."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from array import array
from pathlib import Path

import pysam


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "materialize_latest_tagged_subset.py"
SPEC = importlib.util.spec_from_file_location("materialize_latest_tagged_subset", SCRIPT)
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


def make_alignment(name: str, start: int, hp: str, length: int = 200) -> pysam.AlignedSegment:
    alignment = pysam.AlignedSegment()
    alignment.query_name = name
    alignment.query_sequence = "C" * length
    alignment.flag = 0
    alignment.reference_id = 0
    alignment.reference_start = start
    alignment.mapping_quality = 60
    alignment.cigarstring = f"{length}M"
    alignment.query_qualities = pysam.qualitystring_to_array("I" * length)
    alignment.set_tag("HP", hp, value_type="Z")
    alignment.set_tag("PS", 1, value_type="i")
    alignment.set_tag("MM", "C+m,0;", value_type="Z")
    alignment.set_tag("ML", array("B", [250]))
    return alignment


class MaterializationContractTest(unittest.TestCase):
    def test_latest_tags_replace_embedded_tags_and_duplicate_window_exposure_collapses(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}]}
            with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
                bam.write(make_alignment("read_a", 100, "1", length=700))
                bam.write(make_alignment("read_b", 500, "2"))
            pysam.index(str(raw))
            with pysam.AlignmentFile(str(raw), "rb") as bam:
                records = list(bam.fetch(until_eof=True))

            sidecar_plain = root / "tags.tsv"
            with sidecar_plain.open("w", encoding="ascii") as handle:
                handle.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
                latest = [("1-1", "101"), ("2-1", "202")]
                for record, (hp, ps) in zip(records, latest):
                    digest = hashlib.blake2b((record.cigarstring or "*").encode(), digest_size=8).hexdigest()
                    handle.write(
                        f"{record.reference_name}\t{record.reference_start}\t{record.reference_end}\t"
                        f"{record.query_name}\t{record.flag}\t{record.mapping_quality}\t{digest}\t{hp}\t{ps}\n"
                    )
            sidecar = root / "tags.tsv.gz"
            pysam.tabix_compress(str(sidecar_plain), str(sidecar), force=True)
            pysam.tabix_index(str(sidecar), seq_col=0, start_col=1, end_col=2, zerobased=True, force=True)
            bed = root / "windows.bed"
            bed.write_text("chr1\t50\t400\nchr1\t450\t800\n", encoding="ascii")
            workspace = root / "workspace"
            entry = {
                "sample": "SYNTH",
                "raw_alignment": {"path": str(raw), "index_path": str(raw) + ".bai"},
                "latest_read_tag_sidecar": {"path": str(sidecar), "index_path": str(sidecar) + ".tbi"},
                "window_bed": {"path": str(bed)},
            }
            receipt = M.materialize_sample(entry, str(workspace), compression_threads=1)
            self.assertTrue(receipt["pass"])
            self.assertEqual(receipt["diagnostics"]["sidecar_missing"], 0)
            self.assertEqual(receipt["diagnostics"]["written_unique_alignments"], 2)
            self.assertEqual(receipt["diagnostics"]["duplicate_identity_collapsed"], 1)
            with pysam.AlignmentFile(receipt["output_bam"], "rb") as bam:
                output = list(bam.fetch(until_eof=True))
            self.assertEqual([record.get_tag("HP") for record in output], ["1-1", "2-1"])
            self.assertEqual([record.get_tag("PS") for record in output], [101, 202])
            self.assertTrue(all(record.has_tag("MM") and record.has_tag("ML") for record in output))

    def test_conflicting_sidecar_identity_is_rejected(self) -> None:
        key = ("q", "chr1", 10, 20, 0, "abc")
        values = {key: ("1", "1")}
        self.assertNotEqual(values[key], ("2", "2"))

    def test_nonanalytic_aux_difference_collapses(self) -> None:
        first = make_alignment("read", 100, "1")
        first.set_tag("cm", 10)
        second = make_alignment("read", 100, "1")
        second.set_tag("cm", 11)
        selected, reason, differences = M.resolve_duplicate(first, second)
        self.assertIs(selected, first)
        self.assertEqual(reason, "nonanalytic_aux_tag_variation")
        self.assertEqual(differences, ["cm"])

    def test_complete_methylation_payload_is_preferred(self) -> None:
        complete = make_alignment("read", 100, "1")
        incomplete = make_alignment("read", 100, "1")
        incomplete.set_tag("MM", None)
        incomplete.set_tag("ML", None)
        selected, reason, differences = M.resolve_duplicate(incomplete, complete)
        self.assertIs(selected, complete)
        self.assertEqual(reason, "prefer_complete_methylation_payload")
        self.assertEqual(differences, ["ML", "MM"])

    def test_conflicting_complete_methylation_payload_is_rejected(self) -> None:
        first = make_alignment("read", 100, "1")
        second = make_alignment("read", 100, "1")
        second.set_tag("MM", "C+m,1;", value_type="Z")
        with self.assertRaisesRegex(RuntimeError, "material duplicate conflict"):
            M.resolve_duplicate(first, second)


if __name__ == "__main__":
    unittest.main()
