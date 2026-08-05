#!/usr/bin/env python3
"""Synthetic end-to-end tests for the lossless ClairS-to-sSNV site ledger."""

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pysam


SCRIPT = Path(__file__).resolve().parent / "build_ssnv_site_ledger.py"


def write_vcf(path, records):
    header = [
        "##fileformat=VCFv4.2",
        "##source=ClairS synthetic contract fixture",
        "##contig=<ID=chr1,length=10000>",
        "##contig=<ID=chrX,length=10000>",
        "##FILTER=<ID=LowQual,Description=synthetic non-PASS>",
        "##FILTER=<ID=LPSLowQual,Description=synthetic LongPhase-S recalibration exclusion>",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=synthetic AF>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL",
    ]
    lines = list(header)
    for chrom, pos, ref, alt, filt in records:
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\tAF=0.25\tGT\t0/1\t0/0")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class SiteLedgerContractTest(unittest.TestCase):
    def test_every_raw_record_has_one_auditable_disposition(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = [
                ("chr1", 100, "A", "C", "PASS"),
                ("chr1", 150, "G", "T", "PASS"),
                ("chr1", 1000, "C", "A", "PASS"),
                ("chr1", 2000, "T", "G", "PASS"),
                ("chr1", 2010, "A", "G", "PASS"),
                ("chr1", 2020, "C", "T", "PASS"),
                ("chr1", 5000, "A", "T", "LowQual"),
                ("chr1", 6000, "A", "AT", "PASS"),
                ("chr1", 7000, "T", "C", "PASS"),
                ("chrX", 100, "G", "A", "PASS"),
            ]
            raw = root / "raw.vcf"
            caller_pass = root / "pass.vcf"
            recalibrated = root / "recalibrated.vcf"
            write_vcf(raw, records)
            longphase_records = [record for record in records if record[-1] == "PASS"]
            write_vcf(caller_pass, longphase_records)
            recalibrated_records = [
                (*record[:-1], "LPSLowQual" if record[1] == 7000 else "PASS")
                for record in longphase_records
            ]
            write_vcf(recalibrated, recalibrated_records)
            tree_input = root / "longphase_pass.vcf"
            write_vcf(tree_input, [record for record in recalibrated_records if record[-1] == "PASS"])
            mlhp = {
                "schema_version": "2.0",
                "params": {"TIER_R": 50, "MAX_SNV": 2},
                "groups": [{
                    "chrom": "chr1",
                    "positions": [100, 150],
                    "col_coverage": {"100": [7, 3], "150": [6, 4]},
                    "col_coverage_by_hp": {"1": {"100": [3, 2], "150": [2, 3]}},
                    "raw_HP_counts": {"1": 5, "1-1": 2, ".": 3},
                    "raw_HP_with_PS_counts": {"1": 5, "1-1": 2},
                    "n_unique_phase_sets": 2,
                    "phase_set_counts": {"100": 5, "200": 2},
                    "phase_set_HP_counts": {"100|1": 5, "200|1-1": 2},
                }],
            }
            mlhp_path = root / "mlhp_part_1.json"
            mlhp_path.write_text(json.dumps(mlhp), encoding="utf-8")
            output_tsv = root / "ledger.tsv.gz"
            output_summary = root / "summary.json"
            proc = subprocess.run(
                [sys.executable, str(SCRIPT), "--sample", "TEST", "--caller-raw-vcf", str(raw),
                 "--longphase-input-vcf", str(caller_pass), "--tree-input-vcf", str(tree_input),
                 "--recalibrated-vcf", str(recalibrated),
                 "--mlhp-glob", str(mlhp_path), "--output-tsv-gz", str(output_tsv),
                 "--output-summary", str(output_summary)],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
            summary = json.loads(output_summary.read_text(encoding="utf-8"))
            self.assertTrue(summary["pass"])
            self.assertTrue(Path(summary["output_index"]).is_file())
            self.assertEqual(summary["raw_clairs_records"], 10)
            self.assertEqual(summary["longphase_input_records"], 9)
            self.assertEqual(summary["longphase_recalibrated_records"], 9)
            self.assertEqual(summary["tree_input_records"], 8)
            self.assertEqual(summary["branch_counts"], {
                "retained": 2,
                "positional_singleton": 1,
                "read_unsupported": 2,
                "max_snv_excluded": 1,
                "unsupported_non_biallelic": 1,
                "out_of_scope_non_autosomal": 1,
                "excluded_before_longphase_nonPASS": 1,
                "excluded_by_longphase_filter": 1,
            })
            with gzip.open(output_tsv, "rt", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 10)
            self.assertEqual(len({(row["chrom"], row["pos"], row["ref"], row["alt"]) for row in rows}), 10)
            retained = next(row for row in rows if row["pos"] == "100")
            self.assertEqual(retained["ssnv_branch"], "retained")
            self.assertEqual(json.loads(retained["raw_HP_counts_json"])["1-1"], 2)
            self.assertEqual(json.loads(retained["phase_set_counts_json"]), {"100": 5, "200": 2})
            self.assertEqual(json.loads(retained["phase_set_HP_counts_json"])["200|1-1"], 2)
            self.assertEqual(json.loads(retained["sample_values_json"])["TUMOR"]["GT"], [0, 1])
            self.assertEqual(json.loads(retained["longphase_sample_values_json"])["TUMOR"]["GT"], [0, 1])
            self.assertEqual(json.loads(retained["longphase_info_json"])["AF"], [0.25])
            nonpass = next(row for row in rows if row["pos"] == "5000")
            self.assertEqual(nonpass["longphase_recalibrated_filter"], "NOT_INPUT")
            self.assertEqual(nonpass["longphase_filter_transition"], "LowQual->NOT_INPUT")
            lps_filtered = next(row for row in rows if row["pos"] == "7000")
            self.assertEqual(lps_filtered["longphase_recalibrated_filter"], "LPSLowQual")
            self.assertEqual(lps_filtered["longphase_filter_transition"], "PASS->LPSLowQual")
            self.assertEqual(lps_filtered["ssnv_branch"], "excluded_by_longphase_filter")
            with pysam.TabixFile(str(output_tsv)) as indexed:
                fetched = list(indexed.fetch("chr1", 99, 100))
            self.assertEqual(len(fetched), 1)
            self.assertIn("\t100\t", fetched[0])

    def test_raw_all_input_tracks_rescue_remove_and_stable_transitions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_records = [
                ("chr1", 100, "A", "C", "LowQual"),
                ("chr1", 150, "G", "T", "PASS"),
                ("chr1", 200, "C", "A", "PASS"),
                ("chrX", 100, "T", "G", "LowQual"),
            ]
            recalibrated_records = [
                ("chr1", 100, "A", "C", "PASS"),
                ("chr1", 150, "G", "T", "LPSLowQual"),
                ("chr1", 200, "C", "A", "PASS"),
                ("chrX", 100, "T", "G", "LPSLowQual"),
            ]
            raw = root / "raw-all.vcf"
            recalibrated = root / "recalibrated.vcf"
            tree_input = root / "tree.vcf"
            write_vcf(raw, raw_records)
            write_vcf(recalibrated, recalibrated_records)
            write_vcf(tree_input, [recalibrated_records[0], recalibrated_records[2]])
            mlhp = {
                "schema_version": "2.0",
                "params": {"TIER_R": 500, "MAX_SNV": 2},
                "groups": [{
                    "chrom": "chr1",
                    "positions": [100, 200],
                    "col_coverage": {"100": [8, 2], "200": [7, 3]},
                    "col_coverage_by_hp": {"1": {"100": [4, 1], "200": [3, 2]}},
                    "raw_HP_counts": {"1": 5},
                    "raw_HP_with_PS_counts": {"1": 5},
                    "n_unique_phase_sets": 1,
                    "phase_set_counts": {"100": 5},
                    "phase_set_HP_counts": {"100|1": 5},
                }],
            }
            mlhp_path = root / "mlhp_part_1.json"
            mlhp_path.write_text(json.dumps(mlhp), encoding="utf-8")
            output_tsv = root / "ledger.tsv.gz"
            output_summary = root / "summary.json"
            proc = subprocess.run(
                [sys.executable, str(SCRIPT), "--sample", "RAW_ALL", "--caller-raw-vcf", str(raw),
                 "--longphase-input-vcf", str(raw), "--longphase-input-contract", "clairs_raw_all",
                 "--tree-input-vcf", str(tree_input), "--recalibrated-vcf", str(recalibrated),
                 "--mlhp-glob", str(mlhp_path), "--output-tsv-gz", str(output_tsv),
                 "--output-summary", str(output_summary)],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
            summary = json.loads(output_summary.read_text(encoding="utf-8"))
            self.assertTrue(summary["pass"])
            self.assertEqual(summary["longphase_input_contract"], "clairs_raw_all")
            self.assertTrue(summary["checks"]["raw_all_equals_longphase_input"])
            self.assertNotIn("raw_PASS_equals_longphase_input", summary["checks"])
            self.assertEqual(summary["branch_counts"], {
                "retained": 2,
                "excluded_by_longphase_filter": 2,
            })
            self.assertEqual(summary["filter_transition_counts"], {
                "LowQual->LPSLowQual": 1,
                "LowQual->PASS": 1,
                "PASS->LPSLowQual": 1,
                "PASS->PASS": 1,
            })
            with gzip.open(output_tsv, "rt", encoding="utf-8") as handle:
                rows = {int(row["pos"]): row for row in csv.DictReader(handle, delimiter="\t")
                        if row["chrom"] == "chr1"}
            self.assertEqual(rows[100]["ssnv_branch"], "retained")
            self.assertEqual(rows[100]["longphase_filter_transition"], "LowQual->PASS")
            self.assertEqual(rows[150]["ssnv_branch"], "excluded_by_longphase_filter")
            self.assertEqual(rows[150]["longphase_filter_transition"], "PASS->LPSLowQual")

    def test_raw_all_input_with_clairs_pass_tree_checks_caller_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_records = [
                ("chr1", 100, "A", "C", "LowQual"),
                ("chr1", 150, "G", "T", "PASS"),
                ("chr1", 200, "C", "A", "PASS"),
                ("chrX", 100, "T", "G", "LowQual"),
            ]
            recalibrated_records = [
                ("chr1", 100, "A", "C", "PASS"),
                ("chr1", 150, "G", "T", "LPSLowQual"),
                ("chr1", 200, "C", "A", "PASS"),
                ("chrX", 100, "T", "G", "LPSLowQual"),
            ]
            raw = root / "raw-all.vcf"
            recalibrated = root / "recalibrated.vcf"
            tree_input = root / "clairs-pass.vcf"
            write_vcf(raw, raw_records)
            write_vcf(recalibrated, recalibrated_records)
            write_vcf(tree_input, [raw_records[1], raw_records[2]])
            mlhp = {
                "schema_version": "2.0",
                "params": {"TIER_R": 500, "MAX_SNV": 2},
                "groups": [{
                    "chrom": "chr1",
                    "positions": [150, 200],
                    "col_coverage": {"150": [8, 2], "200": [7, 3]},
                    "col_coverage_by_hp": {"1": {"150": [4, 1], "200": [3, 2]}},
                    "raw_HP_counts": {"1": 5},
                    "raw_HP_with_PS_counts": {"1": 5},
                    "n_unique_phase_sets": 1,
                    "phase_set_counts": {"100": 5},
                    "phase_set_HP_counts": {"100|1": 5},
                }],
            }
            mlhp_path = root / "mlhp_part_1.json"
            mlhp_path.write_text(json.dumps(mlhp), encoding="utf-8")
            output_tsv = root / "ledger.tsv.gz"
            output_summary = root / "summary.json"
            proc = subprocess.run(
                [sys.executable, str(SCRIPT), "--sample", "RAW_ALL_CLAIRS_PASS",
                 "--caller-raw-vcf", str(raw), "--longphase-input-vcf", str(raw),
                 "--longphase-input-contract", "clairs_raw_all",
                 "--tree-contract", "clairs_PASS_input", "--tree-input-vcf", str(tree_input),
                 "--recalibrated-vcf", str(recalibrated), "--mlhp-glob", str(mlhp_path),
                 "--output-tsv-gz", str(output_tsv), "--output-summary", str(output_summary)],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
            summary = json.loads(output_summary.read_text(encoding="utf-8"))
            self.assertTrue(summary["pass"])
            self.assertTrue(summary["checks"]["raw_all_equals_longphase_input"])
            self.assertTrue(summary["checks"]["caller_raw_PASS_equals_tree_input"])
            self.assertNotIn("longphase_input_equals_tree_input", summary["checks"])
            self.assertEqual(summary["tree_input_records"], 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
