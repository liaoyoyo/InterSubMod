#!/usr/bin/env python3
"""Synthetic contract tests for exact BAM-alignment to external HP/PS sidecar joins."""

import argparse
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

import pysam


SCRIPT_DIR = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("sm_linkage_genomewide", SCRIPT_DIR / "sm_linkage_genomewide.py")
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)

EXPECTED_CASES = {
    "test_exact_sidecar_matches_direct_tagged_bam": (
        "same-QNAME primary and supplementary alignments remain distinct and direct BAM HP/PS "
        "equals joined sidecar HP/PS"
    ),
    "test_missing_alignment_is_detected": (
        "an omitted coordinate identity increments sidecar_missing and cannot silently supply HP/PS"
    ),
    "test_duplicate_identity_conflict_is_detected": (
        "different HP/PS values for one coordinate identity increment sidecar_conflicts and are not consumed"
    ),
    "test_redundant_identical_sidecar_row_collapses_without_conflict": (
        "repeated coordinate identities with identical HP/PS collapse to one read-tag subject without extra exposure"
    ),
    "test_duplicate_bam_identity_with_conflicting_allele_is_excluded": (
        "one coordinate identity yielding different alleles at the same sSNV is excluded and counted as an allele conflict"
    ),
}


def make_alignment(name, start, flag, sequence, hp=None, ps=None):
    alignment = pysam.AlignedSegment()
    alignment.query_name = name
    alignment.query_sequence = sequence
    alignment.flag = flag
    alignment.reference_id = 0
    alignment.reference_start = start
    alignment.mapping_quality = 60
    alignment.cigarstring = f"{len(sequence)}M"
    alignment.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    if hp is not None:
        alignment.set_tag("HP", hp, value_type="Z")
    if ps is not None:
        alignment.set_tag("PS", ps, value_type="i")
    return alignment


def write_bam(path, tagged):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    records = [
        make_alignment("shared_qname", 100, 0, "CAAAAAAAAA", "1-1" if tagged else None, 101 if tagged else None),
        make_alignment("shared_qname", 200, 2048, "GAAAAAAAAA", "2-1" if tagged else None, 202 if tagged else None),
    ]
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for record in records:
            bam.write(record)
    pysam.index(str(path))
    # Re-read the records so each fixture is associated with the BAM header.
    # Alignment identity deliberately includes reference_name, which is not
    # available on a detached AlignedSegment created without a header.
    with pysam.AlignmentFile(str(path), "rb") as bam:
        return list(bam.fetch(until_eof=True))


def write_allele_conflict_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    records = [
        make_alignment("duplicate_identity", 100, 0, "CAAAAAAAAA"),
        make_alignment("duplicate_identity", 100, 0, "AAAAAAAAAA"),
    ]
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for record in records:
            bam.write(record)
    pysam.index(str(path))
    with pysam.AlignmentFile(str(path), "rb") as bam:
        return list(bam.fetch(until_eof=True))


def write_sidecar(
    path, records, overrides=None, omit=None, duplicate_conflict=None, duplicate_identical=None
):
    overrides = overrides or {}
    omit = omit or set()
    plain = path.with_suffix("")
    with plain.open("w", encoding="utf-8") as handle:
        handle.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
        for record in records:
            key = M._alignment_key(record)
            if key in omit:
                continue
            if key in overrides:
                hp, ps = overrides[key]
            else:
                hp, ps = record.get_tag("HP"), str(record.get_tag("PS"))
            fields = [
                record.reference_name,
                str(record.reference_start),
                str(record.reference_end),
                record.query_name,
                str(record.flag),
                str(record.mapping_quality),
                M._cigar_digest(record),
                str(hp),
                str(ps),
            ]
            handle.write("\t".join(fields) + "\n")
            if duplicate_identical == key:
                handle.write("\t".join(fields) + "\n")
            if duplicate_conflict == key:
                fields[-2:] = ["4", "999"]
                handle.write("\t".join(fields) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), seq_col=0, start_col=1, end_col=2, zerobased=True, force=True)


def close_sidecar():
    if M._TAG_TABIX is not None:
        M._TAG_TABIX.close()
    M._TAG_TABIX = None


def observe(bam_path, sidecar=""):
    close_sidecar()
    M.READ_TAG_SIDECAR = str(sidecar) if sidecar else ""
    group = [(101, "A", "C", "TP"), (201, "A", "G", "TP")]
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        calls, hp, ps = M.per_read_calls(bam, "chr1", group)
    diagnostics = M.last_tag_diagnostics()
    normalized_ps = {key: str(value) for key, value in ps.items()}
    return calls, hp, normalized_ps, diagnostics


class ReadTagSidecarContractTest(unittest.TestCase):
    def tearDown(self):
        close_sidecar()
        M.READ_TAG_SIDECAR = ""

    def test_exact_sidecar_matches_direct_tagged_bam(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tagged = root / "tagged.bam"
            raw = root / "raw.bam"
            records = write_bam(tagged, tagged=True)
            write_bam(raw, tagged=False)
            sidecar = root / "tags.tsv.gz"
            write_sidecar(sidecar, records)

            direct = observe(tagged)
            joined = observe(raw, sidecar)
            self.assertEqual(direct[:3], joined[:3])
            self.assertEqual(joined[3]["alignment_group_exposures"], 2)
            self.assertEqual(joined[3]["sidecar_exact_matches"], 2)
            self.assertEqual(joined[3]["sidecar_missing"], 0)
            self.assertEqual(joined[3]["sidecar_conflicts"], 0)
            self.assertEqual(set(joined[1].values()), {"1-1", "2-1"})

    def test_missing_alignment_is_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tagged = root / "tagged.bam"
            raw = root / "raw.bam"
            records = write_bam(tagged, tagged=True)
            write_bam(raw, tagged=False)
            sidecar = root / "tags.tsv.gz"
            write_sidecar(sidecar, records, omit={M._alignment_key(records[1])})

            _calls, hp, _ps, diagnostics = observe(raw, sidecar)
            self.assertEqual(diagnostics["sidecar_exact_matches"], 1)
            self.assertEqual(diagnostics["sidecar_missing"], 1)
            self.assertEqual(len(hp), 1)

    def test_duplicate_identity_conflict_is_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tagged = root / "tagged.bam"
            raw = root / "raw.bam"
            records = write_bam(tagged, tagged=True)
            write_bam(raw, tagged=False)
            conflict_key = M._alignment_key(records[0])
            sidecar = root / "tags.tsv.gz"
            write_sidecar(sidecar, records, duplicate_conflict=conflict_key)

            _calls, hp, _ps, diagnostics = observe(raw, sidecar)
            self.assertEqual(diagnostics["sidecar_conflicts"], 1)
            self.assertNotIn(conflict_key, hp)

    def test_redundant_identical_sidecar_row_collapses_without_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tagged = root / "tagged.bam"
            raw = root / "raw.bam"
            records = write_bam(tagged, tagged=True)
            write_bam(raw, tagged=False)
            redundant_key = M._alignment_key(records[0])
            sidecar = root / "tags.tsv.gz"
            write_sidecar(sidecar, records, duplicate_identical=redundant_key)

            _calls, hp, _ps, diagnostics = observe(raw, sidecar)
            self.assertEqual(diagnostics["alignment_group_exposures"], 2)
            self.assertEqual(diagnostics["sidecar_exact_matches"], 2)
            self.assertEqual(diagnostics["sidecar_conflicts"], 0)
            self.assertEqual(diagnostics["alignment_identity_allele_conflicts"], 0)
            self.assertEqual(len(hp), 2)

    def test_duplicate_bam_identity_with_conflicting_allele_is_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.bam"
            records = write_allele_conflict_bam(raw)
            key = M._alignment_key(records[0])
            sidecar = root / "tags.tsv.gz"
            write_sidecar(
                sidecar,
                records,
                overrides={key: ("1", "100")},
                duplicate_identical=key,
            )

            calls, hp, _ps, diagnostics = observe(raw, sidecar)
            self.assertEqual(diagnostics["alignment_identity_allele_conflicts"], 1)
            self.assertNotIn(key, calls)
            self.assertNotIn(key, hp)

def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_new_text(path, text):
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--test-log", type=Path)
    args = parser.parse_args()
    if bool(args.receipt) != bool(args.test_log):
        parser.error("--receipt and --test-log must be provided together")

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ReadTagSidecarContractTest)
    observed_cases = {test._testMethodName for test in suite}
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    log_text = stream.getvalue()
    print(log_text, end="")

    exact_contract = observed_cases == set(EXPECTED_CASES) and result.testsRun == len(EXPECTED_CASES)
    success = result.wasSuccessful() and exact_contract
    if args.test_log:
        write_new_text(args.test_log.resolve(), log_text)
    if args.receipt and success:
        test_source = Path(__file__).resolve()
        consumer_source = (SCRIPT_DIR / "sm_linkage_genomewide.py").resolve()
        test_log = args.test_log.resolve()
        receipt = {
            "schema_version": "1.0",
            "scope": "SYNTHETIC_COORDINATE_JOIN_V1_CONTRACT",
            "claim_limit": "join edge-case behavior only; not production-data validation",
            "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "command_argv": [sys.executable, *sys.argv],
            "test_source": {"path": str(test_source), "sha256": file_sha256(test_source)},
            "consumer_source": {"path": str(consumer_source), "sha256": file_sha256(consumer_source)},
            "test_log": {"path": str(test_log), "sha256": file_sha256(test_log)},
            "tests": [
                {"name": name, "assertion": assertion, "pass": True}
                for name, assertion in EXPECTED_CASES.items()
            ],
            "tests_run": result.testsRun,
            "failures": len(result.failures),
            "errors": len(result.errors),
            "exit_code": 0,
            "pass": True,
        }
        write_new_text(
            args.receipt.resolve(),
            json.dumps(receipt, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        )
    if not exact_contract:
        print(
            f"expected exact join cases {sorted(EXPECTED_CASES)}, observed {sorted(observed_cases)}",
            file=sys.stderr,
        )
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
