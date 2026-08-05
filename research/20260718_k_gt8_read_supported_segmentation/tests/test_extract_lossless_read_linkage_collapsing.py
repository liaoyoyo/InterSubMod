from pathlib import Path
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from extract_lossless_read_linkage_collapsing import (  # noqa: E402
    SidecarTags,
    alignment_id,
    extract_one,
    parse_args,
    target_observation_digest,
)


BASE_TAGS = SidecarTags(hp="1", ps="42001")
BASE_INDICES = (3, 7, 11)
BASE_CODES = ("R", "A", "R")
BASE_QUALITIES = (31, 42, 27)
BASE_MAPQ = 60


def test_cli_accepts_explicit_samtools_executable(tmp_path):
    samtools = tmp_path / "samtools"
    manifest = tmp_path / "manifest.json"
    args = parse_args(
        [
            "--manifest",
            str(manifest),
            "--sample",
            "HCC1395",
            "--chrom",
            "chr22",
            "--output-dir",
            str(tmp_path / "output"),
            "--samtools",
            str(samtools),
        ]
    )
    assert args.samtools == str(samtools)


def observation_digest(
    *,
    tags=BASE_TAGS,
    indices=BASE_INDICES,
    codes=BASE_CODES,
    qualities=BASE_QUALITIES,
    mapq=BASE_MAPQ,
):
    return target_observation_digest(
        tags=tags,
        indices=indices,
        codes=codes,
        qualities=qualities,
        mapq=mapq,
    )


def test_identical_target_evidence_meets_duplicate_collapse_predicate():
    canonical_key = (
        "read-001",
        "chr1",
        100,
        250,
        0,
        "0123456789abcdef",
    )
    identity = bytes.fromhex(alignment_id(canonical_key))
    seen = {identity: observation_digest()}

    independently_recomputed = observation_digest(
        tags=SidecarTags(hp="1", ps="42001"),
        indices=[3, 7, 11],
        codes=["R", "A", "R"],
        qualities=[31, 42, 27],
        mapq=60,
    )

    assert seen[identity] == independently_recomputed


@pytest.mark.parametrize(
    ("changed_field", "overrides"),
    [
        ("HP", {"tags": SidecarTags(hp="2", ps="42001")}),
        ("PS", {"tags": SidecarTags(hp="1", ps="42002")}),
        ("site indices", {"indices": (3, 8, 11)}),
        ("calls", {"codes": ("R", "R", "R")}),
        ("base qualities", {"qualities": (31, 41, 27)}),
        ("MAPQ", {"mapq": 59}),
    ],
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_any_target_evidence_difference_prevents_duplicate_collapse(
    changed_field, overrides
):
    assert observation_digest(**overrides) != observation_digest(), changed_field


@pytest.mark.parametrize(
    "policy",
    [
        None,
        "keep_all_rows",
        "collapse_redundant_rows_without_verification",
    ],
)
def test_extract_one_rejects_unsupported_manifest_duplicate_policy_before_io(
    tmp_path, policy
):
    missing = tmp_path / "must-not-be-opened"
    entry = {
        "sample": "HCC1395",
        "alignment_payload": {"path": str(missing.with_suffix(".bam"))},
        "somatic": {"tree_vcf": {"path": str(missing.with_suffix(".vcf.gz"))}},
        "read_tags": {
            "sidecar": {"path": str(missing.with_suffix(".bed.gz"))},
            "index": {"path": str(missing.with_suffix(".bed.gz.tbi"))},
            "duplicate_identity_policy": policy,
        },
    }
    output_dir = tmp_path / "output"

    with pytest.raises(
        RuntimeError,
        match="unsupported read-tag duplicate identity policy",
    ):
        extract_one(
            entry,
            "chr1",
            output_dir,
            mapq_min=20,
            baseq_min=20,
            thresholds=(1, 2, 3, 5),
            samtools_threads=1,
        )

    assert not output_dir.exists()
