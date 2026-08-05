from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pysam
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "latest_tag_join.py"
SPEC = importlib.util.spec_from_file_location("latest_tag_join", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

AUDIT_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_latest_tag_join_against_fp.py"
AUDIT_SPEC = importlib.util.spec_from_file_location("audit_latest_tag_join_against_fp", AUDIT_SCRIPT)
assert AUDIT_SPEC is not None and AUDIT_SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(AUDIT_SPEC)
sys.modules[AUDIT_SPEC.name] = AUDIT
AUDIT_SPEC.loader.exec_module(AUDIT)


def make_sidecar(tmp_path: Path, rows: list[str]) -> Path:
    plain = tmp_path / "tags.tsv"
    plain.write_text(
        "#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n" + "\n".join(rows) + "\n",
        encoding="ascii",
    )
    compressed = tmp_path / "tags.tsv.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), seq_col=0, start_col=1, end_col=2, zerobased=True, force=True)
    return compressed


def make_reads(tmp_path: Path, hp: str = "2") -> Path:
    path = tmp_path / "reads.tsv"
    path.write_text(
        "read_id\tread_name\tchr\tstart\tend\tmapq\thp\talt_support\tis_tumor\tstrand\n"
        f"0\tread-a\tchr1\t90\t150\t60\t{hp}\tALT\t1\t+\n",
        encoding="ascii",
    )
    return path


def test_latest_hp_ps_join_and_unphased_normalization(tmp_path: Path) -> None:
    sidecar = make_sidecar(tmp_path, ["chr1\t90\t150\tread-a\t0\t60\tcigar-a\t2\t101"])
    with pysam.TabixFile(str(sidecar)) as tabix:
        lookup = MODULE.fetch_site_lookup(tabix, "chr1", 100)
    joined = MODULE.join_read_rows(make_reads(tmp_path), lookup)
    assert joined[0][1] == MODULE.LatestTags(hp="2", ps=101)
    assert joined[0][2] == 1
    assert MODULE.normalize_hp(".") == "0"
    assert MODULE.normalize_ps(".") is None


def test_conflicting_projected_tags_are_rejected(tmp_path: Path) -> None:
    sidecar = make_sidecar(
        tmp_path,
        [
            "chr1\t90\t150\tread-a\t0\t60\tcigar-a\t1\t101",
            "chr1\t90\t150\tread-a\t0\t60\tcigar-b\t2\t102",
        ],
    )
    with pysam.TabixFile(str(sidecar)) as tabix, pytest.raises(RuntimeError, match="Ambiguous HP/PS"):
        MODULE.fetch_site_lookup(tabix, "chr1", 100)


def test_supplementary_conflict_is_excluded_like_intersubmod(tmp_path: Path) -> None:
    sidecar = make_sidecar(
        tmp_path,
        [
            "chr1\t90\t150\tread-a\t0\t60\tcigar-a\t1\t101",
            "chr1\t90\t150\tread-a\t2048\t60\tcigar-b\t2\t102",
        ],
    )
    with pysam.TabixFile(str(sidecar)) as tabix:
        lookup = MODULE.fetch_site_lookup(tabix, "chr1", 100)
    assert lookup.rows_fetched == 2
    assert lookup.rows_eligible == 1
    assert MODULE.join_read_rows(make_reads(tmp_path, hp="1"), lookup)[0][1].hp == "1"


def test_missing_read_identity_is_rejected(tmp_path: Path) -> None:
    sidecar = make_sidecar(tmp_path, ["chr1\t90\t150\tother-read\t0\t60\tcigar-a\t2\t101"])
    with pysam.TabixFile(str(sidecar)) as tabix:
        lookup = MODULE.fetch_site_lookup(tabix, "chr1", 100)
    with pytest.raises(KeyError, match="missing reads.tsv identity"):
        MODULE.join_read_rows(make_reads(tmp_path), lookup)


def test_fp_regression_audit_is_hp_only_when_reads_tsv_has_no_ps(tmp_path: Path) -> None:
    sidecar = make_sidecar(tmp_path, ["chr1\t90\t150\tread-a\t0\t60\tcigar-a\t2\t101"])
    fp_root = tmp_path / "fp"
    reads = fp_root / "S1" / "callset" / "chr1" / "chr1_100" / "window" / "reads" / "reads.tsv"
    reads.parent.mkdir(parents=True)
    reads.write_text(
        "read_id\tread_name\tchr\tstart\tend\tmapq\thp\talt_support\tis_tumor\tstrand\n"
        "0\tread-a\tchr1\t90\t150\t60\t2\tALT\t1\t+\n",
        encoding="ascii",
    )

    result = AUDIT.audit_sample(
        "S1",
        str(fp_root),
        str(sidecar),
        str(sidecar) + ".tbi",
    )

    assert result["pass"] is True
    assert result["audit_scope"] == "HP_ONLY"
    assert result["counts"]["hp_matches"] == 1
    assert result["counts"]["reads_tsv_sites_without_ps_column"] == 1
    assert result["ps_evidence"]["status"] == "NOT_EVALUATED_HP_ONLY_AUDIT"
    assert result["ps_evidence"]["downstream_full_hp_ps_join_required"] is True
    assert "not_ps_completeness" in result["pass_semantics"]
