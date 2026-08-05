from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "build_pattern_census.py"
)
SPEC = importlib.util.spec_from_file_location("build_pattern_census", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
census = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = census
SPEC.loader.exec_module(census)


SPARSE_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "qname_sha256",
    "read_group",
    "alignment_id",
    "start0",
    "end0",
    "flag",
    "mapq",
    "strand",
    "hp_raw",
    "hp_family",
    "phase_set",
    "site_indices",
    "positions1",
    "call_codes",
    "base_qualities",
    "n_sites_in_span",
    "n_fixed_ra",
    "n_alt",
)


def digest(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def topology_row(
    *,
    region_id: str,
    unit_id: str,
    phase_set: str,
    hp_family: str,
    active_positions: list[int],
) -> dict[str, object]:
    return {
        "schema_name": "intersubmod.exact_ps_cpp_topology_af.unit",
        "schema_version": "1.0.0",
        "sample": "H1437",
        "chrom": "chr1",
        "region_id": region_id,
        "unit_id": unit_id,
        "phase_set": phase_set,
        "hp_family": hp_family,
        "active_positions": active_positions,
        "active_bit_count": len(active_positions),
    }


def sparse_row(
    *,
    label: str,
    hp_raw: str,
    hp_family: str,
    phase_set: str,
    positions: tuple[int, ...],
    codes: str,
) -> dict[str, str]:
    return {
        "dataset": "H1437",
        "chrom": "chr1",
        "molecule_id": digest(f"molecule:{label}"),
        "qname_sha256": digest(f"qname:{label}"),
        "read_group": f"RG-{label[-1]}",
        "alignment_id": digest(f"alignment:{label}"),
        "start0": str(min(positions) - 20),
        "end0": str(max(positions) + 20),
        "flag": "0",
        "mapq": "60",
        "strand": "+" if int(digest(label)[0], 16) % 2 == 0 else "-",
        "hp_raw": hp_raw,
        "hp_family": hp_family,
        "phase_set": phase_set,
        "site_indices": ",".join(str(index) for index in range(len(positions))),
        "positions1": ",".join(str(position) for position in positions),
        "call_codes": codes,
        "base_qualities": ",".join("30" for _ in positions),
        "n_sites_in_span": str(len(positions)),
        "n_fixed_ra": str(sum(code in {"R", "A"} for code in codes)),
        "n_alt": str(codes.count("A")),
    }


def write_inputs(tmp_path: Path, *, malformed: bool = False) -> tuple[Path, Path]:
    topology_root = tmp_path / "topology"
    topology_path = topology_root / "samples" / "H1437" / "H1437.topology.jsonl"
    topology_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        topology_row(
            region_id="chr1|PS=P1|HP=1|UPAIR:B0001",
            unit_id="UPAIR",
            phase_set="P1",
            hp_family="1",
            active_positions=[100, 200],
        ),
        topology_row(
            region_id="chr1|PS=P1|HP=1|UPARTIAL:B0001",
            unit_id="UPARTIAL",
            phase_set="P1",
            hp_family="1",
            active_positions=[100, 200, 300],
        ),
        topology_row(
            region_id="chr1|PS=P2|HP=2|UTRIPLE:B0001",
            unit_id="UTRIPLE",
            phase_set="P2",
            hp_family="2",
            active_positions=[400, 500, 600],
        ),
        topology_row(
            region_id="chr1|PS=P3|HP=1|UZERO:B0001",
            unit_id="UZERO",
            phase_set="P3",
            hp_family="1",
            active_positions=[],
        ),
    ]
    topology_path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )

    sparse_root = tmp_path / "production"
    sparse_path = (
        sparse_root
        / "samples"
        / "H1437"
        / "chromosomes"
        / "chr1"
        / "extraction"
        / "H1437.chr1.molecule_sparse_calls.tsv.gz"
    )
    sparse_path.parent.mkdir(parents=True, exist_ok=True)
    sparse_rows: list[dict[str, str]] = []

    for index in range(60):
        if index < 10:
            codes = "RA"
        elif index < 20:
            codes = "AA"
        elif index < 40:
            codes = "RR"
        else:
            codes = "RL"
        sparse_rows.append(
            sparse_row(
                label=f"raw11-{index:02d}",
                hp_raw="1-1",
                hp_family="1",
                phase_set="P1",
                positions=(100, 200),
                codes=codes,
            )
        )
    for index in range(40):
        sparse_rows.append(
            sparse_row(
                label=f"raw1-{index:02d}",
                hp_raw="1",
                hp_family="1",
                phase_set="P1",
                positions=(100, 200),
                codes="RR",
            )
        )
    for index in range(40):
        if index < 8:
            codes = "RRA"
        elif index < 16:
            codes = "AAA"
        else:
            codes = "RRR"
        sparse_rows.append(
            sparse_row(
                label=f"raw21-{index:02d}",
                hp_raw="2-1",
                hp_family="2",
                phase_set="P2",
                positions=(400, 500, 600),
                codes=codes,
            )
        )
    sparse_rows.reverse()
    if malformed:
        sparse_rows[0]["call_codes"] = "R"

    with gzip.open(sparse_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SPARSE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(sparse_rows)
    return topology_root, sparse_root


def run_census(
    tmp_path: Path, output_name: str, *, malformed: bool = False
) -> Path:
    topology_root, sparse_root = write_inputs(tmp_path, malformed=malformed)
    output_dir = tmp_path / output_name
    census.build_pattern_census(
        topology_root=topology_root,
        hcc1395_sparse_root=tmp_path / "unused-hcc1395",
        production_sparse_root=sparse_root,
        output_dir=output_dir,
        samples=("H1437",),
        chromosomes=("chr1",),
        sort_binary="sort",
    )
    return output_dir


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def recursive_digests(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_exact_raw_hp_counts_partial_states_and_formal_gates(tmp_path: Path) -> None:
    output_dir = run_census(tmp_path, "out")
    rows = read_tsv(output_dir / "pattern_counts.tsv")
    by_key = {(row["region_id"], row["hp_raw"]): row for row in rows}

    pair_11 = by_key[("chr1|PS=P1|HP=1|UPAIR:B0001", "1-1")]
    assert pair_11["n_total"] == "60"
    assert pair_11["n_complete"] == "40"
    assert pair_11["n_partial"] == "20"
    assert json.loads(pair_11["state_count_json"]) == {
        "AA": 10,
        "RA": 10,
        "RR": 20,
    }
    assert json.loads(pair_11["partial_state_count_json"]) == {"RX": 20}
    assert pair_11["formal_n5"] == "true"
    assert pair_11["formal_n8"] == "true"
    assert pair_11["formal_n10"] == "true"
    assert pair_11["pair_full4"] == "false"
    assert pair_11["k_ge_3"] == "false"

    pair_1 = by_key[("chr1|PS=P1|HP=1|UPAIR:B0001", "1")]
    assert pair_1["n_total"] == "40"
    assert json.loads(pair_1["state_count_json"]) == {"RR": 40}
    assert pair_1["formal_n5"] == "false"
    assert "1-1" != "1"

    triple = by_key[("chr1|PS=P2|HP=2|UTRIPLE:B0001", "2-1")]
    assert triple["n_total"] == "40"
    assert triple["k_ge_3"] == "true"
    assert triple["formal_n5"] == "true"
    assert triple["formal_n8"] == "true"
    assert triple["formal_n10"] == "false"

    partial_triple = by_key[("chr1|PS=P1|HP=1|UPARTIAL:B0001", "1-1")]
    assert partial_triple["n_complete"] == "0"
    assert json.loads(partial_triple["partial_state_count_json"]) == {
        "AAX": 10,
        "RAX": 10,
        "RRX": 20,
        "RXX": 20,
    }
    assert partial_triple["formal_n5"] == "false"

    marker_rows = read_tsv(output_dir / "marker_universe.tsv")
    assert len(marker_rows) == 8
    assert not any(row["region_id"].endswith("UZERO:B0001") for row in marker_rows)


def test_candidate_shard_sorted_contract_and_deterministic_bytes(
    tmp_path: Path,
) -> None:
    first = run_census(tmp_path, "out-first")
    second = run_census(tmp_path, "out-second")
    assert recursive_digests(first) == recursive_digests(second)

    manifest = read_tsv(first / "candidate_read_join.manifest.tsv")
    assert len(manifest) == 1
    shard = first / manifest[0]["relative_path"]
    with gzip.open(shard, "rt", encoding="utf-8", newline="") as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "dataset",
        "chrom",
        "region_id",
        "unit_id",
        "phase_set",
        "hp_family",
        "hp_raw",
        "qname_sha256",
        "read_group",
        "strand",
        "mapq",
        "start0",
        "end0",
        "active_positions",
        "pattern",
        "complete_pattern",
        "n_active_bits",
    }
    assert required.issubset(candidate_rows[0])
    primary_keys = [
        (
            row["dataset"],
            row["chrom"],
            row["region_id"],
            row["hp_raw"],
            row["qname_sha256"],
        )
        for row in candidate_rows
    ]
    assert primary_keys == sorted(primary_keys)
    assert {row["complete_pattern"] for row in candidate_rows} == {"true", "false"}
    assert {row["hp_raw"] for row in candidate_rows} == {"1", "1-1", "2-1"}

    receipt = json.loads((first / "pattern_census.receipt.json").read_text())
    assert receipt["all_pass"] is True
    assert receipt["inputs"]["sparse_tsv_count"] == 1
    assert receipt["contracts"]["analysis_stratum"] == "exact_hp_raw"


def test_malformed_sparse_call_vector_fails_closed(tmp_path: Path) -> None:
    topology_root, sparse_root = write_inputs(tmp_path, malformed=True)
    try:
        census.build_pattern_census(
            topology_root=topology_root,
            hcc1395_sparse_root=tmp_path / "unused-hcc1395",
            production_sparse_root=sparse_root,
            output_dir=tmp_path / "bad-output",
            samples=("H1437",),
            chromosomes=("chr1",),
            sort_binary="sort",
        )
    except census.CensusContractError as exc:
        assert "length mismatch" in str(exc)
    else:
        raise AssertionError("malformed sparse vector did not fail closed")
    assert not (tmp_path / "bad-output").exists()


def test_cli_scope_resolution_supports_subset_and_all7() -> None:
    assert census.resolve_samples(all7=True, sample_values=None) == census.DATASETS
    assert census.resolve_samples(
        all7=False, sample_values=["H1437,HCC1395"]
    ) == ("H1437", "HCC1395")
    assert census.resolve_chromosomes(["chr22,chr1"]) == ("chr1", "chr22")
    try:
        census.resolve_samples(all7=False, sample_values=None)
    except census.CensusContractError as exc:
        assert "choose --all7" in str(exc)
    else:
        raise AssertionError("missing scope selection did not fail closed")


def test_formal_gate_uses_complete_n_and_pair_full4_uses_n5() -> None:
    region = census.TopologyRegion(
        dataset="H1437",
        chrom="chr1",
        region_id="PAIR",
        unit_id="U",
        phase_set="P",
        hp_family="1",
        active_positions=(100, 200),
    )
    partial_inflated = census.PatternAccumulator(
        region=region,
        hp_raw="1-1",
        n_total=60,
        complete=census.Counter({"RR": 10, "AA": 10}),
        partial=census.Counter({"RX": 40}),
    )
    assert partial_inflated.as_row()["formal_n5"] == "false"

    full4 = census.PatternAccumulator(
        region=region,
        hp_raw="1-1",
        n_total=40,
        complete=census.Counter({"RR": 25, "RA": 5, "AR": 5, "AA": 5}),
    )
    row = full4.as_row()
    assert row["formal_n5"] == "true"
    assert row["pair_full4"] == "true"


class PatternCensusContractTest(unittest.TestCase):
    def test_exact_raw_hp_counts_partial_states_and_formal_gates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            test_exact_raw_hp_counts_partial_states_and_formal_gates(
                Path(temporary)
            )

    def test_candidate_shard_sorted_contract_and_deterministic_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            test_candidate_shard_sorted_contract_and_deterministic_bytes(
                Path(temporary)
            )

    def test_malformed_sparse_call_vector_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            test_malformed_sparse_call_vector_fails_closed(Path(temporary))

    def test_cli_scope_resolution_supports_subset_and_all7(self) -> None:
        test_cli_scope_resolution_supports_subset_and_all7()

    def test_formal_gate_uses_complete_n_and_pair_full4_uses_n5(self) -> None:
        test_formal_gate_uses_complete_n_and_pair_full4_uses_n5()


if __name__ == "__main__":
    unittest.main()
