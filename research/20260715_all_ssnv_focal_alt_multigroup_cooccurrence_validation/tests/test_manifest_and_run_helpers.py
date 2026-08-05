from __future__ import annotations

import importlib.util
from pathlib import Path


TOPIC = Path(__file__).resolve().parents[1]


def load(name: str):
    path = TOPIC / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MANIFEST = load("prepare_all_ssnv_manifest")
RUNNER = load("run_all_ssnv_intersubmod")


def test_dataset_contract_is_fixed() -> None:
    assert MANIFEST.DATASETS == RUNNER.DATASETS
    assert len(MANIFEST.DATASETS) == 7
    assert len(set(MANIFEST.DATASETS)) == 7


def test_validate_index_accepts_explicit_index(tmp_path: Path) -> None:
    data = tmp_path / "input.bam"
    index = tmp_path / "custom.bai"
    data.write_bytes(b"bam")
    index.write_bytes(b"index")
    assert MANIFEST.validate_index(data, index) == index


def test_sha256_is_stable(tmp_path: Path) -> None:
    path = tmp_path / "value.txt"
    path.write_text("InterSubMod\n", encoding="ascii")
    assert MANIFEST.sha256(path) == RUNNER.sha256(path)


def make_site_artifacts(root: Path, position: int) -> None:
    region = root / "sample" / "chr1" / f"chr1_{position}" / f"chr1_{position - 5}_{position + 5}"
    for relative in ("reads/reads.tsv", "methylation/methylation.csv", "distance/BERNOULLI/matrix.csv"):
        path = region / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("header\n", encoding="ascii")


def test_validate_output_requires_clean_log_and_tsv_status(tmp_path: Path) -> None:
    for position in (10, 20):
        make_site_artifacts(tmp_path, position)
    (tmp_path / "significance_summary.csv").write_text("site\nchr1:10\n", encoding="ascii")
    (tmp_path / "significance_statistics.txt").write_text("ok\n", encoding="ascii")
    (tmp_path / "region_stratification_status.tsv").write_text(
        "status\treason\nNOT_APPLICABLE_TUMOR_ONLY\tNORMAL_BAM_REQUIRED\n", encoding="ascii"
    )
    (tmp_path / "inter_sub_mod.log").write_text("[INFO] complete\n", encoding="ascii")
    assert RUNNER.validate_output(tmp_path, 2)["pass"] is True
    (tmp_path / "inter_sub_mod.log").write_text(
        "[ERROR] Region 1 (chr1:20) failed: test\n", encoding="ascii"
    )
    validation = RUNNER.validate_output(tmp_path, 2)
    assert validation["pass"] is False
    assert validation["region_failures_in_log"] == 1
