from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = (
    REPO_ROOT
    / "research/20260723_production_exact_ps_strict_read_linkage/scripts"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


QA = load_module("all7_region_qa_hardening", SCRIPT_ROOT / "qa_all7_region_html.py")
FINALIZER = load_module(
    "all7_region_finalizer",
    SCRIPT_ROOT / "finalize_all7_region_bundle.py",
)


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )


def write_sidecar(sidecar: Path, target: Path) -> None:
    sidecar.write_text(
        f"{FINALIZER.sha256_path(target)}  {target.name}\n",
        encoding="ascii",
    )


def create_bundle_fixture(tmp_path: Path) -> dict[str, Path]:
    root = tmp_path / "bundle"
    data_dir = root / "data"
    query_dir = root / "queries"
    qa_dir = root / "visual_qa"
    data_dir.mkdir(parents=True)
    query_dir.mkdir()
    qa_dir.mkdir()

    dataset_rows = [
        {"dataset": dataset, "all_pass": True}
        for dataset in FINALIZER.CANONICAL_DATASETS
    ]
    chromosome_rows = [
        {"dataset": dataset, "chrom": chrom, "all_pass": True}
        for dataset in FINALIZER.CANONICAL_DATASETS
        for chrom in FINALIZER.AUTOSOMES
    ]
    data = {
        "schema_name": "intersubmod.strict_region_all7_report_data",
        "schema_version": "1.1.0",
        "all_pass": True,
        "scope": {
            "task_type": "B_comprehensive_validation",
            "datasets": list(FINALIZER.CANONICAL_DATASETS),
            "chromosomes": list(FINALIZER.AUTOSOMES),
        },
        "checks": {
            "fixture_data_pass": True,
            "all_154_extraction_receipt_identities_verified": True,
        },
        "datasets": dataset_rows,
        "chromosomes": chromosome_rows,
        "inputs": {
            dataset: {
                "extraction_receipts": {
                    chrom: {
                        "path": f"/fixture/{dataset}/{chrom}/extraction/receipt.json",
                        "size_bytes": 1,
                        "sha256": "1" * 64,
                    }
                    for chrom in FINALIZER.AUTOSOMES
                }
            }
            for dataset in FINALIZER.CANONICAL_DATASETS
        },
    }
    data_path = data_dir / "all7_report_data.json"
    write_json(data_path, data)
    for name in FINALIZER.DATA_OUTPUT_NAMES - {data_path.name}:
        (data_dir / name).write_text("fixture\n", encoding="utf-8")
    data_receipt_path = data_dir / "receipt.json"
    write_json(
        data_receipt_path,
        {
            "schema_name": "intersubmod.strict_region_all7_report_data_receipt",
            "schema_version": "1.1.0",
            "all_pass": True,
            "checks": data["checks"],
            "outputs": {
                name: FINALIZER.file_identity(data_dir / name)
                for name in sorted(FINALIZER.DATA_OUTPUT_NAMES)
            },
        },
    )
    data_sidecar_path = data_dir / "receipt.json.sha256"
    write_sidecar(data_sidecar_path, data_receipt_path)

    artifact = {
        "surface": "report",
        "manifest": {
            "surface": "report",
            "title": "Fixture all-7 report",
            "charts": [{"id": "chart_a"}],
            "tables": [{"id": "table_a"}],
            "blocks": [
                {"type": "chart", "chartId": "chart_a"},
                {"type": "table", "tableId": "table_a"},
            ],
        },
        "snapshot": {
            "status": "ready",
            "datasets": {
                "dataset_summary": dataset_rows,
                "chromosome_detail": chromosome_rows,
            },
        },
    }
    artifact_path = root / "artifact.json"
    write_json(artifact_path, artifact)
    database_path = data_dir / "report.sqlite"
    database_path.write_bytes(b"fixture sqlite")
    for dataset_id in artifact["snapshot"]["datasets"]:
        (query_dir / f"{dataset_id}.sql").write_text(
            f'SELECT * FROM "{dataset_id}";\n',
            encoding="utf-8",
        )
    artifact_receipt_path = root / "artifact_receipt.json"
    artifact_outputs = {
        artifact_path.name: FINALIZER.file_identity(artifact_path),
        "data/report.sqlite": FINALIZER.file_identity(database_path),
        **{
            f"queries/{path.name}": FINALIZER.file_identity(path)
            for path in sorted(query_dir.glob("*.sql"))
        },
    }
    write_json(
        artifact_receipt_path,
        {
            "schema_name": "intersubmod.strict_region_all7_artifact_receipt",
            "schema_version": "1.0.0",
            "all_pass": True,
            "scope": {
                "task_type": "B_comprehensive_validation",
                "datasets": list(FINALIZER.CANONICAL_DATASETS),
                "chromosomes": list(FINALIZER.AUTOSOMES),
            },
            "checks": {"fixture_artifact_pass": True},
            "inputs": {
                "data": FINALIZER.file_identity(data_path),
                "data_receipt": FINALIZER.file_identity(data_receipt_path),
                "data_receipt_sidecar": FINALIZER.file_identity(
                    data_sidecar_path
                ),
            },
            "outputs": artifact_outputs,
        },
    )
    artifact_sidecar_path = root / "artifact_receipt.json.sha256"
    write_sidecar(artifact_sidecar_path, artifact_receipt_path)

    html_path = root / "report.html"
    html_path.write_text("<!doctype html><title>Fixture</title>\n", encoding="utf-8")
    browser_path = root / "fixture-chromium"
    browser_path.write_bytes(b"fixture browser")
    delivery_receipt_path = root / "delivery_receipt.json"
    write_json(
        delivery_receipt_path,
        {
            "version": 1,
            "ok": True,
            "input": {
                "path": str(artifact_path.resolve()),
                "bytes": artifact_path.stat().st_size,
                "sha256": FINALIZER.sha256_path(artifact_path),
            },
            "output": {
                "path": str(html_path.resolve()),
                "bytes": html_path.stat().st_size,
                "sha256": FINALIZER.sha256_path(html_path),
            },
            "browser": {
                "path": str(browser_path.resolve()),
                "bytes": browser_path.stat().st_size,
                "sha256": FINALIZER.sha256_path(browser_path),
            },
            "stages": {
                "canonicalValidation": "passed",
                "package": "passed",
                "staticSvgExtraction": "passed",
                "officialBrowserVerification": "passed",
                "atomicOutputRename": "passed",
            },
            "staticSvg": {
                "chartIds": ["chart_a"],
                "count": 1,
                "themeVariants": 2,
                "embeddedThemeVariants": 2,
            },
            "verification": {
                "counts": {"charts": 1},
                "sourceDialog": "passed",
                "viewports": [1440, 390],
            },
        },
    )

    evidence_files = {}
    for name in FINALIZER.EXPECTED_QA_EVIDENCE:
        path = qa_dir / name
        path.write_bytes(f"fixture {name}".encode())
        evidence_files[name] = FINALIZER.file_identity(path)
    qa_receipt_path = qa_dir / "visual_qa_receipt.json"
    qa_viewports = []
    for mode, label, width, height in sorted(FINALIZER.EXPECTED_QA_MODES):
        qa_viewports.append(
            {
                "mode": mode,
                "viewport": {
                    "label": label,
                    "width": width,
                    "height": height,
                },
                "all_pass": True,
                "checks": {"fixture_viewport_pass": True},
                "metrics": {"horizontalOverflowPx": 0},
            }
        )
    write_json(
        qa_receipt_path,
        {
            "schema_name": "intersubmod.strict_region_html_visual_qa",
            "schema_version": "2.0.0",
            "all_pass": True,
            "inputs": {
                "artifact": FINALIZER.file_identity(artifact_path),
                "html": FINALIZER.file_identity(html_path),
            },
            "browser": {
                "engine": "chromium",
                "executable": FINALIZER.file_identity(browser_path),
            },
            "expected": {
                "chart_ids": ["chart_a"],
                "table_ids": ["table_a"],
                "visible_charts": 1,
                "visible_tables": 1,
                "static_svg_variants": 2,
                "modes": ["javascript_enabled", "javascript_disabled"],
            },
            "evidence": {
                "all_expected_files_present": True,
                "files": evidence_files,
            },
            "viewports": qa_viewports,
        },
    )
    return {
        "root": root,
        "data": data_path,
        "data_receipt": data_receipt_path,
        "data_sidecar": data_sidecar_path,
        "artifact": artifact_path,
        "artifact_receipt": artifact_receipt_path,
        "artifact_sidecar": artifact_sidecar_path,
        "html": html_path,
        "delivery_receipt": delivery_receipt_path,
        "visual_qa_receipt": qa_receipt_path,
        "ready": root / "READY.json",
    }


def finalizer_args(paths: dict[str, Path]) -> list[str]:
    return [
        "--data",
        str(paths["data"]),
        "--data-receipt",
        str(paths["data_receipt"]),
        "--data-sidecar",
        str(paths["data_sidecar"]),
        "--artifact",
        str(paths["artifact"]),
        "--artifact-receipt",
        str(paths["artifact_receipt"]),
        "--artifact-sidecar",
        str(paths["artifact_sidecar"]),
        "--html",
        str(paths["html"]),
        "--delivery-receipt",
        str(paths["delivery_receipt"]),
        "--visual-qa-receipt",
        str(paths["visual_qa_receipt"]),
        "--output",
        str(paths["ready"]),
    ]


def test_visual_id_contract_rejects_count_preserving_id_swap():
    valid = {
        "manifest": {
            "charts": [{"id": "chart_a"}, {"id": "chart_b"}],
            "tables": [{"id": "table_a"}],
            "blocks": [
                {"type": "chart", "chartId": "chart_a"},
                {"type": "chart", "chartId": "chart_b"},
                {"type": "table", "tableId": "table_a"},
            ],
        }
    }
    assert QA.expected_visual_ids(valid) == (
        ["chart_a", "chart_b"],
        ["table_a"],
    )
    valid["manifest"]["blocks"][1]["chartId"] = "chart_c"
    with pytest.raises(ValueError, match="do not exactly match"):
        QA.expected_visual_ids(valid)


def test_finalizer_publishes_ready_last_and_refuses_overwrite(tmp_path: Path):
    paths = create_bundle_fixture(tmp_path)
    assert FINALIZER.main(finalizer_args(paths)) == 0
    assert paths["ready"].is_file()
    ready_sidecar = paths["ready"].with_name("READY.json.sha256")
    FINALIZER.verify_sidecar(ready_sidecar, paths["ready"])
    ready = FINALIZER.read_json(paths["ready"])
    assert ready["all_pass"] is True
    assert ready["scope"]["dataset_chromosome_rows"] == 154
    assert ready["visual_contract"]["chart_ids"] == ["chart_a"]
    with pytest.raises(ValueError, match="refusing to overwrite"):
        FINALIZER.main(finalizer_args(paths))


def test_finalizer_rejects_tampered_html_without_ready(tmp_path: Path):
    paths = create_bundle_fixture(tmp_path)
    paths["html"].write_text("<!doctype html><title>Tampered</title>\n", encoding="utf-8")
    with pytest.raises(ValueError, match="file identity mismatch"):
        FINALIZER.main(finalizer_args(paths))
    assert not paths["ready"].exists()
    assert not paths["ready"].with_name("READY.json.sha256").exists()


def test_finalizer_requires_complete_javascript_disabled_qa_grid(tmp_path: Path):
    paths = create_bundle_fixture(tmp_path)
    receipt = FINALIZER.read_json(paths["visual_qa_receipt"])
    receipt["viewports"] = [
        row
        for row in receipt["viewports"]
        if not (
            row["mode"] == "javascript_disabled"
            and row["viewport"]["label"] == "mobile"
        )
    ]
    write_json(paths["visual_qa_receipt"], receipt)
    with pytest.raises(ValueError, match="grid is incomplete"):
        FINALIZER.main(finalizer_args(paths))
    assert not paths["ready"].exists()
    assert not paths["ready"].with_name("READY.json.sha256").exists()


def test_finalizer_rejects_missing_extraction_receipt_grid(tmp_path: Path):
    paths = create_bundle_fixture(tmp_path)
    data = FINALIZER.read_json(paths["data"])
    del data["inputs"]["HCC1395"]["extraction_receipts"]["chr1"]
    write_json(paths["data"], data)
    with pytest.raises(ValueError, match="not an exact chr1–22 grid"):
        FINALIZER.main(finalizer_args(paths))
    assert not paths["ready"].exists()
    assert not paths["ready"].with_name("READY.json.sha256").exists()
