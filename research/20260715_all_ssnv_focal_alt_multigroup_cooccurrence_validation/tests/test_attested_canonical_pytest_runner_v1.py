from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import stat
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    TOPIC_ROOT / "audit_notes" / "run_attested_canonical_pytest_v5.py"
)


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "run_attested_canonical_pytest_v5_test",
        MODULE_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load_module()


def test_command_requires_bound_python_fd(module: ModuleType) -> None:
    python_path = str(module.PYTHON.resolve())
    command = module.canonical_command(python_path)
    assert command[0] == python_path
    assert str(module.TEST_ROOT) in command
    assert f"--junitxml={module.RAW_XML_PATH}" in command
    assert module.BOUND_TEST_DIR.parent == module.TOPIC_ROOT
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "pass_fds=(runtime_fd, *(binding[0] for binding in test_bindings))" in source
    assert 'pytest_executable = f"/proc/self/fd/{runtime_fd}"' in source
    assert "executable=pytest_executable" in source
    assert "-c" in command
    assert "/dev/null" in command
    assert "--noconftest" in command
    assert command[2:4] == [
        "-X",
        f"pycache_prefix={module.PYTHON_CACHE_DIR}",
    ]


def test_reader_detects_path_replacement(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.py"
    displaced = tmp_path / "source.displaced.py"
    source.write_bytes(b"trusted = True\n")
    source.chmod(0o444)
    with module.BoundArtifactReader() as reader:
        reader.open(source, include_mode=True)
        os.rename(source, displaced)
        source.write_bytes(b"trusted = False\n")
        source.chmod(0o444)
        with pytest.raises(module.TestRunError, match="binding changed"):
            reader.require_paths_still_bound()


def test_trust_modules_execute_from_bound_bytes(module: ModuleType) -> None:
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "execute_module_from_bytes(" in source
    assert "release_bytes" in source
    assert "evidence_bytes" in source
    assert "importlib" not in source
    assert "spec_from_file_location" not in source


def test_writer_is_exclusive_readonly(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    output = tmp_path / "evidence.txt"
    record = module.write_exclusive_readonly(output, b"evidence\n")
    assert record["mode"] == "0o444"
    assert stat.S_IMODE(output.stat().st_mode) == 0o444
    with pytest.raises(FileExistsError):
        module.write_exclusive_readonly(output, b"replacement\n")


def test_manifest_records_raw_and_canonical_junit(module: ModuleType) -> None:
    artifact = {
        "path": "/tmp/artifact",
        "size_bytes": 1,
        "sha256": "a" * 64,
        "mode": "0o444",
    }

    class Evidence:
        @staticmethod
        def source_set_digest(_sources: object) -> str:
            return "b" * 64

    payload = module.build_manifest_payload(
        run_id="45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        started_at="2026-07-19T00:00:00+00:00",
        finished_at="2026-07-19T00:01:00+00:00",
        git_head="c" * 40,
        sources={"source": artifact},
        test_sources=[artifact],
        support={"support": artifact},
        runtime=artifact,
        git_record=artifact,
        public_key=artifact,
        junit=artifact,
        counts={
            "tests": 1,
            "failures": 0,
            "errors": 0,
            "skipped": 0,
        },
        stdout_record=artifact,
        stderr_record=artifact,
        raw_junit_record=artifact,
        evidence_module=Evidence,
        command=module.canonical_command("/tmp/python"),
        entrypoint={
            "schema_name": "intersubmod.bound_python_entrypoint",
            "schema_version": "1.0.0",
            "argv0": "/proc/self/fd/32",
            "script": artifact,
            "python": artifact,
            "isolated_mode": True,
            "no_user_site": True,
            "environment": module.ENVIRONMENT,
        },
        pytest_process={
            "argv0": "/tmp/python",
            "executable": "/proc/self/fd/31",
            "python": artifact,
        },
        inventory={"count": 1, "sha256": "d" * 64},
    )
    assert payload["junit"]["artifact"] == artifact
    assert payload["junit"]["testcase_inventory"]["count"] == 1
    assert payload["captured_output"]["raw_junit"] == artifact
    assert payload["runtime"]["pytest_process"]["executable"] == "/proc/self/fd/31"
    assert payload["checks"]["all_outputs_no_clobber_mode_0444"] is True


def test_formal_output_slots_include_raw_evidence_directory(
    module: ModuleType,
) -> None:
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "RAW_OUTPUT_DIR" in source
    assert "os.mkdir(RAW_OUTPUT_DIR, mode=0o700)" in source
    assert "os.mkdir(PYTHON_CACHE_DIR, mode=0o700)" in source
    assert 'f"pycache_prefix={PYTHON_CACHE_DIR}"' in source
    assert "BOUND_TEST_DIR / path.name" in source
    assert "write_exclusive_readonly(XML_PATH, xml_bytes)" in source
