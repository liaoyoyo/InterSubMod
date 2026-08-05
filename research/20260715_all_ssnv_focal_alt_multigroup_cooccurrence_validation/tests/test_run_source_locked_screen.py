from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "run_source_locked_screen.py"
)
SPEC = importlib.util.spec_from_file_location("run_source_locked_screen", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


FAKE_PRODUCER = r'''#!/usr/bin/env python3
import argparse
import hashlib
import json
import sys
from pathlib import Path

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def artifact(path):
    path = path.resolve()
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": digest(path)}

parser = argparse.ArgumentParser()
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument("--mutate", type=Path)
args = parser.parse_args()
args.output_dir.mkdir()
if args.mutate:
    args.mutate.write_text("changed\n", encoding="utf-8")
summary = args.output_dir / "all_ssnv_summary.json"
summary.write_text("{}\n", encoding="utf-8")
receipt = {
    "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
    "status": "EXECUTION_PASS",
    "pass": True,
    "command": sys.argv,
    "execution": {"selected_samples": ["S1"]},
    "source_code": {"analyzer": artifact(Path(__file__))},
}
(args.output_dir / "run_manifest.json").write_text(json.dumps(receipt) + "\n", encoding="utf-8")
'''


def make_args(
    tmp_path: Path, producer: Path, output_dir: Path, *extra: str
):
    return MODULE.parse_args(
        [
            "--producer",
            str(producer),
            "--output-dir",
            str(output_dir),
            "--workdir",
            str(tmp_path),
            "--",
            sys.executable,
            str(producer),
            "--output-dir",
            str(output_dir),
            *extra,
        ]
    )


def test_source_lock_receipt_binds_start_end_child_and_thread_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    for key in MODULE.THREAD_ENV_KEYS:
        monkeypatch.setenv(key, "1")
    producer = tmp_path / "producer.py"
    producer.write_text(FAKE_PRODUCER, encoding="utf-8")
    output_dir = tmp_path / "output"
    receipt_path = MODULE.run(make_args(tmp_path, producer, output_dir))
    payload = json.loads(receipt_path.read_text(encoding="utf-8"))
    key = str(producer.resolve())
    assert payload["pass"] is True
    assert payload["source_identity_before"] == payload["source_identity_after"]
    assert payload["source_identity_before"][key]["sha256"] == MODULE.sha256(producer)
    assert payload["computational_thread_environment"] == {
        name: "1" for name in MODULE.THREAD_ENV_KEYS
    }
    assert payload["checks"] == {
        "child_exit_zero": True,
        "all_locked_sources_unchanged": True,
        "wrapper_unchanged": True,
        "child_receipt_pass": True,
        "child_receipt_command_exact": True,
        "child_receipt_analyzer_matches_start_lock": True,
    }


def test_source_lock_rejects_locked_dependency_mutation(tmp_path: Path) -> None:
    producer = tmp_path / "producer.py"
    producer.write_text(FAKE_PRODUCER, encoding="utf-8")
    dependency = tmp_path / "dependency.py"
    dependency.write_text("original\n", encoding="utf-8")
    output_dir = tmp_path / "output"
    args = make_args(
        tmp_path,
        producer,
        output_dir,
        "--mutate",
        str(dependency),
    )
    args.locked_file = [dependency]
    with pytest.raises(RuntimeError, match="locked source file changed"):
        MODULE.run(args)
    assert not (output_dir / "source_lock_receipt.json").exists()


def test_source_lock_refuses_existing_output(tmp_path: Path) -> None:
    producer = tmp_path / "producer.py"
    producer.write_text(FAKE_PRODUCER, encoding="utf-8")
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    with pytest.raises(FileExistsError, match="Refusing to reuse"):
        MODULE.run(make_args(tmp_path, producer, output_dir))
