from __future__ import annotations

import ast
import hashlib
import importlib.util
import inspect
import copy
from collections import Counter
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import signal
import shlex
import subprocess
import sys

import pytest


for thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ[thread_env] = "1"


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
RESULT_ROOT = TOPIC_ROOT / "results"
FAILED_V14_VERIFICATION_RECEIPT = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260723_v14_signed_authority_r_metadata_only_relation_schema_mismatch"
    / "tumor_ref_source_receipt_promotion_verification.recovery.v14.json"
)
FAILED_V28_FINAL_DATASET = (
    AUDIT_ROOT
    / "failed_formal_runs"
    / "20260724_v28_signed_dataset_c_report_metric_key_order_mismatch"
    / "observed_downstream_outputs"
    / "all_ssnv_final_report_dataset_v5_m2v5_source_attested"
    / "final_report_dataset.json"
)
INPUT_MANIFEST = RESULT_ROOT / "all_ssnv_input_manifest.json"


def test_v30_regression_suite_executes_bound_source_fd() -> None:
    token = os.environ.get("INTERSUBMOD_REGRESSION_TEST_FD")
    canonical = AUDIT_ROOT / "schema_recovery_tests" / "test_tumor_ref_schema_recovery_v30.py"
    if token is None:
        probe_source = (
            AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v30.py"
        ).read_text(encoding="utf-8")
        assert "INTERSUBMOD_REGRESSION_TEST_FD" in probe_source
        assert 'f"/proc/self/fd/{test_fd}"' in probe_source
        assert "pass_fds=(python_fd, test_fd)" in probe_source
        return

    descriptor = int(token)
    expected_path = Path(os.environ["INTERSUBMOD_REGRESSION_TEST_CANONICAL_PATH"])
    opened = os.fstat(descriptor)
    live = os.stat(expected_path, follow_symlinks=False)
    source_token = f"/proc/self/fd/{descriptor}"
    cmdline = [
        os.fsdecode(value)
        for value in Path("/proc/self/cmdline").read_bytes().split(b"\0")
        if value
    ]
    assert expected_path == canonical
    assert Path(__file__).resolve(strict=True) == canonical
    assert source_token in cmdline
    assert sys.executable == os.environ["INTERSUBMOD_REGRESSION_PYTHON_ALIAS"]
    assert globals().get("_INTERSUBMOD_BOUND_SOURCE_FD") == descriptor
    assert (opened.st_dev, opened.st_ino, opened.st_size) == (
        live.st_dev,
        live.st_ino,
        live.st_size,
    )


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    return module


def basic_artifact(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    observed = os.stat(path, follow_symlinks=False)
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(observed.st_mode & 0o7777),
    }


def legacy_stat_artifact(path: Path) -> dict[str, object]:
    record = basic_artifact(path)
    observed = os.stat(path, follow_symlinks=False)
    record.update(
        {
            "device": observed.st_dev,
            "inode": observed.st_ino,
            "mtime_ns": observed.st_mtime_ns,
            "ctime_ns": observed.st_ctime_ns,
        }
    )
    return record


def generate_ed25519_keypair(openssl: Path, tmp_path: Path, stem: str):
    private = tmp_path / f"{stem}.private.pem"
    public = tmp_path / f"{stem}.public.pem"
    subprocess.run(
        [str(openssl), "genpkey", "-algorithm", "ED25519", "-out", str(private)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            str(openssl),
            "pkey",
            "-in",
            str(private),
            "-pubout",
            "-out",
            str(public),
        ],
        check=True,
        capture_output=True,
    )
    private.chmod(0o400)
    public.chmod(0o444)
    return private, public


def require_missing_review_path(error: FileNotFoundError, review_paths) -> None:
    assert error.filename is not None
    assert Path(error.filename) in set(review_paths)


class NoopWriterLease:
    def recheck(self):
        return None

    def close(self):
        return None


def configure_temporary_absence_contract(builder, validator, monkeypatch, tmp_path):
    result_root = tmp_path / "results"
    review_root = tmp_path / "reviews"
    workspace_root = tmp_path / "workspace"
    for path in (result_root, review_root, workspace_root):
        path.mkdir()
    bundle = result_root / "authority.v30.bundle"
    forbidden = result_root / "forbidden-output.json"
    monkeypatch.setattr(validator, "RESULT_ROOT", result_root)
    monkeypatch.setattr(validator, "REVIEW_ROOT", review_root)
    monkeypatch.setattr(validator, "WORKSPACE_ROOT", workspace_root)
    monkeypatch.setattr(validator, "AUTHORITY_BUNDLE", bundle)
    monkeypatch.setattr(validator, "AUTHORITY", bundle / "authority.json")
    monkeypatch.setattr(
        validator, "AUTHORITY_SIGNATURE", bundle / "authority.ed25519.sig"
    )
    monkeypatch.setattr(validator, "AUTHORITY_COMMIT", bundle / "commit.json")
    monkeypatch.setattr(validator, "CEREMONY_FORBIDDEN_OUTPUT_SLOTS", (forbidden,))
    monkeypatch.setattr(validator, "EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT", 1)
    monkeypatch.setattr(
        validator,
        "CEREMONY_STAGING_PATTERNS",
        (".tumor_ref_promotion_schema_recovery_authority.v30.bundle.staging.*",),
    )
    monkeypatch.setattr(
        validator,
        "_validate_ceremony_absence_contract",
        lambda: {"forbidden_output_slot_count": 1, "pass": True},
    )
    leases = builder.open_absence_directory_leases(validator)
    return bundle, forbidden, leases


@pytest.fixture(scope="module")
def modules():
    return {
        "validator": load_module(
            "schema_recovery_validator_v30_test",
            AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v30.py",
        ),
        "verifier": load_module(
            "schema_recovery_verifier_v30_test",
            AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v30.py",
        ),
        "replayer": load_module(
            "schema_recovery_replayer_v30_test",
            AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v30.py",
        ),
        "continuation": load_module(
            "schema_recovery_continuation_v30_test",
            AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v30.py",
        ),
        "builder": load_module(
            "schema_recovery_builder_v30_test",
            AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v30.py",
        ),
        "probe": load_module(
            "schema_recovery_probe_v30_test",
            AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v30.py",
        ),
        "release_authority": load_module(
            "schema_recovery_release_authority_test",
            TOPIC_ROOT / "scripts" / "release_source_authority.py",
        ),
        "dataset_builder": load_module(
            "schema_recovery_final_dataset_builder_test",
            AUDIT_ROOT / "build_all_ssnv_final_report_dataset_schema_recovery_v29.py",
        ),
        "report_builder": load_module(
            "schema_recovery_report_builder_test",
            AUDIT_ROOT / "build_all_ssnv_report_artifact_schema_recovery_v29.py",
        ),
        "finalizer": load_module(
            "schema_recovery_result_report_finalizer_v30_test",
            AUDIT_ROOT / "finalize_task_b_result_release_schema_recovery_v30.py",
        ),
    }


def v30_runner_prefix(replayer) -> bytes:
    runner_lines = replayer.COMPLETION_RUNNER.read_bytes().splitlines(keepends=True)
    assert len(runner_lines) >= replayer.EXPECTED_RUNNER_LINE_COUNT_MINIMUM
    assert runner_lines[replayer.RUNNER_FIRST_EXCLUDED_LINE - 1].strip() == b""
    return b"".join(runner_lines[: replayer.RUNNER_INCLUDED_LINE_END])


def bind_v30_archive_rebase_inputs(replayer):
    reader = replayer.BoundArtifactReader()
    try:
        _, result_archive_data, result_archive_record = reader.open(
            replayer.FAILED_V28_RESULT_KEY_ARCHIVE_RECORD
        )
        _, _, result_public_record = reader.open(
            replayer.FAILED_V28_RESULT_PUBLIC_KEY
        )
        _, report_archive_data, report_archive_record = reader.open(
            replayer.FAILED_V28_REPORT_KEY_ARCHIVE_RECORD
        )
        _, _, report_public_record = reader.open(
            replayer.FAILED_V28_REPORT_PUBLIC_KEY
        )
        _, report_private_record = reader.open_metadata_with_size(
            replayer.FAILED_V28_REPORT_PRIVATE_KEY
        )
        report_private_link_count = reader.stat_for(
            replayer.FAILED_V28_REPORT_PRIVATE_KEY
        ).st_nlink
    except BaseException:
        reader.close()
        raise
    return reader, {
        "result_archive_data": result_archive_data,
        "result_archive_record": result_archive_record,
        "result_public_record": result_public_record,
        "report_archive_data": report_archive_data,
        "report_archive_record": report_archive_record,
        "report_public_record": report_public_record,
        "report_private_record": report_private_record,
        "report_private_link_count": report_private_link_count,
    }


def test_v30_archive_rebase_is_exact_and_deterministic(modules):
    replayer = modules["replayer"]
    original = v30_runner_prefix(replayer)

    transformed, contract = replayer.archive_rebase_runner_prefix(original)

    assert hashlib.sha256(original).hexdigest() == replayer.EXPECTED_RUNNER_PREFIX_SHA256
    assert hashlib.sha256(transformed).hexdigest() == (
        replayer.EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256
    )
    assert contract["contract"] == "deterministic_three_path_archive_rebase_v1"
    assert contract["mapping_count"] == len(replayer.RUNNER_ARCHIVE_REBASE) == 3
    assert contract["private_key_bytes_read"] is False
    assert contract["host_live_key_paths_created"] is False
    assert contract["pass"] is True
    assert len(original.splitlines(keepends=True)) == 358
    assert len(transformed.splitlines(keepends=True)) == 358
    for mapping, (old_path, archive_path, role) in zip(
        contract["mappings"], replayer.RUNNER_ARCHIVE_REBASE, strict=True
    ):
        assert original.count(old_path) == 1
        assert original.count(archive_path) == 0
        assert transformed.count(old_path) == 0
        assert transformed.count(archive_path) == 1
        assert mapping == {
            "role": role,
            "original_path": os.fsdecode(old_path),
            "archive_path": os.fsdecode(archive_path),
            "replacement_count": 1,
        }


def test_v30_fixture_and_canonical_source_bindings_are_exact(modules):
    expected = {
        "validator": AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v30.py",
        "verifier": AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_recovery_v30.py",
        "replayer": AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v30.py",
        "continuation": (
            AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_recovery_v30.py"
        ),
        "builder": AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v30.py",
        "probe": AUDIT_ROOT / "probe_tumor_ref_schema_recovery_sources_v30.py",
    }
    assert {
        role: Path(modules[role].__file__).resolve(strict=True)
        for role in expected
    } == expected
    assert modules["validator"].AUTHORITY_ID == (
        "20260724_tumor_ref_schema_recovery_v30"
    )
    assert modules["verifier"].VERIFIER == expected["verifier"]
    assert modules["replayer"].REPLAYER == expected["replayer"]
    assert modules["continuation"].CONTINUATION_RUNNER == expected["continuation"]
    assert modules["builder"].CEREMONY_BUILDER == expected["builder"]
    assert modules["builder"].VALIDATOR == expected["validator"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("original_digest", "Original runner prefix identity drift"),
        ("missing_original", "cardinality drift"),
        ("duplicate_original", "cardinality drift"),
        ("preexisting_archive", "cardinality drift"),
        ("transformed_digest", "Archive-rebased runner prefix identity drift"),
    ),
)
def test_v30_archive_rebase_mutations_fail_closed(
    modules, monkeypatch, mutation, message
):
    replayer = modules["replayer"]
    original = v30_runner_prefix(replayer)
    old_path, archive_path, _ = replayer.RUNNER_ARCHIVE_REBASE[0]
    mutated = original

    if mutation == "original_digest":
        mutated = bytes((original[0] ^ 1,)) + original[1:]
    elif mutation == "missing_original":
        mutated = original.replace(old_path, b"/missing/failed-v28-key.pem", 1)
    elif mutation == "duplicate_original":
        mutated = original + b"\n" + old_path
    elif mutation == "preexisting_archive":
        mutated = original + b"\n" + archive_path
    elif mutation == "transformed_digest":
        monkeypatch.setattr(
            replayer, "EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256", "0" * 64
        )
    else:
        raise AssertionError(f"Unhandled mutation: {mutation}")

    if mutation not in {"original_digest", "transformed_digest"}:
        monkeypatch.setattr(
            replayer,
            "EXPECTED_RUNNER_PREFIX_SHA256",
            hashlib.sha256(mutated).hexdigest(),
        )
    with pytest.raises(replayer.GateError, match=message):
        replayer.archive_rebase_runner_prefix(mutated)


def test_v30_archive_rebase_inputs_bind_real_archives_without_private_bytes(modules):
    replayer = modules["replayer"]
    reader, inputs = bind_v30_archive_rebase_inputs(replayer)
    try:
        observed = replayer.validate_archive_rebase_inputs(**inputs)
        private_entry = reader._opened[
            replayer.FAILED_V28_REPORT_PRIVATE_KEY.resolve(strict=True)
        ]

        assert private_entry[1] is None
        assert observed["private_key_bytes_read"] is False
        assert observed["key_reuse_forbidden_preserved"] is True
        assert observed["archive_paths_remain_canonical"] is True
        assert observed["report_private_key_metadata_only"] == {
            "path": str(replayer.FAILED_V28_REPORT_PRIVATE_KEY),
            "mode": "0o400",
            "size_bytes": 290,
        }
        assert observed["report_private_key_link_count"] == 1
        assert observed["pass"] is True
    finally:
        reader.close()


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("result_archive_hash", "artifact identity drift"),
        ("result_public_mode", "artifact identity drift"),
        ("report_private_size", "artifact identity drift"),
        ("report_private_link_count", "artifact identity drift"),
        ("key_reuse_forbidden", "archive record contract drift"),
        ("private_key_bytes_read", "archive record contract drift"),
    ),
)
def test_v30_archive_rebase_input_mutations_fail_closed(
    modules, monkeypatch, mutation, message
):
    replayer = modules["replayer"]
    reader, bound_inputs = bind_v30_archive_rebase_inputs(replayer)
    inputs = copy.deepcopy(bound_inputs)
    try:
        if mutation == "result_archive_hash":
            inputs["result_archive_record"]["sha256"] = "0" * 64
        elif mutation == "result_public_mode":
            inputs["result_public_record"]["mode"] = "0o644"
        elif mutation == "report_private_size":
            inputs["report_private_record"]["size_bytes"] += 1
        elif mutation == "report_private_link_count":
            inputs["report_private_link_count"] += 1
        elif mutation in {"key_reuse_forbidden", "private_key_bytes_read"}:
            payload = json.loads(inputs["report_archive_data"])
            if mutation == "key_reuse_forbidden":
                payload["key_reuse_forbidden"] = False
            else:
                payload["private_key_bytes_read"] = True
            mutated_data = replayer.encode_json(payload)
            mutated_sha256 = hashlib.sha256(mutated_data).hexdigest()
            inputs["report_archive_data"] = mutated_data
            inputs["report_archive_record"]["sha256"] = mutated_sha256
            monkeypatch.setattr(
                replayer,
                "EXPECTED_FAILED_V28_REPORT_ARCHIVE_RECORD_SHA256",
                mutated_sha256,
            )
        else:
            raise AssertionError(f"Unhandled mutation: {mutation}")

        with pytest.raises(replayer.GateError, match=message):
            replayer.validate_archive_rebase_inputs(**inputs)
    finally:
        reader.close()


def test_v30_real_archived_state_raw_fails_rebased_passes_without_outputs(modules):
    replayer = modules["replayer"]
    original = v30_runner_prefix(replayer)
    transformed, _ = replayer.archive_rebase_runner_prefix(original)
    protected_paths = tuple(
        dict.fromkeys(
            (
                *replayer.FAILED_V28_LIVE_KEY_PATHS,
                *replayer.DOWNSTREAM_OUTPUT_SLOTS,
                replayer.OUTPUT_LOG,
                replayer.OUTPUT_RECEIPT,
                replayer.OUTPUT_SUCCESS_WITNESS,
            )
        )
    )
    before = {str(path): not os.path.lexists(path) for path in protected_paths}
    assert all(before.values()), {
        path: absent for path, absent in before.items() if not absent
    }

    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    bash_fd = os.open(replayer.BASH, flags)
    try:
        command = [f"/proc/self/fd/{bash_fd}", "-s"]
        raw = subprocess.run(
            command,
            input=original,
            cwd=REPO_ROOT,
            env=replayer.EXPECTED_ENVIRONMENT,
            pass_fds=(bash_fd,),
            capture_output=True,
            check=False,
            timeout=60,
        )
        rebased = subprocess.run(
            command,
            input=transformed,
            cwd=REPO_ROOT,
            env=replayer.EXPECTED_ENVIRONMENT,
            pass_fds=(bash_fd,),
            capture_output=True,
            check=False,
            timeout=60,
        )
    finally:
        os.close(bash_fd)

    expected_missing = os.fsdecode(replayer.RUNNER_ARCHIVE_REBASE[0][0])
    assert raw.returncode == 1
    assert raw.stdout == b""
    assert raw.stderr == f"Required file is missing or empty: {expected_missing}\n".encode()
    assert rebased.returncode == 0
    assert rebased.stdout == b""
    assert rebased.stderr == b""
    after = {str(path): not os.path.lexists(path) for path in protected_paths}
    assert after == before
    assert all(after.values())


def test_v30_worker_subprocess_receives_only_archive_rebased_stdin(modules):
    replayer = modules["replayer"]
    worker_tree = ast.parse(inspect.getsource(replayer.worker_main))
    run_calls = []
    transform_assignments = []
    for node in ast.walk(worker_tree):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            targets = {
                element.id
                for target in node.targets
                if isinstance(target, (ast.Name, ast.Tuple))
                for element in (
                    target.elts if isinstance(target, ast.Tuple) else (target,)
                )
                if isinstance(element, ast.Name)
            }
            if (
                "archive_rebased_runner_prefix" in targets
                and isinstance(node.value.func, ast.Name)
                and node.value.func.id == "archive_rebase_runner_prefix"
            ):
                transform_assignments.append(node)
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "subprocess"
            and node.func.attr == "run"
        ):
            input_values = [
                keyword.value for keyword in node.keywords if keyword.arg == "input"
            ]
            if input_values:
                run_calls.append((node, input_values))

    assert len(transform_assignments) == 1
    assert len(run_calls) == 1
    replay_call, input_values = run_calls[0]
    assert len(input_values) == 1
    assert isinstance(input_values[0], ast.Name)
    assert input_values[0].id == "archive_rebased_runner_prefix"
    assert transform_assignments[0].lineno < replay_call.lineno


@pytest.mark.parametrize(
    "name",
    (
        "EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256",
        "FAILED_V28_RESULT_KEY_ARCHIVE",
        "FAILED_V28_REPORT_KEY_ARCHIVE",
        "FAILED_V28_RESULT_PUBLIC_KEY",
        "FAILED_V28_RESULT_KEY_ARCHIVE_RECORD",
        "FAILED_V28_REPORT_PUBLIC_KEY",
        "FAILED_V28_REPORT_PRIVATE_KEY",
        "FAILED_V28_REPORT_KEY_ARCHIVE_RECORD",
        "EXPECTED_FAILED_V28_RESULT_PUBLIC_KEY_SHA256",
        "EXPECTED_FAILED_V28_REPORT_PUBLIC_KEY_SHA256",
        "EXPECTED_FAILED_V28_RESULT_ARCHIVE_RECORD_SHA256",
        "EXPECTED_FAILED_V28_REPORT_ARCHIVE_RECORD_SHA256",
        "RUNNER_ARCHIVE_REBASE",
        "FAILED_V28_LIVE_KEY_PATHS",
    ),
)
def test_v30_archive_rebase_critical_global_mutations_are_detected(
    modules, monkeypatch, name
):
    replayer = modules["replayer"]
    shadow = dict(vars(replayer))
    original = getattr(replayer, name)
    if isinstance(original, Path):
        mutated = Path(f"{original}.mutation")
    elif isinstance(original, str):
        mutated = f"{original}.mutation"
    elif isinstance(original, tuple):
        mutated = original + (original[0],)
    else:
        raise AssertionError(f"Unsupported critical-global type for {name}")

    assert replayer.critical_globals_equal(vars(replayer), shadow)
    monkeypatch.setattr(replayer, name, mutated)
    assert not replayer.critical_globals_equal(vars(replayer), shadow)


def test_v30_claims_never_call_transformed_bytes_an_exact_original_replay(modules):
    replayer = modules["replayer"]
    worker_tree = ast.parse(inspect.getsource(replayer.worker_main))
    claims_node = None
    for node in ast.walk(worker_tree):
        if not isinstance(node, ast.Assign) or not any(
            isinstance(target, ast.Name) and target.id == "receipt_payload"
            for target in node.targets
        ):
            continue
        assert isinstance(node.value, ast.Dict)
        for key, value in zip(node.value.keys, node.value.values, strict=True):
            if isinstance(key, ast.Constant) and key.value == "claims":
                claims_node = value
                break
    assert isinstance(claims_node, ast.Dict)
    claims = {
        key.value: value.value
        for key, value in zip(claims_node.keys, claims_node.values, strict=True)
        if isinstance(key, ast.Constant)
        and isinstance(key.value, str)
        and isinstance(value, ast.Constant)
        and isinstance(value.value, bool)
    }

    positive_line_replay_claims = {
        key
        for key, value in claims.items()
        if value is True and "runner_lines_1_358" in key
    }
    assert positive_line_replay_claims
    assert all(
        "archive_rebased" in key for key in positive_line_replay_claims
    ), positive_line_replay_claims
    assert claims["runner_prefix_archive_rebased_three_paths"] is True
    assert claims["original_completion_runner_end_to_end_pass_claimed"] is False
    assert claims["shell_pathname_gate_result_is_authoritative"] is False
    assert claims["downstream_authorized_by_this_receipt"] is False
    assert claims["runner_line_359_or_later_executed"] is False


def test_v30_human_claim_surfaces_qualify_archive_rebased_execution(modules):
    replayer = modules["replayer"]
    completed = subprocess.CompletedProcess([], 0, b"", b"")
    log = replayer.build_log(
        completed,
        completed,
        replayer.EXPECTED_RUNNER_PREFIX_SHA256,
        replayer.EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256,
    ).decode("utf-8")
    terminal_claim = next(
        line for line in log.splitlines() if line.startswith("PASS EVIDENCE ONLY:")
    )
    unqualified = {}
    module_claim = inspect.getdoc(replayer) or ""
    if "archive-rebased" not in module_claim.lower():
        unqualified["module_docstring"] = module_claim
    if "archive-rebased" not in terminal_claim.lower():
        unqualified["terminal_log_claim"] = terminal_claim

    assert (
        f"Runner prefix SHA-256: {replayer.EXPECTED_RUNNER_PREFIX_SHA256}" in log
    )
    assert (
        "Archive-rebased runner prefix SHA-256: "
        f"{replayer.EXPECTED_ARCHIVE_REBASED_RUNNER_PREFIX_SHA256}"
    ) in log
    assert "does not authorize downstream execution" in log
    assert not unqualified, unqualified


def load_v28_metric_regression_fixture():
    final = json.loads(FAILED_V28_FINAL_DATASET.read_text(encoding="utf-8"))
    manifest = json.loads(INPUT_MANIFEST.read_text(encoding="utf-8"))
    return final, manifest


def test_v29_report_metrics_accept_canonical_sorted_json_key_order(modules):
    report_builder = modules["report_builder"]
    final, manifest = load_v28_metric_regression_fixture()

    canonical_bytes = json.dumps(
        final,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    reparsed = json.loads(canonical_bytes)

    assert tuple(reparsed["funnel_metrics"]["per_sample"]) != report_builder.EXPECTED_DATASETS
    metrics = report_builder.validate_metrics(reparsed, manifest)
    assert tuple(metrics["per_sample"]) == report_builder.EXPECTED_DATASETS
    assert metrics["pooled"]["M1"]["numerator"] == 102_842
    assert metrics["pooled"]["M2"]["numerator"] == 919
    assert metrics["pooled"]["G1"]["numerator"] == 7
    assert metrics["pooled"]["G2"]["numerator"] == 0


@pytest.mark.parametrize(
    ("field", "missing_key", "extra_key"),
    (
        ("per_sample", "HCC1395", "EXTRA_SAMPLE"),
        ("truth_strata", "TP", "EXTRA_TRUTH"),
        ("sample_by_truth", "HCC1395|TP", "EXTRA_SAMPLE|TP"),
    ),
)
def test_v29_report_metrics_reject_missing_or_extra_keys(
    modules, field, missing_key, extra_key
):
    report_builder = modules["report_builder"]
    final, manifest = load_v28_metric_regression_fixture()

    missing = copy.deepcopy(final)
    missing["funnel_metrics"][field].pop(missing_key)
    with pytest.raises(report_builder.ReportContractError):
        report_builder.validate_metrics(missing, manifest)

    extra = copy.deepcopy(final)
    exemplar = next(iter(extra["funnel_metrics"][field].values()))
    extra["funnel_metrics"][field][extra_key] = copy.deepcopy(exemplar)
    with pytest.raises(report_builder.ReportContractError):
        report_builder.validate_metrics(extra, manifest)


def tumor_ref_primary_audit_reference(dataset_builder):
    receipt_path = (
        dataset_builder.CANONICAL_TASK_B_PATHS["tumor_ref_dir"] / "run_manifest.json"
    )
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    return receipt["inputs"]["primary_artifact_audit_pre"]


def test_v29_real_tumor_ref_v1_v6_primary_audit_lineage_is_exact(modules):
    dataset_builder = modules["dataset_builder"]
    legacy, lineage = dataset_builder.validate_tumor_ref_primary_audit_lineage(
        tumor_ref_primary_audit_reference(dataset_builder),
        dataset_builder.CANONICAL_TASK_B_PATHS["primary_artifact_audit_pre"],
    )
    assert legacy["schema_version"] == "1.0.0"
    assert lineage["same_data_plane_identity"] is True
    assert lineage["current_v6_is_formal_release_gate"] is True
    assert lineage["counts"] == {
        "stable_sites": 102_842,
        "assignment_records": 102_842,
        "primary_artifacts_expected": 308_526,
        "primary_artifacts_verified": 308_526,
    }
    assert lineage["artifact_set_sha256"] == (
        "195e77d571576f37debf1627edb57f9dc7edb935849bd0ae6e29b323b380b2ca"
    )
    assert lineage["pass"] is True


def test_v29_tumor_ref_lineage_uses_one_bound_fd_snapshot_per_audit(
    modules, monkeypatch
):
    dataset_builder = modules["dataset_builder"]
    reference = tumor_ref_primary_audit_reference(dataset_builder)
    legacy_path = dataset_builder.LEGACY_TUMOR_REF_PRIMARY_AUDIT.resolve()
    current_path = dataset_builder.CANONICAL_TASK_B_PATHS[
        "primary_artifact_audit_pre"
    ].resolve()
    protected = {legacy_path, current_path}
    real_open = dataset_builder.os.open
    real_read_text = Path.read_text
    opened = Counter()

    def counted_open(path, *args, **kwargs):
        resolved = Path(path).resolve()
        if resolved in protected:
            opened[resolved] += 1
        return real_open(path, *args, **kwargs)

    def reject_path_read_text(path, *args, **kwargs):
        if path.resolve() in protected:
            raise AssertionError("lineage audit must parse descriptor-bound bytes")
        return real_read_text(path, *args, **kwargs)

    monkeypatch.setattr(dataset_builder.os, "open", counted_open)
    monkeypatch.setattr(Path, "read_text", reject_path_read_text)
    _, lineage = dataset_builder.validate_tumor_ref_primary_audit_lineage(
        reference, current_path
    )
    assert opened == Counter({legacy_path: 1, current_path: 1})
    assert lineage["legacy_audit"] == reference
    assert lineage["same_data_plane_identity"] is True


def test_v29_legacy_tumor_ref_audit_is_in_gate_and_mutation_inventory(modules):
    continuation = modules["continuation"]
    legacy = continuation.LEGACY_PRIMARY_PRE.resolve()
    assert continuation.GATE_INPUT_PATHS["legacy_primary_pre"] == legacy
    gate_source = inspect.getsource(continuation.validate_gate_inputs)
    run_source = inspect.getsource(continuation.run)
    assert 'gates.open(\n                path,' in gate_source
    assert 'gate_validation["records"]' in run_source
    assert "PathMutationSentinel(" in run_source


def test_v29_historical_signed_runtime_projection_is_exact(modules):
    continuation = modules["continuation"]
    current_roles = set(continuation.REVIEWED_RUNTIME_SOURCE_CONTRACTS)
    historical_roles = set(continuation.HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES)
    recovery_only_roles = set(continuation.RECOVERY_ONLY_RUNTIME_SOURCE_ROLES)
    assert historical_roles.isdisjoint(recovery_only_roles)
    assert historical_roles | recovery_only_roles == current_roles
    records = {
        role: {"role": role, "sha256": hashlib.sha256(role.encode()).hexdigest()}
        for role in current_roles
    }
    projection = continuation.historical_signed_runtime_projection(records)
    assert set(projection) == historical_roles
    assert not set(projection) & recovery_only_roles
    assert set(records) == current_roles
    reordered = dict(reversed(list(records.items())))
    assert continuation.historical_signed_runtime_projection(reordered) == projection


@pytest.mark.parametrize(
    "mutation", ["historical_overlap", "missing_recovery_role", "unknown_historical_role"]
)
def test_v29_historical_runtime_role_partition_definition_mutations_fail_closed(
    modules, monkeypatch, mutation
):
    continuation = modules["continuation"]
    historical = set(continuation.HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES)
    recovery_only = set(continuation.RECOVERY_ONLY_RUNTIME_SOURCE_ROLES)
    if mutation == "historical_overlap":
        historical.add(next(iter(recovery_only)))
    elif mutation == "missing_recovery_role":
        recovery_only.pop()
    else:
        historical.add("unreviewed_historical_runtime")
    monkeypatch.setattr(
        continuation, "HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES", frozenset(historical)
    )
    monkeypatch.setattr(
        continuation, "RECOVERY_ONLY_RUNTIME_SOURCE_ROLES", frozenset(recovery_only)
    )
    records = {
        role: {"role": role}
        for role in continuation.REVIEWED_RUNTIME_SOURCE_CONTRACTS
    }
    with pytest.raises(
        continuation.ContinuationError,
        match="Historical/current runtime role partition drift",
    ):
        continuation.historical_signed_runtime_projection(records)


@pytest.mark.parametrize("mutation", ["omit", "inject", "alias"])
def test_v29_historical_runtime_partition_mutations_fail_closed(
    modules, mutation
):
    continuation = modules["continuation"]
    records = {
        role: {"role": role}
        for role in continuation.REVIEWED_RUNTIME_SOURCE_CONTRACTS
    }
    if mutation == "omit":
        records.pop(next(iter(continuation.HISTORICAL_SIGNED_RUNTIME_SOURCE_ROLES)))
    elif mutation == "inject":
        records["unreviewed_runtime"] = {"role": "unreviewed_runtime"}
    else:
        role = next(iter(continuation.RECOVERY_ONLY_RUNTIME_SOURCE_ROLES))
        records[f"alias_{role}"] = records.pop(role)
    with pytest.raises(
        continuation.ContinuationError,
        match="Historical/current runtime role partition drift",
    ):
        continuation.historical_signed_runtime_projection(records)


def test_v29_signed_projection_is_used_only_for_historical_payloads(modules):
    continuation = modules["continuation"]
    promotion_source = inspect.getsource(continuation.validate_promotion_chain)
    gate_source = inspect.getsource(continuation.validate_gate_inputs)
    terminal_source = inspect.getsource(continuation.validate_terminal_releases)
    assert "historical_signed_runtime_projection" in promotion_source
    assert "**historical_runtime_sources" in promotion_source
    assert 'chain["reviewed_runtime_sources"]' not in promotion_source
    assert "REVIEWED_RUNTIME_SOURCE_CONTRACTS.items()" in gate_source
    assert "set(reviewed_runtime_sources)" in terminal_source


def test_v29_current_reviewed_runtime_gate_registry_is_exact(modules):
    continuation = modules["continuation"]
    observed = continuation.validate_reviewed_runtime_gate_registry()
    expected = {
        role: contract["path"]
        for role, contract in continuation.REVIEWED_RUNTIME_SOURCE_CONTRACTS.items()
    }
    assert observed == expected
    assert len(expected) == 14
    assert "final_release_finalizer" not in continuation.GATE_INPUT_PATHS
    for role, contract in continuation.REVIEWED_RUNTIME_SOURCE_CONTRACTS.items():
        path = contract["path"]
        assert path == expected[role]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == contract["sha256"]


@pytest.mark.parametrize("mutation", ["missing", "alias", "swap", "duplicate"])
def test_v29_current_reviewed_runtime_gate_registry_mutations_fail_closed(
    modules, monkeypatch, mutation
):
    continuation = modules["continuation"]
    registry = dict(continuation.GATE_INPUT_PATHS)
    roles = sorted(continuation.RECOVERY_ONLY_RUNTIME_SOURCE_ROLES)
    if mutation == "missing":
        registry.pop(roles[0])
    elif mutation == "alias":
        registry["final_release_finalizer"] = registry[roles[1]]
    elif mutation == "swap":
        registry[roles[0]], registry[roles[1]] = (
            registry[roles[1]],
            registry[roles[0]],
        )
    else:
        registry["duplicate_recovery_builder"] = registry[roles[2]]
    monkeypatch.setattr(continuation, "GATE_INPUT_PATHS", registry)
    with pytest.raises(
        continuation.ContinuationError,
        match="Reviewed runtime gate path registry drift",
    ):
        continuation.validate_reviewed_runtime_gate_registry()


def test_v29_tumor_ref_primary_audit_rejects_noncanonical_legacy_anchor(modules):
    dataset_builder = modules["dataset_builder"]
    reference = dict(tumor_ref_primary_audit_reference(dataset_builder))
    reference["path"] = str(
        dataset_builder.CANONICAL_TASK_B_PATHS["primary_artifact_audit_pre"]
    )
    with pytest.raises(dataset_builder.ContractError, match="sole authorized lineage anchor"):
        dataset_builder.validate_tumor_ref_primary_audit_lineage(
            reference,
            dataset_builder.CANONICAL_TASK_B_PATHS["primary_artifact_audit_pre"],
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("current_not_formal", "not formal-release eligible"),
        ("count_drift", "counts differ"),
        ("artifact_set_drift", "artifact-set SHA-256 differs"),
        ("data_plane_drift", "audited data-plane identity differs"),
    ],
)
def test_v29_tumor_ref_primary_audit_lineage_mutations_fail_closed(
    modules, monkeypatch, mutation, message
):
    dataset_builder = modules["dataset_builder"]
    real_load_bound = dataset_builder.load_bound_json_artifact
    legacy_path = dataset_builder.LEGACY_TUMOR_REF_PRIMARY_AUDIT.resolve()
    current_path = dataset_builder.CANONICAL_TASK_B_PATHS[
        "primary_artifact_audit_pre"
    ].resolve()

    def mutated_load_bound(path, label, *, declared=None):
        payload, record = real_load_bound(path, label, declared=declared)
        payload = copy.deepcopy(payload)
        resolved = Path(path).resolve()
        if mutation == "current_not_formal" and resolved == current_path:
            payload["formal_task_b_release_eligible"] = False
        elif mutation == "count_drift" and resolved == legacy_path:
            payload["counts"]["stable_sites"] -= 1
        elif mutation == "artifact_set_drift" and resolved == current_path:
            payload["verification"]["artifact_set_sha256"] = "0" * 64
        elif mutation == "data_plane_drift" and resolved == current_path:
            payload["inputs"]["site_results"]["mtime_ns"] -= 1
        return payload, record

    monkeypatch.setattr(
        dataset_builder, "load_bound_json_artifact", mutated_load_bound
    )
    monkeypatch.setattr(
        dataset_builder,
        "_audit_data_artifact",
        lambda reference, *, label: dict(reference),
    )
    with pytest.raises(dataset_builder.ContractError, match=message):
        dataset_builder.validate_tumor_ref_primary_audit_lineage(
            tumor_ref_primary_audit_reference(dataset_builder), current_path
        )


def recovered_window_payload(started, finished):
    return {
        "started_at_utc": started,
        "finished_at_utc": finished,
        "created_at_utc": finished,
    }


def test_v29_recovered_primary_audit_split_window_accepts_exact_epochs(modules):
    dataset_builder = modules["dataset_builder"]
    observed = dataset_builder.validate_recovered_primary_artifact_audit_window(
        recovered_window_payload("2026-07-01T00:00:00+00:00", "2026-07-01T01:00:00+00:00"),
        recovered_window_payload("2026-07-01T03:00:00+00:00", "2026-07-01T04:00:00+00:00"),
        recovered_window_payload("2026-07-01T08:00:00+00:00", "2026-07-01T09:00:00+00:00"),
        {
            "tumor-REF": {
                "started_at_utc": "2026-07-01T01:30:00+00:00",
                "finished_at_utc": "2026-07-01T02:30:00+00:00",
            },
            "cooccurrence-preflight": {
                "started_at_utc": "2026-07-01T04:30:00+00:00",
                "finished_at_utc": "2026-07-01T05:00:00+00:00",
            },
            "cooccurrence": {
                "started_at_utc": "2026-07-01T05:30:00+00:00",
                "finished_at_utc": "2026-07-01T06:00:00+00:00",
            },
        },
    )
    assert observed["tumor_ref_uses_legacy_epoch_only"] is True
    assert observed["all_other_consumers_use_current_v6_epoch"] is True
    assert observed["pass"] is True


@pytest.mark.parametrize("mutation", ["missing_required_with_extra", "tumor_ref_late", "consumer_early"])
def test_v29_recovered_primary_audit_split_window_mutations_fail_closed(
    modules, mutation
):
    dataset_builder = modules["dataset_builder"]
    producers = {
        "tumor-REF": {
            "started_at_utc": "2026-07-01T01:30:00+00:00",
            "finished_at_utc": "2026-07-01T02:30:00+00:00",
        },
        "cooccurrence-preflight": {
            "started_at_utc": "2026-07-01T04:30:00+00:00",
            "finished_at_utc": "2026-07-01T05:00:00+00:00",
        },
        "cooccurrence": {
            "started_at_utc": "2026-07-01T05:30:00+00:00",
            "finished_at_utc": "2026-07-01T06:00:00+00:00",
        },
    }
    if mutation == "missing_required_with_extra":
        producers.pop("cooccurrence")
        producers["unrelated-extra"] = {
            "started_at_utc": "2026-07-01T06:00:00+00:00",
            "finished_at_utc": "2026-07-01T06:30:00+00:00",
        }
    elif mutation == "tumor_ref_late":
        producers["tumor-REF"]["finished_at_utc"] = "2026-07-01T03:30:00+00:00"
    else:
        producers["cooccurrence"]["started_at_utc"] = "2026-07-01T03:30:00+00:00"
    with pytest.raises(dataset_builder.ContractError):
        dataset_builder.validate_recovered_primary_artifact_audit_window(
            recovered_window_payload("2026-07-01T00:00:00+00:00", "2026-07-01T01:00:00+00:00"),
            recovered_window_payload("2026-07-01T03:00:00+00:00", "2026-07-01T04:00:00+00:00"),
            recovered_window_payload("2026-07-01T08:00:00+00:00", "2026-07-01T09:00:00+00:00"),
            producers,
        )


@pytest.mark.parametrize(
    "present_indices",
    [(), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)],
)
def test_v29_review_evidence_all_states_are_explicit(
    modules, monkeypatch, present_indices
):
    probe = modules["probe"]
    present = {probe.REVIEW_EVIDENCE_PATHS[index] for index in present_indices}
    monkeypatch.setattr(probe.os.path, "lexists", lambda path: path in present)
    if 0 < len(present) < len(probe.REVIEW_EVIDENCE_PATHS):
        with pytest.raises(RuntimeError, match="partially present"):
            probe.classify_review_evidence(probe.REVIEW_EVIDENCE_PATHS)
    else:
        state, label = probe.classify_review_evidence(probe.REVIEW_EVIDENCE_PATHS)
        assert sum(state.values()) == len(present)
        assert label == ("all_present" if present else "all_absent")


@pytest.mark.parametrize("slot_index", range(424))
def test_v29_every_canonical_forbidden_slot_fails_closed(
    modules, monkeypatch, tmp_path, slot_index
):
    builder = modules["builder"]
    validator = modules["validator"]
    original = validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS[slot_index]
    original_parents = {
        validator.RESULT_ROOT: "results",
        validator.REVIEW_ROOT: "reviews",
        validator.WORKSPACE_ROOT: "workspace",
    }
    roots = {name: tmp_path / name for name in original_parents.values()}
    for root in roots.values():
        root.mkdir()
    translated = roots[original_parents[original.parent]] / original.name
    monkeypatch.setattr(validator, "RESULT_ROOT", roots["results"])
    monkeypatch.setattr(validator, "REVIEW_ROOT", roots["reviews"])
    monkeypatch.setattr(validator, "WORKSPACE_ROOT", roots["workspace"])
    monkeypatch.setattr(validator, "CEREMONY_FORBIDDEN_OUTPUT_SLOTS", (translated,))
    monkeypatch.setattr(validator, "EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT", 1)
    monkeypatch.setattr(validator, "CEREMONY_STAGING_PATTERNS", ())
    monkeypatch.setattr(
        validator,
        "_validate_ceremony_absence_contract",
        lambda: {"forbidden_output_slot_count": 1, "pass": True},
    )
    leases = builder.open_absence_directory_leases(validator)
    translated.write_bytes(b"occupied")
    try:
        with pytest.raises(RuntimeError, match="forbidden namespace"):
            builder.require_ceremony_absence(validator, leases)
    finally:
        for lease in leases.values():
            lease.close()


def test_v29_unexpected_staging_after_probe_fails_closed(modules, monkeypatch, tmp_path):
    builder = modules["builder"]
    validator = modules["validator"]
    _, _, leases = configure_temporary_absence_contract(
        builder, validator, monkeypatch, tmp_path
    )
    unexpected = validator.RESULT_ROOT / (
        ".tumor_ref_promotion_schema_recovery_authority.v30.bundle.staging.injected"
    )
    unexpected.mkdir()
    try:
        with pytest.raises(RuntimeError, match="forbidden namespace"):
            builder.require_ceremony_absence(validator, leases)
    finally:
        for lease in leases.values():
            lease.close()


@pytest.mark.parametrize("pattern_index", range(28))
def test_v29_every_ceremony_staging_pattern_fails_closed(
    modules, monkeypatch, tmp_path, pattern_index
):
    builder = modules["builder"]
    validator = modules["validator"]
    pattern = validator.CEREMONY_STAGING_PATTERNS[pattern_index]
    _, _, leases = configure_temporary_absence_contract(
        builder, validator, monkeypatch, tmp_path
    )
    monkeypatch.setattr(validator, "CEREMONY_STAGING_PATTERNS", (pattern,))
    unexpected = validator.RESULT_ROOT / pattern.replace("*", "injected")
    unexpected.mkdir()
    try:
        with pytest.raises(RuntimeError, match="forbidden namespace"):
            builder.require_ceremony_absence(validator, leases)
    finally:
        for lease in leases.values():
            lease.close()


def test_v29_clean_namespace_and_exact_active_stage_pass(modules, monkeypatch, tmp_path):
    builder = modules["builder"]
    validator = modules["validator"]
    _, _, leases = configure_temporary_absence_contract(
        builder, validator, monkeypatch, tmp_path
    )
    stage_name = (
        ".tumor_ref_promotion_schema_recovery_authority.v30.bundle.staging.allowed"
    )
    try:
        assert builder.require_ceremony_absence(validator, leases)["pass"] is True
        (validator.RESULT_ROOT / stage_name).mkdir()
        result = builder.require_ceremony_absence(
            validator, leases, allowed_stage_name=stage_name
        )
        assert result["pass"] is True
        assert result["allowed_stage_name"] == stage_name
    finally:
        for lease in leases.values():
            lease.close()


@pytest.mark.parametrize("mutation", ["persist", "move_after_create"])
def test_v29_midscan_directory_generation_drift_fails_closed(
    modules, monkeypatch, tmp_path, mutation
):
    builder = modules["builder"]
    validator = modules["validator"]
    _, forbidden, leases = configure_temporary_absence_contract(
        builder, validator, monkeypatch, tmp_path
    )
    real_stat = builder.os.stat
    injected = False

    def stat_with_early_slot_injection(path, *args, **kwargs):
        nonlocal injected
        try:
            return real_stat(path, *args, **kwargs)
        except FileNotFoundError:
            if (
                not injected
                and path == forbidden.name
                and kwargs.get("dir_fd")
                == leases[validator.RESULT_ROOT].descriptor
            ):
                injected = True
                forbidden.write_bytes(b"injected after the first slot check")
                if mutation == "move_after_create":
                    forbidden.rename(tmp_path / "moved-after-first-scan.json")
            raise

    monkeypatch.setattr(builder.os, "stat", stat_with_early_slot_injection)
    try:
        with pytest.raises(
            RuntimeError,
            match="generation drift|forbidden namespace|mutation event",
        ):
            builder.require_ceremony_absence(validator, leases)
        assert injected is True
    finally:
        for lease in leases.values():
            lease.close()


def test_v29_terminal_absence_call_order_is_structurally_guarded(modules):
    builder = modules["builder"]
    assert "require_ceremony_absence" in builder.terminal_recheck.__code__.co_names
    assert "require_ceremony_absence" in builder.publish_staged_bundle.__code__.co_names
    assert "DirectoryMutationWatch" in builder.require_ceremony_absence.__code__.co_names
    source = inspect.getsource(builder.publish_staged_bundle)
    external = source.index("recheck_staged_bundle")
    final_absence = source.index("require_ceremony_absence", external)
    terminal_stage = source.index("terminal_stage_identity_recheck", final_absence)
    atomic_publish = source.index(
        "atomic_publish_after_terminal_contract", terminal_stage
    )
    assert external < final_absence < terminal_stage < atomic_publish
    assert source.count("require_ceremony_absence") == 1
    terminal_source = inspect.getsource(builder.atomic_publish_after_terminal_contract)
    last_parent_token = terminal_source.rindex("live_parent_identity")
    rename_syscall = terminal_source.rindex("_RENAMEAT2")
    assert last_parent_token < rename_syscall
    assert "rename_noreplace(" not in terminal_source
    assert "_acquire_recovery_output_writer_lease" in inspect.getsource(
        modules["verifier"]._verify
    )
    assert "acquire_recovery_output_writer_lease" in inspect.getsource(
        modules["replayer"].worker_main
    )
    continuation_main_source = inspect.getsource(modules["continuation"].main)
    assert "acquire_recovery_output_writer_lease" in continuation_main_source
    assert "should_record_continuation_incident" in continuation_main_source


@pytest.mark.parametrize(
    "mode,lease_acquired,expected_incidents",
    [
        ("--supervise-and-sign", False, 0),
        ("--supervise-and-sign", True, 1),
        ("--supervised-child", False, 0),
        ("--invalid-mode", False, 0),
    ],
)
def test_v29_only_lease_holding_supervisor_failure_records_incident(
    modules, monkeypatch, mode, lease_acquired, expected_incidents
):
    continuation = modules["continuation"]
    incident_calls = 0

    def acquire_or_fail():
        if lease_acquired:
            return 123
        raise continuation.ContinuationError("synthetic writer lease refusal")

    def fail_execution():
        raise continuation.ContinuationError("synthetic continuation failure")

    def record_incident(_error):
        nonlocal incident_calls
        incident_calls += 1

    def exit_process(code):
        raise SystemExit(code)

    monkeypatch.setattr(
        continuation.sys, "argv", [str(continuation.CONTINUATION_RUNNER), mode]
    )
    monkeypatch.setattr(
        continuation, "acquire_recovery_output_writer_lease", acquire_or_fail
    )
    monkeypatch.setattr(continuation, "supervise_and_sign", fail_execution)
    monkeypatch.setattr(continuation, "run", fail_execution)
    monkeypatch.setattr(continuation, "record_continuation_incident", record_incident)
    monkeypatch.setattr(continuation.os, "write", lambda *args, **kwargs: 0)
    monkeypatch.setattr(continuation.os, "_exit", exit_process)
    with pytest.raises(SystemExit) as error:
        continuation.main()
    assert error.value.code == 1
    assert incident_calls == expected_incidents


def historical_identity() -> dict[str, object]:
    payload = json.loads(
        (
            RESULT_ROOT
            / "tumor_ref_source_receipt_promotion_authorization.v3.json"
        ).read_text(encoding="utf-8")
    )
    return payload["evidence"]["focal_source_identity_transition"]["during_execution"]


def signed_chain_payload(name: str) -> dict[str, object]:
    return json.loads((RESULT_ROOT / name).read_text(encoding="utf-8"))


def test_historical_record_is_exact_eight_keys_without_link_count(modules):
    verifier = modules["verifier"]
    record = historical_identity()
    assert set(record) == set(verifier.HISTORICAL_EXECUTION_ARTIFACT_KEYS)
    assert len(record) == 8
    assert "link_count" not in record
    verifier._require_historical_execution_artifact(record, "historical record")


@pytest.mark.parametrize("mutation", ["missing", "extra_link_count", "extra_other"])
def test_historical_schema_mutations_fail_closed(modules, mutation):
    verifier = modules["verifier"]
    record = dict(historical_identity())
    if mutation == "missing":
        record.pop("ctime_ns")
    elif mutation == "extra_link_count":
        record["link_count"] = 1
    else:
        record["unexpected"] = "forbidden"
    with pytest.raises(verifier.VerificationError):
        verifier._require_historical_execution_artifact(record, mutation)


def test_current_source_identity_still_has_nine_keys_and_one_link(modules):
    verifier = modules["verifier"]
    source = RESULT_ROOT / "tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json"
    data = source.read_bytes()
    observed = os.stat(source, follow_symlinks=False)
    basic = {
        "path": str(source),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(observed.st_mode & 0o7777),
    }
    current = verifier._stat_artifact(basic, observed)
    assert set(current) == set(verifier.STAT_ARTIFACT_KEYS)
    assert len(current) == 9
    assert current["link_count"] == 1


@pytest.mark.parametrize(
    "name",
    [
        "tumor_ref_source_receipt_promotion_authorization.v3.json",
        "tumor_ref_source_receipt_promotion.v3.json",
    ],
)
def test_transition_aware_walker_accepts_signed_chain_and_binds_current(modules, name):
    replayer = modules["replayer"]
    payload = signed_chain_payload(name)
    current = payload["evidence"]["focal_source_identity_transition"]["current"]
    reader = replayer.BoundArtifactReader()
    try:
        relations = replayer.recursively_validate_artifact_relations(
            payload, reader, require_focal_transition=True
        )
    finally:
        reader.close()
    assert current["path"] in relations
    assert relations[current["path"]]["mode"] == "0o444"


@pytest.mark.parametrize(
    "mutation",
    [
        "historical_missing_key",
        "historical_extra_link_count",
        "historical_mode",
        "current_mode",
        "current_link_count",
        "current_ctime_equals_historical",
        "current_ctime_precedes_historical",
        "unchanged_fields",
        "observed_differences",
        "interpretation",
        "moved_context",
    ],
)
def test_transition_aware_walker_mutations_fail_closed(modules, mutation):
    replayer = modules["replayer"]
    payload = copy.deepcopy(
        signed_chain_payload("tumor_ref_source_receipt_promotion_authorization.v3.json")
    )
    transition = payload["evidence"]["focal_source_identity_transition"]
    if mutation == "historical_missing_key":
        transition["during_execution"].pop("ctime_ns")
    elif mutation == "historical_extra_link_count":
        transition["during_execution"]["link_count"] = 1
    elif mutation == "historical_mode":
        transition["during_execution"]["mode"] = "0o444"
    elif mutation == "current_mode":
        transition["current"]["mode"] = "0o664"
    elif mutation == "current_link_count":
        transition["current"]["link_count"] = 2
    elif mutation == "current_ctime_equals_historical":
        transition["current"]["ctime_ns"] = transition["during_execution"]["ctime_ns"]
    elif mutation == "current_ctime_precedes_historical":
        transition["current"]["ctime_ns"] = transition["during_execution"]["ctime_ns"] - 1
        transition["observed_differences"]["ctime_ns"]["current"] = transition[
            "current"
        ]["ctime_ns"]
    elif mutation == "unchanged_fields":
        transition["unchanged_fields"][0] = "ctime_ns"
    elif mutation == "observed_differences":
        transition["observed_differences"]["mode"]["during"] = "0o444"
    elif mutation == "interpretation":
        transition["interpretation"] = "historical mode ignored"
    elif mutation == "moved_context":
        payload["moved_transition"] = payload["evidence"].pop(
            "focal_source_identity_transition"
        )
    reader = replayer.BoundArtifactReader()
    try:
        with pytest.raises(replayer.GateError):
            replayer.recursively_validate_artifact_relations(
                payload, reader, require_focal_transition=True
            )
    finally:
        reader.close()


def test_transition_is_rejected_in_payload_without_transition_contract(modules):
    replayer = modules["replayer"]
    payload = signed_chain_payload(
        "tumor_ref_source_receipt_promotion_authorization.v3.json"
    )
    reader = replayer.BoundArtifactReader()
    try:
        with pytest.raises(replayer.GateError):
            replayer.recursively_validate_artifact_relations(payload, reader)
    finally:
        reader.close()


def v14_terminal_key_state_payload() -> dict[str, object]:
    receipt = json.loads(FAILED_V14_VERIFICATION_RECEIPT.read_text(encoding="utf-8"))
    return {
        "terminal_key_state": receipt["schema_recovery_authority"][
            "terminal_key_state"
        ]
    }


def test_v29_actual_v14_private_key_state_canary_and_fix(modules):
    payload = v14_terminal_key_state_payload()
    legacy_path = Path(
        payload["terminal_key_state"]["legacy_private_key"]["path"]
    )
    recovery_path = Path(
        payload["terminal_key_state"]["recovery_private_key"]["path"]
    )
    v15 = load_module(
        "schema_recovery_replayer_v15_metadata_relation_canary",
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v14.py",
    )
    v15_reader = v15.BoundArtifactReader()
    try:
        v15_reader.open_metadata(legacy_path)
        with pytest.raises(
            v15.GateError,
            match="Metadata-only artifact cannot be reopened for bytes",
        ):
            v15.recursively_validate_artifact_relations(payload, v15_reader)
    finally:
        v15_reader.close()

    v29 = modules["replayer"]
    v29_reader = v29.BoundArtifactReader()
    try:
        v29_reader.open_metadata(legacy_path)
        relations = v29.recursively_validate_artifact_relations(payload, v29_reader)
        assert relations[str(legacy_path)]["size_bytes"] == 119
        assert relations[str(recovery_path)]["size_bytes"] == 119
        assert v29_reader._opened[legacy_path.resolve(strict=True)][1] is None
        assert v29_reader._opened[recovery_path.resolve(strict=True)][1] is None
    finally:
        v29_reader.close()


@pytest.mark.parametrize(
    "mutation",
    ["missing_state", "extra_field", "wrong_mode", "boolean_size", "empty_state"],
)
def test_v29_metadata_plus_size_schema_mutations_fail_closed(modules, mutation):
    replayer = modules["replayer"]
    payload = v14_terminal_key_state_payload()
    record = payload["terminal_key_state"]["legacy_private_key"]
    if mutation == "missing_state":
        record.pop("state")
    elif mutation == "extra_field":
        record["unexpected"] = "forbidden"
    elif mutation == "wrong_mode":
        record["mode"] = "0o444"
    elif mutation == "boolean_size":
        record["size_bytes"] = True
    else:
        record["state"] = ""
    reader = replayer.BoundArtifactReader()
    try:
        with pytest.raises(replayer.GateError, match="metadata-plus-size schema drift"):
            replayer.recursively_validate_artifact_relations(payload, reader)
    finally:
        reader.close()


def test_v29_replayer_metadata_plus_size_conflicting_state_fails_closed(modules):
    replayer = modules["replayer"]
    payload = v14_terminal_key_state_payload()
    reader = replayer.BoundArtifactReader()
    try:
        replayer.recursively_validate_artifact_relations(payload, reader)
        changed = json.loads(json.dumps(payload))
        changed["terminal_key_state"]["legacy_private_key"]["state"] = (
            "CONFLICTING_STATE"
        )
        with pytest.raises(replayer.GateError, match="Conflicting declared relation"):
            replayer.recursively_validate_artifact_relations(changed, reader)
    finally:
        reader.close()


@pytest.mark.parametrize("metadata_first", [False, True])
def test_v29_replayer_full_metadata_reclassification_is_order_independent(
    modules, tmp_path, metadata_first
):
    replayer = modules["replayer"]
    target = tmp_path / "one-time-key.pem"
    target.write_bytes(b"private test material\n")
    target.chmod(0o400)
    full = basic_artifact(target)
    metadata = {
        "path": str(target),
        "size_bytes": target.stat().st_size,
        "mode": "0o400",
        "state": "resident",
    }
    payload = [metadata, full] if metadata_first else [full, metadata]
    reader = replayer.BoundArtifactReader()
    try:
        with pytest.raises(replayer.GateError, match="Conflicting declared relation"):
            replayer.recursively_validate_artifact_relations(payload, reader)
    finally:
        reader.close()


def test_v29_replayer_identical_metadata_plus_size_repeat_is_idempotent(
    modules, tmp_path
):
    replayer = modules["replayer"]
    target = tmp_path / "one-time-key.pem"
    target.write_bytes(b"private test material\n")
    target.chmod(0o400)
    record = {
        "path": str(target),
        "size_bytes": target.stat().st_size,
        "mode": "0o400",
        "state": "resident",
    }
    reader = replayer.BoundArtifactReader()
    try:
        relations = replayer.recursively_validate_artifact_relations(
            [dict(record), dict(record)], reader
        )
        assert relations[str(target)]["size_bytes"] == target.stat().st_size
        assert len(reader._declared_relations) == 1
        assert reader._opened[target.resolve(strict=True)][1] is None
    finally:
        reader.close()


def test_v29_replayer_declared_relation_registry_is_per_reader(modules, tmp_path):
    replayer = modules["replayer"]
    target = tmp_path / "one-time-key.pem"
    target.write_bytes(b"private test material\n")
    target.chmod(0o400)
    metadata = {
        "path": str(target),
        "size_bytes": target.stat().st_size,
        "mode": "0o400",
        "state": "resident",
    }
    full = basic_artifact(target)
    metadata_reader = replayer.BoundArtifactReader()
    content_reader = replayer.BoundArtifactReader()
    try:
        replayer.recursively_validate_artifact_relations(metadata, metadata_reader)
        relations = replayer.recursively_validate_artifact_relations(
            full, content_reader
        )
        assert relations[str(target)] == full
        assert metadata_reader._declared_relations != content_reader._declared_relations
    finally:
        metadata_reader.close()
        content_reader.close()


def test_v29_continuation_consumes_actual_v14_metadata_plus_size_without_bytes(
    modules, monkeypatch
):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    payload = v14_terminal_key_state_payload()
    private_paths = {
        Path(record["path"])
        for role, record in payload["terminal_key_state"].items()
        if role.endswith("private_key")
        if isinstance(record, dict) and "path" in record
    }
    original_open = reader.open

    def guarded_open(path, *args, **kwargs):
        if Path(path) in private_paths:
            raise AssertionError("private key byte reader was called")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(reader, "open", guarded_open)
    try:
        relations = continuation.bind_declared_relations(
            payload, reader, historical_mode_records=set()
        )
        assert {relations[str(path)]["size_bytes"] for path in private_paths} == {119}
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


def test_v29_mandatory_gate_executables_exist_and_match_frozen_contract(modules):
    continuation = modules["continuation"]
    contracts = continuation.MANDATORY_GATE_EXECUTABLE_CONTRACTS
    assert set(contracts) == {
        "primary_python_runtime",
        "qa_python_runtime",
        "node",
        "qa_chromium",
    }
    assert continuation.NODE == Path(
        "/bip7_disk/liaoyoyo2001/.nvm/versions/node/v22.22.1/bin/node"
    )
    for role, contract in contracts.items():
        path = contract["path"]
        digest = hashlib.sha256()
        assert continuation.GATE_INPUT_PATHS[role] == path
        assert path.is_file()
        assert path.resolve(strict=True) == path
        assert path.stat().st_mode & 0o7777 == int(contract["mode"], 8)
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(chunk)
        assert digest.hexdigest() == contract["sha256"]


def test_v29_bound_node_executes_pinned_fd_bytes(modules):
    continuation = modules["continuation"]
    gates = continuation.StreamBoundArtifactSet()
    descriptor, record, _ = gates.open(
        continuation.NODE,
        expected_mode="0o755",
        expected_link_count=1,
    )
    probe = (
        "const fs=require('fs');"
        "const fd=Number(process.argv[1]);"
        "const st=fs.fstatSync(fd,{bigint:false});"
        "console.log(JSON.stringify({execPath:process.execPath,argv0:process.argv0,"
        "fdDevice:st.dev,fdInode:st.ino}));"
    )
    try:
        result = subprocess.run(
            [str(continuation.NODE), "-e", probe, str(descriptor)],
            executable=f"/proc/self/fd/{descriptor}",
            pass_fds=(descriptor,),
            env=continuation.EXPECTED_ENVIRONMENT,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
        observed = json.loads(result.stdout)
        opened = os.fstat(descriptor)
        assert observed == {
            "execPath": str(continuation.NODE),
            "argv0": str(continuation.NODE),
            "fdDevice": opened.st_dev,
            "fdInode": opened.st_ino,
        }
        assert record["sha256"] == continuation.EXPECTED_NODE_SHA256
    finally:
        gates.close()


def test_v29_bound_chromium_executes_pinned_fd_bytes(modules):
    continuation = modules["continuation"]
    gates = continuation.StreamBoundArtifactSet()
    descriptor, record, _ = gates.open(
        continuation.QA_CHROMIUM,
        expected_mode="0o755",
        expected_link_count=1,
    )
    try:
        result = subprocess.run(
            [str(continuation.QA_CHROMIUM), "--version"],
            executable=f"/proc/self/fd/{descriptor}",
            pass_fds=(descriptor,),
            env=continuation.EXPECTED_ENVIRONMENT,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
        assert result.stdout.decode("utf-8").strip().startswith(
            "Google Chrome for Testing "
        )
        opened = os.fstat(descriptor)
        live = os.stat(continuation.QA_CHROMIUM, follow_symlinks=False)
        assert (opened.st_dev, opened.st_ino) == (live.st_dev, live.st_ino)
        assert record["sha256"] == continuation.EXPECTED_QA_CHROMIUM_SHA256
    finally:
        gates.close()


@pytest.mark.parametrize("mutation", ["missing_node", "substituted_chromium"])
def test_v29_missing_or_substituted_mandatory_runtime_fails_closed(
    modules, mutation
):
    continuation = modules["continuation"]
    records = {
        contract["path"]: basic_artifact(contract["path"])
        for contract in continuation.MANDATORY_GATE_EXECUTABLE_CONTRACTS.values()
    }
    if mutation == "missing_node":
        records.pop(continuation.NODE)
    else:
        records[continuation.QA_CHROMIUM] = {
            **records[continuation.QA_CHROMIUM],
            "sha256": "0" * 64,
        }

    class FrozenGateRecords:
        def record_for(self, path):
            return dict(records[path])

    with pytest.raises(
        continuation.ContinuationError,
        match="Mandatory gate executable",
    ):
        continuation.validate_mandatory_gate_executable_records(FrozenGateRecords())


@pytest.mark.parametrize(
    "mutation",
    ["missing_state", "extra_field", "wrong_mode", "boolean_size", "empty_state"],
)
def test_v29_continuation_metadata_plus_size_mutations_fail_closed(
    modules, mutation
):
    continuation = modules["continuation"]
    payload = v14_terminal_key_state_payload()
    record = payload["terminal_key_state"]["legacy_private_key"]
    if mutation == "missing_state":
        record.pop("state")
    elif mutation == "extra_field":
        record["unexpected"] = "forbidden"
    elif mutation == "wrong_mode":
        record["mode"] = "0o444"
    elif mutation == "boolean_size":
        record["size_bytes"] = True
    else:
        record["state"] = ""
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        with pytest.raises(continuation.ContinuationError):
            continuation.bind_declared_relations(
                payload, reader, historical_mode_records=set()
            )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


def test_v29_continuation_metadata_plus_size_conflicting_state_fails_closed(
    modules
):
    continuation = modules["continuation"]
    payload = v14_terminal_key_state_payload()
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        continuation.bind_declared_relations(
            payload, reader, historical_mode_records=set()
        )
        changed = json.loads(json.dumps(payload))
        changed["terminal_key_state"]["legacy_private_key"]["state"] = (
            "CONFLICTING_STATE"
        )
        with pytest.raises(continuation.ContinuationError, match="Conflicting"):
            continuation.bind_declared_relations(
                changed, reader, historical_mode_records=set()
            )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


def collect_records_for_path(value, path):
    records = []
    if isinstance(value, dict):
        if value.get("path") == str(path):
            records.append(value)
        for child in value.values():
            records.extend(collect_records_for_path(child, path))
    elif isinstance(value, list):
        for child in value:
            records.extend(collect_records_for_path(child, path))
    return records


def test_v29_continuation_accepts_actual_failed_v17_replay_metadata_enrichment(
    modules
):
    continuation = modules["continuation"]
    archived_replay = (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v17_signed_authority_c_metadata_enrichment_relation_conflict"
        / "m2v5_runner_only_gate_replay.recovery.v17.json"
    )
    payload = json.loads(archived_replay.read_text(encoding="utf-8"))
    legacy_key = Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
        "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem"
    )
    records = collect_records_for_path(payload, legacy_key)
    assert sorted(sorted(record) for record in records) == sorted(
        [
            ["mode", "path"],
            ["mode", "path", "size_bytes", "state"],
            ["must_remain_unretired", "path", "required_mode"],
        ]
    )
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        relations = continuation.bind_declared_relations(
            records, reader, historical_mode_records=set()
        )
        assert relations[str(legacy_key)] == {
            "path": str(legacy_key),
            "mode": "0o400",
            "size_bytes": 119,
        }
        kind, declared = continuation._DECLARED_RELATION_REGISTRY[reader][
            str(legacy_key)
        ]
        assert kind == "metadata_plus_size"
        assert declared["state"] == "LEGACY_SIGNED_CONTRACT_UNUSED_NOT_RETIRED"
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize("metadata_first", [False, True])
def test_v29_continuation_metadata_enrichment_is_order_independent(
    modules, tmp_path, metadata_first
):
    continuation = modules["continuation"]
    target = tmp_path / "one-time-key.pem"
    target.write_bytes(b"private test material\n")
    target.chmod(0o400)
    metadata = {"path": str(target), "mode": "0o400"}
    enriched = {
        "path": str(target),
        "mode": "0o400",
        "size_bytes": target.stat().st_size,
        "state": "resident",
    }
    payload = [metadata, enriched] if metadata_first else [enriched, metadata]
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        relations = continuation.bind_declared_relations(
            payload, reader, historical_mode_records=set()
        )
        assert relations[str(target)]["size_bytes"] == target.stat().st_size
        kind, declared = continuation._DECLARED_RELATION_REGISTRY[reader][
            str(target)
        ]
        assert kind == "metadata_plus_size"
        assert declared == enriched
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize(
    ("key_path", "state"),
    [
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260722_m2v5_terminal_v2/ed25519_private_one_time_resident.pem",
            "LEGACY_SIGNED_CONTRACT_UNUSED_NOT_RETIRED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v6/ed25519_private_one_time_resident.pem",
            "FAILED_V16_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v7/ed25519_private_one_time_resident.pem",
            "FAILED_V17_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v8/ed25519_private_one_time_resident.pem",
            "FAILED_V18_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v9/ed25519_private_one_time_resident.pem",
            "FAILED_V19_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v10/ed25519_private_one_time_resident.pem",
            "REJECTED_V20_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v11/ed25519_private_one_time_resident.pem",
            "FAILED_V21_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v12/ed25519_private_one_time_resident.pem",
            "REJECTED_V22_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260723_m2v5_terminal_v13/ed25519_private_one_time_resident.pem",
            "REJECTED_V23_UNUSED_NOT_RETIRED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260723_m2v5_terminal_v15_unused_after_signed_v25_c_failure_01/"
            "ed25519_private_one_time_resident.pem",
            "FAILED_V25_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v16_unused_after_signed_v26_c_failure_01/"
            "ed25519_private_one_time_resident.pem",
            "FAILED_V26_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "archive/20260724_m2v5_terminal_v19_unused_after_signed_v29_r_failure_01/"
            "ed25519_private_one_time_resident.pem",
            "FAILED_V29_UNUSED_ARCHIVED_MUST_NOT_BE_REUSED",
        ),
        (
            "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
            "20260724_m2v5_terminal_v21/ed25519_private_one_time_resident.pem",
            "RECOVERY_V30_READY_FOR_SINGLE_SIGNATURE",
        ),
    ],
)
def test_v29_continuation_actual_terminal_key_metadata_enrichment(
    modules, key_path, state
):
    continuation = modules["continuation"]
    path = Path(key_path)
    payload = [
        {"path": str(path), "mode": "0o400"},
        {
            "path": str(path),
            "mode": "0o400",
            "size_bytes": 119,
            "state": state,
        },
    ]
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        relations = continuation.bind_declared_relations(
            payload, reader, historical_mode_records=set()
        )
        assert relations[str(path)]["size_bytes"] == 119
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize(
    ("field", "value"),
    [("size_bytes", 120), ("state", "CONFLICTING_STATE")],
)
def test_v29_continuation_enriched_metadata_conflict_still_fails_closed(
    modules, tmp_path, field, value
):
    continuation = modules["continuation"]
    target = tmp_path / "one-time-key.pem"
    target.write_bytes(b"private test material\n")
    target.chmod(0o400)
    metadata = {"path": str(target), "mode": "0o400"}
    enriched = {
        "path": str(target),
        "mode": "0o400",
        "size_bytes": target.stat().st_size,
        "state": "resident",
    }
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        continuation.bind_declared_relations(
            [metadata, enriched], reader, historical_mode_records=set()
        )
        conflicting = dict(enriched)
        conflicting[field] = value
        with pytest.raises(continuation.ContinuationError, match="Conflicting"):
            continuation.bind_declared_relations(
                conflicting, reader, historical_mode_records=set()
            )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


def test_v29_continuation_declared_relation_registry_is_per_reader(
    modules, tmp_path
):
    continuation = modules["continuation"]
    target = tmp_path / "one-time-key.pem"
    target.write_bytes(b"private test material\n")
    target.chmod(0o400)
    metadata = {
        "path": str(target),
        "size_bytes": target.stat().st_size,
        "mode": "0o400",
        "state": "resident",
    }
    full = basic_artifact(target)
    metadata_reader = modules["release_authority"].BoundArtifactReader()
    content_reader = modules["release_authority"].BoundArtifactReader()
    try:
        continuation.bind_declared_relations(
            metadata, metadata_reader, historical_mode_records=set()
        )
        relations = continuation.bind_declared_relations(
            full, content_reader, historical_mode_records=set()
        )
        assert relations[str(target)] == full
        assert len(continuation._DECLARED_RELATION_REGISTRY) == 2
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        metadata_reader.close()
        content_reader.close()


def executable_alias_record(alias: Path, target: Path) -> dict[str, object]:
    data = target.read_bytes()
    return {
        "path": str(alias),
        "resolved_path": str(target),
        "symlink_target": target.name,
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(target.stat().st_mode & 0o7777),
    }


def close_executable_alias_leases(continuation) -> None:
    for binding in continuation._EXECUTABLE_ALIAS_LEASES.values():
        os.close(binding.alias_fd)
        binding.mutation_watch.close()
    continuation._EXECUTABLE_ALIAS_LEASES.clear()
    continuation._DECLARED_RELATION_REGISTRY.clear()


def make_executable_alias(tmp_path: Path) -> tuple[Path, Path, dict[str, object]]:
    target = tmp_path / "python3.11"
    target.write_bytes(b"reviewed python runtime\n")
    target.chmod(0o755)
    alias = tmp_path / "python"
    alias.symlink_to(target.name)
    return alias, target, executable_alias_record(alias, target)


def test_v29_declared_relation_binds_executable_alias_and_target(modules, tmp_path):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    alias, _, record = make_executable_alias(tmp_path)
    try:
        relations = continuation.bind_declared_relations(
            {"runtime": {"python": record}},
            reader,
            historical_mode_records=set(),
        )
        assert relations[str(alias)] == record
        continuation.require_bootstrap_stable()
    finally:
        close_executable_alias_leases(continuation)
        reader.close()


def test_v29_replayer_binds_expected_executable_alias_target(modules, tmp_path):
    replayer = modules["replayer"]
    reader = replayer.BoundArtifactReader()
    alias, target, record = make_executable_alias(tmp_path)
    try:
        _, data, observed = reader.open_executable_alias(
            alias,
            expected_link_text=target.name,
            expected_resolved_path=target,
            expected_target_size_bytes=record["size_bytes"],
            expected_target_sha256=record["sha256"],
            expected_target_mode=record["mode"],
        )
        assert data == target.read_bytes()
        assert observed == record
        reader.require_paths_still_bound()
    finally:
        reader.close()


@pytest.mark.parametrize("mutation", ["size_bytes", "sha256", "mode"])
def test_v29_replayer_rejects_executable_target_drift_before_use(
    modules, tmp_path, mutation
):
    replayer = modules["replayer"]
    reader = replayer.BoundArtifactReader()
    alias, target, record = make_executable_alias(tmp_path)
    expected = dict(record)
    if mutation == "size_bytes":
        expected[mutation] += 1
    elif mutation == "sha256":
        expected[mutation] = "0" * 64
    else:
        expected[mutation] = "0o700"
    try:
        with pytest.raises(replayer.GateError, match="target identity drift"):
            reader.open_executable_alias(
                alias,
                expected_link_text=target.name,
                expected_resolved_path=target,
                expected_target_size_bytes=expected["size_bytes"],
                expected_target_sha256=expected["sha256"],
                expected_target_mode=expected["mode"],
            )
    finally:
        reader.close()


def test_v29_replayer_transient_target_mutation_is_permanent_failure(
    modules, tmp_path
):
    replayer = modules["replayer"]
    reader = replayer.BoundArtifactReader()
    alias, target, record = make_executable_alias(tmp_path)
    original = target.read_bytes()
    try:
        reader.open_executable_alias(
            alias,
            expected_link_text=target.name,
            expected_resolved_path=target,
            expected_target_size_bytes=record["size_bytes"],
            expected_target_sha256=record["sha256"],
            expected_target_mode=record["mode"],
        )
        target.write_bytes(b"transient runtime bytes\n")
        target.write_bytes(original)
        with pytest.raises(replayer.GateError, match="mutation event"):
            reader.require_paths_still_bound()
        with pytest.raises(replayer.GateError, match="mutation event"):
            reader.require_paths_still_bound()
    finally:
        reader.close()


def test_v29_terminal_prelink_failure_leaves_target_absent_and_restores_signals(
    modules, tmp_path
):
    replayer = modules["replayer"]
    target = tmp_path / "terminal-receipt.json"
    signal_mask_before = signal.pthread_sigmask(signal.SIG_BLOCK, set())

    def fail_terminal_gate() -> None:
        raise replayer.GateError("injected terminal gate failure")

    with pytest.raises(replayer.GateError, match="injected terminal gate failure"):
        replayer.write_exclusive_readonly(
            target,
            b'{"pass":true}\n',
            terminal_commit=fail_terminal_gate,
        )
    signal_mask_after = signal.pthread_sigmask(signal.SIG_BLOCK, set())
    assert signal_mask_after == signal_mask_before
    assert not target.exists()


def test_v29_terminal_success_links_then_exits_zero_without_returning(
    modules, tmp_path
):
    replayer = modules["replayer"]
    source = Path(replayer.__file__).resolve()
    target = tmp_path / "terminal-receipt.json"
    payload = b'{"pass":true}\n'
    program = f"""
import importlib.util
from pathlib import Path
spec = importlib.util.spec_from_file_location("terminal_replayer", {str(source)!r})
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.write_exclusive_readonly(
    Path({str(target)!r}),
    {payload!r},
    terminal_commit=lambda: None,
)
raise SystemExit(91)
"""
    completed = subprocess.run(
        [str(replayer.PYTHON_TARGET), "-I", "-B", "-c", program],
        check=False,
        capture_output=True,
    )
    assert completed.returncode == 0
    assert completed.stdout == b""
    assert completed.stderr == b""
    assert target.read_bytes() == payload
    assert target.stat().st_mode & 0o7777 == 0o444


@pytest.mark.parametrize(
    "wait_status",
    [7 << 8, int(signal.SIGTERM), True],
)
def test_v29_replayer_worker_nonzero_signal_or_bool_wait_fails_closed(
    modules, wait_status
):
    replayer = modules["replayer"]
    with pytest.raises(replayer.GateError, match="normal exit zero"):
        replayer.require_worker_normal_exit(123, 123, wait_status)


def valid_replay_success_witness(continuation, source, log, receipt):
    return {
        "schema_name": continuation.RUNNER_GATE_SUCCESS_WITNESS_SCHEMA,
        "schema_version": "1.0.0",
        "created_at_utc": "2026-07-23T12:00:01+00:00",
        "task_type": "B_comprehensive_validation",
        "status": "WORKER_WAITPID_NORMAL_EXIT_ZERO_ATTESTED",
        "scope": "completion_runner_physical_lines_1_358_only",
        "canonical_command": continuation.replay_command(),
        "clean_environment": continuation.EXPECTED_ENVIRONMENT,
        "supervisor_source": source,
        "worker_wait": {
            "pid": 123,
            "raw_wait_status": 0,
            "normal_exit": True,
            "exit_code": 0,
            "term_signal": None,
        },
        "worker_outputs": {"log": log, "receipt": receipt},
        "checks": continuation.REPLAY_SUCCESS_WITNESS_CHECKS,
        "pass": True,
        "pass_semantics": continuation.REPLAY_SUCCESS_WITNESS_PASS_SEMANTICS,
    }


def test_v29_replay_success_witness_accepts_exact_normal_exit_zero(modules):
    continuation = modules["continuation"]
    source = {
        "path": "/tmp/replayer.py",
        "size_bytes": 1,
        "sha256": "1" * 64,
        "mode": "0o444",
    }
    log = {
        "path": "/tmp/replay.log",
        "size_bytes": 2,
        "sha256": "2" * 64,
        "mode": "0o444",
    }
    receipt = {
        "path": "/tmp/replay.json",
        "size_bytes": 3,
        "sha256": "3" * 64,
        "mode": "0o444",
    }
    witness = valid_replay_success_witness(continuation, source, log, receipt)
    observed = continuation.validate_replay_success_witness(
        witness,
        chain={"replay_record": source},
        receipt_record=receipt,
        log_record=log,
        replay_time=datetime(2026, 7, 23, 12, 0, tzinfo=timezone.utc),
    )
    assert observed["worker_wait"]["exit_code"] == 0


@pytest.mark.parametrize(
    "mutation", ["missing_key", "extra_key", "bool_exit", "receipt_mismatch"]
)
def test_v29_replay_success_witness_schema_and_identity_fail_closed(
    modules, mutation
):
    continuation = modules["continuation"]
    source = {
        "path": "/tmp/replayer.py",
        "size_bytes": 1,
        "sha256": "1" * 64,
        "mode": "0o444",
    }
    log = {
        "path": "/tmp/replay.log",
        "size_bytes": 2,
        "sha256": "2" * 64,
        "mode": "0o444",
    }
    receipt = {
        "path": "/tmp/replay.json",
        "size_bytes": 3,
        "sha256": "3" * 64,
        "mode": "0o444",
    }
    witness = valid_replay_success_witness(continuation, source, log, receipt)
    if mutation == "missing_key":
        witness.pop("checks")
    elif mutation == "extra_key":
        witness["unexpected"] = True
    elif mutation == "bool_exit":
        witness["worker_wait"]["exit_code"] = False
    else:
        witness["worker_outputs"]["receipt"] = dict(receipt, sha256="4" * 64)
    with pytest.raises(continuation.ContinuationError):
        continuation.validate_replay_success_witness(
            witness,
            chain={"replay_record": source},
            receipt_record=receipt,
            log_record=log,
            replay_time=datetime(2026, 7, 23, 12, 0, tzinfo=timezone.utc),
        )


def test_v29_actual_promotion_payloads_accept_exact_legacy_eight_key_stat(modules):
    continuation = modules["continuation"]
    authorization = json.loads(continuation.AUTHORIZATION.read_text(encoding="utf-8"))
    completion = json.loads(continuation.COMPLETION.read_text(encoding="utf-8"))
    historical_modes = {
        (
            (
                "authorization_signature_contract",
                "private_key_lifecycle",
                "pre_signature",
            ),
            str(continuation.AUTHORIZATION_PRIVATE_KEY),
            "0o400",
        ),
        (
            (
                "completion_signature_contract",
                "private_key_lifecycle",
                "pre_signature",
            ),
            str(continuation.COMPLETION_PRIVATE_KEY),
            "0o400",
        ),
        (
            ("continuation_gate", "terminal_signature_contract", "private_key"),
            str(continuation.CONTINUATION_PRIVATE_KEY),
            "0o400",
        ),
        (
            ("signature_contract", "private_key_lifecycle", "pre_signature"),
            str(continuation.COMPLETION_PRIVATE_KEY),
            "0o400",
        ),
    }
    focal_during = (
        authorization.get("evidence", {})
        .get("focal_source_identity_transition", {})
        .get("during_execution")
    )
    if isinstance(focal_during, dict):
        historical_modes.add(
            (
                (
                    "evidence",
                    "focal_source_identity_transition",
                    "during_execution",
                ),
                str(focal_during.get("path")),
                str(focal_during.get("mode")),
            )
        )
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        relations = continuation.bind_declared_relations(
            authorization,
            reader,
            historical_mode_records=historical_modes,
        )
        continuation.bind_declared_relations(
            completion,
            reader,
            historical_mode_records=historical_modes,
            relations=relations,
        )
        legacy_path = str(
            REPO_ROOT
            / "research/20260715_single_fp_alt_multicluster_subclone_validation/"
            "scripts/focal_alt_cluster_lib.py"
        )
        assert relations[legacy_path] == basic_artifact(Path(legacy_path))
    finally:
        close_executable_alias_leases(continuation)
        reader.close()

    relocated_reader = modules["release_authority"].BoundArtifactReader()
    try:
        with pytest.raises(continuation.ContinuationError):
            continuation.bind_declared_relations(
                {"relocated_historical_record": focal_during},
                relocated_reader,
                historical_mode_records=historical_modes,
            )
    finally:
        close_executable_alias_leases(continuation)
        relocated_reader.close()


@pytest.mark.parametrize("order", ["basic_first", "legacy_first"])
def test_v29_basic_and_legacy_stat_are_compatible_content_enrichment(
    modules, tmp_path, order
):
    continuation = modules["continuation"]
    target = tmp_path / "artifact.json"
    target.write_text("{}\n", encoding="utf-8")
    target.chmod(0o444)
    basic = basic_artifact(target)
    legacy = legacy_stat_artifact(target)
    sequence = [basic, legacy] if order == "basic_first" else [legacy, basic]
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        relations = {}
        for record in sequence:
            continuation.bind_declared_relations(
                {"artifact": record},
                reader,
                historical_mode_records=set(),
                relations=relations,
            )
        assert relations[str(target)] == basic
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize(
    "mutation",
    ["missing_ctime", "extra_field", "boolean_inode", "shared_sha_conflict"],
)
def test_v29_legacy_stat_near_misses_fail_closed(modules, tmp_path, mutation):
    continuation = modules["continuation"]
    target = tmp_path / "artifact.json"
    target.write_text("{}\n", encoding="utf-8")
    target.chmod(0o444)
    basic = basic_artifact(target)
    record = legacy_stat_artifact(target)
    if mutation == "missing_ctime":
        record.pop("ctime_ns")
    elif mutation == "extra_field":
        record["unexpected"] = True
    elif mutation == "boolean_inode":
        record["inode"] = True
    else:
        record["sha256"] = "0" * 64
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        if mutation == "shared_sha_conflict":
            continuation.bind_declared_relations(
                {"artifact": basic}, reader, historical_mode_records=set()
            )
            with pytest.raises(continuation.ContinuationError, match="Conflicting"):
                continuation.bind_declared_relations(
                    {"artifact": record}, reader, historical_mode_records=set()
                )
        else:
            with pytest.raises(continuation.ContinuationError):
                continuation.bind_declared_relations(
                    {"artifact": record}, reader, historical_mode_records=set()
                )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize("schema", ["basic", "stat", "legacy_stat"])
def test_v29_declared_artifact_schema_rejects_extra_fields(modules, tmp_path, schema):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    target = tmp_path / "artifact.json"
    target.write_text("{}\n", encoding="utf-8")
    target.chmod(0o444)
    record = basic_artifact(target)
    if schema in {"stat", "legacy_stat"}:
        observed = target.stat()
        record.update(
            {
                "device": observed.st_dev,
                "inode": observed.st_ino,
                "mtime_ns": observed.st_mtime_ns,
                "ctime_ns": observed.st_ctime_ns,
            }
        )
        if schema == "stat":
            record["link_count"] = observed.st_nlink
    record["unexpected"] = True
    try:
        with pytest.raises(
            continuation.ContinuationError,
            match="Ambiguous identity-like relation schema",
        ):
            continuation.bind_declared_relations(
                {"artifact": record}, reader, historical_mode_records=set()
            )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize(
    "record",
    [
        {"path": "/tmp/incomplete", "size_bytes": 1},
        {"path": "/tmp/incomplete", "sha256": "0" * 64},
        {"path": "/tmp/incomplete", "resolved_path": "/tmp/target"},
        {"path": "/tmp/incomplete", "mode": "0o444", "unexpected": True},
    ],
)
def test_v29_partial_identity_like_records_fail_closed(modules, record):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    try:
        with pytest.raises(
            continuation.ContinuationError,
            match="Ambiguous identity-like relation schema",
        ):
            continuation.bind_declared_relations(
                {"artifact": record}, reader, historical_mode_records=set()
            )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize("kind", ["basic_size", "stat_device"])
def test_v29_identity_integer_fields_reject_bool(modules, tmp_path, kind):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    target = tmp_path / "artifact.json"
    target.write_text("{}\n", encoding="utf-8")
    target.chmod(0o444)
    record = basic_artifact(target)
    if kind == "basic_size":
        record["size_bytes"] = True
    else:
        observed = target.stat()
        record.update(
            {
                "device": True,
                "inode": observed.st_ino,
                "link_count": observed.st_nlink,
                "mtime_ns": observed.st_mtime_ns,
                "ctime_ns": observed.st_ctime_ns,
            }
        )
    try:
        with pytest.raises(continuation.ContinuationError, match="type drift"):
            continuation.bind_declared_relations(
                {"artifact": record}, reader, historical_mode_records=set()
            )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


def test_v29_same_path_schema_downgrade_fails_closed(modules, tmp_path):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    target = tmp_path / "artifact.json"
    target.write_text("{}\n", encoding="utf-8")
    target.chmod(0o444)
    basic = basic_artifact(target)
    metadata = {"path": str(target), "mode": basic["mode"]}
    try:
        continuation.bind_declared_relations(
            {"artifact": basic}, reader, historical_mode_records=set()
        )
        with pytest.raises(continuation.ContinuationError, match="Conflicting"):
            continuation.bind_declared_relations(
                {"artifact": metadata}, reader, historical_mode_records=set()
            )
    finally:
        continuation._DECLARED_RELATION_REGISTRY.clear()
        reader.close()


@pytest.mark.parametrize(
    "location", ["runtime", "promotion_fd", "completion_contract"]
)
def test_v29_replay_python_alias_cannot_downgrade_to_basic(
    modules, monkeypatch, tmp_path, location
):
    continuation = modules["continuation"]
    alias, target, alias_record = make_executable_alias(tmp_path)
    target_record = basic_artifact(target)
    monkeypatch.setattr(continuation, "PYTHON", alias)
    monkeypatch.setattr(continuation, "PRIMARY_PYTHON_RUNTIME", target)
    records = {
        "runtime": dict(alias_record),
        "promotion_fd": dict(alias_record),
        "completion_contract": dict(alias_record),
    }
    downgraded = {
        key: value
        for key, value in target_record.items()
    }
    records[location] = downgraded
    assert alias != target
    with pytest.raises(continuation.ContinuationError, match="Python alias drift"):
        continuation.require_replay_python_alias_consistency(
            target_record,
            records["runtime"],
            records["promotion_fd"],
            records["completion_contract"],
        )


def test_v29_duplicate_executable_alias_relation_is_idempotent(modules, tmp_path):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    alias, _, record = make_executable_alias(tmp_path)
    try:
        first = continuation.bind_declared_relations(
            {"runtime": {"python": record}}, reader, historical_mode_records=set()
        )
        second = continuation.bind_declared_relations(
            {"runtime": {"python": dict(record)}},
            reader,
            historical_mode_records=set(),
        )
        assert first == second == {str(alias): record}
        assert len(continuation._EXECUTABLE_ALIAS_LEASES) == 1
    finally:
        close_executable_alias_leases(continuation)
        reader.close()


def test_v29_conflicting_executable_alias_relation_fails_closed(modules, tmp_path):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    _, _, record = make_executable_alias(tmp_path)
    try:
        continuation.bind_declared_relations(
            {"runtime": {"python": record}}, reader, historical_mode_records=set()
        )
        conflicting = dict(record)
        conflicting["sha256"] = "0" * 64
        with pytest.raises(continuation.ContinuationError, match="Conflicting"):
            continuation.bind_declared_relations(
                {"runtime": {"python": conflicting}},
                reader,
                historical_mode_records=set(),
            )
    finally:
        close_executable_alias_leases(continuation)
        reader.close()


def test_v29_stream_bound_alias_accepts_exact_relative_canonical_target(
    modules, tmp_path
):
    continuation = modules["continuation"]
    target = tmp_path / "all_ssnv_tumor_ref_control_summary.json"
    alias = tmp_path / "summary.json"
    target.write_bytes(b'{"pass":true}\n')
    alias.symlink_to(target.name)
    artifacts = continuation.StreamBoundArtifactSet()
    try:
        artifacts.open(target, expected_link_count=1)
        record = artifacts.open_alias(
            alias, target, expected_link_text=target.name
        )
        assert artifacts.descriptor_count == 2
        assert record == artifacts.alias_record_for(alias)
        assert record["path"] == str(alias)
        assert record["link_text"] == target.name
        assert record["resolved_target"] == str(target)
        assert continuation.GATE_INPUT_PATHS["tumor_ref_summary"] == (
            continuation.TUMOR_REF_SUMMARY
        )
        assert continuation.TUMOR_REF_SUMMARY_ALIAS != continuation.TUMOR_REF_SUMMARY
        source = continuation.CONTINUATION_RUNNER.read_text(encoding="utf-8")
        assert source.count(
            "all_event_directories={PORTABLE_ASSET_ROOT, TUMOR_REF}"
        ) == 3
        assert source.count(
            "str(path) for path in (PORTABLE_ASSET_ROOT, TUMOR_REF)"
        ) == 2
        artifacts.require_paths_still_bound()
    finally:
        artifacts.close()


def test_v29_stream_bound_alias_requires_target_bound_first(modules, tmp_path):
    continuation = modules["continuation"]
    target = tmp_path / "target.json"
    alias = tmp_path / "summary.json"
    target.write_bytes(b"{}\n")
    alias.symlink_to(target.name)
    artifacts = continuation.StreamBoundArtifactSet()
    try:
        with pytest.raises(continuation.ContinuationError, match="bound first"):
            artifacts.open_alias(alias, target, expected_link_text=target.name)
    finally:
        artifacts.close()


def test_v29_stream_bound_alias_rejects_wrong_target(modules, tmp_path):
    continuation = modules["continuation"]
    target = tmp_path / "target.json"
    decoy = tmp_path / "decoy.json"
    alias = tmp_path / "summary.json"
    target.write_bytes(b'{"source":"target"}\n')
    decoy.write_bytes(b'{"source":"decoy"}\n')
    alias.symlink_to(decoy.name)
    artifacts = continuation.StreamBoundArtifactSet()
    try:
        artifacts.open(target, expected_link_count=1)
        with pytest.raises(continuation.ContinuationError, match="target contract"):
            artifacts.open_alias(alias, target, expected_link_text=target.name)
    finally:
        artifacts.close()


def test_v29_stream_bound_alias_rejects_absolute_link_text(modules, tmp_path):
    continuation = modules["continuation"]
    target = tmp_path / "target.json"
    alias = tmp_path / "summary.json"
    target.write_bytes(b"{}\n")
    alias.symlink_to(target)
    artifacts = continuation.StreamBoundArtifactSet()
    try:
        artifacts.open(target, expected_link_count=1)
        with pytest.raises(continuation.ContinuationError, match="identity drift"):
            artifacts.open_alias(alias, target, expected_link_text=target.name)
    finally:
        artifacts.close()


def test_v29_stream_bound_alias_detects_alias_replacement(modules, tmp_path):
    continuation = modules["continuation"]
    target = tmp_path / "target.json"
    alias = tmp_path / "summary.json"
    target.write_bytes(b"{}\n")
    alias.symlink_to(target.name)
    artifacts = continuation.StreamBoundArtifactSet()
    try:
        artifacts.open(target, expected_link_count=1)
        artifacts.open_alias(alias, target, expected_link_text=target.name)
        alias.unlink()
        alias.symlink_to(target.name)
        with pytest.raises(continuation.ContinuationError, match="alias drifted"):
            artifacts.require_paths_still_bound()
    finally:
        artifacts.close()


def test_v29_stream_bound_alias_detects_canonical_target_replacement(
    modules, tmp_path
):
    continuation = modules["continuation"]
    target = tmp_path / "target.json"
    replacement = tmp_path / "replacement.json"
    alias = tmp_path / "summary.json"
    target.write_bytes(b'{"version":1}\n')
    replacement.write_bytes(b'{"version":2}\n')
    alias.symlink_to(target.name)
    artifacts = continuation.StreamBoundArtifactSet()
    try:
        artifacts.open(target, expected_link_count=1)
        artifacts.open_alias(alias, target, expected_link_text=target.name)
        replacement.replace(target)
        with pytest.raises(continuation.ContinuationError, match="input drifted"):
            artifacts.require_paths_still_bound()
    finally:
        artifacts.close()


def test_v29_tumor_ref_directory_watch_detects_transient_alias_restore(
    modules, tmp_path
):
    continuation = modules["continuation"]
    target = tmp_path / "all_ssnv_tumor_ref_control_summary.json"
    alias = tmp_path / "summary.json"
    target.write_bytes(b'{"pass":true}\n')
    alias.symlink_to(target.name)
    sentinel = continuation.PathMutationSentinel(
        {target}, all_event_directories={tmp_path}
    )
    try:
        sentinel.assert_clean()
        alias.unlink()
        alias.symlink_to(target.name)
        assert os.readlink(alias) == target.name
        assert alias.resolve(strict=True) == target
        with pytest.raises(
            continuation.ContinuationError,
            match="Protected path mutation was observed",
        ):
            sentinel.assert_clean()
    finally:
        os.close(sentinel._fd)


def test_v29_mode000_file_uses_named_parent_watch_and_detects_chmod(
    modules, tmp_path
):
    continuation = modules["continuation"]
    retired_key = tmp_path / "retired_private.pem"
    retired_key.write_bytes(b"retired-key-metadata-only\n")
    retired_key.chmod(0o000)
    sentinel = continuation.PathMutationSentinel({retired_key})
    try:
        expected_hash = hashlib.sha256(
            (str(retired_key) + "\n").encode("utf-8")
        ).hexdigest()
        assert sentinel.protected_path_count == 1
        assert sentinel.direct_file_watch_count == 0
        assert sentinel.parent_only_mode000_path_count == 1
        assert sentinel.parent_only_mode000_path_set_sha256 == expected_hash
        assert sentinel.setup_identity_recheck_count == 1
        assert sentinel.setup_identity_recheck_path_set_sha256 == expected_hash
        sentinel.assert_clean()
        retired_key.chmod(0o600)
        with pytest.raises(
            continuation.ContinuationError,
            match="Protected path mutation was observed",
        ):
            sentinel.assert_clean()
    finally:
        retired_key.chmod(0o600)
        os.close(sentinel._fd)


def test_v29_mode000_named_parent_watch_detects_replacement(modules, tmp_path):
    continuation = modules["continuation"]
    retired_key = tmp_path / "retired_private.pem"
    replacement = tmp_path / "replacement.pem"
    retired_key.write_bytes(b"original\n")
    replacement.write_bytes(b"replacement\n")
    retired_key.chmod(0o000)
    replacement.chmod(0o000)
    sentinel = continuation.PathMutationSentinel({retired_key})
    try:
        sentinel.assert_clean()
        replacement.replace(retired_key)
        with pytest.raises(
            continuation.ContinuationError,
            match="Protected path mutation was observed",
        ):
            sentinel.assert_clean()
    finally:
        retired_key.chmod(0o600)
        os.close(sentinel._fd)


@pytest.mark.parametrize("mutation", ["chmod_restore", "replacement_restore"])
def test_v29_mode000_setup_window_transient_mutation_fails_closed(
    modules, monkeypatch, tmp_path, mutation
):
    continuation = modules["continuation"]
    retired_key = tmp_path / "retired_private.pem"
    retired_key.write_bytes(b"original\n")
    retired_key.chmod(0o000)
    real_same_stat = continuation.same_stat
    injected = False

    def mutate_then_compare(left, right):
        nonlocal injected
        if not injected:
            injected = True
            if mutation == "chmod_restore":
                retired_key.chmod(0o400)
                retired_key.chmod(0o000)
            else:
                replacement = tmp_path / "replacement.pem"
                held_original = tmp_path / "held-original.pem"
                displaced_replacement = tmp_path / "displaced-replacement.pem"
                replacement.write_bytes(b"replacement\n")
                replacement.chmod(0o000)
                retired_key.rename(held_original)
                replacement.rename(retired_key)
                retired_key.rename(displaced_replacement)
                held_original.rename(retired_key)
        return real_same_stat(left, right)

    monkeypatch.setattr(continuation, "same_stat", mutate_then_compare)
    try:
        with pytest.raises(
            continuation.ContinuationError,
            match="Protected path mutation was observed",
        ):
            continuation.PathMutationSentinel({retired_key})
        assert injected is True
    finally:
        retired_key.chmod(0o600)


def test_v29_parent_only_fallback_rejects_other_modes_and_errors(modules):
    sentinel = modules["continuation"].PathMutationSentinel
    continuation = modules["continuation"]
    assert sentinel._allows_parent_only_mode000_fallback(
        error_number=continuation.errno.EACCES,
        file_mode=0o000,
    )
    assert not sentinel._allows_parent_only_mode000_fallback(
        error_number=continuation.errno.EACCES,
        file_mode=0o400,
    )
    assert not sentinel._allows_parent_only_mode000_fallback(
        error_number=continuation.errno.EPERM,
        file_mode=0o000,
    )


def test_v29_terminal_key_rotation_preserves_legacy_and_pins_fresh_key(modules):
    validator = modules["validator"]
    continuation = modules["continuation"]
    contract, state = validator._validate_terminal_key_rotation(
        expected_recovery_private_mode="0o400"
    )
    assert state["legacy_private_key"]["mode"] == "0o400"
    assert state["failed_v16_private_key"]["mode"] == "0o400"
    assert state["failed_v17_private_key"]["mode"] == "0o400"
    assert state["failed_v18_private_key"]["mode"] == "0o400"
    assert state["failed_v19_private_key"]["mode"] == "0o400"
    assert state["rejected_v20_private_key"]["mode"] == "0o400"
    assert state["failed_v21_private_key"]["mode"] == "0o400"
    assert state["rejected_v22_private_key"]["mode"] == "0o400"
    assert state["rejected_v23_private_key"]["mode"] == "0o400"
    assert state["rejected_v24_private_key"]["mode"] == "0o400"
    assert state["failed_v25_private_key"]["mode"] == "0o400"
    assert state["failed_v26_private_key"]["mode"] == "0o400"
    assert state["rejected_v27_private_key"]["mode"] == "0o400"
    assert state["failed_v28_private_key"]["mode"] == "0o400"
    assert state["recovery_private_key"]["mode"] == "0o400"
    assert state["legacy_public_key"]["sha256"] == (
        continuation.EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["recovery_public_key"]["sha256"] == (
        continuation.EXPECTED_RECOVERY_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v16_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V16_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v17_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V17_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v18_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V18_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v19_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V19_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["rejected_v20_public_key"]["sha256"] == (
        continuation.EXPECTED_REJECTED_V20_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v21_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V21_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["rejected_v22_public_key"]["sha256"] == (
        continuation.EXPECTED_REJECTED_V22_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["rejected_v23_public_key"]["sha256"] == (
        continuation.EXPECTED_REJECTED_V23_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["rejected_v24_public_key"]["sha256"] == (
        continuation.EXPECTED_REJECTED_V24_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v25_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V25_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v26_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V26_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["rejected_v27_public_key"]["sha256"] == (
        continuation.EXPECTED_REJECTED_V27_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert state["failed_v28_public_key"]["sha256"] == (
        validator.FAILED_V28_KEY_ARCHIVES["terminal_v18"]["public_sha256"]
    )
    assert state["failed_v29_private_key"]["mode"] == "0o400"
    assert state["failed_v29_public_key"]["sha256"] == (
        continuation.EXPECTED_FAILED_V29_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert len(
        {
            state["legacy_public_key"]["sha256"],
            state["failed_v16_public_key"]["sha256"],
            state["failed_v17_public_key"]["sha256"],
            state["failed_v18_public_key"]["sha256"],
            state["failed_v19_public_key"]["sha256"],
            state["rejected_v20_public_key"]["sha256"],
            state["failed_v21_public_key"]["sha256"],
            state["rejected_v22_public_key"]["sha256"],
            state["rejected_v23_public_key"]["sha256"],
            state["rejected_v24_public_key"]["sha256"],
            state["failed_v25_public_key"]["sha256"],
            state["failed_v26_public_key"]["sha256"],
            state["rejected_v27_public_key"]["sha256"],
            state["failed_v28_public_key"]["sha256"],
            state["failed_v29_public_key"]["sha256"],
            state["recovery_public_key"]["sha256"],
        }
    ) == 16
    assert contract["legacy_signed_contract"]["status"] == "PRESERVED_NOT_EXECUTED"
    assert contract["failed_v16_contract"]["status"] == (
        "FAILED_BEFORE_C_CHILD_UNUSED_KEY_QUARANTINED"
    )
    assert contract["failed_v17_contract"]["status"] == (
        "FAILED_AFTER_FRESH_V_BEFORE_DOWNSTREAM_UNUSED_KEY_QUARANTINED"
    )
    assert contract["failed_v18_contract"]["status"] == (
        "FAILED_AFTER_FRESH_V_AND_R_BEFORE_DOWNSTREAM_UNUSED_KEY_QUARANTINED"
    )
    assert contract["failed_v19_contract"]["status"] == (
        "FAILED_AFTER_FRESH_V_AND_R_BEFORE_DOWNSTREAM_MODE000_INOTIFY_"
        "EACCES_UNUSED_KEY_QUARANTINED"
    )
    assert contract["rejected_v20_contract"]["status"] == (
        "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED"
    )
    assert contract["failed_v21_contract"]["status"] == (
        "FAILED_AFTER_INTERMEDIATE_DOWNSTREAM_OUTPUTS_BEFORE_FINAL_DATASET_"
        "UNUSED_KEY_QUARANTINED"
    )
    assert contract["rejected_v22_contract"]["status"] == (
        "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED"
    )
    assert contract["rejected_v23_contract"]["status"] == (
        "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED"
    )
    assert contract["rejected_v24_contract"]["status"] == (
        "REJECTED_PRE_AUTHORITY_UNUSED_KEY_QUARANTINED"
    )
    assert contract["failed_v25_contract"]["status"] == (
        "SIGNED_AUTHORITY_C_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED"
    )
    assert contract["failed_v26_contract"]["status"] == (
        "SIGNED_AUTHORITY_C_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED"
    )
    assert contract["rejected_v27_contract"]["status"] == (
        "REJECTED_PRE_AUTHORITY_UNUSED_TERMINAL_KEY_ARCHIVED"
    )
    assert contract["failed_v28_contract"]["status"] == (
        "SIGNED_AUTHORITY_C_REPORT_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED"
    )
    assert contract["failed_v29_contract"]["status"] == (
        "SIGNED_AUTHORITY_R_FAILED_UNUSED_TERMINAL_KEY_ARCHIVED"
    )
    assert contract["recovery_v30_contract"]["status"] == (
        "AUTHORIZED_BY_SIGNED_V30_RECOVERY_AUTHORITY_ONLY"
    )


def test_v29_v_and_r_validate_legacy_key_while_c_separates_recovery_key(modules):
    validator = modules["validator"]
    verifier = modules["verifier"]
    replayer = modules["replayer"]
    continuation = modules["continuation"]
    builder = modules["builder"]
    validator_sha256 = hashlib.sha256(validator.VALIDATOR.read_bytes()).hexdigest()
    assert validator_sha256 != (
        "54716f6e3da56f9b47d60972ec2ed42469412a62c789026abd698890fe94171a"
    )
    assert {
        verifier.EXPECTED_RECOVERY_AUTHORITY_VALIDATOR_SHA256,
        replayer.EXPECTED_SCHEMA_RECOVERY_VALIDATOR_SHA256,
        continuation.EXPECTED_SCHEMA_RECOVERY_VALIDATOR_SHA256,
        builder.EXPECTED_VALIDATOR_SHA256,
    } == {validator_sha256}
    assert verifier.CONTINUATION_PUBLIC_KEY == continuation.CONTINUATION_PUBLIC_KEY
    assert replayer.CONTINUATION_PUBLIC_KEY == continuation.CONTINUATION_PUBLIC_KEY
    assert verifier.EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256 == (
        continuation.EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert replayer.EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256 == (
        continuation.EXPECTED_CONTINUATION_PUBLIC_KEY_SHA256
    )
    assert continuation.CONTINUATION_PUBLIC_KEY != (
        continuation.RECOVERY_CONTINUATION_PUBLIC_KEY
    )
    assert continuation.CONTINUATION_PRIVATE_KEY.stat().st_mode & 0o7777 == 0o400
    assert (
        continuation.RECOVERY_CONTINUATION_PRIVATE_KEY.stat().st_mode & 0o7777
        == 0o400
    )


def test_v29_executable_alias_transient_target_mutation_fails_closed(
    modules, tmp_path
):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    _, target, record = make_executable_alias(tmp_path)
    original = target.read_bytes()
    try:
        continuation.bind_declared_relations(
            {"runtime": {"python": record}}, reader, historical_mode_records=set()
        )
        target.write_bytes(b"temporary replacement\n")
        target.write_bytes(original)
        with pytest.raises(continuation.ContinuationError, match="alias"):
            continuation.require_bootstrap_stable()
    finally:
        close_executable_alias_leases(continuation)
        reader.close()


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_resolved_path",
        "missing_symlink_target",
        "extra_field",
        "bool_size",
        "noncanonical_hash",
    ],
)
def test_v29_executable_alias_schema_ambiguity_fails_closed(
    modules, tmp_path, mutation
):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    _, _, record = make_executable_alias(tmp_path)
    if mutation == "missing_resolved_path":
        record.pop("resolved_path")
    elif mutation == "missing_symlink_target":
        record.pop("symlink_target")
    elif mutation == "extra_field":
        record["unexpected"] = True
    elif mutation == "bool_size":
        record["size_bytes"] = True
    else:
        record["sha256"] = "not-a-sha256"
    try:
        with pytest.raises(continuation.ContinuationError):
            continuation.bind_declared_relations(
                {"runtime": {"python": record}},
                reader,
                historical_mode_records=set(),
            )
    finally:
        close_executable_alias_leases(continuation)
        reader.close()


@pytest.mark.parametrize(
    "mutation", ["wrong_link_text", "regular_file_spoof", "target_bytes"]
)
def test_v29_executable_alias_initial_identity_drift_fails_closed(
    modules, tmp_path, mutation
):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    alias, target, record = make_executable_alias(tmp_path)
    if mutation == "wrong_link_text":
        record["symlink_target"] = "python3.10"
    elif mutation == "regular_file_spoof":
        alias.rename(tmp_path / "python.original-link")
        alias.write_bytes(target.read_bytes())
        alias.chmod(0o755)
    else:
        record["sha256"] = "0" * 64
    try:
        with pytest.raises(continuation.ContinuationError):
            continuation.bind_declared_relations(
                {"runtime": {"python": record}},
                reader,
                historical_mode_records=set(),
            )
    finally:
        close_executable_alias_leases(continuation)
        reader.close()


@pytest.mark.parametrize("mutation", ["alias_replaced", "target_replaced"])
def test_v29_executable_alias_post_bind_replacement_fails_closed(
    modules, tmp_path, mutation
):
    continuation = modules["continuation"]
    reader = modules["release_authority"].BoundArtifactReader()
    alias, target, record = make_executable_alias(tmp_path)
    try:
        continuation.bind_declared_relations(
            {"runtime": {"python": record}},
            reader,
            historical_mode_records=set(),
        )
        if mutation == "alias_replaced":
            alias.rename(tmp_path / "python.bound-link")
            alias.symlink_to(target.name)
        else:
            target.rename(tmp_path / "python3.11.bound-target")
            target.write_bytes(b"replacement python runtime\n")
            target.chmod(0o755)
        with pytest.raises(continuation.ContinuationError, match="alias"):
            continuation.require_bootstrap_stable()
    finally:
        close_executable_alias_leases(continuation)
        reader.close()


def test_authorization_and_promotion_transitions_must_match(modules):
    replayer = modules["replayer"]
    authorization = signed_chain_payload(
        "tumor_ref_source_receipt_promotion_authorization.v3.json"
    )
    promotion = signed_chain_payload("tumor_ref_source_receipt_promotion.v3.json")
    replayer.require_signed_focal_transitions_equal(authorization, promotion)
    promotion["evidence"]["focal_source_identity_transition"]["interpretation"] += (
        " drift"
    )
    with pytest.raises(replayer.GateError):
        replayer.require_signed_focal_transitions_equal(authorization, promotion)


def test_v1_walker_reproduces_the_recorded_mode_failure():
    replayer = load_module(
        "schema_recovery_replayer_v1_failure_canary",
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_recovery_v1.py",
    )
    payload = signed_chain_payload(
        "tumor_ref_source_receipt_promotion_authorization.v3.json"
    )
    reader = replayer.BoundArtifactReader()
    try:
        with pytest.raises(replayer.GateError, match="field=mode"):
            replayer.recursively_validate_artifact_relations(payload, reader)
    finally:
        reader.close()


def test_validator_static_contract_binds_legacy_chain(modules, monkeypatch):
    validator = modules["validator"]
    verifier = modules["verifier"]
    replayer = modules["replayer"]
    result = validator.static_contract_probe()
    assert result["pass"] is True
    assert result["historical_execution_identity_keys"] == validator.HISTORICAL_EXECUTION_KEYS
    assert result["current_identity_keys"] == validator.CURRENT_IDENTITY_KEYS
    assert "link_count" not in result["historical_execution_identity_keys"]
    assert "link_count" in result["current_identity_keys"]
    assert result["prior_recovery_chain"]["authority_signature_verified"] is True
    assert result["prior_recovery_chain"]["runner_output_created"] is False
    failed_v9 = result["prior_failed_signed_recovery"]["v9"]
    failed_v10 = result["prior_failed_signed_recovery"]["v10"]
    failed_v11 = result["prior_failed_signed_recovery"]["v11"]
    failed_v12 = result["prior_failed_signed_recovery"]["v12"]
    failed_v14 = result["prior_failed_signed_recovery"]["v14"]
    failed_v16 = result["prior_failed_signed_recovery"]["v16"]
    failed_v17 = result["prior_failed_signed_recovery"]["v17"]
    failed_v18 = result["prior_failed_signed_recovery"]["v18"]
    failed_v19 = result["prior_failed_signed_recovery"]["v19"]
    failed_v28 = result["prior_failed_signed_recovery"]["v28"]
    rejected_v19_round1 = result["rejected_unsigned_generations"]["v19_round1"]
    rejected_v29_round1 = result["rejected_unsigned_generations"]["v29_round1"]
    assert failed_v9["authority_signature_verified"] is True
    assert failed_v9["atomic_commit_verified"] is True
    assert failed_v9["runtime_outputs_created"] is False
    assert failed_v10["authority_signature_verified"] is True
    assert failed_v10["atomic_commit_verified"] is True
    assert failed_v10["verification_receipt_created"] is True
    assert failed_v10["replay_receipt_and_log_created"] is True
    assert failed_v10["continuation_outputs_created"] is False
    assert failed_v11["basic_to_full_stat_mismatch_reproduced"] is True
    assert failed_v11["continuation_outputs_created"] is False
    assert failed_v12["executable_alias_relation_mismatch_reproduced"] is True
    assert failed_v12["continuation_outputs_created"] is False
    assert failed_v14["metadata_only_relation_schema_mismatch_reproduced"] is True
    assert failed_v14["replay_outputs_created"] is False
    assert failed_v14["continuation_outputs_created"] is False
    assert failed_v16["authority_signature_verified"] is True
    assert failed_v16["atomic_commit_verified"] is True
    assert failed_v16["verification_receipt_created"] is True
    assert failed_v16["replay_receipt_and_log_created"] is True
    assert failed_v16["replay_success_witness_created"] is True
    assert failed_v16["legacy_eight_key_stat_schema_mismatch_reproduced"] is True
    assert failed_v16["continuation_outputs_created"] is False
    assert failed_v17["authority_signature_verified"] is True
    assert failed_v17["atomic_commit_verified"] is True
    assert failed_v17["verification_receipt_created"] is True
    assert failed_v17["replay_receipt_and_log_created"] is True
    assert failed_v17["replay_success_witness_created"] is True
    assert failed_v17["continuation_child_started"] is True
    assert failed_v17["fresh_verifier_passed"] is True
    assert failed_v17["metadata_enrichment_conflict_reproduced"] is True
    assert failed_v17["syscall_trace_evidence_bound"] is True
    assert failed_v17["continuation_outputs_created"] is False
    assert failed_v18["authority_signature_verified"] is True
    assert failed_v18["atomic_commit_verified"] is True
    assert failed_v18["verification_receipt_created"] is True
    assert failed_v18["replay_receipt_and_log_created"] is True
    assert failed_v18["replay_success_witness_created"] is True
    assert failed_v18["continuation_child_started"] is True
    assert failed_v18["fresh_verifier_passed"] is True
    assert failed_v18["tumor_ref_summary_alias_noncanonical_reproduced"] is True
    assert failed_v18["syscall_trace_evidence_bound"] is True
    assert failed_v18["continuation_outputs_created"] is False
    assert failed_v19["authority_signature_verified"] is True
    assert failed_v19["atomic_commit_verified"] is True
    assert failed_v19["verification_receipt_created"] is True
    assert failed_v19["replay_receipt_and_log_created"] is True
    assert failed_v19["replay_success_witness_created"] is True
    assert failed_v19["continuation_child_started"] is True
    assert failed_v19["fresh_verifier_passed"] is True
    assert failed_v19["mode000_inotify_eacces_reproduced"] is True
    assert failed_v19["syscall_trace_evidence_bound"] is True
    assert failed_v19["continuation_outputs_created"] is False
    assert failed_v28["authority_signature_verified"] is True
    assert failed_v28["dataset_release_signature_verified"] is True
    assert failed_v28["archived_inventory_count"] == 742
    assert failed_v28["failure"]["stage"] == "C_REPORT_BUILD_AFTER_SIGNED_DATASET"
    assert failed_v28["failure"]["scientific_payload_changed"] is False
    assert failed_v28["pass"] is True
    assert rejected_v19_round1["authority_created"] is False
    assert rejected_v19_round1["signature_created"] is False
    assert rejected_v19_round1["strictest_review_wins"] is True
    assert rejected_v19_round1["findings_corrected_in_active_candidate"] is True
    assert rejected_v29_round1["payload_inventory_count"] == 20
    assert rejected_v29_round1["mendel_initial_verdict"] == "APPROVE_SUPERSEDED"
    assert rejected_v29_round1["mendel_corrected_verdict"] == "REQUEST_CHANGES"
    assert rejected_v29_round1["nash_verdict"] == "REQUEST_CHANGES"
    assert rejected_v29_round1["strictest_reproducible_review_wins"] is True
    assert rejected_v29_round1["authority_created_before_correction"] is False
    assert rejected_v29_round1["any_key_consumed_before_correction"] is False

    recovery_sources = validator._records(validator.RECOVERY_SOURCE_PATHS)
    evidence = validator._build_validation_evidence(
        expected_runtime_role="continuation_verifier",
        recovery_sources=recovery_sources,
        authority_record={},
        signature_record={},
            commit_record={},
            public_record={},
            expected_private={"path": str(validator.PRIVATE_KEY), "mode": "0o0"},
            terminal_key_rotation=result["terminal_key_rotation"],
            terminal_key_state=result["terminal_key_state"],
        )
    assert evidence["authority_validator"] == recovery_sources["authority_validator"]
    reader = replayer.BoundArtifactReader()
    try:
        relations = replayer.recursively_validate_artifact_relations(
            {"authority_validator": evidence["authority_validator"]}, reader
        )
    finally:
        reader.close()
    assert str(validator.VALIDATOR) in relations

    public_key = validator._records({"public_key": validator.PUBLIC_KEY})["public_key"]
    prior_recovery_chain = validator._validate_prior_recovery_chain(
        validator._records(validator.PRIOR_RECOVERY_CHAIN_PATHS)
    )
    rejected_generations = validator._validate_rejected_generations()
    prior_failed_signed_recovery = validator._validate_prior_failed_signed_recoveries()
    fresh_key_bootstrap = validator._validate_v30_key_bootstrap()
    terminal_key_rotation, _ = validator._validate_terminal_key_rotation(
        expected_recovery_private_mode="0o400"
    )
    review = {
        "schema_name": validator.REVIEW_SCHEMA_NAME,
        "schema_version": validator.SCHEMA_VERSION,
        "reviewer": validator.EXPECTED_REVIEWERS["mendel"],
        "reviewer_agent_id": validator.EXPECTED_REVIEWER_AGENT_IDS["mendel"],
        "attribution": {
            "type": "orchestrator_recorded_transport",
            "transport": "multi_agent_v1",
            "transport_id": validator.EXPECTED_REVIEWER_AGENT_IDS["mendel"],
            "cryptographic_reviewer_authorship_proven": False,
            "semantics": validator.REVIEW_ATTRIBUTION_SEMANTICS,
        },
        "review_completed_at_utc": "2026-07-22T18:30:00Z",
        "verdict": "APPROVE",
        "reviewed_source_set_sha256": validator._json_sha256(recovery_sources),
        "legacy_source_set_sha256": validator._json_sha256(
            validator._records(validator.LEGACY_SOURCE_PATHS)
        ),
        "prior_recovery_chain_sha256": validator._json_sha256(prior_recovery_chain),
        "rejected_generations_sha256": validator._json_sha256(
            rejected_generations
        ),
        "prior_failed_signed_recovery_sha256": validator._json_sha256(
            prior_failed_signed_recovery
        ),
        "fresh_key_bootstrap_sha256": validator._json_sha256(
            fresh_key_bootstrap
        ),
        "trusted_recovery_public_key": public_key,
        "scope_sha256": validator._json_sha256(validator.RECOVERY_SCOPE),
        "terminal_key_rotation_sha256": validator._json_sha256(
            terminal_key_rotation
        ),
        "readonly_probe": {
            "command": validator.READONLY_PROBE_COMMAND,
            "exit_code": 0,
            "no_output_writes": True,
            "status": "PASS",
            "forbidden_output_slots_checked": (
                validator.EXPECTED_CEREMONY_FORBIDDEN_OUTPUT_SLOT_COUNT
            ),
            "regression_summary": validator.EXPECTED_REGRESSION_SUMMARY,
        },
        "high_findings": [],
        "medium_findings": [],
        "low_findings": [],
        "unresolved_conditions": [],
        "summary": "Regression review approves the exact v30 source set.",
        "pass": True,
    }
    legacy_sources = validator._records(validator.LEGACY_SOURCE_PATHS)
    assert (
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
        == validator.EXPECTED_REVIEWER_AGENT_IDS["mendel"]
    )
    review_without_prior_failed_digest = copy.deepcopy(review)
    review_without_prior_failed_digest.pop("prior_failed_signed_recovery_sha256")
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review_without_prior_failed_digest,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    review_with_wrong_prior_failed_digest = copy.deepcopy(review)
    review_with_wrong_prior_failed_digest["prior_failed_signed_recovery_sha256"] = (
        "0" * 64
    )
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review_with_wrong_prior_failed_digest,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    mutated_prior_failed_signed_recovery = copy.deepcopy(prior_failed_signed_recovery)
    mutated_prior_failed_signed_recovery["v9"]["pass"] = False
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            mutated_prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    review["reviewer_agent_id"] = "stale-reviewer-id"
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    review["reviewer_agent_id"] = validator.EXPECTED_REVIEWER_AGENT_IDS["mendel"]
    review["attribution"]["cryptographic_reviewer_authorship_proven"] = True
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    review["attribution"]["cryptographic_reviewer_authorship_proven"] = False
    review["readonly_probe"]["command"] = ["/usr/bin/true"]
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    review["readonly_probe"]["command"] = list(validator.READONLY_PROBE_COMMAND)
    review["readonly_probe"]["exit_code"] = False
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    review["readonly_probe"]["exit_code"] = 0
    review["readonly_probe"]["no_output_writes"] = 0
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    review["readonly_probe"]["no_output_writes"] = True
    review["scope_sha256"] = "0" * 64
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_review(
            "mendel",
            review,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            prior_failed_signed_recovery,
            fresh_key_bootstrap,
            public_key,
            terminal_key_rotation,
        )
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._load_json(b'{"pass":true,"pass":true}', "duplicate-key canary")
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._load_json(b'{"value":NaN}', "non-finite canary")
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._load_json(b'{"value":1e400}', "overflow-float canary")

    private_record, private_fd = validator._open_retired_private_key_bound(
        verifier.AUTHORIZATION_PRIVATE_KEY
    )
    assert private_record == {
        "path": str(verifier.AUTHORIZATION_PRIVATE_KEY),
        "mode": "0o0",
    }
    real_fstat = validator.os.fstat

    def drift_private_mode(descriptor):
        observed = real_fstat(descriptor)
        if descriptor != private_fd:
            return observed
        values = list(observed)
        values[0] = (observed.st_mode & ~0o777) | 0o400
        return os.stat_result(values)

    monkeypatch.setattr(validator.os, "fstat", drift_private_mode)
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._require_leases()


def real_v10_verification_contract(modules):
    continuation = modules["continuation"]
    failed_root = (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v10_signed_authority_c_verification_incident_projection_mismatch"
    )
    payload = json.loads(
        (failed_root / "tumor_ref_source_receipt_promotion_verification.recovery.v10.json")
        .read_text(encoding="ascii")
    )
    payload["schema_recovery_authority"]["authority_id"] = (
        "20260724_tumor_ref_schema_recovery_v30"
    )
    authorization = signed_chain_payload(
        "tumor_ref_source_receipt_promotion_authorization.v3.json"
    )
    prepare = payload["producer_exit_attestations"]["prepare"]
    promote = payload["producer_exit_attestations"]["promote"]
    chain = {
        "authorization": authorization,
        "authorization_record": payload["authorization"]["artifact"],
        "authorization_signature_record": payload["authorization"]["signature"],
        "authorization_public_key_record": prepare["public_key"],
        "authorization_private_key_record": payload["private_key_retirement"][
            "authorization"
        ],
        "completion_record": payload["completion"]["artifact"],
        "completion_signature_record": payload["completion"]["signature"],
        "completion_public_key_record": promote["public_key"],
        "completion_private_key_record": payload["private_key_retirement"][
            "completion"
        ],
        "prepare_exit_attestation_record": prepare["artifact"],
        "promote_exit_attestation_record": promote["artifact"],
        "source_record": payload["source_receipt"],
        "canonical_record": payload["canonical_receipt"],
        "trusted_keys": payload["trusted_signing_keys"],
        "authority_record": payload["release_source_authority_validator"],
        "schema_recovery_authority": payload["schema_recovery_authority"],
        "verifier_record": payload["schema_recovery_authority"]["runtime_source"],
    }
    return continuation, payload, authorization, chain


def test_v29_continuation_accepts_real_v10_verifier_incident_projection(modules):
    continuation, payload, authorization, chain = real_v10_verification_contract(
        modules
    )
    archive = authorization["evidence"]["failed_attempt_archive"]
    assert set(archive) == {"receipt", "archived_log", "archived_cache"}
    assert payload["historical_incident_disclosure"]["archive_receipt"] == (
        archive["receipt"]
    )
    continuation.validate_verification_payload(
        payload,
        expected_mode="verify-and-record",
        chain=chain,
    )


@pytest.mark.parametrize("mutation", ["full_archive_object", "receipt_hash"])
def test_v29_continuation_rejects_incident_projection_drift(modules, mutation):
    continuation, payload, authorization, chain = real_v10_verification_contract(
        modules
    )
    payload = copy.deepcopy(payload)
    if mutation == "full_archive_object":
        payload["historical_incident_disclosure"]["archive_receipt"] = (
            authorization["evidence"]["failed_attempt_archive"]
        )
    else:
        payload["historical_incident_disclosure"]["archive_receipt"]["sha256"] = (
            "0" * 64
        )
    with pytest.raises(continuation.ContinuationError, match="payload contract drift"):
        continuation.validate_verification_payload(
            payload,
            expected_mode="verify-and-record",
            chain=chain,
        )


def test_v29_preserves_canonical_v10_scientific_inputs(modules):
    expected_key_parent = "20260722_all_ssnv_v10_strict_command_parity_bootstrap"
    verifier_key = modules["verifier"].EXPECTED_V7_PUBLIC_KEY_PATH
    for public_key in (
        verifier_key,
        modules["replayer"].V7_PUBLIC_KEY,
        modules["continuation"].V7_PUBLIC_KEY,
    ):
        assert public_key.parent.name == expected_key_parent
        assert public_key.is_file()

    screen = modules["continuation"].SCREEN
    assert screen.name == (
        "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full"
    )
    assert screen.is_dir()


def test_verifier_uses_distinct_legacy_and_recovery_receipts(modules):
    verifier = modules["verifier"]
    assert verifier.AUTHORIZED_VERIFIER.name == "verify_tumor_ref_receipt_promotion_v2.py"
    assert verifier.VERIFIER.name == "verify_tumor_ref_receipt_promotion_recovery_v30.py"
    assert verifier.VERIFICATION_RECEIPT != verifier.RECOVERY_VERIFICATION_RECEIPT
    assert verifier.RECOVERY_VERIFICATION_SCHEMA.endswith(".recovery")
    assert verifier.RECOVERY_SCHEMA_VERSION == "1.0.0"
    assert verifier.RECOVERY_VERIFICATION_RECEIPT.name.endswith("recovery.v30.json")


def test_v29_failed_v9_signed_generation_is_cryptographically_bound(modules):
    validator = modules["validator"]
    result = validator._validate_failed_v9_signed_recovery()
    assert result["pass"] is True
    assert result["authority_signature_verified"] is True
    assert result["atomic_commit_verified"] is True
    assert result["historical_authorized_output_slot_count"] == 17
    assert result["current_recovery_output_slot_count"] == 23


def test_v29_failed_v10_signed_generation_is_cryptographically_bound(modules):
    validator = modules["validator"]
    result = validator._validate_failed_v10_signed_recovery()
    assert result["pass"] is True
    assert result["authority_signature_verified"] is True
    assert result["atomic_commit_verified"] is True
    assert result["verification_receipt_created"] is True
    assert result["replay_receipt_and_log_created"] is True
    assert result["continuation_outputs_created"] is False
    assert result["canonical_downstream_outputs_created"] is False


def test_v29_failed_v25_signed_generation_and_lineage_defect_are_bound(modules):
    validator = modules["validator"]
    result = validator._validate_failed_v25_signed_recovery()
    assert result["pass"] is True
    assert result["failure"] == {
        "stage": "C_FINAL_DATASET_INTEGRATION",
        "contract": "tumor_ref_v1_declared_vs_v6_canonical_audit_path",
        "scientific_payload_changed": False,
        "terminal_signature_created": False,
    }
    assert result["authority_retired_private_key"]["mode"] == "0o0"
    assert result["terminal_unused_private_key"]["mode"] == "0o400"
    assert set(validator.FAILED_V25_ORIGINAL_OUTPUT_SLOTS).issubset(
        validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS
    )
    assert not any(
        os.path.lexists(path) for path in validator.FAILED_V25_ORIGINAL_OUTPUT_SLOTS
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "transaction_id",
        "member_set",
        "authority_digest",
        "retired_private_key",
    ],
)
def test_v29_failed_v25_atomic_commit_mutations_fail_closed(
    modules, monkeypatch, mutation
):
    validator = modules["validator"]
    real_load = validator._load_json

    def mutated_load(data, label):
        payload = copy.deepcopy(real_load(data, label))
        if label != "failed v25 authority commit":
            return payload
        if mutation == "transaction_id":
            payload["transaction_id"] = "0" * 32
        elif mutation == "member_set":
            payload["members"] = ["authority.json", "commit.json"]
        elif mutation == "authority_digest":
            payload["authority"]["sha256"] = "0" * 64
        else:
            payload["retired_private_key"]["path"] += ".substituted"
        return payload

    monkeypatch.setattr(validator, "_load_json", mutated_load)
    with pytest.raises(
        validator.RecoveryAuthorityError,
        match="Failed v25 signed-generation contract drift",
    ):
        validator._validate_failed_v25_signed_recovery()


def test_v29_failed_v11_signed_generation_and_repair_are_bound(modules):
    validator = modules["validator"]
    replayer = modules["replayer"]
    builder = modules["builder"]
    result = validator._validate_failed_v11_signed_recovery()
    assert result["pass"] is True
    assert result["authority_signature_verified"] is True
    assert result["atomic_commit_verified"] is True
    assert result["basic_to_full_stat_mismatch_reproduced"] is True
    assert result["continuation_outputs_created"] is False
    assert result["canonical_downstream_outputs_created"] is False
    assert len(validator.FAILED_V11_OUTPUT_SLOTS) == 13
    assert set(validator.FAILED_V11_OUTPUT_SLOTS).issubset(
        set(validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS)
    )
    assert not any(os.path.lexists(path) for path in validator.FAILED_V11_OUTPUT_SLOTS)
    source = inspect.getsource(replayer)
    assert '"historical_source_receipt": source_full_record' in source
    assert '"canonical_source_receipt": canonical_full_record' in source
    required = builder.required_module_paths(validator)
    assert {
        validator.FAILED_V11_EVIDENCE,
        validator.FAILED_V11_PUBLIC_KEY,
        validator.FAILED_V11_PRIVATE_KEY,
        validator.FAILED_V11_CONTINUATION_PUBLIC_KEY,
        validator.FAILED_V11_CONTINUATION_PRIVATE_KEY,
        *validator.FAILED_V11_SOURCE_PATHS.values(),
        *validator.FAILED_V11_ARCHIVED_ARTIFACT_PATHS.values(),
    }.issubset(required)


def real_v11_replay_source_records():
    failed_root = (
        AUDIT_ROOT
        / "failed_formal_runs"
        / "20260723_v11_signed_authority_c_replay_source_stat_projection_mismatch"
    )
    verification = json.loads(
        (failed_root / "tumor_ref_source_receipt_promotion_verification.recovery.v11.json")
        .read_text(encoding="ascii")
    )
    replay = json.loads(
        (failed_root / "m2v5_runner_only_gate_replay.recovery.v11.json").read_text(
            encoding="ascii"
        )
    )
    return verification, replay["promotion_trust_chain"]


@pytest.mark.parametrize(
    "field,verification_field",
    [
        ("historical_source_receipt", "source_receipt"),
        ("canonical_source_receipt", "canonical_receipt"),
    ],
)
def test_v29_real_v11_replay_exposes_each_basic_to_full_stat_mismatch(
    modules, field, verification_field
):
    continuation = modules["continuation"]
    verification, trust_chain = real_v11_replay_source_records()
    repaired = copy.deepcopy(trust_chain)
    repaired["historical_source_receipt"] = verification["source_receipt"]
    repaired["canonical_source_receipt"] = verification["canonical_receipt"]
    observed = trust_chain[field]
    expected = verification[verification_field]
    assert set(observed) == set(continuation.BASIC_ARTIFACT_KEYS)
    assert set(expected) == set(continuation.STAT_ARTIFACT_KEYS)
    assert set(expected) - set(observed) == {
        "device",
        "inode",
        "link_count",
        "mtime_ns",
        "ctime_ns",
    }
    broken = copy.deepcopy(repaired)
    broken[field] = observed
    with pytest.raises(continuation.ContinuationError, match=f"{field} full-stat record"):
        continuation.validate_replay_source_identity_contract(
            broken,
            verification["source_receipt"],
            verification["canonical_receipt"],
        )
    continuation.validate_replay_source_identity_contract(
        repaired,
        verification["source_receipt"],
        verification["canonical_receipt"],
    )


@pytest.mark.parametrize(
    "mutation", ["missing_stat_field", "extra_field", "wrong_stat_value"]
)
def test_v29_replay_full_stat_mutations_fail_closed(modules, mutation):
    continuation = modules["continuation"]
    verification, trust_chain = real_v11_replay_source_records()
    repaired = copy.deepcopy(trust_chain)
    repaired["historical_source_receipt"] = copy.deepcopy(verification["source_receipt"])
    repaired["canonical_source_receipt"] = copy.deepcopy(
        verification["canonical_receipt"]
    )
    if mutation == "missing_stat_field":
        repaired["historical_source_receipt"].pop("device")
        expected_label = "historical_source_receipt"
    elif mutation == "extra_field":
        repaired["canonical_source_receipt"]["unexpected"] = 1
        expected_label = "canonical_source_receipt"
    else:
        repaired["historical_source_receipt"]["device"] += 1
        expected_label = "historical_source_receipt"
    expected_error = (
        f"{expected_label} live identity drift"
        if mutation == "wrong_stat_value"
        else f"{expected_label} full-stat record"
    )
    with pytest.raises(
        continuation.ContinuationError,
        match=expected_error,
    ):
        continuation.validate_replay_source_identity_contract(
            repaired,
            verification["source_receipt"],
            verification["canonical_receipt"],
        )


def test_v29_continuation_retires_key_before_signing_and_signature_publication(
    modules,
    monkeypatch,
    tmp_path,
):
    continuation = modules["continuation"]
    source = inspect.getsource(continuation.sign_exit_attestation)
    staging = source.index(
        "private_signing_fd = stage_private_key_in_sealed_memfd(private_fd)"
    )
    retirement = source.index("os.fchmod(private_fd, 0o000)")
    signing = source.index("result = subprocess.run(")
    signing_key_close = source.index("os.close(private_signing_fd)")
    publication = source.index("link_fd_no_replace(")
    assert source.count("os.fchmod(private_fd, 0o000)") == 1
    assert staging < retirement < signing < signing_key_close < publication
    assert 'f"/proc/self/fd/{private_signing_fd}"' in source
    helper = inspect.getsource(continuation.stage_private_key_in_sealed_memfd)
    assert "LINUX_F_ADD_SEALS" in helper
    assert "LINUX_F_GET_SEALS" in helper
    assert "staged.st_nlink != 0" in helper
    assert "stat.S_IMODE(staged.st_mode) != 0o400" in helper
    assert (
        continuation.TERMINAL_GOVERNANCE[
            "published_terminal_signature_is_provisional_until_independent_verifier"
        ]
        is True
    )
    assert (
        continuation.EXIT_ATTESTATION_CHECKS[
            "post_link_supervisor_failure_requires_continuation_incident"
        ]
        is True
    )
    witness_commit_source = inspect.getsource(
        continuation.commit_terminal_success_witness
    )
    final_link = witness_commit_source.rindex("link_fd_no_replace(")
    no_return_exit = witness_commit_source.rindex("os._exit(0)")
    assert final_link < no_return_exit
    assert "os.fsync" not in witness_commit_source[final_link:]
    assert "precommit()" not in witness_commit_source[final_link:]
    assert "pthread_sigmask" in witness_commit_source[:final_link]

    private_path, public_path = generate_ed25519_keypair(
        continuation.OPENSSL, tmp_path, "continuation"
    )
    message_path = tmp_path / "continuation-attestation.json"
    signature_path = tmp_path / "continuation-attestation.ed25519.sig"
    message_path.write_bytes(b'{"supervised_child_exit":0}\n')
    private_fd = os.open(private_path, os.O_RDONLY | os.O_CLOEXEC)
    message_fd = os.open(message_path, os.O_RDONLY | os.O_CLOEXEC)
    staged_private_fd = -1
    try:
        staged_private_fd = continuation.stage_private_key_in_sealed_memfd(
            private_fd
        )
        staged = os.fstat(staged_private_fd)
        assert staged.st_nlink == 0
        assert staged.st_mode & 0o7777 == 0o400
        assert (
            fcntl.fcntl(staged_private_fd, continuation.LINUX_F_GET_SEALS)
            == continuation.supervisor_capability_seal_mask()
        )
        os.fchmod(private_fd, 0o000)
        os.fsync(private_fd)
        assert os.fstat(private_fd).st_mode & 0o7777 == 0
        result = subprocess.run(
            [
                str(continuation.OPENSSL),
                "pkeyutl",
                "-sign",
                "-rawin",
                "-inkey",
                f"/proc/self/fd/{staged_private_fd}",
                "-in",
                f"/proc/self/fd/{message_fd}",
            ],
            pass_fds=(staged_private_fd, message_fd),
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0
        assert result.stderr == b""
        assert len(result.stdout) == 64
        signature_path.write_bytes(result.stdout)
        verified = subprocess.run(
            [
                str(continuation.OPENSSL),
                "pkeyutl",
                "-verify",
                "-pubin",
                "-inkey",
                str(public_path),
                "-rawin",
                "-in",
                str(message_path),
                "-sigfile",
                str(signature_path),
            ],
            capture_output=True,
            check=False,
        )
        assert verified.returncode == 0
    finally:
        if staged_private_fd >= 0:
            os.close(staged_private_fd)
        os.close(message_fd)
        os.close(private_fd)

    present = {continuation.CONTINUATION_SIGNATURE}
    monkeypatch.setattr(
        continuation.os.path, "lexists", lambda path: path in present
    )

    def failed_incident_commit(*args, **kwargs):
        raise OSError("synthetic incident publication failure")

    monkeypatch.setattr(
        continuation, "commit_terminal_receipt", failed_incident_commit
    )
    continuation.record_continuation_incident(RuntimeError("post-link failure"))
    with pytest.raises(
        continuation.ContinuationError,
        match="success witness is absent",
    ):
        continuation.require_success_witness_state(required=True)
    continuation.require_success_witness_state(required=False)
    present.add(continuation.CONTINUATION_SUCCESS_WITNESS)
    continuation.require_success_witness_state(required=True)
    with pytest.raises(
        continuation.ContinuationError,
        match="already exists before prewitness",
    ):
        continuation.require_success_witness_state(required=False)


def test_v29_signed_crosslinks_use_exact_historical_17_slots(modules):
    verifier = modules["verifier"]
    validator = modules["validator"]
    authorization = signed_chain_payload(
        "tumor_ref_source_receipt_promotion_authorization.v3.json"
    )
    completion = signed_chain_payload("tumor_ref_source_receipt_promotion.v3.json")
    authorized = [str(path) for path in verifier.AUTHORIZED_DOWNSTREAM_OUTPUTS]
    current = [str(path) for path in verifier.DOWNSTREAM_OUTPUTS]
    assert len(authorized) == 17
    assert len(current) == 25
    assert validator.RECOVERY_SCOPE["current_recovery_output_slot_count"] == len(
        verifier.DOWNSTREAM_OUTPUTS
    )
    assert validator.RECOVERY_SCOPE["ceremony_forbidden_output_slot_count"] == len(
        validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS
    )
    assert validator.RECOVERY_SCOPE["ceremony_staging_pattern_count"] == len(
        validator.CEREMONY_STAGING_PATTERNS
    )
    assert validator.RECOVERY_SCOPE["rejected_unsigned_generations"][-3:] == [
        "v24",
        "v27",
        "v29_round1",
    ]
    assert validator.AUTHORITY_CHECKS[
        "all_439_forbidden_slots_have_occupied_state_regressions"
    ] is True
    assert validator.AUTHORITY_CHECKS[
        "all_29_staging_patterns_have_occupied_state_regressions"
    ] is True
    assert (
        validator.RECOVERY_SCOPE["ceremony_forbidden_output_slots_sha256"]
        == validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS_SHA256
    )
    assert (
        validator.RECOVERY_SCOPE["ceremony_staging_patterns_sha256"]
        == validator.CEREMONY_STAGING_PATTERNS_SHA256
    )
    assert {
        "exact_typed_relation_schema_and_per_reader_path_registry",
        "replay_process_lifetime_alias_and_target_inotify_monitor",
        "runner_supervisor_waitpid_normal_exit_zero_success_witness",
        "persistent_final_absence_inotify_watch_through_pre_rename_check",
    }.issubset(set(validator.RECOVERY_SCOPE["allowed_changes"]))
    assert current[:17] == [str(path) for path in verifier.CURRENT_DOWNSTREAM_OUTPUTS]
    assert current[:17] == authorized
    assert authorization["downstream_output_absence"] == authorized
    assert completion["governance"]["downstream_output_slots_absent_at_promotion"] == (
        authorized
    )


def test_v29_verifier_separates_crosslink_and_live_absence_in_source(modules):
    verifier = modules["verifier"]
    source = inspect.getsource(verifier._verify)
    assert source.count("for path in AUTHORIZED_DOWNSTREAM_OUTPUTS") == 2
    assert source.count("_require_downstream_output_state(mode)") == 3
    state_source = inspect.getsource(verifier._require_downstream_output_state)
    assert "RECOVERY_REPLAY_OUTPUTS" in state_source
    assert "DOWNSTREAM_OUTPUTS" in state_source


def test_v29_downstream_output_stage_state_machine(modules, monkeypatch):
    verifier = modules["verifier"]
    present: set[Path] = set()
    monkeypatch.setattr(verifier.os.path, "lexists", lambda path: path in present)

    verifier._require_downstream_output_state("--verify-and-record")
    present.add(verifier.DOWNSTREAM_OUTPUTS[0])
    with pytest.raises(verifier.VerificationError, match="occupied_forbidden"):
        verifier._require_downstream_output_state("--verify-and-record")

    present.clear()
    verifier._require_downstream_output_state("--verify")
    present.add(verifier.RECOVERY_REPLAY_OUTPUTS[0])
    with pytest.raises(verifier.VerificationError, match="partial_replay_evidence"):
        verifier._require_downstream_output_state("--verify")
    present.add(verifier.RECOVERY_REPLAY_OUTPUTS[1])
    with pytest.raises(verifier.VerificationError, match="partial_replay_evidence"):
        verifier._require_downstream_output_state("--verify")
    present.add(verifier.RECOVERY_REPLAY_OUTPUTS[2])
    verifier._require_downstream_output_state("--verify")
    present.add(verifier.CURRENT_DOWNSTREAM_OUTPUTS[0])
    with pytest.raises(verifier.VerificationError, match="occupied_forbidden"):
        verifier._require_downstream_output_state("--verify")
    with pytest.raises(verifier.VerificationError, match="Unknown"):
        verifier._require_downstream_output_state("--invalid")


def test_v29_replayer_verifies_before_receipt_and_parent_witness(modules):
    replayer = modules["replayer"]
    source = inspect.getsource(replayer.worker_main)
    preflight_absence = source.index(
        'require_absent(OUTPUT_LOG, "runner gate log")'
    )
    verifier_call = source.index("verifier_result = subprocess.run(")
    log_commit = source.index("log_record = write_exclusive_readonly(")
    receipt_commit = source.rindex("write_exclusive_readonly(")
    assert preflight_absence < verifier_call < log_commit < receipt_commit
    supervisor_source = inspect.getsource(replayer.supervise_worker)
    waitpid_check = supervisor_source.index("require_worker_normal_exit(")
    witness_commit = supervisor_source.rindex("write_exclusive_readonly(")
    assert waitpid_check < witness_commit
    assert "return supervise_worker()" in inspect.getsource(replayer.main)


def test_v29_failed_v9_original_slots_are_forbidden_and_absent(modules):
    validator = modules["validator"]
    forbidden = set(validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS)
    assert len(validator.FAILED_V9_OUTPUT_SLOTS) == 13
    assert set(validator.FAILED_V9_OUTPUT_SLOTS).issubset(forbidden)
    assert not any(os.path.lexists(path) for path in validator.FAILED_V9_OUTPUT_SLOTS)


def test_v29_failed_v10_original_slots_are_forbidden_and_absent(modules):
    validator = modules["validator"]
    forbidden = set(validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS)
    assert len(validator.FAILED_V10_OUTPUT_SLOTS) == 13
    assert set(validator.FAILED_V10_OUTPUT_SLOTS).issubset(forbidden)
    assert not any(os.path.lexists(path) for path in validator.FAILED_V10_OUTPUT_SLOTS)


def test_v29_failed_v26_signed_generation_is_bound_and_quarantined(modules):
    validator = modules["validator"]
    observed = validator._validate_failed_v26_signed_recovery()
    assert observed["pass"] is True
    assert observed["failure"] == {
        "stage": "C_PRE_DOWNSTREAM_VALIDATE_PROMOTION_CHAIN",
        "contract": "historical_signed_runtime_projection_exact_11_roles",
        "current_reviewed_runtime_role_count": 14,
        "scientific_payload_changed": False,
        "terminal_signature_created": False,
    }
    assert set(validator.FAILED_V26_ORIGINAL_OUTPUT_SLOTS).issubset(
        set(validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS)
    )
    assert not any(
        os.path.lexists(path) for path in validator.FAILED_V26_ORIGINAL_OUTPUT_SLOTS
    )


def test_v29_failed_v28_signed_dataset_generation_is_bound_and_quarantined(modules):
    validator = modules["validator"]
    observed = validator._validate_failed_v28_signed_recovery()
    assert observed["pass"] is True
    assert observed["archived_inventory_count"] == 742
    assert observed["authority_signature_verified"] is True
    assert observed["dataset_release_signature_verified"] is True
    assert observed["failure"] == {
        "stage": "C_REPORT_BUILD_AFTER_SIGNED_DATASET",
        "contract": "json_mapping_exact_key_set_not_insertion_order",
        "scientific_payload_changed": False,
        "terminal_signature_created": False,
        "report_signature_created": False,
    }
    assert {
        role: record["state"] for role, record in observed["key_archives"].items()
    } == {
        "authority_v28": "RETIRED_AFTER_SINGLE_SIGNED_AUTHORITY_THAT_FAILED_AT_C",
        "terminal_v18": "UNUSED_NO_TERMINAL_SIGNATURE_RETIRED_FROM_REUSE",
        "result_v5": "CONSUMED_DATASET_SIGNATURE_PROVISIONAL_AFTER_C_FAILURE",
        "report_v5": "UNUSED_REPORT_KEY_RETIRED_FROM_FAILED_GENERATION",
    }
    assert set(validator.FAILED_V28_ORIGINAL_OUTPUT_SLOTS).issubset(
        set(validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS)
    )
    assert not any(
        os.path.lexists(path) for path in validator.FAILED_V28_ORIGINAL_OUTPUT_SLOTS
    )


def test_v30_probe_and_validator_guard_identical_439_slots(modules):
    validator = modules["validator"]
    probe = modules["probe"]
    assert len(validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS) == 439
    assert set(probe.FORBIDDEN_OUTPUT_SLOTS) == set(
        validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS
    )
    summary = probe.EXPECTED_REGRESSION_SUMMARY
    count = int(summary.split()[0])
    assert probe.parse_exact_pytest_summary(
        f"{summary} in 7.47s\n".encode("ascii")
    ) == summary
    for invalid in (
        f"{count + 1000} passed in 7.47s\n".encode("ascii"),
        f"{summary}, 1 xfailed in 7.47s\n".encode("ascii"),
        f"{summary} in 7.47s\n1 warning in 0.01s\n".encode("ascii"),
    ):
        assert probe.parse_exact_pytest_summary(invalid) is None


def test_v29_rejected_v27_v24_v23_v22_and_v21_slots_are_bound_and_watched(modules):
    validator = modules["validator"]
    observed = validator._validate_rejected_v22_pre_authority()
    observed_v23 = validator._validate_rejected_v23_pre_authority()
    observed_v24 = validator._validate_rejected_v24_pre_authority()
    observed_v27 = validator._validate_rejected_v27_pre_authority()
    forbidden = set(validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS)
    assert observed["authority_created"] is False
    assert observed["strictest_review_wins"] is True
    assert observed["review_authorship_cryptographically_proven"] is False
    assert observed["outputs_remain_absent"] is True
    assert observed["pass"] is True
    assert observed_v23["authority_created"] is False
    assert observed_v23["strictest_review_wins"] is True
    assert observed_v23["review_authorship_cryptographically_proven"] is False
    assert observed_v23["outputs_remain_absent"] is True
    assert observed_v23["pass"] is True
    assert observed_v24["authority_created"] is False
    assert observed_v24["strictest_review_wins"] is True
    assert observed_v24["review_authorship_cryptographically_proven"] is False
    assert observed_v24["outputs_remain_absent"] is True
    assert observed_v24["pass"] is True
    assert observed_v27["authority_created"] is False
    assert observed_v27["strictest_reproducible_finding_wins"] is True
    assert observed_v27["review_authorship_cryptographically_proven"] is False
    assert observed_v27["outputs_remain_absent"] is True
    assert observed_v27["pass"] is True
    assert set(validator.REJECTED_V22_OUTPUT_SLOTS).issubset(forbidden)
    assert set(validator.REJECTED_V23_OUTPUT_SLOTS).issubset(forbidden)
    assert set(validator.REJECTED_V24_OUTPUT_SLOTS).issubset(forbidden)
    assert set(validator.REJECTED_V27_OUTPUT_SLOTS).issubset(forbidden)
    assert set(validator.FAILED_V21_ORIGINAL_OBSERVED_OUTPUT_SLOTS).issubset(
        forbidden
    )


def test_v29_rejected_v13_pre_authority_evidence_is_bound(modules):
    validator = modules["validator"]
    probe = modules["probe"]
    observed = validator._validate_rejected_v13_pre_authority()
    assert observed["authority_created"] is False
    assert observed["formal_reviews_published"] is False
    assert observed["signatures_created"] is False
    assert observed["outputs_remain_absent"] is True
    assert observed["blocking_finding_count"] == 2
    assert observed["pass"] is True
    assert set(probe.FORBIDDEN_OUTPUT_SLOTS) == set(
        validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS
    )


def test_v29_rejected_v15_pre_authority_evidence_is_bound(modules):
    validator = modules["validator"]
    observed = validator._validate_rejected_v15_pre_authority()
    assert observed["authority_created"] is False
    assert observed["formal_reviews_published"] is False
    assert observed["signatures_created"] is False
    assert observed["outputs_remain_absent"] is True
    assert observed["strictest_review_wins"] is True
    assert observed["pass"] is True
    assert validator.RECOVERY_SCOPE["rejected_v15_authority_created"] is False
    assert validator.RECOVERY_SCOPE["rejected_v15_formal_reviews_published"] is False
    assert validator.RECOVERY_SCOPE[
        "rejected_v15_keys_archived_unused_never_signed"
    ] is True


def test_v19_round1_request_changes_are_archived_and_bound(modules):
    validator = modules["validator"]
    observed = validator._validate_rejected_v19_round1_pre_authority()
    assert observed["authority_created"] is False
    assert observed["signature_created"] is False
    assert observed["strictest_review_wins"] is True
    assert observed["findings_corrected_in_active_candidate"] is True
    assert observed["pass"] is True
    assert set(observed["candidate_reviews"]) == {"mendel", "nash"}
    assert observed["sources"]["downstream_continuation"]["sha256"] == (
        "ade2f9b9331b709487e06ca9a72bf46a248f7dc196edc7f01ec655d2472ccc2f"
    )
    assert hashlib.sha256(validator.RECOVERY_CONTINUATION.read_bytes()).hexdigest() != (
        observed["sources"]["downstream_continuation"]["sha256"]
    )
    assert validator.RECOVERY_SCOPE[
        "rejected_v19_round1_findings_corrected_in_active_candidate"
    ] is True


def test_v29_round1_terminal_projection_rejection_is_archived_and_bound(modules):
    validator = modules["validator"]
    continuation = modules["continuation"]
    probe = modules["probe"]
    observed = validator._validate_rejected_v29_round1_pre_authority()
    assert observed["payload_inventory_count"] == 20
    assert observed["mendel_initial_verdict"] == "APPROVE_SUPERSEDED"
    assert observed["mendel_corrected_verdict"] == "REQUEST_CHANGES"
    assert observed["nash_verdict"] == "REQUEST_CHANGES"
    assert observed["external_claude_verdict"] == "APPROVE_SUPERSEDED"
    assert observed["strictest_reproducible_review_wins"] is True
    assert observed["authority_created_before_correction"] is False
    assert observed["any_key_consumed_before_correction"] is False
    assert observed["scientific_payload_changed"] is False
    assert observed["claim_ceiling_changed"] is False
    assert observed["pass"] is True
    probe_observed = probe.validate_rejected_v29_round1_archive()
    assert probe_observed["payload_inventory_count"] == 20
    assert probe_observed["strictest_reproducible_review_wins"] is True
    assert probe_observed["pass"] is True
    archived_continuation = observed["payload_inventory"][
        "source_snapshot/continue_m2v5_after_tumor_ref_promotion_recovery_v29.py"
    ]
    assert archived_continuation["sha256"] == (
        "e6c4180a8c9fb3e7a3dd87e767c9cc6f12b6e76706ba75edfd4c71bec7094294"
    )
    assert hashlib.sha256(validator.RECOVERY_CONTINUATION.read_bytes()).hexdigest() != (
        archived_continuation["sha256"]
    )
    assert continuation.FAILED_V28_CONTINUATION_PRIVATE_KEY in (
        continuation.PRIVATE_KEY_PATHS
    )
    assert validator.RECOVERY_SCOPE[
        "rejected_v29_round1_terminal_projection_gap_reproduced"
    ] is True
    assert validator.AUTHORITY_CHECKS[
        "rejected_v29_round1_findings_corrected_in_active_candidate"
    ] is True


def test_v29_signed_semantics_and_fresh_reviewer_provenance_are_exact(modules):
    validator = modules["validator"]
    assert validator.AUTHORITY_STATUS == (
        "APPROVED_FOR_TUMOR_REF_V1_V6_AUDIT_LINEAGE_SPLIT_CHRONOLOGY_"
        "AND_HISTORICAL_CURRENT_RUNTIME_PROJECTION_SEPARATION_"
        "AND_BOUND_FINAL_DATASET_REPORT_RECOVERY_V30"
    )
    assert validator.HISTORICAL_V19_V21_AUTHORITY_STATUS == (
        "APPROVED_FOR_TRANSITION_ALIAS_METADATA_PLUS_SIZE_"
        "AND_TERMINAL_KEY_RECOVERY_ONLY"
    )
    assert validator.AUTHORITY_STATUS != validator.HISTORICAL_V19_V21_AUTHORITY_STATUS
    assert "exact metadata-plus-size relation validation" in (
        validator.AUTHORITY_PASS_SEMANTICS
    )
    assert "exact legacy eight-key stat relation" in (
        validator.AUTHORITY_PASS_SEMANTICS
    )
    assert "fresh v21 key for recovery-v30" in validator.AUTHORITY_PASS_SEMANTICS
    assert validator.AUTHORITY_CHECKS[
        "recovery_outputs_use_v30_paths_without_overwriting_v1_through_v29"
    ] is True
    assert "not cryptographic signatures" in validator.REVIEW_ATTRIBUTION_SEMANTICS
    assert validator.AUTHORITY_CHECKS[
        "reviewer_cryptographic_authorship_not_claimed"
    ] is True
    assert len(set(validator.EXPECTED_REVIEWER_AGENT_IDS.values())) == 3
    assert validator.EXPECTED_REVIEWER_AGENT_IDS == {
        "mendel": "019f929e-b3db-7c81-b483-ca9c4bcdf155",
        "nash": "019f929f-3d98-7390-bec5-294f49ae2c56",
        "external_claude_opus": None,
    }


def test_replayer_legacy_bindings_and_runtime_paths_are_separate(modules):
    replayer = modules["replayer"]
    assert replayer.expected_command() != replayer.authorized_replayer_command()
    assert replayer.verifier_command() != replayer.authorized_verifier_command()
    assert (
        replayer.downstream_continuation_command()
        != replayer.authorized_downstream_continuation_command()
    )
    assert replayer.PROMOTION_VERIFICATION_RECEIPT != (
        replayer.AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT
    )
    assert replayer.OUTPUT_RECEIPT != replayer.AUTHORIZED_OUTPUT_RECEIPT


def test_continuation_legacy_bindings_and_runtime_paths_are_separate(modules):
    validator = modules["validator"]
    continuation = modules["continuation"]
    assert continuation.canonical_command() != continuation.authorized_canonical_command()
    assert continuation.verifier_command() != continuation.authorized_verifier_command()
    assert continuation.replay_command() != continuation.authorized_replay_command()
    assert continuation.CONTINUATION_RECEIPT != continuation.AUTHORIZED_CONTINUATION_RECEIPT
    assert continuation.CONTINUATION_EXIT_ATTESTATION != (
        continuation.AUTHORIZED_CONTINUATION_EXIT_ATTESTATION
    )
    assert continuation.PROMOTION_VERIFICATION_RECEIPT != (
        continuation.AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT
    )
    assert continuation.RUNNER_GATE_RECEIPT != continuation.AUTHORIZED_RUNNER_GATE_RECEIPT

    authorization = json.loads(continuation.AUTHORIZATION.read_text(encoding="utf-8"))
    completion = json.loads(continuation.COMPLETION.read_text(encoding="utf-8"))
    assert authorization["commands"]["continuation_verify"] == (
        continuation.authorized_verifier_command()
    )
    assert authorization["commands"]["runner_gate_replay"] == (
        continuation.authorized_replay_command()
    )
    assert authorization["commands"]["downstream_continuation"] == (
        continuation.authorized_canonical_command()
    )
    assert authorization["downstream_output_absence"] == [
        str(path) for path in continuation.AUTHORIZED_DOWNSTREAM_OUTPUT_SLOTS
    ]
    signed_gate = authorization["continuation_gate"]
    assert signed_gate["verification_receipt"] == str(
        continuation.AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT
    )
    assert signed_gate["runner_gate_receipt"] == str(
        continuation.AUTHORIZED_RUNNER_GATE_RECEIPT
    )
    assert signed_gate["verifier"]["path"] == str(
        continuation.AUTHORIZED_CONTINUATION_VERIFIER
    )
    assert signed_gate["runner_gate_replay"]["path"] == str(
        continuation.AUTHORIZED_RUNNER_GATE_REPLAY
    )
    assert signed_gate["downstream_continuation"]["path"] == str(
        continuation.AUTHORIZED_CONTINUATION_RUNNER
    )
    assert completion["governance"]["continuation_verification_receipt"] == str(
        continuation.AUTHORIZED_PROMOTION_VERIFICATION_RECEIPT
    )
    assert completion["governance"]["runner_gate_receipt"] == str(
        continuation.AUTHORIZED_RUNNER_GATE_RECEIPT
    )

    active_self = basic_artifact(continuation.CONTINUATION_RUNNER)
    authorized_self = basic_artifact(continuation.AUTHORIZED_CONTINUATION_RUNNER)
    authorized_public = authorization["trusted_signing_keys"][
        "continuation_public_key"
    ]
    failed_v16_public = basic_artifact(
        continuation.FAILED_V16_CONTINUATION_PUBLIC_KEY
    )
    failed_v17_public = basic_artifact(
        continuation.FAILED_V17_CONTINUATION_PUBLIC_KEY
    )
    failed_v18_public = basic_artifact(
        continuation.FAILED_V18_CONTINUATION_PUBLIC_KEY
    )
    failed_v19_public = basic_artifact(
        continuation.FAILED_V19_CONTINUATION_PUBLIC_KEY
    )
    rejected_v20_public = basic_artifact(
        continuation.REJECTED_V20_CONTINUATION_PUBLIC_KEY
    )
    failed_v21_public = basic_artifact(
        continuation.FAILED_V21_CONTINUATION_PUBLIC_KEY
    )
    rejected_v22_public = basic_artifact(
        continuation.REJECTED_V22_CONTINUATION_PUBLIC_KEY
    )
    rejected_v23_public = basic_artifact(
        continuation.REJECTED_V23_CONTINUATION_PUBLIC_KEY
    )
    rejected_v24_public = basic_artifact(
        continuation.REJECTED_V24_CONTINUATION_PUBLIC_KEY
    )
    failed_v25_public = basic_artifact(
        continuation.FAILED_V25_CONTINUATION_PUBLIC_KEY
    )
    failed_v26_public = basic_artifact(
        continuation.FAILED_V26_CONTINUATION_PUBLIC_KEY
    )
    rejected_v27_public = basic_artifact(
        continuation.REJECTED_V27_CONTINUATION_PUBLIC_KEY
    )
    failed_v28_public = basic_artifact(
        continuation.FAILED_V28_CONTINUATION_PUBLIC_KEY
    )
    failed_v29_public = basic_artifact(
        continuation.FAILED_V29_CONTINUATION_PUBLIC_KEY
    )
    recovery_public = basic_artifact(continuation.RECOVERY_CONTINUATION_PUBLIC_KEY)
    assert authorization["reviewed_sources"]["downstream_continuation"] == authorized_self
    terminal_key_rotation, terminal_key_state = validator._validate_terminal_key_rotation(
        expected_recovery_private_mode="0o400"
    )
    assert set(terminal_key_rotation) == continuation.TERMINAL_KEY_ROTATION_KEYS
    assert set(terminal_key_state) == continuation.TERMINAL_KEY_STATE_KEYS
    assert terminal_key_state["failed_v28_public_key"] == failed_v28_public
    recovery_authority = {
        "authority_id": "20260724_tumor_ref_schema_recovery_v30",
        "runtime_role": "downstream_continuation",
        "runtime_source": active_self,
        "terminal_key_rotation": terminal_key_rotation,
        "terminal_key_state": terminal_key_state,
        "pass": True,
    }
    continuation.validate_signed_terminal_precommit_contract(
        authorization=authorization,
        completion=completion,
        active_self_record=active_self,
        authorized_self_record=authorized_self,
        authorized_public_record=authorized_public,
        failed_v16_public_record=failed_v16_public,
        failed_v17_public_record=failed_v17_public,
        failed_v18_public_record=failed_v18_public,
        failed_v19_public_record=failed_v19_public,
        rejected_v20_public_record=rejected_v20_public,
        failed_v21_public_record=failed_v21_public,
        rejected_v22_public_record=rejected_v22_public,
        rejected_v23_public_record=rejected_v23_public,
        rejected_v24_public_record=rejected_v24_public,
        failed_v25_public_record=failed_v25_public,
        failed_v26_public_record=failed_v26_public,
        rejected_v27_public_record=rejected_v27_public,
        failed_v28_public_record=failed_v28_public,
        failed_v29_public_record=failed_v29_public,
        recovery_public_record=recovery_public,
        canonical_record=completion["canonical_receipt"],
        source_record=completion["source_receipt"],
        recovery_authority=recovery_authority,
    )
    recovery_authority["runtime_source"] = authorized_self
    with pytest.raises(continuation.ContinuationError):
        continuation.validate_signed_terminal_precommit_contract(
            authorization=authorization,
            completion=completion,
            active_self_record=active_self,
            authorized_self_record=authorized_self,
            authorized_public_record=authorized_public,
            failed_v16_public_record=failed_v16_public,
            failed_v17_public_record=failed_v17_public,
            failed_v18_public_record=failed_v18_public,
            failed_v19_public_record=failed_v19_public,
            rejected_v20_public_record=rejected_v20_public,
                failed_v21_public_record=failed_v21_public,
                rejected_v22_public_record=rejected_v22_public,
                rejected_v23_public_record=rejected_v23_public,
                rejected_v24_public_record=rejected_v24_public,
                failed_v25_public_record=failed_v25_public,
                failed_v26_public_record=failed_v26_public,
                rejected_v27_public_record=rejected_v27_public,
                failed_v28_public_record=failed_v28_public,
                failed_v29_public_record=failed_v29_public,
                recovery_public_record=recovery_public,
            canonical_record=completion["canonical_receipt"],
            source_record=completion["source_receipt"],
            recovery_authority=recovery_authority,
        )


def test_replayer_and_continuation_guard_the_same_downstream_slots(modules):
    replayer = modules["replayer"]
    continuation = modules["continuation"]
    assert tuple(replayer.DOWNSTREAM_OUTPUT_SLOTS) == tuple(
        continuation.DOWNSTREAM_OUTPUT_SLOTS
    )


def test_v29_continuation_writer_lease_binds_active_public_key(modules):
    continuation = modules["continuation"]
    validator = modules["validator"]
    writer_key = validator.PUBLIC_KEY
    active_data = writer_key.read_bytes()
    active_stat = writer_key.stat()
    assert hashlib.sha256(active_data).hexdigest() == (
        validator.EXPECTED_PUBLIC_KEY_SHA256
    )
    assert active_stat.st_nlink == 1
    assert oct(active_stat.st_mode & 0o7777) == "0o444"

    try:
        descriptor = continuation.acquire_recovery_output_writer_lease()
    except continuation.ContinuationError as error:
        assert str(error) == (
            "Another cooperating recovery output writer holds the v30 lease"
        )
        contender = os.open(
            writer_key,
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW,
        )
        try:
            with pytest.raises(BlockingIOError):
                fcntl.flock(contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
        finally:
            os.close(contender)
        return
    try:
        observed = os.fstat(descriptor)
        assert observed.st_nlink == 1
        assert oct(observed.st_mode & 0o7777) == "0o444"
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def test_v29_continuation_preflight_matches_authorized_runner(modules):
    continuation = modules["continuation"]
    runner_source = continuation.COMPLETION_RUNNER.read_text(encoding="utf-8")
    assert continuation.COOCCURRENCE_PREFLIGHT.name == (
        "cooccurrence_task_contract_preflight.v9_command_parity_full_runtime.json"
    )
    assert continuation.COOCCURRENCE_PREFLIGHT.is_file()
    assert continuation.COOCCURRENCE_PREFLIGHT.name in runner_source
    assert hashlib.sha256(continuation.COOCCURRENCE_PREFLIGHT.read_bytes()).hexdigest() == (
        "e0826f0691514da04b84b7f45bd7e56a135aecff7d9a818bdb5795407b3a0591"
    )


def test_v2_builder_second_read_canary_reproduces_h1(monkeypatch):
    builder = load_module(
        "schema_recovery_v2_builder_h1_canary",
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v2.py",
    )
    real_read_bytes = Path.read_bytes
    calls = 0

    def counted_read_bytes(path):
        nonlocal calls
        if path == builder.VALIDATOR:
            calls += 1
        return real_read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", counted_read_bytes)
    builder.load_validator()
    assert calls == 2


def test_v29_builder_validator_single_open_uses_bound_bytes(modules, monkeypatch):
    builder = modules["builder"]
    real_open = builder.os.open
    real_read_bytes = Path.read_bytes
    opens = 0

    def counted_open(path, *args, **kwargs):
        nonlocal opens
        if Path(path) == builder.VALIDATOR:
            opens += 1
        return real_open(path, *args, **kwargs)

    def forbid_validator_path_read(path):
        if path == builder.VALIDATOR:
            raise AssertionError("validator Path.read_bytes() is forbidden")
        return real_read_bytes(path)

    monkeypatch.setattr(builder.os, "open", counted_open)
    monkeypatch.setattr(Path, "read_bytes", forbid_validator_path_read)
    validator, lease = builder.load_validator()
    try:
        assert opens == 1
        assert lease.sha256 == builder.EXPECTED_VALIDATOR_SHA256
        assert validator.AUTHORITY_ID == "20260724_tumor_ref_schema_recovery_v30"
    finally:
        lease.close()


def test_v29_copied_builder_fails_before_any_output(modules, tmp_path):
    validator = modules["validator"]
    builder = modules["builder"]
    copied = tmp_path / "copied_builder_v30.py"
    copied.write_bytes(builder.CEREMONY_BUILDER.read_bytes())
    copied.chmod(0o444)
    before = {
        str(path): os.path.lexists(path)
        for path in validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS
    }
    result = subprocess.run(
        [
            str(validator.PYTHON),
            "-I",
            "-B",
            str(copied),
            "--preflight",
        ],
        cwd=REPO_ROOT,
        env={
            "PATH": "/usr/bin:/bin",
            "HOME": "/tmp",
            "LANG": "C.UTF-8",
            "LC_ALL": "C.UTF-8",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "PYTHONDONTWRITEBYTECODE": "1",
        },
        capture_output=True,
        check=False,
        timeout=60,
    )
    after = {
        str(path): os.path.lexists(path)
        for path in validator.CEREMONY_FORBIDDEN_OUTPUT_SLOTS
    }
    assert result.returncode != 0
    assert b"exact canonical reviewed path" in result.stderr
    assert before == after


@pytest.mark.parametrize("mutation", ["mode", "content"])
def test_v29_bound_input_terminal_drift_fails_closed(modules, tmp_path, mutation):
    builder = modules["builder"]
    path = tmp_path / "bound-input.bin"
    path.write_bytes(b"AAAA")
    writer = os.open(path, os.O_RDWR)
    os.chmod(path, 0o444)
    lease = builder.BoundInput.open(path, expected_mode="0o444")
    try:
        if mutation == "mode":
            os.chmod(path, 0o644)
        else:
            os.pwrite(writer, b"B", 0)
            os.fsync(writer)
        with pytest.raises(RuntimeError, match="terminal recheck"):
            lease.recheck(expected_mode="0o444")
    finally:
        lease.close()
        os.close(writer)


def test_v29_private_retirement_recheck_requires_same_fd_and_forward_ctime(
    modules, tmp_path
):
    builder = modules["builder"]
    path = tmp_path / "one-time-key.pem"
    path.write_bytes(b"private-key-test-bytes")
    os.chmod(path, 0o400)
    lease = builder.BoundInput.open(path, expected_mode="0o400")
    try:
        builder.retire_private_key(lease)
        lease.recheck(expected_mode="0o0", allow_retirement=True)
        with pytest.raises(RuntimeError, match="terminal recheck"):
            lease.recheck(expected_mode="0o400")
    finally:
        lease.close()


def test_v29_writer_lease_is_exclusive(modules, monkeypatch, tmp_path):
    builder = modules["builder"]
    validator = modules["validator"]
    public_key = tmp_path / "ed25519_public.pem"
    public_key.write_bytes(b"v30 cooperating writer lease")
    public_key.chmod(0o444)
    monkeypatch.setattr(validator, "PUBLIC_KEY", public_key)
    monkeypatch.setattr(
        validator,
        "EXPECTED_PUBLIC_KEY_SHA256",
        hashlib.sha256(public_key.read_bytes()).hexdigest(),
    )
    first = builder.WriterLease.acquire(validator)
    try:
        with pytest.raises(RuntimeError, match="holds the v30 lease"):
            builder.WriterLease.acquire(validator)
    finally:
        first.close()
    second = builder.WriterLease.acquire(validator)
    second.close()


@pytest.mark.parametrize("mutated_parent_name", ["REVIEW_ROOT", "WORKSPACE_ROOT"])
def test_v29_terminal_publish_rechecks_every_absence_parent_token(
    modules, monkeypatch, tmp_path, mutated_parent_name
):
    builder = modules["builder"]
    validator = modules["validator"]
    bundle, _, leases = configure_temporary_absence_contract(
        builder, validator, monkeypatch, tmp_path
    )
    parent = leases[validator.RESULT_ROOT]
    stage_name = ".authority.v30.all-parent-token.staging"
    os.mkdir(stage_name, 0o700, dir_fd=parent.descriptor)
    stage_fd = os.open(
        stage_name,
        os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY,
        dir_fd=parent.descriptor,
    )
    authority_data = b'{"authority":true}\n'
    signature = b"S" * 64
    transaction_id = "b" * 32
    data = {
        "authority.json": authority_data,
        "authority.ed25519.sig": signature,
        "commit.json": builder.commit_data(
            validator, transaction_id, authority_data, signature
        ),
    }
    member_fds = {
        name: builder.write_member(stage_fd, name, value)
        for name, value in data.items()
    }
    mutation_watch = builder.DirectoryMutationWatch(leases)
    try:
        terminal_tokens = builder.terminal_stage_identity_recheck(
            validator,
            parent,
            stage_name,
            stage_fd,
            member_fds,
            transaction_id,
            authority_data,
            signature,
        )
        absence = builder.require_ceremony_absence(
            validator,
            leases,
            allowed_stage_name=stage_name,
            mutation_watch=mutation_watch,
        )
        mutated_parent = getattr(validator, mutated_parent_name)
        (mutated_parent / "late-forbidden-output.json").write_bytes(b"late")
        with pytest.raises(
            RuntimeError,
            match="absence-parent token drift|mutation event observed",
        ):
            builder.atomic_publish_after_terminal_contract(
                parent,
                leases,
                mutation_watch,
                absence["directory_generation_tokens"],
                stage_name,
                stage_fd,
                member_fds,
                bundle.name,
                terminal_tokens,
                lambda: None,
            )
        assert not bundle.exists()
    finally:
        mutation_watch.close()
        for descriptor in member_fds.values():
            os.close(descriptor)
        os.close(stage_fd)
        for lease in leases.values():
            lease.close()


@pytest.mark.parametrize("occupied", [False, True])
def test_v29_rename_noreplace_exposes_only_complete_directory(
    modules, tmp_path, occupied
):
    builder = modules["builder"]
    parent = builder.DirectoryLease(tmp_path)
    stage = ".authority.staging"
    final = "authority.bundle"
    stage_path = tmp_path / stage
    stage_path.mkdir(mode=0o700)
    for name in builder.MEMBERS:
        (stage_path / name).write_bytes(name.encode("ascii"))
    if occupied:
        (tmp_path / final).mkdir()
    try:
        if occupied:
            with pytest.raises(OSError):
                builder.rename_noreplace(parent.descriptor, stage, final)
            assert stage_path.is_dir()
            assert (tmp_path / final).is_dir()
        else:
            builder.rename_noreplace(parent.descriptor, stage, final)
            assert not stage_path.exists()
            assert set(path.name for path in (tmp_path / final).iterdir()) == set(
                builder.MEMBERS
            )
    finally:
        parent.close()


def test_v29_atomic_commit_binds_authority_signature_and_transaction(modules):
    builder = modules["builder"]
    validator = modules["validator"]
    transaction_id = "a" * 32
    authority_data = b'{"pass":true}\n'
    signature = b"S" * 64
    commit = builder.commit_data(
        validator, transaction_id, authority_data, signature
    )
    authority_record = {
        "path": str(validator.AUTHORITY),
        "size_bytes": len(authority_data),
        "sha256": hashlib.sha256(authority_data).hexdigest(),
        "mode": "0o444",
    }
    signature_record = {
        "path": str(validator.AUTHORITY_SIGNATURE),
        "size_bytes": len(signature),
        "sha256": hashlib.sha256(signature).hexdigest(),
        "mode": "0o444",
    }
    private = {"path": str(validator.PRIVATE_KEY), "mode": "0o0"}
    validator._validate_bundle_commit(
        commit, authority_record, signature_record, private
    )
    drift = json.loads(commit)
    drift["authority"]["sha256"] = "0" * 64
    with pytest.raises(validator.RecoveryAuthorityError):
        validator._validate_bundle_commit(
            json.dumps(drift).encode("ascii"),
            authority_record,
            signature_record,
            private,
        )


def test_v3_full_builder_reopens_validator_canary(monkeypatch):
    builder = load_module(
        "schema_recovery_v3_builder_second_open_canary",
        AUDIT_ROOT / "build_tumor_ref_schema_recovery_authority_v3.py",
    )
    real_open = builder.os.open
    opens = 0

    def counted_open(path, *args, **kwargs):
        nonlocal opens
        if Path(path) == builder.VALIDATOR:
            opens += 1
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr(builder.os, "open", counted_open)
    validator, lease = builder.load_validator()
    try:
        with pytest.raises(FileNotFoundError):
            builder.build_payloads(validator, "a" * 32)
        assert opens == 2
    finally:
        lease.close()


def test_v29_full_builder_reuses_bootstrap_validator_record(modules, monkeypatch):
    builder = modules["builder"]
    validator = modules["validator"]
    builder_source = builder.BoundInput.open(
        builder.CEREMONY_BUILDER, expected_mode="0o444"
    )
    bootstrap = builder.BoundInput.open(builder.VALIDATOR, expected_mode="0o444")
    real_records = validator._records
    observed_calls = 0

    def guarded_records(paths):
        nonlocal observed_calls
        observed_calls += 1
        assert builder.VALIDATOR not in set(paths.values())
        return real_records(paths)

    monkeypatch.setattr(validator, "_records", guarded_records)
    review_state = {
        role: os.path.lexists(path) for role, path in validator.REVIEW_PATHS.items()
    }
    assert all(review_state.values()) or not any(review_state.values())
    fd_lease_count = len(validator._FD_LEASES)
    directory_lease_count = len(validator._DIRECTORY_LEASES)
    try:
        expected_probe_sources = {
            "authority_validator": bootstrap.record(),
            **validator._records(
                {
                    role: path
                    for role, path in validator.RECOVERY_SOURCE_PATHS.items()
                    if role != "authority_validator"
                }
            ),
        }
        builder_binding = {
            "source": builder_source.record(),
            "canonical_file": str(builder.CEREMONY_BUILDER),
            "sys_argv0": str(builder.CEREMONY_BUILDER),
            "proc_cmdline_contains_exact_canonical_script_once": True,
            "bound_before_validator_load": True,
            "same_path_inode_as_bound_source": True,
            "threat_model": (
                "trusted_same_uid_account_no_malicious_runtime_code_injection"
            ),
            "pass": True,
        }
        probe_execution = {
            "probe_source": expected_probe_sources["readonly_probe"],
            "python_runtime": validator.EXPECTED_PYTHON_RUNTIME,
            "launch": "bound_python_fd_exec_a_canonical_alias_bound_probe_source_fd",
            "regression_test_source": {
                "source": expected_probe_sources["regression_tests"],
                "python_runtime": validator.EXPECTED_PYTHON_RUNTIME,
                "execution": (
                    "pytest_preloaded_module_compiled_from_bound_source_fd_"
                    "via_bound_python_fd"
                ),
                "canonical_python_argv0": str(validator.PYTHON),
                "pass": True,
            },
            "pass": True,
        }
        try:
            reviews, authority_data, public_fd = builder.build_payloads(
                validator,
                builder_source,
                bootstrap,
                "b" * 32,
                builder_binding=builder_binding,
                probe_execution=probe_execution,
                expected_probe_sources=expected_probe_sources,
            )
        except FileNotFoundError as error:
            assert not any(review_state.values())
            require_missing_review_path(error, validator.REVIEW_PATHS.values())
        else:
            assert all(review_state.values())
            assert set(reviews) == set(validator.REVIEW_PATHS)
            assert authority_data
            assert public_fd >= 0
        assert observed_calls > 0
    finally:
        for _, descriptor, _, _ in validator._FD_LEASES[fd_lease_count:]:
            os.close(descriptor)
        del validator._FD_LEASES[fd_lease_count:]
        for _, descriptor, _ in validator._DIRECTORY_LEASES[directory_lease_count:]:
            os.close(descriptor)
        del validator._DIRECTORY_LEASES[directory_lease_count:]
        bootstrap.close()
        builder_source.close()


def test_v29_nonreview_FileNotFound_is_not_classified_as_review_absence(modules):
    validator = modules["validator"]
    error = FileNotFoundError(2, "missing", str(validator.REJECTED_V7_EVIDENCE))
    with pytest.raises(AssertionError):
        require_missing_review_path(error, validator.REVIEW_PATHS.values())


def test_v29_sealed_memfd_uses_linux_flag_fallbacks(modules):
    builder = modules["builder"]
    assert builder.LINUX_MFD_CLOEXEC == getattr(os, "MFD_CLOEXEC", 0x0001)
    assert builder.LINUX_MFD_ALLOW_SEALING == getattr(
        os, "MFD_ALLOW_SEALING", 0x0002
    )
    descriptor = builder.sealed_memfd("schema_recovery_flag_test", b"sealed")
    try:
        assert os.pread(descriptor, 6, 0) == b"sealed"
        expected_seals = (
            builder.LINUX_F_SEAL_SEAL
            | builder.LINUX_F_SEAL_SHRINK
            | builder.LINUX_F_SEAL_GROW
            | builder.LINUX_F_SEAL_WRITE
        )
        assert (
            fcntl.fcntl(descriptor, builder.LINUX_F_GET_SEALS) == expected_seals
        )
    finally:
        os.close(descriptor)


def test_v29_success_witness_commit_is_real_no_return_link(
    modules, monkeypatch, tmp_path
):
    continuation = modules["continuation"]
    witness = tmp_path / "terminal-success-witness.json"
    monkeypatch.setattr(continuation, "CONTINUATION_SUCCESS_WITNESS", witness)
    payload = {"schema_name": "synthetic.success", "pass": True}
    child = os.fork()
    if child == 0:
        try:
            continuation.commit_terminal_success_witness(payload, lambda: None)
        except BaseException:
            os._exit(93)
    _, status = os.waitpid(child, 0)
    assert os.waitstatus_to_exitcode(status) == 0
    assert json.loads(witness.read_text(encoding="ascii")) == payload
    observed = witness.stat()
    assert observed.st_nlink == 1
    assert oct(observed.st_mode & 0o7777) == "0o444"


def test_v29_final_verifier_requires_released_writer_lease_and_witness(
    modules, monkeypatch, tmp_path
):
    continuation = modules["continuation"]
    public_key = tmp_path / "ed25519_public.pem"
    public_key.write_bytes(b"v29 continuation terminal lease")
    public_key.chmod(0o444)
    monkeypatch.setattr(
        continuation, "RECOVERY_OUTPUT_LEASE_PUBLIC_KEY", public_key
    )
    monkeypatch.setattr(
        continuation,
        "EXPECTED_RECOVERY_OUTPUT_LEASE_PUBLIC_KEY_SHA256",
        hashlib.sha256(public_key.read_bytes()).hexdigest(),
    )

    writer = continuation.acquire_recovery_output_writer_lease()
    try:
        with pytest.raises(
            continuation.ContinuationError,
            match="writer is still active",
        ):
            continuation.acquire_recovery_output_verifier_lease()
    finally:
        os.close(writer)

    verifier = continuation.acquire_recovery_output_verifier_lease()
    try:
        with pytest.raises(
            continuation.ContinuationError,
            match="holds the v30 lease",
        ):
            continuation.acquire_recovery_output_writer_lease()
    finally:
        os.close(verifier)

    retained = os.open(public_key, os.O_RDONLY | os.O_CLOEXEC)
    monkeypatch.setattr(continuation, "require_clean_runtime", lambda **kwargs: None)
    monkeypatch.setattr(
        continuation, "acquire_recovery_output_verifier_lease", lambda: retained
    )
    monkeypatch.setattr(continuation, "require_incident_absent", lambda: None)
    monkeypatch.setattr(continuation.os.path, "lexists", lambda _path: False)
    try:
        with pytest.raises(
            continuation.ContinuationError,
            match="success witness is absent",
        ):
            continuation.verify_signed_terminal(require_success_witness=True)
    finally:
        continuation._VERIFIER_LEASE_FDS.remove(retained)
        os.close(retained)


def test_v29_authority_retires_persistent_key_before_real_memfd_signing(
    modules, tmp_path
):
    builder = modules["builder"]
    validator = modules["validator"]
    private_path, public_path = generate_ed25519_keypair(
        validator.OPENSSL, tmp_path, "authority"
    )
    private = builder.BoundInput.open(private_path, expected_mode="0o400")
    public_fd = os.open(public_path, os.O_RDONLY | os.O_CLOEXEC)
    signing_key_fd = -1
    try:
        signing_key_fd = builder.stage_signing_key(private)
        staged = os.fstat(signing_key_fd)
        assert staged.st_nlink == 0
        assert staged.st_mode & 0o7777 == 0o400
        builder.retire_private_key(private)
        assert private_path.stat().st_mode & 0o7777 == 0
        signature = builder.sign_and_verify(
            validator, b'{"authority":true}\n', signing_key_fd, public_fd
        )
        assert len(signature) == 64
    finally:
        if signing_key_fd >= 0:
            os.close(signing_key_fd)
        private.close()
        os.close(public_fd)


@pytest.mark.parametrize("exit_type", [KeyboardInterrupt, SystemExit])
def test_v29_publish_retires_key_on_post_sign_BaseException(
    modules, monkeypatch, tmp_path, exit_type
):
    builder = modules["builder"]
    validator = modules["validator"]
    private_path = tmp_path / "private.pem"
    private_path.write_bytes(b"unused one-time key")
    private_path.chmod(0o400)
    private = builder.BoundInput.open(private_path, expected_mode="0o400")
    parent = builder.DirectoryLease(tmp_path)
    state = {
        "private": private,
        "writer": NoopWriterLease(),
        "parent": parent,
        "absence_parents": {tmp_path: parent},
        "authority_data": b"authority",
        "public_fd": -1,
        "transaction_id": "d" * 32,
        "reviews": {},
        "probe": {"pass": True},
    }

    def interrupted_sign(*args):
        raise exit_type()

    monkeypatch.setattr(builder, "preflight", lambda *args: state)
    monkeypatch.setattr(builder, "terminal_recheck", lambda *args, **kwargs: None)
    monkeypatch.setattr(builder, "sign_and_verify", interrupted_sign)
    with pytest.raises(exit_type):
        builder.publish(validator, None, None, {})
    assert private_path.stat().st_mode & 0o7777 == 0


def test_v29_partial_retirement_is_retried_after_BaseException(
    modules, monkeypatch, tmp_path
):
    builder = modules["builder"]
    validator = modules["validator"]
    private_path = tmp_path / "private.pem"
    private_path.write_bytes(b"unused one-time key")
    private_path.chmod(0o400)
    private = builder.BoundInput.open(private_path, expected_mode="0o400")
    parent = builder.DirectoryLease(tmp_path)
    state = {
        "private": private,
        "writer": NoopWriterLease(),
        "parent": parent,
        "absence_parents": {tmp_path: parent},
        "authority_data": b"authority",
        "public_fd": -1,
        "transaction_id": "f" * 32,
        "reviews": {},
        "probe": {"pass": True},
    }
    real_retire = builder.retire_private_key
    retirement_calls = 0

    def forbidden_sign(*args):
        pytest.fail("signing must not run after interrupted persistent-key retirement")

    def interrupted_retirement(bound):
        nonlocal retirement_calls
        retirement_calls += 1
        if retirement_calls == 1:
            os.fchmod(bound.descriptor, 0)
            raise KeyboardInterrupt()
        real_retire(bound)

    monkeypatch.setattr(builder, "preflight", lambda *args: state)
    monkeypatch.setattr(builder, "terminal_recheck", lambda *args, **kwargs: None)
    monkeypatch.setattr(builder, "sign_and_verify", forbidden_sign)
    monkeypatch.setattr(builder, "stage_bundle", lambda *args: ("unused", -1, {}))
    monkeypatch.setattr(builder, "retire_private_key", interrupted_retirement)
    with pytest.raises(KeyboardInterrupt):
        builder.publish(validator, None, None, {})
    assert retirement_calls == 2
    assert private_path.stat().st_mode & 0o7777 == 0


def test_v29_publish_never_calls_runtime_validator_after_atomic_publish(modules):
    builder = modules["builder"]
    names = set(builder.publish.__code__.co_names)
    assert "validate_recovery_authority" not in names
    assert "publish_staged_bundle" in names
    assert "published_record" in names


def test_v29_forbidden_slot_before_sign_aborts_without_burning_key(
    modules, monkeypatch, tmp_path
):
    builder = modules["builder"]
    validator = modules["validator"]
    bundle, forbidden, absence_parents = configure_temporary_absence_contract(
        builder, validator, monkeypatch, tmp_path
    )
    private_path = tmp_path / "private-pre-sign.pem"
    private_path.write_bytes(b"unused one-time key")
    private_path.chmod(0o400)
    monkeypatch.setattr(validator, "PRIVATE_KEY", private_path)
    private = builder.BoundInput.open(private_path, expected_mode="0o400")
    parent = absence_parents[validator.RESULT_ROOT]
    state = {
        "private": private,
        "writer": NoopWriterLease(),
        "parent": parent,
        "absence_parents": absence_parents,
        "authority_data": b'{"authority":true}\n',
        "public_fd": -1,
        "transaction_id": "7" * 32,
        "reviews": {},
        "probe": {"pass": True},
    }

    class StableValidatorLease:
        def recheck(self, **kwargs):
            return None

    sign_calls = 0

    def preflight_with_injection(*args):
        forbidden.write_bytes(b"appeared after build_payloads")
        return state

    def forbidden_sign(*args):
        nonlocal sign_calls
        sign_calls += 1
        pytest.fail("signing must not follow a terminal absence failure")

    monkeypatch.setattr(validator, "_require_leases", lambda **kwargs: None)
    monkeypatch.setattr(builder, "preflight", preflight_with_injection)
    monkeypatch.setattr(builder, "sign_and_verify", forbidden_sign)
    with pytest.raises(RuntimeError, match="forbidden namespace"):
        builder.publish(
            validator, StableValidatorLease(), StableValidatorLease(), {}
        )
    assert sign_calls == 0
    assert private_path.stat().st_mode & 0o7777 == 0o400
    assert not bundle.exists()


def test_v29_forbidden_slot_after_external_verification_blocks_rename_and_retires_key(
    modules, monkeypatch, tmp_path
):
    builder = modules["builder"]
    validator = modules["validator"]
    bundle, forbidden, absence_parents = configure_temporary_absence_contract(
        builder, validator, monkeypatch, tmp_path
    )
    private_path = tmp_path / "private-post-sign.pem"
    private_path.write_bytes(b"unused one-time key")
    private_path.chmod(0o400)
    monkeypatch.setattr(validator, "PRIVATE_KEY", private_path)
    private = builder.BoundInput.open(private_path, expected_mode="0o400")
    parent = absence_parents[validator.RESULT_ROOT]
    state = {
        "private": private,
        "writer": NoopWriterLease(),
        "parent": parent,
        "absence_parents": absence_parents,
        "authority_data": b'{"authority":true}\n',
        "public_fd": -1,
        "transaction_id": "8" * 32,
        "reviews": {},
        "probe": {"pass": True},
    }

    class StableValidatorLease:
        def recheck(self, **kwargs):
            return None

    publish_calls = 0

    def successful_sign(*args):
        return b"S" * 64

    def inject_after_external_verification(*args, **kwargs):
        forbidden.write_bytes(b"appeared after external verification")

    def forbidden_publish(*args):
        nonlocal publish_calls
        publish_calls += 1
        pytest.fail("rename must not follow a terminal absence failure")

    monkeypatch.setattr(validator, "_require_leases", lambda **kwargs: None)
    monkeypatch.setattr(builder, "preflight", lambda *args: state)
    monkeypatch.setattr(builder, "sign_and_verify", successful_sign)
    monkeypatch.setattr(
        builder, "verify_signature_fds", inject_after_external_verification
    )
    monkeypatch.setattr(
        builder, "atomic_publish_after_terminal_contract", forbidden_publish
    )
    with pytest.raises(RuntimeError, match="forbidden namespace"):
        builder.publish(
            validator, StableValidatorLease(), StableValidatorLease(), {}
        )
    assert publish_calls == 0
    assert private_path.stat().st_mode & 0o7777 == 0
    assert not bundle.exists()


def test_v29_successful_publish_uses_retained_fds_without_runtime_validator(
    modules, monkeypatch, tmp_path
):
    builder = modules["builder"]
    validator = modules["validator"]
    bundle = tmp_path / "authority.v30.bundle"
    private_path = tmp_path / "private.pem"
    private_path.write_bytes(b"unused one-time key")
    private_path.chmod(0o400)
    private = builder.BoundInput.open(private_path, expected_mode="0o400")
    parent = builder.DirectoryLease(tmp_path)
    state = {
        "private": private,
        "writer": NoopWriterLease(),
        "parent": parent,
        "absence_parents": {tmp_path: parent},
        "authority_data": b'{"authority":true}\n',
        "public_fd": -1,
        "transaction_id": "9" * 32,
        "reviews": {"mendel": {}, "nash": {}, "external_claude_opus": {}},
        "probe": {"pass": True},
    }
    runtime_validator_calls = 0
    signal_mask_before = signal.pthread_sigmask(signal.SIG_BLOCK, set())
    final_recheck_signal_mask = None
    real_terminal_stage_identity_recheck = builder.terminal_stage_identity_recheck

    def successful_sign(*args):
        return b"S" * 64

    def forbidden_runtime_validator(*args, **kwargs):
        nonlocal runtime_validator_calls
        runtime_validator_calls += 1
        pytest.fail("published bundle must not be reopened by the runtime validator")

    def observe_terminal_stage_identity_recheck(*args, **kwargs):
        nonlocal final_recheck_signal_mask
        if args[2] == bundle.name:
            final_recheck_signal_mask = signal.pthread_sigmask(signal.SIG_BLOCK, set())
        return real_terminal_stage_identity_recheck(*args, **kwargs)

    def stable_ceremony_absence(module, absence_parents, **kwargs):
        return {
            "directory_generation_tokens": {
                str(path): builder.full_identity(os.fstat(lease.descriptor))
                for path, lease in absence_parents.items()
            }
        }

    monkeypatch.setattr(validator, "AUTHORITY_BUNDLE", bundle)
    monkeypatch.setattr(validator, "AUTHORITY", bundle / "authority.json")
    monkeypatch.setattr(
        validator, "AUTHORITY_SIGNATURE", bundle / "authority.ed25519.sig"
    )
    monkeypatch.setattr(validator, "AUTHORITY_COMMIT", bundle / "commit.json")
    monkeypatch.setattr(validator, "PRIVATE_KEY", private_path)
    monkeypatch.setattr(builder, "preflight", lambda *args: state)
    monkeypatch.setattr(builder, "terminal_recheck", lambda *args, **kwargs: None)
    monkeypatch.setattr(builder, "require_ceremony_absence", stable_ceremony_absence)
    monkeypatch.setattr(builder, "sign_and_verify", successful_sign)
    monkeypatch.setattr(builder, "verify_signature_fds", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        builder,
        "terminal_stage_identity_recheck",
        observe_terminal_stage_identity_recheck,
    )
    monkeypatch.setattr(
        validator, "validate_recovery_authority", forbidden_runtime_validator
    )

    result = builder.publish(validator, None, None, {})
    signal_mask_after = signal.pthread_sigmask(signal.SIG_BLOCK, set())

    assert result["pass"] is True
    assert result["publication_validation"] == {
        "basis": "retained_stage_directory_and_member_fds",
        "validator_path_reopened": False,
        "terminal_stage_recheck_after_external_verification": True,
        "pass": True,
    }
    assert runtime_validator_calls == 0
    assert final_recheck_signal_mask is not None
    assert signal.SIGTERM in final_recheck_signal_mask
    assert signal_mask_after == signal_mask_before
    assert private_path.stat().st_mode & 0o7777 == 0
    assert bundle.is_dir()
    assert {path.name for path in bundle.iterdir()} == set(builder.MEMBERS)


@pytest.mark.parametrize("mutation", ["member", "directory"])
def test_v29_pre_rename_recheck_rejects_staged_drift(
    modules, monkeypatch, tmp_path, mutation
):
    builder = modules["builder"]
    validator = modules["validator"]
    parent = builder.DirectoryLease(tmp_path)
    stage_name = ".authority.v30.staging"
    os.mkdir(stage_name, 0o700, dir_fd=parent.descriptor)
    stage_fd = os.open(
        stage_name,
        os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY,
        dir_fd=parent.descriptor,
    )
    authority_data = b'{"authority":true}\n'
    signature = b"S" * 64
    transaction_id = "c" * 32
    data = {
        "authority.json": authority_data,
        "authority.ed25519.sig": signature,
        "commit.json": builder.commit_data(
            validator, transaction_id, authority_data, signature
        ),
    }
    member_fds = {
        name: builder.write_member(stage_fd, name, value)
        for name, value in data.items()
    }
    monkeypatch.setattr(builder, "verify_signature_fds", lambda *args, **kwargs: None)
    try:
        if mutation == "member":
            os.pwrite(member_fds["authority.json"], b"X", 0)
            os.fsync(member_fds["authority.json"])
        else:
            os.rename(
                stage_name,
                f"{stage_name}.moved",
                src_dir_fd=parent.descriptor,
                dst_dir_fd=parent.descriptor,
            )
            os.mkdir(stage_name, 0o700, dir_fd=parent.descriptor)
        with pytest.raises(RuntimeError, match="staging|staged member"):
            builder.recheck_staged_bundle(
                validator,
                parent,
                stage_name,
                stage_fd,
                member_fds,
                transaction_id,
                authority_data,
                signature,
                -1,
            )
    finally:
        for descriptor in member_fds.values():
            os.close(descriptor)
        os.close(stage_fd)
        parent.close()


@pytest.mark.parametrize("mutation", ["member", "directory"])
def test_v29_terminal_recheck_rejects_drift_after_external_verification(
    modules, monkeypatch, tmp_path, mutation
):
    builder = modules["builder"]
    validator = modules["validator"]
    parent = builder.DirectoryLease(tmp_path)
    stage_name = ".authority.v30.post-verification.staging"
    os.mkdir(stage_name, 0o700, dir_fd=parent.descriptor)
    stage_fd = os.open(
        stage_name,
        os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY,
        dir_fd=parent.descriptor,
    )
    authority_data = b'{"authority":true}\n'
    signature = b"S" * 64
    transaction_id = "e" * 32
    data = {
        "authority.json": authority_data,
        "authority.ed25519.sig": signature,
        "commit.json": builder.commit_data(
            validator, transaction_id, authority_data, signature
        ),
    }
    member_fds = {
        name: builder.write_member(stage_fd, name, value)
        for name, value in data.items()
    }

    def mutate_after_initial_identity_check(*args, **kwargs):
        if mutation == "member":
            os.pwrite(member_fds["authority.json"], b"X", 0)
            os.fsync(member_fds["authority.json"])
        else:
            os.rename(
                stage_name,
                f"{stage_name}.moved",
                src_dir_fd=parent.descriptor,
                dst_dir_fd=parent.descriptor,
            )
            os.mkdir(stage_name, 0o700, dir_fd=parent.descriptor)

    monkeypatch.setattr(builder, "verify_signature_fds", mutate_after_initial_identity_check)
    monkeypatch.setattr(builder, "require_ceremony_absence", lambda *args, **kwargs: {})
    monkeypatch.setattr(
        builder,
        "atomic_publish_after_terminal_contract",
        lambda *args: pytest.fail("rename must not follow terminal identity drift"),
    )
    try:
        with pytest.raises(RuntimeError, match="staging|staged member"):
            builder.publish_staged_bundle(
                validator,
                parent,
                {},
                stage_name,
                stage_fd,
                member_fds,
                transaction_id,
                authority_data,
                signature,
                -1,
            )
    finally:
        for descriptor in member_fds.values():
            os.close(descriptor)
        os.close(stage_fd)
        parent.close()


@pytest.mark.parametrize(
    ("primary_fd", "qa_fd"),
    [(2, 9), (9, 2), (9, 9), ("9", 10)],
)
def test_v29_runtime_fd_bridge_rejects_invalid_descriptor_pairs(
    modules, primary_fd, qa_fd
):
    continuation = modules["continuation"]
    with pytest.raises(
        continuation.ContinuationError,
        match="Runtime Python descriptor bridge FD drift",
    ):
        continuation.render_runtime_fd_bridge(primary_fd, qa_fd)


def test_v29_composed_runner_binds_wrappers_canonical_argv0_and_fd_paths(modules):
    continuation = modules["continuation"]
    finalizer = modules["finalizer"]
    probe = modules["probe"]
    primary_fd = os.open(continuation.PRIMARY_PYTHON_RUNTIME, os.O_RDONLY)
    qa_fd = os.open(continuation.QA_PYTHON_RUNTIME, os.O_RDONLY)
    try:
        script, contract = continuation.compose_downstream_script(
            continuation.COMPLETION_RUNNER.read_bytes(),
            primary_python_fd=primary_fd,
            qa_python_fd=qa_fd,
            primary_python_record={"path": str(continuation.PRIMARY_PYTHON_RUNTIME)},
            qa_python_record={"path": str(continuation.QA_PYTHON_RUNTIME)},
        )
    finally:
        os.close(primary_fd)
        os.close(qa_fd)

    wrappers = contract["runtime_fd_canonical_argv0_wrappers"]
    assert wrappers["primary_python"]["wrapper_function"] == (
        continuation.PRIMARY_PYTHON_WRAPPER
    )
    assert wrappers["primary_python"]["canonical_argv0"] == str(
        continuation.PYTHON
    )
    assert wrappers["primary_python"]["bound_executable_fd_path"] == (
        f"/proc/self/fd/{primary_fd}"
    )
    assert wrappers["qa_python"]["wrapper_function"] == continuation.QA_PYTHON_WRAPPER
    assert wrappers["qa_python"]["canonical_argv0"] == str(continuation.QA_PYTHON)
    assert wrappers["qa_python"]["bound_executable_fd_path"] == (
        f"/proc/self/fd/{qa_fd}"
    )
    assert b'readonly PYTHON="run_bound_primary_python"' in script
    assert b'readonly QA_PYTHON="run_bound_qa_python"' in script
    assert contract["rendered_runtime_fd_bridge"]["placeholder_count"] == 0
    assert contract["rendered_runtime_fd_bridge"]["execution_semantics"] == (
        "bound_fd_bytes_with_canonical_argv0"
    )
    expected_rewrites = {
        "${SCRIPT_ROOT}/build_all_ssnv_final_report_dataset.py": str(
            continuation.FINAL_DATASET_BUILDER
        ),
        "${SCRIPT_ROOT}/build_all_ssnv_report_artifact.py": str(
            continuation.REPORT_BUILDER
        ),
        "${SCRIPT_ROOT}/finalize_task_b_result_release.py": str(
            continuation.FINAL_RELEASE_FINALIZER
        ),
    }
    assert {
        source: record["executed_path"]
        for source, record in contract["recovery_script_rewrites"].items()
    } == expected_rewrites
    for source, executed in expected_rewrites.items():
        assert source.encode("ascii") not in script
        assert executed.encode("ascii") in script
    expected_key_rewrites = {
        role: (source, replacement)
        for source, replacement, role in continuation.RUNNER_CURRENT_RELEASE_KEY_REBASE
    }
    assert set(contract["current_release_key_rewrites"]) == set(
        expected_key_rewrites
    )
    for role, (source, replacement) in expected_key_rewrites.items():
        record = contract["current_release_key_rewrites"][role]
        assert record == {
            "source": source,
            "executed_value": replacement,
            "source_occurrence_count": 1,
            "executed_occurrence_count": 1,
        }
        assert source.encode("ascii") not in script
        assert replacement.encode("ascii") in script
    assert finalizer.RESULT_PUBLIC_KEY == continuation.FINAL_RELEASE_PUBLIC_KEY
    assert finalizer.RESULT_PRIVATE_KEY == continuation.FINAL_RELEASE_PRIVATE_KEY
    assert finalizer.RESULT_PUBLIC_KEY_SHA256 == (
        continuation.EXPECTED_FINAL_RELEASE_PUBLIC_KEY_SHA256
    )
    assert finalizer.REPORT_PUBLIC_KEY == continuation.REPORT_RELEASE_PUBLIC_KEY
    assert finalizer.REPORT_PRIVATE_KEY == continuation.REPORT_RELEASE_PRIVATE_KEY
    assert finalizer.REPORT_PUBLIC_KEY_SHA256 == (
        continuation.EXPECTED_REPORT_RELEASE_PUBLIC_KEY_SHA256
    )
    assert finalizer.RELEASE_RECEIPT == continuation.FINAL_DATASET_RELEASE_RECEIPT
    assert finalizer.RELEASE_SIGNATURE == continuation.FINAL_DATASET_RELEASE_SIGNATURE
    assert finalizer.REPORT_RELEASE_RECEIPT == (
        continuation.FINAL_REPORT_RELEASE_RECEIPT
    )
    assert finalizer.REPORT_RELEASE_SIGNATURE == (
        continuation.FINAL_REPORT_RELEASE_SIGNATURE
    )
    assert probe.EXPECTED_SOURCES["recovery_result_report_finalizer"][:3] == (
        continuation.FINAL_RELEASE_FINALIZER,
        33_422,
        "b9645cf4f57653f078357421592cef19a168db8389ee68ae90898b9c8c63d318",
    )
    executed_prefix_size = contract["executed_runtime_bound_prefix"]["size_bytes"]
    prefix_result = subprocess.run(
        [str(continuation.BASH), "-s"],
        input=script[:executed_prefix_size],
        env=continuation.EXPECTED_ENVIRONMENT,
        capture_output=True,
        check=False,
    )
    assert prefix_result.returncode == 0, prefix_result.stderr.decode(
        "utf-8", errors="replace"
    )
    assert prefix_result.stderr == b""
    assert continuation.FINAL_DATASET.name.encode("ascii") in script
    assert continuation.FINAL_REPORT.name.encode("ascii") in script


def test_v29_bound_python_wrappers_execute_fd_bytes_with_canonical_argv0(modules):
    continuation = modules["continuation"]
    primary_fd = os.open(continuation.PRIMARY_PYTHON_RUNTIME, os.O_RDONLY)
    qa_fd = os.open(continuation.QA_PYTHON_RUNTIME, os.O_RDONLY)
    probe = (
        "import json,os,sys;from pathlib import Path;fd=int(sys.argv[1]);"
        "print(json.dumps({\"sys_executable\":sys.executable,"
        "\"cmdline_argv0\":Path(\"/proc/self/cmdline\").read_bytes().split(b\"\\0\",1)[0].decode(),"
        "\"fd_device\":os.fstat(fd).st_dev,\"fd_inode\":os.fstat(fd).st_ino},sort_keys=True))"
    )
    rendered = continuation.render_runtime_fd_bridge(primary_fd, qa_fd)
    script = (
        b"set -Eeuo pipefail\n"
        + rendered
        + (
            f'"{continuation.PRIMARY_PYTHON_WRAPPER}" -I -B -c '
            f"{shlex.quote(probe)} {primary_fd}\n"
            f'"{continuation.QA_PYTHON_WRAPPER}" -I -B -c '
            f"{shlex.quote(probe)} {qa_fd}\n"
        ).encode("ascii")
    )
    try:
        result = subprocess.run(
            [str(continuation.BASH), "-s"],
            input=script,
            env=continuation.EXPECTED_ENVIRONMENT,
            pass_fds=(primary_fd, qa_fd),
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr.decode("utf-8", errors="replace")
        rows = [json.loads(line) for line in result.stdout.decode("ascii").splitlines()]
        assert len(rows) == 2
        expected = (
            (continuation.PYTHON, continuation.PRIMARY_PYTHON_RUNTIME, primary_fd),
            (continuation.QA_PYTHON, continuation.QA_PYTHON_RUNTIME, qa_fd),
        )
        for row, (canonical, runtime, descriptor) in zip(rows, expected, strict=True):
            runtime_stat = os.stat(runtime, follow_symlinks=False)
            assert row == {
                "sys_executable": str(canonical),
                "cmdline_argv0": str(canonical),
                "fd_device": runtime_stat.st_dev,
                "fd_inode": runtime_stat.st_ino,
            }
            assert os.fstat(descriptor).st_ino == runtime_stat.st_ino
    finally:
        os.close(primary_fd)
        os.close(qa_fd)
