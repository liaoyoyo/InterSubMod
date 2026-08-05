#!/usr/bin/env python3
"""Read-only runtime evidence for external review of the frozen release sources."""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
from typing import Any


sys.dont_write_bytecode = True


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
AUTHORITY_SOURCE = TOPIC_ROOT / "scripts" / "release_source_authority.py"
PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")

EXPECTED_SOURCES = {
    "P": (
        AUDIT_ROOT / "promote_tumor_ref_source_receipt_v2.py",
        157_669,
        "02fb9039b362fa261619b2236ddb764db278b23ae4467472fe2caa106770e06c",
    ),
    "V": (
        AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v2.py",
        120_163,
        "03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8",
    ),
    "R": (
        AUDIT_ROOT / "replay_m2v5_runner_only_gates_v1.py",
        118_987,
        "10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694",
    ),
    "C": (
        AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_v1.py",
        277_381,
        "f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd",
    ),
    "release_source_authority": (
        AUTHORITY_SOURCE,
        36_728,
        "36621dc735d17cf3322e34aeb5d72247869935e7b87306cdb600ff6d522c52af",
    ),
}


def artifact_record(path: Path, expected_mode: int | None = 0o444) -> dict[str, Any]:
    canonical = path.resolve(strict=True)
    data = canonical.read_bytes()
    opened = canonical.stat()
    if (
        canonical != path
        or not stat.S_ISREG(opened.st_mode)
        or (
            expected_mode is not None
            and stat.S_IMODE(opened.st_mode) != expected_mode
        )
        or opened.st_nlink != 1
    ):
        raise RuntimeError(f"Frozen artifact contract drift: {path}")
    return {
        "path": str(canonical),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(stat.S_IMODE(opened.st_mode)),
    }


def load_module(tag: str, path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(f"formal_probe_{tag}", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot load frozen module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def normalize_runtime_contracts(module: Any, tag: str) -> dict[str, dict[str, str]]:
    if tag == "P":
        return {
            role: {
                "path": str(module.REVIEWED_SOURCE_ROLES[role]),
                "sha256": contract[0],
                "mode": contract[1],
            }
            for role, contract in module.REVIEWED_SOURCE_RUNTIME_CONTRACTS.items()
        }
    normalized = {}
    for role, contract in module.REVIEWED_RUNTIME_SOURCE_CONTRACTS.items():
        if isinstance(contract, tuple):
            path, sha256, mode = contract
        else:
            path = contract["path"]
            sha256 = contract["sha256"]
            mode = contract["mode"]
        normalized[role] = {
            "path": str(path),
            "sha256": sha256,
            "mode": mode,
        }
    return normalized


def recursive_source_self_check(tag: str, path: Path, function_name: str) -> Any:
    code = """import json,sys
from pathlib import Path
p=Path(sys.argv[1]).resolve()
function_name=sys.argv[2]
data=p.read_bytes()
namespace={\"__name__\":\"formal_probe_shadow\",\"__file__\":str(p),\"__package__\":\"\"}
result={}
def trace(frame,event,arg):
    if (event==\"return\" and frame.f_code.co_name==\"<module>\"
            and frame.f_code.co_filename==str(p)
            and function_name in namespace and not result):
        observed=namespace[function_name](data)
        result[\"value\"]=\"PASS\" if observed is None else observed
    return trace
sys.settrace(trace)
exec(compile(data,str(p),\"exec\"),namespace)
sys.settrace(None)
print(json.dumps(result,sort_keys=True))
"""
    result = subprocess.run(
        [str(PYTHON), "-I", "-B", "-c", code, str(path), function_name],
        cwd=REPO_ROOT,
        env={"HOME": str(Path.home()), "PATH": "/usr/bin:/bin", "TZ": "UTC", "LC_ALL": "C"},
        capture_output=True,
        check=False,
        timeout=120,
    )
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(
            f"Recursive source self-check failed for {tag}: exit={result.returncode}"
        )
    payload = json.loads(result.stdout)
    if "value" not in payload:
        raise RuntimeError(f"Recursive source self-check produced no result: {tag}")
    return payload["value"]


def producer_wait_fault_probe(module: Any) -> dict[str, Any]:
    def exit_zero(_: dict[str, Any]) -> None:
        return None

    def exit_one(_: dict[str, Any]) -> None:
        os._exit(1)

    def terminate(_: dict[str, Any]) -> None:
        os.kill(os.getpid(), signal.SIGTERM)

    def kill(_: dict[str, Any]) -> None:
        os.kill(os.getpid(), signal.SIGKILL)

    expected = {
        "exit_zero": (True, 0, False, None, 0),
        "exit_one": (True, 1, False, None, 256),
        "sigterm": (False, None, True, signal.SIGTERM, signal.SIGTERM),
        "sigkill": (False, None, True, signal.SIGKILL, signal.SIGKILL),
    }
    observed = {}
    for label, action in (
        ("exit_zero", exit_zero),
        ("exit_one", exit_one),
        ("sigterm", terminate),
        ("sigkill", kill),
    ):
        session, child_wait = module.fork_no_exec_and_wait(
            "--prepare-authorization", action
        )
        module.validate_producer_session(
            session, expected_mode="--prepare-authorization"
        )
        module.validate_child_wait(
            child_wait, expected_pid=session["child_pid"], require_success=False
        )
        actual = (
            child_wait["wifexited"],
            child_wait["exit_status"],
            child_wait["wifsignaled"],
            child_wait["terminating_signal"],
            child_wait["raw_wait_status"],
        )
        if actual != expected[label]:
            raise RuntimeError(f"Producer child wait decoding drift: {label}")
        rejected = False
        try:
            module.validate_child_wait(
                child_wait,
                expected_pid=session["child_pid"],
                require_success=True,
            )
        except module.PromotionError:
            rejected = True
        if label == "exit_zero" and rejected:
            raise RuntimeError("Normal producer exit zero was rejected")
        if label != "exit_zero" and not rejected:
            raise RuntimeError(f"Abnormal producer child was accepted: {label}")
        observed[label] = {
            "raw_wait_status": child_wait["raw_wait_status"],
            "normal_success_accepted": label == "exit_zero",
            "abnormal_success_rejected": label != "exit_zero",
        }
    return observed


def main() -> None:
    records: dict[str, dict[str, Any]] = {}
    modules: dict[str, Any] = {}
    for tag, (path, expected_size, expected_sha256) in EXPECTED_SOURCES.items():
        record = artifact_record(path)
        if (
            record["size_bytes"] != expected_size
            or record["sha256"] != expected_sha256
        ):
            raise RuntimeError(f"Frozen source identity drift: {tag}")
        data = path.read_bytes()
        ast.parse(data, filename=str(path))
        compile(data, str(path), "exec")
        records[tag] = record
        if tag in {"P", "V", "R", "C"}:
            modules[tag] = load_module(tag, path)

    p_module = modules["P"]
    v_module = modules["V"]
    r_module = modules["R"]
    c_module = modules["C"]
    source_self_checks = {
        "P": recursive_source_self_check(
            "P", EXPECTED_SOURCES["P"][0],
            "verify_executing_source_matches_bound_bytes",
        ),
        "V": recursive_source_self_check(
            "V", EXPECTED_SOURCES["V"][0],
            "_require_executing_verifier_matches_bound_source",
        ),
        "R": recursive_source_self_check(
            "R", EXPECTED_SOURCES["R"][0],
            "require_executing_replayer_matches_bound_source",
        ),
        "C": recursive_source_self_check(
            "C", EXPECTED_SOURCES["C"][0], "verify_executing_source"
        ),
    }
    if source_self_checks["P"].get("pass") is not True:
        raise RuntimeError("Promotion recursive source self-check did not pass")
    if source_self_checks["C"].get("pass") is not True:
        raise RuntimeError("Continuation recursive source self-check did not pass")
    if source_self_checks["V"] != "PASS" or source_self_checks["R"] != "PASS":
        raise RuntimeError("Verifier/replayer recursive source self-check did not pass")

    matrix_values = {
        "authorization_checks": [
            module.AUTHORIZATION_CHECKS for module in (p_module, v_module, r_module, c_module)
        ],
        "completion_checks": [
            module.COMPLETION_CHECKS for module in (p_module, v_module, r_module, c_module)
        ],
        "prepare_exit_checks": [
            p_module.PREPARE_EXIT_CHECKS,
            v_module.PREPARE_EXIT_ATTESTATION_CHECKS,
            r_module.PREPARE_EXIT_ATTESTATION_CHECKS,
            c_module.PREPARE_EXIT_CHECKS,
        ],
        "promote_exit_checks": [
            p_module.PROMOTE_EXIT_CHECKS,
            v_module.PROMOTE_EXIT_ATTESTATION_CHECKS,
            r_module.PROMOTE_EXIT_ATTESTATION_CHECKS,
            c_module.PROMOTE_EXIT_CHECKS,
        ],
    }
    for name, values in matrix_values.items():
        if any(value != values[0] for value in values[1:]):
            raise RuntimeError(f"P/V/R/C exact check-set drift: {name}")

    runtime_contracts = {
        tag: normalize_runtime_contracts(modules[tag], tag)
        for tag in ("P", "V", "R", "C")
    }
    if (
        len(runtime_contracts["P"]) != 11
        or any(value != runtime_contracts["P"] for value in runtime_contracts.values())
    ):
        raise RuntimeError("P/V/R/C reviewed 4+11 runtime closure drift")

    prepare_schemas = [
        p_module.PREPARE_EXIT_ATTESTATION_SCHEMA,
        {
            "name": v_module.PREPARE_EXIT_ATTESTATION_SCHEMA,
            "version": v_module.EXIT_ATTESTATION_SCHEMA_VERSION,
        },
        {"name": r_module.PREPARE_EXIT_ATTESTATION_SCHEMA, "version": "1.0.0"},
        c_module.PREPARE_EXIT_ATTESTATION_SCHEMA,
    ]
    promote_schemas = [
        p_module.PROMOTE_EXIT_ATTESTATION_SCHEMA,
        {
            "name": v_module.PROMOTE_EXIT_ATTESTATION_SCHEMA,
            "version": v_module.EXIT_ATTESTATION_SCHEMA_VERSION,
        },
        {"name": r_module.PROMOTE_EXIT_ATTESTATION_SCHEMA, "version": "1.0.0"},
        c_module.PROMOTE_EXIT_ATTESTATION_SCHEMA,
    ]
    if any(value != prepare_schemas[0] for value in prepare_schemas[1:]):
        raise RuntimeError("Prepare exit-attestation schema drift")
    if any(value != promote_schemas[0] for value in promote_schemas[1:]):
        raise RuntimeError("Promote exit-attestation schema drift")

    p_reviews = [
        (
            role,
            slot["reviewer"],
            slot["reviewer_agent_id"],
            str(slot["path"]),
        )
        for role, slot in p_module.REVIEW_SLOTS.items()
    ]
    v_reviews = [
        (role, reviewer, agent_id, str(path))
        for role, reviewer, agent_id, path in v_module.REVIEW_SLOTS
    ]
    if p_reviews != v_reviews or len(p_reviews) != 3:
        raise RuntimeError("Promotion three-review slot contract drift")
    slots = {
        "P": tuple(map(str, p_module.DOWNSTREAM_OUTPUTS)),
        "V": tuple(map(str, v_module.DOWNSTREAM_OUTPUTS)),
        "R": tuple(map(str, r_module.DOWNSTREAM_OUTPUT_SLOTS)),
        "C": tuple(map(str, c_module.DOWNSTREAM_OUTPUT_SLOTS)),
    }
    if len(slots["C"]) != 17 or any(value != slots["C"] for value in slots.values()):
        raise RuntimeError("P/V/R/C downstream slot contract drift")
    for name in (
        "SUPERVISOR_CAPABILITY_PROTOCOL",
        "CONTINUATION_THREAT_MODEL",
        "CONTINUATION_POLICY",
    ):
        values = [getattr(modules[tag], name) for tag in ("P", "V", "R", "C")]
        if any(value != values[0] for value in values):
            raise RuntimeError(f"P/V/R/C contract drift: {name}")

    supervisor_commands = [
        [
            str(p_module.PYTHON),
            "-I",
            "-B",
            str(p_module.DOWNSTREAM_CONTINUATION),
            "--supervise-and-sign",
        ],
        v_module._downstream_continuation_command(),
        r_module.downstream_continuation_command(),
        c_module.canonical_command(),
    ]
    if any(command != supervisor_commands[0] for command in supervisor_commands):
        raise RuntimeError("P/V/R/C supervisor command drift")

    private_modes = {}
    private_headers = {}
    for label, path in (
        ("authorization", p_module.AUTHORIZATION_PRIVATE_KEY),
        ("completion", p_module.COMPLETION_PRIVATE_KEY),
        ("continuation", p_module.CONTINUATION_PRIVATE_KEY),
    ):
        opened = path.stat()
        mode = stat.S_IMODE(opened.st_mode)
        if mode != 0o400 or opened.st_nlink != 1:
            raise RuntimeError(f"One-time private key pre-signature drift: {label}")
        private_modes[label] = oct(mode)
        with path.open("rb") as stream:
            private_headers[label] = stream.readline().decode("ascii").strip()
    if private_headers != {
        "authorization": "-----BEGIN ENCRYPTED PRIVATE KEY-----",
        "completion": "-----BEGIN ENCRYPTED PRIVATE KEY-----",
        "continuation": "-----BEGIN PRIVATE KEY-----",
    }:
        raise RuntimeError("One-time private-key encryption/residency contract drift")
    continuation_public = subprocess.run(
        [
            "/usr/bin/openssl",
            "pkey",
            "-in",
            str(p_module.CONTINUATION_PRIVATE_KEY),
            "-pubout",
        ],
        cwd=REPO_ROOT,
        env={"PATH": "/usr/bin:/bin", "HOME": "/tmp", "LC_ALL": "C"},
        capture_output=True,
        check=False,
        timeout=30,
    )
    if (
        continuation_public.returncode != 0
        or continuation_public.stderr
        or continuation_public.stdout
        != p_module.CONTINUATION_PUBLIC_KEY.read_bytes()
    ):
        raise RuntimeError("Resident continuation private/public key pair drift")

    reviews = [Path(slot["path"]) for slot in p_module.REVIEW_SLOTS.values()]
    promotion = [
        p_module.AUTHORIZATION,
        p_module.PREPARE_EXIT_ATTESTATION,
        p_module.AUTHORIZATION_SIGNATURE,
        p_module.LEGACY_AUTHORIZATION_PAYLOAD_SIGNATURE,
        p_module.CANONICAL_RECEIPT,
        p_module.COMPLETION_RECEIPT,
        p_module.PROMOTE_EXIT_ATTESTATION,
        p_module.COMPLETION_SIGNATURE,
        p_module.LEGACY_COMPLETION_PAYLOAD_SIGNATURE,
        p_module.VERIFICATION_RECEIPT,
        p_module.RUNNER_GATE_RECEIPT,
        p_module.RUNNER_GATE_LOG,
        p_module.PROMOTION_INCIDENT,
        c_module.CONTINUATION_RECEIPT,
        c_module.CONTINUATION_EXIT_ATTESTATION,
        c_module.CONTINUATION_SIGNATURE,
        c_module.CONTINUATION_INCIDENT,
    ]
    occupied = [
        str(path)
        for path in [*reviews, *promotion, *p_module.DOWNSTREAM_OUTPUTS]
        if os.path.lexists(path)
    ]
    if occupied:
        raise RuntimeError(f"Formal output slots are not pristine: {occupied}")

    producer_faults = producer_wait_fault_probe(p_module)
    occupied_after_faults = [
        str(path)
        for path in [*reviews, *promotion, *p_module.DOWNSTREAM_OUTPUTS]
        if os.path.lexists(path)
    ]
    if occupied_after_faults:
        raise RuntimeError(
            f"Producer wait fault probe created formal output: {occupied_after_faults}"
        )

    reservation = c_module.reserve_supervisor_capability_descriptor()
    authority_module, bootstrap_record = c_module.bootstrap_release_authority()
    v7_validation = c_module.validate_v7(authority_module)
    leases = [os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC) for _ in range(100)]
    try:
        if c_module.SUPERVISOR_CAPABILITY_FD in leases:
            raise RuntimeError("FD 198 reservation was consumed by a validation lease")
        python_record = artifact_record(PYTHON.resolve(strict=True), None)
        capability_fd, capability_record = c_module.create_supervisor_capability(
            {
                "self_record": records["C"],
                "python_record": python_record,
            },
            reservation,
        )
        try:
            capability_stat = os.fstat(capability_fd)
            capability_seals = c_module.fcntl.fcntl(
                capability_fd, c_module.LINUX_F_GET_SEALS
            )
            if (
                capability_fd != c_module.SUPERVISOR_CAPABILITY_FD
                or stat.S_IMODE(capability_stat.st_mode) != 0o400
                or capability_stat.st_nlink != 0
                or capability_seals != c_module.supervisor_capability_seal_mask()
                or capability_record.get("pass") is not True
            ):
                raise RuntimeError("Sealed supervisor capability runtime drift")
        finally:
            os.close(capability_fd)
    finally:
        for descriptor in leases:
            os.close(descriptor)

    clean_environment = dict(c_module.EXPECTED_ENVIRONMENT)
    direct_child = subprocess.run(
        c_module.supervised_child_command(),
        cwd=REPO_ROOT,
        env=clean_environment,
        capture_output=True,
        check=False,
        timeout=120,
    )
    if (
        direct_child.returncode != 1
        or b"lacks the inherited supervisor capability descriptor" not in direct_child.stderr
    ):
        raise RuntimeError("Direct supervised-child invocation was not rejected")

    python_fd = os.open(PYTHON, os.O_RDONLY | os.O_CLOEXEC)
    try:
        fd_python_token = f"/proc/self/fd/{python_fd}"
        executable_probe = subprocess.run(
            [
                fd_python_token,
                "-I",
                "-B",
                "-c",
                "import sys; print(sys.executable)",
            ],
            executable=fd_python_token,
            cwd=REPO_ROOT,
            env=clean_environment,
            pass_fds=(python_fd,),
            capture_output=True,
            check=False,
            timeout=30,
        )
        if (
            executable_probe.returncode != 0
            or executable_probe.stderr
            or executable_probe.stdout.decode("ascii").strip() != fd_python_token
            or c_module.expected_finalizer_command(
                "create-dataset", python_launcher=fd_python_token
            )[0]
            != fd_python_token
            or c_module.expected_finalizer_command(
                "create-report", python_launcher=fd_python_token
            )[0]
            != fd_python_token
        ):
            raise RuntimeError("FD-bound finalizer command normalization drift")
    finally:
        os.close(python_fd)

    python_fd = os.open(PYTHON, os.O_RDONLY | os.O_CLOEXEC)
    source_fd = os.open(EXPECTED_SOURCES["C"][0], os.O_RDONLY | os.O_CLOEXEC)
    try:
        fd_signed_verifier = subprocess.run(
            [
                str(PYTHON),
                "-I",
                "-B",
                f"/proc/self/fd/{source_fd}",
                "--verify-signed-terminal",
            ],
            executable=f"/proc/self/fd/{python_fd}",
            cwd=REPO_ROOT,
            env=clean_environment,
            pass_fds=(python_fd, source_fd),
            capture_output=True,
            check=False,
            timeout=120,
        )
    finally:
        os.close(source_fd)
        os.close(python_fd)
    if (
        fd_signed_verifier.returncode != 1
        or b"canonical direct command" in fd_signed_verifier.stderr
        or b"Python launcher token drift" in fd_signed_verifier.stderr
        or b"Unexpected mode" in fd_signed_verifier.stderr
    ):
        raise RuntimeError("FD-bound signed-verifier source launch was not accepted")
    occupied_after = [
        str(path)
        for path in [*reviews, *promotion, *p_module.DOWNSTREAM_OUTPUTS]
        if os.path.lexists(path)
    ]
    if occupied_after:
        raise RuntimeError(f"Read-only probe created formal output: {occupied_after}")

    output = {
        "schema_name": "intersubmod.external_formal_readonly_probe",
        "schema_version": "1.0.0",
        "frozen_sources": records,
        "code_ast_compile_pass": True,
        "release_source_authority_runtime": {
            "source": bootstrap_record,
            "v7_authority_id": v7_validation["authority_id"],
            "v7_source_set_sha256": v7_validation["source_set_sha256"],
            "descriptor_lease": v7_validation["validation_binding_lease"],
            "pass": v7_validation["pass"],
        },
        "contract_matrix": {
            "downstream_slot_count": 17,
            "downstream_slots_equal": True,
            "capability_protocol_equal": True,
            "threat_model_equal": True,
            "policy_equal": True,
            "supervisor_commands_equal": True,
            "authorization_check_count": len(p_module.AUTHORIZATION_CHECKS),
            "completion_check_count": len(p_module.COMPLETION_CHECKS),
            "prepare_exit_check_count": len(p_module.PREPARE_EXIT_CHECKS),
            "promote_exit_check_count": len(p_module.PROMOTE_EXIT_CHECKS),
            "runtime_closure_role_count": len(runtime_contracts["P"]),
            "runtime_closure_equal": True,
            "exit_attestation_schemas_equal": True,
            "review_slot_count": len(p_reviews),
            "review_slots_equal": True,
        },
        "recursive_source_self_checks": source_self_checks,
        "producer_wait_faults": producer_faults,
        "fd_bound_finalizer_command": {
            "sys_executable_matches_inherited_fd_token": True,
            "dataset_command_exact": True,
            "report_command_exact": True,
        },
        "fd_bound_signed_verifier_source": {
            "exit_code_before_promotion": fd_signed_verifier.returncode,
            "runtime_invocation_contract_accepted": True,
            "formal_output_written": False,
        },
        "capability_runtime": {
            "reserved_before_authority_validation": True,
            "validation_lease_count_after_reservation": len(leases),
            "reserved_fd_consumed_by_lease": False,
            "descriptor": c_module.SUPERVISOR_CAPABILITY_FD,
            "mode": "0o400",
            "link_count": 0,
            "seal_mask": c_module.supervisor_capability_seal_mask(),
            "pass": True,
        },
        "direct_child_without_capability": {
            "exit_code": direct_child.returncode,
            "rejected_before_downstream": True,
        },
        "private_key_pre_signature_modes": private_modes,
        "key_material_contract": {
            "authorization_encrypted_live_signer_required": True,
            "completion_encrypted_live_signer_required": True,
            "continuation_unencrypted_resident_one_time_key": True,
            "continuation_private_public_pair_verified": True,
            "continuation_required_post_signature_mode": "0o0",
        },
        "review_slots_absent": True,
        "promotion_slots_absent": True,
        "downstream_slots_absent": True,
        "filesystem_artifacts_written": 0,
        "pass": True,
        "pass_semantics": "read_only_release_integrity_runtime_evidence_not_scientific_confirmation",
    }
    print(json.dumps(output, ensure_ascii=True, sort_keys=True))


if __name__ == "__main__":
    main()
