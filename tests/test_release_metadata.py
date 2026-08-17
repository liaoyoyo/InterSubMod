from __future__ import annotations

import json
import os
import re
from pathlib import Path

import yaml
from jsonschema.validators import validator_for


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(os.environ.get("INTERSUBMOD_TEST_SOURCE_ROOT", REPO_ROOT)).resolve()
MARKDOWN_METADATA = (
    REPO_ROOT / "CONTRIBUTING.md",
    REPO_ROOT / "SECURITY.md",
    REPO_ROOT / "CHANGELOG.md",
)
BOUNDARY_PATTERNS = (
    r"zero\s+confirmed\s+cellular\s+subclones",
    r"zero\s+confirmed\s+linear\s+ancestry",
    r"88\.2579",
    r"P3/P4/P5/P7/P8",
    r"non[- ]production",
    r"non[- ]clinical",
)
SCHEMA = REPO_ROOT / (
    "docs/handoff/20260813_完整研究資料與軟體交接_01/"
    "schemas/research-handoff-release-manifest.schema.json"
)
AGGREGATE_SCHEMA = REPO_ROOT / (
    "docs/handoff/20260813_完整研究資料與軟體交接_01/"
    "schemas/aggregate-release-acceptance.schema.json"
)
LONGLINEAGE_SCHEMA = REPO_ROOT / (
    "docs/handoff/20260813_完整研究資料與軟體交接_01/"
    "schemas/longlineage-preview-safety-receipt.schema.json"
)
COMMON_GATE_SCHEMA = REPO_ROOT / (
    "docs/handoff/20260813_完整研究資料與軟體交接_01/"
    "schemas/gate-acceptance-receipt.schema.json"
)
READER_SCHEMA = REPO_ROOT / (
    "docs/handoff/20260813_完整研究資料與軟體交接_01/"
    "schemas/reader-acceptance.schema.json"
)
LONG_LINEAGE_COMMIT = "b9aaa12a11fa00606bd174dabd0f172a5d112359"
AUTHORITY_SHA256 = "a88afc589206b7b3dab4d17c1d7b02a6cf20b125d847409649577c02abf0bfa0"
DENOMINATOR_SHA256 = "a41f726cbf66f22b7e95ddb44216a08bf480d12961758ebd5b1ab6e3e61db6a4"
REQUIRED_GATES = {
    "authority_19_of_19_replay",
    "registry_and_handoff_package",
    "repo_hygiene",
    "tracked_large_asset_policy",
    "clean_release_build",
    "ctest_complete",
    "python_tests_complete",
    "tiny_synthetic_e2e",
    "bip7_fresh_clone_acceptance",
    "bip8_fresh_clone_acceptance",
    "html_browser_qa_84",
    "docs_claim_link_validation",
    "public_claim_registry_158",
    "github_publication_commit_alignment",
    "github_branch_protection_required_ci",
    "full_history_secret_scan",
    "release_asset_sha256_roundtrip",
    "reader_ai_six_question_acceptance",
    "intersubmod_project_license_confirmation",
    "longlineage_preview_safety",
}


def normalized_prose(text: str) -> str:
    return re.sub(r"[*_`]+", "", text)


def test_citation_cff_is_bounded_and_omits_unverified_publication_metadata() -> None:
    citation = yaml.safe_load((REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    assert citation["cff-version"] == "1.2.0"
    assert citation["type"] == "software"
    assert citation["title"] == "InterSubMod"
    assert citation["authors"] == [{"name": "InterSubMod contributors"}]
    assert "license" not in citation
    assert not {"doi", "date-released", "version"}.intersection(citation)
    prose = " ".join(str(value) for value in citation.values())
    for pattern in BOUNDARY_PATTERNS:
        assert re.search(pattern, prose, flags=re.IGNORECASE), pattern
    assert "project-level SPDX choice still requires maintainer confirmation" in prose
    assert "not functions of the inter_sub_mod executable" in prose
    assert "stratification/annotation only" in prose


def test_release_markdown_boundaries_and_local_links() -> None:
    markdown_link = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")
    for source in MARKDOWN_METADATA:
        text = source.read_text(encoding="utf-8")
        assert text.startswith("<!--\n"), source
        prose = normalized_prose(text)
        for pattern in BOUNDARY_PATTERNS:
            assert re.search(pattern, prose, flags=re.IGNORECASE), (source, pattern)
        for raw_target in markdown_link.findall(text):
            target = raw_target.split("#", 1)[0]
            if not target or re.match(r"^[a-z][a-z0-9+.-]*://", target):
                continue
            candidate = (source.parent / target).resolve()
            if not candidate.exists() and source.parent == REPO_ROOT:
                candidate = (SOURCE_ROOT / target).resolve()
            assert candidate.exists(), (source, target)


def test_release_schema_can_never_promote_longlineage_scientific_phases() -> None:
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    gate_properties = schema["properties"]["longlineage"]["properties"]["gates"]["properties"]
    assert set(gate_properties) == {"P3", "P4", "P5", "P7", "P8"}
    assert {gate: contract["const"] for gate, contract in gate_properties.items()} == {
        "P3": "BLOCKED",
        "P4": "BLOCKED",
        "P5": "BLOCKED",
        "P7": "BLOCKED",
        "P8": "BLOCKED",
    }
    longlineage = schema["properties"]["longlineage"]["properties"]
    assert longlineage["candidate_commit"] == {"const": LONG_LINEAGE_COMMIT}
    assert longlineage["production_ready"] == {"const": False}
    assert "RESEARCH_PREVIEW_SAFETY_VERIFIED" in longlineage["preview_status"]["enum"]

    authority = schema["properties"]["authority"]["properties"]
    assert authority["manifest_sha256"] == {"const": AUTHORITY_SHA256}
    assert authority["denominator_registry_sha256"] == {"const": DENOMINATOR_SHA256}


def test_aggregate_acceptance_schema_has_exact_closed_gate_set() -> None:
    schema = json.loads(AGGREGATE_SCHEMA.read_text(encoding="utf-8"))
    assert schema["additionalProperties"] is False
    gates = schema["properties"]["gates"]
    assert gates["additionalProperties"] is False
    assert set(gates["required"]) == REQUIRED_GATES
    assert set(gates["properties"]) == REQUIRED_GATES
    gate_contract = schema["$defs"]["gate"]
    assert gate_contract["additionalProperties"] is False
    assert gate_contract["properties"]["status"] == {"const": "PASS"}
    assert gate_contract["properties"]["evidence"]["minItems"] == 1


def test_longlineage_safety_schema_is_closed_and_pins_candidate() -> None:
    schema = json.loads(LONGLINEAGE_SCHEMA.read_text(encoding="utf-8"))
    assert schema["additionalProperties"] is False
    assert schema["properties"]["longlineage_candidate_commit"] == {
        "const": LONG_LINEAGE_COMMIT
    }
    gates = schema["properties"]["preview_safety_gates"]
    assert gates["additionalProperties"] is False
    assert set(gates["required"]) == {"license", "source_origin", "dependencies", "public_safety"}
    assert set(schema["properties"]["phase_gates"]["properties"]) == {
        "P3",
        "P4",
        "P5",
        "P7",
        "P8",
    }
    assert "safety_stacked_commit" in schema["required"]


def test_common_gate_schema_requires_raw_log_binding_and_reader_agent() -> None:
    schema = json.loads(COMMON_GATE_SCHEMA.read_text(encoding="utf-8"))
    execution = schema["properties"]["execution"]
    assert {"log_uri", "log_sha256"}.issubset(execution["required"])
    assert execution["properties"]["log_uri"]["pattern"].startswith("^release-assets/")
    reader_condition = next(
        condition
        for condition in schema["allOf"]
        if condition["if"]["properties"]["gate_id"].get("const")
        == "reader_ai_six_question_acceptance"
    )
    assert reader_condition["then"]["properties"]["execution"]["properties"]["kind"] == {
        "const": "READER_AGENT"
    }


def test_all_release_metadata_schemas_are_valid_draft_2020_12() -> None:
    for path in (
        SCHEMA,
        AGGREGATE_SCHEMA,
        LONGLINEAGE_SCHEMA,
        COMMON_GATE_SCHEMA,
        READER_SCHEMA,
    ):
        schema = json.loads(path.read_text(encoding="utf-8"))
        validator_class = validator_for(schema)
        validator_class.check_schema(schema)
