#!/usr/bin/env python3
"""Fail-closed validation for the complete public-claim remediation registry."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any

from claim_registry_contract import (
    ALLOWED_CURRENT_VERDICTS,
    C066_CANONICAL_WORDING,
    EXPECTED_GATES,
    EXPECTED_INVENTORY_SHA256,
    EXPECTED_STATUS_SEMANTICS,
    GENERATED_AT,
    REGISTRY_ID,
    REGISTRY_SCHEMA_NAME,
    REGISTRY_SCHEMA_VERSION,
    TASK_TYPE,
    build_about_snapshot_manifest,
    build_counts,
    build_public_source_manifest,
    build_receipt_payload,
    canonical_json_bytes,
    canonical_object_sha256,
    sha256_bytes,
    sha256_path,
)
from github_about_evidence import validate_about_receipt
from validate_public_p0_claims import validate_registry as validate_p0_registry


ROOT = Path(__file__).resolve().parents[3]
PROJECT = ROOT / "research/20260813_public_docs_p0_correction"
DEFAULT_REGISTRY = PROJECT / "claim_remediation_registry.json"
DEFAULT_RECEIPT = PROJECT / "claim_remediation_build_receipt.json"
INVENTORY = ROOT / "research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv"
SOURCE_SCOPE = ROOT / "research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv"
DENOMINATORS = ROOT / "docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv"
P0_REGISTRY = PROJECT / "scripts/p0_claim_registry.json"
ABOUT_RECEIPT = PROJECT / "github_about_c108_receipt.json"
ALLOWED_VERDICTS = set(ALLOWED_CURRENT_VERDICTS)
ALLOWED_SOURCE = {
    "FROZEN_EVIDENCE_CONFIRMED",
    "FROZEN_CONFIRMED_WITH_LIMITS",
    "SOURCE_READY",
    "SOURCE_EDITED_REVIEW_REQUIRED",
    "EXTERNAL_ONLY_BLOCKED",
    "EXTERNAL_LIVE_CONFIRMED_WITH_LIMITS",
    "VALIDATED_DERIVED_WITH_LIMITS",
}
CURRENT_STATUS_OVERRIDES = {
    "C047": ("CONFIRMED_WITH_LIMITS", "VALIDATED_DERIVED_WITH_LIMITS"),
    "C048": ("CONFIRMED_WITH_LIMITS", "VALIDATED_DERIVED_WITH_LIMITS"),
    "C089": ("CONFIRMED_WITH_LIMITS", "VALIDATED_DERIVED_WITH_LIMITS"),
    "C114": ("CONFIRMED_WITH_LIMITS", "FROZEN_CONFIRMED_WITH_LIMITS"),
}
LOCAL_PATHS = {
    "A001": "EXTERNAL:GitHub About",
    "A002": "README.md",
    "A003": "README.md",
    "A004": "README.zh-TW.md",
    "A005": "QUICKSTART.md",
    "A006": "README_PROJECT_SUMMARY.md",
    "M001": "README.zh-TW.md",
}
EXPECTED_TOP_LEVEL_KEYS = {
    "schema_name",
    "schema_version",
    "registry_id",
    "generated_at",
    "task_type",
    "output_path",
    "build_receipt",
    "source_inventory",
    "source_inventory_sha256",
    "source_scope",
    "source_scope_sha256",
    "p0_guard_registry",
    "p0_guard_registry_sha256",
    "p0_guard_validation_summary",
    "github_about_receipt",
    "github_about_receipt_sha256",
    "github_about_snapshot_manifest",
    "github_about_snapshot_manifest_sha256",
    "public_source_manifest",
    "public_source_manifest_sha256",
    "allowed_current_verdicts",
    "status_semantics",
    "gates",
    "counts",
    "claims",
}
EXPECTED_DENOMINATORS = {
    "ssnv_dataset_records": (469849, 469849, "100.0000"),
    "strict_components": (255752, 255752, "100.0000"),
    "k1_linkage_abstain": (170131, 255752, "66.5219"),
    "ranked_complete": (71955, 75224, "95.6543"),
    "one_rooted_unlabeled_topology": (63506, 71955, "88.2579"),
    "methyl_evaluable_units": (811, 1045, "77.6077"),
    "methyl_robust_associations": (3, 811, "0.3699"),
}
EXPECTED_REMEDIATED_CLAIMS = {
    "C022": (
        "For 63,506 of 71,955 ranked-complete units, the frozen model's AF-best candidate set "
        "has one rooted-unlabeled topology signature within that unit. This is a "
        "model-conditional within-unit graph-shape summary, not evidence that units share one "
        "topology, nor a unique biological tree, accuracy, or prevalence."
    ),
    "C023": (
        "For 63,506 / 71,955 = 88.2579% of ranked-complete units, the frozen model's AF-best "
        "candidate set has one rooted-unlabeled topology signature within that unit. This is "
        "not evidence that units share one topology, nor biological correctness, accuracy, or "
        "prevalence."
    ),
    "C024": "Confirmed cellular subclone count = 0; no cellular subclone has been confirmed.",
    "C043": (
        "InterSubMod does not write BAM. LongLineage private baseline/main snapshot 5daf50f "
        "has no writer; private public-preview candidate b9aaa12 includes longlineage-tag-bam "
        "but remains NOT_READY/non-production with P3/P4/P5/P7/P8 blocked."
    ),
    "C047": (
        "The 2026-08-13 inventory identifies 14 historical PRE-FIX canonical tagged BAMs "
        "(7 paired_full + 7 paired_pileup) totaling 3,709,322,840,333 bytes. The seven "
        "paired_full objects total 1,840,983,466,353 bytes (1.6743647 TiB) and match a "
        "2026-07-11 first/middle/last-1MiB sampled-identity receipt. Full-file SHA-256 was not "
        "computed, so this remains a PARTIAL/NON_FINAL storage observation, not frozen authority."
    ),
    "C048": (
        "Exact stat sizes yield 294.2669× for seven historical paired_full BAMs divided by seven "
        "later frozen sidecars. Because producer runs/revisions and full content equivalence were "
        "not established, this remains PARTIAL/NON_FINAL and is only a cross-generation "
        "storage-footprint quotient, not causal "
        "compression, verified storage reduction, or byte-equivalent replacement; 287× is invalid."
    ),
    "C061": (
        "Historical revision 9098f11 used a 99-permutation StructureTest production "
        "override while label-first/config defaults were 999. Frozen release baseline ddd8909a "
        "uses 999 permutations with add-one p-floor 0.001. The 73afaeac-dirty snapshot is a "
        "historical byte-equivalent IN_PROGRESS/PARTIAL audit only and is excluded from the release source."
    ),
    "C062": (
        "A 99-permutation add-one grid has minimum attainable p approximately 0.01; this "
        "describes the historical 9098f11 StructureTest path, not the current 999-permutation "
        "frozen ddd8909a research contract."
    ),
    "C063": (
        "Historical StructureTest used a coarse three-bucket PERMDISP lookup. Current core "
        "uses an analytic F-distribution upper-tail p-value; neither version establishes "
        "cluster or cellular truth."
    ),
    "C066": C066_CANONICAL_WORDING,
    "C078": (
        "Local exact-PS×HP somatic molecular co-occurrence is non-circular only with "
        "respect to methylation-derived labels; it still depends on variant calling, "
        "basecalling, alignment and haplotag assumptions and does not establish cellular identity."
    ),
    "C087": (
        "All named checks in the frozen cohort receipt passed for the cited frozen artifacts; "
        "this is neither current-source certification nor a production/release gate PASS."
    ),
    "C089": (
        "Seven historical PRE-FIX paired_full BAMs are identified by exact paths/bytes/mtimes and "
        "a 7/7 sampled-identity replay; seven later frozen sidecars total 6,256,168,164 bytes. "
        "Their exact-size quotient is 294.2669×, but absent full-file hashes and generation-level "
        "producer/content equivalence it remains PARTIAL/NON_FINAL cross-generation footprint "
        "evidence, not a fully verified frozen storage snapshot or compression claim."
    ),
    "C093": (
        "Observed chr2:18M patterns support multiple local read-state candidates in one "
        "HCC1395 biological case; no minimum biological or complete candidate count is established."
    ),
    "C114": (
        "The frozen canonical/cohort_receipt.json records technical_all_pass=true, while "
        "summary/all7_summary.json records validation_evidence_eligible=false; these are not "
        "top-level authority_manifest.json fields. This names the exact frozen artifacts and "
        "does not imply current-source certification or production/release readiness."
    ),
    "C096": (
        "Within shared-coverage reads, pos5 ALT was observed only with pos3 ALT (HKU=4, "
        "DORADO=2), a provisional local nesting pattern that does not establish parent/child, "
        "order, ancestry, or lineage."
    ),
    "C126": (
        "The normal-ASM-screen-negative subset supports a bounded, mutation-conditioned "
        "methylation association; it does not establish acquired, clone-specific, or "
        "lineage-associated methylation."
    ),
    "C152": (
        "SEQC2 LOH annotation plus HP imbalance is compatible with LOH context; it may "
        "reduce but does not remove a simple two-parent explanation and does not confirm "
        "mechanism. CN, tag/alignment/sequence error and other alternatives remain."
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--simulate-drop-claim")
    parser.add_argument("--simulate-duplicate-claim")
    parser.add_argument("--simulate-illegal-verdict", action="store_true")
    parser.add_argument("--simulate-missing-evidence", action="store_true")
    parser.add_argument("--simulate-source-live-collapse", action="store_true")
    parser.add_argument("--simulate-denominator-drift", action="store_true")
    parser.add_argument("--simulate-capability-commit-drift", action="store_true")
    parser.add_argument("--simulate-unverified-upgrade", action="store_true")
    parser.add_argument("--simulate-p0-guard-hash-drift", action="store_true")
    parser.add_argument("--simulate-about-receipt-hash-drift", action="store_true")
    parser.add_argument("--simulate-schema-drift", action="store_true")
    parser.add_argument("--simulate-generated-at-drift", action="store_true")
    parser.add_argument("--simulate-task-type-drift", action="store_true")
    parser.add_argument("--simulate-source-path-drift", action="store_true")
    parser.add_argument("--simulate-allowed-verdicts-drift", action="store_true")
    parser.add_argument("--simulate-counts-drift", action="store_true")
    parser.add_argument("--simulate-gates-drift", action="store_true")
    parser.add_argument("--simulate-p0-replay-drift", action="store_true")
    parser.add_argument("--simulate-public-manifest-drift", action="store_true")
    parser.add_argument("--simulate-about-snapshot-drift", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_source_map() -> dict[str, str]:
    result = dict(LOCAL_PATHS)
    for row in read_tsv(SOURCE_SCOPE):
        artifact_id = row["artifact_id"]
        artifact = row["artifact"]
        if artifact_id.startswith("W"):
            result[artifact_id] = f"docs/wiki/{artifact}"
        elif artifact_id.startswith("P"):
            result[artifact_id] = f"docs/explain/{artifact}"
    return result


def occurrence_targets(value: str, source_map: dict[str, str]) -> list[str]:
    ids = [item.split(":", 1)[0] for item in value.split("|") if item]
    return sorted({source_map.get(item, f"UNMAPPED:{item}") for item in ids})


def _apply_simulations(registry: dict[str, Any], args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.simulate_p0_guard_hash_drift:
        registry["p0_guard_registry_sha256"] = "0" * 64
    if args.simulate_about_receipt_hash_drift:
        registry["github_about_receipt_sha256"] = "0" * 64
    if args.simulate_schema_drift:
        registry["schema_version"] = "999.0.0"
    if args.simulate_generated_at_drift:
        registry["generated_at"] = "2026-08-13T00:00:01+08:00"
    if args.simulate_task_type_drift:
        registry["task_type"] = ["A_EXPLORATORY_PILOT"]
    if args.simulate_source_path_drift:
        registry["source_inventory"] = "README.md"
    if args.simulate_allowed_verdicts_drift:
        registry["allowed_current_verdicts"] = [*ALLOWED_CURRENT_VERDICTS, "RESOLVED"]
    if args.simulate_counts_drift and isinstance(registry.get("counts"), dict):
        registry["counts"]["claims"] = 157
    if args.simulate_gates_drift and isinstance(registry.get("gates"), dict):
        registry["gates"]["PUBLICATION_READY"]["status"] = "PASS"
    if args.simulate_p0_replay_drift and isinstance(registry.get("p0_guard_validation_summary"), dict):
        registry["p0_guard_validation_summary"]["errors"] = 1
    if args.simulate_public_manifest_drift and isinstance(registry.get("public_source_manifest"), dict):
        registry["public_source_manifest"]["files"][0]["sha256"] = "0" * 64
    if args.simulate_about_snapshot_drift and isinstance(registry.get("github_about_snapshot_manifest"), dict):
        registry["github_about_snapshot_manifest"]["files"][0]["sha256"] = "0" * 64

    claims = registry.get("claims")
    if not isinstance(claims, list):
        return []
    if args.simulate_drop_claim:
        claims = [item for item in claims if item.get("claim_id") != args.simulate_drop_claim]
        registry["claims"] = claims
    if args.simulate_duplicate_claim and claims:
        match = next((item for item in claims if item.get("claim_id") == args.simulate_duplicate_claim), claims[0])
        claims.append(dict(match))
    if args.simulate_illegal_verdict and claims:
        claims[0]["current_verdict"] = "RESOLVED"
    if args.simulate_missing_evidence and claims:
        claims[0]["inventory_evidence"] = ""
    if args.simulate_source_live_collapse and claims:
        match = next((item for item in claims if item.get("source_status") == "SOURCE_READY"), None)
        if match:
            match["live_status"] = "CONFIRMED"
    if args.simulate_unverified_upgrade and claims:
        match = next((item for item in claims if item.get("source_status") == "SOURCE_EDITED_REVIEW_REQUIRED"), None)
        if match:
            match["current_verdict"] = "CONFIRMED_WITH_LIMITS"
    return claims


def validate(registry: dict[str, Any], args: argparse.Namespace) -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    claims = _apply_simulations(registry, args)

    actual_keys = set(registry)
    missing_keys = sorted(EXPECTED_TOP_LEVEL_KEYS - actual_keys)
    unknown_keys = sorted(actual_keys - EXPECTED_TOP_LEVEL_KEYS)
    if missing_keys:
        errors.append("missing top-level keys: " + ", ".join(missing_keys))
    if unknown_keys:
        errors.append("unknown top-level keys: " + ", ".join(unknown_keys))
    exact_top_level = {
        "schema_name": REGISTRY_SCHEMA_NAME,
        "schema_version": REGISTRY_SCHEMA_VERSION,
        "registry_id": REGISTRY_ID,
        "generated_at": GENERATED_AT,
        "task_type": TASK_TYPE,
        "output_path": str(DEFAULT_REGISTRY.relative_to(ROOT)),
        "build_receipt": str(DEFAULT_RECEIPT.relative_to(ROOT)),
        "source_inventory": str(INVENTORY.relative_to(ROOT)),
        "source_scope": str(SOURCE_SCOPE.relative_to(ROOT)),
        "p0_guard_registry": str(P0_REGISTRY.relative_to(ROOT)),
        "github_about_receipt": str(ABOUT_RECEIPT.relative_to(ROOT)),
        "allowed_current_verdicts": ALLOWED_CURRENT_VERDICTS,
        "status_semantics": EXPECTED_STATUS_SEMANTICS,
    }
    for key, expected in exact_top_level.items():
        if registry.get(key) != expected:
            errors.append(f"top-level {key} mismatch: expected {expected!r}, got {registry.get(key)!r}")

    inventory = read_tsv(INVENTORY)
    inventory_ids = [row["claim_id"] for row in inventory]
    inventory_sha = sha256_path(INVENTORY)
    if inventory_sha != EXPECTED_INVENTORY_SHA256:
        errors.append(f"frozen inventory hash drift: expected {EXPECTED_INVENTORY_SHA256}, got {inventory_sha}")
    if registry.get("source_inventory_sha256") != EXPECTED_INVENTORY_SHA256:
        errors.append("registry source_inventory_sha256 does not pin the frozen inventory")
    if registry.get("source_scope_sha256") != sha256_path(SOURCE_SCOPE):
        errors.append("registry source_scope_sha256 does not pin the source-to-target mapping")
    if len(inventory_ids) != 158 or len(set(inventory_ids)) != 158:
        errors.append("frozen inventory itself must contain exactly 158 unique claim IDs")

    p0_bytes = b""
    try:
        p0_bytes = P0_REGISTRY.read_bytes()
        p0 = json.loads(p0_bytes.decode("utf-8"))
        p0_errors, p0_summary = validate_p0_registry(p0, ROOT)
    except (OSError, UnicodeError, json.JSONDecodeError, TypeError, ValueError) as error:
        p0, p0_errors, p0_summary = {}, [f"P0 registry cannot be replayed: {error}"], {}
    errors.extend(f"P0 guard: {item}" for item in p0_errors)
    if not p0_bytes or registry.get("p0_guard_registry_sha256") != sha256_bytes(p0_bytes):
        errors.append("generated registry does not pin the current authoritative P0 guard registry hash")
    if registry.get("p0_guard_validation_summary") != p0_summary:
        errors.append("generated P0 validation summary does not exactly match a fresh replay")
    p0_claims = p0.get("claims", []) if isinstance(p0, dict) else []
    p0_map = {
        item.get("claim_id"): item
        for item in p0_claims
        if isinstance(item, dict) and isinstance(item.get("claim_id"), str)
    }

    about_receipt, about_errors = validate_about_receipt(ROOT, ABOUT_RECEIPT)
    errors.extend(f"GitHub About: {item}" for item in about_errors)
    about_sha = sha256_path(ABOUT_RECEIPT)
    if registry.get("github_about_receipt_sha256") != about_sha:
        errors.append("generated registry does not pin the C108 GitHub About receipt hash")
    try:
        about_manifest = build_about_snapshot_manifest(ROOT, about_receipt)
    except (OSError, ValueError, TypeError) as error:
        about_manifest = {}
        errors.append(f"GitHub About snapshot manifest cannot be rebuilt: {error}")
    if registry.get("github_about_snapshot_manifest") != about_manifest:
        errors.append("GitHub About snapshot manifest differs from re-read snapshot bytes")
    if registry.get("github_about_snapshot_manifest_sha256") != canonical_object_sha256(about_manifest):
        errors.append("GitHub About snapshot manifest hash mismatch")

    try:
        public_manifest = build_public_source_manifest(ROOT, p0)
    except (OSError, ValueError, TypeError, RuntimeError) as error:
        public_manifest = {}
        errors.append(f"public source byte/git manifest cannot be rebuilt: {error}")
    if registry.get("public_source_manifest") != public_manifest:
        errors.append("public source manifest differs from exact working/index/HEAD bytes")
    if registry.get("public_source_manifest_sha256") != canonical_object_sha256(public_manifest):
        errors.append("public source manifest hash mismatch")

    if not isinstance(registry.get("claims"), list):
        errors.append("claims must be a list")
        claims = []
    ids = [
        item.get("claim_id")
        for item in claims
        if isinstance(item, dict) and isinstance(item.get("claim_id"), str)
    ]
    if len(ids) != len(claims):
        errors.append("every claim must be an object with a string claim_id")
    duplicates = sorted(str(key) for key, value in Counter(ids).items() if value > 1)
    if duplicates:
        errors.append("duplicate claim IDs: " + ", ".join(duplicates))
    missing = sorted(set(inventory_ids) - set(ids))
    extra = sorted(str(item) for item in set(ids) - set(inventory_ids))
    if missing:
        errors.append("missing frozen claim IDs: " + ", ".join(missing))
    if extra:
        errors.append("extra claim IDs: " + ", ".join(extra))
    if len(claims) != 158:
        errors.append(f"claim count must be 158, got {len(claims)}")

    inventory_map = {row["claim_id"]: row for row in inventory}
    source_map = build_source_map()
    for item in claims:
        if not isinstance(item, dict):
            errors.append("every claim must be an object")
            continue
        claim_id = item.get("claim_id")
        row = inventory_map.get(claim_id)
        if item.get("current_verdict") not in ALLOWED_VERDICTS:
            errors.append(f"{claim_id}: illegal current_verdict {item.get('current_verdict')!r}")
        if item.get("source_status") not in ALLOWED_SOURCE:
            errors.append(f"{claim_id}: illegal source_status {item.get('source_status')!r}")
        for key in ("inventory_evidence", "inventory_proposition", "remediated_claim"):
            if not isinstance(item.get(key), str) or not item[key].strip():
                errors.append(f"{claim_id}: {key} is required")
        if not isinstance(item.get("targets"), list) or not item["targets"]:
            errors.append(f"{claim_id}: at least one exact target is required")
        if item.get("verdict_applies_to") != "remediated_claim":
            errors.append(f"{claim_id}: current_verdict must explicitly apply to remediated_claim")
        if (
            claim_id in EXPECTED_REMEDIATED_CLAIMS
            and item.get("remediated_claim") != EXPECTED_REMEDIATED_CLAIMS[claim_id]
        ):
            errors.append(f"{claim_id}: bounded canonical remediated claim drifted")
        if row:
            frozen_fields = {
                "category": row["category"],
                "inventory_proposition": row["proposition"],
                "priority": row["priority"],
                "inventory_verdict": row["verdict"],
                "occurrences": row["occurrences"],
                "inventory_evidence": row["evidence"],
                "inventory_minimum_rewrite": row["minimum_rewrite"],
            }
            for key, expected in frozen_fields.items():
                if item.get(key) != expected:
                    errors.append(f"{claim_id}: {key} drifted from frozen inventory")
            expected_targets = occurrence_targets(row["occurrences"], source_map)
            if item.get("targets") != expected_targets or any(
                str(target).startswith("UNMAPPED:") for target in expected_targets
            ):
                errors.append(f"{claim_id}: targets drifted or contain an unmapped source artifact")
        if (
            item.get("source_status") == "SOURCE_EDITED_REVIEW_REQUIRED"
            and item.get("current_verdict") != "UNVERIFIED"
        ):
            errors.append(f"{claim_id}: unguarded P1/P2 source edit must remain UNVERIFIED")
        expected_live = (
            "CONFIRMED_WITH_LIMITS_AFTER_REFETCH"
            if claim_id == "C108"
            else "UNVERIFIED_AFTER_20260812_LOCKED_SNAPSHOT"
        )
        if item.get("live_status") != expected_live:
            errors.append(f"{claim_id}: expected live_status {expected_live!r}, got {item.get('live_status')!r}")
        p0_item = p0_map.get(claim_id)
        if row and not p0_item:
            frozen_verdict = row["verdict"]
            expected = CURRENT_STATUS_OVERRIDES.get(
                str(claim_id),
                {
                    "CONFIRMED": ("CONFIRMED", "FROZEN_EVIDENCE_CONFIRMED"),
                    "CONFIRMED_WITH_LIMITS": ("CONFIRMED_WITH_LIMITS", "FROZEN_CONFIRMED_WITH_LIMITS"),
                }.get(frozen_verdict, ("UNVERIFIED", "SOURCE_EDITED_REVIEW_REQUIRED")),
            )
            actual = (item.get("current_verdict"), item.get("source_status"))
            if actual != expected:
                errors.append(f"{claim_id}: fail-closed verdict/status mapping expected {expected}, got {actual}")
        if p0_item and p0_item.get("disposition") != "external_action_required":
            if item.get("source_status") != "SOURCE_READY":
                errors.append(f"{claim_id}: guarded local P0 must be SOURCE_READY")
            if item.get("p0_guard_registry") != "scripts/p0_claim_registry.json":
                errors.append(f"{claim_id}: P0 claim must point to the authoritative guard registry")
            if item.get("p0_local_check_count") != len(p0_item.get("checks", [])):
                errors.append(f"{claim_id}: generated P0 local check count drifted from guard registry")
        if claim_id == "C108":
            if (item.get("current_verdict"), item.get("source_status")) != (
                "CONFIRMED_WITH_LIMITS",
                "EXTERNAL_LIVE_CONFIRMED_WITH_LIMITS",
            ):
                errors.append("C108 must record the bounded external live correction")
            if item.get("live_receipt") != str(ABOUT_RECEIPT.relative_to(ROOT)):
                errors.append("C108 must point to the checked-in live receipt")
            if item.get("live_receipt_sha256") != about_sha:
                errors.append("C108 live receipt hash mismatch")

    valid_claim_objects = [item for item in claims if isinstance(item, dict)]
    expected_counts = build_counts(valid_claim_objects)
    if registry.get("counts") != expected_counts:
        errors.append("top-level counts do not exactly match claim rows")
    if registry.get("gates") != EXPECTED_GATES:
        errors.append("top-level gates differ from the fixed fail-closed gate contract")

    expected_denominators = dict(EXPECTED_DENOMINATORS)
    if args.simulate_denominator_drift:
        expected_denominators["one_rooted_unlabeled_topology"] = (63507, 71955, "88.2579")
    denom = {row["metric"]: row for row in read_tsv(DENOMINATORS)}
    for metric, expected in expected_denominators.items():
        row = denom.get(metric)
        actual = None if row is None else (int(row["numerator"]), int(row["denominator"]), row["percent"])
        if actual != expected:
            errors.append(f"denominator drift for {metric}: expected {expected}, got {actual}")

    capability_text = "\n".join(
        (ROOT / path).read_text(encoding="utf-8")
        for path in ("README.md", "README.zh-TW.md", "docs/wiki/System-Overview.md")
    )
    if args.simulate_capability_commit_drift:
        capability_text = capability_text.replace("b9aaa12", "BADCOMMIT")
    for anchor in ("b9aaa12", "5daf50f", "longlineage-tag-bam", "199 columns", "VerificationSchemaVersion=2"):
        if anchor not in capability_text:
            errors.append(f"capability/schema anchor missing: {anchor}")
    if not re.search(r"Confirmed cellular subclones\*\* \| \*\*0", capability_text):
        errors.append("README must retain confirmed cellular subclones = 0")

    cross_guards = p0.get("cross_document_guards") if isinstance(p0, dict) else None
    required_paths = cross_guards.get("required_tracked_paths", []) if isinstance(cross_guards, dict) else []
    executable_paths = (
        set(cross_guards.get("required_executable_paths", []))
        if isinstance(cross_guards, dict)
        else set()
    )
    if not required_paths:
        errors.append("authoritative P0 registry must declare tracked tiny-fixture shipment paths")
    for relative in required_paths:
        if not (ROOT / relative).is_file():
            errors.append(f"documented tiny fixture path missing: {relative}")
    for relative in executable_paths:
        if not os.access(ROOT / relative, os.X_OK):
            errors.append(f"documented tiny fixture shell entrypoint is not executable: {relative}")

    summary = {
        "claims": len(claims),
        "unique_claim_ids": len(set(ids)),
        "current_verdicts": dict(
            sorted(Counter(str(item.get("current_verdict")) for item in valid_claim_objects).items())
        ),
        "source_statuses": dict(
            sorted(Counter(str(item.get("source_status")) for item in valid_claim_objects).items())
        ),
        "p0_guard": p0_summary,
        "github_about_receipt_sha256": about_sha,
        "public_source_files": public_manifest.get("counts", {}).get("files"),
        "github_about_snapshot_files": about_manifest.get("counts", {}).get("files"),
        "denominator_checks": len(expected_denominators),
        "tiny_fixture_paths": len(required_paths),
        "errors": len(errors),
        "verdict": "PASS" if not errors else "FAIL",
    }
    return errors, summary


def validate_build_receipt(
    registry: dict[str, Any], registry_bytes: bytes, receipt: dict[str, Any]
) -> list[str]:
    try:
        expected = build_receipt_payload(registry, sha256_bytes(registry_bytes))
    except (KeyError, TypeError, ValueError) as error:
        return [f"build receipt cannot be derived from registry: {error}"]
    if receipt != expected:
        return ["build receipt does not exactly bind output/source/P0/About/manifests/counts/gates"]
    return []


def main() -> int:
    args = parse_args()
    errors: list[str] = []
    try:
        registry_bytes = args.registry.read_bytes()
        registry = json.loads(registry_bytes.decode("utf-8"))
        if not isinstance(registry, dict):
            raise TypeError("registry root must be an object")
    except (OSError, UnicodeError, json.JSONDecodeError, TypeError) as error:
        print(json.dumps({"errors": 1, "verdict": "FAIL"}, indent=2))
        print(f"ERROR: registry cannot be parsed: {error}")
        return 1
    if registry_bytes != canonical_json_bytes(registry):
        errors.append("registry bytes are not canonical UTF-8 JSON")

    validation_errors, summary = validate(registry, args)
    errors.extend(validation_errors)
    try:
        receipt_bytes = args.receipt.read_bytes()
        receipt = json.loads(receipt_bytes.decode("utf-8"))
        if not isinstance(receipt, dict):
            raise TypeError("receipt root must be an object")
        if receipt_bytes != canonical_json_bytes(receipt):
            errors.append("build receipt bytes are not canonical UTF-8 JSON")
        errors.extend(validate_build_receipt(registry, registry_bytes, receipt))
    except (OSError, UnicodeError, json.JSONDecodeError, TypeError) as error:
        errors.append(f"build receipt cannot be parsed: {error}")

    summary["errors"] = len(errors)
    summary["verdict"] = "PASS" if not errors else "FAIL"
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    for error in errors:
        print(f"ERROR: {error}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
