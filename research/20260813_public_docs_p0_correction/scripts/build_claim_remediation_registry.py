#!/usr/bin/env python3
"""Build the complete 158-claim remediation registry from the frozen audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

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
    atomic_write_json,
    build_about_snapshot_manifest,
    build_counts,
    build_public_source_manifest,
    build_receipt_payload,
    canonical_json_bytes,
    canonical_object_sha256,
    sha256_bytes,
    sha256_path,
)
from github_about_evidence import EXPECTED_DESCRIPTION, validate_about_receipt
from validate_public_p0_claims import validate_registry as validate_p0_registry


ROOT = Path(__file__).resolve().parents[3]
PROJECT = ROOT / "research/20260813_public_docs_p0_correction"
INVENTORY = ROOT / "research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv"
SOURCE_SCOPE = ROOT / "research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv"
P0_REGISTRY = PROJECT / "scripts/p0_claim_registry.json"
OUTPUT = PROJECT / "claim_remediation_registry.json"
RECEIPT = PROJECT / "claim_remediation_build_receipt.json"
ABOUT_RECEIPT = PROJECT / "github_about_c108_receipt.json"

VERDICT_MAP = {
    "CONFIRMED": ("CONFIRMED", "FROZEN_EVIDENCE_CONFIRMED"),
    "CONFIRMED_WITH_LIMITS": ("CONFIRMED_WITH_LIMITS", "FROZEN_CONFIRMED_WITH_LIMITS"),
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

# Canonical bounded statements for claims whose frozen audit proposition or minimum
# rewrite is itself too broad.  current_verdict applies to these statements, never
# to the historical wording that triggered the audit.
CANONICAL_REMEDIATED_CLAIMS = {
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


def build() -> dict[str, object]:
    inventory = read_tsv(INVENTORY)
    inventory_sha256 = sha256_path(INVENTORY)
    inventory_ids = [row.get("claim_id") for row in inventory]
    if inventory_sha256 != EXPECTED_INVENTORY_SHA256:
        raise RuntimeError(
            f"frozen claim inventory hash drift: expected {EXPECTED_INVENTORY_SHA256}, got {inventory_sha256}"
        )
    if len(inventory_ids) != 158 or len(set(inventory_ids)) != 158:
        raise RuntimeError("frozen claim inventory must contain exactly 158 unique claim IDs")
    p0_bytes = P0_REGISTRY.read_bytes()
    p0 = json.loads(p0_bytes.decode("utf-8"))
    p0_sha256 = sha256_bytes(p0_bytes)
    p0_errors, p0_summary = validate_p0_registry(p0, ROOT)
    if p0_errors:
        formatted = "\n".join(f"- {item}" for item in p0_errors)
        raise RuntimeError(f"authoritative P0 guard registry failed validation:\n{formatted}")
    p0_map = {item["claim_id"]: item for item in p0["claims"]}
    about_receipt, about_errors = validate_about_receipt(ROOT, ABOUT_RECEIPT)
    if about_errors:
        formatted = "\n".join(f"- {item}" for item in about_errors)
        raise RuntimeError(f"GitHub About evidence failed semantic/snapshot validation:\n{formatted}")
    about_receipt_sha256 = sha256_path(ABOUT_RECEIPT)
    public_source_manifest = build_public_source_manifest(ROOT, p0)
    about_snapshot_manifest = build_about_snapshot_manifest(ROOT, about_receipt)
    source_map = build_source_map()
    claims: list[dict[str, object]] = []

    for row in inventory:
        claim_id = row["claim_id"]
        if claim_id in p0_map:
            p0_item = p0_map[claim_id]
            if p0_item["disposition"] == "external_action_required":
                current_verdict = "CONFIRMED_WITH_LIMITS"
                source_status = "EXTERNAL_LIVE_CONFIRMED_WITH_LIMITS"
            else:
                current_verdict = "CONFIRMED_WITH_LIMITS"
                source_status = "SOURCE_READY"
        elif claim_id in CURRENT_STATUS_OVERRIDES:
            current_verdict, source_status = CURRENT_STATUS_OVERRIDES[claim_id]
        elif row["verdict"] in VERDICT_MAP:
            current_verdict, source_status = VERDICT_MAP[row["verdict"]]
        else:
            current_verdict = "UNVERIFIED"
            source_status = "SOURCE_EDITED_REVIEW_REQUIRED"

        if claim_id in CANONICAL_REMEDIATED_CLAIMS:
            remediated_claim = CANONICAL_REMEDIATED_CLAIMS[claim_id]
        elif source_status == "SOURCE_READY":
            remediated_claim = f"Correction contract: {row['minimum_rewrite']}"
        elif source_status in {"FROZEN_EVIDENCE_CONFIRMED", "FROZEN_CONFIRMED_WITH_LIMITS"}:
            remediated_claim = row["proposition"]
        elif source_status == "EXTERNAL_LIVE_CONFIRMED_WITH_LIMITS":
            remediated_claim = EXPECTED_DESCRIPTION
        elif source_status == "SOURCE_EDITED_REVIEW_REQUIRED":
            remediated_claim = f"Proposed source correction under review: {row['minimum_rewrite']}"
        else:
            remediated_claim = None

        item: dict[str, object] = {
            "claim_id": claim_id,
            "category": row["category"],
            "inventory_proposition": row["proposition"],
            "remediated_claim": remediated_claim,
            "verdict_applies_to": "remediated_claim",
            "priority": row["priority"],
            "inventory_verdict": row["verdict"],
            "current_verdict": current_verdict,
            "source_status": source_status,
            "live_status": (
                about_receipt["live_status"]
                if claim_id == "C108"
                else "UNVERIFIED_AFTER_20260812_LOCKED_SNAPSHOT"
            ),
            "occurrences": row["occurrences"],
            "inventory_evidence": row["evidence"],
            "inventory_minimum_rewrite": row["minimum_rewrite"],
            "targets": occurrence_targets(row["occurrences"], source_map),
        }
        if claim_id in p0_map:
            item["p0_guard_registry"] = "scripts/p0_claim_registry.json"
            item["p0_local_check_count"] = len(p0_map[claim_id].get("checks", []))
            item["external_actions"] = p0_map[claim_id]["external_actions"]
        if claim_id == "C108":
            item["live_receipt"] = str(ABOUT_RECEIPT.relative_to(ROOT))
            item["live_receipt_sha256"] = about_receipt_sha256
        claims.append(item)

    return {
        "schema_name": REGISTRY_SCHEMA_NAME,
        "schema_version": REGISTRY_SCHEMA_VERSION,
        "registry_id": REGISTRY_ID,
        "generated_at": GENERATED_AT,
        "task_type": list(TASK_TYPE),
        "output_path": str(OUTPUT.relative_to(ROOT)),
        "build_receipt": str(RECEIPT.relative_to(ROOT)),
        "source_inventory": str(INVENTORY.relative_to(ROOT)),
        "source_inventory_sha256": inventory_sha256,
        "source_scope": str(SOURCE_SCOPE.relative_to(ROOT)),
        "source_scope_sha256": sha256_path(SOURCE_SCOPE),
        "p0_guard_registry": str(P0_REGISTRY.relative_to(ROOT)),
        "p0_guard_registry_sha256": p0_sha256,
        "p0_guard_validation_summary": p0_summary,
        "github_about_receipt": str(ABOUT_RECEIPT.relative_to(ROOT)),
        "github_about_receipt_sha256": about_receipt_sha256,
        "github_about_snapshot_manifest": about_snapshot_manifest,
        "github_about_snapshot_manifest_sha256": canonical_object_sha256(about_snapshot_manifest),
        "public_source_manifest": public_source_manifest,
        "public_source_manifest_sha256": canonical_object_sha256(public_source_manifest),
        "allowed_current_verdicts": list(ALLOWED_CURRENT_VERDICTS),
        "status_semantics": EXPECTED_STATUS_SEMANTICS,
        "gates": EXPECTED_GATES,
        "counts": build_counts(claims),
        "claims": claims,
    }


def main() -> int:
    registry = build()
    output_bytes = canonical_json_bytes(registry)
    receipt = build_receipt_payload(registry, sha256_bytes(output_bytes))
    atomic_write_json(OUTPUT, registry)
    atomic_write_json(RECEIPT, receipt)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
