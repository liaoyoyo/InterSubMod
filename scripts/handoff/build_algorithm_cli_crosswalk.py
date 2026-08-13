#!/usr/bin/env python3
"""Build the commit/hash-bound 35-row algorithm and CLI claim crosswalk."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "docs/handoff/20260813_完整研究資料與軟體交接_01"
MATRIX_RELATIVE = "research/20260812_intersubmod_github_public_docs_full_validation/algorithm_cli_matrix.tsv"
CLAIM_INVENTORY_RELATIVE = "research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv"
GUARD_RELATIVE = "research/20260813_public_docs_p0_correction/scripts/p0_claim_registry.json"
OUTPUT = PACKAGE / "evidence/algorithm_cli_claim_crosswalk.json"
VALIDATOR = "scripts/handoff/validate_algorithm_cli_crosswalk.py"
ASSESSED_GIT_COMMIT = "e83437ab9775d2cd3e27eb4e57f89f7d38126023"
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
EXPECTED_COLUMNS = [
    "claim_id",
    "public_claim",
    "occurrence",
    "source_symbol",
    "runtime_or_test",
    "verdict",
    "minimum_rewrite",
    "caveat",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_row_hash(row: dict[str, str]) -> str:
    payload = ("\t".join(row[column] for column in EXPECTED_COLUMNS) + "\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def snapshot_tree_hash(bindings: list[dict[str, str]]) -> str:
    payload = json.dumps(bindings, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def commit_path_sha256(commit: str, relative: str) -> str | None:
    probe = subprocess.run(
        ["git", "-C", str(ROOT), "show", f"{commit}:{relative}"],
        capture_output=True,
    )
    if probe.returncode != 0:
        return None
    return hashlib.sha256(probe.stdout).hexdigest()


def row(
    bounded: str,
    disposition: str,
    related: list[str],
    correction: str,
    targets: list[str],
    correction_note: str,
    guard_status: str,
    guard_ids: list[str],
    evidence_status: str,
    evidence_extra: list[str],
    ceiling: str,
    limits: list[str],
    release_gate: str,
    next_action: str,
) -> dict[str, object]:
    return {
        "bounded_current_claim": bounded,
        "current_disposition": disposition,
        "related_claim_ids": related,
        "source_correction": {"status": correction, "target_paths": targets, "note": correction_note},
        "guard": {"status": guard_status, "guard_ids": guard_ids, "validator": VALIDATOR},
        "evidence": {"status": evidence_status, "references": evidence_extra, "claim_ceiling": ceiling},
        "known_limits": limits,
        "release_gate": release_gate,
        "next_action": next_action,
    }


POLICY = {
    "ALG-001": row(
        "The normal Config/ArgParser path uses 16 threads when --threads is omitted; the current help string still says Default: 1.",
        "CONFIRMED_WITH_LIMITS", ["C052"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md"],
        "Public docs disclose the help/runtime mismatch instead of treating help text as authority.",
        "MISSING", [], "VERIFIED_WITH_LIMITS", [],
        "Applies to the normal Config plus ArgParser entry path at the assessed source snapshot.",
        ["The CLI help text itself remains stale and no dedicated cross-document guard covers this mismatch."],
        "SOURCE_READY_WITH_LIMITS", "Correct the help string or add a generated CLI-default conformance test before release.",
    ),
    "ALG-002": row(
        "OpenMP parallelism is primarily across regions; per-region distance work is configured to one worker in the inspected path.",
        "CONFIRMED_WITH_LIMITS", ["C051"], "LOCAL_SOURCE_CORRECTED",
        ["README_PROJECT_SUMMARY.md", "docs/wiki/System-Overview.md"],
        "Current summary scopes parallelism to regions.", "MISSING", [], "VERIFIED_WITH_LIMITS", [],
        "Source and selected tests support wiring, not a throughput or scaling guarantee.",
        ["No controlled cross-core scaling receipt is included."], "SOURCE_READY_WITH_LIMITS",
        "Add a commit- and hardware-pinned 1/4/8/16/32-thread benchmark if performance is cited.",
    ),
    "ALG-003": row(
        "No reproducible controlled 32-core scaling benchmark is available, so approximately linear scaling is not a current claim.",
        "CONFIRMED_WITH_LIMITS", ["C157"], "LOCAL_SOURCE_CORRECTED", ["README_PROJECT_SUMMARY.md"],
        "The unqualified scaling claim was replaced by an explicit evidence gap.", "MISSING", [], "HISTORICAL_ONLY", [],
        "This confirms the documentation/evidence boundary, not a measured scaling curve.",
        ["The historical estimate used four threads and 32 SNVs rather than a controlled core-count series."],
        "SOURCE_READY_WITH_LIMITS", "Keep the performance claim absent until a controlled scaling receipt is published.",
    ),
    "ALG-004": row(
        "The historical 222 ms mean is configuration-specific and does not establish an unqualified less-than-300-ms per-region runtime.",
        "CONFIRMED_WITH_LIMITS", ["C157"], "LOCAL_SOURCE_CORRECTED", ["README_PROJECT_SUMMARY.md"],
        "The universal runtime claim was removed.", "MISSING", [], "HISTORICAL_ONLY", [],
        "Historical timing and one slower counterexample bound the claim; neither is a current performance benchmark.",
        ["Hardware, I/O cache, window, threads and region distribution are not controlled across observations."],
        "SOURCE_READY_WITH_LIMITS", "Publish a versioned benchmark distribution before stating a runtime target.",
    ),
    "ALG-005": row(
        "The effective no-argument default is NHD; explicit NHD/BERNOULLI selections are honored, and repeated selections are retained.",
        "CONFIRMED_WITH_LIMITS", ["C053", "C143", "C144"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md", "docs/wiki/How-to-Run.md", "docs/explain/12_intersubmod-io.standalone.html", "docs/explain/16_how-to-run.standalone.html"],
        "Current public-source wording distinguishes omission from explicit metric selection.",
        "PASS", ["G_NHD_CLI_DEFAULT_SEMANTICS", "P0:C143", "P0:C144"], "VERIFIED_WITH_LIMITS", [],
        "This is parser behavior at the assessed commit, not a promise for other entry points or revisions.",
        ["Config.hpp still initializes BERNOULLI before the normal parser replaces the list."],
        "SOURCE_READY_WITH_LIMITS", "Keep the fail-closed semantic guard in the release validation suite.",
    ),
    "ALG-006": row(
        "Repeatable distance metrics are accepted; the first drives clustering and later metrics add output distance matrices only.",
        "CONFIRMED_WITH_LIMITS", ["C055"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md"],
        "Current README and Wiki sources state the order-sensitive first-metric behavior and repeatable outputs.",
        "PASS", ["G_NHD_CLI_DEFAULT_SEMANTICS"], "SOURCE_INSPECTION_ONLY", [],
        "The source loop establishes behavior, but no dedicated multi-metric runtime receipt is embedded.",
        ["Changing metric order can change clustering."], "SOURCE_READY_WITH_LIMITS",
        "Add a two-metric E2E assertion covering matrix outputs and first-metric clustering.",
    ),
    "ALG-007": row(
        "At the assessed source, inter_sub_mod --help returns exit code 1, which wrappers must not treat as a normal success code.",
        "CONFIRMED_WITH_LIMITS", [], "NO_CHANGE_REQUIRED", ["docs/wiki/InterSubMod-Engine.md"],
        "The nonstandard help exit is already disclosed.", "MISSING", [], "VERIFIED_WITH_LIMITS", [],
        "The observation is revision-specific and is not a POSIX compatibility guarantee.",
        ["No dedicated CI guard pins the help exit code."], "SOURCE_READY_WITH_LIMITS",
        "Either normalize help to exit 0 or add a revision-specific CLI test.",
    ),
    "ALG-008": row(
        "CLI option counts are revision-scoped and must be generated from the exact binary/help output rather than treated as a stable API count.",
        "CONFIRMED_WITH_LIMITS", [], "LOCAL_SOURCE_CORRECTED",
        ["docs/wiki/System-Overview.md", "docs/wiki/InterSubMod-Engine.md", "docs/explain/11_system-map-overview.standalone.html", "docs/explain/12_intersubmod-io.standalone.html"],
        "Current source pins the 29-option census to ddd8909a and labels the count branch-scoped.",
        "PASS", ["G_CLI_OPTION_COUNT_BRANCH_SCOPED"], "HISTORICAL_ONLY", [],
        "This bounds the count claim; it does not establish the option count of another branch.",
        ["Feature commit 73afaea and public baseline snapshots expose different interfaces; the underlying generated option census is not embedded here."],
        "SOURCE_READY_WITH_LIMITS", "Regenerate the option table from the exact release binary and retain the branch-scope guard.",
    ),
    "ALG-009": row(
        "At the assessed source, pmd_bed_path, pmd_gating and min_site_coverage have declarations but no repository-wide runtime use.",
        "CONFIRMED_WITH_LIMITS", ["C056"], "NO_CHANGE_REQUIRED", ["docs/wiki/InterSubMod-Engine.md"],
        "The Wiki labels these fields dormant rather than implemented method features.",
        "MISSING", [], "SOURCE_INSPECTION_ONLY", [],
        "The absence statement is pinned to the assessed source snapshot.",
        ["A later revision can wire these fields without changing their declarations."], "SOURCE_READY_WITH_LIMITS",
        "Add symbol-use and CLI coverage tests if the fields become active.",
    ),
    "ALG-010": row(
        "In the canonical RegionProcessor path, reads lacking required MM/ML modified-base tags are excluded from methylation analysis.",
        "CONFIRMED", ["C050"], "NO_CHANGE_REQUIRED",
        ["pipeline/README.md", "docs/wiki/InterSubMod-Engine.md", "docs/explain/03_methylation-read-filter.standalone.html"],
        "Current public sources state the MM/ML requirement.", "NOT_REQUIRED", [], "VERIFIED", [],
        "Confirmed for the canonical RegionProcessor path and tested missing-MM behavior.", [], "SOURCE_READY",
        "Preserve missing-tag regression coverage.",
    ),
    "ALG-011": row(
        "The parser consumes ML:B:C and extracts C+m blocks, but exhaustive support for every SAM MM grammar variant is not established.",
        "CONFIRMED_WITH_LIMITS", ["C050"], "LOCAL_SOURCE_CORRECTED",
        ["README_PROJECT_SUMMARY.md", "docs/explain/03_methylation-read-filter.standalone.html"],
        "The completeness claim is scoped to covered fixtures and test cases.", "MISSING", [], "VERIFIED_WITH_LIMITS", [],
        "Verified fixtures do not cover all MM grammar and modification combinations.",
        ["Other modification blocks are advanced over rather than modeled."], "SOURCE_READY_WITH_LIMITS",
        "Add named fixtures for every grammar variant before using complete-support language.",
    ),
    "ALG-012": row(
        "HP parsing has direct tests; PS and lineage auxiliary extraction are visible in source but lack equivalent direct extraction tests.",
        "CONFIRMED_WITH_LIMITS", ["C069"], "NO_CHANGE_REQUIRED",
        ["pipeline/README.md", "docs/handoff/20260813_完整研究資料與軟體交接_01/20260813_軟體輸入輸出與研究流程_01.md"],
        "The handoff distinguishes supported HP/ALT use from private lineage-tag preview fields.",
        "MISSING", [], "SOURCE_INSPECTION_ONLY", [],
        "Only HP integer/string forms have the direct tests cited by the audit.",
        ["An orphan lv without ls is discarded; germline-hp-only narrows accepted values."],
        "SOURCE_READY_WITH_LIMITS", "Add direct PS and lineage-aux extraction tests before claiming equivalent coverage.",
    ),
    "ALG-013": row(
        "Repeatable lineage-axis syntax is accepted, but the inspected evaluation loop processes only the first lineage axis.",
        "CONFIRMED_WITH_LIMITS", [], "LOCAL_SOURCE_CORRECTED", ["pipeline/README.md"],
        "The pipeline documentation now says only the first lineage axis is evaluated and denies independent multi-axis evaluation.",
        "PASS", ["G_LINEAGE_AXIS_FIRST_ONLY"], "SOURCE_INSPECTION_ONLY", [],
        "The conclusion applies to the inspected feature-source loop, which breaks after the first axis.",
        ["The lineage-axis interface is preview-only and excluded from the supported release baseline."],
        "SOURCE_READY_WITH_LIMITS", "Preserve the first-axis guard and add a direct runtime test before claiming broader support.",
    ),
    "ALG-014": row(
        "Per-read observations are sparse-like during collection, but the final methylation matrix is a dense rectangular read-by-CpG vector.",
        "CONFIRMED_WITH_LIMITS", [], "LOCAL_SOURCE_CORRECTED", ["README_PROJECT_SUMMARY.md"],
        "The summary now distinguishes temporary observations from the final representation.",
        "PASS", ["G_MATRIX_REPRESENTATION_DENSE_FINAL"], "SOURCE_INSPECTION_ONLY", [],
        "This is an implementation representation statement, not a measured memory benchmark.",
        ["Memory grows with the product of read and CpG counts."], "SOURCE_READY_WITH_LIMITS",
        "Add a representation regression or memory bound before publishing scalability estimates.",
    ),
    "ALG-015": row(
        "methylation.csv starts with a matrix row index; join it positionally to reads.tsv read_id, whose read_name is the BAM QNAME.",
        "CONFIRMED_WITH_LIMITS", ["C058"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md", "docs/explain/12_intersubmod-io.standalone.html"],
        "Current docs warn that the first column is not a read name.", "MISSING", [], "SOURCE_INSPECTION_ONLY", [],
        "The positional join holds for canonical producer order only.",
        ["Downstream sorting or editing can silently break positional alignment."], "SOURCE_READY_WITH_LIMITS",
        "Add a checked read key to methylation.csv or validate row alignment explicitly.",
    ),
    "ALG-016": row(
        "tree.nwk leaves are reads labeled by read_name and represent methylation-distance clustering, not clone phylogeny.",
        "CONFIRMED", ["C057"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md", "docs/handoff/20260813_完整研究資料與軟體交接_01/20260813_軟體輸入輸出與研究流程_01.md"],
        "All current entry paths state the read-dendrogram ceiling.", "NOT_REQUIRED", [], "VERIFIED", [],
        "Confirmed for the assessed implementation; invalid-distance reads can be removed before construction.", [],
        "SOURCE_READY", "Preserve leaf-label and no-cellular-phylogeny tests in E2E validation.",
    ),
    "ALG-017": row(
        "At frozen release baseline ddd8909a, significance_summary.csv has 199 columns with VerificationSchemaVersion=2 and RegionStratificationSchemaVersion=1.",
        "CONFIRMED_WITH_LIMITS", ["C059"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md", "docs/explain/12_intersubmod-io.standalone.html"],
        "Current source replaces the stale 193-column statement with the baseline-pinned 199-column contract.",
        "PASS", ["G_SCHEMA_CURRENT_199", "P0:C059"], "VERIFIED_WITH_LIMITS", [],
        "The count is commit-specific; there is still no single whole-file layout-version field.",
        ["Historical outputs have other column counts."], "SOURCE_READY_WITH_LIMITS",
        "Continue binding consumers by header name and producer commit.",
    ),
    "ALG-018": row(
        "Consumers must parse significance_summary.csv by column name and verify required component schema fields rather than fixed numeric positions.",
        "CONFIRMED_WITH_LIMITS", ["C059"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md", "docs/explain/12_intersubmod-io.standalone.html"],
        "Current docs require name-based parsing and commit/schema provenance.",
        "PASS", ["G_SCHEMA_CURRENT_199", "P0:C059"], "VERIFIED_WITH_LIMITS", [],
        "Name binding prevents positional drift but does not prove semantic compatibility.",
        ["No single whole-file semantic version exists."], "SOURCE_READY_WITH_LIMITS",
        "Keep required-column/schema validation in downstream consumers.",
    ),
    "ALG-019": row(
        "subclone_structure.txt is a legacy-named deprecation/region-stratification artifact and does not report inferred cellular subclones.",
        "CONFIRMED_WITH_LIMITS", ["C024"], "LOCAL_SOURCE_CORRECTED",
        ["docs/wiki/InterSubMod-Engine.md", "docs/explain/12_intersubmod-io.standalone.html", "docs/handoff/20260813_完整研究資料與軟體交接_01/20260813_研究結論時間與Finality_01.md"],
        "Current Wiki and standalone I/O tables identify the file as a legacy-named deprecation/region-stratification artifact and deny cellular-subclone inference.",
        "PASS", ["G_SUBCLONE_STRUCTURE_DEPRECATION", "G_CONFIRMED_CLAIM_CEILING_ZERO"], "VERIFIED_WITH_LIMITS", [],
        "Source shows the deprecation stub; the filename itself remains misleading.",
        ["The legacy filename remains for compatibility and must never be interpreted by name alone."],
        "SOURCE_READY_WITH_LIMITS", "Preserve the file-specific deprecation and zero-cellular-claim guards.",
    ),
    "ALG-020": row(
        "Allele-specific CN, purity/CCF and SAVANA-segment integration are absent; optional LOH BED annotation and coverage-proxy categories exist.",
        "CONFIRMED_WITH_LIMITS", ["C068"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/InterSubMod-Engine.md", "docs/handoff/20260813_完整研究資料與軟體交接_01/20260813_研究結論時間與Finality_01.md"],
        "Current docs distinguish optional LOH annotation/stratification from missing CN-adjusted biological inference.",
        "PASS", ["G_CN_LOH_INFERENCE_BOUNDARY"], "VERIFIED_WITH_LIMITS", [],
        "Coverage-proxy categories are not copy-number inference.",
        ["CN loss, cnLOH, multiplicity, purity and CCF remain unresolved."], "SOURCE_READY_WITH_LIMITS",
        "Preserve the CN/LOH inference-boundary guard until corrected inference is implemented and validated.",
    ),
    "ALG-021": row(
        "The audited workstation generator exits 3 and refuses output when one of its enumerated required metrics is absent.",
        "CONFIRMED", ["C073"], "NO_CHANGE_REQUIRED",
        ["README.md", "docs/wiki/Analysis-and-Presentation.md"],
        "Current docs scope refuse-on-missing to declared required fields.", "NOT_REQUIRED", [], "VERIFIED", [],
        "The result proves presence/non-null enforcement only.", [], "SOURCE_READY",
        "Retain the negative missing-required-field test.",
    ),
    "ALG-022": row(
        "Refuse-on-missing does not validate value truth, denominator correctness, source existence or scientific interpretation.",
        "CONFIRMED", ["C133"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "docs/wiki/Analysis-and-Presentation.md"],
        "Current docs explicitly deny the former structurally-impossible-fabrication claim.",
        "NOT_REQUIRED", [], "VERIFIED", [],
        "The counterexample establishes the validation boundary, not that published values were fabricated.", [],
        "SOURCE_READY", "Preserve source/value adversarial tests and keep the bounded wording.",
    ),
    "ALG-023": row(
        "The historical 334-source-reference/111-asset audit cannot be independently replayed from this repository because its inventories and receipt are absent.",
        "UNVERIFIED", ["C134", "C150"], "LOCAL_SOURCE_CORRECTED",
        ["docs/wiki/Home.md", "docs/wiki/System-Overview.md", "docs/explain/11_system-map-overview.standalone.html"],
        "Current public-source pages retain the historical self-report only with an adjacent ALG-023 UNVERIFIED warning.",
        "PASS", ["G_ALG023_HISTORICAL_AUDIT_UNVERIFIED"], "MISSING", [],
        "Absence of a receipt does not prove the historical audit did not happen; it prevents independent reproduction and current blanket guarantees.",
        ["No versioned 334/111 inventories, command log, source hashes or replayable result table were found."],
        "BLOCKED", "Publish commit-pinned 334/111 inventories and a replay receipt, or remove the counts from current public surfaces.",
    ),
    "ALG-024": row(
        "InterSubMod does not write BAM; LongLineage private candidate b9aaa12 contains longlineage-tag-bam, while the private baseline/main snapshot 5daf50f does not.",
        "CONFIRMED_WITH_LIMITS", ["C043", "C072"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/System-Overview.md", "docs/wiki/Upstream-and-Data.md", "docs/explain/11_system-map-overview.standalone.html", "docs/explain/14_upstream-data.standalone.html"],
        "Current source uses revision-scoped private-preview terminology.",
        "PASS", ["P0:C043"], "VERIFIED_WITH_LIMITS", [],
        "Feature presence at b9aaa12 is not a public-main or production capability.",
        ["LongLineage remains NOT_READY/non-production with P3/P4/P5/P7/P8 blocked."],
        "SOURCE_READY_WITH_LIMITS", "Keep immutable revision pins and verify live public visibility separately.",
    ),
    "ALG-025": row(
        "InterSubMod reads HP/PS from BAM aux tags; LongLineage/exact-PS use sidecars for candidate reconstruction. These are parallel provenance branches, not a direct engine-to-engine runtime interface.",
        "CONFIRMED_WITH_LIMITS", ["C044", "C049", "C069", "C072"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/System-Overview.md", "docs/wiki/Upstream-and-Data.md", "docs/explain/11_system-map-overview.standalone.html", "docs/explain/14_upstream-data.standalone.html"],
        "Current source separates the BAM-tag path from sidecar-based reconstruction and denies an automatic bridge.",
        "PASS", ["G_HPPS_REPRESENTATION_BOUNDARY"], "SOURCE_INSPECTION_ONLY", [],
        "Source inspection establishes representations and consumers, not an end-to-end supported cross-engine workflow.",
        ["No tracked top-level script automatically invokes both engines end to end; tagged-BAM capability is private-preview only."],
        "SOURCE_READY_WITH_LIMITS", "Preserve the boundary guard and add an E2E bridge receipt only if such a workflow becomes supported.",
    ),
    "ALG-026": row(
        "The formal production LongLineage dataset-gate endpoint has zero topology units; separate exact-PS research artifacts do not change that endpoint verdict.",
        "CONFIRMED_WITH_LIMITS", ["C031", "C039"], "NO_CHANGE_REQUIRED",
        ["docs/wiki/LongLineage-Engine.md", "docs/handoff/20260813_完整研究資料與軟體交接_01/evidence/longlineage_capability_matrix.md"],
        "Current handoff separates formal production gates from exploratory exact-PS evidence.",
        "MISSING", [], "VERIFIED_WITH_LIMITS", [],
        "Zero is scoped to the formal dataset-gate output, not all topology evidence.",
        ["LongLineage preview phases P3/P4/P5/P7/P8 remain blocked."], "SOURCE_READY_WITH_LIMITS",
        "Retain endpoint- and revision-specific wording.",
    ),
    "ALG-027": row(
        "Seven of seven frozen technical datasets have strict read-linkage sidecars with technical_all_pass=true, while validation_evidence_eligible=false.",
        "CONFIRMED_WITH_LIMITS", ["C025", "C045"], "NO_CHANGE_REQUIRED",
        ["README.md", "README.zh-TW.md", "docs/wiki/System-Overview.md", "docs/explain/11_system-map-overview.standalone.html"],
        "Current docs distinguish technical completion from biological eligibility.",
        "MISSING", [], "VERIFIED_WITH_LIMITS", [],
        "Technical receipt success is not biological validation eligibility.",
        ["The claim is frozen-snapshot and cohort-receipt scoped."], "SOURCE_READY_WITH_LIMITS",
        "Keep exact cohort receipt paths and both machine fields adjacent.",
    ),
    "ALG-028": row(
        "The historical LongLineage run_sample.sh invocation omitted required partition-root, topology and sidecar inputs and is not a supported executable example.",
        "CONFIRMED_WITH_LIMITS", ["C149"], "LOCAL_SOURCE_CORRECTED", ["pipeline/README.md"],
        "The historical pipeline is now marked SUPERSEDED/HISTORICAL and no longer offers that command as an execution entry.",
        "PASS", ["P0:C149", "G_HISTORICAL_PIPELINE_REDIRECT"], "HISTORICAL_ONLY", [],
        "The negative execution result pertains to the inspected private feature script.",
        ["A future public script can have a different interface."], "SOURCE_READY_WITH_LIMITS",
        "Do not restore an example until a public-clone smoke receipt covers all required arguments.",
    ),
    "ALG-029": row(
        "In the inspected private feature script, the default merged output is HCC1395.lineage_tagged.bam; per-chromosome outputs require --split-by-chrom.",
        "CONFIRMED_WITH_LIMITS", [], "LOCAL_SOURCE_CORRECTED", ["pipeline/README.md"],
        "The superseded historical pipeline no longer promises a per-chromosome default output.",
        "MISSING", [], "SOURCE_INSPECTION_ONLY", [],
        "Output naming is private-feature behavior and not a public production contract.",
        ["LongLineage preview is NOT_READY and production phases remain blocked."], "SOURCE_READY_WITH_LIMITS",
        "Document output naming only in the revision-pinned private preview until publication.",
    ),
    "ALG-030": row(
        "Historical LongLineage documentation records whole-genome tagged-BAM write counts, but the immutable raw receipt/log is not present in this handoff.",
        "UNVERIFIED", [], "NOT_APPLICABLE", [],
        "The counts are not promoted into current InterSubMod public-source claims.",
        "MISSING", [], "HISTORICAL_ONLY", [],
        "A secondary documentation record is insufficient to independently validate the exact counts.",
        ["Generating commit and immutable raw machine log/hash are missing."], "BLOCKED",
        "Publish the raw receipt/log hashes and exact producer commit before citing the counts as verified.",
    ),
    "ALG-031": row(
        "Current outputs establish read clustering, region stratification and bounded association, not cellular clones, clone ancestry or cellular prevalence.",
        "CONFIRMED", ["C024", "C057", "C070"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/System-Overview.md", "docs/handoff/20260813_完整研究資料與軟體交接_01/20260813_研究結論時間與Finality_01.md"],
        "All primary handoff entry paths enforce zero confirmed cellular subclones and zero confirmed linear ancestry.",
        "PASS", ["G_CONFIRMED_CLAIM_CEILING_ZERO"], "VERIFIED", [],
        "The claim is a scientific ceiling, not evidence that true subclones do not exist.", [], "SOURCE_READY",
        "Do not promote molecule-level or model-selected results without independent cellular truth.",
    ),
    "ALG-032": row(
        "technical_all_pass belongs to canonical/cohort_receipt.json and validation_evidence_eligible belongs to summary/all7_summary.json; neither is a top-level authority-manifest field.",
        "CONFIRMED", ["C025"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "docs/wiki/System-Overview.md", "docs/explain/11_system-map-overview.standalone.html"],
        "Current docs name the correct containers and deny authority-manifest attribution.",
        "NOT_REQUIRED", [], "VERIFIED", [],
        "The values and paths are frozen-snapshot scoped.", [], "SOURCE_READY",
        "Keep container paths adjacent to the fields in all future summaries.",
    ),
    "ALG-033": row(
        "The assessed InterSubMod baseline builds with CMake 3.14, C++17, HTSlib and declared Python dependencies on the validated bip7 environment.",
        "CONFIRMED_WITH_LIMITS", ["C051"], "NO_CHANGE_REQUIRED",
        ["README.md", "QUICKSTART.md", "docs/wiki/How-to-Run.md"],
        "Current docs state dependencies and separate machine-specific validation from portability claims.",
        "MISSING", [], "VERIFIED_WITH_LIMITS", [],
        "This is one-machine evidence, not a complete dependency-range or cross-platform matrix.",
        ["bip8 fresh-clone validation remains separately gated."], "SOURCE_READY_WITH_LIMITS",
        "Publish exact tool hashes and bip8 receipt before claiming dual-machine portability.",
    ),
    "ALG-034": row(
        "The exact-commit clean InterSubMod test suite must report exit 0 and zero failures; test and suite counts are generated dynamically rather than hard-coded.",
        "CONFIRMED_WITH_LIMITS", ["C145"], "LOCAL_SOURCE_CORRECTED",
        ["README.md", "README.zh-TW.md", "QUICKSTART.md", "docs/wiki/How-to-Run.md"],
        "Current docs removed the stale fixed 265/38 and 270/39 counts.",
        "PASS", ["G_README_DYNAMIC_TEST_COUNTS"], "VERIFIED_WITH_LIMITS",
        ["research/20260813_complete_research_handoff/receipts/bip7_portable_validation_587aeb67.json"],
        "A passing software suite does not validate biological claims and must be bound to the exact tested commit.",
        ["Hosted CI and bip8 remain separate gates."], "SOURCE_READY_WITH_LIMITS",
        "Regenerate the count and receipt on the final immutable release candidate.",
    ),
    "ALG-035": row(
        "The tracked tiny synthetic DEMO completes a one-SNV run and validates a 199-column schema plus read-level tree semantics.",
        "CONFIRMED_WITH_LIMITS", ["C059"], "LOCAL_SOURCE_CORRECTED",
        ["QUICKSTART.md", "docs/handoff/20260813_完整研究資料與軟體交接_01/20260813_bip7_bip8操作與驗證_01.md"],
        "The public run path is explicitly labeled DEMO and schema-only evidence.",
        "PASS", ["G_FIXTURE_DOCUMENTED_AS_TRACKED_DEMO", "G_SCHEMA_CURRENT_199"], "VERIFIED_WITH_LIMITS",
        ["research/20260813_complete_research_handoff/receipts/bip7_tiny_e2e_fb806d9b.json"],
        "The one-region synthetic run is not a cohort performance or biological benchmark.",
        ["Fixture inputs are synthetic and intentionally tiny."], "SOURCE_READY_WITH_LIMITS",
        "Replay on the final candidate and keep DEMO labeling adjacent to every result.",
    ),
}


def build(assessed_git_commit: str = ASSESSED_GIT_COMMIT) -> dict[str, object]:
    if not COMMIT_RE.fullmatch(assessed_git_commit):
        raise ValueError("assessed_git_commit must be a full lowercase Git SHA")
    commit_probe = subprocess.run(
        ["git", "-C", str(ROOT), "cat-file", "-e", f"{assessed_git_commit}^{{commit}}"],
        capture_output=True,
        text=True,
    )
    if commit_probe.returncode != 0:
        raise ValueError(f"assessed_git_commit is unavailable: {assessed_git_commit}")

    matrix_path = ROOT / MATRIX_RELATIVE
    claim_inventory_path = ROOT / CLAIM_INVENTORY_RELATIVE
    guard_path = ROOT / GUARD_RELATIVE
    with matrix_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError(f"unexpected matrix columns: {reader.fieldnames!r}")
        matrix_rows = list(reader)
    expected_ids = [f"ALG-{index:03d}" for index in range(1, 36)]
    matrix_ids = [source["claim_id"] for source in matrix_rows]
    if matrix_ids != expected_ids or sorted(POLICY) != expected_ids:
        raise ValueError("source matrix and policy must both cover ordered ALG-001..ALG-035 exactly once")
    with claim_inventory_path.open(encoding="utf-8", newline="") as handle:
        claim_inventory_rows = list(csv.DictReader(handle, delimiter="\t"))
    claim_inventory_ids = [source["claim_id"] for source in claim_inventory_rows]
    if claim_inventory_ids != [f"C{index:03d}" for index in range(1, 159)]:
        raise ValueError("claim inventory must contain ordered C001..C158 exactly once")

    records: list[dict[str, object]] = []
    snapshot_paths: set[str] = set()
    for source in matrix_rows:
        policy = POLICY[source["claim_id"]]
        unknown_related = set(policy["related_claim_ids"]) - set(claim_inventory_ids)
        if unknown_related:
            raise ValueError(f"{source['claim_id']} has unknown related claim IDs: {sorted(unknown_related)!r}")
        snapshot_paths.update(policy["source_correction"]["target_paths"])
        evidence_refs = [source["source_symbol"], source["runtime_or_test"], *policy["evidence"]["references"]]
        evidence_refs = list(dict.fromkeys(reference for reference in evidence_refs if reference))
        record = {
            "algorithm_claim_id": source["claim_id"],
            "source_row_sha256": canonical_row_hash(source),
            "audit_verdict": source["verdict"],
            "original_public_claim": source["public_claim"],
            **policy,
        }
        record["evidence"] = {**policy["evidence"], "references": evidence_refs}
        records.append(record)

    for path in snapshot_paths:
        if not (ROOT / path).is_file():
            raise FileNotFoundError(f"source snapshot path is missing: {path}")
    source_snapshot = [
        {"path": path, "sha256": sha256_file(ROOT / path)} for path in sorted(snapshot_paths)
    ]
    disposition_counts = Counter(record["current_disposition"] for record in records)
    correction_counts = Counter(record["source_correction"]["status"] for record in records)
    guard_counts = Counter(record["guard"]["status"] for record in records)
    evidence_counts = Counter(record["evidence"]["status"] for record in records)
    commit = assessed_git_commit
    matrix_binding = {"path": MATRIX_RELATIVE, "sha256": sha256_file(matrix_path), "rows": len(matrix_rows)}
    claim_inventory_binding = {
        "path": CLAIM_INVENTORY_RELATIVE,
        "sha256": sha256_file(claim_inventory_path),
        "rows": len(claim_inventory_rows),
    }
    guard_binding = {"path": GUARD_RELATIVE, "sha256": sha256_file(guard_path)}
    commit_pinned = all(
        commit_path_sha256(commit, binding["path"]) == binding["sha256"]
        for binding in [matrix_binding, claim_inventory_binding, guard_binding, *source_snapshot]
    )
    return {
        "schema_name": "intersubmod.algorithm_cli_claim_crosswalk",
        "schema_version": "1.0.0",
        "crosswalk_id": "intersubmod-algorithm-cli-claims-20260813",
        "assessed_at": "2026-08-13T20:00:00+08:00",
        "assessed_git_commit": commit,
        "source_state": "COMMIT_PINNED" if commit_pinned else "WORKTREE_HASH_BOUND_PENDING_COMMIT",
        "source_matrix": matrix_binding,
        "claim_inventory": claim_inventory_binding,
        "guard_registry": guard_binding,
        "source_snapshot": source_snapshot,
        "source_snapshot_tree_sha256": snapshot_tree_hash(source_snapshot),
        "disposition_semantics": {
            "applies_to": "bounded_current_claim_not_original_public_claim",
            "CONFIRMED": "The bounded current claim has direct verified evidence; this never revives a contradicted original claim.",
            "CONFIRMED_WITH_LIMITS": "The bounded current claim is supported only within its explicit revision, path, fixture or evidence ceiling.",
            "UNVERIFIED": "The claim lacks a replayable or sufficiently direct receipt and remains fail-closed.",
        },
        "counts": {
            "records": len(records),
            "related_claim_links": sum(len(record["related_claim_ids"]) for record in records),
            "records_with_related_claim_ids": sum(bool(record["related_claim_ids"]) for record in records),
            "records_without_related_claim_ids": sum(not record["related_claim_ids"] for record in records),
            "by_current_disposition": dict(sorted(disposition_counts.items())),
            "by_source_correction": dict(sorted(correction_counts.items())),
            "by_guard_status": dict(sorted(guard_counts.items())),
            "by_evidence_status": dict(sorted(evidence_counts.items())),
        },
        "gates": {
            "CROSSWALK_COMPLETENESS": {
                "status": "PASS",
                "reason": "ALG-001 through ALG-035 are present exactly once and hash-bound to the 2026-08-12 source rows."
            },
            "PUBLICATION_ASSET_READY": {
                "status": "PASS",
                "reason": "The crosswalk is publishable as a bounded, hash-bound inventory; this does not prove live default-branch, Wiki or Pages publication/refetch."
            },
            "RELEASE_ASSET_READY": {
                "status": "PASS",
                "reason": "The crosswalk is acceptable as a research-snapshot asset because ALG-023 and ALG-030 remain conspicuously UNVERIFIED and row-level BLOCKED; this is not aggregate repository release approval."
            },
        },
        "records": records,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--assessed-git-commit",
        default=ASSESSED_GIT_COMMIT,
        help="immutable source commit whose bound inputs are being assessed",
    )
    arguments = parser.parse_args(argv)
    output = build(arguments.assessed_git_commit)
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT.relative_to(ROOT)),
        "sha256": sha256_file(OUTPUT),
        "records": output["counts"]["records"],
        "by_current_disposition": output["counts"]["by_current_disposition"],
        "source_state": output["source_state"],
        "release_asset_ready": output["gates"]["RELEASE_ASSET_READY"]["status"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
