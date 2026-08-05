#!/usr/bin/env python3
"""Append the all-7 strict topology completion-boundary correction once."""

from __future__ import annotations

import json
from pathlib import Path


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
LEDGER = REPO / "research/autoresearch/evidence_ledger.jsonl"
CYCLE_ID = "20260723_exact_ps_strict_read_linkage_all7_completion_boundary_r1"
WORDING_CYCLE_ID = (
    "20260723_exact_ps_strict_read_linkage_all7_completion_boundary_r2"
)

ENTRY = {
    "cycle_id": CYCLE_ID,
    "hypothesis_id": (
        "STRICT_W_DOES_NOT_IMPLY_PRODUCTION_DIRECTED_TOPOLOGY_OR_CELLULAR_CLONE_TREE"
    ),
    "hypothesis": (
        "在所有7個technical datasets重新盤點最新receipts後，是否可將exact-PS strict "
        "read-linkage完成狀態與directed topology、clone count、exact parent-child及"
        "fusion tree狀態明確分層"
    ),
    "pipeline_track": "exact_ps_hp_strict_endpoint_read_linkage",
    "type": "comprehensive_completion_boundary_audit_and_html_correction",
    "tier": 3,
    "tier_used": 3,
    "task_type": "B_comprehensive_validation",
    "goals": ["G1", "G3", "G4", "G5"],
    "datasets_tested": [
        "HCC1395",
        "HCC1395_DORADO",
        "COLO829",
        "H1437",
        "H2009",
        "HCC1937",
        "HCC1954",
    ],
    "scale": (
        "7_technical_datasets_6_biological_cell_lines_chr1_22_"
        "154_dataset_chromosome_records"
    ),
    "corrects": "20260723_exact_ps_strict_read_linkage_all7",
    "verdict": (
        "L1_STRICT_READ_LINKAGE_COMPLETE_7_OF_7__"
        "PRODUCTION_STRICT_DIRECTED_TOPOLOGY_COMPLETE_0_OF_7__"
        "CELLULAR_CLONE_PARENT_CHILD_FUSION_VALIDATED_0_OF_7"
    ),
    "human_decision": (
        "report_completion_layers_separately;do_not_promote_undirected_W_or_"
        "legacy_candidate_trees_to_current_directed_topology_or_cellular_clone_tree"
    ),
    "key_observations": (
        "154/154 extraction receipts,154/154 strict receipts and7/7 summaries remain "
        "PASS with no numerical drift: S=469849,active loci=342374,memberships=613480,"
        "components=255752,k1=170131,W=85621,edges=1197530. Bounded inventory found no "
        "full-scope eligible v4 strict topology receipt. Latest HCC1395 v4 receipt stops "
        "at partition and topology=null. Earlier HCC1395 exact-PS topology is a non-binding "
        "technical pilot whose upstream receipt is exploratory/PARTIAL and "
        "validation_evidence_eligible=false; it is not bound to current strict W. DORADO "
        "and the other five datasets have no production strict topology receipt. Legacy "
        "all-7 candidate census uses coordinate/50kb-era grouping and contains5623 mixed-PS "
        "regions. Therefore current strict directed topology=0/7; clone count,exact cellular "
        "parent-child and cross-HP/cross-technical fused trees=0/7."
    ),
    "caveats": (
        "NOT_RUN is an absence claim bounded to the two declared synthesis search roots and "
        "v4 strict receipt schema at audit time. Endpoint edges are undirected read-linkage "
        "evidence,not evolutionary parent-child edges. HCC1395 and HCC1395_DORADO are "
        "technical products of one biological cell line. Earlier132 chromosome artifacts "
        "and HCC1954 use different builder SHA with data-specific no-trigger equivalence only."
    ),
    "artifacts_path": (
        "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01/"
    ),
    "html_report": (
        "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01.html"
    ),
    "topology_completion_audit": (
        "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01/topology_status_audit/"
        "all7_topology_completion_audit.json"
    ),
    "topology_completion_tsv": (
        "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01/data/"
        "topology_completion_status.tsv"
    ),
    "ready_receipt": (
        "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01/READY.json"
    ),
    "output_html_sha256": (
        "39fbfc71bee0b8123a9322d74c44866fbc398de714144f57941f494049922b81"
    ),
    "data_sha256": (
        "7c8fdaf0ae57a55c60b166ae4a34c8322de5efbd09e6e2483550b134d15b56f0"
    ),
    "artifact_sha256": (
        "80e66bfa300bb3533dd0d40e3bf05ac59b335dc42671ec63ecf62b5c5d868028"
    ),
    "ready_sha256": (
        "e7e8769cb80a849af45cbe3edeb4ca991c2802a29aef66046a0942cf27ca9942"
    ),
    "topology_audit_sha256": (
        "84f23c68f6225474bca7106f2d3100a8601b85c7c94cdaf1febfb3b2499aeaf3"
    ),
    "verification": (
        "25/25 report-data checks;Python32 tests;C++ strict graph PASS;artifact "
        "validator ok=true with18 datasets/19 sources;portable renderer11 charts/"
        "22 SVG variants/8 tables;official1440/390 QA and independent JS on/off "
        "desktop/mobile QA all PASS with0px overflow"
    ),
    "reviewer": (
        "root_plus_three_read_only_inventory_agents_plus_independent_final_bundle_"
        "agent_plus_data_analytics_artifact_validator_plus_official_browser_verifier_"
        "plus_independent_playwright_visual_QA"
    ),
    "next_action": (
        "Before claiming any directed topology, execute the v4 strict topology stage "
        "on a preregistered full scope and publish eligible per-sample receipts. Then "
        "preregister nonduplicated VAF/read-pattern likelihood ranking and separately "
        "validate cellular clone count,parent-child lineage and cross-HP fusion with "
        "CN/purity/CCF/methylation or orthogonal truth."
    ),
    "identified_issues": [
        "strict topology production not run",
        "HCC1395 pilot is partial and non-binding",
        "DORADO topology not run",
        "legacy all7 topology is incompatible with exact-PS strict W",
        "endpoint linkage is undirected",
        "clone count and cellular parent-child unvalidated",
        "cross-HP fusion unvalidated",
        "mixed builder SHA has data-specific equivalence only",
    ],
    "operator": "codex-gpt5+four_collaborating_agents",
    "timestamp": "2026-07-23T19:05:29+08:00",
}

WORDING_ENTRY = {
    "cycle_id": WORDING_CYCLE_ID,
    "hypothesis_id": "STRICT_TOPOLOGY_COMPLETION_REPORTING_UNIT_WORDING_CORRECTION",
    "hypothesis": (
        "完成層級HTML是否應將舊candidate census精確寫成7 technical datasets / "
        "6 biological cell lines，而非容易誤讀的7-sample"
    ),
    "pipeline_track": "exact_ps_hp_strict_endpoint_read_linkage",
    "type": "delivery_wording_correction_with_full_bundle_revalidation",
    "tier": 3,
    "tier_used": 3,
    "task_type": "B_comprehensive_validation",
    "goals": ["G4", "G5"],
    "corrects": CYCLE_ID,
    "verdict": (
        "PASS_7_DATASET_6_CELL_LINE_WORDING_CORRECTED__"
        "COMPLETION_COUNTS_AND_ALL_NUMERICAL_RESULTS_UNCHANGED"
    ),
    "human_decision": (
        "use_7_dataset_6_cell_line_for_technical_product_scope;"
        "retain_L1_7_of_7_and_strict_topology_clone_fusion_0_of_7"
    ),
    "key_observations": (
        "Independent final bundle audit passed all identities,numbers,status rows and "
        "browser checks,with one non-blocking wording finding: legacy '7-sample' could "
        "misstate7 technical datasets as7 independent biological samples. HTML was "
        "corrected to '7-dataset / 6-cell-line candidate-tree census' and fully rebuilt. "
        "All analytical data are unchanged: W=85621;L1=7/7;production strict directed "
        "topology=0/7;clone,parent-child,fusion=0/7."
    ),
    "caveats": (
        "This correction changes reader-facing wording and artifact/HTML/READY hashes "
        "only; it does not add topology execution or biological validation."
    ),
    "artifacts_path": (
        "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01/"
    ),
    "html_report": (
        "InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01/"
        "20260723_exactPS嚴格ReadLinkage全資料集報告_01.html"
    ),
    "output_html_sha256": (
        "279795f5463063266921a18bd4c3885bf2ecd468a40d0a3fb1f07ca6f5feb516"
    ),
    "data_sha256": (
        "7c8fdaf0ae57a55c60b166ae4a34c8322de5efbd09e6e2483550b134d15b56f0"
    ),
    "artifact_sha256": (
        "042ada37bdf4e0543a173c9f4972fc931869a22e8118169d709fa62e4ab1c639"
    ),
    "ready_sha256": (
        "1c95ed3c6509e9477b2ab38dc380d16e8072e0ba0ce8e17f82842dee7a0dd2ee"
    ),
    "verification": (
        "artifact validator ok=true;portable renderer11 charts/22 SVG variants/"
        "8 tables;official1440/390 and independent JS on/off browser QA all PASS;"
        "READY sidecar and bundle identities PASS"
    ),
    "reviewer": "independent_final_completion_bundle_audit_plus_root_revalidation",
    "next_action": (
        "Run full-scope v4 strict topology before making any directed topology or "
        "cellular clone/parent/fusion claim."
    ),
    "identified_issues": [
        "reader-facing technical-dataset versus biological-sample terminology",
        "strict directed topology remains not run",
    ],
    "operator": "codex-gpt5+independent_final_bundle_agent",
    "timestamp": "2026-07-23T19:11:00+08:00",
}


def main() -> int:
    entries = [
        json.loads(line)
        for line in LEDGER.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    existing = {entry.get("cycle_id") for entry in entries}
    pending = [
        entry
        for entry in (ENTRY, WORDING_ENTRY)
        if entry["cycle_id"] not in existing
    ]
    if not pending:
        print(f"already_present={CYCLE_ID},{WORDING_CYCLE_ID}")
        return 0
    with LEDGER.open("a", encoding="utf-8") as handle:
        for entry in pending:
            handle.write(
                json.dumps(entry, ensure_ascii=False, separators=(",", ":"))
                + "\n"
            )
    print(
        f"appended={','.join(entry['cycle_id'] for entry in pending)} "
        f"ledger={LEDGER}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
