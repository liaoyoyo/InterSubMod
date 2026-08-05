#!/usr/bin/env python3
"""Append the finalized singleton ALT methyl-substructure result to the ledger."""

from __future__ import annotations

import json
from pathlib import Path


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
LEDGER = REPO / "research/autoresearch/evidence_ledger.jsonl"
CYCLE_ID = "20260718_singleton_alt_methyl_substructure_validation"
CORRECTION_CYCLE_ID = "20260718_singleton_alt_methyl_substructure_validation_heatmap_r1"

ENTRY = {
    "cycle_id": CYCLE_ID,
    "hypothesis_id": "POSITIONAL_SINGLETON_FOCAL_ALT_METHYL_SUBSTRUCTURE_AND_CLONE_LABEL",
    "hypothesis": (
        "完整 HCC1395 positional-singleton sSNV 母體中，focal-ALT reads 是否存在可重現的甲基子結構，"
        "且現有證據是否足以標記特殊 cellular clone/subclone"
    ),
    "pipeline_track": "longphase_s_pass_focal_alt_read_level_epigenetic_heterogeneity",
    "type": "comprehensive_validation_with_portable_html_and_read_distance_heatmaps",
    "tier": 3,
    "tier_used": 3,
    "task_type": "B_comprehensive_validation",
    "goals": ["G3", "G4", "G5"],
    "datasets_tested": [
        "COLO829",
        "H1437",
        "H2009",
        "HCC1395",
        "HCC1395_DORADO",
        "HCC1937",
        "HCC1954",
    ],
    "scale": "7_datasets_50432_singleton_dataset_sites_HCC1395_8279_full_sites",
    "verdict": (
        "PASS_COMPLETE_SINGLETON_CENSUS_AND_TWO_RESIDUAL_EPIGENETIC_PARTITION_CANDIDATES_"
        "CELLULAR_CLONE_SUBCLONE_NO_GO"
    ),
    "human_decision": (
        "use_M1_734_as_screen_and_M2_2_as_orthogonal_followup_candidates_only;"
        "do_not_label_as_confirmed_cellular_clone_subclone_or_lineage"
    ),
    "key_observations": (
        "HCC1395 LongPhase-S recalibrated PASS autosomal biallelic sSNV=79687; 50kb positional components=16501, "
        "singleton sites=8279, multilocus components=8222/71408 sites. M1 evaluable=8074 and stable multigroup="
        "734/8279=8.8658% (734/8074=9.0909%). M2 measured-axis clear residual partition=2/8279=0.02416%="
        "24.16 per100k; M2 FAIL=0, NOT_EVALUABLE=732, NOT_RUN=7545. The two HCC1395 loci are "
        "chr14:86272476 A>T (86/22 ALT reads; between/within distance ratio1.8919; global high/low methyl partition) "
        "and chr22:47466517 A>G (88/21; ratio1.5668; CpG-pattern partition). Exact joins cover217/217 focal-ALT core "
        "reads; both pass HP/strand/geometry/MAPQ/CpG-called measured-axis guardrails. HCC1395 versus DORADO has "
        "7484 exact common singleton loci and M1 state agreement90.43%, but exact M2-PASS overlap=0. Seven-dataset "
        "singleton total=50432, M1 flags=5961, M2 PASS=30/FAIL=18/NOT_EVALUABLE=5913/NOT_RUN=44471. "
        "Confirmed cellular clone/subclone=0 because required biological gates were not run. Portable HTML has "
        "39 blocks,14 charts,8 tables,4 metrics; canonical and independent browser verification pass at1440/390."
    ),
    "caveats": (
        "8279 means 8279 component-size-1 dataset-sites, not 8279 variants in one region. M2 operational yield is "
        "not biological prevalence; NOT_EVALUABLE and NOT_RUN are not negatives. A singleton has no local second "
        "sSNV, so mutation co-membership, clone number, parent-child order and a unique evolutionary tree are "
        "unidentifiable. VAF is focal allele-burden context only and was not used to discover methyl groups or "
        "infer ancestry. G1/G2/formal R1, matched-normal methylation, CN/purity/CCF and cellular lineage are NOT_RUN. "
        "HCC1395 and DORADO are the same cell line but not independent biological replication, and exact M2-PASS "
        "replication is zero. Gene/drug annotations are prioritization context, not functional validation."
    ),
    "artifacts_path": "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/",
    "html_report": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/"
        "20260718_HCC1395_singleton_ALT甲基子結構驗證報告_01.html"
    ),
    "markdown_report": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/"
        "20260718_HCC1395_singleton_ALT甲基子結構驗證報告_01.md"
    ),
    "site_audit": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/"
        "hcc1395_singleton_site_audit.tsv.gz"
    ),
    "analysis_receipt": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/"
        "validation_receipt.json"
    ),
    "html_qa_receipt": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/"
        "html_qa_receipt.json"
    ),
    "output_html_sha256": "fccd7846dc2a38e2b015fe3dc86cd47cc01b69942161bdd46c3bda469c61784c",
    "output_markdown_sha256": "7c4758964075b341cdc68b8bbd6c6ce5305c1d30cba498fddc1bceee503346b6",
    "site_audit_sha256": "45be5c5f168ed4726c9c5c4134528f71b0e152697b5704eb13ed8eceaf39e783",
    "analysis_receipt_sha256": "ec93a24034ef9815529d409a7594d3f174c8746abe28efd115e3167bd2551962",
    "html_qa_receipt_sha256": "3c8f9b4465f9f9b01a21eea3478ac52e2830cb2d18b7c06f3859ef67b08dc5f1",
    "reviewer": (
        "root_full_recount_plus_exact_HP_read_matrix_auditor_plus_span_grid_auditor_plus_"
        "claim_redteam_plus_canonical_artifact_validator_plus_independent_browser_QA"
    ),
    "next_action": (
        "For the two M2 candidates, preregister matched-normal/tumor-REF methyl specificity and an independent "
        "second somatic marker or single-cell/colony/multi-region truth; add allele-specific CN, purity and CCF. "
        "Only after those gates may labels advance from residual epigenetic partition candidate to molecular "
        "haplotype or cellular clone/subclone."
    ),
    "identified_issues": [
        "singleton lacks local second genetic marker",
        "M2 operational yield is not prevalence",
        "most M1 flags are M2 NOT_EVALUABLE",
        "matched-normal and tumor-REF specificity not run",
        "CN purity CCF and lineage not run",
        "VAF cannot establish parent-child order",
        "no exact HCC1395/DORADO M2-PASS replication",
        "same-cell-line technical sources are not biological replication",
        "gene and drug context is not functional validation",
    ],
    "operator": "codex-gpt5+three_collaborating_agents",
    "timestamp": "2026-07-19T00:38:19+08:00",
}

CORRECTION_ENTRY = {
    "cycle_id": CORRECTION_CYCLE_ID,
    "hypothesis_id": "POSITIONAL_SINGLETON_HEATMAP_RENDERING_DELIVERY_CORRECTION",
    "hypothesis": (
        "正式 portable HTML 是否把兩個 HCC1395 M2 loci 的 distance、shared-CpG 與 methylation "
        "矩陣真正呈現為可見彩色 heatmap，而非 native renderer 的數字表 fallback"
    ),
    "pipeline_track": "longphase_s_pass_focal_alt_read_level_epigenetic_heterogeneity",
    "type": "delivery_correction_with_dom_and_visual_qa",
    "tier": 3,
    "tier_used": 3,
    "task_type": "B_comprehensive_validation",
    "goals": ["G3", "G5"],
    "corrects": CYCLE_ID,
    "verdict": "PASS_SIX_COLORED_HEATMAPS_RENDERED_AND_BROWSER_VERIFIED",
    "human_decision": (
        "use_the_r1_HTML_as_the_final_delivery;supersede_the_initial_native_heatmap_table_fallback"
    ),
    "key_observations": (
        "Independent red-team found that six native heatmap chart blocks were exported as numeric tables without "
        "colored cells. They were replaced with six schema-supported, script-free HTML heatmaps: distance, "
        "shared-CpG and methylation for each of chr14:86272476 and chr22:47466517. Each heatmap contains64 cells, "
        "explicit color scale, Group A/B boundaries, values, NA semantics and accessibility labels. Playwright DOM "
        "QA passes6/6 at1440 and390px; unique rendered colors range11-28, iframe geometry is nonzero with no clipping, "
        "and parent horizontal overflow is absent. Canonical and independent portable verifiers pass39 blocks, "
        "8 native charts,6 HTML heatmaps,8 tables and4 metrics. Analytical counts are unchanged: HCC1395 singleton "
        "8279, M1=734, M2 PASS=2, established cellular clone/subclone=0 because required gates were not run."
    ),
    "caveats": (
        "The heatmaps display four deterministic representative reads per group for readability; all reported "
        "distance effects and M2 decisions use all108/109 focal-ALT core reads. Colored rendering improves visual "
        "delivery but adds no biological validation and does not upgrade the clone/subclone claim ceiling."
    ),
    "artifacts_path": "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/",
    "html_report": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/"
        "20260718_HCC1395_singleton_ALT甲基子結構驗證報告_01.html"
    ),
    "analysis_receipt": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/"
        "validation_receipt.json"
    ),
    "html_qa_receipt": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/"
        "html_qa_receipt.json"
    ),
    "heatmap_rendering_qa": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/"
        "heatmap_rendering_qa.json"
    ),
    "heatmap_screenshots": (
        "InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/heatmap_qa/"
    ),
    "output_html_sha256": "13307dd73ee26856e28e699c764de7281cb304ac5fa022baa781da30b0fde3e9",
    "artifact_sha256": "41938091fced89ebb719c10361e5963463bf299788adec72983098a691f1b5bb",
    "site_audit_sha256": "cfa8d15c522f789411878f435b0827cb5074954b55a0f30342302017d12387ce",
    "analysis_receipt_sha256": "586181d0898640d8ec47f43e9b0d342080a4ad8a29c7d940656d131a27e999d8",
    "html_qa_receipt_sha256": "649546a5b8f173afede65b7b7eea69a3a7250ea3eb0f37c0357bae9105265a99",
    "heatmap_rendering_qa_sha256": "5345c5553c5e972e836ca1230d9c3b81658d269f2283693ac8193a03f0713265",
    "reviewer": (
        "root_plus_final_claim_redteam_plus_canonical_artifact_validator_plus_independent_portable_verifier_"
        "plus_playwright_frame_level_heatmap_DOM_QA_plus_visual_screenshot_inspection"
    ),
    "next_action": (
        "Use the corrected HTML for interpretation. Biological upgrade still requires G1/G2 genetic co-membership, "
        "matched-normal/tumor-REF methyl controls, CN/purity/CCF and orthogonal cellular truth."
    ),
    "identified_issues": [
        "native portable heatmap exported as numeric-table fallback",
        "generic chart-count QA did not establish colored heatmap rendering",
        "representative-read visualization is not the full analytical denominator",
        "visual correction does not validate cellular clone identity",
    ],
    "operator": "codex-gpt5+claim_redteam",
    "timestamp": "2026-07-19T00:51:57+08:00",
}


def main() -> int:
    lines = LEDGER.read_text(encoding="utf-8").splitlines()
    entries = [json.loads(line) for line in lines if line.strip()]
    existing = {entry.get("cycle_id") for entry in entries}
    pending = [entry for entry in (ENTRY, CORRECTION_ENTRY) if entry["cycle_id"] not in existing]
    if not pending:
        print(f"already_present={CYCLE_ID},{CORRECTION_CYCLE_ID}")
        return 0
    with LEDGER.open("a", encoding="utf-8") as handle:
        for entry in pending:
            handle.write(json.dumps(entry, ensure_ascii=False, separators=(",", ":")) + "\n")
    print(f"appended={','.join(entry['cycle_id'] for entry in pending)} ledger={LEDGER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
