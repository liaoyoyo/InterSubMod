#!/usr/bin/env python3
"""Validate the final HCC1395 coarse-topology report delivery."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
TOPIC = ROOT / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation"
REPORT = ROOT / (
    "docs/reports/in_progress/2026/07/"
    "20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01"
)
STEM = "20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact_path = REPORT / "artifact.json"
    html_path = REPORT / f"{STEM}.html"
    receipt_path = REPORT / "qa_receipt.json"
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    analysis = json.loads((TOPIC / "data/topology_pair_analysis.json").read_text(encoding="utf-8"))
    annotation = read_tsv(TOPIC / "data/hcc1395_annotation_reproducibility.tsv")
    summary = read_tsv(TOPIC / "data/coarse_topology_all_dataset_summary.tsv")
    vaf_metrics = read_tsv(TOPIC / "data/hcc1395_vaf_pair_metrics.tsv")
    vaf_summary = json.loads((TOPIC / "data/hcc1395_vaf_pair_summary.json").read_text(encoding="utf-8"))
    vaf_checks = read_tsv(TOPIC / "data/hcc1395_vaf_pair_checks.tsv")
    html = html_path.read_text(encoding="utf-8")

    checks: list[tuple[str, bool]] = []
    checks.append(("seven dataset rows", len(summary) == 7))
    checks.append(("47,377 primary regions", sum(int(row["primary_regions"]) for row in summary) == 47_377))
    checks.append(("39,885 complete regions", sum(int(row["complete_regions"]) for row in summary) == 39_885))
    checks.append(("five-class conservation", all(
        sum(int(row[field]) for field in [
            "no_within_hp_relation", "sister_only", "direct_only",
            "sister_and_direct", "topology_multiple_unresolved",
        ]) == int(row["complete_regions"])
        for row in summary
    )))
    checks.append(("topology checks 64/64", analysis["validation"]["checks"] == 64 == analysis["validation"]["passed"]))
    exact = next(row for row in analysis["hcc1395_pair_metrics"] if row["scenario"] == "exact_coordinate")
    checks.extend([
        ("exact matches 6,252", exact["matched_all"] == 6_252),
        ("complete-both 5,720", exact["complete_both"] == 5_720),
        ("agreements 3,969", round(exact["raw_agreement"] * exact["complete_both"]) == 3_969),
        ("coarse agreement", abs(exact["raw_agreement"] - 3_969 / 5_720) < 1e-12),
        ("kappa 0.497329", abs(exact["cohen_kappa"] - 0.4973289) < 1e-6),
        ("permutation p=1/5001", abs(exact["permutation_null"]["agreement_p_ge"] - 1 / 5_001) < 1e-12),
    ])
    annotation_features = [row for row in annotation if row["feature"] != "ALL"]
    checks.append(("seven annotation features", len(annotation_features) == 7))
    checks.append(("all annotation p>0.05", all(float(row["permutation_p_two_sided"]) > 0.05 for row in annotation_features)))

    vaf_lookup = {
        (row["metric_family"], row["population"], row["comparison"]): row
        for row in vaf_metrics
    }
    shape_swap = vaf_lookup[("selected_single_shape", "both_shapes_require_vaf_from_original_Topo_gt1", "HP1_HP2_swap_tolerant")]
    exact_swap = vaf_lookup[("vaf_unique_exact_first_candidate", "both_sides_actually_use_vaf_ranking_and_unique", "HP1_HP2_swap_tolerant")]
    checks.extend([
        ("VAF checks 20/20", len(vaf_checks) == 20 and all(row["pass"] == "True" for row in vaf_checks)),
        ("VAF summary 20/20", vaf_summary["validation"]["passed"] == vaf_summary["validation"]["total"] == 20),
        ("HCC shape-selectable 6,798", vaf_summary["sample_conservation"]["HCC1395"]["structural_or_vaf_shape_selectable"] == 6_798),
        ("DORADO shape-selectable 6,082", vaf_summary["sample_conservation"]["HCC1395_DORADO"]["structural_or_vaf_shape_selectable"] == 6_082),
        ("VAF shape-rescue swap 1,624/2,089", int(shape_swap["numerator"]) == 1_624 and int(shape_swap["denominator"]) == 2_089),
        ("VAF unique-exact swap 949/2,543", int(exact_swap["numerator"]) == 949 and int(exact_swap["denominator"]) == 2_543),
    ])

    manifest = artifact["manifest"]
    artifact_text = json.dumps(artifact, ensure_ascii=False)
    effective_blocks = sum(
        len(block.get("cardIds", [])) if block["type"] == "metric-strip" else 1
        for block in manifest["blocks"]
    )
    checks.extend([
        ("49 rendered artifact blocks", effective_blocks == 49),
        ("six artifact charts", len(manifest["charts"]) == 6),
        ("19 artifact tables", len(manifest["tables"]) == 19),
        ("four metric cards", len(manifest["cards"]) == 4),
        ("partial status", artifact["snapshot"]["status"] == "partial"),
        ("no false source-query blockers", artifact["snapshot"]["accessIssues"] == []),
        ("header compatibility style", "hcc1395-portable-header-scrollbar-compat" in html),
        ("correct TVD text", "20.03 percentage points" in artifact_text and "0.20 pp" not in artifact_text),
        ("same biological sample caveat", "同一生物 cell line" in html),
        ("VAF inference caveat", "VAF-supported 推測" in artifact_text and "不是 posterior" in artifact_text),
        ("VAF shape result", "77.74%" in artifact_text),
        ("VAF exact result", "37.32%" in artifact_text),
        ("scientific no-go", "SCIENTIFIC NO-GO" in html),
        ("latest freeze cutoff", "2026-07-12T02:23:35+08:00" in artifact_text),
        ("producer seven of seven", "clean-v3 producer 7/7 PASS" in artifact_text),
        ("closeout failure disclosed", "E_SIDECAR_VALIDATION" in artifact_text),
        ("continuation failure disclosed", "continuation=FAILED" in artifact_text),
        ("superseded six-of-seven absent", "clean-v3 producer 6/7 PASS" not in artifact_text),
    ])
    for canonical_name in ["無 HP 內關係", "姐妹 only", "直系 only", "姐妹＋直系", "Topo>1 未定"]:
        checks.append((f"canonical class: {canonical_name}", canonical_name in artifact_text))

    for key in [
        "html", "artifact", "markdown", "analysis_json", "annotation_json",
        "vaf_pair_metrics", "vaf_pair_summary", "vaf_pair_regions", "vaf_pair_checks",
        "latest_pipeline_status", "qa_desktop", "qa_mobile",
    ]:
        record = receipt["outputs"][key]
        path = ROOT.parent / record["path"]
        checks.append((f"receipt bytes: {key}", path.stat().st_size == record["bytes"]))
        checks.append((f"receipt sha256: {key}", sha256(path) == record["sha256"]))

    failed = [label for label, passed in checks if not passed]
    print(f"INPUT artifact -> {artifact_path.relative_to(ROOT)}")
    print(f"INPUT analysis -> {(TOPIC / 'data/topology_pair_analysis.json').relative_to(ROOT)}")
    print(f"INPUT annotation -> {(TOPIC / 'data/hcc1395_annotation_reproducibility.tsv').relative_to(ROOT)}")
    print(f"INPUT VAF pair -> {(TOPIC / 'data/hcc1395_vaf_pair_summary.json').relative_to(ROOT)}")
    print(f"INPUT receipt -> {receipt_path.relative_to(ROOT)}")
    print("OUTPUT -> no file; validation-only")
    print(f"CHECKS -> {len(checks) - len(failed)}/{len(checks)} PASS")
    if failed:
        print("FAIL -> " + " | ".join(failed))
        return 1
    print("RESULT -> PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
