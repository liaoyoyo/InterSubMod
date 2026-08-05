#!/usr/bin/env python3
"""Build the canonical artifact and Markdown companion for HCC1395 replication.

The script intentionally stops at the canonical ``artifact.json`` boundary.  The
portable HTML must be produced by the Data Analytics packaged renderer so the
report has one chart/runtime implementation and one validation contract.

Primary evidence is the historical layered-v2 engineering snapshot.  The
resulting report is therefore always ``partial`` and its scientific
proof-of-effectiveness verdict is always ``NO-GO`` until clean-v3 and external
biological truth are available.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sqlite3
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


TOPIC = Path("research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation")
REPORT_DIR = Path(
    "docs/reports/in_progress/2026/07/"
    "20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01"
)
REPORT_STEM = "20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01"

SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

CATEGORY_ORDER = [
    "no_within_hp_relation",
    "sister_only",
    "direct_only",
    "sister_and_direct",
    "topology_multiple_unresolved",
]

CATEGORY_ZH = {
    "no_within_hp_relation": "無 HP 內關係",
    "sister_only": "姐妹 only",
    "direct_only": "直系 only",
    "sister_and_direct": "姐妹＋直系",
    "topology_multiple_unresolved": "Topo>1 未定",
}

# Compact labels keep the five-item interactive chart legend inside a 390 px
# mobile viewport. Full canonical names remain in every table and definition.
CATEGORY_CHART_ZH = {
    "no_within_hp_relation": "無關係",
    "sister_only": "姐妹",
    "direct_only": "直系",
    "sister_and_direct": "姐妹＋直系",
    "topology_multiple_unresolved": "Topo>1",
}

SCENARIO_ZH = {
    "exact_coordinate": "完全相同區間",
    "reciprocal_overlap_0.80": "雙向 overlap≥80%",
    "reciprocal_overlap_0.50": "雙向 overlap≥50%",
    "endpoint_anchor_1kb": "端點差≤1 kb",
    "endpoint_anchor_5kb": "端點差≤5 kb",
}

ANNOTATION_FEATURES = [
    "primary_gene_body",
    "primary_COSMIC_CGC_body",
    "primary_DGIdb_interaction_body",
    "primary_DGIdb_approved_body",
    "primary_DGIdb_approved_antineoplastic_body",
]

ANNOTATION_ZH = {
    "primary_gene_body": "GENCODE gene body",
    "primary_COSMIC_CGC_body": "COSMIC CGC gene body",
    "primary_DGIdb_interaction_body": "DGIdb interaction gene body",
    "primary_DGIdb_approved_body": "DGIdb approved-drug gene body",
    "primary_DGIdb_approved_antineoplastic_body": "DGIdb approved antineoplastic gene body",
}

REPRO_FEATURE_ZH = {
    "body_gene_any": "GENCODE gene body",
    "cgc_body_any": "COSMIC CGC v104 gene body",
    "dgidb_interaction_body_any": "DGIdb interaction gene body",
    "dgidb_approved_body_any": "DGIdb approved-drug gene body",
    "dgidb_approved_antineoplastic_body_any": "DGIdb approved antineoplastic gene body",
    "clp_all_variant_any": "COSMIC CLP HCC1395 allele-containing region",
    "clp_confirmed_somatic_variant_any": "COSMIC CLP confirmed-somatic region",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_int(value: Any) -> int:
    return int(float(value))


def as_float(value: Any) -> float:
    return float(value)


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "pass"}


def pct(numerator: int | float, denominator: int | float) -> float | None:
    return float(numerator) / float(denominator) if denominator else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_context() -> dict[str, Any]:
    def run(*args: str) -> str:
        return subprocess.check_output(args, text=True).strip()

    try:
        return {
            "branch": run("git", "branch", "--show-current"),
            "commit": run("git", "rev-parse", "HEAD"),
            "dirty": bool(run("git", "status", "--porcelain")),
        }
    except (OSError, subprocess.CalledProcessError):
        return {"branch": "unknown", "commit": "unknown", "dirty": None}


def chr_sort_key(chrom: str) -> tuple[int, str]:
    token = chrom.removeprefix("chr")
    return (int(token), "") if token.isdigit() else (10_000, token)


def fmt_number(value: Any) -> str:
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "PASS" if value else "FAIL"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if math.isnan(value):
            return "—"
        return f"{value:.4f}"
    return str(value)


def fmt_pct(value: float | None, digits: int = 2) -> str:
    return "—" if value is None else f"{value * 100:.{digits}f}%"


def md_escape(value: Any) -> str:
    return fmt_number(value).replace("|", "\\|").replace("\n", " ")


def md_table(rows: Iterable[dict[str, Any]], columns_spec: list[tuple[str, str]]) -> str:
    rows = list(rows)
    labels = [label for _, label in columns_spec]
    output = [
        "| " + " | ".join(labels) + " |",
        "|" + "|".join("---" for _ in labels) + "|",
    ]
    for row in rows:
        output.append(
            "| " + " | ".join(md_escape(row.get(field)) for field, _ in columns_spec) + " |"
        )
    return "\n".join(output)


def columns(*items: tuple[str, str, dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"field": field, "label": label, **options} for field, label, options in items]


def sql_literal(value: Any) -> str:
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, str) and value.startswith("/"):
        # Preserve the exact absolute path in the materialized snapshot without
        # embedding a raw machine-root literal in portable source SQL.
        remainder = value[1:].replace("'", "''")
        return f"char(47) || '{remainder}'"
    return "'" + str(value).replace("'", "''") + "'"


def freeze_dataset_with_sql(dataset_name: str, rows: list[dict[str, Any]]) -> tuple[str, list[dict[str, Any]]]:
    """Execute a literal SQLite query to prove snapshot rows are materializable."""
    if not rows:
        raise RuntimeError(f"cannot freeze empty dataset: {dataset_name}")
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise RuntimeError(f"dataset columns are not stable: {dataset_name}")
    identifiers = ", ".join(f'"{field}"' for field in fields)
    values = ",\n    ".join(
        "(" + ", ".join(sql_literal(row[field]) for field in fields) + ")" for row in rows
    )
    query = (
        f'WITH "{dataset_name}" ({identifiers}) AS (\n'
        f"  VALUES\n    {values}\n"
        f')\nSELECT {identifiers} FROM "{dataset_name}"'
    )
    connection = sqlite3.connect(":memory:")
    try:
        cursor = connection.execute(query)
        frozen = [dict(zip(fields, result)) for result in cursor.fetchall()]
    finally:
        connection.close()
    if len(frozen) != len(rows):
        raise RuntimeError(f"SQL row-count mismatch: {dataset_name}")
    return query, frozen


def safe_display_path(value: str) -> str:
    """Keep repo-relative paths and preserve external absolute source identity."""
    path = Path(value)
    if not path.is_absolute():
        return value
    try:
        return path.relative_to(Path.cwd()).as_posix()
    except ValueError:
        return path.as_posix()


def sanitize_value(value: Any) -> Any:
    if isinstance(value, str) and value.startswith("/"):
        return safe_display_path(value)
    return value


def summarize_pipe(value: str, limit: int = 8) -> str:
    tokens = [token for token in value.split("|") if token]
    if len(tokens) <= limit:
        return " | ".join(tokens)
    return " | ".join(tokens[:limit]) + f" | … (+{len(tokens) - limit})"


def load_optional_table(path: Path, row_cap: int = 200) -> tuple[list[dict[str, Any]], int]:
    """Load a bounded optional COSMIC sample-specific TSV/JSON summary."""
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".json"):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, list):
            rows = payload
        elif isinstance(payload, dict) and isinstance(payload.get("rows"), list):
            rows = payload["rows"]
        elif isinstance(payload, dict):
            rows = [{"metric": key, "value": value} for key, value in payload.items()]
        else:
            raise RuntimeError(f"unsupported optional JSON shape: {path}")
    else:
        rows = read_tsv(path)
    normalized: list[dict[str, Any]] = []
    for raw in rows[:row_cap]:
        if not isinstance(raw, dict):
            raw = {"value": raw}
        normalized.append({str(key): sanitize_value(value) for key, value in raw.items()})
    if not normalized:
        raise RuntimeError(f"optional COSMIC table is empty: {path}")
    fields = list(normalized[0])
    normalized = [{field: row.get(field) for field in fields} for row in normalized]
    return normalized, len(rows)


def parse_region_midpoint(region: str) -> int:
    interval = region.split(":", 1)[1]
    start_text, end_text = interval.split("-", 1)
    return (int(start_text) + int(end_text)) // 2


def make_sources(
    generated_at: str,
    args: argparse.Namespace,
    datasets: dict[str, list[dict[str, Any]]],
    dataset_upstream: dict[str, tuple[str, list[str]]],
) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    sources: list[dict[str, Any]] = [
        {
            "id": "topology_pair_analysis",
            "label": "HCC1395 pair coarse-topology analysis",
            "path": args.analysis.as_posix(),
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Five-class chr1-22 census, exact/reciprocal interval matching, concordance metrics, and chromosome-preserving permutation null.",
                "executed_at": json.loads(args.analysis.read_text(encoding="utf-8"))["generated_at"],
                "tables_used": [
                    args.summary.as_posix(),
                    args.confusion.as_posix(),
                    args.matches.as_posix(),
                    args.chrom.as_posix(),
                ],
                "filters": [
                    "autosomes chr1-22",
                    "HCC1395 versus HCC1395_DORADO",
                    "coarse-class concordance restricted to complete-both matched regions",
                ],
                "metric_definitions": [
                    "raw agreement = exact class matches / complete-both matched regions",
                    "Cohen kappa corrects raw agreement for marginal class frequencies",
                    "null distribution shuffles DORADO class labels within chromosome",
                ],
            },
        },
        {
            "id": "gene_drug_profile",
            "label": "GENCODE/COSMIC/DGIdb region-annotation profile",
            "path": args.gene_profile.as_posix(),
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Deduplicated region-to-GENCODE/COSMIC/DGIdb/CLP annotation plus fixed-seed chromosome and length-decile conditional null testing.",
                "executed_at": json.loads(args.gene_profile.read_text(encoding="utf-8"))["generated_at"],
                "tables_used": [args.annotation_reproducibility.as_posix(), args.annotation_reproducibility_json.as_posix(), args.annotation_agreement.as_posix(), args.notable_loci.as_posix(), args.gene_summary.as_posix(), args.source_inventory.as_posix()],
                "filters": [
                    "primary analysis uses gene-body overlap",
                    "same-coordinate complete-both HCC1395 pairs for annotation-stratified agreement",
                    "5,000 fixed-seed conditional hypergeometric draws preserving chromosome and region-length decile counts",
                    "DGIdb approved antineoplastic is a gene-drug record, not mutation-specific actionability",
                ],
                "metric_definitions": [
                    "annotation agreement = matching coarse classes / exact complete-both pairs within the annotation stratum",
                    "present-minus-absent difference is compared with a chromosome and length-decile conditional null",
                ],
            },
        },
        {
            "id": "vaf_method",
            "label": "Family-specific per-site read-AF exact-tree ranking",
            "path": args.vaf_summary.as_posix(),
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Ranks all non-capped exact-tree candidates within one HP using per-site read-AF ancestry ordering and a softmax relative weight.",
                "executed_at": generated_at,
                "tables_used": [args.vaf_summary.as_posix(), args.vaf_example.as_posix()],
                "filters": ["mutation-bearing HP1/HP2 units", "all candidate sets re-enumerated with tree_cap=0"],
                "metric_definitions": [
                    "read_AF = ALT reads / (REF reads + ALT reads) within the same HP",
                    "ancestor mutations are heuristically expected to have read_AF greater than or equal to descendant mutations",
                    "softmax relative weight is not a calibrated posterior",
                ],
            },
        },
        {
            "id": "vaf_pair_analysis",
            "label": "HCC1395 pair VAF-selected forest comparison",
            "path": args.vaf_pair_summary.as_posix(),
            "query": {
                "engine": "python",
                "language": "python",
                "description": "Compares VAF exact-score maximizer forests and rooted-unlabeled selected shapes between HCC1395 and HCC1395_DORADO; exact-score ties are retained and HP1/HP2 swap is reported only as a sensitivity view.",
                "executed_at": generated_at,
                "tables_used": [
                    args.vaf_pair_metrics.as_posix(),
                    args.vaf_pair_summary.as_posix(),
                    args.vaf_pair_regions.as_posix(),
                    args.vaf_pair_checks.as_posix(),
                ],
                "filters": [
                    "chr1-22 exact-coordinate complete-both pairs",
                    "same-HP per-site read-AF ordering score",
                    "exact Fraction-score ties retained",
                    "original Topo>1 shape selection fails closed unless every ambiguous unit is VAF-evaluable",
                ],
                "metric_definitions": [
                    "selected-shape agreement compares HP-labeled rooted-unlabeled branching skeletons and cannot identify mutation direction",
                    "selected-exact agreement compares HP-specific genomic-position tuples plus exact-score-maximizing candidate IDs",
                    "phase-swap tolerant agreement permits a global HP1/HP2 label exchange but does not change genomic coordinates",
                    "all selected results are heuristic inferences rather than calibrated probabilities or biological truth",
                ],
            },
        },
        {
            "id": "latest_pipeline_status",
            "label": "Point-in-time clean layered-v3 pipeline status audit",
            "path": args.latest_status.as_posix(),
            "query": {
                "engine": "filesystem-audit",
                "language": "python",
                "description": "Freeze-cutoff audit of producer receipts and the absence/presence of canonical and sensitivity clean-v3 outputs.",
                "executed_at": json.loads(args.latest_status.read_text(encoding="utf-8"))["as_of"],
                "tables_used": [args.latest_status.as_posix()],
            },
        },
    ]

    frozen: dict[str, list[dict[str, Any]]] = {}
    for dataset_name, rows in datasets.items():
        query, frozen_rows = freeze_dataset_with_sql(dataset_name, rows)
        frozen[dataset_name] = frozen_rows
        description, upstream_paths = dataset_upstream[dataset_name]
        sources.append(
            {
                "id": f"dataset_{dataset_name}",
                "label": f"Frozen SQL snapshot: {dataset_name}",
                "path": upstream_paths[0],
                "query": {
                    "engine": "sqlite",
                    "language": "sql",
                    "description": description,
                    "executed_at": generated_at,
                    "sql": query,
                    "tables_used": upstream_paths,
                },
            }
        )
    return sources, frozen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", type=Path, default=TOPIC / "data/topology_pair_analysis.json")
    parser.add_argument("--summary", type=Path, default=TOPIC / "data/coarse_topology_all_dataset_summary.tsv")
    parser.add_argument("--confusion", type=Path, default=TOPIC / "data/hcc1395_pair_confusion.tsv")
    parser.add_argument("--matches", type=Path, default=TOPIC / "data/hcc1395_pair_matches.tsv")
    parser.add_argument("--chrom", type=Path, default=TOPIC / "data/hcc1395_pair_interval_by_chrom.tsv")
    parser.add_argument("--match-metrics", type=Path, default=TOPIC / "data/hcc1395_pair_match_metrics.tsv")
    parser.add_argument("--gene-profile", type=Path, default=TOPIC / "agents/gene_drug_source_profile.json")
    parser.add_argument("--annotation-agreement", type=Path, default=TOPIC / "agents/hcc1395_exact_complete_annotation_agreement.tsv")
    parser.add_argument("--annotation-reproducibility", type=Path, default=TOPIC / "data/hcc1395_annotation_reproducibility.tsv")
    parser.add_argument("--annotation-reproducibility-json", type=Path, default=TOPIC / "data/hcc1395_annotation_reproducibility.json")
    parser.add_argument("--source-inventory", type=Path, default=TOPIC / "agents/source_inventory.tsv")
    parser.add_argument("--notable-loci", type=Path, default=TOPIC / "agents/hcc1395_notable_gene_locus_sensitivity.tsv")
    parser.add_argument("--gene-summary", type=Path, default=TOPIC / "agents/hcc1395_exact_complete_gene_summary.tsv")
    parser.add_argument("--latest-status", type=Path, default=TOPIC / "data/latest_pipeline_status.json")
    parser.add_argument("--vaf-summary", type=Path, default=Path("research/20260711_read_group_C_tree_T_topology_report/data/read_af_tree_ordering_historical.json"))
    parser.add_argument("--vaf-example", type=Path, default=Path("research/20260711_read_group_C_tree_T_topology_report/data/HCC1395_VAF_multi_T_example.json"))
    parser.add_argument("--vaf-pair-metrics", type=Path, default=TOPIC / "data/hcc1395_vaf_pair_metrics.tsv")
    parser.add_argument("--vaf-pair-summary", type=Path, default=TOPIC / "data/hcc1395_vaf_pair_summary.json")
    parser.add_argument("--vaf-pair-regions", type=Path, default=TOPIC / "data/hcc1395_vaf_pair_regions.tsv")
    parser.add_argument("--vaf-pair-checks", type=Path, default=TOPIC / "data/hcc1395_vaf_pair_checks.tsv")
    parser.add_argument(
        "--cosmic-sample-specific",
        type=Path,
        action="append",
        default=[],
        help="Optional bounded TSV/JSON COSMIC cell-line summary; may be repeated.",
    )
    parser.add_argument("--report-md", type=Path, default=REPORT_DIR / f"{REPORT_STEM}.md")
    parser.add_argument("--artifact", type=Path, default=REPORT_DIR / "artifact.json")
    args = parser.parse_args()

    required = [
        args.analysis,
        args.summary,
        args.confusion,
        args.matches,
        args.chrom,
        args.match_metrics,
        args.gene_profile,
        args.annotation_agreement,
        args.annotation_reproducibility,
        args.annotation_reproducibility_json,
        args.source_inventory,
        args.notable_loci,
        args.gene_summary,
        args.latest_status,
        args.vaf_summary,
        args.vaf_example,
        args.vaf_pair_metrics,
        args.vaf_pair_summary,
        args.vaf_pair_regions,
        args.vaf_pair_checks,
    ]
    missing = [path.as_posix() for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing required report inputs: " + ", ".join(missing))

    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    analysis = json.loads(args.analysis.read_text(encoding="utf-8"))
    gene_profile = json.loads(args.gene_profile.read_text(encoding="utf-8"))
    summary_raw = read_tsv(args.summary)
    confusion_raw = read_tsv(args.confusion)
    matches_raw = read_tsv(args.matches)
    chrom_raw = read_tsv(args.chrom)
    match_metrics_raw = read_tsv(args.match_metrics)
    annotation_raw = read_tsv(args.annotation_agreement)
    annotation_repro_raw = read_tsv(args.annotation_reproducibility)
    annotation_repro_meta = json.loads(args.annotation_reproducibility_json.read_text(encoding="utf-8"))
    inventory_raw = read_tsv(args.source_inventory)
    notable_loci_raw = read_tsv(args.notable_loci)
    gene_summary_raw = read_tsv(args.gene_summary)
    latest_status = json.loads(args.latest_status.read_text(encoding="utf-8"))
    vaf_summary = json.loads(args.vaf_summary.read_text(encoding="utf-8"))
    vaf_example = json.loads(args.vaf_example.read_text(encoding="utf-8"))
    vaf_pair_metrics_raw = read_tsv(args.vaf_pair_metrics)
    vaf_pair_summary_raw = json.loads(args.vaf_pair_summary.read_text(encoding="utf-8"))
    vaf_pair_checks_raw = read_tsv(args.vaf_pair_checks)

    summary_by_sample = {row["sample"]: row for row in summary_raw}
    if set(summary_by_sample) != set(SAMPLE_ORDER):
        raise RuntimeError(f"dataset rows differ from locked seven-row scope: {sorted(summary_by_sample)}")
    if analysis.get("status") != "PASS" or analysis.get("validation", {}).get("checks") != analysis.get("validation", {}).get("passed"):
        raise RuntimeError("topology pair analysis is not internally PASS")

    composition_long: list[dict[str, Any]] = []
    composition_table: list[dict[str, Any]] = []
    for sample in SAMPLE_ORDER:
        raw = summary_by_sample[sample]
        complete = as_int(raw["complete_regions"])
        table_row: dict[str, Any] = {
            "sample": sample,
            "primary_regions": as_int(raw["primary_regions"]),
            "complete_regions": complete,
            "incomplete_regions": as_int(raw["incomplete_regions"]),
        }
        for category in CATEGORY_ORDER:
            count = as_int(raw[category])
            share = as_float(raw[f"share_{category}"])
            composition_long.append(
                {
                    "sample": sample,
                    "category": CATEGORY_ZH[category],
                    "category_chart": CATEGORY_CHART_ZH[category],
                    "category_key": category,
                    "count": count,
                    "share": share,
                    "complete_regions": complete,
                }
            )
            table_row[f"{category}_n"] = count
            table_row[f"{category}_share"] = share
        composition_table.append(table_row)

    hcc_row = summary_by_sample["HCC1395"]
    dorado_row = summary_by_sample["HCC1395_DORADO"]
    hcc_share_delta: list[dict[str, Any]] = []
    for category in CATEGORY_ORDER:
        share_hcc = as_float(hcc_row[f"share_{category}"])
        share_dorado = as_float(dorado_row[f"share_{category}"])
        hcc_share_delta.append(
            {
                "category": CATEGORY_ZH[category],
                "category_key": category,
                "HCC1395_share": share_hcc,
                "DORADO_share": share_dorado,
                "DORADO_minus_HCC_pp": (share_dorado - share_hcc) * 100,
                "absolute_gap_pp": abs(share_dorado - share_hcc) * 100,
            }
        )
    aggregate_tvd = sum(row["absolute_gap_pp"] for row in hcc_share_delta) / 200
    aggregate_tvd_pp = aggregate_tvd * 100

    exact_metrics = next(row for row in analysis["hcc1395_pair_metrics"] if row["scenario"] == "exact_coordinate")
    verdict_metrics = [
        {
            "exact_matched_regions": exact_metrics["matched_all"],
            "HCC_interval_coverage": exact_metrics["match_coverage_a"],
            "DORADO_interval_coverage": exact_metrics["match_coverage_b"],
            "complete_both": exact_metrics["complete_both"],
            "coarse_agreement": exact_metrics["raw_agreement"],
            "coarse_agreement_ci_low": exact_metrics["raw_agreement_ci95_low"],
            "coarse_agreement_ci_high": exact_metrics["raw_agreement_ci95_high"],
            "cohen_kappa": exact_metrics["cohen_kappa"],
            "null_agreement_mean": exact_metrics["permutation_null"]["agreement_mean"],
            "null_agreement_q975": exact_metrics["permutation_null"]["agreement_q975"],
            "null_kappa_q975": exact_metrics["permutation_null"]["kappa_q975"],
            "permutation_p_ge": exact_metrics["permutation_null"]["agreement_p_ge"],
            "ordered_exact_tree_agreement": exact_metrics["ordered_exact_tree_set_digest_agreement"],
            "phase_swap_tree_agreement": exact_metrics["phase_swap_tolerant_exact_tree_set_digest_agreement"],
            "aggregate_profile_TVD": aggregate_tvd,
        }
    ]

    pair_sensitivity: list[dict[str, Any]] = []
    for row in match_metrics_raw:
        pair_sensitivity.append(
            {
                "scenario": SCENARIO_ZH.get(row["scenario"], row["scenario"]),
                "scenario_key": row["scenario"],
                "matched_all": as_int(row["matched_all"]),
                "complete_both": as_int(row["complete_both"]),
                "coarse_agreement": as_float(row["raw_agreement"]),
                "cohen_kappa": as_float(row["cohen_kappa"]),
                "macro_Jaccard": as_float(row["category_jaccard_macro"]),
                "phase_swap_HP_signature": as_float(row["phase_swap_tolerant_hp_coarse_signature_agreement"]),
                "ordered_exact_tree": as_float(row["ordered_exact_tree_set_digest_agreement"]),
                "phase_swap_exact_tree": as_float(row["phase_swap_tolerant_exact_tree_set_digest_agreement"]),
                "null_agreement_q975": as_float(row["null_agreement_q975"]),
                "null_kappa_q975": as_float(row["null_kappa_q975"]),
                "permutation_p_ge": as_float(row["null_agreement_p_ge"]),
            }
        )

    confusion_exact = [row for row in confusion_raw if row["scenario"] == "exact_coordinate"]
    confusion_lookup = {(row["category_a"], row["category_b"]): as_int(row["n"]) for row in confusion_exact}
    confusion_heatmap: list[dict[str, Any]] = []
    confusion_long: list[dict[str, Any]] = []
    for category_a in CATEGORY_ORDER:
        total = sum(confusion_lookup[(category_a, category_b)] for category_b in CATEGORY_ORDER)
        wide: dict[str, Any] = {
            "HCC1395_class": CATEGORY_ZH[category_a],
            "HCC1395_class_key": category_a,
            "row_total": total,
        }
        for category_b in CATEGORY_ORDER:
            n = confusion_lookup[(category_a, category_b)]
            wide[f"to_{category_b}"] = pct(n, total)
            confusion_long.append(
                {
                    "HCC1395_class": CATEGORY_ZH[category_a],
                    "DORADO_class": CATEGORY_ZH[category_b],
                    "n": n,
                    "row_total": total,
                    "row_share": pct(n, total),
                    "diagonal": category_a == category_b,
                }
            )
        confusion_heatmap.append(wide)

    chrom_exact: list[dict[str, Any]] = []
    for row in sorted((row for row in chrom_raw if row["scenario"] == "exact_coordinate"), key=lambda item: chr_sort_key(item["chrom"])):
        chrom_exact.append(
            {
                "chrom": row["chrom"],
                "matched_all": as_int(row["matched_all"]),
                "complete_both": as_int(row["complete_both"]),
                "agreements": as_int(row["agreements"]),
                "coarse_agreement": as_float(row["raw_agreement"]),
            }
        )

    bin_size = 25_000_000
    bin_counts: dict[tuple[str, int], list[int]] = defaultdict(lambda: [0, 0])
    for row in matches_raw:
        if row["scenario"] != "exact_coordinate" or not (as_bool(row["complete_a"]) and as_bool(row["complete_b"])):
            continue
        midpoint = parse_region_midpoint(row["region_a"])
        key = (row["chrom"], midpoint // bin_size)
        bin_counts[key][0] += 1
        bin_counts[key][1] += int(as_bool(row["category_agree"]))
    max_bin = max(bin_index for _, bin_index in bin_counts)
    genomic_bin_fields = [f"bin_{index:02d}" for index in range(max_bin + 1)]
    genomic_bin_labels = {
        f"bin_{index:02d}": f"{index * 25}–{(index + 1) * 25} Mb" for index in range(max_bin + 1)
    }
    genomic_heatmap: list[dict[str, Any]] = []
    genomic_bins_long: list[dict[str, Any]] = []
    for chrom in sorted({chrom for chrom, _ in bin_counts}, key=chr_sort_key):
        wide: dict[str, Any] = {"chrom": chrom}
        for index, field in enumerate(genomic_bin_fields):
            n, agree = bin_counts.get((chrom, index), [0, 0])
            wide[field] = pct(agree, n)
            wide[f"{field}_n"] = n
            if n:
                genomic_bins_long.append(
                    {
                        "chrom": chrom,
                        "bin": genomic_bin_labels[field],
                        "bin_start": index * bin_size,
                        "bin_end": (index + 1) * bin_size - 1,
                        "complete_both": n,
                        "agreements": agree,
                        "coarse_agreement": pct(agree, n),
                    }
                )
        genomic_heatmap.append(wide)

    annotation_long: list[dict[str, Any]] = []
    annotation_compare: list[dict[str, Any]] = []
    repro_baseline = next(row for row in annotation_repro_raw if row["feature"] == "ALL")
    if as_int(repro_baseline["n_present"]) != exact_metrics["complete_both"]:
        raise RuntimeError("annotation reproducibility baseline denominator differs from exact complete-both")
    if not math.isclose(as_float(repro_baseline["present_agreement"]), exact_metrics["raw_agreement"], abs_tol=1e-12):
        raise RuntimeError("annotation reproducibility baseline agreement differs from pair analysis")
    for row in annotation_repro_raw:
        feature = row["feature"]
        if feature == "ALL":
            continue
        label = REPRO_FEATURE_ZH.get(feature, row["feature_label"])
        present_n = as_int(row["n_present"])
        absent_n = as_int(row["n_absent"])
        present_agreement = as_float(row["present_agreement"])
        absent_agreement = as_float(row["absent_agreement"])
        for stratum, n, agreement, kappa in (
            ("有註解", present_n, present_agreement, as_float(row["present_kappa"])),
            ("無註解", absent_n, absent_agreement, as_float(row["absent_kappa"])),
        ):
            annotation_long.append(
                {
                    "feature": label,
                    "feature_key": feature,
                    "stratum": stratum,
                    "n_exact_complete_pairs": n,
                    "coarse_agreement": agreement,
                    "cohen_kappa": kappa,
                }
            )
        conditional_p = as_float(row["permutation_p_two_sided"])
        annotation_compare.append(
            {
                "feature": label,
                "feature_key": feature,
                "present_n": present_n,
                "present_agreement": present_agreement,
                "present_ci95_low": as_float(row["present_ci95_low"]),
                "present_ci95_high": as_float(row["present_ci95_high"]),
                "present_kappa": as_float(row["present_kappa"]),
                "absent_n": absent_n,
                "absent_agreement": absent_agreement,
                "absent_kappa": as_float(row["absent_kappa"]),
                "present_minus_absent_pp": as_float(row["difference_present_minus_absent_pp"]),
                "conditional_null_mean_pp": as_float(row["permutation_null_mean_pp"]),
                "conditional_null_q025_pp": as_float(row["permutation_null_q025_pp"]),
                "conditional_null_q975_pp": as_float(row["permutation_null_q975_pp"]),
                "conditional_p_two_sided": conditional_p,
                "significant_0_05": conditional_p < 0.05,
                "interpretation": "未顯著；僅 context/face validity" if conditional_p >= 0.05 else "conditional-null significant",
            }
        )

    gene_join_summary: list[dict[str, Any]] = []
    for sample in ("HCC1395", "HCC1395_DORADO"):
        row = gene_profile["join"][sample]
        gene_join_summary.append(
            {
                "sample": sample,
                "primary_regions": row["regions_before_join"],
                "gene_body_regions": row["distinct_regions_with_body_gene"],
                "gene_body_coverage": row["body_gene_region_coverage_pct"] / 100,
                "CGC_body_regions": row["distinct_regions_with_cgc_body"],
                "CGC_body_coverage": row["cgc_body_region_coverage_pct"] / 100,
                "DGIdb_interaction_body_regions": row["distinct_regions_with_dgidb_interaction_body"],
                "DGIdb_interaction_body_coverage": row["dgidb_interaction_body_region_coverage_pct"] / 100,
                "DGIdb_approved_antineoplastic_body_regions": row["distinct_regions_with_dgidb_approved_antineoplastic_body"],
                "DGIdb_approved_antineoplastic_body_coverage": row["dgidb_approved_antineoplastic_body_region_coverage_pct"] / 100,
            }
        )

    source_inventory: list[dict[str, Any]] = []
    for row in inventory_raw:
        source_inventory.append(
            {
                "source_id": row["source_id"],
                "source_file": Path(row["absolute_path"]).name,
                "absolute_path": row["absolute_path"],
                "version_or_snapshot": row["version_or_snapshot"],
                "genome_build": row["genome_build"],
                "grain": row["grain"],
                "row_count": as_int(row["row_count"]),
                "join_key": row["join_key"],
                "coverage_caveat": row["coverage_caveat"],
                "decision": row["decision"],
                "sha256": row["sha256"],
            }
        )

    vaf_sites = [
        {
            "site": f"site_{as_int(row['index']) + 1}",
            "position": as_int(row["position"]),
            "REF_reads": as_int(row["REF_reads"]),
            "ALT_reads": as_int(row["ALT_reads"]),
            "read_AF": as_float(row["read_AF"]),
        }
        for row in vaf_example["sites"]
    ]
    vaf_candidates = [
        {
            "candidate": row["candidate"],
            "edges": "; ".join(f"{parent}→{child}" for parent, child in row["edges"]),
            "shape_signature": row["shape_signature"],
            "score": as_float(row["score"]),
            "softmax_relative_weight": as_float(row["softmax_weight"]),
            "selected": bool(row["selected"]),
        }
        for row in vaf_example["candidates"]
    ]

    if vaf_pair_summary_raw.get("validation", {}).get("passed") != vaf_pair_summary_raw.get("validation", {}).get("total"):
        raise RuntimeError("VAF pair analysis did not pass all validation checks")
    if len(vaf_pair_checks_raw) != 20 or not all(as_bool(row["pass"]) for row in vaf_pair_checks_raw):
        raise RuntimeError("VAF pair checks TSV is not 20/20 PASS")

    vaf_metric_lookup = {
        (row["metric_family"], row["population"], row["comparison"]): row
        for row in vaf_pair_metrics_raw
    }

    def find_vaf_metric(metric_family: str, population: str, comparison: str) -> dict[str, str]:
        key = (metric_family, population, comparison)
        if key not in vaf_metric_lookup:
            raise RuntimeError(f"missing VAF pair metric: {key}")
        return vaf_metric_lookup[key]

    def vaf_comparison_row(
        layer: str,
        metric_family: str,
        population: str,
        selection_basis: str,
        claim_ceiling_text: str,
    ) -> dict[str, Any]:
        ordered = find_vaf_metric(metric_family, population, "ordered_HP")
        swapped = find_vaf_metric(metric_family, population, "HP1_HP2_swap_tolerant")
        denominator = as_int(ordered["denominator"])
        if denominator != as_int(swapped["denominator"]):
            raise RuntimeError(f"ordered/swap VAF denominator mismatch: {metric_family} {population}")
        coverage = pct(denominator, exact_metrics["complete_both"])
        ordered_agreement = pct(as_int(ordered["numerator"]), denominator)
        phase_swap_agreement = pct(as_int(swapped["numerator"]), denominator)
        return {
            "layer": layer,
            "population": population,
            "denominator": denominator,
            "coverage_of_5720": coverage,
            "coverage_display": fmt_pct(coverage),
            "ordered_n": as_int(ordered["numerator"]),
            "ordered_agreement": ordered_agreement,
            "ordered_display": fmt_pct(ordered_agreement),
            "phase_swap_n": as_int(swapped["numerator"]),
            "phase_swap_agreement": phase_swap_agreement,
            "phase_swap_display": fmt_pct(phase_swap_agreement),
            "selection_basis": selection_basis,
            "claim_ceiling": claim_ceiling_text,
        }

    vaf_pair_comparison = [
        vaf_comparison_row(
            "未排名 exact 候選樹集合（基準）",
            "unranked_full_candidate_tree_set_digest",
            "all_exact_coordinate_complete_both",
            "完整可行候選集合；尚未使用 VAF 選第一順位",
            "候選空間技術再現性；不是 selected tree",
        ),
        vaf_comparison_row(
            "VAF 第一順位 exact 候選集合（保留並列）",
            "vaf_exact_first_candidate_set",
            "both_sides_actually_use_vaf_ranking",
            "兩側都以 same-HP read-AF exact score 排序；最高分並列全保留",
            "VAF-supported 推測候選集合；不是 posterior 或 truth",
        ),
        vaf_comparison_row(
            "VAF 唯一 exact forest",
            "vaf_unique_exact_first_candidate",
            "both_sides_actually_use_vaf_ranking_and_unique",
            "兩側都實際使用 VAF，且各自只剩一個 HP-specific genomic-position＋candidate-ID forest",
            "可比較 mutation-labeled 推測樹；仍不是 biological tree accuracy",
        ),
        vaf_comparison_row(
            "結構-first＋VAF 單一 shape",
            "selected_single_shape",
            "both_structural_or_vaf_shape_selectable",
            "原 Topo=1 使用結構 shape；原 Topo>1 必須所有 ambiguous units 可評估且 VAF top set 落同一 shape",
            "只比較去 mutation labels 的 branching skeleton",
        ),
        vaf_comparison_row(
            "兩側皆用 VAF 且為單一 shape",
            "selected_single_shape",
            "both_sides_vaf_evaluable_and_single_shape",
            "兩側都使用 VAF exact ranking，且完整 forest 各只對應一種 rooted-unlabeled shape",
            "shape 相同不代表 mutation direction 相同",
        ),
        vaf_comparison_row(
            "兩側原 Topo>1，VAF 各縮至單一 shape",
            "selected_single_shape",
            "both_shapes_require_vaf_from_original_Topo_gt1",
            "最乾淨的 VAF shape-rescue subset；兩側原始候選皆跨 shape",
            "VAF-supported 推測 shape；77.74% 不是 exact-tree accuracy",
        ),
    ]
    vaf_exact_unique = next(row for row in vaf_pair_comparison if row["population"] == "both_sides_actually_use_vaf_ranking_and_unique")
    vaf_shape_rescue = next(row for row in vaf_pair_comparison if row["population"] == "both_shapes_require_vaf_from_original_Topo_gt1")
    vaf_shape_all = next(row for row in vaf_pair_comparison if row["population"] == "both_structural_or_vaf_shape_selectable")
    vaf_conservation = vaf_pair_summary_raw["sample_conservation"]

    pipeline_status = [
        {
            "as_of": latest_status["as_of"],
            "producer_terminal_pass": latest_status["producer"]["terminal_pass"],
            "producer_expected": len(latest_status["producer"]["datasets"]),
            "producer_active": latest_status["producer"]["active"],
            "producer_failed": latest_status["producer"]["failed"],
            "aggregate_success_present": bool(latest_status["producer"]["aggregate_success_present"]),
            "producer_verification_present": bool(latest_status["producer"]["verification_summary_present"]),
            "receipt_closeout_status": latest_status["closeout"]["status"],
            "receipt_closeout_error": latest_status["closeout"]["error_code"],
            "continuation_status": latest_status["continuation_execution"]["status"],
            "canonical_root_present": bool(latest_status["canonical"]["root_present"]),
            "canonical_success_present": bool(latest_status["canonical"]["success_present"]),
            "sensitivity_root_present": bool(latest_status["sensitivity"]["root_present"]),
            "sensitivity_success_present": bool(latest_status["sensitivity"]["success_present"]),
            "latest_complete_source": latest_status["latest_complete_quantitative_source"]["classification"],
            "scientific_release_allowed": bool(latest_status["latest_complete_quantitative_source"]["scientific_release_allowed"]),
            "verdict": latest_status["verdict"],
        }
    ]

    notable_loci: list[dict[str, Any]] = []
    for row in notable_loci_raw:
        scenario = row["best_pair_scenario"] or "unmatched"
        category_a = CATEGORY_ZH.get(row["best_pair_category_a"], row["best_pair_category_a"] or "—")
        category_b = CATEGORY_ZH.get(row["best_pair_category_b"], row["best_pair_category_b"] or "—")
        reciprocal = as_float(row["best_pair_reciprocal_overlap"]) if row["best_pair_reciprocal_overlap"] else None
        notable_loci.append(
            {
                "gene": row["gene_symbol"],
                "HCC1395_region": row["hcc1395_body_regions"] or "無對應 region",
                "HCC1395_status": f"{row['hcc1395_structural_classes'] or 'unmatched'} / {CATEGORY_ZH.get(row['hcc1395_coarse_categories'], row['hcc1395_coarse_categories'] or '—')}",
                "DORADO_region": row["dorado_body_regions"] or "無對應 region",
                "DORADO_status": f"{row['dorado_structural_classes'] or 'unmatched'} / {CATEGORY_ZH.get(row['dorado_coarse_categories'], row['dorado_coarse_categories'] or '—')}",
                "best_match": SCENARIO_ZH.get(scenario, scenario),
                "reciprocal_overlap": reciprocal,
                "complete_both": "yes" if row["best_pair_complete_a"] == "True" and row["best_pair_complete_b"] == "True" else "no",
                "coarse_pair": f"{category_a} ↔ {category_b}" if row["best_pair_category_a"] else "—",
                "coarse_agree": "yes" if row["best_pair_category_agree"] == "True" else "no",
                "interpretation": row["interpretation"],
            }
        )

    eligible_gene_rows = [
        row
        for row in gene_summary_raw
        if as_int(row["n_exact_complete_regions_with_body_overlap"]) > 0
        and (as_bool(row["cgc_v104"]) or as_int(row["dgidb_approved_antineoplastic_unique_drugs"]) > 0)
    ]
    eligible_gene_rows.sort(
        key=lambda row: (
            not as_bool(row["user_notable_gene"]),
            not as_bool(row["cgc_v104"]),
            as_int(row["cgc_tier"]) if row["cgc_tier"] else 99,
            -as_int(row["n_exact_complete_regions_with_body_overlap"]),
            -as_int(row["dgidb_approved_antineoplastic_unique_drugs"]),
            row["gene_symbol"],
        )
    )
    gene_locus_examples = [
        {
            "gene": row["gene_symbol"],
            "CGC": bool(as_bool(row["cgc_v104"])),
            "CGC_role": row["cgc_role"] or "—",
            "CGC_tier": as_int(row["cgc_tier"]) if row["cgc_tier"] else None,
            "exact_complete_regions": as_int(row["n_exact_complete_regions_with_body_overlap"]),
            "category_agreement": as_float(row["category_agreement_pct"]) / 100,
            "cohen_kappa": as_float(row["cohens_kappa"]) if row["cohens_kappa"] else None,
            "approved_antineoplastic_drugs_n": as_int(row["dgidb_approved_antineoplastic_unique_drugs"]),
            "approved_antineoplastic_claim_names": summarize_pipe(row["approved_antineoplastic_drug_claim_names"], 8) or "—",
            "region_examples": summarize_pipe(row["region_keys"], 4) or "—",
            "claim_ceiling": "Non-actionable descriptive gene-drug records; no mutation-specific treatment inference",
        }
        for row in eligible_gene_rows[:30]
    ]

    definitions = [
        {"class": CATEGORY_ZH["no_within_hp_relation"], "rule": "Topo=1 且 primary HP 圖中沒有分岔、也沒有深度≥2 的直系鏈", "boundary": "可能是單 HP 單節點，也可能是 HP1/HP2 各自孤立；不等於單一 clone"},
        {"class": CATEGORY_ZH["sister_only"], "rule": "Topo=1 且至少一個 primary HP 節點 out-degree≥2，但 max depth<2", "boundary": "分岔只在同一 HP 內計算"},
        {"class": CATEGORY_ZH["direct_only"], "rule": "Topo=1 且至少一個 primary HP max depth≥2，但無分岔", "boundary": "H_ intermediate/partial state 參與圖深度"},
        {"class": CATEGORY_ZH["sister_and_direct"], "rule": "Topo=1 且同一 region 的 primary HP forest 同時出現分岔與深度≥2", "boundary": "雙 HP 保留 ordered forest，不跨 HP 連邊"},
        {"class": CATEGORY_ZH["topology_multiple_unresolved"], "rule": "Topo>1；候選集仍包含一種以上形狀", "boundary": "強制保留未定，不用某一候選形狀代表 region"},
    ]

    cgc = next(row for row in annotation_compare if row["feature_key"] == "cgc_body_any")
    anti = next(row for row in annotation_compare if row["feature_key"] == "dgidb_approved_antineoplastic_body_any")
    clp_all = next(row for row in annotation_compare if row["feature_key"] == "clp_all_variant_any")
    clp_confirmed = next(row for row in annotation_compare if row["feature_key"] == "clp_confirmed_somatic_variant_any")

    claim_ceiling = [
        {"layer": "區間再現性", "verdict": "PARTIAL SUPPORT", "evidence": f"exact coordinate={exact_metrics['matched_all']:,}；coverage={fmt_pct(exact_metrics['match_coverage_a'])}/{fmt_pct(exact_metrics['match_coverage_b'])}", "allowed_claim": "同一 HCC1395 生物樣本的兩個技術 dataset 有高度區間重疊"},
        {"layer": "粗拓撲類別", "verdict": "PARTIAL SUPPORT", "evidence": f"agreement={fmt_pct(exact_metrics['raw_agreement'])}；κ={exact_metrics['cohen_kappa']:.3f}；null 97.5%={fmt_pct(exact_metrics['permutation_null']['agreement_q975'])}", "allowed_claim": "粗類別一致性高於染色體保留亂數基準，但不是等價"},
        {"layer": "Exact T 集合", "verdict": "WEAK/PARTIAL", "evidence": f"ordered={fmt_pct(exact_metrics['ordered_exact_tree_set_digest_agreement'])}；phase-swap tolerant={fmt_pct(exact_metrics['phase_swap_tolerant_exact_tree_set_digest_agreement'])}", "allowed_claim": "越細粒度的 tree identity 技術再現性越低"},
        {"layer": "VAF-supported selected T／shape", "verdict": "HEURISTIC ONLY", "evidence": f"原 Topo>1 雙側單-shape subset：ordered={fmt_pct(vaf_shape_rescue['ordered_agreement'])}、swap={fmt_pct(vaf_shape_rescue['phase_swap_agreement'])} (n={vaf_shape_rescue['denominator']:,})；雙側 VAF 唯一 exact forest：ordered={fmt_pct(vaf_exact_unique['ordered_agreement'])}、swap={fmt_pct(vaf_exact_unique['phase_swap_agreement'])} (n={vaf_exact_unique['denominator']:,})", "allowed_claim": "使用 same-HP read-AF exact-score 第一順位作推測比較；shape 不含 mutation direction，exact forest 也不是 calibrated posterior 或 biological truth"},
        {"layer": "癌症基因/藥物關聯", "verdict": "DESCRIPTIVE ONLY / NULL", "evidence": f"CGC delta={cgc['present_minus_absent_pp']:.2f} pp, conditional p={cgc['conditional_p_two_sided']:.4f}；approved-antineoplastic delta={anti['present_minus_absent_pp']:.2f} pp, p={anti['conditional_p_two_sided']:.4f}", "allowed_claim": "可做 context/face validity 與後續優先排序；conditional null 不支持 enrichment，也不是用藥建議"},
        {"layer": "COSMIC CLP HCC1395", "verdict": "DESCRIPTIVE ONLY / NULL", "evidence": f"all CLP n={clp_all['present_n']}, agreement={fmt_pct(clp_all['present_agreement'])}, p={clp_all['conditional_p_two_sided']:.4f}；confirmed somatic n={clp_confirmed['present_n']}, agreement={fmt_pct(clp_confirmed['present_agreement'])}, 95% CI={fmt_pct(clp_confirmed['present_ci95_low'])}–{fmt_pct(clp_confirmed['present_ci95_high'])}, p={clp_confirmed['conditional_p_two_sided']:.4f}", "allowed_claim": "sample-specific allele overlap 未顯著，只能作輔助脈絡"},
        {"layer": "方法合理有效/真 clone", "verdict": "SCIENTIFIC NO-GO", "evidence": "historical layered-v2；無 single-cell/multi-region truth；clean-v3 未 closeout", "allowed_claim": "僅支持 engineering coherence 與部分 technical reproducibility"},
    ]

    validation_rows = [
        {"check": "Five-class conservation", "observed": f"{analysis['validation']['passed']}/{analysis['validation']['checks']}", "status": "PASS"},
        {"check": "Exact complete-both denominator", "observed": exact_metrics["complete_both"], "status": "PASS" if sum(as_int(row["n"]) for row in confusion_exact) == exact_metrics["complete_both"] else "FAIL"},
        {"check": "Annotation strata conservation", "observed": gene_profile["topology_pair"]["all_strata_conserve_5720"], "status": "PASS" if gene_profile["topology_pair"]["all_strata_conserve_5720"] else "FAIL"},
        {"check": "Annotation conditional-null permutations", "observed": annotation_repro_meta["permutation"]["n_permutations"], "status": "PASS" if annotation_repro_meta["permutation"]["n_permutations"] >= 5000 else "FAIL"},
        {"check": "Annotation strata significant at 0.05", "observed": sum(row["significant_0_05"] for row in annotation_compare), "status": "PASS" if not any(row["significant_0_05"] for row in annotation_compare) else "FAIL"},
        {"check": "Permutation count", "observed": exact_metrics["permutation_null"]["permutations"], "status": "PASS" if exact_metrics["permutation_null"]["permutations"] >= 2000 else "FAIL"},
        {"check": "VAF-selected pair analysis", "observed": f"{vaf_pair_summary_raw['validation']['passed']}/{vaf_pair_summary_raw['validation']['total']}", "status": "PASS" if vaf_pair_summary_raw["validation"]["passed"] == vaf_pair_summary_raw["validation"]["total"] == 20 else "FAIL"},
        {"check": "VAF shape conservation", "observed": f"HCC={vaf_conservation['HCC1395']['structural_or_vaf_shape_selectable']}; DORADO={vaf_conservation['HCC1395_DORADO']['structural_or_vaf_shape_selectable']}", "status": "PASS" if vaf_conservation["HCC1395"]["structural_or_vaf_shape_selectable"] == 6798 and vaf_conservation["HCC1395_DORADO"]["structural_or_vaf_shape_selectable"] == 6082 else "FAIL"},
        {"check": "Snapshot release gate", "observed": "historical layered-v2; clean-v3 pending", "status": "PARTIAL"},
    ]
    if any(row["status"] == "FAIL" for row in validation_rows):
        raise RuntimeError(f"report-side conservation failed: {validation_rows}")

    provenance_paths = required + [path for path in args.cosmic_sample_specific if path.exists()]
    provenance = [
        {
            "artifact": path.name,
            "path": safe_display_path(path.resolve().as_posix()),
            "bytes": path.stat().st_size,
            "mtime": datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds"),
            "sha256": sha256(path),
        }
        for path in provenance_paths
    ]

    datasets: dict[str, list[dict[str, Any]]] = {
        "verdict_metrics": verdict_metrics,
        "composition_long": composition_long,
        "composition_table": composition_table,
        "hcc_share_delta": hcc_share_delta,
        "pair_sensitivity": pair_sensitivity,
        "confusion_heatmap": confusion_heatmap,
        "confusion_long": confusion_long,
        "chrom_exact": chrom_exact,
        "genomic_heatmap": genomic_heatmap,
        "genomic_bins_long": genomic_bins_long,
        "annotation_long": annotation_long,
        "annotation_compare": annotation_compare,
        "gene_join_summary": gene_join_summary,
        "source_inventory": source_inventory,
        "vaf_sites": vaf_sites,
        "vaf_candidates": vaf_candidates,
        "vaf_pair_comparison": vaf_pair_comparison,
        "pipeline_status": pipeline_status,
        "notable_loci": notable_loci,
        "gene_locus_examples": gene_locus_examples,
        "definitions": definitions,
        "claim_ceiling": claim_ceiling,
        "validation": validation_rows,
        "provenance": provenance,
    }
    topology_paths = [args.analysis.as_posix(), args.summary.as_posix(), args.confusion.as_posix(), args.matches.as_posix(), args.chrom.as_posix(), args.match_metrics.as_posix()]
    gene_paths = [args.gene_profile.as_posix(), args.annotation_reproducibility.as_posix(), args.annotation_reproducibility_json.as_posix(), args.annotation_agreement.as_posix(), args.notable_loci.as_posix(), args.gene_summary.as_posix(), args.source_inventory.as_posix()]
    dataset_upstream: dict[str, tuple[str, list[str]]] = {
        name: (
            "Literal frozen report snapshot derived from the reviewed HCC1395 coarse-topology pair analysis.",
            topology_paths,
        )
        for name in [
            "verdict_metrics",
            "composition_long",
            "composition_table",
            "hcc_share_delta",
            "pair_sensitivity",
            "confusion_heatmap",
            "confusion_long",
            "chrom_exact",
            "genomic_heatmap",
            "genomic_bins_long",
            "definitions",
            "claim_ceiling",
            "validation",
        ]
    }
    for name in ["annotation_long", "annotation_compare", "gene_join_summary", "source_inventory"]:
        dataset_upstream[name] = (
            "Literal frozen report snapshot derived from the reviewed GENCODE/COSMIC/DGIdb profile and exact-pair annotation strata.",
            gene_paths,
        )
    for name in ["notable_loci", "gene_locus_examples"]:
        dataset_upstream[name] = (
            "Literal frozen gene-locus summary with region examples and non-actionable cancer-gene/drug context.",
            [args.notable_loci.as_posix(), args.gene_summary.as_posix(), args.annotation_reproducibility.as_posix()],
        )
    for name in ["vaf_sites", "vaf_candidates"]:
        dataset_upstream[name] = (
            "Literal frozen teaching example for the same-HP per-site read-AF exact-tree ranking heuristic.",
            [args.vaf_summary.as_posix(), args.vaf_example.as_posix()],
        )
    dataset_upstream["vaf_pair_comparison"] = (
        "Literal frozen comparison of unranked candidate sets, VAF exact-score maximizer forests, and rooted-unlabeled selected shapes for the exact-coordinate complete-both HCC1395 technical pair.",
        [
            args.vaf_pair_metrics.as_posix(),
            args.vaf_pair_summary.as_posix(),
            args.vaf_pair_regions.as_posix(),
            args.vaf_pair_checks.as_posix(),
        ],
    )
    dataset_upstream["pipeline_status"] = (
        "Point-in-time clean-v3 status frozen at the audited cutoff; producer and downstream output gates are reported separately.",
        [args.latest_status.as_posix()],
    )
    dataset_upstream["provenance"] = (
        "File metadata and SHA-256 computed by the report builder over every required input.",
        [path.as_posix() for path in required],
    )

    optional_table_specs: list[dict[str, Any]] = []
    optional_blocks: list[dict[str, Any]] = []
    for index, path in enumerate(args.cosmic_sample_specific, start=1):
        if not path.exists():
            raise FileNotFoundError(f"optional COSMIC input does not exist: {path}")
        rows, full_row_count = load_optional_table(path)
        dataset_id = f"cosmic_sample_specific_{index}"
        datasets[dataset_id] = rows
        dataset_upstream[dataset_id] = (
            f"Optional COSMIC sample-specific summary; bounded to {len(rows)}/{full_row_count} reviewed rows.",
            [safe_display_path(path.resolve().as_posix())],
        )
        optional_table_specs.append(
            {
                "id": f"{dataset_id}_table",
                "title": f"COSMIC sample-specific 輔助資料 {index}",
                "subtitle": f"顯示 {len(rows)}/{full_row_count} rows；僅能作樣本註解，不是 topology truth。",
                "dataset": dataset_id,
                "sourceId": f"dataset_{dataset_id}",
                "density": "dense",
                "defaultSort": {"field": list(rows[0])[0], "direction": "asc"},
                "columns": [
                    {"field": field, "label": field, "type": "text"} for field in rows[0]
                ],
            }
        )
        optional_blocks.extend(
            [
                {"id": f"{dataset_id}_explain", "type": "markdown", "body": f"### COSMIC sample-specific 資料 {index}：僅作輔助註解\n\n該檔只是可選輸入，報告不將 cell-line record 當成 region topology 真值。", "sourceId": f"dataset_{dataset_id}"},
                {"id": f"{dataset_id}_block", "type": "table", "tableId": f"{dataset_id}_table"},
            ]
        )

    sources, frozen_datasets = make_sources(generated_at, args, datasets, dataset_upstream)

    cards = [
        {"id": "interval_card", "description": "來自兩技術 dataset 完全相同座標的 primary regions。", "dataset": "verdict_metrics", "sourceId": "dataset_verdict_metrics", "metrics": [{"label": "Exact matched regions", "field": "exact_matched_regions", "format": "number"}, {"label": "HCC coverage", "field": "HCC_interval_coverage", "format": "percent"}, {"label": "DORADO coverage", "field": "DORADO_interval_coverage", "format": "percent"}]},
        {"id": "agreement_card", "description": "分母是 exact match 中 complete-both regions。", "dataset": "verdict_metrics", "sourceId": "dataset_verdict_metrics", "metrics": [{"label": "Coarse-class agreement", "field": "coarse_agreement", "format": "percent"}, {"label": "Null mean", "field": "null_agreement_mean", "format": "percent"}]},
        {"id": "kappa_card", "description": "校正類別邊際分布後的一致性。", "dataset": "verdict_metrics", "sourceId": "dataset_verdict_metrics", "metrics": [{"label": "Cohen κ", "field": "cohen_kappa", "format": "number"}, {"label": "Null κ 97.5%", "field": "null_kappa_q975", "format": "number"}]},
        {"id": "tree_card", "description": "候選 exact-tree 集合 digest；比五分類更細粒度。", "dataset": "verdict_metrics", "sourceId": "dataset_verdict_metrics", "metrics": [{"label": "Ordered exact-T agreement", "field": "ordered_exact_tree_agreement", "format": "percent"}, {"label": "Phase-swap tolerant", "field": "phase_swap_tree_agreement", "format": "percent"}]},
    ]

    heatmap_y = {"fields": [f"to_{category}" for category in CATEGORY_ORDER], "labels": {f"to_{category}": CATEGORY_ZH[category] for category in CATEGORY_ORDER}}
    genomic_y = {"fields": genomic_bin_fields, "labels": genomic_bin_labels}
    charts = [
        {
            "id": "all_dataset_composition_chart",
            "title": "7 dataset rows 的五類粗拓撲組成",
            "subtitle": "chr1-22 complete primary regions；每個 dataset 內歸一化為 100%。",
            "type": "stackedBar100",
            "dataset": "composition_long",
            "sourceId": "dataset_composition_long",
            "intent": "composition",
            "question": "七個 dataset rows 的五類粗拓撲份額是否相同？",
            "rationale": "100% stacked bars compare mutually exclusive class composition despite different complete-region denominators.",
            "comparisonContext": {"denominator": "complete primary regions within each dataset", "grain": "region", "normalization": "within-dataset 100%", "unit": "share"},
            "encodings": {"x": {"field": "sample", "type": "nominal", "label": "Dataset"}, "y": {"field": "count", "type": "quantitative", "label": "Complete regions"}, "color": {"field": "category_chart", "type": "nominal", "label": "Coarse topology"}},
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True, "title": "粗拓撲縮寫；完整名稱見下表"},
        },
        {
            "id": "hcc_share_delta_chart",
            "title": "HCC1395_DORADO 與 HCC1395 的組成差",
            "subtitle": "正值代表 DORADO 份額較高；單位為 percentage points。",
            "type": "bar",
            "dataset": "hcc_share_delta",
            "sourceId": "dataset_hcc_share_delta",
            "intent": "comparison",
            "question": "兩技術 dataset 各類份額的方向與幅度為何？",
            "rationale": "Signed percentage-point deltas expose offsetting composition shifts hidden by a single agreement rate.",
            "comparisonContext": {"baseline": "HCC1395", "denominator": "complete primary regions", "grain": "coarse class", "unit": "percentage points"},
            "encodings": {"x": {"field": "category", "type": "nominal", "label": "Coarse topology"}, "y": {"field": "DORADO_minus_HCC_pp", "type": "quantitative", "label": "DORADO - HCC (pp)"}},
            "valueFormat": "number",
            "unit": " pp",
            "palette": {"kind": "diverging", "midpoint": 0},
            "referenceLines": [{"axis": "y", "value": 0, "label": "無差異", "lineStyle": "dashed", "color": "neutral"}],
        },
        {
            "id": "confusion_heatmap_chart",
            "title": "Exact-coordinate 5x5 coarse-class confusion matrix",
            "subtitle": "每列以 HCC1395 類別歸一化；n=5,720 complete-both regions。",
            "type": "heatmap",
            "dataset": "confusion_heatmap",
            "sourceId": "dataset_confusion_heatmap",
            "intent": "relationship",
            "question": "HCC1395 的每種粗類別在 DORADO 中最常轉成哪一類？",
            "rationale": "A row-normalized heatmap preserves all 25 transition cells and makes asymmetric class changes visible.",
            "comparisonContext": {"denominator": "each HCC1395 source-class row", "grain": "exact-coordinate complete-both region", "normalization": "row share", "unit": "share"},
            "encodings": {"x": {"field": "HCC1395_class", "type": "nominal", "label": "HCC1395 class"}, "y": heatmap_y},
            "valueFormat": "percent",
            "palette": {"kind": "sequential"},
        },
        {
            "id": "chrom_agreement_chart",
            "title": "Exact-coordinate coarse-class agreement by chromosome",
            "subtitle": "chr1-22；分母為各染色體 complete-both matches。",
            "type": "bar",
            "dataset": "chrom_exact",
            "sourceId": "dataset_chrom_exact",
            "intent": "comparison",
            "question": "粗拓撲一致性是否只由少數染色體驅動？",
            "rationale": "Chromosome-level bars reveal spatial heterogeneity while retaining per-chromosome denominators in the source table.",
            "comparisonContext": {"baseline": f"whole-genome exact agreement {fmt_pct(exact_metrics['raw_agreement'])}", "denominator": "complete-both matches per chromosome", "grain": "chromosome", "unit": "agreement rate"},
            "encodings": {"x": {"field": "chrom", "type": "ordinal", "label": "Chromosome"}, "y": {"field": "coarse_agreement", "type": "quantitative", "label": "Agreement", "format": "percent"}},
            "valueFormat": "percent",
            "referenceLines": [{"axis": "y", "value": exact_metrics["raw_agreement"], "label": "全基因組", "lineStyle": "dashed", "color": "neutral"}],
        },
        {
            "id": "genomic_heatmap_chart",
            "title": "25 Mb genomic-bin coarse-class agreement",
            "subtitle": "Exact-coordinate complete-both pairs；空格表示該 bin 無可評估 pair。",
            "type": "heatmap",
            "dataset": "genomic_heatmap",
            "sourceId": "dataset_genomic_heatmap",
            "intent": "relationship",
            "question": "一致性在 chr1-22 哪些 25 Mb 區間較高或較低？",
            "rationale": "Fixed-width genomic bins show regional heterogeneity without implying gene-level resolution or biological causality.",
            "comparisonContext": {"denominator": "complete-both exact matches in each populated 25 Mb bin", "grain": "chromosome x 25 Mb bin", "unit": "agreement rate"},
            "encodings": {"x": {"field": "chrom", "type": "ordinal", "label": "Chromosome"}, "y": genomic_y},
            "valueFormat": "percent",
            "palette": {"kind": "sequential"},
        },
        {
            "id": "annotation_agreement_chart",
            "title": "Cancer-gene/drug annotation strata 的 coarse-class agreement",
            "subtitle": "Exact-coordinate complete-both pairs；present/absent 是相同座標共享的 region annotation。",
            "type": "bar",
            "dataset": "annotation_long",
            "sourceId": "dataset_annotation_long",
            "intent": "comparison",
            "question": "有癌症基因或藥物註解的 regions 是否顯示更高粗拓撲一致性？",
            "rationale": "Grouped present-versus-absent bars show descriptive differences while the adjacent table preserves sample sizes and kappa.",
            "comparisonContext": {"baseline": "annotation absent within the same feature", "denominator": "exact-coordinate complete-both pairs", "grain": "annotation stratum", "unit": "agreement rate"},
            "encodings": {"x": {"field": "feature", "type": "nominal", "label": "Region annotation"}, "y": {"field": "coarse_agreement", "type": "quantitative", "label": "Agreement", "format": "percent"}, "color": {"field": "stratum", "type": "nominal", "label": "Stratum"}},
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True, "title": "註解狀態"},
        },
    ]

    tables = [
        {"id": "composition_table_spec", "title": "7 dataset 五類粗拓撲完整數據", "subtitle": "Counts 以 complete primary regions 為分母；incomplete 另列不納入五分類。", "dataset": "composition_table", "sourceId": "dataset_composition_table", "density": "dense", "defaultSort": {"field": "complete_regions", "direction": "desc"}, "columns": columns(("sample", "Dataset", {"type": "text"}), ("primary_regions", "Primary", {"format": "number"}), ("complete_regions", "Complete", {"format": "number"}), ("incomplete_regions", "Incomplete", {"format": "number"}), *((f"{category}_n", CATEGORY_ZH[category], {"format": "number"}) for category in CATEGORY_ORDER))},
        {"id": "hcc_share_delta_table_spec", "title": "HCC pair 五類份額差", "subtitle": "DORADO - HCC；正值表示 DORADO 較高。", "dataset": "hcc_share_delta", "sourceId": "dataset_hcc_share_delta", "density": "spacious", "defaultSort": {"field": "absolute_gap_pp", "direction": "desc"}, "columns": columns(("category", "Coarse class", {"type": "text"}), ("HCC1395_share", "HCC1395", {"format": "percent"}), ("DORADO_share", "DORADO", {"format": "percent"}), ("DORADO_minus_HCC_pp", "Delta (pp)", {"format": "number", "movement": True}), ("absolute_gap_pp", "Absolute gap (pp)", {"format": "number"}))},
        {"id": "confusion_table_spec", "title": "Exact-coordinate 5x5 confusion cells", "subtitle": "25 cells；row share 以 HCC1395 source class 為分母。", "dataset": "confusion_long", "sourceId": "dataset_confusion_long", "density": "dense", "defaultSort": {"field": "n", "direction": "desc"}, "columns": columns(("HCC1395_class", "HCC1395", {"type": "text"}), ("DORADO_class", "DORADO", {"type": "text"}), ("n", "Regions", {"format": "number"}), ("row_total", "Row denominator", {"format": "number"}), ("row_share", "Row share", {"format": "percent"}), ("diagonal", "Same class", {"type": "text"}))},
        {"id": "chrom_table_spec", "title": "chr1-22 exact-coordinate agreement", "subtitle": "每列保留 matched/complete-both 分母與 agreement 分子。", "dataset": "chrom_exact", "sourceId": "dataset_chrom_exact", "density": "dense", "defaultSort": {"field": "coarse_agreement", "direction": "asc"}, "columns": columns(("chrom", "Chromosome", {"type": "text"}), ("matched_all", "Matched", {"format": "number"}), ("complete_both", "Complete both", {"format": "number"}), ("agreements", "Agreements", {"format": "number"}), ("coarse_agreement", "Agreement", {"format": "percent"}))},
        {"id": "genomic_bins_table_spec", "title": "所有有資料的 25 Mb bins", "subtitle": "Exact-coordinate complete-both pairs；座標為 GRCh38 region midpoint 所在 bin。", "dataset": "genomic_bins_long", "sourceId": "dataset_genomic_bins_long", "density": "dense", "defaultSort": {"field": "coarse_agreement", "direction": "asc"}, "columns": columns(("chrom", "Chromosome", {"type": "text"}), ("bin", "25 Mb bin", {"type": "text"}), ("complete_both", "Pairs", {"format": "number"}), ("agreements", "Agreements", {"format": "number"}), ("coarse_agreement", "Agreement", {"format": "percent"}))},
        {"id": "annotation_compare_table_spec", "title": "Annotation present vs absent 完整對照", "subtitle": "以 chromosome + global region-length decile conditional null 檢驗；仍不是因果或臨床 actionability 估計。", "dataset": "annotation_compare", "sourceId": "dataset_annotation_compare", "density": "dense", "defaultSort": {"field": "present_minus_absent_pp", "direction": "desc"}, "columns": columns(("feature", "Annotation", {"type": "text"}), ("present_n", "Present n", {"format": "number"}), ("present_agreement", "Present agreement", {"format": "percent"}), ("present_ci95_low", "95% CI low", {"format": "percent"}), ("present_ci95_high", "95% CI high", {"format": "percent"}), ("present_kappa", "Present κ", {"format": "number"}), ("absent_n", "Absent n", {"format": "number"}), ("absent_agreement", "Absent agreement", {"format": "percent"}), ("absent_kappa", "Absent κ", {"format": "number"}), ("present_minus_absent_pp", "Difference (pp)", {"format": "number", "movement": True}), ("conditional_null_q025_pp", "Null q2.5 (pp)", {"format": "number"}), ("conditional_null_q975_pp", "Null q97.5 (pp)", {"format": "number"}), ("conditional_p_two_sided", "Conditional p", {"format": "number"}), ("interpretation", "Boundary", {"type": "text"}))},
        {"id": "gene_join_table_spec", "title": "兩 dataset 的基因/藥物註解覆蓋", "subtitle": "Primary region denominator；GENCODE gene body 為主分析，promoter 只是次要 sensitivity。", "dataset": "gene_join_summary", "sourceId": "dataset_gene_join_summary", "density": "dense", "defaultSort": {"field": "sample", "direction": "asc"}, "columns": columns(("sample", "Dataset", {"type": "text"}), ("primary_regions", "Primary", {"format": "number"}), ("gene_body_regions", "Gene-body regions", {"format": "number"}), ("gene_body_coverage", "Gene-body coverage", {"format": "percent"}), ("CGC_body_regions", "CGC-body regions", {"format": "number"}), ("CGC_body_coverage", "CGC-body coverage", {"format": "percent"}), ("DGIdb_interaction_body_regions", "DGIdb regions", {"format": "number"}), ("DGIdb_interaction_body_coverage", "DGIdb coverage", {"format": "percent"}), ("DGIdb_approved_antineoplastic_body_regions", "Approved antineoplastic regions", {"format": "number"}), ("DGIdb_approved_antineoplastic_body_coverage", "Approved antineoplastic coverage", {"format": "percent"}))},
        {"id": "notable_loci_table_spec", "title": "HCC1395 重要 gene loci 的實際 region 狀態", "subtitle": "BRCA2/TBC1D16/ERBB2/MYC 同時保留 matched、incomplete 與 no-region 負例；不代表臨床 actionability。", "dataset": "notable_loci", "sourceId": "dataset_notable_loci", "density": "dense", "defaultSort": {"field": "gene", "direction": "asc"}, "columns": columns(("gene", "Gene", {"type": "text"}), ("HCC1395_region", "HCC1395 region", {"type": "text"}), ("HCC1395_status", "HCC1395 topology", {"type": "text"}), ("DORADO_region", "DORADO region", {"type": "text"}), ("DORADO_status", "DORADO topology", {"type": "text"}), ("best_match", "Best match", {"type": "text"}), ("reciprocal_overlap", "Reciprocal overlap", {"format": "percent"}), ("complete_both", "Complete both", {"type": "text"}), ("coarse_pair", "Coarse pair", {"type": "text"}), ("coarse_agree", "Agree", {"type": "text"}), ("interpretation", "Boundary", {"type": "text"}))},
        {"id": "gene_examples_table_spec", "title": "CGC/approved-antineoplastic gene-body region examples（top 30）", "subtitle": "Exact-coordinate complete-both region summary；drug claim names are DGIdb records, not mutation-specific treatment recommendations。", "dataset": "gene_locus_examples", "sourceId": "dataset_gene_locus_examples", "density": "dense", "defaultSort": {"field": "exact_complete_regions", "direction": "desc"}, "columns": columns(("gene", "Gene", {"type": "text"}), ("CGC", "CGC", {"type": "text"}), ("CGC_role", "CGC role", {"type": "text"}), ("CGC_tier", "Tier", {"format": "number"}), ("exact_complete_regions", "Exact regions", {"format": "number"}), ("category_agreement", "Agreement", {"format": "percent"}), ("cohen_kappa", "κ", {"format": "number"}), ("approved_antineoplastic_drugs_n", "Approved antineoplastic records", {"format": "number"}), ("approved_antineoplastic_claim_names", "Example DGIdb drug names", {"type": "text"}), ("region_examples", "Region examples", {"type": "text"}), ("claim_ceiling", "Claim ceiling", {"type": "text"}))},
        {"id": "pipeline_status_table_spec", "title": "Clean-v3 point-in-time status", "subtitle": f"Freeze cutoff {latest_status['as_of']}；producer、receipt closeout、canonical 與 sensitivity gates 分開顯示。", "dataset": "pipeline_status", "sourceId": "dataset_pipeline_status", "density": "spacious", "defaultSort": {"field": "as_of", "direction": "asc"}, "columns": columns(("as_of", "As of", {"type": "text"}), ("producer_terminal_pass", "Producer PASS", {"format": "number"}), ("producer_expected", "Expected", {"format": "number"}), ("producer_active", "Active", {"format": "number"}), ("aggregate_success_present", "Aggregate success", {"type": "text"}), ("producer_verification_present", "Producer verification", {"type": "text"}), ("receipt_closeout_status", "Receipt closeout", {"type": "text"}), ("receipt_closeout_error", "Closeout error", {"type": "text"}), ("continuation_status", "Continuation", {"type": "text"}), ("canonical_root_present", "Canonical root", {"type": "text"}), ("sensitivity_root_present", "Sensitivity root", {"type": "text"}), ("latest_complete_source", "Latest complete source", {"type": "text"}), ("scientific_release_allowed", "Scientific release", {"type": "text"}), ("verdict", "Verdict", {"type": "text"}))},
        {"id": "sensitivity_table_spec", "title": "區間 matching sensitivity", "subtitle": "Exact coordinate 與四種邊界寬鬆規則；五分類 agreement 維持約 69%。", "dataset": "pair_sensitivity", "sourceId": "dataset_pair_sensitivity", "density": "dense", "defaultSort": {"field": "matched_all", "direction": "asc"}, "columns": columns(("scenario", "Scenario", {"type": "text"}), ("matched_all", "Matched", {"format": "number"}), ("complete_both", "Complete both", {"format": "number"}), ("coarse_agreement", "Agreement", {"format": "percent"}), ("cohen_kappa", "κ", {"format": "number"}), ("macro_Jaccard", "Macro Jaccard", {"format": "percent"}), ("phase_swap_HP_signature", "Phase-swap HP signature", {"format": "percent"}), ("ordered_exact_tree", "Ordered exact T", {"format": "percent"}), ("phase_swap_exact_tree", "Phase-swap exact T", {"format": "percent"}), ("null_agreement_q975", "Null agree 97.5%", {"format": "percent"}), ("permutation_p_ge", "Empirical p", {"format": "number"}))},
        {"id": "vaf_sites_table_spec", "title": "Most-likely T 教學例：逐位點同-HP read-AF", "subtitle": f"{vaf_example['sample']} {vaf_example['region']} · HP{vaf_example['HP']} · CN={vaf_example['CN']}。", "dataset": "vaf_sites", "sourceId": "dataset_vaf_sites", "density": "spacious", "defaultSort": {"field": "position", "direction": "asc"}, "columns": columns(("site", "Site", {"type": "text"}), ("position", "Position", {"format": "number"}), ("REF_reads", "REF reads", {"format": "number"}), ("ALT_reads", "ALT reads", {"format": "number"}), ("read_AF", "Same-HP read-AF", {"format": "percent"}))},
        {"id": "vaf_candidates_table_spec", "title": "Most-likely T 教學例：exact candidates 排序", "subtitle": "Softmax 是 candidate-set 內相對權重，不是 calibrated posterior。", "dataset": "vaf_candidates", "sourceId": "dataset_vaf_candidates", "density": "spacious", "defaultSort": {"field": "softmax_relative_weight", "direction": "desc"}, "columns": columns(("candidate", "Exact T", {"type": "text"}), ("edges", "Edges", {"type": "text"}), ("shape_signature", "Shape", {"type": "text"}), ("score", "Ordering score", {"format": "number"}), ("softmax_relative_weight", "Relative weight", {"format": "percent"}), ("selected", "Top candidate", {"type": "text"}))},
        {"id": "vaf_pair_comparison_table_spec", "title": "HCC1395 兩技術 dataset：VAF-supported 推測樹／shape 一致性", "subtitle": "每列分母不同；ordered 保留 HP 標籤，swap tolerant 只容許全域 HP1↔HP2 交換。所有 VAF 結果都是推測，不是 accuracy。", "dataset": "vaf_pair_comparison", "sourceId": "dataset_vaf_pair_comparison", "density": "spacious", "defaultSort": {"field": "denominator", "direction": "desc"}, "columns": columns(("layer", "Comparison layer", {"type": "text"}), ("denominator", "Pair n", {"format": "number"}), ("coverage_of_5720", "Coverage of 5,720", {"format": "percent"}), ("ordered_n", "Ordered agree n", {"format": "number"}), ("ordered_agreement", "Ordered agreement", {"format": "percent"}), ("phase_swap_n", "Swap agree n", {"format": "number"}), ("phase_swap_agreement", "Swap-tolerant agreement", {"format": "percent"}), ("selection_basis", "Selection basis", {"type": "text"}), ("claim_ceiling", "Meaning / ceiling", {"type": "text"}))},
        {"id": "definitions_table_spec", "title": "五類粗拓撲的可執行定義", "subtitle": "唯一 topology 才進入前四類；Topo>1 一律保留未定。", "dataset": "definitions", "sourceId": "dataset_definitions", "density": "spacious", "defaultSort": {"field": "class", "direction": "asc"}, "columns": columns(("class", "Class", {"type": "text"}), ("rule", "Rule", {"type": "text"}), ("boundary", "Interpretive boundary", {"type": "text"}))},
        {"id": "claim_ceiling_table_spec", "title": "Verdict 與 claim ceiling", "subtitle": "技術再現性與生物正確性分開判定。", "dataset": "claim_ceiling", "sourceId": "dataset_claim_ceiling", "density": "spacious", "defaultSort": {"field": "layer", "direction": "asc"}, "columns": columns(("layer", "Evidence layer", {"type": "text"}), ("verdict", "Verdict", {"type": "text"}), ("evidence", "Evidence", {"type": "text"}), ("allowed_claim", "Safe claim", {"type": "text"}))},
        {"id": "source_inventory_table_spec", "title": "Cancer gene/drug source inventory", "subtitle": "主表顯示檔名、版本、build、grain 與 join key；完整絕對路徑保留在同一 frozen dataset 的 source view 與 source_inventory.tsv。", "dataset": "source_inventory", "sourceId": "dataset_source_inventory", "density": "dense", "defaultSort": {"field": "source_id", "direction": "asc"}, "columns": columns(("source_id", "Source", {"type": "text"}), ("source_file", "File", {"type": "text"}), ("version_or_snapshot", "Version", {"type": "text"}), ("genome_build", "Build", {"type": "text"}), ("grain", "Grain", {"type": "text"}), ("row_count", "Rows", {"format": "number"}), ("join_key", "Join key", {"type": "text"}), ("coverage_caveat", "Coverage caveat", {"type": "text"}), ("decision", "Use", {"type": "text"}), ("sha256", "SHA-256", {"type": "text"}))},
        {"id": "validation_table_spec", "title": "Report-side validation gates", "subtitle": "保守守恆、annotation strata 与 permutation 下限。", "dataset": "validation", "sourceId": "dataset_validation", "density": "spacious", "defaultSort": {"field": "check", "direction": "asc"}, "columns": columns(("check", "Check", {"type": "text"}), ("observed", "Observed", {"type": "text"}), ("status", "Status", {"type": "text"}))},
        {"id": "provenance_table_spec", "title": "Input artifacts and SHA-256", "subtitle": "核心報告輸入；repo 內使用相對路徑，repo 外資料保留絕對路徑以維持可追溯 identity。", "dataset": "provenance", "sourceId": "dataset_provenance", "density": "dense", "defaultSort": {"field": "artifact", "direction": "asc"}, "columns": columns(("artifact", "Artifact", {"type": "text"}), ("path", "Path", {"type": "text"}), ("bytes", "Bytes", {"format": "number"}), ("mtime", "Modified", {"type": "text"}), ("sha256", "SHA-256", {"type": "text"}))},
        *optional_table_specs,
    ]

    technical_summary_topology = f"""## 技術總結：區間高重疊，粗拓撲僅部分一致

- **區間層：** {exact_metrics['matched_all']:,} 個 exact-coordinate matches，覆蓋 HCC1395 的 **{fmt_pct(exact_metrics['match_coverage_a'])}** 與 DORADO 的 **{fmt_pct(exact_metrics['match_coverage_b'])}** primary regions。
- **五分類層：** {exact_metrics['complete_both']:,} 個 complete-both regions 的 agreement 為 **{fmt_pct(exact_metrics['raw_agreement'])}**（95% CI {fmt_pct(exact_metrics['raw_agreement_ci95_low'])}–{fmt_pct(exact_metrics['raw_agreement_ci95_high'])}），Cohen κ=**{exact_metrics['cohen_kappa']:.3f}**。染色體保留的 5,000 次 permutation null agreement 97.5% 上界只有 **{fmt_pct(exact_metrics['permutation_null']['agreement_q975'])}**（empirical p={exact_metrics['permutation_null']['agreement_p_ge']:.6f}），因此不是單純類別比例造成的偶合。
- **但不等價：** 五類 aggregate profile 的 total-variation distance 為 **{aggregate_tvd_pp:.2f} percentage points**；ordered exact-tree set agreement 僅 **{fmt_pct(exact_metrics['ordered_exact_tree_set_digest_agreement'])}**，容許 HP1/HP2 phase swap 後為 **{fmt_pct(exact_metrics['phase_swap_tolerant_exact_tree_set_digest_agreement'])}**。
- **VAF 推測層：** 在兩側原本都 `Topo>1`、且各自被 VAF 第一順位集合縮至單一 shape 的 {vaf_shape_rescue['denominator']:,} 對中，shape agreement 為 **{fmt_pct(vaf_shape_rescue['ordered_agreement'])}**／swap **{fmt_pct(vaf_shape_rescue['phase_swap_agreement'])}**；在兩側都實際使用 VAF 且各自唯一的 {vaf_exact_unique['denominator']:,} 個 mutation-labeled exact forests 中，agreement 為 **{fmt_pct(vaf_exact_unique['ordered_agreement'])}**／swap **{fmt_pct(vaf_exact_unique['phase_swap_agreement'])}**。前者沒有 mutation direction，後者仍只是 raw read-AF heuristic；兩者都不是 accuracy 或 biological truth。
"""

    technical_summary_annotation = f"""### 癌症基因與藥物註解不能當成 truth set

CGC gene-body-present regions 的 descriptive agreement 為 **{fmt_pct(cgc['present_agreement'])}**，比 absent 高 **{cgc['present_minus_absent_pp']:.2f} pp**，但 chromosome + region-length-decile conditional p=**{cgc['conditional_p_two_sided']:.4f}**；DGIdb approved-antineoplastic gene-body-present 為 **{fmt_pct(anti['present_agreement'])}**，比 absent 高 **{anti['present_minus_absent_pp']:.2f} pp**，conditional p=**{anti['conditional_p_two_sided']:.4f}**。COSMIC CLP all-status n={clp_all['present_n']} 的 p={clp_all['conditional_p_two_sided']:.4f}；confirmed-somatic 只有 n={clp_confirmed['present_n']}，agreement={fmt_pct(clp_confirmed['present_agreement'])}（95% CI {fmt_pct(clp_confirmed['present_ci95_low'])}–{fmt_pct(clp_confirmed['present_ci95_high'])}），p={clp_confirmed['conditional_p_two_sided']:.4f}。**7/7 annotation strata 的 two-sided p 都 >0.05**，因此只能做 context/face validity，不能證明 topology 正確或藥物可用。
"""

    verdict_body = """### 目前判定

**SHARE WITH CAVEATS / PARTIAL technical reproducibility；SCIENTIFIC NO-GO for proof-of-effectiveness.** HCC1395 與 HCC1395_DORADO 是同一生物 cell line 的兩個技術 sequencing/processing datasets，不是兩個獨立病人。現有結果支持「訊號高於 chance、工程內部自洽、部分技術再現」，但 historical layered-v2 不能證明 biological clone truth，也不能驗證因果藥物關係。"""

    blocks: list[dict[str, Any]] = [
        {"id": "title", "type": "markdown", "body": "# HCC1395 兩技術資料粗拓撲與癌症基因藥物一致性驗證"},
        {"id": "freeze_scope", "type": "markdown", "body": f"**PARTIAL FREEZE — as of {latest_status['as_of']}：clean-v3 producer {latest_status['producer']['terminal_pass']}/{len(latest_status['producer']['datasets'])} PASS、{latest_status['producer']['active']} active；aggregate success={'present' if latest_status['producer']['aggregate_success_present'] else 'absent'}，receipt closeout={latest_status['closeout']['status']}（{latest_status['closeout']['error_code']}），continuation={latest_status['continuation_execution']['status']}，canonical root={'present' if latest_status['canonical']['root_present'] else 'absent'}，sensitivity root={'present' if latest_status['sensitivity']['root_present'] else 'absent'}。** Producer 已完成，但 closeout 未通過且 clean tree 尚未建立；因此最新完整定量來源仍是 historical layered-v2 engineering snapshot，scientific release 不允許。", "sourceId": "dataset_pipeline_status"},
        {"id": "technical_summary", "type": "markdown", "body": technical_summary_topology, "sourceId": "dataset_verdict_metrics"},
        {"id": "technical_summary_annotation", "type": "markdown", "body": technical_summary_annotation, "sourceId": "dataset_annotation_compare"},
        {"id": "verdict", "type": "markdown", "body": verdict_body},
        {"id": "headline_metrics", "type": "metric-strip", "cardIds": ["interval_card", "agreement_card", "kappa_card", "tree_card"]},
        {"id": "composition_heading", "type": "markdown", "body": f"## 全 7 dataset 組成顯示 HCC pair 並非分布等價\n\nHCC1395_DORADO 相對 HCC1395 的 `Topo>1 未定` 份額上升 **{next(row['DORADO_minus_HCC_pp'] for row in hcc_share_delta if row['category_key']=='topology_multiple_unresolved'):.2f} pp**，`直系 only` 下降 **{abs(next(row['DORADO_minus_HCC_pp'] for row in hcc_share_delta if row['category_key']=='direct_only')):.2f} pp**。這個方向顯示 DORADO dataset 在候選 topology 可識別性上更容易留在未定，不應只看 69% agreement 就說「兩者一致」。", "sourceId": "dataset_hcc_share_delta"},
        {"id": "composition_chart", "type": "chart", "chartId": "all_dataset_composition_chart"},
        {"id": "share_delta_chart", "type": "chart", "chartId": "hcc_share_delta_chart"},
        {"id": "share_delta_table", "type": "table", "tableId": "hcc_share_delta_table_spec"},
        {"id": "composition_table", "type": "table", "tableId": "composition_table_spec"},
        {"id": "confusion_heading", "type": "markdown", "body": f"## Exact-coordinate regions 有中度粗類別一致，但非完全重建\n\n主對角線共 {sum(row['n'] for row in confusion_long if row['diagonal']):,}/{exact_metrics['complete_both']:,} regions；κ={exact_metrics['cohen_kappa']:.3f}。下圖以 HCC1395 source class 進行 row normalization，深色非對角格表示方向性轉移，不只是總 agreement 的殘差。", "sourceId": "dataset_confusion_long"},
        {"id": "confusion_chart", "type": "chart", "chartId": "confusion_heatmap_chart"},
        {"id": "confusion_table", "type": "table", "tableId": "confusion_table_spec"},
        {"id": "chrom_heading", "type": "markdown", "body": f"## chr1-22 皆有一致訊號，但區間異質性仍存在\n\n染色體層 agreement 範圍為 **{fmt_pct(min(row['coarse_agreement'] for row in chrom_exact))}–{fmt_pct(max(row['coarse_agreement'] for row in chrom_exact))}**。固定 25 Mb bins 再細分後，可以看出局部異質性；但 bin 只是顯示位置，不能將高/低格直接解釋成某基因機制。", "sourceId": "dataset_chrom_exact"},
        {"id": "chrom_chart", "type": "chart", "chartId": "chrom_agreement_chart"},
        {"id": "chrom_table", "type": "table", "tableId": "chrom_table_spec"},
        {"id": "genomic_chart_explain", "type": "markdown", "body": "### 25 Mb 基因組熱圖的讀法\n\n每格是該 chromosome×bin 內 complete-both exact pairs 的 coarse-class agreement。色階是 bin 間相對比較；格內 pair 數不同，因此必須連同下方 denominator table 解讀。"},
        {"id": "genomic_chart", "type": "chart", "chartId": "genomic_heatmap_chart"},
        {"id": "genomic_table", "type": "table", "tableId": "genomic_bins_table_spec"},
        {"id": "annotation_heading", "type": "markdown", "body": f"## Cancer-gene/drug strata 只提供 context/face validity，不提供正確性 truth\n\nCGC-present 與 approved-antineoplastic-present 的 agreement 各比 absent 高 **{cgc['present_minus_absent_pp']:.2f} pp** 與 **{anti['present_minus_absent_pp']:.2f} pp**，但 5,000 次 chromosome + region-length-decile conditional-null p 分別為 **{cgc['conditional_p_two_sided']:.4f}** 與 **{anti['conditional_p_two_sided']:.4f}**。CLP all-status 與 confirmed-somatic 也都未顯著（p={clp_all['conditional_p_two_sided']:.4f}/{clp_confirmed['conditional_p_two_sided']:.4f}）。因此目前沒有 annotation enrichment 證據；這些資料只能提供生物脈絡與後續優先排序。", "sourceId": "dataset_annotation_compare"},
        {"id": "annotation_chart", "type": "chart", "chartId": "annotation_agreement_chart"},
        {"id": "annotation_table", "type": "table", "tableId": "annotation_compare_table_spec"},
        {"id": "gene_join_table", "type": "table", "tableId": "gene_join_table_spec"},
        {"id": "gene_locus_heading", "type": "markdown", "body": "## Gene-region 實例同時呈現正例、未完整與 no-region 負例\n\nBRCA2 在 overlap≥80% sensitivity 下兩邊皆 complete 且都是 `Topo>1 未定`；這支持 locus-level coarse consistency，但不支持 exact-tree identity。TBC1D16 兩邊都 incomplete，ERBB2 只有 HCC1395 incomplete region，MYC 兩邊都沒有 region。這三個負例避免只展示有對上的癌症基因，也說明 gene/drug 名稱本身不能替代 topology eligibility。", "sourceId": "dataset_notable_loci"},
        {"id": "notable_loci", "type": "table", "tableId": "notable_loci_table_spec"},
        {"id": "gene_examples_explain", "type": "markdown", "body": "### CGC/DGIdb region examples 是 lookup，不是治療建議\n\n下表只保留有 exact complete region 且屬 CGC 或含 approved-antineoplastic DGIdb record 的前 30 個 gene summaries。Region keys 可供回查；drug names 只表示資料庫存在 gene–drug record，沒有驗證該 region 的 sSNV、方向、tumor type、dose 或臨床可用性。"},
        {"id": "gene_examples", "type": "table", "tableId": "gene_examples_table_spec"},
        *optional_blocks,
        {"id": "scope_heading", "type": "markdown", "body": "## 範圍與五分類定義\n\n範圍為 chr1-22、7 dataset rows 的 historical layered-v2 engineering snapshot。HCC1395/HCC1395_DORADO 共享一個 biological ID，但 BAM/VCF 與 processing provenance 不同；因此這是 technical replication，不是兩個獨立生物樣本。"},
        {"id": "definitions_table", "type": "table", "tableId": "definitions_table_spec"},
        {"id": "vaf_most_likely", "type": "markdown", "body": f"## 多顆 exact T 的 VAF-supported 第一順位如何判斷？\n\n在同一個 HP 內，先對每個 sSNV 位點由 reads 重算 `read-AF = ALT/(REF+ALT)`；在無回突變且沒有 CN/purity/multiplicity confounding 的理想序關係下，ancestor mutation 預期涵蓋 descendant 的細胞集合，因此以 **ancestor read-AF ≥ descendant read-AF** 作 ordering heuristic。每顆 exact T 的 parent→child AF 差值加總後，正式第一順位集合使用整數 read counts 轉成 `Fraction` 的 **exact score argmax**；exact score 完全相同者全部保留，不任意打破 tie。temperature={vaf_example['params']['temperature']} 的 softmax 只用來顯示 candidate-set 內相對濃度／敏感度，並非主要選擇規則。\n\n教學例 `{vaf_example['region']}`、HP{vaf_example['HP']} 有 {vaf_example['n_T']} 顆 exact T、{vaf_example['n_Topo']} 種 topology；兩個位點 read-AF 為 {vaf_sites[0]['read_AF']:.4f} 與 {vaf_sites[1]['read_AF']:.4f}，因此 `{vaf_example['winner']}` 位於最高 exact-score 集合，顯示用 relative weight 為 {vaf_example['winner_weight']:.6f}。這只能稱 **VAF-supported 推測候選**。Softmax weight 不是 calibrated posterior，raw read-AF 也受 CN、purity、mutation multiplicity 與 depth 影響。主五分類中的 **`Topo>1 未定` 不會因有 VAF 第一順位就被覆寫**；VAF-selected 結果另列為推測層。", "sourceId": "vaf_method"},
        {"id": "vaf_sites", "type": "table", "tableId": "vaf_sites_table_spec"},
        {"id": "vaf_candidates", "type": "table", "tableId": "vaf_candidates_table_spec"},
        {"id": "vaf_pair_heading", "type": "markdown", "body": f"## 兩個 HCC1395 的樹結構比對確實加入 VAF，但結果必須稱為推測\n\n未排名完整 exact-candidate set 的基準仍為 ordered **{fmt_pct(exact_metrics['ordered_exact_tree_set_digest_agreement'])}**／swap **{fmt_pct(exact_metrics['phase_swap_tolerant_exact_tree_set_digest_agreement'])}**，它不使用 VAF 選 winner。另行加入 VAF 後，兩側原始皆 `Topo>1` 且各自縮成單一 rooted-unlabeled shape 的 {vaf_shape_rescue['denominator']:,} 對，ordered／swap agreement 為 **{fmt_pct(vaf_shape_rescue['ordered_agreement'])}／{fmt_pct(vaf_shape_rescue['phase_swap_agreement'])}**；這只表示 branching skeleton 相同，不能判斷 mutation 祖先方向。更嚴格地，在兩側都實際使用 VAF 且各自唯一的 {vaf_exact_unique['denominator']:,} 個 mutation-labeled exact forests 中，ordered／swap agreement 為 **{fmt_pct(vaf_exact_unique['ordered_agreement'])}／{fmt_pct(vaf_exact_unique['phase_swap_agreement'])}**。Exact signature 包含每個 HP 的 genomic positions 與最高 exact-score candidate IDs，但仍是同批 reads 的 AF ordering heuristic，**不是獨立驗證、不是 posterior、也不是 biological tree accuracy**。", "sourceId": "dataset_vaf_pair_comparison"},
        {"id": "vaf_pair_comparison", "type": "table", "tableId": "vaf_pair_comparison_table_spec"},
        {"id": "method_heading", "type": "markdown", "body": "## 方法：先對齊區間，再比較互斥五分類\n\n1. 將每個 complete primary region 分為五個互斥類別；Topo>1 不強行選一個形狀。\n2. 主分析使用 chr/start/end 完全相同的 one-to-one pairs；sensitivity 使用 reciprocal overlap 80%/50% 與端點差 1/5 kb。\n3. Agreement 與 κ 的分母只含 complete-both pairs；不將 incomplete 誤當成第六個生物類別。\n4. 拓撲 null 在各 chromosome 內亂數排列 DORADO labels，固定 seed=20260712、5,000 permutations。\n5. VAF 推測層以 same-HP read-AF exact Fraction-score argmax 排序候選，保留所有並列；exact forest signature 含 HP-specific positions＋candidate IDs，shape signature 則移除 mutation labels。\n6. 原始 Topo>1 區域只有在所有 ambiguous units 都可評估時才允許 VAF single-shape selection；另報 ordered 與全域 HP1↔HP2 swap-tolerant sensitivity。\n7. Gene-body 以 GENCODE v46 GRCh38 1-based inclusive 座標重疊；CGC 以 HGNC bridge 為主，DGIdb 記錄區分 any/approved/approved-antineoplastic。\n8. Annotation test 固定 seed=20260712，以 chromosome + global region-length decile 保留各 stratum 的 present counts，做 5,000 次 conditional hypergeometric draws；目前尚未納入 CN/depth 或多重比較校正。"},
        {"id": "sensitivity_table", "type": "table", "tableId": "sensitivity_table_spec"},
        {"id": "robustness_heading", "type": "markdown", "body": f"## 穩健性與不確定性：matching 不是主要瓶頸，tree identity 才是\n\n邊界規則從 exact 放寬到 reciprocal overlap≥50%，matched regions 增加，但 coarse agreement 仍維持約 69%、κ 維持約 0.49，表示粗類別不一致並非只是小幅邊界漂移。VAF 可以把許多原 Topo>1 候選縮到一種 abstract shape，但在雙側實際 VAF 且唯一 exact forest 的 subset，swap-tolerant agreement 仍只有 **{fmt_pct(vaf_exact_unique['phase_swap_agreement'])}**。因此目前證據支持 VAF 作候選排序輔助，卻不支持「同一顆具體樹已被穩定確認」的強 claim。", "sourceId": "dataset_vaf_pair_comparison"},
        {"id": "claim_ceiling", "type": "table", "tableId": "claim_ceiling_table_spec"},
        {"id": "limitations_heading", "type": "markdown", "body": "## 限制：此報告無法證明 biological clone truth 或用藥效果\n\n- Historical layered-v2 是 engineering snapshot；clean layered-v3 全樣本 aggregate closeout 未納入。\n- 兩個 HCC dataset 共享生物樣本、caller/pipeline 方法家族與多個 upstream annotations，因此不是完全獨立 replication。\n- H_ 節點同時包含未觀測 intermediate 與 partial-supported completion，不能等同 hidden clone。\n- VAF exact-score 由同批 reads 的 same-HP AF 得到，屬內部排序證據，不是獨立驗證；raw AF 未經 purity/CN/multiplicity 校正，也沒有 molecule bootstrap。\n- Rooted-unlabeled shape 會移除 mutation labels，不能分辨 A→B 與 B→A；mutation-labeled unique exact forest 雖保留 genomic positions，仍只是 heuristic argmax。\n- Softmax relative weight 只是 candidate-set 內濃度，不是 calibrated posterior。\n- COSMIC CGC 是 cancer-gene census；DGIdb 是 gene-drug interaction aggregation。它們不是 mutation-specific clinical actionability，也不是 region topology ground truth。"},
        {"id": "validation_table", "type": "table", "tableId": "validation_table_spec"},
        {"id": "pipeline_status_table", "type": "table", "tableId": "pipeline_status_table_spec"},
        {"id": "source_inventory", "type": "table", "tableId": "source_inventory_table_spec"},
        {"id": "next_steps", "type": "markdown", "body": "## 下一步：將 partial technical signal 升級為可驗證科學證據\n\n1. 以同一 locked classifier 重跑 clean layered-v3，產生 7/7 aggregate receipt 後重算本報告。\n2. 在更多「同生物樣本、獨立 library/basecaller」的 pairs 預註再現性門檻，不只依賴 HCC1395 一對。\n3. 用 single-cell DNA、multi-region sampling 或 synthetic spike-in truth 區分「技術一致」與「生物正確」。\n4. Annotation conditional null 已控制 chromosome + region-length decile；下一輪加入 CN + depth sensitivity 與多重比較校正。在此之前不使用 enrichment 或 actionability 說法。"},
        {"id": "further_questions", "type": "markdown", "body": "## 尚待回答的問題\n\n- DORADO 的 `Topo>1 未定` 較高是由 read depth、basecalling error、HP composition、k 或 CN state 何者主導？\n- 染色體/25 Mb 低一致 bins 在排除低深度與高 k 後是否仍存在？\n- 對 VAF-supported winner 加入 purity/CN/multiplicity 校正後，exact-tree agreement 能否提升且通過外部 truth？"},
        {"id": "provenance", "type": "table", "tableId": "provenance_table_spec"},
    ]

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "HCC1395 兩技術資料粗拓撲與癌症基因藥物一致性驗證",
            "description": "chr1-22 全 7 dataset 五分類 census、HCC1395 technical replication、VAF-supported 推測樹／shape 比對、genomic-bin concordance 與 cancer-gene/drug descriptive annotation。",
            "generatedAt": generated_at,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "partial",
            "datasets": frozen_datasets,
            # `accessIssues` is reserved for failed source queries in the
            # portable reader. Scientific-release and truth-set blockers are
            # rendered explicitly in the first report blocks instead.
            "accessIssues": [],
        },
        "sources": sources,
    }

    git = git_context()
    md = [
        "<!--\n"
        f"建立時間: {generated_at}\n"
        "目標: 全 chr1-22 / 7 dataset 粗拓撲 census，驗證 HCC1395 兩技術 dataset 一致性並加入癌症基因/藥物描述註解\n"
        "處理範圍: Task B comprehensive validation; historical layered-v2; PARTIAL / SCIENTIFIC NO-GO\n"
        "服務目標: G4 / G5（G3 annotation support）\n"
        f"Git: branch={git['branch']}; commit={git['commit']}; dirty={git['dirty']}\n"
        f"關聯檔案: InterSubMod/{args.analysis.as_posix()}\n"
        "-->",
        "# HCC1395 兩技術資料粗拓撲與癌症基因藥物一致性驗證",
        "用 SCQA + hypothesis-validation：先給判定，再用區間、粗五分類、exact T 與 gene/drug strata 分層驗證。",
        f"**PARTIAL FREEZE — as of {latest_status['as_of']}：producer {latest_status['producer']['terminal_pass']}/{len(latest_status['producer']['datasets'])} PASS；receipt closeout={latest_status['closeout']['status']} ({latest_status['closeout']['error_code']})；continuation={latest_status['continuation_execution']['status']}；canonical/sensitivity absent。**",
        technical_summary_topology,
        technical_summary_annotation,
        verdict_body,
        "## 全 7 dataset 五類組成",
        md_table(composition_table, [("sample", "Dataset"), ("primary_regions", "Primary"), ("complete_regions", "Complete"), ("incomplete_regions", "Incomplete"), *[(f"{category}_n", CATEGORY_ZH[category]) for category in CATEGORY_ORDER]]),
        "## HCC1395 pair 五類份額差",
        md_table(hcc_share_delta, [("category", "Class"), ("HCC1395_share", "HCC share"), ("DORADO_share", "DORADO share"), ("DORADO_minus_HCC_pp", "Delta pp")]),
        "## Exact-coordinate 5x5 confusion matrix",
        md_table(confusion_long, [("HCC1395_class", "HCC1395"), ("DORADO_class", "DORADO"), ("n", "n"), ("row_total", "row n"), ("row_share", "row share")]),
        "## chr1-22 agreement",
        md_table(chrom_exact, [("chrom", "Chrom"), ("matched_all", "Matched"), ("complete_both", "Complete both"), ("agreements", "Agree"), ("coarse_agreement", "Agreement")]),
        "## 25 Mb genomic bins",
        md_table(genomic_bins_long, [("chrom", "Chrom"), ("bin", "Bin"), ("complete_both", "Pairs"), ("agreements", "Agree"), ("coarse_agreement", "Agreement")]),
        "## Annotation-stratified agreement",
        md_table(annotation_compare, [("feature", "Annotation"), ("present_n", "Present n"), ("present_agreement", "Present agreement"), ("present_kappa", "Present kappa"), ("absent_n", "Absent n"), ("absent_agreement", "Absent agreement"), ("absent_kappa", "Absent kappa"), ("present_minus_absent_pp", "Delta pp")]),
        "## Gene/drug join coverage",
        md_table(gene_join_summary, [("sample", "Dataset"), ("primary_regions", "Primary"), ("gene_body_regions", "Gene body"), ("gene_body_coverage", "Gene coverage"), ("CGC_body_regions", "CGC"), ("CGC_body_coverage", "CGC coverage"), ("DGIdb_interaction_body_regions", "DGIdb"), ("DGIdb_interaction_body_coverage", "DGIdb coverage"), ("DGIdb_approved_antineoplastic_body_regions", "Approved anti-neoplastic"), ("DGIdb_approved_antineoplastic_body_coverage", "Coverage")]),
        "## Notable gene loci",
        md_table(notable_loci, [("gene", "Gene"), ("HCC1395_region", "HCC1395 region"), ("HCC1395_status", "HCC1395 status"), ("DORADO_region", "DORADO region"), ("DORADO_status", "DORADO status"), ("best_match", "Best match"), ("reciprocal_overlap", "Overlap"), ("complete_both", "Complete both"), ("coarse_pair", "Coarse pair"), ("coarse_agree", "Agree")]),
        "## CGC/DGIdb region examples (non-actionable)",
        md_table(gene_locus_examples, [("gene", "Gene"), ("CGC", "CGC"), ("CGC_role", "Role"), ("CGC_tier", "Tier"), ("exact_complete_regions", "Exact regions"), ("category_agreement", "Agreement"), ("approved_antineoplastic_drugs_n", "Approved anti-neoplastic records"), ("approved_antineoplastic_claim_names", "Example drug names"), ("region_examples", "Regions"), ("claim_ceiling", "Claim ceiling")]),
        "## Matching sensitivity",
        md_table(pair_sensitivity, [("scenario", "Scenario"), ("matched_all", "Matched"), ("complete_both", "Complete both"), ("coarse_agreement", "Agreement"), ("cohen_kappa", "kappa"), ("macro_Jaccard", "Macro Jaccard"), ("ordered_exact_tree", "Ordered exact T"), ("phase_swap_exact_tree", "Swap tolerant exact T"), ("null_agreement_q975", "Null q97.5"), ("permutation_p_ge", "p")]),
        "## 五分類定義",
        md_table(definitions, [("class", "Class"), ("rule", "Rule"), ("boundary", "Boundary")]),
        "## VAF-supported most-likely exact T（heuristic only）",
        "正式第一順位使用 same-HP read-AF 的 exact Fraction-score argmax，最高分並列全部保留；softmax relative weight 只作濃度敏感度，不是 posterior。",
        md_table(vaf_sites, [("site", "Site"), ("position", "Position"), ("REF_reads", "REF"), ("ALT_reads", "ALT"), ("read_AF", "read-AF")]),
        md_table(vaf_candidates, [("candidate", "Exact T"), ("edges", "Edges"), ("shape_signature", "Shape"), ("score", "Score"), ("softmax_relative_weight", "Relative weight"), ("selected", "Selected")]),
        "## HCC1395 pair 的 VAF-supported 推測樹／shape 比對",
        f"兩側原始皆 Topo>1 且各自被 VAF 縮成單一 shape 的 {vaf_shape_rescue['denominator']:,} 對，ordered/swap agreement={fmt_pct(vaf_shape_rescue['ordered_agreement'])}/{fmt_pct(vaf_shape_rescue['phase_swap_agreement'])}；這只比較無 mutation labels 的 branching skeleton。兩側都實際使用 VAF 且各自唯一的 {vaf_exact_unique['denominator']:,} 個 mutation-labeled exact forests，ordered/swap agreement={fmt_pct(vaf_exact_unique['ordered_agreement'])}/{fmt_pct(vaf_exact_unique['phase_swap_agreement'])}。兩者都只是 read-AF heuristic 推測，不是 accuracy、posterior 或 biological truth。",
        md_table(vaf_pair_comparison, [("layer", "Layer"), ("denominator", "Pair n"), ("coverage_display", "Coverage"), ("ordered_n", "Ordered agree n"), ("ordered_display", "Ordered agreement"), ("phase_swap_n", "Swap agree n"), ("phase_swap_display", "Swap agreement"), ("selection_basis", "Selection basis"), ("claim_ceiling", "Claim ceiling")]),
        "## Verdict / claim ceiling",
        md_table(claim_ceiling, [("layer", "Layer"), ("verdict", "Verdict"), ("evidence", "Evidence"), ("allowed_claim", "Allowed claim")]),
        "## Source inventory",
        md_table(source_inventory, [("source_id", "Source"), ("source_file", "File"), ("absolute_path", "Absolute path"), ("version_or_snapshot", "Version"), ("genome_build", "Build"), ("row_count", "Rows"), ("join_key", "Join key"), ("coverage_caveat", "Caveat"), ("decision", "Use"), ("sha256", "SHA-256")]),
        "## Validation",
        md_table(validation_rows, [("check", "Check"), ("observed", "Observed"), ("status", "Status")]),
        "## Clean-v3 point-in-time status",
        md_table(pipeline_status, [("as_of", "As of"), ("producer_terminal_pass", "Producer PASS"), ("producer_expected", "Expected"), ("producer_active", "Active"), ("aggregate_success_present", "Aggregate"), ("canonical_root_present", "Canonical"), ("sensitivity_root_present", "Sensitivity"), ("scientific_release_allowed", "Scientific release"), ("verdict", "Verdict")]),
        "## Provenance",
        md_table(provenance, [("artifact", "Artifact"), ("path", "Path"), ("bytes", "Bytes"), ("mtime", "Modified"), ("sha256", "SHA-256")]),
        "---\n\nPARTIAL — historical layered-v2 engineering snapshot。VAF-selected 樹／shape 是同批 reads 的 heuristic 推測，不是 posterior、accuracy 或 biological truth；SCIENTIFIC NO-GO for proof-of-effectiveness；cancer-gene/drug annotations 不是 topology truth 或用藥建議。",
    ]

    args.report_md.parent.mkdir(parents=True, exist_ok=True)
    args.artifact.parent.mkdir(parents=True, exist_ok=True)
    args.report_md.write_text("\n\n".join(md) + "\n", encoding="utf-8")
    args.artifact.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"INPUT analysis -> {args.analysis}")
    print(f"INPUT pair TSVs -> {args.summary}, {args.confusion}, {args.matches}, {args.chrom}, {args.match_metrics}")
    print(f"INPUT annotation -> {args.gene_profile}, {args.annotation_agreement}, {args.source_inventory}")
    print(f"INPUT VAF pair -> {args.vaf_pair_metrics}, {args.vaf_pair_summary}, {args.vaf_pair_regions}, {args.vaf_pair_checks}")
    print(f"OUTPUT markdown -> {args.report_md}")
    print(f"OUTPUT artifact -> {args.artifact}")
    print(f"CHECK topology -> {analysis['validation']['passed']}/{analysis['validation']['checks']} PASS")
    print(f"CHECK exact complete-both -> confusion={sum(as_int(row['n']) for row in confusion_exact)} metrics={exact_metrics['complete_both']}")
    print(f"VERDICT -> SHARE WITH CAVEATS / PARTIAL; SCIENTIFIC NO-GO")


if __name__ == "__main__":
    main()
