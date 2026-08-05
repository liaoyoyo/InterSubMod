from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import importlib.util
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "build_pattern_methyl_report.py"
)
SPEC = importlib.util.spec_from_file_location("build_pattern_methyl_report", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
report = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = report
SPEC.loader.exec_module(report)


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.endswith(".gz"):
        handle_context = gzip.open(
            path, mode="wt", encoding="utf-8", newline=""
        )
    else:
        handle_context = path.open(mode="wt", encoding="utf-8", newline="")
    with handle_context as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if path.name.endswith(".gz"):
        handle_context = gzip.open(path, mode="rt", encoding="utf-8", newline="")
    else:
        handle_context = path.open(mode="rt", encoding="utf-8", newline="")
    with handle_context as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_row(
    *,
    dataset: str,
    chrom: str,
    region: str,
    unit: str,
    hp_raw: str,
    positions: list[int],
    states: dict[str, int],
    partial: dict[str, int] | None = None,
) -> dict[str, str]:
    partial = partial or {}
    row = {field: "" for field in sorted(report.PATTERN_COUNT_REQUIRED)}
    row.update(
        {
            "dataset": dataset,
            "chrom": chrom,
            "region_id": region,
            "unit_id": unit,
            "phase_set": f"PS-{unit}",
            "hp_family": hp_raw.split("-")[0],
            "hp_raw": hp_raw,
            "active_positions": ",".join(str(value) for value in positions),
            "n_active_bits": str(len(positions)),
            "n_complete": str(sum(states.values())),
            "n_partial": str(sum(partial.values())),
            "state_count_json": json.dumps(states, separators=(",", ":")),
            "partial_state_count_json": json.dumps(partial, separators=(",", ":")),
            "formal_n5": "true",
        }
    )
    return row


def evidence_row(
    count: dict[str, str],
    *,
    assessment: str,
    r2: float | None,
    q_by: float | None,
    best_pair: str,
    hamming: int,
    relation: str,
    evaluation_status: str = "EVALUABLE",
    invalid_reason: str = "",
    permanova_p: float = 0.002,
    permutations_requested: int = 49_999,
) -> dict[str, str]:
    row = {field: "" for field in sorted(report.EVIDENCE_REQUIRED)}
    row.update(
        {
            "schema_version": "1.0.0",
            "dataset": count["dataset"],
            "chrom": count["chrom"],
            "region_id": count["region_id"],
            "unit_id": count["unit_id"],
            "phase_set": count["phase_set"],
            "hp_family": count["hp_family"],
            "hp_raw": count["hp_raw"],
            "active_positions": count["active_positions"],
            "n_active_bits": count["n_active_bits"],
            "pair_full4": "false",
            "k_ge_3": "true" if int(count["n_active_bits"]) >= 3 else "false",
            "input_n_complete": count["n_complete"],
            "input_state_counts_json": count["state_count_json"],
            "analysis_n": count["n_complete"] if evaluation_status == "EVALUABLE" else "",
            "analysis_state_counts_json": (
                count["state_count_json"] if evaluation_status == "EVALUABLE" else "{}"
            ),
            "n_common_cpg": "4" if evaluation_status == "EVALUABLE" else "",
            "n_distal_cpg": "2" if evaluation_status == "EVALUABLE" else "",
            "permanova_r2": "" if r2 is None else str(r2),
            "permanova_p": str(permanova_p) if r2 is not None else "",
            "permanova_permutations_requested": (
                str(permutations_requested)
                if evaluation_status == "EVALUABLE"
                else ""
            ),
            "permdisp_p": (
                "0.01" if assessment == "CONFOUNDED" else "0.30"
            )
            if evaluation_status == "EVALUABLE"
            else "",
            "best_pair": best_pair if evaluation_status == "EVALUABLE" else "",
            "best_pair_hamming": str(hamming)
            if evaluation_status == "EVALUABLE"
            else "",
            "best_pair_distance_contrast": "0.22"
            if evaluation_status == "EVALUABLE"
            else "",
            "best_pair_standardized_effect": "0.80"
            if evaluation_status == "EVALUABLE"
            else "",
            "best_pair_topology_relation": relation
            if evaluation_status == "EVALUABLE"
            else "",
            "max_geometry_smd": (
                "0.72" if assessment == "CONFOUNDED" else "0.20"
            )
            if evaluation_status == "EVALUABLE"
            else "",
            "all_states_n8": "true",
            "all_states_n10": "true",
            "equal_n_retention": "0.80"
            if evaluation_status == "EVALUABLE"
            else "",
            "rarefaction_retention": "0.85"
            if evaluation_status == "EVALUABLE"
            else "",
            "distal_retention": "0.70"
            if evaluation_status == "EVALUABLE"
            else "",
            "multiplicity_family": "CONFIRMATORY_FULL4_OR_LONG",
            "q_by": "" if q_by is None else str(q_by),
            "p_holm": "0.02" if q_by is not None else "",
            "assessment": assessment,
            "evaluation_status": evaluation_status,
            "invalid_reason": invalid_reason,
        }
    )
    return row


def detail_row(
    count: dict[str, str],
    assessment: str,
    *,
    first: str,
    second: str,
    hamming: int,
    relation: str,
) -> dict[str, Any]:
    states = json.loads(count["state_count_json"])
    read_patterns = [
        state
        for state in sorted(states)
        for _ in range(states[state])
    ]
    read_order = [
        {
            "pattern": pattern,
            "qname_sha256": hashlib.sha256(
                f"{count['region_id']}:{index}".encode("utf-8")
            ).hexdigest(),
            "read_group": "RG0",
            "strand": "+" if index % 2 == 0 else "-",
        }
        for index, pattern in enumerate(read_patterns)
    ]
    distance_matrix = [
        [
            0.0
            if first_index == second_index
            else (
                0.10 + 0.01 * ((first_index + second_index) % 3)
                if read_patterns[first_index] == read_patterns[second_index]
                else 0.70 + 0.01 * ((first_index + second_index) % 3)
            )
            for second_index in range(len(read_patterns))
        ]
        for first_index in range(len(read_patterns))
    ]
    return {
        "schema_name": "intersubmod.pattern_methyl_evidence.detail",
        "schema_version": "1.0.0",
        "dataset": count["dataset"],
        "chrom": count["chrom"],
        "region_id": count["region_id"],
        "hp_raw": count["hp_raw"],
        "active_positions": [
            int(value) for value in count["active_positions"].split(",")
        ],
        "cpg_positions": [80, 120, 180, 220],
        "state_counts": states,
        "state_mean_profiles": {
            state: [0.10 + index * 0.25, 0.20, 0.75, 0.90 - index * 0.2]
            for index, state in enumerate(sorted(states))
        },
        "pairwise_effects": [
            {
                "first": first,
                "second": second,
                "hamming": hamming,
                "between_mean": 0.72,
                "pooled_within_mean": 0.31,
                "distance_contrast": 0.41,
                "standardized_effect": 0.90,
                "topology_relation": relation,
            }
        ],
        "topology": {
            "best_tree_unique": relation == "HAMMING1_GLOBAL_BEST_UNANIMOUS",
            "best_tree_tie_count": 1,
            "representative_best_edges": [
                {
                    "parent_vertex": 0,
                    "child_vertex": 1,
                    "parent_label": first,
                    "child_label": second,
                    "acquired_position": 180,
                }
            ],
        },
        "balanced_n_per_state": 20,
        "assessment": assessment,
        "read_order": read_order,
        "distance_matrix": distance_matrix,
    }


class SyntheticInputs:
    def __init__(self, root: Path):
        self.root = root
        self.census_receipt = root / "census_receipt.json"
        self.pattern_counts = root / "pattern_counts.tsv"
        self.catalog_summary = root / "artifact_catalog.json"
        self.catalog_tsv = root / "artifact_catalog.tsv"
        self.bernoulli_parity_summary = root / "bernoulli_parity_summary.json"
        self.evidence = root / "pattern_methyl_evidence.v1.tsv.gz"
        self.details = root / "pattern_methyl_details.v1.jsonl.gz"
        self.analysis_summary = root / "analysis_summary.json"
        self._write()

    def _write(self) -> None:
        counts = [
            count_row(
                dataset="H2009",
                chrom="chr5",
                region="chr5|PS=A|HP=2|US1:B0001",
                unit="US1",
                hp_raw="2-1",
                positions=[18096980, 18100000],
                states={"AA": 25, "AR": 20},
                partial={"AX": 5},
            ),
            count_row(
                dataset="HCC1395",
                chrom="chr22",
                region="chr22|PS=B|HP=1|US2:B0001",
                unit="US2",
                hp_raw="1-1",
                positions=[46257699, 46259000],
                states={"RA": 25, "AA": 20},
            ),
            count_row(
                dataset="H1437",
                chrom="chr1",
                region="chr1|PS=C|HP=1|UEXT:B0001",
                unit="UEXT",
                hp_raw="1-1",
                positions=[100, 200, 300],
                states={"RRR": 24, "AAA": 21},
            ),
            count_row(
                dataset="COLO829",
                chrom="chr2",
                region="chr2|PS=D|HP=2|UNULL:B0001",
                unit="UNULL",
                hp_raw="2",
                positions=[500, 600],
                states={"RR": 22, "AA": 20},
            ),
        ]
        write_tsv(
            self.pattern_counts,
            sorted(report.PATTERN_COUNT_REQUIRED),
            counts,
        )
        self.census_receipt.write_text(
            json.dumps(
                {
                    "schema_name": "intersubmod.exact_raw_hp_pattern_census",
                    "schema_version": "1.0.0",
                    "all_pass": True,
                    "scope": {
                        "task_type": "B_comprehensive_validation",
                        "datasets": [
                            "COLO829",
                            "H1437",
                            "H2009",
                            "HCC1395",
                            "HCC1395_DORADO",
                            "HCC1937",
                            "HCC1954",
                        ],
                        "chromosomes": [f"chr{index}" for index in range(1, 23)],
                    },
                    "counts": {
                        "pattern_count_rows": 4,
                        "formal_n5": 4,
                    },
                    "outputs": {
                        "pattern_counts": {
                            "sha256": sha256_path(self.pattern_counts),
                            "size_bytes": self.pattern_counts.stat().st_size,
                        }
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

        catalog_rows = [
            {
                "dataset": count["dataset"],
                "chrom": count["chrom"],
                "position1": str(position),
                "status": "PASS",
            }
            for count in counts
            for position in (
                int(value) for value in count["active_positions"].split(",")
            )
        ]
        write_tsv(self.catalog_tsv, sorted(report.CATALOG_REQUIRED), catalog_rows)
        self.catalog_summary.write_text(
            json.dumps(
                {
                    "schema_name": (
                        "intersubmod.multisite_pattern_methyl_artifact_catalog"
                    ),
                    "schema_version": "1.0.0",
                    "pass": True,
                    "status": "PASS",
                    "summary": {
                        "markers_total": len(catalog_rows),
                        "markers_pass": len(catalog_rows),
                        "markers_fail": 0,
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        self.bernoulli_parity_summary.write_text(
            json.dumps(
                {
                    "schema_name": "intersubmod.bernoulli_artifact_parity",
                    "schema_version": "1.0.0",
                    "all_pass": True,
                    "contract": {
                        "scope": "every formal marker",
                        "read_selection": (
                            "lowest SHA256(dataset,chrom,position1,read_id)"
                        ),
                        "max_reads_per_marker": 16,
                        "all_cpgs_retained": True,
                        "min_common_cpg": 3,
                        "absolute_tolerance": 0.0001,
                    },
                    "counts": {
                        "markers_total": len(catalog_rows),
                        "markers_pass": len(catalog_rows),
                        "markers_fail": 0,
                        "pair_cells_checked": 24,
                        "invalid_mask_mismatches": 0,
                    },
                    "max_absolute_error": 4.2e-8,
                    "inputs": {
                        "artifact_catalog": {
                            "path": str(self.catalog_tsv.resolve()),
                            "sha256": sha256_path(self.catalog_tsv),
                        }
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

        evidence_rows = [
            evidence_row(
                counts[0],
                assessment="CONFOUNDED",
                r2=0.19,
                q_by=0.03,
                best_pair="AA|AR",
                hamming=1,
                relation="HAMMING1_NOT_UNANIMOUS",
            ),
            evidence_row(
                counts[1],
                assessment="EVALUABLE_NO_ROBUST_ASSOCIATION",
                r2=0.04,
                q_by=0.72,
                best_pair="RA|AA",
                hamming=1,
                relation="HAMMING1_NOT_IN_GLOBAL_BEST",
                permanova_p=1.0 / 50_000.0,
            ),
            evidence_row(
                counts[2],
                assessment="ROBUST_ASSOCIATION",
                r2=0.271,
                q_by=0.004,
                best_pair="RRR|AAA",
                hamming=3,
                relation="PAIR_BAND_ONLY_HAMMING_GT1",
            ),
            evidence_row(
                counts[3],
                assessment="NOT_EVALUABLE",
                r2=None,
                q_by=None,
                best_pair="",
                hamming=0,
                relation="",
                evaluation_status="NOT_EVALUABLE",
                invalid_reason="insufficient_common_cpg",
            ),
        ]
        write_tsv(
            self.evidence,
            sorted(report.EVIDENCE_REQUIRED),
            evidence_rows,
        )
        details = [
            detail_row(
                counts[0],
                "CONFOUNDED",
                first="AA",
                second="AR",
                hamming=1,
                relation="HAMMING1_NOT_UNANIMOUS",
            ),
            detail_row(
                counts[2],
                "ROBUST_ASSOCIATION",
                first="RRR",
                second="AAA",
                hamming=3,
                relation="PAIR_BAND_ONLY_HAMMING_GT1",
            ),
        ]
        with gzip.open(self.details, "wt", encoding="utf-8", newline="") as handle:
            for detail in details:
                handle.write(
                    json.dumps(detail, ensure_ascii=False, separators=(",", ":"))
                    + "\n"
                )
        self.analysis_summary.write_text(
            json.dumps(
                {
                    "schema_name": (
                        "intersubmod.pattern_methyl_evidence.summary"
                    ),
                    "schema_version": "1.0.0",
                    "claim_ceiling": "input claim text is not copied verbatim",
                    "counts": {
                        "units_total": 4,
                        "units_evaluable": 3,
                        "detail_records": 2,
                        "assessment": {
                            "CONFOUNDED": 1,
                            "EVALUABLE_NO_ROBUST_ASSOCIATION": 1,
                            "NOT_EVALUABLE": 1,
                            "ROBUST_ASSOCIATION": 1,
                        },
                    },
                    "outputs": {
                        "evidence": {
                            "sha256": sha256_path(self.evidence),
                            "size_bytes": self.evidence.stat().st_size,
                        },
                        "details": {
                            "sha256": sha256_path(self.details),
                            "size_bytes": self.details.stat().st_size,
                        },
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    def rewrite_evidence(
        self, mutate: Any
    ) -> None:
        rows = read_tsv(self.evidence)
        mutate(rows)
        write_tsv(self.evidence, sorted(report.EVIDENCE_REQUIRED), rows)
        summary = json.loads(self.analysis_summary.read_text(encoding="utf-8"))
        summary["outputs"]["evidence"] = {
            "sha256": sha256_path(self.evidence),
            "size_bytes": self.evidence.stat().st_size,
        }
        self.analysis_summary.write_text(
            json.dumps(summary) + "\n", encoding="utf-8"
        )

    def rewrite_details(
        self, mutate: Any
    ) -> None:
        with gzip.open(self.details, "rt", encoding="utf-8") as handle:
            rows = [json.loads(line) for line in handle if line.strip()]
        mutate(rows)
        with gzip.open(self.details, "wt", encoding="utf-8", newline="") as handle:
            for row in rows:
                handle.write(
                    json.dumps(row, ensure_ascii=False, separators=(",", ":"))
                    + "\n"
                )
        summary = json.loads(self.analysis_summary.read_text(encoding="utf-8"))
        summary["outputs"]["details"] = {
            "sha256": sha256_path(self.details),
            "size_bytes": self.details.stat().st_size,
        }
        self.analysis_summary.write_text(
            json.dumps(summary) + "\n", encoding="utf-8"
        )

    def rewrite_catalog(self, mutate: Any) -> None:
        rows = read_tsv(self.catalog_tsv)
        mutate(rows)
        write_tsv(self.catalog_tsv, sorted(report.CATALOG_REQUIRED), rows)
        parity = json.loads(
            self.bernoulli_parity_summary.read_text(encoding="utf-8")
        )
        parity["inputs"]["artifact_catalog"]["sha256"] = sha256_path(
            self.catalog_tsv
        )
        self.bernoulli_parity_summary.write_text(
            json.dumps(parity) + "\n", encoding="utf-8"
        )

    def build_data(self) -> dict[str, Any]:
        return report.build_report_data(
            census_receipt_path=self.census_receipt,
            pattern_counts_path=self.pattern_counts,
            catalog_summary_path=self.catalog_summary,
            catalog_tsv_path=self.catalog_tsv,
            bernoulli_parity_summary_path=self.bernoulli_parity_summary,
            evidence_path=self.evidence,
            details_path=self.details,
            analysis_summary_path=self.analysis_summary,
        )


class PatternMethylReportTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.inputs = SyntheticInputs(self.root / "inputs")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_builds_self_contained_four_layer_report_with_injected_values(self) -> None:
        output_html = self.root / "out" / "report.html"
        output_json = self.root / "out" / "report_data.json"
        data = self.inputs.build_data()
        report.write_outputs(data, output_html, output_json)
        html = output_html.read_text(encoding="utf-8")

        self.assertNotIn("https://", html)
        self.assertNotIn("http://", html)
        self.assertNotIn("<script src=", html)
        self.assertNotIn("<link rel=", html)
        self.assertNotIn("lineage", html.lower())
        self.assertNotIn("ancestry", html.lower())
        self.assertIn("X pattern 僅為 subcube evidence", html)
        self.assertIn("Hamming&gt;1 僅畫 pair band", html)
        self.assertIn("@media (max-width:700px)", html)
        self.assertIn("@media print", html)
        self.assertIn('data-testid="claim-ceiling"', html)
        self.assertIn('data-testid="direct-answer"', html)
        self.assertIn('data-testid="strongest-secondary"', html)
        self.assertIn('data-testid="task-scope"', html)
        self.assertIn('data-testid="bernoulli-parity"', html)
        for test_id in (
            "filter-dataset",
            "filter-chrom",
            "filter-hp",
            "filter-pattern",
            "filter-assessment",
            "aggregate-layer",
            "canonical-layer",
            "extreme-layer",
            "explained-layer",
            "profile-heatmap",
            "read-cluster-heatmap",
            "read-cluster-canvas",
            "state-pair-matrix",
            "topology-relation",
        ):
            self.assertIn(f'data-testid="{test_id}"', html)
        positions = [
            html.index('data-testid="aggregate-layer"'),
            html.index('data-testid="canonical-layer"'),
            html.index('data-testid="extreme-layer"'),
            html.index('data-testid="explained-layer"'),
        ]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("18,096,980", html)
        self.assertIn("46,257,699", html)
        self.assertIn("0.271", html)
        self.assertIn("0.004", html)

        match = re.search(
            r'<script id="report-data" type="application/json">(.*?)</script>',
            html,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(match)
        embedded = json.loads(match.group(1))
        external = json.loads(output_json.read_text(encoding="utf-8"))
        self.assertEqual(embedded, external)
        self.assertEqual(external["aggregate"]["markers_pass"], 9)
        self.assertEqual(external["aggregate"]["analysis_units"], 4)
        self.assertEqual(
            external["aggregate"]["full4_units"],
            sum(bool(row["pair_full4"]) for row in external["cases"]),
        )
        self.assertEqual(
            external["aggregate"]["k_ge_3_units"],
            sum(bool(row["k_ge_3"]) for row in external["cases"]),
        )
        self.assertIn("formal full-four", external["direct_answer"])
        self.assertIn("secondary", external["direct_answer"])
        for label in (
            "全域 Formal full-four",
            "全域 k≥3 formal units",
            "全域 Artifact markers",
            "全域 Partial read projections",
        ):
            self.assertIn(label, html)
        self.assertEqual(
            set(external["aggregate"]["dataset_formal_unit_counts"]),
            set(external["scope"]["datasets"]),
        )
        self.assertIn("ordinalMarkerX", html)
        self.assertIn("outside displayed CpG span", html)
        self.assertEqual(external["aggregate"]["resolution_limited_units"], 1)
        self.assertEqual(external["bernoulli_parity"]["pair_cells_checked"], 24)
        self.assertEqual(
            external["bernoulli_parity"]["invalid_mask_mismatches"], 0
        )
        self.assertEqual(
            external["bernoulli_parity"]["max_absolute_error"], 4.2e-8
        )
        self.assertEqual(
            external["sources"]["bernoulli_parity_summary"]["sha256"],
            sha256_path(self.inputs.bernoulli_parity_summary),
        )
        resolution_limited = [
            row for row in external["cases"] if row["resolution_limited"]
        ]
        self.assertEqual(len(resolution_limited), 1)
        self.assertEqual(resolution_limited[0]["dataset"], "HCC1395")
        self.assertEqual(
            resolution_limited[0]["assessment"],
            "EVALUABLE_NO_ROBUST_ASSOCIATION",
        )
        self.assertIn('data-testid="resolution-limited-flag"', html)
        self.assertEqual(
            external["aggregate"]["assessment_counts"]["ROBUST_ASSOCIATION"],
            1,
        )
        self.assertTrue(external["scope"]["full_task_b_scope"])
        self.assertTrue(
            all(sentinel["found_in_pattern_counts"] for sentinel in external["sentinels"])
        )
        clusters = [
            row["detail"]["read_cluster"]
            for row in external["cases"]
            if row["detail"] is not None
        ]
        self.assertEqual(external["aggregate"]["read_cluster_records"], 2)
        self.assertEqual(external["aggregate"]["read_cluster_unavailable_records"], 2)
        self.assertEqual(
            external["aggregate"]["read_cluster_records"]
            + external["aggregate"]["read_cluster_source_null_records"]
            + external["aggregate"]["read_cluster_evaluable_no_detail_records"]
            + external["aggregate"]["read_cluster_non_evaluable_records"],
            external["aggregate"]["analysis_units"],
        )
        self.assertEqual(external["aggregate"]["read_cluster_reads_total"], 90)
        self.assertEqual(external["aggregate"]["read_cluster_cells_total"], 4050)
        self.assertTrue(all(cluster["available"] for cluster in clusters))
        for cluster in clusters:
            matrix_bytes = base64.b64decode(cluster["matrix_u8_b64"])
            self.assertEqual(
                len(matrix_bytes), cluster["n_reads"] * cluster["n_reads"]
            )
            self.assertEqual(
                hashlib.sha256(matrix_bytes).hexdigest(),
                cluster["matrix_u8_sha256"],
            )
            self.assertEqual(
                len(base64.b64decode(cluster["pattern_codes_u8_b64"])),
                cluster["n_reads"],
            )
            self.assertEqual(
                len(cluster["dendrogram_merges"]), cluster["n_reads"] - 1
            )
            self.assertFalse(cluster["order_uses_pattern_labels"])
            self.assertEqual(
                cluster["order_tie_breaker"],
                "SOURCE_LEAF_ORDINAL_FOR_EXACT_DISTANCE_TIES",
            )
            self.assertTrue(cluster["read_group_levels_pseudonymized"])
            self.assertEqual(cluster["read_group_levels"], ["1"])
        serialized_output = output_json.read_text(encoding="utf-8")
        self.assertNotIn("distance_matrix", serialized_output)
        self.assertNotIn("qname_sha256", serialized_output)
        self.assertNotIn('"RG0"', serialized_output)
        self.assertNotIn("實線箭頭", html)
        self.assertNotIn('<path d="M490 90', html)
        self.assertIn("undirected", html)
        self.assertNotIn("representative_best_edges", output_json.read_text(encoding="utf-8"))
        self.assertNotIn("parent_label", output_json.read_text(encoding="utf-8"))
        self.assertIn("未通過既定混雜 gate", external["technical_headline"])
        self.assertNotIn("由已量測混雜解釋", external["technical_headline"])
        self.assertNotIn("受已量測混雜影響", external["technical_headline"])

    def test_fails_closed_when_analysis_assessment_counts_drift(self) -> None:
        payload = json.loads(
            self.inputs.analysis_summary.read_text(encoding="utf-8")
        )
        payload["counts"]["assessment"]["ROBUST_ASSOCIATION"] = 2
        self.inputs.analysis_summary.write_text(
            json.dumps(payload) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(
            report.ReportContractError, "assessment counts"
        ):
            self.inputs.build_data()

    def test_accepts_hash_bound_combined_analysis_summary(self) -> None:
        payload = json.loads(
            self.inputs.analysis_summary.read_text(encoding="utf-8")
        )
        payload["schema_name"] = (
            "intersubmod.pattern_methyl_evidence.combined_summary"
        )
        payload["counts"]["confirmatory_input_units"] = 3
        payload["counts"]["secondary_input_units"] = 1
        self.inputs.analysis_summary.write_text(
            json.dumps(payload) + "\n", encoding="utf-8"
        )
        data = self.inputs.build_data()
        self.assertEqual(data["aggregate"]["analysis_units"], 4)

    def test_resolution_limited_uses_frozen_adaptive_boundary(self) -> None:
        base = {
            "family": report.RESOLUTION_LIMITED_FAMILY,
            "permutations_requested": 49_999,
            "p_value": 1.0 / 50_000.0,
            "q_by": 0.0500001,
        }
        self.assertTrue(report.is_resolution_limited(**base))
        for changed in (
            {"q_by": 0.05},
            {"p_value": 2.0 / 50_000.0},
            {"permutations_requested": 999},
            {"family": "EXPLORATORY_OTHER"},
        ):
            candidate = {**base, **changed}
            self.assertFalse(report.is_resolution_limited(**candidate))

    def test_fails_closed_when_census_hash_binding_drifts(self) -> None:
        with self.inputs.pattern_counts.open("a", encoding="utf-8") as handle:
            handle.write("\n")
        with self.assertRaisesRegex(report.ReportContractError, "SHA mismatch"):
            self.inputs.build_data()

    def test_fails_closed_on_every_formal_evidence_identity_drift(self) -> None:
        mutations = (
            ("unit_id", "UNIT-DRIFT", "unit_id mismatch"),
            ("phase_set", "PS-DRIFT", "phase_set mismatch"),
            ("hp_family", "9", "hp_family mismatch"),
            (
                "active_positions",
                "18096980,18100001",
                "active_positions mismatch",
            ),
            ("n_active_bits", "3", "n_active_bits mismatch"),
            ("input_n_complete", "44", "input_n_complete mismatch"),
            (
                "input_state_counts_json",
                '{"AA":24,"AR":21}',
                "input_state_counts mismatch",
            ),
        )
        for index, (field, replacement, expected) in enumerate(mutations):
            with self.subTest(field=field):
                inputs = SyntheticInputs(self.root / f"identity-{index}")

                def mutate(rows: list[dict[str, str]]) -> None:
                    rows[0][field] = replacement

                inputs.rewrite_evidence(mutate)
                with self.assertRaisesRegex(
                    report.ReportContractError, expected
                ):
                    inputs.build_data()

    def test_fails_closed_when_catalog_has_same_count_but_wrong_marker_set(
        self,
    ) -> None:
        def mutate(rows: list[dict[str, str]]) -> None:
            rows[0]["position1"] = "18096981"

        self.inputs.rewrite_catalog(mutate)
        with self.assertRaisesRegex(
            report.ReportContractError, "does not exactly match"
        ):
            self.inputs.build_data()

    def test_fails_closed_on_detail_state_profile_and_pair_drift(self) -> None:
        mutations = (
            (
                "active_positions",
                lambda rows: rows[0].update(
                    {"active_positions": [18096980, 18100001]}
                ),
                "active_positions mismatch",
            ),
            (
                "state_counts",
                lambda rows: rows[0]["state_counts"].update({"AA": 24}),
                "state_counts mismatch",
            ),
            (
                "profile_keys",
                lambda rows: rows[0]["state_mean_profiles"].pop("AR"),
                "profile/state key mismatch",
            ),
            (
                "profile_width",
                lambda rows: rows[0]["state_mean_profiles"]["AA"].pop(),
                "profile width mismatch",
            ),
            (
                "pair_hamming",
                lambda rows: rows[0]["pairwise_effects"][0].update(
                    {"hamming": 2}
                ),
                "pair hamming mismatch",
            ),
            (
                "pair_relation",
                lambda rows: rows[0]["pairwise_effects"][0].update(
                    {"topology_relation": "HAMMING1_TOPOLOGY_UNAVAILABLE"}
                ),
                "best-pair hamming/relation mismatch",
            ),
        )
        for index, (label, mutate, expected) in enumerate(mutations):
            with self.subTest(label=label):
                inputs = SyntheticInputs(self.root / f"detail-{index}")
                inputs.rewrite_details(mutate)
                with self.assertRaisesRegex(
                    report.ReportContractError, expected
                ):
                    inputs.build_data()

    def test_read_cluster_is_deterministic_and_pattern_label_independent(
        self,
    ) -> None:
        first = self.inputs.build_data()
        second = self.inputs.build_data()
        first_clusters = [
            row["detail"]["read_cluster"]
            for row in first["cases"]
            if row["detail"] is not None
        ]
        second_clusters = [
            row["detail"]["read_cluster"]
            for row in second["cases"]
            if row["detail"] is not None
        ]
        self.assertEqual(first_clusters, second_clusters)

        matrix = [
            [0.0, 0.1, 0.7, 0.8],
            [0.1, 0.0, 0.75, 0.72],
            [0.7, 0.75, 0.0, 0.12],
            [0.8, 0.72, 0.12, 0.0],
        ]
        case = {
            "dataset": "TEST",
            "chrom": "chr1",
            "region_id": "U1",
            "hp_raw": "1-1",
            "analysis_n": 4,
            "n_active_bits": 2,
        }
        qnames = [
            hashlib.sha256(f"read-{index}".encode()).hexdigest()
            for index in range(4)
        ]

        def payload(patterns: list[str]) -> dict[str, Any]:
            return {
                "read_order": [
                    {
                        "qname_sha256": qnames[index],
                        "pattern": pattern,
                        "read_group": "RG0",
                        "strand": "+" if index % 2 == 0 else "-",
                    }
                    for index, pattern in enumerate(patterns)
                ],
                "distance_matrix": matrix,
            }

        original = report.normalize_read_cluster(
            payload(["AA", "AA", "AR", "AR"]),
            case=case,
            state_counts={"AA": 2, "AR": 2},
        )
        relabeled = report.normalize_read_cluster(
            payload(["RR", "RR", "RA", "RA"]),
            case=case,
            state_counts={"RR": 2, "RA": 2},
        )
        self.assertEqual(
            original["dendrogram_merges"], relabeled["dendrogram_merges"]
        )
        self.assertEqual(original["matrix_u8_b64"], relabeled["matrix_u8_b64"])
        self.assertFalse(original["order_uses_pattern_labels"])
        self.assertNotEqual(original["pattern_levels"], relabeled["pattern_levels"])

    def test_read_cluster_matrix_contract_fails_closed(self) -> None:
        mutations = (
            (
                "asymmetric",
                lambda rows: rows[0]["distance_matrix"][0].__setitem__(1, 0.44),
                "asymmetric",
            ),
            (
                "non_zero_diagonal",
                lambda rows: rows[0]["distance_matrix"][0].__setitem__(0, 0.1),
                "non-zero diagonal",
            ),
            (
                "out_of_range",
                lambda rows: rows[0]["distance_matrix"][0].__setitem__(1, 1.1),
                "within 0..1",
            ),
            (
                "non_square",
                lambda rows: rows[0]["distance_matrix"][0].pop(),
                "width mismatch",
            ),
            (
                "duplicate_digest",
                lambda rows: rows[0]["read_order"][1].update(
                    {
                        "qname_sha256": rows[0]["read_order"][0][
                            "qname_sha256"
                        ]
                    }
                ),
                "duplicate qname digest",
            ),
            (
                "read_pattern_mass",
                lambda rows: rows[0]["read_order"][0].update({"pattern": "AR"}),
                "read_order/state count mismatch",
            ),
        )
        for index, (label, mutate, expected) in enumerate(mutations):
            with self.subTest(label=label):
                inputs = SyntheticInputs(self.root / f"cluster-contract-{index}")
                inputs.rewrite_details(mutate)
                with self.assertRaisesRegex(
                    report.ReportContractError, expected
                ):
                    inputs.build_data()

    def test_source_null_matrix_above_limit_is_explicit(self) -> None:
        patterns = ["AA"] * 81 + ["AR"] * 80
        payload = {
            "read_order": [
                {
                    "qname_sha256": hashlib.sha256(
                        f"large-{index}".encode()
                    ).hexdigest(),
                    "pattern": pattern,
                    "read_group": "RG0",
                    "strand": "+",
                }
                for index, pattern in enumerate(patterns)
            ],
            "distance_matrix": None,
        }
        case = {
            "dataset": "TEST",
            "chrom": "chr2",
            "region_id": "U-LARGE",
            "hp_raw": "2-1",
            "analysis_n": 161,
            "n_active_bits": 2,
        }
        cluster = report.normalize_read_cluster(
            payload,
            case=case,
            state_counts={"AA": 81, "AR": 80},
        )
        self.assertFalse(cluster["available"])
        self.assertEqual(cluster["source_status"], "SOURCE_NULL_N_GT_160")
        self.assertEqual(cluster["n_reads"], 161)

    def test_fails_closed_when_evidence_best_pair_metadata_drifts_from_detail(
        self,
    ) -> None:
        for index, (field, value) in enumerate(
            (
                ("best_pair_hamming", "2"),
                (
                    "best_pair_topology_relation",
                    "HAMMING1_TOPOLOGY_UNAVAILABLE",
                ),
            )
        ):
            with self.subTest(field=field):
                inputs = SyntheticInputs(self.root / f"best-pair-{index}")

                def mutate(rows: list[dict[str, str]]) -> None:
                    rows[0][field] = value

                inputs.rewrite_evidence(mutate)
                with self.assertRaisesRegex(
                    report.ReportContractError,
                    "best-pair hamming/relation mismatch",
                ):
                    inputs.build_data()

    def test_sentinel_prefers_fixed_exact_hp_and_discloses_truncation(
        self,
    ) -> None:
        target = report.SENTINEL_TARGETS[0]
        matches: list[dict[str, Any]] = []
        for index, hp_raw in enumerate(("1", "2", "2-2", "1-1", "2-1")):
            matches.append(
                {
                    "dataset": target["dataset"],
                    "chrom": target["chrom"],
                    "region_id": f"sentinel-{index}",
                    "unit_id": f"U{index}",
                    "phase_set": f"PS{index}",
                    "hp_family": hp_raw.split("-")[0],
                    "hp_raw": hp_raw,
                    "active_positions": [target["position"], target["position"] + 1],
                    "n_active_bits": 2,
                    "n_complete": 100 - index,
                    "n_partial": 0,
                    "state_counts": {"AA": 50, "RR": 50 - index},
                    "partial_state_counts": {},
                    "formal_n5": hp_raw != "2-1",
                }
            )
        sentinels = report.build_sentinels(matches, [])
        sentinel = next(row for row in sentinels if row["id"] == target["id"])
        self.assertEqual(sentinel["matches"][0]["hp_raw"], target["hp_raw"])
        self.assertTrue(sentinel["matches_truncated"])
        self.assertEqual(sentinel["displayed_match_count"], 3)
        self.assertEqual(sentinel["truncated_match_count"], 2)
        self.assertEqual(sentinel["preferred_hp_match_count"], 1)
        self.assertIn(
            'data-testid="sentinel-truncation"', report.HTML_TEMPLATE
        )

    def test_screen_reader_profile_table_is_layout_contained(self) -> None:
        self.assertIn(
            "contain:strict!important;table-layout:fixed!important",
            report.HTML_TEMPLATE,
        )
        self.assertIn("scroll-margin-top:132px", report.HTML_TEMPLATE)

    def test_fails_closed_when_parity_marker_count_drifts(self) -> None:
        payload = json.loads(
            self.inputs.bernoulli_parity_summary.read_text(encoding="utf-8")
        )
        payload["counts"]["markers_total"] = 3
        payload["counts"]["markers_pass"] = 3
        self.inputs.bernoulli_parity_summary.write_text(
            json.dumps(payload) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(
            report.ReportContractError, "does not match catalog"
        ):
            self.inputs.build_data()

    def test_fails_closed_when_parity_is_not_all_pass(self) -> None:
        payload = json.loads(
            self.inputs.bernoulli_parity_summary.read_text(encoding="utf-8")
        )
        payload["all_pass"] = False
        self.inputs.bernoulli_parity_summary.write_text(
            json.dumps(payload) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(
            report.ReportContractError, "not all_pass"
        ):
            self.inputs.build_data()

    def test_rejects_non_x_state_in_partial_subcube_evidence(self) -> None:
        with self.inputs.pattern_counts.open(
            "rt", encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows[0]["partial_state_count_json"] = '{"AA":5}'
        write_tsv(
            self.inputs.pattern_counts,
            sorted(report.PATTERN_COUNT_REQUIRED),
            rows,
        )
        receipt = json.loads(
            self.inputs.census_receipt.read_text(encoding="utf-8")
        )
        receipt["outputs"]["pattern_counts"]["sha256"] = sha256_path(
            self.inputs.pattern_counts
        )
        receipt["outputs"]["pattern_counts"]["size_bytes"] = (
            self.inputs.pattern_counts.stat().st_size
        )
        self.inputs.census_receipt.write_text(
            json.dumps(receipt) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(report.ReportContractError, "non-X state"):
            self.inputs.build_data()

    def test_refuses_to_overwrite_outputs(self) -> None:
        output_html = self.root / "report.html"
        output_json = self.root / "report.json"
        data = self.inputs.build_data()
        report.write_outputs(data, output_html, output_json)
        with self.assertRaisesRegex(report.ReportContractError, "overwrite"):
            report.write_outputs(data, output_html, output_json)


if __name__ == "__main__":
    unittest.main()
