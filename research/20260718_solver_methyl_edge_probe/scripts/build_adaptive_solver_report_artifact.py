#!/usr/bin/env python3
"""Build a canonical Data Analytics report artifact from the adaptive receipt."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


REPO_ROOT = Path(__file__).resolve().parents[3]
PROBE_ROOT = REPO_ROOT / "research" / "20260718_solver_methyl_edge_probe"
TITLE = "Hypercube solver：理論縮減與安全放寬評估"


def source(
    source_id: str,
    label: str,
    path: str,
    description: str,
    metric_definitions: List[str] | None = None,
) -> Dict[str, Any]:
    item: Dict[str, Any] = {
        "id": source_id,
        "label": label,
        "path": path,
        "query": {
            "engine": "verified-file",
            "description": description,
        },
    }
    if metric_definitions:
        item["query"]["metric_definitions"] = metric_definitions
    return item


def sql_literal(value: Any) -> str:
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, float)):
        return repr(value)
    return "'" + str(value).replace("'", "''") + "'"


def values_sql_source(
    source_id: str,
    label: str,
    rows: List[Dict[str, Any]],
    columns: List[str],
    description: str,
    metric_definitions: List[str],
) -> Dict[str, Any]:
    values = ",".join(
        "(" + ",".join(sql_literal(row[column]) for column in columns) + ")"
        for row in rows
    )
    column_sql = ",".join(columns)
    sql = (
        f"WITH rows({column_sql}) AS (VALUES {values}) "
        f"SELECT {column_sql} FROM rows"
    )
    return {
        "id": source_id,
        "label": label,
        "query": {
            "engine": "sqlite",
            "language": "sql",
            "sql": sql,
            "description": description,
            "metric_definitions": metric_definitions,
        },
    }


def build_artifact(receipt: Dict[str, Any]) -> Dict[str, Any]:
    generated_at = datetime.now(timezone.utc).isoformat()
    adaptive_rows = receipt["theory"]["adaptive_gate"]["rows"]
    fixture_rows = receipt["bounded_real_fixture_evidence"]["cases"]
    sample_rows = receipt["canonical_v5_context"]["sample_incomplete_rows"]
    tail = receipt["m2_frozen_v2_tail_diagnostic"]
    historical = receipt["historical_layered_v2_vaf_context"]
    current = receipt["canonical_v5_context"]

    sources = [
        source(
            "adaptive_receipt",
            "Adaptive solver limit assessment receipt",
            "InterSubMod/research/20260718_solver_methyl_edge_probe/results/adaptive_solver_limit_assessment_receipt.json",
            "公式、兩個bounded real fixtures、current-v5、historical VAF與M2 tail diagnostic的守恆評估；checks.all_pass=true。",
            [
                "proxy operations = 3^q*2^m + 2^q*m*2^(m-1)",
                "dense bytes = 2^q*2^m*16 in the planning gate",
                "proxy values are not measured runtime",
            ],
        ),
        source(
            "solver_probe",
            "Exact B&B bounded solver receipt",
            "InterSubMod/research/20260718_solver_methyl_edge_probe/results/solver_probe_receipt.json",
            "H2009_M31與COLO829_M31的B&B、persistent backend與完整optimal-set digest比較。",
            [
                "single-run timing is exploratory",
                "complete-set digest equality is the correctness evidence",
            ],
        ),
        source(
            "canonical_v5",
            "Current LongPhase-S PASS v5 topology summary",
            "InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json",
            "7 technical datasets／6 biological samples、chr1-22 current canonical summary；all_pass=true。",
            [
                "incomplete fraction = incomplete_regions / W_primary",
                "partial-only fraction = primary_units_partial_only / primary_units",
            ],
        ),
        source(
            "historical_vaf",
            "Historical layered-v2 VAF exact-top census",
            "InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_census.json",
            "Historical engineering snapshot；只用於ranking endpoint設計，不能與current-v5分母合併。",
            [
                "one-topology first-rank = unique_first + tied_first_same_topology",
                "selection is not independent biological confirmation",
            ],
        ),
        source(
            "m2_tail",
            "M2 frozen-v2 formal pilot NO-GO record",
            "InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md",
            "8小時timeout留下的diagnostic-only runtime；沒有完整ranking receipt。",
            [
                "7,058 is a parameter-grid unit-instance count, not unique biological regions",
                "tail-time shares are performance diagnostics only",
            ],
        ),
    ]

    tail_rows = [
        {
            "measure": "Tail units占全部units",
            "share": tail["tail_unit_fraction"],
            "numerator": tail["tail_candidate_only_units"],
            "denominator": tail["total_unit_instances"],
        },
        {
            "measure": "Tail units占candidate time",
            "share": tail["tail_candidate_time_fraction"],
            "numerator": tail["tail_candidate_seconds"],
            "denominator": tail["total_candidate_seconds"],
        },
    ]

    adaptive_dataset = [
        {
            "q": row["q"],
            "maximum_m": row["maximum_m"],
            "proxy_ops": row["proxy_ops_at_maximum_m"],
            "dense_mib": row["dense_mib_at_maximum_m"],
        }
        for row in adaptive_rows
    ]

    incomplete_dataset = [
        {
            "sample": row["sample"],
            "W_primary": row["W_primary"],
            "complete_regions": row["complete_regions"],
            "incomplete_regions": row["incomplete_regions"],
            "incomplete_fraction": row["incomplete_fraction"],
            "share_of_all_incomplete": row["share_of_all_incomplete"],
        }
        for row in sample_rows
    ]

    fixture_dataset = [
        {
            "case_id": row["case_id"],
            "raw_k": row["raw_k"],
            "effective_m": row["effective_m"],
            "input_groups": row["input_groups"],
            "reduced_groups": row["reduced_groups"],
            "naive_layers": row["naive_subset_layers_through_h_star"],
            "bnb_visited": row["bnb_visited_states"],
            "search_reduction": row["bnb_search_state_reduction_fraction"],
            "optimal_sets": row["optimal_vertex_sets"],
            "current_seconds": row["current_scipy_seconds"],
            "bnb_seconds": row["bnb_seconds"],
            "wall_ratio": row["single_run_wall_ratio_current_over_bnb"],
        }
        for row in fixture_rows
    ]

    policy_rows = [
        {
            "order": 1,
            "scope": "Production default",
            "limit": "effective m≤12",
            "status": "維持",
            "reason": "尚未完成33-tail與跨資料集dual-backend驗證",
        },
        {
            "order": 2,
            "scope": "Research pilot 1",
            "limit": "q≤6、m≤15",
            "status": "建議先試",
            "reason": "50M proxy gate內；先量actual RSS與p95/p99",
        },
        {
            "order": 3,
            "scope": "Research pilot 2",
            "limit": "q≤4、m≤17",
            "status": "條件式",
            "reason": "需pilot 1、packed memory與33-tail先PASS",
        },
        {
            "order": 4,
            "scope": "Small-q theory",
            "limit": "q≤2、m≤19",
            "status": "未核准",
            "reason": "只是operation/memory proxy，不是runtime證據",
        },
        {
            "order": 5,
            "scope": "Candidate family",
            "limit": "max_sets=256",
            "status": "不盲目提高",
            "reason": "改用objective／uniqueness／compressed family／full enumeration模式",
        },
    ]

    version_rows = [
        {
            "version": "Current canonical v5",
            "complete_regions": current["complete_regions"],
            "topology_multiple": current["topology_multiple_exact_multiple"],
            "vaf_ranked": "尚未以M2全量重算",
            "allowed_use": "current workload與stress targeting",
        },
        {
            "version": "Historical layered-v2",
            "complete_regions": historical["complete_regions"],
            "topology_multiple": 17909,
            "vaf_ranked": (
                f"{historical['first_rank_one_topology']}/"
                f"{historical['evaluable_ambiguous_regions']} first-rank one Topo"
            ),
            "allowed_use": "ranking endpoint設計；不可作current成功率",
        },
    ]

    headline = [
        {
            "production_m_cap": 12,
            "pilot_m_cap": 15,
            "pilot_q_cap": 6,
            "tail_unit_share": tail["tail_unit_fraction"],
            "tail_time_share": tail["tail_candidate_time_fraction"],
            "current_incomplete_share": current["incomplete_fraction"],
        }
    ]

    sources.extend(
        [
            values_sql_source(
                "headline_sql",
                "Reviewed adaptive-solver headline transformation",
                headline,
                [
                    "production_m_cap",
                    "pilot_m_cap",
                    "pilot_q_cap",
                    "tail_unit_share",
                    "tail_time_share",
                    "current_incomplete_share",
                ],
                "將adaptive receipt的已核對headline值物化為一列。",
                [
                    "tail shares are fractions on a 0-1 scale",
                    "production and pilot caps are policy values, not observed biological metrics",
                ],
            ),
            values_sql_source(
                "tail_concentration_sql",
                "Reviewed M2 tail-concentration transformation",
                tail_rows,
                ["measure", "share", "numerator", "denominator"],
                "將M2 frozen-v2 diagnostic-only unit與candidate-time占比物化為圖表rows。",
                [
                    "unit share = 33 / 7,058",
                    "candidate-time share = 14,715.729568 / 15,304.384699",
                ],
            ),
            values_sql_source(
                "adaptive_gate_sql",
                "Reviewed small-q adaptive-gate transformation",
                adaptive_dataset,
                ["q", "maximum_m", "proxy_ops", "dense_mib"],
                "依50M proxy operations與512 MiB planning gate物化q=1..10最大m。",
                [
                    "proxy_ops = 3^q*2^m + 2^q*m*2^(m-1)",
                    "dense_mib assumes 16 bytes per dense cell",
                ],
            ),
            values_sql_source(
                "incomplete_by_sample_sql",
                "Reviewed current-v5 incomplete-by-sample transformation",
                incomplete_dataset,
                [
                    "sample",
                    "W_primary",
                    "complete_regions",
                    "incomplete_regions",
                    "incomplete_fraction",
                    "share_of_all_incomplete",
                ],
                "從current-v5 canonical summary物化每technical dataset的incomplete rate與count。",
                [
                    "incomplete_fraction = incomplete_regions / W_primary",
                    "share_of_all_incomplete = incomplete_regions / 7,975",
                ],
            ),
            values_sql_source(
                "fixture_reduction_sql",
                "Reviewed bounded-fixture reduction transformation",
                fixture_dataset,
                [
                    "case_id",
                    "raw_k",
                    "effective_m",
                    "input_groups",
                    "reduced_groups",
                    "naive_layers",
                    "bnb_visited",
                    "search_reduction",
                    "optimal_sets",
                    "current_seconds",
                    "bnb_seconds",
                    "wall_ratio",
                ],
                "從solver receipts物化兩個bounded real fixtures的exact search與single-run timing。",
                [
                    "search_reduction = 1 - bnb_visited / naive_layers",
                    "wall_ratio = current SciPy seconds / B&B seconds",
                    "timing is exploratory; digest equality is the correctness gate",
                ],
            ),
            values_sql_source(
                "limit_policy_sql",
                "Reviewed solver-limit policy transformation",
                policy_rows,
                ["order", "scope", "limit", "status", "reason"],
                "把量化評估轉為production與research分層政策。",
                [
                    "policy rows are recommendations, not measured outcomes",
                    "any capped or deadline-limited enumeration remains complete=false",
                ],
            ),
            values_sql_source(
                "version_separation_sql",
                "Reviewed current/historical denominator separation",
                version_rows,
                [
                    "version",
                    "complete_regions",
                    "topology_multiple",
                    "vaf_ranked",
                    "allowed_use",
                ],
                "將current-v5與historical layered-v2分母、VAF資格與允許用途並列。",
                [
                    "historical VAF results must not be combined with current-v5 denominators",
                    "topology_multiple means T>1 and Topo>1 within each version contract",
                ],
            ),
        ]
    )

    cards = [
        {
            "id": "production_cap",
            "description": "Current production exact contract；在dual-pilot前不變。",
            "dataset": "headline",
            "sourceId": "headline_sql",
            "metrics": [
                {
                    "label": "Production effective m上限",
                    "field": "production_m_cap",
                    "format": "number",
                }
            ],
        },
        {
            "id": "first_pilot_cap",
            "description": "第一階段research relaxation，必須同時通過q與resource gate。",
            "dataset": "headline",
            "sourceId": "headline_sql",
            "metrics": [
                {
                    "label": "Pilot effective m上限",
                    "field": "pilot_m_cap",
                    "format": "number",
                },
                {
                    "label": "同時要求 q≤",
                    "field": "pilot_q_cap",
                    "format": "number",
                },
            ],
        },
        {
            "id": "tail_unit_share",
            "description": "M2 frozen-v2 chr6 diagnostic中的長尾unit占比。",
            "dataset": "headline",
            "sourceId": "headline_sql",
            "metrics": [
                {
                    "label": "Tail unit share",
                    "field": "tail_unit_share",
                    "format": "percent",
                }
            ],
        },
        {
            "id": "tail_time_share",
            "description": "同一diagnostic中長尾units占candidate-generation時間。",
            "dataset": "headline",
            "sourceId": "headline_sql",
            "metrics": [
                {
                    "label": "Tail candidate-time share",
                    "field": "tail_time_share",
                    "format": "percent",
                }
            ],
        },
    ]

    charts = [
        {
            "id": "tail_concentration_chart",
            "title": "M2 長尾 unit 與 candidate-generation time 占比",
            "subtitle": "HCC1395_DORADO chr6、structural-minread=1；7,058 parameter-grid unit instances，diagnostic only。",
            "intent": "comparison",
            "question": "少量長尾units是否主導candidate-generation時間？",
            "rationale": "兩個百分比直接呈現workload concentration；不將unit share與time share誤作同一分母。",
            "comparisonContext": {
                "baseline": "all 7,058 unit instances / all 15,304.385 candidate seconds",
                "denominator": "各bar使用自身明示分母",
                "grain": "performance diagnostic share",
                "unit": "fraction",
            },
            "type": "horizontalBar",
            "dataset": "tail_concentration",
            "sourceId": "tail_concentration_sql",
            "encodings": {
                "x": {
                    "field": "measure",
                    "type": "ordinal",
                    "label": "Measure",
                },
                "y": {
                    "field": "share",
                    "type": "quantitative",
                    "label": "Share",
                    "format": "percent",
                },
                "tooltip": [
                    {
                        "field": "numerator",
                        "type": "quantitative",
                        "label": "Numerator",
                        "format": "number",
                    },
                    {
                        "field": "denominator",
                        "type": "quantitative",
                        "label": "Denominator",
                        "format": "number",
                    },
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "labels": {"values": "all"},
        },
        {
            "id": "adaptive_gate_chart",
            "title": "Small-q subset-DP proxy gate：q 與最大 effective m",
            "subtitle": "50M proxy operations、512 MiB、16 B/cell；規劃門檻，不是runtime保證。",
            "intent": "comparison",
            "question": "不同reduced q下，哪些effective m值得進research pilot？",
            "rationale": "Bar chart顯示q增加時允許m下降的離散關係。",
            "comparisonContext": {
                "baseline": "fixed planning resource gate",
                "denominator": "one maximum m for each q",
                "grain": "theoretical planning row",
                "unit": "effective dimensions",
            },
            "type": "bar",
            "dataset": "adaptive_gate",
            "sourceId": "adaptive_gate_sql",
            "encodings": {
                "x": {"field": "q", "type": "ordinal", "label": "Reduced q"},
                "y": {
                    "field": "maximum_m",
                    "type": "quantitative",
                    "label": "Maximum effective m",
                    "format": "number",
                },
                "tooltip": [
                    {
                        "field": "proxy_ops",
                        "type": "quantitative",
                        "label": "Proxy operations",
                        "format": "number",
                    },
                    {
                        "field": "dense_mib",
                        "type": "quantitative",
                        "label": "Dense memory (MiB)",
                        "format": "number",
                    },
                ],
            },
            "valueFormat": "number",
            "layout": "full",
            "labels": {"values": "all"},
        },
        {
            "id": "incomplete_by_sample_chart",
            "title": "Current-v5 incomplete rate by technical dataset",
            "subtitle": "LongPhase-S recalibrated FILTER=PASS、chr1-22；incomplete/W_primary。",
            "intent": "comparison",
            "question": "哪兩個technical datasets最適合做solver stress panel？",
            "rationale": "排序horizontal bars同時保留rate與count tooltip，避免只看樣本量。",
            "comparisonContext": {
                "baseline": "each dataset W_primary",
                "denominator": "dataset-specific W_primary",
                "grain": "technical dataset × regional primary output",
                "unit": "fraction",
            },
            "type": "horizontalBar",
            "dataset": "incomplete_by_sample",
            "sourceId": "incomplete_by_sample_sql",
            "encodings": {
                "x": {
                    "field": "sample",
                    "type": "nominal",
                    "label": "Technical dataset",
                },
                "y": {
                    "field": "incomplete_fraction",
                    "type": "quantitative",
                    "label": "Incomplete rate",
                    "format": "percent",
                },
                "tooltip": [
                    {
                        "field": "incomplete_regions",
                        "type": "quantitative",
                        "label": "Incomplete regions",
                        "format": "number",
                    },
                    {
                        "field": "W_primary",
                        "type": "quantitative",
                        "label": "W_primary",
                        "format": "number",
                    },
                    {
                        "field": "share_of_all_incomplete",
                        "type": "quantitative",
                        "label": "Share of all incomplete",
                        "format": "percent",
                    },
                ],
            },
            "valueFormat": "percent",
            "layout": "full",
            "labels": {"values": "all"},
        },
    ]

    tables = [
        {
            "id": "fixture_reduction_table",
            "title": "兩個 bounded real fixtures 的 exact search縮減",
            "subtitle": "Objective、complete flag與完整optimal-set digest皆一致；timing為single-run exploratory。",
            "dataset": "fixture_reduction",
            "sourceId": "fixture_reduction_sql",
            "layout": "full",
            "defaultSort": {"field": "search_reduction", "direction": "desc"},
            "columns": [
                {"field": "case_id", "label": "Case", "type": "text"},
                {"field": "raw_k", "label": "Raw k", "type": "number"},
                {"field": "effective_m", "label": "Effective m", "type": "number"},
                {"field": "input_groups", "label": "Input groups", "type": "number"},
                {
                    "field": "reduced_groups",
                    "label": "Reduced groups",
                    "type": "number",
                },
                {"field": "naive_layers", "label": "Naive layers", "type": "number"},
                {"field": "bnb_visited", "label": "B&B visited", "type": "number"},
                {
                    "field": "search_reduction",
                    "label": "Search-state reduction",
                    "type": "percent",
                },
                {"field": "optimal_sets", "label": "Optimal V", "type": "number"},
                {"field": "wall_ratio", "label": "Wall ratio", "type": "number"},
            ],
        },
        {
            "id": "limit_policy_table",
            "title": "建議的限制政策",
            "subtitle": "Production與research pilot分層；理論gate不直接升級為production contract。",
            "dataset": "limit_policy",
            "sourceId": "limit_policy_sql",
            "layout": "full",
            "defaultSort": {"field": "order", "direction": "asc"},
            "columns": [
                {"field": "order", "label": "#", "type": "number"},
                {"field": "scope", "label": "Scope", "type": "text"},
                {"field": "limit", "label": "Limit", "type": "text"},
                {"field": "status", "label": "Status", "type": "text"},
                {"field": "reason", "label": "Reason", "type": "text"},
            ],
        },
        {
            "id": "version_separation_table",
            "title": "Current-v5 與 historical VAF 分母必須分開",
            "subtitle": "兩版本可回答不同問題；不得產生混合的最新成功率。",
            "dataset": "version_separation",
            "sourceId": "version_separation_sql",
            "layout": "full",
            "defaultSort": {"field": "complete_regions", "direction": "desc"},
            "columns": [
                {"field": "version", "label": "Version", "type": "text"},
                {
                    "field": "complete_regions",
                    "label": "Complete regions",
                    "type": "number",
                },
                {
                    "field": "topology_multiple",
                    "label": "T>1/Topo>1",
                    "type": "number",
                },
                {"field": "vaf_ranked", "label": "VAF status", "type": "text"},
                {"field": "allowed_use", "label": "Allowed use", "type": "text"},
            ],
        },
    ]

    blocks = [
        {"id": "title", "type": "markdown", "layout": "full", "body": f"# {TITLE}"},
        {
            "id": "technical_summary",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 技術結論：可做自適應放寬，但production先不改\n\n"
                "**可以放寬raw k，不能只把effective m上限改大。**先exact壓縮active bits與groups，"
                "再依effective m、reduced q、預估資源與optimal family大小V分流。"
                "Production維持m≤12；第一階段只pilot q≤6、m≤15；通過33-tail、H2009／H1437與"
                "packed-memory gate後，才考慮q≤4、m≤17。q≤2、m≤19只是理論proxy。"
                "\n\n**最大實際headroom集中在少數長尾。**33/7,058（0.468%）diagnostic units"
                "消耗96.154% candidate-generation time；應優先路由長尾，而非全面重寫。"
                "\n\n**演算法加速不等於生物拓撲確認。**Historical VAF的92.18%只能作ranking endpoint設計；"
                "current-v5與historical layered-v2分母必須分開。"
            ),
        },
        {
            "id": "headline_metrics",
            "type": "metric-strip",
            "layout": "full",
            "cardIds": [
                "production_cap",
                "first_pilot_cap",
                "tail_unit_share",
                "tail_time_share",
            ],
        },
        {
            "id": "tail_narrative",
            "type": "markdown",
            "layout": "full",
            "sourceId": "m2_tail",
            "body": (
                "## 0.468% 的長尾 units 主導 96.154% candidate-generation time\n\n"
                "圖中兩個bar使用不同但明示的分母：unit share是33/7,058；time share是"
                "14,715.730/15,304.385秒。這支持選擇性router與per-unit total deadline。"
                "若把tail cost理想化為0，candidate stage的純數學headroom約26.0×；這不是新backend預測。"
            ),
        },
        {
            "id": "tail_chart_block",
            "type": "chart",
            "layout": "full",
            "chartId": "tail_concentration_chart",
        },
        {
            "id": "adaptive_narrative",
            "type": "markdown",
            "layout": "full",
            "sourceId": "adaptive_receipt",
            "body": (
                "## q 比 raw k 更適合決定能否放寬 effective m\n\n"
                "Subset-DP規劃式同時對q與m指數成長。50M proxy gate下，q≤6可到m=15，"
                "q≤4可到m=17；但16 B/cell是樂觀packed假設，且objective proof不包含V個答案的完整輸出。"
                "因此圖是research routing map，不是production SLA。"
            ),
        },
        {
            "id": "adaptive_chart_block",
            "type": "chart",
            "layout": "full",
            "chartId": "adaptive_gate_chart",
        },
        {
            "id": "fixture_narrative",
            "type": "markdown",
            "layout": "full",
            "sourceId": "solver_probe",
            "body": (
                "## Exact B&B在兩個stress fixtures少走99.6%以上search states\n\n"
                "H2009與COLO829從206,368個naive subset layers降到814與531個visited states，"
                "仍完整輸出242／216個optimal sets並通過digest比對。這證明值得擴大pilot；"
                "兩例single-run的62.5–72.5× wall比值不能外推全樣本。"
            ),
        },
        {
            "id": "fixture_table_block",
            "type": "table",
            "layout": "full",
            "tableId": "fixture_reduction_table",
        },
        {
            "id": "sample_narrative",
            "type": "markdown",
            "layout": "full",
            "sourceId": "canonical_v5",
            "body": (
                "## H2009與H1437占current-v5全部incomplete的68.89%\n\n"
                "H2009 incomplete rate為39.05%，H1437為19.64%；兩者合計5,494/7,975。"
                "這支持把H2009列為第一個跨資料集stress panel、H1437為第二個。"
                "Current-v5 incomplete屬舊candidate/cap語意，只能作抽樣優先序，不能預測M2 timeout率。"
            ),
        },
        {
            "id": "sample_chart_block",
            "type": "chart",
            "layout": "full",
            "chartId": "incomplete_by_sample_chart",
        },
        {
            "id": "limit_policy_narrative",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 限制應改成路由政策，不是單一常數\n\n"
                "Raw k先做active-bit compression；effective m、q與resource prediction共同決定backend。"
                "Candidate cap不宜盲目提高，應拆成objective、uniqueness、exact count／compressed family及"
                "full enumeration模式。任何cap、deadline或BDD width限制都必須保留complete=false。"
            ),
        },
        {
            "id": "limit_policy_table_block",
            "type": "table",
            "layout": "full",
            "tableId": "limit_policy_table",
        },
        {
            "id": "version_narrative",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## Historical VAF可指引第三層排序，但不能填入current-v5主表\n\n"
                "Historical layered-v2的25,978/28,183=92.18%代表first-rank set只含一種Topo，"
                "不是current-v5或M2的最新成功率。2,205個跨Topo ties可作第三層證據target；"
                "不可評估者應abstain。另應將historical full-ALT group數與current candidate tree count"
                "分別命名為C_full_alt_groups與T_joint_candidates，停止共用C_region。"
            ),
        },
        {
            "id": "version_table_block",
            "type": "table",
            "layout": "full",
            "tableId": "version_separation_table",
        },
        {
            "id": "next_steps",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 下一步：先證明等價，再決定是否升級production\n\n"
                "1. 建current route census：raw_k、effective_m、q、h*、V、complete、wall與RSS。\n"
                "2. 用m≤5／可承受m≤6 exhaustive oracle驗證subset-DP objective 0 mismatch。\n"
                "3. 重播33-tail＋H2009／H1437；complete cases的sorted vertex-set digest必須0 mismatch。\n"
                "4. 先啟用q≤6、m≤15 research route；只有correctness、RSS、p95/p99全部過關才試m≤17。\n"
                "5. Strict infinite-sites維持CN/LOH-gated sensitivity；VAF、甲基與CCF只在完整候選後排序。"
            ),
        },
        {
            "id": "further_questions",
            "type": "markdown",
            "layout": "full",
            "body": (
                "## 尚未回答、會改變決策的問題\n\n"
                "- Current-v5／M2全量的effective m與reduced q實際分布為何？\n"
                "- 33-tail中多少是V巨大、多少是單次objective困難、多少可被ZDD壓縮？\n"
                "- Packed DP的actual RSS與traceback成本是否符合16 B/cell規劃？\n"
                "- Current-v5重跑VAF／likelihood後，historical 2,205個跨Topo tie target是否仍存在？"
            ),
        },
    ]

    manifest = {
        "version": 1,
        "surface": "report",
        "title": TITLE,
        "description": "Boolean-hypercube exact solver的理論縮減、實測長尾與自適應放寬決策。",
        "generatedAt": generated_at,
        "blocks": blocks,
        "cards": cards,
        "charts": charts,
        "tables": tables,
        "sources": sources,
    }
    snapshot = {
        "version": 1,
        "status": "ready",
        "generatedAt": generated_at,
        "datasets": {
            "headline": headline,
            "tail_concentration": tail_rows,
            "adaptive_gate": adaptive_dataset,
            "incomplete_by_sample": incomplete_dataset,
            "fixture_reduction": fixture_dataset,
            "limit_policy": policy_rows,
            "version_separation": version_rows,
        },
    }
    return {
        "surface": "report",
        "manifest": manifest,
        "snapshot": snapshot,
        "sources": sources,
        "package_info": {
            "scope": "PARTIAL_EXPLORATORY_NOT_PRODUCTION_VALIDATION",
            "audience": "InterSubMod PI and method developers",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--receipt",
        type=Path,
        default=PROBE_ROOT / "results" / "adaptive_solver_limit_assessment_receipt.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROBE_ROOT / "results" / "adaptive_solver_report_artifact.json",
    )
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text(encoding="utf-8"))
    if not receipt["checks"]["all_pass"]:
        raise RuntimeError("Adaptive solver receipt did not pass all checks")
    artifact = build_artifact(receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "PASS"
        f" blocks={len(artifact['manifest']['blocks'])}"
        f" charts={len(artifact['manifest']['charts'])}"
        f" tables={len(artifact['manifest']['tables'])}"
        f" datasets={len(artifact['snapshot']['datasets'])}"
        f" output={args.output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
