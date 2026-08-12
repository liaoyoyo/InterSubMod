#!/usr/bin/env python3
"""Build the canonical portable-report artifact from audited CSV/JSON outputs."""
from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUTPUT = ROOT / "artifact.json"


def read_csv(name: str) -> list[dict]:
    with (RESULTS / name).open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def num(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def integer(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def ratio_from_pct(value) -> float:
    return num(value) / 100.0


def main() -> int:
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    overview_raw = read_csv("bundle_overview.csv")
    coverage_raw = read_csv("methylation_coverage_by_k.csv")
    axis_raw = read_csv("methylation_axis_metrics.csv")
    cohort_raw = read_csv("cohort_topology_metrics.csv")
    audit = json.loads((RESULTS / "metrics_audit.json").read_text(encoding="utf-8"))
    legacy = json.loads((RESULTS / "legacy_browser_crosswalk.json").read_text(encoding="utf-8"))
    visual = json.loads(
        (RESULTS / "legacy_browser_visual_audit.generated.json").read_text(encoding="utf-8")
    )

    overview = []
    for row in overview_raw:
        overview.append({
            "bundle": row["bundle"],
            "built_at": row["built_at"],
            "status": row["validation_status_recomputed"],
            "regions": integer(row["regions"]),
            "tree_coverage": ratio_from_pct(row["tree_coverage_pct_all_regions"]),
            "unique_among_tree": ratio_from_pct(row["unique_pct_among_tree"]),
            "distinct_ssnv": integer(row["distinct_ssnv"]),
            "ism_linked": integer(row["methyl_topology_linked"]),
            "ism_linkage": ratio_from_pct(row["methyl_topology_linkage_pct"]),
            "selfcheck_fail": integer(row["selfcheck_fail"]),
            "selfcheck_skip": integer(row["selfcheck_skip"]),
            "open_selfcheck": integer(row["selfcheck_fail"]) + integer(row["selfcheck_skip"]),
            "panel_actual_bytes": integer(row["panel_all_actual_bytes"]),
            "panel_unreported_bytes": integer(row["panel_unreported_bytes"]),
        })

    v1 = next(row for row in overview if row["bundle"] == "v1")
    headline = [{
        "v1_regions": v1["regions"],
        "v1_ism_linkage": v1["ism_linkage"],
        "v1_open_selfcheck": v1["selfcheck_fail"] + v1["selfcheck_skip"],
        "topology_datasets": len(cohort_raw),
        "biological_samples": len({r["biological_id"] for r in cohort_raw}),
    }]

    k_order = {str(i): i for i in range(1, 8)} | {"8+": 8}
    coverage = [{
        "bundle": row["bundle"],
        "k_bin": row["k_bin"],
        "k_order": k_order[row["k_bin"]],
        "linked": integer(row["numerator"]),
        "denominator": integer(row["denominator"]),
        "coverage": ratio_from_pct(row["percent"]),
    } for row in coverage_raw]
    coverage.sort(key=lambda r: (r["k_order"], r["bundle"]))

    axis = []
    for row in axis_raw:
        if row["bundle"] != "v1" or row["scope"] != "topology_linked_distinct_ssnv":
            continue
        axis.append({
            "axis": row["axis"],
            "tested_n": integer(row["tested_n"]),
            "raw_sig_n": integer(row["raw_p_le_0_05_n"]),
            "raw_sig_rate": ratio_from_pct(row["raw_p_le_0_05_pct"]),
            "bh_sig_n": integer(row["bh_fdr_q_le_0_05_n"]),
            "bh_sig_rate": ratio_from_pct(row["bh_fdr_q_le_0_05_pct"]),
            "effect_field": row["source_effect_field"],
            "effect_unit": row["effect_unit"],
            "interpretation_gate": row["interpretation_gate"],
        })

    cohort = []
    cohort_chart = []
    for row in cohort_raw:
        item = {
            "sample": row["sample"],
            "biological_id": row["biological_id"],
            "technical_replicate": row["technical_replicate"].lower() == "true",
            "regions": integer(row["region_n"]),
            "distinct_ssnv": integer(row["distinct_ssnv_n"]),
            "tree_coverage": ratio_from_pct(row["tree_pct_all_regions"]),
            "unique_among_tree": ratio_from_pct(row["unique_pct_among_tree"]),
            "hash_match": row["hash_match"].lower() == "true",
            "sample_identity_match": row["sample_identity_all_match"].lower() == "true",
        }
        cohort.append(item)
        cohort_chart.extend([
            {"sample": item["sample"], "metric": "tree coverage", "value": item["tree_coverage"]},
            {"sample": item["sample"], "metric": "unique among tree", "value": item["unique_among_tree"]},
        ])

    lca_summary = audit["lca_ab_summary"]
    lca = [{
        "shared_chromosomes": lca_summary["shared_chrom_n"],
        "same_in_bam": lca_summary["same_in_bam_n"],
        "same_threads": lca_summary["same_threads_n"],
        "pre_lv": lca_summary["pre_lv_written"],
        "post_lv": lca_summary["post_lv_written"],
        "lca_resolved": lca_summary["lca_resolved"],
        "net_new_lv": lca_summary["net_new_lv_written"],
        "causal_gate": lca_summary["causal_ab_gate"],
    }]

    legacy_summary = [{
        "status": legacy["status"],
        "legacy_loci": legacy["coordinate_crosswalk"]["legacy_unique_loci"],
        "current_loci": legacy["coordinate_crosswalk"]["current_unique_loci"],
        "shared_loci": legacy["coordinate_crosswalk"]["shared_loci"],
        "coordinate_jaccard": legacy["coordinate_crosswalk"]["coordinate_jaccard"],
        "legacy_regions": legacy["legacy_funnel"][0]["n"],
        "legacy_candidates": legacy["legacy_funnel"][2]["n"],
        "legacy_displayed": legacy["legacy_funnel"][3]["n"],
        "current_axis_untested_16_available": legacy[
            "current_axis_encoding_provenance"
        ]["axis_untested_code_16_available"],
        "visual_evidence": visual["evidence_status"],
        "desktop_overflow": visual["current_desktop"]["widths"]["horizontalOverflow"],
        "mobile_overflow": visual["current_mobile"]["widths"]["horizontalOverflow"],
        "desktop_errors": len(visual["current_desktop"]["browserErrors"]),
        "mobile_errors": len(visual["current_mobile"]["browserErrors"]),
        "selfcheck_status": visual["current_mobile"]["selfcheck"]["status"],
        "methyl_qa_chrom": visual["current_desktop"]["selectedLocus"]["chrom"],
        "methyl_qa_pos": visual["current_desktop"]["selectedLocus"]["pos"],
        "methyl_qa_axis_code": visual["current_desktop"]["selectedLocus"]["axisCode"],
        "methyl_qa_has_summary": "此位點沒有 ISM 資料" not in (
            visual["current_desktop"]["methylSectionText"] or ""
        ),
    }]

    priorities = [
        {"priority": "P0", "item": "SEQC2 truth-set / HC BED / som.py 未整合", "status": "BLOCKED",
         "acceptance": "同一 callset 輸出 TP/FP/FN、precision、recall、F1，並保存 benchmark receipt"},
        {"priority": "P0", "item": "v1 ISM source root 已不存在且 receipt 無內容 hash", "status": "BLOCKED",
         "acceptance": "以可尋址 summary/run_params 重建；若無法，永久標示 legacy observation"},
        {"priority": "P0", "item": "LCA pre/post 不是乾淨 A/B", "status": "BLOCKED",
         "acceptance": "22 chromosome 的 in_bam、threads 與非 LCA identity 全相同"},
        {"priority": "P1", "item": "多樣本 sample identity fail-closed", "status": "DONE (source)",
         "acceptance": "COLO829 不借用 HCC1395 ISM/LCA/lineage；43 contract tests 通過"},
        {"priority": "P1", "item": "receipt v2 validation/provenance gates", "status": "DONE (source)",
         "acceptance": "selfcheck、sample/schema、dirty source、file hash 與 truth-set 分層記錄"},
        {"priority": "P1", "item": "typed multi-sample metric tables", "status": "PARTIAL",
         "acceptance": "axis_result 保存 p/q/effect/group_n/valid/missing_reason/family"},
        {"priority": "P1", "item": "legacy analytic manifest / deterministic case selection", "status": "BLOCKED",
         "acceptance": "移除 /tmp 上游依賴；保存 argv/commit/hash 與 18→14 inclusion/exclusion reason"},
        {"priority": "P1", "item": "多樣本 ISM/lineage drilldown", "status": "NOT STARTED",
         "acceptance": "6 biological samples 同 window/metric contract，逐樣本分母 + median/IQR"},
        {"priority": "P2", "item": "mobile sticky/overflow 與 UI 誠實度", "status": "DONE (source)",
         "acceptance": "390 px 無 body overflow/重疊、無 filterMap 不宣稱可點、browser errors=0"},
    ]

    def snapshot_source(source_id: str, dataset: str, sql: str, description: str) -> dict:
        """Canonical query selecting reviewed rows from one bounded snapshot dataset."""
        return {
            "id": source_id,
            "label": f"Reviewed snapshot: {dataset}",
            "query": {
                "engine": "portable artifact snapshot",
                "language": "sql",
                "sql": sql,
                "description": description,
                "tables_used": [dataset],
                "executed_at": generated,
            },
        }

    sources = [
        {"id": "audit_metrics", "label": "HCC1395 drilldown audit metrics",
         "path": "research/20260813_hcc1395_drilldown_validation/results/metrics_audit.json"},
        {"id": "contract_tests", "label": "Drilldown regression contracts (43 passed)",
         "path": "tests/test_drilldown_multisample_contract.py"},
        {"id": "source_fixes", "label": "Current drilldown generator source",
         "path": "scripts/build_drilldown_dashboard.py"},
        {"id": "legacy_method", "label": "Legacy methyl browser method audit",
         "path": "research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md"},
        {"id": "legacy_crosswalk", "label": "Legacy/current coordinate crosswalk",
         "path": "research/20260813_hcc1395_drilldown_validation/results/legacy_browser_crosswalk.json"},
        {"id": "visual_receipt", "label": "Direct-generated browser QA receipt",
         "path": "research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.generated.json"},
        {"id": "claude_round4", "label": "Claude Code Round 4 read-only adversarial review",
         "path": "research/20260813_hcc1395_drilldown_validation/claude_code_round4_review.md"},
        snapshot_source("q_headline", "headline", "SELECT * FROM headline;",
                        "Headline metrics selected from the reviewed bounded snapshot."),
        snapshot_source("q_overview", "overview", "SELECT * FROM overview ORDER BY bundle;",
                        "Exact v1/v3 bundle metrics."),
        snapshot_source("q_coverage", "coverage",
                        "SELECT * FROM coverage ORDER BY k_order, bundle;",
                        "ISM linkage numerator, denominator, and rate by bundle and active-k bin."),
        snapshot_source("q_axis", "axis", "SELECT * FROM axis ORDER BY axis;",
                        "Topology-linked v1 methylation-axis audit rows."),
        snapshot_source("q_cohort", "cohort", "SELECT * FROM cohort ORDER BY sample;",
                        "Seven-dataset topology inventory with per-dataset denominators."),
        snapshot_source("q_cohort_chart", "cohort_chart",
                        "SELECT * FROM cohort_chart ORDER BY sample, metric;",
                        "Long-form tree coverage and uniqueness rates for the cohort chart."),
        snapshot_source("q_priorities", "priorities",
                        "SELECT * FROM priorities ORDER BY priority, item;",
                        "Prioritized remediation register assembled from the verified audit."),
    ]

    cards = [
        {"id": "c_regions", "dataset": "headline", "sourceId": "q_headline",
         "description": "v1 的全染色體 region 母體。",
         "metrics": [{"label": "v1 regions", "field": "v1_regions", "format": "number"}]},
        {"id": "c_ism", "dataset": "headline", "sourceId": "q_headline",
         "description": "distinct topology sSNV 中有 v1 summary row 的比例。",
         "metrics": [{"label": "v1 ISM linkage", "field": "v1_ism_linkage", "format": "percent"}]},
        {"id": "c_selfcheck", "dataset": "headline", "sourceId": "q_headline",
         "description": "原 v1 bundle 的 FAIL + SKIP；不是通過數。",
         "metrics": [{"label": "open selfcheck items", "field": "v1_open_selfcheck", "format": "number"}]},
        {"id": "c_cohort", "dataset": "headline", "sourceId": "q_headline",
         "description": "7 datasets 中 HCC1395_DORADO 是 technical replicate，故 biological n=6。",
         "metrics": [
             {"label": "topology datasets", "field": "topology_datasets", "format": "number"},
             {"label": "biological samples", "field": "biological_samples", "format": "number"},
         ]},
    ]

    charts = [
        {
            "id": "ch_coverage_k", "title": "ISM linkage 隨 active k 上升而下降",
            "subtitle": "每點保留 linked / denominator；v1 與 v3 的 window 不同，不是 controlled A/B。",
            "type": "line", "intent": "trend", "dataset": "coverage", "sourceId": "q_coverage",
            "question": "甲基 summary 的可得性是否與 topology 複雜度無關？",
            "rationale": "若 linkage 隨 k 系統性變化，methylation 子集不能當成母體的隨機樣本。",
            "comparisonContext": {"grain": "bundle × active-k bin", "denominator": "distinct topology sSNV", "unit": "rate"},
            "encodings": {
                "x": {"field": "k_bin", "type": "ordinal", "label": "active k"},
                "y": {"field": "coverage", "type": "quantitative", "format": "percent", "label": "ISM linkage"},
                "color": {"field": "bundle", "type": "nominal", "label": "bundle"},
                "tooltip": [
                    {"field": "linked", "type": "quantitative", "label": "linked"},
                    {"field": "denominator", "type": "quantitative", "label": "denominator"},
                ],
            },
            "maxRows": 20,
        },
    ]

    tables = [
        {"id": "t_bundle", "title": "v1 / v3 bundle exact metrics", "dataset": "overview",
         "sourceId": "q_overview", "density": "dense", "columns": [
             {"field": "bundle", "label": "bundle", "type": "text"},
             {"field": "status", "label": "validation", "type": "text"},
             {"field": "regions", "label": "regions", "format": "number"},
             {"field": "tree_coverage", "label": "tree coverage", "format": "percent"},
             {"field": "unique_among_tree", "label": "unique among tree", "format": "percent"},
             {"field": "ism_linkage", "label": "ISM linkage", "format": "percent"},
             {"field": "open_selfcheck", "label": "FAIL + SKIP", "format": "number"},
         ]},
        {"id": "t_axis", "title": "v1 topology-linked methylation axis audit", "dataset": "axis",
         "sourceId": "q_axis", "density": "dense", "columns": [
             {"field": "axis", "label": "axis", "type": "text"},
             {"field": "tested_n", "label": "tested n", "format": "number"},
             {"field": "raw_sig_rate", "label": "raw p≤.05", "format": "percent"},
             {"field": "bh_sig_rate", "label": "recomputed BH q≤.05", "format": "percent"},
             {"field": "interpretation_gate", "label": "gate", "type": "text"},
         ]},
        {"id": "t_cohort", "title": "Seven-dataset topology inventory", "dataset": "cohort",
         "sourceId": "q_cohort", "density": "dense", "defaultSort": {"field": "sample", "direction": "asc"},
         "columns": [
             {"field": "sample", "label": "dataset", "type": "text"},
             {"field": "biological_id", "label": "biological sample", "type": "text"},
             {"field": "technical_replicate", "label": "technical replicate", "type": "text"},
             {"field": "regions", "label": "regions", "format": "number"},
             {"field": "distinct_ssnv", "label": "distinct sSNV", "format": "number"},
             {"field": "tree_coverage", "label": "tree coverage", "format": "percent"},
             {"field": "unique_among_tree", "label": "unique among tree", "format": "percent"},
         ]},
        {"id": "t_priorities", "title": "Prioritized remediation register", "dataset": "priorities",
         "sourceId": "q_priorities", "density": "spacious", "columns": [
             {"field": "priority", "label": "priority", "type": "text"},
             {"field": "item", "label": "item", "type": "text"},
             {"field": "status", "label": "current state", "type": "text"},
         ]},
    ]

    blocks = [
        {"id": "b_summary", "type": "markdown", "sourceId": "audit_metrics", "body":
         "## 技術摘要 — 兩個 bundle 都是 observation，不是 validation\n\n"
         "**v1 與 v3 均為 BLOCKED。** v1 的 11,590 regions、19,849 distinct sSNV 與 81.591% ISM linkage 可作單樣本描述；但原 v1 自檢有 1 FAIL + 1 SKIP、ISM source root 已不存在、沒有 SEQC2 truth-set benchmark，LCA pre/post 也不是乾淨 A/B。v3 多耗約 2.9 GB，ISM linked sSNV 只增加 109，仍未解除上述 claim gate。\n\n"
         "**決策：保留 v1/v3 原檔做 legacy observation；不要 cite、不要覆寫。** 後續只從已修補的產生器建立新 candidate，通過來源、身分、自檢與 truth-set 四層 gate 後才升級。"},
        {"id": "b_metrics", "type": "metric-strip", "cardIds": ["c_regions", "c_ism", "c_selfcheck", "c_cohort"]},
        {"id": "b_bundle_intro", "type": "markdown", "sourceId": "audit_metrics", "body":
         "## v3 的容量成本沒有換到相稱的驗證力\n\n"
         "v3 的主要新增是 ±5 kb ISM 與更大的 IGV/面板資產；它把 topology-linked methylation sSNV 從 16,195 提到 16,304（+109），但 window 改變，因此不是 controlled A/B。兩版 receipt 的 panel bytes 都只記 base PNG，分別漏報 82,516,832 與 144,306,817 bytes。"},
        {"id": "b_bundle_table", "type": "table", "tableId": "t_bundle"},
        {"id": "b_coverage_intro", "type": "markdown", "sourceId": "audit_metrics", "body":
         "## 甲基可得性不是隨機缺失\n\n"
         "v1 linkage 從 k=1 的 85.42% 降到 k≥8 的 41.05%；v3 只有小幅改善。這表示高複雜度 topology 較容易缺甲基 summary，任何 k × methylation 結論都必須保留 linked / denominator，不能把 available subset 當成全體。"},
        {"id": "b_coverage_chart", "type": "chart", "chartId": "ch_coverage_k"},
        {"id": "b_coverage_after", "type": "markdown", "sourceId": "audit_metrics", "body":
         "**解讀：** 這張圖支持 selection-bias 警告，不支持生物機制。要比較樣本，必須先凍結 window、metric、valid/missing 規則，再逐樣本計算 coverage；不得把未測得當作不顯著。"},
        {"id": "b_axis_intro", "type": "markdown", "sourceId": "audit_metrics", "body":
         "## p 值、FDR、effect 與 missingness 必須拆開\n\n"
         "下表的 BH q-value 是本次依保存的 raw p、在 axis × scope 內重新計算，只能作 exploratory audit；artifact 沒有記原始 multiplicity provenance。cluster 軸 11,545/11,546 topology-linked tests 顯著，因分群與檢定使用同一甲基距離，屬 double-dipping，不能作 evidence。其他軸仍欠一致的 group n、valid 與 effect 定義。"},
        {"id": "b_axis_table", "type": "table", "tableId": "t_axis"},
        {"id": "b_lca", "type": "markdown", "sourceId": "audit_metrics", "body":
         "## 4.969× 是描述值，不是 LCA 因果效果\n\n"
         f"共同 22 條染色體的 lv_written 為 {lca[0]['pre_lv']:,} → {lca[0]['post_lv']:,}，但只有 {lca[0]['same_in_bam']}/22 的 in_bam 相同，threads 為 {lca[0]['same_threads']}/22 相同。threads 本身不是生物 confounder，但其差異破壞『唯一變項只有 LCA』的 controlled-comparison 前提。即使非 LCA stats 逐染色體相等，input identity gate 仍 FAIL；341,482 net-new lv 與 341,819 lca_resolved 也不是同一量。"},
        {"id": "b_cohort_intro", "type": "markdown", "sourceId": "audit_metrics", "body":
         "## Upstream 已有七套資料，但 downstream 還不是多樣本產品\n\n"
         "七套 topology + MLHP 的 schema、model、AF basis 與 receipt hash 可比；其中 HCC1395_DORADO 是 technical replicate，所以 biological n=6。tree coverage 為 64.17%–81.71%，unique-among-tree 為 35.29%–90.58%，異質性大，不可 pooled。`drilldown_out` 仍只有 HCC1395 v1/v3，尚無六 biological samples 的 ISM/lineage bundle。"},
        {"id": "b_cohort_after", "type": "markdown", "sourceId": "audit_metrics", "body":
         "**下一個合理單位是 per-sample macro，不是 pooled loci。** 建議固定 biological replicate policy，報每樣本 numerator/denominator，再用 median + IQR；DORADO 留作 reproducibility 對照，不算第七個獨立生物樣本。"},
        {"id": "b_cohort_table", "type": "table", "tableId": "t_cohort"},
        {"id": "b_legacy", "type": "markdown", "sourceId": "legacy_crosswalk", "body":
         "## 2026-07-26 standalone：搬閱讀骨架，不搬 A/B claim\n\n"
         f"legacy/current coordinate loci 為 {legacy_summary[0]['legacy_loci']:,}/{legacy_summary[0]['current_loci']:,}，"
         f"只共享 {legacy_summary[0]['shared_loci']:,}（coordinate Jaccard {legacy_summary[0]['coordinate_jaccard']:.3%}）。"
         "legacy B 有 88.0% 亦達 allele-high 門檻；472 candidates 中 68 含 FP，且 18 個 A≥3 eligible regions 只顯示 14 個，沒有 deterministic tie-break。"
         "因此只移植 claim ceiling、population funnel、progressive disclosure 與案例圖共置；不映射 A/B taxonomy、raw linkage 或 selection prevalence。"
         "crosswalk 讀 immutable v1 axis codes 0–8；後加的 AXIS_UNTESTED=16 不會回填舊 bundle。"},
        {"id": "b_visual", "type": "markdown", "sourceId": "visual_receipt", "body":
         "## Source-complete direct-generated browser QA\n\n"
         f"machine receipt 為 `{legacy_summary[0]['visual_evidence']}`；desktop/mobile horizontal overflow 都是 false，"
         f"browser errors 分別為 {legacy_summary[0]['desktop_errors']}/{legacy_summary[0]['mobile_errors']}。"
         "頁面如實顯示 2 FAIL / 0 SKIP → BLOCKED，並提供 4 個 co-occurrence tables、106 個 denominator labels、0 個 fake cells。"
         f"甲基 detail 另機器選出 {legacy_summary[0]['methyl_qa_chrom']}:"
         f"{legacy_summary[0]['methyl_qa_pos']:,}（axis code {legacy_summary[0]['methyl_qa_axis_code']}），"
         "呈現 105 reads、371 CpG、raw p/effect 與未檢定第三態。"
         "此 staging 使用完整 v3 summary，但 `igv=skip, panels=0`；它驗 IA/runtime，不驗完整影像 delivery。"},
        {"id": "b_review", "type": "markdown", "sourceId": "claude_round4", "body":
         "## Claude Code Round 4 independent review\n\n"
         "唯讀 reviewer fresh 重跑 43 個 pytest 與五個 JavaScript syntax checks，工程/資料產品層判 `ACCEPT`；"
         "legacy taxonomy/linkage/selection 反證、sample fail-closed、receipt、LCA gate、crosswalk 與 machine receipt 均可獨立重算。"
         "科學層仍是 `BLOCKED`、多樣本仍是 `PARTIAL`；此 ACCEPT 不回溯升級 immutable v1/v3 或 legacy standalone。"},
        {"id": "b_fixes", "type": "markdown", "sourceId": "contract_tests", "body":
         "## 產生器已完成的安全修補\n\n"
         "目前 source 已加入 sample-aware routing 與 topology/MLHP/ISM/LCA/lineage/strict-edge identity checks；COLO829 不再借用 HCC1395 extension。receipt v2 分列 generator commit/dirty/hash、command、capabilities、inputs、output inventory 與 validation gates；自檢標題、假互動、硬編碼 CN claim、mobile sticky/overflow 也已修正。\n\n"
         "回歸命令 `python3 -m pytest tests/test_drilldown_contract.py tests/test_drilldown_multisample_contract.py -q` 實得 **43 passed**。這證明 contract 修補通過，不等於既有 v1/v3 bundle 被重新驗證。"},
        {"id": "b_priorities", "type": "table", "tableId": "t_priorities"},
        {"id": "b_method", "type": "markdown", "sourceId": "audit_metrics", "body":
         "## 指標定義與方法\n\n"
         "- distinct sSNV：topology `active_positions` 的 unique (chrom, 1-based position)。\n"
         "- tree coverage：有 representative tree 的 regions / all regions。\n"
         "- unique among tree：unique-best-tree regions / regions with a tree。\n"
         "- ISM linkage：有 significance summary row 的 distinct topology sSNV / distinct topology sSNV。\n"
         "- BH-FDR：本次僅在單一 bundle × axis × declared scope 內重算；跨樣本不 pool p。\n"
         "- strict-edge C12：同一 construction 的 round-trip consistency；不是 read-linkage 的獨立證明。"},
        {"id": "b_limits", "type": "markdown", "body":
         "## 限制與 robustness\n\n"
         "**本報告完成的是全量 bundle engineering audit，不是全基因組 truth validation。** v1/v3 的所有現存檔案與 22 chromosomes 都納入，但 v1 ISM 上游已遺失；KDE-corrected 狀態、原始 multiple-testing family 與 truth-set performance 無法回溯。截圖只能驗 UI，不會提升數據 claim。多樣本部分完成七套 topology inventory；多樣本 ISM/lineage 仍屬 PARTIAL scope。"},
        {"id": "b_next", "type": "markdown", "body":
         "## 建議執行順序\n\n"
         "1. 將 v1/v3 標為 `QUARANTINE / DO-NOT-CITE`，保留唯讀 legacy observation。\n"
         "2. 以乾淨 commit 重建一個 HCC1395 candidate，要求 receipt validation gates 全過。\n"
         "3. 整合 SEQC2 HC BED + som.py，新增 TP/FP/FN/precision/recall/F1 與 benchmark receipt。\n"
         "4. 重做同 input、同 threads 的 LCA pre/post controlled A/B。\n"
         "5. 凍結 multi-sample contract 後擴到 6 biological samples；先輸出逐樣本表，再做 median/IQR。"},
        {"id": "b_questions", "type": "markdown", "body":
         "## 仍會改變決策的問題\n\n"
         "- 要把哪一個 caller/callset 定為 SEQC2 truth benchmark 的主輸入？\n"
         "- v1 缺失的 ISM source 是否有可驗 hash 的 archive，還是應永久降為 legacy？\n"
         "- 多樣本 phase 要以 ±1 kb 還是 ±5 kb 為 frozen window？兩者不可直接混合。\n"
         "- HCC1395_DORADO 要只做 technical reproducibility，還是另設 platform-stratified analysis？"},
    ]

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "HCC1395 drilldown 驗證稽核",
            "description": "完整 bundle 工程稽核、legacy methyl browser 方法/crosswalk、direct-generated 視覺 QA、Claude Code 交叉審查與多樣本擴展評估。",
            "generatedAt": generated,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated,
            "status": "blocked",
            "datasets": {
                "headline": headline,
                "overview": overview,
                "coverage": coverage,
                "axis": axis,
                "cohort": cohort,
                "cohort_chart": cohort_chart,
                "lca": lca,
                "legacy_summary": legacy_summary,
                "priorities": priorities,
            },
            "accessIssues": [
                {"id": "truth_missing", "scope": "scientific validity", "sourceId": "audit_metrics",
                 "message": "SEQC2 truth calls / HC BED / som.py benchmark are not present in either bundle."},
                {"id": "v1_source_missing", "scope": "reproducibility", "sourceId": "audit_metrics",
                 "message": "The v1 ISM source root no longer exists and its legacy receipt did not content-address the directory."},
                {"id": "multisample_downstream_missing", "scope": "G4 multi-sample validation", "sourceId": "audit_metrics",
                 "message": "Seven topology datasets exist, but comparable multi-sample ISM/lineage drilldown bundles do not."},
                {"id": "legacy_manifest_missing", "scope": "legacy reproducibility", "sourceId": "legacy_method",
                 "message": "The legacy builder uses /tmp intermediates and does not encode a deterministic 18-to-14 case-selection contract."},
                {"id": "full_assets_not_evaluated", "scope": "asset delivery", "sourceId": "visual_receipt",
                 "message": "Direct-generated QA used igv=skip and panels=0; full image-asset delivery remains NOT EVALUATED."},
            ],
        },
        "sources": sources,
    }
    OUTPUT.write_text(json.dumps(artifact, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "blocks": len(blocks), "cards": len(cards), "charts": len(charts),
        "tables": len(tables), "datasets": {k: len(v) for k, v in artifact["snapshot"]["datasets"].items()},
        "status": artifact["snapshot"]["status"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
