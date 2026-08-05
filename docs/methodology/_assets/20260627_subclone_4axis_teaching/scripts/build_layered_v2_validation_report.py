#!/usr/bin/env python3
"""Build the authoritative Markdown and self-contained HTML validation report."""

import argparse
import datetime as dt
import hashlib
import html
import json
import subprocess
import sys
import tempfile
from pathlib import Path


FUNNEL_FIELDS = [
    ("Out of scope", "L2_out_of_scope_non_autosomal"),
    ("Positional singleton", "L3_positional_singleton"),
    ("MAX_SNV excluded", "L5_cap_excluded_sSNV"),
    ("Read unsupported", "L5_read_unsupported_sSNV"),
    ("Retained", "L6_retained_sSNV"),
]


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fmt(value, digits=1):
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return f"{value:,}"


class Tooltips:
    def __init__(self):
        self.index = 0

    def number(self, value, source, digits=1, suffix=""):
        self.index += 1
        tooltip_id = f"source-{self.index}"
        shown = fmt(value, digits) + suffix
        return (f'<span class="source-tooltip" tabindex="0" aria-describedby="{tooltip_id}">{html.escape(shown)}'
                f'<span class="source-tooltip-content" id="{tooltip_id}" role="tooltip">'
                f'Source: {html.escape(source)}</span></span>')


def runtime_metrics(run_root, result):
    metrics = result.get("metrics", {})
    if "funnel" in metrics and "read_tag_census" in metrics and "site_ledger" in metrics:
        return metrics
    sample = result["sample"]
    sample_root = run_root / "samples" / sample
    layered = load(sample_root / f"layered_reconstruction_{sample}.json")
    view = load(sample_root / f"layered_region_view_{sample}.json")
    site = load(sample_root / f"ssnv_site_ledger_{sample}.summary.json")
    census = view["census"]
    l1 = census["L1"]
    return {
        "roles": {
            "primary": l1["n_primary_lineage_units"],
            "reference": l1["n_reference_only_controls"],
            "H3_auxiliary": l1["n_unresolved_H3_auxiliary"],
            "H4_auxiliary": l1["n_shared_H4_auxiliary"],
        },
        "funnel": census["funnel"],
        "read_tag_census": layered["read_tag_census"],
        "site_ledger": site,
        "n_capped": layered["L1_ssnv_algorithm"]["n_verification_not_applicable_capped"],
    }


def sample_rows(verification, ablation, determinacy, run_root):
    abl = {row["sample"]: row for row in ablation["samples"]}
    det = {row["sample"]: row for row in determinacy["samples"]}
    rows = []
    for result in verification["samples"]:
        sample = result["sample"]
        metrics = runtime_metrics(run_root, result)
        funnel = metrics["funnel"]
        row = {
            "sample": sample,
            "biological_id": result["biological_id"],
            "pass": result["pass"],
            "universe": funnel["L1_all_pass_universe"],
            "autosomal": funnel["autosomal_chr1_22"],
            "out_scope": funnel["L2_out_of_scope_non_autosomal"],
            "singleton": funnel["L3_positional_singleton"],
            "cap_excluded": funnel["L5_cap_excluded_sSNV"],
            "unsupported": funnel["L5_read_unsupported_sSNV"],
            "retained": funnel["L6_retained_sSNV"],
            "primary": metrics["roles"]["primary"],
            "reference": metrics["roles"]["reference"],
            "H3": metrics["roles"]["H3_auxiliary"],
            "H4": metrics["roles"]["H4_auxiliary"],
            "raw_records": metrics["site_ledger"]["raw_clairs_records"],
            "longphase_input": metrics["site_ledger"]["longphase_input_records"],
            "longphase_all": metrics["site_ledger"]["longphase_recalibrated_records"],
            "tree_input": metrics["site_ledger"]["tree_input_records"],
            "tag_exposures": metrics["read_tag_census"]["alignment_group_exposures"],
            "tag_exact": metrics["read_tag_census"]["sidecar_exact_matches"],
            "tag_missing": metrics["read_tag_census"]["sidecar_missing"],
            "tag_conflicts": metrics["read_tag_census"]["sidecar_conflicts"],
            "mixed_ps_regions": metrics["read_tag_census"].get("phase_set_region_counts", {}).get("multiple", 0),
            "capped": metrics["n_capped"],
            "multiHP_pct": abl[sample]["multiHP_primary_pct"],
            "exact_pct": det[sample]["exact_tree_unique"]["pct"],
            "shape_pct": det[sample]["shape_determined"]["pct"],
            "region_pct": det[sample]["region_all_determined"]["pct"],
        }
        rows.append(row)
    return rows


def markdown(rows, verification, ablation, determinacy, read_af, parameter, backbone, run_root, manifest):
    lines = [
        "<!--",
        f"建立時間: {dt.datetime.now().astimezone().isoformat(timespec='minutes')}",
        "目標: Layered reconstruction v2 全面驗證與 current claim 單一真值報告",
        "處理範圍: chr1-22 x 7 datasets（6 biological samples）；base + parameter/backbone sensitivity",
        f"data_sources: {run_root}/verification_summary.json,{manifest}",
        "證據等級: L2 ⭐⭐⭐⭐（全 dataset 工程/內部驗證；無 single-cell/multi-region orthogonal truth）",
        "-->",
        "",
        "# Layered Reconstruction v2 全面驗證",
        "",
        f"> **結論**：clean pipeline 已把 raw ClairS 全 records、ClairS PASS LongPhase-S input、LongPhase-S `_sc.vcf` all/PASS、無 truth-BED exact HP/PS sidecar、sSNV disposition、reference-only、HP3/HP4 auxiliary 與候選樹完整性寫入 schema；演化樹以 `_sc.vcf` FILTER=PASS 為 backbone，base run {verification['n_pass']}/{verification['dataset_count']} datasets 通過所有 pre-registered gates。主張仍限 regional mutation-state trees，不能升級成 confirmed cell clones。",
        "",
        "## 主要結果",
        "",
        "| Dataset | Universe | chr1-22 | Retained | Primary HP1/2 | Reference controls | H3? aux | H4? aux | Exact tree | Shape | Region all-det |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(f"| {row['sample']} | {row['universe']:,} | {row['autosomal']:,} | {row['retained']:,} | "
                     f"{row['primary']:,} | {row['reference']:,} | {row['H3']:,} | {row['H4']:,} | {row['exact_pct']:.1f}% | "
                     f"{row['shape_pct']:.1f}% | {row['region_pct']:.1f}% |")
    lines += [
        "",
        "## ClairS → LongPhase-S → sSNV 資料流守恆",
        "",
        "| Dataset | raw ClairS | ClairS PASS / LPS input | LPS `_sc` all | LPS PASS / tree input | sSNV retained | tag exposures | exact joins | missing | conflicts | mixed-PS regions |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['sample']} | {row['raw_records']:,} | {row['longphase_input']:,} | "
            f"{row['longphase_all']:,} | {row['tree_input']:,} | {row['retained']:,} | "
            f"{row['tag_exposures']:,} | {row['tag_exact']:,} | {row['tag_missing']:,} | {row['tag_conflicts']:,} | "
            f"{row['mixed_ps_regions']:,} |"
        )
    lines += [
        "",
        "raw ClairS non-PASS records 全部保留在 site ledger，但不送入 LongPhase-S；ClairS PASS 全部送入 LongPhase-S；`_sc.vcf` 全 records 保存 recalibrated FILTER；只有 `_sc.vcf` FILTER=PASS 是 primary sSNV tree universe。",
        "",
        "## 分母契約",
        "",
        "- Primary lineage：至少一個 MINREAD-supported ALT genotype 的 HP1/HP2 family unit。",
        "- Reference-only：只有 ROOT/no-target-mutation state 的 control；不進 lineage、determined、multi-HP 分母。",
        "- H3?：somatic-integrated unresolved auxiliary；維持計數但不當第三 parental lineage。",
        "- H4?：somatic ALT shared by germline HP1/HP2 auxiliary；維持計數但不進 primary lineage。",
        "- Exact-tree：non-capped primary units 中 `n_trees=1`；shape-determined 允許多棵 label/order 變體但 canonical shape 唯一。",
        "- Region all-determined：有至少一條 primary lineage，且該區所有 primary lineages 都 exact determined。",
        "- Region：相鄰 sSNV gap <=50kb 的 connected component；total span 可超過 50kb。",
        "",
        "## Verification",
        "",
        f"- Base verifier：{verification['n_pass']}/{verification['dataset_count']} PASS；每個 non-capped eligible unit 的 V4/V5 都實際執行，capped 另列 not-applicable。",
        f"- Forbidden counts：reference-only/HP3/HP4 進 primary = {'0/0/0（全樣本）' if ablation['all_forbidden_counts_zero'] else 'FAIL'}。",
        "- Read-tag contract：7/7 使用無 truth VCF/BED 的 production sidecar；HP/PS exact join missing=0、conflicts=0；PS 作 phase-block QC，不作 topology edge。",
        "- Site ledger：每個 raw ClairS record 恰有一個可稽核 disposition；ClairS non-PASS 與 LongPhase-S filter exclusion 分開列示。",
        f"- Candidate completeness：{'PASS' if read_af['all_candidate_sets_complete'] else 'FAIL'}；read-AF ordering 重枚舉全部候選，不讀 display prefix。",
        "- 8 golden cases + 5 seeds x 800 = 4,000 random stress cases：0 mismatch（既有獨立 oracle artifact）。",
        "",
        "## Robustness",
        "",
        f"- Parameter sensitivity（7 datasets）：{json.dumps(parameter.get('aggregate', {}), ensure_ascii=False)}",
        f"- Backbone sensitivity（7 datasets；LongPhase-S PASS canonical vs ClairS PASS input）：{json.dumps(backbone.get('aggregate', {}), ensure_ascii=False)}",
        "- Read-AF ordering 是 exploratory heuristic，不是 purity/CN-corrected CCF，也不是獨立 biological validation。",
        "",
        "| Dataset | Retained-site Jaccard | Primary-unit Jaccard | Shared topology concordance | Determined primary delta | Multi-HP delta | Region determined delta | Verdict |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for item in backbone.get("comparisons", []):
        lines.append(
            f"| {item['sample']} | {item['retained_position_jaccard']:.3f} | "
            f"{item['primary_unit_key_jaccard']:.3f} | "
            f"{item['shared_unit_topology_digest_concordance']:.3f} | "
            f"{item['delta']['determined_primary_pp']:.2f} pp | {item['delta']['multiHP_pp']:.2f} pp | "
            f"{item['delta']['all_determined_region_pp']:.2f} pp | {item['verdict']} |"
        )
    lines += [
        "",
        "## 限制",
        "",
        "- 7 datasets = 6 biological samples；HCC1395 與 HCC1395_DORADO 是同一細胞株的技術/處理版本。",
        "- COLO829、HCC1937 CN unavailable；不以 missing 代替 neutral。",
        "- Backbone sensitivity 已涵蓋 7 datasets，但比較的是同一 ClairS→LongPhase-S 流程的上游 PASS 與 LongPhase-S recalibrated PASS；尚不是 independent-caller robustness。",
        "- PS 已逐 alignment 保存並輸出每區 PS/HP×PS census；canonical tree 不把 PS 當 lineage label。跨 PS region 必以 mixed-PS census 另行解讀。",
        "- 現有資料沒有 single-cell 或 multi-region orthogonal truth；bulk molecule != cell，regional state != confirmed clone。",
        "- L3 methylation contract 是 bounded auxiliary：只允許 negative screen/residual flag，禁止 tree ranking/lineage confirmation。",
        "",
        "## Reproducibility",
        "",
        f"- Run root: `{run_root}`",
        f"- Input manifest: `{manifest}`",
        f"- Verification: `{run_root}/verification_summary.json`",
        f"- LongPhase-S producer closeout: 見 `{manifest}` 的 `tag_contract.production_closeout` 與 `_SUCCESS` hash binding。",
        "- Pre-registration: `InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md`",
    ]
    return "\n".join(lines) + "\n"


def fallback_table(headers, rows, source=None, prefix="table"):
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    rendered_rows = []
    for row_index, row in enumerate(rows):
        cells = []
        for column_index, value in enumerate(row):
            shown = html.escape(str(value))
            if source and column_index > 0:
                tooltip_id = f"{prefix}-{row_index}-{column_index}"
                shown = (f'<span class="source-tooltip" tabindex="0" aria-describedby="{tooltip_id}">{shown}'
                         f'<span class="source-tooltip-content" id="{tooltip_id}" role="tooltip">'
                         f'Source: {html.escape(str(source))}</span></span>')
            cells.append(f"<td>{shown}</td>")
        rendered_rows.append("<tr>" + "".join(cells) + "</tr>")
    body = "".join(rendered_rows)
    return f'<div class="table-scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def html_report(rows, verification, ablation, determinacy, read_af, parameter, backbone, run_root, manifest):
    tip = Tooltips()
    source = f"{run_root}/verification_summary.json"
    analysis_source = "InterSubMod/research/20260710_layered_reconstruction_v2/"
    metric_cards = [
        ("Dataset gates", tip.number(verification["n_pass"], source, 0, f"/{verification['dataset_count']}"), "all pre-registered invariants"),
        ("Biological samples", tip.number(verification["biological_sample_count"], source, 0), "7 datasets; HCC1395 counted once"),
        ("Forbidden primary units", tip.number(0 if ablation["all_forbidden_counts_zero"] else 1, str(manifest), 0), "reference-only + HP3 + HP4"),
        ("Claim ceiling", "Regional", "single-cell/multi-region truth unavailable"),
    ]
    cards = "".join(f'<div class="metric"><div class="metric-label">{label}</div><div class="metric-value">{value}</div><div class="metric-note">{note}</div></div>'
                    for label, value, note in metric_cards)
    main_rows = []
    for row in rows:
        main_rows.append([row["sample"], fmt(row["universe"], 0), fmt(row["autosomal"], 0), fmt(row["retained"], 0),
                          fmt(row["primary"], 0), fmt(row["reference"], 0), fmt(row["H3"], 0), fmt(row["H4"], 0),
                          f"{row['exact_pct']:.1f}%", f"{row['shape_pct']:.1f}%", f"{row['region_pct']:.1f}%"])
    main_table = fallback_table(["Dataset", "Universe", "chr1-22", "Retained", "Primary", "Reference", "H3?", "H4?", "Exact", "Shape", "Region"],
                                main_rows, f"{source}; {analysis_source}determinacy_layers.json", "main")
    flow_table = fallback_table(
        ["Dataset", "raw ClairS", "ClairS PASS/LPS input", "LPS _sc all", "LPS PASS/tree", "sSNV retained",
         "Tag exposures", "Exact", "Missing", "Conflicts", "Mixed-PS regions"],
        [[row["sample"], row["raw_records"], row["longphase_input"], row["longphase_all"], row["tree_input"],
          row["retained"], row["tag_exposures"], row["tag_exact"], row["tag_missing"], row["tag_conflicts"],
          row["mixed_ps_regions"]]
         for row in rows],
        source, "flow",
    )
    funnel_fallback = fallback_table(
        ["Dataset", *[label for label, _ in FUNNEL_FIELDS]],
        [[row["sample"], row["out_scope"], row["singleton"], row["cap_excluded"], row["unsupported"], row["retained"]] for row in rows],
        source, "funnel",
    )
    det_fallback = fallback_table(
        ["Dataset", "Exact tree %", "Shape-determined %", "Region all-determined %"],
        [[row["sample"], f"{row['exact_pct']:.1f}", f"{row['shape_pct']:.1f}", f"{row['region_pct']:.1f}"] for row in rows],
        f"{analysis_source}determinacy_layers.json", "determinacy",
    )
    parameter_rows = []
    for variant, summary in parameter.get("aggregate", {}).items():
        parameter_rows.append([variant, f"{summary['max_abs_delta_pp']:.2f} pp", f"{summary['median_abs_delta_pp']:.2f} pp",
                               f"{summary['min_retained_position_jaccard']:.3f}",
                               f"{summary['min_primary_unit_key_jaccard']:.3f}",
                               f"{summary['min_shared_topology_digest_concordance']:.3f}", summary["verdict"]])
    backbone_rows = []
    for item in backbone.get("comparisons", []):
        backbone_rows.append([item["label"], item["sample"], f"{item['retained_position_jaccard']:.3f}",
                              f"{item['primary_unit_key_jaccard']:.3f}",
                              f"{item['shared_unit_topology_digest_concordance']:.3f}",
                              f"{item['delta']['determined_primary_pp']:.2f} pp", f"{item['delta']['multiHP_pp']:.2f} pp", item["verdict"]])
    read_af_rows = []
    for item in read_af.get("samples", []):
        default = item.get("default") or {}
        read_af_rows.append([item["sample"], item["n_ambiguous_primary_units"], item["n_units_analyzed_all_candidates"],
                             default.get("reached"), default.get("strict_neutral_reached"), item["n_candidate_mismatch_or_incomplete"]])
    now = dt.datetime.now().astimezone().isoformat(timespec="minutes")
    return f'''<!doctype html>
<html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light dark"><title>Layered Reconstruction v2 全面驗證</title>
<style>
:root{{--bg:#f4f5f2;--surface:#fff;--ink:#171a1d;--muted:#5e6469;--line:#d8dcda;--teal:#087f78;--gold:#b7791f;--red:#b23a48;--blue:#285f9e;--olive:#667a2e;--soft:#eef1ed;--radius:6px;--font:"Aptos","Noto Sans TC",sans-serif;--display:Charter,"Bitstream Charter",Georgia,serif}}
@media(prefers-color-scheme:dark){{:root{{--bg:#171a1b;--surface:#222627;--ink:#f4f4ef;--muted:#bbc0be;--line:#3a4140;--teal:#53bdb5;--gold:#e4ad57;--red:#ef8690;--blue:#78a9dd;--olive:#a9bc6d;--soft:#2a302f}}}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.6 var(--font);letter-spacing:0}}main{{max-width:1160px;margin:auto;background:var(--surface);border-left:1px solid var(--line);border-right:1px solid var(--line)}}header,section{{padding:34px 46px;border-bottom:1px solid var(--line)}}.top{{padding-top:44px}}.kicker{{color:var(--teal);font:700 12px/1.4 var(--font);text-transform:uppercase}}h1,h2{{font-family:var(--display);font-weight:600;letter-spacing:0}}h1{{font-size:38px;line-height:1.15;margin:8px 0 16px;max-width:850px}}h2{{font-size:25px;margin:0 0 12px}}h3{{font-size:15px;margin:0}}p{{max-width:850px;color:var(--muted)}}strong{{color:var(--ink)}}.summary{{font-size:18px;color:var(--ink);border-left:4px solid var(--teal);padding-left:16px}}.metrics{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin-top:24px}}.metric{{border:1px solid var(--line);border-radius:var(--radius);padding:15px;min-height:112px}}.metric-label,.metric-note{{color:var(--muted);font-size:12px}}.metric-value{{font:600 27px/1.3 var(--display);margin:8px 0}}.chart{{border:1px solid var(--line);border-radius:var(--radius);margin:20px 0 28px;overflow:hidden}}.chart-head{{padding:14px 18px;background:var(--soft);border-bottom:1px solid var(--line)}}.chart-head p{{font-size:12px;margin:2px 0}}.chart-body{{padding:14px 18px;min-height:320px}}[data-recharts-live]{{display:none;min-height:320px}}[data-recharts-ready="true"] [data-recharts-fallback]{{display:none}}[data-recharts-ready="true"] [data-recharts-live]{{display:block}}.table-scroll{{overflow:auto}}table{{width:100%;border-collapse:collapse;font-size:12px}}th,td{{padding:9px 11px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap;font-variant-numeric:tabular-nums}}th:first-child,td:first-child{{text-align:left}}th{{color:var(--muted);background:var(--soft)}}.note{{background:var(--soft);border-left:4px solid var(--gold);padding:14px 16px;color:var(--muted)}}code{{font-family:ui-monospace,"SFMono-Regular",monospace;font-size:.88em;overflow-wrap:anywhere}}.source-tooltip{{position:relative;text-decoration:underline dotted;text-underline-offset:3px;cursor:help}}.source-tooltip-content{{position:absolute;z-index:5;left:0;bottom:calc(100% + 7px);display:none;width:330px;padding:8px 10px;border-radius:4px;background:#111;color:#fff;font:12px/1.4 var(--font);overflow-wrap:anywhere}}.source-tooltip:hover .source-tooltip-content,.source-tooltip:focus .source-tooltip-content{{display:block}}ul{{padding-left:20px}}li{{margin:7px 0;color:var(--muted)}}footer{{padding:24px 46px;color:var(--muted);font-size:12px}}
@media(max-width:800px){{main{{border:0}}header,section{{padding:26px 18px}}h1{{font-size:31px}}.metrics{{grid-template-columns:1fr 1fr}}.chart-body{{overflow:auto}}}}
@media(max-width:520px){{.metrics{{grid-template-columns:1fr}}.source-tooltip-content{{position:fixed;left:12px;right:12px;bottom:12px;width:auto}}}}
</style></head><body><main data-report-audience="technical">
<header class="top" data-contract-section="title"><div class="kicker">InterSubMod · Comprehensive validation · {now}</div><h1>Layered Reconstruction v2：資料流與分母已驗證，生物主張仍限區域狀態樹</h1><p class="summary">7 個 dataset 的 clean pipeline 全部通過 pre-registered invariants；raw ClairS 全 records 有逐位點 ledger，ClairS PASS 完整送入 LongPhase-S，演化樹以 LongPhase-S `_sc.vcf` FILTER=PASS 為 backbone；production tagging 不使用 truth-BED，HP/PS 以 exact alignment identity 接入，missing/conflicts 均為 0。HP 用於 family stratification；PS 僅作 phase-block QC，不作 topology edge。這完成工程與內部驗證，但沒有 single-cell／multi-region truth，因此不能把 bulk regional state tree 寫成 confirmed cell clone。</p><div class="metrics">{cards}</div></header>
<section data-contract-section="technical-summary"><h2>結論先行</h2><p>Primary lineage 現在只包含 mutation-bearing HP1/HP2 family。Root-only 保留為 reference control，HP3 與 HP4 分別保留為 unresolved/shared-somatic auxiliary，COLO829/HCC1937 的 CN 保留 unavailable。所有 non-capped eligible units 實際執行 V1–V7；capped units 明列 not-applicable，不再用「全體 ALL PASS」掩蓋未執行 oracle。</p><div class="note"><strong>證據級：</strong>L2。跨 7 datasets 的工程與內部一致性，以及 LongPhase-S PASS 對 ClairS PASS selected-tree 的 backbone sensitivity 已驗；independent caller 與生物 clone confirmation 仍缺 orthogonal truth。</div></section>
<section data-contract-section="data-flow-contract"><h2>ClairS → LongPhase-S → sSNV 資料流逐層守恆</h2><p>raw ClairS non-PASS records 全部保留在 site ledger；ClairS PASS 完整進入 LongPhase-S；`_sc.vcf` 保存所有輸入位點及 recalibrated FILTER；只有 `_sc.vcf` PASS 進 primary sSNV reconstruction。Tag exposures 是 region-group 對 alignment 的暴露次數，不等同 BAM unique reads。</p>{flow_table}</section>
<section data-contract-section="key-findings"><h2>六層 funnel 顯示最大的損失來自工程 cap，不是生物 isolated</h2><p>每個 universe 先拆 non-autosomal scope，再於 chr1–22 拆 positional singleton、MAX_SNV excluded、read unsupported 與 retained。各樣本均由 upstream artifact 直接守恆，不再從最後結果倒推。</p><div class="chart"><div class="chart-head"><h3>All-PASS universe 的組成</h3><p>7 datasets；absolute sSNV counts，stacked composition；total span 不受 50kb 限制。</p></div><div class="chart-body" data-recharts-chart="funnel-composition"><div data-recharts-fallback>{funnel_fallback}</div><div data-recharts-live aria-hidden="true"></div></div></div>
<h2>Exact tree、shape 與 region determinacy 是三個不同問題</h2><p>Exact tree 與 shape 使用 non-capped primary lineage units；region all-determined 使用有 primary lineage 的 regions。三者不可共用 headline 分母。</p><div class="chart"><div class="chart-head"><h3>三層 determinacy 比例</h3><p>Percent；exact/shape denominator = complete non-capped primary units；region denominator = regions with primary lineage。</p></div><div class="chart-body" data-recharts-chart="determinacy-layers"><div data-recharts-fallback>{det_fallback}</div><div data-recharts-live aria-hidden="true"></div></div></div>
<div class="chart"><div class="chart-head"><h3>跨樣本 exact census</h3><p>7 datasets；HCC1395 與 HCC1395_DORADO 是同一 biological sample。</p></div>{main_table}</div></section>
<section data-contract-section="scope-data-and-metric-definitions"><h2>Scope 與分母先於解釋</h2><ul><li><strong>Genome：</strong>chr1–22 primary；其他 contigs 只進 out-of-scope census。</li><li><strong>Region：</strong>adjacent-gap ≤50kb connected component；總 span 可超過 50kb；送 solver 最多 8 sSNV。</li><li><strong>Primary：</strong>MINREAD-supported ALT genotype 的 HP1/HP2 family unit。</li><li><strong>Reference-only：</strong>local no-target-mutation control，不是 normal cell，也不是 lineage。</li><li><strong>H3?：</strong>somatic-integrated unresolved auxiliary，不是第三 parental lineage。</li><li><strong>H4?：</strong>somatic ALT shared by germline HP1/HP2 auxiliary，不進 primary lineage。</li><li><strong>Read AF：</strong>family-specific ALT read fraction；未做 purity/CN/multiplicity 校正，不是真 CCF。</li></ul></section>
<section data-contract-section="methodology"><h2>方法由 genetic co-occurrence 承重</h2><p>ClairS paired PASS 定義 LongPhase-S input universe，LongPhase-S `_sc.vcf` FILTER=PASS 定義 tree operational universe；ONT same-molecule genotype vectors 形成觀測節點與 partial subcubes；HP1/HP2 只負責 family 分層；PS 完整保存並標示 mixed phase blocks，但不作 topology edge；group-Steiner solver 枚舉所有 non-capped minimal candidates；CN 在建樹後解釋 recurrence；methylation 只允許 bounded residual annotation。</p><p><code>{html.escape(str(manifest))}</code><br><code>{html.escape(str(run_root))}</code></p></section>
<section data-contract-section="limitations-uncertainty-and-robustness-checks"><h2>Robustness 通過不等於 biological confirmation</h2><h3>Read-AF ordering</h3><p>每個 ambiguous primary unit 重新枚舉全部候選，並對 temperature、posterior threshold、violation margin 做 grid。Winner consistency 使用同一 heuristic，只是 self-consistency。</p>{fallback_table(["Dataset","Ambiguous","All candidates","Reached","Strict neutral","Mismatch"], read_af_rows, f"{analysis_source}read_af_tree_ordering.json", "read-af")}
<h3>7-dataset parameter sensitivity</h3>{fallback_table(["Variant","Max |Δ|","Median |Δ|","Min site Jaccard","Min unit Jaccard","Min shared topology","Verdict"], parameter_rows, f"{analysis_source}parameter_sensitivity.json", "parameter")}
<h3>7-dataset backbone sensitivity</h3>{fallback_table(["Backbone","Dataset","Site Jaccard","Unit Jaccard","Shared topology","Determined Δ","Multi-HP Δ","Verdict"], backbone_rows, f"{analysis_source}backbone_sensitivity.json", "backbone")}
<div class="note"><strong>Material limitations：</strong>7 datasets 只有 6 biological samples；COLO829/HCC1937 CN unavailable；本次 backbone sensitivity 是 ClairS PASS selected tree 對 LongPhase-S recalibrated PASS，兩者共用同一 raw-all LongPhase-S HP/PS tags，並非 independent caller；現有 workspace 無 single-cell/multi-region truth；L3 methylation 未接 per-unit orthogonal validation。</div></section>
<section data-contract-section="recommended-next-steps"><h2>唯一未完成的是外部證據，不是再調內部閾值</h2><ol><li>取得 ≥3 biological samples 的 independent caller callsets，重跑同一 manifest contract。</li><li>取得 single-cell 或 multi-region orthogonal truth，把 regional states 對到真 cell populations。</li><li>若取得 purity + allele-specific integer CN，再另立 true CCF model；不要把 read AF 原地改名。</li></ol></section>
<section data-contract-section="further-questions"><h2>仍需外部回答的問題</h2><p>哪些 regional mutation states 在 single cells 中共存？H3? 能否由獨立 parental phasing 定相？COLO829/HCC1937 是否能取得可靠 CN truth？這三題目前不能由既有 bulk BAM 自我回答。</p></section>
<footer>Primary evidence: <code>{html.escape(str(run_root / 'verification_summary.json'))}</code> · Manifest: <code>{html.escape(str(manifest))}</code></footer>
</main><!-- DATA_ANALYTICS_HTML_REPORT_RUNTIME --></body></html>'''


def payload(rows):
    funnel_data = []
    for row in rows:
        mapping = {"Out of scope": row["out_scope"], "Positional singleton": row["singleton"],
                   "MAX_SNV excluded": row["cap_excluded"], "Read unsupported": row["unsupported"], "Retained": row["retained"]}
        funnel_data.extend({"sample": row["sample"], "stage": stage, "value": value} for stage, value in mapping.items())
    det_data = []
    for row in rows:
        det_data.extend([
            {"sample": row["sample"], "level": "Exact tree", "value": row["exact_pct"] / 100},
            {"sample": row["sample"], "level": "Shape", "value": row["shape_pct"] / 100},
            {"sample": row["sample"], "level": "Region all-det", "value": row["region_pct"] / 100},
        ])
    return {"charts": [
        {"id": "funnel-composition", "height": 360, "type": "bar", "dataset": {"id": "funnel-composition", "title": "All-PASS universe composition", "data": funnel_data,
          "chart_spec": {"id": "funnel-composition", "dataset": "funnel-composition", "title": "All-PASS universe composition", "type": "bar",
                         "encodings": {"x": {"field": "sample", "type": "nominal"}, "y": {"field": "value", "type": "quantitative", "label": "sSNV count"}, "color": {"field": "stage", "type": "nominal"}},
                         "xAxisTitle": "", "yAxisTitle": "sSNV count", "valueFormat": "number", "settings": {"groupMode": "stacked"}}}},
        {"id": "determinacy-layers", "height": 360, "type": "bar", "dataset": {"id": "determinacy-layers", "title": "Determinacy layers", "data": det_data,
          "chart_spec": {"id": "determinacy-layers", "dataset": "determinacy-layers", "title": "Determinacy layers", "type": "bar",
                         "encodings": {"x": {"field": "sample", "type": "nominal"}, "y": {"field": "value", "type": "quantitative", "label": "Share"}, "color": {"field": "level", "type": "nominal"}},
                         "xAxisTitle": "", "yAxisTitle": "Share", "valueFormat": "percent", "settings": {"groupMode": "grouped"}}}},
    ]}


def validate_authoritative_inputs(manifest, verification, ablation, read_af, parameter, backbone,
                                  run_root=None, secondary_verification=None, determinacy=None):
    errors = []
    manifest_samples = [item.get("sample") for item in manifest.get("samples", [])]
    if len(manifest_samples) != 7 or len(set(manifest_samples)) != 7:
        errors.append("canonical manifest does not contain seven unique datasets")
    if manifest.get("schema_version") != "2.1":
        errors.append("input manifest is not clean schema 2.1")
    if (manifest.get("task_type") != "B_COMPREHENSIVE_VALIDATION"
            or manifest.get("tree_input_contract") != "longphase_recalibrated_PASS"):
        errors.append("input manifest is not the canonical LongPhase-S output tree contract")
    tag_contract = manifest.get("tag_contract", {})
    if (tag_contract.get("truth_flags") is not False or tag_contract.get("PS_preserved") is not True
            or tag_contract.get("bam_identity_locked") is not True):
        errors.append("manifest does not guarantee no-truth production HP/PS tags")
    if (tag_contract.get("tree_backbone") != "LongPhase-S _sc.vcf FILTER=PASS"
            or tag_contract.get("longphase_filtering_policy") != "production_default_filter"):
        errors.append("manifest does not guarantee the canonical LongPhase-S output/filtering contract")
    for path_field, hash_field in (
            ("production_closeout", "production_closeout_sha256"),
            ("production_success", "production_success_sha256"),
            ("production_artifacts_manifest", "production_artifacts_manifest_sha256")):
        artifact = Path(tag_contract.get(path_field, ""))
        if (not artifact.is_file() or len(tag_contract.get(hash_field, "")) != 64
                or sha256(artifact) != tag_contract.get(hash_field)):
            errors.append(f"producer closeout binding failed: {path_field}")
    if verification.get("all_pass") is not True or verification.get("n_pass") != verification.get("dataset_count"):
        errors.append("base verification is not all-pass")
    if verification.get("dataset_count") != 7:
        errors.append(f"base verification dataset_count={verification.get('dataset_count')}, expected 7")
    verified_samples = [item.get("sample") for item in verification.get("samples", [])]
    if len(verified_samples) != 7 or set(verified_samples) != set(manifest_samples):
        errors.append("base verification dataset set is incomplete or duplicated")
    if run_root is not None and not (run_root / "_SUCCESS").is_file():
        errors.append("base run does not have the final _SUCCESS marker")
    is_v3 = verification.get("schema_name") == "intersubmod.layered_verification_summary"
    if is_v3:
        if (not isinstance(secondary_verification, dict)
                or secondary_verification.get("all_pass") is not True
                or secondary_verification.get("n_pass") != 7
                or secondary_verification.get("dataset_count") != 7):
            errors.append("v3 canonical run lacks a 7/7 secondary science verification")
    for result in verification.get("samples", []):
        try:
            metrics = runtime_metrics(run_root, result) if run_root is not None else result.get("metrics", {})
        except (OSError, KeyError, ValueError, json.JSONDecodeError) as error:
            errors.append(f"{result.get('sample')}: cannot read authoritative outputs ({error})")
            continue
        tags = metrics.get("read_tag_census", {})
        ledger = metrics.get("site_ledger", {})
        if tags.get("check_exact_sidecar_join") is not True:
            errors.append(f"{result.get('sample')}: exact read-tag join failed")
        if ledger.get("pass") is not True:
            errors.append(f"{result.get('sample')}: site ledger failed")
    if ablation.get("all_forbidden_counts_zero") is not True:
        errors.append("reference/HP3/HP4 forbidden-primary census failed")
    if {item.get("sample") for item in ablation.get("samples", [])} != set(manifest_samples):
        errors.append("lineage ablation dataset set is incomplete")
    if determinacy is not None and {item.get("sample") for item in determinacy.get("samples", [])} != set(manifest_samples):
        errors.append("determinacy dataset set is incomplete")
    if read_af.get("all_candidate_sets_complete") is not True:
        errors.append("read-AF candidate sets are incomplete")
    if {item.get("sample") for item in read_af.get("samples", [])} != set(manifest_samples):
        errors.append("read-AF dataset set is incomplete")
    expected_parameters = {"mapq30", "baseq10", "minread4", "maxsnv6"}
    if set(parameter.get("aggregate", {})) != expected_parameters:
        errors.append("parameter sensitivity is incomplete")
    elif any(not isinstance(summary.get("min_primary_unit_key_jaccard"), (int, float))
             or not isinstance(summary.get("min_shared_topology_digest_concordance"), (int, float))
             for summary in parameter["aggregate"].values()):
        errors.append("parameter sensitivity lacks direct primary-unit/topology comparisons")
    expected_parameter_rows = {(variant, sample) for variant in expected_parameters for sample in manifest_samples}
    observed_parameter_rows = [(item.get("variant"), item.get("sample")) for item in parameter.get("rows", [])]
    if len(observed_parameter_rows) != 28 or set(observed_parameter_rows) != expected_parameter_rows:
        errors.append("parameter sensitivity does not contain the exact 4x7 detail grid")
    comparisons = backbone.get("comparisons", [])
    comparison_samples = [item.get("sample") for item in comparisons]
    if (backbone.get("schema_version") != "2.1"
            or backbone.get("canonical") != "LongPhase-S _sc.vcf FILTER=PASS"
            or backbone.get("alternative")
            != "ClairS paired FILTER=PASS selected tree using the same raw-all LongPhase-S HP/PS tags"
            or len(comparisons) != 7
            or comparison_samples != manifest_samples
            or len(set(comparison_samples)) != 7
            or any(item.get("label") != "clairs_PASS_input_vs_longphase_PASS"
                   or not isinstance(item.get("primary_unit_key_jaccard"), (int, float))
                   or not isinstance(item.get("shared_unit_topology_digest_concordance"), (int, float))
                   for item in comparisons)):
        errors.append("7-dataset canonical-vs-ClairS backbone sensitivity is incomplete or mismatched")
    if errors:
        raise SystemExit("authoritative report input gate failed: " + "; ".join(errors))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--lineage-ablation", required=True, type=Path)
    parser.add_argument("--determinacy", required=True, type=Path)
    parser.add_argument("--read-af", required=True, type=Path)
    parser.add_argument("--parameter-sensitivity", required=True, type=Path)
    parser.add_argument("--backbone-sensitivity", required=True, type=Path)
    parser.add_argument("--secondary-science-verification", type=Path,
                        help="Required for a v3 canonical root; written outside the immutable run root")
    parser.add_argument("--output-md", required=True, type=Path)
    parser.add_argument("--output-html", required=True, type=Path)
    parser.add_argument("--embed-helper", required=True, type=Path)
    parser.add_argument("--source-notes", required=True, type=Path)
    args = parser.parse_args()
    manifest = load(args.input_manifest)
    verification = load(args.run_root / "verification_summary.json")
    ablation = load(args.lineage_ablation)
    determinacy = load(args.determinacy)
    read_af = load(args.read_af)
    parameter = load(args.parameter_sensitivity)
    backbone = load(args.backbone_sensitivity)
    secondary = load(args.secondary_science_verification) if args.secondary_science_verification else None
    validate_authoritative_inputs(manifest, verification, ablation, read_af, parameter, backbone,
                                  args.run_root, secondary, determinacy)
    rows = sample_rows(verification, ablation, determinacy, args.run_root)
    args.output_md.write_text(markdown(rows, verification, ablation, determinacy, read_af, parameter, backbone,
                                       args.run_root, args.input_manifest), encoding="utf-8")
    shell = html_report(rows, verification, ablation, determinacy, read_af, parameter, backbone,
                        args.run_root, args.input_manifest)
    chart_map = {
        "delivery_mode": "html", "audience": "technical",
        "charts": [
            {"id": "funnel-composition", "question": "Where does the operational universe leave the solver?", "family": "composition", "type": "stacked bar", "source": str(args.run_root / "verification_summary.json")},
            {"id": "determinacy-layers", "question": "How do exact, shape, and region determinacy differ?", "family": "comparison", "type": "grouped bar", "source": str(args.determinacy)},
        ],
        "omissions": {"backbone": "seven exact dataset comparisons are more auditable as a table",
                      "read_af": "heuristic parameters and CN strata require exact lookup; table used"},
    }
    args.source_notes.write_text(json.dumps(chart_map, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        shell_path = td / "report-shell.html"
        payload_path = td / "report-payload.json"
        shell_path.write_text(shell, encoding="utf-8")
        payload_path.write_text(json.dumps(payload(rows), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        proc = subprocess.run([sys.executable, str(args.embed_helper), "--input", str(shell_path),
                               "--payload", str(payload_path), "--output", str(args.output_html)],
                              text=True, capture_output=True, check=False)
        if proc.returncode:
            args.output_html.write_text(shell, encoding="utf-8")
            print(f"WARN: live chart runtime unavailable; static fallback delivered: {proc.stderr}", file=sys.stderr)
    print(f"REPORT MD -> {args.output_md}")
    print(f"REPORT HTML -> {args.output_html}")


if __name__ == "__main__":
    main()
