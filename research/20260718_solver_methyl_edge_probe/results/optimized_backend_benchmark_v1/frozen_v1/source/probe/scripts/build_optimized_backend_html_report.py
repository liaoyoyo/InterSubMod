#!/usr/bin/env python3
"""Build a standalone, source-bound HTML audit report for the optimized probe."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import pathlib
import platform
import time
from string import Template
from typing import Any


REPO = pathlib.Path("/big7_disk/liaoyoyo2001/InterSubMod")
PROBE_ROOT = REPO / "research/20260718_solver_methyl_edge_probe"
SOURCE_PATHS = {
    "probe": PROBE_ROOT / "scripts/solver_probe.py",
    "optimized_backend": PROBE_ROOT / "scripts/optimized_hypercube_backend.py",
    "current_solver": (
        REPO
        / "research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py"
    ),
    "fixtures": PROBE_ROOT / "tests/fixtures/real_units.json",
    "benchmark_harness": PROBE_ROOT / "scripts/run_optimized_backend_benchmark.py",
    "optimized_tests": PROBE_ROOT / "tests/test_optimized_hypercube_backend.py",
    "oracle_script": PROBE_ROOT / "scripts/verify_optimized_backend_oracles.py",
}
PLAN_PATH = PROBE_ROOT / "20260718_Hypercube邊與subcube改良研究計畫_01.md"
IMPLEMENTATION_NOTES_PATH = PROBE_ROOT / "implementation-notes.md"


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def seconds(value: float) -> str:
    if value < 0.001:
        return f"{value * 1_000_000:.1f} µs"
    if value < 1:
        return f"{value * 1000:.3f} ms"
    return f"{value:.3f} s"


def short_hash(value: str) -> str:
    return f"{value[:12]}…{value[-8:]}"


def assert_receipts(
    benchmark: dict[str, Any],
    oracle: dict[str, Any],
    benchmark_path: pathlib.Path,
    oracle_path: pathlib.Path,
) -> None:
    if benchmark.get("schema") != "intersubmod.optimized_backend_benchmark.suite.v1":
        raise ValueError("unexpected benchmark schema")
    if oracle.get("schema") != "intersubmod.optimized_backend_oracles.v1":
        raise ValueError("unexpected oracle schema")
    if benchmark["promotion_gate"]["overall"] != "PASS_FOR_BOUNDED_DUAL_PILOT":
        raise ValueError("benchmark receipt is not a bounded PASS")
    if not oracle["summary"]["all_pass"] or oracle["summary"]["total_mismatches"]:
        raise ValueError("oracle receipt is not mismatch-free")
    for comparison in benchmark["comparisons"]:
        if not comparison["exact_match"] or not comparison["repeat_stable"]:
            raise ValueError(f"comparison is not exact/stable: {comparison['case_id']}")
    expected = benchmark["source_sha256"]
    for key, path in SOURCE_PATHS.items():
        observed = sha256_file(path)
        if observed != expected[key]:
            raise ValueError(
                f"source drift for {key}: receipt={expected[key]} current={observed}"
            )
    if sha256_file(oracle_path) != expected["oracle_receipt"]:
        raise ValueError("oracle receipt hash differs from benchmark binding")
    if benchmark_path.stat().st_size == 0:
        raise ValueError("empty benchmark receipt")


def comparison_rows(comparisons: list[dict[str, Any]]) -> str:
    rows = []
    for row in comparisons:
        current = row["current"]["solver_elapsed_seconds"]
        optimized = row["optimized"]["solver_elapsed_seconds"]
        rss_current = row["current"]["ru_maxrss_kib"]
        rss_new = row["optimized"]["ru_maxrss_kib"]
        rows.append(
            "<tr>"
            f"<th scope='row'>{esc(row['case_id'])}</th>"
            f"<td>{row['h_star']}</td>"
            f"<td>{row['optimal_family_count']:,}</td>"
            f"<td><code>{esc(short_hash(row['canonical_family_digest']))}</code></td>"
            f"<td>{seconds(current['median'])}<small>p95 {seconds(current['p95'])}</small></td>"
            f"<td>{seconds(optimized['median'])}<small>p95 {seconds(optimized['p95'])}</small></td>"
            f"<td class='number-strong'>{row['median_speedup_solver']:.2f}×</td>"
            f"<td>{rss_current['median']/1024:.1f} → {rss_new['median']/1024:.1f} MiB"
            f"<small>{row['median_rss_ratio_optimized_over_current']:.3f}×</small></td>"
            "<td><span class='status pass'>EXACT</span></td>"
            "</tr>"
        )
    return "\n".join(rows)


def timing_bars(comparisons: list[dict[str, Any]]) -> str:
    cards = []
    for row in comparisons:
        current = row["current"]["solver_elapsed_seconds"]["median"]
        optimized = row["optimized"]["solver_elapsed_seconds"]["median"]
        new_width = max(1.5, 100.0 * optimized / current)
        cards.append(
            "<article class='timing-card'>"
            f"<div class='timing-head'><h3>{esc(row['case_id'])}</h3>"
            f"<strong>{row['median_speedup_solver']:.2f}×</strong></div>"
            "<div class='bar-label'><span>Current SciPy MILP</span>"
            f"<span>{seconds(current)}</span></div>"
            "<div class='track'><span class='bar current' style='width:100%'></span></div>"
            "<div class='bar-label'><span>DP + bitset B&amp;B</span>"
            f"<span>{seconds(optimized)}</span></div>"
            f"<div class='track'><span class='bar optimized' style='width:{new_width:.3f}%'></span></div>"
            f"<p>完整 family：{row['optimal_family_count']:,} 組；digest 逐組一致。</p>"
            "</article>"
        )
    return "\n".join(cards)


def source_rows(
    benchmark: dict[str, Any],
    benchmark_path: pathlib.Path,
    oracle_path: pathlib.Path,
) -> str:
    path_map = {
        **SOURCE_PATHS,
        "oracle_receipt": oracle_path,
        "benchmark_receipt": benchmark_path,
    }
    expected = {
        **benchmark["source_sha256"],
        "benchmark_receipt": sha256_file(benchmark_path),
    }
    rows = []
    for key, path in path_map.items():
        value = expected[key]
        relative = path.relative_to(REPO)
        rows.append(
            "<tr>"
            f"<th scope='row'>{esc(key)}</th>"
            f"<td><code>{esc(str(relative))}</code></td>"
            f"<td><code title='{esc(value)}'>{esc(short_hash(value))}</code></td>"
            f"<td>{path.stat().st_size:,} B</td>"
            "</tr>"
        )
    return "\n".join(rows)


def build_html(
    benchmark: dict[str, Any],
    oracle: dict[str, Any],
    benchmark_path: pathlib.Path,
    oracle_path: pathlib.Path,
) -> str:
    comparisons = benchmark["comparisons"]
    bench_receipt_sha = sha256_file(benchmark_path)
    oracle_receipt_sha = sha256_file(oracle_path)
    plan_sha = sha256_file(PLAN_PATH)
    notes_sha = sha256_file(IMPLEMENTATION_NOTES_PATH)
    generated = time.strftime("%Y-%m-%d %H:%M:%S %z")

    commands = """# Finite-domain oracle
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
python3 InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/verify_optimized_backend_oracles.py \\
  --output InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/oracle_receipt.json \\
  --seeded-k4-cases 300

# Five-repeat bounded comparison
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
python3 InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/run_optimized_backend_benchmark.py \\
  --output-dir InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/runs \\
  --receipt InterSubMod/research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/benchmark_receipt.json \\
  --import-current-dir /tmp/intersubmod_solver_baseline_20260718 \\
  --repeats 5 --warmups 1 --time-limit 30 --q-max 8

# Complete probe regression
INTERSUBMOD_HIGHSPY_PATH=/tmp/intersubmod_highspy_1_15_1 \\
python3 -m unittest discover \\
  -s InterSubMod/research/20260718_solver_methyl_edge_probe/tests \\
  -p 'test_*.py' -v"""

    benchmark_table = comparison_rows(comparisons)
    bars = timing_bars(comparisons)
    sources = source_rows(benchmark, benchmark_path, oracle_path)
    receipt_json = json.dumps(
        {
            "benchmark_receipt_sha256": bench_receipt_sha,
            "oracle_receipt_sha256": oracle_receipt_sha,
            "objective_status": benchmark["certificates"]["objective"]["status"],
            "family_status": benchmark["certificates"]["candidate_family"]["status"],
            "oracle_cases": oracle["summary"]["total_structural_cases"],
            "oracle_mismatches": oracle["summary"]["total_mismatches"],
        },
        ensure_ascii=False,
    ).replace("<", "\\u003c")

    template = Template(r"""<!doctype html>
<html lang="zh-Hant" data-theme="paper">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="light dark">
  <meta name="description" content="InterSubMod Hypercube optimized backend bounded implementation, exactness and performance audit">
  <title>Hypercube Optimized Backend｜完整流程、正確性與效能稽核</title>
  <style>
    :root {
      --paper:#f3efe4; --paper-2:#e9e1cf; --ink:#17252d; --muted:#58666d;
      --navy:#173f52; --navy-2:#0e2a37; --rust:#bd4f2f; --amber:#dba646;
      --sage:#507565; --line:#b9ae98; --white:#fffdf7; --shadow:rgba(23,37,45,.13);
      --serif:"Iowan Old Style","Noto Serif TC","Songti TC",Georgia,serif;
      --sans:"Avenir Next","Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif;
      --mono:"IBM Plex Mono","SFMono-Regular","Cascadia Code","Noto Sans Mono",monospace;
    }
    html[data-theme="night"] {
      --paper:#111b20; --paper-2:#19272d; --ink:#e9e5d8; --muted:#aeb9b9;
      --navy:#78b7c6; --navy-2:#d5edf0; --rust:#ef8c69; --amber:#f0c36a;
      --sage:#8fc3aa; --line:#43565a; --white:#17242a; --shadow:rgba(0,0,0,.3);
    }
    * { box-sizing:border-box; }
    html { scroll-behavior:smooth; }
    body {
      margin:0; color:var(--ink); background:
        linear-gradient(90deg,transparent 0 47px,rgba(189,79,47,.16) 48px 49px,transparent 50px),
        repeating-linear-gradient(0deg,transparent 0 31px,rgba(23,63,82,.045) 32px 33px),
        var(--paper);
      font-family:var(--sans); line-height:1.72;
    }
    a { color:var(--navy); text-underline-offset:3px; }
    a:focus-visible,button:focus-visible,summary:focus-visible { outline:3px solid var(--amber); outline-offset:3px; }
    .skip { position:fixed; left:1rem; top:-8rem; z-index:100; background:var(--ink); color:var(--paper); padding:.7rem 1rem; }
    .skip:focus { top:1rem; }
    .scope-ribbon {
      position:sticky; top:0; z-index:40; padding:.52rem clamp(1rem,4vw,4rem);
      display:flex; gap:1rem; justify-content:space-between; align-items:center;
      background:var(--rust); color:#fff9eb; font:700 .78rem/1.3 var(--mono); letter-spacing:.08em;
      box-shadow:0 3px 20px var(--shadow);
    }
    .scope-ribbon button,.toolbar button {
      border:1px solid currentColor; color:inherit; background:transparent; border-radius:999px;
      padding:.34rem .72rem; cursor:pointer; font:inherit;
    }
    header.hero {
      min-height:72vh; display:grid; grid-template-columns:minmax(0,1.55fr) minmax(250px,.7fr);
      gap:clamp(2rem,6vw,7rem); align-items:end; padding:clamp(4rem,9vw,9rem) clamp(2rem,8vw,9rem) 4rem;
      border-bottom:2px solid var(--ink); position:relative; overflow:hidden;
    }
    header.hero::after {
      content:"Qₘ"; position:absolute; right:-.03em; top:-.2em; font:900 clamp(12rem,36vw,36rem)/1 var(--serif);
      color:var(--navy); opacity:.045; pointer-events:none;
    }
    .eyebrow { font:700 .78rem/1.3 var(--mono); letter-spacing:.16em; color:var(--rust); text-transform:uppercase; }
    h1,h2,h3 { font-family:var(--serif); text-wrap:balance; }
    h1 { font-size:clamp(3rem,7.5vw,7.6rem); line-height:.91; letter-spacing:-.055em; margin:.45rem 0 1.5rem; max-width:12ch; }
    h1 em { color:var(--rust); font-style:italic; }
    .lede { font:500 clamp(1.05rem,1.8vw,1.4rem)/1.65 var(--serif); max-width:54rem; }
    .verdict-plate { border:2px solid var(--ink); background:var(--white); padding:1.4rem; box-shadow:12px 12px 0 var(--navy); position:relative; z-index:1; }
    .verdict-plate .stamp { display:inline-block; transform:rotate(-2deg); border:3px double var(--sage); color:var(--sage); padding:.25rem .65rem; font:800 .78rem var(--mono); }
    .verdict-plate strong { display:block; font:700 2.35rem/1 var(--serif); margin:1.2rem 0 .4rem; }
    .hero-meta { margin:1.4rem 0 0; padding:0; list-style:none; font:.78rem/1.65 var(--mono); color:var(--muted); }
    .layout { display:grid; grid-template-columns:260px minmax(0,1fr); max-width:1560px; margin:auto; }
    nav.toc { position:sticky; top:3rem; align-self:start; height:calc(100vh - 3rem); padding:3rem 1.5rem 2rem 2rem; border-right:1px solid var(--line); overflow:auto; }
    nav.toc strong { font:800 .72rem var(--mono); letter-spacing:.14em; color:var(--rust); }
    nav.toc ol { list-style:none; padding:0; margin:1rem 0; counter-reset:nav; }
    nav.toc li { counter-increment:nav; margin:.24rem 0; }
    nav.toc a { display:block; padding:.45rem .55rem; border-left:3px solid transparent; color:var(--muted); text-decoration:none; font-size:.87rem; }
    nav.toc a::before { content:counter(nav,decimal-leading-zero) " / "; font-family:var(--mono); color:var(--line); }
    nav.toc a.active,nav.toc a:hover { color:var(--ink); border-left-color:var(--rust); background:var(--paper-2); }
    .toolbar { display:flex; gap:.5rem; flex-wrap:wrap; }
    .toolbar button { color:var(--ink); font:.72rem var(--mono); }
    main { min-width:0; padding:0 clamp(1.4rem,5vw,6rem) 7rem; }
    section { padding:5.5rem 0; border-bottom:1px solid var(--line); scroll-margin-top:4rem; }
    .section-kicker { font:800 .74rem var(--mono); color:var(--rust); letter-spacing:.14em; text-transform:uppercase; }
    h2 { font-size:clamp(2.2rem,4.6vw,4.8rem); line-height:1; letter-spacing:-.035em; margin:.45rem 0 2rem; max-width:15ch; }
    h3 { font-size:1.55rem; line-height:1.2; }
    .answer {
      display:grid; grid-template-columns:minmax(0,1.6fr) minmax(230px,.65fr); gap:2rem;
      padding:2rem; border:2px solid var(--ink); background:var(--white); box-shadow:9px 9px 0 var(--amber);
    }
    .answer p:first-child { font:700 clamp(1.35rem,2.5vw,2.15rem)/1.28 var(--serif); margin:0; }
    .answer aside { border-left:1px solid var(--line); padding-left:1.5rem; }
    .metric-grid,.certificate-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:1rem; margin:2rem 0; }
    .metric,.certificate { background:var(--white); border:1px solid var(--line); padding:1.25rem; min-height:150px; }
    .metric span,.certificate span { font:700 .67rem var(--mono); color:var(--muted); letter-spacing:.08em; text-transform:uppercase; }
    .metric strong { display:block; font:800 clamp(2rem,4vw,3.8rem)/1 var(--serif); margin:.7rem 0; }
    .certificate { border-top:5px solid var(--sage); }
    .certificate.warn { border-top-color:var(--amber); }
    .certificate.stop { border-top-color:var(--rust); }
    .certificate strong { display:block; margin:.65rem 0; font-size:1.05rem; }
    .status { display:inline-flex; align-items:center; border-radius:999px; padding:.22rem .58rem; font:800 .67rem var(--mono); letter-spacing:.05em; }
    .status.pass { background:rgba(80,117,101,.16); color:var(--sage); }
    .status.warn { background:rgba(219,166,70,.19); color:#7b5610; }
    .status.stop { background:rgba(189,79,47,.14); color:var(--rust); }
    .flow { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:1px; background:var(--line); border:1px solid var(--line); margin:2.5rem 0; }
    .flow article { background:var(--white); padding:1.35rem; position:relative; min-height:190px; }
    .flow article::before { content:attr(data-step); font:900 3rem/1 var(--serif); color:var(--rust); opacity:.38; }
    .flow h3 { margin:.5rem 0; }
    .flow code { color:var(--navy); }
    .decision-table,.data-table { width:100%; border-collapse:collapse; margin:1.6rem 0; background:var(--white); font-size:.88rem; }
    .decision-table th,.decision-table td,.data-table th,.data-table td { border:1px solid var(--line); padding:.75rem .8rem; text-align:left; vertical-align:top; }
    .decision-table thead th,.data-table thead th { background:var(--navy-2); color:var(--paper); font:700 .72rem var(--mono); letter-spacing:.05em; position:sticky; top:2.1rem; }
    html[data-theme="night"] .decision-table thead th,html[data-theme="night"] .data-table thead th { color:#111b20; }
    .data-table td small { display:block; color:var(--muted); margin-top:.15rem; }
    .number-strong { font:800 1.15rem var(--mono); color:var(--rust); }
    code,pre { font-family:var(--mono); }
    code { overflow-wrap:anywhere; }
    pre { background:var(--navy-2); color:#e8f0ed; padding:1.2rem; border-radius:2px; overflow:auto; font-size:.76rem; line-height:1.6; }
    .formula { display:block; padding:1.1rem 1.4rem; background:var(--paper-2); border-left:5px solid var(--navy); font:600 1rem/1.65 var(--mono); overflow:auto; }
    .callout { border-left:6px solid var(--amber); background:var(--white); padding:1rem 1.3rem; margin:1.5rem 0; }
    .callout.stop { border-left-color:var(--rust); }
    .timing-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:1rem; }
    .timing-card { background:var(--white); border:1px solid var(--line); padding:1.25rem; }
    .timing-head { display:flex; justify-content:space-between; align-items:baseline; }
    .timing-head h3 { margin:.15rem 0 1rem; }
    .timing-head strong { color:var(--rust); font:900 1.5rem var(--mono); }
    .bar-label { display:flex; justify-content:space-between; gap:1rem; font:.7rem var(--mono); margin:.7rem 0 .25rem; }
    .track { height:10px; background:var(--paper-2); overflow:hidden; }
    .bar { display:block; height:100%; min-width:3px; transform-origin:left; animation:grow .9s ease-out both; }
    .bar.current { background:var(--line); }
    .bar.optimized { background:var(--sage); }
    @keyframes grow { from { transform:scaleX(0); } }
    details { border:1px solid var(--line); background:var(--white); margin:.7rem 0; }
    summary { padding:1rem 1.2rem; cursor:pointer; font-weight:700; }
    details > div { padding:0 1.2rem 1.2rem; }
    .two-col { display:grid; grid-template-columns:1fr 1fr; gap:1.5rem; }
    .checklist { list-style:none; padding:0; }
    .checklist li { padding:.55rem 0 .55rem 2rem; border-bottom:1px dotted var(--line); position:relative; }
    .checklist li::before { content:"✓"; position:absolute; left:.35rem; color:var(--sage); font-weight:900; }
    .checklist li.stop::before { content:"×"; color:var(--rust); }
    .hash-note { font:.72rem var(--mono); color:var(--muted); overflow-wrap:anywhere; }
    footer { padding:3rem clamp(2rem,8vw,8rem); background:var(--navy-2); color:#e7ebe3; }
    footer a { color:#dceecf; }
    .print-only { display:none; }
    @media (max-width:1050px) {
      header.hero { grid-template-columns:1fr; min-height:auto; }
      .layout { grid-template-columns:1fr; }
      nav.toc { position:relative; top:auto; height:auto; border-right:0; border-bottom:1px solid var(--line); padding:1.5rem; }
      nav.toc ol { columns:2; }
      .metric-grid,.certificate-grid { grid-template-columns:repeat(2,1fr); }
      .flow,.timing-grid { grid-template-columns:1fr 1fr; }
    }
    @media (max-width:680px) {
      body { background:var(--paper); }
      .scope-ribbon { align-items:flex-start; }
      header.hero { padding:4rem 1.25rem 2.5rem; }
      h1 { font-size:3.25rem; }
      main { padding:0 1.1rem 4rem; }
      section { padding:4rem 0; }
      nav.toc ol { columns:1; }
      .answer,.two-col { grid-template-columns:1fr; }
      .answer aside { border-left:0; border-top:1px solid var(--line); padding:1rem 0 0; }
      .metric-grid,.certificate-grid,.flow,.timing-grid { grid-template-columns:1fr; }
      .table-wrap { overflow-x:auto; }
    }
    @media (prefers-reduced-motion:reduce) {
      html { scroll-behavior:auto; }
      *,*::before,*::after { animation:none!important; transition:none!important; }
    }
    @media print {
      .scope-ribbon,nav.toc,.toolbar { display:none!important; }
      .layout { display:block; }
      header.hero { min-height:auto; padding:1cm; }
      main { padding:0; }
      section { break-inside:avoid; padding:1cm 0; }
      body { background:white; color:black; font-size:10pt; }
      .print-only { display:block; }
      a { color:black; text-decoration:none; }
    }
  </style>
</head>
<body>
  <a class="skip" href="#main">跳到主要內容</a>
  <div class="scope-ribbon">
    <span>PARTIAL · BOUNDED PROBE · NON-PRODUCTION · production_claim_allowed=false</span>
    <button id="themeToggle" type="button" aria-label="切換明暗主題">切換主題</button>
  </div>

  <header class="hero">
    <div>
      <p class="eyebrow">InterSubMod / Exact structural solver audit / 2026-07-18</p>
      <h1>從反覆求解到<em>一次證明</em></h1>
      <p class="lede">Bitset obligation B&amp;B、動態 subcube antichain、small-q subset DP 與 fixed-N parent mapping 的完整流程、判斷細節、正確性證據與本機效能比較。</p>
    </div>
    <aside class="verdict-plate" aria-label="總判定">
      <span class="stamp">PASS FOR BOUNDED DUAL PILOT</span>
      <strong>結果相同；局部明顯加速</strong>
      <p>三個完整案例的 <code>h*</code>、candidate 數與 canonical family digest 全相同。這允許下一階段 stress panel，不允許替換 canonical 或啟動 production full run。</p>
      <ul class="hero-meta">
        <li>Benchmark receipt: $benchmark_receipt_short</li>
        <li>Oracle receipt: $oracle_receipt_short</li>
        <li>Generated: $generated</li>
      </ul>
    </aside>
  </header>

  <div class="layout">
    <nav class="toc" aria-label="報告目錄">
      <strong>DOCUMENT MAP</strong>
      <ol>
        <li><a href="#verdict">直答與四證書</a></li>
        <li><a href="#scope">範圍與 provenance</a></li>
        <li><a href="#flow">完整處理流程</a></li>
        <li><a href="#methods">基本判斷與方法</a></li>
        <li><a href="#before-after">之前與現在</a></li>
        <li><a href="#results">結果與時間</a></li>
        <li><a href="#oracles">Oracle 與 fail-closed</a></li>
        <li><a href="#edge-biology">Edge／chain／biology</a></li>
        <li><a href="#files">程式與命令</a></li>
        <li><a href="#gate">限制與 promotion gate</a></li>
      </ol>
      <div class="toolbar">
        <button id="expandAll" type="button">展開附錄</button>
        <button id="collapseAll" type="button">收合附錄</button>
        <button type="button" onclick="window.print()">列印</button>
      </div>
    </nav>

    <main id="main">
      <section id="verdict">
        <p class="section-kicker">01 / Executive verdict</p>
        <h2>可以說什麼，不能說什麼</h2>
        <div class="answer">
          <div>
            <p>可以確認：在 AAAA、H2009_M31、COLO829_M31 三個 bounded complete cases，新舊輸出逐組相同；新 backend 的本機 solver median 快 $speed_range，isolated-process peak RSS 為 current 的 $rss_range。</p>
            <p>不能確認：33-tail、全樣本、production SLA、真實 parent ancestry 或生物準確率。這些沒有被本輪數據覆蓋。</p>
          </div>
          <aside>
            <span class="status pass">影響：高</span>
            <span class="status pass">structural 信心：高</span>
            <span class="status warn">production 信心：待驗證</span>
            <p class="hash-note">Task classification：A exploratory implementation，附 B comprehensive-style exactness gate；服務 G4/G5。</p>
          </aside>
        </div>
        <div class="certificate-grid">
          <article class="certificate"><span>Objective certificate</span><strong>PROVEN_OPTIMAL</strong><p>DP 或完整 B&amp;B 證明 minimum-extra <code>h*</code>。</p></article>
          <article class="certificate"><span>Candidate family</span><strong>COMPLETE / 3 cases</strong><p>完整 family count 與 digest 一致；30/30 timing runs complete。</p></article>
          <article class="certificate warn"><span>Edge certificate</span><strong>ANALYTIC ORACLE PASS</strong><p>固定 N 的 additive score 正確；real edge evidence 未測。</p></article>
          <article class="certificate stop"><span>Biology certificate</span><strong>NOT TESTED / GATED</strong><p>M1只作CN/LOH-clean sensitivity；M2維持 unresolved。</p></article>
        </div>
        <div class="metric-grid">
          <article class="metric"><span>Finite structural cases</span><strong>$oracle_cases</strong><p>DP/B&amp;B 對 brute-force mismatch = 0。</p></article>
          <article class="metric"><span>Five-repeat raw runs</span><strong>30/30</strong><p>Current + optimized 全部完整。</p></article>
          <article class="metric"><span>Digest mismatches</span><strong>0</strong><p>三案每次 repeat 均穩定。</p></article>
          <article class="metric"><span>Incomplete ranked</span><strong>0</strong><p>cap/deadline 一律阻擋 winner。</p></article>
        </div>
      </section>

      <section id="scope">
        <p class="section-kicker">02 / Scope & provenance</p>
        <h2>先固定證據，再談速度</h2>
        <div class="two-col">
          <div>
            <h3>本輪實際涵蓋</h3>
            <ul class="checklist">
              <li>AAAA <code>k=4</code> exhaustive family（24組）。</li>
              <li>H2009_M31（<code>h*=5, V=242</code>）。</li>
              <li>COLO829_M31（<code>h*=5, V=216</code>）。</li>
              <li><code>k≤3</code> exhaustive + 300 seeded <code>k=4</code> + noncontiguous <code>k=12</code>。</li>
            </ul>
          </div>
          <div>
            <h3>本輪明確未涵蓋</h3>
            <ul class="checklist">
              <li class="stop">33 個歷史 solver tails。</li>
              <li class="stop">H2009/H1437 stress panel與全7 technical datasets。</li>
              <li class="stop">canonical router、production tagged BAM或full run。</li>
              <li class="stop">真實 methylation edge、CN-corrected CCF或ancestry truth。</li>
            </ul>
          </div>
        </div>
        <div class="callout">
          <strong>共享 workspace drift 如何處理：</strong>研究計畫曾由外部並行 session 更新，因此 plan hash只列為 mutable planning context；可執行證據由 solver、fixture、tests、oracle與receipts的byte hashes獨立綁定。
        </div>
        <div class="callout">
          <strong>Frozen bounded evidence：</strong>本報告完成瀏覽器QA後，連同source、fixture、
          current baseline raw records、optimized raw records、receipts與QA evidence，以
          no-clobber方式封存在
          <code>results/optimized_backend_benchmark_v1/frozen_v1/</code>；
          <code>manifest.json.sha256</code>是外部完整性sidecar。這不等於未來33-tail輸入也已freeze。
        </div>
        <div class="table-wrap">
          <table class="data-table">
            <thead><tr><th>Role</th><th>Path</th><th>SHA-256</th><th>Bytes</th></tr></thead>
            <tbody>$source_rows</tbody>
          </table>
        </div>
        <p class="hash-note">Mutable plan snapshot：<code>$plan_sha</code><br>Implementation notes snapshot：<code>$notes_sha</code><br>OS：$platform</p>
      </section>

      <section id="flow">
        <p class="section-kicker">03 / End-to-end flow</p>
        <h2>從 R/A/X 到可排名 family</h2>
        <div class="flow">
          <article data-step="01"><h3>Freeze contract</h3><p>綁定input、source SHA、thread env、deadline與<code>max_sets</code>；未綁定即不比較。</p></article>
          <article data-step="02"><h3>Active-bit compression</h3><p>raw <code>k → effective m</code>；只保留真正出現ALT的維度，不混淆raw座標與dense index。</p></article>
          <article data-step="03"><h3>Exact group reduction</h3><p>duplicate、mandatory-hit、singleton、subset dominance與downward closure；scoring counts不刪。</p></article>
          <article data-step="04"><h3>Subset DP</h3><p><code>q≤8</code>且cells/bytes/ops gate通過時，一次證明<code>h*</code>；只發objective certificate。</p></article>
          <article data-step="05"><h3>Bitset B&amp;B</h3><p>動態group/parent obligations、antichain、singleton propagation、MRV與safe LB；列完整optimal family。</p></article>
          <article data-step="06"><h3>Family gate</h3><p>只有<code>family_complete=true</code>才允許read-likelihood ranking；cap/deadline立即abstain。</p></article>
          <article data-step="07"><h3>Fixed-N parents</h3><p>edge-local additive分數以<code>O(E_N)</code>一次算best/ties/posterior；不展開Cartesian product。</p></article>
          <article data-step="08"><h3>Unary projection</h3><p>無證據chain只作multi-mutation display equivalence，順序數<code>n!</code>；保留原candidate。</p></article>
          <article data-step="09"><h3>Biology gate</h3><p>M1要求1+1 allele-specific CN與完整LOH/CNA/QC；M2未實作則unresolved，不硬猜。</p></article>
        </div>
      </section>

      <section id="methods">
        <p class="section-kicker">04 / Decisions at the lowest level</p>
        <h2>每一個判斷條件</h2>
        <div class="table-wrap">
          <table class="decision-table">
            <thead><tr><th>判斷</th><th>True 時</th><th>False／不明時</th><th>Exactness理由</th></tr></thead>
            <tbody>
              <tr><th>Group 已被 required vertex 命中？</th><td>移除structural row</td><td>保留</td><td>required 永遠selected；不影響likelihood counts。</td></tr>
              <tr><th><code>Gₐ ⊆ Gᵦ</code>？</th><td>移除superset <code>Gᵦ</code></td><td>兩者保留</td><td>命中subset必然命中superset；反向移除會錯。</td></tr>
              <tr><th>Obligation domain size = 1？</th><td>強制加入唯一vertex</td><td>交給MRV分枝</td><td>任何可行解都必選它；非mandatory仍計cost 1。</td></tr>
              <tr><th>Selected child有selected predecessor？</th><td>parent obligation已滿足</td><td>新增<code>Pred(v)\excluded</code> obligation</td><td>rank每步下降，local parent條件即root connectivity。</td></tr>
              <tr><th>某candidate仍可避開excluded走到root？</th><td>保留domain</td><td>從domain移除；空則infeasible</td><td>只刪除未來不可能root-connected的state。</td></tr>
              <tr><th><code>LB &gt; best</code>？</th><td>剪枝</td><td>繼續；<code>LB==best</code>不可剪</td><td>保留所有同分optimal ties。</td></tr>
              <tr><th>DP route cells/bytes/ops超標？</th><td>拒絕DP route</td><td>配置table並求objective</td><td>resource gate不截terminal；可退回exact B&amp;B。</td></tr>
              <tr><th>Family完整？</th><td>允許每個distinct N做一次ranking</td><td><strong>ABSTAIN</strong></td><td>candidate prefix不能代表全部minimum family。</td></tr>
              <tr><th>Edge score可分解為child-local sum？</th><td>parent mapping解析求best</td><td>拒絕analytic route</td><td>degree/path/sibling coupling必須joint solver。</td></tr>
              <tr><th>M1所有CN/LOH/QC gate通過，且每個non-root state有selected predecessor？</th><td>strict只作sensitivity</td><td><code>ABSTAIN_CN_LOH_GATE</code>或<code>STRICT_INFEASIBLE</code></td><td>2+0 cnLOH、loss、amplification與不連root topology皆不可當strict-compatible。</td></tr>
            </tbody>
          </table>
        </div>

        <details open><summary>Bitset obligation B&amp;B：完整狀態與分枝</summary><div>
          <p class="formula">State = (selected_bits, excluded_bits)<br>Q = uncovered groups ∪ orphan predecessor domains<br>branch chosen D as: include c₁; exclude c₁+include c₂; …</p>
          <p>Sibling exclusion將解空間分成互斥區塊；memo key保留selected與excluded，不以「未滿足domain看起來相同」錯誤合併。Lower bound取disjoint packing、coverage與root-connection的最大值，不能相加，因同一vertex可能同時滿足多個義務。</p>
        </div></details>

        <details><summary>Small-q DP：為何merge不用再扣junction成本</summary><div>
          <p class="formula">DP[S,v] = 從 v 往下覆蓋 terminals S 的最小成本，且不含 c(v)<br>merge: DP[A,v] + DP[S\A,v]<br>edge relax: c(child) + DP[S,child]</p>
          <p>因為<code>v</code>本身成本從未放入兩個subproblem，merge不會重複支付；parent→child時才支付child。27,054個exhaustive cases與shared-connector反例均支持此 convention。</p>
        </div></details>

        <details><summary>Dynamic antichain：只做existential structural reduction</summary><div>
          <p>若當下domain <code>D₁ ⊆ D₂</code>，命中D₁即同時命中D₂，所以branching只需保留inclusion-minimal domains。這不代表兩條reads、兩個child identities或edge evidence可合併；原始groups/counts仍保留在structural以外的層。</p>
        </div></details>
      </section>

      <section id="before-after">
        <p class="section-kicker">05 / Before → After → Impact</p>
        <h2>改變的是求解方式，不是 estimand</h2>
        <div class="table-wrap">
          <table class="decision-table">
            <thead><tr><th>面向</th><th>之前</th><th>現在</th><th>影響</th></tr></thead>
            <tbody>
              <tr><th>Objective proof</th><td>Current MILP先求一解；舊one-pass B&amp;B遇cap時可能只有incumbent。</td><td>small-q DP先獨立證明<code>h*</code>；證書欄位分離。</td><td>避免把cost 3 incumbent誤稱真<code>h*=2</code>。</td></tr>
              <tr><th>All-optimal family</th><td>Current需約<code>V+1</code>次MILP solves。</td><td>Bitset B&amp;B一次traversal同時列family。</td><td>兩個M31由十幾/二十多秒降至數百毫秒。</td></tr>
              <tr><th>Constraint representation</th><td>Python sets/frozensets與重複domain掃描。</td><td>integer bitset、dynamic antichain、dense-index mapping。</td><td>AND/popcount為主要操作；noncontiguous raw bits不截斷。</td></tr>
              <tr><th>Resource policy</th><td>Current 30秒是每個MILP solve，unit可累積很多次。</td><td>30秒包prepare＋DP＋B&amp;B；DP另有cells/bytes/ops gate。</td><td>更接近per-unit fail-closed，但尚非production worker SLA。</td></tr>
              <tr><th>Parent edges</th><td>已有tree count；read likelihood只依N。</td><td>補best/tie/posterior解析mapping。</td><td>可避免parent Cartesian展開；不創造新read evidence。</td></tr>
              <tr><th>Unsupported unary chain</th><td>可能展示一條任意完整順序。</td><td>multi-mutation edge＋<code>n!</code>未定序數。</td><td>降低假精確；仍是projection-only model boundary。</td></tr>
              <tr><th>Infinite sites</th><td>容易被誤當普遍加速限制。</td><td>只在1+1 CN與全QC通過時作M1 sensitivity。</td><td>CNA/LOH區域維持M0或M2 unresolved。</td></tr>
            </tbody>
          </table>
        </div>
        <div class="callout stop"><strong>沒有改的資料：</strong>BAM reads、R/A/X pattern counts、methylation、VAF/CCF、CN/LOH calls與primary likelihood都沒有重算。本輪證明的是structural候選求解等價與本機效能，不是F1或biology改善。</div>
      </section>

      <section id="results">
        <p class="section-kicker">06 / Exactness & local performance</p>
        <h2>結果相同，時間明顯縮短</h2>
        <div class="timing-grid">$timing_bars</div>
        <div class="table-wrap">
          <table class="data-table">
            <thead><tr><th>Case</th><th>h*</th><th>V</th><th>Canonical digest</th><th>Current median</th><th>New median</th><th>Solver speedup</th><th>Peak RSS median</th><th>Gate</th></tr></thead>
            <tbody>$benchmark_rows</tbody>
          </table>
        </div>
        <div class="callout">
          <strong>「結果相同」的精確定義：</strong>同一input、同一M0 recurrence-allowed objective下，<code>h*</code>相同、optimal family count相同，且將每個vertex set排序後個別hash、再組成family hash的canonical digest相同。不是只比較第一組解，也不是只比較candidate數。
        </div>
        <p class="hash-note">Timing：每case先warmup 1次，再正式5次；threads固定1。Current baseline與new run為順序執行、非interleaved。RSS是隔離process peak，不是production pipeline增量。</p>
      </section>

      <section id="oracles">
        <p class="section-kicker">07 / Independent oracles & negative tests</p>
        <h2>如何知道沒有錯剪枝</h2>
        <div class="metric-grid">
          <article class="metric"><span>Exhaustive k≤3</span><strong>27,054</strong><p>mandatory subsets × 0–2 subcube groups。</p></article>
          <article class="metric"><span>Seeded k=4</span><strong>300</strong><p>Objective與完整family對brute force。</p></article>
          <article class="metric"><span>Noncontiguous k=12</span><strong>1</strong><p>raw states 0/1/2048/2049 mapping正確。</p></article>
          <article class="metric"><span>Total mismatch</span><strong>0</strong><p>合計 $oracle_cases cases。</p></article>
        </div>
        <div class="two-col">
          <div>
            <h3>正向 controls</h3>
            <ul class="checklist">
              <li><code>AAA</code>保留6條optimal chains。</li>
              <li><code>AAAA</code>保留24組family。</li>
              <li><code>AX + XA</code>聯合group coverage為3組。</li>
              <li>Shared connector的DP node cost只付一次。</li>
              <li>Diamond equal edge score保留2個best parents。</li>
            </ul>
          </div>
          <div>
            <h3>Fail-closed controls</h3>
            <ul class="checklist">
              <li><code>max_sets=1</code>且V=24：family incomplete、ranking=false。</li>
              <li>total deadline先到：objective不冒充已證明。</li>
              <li>無incumbent：<code>NO_FEASIBLE_CERTIFICATE_*</code>。</li>
              <li>missing/nonfinite edge score與非法tolerance直接拒絕。</li>
              <li>單edge或跨child累加overflow直接拒絕。</li>
              <li>M1缺root、state越界、不連root或CN矛盾直接reject/infeasible/abstain。</li>
            </ul>
          </div>
        </div>
      </section>

      <section id="edge-biology">
        <p class="section-kicker">08 / Edge, unary chain & biology</p>
        <h2>計算得出，不代表生物可識別</h2>
        <h3>Fixed-N parent mapping</h3>
        <p class="formula">T(N) = ∏ |Pred_N(v)|<br>S*(N) = Σ maxₚ s(p,v)<br>best-tree count = ∏ |argmaxₚ s(p,v)|</p>
        <p>只在score是child-local additive時成立。若加入degree cap、sibling competition、CCF sum rule、strict gain-once或path-level methylation consistency，就必須拒絕解析route並改用joint edge solver。</p>

        <h3>Unary hidden chain</h3>
        <p>以固定parent tree為單位，若中間connector沒有read/evidence且每點只有一個child，可顯示為一條<code>MULTI_MUTATION_EDGE_EQUIVALENCE</code>；三個mutations時未定序數是<code>3!=6</code>。原始candidate與objective仍保留，報告不宣稱這是跨所有candidate的壓縮family或objective-preserving新模型。</p>

        <h3>M0／M1／M2</h3>
        <div class="table-wrap">
          <table class="decision-table">
            <thead><tr><th>Mode</th><th>允許條件</th><th>目前輸出</th><th>可說的上限</th></tr></thead>
            <tbody>
              <tr><th>M0 recurrence-allowed</th><td>所有structural units</td><td>Baseline exact family</td><td>rooted monotone molecule-state topology</td></tr>
              <tr><th>M1 strict infinite-sites</th><td>allele-specific 1+1 CN；無clonal/subclonal LOH/deletion/amplification/WGD/duplicate copy；read與somatic QC全過；selected states連根</td><td>STRICT_COMPATIBLE / INFEASIBLE / ABSTAIN</td><td>strict-perfect-compatible molecule-state sensitivity</td></tr>
              <tr><th>M2 loss-supported Dollo</th><td>需明確copy/loss state與evidence</td><td>UNRESOLVED_NOT_IMPLEMENTED</td><td>不可發表biology claim</td></tr>
            </tbody>
          </table>
        </div>
      </section>

      <section id="files">
        <p class="section-kicker">09 / Code map & reproducibility</p>
        <h2>改了哪些檔，怎麼重播</h2>
        <div class="table-wrap">
          <table class="decision-table">
            <thead><tr><th>檔案</th><th>角色</th><th>變動</th></tr></thead>
            <tbody>
              <tr><th><code>scripts/optimized_hypercube_backend.py</code></th><td>DP、bitset B&amp;B、parent mapping、chain projection、M1/M2 gates</td><td>新增；isolated research backend</td></tr>
              <tr><th><code>tests/test_optimized_hypercube_backend.py</code></th><td>structural/edge/biology/failure regression</td><td>新增</td></tr>
              <tr><th><code>scripts/verify_optimized_backend_oracles.py</code></th><td>27,355-case finite oracle</td><td>新增</td></tr>
              <tr><th><code>scripts/run_optimized_backend_benchmark.py</code></th><td>isolated repeats、digest與receipt</td><td>新增</td></tr>
              <tr><th><code>tests/test_solver_probe.py</code></th><td>prepared-base build diagnostic</td><td>stale 25-build assertion修為current 1</td></tr>
              <tr><th><code>implementation-notes.md</code></th><td>living decisions／deviations／commands</td><td>追加本輪證據</td></tr>
              <tr><th><code>scripts/build_optimized_backend_html_report.py</code></th><td>source-bound standalone HTML</td><td>新增</td></tr>
              <tr><th><code>scripts/freeze_optimized_backend_evidence.py</code></th><td>no-clobber copy、manifest、sidecars、read-only verify</td><td>新增</td></tr>
              <tr><th><code>hypercube_exact.py</code></th><td>current comparison backend</td><td><strong>未修改</strong></td></tr>
            </tbody>
          </table>
        </div>
        <details open><summary>完整重播命令</summary><div><pre><code>$commands</code></pre></div></details>
        <p class="hash-note">最新 regression：本 probe 33/33 PASS；未修改的 current
        <code>test_hypercube_exact.py</code> 27/27 PASS。</p>
        <p>主要輸出：</p>
        <ul>
          <li><code>results/optimized_backend_benchmark_v1/benchmark_receipt.json</code></li>
          <li><code>results/optimized_backend_benchmark_v1/oracle_receipt.json</code></li>
          <li><code>results/optimized_backend_benchmark_v1/runs/optimized_dp_bitset_bnb/</code></li>
        </ul>
      </section>

      <section id="gate">
        <p class="section-kicker">10 / Promotion decision</p>
        <h2>下一步被授權，但 production 未被授權</h2>
        <div class="two-col">
          <div>
            <h3>本輪已通過</h3>
            <ul class="checklist">
              <li>Complete controls：digest mismatch 0。</li>
              <li>Objective/family certificate分離。</li>
              <li>Incomplete family不產生winner。</li>
              <li>DP resource gate與total-unit deadline。</li>
              <li>Local timing與RSS repeat evidence。</li>
            </ul>
          </div>
          <div>
            <h3>Promotion 前仍需</h3>
            <ul class="checklist">
              <li class="stop">33-tail／stress panel 各自的新 frozen input snapshot與sidecar。</li>
              <li class="stop">33-tail + H2009/H1437 stress panel。</li>
              <li class="stop">Total wall ≥5×、p95不劣、RSS≤1.5× current gate。</li>
              <li class="stop">Production feature flag、fallback與release QA。</li>
              <li class="stop">M2 copy/loss model若要進CNA/LOH區域。</li>
            </ul>
          </div>
        </div>
        <div class="callout stop">
          <strong>最終操作判定：</strong>不可變bounded evidence bundle完成後，本輪只授權
          「以另一份 frozen input contract 跑33-tail／H2009／H1437 dual stress」。
          禁止直接替換canonical、啟動新的production全量run，或宣稱生物ancestry改善。
        </div>
        <details><summary>狀態字典</summary><div>
          <ul>
            <li><code>OPTIMAL_VALUE_CERTIFIED</code>：只證明h*，不代表family完整。</li>
            <li><code>CANDIDATE_SET_COMPLETE</code>：所有minimum vertex sets已列完。</li>
            <li><code>CANDIDATE_SET_INCOMPLETE_*</code>：h*可已知，但family不可排名。</li>
            <li><code>FEASIBLE_UNPROVEN_*</code>：有incumbent，但h*未證明。</li>
            <li><code>NO_FEASIBLE_CERTIFICATE_*</code>：限制先到，連feasible incumbent都沒有。</li>
            <li><code>ABSTAIN_CN_LOH_GATE</code>：M1輸入資格不完整或有矛盾證據。</li>
          </ul>
        </div></details>
      </section>
    </main>
  </div>

  <footer>
    <p><strong>InterSubMod bounded solver audit</strong> · generated from source-bound JSON receipts.</p>
    <p class="hash-note">Benchmark receipt SHA-256: $benchmark_receipt_sha<br>Oracle receipt SHA-256: $oracle_receipt_sha</p>
    <p>這份報告是PARTIAL工程證據，不是production validation或biological truth report。</p>
  </footer>

  <script type="application/json" id="report-evidence">$receipt_json</script>
  <script>
    (function () {
      const root = document.documentElement;
      const toggle = document.getElementById("themeToggle");
      toggle.addEventListener("click", function () {
        root.dataset.theme = root.dataset.theme === "night" ? "paper" : "night";
      });
      document.getElementById("expandAll").addEventListener("click", function () {
        document.querySelectorAll("details").forEach(function (node) { node.open = true; });
      });
      document.getElementById("collapseAll").addEventListener("click", function () {
        document.querySelectorAll("details").forEach(function (node) { node.open = false; });
      });
      const links = Array.from(document.querySelectorAll("nav.toc a"));
      const sections = links.map(function (link) { return document.querySelector(link.getAttribute("href")); });
      if ("IntersectionObserver" in window) {
        const observer = new IntersectionObserver(function (entries) {
          entries.forEach(function (entry) {
            if (!entry.isIntersecting) return;
            links.forEach(function (link) { link.classList.remove("active"); });
            const active = links.find(function (link) { return link.getAttribute("href") === "#" + entry.target.id; });
            if (active) active.classList.add("active");
          });
        }, { rootMargin:"-20% 0px -65% 0px", threshold:0 });
        sections.forEach(function (section) { if (section) observer.observe(section); });
      }
    }());
  </script>
</body>
</html>
""")
    speedups = [row["median_speedup_solver"] for row in comparisons]
    rss_ratios = [
        row["median_rss_ratio_optimized_over_current"] for row in comparisons
    ]
    return template.substitute(
        benchmark_receipt_short=esc(short_hash(bench_receipt_sha)),
        oracle_receipt_short=esc(short_hash(oracle_receipt_sha)),
        benchmark_receipt_sha=esc(bench_receipt_sha),
        oracle_receipt_sha=esc(oracle_receipt_sha),
        generated=esc(generated),
        speed_range=esc(f"{min(speedups):.2f}×–{max(speedups):.2f}×"),
        rss_range=esc(f"{min(rss_ratios):.3f}×–{max(rss_ratios):.3f}×"),
        oracle_cases=f"{oracle['summary']['total_structural_cases']:,}",
        source_rows=sources,
        plan_sha=esc(plan_sha),
        notes_sha=esc(notes_sha),
        platform=esc(platform.platform()),
        timing_bars=bars,
        benchmark_rows=benchmark_table,
        commands=esc(commands),
        receipt_json=receipt_json,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", type=pathlib.Path, required=True)
    parser.add_argument("--oracle", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    benchmark_path = args.benchmark.resolve()
    oracle_path = args.oracle.resolve()
    benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
    oracle = json.loads(oracle_path.read_text(encoding="utf-8"))
    assert_receipts(benchmark, oracle, benchmark_path, oracle_path)
    document = build_html(benchmark, oracle, benchmark_path, oracle_path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(document, encoding="utf-8")
    print(
        json.dumps(
            {
                "status": "PASS",
                "output": str(args.output.resolve()),
                "bytes": args.output.stat().st_size,
                "sha256": sha256_file(args.output),
                "sections": 10,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
