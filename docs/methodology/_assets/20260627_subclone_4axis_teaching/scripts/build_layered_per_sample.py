#!/usr/bin/env python3
"""
Build the portable, per-sample layered reconstruction workstation.

The landing page is always generated from each sample's ``census`` object. Use
``--index-only`` to refresh the cohort landing page without rebuilding the large
sample HTML files.

Inputs:
  - layered_region_view_{sample}.json for every configured sample
Outputs:
  - layered_workstation/{sample}.html (default mode only)
  - layered_workstation/index.html
"""

import argparse
import html
import json
import os
import subprocess
from datetime import datetime


HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "build_layered_workstation.py")
PY = "python3"
OUTDIR = os.path.normpath(os.path.join(HERE, "..", "..", "layered_workstation"))
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"

# Sample order is an interface choice; every displayed measurement comes from census.
SAMPLES = [
    (
        "HCC1395",
        "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/"
        "20260618_subcluster_pilot/layered_region_view_HCC1395.json",
        "SEQC2 truth · primary 5 kHz run",
    ),
    ("COLO829", f"{MSROOT}/COLO829/layered_region_view_COLO829.json", "melanoma · external truth"),
    ("H1437", f"{MSROOT}/H1437/layered_region_view_H1437.json", "LUAD"),
    ("H2009", f"{MSROOT}/H2009/layered_region_view_H2009.json", "LUAD"),
    (
        "HCC1395_DORADO",
        f"{MSROOT}/HCC1395_DORADO/layered_region_view_HCC1395_DORADO.json",
        "same cell line · Dorado basecaller",
    ),
    ("HCC1937", f"{MSROOT}/HCC1937/layered_region_view_HCC1937.json", "breast"),
    ("HCC1954", f"{MSROOT}/HCC1954/layered_region_view_HCC1954.json", "breast"),
]


def human(n_bytes):
    """Format a file size for display."""
    value = float(n_bytes)
    for unit in ("B", "KB", "MB"):
        if value < 1024:
            return f"{value:.0f} {unit}"
        value /= 1024
    return f"{value:.1f} GB"


def local_time(timestamp):
    """Render a local timestamp with timezone for provenance."""
    return datetime.fromtimestamp(timestamp).astimezone().isoformat(timespec="minutes")


def escaped(value):
    return html.escape(str(value), quote=True)


def repo_display_path(value):
    """Use the repository prefix required for displayed Markdown paths."""
    path = str(value)
    return f"InterSubMod/{path}" if path.startswith("docs/") else path


def collect_rows(build_samples):
    """Build sample pages if requested and collect landing-page rows."""
    rows = []
    for name, region_view, note in SAMPLES:
        output = os.path.join(OUTDIR, name + ".html")
        row = {
            "name": name,
            "note": note,
            "region_view": region_view,
            "output": output,
            "source_ready": os.path.exists(region_view),
            "page_ready": os.path.exists(output),
            "build_error": False,
        }
        if not row["source_ready"]:
            rows.append(row)
            continue

        if build_samples:
            env = dict(os.environ, SM_RV=region_view, SM_OUT=output)
            result = subprocess.run(
                [PY, BUILD],
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if result.returncode != 0:
                row["build_error"] = True
                print(f"[{name}] BUILD FAIL:\n{result.stdout.decode(errors='replace')[:500]}")
                rows.append(row)
                continue
            row["page_ready"] = True

        with open(region_view, encoding="utf-8") as handle:
            payload = json.load(handle)
        census = payload["census"]
        determinacy = census["region_determinacy"]
        multiplicity = census["hp_multiplicity"]
        layer_one = census["L1"]
        n_lineages = layer_one["n_lineage_units"]

        row.update(
            {
                "ssnv": census.get("U1_sSNV_somatic_total", 0),
                "n_regions": census["n_regions"],
                "multi_hp": multiplicity.get("2", 0),
                "determined": determinacy.get("all_determined", 0),
                "ambiguous": determinacy.get("has_ambiguous", 0),
                "capped": determinacy.get("has_capped", 0),
                "recurrence": determinacy.get("has_recurrence", 0),
                "avg_trees": (
                    census.get("U5_trees", {}).get("sum_ntrees_noncapped", 0) / n_lineages
                    if n_lineages
                    else 0
                ),
                "v1_v7_pass": layer_one.get("all_V1V7_pass", False),
                "l3_ready": bool(census.get("L3")),
                "backbone": census.get("U1_backbone_source", "not recorded"),
                "data_model": payload.get("data_model_spec", "not recorded"),
                "source_mtime": os.path.getmtime(region_view),
                "size": os.path.getsize(output) if os.path.exists(output) else 0,
            }
        )
        rows.append(row)
        if row["page_ready"]:
            print(f"[{name}] {'INDEX' if not build_samples else 'OK'} {human(row['size'])}")
        else:
            print(f"[{name}] census ready; sample HTML is not built")
    return rows


def percent_cell(value, total, color, label):
    percentage = 100 * value / total if total else 0
    return (
        f'<td class="numeric" data-label="{escaped(label)}">'
        f'<span class="metric-value" style="color:{color}">{percentage:.1f}%</span>'
        f'<span class="metric-count">{value:,}</span>'
        f'<span class="metric-track" aria-hidden="true">'
        f'<span style="width:{percentage:.1f}%;background:{color}"></span></span></td>'
    )


def sample_rows_html(rows):
    rendered = []
    for row in rows:
        name = escaped(row["name"])
        note = escaped(row["note"])
        if not row["source_ready"]:
            rendered.append(
                '<tr class="unavailable">'
                f'<td class="sample-cell"><strong>{name}</strong><small>{note}</small></td>'
                '<td colspan="9"><span class="status status-muted">UNAVAILABLE</span> '
                'No layered region-view census was found.</td><td aria-label="Open">—</td></tr>'
            )
            continue
        if row["build_error"]:
            rendered.append(
                '<tr class="unavailable">'
                f'<td class="sample-cell"><strong>{name}</strong><small>{note}</small></td>'
                '<td colspan="9"><span class="status status-review">BUILD ERROR</span> '
                'The census exists, but the sample page could not be refreshed.</td><td aria-label="Open">—</td></tr>'
            )
            continue

        if row["page_ready"]:
            name_markup = f'<a href="{name}.html"><strong>{name}</strong></a>'
            action = f'<a class="button" href="{name}.html">Open sample<span aria-hidden="true"> →</span></a>'
        else:
            name_markup = f'<strong>{name}</strong>'
            action = '<span class="status status-muted">NOT BUILT</span>'
        validation = "PASS" if row["v1_v7_pass"] else "REVIEW"
        validation_class = "status-pass" if row["v1_v7_pass"] else "status-review"
        l3_label = "AVAILABLE" if row["l3_ready"] else "L3 PENDING"
        l3_class = "status-pass" if row["l3_ready"] else "status-pending"
        rendered.append(
            '<tr>'
            f'<td class="sample-cell">{name_markup}<small>{note}</small></td>'
            f'<td class="numeric" data-label="Somatic sSNV">{row["ssnv"]:,}</td>'
            f'<td class="numeric" data-label="Eligible regions">{row["n_regions"]:,}</td>'
            + percent_cell(row["multi_hp"], row["n_regions"], "#b42373", "HP multiplicity equals two")
            + percent_cell(row["determined"], row["n_regions"], "#117a56", "All determined")
            + percent_cell(row["ambiguous"], row["n_regions"], "#b35c00", "Ambiguous")
            + percent_cell(row["capped"], row["n_regions"], "#7551a5", "Capped")
            + f'<td class="numeric" data-label="Recurrence">{row["recurrence"]:,}</td>'
            + f'<td class="numeric" data-label="Average trees">{row["avg_trees"]:.2f}</td>'
            + f'<td data-label="Validation and file"><span class="status {validation_class}">V1–V7 {validation}</span>'
            + f'<span class="status {l3_class}">{l3_label}</span><small>{human(row["size"]) if row["size"] else "—"}</small></td>'
            + f'<td data-label="Open">{action}</td></tr>'
        )
    return "".join(rendered)


def cohort_metadata(rows):
    available = [row for row in rows if row.get("n_regions") is not None]
    complete = [row for row in available if row["page_ready"]]
    backbones = sorted({row["backbone"] for row in available})
    data_models = sorted({repo_display_path(row["data_model"]) for row in available})
    l3_count = sum(1 for row in available if row["l3_ready"])
    latest_source = max((row["source_mtime"] for row in available), default=None)
    return {
        "available": available,
        "complete": complete,
        "backbone": " · ".join(backbones) if backbones else "not available",
        "data_model": " · ".join(data_models) if data_models else "not available",
        "l3_count": l3_count,
        "source_snapshot": local_time(latest_source) if latest_source is not None else "not available",
        "generated_at": datetime.now().astimezone().isoformat(timespec="minutes"),
    }


def build_index(rows):
    meta = cohort_metadata(rows)
    n_configured = len(rows)
    n_available = len(meta["available"])
    n_complete = len(meta["complete"])
    l3_status = (
        "AVAILABLE"
        if n_available and meta["l3_count"] == n_available
        else "PARTIAL"
        if meta["l3_count"]
        else "PENDING"
    )
    table_rows = sample_rows_html(rows)
    backbone = escaped(meta["backbone"])
    data_model = escaped(meta["data_model"])
    source_snapshot = escaped(meta["source_snapshot"])
    generated_at = escaped(meta["generated_at"])

    return f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<title>Layered reconstruction · cohort console</title>
<style>
:root{{--ink:#17212b;--muted:#596775;--line:#d7dee6;--paper:#f4f7f9;--panel:#ffffff;--blue:#1769aa;--hp1:#1971c2;--hp2:#b42373;--hp3:#087f8c;--good:#117a56;--warn:#9a5200;--pending:#6d4b8f;--shadow:0 12px 34px rgba(25,44,64,.07)}}
*{{box-sizing:border-box}}
html{{background:var(--paper);scroll-behavior:smooth}}
body{{margin:0;overflow-x:hidden;color:var(--ink);background:linear-gradient(180deg,#eef3f7 0,#f7f9fb 360px,var(--paper) 100%);font-family:"Aptos","Noto Sans TC","Microsoft JhengHei",sans-serif;line-height:1.55}}
a{{color:var(--blue);text-decoration-thickness:1px;text-underline-offset:3px}}
a:hover{{text-decoration-thickness:2px}}
a:focus-visible,button:focus-visible{{outline:3px solid #f1a93b;outline-offset:3px}}
.skip-link{{position:fixed;left:12px;top:8px;z-index:30;transform:translateY(-160%);background:#fff;color:#111;padding:8px 12px;border:2px solid #111}}
.skip-link:focus{{transform:none}}
.wrap{{width:min(100%,1280px);margin:0 auto;padding:24px clamp(14px,3vw,38px) 46px}}
.hero{{position:relative;overflow:hidden;border:1px solid #cbd5df;border-top:4px solid var(--blue);background:rgba(255,255,255,.96);box-shadow:var(--shadow)}}
.hero:before{{content:"";position:absolute;inset:0 0 auto auto;width:220px;height:100%;pointer-events:none;background:repeating-linear-gradient(90deg,transparent 0 15px,rgba(23,105,170,.05) 15px 16px)}}
.hero-main{{position:relative;padding:clamp(22px,4vw,42px)}}
.eyebrow{{display:flex;align-items:center;gap:10px;margin:0 0 16px;color:#365064;font:700 11px/1.2 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;letter-spacing:.12em;text-transform:uppercase}}
.signal{{width:9px;height:9px;border-radius:50%;background:var(--good);box-shadow:0 0 0 4px rgba(17,122,86,.13)}}
h1{{max-width:880px;margin:0;font-size:clamp(29px,4.6vw,54px);font-weight:650;letter-spacing:-.035em;line-height:1.04}}
.lede{{max-width:850px;margin:18px 0 0;color:#3e5060;font-size:clamp(15px,1.8vw,19px)}}
.hp-rail{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));height:5px;margin-top:28px}}
.hp-rail span:nth-child(1){{background:var(--hp1)}}.hp-rail span:nth-child(2){{background:var(--hp2)}}.hp-rail span:nth-child(3){{background:var(--hp3)}}.hp-rail span:nth-child(4){{background:#8995a1}}
.provenance{{display:grid;grid-template-columns:1.45fr 1fr 1fr 1fr;border-top:1px solid var(--line);background:#f8fafb}}
.provenance div{{min-width:0;padding:12px 16px;border-right:1px solid var(--line)}}
.provenance div:last-child{{border-right:0}}
.provenance dt{{margin:0 0 3px;color:var(--muted);font:700 10px/1.2 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;letter-spacing:.08em;text-transform:uppercase}}
.provenance dd{{margin:0;overflow-wrap:anywhere;font-size:12px;font-weight:650}}
.section{{margin-top:28px}}
.section-head{{display:flex;align-items:flex-end;justify-content:space-between;gap:16px;margin-bottom:10px}}
.section-kicker{{margin:0 0 3px;color:var(--blue);font:700 10.5px/1.2 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;letter-spacing:.1em;text-transform:uppercase}}
h2{{margin:0;font-size:clamp(20px,2.5vw,28px);letter-spacing:-.02em}}
.section-note{{max-width:580px;margin:0;color:var(--muted);font-size:12.5px}}
.claim-grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}}
.claim{{min-height:184px;padding:18px;border:1px solid var(--line);background:var(--panel);box-shadow:0 5px 18px rgba(25,44,64,.035)}}
.claim:nth-child(1){{border-top:4px solid var(--hp1)}}.claim:nth-child(2){{border-top:4px solid var(--hp3)}}.claim:nth-child(3){{border-top:4px solid var(--warn)}}
.claim-label{{display:block;margin-bottom:18px;color:var(--muted);font:700 10px/1.2 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;letter-spacing:.11em;text-transform:uppercase}}
.claim h3{{margin:0 0 8px;font-size:17px;line-height:1.25}}
.claim p{{margin:0;color:#4a5966;font-size:13.5px}}
.instrument-strip{{display:grid;grid-template-columns:1fr 1fr 1fr;margin-top:12px;border:1px solid var(--line);background:#17212b;color:#edf4f8}}
.readout{{min-width:0;padding:15px 18px;border-right:1px solid #35424e}}
.readout:last-child{{border-right:0}}
.readout span{{display:block;color:#aebdca;font:700 10px/1.2 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;letter-spacing:.09em;text-transform:uppercase}}
.readout strong{{display:block;margin-top:7px;overflow-wrap:anywhere;font:650 18px/1.2 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;font-variant-numeric:tabular-nums}}
.readout.pending strong{{color:#dfc7f2}}
.model-line{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));border:1px solid var(--line);background:#fff}}
.layer{{position:relative;min-width:0;padding:16px 18px;border-right:1px solid var(--line)}}
.layer:last-child{{border-right:0}}
.layer:not(:last-child):after{{content:"→";position:absolute;z-index:2;right:-10px;top:50%;transform:translateY(-50%);display:grid;place-items:center;width:19px;height:19px;border:1px solid var(--line);border-radius:50%;background:#fff;color:var(--muted);font-size:11px}}
.layer-code{{color:var(--blue);font:750 11px/1 ui-monospace,"SFMono-Regular","Cascadia Code",monospace}}
.layer strong{{display:block;margin:7px 0 3px;font-size:13px}}
.layer small{{display:block;color:var(--muted);font-size:11.5px}}
.layer.is-pending{{background:#faf7fc}}
.layer.is-pending .layer-code{{color:var(--pending)}}
.table-shell{{max-width:100%;border:1px solid #cbd5df;background:#fff;box-shadow:var(--shadow)}}
.scroll-cue{{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:9px 12px;border-bottom:1px solid var(--line);background:#f8fafb;color:var(--muted);font-size:11.5px}}
.table-scroll{{width:100%;max-width:100%;overflow-x:auto;overscroll-behavior-inline:contain;-webkit-overflow-scrolling:touch}}
.table-scroll:focus-visible{{outline:3px solid #f1a93b;outline-offset:-3px}}
table{{width:100%;min-width:1160px;border-collapse:collapse;font-size:12.5px}}
th,td{{padding:11px 10px;text-align:left;vertical-align:middle;border-bottom:1px solid #e5eaf0}}
thead th{{position:sticky;top:0;z-index:2;background:#eaf0f4;color:#344656;font:750 10.5px/1.25 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;letter-spacing:.025em}}
tbody tr:last-child td{{border-bottom:0}}
tbody tr:hover td{{background:#f8fbfd}}
th:first-child,td:first-child{{position:sticky;left:0;z-index:1}}
thead th:first-child{{z-index:3;background:#eaf0f4}}
tbody td:first-child{{background:#fff;box-shadow:6px 0 10px -10px #17212b}}
tbody tr:hover td:first-child{{background:#f8fbfd}}
.sample-cell{{min-width:180px}}
.sample-cell strong{{display:block;font-size:13px}}
.sample-cell small,td small{{display:block;margin-top:3px;color:var(--muted);font-size:11px}}
.numeric{{text-align:right;font-family:ui-monospace,"SFMono-Regular","Cascadia Code",monospace;font-variant-numeric:tabular-nums}}
.metric-value{{display:block;font-weight:800}}
.metric-count{{display:block;margin-top:1px;color:var(--muted);font-size:10.5px}}
.metric-track{{display:block;width:72px;height:3px;margin:6px 0 0 auto;overflow:hidden;background:#e3e9ee}}
.metric-track span{{display:block;height:100%}}
.status{{display:inline-block;margin:0 5px 4px 0;padding:3px 6px;border:1px solid;border-radius:2px;font:750 9.5px/1.2 ui-monospace,"SFMono-Regular","Cascadia Code",monospace;letter-spacing:.035em;white-space:nowrap}}
.status-pass{{border-color:#93c9b7;color:#0b6848;background:#edf8f4}}.status-review{{border-color:#e0b879;color:#7a4000;background:#fff8ec}}.status-pending{{border-color:#c9b1dc;color:#62417e;background:#faf6fd}}.status-muted{{border-color:#c9d0d6;color:#55616c;background:#f4f6f8}}
.button{{display:inline-block;padding:7px 10px;border:1px solid var(--blue);background:var(--blue);color:#fff;text-decoration:none;font-size:11.5px;font-weight:750;white-space:nowrap}}
.button:hover{{background:#0f578f;text-decoration:none}}
.unavailable{{color:var(--muted)}}
.definition-grid{{display:grid;grid-template-columns:1.2fr 1fr;gap:12px}}
.definition-card{{border:1px solid var(--line);background:#fff;padding:18px}}
.definition-card h3{{margin:0 0 10px;font-size:14px}}
.definition-card p,.definition-card li{{color:#465663;font-size:12.5px}}
.definition-card ul{{margin:0;padding-left:18px}}
.definition-card li+li{{margin-top:6px}}
.rule{{border-left:4px solid var(--warn)}}
.mono{{font-family:ui-monospace,"SFMono-Regular","Cascadia Code",monospace;overflow-wrap:anywhere}}
footer{{margin-top:25px;padding-top:14px;border-top:1px solid var(--line);color:var(--muted);font-size:11.5px}}
@media(max-width:820px){{.provenance{{grid-template-columns:1fr 1fr}}.provenance div:nth-child(2){{border-right:0}}.provenance div:nth-child(-n+2){{border-bottom:1px solid var(--line)}}.claim-grid,.definition-grid{{grid-template-columns:1fr}}.claim{{min-height:0}}.instrument-strip{{grid-template-columns:1fr}}.readout{{border-right:0;border-bottom:1px solid #35424e}}.readout:last-child{{border-bottom:0}}.model-line{{grid-template-columns:1fr 1fr}}.layer:nth-child(2){{border-right:0}}.layer:nth-child(-n+2){{border-bottom:1px solid var(--line)}}.layer:nth-child(2):after{{display:none}}}}
@media(max-width:520px){{.wrap{{padding:14px 12px 32px}}.hero-main{{padding:22px 18px}}.provenance{{grid-template-columns:1fr}}.provenance div{{border-right:0;border-bottom:1px solid var(--line)}}.provenance div:last-child{{border-bottom:0}}.section-head{{display:block}}.section-note{{margin-top:6px}}.model-line{{grid-template-columns:1fr}}.layer{{border-right:0;border-bottom:1px solid var(--line)}}.layer:last-child{{border-bottom:0}}.layer:after{{display:none!important}}.scroll-cue span:last-child{{display:none}}}}
@media(prefers-reduced-motion:reduce){{html{{scroll-behavior:auto}}}}
@media print{{body{{background:#fff}}.wrap{{max-width:none;padding:0}}.hero,.table-shell{{box-shadow:none}}.table-scroll{{overflow:visible}}table{{min-width:0;font-size:8.5px}}.scroll-cue,.button{{display:none}}}}
</style>
</head>
<body>
<a class="skip-link" href="#cohort-table">Skip to cohort table</a>
<main class="wrap">
  <header class="hero">
    <div class="hero-main">
      <p class="eyebrow"><span class="signal" aria-hidden="true"></span>Current · canonical research interface</p>
      <h1>Layered reconstruction<br>cohort console</h1>
      <p class="lede">以 germline HP family 為邊界，逐家族枚舉所有相容的最小 sSNV 狀態樹；唯一、歧義與計算上限都保留為可稽核結果。</p>
      <div class="hp-rail" aria-label="HP family palette"><span title="HP1"></span><span title="HP2"></span><span title="HP3"></span><span title="Unassigned"></span></div>
    </div>
    <dl class="provenance">
      <div><dt>Data version / backbone</dt><dd>{backbone}</dd></div>
      <div><dt>Data snapshot</dt><dd>{source_snapshot}</dd></div>
      <div><dt>Page generated</dt><dd>{generated_at}</dd></div>
      <div><dt>Model specification</dt><dd class="mono">{data_model}</dd></div>
    </dl>
  </header>

  <section class="section" aria-labelledby="claim-title">
    <div class="section-head"><div><p class="section-kicker">Interpretation boundary</p><h2 id="claim-title">先理解這個工作站能回答到哪裡</h2></div><p class="section-note">每張圖都應回到「觀察、推論、限制」三層；不能唯一決定不是失敗，而是資料可識別性的結果。</p></div>
    <div class="claim-grid">
      <article class="claim"><span class="claim-label">01 · 回答什麼</span><h3>哪些區域樹被資料允許？</h3><p>對每個 eligible multilocus region，分 HP family 回答相容最小樹是否唯一、存在多個等價解，或因 dense state space 被 capped。</p></article>
      <article class="claim"><span class="claim-label">02 · 證據是什麼</span><h3>read-level sSNV 共現狀態</h3><p>L0 HP family 與 L1 單分子共現是樹骨幹；solver 枚舉全部最小解。L2 CN 只在有資料時提供事後裁決，V1–V7 記錄內部一致性。</p></article>
      <article class="claim"><span class="claim-label">03 · 不能宣稱什麼</span><h3>不是已確認的生物 subclone</h3><p>不能由此單獨宣稱全腫瘤 clone phylogeny、細胞比例或正交確認。L3 methylation 目前 pending，且不參與樹的排序或選擇。</p></article>
    </div>
    <div class="instrument-strip" aria-label="Cohort status readout">
      <div class="readout"><span>Sample workstations</span><strong>{n_complete} / {n_configured}</strong></div>
      <div class="readout"><span>Census inputs available</span><strong>{n_available} / {n_configured}</strong></div>
      <div class="readout pending"><span>L3 methylation census</span><strong>{l3_status} · {meta["l3_count"]} / {n_available}</strong></div>
    </div>
  </section>

  <section class="section" aria-labelledby="model-title">
    <div class="section-head"><div><p class="section-kicker">Evidence stack</p><h2 id="model-title">L0 → L3 分層軌跡</h2></div><p class="section-note">遺傳證據承重；後續層只能在可用時追加註解或裁決，不回頭重寫觀察到的 read state。</p></div>
    <div class="model-line">
      <div class="layer"><span class="layer-code">L0</span><strong>HP family</strong><small>先拆 germline lineage，避免 pooled allelic/clonal 混淆。</small></div>
      <div class="layer"><span class="layer-code">L1</span><strong>sSNV enumeration</strong><small>用 read-level 共現狀態枚舉所有最小樹。</small></div>
      <div class="layer"><span class="layer-code">L2</span><strong>Copy-number context</strong><small>資料可用時裁決 recurrence；不可用即保留 unavailable。</small></div>
      <div class="layer is-pending"><span class="layer-code">L3 · {l3_status}</span><strong>Methylation auxiliary</strong><small>目前不 rank、不選樹；待有效 census 才啟用展示。</small></div>
    </div>
  </section>

  <section class="section" id="cohort-table" aria-labelledby="cohort-title">
    <div class="section-head"><div><p class="section-kicker">Cross-sample inventory</p><h2 id="cohort-title">跨樣本 census</h2></div><p class="section-note">比例的分母都是各樣本 eligible region；絕對數僅作資料盤點，不能直接當生物差異。</p></div>
    <div class="table-shell">
      <div class="scroll-cue"><span>表格可水平捲動；樣本欄固定。</span><span>Region-level denominator</span></div>
      <div class="table-scroll" role="region" aria-label="Cross-sample layered reconstruction census" tabindex="0">
        <table>
          <thead><tr><th>Sample</th><th class="numeric">Somatic<br>sSNV</th><th class="numeric">Eligible<br>regions</th><th class="numeric">HP multiplicity<br>= 2, %</th><th class="numeric">All<br>determined %</th><th class="numeric">Has<br>ambiguous %</th><th class="numeric">Has<br>capped %</th><th class="numeric">Recurrence</th><th class="numeric">Avg trees /<br>lineage</th><th>Validation · L3 · file</th><th>Inspect</th></tr></thead>
          <tbody>{table_rows}</tbody>
        </table>
      </div>
    </div>
  </section>

  <section class="section" aria-labelledby="reading-title">
    <div class="section-head"><div><p class="section-kicker">Reading protocol</p><h2 id="reading-title">比較前先鎖定分母與證據層</h2></div></div>
    <div class="definition-grid">
      <article class="definition-card">
        <h3>欄位讀法</h3>
        <ul>
          <li><strong>HP multiplicity = 2</strong>：目前 census 記錄兩個 germline family units 的 region 比例；可能包含 root-only control，<strong>不等同 both-mutation-bearing</strong>。</li>
          <li><strong>All determined</strong>：該 region 的所有 census lineages 都只有一棵最小樹；不表示每個 lineage 都有 mutation edge。</li>
          <li><strong>Has ambiguous</strong>：至少一個 lineage 有多棵等價最小樹；完整保留解集合。</li>
          <li><strong>Has capped</strong>：隱藏祖先空間超過列舉上限；不是唯一樹證據。</li>
        </ul>
      </article>
      <article class="definition-card rule">
        <h3>不可跨用的分母</h3>
        <p><strong>Region %、lineage-unit % 與舊 pooled % 不可互換。</strong>工作站首頁只呈現 region-level 比例；每個樣本頁再展開 U1–U6 單位與 V1–V7 檢核。</p>
        <p class="mono">InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md</p>
      </article>
    </div>
  </section>

  <footer>Offline portable artifact · zero external dependencies · index regenerated from layered region-view census with <span class="mono">python3 build_layered_per_sample.py --index-only</span>.</footer>
</main>
</body>
</html>"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--index-only",
        action="store_true",
        help="Regenerate index.html from census inputs without rebuilding sample HTML files.",
    )
    args = parser.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    rows = collect_rows(build_samples=not args.index_only)
    index_path = os.path.join(OUTDIR, "index.html")
    with open(index_path, "w", encoding="utf-8") as handle:
        handle.write(build_index(rows))
    complete = sum(1 for row in rows if row["source_ready"] and row["page_ready"] and not row["build_error"])
    print(f"OK index.html + {complete}/{len(rows)} sample pages -> {OUTDIR}")


if __name__ == "__main__":
    main()
